#!/usr/bin/env python3
"""
Unified pangenome middleware.

This is the single in-process API a service layer (e.g. an HTTP server for the
UCSC Genome Browser tool) should call. It combines the two engines:

  * coordinate translation + anchor building via liftover_ext.Index (in-process)
  * read mapping + anchor-driven surjection via `vg giraffe-server`
    (a long-lived subprocess)

It is kept at parity with the reference orchestration in
`smoke_test_middleware.py`, in particular the fast **anchor** surjection
pipeline (map -> build_surject_anchors -> surject_with_anchors), which surjects
onto ANY haplotype at query time without pre-indexing it.

Preset note: the long-read / BLAT mapping preset is applied inside the giraffe
engine itself (apply_blat_preset in giraffe_engine.cpp); the middleware does not
pass `--parameter-preset`. The only mapping knob exposed here is
`max_multimaps` (0 = the engine's built-in BLAT default of 100; 1 = single best
alignment per read, recommended for the one-line-per-read UI). Pass the
LONG-READ minimizer + zipcodes to match the engine's preset.
"""
from __future__ import annotations

import threading
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from .giraffe_server_middleware import (
    FastqRead,
    GiraffeServerConfig,
    GiraffeServerMiddleware,
)

# `liftover_ext` is the native pybind11 extension. It is imported lazily inside
# __init__ (not at module import) so this module — and the HTTP API / stub built
# on it — can be imported on hosts that don't have the extension built.


@dataclass
class CoordinateIndexPaths:
    gbz_path: str
    ri_path: str
    tags_path: str
    gbwt_ri_path: str
    table1_path: str
    table2_path: str = ""      # optional: empty = table-free translation


# Result of building anchors for one graph alignment onto a target haplotype:
# (anchors, target_path_length, status).
AnchorBuild = Tuple[list, int, str]


def _fold_to_intervals(raw) -> List[Dict[str, Any]]:
    """Fold liftover_ext's per-base source→target correspondences into target
    intervals with strand.

    Each raw point (a liftover_ext.TranslatedInterval) carries:
      - .haplotype : the full target CONTIG path it landed on (e.g. "H#2#CM09.1")
      - .start     : the SOURCE coordinate (contig-local)
      - .end       : the TARGET coordinate (contig-local)
    (Its own .strand is a stubbed '+', so we ignore it and derive strand here.)

    trace_coordinates_gbwt emits one point per source base, so points on a
    colinear stretch step target by ±1. We group by target contig and split
    each contig's points (in source order) into maximal runs of a single target
    direction; each run becomes one interval:
      - start/end : min/max TARGET coordinate of the run (half-open, +1 on end)
      - strand    : '+' if target rises with source, '-' if it falls
    A direction reversal (e.g. an inversion) or a different contig naturally
    breaks into separate pieces — this is what yields the agreed 0..N intervals.
    """
    by_contig: Dict[str, List[Tuple[int, int]]] = {}
    for p in raw:
        by_contig.setdefault(p.haplotype, []).append((int(p.start), int(p.end)))

    out: List[Dict[str, Any]] = []
    for contig, pts in by_contig.items():
        pts.sort()  # source ascending (ties: target ascending)
        i, n = 0, len(pts)
        while i < n:
            j = i
            direction = 0  # +1 target rising, -1 falling, 0 undecided
            while j + 1 < n:
                dt = pts[j + 1][1] - pts[j][1]
                if dt == 0:
                    j += 1  # plateau (e.g. source insertion): stay in the run
                    continue
                step = 1 if dt > 0 else -1
                if direction == 0:
                    direction = step
                if step != direction:
                    break  # target reversed → end this run, start a new piece
                j += 1
            run_targets = [t for _s, t in pts[i:j + 1]]
            lo, hi = min(run_targets), max(run_targets)
            out.append({
                "haplotype": contig,
                "start": lo,
                "end": hi + 1,
                "strand": "-" if direction < 0 else "+",
            })
            i = j + 1
    return out


def _fold_to_blocks(raw) -> List[Dict[str, Any]]:
    """Fold per-base correspondences into BLOCK-LEVEL alignment.

    Where _fold_to_intervals() reports one min/max span per direction-run — so an
    indel inside a region is spanned straight across, and exon structure is lost —
    this emits every maximal colinear block, each carrying BOTH sides:

        {haplotype, source_start, source_end, target_start, target_end, strand}

    A block ends wherever the 1:1 correspondence breaks:
      * source advances by more than 1  -> bases with no target (deletion)
      * target advances by more than 1  -> bases inserted in the target
      * target direction flips          -> inversion
      * different target contig
    So a gene lifts over as its exon/indel structure rather than one fused span,
    and one call per (region, target) replaces one call per (feature, target).
    All coordinates are 0-based half-open; target_start < target_end always, with
    orientation carried by `strand`.
    """
    by_contig: Dict[str, List[Tuple[int, int]]] = {}
    for p in raw:
        by_contig.setdefault(p.haplotype, []).append((int(p.start), int(p.end)))

    out: List[Dict[str, Any]] = []
    for contig, pts in by_contig.items():
        pts.sort()
        i, n = 0, len(pts)
        while i < n:
            s0, t0 = pts[i]
            j = i
            direction = 0
            while j + 1 < n:
                s_prev, t_prev = pts[j]
                s_next, t_next = pts[j + 1]
                if s_next - s_prev != 1:
                    break                      # unmapped source bases
                dt = t_next - t_prev
                if dt == 1:
                    step = 1
                elif dt == -1:
                    step = -1
                else:
                    break                      # jump on the target side
                if direction == 0:
                    direction = step
                elif step != direction:
                    break                      # orientation change
                j += 1
            s1, t1 = pts[j]
            tgt_lo, tgt_hi = (t0, t1) if direction >= 0 else (t1, t0)
            out.append({
                "haplotype": contig,
                "source_start": s0, "source_end": s1 + 1,
                "target_start": tgt_lo, "target_end": tgt_hi + 1,
                "strand": "-" if direction < 0 else "+",
            })
            i = j + 1
    out.sort(key=lambda b: (b["haplotype"], b["source_start"]))
    return out


class PangenomeMiddleware:
    """
    Unified middleware:
      - coordinate translation + anchor building via liftover_ext.Index (in-process)
      - read mapping + anchor-driven surjection via vg giraffe-server (subprocess)

    Typical lifecycle:
        mw = PangenomeMiddleware.from_paths(vg_binary=..., gbz=..., ...)
        mw.wait_until_ready()                 # block until indexes are loaded
        lines = mw.surject_sequences([...], target="CHM13#0#chr10")
        mw.close()
    """

    def __init__(self, coord_paths: CoordinateIndexPaths, giraffe_cfg: GiraffeServerConfig) -> None:
        import liftover_ext  # lazy: only the real (indexed) path needs the extension
        self._coord = liftover_ext.Index()
        self._coord.load(
            coord_paths.gbz_path,
            coord_paths.ri_path,
            coord_paths.tags_path,
            coord_paths.gbwt_ri_path,
            coord_paths.table1_path,
            coord_paths.table2_path,
        )
        # Serializes coordinate-index queries: the liftover_ext query path is not
        # guaranteed re-entrant (lazy index inits, const_cast reads), and the API
        # can call translate concurrently with mapping. Translation is fast and
        # low-volume, so a lock is cheap insurance.
        self._coord_lock = threading.Lock()
        self._giraffe = GiraffeServerMiddleware(giraffe_cfg)
        self._giraffe.start()

    @classmethod
    def from_paths(
        cls,
        *,
        vg_binary: str,
        gbz: str,
        minimizer: str,
        distance: str,
        zipcodes: str,
        ri: str,
        tags: str,
        gbwt_ri: str,
        t1: str,
        t2: str = "",
        threads: int = 8,
        max_multimaps: int = 1,
        batch_size: int = 256,
        output_timeout_s: float = 120.0,
    ) -> "PangenomeMiddleware":
        """Build both index configs from paths with correct defaults.

        Pass the LONG-READ minimizer + zipcodes (they must match the engine's
        long-read/BLAT preset). `max_multimaps=1` returns a single best
        alignment per read (recommended for the UI); set 0 to fall back to the
        engine's BLAT default of 100. No `--surject-target` is pre-indexed: the
        anchor pipeline surjects onto any haplotype at query time.
        """
        coord_paths = CoordinateIndexPaths(
            gbz_path=gbz,
            ri_path=ri,
            tags_path=tags,
            gbwt_ri_path=gbwt_ri,
            table1_path=t1,
            table2_path=t2,
        )
        giraffe_cfg = GiraffeServerConfig(
            vg_binary=vg_binary,
            gbz_path=gbz,
            minimizer_path=minimizer,
            distance_path=distance,
            zipcode_path=zipcodes,
            threads=threads,
            max_multimaps=max_multimaps,
            batch_size=batch_size,
            output_timeout_s=output_timeout_s,
            surject_target_paths=(),
        )
        return cls(coord_paths, giraffe_cfg)

    def set_call_timeout(self, seconds: float) -> None:
        """Set the per-giraffe-call output timeout (used to bound job wall time)."""
        self._giraffe.cfg.output_timeout_s = seconds

    # ------------------------------------------------------------------ #
    # Lifecycle
    # ------------------------------------------------------------------ #

    def wait_until_ready(self) -> None:
        """Block until giraffe-server has finished loading its indexes. The
        coordinate index is already loaded synchronously in __init__."""
        self._giraffe.wait_until_ready()

    def is_running(self) -> bool:
        return self._giraffe.is_running()

    def close(self) -> None:
        self._giraffe.stop()

    # ------------------------------------------------------------------ #
    # Coordinate translation (liftover_ext, in-process)
    # ------------------------------------------------------------------ #

    def translate(self, src_haplotype: str, start: int, end: int, tgt_haplotype: str):
        """Raw per-base source→target correspondences (liftover_ext.Index.translate).
        Prefer translate_intervals() for the folded target-interval form."""
        with self._coord_lock:
            return self._coord.translate(src_haplotype, start, end, tgt_haplotype)

    def translate_intervals(
        self, src: str, start: int, end: int, tgts: Sequence[str],
        timeout_ms: float = 5000.0,
        warnings: Optional[List[Dict[str, Any]]] = None,
        blocks: bool = False,
    ) -> List[Dict[str, Any]]:
        """Translate a source contig interval to one or more target haplotypes,
        returning folded target INTERVALS with strand.

        `src` is a full contig path (coordinates live on a contig); each entry of
        `tgts` is a target haplotype ("H#1") or contig ("H#1#chrX"). Each returned
        dict is {haplotype (full target contig), start, end (0-based half-open),
        strand ('+'/'-')}. Zero intervals for a target means the region does not
        exist there. Results across targets are concatenated; each interval is
        self-identifying via its contig name.
        """
        out: List[Dict[str, Any]] = []
        checked = getattr(self._coord, "translate_checked", None)
        with self._coord_lock:
            for tgt in tgts:
                if checked is not None and timeout_ms and timeout_ms > 0:
                    # Per-target deadline: one pathological haplotype is dropped
                    # with a warning instead of blocking the whole query.
                    run = checked(src, int(start), int(end), tgt, float(timeout_ms))
                    raw = run.intervals
                    if run.timed_out and warnings is not None:
                        warnings.append({
                            "haplotype": tgt,
                            "reason": "timeout",
                            "elapsed_ms": round(float(run.elapsed_ms), 1),
                            "message": (f"Translation to {tgt} exceeded "
                                        f"{int(timeout_ms)} ms and was stopped; "
                                        "results for it may be incomplete."),
                        })
                else:
                    raw = self._coord.translate(src, int(start), int(end), tgt)
                out.extend(_fold_to_blocks(raw) if blocks else _fold_to_intervals(raw))
        return out

    def translatable_haplotypes_scored(self, src: str, start: int, end: int,
                                       min_coverage: float = 0.0,
                                       max_nodes: int = 0) -> List[Dict[str, Any]]:
        """Haplotypes a source interval can reach, each with a 0-100 coverage
        score. Works without Table 2. Returns [] on an extension that predates
        the call, so callers can treat it as optional."""
        with self._coord_lock:
            fn = getattr(self._coord, "translatable_haplotypes_scored", None)
            if fn is None:
                return []
            rows = fn(src, int(start), int(end), float(min_coverage), int(max_nodes))
        return [{"haplotype": r.haplotype,
                 "coverage": round(float(r.coverage), 2),
                 "covered_bp": int(r.covered_bp)} for r in rows]

    def translatable_haplotypes(self, src: str, start: int, end: int) -> List[str]:
        """Discovery: the target haplotypes a source contig interval CAN translate
        to (names only, no coordinates).

        Uses the native liftover_ext.Index.translatable_haplotypes (a cheap
        Table-2 overlap check) when the loaded extension provides it. Falls back
        to probing every haplotype with translate() when it doesn't — correct but
        heavier — so this works before the extension is rebuilt and upgrades to
        the fast path automatically once it is.
        """
        with self._coord_lock:
            native = getattr(self._coord, "translatable_haplotypes", None)
            if native is not None:
                return list(native(src, int(start), int(end)))
            # Fallback: a haplotype is reachable iff translate() yields anything.
            names = self._coord.get_haplotype_names()
            out = [h for h in names
                   if self._coord.translate(src, int(start), int(end), h)]
            out.sort()
            return out

    def haplotype_coverage(self, gaf: str, min_coverage: float = 0.0,
                           include_zero: bool = False) -> List[Dict[str, Any]]:
        """Score every haplotype by how much of one alignment it accounts for.

        Returns [{haplotype, coverage, covered_bp}, ...] sorted by descending
        coverage, where `coverage` is a 0-100 percentage of the alignment's
        aligned bases lying on nodes that haplotype also visits.

        This is a graded companion to the engine's "carried by" list: that list
        holds only haplotypes threading the read's exact allele path, so a
        haplotype differing at a single variant is absent from it entirely,
        while here it scores near 100.

        `min_coverage=0` (the default) reports every haplotype sharing any node,
        however partial; `include_zero` additionally lists haplotypes sharing no
        node at all, scored 0, giving a complete table of every haplotype.

        Returns [] if the loaded extension predates this call, so callers can
        treat the field as optional rather than version-gating.
        """
        with self._coord_lock:
            fn = getattr(self._coord, "haplotype_coverage", None)
            if fn is None:
                return []
            try:
                rows = fn(gaf, float(min_coverage), bool(include_zero))
            except TypeError:
                # Extension built before include_zero existed.
                rows = fn(gaf, float(min_coverage))
        return [{"haplotype": r.haplotype,
                 "coverage": round(float(r.coverage), 2),
                 "covered_bp": int(r.covered_bp)} for r in rows]

    def get_haplotype_names(self) -> List[str]:
        """All haplotype/path names known to the coordinate index."""
        with self._coord_lock:
            return self._coord.get_haplotype_names()

    # ------------------------------------------------------------------ #
    # Mapping (giraffe-server)
    # ------------------------------------------------------------------ #

    def map_reads(
        self,
        reads: Iterable[FastqRead],
        surject_target: Optional[str] = None,
    ) -> List[List[str]]:
        """Map reads to the graph. Each result GAF line already carries the
        haplotype tags (hp/hn/hb/hf/hc/hl/hr/ht/hq/hv/hs/hm) the engine emits.
        Passing `surject_target` uses the engine's single-call surjection
        (requires the target pre-indexed via --surject-target); for arbitrary
        haplotypes use the anchor pipeline below instead."""
        return self._giraffe.map_reads(reads, surject_target=surject_target)

    def map_sequences(
        self,
        sequences: Iterable[str],
        surject_target: Optional[str] = None,
    ) -> List[List[str]]:
        return self._giraffe.map_sequences(sequences, surject_target=surject_target)

    # ------------------------------------------------------------------ #
    # Anchor-driven surjection (the fast path — surject onto ANY haplotype)
    # ------------------------------------------------------------------ #

    def build_surject_anchors(self, graph_alignment_gaf: str, target: str) -> AnchorBuild:
        """Build surjection anchors for one graph-alignment GAF line onto
        `target`. Returns (anchors, target_path_length, status). Uses
        build_surject_anchors_full when available so the cached target path
        length is returned (lets the server skip an O(path) re-walk)."""
        if hasattr(self._coord, "build_surject_anchors_full"):
            res = self._coord.build_surject_anchors_full(graph_alignment_gaf, target)
            return list(res.anchors), int(res.target_path_length), str(getattr(res, "status", "ok"))
        anchors = self._coord.build_surject_anchors(graph_alignment_gaf, target)
        return list(anchors), 0, "ok"

    def surject_with_anchors(
        self,
        graph_alignment_gaf: str,
        anchors,
        target: str,
        target_path_length: int = 0,
    ) -> List[str]:
        """Surject one graph alignment onto `target` using pre-computed anchors
        (SURJECT_WITH_ANCHORS on giraffe-server). Returns the GAF line(s) with
        surjection tags (sj/sn/sp/sr/ss/sm/sc/an/ap) appended."""
        return self._giraffe.surject_with_anchors(
            graph_alignment_gaf, anchors, target, target_path_length=target_path_length
        )

    def surject_sequence(
        self,
        sequence: str,
        target: str,
        name: Optional[str] = None,
        quality: Optional[str] = None,
    ) -> List[str]:
        """Full 3-step anchor pipeline for ONE sequence onto `target`:
            1. map to the graph (graph-only)       -> giraffe-server
            2. build anchors per alignment          -> liftover_ext
            3. surject each with those anchors      -> giraffe-server
        Returns one GAF line per (multimapped) alignment, each carrying both the
        haplotype tags (from step 1) and the surjection tags (from step 3).
        An alignment with no anchors is returned unchanged with sj:Z:no_anchors.
        """
        seq = sequence.strip().upper()
        if not seq:
            return []
        nm = name or f"q_{seq[:8]}"
        qual = quality if quality else ("I" * len(seq))

        # Step 1: graph-only mapping (a 3-tuple -> no surjection target).
        out = self._giraffe.map_reads([(nm, seq, qual)])
        graph_alignments = out[0] if out else []

        lines: List[str] = []
        for gaf in graph_alignments:
            anchors, path_len, _status = self.build_surject_anchors(gaf, target)
            if not anchors:
                lines.append(gaf + "\tsj:Z:no_anchors")
                continue
            lines.extend(
                self._giraffe.surject_with_anchors(
                    gaf, anchors, target, target_path_length=path_len
                )
            )
        return lines

    def surject_sequences(
        self,
        sequences: Sequence[str],
        target: str,
        names: Optional[Sequence[str]] = None,
    ) -> List[List[str]]:
        """Batched full anchor pipeline for many sequences onto `target`.

        Mapping (step 1) is done in a single batched giraffe-server call for
        throughput; anchor building (step 2) and surjection (step 3) run per
        alignment. Returns one list of GAF lines per input sequence, in input
        order. Sequences must be non-empty (validate/cap upstream — e.g. the UI
        limits to 50 sequences per request).
        """
        reads: List[FastqRead] = []
        for i, s in enumerate(sequences):
            su = s.strip().upper()
            if not su:
                raise ValueError(f"sequence at index {i} is empty")
            nm = names[i] if names is not None else f"q{i}_{su[:8]}"
            reads.append((nm, su, "I" * len(su)))
        if not reads:
            return []

        # Step 1: one batched graph-only mapping call (1:1 with reads).
        mapped = self._giraffe.map_reads(reads)

        results: List[List[str]] = []
        for graph_alignments in mapped:
            per_read: List[str] = []
            for gaf in graph_alignments:
                anchors, path_len, _status = self.build_surject_anchors(gaf, target)
                if not anchors:
                    per_read.append(gaf + "\tsj:Z:no_anchors")
                    continue
                per_read.extend(
                    self._giraffe.surject_with_anchors(
                        gaf, anchors, target, target_path_length=path_len
                    )
                )
            results.append(per_read)
        return results
