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

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple

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
    table2_path: str


# Result of building anchors for one graph alignment onto a target haplotype:
# (anchors, target_path_length, status).
AnchorBuild = Tuple[list, int, str]


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
        t2: str,
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
        """Translate an interval on one haplotype/contig to another."""
        return self._coord.translate(src_haplotype, start, end, tgt_haplotype)

    def get_haplotype_names(self) -> List[str]:
        """All haplotype/path names known to the coordinate index."""
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
