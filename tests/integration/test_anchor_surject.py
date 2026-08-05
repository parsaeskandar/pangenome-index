#!/usr/bin/env python3
"""
test_anchor_surject.py — comprehensive comparison of anchor-based surjection
against ordinary (ReferencePathOverlay) surjection on a running pangenome
server.

For each sampled read it:
  1. maps + surjects via giraffe-server's built-in path (ReferencePathOverlay)
     in a single map_reads call → "regular" surjection,
  2. takes that same graph alignment, builds anchors with
     liftover_ext.Index.build_surject_anchors_full, and surjects via the
     SURJECT_WITH_ANCHORS stdin command → "anchor" surjection,
  3. compares the two surjections field-by-field
     (sj status, sn path, sp position, sr strand, ss score, sc CIGAR).

Reads are sampled from a newline-separated sequences file (one haplotype
sequence per line, no headers/names). For the regular path to actually
surject (not report INCOMPATIBLE), reads should come from the target
haplotype — so the most informative run uses a --sequences file that is the
target haplotype's own sequence, and --surject-target naming that path. In
that setup --ground-truth also checks each surjected position against the
known sample offset.

Example (chrM):
  WORK=/path/to/chrM
  python3 test_anchor_surject.py \
    --vg ../vg --gbz $WORK/graph.gbz \
    --minimizer $WORK/index.shortread.withzip.min \
    --dist $WORK/index.dist --zipcodes $WORK/index.shortread.zipcodes \
    --ri $WORK/rlbwt_rindex.ri --tags $WORK/sampled.tags \
    --gbwt-ri $WORK/gbwt_fastlocate.ri --t1 $WORK/output.t1 --t2 $WORK/output.t2 \
    --sequences $WORK/HG002_2_chrM.seq.txt \
    --surject-target 'HG002#2#2chrM#0' \
    --num-reads 200 --read-length 150 --ground-truth
"""
from __future__ import annotations

# Repo root (holds liftover_ext.so + the middleware package) on sys.path, so
# this script runs from anywhere after being moved under tests/integration/.
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import argparse
import importlib.util
import os
import random
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional

# Load middleware/giraffe_server_middleware.py directly rather than via
# `from middleware import ...`, because the package __init__ eagerly imports
# pangenome_middleware → liftover_ext (the compiled extension). Loading the
# file directly keeps --help / argument parsing working without the native
# library present; liftover_ext is imported lazily in main() when actually
# running. (giraffe_server_middleware.py itself has no native dependency.)
def _load_giraffe_middleware():
    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    path = os.path.join(repo_root, "middleware", "giraffe_server_middleware.py")
    spec = importlib.util.spec_from_file_location("giraffe_server_middleware", path)
    mod = importlib.util.module_from_spec(spec)
    # Register before exec so @dataclass (which looks up cls.__module__ in
    # sys.modules) resolves correctly.
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


_gm = _load_giraffe_middleware()
GiraffeServerConfig = _gm.GiraffeServerConfig
GiraffeServerMiddleware = _gm.GiraffeServerMiddleware

# Surjection-related GAF optional tags appended by both code paths. Stripped
# from a regular GAF to recover the clean graph alignment, and parsed to
# compare the two surjections. (an/ap/se are anchor-path extras.)
_SURJ_TAG_KEYS = {"sj", "sn", "sp", "sr", "ss", "sm", "sc", "an", "ap", "se"}

_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_RC)[::-1]


# ──────────────────────────────────────────────────────────────────────────
# Arguments
# ──────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compare anchor-based vs ordinary surjection on a pangenome server.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # giraffe-server indexes
    p.add_argument("--vg", required=True, help="Path to vg binary (with giraffe-server + SURJECT_WITH_ANCHORS)")
    p.add_argument("--gbz", required=True)
    p.add_argument("--minimizer", required=True)
    p.add_argument("--dist", required=True)
    p.add_argument("--zipcodes", required=True)
    # coordinate (liftover) indexes
    p.add_argument("--ri", required=True, help="RLBWT r-index (.ri)")
    p.add_argument("--tags", required=True, help="Sampled tag array (.tags)")
    p.add_argument("--gbwt-ri", required=True, help="GBWT FastLocate r-index (.ri)")
    p.add_argument("--t1", required=True, help="Translation table 1 (.t1)")
    p.add_argument("--t2", required=True, help="Translation table 2 (.t2)")
    # reads
    p.add_argument("--sequences", required=True,
                   help="Newline-separated haplotype sequences (one per line, no headers).")
    p.add_argument("--surject-target", required=True,
                   help="Path/haplotype name to surject onto (both methods use this).")
    # sampling
    p.add_argument("--num-reads", type=int, default=200)
    p.add_argument("--read-length", type=int, default=150)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--revcomp-fraction", type=float, default=0.0,
                   help="Fraction of reads to sample as reverse complement [0..1].")
    p.add_argument("--max-n-fraction", type=float, default=0.1,
                   help="Skip windows with more than this fraction of N bases.")
    # source-file memory guards (whole-genome sequence files can be huge)
    p.add_argument("--max-source-seqs", type=int, default=1000,
                   help="Read at most this many sequence lines from --sequences.")
    p.add_argument("--max-total-bytes", type=int, default=200_000_000,
                   help="Stop reading --sequences after this many sequence bytes.")
    # server config
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--batch-size", type=int, default=64)
    p.add_argument("--max-multimaps", type=int, default=0)
    # comparison / reporting
    p.add_argument("--ground-truth", action="store_true",
                   help="Treat the sample start offset as the expected surjected "
                        "position (reads must be sampled from the target's own sequence).")
    p.add_argument("--pos-tolerance", type=int, default=0,
                   help="Max |position difference| still counted as a position match.")
    p.add_argument("--min-agreement", type=float, default=0.95,
                   help="Minimum position-agreement rate among both-ok primaries to PASS.")
    p.add_argument("--out-tsv", default="",
                   help="Optional path to write per-read results as TSV.")
    p.add_argument("--max-mismatches-shown", type=int, default=20)
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


# ──────────────────────────────────────────────────────────────────────────
# Sequence loading + read sampling
# ──────────────────────────────────────────────────────────────────────────

def load_sequences(path: str, max_seqs: int, max_total_bytes: int) -> List[str]:
    seqs: List[str] = []
    total = 0
    with open(path) as fh:
        for line in fh:
            s = line.strip().upper()
            if not s:
                continue
            seqs.append(s)
            total += len(s)
            if len(seqs) >= max_seqs or total >= max_total_bytes:
                break
    return seqs


@dataclass
class SampledRead:
    name: str
    seq: str
    source_idx: int
    source_start: int
    is_revcomp: bool


def sample_reads(seqs: List[str], n: int, read_len: int,
                 revcomp_fraction: float, max_n_fraction: float,
                 rng: random.Random) -> List[SampledRead]:
    # Weight sequence choice by length so longer contigs get proportionally
    # more reads (uniform coverage over bases, not over sequences).
    eligible = [(i, len(s)) for i, s in enumerate(seqs) if len(s) > 0]
    if not eligible:
        return []
    weights = [w for _, w in eligible]
    reads: List[SampledRead] = []
    attempts = 0
    max_attempts = n * 20
    while len(reads) < n and attempts < max_attempts:
        attempts += 1
        (src_idx, src_len) = rng.choices(eligible, weights=weights, k=1)[0]
        s = seqs[src_idx]
        if src_len <= read_len:
            start, window = 0, s
        else:
            start = rng.randint(0, src_len - read_len)
            window = s[start:start + read_len]
        # Reject windows that are mostly N (won't map cleanly).
        n_count = window.count("N")
        if n_count > max_n_fraction * len(window):
            continue
        is_rc = rng.random() < revcomp_fraction
        seq = revcomp(window) if is_rc else window
        reads.append(SampledRead(
            name=f"r{len(reads)}_{src_idx}_{start}{'_rc' if is_rc else ''}",
            seq=seq, source_idx=src_idx, source_start=start, is_revcomp=is_rc,
        ))
    return reads


# ──────────────────────────────────────────────────────────────────────────
# GAF tag helpers
# ──────────────────────────────────────────────────────────────────────────

def parse_surj_tags(gaf_line: str) -> Dict[str, str]:
    """Return {tag_key: value} for the surjection tags present on a GAF line."""
    out: Dict[str, str] = {}
    for tok in gaf_line.split("\t"):
        # tag form: KEY:TYPE:VALUE
        if len(tok) >= 5 and tok[2] == ":" and tok[4] == ":":
            key = tok[:2]
            if key in _SURJ_TAG_KEYS:
                out[key] = tok[5:]
    return out


def strip_surj_tags(gaf_line: str) -> str:
    """Drop surjection-related optional tags, recovering the graph alignment."""
    kept = []
    for tok in gaf_line.split("\t"):
        if len(tok) >= 5 and tok[2] == ":" and tok[4] == ":" and tok[:2] in _SURJ_TAG_KEYS:
            continue
        kept.append(tok)
    return "\t".join(kept)


@dataclass
class Surjection:
    status: str = "missing"
    path_name: str = ""
    position: Optional[int] = None
    strand: Optional[int] = None
    score: Optional[int] = None
    mapq: Optional[int] = None
    cigar: str = ""

    @property
    def ok(self) -> bool:
        return self.status == "ok"

    @classmethod
    def from_gaf(cls, gaf_line: str) -> "Surjection":
        t = parse_surj_tags(gaf_line)
        def _int(k):
            v = t.get(k)
            try:
                return int(v) if v is not None else None
            except ValueError:
                return None
        return cls(
            status=t.get("sj", "missing"),
            path_name=t.get("sn", ""),
            position=_int("sp"),
            strand=_int("sr"),
            score=_int("ss"),
            mapq=_int("sm"),
            cigar=t.get("sc", ""),
        )


# ──────────────────────────────────────────────────────────────────────────
# Per-read comparison
# ──────────────────────────────────────────────────────────────────────────

@dataclass
class ReadResult:
    read: SampledRead
    mapped: bool = False
    regular: Surjection = field(default_factory=Surjection)
    anchor: Surjection = field(default_factory=Surjection)
    n_anchors: int = 0
    anchor_build_status: str = ""
    error: str = ""

    @property
    def both_ok(self) -> bool:
        return self.regular.ok and self.anchor.ok


def run_one(mw: GiraffeServerMiddleware, coord_index, read: SampledRead,
            target: str) -> ReadResult:
    rr = ReadResult(read=read)
    qual = "I" * len(read.seq)

    # Regular: map + surject in one call. The primary (best) alignment is [0].
    try:
        reg_out = mw.map_reads([(read.name, read.seq, qual, target)])
    except Exception as exc:  # noqa: BLE001
        rr.error = f"regular map_reads: {exc}"
        return rr
    reg_alignments = reg_out[0] if reg_out else []
    if not reg_alignments:
        return rr  # unmapped
    rr.mapped = True
    primary_reg_gaf = reg_alignments[0]
    rr.regular = Surjection.from_gaf(primary_reg_gaf)

    # Recover the clean graph alignment from the regular primary, then anchor.
    graph_gaf = strip_surj_tags(primary_reg_gaf)
    try:
        build = coord_index.build_surject_anchors_full(graph_gaf, target)
        rr.anchor_build_status = build.status
        anchors = list(build.anchors)
        rr.n_anchors = len(anchors)
        path_len = build.target_path_length
    except Exception as exc:  # noqa: BLE001
        rr.error = f"build_surject_anchors_full: {exc}"
        return rr

    if not anchors:
        # No anchors → anchor surjection can't run; leave anchor.status as the
        # build status so the report distinguishes this from a surject failure.
        rr.anchor = Surjection(status=f"no_anchors({build.status})")
        return rr

    try:
        anc_lines = mw.surject_with_anchors(
            graph_gaf, anchors, target,
            target_path_length=path_len, read_name=read.name + "_anc",
        )
    except Exception as exc:  # noqa: BLE001
        rr.error = f"surject_with_anchors: {exc}"
        return rr
    if anc_lines:
        rr.anchor = Surjection.from_gaf(anc_lines[0])
    return rr


# ──────────────────────────────────────────────────────────────────────────
# Aggregation + reporting
# ──────────────────────────────────────────────────────────────────────────

def pos_match(a: Optional[int], b: Optional[int], tol: int) -> bool:
    return a is not None and b is not None and abs(a - b) <= tol


def report(results: List[ReadResult], args: argparse.Namespace) -> int:
    n_total = len(results)
    n_mapped = sum(1 for r in results if r.mapped)
    n_errors = sum(1 for r in results if r.error)
    both_ok = [r for r in results if r.both_ok]
    reg_ok_only = [r for r in results if r.regular.ok and not r.anchor.ok]
    anc_ok_only = [r for r in results if r.anchor.ok and not r.regular.ok]
    neither_ok = [r for r in results if r.mapped and not r.regular.ok and not r.anchor.ok]

    print("=" * 70)
    print("Anchor vs Regular surjection comparison")
    print("=" * 70)
    print(f"target            : {args.surject_target}")
    print(f"reads tested      : {n_total}  "
          f"(read_length={args.read_length}, revcomp={args.revcomp_fraction:.0%})")
    print(f"mapped (>=1 aln)  : {n_mapped} / {n_total}")
    if n_errors:
        print(f"per-read errors   : {n_errors}  (see details below)")
    print()
    print("Surjection outcomes (primary alignment):")
    print(f"  both ok                 : {len(both_ok)}")
    print(f"  regular ok, anchor not  : {len(reg_ok_only)}")
    print(f"  anchor ok, regular not  : {len(anc_ok_only)}")
    print(f"  neither ok              : {len(neither_ok)}")
    print()

    # Agreement among both-ok primaries.
    sn_match = sum(1 for r in both_ok if r.regular.path_name == r.anchor.path_name)
    sp_match = sum(1 for r in both_ok if pos_match(r.regular.position, r.anchor.position, args.pos_tolerance))
    sr_match = sum(1 for r in both_ok if r.regular.strand == r.anchor.strand)
    ss_match = sum(1 for r in both_ok if r.regular.score == r.anchor.score)
    sc_match = sum(1 for r in both_ok if r.regular.cigar == r.anchor.cigar)
    n_both = len(both_ok)

    def pct(x: int) -> str:
        return f"{(100.0 * x / n_both):.1f}%" if n_both else "n/a"

    print(f"Agreement among both-ok primaries ({n_both}):")
    print(f"  path name (sn) match : {sn_match} / {n_both}  ({pct(sn_match)})")
    print(f"  position  (sp) match : {sp_match} / {n_both}  ({pct(sp_match)})  (tol={args.pos_tolerance})")
    print(f"  strand    (sr) match : {sr_match} / {n_both}  ({pct(sr_match)})")
    print(f"  score     (ss) match : {ss_match} / {n_both}  ({pct(ss_match)})")
    print(f"  CIGAR     (sc) match : {sc_match} / {n_both}  ({pct(sc_match)})")

    # |sp| difference histogram.
    buckets = {"0": 0, "1-5": 0, "6-50": 0, "51+": 0, "na": 0}
    for r in both_ok:
        if r.regular.position is None or r.anchor.position is None:
            buckets["na"] += 1
            continue
        d = abs(r.regular.position - r.anchor.position)
        if d == 0:
            buckets["0"] += 1
        elif d <= 5:
            buckets["1-5"] += 1
        elif d <= 50:
            buckets["6-50"] += 1
        else:
            buckets["51+"] += 1
    print("  |sp diff| histogram  : " +
          "  ".join(f"{k}:{v}" for k, v in buckets.items()))
    print()

    # Ground-truth check (expects reads sampled from the target's own sequence).
    gt_lines = []
    if args.ground_truth:
        reg_gt = anc_gt = gt_n = 0
        for r in both_ok:
            exp = r.read.source_start
            gt_n += 1
            if pos_match(r.regular.position, exp, args.pos_tolerance):
                reg_gt += 1
            if pos_match(r.anchor.position, exp, args.pos_tolerance):
                anc_gt += 1
        gt_lines.append(f"Ground truth (expected pos = sample offset), among {gt_n} both-ok:")
        gt_lines.append(f"  regular within tol : {reg_gt} / {gt_n}")
        gt_lines.append(f"  anchor  within tol : {anc_gt} / {gt_n}")
        gt_lines.append("  NOTE: meaningful only if --sequences IS the target's sequence "
                        "and the target is a single (non-subpath-offset) path.")
        print("\n".join(gt_lines))
        print()

    # Mismatch listing (both-ok but disagree on position or CIGAR).
    mismatches = [r for r in both_ok
                  if not pos_match(r.regular.position, r.anchor.position, args.pos_tolerance)
                  or r.regular.cigar != r.anchor.cigar]
    if mismatches:
        print(f"Mismatches among both-ok ({len(mismatches)}; showing up to {args.max_mismatches_shown}):")
        for r in mismatches[:args.max_mismatches_shown]:
            dpos = ("na" if r.regular.position is None or r.anchor.position is None
                    else str(r.anchor.position - r.regular.position))
            print(f"  [{r.read.name}] "
                  f"sn reg={r.regular.path_name} anc={r.anchor.path_name} | "
                  f"sp reg={r.regular.position} anc={r.anchor.position} (Δ{dpos}) | "
                  f"sr reg={r.regular.strand} anc={r.anchor.strand} | "
                  f"ss reg={r.regular.score} anc={r.anchor.score} | "
                  f"sc reg={r.regular.cigar} anc={r.anchor.cigar} | "
                  f"n_anchors={r.n_anchors}")
        print()

    # Cases where the methods disagree on whether surjection succeeded.
    if reg_ok_only or anc_ok_only:
        print(f"Status disagreements (showing up to {args.max_mismatches_shown}):")
        for r in (reg_ok_only + anc_ok_only)[:args.max_mismatches_shown]:
            print(f"  [{r.read.name}] regular={r.regular.status} anchor={r.anchor.status} "
                  f"(build={r.anchor_build_status}, n_anchors={r.n_anchors})"
                  + (f" err={r.error}" if r.error else ""))
        print()

    if n_errors and args.verbose:
        print("Per-read errors:")
        for r in results:
            if r.error:
                print(f"  [{r.read.name}] {r.error}")
        print()

    # Optional TSV dump.
    if args.out_tsv:
        _write_tsv(results, args.out_tsv, args.pos_tolerance)
        print(f"Per-read results written to {args.out_tsv}")
        print()

    # PASS/FAIL: gated on position agreement among reads both methods placed.
    # Status disagreements (one method surjects, the other doesn't) are
    # reported but NOT failed on — they can be legitimate (e.g. the regular
    # path's haplotype-compatibility check rejects a read the anchor path can
    # still place). The point of this test is: when both place a read, do
    # they place it in the same spot?
    pass_pos = (sp_match / n_both) >= args.min_agreement if n_both else False
    ok = bool(n_both) and pass_pos
    verdict = "PASS" if ok else "FAIL"
    reason = ""
    if not n_both:
        reason = " (no reads where both methods surjected — check target/sequences)"
    elif not pass_pos:
        reason = f" (position agreement {pct(sp_match)} < required {args.min_agreement:.0%})"
    print(f"{verdict}{reason}")
    return 0 if ok else 1


def _write_tsv(results: List[ReadResult], path: str, tol: int) -> None:
    cols = [
        "read", "source_idx", "source_start", "revcomp", "mapped",
        "reg_status", "anc_status", "anchor_build_status", "n_anchors",
        "reg_sn", "anc_sn", "reg_sp", "anc_sp", "sp_diff",
        "reg_sr", "anc_sr", "reg_ss", "anc_ss",
        "reg_sc", "anc_sc", "pos_match", "cigar_match", "error",
    ]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in results:
            sp_diff = ("" if r.regular.position is None or r.anchor.position is None
                       else str(r.anchor.position - r.regular.position))
            row = [
                r.read.name, r.read.source_idx, r.read.source_start,
                int(r.read.is_revcomp), int(r.mapped),
                r.regular.status, r.anchor.status, r.anchor_build_status, r.n_anchors,
                r.regular.path_name, r.anchor.path_name,
                r.regular.position, r.anchor.position, sp_diff,
                r.regular.strand, r.anchor.strand,
                r.regular.score, r.anchor.score,
                r.regular.cigar, r.anchor.cigar,
                int(pos_match(r.regular.position, r.anchor.position, tol)),
                int(r.regular.cigar == r.anchor.cigar),
                r.error,
            ]
            fh.write("\t".join("" if c is None else str(c) for c in row) + "\n")


# ──────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────

def main() -> int:
    args = parse_args()
    rng = random.Random(args.seed)

    print("Loading sequences…", file=sys.stderr)
    seqs = load_sequences(args.sequences, args.max_source_seqs, args.max_total_bytes)
    if not seqs:
        print(f"ERROR: no sequences read from {args.sequences}", file=sys.stderr)
        return 2
    total_bp = sum(len(s) for s in seqs)
    print(f"  {len(seqs)} sequence(s), {total_bp} bp total", file=sys.stderr)

    reads = sample_reads(seqs, args.num_reads, args.read_length,
                         args.revcomp_fraction, args.max_n_fraction, rng)
    if not reads:
        print("ERROR: could not sample any reads (sequences too short / too many Ns?)",
              file=sys.stderr)
        return 2
    print(f"  sampled {len(reads)} reads", file=sys.stderr)

    cfg = GiraffeServerConfig(
        vg_binary=args.vg,
        gbz_path=args.gbz,
        minimizer_path=args.minimizer,
        distance_path=args.dist,
        zipcode_path=args.zipcodes,
        threads=args.threads,
        max_multimaps=args.max_multimaps,
        batch_size=args.batch_size,
        surject_target_paths=[args.surject_target],
    )
    mw = GiraffeServerMiddleware(cfg)
    mw.start()

    try:
        print("Waiting for giraffe-server to finish loading indexes…", file=sys.stderr)
        mw.wait_until_ready()
        print("giraffe-server is ready.", file=sys.stderr)

        import liftover_ext
        coord_index = liftover_ext.Index()
        print("Loading coordinate index…", file=sys.stderr)
        coord_index.load(args.gbz, args.ri, args.tags, args.gbwt_ri, args.t1, args.t2)
        print("Coordinate index loaded.", file=sys.stderr)

        results: List[ReadResult] = []
        for i, read in enumerate(reads):
            rr = run_one(mw, coord_index, read, args.surject_target)
            results.append(rr)
            if args.verbose:
                print(f"[{i}] {read.name} mapped={rr.mapped} "
                      f"reg={rr.regular.status}/{rr.regular.position} "
                      f"anc={rr.anchor.status}/{rr.anchor.position} "
                      f"n_anchors={rr.n_anchors}", file=sys.stderr)
            elif (i + 1) % 25 == 0:
                print(f"  …{i + 1}/{len(reads)} reads processed", file=sys.stderr)
    finally:
        mw.stop()

    return report(results, args)


if __name__ == "__main__":
    raise SystemExit(main())
