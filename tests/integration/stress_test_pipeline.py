#!/usr/bin/env python3
"""
stress_test_pipeline.py — end-to-end stress / demo test for the pangenome
server pipeline:  giraffe mapping  →  anchor build  →  surjection onto ANY
haplotype.

WHY THIS IS CONVINCING (the unfakeable test)
  Every read is sampled from a KNOWN base offset on a KNOWN haplotype. We then
  optionally reverse-complement it and inject sequencing errors, map it to the
  pangenome with no target, and surject it back onto its source haplotype using
  pre-computed anchors. A correct pipeline MUST place it at the exact original
  base — there is no way to fake this: the ground truth is the sampling offset.

WHAT IT DEMONSTRATES
  • correctness at scale          — thousands of reads, % landing on the exact base
  • the unique capability         — surjection onto NON-reference haplotypes,
                                    which `vg giraffe --ref-paths` cannot target
  • robustness                    — forward AND reverse strand; clean reads AND
                                    reads with substitution errors; multiple
                                    read lengths
  • coverage                      — many distinct haplotypes and chromosomes
  • performance                   — per-stage timing and throughput

Reuses the proven pipeline from benchmark_surject_anchors.py (map → build
anchors → surject, with the same ground-truth position check) and adds error
injection, a varied read mix, category breakdowns, and a demo-quality report.

Example:
  python3 stress_test_pipeline.py \
    --vg ../1.server/Giraffe_server/bin/vg --gbz $W/graph.gbz \
    --minimizer $W/...longread.withzip.min --dist $W/...dist \
    --zipcodes $W/...longread.zipcodes \
    --ri $W/rlbwt_rindex.ri --tags $W/sampled.tags --gbwt-ri $W/fastlocate.ri \
    --t1 $W/table1.t1 --t2 $W/table2.t2 \
    --sequences $W/whole_genome_info --names $W/path_names.txt --both-orientations \
    --read-lengths 150,1000,5000 --reads-per-length 400 \
    --revcomp-fraction 0.5 --error-fraction 0.5 --error-rate 0.01 \
    --reads-cache stress.reads.tsv --out-tsv stress.results.tsv
"""
from __future__ import annotations

# Repo root (holds liftover_ext.so + the middleware package) on sys.path, so
# this script runs from anywhere after being moved under tests/integration/.
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import argparse
import os
import random
import re
import sys
from typing import Dict, List, Optional, Tuple

from benchmark_surject_anchors import (
    BenchResult,
    GiraffeServerConfig,
    GiraffeServerMiddleware,
    Read,
    base_target,
    load_from_fasta,
    load_from_pair,
    load_reads,
    revcomp,
    run_one,
    save_reads,
    _fmt,
    timing_stats,
)


# ──────────────────────────────────────────────────────────────────────────
# Arguments
# ──────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="End-to-end stress/demo test for giraffe mapping + surjection onto any haplotype.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # giraffe-server indexes
    p.add_argument("--vg", required=True)
    p.add_argument("--gbz", required=True)
    p.add_argument("--minimizer", required=True)
    p.add_argument("--dist", required=True)
    p.add_argument("--zipcodes", required=True)
    # coordinate (liftover) indexes
    p.add_argument("--ri", required=True)
    p.add_argument("--tags", required=True)
    p.add_argument("--gbwt-ri", required=True)
    p.add_argument("--t1", required=True)
    p.add_argument("--t2", required=True)
    # name↔sequence input
    p.add_argument("--fasta", help="vg paths --extract-fasta output.")
    p.add_argument("--sequences", help="gbz_extract sequences file.")
    p.add_argument("--names", help="vg paths -L output (path names).")
    p.add_argument("--both-orientations", action="store_true",
                   help="sequences file built with gbz_extract -b (2 lines/path).")
    p.add_argument("--name-filter", default="",
                   help="Only draw reads from paths whose name matches this regex "
                        "(e.g. 'chr' or '#CM0' to favor chromosomes over scaffolds).")
    p.add_argument("--max-haplotypes", type=int, default=40)
    p.add_argument("--max-seq-len", type=int, default=20_000_000)
    p.add_argument("--target-name-fields", type=int, default=3)
    # read mix
    p.add_argument("--read-lengths", default="150,1000,5000",
                   help="Comma-separated read lengths to test (each gets --reads-per-length reads).")
    p.add_argument("--reads-per-length", type=int, default=300)
    p.add_argument("--revcomp-fraction", type=float, default=0.5,
                   help="Fraction of reads sampled on the reverse strand.")
    p.add_argument("--error-fraction", type=float, default=0.5,
                   help="Fraction of reads that get substitution errors injected.")
    p.add_argument("--error-rate", type=float, default=0.01,
                   help="Per-base substitution rate for error reads (e.g. 0.01 = 1%% per base).")
    p.add_argument("--max-n-fraction", type=float, default=0.05)
    p.add_argument("--reference-samples", default="CHM13,GRCh38",
                   help="Comma-separated sample names treated as 'reference' for the "
                        "reference vs non-reference breakdown.")
    p.add_argument("--seed", type=int, default=0)
    # server config
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--batch-size", type=int, default=64)
    p.add_argument("--max-multimaps", type=int, default=1,
                   help="Max alignments per read. Default 1: only the primary is surjected, "
                        "so reporting many just ships huge GAF lines and risks framing desync.")
    p.add_argument("--timeout", type=float, default=120.0)
    # correctness / reporting
    p.add_argument("--pos-tolerance", type=int, default=0,
                   help="Max |surjected pos - sample offset| counted as exact. 0 = exact base.")
    p.add_argument("--min-exact", type=float, default=0.98,
                   help="Min exact-position rate (among surjected) to PASS.")
    p.add_argument("--reads-cache", default="",
                   help="Cache the generated reads here; reused on later runs to skip "
                        "the (slow) sequence-file streaming. Delete to regenerate.")
    p.add_argument("--out-tsv", default="")
    p.add_argument("--warmup", type=int, default=10)
    p.add_argument("--max-mismatches-shown", type=int, default=15)
    p.add_argument("--progress-every", type=int, default=100)
    return p.parse_args()


# ──────────────────────────────────────────────────────────────────────────
# Read generation with error injection
# ──────────────────────────────────────────────────────────────────────────

_BASES = "ACGT"


def inject_substitutions(seq: str, rate: float, rng: random.Random) -> Tuple[str, int]:
    """Substitute bases at the given per-base rate (never to the same base).
    Substitutions don't shift coordinates, so the ground-truth offset is
    preserved exactly — the surjected position must still equal the offset."""
    if rate <= 0:
        return seq, 0
    out = list(seq)
    n = 0
    for i, c in enumerate(out):
        if c in _BASES and rng.random() < rate:
            out[i] = rng.choice([b for b in _BASES if b != c])
            n += 1
    return "".join(out), n


def generate_reads(named_seqs: List[Tuple[str, str]], lengths: List[int],
                   reads_per_length: int, revcomp_fraction: float,
                   error_fraction: float, error_rate: float,
                   max_n_fraction: float, target_fields: int,
                   rng: random.Random) -> List[Read]:
    """Build a varied read set: for each length, sample reads from known offsets,
    half reverse-complemented and a fraction given substitution errors. The read
    name encodes the truth (offset, length, strand, #errors) so category and
    ground truth survive caching."""
    eligible = [(nm, s) for nm, s in named_seqs if len(s) > 0]
    if not eligible:
        return []
    weights = [len(s) for _, s in eligible]
    reads: List[Read] = []
    idx = 0
    for L in lengths:
        made = 0
        attempts = 0
        cap = reads_per_length * 40
        while made < reads_per_length and attempts < cap:
            attempts += 1
            nm, s = rng.choices(eligible, weights=weights, k=1)[0]
            if len(s) < L:
                continue
            start = rng.randint(0, len(s) - L)
            window = s[start:start + L]
            if window.count("N") > max_n_fraction * L:
                continue
            is_rc = rng.random() < revcomp_fraction
            n_err = 0
            seq = window
            if rng.random() < error_fraction and error_rate > 0:
                seq, n_err = inject_substitutions(seq, error_rate, rng)
            if is_rc:
                seq = revcomp(seq)
            name = (f"s{idx}_{start}_{L}"
                    + ("_rc" if is_rc else "")
                    + (f"_e{n_err}" if n_err else ""))
            reads.append(Read(
                name=name, seq=seq, source_name=nm,
                target=base_target(nm, target_fields),
                offset=start, is_revcomp=is_rc,
            ))
            idx += 1
            made += 1
    rng.shuffle(reads)  # interleave lengths/strands so warmup isn't all one kind
    return reads


# ──────────────────────────────────────────────────────────────────────────
# Category derived from a Read (recoverable after caching)
# ──────────────────────────────────────────────────────────────────────────

class Category:
    __slots__ = ("length", "strand", "errored", "n_errors", "sample", "is_reference")

    def __init__(self, read: Read, reference_samples) -> None:
        self.length = len(read.seq)
        self.strand = "-" if read.is_revcomp else "+"
        m = re.search(r"_e(\d+)$", read.name)
        self.n_errors = int(m.group(1)) if m else 0
        self.errored = self.n_errors > 0
        self.sample = read.source_name.split("#")[0]
        self.is_reference = self.sample in reference_samples


# ──────────────────────────────────────────────────────────────────────────
# Reporting
# ──────────────────────────────────────────────────────────────────────────

def _pct(x: int, n: int) -> str:
    return f"{(100.0 * x / n):.1f}%" if n else "n/a"


def _breakdown(rows: List[Tuple[Read, Category, BenchResult]],
               key_fn, tol: int) -> List[Tuple[str, int, int]]:
    """Group surjected reads by key_fn(category) → (label, n_surjected, n_exact)."""
    groups: Dict[str, List[int]] = {}
    for _read, cat, rr in rows:
        if not rr.surj.ok:
            continue
        k = str(key_fn(cat))
        g = groups.setdefault(k, [0, 0])
        g[0] += 1
        if rr.is_pos_correct(tol):
            g[1] += 1
    # sort numerically when keys look like ints, else lexically
    def sort_key(item):
        try:
            return (0, int(item[0]))
        except ValueError:
            return (1, item[0])
    return [(k, v[0], v[1]) for k, v in sorted(groups.items(), key=sort_key)]


def _fmt_breakdown(groups: List[Tuple[str, int, int]]) -> str:
    return "   ".join(f"{k}: {_pct(ex, n)} ({ex}/{n})" for k, n, ex in groups)


def report(rows: List[Tuple[Read, Category, BenchResult]], args: argparse.Namespace,
           n_paths: int, ref_samples) -> int:
    tol = args.pos_tolerance
    total = len(rows)
    mapped = [t for t in rows if t[2].mapped]
    surjected = [t for t in rows if t[2].surj.ok]
    nm = len(mapped)
    ns = len(surjected)

    exact = sum(1 for _r, _c, rr in surjected if rr.is_pos_correct(tol))
    strand_ok = sum(1 for _r, _c, rr in surjected if rr.is_strand_correct())
    target_ok = sum(1 for _r, _c, rr in surjected
                    if rr.surj.path_name and
                    rr.surj.path_name.split("#")[:max(args.target_name_fields, 1)]
                    == _r.target.split("#")[:max(args.target_name_fields, 1)])

    # coverage / capability
    target_samples = {c.sample for _r, c, rr in surjected}
    nonref_samples = {s for s in target_samples if s not in ref_samples}
    distinct_targets = {rr.surj.path_name for _r, _c, rr in surjected if rr.surj.path_name}
    nonref_surj = sum(1 for _r, c, rr in surjected if not c.is_reference)

    bar = "=" * 72
    print(bar)
    print("  PANGENOME SERVER — END-TO-END STRESS TEST")
    print("  giraffe mapping  →  anchor build  →  surjection onto ANY haplotype")
    print(bar)
    print(f"  graph        : {os.path.basename(args.gbz)}")
    print(f"  drawn from   : {n_paths} source haplotype paths")
    lengths = ",".join(str(x) for x in sorted({c.length for _r, c, _rr in rows}))
    print(f"  reads        : {total}   lengths(bp): {lengths}   "
          f"rev-strand: {args.revcomp_fraction:.0%}   "
          f"with-errors: {args.error_fraction:.0%} @ {args.error_rate:.1%}/base")
    print("  method       : each read sampled at a KNOWN base on a KNOWN haplotype,")
    print("                 mapped to the graph, surjected back — must land on that base.")
    print()

    print("  PIPELINE OUTCOMES")
    print(f"    mapped to graph        : {nm} / {total}   ({_pct(nm, total)})")
    print(f"    surjected onto target  : {ns} / {nm}   ({_pct(ns, nm)} of mapped)")
    print()

    print(f"  CORRECTNESS  (over {ns} surjected reads; tolerance = {tol} bp)")
    print(f"    landed on the EXACT base : {exact} / {ns}   ({_pct(exact, ns)})")
    print(f"    strand correct           : {strand_ok} / {ns}   ({_pct(strand_ok, ns)})")
    print(f"    target haplotype correct : {target_ok} / {ns}   ({_pct(target_ok, ns)})")
    print()
    print("    exact-base rate by read length :  "
          + _fmt_breakdown(_breakdown(rows, lambda c: c.length, tol)))
    print("    exact-base rate by strand      :  "
          + _fmt_breakdown(_breakdown(rows, lambda c: ("forward" if c.strand == "+" else "reverse"), tol)))
    print("    exact-base rate by read type   :  "
          + _fmt_breakdown(_breakdown(rows, lambda c: ("error" if c.errored else "clean"), tol)))
    print("    exact-base rate by target      :  "
          + _fmt_breakdown(_breakdown(rows, lambda c: ("reference" if c.is_reference else "non-reference"), tol)))
    print()

    print("  CAPABILITY  (what reference-only surjection cannot do)")
    print(f"    reads surjected onto NON-reference haplotypes : {nonref_surj}")
    print(f"    distinct haplotype samples targeted           : {len(target_samples)} "
          f"({len(nonref_samples)} non-reference)")
    print(f"    distinct target paths (haplotype contigs)     : {len(distinct_targets)}")
    print("    → `vg giraffe --ref-paths` can only surject onto reference paths;")
    print("      this pipeline targets any of them, chosen at query time.")
    print()

    # performance (exclude warmup)
    timed = [rr for _r, _c, rr in rows if rr.mapped and not rr.error]
    warm = timed[args.warmup:] if len(timed) > args.warmup else timed
    if warm:
        print(f"  PERFORMANCE  (serial, single pipeline; {len(warm)} reads, "
              f"{min(args.warmup, len(timed))} warmup excluded)")
        print(f"    map           : {_fmt(timing_stats([r.map_ms for r in warm]))}")
        print(f"    build_anchors : {_fmt(timing_stats([r.build_ms for r in warm if r.build_ms > 0]))}")
        print(f"    surject       : {_fmt(timing_stats([r.surject_ms for r in warm if r.surject_ms > 0]))}")
        print(f"    total / read  : {_fmt(timing_stats([r.total_ms for r in warm]))}")
        tot = sum(r.total_ms for r in warm)
        if tot > 0:
            print(f"    throughput    : {1000.0 * len(warm) / tot:.1f} reads/s")
        print()

    # show any reads that surjected but to the wrong base (should be ~none)
    wrong = [(r, c, rr) for r, c, rr in surjected if not rr.is_pos_correct(tol)]
    if wrong:
        print(f"  POSITION MISMATCHES ({len(wrong)}; showing up to {args.max_mismatches_shown})")
        for r, c, rr in wrong[:args.max_mismatches_shown]:
            d = "na" if rr.surj.position is None else str(rr.surj.position - r.offset)
            print(f"    [{r.name}] {r.source_name}  expect={r.offset} got={rr.surj.position} "
                  f"(Δ{d}) len={c.length} {c.strand} err={c.n_errors} cigar={rr.surj.cigar}")
        print()

    if args.out_tsv:
        _write_tsv(rows, args.out_tsv, tol)
        print(f"  per-read results → {args.out_tsv}")
        print()

    exact_rate = (exact / ns) if ns else 0.0
    strand_rate = (strand_ok / ns) if ns else 0.0
    target_rate = (target_ok / ns) if ns else 0.0
    ok = (ns > 0 and exact_rate >= args.min_exact
          and strand_rate >= 0.999 and target_rate >= 0.999)
    print(bar)
    if ok:
        print(f"  VERDICT: PASS — {_pct(exact, ns)} of surjected reads landed on the exact "
              f"base, across {len(target_samples)} haplotypes.")
    else:
        why = []
        if not ns:
            why.append("no reads surjected")
        if ns and exact_rate < args.min_exact:
            why.append(f"exact-base {exact_rate:.1%} < {args.min_exact:.0%}")
        if ns and strand_rate < 0.999:
            why.append(f"strand {strand_rate:.1%}")
        if ns and target_rate < 0.999:
            why.append(f"target {target_rate:.1%}")
        print(f"  VERDICT: FAIL — {'; '.join(why)}")
    print(bar)
    return 0 if ok else 1


def _write_tsv(rows: List[Tuple[Read, Category, BenchResult]], path: str, tol: int) -> None:
    cols = ["read", "source_name", "target", "expect_offset", "length", "strand",
            "n_errors", "is_reference", "mapped", "surj_status", "surj_path",
            "surj_pos", "pos_diff", "pos_exact", "strand_ok", "surj_cigar",
            "map_ms", "build_ms", "surject_ms"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r, c, rr in rows:
            d = "" if rr.surj.position is None else str(rr.surj.position - r.offset)
            row = [
                r.name, r.source_name, r.target, r.offset, c.length, c.strand,
                c.n_errors, int(c.is_reference), int(rr.mapped), rr.surj.status,
                rr.surj.path_name, "" if rr.surj.position is None else rr.surj.position,
                d, int(rr.is_pos_correct(tol)), int(rr.is_strand_correct()),
                rr.surj.cigar, f"{rr.map_ms:.2f}", f"{rr.build_ms:.2f}", f"{rr.surject_ms:.2f}",
            ]
            fh.write("\t".join(str(x) for x in row) + "\n")


# ──────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────

def main() -> int:
    args = parse_args()
    rng = random.Random(args.seed)
    lengths = [int(x) for x in args.read_lengths.split(",") if x.strip()]
    ref_samples = {s for s in args.reference_samples.split(",") if s.strip()}

    # Generate (or load cached) reads.
    reads: Optional[List[Read]] = None
    if args.reads_cache and os.path.exists(args.reads_cache):
        print(f"Loading cached reads from {args.reads_cache}…", file=sys.stderr)
        reads = load_reads(args.reads_cache) or None

    if reads is None:
        if args.fasta:
            print("Loading sequences from FASTA…", file=sys.stderr)
            named = load_from_fasta(args.fasta, args.name_filter, args.max_haplotypes,
                                    args.max_seq_len, rng)
        elif args.sequences and args.names:
            print("Loading sequences (gbz_extract + names)…", file=sys.stderr)
            with open(args.names) as fh:
                names = [ln.strip() for ln in fh if ln.strip()]
            named = load_from_pair(args.sequences, names, args.both_orientations,
                                   args.name_filter, args.max_haplotypes,
                                   args.max_seq_len, rng)
        else:
            print("ERROR: provide --fasta, or both --sequences and --names.", file=sys.stderr)
            return 2
        if not named:
            print("ERROR: no source sequences (check --name-filter / inputs).", file=sys.stderr)
            return 2
        print(f"  {len(named)} source path(s); {sum(len(s) for _, s in named)} bp", file=sys.stderr)
        reads = generate_reads(named, lengths, args.reads_per_length,
                               args.revcomp_fraction, args.error_fraction,
                               args.error_rate, args.max_n_fraction,
                               args.target_name_fields, rng)
        if not reads:
            print("ERROR: could not generate reads (sequences shorter than read lengths?).",
                  file=sys.stderr)
            return 2
        if args.reads_cache:
            save_reads(reads, args.reads_cache)
            print(f"  cached {len(reads)} reads to {args.reads_cache}", file=sys.stderr)

    n_source_paths = len({r.source_name for r in reads})
    print(f"  {len(reads)} reads from {n_source_paths} source paths", file=sys.stderr)

    cfg = GiraffeServerConfig(
        vg_binary=args.vg, gbz_path=args.gbz, minimizer_path=args.minimizer,
        distance_path=args.dist, zipcode_path=args.zipcodes,
        threads=args.threads, max_multimaps=args.max_multimaps,
        batch_size=args.batch_size, output_timeout_s=args.timeout,
        surject_target_paths=[],  # anchor path needs no pre-indexed targets
    )
    mw = GiraffeServerMiddleware(cfg)
    mw.start()
    rows: List[Tuple[Read, Category, BenchResult]] = []
    try:
        print("Waiting for giraffe-server…", file=sys.stderr)
        mw.wait_until_ready()
        print("giraffe-server ready.", file=sys.stderr)

        import liftover_ext
        coord = liftover_ext.Index()
        print("Loading coordinate index…", file=sys.stderr)
        coord.load(args.gbz, args.ri, args.tags, args.gbwt_ri, args.t1, args.t2)
        print("Coordinate index loaded. Running stress test…", file=sys.stderr)

        for i, read in enumerate(reads):
            rr = run_one(mw, coord, read)
            rows.append((read, Category(read, ref_samples), rr))
            if (i + 1) % args.progress_every == 0:
                ok = sum(1 for _r, _c, x in rows if x.surj.ok and x.is_pos_correct(args.pos_tolerance))
                print(f"  …{i + 1}/{len(reads)}  exact-so-far={ok}", file=sys.stderr)
    finally:
        mw.stop()

    print(file=sys.stderr)
    return report(rows, args, n_source_paths, ref_samples)


if __name__ == "__main__":
    raise SystemExit(main())
