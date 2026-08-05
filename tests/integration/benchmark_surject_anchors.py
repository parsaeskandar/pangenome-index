#!/usr/bin/env python3
"""
benchmark_surject_anchors.py — ground-truth correctness + timing benchmark for
the giraffe + surject-with-anchors pipeline.

Idea (step 1): sample reads from known haplotype sequences at known offsets,
then surject each read back onto the SAME haplotype it came from. Because we
know the sample offset, we know the answer: the surjected position must equal
that offset. No second tool needed — the haplotype sequence is its own ground
truth. The anchor pipeline is the right fit here because it needs no
pre-indexing, so every read can target a different haplotype at query time.

For each read it times three stages and checks correctness:
  1. map           — giraffe-server graph mapping (no surjection)   [map_ms]
  2. build_anchors — liftover_ext.build_surject_anchors_full         [build_ms]
  3. surject       — SURJECT_WITH_ANCHORS on giraffe-server          [surject_ms]
correctness = surjected position within --pos-tolerance of the sample offset,
on the expected strand, on the expected path.

Inputs (two ways to provide name↔sequence pairs):
  A) --fasta FILE        : `vg paths --extract-fasta` output (>NAME\\nSEQ).
                           Intrinsically paired; recommended (no ordering risk).
  B) --sequences FILE --names FILE [--both-orientations]
                         : gbz_extract output + `vg paths -L` output. Pairs by
                           order: names line i ↔ sequences line (2i if -b else i).

The anchor target for a path named "A#B#C#D" is its T1 base name = the first
--target-name-fields (#-separated) fields, default 3 → "A#B#C" (what
Translation Table 1 stores). Override with --target-name-fields.

Example:
  $VG paths --extract-fasta -x $W/graph.gbz -S HG002 > hg002.fa     # subset
  python3 benchmark_surject_anchors.py \\
    --vg ../vg --gbz $W/graph.gbz \\
    --minimizer $W/...min --dist $W/...dist --zipcodes $W/...zipcodes \\
    --ri $W/rlbwt_rindex.ri --tags $W/sampled.tags --gbwt-ri $W/fastlocate.ri \\
    --t1 $W/table1.t1 --t2 $W/table2.t2 \\
    --fasta hg002.fa --num-reads 500 --read-length 150 \\
    --max-haplotypes 20 --warmup 10 --out-tsv bench.tsv
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
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# Reuse the middleware loader + GAF helpers from the comparison test (no native
# dependency at import time; liftover_ext is imported lazily in main()).
from test_anchor_surject import (
    GiraffeServerConfig,
    GiraffeServerMiddleware,
    Surjection,
    revcomp,
    strip_surj_tags,
)


# ──────────────────────────────────────────────────────────────────────────
# Arguments
# ──────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Ground-truth correctness + timing benchmark for giraffe + surject-with-anchors.",
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
    # name↔sequence input (choose A or B)
    p.add_argument("--fasta", help="A) vg paths --extract-fasta output (recommended).")
    p.add_argument("--sequences", help="B) gbz_extract sequences file (one seq per line).")
    p.add_argument("--names", help="B) vg paths -L output (one path name per line).")
    p.add_argument("--both-orientations", action="store_true",
                   help="B) sequences file was made with gbz_extract -b (2 lines/path).")
    # which haplotypes/paths to draw reads from
    p.add_argument("--max-haplotypes", type=int, default=20,
                   help="Sample reads from at most this many paths (randomly chosen).")
    p.add_argument("--name-filter", default="",
                   help="Only consider paths whose name matches this regex (e.g. '^HG002' or 'chr20').")
    p.add_argument("--max-seq-len", type=int, default=10_000_000,
                   help="Truncate each path sequence to this length (memory guard).")
    p.add_argument("--target-name-fields", type=int, default=3,
                   help="Use the first N '#'-fields of the path name as the anchor target "
                        "(T1 base name = sample#phase#contig = 3). 0 = use the full name.")
    # read sampling
    p.add_argument("--num-reads", type=int, default=500,
                   help="Reads to sample. With --read-lengths this is PER length.")
    p.add_argument("--read-length", type=int, default=150,
                   help="Single read length (used only when --read-lengths is not given).")
    p.add_argument("--read-lengths", default="",
                   help="Comma-separated read lengths to benchmark together in ONE sampling "
                        "pass, e.g. 150,500,1000,5000,10000,20000,50000,100000. Overrides "
                        "--read-length; --num-reads then applies per length. Results are "
                        "broken down by length. Only sequences at least as long as a given "
                        "length contribute reads of that length.")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--revcomp-fraction", type=float, default=0.0,
                   help="Fraction of reads sampled as reverse complement (reverse-strand "
                        "surjection is supported).")
    p.add_argument("--max-n-fraction", type=float, default=0.1)
    # server config
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--batch-size", type=int, default=64)
    p.add_argument("--max-multimaps", type=int, default=1,
                   help="Max alignments giraffe reports per read. Default 1: we only "
                        "surject the primary, so reporting up to 100 (the BLAT default) "
                        "just ships ~100× more huge GAF lines per frame — slower and a "
                        "framing-desync risk. Set 0 for the engine default, or higher to "
                        "inspect multimaps.")
    p.add_argument("--timeout", type=float, default=60.0,
                   help="Per-call giraffe-server response timeout (s). If "
                        "surject_with_anchors hangs, a low value here surfaces "
                        "a vg binary that lacks the SURJECT_WITH_ANCHORS handler.")
    # correctness + timing + reporting
    p.add_argument("--pos-tolerance", type=int, default=0,
                   help="Max |surjected pos - sample offset| still counted correct.")
    p.add_argument("--warmup", type=int, default=10,
                   help="Exclude the first N reads from timing aggregates (server warmup).")
    p.add_argument("--min-correct", type=float, default=0.90,
                   help="Minimum position-correct rate (over surjected reads) to PASS.")
    p.add_argument("--out-tsv", default="")
    p.add_argument("--reads-cache", default="",
                   help="Cache the sampled reads here. If the file exists, reads are "
                        "loaded from it and the (slow) --sequences/--fasta streaming is "
                        "skipped entirely — ideal when rerunning to measure surjection "
                        "while keeping the read set fixed. Delete the file to resample.")
    p.add_argument("--max-mismatches-shown", type=int, default=20)
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


# ──────────────────────────────────────────────────────────────────────────
# Name ↔ sequence loading (subset, memory-bounded)
# ──────────────────────────────────────────────────────────────────────────

def _select_indices(names: List[str], name_filter: str, max_hap: int,
                    rng: random.Random) -> List[int]:
    if name_filter:
        pat = re.compile(name_filter)
        idx = [i for i, n in enumerate(names) if pat.search(n)]
    else:
        idx = list(range(len(names)))
    if not idx:
        return []
    if len(idx) > max_hap:
        idx = rng.sample(idx, max_hap)
    return sorted(idx)


def load_from_fasta(path: str, name_filter: str, max_hap: int,
                    max_seq_len: int, rng: random.Random) -> List[Tuple[str, str]]:
    """Stream a FASTA; keep a random subset of up-to-max_hap records (reservoir
    sample over headers so we don't hold the whole file). Returns [(name, seq)]."""
    pat = re.compile(name_filter) if name_filter else None
    # First pass: collect matching headers + their record index (cheap).
    headers: List[Tuple[int, str]] = []
    with open(path) as fh:
        rec = -1
        for line in fh:
            if line.startswith(">"):
                rec += 1
                name = line[1:].strip().split()[0] if line[1:].strip() else ""
                if name and (pat is None or pat.search(name)):
                    headers.append((rec, name))
    if not headers:
        return []
    chosen = headers if len(headers) <= max_hap else rng.sample(headers, max_hap)
    want = {rec: name for rec, name in chosen}
    # Second pass: capture the sequences for chosen records only.
    out: List[Tuple[str, str]] = []
    with open(path) as fh:
        rec = -1
        cur_name: Optional[str] = None
        cur: List[str] = []
        def flush():
            if cur_name is not None and rec in want:
                seq = "".join(cur).upper()[:max_seq_len]
                out.append((want[rec], seq))
        for line in fh:
            if line.startswith(">"):
                flush()
                rec += 1
                cur_name = line[1:].strip().split()[0] if line[1:].strip() else None
                cur = []
            else:
                if rec in want:
                    cur.append(line.strip())
        flush()
    return out


def load_from_pair(seq_path: str, names: List[str], both_orient: bool,
                   name_filter: str, max_hap: int, max_seq_len: int,
                   rng: random.Random) -> List[Tuple[str, str]]:
    """gbz_extract sequences + vg-paths names, paired by order. names line i ↔
    sequence line (2i if -b else i)."""
    sel = _select_indices(names, name_filter, max_hap, rng)
    if not sel:
        return []
    want_line = {(2 * i if both_orient else i): i for i in sel}
    captured: Dict[int, str] = {}
    with open(seq_path) as fh:
        for ln, line in enumerate(fh):
            if ln in want_line:
                captured[want_line[ln]] = line.strip().upper()[:max_seq_len]
                if len(captured) == len(want_line):
                    break
    return [(names[i], captured[i]) for i in sel if i in captured and captured[i]]


# ──────────────────────────────────────────────────────────────────────────
# Read sampling
# ──────────────────────────────────────────────────────────────────────────

@dataclass
class Read:
    name: str
    seq: str
    source_name: str   # full path name the read was drawn from
    target: str        # anchor target (T1 base name)
    offset: int        # 0-based start within the source path sequence
    is_revcomp: bool


def base_target(path_name: str, fields: int) -> str:
    if fields <= 0:
        return path_name
    parts = path_name.split("#")
    return "#".join(parts[:fields]) if len(parts) > fields else path_name


def sample_reads(named_seqs: List[Tuple[str, str]], n: int, read_len: int,
                 target_fields: int, revcomp_fraction: float,
                 max_n_fraction: float, rng: random.Random) -> List[Read]:
    eligible = [(nm, s) for nm, s in named_seqs if len(s) > 0]
    if not eligible:
        return []
    weights = [len(s) for _, s in eligible]
    reads: List[Read] = []
    attempts = 0
    max_attempts = n * 20
    while len(reads) < n and attempts < max_attempts:
        attempts += 1
        nm, s = rng.choices(eligible, weights=weights, k=1)[0]
        if len(s) <= read_len:
            start, window = 0, s
        else:
            start = rng.randint(0, len(s) - read_len)
            window = s[start:start + read_len]
        if window.count("N") > max_n_fraction * len(window):
            continue
        is_rc = rng.random() < revcomp_fraction
        seq = revcomp(window) if is_rc else window
        reads.append(Read(
            name=f"r{len(reads)}_{start}{'_rc' if is_rc else ''}",
            seq=seq, source_name=nm,
            target=base_target(nm, target_fields),
            offset=start, is_revcomp=is_rc,
        ))
    return reads


def sample_reads_multi(named_seqs: List[Tuple[str, str]], lengths: List[int],
                       n_per_length: int, target_fields: int,
                       revcomp_fraction: float, max_n_fraction: float,
                       rng: random.Random) -> List[Read]:
    """Sample n_per_length reads for EACH length, from the already-loaded
    sequences (one file pass already happened in load_*). Only sequences at
    least as long as a given length contribute reads of that length, so the
    read length equals len(seq) exactly — the report buckets cleanly by it.
    Read names are globally unique and encode the length (`_L<len>`)."""
    eligible = [(nm, s) for nm, s in named_seqs if len(s) > 0]
    if not eligible:
        return []
    reads: List[Read] = []
    idx = 0
    for L in lengths:
        pool = [(nm, s) for nm, s in eligible if len(s) >= L]
        if not pool:
            print(f"  WARNING: no source sequence >= {L} bp; skipping length {L}.",
                  file=sys.stderr)
            continue
        weights = [len(s) for _, s in pool]
        made = 0
        attempts = 0
        cap = n_per_length * 40
        while made < n_per_length and attempts < cap:
            attempts += 1
            nm, s = rng.choices(pool, weights=weights, k=1)[0]
            start = 0 if len(s) == L else rng.randint(0, len(s) - L)
            window = s[start:start + L]
            if window.count("N") > max_n_fraction * L:
                continue
            is_rc = rng.random() < revcomp_fraction
            seq = revcomp(window) if is_rc else window
            reads.append(Read(
                name=f"r{idx}_{start}_L{L}{'_rc' if is_rc else ''}",
                seq=seq, source_name=nm,
                target=base_target(nm, target_fields),
                offset=start, is_revcomp=is_rc,
            ))
            idx += 1
            made += 1
        if made < n_per_length:
            print(f"  WARNING: only sampled {made}/{n_per_length} reads of length {L}.",
                  file=sys.stderr)
    rng.shuffle(reads)  # interleave lengths so warmup isn't all one size
    return reads


def save_reads(reads: List[Read], path: str) -> None:
    """Persist sampled reads (with their ground-truth metadata) as TSV so a
    later run can skip the expensive sequence-file streaming."""
    with open(path, "w") as fh:
        fh.write("# name\tseq\tsource_name\ttarget\toffset\tis_revcomp\n")
        for r in reads:
            fh.write(f"{r.name}\t{r.seq}\t{r.source_name}\t{r.target}\t"
                     f"{r.offset}\t{int(r.is_revcomp)}\n")


def load_reads(path: str) -> List[Read]:
    """Load reads previously written by save_reads()."""
    out: List[Read] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 6:
                continue
            name, seq, source_name, target, offset, rc = parts
            out.append(Read(name=name, seq=seq, source_name=source_name,
                            target=target, offset=int(offset),
                            is_revcomp=bool(int(rc))))
    return out


# ──────────────────────────────────────────────────────────────────────────
# Per-read pipeline + timing
# ──────────────────────────────────────────────────────────────────────────

@dataclass
class BenchResult:
    read: Read
    mapped: bool = False
    map_ms: float = 0.0
    build_ms: float = 0.0
    surject_ms: float = 0.0
    n_anchors: int = 0
    build_status: str = ""
    surj: Surjection = field(default_factory=Surjection)
    error: str = ""
    # find_sequences_for_tag LF-cost diagnostics (from build_surject_anchors_full).
    find_seq_calls: int = 0
    find_seq_runs: int = 0
    find_seq_lf_steps: int = 0
    find_seq_visits: int = 0
    last_run_nav_steps: int = 0
    last_run_length: int = 0
    # target-path GBWT LF walk diagnostics (the other candidate bottleneck).
    walk_lf_steps: int = 0
    walk_span: int = 0
    first_anchor_base: int = 0
    last_anchor_base: int = 0
    # wall-clock attribution of build_ms (the definitive breakdown).
    find_seq_ms: float = 0.0
    decompress_sa_ms: float = 0.0
    walk_ms: float = 0.0
    decompress_sa_calls: int = 0
    decompress_sa_entries: int = 0
    n_target_subpaths: int = 0
    n_source_mappings: int = 0

    @property
    def total_ms(self) -> float:
        return self.map_ms + self.build_ms + self.surject_ms

    def is_pos_correct(self, tol: int) -> bool:
        return (self.surj.ok and self.surj.position is not None
                and abs(self.surj.position - self.read.offset) <= tol)

    def is_strand_correct(self) -> bool:
        if not self.surj.ok or self.surj.strand is None:
            return False
        return self.surj.strand == (1 if self.read.is_revcomp else 0)


def run_one(mw: GiraffeServerMiddleware, coord, read: Read) -> BenchResult:
    rr = BenchResult(read=read)
    qual = "I" * len(read.seq)

    t0 = time.perf_counter()
    try:
        out = mw.map_reads([(read.name, read.seq, qual)])
    except Exception as exc:  # noqa: BLE001
        rr.error = f"map: {exc}"
        return rr
    rr.map_ms = (time.perf_counter() - t0) * 1e3
    alns = out[0] if out else []
    if not alns:
        return rr
    rr.mapped = True
    graph_gaf = strip_surj_tags(alns[0])  # clean graph alignment (primary)

    t0 = time.perf_counter()
    try:
        build = coord.build_surject_anchors_full(graph_gaf, read.target)
    except Exception as exc:  # noqa: BLE001
        rr.error = f"build: {exc}"
        return rr
    rr.build_ms = (time.perf_counter() - t0) * 1e3
    rr.build_status = build.status
    anchors = list(build.anchors)
    rr.n_anchors = len(anchors)
    path_len = build.target_path_length
    # find_sequences_for_tag LF-cost diagnostics.
    rr.find_seq_calls     = build.find_seq_calls
    rr.find_seq_runs      = build.find_seq_runs
    rr.find_seq_lf_steps  = build.find_seq_lf_steps
    rr.find_seq_visits    = build.find_seq_visits
    rr.last_run_nav_steps = build.last_run_nav_steps
    rr.last_run_length    = build.last_run_length
    rr.walk_lf_steps      = build.walk_lf_steps
    rr.walk_span          = build.walk_span
    rr.first_anchor_base  = build.first_anchor_base
    rr.last_anchor_base   = build.last_anchor_base
    rr.find_seq_ms          = build.find_seq_ms
    rr.decompress_sa_ms     = build.decompress_sa_ms
    rr.walk_ms              = build.walk_ms
    rr.decompress_sa_calls  = build.decompress_sa_calls
    rr.decompress_sa_entries = build.decompress_sa_entries
    rr.n_target_subpaths    = build.n_target_subpaths
    rr.n_source_mappings    = build.n_source_mappings
    if not anchors:
        rr.surj = Surjection(status=f"no_anchors({build.status})")
        return rr

    t0 = time.perf_counter()
    try:
        lines = mw.surject_with_anchors(graph_gaf, anchors, read.target,
                                        target_path_length=path_len,
                                        read_name=read.name + "_anc")
    except Exception as exc:  # noqa: BLE001
        rr.error = f"surject: {exc}"
        return rr
    rr.surject_ms = (time.perf_counter() - t0) * 1e3
    if lines:
        rr.surj = Surjection.from_gaf(lines[0])
    return rr


# ──────────────────────────────────────────────────────────────────────────
# Stats + reporting
# ──────────────────────────────────────────────────────────────────────────

def _pct(sorted_vals: List[float], q: float) -> float:
    if not sorted_vals:
        return 0.0
    k = max(0, min(len(sorted_vals) - 1, int(round(q * (len(sorted_vals) - 1)))))
    return sorted_vals[k]


def timing_stats(vals: List[float]) -> Dict[str, float]:
    if not vals:
        return {"n": 0, "mean": 0, "p50": 0, "p90": 0, "p99": 0, "min": 0, "max": 0}
    sv = sorted(vals)
    return {
        "n": len(sv),
        "mean": sum(sv) / len(sv),
        "p50": _pct(sv, 0.50),
        "p90": _pct(sv, 0.90),
        "p99": _pct(sv, 0.99),
        "min": sv[0],
        "max": sv[-1],
    }


def _fmt(s: Dict[str, float]) -> str:
    return (f"mean={s['mean']:.2f}  p50={s['p50']:.2f}  p90={s['p90']:.2f}  "
            f"p99={s['p99']:.2f}  min={s['min']:.2f}  max={s['max']:.2f}  (ms)")


def report(results: List[BenchResult], args: argparse.Namespace,
           n_paths: int) -> int:
    n_total = len(results)
    mapped = [r for r in results if r.mapped]
    surjected = [r for r in results if r.surj.ok]
    n_errors = sum(1 for r in results if r.error)

    print("=" * 70)
    print("giraffe + surject-with-anchors — ground-truth benchmark")
    print("=" * 70)
    print(f"paths sampled from : {n_paths}")
    print(f"reads tested       : {n_total}  "
          f"(read_length={args.read_length}, revcomp={args.revcomp_fraction:.0%}, "
          f"target=first {args.target_name_fields or 'all'} name fields)")
    print(f"mapped             : {len(mapped)} / {n_total}")
    print(f"surjected (sj:ok)  : {len(surjected)} / {n_total}")
    if n_errors:
        print(f"per-read errors    : {n_errors}")
    print()

    # Correctness (over surjected reads — that's where a position exists).
    pos_ok = sum(1 for r in surjected if r.is_pos_correct(args.pos_tolerance))
    strand_ok = sum(1 for r in surjected if r.is_strand_correct())
    path_ok = sum(1 for r in surjected if r.surj.path_name and
                  r.surj.path_name.split("#")[:max(args.target_name_fields, 1)]
                  == r.read.target.split("#")[:max(args.target_name_fields, 1)])
    n_s = len(surjected)

    def pc(x):
        return f"{(100.0 * x / n_s):.1f}%" if n_s else "n/a"

    print(f"Correctness (over {n_s} surjected reads, pos tol={args.pos_tolerance}):")
    print(f"  position correct : {pos_ok} / {n_s}  ({pc(pos_ok)})")
    print(f"  strand correct   : {strand_ok} / {n_s}  ({pc(strand_ok)})")
    print(f"  target path match: {path_ok} / {n_s}  ({pc(path_ok)})")

    # |pos - offset| histogram.
    buckets = {"0": 0, "1-5": 0, "6-50": 0, "51+": 0}
    for r in surjected:
        if r.surj.position is None:
            continue
        d = abs(r.surj.position - r.read.offset)
        k = "0" if d == 0 else "1-5" if d <= 5 else "6-50" if d <= 50 else "51+"
        buckets[k] += 1
    print("  |pos diff| hist  : " + "  ".join(f"{k}:{v}" for k, v in buckets.items()))

    # build-status distribution (why some reads produced no anchors).
    bs: Dict[str, int] = {}
    for r in mapped:
        bs[r.build_status or "(none)"] = bs.get(r.build_status or "(none)", 0) + 1
    print("  build status     : " + "  ".join(f"{k}:{v}" for k, v in sorted(bs.items())))
    na = [r.n_anchors for r in surjected if r.n_anchors]
    if na:
        print(f"  anchors/read     : mean={sum(na)/len(na):.1f}  min={min(na)}  max={max(na)}")
    print()

    # Timing (exclude warmup reads from aggregates).
    timed = [r for r in mapped if not r.error]
    warm = timed[args.warmup:] if len(timed) > args.warmup else timed
    print(f"Timing over {len(warm)} reads (excluding {min(args.warmup, len(timed))} warmup):")
    print(f"  map           : {_fmt(timing_stats([r.map_ms for r in warm]))}")
    print(f"  build_anchors : {_fmt(timing_stats([r.build_ms for r in warm if r.build_ms > 0]))}")
    print(f"  surject       : {_fmt(timing_stats([r.surject_ms for r in warm if r.surject_ms > 0]))}")
    print(f"  total/read    : {_fmt(timing_stats([r.total_ms for r in warm]))}")
    tot = sum(r.total_ms for r in warm)
    if tot > 0:
        print(f"  throughput    : {1000.0 * len(warm) / tot:.1f} reads/s (serial, single pipeline)")
    print()

    # Per read-length breakdown — how the pipeline scales with read length.
    # (Bucketed over the same post-warmup reads; p50 ms per stage.)
    lengths_present = sorted({len(r.read.seq) for r in warm})
    if len(lengths_present) > 1:
        print("Per read-length breakdown (p50 ms; correctness over surjected):")
        print(f"  {'length':>8} {'reads':>6} {'surj':>5} {'exact':>7} "
              f"{'map':>7} {'build':>7} {'surj':>8} {'total':>9} {'reads/s':>8}")
        for L in lengths_present:
            g = [r for r in warm if len(r.read.seq) == L]
            gs = [r for r in g if r.surj.ok]
            gx = sum(1 for r in gs if r.is_pos_correct(args.pos_tolerance))
            mp = timing_stats([r.map_ms for r in g])["p50"]
            bd = timing_stats([r.build_ms for r in g if r.build_ms > 0])["p50"]
            sj = timing_stats([r.surject_ms for r in gs if r.surject_ms > 0])["p50"]
            tt = timing_stats([r.total_ms for r in g])["p50"]
            gtot = sum(r.total_ms for r in g)
            rps = (1000.0 * len(g) / gtot) if gtot > 0 else 0.0
            exact_str = f"{(100.0*gx/len(gs)):.0f}%" if gs else "n/a"
            print(f"  {L:>8} {len(g):>6} {len(gs):>5} {exact_str:>7} "
                  f"{mp:>7.1f} {bd:>7.1f} {sj:>8.1f} {tt:>9.1f} {rps:>8.1f}")
        print()

    # Where the build time goes: find_sequences_for_tag LF vs. the target-path
    # walk LF. The one whose ns-per-step is ~constant across reads is the cost.
    measured = [r for r in warm if r.build_ms > 0]
    if measured:
        fs = [r.find_seq_lf_steps for r in measured]
        wk = [r.walk_lf_steps for r in measured]
        vis = [r.find_seq_visits for r in measured]
        spans = [r.walk_span for r in measured]
        print(f"LF cost breakdown over {len(measured)} reads:")
        print(f"  find_seq lf   : {_fmt(timing_stats([float(x) for x in fs]))}")
        print(f"  walk lf       : {_fmt(timing_stats([float(x) for x in wk]))}")
        print(f"  walk span(bp) : {_fmt(timing_stats([float(x) for x in spans]))}")
        print(f"  find_seq vis  : {_fmt(timing_stats([float(x) for x in vis]))}")
        calls = [r.find_seq_calls for r in measured]
        print(f"  find_seq calls/read : mean={sum(calls)/len(calls):.1f}  "
              f"min={min(calls)}  max={max(calls)}")

        # Wall-clock attribution — the definitive answer to "where does build go".
        print(f"  time(ms) find_seq    : {_fmt(timing_stats([r.find_seq_ms for r in measured]))}")
        print(f"  time(ms) decompressSA: {_fmt(timing_stats([r.decompress_sa_ms for r in measured]))}")
        print(f"  time(ms) walk        : {_fmt(timing_stats([r.walk_ms for r in measured]))}")
        ent = [r.decompress_sa_entries for r in measured]
        print(f"  decompressSA entries : {_fmt(timing_stats([float(x) for x in ent]))}")
        sum_fs = sum(r.find_seq_ms for r in measured)
        sum_ds = sum(r.decompress_sa_ms for r in measured)
        sum_wk = sum(r.walk_ms for r in measured)
        sum_build = sum(r.build_ms for r in measured)
        if sum_build > 0:
            print(f"  share of total build : find_seq={100*sum_fs/sum_build:.0f}%  "
                  f"decompressSA={100*sum_ds/sum_build:.0f}%  walk={100*sum_wk/sum_build:.0f}%  "
                  f"(accounted={100*(sum_fs+sum_ds+sum_wk)/sum_build:.0f}%)")
        slow = max(measured, key=lambda r: r.build_ms)
        fast = min(measured, key=lambda r: r.build_ms)
        print(f"  slowest: build={slow.build_ms:.0f}ms  find_seq={slow.find_seq_ms:.0f}  "
              f"decompressSA={slow.decompress_sa_ms:.0f}  walk={slow.walk_ms:.0f}  "
              f"ds_entries={slow.decompress_sa_entries}")
        print(f"  fastest: build={fast.build_ms:.0f}ms  find_seq={fast.find_seq_ms:.0f}  "
              f"decompressSA={fast.decompress_sa_ms:.0f}  walk={fast.walk_ms:.0f}  "
              f"ds_entries={fast.decompress_sa_entries}")
    print()

    # Mismatches: surjected but wrong position.
    wrong = [r for r in surjected if not r.is_pos_correct(args.pos_tolerance)]
    if wrong:
        print(f"Position-wrong reads ({len(wrong)}; showing up to {args.max_mismatches_shown}):")
        for r in wrong[:args.max_mismatches_shown]:
            d = ("na" if r.surj.position is None
                 else str(r.surj.position - r.read.offset))
            print(f"  [{r.read.name}] src={r.read.source_name} tgt={r.read.target} "
                  f"expect_pos={r.read.offset} got_pos={r.surj.position} (Δ{d}) "
                  f"sn={r.surj.path_name} sr={r.surj.strand}/{1 if r.read.is_revcomp else 0} "
                  f"sc={r.surj.cigar} n_anchors={r.n_anchors}")
        print()

    if n_errors and args.verbose:
        print("Errors:")
        for r in results:
            if r.error:
                print(f"  [{r.read.name}] {r.error}")
        print()

    if args.out_tsv:
        _write_tsv(results, args.out_tsv, args.pos_tolerance)
        print(f"Per-read results → {args.out_tsv}")
        print()

    correct_rate = (pos_ok / n_s) if n_s else 0.0
    ok = bool(n_s) and correct_rate >= args.min_correct
    verdict = "PASS" if ok else "FAIL"
    reason = ""
    if not n_s:
        reason = " (no reads surjected — check target-name-fields / index pairing)"
    elif not ok:
        reason = f" (position-correct {pc(pos_ok)} < required {args.min_correct:.0%})"
    print(f"{verdict}{reason}")
    return 0 if ok else 1


def _write_tsv(results: List[BenchResult], path: str, tol: int) -> None:
    cols = ["read", "source_name", "target", "offset", "revcomp", "mapped",
            "build_status", "n_anchors", "surj_status", "surj_sn", "surj_sp",
            "pos_diff", "surj_sr", "expect_sr", "surj_sc",
            "pos_correct", "map_ms", "build_ms", "surject_ms", "total_ms", "error"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in results:
            pos_diff = ("" if r.surj.position is None
                        else str(r.surj.position - r.read.offset))
            row = [
                r.read.name, r.read.source_name, r.read.target, r.read.offset,
                int(r.read.is_revcomp), int(r.mapped),
                r.build_status, r.n_anchors, r.surj.status, r.surj.path_name,
                "" if r.surj.position is None else r.surj.position, pos_diff,
                "" if r.surj.strand is None else r.surj.strand,
                1 if r.read.is_revcomp else 0, r.surj.cigar,
                int(r.is_pos_correct(tol)),
                f"{r.map_ms:.3f}", f"{r.build_ms:.3f}", f"{r.surject_ms:.3f}",
                f"{r.total_ms:.3f}", r.error,
            ]
            fh.write("\t".join(str(c) for c in row) + "\n")


# ──────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────

def main() -> int:
    args = parse_args()
    rng = random.Random(args.seed)

    # Fast path: reuse a cached read set and skip the (slow) sequence-file
    # streaming entirely. Ideal when rerunning to measure surjection.
    reads: Optional[List[Read]] = None
    if args.reads_cache and os.path.exists(args.reads_cache):
        print(f"Loading cached reads from {args.reads_cache} "
              f"(skipping sequence streaming)…", file=sys.stderr)
        reads = load_reads(args.reads_cache)
        if not reads:
            print(f"WARNING: {args.reads_cache} held no reads; re-sampling.",
                  file=sys.stderr)
            reads = None

    if reads is None:
        # Load name↔sequence pairs (subset) — this is the expensive part for a
        # whole-genome sequence file, so we only do it on a cache miss.
        if args.fasta:
            print("Loading sequences from FASTA…", file=sys.stderr)
            named = load_from_fasta(args.fasta, args.name_filter, args.max_haplotypes,
                                    args.max_seq_len, rng)
        elif args.sequences and args.names:
            print("Loading sequences from gbz_extract + names…", file=sys.stderr)
            with open(args.names) as fh:
                names = [ln.strip() for ln in fh if ln.strip()]
            named = load_from_pair(args.sequences, names, args.both_orientations,
                                   args.name_filter, args.max_haplotypes,
                                   args.max_seq_len, rng)
        else:
            print("ERROR: provide --fasta, or both --sequences and --names.", file=sys.stderr)
            return 2

        if not named:
            print("ERROR: no sequences selected (check --name-filter / inputs).", file=sys.stderr)
            return 2
        print(f"  {len(named)} path(s) selected; "
              f"{sum(len(s) for _, s in named)} bp", file=sys.stderr)

        if args.read_lengths.strip():
            lengths = [int(x) for x in args.read_lengths.split(",") if x.strip()]
            print(f"  sampling {args.num_reads} reads each for lengths {lengths}…",
                  file=sys.stderr)
            reads = sample_reads_multi(named, lengths, args.num_reads,
                                       args.target_name_fields, args.revcomp_fraction,
                                       args.max_n_fraction, rng)
        else:
            reads = sample_reads(named, args.num_reads, args.read_length,
                                 args.target_name_fields, args.revcomp_fraction,
                                 args.max_n_fraction, rng)
        if not reads:
            print("ERROR: could not sample reads.", file=sys.stderr)
            return 2
        if args.reads_cache:
            save_reads(reads, args.reads_cache)
            print(f"  cached {len(reads)} reads to {args.reads_cache} "
                  f"(reuse with the same --reads-cache to skip streaming)", file=sys.stderr)

    print(f"  {len(reads)} reads from {len(set(r.source_name for r in reads))} paths",
          file=sys.stderr)

    cfg = GiraffeServerConfig(
        vg_binary=args.vg, gbz_path=args.gbz, minimizer_path=args.minimizer,
        distance_path=args.dist, zipcode_path=args.zipcodes,
        threads=args.threads, max_multimaps=args.max_multimaps,
        batch_size=args.batch_size, output_timeout_s=args.timeout,
        surject_target_paths=[],  # anchor path needs NO pre-indexed targets
    )
    mw = GiraffeServerMiddleware(cfg)
    mw.start()
    try:
        print("Waiting for giraffe-server…", file=sys.stderr)
        mw.wait_until_ready()
        print("giraffe-server ready.", file=sys.stderr)

        import liftover_ext
        coord = liftover_ext.Index()
        print("Loading coordinate index…", file=sys.stderr)
        coord.load(args.gbz, args.ri, args.tags, args.gbwt_ri, args.t1, args.t2)
        print("Coordinate index loaded.", file=sys.stderr)

        results: List[BenchResult] = []
        for i, read in enumerate(reads):
            rr = run_one(mw, coord, read)
            results.append(rr)
            # Always show the first few reads (so a stall/slow stage is visible
            # immediately, not after 50 reads), then every 25, plus --verbose.
            if args.verbose or i < 3 or (i + 1) % 25 == 0:
                print(f"[{i}] {read.name} src={read.source_name} tgt={read.target} "
                      f"mapped={rr.mapped} build={rr.build_status} "
                      f"sj={rr.surj.status} sp={rr.surj.position} exp={read.offset} "
                      f"map={rr.map_ms:.0f} build={rr.build_ms:.0f} surj={rr.surject_ms:.0f} ms"
                      + (f"  ERR={rr.error}" if rr.error else ""),
                      file=sys.stderr)
                # find_sequences_for_tag LF-cost breakdown: confirms whether
                # build time tracks the node-enumeration LF steps (the suspected
                # bottleneck). last_run = the final tag run ("vector") iterated:
                # nav LF steps to reach it + how many positions it spans.
                print(f"      find_seq: calls={rr.find_seq_calls} runs={rr.find_seq_runs} "
                      f"lf_steps={rr.find_seq_lf_steps} visits={rr.find_seq_visits} | "
                      f"last_run: nav_lf={rr.last_run_nav_steps} len={rr.last_run_length}",
                      file=sys.stderr)
                # The target-path walk: walk_lf_steps is the suspected real cost;
                # walk_span shows how far apart the chosen anchors landed (a huge
                # span ⇒ boundary nodes recurred and the walk crossed the chromosome).
                print(f"      walk:     lf_steps={rr.walk_lf_steps} span={rr.walk_span} "
                      f"anchors=[{rr.first_anchor_base},{rr.last_anchor_base}]",
                      file=sys.stderr)
                # Wall-clock attribution: which phase actually eats build_ms.
                # decompressSA entries shows how much the GBWT call enumerates.
                accounted = rr.find_seq_ms + rr.decompress_sa_ms + rr.walk_ms
                print(f"      time(ms): find_seq={rr.find_seq_ms:.0f} "
                      f"decompressSA={rr.decompress_sa_ms:.0f} "
                      f"(calls={rr.decompress_sa_calls} entries={rr.decompress_sa_entries}) "
                      f"walk={rr.walk_ms:.0f} | accounted={accounted:.0f}/{rr.build_ms:.0f}",
                      file=sys.stderr)
                # The call explosion: decompressSA calls ≈ subpaths × 2 × mappings.
                exp = rr.n_target_subpaths * 2 * rr.n_source_mappings
                print(f"      cause:    subpaths={rr.n_target_subpaths} "
                      f"source_mappings={rr.n_source_mappings} "
                      f"| subpaths×2×mappings={exp} vs decompressSA_calls={rr.decompress_sa_calls}",
                      file=sys.stderr)
            # Fast-fail: if the very first read can't surject because the server
            # never answered SURJECT_WITH_ANCHORS, don't grind through 500 reads.
            if i == 0 and rr.error.startswith("surject:") and "Timed out" in rr.error:
                print("\nERROR: surject_with_anchors timed out on the first read. "
                      "Your vg binary likely lacks the SURJECT_WITH_ANCHORS handler "
                      "(rebuild vg from Giraffe_server), or raise --timeout if the "
                      "server is just slow.\n", file=sys.stderr)
                break
    finally:
        mw.stop()

    return report(results, args, n_paths=len(named))


if __name__ == "__main__":
    raise SystemExit(main())
