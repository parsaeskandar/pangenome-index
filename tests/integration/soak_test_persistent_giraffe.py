#!/usr/bin/env python3
from __future__ import annotations

# Repo root (holds liftover_ext.so + the middleware package) on sys.path, so
# this script runs from anywhere after being moved under tests/integration/.
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import argparse
import os
import random
import time
from typing import List

from middleware.giraffe_server_middleware import GiraffeServerConfig, GiraffeServerMiddleware


DNA = "ACGT"


def random_seq(n: int, rng: random.Random) -> str:
    return "".join(rng.choice(DNA) for _ in range(n))


def read_rss_kb(pid: int) -> int:
    status_path = f"/proc/{pid}/status"
    try:
        with open(status_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    parts = line.split()
                    return int(parts[1])
    except FileNotFoundError:
        return -1
    return -1


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Soak test: persistent in-memory giraffe middleware")
    p.add_argument("--vg", required=True)
    p.add_argument("--gbz", required=True)
    p.add_argument("--minimizer", required=True)
    p.add_argument("--dist", required=True)
    p.add_argument("--zipcodes", required=True)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--batch-size", type=int, default=64)
    p.add_argument("--max-multimaps", type=int, default=1)
    p.add_argument("--rounds", type=int, default=500, help="Number of query rounds")
    p.add_argument("--reads-per-round", type=int, default=8)
    p.add_argument("--read-len", type=int, default=50)
    p.add_argument("--seed", type=int, default=13)
    p.add_argument("--sleep-ms", type=int, default=0, help="Optional delay between rounds")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    cfg = GiraffeServerConfig(
        vg_binary=args.vg,
        gbz_path=args.gbz,
        minimizer_path=args.minimizer,
        distance_path=args.dist,
        zipcode_path=args.zipcodes,
        threads=args.threads,
        max_multimaps=args.max_multimaps,
        batch_size=args.batch_size,
    )
    rng = random.Random(args.seed)

    mw = GiraffeServerMiddleware(cfg)
    t0 = time.time()
    failures = 0
    total_reads = 0

    try:
        mw.start()
        if not mw.is_running() or mw.pid is None:
            print("FAIL: middleware process did not start")
            return 1
        pid0 = mw.pid
        rss0 = read_rss_kb(pid0)
        print(f"Started vg giraffe-server pid={pid0} rss_kb={rss0}")

        for r in range(1, args.rounds + 1):
            reads: List[tuple[str, str, str]] = []
            for i in range(args.reads_per_round):
                seq = random_seq(args.read_len, rng)
                reads.append((f"round{r}_read{i}", seq, "I" * len(seq)))

            out = mw.map_reads(reads)
            total_reads += len(reads)

            if mw.pid != pid0 or not mw.is_running():
                print(f"FAIL: process changed or exited at round {r} (pid now {mw.pid})")
                return 1
            if len(out) != len(reads):
                print(f"FAIL: output size mismatch at round {r}: got {len(out)} expected {len(reads)}")
                return 1

            round_ok = all(len(x) >= 1 for x in out)
            if not round_ok:
                failures += 1

            if r % 25 == 0 or r == 1 or r == args.rounds:
                rss = read_rss_kb(pid0)
                elapsed = time.time() - t0
                rps = total_reads / max(elapsed, 1e-6)
                print(
                    f"round={r}/{args.rounds} pid={pid0} rss_kb={rss} "
                    f"total_reads={total_reads} reads_per_sec={rps:.2f} failures={failures}"
                )

            if args.sleep_ms > 0:
                time.sleep(args.sleep_ms / 1000.0)

    finally:
        mw.stop()

    elapsed = time.time() - t0
    print(
        f"PASS: completed rounds={args.rounds} total_reads={total_reads} "
        f"elapsed_s={elapsed:.2f} failures={failures}"
    )
    return 0 if failures == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
