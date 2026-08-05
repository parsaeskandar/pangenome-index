#!/usr/bin/env python3
"""
Integration-style test for the Python coordinate-translation extension.

This script intentionally uses plain assertions so it can run without pytest:
    ./.venv312/bin/python tests/test_python_coordinate_translation.py
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

# Ensure the repo root (where liftover_ext*.so is produced) is importable.
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import liftover_ext


def expect_exception(fn, exc_type, label: str) -> None:
    try:
        fn()
    except exc_type:
        return
    except Exception as err:  # pragma: no cover - defensive path for debugging
        raise AssertionError(
            f"{label}: expected {exc_type.__name__}, got {type(err).__name__}: {err}"
        ) from err
    raise AssertionError(f"{label}: expected {exc_type.__name__}, but no exception was raised")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Test Python coordinate translation extension")
    parser.add_argument(
        "--dataset",
        default="coord_translation_tests/test0_target_loop",
        help="Dataset directory containing graph.gbz, rlbwt_rindex.ri, sampled.tags, "
        "fastlocate.ri, output.t1, output.t2",
    )
    return parser.parse_args()


def build_paths(dataset_dir: Path) -> dict[str, Path]:
    paths = {
        "gbz": dataset_dir / "graph.gbz",
        "ri": dataset_dir / "rlbwt_rindex.ri",
        "tags": dataset_dir / "sampled.tags",
        "gbwt_ri": dataset_dir / "fastlocate.ri",
        "t1": dataset_dir / "output.t1",
        "t2": dataset_dir / "output.t2",
    }
    missing = [str(p) for p in paths.values() if not p.exists()]
    assert not missing, f"Missing required dataset files: {missing}"
    return paths


def main() -> None:
    args = parse_args()
    dataset_dir = Path(args.dataset).resolve()
    assert dataset_dir.exists(), f"Dataset directory does not exist: {dataset_dir}"
    paths = build_paths(dataset_dir)

    # 1) Error path: translate() before load() should fail.
    fresh_idx = liftover_ext.Index()
    expect_exception(
        lambda: fresh_idx.translate("_gbwt_ref#4294967295", 0, 10, "_gbwt_ref#4294967295"),
        RuntimeError,
        "translate before load",
    )

    # 2) Happy path: load index and inspect haplotypes.
    idx = liftover_ext.Index()
    idx.load(
        str(paths["gbz"]),
        str(paths["ri"]),
        str(paths["tags"]),
        str(paths["gbwt_ri"]),
        str(paths["t1"]),
        str(paths["t2"]),
    )
    haplotypes = idx.get_haplotype_names()
    assert isinstance(haplotypes, list), "get_haplotype_names() must return list"
    assert len(haplotypes) > 0, "Expected at least one haplotype after load()"
    src = haplotypes[0]
    tgt = haplotypes[0]

    # 3) Input validation paths.
    expect_exception(
        lambda: idx.translate(src, 100, 99, tgt),
        ValueError,
        "invalid interval start > end",
    )
    expect_exception(
        lambda: idx.translate(src, 0, 10_000_001, tgt),
        ValueError,
        "interval too long",
    )
    expect_exception(
        lambda: idx.translate("__missing_source__", 0, 10, tgt),
        ValueError,
        "missing source haplotype",
    )

    # 4) Translation API behavior for a handful of intervals.
    # Result cardinality depends on dataset and may be 0. We assert type/shape.
    query_intervals = [(0, 10), (0, 100), (100, 200), (0, 1000)]
    total_results = 0
    for start, end in query_intervals:
        result = idx.translate(src, start, end, tgt)
        assert isinstance(result, list), "translate() must return a list"
        total_results += len(result)
        for rec in result:
            assert isinstance(rec.haplotype, str), "TranslatedInterval.haplotype must be str"
            assert isinstance(rec.start, int), "TranslatedInterval.start must be int"
            assert isinstance(rec.end, int), "TranslatedInterval.end must be int"
            assert rec.strand in {"+", "-"}, "TranslatedInterval.strand must be '+' or '-'"

    print("PASS: Python coordinate translation integration test")
    print(f"dataset={dataset_dir}")
    print(f"haplotypes={len(haplotypes)} first={haplotypes[:3]}")
    print(f"intervals_tested={len(query_intervals)} translated_records={total_results}")


if __name__ == "__main__":
    main()

