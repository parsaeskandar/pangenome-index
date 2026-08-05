"""Concurrency tests for the multiplexed GiraffeServerMiddleware.

These use a mock `vg giraffe-server` (tests/fixtures/mock_giraffe_server.py) so
they exercise the dispatcher + demuxer — coalescing, per-name routing, no
cross-talk between concurrent callers — without needing a real vg binary or any
indexes. This is where the concurrent code is validated; the real stack is
validated separately on the server.
"""
from __future__ import annotations

import os
import stat
import sys
import threading
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from middleware.giraffe_server_middleware import (  # noqa: E402
    GiraffeServerConfig,
    GiraffeServerMiddleware,
)

MOCK = REPO_ROOT / "tests" / "fixtures" / "mock_giraffe_server.py"


class FakeAnchor:
    """Minimal stand-in with the fields surject_with_anchors serializes."""
    gbwt_edge_begin_node = 1
    gbwt_edge_begin_offset = 0
    gbwt_edge_end_node = 2
    gbwt_edge_end_offset = 0
    path_offset_step_begin = 0
    path_offset_step_end = 1
    read_begin_offset = 0
    read_end_offset = 1
    source_mapping_begin = 0
    source_mapping_end = 1


@pytest.fixture
def mw(tmp_path):
    # Make the mock executable so Popen can run it as the "vg" binary.
    MOCK.chmod(MOCK.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    # Dummy index files so start()'s existence checks pass.
    paths = {}
    for name in ("gbz", "min", "dist", "zip"):
        p = tmp_path / f"{name}.idx"
        p.write_text("x")
        paths[name] = str(p)
    cfg = GiraffeServerConfig(
        vg_binary=str(MOCK),
        gbz_path=paths["gbz"],
        minimizer_path=paths["min"],
        distance_path=paths["dist"],
        zipcode_path=paths["zip"],
        threads=4,
        batch_size=256,
        output_timeout_s=30.0,
    )
    m = GiraffeServerMiddleware(cfg)
    m.start()
    m.wait_until_ready()
    yield m
    m.stop()


def test_single_batch_roundtrips_in_order(mw):
    reads = [(f"r{i}", f"AAAA{i}CCCC", "I" * len(f"AAAA{i}CCCC")) for i in range(6)]
    out = mw.map_reads(reads)
    assert len(out) == 6
    for i in range(6):
        assert len(out[i]) == 1
        assert f"gaf:AAAA{i}CCCC" in out[i][0]


def test_concurrent_map_has_no_crosstalk(mw):
    """Many threads mapping at once: each must get back exactly its own reads."""
    n_threads, n_reads = 12, 8
    results = {}
    errors = []

    def worker(t):
        try:
            reads = [(f"t{t}r{i}", f"SEQ{t}_{i}", "I" * len(f"SEQ{t}_{i}"))
                     for i in range(n_reads)]
            results[t] = mw.map_reads(reads)
        except Exception as exc:  # noqa: BLE001
            errors.append((t, repr(exc)))

    threads = [threading.Thread(target=worker, args=(t,)) for t in range(n_threads)]
    for th in threads:
        th.start()
    for th in threads:
        th.join()

    assert not errors, errors
    for t in range(n_threads):
        out = results[t]
        assert len(out) == n_reads, (t, len(out))
        for i in range(n_reads):
            assert len(out[i]) == 1
            # Read i of thread t must carry thread t's marker — never another's.
            assert f"gaf:SEQ{t}_{i}" in out[i][0], (t, i, out[i])


def test_concurrent_surject_routes_by_name(mw):
    n_threads = 12
    results = {}
    errors = []

    def worker(t):
        try:
            gaf = f"readT{t}\t100\t0\t100\t+\tpathX\t200\t0\t100\t100\t100\t60"
            results[t] = mw.surject_with_anchors(
                gaf, [FakeAnchor()], "TARGET#0#chr1", target_path_length=200)
        except Exception as exc:  # noqa: BLE001
            errors.append((t, repr(exc)))

    threads = [threading.Thread(target=worker, args=(t,)) for t in range(n_threads)]
    for th in threads:
        th.start()
    for th in threads:
        th.join()

    assert not errors, errors
    for t in range(n_threads):
        out = results[t]
        assert out, (t, out)
        assert any(f"surj:readT{t}" in line for line in out), (t, out)


def test_mixed_map_and_surject_concurrently(mw):
    """Map workers and surject workers hitting the one subprocess together."""
    results = {}
    errors = []

    def map_worker(t):
        try:
            reads = [(f"m{t}_{i}", f"MAP{t}_{i}", "I" * len(f"MAP{t}_{i}")) for i in range(4)]
            results[("map", t)] = mw.map_reads(reads)
        except Exception as exc:  # noqa: BLE001
            errors.append(("map", t, repr(exc)))

    def surj_worker(t):
        try:
            gaf = f"sread{t}\t50\t0\t50\t+\tp\t100\t0\t50\t50\t50\t60"
            results[("surj", t)] = mw.surject_with_anchors(gaf, [FakeAnchor()], "TGT")
        except Exception as exc:  # noqa: BLE001
            errors.append(("surj", t, repr(exc)))

    threads = []
    for t in range(6):
        threads.append(threading.Thread(target=map_worker, args=(t,)))
        threads.append(threading.Thread(target=surj_worker, args=(t,)))
    for th in threads:
        th.start()
    for th in threads:
        th.join()

    assert not errors, errors
    for t in range(6):
        m = results[("map", t)]
        assert len(m) == 4
        for i in range(4):
            assert f"gaf:MAP{t}_{i}" in m[i][0]
        s = results[("surj", t)]
        assert any(f"surj:sread{t}" in line for line in s)
