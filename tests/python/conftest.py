"""Shared pytest fixtures for the pangenome-index Python test suite.

The suite is built around a *coordinate fixture*: a directory holding a graph
(`graph.gbz`) and its coordinate indexes (`rlbwt_rindex.ri`, `sampled.tags`,
`fastlocate.ri`, `output.t1`, `output.t2`). A loaded `liftover_ext.Index` over
that fixture is exposed as the `coord_index` fixture.

By default the tiny in-repo fixture ``coord_translation_tests/test0_target_loop``
is used. It is *degenerate* (a single generic ``_gbwt_ref`` path) — enough for
the API / validation contract tests, but its self-translation is empty, so the
``@pytest.mark.oracle`` tests self-skip on it. Point

    PANGENOME_TEST_FIXTURE=/path/to/fixture_dir

at a real multi-haplotype fixture (e.g. one produced by
``tests/fixtures/build_fixture_indexes.sh`` on chrM) to exercise the oracle tests
too.

The ``liftover_ext`` extension is arch/ABI-specific: it must be built for the
exact Python running pytest. On this repo's dev Mac that means an arm64 Python
3.12 (Homebrew), not the x86_64 conda Python. See tests/README.md.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import shutil
import subprocess
from typing import Optional

import pytest

# Repo root = two levels up from this file (tests/python/conftest.py). Put it on
# sys.path so `import liftover_ext` (built at the repo root) and
# `from middleware... import ...` both resolve.
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# Canonical coordinate-index file names inside a fixture directory. These match
# both the in-repo fixtures and the output of build_fixture_indexes.sh.
COORD_FILES = {
    "gbz": "graph.gbz",
    "ri": "rlbwt_rindex.ri",
    "tags": "sampled.tags",
    "gbwt_ri": "fastlocate.ri",
    "t1": "output.t1",
    "t2": "output.t2",
}

DEFAULT_FIXTURE = REPO_ROOT / "coord_translation_tests" / "test0_target_loop"

# The committed small real graph (34 K, 44 haplotypes, reference + accession
# contig names). Its coordinate indexes are built ON DEMAND into this ignored
# cache — never committed (they're large and index-format-version coupled).
CHRM_GRAPH = REPO_ROOT / "test_data" / "chrM" / "chrM.gbz"
CHRM_CACHE = REPO_ROOT / "tests" / "fixtures" / "chrM"        # ignored (tests/fixtures/*/)
FIXTURE_BUILDER = REPO_ROOT / "tests" / "fixtures" / "build_fixture_indexes.sh"


def _cache_ready(d: Path) -> bool:
    """True if every coordinate index is present and newer than the graph."""
    if not d.is_dir():
        return False
    files = [d / name for name in COORD_FILES.values()]
    if not all(f.exists() for f in files):
        return False
    graph_mtime = CHRM_GRAPH.stat().st_mtime if CHRM_GRAPH.exists() else 0.0
    return min(f.stat().st_mtime for f in files) >= graph_mtime


def _build_tools_available() -> bool:
    vg = os.environ.get("VG") or shutil.which("vg")
    have_rlbwt = shutil.which("gbz_extract") and shutil.which("grlbwt-cli")
    have_bins = (REPO_ROOT / "bin" / "build_tags").exists()
    return bool(vg and have_rlbwt and have_bins
                and FIXTURE_BUILDER.exists() and CHRM_GRAPH.exists())


def _ensure_chrm_built() -> Optional[Path]:
    """Build the chrM coordinate indexes once into the ignored cache and return
    it, or None if the build toolchain isn't available (so callers fall back)."""
    if _cache_ready(CHRM_CACHE):
        return CHRM_CACHE
    if not _build_tools_available():
        return None
    env = dict(os.environ)
    env.setdefault("BIN", str(REPO_ROOT / "bin"))
    env.setdefault("VG", shutil.which("vg") or "")
    CHRM_CACHE.mkdir(parents=True, exist_ok=True)
    try:
        subprocess.run(
            ["bash", str(FIXTURE_BUILDER), str(CHRM_GRAPH), str(CHRM_CACHE), "--coord-only"],
            check=True, env=env, timeout=1800,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
    except Exception:
        return None
    return CHRM_CACHE if _cache_ready(CHRM_CACHE) else None


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "oracle: needs a real multi-haplotype fixture with actual homology; "
        "self-skips on the degenerate default fixture",
    )
    config.addinivalue_line(
        "markers", "needs_giraffe: needs the giraffe indexes (dist/min/zipcodes) + a vg binary"
    )
    config.addinivalue_line("markers", "slow: slow / large-graph test")


def _fixture_dir() -> Path:
    # Explicit override wins.
    env = os.environ.get("PANGENOME_TEST_FIXTURE")
    if env:
        return Path(env).expanduser().resolve()
    # Otherwise use the real chrM fixture when the build toolchain is present
    # (its indexes are built once into the ignored cache); fall back to the
    # tiny degenerate in-repo fixture so contract tests still run everywhere.
    built = _ensure_chrm_built()
    return built if built is not None else DEFAULT_FIXTURE


@pytest.fixture(scope="session")
def liftover():
    """The native coordinate-translation extension, or skip if it can't load."""
    try:
        import liftover_ext  # noqa: WPS433 (import inside fixture is intentional)
    except ImportError as exc:
        pytest.skip(
            f"liftover_ext not importable ({exc}); build the extension for THIS "
            f"Python ({sys.executable}, {sys.version.split()[0]})"
        )
    return liftover_ext


@pytest.fixture(scope="session")
def fixture_dir() -> Path:
    d = _fixture_dir()
    if not d.is_dir():
        pytest.skip(f"fixture dir not found: {d} (set PANGENOME_TEST_FIXTURE)")
    return d


@pytest.fixture(scope="session")
def fixture_paths(fixture_dir):
    paths = {key: fixture_dir / name for key, name in COORD_FILES.items()}
    missing = [str(p) for p in paths.values() if not p.exists()]
    if missing:
        pytest.skip(f"fixture {fixture_dir} is missing coord files: {missing}")
    return paths


@pytest.fixture(scope="session")
def coord_index(liftover, fixture_paths):
    """A loaded liftover_ext.Index over the configured fixture (loaded once)."""
    idx = liftover.Index()
    idx.load(
        str(fixture_paths["gbz"]),
        str(fixture_paths["ri"]),
        str(fixture_paths["tags"]),
        str(fixture_paths["gbwt_ri"]),
        str(fixture_paths["t1"]),
        str(fixture_paths["t2"]),
    )
    return idx


@pytest.fixture(scope="session")
def haplotype_names(coord_index):
    names = coord_index.get_haplotype_names()
    assert isinstance(names, list) and names, "fixture reports no haplotypes"
    return names
