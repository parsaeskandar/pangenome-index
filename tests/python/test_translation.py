"""Coordinate-translation tests (liftover_ext.Index).

Two tiers:

* **contract** tests run on any loadable fixture (including the degenerate
  in-repo default). They pin the API shape, input validation, and determinism —
  the things that must hold regardless of graph content.
* **oracle** tests (``@pytest.mark.oracle``) need a real multi-haplotype fixture
  with actual homology. They discover a translatable interval at runtime and
  self-skip if the fixture yields none, so they are harmless on the degenerate
  default and meaningful on a chrM/chrY fixture.

Run:  pytest tests/python -q
Oracle:  PANGENOME_TEST_FIXTURE=/path/to/chrM_fixture pytest tests/python -q
"""
from __future__ import annotations

import pytest


def _tuples(records):
    """Normalize TranslatedInterval records to comparable tuples."""
    return [(r.haplotype, r.start, r.end, r.strand) for r in records]


# --------------------------------------------------------------------------- #
# Contract: error / validation paths
# --------------------------------------------------------------------------- #

def test_translate_before_load_raises(liftover):
    """Using the index before load() is a usage error, not a crash."""
    idx = liftover.Index()
    with pytest.raises(RuntimeError):
        idx.translate("_gbwt_ref#4294967295", 0, 10, "_gbwt_ref#4294967295")


def test_get_haplotype_names_is_nonempty_str_list(haplotype_names):
    assert isinstance(haplotype_names, list)
    assert haplotype_names
    assert all(isinstance(n, str) and n for n in haplotype_names)


def test_get_haplotype_names_is_deterministic(coord_index):
    assert coord_index.get_haplotype_names() == coord_index.get_haplotype_names()


def test_translate_returns_list_of_typed_intervals(coord_index, haplotype_names):
    src = haplotype_names[0]
    result = coord_index.translate(src, 0, 50, src)
    assert isinstance(result, list)
    for rec in result:
        assert isinstance(rec.haplotype, str) and rec.haplotype
        assert isinstance(rec.start, int)
        assert isinstance(rec.end, int)
        assert rec.start <= rec.end
        assert rec.strand in ("+", "-")


@pytest.mark.parametrize("start,end", [(100, 99), (50, 0)])
def test_translate_rejects_inverted_interval(coord_index, haplotype_names, start, end):
    with pytest.raises(ValueError):
        coord_index.translate(haplotype_names[0], start, end, haplotype_names[0])


def test_translate_rejects_overlong_interval(coord_index, haplotype_names):
    # The binding rejects intervals longer than MAX_INTERVAL_LENGTH (1e7).
    with pytest.raises(ValueError):
        coord_index.translate(haplotype_names[0], 0, 10_000_001, haplotype_names[0])


def test_translate_rejects_missing_source(coord_index, haplotype_names):
    with pytest.raises(ValueError):
        coord_index.translate("__no_such_haplotype__", 0, 10, haplotype_names[0])


def test_translate_is_deterministic(coord_index, haplotype_names):
    src = haplotype_names[0]
    first = coord_index.translate(src, 0, 100, src)
    second = coord_index.translate(src, 0, 100, src)
    assert _tuples(first) == _tuples(second)


# --------------------------------------------------------------------------- #
# Oracle: real homology required (self-skips on the degenerate fixture)
# --------------------------------------------------------------------------- #

def _first_translatable(idx, src, tgt, span=200, step=500, max_offset=200_000):
    """Scan offsets on ``src`` for the first [off, off+span) that yields >=1
    translation to ``tgt``. Returns (start, end, records) or None."""
    off = 0
    while off <= max_offset:
        try:
            res = idx.translate(src, off, off + span, tgt)
        except ValueError:
            res = []
        if res:
            return off, off + span, res
        off += step
    return None


@pytest.mark.oracle
def test_identity_translation_preserves_interval(coord_index, haplotype_names):
    """translate(A, i, j, A) must return the same interval on A, forward strand.

    Self-oracle: the source haplotype's own sequence is ground truth, so an
    identity lift is exact by construction.
    """
    src = haplotype_names[0]
    found = _first_translatable(coord_index, src, src)
    if not found:
        pytest.skip(f"no self-translatable interval on {src!r} (degenerate fixture?)")
    start, end, res = found
    same = [r for r in res if r.haplotype == src]
    assert same, f"identity produced no record on {src!r}: {_tuples(res)}"
    assert any(r.start == start and r.end == end and r.strand == "+" for r in same), (
        f"identity did not preserve [{start},{end}) on {src!r}: "
        f"{[(r.start, r.end, r.strand) for r in same]}"
    )


@pytest.mark.oracle
def test_translation_round_trips(coord_index, haplotype_names):
    """translate(A->B) then translate(B->A) must recover the original interval
    (exact for SNP-only regions, within an indel-sized tolerance elsewhere)."""
    if len(haplotype_names) < 2:
        pytest.skip("round-trip needs >=2 haplotypes")

    src = haplotype_names[0]
    found = None
    for tgt in haplotype_names[1:]:
        found = _first_translatable(coord_index, src, tgt)
        if found:
            break
    if not found:
        pytest.skip(f"no cross-haplotype homology found from {src!r}")

    start, end, forward = found
    tolerance = 50  # allow for indels between the two haplotypes
    recovered = False
    for fwd in forward:
        back = coord_index.translate(fwd.haplotype, fwd.start, fwd.end, src)
        if any(
            b.haplotype == src
            and abs(b.start - start) <= tolerance
            and abs(b.end - end) <= tolerance
            for b in back
        ):
            recovered = True
            break
    assert recovered, (
        f"round-trip {src} -> {forward[0].haplotype} -> {src} did not recover "
        f"[{start},{end}) within {tolerance}bp"
    )
