# Test suite

Correctness tests for the pangenome mapping tool: coordinate translation,
anchor building, surjection, and the mapping/haplotype pipeline.

Every test asserts against an **independent oracle** (see the test plan), not
just "it ran": the haplotype sequence for translation, a known sampled-read
locus for mapping/surjection, stock `vg surject` for the differential check, and
round-trip invariants. Tests that only need code (no index data) are pure unit
tests; the rest run against a *fixture* graph with its indexes built.

## Layout

| Where | What | Harness |
|-------|------|---------|
| `Giraffe_server/src/unittest/` | C++ component tests (anchor-backed graph, surjector, haplotype assigner) | Catch2, run via `vg test` |
| `tests/*.cpp` | C++ tests for the index tools (anchor builder, r-index, b+ tree) | standalone, `tests/makefile` |
| `tests/python/` | Python tests (translation, middleware multiplexer, …) | pytest |
| `tests/integration/` | standalone end-to-end / benchmark drivers (run manually with a real graph) | `python …` |
| `tests/fixtures/` | fixture-builder script + mock giraffe-server; built fixtures land in ignored subdirs here | — |

## C++ component tests (Catch2)

These build into the `vg` binary and run through its test subcommand. After
building `vg` (see `Giraffe_server/README.md`):

```bash
# all anchor-backed position-graph tests, including the two regression cases
vg test "[anchor_backed_graph]"

# just the regressions (O(gap) walk; target-only for_each_step_on_handle)
vg test "[regression]"
```

Regression coverage added here:
- `get_position_of_step` walks **O(gap)** backward to the nearest known step and
  backfills — not O(path) forward from the start (the old O(path²) behavior).
- `for_each_step_on_handle` returns **only the target path's** steps over the
  read region — not every haplotype's steps on the node (the 12 s surjection).

## Python tests (pytest)

### The extension must match your Python

`liftover_ext` is a native (pybind11) extension; import fails with an
"incompatible architecture" / ABI error unless pytest runs under the *exact*
Python the `.so` was built for. On the dev Mac that is an **arm64 Python 3.12**
(Homebrew), not the x86_64 conda Python:

```bash
/opt/homebrew/bin/python3.12 -m venv .venv-test312
.venv-test312/bin/python -m pip install pytest
.venv-test312/bin/python -m pytest tests/python -q
```

On the server, use the same interpreter the service runs under (the one
`liftover_ext` was built against).

### Fixtures

Tests resolve a *coordinate fixture* — a directory with `graph.gbz` plus its
coordinate indexes (`rlbwt_rindex.ri`, `sampled.tags`, `fastlocate.ri`,
`output.t1`, `output.t2`).

- **Default (no setup):** the tiny in-repo `coord_translation_tests/test0_target_loop`.
  It is degenerate (one generic `_gbwt_ref` path), so the API/validation
  **contract** tests run and the `@pytest.mark.oracle` tests self-skip.
- **Real fixture:** build one and point the suite at it:

  ```bash
  # needs: vg, gbz_extract, grlbwt-cli on PATH; pangenome-index bins in ./bin
  tests/fixtures/build_fixture_indexes.sh test_data/chrM/chrM.gbz tests/fixtures/chrM
  PANGENOME_TEST_FIXTURE=tests/fixtures/chrM .venv-test312/bin/python -m pytest tests/python -q
  ```

  Pass `--coord-only` to `build_fixture_indexes.sh` to skip the giraffe indexes
  (dist/min/zipcodes) when you only need translation tests.

### Markers

- `oracle` — needs a real multi-haplotype fixture with homology; self-skips otherwise.
- `needs_giraffe` — needs the giraffe indexes + a `vg` binary.
- `slow` — large-graph / long-running.

```bash
.venv-test312/bin/python -m pytest tests/python -m "not slow" -q
```

## Status

Implemented:
- C++ `AnchorBackedPositionGraph` unit + regression tests.
- Python coordinate-translation contract tests + oracle scaffold (round-trip,
  identity), fixture resolution, and the fixture builder.

Next (see the test plan): sampled-read → surject oracle + stock-`vg` differential,
haplotype-assigner unit tests, anchor-builder unit tests, the engineered tiny
fixture, index-construction tests, and middleware/API completion.
