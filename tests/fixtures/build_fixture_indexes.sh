#!/usr/bin/env bash
#
# Build the full index set for a test-fixture graph.
#
# The pytest suite (tests/python/) and the middleware both consume a fixed set
# of index files with canonical names. This script builds them from a single
# .gbz so a small real graph (e.g. chrM or chrY) can be turned into a complete,
# self-contained test fixture with one command.
#
# Usage:
#   tests/fixtures/build_fixture_indexes.sh <graph.gbz> <output_dir> [--coord-only]
#
# Environment:
#   BIN             pangenome-index bin dir (build_tags, build_rindex, ...).
#                   Default: <repo>/bin
#   VG              path to the `vg` binary. Default: `vg` on PATH.
#   THREADS         thread count for the parallel steps. Default: 4.
#   MINIMIZER_ARGS  extra args to `vg minimizer` (e.g. "-k 29 -w 11" for the
#                   short-read preset). Default: vg's own defaults.
#
# Produces, in <output_dir>:
#   graph.gbz          copy of the input graph
#   index.dist         distance index            (giraffe)      [skipped with --coord-only]
#   index.min          minimizer index           (giraffe)      [skipped with --coord-only]
#   index.zipcodes     zipcode index             (giraffe)      [skipped with --coord-only]
#   rlbwt_rindex.ri    RLBWT r-index             (coordinate translation)
#   sampled.tags       sampled tag array         (coordinate translation)
#   fastlocate.ri      GBWT FastLocate r-index   (coordinate translation)
#   output.t1          translation table 1
#   output.t2          translation table 2
#
# The coordinate chain (ri/tags/fastlocate/t1/t2) is enough for the translation
# tests. The giraffe chain (dist/min/zipcodes) is additionally needed for the
# mapping + surjection tests; pass --coord-only to skip it.
#
# NOTE on the minimizer preset: the minimizer index must be built with the same
# k/w the mapper expects. The middleware documents "pass the LONG-READ minimizer
# + zipcodes." If your production graphs use a specific preset, set MINIMIZER_ARGS
# to match; otherwise vg's defaults are used (fine for the small fixture, where
# the tests only require that sampled reads map back).

set -euo pipefail

if [[ $# -lt 2 ]]; then
    grep '^#' "$0" | sed 's/^# \{0,1\}//'   # print this header as help
    exit 2
fi

GBZ_IN="$1"
OUT="$2"
COORD_ONLY=0
[[ "${3:-}" == "--coord-only" ]] && COORD_ONLY=1

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BIN="${BIN:-$REPO_ROOT/bin}"
VG="${VG:-vg}"
THREADS="${THREADS:-4}"
MINIMIZER_ARGS="${MINIMIZER_ARGS:-}"

[[ -s "$GBZ_IN" ]] || { echo "ERROR: input graph not found: $GBZ_IN" >&2; exit 1; }
mkdir -p "$OUT"
TMP="$OUT/tmp"; mkdir -p "$TMP"

# Tool presence check (fail early with a clear message).
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found on PATH" >&2; exit 1; }; }
for t in "$BIN/build_tags" "$BIN/build_rindex" "$BIN/convert_tags" \
         "$BIN/build_sampled_tags" "$BIN/build_translation_tables"; do
    [[ -x "$t" ]] || { echo "ERROR: missing pangenome-index tool: $t (set BIN=)" >&2; exit 1; }
done
need gbz_extract
need grlbwt-cli
if [[ $COORD_ONLY -eq 0 ]]; then need "$VG"; fi

echo "=== [0/8] stage graph.gbz ==="
GBZ="$OUT/graph.gbz"
cp -f "$GBZ_IN" "$GBZ"

echo "=== [1/8] gbz_extract -> RL-BWT info ==="
INFO="$OUT/graph_info"
gbz_extract -t "$THREADS" -b "$GBZ" > "$INFO"

echo "=== [2/8] grlbwt-cli -> ${INFO}.rl_bwt ==="
grlbwt-cli -t "$THREADS" -T "$TMP" "$INFO"
RLBWT="${INFO}.rl_bwt"
[[ -s "$RLBWT" ]] || { echo "ERROR: grlbwt-cli did not produce $RLBWT" >&2; exit 1; }

echo "=== [3/8] build_tags ==="
"$BIN/build_tags" "$GBZ" "$RLBWT" "$OUT/raw.tags"

echo "=== [4/8] build_rindex -> rlbwt_rindex.ri ==="
"$BIN/build_rindex" "$RLBWT" > "$OUT/rlbwt_rindex.ri"

echo "=== [5/8] convert_tags + build_sampled_tags -> sampled.tags ==="
"$BIN/convert_tags" "$OUT/raw.tags" "$OUT/compact.tags"
"$BIN/build_sampled_tags" "$OUT/compact.tags" "$OUT/sampled.tags"

echo "=== [6/8] vg gbwt -> fastlocate.ri (GBWT FastLocate) ==="
# The GBWT FastLocate r-index is built from the GBZ's embedded GBWT.
# (Flag names are stable across recent vg; adjust if your vg differs.)
"${VG:-vg}" gbwt -Z "$GBZ" -r "$OUT/fastlocate.ri"

echo "=== [7/8] build_translation_tables -> output.t1 / output.t2 ==="
"$BIN/build_translation_tables" "$GBZ" "$OUT/sampled.tags" \
    --table1 "$OUT/output.t1" --table2 "$OUT/output.t2"

if [[ $COORD_ONLY -eq 1 ]]; then
    echo "=== done (coordinate indexes only) ==="
    ls -la "$OUT"
    exit 0
fi

echo "=== [8/8] giraffe indexes: distance / minimizer / zipcodes ==="
"$VG" index -t "$THREADS" -j "$OUT/index.dist" "$GBZ"
# shellcheck disable=SC2086
"$VG" minimizer -t "$THREADS" -d "$OUT/index.dist" -z "$OUT/index.zipcodes" \
    $MINIMIZER_ARGS -o "$OUT/index.min" "$GBZ"

echo "=== done ==="
ls -la "$OUT"
