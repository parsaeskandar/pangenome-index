#!/bin/bash
#SBATCH --job-name=merge_chr_tags
#SBATCH --partition=long
#SBATCH --mail-user=seeskand@ucsc.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=700gb
#SBATCH --time=24:00:00
#SBATCH --output=%x.%j.log
#SBATCH --no-requeue
#
# Stage 3 of the per-chromosome tag pipeline.
# Merge every per-chunk tag array (built by build_chr_tags_array.sh) into the
# whole-graph compressed tag array using the WHOLE-graph gbz + r-index, then
# build the sampled tags the coordinate translator / server load.
#
# merge_tags derives the per-chunk BWT offsets from the whole-graph gbz + r-index
# (it does NOT rely on chunk file order), so TAGS_DIR must contain exactly the
# per-chromosome .tags files and nothing else.
#
set -euo pipefail
source ~/.bashrc

WORK="${WORK:-/private/groups/cgl/seeskand/graph_index/13.shloka_graph}"
GBZ="${GBZ:-$WORK/E821-16-sampled.gbz}"
RINDEX="${RINDEX:-$WORK/rlbwt_rindex.ri}"          # whole-graph build_rindex output
CHUNKDIR="${CHUNKDIR:-$WORK/chromosome_tags}"
TAGS_DIR="${TAGS_DIR:-$CHUNKDIR/tags_dir}"
BIN="${BIN:-/private/groups/cgl/seeskand/pangenome-index/bin}"
# --in-memory is fast but RAM-heavy; set IN_MEMORY="" to stream if it OOMs.
IN_MEMORY="${IN_MEMORY:---in-memory}"

for f in "$GBZ" "$RINDEX"; do
  [[ -s "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done
if ! ls "$TAGS_DIR"/*.tags >/dev/null 2>&1; then
  echo "ERROR: no .tags files in $TAGS_DIR" >&2; exit 1
fi

cd "$WORK"
# merge_tags opens its output in append mode; clear any partial prior run.
rm -f "$WORK"/whole_genome_tag_array_compressed.tags*

echo "Per-chromosome tags to merge:"
ls -la "$TAGS_DIR"/*.tags

echo "=== merge_tags ==="
/usr/bin/time -v "$BIN/merge_tags" $IN_MEMORY "$GBZ" "$RINDEX" "$TAGS_DIR"
# -> writes ./whole_genome_tag_array_compressed.tags in $WORK

echo "=== build_sampled_tags ==="
/usr/bin/time -v "$BIN/build_sampled_tags" \
    "$WORK/whole_genome_tag_array_compressed.tags" "$WORK/sampled.tags"

echo "=== done ==="
ls -la "$WORK/sampled.tags"
