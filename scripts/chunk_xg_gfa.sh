#!/bin/bash
#SBATCH --job-name=chunk_xg_gfa
#SBATCH --partition=long
#SBATCH --mail-user=seeskand@ucsc.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=800gb
#SBATCH --time=48:00:00
#SBATCH --output=%x.%j.log
#SBATCH --no-requeue
#
# Stage 1 of the per-chromosome tag pipeline (correct chunking).
#
#   GBZ --(vg convert -x)--> XG --(vg chunk -C -O gfa)--> one GFA per component
#
# Converting to XG first turns every GBWT haplotype thread into a path, so each
# per-component GFA carries the FULL haplotype set. (Chunking the GBZ directly
# with `vg chunk -C` yields plain graphs that drop the haplotypes -> the
# per-chunk gbz comes out empty -> grlbwt crashes / tags are wrong.)
#
# build_chr_tags_array.sh then rebuilds each chunk's GBZ from its GFA with
# `vg gbwt -G <gfa> --max-node 0 -g <gbz>` (original node ids preserved so the
# per-chunk tags line up with the whole-graph r-index at merge time).
#
set -euo pipefail
source ~/.bashrc

WORK="${WORK:-/private/groups/cgl/seeskand/graph_index/13.shloka_graph}"
GBZ="${GBZ:-$WORK/E821-16-sampled.gbz}"
CHUNKDIR="${CHUNKDIR:-$WORK/chromosome_tags}"
# Use the vg that produced the working d46 chunks (has the convert/chunk/gbwt
# behavior we rely on); override with VG=... if it lives elsewhere.
VG="${VG:-$HOME/software/vg_binaries/vg}"
T="${SLURM_CPUS_PER_TASK:-32}"

[[ -s "$GBZ" ]] || { echo "ERROR: missing $GBZ" >&2; exit 1; }
mkdir -p "$CHUNKDIR"; cd "$CHUNKDIR"

# 1) GBZ -> XG (materializes all haplotype threads as paths).
if [[ ! -s "$CHUNKDIR/graph.xg" ]]; then
  /usr/bin/time -v "$VG" convert -x "$GBZ" -t "$T" > "$CHUNKDIR/graph.xg.tmp"
  mv "$CHUNKDIR/graph.xg.tmp" "$CHUNKDIR/graph.xg"
fi

# 2) Chunk by connected component -> one GFA per component.
/usr/bin/time -v "$VG" chunk -C -O gfa -t "$T" -x "$CHUNKDIR/graph.xg" -b "$CHUNKDIR/chunk"

# 3) Manifest of GFA chunks (absolute paths). A fresh CHUNKDIR contains only the
#    chunk GFAs plus graph.xg, so globbing *.gfa is safe.
ls -1 "$CHUNKDIR"/*.gfa 2>/dev/null | sort -V | xargs -r readlink -f > "$CHUNKDIR/chunk_manifest.txt"
n="$(wc -l < "$CHUNKDIR/chunk_manifest.txt" | tr -d ' ')"
if [[ "$n" -eq 0 ]]; then
  echo "ERROR: no *.gfa chunks produced in $CHUNKDIR" >&2
  echo "What vg chunk wrote:"; ls -la "$CHUNKDIR"
  exit 1
fi
echo "Wrote $n chunk GFA path(s) -> $CHUNKDIR/chunk_manifest.txt"
