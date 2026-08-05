#!/bin/bash
#SBATCH --job-name=chr_tags
#SBATCH --partition=long
#SBATCH --mail-user=seeskand@ucsc.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=300gb
#SBATCH --time=24:00:00
#SBATCH --output=%x.%A_%a.log
#SBATCH --no-requeue
#SBATCH --array=1-30            # override on the command line: --array=1-<#chunks>
#
# Stage 2 of the per-chromosome tag pipeline.
# One array task = one chunk from make_gbz_chunks.sh. Each task builds that
# chunk's r-index and tag array; all per-chunk .tags land in one TAGS_DIR that
# stage 3 (merge_chr_tags.sh) then merges into the whole-graph sampled.tags.
#
# Per-chunk sizes are far below the 2^32 counters and the memory that made the
# whole-graph build_tags OOM/segfault, so this is the scalable route.
#
set -euo pipefail
source ~/.bashrc

WORK="${WORK:-/private/groups/cgl/seeskand/graph_index/13.shloka_graph}"
CHUNKDIR="${CHUNKDIR:-$WORK/chromosome_tags}"
MANIFEST="${MANIFEST:-$CHUNKDIR/chunk_manifest.txt}"
TAGS_DIR="${TAGS_DIR:-$CHUNKDIR/tags_dir}"

VG="${VG:-$HOME/software/vg_binaries/vg}"   # the vg used to build the working d46 chunks
GBZ_EXTRACT="${GBZ_EXTRACT:-/private/groups/cgl/seeskand/1.server/Giraffe_server/deps/gbwtgraph/bin/gbz_extract}"
GRLBWT_CLI="${GRLBWT_CLI:-/private/groups/cgl/seeskand/software/grlBWT/build/grlbwt-cli}"
BIN="${BIN:-/private/groups/cgl/seeskand/pangenome-index/bin}"
T="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "$TAGS_DIR"

CHUNK="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST" || true)"
if [[ -z "${CHUNK}" ]]; then
  echo "[task ${SLURM_ARRAY_TASK_ID}] no manifest line ${SLURM_ARRAY_TASK_ID}; nothing to do."
  exit 0
fi

name="$(basename "$CHUNK")"; name="${name%.gfa}"; name="${name%.gbz}"; name="${name%.pg}"; name="${name%.vg}"
d="$CHUNKDIR/work_$name"; mkdir -p "$d/tmp"; cd "$d"
echo "[task ${SLURM_ARRAY_TASK_ID}] chunk=$CHUNK name=$name pwd=$PWD"

# Resumable: skip if this chunk's tags are already done.
if [[ -s "$TAGS_DIR/$name.tags" ]]; then
  echo "[task] tags already exist, skipping: $TAGS_DIR/$name.tags"
  exit 0
fi

# Ensure a GBZ (WITH the full haplotype set) for this chunk.
#   .gfa  -> rebuild with `vg gbwt -G --max-node 0` (chunk_xg_gfa.sh output).
#            The GFA (from the XG chunking) carries every haplotype as a path,
#            so the resulting GBZ has real sequences; --max-node 0 keeps the
#            original node ids so the tags line up with the whole-graph r-index.
#   .gbz  -> use directly.
#   .vg/.pg -> legacy fallback (embedded-path extraction; can drop haplotypes).
if [[ "$CHUNK" == *.gbz ]]; then
  CGBZ="$CHUNK"
elif [[ "$CHUNK" == *.gfa ]]; then
  CGBZ="$d/$name.gbz"
  [[ -s "$CGBZ" ]] || /usr/bin/time -v "$VG" gbwt -p -G "$CHUNK" --max-node 0 -g "$CGBZ"
else
  CGBZ="$d/$name.gbz"
  [[ -s "$CGBZ" ]] || /usr/bin/time -v "$VG" gbwt -p -x "$CHUNK" -g "$CGBZ" --gbz-format -E
fi

# Per-chunk RL-BWT (gbz_extract -> grlbwt-cli).
info="$d/${name}_info"
if [[ ! -s "${info}.rl_bwt" ]]; then
  /usr/bin/time -v "$GBZ_EXTRACT" -t "$T" -b -p "$CGBZ" > "${info}.tmp"
  mv "${info}.tmp" "$info"
  # A degenerate connected component (no haplotype threads) extracts to a
  # (near-)empty info; grlbwt then reads 0 strings and aborts with a SIZE_MAX
  # symbol count. Such a component has no sequence to tag, so skip it — but log
  # loudly with the size so a wrongly-empty BIG chunk (a real conversion
  # failure, not a degenerate component) is obvious in the logs.
  info_bytes="$(stat -c%s "$info" 2>/dev/null || echo 0)"
  if [[ "${info_bytes}" -lt 1000 ]]; then
    echo "[task] SKIP: chunk '$name' RL-BWT input is empty/degenerate (${info_bytes} bytes)" \
         "— no haplotypes in this component to tag."
    exit 0
  fi
  /usr/bin/time -v "$GRLBWT_CLI" -t "$T" -T "$d/tmp" "$info"
fi

# Per-chunk tags -> shared TAGS_DIR (atomic temp+mv so a killed task never
# leaves a half-written .tags that merge_tags would pick up).
/usr/bin/time -v "$BIN/build_tags" -t "$T" "$CGBZ" "${info}.rl_bwt" "$TAGS_DIR/$name.tags.tmp"
mv "$TAGS_DIR/$name.tags.tmp" "$TAGS_DIR/$name.tags"
echo "[task] done: $TAGS_DIR/$name.tags"
