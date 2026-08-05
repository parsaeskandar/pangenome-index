#!/usr/bin/env bash
#SBATCH --job-name=mk_gbz_chunks
#SBATCH --partition=long
#SBATCH --mail-user=seeskand@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=800gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:0
#SBATCH --output=%x.%j.log
#SBATCH --time=80:00:00
#SBATCH --no-requeue
#SBATCH --exclude=phoenix-00
#
# Create chunks from a whole-genome GBZ using a single vg chunk process (no vg loops).
#
# Defaults (HPRC d46 CHM13 + your index layout); override with MAIN_GBZ / OUT_DIR if needed:
#   Graph:  .../hprc-v2.0-mc-chm13.d46.gbz
#   Output: .../graph_index/12.hprc-v2.0-mc-chm13.d46/chromosome_tags
#
# Local:
#   bash scripts/make_gbz_chunks.sh
#
# Slurm (defaults baked in; only export VG_BIN / PATH if vg is not on the default PATH):
#   sbatch /path/to/pangenome-index/scripts/make_gbz_chunks.sh
#
# Optional env:
#   PREFIX     basename prefix for -b (default: chunks)
#   THREADS    vg -t (default: nproc or 16)
#   VG_BIN     default: vg
#   CONTIG     with --gbz only: passed as --contig "$CONTIG"
#
# Chunking (exactly one vg invocation):
#   - vg v1.73+ (e.g. "Ducky"): `vg chunk --help` lists --gbz → uses vg chunk --gbz (multiple .gbz).
#   - Older vg without --gbz: falls back to vg chunk -C (multiple .pg/.vg under OUT_DIR).
# Put vg on PATH before sbatch, or set VG_BIN to the full binary (avoid sourcing ~/.bashrc here).
#
set -euo pipefail

MAIN_GBZ="${MAIN_GBZ:-/private/groups/cgl/hprc-graphs/hprc-v2.0-feb28/hprc-v2.0-mc-chm13/hprc-v2.0-mc-chm13.d46.gbz}"
OUT_DIR="${OUT_DIR:-/private/groups/cgl/seeskand/graph_index/12.hprc-v2.0-mc-chm13.d46/chromosome_tags}"
if [[ ! -f "${MAIN_GBZ}" ]]; then
  echo "ERROR: MAIN_GBZ is not a file: ${MAIN_GBZ}" >&2
  exit 1
fi
PREFIX="${PREFIX:-chunks}"
VG_BIN="${VG_BIN:-vg}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || echo 16)}}"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

base="${PWD}/${PREFIX}"

if "${VG_BIN}" chunk --help 2>&1 | grep -q -- '--gbz'; then
  cmd=( "${VG_BIN}" chunk --gbz -x "${MAIN_GBZ}" -b "${base}" -t "${THREADS}" )
  if [[ -n "${CONTIG:-}" ]]; then
    cmd+=( --contig "${CONTIG}" )
  fi
  "${cmd[@]}"
else
  "${VG_BIN}" chunk -x "${MAIN_GBZ}" -C -b "${base}_cc" -O pg -t "${THREADS}"
fi

manifest="${PWD}/chunk_manifest.txt"
: > "${manifest}.new"
shopt -s nullglob
if compgen -G "*.gbz" > /dev/null; then
  for g in *.gbz; do
    readlink -f "$g"
  done | sort -V >> "${manifest}.new"
else
  for g in "${PREFIX}_cc"*.*; do
    [[ -f "$g" ]] || continue
    readlink -f "$g"
  done | sort -V >> "${manifest}.new"
fi
shopt -u nullglob

if [[ ! -s "${manifest}.new" ]]; then
  echo "ERROR: no chunk files (*.gbz or ${PREFIX}_cc*) found in ${OUT_DIR}" >&2
  exit 1
fi
mv "${manifest}.new" "${manifest}"

echo "Wrote $(wc -l < "${manifest}" | tr -d ' ') paths -> ${manifest}"
