#!/bin/bash
#SBATCH --job-name=tags_d46
#SBATCH --partition=long
#SBATCH --mail-user=seeskand@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=300gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:0
#SBATCH --output=%x.%j.log
#SBATCH --time=20:00:00
#SBATCH --array=1-26
#SBATCH --no-requeue
#SBATCH --exclude=phoenix-00
#
# Build per-chunk tags from EXISTING chunk graph files only.
# This script does NOT run any chunking.
#
# If build_tags crashes loading GBZ with "Expected version 1, got version 2", rebuild tools
# with a gbwtgraph that supports GBZ v2: from repo root run `make gbwtgraph-lib` then `make`
# (requires SDSL_DIR / gbwt / handlegraph same as your pangenome-index build).
#
# Per array task:
#   1) Pick one chunk file from manifest (chunks_cc_*.vg by default).
#   2) Move chunk into its own WORKDIR/<chunk_name>/ folder.
#   3) Convert chunk graph to GBZ.
#   4) Build rlbwt via gbz_extract.
#   5) Build tags via build_tags.
#
set -euo pipefail
source ~/.bashrc

DATA_DIR="${DATA_DIR:-/private/groups/cgl/seeskand/slurm-jobs/data}"
WORKDIR="${WORKDIR:-/private/groups/cgl/seeskand/graph_index/12.hprc-v2.0-mc-chm13.d46/chromosome_tags}"
CHUNK_STAGING="${CHUNK_STAGING:-${WORKDIR}/.staging}"
EXISTING_CHUNK_DIR="${EXISTING_CHUNK_DIR:-${WORKDIR}}"
EXISTING_CHUNK_GLOB="${EXISTING_CHUNK_GLOB:-chunks_cc_*.vg}"

VG_BIN="${VG_BIN:-vg}"
GBZ_EXTRACT="${GBZ_EXTRACT:-/private/groups/cgl/seeskand/1.server/Giraffe_server/deps/gbwtgraph/bin/gbz_extract}"
GRLBWT_CLI="${GRLBWT_CLI:-/private/groups/cgl/seeskand/software/grlBWT/build/grlbwt-cli}"
BUILD_TAGS="${BUILD_TAGS:-/private/groups/cgl/seeskand/pangenome-index/bin/build_tags}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"

export PATH=/private/groups/cgl/seeskand/vg/vg/bin:/private/groups/cgl/seeskand/vg/vg/scripts:${PATH}

CHUNK_MANIFEST="${CHUNK_STAGING}/chunk_manifest.txt"
PREP_DONE="${CHUNK_STAGING}/.prep.done"
PREP_LOCK="${CHUNK_STAGING}/.prep.lock"

prepare_manifest() {
  mkdir -p "${CHUNK_STAGING}"
  : > "${CHUNK_MANIFEST}.tmp"
  # Build manifest from chunk IDs so the script can resume even after chunks
  # have already been moved into per-chunk directories.
  local f d id
  local -a ids=()
  declare -A seen=()
  shopt -s nullglob
  for f in "${EXISTING_CHUNK_DIR}"/${EXISTING_CHUNK_GLOB}; do
    id="$(basename "${f}")"
    id="${id%.vg}"
    id="${id%.pg}"
    id="${id%.gbz}"
    if [[ -n "${id}" && -z "${seen[$id]+x}" ]]; then
      seen["$id"]=1
      ids+=("$id")
    fi
  done
  for d in "${WORKDIR}"/chunks_cc_*; do
    [[ -d "${d}" ]] || continue
    id="$(basename "${d}")"
    if [[ -n "${id}" && -z "${seen[$id]+x}" ]]; then
      seen["$id"]=1
      ids+=("$id")
    fi
  done
  shopt -u nullglob
  if [[ "${#ids[@]}" -eq 0 ]]; then
    echo "[prep] ERROR: no chunk IDs found in ${EXISTING_CHUNK_DIR}/${EXISTING_CHUNK_GLOB} or ${WORKDIR}/chunks_cc_*" >&2
    exit 1
  fi
  printf '%s\n' "${ids[@]}" | sort -V > "${CHUNK_MANIFEST}.tmp"
  mv "${CHUNK_MANIFEST}.tmp" "${CHUNK_MANIFEST}"
  local cnt
  cnt="$(wc -l < "${CHUNK_MANIFEST}" | tr -d ' ')"
  echo "[prep] Manifest has ${cnt} chunk ID(s): ${CHUNK_MANIFEST}"
}

prepare_manifest_locked() {
  mkdir -p "${CHUNK_STAGING}"
  exec {lock_fd}>"${PREP_LOCK}"
  flock "${lock_fd}"
  if [[ ! -f "${PREP_DONE}" ]]; then
    echo "[prep] Building manifest under lock (job=${SLURM_JOB_ID:-local} task=${SLURM_ARRAY_TASK_ID:-1})"
    prepare_manifest
    touch "${PREP_DONE}"
  else
    echo "[prep] Reusing existing manifest (${PREP_DONE} present)."
  fi
  flock -u "${lock_fd}"
  exec {lock_fd}>&-
}

wait_manifest() {
  local waited=0
  while [[ ! -f "${PREP_DONE}" ]]; do
    sleep 10
    waited=$((waited + 10))
    if [[ "${waited}" -ge 864000 ]]; then
      echo "ERROR: timed out waiting for ${PREP_DONE}" >&2
      exit 1
    fi
    echo "[wait] ${waited}s waiting for manifest..."
  done
}

prepare_manifest_locked
wait_manifest

CHUNK_ID="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CHUNK_MANIFEST}")"
if [[ -z "${CHUNK_ID}" ]]; then
  echo "[task] No chunk on manifest line ${SLURM_ARRAY_TASK_ID} (manifest has $(wc -l < "${CHUNK_MANIFEST}" | tr -d ' ') lines); exiting 0."
  exit 0
fi

WORK_ID="${CHUNK_ID}"
WORK_ID="$(echo "${WORK_ID}" | sed 's/[#\/]/_/g')"

mkdir -p "${WORKDIR}/${WORK_ID}"
cd "${WORKDIR}/${WORK_ID}"
echo "[task] array=${SLURM_ARRAY_TASK_ID} chunk_id=${CHUNK_ID} work_id=${WORK_ID} pwd=${PWD}"

CHUNK_GBZ="${PWD}/${WORK_ID}.gbz"
if [[ ! -s "${CHUNK_GBZ}" ]]; then
  LOCAL_CHUNK_GRAPH="${PWD}/${WORK_ID}.vg"
  if [[ ! -f "${LOCAL_CHUNK_GRAPH}" ]]; then
    SOURCE_CHUNK_GRAPH="${EXISTING_CHUNK_DIR}/${WORK_ID}.vg"
    if [[ -f "${SOURCE_CHUNK_GRAPH}" ]]; then
      # Resume-safe: only move if still at source location.
      mv "${SOURCE_CHUNK_GRAPH}" "${LOCAL_CHUNK_GRAPH}"
    fi
  fi
  if [[ ! -f "${LOCAL_CHUNK_GRAPH}" ]]; then
    echo "ERROR: no chunk graph found for ${WORK_ID} (expected ${PWD}/${WORK_ID}.vg or ${EXISTING_CHUNK_DIR}/${WORK_ID}.vg)" >&2
    exit 1
  fi
  /usr/bin/time -v "${VG_BIN}" gbwt -p -x "${LOCAL_CHUNK_GRAPH}" -g "${CHUNK_GBZ}" --gbz-format -E
else
  echo "[task] Reusing existing GBZ: ${CHUNK_GBZ}"
fi

if [[ -s "${WORK_ID}.tags" ]]; then
  echo "[task] Tags already exist, skipping: ${PWD}/${WORK_ID}.tags"
  exit 0
fi

if [[ ! -s "${WORK_ID}_info.rl_bwt" ]]; then
  /usr/bin/time -v "${GBZ_EXTRACT}" -t "${THREADS}" -b -p "${CHUNK_GBZ}" > "${WORK_ID}_info.tmp"
  mv "${WORK_ID}_info.tmp" "${WORK_ID}_info"
  TMP_GRLBWT_DIR="${PWD}/tmp_grlbwt_${WORK_ID}"
  mkdir -p "${TMP_GRLBWT_DIR}"
  /usr/bin/time -v "${GRLBWT_CLI}" -t "${THREADS}" -T "${TMP_GRLBWT_DIR}" "${WORK_ID}_info"
else
  echo "[task] Reusing existing rlbwt: ${PWD}/${WORK_ID}_info.rl_bwt"
fi

/usr/bin/time -v "${BUILD_TAGS}" "${CHUNK_GBZ}" "${WORK_ID}_info.rl_bwt" "${WORK_ID}.tags"

echo "[task] done: ${PWD}/${WORK_ID}.tags"
