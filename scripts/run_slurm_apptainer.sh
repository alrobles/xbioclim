#!/usr/bin/env bash
#SBATCH --job-name=xbioclim
#SBATCH --output=logs/xbioclim_%A_%a.out
#SBATCH --error=logs/xbioclim_%A_%a.err
#SBATCH --array=0-359          # one task per degree of longitude (-180..179)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
# GPU nodes — uncomment the two lines below if using a CUDA-enabled image:
# #SBATCH --partition=gpu
# #SBATCH --gres=gpu:1

set -euo pipefail

# Some clusters require loading the apptainer module; uncomment if needed:
# module load apptainer

# ---------- Configuration ----------
# Root directory containing the input CHELSA monthly rasters.
# Expected layout: ${INPUT_DIR}/${VAR}_${MM}.tif  (e.g. tas_01.tif)
INPUT_DIR="${INPUT_DIR:?Set INPUT_DIR to the directory containing monthly rasters}"
OUTPUT_DIR="${OUTPUT_DIR:?Set OUTPUT_DIR to the desired output directory}"
SIF="${SIF:-$HOME/xbioclim.sif}"
TILE_SIZE="${TILE_SIZE:-512}"

# Guard: make sure the SIF image exists
if [[ ! -f "${SIF}" ]]; then
    echo "ERROR: Apptainer image not found at ${SIF}" >&2
    echo "Build it with:  sudo apptainer build xbioclim.sif apptainer/xbioclim.def" >&2
    exit 1
fi

# ---------- Derived values ----------
LON_IDX=${SLURM_ARRAY_TASK_ID}          # 0..359
LON=$(( LON_IDX - 180 ))                # -180..179

STRIP_DIR="${OUTPUT_DIR}/lon_${LON}"
mkdir -p "${STRIP_DIR}" "logs"

# Build file lists (12 months per variable)
make_files() {
    local var=$1
    local files=()
    for m in $(seq -w 1 12); do
        files+=("${INPUT_DIR}/${var}_${m}.tif")
    done
    echo "${files[@]}"
}

echo "[$(date)] Starting longitude strip ${LON} (task ${LON_IDX})"

apptainer exec \
    --bind "${INPUT_DIR}:${INPUT_DIR}" \
    --bind "${OUTPUT_DIR}:${OUTPUT_DIR}" \
    "${SIF}" \
    xbioclim_cli \
        --tas    $(make_files tas) \
        --tasmax $(make_files tasmax) \
        --tasmin $(make_files tasmin) \
        --pr     $(make_files pr) \
        --outdir "${STRIP_DIR}" \
        --prefix "bio" \
        --tile-size "${TILE_SIZE}"

echo "[$(date)] Finished longitude strip ${LON}"
