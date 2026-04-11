#!/usr/bin/env bash
# validate_cog_outputs.sh — validate xbioclim output TIFFs are valid COGs
#
# Usage: validate_cog_outputs.sh <output_dir>
#
# For every *.tif file in <output_dir> this script:
#   1. Runs `gdalinfo -checksum` and asserts it exits successfully.
#   2. Runs scripts/check_cog.py to assert the file is a valid COG.
#
# Exit code: 0 if all files pass, 1 if any file fails.

set -euo pipefail

OUTPUT_DIR="${1:-.}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COG_CHECKER="${SCRIPT_DIR}/check_cog.py"

if [ ! -f "${COG_CHECKER}" ]; then
    echo "ERROR: COG checker not found at ${COG_CHECKER}" >&2
    exit 1
fi

shopt -s nullglob
tif_files=("${OUTPUT_DIR}"/*.tif)
shopt -u nullglob

if [ "${#tif_files[@]}" -eq 0 ]; then
    echo "ERROR: No .tif files found in '${OUTPUT_DIR}'" >&2
    exit 1
fi

FAILED=0

for f in "${tif_files[@]}"; do
    name=$(basename "$f")
    echo "=== Validating: ${name} ==="

    # 1. gdalinfo -checksum must succeed and produce non-empty output.
    if ! gdalinfo -checksum "$f"; then
        echo "FAIL [${name}]: gdalinfo -checksum returned a non-zero exit code." >&2
        FAILED=1
        continue
    fi

    # 2. Validate COG structure.
    if ! python3 "${COG_CHECKER}" "$f"; then
        echo "FAIL [${name}]: COG structure validation failed." >&2
        FAILED=1
    fi
done

if [ "${FAILED}" -ne 0 ]; then
    echo ""
    echo "ERROR: One or more output files failed COG + checksum validation." >&2
    exit 1
fi

echo ""
echo "SUCCESS: All ${#tif_files[@]} output file(s) passed COG + checksum validation."
