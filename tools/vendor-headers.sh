#!/usr/bin/env bash
# tools/vendor-headers.sh
#
# Downloads xtl and xtensor header-only libraries and places them in
# inst/include/ for vendoring with the R package.
#
# Usage (from repository root):
#   bash tools/vendor-headers.sh
#
# Requirements: curl, tar

set -euo pipefail

XTL_VERSION="0.7.7"
XTENSOR_VERSION="0.25.0"

INST_INCLUDE="$(cd "$(dirname "$0")/.." && pwd)/inst/include"

echo "==> Vendoring headers into ${INST_INCLUDE}"

# --------------------------------------------------------------------------
# xtl
# --------------------------------------------------------------------------
XTL_URL="https://github.com/xtensor-stack/xtl/archive/refs/tags/${XTL_VERSION}.tar.gz"
XTL_DIR="${INST_INCLUDE}/xtl"

echo ""
echo "--> Downloading xtl ${XTL_VERSION}..."
mkdir -p "${XTL_DIR}"
curl -fsSL "${XTL_URL}" \
  | tar -xz --strip-components=2 -C "${XTL_DIR}" \
      "xtl-${XTL_VERSION}/include/xtl"

echo "    xtl headers installed to ${XTL_DIR}"

# --------------------------------------------------------------------------
# xtensor
# --------------------------------------------------------------------------
XTENSOR_URL="https://github.com/xtensor-stack/xtensor/archive/refs/tags/${XTENSOR_VERSION}.tar.gz"
XTENSOR_DIR="${INST_INCLUDE}/xtensor"

echo ""
echo "--> Downloading xtensor ${XTENSOR_VERSION}..."
mkdir -p "${XTENSOR_DIR}"
curl -fsSL "${XTENSOR_URL}" \
  | tar -xz --strip-components=2 -C "${XTENSOR_DIR}" \
      "xtensor-${XTENSOR_VERSION}/include/xtensor"

echo "    xtensor headers installed to ${XTENSOR_DIR}"

echo ""
echo "==> Done. Header counts:"
echo "    xtl:     $(find "${XTL_DIR}"     -name '*.hpp' | wc -l) files"
echo "    xtensor: $(find "${XTENSOR_DIR}" -name '*.hpp' | wc -l) files"
