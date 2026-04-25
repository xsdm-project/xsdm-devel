#!/usr/bin/env bash
set -euo pipefail

# --- Config (edit if you want different refs) ---
XTENSOR_REPO="${XTENSOR_REPO:-https://github.com/alrobles/xtensor.git}"
XTL_REPO="${XTL_REPO:-https://github.com/xtensor-stack/xtl.git}"
XSIMD_REPO="${XSIMD_REPO:-https://github.com/xtensor-stack/xsimd.git}"

# Pin versions/commits for reproducibility
XTENSOR_REF="${XTENSOR_REF:-master}"
XTL_REF="${XTL_REF:-master}"
XSIMD_REF="${XSIMD_REF:-master}"

# Where to vendor
VENDOR_DIR="src/vendor"

# Temp workspace
TMP_DIR="$(mktemp -d)"
cleanup() { rm -rf "$TMP_DIR"; }
trap cleanup EXIT

echo "Vendoring into: ${VENDOR_DIR}"
mkdir -p "${VENDOR_DIR}"

clone_and_checkout () {
  local repo="$1"
  local ref="$2"
  local dest="$3"
  git clone --depth 1 --branch "$ref" "$repo" "$dest" 2>/dev/null || {
    # fallback if ref isn't a branch/tag (e.g. commit SHA)
    git clone "$repo" "$dest"
    (cd "$dest" && git checkout "$ref")
  }
}

echo "Cloning xtl..."
clone_and_checkout "$XTL_REPO" "$XTL_REF" "${TMP_DIR}/xtl"

echo "Cloning xsimd..."
clone_and_checkout "$XSIMD_REPO" "$XSIMD_REF" "${TMP_DIR}/xsimd"

echo "Cloning xtensor..."
clone_and_checkout "$XTENSOR_REPO" "$XTENSOR_REF" "${TMP_DIR}/xtensor"

# Clean existing vendor dirs (optional; comment out if you prefer manual)
rm -rf "${VENDOR_DIR}/xtl" "${VENDOR_DIR}/xsimd" "${VENDOR_DIR}/xtensor"

mkdir -p "${VENDOR_DIR}/xtl" "${VENDOR_DIR}/xsimd" "${VENDOR_DIR}/xtensor"

echo "Copying headers..."
# xtl headers live in xtl/include/xtl
rsync -a --delete "${TMP_DIR}/xtl/include/xtl" "${VENDOR_DIR}/xtl/"
# xsimd headers live in xsimd/include/xsimd
rsync -a --delete "${TMP_DIR}/xsimd/include/xsimd" "${VENDOR_DIR}/xsimd/"
# xtensor headers live in xtensor/include/xtensor
rsync -a --delete "${TMP_DIR}/xtensor/include/xtensor" "${VENDOR_DIR}/xtensor/"

echo "Writing vendor metadata..."
cat > "${VENDOR_DIR}/VENDORED_XTENSOR.md" <<EOF
# Vendored headers

This directory vendors header-only dependencies used by the C++ code in this package.

- xtensor: ${XTENSOR_REPO} @ ${XTENSOR_REF}
- xtl: ${XTL_REPO} @ ${XTL_REF}
- xsimd: ${XSIMD_REPO} @ ${XSIMD_REF}

Vendored by: scripts/vendor-xtensor.sh
EOF

echo "Done. Next steps:"
echo "  1) Update src/Makevars to add -Isrc/vendor/xtensor -Isrc/vendor/xtl -Isrc/vendor/xsimd"
echo "  2) Run R CMD check"
echo "  3) Commit + push"