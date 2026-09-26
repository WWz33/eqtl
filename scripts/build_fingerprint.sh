#!/usr/bin/env bash
# Print what a build actually linked, to be recorded before an A/B comparison.
#
# Same source linked against a different GSL/BLAS resolves to a numerically
# different binary — perm>0 cases diverge because BLAS takes part in the
# eigendecomposition and the dgemms — so two builds must be fingerprinted
# before their outputs are compared, or the comparison measures the build
# environment instead of the change.
#
# Usage: scripts/build_fingerprint.sh [binary]   (default: <repo>/eqtl)
# Expected on this host: libgsl.so.28 and libmkl_rt.so.3 out of
#   /sri/home/sri2025201067/miniforge3/envs/eqtl/lib
# which is what `make GSL_PREFIX=/sri/home/sri2025201067/miniforge3/envs/eqtl -j`
# produces (the same prefix an activated env supplies via CONDA_PREFIX).
# Anything else (base env's libgsl.so.25 + libopenblas.so.0, say) has to be
# recorded as such next to the numbers it produced.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="${1:-$ROOT/eqtl}"
[[ -x "$BIN" ]] || { echo "[E] no binary: $BIN" >&2; exit 1; }

echo "binary: $BIN"
echo "md5:    $(md5sum "$BIN" | cut -d' ' -f1)"
echo "built:  $(stat -c '%y' "$BIN" | cut -d. -f1)  size $(stat -c %s "$BIN")"
echo "libs:"
ldd "$BIN" | grep -Ei "gsl|blas|mkl|lapack" | sed 's/^/  /' || echo "  (none)"
