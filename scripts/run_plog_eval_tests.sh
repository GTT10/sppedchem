#!/usr/bin/env bash
#
# SpeedCHEM PLOG rate/Jacobian test harness (Intel ifx only).
#
# The suite covers:
#
#   1. Closed-form PLOG rates and derivatives on hand-built tables.
#   2. Extreme Arrhenius parameters whose direct A*T**b intermediate would
#      overflow although the complete rate remains finite.
#   3. Full-RHS equivalence between plain PRF and a PLOG-equivalent PRF
#      mechanism, using a normalized physical mass-fraction state.
#   4. Complete constant-volume analytic Jacobian vs central differences.
#   5. Actual integration with the default LSODESJAC solver.
#   6. Two-rank MPI broadcast of the packed PLOG representation.
#
# Usage:  scripts/run_plog_eval_tests.sh
# Exit:   0 = all checks pass, non-zero otherwise.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root_dir=$(cd -- "$script_dir/.." && pwd)
cd "$root_dir"

TAG=ifx
FC=mpiifx
DRV_FLAGS=(-extend-source 132 -module "$TAG/mod" -O2)

if ! command -v "$FC" >/dev/null 2>&1; then
  echo "run_plog_eval_tests.sh: required compiler wrapper '$FC' not found in PATH" >&2
  exit 2
fi

echo "======================================================================"
echo " SpeedCHEM PLOG rate-evaluation test harness (ifx)"
echo "======================================================================"

echo "[1/7] Building library and PLOG test drivers..."
make FC="$FC" -j1 >/dev/null
LIB="$TAG/libSpeedCHEM64.a"
[[ -f "$LIB" ]] || { echo "run_plog_eval_tests.sh: library not produced ($LIB)" >&2; exit 1; }

build_drv() {  # $1 = source basename (without .f90)
  "$FC" "${DRV_FLAGS[@]}" -c "test/$1.f90" -o "$TAG/build/$1.o"
  "$FC" "$TAG/build/$1.o" -Wl,--start-group "$LIB" -Wl,--end-group -o "$TAG/$1"
}
build_drv driver_plog_eval
build_drv driver_plog_extreme
build_drv driver_rhs_probe
build_drv driver_plog_jacobian
build_drv driver_plog_mpi
build_drv driver_plog_integrate

fails=0

# ---- 1. Unit test -----------------------------------------------------------
echo "[2/7] PLOG rate and derivative unit test..."
echo "----------------------------------------------------------------------"
set +e
"$TAG/driver_plog_eval"
rc=$?
set -e
if [[ $rc -ne 0 ]]; then
  echo "    unit test FAILED (exit $rc)"
  fails=$((fails+1))
fi

# ---- 2. Extreme-range numerical stability ---------------------------------
echo "[3/7] Extreme grouped-PLOG overflow/underflow guard..."
echo "----------------------------------------------------------------------"
set +e
"$TAG/driver_plog_extreme"
rc=$?
set -e
if [[ $rc -ne 0 ]]; then
  echo "    extreme-range PLOG test FAILED (exit $rc)"
  fails=$((fails+1))
fi

# ---- 3. RHS equivalence -----------------------------------------------------
echo "[4/7] RHS equivalence: plain PRF vs PLOG-PRF (identical nodes)..."
echo "----------------------------------------------------------------------"
prf_dir="$root_dir/test/data/"
plog_dir="$root_dir/test/data_plog_prf/"
tmp_prf=$(mktemp)
tmp_plog=$(mktemp)
trap 'rm -f "$tmp_prf" "$tmp_plog"' EXIT

for d in "$prf_dir" "$plog_dir"; do
  rm -f "${d}"cklink "${d}"chem.bin "${d}"SpeedCHEM.* "${d}"dat.* "${d}"chem.out 2>/dev/null || true
done

SC_MECHDIR="$prf_dir"  "$TAG/driver_rhs_probe" 2>/dev/null | grep -E '^#|^[0-9]' > "$tmp_prf"
SC_MECHDIR="$plog_dir" "$TAG/driver_rhs_probe" 2>/dev/null | grep -E '^#|^[0-9]' > "$tmp_plog"

if [[ ! -s "$tmp_prf" || ! -s "$tmp_plog" ]]; then
  echo "    FAILED: one of the RHS probes produced no output"
  fails=$((fails+1))
elif diff "$tmp_prf" "$tmp_plog" >/dev/null; then
  echo "    RHS dumps are byte-identical (PLOG reproduces the Arrhenius rate):"
  head -4 "$tmp_prf" | sed 's/^/    | /'
  echo "    | ... (all $(wc -l < "$tmp_prf") lines match)"
else
  echo "    FAILED: RHS dumps differ"
  diff "$tmp_prf" "$tmp_plog" | head -10 | sed 's/^/    /'
  fails=$((fails+1))
fi

# ---- 4. Complete analytic Jacobian -----------------------------------------
echo "[5/7] PLOG analytic Jacobian vs central differences..."
echo "----------------------------------------------------------------------"
set +e
"$TAG/driver_plog_jacobian"
rc=$?
set -e
if [[ $rc -ne 0 ]]; then
  echo "    analytic Jacobian test FAILED (exit $rc)"
  fails=$((fails+1))
fi

# ---- 5. Integration ---------------------------------------------------------
echo "[6/7] PLOG integration with default LSODESJAC..."
echo "----------------------------------------------------------------------"
set +e
"$TAG/driver_plog_integrate"
rc=$?
set -e
if [[ $rc -ne 0 ]]; then
  echo "    PLOG integration test FAILED (exit $rc)"
  fails=$((fails+1))
fi

# ---- 6. MPI -----------------------------------------------------------------
echo "[7/7] PLOG MPI broadcast..."
echo "----------------------------------------------------------------------"
set +e
mpiexec -n 2 "$TAG/driver_plog_mpi"
rc=$?
set -e
if [[ $rc -ne 0 ]]; then
  echo "    MPI broadcast test FAILED (exit $rc)"
  fails=$((fails+1))
fi

echo "======================================================================"
if [[ $fails -eq 0 ]]; then
  echo " RESULT: PASS — PLOG rates, numerical stability, RHS, Jacobian, integration, and MPI"
  echo "======================================================================"
  exit 0
else
  echo " RESULT: FAIL — $fails PLOG evaluation check(s) failed"
  echo "======================================================================"
  exit 1
fi
