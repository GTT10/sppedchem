#!/usr/bin/env bash
#
# SpeedCHEM PLOG rate-evaluation test harness (Intel ifx only) -- stage 2.
#
# Two checks, both independent of the (separately-flaky) stiff ODE solvers:
#
#   1. Unit test (test/driver_plog_eval.f90): calls reacpar::plog_kinf_eval
#      directly on a hand-built 3-node PLOG reaction and compares the
#      interpolated rate constant against an independent closed-form
#      reference at nodes, geometric-mean pressures, out-of-range (clamp),
#      and several temperatures.
#
#   2. RHS equivalence (test/driver_rhs_probe.f90): evaluates the full
#      constant-volume RHS d(state)/dt once for
#         - test/data/          (plain PRF, no PLOG)
#         - test/data_plog_prf/ (PRF with reaction "O+H2=H+OH" replaced by a
#                                3-node PLOG whose nodes reproduce that same
#                                Arrhenius rate at every pressure)
#      and asserts the two derivative dumps are byte-identical. This proves
#      plog_kinf_eval is wired into the real RHS (mass_action) and, when the
#      PLOG table encodes the original Arrhenius rate, reproduces it exactly.
#
# Usage:  scripts/run_plog_eval_tests.sh
# Exit:   0 = both checks pass, non-zero otherwise.

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

echo "[1/3] Building library and PLOG test drivers..."
make FC="$FC" -j1 >/dev/null
LIB="$TAG/libSpeedCHEM64.a"
[[ -f "$LIB" ]] || { echo "run_plog_eval_tests.sh: library not produced ($LIB)" >&2; exit 1; }

build_drv() {  # $1 = source basename (without .f90)
  "$FC" "${DRV_FLAGS[@]}" -c "test/$1.f90" -o "$TAG/build/$1.o"
  "$FC" "$TAG/build/$1.o" "$LIB" -o "$TAG/$1"
}
build_drv driver_plog_eval
build_drv driver_rhs_probe

fails=0

# ---- 1. Unit test -----------------------------------------------------------
echo "[2/3] PLOG rate unit test (plog_kinf_eval vs closed form)..."
echo "----------------------------------------------------------------------"
set +e
"$TAG/driver_plog_eval"
rc=$?
set -e
if [[ $rc -ne 0 ]]; then
  echo "    unit test FAILED (exit $rc)"
  fails=$((fails+1))
fi

# ---- 2. RHS equivalence -----------------------------------------------------
echo "[3/3] RHS equivalence: plain PRF vs PLOG-PRF (identical nodes)..."
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

echo "======================================================================"
if [[ $fails -eq 0 ]]; then
  echo " RESULT: PASS — PLOG rates evaluate correctly and feed the RHS"
  echo "======================================================================"
  exit 0
else
  echo " RESULT: FAIL — $fails PLOG evaluation check(s) failed"
  echo "======================================================================"
  exit 1
fi
