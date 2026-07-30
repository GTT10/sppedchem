#!/usr/bin/env bash
#
# SpeedCHEM negative-test harness (Intel ifx only).
#
# Complements scripts/run_tests.sh (which checks a mechanism DOES integrate to
# ignition). This one checks that mechanisms using reaction features SpeedCHEM
# cannot evaluate are REFUSED (fail-closed) rather than silently mis-integrated.
#
# A negative case passes when the driver:
#   1. exits NON-ZERO (the fail-closed path uses `error stop 1`), and
#   2. prints the expected "does not evaluate" rejection message naming the
#      feature.
# It FAILS the harness if the driver exits 0 (feature silently accepted) or the
# expected message is missing.
#
# Usage:
#   scripts/run_neg_tests.sh
#
# Exit code: 0 = all negative cases refused as expected, non-zero otherwise.

set -euo pipefail

# ---- Locations --------------------------------------------------------------
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root_dir=$(cd -- "$script_dir/.." && pwd)
cd "$root_dir"

# ---- Compiler ---------------------------------------------------------------
TAG=ifx
FC=mpiifx
DRV_FLAGS=(-extend-source 132 -module "$TAG/mod" -O2)

if ! command -v "$FC" >/dev/null 2>&1; then
  echo "run_neg_tests.sh: required compiler wrapper '$FC' not found in PATH" >&2
  exit 2
fi

echo "======================================================================"
echo " SpeedCHEM negative-test harness (ifx)"
echo "   compiler : $FC"
echo "======================================================================"

# ---- 1. Build library + driver (same as run_tests.sh) -----------------------
echo "[1/2] Building library and smoke-test driver..."
make FC="$FC" -j1 >/dev/null
LIB="$TAG/libSpeedCHEM64.a"
[[ -f "$LIB" ]] || { echo "run_neg_tests.sh: library not produced ($LIB)" >&2; exit 1; }
DRV_BIN="$TAG/driver_smoke"
"$FC" "${DRV_FLAGS[@]}" -c test/driver_smoke.f90 -o "$TAG/build/driver_smoke.o"
"$FC" "$TAG/build/driver_smoke.o" "$LIB" -o "$DRV_BIN"

# ---- 2. Negative cases ------------------------------------------------------
# Each case: "<mechanism dir>|<substring the rejection message must contain>"
neg_cases=(
  "test/data_neg_ford/|modified species order (FORD/RORD)"
  # strict_chemkin: duplicate pressure within one PLOG reaction is
  # rejected by CKINTP (error stop) before a usable cklink is written.
  "test/data_neg_plog_dup/|not strictly ascending"
  "test/data_neg_plog_a/|must be finite and positive"
  "test/data_neg_plog_rev/|also uses REV/third-body/falloff syntax"
  "test/data_neg_plog_thirdbody/|also uses REV/third-body/falloff syntax"
)

echo "[2/2] Running negative cases..."
echo "----------------------------------------------------------------------"

fails=0
tmp_root=$(mktemp -d)
trap 'rm -rf -- "$tmp_root"' EXIT
for case in "${neg_cases[@]}"; do
  mech="${case%%|*}"
  want="${case#*|}"
  source_dir="$root_dir/$mech"
  mechdir="$source_dir"
  if [[ ! -f "${source_dir}therm.dat" ]]; then
    staged="$tmp_root/$(basename "${source_dir%/}")"
    mkdir -p "$staged"
    cp "${source_dir}chem.inp" "$staged/chem.inp"
    cp "$root_dir/test/data_plog/therm.dat" "$staged/therm.dat"
    mechdir="$staged/"
  fi

  # Clean any stale generated files so CKINTP/SCcklink run fresh.
  rm -f "${mechdir}"cklink "${mechdir}"chem.bin "${mechdir}"SpeedCHEM.* \
        "${mechdir}"dat.* "${mechdir}"chem.out 2>/dev/null || true

  echo "  case: $mech"
  set +e
  out=$(SC_MECHDIR="$mechdir" "$DRV_BIN" 2>&1)
  rc=$?
  set -e

  ok=1
  if [[ $rc -eq 0 ]]; then
    echo "    FAIL: driver exited 0 (feature was NOT refused)"
    ok=0
  fi
  if ! grep -qF "$want" <<<"$out"; then
    echo "    FAIL: expected rejection message not found:"
    echo "          want substring: $want"
    ok=0
  fi

  if [[ $ok -eq 1 ]]; then
    echo "    PASS: refused (exit $rc) with expected message"
  else
    echo "    --- driver output ---"
    sed 's/^/    | /' <<<"$out"
    fails=$((fails+1))
  fi
done

echo "======================================================================"
if [[ $fails -eq 0 ]]; then
  echo " RESULT: PASS — all negative cases refused as expected"
  echo "======================================================================"
  exit 0
else
  echo " RESULT: FAIL — $fails negative case(s) not refused as expected"
  echo "======================================================================"
  exit 1
fi
