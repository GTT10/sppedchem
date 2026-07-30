#!/usr/bin/env bash
#
# SpeedCHEM PLOG data-plumbing and unit-conversion test harness.
#
# Exercises the PLOG parse -> cklink v2 write -> read -> packed arrays
# path and the REACTIONS MOLECULES conversion. It builds
# the parse-only driver (test/driver_plog.f90), runs it on the bundled
# PLOG mechanism (test/data_plog/), and asserts the
# canonical PLOG dump matches the golden file test/data_plog/plog_expected.txt.
#
# Unsupported PLOG combinations and malformed input are covered separately
# by run_neg_tests.sh.
#
# Usage:   scripts/run_plog_tests.sh
# Exit:    0 = dump matches golden, non-zero otherwise.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root_dir=$(cd -- "$script_dir/.." && pwd)
cd "$root_dir"

TAG=ifx
FC=mpiifx
DRV_FLAGS=(-extend-source 132 -module "$TAG/mod" -O2)

if ! command -v "$FC" >/dev/null 2>&1; then
  echo "run_plog_tests.sh: required compiler wrapper '$FC' not found in PATH" >&2
  exit 2
fi

echo "======================================================================"
echo " SpeedCHEM PLOG plumbing test harness (ifx)"
echo "   compiler : $FC"
echo "======================================================================"

# ---- 1. Build library + parse-only PLOG driver ------------------------------
echo "[1/4] Building library and PLOG parse-only driver..."
make FC="$FC" -j1 >/dev/null
LIB="$TAG/libSpeedCHEM64.a"
[[ -f "$LIB" ]] || { echo "run_plog_tests.sh: library not produced ($LIB)" >&2; exit 1; }
DRV_BIN="$TAG/driver_plog"
"$FC" "${DRV_FLAGS[@]}" -c test/driver_plog.f90 -o "$TAG/build/driver_plog.o"
"$FC" "$TAG/build/driver_plog.o" "$LIB" -o "$DRV_BIN"

# ---- 2. Parse + link + dump -------------------------------------------------
mech="test/data_plog/"
mechdir="$root_dir/$mech"
golden="${mechdir}plog_expected.txt"
[[ -f "$golden" ]] || { echo "run_plog_tests.sh: missing golden file $golden" >&2; exit 1; }

echo "[2/4] Parsing + linking + dumping PLOG data ($mech)..."
# Force a fresh CKINTP + SCcklink pass.
rm -f "${mechdir}"cklink "${mechdir}"chem.bin "${mechdir}"SpeedCHEM.* \
      "${mechdir}"dat.* "${mechdir}"chem.out 2>/dev/null || true

got=$(SC_MECHDIR="$mechdir" "$DRV_BIN" 2>/dev/null \
      | grep -E '^# PLOG|^n_plog|^[0-9]')

# ---- 3. Compare against golden ---------------------------------------------
echo "[3/4] Comparing canonical dump against golden..."
echo "----------------------------------------------------------------------"
if diff <(printf '%s\n' "$got") "$golden" >/dev/null; then
  echo " canonical PLOG dump matches golden:"
  printf '%s\n' "$got" | sed 's/^/    | /'
else
  echo " MISMATCH between dump and golden:"
  diff <(printf '%s\n' "$got") "$golden" | sed 's/^/    /' || true
  echo "======================================================================"
  echo " RESULT: FAIL — PLOG dump does not match golden"
  echo "======================================================================"
  exit 1
fi

echo "[4/4] Verifying REACTIONS MOLECULES conversion..."
tmp_units=$(mktemp -d)
trap 'rm -rf -- "$tmp_units"' EXIT
cp test/data_plog_molecules/chem.inp "$tmp_units/chem.inp"
cp test/data_plog/therm.dat "$tmp_units/therm.dat"
units_dump=$(SC_MECHDIR="$tmp_units/" "$DRV_BIN" 2>/dev/null |
             awk '/^[0-9]/{print $4; exit}')
awk -v got="$units_dump" 'BEGIN {
  ref=1.0e12*6.0221367e23;
  err=(got-ref)/ref; if (err<0) err=-err;
  if (err>1.0e-12) exit 1
}' || {
  echo " PLOG MOLECULES conversion mismatch: got $units_dump" >&2
  exit 1
}
echo " converted PLOG A = $units_dump"
echo "======================================================================"
echo " RESULT: PASS — PLOG cklink round-trip and units are correct"
echo "======================================================================"
