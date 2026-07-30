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
echo "[1/7] Building library and PLOG link/read drivers..."
make FC="$FC" -j1 >/dev/null
LIB="$TAG/libSpeedCHEM64.a"
[[ -f "$LIB" ]] || { echo "run_plog_tests.sh: library not produced ($LIB)" >&2; exit 1; }
DRV_BIN="$TAG/driver_plog"
"$FC" "${DRV_FLAGS[@]}" -c test/driver_plog.f90 -o "$TAG/build/driver_plog.o"
"$FC" "$TAG/build/driver_plog.o" "$LIB" -o "$DRV_BIN"
READ_BIN="$TAG/driver_cklink_read"
"$FC" "${DRV_FLAGS[@]}" -c test/driver_cklink_read.f90 \
  -o "$TAG/build/driver_cklink_read.o"
"$FC" "$TAG/build/driver_cklink_read.o" "$LIB" -o "$READ_BIN"

# ---- 2. Parse + link + dump -------------------------------------------------
mech="test/data_plog/"
mechdir="$root_dir/$mech"
golden="${mechdir}plog_expected.txt"
[[ -f "$golden" ]] || { echo "run_plog_tests.sh: missing golden file $golden" >&2; exit 1; }

echo "[2/7] Parsing + linking + dumping PLOG data ($mech)..."
# Force a fresh CKINTP + SCcklink pass.
rm -f "${mechdir}"cklink "${mechdir}"chem.bin "${mechdir}"SpeedCHEM.* \
      "${mechdir}"dat.* "${mechdir}"chem.out 2>/dev/null || true

got=$(SC_MECHDIR="$mechdir" "$DRV_BIN" 2>/dev/null \
      | grep -E '^# PLOG|^n_plog|^[0-9]')

# ---- 3. Compare against golden ---------------------------------------------
echo "[3/7] Comparing canonical dump against golden..."
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

echo "[4/7] Verifying grouped same-pressure PLOG round-trip..."
tmp_grouped=$(mktemp -d)
cp test/data_plog_grouped/chem.inp "$tmp_grouped/chem.inp"
cp test/data_plog/therm.dat "$tmp_grouped/therm.dat"
grouped_dump=$(SC_MECHDIR="$tmp_grouped/" "$DRV_BIN" 2>/dev/null |
               grep -E '^n_plog|^[0-9]')
grouped_entries=$(grep -c '^[0-9]' <<<"$grouped_dump")
grouped_unique_pressures=$(awk '/^[0-9]/{print $3}' <<<"$grouped_dump" |
                           sort -u | wc -l)
if [[ "$grouped_entries" -ne 4 || "$grouped_unique_pressures" -ne 3 ]]; then
  echo " grouped PLOG entries were not preserved:" >&2
  printf '%s\n' "$grouped_dump" >&2
  exit 1
fi
echo " grouped PLOG: 4 Arrhenius entries at 3 pressure nodes"

echo "[5/7] Verifying REACTIONS MOLECULES conversion..."
tmp_units=$(mktemp -d)
trap 'rm -rf -- "$tmp_grouped" "$tmp_units"' EXIT
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

echo "[6/7] Rejecting an old/foreign cklink without v2 magic..."
tmp_old=$(mktemp -d)
printf 'legacy cklink\n' >"$tmp_old/cklink"
set +e
old_out=$(SC_MECHDIR="$tmp_old/" "$READ_BIN" 2>&1)
old_rc=$?
set -e
if [[ $old_rc -eq 0 ]] || ! grep -q 'missing magic record' <<<"$old_out"; then
  echo " old/foreign cklink was not rejected cleanly" >&2
  printf '%s\n' "$old_out" >&2
  exit 1
fi
echo " old/foreign cklink rejected (exit $old_rc)"

echo "[7/7] Rejecting a truncated PLOG section..."
tmp_truncated=$(mktemp -d)
cp "${mechdir}cklink" "$tmp_truncated/cklink"
size=$(stat -c %s "$tmp_truncated/cklink")
truncate -s "$((size-4))" "$tmp_truncated/cklink"
set +e
truncated_out=$(SC_MECHDIR="$tmp_truncated/" "$READ_BIN" 2>&1)
truncated_rc=$?
set -e
if [[ $truncated_rc -eq 0 ]] || \
   ! grep -q 'PLOG checksum record missing' <<<"$truncated_out"; then
  echo " truncated PLOG section was not rejected cleanly" >&2
  printf '%s\n' "$truncated_out" >&2
  exit 1
fi
echo " truncated PLOG section rejected (exit $truncated_rc)"

rm -rf -- "$tmp_old" "$tmp_truncated"
echo "======================================================================"
echo " RESULT: PASS — PLOG cklink round-trip, units, and rejection are correct"
echo "======================================================================"
