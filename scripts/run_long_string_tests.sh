#!/usr/bin/env bash
#
# Regression test for the CHEMKIN/NASA-7 string boundaries:
# 18-character species identifiers and mechanism records longer than 80 chars.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root_dir=$(cd -- "$script_dir/.." && pwd)
cd "$root_dir"

TAG=ifx
FC=mpiifx
DRV_FLAGS=(-extend-source 132 -module "$TAG/mod" -O2)

if ! command -v "$FC" >/dev/null 2>&1; then
  echo "run_long_string_tests.sh: required compiler wrapper '$FC' not found" >&2
  exit 2
fi
if ! command -v mpiexec >/dev/null 2>&1; then
  echo "run_long_string_tests.sh: required launcher 'mpiexec' not found" >&2
  exit 2
fi

fixture="$root_dir/test/data_long_strings"
species_width=$(awk '/^LONGSPECIES/{print length($1); exit}' "$fixture/chem.inp")
record_width=$(awk '{if (length > max) max=length} END{print max}' "$fixture/chem.inp")
if [[ "$species_width" -ne 18 || "$record_width" -le 80 ]]; then
  echo "run_long_string_tests.sh: fixture does not exercise both boundaries" >&2
  exit 1
fi

tmp_mech=$(mktemp -d)
trap 'rm -rf -- "$tmp_mech"' EXIT
cp "$fixture/chem.inp" "$fixture/therm.dat" "$tmp_mech/"

echo "======================================================================"
echo " SpeedCHEM long chemistry string regression"
echo "   species width : $species_width"
echo "   record width  : $record_width"
echo "======================================================================"

echo "[1/4] Building library..."
make FC="$FC" -j1 >/dev/null
LIB="$TAG/libSpeedCHEM64.a"

echo "[2/4] Building regression driver..."
"$FC" "${DRV_FLAGS[@]}" -c test/driver_long_strings.f90 \
  -o "$TAG/build/driver_long_strings.o"
"$FC" "$TAG/build/driver_long_strings.o" \
  -Wl,--start-group "$LIB" -Wl,--end-group -o "$TAG/driver_long_strings"

echo "[3/4] Parsing and round-tripping on one MPI rank..."
single_out=$(SC_MECHDIR="$tmp_mech/" "$TAG/driver_long_strings")
printf '%s\n' "$single_out"
grep -q '^RESULT: PASS ' <<<"$single_out"

rm -f "$tmp_mech/cklink" "$tmp_mech/chem.out" "$tmp_mech/SpeedCHEM.out" \
      "$tmp_mech"/dat.*

echo "[4/4] Broadcasting and validating on two MPI ranks..."
mpi_out=$(SC_MECHDIR="$tmp_mech/" mpiexec -n 2 "$TAG/driver_long_strings")
printf '%s\n' "$mpi_out"
grep -q '^RESULT: PASS ' <<<"$mpi_out"

echo "======================================================================"
echo " RESULT: PASS - long chemistry strings are preserved"
echo "======================================================================"
