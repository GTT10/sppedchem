#!/usr/bin/env bash
#
# Verify that a process can load mechanism A, finalize it, load a PLOG
# mechanism B, finalize it, and reload A without stale global state.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root_dir=$(cd -- "$script_dir/.." && pwd)
cd "$root_dir"

TAG=ifx
FC=mpiifx
DRV_FLAGS=(-extend-source 132 -module "$TAG/mod" -O2)

if ! command -v "$FC" >/dev/null 2>&1; then
  echo "run_reload_tests.sh: required compiler wrapper '$FC' not found" >&2
  exit 2
fi

stage_a=$(mktemp -d)
stage_b=$(mktemp -d)
trap 'rm -rf -- "$stage_a" "$stage_b"' EXIT

echo "======================================================================"
echo " SpeedCHEM mechanism reload test (A -> PLOG B -> A)"
echo "======================================================================"

echo "[1/5] Building library and drivers..."
make FC="$FC" -j1 >/dev/null
LIB="$TAG/libSpeedCHEM64.a"

"$FC" "${DRV_FLAGS[@]}" -c test/driver_plog.f90 \
  -o "$TAG/build/driver_plog_reload_link.o"
"$FC" "$TAG/build/driver_plog_reload_link.o" "$LIB" \
  -o "$TAG/driver_plog_reload_link"

"$FC" "${DRV_FLAGS[@]}" -c test/driver_mechanism_reload.f90 \
  -o "$TAG/build/driver_mechanism_reload.o"
"$FC" "$TAG/build/driver_mechanism_reload.o" "$LIB" \
  -o "$TAG/driver_mechanism_reload"

echo "[2/5] Staging and linking both mechanisms..."
cp test/data/chem.inp test/data/therm.dat "$stage_a/"
cp test/data_plog/chem.inp test/data_plog/therm.dat "$stage_b/"
SC_MECHDIR="$stage_a/" "$TAG/driver_plog_reload_link" >/dev/null
SC_MECHDIR="$stage_b/" "$TAG/driver_plog_reload_link" >/dev/null

echo "[3/5] Reloading with LSODESJAC workspaces..."
reload_output=$(SC_SOLVER=LSODESJAC \
  "$TAG/driver_mechanism_reload" "$stage_a" "$stage_b")

echo "[4/5] Reloading with pointer-owning VODESJAC options..."
vodes_output=$(SC_SOLVER=VODESJAC \
  "$TAG/driver_mechanism_reload" "$stage_a" "$stage_b")

echo "[5/5] Checking exact A reload and complete finalization..."
result=$(grep '^RESULT: PASS' <<<"$reload_output")
vodes_result=$(grep '^RESULT: PASS' <<<"$vodes_output")
if [[ -z "$result" || -z "$vodes_result" ]]; then
  printf '%s\n' "$reload_output" >&2
  printf '%s\n' "$vodes_output" >&2
  echo "run_reload_tests.sh: reload driver did not report PASS" >&2
  exit 1
fi
echo "LSODESJAC: $result"
echo "VODESJAC:  $vodes_result"
echo "======================================================================"
echo " RESULT: PASS - mechanism lifecycle is reusable in one process"
echo "======================================================================"
