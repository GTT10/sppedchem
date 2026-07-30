#!/usr/bin/env bash
#
# Build an isolated OpenMP + bounds-checking library and verify that the
# mechanism lifecycle remains reusable with two threads. The canonical ifx/
# build is not touched.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root_dir=$(cd -- "$script_dir/.." && pwd)
cd "$root_dir"

FC=mpiifx
if ! command -v "$FC" >/dev/null 2>&1; then
  echo "run_openmp_reload_tests.sh: '$FC' not found" >&2
  exit 2
fi

test_root=$(mktemp -d)
trap 'rm -rf -- "$test_root"' EXIT

omp_out="$test_root/ifx"
build_log="$test_root/build.log"
fflags="-c -extend-source 132 -module $omp_out/mod -O0 -g -check bounds -traceback -qopenmp"
driver_flags=(-extend-source 132 -module "$omp_out/mod" -O0 -g \
              -check bounds -traceback -qopenmp)

echo "======================================================================"
echo " SpeedCHEM OpenMP bounds/reload test"
echo "======================================================================"

echo "[1/5] Building isolated OpenMP library with bounds checking..."
if ! make OUTDIR="$omp_out" FC="$FC" FFLAGS="$fflags" -j1 \
     >"$build_log" 2>&1; then
  tail -80 "$build_log" >&2
  exit 1
fi

echo "[2/5] Building link and lifecycle drivers..."
"$FC" "${driver_flags[@]}" -c test/driver_plog.f90 \
  -o "$test_root/driver_plog.o"
"$FC" -qopenmp "$test_root/driver_plog.o" \
  "$omp_out/libSpeedCHEM64.a" -o "$test_root/driver_plog"
"$FC" "${driver_flags[@]}" -c test/driver_mechanism_reload.f90 \
  -o "$test_root/driver_reload.o"
"$FC" -qopenmp "$test_root/driver_reload.o" \
  "$omp_out/libSpeedCHEM64.a" -o "$test_root/driver_reload"

echo "[3/5] Staging and linking non-PLOG A and PLOG B..."
mkdir "$test_root/a" "$test_root/b"
cp test/data/chem.inp test/data/therm.dat "$test_root/a/"
cp test/data_plog/chem.inp test/data_plog/therm.dat "$test_root/b/"
SC_MECHDIR="$test_root/a/" "$test_root/driver_plog" >/dev/null
SC_MECHDIR="$test_root/b/" "$test_root/driver_plog" >/dev/null

echo "[4/5] Running A -> B -> A with two OpenMP threads..."
OMP_NUM_THREADS=2 SC_SOLVER=LSODESJAC \
  "$test_root/driver_reload" "$test_root/a" "$test_root/b" \
  >"$test_root/lsodes.log" 2>"$test_root/lsodes.err"
OMP_NUM_THREADS=2 SC_SOLVER=VODESJAC \
  "$test_root/driver_reload" "$test_root/a" "$test_root/b" \
  >"$test_root/vodes.log" 2>"$test_root/vodes.err"

echo "[5/5] Checking bounds diagnostics and exact reload results..."
if grep -qi 'severe' "$test_root/lsodes.err" "$test_root/vodes.err"; then
  tail -80 "$test_root/lsodes.err" "$test_root/vodes.err" >&2
  exit 1
fi
lsodes_result=$(grep '^RESULT: PASS' "$test_root/lsodes.log")
vodes_result=$(grep '^RESULT: PASS' "$test_root/vodes.log")
echo "LSODESJAC: $lsodes_result"
echo "VODESJAC:  $vodes_result"
echo "======================================================================"
echo " RESULT: PASS - OpenMP build, bounds checks, and 2-thread reload"
echo "======================================================================"
