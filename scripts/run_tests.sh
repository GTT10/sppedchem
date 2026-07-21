#!/usr/bin/env bash
#
# SpeedCHEM build-and-test harness.
#
#   1. Builds libSpeedCHEM64.a for the chosen compiler (via the top-level Makefile).
#   2. Compiles and links test/driver_smoke.f90 against it.
#   3. Runs the smoke test against a real mechanism directory and gates on its exit code.
#
# The library uses MPI (SCbroadcast.f90 / sparse_MPI.f), so the driver is compiled and
# linked with the MPI compiler wrapper (mpiifx / mpif90), matching how the library is built.
#
# Usage:
#   scripts/run_tests.sh                 # ifx (default; mpiifx) — matches "ifx-first" policy
#   scripts/run_tests.sh gfortran        # gfortran (mpif90)
#   scripts/run_tests.sh ifx             # explicit ifx
#   SC_MECHDIR=/path/to/data/ scripts/run_tests.sh   # override mechanism directory
#
# Exit code: 0 = all passed, non-zero = build or test failure.

set -euo pipefail

# ---- Locations --------------------------------------------------------------------------
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root_dir=$(cd -- "$script_dir/.." && pwd)
cd "$root_dir"

# ---- Compiler selection -----------------------------------------------------------------
COMPILER=${1:-ifx}

case "$COMPILER" in
  ifx)
    TAG=ifx
    LIB_FC=mpiifx
    DRV_FC=mpiifx
    DRV_FLAGS=(-extend-source 132 -module "$TAG/mod" -O2)
    ;;
  gfortran|gfort)
    TAG=gfortran
    LIB_FC=mpif90
    DRV_FC=mpif90
    DRV_FLAGS=(-ffree-line-length-none -I"$TAG/mod" -J"$TAG/mod" -O2 -fallow-argument-mismatch)
    ;;
  *)
    echo "run_tests.sh: unknown compiler '$COMPILER' (use 'ifx' or 'gfortran')" >&2
    exit 2
    ;;
esac

if ! command -v "$LIB_FC" >/dev/null 2>&1; then
  echo "run_tests.sh: required compiler wrapper '$LIB_FC' not found in PATH" >&2
  exit 2
fi

# ---- Mechanism directory ----------------------------------------------------------------
# Default is the bundled test mechanism (test/data/); override with SC_MECHDIR to point at
# another dir (e.g. production data). Must end with a separator.
export SC_MECHDIR=${SC_MECHDIR:-$root_dir/test/data/}
case "$SC_MECHDIR" in
  */) ;;
  *) SC_MECHDIR="$SC_MECHDIR/" ;;
esac

echo "======================================================================"
echo " SpeedCHEM test harness"
echo "   compiler : $COMPILER  (lib: $LIB_FC, driver: $DRV_FC)"
echo "   mechdir  : $SC_MECHDIR"
echo "======================================================================"

# Sanity-check the input data up front so a missing mount fails clearly, not deep in Fortran I/O.
for f in chem.inp therm.dat; do
  if [[ ! -r "${SC_MECHDIR}${f}" ]]; then
    echo "run_tests.sh: cannot read ${SC_MECHDIR}${f}" >&2
    echo "  (set SC_MECHDIR to a directory containing chem.inp and therm.dat)" >&2
    exit 2
  fi
done

# ---- 1. Build the library ---------------------------------------------------------------
echo "[1/3] Building libSpeedCHEM64.a ($TAG)..."
make FC="$LIB_FC" COMPILER_TAG="$TAG" -j1

LIB="$TAG/libSpeedCHEM64.a"
[[ -f "$LIB" ]] || { echo "run_tests.sh: library not produced ($LIB)" >&2; exit 1; }

# ---- 2. Build the driver ----------------------------------------------------------------
echo "[2/3] Compiling and linking smoke-test driver..."
DRV_BIN="$TAG/driver_smoke"
"$DRV_FC" "${DRV_FLAGS[@]}" -c test/driver_smoke.f90 -o "$TAG/build/driver_smoke.o"
"$DRV_FC" "$TAG/build/driver_smoke.o" "$LIB" -o "$DRV_BIN"

# ---- 3. Run ------------------------------------------------------------------------------
echo "[3/3] Running smoke test..."
echo "----------------------------------------------------------------------"
# The driver writes SpeedCHEM.out etc. into mechdir; run from there so relative outputs land
# alongside the data (and so any bare-filename I/O resolves as the production code expects).
if "$root_dir/$DRV_BIN"; then
  echo "======================================================================"
  echo " RESULT: PASS ($COMPILER)"
  echo "======================================================================"
  exit 0
else
  rc=$?
  echo "======================================================================"
  echo " RESULT: FAIL ($COMPILER) — driver exit code $rc"
  echo "======================================================================"
  exit "$rc"
fi
