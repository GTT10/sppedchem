#!/usr/bin/env bash
#
# Stage-4 integration test against the repo's real LLNL CFD-270 PLOG
# mechanism. This is intentionally separate from the self-contained fast
# suite because it depends on ~/projects/mechanism and Cantera.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root_dir=$(cd -- "$script_dir/.." && pwd)
mechanism_repo=${PLOG_MECH_REPO:-"$HOME/projects/mechanism"}
source_dir="$mechanism_repo/results/zerork_yaml2ck/llnl_hxn_amn_cfd_lumped270"
source_mech="$source_dir/llnl_hxn_amn_cfd_lumped270.mech.txt"
source_therm="$source_dir/llnl_hxn_amn_cfd_lumped270.therm.repaired.txt"
cantera_python=${PLOG_CANTERA_PYTHON:-"$mechanism_repo/.venv-cantera/bin/python"}
ck2yaml=${PLOG_CK2YAML:-"$mechanism_repo/.venv-cantera/bin/ck2yaml"}

for required in "$source_mech" "$source_therm" "$cantera_python" "$ck2yaml"; do
  [[ -f "$required" ]] || {
    echo "run_plog_real_mechanism.sh: missing $required" >&2
    exit 2
  }
done

cd "$root_dir"
stage_dir=$(mktemp -d)
trap 'rm -rf -- "$stage_dir"' EXIT
cp "$source_mech" "$stage_dir/chem.inp"
cp "$source_therm" "$stage_dir/therm.dat"

echo "======================================================================"
echo " SpeedCHEM real PLOG mechanism integration (LLNL CFD-270)"
echo "======================================================================"

echo "[1/6] Building library and real-mechanism drivers..."
make FC=mpiifx -j1 >/dev/null
mpiifx -extend-source 132 -module ifx/mod -O2 \
  -c test/driver_plog.f90 -o "$stage_dir/driver_plog.o"
mpiifx "$stage_dir/driver_plog.o" ifx/libSpeedCHEM64.a \
  -o "$stage_dir/driver_plog"
mpiifx -extend-source 132 -module ifx/mod -O2 \
  -c test/driver_plog_real_history.f90 -o "$stage_dir/driver_history.o"
mpiifx "$stage_dir/driver_history.o" \
  -Wl,--start-group ifx/libSpeedCHEM64.a -Wl,--end-group \
  -o "$stage_dir/driver_history"
mpiifx -extend-source 132 -module ifx/mod -O2 \
  -c test/driver_plog_mpi_real.f90 -o "$stage_dir/driver_mpi.o"
mpiifx "$stage_dir/driver_mpi.o" \
  -Wl,--start-group ifx/libSpeedCHEM64.a -Wl,--end-group \
  -o "$stage_dir/driver_mpi"

echo "[2/6] Parsing the unmodified 272-species / 1715-reaction mechanism..."
SC_MECHDIR="$stage_dir/" "$stage_dir/driver_plog" >"$stage_dir/parse.log"
grep -E '^# linked: 272 species, 1715 reactions$' "$stage_dir/parse.log"
grep -E '^n_plog_reactions 212$' "$stage_dir/parse.log"
grep -E '^n_plog_nodes 1352$' "$stage_dir/parse.log"

echo "[3/6] Converting the exact staged CHEMKIN pair for Cantera..."
"$ck2yaml" --input="$stage_dir/chem.inp" --thermo="$stage_dir/therm.dat" \
  --output="$stage_dir/reference.yaml" --permissive \
  >"$stage_dir/ck2yaml.log" 2>&1

echo "[4/6] Running numerical-Jacobian ignition history..."
SC_MECHDIR="$stage_dir/" SC_SOLVER=LSODES SC_NOUT="${SC_REAL_NOUT:-240}" \
  OMP_NUM_THREADS=1 "$stage_dir/driver_history" >"$stage_dir/numeric.log"
grep '^SUMMARY' "$stage_dir/numeric.log"

analytic_args=()
if [[ ${SC_RUN_HEAVY_ANALYTIC:-0} == 1 ]]; then
  echo "      Running heavy analytical-Jacobian cross-check..."
  SC_MECHDIR="$stage_dir/" SC_SOLVER=LSODESJAC \
    SC_NOUT="${SC_REAL_ANALYTIC_NOUT:-12}" OMP_NUM_THREADS=1 \
    "$stage_dir/driver_history" >"$stage_dir/analytic.log"
  grep '^SUMMARY' "$stage_dir/analytic.log"
  analytic_args=(--analytic "$stage_dir/analytic.log")
fi

echo "[5/6] Comparing IDT, dT/dt, T/P history, and major species..."
PYTHONWARNINGS=ignore "$cantera_python" test/compare_real_plog.py \
  --mechanism "$stage_dir/reference.yaml" \
  --numeric "$stage_dir/numeric.log" "${analytic_args[@]}"

echo "[6/6] Checking real PLOG rates and Jacobian on 2 MPI ranks..."
SC_MECHDIR="$stage_dir/" SC_SOLVER=LSODESJAC OMP_NUM_THREADS=1 \
  mpiexec -n 2 "$stage_dir/driver_mpi"

echo "======================================================================"
echo " RESULT: PASS — LLNL CFD-270 PLOG integration is independently verified"
echo "======================================================================"
