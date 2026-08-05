#!/usr/bin/env bash
#
# End-to-end validation against a public, published PLOG mechanism:
# C3MechV4.0.1 MID 2900 (H2/O2 core, CC BY 4.0), pinned to an immutable
# upstream commit. The exact same CHEMKIN pair is evaluated by SpeedCHEM and
# converted by Cantera for independent rate, RHS, ignition-history, solver,
# and MPI comparisons.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root_dir=$(cd -- "$script_dir/.." && pwd)
cd "$root_dir"

FC=${PLOG_FC:-mpiifx}
cantera_python=${PLOG_CANTERA_PYTHON:-python3}
ck2yaml=${PLOG_CK2YAML:-ck2yaml}

source_repo=C3Mech/C3Mech
source_commit=e9cae8ac06198234300adb8f99e1b3c2e8f65f19
source_root="https://raw.githubusercontent.com/${source_repo}/${source_commit}"
source_mech_path=PRECOMPILED/C0/C0/Cantera/C3MechV4.0.1_2900_C0_Cantera.CKI
source_therm_path=PRECOMPILED/C0/C0/Chemkin/C3MechV4.0.1_2900_C0.THERM
source_yaml_path=PRECOMPILED/C0/C0/Cantera/C3MechV4.0.1_2900_C0_Cantera.yaml
source_mech_blob=a91ab6d11018a4b741680fe420037d35f4f23da5
source_therm_blob=7533a72d9e1c75cb65c32585ffaca5a3ac3481b4
source_yaml_blob=5cb671203bd6a381bd9a4f4bec1ffb990d40c366

for command in "$FC" mpiexec curl git "$cantera_python" "$ck2yaml"; do
  if ! command -v "$command" >/dev/null 2>&1; then
    echo "run_plog_real_mechanism.sh: required command '$command' not found" >&2
    exit 2
  fi
done

"$cantera_python" - <<'PY'
import cantera as ct
print(f"Cantera {ct.__version__}")
PY

stage_dir=$(mktemp -d)
output_dir=${PLOG_OUTPUT_DIR:-"$root_dir/ifx/real-plog-results"}
mkdir -p "$output_dir"
trap 'rm -rf -- "$stage_dir"' EXIT

fetch_pinned() {
  local relative_path=$1
  local output=$2
  local expected_blob=$3
  curl --fail --location --silent --show-error --retry 4 --retry-all-errors \
    "$source_root/$relative_path" -o "$output"
  local actual_blob
  actual_blob=$(git hash-object "$output")
  if [[ "$actual_blob" != "$expected_blob" ]]; then
    echo "run_plog_real_mechanism.sh: blob mismatch for $relative_path" >&2
    echo "  expected: $expected_blob" >&2
    echo "  actual:   $actual_blob" >&2
    exit 1
  fi
}

fetch_pinned "$source_mech_path"  "$stage_dir/chem.inp"       "$source_mech_blob"
fetch_pinned "$source_therm_path" "$stage_dir/therm.dat"      "$source_therm_blob"
fetch_pinned "$source_yaml_path"  "$stage_dir/published.yaml" "$source_yaml_blob"

cat >"$output_dir/SOURCE.txt" <<EOF
repository=$source_repo
commit=$source_commit
chemkin_path=$source_mech_path
chemkin_blob=$source_mech_blob
thermo_path=$source_therm_path
thermo_blob=$source_therm_blob
published_yaml_path=$source_yaml_path
published_yaml_blob=$source_yaml_blob
EOF

echo "======================================================================"
echo " SpeedCHEM public real-PLOG validation"
echo "   source   : $source_repo@$source_commit"
echo "   model    : C3MechV4.0.1 MID 2900"
echo "   chemistry: H2/O2/Ar with 0.01--300 atm PLOG fit"
echo "======================================================================"

read -r expected_ns expected_nr expected_plog expected_terms < <(
  "$cantera_python" - "$stage_dir/published.yaml" <<'PY'
import sys
import cantera as ct

gas = ct.Solution(sys.argv[1])
plog = [r for r in gas.reactions()
        if getattr(r, "reaction_type", "") == "pressure-dependent-Arrhenius"
        or type(r.rate).__name__ == "PlogRate"]
terms = sum(len(r.rate.rates) for r in plog)
print(gas.n_species, gas.n_reactions, len(plog), terms)
PY
)
printf 'Published Cantera metadata: ns=%s nr=%s PLOG=%s terms=%s\n' \
  "$expected_ns" "$expected_nr" "$expected_plog" "$expected_terms"

TAG=ifx
LIB="$TAG/libSpeedCHEM64.a"
DRV_FLAGS=(-extend-source 132 -module "$TAG/mod" -O2)

echo "[1/8] Building SpeedCHEM and public-mechanism drivers..."
make FC="$FC" -j1 >/dev/null
[[ -f "$LIB" ]] || { echo "library not produced: $LIB" >&2; exit 1; }

build_driver() {
  local name=$1
  "$FC" "${DRV_FLAGS[@]}" -c "test/$name.f90" -o "$stage_dir/$name.o"
  "$FC" "$stage_dir/$name.o" -Wl,--start-group "$LIB" \
    -Wl,--end-group -o "$stage_dir/$name"
}

build_driver driver_plog
build_driver driver_plog_real_rates
build_driver driver_plog_real_state_grid
build_driver driver_plog_real_history
build_driver driver_plog_mpi_real

echo "[2/8] Parsing and round-tripping the pinned CHEMKIN mechanism..."
SC_MECHDIR="$stage_dir/" "$stage_dir/driver_plog" >"$stage_dir/parse.log"
cat "$stage_dir/parse.log" | grep -E '^# linked:|^n_plog_reactions|^n_plog_nodes'
linked_ns=$(awk '/^# linked:/{print $3; exit}' "$stage_dir/parse.log")
linked_nr=$(awk '/^# linked:/{print $5; exit}' "$stage_dir/parse.log")
linked_plog=$(awk '/^n_plog_reactions/{print $2; exit}' "$stage_dir/parse.log")
linked_terms=$(awk '/^n_plog_nodes/{print $2; exit}' "$stage_dir/parse.log")
if [[ "$linked_ns" != "$expected_ns" || "$linked_nr" != "$expected_nr" || \
      "$linked_plog" != "$expected_plog" || "$linked_terms" != "$expected_terms" ]]; then
  echo "SpeedCHEM/Cantera mechanism topology mismatch" >&2
  echo "  SpeedCHEM: $linked_ns $linked_nr $linked_plog $linked_terms" >&2
  echo "  Cantera:   $expected_ns $expected_nr $expected_plog $expected_terms" >&2
  exit 1
fi

echo "[3/8] Regenerating Cantera YAML from the exact staged CHEMKIN pair..."
"$ck2yaml" --input="$stage_dir/chem.inp" --thermo="$stage_dir/therm.dat" \
  --output="$stage_dir/reference.yaml" --permissive \
  >"$stage_dir/ck2yaml.log" 2>&1
"$cantera_python" - "$stage_dir/reference.yaml" <<'PY'
import sys
import cantera as ct

gas = ct.Solution(sys.argv[1])
print(f"Regenerated Cantera mechanism: {gas.n_species} species, {gas.n_reactions} reactions")
PY

echo "[4/8] Comparing every real PLOG rate over T/P nodes, intervals, and clamps..."
SC_MECHDIR="$stage_dir/" "$stage_dir/driver_plog_real_rates" \
  >"$stage_dir/rates.log"
grep '^# REAL_PLOG_RATES' "$stage_dir/rates.log"

echo "[5/8] Evaluating the complete constant-volume RHS over a 33-state grid..."
SC_MECHDIR="$stage_dir/" "$stage_dir/driver_plog_real_state_grid" \
  >"$stage_dir/states.log"
printf 'State rows: %s\n' "$(grep -c '^STATE,' "$stage_dir/states.log")"

echo "[6/8] Integrating H2/O2/Ar ignition with numeric and analytic Jacobians..."
real_t0=${SC_REAL_T0:-1200.0}
real_p0=${SC_REAL_P0:-1013250.0}
real_tend=${SC_REAL_TEND:-0.002}
real_nout=${SC_REAL_NOUT:-400}
common_history_env=(
  SC_MECHDIR="$stage_dir/"
  SC_T0="$real_t0"
  SC_P0="$real_p0"
  SC_TEND="$real_tend"
  SC_NOUT="$real_nout"
  OMP_NUM_THREADS=1
)
env "${common_history_env[@]}" SC_SOLVER=LSODES \
  "$stage_dir/driver_plog_real_history" >"$stage_dir/numeric.log"
env "${common_history_env[@]}" SC_SOLVER=LSODESJAC \
  "$stage_dir/driver_plog_real_history" >"$stage_dir/analytic.log"
grep '^SUMMARY' "$stage_dir/numeric.log"
grep '^SUMMARY' "$stage_dir/analytic.log"

echo "[7/8] Comparing rates, full RHS, ignition history, and solver paths with Cantera..."
PYTHONWARNINGS=ignore "$cantera_python" test/compare_real_plog.py \
  --mechanism "$stage_dir/reference.yaml" \
  --published "$stage_dir/published.yaml" \
  --rates "$stage_dir/rates.log" \
  --states "$stage_dir/states.log" \
  --numeric "$stage_dir/numeric.log" \
  --analytic "$stage_dir/analytic.log" | tee "$stage_dir/comparison.log"

echo "[8/8] Checking PLOG rate, RHS, and analytic Jacobian on two MPI ranks..."
SC_MECHDIR="$stage_dir/" SC_SOLVER=LSODESJAC OMP_NUM_THREADS=1 \
  mpiexec -n 2 "$stage_dir/driver_plog_mpi_real" | tee "$stage_dir/mpi.log"

cp "$stage_dir/parse.log" "$stage_dir/ck2yaml.log" \
   "$stage_dir/rates.log" "$stage_dir/states.log" \
   "$stage_dir/numeric.log" "$stage_dir/analytic.log" \
   "$stage_dir/comparison.log" "$stage_dir/mpi.log" \
   "$output_dir/"

echo "======================================================================"
echo " RESULT: PASS - public C3Mech PLOG mechanism independently verified"
echo " Results: $output_dir"
echo "======================================================================"
