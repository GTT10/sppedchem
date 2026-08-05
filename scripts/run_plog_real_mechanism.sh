#!/usr/bin/env bash
# End-to-end validation against public C3MechV4.0.1 MID 2900.
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
  command -v "$command" >/dev/null 2>&1 || {
    echo "run_plog_real_mechanism.sh: missing command '$command'" >&2
    exit 2
  }
done

"$cantera_python" - <<'PY'
import cantera as ct
print(f"Cantera {ct.__version__}")
PY

stage_dir=$(mktemp -d)
output_dir=${PLOG_OUTPUT_DIR:-"$root_dir/ifx/real-plog-results"}
mkdir -p "$output_dir"

cleanup() {
  local rc=$?
  trap - EXIT
  if [[ $rc -ne 0 ]]; then
    echo "--- public PLOG failure diagnostics ---" >&2
    for file in "$stage_dir"/*.log; do
      [[ -f "$file" ]] || continue
      cp "$file" "$output_dir/" 2>/dev/null || true
      echo "--- $(basename "$file") ---" >&2
      tail -160 "$file" >&2 || true
    done
  fi
  rm -rf -- "$stage_dir"
  exit "$rc"
}
trap cleanup EXIT

run_logged() {
  local logfile=$1
  shift
  if ! "$@" >"$logfile" 2>&1; then
    cat "$logfile" >&2
    return 1
  fi
}

fetch_pinned() {
  local relative_path=$1 output=$2 expected_blob=$3
  curl --fail --location --silent --show-error --retry 4 --retry-all-errors \
    "$source_root/$relative_path" -o "$output"
  local actual_blob
  actual_blob=$(git hash-object "$output")
  [[ "$actual_blob" == "$expected_blob" ]] || {
    printf 'blob mismatch for %s\nexpected=%s\nactual=%s\n' \
      "$relative_path" "$expected_blob" "$actual_blob" >&2
    exit 1
  }
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
echo " source: $source_repo@$source_commit"
echo " model : C3MechV4.0.1 MID 2900 (H2/O2/Ar)"
echo "======================================================================"

read -r expected_ns expected_nr expected_plog expected_terms < <(
  "$cantera_python" - "$stage_dir/published.yaml" <<'PY'
import sys
import cantera as ct
gas = ct.Solution(sys.argv[1])
plog = [r for r in gas.reactions()
        if getattr(r, "reaction_type", "") == "pressure-dependent-Arrhenius"
        or type(r.rate).__name__ == "PlogRate"]
print(gas.n_species, gas.n_reactions, len(plog),
      sum(len(r.rate.rates) for r in plog))
PY
)
printf 'Published metadata: ns=%s nr=%s PLOG=%s terms=%s\n' \
  "$expected_ns" "$expected_nr" "$expected_plog" "$expected_terms"

TAG=ifx
LIB="$TAG/libSpeedCHEM64.a"
DRV_FLAGS=(-extend-source 132 -module "$TAG/mod" -O2)
build_driver() {
  local name=$1
  "$FC" "${DRV_FLAGS[@]}" -c "test/$name.f90" -o "$stage_dir/$name.o"
  "$FC" "$stage_dir/$name.o" -Wl,--start-group "$LIB" \
    -Wl,--end-group -o "$stage_dir/$name"
}

echo "[1/9] Build library and drivers"
make FC="$FC" -j1 >/dev/null
for driver in driver_plog driver_plog_real_rates driver_plog_real_debug \
              driver_plog_real_state_grid driver_plog_real_history \
              driver_plog_mpi_real; do
  build_driver "$driver"
done

echo "[2/9] Parse and round-trip pinned CHEMKIN"
run_logged "$stage_dir/parse.log" env SC_MECHDIR="$stage_dir/" \
  "$stage_dir/driver_plog"
grep -E '^# linked:|^n_plog_reactions|^n_plog_nodes' "$stage_dir/parse.log"
linked_ns=$(awk '/^# linked:/{print $3; exit}' "$stage_dir/parse.log")
linked_nr=$(awk '/^# linked:/{print $5; exit}' "$stage_dir/parse.log")
linked_plog=$(awk '/^n_plog_reactions/{print $2; exit}' "$stage_dir/parse.log")
linked_terms=$(awk '/^n_plog_nodes/{print $2; exit}' "$stage_dir/parse.log")
[[ "$linked_ns $linked_nr $linked_plog $linked_terms" == \
   "$expected_ns $expected_nr $expected_plog $expected_terms" ]] || {
  echo "SpeedCHEM/Cantera topology mismatch" >&2
  exit 1
}

echo "[3/9] Convert exact CHEMKIN pair with Cantera"
run_logged "$stage_dir/ck2yaml.log" "$ck2yaml" \
  --input="$stage_dir/chem.inp" --thermo="$stage_dir/therm.dat" \
  --output="$stage_dir/reference.yaml" --permissive
"$cantera_python" - "$stage_dir/reference.yaml" <<'PY'
import sys
import cantera as ct
gas = ct.Solution(sys.argv[1])
print(f"Regenerated: {gas.n_species} species, {gas.n_reactions} reactions")
PY

echo "[4/9] Evaluate PLOG rate grid"
run_logged "$stage_dir/rates.log" env SC_MECHDIR="$stage_dir/" \
  "$stage_dir/driver_plog_real_rates"
grep '^# REAL_PLOG_RATES' "$stage_dir/rates.log"

echo "[5/9] Decompose the real PLOG forward/reverse reaction path"
run_logged "$stage_dir/debug.log" env SC_MECHDIR="$stage_dir/" \
  "$stage_dir/driver_plog_real_debug"
grep '^DEBUG' "$stage_dir/debug.log"

echo "[6/9] Evaluate complete constant-volume RHS grid"
run_logged "$stage_dir/states.log" env SC_MECHDIR="$stage_dir/" \
  "$stage_dir/driver_plog_real_state_grid"
printf 'State rows: %s\n' "$(grep -c '^STATE,' "$stage_dir/states.log")"

echo "[7/9] Integrate numeric and analytic-Jacobian histories"
real_t0=${SC_REAL_T0:-1200.0}
real_p0=${SC_REAL_P0:-1013250.0}
real_tend=${SC_REAL_TEND:-0.002}
real_nout=${SC_REAL_NOUT:-400}
common_env=(SC_MECHDIR="$stage_dir/" SC_T0="$real_t0" SC_P0="$real_p0"
            SC_TEND="$real_tend" SC_NOUT="$real_nout" OMP_NUM_THREADS=1)
run_logged "$stage_dir/numeric.log" env "${common_env[@]}" SC_SOLVER=LSODES \
  "$stage_dir/driver_plog_real_history"
run_logged "$stage_dir/analytic.log" env "${common_env[@]}" SC_SOLVER=LSODESJAC \
  "$stage_dir/driver_plog_real_history"
grep '^SUMMARY' "$stage_dir/numeric.log"
grep '^SUMMARY' "$stage_dir/analytic.log"

echo "[8/9] Compare PLOG rates/ROP and solver histories with Cantera"
run_logged "$stage_dir/comparison.log" env PYTHONWARNINGS=ignore \
  "$cantera_python" test/compare_real_plog.py \
  --mechanism "$stage_dir/reference.yaml" --published "$stage_dir/published.yaml" \
  --rates "$stage_dir/rates.log" --debug "$stage_dir/debug.log" \
  --states "$stage_dir/states.log" --numeric "$stage_dir/numeric.log" \
  --analytic "$stage_dir/analytic.log"
cat "$stage_dir/comparison.log"

echo "[9/9] Check two-rank MPI rate/RHS/Jacobian equality"
run_logged "$stage_dir/mpi.log" env SC_MECHDIR="$stage_dir/" \
  SC_SOLVER=LSODESJAC OMP_NUM_THREADS=1 mpiexec -n 2 \
  "$stage_dir/driver_plog_mpi_real"
cat "$stage_dir/mpi.log"

cp "$stage_dir"/*.log "$output_dir/"
echo "======================================================================"
echo " RESULT: PASS - public C3Mech PLOG independently verified"
echo " Results: $output_dir"
echo "======================================================================"
