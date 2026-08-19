#!/usr/bin/env bash
# Regression test for issue #5: crossing the fixed CKINTP IDIM limit must
# report the reaction-capacity error immediately, before later PLOG records
# can be attached to reaction IDIM, and must leave no cklink behind.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root_dir=$(cd -- "$script_dir/.." && pwd)
cd "$root_dir"

TAG=ifx
FC=mpiifx
DRV_FLAGS=(-extend-source 132 -module "$TAG/mod" -O2)

if ! command -v "$FC" >/dev/null 2>&1; then
  echo "run_capacity_tests.sh: required compiler wrapper '$FC' not found in PATH" >&2
  exit 2
fi

echo "======================================================================"
echo " SpeedCHEM CKINTP reaction-capacity regression"
echo "   compiler : $FC"
echo "======================================================================"

make FC="$FC" -j1 >/dev/null
LIB="$TAG/libSpeedCHEM64.a"
[[ -f "$LIB" ]] || { echo "library not produced: $LIB" >&2; exit 1; }
DRV_BIN="$TAG/driver_capacity"
"$FC" "${DRV_FLAGS[@]}" -c test/driver_smoke.f90 -o "$TAG/build/driver_capacity.o"
"$FC" "$TAG/build/driver_capacity.o" "$LIB" -o "$DRV_BIN"

tmp_root=$(mktemp -d)
trap 'rm -rf -- "$tmp_root"' EXIT
mechdir="$tmp_root/reaction_capacity"
mkdir -p "$mechdir"

cat > "$mechdir/chem.inp" <<'EOF'
ELEMENTS
H O
END
SPECIES
H2 H O OH
END
REACTIONS CAL/MOLE
EOF

# IDIM is 9000. The next reaction must cause an immediate primary failure.
for ((i=1; i<=9001; i++)); do
  printf '%s\n' 'H2+O=H+OH  1.000E+13  0.0  0.0' >> "$mechdir/chem.inp"
done

# With the old behavior these lines were incorrectly assigned to reaction
# 9000, replacing the capacity diagnostic with a PLOG-order error.
cat >> "$mechdir/chem.inp" <<'EOF'
PLOG / 1.0 1.000E+13 0.0 0.0 /
PLOG / 0.1 1.000E+13 0.0 0.0 /
END
EOF

cp "$root_dir/test/data_plog/therm.dat" "$mechdir/therm.dat"

set +e
out=$(SC_MECHDIR="$mechdir/" "$DRV_BIN" 2>&1)
rc=$?
set -e

want="reaction count exceeds CKINTP capacity IDIM=9000"
if [[ $rc -eq 0 ]]; then
  echo "FAIL: driver exited 0 after exceeding IDIM"
  printf '%s\n' "$out"
  exit 1
fi
if ! grep -qF "$want" <<<"$out"; then
  echo "FAIL: primary reaction-capacity diagnostic not found"
  printf '%s\n' "$out"
  exit 1
fi
if grep -qF "PLOG pressures for reaction" <<<"$out"; then
  echo "FAIL: secondary PLOG-order diagnostic masked the capacity error"
  printf '%s\n' "$out"
  exit 1
fi
if [[ -e "$mechdir/cklink" ]]; then
  echo "FAIL: cklink exists after failed parse"
  exit 1
fi

echo "RESULT: PASS - overflow refused before PLOG parsing; no cklink remains"
echo "======================================================================"
