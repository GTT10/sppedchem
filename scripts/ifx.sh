#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$script_dir"

SRCDIR=../src
SCLIB=${SCLIB:-../ifx/libSpeedCHEM64.a}
MODULEDIR=../ifx/mod
LOG_FILE=${LOG_FILE:-../ifx/compile_ifx.log} # ログファイル名を変更

mkdir -p ../ifx
mkdir -p "$MODULEDIR"

# Intel MPIラッパー(mpiifx)を使用するため, 手動のMPIパス探索は不要である

FCFLAGS=(
  -c
  -m64
  -extend-source 132         # -ffixed-line-length-none の代替
  -convert big_endian        # -fconvert=big-endian の代替
  -warn uninitialized        # -Wuninitialized の代替
  -g                         # フルデバッグシンボル
  -O3                        # 最適化レベル3
  -xHost                     # 実行ホストのCPUアーキテクチャへの最適化(推奨)
  -fno-omit-frame-pointer
  -traceback                 # -fbacktrace の代替
  -module "$MODULEDIR"       # -J の代替
  # -fallow-argument-mismatch に相当するオプションはifxでは通常不要である.
  # 厳密な型チェックでエラーが出る場合は -diag-disable=8284 等で個別に警告を抑止する.
)

SRCS=(
  "$SRCDIR/working_precision.f90"
  "$SRCDIR/SCutilities.f"
  "$SRCDIR/SCsparse_definitions.f90"
  "$SRCDIR/SCsparse.f"
  "$SRCDIR/sparse_MPI.f"
  "$SRCDIR/dvode_f90_m.f90"
  "$SRCDIR/SCmodule.f"
  "$SRCDIR/chemkin_module.f90"
  "$SRCDIR/SCconV.f90"
  "$SRCDIR/SCsetup.f"
  "$SRCDIR/SCallocate.f"
  "$SRCDIR/SCcklink.f90"
  "$SRCDIR/gam.f"
  "$SRCDIR/gamsub.f"
  "$SRCDIR/opkdmain.f"
  "$SRCDIR/opkda1.f"
  "$SRCDIR/opkda2.f"
  "$SRCDIR/ddaspk.f"
  "$SRCDIR/rodas.f"
  "$SRCDIR/vode.f"
  "$SRCDIR/MEBDFSO.f"
  "$SRCDIR/rowmap.f"
  "$SRCDIR/radau5.f"
  "$SRCDIR/radaua.f"
  "$SRCDIR/chemistry_input.f90"
  "$SRCDIR/SCbroadcast.f90"
)

rm -f "$SCLIB" ./*.mod

printf 'Compiling to generate %s...\n' "$SCLIB"

{
  # gfortranの代わりにmpiifxを使用する
  mpiifx "${FCFLAGS[@]}" -c "${SRCS[@]}"
  ar cr "$SCLIB" ./*.o
  ranlib "$SCLIB"
} >"$LOG_FILE" 2>&1

rm -f ./*.o

if [[ -f "$SCLIB" ]]; then
  printf 'Compilation finished successfully (see %s).\n' "$LOG_FILE"
else
  printf 'Compilation failed. Check %s for details.\n' "$LOG_FILE" >&2
  exit 1
fi