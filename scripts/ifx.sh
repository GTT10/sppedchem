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
  -extend-source 132         # 固定形式の行長制限を緩和
  -convert big_endian        # ビッグエンディアンのバイナリ入出力
  -warn uninitialized        # 未初期化変数の警告
  -g                         # フルデバッグシンボル
  -O3                        # 最適化レベル3
  -xHost                     # 実行ホストのCPUアーキテクチャへの最適化(推奨)
  -fno-omit-frame-pointer
  -traceback                 # 実行時エラーでバックトレース
  -module "$MODULEDIR"       # モジュール(.mod)出力先
  # 厳密な型チェックでエラーが出る場合は -diag-disable=8284 等で個別に警告を抑止する.
)

SRCS=(
  "$SRCDIR/working_precision.f90"
  "$SRCDIR/SCstring_limits.f90"
  "$SRCDIR/SCutilities.f90"
  "$SRCDIR/SCsparse_definitions.f90"
  "$SRCDIR/SCsparse.f90"
  "$SRCDIR/sparse_MPI.f90"
  "$SRCDIR/dvode_f90_m.f90"
  "$SRCDIR/SCmodule.f90"
  "$SRCDIR/chemkin_module.f90"
  "$SRCDIR/SCconV.f90"
  "$SRCDIR/SCsetup.f90"
  "$SRCDIR/SCallocate.f90"
  "$SRCDIR/SCcklink.f90"
  "$SRCDIR/gam.f90"
  "$SRCDIR/gamsub.f90"
  "$SRCDIR/opkdmain.f90"
  "$SRCDIR/opkda1.f90"
  "$SRCDIR/opkda2.f90"
  "$SRCDIR/ddaspk.f90"
  "$SRCDIR/rodas.f90"
  "$SRCDIR/vode.f90"
  "$SRCDIR/MEBDFSO.f90"
  "$SRCDIR/rowmap.f90"
  "$SRCDIR/radau5.f90"
  "$SRCDIR/radaua.f90"
  "$SRCDIR/chemistry_input.f90"
  "$SRCDIR/SCbroadcast.f90"
)

rm -f "$SCLIB" ./*.mod

printf 'Compiling to generate %s...\n' "$SCLIB"

{
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
