#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$script_dir"

SRCDIR=../src
SCLIB=${SCLIB:-../gfortran/libSpeedCHEM64.a}
MODULEDIR=../gfortran/mod
LOG_FILE=${LOG_FILE:-../gfortran/compile_gft.log}

mkdir -p ../gfortran
mkdir -p "$MODULEDIR"

if [[ -n "${GNUHOME:-}" ]]; then
  export PATH="$GNUHOME/bin:$PATH"
fi

detect_mpi_include_dir() {
  if [[ -n "${MPIHOME:-}" && -d "$MPIHOME/include" ]]; then
    printf '%s\n' "$MPIHOME/include"
    return 0
  fi

  if command -v mpif90 >/dev/null 2>&1; then
    if incdirs=$(mpif90 -showme:incdirs 2>/dev/null); then
      for dir in $incdirs; do
        if [[ -d "$dir" ]]; then
          printf '%s\n' "$dir"
          return 0
        fi
      done
    fi

    if show=$(mpif90 -show 2>/dev/null); then
      for token in $show; do
        if [[ $token == -I* ]]; then
          dir=${token#-I}
          if [[ -d "$dir" ]]; then
            printf '%s\n' "$dir"
            return 0
          fi
        fi
      done
    fi
  fi

  return 1
}

MPI_INC_DEFAULT=""
if auto_inc=$(detect_mpi_include_dir); then
  MPI_INC_DEFAULT="$auto_inc"
fi
MPI_INC=${MPI_INC:-$MPI_INC_DEFAULT}

FCFLAGS=(
  -c
  -malign-double
  -m64
  -ffixed-line-length-none
  -ffree-line-length-none
  -fconvert=big-endian
  -Wuninitialized
  -g                          # Full debug symbols for valgrind
  -O3                         # Optimization level 3 for performance
#  -fcheck=all                 # Runtime checks (bounds, pointers, etc.)
  -fno-omit-frame-pointer     # Keep frame pointers for better stack traces
  -fbacktrace                 # Enable backtrace on runtime errors
  -fdump-core                 # Generate core dump on crash
  "-J$MODULEDIR"              # Module output directory
  -fallow-argument-mismatch  # Allow argument mismatch (for some libraries)
)

if [[ -n "$MPI_INC" ]]; then
  FCFLAGS+=("-I$MPI_INC")
fi

if [[ -n "$MPI_INC" ]]; then
  printf 'Using MPI include directory: %s\n' "$MPI_INC"
fi

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
  gfortran "${FCFLAGS[@]}" -c "${SRCS[@]}"
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
