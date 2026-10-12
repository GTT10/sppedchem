#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd -- "$script_dir/.." && pwd)
cd "$repo_dir"

log_file=${LOG_FILE:-ifx/compile_ifx.log}
if [[ "$log_file" != /* ]]; then
  log_file="$repo_dir/$log_file"
fi
mkdir -p -- "$(dirname -- "$log_file")"

fcflags=(
  -c
  -m64
  -extend-source 132
  -convert big_endian
  -warn uninitialized
  -g
  -O3
  -xHost
  -fno-omit-frame-pointer
  -traceback
  -module ifx/mod
)

printf 'Building SpeedCHEM with %s...\n' "${FC:-mpiifx}"

# Rebuild all objects so a previous build cannot retain different flags.
# The Makefile supplies the source list and module dependency order.
make -B all "FC=${FC:-mpiifx}" "FFLAGS=${fcflags[*]}" "$@" \
  2>&1 | tee "$log_file"

printf 'Build finished successfully (log: %s).\n' "$log_file"
