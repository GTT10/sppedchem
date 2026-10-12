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

printf 'Building SpeedCHEM...\n'

# Rebuild all objects so a previous build cannot retain different flags.
# The Makefile supplies compiler options, sources, and module dependency order.
make -B all "$@" \
  2>&1 | tee "$log_file"

printf 'Build finished successfully (log: %s).\n' "$log_file"
