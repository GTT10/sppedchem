#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
make FC=mpiifx -j1
mpiifx -O2 -extend-source 132 -I ifx/mod test/driver_sparse_row_replace.f90 ifx/libSpeedCHEM64.a -o ifx/driver_sparse_row_replace
ifx/driver_sparse_row_replace
