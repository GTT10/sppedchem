# SpeedCHEM regression tests

The canonical build and test entry points use Intel `ifx` through the
`mpiifx` wrapper.

```bash
make test-plog       # PLOG parse/cklink, rates, RHS, Jacobian, integration, MPI, negatives
make test            # smoke + PLOG + reload + long-string regressions
make test-openmp     # isolated OpenMP/bounds/reload regression
make test-all        # all bundled self-contained tests
```

`make test-real-plog` is intentionally separate. It requires the external
LLNL CFD-270 mechanism repository and a Cantera environment configured by
`scripts/run_plog_real_mechanism.sh`.

The GitHub Actions workflow in `.github/workflows/ci.yml` installs the Intel
Fortran compiler and Intel MPI from Intel's APT repository, then runs
`make test-all` on pull requests.
