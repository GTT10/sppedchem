# SpeedCHEM

[![CI](https://github.com/GTT10/sppedchem/actions/workflows/ci.yml/badge.svg)](https://github.com/GTT10/sppedchem/actions/workflows/ci.yml)

Fortran library for stiff gas-phase chemical kinetics with sparse analytical
Jacobians, constant-volume integration, CHEMKIN/NASA input, PLOG, and MPI.

This is a community-maintained version of [Federico Perini’s SpeedCHEM](https://www.federicoperini.info/speedchem), not an official upstream release.

## Build

Requires Linux, Intel oneAPI `ifx`/`mpiifx`, Intel MPI, and GNU Make.

```bash
source /opt/intel/oneapi/setvars.sh
make
```

Outputs: `ifx/libSpeedCHEM64.a` and `ifx/mod/*.mod`.

## Test

```bash
make test-all          # Bundled regressions, MPI, and OpenMP
```

For the pinned public C3Mech comparison against Cantera:

```bash
python -m pip install "cantera==3.2.0"
make test-real-plog
```

See [testing details](docs/testing.md) and [supported reaction forms](CLAUDE.md#plog-support--cklink-v2).
Unsupported forms are rejected.

## Use

1. Set `chemistry_setup::mechdir` to the mechanism directory.
2. Call `chemistry_input` to initialize.
3. Call `chemistry_ODE_integrate` for constant-volume integration.
4. Call `chemistry_finalize` before loading another mechanism.

Examples: [test drivers](test/). Internals: [CLAUDE.md](CLAUDE.md).

## Citation and license

Cite the [SpeedCHEM methods paper](https://doi.org/10.1021/ef300747n).
The [direct/Krylov solver paper](https://doi.org/10.1016/j.combustflame.2013.11.017)
and mechanism references are listed in [CITATION.cff](CITATION.cff).

SpeedCHEM-authored code: [GPL-3.0-or-later](LICENSE).
Third-party code and test-data terms, including unresolved mechanism redistribution:
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
Contributions: [CONTRIBUTING.md](CONTRIBUTING.md).
