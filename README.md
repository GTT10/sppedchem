# SpeedCHEM

SpeedCHEM is a Fortran library for gas-phase chemical kinetics.
It integrates stiff species and temperature equations in
constant-volume reactors and can be embedded in combustion CFD codes.

Sparse analytical Jacobians reduce the cost of solving large
reaction mechanisms.

## Features

- CHEMKIN mechanism and NASA thermodynamic input
- Sparse analytical Jacobians
- Multiple stiff ODE/DAE solvers
- Pressure-dependent reactions, including PLOG
- MPI support for mechanism initialization

See [supported input](#supported-input)
for supported reaction forms and limitations.

## Usage

Set the mechanism directory, initialize the library, and call
the integration routine from your application:

1. Set `chemistry_setup::mechdir`.
2. Call `chemistry_input`.
3. Call `chemistry_ODE_integrate`.
4. Call `chemistry_finalize` when finished.

## Supported input

CHEMKIN mechanism files and NASA-7 thermodynamic records support
18-character species names and 256-character reaction records.
Thermodynamic lookup accepts full-width identifiers with adjacent dates
and compact headers for short species names.

PLOG requires positive Arrhenius A factors and non-decreasing pressures.
Terms at the same pressure are summed before logarithmic pressure
interpolation; pressures outside the table use the nearest endpoint.
PLOG combined with explicit `REV`, third-body/falloff, `FORD/RORD`, or
non-integer stoichiometry is rejected.

The analytical Jacobian includes pressure dependence through temperature
and composition. PLOG disables simplified sparsity.

Activation energies use `8.31446261815324 J/mol/K` and `4.184 J/cal`.
These conversions change rates relative to older versions that used
inconsistent constants. Regenerate `cklink` when updating from those
versions or changing compiler flags or binary byte order.

## About this version

This repository maintains Federico Perini's SpeedCHEM,
with added PLOG support and improvements to mechanism loading
and validation.

It is not an official upstream release.

## Building

The supplied [Makefile](Makefile) uses Intel oneAPI and Intel MPI.
With the oneAPI environment initialized, build with:

    make

Outputs are `ifx/libSpeedCHEM64.a` and `ifx/mod/*.mod`.
The Makefile lists sources in module dependency order and serializes
compilation even with `make -j`.

For an optimized build with debug symbols and a saved compilation log:

    ./scripts/ifx.sh

The script rebuilds through the Makefile using `-O3`, `-xHost`, and
big-endian binary I/O. The library is optimized for the build host's CPU;
regenerate `cklink` when switching from a build with different byte order.
The log is saved to `ifx/compile_ifx.log` (`LOG_FILE` can override it).

## Source layout

| File | Responsibility |
| --- | --- |
| `SCmodule.f90` | Mechanism state, thermodynamics, kinetics, and solver workspace |
| `chemkin_module.f90`, `SCcklink.f90` | CHEMKIN parsing and binary mechanism loading |
| `SCconV.f90` | Constant-volume equations and analytical Jacobian |
| `chemistry_input.f90` | Initialization, solver selection, and integration |
| `SCallocate.f90` | Allocation and mechanism finalization |
| `SCsparse*.f90`, `sparse_MPI.f90` | Sparse storage, algebra, and broadcast |
| `SCbroadcast.f90` | Mechanism broadcast to MPI workers |

Module names can differ from filenames. Keep `working_precision.f90`
first in the compilation order; `radau5s.f90` is excluded from the build.
Preserve source encodings and author notices when editing legacy files.
The current `cklink` format uses magic `SCLKv2  ` and schema 3; changes
to its binary layout must update both `CK_SCHEMA` and `CK_SCHEMA_EXPECT`.

## References and license

Original project: https://www.federicoperini.info/speedchem

Citation: [CITATION.cff](CITATION.cff)
SpeedCHEM-authored code: [GPL-3.0-or-later](LICENSE).
Bundled solver author credits, license notices, and provenance notes:
[NOTICE](NOTICE).
