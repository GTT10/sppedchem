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

See [supported input](docs/development.md#supported-input)
for supported reaction forms and limitations.

## Usage

Set the mechanism directory, initialize the library, and call
the integration routine from your application:

1. Set `chemistry_setup::mechdir`.
2. Call `chemistry_input`.
3. Call `chemistry_ODE_integrate`.
4. Call `chemistry_finalize` when finished.

See [example drivers](test/) for complete calling examples.

## About this version

This repository maintains Federico Perini's SpeedCHEM,
with added PLOG support and improvements to mechanism loading
and validation.

It is not an official upstream release.

## Building

The supplied Makefile uses Intel oneAPI and Intel MPI:

    make

See the [development guide](docs/development.md)
for build and test details.

## References and license

Original project: https://www.federicoperini.info/speedchem

Citation: [CITATION.cff](CITATION.cff)
License: [LICENSE](LICENSE)
Third-party notices: [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)
