# SpeedCHEM

[![CI](https://github.com/GTT10/sppedchem/actions/workflows/ci.yml/badge.svg)](https://github.com/GTT10/sppedchem/actions/workflows/ci.yml)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue.svg)](LICENSE)

SpeedCHEM is a Fortran library for integrating stiff gas-phase chemical-kinetics systems with sparse analytical Jacobians. It can be linked into multidimensional combustion codes or used to build stand-alone constant-volume reactor applications.

This repository is a **community-maintained modernization of Federico Perini's SpeedCHEM distribution**. It is not the official upstream site and is not presented as an official release by the original author. The original project page remains the authoritative source for the historical SpeedCHEM release and documentation.

## Current scope

The library currently provides:

- CHEMKIN-style ASCII mechanism and NASA/JANAF thermodynamic-data parsing;
- the SpeedCHEM sparse thermochemistry and analytical-Jacobian implementation;
- constant-volume kinetics integration through several stiff ODE/DAE solvers;
- MPI broadcast of parsed mechanisms and sparse structures;
- CHEMKIN `PLOG / P A b E /` support, including grouped same-pressure terms, pressure interpolation, analytical Jacobian coupling, and cklink serialization;
- safe mechanism finalization and reload in a single process;
- regression coverage for parsing, rates, Jacobians, integration, MPI, OpenMP, and long CHEMKIN identifiers.

Unsupported reaction forms are rejected rather than silently approximated. See `docs/testing.md` and the negative-test fixtures for the current compatibility boundary.

## Requirements

The maintained build currently targets Linux with:

- Intel oneAPI Fortran compiler (`ifx`), invoked through `mpiifx`;
- Intel MPI development tools;
- GNU Make.

Python 3 and Cantera 3.2 are required only for the external cross-implementation validation target.

## Build

Initialize the Intel oneAPI environment, then build the static library:

```bash
source /opt/intel/oneapi/setvars.sh
make
```

The build produces:

```text
ifx/libSpeedCHEM64.a
ifx/mod/*.mod
ifx/build/*.o
```

The source compilation order is intentionally serialized because the legacy modules do not yet have an automatically generated dependency graph. The top-level `Makefile` is the canonical build entry point.

## Test

Run all self-contained regression tests:

```bash
make test-all
```

Run the pinned public-mechanism comparison against Cantera:

```bash
python -m pip install "cantera==3.2.0"
make test-real-plog
```

The external test downloads a fixed C3Mech revision, verifies the downloaded Git blobs, converts the exact CHEMKIN/thermodynamic pair with Cantera, and compares PLOG rate coefficients, rates of progress, finite constant-volume evaluations, numerical-versus-analytical solver histories, and two-rank MPI results.

## Using the library

SpeedCHEM is a library rather than a complete end-user executable. A host application normally:

1. sets the mechanism directory through `chemistry_setup::mechdir`;
2. calls `chemistry_input` to parse and initialize the mechanism and solver state;
3. calls `chemistry_ODE_integrate` for each constant-volume integration;
4. calls `chemistry_finalize` before loading a different mechanism in the same process.

The test drivers under `test/` provide compact integration examples. `CLAUDE.md` contains the current module map, build-order constraints, and implementation notes for repository maintainers.

## Original project and attribution

SpeedCHEM was originally developed by Federico Perini. The original source headers identify:

```text
Copyright (C) 2010-2013 Federico Perini
```

Original project page: <https://www.federicoperini.info/speedchem>

This maintained tree preserves the original notices and records later modifications in Git history. Do not remove file-level author, copyright, provenance, or license notices.

## Citation

The original documentation requests citation of the SpeedCHEM numerical-method papers. The main reference is:

> F. Perini, E. Galligani, and R. D. Reitz, “An Analytical Jacobian Approach to Sparse Reaction Kinetics for Computationally Efficient Combustion Modeling with Large Reaction Mechanisms,” *Energy & Fuels*, 26(8), 4804–4822, 2012. DOI: 10.1021/ef300747n.

For direct/Krylov sparse-solver work, also cite:

> F. Perini, E. Galligani, and R. D. Reitz, “A Study of Direct and Krylov Iterative Sparse Solver Techniques to Approach Linear Scaling of the Integration of Chemical Kinetics with Detailed Combustion Mechanisms,” *Combustion and Flame*, 161(5), 1180–1195, 2014. DOI: 10.1016/j.combustflame.2013.11.017.

Machine-readable citation metadata is provided in `CITATION.cff`.

## License and third-party code

The SpeedCHEM-authored portions are distributed under the **GNU General Public License, version 3 or any later version** (`GPL-3.0-or-later`), consistent with the original source-file notices. See `LICENSE`.

This repository also contains numerical solver and legacy interface code with separate authorship and provenance. The root license does not erase file-level notices or independently relicense third-party material. See `THIRD_PARTY_NOTICES.md` before redistributing the source or a linked binary.

Because GPL-3.0-or-later is a strong copyleft license, distribution of linked applications can create source-code and licensing obligations. Consult the license text and obtain appropriate legal review for commercial or externally distributed products.

## Contributing

See `CONTRIBUTING.md`. Changes must preserve provenance notices, document newly introduced third-party material, and pass the relevant regression targets before review.
