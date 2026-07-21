# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

SpeedCHEM is a Fortran library for chemical-kinetics ODE/DAE problems (combustion chemistry). The build produces a static archive `libSpeedCHEM64.a` — there is no executable or main program here; this repo is a *library* that a host application links against. Original author: Federico Perini (GPLv3), circa 2010–2013; this tree adds a modern Makefile, `fpm.toml`, and dual gfortran/Intel-ifx build support.

## Building

The library is built two ways depending on the compiler. Output is segregated per compiler under `gfortran/` or `ifx/` (`<tag>/libSpeedCHEM64.a`, `<tag>/mod/*.mod`, `<tag>/build/*.o`).

```bash
make                       # gfortran -> gfortran/libSpeedCHEM64.a
make FC=mpif90             # MPI gfortran wrapper (still tagged 'gfortran')
make FC=mpiifx             # Intel ifx via MPI wrapper -> ifx/libSpeedCHEM64.a
make clean                 # remove only the current COMPILER_TAG output dir
make distclean             # remove gfortran/ ifx/ build/
```

`COMPILER_TAG` is auto-derived from `FC` (anything containing `ifx`/`ifort` → `ifx`, else `gfortran`) and can be overridden: `make COMPILER_TAG=ifx FC=mpiifx`.

Alternative standalone build scripts (do the same compile without make): `scripts/compile_gfort64opt.sh` (gfortran) and `scripts/ifx.sh` (mpiifx). These use `-O3 -g` with debug/backtrace flags; the Makefile uses `-O2`.

An `fpm.toml` exists but the Makefile/scripts are the canonical build. `fpm` is risky here: it auto-globs `src/` and would try to compile `radau5s.f`, which is deliberately excluded from the real build (see below).

### Critical: compile order is fixed and must not change

The source files have inter-module dependencies but there is **no dependency-graph build** — modules are compiled in one hard-coded serial order. This exact order is duplicated in three places and must be kept in sync: the `SRCS` list in [Makefile](Makefile), `scripts/compile_gfort64opt.sh`, and `scripts/ifx.sh`. Key constraints:

- `working_precision.f90` **must** be first (every other module `use`s it).
- The Makefile builds strictly serially even under `-j` (it chains each `.o` as a prerequisite of the next); do not "parallelize" it by removing that chain.

### Compiler flags that matter

- **`-fconvert=big-endian` (gfortran) / `-convert big_endian` (ifx)**: the library reads/writes big-endian binary data. Changing this breaks I/O compatibility with data files.
- **`-fallow-argument-mismatch`** (gfortran) is required — legacy solver code (LSODE/VODE families, etc.) passes mismatched argument types; without it the build fails. ifx generally does not need an equivalent.
- Fixed-form files (`.f`) use extended line length (`-ffixed-line-length-none` / `-extend-source 132`).
- `.f` = fixed-form legacy solvers/utilities; `.f90` = free-form (mostly the newer SpeedCHEM-authored modules).

### Encoding gotcha

Several `.f`/`.f90` files are **ISO-8859 / extended-ASCII, not UTF-8** (e.g. `SCconV.f90`, `SCutilities.f`). Plain `grep -E` may silently fail to match inside them — use `grep -a` when searching source. Do not "fix" the encoding blindly; the compilers accept these files.

## Architecture

### Module-name ≠ file-name

The most important thing to know: **module names do not match file names**, and several files each contain *multiple* modules. To find a module, search rather than guess the file. Map of the SpeedCHEM-authored modules:

| Module (`use ...`)        | File                          |
|---------------------------|-------------------------------|
| `working_precision`       | `working_precision.f90`       |
| `utilities`               | `SCutilities.f`               |
| `sparse_definitions`      | `SCsparse_definitions.f90`    |
| `sparse_algebra`          | `SCsparse.f`                  |
| `sparse_MPI`              | `sparse_MPI.f`                |
| `chemistry_setup`, `sparse_chemistry`, `speedchem`, `find_mod`, `universal_constants`, `SCthermodata`, `SCspeciesthermo`, `SCmixturethermo`, `troepar`, `reacpar`, `kinetics_mod`, `ode_solver` | **all in `SCmodule.f`** |
| `chemkinII`, `chemkinII_interpreter`, `chemkin_kiva` | **all in `chemkin_module.f90`** |
| `speedchem_conV`          | `SCconV.f90`                  |
| `radau_sparse`            | `radau_sparse.f90`            |
| `DVODE_F90_M`             | `dvode_f90_m.f90`             |

`SCmodule.f` is the core: it holds the mechanism data structures, thermodynamics, reaction-rate kinetics, and the `speedchem` public API plus `ode_solver` workspace state — roughly a dozen tightly-coupled modules in one ~186 KB file.

### Layers, roughly bottom-up

1. **Precision** — `working_precision` defines `dp` (real kind 8) used everywhere.
2. **Sparse infrastructure** — `sparse_definitions` (types), `sparse_algebra` (ops), `sparse_MPI` (MPI broadcast of sparse structures). The library's whole point is a **sparse analytical Jacobian** for stiff chemistry, so sparse types are pervasive.
3. **Chemistry core** (`SCmodule.f`) — thermo (NASA polynomials), reaction kinetics, the mechanism/setup state, and `speedchem` public routines.
4. **CHEMKIN I/O** (`chemkin_module.f90`) — parses CHEMKIN-format mechanism + thermo files (`chemkinII_interpreter`) into internal structures; `SCcklink.f90` links a CHEMKIN mechanism into SpeedCHEM state.
5. **RHS / Jacobian callbacks** (`SCconV.f90`, `speedchem_conV`) — the constant-volume derivative and Jacobian functions that the ODE solvers call back into.
6. **Solver drivers** (the many `.f` files: `vode.f`, `opkd*.f` = ODEPACK/LSODE family, `ddaspk.f`, `radau5.f`/`radaua.f`/`radau_sparse.f90`, `rodas.f`, `rowmap.f`, `MEBDFSO.f`, `gam.f`/`gamsub.f`, plus `dvode_f90_m.f90`). These are largely third-party stiff ODE/DAE integrators, kept in legacy fixed form.
7. **Top-level entry** — `chemistry_input.f90`: `chemistry_input` (read problem), `ODE_solver_speedchem_init` (set up workspaces/tolerances), and `chemistry_ODE_integrate` (the main integration driver).

### Solver selection is a string dispatch

`chemistry_ODE_integrate` in `chemistry_input.f90` picks the integrator with `select case (solver)` on a character string: `"VODE"`, `"VODES"`, `"VODEJAC"`, `"VODESJAC"`, `"LSODE"/"LSODES"/"LSODA"` (+`JAC` variants), `"RADAU5"`, `"RODAS"`, `"ROWMAP"`, etc. The `...JAC` / `...S` suffixes select analytical vs. numerical Jacobian and dense vs. sparse. When adding or wiring a solver, this `select case` is the integration point, and the corresponding driver file must be in the compile list.

### MPI

MPI is used only in `sparse_MPI.f` and `SCbroadcast.f90` (broadcasting the parsed mechanism / sparse structures to worker ranks). Building with a plain compiler (`gfortran`) still works; the MPI paths compile because the wrappers (`mpif90`/`mpiifx`) provide `mpif.h`. The gfortran script auto-detects the MPI include dir (`MPIHOME`, or via `mpif90 -showme`).

### Files present but not built

- `radau5s.f` — a sparse RADAU5 variant, **not** in any build list (superseded by `radau_sparse.f90`). Do not add it to the compile order.
- `gamparam.dat` — an `INCLUDE`d parameter file for `gamsub.f`, not compiled directly.
- `scripts/Makefile`, `scripts/compile_gfort*.bat`, `scripts/SC_gather_results.m` — legacy/aux (old in-tree Makefile fragment, Windows batch builds, a MATLAB results-gathering script). The top-level `Makefile` supersedes `scripts/Makefile`.
