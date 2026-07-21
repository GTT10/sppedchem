# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

SpeedCHEM is a Fortran library for chemical-kinetics ODE/DAE problems (combustion chemistry). The build produces a static archive `libSpeedCHEM64.a` — there is no executable or main program here; this repo is a *library* that a host application links against. Original author: Federico Perini (GPLv3), circa 2010–2013. **This tree targets Intel ifx only**, built via the MPI wrapper `mpiifx`.

## Building

Output goes under `ifx/` (`ifx/libSpeedCHEM64.a`, `ifx/mod/*.mod`, `ifx/build/*.o`).

```bash
make                       # -> ifx/libSpeedCHEM64.a  (FC defaults to mpiifx)
make FC=mpiifx             # explicit (same as default)
make clean                 # remove the ifx/ output dir
```

Alternative standalone build script (same compile without make): `scripts/ifx.sh`. It uses `-O3 -g` with debug/backtrace flags; the Makefile uses `-O2`.

An `fpm.toml` exists but the Makefile/scripts are the canonical build. `fpm` is risky here: it auto-globs `src/` and would try to compile `radau5s.f90`, which is deliberately excluded from the real build (see below).

## Testing

There is one end-to-end smoke test that drives the whole pipeline and gates on the result. Since the repo is a library with no main program, the test provides its own driver.

```bash
scripts/run_tests.sh                             # build + smoke test (test/data/)
SC_MECHDIR=/path/to/data/ scripts/run_tests.sh   # run against a different mechanism dir
```

`run_tests.sh` builds `libSpeedCHEM64.a`, compiles+links [test/driver_smoke.f90](test/driver_smoke.f90) via the MPI wrapper `mpiifx` (the library uses MPI in `sparse_MPI.f90`/`SCbroadcast.f90`, though the driver itself never calls `MPI_INIT`), runs it, and exits non-zero on failure.

What the driver does: sets `mechdir`, calls `chemistry_input` (the full setup orchestrator: CKINTP → cklink → SCcklink → sparse setup → `ODE_solver_speedchem_init`), seeds a **hot stoichiometric n-hexadecane/air state** (looked up by species name via the public `specie`/`ns` from `speedchem`), integrates constant-volume with `chemistry_ODE_integrate`, and asserts the solution is finite **and the temperature rose >100 K (autoignition)**. Baseline: the bundled mechanism links to 51 species / 153 reactions and ignites 1400 K → ~1798 K.

Test mechanism lives in [test/data/](test/data/) (`chem.inp` = a 51-species PRF n-hexadecane/iso-cetane skeletal mechanism, `therm.dat`). Only those two inputs are tracked; the run writes `cklink`, `SpeedCHEM.out`, `chem.out`, `dat.*` etc. into that dir and they are gitignored. Use this test as the regression gate for any refactor.

### Negative tests (fail-closed guards)

`scripts/run_neg_tests.sh` is the complement to `run_tests.sh`: it checks that mechanisms using reaction features SpeedCHEM **cannot** evaluate are *refused* (fail-closed), not silently mis-integrated. SpeedCHEM's most dangerous failure mode is "the interpreter accepts a keyword, the mechanism links and integrates, but the rate/Jacobian code ignores the feature" — so `SCcklink` counts each unsupported feature from the cklink header and, when present, prints a rejection naming the keyword + reaction count and `error stop 1`s. A negative case passes when the driver exits **non-zero** AND emits the expected message. (Note: plain `stop` is exit 0 in Fortran and would read as success; fail-closed guards use `error stop`.) Bundled cases: [test/data_neg_ford/](test/data_neg_ford/) — the PRF mechanism with one added `FORD` line; [test/data_plog/](test/data_plog/) — a PLOG mechanism, refused at integration because PLOG rate evaluation isn't implemented yet (stage 1 parses/stores but does not integrate PLOG); [test/data_neg_plog_dup/](test/data_neg_plog_dup/) — a PLOG reaction with a duplicate pressure, rejected by CKINTP under the `strict_chemkin` dialect. See each dir's README. Add a case by dropping a mechanism dir under `test/` and a `"dir|expected-substring"` line in the harness's `neg_cases`.

### PLOG data plumbing + cklink v2 (stage 1)

PLOG (pressure-dependent Arrhenius, `PLOG / P A b E /`) is parsed, stored, and round-tripped (stage 1) **and its forward rates are now evaluated** (stage 2 — see the next subsection). The path: `CKAUXL` recognises the keyword and calls into a dedicated **`plog_collect`** module (in `chemkin_module.f90`, before `chemkinII_interpreter`) that accumulates PLOG lines dynamically — no fixed `MAXPLOG` array. `CPREAC` applies the same `EFAC` activation-energy conversion to PLOG `E` as to `PAR(3,*)` (so PLOG units can't drift). `plog_finalize` packs them into per-reaction CSR arrays (`plog_reaction`/`plog_node_ptr`/`plog_logP`/`plog_A`/`plog_b`/`plog_EoverR`, defined in the `reacpar` module in `SCmodule.f90`) and enforces the `strict_chemkin` dialect (pressures strictly ascending per reaction; duplicate/out-of-order → `error stop 1`).

The linking file is now **cklink v2**: a leading `WRITE (LINC) CK_MAGIC, CK_SCHEMA` record (magic `'SCLKv2  '`, schema `2`) precedes the legacy `VERS,PREC,KERR` header, and a PLOG section (counts → CSR arrays → integer checksum) is appended after the reaction data. `SCcklink` verifies the magic/schema (fail-closed on mismatch — this replaced the stage-0 `VERS`-string warning), then reads the PLOG section into the `reacpar` arrays and sets a per-reaction `rate_form` tag (`RATE_ARRHENIUS`/`RATE_PLOG`/…). **A mechanism with no PLOG writes a count-0 section and is byte-for-byte numerically identical downstream — the PRF baseline stays 1797.681 K.** Bumping the on-disk layout means bumping `CK_SCHEMA` in `chemkin_module.f90` and `CK_SCHEMA_EXPECT` in `SCcklink.f90` together.

Because SCcklink rejects every unsupported optional section (`FORD`/`JAN`/…) with `error stop` *before* the PLOG read, the reader can rely on the file being positioned exactly at the PLOG section after the REV/FAL/THB reads — so the trailing PLOG section is read by a plain forward read, no rewind/skip bookkeeping.

`scripts/run_plog_tests.sh` is the **stage-1 plumbing** harness: it builds the parse-only driver [test/driver_plog.f90](test/driver_plog.f90) (which calls only `ckintp` + `SCcklink`, never integrates), runs it with `SC_PLOG_DUMP=1` (an env flag, checked by `reacpar::plog_dump_only`, that lets the parse-only path run), and diffs the canonical PLOG dump against the golden `test/data_plog/plog_expected.txt`. Regenerate the golden only when the packed representation legitimately changes.

### PLOG rate evaluation (stage 2)

PLOG **forward rate constants are evaluated** by `reacpar::plog_kinf_eval(Ta, P_pa, kinf)`, called from `mass_action` (in `SCmodule.f90`) right after the Arrhenius `kinf` is formed and before `kf = kinf`. For each PLOG reaction it overwrites `kinf(r)` with the pressure-interpolated value: per node `k_i = A_i·exp(b_i·lnT − (E/R)_i/T)` (units match `Ainf`; `E/R` is in K so no `Rcal`), then **log-linear interpolation of ln k in ln P** between the two bracketing pressure nodes, with **nearest-endpoint clamping** outside the node range (`s=0`). Pressure comes from `SCP` [Pa] (the current state pressure; PLOG nodes store `ln(P[Pa])`). A mechanism with no PLOG reactions never calls it, so non-PLOG numerics are unchanged (PRF stays 1797.681 K).

**Analytic Jacobian is stage 3, so it is not yet PLOG-aware.** `chemistry_ODE_integrate` therefore **refuses any `...JAC` solver when PLOG is present** (`index(solver,"JAC")>0` → `error stop 1`) — integrating with the analytic Jacobian (`constV_jac_sparse`, unchanged) would use a wrong Jacobian. Numeric-Jacobian solvers finite-difference the RHS (which includes PLOG) and are fine in principle. **Caveat:** the stiff numeric solvers (LSODES/VODE without `JAC`) fail to converge in this build *even on plain PRF* — a pre-existing limitation unrelated to PLOG — so the default working solver is `LSODESJAC`.

`scripts/run_plog_eval_tests.sh` is the **stage-2** harness, both checks independent of the flaky stiff solvers: (1) a unit test [test/driver_plog_eval.f90](test/driver_plog_eval.f90) that compares `plog_kinf_eval` against an independent closed-form reference at nodes / geometric-mean pressures / out-of-range / multiple T (machine-precision match, cross-checked against Python); (2) an **RHS-equivalence** check [test/driver_rhs_probe.f90](test/driver_rhs_probe.f90) that evaluates the constant-volume RHS once for plain PRF and for [test/data_plog_prf/](test/data_plog_prf/) (PRF with `O+H2=H+OH` PLOG-ified using nodes that reproduce its Arrhenius rate) and asserts the two derivative dumps are **byte-identical** — proving PLOG feeds the real RHS correctly. The smoke driver gained an `SC_SOLVER` env override for choosing a non-default solver.

### Critical: compile order is fixed and must not change

The source files have inter-module dependencies but there is **no dependency-graph build** — modules are compiled in one hard-coded serial order. This exact order is duplicated in two places and must be kept in sync: the `SRCS` list in [Makefile](Makefile) and `scripts/ifx.sh`. Key constraints:

- `working_precision.f90` **must** be first (every other module `use`s it).
- The Makefile builds strictly serially even under `-j` (it chains each `.o` as a prerequisite of the next); do not "parallelize" it by removing that chain.

### Source form

All sources are now **free-form `.f90`** (the original fixed-form `.f` files were machine-converted with `findent`; see the "convert all fixed-form sources" commit). `-extend-source 132` remains in the flags to allow long free-form lines. ifx tolerates the legacy argument-type mismatches in the LSODE/VODE solver families without a special flag.

Note when touching the ex-fixed-form solvers: a few relied on fixed-form quirks that free form makes fatal — blanks inside tokens (`1 0` meaning `10`, `2.0  D0` meaning `2.0D0`) and tab-indented lines. These were fixed during conversion, but keep it in mind if you ever re-import upstream fixed-form code.

### Compiler flags that matter

- **`-convert big_endian`** (in `scripts/ifx.sh`): the library reads/writes big-endian binary data. Changing this breaks I/O compatibility with data files. (The Makefile currently omits this.)

### Encoding gotcha

Several sources are **ISO-8859 / extended-ASCII, not UTF-8** (e.g. `SCconV.f90`, `SCutilities.f90`). Plain `grep -E` may silently fail to match inside them — use `grep -a` when searching source. Do not "fix" the encoding blindly; ifx accepts these files.

## Architecture

### Module-name ≠ file-name

The most important thing to know: **module names do not match file names**, and several files each contain *multiple* modules. To find a module, search rather than guess the file. Map of the SpeedCHEM-authored modules:

| Module (`use ...`)        | File                          |
|---------------------------|-------------------------------|
| `working_precision`       | `working_precision.f90`       |
| `utilities`               | `SCutilities.f90`               |
| `sparse_definitions`      | `SCsparse_definitions.f90`    |
| `sparse_algebra`          | `SCsparse.f90`                  |
| `sparse_MPI`              | `sparse_MPI.f90`                |
| `chemistry_setup`, `sparse_chemistry`, `speedchem`, `find_mod`, `universal_constants`, `SCthermodata`, `SCspeciesthermo`, `SCmixturethermo`, `troepar`, `reacpar`, `kinetics_mod`, `ode_solver` | **all in `SCmodule.f90`** |
| `chemkinII`, `chemkinII_interpreter`, `chemkin_kiva` | **all in `chemkin_module.f90`** |
| `speedchem_conV`          | `SCconV.f90`                  |
| `radau_sparse`            | `radau_sparse.f90`            |
| `DVODE_F90_M`             | `dvode_f90_m.f90`             |

`SCmodule.f90` is the core: it holds the mechanism data structures, thermodynamics, reaction-rate kinetics, and the `speedchem` public API plus `ode_solver` workspace state — roughly a dozen tightly-coupled modules in one ~186 KB file.

### Layers, roughly bottom-up

1. **Precision** — `working_precision` defines `dp` (real kind 8) used everywhere.
2. **Sparse infrastructure** — `sparse_definitions` (types), `sparse_algebra` (ops), `sparse_MPI` (MPI broadcast of sparse structures). The library's whole point is a **sparse analytical Jacobian** for stiff chemistry, so sparse types are pervasive.
3. **Chemistry core** (`SCmodule.f90`) — thermo (NASA polynomials), reaction kinetics, the mechanism/setup state, and `speedchem` public routines.
4. **CHEMKIN I/O** (`chemkin_module.f90`) — parses CHEMKIN-format mechanism + thermo files (`chemkinII_interpreter`) into internal structures; `SCcklink.f90` links a CHEMKIN mechanism into SpeedCHEM state.
5. **RHS / Jacobian callbacks** (`SCconV.f90`, `speedchem_conV`) — the constant-volume derivative and Jacobian functions that the ODE solvers call back into.
6. **Solver drivers** (`vode.f90`, `opkd*.f90` = ODEPACK/LSODE family, `ddaspk.f90`, `radau5.f90`/`radaua.f90`/`radau_sparse.f90`, `rodas.f90`, `rowmap.f90`, `MEBDFSO.f90`, `gam.f90`/`gamsub.f90`, plus `dvode_f90_m.f90`). These are largely third-party stiff ODE/DAE integrators (originally fixed-form, now converted to free form).
7. **Top-level entry** — `chemistry_input.f90`: `chemistry_input` (read problem), `ODE_solver_speedchem_init` (set up workspaces/tolerances), and `chemistry_ODE_integrate` (the main integration driver).

### Solver selection is a string dispatch

`chemistry_ODE_integrate` in `chemistry_input.f90` picks the integrator with `select case (solver)` on a character string: `"VODE"`, `"VODES"`, `"VODEJAC"`, `"VODESJAC"`, `"LSODE"/"LSODES"/"LSODA"` (+`JAC` variants), `"RADAU5"`, `"RODAS"`, `"ROWMAP"`, etc. The `...JAC` / `...S` suffixes select analytical vs. numerical Jacobian and dense vs. sparse. When adding or wiring a solver, this `select case` is the integration point, and the corresponding driver file must be in the compile list.

### MPI

MPI is used only in `sparse_MPI.f90` and `SCbroadcast.f90` (broadcasting the parsed mechanism / sparse structures to worker ranks). Everything is built through the MPI wrapper `mpiifx`, which supplies `mpif.h`; the smoke-test driver never calls `MPI_INIT` and runs single-process.

### Files present but not built

- `radau5s.f90` — a sparse RADAU5 variant, **not** in any build list (superseded by `radau_sparse.f90`). Do not add it to the compile order.
- `gamparam.dat` — an `INCLUDE`d parameter file for `gamsub.f90`, not compiled directly.
- `scripts/Makefile`, `scripts/SC_gather_results.m` — legacy/aux (old in-tree Makefile fragment; a MATLAB results-gathering script). The top-level `Makefile` supersedes `scripts/Makefile`.
