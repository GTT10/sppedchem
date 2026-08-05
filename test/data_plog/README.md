# PLOG end-to-end test mechanism (cklink v2)

A minimal, **illustrative** (not physically calibrated) mechanism that
contains one `PLOG` (pressure-dependent Arrhenius) reaction, used for
binary plumbing, rate, analytic Jacobian, integration, and MPI tests:

```
chem.inp + therm.dat --[CKINTP]--> cklink(v2) --[SCcklink]--> reacpar PLOG arrays
                     --[plog_dump_canonical]--> canonical text
```

Species and thermo are reused from the bundled PRF `therm.dat` (copied
here), so no new thermodynamic data is needed.

## What it checks

- `PLOG / P A b E /` lines are recognised by the CKINTP interpreter
  (`CKAUXL`) and collected in the `plog_collect` module.
- Pressure is converted atm → `ln(P[Pa])`; activation energy goes through
  the **same** `EFAC` unit conversion as every other Arrhenius `E`
  (here cal/mole → `E/R [K]`), so PLOG cannot drift from the main rate
  units (cf. the KJOU 4× bug fixed in stage 0).
- cklink **v2** (magic `SCLKv2` + schema 3) writes the packed PLOG
  section and the full 18-character CHEMKIN/NASA-7 species field;
  `SCcklink` reads it back into the `reacpar` arrays and verifies the
  section checksum.
- The PLOG reaction is the **3rd** reaction (not reaction 1), so the
  reaction-index map in the packed arrays is genuinely exercised.
- Adjacent entries at the same pressure are legal grouped terms. They
  survive the round trip and are summed in log space before pressure
  interpolation.

## Files

- `chem.inp` — 2 plain Arrhenius reactions, one 3-node PLOG reaction, one
  more plain Arrhenius.
- `therm.dat` — copy of the PRF thermo database.
- `plog_expected.txt` — **golden** canonical dump. `scripts/run_plog_tests.sh`
  regenerates the dump and diffs it against this file.

Generated outputs (`cklink`, `chem.bin`, `SpeedCHEM.*`, `dat.*`,
`chem.out`) are git-ignored.

## How it is driven

- `make test-plog`: reproducible entry point for the bundled PLOG plumbing,
  rate/Jacobian/integration/MPI, and negative tests.
- `scripts/run_plog_tests.sh`: builds `test/driver_plog.f90` and asserts
  the canonical dump equals `plog_expected.txt`; it also checks grouped
  same-pressure terms and `REACTIONS MOLECULES` conversion.
- `scripts/run_plog_eval_tests.sh`: checks ordinary and extreme-range rate
  values and derivatives, a normalized full-RHS equivalence case, the
  analytic Jacobian, default-`LSODESJAC` integration, and two-rank MPI
  broadcast.
- `scripts/run_neg_tests.sh`: rejects malformed PLOG input and unsupported
  PLOG combinations such as explicit `REV` and third-body/falloff syntax.
  Duplicate pressure entries are not rejected; they are standard grouped
  PLOG terms and are covered by the positive fixture in
  `test/data_plog_grouped/`.
