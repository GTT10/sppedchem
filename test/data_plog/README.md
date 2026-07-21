# PLOG plumbing test mechanism (cklink v2, stage 1)

A minimal, **illustrative** (not physically calibrated) mechanism that
contains one `PLOG` (pressure-dependent Arrhenius) reaction, used to test
the stage-1 PLOG data-plumbing path:

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
- cklink **v2** (magic `SCLKv2` + schema 2) writes the packed PLOG
  section; `SCcklink` reads it back into the `reacpar` CSR arrays and
  verifies the section checksum.
- The PLOG reaction is the **3rd** reaction (not reaction 1), so the
  reaction-index map in the packed arrays is genuinely exercised.

## Files

- `chem.inp` — 2 plain Arrhenius reactions, one 3-node PLOG reaction, one
  more plain Arrhenius.
- `therm.dat` — copy of the PRF thermo database.
- `plog_expected.txt` — **golden** canonical dump. `scripts/run_plog_tests.sh`
  regenerates the dump and diffs it against this file.

Generated outputs (`cklink`, `chem.bin`, `SpeedCHEM.*`, `dat.*`,
`chem.out`) are git-ignored.

## How it is driven

- `scripts/run_plog_tests.sh` (positive): builds `test/driver_plog.f90`,
  runs it with `SC_PLOG_DUMP=1` (which makes `SCcklink` skip its PLOG
  fail-closed guard), and asserts the dump equals `plog_expected.txt`.
- `scripts/run_neg_tests.sh` (negative): runs the **integrating**
  `driver_smoke` on this same mechanism *without* `SC_PLOG_DUMP`, and
  asserts it is refused (exit non-zero, message contains `PLOG`). Stage 1
  parses and stores PLOG but does not evaluate PLOG rates, so integrating
  must fail closed rather than silently ignore the pressure dependence.
  That guard is lifted in stage 2 when PLOG rates are implemented.
