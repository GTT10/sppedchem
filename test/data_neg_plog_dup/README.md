# Negative-test mechanism: duplicate PLOG pressure (strict_chemkin)

The PLOG plumbing mechanism (`test/data_plog/`) with the PLOG reaction
changed so **two PLOG lines share the same pressure** (1.0 atm):

```
H + O2 = HO2             1.000E+12   0.00     500.0
    PLOG /  1.0   1.000E+12   0.00     500.0 /
    PLOG /  1.0   2.000E+12   0.00     600.0 /   <-- duplicate pressure
```

## What it tests

SpeedCHEM's default PLOG dialect is `strict_chemkin`, which forbids two
PLOG lines at the same pressure within one reaction (the current Ansys
Chemkin rule; use `DUPLICATE` reactions instead). `plog_finalize`
(in the `plog_collect` module) enforces strictly-ascending pressures per
reaction and, on a violation, prints

```
ERROR...PLOG pressures for reaction N are not strictly ascending
        (duplicate or out-of-order pressure; strict_chemkin)
```

and **`error stop 1`**s. It uses `error stop`, not a plain `stop` — a
Fortran `stop` is exit 0 and would read as success to the harness; only
`error stop` makes the rejection detectable by exit code (the same
fail-closed lesson as stage 0).

`scripts/run_neg_tests.sh` drives this mechanism (via `driver_smoke`) and
asserts the run exits non-zero **and** emits the `not strictly ascending`
message. If the strict-dialect check regresses, the mechanism would link
(and, once PLOG rates exist, integrate) with a silently mis-ordered
pressure table.

Only `chem.inp` (+ `therm.dat`, copied from `test/data_plog/`) is
meaningful; generated outputs are git-ignored.
