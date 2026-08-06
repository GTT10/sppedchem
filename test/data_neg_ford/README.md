# Negative-test mechanism: modified species order (FORD)

This is the bundled 51-species PRF mechanism (`test/data/`) with **one line
added**: a `FORD` (forward reaction order) override on the first reaction:

```
IC16H34+O2 = IC16H33+HO2 6.000E+13 0.00 46000.0
 REV / 1.000E+12 0.00 0.0 /
 FORD / IC16H34 0.5 /          <-- added
```

## Provenance

This is a locally modified negative-test fixture derived from the Fan–Jia–Chang–Xie PRF mechanism described in `test/data/README.md`, not an original author release. Its publication attribution and unresolved redistribution status are recorded there and in `THIRD_PARTY_NOTICES.md`. `therm.dat` is copied unchanged from `test/data/`.

## What it tests

`FORD`/`RORD` (modified species order) is parsed by the CKINTP interpreter and
written into `cklink`, but SpeedCHEM's rate/Jacobian code does **not** evaluate
it (species orders are truncated to integers). Previously such a mechanism would
link and integrate while silently ignoring the override — the most dangerous
failure mode.

`SCcklink` now fails closed: it counts `NORD` from the (corrected) cklink header
and, when `nord > 0`, prints a rejection naming the feature and reaction count
and `error stop 1`s.

`scripts/run_neg_tests.sh` drives this mechanism and asserts the run exits
non-zero **and** emits the expected `modified species order (FORD/RORD)`
message. If someone regresses the fail-closed guard (or the header read), the
driver would either integrate to ignition (exit 0) or omit the message, and the
negative harness would fail.

Only `chem.inp` (+ `therm.dat`, copied from `test/data/`) is meaningful here;
generated run outputs (`cklink`, `SpeedCHEM.out`, `chem.out`, `dat.*`, …) are
git-ignored.
