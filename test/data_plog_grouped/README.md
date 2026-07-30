# Positive-test mechanism: grouped same-pressure PLOG terms

The PLOG reaction has two Arrhenius entries at 1 atm. Standard
CHEMKIN/Cantera semantics sum those two rates at the pressure node before
logarithmic pressure interpolation.

`scripts/run_plog_tests.sh` verifies that both entries survive the
CHEMKIN -> cklink v2 -> SpeedCHEM round trip. The evaluator-level sum,
temperature derivative, and pressure slope are checked independently by
`test/driver_plog_eval.f90`.

Only `chem.inp` is stored here. The harness stages the shared test
`therm.dat` in a temporary directory.
