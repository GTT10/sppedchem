# Invalid PLOG A test

The PLOG pre-exponential factor is zero. SpeedCHEM must reject it instead
of taking `log(0)` in the rate evaluator.
