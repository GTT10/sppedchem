# PLOG-on-PRF regression mechanism

The bundled 51-species PRF n-hexadecane/iso-cetane mechanism
(`test/data/`) with **one reaction converted to PLOG form**:

```
O+H2 = H+OH 0.508E+05 2.67 0.629E+04
    PLOG /  0.1   0.508E+05  2.67  0.629E+04 /
    PLOG /  1.0   0.508E+05  2.67  0.629E+04 /
    PLOG / 100.0  0.508E+05  2.67  0.629E+04 /
```

All three PLOG nodes carry the **same Arrhenius coefficients as the base
reaction line**, so the log-linear PLOG interpolation returns that exact
Arrhenius rate at *every* pressure. The mechanism is therefore
numerically identical to plain PRF — but it now exercises the PLOG
evaluation path (`reacpar::plog_kinf_eval`, called from `mass_action`).

Because it derives from the working PRF mechanism, it links and builds
the sparse chemistry cleanly (unlike a hand-built toy mechanism).

## Provenance

This is a locally modified regression fixture derived from the Fan–Jia–Chang–Xie PRF mechanism described in `test/data/README.md`, not an original author release. Its publication attribution and unresolved redistribution status are recorded there and in `THIRD_PARTY_NOTICES.md`. `therm.dat` is copied unchanged from `test/data/`.

## What it tests

- **RHS equivalence** (`scripts/run_plog_eval_tests.sh`, via
  `test/driver_rhs_probe.f90`): the constant-volume derivative
  d(state)/dt for this mechanism is **byte-identical** to plain PRF's.
  That proves PLOG is wired into the real RHS and reproduces the
  Arrhenius rate it stands in for.

- **Analytic-Jacobian regression** (`scripts/run_plog_eval_tests.sh`):
  the default `LSODESJAC` path is supported. A separate compact test
  compares every analytic Jacobian column with central differences.

Only `chem.inp` (+ `therm.dat`, copied from `test/data/`) is meaningful;
generated outputs are git-ignored.
