# PLOG-on-PRF test mechanism (stage 2)

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

## What it tests

- **RHS equivalence** (`scripts/run_plog_eval_tests.sh`, via
  `test/driver_rhs_probe.f90`): the constant-volume derivative
  d(state)/dt for this mechanism is **byte-identical** to plain PRF's.
  That proves PLOG is wired into the real RHS and reproduces the
  Arrhenius rate it stands in for.

- **Analytic-Jacobian guard** (`scripts/run_neg_tests.sh`): run with the
  default solver (`LSODESJAC`, analytic Jacobian), `chemistry_ODE_integrate`
  **refuses** it — PLOG forward rates are evaluated (stage 2) but the
  analytic PLOG Jacobian is stage 3, so a `...JAC` solver would use an
  inconsistent Jacobian. The rejection names the solver and points to a
  numeric-Jacobian solver. To integrate it, set `SC_SOLVER` to a
  numeric-Jacobian solver (note: the stiff numeric solvers in this build
  are independently flaky even on plain PRF, so RHS equivalence — not an
  ignition run — is the stage-2 integration proof).

Only `chem.inp` (+ `therm.dat`, copied from `test/data/`) is meaningful;
generated outputs are git-ignored.
