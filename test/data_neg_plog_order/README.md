# Negative-test mechanism: out-of-order PLOG pressures

PLOG pressure entries must be non-decreasing within a reaction. This
fixture decreases from 1.0 atm to 0.1 atm and must fail closed during
linking with an `out of order` diagnostic.
