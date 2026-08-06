# Bundled PRF smoke-test mechanism

This directory contains the mechanism used by the self-contained ignition smoke test. It is test data, not an original SpeedCHEM component and not a new mechanism released by this repository.

## Scientific provenance

The `chem.inp` header attributes the n-hexadecane/iso-cetane skeletal primary-reference-fuel mechanism to:

> W. Fan, M. Jia, Y. Chang, and M. Xie, “Understanding the Relationship between Cetane Number and the Ignition Delay in Shock Tubes for Different Fuels Based on a Skeletal Primary Reference Fuel (n-Hexadecane/Iso-cetane) Mechanism,” *Energy & Fuels*, 29(5), 3413–3427, 2015. DOI: 10.1021/ef5028185.

The publication identifies the CHEMKIN mechanism and thermodynamic properties as Supporting Information.

## Local adaptation

The tracked mechanism is not represented as a byte-identical copy of the authors' Supporting Information. It contains SpeedCHEM-oriented adaptations, including comments and pseudo-fragment handling used to satisfy participant limits. Git history records later test-maintenance changes.

Derived fixtures are located in:

- `test/data_plog_prf/`: one Arrhenius reaction converted to an equivalent PLOG representation;
- `test/data_neg_ford/`: one `FORD` line added to exercise fail-closed behavior.

Both derived directories copy `therm.dat` from this directory.

## Redistribution status

No explicit redistribution or modification license for these tracked mechanism and thermodynamic-data files was found in this repository during the 2026 provenance review. Free online availability as Supporting Information does not by itself establish redistribution permission.

Preserve the attribution above, do not describe these files as an original release by the paper authors, and consult `THIRD_PARTY_NOTICES.md` and issue #3 before redistributing the repository or these fixtures outside the research group.
