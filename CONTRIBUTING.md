# Contributing

## Scope

This repository maintains a research-oriented modernization of SpeedCHEM. Changes should preserve numerical behavior unless a pull request explicitly documents and validates an intended behavior change.

## Before submitting a change

1. Create a focused branch from `master`.
2. Preserve all existing copyright, license, author, contract, and provenance notices.
3. Do not add third-party source unless its redistribution terms are known and recorded in `THIRD_PARTY_NOTICES.md`.
4. Add or update regression coverage for behavior changes.
5. Run the narrowest relevant tests and `make test-all` before requesting merge.

For PLOG or kinetics changes, also run `make test-real-plog` with Cantera available.

## Licensing of contributions

Unless a file states otherwise, contributions to SpeedCHEM-maintained files are submitted under `GPL-3.0-or-later`, consistent with the original SpeedCHEM source notices. Contributors must have the right to submit their changes under those terms.

Do not add a GPL SPDX identifier to imported third-party files merely because they are compiled into this repository. Retain their original notices and document their provenance separately.

## Source conventions

- The canonical compiler is Intel `ifx` through `mpiifx`.
- Keep the source order in `Makefile` and `scripts/ifx.sh` synchronized.
- Do not add `src/radau5s.f90` to the canonical build.
- Keep mechanism parsers fail-closed for unsupported reaction forms.
- Avoid changing source encodings without a dedicated migration and regression plan.

## Pull-request description

A pull request should state:

- what changed and why;
- affected numerical or file-format behavior;
- tests run and their results;
- any new external source, data, license, or citation requirement;
- compatibility or migration implications.
