# Development

## Build

Use the root [Makefile](../Makefile) with Intel `mpiifx`. Its `SRCS` list is
the single source of compilation order; dependencies are serialized even with
`make -j`. Keep `working_precision.f90` first and exclude `radau5s.f90`.

Sources use free-form Fortran, with some legacy extended-ASCII text. Preserve
source encodings and third-party notices. Build flags and unformatted I/O byte
order affect binary compatibility; regenerate `cklink` after changing them.

## Code map

| File | Responsibility |
| --- | --- |
| `SCmodule.f90` | Mechanism state, thermo, kinetics, sparse chemistry, solver workspace |
| `chemkin_module.f90`, `SCcklink.f90` | CHEMKIN parsing and binary mechanism loading |
| `SCconV.f90` | Constant-volume RHS and analytical Jacobian |
| `chemistry_input.f90` | Initialization, solver dispatch, integration |
| `SCallocate.f90` | Allocation and mechanism finalization |
| `SCsparse*.f90`, `sparse_MPI.f90` | Sparse storage, algebra, broadcast |
| `SCbroadcast.f90` | Mechanism broadcast to MPI workers |

Module names can differ from filenames; several files contain multiple modules.
Search definitions before editing. Imported ODE/DAE solvers have separate
[provenance notices](../THIRD_PARTY_NOTICES.md).

## Supported input

CHEMKIN mechanism files and NASA-7 thermo records support 18-character species
names and 256-character reaction records. Thermo lookup accepts both a full
18-column identifier with an adjacent date and compact short-name headers.

PLOG requires positive A factors and non-decreasing pressures. Same-pressure
terms are summed before log-pressure interpolation; pressures outside the table
use the nearest endpoint. PLOG combined with explicit `REV`, third-body/falloff,
`FORD/RORD`, or non-integer stoichiometry is rejected.

Activation energies use a shared gas constant (`8.31446261815324 J/mol/K`)
and `4.184 J/cal`. Regenerate `cklink` after updating from older versions;
these corrections change rates relative to the old inconsistent conversions.

The analytical Jacobian includes pressure dependence through temperature and
composition. PLOG disables simplified sparsity. Call `chemistry_finalize`
before loading another mechanism in the same process.

The `cklink` magic is `SCLKv2  ` and its current schema is 3. Change both
`CK_SCHEMA` and `CK_SCHEMA_EXPECT` if the binary layout changes.

## Testing

| Target | Coverage |
| --- | --- |
| `make test-smoke` | Bundled PRF ignition |
| `make test-sparse` | Sparse row replacement and subsequent collider columns |
| `make test-constants` | Energy conversions in five units and explicit zero efficiencies |
| `make test-plog` | Parse/cklink, rates, RHS, Jacobian, integration, MPI, rejected inputs, capacity |
| `make test-lifecycle` | Mechanism reload and compact/full-width string handling |
| `make test` | Smoke, sparse, constants, PLOG, and lifecycle tests |
| `make test-openmp` | Isolated bounds-checking build and two-thread reload |
| `make test-all` | All bundled tests |
| `make test-real-plog` | Pinned public C3Mech comparison with Cantera 3.2 |

The real-PLOG test downloads C3MechV4.0.1 MID 2900 (H2/O2/Ar), verifies pinned
Git blobs, converts the exact CHEMKIN/thermo pair, and compares rates, production
rates, constant-volume RHS/history, and two-rank MPI results. Source revisions
and paths are in [the runner](../scripts/run_plog_real_mechanism.sh).

CI runs `test-all`, `test-real-plog`, and checks that tracked files remain clean.
For local changes, start with the affected tests; run the full suite before
merging code or build changes. Documentation-only changes need diff and link
checks. Test mechanism attribution and redistribution terms are in
[THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md).
