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

## Change validation

CI builds the library with Intel oneAPI and checks that tracked files remain
clean. For numerical changes, record the mechanisms, conditions, and reference
results used for validation. Mechanism data and test drivers are not bundled
with this source distribution. Documentation-only changes need diff and link
checks.
