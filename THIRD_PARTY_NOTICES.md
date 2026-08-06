# Third-party notices and provenance

This file records authorship and provenance visible in the source tree. It is not a legal opinion and does not replace the license text or notices embedded in individual files. File-level notices remain controlling for their respective material.

## SpeedCHEM core

The SpeedCHEM-authored source files identify:

```text
SpeedCHEM - A fast, portable Fortran library for Chemical Kinetics problems
Copyright (C) 2010-2013 Federico Perini
```

Those files state that they may be redistributed and modified under the GNU General Public License, version 3 or, at the recipient's option, any later version. The repository-level copy of that license is in `LICENSE`.

## RADAU5 and RODAS family

Relevant files include:

- `src/radau5.f90`
- `src/radaua.f90`
- `src/rodas.f90`

The source identifies Ernst Hairer and Gerhard Wanner as the authors of RADAU5/RODAS and cites *Solving Ordinary Differential Equations II*. The original code distribution is accompanied by the following notice:

```text
Copyright (c) 2004, Ernst Hairer

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
```

The repository preserves the in-source author and version notices. The notice above must also be retained when those files or resulting binaries are redistributed.

## ODEPACK and LLNL solver lineage

Relevant files include:

- `src/opkdmain.f90`
- `src/opkda1.f90`
- `src/opkda2.f90`
- `src/vode.f90`
- `src/ddaspk.f90`

The source headers identify Alan C. Hindmarsh and other Lawrence Livermore National Laboratory contributors. LLNL describes ODEPACK as public-domain software. `vode.f90` additionally records later SpeedCHEM-specific modifications by Federico Perini. `ddaspk.f90` records work performed under a U.S. Department of Energy contract.

Preserve all author, laboratory, contract, reference, and modification notices in these files. This notice does not attempt to assign a new license to material whose source file does not contain an explicit standalone grant.

## DVODE_F90 and sparse linear-algebra lineage

`src/dvode_f90_m.f90` identifies G. D. Byrne and S. Thompson as the DVODE_F90 maintainers and records contributions or source lineage from VODE, LSODAR, LSODES, LSOD28/MA28, and JACSP.

No standalone redistribution license for this complete combined file was found in this repository during the 2026 metadata review. Before redistributing this file or binaries that materially depend on it outside the research group, verify the applicable upstream terms, particularly for the embedded sparse linear-algebra lineage. The original provenance comments must remain intact.

## MEBDFSO

`src/MEBDFSO.f90` identifies T. J. Abdulla and J. R. Cash as authors and credits Torsten Hennig for a Fortran 90 version. No explicit standalone redistribution license was found in this repository during the 2026 metadata review. Preserve the source notice and verify upstream redistribution terms before external distribution.

## GAM

`src/gam.f90`, `src/gamsub.f90`, and `src/gamparam.dat` identify F. Iavernaro and F. Mazzia as authors. The source also states that some routines, comments, and implementation techniques were imported from RADAU5. No explicit standalone redistribution license for the complete GAM source was found in this repository during the 2026 metadata review. Preserve the source notices and verify upstream terms before external distribution.

## ROWMAP

`src/rowmap.f90` identifies H. Podhaisky, R. Weiner, and B. A. Schmitt as authors. No explicit standalone redistribution license was found in this repository during the 2026 metadata review. Preserve the source notice and verify upstream terms before external distribution.

## CHEMKIN/KIVA-derived interface

`src/chemkin_module.f90` identifies the code as a CHEMKIN-II interpreter/runtime taken from a KIVA4 chemistry version and modified by Federico Perini. The repository does not currently contain a separate upstream license or provenance package for that imported baseline. Preserve the existing header and verify the original KIVA/CHEMKIN redistribution terms before external distribution.

## ICSPLODE / sparse RADAU implementation

`src/radau_sparse.f90` identifies Emanuele Galligani and Federico Perini as authors. It is maintained as part of this SpeedCHEM tree; its original authorship notice must remain intact.

## Test-only external mechanism

`scripts/run_plog_real_mechanism.sh` downloads a pinned C3Mech mechanism and thermodynamic data at test time. Those files are not vendored into this repository. Their own upstream notices and terms apply to the downloaded material.

## Open provenance items

The repository deliberately records unresolved provenance instead of assuming that the root GPL notice relicenses every imported solver. See the project issue tracker for the current audit checklist. Until each item is resolved, internal research use and source preservation are lower-risk than redistributing a compiled library or linked application to third parties.
