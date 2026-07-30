!   ***********************************************************************************************
!   **                                                                                           **
!   **         |\   -  -.   ./                                                                   **
!   **         | \./ \/ | ./ /     _____                     __________  __________  ___         **
!   **       __|         /  /     / ___/____  ___  ___  ____/ / ____/ / / / ____/  |/  /         **
!   **       \    .        /.-/   \__ \/ __ \/ _ \/ _ \/ __  / /   / /_/ / __/ / /|_/ /          **
!   **        \   |\.|\/|    /   ___/ / /_/ /  __/  __/ /_/ / /___/ __  / /___/ /  / /           **
!   **         \__\     /___/   /____/ .___/\___/\___/\__,_/\____/_/ /_/_____/_/  /_/            **
!   **                              /_/                                                          **
!   **                                                                                           **
!   **       SpeedCHEM - A fast, portable Fortran library for Chemical Kinetics problems         **
!   **                          http://www.federicoperini.info/speedchem                         **
!   **                 Copyright (C) 2010-2013 Federico Perini. Version 1.0                      **
!   **                                                                                           **
!   ***********************************************************************************************
!   **                                                                                           **
!   ** License                                                                                   **
!   ** This file is part of SpeedCHEM.                                                           **
!   **                                                                                           **
!   ** SpeedCHEM is free software: you can redistribute it and/or modify                         **
!   ** it under the terms of the GNU General Public License as                                   **
!   ** published by the Free Software Foundation, either version 3 of the                        **
!   ** License, or any later version.                                                            **
!   **                                                                                           **
!   ** SpeedCHEM is distributed in the hope that it will be useful,                              **
!   ** but WITHOUT ANY WARRANTY; without even the implied warranty of                            **
!   ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                             **
!   ** GNU General Public License for more details.                                              **
!   **                                                                                           **
!   ** You should have received a copy of the GNU General Public License                         **
!   ** along with SpeedCHEM. If not, see <http://www.gnu.org/licenses/>.                         **
!   **                                                                                           **
!   ** Please send comments, requests, bug reports to federico.perini@gmail.com                  **
!   **                                                                                           **
!   ***********************************************************************************************

!     *****************************************************************
!     **                                                             **
!     **                  SpeedCHEM - CHEMKIN link                   **
!     **                                                             **
!     **   Setup mechanism properties from CHEMKIN-II binary file    **
!     **             into the speedchem working arrays               **
!     **                                                             **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: thursday, 27/04/2012                         **
!     **                                                             **
!     *****************************************************************

      subroutine SCcklink
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

      use working_precision
!ck2015      use chemistry_setup, only: mechanism
      use chemistry_setup, only: mechanism, mechdir
      use speedchem,   only: nel, ns, nr, neq, uneq, species,         &
                             reactions, elMW, uelMW, elementi,        &
                             SCMW, uMW, specie, Arev, brev, Erev,     &
                             A0, b0, E0, Ainf, binf, Einf,            &
                             aTROE, T1TROE, T2TROE, T3TROE, TROE,     &
                             HIGH, LOW, reversibile, third_body,      &
                             inotrev, THREE, iTHREE, nTHREE,          &
                             nudiff, REV, third_body_beta, nnotrev,   &
                             isumnudiff, nTB, iTB
      use SCthermodata, only: aL, bL, cL, dL, eL, fL, gL, tsw,        &
                              aH, bH, cH, dH, eH, fH, gH, EM

      use troepar,      only: todotroe, aT2, uT1T2, T2T2, uT3T2,      &
                              zeroT2

      use reacpar,      only: Arrhreac, Lindreac, Troereac, Revreac,   &
                              rate_form, RATE_ARRHENIUS, RATE_PLOG,     &
                              n_plog_reactions, n_plog_nodes,           &
                              n_plog_terms,                             &
                              plog_reaction, plog_node_ptr, plog_logP,  &
                              plog_term_ptr, plog_A, plog_b, plog_EoverR

      use kinetics_mod, only: ip, jp, ir, jr, v_stoich_p, v_stoich_r, &
                              indice_p, indice_r, ijp, ijr,           &
                              iv_stoich_p, i2D1p, i2D2p, i1p, i2p,    &
                              iv_stoich_r, i2D1r, i2D2r, i1r, i2r

      use find_mod,     only: find_indices, find_indices2D, indices,  &
                              i2D1, i2D2
      use sparse_chemistry, only: stoich_r_sp, stoich_p_sp,           &
                                  nudiff_sparse, third_body_sp,       &
                                  tb_beta_sp

      use sparse_algebra
      use universal_constants, only: zero, one, R, joule_to_cal


!     Note: for consistency with the CHEMKIN the universal gas
!           constant in calorie units has been taken to be
!           1.987 cal / mol / K instead than the exact value
!      use universal_constants, only: Rcal

      implicit none

      real (dp)        :: Rcal!, parameter :: Rcal = 1.987d0


!     ** Variable declaration *****************************************


      integer,          parameter :: linc   = 51,                     &
                                     lout   = 50
!ck2015      character(len=*), parameter :: cklink = 'cklink',               &
!                                     cklin2 = 'chem.bin',             &
!                                     output = 'SpeedCHEM.out'
      character(len=255) :: cklink,cklin2,output

      character(len=*), parameter ::                                  &
         fmt_nofile = "(' No chemkin linking file found; ')",         &
         fmt_nofil2 = "(' Either cklink or chem.bin required for',    &
                        ' chemistry mechanism input. Exiting. ')",    &
         fmt_head   = "(72('-'))",                                    &
         fmt_head2  = "(23x,'CHEMKIN to SpeedCHEM link',23x)",        &
         fmt_fname  = "(' Reading mechanism input from: ',A20)",      &
         fmt_info   = "(' CKII version: ',A4,'; precision: ',A10)",   &
         fmt_nel    = "(' Number of elements : ',I5)",                &
         fmt_nr     = "(' Number of reactions: ',I5)",                &
         fmt_ns     = "(' Number of species  : ',I5)",                &
         fmt_lt     = "(' Landau-Teller reactions (',I5,')'"//        &
                       "' are not supported in this code. ')",        &
         fmt_chrg   = "(' Reactions with charge not supported in',"// &
                       "' this code. ')",                             &
         fmt_tmp    = "(' Thermo data needs 3 temperature bounds',"// &
                        "' and 7 JANAF coefficient;'/' found',"//     &
                       "' ntemp = ',I2,' and ncoefs = ',I2)",         &
         fmt_eldata = "(' Element   Atomic weight (AMU)',/,"//        &
                       "' -----------------------------')",           &
         fmt_elit   = "(' ',I2,'. ',A2,8X,F10.5)",                    &
         fmt_spec   = "(' Species and atomic weights [g/mol] ',/,"//  &
                       "' ',40('-'))",                                &
         fmt_spdet  = "(' ',I5,'. ',A18,2X,F9.5)",                    &
         fmt_thhead = "(' Thermodynamic database data for species',"//&
                         "' ',87('-'))",                              &
         fmt_thl1   = "(4(1X,1PE21.12))",                             &
         fmt_sri    = "(' ',I3,' SRI reactions not supported.')",     &
         fmt_erfl   = "(' Unrecognized falloff reaction type: ',I3)", &
         fmt_prec   = "(' Unsupported real data precision ')",        &
         fmt_unsup  = "(' ERROR: ',I5,' reaction(s) use ',A,', which',"//&
                       "' SpeedCHEM does not evaluate. The CKINTP',"//    &
                       "' interpreter accepts these keywords but the',"// &
                       "' rate/Jacobian code ignores them, which would',"//&
                       "' silently produce a wrong mechanism.',"//        &
                       "' Aborting.')"

!     ckfile receives adjustl(cklink)/adjustl(cklin2), which are len=255.
!     At len=80 long mechdir paths were silently truncated. Match them.
      character(len=255)                           :: ckfile
      character(len=16)                            :: ckvers, ckprec
      character(len=16), dimension(:), allocatable :: tmpchar

      logical :: ispresent, ispresent2, kerr, invalid_plog
      integer :: i, j, k, idummy1, idummy2
      integer :: liwork, lrwork, lcwork, mm, kk, ii, maxsp, maxtb,    &
                 maxtp, nthcf, nipar, nitar, nifar, nrv, nfl, nthb,   &
                 nlt, nrl, nw, nchrg, max_nreacs, max_nprods
!ck2015 for real stoichometric coefficients
!     nei=NEIM, nja=NJAR, njan=NJAN, nif1=NF1R, nf1=NFT1, nex=NEXC,
!     nsto=NRNU, nord=NORD, maxord=MAXORD  (see CKINTP header WRITE order)
      integer :: nei,nja, njan, nf1, nif1, nex, nsto, nord, maxord
      real(dp) :: ckmin
!ck2015 cklink v2 leading record: magic string + integer schema version.
!     Must match CKINTP's CK_MAGIC/CK_SCHEMA (chemkin_module.f90).
      character(len=8) :: ck_magic_in
      integer          :: ck_schema_in
      character(len=8), parameter :: CK_MAGIC_EXPECT  = 'SCLKv2  '
      integer,          parameter :: CK_SCHEMA_EXPECT = 2
!     PLOG section scratch + checksum
      integer :: plog_chk_in, plog_chk_calc, ipl, inode

      logical,          dimension(:)  , allocatable :: reac_cond

      integer,          dimension(:)  , allocatable :: idummy0, idummy
      integer,          dimension(:)  , allocatable :: ire, itype
      integer,          dimension(:),   allocatable :: inthb, imole
      integer,          dimension(:),   allocatable :: count_reacs
      integer,          dimension(:,:), allocatable :: itmpst,  itmpsp
!ck2015 for real stoichometric coefficients
      real (dp)       , dimension(:,:), allocatable :: rtmpst
      integer,          dimension(:),   allocatable :: tmpip, tmpjp
      real (dp)       , dimension(:),   allocatable :: tmpv
      real (8)        , dimension(:),   allocatable :: tmpMW, tmp_ELMW
      real (8)        , dimension(:),   allocatable :: taL,tbL,tcL,   &
                                                       tdL,teL,tfL,tgL
      real (8)        , dimension(:),   allocatable :: taH,tbH,tcH,   &
                                                       tdH,teH,tfH,tgH
      real (8)        , dimension(:)  , allocatable :: tmpA, tmpb,    &
                                                       tmpE
      real (8)        , dimension(:,:), allocatable :: tmptemp, tmprea
      real (8)        , dimension(:,:), allocatable :: tmpfar, tmpmol
      real (dp)       , dimension(:,:), allocatable :: third_body_dense

      type(sparseint)                               :: tmp_isp

!     *****************************************************************

!ck2015 mechdir add 3 lines
      cklink = trim(mechdir)//"cklink"
      cklin2 = trim(mechdir)//"chem.bin"
      output = trim(mechdir)//"SpeedCHEM.out"
!     Set physical gas constant
!     [FP] For compatibility with different versions of the CHEMKIN
!     reaction mechanism interpreter, one might use the approximate
!     Rcal = 1.987 kcal/mol/K value instead than the more accurate
!     ratio between the gas constant and the conversion factor
!      Rcal = 8.314_dp/4.184_dp
      Rcal = 1.987_dp

!     ** Open output file for mechanism processing
      open(unit=lout,file=output,status='unknown')

      write(lout,fmt_head )
      write(lout,fmt_head2)
      write(lout,fmt_head )
      write(lout,        *)

!     ** Check presence of chemkin linking file
      inquire(file = cklink,exist = ispresent )
      if (ispresent ) ckfile = adjustl(cklink)
      inquire(file = cklin2,exist = ispresent2)
      if (ispresent2) ckfile = adjustl(cklin2)

!     ** Halt on missing file
      if ((.not.ispresent) .and. (.not.ispresent2)) then
!ck2015         write(*   ,fmt_nofile)
!ck2015         write(lout,fmt_nofile)
         write(*   ,fmt_nofile) trim(cklink)," ", trim(cklin2)
         write(lout,fmt_nofile) trim(cklink)," ", trim(cklin2)
         stop
      endif

!      write(*   , fmt_fname)trim(adjustl(ckfile))
      write(lout, fmt_fname)trim(adjustl(ckfile))
      write(lout, *        )

!     Open the linking file for reading. Use action='read'/status='old'
!     and capture iostat so a missing/unreadable file fails closed with a
!     message instead of silently proceeding (or creating an empty file,
!     as status='unknown' without action='read' could).
      open(unit=linc, file = ckfile,  form   = 'unformatted', &
                                      status = 'old',         &
                                      action = 'read',        &
                                      iostat = idummy1        )
      if (idummy1 /= 0) then
         write(*   , "(' ERROR: cannot open chemkin linking file: ',A,"//&
                     "' (iostat=',I0,'). Aborting.')") &
                     trim(adjustl(ckfile)), idummy1
         write(lout, "(' ERROR: cannot open chemkin linking file: ',A,"//&
                     "' (iostat=',I0,'). Aborting.')") &
                     trim(adjustl(ckfile)), idummy1
         error stop 1
      endif
      rewind linc

!     ** PRELIMINARY PROBLEM DIMENSIONS AND DATA ***********************

!     Mechanism label
      mechanism = trim(adjustl("Mechanism imported from "//ckfile))

!     ** cklink v2 magic + schema (fail-closed) **********************
!     CKINTP writes a leading record {magic, schema} before the legacy
!     {VERS,PREC,KERR} header. This positively identifies the on-disk
!     format and version, replacing the fragile VERS-string heuristic.
!     A missing/wrong magic means an old-format or foreign cklink whose
!     record layout we must not assume -> refuse (regenerate).
      read (linc, iostat=idummy1) ck_magic_in, ck_schema_in
      if (idummy1 /= 0) then
         write(*   , "(' ERROR: cklink is not a SpeedCHEM v2 linking',"//   &
                     "' file (missing magic record). Regenerate cklink',"//&
                     "' from chem.inp/therm.dat. Aborting.')")
         write(lout, "(' ERROR: cklink is not a SpeedCHEM v2 linking',"//   &
                     "' file (missing magic record). Aborting.')")
         error stop 1
      endif
      if (ck_magic_in /= CK_MAGIC_EXPECT .or.                          &
          ck_schema_in /= CK_SCHEMA_EXPECT) then
         write(*   , "(' ERROR: cklink magic/schema mismatch (got [',A,"//&
                     "'] v',I0,', expected [',A,'] v',I0,'). Regenerate',"//&
                     "' cklink. Aborting.')") ck_magic_in, ck_schema_in,  &
                     CK_MAGIC_EXPECT, CK_SCHEMA_EXPECT
         write(lout, "(' ERROR: cklink magic/schema mismatch (got [',A,"//&
                     "'] v',I0,', expected [',A,'] v',I0,'). Aborting.')")&
                     ck_magic_in, ck_schema_in, CK_MAGIC_EXPECT,          &
                     CK_SCHEMA_EXPECT
         error stop 1
      endif

!     File header: Chemkin version, machine precision and check OK
      read (linc) ckvers, ckprec, kerr
      write(lout, fmt_info)ckvers,ckprec

!     ** Fail-closed on a linking file the interpreter flagged as bad ***
!     CKINTP writes the header (VERS, PREC, KERR) even when the mechanism
!     input had errors, then STOPs -- leaving a header-only, truncated
!     cklink behind. Reading past the header here would give a garbled
!     mechanism, so refuse when KERR is set.
      if (kerr) then
         write(*   , "(' ERROR: chemkin linking file reports interpreter',"//&
                     "' errors (KERR=.true.). Regenerate cklink from',"//     &
                     "' chem.inp/therm.dat. Aborting.')")
         write(lout, "(' ERROR: chemkin linking file reports interpreter',"//&
                     "' errors (KERR=.true.). Regenerate cklink from',"//     &
                     "' chem.inp/therm.dat. Aborting.')")
         error stop 1
      endif

!     ** Inner CHEMKIN version (informational) ***********************
!     The magic/schema record above is now the authoritative format
!     gate; ckvers is the inner CHEMKIN-II data version ('3.1'), kept
!     for the log only.
      if (trim(adjustl(ckvers)) /= "3.1") then
         write(lout, "(' NOTE: inner chemkin data version ',A,"//         &
                     "' (expected 3.1).')") trim(adjustl(ckvers))
      endif

!     ** Actions depending on data precision ***************************
      if (ckprec(1:6) /= "DOUBLE") then
         write(lout, fmt_prec)
         stop
      endif


!     Read integer problem dimensions
!     len*work   = lengths of working arrays: integer, dp, character
!     mm, kk, ii = number of elements, species, reactions
!     maxsp      = max number of species per reaction (reac+prod)
!     maxtb      = max number of enhanced third-vody efficiencies
!     maxtp      = max number of temperature bounds for thermo
!                  (usually, maxtp=3, i.e. 2 temperature intervals)
!     nthcf      = number of polynomial coefficients for Cp
!     nipar      = number of parameters in Arrhenius expression (=3)
!     nitar      = number of parameters for landau-teller reactions
!     nifar      = number of parameters for falloff reactions
!     nrv        = number of reactions with explicit reverse rates
!     nfl        = number of falloff reactions
!     nthb       = number of third-body reactions
!     nlt        = number of landau-teller reactions
!     nrl        = number of landau-teller reactions with explicit
!                  reverse parameters
!     nw         = number of reactions with radiation wavelength enhanc.
!     nchrg      = number of species with non-zero charge
      read (linc) liwork, lrwork, lcwork, mm, kk, ii, maxsp, maxtb,   &
                  maxtp, nthcf, nipar, nitar, nifar, nrv, nfl, nthb,  &
!ck2015                  nlt, nrl, nw, nchrg
!     Must mirror the CKINTP header WRITE exactly (chemkin_module.f90):
!       ... NEIM, NJAR, NJAN, NF1R, NFT1, NEXC, NRNU, NORD, MAXORD, CKMIN
!     The original read collapsed NF1R and NFT1 into one `nf1` and dropped
!     NORD/MAXORD/CKMIN, so modified species orders (FORD/RORD) could not
!     even be counted, let alone rejected. nif1=NF1R, nf1=NFT1.
               nlt, nrl, nw, nchrg, nei, nja, njan, nif1, nf1, nex,   &
                  nsto, nord, maxord, ckmin
!     Assign dimensions to speedchem data
      nel = mm
      ns  = kk
      nr  = ii

!     Allocate species and reaction indices
      allocate(species(ns), reactions(nr))
      species   = [(j,j=1,ns)]
      reactions = [(j,j=1,nr)]

!     Number of detailed chemistry problem equations
      neq  = ns + 1
      uneq = 1.e0_dp / neq

      write(lout, *      )
      write(lout, fmt_nel)nel
      write(lout, fmt_ns )ns
      write(lout, fmt_nr )nr
      write(lout, *      )

!     Check for presence of reaction types not supported by this code
      if (nlt > 0 .or. nrl > 0) then
         write(*   , fmt_lt)nlt
         write(lout, fmt_lt)nlt
         stop
      endif
      if (nw > 0 .or. nchrg > 0) then
         write(*   , fmt_chrg)
         write(lout, fmt_chrg)
         stop
      endif
!     Fail-closed on reaction features that CKINTP parses into cklink but
!     that this reader/rate code does NOT convert into an evaluated rate
!     form. Previously these records were silently skipped, so a mechanism
!     using them would link and integrate while quietly ignoring the
!     feature. Refuse instead, naming the keyword and reaction count.
      if (nei > 0) then
         write(*   , fmt_unsup) nei, 'electron-impact temperature dep. (EIM/TDEP)'
         write(lout, fmt_unsup) nei, 'electron-impact temperature dep. (EIM/TDEP)'
         error stop 1
      endif
      if (njan > 0) then
         write(*   , fmt_unsup) njan, 'Jannev-Langer-Evans-Post (JAN)'
         write(lout, fmt_unsup) njan, 'Jannev-Langer-Evans-Post (JAN)'
         error stop 1
      endif
      if (nf1 > 0) then
         write(*   , fmt_unsup) nf1, 'fit #1 (FIT1)'
         write(lout, fmt_unsup) nf1, 'fit #1 (FIT1)'
         error stop 1
      endif
      if (nex > 0) then
         write(*   , fmt_unsup) nex, 'excitation energy loss (EXCI)'
         write(lout, fmt_unsup) nex, 'excitation energy loss (EXCI)'
         error stop 1
      endif
      if (nord > 0) then
         write(*   , fmt_unsup) nord, 'modified species order (FORD/RORD)'
         write(lout, fmt_unsup) nord, 'modified species order (FORD/RORD)'
         error stop 1
      endif
      if (nsto > 0) then
         write(*   , fmt_unsup) nsto, 'real (non-integer) stoichiometry (REAL)'
         write(lout, fmt_unsup) nsto, 'real (non-integer) stoichiometry (REAL)'
         error stop 1
      endif

      if (maxtp /= 3 .or. nthcf /= 5) then
         write(*   , fmt_tmp)maxtp, nthcf
         write(lout, fmt_tmp)maxtp, nthcf
         stop
      endif

!     Allocate mechanism arrays
      call SCallocate


!     ** ELEMENT-RELATED DATA ******************************************

!     Read element names and atomic weights
      allocate( tmpchar(nel), tmp_elMW(nel) )

      read(linc)(tmpchar(i), tmp_elMW(i), i = 1,nel)
      ! Convert into current working format
      elMW  = real(tmp_elMW, kind=dp)
      uelMW = one/elMW

      write(lout, *)
      write(lout, fmt_eldata)

      ! Element names
      do i = 1, nel
         elementi(i) = tmpchar(i)(1:2)
         write(lout, fmt_elit)i, elementi(i), elMW(i)
      end do

      write(lout, *)
      deallocate(tmpchar, tmp_elMW)

!     ** SPECIES-RELATED DATA ******************************************

!     Read species names
      allocate(EM(ns,nel))
      allocate(tmpchar(ns), idummy(ns), tmptemp(maxtp,ns))

      allocate(tmpMW(ns))
      allocate(taL(ns),tbL(ns),tcL(ns),tdL(ns),teL(ns),tfL(ns),tgL(ns))
      allocate(taH(ns),tbH(ns),tcH(ns),tdH(ns),teH(ns),tfH(ns),tgH(ns))

      read (linc) (tmpchar(i),               & ! species name
                  (EM(i,j),j=1,nel),         & ! species data in EM
                   idummy1,                  & ! phase (gas only)
                   idummy2,                  & ! charge
                   tmpMW(i),                 & ! molecular weight [g/mol]
                   idummy(i),                & ! Num. of temperature pts
                  (tmptemp(k,i), k=1,maxtp), & ! Thermo temp. bounds [K]
                   taL(i), tbL(i), tcL(i),   & ! Thermo coefficients
                   tdL(i), teL(i), tfL(i),   &
                   tgL(i), taH(i), tbH(i),   &
                   tcH(i), tdH(i), teH(i),   &
                   tfH(i), tgH(i),           &
                   i = 1, ns )                 ! Main loop is on species
!
!      write(*,*)  (tmpchar(i),               & ! species name
!                  (EM(i,j),j=1,nel),         & ! species data in EM
!                   idummy1,                  & ! phase (gas only)
!                   idummy2,                  & ! charge
!                   tmpMW(i),                 & ! molecular weight [g/mol]
!                   idummy(i),                & ! Num. of temperature pts
!                  (tmptemp(k,i), k=1,maxtp), & ! Thermo temp. bounds [K]
!                   taL(i), tbL(i), tcL(i),   & ! Thermo coefficients
!                   tdL(i), teL(i), tfL(i),   &
!                   tgL(i), taH(i), tbH(i),   &
!                   tcH(i), tdH(i), teH(i),   &
!                   tfH(i), tgH(i),           &
!                   i = 1, ns )

!     Convert data into current working format
      SCMW = real(tmpMW, kind = dp)
      aL   = real(taL  , kind = dp)
      bL   = real(tbL  , kind = dp)
      cL   = real(tcL  , kind = dp)
      dL   = real(tdL  , kind = dp)
      eL   = real(teL  , kind = dp)
      fL   = real(tfL  , kind = dp)
      gL   = real(tgL  , kind = dp)
      aH   = real(taH  , kind = dp)
      bH   = real(tbH  , kind = dp)
      cH   = real(tcH  , kind = dp)
      dH   = real(tdH  , kind = dp)
      eH   = real(teH  , kind = dp)
      fH   = real(tfH  , kind = dp)
      gH   = real(tgH  , kind = dp)



      write(lout, *       )
      write(lout, fmt_spec)
      species_fix: do i = 1, ns

!        Fix species names
         specie(i) = trim(adjustl(tmpchar(i)))
         write(lout, fmt_spdet)i, specie(i), SCMW(i)

!        Thermodynamic database switch temperature
!        NB! Check if temperatures are always in ascending order!
         tsw(i) = real(tmptemp(2,i), kind = dp)

      end do species_fix
      write(lout, *        )

!     Reciprocal of species molecular weights [mol/g]
      uMW = one / SCMW

 !    Print thermodynamic database in SpeedCHEM format (SCthermo.dat)
 !    to file lout
      write(lout, fmt_thhead)
      write(lout, *         )'nrecords ',ns
      thermo_database: do i = 1, ns

        write(lout, *       ) specie(i)
        write(lout, *       ) 'HI TEMP'
        write(lout, fmt_thl1) aH(i),bH(i),cH(i),dH(i)
        write(lout, fmt_thl1) eH(i),fH(i),gH(i)
        write(lout, *       ) 'TSWITCH', tsw(i)
        write(lout, fmt_thl1) aL(i),bL(i),cL(i),dL(i)
        write(lout, fmt_thl1) eL(i),fL(i),gL(i)

      end do thermo_database

      write(lout, *         )



      deallocate(tmpchar, tmptemp)
      deallocate(tmpMW)
      deallocate(taL, tbL, tcL, tdL, teL, tfL, tgL)
      deallocate(taH, tbH, tcH, tdH, teH, tfH, tgH)

!     Check that all of the thermodynamic parameter values are stored
!     according to the standard JANAF format with 2 temperature
!     intervals only
      if (any(idummy /= 3)) then
         stop
      endif

      deallocate(idummy)

!     Write element structure to output file


!     ** REACTION-RELATED DATA *****************************************

      allocate( idummy0(nr), idummy(nr), tmprea(nipar,nr),             &
                itmpst(maxsp,nr), itmpsp(maxsp,nr) )


      tmprea = zero

!     Read main reaction arrays
      read(linc)( idummy0(k),                     & ! Number of species involved (irreversible if < 0)
                  idummy (k),                     & ! Number of reactants
                 (tmprea(i,k), i = 1, nipar),     & ! Reaction main Arrhenius parameters
                 (itmpst(j,k),                    & ! Stoichiometric coefficients of species involved in this reaction
                  itmpsp(j,k), j = 1, maxsp),     & ! Species indexes for the stoichiometric coefs
                               k = 1, nr    )

!     Assign reactions_related arrays
!     NB it appears that Chemkin stores the high pressure limit values in
!     this array when a falloff reaction is present, however,
!     need to check this
      Ainf = real(tmprea(1,:), kind = dp)
      binf = real(tmprea(2,:), kind = dp)
      Einf = real(tmprea(3,:), kind = dp) * Rcal

!     Reaction reversibility indices (to be updated later with explicit
!     reverse reactions data!)
      where (idummy0 < 0)
         reversibile = .false.
      elsewhere
         reversibile = .true.
      end where

      deallocate(idummy0, idummy)

      assign_stoichiometric_coefficients: do k = 1, nr

        do j = 1, maxsp

!          A reactant is present
           if (itmpst(j,k) < 0) then
              call add_value(stoich_r_sp,k,itmpsp(j,k),-real(itmpst(j,k),dp))
!          A product is present
           elseif (itmpst(j,k) > 0) then
              call add_value(stoich_p_sp,k,itmpsp(j,k), real(itmpst(j,k),dp))
           endif

        end do

      end do assign_stoichiometric_coefficients
!ck2015s for real stoichometric coefficients
      if(nsto > 0) then
!        FAIL-CLOSED: real (non-integer) stoichiometric coefficients are read
!        here but the downstream rate/Jacobian code truncates them via int()
!        (see iv_stoich_p/r assignment below and mass_action in SCmodule),
!        so they are NOT correctly evaluated. Refuse rather than silently
!        compute a wrong mechanism.
         write(*   , "(' ERROR: ',I5,' reactions use real (non-integer)',"// &
                     "' stoichiometric coefficients, which SpeedCHEM does',"//&
                     "' not yet evaluate correctly. Aborting.')") nsto
         write(lout, "(' ERROR: ',I5,' reactions use real (non-integer)',"// &
                     "' stoichiometric coefficients, which SpeedCHEM does',"//&
                     "' not yet evaluate correctly. Aborting.')") nsto
         error stop 1
      allocate( idummy0(nsto), rtmpst(maxsp,nsto) )
         if(nrv > 0) read(linc)
         if(nfl > 0) read(linc)
         if(nthb > 0) read(linc)
         if(nlt > 0) read(linc)
         if(nrl > 0) read(linc)
         if(nw > 0) read(linc)
         if(nei > 0) read(linc)
         if(njan > 0) read(linc)
         if(nf1 > 0) read(linc)
         if(nex > 0) read(linc)
         read(linc) ( idummy0(k),  &
                     (rtmpst(j,k), j = 1, maxsp),  &
                                   k = 1, nsto)
!        NOTE: rtmpst is packed as (maxsp, nsto); its 2nd index is the packed
!        record index i (1..nsto), NOT the global reaction number k. The
!        original code indexed rtmpst(j,k) with k=idummy0(i), which reads the
!        wrong record and can go out of bounds when k>nsto. Fixed to (j,i).
      assign_real_stoichiometric_coefficients: do i = 1, nsto
        k = idummy0(i)
        do j = 1, maxsp
!          A reactant is present
           if (rtmpst(j,i) < 0) then
               call add_value(stoich_r_sp,k,itmpsp(j,k),-rtmpst(j,i))
!          A product is present
           elseif (rtmpst(j,i) > 0) then
              call add_value(stoich_p_sp,k,itmpsp(j,k), rtmpst(j,i))
           endif
        end do
      end do assign_real_stoichiometric_coefficients

      rewind linc
      read(linc)
      read(linc)
      read(linc)
      read(linc)
      read(linc)

      deallocate( idummy0, rtmpst )

      endif
!ck2015e

!     Fix sparse matrix representation for unactive species
      stoich_r_sp%nc = ns
      stoich_p_sp%nc = ns

!     Assign nudiff
!      nudiff_sparse  = sparse_sum( 1e0_dp, stoich_p_sp,&
!                                  -1e0_dp, stoich_r_sp )
      nudiff_sparse  = stoich_p_sp - stoich_r_sp

!     Count maximum reaction dimensions (products or reactants)
      allocate(count_reacs(nr))

      call sparse_internal_count(stoich_r_sp,count_reacs,dim=2)
      max_nreacs = maxval(count_reacs)

      call sparse_internal_count(stoich_p_sp,count_reacs,dim=2)
      max_nprods = maxval(count_reacs)

      deallocate(count_reacs)

!     Initialising nudiff
!      nudiff = stoich_p - stoich_r
      tmp_isp = nudiff_sparse
      if (.not.allocated(isumnudiff)) then
         allocate(isumnudiff(nr))
         call sparseint_internal_sum(tmp_isp,isumnudiff,2)
!         isumnudiff = int(sum(nudiff, dim = 2))
      endif

      call spdeallocate(tmp_isp)
      deallocate  (itmpst, itmpsp, tmprea)

!
!!     PERTURBATION FACTOR............ to be understood
!
!      DO 10 I = 1, II
!         RCKWRK(NcCO + (I-1)*(NPAR+1) + NPAR) = 1.0
!   10 CONTINUE
!

!      Reactions with explicit reverse Arrhenius parameters
!      ire = reaction number
!      tmpA = is the pre-exponential factor (mole-cm-sec-K),
!      tmpb = is the temperature exponent, and
!      tmpE = is the activation energy (Kelvins).
       REV  = .false.
       Arev = 0e0_dp
       brev = 0e0_dp
       Erev = 0e0_dp

       explicit_reverse_arrhenius: if (nrv > 0) then
          allocate (ire(nrv), tmpA(nrv), tmpE(nrv), tmpb(nrv))

!         Read data from binary file
          read(linc) (ire(k), tmpA(k), tmpb(k), tmpE(k), k = 1, nrv)

!         Assign data to SpeedCHEM
          Arev(ire) = real(tmpA, kind = dp)
          brev(ire) = real(tmpb, kind = dp)
          Erev(ire) = real(tmpE, kind = dp) * Rcal! Conversion from [K] to [cal/mole]

!         Check about masked irreversible reactions that are treated
!         as reversible ones, but have explicit reverse coefficients
!         set to zero
          where (Arev(ire) == 0e0_dp) reversibile(ire) = .false.

          deallocate(ire, tmpA, tmpE, tmpb)
       end if explicit_reverse_arrhenius

!     Flag for explicit reversible reaction rates
      REV = Arev /= 0e0_dp
!      REV = (Arev/=0e0_dp .or. brev/=0e0_dp .or. Erev/=0e0_dp)


!     Packed array of non reversible reactions
      allocate(reac_cond(nr))
      if (.not.allocated(inotrev)) then
         reac_cond = .NOT.reversibile
         nnotrev = count(reac_cond)
         if (nnotrev > 0) then
            allocate(inotrev(nnotrev))
            inotrev = pack(reactions, reac_cond)
         endif
      endif


!      Falloff reactions
!      reaction type: 1 for 3-parameter Lindemann Form
!                     2 for 6- or 8-parameter SRI Form (not supported)s
!                     3 for 6-parameter Troe Form
!                     4 for 7-parameter Troe form
!      Molecularity:  0 for total gas concentration
!                     K for concentration of the species K only

       A0   = 0e0_dp
       b0   = 0e0_dp
       E0   = 0e0_dp

       HIGH = .false.
       LOW  = .false.
       TROE = .false.

!      Assign falloff reaction data depending on type of reaction
       aTROE  = 0e0_dp
       T3TROE = 0e0_dp
       T1TROE = 0e0_dp
       T2TROE = 0e0_dp
!       third_body(:,:) = 0e0_dp
!       third_body_beta = 0e0_dp
       call allocate(nr,     ns, 0, third_body_sp)



       falloff_reactions: if (nfl > 0) then

          allocate(ire(nfl), itype(nfl), imole(nfl), tmpfar(nifar,nfl) )

!         Read data from binary file
          read(linc) (ire(k),                   & ! Falloff reaction indx
                      itype(k),                 & ! Reaction falloff type
                      imole(k),                 & ! Molecularity
                     (tmpfar(j,k), j=1, nifar), & ! parameter values
                                   k=1, nfl     )

!         Halt on unsupported reactions
          if (any(itype == 2)) then
             write(*, fmt_sri)count(itype==2)
             stop
          endif

!         Logical switch for presence of low pressure limit
          LOW(ire) = .true.

          if (nfl > 0) then
          do i = 1, nfl

!            Molecularity data
!             if (imole(i) == 0) third_body(ire(i),:) = 1e0_dp
!             if (imole(i) /= 0) third_body(ire(i),imole(i)) = 1e0_dp

             if (imole(i) == 0) &
             call add_line(third_body_sp,ire(i),[(1e0_dp,j=1,ns)])

             if (imole(i) /= 0) &
             call add_value(third_body_sp,ire(i),imole(i),1e0_dp)

             select case (itype(i))

                case (1) ! LINDEMANN form, 3 Arrhenius parameters

                    A0(ire(i)) = real(tmpfar(1,i), dp)
                    b0(ire(i)) = real(tmpfar(2,i), dp)
                    E0(ire(i)) = real(tmpfar(3,i), dp) * Rcal

                case (3) ! TROE form with 3 parameters + 3 Arrhenius

                    A0(ire(i)) = real(tmpfar(1,i), dp)
                    b0(ire(i)) = real(tmpfar(2,i), dp)
                    E0(ire(i)) = real(tmpfar(3,i), dp) * Rcal

                    TROE  (ire(i)) = .TRUE.
                    aTROE (ire(i)) = real(tmpfar(4,i), dp)
                    T3TROE(ire(i)) = real(tmpfar(5,i), dp)
                    T1TROE(ire(i)) = real(tmpfar(6,i), dp)


                case (4) ! TROE form with 4 parameters + 3 Arrhenius

                    A0(ire(i)) = real(tmpfar(1,i), dp)
                    b0(ire(i)) = real(tmpfar(2,i), dp)
                    E0(ire(i)) = real(tmpfar(3,i), dp) * Rcal

                    TROE  (ire(i)) = .TRUE.
                    aTROE (ire(i)) = real(tmpfar(4,i), dp)
                    T3TROE(ire(i)) = real(tmpfar(5,i), dp)
                    T1TROE(ire(i)) = real(tmpfar(6,i), dp)

                    ! Prevent T2TROE from being zero for CHEMKIN compatibility
                    ! (it can happen that someone wants T2TROE == 0.0)
                    if (abs(real(tmpfar(7,i),dp))<tiny(zero)) then
                       T2TROE(ire(i)) = 1.0e-99_dp
                    else
                       T2TROE(ire(i)) = real(tmpfar(7,i), dp)
                    endif


                case default

                   write(*, fmt_erfl)itype(i)
                   stop

             end select

          end do
          endif

          deallocate(ire, itype, imole, tmpfar)

       end if falloff_reactions


!     Allocating variables for TROE pressure-dependent reactions
      if (.not.allocated(todotroe) ) then

        allocate(todotroe(count(TROE)))
        allocate(aT2     (count(TROE)))
        allocate(uT1T2   (count(TROE)))
        allocate(T2T2    (count(TROE)))
        allocate(uT3T2   (count(TROE)))
!        allocate(expuT3T2(count(TROE)))
!        allocate(expuT1T2(count(TROE)))

        todotroe = pack(reactions, TROE)

        aT2      = aTROE(todotroe)
        uT1T2    = one/T1TROE(todotroe)
        T2T2     = T2TROE(todotroe)
        uT3T2    = one/T3TROE(todotroe)
!        expuT1T2 = exp(uT1T2)
!        expuT3T2 = exp(uT3T2)


        allocate( zeroT2(count(T2T2 == zero)) )
        zeroT2 = pack([(j,j=1,count(TROE))], T2T2 == zero)

       endif



!      Assign third body species enhanced coefficients
       third_body_data: if (nthb > 0) then

           allocate(ire(nthb), inthb(nthb), itmpsp(maxtb,nthb), tmpmol(maxtb,nthb))
           itmpsp = 0
           tmpmol = 0e0_dp

!          Read third-body related data
           read (linc) ( ire(i),        &  ! Number of reaction with third bodies
                         inthb(i),      &  ! Number of enhanced efficiencies
                        (itmpsp(k,i),   &  ! Index of enhancing species
                         tmpmol(k,i),   &  ! Enhanced molecularity coefficient
                         k = 1, maxtb), &
                         i = 1, nthb )

!          Assign third-body reaction data
           tb_assign: do i = 1, nthb

!               third_body(ire(i),:) = 1e0_dp

               call add_line(third_body_sp,ire(i),[(one,j=1,ns)])

               if (inthb(i) > 0) then
                  do j = 1, inthb(i)
                     call add_value(third_body_sp,ire(i),itmpsp(j,i),real(tmpmol(j,i), dp))
                  end do
               endif

           end do tb_assign

!          Logical mask of third-body reactions
           THREE(ire) = .true.

           deallocate(ire, inthb, itmpsp, tmpmol)

       endif third_body_data

!      ** cklink v2 PLOG SECTION (fail-closed on truncation) ***********
!      Positioned here because every unsupported optional section
!      (LAN/RLT/WL/EIM/JAN/FT1/EXC/RNU/ORD) is rejected earlier with
!      error stop, so for any mechanism that reaches this point only
!      REV/FAL/THB were written and consumed above -- the file is now at
!      the PLOG counts record (see CKINTP write order). Layout:
!        counts -> [reaction map] -> [node arrays] -> checksum.
       read(linc, iostat=idummy1) n_plog_reactions, n_plog_nodes
       if (idummy1 /= 0) then
          write(*   , "(' ERROR: cklink v2 PLOG section missing/truncated',"//&
                      "' (counts record). Regenerate cklink. Aborting.')")
          write(lout, "(' ERROR: cklink v2 PLOG section missing/truncated',"//&
                      "' (counts record). Aborting.')")
          error stop 1
       endif
       if (n_plog_reactions < 0 .or. n_plog_nodes < 0 .or.             &
           (n_plog_reactions == 0 .and. n_plog_nodes /= 0) .or.        &
           (n_plog_reactions > 0 .and.                                &
            n_plog_nodes < n_plog_reactions)) then
          write(*   , "(' ERROR: invalid cklink v2 PLOG counts: ',"//  &
                      "'reactions=',I0,', nodes=',I0,'. Aborting.')")   &
                      n_plog_reactions, n_plog_nodes
          write(lout, "(' ERROR: invalid cklink v2 PLOG counts.')")
          error stop 1
       endif
!      cklink v2 stores one entry per Arrhenius term. Adjacent entries
!      may share a pressure; the evaluator groups and sums them.
       n_plog_terms = n_plog_nodes

!      Defensive: make SCcklink re-entrant (a second call in the same
!      process must not hit "already allocated").
       if (allocated(plog_reaction))  deallocate(plog_reaction)
       if (allocated(plog_node_ptr))  deallocate(plog_node_ptr)
       if (allocated(plog_term_ptr))  deallocate(plog_term_ptr)
       if (allocated(plog_logP))      deallocate(plog_logP)
       if (allocated(plog_A))         deallocate(plog_A)
       if (allocated(plog_b))         deallocate(plog_b)
       if (allocated(plog_EoverR))    deallocate(plog_EoverR)

       if (n_plog_reactions > 0) then
          allocate(plog_reaction(n_plog_reactions),                    &
                   plog_node_ptr(0:n_plog_reactions))
          allocate(plog_logP(n_plog_nodes), plog_A(n_plog_nodes),      &
                   plog_b(n_plog_nodes), plog_EoverR(n_plog_nodes))
!         v2 entry-level storage: term_ptr(k)=k. Equal-pressure entries
!         are grouped without changing the backward-compatible layout.
          allocate(plog_term_ptr(0:n_plog_nodes))
          plog_term_ptr = [(ipl, ipl = 0, n_plog_nodes)]

          read(linc, iostat=idummy1)                                   &
               (plog_reaction(ipl), ipl = 1, n_plog_reactions)
          if (idummy1==0) read(linc, iostat=idummy1)                   &
               (plog_node_ptr(ipl), ipl = 0, n_plog_reactions)
          if (idummy1==0) read(linc, iostat=idummy1)                   &
               (plog_logP(ipl),   ipl = 1, n_plog_nodes)
          if (idummy1==0) read(linc, iostat=idummy1)                   &
               (plog_A(ipl),      ipl = 1, n_plog_nodes)
          if (idummy1==0) read(linc, iostat=idummy1)                   &
               (plog_b(ipl),      ipl = 1, n_plog_nodes)
          if (idummy1==0) read(linc, iostat=idummy1)                   &
               (plog_EoverR(ipl), ipl = 1, n_plog_nodes)
          if (idummy1 /= 0) then
             write(*   , "(' ERROR: cklink v2 PLOG arrays truncated.',"//  &
                         "' Regenerate cklink. Aborting.')")
             write(lout, "(' ERROR: cklink v2 PLOG arrays truncated.',"//  &
                         "' Aborting.')")
             error stop 1
          endif

!         Validate the complete packed topology before any array is used.
!         This turns corrupt-but-readable cklink records into a controlled
!         failure instead of an out-of-bounds access in the evaluator.
          if (plog_node_ptr(0) /= 0 .or.                               &
              plog_node_ptr(n_plog_reactions) /= n_plog_nodes) then
             write(*   , "(' ERROR: invalid cklink v2 PLOG node',"//  &
                         "' pointers. Regenerate cklink. Aborting.')")
             write(lout, "(' ERROR: invalid cklink v2 PLOG node pointers.')")
             error stop 1
          endif
          do ipl = 1, n_plog_reactions
             invalid_plog = plog_reaction(ipl) < 1 .or.                &
                            plog_reaction(ipl) > nr .or.                &
                            plog_node_ptr(ipl) <= plog_node_ptr(ipl-1)
             if (ipl > 1) invalid_plog = invalid_plog .or.             &
                  plog_reaction(ipl) <= plog_reaction(ipl-1)
             if (invalid_plog) then
                write(*   , "(' ERROR: invalid cklink v2 PLOG reaction',"//&
                            "' map/pointers at packed reaction ',I0,"//   &
                            "'. Aborting.')") ipl
                write(lout, "(' ERROR: invalid cklink v2 PLOG reaction map.')")
                error stop 1
             endif
             do inode = plog_node_ptr(ipl-1)+1, plog_node_ptr(ipl)
                invalid_plog = .not. ieee_is_finite(plog_logP(inode)) .or. &
                               .not. ieee_is_finite(plog_A(inode)) .or.    &
                               .not. ieee_is_finite(plog_b(inode)) .or.    &
                               .not. ieee_is_finite(plog_EoverR(inode)) .or.&
                               .not. plog_A(inode) > zero
                if (inode > plog_node_ptr(ipl-1)+1)                     &
                   invalid_plog = invalid_plog .or.                     &
                      plog_logP(inode) < plog_logP(inode-1) - 1.0e-9_dp
                if (invalid_plog) then
                   write(*   , "(' ERROR: invalid cklink v2 PLOG node ',"//&
                               "I0,' for reaction ',I0,'. Aborting.')")  &
                               inode, plog_reaction(ipl)
                   write(lout, "(' ERROR: invalid cklink v2 PLOG node.')")
                   error stop 1
                endif
             enddo
          enddo
       else
!         No PLOG reactions: leave arrays unallocated (size 0 conceptually).
          allocate(plog_reaction(0), plog_node_ptr(0:0))
          plog_node_ptr(0) = 0
          allocate(plog_logP(0), plog_A(0), plog_b(0), plog_EoverR(0))
          allocate(plog_term_ptr(0:0)); plog_term_ptr(0) = 0
       endif

!      Section checksum (matches CKINTP): reactions + nodes + last ptr.
       read(linc, iostat=idummy1) plog_chk_in
       if (idummy1 /= 0) then
          write(*   , "(' ERROR: cklink v2 PLOG checksum record missing.',"//&
                      "' Regenerate cklink. Aborting.')")
          write(lout, "(' ERROR: cklink v2 PLOG checksum record missing.')")
          error stop 1
       endif
       plog_chk_calc = n_plog_reactions + n_plog_nodes
       if (n_plog_reactions > 0)                                       &
          plog_chk_calc = plog_chk_calc + plog_node_ptr(n_plog_reactions)
       if (plog_chk_in /= plog_chk_calc) then
          write(*   , "(' ERROR: cklink v2 PLOG checksum mismatch (got',"//&
                      "' ',I0,', expected ',I0,'). Corrupt/misaligned',"//&
                      "' cklink. Aborting.')") plog_chk_in, plog_chk_calc
          write(lout, "(' ERROR: cklink v2 PLOG checksum mismatch.')")
          error stop 1
       endif

!      Standard gas-phase PLOG cannot be combined with an explicit REV
!      rate, falloff syntax, TROE, or another pressure form.
!      These combinations require different mathematics and must never
!      fall through to the ordinary PLOG evaluator.
       do ipl = 1, n_plog_reactions
          k = plog_reaction(ipl)
          if (REV(k) .or. HIGH(k) .or. LOW(k) .or. TROE(k)) then
             write(*   , "(' ERROR: PLOG reaction ',I0,' has',"//      &
                         "' unsupported REV/HIGH/LOW/TROE flags: ',"// &
                         "4(L1,1X))")                                   &
                         k, REV(k), HIGH(k), LOW(k), TROE(k)
             write(lout, "(' ERROR: unsupported PLOG combination at',"//&
                          "' reaction ',I0,'.')") k
             error stop 1
          endif
       enddo

!      CKINTP's legacy linking layout classifies PLOG reactions in the
!      third-body table even when the CHEMKIN equation contains no +M.
!      PLOG already includes pressure dependence in k(T,P); retaining
!      those rows would multiply the rate by [M] a second time. Rebuild
!      the matrix without PLOG rows. Explicit +M/falloff was rejected
!      while parsing, so every removed row is the legacy PLOG marker.
       if (n_plog_reactions > 0) then
          allocate(third_body_dense(nr,ns))
          third_body_dense = zero
          do i = 1, nr
             if (any(plog_reaction == i)) cycle
             do j = 1, ns
                third_body_dense(i,j) = sparse_value(third_body_sp,i,j)
             enddo
          enddo
          call deallocate(third_body_sp)
          call allocate(nr, ns, 0, third_body_sp)
          do i = 1, nr
             do j = 1, ns
                if (third_body_dense(i,j) /= zero)                     &
                   call add_value(third_body_sp,i,j,third_body_dense(i,j))
             enddo
          enddo
          deallocate(third_body_dense)
       endif

!      Build the per-reaction rate-form tag. Legacy code paths still use
!      Arrhreac/Lindreac/Troereac; rate_form is stage-1 plumbing that
!      records which reactions are PLOG (RATE_PLOG). It is populated for
!      every reaction so stage 2+ can dispatch on it. A no-PLOG mechanism
!      leaves it all-zero and unused, so nothing numeric changes.
       if (.not. allocated(rate_form)) allocate(rate_form(nr))
       rate_form = RATE_ARRHENIUS
       do ipl = 1, n_plog_reactions
          rate_form(plog_reaction(ipl)) = RATE_PLOG
       end do

       if (n_plog_reactions > 0) then
!         PLOG forward rates and exact constant-volume analytic
!         derivatives are evaluated by reacpar::plog_kinf_eval.
          write(lout, "(' PLOG: ',I0,' pressure-dependent reaction(s),',"//&
                      "' ',I0,' pressure node(s) loaded and evaluable',"//  &
                      "' (cklink v2, analytic Jacobian enabled).')")        &
                      n_plog_reactions, n_plog_nodes
       endif

!      Fix third-body reactions sparse matrix rows



!      Count number of third-body reactions
       allocate(ire(nr))
       call sparse_internal_count(third_body_sp,ire,dim=2)
       where ( ire > 0 )
          THREE = .true.
       elsewhere
          THREE = .false.
       end where
       deallocate(ire)

!      Allocate THREE indices
       nTHREE = count(THREE)
       allocate(iTHREE(nTHREE))
       iTHREE = pack(reactions, THREE)


!      Alternative formulation for third-body matrix
!       third_body_beta = 1e0_dp - transpose(third_body)
!      Check if the formualtion adopted is equal to the following:
!       third_body_beta = 0e0_dp
       call allocate(nTHREE, ns, 0, tb_beta_sp   )
       do j = 1, nTHREE
          i = iTHREE(j)
!          where(third_body     (i,:)/= 1e0_dp) &
!                third_body_beta(:,i) = 1e0_dp - third_body(i,:)


          do k = 1, ns
             call       add_value   (tb_beta_sp,    j, k, &
                  one - sparse_value(third_body_sp, i, k) )

          end do

       end do

!      Fix number of columns
       tb_beta_sp%nc = ns

!      Logical index for simple third_body reactions

       nTB = count(THREE .and. (.not.(HIGH.or.LOW)))
       allocate(iTB(nTB))
       iTB = pack(reactions, THREE.and.(.not.(HIGH.or.LOW)))

!      Landau-teller reactions are not supported
       if (nlt  > 0) stop
       if (nrl  > 0) stop

!      Also radiation wavelength is not supported
       if (nw   > 0) stop

!     Close mechanism file and output transcript
      close(linc)
      close(lout)



!     ******************************************************************
!     Allocating types of reactions in module reacpar


      reac_cond = (.not.LOW.and..not.HIGH.and..not.TROE)
      allocate(Arrhreac(count(reac_cond)))
      Arrhreac = pack(reactions, reac_cond)

      reac_cond = (LOW.or.HIGH).and.(TROE)
      allocate(Troereac(count(reac_cond)))
      Troereac = pack(reactions, reac_cond)

      reac_cond = (LOW.or.HIGH).and.(.not.TROE)
      allocate(Lindreac(count(reac_cond)))
      Lindreac = pack(reactions, reac_cond)

      allocate(Revreac(count(REV)))
      Revreac  = pack(reactions, REV)


!     PRODUCTS ********************************************************
!     (data is only needed for reversible reactions)


!     jp, ip = row and column indexes (in stoich_p matrix) referring
!     to non-zero stoichiometric coefficients of species active in
!     reversible reactions

!     The arrays are ordered to speedup the calculation of species
!     concentration productories that contribute to the reaction
!     rate of progress variable

!      allocate(stoichprev(nr,ns))

      allocate(ip(stoich_p_sp%n),jp(stoich_p_sp%n),v_stoich_p(stoich_p_sp%n))
      call extract_rowcol_indices_columnwise(stoich_p_sp,ip,jp,v_stoich_p)

!     NB: not all the reactions are reversible! save time by
!         reducing the mass action productories arrays to the
!         reversible reactions only (index is ip)
      if (count(reversibile(ip)) < size(ip)) then
        allocate(tmpip(count(reversibile(ip))))
        allocate(tmpjp(count(reversibile(ip))))
        allocate(tmpv (count(reversibile(ip))))

        tmpip = pack(ip,         reversibile(ip))
        tmpjp = pack(jp,         reversibile(ip))
        tmpv  = pack(v_stoich_p, reversibile(ip))

        deallocate(ip,jp,v_stoich_p)
        allocate(ip(size(tmpip)), jp(size(tmpjp)), v_stoich_p(size(tmpip)))
        ip         = tmpip
        jp         = tmpjp
        v_stoich_p = tmpv

        deallocate(tmpip, tmpjp, tmpv)

      endif

!      stoichprev = stoich_p>tiny(0.0d0) .and. spread(reversibile, 2, ns)
!      call find_indices2D(stoichprev)
!      allocate(ip(size(i2D1)),jp(size(i2D2)))
!
!      ip = i2D1
!      jp = i2D2


!      deallocate(stoichprev)
      allocate(indice_p(nr,max_nprods))
      indice_p(:,:) = 0

      if (size(ip)>0) then
         do i=1,nr
           call find_indices(ip.eq.i)
           if (size(indices)>0)indice_p(i,1:size(indices)) = indices(:)
         end do
      endif

      call find_indices2D(indice_p > 0)
      allocate(i1p(nr,max_nprods))
      if (count(indice_p>0)>0)allocate(i2p(size(i2D1)))

      i1p = 0
      where (indice_p>0) i1p = 1


      if (allocated(i2p)) then

      i2p = pack(indice_p,indice_p>0)

!     Reorder v_stoich_p
      allocate(iv_stoich_p(size(v_stoich_p)),ijp(size(v_stoich_p)))

      do j = 1, size(v_stoich_p)
        iv_stoich_p(j) = int(v_stoich_p(i2p(j)))
        ijp(j) = jp(i2p(j))
      end do

      endif

      if (.not.allocated(i2D1p)) then
      call find_indices2D(i1p.eq.1)
      allocate(i2D1p(size(i2D1)))
      i2D1p = i2D1
      allocate(i2D2p(size(i2D2)))
      i2D2p = i2D2
      endif

!     REACTANTS *******************************************************


      allocate(ir(stoich_r_sp%n),jr(stoich_r_sp%n))
      allocate(v_stoich_r(size(ir)))
      call extract_rowcol_indices_columnwise(stoich_r_sp,ir,jr,v_stoich_r)


!      v_stoich_r = pack(stoich_r,stoich_r>0)

      allocate(indice_r(nr,max_nreacs))
      indice_r(:,:) = 0

      do i=1,stoich_r_sp%nr
         call find_indices(ir.eq.i)
         indice_r(i,1:size(indices)) = indices(:)
      end do

      call find_indices2D(indice_r > 0)
      allocate(i1r(nr,max_nreacs))
      allocate(i2r(size(i2D1)))

      i1r(:,:) = 0
      do i=1,size(indice_r,1)
         do j=1,size(indice_r,2)
          if (indice_r(i,j).gt.0) i1r(i,j) = 1
         end do
      end do

      i2r = pack(indice_r,indice_r>0)

!     Reorder v_stoich_r
      allocate(iv_stoich_r(size(v_stoich_r)),ijr(size(v_stoich_r)))

      do j = 1, size(v_stoich_r)
        iv_stoich_r(j) = int(v_stoich_r(i2r(j)))
        ijr(j) = jr(i2r(j))
      end do



      if (.not.allocated(i2D1r)) then
      call find_indices2D(i1r.eq.1)
      allocate(i2D1r(size(i2D1)))
      i2D1r = i2D1
      allocate(i2D2r(size(i2D2)))
      i2D2r = i2D2
      endif

      deallocate(i2D1, i2D2)

!     *****************************************************************
      end subroutine SCcklink
