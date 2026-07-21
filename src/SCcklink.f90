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

      use reacpar,      only: Arrhreac, Lindreac, Troereac, Revreac

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
         fmt_prec   = "(' Unsupported real data precision ')"

      character(len=80)                            :: ckfile
      character(len=16)                            :: ckvers, ckprec
      character(len=16), dimension(:), allocatable :: tmpchar

      logical :: ispresent, ispresent2, kerr
      integer :: i, j, k, idummy1, idummy2
      integer :: liwork, lrwork, lcwork, mm, kk, ii, maxsp, maxtb,    &
                 maxtp, nthcf, nipar, nitar, nifar, nrv, nfl, nthb,   &
                 nlt, nrl, nw, nchrg, max_nreacs, max_nprods
!ck2015 for real stoichometric coefficients
      integer :: nei,nja, njan, nf1, nif1, nex, nsto

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

      open(unit=linc, file = ckfile,  form   = 'unformatted', &
                                      status = 'unknown'      )
      rewind linc

!     ** PRELIMINARY PROBLEM DIMENSIONS AND DATA ***********************

!     Mechanism label
      mechanism = trim(adjustl("Mechanism imported from "//ckfile))

!     File header: Chemkin version, machine precision and check OK
      read (linc) ckvers, ckprec, kerr
      write(lout, fmt_info)ckvers,ckprec

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
               nlt, nrl, nw, nchrg, nei, nja, njan, nf1, nf1, nex, nsto
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
      assign_real_stoichiometric_coefficients: do i = 1, nsto
        k = idummy0(i)
        do j = 1, maxsp
!          A reactant is present
           if (rtmpst(j,k) < 0) then
               call add_value(stoich_r_sp,k,itmpsp(j,k),-rtmpst(j,k))
!          A product is present
           elseif (rtmpst(j,k) > 0) then
              call add_value(stoich_p_sp,k,itmpsp(j,k), rtmpst(j,k))
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