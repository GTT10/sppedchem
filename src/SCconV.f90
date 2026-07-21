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
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **          Constant-volume reactor chemistry module           **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: wednesday, 24/04/2012                        **
!     **                                                             **
!     *****************************************************************

      module speedchem_conV

      use working_precision

      implicit none
      private

!     ** Public data **************************************************
!      public :: conV_fun
!      public :: conV_jac

      public :: SC_conV
      public :: SC_conV_DAE
      public :: SC_conV_MEBDF
      public :: SC_conV_VODES

      public :: constV_jac_sparse
      public :: constV_jac_LSODES
      public :: constV_jac_MEBDF
      public :: constV_jac_VODES
      public :: constV_jac_DASPK
      public :: constV_jac_DASPK_sp
      public :: constV_jacv_ROWMAP
      public :: constV_jac
      public :: constV_jac_GAM
      public :: constV_jac_RADAUS
      public :: constV_jac_RADAU

      interface conV_fun
         module procedure SC_conV
         module procedure SC_conV_DAE
         module procedure SC_conV_MEBDF
      end interface conV_fun

      interface conV_jac
         module procedure constV_jac_sparse
         module procedure constV_jac_LSODES
         module procedure constV_jac_MEBDF
         module procedure constV_jac_VODES
!         module procedure constV_jac_RADAU
         module procedure constV_jac_DASPK
         module procedure constV_jac_DASPK_sp
         module procedure constV_jacv_ROWMAP
         module procedure constV_jac
         module procedure constV_jac_GAM
         module procedure constV_jac_RADAUS
      end interface conV_jac


!     *****************************************************************
      contains

!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **     Differential equation for constant-volume   reactor     **
!     **                                                             **
!     ** INPUT DATA                                                  **
!     ** yin   = array of the unknowns:                              **
!     **         y(1)           = T [K]                              **
!     **         y(2:nspecie+1) = Y1, Y2, ..., Y_n [-]               **
!     **                                                             **
!     ** OUTPUT DATA                                                 **
!     ** dydt  = unknowns rate of change [K/s],[mol/cm^3/s]          **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: wedesday, 03/12/2011                         **
!     **   Usage, copying, distribution of the proprietary routines  **
!     **   or any modified versions of them are not allowed without  **
!     **   approval of the author.                                   **
!     **                                                             **
!     *****************************************************************
         subroutine SC_conV(neq,time,yin,dyindt)

         use speedchem,       only: species_permutations,            &
                                    species_inverse_permutations
         use chemistry_setup, only: accurate_scthermo,               &
     &                              Temp_HIlim, Temp_LOlim,          &
     &                              permutate_species
         use SCthermodata,    only: temperature_array, table_indices,&
                                    use_table
         use SCmixturethermo, only: cvmas, SCrho, SCP, pressurerhoT
         use SCspeciesthermo, only: int_energy
         use kinetics_mod,    only: SCdwdt, reaction_rates
         use speedchem,       only: SCMW
         use universal_constants, only: zero, kilo
         use reacpar, only: uequilC
         use omp_lib

         implicit none

         integer,          intent(in)               :: neq
         real (dp)       , intent(in)               :: time
         real (dp)       , intent(in),     target   :: yin   (neq)
         real (dp)       , intent(out),    target   :: dyindt(neq)

         real (dp)       ,                 pointer  :: T, dTdt
         real (dp)       , dimension(:),   pointer  :: Y, dYdt

         logical,          dimension(neq-1)         :: nonzeroes
         real (dp)       , dimension(neq-1)         :: dwdt, Umol, Yp
         real (dp)       , dimension(6)             :: Ta
         real (dp)                                  :: cv, frac, urho

         integer                                    :: iT

         integer :: j, ir


!     *****************************************************************

!         Pointer associations
          if (permutate_species) then
             T => yin(neq)
             Y => yin(1:neq-1)

             ! revert to original order
             Yp = yin(species_inverse_permutations(2:neq))

             dTdt => dyindt(neq)
             dYdt => dyindt(1:neq-1)
          else
             T => yin(1)
             Y => yin(2:neq)

             Yp = Y

             dTdt => dyindt(1)
             dYdt => dyindt(2:neq)
          endif

!         Error handling for NaNs or negatives arising from ODE solver
          where ( abs(Yp) < small ) Yp = sign(small,Yp)
          if (.not. T > zero) return


!         NB: problem constant, density SCrho [kg/m3],
!             has to be set prior to calling this subroutine!

!         Compute temperatures declination
          Ta = temperature_array(T)

!         Compute system pressure
          SCP = pressurerhoT(T,Yp)

          if (.not. use_table(T)) then

!         Compute mixture constant volume specific heat [J/kg/K]
          cv = cvmas(Ta,Yp)

!         Compute species internal energies in moles [J/mol/K]
          Umol = int_energy(Ta)

!         Compute species production rates dwdt [mol/cm3/s] at
!         constant volume
          dwdt = SCdwdt(Ta,Yp)

          else

!         Array for reducing number of species computations
          call table_indices(T, iT, frac)
!          nonzeroes = ( Yp /= 0.e0_dp )

!         Compute mixture constant volume specific heat [J/kg/K]
          cv = cvmas(Ta,Yp,iT,frac)

!         Compute species internal energies in moles [J/mol/K]
          Umol = int_energy(Ta,iT,frac)

!         Compute species production rates dwdt [mol/cm3/s] at
!         constant volume
          dwdt = SCdwdt(Ta,Yp,iT,frac)

          endif

!         Compute reciprocal density [cm3/kg]
          urho = 1.d3/SCrho

!         Computing mass fraction rate of change imposing mass consv.
          dYdt = urho * dwdt * SCMW

!         Computing energy conservation in the system [K/s]
          dTdt = - kilo * sum(Umol * dwdt) * urho / cv

          if (permutate_species) dYdt = dYdt(species_permutations(1:neq-1)-1)

         end subroutine SC_conV
!     *****************************************************************

!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **               Differential-Algebraic Equation               **
!     **        function for constant-volume reactor problem         **
!     **                                                             **
!     **      G(t, y, yprime) = yprime - f(t, y) = residual (=0)     **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: thursday, 24/05/2012                         **
!     **   Usage, copying, distribution of the proprietary routines  **
!     **   or any modified versions of them are not allowed without  **
!     **   approval of the author.                                   **
!     **                                                             **
!     *****************************************************************
         subroutine SC_conV_DAE(time,y,yprime,CJ,delta,ires)

         use speedchem, only: neq

         implicit none

         real (dp)                           , intent(inout) :: time
         real (dp)       , dimension(neq)    , intent(inout) :: y
         real (dp)       , dimension(neq)    , intent(inout) :: yprime
         real (dp)       ,                     intent(inout) :: cj
         real (dp)       , dimension(neq)    , intent(out)   :: delta
         integer,                              intent(inout) :: ires
         integer                              :: j
         real (dp)       , dimension(size(y)) :: dydt


!     *****************************************************************

!           Call derivative routine
            call SC_conV(neq,time,y,dydt)

!           Compute DAE residual
            delta = yprime - dydt

!           Successful step: ires = 0
            ires = 0

         end subroutine SC_conV_DAE
!     *****************************************************************

         subroutine SC_conV_VODES(NEQ,T,Y,YDOT)

         implicit none

         integer, parameter :: WP = KIND(1.0D0)
         integer                             , intent(in ) :: neq
         real (wp)                           , intent(in ) :: T
         real (wp)       , dimension(neq)    , intent(in ) :: Y
         real (wp)       , dimension(neq)    , intent(out) :: YDOT
         real (dp)       , dimension(neq)                  :: YDOT_DP

!     *****************************************************************
!        Call derivative routine
         call SC_conV(NEQ,real(T,dp),real(Y,dp),YDOT_DP)
         YDOT = dble(YDOT_DP)

         end subroutine SC_conV_VODES
!     *****************************************************************



!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **                         ODE system                          **
!     **        function for constant-volume reactor problem         **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: thursday, 24/05/2012                         **
!     **   Usage, copying, distribution of the proprietary routines  **
!     **   or any modified versions of them are not allowed without  **
!     **   approval of the author.                                   **
!     **                                                             **
!     *****************************************************************
         subroutine SC_conV_MEBDF(neq,time,y,dydt,ipar,rpar,ierr)

         implicit none

         integer                             , intent(in   ) :: neq
         real (dp)                           , intent(in   ) :: time
         real (dp)       , dimension(neq)    , intent(in   ) :: y
         real (dp)       , dimension(neq)    , intent(out  ) :: dydt
         real (dp)       , dimension(*)      , intent(inout) :: rpar
         integer         , dimension(*)      , intent(inout) :: ipar
         integer,                              intent(out  ) :: ierr
         integer                              :: j
!     *****************************************************************
!           Call derivative routine
            call SC_conV(neq,time,y,dydt)

!           Successful step: ierr = 0
            ierr = 0

         end subroutine SC_conV_MEBDF
!     *****************************************************************


!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **   Computing chemistry jacobian for constant-volume reactor  **
!     **                                                             **
!     **   constV_jac = [df1/dT  df1/dY1   .... df1/dY_n ]           **
!     **                | ...                            |           **
!     **                [dfm/dT  dfm/dY1   .... dfm/dY_n ]           **
!     **                                                             **
!     **                                                             **
!     ** INPUT DATA                                                  **
!     ** yin   = array of the unknowns:                              **
!     **         y(1)           = T [K]                              **
!     **         y(2:nspecie+1) = Y1, Y2, ..., Y_n [-]               **
!     **                                                             **
!     ** OUTPUT DATA                                                 **
!     ** Jacobian matrix in sparse form (JAC_sparse)                 **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: wedesday, 19/04/2012                         **
!     **   Usage, copying, distribution of the proprietary routines  **
!     **   or any modified versions of them are not allowed without  **
!     **   approval of the author.                                   **
!     **                                                             **
!     *****************************************************************

      subroutine constV_jac_sparse(neq,time,yin)
      use speedchem,       only : nr,ns,nel,SCMW,HIGH,LOW,nudiff,TROE,    &
                                  A0,Ainf,b0,binf,E0,Einf,iTHREE,iTB,     &
                                  sumnudiff,inotrev,nnotrev,              &
                                  ARi,bRi,ERi,REV,uMW, reversibile,       &
                                  nTHREE,nTB, sparse_jac, njac,           &
                                  stoich_r_pack,                          &
                                  stoich_p_pack, i_stoich_r, i_stoich_p,  &
                                  n_stoich_r, n_stoich_p, n_nudiff,       &
                                  i_nudiff, nudiff_pack, species_permutations,&
                                  species_inverse_permutations

      use chemistry_setup, only: accurate_scthermo, simplified_for_sparsity,&
                                  Temp_table_accuracy, Temp_LOlim,        &
                                  Temp_HIlim, tab_nsteps,                 &
                                  rec_Ttable_accuracy,                    &
                                  save_thermal_parameters, solver,        &
                                  permutate_species

      use SCthermodata, only : aL, bL, cL, dL, eL, fL, gL, tsw,           &
                               aH, bH, cH, dH, eH, fH, gH,                &
                               tab_CpsuR, tab_HsuRT, tab_SsuR,            &
                               tab_uKc, tab_dKcdT, use_table,             &
                               tab_dGdT, tab_dCvdT, interp4_coefs,        &
                               temperature_array, table_indices
      use troepar,      only : todotroe, aT2, uT1T2, T2T2, uT3T2,         &
                               zeroT2, tab_log10Fcent, itbALL, ntbALL,    &
                               iTBTROE, itbSIMP, itbFALL, itbLIND,        &
                               iTROEiTBALL, iSIMPitbALL, iFALLitbALL,     &
                               iLINDitbALL, ntbFALL, ntbLIND, ntbSIMP,    &
                               ntbTROE, iTROE4, troe_logfac, iTROEitbFALL,&
                               iLINDitbFALL
      use kinetics_mod, only : ijr, i2r, iv_stoich_r, i2D1r,              &
                               ijp, i2p, iv_stoich_p, i2D1p,              &
                               tab_k0,  tab_kinf,   iA0,                  &
                               tab_dk0dT, tab_dkinfdT,                    &
                               reaction_rates_and_derivative,             &
                               explicit_rev_reaction_rates_deriv,         &
                               reaction_rates_derivatives,                &
                               save_k0, save_kinf, lsavek,                &
                               mass_action_productories
      use reacpar,      only : Troereac, nEQREV, iEQREV, nXREV, iXREV,    &
                               uequilC_and_derivative, nTREV, iTREV,      &
                               is_beta_pack, tb_beta_pack, n_tb_beta
      use SCspeciesthermo, only : CvuRmol, int_energy, dCv_dT
      use SCmixturethermo, only : SCP,SCrho,cp,cv, cvmas, molar_volumes
      use sparse_chemistry, only: sparse, nudiffT_sparse,                 &
                                  dq_dY_sparse, sparse_allocate,          &
                                  sparse_to_dense,JACYY_sparse,           &
                                  nudiffT_molarv_sparse, nudiffT_sparse,  &
                                  JACYYT_sparse, JAC_sparse, JACT_sparse, &
                                  sparse_transpose, sparse_symbolic_mm,   &
                                  dq_dY_T_sparse,stoich_r_sp,stoich_p_sp, &
                                  extract_sparse_index,                   &
                                  sparse_2_matmul,identity,               &
                                  sparse_transpose_valuesonly,            &
                                  tb_beta_sp, stoich_p_eff_sp,            &
                                  stoich_r_eff_sp
      use sparse_algebra, only: dense_to_sparse, sparse_row_prod,         &
                                sparse_value, sparse_nullify_general,     &
                                sparse_col_prod, sparse_partial_sum,      &
                                sparse_col_prod_valonly,                  &
                                sparse_row_prod_valonly,sparse_sum,       &
                                print_sparse_to_file, &
                                sparse_compress, sparse_sort,             &
                                sparse_matmulT, sparse_square_permutation,&
                                add_value, add_line

      use sparse_definitions

      use universal_constants, only: Patm, uPatm, R, uR, Rcal, uRcal,     &
                                     ln10, uln10, u2, u3, u6, one, two,   &
                                     ten, zero, milli, kilo, mega

      implicit none

!     ** Data input-output ********************************************
      integer,                          intent(in)         :: neq
      real (dp)       ,                 intent(in)         :: time
      real (dp)       , dimension(neq), intent(in), target :: yin

      real (dp)       , dimension(:,:), allocatable        :: JAC

!     *****************************************************************

      real (dp)       ,               pointer :: T
      real (dp)       , dimension(:), pointer :: Y

      real (dp)       , dimension(:,:), allocatable :: JACYY
      real (dp)       , dimension(ns)    :: JACTY

      real (dp)       , dimension(nr)    :: qf, qb

      real (dp)       , dimension(6), target :: Ta
      real (dp)        :: dTdt, uYj

      real (dp)       , dimension(neq-1) :: X,C,dwdt,dYdt,molar_volume


!     SCthermo variables
      real (dp)                       :: Ctot, ucv
      real (dp)       , dimension(ns) :: CpsuR,CpsuRH,HsuRT,SsuR

      real (dp)       , dimension(ns) :: Cpmol,Cvmol,Hmol,Umol,           &
                                                Smol,cpm,cvm,hm

!     SCthermomix variables
      real (dp)        :: MWm, Rm, h, e

!     Production rate variables
      real (dp)       , dimension(nXREV) :: kbXrev, dkbXrevdT
      real (dp)       , dimension(nEQREV,2), target  :: uKc_dKcdT
      real (dp)       , dimension(:),        pointer :: uKc, dKcdT
      real (dp)       , dimension(nr) :: kf
      real (dp)       , dimension(nr) :: dG0, uKp, kb
      real (dp)       , dimension(nr) :: q, qeff, qefff, qeffb
      real (dp)        :: faceq, ufaceq
      integer :: iel , ii, iii, iT, ivsp, ivsr, j0, jf

!     Request for accurate estimation of thermodynamic properties
      logical :: compute_accurate_thermo

      logical, parameter :: parabolic = .true.

!     deltaG0
      real (dp)       , dimension(ns) :: gibbs

!     k_forward
      real (dp)       , dimension(nr) :: kinf



!     TREE BODY REACTIONS
      real (dp)       , dimension(ntbALL ) :: M,uFTLALL,           &
                                              Pr_dlogFT_dlodPr_uMi,&
                                              tenlogFall,          &
                                              prod_rate_consts
      real (dp)       , dimension(ntbFALL) :: Pr, ukinfpmk0
      real (dp)       , dimension(ntbTROE) :: Pr2,Fcent,           &
                                              log10Fcent,          &
                                              log10Pr,             &
                                              ctroe,ntroe,logF,    &
                                              k_uni,fattore,       &
                                              ten_pow_logF

      real (dp)       , parameter :: dtroe = 0.14e0_dp

!     constV_energy
      real (dp)        :: urho, coefs(5)

      integer :: j,i,k,isp,ire,idqdY
      real (dp)        :: frac,tenthT

!     PART 1 - DERIVATIVES WITH RESPECT TO TEMPERATURE ****************
      real (dp)                              :: dT_dT
      real (dp)       , dimension(nr)        :: dkinfdT, dkbdT
      real (dp)       , dimension(ntbFALL)   :: dk0dTs, k0s
      real (dp)       , dimension(ntbFALL)   :: dk0dT, k0, k0M
      real (dp)       , dimension(nr)        :: dq_dT
      real (dp)       , dimension(nr)        :: prod_f, prod_b

      real (dp)       , dimension(neq-1)     :: dCvdT
      real (dp)       , dimension(ns)        :: dY_dT

      real (dp)       , dimension(ntbFALL)   :: FTL, uFTL, dFTL_dT
      real (dp)       , dimension(ntbTROE)   :: dlogFtroe_dlogPr, &
                                                dPr_dT, dlogP_dT, &
                                                logPpc

!     PART 2 - SPECIES DERIVATIVES WITH RESPECT TO SPECIES ************

      real (dp)       , dimension(neq-1)     :: uY, dCtot_dY

!     Factor for Three-body reactions
      real (dp)       , dimension(:,:), allocatable :: dMeff_dY
      real (dp)                                     :: C_pow_nu
      real (dp)       , dimension(neq-1)            :: dlogP_dY
      real (dp)       , dimension(neq-1), target    :: YY

!     PART 3 - TEMPERATURE DERIVATIVES WITH RESPECT TO SPECIES ********
      real (dp)       , dimension(neq-1)         :: JACYYT_UuMW
      type(sparse)                               :: dMeff_dY_sp,      &
                                                    dFTL_dY_sp,       &
                                                    tmp_sp

!     *****************************************************************

!     ** Preliminary data assignments *********************************

!     Dimension checks
      if (sparse_jac.and.(.not.allocated(JAC_sparse%A))) then
         write(*,*)'Error in sparse jacobian dimension, constV_jac: '
         write(*,*)'n_jac = ',JAC_sparse%n

         stop
      endif

!     Security check for temporarily empty/wrong items
      if (any(yin == yin + one)) then
         JAC_sparse%A = small
         return
      endif

!     Assigning temperature and species mass fractions

      if (permutate_species) then
         YY =  yin(species_inverse_permutations(2:neq))
         T  => yin(neq)
      else
         YY =  yin(2:neq)
         T  => yin(1)
      endif

!     Guarantee that all the species concentrations are non zero!
!     Needed for not having discontinuities in derivative formulations
      Y  => YY
      where ( abs(Y) < small ) Y = sign(small,Y)
      uY = one/Y

!     Error check on temperature (.not.T>0.e0_dp) means T is negative
!     or NaN; the same for pressure
      if (.not. T   > zero) return
      if (.not. SCP > zero) return

!     Flag to compute using accurate/tabulated thermodynamic data
      compute_accurate_thermo = .not.use_table(T)

!     Assign Jacobian sub-parts: ÚÄÄÄÄÂÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
!                                ³dTdT³        JACTY          ³
!                                ÃÄÄÄÄÅÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ´
!                                ³    ³                       ³
!                                ³ d  ³                       ³
!                                ³ Y  ³                       ³
!                                ³ _  ³                       ³
!                                ³ d  ³        JACYY          ³
!                                ³ T  ³                       ³
!                                ³    ³                       ³
!                                ³    ³                       ³
!                                ³    ³                       ³
!                                ÀÄÄÄÄÁÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ

!     ** Load thermodynamic temperature-dependent parameters **************
!     Temperatures and beyond
      Ta     = temperature_array(T)

!     Tabulated / computed NASA thermodynamic properties
      if (compute_accurate_thermo) then

         dCvdT      = dCv_dT    (Ta)
         Umol       = int_energy(Ta)
         if (ntbTROE>0) log10Fcent = troe_logfac(Ta)
         call CvuRmol(Ta,Cvmol)
         call reaction_rates_and_derivative(Ta, k0, dk0dT, kinf, dkinfdT)

      else

!        Compute temperature-dependent tabulation parameters
         call table_indices(T, iT, frac)
         coefs = interp4_coefs(frac)

         dCvdT      = dCv_dT    (Ta,iT,frac)
         Umol       = int_energy(Ta,iT,frac)
         if (ntbTROE>0) log10Fcent = troe_logfac(Ta,iT,frac)
         call CvuRmol(Ta,Cvmol,iT,frac)
         call reaction_rates_and_derivative(Ta, k0, dk0dT, kinf, dkinfdT, iT, frac)

      endif

!     Specific heat [J/mol/K]
      Cvmol = Cvmol * R

!     ** Average mixture properties ************************************
!     constant volume specific heats [J/kg/K]
      cv  = kilo * sum(Y * Cvmol * uMW)
      ucv = one / cv

!     Average mixture concentration Ctot [mol/cm3] and derivative
!     Species concentrations [mol/cm3]
      dCtot_dY = milli * SCrho * uMW
      C        = dCtot_dY * Y
      Ctot     = sum(C)

!     *****************************************************************
!                      PRODUCTION RATE [mol/cm3/s]
!     *****************************************************************

!     Computing forward reaction rate constants
      kf = kinf

!     Computing effective molecularity of the reactions
      if (ntbALL > 0) M = Ctot - (tb_beta_sp*C)

!     Computing reaction reduced pressure value for falloff reactions
      if (ntbFALL > 0) then
         k0M        = k0 * M(iFALLitbALL)
         ukinfpmk0  = one / (kinf(itbFALL) + k0M)
         Pr         = k0M * ukinfpmk0
      endif

!     Computing reaction rate constants according to Lindemann's form
      if (ntbLIND > 0) FTL(iLINDitbFALL) = Pr(iLINDitbFALL)

!     Computing reaction rate constants according to Troe's form
      troefactors: if (ntbTROE > 0) then

        Pr2 = Pr(iTROEiTBFALL)

!       ** Computing troe centering factor
         log10Pr    = log10(Pr2/(one-Pr2))

!        ** Troe model parameters
         ctroe = -0.40_dp - 0.67_dp * log10Fcent
         ntroe =  0.75_dp - 1.27_dp * log10Fcent

         fattore = ((log10Pr+ctroe)/(ntroe-dtroe*(log10Pr+ctroe)))
         logF = log10Fcent/(one + fattore*fattore)

         ten_pow_logF = ten ** logF

         FTL(iTROEitbFALL) = Pr2 * ten_pow_logF

      endif troefactors

      if (ntbFALL > 0) then
         kf(itbFALL) = kf(itbFALL) * FTL
!        Reciprocal of pressure-dependent reaction factor
         uFTL = one/FTL
      endif

!     *****************************************************************

!     Computing backward reaction rate constants based on the
!     equilibrium assumption
      if (nnotrev > 0) then
         kb   (inotrev) = zero
         dkbdT(inotrev) = zero
      endif

!     ** Computing EQUILIBRIUM CONSTANTS ******************************
!        and derivative with respect to temperature,
!        for reversible reactions with equilibrium-driven backw rate
      if (nEQREV > 0) then

        if (compute_accurate_thermo) then
           uKc_dKcdT = uequilC_and_derivative(Ta)
        else
           uKc_dKcdT = uequilC_and_derivative(Ta,iT,frac)
        endif

        uKc   => uKc_dKcdT(:,1)
        dKcdT => uKc_dKcdT(:,2)

        kb   (iEQREV) = uKc * kf(iEQREV)
        dkbdT(iEQREV) = uKc * (dkinfdT(iEQREV) - kinf(iEQREV)*uKc*dKcdT)
      endif

!     Explicit reverse reaction rates in mechanism input
      if (nXREV > 0) then

        if (compute_accurate_thermo) then
            call explicit_rev_reaction_rates_deriv &
                              (Ta, kbXrev, dkbXrevdT)
        else

            call explicit_rev_reaction_rates_deriv &
                              (Ta, kbXrev, dkbXrevdT, iT, frac)

        endif

        kb(iXREV)    = kbXrev
        dkbdT(iXREV) = dkbXrevdT

      endif

!     Computing reaction progress variable (KINETICS) *****************
      call mass_action_productories(C,prod_f,prod_b)

!     Final forward/backward progress variables
      qf = prod_f*kf
      qb = prod_b*kb

      qefff = qf
      qeffb = qb
      if (ntbSIMP > 0) then
         qefff(itbSIMP) = qf(itbSIMP) * M(iSIMPitbALL)
         qeffb(itbSIMP) = qb(itbSIMP) * M(iSIMPitbALL)
      endif

!     Total rates of change of progress variable
      q    = qf    - qb

!     Effective progress variables in presence of three body reactions
      qeff = qefff - qeffb

!     *****************************************************************

!     Sparse multiplication
      dwdt = nudiffT_sparse*qeff

!     Reciprocal density [m3/kg]
      urho = one/SCrho

!     ** END OF PRODUCTION RATE ***************************************

!     Computing mass fraction rate of change imposing mass consv. [1/s]
      molar_volume = kilo * SCMW * urho    ! [cm3/mol]
      dYdt = dwdt * molar_volume

!     Computing energy conservation in the system [K/s]
      dTdt = - mega * sum(Umol*dwdt) * urho * ucv

!     *****************************************************************
!     *****************************************************************
!            S T A R T   J A C O B I A N   A S S I G N M E N T
!     *****************************************************************
!     *****************************************************************

!     ** PART 1 - DERIVATIVES WITH RESPECT TO TEMPERATURE *************

!     ** Derivatives of reaction progress variables with respect to
!     ** temperature [mol/cm3/s/K]

!     Simple, arrhenius reaction rates involved
      dq_dT = prod_f * dkinfdT - prod_b * dkbdT

!     Simple third-body reactions have dMeff_dT = 0 always
      if (ntbSIMP > 0) dq_dT(itbSIMP) = dq_dT(itbSIMP) * M(iSIMPitbALL)

!     Updating rates for LINDEMANN reactions (also needed to build TROE parameters)
      if (ntbFALL > 0) dFTL_dT = M(iFALLitbALL) * ukinfpmk0**2 * &
                               ( kinf(itbFALL) * dk0dT - dkinfdt(itbFALL) * k0 )


!     Updating rates for TROE reactions (array operations)
      TROE_dFTLdT: if (ntbTROE > 0) then

         logPpc           = log10Pr + ctroe
         dlogFtroe_dlogPr = - two * ntroe * log10Fcent * logPpc * (ntroe - dtroe*logPpc) &
                          / ( ntroe**2 - two * ntroe * dtroe * logPpc                    &
                          + (dtroe**2 + one) * logPpc**2 ) **2
         dlogP_dT         = (dk0dT  (iTROEitbFALL)/k0  (iTROEitbFALL)                    &
                            -dkinfdT(itbTROE     )/kinf(itbTROE     )) * uln10
         dFTL_dT(iTROEitbFALL) = ten_pow_logF * (   dFTL_dT(iTROEitbFALL)                &
                                            + Pr2 * ln10 * dlogFtroe_dlogPr * dlogP_dT )

      endif TROE_dFTLdT

!     ** Derivative of rates of progress variable w/ respect to temperature
!        for falloff reactions
!       NB Pressure-dependent reactions not surely have equilibrium-based computatin
!          of backward reaction rate, but formula is OK because
!          in case of explicit reverse reaction rate
!          kb has already been included the TROE/lindemann correction factor FTL
      if (ntbFALL > 0)                                               &
         dq_dT(itbFALL) = dq_dT(itbFALL) * FTL                       &
                        + dFTL_dT *                                  &
                        (  kinf(itbFALL) * prod_f(itbFALL)           &
                          -kb  (itbFALL) * prod_b(itbFALL) * uFTL )


!     ** Derivative of mass fraction variations with respect to temperature [1/s/K]
      if (.not.sparse_jac) call molar_volumes
      dY_dT = nudiffT_molarv_sparse * dq_dT
!     Set a safety-related very small value to zeroes
      where (abs(dY_dT) < small) dY_dT = sign(small, dY_dT)

!     ** Derivative of temperature variation dT/dt with respect to
!     ** temperature [1/s]
      dT_dT = - ucv * kilo * (dTdt * sum(Y*uMW*dCvdT) +              &
                              sum(uMW * (Cvmol*dYdt + Umol*dY_dT) ) )

!     ** PART 2 - SPECIES DERIVATIVES WITH RESPECT TO SPECIES *********

!     Compute dMeff/dY [mol/cm3]
!     If not simplified for sparsity, this matrix is full; otherwise,
!     it is usually very sparse.
      if (.not.simplified_for_sparsity) then

         call allocate(ntbALL, ns, ns*nTBALL, dMeff_dY_sp)

!        Assign constant term
         dMeff_dY_sp%IA(1)      = 1
         do j = 1, nTBALL
            dMeff_dY_sp%IA(j+1)             = 1 + j*ns
            dMeff_dY_sp%JA(1+(j-1)*ns:j*ns) = [(i, i=1, ns)]
            do i = 1+(j-1)*ns, j*ns
               ii = i - (j-1)*ns
               dMeff_dY_sp%A (i) = dCtot_dY(ii) * (one - sparse_value(tb_beta_sp,j,ii))
            end do
         end do

      else

         dMeff_dY_sp = sparse_row_prod(tb_beta_sp, -dCtot_dY)

      endif


!     ** Derivatives of falloff reactions enhancement factor with respect
!        to species mass fractions, dFTL_dY(ntbALL, ns)

!        NB: dFTL_dY matrix = dMeff_dY at simple third-body reactions!
         if (ntbALL > 0) then

!        A) Lindemann reactions (also brick for TROE ones)
         prod_rate_consts(iSIMPitbALL) = zero
         prod_rate_consts(iFALLitbALL) = k0*ukinfpmk0**2*kinf(itbFALL)
         dFTL_dY_sp                    = sparse_col_prod(dMeff_dY_sp, prod_rate_consts)

!        B) Troe reactions
         Pr_dlogFT_dlodPr_uMi              = one
         Pr_dlogFT_dlodPr_uMi(iTROEitbALL) = Pr_dlogFT_dlodPr_uMi(iTROEitbALL) &
                                           * Pr2 * dlogFtroe_dlogPr /M(iTROEitbALL)

         tmp_sp = dFTL_dY_sp + sparse_col_prod(dMeff_dY_sp, Pr_dlogFT_dlodPr_uMi)
         dFTL_dY_sp = tmp_sp

         tenlogFall              = q(itbALL)
         tenlogFall(iTROEitbALL) = tenlogFall(iTROEitbALL) * ten_pow_logF
         tenlogFall(iFALLitbALL) = tenlogFall(iFALLitbALL) * uFTL
         call sparse_col_prod_valonly(dFTL_dY_sp, tenlogFall)

         endif


!     ** Derivatives of rates of progress variable with respect to the species
!        mass fractions, dq_dY(nr, ns)
      dq_dY_sparse = sparse_col_prod(stoich_r_eff_sp, qefff)

      if (nnotrev < nr) dq_dY_sparse = dq_dY_sparse - sparse_col_prod(stoich_p_eff_sp, qeffb)

      call sparse_row_prod_valonly(dq_dY_sparse, uY)

      tmp_sp = sparse_partial_sum(dq_dY_sparse,dFTL_dY_sp,itbALL)
      dq_dY_sparse = tmp_sp

!     RETRIEVING SPARSE FORMULATION FOR JACYY_SPARSE MATRIX
      if (.not.sparse_jac) call sparse_symbolic_mm(nudiffT_molarv_sparse, dq_dY_sparse, JACYY_sparse)
      call sparse_2_matmul(nudiffT_molarv_sparse, dq_dY_sparse, JACYY_sparse)

!     ** PART 3 - TEMPERATURE DERIVATIVES WITH RESPECT TO SPECIES *****
!     JACTY = d(dT/dt)/dY [K/s]
      JACYYT_UuMW = sparse_matmulT(JACYY_sparse,Umol*uMW)
      JACTY = - kilo * ucv * (dTdt * Cvmol * uMW + JACYYT_UuMW)
!     ** Set a safety-related very small value to zeroes
      where (abs(JACTY) < small) JACTY = sign(small, JACTY)



!     Compute sparse jacobian structure, if needed
      if (.not.sparse_jac) then

         call add_line (JAC_sparse,1,[dT_dT,JACTY])
         do j = 2, neq
            call add_value(JAC_sparse,j,1,dY_dT(j-1))
            do i = 2, neq
               call add_value(JAC_sparse,i,j,sparse_value(JACYY_sparse,i-1,j-1))
            end do
         end do

         if (permutate_species) JAC_sparse = sparse_square_permutation(JAC_sparse,species_permutations)

      end if


!     ASSEMBLING JACOBIAN IN SPARSE FORM
!     When the (permutate_species) option is selected, temperature is
!     the last of the unknowns; when not, it is the first.
      jacobian_assemble: if (.not.permutate_species) then

!     Designating first row of sparse jacobian matrix
      JAC_sparse%A(1)     = dT_dT
      JAC_sparse%A(2:neq) = JACTY

!     Designating species rows of sparse jacobian matrix
      do j = 1, ns

!       1st item in the row is temperature derivative
        i = JAC_sparse%IA(j+1)

        if (i/=JAC_sparse%IA(j+2)) then

        JAC_sparse%A (i) = dY_dT(j)
        JAC_sparse%JA(i) = 1

!       Other items are the species derivatives
        JAC_sparse  %A(i+1               :JAC_sparse%IA  (j+2)-1) = &
        JACYY_sparse%A(JACYY_sparse%IA(j):JACYY_sparse%IA(j+1)-1)

!       Need to update columns, as symbolic matrix matrix multiply routine
!       does change the column ordering in each row for faster computation
        JAC_sparse  %JA(i+1               :JAC_sparse%IA  (j+2)-1) = &
        JACYY_sparse%JA(JACYY_sparse%IA(j):JACYY_sparse%IA(j+1)-1) + 1

        endif

      end do

      else ! jacobian_assemble


!        Designating species rows of sparse jacobian matrix
         do j = 1, ns

!          Current row in permutated order
           ii = species_permutations(j)-1

!          1st item in the row is temperature derivative
           i = JAC_sparse%IA(j)

           if (i/=JAC_sparse%IA(j+1)) then

!          Other items are the species derivatives
           JAC_sparse  %A(i                 :JAC_sparse%IA  (j+1)-2) = &
           JACYY_sparse%A(JACYY_sparse%IA(ii):JACYY_sparse%IA(ii+1)-1)

!          Need to update columns, as symbolic matrix matrix multiply routine
!          does change the column ordering in each row for faster computation
           JAC_sparse%JA(i                   :JAC_sparse%IA  (j+1)-2) =   &
           species_inverse_permutations(1+JACYY_sparse%JA(JACYY_sparse%IA(ii):JACYY_sparse%IA(ii+1)-1))

!          Last column in the row is the temperature derivative of
!          current species
           JAC_sparse%A(JAC_sparse%IA(j+1)-1) = dY_dT(species_permutations(j)-1)
           JAC_sparse%JA(JAC_sparse%IA(j+1)-1) = neq

           endif

         end do


!        Designating last row of sparse jacobian matrix
         i = JAC_sparse%IA(ns+1)
         JAC_sparse%A(i:i+ns-1) = JACTY(species_permutations(1:ns)-1)
         JAC_sparse%A(i+ns)     = dT_dT
         JAC_sparse%JA(i:i+ns) = [(j,j=1,neq)]


      endif jacobian_assemble


      where(abs(JAC_sparse%A(1:JAC_sparse%n)) < small) &
      JAC_sparse%A(1:JAC_sparse%n) = sign(small,JAC_sparse%A(1:JAC_sparse%n))

!     Clear sparse memory storage allocation
      call sparse_nullify_general(dMeff_dY_sp)
      call sparse_nullify_general(dFTL_dY_sp )
      call sparse_nullify_general(tmp_sp     )

 900  format(' error in thermodynamic data calculation: T = ',F7.2,' K')
      end subroutine constV_jac_sparse








!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **   Interface for chemistry Jacobian computation by LSODES    **
!     **                    stiff ODE integrator                     **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: sunday, 20/11/2011                           **
!     **   Usage, copying, distribution of the proprietary routines  **
!     **   or any modified versions of them are not allowed without  **
!     **   approval of the author.                                   **
!     **                                                             **
!     *****************************************************************

      subroutine constV_jac_LSODES(neq,time,yin,J,IAN,JAN,PDJ)

      use speedchem,        only: species_permutations
      use sparse_chemistry, only: JAC_sparse, JACT_sparse
      use sparse_algebra,   only: sparse_transpose,            &
                                  sparse_transpose_valuesonly, &
                                  sparse, sparse_square_permutation, &
                                  sparse_transpose_permutations

      implicit none

!     ** Data input-output ********************************************
      integer,                          intent(in)         :: neq
      real (dp)       ,                 intent(in)         :: time
      real (dp)       , dimension(neq), intent(inout)      :: yin
      integer,                          intent(in)         :: J

      integer, dimension(:),            intent(inout)      :: IAN, JAN

      real (dp)       , dimension(neq), intent(out)        :: PDJ

      integer :: i, ii

      if (J == 1) then
!        Interface to transposed jacobian call
         call constV_jac_sparse(neq,time,yin)

!        Store Jacobian transposed matrix (faster conversion to dense format,
!        if items are filled columnwise!)

         if (JACT_sparse%n /= JAC_sparse%n) then
           JACT_sparse = sparse_transpose(JAC_sparse)
         else
           call sparse_transpose_valuesonly(JAC_sparse,JACT_sparse%A)
         endif

!          JACT_sparse = sparse_square_permutation(JACT_sparse,[1, species_permutations+1])

!         [fp]: no need to set zeroes as the LSODES solver already does it
!               and only takes values that match the current sparsity pattern
!         PDJ(1:neq) = 0.e0_dp

!      else
!         [fp]: no need to set zeroes as the LSODES solver already does it
!               and only takes values that match the current sparsity pattern
!         do i = JACT_sparse%IA(J-1),JACT_sparse%IA(J)-1
!            PDJ(JACT_sparse%JA(i)) = 0.e0_dp
!         end do

      endif

      do i = JACT_sparse%IA(J),JACT_sparse%IA(J+1)-1
         PDJ(JACT_sparse%JA(i)) = JACT_sparse%A(i)
      end do


      end subroutine constV_jac_LSODES


!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **   Interface for chemistry Jacobian computation by MEBDF     **
!     **                    stiff ODE integrator                     **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: sunday, 05/08/2012                           **
!     **   Usage, copying, distribution of the proprietary routines  **
!     **   or any modified versions of them are not allowed without  **
!     **   approval of the author.                                   **
!     **                                                             **
!     *****************************************************************

      subroutine constV_jac_MEBDF(neq,time,yin,J,PDJ,IPAR,RPAR,IERR)

!      use speedchem, only: LSODES_JAC
      use speedchem,        only: species_permutations
      use sparse_chemistry, only: JAC_sparse, JACT_sparse
      use sparse_algebra,   only: sparse_transpose,                  &
                                  sparse_transpose_valuesonly,       &
                                  sparse, sparse_square_permutation, &
                                  sparse_transpose_permutations,     &
                                  identity
      use sparse_definitions

      implicit none

!     ** Data input-output ********************************************
      integer,                          intent(in)         :: neq
      real (dp)       ,                 intent(in)         :: time
      real (dp)       , dimension(neq), intent(inout)      :: yin
      integer,                          intent(in)         :: J

      integer,          dimension(:),   intent(inout)      :: IPAR
      real (dp)       , dimension(:),   intent(inout)      :: RPAR
      integer,                          intent(out)        :: IERR

      real (dp)       , dimension(neq), intent(out)        :: PDJ

      type(sparse)                                         :: JAC_MEBDF

      integer :: i, ii

      if (J == 1) then
!        Interface to transposed jacobian call
         call constV_jac_sparse(neq,time,yin)

         JAC_MEBDF = JAC_sparse + identity(neq,tiny(0.e0_dp))

!        Store Jacobian transposed matrix (faster conversion to dense
!        format, if items are filled columnwise!)

         if (JACT_sparse%n /= JAC_MEBDF%n) then
           JACT_sparse = sparse_transpose(JAC_MEBDF)
         else
           call sparse_transpose_valuesonly(JAC_MEBDF,JACT_sparse%A)
         endif

         call deallocate(JAC_MEBDF)

!         [fp]: no need to set zeroes as the LSODES solver already does it
!               and only takes values that match the current sparsity pattern
!         PDJ(1:neq) = 0.e0_dp

!      else
!         [fp]: no need to set zeroes as the LSODES solver already does it
!               and only takes values that match the current sparsity pattern
!         do i = JACT_sparse%IA(J-1),JACT_sparse%IA(J)-1
!            PDJ(JACT_sparse%JA(i)) = 0.e0_dp
!         end do

      endif

      do i = JACT_sparse%IA(J),JACT_sparse%IA(J+1)-1
         PDJ(JACT_sparse%JA(i)) = JACT_sparse%A(i)
      end do

      ! Set error flag to operation completed successfully
      IERR = 0

      end subroutine constV_jac_MEBDF




!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **   Interface for chemistry Jacobian in sparse format         **
!     **              for VODE F90 sparse ODE solver                 **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: tuesday, 22/11/2011                          **
!     **   Usage, copying, distribution of the proprietary routines  **
!     **   or any modified versions of them are not allowed without  **
!     **   approval of the author.                                   **
!     **                                                             **
!     *****************************************************************
      subroutine constV_jac_VODES(neq,time,yin,IA,JA,NZ,PD)

      use sparse_chemistry, only: JAC_sparse, JACT_sparse
      use sparse_algebra

      use ode_solver, only: ncJAC, nIA, nJA, nPD, nVF90JAC, JACT_VF90

      implicit none

!     ** Data input-output ********************************************
      integer,                          intent(in)           :: neq
      double precision,                 intent(in)           :: time
      integer         ,                 intent(inout)        :: NZ
      double precision, dimension(neq), intent(inout)        :: yin
      integer         , dimension(nIA),        intent(inout) :: IA
      integer         , dimension(nJA),        intent(inout) :: JA
      double precision, dimension(nPD),        intent(inout) :: PD

      integer :: i


!        Do assignments
!        NB: Row indexing of sparse(A) is equal to column indexing of
!            sparse(transpose(A))!!
         call constV_jac_sparse(neq,real(time, dp),real(yin, dp))

         if (JACT_sparse%n /= JAC_sparse%n) then
           JACT_sparse = sparse_transpose(JAC_sparse)
         else
           call sparse_transpose_valuesonly(JAC_sparse,JACT_sparse%A)
         endif

         JACT_VF90 = JACT_sparse + identity(neq,tiny(0.e0_dp))
         nz                   = JACT_VF90%n
         IA(1:JACT_VF90%nr+1) = JACT_VF90%IA(1:JACT_VF90%nr+1)
         JA(1:JACT_VF90%n)    = JACT_VF90%JA(1:JACT_VF90%n)
         PD(1:JACT_VF90%n )   = JACT_VF90%A (1:JACT_VF90%n)


         call deallocate(JACT_VF90)


      end subroutine constV_jac_VODES


!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **   Interface for chemistry Jacobian computation by RADAU     **
!     **                    stiff ODE integrator                     **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: sunday, 20/11/2011                           **
!     **   Usage, copying, distribution of the proprietary routines  **
!     **   or any modified versions of them are not allowed without  **
!     **   approval of the author.                                   **
!     **                                                             **
!     *****************************************************************

      subroutine constV_jac_RADAU(neq,time,yin,JAC,nrowJ)
      implicit none

!     ** Data input-output ********************************************
      integer,                          intent(in)         :: neq
      real (dp)       ,                 intent(in)         :: time
      real (dp)       , dimension(neq), intent(inout)      :: yin

      integer,                              intent(inout) :: nrowJ
      real (dp)       , dimension(neq,neq), intent(out)   :: JAC

      integer :: nrowL, ncolU
      data nrowL/0/, ncolU/0/


!     Interface to standard jacobian call
      call constV_jac(neq,time,yin,nrowL,ncolU,JAC,nrowJ)

      end subroutine constV_jac_RADAU


!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **   Interface for chemistry Jacobian computation by DASPK     **
!     **                    stiff ODE integrator                     **
!     **                                                             **
!     **   The output matrix is: A = dG/dY + cj * dG(i)/dYprime(j)   **
!     **                                                             **
!     **   Provided that G(t,y,yprime) = yprime - f(t,y)             **
!     **   we have: A = - JAC(t,y) + cj * identity(neq)              **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: thursday, 24/05/2012                         **
!     **   Usage, copying, distribution of the proprietary routines  **
!     **   or any modified versions of them are not allowed without  **
!     **   approval of the author.                                   **
!     **                                                             **
!     *****************************************************************

      subroutine constV_jac_DASPK(time,yin,yprime,A,cj)

      use speedchem,        only: neq
      use sparse_chemistry, only: JAC_sparse, JACT_sparse
      use sparse_algebra
      use sparse_definitions


      implicit none

!     ** Data input-output ********************************************
      real (dp)       ,                 intent(inout)      :: time
      real (dp)       , dimension(neq), intent(inout)      :: yin, &
                                                              yprime

!     NB: JAC is already initialised to zeroes by DASPK
      real (dp)       , dimension(neq,neq), intent(inout)   :: A

      real (dp)       ,                     intent(inout)   :: cj
      type(sparse) :: ineq, sparseA
      logical      :: tr

!     Interface to standard jacobian call
      call constV_jac_sparse(neq,time,yin)

!     Convert into dense format
      if (JACT_sparse%n /= JAC_sparse%n) then
        JACT_sparse = sparse_transpose(JAC_sparse)
      else
        call sparse_transpose_valuesonly(JAC_sparse,JACT_sparse%A)
      endif

      sparseA = identity(neq,cj) - JACT_sparse
      call sparse_to_dense_columnwise(.true.,sparseA,A)

      end subroutine constV_jac_DASPK

!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **   Interface for chemistry Jacobian computation by DASPK     **
!     **                    stiff ODE integrator                     **
!     **                                                             **
!     **   The output matrix is: A = dG/dY + cj * dG(i)/dYprime(j)   **
!     **                                                             **
!     **   Provided that G(t,y,yprime) = yprime - f(t,y)             **
!     **   we have: A = - JAC(t,y) + cj * identity(neq)              **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: monday,   30/07/2012                         **
!     **   Usage, copying, distribution of the proprietary routines  **
!     **   or any modified versions of them are not allowed without  **
!     **   approval of the author.                                   **
!     **                                                             **
!     *****************************************************************

      subroutine constV_jac_DASPK_sp(time,yin,yprime,A,cj)

      use speedchem,        only: neq
      use sparse_chemistry, only: JAC_sparse
      use sparse_algebra
      use sparse_definitions


      implicit none

!     ** Data input-output ********************************************
      real (dp)       ,                 intent(inout)      :: time
      real (dp)       , dimension(neq), intent(inout)      :: yin, &
                                                              yprime

!     NB: JAC is already initialised to zeroes by DASPK
      type(sparse),                         intent(inout)   :: A

      real (dp)       ,                     intent(inout)   :: cj
      type(sparse) :: ineq
      logical      :: tr

!     Interface to standard jacobian call
      call constV_jac_sparse(neq,time,yin)

!     Output in sparse form
      A = identity(neq,cj) - JAC_sparse

      end subroutine constV_jac_DASPK_sp




!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **   Interface for chemistry Jacobian * vector product         **
!     **     computation by the ROWMAP stiff ODE integrator          **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: tuesday, 22/11/2011                          **
!     **   Usage, copying, distribution of the proprietary routines  **
!     **   or any modified versions of them are not allowed without  **
!     **   approval of the author.                                   **
!     **                                                             **
!     *****************************************************************
      subroutine constV_jacv_ROWMAP(neq,time,yin,vin,jacv)

      use sparse_chemistry, only: JAC_sparse
      use sparse_algebra

      implicit none

!     ** Data input-output ********************************************
      integer,                          intent(in)         :: neq
      real (dp)       ,                 intent(in)         :: time
      real (dp)       , dimension(neq), intent(inout)      :: yin
      real (dp)       , dimension(neq), intent(in   )      :: vin
      real (dp)       , dimension(neq), intent(inout)      :: jacv

!      real (dp)       , dimension(neq,neq), target :: JAC
!      real (dp)       , dimension(neq,neq)   :: TJAC

      integer :: i, j, nrowL = 0, ncolU = 0, nrowJ = 0
      integer :: diagJ =0, moveB = 1

!     Interface to standard jacobian call
      call constV_jac_sparse(neq,time,yin)!,nrowL,ncolU,JAC,nrowJ)
!      TJAC = transpose(JAC)

      JACV = JAC_sparse * vin

!     Compute jacobian - vector product
!      do j = 1, neq
!         jacv(j) = sum(TJAC(:,j) * vin,mask = TJAC(:,j)/=0.e0_dp)
!      end do

      end subroutine constV_jacv_ROWMAP



!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **   Computing chemistry jacobian for constant-volume reactor  **
!     **                                                             **
!     **   constV_jac = [df1/dT  df1/dY1   .... df1/dY_n ]           **
!     **                | ...                            |           **
!     **                [dfm/dT  dfm/dY1   .... dfm/dY_n ]           **
!     **                                                             **
!     **   Dense matrix version interface for traditional solvers    **
!     **                                                             **
!     ** INPUT DATA                                                  **
!     ** yin   = array of the unknowns:                              **
!     **         y(1)           = T [K]                              **
!     **         y(2:nspecie+1) = Y1, Y2, ..., Y_n [-]               **
!     **                                                             **
!     ** OUTPUT DATA                                                 **
!     ** dydt  = unknowns rate of change [K/s],[mol/cm^3/s]          **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: wedesday, 26/05/2010                         **
!     **   Usage, copying, distribution of the proprietary routines  **
!     **   or any modified versions of them are not allowed without  **
!     **   approval of the author.                                   **
!     **                                                             **
!     *****************************************************************


      subroutine constV_jac(neq,time,yin,nrowL,ncolU,JAC,nrowJ)

      use sparse_chemistry, only: JAC_sparse, JACT_sparse
      use sparse_algebra

      implicit none

!     ** Data input-output ********************************************
      integer,                          intent(in)         :: neq
      real (dp)       ,                 intent(in)         :: time
      real (dp)       , dimension(neq), intent(in), target :: yin

      integer,                              intent(inout) :: nrowL, &
                                                             ncolU
      integer,                              intent(inout) :: nrowJ
      real (dp)       , dimension(neq,neq), intent(out), target :: JAC

!     *****************************************************************

      integer :: diagJ = 0, moveB = 1
      logical :: tr

!     Call jacobian matrix in sparse form (updates JAC_sparse)
      call constV_jac_sparse(neq,time,yin)


!     Convert into dense format
      if (.not.(JACT_sparse%n > 0)) then
         tr = .false.
         call sparse_to_dense_columnwise(tr,JAC_sparse,JAC)
      else
         tr = .true.
         call sparse_to_dense_columnwise(tr,JACT_sparse,JAC)
      endif



      end subroutine constV_jac




!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **   Computing chemistry jacobian for constant-volume reactor  **
!     **   Dense matrix version interface for GAM solvers            **
!     **                                                             **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: wedesday, 28/02/2012                         **
!     **   Usage, copying, distribution of the proprietary routines  **
!     **   or any modified versions of them are not allowed without  **
!     **   approval of the author.                                   **
!     **                                                             **
!     *****************************************************************


      subroutine constV_jac_GAM(neq,time,yin,JAC,nrowL,RPAR,IPAR)

      use sparse_chemistry, only: JAC_sparse, JACT_sparse
      use sparse_algebra

      implicit none

!     ** Data input-output ********************************************
      integer,                          intent(in)         :: neq
      real (dp)       ,                 intent(in)         :: time
      real (dp)       , dimension(neq), intent(in), target :: yin

      integer,                              intent(inout) :: nrowL

      real (dp)       , dimension(neq,neq), intent(out), target :: JAC
      real (dp)       , dimension(:), intent(in), optional :: RPAR
      integer,          dimension(:), intent(in), optional :: IPAR

!     *****************************************************************

      integer :: diagJ = 0, moveB = 1
      logical :: tr

!     Call jacobian matrix in sparse form (updates JAC_sparse)
      call constV_jac_sparse(neq,time,yin)

!     Convert into dense format
      if (.not.(JACT_sparse%n > 0)) then
         tr = .false.
         call sparse_to_dense_columnwise(tr,JAC_sparse,JAC)
      else
         tr = .true.
         call sparse_to_dense_columnwise(tr,JACT_sparse,JAC)
      endif



      end subroutine constV_jac_GAM

!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **   Interface for chemistry Jacobian computation by RADAU     **
!     **            stiff ODE integrator in sparse form              **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: friday, 18/05/2012                           **
!     **   Usage, copying, distribution of the proprietary routines  **
!     **   or any modified versions of them are not allowed without  **
!     **   approval of the author.                                   **
!     **                                                             **
!     *****************************************************************

      subroutine constV_jac_RADAUS(neq,time,yin,sparsejac)

      use sparse_chemistry, only: JAC_sparse
      use sparse_algebra

      implicit none

!     ** Data input-output ********************************************
      integer,                          intent(in)         :: neq
      real (dp)       ,                 intent(in)         :: time
      real (dp)       , dimension(neq), intent(inout)      :: yin

      type(sparse),                     intent(out)        :: sparsejac


!     Interface to standard jacobian call
      call constV_jac_sparse(neq,time,yin)

!     Assign output array (DIAGONAL ELEMENTS MUST BE NONZERO!)
      sparsejac = JAC_sparse 

      end subroutine constV_jac_RADAUS



      end module speedchem_conV

!     *****************************************************************

!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **     Interface subroutine for integrating the constant       **
!     **                volume initial value problem                 **
!     **                                                             **
!     **   Author:      Federico Perini                              **
!     **   Last update: tuesday, 01/05/2012                          **
!     **   Usage, copying, distribution of the proprietary routines  **
!     **   or any modified versions of them are not allowed without  **
!     **   approval of the author.                                   **
!     **                                                             **
!     *****************************************************************
      subroutine conV_integrate(unk0,rho0,time0,timef,unkf)

      use working_precision
      use speedchem,       only: ns, neq, species_permutations,       &
                                 species_inverse_permutations
      use chemistry_setup, only: permutate_species
      use chemkin_kiva,    only: CKp, CKrho, intwork, reawork,        &
                                 ckII_conV, ckII_integrate,           &
                                 lenrwk, leniwk
      use chemkinII,       only: ckrhoy, ckpy
      use scmixturethermo, only: SCp, SCrho, rhoY, molar_volumes,     &
                                 pressurerhoT
      use chemistry_setup, only: use_speedchem, analytical_jac
      use ode_solver,      only: rtol, atol
      use omp_lib
      implicit none

      real (dp)       , dimension(neq), intent(in)  :: unk0
      real (dp)       ,                 intent(in)  :: rho0, time0,   &
                                                       timef
      real (dp)       , dimension(neq), intent(out), target :: unkf
      integer                                       :: n

      real (dp)       ,               pointer :: T
      real (dp)       , dimension(:), pointer :: Y

!     *****************************************************************

!     Initialize array
      unkf = unk0

!     Get current thread number
      n = 1!omp_get_thread_num() + 1

!     Pressure value in p0 = [Pa];
!     Initialising average mixture density **************************
        if (use_speedchem) then

!         Density value in SCrho = [kg/m3]; pressure SCP = [Pa];
          SCrho = rho0

!         Update molar volumes matrix for jacobian computations
          if (analytical_jac) call molar_volumes

!         Integrate and retrieve output unkf
          call chemistry_ODE_integrate(neq,rtol,atol,time0,timef,unkf)

        else

!         Insert here the link to your alternative solver    
          stop '[conV_integrate] No alternative solver implemented yet.'

        endif


      end subroutine conV_integrate














