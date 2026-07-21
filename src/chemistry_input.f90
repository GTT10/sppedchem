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
!     **                         SpeedCHEM                           **
!     **                                                             **
!     **   Chemistry ODE integration according to user solver choice **
!     **                                                             **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: tuesday, 22/11/2011                          **
!     **                                                             **
!     *****************************************************************

      subroutine chemistry_ODE_integrate(neq,rtol,atol,t0,tf,yin)

      use working_precision,      only: dp
      use speedchem,              only: nr, njac, ns,                &
                                        species_permutations,        &
                                        species_inverse_permutations
      use ode_solver,             only: rwork, iwork, ncJAC, ncCONV, &
                                        nsteps, nLUDEC, nNewton, hs, &
                                        dummy, ifcn, ijac, imas,     &
                                        iopt, vf90_opts, istate,     &
                                        iout, itask, itol, istats,   &
                                        rstats, maxnsteps, liw, lrw, &
                                        method, mumas, mujac, mljac, &
                                        mlmas, idfx, daspkinfo, maxk
      use chemistry_setup,        only: solver, permutate_species,   &
                                        accurate_scthermo
      use reacpar,                only: n_plog_reactions
      use sparse_chemistry,       only: JAC_sparse, JACT_sparse
      use dvode_f90_m,            only: dvode_f90, vode_opts,    &
                                        USERSETS_IAJA, set_opts, &
                                        get_stats
      use universal_constants,    only: ten
      use sparse_algebra,         only: print_sparse_to_file

      use speedchem_conV

      implicit none


!     ** Problem - related variables
      integer,                   intent(in)    :: neq
      real (dp),                 intent(in)    :: rtol
      real (dp), dimension(neq), intent(in)    :: atol
      real (dp),                 intent(inout) :: t0, tf
      real (dp), dimension(neq), intent(inout) :: yin

      double precision                  :: dt0, dtf
      double precision, dimension(neq)  :: dyin


      real (dp), dimension(neq) :: yin0, yprime
      real (dp)                 :: rpar, rto2, h0, t00
!     Saved initial time so a retry restarts from the original t0. Most
!     solvers take t0 as intent(inout) and advance it, so without this the
!     next attempt would integrate a shorter interval [advanced_t0, tf].
      real (dp)                 :: t0_init
      integer                   :: ipar, it, ierr
	  integer, parameter        :: nit = 5
	  integer, dimension(20)    :: info
	  logical                   :: original_tabulation

!     ** Output formats
      character(len=*), parameter ::                                   &
       fmt_ersol = "(' chemistry solver choice not supported: ',A15)", &
       fmt_erint = "(' error in speedchem chemistry ODE integration')",&
	   fmt_nit   = "(' ODE integration completed in ',I2,' attempts')",&
	   fmt_time  = "(' Warning: integration time is not positive ')"



!     State-vector contract check. The unknown array is
!         yin(1)      = T [K]
!         yin(2:neq)  = species mass fractions Y_1..Y_ns [-]   (NOT molar
!                       concentrations; see SC_conV in SCconV.f90)
!     so neq must be exactly ns+1. A mismatch means the caller sized the
!     state wrong; integrating anyway would read/write out of bounds or
!     silently drop a species, so refuse up front.
      if (neq /= ns + 1) then
         write(*,"(' ERROR chemistry_ODE_integrate: neq (',I0,') /= ns+1"//&
                  " (',I0,'). The state must be [T, Y_1..Y_ns].')") neq, ns+1
         error stop 1
      endif

!     PLOG + analytic Jacobian guard (stage 2 = numeric Jacobian only).
!     PLOG forward rates ARE evaluated in the RHS (mass_action), so
!     numeric-Jacobian solvers integrate PLOG correctly. The analytic
!     Jacobian (constV_jac_sparse) does NOT yet include the PLOG
!     temperature/pressure derivative (stage 3), so a "...JAC" solver
!     would silently use a wrong Jacobian. Refuse rather than integrate
!     with an inconsistent Jacobian. (The check is on the solver name
!     containing "JAC"; the numeric sparse "...S" solvers are fine.)
      if (n_plog_reactions > 0 .and. index(solver, "JAC") > 0) then
         write(*,"(' ERROR chemistry_ODE_integrate: solver ',A,' uses the"//&
                  " analytic Jacobian, but ',I0,' PLOG reaction(s) are"//    &
                  " present and the analytic PLOG Jacobian is not"//         &
                  " implemented yet (stage 3). Use a numeric-Jacobian"//     &
                  " solver (e.g. LSODES/VODE without the JAC suffix).')")     &
                  trim(solver), n_plog_reactions
         error stop 1
      endif

!     Initialization
      istate = 1

!     Preliminary check on integration time
      if ( tf-t0<=tiny(0.d0) ) write(*,fmt_time)

!     Assign system unknowns array: if the species permutation flag
!     is adopted, the order in which species are stored in the reaction
!     mechanism is different than the one it is input by the user.
!     In this case, the species array provided by the user has to be
!     reordered to the new internal species order that reduces fillin
!     during Jacobian matrix decomposition.
      if (permutate_species) then
         yin0 = yin(species_permutations)
      else
         yin0          = yin
      endif

!     Save the initial time so each retry can restart the interval from it
      t0_init = t0

	  it     = 0

!     If the first integration fails, temperature tabulation is suppressed
!     for more accurate solution. Original flag is stored
      original_tabulation = accurate_scthermo


	  integration_attempts: do while (it < nit .and. istate<2)

	     yin     = yin0
		 t0      = t0_init
		 it      = it + 1
		 istate  = 1
		 rto2    = rtol * ten**(it-1)

!        Deactivate tabulation if first attempt has failed
		 if (it > 1) accurate_scthermo = .true.

!     ** MAIN SOLVER CALL **********************************************

!        ** Calling solution for chemistry ODE system within the cell **
         select case (solver)

         case ("VODE")

!           Arbitrary number of steps
            iwork(6) = maxnsteps

            call DVODE (SC_conV,neq,yin,t0,tf,                       &
               ITOL,RTO2,ATOL,ITASK,ISTATE,IOPT,RWORK,LRW,IWORK,      &
               LIW,dummy,method)

            ncJAC   = ncJAC   + IWORK(13)
            ncCONV  = ncCONV  + IWORK(12)
            nsteps  = nsteps  + IWORK(11)
            nLUDEC  = nLUDEC  + IWORK(19)
            nNewton = nNewton + IWORK(20)

         case ("VODEJAC")

!           Arbitrary number of steps
            iwork(6) = maxnsteps

            call DVODE (SC_conV,neq,yin,t0,tf,                        &
               ITOL,RTO2,ATOL,ITASK,ISTATE,IOPT,RWORK,LRW,IWORK,      &
               LIW,constV_jac,method)

            ncJAC   = ncJAC   + IWORK(13)
            ncCONV  = ncCONV  + IWORK(12)
            nsteps  = nsteps  + IWORK(11)
            nLUDEC  = nLUDEC  + IWORK(19)
            nNewton = nNewton + IWORK(20)

         case ("VODESJAC")

           dyin = dble(yin)
           dt0  = dble(t0)
           dtf  = dble(tf)

           call DVODE_F90(SC_conV_VODES,neq,dyin,dt0,dtf,ITASK,ISTATE, &
                          vf90_opts,J_FCN=constV_JAC_VODES)

           t0   = real(dt0, dp)
           tf   = real(dtf, dp)
           yin  = real(dyin, dp)

           call GET_STATS(RSTATS, ISTATS)
           ncJAC   = ncJAC   + ISTATS(13)
           ncCONV  = ncCONV  + ISTATS(12)
           nsteps  = nsteps  + ISTATS(11)
           nLUDEC  = nLUDEC  + ISTATS(19)
           nNewton = nNewton + ISTATS(20)


         case ("VODES")

           dyin = dble(yin)
           dt0  = dble(t0)
           dtf  = dble(tf)
           call DVODE_F90(SC_conV_VODES,neq,dyin,dt0,dtf,ITASK,ISTATE,      &
                          vf90_opts)

           t0   = real(dt0, dp)
           tf   = real(dtf, dp)
           yin  = real(dyin, dp)

           call GET_STATS(RSTATS, ISTATS)
           ncJAC   = ncJAC   + ISTATS(13)
           ncCONV  = ncCONV  + ISTATS(12)
           nsteps  = nsteps  + ISTATS(11)
           nLUDEC  = nLUDEC  + ISTATS(19)
           nNewton = nNewton + ISTATS(20)

         case ("LSODES")

!           Arbitrary number of steps
            iwork(6) = maxnsteps


            call DLSODES (SC_conV,neq,yin,t0,tf,                       &
                  ITOL,RTO2,ATOL,ITASK,                                &
                  ISTATE, IOPT, RWORK, LRW, IWORK, LIW, dummy, method)

            ncJAC   = ncJAC  + IWORK(13)
            ncCONV  = ncCONV + IWORK(12)
            nsteps  = nsteps + IWORK(11)
            nLUDEC  = nLUDEC + IWORK(21)


         case ("LSODESJAC")

!           Arbitrary number of steps
            iwork(6) = maxnsteps

!           Input jacobian matrix sparsity structure
            iwork(30+1:30+neq+1) = JAC_sparse%IA(1:neq+1)
            iwork(31+neq+1:31+neq+njac) = JAC_sparse%JA(1:JAC_sparse%n)

            call DLSODES (SC_conV,neq,yin,t0,tf,                       &
                  ITOL,RTO2,ATOL,ITASK, ISTATE, IOPT,                  &
                  RWORK, LRW, IWORK, LIW, constV_JAC_LSODES, method)

            ncJAC   = ncJAC  + IWORK(13)
            ncCONV  = ncCONV + IWORK(12)
            nsteps  = nsteps + IWORK(11)
            nLUDEC  = nLUDEC + IWORK(21)

         case ("LSODE")

            iwork(6) = maxnsteps
            call DLSODE  (SC_conV,neq,yin,t0,tf,                       &
                  ITOL,RTO2,ATOL,ITASK,                                &
                  ISTATE, IOPT, RWORK, LRW, IWORK, LIW, dummy, method)

            ncJAC   = ncJAC  + IWORK(13)
            ncCONV  = ncCONV + IWORK(12)
            nsteps  = nsteps + IWORK(11)
            nLUDEC  = nLUDEC + IWORK(13)


         case ("LSODEJAC")

            iwork(6) = maxnsteps
            call DLSODE  (SC_conV,neq,yin,t0,tf,                       &
                  ITOL,RTO2,ATOL,ITASK, ISTATE,                        &
                  IOPT, RWORK, LRW, IWORK, LIW, constV_jac, method)

            ncJAC   = ncJAC + IWORK(13)
            ncCONV  = ncCONV + IWORK(12)
            nsteps  = nsteps + IWORK(11)


         case ("LSODA")

            call DLSODA  (SC_conV,neq,yin,t0,tf,                       &
                  ITOL,RTO2,ATOL,ITASK, istate,                        &
                  IOPT, RWORK, LRW, IWORK, LIW, constV_jac, method)

         case ("LSODAJAC")

            call DLSODA  (SC_conV,neq,yin,t0,tf,                       &
                  ITOL,RTO2,ATOL,ITASK, ISTATE,                        &
                  IOPT, RWORK, LRW, IWORK, LIW, constV_jac, method)

         case ("RADAU5")


            call RADAU5 (neq,SC_conV,t0,yin,tf,                        &
                        hs,spread(RTO2,1,neq),ATOL,ITOL,               &
                        dummy ,IJAC, neq, MUJAC,                       &
                        dummy ,IMAS, neq, MUMAS,                       &
                        dummy ,IOUT,                                   &
                        RWORK,LRW,IWORK,LIW,RPAR,IPAR,istate)

            if (istate==1) istate = 100

            ncJAC   = ncJAC + IWORK(15)
!           Here function evaluations for constructing jacobian aren't counted
            ncCONV  = ncCONV + IWORK(14)
            nsteps  = nsteps + IWORK(16)
            nLUDEC  = nLUDEC + IWORK(19)

         case ("RADAU5JAC")


            hs = tf - t0
            call RADAU5 (neq,SC_conV,t0,yin,tf,                        &
                        hs,spread(RTO2,1,neq),ATOL,ITOL,               &
                        constV_jac_RADAUS,IJAC, mljac, MUJAC,          &
                        dummy ,IMAS, neq, MUMAS,                       &
                        dummy ,IOUT,                                   &
                        RWORK,LRW,IWORK,LIW,RPAR,IPAR,istate)

            if (istate==1) istate = 100

            ncJAC   = ncJAC + IWORK(15)
!           Here function evaluations for constructing jacobian aren't counted
            ncCONV  = ncCONV + IWORK(14)
            nsteps  = nsteps + IWORK(16)
            nLUDEC  = nLUDEC + IWORK(19)

         case ("RODAS")

            iwork(1) = maxnsteps

            call RODAS (neq,SC_conV,itask,t0,yin,tf,                   &
                        hs,spread(RTO2,1,neq),ATOL,ITOL,               &
                        dummy ,IJAC, mljac, MUJAC,                     &
                        dummy ,IDFX,                                   &
                        dummy ,IMAS, mlmas, MUMAS,                     &
                        dummy ,IOUT,                                   &
                        RWORK,LRW,IWORK,LIW,RPAR,IPAR,istate)

            if (istate==1) istate = 100

            ncJAC   = ncJAC + IWORK(15)
!           Here function evaluations for constructing jacobian aren't counted
            ncCONV  = ncCONV + IWORK(14)
            nsteps  = nsteps + IWORK(16)
            nLUDEC  = nLUDEC + IWORK(19)

         case ("RODASJAC")

            iwork(1) = maxnsteps

            call RODAS (neq,SC_conV,itask,t0,yin,tf,                   &
                        hs,spread(RTO2,1,neq),ATOL,ITOL,               &
                        constV_jac_RADAU ,IJAC, mljac, MUJAC,          &
                        dummy ,IDFX,                                   &
                        dummy ,IMAS, mlmas, MUMAS,                     &
                        dummy ,IOUT,                                   &
                        RWORK,LRW,IWORK,LIW,RPAR,IPAR,istate)

            if (istate==1) istate = 100

            ncJAC   = ncJAC + IWORK(15)
!           Here function evaluations for constructing jacobian aren't counted
            ncCONV  = ncCONV + IWORK(14)
            nsteps  = nsteps + IWORK(16)
            nLUDEC  = nLUDEC + IWORK(19)

         case ("ROWMAP")

            call ROWMAP(neq,SC_conV,ifcn,t0,yin,                       &
                        tf,hs,rto2,atol(1),itol,                       &
                        dummy,ijac,dummy,ifcn,dummy,iout,rwork,        &
                        lrw,iwork,liw,rpar,ipar,istate)

            if (istate==1) istate = 100

            ncJAC   = ncJAC + IWORK(8)
!           Here function evaluations for constructing jacobian aren't counted
            ncCONV  = ncCONV + IWORK(7)
            nsteps  = nsteps + IWORK(5)

         case ("ROWMAPJAC")

            call ROWMAP(neq,SC_conV,ifcn,t0,yin,                       &
                        tf,hs,rto2,atol(1),itol,                       &
                        constV_jacv_ROWMAP,ijac,dummy,ifcn,dummy,iout, &
                        rwork,lrw,iwork,liw,rpar,ipar,istate)

            if (istate==1) istate = 100

            ncJAC   = ncJAC + IWORK(8)
!           Here function evaluations for constructing jacobian aren't counted
            ncCONV  = ncCONV + IWORK(7)
            nsteps  = nsteps + IWORK(5)


         case ("GAM")

            iwork(2) = maxnsteps

            call GAM(neq,SC_conV,t0,yin,tf,hs,rto2,atol(1),   &
                     itol, dummy, ijac, mljac, mujac, dummy,  &
                     iout, rwork, lrw, iwork, liw, rpar, ipar,&
                     istate)

            if (istate==1)istate = 2

            ncJAC   = ncJAC + IWORK(11)
            ncCONV  = ncCONV + IWORK(10)
            nsteps  = nsteps + sum(IWORK(12:15))
            nLUDEC  = nLUDEC + IWORK(24)


         case ("GAMJAC")

            iwork(2) = maxnsteps

            call GAM(neq,SC_conV,t0,yin,tf,hs,rto2,atol(1),   &
                     itol, constV_jac_GAM, ijac, mljac, mujac,&
                     dummy, iout, rwork, lrw, iwork, liw,     &
                     rpar, ipar, istate)

            if (istate==1)istate = 2

            ncJAC   = ncJAC + IWORK(11)
            ncCONV  = ncCONV + IWORK(10)
            nsteps  = nsteps + sum(IWORK(12:15))
            nLUDEC  = nLUDEC + IWORK(24)

         case ("DASPK")

            iwork(1:liw) = 0
            info = daspkinfo

            call SC_conV(neq,t0,yin,yprime)
            call DDASPK(SC_conV_DAE,neq,t0,yin,yprime,tf,     &
                        info,spread(rto2,1,neq),atol,         &
                        istate, rwork, lrw, iwork, liw,       &
                        rpar, ipar, dummy, dummy)

            ncJAC   = ncJAC  + IWORK(13)
            ncCONV  = ncCONV + IWORK(12) + 1
            nsteps  = nsteps + IWORK(11)
            nLUDEC  = nLUDEC + IWORK(21)

         case ("DASPKJAC")

            iwork(1:liw) = 0
            info = daspkinfo

            call SC_conV(neq,t0,yin,yprime)
            call DDASPK(SC_conV_DAE,neq,t0,yin,yprime,tf,     &
                        info,spread(rto2,1,neq),atol,         &
                        istate, rwork, lrw, iwork, liw,       &
                        rpar, ipar, constV_jac_DASPK, dummy)

            ncJAC   = ncJAC  + IWORK(13)
            ncCONV  = ncCONV + IWORK(12) + 1
            nsteps  = nsteps + IWORK(11)
            nLUDEC  = nLUDEC + IWORK(21)


         case ("MEBDF")

            write(*,*)' MEBDF solver only supported w/ '//    &
                      ' analytical Jacobian choice '
            stop

         case ("MEBDFJAC")

            ierr = 0

            iwork(11) = maxnsteps
            iwork(1:10) = 0

            t00 = t0
            h0  = min(hs, tf-t00)

            call MEBDFSO(neq,t00,h0,yin,tf,tf,                &
                         method,istate,6,lrw,rwork,liw,       &
                         iwork,maxk,itol,rto2,atol(1),        &
                         SC_conV_MEBDF,constV_jac_MEBDF,      &
                         ijac,njac+neq,ipar,rpar,ierr)


            istate = 0

            call MEBDFSO(neq,t00,h0,yin,tf,tf,                &
                         method,istate,6,lrw,rwork,liw,       &
                         iwork,maxk,itol,rto2,atol(1),        &
                         SC_conV_MEBDF,constV_jac_MEBDF,      &
                         ijac,njac+neq,ipar,rpar,ierr)

            ncJAC   = ncJAC  + IWORK( 5)
            ncCONV  = ncCONV + IWORK( 4)
            nsteps  = nsteps + IWORK( 2)
            nLUDEC  = nLUDEC + IWORK( 6)

!            write(*,*)IWORK(1:12)

            ! IF Istate=0, normal output
            if (istate>=0) istate = 2
            if (ierr/=0  ) istate = -1

         case default

            write(*,fmt_ersol)solver
            stop

         end select

      end do integration_attempts


!     Restore the original tabulation flag
      if (it > 1) accurate_scthermo = original_tabulation

!     Restore the user's species order if permutation during integration
!     had been adopted
      if (permutate_species) then
         yin  = yin(species_inverse_permutations)
      endif


!     ** END OF MAIN SOLVER CALL ***************************************
      if (it > 1) write(*, fmt_nit)it

      if (istate < 2) then
         write(*, fmt_erint)
         write(*,*)'yin0',yin0(1:neq)
         stop
      endif

      end subroutine chemistry_ODE_integrate

!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **     Initialization of ODE solver runtime parameters         **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: 24/04/2012                                   **
!     **   Usage, copying, distribution of the proprietary routines  **
!     **   or any modified versions of them are not allowed without  **
!     **   approval of the author.                                   **
!     **                                                             **
!     *****************************************************************
      subroutine ODE_solver_speedchem_init

      use working_precision,only: dp
      use speedchem,        only: neq, nr, njac
      use chemistry_setup,  only: solver, TOLR, YTOLA, TTOLA,         &
                                  separate_tols
      use sparse_chemistry, only: JAC_sparse, JACT_sparse
      use ode_solver,       only: lrw, liw, itol, iopt, itask, istate,&
                                  maxnsteps, method, vf90_opts, ijac, &
                                  rtol, atol, mujac, mumas, mljac,    &
                                  iout, imas, hs, ifcn, maxk,         &
                                  rwork, iwork, idfx, mlmas,          &
                                  daspkinfo,nIA, nJA, nPD, nVF90JAC
      use universal_constants, only: milli
      use sparse_definitions
      use sparse_algebra,   only: sparse_transpose, identity
      use dvode_f90_m,      only: usersets_iaja, set_opts
      use utilities,        only: force_allocate

      implicit none

      integer      :: i, ipar, ierr = 0, j = 0
      external     :: SC_conV_MEBDF, constV_JAC_MEBDF
      real (dp)    :: rpar
      type(sparse) :: vf90_jac

      character(len=*), parameter ::                                  &
        fmt_ersol = "(' ODE solver choice not supported: ',A15)",     &
        fmt_erint = "(' chemistry still not initialised. ')",         &
        fmt_erjac = "(' chemistry jacobian still not initialised. ')",&
        fmt_eline = "(' solver unsuitable for Jacobian with empty')", &
        fmt_elin2 = "(' lines, switching to VODE ')"


!       Problem initialisation check
        if (.not.neq>0) then
           write(*,fmt_erint)
           stop
        endif

!       ** Solver-independent runtime parameters and tolerances
        maxnsteps = 50000

        call force_allocate(atol,neq)
                           atol    = YTOLA
                           rtol    =  TOLR
        if (separate_tols) atol(1) = TTOLA



!       ** Definition of integrator-specific parameters and work array**
!       ** dimensions                                                 **
        integrator: select case (solver)

!         ****   VODE   ************************************************
          case ("VODE")

            lrw    = 22 + 9*NEQ + 2*NEQ**2
            liw    = 30 + neq
            method = 22
            itol   = 2
            iopt   = 1
            itask  = 1

!         ****   VODE + JACOBIAN   *************************************
          case ("VODEJAC")

            lrw    = 22 + 9*NEQ + 2*neq**2
            liw    = 30 + neq
            method = 21
            itol   = 2
            iopt   = 1
            itask  = 1

!         ****   VODE FORTRAN 90 (SPARSE) + JACOBIAN   *****************
          case ("VODESJAC")

            itask    = 1
            method   = 26

            vf90_opts = SET_OPTS( SPARSE_J               = .true.,     &
                                  METHOD_FLAG            = method,     &
                                  MXSTEP                 = maxnsteps,  &
                                  NZSWAG                 = njac+neq,   &
                                  ABSERR_VECTOR          = dble(ATOL), &
                                  RELERR                 = dble(RTOL), &
                                  MA28_ELBOW_ROOM        = 2,          &
                                  MC19_SCALING           = .false.,    &
                                  MA28_MESSAGES          = .false.,    &
                                  MA28_EPS               = 1.d-16,     &
                                  MA28_RPS               = .false.,    &
                                  USER_SUPPLIED_SPARSITY = .true.,     &
                                  USER_SUPPLIED_JACOBIAN = .true.,     &
                                  SAVE_JACOBIAN          = .true.      )

!           Prepare sparse jacobian features
            if (.not.JAC_sparse%n > 0) then
               write(*,fmt_erjac)
               stop
            else

!              Compute transposed Jacobian matrix
               JACT_sparse = sparse_transpose(JAC_sparse)

!              Add diagonal elements (required by VODE_F90)
               vf90_jac = JACT_sparse + identity(neq)

               call USERSETS_IAJA   (  IAUSER  = vf90_jac%IA,   &
                                       NIAUSER = vf90_jac%nr+1, &
                                       JAUSER  = vf90_jac%JA,   &
                                       NJAUSER = vf90_jac%n     )

!              Store array dimensions for constV_jac_VODES
               nIA = vf90_jac%nr+1
               nJA = vf90_jac%n
               nPD = vf90_jac%n
               nVF90JAC = vf90_jac%n


               call deallocate(vf90_jac)
            endif

!         ****   VODE FORTRAN 90 (SPARSE)   ****************************
          case ("VODES")

            itask  = 1
            method = 22

            vf90_opts = set_opts( SPARSE_J   = .true.,              &
                                  USER_SUPPLIED_SPARSITY = .false., &
                                  USER_SUPPLIED_JACOBIAN = .false., &
                                  METHOD_FLAG          = method,    &
                                  MXSTEP               = maxnsteps, &
                                  ABSERR_VECTOR        = dble(ATOL),&
                                  RELERR               = dble(RTOL),&
                                  MA28_ELBOW_ROOM      = 2,         &
                                  MC19_SCALING         = .false.,   &
                                  MA28_MESSAGES        = .false.,   &
                                  MA28_EPS             = 1.d-04,    &
                                  MA28_RPS             = .true.,    &
                                  SAVE_JACOBIAN        = .true.     )

!         ****   LSODES   **********************************************
          case ("LSODES")

            lrw    = 20 + 21 * neq + 4*njac
            liw    = 30
            itol   = 2
            method = 222
            iopt   = 1
            itask  = 1

!         ****   LSODES + JACOBIAN   ***********************************
          case ("LSODESJAC")

            lrw    = 20 + 21 * neq + 4*njac
            liw    = 31 + neq + njac
            itol   = 2
            method = 21
            iopt   = 1
            itask  = 1

!         ****   LSODE   ***********************************************
          case ("LSODE")

            lrw    = 22 +  9*NEQ + NEQ**2
            liw    = 20 + NEQ
            method = 22
            itol   = 2
            iopt   = 1
            itask  = 1

!         ****   LSODE + JACOBIAN   ************************************
          case ("LSODEJAC")

            lrw    = 22 +  9*NEQ + NEQ**2
            liw    = 20 + NEQ
            method = 21
            itol   = 2
            iopt   = 1
            itask  = 1

!         ****   LSODA   ***********************************************
          case ("LSODA")

            lrw    = 22 +  9*NEQ + NEQ**2
            liw    = 20 + NEQ
            method = 2
            itol   = 2
            iopt   = 1
            itask  = 1

!         ****   LSODA + JACOBIAN   ************************************
          case ("LSODAJAC")

            lrw    = 22 +  9*NEQ + NEQ**2
            liw    = 20 + NEQ
            method = 2
            itol   = 2
            iopt   = 1
            itask  = 1

!         ****   RADAU5   **********************************************
          case ("RADAU5")

            lrw   = 4*neq**2 + 12*neq + 20
            liw   = 3*neq + 20
            itol  = 1

            hs    = 1.d-3!min(1.d-3, tf-t0)
            ijac  = 0
            imas  = 0
            iout  = 0
            mljac = neq
            mujac = neq
            mumas = neq

!         ****   RODAS    **********************************************
          case ("RODAS")

            itask = 0 !ODE system is autonomous
            hs    = 1.d-3!min(1.d-3, tf-t0)
            itol  = 1

            mljac = neq
            mujac = neq
            mlmas = neq
            mumas = neq

            ijac  = 0
            idfx  = 0
            iout  = 0
            imas  = 0

            lrw   = 4*neq**2 + 12*neq + 20
            liw   = neq + 20

!         ****   RADAU5 + JACOBIAN   ***********************************
          case ("RADAU5JAC")

            lrw   = 4*neq**2 + 12*neq + 20
            liw   = 3*neq + 20
            itol  = 1

            hs    = 1.d-3
            ijac  = 1
            imas  = 0
            iout  = 0
            mljac = neq
            mujac = neq
            mumas = neq


!           Prepare allocation for sparse jacobian exploitation
!            write(*,*)'Allocate sparse algebra solver...'
!            call sparse_radau5_algebra_allocation(neq,njac)
!            write(*,*)'done.'

!         ****   RODAS + JACOBIAN   ************************************
          case ("RODASJAC")

            itask = 0 !ODE system is autonomous
            hs    = milli!min(1.d-3, tf-t0)
            itol  = 1

            mljac = neq
            mujac = neq
            mlmas = neq
            mumas = neq

            ijac  = 1
            idfx  = 0
            iout  = 0
            imas  = 0

            lrw   = 4*neq**2 + 12*neq + 20
            liw   = neq + 20

!         ****   ROWMAP   **********************************************
          case ("ROWMAP")

           hs    = milli
           ifcn  = 0  ! system is autonomous
           itol  = 0  ! tolerance constraints are scalars
           ijac  = 0  ! jacobian computed internally
           iout  = 0  ! no intermediate output
           maxk  = 70 ! maximum krylov dimension
           lrw   = 10+neq*(maxk+11)+maxk*(maxk+4)
           liw   = maxk + 20

!         ****   MEBDF    **********************************************
          case ("MEBDFJAC")

           hs     = milli
           method = 25
           maxk   = 7 ! Maximum integration order
           ijac   = 2*(8*neq+2+2*(njac+neq)) ! Space for sparse Yale package
           lrw    = 33*neq + 2*(njac+neq) + ijac + 3
           liw    = 6*neq + 2*(njac+neq) + 15
           itol   = 2 ! RTOL=scalar, ATOL=scalar


!         ****   ROWMAP + JACOBIAN   ***********************************
          case ("ROWMAPJAC")

           hs    = milli
           ifcn  = 0  ! system is autonomous
           itol  = 0  ! tolerance constraints are scalars
           ijac  = 1  ! jacobian provided analytically
           iout  = 0  ! no intermediate output
           maxk  = 70 ! maximum krylov dimension
           lrw   = 10+neq*(maxk+11)+maxk*(maxk+4)
           liw   = maxk + 20

!         ****   GAM   *************************************************
          case ("GAM")

            hs   = milli
            itol = 0
            ijac = 0
            mljac = neq
            mujac = neq
            iout = 0
            lrw  = 2*neq**2 + 42*neq + 18
            liw  = 24 + neq

!         ****   GAM + JACOBIAN   **************************************
          case ("GAMJAC")

            hs   = milli
            itol = 0
            ijac = 1
            mljac = neq
            mujac = neq
            iout = 0
            lrw  = 2*neq**2 + 42*neq + 18
            liw  = 24 + neq

!         ****   DASPK   ***********************************************
          case ("DASPK")


            liw = 40 + neq
            lrw = 50 + 9*neq + neq**2

            daspkinfo(1)  = 0 ! It is initialisation
            daspkinfo(2)  = 1 ! RTOL, ATOL are arrays
            daspkinfo(3)  = 0 ! no intermediate output
            daspkinfo(4)  = 0
            daspkinfo(5)  = 0 ! Finite-differences jacobian
            daspkinfo(6)  = 0 ! Jacobian is not banded
            daspkinfo(7)  = 0 ! No maximum stepsize given
            daspkinfo(8)  = 0 ! Initial step computed automatically
            daspkinfo(9)  = 0 ! Use default max order (=5th)
            daspkinfo(10) = 0 ! Do not check for negative items
            daspkinfo(11) = 0 ! Initial Y and Y' are consistent
            daspkinfo(12) = 0 ! Do not use Krylov techniques
            daspkinfo(13:20) = 0 ! Used only on Krylov cases

!         ****   DASPK + JACOBIAN (no KRYLOV)   ************************
          case ("DASPKJAC")

            liw = 40 + neq
            lrw = 50 + 9*neq + neq**2

            daspkinfo(1)  = 0 ! It is initialisation
            daspkinfo(2)  = 1 ! RTOL, ATOL are arrays
            daspkinfo(3)  = 0 ! no intermediate output
            daspkinfo(4)  = 0
            daspkinfo(5)  = 1 ! Finite-differences jacobian
            daspkinfo(6)  = 0 ! Jacobian is not banded
            daspkinfo(7)  = 0 ! No maximum stepsize given
            daspkinfo(8)  = 0 ! Initial step computed automatically
            daspkinfo(9)  = 0 ! Use default max order (=5th)
            daspkinfo(10) = 0 ! Do not check for negative items
            daspkinfo(11) = 0 ! Initial Y and Y' are consistent
            daspkinfo(12) = 2 ! Use SPARSE jacobian
            daspkinfo(13:20) = 0 ! Used only on Krylov cases

          case default

             write(*,fmt_ersol)solver
             stop

        end select integrator

!       Allocate solver working arrays
        call force_allocate(rwork, lrw)
        call force_allocate(iwork, liw)
        rwork = 0.d0
        iwork = 0

        end subroutine ODE_solver_speedchem_init




!     *****************************************************************
!     **                                                             **
!     **                         SpeedCHEM                           **
!     **                                                             **
!     **                Combustion mechanism input                   **
!     **                                                             **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: tuesday, 15/11/2011                          **
!     **   Usage, copying, distribution of the proprietary routines  **
!     **   or any modified versions of them are not allowed without  **
!     **   approval of the author.                                   **
!     **                                                             **
!     *****************************************************************

      subroutine chemistry_input

      use working_precision,only: dp
      use chemistry_setup,  only: solver_setup, mechanism, use_speedchem, &
                                  accurate_scthermo,                      &
                                  check_reaction_mechanism,               &
!ck2015                                  permutate_species
                                  permutate_species, mechdir
      use speedchem,        only: ns, nr, neq, scmechanism_to_file,       &
                                  init_stoich_indices,                    &
                                  scmechanism_to_chemkin,                 &
                                  species_permutations,                   &
                                  permutate_mechanism_species,            &
                                  species_inverse_permutations,           &
                                  nnotrev, inotrev, Ainf,                 &
                                  check_duplicate_reactions,              &
                                  compute_optimal_ordering
      use SCthermodata,     only: tabulate_SCthermo, store_thermo_coeffs, &
                                  element_matrix, element_mass_fracs,     &
                                  check_atom_conservation,                &
                                  SCthermo_to_janaf, permutate_thermo,    &
                                  check_janaf_polynomials
      use sparse_chemistry, only: permutate_sparse_matrices
      use SCmixturethermo,  only: SCP, SCrho
      use reacpar,          only: init_reaction_indices, &
                                  tabulate_equilibrium
      use troepar,          only: init_thirdbody_indices, &
                                  tabulate_troepars
      use kinetics_mod,     only: tabulate_kinetics, permutate_kinetics
      use sparse_chemistry, only: sparse_chemistry_setup, JAC_sparse
      use sparse_algebra,   only: print_sparse_to_file, dense_to_sparse
      use ode_solver,       only: set_parallel
      use speedchem_conV,   only: constV_jac_sparse
      use universal_constants, only: one, two, ten, atm_to_Pa
!ck2014 add 3 lines
      use chemkinii
      use chemkinii_interpreter
      use chemkin_kiva

      implicit none

!     Variables related to chemistry jacobian sparsity evaluation
      real (dp)        :: dummyt = 0.e0_dp
      real (kind = 4)  :: t0, t1
      integer          :: i, j, dummynJ=0, dummyn1=0, dummyn2=0
      logical          :: present, present_CK1, present_CK2, present_CKascii,    &
                          present_CKtherm, ckinp_run = .false.
      real (dp), dimension(:)  , allocatable :: dummyy, dummydyindt
      real (dp), dimension(:,:), allocatable :: dummyjac


      character(len=*), parameter :: &
        fmt_bar = "(' ------------------------------------------------------')", &
        fmt_lnk = "(' chemistry link completed: ',I4,' reactions, ',I4, "        &
                    //"' species')", &
        fmt_ttl = "('       SpeedCHEM - Constant Volume Reactor             ')", &
        fmt_lab = "(' Mechanism: ',A80)",                                        &
        fmt_ckintp = "(' calling CHEMKIN-II interpreter - therm.dat chem.inp')", &
        fmt_chkck  = "(' checking CHEMKIN-II linking file...')",                 &
        fmt_chkok  = "(' no errors found in CHEMKIN-II linking file.')",         &
        fmt_chkno  = "(' errors found in chem.out, exiting. ')",                 &
        fmt_tabt   = "(' tabulation time: ',f7.3,' s.')",                        &
        fmt_parst  = "(' parsing    time: ',f7.3,' s.')",                        &
        fmt_parse  = "(' parsing reaction mechanism... ')"


      write(*,*)
      write(*,fmt_bar)
      write(*,fmt_ttl)
      write(*,fmt_bar)
      write(*,*)

!  allocate(species_permutations(161),species_inverse_permutations(161))
!      open(unit=123,file='ordine2.dat')
!      read(123,*)(species_permutations(j),j=1,160)
!      close(123)
!
!      species_permutations(1:160) = species_permutations(1:160) + 1
!      species_permutations(161)   = 1

!

!      species_inverse_permutations(species_permutations) = [(j,j=1,161)]


!     Load chemistry solution setup
      call solver_setup

!     Choose number of threads to be used
      call set_parallel(0)

!     Load mechanism and thermodynamic data.
!     Use CHEMKIN linking file if possible, in binary or
!     in ASCII forms; otherwise, look for data in
!     "SpeedCHEM" ascii format
!ck2015s      inquire(file='cklink'  ,  exist = present_CK1)
!      inquire(file='chem.bin',  exist = present_CK2)
!      inquire(file='chem.inp',  exist = present_CKascii)
!      inquire(file='therm.dat', exist = present_CKtherm)
      inquire(file=trim(mechdir)//"cklink"  ,  exist = present_CK1)
      inquire(file=trim(mechdir)//"chem.bin",  exist = present_CK2)
      inquire(file=trim(mechdir)//"chem.inp",  exist = present_CKascii)
      inquire(file=trim(mechdir)//"therm.dat", exist = present_CKtherm)
!ck2015e
!     Parse reaction mechanism
      call cpu_time(t0)
      write(*,fmt_parse)
      if (present_CK1 .or. present_CK2) then
         call SCcklink
      elseif (present_CKascii .and. present_CKtherm) then
         call ckintp
         call SCcklink
         ckinp_run = .true.
      else
         call SCsetup
      endif

!     Initialize sparse chemistry algebra
!      call init_stoich_indices
      call init_reaction_indices
      call init_thirdbody_indices
!      if (accurate_scthermo) call store_thermo_coeffs
!     Compute mechanism's sparse matrices
      call sparse_chemistry_setup(inotrev,nnotrev,Ainf)

!     Compute stoichiometry: species and element invariants
      call element_matrix
!      call stoich_afr
      call element_mass_fracs

      if (check_reaction_mechanism) call check_atom_conservation
      call check_duplicate_reactions
      call check_janaf_polynomials
      close(820)

      call cpu_time(t1)
      write(*,fmt_parst)t1-t0

!     Finally, print mechanism information to file
      call SCmechanism_to_file
      call SCmechanism_to_chemkin
      call SCthermo_to_janaf

!     Deallocate unuseful storage
      call SCdeallocate

!     Tabulate thermodynamic, equilibrium and kinetics functions for
!     faster computation

      if (.not.accurate_scthermo) then
         call cpu_time(t0)
         call tabulate_SCthermo
         call tabulate_equilibrium
         call tabulate_troepars
         call tabulate_kinetics
         call cpu_time(t1)
         write(*,fmt_tabt)t1-t0
      endif

!     Evaluate chemistry jacobian sparsity pattern: in case OpenMP execution
!     is chosen, initialisation has to be done per each thread, in order to
!     allocate the sparse arrays for jacobian computation
      write(*,*)
      write(*,"(1X,A80)")mechanism
      write(*,*)
      allocate( dummyy(neq) )

      if (.not.permutate_species) then
         dummyy = [1500.0_dp,(one/real(ns, dp),i=1,ns)]
      else
         dummyy = [(one/real(ns, dp),i=1,ns), 1500.0_dp]
      endif

!$OMP PARALLEL
!$OMP CRITICAL
      SCP   = ten * atm_to_pa ! [Pa]
      SCrho = two             ! [kg/m3]

      call constV_jac_sparse(neq,dummyt,dummyy)
      call jac_sparsity(neq)
!     Initialize ODE solver parameters
!     (needs to be done after call to the jacobian routine!)
      call ODE_solver_speedchem_init
!$OMP END CRITICAL
!$OMP END PARALLEL

!     Introduce optimal permutation
      call compute_optimal_ordering

!     Deallocate unuseful dense matrix storage
      deallocate(dummyy)


!     ** CHEMKIN-II related initialization
      use_chemkin: if (.not.use_speedchem) then

         write(*,*)
         write(*,fmt_ckintp)
         if (.not.ckinp_run.and.(.not.(present_CK1 .or. present_CK2))) then

            call ckintp
            write(*,fmt_chkck)

!           Check output file
!ck2014            inquire (file = chemdat, exist = present)
            inquire (file = trim(mechdir)//"chem.out", exist = present)

            if (present) then
              write(*,fmt_chkok)

            else
              write(*,fmt_chkno)
              stop
            endif

         endif

!        Initialize chemkin calculations
         call chemkin_initialize

      endif use_chemkin


      write(*,*)
      write(*,fmt_lnk) nr,ns
      write(*,fmt_bar)
      write(*,*)



      end subroutine chemistry_input