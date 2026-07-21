!     ******************************************************************
!     **                                                              **
!     **                           ICSPLODE                           **
!     **        Implicit Chemistry SParse soLver for ODE systems      **
!     **                                                              **
!     **   A Fortran 2008 implementation of a sparse, implicit Runge  **
!     **   Kutta of order 5 solver for stiff systems of Ordinary      **
!     **   Differential Equations                                     **
!     **                                                              **
!     **   Authors:     Emanuele Galligani, Federico Perini           **
!     **   Last update: wednesday, 29/02/2012                         **
!     **                                                              **
!     ******************************************************************

      module radau_sparse

         use sparse_algebra, only: sparse
         implicit none
         private

!        ***************************************************************
!        **  SOLVER WORKING PRECISION (default: double precision)     **
!        ***************************************************************
         integer, parameter :: dp = KIND(0.d0)

!        ***************************************************************
!        **  SOLVER'S DERIVED DATA TYPES                              **
!        ***************************************************************

!        1) Derived type for integration timestep
         type step
            sequence
!           step count
            integer   :: n
!           step width [s]
            real (dp) :: h
!           logical flags
            logical   :: isfirst
            logical   :: islast
            logical   :: accepted
         end type step

!        ***************************************************************
!        **  COMMON SIMULATION PARAMETERS AND OPTIONS                 **
!        ***************************************************************

!        Number of equations and its powers
         integer :: neq, neq2, neq3, neq4

!        Jacobian matrix evaluation: numerical/analytical
         logical :: analytical_jac = .false.
         logical :: sparse_jac     = .false.

!        ***************************************************************
!        **  WORKING ARRAYS AND COMMON VARIABLES                      **
!        ***************************************************************

!        Jacobian matrix storage
         real (dp), dimension(:,:), allocatable :: JAC
         type (sparse)                          :: JACS

!        ** Time advancement and control *******************************
         type (step) :: dt, dtold, dtnew
         real (dp)   :: t, tstart, tend, ttoend
         real (dp)   :: hmin, hmax


!        ***************************************************************
!        **  INTEGRATION STATISTICS                                   **
!        ***************************************************************

!        Number of jacobian evaluations
         integer :: njac = 0

!        Accepted and rejected integration steps
         integer :: naccst = 0, nrejst = 0


!        ***************************************************************
!        **  NUMBER CONSTANTS                                         **
!        ***************************************************************
         real (dp), parameter :: zero  = 0.00000000000000000000d0
         real (dp), parameter :: one   = 1.00000000000000000000d0
         real (dp), parameter :: three = 3.00000000000000000000d0
         real (dp), parameter :: four  = 4.00000000000000000000d0
         real (dp), parameter :: five  = 5.00000000000000000000d0
         real (dp), parameter :: six   = 6.00000000000000000000d0
         real (dp), parameter :: seven = 7.00000000000000000000d0
         real (dp), parameter :: ten   = 1.00000000000000000000d1
         real (dp), parameter :: tenth = 0.10000000000000000000d0
         real (dp), parameter :: third = one/three
         real (dp), parameter :: sixth = one/six
         real (dp), parameter :: sq3   = dsqrt(three)
         real (dp), parameter :: sq6   = dsqrt(six  )

!        ***************************************************************
!        **  MACHINE ROUNDING                                         **
!        ***************************************************************

         real (dp), parameter :: small = 1.d-19



!        ***************************************************************
!        **  PARAMETER VALUES                                         **
!        ***************************************************************

         character(len=10), parameter :: solname = "ICSPLODE"

         real (dp), parameter :: c1    = tenth * (four - sq6)
         real (dp), parameter :: c2    = tenth * (four + sq6)
         real (dp), parameter :: c1m1  = c1 - one
         real (dp), parameter :: c2m1  = c2 - one
         real (dp), parameter :: c1mc2 = c1 - c2
         real (dp), parameter :: dd1   = - third * (seven * sq6 + 13.d0)
         real (dp), parameter :: dd2   =   third * (seven * sq6 - 13.d0)
         real (dp), parameter :: dd3   = - third
         real (dp), parameter :: u1    = one/(third * tenth * (6.0d0 + &
                                         81.d0**third - 9.d0**third))
         real (dp), parameter :: alp   = sixth * tenth * (12.d0 -      &
                                         81.d0**third + 9.d0**third)
         real (dp), parameter :: bet   = sixth * tenth * sq3 * (       &
                                         81.d0**third + 9.d0**third)
         real (dp), parameter :: cno   = alp**2 + bet**2
         real (dp), parameter :: alph  = alp/cno
         real (dp), parameter :: beta  = bet/cno
         real (dp), parameter :: T11   =  9.12323948708929427920D-02
         real (dp), parameter :: T12   = -0.14125529502095420843D0
         real (dp), parameter :: T13   = -3.00291941051474244920D-02
         real (dp), parameter :: T21   =  0.24171793270710701896D0
         real (dp), parameter :: T22   =  0.20412935229379993199D0
         real (dp), parameter :: T23   =  0.38294211275726193779D0
         real (dp), parameter :: T31   =  0.96604818261509293619D0
         real (dp), parameter :: TI11  =  4.32557989006315535100D0
         real (dp), parameter :: TI12  =  0.33919925181580986954D0
         real (dp), parameter :: TI13  =  0.54177053993587487119D0
         real (dp), parameter :: TI21  = -4.17871859155190472730D0
         real (dp), parameter :: TI22  = -0.32768282076106238708D0
         real (dp), parameter :: TI23  =  0.47662355450055045196D0
         real (dp), parameter :: TI31  = -0.50287263494578687595D0
         real (dp), parameter :: TI32  =  2.57192694985560542920D0
         real (dp), parameter :: TI33  = -0.59603920482822492497D0


         contains

!        **************************************************************
!        **  MAIN SOLVER ROUTINE                                     **
!        **************************************************************
         subroutine ode_solve
         implicit none

         end subroutine ode_solve


!        **************************************************************
!        **  READ PROBLEM SETUP                                      **
!        **************************************************************
         subroutine solver_setup
         implicit none

         end subroutine solver_setup



!        **************************************************************
!        **  ICSPLODE PROBLEM INITIALISATION                         **
!        **************************************************************
         subroutine solver_init
         implicit none

         character(len=*), parameter ::                               &
           fmt_erneg = "(1x,A10,': negative integration interval not',&
                         ' allowed in this version. Exiting')"


!          Number of equations-related parameters
           neq2 = neq  * neq
           neq3 = neq2 * neq
           neq4 = neq3 * neq

!          Initialize Jacobian matrix allocation
           if (sparse_jac) then

!             Compute Jacobian matrix structure via finite diffs
!             TO BE DONEEEEEEEEEEEE!!!!!!!!!!!!!!!!
!             Allocate Jacobian sparse matrix representation
!             TO BE DONEEEEEEEEEEEE!!!!!!!!!!!!!!!!

           else
              if (allocated(JAC)) deallocate(JAC)
              allocate(JAC(neq,neq))
           end if

!          Set initial time
           t = tstart

!          Check integration interval
           if (tend-tstart<0) call ode_exit(1)

!          Interval to be integrated
           ttoend = tend - t

!          Minimum and maximum timesteps
           hmin = ten * small

!          Possible extension: user-defined maximum timestep
           hmax = ttoend

         end subroutine solver_init

!        **************************************************************
!        **  CHOOSE INITIAL STEPSIZE                                 **
!        **************************************************************
         subroutine step_init
         implicit none

!           Initialize first timestep
            dt%n        = 1
            dt%accepted = .true.
            dt%isfirst  = .true.
            dt%islast   = .false.

            dt%h   = min(hmax, ttoend)

!           Check if it's the last one
            if (dt%h * 1.000001 > ttoend) then
                dt%h      = ttoend
                dt%islast = .true.
            endif

!           Store tentative step information
            dtold = dt


         end subroutine step_init

!        **************************************************************
!        **  STEP SIZE SELECTION                                     **
!        **************************************************************
         subroutine step_size
         implicit none


         end subroutine step_size

!        **************************************************************
!        **  COMPUTE INTEGRATION STEP ATTEMPT                        **
!        **************************************************************
         subroutine step_integrate
         implicit none

!           Inserire calcolo dello Jacobiano



         end subroutine step_integrate

!        **************************************************************
!        **  ERROR ESTIMATION                                        **
!        **************************************************************
         subroutine error_estimate
         implicit none


         end subroutine error_estimate


!        **************************************************************
!        **  JACOBIAN MATRIX EVALUATION INTERFACE                    **
!        **************************************************************
         subroutine jacobian(odeF, jacF, t, Y, dYdt0, RPAR, IPAR)
         implicit none

         real (dp),                         intent(in)    :: t
         real (dp), dimension(neq),         intent(inout) :: Y
         real (dp), dimension(neq),         intent(in)    :: dYdt0
         real (dp), dimension(:), optional, intent(in)    :: rpar
         integer,   dimension(:), optional, intent(in)    :: ipar
         external                                         :: odeF,jacF

           if (.not. analytical_jac) then
              call jac_finite_diff(odeF, t, Y, dYdt0, RPAR, IPAR)
           else
              if (.not. sparse_jac) then
                 call jac_analytical(jacF, t, Y, RPAR, IPAR)
              else
                 call jac_sp_analytical(jacF, t, Y, RPAR, IPAR)
              endif
           end if

         end subroutine jacobian



!        **************************************************************
!        **  JACOBIAN MATRIX COMPUTATION THROUGH FINITE DIFFERENCES  **
!        **************************************************************
         subroutine jac_finite_diff(odeF, t, Y, dYdt0, RPAR, IPAR)
         implicit none


         real (dp),                         intent(in)    :: t
         real (dp), dimension(neq),         intent(inout) :: Y
         real (dp), dimension(neq),         intent(in)    :: dYdt0
         real (dp), dimension(:), optional, intent(in)    :: rpar
         integer,   dimension(:), optional, intent(in)    :: ipar
!         real (dp), dimension(neq,neq)                    :: jac
!         external                                         :: odeF

         real (dp)                 :: deltaY, Ybk
         real (dp), dimension(neq) :: dYdt
         integer                   :: i


            perturbate: do i = 1, neq

!              Back-up current value of Y(i)
               Ybk    = Y(i)

!              Evaluate finite increment of i-th variable
               deltaY = max(small, sqrt(abs(Ybk)) )
               Y(i)   = Ybk + deltaY

!              Evaluate function at current perturbated Y array
               call odeF(neq,t,Y,dYdt,RPAR,IPAR)

!              Compute finite differences and store
               jac(:, i) = (dYdt - dYdt0)/deltaY

!              Restore Y(i)
               Y(i)   = Ybk

            end do perturbate


         end subroutine jac_finite_diff

!        **************************************************************
!        **  FULL ANALYTICAL JACOBIAN COMPUTATION THROUGH jacF       **
!        **************************************************************
         subroutine jac_analytical(jacF, t, Y, RPAR, IPAR)
         implicit none

         real (dp),                         intent(in)    :: t
         real (dp), dimension(neq),         intent(inout) :: Y
         real (dp), dimension(:), optional, intent(in)    :: rpar
         integer,   dimension(:), optional, intent(in)    :: ipar
         real (dp), dimension(neq,neq)                    :: jac
!         external                                         :: jacF

!           Call to the external routine for analytical jacobian
            call jacF(neq,t,Y,jac,rpar,ipar)

         end subroutine jac_analytical

!        **************************************************************
!        **  SPARSE ANALYTICAL JACOBIAN COMPUTATION THROUGH jacF     **
!        **************************************************************
         subroutine jac_sp_analytical(jacF, t, Y, RPAR, IPAR)
         implicit none

         real (dp),                         intent(in)    :: t
         real (dp), dimension(neq),         intent(inout) :: Y
         real (dp), dimension(:), optional, intent(in)    :: rpar
         integer,   dimension(:), optional, intent(in)    :: ipar

!           Call to the external routine for analytical jacobian
            call jacF(neq,t,Y,jacs,rpar,ipar)

         end subroutine jac_sp_analytical


!        **************************************************************
!        **  EXIT ON ERROR                                           **
!        **************************************************************
         subroutine ode_exit(ier)
         implicit none

         integer, intent(in) :: ier

         character(len=*), parameter ::                               &
           fmt_erneg = "(1x,A10,': negative integration interval not',&
                         ' allowed in this version. Exiting')"

           if (ier == 1) write(*,fmt_erneg)solname


           stop


         end subroutine ode_exit



!     ** End of module ************************************************
      end module radau_sparse
