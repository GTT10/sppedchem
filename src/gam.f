       SUBROUTINE   GAM(R,FCN,T0,Y0,TEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JAC ,IJAC,MLJAC,MUJAC,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C
C     VERSION: NOVEMBER 25, 1999
C
C     PURPOSE: THE GAM CODE NUMERICALLY SOLVES A (POSSIBLY STIFF)
C              SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL
C              EQUATIONS IN THE FORM  Y'=F(X,Y), WITH A GIVEN
C              INITIAL CONDITION.
C
C     AUTHORS: F. IAVERNARO AND F. MAZZIA
C              UNIVERSITA' DEGLI STUDI DI BARI,
C              DIPARTIMENTO DI MATEMATICA
C              VIA ORABONA 4, 70125 BARI, ITALY
C              E-MAIL:  LABOR@ALPHAMATH.DM.UNIBA.IT
C
C     METHODS: THE METHODS USED ARE IN THE CLASS OF BOUNDARY VALUE
C              METHODS (BVMs), NAMELY THE GENERALIZED ADAMS METHODS
C              (GAMs) OF ORDER 3-5-7-9 WITH STEP SIZE CONTROL.
C
C  REFERENCES: L.BRUGNANO, D.TRIGIANTE,  Solving Differential Problems
C              by Multistep Initial and Boundary Value Methods,
C              Gordon & Breach.
C
C              F.IAVERNARO, F.MAZZIA,  Block-Boundary Value Methods
C              for the solution of Ordinary Differential Equations,
C              (submitted).
C
C              F.IAVERNARO, F.MAZZIA,  Solving Ordinary Differential
C              Equations by Generalized Adams Methods: properties and
C              implementation techniques,
C              proceedings of NUMDIFF8, Appl. Numer. Math. (to appear).
C
C DESCRIPTION: THE GAM CODE CONSISTS OF THREE FILES:
C            - gam.f CONTAINS THE MAIN SUBROUTINES THAT IMPLEMENT THE
C              INTEGRATION PROCEDURE;
C            - subgam.f  CONTAINS THE ADDITIONAL LINEAR ALGEBRA ROUTINES
C              REQUIRED BY GAM.F PLUS SOME OTHER SUBROUTINES PROPER OF
C              THE USED METHODS;
C            - param.dat CONTAINS THE DEFINITION OF ALL THE COSTANTS.
C              THIS FILE MUST BE INSERTED IN THE SAME DIRECTORY OF THE MAIN 
C              PROGRAM; IT IS INCLUDED IN ALL THE SUBROUTINE WITH THE
C              FORTRAN 77 INSTRUCTION
C              INCLUDE "param.dat" .
C
C    COMMENTS: THE PHILOSOFY AND THE STYLE USED IN WRITING THE CODE ARE VERY
C              SIMILAR TO THOSE CHARACTERIZING THE CODE  RADAU5.
C              INDEED THE AUTHORS IMPORTED FROM RADAU5 SOME SUBROUTINES,
C              COMMENTS AND IMPLEMENTATION THECNIQES LEAVING UNCHANGED
C              THE NAME AND THE MEANING OF  A NUMBER OF VARIABLES.
C              THE AUTHORS ARE VERY GRATEFUL TO ANYONE USES THE CODE AND
C              WOULD APPRECIATE ANY CRITICISM AND REMARKS ON HOW IT PERFORMS.
C
C
C
C -------------------------------------------------------------------------
C     INPUT PARAMETERS
C -------------------------------------------------------------------------
C     R           DIMENSION OF THE SYSTEM
C
C     FCN         NAME (EXTERNAL) OF THE SUBROUTINE COMPUTING THE
C                 VALUE OF F(T,Y):
C                    SUBROUTINE FCN(R,T,Y,F,RPAR,IPAR)
C                    DOUBLE PRECISION X,Y(R),F(R)
C                    F(1)=...   ETC.
C                 (RPAR, IPAR    SEE BELOW)
C
C     T0          INITIAL T-VALUE
C
C     Y0          INITIAL VALUES FOR Y
C
C     TEND        FINAL T-VALUE (TEND-T MUST BE POSITIVE)
C
C     H           INITIAL STEP SIZE GUESS;
C                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT,
C                 H=1.D0/(NORM OF F'), USUALLY 1.D-3 OR 1.D-5, IS GOOD.
C                 THIS CHOICE IS NOT VERY IMPORTANT, THE STEP SIZE IS
C                 QUICKLY ADAPTED. (IF H=0.D0, THE CODE PUTS H=1.D-6).
C
C     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
C                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH R.
C
C     ITOL        SWITCH FOR RTOL AND ATOL:
C                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
C                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
C                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL
C                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
C                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
C                     RTOL(I)*ABS(Y(I))+ATOL(I).
C
C     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
C                 THE PARTIAL DERIVATIVES OF F(T,Y) WITH RESPECT TO Y
C                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1; SUPPLY
C                 A DUMMY SUBROUTINE IN THE CASE IJAC=0).
C                 FOR IJAC=1, THIS SUBROUTINE MUST HAVE THE FORM
C                    SUBROUTINE JAC(R,T,Y,DFY,LDFY,RPAR,IPAR)
C                    DOUBLE PRECISION T,Y(R),DFY(LDFY,R)
C                    DFY(1,1)= ...
C                 LDFY, THE COLUMN-LENGTH OF THE ARRAY, IS
C                 FURNISHED BY THE CALLING PROGRAM.
C                 IF (MLJAC.EQ.R) THE JACOBIAN IS SUPPOSED TO
C                    BE FULL AND THE PARTIAL DERIVATIVES ARE
C                    STORED IN DFY AS
C                       DFY(I,J) = PARTIAL F(I) / PARTIAL Y(J)
C                 ELSE, THE JACOBIAN IS TAKEN AS BANDED AND
C                    THE PARTIAL DERIVATIVES ARE STORED
C                    DIAGONAL-WISE AS
C                       DFY(I-J+MUJAC+1,J) = PARTIAL F(I) / PARTIAL Y(J).
C
C     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN:
C                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE
C                       DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED.
C                    IJAC=1: JACOBIAN IS SUPPLIED BY SUBROUTINE JAC.
C
C     MLJAC       SWITCH FOR THE BANDED STRUCTURE OF THE JACOBIAN:
C                    MLJAC=R: JACOBIAN IS A FULL MATRIX. THE LINEAR
C                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
C                       0<=MLJAC<R: MLJAC IS THE LOWER BANDWITH OF JACOBIAN
C                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
C                       THE MAIN DIAGONAL).
C
C     MUJAC       UPPER BANDWITH OF JACOBIAN  MATRIX (>= NUMBER OF NON-
C                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
C                 NEED NOT BE DEFINED IF MLJAC=R.
C
C
C     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
C                 NUMERICAL SOLUTION DURING INTEGRATION.
C                 IF IOUT=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
C                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0.
C                 IT MUST HAVE THE FORM
C                    SUBROUTINE SOLOUT(R,TP,YP,FF,NT,DBLK,ORD,RPAR,IPAR,IRTRN)
C                    INTEGER R, DBLK, ORD, IPAR(*), IRTRN, NT
C                    DOUBLE PRECISION TP(*), YP(R,*), RPAR(*), FF(R,*)
C                    ....
C                 SOLOUT FURNISHES THE SOLUTION "YP" AT THE
C                    GRID-POINTS "TP(*)".
C                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
C                    IS SET <0, GAM  RETURNS TO THE CALLING PROGRAM.
C
C                 CONTINUOUS OUTPUT:
C                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
C                 FOR THE INTERVAL [TP(1),TP(DBLK+1)] IS AVAILABLE THROUGH
C                 THE FUNCTION
C                        >>>   CONTR(I,R,T,TP,FF,DBLK,NT)   <<<
C                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH
C                 COMPONENT OF THE SOLUTION AT THE POINT T. THE VALUE
C                 T SHOULD LIE IN THE INTERVAL [TP(1),TP(DBLK+1)] ON
C                 WHICH THE SOLUTION IS COMPUTED AT CURRENT STEP.
C                 DO NOT CHANGE THE ENTRIES OF FF(R,*) and NT, IF THE
C                 DENSE OUTPUT FUNCTION IS USED.
C
C     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT:
C                    IOUT=0: SUBROUTINE IS NEVER CALLED
C                    IOUT=1: SUBROUTINE IS AVAILABLE FOR OUTPUT.
C
C     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
C                 WORK(1), WORK(2),.., WORK(19) SERVE AS PARAMETERS
C                 FOR THE CODE. FOR STANDARD USE OF THE CODE
C                 WORK(1),..,WORK(18) MUST BE SET TO ZERO BEFORE
C                 CALLING. SEE BELOW FOR A MORE SOPHISTICATED USE.
C                 WORK(19),..,WORK(LWORK) SERVE AS WORKING SPACE
C                 FOR ALL VECTORS AND MATRICES.
C                 "LWORK" MUST BE AT LEAST:
C
C                 IN THE USUAL CASE WHERE THE JACOBIAN IS FULL
C                 STORAGE REQUIREMENT IS
C                             LWORK = 2*R*R+42*R+18.
C
C                 IN THE CASE WHERE THE JACOBIAN IS SPARSE
C                 STORAGE REQUIREMENT IS
C                             LWORK = (3*MLJAC+2*MUJAC+40)*R+18
C
C
C     LWORK       DECLARED LENGTH OF ARRAY "WORK".
C
C     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK".
C                 IWORK(1),IWORK(2),...,IWORK(24) SERVE AS PARAMETERS
C                 FOR THE CODE. FOR STANDARD USE, SET IWORK(1),..,
C                 IWORK(8) TO ZERO BEFORE CALLING.
C                 IWORK(10),...,IWORK(LIWORK) SERVE AS WORKING AREA.
C                 "LIWORK" MUST BE AT LEAST
C
C                          LIWORK = 24 + R
C
C     LIWORK      DECLARED LENGTH OF ARRAY "IWORK".
C
C     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH
C                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING
C                 PROGRAM AND THE FCN, JAC SUBROUTINES.
C
C     SEE THE EXAMPLE BELOW FOR A PRACTICAL EXPLANATION ON THE USE OF
C     SOME OF THE LISTED VARIABLES AND SUBROUTINES.
C
C -------------------------------------------------------------------------
C     EXAMPLE PROBLEM.
C---------------------------------------------------------------------------
c
c the following is a simple example problem, with the coding
c needed for its solution by GAM.  the problem is from chemical
c kinetics, and consists of the following three rate equations..
c     dy1/dt = -.04*y1 + 1.e4*y2*y3
c     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
c     dy3/dt = 3.e7*y2**2
c on the interval from t = 0.0 to t = 4.e10, with initial conditions
c y1 = 1.0, y2 = y3 = 0.  the problem is stiff.
c
c the following coding solves this problem with GAM,
c printing results at t = .4, 4., ..., 4.e10.  it uses
c itol = 1 and atol much smaller for y2 than y1 or y3 because
c y2 has much smaller values.
c at the end of the run, statistical quantities of interest are
c printed (see optional outputs in the full description below).
c
c      implicit none
c      integer     neq, lwork, liwork
c      parameter ( neq=3, lwork=2*neq*neq+42*neq+18, liwork=neq+24)
c      double precision atol(neq), rtol(neq),work(lwork),t,tout,y(neq)
c      double precision h, rpar
c      integer  i, iwork(liwork), itol, iout, nsteps, naccept
c      integer  mljac, mujac, ijac, ipar, idid
c      external feval, jeval, solout
c
c      y(1) = 1.0d0
c      y(2) = 0.0d0
c      y(3) = 0.0d0
c      t = 0.0d0
c      tout = 4.0d10
c      iout = 1
c      ijac  = 0
c      mljac = neq
c      mujac = neq
c      h = 1d-6
c
c      do i = 1,18
c         work(i) = 0.0d0
c      enddo
c      do i = 1,24
c         iwork(i) = 0.0d0
c      enddo
c
c      itol = 1
c      rtol(1) = 1.0d-5
c      rtol(2) = 1.0d-5
c      rtol(3) = 1.0d-5
c      atol(1) = 1.0d-5
c      atol(2) = 1.0d-8
c      atol(3) = 1.0d-5
c
c      rpar = 0.4d0
c      call GAM(neq,feval,t,y,tout,h,rtol,atol,itol,
c     &                  jeval ,ijac,mljac,mujac,solout,iout,
c     &                  work,lwork,iwork,liwork,rpar,ipar,idid)
c
c 50         format(7h at t =,e12.4,6h   y =,3e14.6)
c      write(6,*)
c      write(6,50) t, y(1), y(2), y(3)
c
c      nsteps = 0
c      do i=12,23
c         nsteps = nsteps + iwork(i)
c      end do
c      naccept = iwork(12)+iwork(13)+iwork(14)+iwork(15)
c      write(6,41) nsteps,naccept,iwork(10),iwork(11),iwork(24)
c   41 format(  ' # steps  =         ',i8,/,
c     +         ' # accept =         ',i8,/,
c     +         ' # f-eval =         ',i8,/,
c     +         ' # Jac-eval =       ',i8,/,
c     +         ' # LU-decomp =      ',i8/)
c
c      stop
c      end
c
c      subroutine feval(neq, t, y, ydot, rpar, ipar)
c      double precision t, y(3), ydot(3), rpar
c      integer ipar
c      ydot(1) = -.04d0*y(1) + 1.0d4*y(2)*y(3)
c      ydot(3) = 3.0d7*y(2)*y(2)
c      ydot(2) = -ydot(1) - ydot(3)
c      return
c      end
c      subroutine jeval(neqn,t,y,jac,ldim,rpar,ipar)
c      double precision t,y(neqn),jac(ldim,neqn),rpar
c      integer neqn,ldim,ipar
c      return
c      end
c
c      subroutine solout(r,tp,yp,ff,nt1,dblk,ord,rpar,ipar,irtrn)
c      implicit none
c      integer r,dblk,ord,ipar(*),irtrn,nt1
c      double precision tp(*),yp(r,*),rpar(*),ff(r,*),y(3),contr,t
c             t = rpar(1)
c             if ( (tp(1).le.t).and.(t.lt.tp(dblk+1)) ) then
c               y(1) = contr(1,r,t,tp,ff,dblk,nt1)
c               y(2) = contr(2,r,t,tp,ff,dblk,nt1)
c               y(3) = contr(3,r,t,tp,ff,dblk,nt1)
c               write(6,50) t, y(1), y(2), y(3)
c               rpar(1) = rpar(1)*10d0
c             endif
c 50         format(7h at t =,e12.4,6h   y =,3e14.6)
c      return
c      end
c
C -------------------------------------------------------------------------
C     SOPHISTICATED SETTING OF PARAMETERS
C -------------------------------------------------------------------------
C              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK
C              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),...
C              AS WELL AS IWORK(1),... DIFFERENT FROM ZERO.
C              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES:
C
C    IWORK(1)  NOT USED
C
C    IWORK(2)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
C              THE DEFAULT VALUE (FOR IWORK(2)=0) IS 100000.
C
C    IWORK(3)  ORDMIN, 3 <= ORDMIN <= 9,
C
C    IWORK(4)  ORDMAX, ORDMIN <= ORDMAX <= 9
C
C    IWORK(5)  THE MAXIMUM NUMBER OF SPLITTING-NEWTON ITERATIONS FOR THE
C              SOLUTION OF THE IMPLICIT SYSTEM IN EACH STEP FOR ORDER 3.
C              THE DEFAULT VALUE (FOR IWORK(5)=0) IS 10.
C
C    IWORK(6)  THE MAXIMUM NUMBER OF SPLITTING-NEWTON ITERATION FOR
C              ORDER 5, THE DEFAULT VALUE (FOR IWORK(6)=0) IS 18.
C
C    IWORK(7)  THE MAXIMUM NUMBER OF SPLITTING-NEWTON ITERATION FOR
C              ORDER 7, THE DEFAULT VALUE (FOR IWORK(7)=0) IS 26.
C
C    IWORK(8)  THE MAXIMUM NUMBER OF SPLITTING-NEWTON ITERATION FOR
C              ORDER 9, THE DEFAULT VALUE (FOR IWORK(5)=0) IS 36.
C
C
C    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.
C
C    WORK(2)   HMAX  MAXIMAL STEP SIZE, DEFAULT TEND-T0.
C
C    WORK(3)   THET DECIDE WHETHER THE JACOBIAN SHOULD BE RECOMPUTED
C
C    WORK(4)   FACNEWT:  stopping criterion for splitting-Newton method
C                 for small values of min(abs(y_i)) and min(abs(f_j)).
C
C    WORK(5)   TETAK0 stopping criterium for the splitting-Newton method
C                the iterates must be decreasing by a factor tetak0
C
C    WORK(6)   CS(2): EMPIRICAL COMPUTATIONAL COST FOR ORDER  5 METHOD
C              USED IN THE ORDER VARIATION STRATEGY
C              (DEFAULT WORK(6) = 2.4D0)
C
C    WORK(7)   CS(3): EMPIRICAL COMPUTATIONAL COST FOR ORDER  7 METHOD
C              USED IN THE ORDER VARIATION STRATEGY
C              (DEFAULT WORK(6) = 4.0D0)
C
C    WORK(8)   CS(4): EMPIRICAL COMPUTATIONAL COST FOR ORDER  9 METHOD
C              USED IN THE ORDER VARIATION STRATEGY
C              (DEFAULT WORK(6) = 7.2D0)
C
C    WORK(9)-WORK(10)   FACL-FACR: PARAMETERS FOR STEP SIZE SELECTION
C               THE NEW STEPSIZE IS CHOSEN SUBJECT TO THE RESTRICTION
C               FACL  <=  HNEW/HOLD <= FACR
C               (DEFAULT WORK(9) = 0.12, WORK(10) = 10 )
C
C    WORK(11)  SFDOWN:SAFETY FACTOR IN STEP SIZE PREDICTION
C                  USED FOR THE LOWER ORDER METHOD
C                  (DEFAULT WORK(11) = 20D0)
C
C    WORK(12)  SFUP:SAFETY FACTOR IN STEP SIZE PREDICTION
C                  USED FOR THE UPPER ORDER METHOD
C                  (DEFAULT WORK(12) = 40D0)
C
C    WORK(13)  SFSAME: SAFETY FACTOR IN STEP SIZE PREDICTION
C                  USED FOR THE CURRENT ORDER METHOD
C                  (DEFAULT WORK(13) = 18D0)
C
C    WORK(14)  SF: SAFETY FACTOR IN STEP SIZE PREDICTION
C                  USED FOR THE CURRENT ORDER METHOD WHEN IS
C                  FAILED THE ERROR CONTROL TEST (DEFAULT WORK(14) = 15D0)
C
C
C
C    WORK(15)  FACNEWT stopping criterion for splitting-Newton method ORDER 3
C                  (DEFAULT WORK(15) = 1.0D-3)
C
C    WORK(16)  FACNEWT stopping criterion for splitting-Newton method ORDER 5
C                  (DEFAULT WORK(16) = 9.0D-2)
C
C    WORK(17)  FACNEWT stopping criterion for splitting-Newton method ORDER 7
C                  (DEFAULT WORK(17) = 9.0D-1)
C
C    WORK(18)  FACNEWT stopping criterion for splitting-Newton method ORDER 9
C                  (DEFAULT WORK(18) = 9.9D-1)
C
C -------------------------------------------------------------------------
C     OUTPUT PARAMETERS
C -------------------------------------------------------------------------
C     T0          T-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED
C                 (AFTER SUCCESSFUL RETURN T0=TEND).
C
C     Y(N)        NUMERICAL SOLUTION AT T0
C
C     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
C
C     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
C                   IDID= 1  COMPUTATION SUCCESSFUL,
C                   IDID=-1  INPUT IS NOT CONSISTENT,
C                   IDID=-2  LARGER NMAX IS NEEDED,
C                   IDID=-3  STEP SIZE BECOMES TOO SMALL,
C                   IDID=-4  MATRIX IS REPEATEDLY SINGULAR.
C
C   IWORK(10)  NFCN    NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL
C                      EVALUATION OF THE JACOBIAN ARE NOT COUNTED)
C   IWORK(11)  NJAC    NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY
C                      OR NUMERICALLY)
C   IWORK(12)  NSTEP(1)  NUMBER OF COMPUTED STEPS   ORD 3
C   IWORK(13)  NSTEP(2)  NUMBER OF COMPUTED STEPS   ORD 5
C   IWORK(14)  NSTEP(3)  NUMBER OF COMPUTED STEPS   ORD 7
C   IWORK(15)  NSTEP(4)  NUMBER OF COMPUTED STEPS   ORD 9
C   IWORK(16)  NNEWT(1)  NUMBER OF REJECTED STEPS (DUE TO NEWTON CONVERGENCE) 3
C   IWORK(17)  NNEWT(2)  NUMBER OF REJECTED STEPS (DUE TO NEWTON CONVERGENCE) 5
C   IWORK(18)  NNEWT(3)  NUMBER OF REJECTED STEPS (DUE TO NEWTON CONVERGENCE) 7
C   IWORK(19)  NNEWT(4)  NUMBER OF REJECTED STEPS (DUE TO NEWTON CONVERGENCE) 9
C   IWORK(20)  NERR(1)   NUMBER OF REJECTED STEPS (DUE TO ERROR TEST) 3
C   IWORK(21)  NERR(2)   NUMBER OF REJECTED STEPS (DUE TO ERROR TEST) 5
C   IWORK(22)  NERR(3)   NUMBER OF REJECTED STEPS (DUE TO ERROR TEST) 7
C   IWORK(23)  NERR(4)   NUMBER OF REJECTED STEPS (DUE TO ERROR TEST) 9
C   IWORK(24)  NDEC      NUMBER OF LU-DECOMPOSITIONS
C-----------------------------------------------------------------------
C     DECLARATIONS
C -------------------------------------------------------------------------
       IMPLICIT NONE
C
C   INPUT VARIABLES
C------------------------------------
       INTEGER R, ORDMIN, ORDMAX, IDID, ITOL, IJAC, MLJAC, MUJAC, IOUT,
     &         IPAR(*), ITINT(4), ITMAX, IJOB, NMAX, LDJAC, LDLU
   
       DOUBLE PRECISION TEND, ATOL(*), RTOL(*), RPAR(*), FACNORD(4),
     &                  HMAX, THET, FACNEWT, TETAK0, CS(4), FACL, FACR,
     &                  SFDOWN, SFUP, SFSAME, SF, UROUND

C
C   OUTPUT VARIABLES
C------------------------------------
       INTEGER NDEC, NFCN, NJAC, NSTEP(4), NNEWT(4), NERR(4)
C
C   INPUT/OUTPUT VARIABLES
C------------------------------------
       INTEGER  LIWORK, LWORK, IWORK(LIWORK)
       DOUBLE PRECISION T0, Y0(R), H, WORK(LWORK)

C
C   LOCAL VARIABLES
C------------------------------------
       INTEGER  IEYP, IEFP, IEDN,IEF,IEF1, IEJF0,IELU,ISTORE,
     &          I, IEIPIV, IESC
       LOGICAL  ARRET, JBAND

C
C   EXTERNAL FUNCTIONS
C------------------------------------
       EXTERNAL FCN,JAC, SOLOUT
C -------------------------------------------------------------------------
C     SETTING THE PARAMETERS
C -------------------------------------------------------------------------
      NFCN    =0
      NJAC    =0
      NSTEP(1)=0
      NSTEP(2)=0
      NSTEP(3)=0
      NSTEP(4)=0
      NNEWT(1)=0
      NNEWT(2)=0
      NNEWT(3)=0
      NNEWT(4)=0
      NERR(1) =0
      NERR(2) =0
      NERR(3) =0
      NERR(4) =0
      NDEC    =0
      ARRET   = .FALSE.
C -------- NMAX := THE MAXIMAL NUMBER OF STEPS -----
      IF (IWORK(2).EQ.0) THEN
         NMAX=100000
      ELSE
         NMAX=IWORK(2)
         IF (NMAX.LE.0) THEN
            WRITE(6,*)' WRONG INPUT IWORK(2)=',IWORK(2)
            ARRET=.TRUE.
         END IF
      END IF
C--------- ORDMIN  :=  MINIMAL ORDER
      IF (IWORK(3).EQ.0) THEN
         ORDMIN = 1
      ELSE
         ORDMIN=(IWORK(3)-1)/2
         IF ((ORDMIN.LE.0).OR.(ORDMIN.GT.4)) THEN
            WRITE(6,*)' CURIOUS INPUT IWORK(3)=',IWORK(3)
            ARRET=.TRUE.
         END IF
      END IF
C--------- ORDMAX :=  MAXIMAL ORDER
      IF (IWORK(4).EQ.0) THEN
         ORDMAX = 4
      ELSE
         ORDMAX=(IWORK(4)-1)/2
         IF ((ORDMAX.LE.0).OR.(ORDMAX.GT.4).OR.(ORDMAX.LT.ORDMIN)) THEN
            WRITE(6,*)' CURIOUS INPUT IWORK(4)=',IWORK(4)
            ARRET=.TRUE.
         END IF
      END IF
C -------- ITINT(1) :=  NUMBER OF SPLITTING-NEWTON ITERATIONS ORD 3
      IF (IWORK(5).EQ.0) THEN
         ITINT(1)=12
      ELSE
         ITINT(1)=IWORK(5)
         IF (ITINT(1).LE.0) THEN
            WRITE(6,*)' CURIOUS INPUT IWORK(5)=',IWORK(5)
            ARRET=.TRUE.
         END IF
      END IF
      ITMAX = ITINT(1)
C -------- ITINT(2) :=  NUMBER OF SPLITTING-NEWTON ITERATIONS ORD 5
      IF (IWORK(6).EQ.0) THEN
         ITINT(2)=18
      ELSE
         ITINT(2)=IWORK(6)
         IF (ITINT(2).LT.0) THEN
            WRITE(6,*)' CURIOUS INPUT IWORK(6)=',IWORK(6)
            ARRET=.TRUE.
         END IF
      END IF
C -------- ITINT(3) :=  NUMBER OF SPLITTING-NEWTON ITERATIONS ORD 7
      IF (IWORK(7).EQ.0) THEN
         ITINT(3)= 26
      ELSE
         ITINT(3)=IWORK(7)
         IF (ITINT(3).LT.0) THEN
            WRITE(6,*)' CURIOUS INPUT IWORK(8)=',IWORK(7)
            ARRET=.TRUE.
         END IF
      END IF
C -------- ITINT(4) :=  NUMBER OF SPLITTING-NEWTON ITERATIONS ORD 9
      IF (IWORK(8).EQ.0) THEN
         ITINT(4)= 36
      ELSE
         ITINT(4)=IWORK(8)
         IF (ITINT(4).LT.0) THEN
            WRITE(6,*)' CURIOUS INPUT IWORK(8)=',IWORK(8)
            ARRET=.TRUE.
         END IF
      END IF

C -------- UROUND :=  SMALLEST NUMBER SATISFYING 1.0D0+UROUND>1.0D0
      IF (WORK(1).EQ.0.0D0) THEN
         UROUND=1.0D-16
      ELSE
         UROUND=WORK(1)
         IF (UROUND.LE.1.0D-19.OR.UROUND.GE.1.0D0) THEN
            WRITE(6,*)'COEFFICIENTS HAVE 20 DIGITS, UROUND=',WORK(1)
            ARRET=.TRUE.
         END IF
      END IF
C -------- HMAX := MAXIMAL STEP SIZE
      IF (WORK(2).EQ.0.D0) THEN
         HMAX=TEND-T0
      ELSE
         HMAX=WORK(2)
         IF (HMAX.GT.TEND-T0) THEN
            HMAX=TEND-T0
         END IF
      END IF
C -------- THET  DECIDE WHETHER THE JACOBIAN SHOULD BE RECOMPUTED
      IF (WORK(3).EQ.0.D0) THEN
         THET = 0.005
      ELSE
         THET=WORK(3)
         IF (THET .GT. 1d0) THEN
            WRITE(6,*)' CURIOUS INPUT WORK(3)=',WORK(3)
            ARRET=.TRUE.
         END IF
      END IF
C--------- FACNEWT: STOPPING CRITERION FOR SPLITTING-NEWTON METHOD
C--------           FOR SMALL VALUES OF min(abs(y_i)) and min(abs(f_j))
      IF (WORK(4).EQ.0.D0) THEN
          FACNEWT=5d-3
          FACNEWT=DMAX1(FACNEWT,UROUND/RTOL(1) )
      ELSE
         FACNEWT=WORK(4)
         FACNEWT=DMAX1(FACNEWT,UROUND/RTOL(1) )
         IF (FACNEWT.GE.1.0D0) THEN
            WRITE(6,*)'WRONG INPUT FOR WORK(4) ',WORK(4)
            ARRET=.TRUE.
         END IF
      END IF
C--------- FACNORD(1): STOPPING CRITERION FOR SPLITTING-NEWTON METHOD
C---------             ORDER 3
      IF (WORK(15).EQ.0.D0) THEN
          FACNORD(1) = 1d-3
          FACNORD(1) = DMAX1(FACNORD(1) ,UROUND/RTOL(1) )
      ELSE
         FACNORD(1) = WORK(15)
         FACNORD(1) = DMAX1(FACNORD(1) ,UROUND/RTOL(1) )
         IF (FACNEWT.GE.1.0D0) THEN
            WRITE(6,*)'WRONG INPUT FOR WORK(15) ',WORK(15)
            ARRET=.TRUE.
         END IF
      END IF
C--------- FACNORD(2): STOPPING CRITERION FOR SPLITTING-NEWTON METHOD
C--------              ORDER 5
      IF (WORK(16).EQ.0.D0) THEN
          FACNORD(2) = 9d-2
          FACNORD(2) = DMAX1(FACNORD(2) ,UROUND/RTOL(1) )
      ELSE
         FACNORD(2) = WORK(16)
         FACNORD(2) = DMAX1(FACNORD(2) ,UROUND/RTOL(1) )
         IF (FACNORD(2).GE.1.0D0) THEN
            WRITE(6,*)'WRONG INPUT FOR WORK(16) ',WORK(16)
            ARRET=.TRUE.
         END IF
      END IF
C--------- FACNORD(3): STOPPING CRITERION FOR SPLITTING-NEWTON METHOD
C---------             ORDER 7
      IF (WORK(17).EQ.0.D0) THEN
          FACNORD(3) = 9d-1
          FACNORD(3) = DMAX1(FACNORD(3), UROUND/RTOL(1) )
      ELSE
         FACNORD(3) = WORK(17)
         FACNORD(3) = DMAX1(FACNORD(3), UROUND/RTOL(1) )
         IF (FACNORD(3).GE.1.0D0) THEN
            WRITE(6,*)'WRONG INPUT FOR WORK(17) ',WORK(17)
            ARRET=.TRUE.
         END IF
      END IF
C--------- FACNORD(4): STOPPING CRITERION FOR SPLITTING-NEWTON METHOD
C---------             ORDER 9
      IF (WORK(18).EQ.0.D0) THEN
          FACNORD(4) = 9.9d-1
          FACNORD(4) = DMAX1(FACNORD(4),UROUND/RTOL(1) )
      ELSE
         FACNORD(4) = WORK(18)
         FACNORD(4) = DMAX1(FACNORD(4),UROUND/RTOL(1) )
         IF (FACNORD(4).GE.1.0D0) THEN
            WRITE(6,*)'WRONG INPUT FOR WORK(18) ',WORK(18)
            ARRET=.TRUE.
         END IF
      END IF

C--------- TETAK0: STOPPING CRITERIUM FOR THE SPLITTING-NEWTON METHOD
C---------         THE ERROR IN THE ITERATES MUST BE DECREASING
C---------         BY A FACTOR TETAK0
      IF (WORK(5).EQ.0.D0) THEN
         TETAK0 = 0.9D0
      ELSE
         TETAK0 = WORK(5)
         IF (TETAK0.LE.0.0D0) THEN
            WRITE(6,*)'WRONG INPUT FOR WORK(5) ',WORK(5)
            ARRET=.TRUE.
         END IF
      END IF
      CS(1) = 1.0D0
C--------- CS(2): EMPIRICAL COMPUTATIONAL COST FOR ORDER 5
C---------        USED IN THE ORDER VARIATION STRATEGY.
      IF (WORK(6).EQ.0.D0) THEN
         CS(2) = 2.4D0
      ELSE
         CS(2) = WORK(6)
         IF (CS(2).LE.0.0D0) THEN
            WRITE(6,*)'WRONG INPUT FOR WORK(6) ',WORK(6)
            ARRET=.TRUE.
         END IF
      END IF
C---------  CS(3): EMPIRICAL COMPUTATIONAL COST FOR ORDER 7
C---------         USED IN THE ORDER VARIATION STRATEGY.
      IF (WORK(7).EQ.0.D0) THEN
         CS(3) = 4.0D0
      ELSE
         CS(3) = WORK(7)
         IF (CS(3).LE.0.0D0) THEN
            WRITE(6,*)'WRONG INPUT FOR WORK(7) ',WORK(7)
            ARRET=.TRUE.
         END IF
      END IF
C--------- CS(4): EMPIRICAL COMPUTATIONAL COST FOR ORDER 9
C---------        USED IN THE ORDER VARIATION STRATEGY.
      IF (WORK(8).EQ.0.D0) THEN
         CS(4) =7.2D0
      ELSE
         CS(4) = WORK(8)
         IF (CS(4).LE.0.0D0) THEN
            WRITE(6,*)'WRONG INPUT FOR WORK(8) ',WORK(8)
            ARRET=.TRUE.
         END IF
      END IF
C--------- FACL: PARAMETER FOR STEP SIZE SELECTION
C---------       THE NEW STEPSIZE IS CHOSEN SUBJECT TO THE RESTRICTION
C---------       FACL <= HNEW/HOLD
      IF (WORK(9).EQ.0.D0) THEN
         FACL = 0.12D0
      ELSE
         FACL = WORK(9)
         IF (FACL.LE.0.0D0) THEN
            WRITE(6,*)'WRONG INPUT FOR WORK(9) ',WORK(9)
            ARRET=.TRUE.
         END IF
      END IF
C--------- FACR: PARAMETER FOR STEP SIZE SELECTION
C---------       THE NEW STEPSIZE IS CHOSEN SUBJECT TO THE RESTRICTION
C---------       HNEW/HOLD <= FACR
      IF (WORK(10).EQ.0.D0) THEN
         FACR = 10D0
      ELSE
         FACR = WORK(10)
         IF (FACR.LE.0.0D0) THEN
            WRITE(6,*)'WRONG INPUT FOR WORK(10) ',WORK(10)
            ARRET=.TRUE.
         END IF
      END IF
C--------- SFDOWN: SAFETY FACTOR IN STEP SIZE PREDICTION
C---------         USED FOR THE LOWER ORDER METHOD
      IF (WORK(11).EQ.0.D0) THEN
         SFDOWN = 20.0D0
      ELSE
         SFDOWN = WORK(11)
         IF (SFDOWN.LE.0.0D0) THEN
            WRITE(6,*)'WRONG INPUT FOR WORK(11) ',WORK(11)
            ARRET=.TRUE.
         END IF
      END IF
C--------- SFUP:  SAFETY FACTOR IN STEP SIZE PREDICTION
C---------        USED FOR THE UPPER ORDER METHOD
      IF (WORK(12).EQ.0.D0) THEN
         SFUP = 40.0D0
      ELSE
         SFUP = WORK(12)
         IF (SFUP.LE.0.0D0) THEN
            WRITE(6,*)'WRONG INPUT FOR WORK(12) ',WORK(12)
            ARRET=.TRUE.
         END IF
      END IF
C--------- SFSAME: SAFETY FACTOR IN STEP SIZE PREDICTION
C---------         USED FOR THE CURRENT ORDER METHOD
      IF (WORK(13).EQ.0.D0) THEN
         SFSAME = 18.0D0
      ELSE
         SFSAME = WORK(13)
         IF (SFSAME.LE.0.0D0) THEN
            WRITE(6,*)'WRONG INPUT FOR WORK(13) ',WORK(13)
            ARRET=.TRUE.
         END IF
      END IF
C--------- SF: SAFETY FACTOR IN STEP SIZE PREDICTION
C---------     USED FOR THE CURRENT ORDER METHOD WHEN
C---------     THE ERROR CONTROL TEST fails
      IF (WORK(14).EQ.0.D0) THEN
         SF = 15.0D0
      ELSE
         SF = WORK(14)
         IF (SF.LE.0.0D0) THEN
            WRITE(6,*)'WRONG INPUT FOR WORK(14) ',WORK(14)
            ARRET=.TRUE.
         END IF
      END IF
C -------- CHECK  THE TOLERANCES
      IF (ITOL.EQ.0) THEN
          IF (ATOL(1).LE.0.D0.OR.RTOL(1).LE. UROUND) THEN
              WRITE (6,*) ' TOLERANCES ARE TOO SMALL'
              ARRET=.TRUE.
          END IF
      ELSE
          DO I=1,R
          IF (ATOL(I).LE.0.D0.OR.RTOL(I).LE. UROUND) THEN
              WRITE (6,*) ' TOLERANCES(',I,') ARE TOO SMALL'
              ARRET=.TRUE.
          END IF
          END DO
      END IF
C -------------------------------------------------------------------------
C     COMPUTATION OF ARRAY ENTRIES
C -------------------------------------------------------------------------
C -------- BANDED OR NOT 
      JBAND=MLJAC.LT.R
C -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS
      IF (JBAND) THEN
        LDJAC = MLJAC+MUJAC+1
        LDLU  = MLJAC+LDJAC
      ELSE
        LDJAC = R
        LDLU  = R
      END IF
      IF (JBAND) THEN
         IJOB=2
      ELSE
         IJOB=1
      END IF
C -------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK
      IEYP  = 19
      IEFP  = IEYP  + 10*R
      IEDN  = IEFP  + 10*R
      IEF   = IEDN  + R
      IEF1  = IEF   + 10*R
      IESC  = IEF1  + 10*R
      IEJF0 = IESC  + R
      IELU  = IEJF0 + R*LDJAC
C--------- TOTAL STORAGE REQUIREMENT
      ISTORE = IELU + R*LDLU - 1
      IF(ISTORE.GT.LWORK)THEN
         WRITE(6,*)' INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=',ISTORE
         ARRET=.TRUE.
      END IF
C -------- ENTRY POINTS FOR INTEGER WORKSPACE
      IEIPIV=25
C -------- TOTAL REQUIREMENT
      ISTORE=IEIPIV+R-1
      IF (ISTORE.GT.LIWORK) THEN
         WRITE(6,*)' INSUFF. STORAGE FOR IWORK, MIN. LIWORK=',ISTORE
         ARRET=.TRUE.
      END IF

C -------- WHEN A FAIL HAS OCCURED, GAM RETURNs WITH IDID=-1
      IF (ARRET) THEN
         IDID=-1
         RETURN
      END IF

      DO I=1,4
        NSTEP(I) = 0
        NNEWT(I) = 0
        NERR(I)  = 0
      END DO
      NDEC = 0
      NFCN = 0
      NJAC = 0


C -------------------------------------------------------------------------
C     CALL TO CORE INTEGRATOR
C -------------------------------------------------------------------------
      CALL ETRO(R, FCN, T0, Y0, TEND, HMAX, H, RTOL, ATOL, ITOL,
     &   JAC, IJAC, MLJAC, MUJAC, SOLOUT, IOUT, IDID, NMAX,
     &   UROUND, THET, FACNEWT, FACNORD, TETAK0, CS, FACL, FACR, SFDOWN,
     &   SFUP, SFSAME, SF, ORDMIN, ORDMAX, ITINT, ITMAX,
     &   JBAND, IJOB, LDJAC, LDLU, WORK(IEYP), WORK(IEFP),
     &   WORK(IEDN), WORK(IEF), WORK(IEF1), WORK(IESC),
     &   WORK(IEJF0), WORK(IELU), IWORK(IEIPIV),
     &   NFCN, NJAC, NSTEP, NNEWT, NERR, NDEC, RPAR, IPAR)

      IWORK(10)= NFCN
      IWORK(11)= NJAC
      IWORK(12)= NSTEP(1)
      IWORK(13)= NSTEP(2)
      IWORK(14)= NSTEP(3)
      IWORK(15)= NSTEP(4)
      IWORK(16)= NNEWT(1)
      IWORK(17)= NNEWT(2)
      IWORK(18)= NNEWT(3)
      IWORK(19)= NNEWT(4)
      IWORK(20)= NERR(1)
      IWORK(21)= NERR(2)
      IWORK(22)= NERR(3)
      IWORK(23)= NERR(4)
      IWORK(24)= NDEC

      RETURN
      END
C
C--------- END OF SUBROUTINE GAM
C
C -------------------------------------------------------------------------
C     SUBROUTINE  ETRO (Extended trapezoidal Rules of Odd order,
C                       that is GAMs)
C -------------------------------------------------------------------------
      SUBROUTINE  ETRO(R,FCN,T0,Y0,TEND,HMAX,H,RTOL,ATOL,ITOL,
     &   JAC,IJAC,MLJAC,MUJAC,SOLOUT,IOUT,IDID,
     &   NMAX,UROUND,THET,FACNEWT,FACNORD,TETAK0,CS,FACL,FACR,SFDOWN,
     &   SFUP,SFSAME,SF, ORDMIN,ORDMAX,ITINT,ITMAX,
     &   JBAND,IJOB,LDJAC,LDLU,YP,FP,
     &   DN,F,F1,SCAL, JF0, LU, IPIV,
     &   NFCN,NJAC,NSTEP,NNEWT,NERR,NDEC,RPAR,IPAR)
C -------------------------------------------------------------------------
C     CORE INTEGRATOR FOR GAM
C     PARAMETERS SAME AS IN GAM WITH ADDED WORKSPACE
C -------------------------------------------------------------------------
C     DECLARATIONS
C ----------------------------------------------------------

      use working_precision, only: dp
      use speedchem_conV, only: constV_jac_sparse
      IMPLICIT NONE
C
C   COMMON
C------------------------------------
      COMMON/LINAL/MLLU,MULU,MDIAG
!$OMP THREADPRIVATE(/LINAL/)
C
C   INPUT VARIABLES
C------------------------------------
       INTEGER R, ORDMIN, ORDMAX, IDID, ITOL, IJAC, MLJAC, MUJAC, IOUT,
     &         IPAR(*), ITINT(4), ITMAX, IJOB, NMAX, LDJAC, LDLU
   
       DOUBLE PRECISION TEND, ATOL(*), RTOL(*), RPAR(*), FACNORD(4),
     &                  HMAX, THET, FACNEWT, TETAK0, CS(4), FACL, FACR,
     &                  SFDOWN, SFUP, SFSAME, SF, UROUND

C
C   OUTPUT VARIABLES
C------------------------------------
       INTEGER NDEC, NFCN, NJAC, NSTEP(4), NNEWT(4), NERR(4), IER
C
C   INPUT/OUTPUT VARIABLES
C---- 
       DOUBLE PRECISION T0, Y0(R), H, SCAL(R), YP(R,10), FP(R,10), 
     &                  F(R,9),DN(R), F1(R,9), JF0(LDJAC,R), LU(LDLU,R)

C
C   LOCAL VARIABLES
C------------------------------------
       INTEGER  I, J, IPIV(R), NSING,
     &          FAILNI, FAILEI, ORDOLD, ORD, ORD2, ORDN, 
     &          IT, DBL(4), DBLK, DBLKOLD,
     &          NSTEPS, IRTRN, NT1, MLLU, MULU, MDIAG

       DOUBLE PRECISION ERRV(10), TP(10), T1(10), YSAFE, DELT, 
     &                  THETA, TETAK, TETAKOLD,
     &                  HOLD, ERRNEWT, ERRNEWT0,  ERRNEWT1, ESP,
     &                  ERRUP, ERRSAME, ERRDOWN, RR, RRN, TH, THN, FACN
      LOGICAL  JBAND, CALJAC, NEWJAC, JVAI, TER, EXTRAP
C
C   EXTERNAL FUNCTIONS
C------------------------------------
       EXTERNAL FCN,JAC, SOLOUT
C
C   INCLUDE
C------------------------------------
       INCLUDE "gamparam.dat"
c -------- CONSTANTS
      MLLU=MLJAC
      MULU=MUJAC
      MDIAG=MLLU + MULU +1
C--------- DBL(1:4) := SIZE OF THE COEFFICIENT MATRICES DEFINING THE GAMs

      DBL(1) = 4
      DBL(2) = 6
      DBL(3) = 8
      DBL(4) = 9
      ORD  = ORDMIN
      NSTEPS = 0
C--------- STARTING VALUES FOR NEWTON ITERATION
      DBLK = DBL(ORD)
      H    = MIN( H, ABS(TEND-T0)/DBLK )
      DBLKOLD = DBLK
      CALL FCN(R,T0,Y0,FP(1,1), RPAR,IPAR)
      NFCN = NFCN + 1
C -------- NUMBER OF FAILURES IN THE SPLITTING-NEWTON SCHEME
      FAILNI = 0
C -------- NUMBER OF FAILURES DUE TO THE ERROR TEST
      FAILEI = 0
      NSING  = 0
      ORDOLD = 2*ORD
      HOLD   = 2*H
      CALJAC = .TRUE.
      EXTRAP = .FALSE.
      ITMAX  = ITINT(ORD)
C--------- MAIN LOOP (ADVANCING IN TIME)
100   CONTINUE

C--------- (EVENTUALLY) COMPUTE THE JACOBIAN MATRIX NUMERICALLY
       NEWJAC = .FALSE.
       IF (CALJAC) THEN
          IF (IJAC.EQ.0) THEN
            DO I=1,R
               YSAFE=Y0(I)
               DELT=DSQRT(UROUND*DMAX1(1.D-5,DABS(YSAFE)))
               Y0(I)=YSAFE+DELT
               CALL FCN(R,T0,Y0,F,RPAR,IPAR)
               IF (JBAND) THEN
                 DO J=MAX(1,I-MUJAC),MIN(R,I+MLJAC)
                   JF0(J-I+MUJAC+1,I) = (F(J,1)-FP(J,1))/DELT
                 END DO
               ELSE
                 DO J=1,R
                   JF0(J,I)=(F(J,1)-FP(J,1))/DELT
                 END DO
               END IF
               Y0(I)=YSAFE
            END DO
         ELSE
C -------- COMPUTE JACOBIAN MATRIX ANALYTICALLY
c           CALL JAC(R,T0,Y0,JF0,LDJAC,RPAR,IPAR)
           CALL constV_jac_sparse(R,real(T0,dp),real(Y0, dp))
        END IF
         NJAC = NJAC + 1
         NEWJAC = .TRUE.
       END IF

C--------- FACTORIZE THE ITERATION MATRIX
      IF ((ORDOLD.NE.ORD).OR.(HOLD.NE.H).OR.(NEWJAC) ) THEN
        HOLD = H
        ORDOLD = ORD
        IER = 1

        DO WHILE ( IER .NE. 0)
          CALL  DECLU(R,JF0,H,LDJAC,LU,LDLU,IPIV,ORD,IER,IJOB)
          IF (IER.NE.0) THEN
            NSING = NSING + 1
            IF (NSING.GT.5) THEN
              WRITE(6,*) 'MATRIX IS REPEATEDLY SINGULAR, IER= ',IER
              WRITE(6,900) T0
              IDID=-4
              GOTO 800
            ELSE
              H = H/2D0
            END IF
          END IF
          NDEC = NDEC + 1
        END DO
      END IF

C--------- DEFINE TP AND YP
      IF (EXTRAP) THEN
         T1(1) = T0+H
         DO I=2,DBLK+1
           T1(I) = T1(I-1)+H
         END DO
         CALL INTERP(R,TP,YP,T1,F1,NT1,DBLKOLD,DBLK,T0,Y0,ORD)
      ELSE
         TP(1) = T0
         DO J=1,R
           YP(J,1) = Y0(J)
         END DO
         DO I=2,DBLK+1
           DO J=1,R
             YP(J,I) = Y0(J)
           END DO
           TP(I) = TP(I-1)+H
         END DO
      END IF

C--------- DEFINE SCAL AND FACN

      THN = 1d0
      J    = 0
      IF (ITOL.EQ.0) THEN
          DO I=1,R
            RRN = DABS(Y0(I))
            SCAL(I)=ATOL(1)+RTOL(1)*RRN
            IF (RRN .LT. THN) THEN
               J = I
               THN = RRN
            ENDIF
          END DO
      ELSE
          DO I=1,R
            RRN = DABS( Y0(I) )
            SCAL(I)=ATOL(I)+RTOL(I)*RRN
            IF (RRN .LT. THN) THEN
               J = I
               THN = RRN
            ENDIF
          END DO
      END IF



       FACN = FACNORD(ORD)
       IF (THN .LT. 1d-1) THEN
          IF (DABS(FP(J,1)) .LT. 1d-5) THEN
              FACN = MIN(FACNEWT, FACNORD(ORD) )
          END IF
       END IF


C---------- COMPUTE THE NUMERICAL SOLUTION AT TIMES T1(1)...T1(DBLK)
C---------- DEFINE VARIABLES NEEDED IN THE ITERATION
        ERRNEWT  = FACN+1d0
        ERRNEWT0 = FACN+1d0
C----------
        TETAK    = 1.0D0
        THETA    = 1.0D0
        TETAKOLD = 1.0D0
        ITMAX    = ITINT(ORD)
        IT  = 0
        DO J = 2, DBLK+1
           CALL FCN(R,TP(J),YP(1,J),FP(1,J), RPAR,IPAR)
        END DO
        NFCN = NFCN + DBLK
C
C--------- SPLITTING NEWTON LOOP
C
 300    CONTINUE
          ERRNEWT1 = ERRNEWT0
          ERRNEWT0 = ERRNEWT
          ERRNEWT  = 0D0
C--------- COMPUTE ONE ITERATION FOR THE SELECTED ORDER
          GOTO (101,201,301,401) ORD
101     CALL      TERMNOT3(R,FCN,H,IT,DN, F,FP,YP,TP,NFCN,
     &  ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV, SCAL,IJOB,TER,
     &  RPAR,IPAR)
        GOTO 501
201     CALL      TERMNOT5(R,FCN,H,IT,DN, F,FP,YP,TP,NFCN,
     &  ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV, SCAL,IJOB,TER,
     &  RPAR,IPAR)
        GOTO 501
301     CALL      TERMNOT7(R,FCN,H,IT,DN, F,FP,YP,TP,NFCN,
     &  ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV, SCAL,IJOB,TER,
     &  RPAR,IPAR)
        GOTO 501
401     CALL      TERMNOT9(R,FCN,H,IT,DN, F,FP,YP,TP,NFCN,
     &  ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV, SCAL,IJOB,TER,
     &  RPAR,IPAR)
501     CONTINUE
         IF (TER)  THEN
           ERRNEWT = FACN + 1
           GOTO 999
         END IF
C--------- COMPUTE TETAK, ETAK

           TETAKOLD=TETAK
           TETAK   = ERRNEWT/SQRT(ERRNEWT0*ERRNEWT1)
           IF (IT.LE.2) THEN
              THETA=THET/2d0
           ELSE IF (IT .GT. 2) THEN
              THETA = SQRT(TETAK*TETAKOLD)
           END IF

           IT = IT+1

           JVAI = (IT .LE. ITMAX).AND.(ERRNEWT.GT.FACN).AND.
     &((THETA.LT.TETAK0).OR.(IT.LE.2)).AND. (ERRNEWT.GT.0d0)
           IF (JVAI) GO TO 300
C
C--------- END OF NEWTON LOOP
C
999   CONTINUE
      IF ((ERRNEWT.GT.FACN).OR.(.not.(ERRNEWT.GT.0d0))) THEN
C--------- THE ITERATION DOES NOT CONVERGE
         FAILNI = FAILNI + 1
         NNEWT(ORD) = NNEWT(ORD)+1
C--------- CHOICE OF THE NEW STEPSIZE
         H=H/2d0
         DBLKOLD = DBLK
         EXTRAP = .FALSE.
         IF (FAILNI .EQ. 1) THEN
            CALJAC = .NOT. NEWJAC
         ELSE
            CALJAC = .FALSE.
         END IF
C--------- RETURN TO THE MAIN LOOP
       ELSE
C--------- THE ITERATION CONVERGES
C--------- ERROR ESTIMATION
        CALL  ESTERR(ERRV, ERRSAME, ERRUP, ERRDOWN, FP,
     &     R, H, ORD, DBLK, LU, LDLU,
     &     IPIV, F, F1, SCAL, ORDMAX,ORDMIN,IJOB)
         IF ( ERRSAME .GT. 1D0 ) THEN
           FAILEI = FAILEI + 1
           NERR(ORD) = NERR(ORD) + 1
           CALJAC = (THETA .GT. THET)
C--------- NEW STEPSIZE SELECTION
           ORD2 = 2*ORD
           ESP = 1D0/(ORD2+1D0)
           RRN=DMAX1(FACL,DMIN1(FACR,(SF*ERRSAME)**ESP))
           H = H/RRN
           DBLKOLD = DBLK
           DO I=1,DBLKOLD+1
            DO J=1,R
              F1(J,I) = YP(J,I)
            END DO
           END DO
           call DIFFDIV(TP,F1,R,DBLK,NT1)
           EXTRAP = .TRUE.

C--------- RETURN TO THE MAIN LOOP
         ELSE
C--------- THE STEPSIZE IS ACCEPTED
           NSTEP(ORD) = NSTEP(ORD)+1
c          write(55,*) H
           T0 = TP(DBLK+1)
           DO I=1, R
             Y0(I) = YP(I,DBLK+1)
             FP(I,1) = FP(I,DBLK+1)
           END DO
C--------- NEW STEPSIZE SELECTION
            ORD2 = 2*ORD
            ESP = 1D0/(ORD2+1d0)
            RRN=DMAX1(FACL,DMIN1(FACR,(SFSAME*ERRSAME)**ESP))
            THN=DBL(ORD)/(CS(ORD)*RRN)
            ORDN = ORD
            IF  (ORD.LT.ORDMAX) THEN
              ESP = 1D0/(ORD2+3D0)
              RR=DMAX1(FACL,DMIN1(FACR,(SFUP*ERRUP)**ESP))
              TH=DBL(ORD+1)/(CS(ORD+1)*RR )
              IF (TH .GT. THN ) THEN
                ORDN = ORD + 1
                RRN  = RR
                THN  = TH
              END IF
            END IF

            IF ( ORD.GT.ORDMIN)  THEN
             ESP = 1D0/(ORD2-1d0)
             RR=DMAX1(FACL,DMIN1(FACR,(SFDOWN*ERRDOWN)**ESP))
             TH=DBL(ORD-1)/(CS(ORD-1)*RR )
             IF (TH .GT. THN ) THEN
                ORDN = ORD - 1
                RRN  = RR
             END IF
            END IF
           HOLD = H
c
c
           IF (ORDN.GT.ORD) THEN
               H = MIN(H/RRN, HOLD)
           ELSE
               H = H/RRN
           END IF
           ORDOLD = ORD
           ORD = ORDN
           DBLKOLD = DBLK
           DBLK = DBL(ORD)

           CALJAC = (THETA .GT. THET)
           EXTRAP = .TRUE.

           IF ((FAILNI.NE.0).OR.(FAILEI.NE.0)) THEN
             H = DMIN1( H, HOLD)
           END IF
           IF  (.NOT. CALJAC) THEN
             IF ((H/HOLD.LE.1.1D0 ).AND.(H/HOLD.GE.0.9D0)) THEN
                H = HOLD
             END IF
           END IF
           H = DMIN1( H, DMIN1(HMAX, (TEND-T0)/DBLK) )


           DO I=1,DBLKOLD+1
            DO J=1,R
              F1(J,I) = YP(J,I)
            END DO
           END DO
           CALL DIFFDIV(TP,F1,R,DBLKOLD,NT1)
           EXTRAP = .TRUE.
c
        IF (IOUT.NE.0) THEN
C--------- CALL SOLOUT
          CALL SOLOUT(R,TP,YP,F1,NT1,DBLKOLD,ORDOLD,RPAR,IPAR,IRTRN)
          IF (IRTRN.LT.0) GOTO 800
        END IF
          IF (NSTEPS .EQ. 0) THEN
            FAILNI = 0
            FAILEI = 0
          ELSE
            FAILNI = MAX(FAILNI-1,0)
            FAILEI = MAX(FAILEI-1,0)
          END IF
          NSING  = 0
        END IF
c--------- END IF ERRSAME > 1
      END IF
c--------- END IF ERRNEWT > 1
           NSTEPS = NSTEPS + 1
           IF (0.1d0*DABS(T0-TEND)/DBLK .GT. dabs(T0)*UROUND ) THEN
            IF (0.1d0*DABS(H).LE.DABS(T0)*UROUND) THEN
              WRITE(6,*) ' STEPSIZE TOO SMALL, H=',H
              WRITE(6,900) T0
              IDID=-3
              GOTO 800
            END IF
            IF (NSTEPS.GT.NMAX) THEN
              WRITE(6,*) ' MORE THAN NMAX =',NMAX,'STEPS ARE NEEDED'
              WRITE(6,900) T0
              IDID=-2
              GOTO 800
            END IF

           GOTO 100
C
C---------- END WHILE T0 < T
           ELSE
             H    = HOLD
             IDID = 1
           END IF
 900  FORMAT(' EXIT OF GAM AT T=',E18.4)
 800  RETURN
      END
