SUBROUTINE RADAU5(N,FCN,X,Y,XEND,H,&
&RTOL,ATOL,ITOL,&
&JAC ,IJAC,MLJAC,MUJAC,&
&MAS ,IMAS,MLMAS,MUMAS,&
&SOLOUT,IOUT,&
&WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
! ----------------------------------------------------------
!     NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC)
!     SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS
!                     M*Y'=F(X,Y).
!     THE SYSTEM CAN BE (LINEARLY) IMPLICIT (MASS-MATRIX M .NE. I)
!     OR EXPLICIT (M=I).
!     THE METHOD USED IS AN IMPLICIT RUNGE-KUTTA METHOD (RADAU IIA)
!     OF ORDER 5 WITH STEP SIZE CONTROL AND CONTINUOUS OUTPUT.
!     CF. SECTION IV.8
!
!     AUTHORS: E. HAIRER AND G. WANNER
!              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
!              CH-1211 GENEVE 24, SWITZERLAND
!              E-MAIL:  Ernst.Hairer@math.unige.ch
!                       Gerhard.Wanner@math.unige.ch
!
!     THIS CODE IS PART OF THE BOOK:
!         E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!         EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS 14,
!         SPRINGER-VERLAG 1991, SECOND EDITION 1996.
!
!     VERSION OF JULY 9, 1996
!     (latest small correction: January 18, 2002)
!
!     INPUT PARAMETERS
!     ----------------
!     N           DIMENSION OF THE SYSTEM
!
!     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE
!                 VALUE OF F(X,Y):
!                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR)
!                    DOUBLE PRECISION X,Y(N),F(N)
!                    F(1)=...   ETC.
!                 RPAR, IPAR (SEE BELOW)
!
!     X           INITIAL X-VALUE
!
!     Y(N)        INITIAL VALUES FOR Y
!
!     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)
!
!     H           INITIAL STEP SIZE GUESS;
!                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT,
!                 H=1.D0/(NORM OF F'), USUALLY 1.D-3 OR 1.D-5, IS GOOD.
!                 THIS CHOICE IS NOT VERY IMPORTANT, THE STEP SIZE IS
!                 QUICKLY ADAPTED. (IF H=0.D0, THE CODE PUTS H=1.D-6).
!
!     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
!                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N.
!
!     ITOL        SWITCH FOR RTOL AND ATOL:
!                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
!                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
!                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL
!                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
!                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
!                     RTOL(I)*ABS(Y(I))+ATOL(I).
!
!     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
!                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO Y
!                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1; SUPPLY
!                 A DUMMY SUBROUTINE IN THE CASE IJAC=0).
!                 FOR IJAC=1, THIS SUBROUTINE MUST HAVE THE FORM
!                    SUBROUTINE JAC(N,X,Y,DFY,LDFY,RPAR,IPAR)
!                    DOUBLE PRECISION X,Y(N),DFY(LDFY,N)
!                    DFY(1,1)= ...
!                 LDFY, THE COLUMN-LENGTH OF THE ARRAY, IS
!                 FURNISHED BY THE CALLING PROGRAM.
!                 IF (MLJAC.EQ.N) THE JACOBIAN IS SUPPOSED TO
!                    BE FULL AND THE PARTIAL DERIVATIVES ARE
!                    STORED IN DFY AS
!                       DFY(I,J) = PARTIAL F(I) / PARTIAL Y(J)
!                 ELSE, THE JACOBIAN IS TAKEN AS BANDED AND
!                    THE PARTIAL DERIVATIVES ARE STORED
!                    DIAGONAL-WISE AS
!                       DFY(I-J+MUJAC+1,J) = PARTIAL F(I) / PARTIAL Y(J).
!
!     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN:
!                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE
!                       DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED.
!                    IJAC=1: JACOBIAN IS SUPPLIED BY SUBROUTINE JAC.
!
!     MLJAC       SWITCH FOR THE BANDED STRUCTURE OF THE JACOBIAN:
!                    MLJAC=N: JACOBIAN IS A FULL MATRIX. THE LINEAR
!                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
!                    0<=MLJAC<N: MLJAC IS THE LOWER BANDWITH OF JACOBIAN
!                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
!                       THE MAIN DIAGONAL).
!
!     MUJAC       UPPER BANDWITH OF JACOBIAN  MATRIX (>= NUMBER OF NON-
!                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
!                 NEED NOT BE DEFINED IF MLJAC=N.
!
!     ----   MAS,IMAS,MLMAS, AND MUMAS HAVE ANALOG MEANINGS      -----
!     ----   FOR THE "MASS MATRIX" (THE MATRIX "M" OF SECTION IV.8): -
!
!     MAS         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE MASS-
!                 MATRIX M.
!                 IF IMAS=0, THIS MATRIX IS ASSUMED TO BE THE IDENTITY
!                 MATRIX AND NEEDS NOT TO BE DEFINED;
!                 SUPPLY A DUMMY SUBROUTINE IN THIS CASE.
!                 IF IMAS=1, THE SUBROUTINE MAS IS OF THE FORM
!                    SUBROUTINE MAS(N,AM,LMAS,RPAR,IPAR)
!                    DOUBLE PRECISION AM(LMAS,N)
!                    AM(1,1)= ....
!                    IF (MLMAS.EQ.N) THE MASS-MATRIX IS STORED
!                    AS FULL MATRIX LIKE
!                         AM(I,J) = M(I,J)
!                    ELSE, THE MATRIX IS TAKEN AS BANDED AND STORED
!                    DIAGONAL-WISE AS
!                         AM(I-J+MUMAS+1,J) = M(I,J).
!
!     IMAS       GIVES INFORMATION ON THE MASS-MATRIX:
!                    IMAS=0: M IS SUPPOSED TO BE THE IDENTITY
!                       MATRIX, MAS IS NEVER CALLED.
!                    IMAS=1: MASS-MATRIX  IS SUPPLIED.
!
!     MLMAS       SWITCH FOR THE BANDED STRUCTURE OF THE MASS-MATRIX:
!                    MLMAS=N: THE FULL MATRIX CASE. THE LINEAR
!                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
!                    0<=MLMAS<N: MLMAS IS THE LOWER BANDWITH OF THE
!                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
!                       THE MAIN DIAGONAL).
!                 MLMAS IS SUPPOSED TO BE .LE. MLJAC.
!
!     MUMAS       UPPER BANDWITH OF MASS-MATRIX (>= NUMBER OF NON-
!                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
!                 NEED NOT BE DEFINED IF MLMAS=N.
!                 MUMAS IS SUPPOSED TO BE .LE. MUJAC.
!
!     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
!                 NUMERICAL SOLUTION DURING INTEGRATION.
!                 IF IOUT=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
!                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0.
!                 IT MUST HAVE THE FORM
!                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,
!                                       RPAR,IPAR,IRTRN)
!                    DOUBLE PRECISION X,Y(N),CONT(LRC)
!                    ....
!                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH
!                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS
!                    THE FIRST GRID-POINT).
!                 "XOLD" IS THE PRECEEDING GRID-POINT.
!                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
!                    IS SET <0, RADAU5 RETURNS TO THE CALLING PROGRAM.
!
!          -----  CONTINUOUS OUTPUT: -----
!                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
!                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH
!                 THE FUNCTION
!                        >>>   CONTR5(I,S,CONT,LRC)   <<<
!                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH
!                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE
!                 S SHOULD LIE IN THE INTERVAL [XOLD,X].
!                 DO NOT CHANGE THE ENTRIES OF CONT(LRC), IF THE
!                 DENSE OUTPUT FUNCTION IS USED.
!
!     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT:
!                    IOUT=0: SUBROUTINE IS NEVER CALLED
!                    IOUT=1: SUBROUTINE IS AVAILABLE FOR OUTPUT.
!
!     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
!                 WORK(1), WORK(2),.., WORK(20) SERVE AS PARAMETERS
!                 FOR THE CODE. FOR STANDARD USE OF THE CODE
!                 WORK(1),..,WORK(20) MUST BE SET TO ZERO BEFORE
!                 CALLING. SEE BELOW FOR A MORE SOPHISTICATED USE.
!                 WORK(21),..,WORK(LWORK) SERVE AS WORKING SPACE
!                 FOR ALL VECTORS AND MATRICES.
!                 "LWORK" MUST BE AT LEAST
!                             N*(LJAC+LMAS+3*LE+12)+20
!                 WHERE
!                    LJAC=N              IF MLJAC=N (FULL JACOBIAN)
!                    LJAC=MLJAC+MUJAC+1  IF MLJAC<N (BANDED JAC.)
!                 AND
!                    LMAS=0              IF IMAS=0
!                    LMAS=N              IF IMAS=1 AND MLMAS=N (FULL)
!                    LMAS=MLMAS+MUMAS+1  IF MLMAS<N (BANDED MASS-M.)
!                 AND
!                    LE=N               IF MLJAC=N (FULL JACOBIAN)
!                    LE=2*MLJAC+MUJAC+1 IF MLJAC<N (BANDED JAC.)
!
!                 IN THE USUAL CASE WHERE THE JACOBIAN IS FULL AND THE
!                 MASS-MATRIX IS THE INDENTITY (IMAS=0), THE MINIMUM
!                 STORAGE REQUIREMENT IS
!                             LWORK = 4*N*N+12*N+20.
!                 IF IWORK(9)=M1>0 THEN "LWORK" MUST BE AT LEAST
!                          N*(LJAC+12)+(N-M1)*(LMAS+3*LE)+20
!                 WHERE IN THE DEFINITIONS OF LJAC, LMAS AND LE THE
!                 NUMBER N CAN BE REPLACED BY N-M1.
!
!     LWORK       DECLARED LENGTH OF ARRAY "WORK".
!
!     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK".
!                 IWORK(1),IWORK(2),...,IWORK(20) SERVE AS PARAMETERS
!                 FOR THE CODE. FOR STANDARD USE, SET IWORK(1),..,
!                 IWORK(20) TO ZERO BEFORE CALLING.
!                 IWORK(21),...,IWORK(LIWORK) SERVE AS WORKING AREA.
!                 "LIWORK" MUST BE AT LEAST 3*N+20.
!
!     LIWORK      DECLARED LENGTH OF ARRAY "IWORK".
!
!     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH
!                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING
!                 PROGRAM AND THE FCN, JAC, MAS, SOLOUT SUBROUTINES.
!
! ----------------------------------------------------------------------
!
!     SOPHISTICATED SETTING OF PARAMETERS
!     -----------------------------------
!              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK
!              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),...
!              AS WELL AS IWORK(1),... DIFFERENT FROM ZERO.
!              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES:
!
!    IWORK(1)  IF IWORK(1).NE.0, THE CODE TRANSFORMS THE JACOBIAN
!              MATRIX TO HESSENBERG FORM. THIS IS PARTICULARLY
!              ADVANTAGEOUS FOR LARGE SYSTEMS WITH FULL JACOBIAN.
!              IT DOES NOT WORK FOR BANDED JACOBIAN (MLJAC<N)
!              AND NOT FOR IMPLICIT SYSTEMS (IMAS=1).
!
!    IWORK(2)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
!              THE DEFAULT VALUE (FOR IWORK(2)=0) IS 100000.
!
!    IWORK(3)  THE MAXIMUM NUMBER OF NEWTON ITERATIONS FOR THE
!              SOLUTION OF THE IMPLICIT SYSTEM IN EACH STEP.
!              THE DEFAULT VALUE (FOR IWORK(3)=0) IS 7.
!
!    IWORK(4)  IF IWORK(4).EQ.0 THE EXTRAPOLATED COLLOCATION SOLUTION
!              IS TAKEN AS STARTING VALUE FOR NEWTON'S METHOD.
!              IF IWORK(4).NE.0 ZERO STARTING VALUES ARE USED.
!              THE LATTER IS RECOMMENDED IF NEWTON'S METHOD HAS
!              DIFFICULTIES WITH CONVERGENCE (THIS IS THE CASE WHEN
!              NSTEP IS LARGER THAN NACCPT + NREJCT; SEE OUTPUT PARAM.).
!              DEFAULT IS IWORK(4)=0.
!
!       THE FOLLOWING 3 PARAMETERS ARE IMPORTANT FOR
!       DIFFERENTIAL-ALGEBRAIC SYSTEMS OF INDEX > 1.
!       THE FUNCTION-SUBROUTINE SHOULD BE WRITTEN SUCH THAT
!       THE INDEX 1,2,3 VARIABLES APPEAR IN THIS ORDER.
!       IN ESTIMATING THE ERROR THE INDEX 2 VARIABLES ARE
!       MULTIPLIED BY H, THE INDEX 3 VARIABLES BY H**2.
!
!    IWORK(5)  DIMENSION OF THE INDEX 1 VARIABLES (MUST BE > 0). FOR
!              ODE'S THIS EQUALS THE DIMENSION OF THE SYSTEM.
!              DEFAULT IWORK(5)=N.
!
!    IWORK(6)  DIMENSION OF THE INDEX 2 VARIABLES. DEFAULT IWORK(6)=0.
!
!    IWORK(7)  DIMENSION OF THE INDEX 3 VARIABLES. DEFAULT IWORK(7)=0.
!
!    IWORK(8)  SWITCH FOR STEP SIZE STRATEGY
!              IF IWORK(8).EQ.1  MOD. PREDICTIVE CONTROLLER (GUSTAFSSON)
!              IF IWORK(8).EQ.2  CLASSICAL STEP SIZE CONTROL
!              THE DEFAULT VALUE (FOR IWORK(8)=0) IS IWORK(8)=1.
!              THE CHOICE IWORK(8).EQ.1 SEEMS TO PRODUCE SAFER RESULTS;
!              FOR SIMPLE PROBLEMS, THE CHOICE IWORK(8).EQ.2 PRODUCES
!              OFTEN SLIGHTLY FASTER RUNS
!
!       IF THE DIFFERENTIAL SYSTEM HAS THE SPECIAL STRUCTURE THAT
!            Y(I)' = Y(I+M2)   FOR  I=1,...,M1,
!       WITH M1 A MULTIPLE OF M2, A SUBSTANTIAL GAIN IN COMPUTERTIME
!       CAN BE ACHIEVED BY SETTING THE PARAMETERS IWORK(9) AND IWORK(10).
!       E.G., FOR SECOND ORDER SYSTEMS P'=V, V'=G(P,V), WHERE P AND V ARE
!       VECTORS OF DIMENSION N/2, ONE HAS TO PUT M1=M2=N/2.
!       FOR M1>0 SOME OF THE INPUT PARAMETERS HAVE DIFFERENT MEANINGS:
!       - JAC: ONLY THE ELEMENTS OF THE NON-TRIVIAL PART OF THE
!              JACOBIAN HAVE TO BE STORED
!              IF (MLJAC.EQ.N-M1) THE JACOBIAN IS SUPPOSED TO BE FULL
!                 DFY(I,J) = PARTIAL F(I+M1) / PARTIAL Y(J)
!                FOR I=1,N-M1 AND J=1,N.
!              ELSE, THE JACOBIAN IS BANDED ( M1 = M2 * MM )
!                 DFY(I-J+MUJAC+1,J+K*M2) = PARTIAL F(I+M1) / PARTIAL Y(J+K*M2)
!                FOR I=1,MLJAC+MUJAC+1 AND J=1,M2 AND K=0,MM.
!       - MLJAC: MLJAC=N-M1: IF THE NON-TRIVIAL PART OF THE JACOBIAN IS FULL
!                0<=MLJAC<N-M1: IF THE (MM+1) SUBMATRICES (FOR K=0,MM)
!                     PARTIAL F(I+M1) / PARTIAL Y(J+K*M2),  I,J=1,M2
!                    ARE BANDED, MLJAC IS THE MAXIMAL LOWER BANDWIDTH
!                    OF THESE MM+1 SUBMATRICES
!       - MUJAC: MAXIMAL UPPER BANDWIDTH OF THESE MM+1 SUBMATRICES
!                NEED NOT BE DEFINED IF MLJAC=N-M1
!       - MAS: IF IMAS=0 THIS MATRIX IS ASSUMED TO BE THE IDENTITY AND
!              NEED NOT BE DEFINED. SUPPLY A DUMMY SUBROUTINE IN THIS CASE.
!              IT IS ASSUMED THAT ONLY THE ELEMENTS OF RIGHT LOWER BLOCK OF
!              DIMENSION N-M1 DIFFER FROM THAT OF THE IDENTITY MATRIX.
!              IF (MLMAS.EQ.N-M1) THIS SUBMATRIX IS SUPPOSED TO BE FULL
!                 AM(I,J) = M(I+M1,J+M1)     FOR I=1,N-M1 AND J=1,N-M1.
!              ELSE, THE MASS MATRIX IS BANDED
!                 AM(I-J+MUMAS+1,J) = M(I+M1,J+M1)
!       - MLMAS: MLMAS=N-M1: IF THE NON-TRIVIAL PART OF M IS FULL
!                0<=MLMAS<N-M1: LOWER BANDWIDTH OF THE MASS MATRIX
!       - MUMAS: UPPER BANDWIDTH OF THE MASS MATRIX
!                NEED NOT BE DEFINED IF MLMAS=N-M1
!
!    IWORK(9)  THE VALUE OF M1.  DEFAULT M1=0.
!
!    IWORK(10) THE VALUE OF M2.  DEFAULT M2=M1.
!
! ----------
!
!    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.
!
!    WORK(2)   THE SAFETY FACTOR IN STEP SIZE PREDICTION,
!              DEFAULT 0.9D0.
!
!    WORK(3)   DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED;
!              INCREASE WORK(3), TO 0.1 SAY, WHEN JACOBIAN EVALUATIONS
!              ARE COSTLY. FOR SMALL SYSTEMS WORK(3) SHOULD BE SMALLER
!              (0.001D0, SAY). NEGATIV WORK(3) FORCES THE CODE TO
!              COMPUTE THE JACOBIAN AFTER EVERY ACCEPTED STEP.
!              DEFAULT 0.001D0.
!
!    WORK(4)   STOPPING CRITERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1.
!              SMALLER VALUES OF WORK(4) MAKE THE CODE SLOWER, BUT SAFER.
!              DEFAULT MIN(0.03D0,RTOL(1)**0.5D0)
!
!    WORK(5) AND WORK(6) : IF WORK(5) < HNEW/HOLD < WORK(6), THEN THE
!              STEP SIZE IS NOT CHANGED. THIS SAVES, TOGETHER WITH A
!              LARGE WORK(3), LU-DECOMPOSITIONS AND COMPUTING TIME FOR
!              LARGE SYSTEMS. FOR SMALL SYSTEMS ONE MAY HAVE
!              WORK(5)=1.D0, WORK(6)=1.2D0, FOR LARGE FULL SYSTEMS
!              WORK(5)=0.99D0, WORK(6)=2.D0 MIGHT BE GOOD.
!              DEFAULTS WORK(5)=1.D0, WORK(6)=1.2D0 .
!
!    WORK(7)   MAXIMAL STEP SIZE, DEFAULT XEND-X.
!
!    WORK(8), WORK(9)   PARAMETERS FOR STEP SIZE SELECTION
!              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION
!                 WORK(8) <= HNEW/HOLD <= WORK(9)
!              DEFAULT VALUES: WORK(8)=0.2D0, WORK(9)=8.D0
!
!-----------------------------------------------------------------------
!
!     OUTPUT PARAMETERS
!     -----------------
!     X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED
!                 (AFTER SUCCESSFUL RETURN X=XEND).
!
!     Y(N)        NUMERICAL SOLUTION AT X
!
!     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
!
!     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
!                   IDID= 1  COMPUTATION SUCCESSFUL,
!                   IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT)
!                   IDID=-1  INPUT IS NOT CONSISTENT,
!                   IDID=-2  LARGER NMAX IS NEEDED,
!                   IDID=-3  STEP SIZE BECOMES TOO SMALL,
!                   IDID=-4  MATRIX IS REPEATEDLY SINGULAR.
!
!   IWORK(14)  NFCN    NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL
!                      EVALUATION OF THE JACOBIAN ARE NOT COUNTED)
!   IWORK(15)  NJAC    NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY
!                      OR NUMERICALLY)
!   IWORK(16)  NSTEP   NUMBER OF COMPUTED STEPS
!   IWORK(17)  NACCPT  NUMBER OF ACCEPTED STEPS
!   IWORK(18)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST),
!                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED)
!   IWORK(19)  NDEC    NUMBER OF LU-DECOMPOSITIONS OF BOTH MATRICES
!   IWORK(20)  NSOL    NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS, OF BOTH
!                      SYSTEMS; THE NSTEP FORWARD-BACKWARD SUBSTITUTIONS,
!                      NEEDED FOR STEP SIZE SELECTION, ARE NOT COUNTED
!-----------------------------------------------------------------------
! *** *** *** *** *** *** *** *** *** *** *** *** ***
!          DECLARATIONS
! *** *** *** *** *** *** *** *** *** *** *** *** ***
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)


   DIMENSION Y(N),ATOL(*),RTOL(*),WORK(LWORK),IWORK(LIWORK)
   DIMENSION RPAR(*),IPAR(*)
   LOGICAL IMPLCT,JBAND,ARRET,STARTN,PRED
   EXTERNAL FCN,JAC,MAS,SOLOUT, constV_jac_RADAUS
! *** *** *** *** *** *** ***
!        SETTING THE PARAMETERS
! *** *** *** *** *** *** ***
   NFCN=0
   NJAC=0
   NSTEP=0
   NACCPT=0
   NREJCT=0
   NDEC=0
   NSOL=0
   ARRET=.FALSE.

! -------- UROUND   SMALLEST NUMBER SATISFYING 1.0D0+UROUND>1.0D0
   IF (WORK(1).EQ.0.0D0) THEN
      UROUND=1.0D-19
   ELSE
      UROUND=WORK(1)
      IF (UROUND.LE.1.0D-19.OR.UROUND.GE.1.0D0) THEN
         WRITE(6,*)' COEFFICIENTS HAVE 20 DIGITS, UROUND=',WORK(1)
         ARRET=.TRUE.
      END IF
   END IF
! -------- CHECK AND CHANGE THE TOLERANCES
   EXPM=2.0D0/3.0D0
   IF (ITOL.EQ.0) THEN
      IF (ATOL(1).LE.0.D0.OR.RTOL(1).LE.10.D0*UROUND) THEN
         WRITE (6,*) ' TOLERANCES ARE TOO SMALL'
         ARRET=.TRUE.
      ELSE
         QUOT=ATOL(1)/RTOL(1)
         RTOL(1)=0.1D0*RTOL(1)**EXPM
         ATOL(1)=RTOL(1)*QUOT
      END IF
   ELSE
      DO I=1,N
         IF (ATOL(I).LE.0.D0.OR.RTOL(I).LE.10.D0*UROUND) THEN
            WRITE (6,*) ' TOLERANCES(',I,') ARE TOO SMALL',UROUND
            ARRET=.TRUE.
         ELSE
            QUOT=ATOL(I)/RTOL(I)
            RTOL(I)=0.1D0*RTOL(I)**EXPM
            ATOL(I)=RTOL(I)*QUOT
         END IF
      END DO
   END IF
! -------- NMAX , THE MAXIMAL NUMBER OF STEPS -----
   IF (IWORK(2).EQ.0) THEN
      NMAX=100000
   ELSE
      NMAX=IWORK(2)
      IF (NMAX.LE.0) THEN
         WRITE(6,*)' WRONG INPUT IWORK(2)=',IWORK(2)
         ARRET=.TRUE.
      END IF
   END IF
! -------- NIT    MAXIMAL NUMBER OF NEWTON ITERATIONS
   IF (IWORK(3).EQ.0) THEN
      NIT=7
   ELSE
      NIT=IWORK(3)
      IF (NIT.LE.0) THEN
         WRITE(6,*)' CURIOUS INPUT IWORK(3)=',IWORK(3)
         ARRET=.TRUE.
      END IF
   END IF
! -------- STARTN  SWITCH FOR STARTING VALUES OF NEWTON ITERATIONS
   IF(IWORK(4).EQ.0)THEN
      STARTN=.FALSE.
   ELSE
      STARTN=.TRUE.
   END IF
! -------- PARAMETER FOR DIFFERENTIAL-ALGEBRAIC COMPONENTS
   NIND1=IWORK(5)
   NIND2=IWORK(6)
   NIND3=IWORK(7)
   IF (NIND1.EQ.0) NIND1=N
   IF (NIND1+NIND2+NIND3.NE.N) THEN
      WRITE(6,*)' CURIOUS INPUT FOR IWORK(5,6,7)=',NIND1,NIND2,NIND3
      ARRET=.TRUE.
   END IF
! -------- PRED   STEP SIZE CONTROL
   IF(IWORK(8).LE.1)THEN
      PRED=.TRUE.
   ELSE
      PRED=.FALSE.
   END IF
! -------- PARAMETER FOR SECOND ORDER EQUATIONS
   M1=IWORK(9)
   M2=IWORK(10)
   NM1=N-M1
   IF (M1.EQ.0) M2=N
   IF (M2.EQ.0) M2=M1
   IF (M1.LT.0.OR.M2.LT.0.OR.M1+M2.GT.N) THEN
      WRITE(6,*)' CURIOUS INPUT FOR IWORK(9,10)=',M1,M2
      ARRET=.TRUE.
   END IF
! --------- SAFE     SAFETY FACTOR IN STEP SIZE PREDICTION
   IF (WORK(2).EQ.0.0D0) THEN
      SAFE=0.9D0
   ELSE
      SAFE=WORK(2)
      IF (SAFE.LE.0.001D0.OR.SAFE.GE.1.0D0) THEN
         WRITE(6,*)' CURIOUS INPUT FOR WORK(2)=',WORK(2)
         ARRET=.TRUE.
      END IF
   END IF
! ------ THET     DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED;
   IF (WORK(3).EQ.0.D0) THEN
      THET=0.001D0
   ELSE
      THET=WORK(3)
      IF (THET.GE.1.0D0) THEN
         WRITE(6,*)' CURIOUS INPUT FOR WORK(3)=',WORK(3)
         ARRET=.TRUE.
      END IF
   END IF
! --- FNEWT   STOPPING CRITERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1.
   TOLST=RTOL(1)
   IF (WORK(4).EQ.0.D0) THEN
      FNEWT=MAX(10*UROUND/TOLST,MIN(0.03D0,TOLST**0.5D0))
   ELSE
      FNEWT=WORK(4)
      IF (FNEWT.LE.UROUND/TOLST) THEN
         WRITE(6,*)' CURIOUS INPUT FOR WORK(4)=',WORK(4)
         ARRET=.TRUE.
      END IF
   END IF
! --- QUOT1 AND QUOT2: IF QUOT1 < HNEW/HOLD < QUOT2, STEP SIZE = CONST.
   IF (WORK(5).EQ.0.D0) THEN
      QUOT1=1.D0
   ELSE
      QUOT1=WORK(5)
   END IF
   IF (WORK(6).EQ.0.D0) THEN
      QUOT2=1.2D0
   ELSE
      QUOT2=WORK(6)
   END IF
   IF (QUOT1.GT.1.0D0.OR.QUOT2.LT.1.0D0) THEN
      WRITE(6,*)' CURIOUS INPUT FOR WORK(5,6)=',QUOT1,QUOT2
      ARRET=.TRUE.
   END IF
! -------- MAXIMAL STEP SIZE
   IF (WORK(7).EQ.0.D0) THEN
      HMAX=XEND-X
   ELSE
      HMAX=WORK(7)
   END IF
! -------  FACL,FACR     PARAMETERS FOR STEP SIZE SELECTION
   IF(WORK(8).EQ.0.D0)THEN
      FACL=5.D0
   ELSE
      FACL=1.D0/WORK(8)
   END IF
   IF(WORK(9).EQ.0.D0)THEN
      FACR=1.D0/8.0D0
   ELSE
      FACR=1.D0/WORK(9)
   END IF
   IF (FACL.LT.1.0D0.OR.FACR.GT.1.0D0) THEN
      WRITE(6,*)' CURIOUS INPUT WORK(8,9)=',WORK(8),WORK(9)
      ARRET=.TRUE.
   END IF
! *** *** *** *** *** *** *** *** *** *** *** *** ***
!         COMPUTATION OF ARRAY ENTRIES
! *** *** *** *** *** *** *** *** *** *** *** *** ***
! ---- IMPLICIT, BANDED OR NOT ?
   IMPLCT=IMAS.NE.0
   JBAND=MLJAC.LT.NM1
! -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS ---
! -- JACOBIAN  AND  MATRICES E1, E2
   IF (JBAND) THEN
      LDJAC=MLJAC+MUJAC+1
      LDE1=MLJAC+LDJAC
   ELSE
      MLJAC=NM1
      MUJAC=NM1
      LDJAC=NM1
      LDE1=NM1
   END IF
! -- MASS MATRIX
   IF (IMPLCT) THEN
      IF (MLMAS.NE.NM1) THEN
         LDMAS=MLMAS+MUMAS+1
         IF (JBAND) THEN
            IJOB=4
         ELSE
            IJOB=3
         END IF
      ELSE
         MUMAS=NM1
         LDMAS=NM1
         IJOB=5
      END IF
! ------ BANDWITH OF "MAS" NOT SMALLER THAN BANDWITH OF "JAC"
      IF (MLMAS.GT.MLJAC.OR.MUMAS.GT.MUJAC) THEN
         WRITE (6,*) 'BANDWITH OF "MAS" NOT SMALLER THAN BANDWITH OF&
         & "JAC"'
         ARRET=.TRUE.
      END IF
   ELSE
      LDMAS=0
      IF (JBAND) THEN
         IJOB=2
      ELSE
         IJOB=1
         IF (N.GT.2.AND.IWORK(1).NE.0) IJOB=7
      END IF
   END IF
   LDMAS2=MAX(1,LDMAS)
! ------ HESSENBERG OPTION ONLY FOR EXPLICIT EQU. WITH FULL JACOBIAN
   IF ((IMPLCT.OR.JBAND).AND.IJOB.EQ.7) THEN
      WRITE(6,*)' HESSENBERG OPTION ONLY FOR EXPLICIT EQUATIONS WITH&
      &FULL JACOBIAN'
      ARRET=.TRUE.
   END IF
! ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK -----
   IEZ1=21
   IEZ2=IEZ1+N
   IEZ3=IEZ2+N
   IEY0=IEZ3+N
   IESCAL=IEY0+N
   IEF1=IESCAL+N
   IEF2=IEF1+N
   IEF3=IEF2+N
   IECON=IEF3+N
   IEJAC=IECON+4*N
   IEMAS=IEJAC+N*LDJAC
   IEE1=IEMAS+NM1*LDMAS
   IEE2R=IEE1+NM1*LDE1
   IEE2I=IEE2R+NM1*LDE1
! ------ TOTAL STORAGE REQUIREMENT -----------
   ISTORE=IEE2I+NM1*LDE1-1
   IF(ISTORE.GT.LWORK)THEN
      WRITE(6,*)' INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=',ISTORE
      ARRET=.TRUE.
   END IF
! ------- ENTRY POINTS FOR INTEGER WORKSPACE -----
   IEIP1=21
   IEIP2=IEIP1+NM1
   IEIPH=IEIP2+NM1
! --------- TOTAL REQUIREMENT ---------------
   ISTORE=IEIPH+NM1-1
   IF (ISTORE.GT.LIWORK) THEN
      WRITE(6,*)' INSUFF. STORAGE FOR IWORK, MIN. LIWORK=',ISTORE
      ARRET=.TRUE.
   END IF
! ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1
   IF (ARRET) THEN
      IDID=-1
      RETURN
   END IF
! -------- CALL TO CORE INTEGRATOR ------------
   CALL RADCOR(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,&
   &JAC,IJAC,MLJAC,MUJAC,MAS,MLMAS,MUMAS,SOLOUT,IOUT,IDID,&
   &NMAX,UROUND,SAFE,THET,FNEWT,QUOT1,QUOT2,NIT,IJOB,STARTN,&
   &NIND1,NIND2,NIND3,PRED,FACL,FACR,M1,M2,NM1,&
   &IMPLCT,JBAND,LDJAC,LDE1,LDMAS2,WORK(IEZ1),WORK(IEZ2),&
   &WORK(IEZ3),WORK(IEY0),WORK(IESCAL),WORK(IEF1),WORK(IEF2),&
   &WORK(IEF3),WORK(IEJAC),WORK(IEE1),WORK(IEE2R),WORK(IEE2I),&
   &WORK(IEMAS),IWORK(IEIP1),IWORK(IEIP2),IWORK(IEIPH),&
   &WORK(IECON),NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL,RPAR,IPAR)
   IWORK(14)=NFCN
   IWORK(15)=NJAC
   IWORK(16)=NSTEP
   IWORK(17)=NACCPT
   IWORK(18)=NREJCT
   IWORK(19)=NDEC
   IWORK(20)=NSOL
! -------- RESTORE TOLERANCES
   EXPM=1.0D0/EXPM
   IF (ITOL.EQ.0) THEN
      QUOT=ATOL(1)/RTOL(1)
      RTOL(1)=(10.0D0*RTOL(1))**EXPM
      ATOL(1)=RTOL(1)*QUOT
   ELSE
      DO I=1,N
         QUOT=ATOL(I)/RTOL(I)
         RTOL(I)=(10.0D0*RTOL(I))**EXPM
         ATOL(I)=RTOL(I)*QUOT
      END DO
   END IF
! ----------- RETURN -----------
   RETURN
END
!
!     END OF SUBROUTINE RADAU5
!
! ***********************************************************
!
SUBROUTINE RADCOR(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,&
&JAC,IJAC,MLJAC,MUJAC,MAS,MLMAS,MUMAS,SOLOUT,IOUT,IDID,&
&NMAX,UROUND,SAFE,THET,FNEWT,QUOT1,QUOT2,NIT,IJOB,STARTN,&
&NIND1,NIND2,NIND3,PRED,FACL,FACR,M1,M2,NM1,&
&IMPLCT,BANDED,LDJAC,LDE1,LDMAS,Z1,Z2,Z3,&
&Y0,SCAL,F1,F2,F3,FJAC,E1,E2R,E2I,FMAS,IP1,IP2,IPHES,&
&CONT,NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL,RPAR,IPAR)

   use working_precision, only: dp
   use sparse_algebra, only:sparse,sparse_to_dense,dense_to_sparse,&
   &print_sparse_to_file,&
   &sparse_sum, sparse_block_diagonal,&
   &identity


   use sparse_definitions
   use ode_solver, only: iper1,iiper1,iper2,iiper2,iLU1,iLU2,rLU1,&
   &rLU2,liLU1,lrLU1,liLU2,lrLU2, R5_sys1,&
   &R5_sys2

   use speedchem_conV, only: constV_jac_RADAUS

! ----------------------------------------------------------
!     CORE INTEGRATOR FOR RADAU5
!     PARAMETERS SAME AS IN RADAU5 WITH WORKSPACE ADDED
! ----------------------------------------------------------
!         DECLARATIONS
! ----------------------------------------------------------
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)

   interface

      function perturbed(Y,RTOL,ATOL) result(perturbedY)
         implicit none
         double precision, dimension(:), intent(in) :: Y
         double precision, dimension(:), intent(in) :: RTOL, ATOL
         double precision, dimension(size(Y))       :: perturbedY
      end function perturbed

      subroutine jacf(neq,time,yin,jac)
         use sparse_algebra, only: sparse
         implicit none
         integer,          intent(in) :: neq
         double precision, intent(in) :: time
         double precision, dimension(neq), intent(in)  :: yin
         type(sparse),                     intent(out) :: jac
      end subroutine jacf

      function r5_newton_matrix(n,a,b,g,jac) result(mat)
         use sparse_algebra
         implicit none
         integer,          intent(in) :: n
         double precision, intent(in) :: a, b, g
         type(sparse), intent(in) :: jac
         type(sparse)             :: mat
      end function r5_newton_matrix

      subroutine numerical_factorization(spmat,ip,iip,iLU,rLU)
         use sparse_algebra, only: sparse
         implicit none
         type(sparse), intent(inout) :: spmat
         integer,          dimension(:), intent(inout) :: ip, iip
         integer,          dimension(:), intent(inout) :: iLU
         double precision, dimension(:), intent(inout) :: rLU
      end subroutine numerical_factorization


      function sparse_linear_system(A,b,ip,iip,iLU,rLU) result(x)
         use sparse_algebra, only: sparse
         implicit none
         type(sparse),                   intent(in) :: A
         integer, dimension(:),          intent(inout) :: ip, iip
         integer, dimension(:),          intent(inout) :: iLU
         double precision, dimension(:), intent(inout) :: rLU
         double precision, dimension(:), intent(in) :: b
         double precision, dimension(size(b))       :: x
      end function sparse_linear_system

      function yale_work_dimensions(n,n2) result(lwork)
         implicit none
         integer, intent(in) :: n
         integer, intent(in), optional :: n2
         integer             :: lwork
      end function yale_work_dimensions

      subroutine system_1_prepro(neq,jacmat,gam,ip1,iip1,iLU1,rLU1)
         use sparse_algebra
         implicit none
         integer,                          intent(in) :: neq
         type(sparse),                     intent(in) :: jacmat
         double precision,                 intent(in) :: gam
         integer, dimension(neq),          intent(out):: ip1, iip1
         integer, dimension(:),            intent(inout) :: iLU1
         double precision, dimension(:),   intent(inout) :: rLU1
      end subroutine system_1_prepro

      subroutine system_2_prepro(neq,jacmat,a,b,ip1,iip1,iLU1,rLU1)
         use sparse_algebra
         implicit none
         integer,                          intent(in) :: neq
         type(sparse),                     intent(in) :: jacmat
         double precision,                 intent(in) :: a,b
         integer, dimension(neq),          intent(out):: ip1, iip1
         integer, dimension(:),            intent(inout) :: iLU1
         double precision, dimension(:),   intent(inout) :: rLU1
      end subroutine system_2_prepro

      function system_2_matrix(alfa,beta,jac) result(s2mat)
         use sparse_algebra
         implicit none
         type(sparse),     intent(in) :: jac
         double precision, intent(in) :: alfa, beta
         type(sparse)                 :: s2mat
      end function system_2_matrix

   end interface

   DIMENSION Y(N),Z1(N),Z2(N),Z3(N),Y0(N),SCAL(N),F1(N),F2(N),F3(N)
   DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),CONT(4*N)
   DIMENSION E1(LDE1,NM1),E2R(LDE1,NM1),E2I(LDE1,NM1)
   DIMENSION ATOL(*),RTOL(*),RPAR(*),IPAR(*)
   INTEGER IP1(NM1),IP2(NM1),IPHES(NM1)
   COMMON /CONRA5/NN,NN2,NN3,NN4,XSOL,HSOL,C2M1,C1M1
   COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
!$ OMP THREADPRIVATE(/CONRA5/,/LINAL/)
   LOGICAL REJECT,FIRST,IMPLCT,BANDED,CALJAC,STARTN,CALHES
   LOGICAL INDEX1,INDEX2,INDEX3,LAST,PRED
   EXTERNAL FCN

   DIMENSION E1TMP(LDE1,NM1)

   type(sparse) :: sparseJ, msys1, msys2
   double precision, dimension(n)    :: pertY
   double precision, dimension(LDE1) :: solved_sparse
   double precision, dimension(:), allocatable :: RTO2, ATO2
   double precision, dimension(n)    :: unk1, res1
   double precision, dimension(2*n)  :: unk2, res2
   real (dp) :: Xdp, Ydp(N)

!     Two are the systems of equations to be solved
!c     NB variables are now taken from ode_solver module!
!      integer, dimension(n)             :: iper1, iiper1
!      integer, dimension(2*n)           :: iper2, iiper2
!      integer, dimension(:), allocatable :: iLU1, iLU2
!      integer                            :: liLU1, lrLU1, liLU2, lrLU2
!      double precision, dimension(:), allocatable :: rLU1, rLU2



! *** *** *** *** *** *** ***
!  INITIALISATIONS
! *** *** *** *** *** *** ***
   IER = 0
!      allocate(RTO2(N), ATO2(N))
!      RTO2(1:N) = RTOL(1:N)
!      ATO2(1:N) = ATOL(1:N)



! --------- DUPLIFY N FOR COMMON BLOCK CONT -----
   NN=N
   NN2=2*N
   NN3=3*N
   LRC=4*N
! -------- CHECK THE INDEX OF THE PROBLEM -----
   INDEX1=NIND1.NE.0
   INDEX2=NIND2.NE.0
   INDEX3=NIND3.NE.0
! ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ----------
   IF (IMPLCT) CALL MAS(NM1,FMAS,LDMAS,RPAR,IPAR)
! ---------- CONSTANTS ---------
   SQ6=DSQRT(6.D0)
   C1=(4.D0-SQ6)/10.D0
   C2=(4.D0+SQ6)/10.D0
   C1M1=C1-1.D0
   C2M1=C2-1.D0
   C1MC2=C1-C2
   DD1=-(13.D0+7.D0*SQ6)/3.D0
   DD2=(-13.D0+7.D0*SQ6)/3.D0
   DD3=-1.D0/3.D0
   U1=(6.D0+81.D0**(1.D0/3.D0)-9.D0**(1.D0/3.D0))/30.D0
   ALPH=(12.D0-81.D0**(1.D0/3.D0)+9.D0**(1.D0/3.D0))/60.D0
   BETA=(81.D0**(1.D0/3.D0)+9.D0**(1.D0/3.D0))*DSQRT(3.D0)/60.D0
   CNO=ALPH**2+BETA**2
   U1=1.0D0/U1
   ALPH=ALPH/CNO
   BETA=BETA/CNO
   T11=9.1232394870892942792D-02
   T12=-0.14125529502095420843D0
   T13=-3.0029194105147424492D-02
   T21=0.24171793270710701896D0
   T22=0.20412935229379993199D0
   T23=0.38294211275726193779D0
   T31=0.96604818261509293619D0
   TI11=4.3255798900631553510D0
   TI12=0.33919925181580986954D0
   TI13=0.54177053993587487119D0
   TI21=-4.1787185915519047273D0
   TI22=-0.32768282076106238708D0
   TI23=0.47662355450055045196D0
   TI31=-0.50287263494578687595D0
   TI32=2.5719269498556054292D0
   TI33=-0.59603920482822492497D0


   IF (M1.GT.0) IJOB=IJOB+10
   POSNEG=SIGN(1.D0,XEND-X)
   HMAXN=MIN(ABS(HMAX),ABS(XEND-X))
   IF (ABS(H).LE.10.D0*UROUND) H=1.0D-6
   H=MIN(ABS(H),HMAXN)
   H=SIGN(H,POSNEG)
   HOLD=H
   REJECT=.FALSE.
   FIRST=.TRUE.
   LAST=.FALSE.
   IF ((X+H*1.0001D0-XEND)*POSNEG.GE.0.D0) THEN
      H=XEND-X
      LAST=.TRUE.
   END IF
   HOPT=H
   FACCON=1.D0
   CFAC=SAFE*(1+2*NIT)
   NSING=0
   XOLD=X
   IF (IOUT.NE.0) THEN
      IRTRN=1
      NRSOL=1
      XOSOL=XOLD
      XSOL=X
      DO I=1,N
         CONT(I)=Y(I)
      END DO
      NSOLU=N
      HSOL=HOLD
      CALL SOLOUT(NRSOL,XOSOL,XSOL,Y,CONT,LRC,NSOLU,&
      &RPAR,IPAR,IRTRN)
      IF (IRTRN.LT.0) GOTO 179
   END IF
   MLE=MLJAC
   MUE=MUJAC
   MBJAC=MLJAC+MUJAC+1
   MBB=MLMAS+MUMAS+1
   MDIAG=MLE+MUE+1
   MDIFF=MLE+MUE-MUMAS
   MBDIAG=MUMAS+1
   N2=2*N
   N3=3*N
   IF (ITOL.EQ.0) THEN
      DO I=1,N
         SCAL(I)=ATOL(1)+RTOL(1)*ABS(Y(I))
      END DO
   ELSE
      DO I=1,N
         SCAL(I)=ATOL(I)+RTOL(I)*ABS(Y(I))
      END DO
   END IF
   HHFAC=H
   CALL FCN(N,X,Y,Y0,RPAR,IPAR)
   NFCN=NFCN+1
! --- BASIC INTEGRATION STEP
10 CONTINUE
! *** *** *** *** *** *** ***
!  COMPUTATION OF THE JACOBIAN
! *** *** *** *** *** *** ***
   NJAC=NJAC+1
   IF (IJAC.EQ.0) THEN
! --- COMPUTE JACOBIAN MATRIX NUMERICALLY
      IF (BANDED) THEN
! --- JACOBIAN IS BANDED
         MUJACP=MUJAC+1
         MD=MIN(MBJAC,M2)
         DO MM=1,M1/M2+1
            DO K=1,MD
               J=K+(MM-1)*M2
12             F1(J)=Y(J)
               F2(J)=DSQRT(UROUND*MAX(1.D-5,ABS(Y(J))))
               Y(J)=Y(J)+F2(J)
               J=J+MD
               IF (J.LE.MM*M2) GOTO 12
               CALL FCN(N,X,Y,CONT,RPAR,IPAR)
               J=K+(MM-1)*M2
               J1=K
               LBEG=MAX(1,J1-MUJAC)+M1
14             LEND=MIN(M2,J1+MLJAC)+M1
               Y(J)=F1(J)
               MUJACJ=MUJACP-J1-M1
               DO L=LBEG,LEND
                  FJAC(L+MUJACJ,J)=(CONT(L)-Y0(L))/F2(J)
               END DO
               J=J+MD
               J1=J1+MD
               LBEG=LEND+1
               IF (J.LE.MM*M2) GOTO 14
            END DO
         END DO

         ! Convert jacobian into sparse form
         sparseJ = FJAC(1:N,1:N)

      ELSE
! --- JACOBIAN IS FULL
         DO I=1,N
            YSAFE=Y(I)
            DELT=DSQRT(UROUND*MAX(1.D-5,ABS(YSAFE)))
            Y(I)=YSAFE+DELT
            CALL FCN(N,X,Y,CONT,RPAR,IPAR)
            DO J=M1+1,N
               FJAC(J-M1,I)=(CONT(J)-Y0(J))/DELT
            END DO
            Y(I)=YSAFE
         END DO

         ! Convert jacobian into sparse form
         sparseJ = FJAC(1:N,1:N)

      END IF
   ELSE
! --- COMPUTE JACOBIAN MATRIX ANALYTICALLY
!         CALL JAC(N,X,Y,FJAC,LDJAC,RPAR,IPAR)
!         write(*,*)'Analytical jacobian'
      Xdp = real (X, dp)
      Ydp = real (Y, dp)
      call jac(n,Xdp,Ydp,sparseJ)

   END IF
   CALJAC=.TRUE.
   CALHES=.TRUE.
20 CONTINUE
! --- COMPUTE THE MATRICES E1 AND E2 AND THEIR DECOMPOSITIONS
   if (H/=0.d0) then
      FAC1=U1/H
      ALPHN=ALPH/H
      BETAN=BETA/H
   else
      FAC1=huge(0.d0)
      ALPHN=huge(0.d0)
      BETAN=huge(0.d0)
   endif

!     Compute matrix E1 and decompose
!      CALL DECOMR(N,FJAC,LDJAC,FMAS,LDMAS,MLMAS,MUMAS,
!     &            M1,M2,NM1,FAC1,E1,LDE1,IP1,IER,IJOB,CALHES,IPHES)

!      msys1 =
!      msys2 =


   R5_sys1 = identity(n, real(fac1,dp)) - sparseJ
   R5_sys2 = system_2_matrix(alphn,betan,sparseJ)
!      R5_sys1%A = msys1%A
!      R5_sys2%A = msys2%A

   call sparseLU(R5_sys1)
   call sparseLU(R5_sys2)

!     Compute LU factorization of the two matrices, numerically
!      call numerical_factorization(msys1,iper1,iiper1,iLU1,rLU1)
!      call numerical_factorization(msys2,iper2,iiper2,iLU2,rLU2)

!          Numerical factorization of the r5_newton_matrix
!      call numerical_factorization(sparse_E1,iperm,iinvperm,iLUwork,
!     &                             rLUwork)



!      sparse_E1 = sparse_block_diagonal(E1_L, E1_U)

!      call print_sparse_to_file(sparse_E1,'spE1.dat')


!      call print_sparse_to_file(sparse_E1,'newton_mat.dat')

!      stop



!      IF (IER.NE.0) GOTO 78
!      CALL DECOMC(N,FJAC,LDJAC,FMAS,LDMAS,MLMAS,MUMAS,
!     &            M1,M2,NM1,ALPHN,BETAN,E2R,E2I,LDE1,IP2,IER,IJOB)
   IF (IER.NE.0) GOTO 78
   NDEC=NDEC+1
30 CONTINUE
   NSTEP=NSTEP+1
   IF (NSTEP.GT.NMAX) GOTO 178
   IF (0.1D0*ABS(H).LE.ABS(X)*UROUND.AND.X+H<XEND) GOTO 177
   IF (INDEX2) THEN
      DO I=NIND1+1,NIND1+NIND2
         SCAL(I)=SCAL(I)/HHFAC
      END DO
   END IF
   IF (INDEX3) THEN
      DO I=NIND1+NIND2+1,NIND1+NIND2+NIND3
         SCAL(I)=SCAL(I)/(HHFAC*HHFAC)
      END DO
   END IF
   XPH=X+H
! *** *** *** *** *** *** ***
!  STARTING VALUES FOR NEWTON ITERATION
! *** *** *** *** *** *** ***
   IF (FIRST.OR.STARTN) THEN
      DO I=1,N
         Z1(I)=0.D0
         Z2(I)=0.D0
         Z3(I)=0.D0
         F1(I)=0.D0
         F2(I)=0.D0
         F3(I)=0.D0
      END DO
   ELSE
      C3Q=H/HOLD
      C1Q=C1*C3Q
      C2Q=C2*C3Q
      DO I=1,N
         AK1=CONT(I+N)
         AK2=CONT(I+N2)
         AK3=CONT(I+N3)
         Z1I=C1Q*(AK1+(C1Q-C2M1)*(AK2+(C1Q-C1M1)*AK3))
         Z2I=C2Q*(AK1+(C2Q-C2M1)*(AK2+(C2Q-C1M1)*AK3))
         Z3I=C3Q*(AK1+(C3Q-C2M1)*(AK2+(C3Q-C1M1)*AK3))
         Z1(I)=Z1I
         Z2(I)=Z2I
         Z3(I)=Z3I
         F1(I)=TI11*Z1I+TI12*Z2I+TI13*Z3I
         F2(I)=TI21*Z1I+TI22*Z2I+TI23*Z3I
         F3(I)=TI31*Z1I+TI32*Z2I+TI33*Z3I
      END DO
   END IF
! *** *** *** *** *** *** ***
!  LOOP FOR THE SIMPLIFIED NEWTON ITERATION
! *** *** *** *** *** *** ***
   NEWT=0
   FACCON=MAX(FACCON,UROUND)**0.8D0
   THETA=ABS(THET)
40 CONTINUE
   IF (NEWT.GE.NIT) GOTO 78
! ---     COMPUTE THE RIGHT-HAND SIDE
   DO I=1,N
      CONT(I)=Y(I)+Z1(I)
   END DO
   CALL FCN(N,X+C1*H,CONT,Z1,RPAR,IPAR)
   DO I=1,N
      CONT(I)=Y(I)+Z2(I)
   END DO
   CALL FCN(N,X+C2*H,CONT,Z2,RPAR,IPAR)
   DO I=1,N
      CONT(I)=Y(I)+Z3(I)
   END DO
   CALL FCN(N,XPH,CONT,Z3,RPAR,IPAR)
   NFCN=NFCN+3
! ---     SOLVE THE LINEAR SYSTEMS
   DO I=1,N
      A1=Z1(I)
      A2=Z2(I)
      A3=Z3(I)
      Z1(I)=TI11*A1+TI12*A2+TI13*A3
      Z2(I)=TI21*A1+TI22*A2+TI23*A3
      Z3(I)=TI31*A1+TI32*A2+TI33*A3
   END DO

!          prepare array of the unknowns
   call rad_rhs(N,F1,F2,F3,FAC1,ALPHN,BETAN,Z1,Z2,Z3)

   unk1 = Z1
   unk2 = [Z2,Z3]

!          res1 = msys1.backslash.unk1
!          res2 = msys2.backslash.unk2

   res1 = R5_sys1.backslash.unk1
   res2 = R5_sys2.backslash.unk2

!        CALL SLVRAD(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
!     &          M1,M2,NM1,FAC1,ALPHN,BETAN,E1,E2R,E2I,LDE1,Z1,Z2,Z3,
!     &          F1,F2,F3,CONT,IP1,IP2,IPHES,IER,IJOB)

   Z1 = res1(1:n)
   Z2 = res2(1:n)
   Z3 = res2(n+1:2*n)


   NSOL=NSOL+1
   NEWT=NEWT+1
   DYNO=0.D0
   DO I=1,N
      DENOM=SCAL(I)
      DYNO=DYNO+(Z1(I)/DENOM)**2+(Z2(I)/DENOM)**2&
      &+(Z3(I)/DENOM)**2
   END DO
   DYNO=DSQRT(DYNO/N3)
! ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE
   IF (NEWT.GT.1.AND.NEWT.LT.NIT) THEN
      THQ=DYNO/DYNOLD
      IF (NEWT.EQ.2) THEN
         THETA=THQ
      ELSE
         THETA=SQRT(THQ*THQOLD)
      END IF
      THQOLD=THQ
      IF (THETA.LT.0.99D0) THEN
         FACCON=THETA/(1.0D0-THETA)
         DYTH=FACCON*DYNO*THETA**(NIT-1-NEWT)/FNEWT
         IF (DYTH.GE.1.0D0) THEN
            QNEWT=DMAX1(1.0D-4,DMIN1(20.0D0,DYTH))
            HHFAC=.8D0*QNEWT**(-1.0D0/(4.0D0+NIT-1-NEWT))
            H=HHFAC*H
            REJECT=.TRUE.
            LAST=.FALSE.
            IF (CALJAC) GOTO 20
            GOTO 10
         END IF
      ELSE
         GOTO 78
      END IF
   END IF
   DYNOLD=MAX(DYNO,UROUND)
   DO I=1,N
      F1I=F1(I)+Z1(I)
      F2I=F2(I)+Z2(I)
      F3I=F3(I)+Z3(I)
      F1(I)=F1I
      F2(I)=F2I
      F3(I)=F3I
      Z1(I)=T11*F1I+T12*F2I+T13*F3I
      Z2(I)=T21*F1I+T22*F2I+T23*F3I
      Z3(I)=T31*F1I+    F2I
   END DO
   IF (FACCON*DYNO.GT.FNEWT) GOTO 40 ! End of newton iteration
! --- ERROR ESTIMATION
   CALL ESTRAD (N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,&
   &H,DD1,DD2,DD3,FCN,NFCN,Y0,Y,IJOB,X,M1,M2,NM1,&
   &E1,LDE1,Z1,Z2,Z3,CONT,F1,F2,IP1,IPHES,SCAL,ERR,&
   &FIRST,REJECT,FAC1,RPAR,IPAR)

! --- COMPUTATION OF HNEW
! --- WE REQUIRE .2<=HNEW/H<=8.
   FAC=MIN(SAFE,CFAC/(NEWT+2*NIT))
   QUOT=MAX(FACR,MIN(FACL,ERR**.25D0/FAC))
   HNEW=H/QUOT
! *** *** *** *** *** *** ***
!  IS THE ERROR SMALL ENOUGH ?
! *** *** *** *** *** *** ***
   IF (ERR.LT.1.D0) THEN
! --- STEP IS ACCEPTED
      FIRST=.FALSE.
      NACCPT=NACCPT+1
      IF (PRED) THEN
!       --- PREDICTIVE CONTROLLER OF GUSTAFSSON
         IF (NACCPT.GT.1) THEN
            FACGUS=(HACC/H)*(ERR**2/ERRACC)**0.25D0/SAFE
            FACGUS=MAX(FACR,MIN(FACL,FACGUS))
            QUOT=MAX(QUOT,FACGUS)
            HNEW=H/QUOT
         END IF
         HACC=H
         ERRACC=MAX(1.0D-2,ERR)
      END IF
      XOLD=X
      HOLD=H
      X=XPH
      DO I=1,N
         Y(I)=Y(I)+Z3(I)
         Z2I=Z2(I)
         Z1I=Z1(I)
         CONT(I+N)=(Z2I-Z3(I))/C2M1
         AK=(Z1I-Z2I)/C1MC2
         ACONT3=Z1I/C1
         ACONT3=(AK-ACONT3)/C2
         CONT(I+N2)=(AK-CONT(I+N))/C1M1
         CONT(I+N3)=CONT(I+N2)-ACONT3
      END DO
      IF (ITOL.EQ.0) THEN
         DO I=1,N
            SCAL(I)=ATOL(1)+RTOL(1)*ABS(Y(I))
         END DO
      ELSE
         DO I=1,N
            SCAL(I)=ATOL(I)+RTOL(I)*ABS(Y(I))
         END DO
      END IF
      IF (IOUT.NE.0) THEN
         NRSOL=NACCPT+1
         XSOL=X
         XOSOL=XOLD
         DO I=1,N
            CONT(I)=Y(I)
         END DO
         NSOLU=N
         HSOL=HOLD
         CALL SOLOUT(NRSOL,XOSOL,XSOL,Y,CONT,LRC,NSOLU,&
         &RPAR,IPAR,IRTRN)
         IF (IRTRN.LT.0) GOTO 179
      END IF
      CALJAC=.FALSE.
      IF (LAST) THEN
         H=HOPT
         IDID=1
         call deallocate(msys1)
         call deallocate(msys2)
         call deallocate(sparseJ)
         RETURN
      END IF
      CALL FCN(N,X,Y,Y0,RPAR,IPAR)
      NFCN=NFCN+1
      HNEW=POSNEG*MIN(ABS(HNEW),HMAXN)
      HOPT=HNEW
      HOPT=MIN(H,HNEW)
      IF (REJECT) HNEW=POSNEG*MIN(ABS(HNEW),ABS(H))
      REJECT=.FALSE.
      IF ((X+HNEW/QUOT1-XEND)*POSNEG.GE.0.D0) THEN
         H=XEND-X
         LAST=.TRUE.
      ELSE
         QT=HNEW/H
         HHFAC=H
         IF (THETA.LE.THET.AND.QT.GE.QUOT1.AND.QT.LE.QUOT2) GOTO 30
         H=HNEW
      END IF
      HHFAC=H
      IF (THETA.LE.THET) GOTO 20
      GOTO 10
   ELSE
! --- STEP IS REJECTED
      REJECT=.TRUE.
      LAST=.FALSE.
      IF (FIRST) THEN
         H=H*0.1D0
         HHFAC=0.1D0
      ELSE
         HHFAC=HNEW/H
         H=HNEW
      END IF
      IF (NACCPT.GE.1) NREJCT=NREJCT+1
      IF (CALJAC) GOTO 20
      GOTO 10
   END IF
! --- UNEXPECTED STEP-REJECTION
78 CONTINUE
   IF (IER.NE.0) THEN
      NSING=NSING+1
      IF (NSING.GE.5) GOTO 176
   END IF
   H=H*0.5D0
   HHFAC=0.5D0
   REJECT=.TRUE.
   LAST=.FALSE.
   IF (CALJAC) GOTO 20
   GOTO 10
! --- FAIL EXIT
176 CONTINUE
   WRITE(6,979)X
   WRITE(6,*) ' MATRIX IS REPEATEDLY SINGULAR, IER=',IER
   IDID=-4
   call deallocate(msys1)
   call deallocate(msys2)
   call deallocate(sparseJ)
   RETURN
177 CONTINUE
   WRITE(6,979)X
   WRITE(6,*) ' STEP SIZE T0O SMALL, H=',H
   IDID=-3
   call deallocate(msys1)
   call deallocate(msys2)
   call deallocate(sparseJ)
   RETURN
178 CONTINUE
   WRITE(6,979)X
   WRITE(6,*) ' MORE THAN NMAX =',NMAX,'STEPS ARE NEEDED'
   IDID=-2
   call deallocate(msys1)
   call deallocate(msys2)
   call deallocate(sparseJ)
   RETURN
! --- EXIT CAUSED BY SOLOUT
179 CONTINUE
   WRITE(6,979)X
979 FORMAT(' EXIT OF RADAU5 AT X=',E22.16)
   IDID=2
   call deallocate(msys1)
   call deallocate(msys2)
   call deallocate(sparseJ)
   RETURN
END
!
!     END OF SUBROUTINE RADCOR
!
! ***********************************************************
!
DOUBLE PRECISION FUNCTION CONTR5(I,X,CONT,LRC)
! ----------------------------------------------------------
!     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT. IT PROVIDES AN
!     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT X.
!     IT GIVES THE VALUE OF THE COLLOCATION POLYNOMIAL, DEFINED FOR
!     THE LAST SUCCESSFULLY COMPUTED STEP (BY RADAU5).
! ----------------------------------------------------------
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
   DIMENSION CONT(LRC)
   COMMON /CONRA5/NN,NN2,NN3,NN4,XSOL,HSOL,C2M1,C1M1
!$ OMP THREADPRIVATE(/CONRA5/)
   S=(X-XSOL)/HSOL
   CONTR5=CONT(I)+S*(CONT(I+NN)+(S-C2M1)*(CONT(I+NN2)&
   &+(S-C1M1)*CONT(I+NN3)))
   RETURN
END
!
!     END OF FUNCTION CONTR5
!
! ***********************************************************


!     ******************************************************************
!      Prepare right-hand side columns of the linear system
!     ******************************************************************

subroutine rad_rhs(N,F1,F2,F3,FAC1,ALPHN,BETAN,Z1,Z2,Z3)
   implicit none

   integer, intent(in) :: n

   double precision, dimension(n), intent(in)  :: F1, F2, F3
   double precision,               intent(in)  :: FAC1, ALPHN, BETAN

   double precision, dimension(n), intent(inout) :: Z1, Z2, Z3

   double precision :: s2, s3
   integer          :: i

   do I=1,N
      S2=-F2(I)
      S3=-F3(I)
      Z1(I)=Z1(I)-F1(I)*FAC1
      Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
      Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
   end do

   return
end subroutine rad_rhs

!     ******************************************************************

!     ******************************************************************
!      Build Newton iteration's matrix (dimension: 3n x 3n)
!      N = [ (gamma I - J)                               ]
!          [               (alpha I - J)   (-beta I)     ]
!          [               (beta I)        (alpha I - J) ]
!
!     ******************************************************************

function r5_newton_matrix(n,alpha,beta,gamma,jac) result(mat)

   use working_precision, only: dp
   use sparse_algebra, only: sparse, identity, sparse_sum,&
   &sparse_block_vertical,&
   &sparse_block_horizontal,&
   &sparse_block_diagonal,&
   &print_sparse_to_file,&
   &sparse_nullify_general
   use sparse_definitions

   implicit none

   integer,          intent(in) :: n
   double precision, intent(in) :: alpha, beta, gamma

   type(sparse), intent(in) :: jac
   type(sparse)             :: mat
   type(sparse)             :: gamma_id, alpha_id, beta_id, betam_id,&
   &AIJ, B1, B2, B2L, B2R

!        Build identity matrices
   gamma_id = identity(n,real(gamma,dp))
   alpha_id = identity(n,real(alpha,dp))
   beta_id  = identity(n,real(beta ,dp))
   betam_id = identity(n,real(-beta,dp))

!        Build upper block

   B1  = gamma_id - jac
   call sparse_nullify_general(gamma_id)

   AIJ = alpha_id - jac
   call sparse_nullify_general(alpha_id)

   B2L = sparse_block_vertical(AIJ,  beta_id)
   B2R = sparse_block_vertical(betam_id, AIJ)
   call sparse_nullify_general(AIJ)
   call sparse_nullify_general(betam_id)
   call sparse_nullify_general(beta_id)

   B2  = sparse_block_horizontal(B2L, B2R)
   call sparse_nullify_general(B2R)
   call sparse_nullify_general(B2L)

   mat = sparse_block_diagonal(B1, B2)
   call sparse_nullify_general(B2)
   call sparse_nullify_general(B1)

!         call print_sparse_to_file(B2,'B2.dat')
!         call print_sparse_to_file(B2R,'B2R.dat')
!         call print_sparse_to_file(B2L,'B2L.dat')
!         call print_sparse_to_file(AIJ,'AIJ.dat')

end function r5_newton_matrix

!     ******************************************************************



!     ******************************************************************
!      Preprocessing of Jacobian matrix sparsity data
!     ******************************************************************

subroutine jacobian_prepro(neq,jacf,Y,RTOL,ATOL,a,b,g,iperm,&
&iinvperm,iLUwork,rLUwork)

   use sparse_algebra, only: sparse, sparse_square_permutation,&
   &print_sparse_to_file
   use sparse_definitions
   implicit none

   interface

      function perturbed(Y,RTOL,ATOL) result(perturbedY)
         implicit none
         double precision, dimension(:), intent(in) :: Y
         double precision, dimension(:), intent(in) :: RTOL, ATOL
         double precision, dimension(size(Y))       :: perturbedY
      end function perturbed

      subroutine jacf(neq,time,yin,jac)
         use sparse_algebra, only: sparse
         implicit none
         integer,          intent(in) :: neq
         double precision, intent(in) :: time
         double precision, dimension(neq), intent(in)  :: yin
         type(sparse),                     intent(out) :: jac
      end subroutine jacf

      function r5_newton_matrix(n,a,b,g,jac) result(mat)
         use sparse_algebra
         implicit none
         integer,          intent(in) :: n
         double precision, intent(in) :: a, b, g
         type(sparse), intent(in) :: jac
         type(sparse)             :: mat
      end function r5_newton_matrix


   end interface

   integer,                          intent(in) :: neq
   double precision, dimension(neq), intent(in) :: Y
   double precision, dimension(:),   intent(in) :: RTOL,ATOL
   integer, dimension(3*neq),        intent(out):: iperm, iinvperm
   integer, dimension(:),            intent(inout) :: iLUwork
   double precision, dimension(:),   intent(inout) :: rLUwork

!     Real parameters for RADAU method alpha, beta, u1
!     (assuming step h = 1.d0)
   double precision,                 intent(in) :: a, b, g

   type(sparse) :: jactmp, newtmat

   integer                                      :: lwork, error_flag,&
   &freespace, liwork
   integer,          dimension(:),  allocatable :: iwork
   double precision, dimension(neq)             :: pertY
   double precision, dimension(:),  allocatable :: rwork, z
   double precision                             :: t = 0.d0
   integer      :: j

!     ** Determine Jacobian sparsity structure by adopting a perturbed
!     ** array of the unknowns
   pertY = perturbed(Y,RTOL,ATOL)

!     ** Call jacobian matrix; the sparsity structure is contained
!     ** in jacf%IA, jacf%JA
   call jacf(neq,t,pertY,jactmp)

!     ** Build Newton iteration matrix
   newtmat = r5_newton_matrix(neq,a,b,g,jactmp)
   call deallocate(jactmp)

!     ** Retrieve newton-iterat matrix minimal dimension rows/columns
!     ** ordering for faster and more sparse linear system solution
!     ** involving this matrix

   lwork = 5 * (3*newtmat%nr + 4 * newtmat%n)
   allocate(iwork(lwork))
   iperm    = 0
   iinvperm = 0

   call odrv(newtmat%nr, newtmat%IA, newtmat%col, newtmat%val,&
   &iperm, iinvperm, lwork, iwork, 1, error_flag)

   deallocate(iwork)

   if (error_flag/=0) then
      write(*,*)'Error detected in call to odrv, ierror=',error_flag
      stop
   endif

!     ** Reorder newton iterations matrix elements and perform
!     ** symbolic LU factorisation (rows and columns have same ordering)
!     ** solution array is z (output), that overwrites z
!     ** rhs coefficients array is z (input), intialised as zeroes

   lwork = 2 * (8*(newtmat%nr + 2) + 2*newtmat%n)
   liwork = 2 * lwork

   if (size(iLUwork) < liwork) then
      write(*,*)'iLUwork'
      write(*,*)'Insufficient space for LU matrix factorization'
      write(*,*)'Provided: ',size(iLUwork),' Needed: ',liwork
      stop
   endif

   allocate(iwork(liwork),rwork(lwork))
   iwork   = 0
   rwork   = 0.d0

   allocate(z(newtmat%nr))
   z = 0.d0

   error_flag = 0

   call cdrv(newtmat%nr, iperm, iperm, iinvperm,&
   &newtmat%IA, newtmat%JA, newtmat%A, z, z,&
   &lwork, iwork, rwork, freespace, 5, error_flag)

   if (error_flag/=0 .or. freespace < 0) then
      write(*,*)'cdrv exited with error flag ',error_flag
      write(*,*)'working storage free space: ',freespace
      stop
   endif

   iLUwork(1:liwork)  = iwork
   iLUwork(liwork+1:) = -1

   rLUwork(1:lwork )   = rwork
   rLUwork(lwork+1:)   = 0.d0

   call deallocate(newtmat)

!      subroutine cdrv
!     *     (n, r,c,ic, ia,ja,a, b, z, nsp,isp,rsp,esp, path, flag)
!
!      CALL CDRV (N,IWK(IPR),IWK(IPC),IWK(IPIC),IWK(IPIAN),IWK(IPJAN),
!     1   WK(IPA),WK(IPA),WK(IPA),NSP,IWK(IPISP),WK(IPRSP),IESP,5,IYS)


end subroutine jacobian_prepro

!     ******************************************************************
!       Preprocessing matrices for the solution of 1st system of lin.eq.
!     ******************************************************************

subroutine system_1_prepro(neq,jacmat,gam,ip1,iip1,iLU1,rLU1)

   use sparse_algebra, only: sparse, sparse_square_permutation,&
   &identity, sparse_sum,&
   &sparse_nullify_general,&
   &print_sparse_to_file

   use sparse_definitions
   use working_precision
   implicit none

   interface
      function yale_work_dimensions(n,n2) result(lwork)
         implicit none
         integer, intent(in) :: n
         integer, intent(in), optional :: n2
         integer             :: lwork
      end function yale_work_dimensions
   end interface


   integer,                          intent(in) :: neq
   type(sparse),                     intent(in) :: jacmat
   double precision,                 intent(in) :: gam
   integer, dimension(neq),          intent(out):: ip1, iip1
   integer, dimension(:),            intent(inout), target :: iLU1
   double precision, dimension(:),   intent(inout), target :: rLU1

   type(sparse) :: sys_matrix

   integer                                      :: lwork, error_flag,&
   &freespace, liwork
   integer,          dimension(:),  pointer     :: iwork
   double precision, dimension(neq)             :: pertY
   double precision, dimension(:),  pointer     :: rwork
   double precision, dimension(:),  allocatable ::  z
   integer      :: j

   type(sparse) :: tmp

!     ** Build system matrix
   sys_matrix = identity(neq,real(gam,dp)) - jacmat

!     ** Retrieve matrix minimal dimension rows/columns
!     ** ordering for faster and more sparse linear system solution
!     ** involving this matrix

   lwork = 5 * (3*sys_matrix%nr + 4 * sys_matrix%n)

   iwork => iLU1(1:lwork)

   ip1    = 0
   iip1   = 0

   call odrv(sys_matrix%nr, sys_matrix%IA, sys_matrix%col,&
   &sys_matrix%val, ip1, iip1, lwork, iwork, 1, error_flag)

!      deallocate(iwork)

   if (error_flag/=0) then
      write(*,*)'Error detected in call to odrv, ierror=',error_flag
      stop
   endif

!     ** Reorder newton iterations matrix elements and perform
!     ** symbolic LU factorisation (rows and columns have same ordering)
!     ** solution array is z (output), that overwrites z
!     ** rhs coefficients array is z (input), intialised as zeroes

   lwork = yale_work_dimensions(neq, sys_matrix%n)
   liwork = 2 * lwork

   if (size(iLU1) < liwork) then
      write(*,*)'System 1 prepro'
      write(*,*)'Insufficient space for LU matrix factorization'
      write(*,*)'Provided: ',size(iLU1),' Needed: ',liwork
      stop
   endif

   nullify(iwork)
   iwork  => iLU1(1:liwork)
   rwork  => rLU1(1:lwork )
!      iwork   = 0
!      rwork   = 0.d0

   allocate(z(sys_matrix%nr))
   z = 0.d0

   error_flag = 0

   call cdrv(sys_matrix%nr, ip1, ip1, iip1,&
   &sys_matrix%IA, sys_matrix%JA, sys_matrix%A, z, z,&
   &lwork, iwork, rwork, freespace, 5, error_flag)

   if (error_flag/=0 .or. freespace < 0) then
      write(*,*)'System 1 prepro'
      write(*,*)'cdrv exited with error flag ',error_flag
      write(*,*)'working storage free space: ',freespace
      stop
   endif

   deallocate(z)

   call sparse_nullify_general(sys_matrix)

!      iLU1(1:liwork)  = iwork
!      iLU1(liwork+1:) = -1

!      rLU1(1:lwork )   = rwork
!      rLU1(lwork+1:)   = 0.d0

end subroutine system_1_prepro


!     ******************************************************************
!       Preprocessing matrices for the solution of 1st system of lin.eq.
!     ******************************************************************

subroutine system_2_prepro(neq,jacmat,alfa,bet,ip2,iip2,iLU2,rLU2)

   use sparse_algebra, only: sparse, sparse_square_permutation,&
   &identity, sparse_sum,&
   &sparse_block_vertical,&
   &sparse_block_horizontal,&
   &sparse_nullify_general,&
   &print_sparse_to_file,&
   &print_sparsity_to_file
   use sparse_definitions
   use working_precision, only: dp
   implicit none

   interface
      function yale_work_dimensions(n,n2) result(lwork)
         implicit none
         integer, intent(in) :: n
         integer, intent(in), optional :: n2
         integer             :: lwork
      end function yale_work_dimensions
   end interface


   integer,                          intent(in) :: neq
   type(sparse),                     intent(in) :: jacmat
   double precision,                 intent(in) :: alfa,bet
   integer, dimension(2*neq),        intent(out):: ip2, iip2
   integer, dimension(:),            intent(inout), target :: iLU2
   double precision, dimension(:),   intent(inout), target :: rLU2

   type(sparse) :: m1, m2, partial_left, partial_right, sys_matrix

   integer                                      :: lwork, error_flag,&
   &freespace, liwork
   integer,          dimension(:),  pointer     :: iwork
   double precision, dimension(neq)             :: pertY
   double precision, dimension(:),  pointer     :: rwork
   double precision, dimension(:),  allocatable :: z

   integer      :: j


!     ** Build system matrix
   m1           = identity(neq,real(alfa,dp))  - jacmat
!      call print_sparse_to_file(m1, 'm1.dat')
!
!      m1           = identity(neq,alfa) - jacmat
!      call print_sparse_to_file(m1, 'm1b.dat')
!      stop
   m2           = identity(neq, real(bet,dp))
   partial_left = sparse_block_vertical(m1, m2)
   m2           = identity(neq, real(-bet,dp))
   partial_right= sparse_block_vertical(m2, m1)
   sys_matrix   = sparse_block_horizontal(partial_left,partial_right)

   call sparse_nullify_general(m1)
   call sparse_nullify_general(m2)
   call sparse_nullify_general(partial_left)
   call sparse_nullify_general(partial_right)

!     ** Retrieve matrix minimal dimension rows/columns
!     ** ordering for faster and more sparse linear system solution
!     ** involving this matrix

   lwork = 5 * (3*sys_matrix%nr + 4 * sys_matrix%n)

   iwork => iLU2(1:lwork)

   ip2    = 0
   iip2   = 0

   call odrv(sys_matrix%nr, sys_matrix%IA, sys_matrix%col,&
   &sys_matrix%val, ip2, iip2, lwork, iwork, 1, error_flag)

   nullify(iwork)

   if (error_flag/=0) then
      write(*,*)'Error detected in call to odrv, ierror=',error_flag
      stop
   endif

!     ** Reorder newton iterations matrix elements and perform
!     ** symbolic LU factorisation (rows and columns have same ordering)
!     ** solution array is z (output), that overwrites z
!     ** rhs coefficients array is z (input), intialised as zeroes

   lwork = yale_work_dimensions(2*neq, 2*jacmat%nr + 4*neq)
   liwork = 2 * lwork

   if (size(iLU2) < liwork) then
      write(*,*)'iLU2'
      write(*,*)'Insufficient space for LU matrix factorization'
      write(*,*)'Provided: ',size(iLU2),' Needed: ',liwork
      stop
   endif

!      allocate(iwork(liwork),rwork(lwork))
   iwork => iLU2(1:liwork)
   rwork => rLU2(1:lwork )
!      iwork   = 0
!      rwork   = 0.d0

   allocate(z(sys_matrix%nr))
   z = 0.d0

   error_flag = 0

   call cdrv(sys_matrix%nr, ip2, ip2, iip2,&
   &sys_matrix%IA, sys_matrix%JA, sys_matrix%A, z, z,&
   &lwork, iwork, rwork, freespace, 5, error_flag)

   if (error_flag/=0 .or. freespace < 0) then
      write(*,*)'System 2 prepro'
      write(*,*)'cdrv exited with error flag ',error_flag
      write(*,*)'working storage free space: ',freespace
      stop
   endif

   deallocate(z)

   call sparse_nullify_general(sys_matrix)

!      iLU2(1:liwork)  = iwork
!      iLU2(liwork+1:) = -1

!      rLU2(1:lwork )   = rwork
!      rLU2(lwork+1:)   = 0.d0

end subroutine system_2_prepro


!     ******************************************************************
!      Perturbation of the array of the unknowns by a small amount
!     ******************************************************************

function perturbed(Y,RTOL,ATOL) result(perturbedY)
   implicit none

   interface
      function scale_array(Y,RTOL,ATOL) result(scal)
         implicit none
         double precision, dimension(:), intent(in) :: Y
         double precision, dimension(:), intent(in) :: RTOL, ATOL
         double precision, dimension(size(Y))       :: scal
      end function scale_array
   end interface

   double precision, dimension(:), intent(in) :: Y
   double precision, dimension(:), intent(in) :: RTOL, ATOL
   double precision, dimension(size(Y))       :: perturbedY
   double precision, dimension(size(Y))       :: perturbation

   double precision, dimension(size(Y))       :: Yscale

   integer :: n, j


!     Determine array dimensions
   n = size(Y)

!     Determine scale array
   Yscale = scale_array(Y,RTOL,ATOL)

!     Setup uneven perturbation factor
   perturbation = 1.d0 + 1.d0/(dble([(j,j=1,n)]) + 1.d0)

!     Perturb array y
   perturbedY = Y + perturbation * sign(Yscale,Y)

end function perturbed


!     ******************************************************************
!      Determine the scale array based on tolerance constraints
!      scale(i) = RTOL*ABS(Y(i)) + ATOL(i)
!     ******************************************************************

function scale_array(Y,RTOL,ATOL) result(scal)
   implicit none

   double precision, dimension(:), intent(in) :: Y
   double precision, dimension(:), intent(in) :: RTOL, ATOL

   double precision, dimension(size(Y))       :: scal

   integer :: n


!     Determine array dimensions
   n = size(Y)

   if (size(RTOL)>1 .and. size(RTOL)/=n) then
      write(*,*)'Error: RTOL dimensions not 1 or neq'
      stop
   endif

   if (size(ATOL)>1 .and. size(ATOL)/=n) then
      write(*,*)'Error: ATOL dimensions not 1 or neq'
      stop
   endif

   scal = RTOL * ABS(Y) + ATOL

end function scale_array

!     ******************************************************************
!      Compute numerical factorization of a square sparse matrix

!       IP  = index of row+columns permutations (dimension: nrows)
!       IIP = inverted index of permutations
!       iLUwork = working array processed by subroutine odrv
!     ******************************************************************

subroutine numerical_factorization(spmat,iperm,iinvperm,&
&iLUwork,rLUwork)
   use sparse_algebra, only: sparse,print_sparse_to_file
   implicit none

   interface
      function yale_work_dimensions(n,n2) result(lwork)
         implicit none
         integer, intent(in) :: n
         integer, intent(in), optional :: n2
         integer             :: lwork
      end function yale_work_dimensions
   end interface

   type(sparse), intent(inout) :: spmat
   integer, dimension(:),          intent(inout) :: iperm, iinvperm
   integer, dimension(:),          intent(inout), target :: iLUwork
   double precision, dimension(:), intent(inout), target :: rLUwork

   double precision, dimension(spmat%nr) :: tmp
   integer :: lwork, liwork, freespace, error_flag, j
   integer,          dimension(:), pointer :: iwk
   double precision, dimension(:), pointer :: rwk

!     Dimensional check
   if (spmat%nr /= spmat%nc) then
      write(*,*)'Matrix must be square for LU decomp!'
      write(*,*)'nrows = ',spmat%nr,' ncols = ',spmat%nc
      stop
   endif

   if (size(iperm)/=spmat%nr .or. size(iinvperm)/=spmat%nr) then
      write(*,*)'Wrong pivot array dimensions in LU decomp'
      write(*,*)'nrows = ',spmat%nr,' len(ip)=',size(iperm)
      stop
   endif

   tmp(1:spmat%nr) = 0.d0

   lwork = yale_work_dimensions(spmat%nr,spmat%n)!2 * (8*(spmat%nr + 2) + 2*spmat%n)
   liwork = 2 * lwork


   if (size(iLUwork) < liwork) then
      write(*,*)'iLUwork, numerical_factorization'
      write(*,*)'Insufficient space for LU matrix factorization'
      write(*,*)'Provided: ',size(iLUwork),' Needed: ',liwork
      stop
   endif


   if (size(rLUwork) < lwork) then
      write(*,*)'Insufficient real space for LU matrix factorization'
      write(*,*)'Provided: ',size(rLUwork),' Needed: ',lwork
      stop
   endif

!      allocate(iwk(liwork), rwk(lwork))
   iwk => iLUwork(1:liwork)
   rwk => rLUwork(1:lwork )

   error_flag = 0

!     Call to the CDRV routine for numerical matrix factorization
   call cdrv(spmat%nr, iperm, iperm, iinvperm,&
   &spmat%IA, spmat%JA, spmat%A, tmp, tmp,&
   &lwork, iwk, rwk, freespace, 2, error_flag)

   if (error_flag/=0 .or. freespace < 0) then
      write(*,*)'numerical factorization'
      write(*,*)'cdrv exited with error flag ',error_flag
      write(*,*)'working storage free space: ',freespace
      write(*,*)'working storage provided  : ',lwork
      stop
   endif

!      iLUwork(1:liwork)  = iwk
!      iLUwork(liwork+1:) = -1

!      rLUwork(1:lwork ) = rwk
!      rLUwork(lwork+1:) = 0.d0


!      deallocate(iwk,rwk)

end subroutine numerical_factorization

!     ******************************************************************
!      Solve sparse linear system of equations: spmat x = b
!      NB: need to prior diagonalize the system matrix by call to cdrv

!       IP  = index of row+columns permutations (dimension: nrows)
!       IIP = inverted index of permutations
!       iLUwork = working array processed by subroutine odrv
!     ******************************************************************
function sparse_linear_system(spmat,b,iperm,iinvperm,iLUwork,&
&rLUwork) result(x)
   use sparse_algebra, only: sparse
   implicit none

   interface
      function yale_work_dimensions(n,n2) result(lwork)
         implicit none
         integer, intent(in) :: n
         integer, intent(in), optional :: n2
         integer             :: lwork
      end function yale_work_dimensions
   end interface

   type(sparse),                   intent(in) :: spmat
   integer, dimension(:),          intent(inout) :: iperm, iinvperm
   integer, dimension(:),          intent(inout), target :: iLUwork
   double precision, dimension(:), intent(inout), target :: rLUwork
   double precision, dimension(:), intent(in) :: b
   double precision, dimension(size(b))       :: x

   integer :: lwork, liwork, freespace, error_flag, j
   integer,          dimension(:), pointer :: iwk
   double precision, dimension(:), pointer :: rwk

   if (spmat%nr /= spmat%nc) then
      write(*,*)'Matrix must be square for LU decomp!'
      write(*,*)'nrows = ',spmat%nr,' ncols = ',spmat%nc
      stop
   endif

   if (size(iperm)/=spmat%nr .or. size(iinvperm)/=spmat%nr) then
      write(*,*)'Wrong pivot array dimensions in LU decomp'
      write(*,*)'nrows = ',spmat%nr,' len(ip)=',size(iperm)
      stop
   endif

   if (size(b) /= spmat%nr) then
      write(*,*)'Wrong coefficient array dimensions in lin_solve'
      write(*,*)'size(b)=',size(b), 'mat%nr = ',spmat%nr
      stop
   endif

!      lwork = 2 * (8*(spmat%nr + 2) + 2*spmat%n)
!      liwork = 2 * lwork
   lwork = yale_work_dimensions(spmat%nr,spmat%n)!2 * (8*(spmat%nr + 2) + 2*spmat%n)
   liwork = 2 * lwork

   if (size(iLUwork) < liwork) then
      write(*,*)'iLUwork, sparse_linear_system'
      write(*,*)'Insufficient space for LU matrix factorization'
      write(*,*)'Provided: ',size(iLUwork),' Needed: ',liwork
      stop
   endif

   if (size(rLUwork) < lwork) then
      write(*,*)'Insufficient real space for LU matrix factorization'
      write(*,*)'Provided: ',size(rLUwork),' Needed: ',lwork
      stop
   endif

!      allocate(iwk(liwork), rwk(lwork))
   iwk => iLUwork(1:liwork)
   rwk => rLUwork(1:lwork)

   error_flag = 0

!     Call to the CDRV routine for numerical matrix factorization
   call cdrv(spmat%nr, iperm, iperm, iinvperm,&
   &spmat%IA, spmat%JA, spmat%A, b, x,&
   &lwork, iwk, rwk, freespace, 3, error_flag)


!      write(*,*)error_flag, freespace

   if (error_flag/=0 .or. freespace < 0) then
      write(*,*)'sparse_linear_system'
      write(*,*)'cdrv exited with error flag ',error_flag
      write(*,*)'working storage free space: ',freespace
      stop
   endif

!      iLUwork(1:liwork)  = iwk
!      iLUwork(liwork+1:) = -1
!
!      rLUwork(1:lwork )  = rwk
!      rLUwork(lwork+1:) = 0.d0


end function sparse_linear_system


!     ******************************************************************
!      Compute working array dimensions for Yale sparse package
!      n = number of system unknowns
!      n2 = number of sparse matrix nonzero elements (optional)
!      safety coefficient: the LU factorization can destroy sparsity
!     ******************************************************************

function yale_work_dimensions(n,n2) result(lwork)
   implicit none

   integer, intent(in) :: n
   integer, intent(in), optional :: n2
   integer             :: lwork
   integer             :: safety_coef = 10

   if (present(n2)) then
      lwork = safety_coef * ( 9*(n+2) + 2 * n2   )
   else
      lwork = safety_coef * ( 9*(n+2) + 2 * n**2 )
   endif

end function yale_work_dimensions

!     ******************************************************************
!       Build linear system 2 matrix
!     ******************************************************************
function system_2_matrix(alfa,beta,jac) result(s2mat)
   use working_precision, only: dp
   use sparse_algebra, only: sparse, sparse_sum, identity,&
   &sparse_block_vertical,&
   &sparse_block_horizontal,&
   &print_sparse_to_file
   use sparse_definitions
   implicit none

   type(sparse),     intent(in) :: jac
   double precision, intent(in) :: alfa, beta
   type(sparse)                 :: s2mat

   integer                      :: neq
   type(sparse)                 :: m1, m2, part_l, part_r

   neq     = jac%nr
   m1      = identity(neq, real(alfa,dp)) - jac
   m2      = identity(neq, real(beta,dp))
   part_l  = sparse_block_vertical(m1, m2)
!      call matrix_details(part_l,'partl')

!      call print_sparse_to_file(part_l,'partl.dat')
!      stop
   m2      = identity(neq, real(-beta,dp))
   part_r  = sparse_block_vertical(m2, m1)
   s2mat   = sparse_block_horizontal(part_l, part_r)

   call deallocate(m1)
   call deallocate(m2)
   call deallocate(part_l)
   call deallocate(part_r)

!      call print_sparse_to_file(s2mat,'s2mat.dat')
!      stop

end function system_2_matrix

!     ******************************************************************
!       Prepare sparse algebra matrices allocation for call to
!       sparse RADAU5 solver
!     ******************************************************************

subroutine sparse_radau5_algebra_allocation(neq,njac)

   use working_precision, only: dp
   use sparse_chemistry, only: stoich_r_sp, stoich_p_sp
   use sparse_algebra, only: sparse, identity,&
   &print_sparse_to_file
   use sparse_definitions

   use ode_solver, only: iper1,iiper1,iper2,iiper2,iLU1,iLU2,rLU1,&
   &rLU2,liLU1,lrLU1,liLU2,lrLU2, R5_sys1,&
   &R5_sys2

   use speedchem_conv, only: constV_jac_RADAUS

   implicit none

   interface

      function system_2_matrix(alfa,beta,jac) result(s2mat)
         use sparse_algebra
         implicit none
         type(sparse),     intent(in) :: jac
         double precision, intent(in) :: alfa, beta
         type(sparse)                 :: s2mat
      end function system_2_matrix

      function yale_work_dimensions(n,n2) result(lwork)
         implicit none
         integer, intent(in) :: n
         integer, intent(in), optional :: n2
         integer             :: lwork
      end function yale_work_dimensions

      subroutine system_1_prepro(neq,jacmat,gam,ip1,iip1,iLU1,rLU1)
         use sparse_algebra
         implicit none
         integer,                          intent(in) :: neq
         type(sparse),                     intent(in) :: jacmat
         double precision,                 intent(in) :: gam
         integer, dimension(neq),          intent(out):: ip1, iip1
         integer, dimension(:),            intent(inout) :: iLU1
         double precision, dimension(:),   intent(inout) :: rLU1
      end subroutine system_1_prepro

      subroutine system_2_prepro(neq,jacmat,a,b,ip1,iip1,iLU1,rLU1)
         use sparse_algebra
         implicit none
         integer,                          intent(in) :: neq
         type(sparse),                     intent(in) :: jacmat
         double precision,                 intent(in) :: a,b
         integer, dimension(neq),          intent(out):: ip1, iip1
         integer, dimension(:),            intent(inout) :: iLU1
         double precision, dimension(:),   intent(inout) :: rLU1
      end subroutine system_2_prepro
   end interface

!      integer, parameter  :: dp = kind(0.d0)

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
   real (dp), parameter :: sq3   = 1.73205080756887729353d0
   real (dp), parameter :: sq6   = 2.44948974278317809820d0
   real (dp), parameter :: c1    = tenth * (four - sq6)
   real (dp), parameter :: c2    = tenth * (four + sq6)
   real (dp), parameter :: c1m1  = c1 - one
   real (dp), parameter :: c2m1  = c2 - one
   real (dp), parameter :: c1mc2 = c1 - c2
   real (dp), parameter :: dd1   = - third * (seven * sq6 + 13.d0)
   real (dp), parameter :: dd2   =   third * (seven * sq6 - 13.d0)
   real (dp), parameter :: dd3   = - third
   real (dp), parameter :: u1    = one/(third * tenth * (6.0d0 +&
   &81.d0**third - 9.d0**third))
   real (dp), parameter :: alp   = sixth * tenth * (12.d0 -&
   &81.d0**third + 9.d0**third)
   real (dp), parameter :: bet   = sixth * tenth * sq3 * (&
   &81.d0**third + 9.d0**third)
   real (dp), parameter :: cno   = alp**2 + bet**2
   real (dp), parameter :: alph  = alp/cno
   real (dp), parameter :: beta  = bet/cno


   integer, intent(in)           :: neq
   integer, intent(in), optional :: njac
   integer                       :: j

   type(sparse)                     :: sparseJ, tmpr, tmpp
   real (dp), dimension(neq) :: Y
   integer, dimension(neq-1) :: tmpip

!     No knowledge if arrays had already been allocated in this thread
   if (allocated(iper1 ))deallocate(iper1 )
   if (allocated(iiper1))deallocate(iiper1)
   if (allocated(iper2 ))deallocate(iper2 )
   if (allocated(iiper2))deallocate(iiper2)
   if (allocated(iLU1  ))deallocate(iLU1  )
   if (allocated(iLU2  ))deallocate(iLU2  )
   if (allocated(rLU1  ))deallocate(rLU1  )
   if (allocated(rLU2  ))deallocate(rLU2  )


!     Allocate working arrays for the solution of the two sparse
!     systems of linear equations
!      if (present(njac)) then
!         lrLU1 = yale_work_dimensions(neq,njac + neq)
!      else
!         lrLU1 = yale_work_dimensions(neq)
!      endif
!      liLU1 = 2 * lrLU1
!
!      if (present(njac)) then
!         lrLU2 = yale_work_dimensions(2*neq,2*njac + 4 * neq)
!      else
!         lrLU2 = yale_work_dimensions(2*neq)
!      endif
!      liLU2 = 2 * lrLU2
!
!      write(*,*)'RADAU sparse working array dimensions: ',liLU1,liLU2
!
!      allocate( iLU1(liLU1), rLU1(lrLU1) )
!      allocate( iLU2(liLU2), rLU2(lrLU2) )
!      allocate( iper1(neq)  , iiper1(neq)   )
!      allocate( iper2(2*neq), iiper2(2*neq) )
!
!      iLU1   = 0
!      iLU2   = 0
!      rLU1   = 0.d0
!      rLU2   = 0.d0
!      iper1  = 0
!      iper2  = 0
!      iiper1 = 0
!      iiper2 = 0

!     Preprocess sparse matrices for faster computation of the linear
!     system

!     1) Generate temporary Y for the problem
   Y = [1500._dp,(real(1.d0/real(neq-1,dp),dp),j=1,neq-1)]

!     2) Evaluate jacobian matrix
   call constV_jac_RADAUS(neq,zero,Y,sparseJ)



!     3) Preprocess first linear system
!      call system_1_prepro(neq,sparseJ,u1,iper1,iiper1,iLU1,rLU1)
   R5_sys1 = identity(neq,u1) - sparseJ
!      do j = 2, neq
!         if (iper1(j)>iper1(1)) then
!           tmpip(j-1) = iper1(j)-iper1(1) + 1
!         else
!           tmpip(j-1) = iper1(j)
!         endif
!      end do

!      open(unit=314,file='neworder.dat')
!      write(314,"(1x,I5)")(tmpip(j),j=1,neq-1)
!      close(314)

!      tmpr = sparse_column_permutation(stoich_r_sp,tmpip)
!      tmpp = sparse_column_permutation(stoich_p_sp,tmpip)
!
!      call print_sparse_to_file(tmpr,'tmpr.dat')
!      call print_sparse_to_file(stoich_r_sp,'stoich_r_sp.dat')
!      call print_sparse_to_file(tmpp,'tmpp.dat')
!      call print_sparse_to_file(stoich_p_sp,'stoich_p_sp.dat')
!      stop

!     4) Preprocess second linear system
!      call system_2_prepro(neq,sparseJ,alph,beta,iper2,iiper2,iLU2,rLU2)

   R5_sys2 = system_2_matrix(dble(alph),dble(beta),sparseJ)


end subroutine sparse_radau5_algebra_allocation
