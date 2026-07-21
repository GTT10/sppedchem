!-----------------------------------------------------------------------
!     ADDITIONAL LINEAR ALGEBRA ROUTINES REQUIRED BY GAM
!-----------------------------------------------------------------------
!     VERSION OF AUGUST 20, 1997
!-----------------------------------------------------------------------
!
SUBROUTINE DECLU(R,JF0,H,LDJAC,LU,LDLU,IPIV,ORD,IER,IJOB)

   use working_precision, only: dp
   use sparse_definitions
   use sparse_chemistry, only: JAC_sparse
   use sparse_algebra,   only: identity
   use ode_solver,       only: R5_sys1
   IMPLICIT NONE
!
!   COMMON
!------------------------------------
   COMMON/LINAL/MLLU,MULU,MDIAG
!$ OMP THREADPRIVATE(/LINAL/)
!
!   INPUT VARIABLES
!------------------------------------
   INTEGER R, LDJAC, LDLU, ORD, MLLU, MULU, MDIAG, IJOB
   DOUBLE PRECISION  JF0(LDJAC,1), H
!
!   OUTPUT VARIABLES
!------------------------------------
   INTEGER IER, IPIV(R)
   DOUBLE PRECISION LU(LDLU,1)
!
!   LOCAL VARIABLES
!------------------------------------
   INTEGER I,J
   DOUBLE PRECISION  FAC, L31, L51, L71, L91
   PARAMETER(L31  =  6.411501944628007d-01,&
   &L51  =  6.743555662880509D-01,&
   &L71  =  7.109158294404152D-01,&
   &L91  =  7.440547954061898d-01)
!
!   EXECUTABLE STATEMENTS
!---------------------------------
!
   GOTO (10,20,30,40), ORD
10 FAC = -(L31*H)
   GOTO 50
20 FAC = -(L51*H)
   GOTO 50
30 FAC = -(L71*H)
   GOTO 50
40 FAC = -(L91*H)
50 CONTINUE

   GO TO (1,2) IJOB

1  CONTINUE

! -------- JACOBIAN A FULL MATRIX

!      DO J=1,R
!         DO  I=1,R
!            LU(I,J)= FAC*JF0(I,J)
!         END DO
!         LU(J,J)=LU(J,J)+1d0
!      END DO
!      CALL DEC (R,LDLU,LU,IPIV,IER)

   R5_sys1 = real(FAC,dp)*JAC_sparse + identity(R,real(1.d0-FAC,dp))
   call sparseLU(R5_sys1)

   IER = 0
   RETURN

2  CONTINUE

! -------- JACOBIAN A BAND MATRIX

   DO J=1,R
      DO I=1,MDIAG
         LU(I+MLLU,J)= FAC*JF0(I,J)
      END DO
      LU(MDIAG,J)=LU(MDIAG,J)+1d0
   END DO
   CALL DECB (R,LDLU,LU,MLLU,MULU,IPIV,IER)
   RETURN

END
!
!  SUBROUTINE SOLLU
!
SUBROUTINE SOLLU(R,LU,LDLU,F,IPIV,IJOB)
   IMPLICIT NONE
!
!   COMMON
!------------------------------------
   COMMON/LINAL/MLLU,MULU,MDIAG
!$ OMP THREADPRIVATE(/LINAL/)
!
!   INPUT VARIABLES
!------------------------------------
   INTEGER R, LDLU, IPIV(R), MLLU, MULU, MDIAG, IJOB
   DOUBLE PRECISION  LU(LDLU,1)
!
!   INPUT/OUTPUT VARIABLES
!------------------------------------
   DOUBLE PRECISION F(R)
!
!   EXECUTABLE STATEMENTS
!---------------------------------
!

   GO TO (1,2) IJOB

1  CONTINUE

! -------- JACOBIAN A FULL MATRIX

   CALL SOL (R,LDLU,LU,F,IPIV)

   RETURN

2  CONTINUE


! -------- JACOBIAN A BAND MATRIX

   CALL SOLB (R,LDLU,LU,MLLU,MULU,F,IPIV)
   RETURN

END

!
!  SUBROUTINE NEWTGS
!
SUBROUTINE NEWTGS(R,DBLK,LU,LDLU,IPIV,F,DN,IJOB)
   use ode_solver, only: R5_sys1
   use sparse_definitions
   IMPLICIT NONE
!
!   INPUT VARIABLES
!------------------------------------
   INTEGER R, LDLU, IPIV(R), IJOB, DBLK
   DOUBLE PRECISION  LU(LDLU,1), F(R,DBLK)
!
!   OUTPUT VARIABLES
!------------------------------------
   DOUBLE PRECISION DN(R,DBLK)
!
!   LOCAL VARIABLES
!------------------------------------
   INTEGER I, J
!
!   EXECUTABLE STATEMENTS
!---------------------------------
!
   DO I=1, R
      DN(I,1) = -F(I,1)
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN(1,1),IPIV,IJOB)
   DN(1:R,1) = R5_sys1 .backslash. DN(1:R,1)
   DO J=2,DBLK
      DO I=1, R
         DN(I,J) =  -F(I,J)+DN(I,J-1)
      END DO
      !CALL  SOLLU(R,LU,LDLU,DN(1,J),IPIV,IJOB)
      DN(1:R,J) = R5_sys1 .backslash. DN(1:R,J)
   END DO
   RETURN
END
!
!    SUBROUTINE INTERP
!
SUBROUTINE INTERP(R,TP,YP,T1,F1,NT1,DBLKOLD,DBLK,T0,Y0,ORD)
   IMPLICIT NONE
!
!   INPUT VARIABLES
!------------------------------------
   INTEGER R,  DBLK, DBLKOLD, ORD, NT1
   DOUBLE PRECISION  T0, Y0(R), T1(1), F1(R,1)
!
!   OUTPUT VARIABLES
!------------------------------------
   DOUBLE PRECISION TP(1), YP(R,1)
!
!   LOCAL VARIABLES
!------------------------------------
   INTEGER I, J, N, IT1, NT2
!
!   EXECUTABLE STATEMENTS
!---------------------------------
!
   NT2 = NT1+1
   N = DBLKOLD+1


   DO IT1=2,DBLK+1
      DO I=1,R
         YP(I,IT1) = F1(I,NT1)
         DO J=NT2,N
            YP(I,IT1) = YP(I,IT1)*(T1(IT1-1)-TP(J)) + F1(I,J)
         END DO
      END DO
   END DO

   DO J=1,R
      YP(J,1) = Y0(J)
   END DO

   GOTO (10,20,30,40) ORD
10 CONTINUE
   TP(1) = T0
   TP(2) = T1(1)
   TP(3) = T1(2)
   TP(4) = T1(3)
   TP(5) = T1(4)
   RETURN
20 CONTINUE
   TP(1) = T0
   TP(2) = T1(1)
   TP(3) = T1(2)
   TP(4) = T1(3)
   TP(5) = T1(4)
   TP(6) = T1(5)
   TP(7) = T1(6)
   RETURN
30 CONTINUE
   TP(1) = T0
   TP(2) = T1(1)
   TP(3) = T1(2)
   TP(4) = T1(3)
   TP(5) = T1(4)
   TP(6) = T1(5)
   TP(7) = T1(6)
   TP(8) = T1(7)
   TP(9) = T1(8)
   RETURN
40 CONTINUE
   TP(1) = T0
   TP(2) = T1(1)
   TP(3) = T1(2)
   TP(4) = T1(3)
   TP(5) = T1(4)
   TP(6) = T1(5)
   TP(7) = T1(6)
   TP(8) = T1(7)
   TP(9) = T1(8)
   TP(10) = T1(9)
   RETURN
END
!
!    SUBROUTINE DIFFDIV
!

SUBROUTINE DIFFDIV(TP,YP,R,DBLK,NT1)
   IMPLICIT NONE
!
!   INPUT VARIABLES
!------------------------------------
   INTEGER R,  DBLK
!
!   INPUT/OUTPUT VARIABLES
!------------------------------------
   INTEGER NT1
   DOUBLE PRECISION TP(1), YP(R,1)
!
!   LOCAL VARIABLES
!------------------------------------
   INTEGER I, J, K, N
!
!   EXECUTABLE STATEMENTS
!---------------------------------
!

   N = DBLK+1
   NT1 = 3


   DO J=N-1,NT1,-1
      DO K=1,J
         DO I=1,R
            YP(I,K)= ( YP(I,K)- YP(I,K+1) )/( TP(K)-TP(K+N-J))
         END DO
      END DO
   END DO
   RETURN
END
!
!     FUNCTION  CONTR
!
DOUBLE PRECISION FUNCTION CONTR(I,R,T,TP,FF,DBLK,NT1)
! ----------------------------------------------------------
!     THIS FUNCTION CAN BE USED FOR CONTINUOUS OUTPUT. IT PROVIDES AN
!     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT T.
!     IT GIVES THE VALUE OF THE INTERPOLATION POLYNOMIAL, DEFINED FOR
!     THE LAST SUCCESSFULLY COMPUTED STEP.
! ----------------------------------------------------------
   IMPLICIT NONE
!
!   INPUT VARIABLES
!------------------------------------
   INTEGER R, I, DBLK, NT1
   DOUBLE PRECISION T, TP(1), FF(R,1)
!
!   LOCAL VARIABLES
!------------------------------------
   INTEGER J, N, NT2
   DOUBLE PRECISION YP
!
!   EXECUTABLE STATEMENTS
!---------------------------------
!
   N = DBLK+1
   NT2=NT1+1
   YP = FF(I,NT1)
   DO J=NT2,N
      YP = YP*(T-TP(J)) + FF(I,J)
   END DO
   CONTR = YP
   RETURN
END

!
!  SUBROUTINE TERMNOT3  (ORDER 3)
!
SUBROUTINE  TERMNOT3(R,FCN,H,IT,DN, F1,FP,YP,TP,NFCN,&
&ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV, SCAL,IJOB,TER,&
&RPAR,IPAR)

   use ode_solver, only: R5_sys1
   use sparse_definitions


   IMPLICIT NONE

!
!   INCLUDE
!------------------------------------
   INCLUDE "gamparam.dat"
!
!   INPUT VARIABLES
!------------------------------------
   INTEGER R, IT, IJOB, IPIV(R), LDLU, IPAR(1)
   DOUBLE PRECISION  H, SCAL(R), TP(1), ERRNEWT0,&
   &LU(LDLU,R),RPAR(1)
!
!   INPUT/OUTPUT VARIABLES
!------------------------------------
   INTEGER NFCN
   DOUBLE PRECISION  ERRNEWT, TETAK0, YP(R,1), FP(R,1), F1(R,1),&
   &DN(R)
   LOGICAL TER
!
!   LOCAL VARIABLES
!------------------------------------
   INTEGER  J
   DOUBLE PRECISION  ERRVJ, SUM
!
!   EXTERNAL FUNCTIONS
!------------------------------------

   EXTERNAL FCN

!
!   EXECUTABLE STATEMENTS
!---------------------------------
!
!--------- ONE STEP OF THE ITERATION PROCEDURE
   TER = .FALSE.
!
!
   DO J=1,R
      SUM = B311*FP(J,1)+B312*FP(J,2)+B313*FP(J,3)
      DN(J)=YP(J,2)-YP(J,1)-H*SUM
   END DO
!      CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)

   ERRVJ = 0D0
   DO J=1,R
      YP(J,2)=YP(J,2)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = dsqrt(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF
   CALL FCN(R,TP(2),YP(1,2),F1(1,1), RPAR,IPAR)
   NFCN = NFCN + 1


   DO J=1,R
      SUM = L321*(F1(J,1)-FP(J,2))+B311*FP(J,2)&
      &+B312*FP(J,3)+B313*FP(J,4)
      DN(J)=YP(J,3)-YP(J,2)-H*SUM
   END DO
!      CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,3)=YP(J,3)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF
   CALL FCN(R,TP(3),YP(1,3),F1(1,2), RPAR,IPAR)
   NFCN = NFCN + 1


   DO J=1,R
      SUM = L331*(F1(J,1)-FP(J,2))+L332*(F1(J,2)-FP(J,3))&
      &+B311*FP(J,3)+B312*FP(J,4)+B313*FP(J,5)
      DN(J)=YP(J,4)-YP(J,3)-H*SUM
   END DO
!      CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,4)=YP(J,4)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF
   CALL FCN(R,TP(4),YP(1,4),F1(1,3), RPAR,IPAR)
   NFCN = NFCN + 1

   DO J=1,R
      SUM = L341*(F1(J,1)-FP(J,2))+L342*(F1(J,2)-FP(J,3))&
      &+L343*(F1(J,3)-FP(J,4))&
      &+B313*FP(J,3)+B312*FP(J,4)+B311*FP(J,5)
      DN(J)=YP(J,5)-YP(J,4)-H*SUM
   END DO
!       CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,5)=YP(J,5)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF
   CALL FCN(R,TP(5),YP(1,5),F1(1,4), RPAR,IPAR)
   NFCN = NFCN + 1

   DO J=1,R
      FP(J,2) = F1(J,1)
      FP(J,3) = F1(J,2)
      FP(J,4) = F1(J,3)
      FP(J,5) = F1(J,4)
   END DO


   RETURN
END
!
!  SUBROUTINE TERMNOT5  (ORDER 5)
!
SUBROUTINE  TERMNOT5(R,FCN,H,IT,DN, F1,FP,YP,TP,NFCN,&
&ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV, SCAL,IJOB,TER,&
&RPAR,IPAR)

   use ode_solver, only: R5_sys1
   use sparse_definitions

   IMPLICIT NONE
!
!   INCLUDE
!------------------------------------
   INCLUDE "gamparam.dat"
!
!   INPUT VARIABLES
!------------------------------------
   INTEGER R, IT, IJOB, IPIV(R), LDLU, IPAR(1)
   DOUBLE PRECISION  H, SCAL(R), TP(1), ERRNEWT0,&
   &LU(LDLU,R),RPAR(1)
!
!   INPUT/OUTPUT VARIABLES
!------------------------------------
   INTEGER NFCN
   DOUBLE PRECISION  ERRNEWT, TETAK0, YP(R,1), FP(R,1), F1(R,1),&
   &DN(R)
   LOGICAL TER
!
!   LOCAL VARIABLES
!------------------------------------
   INTEGER  J
   DOUBLE PRECISION  ERRVJ, SUM
!
!   EXTERNAL FUNCTIONS
!------------------------------------

   EXTERNAL FCN

!
!   EXECUTABLE STATEMENTS
!---------------------------------
!
!--------- ONE STEP OF THE ITERATION PROCEDURE

   TER = .FALSE.
   DO J=1,R
      SUM = B511*FP(J,1)+B512*FP(J,2)+B513*FP(J,3)+B514*FP(J,4)&
      &+B515*FP(J,5)
      DN(J)=YP(J,2)-YP(J,1)-H*SUM
   END DO
!              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,2)=YP(J,2)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF
   CALL FCN(R,TP(2),YP(1,2),F1(1,1), RPAR,IPAR)
   NFCN = NFCN + 1

   DO J=1,R
      SUM = L521*(F1(J,1)-FP(J,2))+B521*FP(J,1)+B522*FP(J,2)&
      &+B523*FP(J,3)+B524*FP(J,4)+B525*FP(J,5)
      DN(J)=YP(J,3)-YP(J,2)-H*SUM
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,3)=YP(J,3)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF

   CALL FCN(R,TP(3),YP(1,3),F1(1,2), RPAR,IPAR)
   NFCN = NFCN + 1

   DO J=1,R
      SUM = L531*(F1(J,1)-FP(J,2))+L532*(F1(J,2)-FP(J,3))&
      &+B521*FP(J,2)+B522*FP(J,3)+B523*FP(J,4)&
      &+B524*FP(J,5)+B525*FP(J,6)
      DN(J)=YP(J,4)-YP(J,3)-H*SUM
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,4)=YP(J,4)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF
   CALL FCN(R,TP(4),YP(1,4),F1(1,3), RPAR,IPAR)
   NFCN = NFCN + 1

   DO J=1,R
      SUM = L541*(F1(J,1)-FP(J,2))+L542*(F1(J,2)-FP(J,3))&
      &+L543*(F1(J,3)-FP(J,4))+B521*FP(J,3)+B522*FP(J,4)&
      &+B523*FP(J,5)+B524*FP(J,6)+B525*FP(J,7)
      DN(J)=YP(J,5)-YP(J,4)-H*SUM
   END DO
!              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,5)=YP(J,5)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT.TETAK0)) THEN
      TER = .TRUE.
      RETURN
   END IF

   CALL FCN(R,TP(5),YP(1,5),F1(1,4), RPAR,IPAR)
   NFCN = NFCN + 1

   DO J=1,R
      SUM =  L551*(F1(J,1)-FP(J,2))+L552*(F1(J,2)-FP(J,3))&
      &+L553*(F1(J,3)-FP(J,4))+L554*(F1(J,4)-FP(J,5))&
      &+B525*FP(J,3)+B524*FP(J,4)+B523*FP(J,5)+B522*FP(J,6)+B521*FP(J,7)
      DN(J)=YP(J,6)-YP(J,5)-H*SUM
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,6)=YP(J,6)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0)) THEN
      TER = .TRUE.
      RETURN
   END IF
   CALL FCN(R,TP(6),YP(1,6),F1(1,5), RPAR,IPAR)
   NFCN = NFCN + 1


   DO J=1,R
      SUM =  L561*(F1(J,1)-FP(J,2))+L562*(F1(J,2)-FP(J,3))&
      &+L563*(F1(J,3)-FP(J,4))+L564*(F1(J,4)-FP(J,5))&
      &+L565*(F1(J,5)-FP(J,6))+B515*FP(J,3)&
      &+B514*FP(J,4)+B513*FP(J,5)+B512*FP(J,6)+B511*FP(J,7)
      DN(J)=YP(J,7)-YP(J,6)-H*SUM
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,7)=YP(J,7)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = dsqrt(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   if ((it.ge.1).and.(errnewt/errnewt0 .gt.tetak0)) then
      TER = .TRUE.
      RETURN
   end if
   CALL FCN(R,TP(7),YP(1,7),F1(1,6), RPAR,IPAR)
   NFCN = NFCN + 1



   DO J=1,R
      FP(J,2) = F1(J,1)
      FP(J,3) = F1(J,2)
      FP(J,4) = F1(J,3)
      FP(J,5) = F1(J,4)
      FP(J,6) = F1(J,5)
      FP(J,7) = F1(J,6)
   END DO
   RETURN
END

!
!  SUBROUTINE TERMNOT7  (ORDER 7)
!
SUBROUTINE  TERMNOT7(R,FCN,H,IT,DN, F1,FP,YP,TP,NFCN,&
&ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV, SCAL,IJOB,TER,&
&RPAR,IPAR)

   use ode_solver, only: R5_sys1
   use sparse_definitions
   IMPLICIT NONE
!
!   INCLUDE
!------------------------------------
   INCLUDE "gamparam.dat"
!
!   INPUT VARIABLES
!------------------------------------
   INTEGER R, IT, IJOB, IPIV(R), LDLU, IPAR(1)
   DOUBLE PRECISION  H, SCAL(R), TP(1), ERRNEWT0,&
   &LU(LDLU,R),RPAR(1)
!
!   INPUT/OUTPUT VARIABLES
!------------------------------------
   INTEGER NFCN
   DOUBLE PRECISION  ERRNEWT, TETAK0, YP(R,1), FP(R,10), F1(R,9),&
   &DN(R)
   LOGICAL TER
!
!   LOCAL VARIABLES
!------------------------------------
   INTEGER  J
   DOUBLE PRECISION  ERRVJ, SUM
!
!   EXTERNAL FUNCTIONS
!------------------------------------

   EXTERNAL FCN

!
!   EXECUTABLE STATEMENTS
!---------------------------------
!
!--------- ONE STEP OF THE ITERATION PROCEDURE
   TER = .FALSE.
   DO J=1,R
      SUM= B711*FP(J,1)+B712*FP(J,2)+B713*FP(J,3)&
      &+B714*FP(J,4)+B715*FP(J,5)+B716*FP(J,6)+B717*FP(J,7)
      DN(J) =YP(J,2)-YP(J,1)-H*SUM
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,2)=YP(J,2) - DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF
   CALL FCN(R,TP(2),YP(1,2),F1(1,1), RPAR,IPAR)
   NFCN = NFCN + 1

   DO J=1,R
      SUM = L721*(F1(J,1)-FP(J,2))+B721*FP(J,1)+B722*FP(J,2)+&
      &B723*FP(J,3)+B724*FP(J,4)+B725*FP(J,5)+B726*FP(J,6)+B727*FP(J,7)
      DN(J) = YP(J,3)-YP(J,2)-H*SUM
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,3)=YP(J,3)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF
   CALL FCN(R,TP(3),YP(1,3),F1(1,2), RPAR,IPAR)
   NFCN = NFCN + 1

   DO J=1,R
      SUM = +L731*(F1(J,1)-FP(J,2))+L732*(F1(J,2)-FP(J,3))&
      &+B731*FP(J,1)+B732*FP(J,2)+B733*FP(J,3)&
      &+B734*FP(J,4)+B735*FP(J,5)+B736*FP(J,6)+B737*FP(J,7)
      DN(J) =YP(J,4)-YP(J,3)-H*SUM
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,4)=YP(J,4)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF
   CALL FCN(R,TP(4),YP(1,4),F1(1,3), RPAR,IPAR)
   NFCN = NFCN + 1

   DO J=1,R
      SUM =L741*(F1(J,1)-FP(J,2))+L742*(F1(J,2)-FP(J,3))&
      &+L743*(F1(J,3)-FP(J,4))&
      &+B731*FP(J,2)+B732*FP(J,3)+B733*FP(J,4)&
      &+B734*FP(J,5)+B735*FP(J,6)+B736*FP(J,7)+B737*FP(J,8)
      DN(J) =YP(J,5)-YP(J,4)-H*SUM
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,5)=YP(J,5)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF

   CALL FCN(R,TP(5),YP(1,5),F1(1,4), RPAR,IPAR)
   NFCN = NFCN + 1

   DO J=1,R
      SUM = L751*(F1(J,1)-FP(J,2))+L752*(F1(J,2)-FP(J,3))&
      &+L753*(F1(J,3)-FP(J,4))+L754*(F1(J,4)-FP(J,5))&
      &+B731*FP(J,3)+B732*FP(J,4)+B733*FP(J,5)&
      &+B734*FP(J,6)+B735*FP(J,7)+B736*FP(J,8)+B737*FP(J,9)
      DN(J)=YP(J,6)-YP(J,5)-H*SUM
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,6)=YP(J,6)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   end if
   CALL FCN(R,TP(6),YP(1,6),F1(1,5), RPAR,IPAR)
   NFCN = NFCN + 1

   DO J=1,R
      SUM = L761*(F1(J,1)-FP(J,2))+L762*(F1(J,2)-FP(J,3))&
      &+L763*(F1(J,3)-FP(J,4))+L764*(F1(J,4)-FP(J,5))&
      &+L765*(F1(J,5)-FP(J,6))&
      &+B737*FP(J,3)+B736*FP(J,4)+B735*FP(J,5)&
      &+B734*FP(J,6)+B733*FP(J,7)+B732*FP(J,8)+B731*FP(J,9)
      DN(J) = YP(J,7)-YP(J,6)-H*SUM
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,7)=YP(J,7)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF

   CALL FCN(R,TP(7),YP(1,7),F1(1,6), RPAR,IPAR)
   NFCN = NFCN + 1

   DO J=1,R
      SUM = L771*(F1(J,1)-FP(J,2))+L772*(F1(J,2)-FP(J,3))&
      &+L773*(F1(J,3)-FP(J,4))+L774*(F1(J,4)-FP(J,5))&
      &+L775*(F1(J,5)-FP(J,6))+L776*(F1(J,6)-FP(J,7))&
      &+B727*FP(J,3)+B726*FP(J,4)+B725*FP(J,5)&
      &+B724*FP(J,6)+B723*FP(J,7)+B722*FP(J,8)+B721*FP(J,9)
      DN(J) = YP(J,8)-YP(J,7)-H*SUM
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,8)=YP(J,8)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF
   CALL FCN(R,TP(8),YP(1,8),F1(1,7), RPAR,IPAR)
   NFCN = NFCN + 1


   DO J=1,R
      SUM = L781*(F1(J,1)-FP(J,2))+L782*(F1(J,2)-FP(J,3))&
      &+L783*(F1(J,3)-FP(J,4))+L784*(F1(J,4)-FP(J,5))&
      &+L785*(F1(J,5)-FP(J,6))+L786*(F1(J,6)-FP(J,7))&
      &+L787*(F1(J,7)-FP(J,8))&
      &+B717*FP(J,3)+B716*FP(J,4)+B715*FP(J,5)&
      &+B714*FP(J,6)+B713*FP(J,7)+B712*FP(J,8)+B711*FP(J,9)
      DN(J) = YP(J,9)-YP(J,8)-H*SUM
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,9)=YP(J,9)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF
   CALL FCN(R,TP(9),YP(1,9),F1(1,8), RPAR,IPAR)
   NFCN = NFCN + 1

   DO J=1,R
      FP(J,2) = F1(J,1)
      FP(J,3) = F1(J,2)
      FP(J,4) = F1(J,3)
      FP(J,5) = F1(J,4)
      FP(J,6) = F1(J,5)
      FP(J,7) = F1(J,6)
      FP(J,8) = F1(J,7)
      FP(J,9) = F1(J,8)
   END DO

   RETURN
END
!
!  SUBROUTINE TERMNOT9  (ORDER 9)
!
SUBROUTINE  TERMNOT9(R,FCN,H,IT,DN, F1,FP,YP,TP,NFCN,&
&ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV, SCAL,IJOB,TER,&
&RPAR,IPAR)
   use ode_solver, only: R5_sys1
   use sparse_definitions
   IMPLICIT NONE
!
!   INCLUDE
!------------------------------------
   INCLUDE "gamparam.dat"
!
!   INPUT VARIABLES
!------------------------------------
   INTEGER R, IT, IJOB, IPIV(R), LDLU, IPAR(1)
   DOUBLE PRECISION  H, SCAL(R), TP(1), ERRNEWT0,&
   &LU(LDLU,R),RPAR(1)
!
!   INPUT/OUTPUT VARIABLES
!------------------------------------
   INTEGER NFCN
   DOUBLE PRECISION  ERRNEWT, TETAK0, YP(R,10), FP(R,10), F1(R,9),&
   &DN(R)
   LOGICAL TER
!
!   LOCAL VARIABLES
!------------------------------------
   INTEGER  J
   DOUBLE PRECISION  ERRVJ, SUM
!
!   EXTERNAL FUNCTIONS
!------------------------------------

   EXTERNAL FCN

!
!   EXECUTABLE STATEMENTS
!---------------------------------
!
!--------- ONE STEP OF THE ITERATION PROCEDURE
   TER = .FALSE.

   DO J=1,R
      SUM = B911*FP(J,1)+B912*FP(J,2)+B913*FP(J,3)+B914*FP(J,4)&
      &+B915*FP(J,5)+B916*FP(J,6)+B917*FP(J,7)+B918*FP(J,8)+B919*FP(J,9)
      DN(J) = YP(J,2)-YP(J,1)-H*SUM
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,2)=YP(J,2)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF
   CALL FCN(R,TP(2),YP(1,2),F1(1,1), RPAR,IPAR)
   NFCN = NFCN + 1


   DO J=1,R
      SUM = L921*(F1(J,1)-FP(J,2))+B921*FP(J,1)+B922*FP(J,2)&
      &+B923*FP(J,3)+B924*FP(J,4)+B925*FP(J,5)&
      &+B926*FP(J,6)+B927*FP(J,7)+B928*FP(J,8)+B929*FP(J,9)
      DN(J) = YP(J,3)-YP(J,2)-H*SUM
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,3)=YP(J,3)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF

   CALL FCN(R,TP(3),YP(1,3),F1(1,2), RPAR,IPAR)
   NFCN = NFCN + 1

   DO J=1,R
      SUM = L932*(F1(J,2)-FP(J,3))+B931*FP(J,1)&
      &+B932*FP(J,2)+B933*FP(J,3)+B934*FP(J,4)+B935*FP(J,5)&
      &+B936*FP(J,6)+B937*FP(J,7)+B938*FP(J,8)+B939*FP(J,9)
      DN(J) = YP(J,4)-YP(J,3)-H*SUM
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)DN(1:R) = R5_sys1 .backslash. DN(1:R)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,4)=YP(J,4)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF

   CALL FCN(R,TP(4),YP(1,4),F1(1,3), RPAR,IPAR)
   NFCN = NFCN + 1


   DO J=1,R
      SUM = L943*(F1(J,3)-FP(J,4))+B941*FP(J,1)+B942*FP(J,2)&
      &+B943*FP(J,3)+B944*FP(J,4)+B945*FP(J,5)+B946*FP(J,6)&
      &+B947*FP(J,7)+B948*FP(J,8)+B949*FP(J,9)
      DN(J) =YP(J,5)-YP(J,4)-H*SUM
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,5)=YP(J,5)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF

   CALL FCN(R,TP(5),YP(1,5),F1(1,4), RPAR,IPAR)
   NFCN = NFCN + 1


   DO J=1,R
      SUM = L954*(F1(J,4)-FP(J,5))+B941*FP(J,2)+B942*FP(J,3)&
      &+B943*FP(J,4)+B944*FP(J,5)+B945*FP(J,6)+B946*FP(J,7)&
      &+B947*FP(J,8)+B948*FP(J,9)+B949*FP(J,10)
      DN(J) =YP(J,6)-YP(J,5)-H*SUM
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,6)=YP(J,6)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF

   CALL FCN(R,TP(6),YP(1,6),F1(1,5), RPAR,IPAR)
   NFCN = NFCN + 1

   DO J=1,R
      SUM = L965*(F1(J,5)-FP(J,6))+B949*FP(J,2)+B948*FP(J,3)&
      &+B947*FP(J,4)+B946*FP(J,5)+B945*FP(J,6)+B944*FP(J,7)&
      &+B943*FP(J,8)+B942*FP(J,9)+B941*FP(J,10)
      DN(J) =YP(J,7)-YP(J,6)-H*SUM
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,7)=YP(J,7)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF

   CALL FCN(R,TP(7),YP(1,7),F1(1,6), RPAR,IPAR)
   NFCN = NFCN + 1

   DO J=1,R
      SUM = L976*(F1(J,6)-FP(J,7))+B939*FP(J,2)+B938*FP(J,3)&
      &+B937*FP(J,4)+B936*FP(J,5)+B935*FP(J,6)+B934*FP(J,7)&
      &+B933*FP(J,8)+B932*FP(J,9)+B931*FP(J,10)
      DN(J) =YP(J,8)-YP(J,7)-H*SUM
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,8)=YP(J,8)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF
   CALL FCN(R,TP(8),YP(1,8),F1(1,7), RPAR,IPAR)
   NFCN = NFCN + 1

   DO J=1,R
      SUM = L987*(F1(J,7)-FP(J,8))+B929*FP(J,2)+B928*FP(J,3)&
      &+B927*FP(J,4)+B926*FP(J,5)+B925*FP(J,6)+B924*FP(J,7)&
      &+B923*FP(J,8)+B922*FP(J,9)+B921*FP(J,10)
      DN(J) =YP(J,9)-YP(J,8)-H*SUM
   END DO
   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,9)=YP(J,9)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF
   CALL FCN(R,TP(9),YP(1,9),F1(1,8), RPAR,IPAR)
   NFCN = NFCN + 1

   DO J=1,R
      SUM = L998*(F1(J,8)-FP(J,9))+B919*FP(J,2)+B918*FP(J,3)&
      &+B917*FP(J,4)+B916*FP(J,5)+B915*FP(J,6)+B914*FP(J,7)&
      &+B913*FP(J,8)+B912*FP(J,9)+B911*FP(J,10)
      DN(J) =YP(J,10)-YP(J,9)-H*SUM
   END DO

   !CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
   DN(1:R) = R5_sys1 .backslash. DN(1:R)
   ERRVJ = 0D0
   DO J=1,R
      YP(J,10)=YP(J,10)-DN(J)
      SUM = (DN(J)/SCAL(J))
      ERRVJ =  ERRVJ + SUM*SUM
   END DO
   ERRVJ = DSQRT(ERRVJ/R)
   ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
   IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
      TER = .TRUE.
      RETURN
   END IF
   CALL FCN(R,TP(10),YP(1,10),F1(1,9), RPAR,IPAR)
   NFCN = NFCN + 1

   DO J=1,R
      FP(J,2) = F1(J,1)
      FP(J,3) = F1(J,2)
      FP(J,4) = F1(J,3)
      FP(J,5) = F1(J,4)
      FP(J,6) = F1(J,5)
      FP(J,7) = F1(J,6)
      FP(J,8) = F1(J,7)
      FP(J,9) = F1(J,8)
      FP(J,10) = F1(J,9)
   END DO
   RETURN
END


!
!  SUBROUTINE ESTERR
!  ERRORS ESTIMATION.  ERRSAME: THE CURRENT ORDER
!                         ERRUP: GREATER ORDER (THAN ERRSAME)
!                       ERRDOWN: LOWER ORDER
!
SUBROUTINE  ESTERR(ERRV, ERRSAME, ERRUP, ERRDOWN, FP,&
&R, H, ORD, DBLK, LU, LDLU,&
&IPIV, F, DN,SCAL,ORDMAX,ORDMIN,IJOB)
   use ode_solver, only: R5_sys1
   use sparse_definitions
   IMPLICIT NONE
!
!   INCLUDE
!------------------------------------
   INCLUDE "gamparam.dat"
!
!   INPUT VARIABLES
!------------------------------------
   INTEGER R, ORD,ORDMIN, ORDMAX, DBLK, IJOB, IPIV(R), LDLU
   DOUBLE PRECISION  H, SCAL(R),  LU(LDLU,R)

!
!   INPUT/OUTPUT VARIABLES
!------------------------------------
   DOUBLE PRECISION  ERRV(DBLK), ERRSAME, ERRUP, ERRDOWN,&
   &FP(R,10), F(R,9), DN(R,1)
!
!   LOCAL VARIABLES
!------------------------------------
   INTEGER  I, J, ORDDOWN
   DOUBLE PRECISION ERRVJ,&
   &FP1,FP2,FP3,FP4,FP5,FP6,FP7,FP8,FP9,FP10
!
!   EXECUTABLE STATEMENTS
!---------------------------------
!
!--------- ERRSAME ESTIMATION
   GOTO (10,20,30,40),  ORD
10 CONTINUE
   DO I=1,R
      FP1= FP(I,1)
      FP2= FP(I,2)
      FP3= FP(I,3)
      FP4= FP(I,4)
      FP5= FP(I,5)
      F(I,1) = H*(B3511*FP1+B3512*FP2+B3513*FP3+B3514*FP4+B3515*FP5)
      F(I,2) = H*(B3521*FP1+B3522*FP2+B3523*FP3+B3524*FP4+B3525*FP5)
      F(I,3) = H*(B3531*FP1+B3532*FP2+B3533*FP3+B3534*FP4+B3535*FP5)
      F(I,4) = H*(B3541*FP1+B3542*FP2+B3543*FP3+B3544*FP4+B3545*FP5)
   END DO

   GOTO 50
20 CONTINUE
   DO J=1,R
      FP1= FP(J,1)
      FP2= FP(J,2)
      FP3= FP(J,3)
      FP4= FP(J,4)
      FP5= FP(J,5)
      FP6= FP(J,6)
      FP7= FP(J,7)
      F(J,1) = H*(B5711*FP1+B5712*FP2+B5713*FP3&
      &+B5714*FP4+B5715*FP5+B5716*FP6+B5717*FP7)
      F(J,2) = H*(B5721*FP1+B5722*FP2+B5723*FP3&
      &+B5724*FP4+B5725*FP5+B5726*FP6+B5727*FP7)
      F(J,3) = H*(B5731*FP1+B5732*FP2+B5733*FP3&
      &+B5734*FP4+B5735*FP5+B5736*FP6+B5737*FP7)
      F(J,4) = H*(B5741*FP1+B5742*FP2+B5743*FP3&
      &+B5744*FP4+B5745*FP5+B5746*FP6+B5747*FP7)
      F(J,5) = H*(B5727*FP1+B5726*FP2+B5725*FP3&
      &+B5724*FP4+B5723*FP5+B5722*FP6+B5721*FP7)
      F(J,6) = H*(B5717*FP1+B5716*FP2+B5715*FP3&
      &+B5714*FP4+B5713*FP5+B5712*FP6+B5711*FP7)
   END DO

   GOTO 50
30 CONTINUE
   DO J=1,R
      FP1= FP(J,1)
      FP2= FP(J,2)
      FP3= FP(J,3)
      FP4= FP(J,4)
      FP5= FP(J,5)
      FP6= FP(J,6)
      FP7= FP(J,7)
      FP8= FP(J,8)
      FP9= FP(J,9)
      F(J,1) = H*(B7911*FP1+B7912*FP2+B7913*FP3+B7914*FP4&
      &+B7915*FP5+B7916*FP6+B7917*FP7+B7918*FP8+B7919*FP9)
      F(J,2) = H*(B7921*FP1+B7922*FP2+B7923*FP3+B7924*FP4&
      &+B7925*FP5+B7926*FP6+B7927*FP7+B7928*FP8+B7929*FP9)
      F(J,3) = H*(B7931*FP1+B7932*FP2+B7933*FP3+B7934*FP4&
      &+B7935*FP5+B7936*FP6+B7937*FP7+B7938*FP8+B7939*FP9)
      F(J,4) = H*(B7941*FP1+B7942*FP2+B7943*FP3+B7944*FP4&
      &+B7945*FP5+B7946*FP6+B7947*FP7+B7948*FP8+B7949*FP9)
      F(J,5) = H*(B7951*FP1+B7952*FP2+B7953*FP3+B7954*FP4&
      &+B7955*FP5+B7956*FP6+B7957*FP7+B7958*FP8+B7959*FP9)
      F(J,6) = H*(B7939*FP1+B7938*FP2+B7937*FP3+B7936*FP4&
      &+B7935*FP5+B7934*FP6+B7933*FP7+B7932*FP8+B7931*FP9)
      F(J,7) = H*(B7929*FP1+B7928*FP2+B7927*FP3+B7926*FP4&
      &+B7925*FP5+B7924*FP6+B7923*FP7+B7922*FP8+B7921*FP9)
      F(J,8) = H*(B7919*FP1+B7918*FP2+B7917*FP3+B7916*FP4&
      &+B7915*FP5+B7914*FP6+B7913*FP7+B7912*FP8+B7911*FP9)
   END DO

   GOTO 50
40 CONTINUE
   DO J=1,R
      FP1= FP(J,1)
      FP2= FP(J,2)
      FP3= FP(J,3)
      FP4= FP(J,4)
      FP5= FP(J,5)
      FP6= FP(J,6)
      FP7= FP(J,7)
      FP8= FP(J,8)
      FP9= FP(J,9)
      FP10= FP(J,10)
      F(J,1) = H*(B91011*(FP1-FP10)+B91012*(FP2-FP9)+B91013*(FP3-FP8)&
      &+B91014*(FP4-FP7)+B91015*(FP5-FP6) )
      F(J,2) = H*(B91021*(FP1-FP10)+B91022*(FP2-FP9)+B91023*(FP3-FP8)&
      &+B91024*(FP4-FP7)+B91025*(FP5-FP6) )
      F(J,3) = H*(B91031*(FP1-FP10)+B91032*(FP2-FP9)+B91033*(FP3-FP8)&
      &+B91034*(FP4-FP7)+B91035*(FP5-FP6) )
      F(J,4) = H*(B91041*(FP1-FP10)+B91042*(FP2-FP9)+B91043*(FP3-FP8)&
      &+B91044*(FP4-FP7)+B91045*(FP5-FP6) )
      F(J,5) =  F(J,4)
      F(J,6) = -F(J,4)
      F(J,7) = -F(J,3)
      F(J,8) = -F(J,2)
      F(J,9) = -F(J,1)
   END DO

50 CONTINUE

!--------- A SINGLE SPLITTING-NEWTON ITERATION FOR ERRSAME
   CALL NEWTGS(R,DBLK,LU,LDLU,IPIV,F,DN,IJOB)


!--------- COMPUTE  ERRSAME AND ERRV (VECTOR ERROR)
   ERRSAME = 0D0
   DO J=1,DBLK
      ERRV(J) = 0D0
      DO I=1,R
         FP1 = (DN(I,J)/SCAL(I) )
         ERRV(J) =  ERRV(J)+ FP1*FP1
      END DO
      ERRV(J) = DSQRT(ERRV(J)/R)
      ERRSAME = DMAX1( ERRSAME, ERRV(J) )
   END DO
   ERRSAME = DMAX1(ERRSAME, 1d-15)
   ERRDOWN = 0D0
   ERRUP   = 0D0
   IF ( ERRSAME .LE. 1d0) THEN

      IF (ORD .LT. ORDMAX) THEN
!--------- ERRUP ESTIMATION
         GOTO (11,21,31), ORD
11       CONTINUE
         DO I=1,R
            FP1 =  F(I,1)/CP31
            FP2 =  F(I,2)/CP31
            FP3 =  F(I,3)/CP31
            FP4 = -F(I,4)/CP31
            F(I,1) =  (-FP1 + 2d0*FP2 - FP3)*CP51
            F(I,2) =  (-FP1 + 2d0*FP2 - FP3)*CP52
            F(I,3) = -(-FP2 + 2d0*FP3 - FP4)*CP52
            F(I,4) = -(-FP2 + 2d0*FP3 - FP4)*CP51
         END DO

         GOTO 41
21       CONTINUE
         DO I = 1, R
            FP1 =  F(I,1)/CP51
            FP2 =  F(I,2)/CP52
            FP3 =  F(I,3)/CP52
            FP4 =  F(I,4)/CP52
            FP5 = -F(I,5)/CP52
            FP6 = -F(I,6)/CP51
            F(I,1) =  (-FP1 + 2d0*FP2 - FP3)*CP71
            F(I,2) =  (-FP1 + 2d0*FP2 - FP3)*CP72
            F(I,3) =  (-FP2 + 2d0*FP3 - FP4)*CP73
            F(I,4) = -(-FP3 + 2d0*FP4 - FP5)*CP73
            F(I,5) = -(-FP4 + 2d0*FP5 - FP6)*CP72
            F(I,6) = -(-FP4 + 2d0*FP5 - FP6)*CP71

         END DO
         GOTO 41
31       CONTINUE
         DO I = 1, R
            FP1 =  F(I,1)/CP71
            FP2 =  F(I,2)/CP72
            FP3 =  F(I,3)/CP73
            FP4 =  F(I,4)/CP73
            FP5 =  F(I,5)/CP73
            FP6 = -F(I,6)/CP73
            FP7 = -F(I,7)/CP72
            FP8 = -F(I,8)/CP71
            F(I,1) =  (-FP1 + 2d0*FP2 - FP3)*CP91
            F(I,2) =  (-FP1 + 2d0*FP2 - FP3)*CP92
            F(I,3) =  (-FP2 + 2d0*FP3 - FP4)*CP93
            F(I,4) =  (-FP3 + 2d0*FP4 - FP5)*CP94
            F(I,5) = -(-FP4 + 2d0*FP5 - FP6)*CP94
            F(I,6) = -(-FP5 + 2d0*FP6 - FP7)*CP93
            F(I,7) = -(-FP6 + 2d0*FP7 - FP8)*CP92
            F(I,8) = -(-FP6 + 2d0*FP7 - FP8)*CP91
         END DO

41       CONTINUE

!--------- A SINGLE SPLITTING-NEWTON ITERATION FOR ERRUP

         CALL NEWTGS(R,DBLK,LU,LDLU,IPIV,F,DN,IJOB)


!-------- COMPUTE ERRUP
         ERRUP = 0D0
         DO J=1,DBLK
            ERRVJ = 0D0
            DO I=1,R
               FP1 = (DN(I,J)/SCAL(I) )
               ERRVJ =  ERRVJ + FP1*FP1
            END DO
            ERRVJ = DSQRT(ERRVJ/R)
            ERRUP = DMAX1( ERRUP, ERRVJ )
         END DO
         ERRUP = dmax1(ERRUP, 1d-15)
      END IF
      IF (ORD .GT. ORDMIN) THEN
!--------- ERRDOWN ESTIMATION
         ORDDOWN = ORD-1
         GOTO (13, 23, 33), ORDDOWN
13       CONTINUE
         DO J=1,R
            FP1= FP(J,1)
            FP2= FP(J,2)
            FP3= FP(J,3)
            FP4= FP(J,4)
            FP5= FP(J,5)
            FP6= FP(J,6)
            FP7= FP(J,7)
            F(J,1) = -H*(B3511*FP1+B3512*FP2+B3513*FP3+B3514*FP4+B3515*FP5)
            F(J,2) = -H*(B3521*FP1+B3522*FP2+B3523*FP3+B3524*FP4+B3525*FP5)
            F(J,3) = -H*(B3521*FP2+B3522*FP3+B3523*FP4+B3524*FP5+B3525*FP6)
            F(J,4) = -H*(B3521*FP3+B3522*FP4+B3523*FP5+B3524*FP6+B3525*FP7)
            F(J,5) = -H*(B3531*FP3+B3532*FP4+B3533*FP5+B3534*FP6+B3535*FP7)
            F(J,6) = -H*(B3541*FP3+B3542*FP4+B3543*FP5+B3544*FP6+B3545*FP7)
         END DO

         GOTO 43
23       CONTINUE
         DO J=1,R
            FP1= FP(J,1)
            FP2= FP(J,2)
            FP3= FP(J,3)
            FP4= FP(J,4)
            FP5= FP(J,5)
            FP6= FP(J,6)
            FP7= FP(J,7)
            FP8= FP(J,8)
            FP9= FP(J,9)
            F(J,1) = -H*(B5711*FP1+B5712*FP2+B5713*FP3&
            &+B5714*FP4+B5715*FP5+B5716*FP6+B5717*FP7)
            F(J,2) = -H*(B5721*FP1+B5722*FP2+B5723*FP3&
            &+B5724*FP4+B5725*FP5+B5726*FP6+B5727*FP7)
            F(J,3) = -H*(B5731*FP1+B5732*FP2+B5733*FP3&
            &+B5734*FP4+B5735*FP5+B5736*FP6+B5737*FP7)
            F(J,4) = -H*(B5731*FP2+B5732*FP3+B5733*FP4&
            &+B5734*FP5+B5735*FP6+B5736*FP7+B5737*FP8)
            F(J,5) = -H*(B5731*FP3+B5732*FP4+B5733*FP5&
            &+B5734*FP6+B5735*FP7+B5736*FP8+B5737*FP9)
            F(J,6) = -H*(B5741*FP3+B5742*FP4+B5743*FP5&
            &+B5744*FP6+B5745*FP7+B5746*FP8+B5747*FP9)
            F(J,7) = -H*(B5727*FP3+B5726*FP4+B5725*FP5&
            &+B5724*FP6+B5723*FP7+B5722*FP8+B5721*FP9)
            F(J,8) = -H*(B5717*FP3+B5716*FP4+B5715*FP5&
            &+B5714*FP6+B5713*FP7+B5712*FP8+B5711*FP9)

         END DO

         GOTO 43
33       CONTINUE
         DO J=1,R
            FP1= FP(J,1)
            FP2= FP(J,2)
            FP3= FP(J,3)
            FP4= FP(J,4)
            FP5= FP(J,5)
            FP6= FP(J,6)
            FP7= FP(J,7)
            FP8= FP(J,8)
            FP9= FP(J,9)
            FP10= FP(J,10)
            F(J,1) = H*(B7911*FP1+B7912*FP2+B7913*FP3+B7914*FP4&
            &+B7915*FP5+B7916*FP6+B7917*FP7+B7918*FP8+B7919*FP9)

            F(J,2) = H*(B7921*FP1+B7922*FP2+B7923*FP3+B7924*FP4&
            &+B7925*FP5+B7926*FP6+B7927*FP7+B7928*FP8+B7929*FP9)

            F(J,3) = H*(B7931*FP1+B7932*FP2+B7933*FP3+B7934*FP4&
            &+B7935*FP5+B7936*FP6+B7937*FP7+B7938*FP8+B7939*FP9)

            F(J,4) = H*(B7941*FP1+B7942*FP2+B7943*FP3+B7944*FP4&
            &+B7945*FP5+B7946*FP6+B7947*FP7+B7948*FP8+B7949*FP9)

            F(J,5) = -H*(B7941*FP2+B7942*FP3+B7943*FP4+B7944*FP5&
            &+B7945*FP6+B7946*FP7+B7947*FP8+B7948*FP9+B7949*FP10)


            F(J,6) = -H*(B7951*FP2+B7952*FP3+B7953*FP4+B7954*FP5&
            &+B7955*FP6+B7956*FP7+B7957*FP8+B7958*FP9+B7959*FP10)

            F(J,7) = -H*(B7939*FP2+B7938*FP3+B7937*FP4+B7936*FP5&
            &+B7935*FP6+B7934*FP7+B7933*FP8+B7932*FP9+B7931*FP10)

            F(J,8) = -H*(B7929*FP2+B7928*FP3+B7927*FP4+B7926*FP5&
            &+B7925*FP6+B7924*FP7+B7923*FP8+B7922*FP9+B7921*FP10)

            F(J,9) = -H*(B7919*FP2+B7918*FP3+B7917*FP4+B7916*FP5&
            &+B7915*FP6+B7914*FP7+B7913*FP8+B7912*FP9+B7911*FP10)

         END DO

43       CONTINUE

!--------- A SINGLE SPLITTING-NEWTON ITERATION FOR ERRDOWN
         CALL NEWTGS(R,DBLK,LU,LDLU,IPIV,F,DN,IJOB)


!--------- COMPUTE ERRDOWN
         ERRDOWN = 0D0
         DO J=1,DBLK
            ERRVJ = 0D0
            DO I=1,R
               FP1 = (DN(I,J)/SCAL(I) )
               ERRVJ =  ERRVJ + FP1*FP1
            END DO
            ERRVJ = DSQRT(ERRVJ/R)
            ERRDOWN = DMAX1( ERRDOWN, ERRVJ )
         END DO
         ERRDOWN = dmax1(ERRDOWN, 1d-15)
      END IF
   END IF
   RETURN
END
!
!C     SUBROUTINE DEC
!C
!      SUBROUTINE DEC (N, NDIM, A, IP, IER)
!C VERSION REAL DOUBLE PRECISION
!      INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J
!      DOUBLE PRECISION A,T
!      DIMENSION A(NDIM,N), IP(N)
!C-----------------------------------------------------------------------
!C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION.
!C  INPUT..
!C     N = ORDER OF MATRIX.
!C     NDIM = DECLARED DIMENSION OF ARRAY  A .
!C     A = MATRIX TO BE TRIANGULARIZED.
!C  OUTPUT..
!C     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U .
!C     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
!C     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
!C     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
!C     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
!C           SINGULAR AT STAGE K.
!C  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
!C  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N).
!C  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
!C
!C  REFERENCE..
!C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
!C     C.A.C.M. 15 (1972), P. 274.
!C-----------------------------------------------------------------------
!      IER = 0
!      IP(N) = 1
!      IF (N .EQ. 1) GO TO 70
!      NM1 = N - 1
!      DO 60 K = 1,NM1
!        KP1 = K + 1
!        M = K
!        DO 10 I = KP1,N
!          IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I
! 10     CONTINUE
!        IP(K) = M
!        T = A(M,K)
!        IF (M .EQ. K) GO TO 20
!        IP(N) = -IP(N)
!        A(M,K) = A(K,K)
!        A(K,K) = T
! 20     CONTINUE
!        IF (T .EQ. 0.D0) GO TO 80
!        T = 1.D0/T
!        DO 30 I = KP1,N
! 30       A(I,K) = -A(I,K)*T
!        DO 50 J = KP1,N
!          T = A(M,J)
!          A(M,J) = A(K,J)
!          A(K,J) = T
!          IF (T .EQ. 0.D0) GO TO 45
!          DO 40 I = KP1,N
! 40         A(I,J) = A(I,J) + A(I,K)*T
! 45       CONTINUE
! 50       CONTINUE
! 60     CONTINUE
! 70   K = N
!      IF (A(N,N) .EQ. 0.D0) GO TO 80
!      RETURN
! 80   IER = K
!      IP(N) = 0
!      RETURN
!C----------------------- END OF SUBROUTINE DEC -------------------------
!      END
!
!     SUBROUTINE SOL
!
!      SUBROUTINE SOL (N, NDIM, A, B, IP)
!C VERSION REAL DOUBLE PRECISION
!      INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1
!      DOUBLE PRECISION A,B,T
!      DIMENSION A(NDIM,N), B(N), IP(N)
!C-----------------------------------------------------------------------
!C  SOLUTION OF LINEAR SYSTEM, A*X = B .
!C  INPUT..
!C    N = ORDER OF MATRIX.
!C    NDIM = DECLARED DIMENSION OF ARRAY  A .
!C    A = TRIANGULARIZED MATRIX OBTAINED FROM DEC.
!C    B = RIGHT HAND SIDE VECTOR.
!C    IP = PIVOT VECTOR OBTAINED FROM DEC.
!C  DO NOT USE IF DEC HAS SET IER .NE. 0.
!C  OUTPUT..
!C    B = SOLUTION VECTOR, X .
!C-----------------------------------------------------------------------
!      IF (N .EQ. 1) GO TO 50
!      NM1 = N - 1
!      DO 20 K = 1,NM1
!        KP1 = K + 1
!        M = IP(K)
!        T = B(M)
!        B(M) = B(K)
!        B(K) = T
!        DO 10 I = KP1,N
! 10       B(I) = B(I) + A(I,K)*T
! 20     CONTINUE
!      DO 40 KB = 1,NM1
!        KM1 = N - KB
!        K = KM1 + 1
!        B(K) = B(K)/A(K,K)
!        T = -B(K)
!        DO 30 I = 1,KM1
! 30       B(I) = B(I) + A(I,K)*T
! 40     CONTINUE
! 50   B(1) = B(1)/A(1,1)
!      RETURN
!C----------------------- END OF SUBROUTINE SOL -------------------------
!      END
!
!     SUBROUTINE DECB
!
!      SUBROUTINE DECB (N, NDIM, A, ML, MU, IP, IER)
!      REAL*8 A,T
!      DIMENSION A(NDIM,N), IP(N)
!C-----------------------------------------------------------------------
!C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED
!C  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU
!C  INPUT..
!C     N       ORDER OF THE ORIGINAL MATRIX A.
!C     NDIM    DECLARED DIMENSION OF ARRAY  A.
!C     A       CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS
!C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  A  AND
!C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
!C                ML+1 THROUGH 2*ML+MU+1 OF  A.
!C     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
!C     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
!C  OUTPUT..
!C     A       AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
!C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
!C     IP      INDEX VECTOR OF PIVOT INDICES.
!C     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O .
!C     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE
!C                SINGULAR AT STAGE K.
!C  USE  SOLB  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
!C  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1.
!C  IF IP(N)=O, A IS SINGULAR, SOLB WILL DIVIDE BY ZERO.
!C
!C  REFERENCE..
!C     THIS IS A MODIFICATION OF
!C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
!C     C.A.C.M. 15 (1972), P. 274.
!C-----------------------------------------------------------------------
!      IER = 0
!      IP(N) = 1
!      MD = ML + MU + 1
!      MD1 = MD + 1
!      JU = 0
!      IF (ML .EQ. 0) GO TO 70
!      IF (N .EQ. 1) GO TO 70
!      IF (N .LT. MU+2) GO TO 7
!      DO 5 J = MU+2,N
!      DO 5 I = 1,ML
!  5   A(I,J) = 0.D0
!  7   NM1 = N - 1
!      DO 60 K = 1,NM1
!        KP1 = K + 1
!        M = MD
!        MDL = MIN(ML,N-K) + MD
!        DO 10 I = MD1,MDL
!          IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I
! 10     CONTINUE
!        IP(K) = M + K - MD
!        T = A(M,K)
!        IF (M .EQ. MD) GO TO 20
!        IP(N) = -IP(N)
!        A(M,K) = A(MD,K)
!        A(MD,K) = T
! 20     CONTINUE
!        IF (T .EQ. 0.D0) GO TO 80
!        T = 1.D0/T
!        DO 30 I = MD1,MDL
! 30       A(I,K) = -A(I,K)*T
!        JU = MIN0(MAX0(JU,MU+IP(K)),N)
!        MM = MD
!        IF (JU .LT. KP1) GO TO 55
!        DO 50 J = KP1,JU
!          M = M - 1
!          MM = MM - 1
!          T = A(M,J)
!          IF (M .EQ. MM) GO TO 35
!          A(M,J) = A(MM,J)
!          A(MM,J) = T
! 35       CONTINUE
!          IF (T .EQ. 0.D0) GO TO 45
!          JK = J - K
!          DO 40 I = MD1,MDL
!            IJK = I - JK
! 40         A(IJK,J) = A(IJK,J) + A(I,K)*T
! 45       CONTINUE
! 50       CONTINUE
! 55     CONTINUE
! 60     CONTINUE
! 70   K = N
!      IF (A(MD,N) .EQ. 0.D0) GO TO 80
!      RETURN
! 80   IER = K
!      IP(N) = 0
!      RETURN
!C----------------------- END OF SUBROUTINE DECB ------------------------
!      END
!
!     SUBROUTINE SOLB
!
!      SUBROUTINE SOLB (N, NDIM, A, ML, MU, B, IP)
!      REAL*8 A,B,T
!      DIMENSION A(NDIM,N), B(N), IP(N)
!C-----------------------------------------------------------------------
!C  SOLUTION OF LINEAR SYSTEM, A*X = B .
!C  INPUT..
!C    N      ORDER OF MATRIX A.
!C    NDIM   DECLARED DIMENSION OF ARRAY  A .
!C    A      TRIANGULARIZED MATRIX OBTAINED FROM DECB.
!C    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
!C    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
!C    B      RIGHT HAND SIDE VECTOR.
!C    IP     PIVOT VECTOR OBTAINED FROM DECB.
!C  DO NOT USE IF DECB HAS SET IER .NE. 0.
!C  OUTPUT..
!C    B      SOLUTION VECTOR, X .
!C-----------------------------------------------------------------------
!      MD = ML + MU + 1
!      MD1 = MD + 1
!      MDM = MD - 1
!      NM1 = N - 1
!      IF (ML .EQ. 0) GO TO 25
!      IF (N .EQ. 1) GO TO 50
!      DO 20 K = 1,NM1
!        M = IP(K)
!        T = B(M)
!        B(M) = B(K)
!        B(K) = T
!        MDL = MIN(ML,N-K) + MD
!        DO 10 I = MD1,MDL
!          IMD = I + K - MD
! 10       B(IMD) = B(IMD) + A(I,K)*T
! 20     CONTINUE
! 25   CONTINUE
!      DO 40 KB = 1,NM1
!        K = N + 1 - KB
!        B(K) = B(K)/A(MD,K)
!        T = -B(K)
!        KMD = MD - K
!        LM = MAX0(1,KMD+1)
!        DO 30 I = LM,MDM
!          IMD = I - KMD
! 30       B(IMD) = B(IMD) + A(I,K)*T
! 40     CONTINUE
! 50   B(1) = B(1)/A(MD,1)
!      RETURN
!C----------------------- END OF SUBROUTINE SOLB ------------------------
!      END
