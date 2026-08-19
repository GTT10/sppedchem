!     *****************************************************************
!     **                                                             **
!     **                     KIVA4 - CHEMISTRY                       **
!     **                                                             **
!     **   ChemkinII interpreter and runtime in FORTRAN 2003 format  **
!     **           Taken from version 4.2, 14/9/1993                 **
!     **                                                             **
!     **                                                             **
!     **   Modified by: Federico Perini                              **
!     **   Last update: wedesday, 30/11/2011                         **
!     **                                                             **
!     *****************************************************************

      module chemkinII

      use chemistry_string_limits, only: species_name_len

      implicit none
      public

!     *****************************************************************
!     **     COMMON BLOCKS DECLARATION                               **
!     *****************************************************************

!     ** CKSTRT block *************************************************
!     NMM    - Total number of elements in problem.
!     NKK    - Total number of species in problem.
!     NII    - Total number of reactions in problem.
!     MXSP   - Maximum number of species (reactants plus products)
!              allowed for any reaction;  unless changed in the
!              interpreter, MXSP=6.
!     MXTB   - Maximum number of enhanced third-bodies allowed fo any
!              reaction;  unless changed in the interpreter, MXTB=10.
!     MXTP   - Maximum number of temperatures allowed in fits of
!              thermodynamic properties for any species;  unless
!              changed in the interpreter and the thermodynamic
!              database, MXTP=3.
!     NCP    - Number of polynomial coefficients to fits of CP/R for
!              a species;  unless changed in the interpreter and the
!              thermodynamic database, NCP=5.
!     NCP1   - NCP + 1
!     NCP2   - NCP + 2
!     NCP2T  - Total number of thermodynamic fit coefficients for the
!              species;  unless changed, NCP2T = (MXTP-1)*NCP2 = 14.
!     NPAR   - Number of parameters required in the rate expression
!              for the reactions;  in the current formulation NPAR=3.
!     NLAR   - Number of parameters required for Landau-Teller
!              reactions; NLAR=4.
!     NFAR   - Number of parameters allowed for fall-off reactions;
!              NFAR=8.
!     NLAN   - Total number of Landau-Teller reactions.
!     NFAL   - Total number of fall-off reactions.
!     NREV   - Total number of reactions with reverse parameters.
!     NTHB   - Total number of reactions with third-bodies.
!     NRLT   - Total number of Landau-Teller reactions with reverse
!              parameters.
!     NWL    - Total number of reactions with radiation wavelength
!              enhancement factors.
!
!  STARTING ADDRESSES FOR THE CHARACTER WORK SPACE, CCKWRK.
!
!     IcMM   - Starting address of an array of the NMM element names.
!              CCKWRK(IcMM+M-1) is the name of the Mth element.
!     IcKK   - Starting address of an array of the NKK species names.
!              CCKWRK(icKK+M-1) is the name of the Kth species.
!
!  STARTING ADDRESSES FOR THE INTEGER WORK SPACE, ICKWRK.
!
!     IcNC  - Starting address of an array of the elemental content
!             of the NMM elements in the NKK species.
!             ICKWRK(IcNC+(K-1)*NMM+M-1) is the number of atoms of the
!             Mth element in the Kth species.
!     IcPH  - Starting address of an array of phases of the NKK species.
!             ICKWRK(IcPH+K-1) = -1, the Kth species is a solid
!                              =  0, the Kth species is a gas
!                              = +1, the Kth species is a liquid
!     IcCH  - Starting address of an array of the electronic charges of
!             the NKK species.
!             ICKWRK(IcCH+K-1) = -2, the Kth species has two excess
!                                    electrons.
!     IcNT  - Starting address of an array of the number of temperatures
!             used to fit thermodynamic coefficients for the
!             NKK species.
!             ICKWRK(IcNT+K-1) = N, N temperatures were used in the fit
!                                   for the Kth species.
!     IcNU  - Starting address of a matrix of stoichiometric
!             coefficients of the MXSP species in the NII reactions.
!             ICKWRK(IcNU+(I-1)*MXSP+N-1) is the coefficient of the Nth
!             participant species in the Ith reaction
!     IcNK  - Starting address of a matrix of species index numbers for
!             the MXSP species in the NII reactions.
!             ICKWRK(IcNK+(I-1)*MXSP+N-1) = K, the species number of
!             the Nth participant species in the Ith reaction.
!     IcNS  - Starting address of an array of the total number of
!             participant species for the NII reactions, and the
!             reversibility of the reactions.
!             ICKWRK(IcNS+I-1) = +N, the Ith reaction is reversible
!                                    and has N participant species
!                                    (reactants + products)
!                              = -N, the Ith reaction is irreversible
!                                    and has N participant species
!                                    (reactants + products)
!     IcNR  - Starting address of an array of the number of reactants
!             only for the NII reactions.
!             ICKWRK(IcNR+I-1) is the total number of reactants in the
!             Ith reaction.
!     IcLT  - Starting address of an array of the NLAN reaction numbers
!             for which Landau-Teller parameters have been given.
!             ICKWRK(IcLT+N-1) is the reaction number of the Nth
!             Landau-Teller reaction.
!     IcRL  - Starting address of an array of the NRLT reaction numbers
!             for which reverse Landau-Teller parameters have been
!             given.
!             ICKWRK(IcRL+N-1) is the reaction number of the Nth
!             reaction with reverse Landau-Teller parameters.
!     IcRV  - Starting address of an array of the NREV reaction numbers
!             for which reverse Arhennius coefficients have been given.
!             ICKWRK(IcRV+N-1) is the reaction number of the Nth
!             reaction with reverse coefficients.
!     IcWL  - Starting address of an array of the NWL reactions numbers
!             for which radiation wavelength has been given.
!             ICKWRK(IcWL+N-1) is the reaction number of the Nth
!             reaction with wavelength enhancement.
!     IcFL  - Starting address of an array of the NFAL reaction numbers
!             with fall-off parameters.
!             ICKWRK(IcFL+N-1) is the reaction number of the Nth
!             fall-off reaction.
!     IcFO  - Starting address of an array describing the type of
!             the NFAL fall-off reactions.
!             ICKWRK(IcFO+N-1) is the type of the Nth fall-off
!             reaction: 1 for 3-parameter Lindemann Form
!                       2 for 6- or 8-parameter SRI Form
!                       3 for 6-parameter Troe Form
!                       4 for 7-parameter Troe form
!     IcKF  - Starting address of an array of the third-body species
!             numbers for the NFAL fall-off reactions.
!             ICKWRK(IcKF+N-1) = 0: the concentration of the third-body
!                                   is the total of the concentrations
!                                   of all species in the problem
!                              = K: the concentration of the third-body
!                                   is the concentration of species K.
!     IcTB  - Starting address of an array of reaction numbers for the
!             NTHB third-body reactions.
!             ICKWRK(IcTB+N-1) is the reaction number of the Nth
!             third-body reaction.
!     IcKN  - Starting address of an array of the number of enhanced
!             third bodies for the NTHB third-body reactions.
!             ICKWRK(IcKN+N-1) is the number of enhanced species for
!             the Nth third-body reaction.
!     IcKT  - Starting address of an array of species numbers for the
!             MXTB enhanced 3rd bodies in the NTHB third-body reactions.
!             ICKWRK(IcTB+(N-1)*MXTB+L-1) is the species number of the
!             Lth enhanced species in the Nth third-body reaction.
!
!  STARTING ADDRESSES FOR THE REAL WORK SPACE, RCKWRK.
!
!     NcAW  - Starting address of an array of atomic weights of the
!             NMM elements (gm/mole).
!             RCKWRK(NcAW+M-1) is the atomic weight of element M.
!     NcWT  - Starting address of an array of molecular weights for
!             the NKK species (gm/mole).
!             RCKWRK(NcWT+K-1) is the molecular weight of species K.
!     NcTT  - Starting address of an array of MXTP temperatures used in
!             the fits of thermodynamic properties of the NKK species
!             (Kelvins).
!             RCKWRK(NcTT+(K-1)*MXTP+N-1) is the Nth temperature for the
!             Kth species.
!     NcAA  - Starting address of a three-dimensional array of
!             coefficients for the NCP2 fits to the thermodynamic
!             properties for the NKK species, for (MXTP-1) temperature
!             ranges.
!             RCKWRK(NcAA+(L-1)*NCP2+(K-1)*NCP2T+N-1) = A(N,L,K);
!             A(N,L,K),N=1,NCP2T = polynomial coefficients in the fits
!             for the Kth species and the Lth temperature range, where
!             the total number of temperature ranges for the Kth species
!             is ICKWRK(IcNT+K-1) - 1.
!     NcCO  - Starting address of an array of NPAR Arrhenius parameters
!             for the NII reactions.
!             RCKWRK(NcCO+(I-1)*NPAR+(L-1)) is the Lth parameter of the
!             Ith reaction, where
!                L=1 is the pre-exponential factor (mole-cm-sec-K),
!                L=2 is the temperature exponent, and
!                L=3 is the activation energy (Kelvins).
!     NcRV  - Starting address of an array of NPAR reverse Arrhenius
!             parameters for the NREV reactions.
!             RCKWRK(NcRV+(N-1)*NPAR+(L-1)) is the Lth reverse
!             parameter for the Nth reaction with reverse parameters
!             defined, where
!                L=1 is the pre-exponential factor (mole-cm-sec-K),
!                L=2 is the temperature exponent, and
!                L=3 is the activation energy (Kelvins).
!             The reaction number is ICKWRK(IcRV+N-1).
!     NcLT  - Starting location of an array of the NLAR parameters for
!             the NLAN Landau-Teller reactions.
!             RCKWRK(NcLT+(N-1)*NLAR+(L-1)) is the Lth Landau-Teller
!             parameter for the Nth Landau-Teller reaction, where
!                L=1 is B(I) (Eq. 72) (Kelvins**1/3), and
!                L=2 is C(I) (Eq. 72) (Kelvins**2/3).
!             The reaction number is ICKWRK(IcLT+N-1).
!     NcRL  - Starting location of an array of the NLAR reverse
!             parameters for the NRLT Landau-Teller reactions for which
!             reverse parameters were given.
!             RCKWRK(NcRL+(N-1)*NLAR+(L-1)) is the Lth reverse
!             parameter for the Nth reaction with reverse Landau-Teller
!             parameters, where
!                L=1 is B(I) (Eq. 72) (Kelvins**1/3), and
!                L=2 is C(I) (Eq. 72) (Kelvins**2/3).
!             The reaction number is ICKWRK(IcRL+N-1).
!     NcFL  - Starting location of an array of the NFAR fall-off
!             parameters for the NFL fall-off reactions.
!             RCKWRK(NcFL+(N-1)*NFAR+(L-1)) is the Lth fall-off
!             parameter for the Nth fall-off reaction, where the low
!             pressure limits are defined by
!                L=1 is the pre-exponential factor (mole-cm-sec-K),
!                L=2 is the temperature exponent, and
!                L=3 is the activation energy (Kelvins).
!             Additional parameters define the centering, depending on
!             the type of formulation -
!                Troe: L=4 is the Eq. 68 parameter a,
!                      L=5 is the Eq. 68 parameter T*** (Kelvins),
!                      L=6 is the Eq. 68 parameter T*   (Kelvins), and
!                      L=7 is the Eq. 68 parameter T**  (Kelvins).
!                SRI:  L=4 is the Eq. 69 parameter a,
!                      L=5 is the Eq. 69 parameter b (Kelvins),
!                      L=6 is the Eq. 69 parameter c (kelvins),
!                      L=7 is the Eq. 69 parameter d, and
!                      L=8 is the Eq. 69 parameter e.
!             The reaction number is ICKWRK(IcFL+N-1), and the type
!             of formulation is ICKWRK(IcFO+N-1).
!     NcWL  - Starting location of an array of wavelengths for the NWL
!             wavelength-enhanced reactions.
!             RCKWRK(NcWL+N-1) is the wavelength enhancement (angstrom)
!             for the Nth wavelength-enhanced reaction;
!             the reaction number is ICKWRK(IcWL+N-1).
!     NcKT  - Starting location of an array of MXTB enhancement factors
!             for the NTHB third-body reactions.
!             RCKWRK(NcKT+(N-1)*MXTB+(L-1)) is the enhancement factor
!             for the Lth enhanced species in the Nth third-body
!             reaction;
!             the reaction number is ICKWRK(IcTB+N-1), and the Lth
!             enhanced species index number is
!             ICKWRK(IcKT+(N-1)*MXTB+L-1).
!     NcRU  - RCKWRK(NcRU) is the universal gas constant (ergs/mole-K).
!     NcRC  - RCKWRK(NcRC) is the universal gas constant (cal/mole-K).
!     NcPA  - RCKWRK(NcPA) is the pressure of one standard atmosphere
!             (dynes/cm**2).
!     NcKF  - Starting address of an array of intermediate forward
!             temperature-dependent rates for the II reactions.
!     NcKR  - Starting address of an array of intermediate reverse
!             temperature-dependent rates for the II reactions.
!     NcK1  - Starting addresses of arrays of internal work space
!     NcK2
!     NcK3                  space of length NKK
!     NcK4
!     NcI1  - Starting addresses of arrays of internal work space
!     NcI2
!     NcI3                  space of length NII
!     NcI4
      integer ::   NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,  &
                   NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,  &
                   NTHB, NRLT, NWL,  IcMM, IcKK, IcNC, IcPH, IcCH,  &
                   IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL, IcRV,  &
                   IcWL, IcFL, IcFO, IcKF, IcTB, IcKN, IcKT, NcAW,  &
                   NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL,  &
                   NcKT, NcWL, NcRU, NcRC, NcPA, NcKF, NcKR, NcK1,  &
                   NcK2, NcK3, NcK4, NcI1, NcI2, NcI3, NcI4

!     ** MACH block - Machine constants *******************************
      DOUBLE PRECISION :: SMALL, BIG, EXPARG

!     ** CKCONS block - Model constants *******************************
      character(len=16) :: PREC, VERS
      integer :: LENI, LENR, LENC
      logical :: KERR
!      COMMON /CKCONS/ PREC, VERS, KERR, LENI, LENR, LENC

!     ** CHEMKIN OUTPUT FILE
      character(len=*), parameter :: chemdat = 'chem.out'

!     *****************************************************************
!     **                CHEMKIN SUBROUTINES                          **
!     *****************************************************************
      contains

!  The work arrays contain all the pertinent information about the
!  species and the reaction mechanism.  They also contain some work
!  space needed by various routines for internal manipulations.  If a
!  user wishes to modify a CKLIB subroutine or to write new routines,
!  he will probably want to use the work arrays directly.  The starting
!  adddresses for information stored in the work arrays are found in
!  the labeled common block, COMMON /CKSTRT/, and are explained below.
!
!  COMMON /CKSTRT/ NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,
! 1                NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,
! 2                NTHB, NRLT, NWL,  IcMM, IcKK, IcNC, IcPH, IcCH,
! 3                IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL, IcRV,
! 4                IcWL, IcFL, IcFO, IcKF, IcTB, IcKN, IcKT, NcAW,
! 5                NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL,
! 6                NcKT, NcWL, NcRU, NcRC, NcPA, NcKF, NcKR, NcK1,
! 7                NcK2, NcK3, NcK4, NcI1, NcI2, NcI3, NcI4
!
!  INDEX CONSTANTS.
!

!  The linking file consists of the following binary records:
!   1) Information about the linking file:  VERS, PREC, KERR
!      Where VERS   = character*16 string representing the version
!                     number of the interpreter which created the
!                     the linking file.
!            PREC   = character*16 string representing the machine
!                     precision of the linking file (SINGLE, DOUBLE).
!            KERR   = logical which indicates whether or not
!                    an error occurred in the interpreter input.
!   2) Index constants:
!      LENI, LENR, LENC, NMM,  NKK,  NII,  MXSP, MXTB,
!      MXTP, NCP,  NPAR, NLAR, NFAR, NREV, NFAL, NTHB,
!      NLAN, NRLT, NWL, NCHRG
!      Where LENI = required length of ICKWRK.
!            LENR = required length of RCKWRK.
!            LENC = required length of CCKWRK.
!            NCHRG= total number of species with an electronic
!                   charge not equal to zero.

!  3) Element information:
!     ((CCKWRK(IcMM + M-1),                       !element names
!      RCKWRK(NcAW + M-1)),                       !atomic weights
!      M=1,NMM)

!  4) Species information:
!     ((CCKWRK(IcKK+K-1),                         !species names
!      (ICKWRK(IcNC+(K-1)*NMM+M-1),M=1,MMM),      !composition
!      ICKWRK(IcPH+K-1),                          !phase
!      ICKWRK(IcCH+K-1),                          !charge
!      RCKWRK(NcWT+K-1),                          !molec weight
!      ICKWRK(IcNT+K-1),                          !# of fit temps
!      (RCKWRK(NcTT+(K-1)*MXTP + L-1),L=1,MXTP),  !array of temps
!      ((RCKWRK(NcAA+(L-1)*NCP2+(K-1)*NCP2T+N-1), !fit coeff'nts
!               N=1,NCP2), L=1,(MXTP-1))),
!      K = 1,NKK)

!  5) Reaction information (if NII>0):
!     (ICKWRK(IcNS+I-1),                          !# of species
!      ICKWRK(IcNR+I-1),                          !# of reactants
!      (RCKWRK(NcCO+(I-1)*NPAR+N-1), N=1,NPAR),   !Arr. coefficients
!      (ICKWRK(IcNU+(I-1)*MXSP+N-1),              !stoic coef
!      ICKWRK(IcNK+(I-1)*MXSP+N-1), N=1,MXSP),    !species numbers
!      I = 1,NII)

!  6) Reverse parameter information (if NREV>0):
!     (ICKWRK(IcRV+N-1),                          !reaction numbers
!      (RCKWRK(NcRV+(N-1)*NPAR+L-1),L=1,NPAR),    !reverse coefficients
!      N = 1,NREV)

!  7) Fall-off reaction information (if NFAL>0):
!     (ICKWRK(IcFL+N-1),                          !reaction numbers
!      ICKWRK(IcFO+N-1),                          !fall-off option
!      ICKWRK(IcKF+N-1),                          !3rd-body species
!      (RCKWRK(NcFL+(N-1)*NFAR+L-1),L=1,NFAR),    !fall-off parameters
!      N=1,NFAL)

!  8) Third-body reaction information (if NTHB>0):
!     (ICKWRK(IcTB+N-1),                          !reaction numbers
!      ICKWRK(IcKN+N-1),                          !# of 3rd bodies
!      (ICKWRK(IcKT+(N-1)*MXTB+L-1),              !3rd-body species
!      RCKWRK(NcKT+(N-1)*MXTB+L-1),L=1,MXTB),     !enhancement factors
!      N=1,NTHB)

!  9) Landau-Teller reaction information (if NLAN>0):
!     (ICKWRK(IcLT+N-1),                          !reaction numbers
!      (RCKWRK(NcLT+(N-1)*NLAR+L-1),L=1,NLAR),    !L-T parameters
!      N=1,NLAN)

! 10) Reverse Landau-Teller reaction information (if NRLT>0):
!     (ICKWRK(IcRL+N-1),                          !reaction numbers
!      (RCKWRK(NcRL+(N-1)*NLAR+L-1),L=1,NLAR),    !rev. L-T parameters
!      N=1,NRLT)

! 11) Photon radiation reaction information (if NWL>0):
!     (ICKWRK(IcWL+N-1),                          !reaction numbers
!      RCKWRK(NcWL+N-1),                          !wavelength factor
!      N=1,NWL)

! *********************************************************************

	     SUBROUTINE CKABE  (ICKWRK, RCKWRK, RA, RB, RE)
!	     Returns the Arrhenius coefficients of the reactions;
!	     see Eq. (52).
!
!	     INPUT
!	     ICKWRK - Array of integer workspace.
!	     RCKWRK - Array of real work space.
!
!	  OUTPUT
!	     RA     - Pre-exponential constants for the reactions.
!	     RB     - Temperature dependence exponents for the reactions.
!	     RE     - Activation energies for the reactions.

	      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
	      DIMENSION RA(*), RB(*), RE(*), ICKWRK(*), RCKWRK(*)

	      DO 100 I = 1, NII
	         IND = NcCO + (I-1)*(NPAR+1)
	         RA(I) = RCKWRK(IND)
	         RB(I) = RCKWRK(IND+1)
	         RE(I) = RCKWRK(IND+2)
	  100 CONTINUE

	      RETURN
	      END SUBROUTINE CKABE

!         *************************************************************

      SUBROUTINE CKABML (P, T, X, ICKWRK, RCKWRK, ABML)
!	     Returns the Helmholtz free energy of the mixture in molar units,
!	     given the pressure, temperature, and mole fractions;
!	     see Eq. (46).
!
!	     INPUT
!	     P      - Pressure. [dynes/cm2]
!	     T      - Temperature. [K]
!	     X      - Mole fractions of the species.
!	     ICKWRK - Array of integer workspace.
!	     RCKWRK - Array of real work space.
!	     OUTPUT
!	     ABML   - Mean Helmholtz free energy in molar units. [erg/mol]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION X(*), ICKWRK(*), RCKWRK(*)

      CALL CKSML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKUML (T, ICKWRK, RCKWRK, RCKWRK(NcK2))
      RLNP = RCKWRK(NcRU) * LOG(P / RCKWRK(NcPA))

      ABML = 0.0
      DO 100 K = 1, NKK
         ABML = ABML + X(K) * ( RCKWRK(NcK2 + K - 1) - T *  &
                (RCKWRK(NcK1 + K - 1) - RCKWRK(NcRU)        &
                * LOG(MAX(X(K),SMALL)) - RLNP) )
  100 CONTINUE
      RETURN
      END SUBROUTINE CKABML

!     *****************************************************************

      SUBROUTINE CKABMS (P, T, Y, ICKWRK, RCKWRK, ABMS)
!     Returns the mean Helmholtz free energy of the mixture in
!     mass units, given the pressure, temperature and mass fractions;
!     see Eq. (47).
!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature.
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     ABMS   - Mean Helmholtz free energy in mass units. [erg/g]

      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION Y(*), ICKWRK(*), RCKWRK(*)

      CALL CKSML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKUML (T, ICKWRK, RCKWRK, RCKWRK(NcK2))
      CALL CKYTX (Y, ICKWRK, RCKWRK, RCKWRK(NcK3))

      RLNP = RCKWRK(NcRU) * LOG (P / RCKWRK(NcPA))

      SUM = 0.d0
      DO 100 K = 1, NKK
         SUM = SUM + RCKWRK(NcK3 + K - 1) *                         &
                   ( RCKWRK(NcK2 + K - 1) - T *                     &
                   ( RCKWRK(NcK1 + K - 1) -                         &
                     RCKWRK(NcRU)*                                  &
                     LOG(MAX(RCKWRK(NcK3 + K - 1),SMALL)) - RLNP))
  100 CONTINUE

      CALL CKMMWY (Y, ICKWRK, RCKWRK, WTM)
      ABMS = SUM / WTM
      RETURN
      END SUBROUTINE CKABMS

!     *****************************************************************


      SUBROUTINE CKAML  (T, ICKWRK, RCKWRK, AML)
!     Returns the standard state Helmholtz free energies in molar
!     units;  see Eq. (25).
!     INPUT
!     T      - Temperature.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     AML    - Standard state Helmholtz free energies in molar units
!              for the species. [erg/mol]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), AML(*)

      CALL CKSML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKHML (T, ICKWRK, RCKWRK, RCKWRK(NcK2))

      RUT = T*RCKWRK(NcRU)
      DO 150 K = 1, NKK
         AML(K) = RCKWRK(NcK2 + K - 1) - RUT - T*RCKWRK(NcK1 + K - 1)
150   CONTINUE
      RETURN
      END SUBROUTINE CKAML

!     *****************************************************************

      SUBROUTINE CKAMS  (T, ICKWRK, RCKWRK, AMS)
!     Returns the standard state Helmholtz free energies in mass
!     units;  see Eq. (32).
!     INPUT
!     T      - Temperature.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     AMS    - Standard state Helmholtz free energies in mass units
!              for the species. [erg/g]

      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), AMS(*)

      CALL CKSMS (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKHMS (T, ICKWRK, RCKWRK, RCKWRK(NcK2))

      RUT = T*RCKWRK(NcRU)
      DO 150 K = 1, NKK
         AMS(K) = RCKWRK(NcK2 + K - 1) - RUT/RCKWRK(NcWT + K - 1) &
                                       - T*RCKWRK(NcK1 + K - 1)
150   CONTINUE
      RETURN
      END SUBROUTINE CKAMS


!     *****************************************************************

     SUBROUTINE CKATHM (NDIM1, NDIM2, ICKWRK, RCKWRK, MAXTP, NT, TMP, A)
!     Returns the coefficients of the fits for thermodynamic properties
!     of the species; see Eqns. (19)-(21).

!     INPUT
!     NDIM1  - First dimension of the three-dimensional array of
!              thermodynamic fit coefficients, A; NDIM1 must be at
!              least NCP2, the total number of coefficients for one
!              temperature range.
!     NDIM2  - Second dimension of the three-dimensionalarray of
!              thermodynamic fit coefficients, A; NDIM2 must be at
!              least MXPT-1, the total number of temperature ranges.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     NT     - Number of temperatures used for fitting coefficients of
!              thermodynamic properties for the species.
!     TMP    - Common temperatures dividing the thermodynamic fits for
!              the species. [K]
!                   Dimension TMP(MAXT,*) exactly MAXT for the first
!                   dimension (the maximum number of temperatures
!                   allowed for a species) , and at least KK for the
!                   second dimension (the total number of species)
!     A      - Three dimensional array of fit coefficients to the
!              thermodynamic data for the species.
!              The indicies in  A(N,L,K) mean-
!              N = 1,NN are polynomial coefficients in CP/R
!                CP/R(K)=A(1,L,K) + A(2,L,K)*T + A(3,L,K)*T**2 + ...
!              N = NN+1 is a6 in Eq. (20)
!              N = NN+2 is a7 in Eq. (21)
!              L = 1..MXTP-1 is for each temperature range.
!              K  is  the  species index
!                   Dimension A(NPCP2,NDIM2,*) exactly NPCP2 and MXTP-1
!                   for the first and second dimensions and at least
!                   KK for the third.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION NT(*), TMP(MAXTP,*), A(NDIM1,NDIM2,*), &
                ICKWRK(*), RCKWRK(*)

      DO 100 K = 1, NKK
         NT(K) = ICKWRK(IcNT + K - 1)
  100 CONTINUE

      DO 140 L = 1, MXTP
         DO 140 K = 1, NKK
            TMP(L,K) = RCKWRK(NcTT + (K-1)*MXTP + L - 1)
  140 CONTINUE

      DO 150 K = 1, NKK
         DO 150 L = 1, MXTP-1
            NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
            DO 150 M = 1, NCP2
               A(M, L, K) = RCKWRK(NA1 + M - 1)
150   CONTINUE

      RETURN
      END SUBROUTINE CKATHM

!     *****************************************************************

      SUBROUTINE CKAWT  (ICKWRK, RCKWRK, AWT)
!     Returns the atomic weights of the elements
!     INPUT
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     AWT    - Atomic weights of the elements. [g/mol]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), AWT(*)

      DO 100 M = 1, NMM
         AWT(M) = RCKWRK(NcAW + M - 1)
  100 CONTINUE

      RETURN
      END SUBROUTINE CKAWT

!     *****************************************************************

      SUBROUTINE CKCDC  (T, C, ICKWRK, RCKWRK, CDOT, DDOT)
!     Returns the molar creation and destruction rates of the species
!     given the temperature and molar concentrations;  see Eq. (73).
!     INPUT
!     T      - Temperature.
!     C      - Molar concentrations of the species. [mol/cm3]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     CDOT   - Chemical molar creation rates of the species. [mol/cm3/s]
!     DDOT   - Chemical molar destruction rates of the species.[mol/cm3/s]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), C(*), CDOT(*), DDOT(*)

      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),              &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),          &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,             &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT), &
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),       &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),             &
                   RCKWRK(NcI1))

      DO 50 K = 1, NKK
         RCKWRK(NcK1 + K - 1) = C(K)
   50 CONTINUE

      CALL CKRATX (NII, NKK, MXSP, MXTB, T, RCKWRK(NcK1), ICKWRK(IcNS),  &
                   ICKWRK(IcNU), ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO),     &
                   NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF), NFAR, &
                   RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),       &
                   RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF),             &
                   RCKWRK(NcKR), RCKWRK(NcI1), RCKWRK(NcI2),             &
                   RCKWRK(NcI3))

      CDOT(1:NKK) = 0.d0
      DDOT(1:NKK) = 0.d0

      DO 200 I = 1, NII
         RKF = RCKWRK(NcI1 + I -1)
         RKR = RCKWRK(NcI2 + I -1)
         DO 200 N = 1, 3
            NKR = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
            NUR = ABS(ICKWRK(IcNU + (I-1)*MXSP + N - 1))
            NKP = ICKWRK(IcNK + (I-1)*MXSP + N + 2)
            NUP = ICKWRK(IcNU + (I-1)*MXSP + N + 2)
            IF (NKR /= 0) THEN
               CDOT(NKR) = CDOT(NKR) + RKR*NUR
               DDOT(NKR) = DDOT(NKR) + RKF*NUR
            ENDIF
            IF (NKP /= 0) THEN
               CDOT(NKP) = CDOT(NKP) + RKF*NUP
               DDOT(NKP) = DDOT(NKP) + RKR*NUP
            ENDIF
  200 CONTINUE
      RETURN
      END SUBROUTINE CKCDC


!     *****************************************************************

      SUBROUTINE CKCDXP (P, T, X, ICKWRK, RCKWRK, CDOT, DDOT)
!     Returns the molar creation and destruction rates of the species
!     given pressure, temperature and mole fractions;  see Eq. (73).
!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature.
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     CDOT   - Chemical molar creation rates of the species. [mol/cm3/s]
!     DDOT   - Chemical molar destruction rates of the species.[mol/cm3/s]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

      DIMENSION ICKWRK(*), RCKWRK(*), X(*), CDOT(*), DDOT(*)

      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),              &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),          &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,             &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT), &
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),       &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),             &
                   RCKWRK(NcI1))

      CALL CKXTCP (P, T, X, ICKWRK, RCKWRK, RCKWRK(NcK1))

      CALL CKRATX (NII, NKK, MXSP, MXTB, T, RCKWRK(NcK1), ICKWRK(IcNS),  &
                   ICKWRK(IcNU), ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO),     &
                   NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF), NFAR, &
                   RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),       &
                   RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF),             &
                   RCKWRK(NcKR), RCKWRK(NcI1), RCKWRK(NcI2),             &
                   RCKWRK(NcI3))

      CDOT(1:NKK) = 0.d0
      DDOT(1:NKK) = 0.d0

      DO 200 I = 1, NII
         RKF = RCKWRK(NcI1 + I - 1)
         RKR = RCKWRK(NcI2 + I - 1)
         DO 200 N = 1, 3
            NKR = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
            NUR = ABS(ICKWRK(IcNU + (I-1)*MXSP + N - 1))
            NKP = ICKWRK(IcNK + (I-1)*MXSP + N + 2)
            NUP = ICKWRK(IcNU + (I-1)*MXSP + N + 2)
            IF (NKR .NE. 0) THEN
               CDOT(NKR) = CDOT(NKR) + RKR*NUR
               DDOT(NKR) = DDOT(NKR) + RKF*NUR
            ENDIF
            IF (NKP .NE. 0) THEN
               CDOT(NKP) = CDOT(NKP) + RKF*NUP
               DDOT(NKP) = DDOT(NKP) + RKR*NUP
            ENDIF
  200 CONTINUE
      RETURN
      END SUBROUTINE CKCDXP

!     *****************************************************************

      SUBROUTINE CKCDXR (RHO, T, X, ICKWRK, RCKWRK, CDOT, DDOT)
!     Returns the molar creation and destruction rates of the species
!     given the mass density, temperature and mole fractions;
!     see Eq. (73).
!     INPUT
!     RHO    - Mass density. [g/cm3]
!     T      - Temperature.
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     CDOT   - Chemical molar creation rates of the species. [mol/cm3/s]
!     DDOT   - Chemical molar destruction rates of the species.[mol/cm3/s]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), X(*), CDOT(*), DDOT(*)

      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),              &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),          &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,             &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT), &
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),       &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),             &
                   RCKWRK(NcI1))

      CALL CKXTCR (RHO, T, X, ICKWRK, RCKWRK, RCKWRK(NcK1))

      CALL CKRATX (NII, NKK, MXSP, MXTB, T, RCKWRK(NcK1), ICKWRK(IcNS),  &
                   ICKWRK(IcNU), ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO),     &
                   NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF), NFAR, &
                   RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),       &
                   RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF),             &
                   RCKWRK(NcKR), RCKWRK(NcI1), RCKWRK(NcI2),             &
                   RCKWRK(NcI3))

      CDOT(1:NKK) = 0.d0
      DDOT(1:NKK) = 0.d0

      DO 200 I = 1, NII
         RKF = RCKWRK(NcI1 + I - 1)
         RKR = RCKWRK(NcI2 + I - 1)
         DO 200 N = 1, 3
            NKR = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
            NUR = ABS(ICKWRK(IcNU + (I-1)*MXSP + N - 1))
            NKP = ICKWRK(IcNK + (I-1)*MXSP + N + 2)
            NUP = ICKWRK(IcNU + (I-1)*MXSP + N + 2)
            IF (NKR .NE. 0) THEN
               CDOT(NKR) = CDOT(NKR) + RKR*NUR
               DDOT(NKR) = DDOT(NKR) + RKF*NUR
            ENDIF
            IF (NKP .NE. 0) THEN
               CDOT(NKP) = CDOT(NKP) + RKF*NUP
               DDOT(NKP) = DDOT(NKP) + RKR*NUP
            ENDIF
  200 CONTINUE
      RETURN
      END SUBROUTINE CKCDXR



!     *****************************************************************

      SUBROUTINE CKCDYP (P, T, Y, ICKWRK, RCKWRK, CDOT, DDOT)
!     Returns the molar creation and destruction rates of the species
!     given mass density, temperature and mass fractions;
!     see Eq. (73).
!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature.
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     CDOT   - Chemical molar creation rates of the species. [mol/cm3/s]
!     DDOT   - Chemical molar destruction rates of the species.[mol/cm3/s]

      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

      DIMENSION ICKWRK(*), RCKWRK(*), Y(*), CDOT(*), DDOT(*)

      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),              &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),          &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,             &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT), &
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),       &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),             &
                   RCKWRK(NcI1))

      CALL CKYTCP (P, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK1))

      CALL CKRATX (NII, NKK, MXSP, MXTB, T, RCKWRK(NcK1), ICKWRK(IcNS),  &
                   ICKWRK(IcNU), ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO),     &
                   NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF), NFAR, &
                   RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),       &
                   RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF),             &
                   RCKWRK(NcKR), RCKWRK(NcI1), RCKWRK(NcI2),             &
                   RCKWRK(NcI3))

      DO 100 K = 1, NKK
         CDOT(K) = 0.0
         DDOT(K) = 0.0
  100 CONTINUE
      DO 200 I = 1, NII
         RKF = RCKWRK(NcI1 + I - 1)
         RKR = RCKWRK(NcI2 + I - 1)
         DO 200 N = 1, 3
            NKR = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
            NUR = ABS(ICKWRK(IcNU + (I-1)*MXSP + N - 1))
            NKP = ICKWRK(IcNK + (I-1)*MXSP + N + 2)
            NUP = ICKWRK(IcNU + (I-1)*MXSP + N + 2)
            IF (NKR .NE. 0) THEN
               CDOT(NKR) = CDOT(NKR) + RKR*NUR
               DDOT(NKR) = DDOT(NKR) + RKF*NUR
            ENDIF
            IF (NKP .NE. 0) THEN
               CDOT(NKP) = CDOT(NKP) + RKF*NUP
               DDOT(NKP) = DDOT(NKP) + RKR*NUP
            ENDIF
  200 CONTINUE
      RETURN
      END SUBROUTINE CKCDYP

!     *****************************************************************


      SUBROUTINE CKCDYR (RHO, T, Y, ICKWRK, RCKWRK, CDOT, DDOT)
!     Returns the molar creation and destruction rates of the species
!     given the mass density, temperature and mass fractions;
!     see Eq. (73).
!     INPUT
!     RHO    - Mass density. [g/cm3]
!     T      - Temperature.
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     CDOT   - Chemical molar creation rates of the species.[MOL/CM3/S]
!     DDOT   - Chemical molar destruction rates of the species.[mol/cm3/s]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), Y(*), CDOT(*), DDOT(*)

      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),              &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),          &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,             &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT), &
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),       &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),             &
                   RCKWRK(NcI1))

      CALL CKYTCR (RHO, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK1))

      CALL CKRATX (NII, NKK, MXSP, MXTB, T, RCKWRK(NcK1), ICKWRK(IcNS),  &
                   ICKWRK(IcNU), ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO),     &
                   NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF), NFAR, &
                   RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),       &
                   RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF),             &
                   RCKWRK(NcKR), RCKWRK(NcI1), RCKWRK(NcI2),             &
                   RCKWRK(NcI3))

      DO 100 K = 1, NKK
         CDOT(K) = 0.0
         DDOT(K) = 0.0
  100 CONTINUE
      DO 200 I = 1, NII
         RKF = RCKWRK(NcI1 + I - 1)
         RKR = RCKWRK(NcI2 + I - 1)
         DO 200 N = 1, 3
            NKR = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
            NUR = ABS(ICKWRK(IcNU + (I-1)*MXSP + N - 1))
            NKP = ICKWRK(IcNK + (I-1)*MXSP + N + 2)
            NUP = ICKWRK(IcNU + (I-1)*MXSP + N + 2)
            IF (NKR .NE. 0) THEN
               CDOT(NKR) = CDOT(NKR) + RKR*NUR
               DDOT(NKR) = DDOT(NKR) + RKF*NUR
            ENDIF
            IF (NKP .NE. 0) THEN
               CDOT(NKP) = CDOT(NKP) + RKF*NUP
               DDOT(NKP) = DDOT(NKP) + RKR*NUP
            ENDIF
  200 CONTINUE
      RETURN
      END SUBROUTINE CKCDYR

!     *****************************************************************

      SUBROUTINE CKCHRG (ICKWRK, RCKWRK, KCHARG)
!     Returns the electronic charges of the species
!     INPUT
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     KCHARG - Electronic charges of the species; KCHARG(K)=-2
!              indicates that the Kth species has two excess electrons.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), KCHARG(*)

      DO 100 K = 1, NKK
         KCHARG(K) = ICKWRK(IcCH + K - 1)
  100 CONTINUE
      RETURN
      END SUBROUTINE CKCHRG

!     *****************************************************************

      SUBROUTINE CKCONT (K, Q, ICKWRK, RCKWRK, CIK)
!     Returns the contributions of the reactions to the molar
!     production rate of a species;  see Eqs. (49) and (51).
!     INPUT
!     K      - Integer species number.
!                   Data type - integer scalar
!     Q      - Rates of progress for the reactions. [mol/cm3/s]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     CIK    - Contributions of the reactions to the molar production
!              rate of the Kth species. [mol/cm3/s]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), Q(*), CIK(*)

      CIK(1:NII) = 0.d0

      DO 200 N = 1, MXSP
         DO 200 I = 1, NII
            NK = ICKWRK(IcNK + MXSP*(I-1) + N - 1)
            NC = ICKWRK(IcNU + MXSP*(I-1) + N - 1)
            IF (NK .EQ. K) CIK(I) = CIK(I) + NC*Q(I)
200   CONTINUE
      RETURN
      END SUBROUTINE CKCONT

!     *****************************************************************

      SUBROUTINE CKCOMP (IST, IRAY, II, I)
!     Returns the index of an element of a reference character
!     string array which corresponds to a character string;
!     leading and trailing blanks are ignored.
!     INPUT
!     IST   - A character string.
!     IRAY  - An array of character strings;
!     II    - The length of IRAY.
!     OUTPUT
!     I     - The first integer location in IRAY in which IST
!             corresponds to IRAY(I); if IST is not also an
!             entry in IRAY, I=0.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
      CHARACTER*(*) IST, IRAY(*)

      I = 0

      do N = II, 1, -1
         IS1 = IFIRCH(IST)
         IS2 = ILASCH(IST)
         IR1 = IFIRCH(IRAY(N))
         IR2 = ILASCH(IRAY(N))
         IF ( IS2 >= IS1 .AND. IS2 > 0 .AND.           &
              IR2 >= IR1 .AND. IR2 > 0 .AND.           &
              IST(IS1:IS2).EQ.IRAY(N)(IR1:IR2) ) I = N
      end do



      RETURN
      END SUBROUTINE CKCOMP

!     *****************************************************************

      SUBROUTINE CKCPBL (T, X, ICKWRK, RCKWRK, CPBML)
!     Returns the mean specific heat at constant pressure;
!     see Eq. (33).
!     INPUT
!     T      - Temperature. [K]
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     CPBML  - Mean specific heat at constant pressure in molar units. [erg/mol/K]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

      DIMENSION ICKWRK(*), RCKWRK(*), X(*)

      CALL CKCPML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))

      CPBML = 0.d0

      do K = 1, NKK
         CPBML = CPBML + X(K)*RCKWRK(NcK1 + K - 1)
      end do

      RETURN
      END SUBROUTINE CKCPBL



!     *****************************************************************

      SUBROUTINE CKCPBS (T, Y, ICKWRK, RCKWRK, CPBMS)
!     Returns the mean specific heat at constant pressure;
!     see Eq. (34).
!     INPUT
!     T      - Temperature. [K]
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     CPBMS  - Mean specific heat at constant pressure in mass units. [erg/g/K]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

      DIMENSION ICKWRK(*), RCKWRK(*), Y(*)

      CALL CKCPMS (T, ICKWRK, RCKWRK, RCKWRK(NcK1))

      CPBMS = 0.d0
      DO 100 K = 1, NKK
         CPBMS = CPBMS + Y(K)*RCKWRK(NcK1 + K - 1)
  100 CONTINUE

      RETURN
      END SUBROUTINE CKCPBS


!     *****************************************************************

      SUBROUTINE CKCPML (T, ICKWRK, RCKWRK, CPML)
!     Returns the specific heats at constant pressure in molar units
!     INPUT
!     T      - Temperature. [K]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     CPML   - Specific heats at constant pressure in molar units
!              for the species. [erg/mol/K]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

      DIMENSION ICKWRK(*), RCKWRK(*), CPML(*), TN(10)

      TN(1) = 1.d0
      DO 150 N = 2, NCP
         TN(N) = T**(N-1)
150   CONTINUE

      DO 250 K = 1, NKK
         L = 1
         DO 220 N = 2, ICKWRK(IcNT + K - 1)-1
            TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE

         NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
         CPML(K) = 0.0
         DO 250 N = 1, NCP
            CPML(K) = CPML(K) + RCKWRK(NcRU)*TN(N)*RCKWRK(NA1 + N - 1)
250   CONTINUE

      RETURN
      END SUBROUTINE CKCPML


!     *****************************************************************

      SUBROUTINE CKCPMS (T, ICKWRK, RCKWRK, CPMS)
!     Returns the specific heats at constant pressure in mass units;
!     see Eq. (26).
!     INPUT
!     T      - Temperature.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     CPMS   - Specific heats at constant pressure in mass units
!              for the species. [erg/g/K]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), CPMS(*), TN(10)

      TN(1) = 1.0
      DO 150 N = 2, NCP
         TN(N) = T**(N-1)
150   CONTINUE

      DO 250 K = 1, NKK
         L = 1
         DO 220 N = 2, ICKWRK(IcNT + K - 1)-1
            TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE

         NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
         SUM = 0.0
         DO 240 N = 1, NCP
            SUM = SUM + TN(N)*RCKWRK(NA1 + N - 1)
  240    CONTINUE
         CPMS(K) = RCKWRK(NcRU) * SUM / RCKWRK(NcWT + K - 1)

250   CONTINUE
      RETURN
      END SUBROUTINE CKCPMS

!     *****************************************************************

      SUBROUTINE CKCPOR (T, ICKWRK, RCKWRK, CPOR)
!     Returns the nondimensional specific heats at constant pressure;
!     see Eq. (19).
!     INPUT
!     T      - Temperature.[K]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     CPOR   - Nondimensional specific heats at constant pressure
!              for the species.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

      DIMENSION ICKWRK(*), RCKWRK(*), TN(10), CPOR(*)

      TN(1) = 1.0
      DO 150 N = 2, NCP
         TN(N) = T**(N-1)
150   CONTINUE

      DO 250 K = 1, NKK
         L = 1
         DO 220 N = 2, ICKWRK(IcNT + K - 1)-1
            TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE

         NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
         CPOR(K) = 0.0
         DO 250 N = 1, NCP
            CPOR(K) = CPOR(K) + TN(N)*RCKWRK(NA1 + N - 1)
250   CONTINUE
      RETURN
      END SUBROUTINE CKCPOR

!     *****************************************************************

      SUBROUTINE CKCRAY (LINE, NN, KRAY, LOUT, NDIM, NRAY, NF, KERR)
!     This subroutine is called to parse a character string, LINE,
!     that is composed of several blank-delimited substrings.  Each
!     substring in LINE is compared with an ordered reference array
!     of character strings, KRAY(*).  For each substring in LINE that
!     is also an entry in KRAY(*), the index position in KRAY(*) is
!     returned in the integer array, NRAY(*).  It is expected that
!     each substring in LINE will be found in KRAY(*). If a substring
!     cannot be found in KRAY(*) an error flag will be returned. For
!     example, after reading a line of species names, the subroutine
!     might be called to assign Chemkin species index numbers to the
!     list of species names.  This application is illustrated in the
!     following example:
!
!     input:  LINE    = "OH  N2  NO"
!             KRAY(*) = "H2" "O2" "N2" "H" "O" "N" "OH" "H2O" "NO"
!             NN      = 9, the number of entries in KRAY(*)
!             LOUT    = 6, a logical unit number on which to write
!                       diagnostic messages.
!             NDIM    = 10, the dimension of array NRAY(*)
!     output: NRAY(*) = 7, 3, 9, the index numbers of the entries
!                       in KRAY(*) corresponding to the substrings
!                       in LINE
!             NF      = 3, the number of correspondences found.
!             KERR    = .FALSE.
!     INPUT
!     LINE - A character string.
!                 Data type - CHARACTER*(*)
!     KRAY - An array of character strings; dimension KRAY(*) at
!                 least NN.
!                 Data type - CHARACTER*(*)
!     NN   - Total number of character strings in KRAY
!                 Data type - integer scalar
!     LOUT - Output unit for printed diagnostics
!                 Data type - integer scalar
!     NDIM - Dimension of the integer array NRAY.
!                 Data type - integer scalar
!     OUTPUT
!     NRAY - Index numbers of the elements of KRAY which
!            correspond to the substrings in LINE.
!                 Data type - integer array
!                 Dimension NRAY(*) at least NDIM
!     NF   - Number of correspondences found.
!                 Data type - integer scalar
!     KERR - Error flag; syntax or dimensioning errors will
!            result in KERR=.TRUE.
!                 Data type - logical
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      CHARACTER LINE*(*), KRAY(*)*(*), SUB(80)*80
      DIMENSION NRAY(*)
      LOGICAL KERR, IERR

      KERR = .FALSE.
      NF = 0

      IDIM = 1000
      CALL CKSUBS (LINE, LOUT, IDIM, SUB, NFOUND, IERR)
      IF (IERR) THEN
         KERR = .TRUE.
         WRITE (LOUT,*) ' Error in CKCRAY...'
         RETURN
      ENDIF

      DO 50 N = 1, NFOUND
         CALL CKCOMP (SUB(N), KRAY, NN, K)
         IF (K .LE. 0) THEN
            LT = MAX (ILASCH(SUB(N)), 1)
            WRITE (LOUT,'(A)') &
     &      ' Error in CKCRAY...'//SUB(N)(:LT)//' not found...'
            KERR = .TRUE.
         ELSE
            IF (NF+1 .GT. NDIM) THEN
               WRITE (LOUT,'(A)') &
     &       ' Error in CKCRAY...dimension of NRAY too small...'
               KERR = .TRUE.
            ELSE
               NF = NF + 1
               NRAY(NF) = K
            ENDIF
         ENDIF
   50 CONTINUE
      RETURN
      END SUBROUTINE CKCRAY

!     *****************************************************************

      SUBROUTINE CKCTC  (T, C, ICKWRK, RCKWRK, CDOT, TAU)
!     Returns the molar creation rates and characteristic destruction
!     times of the species given temperature and molar concentrations;
!     see Eqs. (76) and (78).
!     INPUT
!     T      - Temperature. [K]
!     C      - Molar concentrations of the species. [mol/cm3]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     CDOT   - Chemical molar creation rates of the species. [mol/cm3/s]
!     TAU    - Characteristic destruction times of the species. [s]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), C(*), TAU(*), CDOT(*)

      CALL CKCDC (T, C, ICKWRK, RCKWRK, CDOT, RCKWRK(NcK1))
      DO 150 K = 1, NKK
         TAU(K) = C(K) / (RCKWRK(NcK1 + K - 1)+SMALL)
150   CONTINUE
      RETURN
      END SUBROUTINE CKCTC

!     *****************************************************************

      SUBROUTINE CKCTX  (C, ICKWRK, RCKWRK, X)
!     Returns the mole fractions given the molar concentrations;
!     see Eq. (13).
!     INPUT
!     C      - Molar concentrations of the species. [mol/cm3]
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!     OUTPUT
!     X      - Mole fractions of the species.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), C(*), X(*)

      CTOT = 0.0
      DO 100 K = 1, NKK
         CTOT = CTOT + C(K)
  100 CONTINUE
      DO 200 K = 1, NKK
         X(K) = C(K)/CTOT
200   CONTINUE
      RETURN
      END SUBROUTINE CKCTX

!     *****************************************************************

      SUBROUTINE CKCTXP (P, T, X, ICKWRK, RCKWRK, CDOT, TAU)
!     Returns the molar creation rates and characteristic destruction
!     times of the species given the pressure, temperature and mole
!     fractions;  see Eqs. (76) and (78).
!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature. [K]
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     CDOT   - Chemical molar creation rates of the species. [mol/cm3/s]
!     TAU    - Characteristic destruction times of the species. [s]

      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

      DIMENSION ICKWRK(*), RCKWRK(*), X(*), TAU(*), CDOT(*)

      CALL CKCDXP (P, T, X, ICKWRK, RCKWRK, CDOT, RCKWRK(NcK1))
      CALL CKXTCP (P, T, X, ICKWRK, RCKWRK, RCKWRK(NcK2))
      DO 150 K = 1, NKK
         TAU(K) = RCKWRK(NcK2 + K - 1) / (RCKWRK(NcK1 + K - 1)+SMALL)
150   CONTINUE
      RETURN
      END SUBROUTINE CKCTXP

!     *****************************************************************

      SUBROUTINE CKCTXR (RHO, T, X, ICKWRK, RCKWRK, CDOT, TAU)
!     Returns the molar creation rates and characteristic destruction
!     times of the species given the mass density, temperature and
!     mole fractions;  see Eqs. (76) and (78).
!     INPUT
!     RHO    - Mass density. [g/cm3]
!     T      - Temperature. [K]
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     CDOT   - Chemical molar creation rates of the species. [mol/cm3/s]
!     TAU    - Characteristic destruction times of the species. [s]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), X(*), TAU(*), CDOT(*)

      CALL CKCDXR (RHO, T, X, ICKWRK, RCKWRK, CDOT, RCKWRK(NcK1))
      CALL CKXTCR (RHO, T, X, ICKWRK, RCKWRK, RCKWRK(NcK2))
      DO 150 K = 1, NKK
         TAU(K) = RCKWRK(NcK2 + K - 1) / (RCKWRK(NcK1 + K - 1)+SMALL)
150   CONTINUE
      RETURN
      END SUBROUTINE CKCTXR

!     *****************************************************************

      SUBROUTINE CKCTY  (C, ICKWRK, RCKWRK, Y)
!     Returns the mass fractions given the molar concentrations;
!     see Eq. (12).
!     INPUT
!     C      - Molar concentrations of the species[mol/cm3]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     Y      - Mass fractions of the species.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), C(*), Y(*)

      RHO = 0.0
      DO 100 K = 1, NKK
         RHO = RHO + C(K)*RCKWRK(NcWT + K - 1)
  100 CONTINUE

      DO 200 K = 1, NKK
         Y(K) = C(K)*RCKWRK(NcWT + K - 1)/RHO
200   CONTINUE
      RETURN
      END SUBROUTINE CKCTY

!     *****************************************************************

      SUBROUTINE CKCTYP (P, T, Y, ICKWRK, RCKWRK, CDOT, TAU)
!     Returns the molar creation rates and characteristic destruction
!     times of the species given the mass density, temperature and
!     mass fractions;  see Eqs. (76) and (78).
!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature.
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     CDOT   - Chemical molar creation rates of the species.[mol/cm3/s]
!     TAU    - Characteristic destruction times of the species. [s]

      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), Y(*), TAU(*), CDOT(*)
      CALL CKCDYP (P, T, Y, ICKWRK, RCKWRK, CDOT, RCKWRK(NcK1))
      CALL CKYTCP (P, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK2))
      DO 150 K = 1, NKK
         TAU(K) = RCKWRK(NcK2 + K - 1) / (RCKWRK(NcK1 + K - 1)+SMALL)
150   CONTINUE
      RETURN
      END SUBROUTINE CKCTYP

!     *****************************************************************


      SUBROUTINE CKCTYR (RHO, T, Y, ICKWRK, RCKWRK, CDOT, TAU)
!     Returns the molar creation rates and characteristic destruction
!     times of the species given the mass density, temperature and
!     mass fractions;  see Eqs. (76) and (78).
!     INPUT
!     RHO    - Mass density. [g/cm3]
!     T      - Temperature. [K]
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     CDOT   - Chemical molar creation rates of the species [mol/cm3]
!     TAU    - Characteristic destruction times of the species. [s]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

      DIMENSION ICKWRK(*), RCKWRK(*), Y(*), TAU(*), CDOT(*)

      CALL CKCDYR (RHO, T, Y, ICKWRK, RCKWRK, CDOT, RCKWRK(NcK1))
      CALL CKYTCR (RHO, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK2))
      DO 150 K = 1, NKK
         TAU(K) = RCKWRK(NcK2 + K - 1) / (RCKWRK(NcK1 + K - 1)+SMALL)
150   CONTINUE
      RETURN
      END SUBROUTINE CKCTYR

!     *****************************************************************

      SUBROUTINE CKCVBL (T, X, ICKWRK, RCKWRK, CVBML)
!     Returns the mean specific heat at constant volume in molar units;
!     see Eq. (35).
!
!     INPUT
!     T      - Temperature.
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!
!     OUTPUT
!     CVBML  - Mean specific heat at constant volume in molar units.
!                   cgs units - ergs/(mole*K)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), X(*)

      CALL CKCVML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))

      CVBML = 0.0
      DO 100 K = 1, NKK
         CVBML = CVBML + X(K)*RCKWRK(NcK1 + K - 1)
  100 CONTINUE
      RETURN
      END SUBROUTINE CKCVBL


!     *****************************************************************

      SUBROUTINE CKCVBS (T, Y, ICKWRK, RCKWRK, CVBMS)
!     Returns the mean specific heat at constant volume in mass units;
!     see Eq. (36).
!     INPUT
!     T      - Temperature. [K]
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     CVBMS  - Mean specific heat at constant volume in mass units.
!                   cgs units - ergs/(gm*K)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

      DIMENSION ICKWRK(*), RCKWRK(*), Y(*)

      CALL CKCVMS (T, ICKWRK, RCKWRK, RCKWRK(NcK1:NcK1+NKK-1))

      CVBMS = 0.d0
      DO 100 K = 1, NKK
         CVBMS = CVBMS + Y(K)*RCKWRK(NcK1 + K - 1)
  100 CONTINUE
      RETURN
      END SUBROUTINE CKCVBS

!     *****************************************************************

      SUBROUTINE CKCVML (T, ICKWRK, RCKWRK, CVML)
!     Returns the specific heats in constant volume in molar units;
!     see Eq. (22).
!     INPUT
!     T      - Temperature. [K]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     CVML   - Specific heats at constant volume in molar units
!              for the species. [erg/mol/K]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

      DIMENSION ICKWRK(*), RCKWRK(*), CVML(*)

      CALL CKCPML (T, ICKWRK, RCKWRK, CVML)

      DO 150 K = 1, NKK
         CVML(K) = CVML(K) - RCKWRK(NcRU)
150   CONTINUE
      RETURN
      END SUBROUTINE CKCVML

!     *****************************************************************


      SUBROUTINE CKCVMS (T, ICKWRK, RCKWRK, CVMS)
!     Returns the specific heats at constant volume in mass units;
!     see Eq. (29).
!     INPUT
!     T      - Temperature. [K]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     CVMS   - Specific heats at constant volume in mass units
!              for the species. [erg/g/K]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), CVMS(*)

      CALL CKCPMS (T, ICKWRK, RCKWRK, CVMS)

      DO 150 K = 1, NKK
         CVMS(K) = CVMS(K) - RCKWRK(NcRU) / RCKWRK(NcWT + K - 1)
150   CONTINUE
      RETURN
      END SUBROUTINE CKCVMS

!     *****************************************************************


      SUBROUTINE CKEQC  (T, C, ICKWRK, RCKWRK, EQKC)
!     Returns the equilibrium constants of the reactions given
!     temperature and molar concentrations;  see Eq. (54).
!     INPUT
!     T      - Temperature. [K]
!     C      - Molar concentrations of the species. [mol/cm3]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     EQKC   - Equilibrium constants in concentration units
!              for the reactions.
!                   cgs units - (mole/cm**3)**some power, depending on
!                               the reaction

      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), C(*), EQKC(*)

      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),              &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),          &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,             &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT), &
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),       &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),             &
                   RCKWRK(NcI1))

      DO 50 I = 1, NII
         EQKC(I) = RCKWRK(NcI1 + I - 1)
   50 CONTINUE

      RETURN
      END SUBROUTINE CKEQC


!     *****************************************************************

      SUBROUTINE CKEQXP (P, T, X, ICKWRK, RCKWRK, EQKC)
!     Returns the equilibrium constants for the reactions given
!     pressure, temperature and mole fractions;  see Eq. (54).
!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature.
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     EQKC   - Equilibrium constants in concentration units
!              for the reactions.
!                   cgs units - (mole/cm**3)**some power, depending on
!                               the reaction
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), X(*), EQKC(*)

      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),              &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),          &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,             &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT), &
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),       &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),             &
                   RCKWRK(NcI1))

      DO 50 I = 1, NII
         EQKC(I) = RCKWRK(NcI1 + I - 1)
   50 CONTINUE

      RETURN
      END SUBROUTINE CKEQXP

!     *****************************************************************

      SUBROUTINE CKEQXR (RHO, T, X, ICKWRK, RCKWRK, EQKC)
!     Returns the equilibrium constants of the reactions given mass
!     density, temperature and mole fractions;  see Eq. (54).
!     INPUT
!     RHO    - Mass density. [g/cm3]
!     T      - Temperature. [K]
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!
!     OUTPUT
!     EQKC   - Equilibrium constants in concentration units
!              for the reactions.
!                   cgs units - (mole/cm**3)**some power, depending on
!                               the reaction

      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), X(*), EQKC(*)
      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),              &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),          &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,             &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT), &
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),       &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),             &
                   RCKWRK(NcI1))

      DO 50 I = 1, NII
         EQKC(I) = RCKWRK(NcI1 + I - 1)
   50 CONTINUE

      RETURN
      END SUBROUTINE CKEQXR


!     *****************************************************************

      SUBROUTINE CKEQYP (P, T, Y, ICKWRK, RCKWRK, EQKC)
!     Returns the equilibrium constants for the reactions given
!     pressure, temperature and mass fractions;  see Eq. (54).
!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature. [K]
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!
!     OUTPUT
!     EQKC   - Equilibrium constants in concentration units
!              for the reactions.
!                   cgs units - (mole/cm**3)**some power, depending on
!                               the reaction
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), Y(*), EQKC(*)

      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),             &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),         &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,            &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT),&
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),      &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),            &
                   RCKWRK(NcI1))

      DO 50 I = 1, NII
         EQKC(I) = RCKWRK(NcI1 + I - 1)
   50 CONTINUE

      RETURN
      END SUBROUTINE CKEQYP

!     *****************************************************************

      SUBROUTINE CKEQYR (RHO, T, Y, ICKWRK, RCKWRK, EQKC)
!     Returns the equilibrium constants of the reactions given mass
!     density, temperature and mass fractions;  see Eq. (54).
!     INPUT
!     RHO    - Mass density. [g/cm3]
!     T      - Temperature. [K]
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     EQKC   - Equilibrium constants in concentration units
!              for the reactions.
!                   cgs units - (mole/cm**3)**some power, depending on
!                               the reaction
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

      DIMENSION ICKWRK(*), RCKWRK(*), Y(*), EQKC(*)

      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),             &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),         &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,            &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT),&
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),      &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),            &
                   RCKWRK(NcI1))

      DO 50 I = 1, NII
         EQKC(I) = RCKWRK(NcI1 + I - 1)
   50 CONTINUE

      RETURN
      END SUBROUTINE CKEQYR

!     *****************************************************************

      SUBROUTINE CKFAL  (NDIM, ICKWRK, RCKWRK, IFOP, KFAL, FPAR)
!     Returns a set of flags indicating whether a reaction has
!     fall-off behavior and an array of the fall-off
!     parameters.
!     INPUT
!     NDIM   - First dimension of the two dimensional array FPAR;
!              NDIM must be greater than or equal to the maximum
!              number of fall-off parameters, NFAR, which is
!              currently equal to 8.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!
!     OUTPUT
!     IFOP   - Array of flags indicating fall-off behavior:
!               0 - No fall-off behavior
!               1 - fall-off behavior - Lindeman form (3 parameters)
!               2 - fall-off behavior - SRI form      (8 parameters)
!               3 - fall-off behavior - Troe form     (6 parameters)
!               4 - fall-off behavior - Troe form     (7 parameters)
!                   Data type - integer array
!                   Dimension IFOP(*) at least II, the total number
!                   of reactions.
!     KFAL   - Array of flags indicating type of bath-gas
!              concentration to be used in fall-off expressions
!              (see footnote on page 23).
!               0 - Use total concentration of gas mixture
!                    (with the added capability of using enhanced
!                     third body coefficients) (default)
!               K - Use the concentration of species K
!                   Data type - integer array
!                   Dimension KFAL(*) at least II, the total number
!                   of reactions.
!     FPAR   - Matrix of fall-off parameters. The number of fall-off
!                   parameters will depend on the particular
!                   functional form indicated by the IFOP array:
!                   FPAR(1,I), FPAR(2,I), FPAR(3,I) are always the
!                   parameters entered on the LOW auxiliary keyword
!                   line in the CHEMKIN interpretor input file.
!                     FPAR(1,I) = Pre-exponential for low pressure
!                                 limiting rate constant
!                                 cgs units - mole-cm-sec-K
!                     FPAR(2,I) = Temperature dependence exponents
!                                 for the low pressure limiting rate
!                                 constants.
!                     FPAR(3,I) = Activation energy for the low
!                                 pressure limiting rate constant.
!                                 cgs units - Kelvins
!                   Additional FPAR values depend on IFOP:
!                   IFOP(I) = 2:
!                     FPAR(4,I) = a           (See Eqn. (69))
!                     FPAR(5,I) = b (Kelvin)  (See Eqn. (69))
!                     FPAR(6,I) = c (Kelvin)  (See Eqn. (69))
!                     FPAR(7,I) = d           (See Eqn. (69))
!                     FPAR(8,I) = e           (See Eqn. (69))
!                   IFOP(I) = 3:
!                     FPAR(4,I) = a             (See Eqn. (68))
!                     FPAR(5,I) = T*** (Kelvin) (See Eqn. (68))
!                     FPAR(6,I) = T*   (Kelvin) (See Eqn. (68))
!                   IFOP(I) = 4:
!                     FPAR(4,I) = a             (See Eqn. (68))
!                     FPAR(5,I) = T*** (Kelvin) (See Eqn. (68))
!                     FPAR(6,I) = T*   (Kelvin) (See Eqn. (68))
!                     FPAR(7,I) = T**  (Kelvin) (See Eqn. (68))
!                   Data type - real array
!                   Dimension FPAR(NDIM,*) exactly NDIM (at least NFAR,
!                   the maximum number of fall-off parameters
!                   - currently 8) for the first
!                   dimension and at least II for the second, the total
!                   number of reactions).
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

      DIMENSION ICKWRK(*), RCKWRK(*), IFOP(*), KFAL(*), FPAR(NDIM,*)

      DO 100 I = 1, NII
        IFOP(I) = 0
        KFAL(I) = 0
  100 CONTINUE

      DO 150 I = 1, NII
         DO 140 N = 1, NFAR
            FPAR(N,I) = 0.0
  140    CONTINUE
  150 CONTINUE

      DO 250 N = 1, NFAL
        I       = ICKWRK(IcFL + N - 1)
        IFOP(I) = ICKWRK(IcFO + N - 1)
        KFAL(I) = ICKWRK(IcKF + N - 1)
        IF (IFOP(I) .EQ. 1) THEN
          DO 210 L = 1, 3
            FPAR(L,I) = RCKWRK(NcFL + (N-1)*NFAR + L - 1)
  210     CONTINUE
        ELSE IF (IFOP(I) .EQ. 2) THEN
          DO 220 L = 1, 8
            FPAR(L,I) = RCKWRK(NcFL + (N-1)*NFAR + L - 1)
  220     CONTINUE
        ELSE IF (IFOP(I) .EQ. 3) THEN
          DO 230 L = 1, 6
            FPAR(L,I) = RCKWRK(NcFL + (N-1)*NFAR + L - 1)
  230     CONTINUE
        ELSE IF (IFOP(I) .EQ. 4) THEN
          DO 240 L = 1, 7
            FPAR(L,I) = RCKWRK(NcFL + (N-1)*NFAR + L - 1)
  240     CONTINUE
        ENDIF
  250 CONTINUE
      RETURN
      END SUBROUTINE CKFAL

!     *****************************************************************

      SUBROUTINE CKGBML (P, T, X, ICKWRK, RCKWRK, GBML)
!     Returns the mean Gibbs free energy of the mixture in molar units,
!     given the pressure, temperature and mole fractions;
!     see Eq. (44).
!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature. [K]
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     GBML   - Mean Gibbs free energy in molar units. [erg/mol]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), X(*)

      CALL CKSML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKHML (T, ICKWRK, RCKWRK, RCKWRK(NcK2))

      RLNP = RCKWRK(NcRU) * LOG(P / RCKWRK(NcPA))
      GBML = 0.0
      DO 100 K = 1, NKK
         GBML = GBML + X(K) * ( RCKWRK(NcK2 + K - 1) - T *   &
                (RCKWRK(NcK1 + K - 1) - RCKWRK(NcRU) *       &
                 LOG(MAX(X(K),SMALL)) - RLNP))
  100 CONTINUE

      RETURN
      END SUBROUTINE CKGBML

!     *****************************************************************

      SUBROUTINE CKGBMS (P, T, Y, ICKWRK, RCKWRK, GBMS)
!     Returns the mean Gibbs free energy of the mixture in mass units,
!     given the pressure, temperature, and mass fractions;
!     see Eq. (45).
!
!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature. [K]
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!
!     OUTPUT
!     GBMS   - Mean Gibbs free energy in mass units. [erg/g]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), Y(*)

      CALL CKSML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKHML (T, ICKWRK, RCKWRK, RCKWRK(NcK2))
      CALL CKYTX (Y, ICKWRK, RCKWRK, RCKWRK(NcK3))

      RLNP = RCKWRK(NcRU) * LOG(P / RCKWRK(NcPA))
      SUM = 0.0
      DO 100 K = 1, NKK
         SUM = SUM + RCKWRK(NcK3 + K - 1) *                      &
                   ( RCKWRK(NcK2 + K - 1) - T *                  &
                   ( RCKWRK(NcK1 + K - 1) -                      &
                     RCKWRK(NcRU) *                              &
                     LOG(MAX(RCKWRK(NcK3 + K - 1),SMALL)) - RLNP))
  100 CONTINUE

      CALL CKMMWY (Y, ICKWRK, RCKWRK, WTM)
      GBMS = SUM / WTM
      RETURN
      END SUBROUTINE CKGBMS

!     *****************************************************************

      SUBROUTINE CKGML  (T, ICKWRK, RCKWRK, GML)
!     Returns the standard state Gibbs free energies in molar units;
!     see Eq. (24).
!     INPUT
!     T      - Temperature. [K]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!
!     OUTPUT
!     GML    - Standard state gibbs free energies in molar units
!              for the species. [erg/mol]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), GML(*)

      CALL CKHML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKSML (T, ICKWRK, RCKWRK, RCKWRK(NcK2))
      DO 150 K = 1, NKK
         GML(K) = RCKWRK(NcK1 + K - 1) - T*RCKWRK(NcK2 + K - 1)
150   CONTINUE
      RETURN
      END SUBROUTINE CKGML

!     *****************************************************************

      SUBROUTINE CKGMS  (T, ICKWRK, RCKWRK, GMS)
!     Returns the standard state Gibbs free energies in mass units;
!     see Eq. (31).
!     INPUT
!     T      - Temperature. [K]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!
!     OUTPUT
!     GMS    - Standard state Gibbs free energies in mass units
!              for the species. [erg/g]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), GMS(*)
      CALL CKHMS (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKSMS (T, ICKWRK, RCKWRK, RCKWRK(NcK2))
      DO 150 K = 1, NKK
         GMS(K) = RCKWRK(NcK1 + K - 1) - T*RCKWRK(NcK2 + K - 1)
150   CONTINUE
      RETURN
      END SUBROUTINE CKGMS

!     *****************************************************************

      SUBROUTINE CKHBML (T, X, ICKWRK, RCKWRK, HBML)
!     Returns the mean enthalpy of the mixture in molar units;
!     see Eq. (37).
!     INPUT
!     T      - Temperature. [K]
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!
!     OUTPUT
!     HBML   - Mean enthalpy in molar units. [erg/mol]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), X(*)
      CALL CKHML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))

      HBML = 0.0
      DO 100 K = 1, NKK
         HBML = HBML + X(K)*RCKWRK(NcK1 + K - 1)
  100 CONTINUE
      RETURN
      END SUBROUTINE CKHBML

!     *****************************************************************

      SUBROUTINE CKHBMS (T, Y, ICKWRK, RCKWRK, HBMS)
!     Returns the mean enthalpy of the mixture in mass units;
!     see Eq. (38).
!     INPUT
!     T      - Temperature. [K]
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!
!     OUTPUT
!     HBMS   - Mean enthalpy in mass units. [erg/g]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), Y(*)

      CALL CKHMS (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      HBMS = 0.0
      DO 100 K = 1, NKK
         HBMS = HBMS + Y(K)*RCKWRK(NcK1 + K - 1)
  100 CONTINUE
      RETURN
      END SUBROUTINE CKHBMS

!     *****************************************************************

      SUBROUTINE CKHML  (T, ICKWRK, RCKWRK, HML)
!     Returns the enthalpies in molar units
!     INPUT
!     T      - Temperature. [K]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     HML    - Enthalpies in molar units for the species. [erg/mol]

      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), HML(*), TN(10)

      RUT = T*RCKWRK(NcRU)
      TN(1) = 1.0
      DO 150 N = 2, NCP
         TN(N) = T**(N-1)/N
150   CONTINUE

      DO 250 K = 1, NKK
         L = 1
         DO 220 N = 2, ICKWRK(IcNT + K - 1)-1
            TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE

         NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
         SUM = 0.0
         DO 225 N = 1, NCP
            SUM = SUM + TN(N)*RCKWRK(NA1 + N - 1)
  225    CONTINUE
         HML(K) = RUT * (SUM + RCKWRK(NA1 + NCP1 - 1)/T)
250   CONTINUE
      RETURN
      END SUBROUTINE CKHML

!     *****************************************************************

      SUBROUTINE CKHMS  (T, ICKWRK, RCKWRK, HMS)
!     Returns the enthalpies in mass units;  see Eq. (27).
!
!    INPUT
!     T      - Temperature. [K]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!
!     OUTPUT
!     HMS    - Enthalpies in mass units for the species. [erg/g]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

      DIMENSION ICKWRK(*), RCKWRK(*), HMS(*), TN(10)

      RUT = T*RCKWRK(NcRU)
      TN(1)=1.0
      DO 150 N = 2, NCP
         TN(N) = T**(N-1)/N
150   CONTINUE

      DO 250 K = 1, NKK
         L = 1
         DO 220 N = 2, ICKWRK(IcNT + K - 1)-1
            TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE

         NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
         SUM = 0.0
         DO 225 N = 1, NCP
            SUM = SUM + TN(N)*RCKWRK(NA1 + N - 1)
  225    CONTINUE
         HMS(K) = RUT * (SUM + RCKWRK(NA1 + NCP1 - 1)/T)    &
                     / RCKWRK(NcWT + K - 1)
250   CONTINUE
      RETURN
      END SUBROUTINE CKHMS

!     *****************************************************************

      SUBROUTINE CKHORT (T, ICKWRK, RCKWRK, HORT)
!     Returns the nondimensional enthalpies;  see Eq. (20).
!     INPUT
!     T      - Temperature. [K]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!
!     OUTPUT
!     HORT   - Nondimensional enthalpies for the species.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION TN(10), HORT(*), ICKWRK(*), RCKWRK(*)
      DO 150 N = 1, NCP
         TN(N) = T**(N-1)/N
150   CONTINUE

      DO 250 K = 1, NKK
         L = 1
         DO 220 N = 2, ICKWRK(IcNT + K - 1)-1
            TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE

         NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
         SUM = 0.0
         DO 225 N = 1, NCP
            SUM = SUM + TN(N)*RCKWRK(NA1 + N - 1)
  225    CONTINUE
         HORT(K) = SUM + RCKWRK(NA1 + NCP1 - 1)/T
250   CONTINUE
      RETURN
      END SUBROUTINE CKHORT

!     *****************************************************************

      SUBROUTINE CKI2CH (NUM, STR, I, KERR)
!     Returns a character string representation of an integer
!     and the effective length of the string.
!
!     INPUT
!     NUM   - A number to be converted to a character string;
!             the maximum magnitude of NUM is machine-dependent.
!                  Data type - integer scalar.
!
!     OUTPUT
!     STR   - A left-justified character string representing NUM
!                  Data type - CHARACTER*(*)
!     I     - The effective length of the character string
!                  Data type - integer scalar
!     KERR  - Error flag;  character length errors will result in
!             KERR=.TRUE.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

      CHARACTER STR*(*), IST(10)*(1)
      LOGICAL KERR
      DATA IST/'0','1','2','3','4','5','6','7','8','9'/
      BIGI = 2147483647.

      I = 0
      STR = ' '
      ILEN = LEN(STR)
      KERR = .FALSE.
      IF (ILEN.LT.1 .OR. ABS(NUM).GT.BIGI) THEN
         WRITE(*,*) 'ERROR IN INTEGER TO STRING CONVERSION'
         WRITE(*,*) 'ILEN = ',ILEN
         WRITE(*,*) 'NUM  = ',abs(NUM)
         KERR = .TRUE.
         RETURN
      ENDIF

      IF (NUM .EQ. 0) THEN
         STR = '0'
         I = 1
         RETURN
      ELSEIF (NUM .LT. 0) THEN
         STR(1:) = '-'
      ENDIF

      INUM = ABS(NUM)
      NCOL = NINT(LOG10(REAL(INUM))) + 1

      DO 10 J = NCOL, 1, -1
         IDIV = INUM / 10.0**(J-1)
         IF (J.EQ.NCOL .AND. IDIV.EQ.0) GO TO 10
         LT = ILASCH(STR)
         IF (LT .EQ. ILEN) THEN
            STR = ' '
            KERR = .TRUE.
            RETURN
         ENDIF
         STR(LT+1:) = IST(IDIV+1)
         INUM = INUM - IDIV*10.0**(J-1)
   10 CONTINUE
      I = ILASCH(STR)

      RETURN
      END SUBROUTINE CKI2CH


!     *****************************************************************

      SUBROUTINE CKINDX (ICKWRK, RCKWRK, MM, KK, II, NFIT)
!     Returns a group of indices defining the size of the particular
!     reaction mechanism
!
!     INPUT
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!
!     OUTPUT
!     MM     - Total number of elements in mechanism.
!     KK     - Total number of species in mechanism.
!     II     - Total number of reactions in mechanism.
!     NFIT   - number of coefficients in fits to thermodynamic data
!              for one temperature range; NFIT = number of
!              coefficients in polynomial fits to CP/R  +  2.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*)

      MM = NMM
      KK = NKK
      II = NII
      NFIT = NCP2
      RETURN
      END SUBROUTINE CKINDX

!     *****************************************************************

      SUBROUTINE CKINIT (LENIWK, LENRWK, LENCWK, LINC, LOUT, ICKWRK,  &
                         RCKWRK, CCKWRK)
!     Reads the linking file and creates the internal work arrays
!     ICKWRK, CCKWRK, and RCKWORK.  CKINIT must be called before any
!     other CHEMKIN subroutine is called.  The work arrays must then
!     be made available as input to the other CHEMKIN subroutines.

!     INPUT
!     LENIWK - Length of the integer work array, ICKWRK.
!                   Data type - integer scalar
!     LENCWK - Length of the character work array, CCKWRK.
!              The minimum length of CCKWRK(*) is MM + KK.
!                   Data type - integer scalar
!     LENRWK - Length of the real work array, WORK.
!                   Data type - integer scalar
!     LINC  -  Logical file number for the linking file.
!                   Data type - integer scalar
!     LOUT  -  Output file for printed diagnostics.
!                   Data type - integer scalar
!
!     OUTPUT
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!     CCKWRK - Array of character work space.
!                   Data type - CHARACTER*18 array
!                   Dimension CCKWRK(*) at least LENCWK.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

      DIMENSION ICKWRK(*), RCKWRK(*)
      CHARACTER CCKWRK(*)*(*)
      LOGICAL IOK, ROK, COK

      DATA RU,RUC,PA /8.314E7, 1.987, 1.01325E6/

      SMALL = 10.0D0**(-300)
      BIG   = 10.0D0**(+300)
      EXPARG = LOG(BIG)

      WRITE (LOUT,15)
   15 FORMAT (/1X,' CKLIB:  Chemical Kinetics Library',              &
              /1X,'         CHEMKIN-II Version 4.2, September 1993', &
              /1X,'         DOUBLE PRECISION')

      CALL CKLEN (LINC, LOUT, LI, LR, LC)

      IOK = (LENIWK .GE. LI)
      ROK = (LENRWK .GE. LR)
      COK = (LENCWK .GE. LC)
      IF (.NOT.IOK .OR. .NOT.ROK .OR. .NOT.COK) THEN
         IF (.NOT. IOK) WRITE (LOUT, 300) LI
         IF (.NOT. ROK) WRITE (LOUT, 350) LR
         IF (.NOT. COK) WRITE (LOUT, 375) LC
         stop!(1)
      ENDIF

      REWIND LINC
      READ (LINC, ERR=110) VERS, PREC, KERR

      READ (LINC, ERR=110) LENI, LENR, LENC, MM, KK, II,              &
                           MAXSP, MAXTB, MAXTP, NTHCF, NIPAR, NITAR,  &
                           NIFAR, NRV, NFL, NTB, NLT, NRL, NW, NCHRG

      IF (LEN(CCKWRK(1)) .LT. species_name_len) THEN
         WRITE (LOUT,475)
         stop
      ENDIF

      NMM = MM
      NKK = KK
      NII = II
      MXSP = MAXSP
      MXTB = MAXTB
      MXTP = MAXTP
      NCP  = NTHCF
      NCP1 = NTHCF+1
      NCP2 = NTHCF+2
      NCP2T = NCP2*(MAXTP-1)
      NPAR = NIPAR
      NLAR = NITAR
      NFAR = NIFAR
      NTHB = NTB
      NLAN = NLT
      NFAL = NFL
      NREV = NRV
      NRLT = NRL
      NWL  = NW

!             APPORTION work arrays

!            SET  ICKWRK(*)=1  TO FLAG THAT CKINIT HAS BEEN CALLED

      ICKWRK(1) = 1

!             STARTING LOCATIONS OF INTEGER SPACE

! elemental composition of species
      IcNC = 2
! species phase array
      IcPH = IcNC + KK*MM
! species charge array
      IcCH = IcPH + KK
! # of temperatures for fit
      IcNT = IcCH + KK
! stoichiometric coefficients
      IcNU = IcNT + KK
! species numbers for the coefficients
      IcNK = IcNU + MAXSP*II
! # of non-zero coefficients  (<0=reversible, >0=irreversible)
      IcNS = IcNK + MAXSP*II
! # of reactants
      IcNR = IcNS + II
! Landau-Teller reaction numbers
      IcLT = IcNR + II
! Reverse Landau-Teller reactions
      IcRL = IcLT + NLAN
! Fall-off reaction numbers
      IcFL = IcRL + NRLT
! Fall-off option numbers
      IcFO = IcFL + NFAL
! Fall-off enhanced species
      IcKF = IcFO + NFAL
! Third-body reaction numbers
      IcTB = IcKF + NFAL
! number of 3rd bodies for above
      IcKN = IcTB + NTHB
! array of species #'s for above
      IcKT = IcKN + NTHB
! Reverse parameter reaction numbers
      IcRV = IcKT + MAXTB*NTHB
! Radiation wavelength reactions
      IcWL = IcRV + NREV
      ITOT = IcWL + NWL - 1

!             STARTING LOCATIONS OF CHARACTER SPACE

! start of element names
      IcMM = 1
! start of species names
      IcKK = IcMM + MM
      ITOC = IcKK + KK - 1

!             STARTING LOCATIONS OF REAL SPACE

! atomic weights
      NcAW = 1
! molecular weights
      NcWT = NcAW + MM
! temperature fit array for species
      NcTT = NcWT + KK
! thermodynamic coefficients
      NcAA = NcTT + MAXTP*KK
! Arrhenius coefficients (3)
      NcCO = NcAA + (MAXTP-1)*NCP2*KK
! Reverse coefficients
      NcRV = NcCO + (NPAR+1)*II
! Landau-Teller #'s for NLT reactions
      NcLT = NcRV + (NPAR+1)*NREV
! Reverse Landau-Teller #'s
      NcRL = NcLT + NLAR*NLAN
! Fall-off parameters for NFL reactions
      NcFL = NcRL + NLAR*NRLT
! 3rd body coef'nts for NTHB reactions
      NcKT = NcFL + NFAR*NFAL
! wavelength
      NcWL = NcKT + MAXTB*NTHB
! universal gas constant
      NcRU = NcWL + NWL
! universal gas constant in units
      NcRC = NcRU + 1
! pressure of one atmosphere
      NcPA = NcRC + 1
! intermediate temperature-dependent forward rates
      NcKF = NcPA + 1
! intermediate temperature-dependent reverse rates
      NcKR = NcKF + II
! internal work space of length kk
      NcK1 = NcKR + II
!          'ditto'
      NcK2 = NcK1 + KK
!          'ditto'
      NcK3 = NcK2 + KK
!          'ditto'
      NcK4 = NcK3 + KK
      NcI1 = NcK4 + KK
      NcI2 = NcI1 + II
      NcI3 = NcI2 + II
      NcI4 = NcI3 + II
      NTOT = NcI4 + II - 1

!        SET UNIVERSAL CONSTANTS IN CGS UNITS

      RCKWRK(NcRU) = RU
      RCKWRK(NcRC) = RUC
      RCKWRK(NcPA) = PA

!element names, !atomic weights
      READ (LINC,err=111) (CCKWRK(IcMM+M-1), RCKWRK(NcAW+M-1), M=1,MM)

!species names, !composition, !phase, !charge, !molec weight,
!# of fit temps, !array of temps, !fit coeff'nts
      READ (LINC,err=222) (CCKWRK(IcKK+K-1),                   &
           (ICKWRK(IcNC+(K-1)*MM + M-1),M=1,MM),               &
           ICKWRK(IcPH+K-1),                                   &
           ICKWRK(IcCH+K-1),                                   &
           RCKWRK(NcWT+K-1),                                   &
           ICKWRK(IcNT+K-1),                                   &
           (RCKWRK(NcTT+(K-1)*MAXTP + L-1),L=1,MAXTP),         &
           ((RCKWRK(NcAA+(L-1)*NCP2+(K-1)*NCP2T+N-1),          &
           N=1,NCP2), L=1,(MAXTP-1)),    K = 1,KK)

      IF (II .EQ. 0) RETURN

!# spec,reactants, !Arr. coefficients, !stoic coef, !species numbers
      READ (LINC,end=100,err=333)                              &
           (ICKWRK(IcNS+I-1), ICKWRK(IcNR+I-1),                &
            (RCKWRK(NcCO+(I-1)*(NPAR+1)+N-1), N=1,NPAR),       &
            (ICKWRK(IcNU+(I-1)*MAXSP+N-1),                     &
             ICKWRK(IcNK+(I-1)*MAXSP+N-1), N=1,MAXSP),         &
            I = 1,II)

!     PERTURBATION FACTOR

      DO 10 I = 1, II
         RCKWRK(NcCO + (I-1)*(NPAR+1) + NPAR) = 1.0
   10 CONTINUE

      IF (NREV .GT. 0) READ (LINC,err=444)                            &
         (ICKWRK(IcRV+N-1), (RCKWRK(NcRV+(N-1)*(NPAR+1)+L-1),         &
         L=1,NPAR), N = 1,NREV)

      IF (NFAL .GT. 0) READ (LINC,err=555)                            &
         (ICKWRK(IcFL+N-1), ICKWRK(IcFO+N-1), ICKWRK(IcKF+N-1),       &
         (RCKWRK(NcFL+(N-1)*NFAR+L-1),L=1,NFAR),N=1,NFAL)

      IF (NTHB .GT. 0) READ (LINC,err=666)                            &
         (ICKWRK(IcTB+N-1), ICKWRK(IcKN+N-1),                         &
         (ICKWRK(IcKT+(N-1)*MAXTB+L-1),                               &
           RCKWRK(NcKT+(N-1)*MAXTB+L-1),L=1,MAXTB),N=1,NTHB)

      IF (NLAN .GT. 0) READ (LINC,err=777)                            &
         (ICKWRK(IcLT+N-1), (RCKWRK(NcLT+(N-1)*NLAR+L-1),L=1,NLAR),   &
          N=1,NLAN)

      IF (NRLT .GT. 0) READ (LINC,err=888)                            &
         (ICKWRK(IcRL+N-1), (RCKWRK(NcRL+(N-1)*NLAR+L-1),L=1,NLAR),   &
          N=1,NRLT)

      IF (NWL .GT. 0) READ (LINC,err=999)                             &
         (ICKWRK(IcWL+N-1), RCKWRK(NcWL+N-1), N=1,NWL)

  100 CONTINUE
      RETURN

  110 WRITE (LOUT,*) ' Error reading linking file...'
      stop
  111 WRITE (LOUT,*) ' Error reading element data...'
      stop
  222 WRITE (LOUT,*) ' Error reading species data...'
      stop
  333 WRITE (LOUT,*) ' Error reading reaction data...'
      stop
  444 WRITE (LOUT,*) ' Error reading reverse Arrhenius parameters...'
      stop
  555 WRITE (LOUT,*) ' Error reading Fall-off data...'
      stop
  666 WRITE (LOUT,*) ' Error reading third-body data...'
      stop
  777 WRITE (LOUT,*) ' Error reading Landau-Teller data...'
      stop
  888 WRITE (LOUT,*) ' Error reading reverse Landau-Teller data...'
      stop
  999 WRITE (LOUT,*) ' Error reading Wavelength data...'
      stop

  300 FORMAT (10X,'ICKWRK MUST BE DIMENSIONED AT LEAST ',I5)
  350 FORMAT (10X,'RCKWRK MUST BE DIMENSIONED AT LEAST ',I5)
  375 FORMAT (10X,'CCKWRK MUST BE DIMENSIONED AT LEAST ',I5)
  475 FORMAT (10X,'CHARACTER LENGTH OF CCKWRK MUST BE AT LEAST 18 ')
      END SUBROUTINE CKINIT


!     *****************************************************************

      SUBROUTINE CKITR  (ICKWRK, RCKWRK, ITHB, IREV)
!     Returns a set of flags indicating whether the reactions are
!     reversible or whether they contain arbitrary third bodies
!     INPUT
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     ITHB   - Third-body flags for the reactions;
!              ITHB(I)= -1  reaction I is not a third-body reactions
!              ITHB(I)=  0  reaction I is is a third-body reaction with
!                           no enhanced third body efficiencies
!              ITHB(I)=  N  reaction I is a third-body reaction with
!                           N species enhanced third-body efficiencies.
!
!     IREV   - Reversibility flags and number of species
!              (reactants plus products) for reactions.
!              IREV(I)=+N, reversible reaction I has N species
!              IREV(I)=-N, irreversible reaction I has N species
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ITHB(*), IREV(*), ICKWRK(*), RCKWRK(*)

      DO 100 I = 1, NII
         IREV(I) = ICKWRK(IcNS + I - 1)
         ITHB(I) = -1
  100 CONTINUE
      DO 150 N = 1, NTHB
         ITHB(ICKWRK(IcTB + N - 1)) = ICKWRK(IcKN + N - 1)
  150 CONTINUE

      RETURN
      END SUBROUTINE CKITR

!     *****************************************************************

      SUBROUTINE CKKFKR (P, T, X, ICKWRK, RCKWRK, FWDK, REVK)
!     Returns the forward and reverse reaction rates for the
!     reactions given pressure, temperature and mole fractions.
!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature. [K]
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     FWDK   - Forward reaction rates for the reactions.
!                   cgs units - depends on the reaction
!     REVK   - Reverse reaction rates for the reactions.
!                   cgs units - depends on the reaction
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), X(*), FWDK(*), REVK(*)

      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),              &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),          &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,             &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT), &
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),       &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),             &
                   RCKWRK(NcI1))

      CALL CKXTCP (P, T, X, ICKWRK, RCKWRK, RCKWRK(NcK1))

      CALL CKRATX (NII, NKK, MXSP, MXTB, T, RCKWRK(NcK1), ICKWRK(IcNS),  &
                   ICKWRK(IcNU), ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO),     &
                   NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF), NFAR, &
                   RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),       &
                   RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF),             &
                   RCKWRK(NcKR), RCKWRK(NcI1), RCKWRK(NcI2),             &
                   RCKWRK(NcI3))

      DO 200 I = 1, NII
         FWDK(I) = RCKWRK(NcI1 + I - 1)
         REVK(I) = RCKWRK(NcI2 + I - 1)
  200 CONTINUE
      RETURN
      END SUBROUTINE CKKFKR

!     *****************************************************************

      SUBROUTINE CKKFRT (P, T, ICKWRK, RCKWRK, RKFT, RKRT)
!     Returns the forward and reverse reaction rates for the
!     reactions given pressure and temperature.

!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature. [K]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!
!     OUTPUT
!     RKFT   - Forward reaction rates for the reactions.
!                   cgs units - depends on the reaction
!     RKRT   - Reverse reaction rates for the reactions.
!                   cgs units - depends on the reaction

      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), RKFT(*), RKRT(*)

      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),             &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),         &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,            &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT),&
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),      &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),            &
                   RCKWRK(NcI1))

      DO 200 I = 1, NII
         RKFT(I) = RCKWRK(NcKF + I - 1)
         RKRT(I) = RCKWRK(NcKR + I - 1)
  200 CONTINUE

      RETURN
      END SUBROUTINE CKKFRT

!     *****************************************************************

      SUBROUTINE CKLEN (LINC, LOUT, LI, LR, LC)
!     Returns the lengths required for the work arrays.
!     INPUT
!     LINC  -  Logical file number for the linking file.
!     LOUT  -  Output file for printed diagnostics.
!
!     OUTPUT
!     LENI  -  Minimum length required for the integer work array.
!     LENR  -  Minimum length required for the real work array.
!     LENC  -  Minimum length required for the character work array.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

      PARAMETER (NLIST = 14)
      LOGICAL VOK, POK, KERR
      CHARACTER LIST(NLIST)*16

      DATA LIST/'1.9','2.0','2.1','2.2','2.3','2.4','2.5','2.6', &
                '2.7','2.8','2.9','3.0','3.1','3.9'/

      VERS = ' '
      PREC = ' '
      LENI = 0
      LENR = 0
      LENC = 0

      KERR = .FALSE.
      REWIND LINC
      READ (LINC, ERR=999) VERS, PREC, KERR
      VOK = .FALSE.
      DO 5 N = 1, NLIST
         IF (VERS .EQ. LIST(N)) VOK = .TRUE.
    5 CONTINUE

      POK = .FALSE.

      IF (INDEX(PREC, 'DOUB') .GT. 0) POK = .TRUE.

      IF (KERR .OR. (.NOT.POK) .OR. (.NOT.VOK)) THEN
         IF (KERR) THEN
            WRITE (LOUT,'(/A,/A)')                                    &
            ' There is an error in the Chemkin linking file...',      &
            ' Check CHEMKIN INTERPRETER output for error conditions.'
            WRITE(LOUT,*) ' KERR = ',KERR
            WRITE(LOUT,*) ' POK  = ',POK
            WRITE(LOUT,*) ' VOK  = ',VOK
         ENDIF
         IF (.NOT. VOK) THEN
            WRITE (LOUT,'(/A,A)')                                     &
            ' Chemkin linking file is incompatible with Chemkin',     &
            ' Library Version 4.2'
         ENDIF
         IF (.NOT. POK) THEN
            WRITE (LOUT,'(/A,A)')                                     &
            ' Precision of Chemkin linking file does not agree with', &
            ' precision of Chemkin library'
         ENDIF
         stop
      ENDIF

      READ (LINC, ERR=999) LENICK, LENRCK, LENCCK, MM, KK, II,        &
                           MAXSP, MAXTB, MAXTP, NTHCF, NIPAR, NITAR,  &
                           NIFAR, NRV, NFL, NTB, NLT, NRL, NW, NCHRG
      REWIND LINC

      LENI = LENICK
      LENR = LENRCK
      LENC = LENCCK
      LI   = LENI
      LR   = LENR
      LC   = LENC
      RETURN

  999 CONTINUE
      WRITE (LOUT, 50)
   50 FORMAT (' Error reading Chemkin Linking file.')
      stop
      END SUBROUTINE CKLEN

!     *****************************************************************

      SUBROUTINE CKMMWC (C, ICKWRK, RCKWRK, WTM)
!     Returns the mean molecular weight of the gas mixture given the
!     molar concentrations;  see Eq. (5).
!     INPUT
!     C      - Molar concentrations of the species. [mol/cm3]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     OUTPUT
!     WTM    - Mean molecular weight of the species mixture. [g/mol]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION C(*), ICKWRK(*), RCKWRK(*)

      CTOT = 0.0
      DO 100 K = 1, NKK
         CTOT = CTOT + C(K)
  100 CONTINUE

      WTM = 0.0
      DO 200 K = 1, NKK
         WTM = WTM + C(K)*RCKWRK(NcWT + K - 1)
  200 CONTINUE
      WTM = WTM / CTOT

      RETURN
      END SUBROUTINE CKMMWC

!     *****************************************************************


      SUBROUTINE CKMMWX (X, ICKWRK, RCKWRK, WTM)
!     Returns the mean molecular weight of the gas mixture given the
!     mole fractions;  see Eq. (4).
!     INPUT
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     WTM    - Mean molecular weight of the species mixture. [g/mol]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION X(*), ICKWRK(*), RCKWRK(*)
      WTM = 0.0
      DO 100 K = 1, NKK
         WTM = WTM + X(K)*RCKWRK(NcWT + K - 1)
  100 CONTINUE
      RETURN
      END SUBROUTINE CKMMWX

!     *****************************************************************

      SUBROUTINE CKMMWY (Y, ICKWRK, RCKWRK, WTM)
!     Returns the mean molecular weight of the gas mixture given the
!     mass fractions;  see Eq. (3).
!     INPUT
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     WTM    - Mean molecular weight of the species mixture. [g/mol]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION Y(*), ICKWRK(*), RCKWRK(*)
      SUMYOW=0.0
      DO 150 K = 1, NKK
         SUMYOW = SUMYOW + Y(K)/RCKWRK(NcWT + K - 1)
150   CONTINUE
      WTM = 1.0/SUMYOW
      RETURN
      END SUBROUTINE CKMMWY

!     *****************************************************************


      SUBROUTINE CKMXTP (ICKWRK, MAXTP)
!     Returns the maximum number of temperatures used in
!     fitting the thermodynamic properties of the species.
!     INPUT
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!     OUTPUT
!     MXTP   - Maximum number of temperatures used in
!              fitting the thermodynamic properties of
!              the species.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*)

      MAXTP = MXTP
      RETURN
      END SUBROUTINE CKMXTP

!     *****************************************************************


      SUBROUTINE CKNCF  (MDIM, ICKWRK, RCKWRK, NCF)
!     Returns the elemental composition of the species
!     INPUT
!     MDIM   - First dimension of the two-dimensional array NCF;
!              MDIM must be equal to or greater than the number of
!              elements, MM.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     NCF    - Matrix of the elemental composition of the species;
!              NCF(M,K) is the number of atoms of the Mth element
!              in the Kth species.
!                   Data type - integer array
!                   Dimension NCF(MDIM,*) exactly MDIM (at least MM,
!                   the total number of elements in the problem) for
!                   the first dimension and at least KK, the total
!                   number of species, for the second.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), NCF(MDIM,*)
      DO 150 K = 1, NKK
         J = IcNC + (K-1)*NMM
         DO 150 M = 1, NMM
            NCF(M,K) = ICKWRK(J + M - 1)
150   CONTINUE
      RETURN
      END SUBROUTINE CKNCF

!     *****************************************************************


      SUBROUTINE CKNPAR (LINE, NPAR, LOUT, IPAR, ISTART, KERR)
!     This subroutine is called to parse a character string, LINE,
!     that is composed of several blank-delimited substrings.
!     That final segment of LINE containing NPAR substrings is
!     found, beginning in the ISTART column; this segment is
!     then copied into the character string IPAR.  This allows
!     format-free input of combined alpha-numeric data.
!     For example, after reading a line containing alpha-numeric
!     information ending with several numbers, the subroutine
!     might be called to find the segment of the line containing
!     the numbers:

!     input:  LINE*80   = "t1 t2 dt  300.0  3.0E3  50"
!             NPAR      = 3, the number of substrings requested
!             LOUT      = 6, a logical unit number on which to write
!                         diagnostic messages.
!     output: IPAR*80   = "300.0  3.0E3  50"
!             ISTART    = 13, the starting column in LINE of the
!                         NPAR substrings
!             KERR      = .FALSE.
!     INPUT
!     LINE   - A character string.
!                   Data type - CHARACTER*(*)
!     NPAR   - Number of substrings expected.
!                   Data type - integer scalar
!     LOUT   - Output unit for printed diagnostics.
!                   Data type - integer scalar

!     OUTPUT
!     IPAR   - A character string containing only the NPAR substrings.
!                   Data type - CHARACTER*(*)
!     ISTART - The starting location in LINE of the NPAR substrings.
!                   Data type - integer scalar
!     KERR   - Error flag; character length or syntax error will
!              result in KERR = .TRUE.
!                   Date type: logical

!     A '!' will comment out a line, or remainder of the line.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      CHARACTER LINE*(*), IPAR*(*)
      LOGICAL FOUND, KERR

!----------Find Comment String (! signifies comment)

      ILEN = IPPLEN(LINE)
      KERR = .FALSE.

      IF (ILEN.GT.0) THEN
         FOUND = .FALSE.
         N = 0
         DO 40 I = ILEN, 1, -1
            IF (FOUND) THEN
               IF (LINE(I:I).EQ.' ') THEN
                  N = N+1
                  FOUND = .FALSE.
                  IF (N.EQ.NPAR) THEN
                     ISTART = I+1
                     L1 = ILEN - ISTART + 1
                     L2 = LEN(IPAR)
                     IF (L2 .GE. L1) THEN
                        IPAR = LINE(ISTART:ILEN)
                     ELSE
                        WRITE (LOUT,*) &
                     ' Error in CKNPAR...character length too small...'
                        KERR = .TRUE.
                     ENDIF
                     RETURN
                  ENDIF
               ENDIF
            ELSE
               IF (LINE(I:I).NE.' ') FOUND = .TRUE.
            ENDIF
   40    CONTINUE
      ENDIF

      WRITE (LOUT,*) ' Error in CKNPAR...',NPAR,' values not found...'
      KERR = .TRUE.
      RETURN
      END SUBROUTINE CKNPAR

!     *****************************************************************

      SUBROUTINE CKNU   (KDIM, ICKWRK, RCKWRK, NUKI)

!     Returns the stoichiometric coefficients of the reaction
!     mechanism;  see Eq. (50).
!     INPUT
!     KDIM   - First dimension of the two-dimensional array NUKI;
!              KDIM must be greater than or equal to the total
!              number of species, KK.
!                   Data type - integer scalar
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!     OUTPUT
!     NUKI   - Matrix of stoichiometric coefficients for the species
!              in the reactions;  NUKI(K,I) is the stoichiometric
!              coefficient of species K in reaction I.
!                   Data type - integer array
!                   Dimension NUKI(KDIM,*) exactly KDIM (at least KK,
!                   the total number of species) for the first
!                   dimension and at least II for the second, the total
!                   number of reactions.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), NUKI(KDIM,*)

      DO 100 I = 1, NII
         DO 100 K = 1, NKK
            NUKI(K,I) = 0
  100 CONTINUE
      DO 200 N = 1, MXSP
         DO 200 I = 1, NII
            K = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
            NU= ICKWRK(IcNU + (I-1)*MXSP + N -1)
            IF (K .NE. 0) NUKI(K,I) = NUKI(K,I) + NU
200   CONTINUE
      RETURN
      END SUBROUTINE CKNU

!     *****************************************************************

      SUBROUTINE CKNUF   (KDIM, ICKWRK, RCKWRK, NUFKI)
!     Returns the stoichiometric coefficients for the forward
!     reactions in the reaction mechanism.  All stoichiometric
!     coefficients for reactants are defined to be negative, by
!     definition; see Eq. (50).  Note this subroutine is to be
!     contrasted with the subroutine, CKNU, which returns the net
!     stoichiometric coefficients for a reaction.
!     INPUT
!     KDIM   - First dimension of the two-dimensional array NUKI;
!              KDIM must be greater than or equal to the total
!              number of species, KK.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     NUFKI  - Matrix of stoichiometric coefficients for the species
!              in the forward directions of the reactions;
!              NUKI(K,I) is the stoichiometric
!              coefficient of species K in forward direction of
!              reaction I.
!                   Data type - integer array
!                   Dimension NUKI(KDIM,*) exactly KDIM (at least KK,
!                   the total number of species) for the first
!                   dimension and at least II for the second, the total
!                   number of reactions.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), NUFKI(KDIM,*)

      DO 100 I = 1, NII
         DO 100 K = 1, NKK
            NUFKI(K,I) = 0
  100 CONTINUE
      IF (MXSP .EQ. 6) THEN
        DO 200 I = 1, NII
	   K1 = ICKWRK(IcNK + (I-1)*MXSP )
	   NU1= ICKWRK(IcNU + (I-1)*MXSP )
	   IF (K1 .NE. 0) NUFKI(K1,I) = NUFKI(K1,I) + NU1
	   K2 = ICKWRK(IcNK + (I-1)*MXSP + 1)
	   NU2= ICKWRK(IcNU + (I-1)*MXSP + 1)
	   IF (K2 .NE. 0) NUFKI(K2,I) = NUFKI(K2,I) + NU2
	   K3 = ICKWRK(IcNK + (I-1)*MXSP + 2)
	   NU3= ICKWRK(IcNU + (I-1)*MXSP + 2)
	   IF (K3 .NE. 0) NUFKI(K3,I) = NUFKI(K3,I) + NU3
 200    CONTINUE
      ELSE
        DO 300 N = 1, (MXSP/2)
        DO 300 I = 1, NII
	   K = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
	   NU= ICKWRK(IcNU + (I-1)*MXSP + N - 1)
	   IF (K .NE. 0) NUFKI(K,I) = NUFKI(K,I) + NU
 300    CONTINUE
      ENDIF
      RETURN
      END SUBROUTINE CKNUF

!     *****************************************************************

      SUBROUTINE CKPC   (RHO, T, C, ICKWRK, RCKWRK, P)
!     Returns the pressure of the gas mixture given the mass density,
!     temperature and molar concentrations;  see Eq. (2).

!     INPUT
!     RHO    - Mass density. [g/cm3]
!     T      - Temperature. [K]
!     C      - Molar concentrations of the species. [mol/cm3]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     P      - Pressure. [dyn/cm2]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION C(*), ICKWRK(*), RCKWRK(*)

      CTOT = 0.0
      SUM = 0.0
      DO 100 K = 1, NKK
         CTOT = CTOT + C(K)
         SUM  = SUM + C(K)*RCKWRK(NcWT + K - 1)
  100 CONTINUE
      P    = RHO*RCKWRK(NcRU) * T * CTOT / SUM
      RETURN
      END SUBROUTINE CKPC

!     *****************************************************************

      SUBROUTINE CKPHAZ (ICKWRK, RCKWRK, KPHASE)
!     Returns a set of flags indicating phases of the species

!     INPUT
!     ICKWRK - Array of integer workspace.

!     RCKWRK - Array of real work space.

!     OUTPUT
!     KPHASE - Phases of the species;
!              KPHASE(K)=-1  the Kth species is solid
!              KPHASE(K)= 0  the Kth species is gaseous
!              KPHASE(K)=+1  the Kth species is liquid
!                   Data type - integer array
!                   Dimension KPHASE(*) at least KK, the total number of
!                   species.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), KPHASE(*)
      DO 100 K = 1, NKK
         KPHASE(K) = ICKWRK(IcPH + K - 1)
  100 CONTINUE
      RETURN
      END SUBROUTINE CKPHAZ

!     *****************************************************************

      SUBROUTINE CKPNT (LSAVE, LOUT, NPOINT, V, P, LI, LR, LC, IERR)
!     Reads from a binary file information about a Chemkin
!     linking file, pointers for the Chemkin Library, and
!     returns lengths of work arrays.

!     INPUT
!     LSAVE  - Integer input unit for binary data file.
!                   Data type - integer scalar
!     LOUT   - Integer output unit for printed diagnostics.
!                   Data type - integer scalar
!
!     OUTPUT
!     NPOINT - Total number of pointers.
!                   Data type - integer scalar
!     VERS   - Version number of the Chemkin linking file.
!                   Data type - real scalar
!     PREC   - Machine precision of the Chemkin linking file.
!                   Data type - character string
!     LENI   - Minimum length required for the integer work array.
!                   Data type - integer scalar
!     LENR   - Minimum length required for the real work array.
!                   Data type - integer scalar
!     LENC   - Minimum length required for the character work array.
!                   Data type - integer scalar
!     KERR   - Logical error flag.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      LOGICAL IERR
      CHARACTER P*16, V*16

!*****exponent range > +/-300
      SMALL = 10.0D0**(-300)
      BIG   = 10.0D0**(+300)
!*****END exponent range > +/-300
      EXPARG = LOG(BIG)

      KERR = .FALSE.
      READ (LSAVE, ERR=100) NPOINT, VERS, PREC, LENI, LENR, LENC,      &
                      NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,  &
                      NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,  &
                      NTHB, NRLT, NWL,  IcMM, IcKK, IcNC, IcPH, IcCH,  &
                      IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL, IcRV,  &
                      IcWL, IcFL, IcFO, IcKF, IcTB, IcKN, IcKT, NcAW,  &
                      NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL,  &
                      NcKT, NcWL, NcRU, NcRC, NcPA, NcKF, NcKR, NcK1,  &
                      NcK2, NcK3, NcK4, NcI1, NcI2, NcI3, NcI4
      V = VERS
      P = PREC
      LI = LENI
      LR = LENR
      LC = LENC
      IERR = KERR
      RETURN

  100 CONTINUE
      WRITE (LOUT, *) ' Error reading Chemkin linking file data...'
      KERR   = .TRUE.
      IERR   = KERR
      NPOINT = 0
      VERS   = ' '
      V      = VERS
      PREC   = ' '
      P      = PREC
      RETURN
      END SUBROUTINE CKPNT

!     *****************************************************************

      SUBROUTINE CKPX   (RHO, T, X, ICKWRK, RCKWRK, P)
!     Returns the pressure of the gas mixture given the mass density,
!     temperature and mole fractions;  see Eq. (*).
!     INPUT
!     RHO    - Mass density. [g/cm3]
!     T      - Temperature. [K]
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     P      - Pressure. [dyn/cm2]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION X(*), ICKWRK(*), RCKWRK(*)

      SUM = 0.0
      DO 100 K = 1, NKK
         SUM = SUM + X(K)*RCKWRK(NcWT + K - 1)
  100 CONTINUE
      P = RHO * RCKWRK(NcRU) * T / SUM

      RETURN
      END SUBROUTINE CKPX


!     *****************************************************************

      SUBROUTINE CKPY   (RHO, T, Y, ICKWRK, RCKWRK, P)
!     Returns the pressure of the gas mixture given the mass density,
!     temperature and mass fractions;  see Eq. (*).
!     INPUT
!     RHO    - Mass density. [g/cm3]
!     T      - Temperature. [K]
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     P      - Pressure. [dyn/cm2]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION Y(*), ICKWRK(*), RCKWRK(*)
      SUMYOW = 0.0
      DO 150 K = 1, NKK
         SUMYOW = SUMYOW + Y(K)/RCKWRK(NcWT + K - 1)
150   CONTINUE
      P = RHO * RCKWRK(NcRU) * T * SUMYOW
      RETURN
      END SUBROUTINE CKPY

!     *****************************************************************

      SUBROUTINE CKQC   (T, C, ICKWRK, RCKWRK, Q)
!     Returns the rates of progress for the reactions given
!     temperature and molar concentrations;  see Eqs. (51) and (58).
!
!     INPUT
!     T      - Temperature. [K]
!     C      - Molar concentrations of the species. [mol/cm3]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!
!     OUTPUT
!     Q      - Rates of progress for the reactions. [mol/cm3/s]

      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION C(*), ICKWRK(*), RCKWRK(*), Q(*)

      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),              &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),          &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,             &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT), &
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),       &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),             &
                   RCKWRK(NcI1))

      DO 50 K = 1, NKK
         RCKWRK(NcK1 + K - 1) = C(K)
   50 CONTINUE

      CALL CKRATX (NII, NKK, MXSP, MXTB, T, RCKWRK(NcK1), ICKWRK(IcNS),  &
                   ICKWRK(IcNU), ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO),     &
                   NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF), NFAR, &
                   RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),       &
                   RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF),             &
                   RCKWRK(NcKR), RCKWRK(NcI1), RCKWRK(NcI2),             &
                   RCKWRK(NcI3))

      DO 100 I = 1, NII
         Q(I) = RCKWRK(NcI1 + I - 1) - RCKWRK(NcI2 + I - 1)
  100 CONTINUE
      RETURN
      END SUBROUTINE CKQC

!     *****************************************************************

      SUBROUTINE CKQXP  (P, T, X, ICKWRK, RCKWRK, Q)
!     Returns the rates of progress for the reactions given pressure,
!     temperature and mole fractions;  see Eqs. (51) and (58).
!
!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature. [K]
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!
!     OUTPUT
!     Q      - Rates of progress for the reactions. [mol/cm3/s]

      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION X(*), ICKWRK(*), RCKWRK(*), Q(*)

      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),              &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),          &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,             &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT), &
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),       &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),             &
                   RCKWRK(NcI1))

      CALL CKXTCP (P, T, X, ICKWRK, RCKWRK, RCKWRK(NcK1))

      CALL CKRATX (NII, NKK, MXSP, MXTB, T, RCKWRK(NcK1), ICKWRK(IcNS),  &
                   ICKWRK(IcNU), ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO),     &
                   NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF), NFAR, &
                   RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),       &
                   RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF),             &
                   RCKWRK(NcKR), RCKWRK(NcI1), RCKWRK(NcI2),             &
                   RCKWRK(NcI3))

      DO 100 I = 1, NII
         Q(I) = RCKWRK(NcI1 + I - 1) - RCKWRK(NcI2 + I - 1)
  100 CONTINUE
      RETURN
      END SUBROUTINE CKQXP

!     *****************************************************************

      SUBROUTINE CKQXR  (RHO, T, X, ICKWRK, RCKWRK, Q)
!     Returns the rates of progress for the reactions given mass
!     density, temperature and mole fractions;  see Eqs. (51) and (58).
!     INPUT
!     RHO    - Mass density. [g/cm3]
!     T      - Temperature. [K]
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!
!     OUTPUT
!     Q      - Rates of progress for the reactions. [mol/cm3/s]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION X(*), ICKWRK(*), RCKWRK(*), Q(*)
      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),              &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),          &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,             &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT), &
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),       &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),             &
                   RCKWRK(NcI1))

      CALL CKXTCR (RHO, T, X, ICKWRK, RCKWRK, RCKWRK(NcK1))

      CALL CKRATX (NII, NKK, MXSP, MXTB, T, RCKWRK(NcK1), ICKWRK(IcNS),  &
                   ICKWRK(IcNU), ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO),     &
                   NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF), NFAR, &
                   RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),       &
                   RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF),             &
                   RCKWRK(NcKR), RCKWRK(NcI1), RCKWRK(NcI2),             &
                   RCKWRK(NcI3))

      DO 100 I = 1, NII
         Q(I) = RCKWRK(NcI1 + I - 1) - RCKWRK(NcI2 + I - 1)
  100 CONTINUE
      RETURN
      END SUBROUTINE CKQXR

!     *****************************************************************

      SUBROUTINE CKQYP  (P, T, Y, ICKWRK, RCKWRK, Q)
!     Returns the rates of progress for the reactions given pressure,
!     temperature and mass fractions;  see Eqs. (51) and (58).
!
!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature. [K]
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!
!     OUTPUT
!     Q      - Rates of progress for the reactions. [mol/cm3/s]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION Y(*), ICKWRK(*), RCKWRK(*), Q(*)

      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),             &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),         &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,            &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT),&
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),      &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),            &
                   RCKWRK(NcI1))

      CALL CKYTCP (P, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK1))

      CALL CKRATX (NII, NKK, MXSP, MXTB, T, RCKWRK(NcK1), ICKWRK(IcNS), &
                   ICKWRK(IcNU), ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO),    &
                   NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF), NFAR,&
                   RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),      &
                   RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF),            &
                   RCKWRK(NcKR), RCKWRK(NcI1), RCKWRK(NcI2),            &
                   RCKWRK(NcI3))

      DO 100 I = 1, NII
         Q(I) = RCKWRK(NcI1 + I - 1) - RCKWRK(NcI2 + I - 1)
  100 CONTINUE
      RETURN
      END SUBROUTINE CKQYP
!     *****************************************************************

      SUBROUTINE CKQYR  (RHO, T, Y, ICKWRK, RCKWRK, Q)
!     Returns the rates of progress for the reactions given mass
!     density, temperature and mass fractions;  see Eqs. (51) and (58).

!     INPUT
!     RHO    - Mass density. [g/cm3]
!     T      - Temperature. [K]
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     Q      - Rates of progress for the reactions. [mol/cm3/s]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION Y(*), ICKWRK(*), RCKWRK(*), Q(*)

      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),             &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),         &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,            &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT),&
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),      &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),            &
                   RCKWRK(NcI1))

      CALL CKYTCR (RHO, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK1))

      CALL CKRATX (NII, NKK, MXSP, MXTB, T, RCKWRK(NcK1), ICKWRK(IcNS), &
                   ICKWRK(IcNU), ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO),    &
                   NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF), NFAR,&
                   RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),      &
                   RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF),            &
                   RCKWRK(NcKR), RCKWRK(NcI1), RCKWRK(NcI2),            &
                   RCKWRK(NcI3))

      DO 100 I = 1, NII
         Q(I) = RCKWRK(NcI1 + I - 1) - RCKWRK(NcI2 + I - 1)
  100 CONTINUE
      RETURN
      END SUBROUTINE CKQYR

!     *****************************************************************

      SUBROUTINE CKR2CH (RNUM, STR, I, KERR)
!     Returns a character string representation of a real number
!     and the effective length of the string.

!     INPUT
!     RNUM   - A number to be converted to a string.
!              the maximum magnitude of RNUM is machine-dependent.
!                   Data type - real scalar

!     OUTPUT
!     STR   - A left-justified character string representing RNUM,
!             with 5 to 10 characters, depending on the input value.
!             i.e., RNUM=  0.0      returns STR=" 0.00"
!                   RNUM= -10.5     returns STR="-1.05E+01"
!                   RNUM= 1.86E-100 returns in STR=" 1.86E-100"
!                   Data type - CHARACTER*(*);
!                   the minimum length of STR required is 5
!     I     - The effective length of STR
!                   Data type - integer scalar
!     KERR  - Error flag;  character length error will result in
!             KERR=.TRUE.
!                   Data type - logical
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
      CHARACTER STR*(*), INUM*3, IEXP*4
      LOGICAL KERR, IERR

      SMALL = 10.0D0**(-300)
      BIG   = 10.0D0**(+300)

      ILEN = LEN(STR)
      STR = ' '
      KERR = .FALSE.
      I = 0
      IF (ILEN .LT. 5) THEN
         KERR = .TRUE.
         RETURN
      ENDIF

      IF (RNUM .EQ. 0.0) THEN
         STR = ' 0.00'
         I = 5
         RETURN
      ENDIF

!     convert RNUM to a value between 100.0 and 999.0

      IF (RNUM.LT.-BIG.OR.RNUM.GT.BIG .OR.        &
         (RNUM.GT.0.0 .AND. RNUM.LT.SMALL) .OR.   &
         (RNUM.LT.0.0 .AND. RNUM.GT.SMALL)) THEN
         KERR = .TRUE.
         RETURN
      ENDIF

      IF (RNUM .LT. 0) THEN
          VAL = -RNUM
      ELSE
          VAL = RNUM
      ENDIF
      IE  = LOG10(VAL)

   25 CONTINUE
      IF (IE .LT. 0) THEN
         RVAL = VAL * 10.0**(ABS(IE) - 1) * 1000.0
      ELSEIF (IE .GT. 0) THEN
         RVAL = VAL * 10.0**(-IE + 1) * 10.0
      ELSE
         RVAL = VAL * 100.0
      ENDIF
      IF (RVAL.LT.100.0 .OR. RVAL.GE.1000.0) THEN
         IF (RVAL .LT. 100.0) IE = IE - 1
         IF (RVAL .GE. 1000.0)IE = IE + 1
         GO TO 25
      ELSE
         IVAL = NINT (RVAL)
         IF (IVAL .EQ. 1000) THEN
            IVAL = 100
            IF (IE .LE. 0) THEN
               IE = IE - 1
            ELSE
               IE = IE + 1
            ENDIF
         ENDIF
      ENDIF

      CALL CKI2CH (IVAL, INUM, L, IERR)
      LT = 0
      IF (IE.NE.0) THEN
         CALL CKI2CH (ABS(IE), IEXP, LEXP, IERR)
         LT = MAX(LEXP, 2) + 2
      ENDIF
      IERR = IERR.OR.(5+LT .GT. ILEN)
      IF (IERR) THEN
         KERR = .TRUE.
         RETURN
      ENDIF

      IF (RNUM .LT. 0.0) STR(1:) = '-'
      STR(2:) = INUM(:1)//'.'//INUM(2:3)
      IF (IE .NE. 0) THEN
         IF (IE .LT. 0) THEN
            STR(6:) = 'E-'
         ELSEIF (IE .GT. 0) THEN
            STR(6:) = 'E+'
         ENDIF
         IF (LEXP .EQ. 1) THEN
            STR(8:) = '0'//IEXP(:1)
         ELSE
            STR(8:) = IEXP(:LEXP)
         ENDIF
      ENDIF

      I = ILASCH(STR)
      RETURN
      END SUBROUTINE CKR2CH


!     *****************************************************************

      SUBROUTINE CKRAEX (I, RCKWRK, RA)
!     Get/put the Pre-exponential coefficient of the Ith reaction
!     INPUT
!     I      - Reaction number; I > 0 gets RA(I) from RCKWRK
!                               I < 0 puts RA(I) into RCKWRK
!                   Data type - integer scalar
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!     If I < 1:
!     RA     - Pre-exponential coefficient for the Ith reaction.
!                   cgs units - mole-cm-sec-K
!                   Data type - real scalar

!     OUTPUT
!     If I > 1:
!     RA     - Pre-exponential coefficient for Ith reaction.
!                   cgs units - mole-cm-sec-K
!                   Data type - real scalar
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION RCKWRK(*)

      NI = NcCO + (ABS(I)-1)*(NPAR+1)
      IF (I .GT. 0) THEN
         RA = RCKWRK(NI)
      ELSE
         RCKWRK(NI) = RA
      ENDIF
      RETURN
      END SUBROUTINE CKRAEX

!     *****************************************************************

      SUBROUTINE CKRAT  (RCKWRK, ICKWRK, II, KK, MAXSP, MAXTB, RU, PATM, &
                         T, C, NSPEC, NU, NUNK, NPAR, PAR, NREV, IREV,   &
                         RPAR, NFAL, IFAL, IFOP, KFAL, NFAR, FPAR, NLAN, &
                         NLAR, ILAN, PLT, NRLT, IRLT, RPLT, NTHB, ITHB,  &
                         NTBS, AIK, NKTB, SMH, RKF, RKR, EQK, CTB)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

      DIMENSION RCKWRK(*), ICKWRK(*), C(*), NSPEC(*), NU(MAXSP,*),       &
                NUNK(MAXSP,*), PAR(NPAR,*), IREV(*), RPAR(NPAR,*),       &
                ILAN(*), IRLT(*), PLT(NLAR,*), RPLT(NLAR,*),             &
                IFAL(*), IFOP(*), KFAL(*), FPAR(NFAR,*), ITHB(*),        &
                NTBS(*), AIK(MAXTB,*), NKTB(MAXTB,*), SMH(*),            &
                RKF(*), RKR(*), EQK(*), CTB(*)

!     COMMON /MACH/ SMALL,BIG,EXPARG

      ALOGT = LOG(T)

      DO 20 I = 1, II
         CTB(I) = 1.0
         RKF(I) = PAR(1,I) * EXP(PAR(2,I)*ALOGT - PAR(3,I)/T)
   20 CONTINUE

!        Landau-Teller reactions

      DO 25 N = 1, NLAN
         I = ILAN(N)
         TFAC = PLT(1,N)/T**(1.0/3.0) + PLT(2,N)/T**(2.0/3.0)
         RKF(I) = RKF(I) * EXP(TFAC)
   25 CONTINUE

      CALL CKSMH (T, ICKWRK, RCKWRK, SMH)
      DO 50 I = 1, II
          SUMSMH = 0.0
          DO 40 N = 1, MAXSP
             IF (NUNK(N,I).NE.0) SUMSMH=SUMSMH+NU(N,I)*SMH(NUNK(N,I))
   40     CONTINUE

          EQK(I) = EXP(MIN(SUMSMH,EXPARG))
   50 CONTINUE

      PFAC = PATM / (RU*T)
      DO 60 I = 1, II
         NUSUMK = NU(1,I)+NU(2,I)+NU(3,I)+NU(4,I)+NU(5,I)+NU(6,I)
         EQK(I) = EQK(I) * PFAC**NUSUMK

!     RKR=0.0 for irreversible reactions, else RKR=RKF/MAX(EQK,SMALL)

         RKR(I) = 0.0
         IF (NSPEC(I).GT.0) RKR(I) = RKF(I) / MAX(EQK(I),SMALL)
   60 CONTINUE

!     if reverse parameters have been given:

      DO 70 N = 1, NREV
         I = IREV(N)
         RKR(I) = RPAR(1,N) * EXP(RPAR(2,N)*ALOGT - RPAR(3,N)/T)
         EQK(I) = RKF(I)/RKR(I)
   70 CONTINUE

!     if reverse Landau-Teller parameters have been given:

      DO 75 N = 1, NRLT
         I = IRLT(N)
         TFAC = RPLT(1,N)/T**(1.0/3.0) + RPLT(2,N)/T**(2.0/3.0)
         RKR(I) = RKR(I) * EXP(TFAC)
         EQK(I) = RKF(I)/RKR(I)
   75 CONTINUE

!     third-body reactions

      CTOT = 0.0
      DO 10 K = 1, KK
         CTOT = CTOT + C(K)
   10 CONTINUE

      DO 80 N = 1, NTHB
         CTB(ITHB(N)) = CTOT
         DO 80 L = 1, NTBS(N)
            CTB(ITHB(N)) = CTB(ITHB(N)) + (AIK(L,N)-1.0)*C(NKTB(L,N))
   80 CONTINUE

!     If fall-off (pressure dependence):

      DO 90 N = 1, NFAL

!        CONCENTRATION OF THIRD BODY

         IF (KFAL(N) .EQ. 0) THEN
            CTHB = CTB(IFAL(N))
            CTB(IFAL(N)) = 1.0
         ELSE
            CTHB = C(KFAL(N))
         ENDIF

         RKLOW = FPAR(1,N) * EXP(FPAR(2,N)*ALOGT - FPAR(3,N)/T)
         PR = RKLOW*CTHB / RKF(IFAL(N))
         PRLOG = LOG10(MAX(PR,SMALL))

         IF (IFOP(N) .EQ. 1) THEN

!           LINDEMANN FORM

            FC = 1.0

         ELSE

            IF (IFOP(N) .EQ. 2) THEN

!              SRI FORM

               XP = 1.0/(1.0 + PRLOG**2)
               FC = ((FPAR(4,N)*EXP(-FPAR(5,N)/T) + EXP(-T/FPAR(6,N))) &
                    **XP) * FPAR(7,N) * T**FPAR(8,N)

            ELSE

!              6-PARAMETER TROE FORM

               FCENT = (1.0-FPAR(4,N)) * EXP(-T/FPAR(5,N)) &
                     + FPAR(4,N) * EXP(-T/FPAR(6,N))

!              7-PARAMETER TROE FORM

               IF (IFOP(N) .EQ. 4) FCENT = FCENT + EXP(-FPAR(7,N)/T)

               FCLOG = LOG10(MAX(FCENT,SMALL))
               XN    = 0.75 - 1.27*FCLOG
               CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
               FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
               FC = 10.0**FLOG
            ENDIF
         ENDIF
         PCOR = FC * PR/(1.0+PR)
         RKF(IFAL(N)) = RKF(IFAL(N)) * PCOR
         RKR(IFAL(N)) = RKR(IFAL(N)) * PCOR
   90 CONTINUE

!     Multiply by the product of reactants and product of products
!     PAR(4,I) is a perturbation factor

      DO 150 I = 1, II
         RKF(I) = RKF(I)*CTB(I)*C(NUNK(1,I))**IABS(NU(1,I))*PAR(4,I)
         RKR(I) = RKR(I)*CTB(I)*C(NUNK(4,I))**NU(4,I)      *PAR(4,I)
         IF (NUNK(2,I) .NE. 0) THEN
            RKF(I)= RKF(I) * C(NUNK(2,I))**IABS(NU(2,I))
            IF (NUNK(3,I) .NE. 0) &
               RKF(I) = RKF(I) * C(NUNK(3,I))**IABS(NU(3,I))
         ENDIF
         IF (NUNK(5,I) .NE. 0) THEN
            RKR(I) = RKR(I) * C(NUNK(5,I))**NU(5,I)
            IF (NUNK(6,I) .NE. 0) RKR(I) = RKR(I)*C(NUNK(6,I))**NU(6,I)
         ENDIF
  150 CONTINUE

      RETURN
      END SUBROUTINE CKRAT


!     *****************************************************************

      SUBROUTINE CKRATT (RCKWRK, ICKWRK, II, MAXSP, RU, PATM, T, NSPEC,&
                         NU, NUNK, NPAR, PAR, NREV, IREV, RPAR, NLAN,  &
                         NLAR, ILAN, PLT, NRLT, IRLT, RPLT, SMH, RKFT, &
                         RKRT, EQK)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION RCKWRK(*), ICKWRK(*), NSPEC(*), NU(MAXSP,*),          &
                NUNK(MAXSP,*), PAR(NPAR,*), IREV(*), RPAR(NPAR,*),    &
                ILAN(*), IRLT(*), PLT(NLAR,*), RPLT(NLAR,*), SMH(*),  &
                RKFT(*), RKRT(*), EQK(*)


      ALOGT = LOG(T)

      DO 20 I = 1, II
         RKFT(I) = PAR(1,I) * EXP(PAR(2,I)*ALOGT - PAR(3,I)/T)
   20 CONTINUE

!        Landau-Teller reactions

      DO 25 N = 1, NLAN
         I = ILAN(N)
         TFAC = PLT(1,N)/T**(1.0/3.0) + PLT(2,N)/T**(2.0/3.0)
         RKFT(I) = RKFT(I) * EXP(TFAC)
   25 CONTINUE

       CALL CKSMH (T, ICKWRK, RCKWRK, SMH)
       DO 50 I = 1, II
          SUMSMH = 0.0
          DO 40 N = 1, MAXSP
             IF (NUNK(N,I).NE.0) SUMSMH=SUMSMH+NU(N,I)*SMH(NUNK(N,I))
   40     CONTINUE
          EQK(I) = EXP(MIN(SUMSMH,EXPARG))
   50 CONTINUE

      PFAC = PATM / (RU*T)
      DO 60 I = 1, II
         NUSUMK = NU(1,I)+NU(2,I)+NU(3,I)+NU(4,I)+NU(5,I)+NU(6,I)
         EQK(I) = EQK(I) * PFAC**NUSUMK

!     RKRT=0.0 for irreversible reactions, else RKRT=RKFT/MAX(EQK,SMALL)

         RKRT(I) = 0.0
         IF (NSPEC(I).GT.0) RKRT(I) = RKFT(I) / MAX(EQK(I),SMALL)
   60 CONTINUE

!     if reverse parameters have been given:

      DO 70 N = 1, NREV
         I = IREV(N)
         RKRT(I) = RPAR(1,N) * EXP(RPAR(2,N)*ALOGT - RPAR(3,N)/T)
         EQK(I)  = RKFT(I)/RKRT(I)
   70 CONTINUE

!     if reverse Landau-Teller parameters have been given:

      DO 75 N = 1, NRLT
         I = IRLT(N)
         TFAC = RPLT(1,N)/T**(1.0/3.0) + RPLT(2,N)/T**(2.0/3.0)
         RKRT(I) = RKRT(I) * EXP(TFAC)
         EQK(I) = RKFT(I)/RKRT(I)
   75 CONTINUE



      RETURN
      END SUBROUTINE CKRATT


!     *****************************************************************


      SUBROUTINE CKRATX (II, KK, MAXSP, MAXTB, T, C, NSPEC, NU, NUNK,   &
                         NPAR, PAR, NFAL, IFAL, IFOP, KFAL, NFAR, FPAR, &
                         NTHB, ITHB, NTBS, AIK, NKTB, RKFT, RKRT, RKF,  &
                         RKR, CTB)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION C(*), NSPEC(*), NU(MAXSP,*), NUNK(MAXSP,*), PAR(NPAR,*),&
                IFAL(*), IFOP(*), KFAL(*), FPAR(NFAR,*), ITHB(*),       &
                NTBS(*), AIK(MAXTB,*), NKTB(MAXTB,*), RKFT(*),          &
                RKRT(*), RKF(*), RKR(*), CTB(*)

      DO 20 I = 1, II
         CTB(I) = 1.d0
         RKF(I) = RKFT(I)
         RKR(I) = RKRT(I)
   20 CONTINUE

!     third-body reactions

      CTOT = 0.d0
      DO 10 K = 1, KK
         CTOT = CTOT + C(K)
   10 CONTINUE

      DO 80 N = 1, NTHB
         CTB(ITHB(N)) = CTOT
         DO 80 L = 1, NTBS(N)
            CTB(ITHB(N)) = CTB(ITHB(N)) + (AIK(L,N)-1.d0)*C(NKTB(L,N))
   80 CONTINUE

!     If fall-off (pressure correction):

      ALOGT = LOG(T)

      DO 90 N = 1, NFAL

         RKLOW = FPAR(1,N) * EXP(FPAR(2,N)*ALOGT - FPAR(3,N)/T)

!        CONCENTRATION OF THIRD BODY

         IF (KFAL(N) .EQ. 0) THEN
            PR = RKLOW * CTB(IFAL(N)) / RKF(IFAL(N))
            CTB(IFAL(N)) = 1.d0
         ELSE
            PR = RKLOW * C(KFAL(N)) / RKF(IFAL(N))
         ENDIF

         PCOR = PR / (1.d0 + PR)

         IF (IFOP(N) .GT. 1) THEN
            PRLOG = LOG10(MAX(PR,SMALL))

            IF (IFOP(N) .EQ. 2) THEN

!              8-PARAMETER SRI FORM

               XP = 1.d0/(1.d0 + PRLOG**2)
               FC = ((FPAR(4,N)*EXP(-FPAR(5,N)/T) + EXP(-T/FPAR(6,N))) &
                    **XP) * FPAR(7,N) * T**FPAR(8,N)

            ELSE

!              6-PARAMETER TROE FORM

               FCENT = (1.d0-FPAR(4,N)) * EXP(-T/FPAR(5,N)) &
                     +       FPAR(4,N) * EXP(-T/FPAR(6,N))

!              7-PARAMETER TROE FORM

               IF (IFOP(N) .EQ. 4) FCENT = FCENT + EXP(-FPAR(7,N)/T)

               FCLOG = LOG10(MAX(FCENT,SMALL))
               XN    = 0.75d0 - 1.27d0*FCLOG
               CPRLOG= PRLOG - (0.4d0 + 0.67d0*FCLOG)
               FLOG = FCLOG/(1.d0 + (CPRLOG/(XN-0.14d0*CPRLOG))**2)
               FC = 10.d0**FLOG
            ENDIF
            PCOR = FC * PCOR
         ENDIF

         RKF(IFAL(N)) = RKF(IFAL(N)) * PCOR
         RKR(IFAL(N)) = RKR(IFAL(N)) * PCOR


   90 CONTINUE


!     Multiply by the product of reactants and product of products

      DO 150 I = 1, II
         RKF(I) = RKF(I)*CTB(I)*C(NUNK(1,I))**IABS(NU(1,I))
         RKR(I) = RKR(I)*CTB(I)*C(NUNK(4,I))**NU(4,I)
         IF (NUNK(2,I) .NE. 0) THEN
            RKF(I) = RKF(I) * C(NUNK(2,I))**IABS(NU(2,I))
            IF (NUNK(3,I) .NE. 0) &
               RKF(I) = RKF(I) * C(NUNK(3,I))**IABS(NU(3,I))
         ENDIF
         IF (NUNK(5,I) .NE. 0) THEN
            RKR(I) = RKR(I) * C(NUNK(5,I))**NU(5,I)
            IF (NUNK(6,I) .NE. 0) &
               RKR(I) = RKR(I) * C(NUNK(6,I))**NU(6,I)
         ENDIF
  150 CONTINUE



!     Perturbation factor

      DO 160 I = 1, II
         RKF(I) = RKF(I) * PAR(4,I)
         RKR(I) = RKR(I) * PAR(4,I)
  160 CONTINUE




      RETURN
      END SUBROUTINE CKRATX


!     *****************************************************************

      SUBROUTINE CKRDEX (I, RCKWRK, RD)
!     Get/put the perturbation factor of the Ith reaction

!     INPUT
!     I      - Reaction number; I > 0 gets RD(I) from RCKWRK
!                               I < 0 puts RD(I) into RCKWRK
!                   Data type - integer scalar
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.
!     If I < 1:
!     RD     - Perturbation factor for the Ith reaction.
!                   cgs units - mole-cm-sec-K
!                   Data type - real scalar

!     OUTPUT
!     If I > 1:
!     RD     - Perturbation factor for Ith reaction.
!                   cgs units - mole-cm-sec-K
!                   Data type - real scalar
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION RCKWRK(*)

      NI = NcCO + (ABS(I)-1)*(NPAR+1) + NPAR
      IF (I .GT. 0) THEN
         RD = RCKWRK(NI)
      ELSE
         RCKWRK(NI) = RD
      ENDIF
      RETURN
      END SUBROUTINE CKRDEX

!     *****************************************************************

      SUBROUTINE CKRHEX (K, RCKWRK, A6)
!     Returns an array of the sixth thermodynamic polynomial
!     coefficients for a species, or changes their value,
!     depending on the sign of K.

!     INPUT
!      K      - Integer species number; K>0 gets A6(*) from RCKWRK,
!                                       K<0 puts A6(*) into RCKWRK.
!                    Data type - integer scalar
!      RCKWRK - Array of real internal work space.

!     OUTPUT
!      A6     - The array of the 6th thermodynamic polynomial
!               coefficients for the Kth species, over the number
!               of temperature ranges used in fitting thermodynamic
!               properties.
!               Dimension A6(*) at least (MXTP-1), where MXTP is
!               the maximum number of temperatures used for fitting
!               the thermodynamic properties of the species.
!                    Data type - real array
!                    cgs units:  none
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION RCKWRK(*), A6(*)
      DO 100 L = 1, MXTP-1
         NA6 = NCAA + (L-1)*NCP2 + (ABS(K)-1)*NCP2T + NCP
         IF (K .GT. 0) THEN
            A6(L) = RCKWRK(NA6)
         ELSE
            RCKWRK(NA6) = A6(L)
         ENDIF
  100 CONTINUE

      RETURN
      END SUBROUTINE CKRHEX


!     *****************************************************************

      SUBROUTINE CKRHOC (P, T, C, ICKWRK, RCKWRK, RHO)
!     Returns the mass density of the gas mixture given the pressure,
!     temperature and molar concentrations;  see Eq. (2).

!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature. [K]
!     C      - Molar concentrations of the species. [mol/cm3]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     RHO    - Mass density. [g/cm3]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION C(*), ICKWRK(*), RCKWRK(*)
      CTOT = 0.0
      SUM  = 0.0
      DO 100 K = 1, NKK
         CTOT = CTOT + C(K)
         SUM = SUM + C(K)*RCKWRK(NcWT + K - 1)
  100 CONTINUE

      RHO  = SUM * P / (RCKWRK(NcRU)*T*CTOT)
      RETURN
      END SUBROUTINE CKRHOC


!     *****************************************************************

      SUBROUTINE CKRHOX (P, T, X, ICKWRK, RCKWRK, RHO)
!     Returns the mass density of the gas mixture given the pressure,
!     temperature and mole fractions;  see Eq. (2).
!
!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature. [K]
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!
!     OUTPUT
!     RHO    - Mass density. [g/cm3]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION X(*), ICKWRK(*), RCKWRK(*)

      SUM = 0.0
      DO 100 K = 1, NKK
         SUM = SUM + X(K)*RCKWRK(NcWT + K - 1)
  100 CONTINUE

      RHO = SUM * P / (RCKWRK(NcRU)*T)
      RETURN
      END SUBROUTINE CKRHOX

!     *****************************************************************

      SUBROUTINE CKRHOY (P, T, Y, ICKWRK, RCKWRK, RHO)
!     Returns the mass density of the gas mixture given the pressure,
!     temperature and mass fractions;  see Eq. (2).

!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature. [K]
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     RHO    - Mass density. [g/cm3]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), Y(*)

      SUMYOW = 0.d0
      DO 150 K = 1, NKK
         SUMYOW = SUMYOW + Y(K)/RCKWRK(NcWT + K - 1)
150   CONTINUE
      RHO = P/(SUMYOW*T*RCKWRK(NcRU))
      RETURN
      END SUBROUTINE CKRHOY

!     *****************************************************************

      SUBROUTINE CKRP   (ICKWRK, RCKWRK, RU, RUC, PA)
!     Returns universal gas constants and the pressure of one standard
!     atmosphere

!     INPUT
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     RU     - Universal gas constant.
!                   cgs units - 8.314E7 ergs/(mole*K)
!                   Data type - real scalar
!     RUC    - Universal gas constant used only in conjuction with
!              activation energy.
!                   preferred units - 1.987 cal/(mole*K)
!                   Data type - real scalar
!     PA     - Pressure of one standard atmosphere.
!                   cgs units - 1.01325E6 dynes/cm**2
!                   Data type - real scalar
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*)
      RU  = RCKWRK(NcRU)
      RUC = RCKWRK(NcRC)
      PA  = RCKWRK(NcPA)
      RETURN
      END SUBROUTINE CKRP

!     *****************************************************************

      SUBROUTINE CKSAVE (LOUT, LSAVE, ICKWRK, RCKWRK, CCKWRK)
!     Writes to a binary file information about a Chemkin
!     linking file, pointers for the Chemkin Library, and
!     Chemkin work arrays.

!     INPUT
!     LOUT   - Output file for printed diagnostics.
!                   Data type - integer scalar
!     LSAVE  - Integer output unit.
!                   Data type - integer scalar
!     ICKWRK - Array of integer workspace containing integer data.
!                   Data type - integer array
!     RCKWRK - Array of real workspace containing real data.
!                   Data type - real array
!     CCKWRK - Array of character workspace containing character data.
!                   Data type - CHARACTER*18 array
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*)
      CHARACTER CCKWRK(*)*(*)
!     LOGICAL KERR
!     COMMON /CKCONS/ PREC, VERS, KERR, LENI, LENR, LENC

      NPOINT = 63
      WRITE (LSAVE, ERR=999)                                         &
                      NPOINT, VERS,   PREC,   LENI,   LENR,   LENC,  &
                      NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,&
                      NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,&
                      NTHB, NRLT, NWL,  IcMM, IcKK, IcNC, IcPH, IcCH,&
                      IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL, IcRV,&
                      IcWL, IcFL, IcFO, IcKF, IcTB, IcKN, IcKT, NcAW,&
                      NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL,&
                      NcKT, NcWL, NcRU, NcRC, NcPA, NcKF, NcKR, NcK1,&
                      NcK2, NcK3, NcK4, NcI1, NcI2, NcI3, NcI4
      WRITE (LSAVE, ERR=999) (ICKWRK(L), L = 1, LENI)
      WRITE (LSAVE, ERR=999) (RCKWRK(L), L = 1, LENR)
      WRITE (LSAVE, ERR=999) (CCKWRK(L), L = 1, LENC)
      RETURN

  999 CONTINUE
      WRITE (LOUT, *) &
       ' Error writing Chemkin linking file information...'
      KERR = .TRUE.
      RETURN
      END SUBROUTINE CKSAVE

!     *****************************************************************

      SUBROUTINE CKSBML (P, T, X, ICKWRK, RCKWRK, SBML)
!     Returns the mean entropy of the mixture in molar units,
!     given the pressure, temperature and mole fractions;
!     see Eq. (42).

!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature. [K]
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     SBML   - Mean entropy in molar units. [erg/mol/K]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION X(*), ICKWRK(*), RCKWRK(*)
!     COMMON /MACH/ SMALL,BIG,EXPARG

      CALL CKSML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))

      RLNP = RCKWRK(NcRU) * LOG(P / RCKWRK(NcPA))
      SBML = 0.0
      DO 100 K = 1, NKK
         SBML = SBML + X(K) * ( RCKWRK(NcK1 + K - 1) -   &
                RCKWRK(NcRU)*LOG(MAX(X(K),SMALL)) - RLNP )
  100 CONTINUE

      RETURN
      END SUBROUTINE CKSBML


!     *****************************************************************

      SUBROUTINE CKSBMS (P, T, Y, ICKWRK, RCKWRK, SBMS)

!     Returns the mean entropy of the mixture in mass units,
!     given the pressure, temperature and mass fractions;
!     see Eq.(43).

!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature. [K]
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     SBMS   - Mean entropy in mass units. [erg/g/K]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION Y(*), ICKWRK(*), RCKWRK(*)
!     COMMON /MACH/ SMALL,BIG,EXPARG

      CALL CKSML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKYTX (Y, ICKWRK, RCKWRK, RCKWRK(NcK2))

      RLNP = RCKWRK(NcRU) * LOG (P / RCKWRK(NcPA))
      SUM = 0.0
      DO 100 K = 1, NKK
         SUM = SUM + RCKWRK(NcK2 + K - 1) *  &
                   ( RCKWRK(NcK1 + K - 1)    &
                   - RCKWRK(NcRU) *          &
                     LOG(MAX(RCKWRK(NcK2 + K - 1),SMALL)) - RLNP)
  100 CONTINUE

      CALL CKMMWY (Y, ICKWRK, RCKWRK, WTM)
      SBMS = SUM / WTM
      RETURN
      END SUBROUTINE CKSBMS


!     *****************************************************************

      SUBROUTINE CKSMH  (T, ICKWRK, RCKWRK, SMH)
!     Returns the array of entropies minus enthalpies for the species.
!     It is normally not called directly by the user.

!     INPUT
!     T      - Temperature. [K]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     SMH    - Entropy minus enthalpy for the species,
!              SMH(K) = S(K)/R - H(K)/RT.
!                   cgs units - none
!                   Data type - real array
!                   Dimension SMH(*) at least KK, the total number of
!                   species.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), SMH(*), TN(10)
      TN(1) = LOG(T) - 1.0
      DO 150 N = 2, NCP
         TN(N) = T**(N-1)/((N-1)*N)
 150  CONTINUE

      DO 250 K = 1, NKK
         L = 1
         DO 220 N = 2, ICKWRK(IcNT + K - 1)-1
            TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE

         NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
         SUM = 0.0
         DO 225 N = 1, NCP
            SUM = SUM + TN(N)*RCKWRK(NA1 + N - 1)
  225    CONTINUE
         SMH(K) = SUM + RCKWRK(NA1 + NCP2 - 1) &
                      - RCKWRK(NA1 + NCP1 - 1)/T

 250  CONTINUE
      RETURN
      END SUBROUTINE CKSMH

!     *****************************************************************

      SUBROUTINE CKSML  (T, ICKWRK, RCKWRK, SML)
!     Returns the standard state entropies in molar units

!     INPUT
!     T      - Temperature. [K]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     SML    - Standard state entropies in molar units for the species. [erg/mol/K]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), SML(*), TN(10)

      TN(1) = LOG(T)
      DO 150 N = 2, NCP
         TN(N) = T**(N-1)/(N-1)
150   CONTINUE

      DO 250 K = 1, NKK
         L = 1
         DO 220 N = 2, ICKWRK(IcNT + K - 1)-1
            TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE

         NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
         SUM = 0.0
         DO 225 N = 1, NCP
            SUM = SUM + TN(N)*RCKWRK(NA1 + N - 1)
  225    CONTINUE
         SML(K) = RCKWRK(NcRU) * (SUM + RCKWRK(NA1+NCP2-1))
250   CONTINUE
      RETURN
      END SUBROUTINE CKSML

!     *****************************************************************

      SUBROUTINE CKSMS  (T, ICKWRK, RCKWRK, SMS)
!     Returns the standard state entropies in mass units;
!     see Eq. (28).

!     INPUT
!     T      - Temperature. [K]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     SMS    - Standard state entropies in mass units for the species. [erg/g/K]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), SMS(*), TN(10)
      TN(1) = LOG(T)
      DO 150 N = 2, NCP
         TN(N) = T**(N-1)/(N-1)
150   CONTINUE

      DO 250 K = 1, NKK
         L = 1
         DO 220 N = 2, ICKWRK(IcNT + K - 1)-1
            TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE

         NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
         SUM = 0.0
         DO 225 N = 1, NCP
            SUM = SUM + TN(N)*RCKWRK(NA1 + N - 1)
  225    CONTINUE
         SMS(K) = RCKWRK(NcRU) * (SUM+RCKWRK(NA1 + NCP2 - 1)) &
                               / RCKWRK(NcWT + K - 1)
250   CONTINUE
      RETURN
      END SUBROUTINE CKSMS

!     *****************************************************************

      SUBROUTINE CKSNUM (LINE, NEXP, LOUT, KRAY, NN, KNUM, NVAL, &
                         RVAL, KERR)
!     This subroutine is called to parse a character string, LINE,
!     that is composed of several blank-delimited substrings.
!     It is expected that the first substring in LINE is also an
!     entry in a reference array of character strings, KRAY(*), in
!     which case the index position in KRAY(*) is returned as KNUM,
!     otherwise an error flag is returned.  The substrings following
!     the first are expected to represent numbers, and are converted
!     to elements of the array RVAL(*).  If NEXP substrings are not
!     found an error flag will be returned.  This allows format-free
!     input of combined alpha-numeric data.  For example, after
!     reading a line containing a species name followed by several
!     numerical values, the subroutine might be called to find
!     a Chemkin species index and convert the other substrings to
!     real values:

!     input:  LINE    = "N2  1.2"
!             NEXP    = 1, the number of values expected
!             LOUT    = 6, a logical unit number on which to write
!                       diagnostic messages.
!             KRAY(*) = "H2" "O2" "N2" "H" "O" "N" "OH" "H2O" "NO"
!             NN      = 9, the number of entries in KRAY(*)
!     output: KNUM    = 3, the index number of the substring in
!                       KRAY(*) which corresponds to the first
!                       substring in LINE
!             NVAL    = 1, the number of values found in LINE
!                       following the first substring
!             RVAL(*) = 1.200E+00, the substring converted to a number
!             KERR    = .FALSE.
!     INPUT
!     LINE   - A character string.
!                   Data type - CHARACTER*80
!     NEXP   - Number of real values to be found in character string.
!              If NEXP is negative, then ABS(NEXP) values are
!              expected.  However, it is not an error condition,
!              if less values are found.
!                   Data type - integer scalar
!     LOUT   - Output unit for printed diagnostics.
!                   Data type - integer scalar
!     KRAY   - Array of character strings.
!                   Data type - CHARACTER*(*)
!     NN     - Total number of character strings in KRAY.
!                   Data type - integer scalar

!     OUTPUT
!     KNUM   - Index number of character string in array which
!              corresponds to the first substring in LINE.
!                   Data type - integer scalar
!     NVAL   - Number of real values found in LINE.
!                   Data type - integer scalar
!     RVAL   - Array of real values found in LINE.
!                   Data type - real array
!                   Dimension RVAL(*) at least NEXP
!     KERR   - Error flag; syntax or dimensioning error,
!              corresponding string not found, or total of
!              values found is not the number of values expected,
!              will result in KERR = .TRUE.
!                   Data type - logical

!     A '!' will comment out a line, or remainder of the line.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

      CHARACTER LINE*(*), KRAY(*)*(*), ISTR*80
      DIMENSION RVAL(*)
      LOGICAL KERR, IERR

      NVAL = 0
      KERR = .FALSE.
      ILEN = MIN (IPPLEN(LINE), ILASCH(LINE))
      IF (ILEN .LE. 0) RETURN

      I1 = IFIRCH(LINE(:ILEN))
      I3 = INDEX(LINE(I1:ILEN),' ')
      IF (I3 .EQ. 0) I3 = ILEN - I1 + 1
      I2 = I1 + I3
      ISTR = ' '
      ISTR = LINE(I1:I2-1)

      CALL CKCOMP (ISTR, KRAY, NN, KNUM)
      IF (KNUM.EQ.0) THEN
         LT = MAX (ILASCH(ISTR), 1)
         WRITE (LOUT,'(A)') &
         ' Error in CKSNUM...'//ISTR(:LT)//' not found...'
         WRITE (*   ,'(A)') &
         ' Error in CKSNUM...'//ISTR(:LT)//' not found...'
         KERR = .TRUE.
      ENDIF

      ISTR = ' '
      ISTR = LINE(I2:ILEN)
      IF (NEXP .NE. 0) &
            CALL CKXNUM (ISTR, NEXP, LOUT, NVAL, RVAL, IERR)

      RETURN
      END SUBROUTINE CKSNUM

!     *****************************************************************

      SUBROUTINE CKSOR  (T, ICKWRK, RCKWRK, SOR)
!     Returns the nondimensional entropies;  see Eq. (21).

!     INPUT
!     T      - Temperature. [K]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     SOR    - Nondimensional entropies for the species.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION TN(10), SOR(*), ICKWRK(*), RCKWRK(*)
      TN(1) = LOG(T)
      DO 150 N = 2, NCP
         TN(N) = T**(N-1)/(N-1)
150   CONTINUE

      DO 250 K = 1, NKK
         L = 1
         DO 220 N = 2, ICKWRK(IcNT + K - 1)-1
            TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE

         NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
         SUM = 0.0
         DO 225 N = 1, NCP
            SUM = SUM + TN(N)*RCKWRK(NA1 + N - 1)
  225    CONTINUE
         SOR(K) = SUM + RCKWRK(NA1 + NCP2 - 1)
250   CONTINUE
      RETURN
      END SUBROUTINE CKSOR

!     *****************************************************************

      SUBROUTINE CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, KERR)
!     Returns an array of substrings in a character string with blanks
!     as the delimiter

!     INPUT
!     LINE   - A character string.
!     LOUT   - Output unit for printed diagnostics.
!     NDIM   - Dimension of array SUB(*)*(*)

!     OUTPUT
!     SUB    - The character substrings of LINE.
!                   Data type - CHARACTER*(*) array
!                   Dimension SUB(*) at least NDIM
!     NFOUND - Number of substrings found in LINE.
!                   Data type - integer
!     KERR   - Error flag; dimensioning errors will result in
!              KERR = .TRUE.
!                   Data type - logical
!     A '!' will comment out a line, or remainder of the line.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      CHARACTER SUB(*)*(*), LINE*(*)
      LOGICAL KERR
      NFOUND = 0
      ILEN = LEN(SUB(1))

      IEND = 0
      KERR = .FALSE.
   25 CONTINUE

      ISTART = IEND + 1
      DO 100 L = ISTART, IPPLEN(LINE)

         IF (LINE(L:L) .NE. ' ') THEN
            IEND   = INDEX(LINE(L:), ' ')
            IF (IEND .EQ. 0) THEN
               IEND = IPPLEN(LINE)
            ELSE
               IEND = L + IEND - 1
            ENDIF
            IF (IEND-L+1 .GT. ILEN) THEN
               WRITE (LOUT,*) ' Error in CKSUBS...substring too long'
               KERR = .TRUE.
            ELSEIF (NFOUND+1 .GT. NDIM) THEN
               WRITE (LOUT,*) ' Error in CKSUBS...NDIM too small'
               KERR = .TRUE.
            ELSE
               NFOUND = NFOUND + 1
               SUB(NFOUND) = LINE(L:IEND)
            ENDIF
            GO TO 25
         ENDIF

  100 CONTINUE
      RETURN
      END SUBROUTINE CKSUBS

!     *****************************************************************

      SUBROUTINE CKSYME (CCKWRK, LOUT, ENAME, KERR)
!     Returns the character strings of element names.
!     INPUT
!     CCKWRK - Array of character work space.
!     LOUT   - Output unit for printed diagnostics.

!     OUTPUT
!     ENAME  - Element names.
!                   Data type - CHARACTER*(*) array
!                   Dimension ENAME at least MM, the total number of
!                   elements in the problem.
!     KERR   - Error flag; character length error will result in
!              KERR = .TRUE.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      CHARACTER*(*) CCKWRK(*), ENAME(*)
      LOGICAL KERR
      KERR = .FALSE.
      ILEN = LEN(ENAME(1))
      DO 150 M = 1, NMM
         LT = ILASCH(CCKWRK(IcMM+M-1))
         ENAME(M) = ' '
         IF (LT .LE. ILEN) THEN
            ENAME(M) = CCKWRK(IcMM+M-1)
         ELSE
            WRITE (LOUT,'(A)') &
            ' Error in CKSYME...character string length too small '
            KERR = .TRUE.
         ENDIF
150   CONTINUE
      RETURN
      END SUBROUTINE CKSYME

!     *****************************************************************

      SUBROUTINE CKSYMR (I, LOUT, ICKWRK, RCKWRK, CCKWRK, LT, ISTR, &
                         KERR)
!     Returns a character string which describes the Ith reaction,
!     and the effective length of the character string.

!     INPUT
!     I      - Reaction index.
!     LOUT   - Output unit for printed diagnostics.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.
!     CCKWRK - Array of character work space.

!     OUTPUT
!     ISTR   - Character string describing the Ith reaction.
!     LT     - Number of characters in the reaction description.
!     KERR   - Error flag;  character length error will result in
!              KERR=.TRUE.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*),RCKWRK(*)
      CHARACTER CCKWRK(*)*(*), ISTR*(*), IDUM*80
      LOGICAL KERR, IERR

      ISTR = ' '
      ILEN = LEN(ISTR)
      KERR = .FALSE.

      DO 100 J = 1,2
         NS = 0
         DO 50 N = 1, MXSP
            NU = ICKWRK(IcNU + (I-1)*MXSP + N - 1)
            K  = ICKWRK(IcNK + (I-1)*MXSP + N - 1)

            IF (J.EQ.1.AND.NU.LT.0 .OR. J.EQ.2.AND.NU.GT.0) THEN
               NS = NS + 1

               IF (NS .GT. 1) THEN
                  LT = ILASCH(ISTR)
                  IF (LT+1 .GT. ILEN) THEN
                     KERR = .TRUE.
                     ISTR = ' '
                     WRITE (LOUT, 500)
                     RETURN
                  ENDIF
                  ISTR(LT+1:) = '+'
               ENDIF
               CALL CKI2CH (ABS(NU), IDUM, L, IERR)
               IF (IERR) THEN
                  KERR = .TRUE.
                  WRITE (LOUT,*) ' Syntax error in CKSYMR...'
                  ISTR = ' '
                  RETURN
               ENDIF
               IF (ABS(NU) .GT. 1) THEN
                  LT = ILASCH(ISTR)
                  IF (LT+L .GT. ILEN) THEN
                      KERR = .TRUE.
                      ISTR = ' '
                      WRITE (LOUT, 500)
                      RETURN
                  ENDIF
                  ISTR(LT+1:) = IDUM
               ENDIF
               LK = ILASCH(CCKWRK(IcKK+K-1))
               LT = ILASCH(ISTR)
               IF (LT+LK .GT. ILEN) THEN
                  KERR = .TRUE.
                  ISTR = ' '
                  WRITE (LOUT, 500)
                  RETURN
               ENDIF
               ISTR(LT+1:) = CCKWRK(IcKK+K-1)(:LK)
            ENDIF
   50    CONTINUE

         DO 60 N = 1, NFAL
            IF (ICKWRK(IcFL+N-1) .EQ. I) THEN
               LT = ILASCH(ISTR)
               IF (ICKWRK(IcKF+N-1) .EQ. 0) THEN
                  IF (LT+4 .GT. ILEN) THEN
                     KERR = .TRUE.
                     ISTR = ' '
                     WRITE (LOUT, 500)
                     RETURN
                  ENDIF
                  ISTR(LT+1:) = '(+M)'
               ELSE
                  IDUM = ' '
                  IDUM = CCKWRK (IcKK + ICKWRK(IcKF+N-1) - 1)
                  LK = ILASCH(IDUM)
                  IF (LT+LK+3 .GT. ILEN) THEN
                     KERR = .TRUE.
                     ISTR = ' '
                     WRITE (LOUT, 500)
                     RETURN
                  ENDIF
                  ISTR(LT+1:) ='(+'//IDUM(:LK)//')'
               ENDIF
            ENDIF
   60    CONTINUE

         DO 70 N = 1, NTHB
            IF (ICKWRK(IcTB+N-1).EQ.I .AND. INDEX(ISTR,'(+M)').LE.0) THEN
               LT = ILASCH(ISTR)
               IF (LT+2 .GT. ILEN) THEN
                  KERR = .TRUE.
                     ISTR = ' '
                     WRITE (LOUT, 500)
                  RETURN
               ENDIF
               ISTR(LT+1:) = '+M'
            ENDIF
   70    CONTINUE

         DO 80 N = 1, NWL
            IF (ICKWRK(IcWL+N-1) .EQ. I) THEN
               W = RCKWRK(NcWL+N-1)
               LT = ILASCH(ISTR)
               IF (LT+3 .GT. ILEN) THEN
                  KERR = .TRUE.
                  ISTR = ' '
                  WRITE (LOUT, 500)
                  RETURN
               ENDIF
               IF (J.EQ.1.AND.W.LT.0.0 .OR. J.EQ.2.AND.W.GT.0.0) &
                   ISTR(LT+1:) = '+HV'
            ENDIF
   80    CONTINUE

         IF (J.EQ.1) THEN
            LT = ILASCH(ISTR)
            IF (ICKWRK(IcNS+I-1) .LT. 0) THEN
               IF (LT+2 .GT. ILEN) THEN
                  KERR = .TRUE.
                  ISTR = ' '
                  WRITE (LOUT, 500)
                  RETURN
               ENDIF
               ISTR(LT+1:) = '=>'
            ELSE
               IF (LT+3 .GT. ILEN) THEN
                  KERR = .TRUE.
                  ISTR = ' '
                  WRITE (LOUT, 500)
                  RETURN
               ENDIF
               ISTR(LT+1:) = '<=>'
            ENDIF
         ENDIF
  100 CONTINUE
      LT = ILASCH(ISTR)

  500 FORMAT (' Error in CKSYMR...character string length too small')
      RETURN
      END SUBROUTINE CKSYMR

!     *****************************************************************

      SUBROUTINE CKSYMS (CCKWRK, LOUT, KNAME, KERR)
!     Returns the character strings of species names

!     INPUT
!     CCKWRK - Array of character work space.
!     LOUT   - Output unit for printed diagnostics.

!     OUTPUT
!     KNAME  - Species names.
!     KERR   - Error flag; character length errors will result in
!              KERR = .TRUE.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      CHARACTER*(*) CCKWRK(*), KNAME(*)
      LOGICAL KERR
      KERR = .FALSE.
      ILEN = LEN(KNAME(1))
      DO 150 K = 1, NKK
         LT = ILASCH(CCKWRK(IcKK + K - 1))
         KNAME(K) = ' '
         IF (LT .LE. ILEN) THEN
            KNAME(K) = CCKWRK(IcKK+K-1)
         ELSE
            WRITE (LOUT,*) &
            ' Error in CKSYM...character string length too small '
            KERR = .TRUE.
         ENDIF
150   CONTINUE
      RETURN
      END SUBROUTINE CKSYMS

!     *****************************************************************

      SUBROUTINE CKTHB  (KDIM, ICKWRK, RCKWRK, AKI)
!     Returns matrix of enhanced third body coefficients;
!     see Eq. (58).

!     INPUT
!     KDIM   - First dimension of the two dimensional array AKI;
!              KDIM must be greater than or equal to the total
!              number of species, KK.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     AKI    - Matrix of enhanced third body efficiencies of the
!              species in the reactions; AKI(K,I) is the enhanced
!              efficiency of the Kth species in the Ith reaction.
!                   Data type - real array
!                   Dimension AKI(KDIM,*) exactly KDIM (at least KK,
!                   the total number of species) for the first
!                   dimension and at least II for the second, the total
!                   number of reactions.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION AKI(KDIM,*), ICKWRK(*), RCKWRK(*)
      DO 150 I = 1, NII
         DO 140 K = 1, NKK
            AKI(K,I) = 1.0
  140    CONTINUE
  150 CONTINUE

      DO 250 N = 1, NTHB
         I = ICKWRK(IcTB + N - 1)
         DO 250 L = 1, ICKWRK(IcKN + N - 1)
            K  = ICKWRK(IcKT + (N-1)*MXTB + L - 1)
            AK = RCKWRK(NcKT + (N-1)*MXTB + L - 1)
            AKI(K,I) = AK
  250 CONTINUE
      RETURN
      END SUBROUTINE CKTHB


!     *****************************************************************

      SUBROUTINE CKUBML (T, X, ICKWRK, RCKWRK, UBML)
!     Returns the mean internal energy of the mixture in molar units;
!     see Eq. (39).

!     INPUT
!     T      - Temperature. [K]
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     UBML   - Mean internal energy in molar units. [erg/mol]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION X(*), ICKWRK(*), RCKWRK(*)
      CALL CKUML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))

      UBML = 0.0
      DO 100 K = 1, NKK
         UBML = UBML + X(K)*RCKWRK(NcK1 + K - 1)
  100 CONTINUE
      RETURN
      END SUBROUTINE CKUBML

!     *****************************************************************

      SUBROUTINE CKUBMS (T, Y, ICKWRK, RCKWRK, UBMS)
!     Returns the mean internal energy of the mixture in mass units;
!     see Eq. (40).

!     INPUT
!     T      - Temperature. [K]
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     UBMS   - Mean internal energy in mass units. [erg/g]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION Y(*), ICKWRK(*), RCKWRK(*)

      CALL CKUMS (T, ICKWRK, RCKWRK, RCKWRK(NcK1))

      UBMS = 0.0
      DO 100 K = 1, NKK
         UBMS = UBMS + Y(K)*RCKWRK(NcK1 + K - 1)
  100 CONTINUE
      RETURN
      END SUBROUTINE CKUBMS

!     *****************************************************************

      SUBROUTINE CKUML  (T, ICKWRK, RCKWRK, UML)
!     Returns the internal energies in molar units;  see Eq. (23).

!     INPUT
!     T      - Temperature. [K]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     UML    - Internal energies in molar units for the species. [erg/mol]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), UML(*)

      CALL CKHML (T, ICKWRK, RCKWRK, UML)
      RUT = T*RCKWRK(NcRU)
      DO 150 K = 1, NKK
         UML(K) = UML(K) - RUT
150   CONTINUE
      RETURN
      END SUBROUTINE CKUML

!     *****************************************************************

      SUBROUTINE CKUMS  (T, ICKWRK, RCKWRK, UMS)
!     Returns the internal energies in mass units;  see Eq. (30).

!     INPUT
!     T      - Temperature. [K]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     UMS    - Internal energies in mass units for the species. [erg/g]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), UMS(*)

      CALL CKHMS (T, ICKWRK, RCKWRK, UMS)
      RUT = T*RCKWRK(NcRU)
      DO 150 K = 1, NKK
         UMS(K) = UMS(K) - RUT/RCKWRK(NcWT+K-1)
150   CONTINUE
      RETURN
      END SUBROUTINE CKUMS

!     *****************************************************************

      SUBROUTINE CKWL   (ICKWRK, RCKWRK, WL)
!     Returns a set of flags providing information on the wave length
!     of photon radiation

!     INPUT
!     ICKWRK - Array of integer workspace.
!                   Data type - integer array
!                   Dimension ICKWRK(*) at least LENIWK.
!     RCKWRK - Array of real work space.
!                   Data type - real array
!                   Dimension RCKWRK(*) at least LENRWK.

!     OUTPUT
!     WL     - Radiation wavelengths for the reactions.
!              WL(I)= 0.  reaction I does not have radiation as
!                         either a reactant or product
!              WL(I)=-A   reaction I has radiation of wavelength A
!                         as a reactant
!              WL(I)=+A   reaction I has radiation of wavelength A
!                         as a product
!              If A = 1.0 then no wavelength information was given;
!                   cgs units - angstrom
!                   Data type - real array
!                   Dimension WL(*) at least II, the total number of
!                   reactions.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION WL(*), ICKWRK(*), RCKWRK(*)

      DO 100 I = 1, NII
         WL(I) = 0.0
  100 CONTINUE
      DO 150 N = 1, NWL
         WL(ICKWRK(IcWL+N-1)) = RCKWRK(NcWL+N-1)
  150 CONTINUE

      RETURN
      END SUBROUTINE CKWL

!     *****************************************************************

      SUBROUTINE CKWC   (T, C, ICKWRK, RCKWRK, WDOT)
!     Returns the molar production rates of the species given the
!     temperature and molar concentrations;  see Eq. (49).

!     INPUT
!     T      - Temperature. [K]
!     C      - Molar concentrations of the species. [mol/cm3]
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     WDOT   - Chemical molar production rates of the species. [mol/cm3/s]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION C(*), ICKWRK(*), RCKWRK(*), WDOT(*)

      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),             &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),         &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,            &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT),&
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),      &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),            &
                   RCKWRK(NcI1))

      DO 25 K = 1, NKK
         RCKWRK(NcK1 + K - 1) = C(K)
         WDOT(K) = 0.0
   25 CONTINUE

      CALL CKRATX (NII, NKK, MXSP, MXTB, T, RCKWRK(NcK1), ICKWRK(IcNS), &
                   ICKWRK(IcNU), ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO),    &
                   NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF), NFAR,&
                   RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),      &
                   RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF),            &
                   RCKWRK(NcKR), RCKWRK(NcI1), RCKWRK(NcI2),            &
                   RCKWRK(NcI3))

      DO 100 N = 1, MXSP
         DO 100 I = 1, NII
            K = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
            NU= ICKWRK(IcNU + (I-1)*MXSP + N - 1)
            IF (K .NE. 0) WDOT(K) = WDOT(K) + NU*    &
                        (RCKWRK(NcI1 + I - 1) - RCKWRK(NcI2 + I - 1))
  100 CONTINUE
      RETURN
      END SUBROUTINE CKWC

!     *****************************************************************

      SUBROUTINE CKWT   (ICKWRK, RCKWRK, WT)
!     Returns the molecular weights of the species

!     INPUT
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     WT     - Molecular weights of the species. [g/mol]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), WT(*)

      DO 100 K = 1, NKK
         WT(K) = RCKWRK(NcWT + K - 1)
  100 CONTINUE

      RETURN
      END SUBROUTINE CKWT

!     *****************************************************************

      SUBROUTINE CKWXP  (P, T, X, ICKWRK, RCKWRK, WDOT)
!     Returns the molar production rates of the species given the
!     pressure, temperature and mole fractions;  see Eq. (49).

!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature. [K]
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     WDOT   - Chemical molar production rates of the species. [mol/cm3/s]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION X(*), ICKWRK(*), RCKWRK(*), WDOT(*)

      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),             &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),         &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,            &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT),&
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),      &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),            &
                   RCKWRK(NcI1))

      CALL CKXTCP (P, T, X, ICKWRK, RCKWRK, RCKWRK(NcK1))

      CALL CKRATX (NII, NKK, MXSP, MXTB, T, RCKWRK(NcK1), ICKWRK(IcNS), &
                   ICKWRK(IcNU), ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO),    &
                   NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF), NFAR,&
                   RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),      &
                   RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF),            &
                   RCKWRK(NcKR), RCKWRK(NcI1), RCKWRK(NcI2),            &
                   RCKWRK(NcI3))

      DO 50 K = 1, NKK
         WDOT(K) = 0.0
   50 CONTINUE

      DO 100 N = 1, MXSP
         DO 100 I = 1, NII
            K = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
            NU= ICKWRK(IcNU + (I-1)*MXSP + N - 1)
            IF (K .NE. 0) &
               WDOT(K) = WDOT(K)+(RCKWRK(NcI1+I-1)-RCKWRK(NcI2+I-1))*NU
  100 CONTINUE
      RETURN
      END SUBROUTINE CKWXP

!     *****************************************************************

      SUBROUTINE CKWXR  (RHO, T, X, ICKWRK, RCKWRK, WDOT)
!     Returns the molar production rates of the species given the
!     mass density, temperature and mole fractions;  see Eq. (49).

!     INPUT
!     RHO    - Mass density. [g/cm3]

!     T      - Temperature. [K]
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     WDOT   - Chemical molar production rates of the species. [mol/cm3/s]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION X(*), ICKWRK(*), RCKWRK(*), WDOT(*)

      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),             &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),         &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,            &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT),&
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),      &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),            &
                   RCKWRK(NcI1))

      CALL CKXTCR (RHO, T, X, ICKWRK, RCKWRK, RCKWRK(NcK1))

      CALL CKRATX (NII, NKK, MXSP, MXTB, T, RCKWRK(NcK1), ICKWRK(IcNS), &
                   ICKWRK(IcNU), ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO),    &
                   NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF), NFAR,&
                   RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),      &
                   RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF),            &
                   RCKWRK(NcKR), RCKWRK(NcI1), RCKWRK(NcI2),            &
                   RCKWRK(NcI3))

      DO 50 K = 1, NKK
         WDOT(K) = 0.0
   50 CONTINUE
      DO 100 N = 1, MXSP
         DO 100 I = 1, NII
            K = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
            NU= ICKWRK(IcNU + (I-1)*MXSP + N - 1)
            IF (K .NE. 0) &
               WDOT(K) = WDOT(K)+(RCKWRK(NcI1+I-1)-RCKWRK(NcI2+I-1))*NU
  100 CONTINUE
      RETURN
      END SUBROUTINE CKWXR

!     *****************************************************************

      SUBROUTINE CKWYP  (P, T, Y, ICKWRK, RCKWRK, WDOT)
!     Returns the molar production rates of the species given the
!     pressure, temperature and mass fractions;  see Eq. (49).

!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature. [K]
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     WDOT   - Chemical molar production rates of the species. [mol/cm3/s]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), Y(*), WDOT(*)

      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),             &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),         &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,            &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT),&
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),      &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),            &
                   RCKWRK(NcI1))

      CALL CKYTCP (P, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK1))

      CALL CKRATX (NII, NKK, MXSP, MXTB, T, RCKWRK(NcK1), ICKWRK(IcNS), &
                   ICKWRK(IcNU), ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO),    &
                   NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF), NFAR,&
                   RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),      &
                   RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF),            &
                   RCKWRK(NcKR), RCKWRK(NcI1), RCKWRK(NcI2),            &
                   RCKWRK(NcI3))

      DO 50 K = 1, NKK
         WDOT(K) = 0.0
   50 CONTINUE
      DO 100 N = 1, MXSP
         DO 100 I = 1, NII
            K = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
            NU= ICKWRK(IcNU + (I-1)*MXSP + N - 1)
            IF (K .NE. 0) &
               WDOT(K) = WDOT(K)+(RCKWRK(NcI1+I-1)-RCKWRK(NcI2+I-1))*NU
  100 CONTINUE
      RETURN
      END SUBROUTINE CKWYP

!     *****************************************************************

      SUBROUTINE CKWYPK  (P, T, Y, RKFT, RKRT, ICKWRK, RCKWRK, WDOT)
!     Returns the molar production rates of the species given the
!     pressure, temperature and mass fractions;  see Eq. (49).

!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature. [K]
!     Y      - Mass fractions of the species.
!                   cgs units - none
!     RKFT   - Forward reaction rates for the reactions
!                   cgs units - depends on the reaction
!     RKRT   - Referse reaction rates for the reactions
!                   cgs units - depends on the reaction
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     WDOT   - Chemical molar production rates of the species. [mol/cm3/s]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), Y(*), RKFT(*), RKRT(*), WDOT(*)
      CALL CKYTCP (P, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK1))
      DO 25 I = 1, NII
         RCKWRK(NcKF + I - 1) = RKFT(I)
         RCKWRK(NcKR + I - 1) = RKRT(I)
   25 CONTINUE

      CALL CKRATX (NII, NKK, MXSP, MXTB, T, RCKWRK(NcK1), ICKWRK(IcNS), &
                   ICKWRK(IcNU), ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO),    &
                   NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF), NFAR,&
                   RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),      &
                   RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF),            &
                   RCKWRK(NcKR), RCKWRK(NcI1), RCKWRK(NcI2),            &
                   RCKWRK(NcI3))

      DO 50 K = 1, NKK
         WDOT(K) = 0.0
   50 CONTINUE
      DO 100 N = 1, MXSP
         DO 100 I = 1, NII
            K = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
            NU= ICKWRK(IcNU + (I-1)*MXSP + N - 1)
            IF (K .NE. 0) &
               WDOT(K) = WDOT(K)+(RCKWRK(NcI1+I-1)-RCKWRK(NcI2+I-1))*NU
  100 CONTINUE
      RETURN
      END SUBROUTINE CKWYPK

!     *****************************************************************

      SUBROUTINE CKWYR  (RHO, T, Y, ICKWRK, RCKWRK, WDOT)
!     Returns the molar production rates of the species given the
!     mass density, temperature and mass fractions;  see Eq. (49).

!     INPUT
!     RHO    - Mass density. [g/cm3]
!     T      - Temperature. [K]
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     WDOT   - Chemical molar production rates of the species. [mol/cm3/s]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION Y(*), RCKWRK(*), WDOT(*), ICKWRK(*)

      CALL CKRATT (RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),             &
                   RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),         &
                   ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,            &
                   ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT),&
                   RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),      &
                   RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),            &
                   RCKWRK(NcI1))

      CALL CKYTCR (RHO, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK1))

      CALL CKRATX (NII, NKK, MXSP, MXTB, T, RCKWRK(NcK1), ICKWRK(IcNS), &
                   ICKWRK(IcNU), ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO),    &
                   NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF), NFAR,&
                   RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),      &
                   RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF),            &
                   RCKWRK(NcKR), RCKWRK(NcI1), RCKWRK(NcI2),            &
                   RCKWRK(NcI3))

      DO 50 K = 1, NKK
         WDOT(K) = 0.d0
   50 CONTINUE
      DO 100 N = 1, MXSP
         DO 100 I = 1, NII
            K = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
            NU= ICKWRK(IcNU + (I-1)*MXSP + N - 1)
            IF (K .NE. 0) &
               WDOT(K) = WDOT(K)+(RCKWRK(NcI1+I-1)-RCKWRK(NcI2+I-1))*NU
  100 CONTINUE
      RETURN
      END SUBROUTINE CKWYR

!     *****************************************************************

      SUBROUTINE CKXNUM (LINE, NEXP, LOUT, NVAL, RVAL, KERR)
!     This subroutine is called to parse a character string, LINE,
!     that is composed of several blank-delimited substrings.
!     Each substring is expected to represent a number, which
!     is converted to entries in the array of real numbers, RVAL(*).
!     NEXP is the number of values expected, and NVAL is the
!     number of values found.  This allows format-free input of
!     numerical data.  For example:

!     input:  LINE    = " 0.170E+14 0 47780.0"
!             NEXP    = 3, the number of values requested
!             LOUT    = 6, a logical unit number on which to write
!                       diagnostic messages.
!     output: NVAL    = 3, the number of values found
!             RVAL(*) = 1.700E+13, 0.000E+00, 4.778E+04
!             KERR    = .FALSE.

!     INPUT
!     LINE   - A character string.
!     NEXP   - Number of real values to be found in character string.
!              If NEXP is negative, then ABS(NEXP) values are
!              expected.  However, it is not an error condition,
!              if less values are found.
!     LOUT   - Output unit for printed diagnostics.

!     OUTPUT
!     NVAL   - Number of real values found in character string.
!     RVAL   - Array of real values found.
!     KERR   - Error flag;  syntax or dimensioning error results
!              in KERR = .TRUE.
!     A '!' will comment out a line, or remainder of the line.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      CHARACTER LINE*(*), ITEMP*80
      DIMENSION RVAL(*), RTEMP(80)
      LOGICAL KERR

!----------Find Comment String (! signifies comment)

      ILEN = IPPLEN(LINE)
      NVAL = 0
      KERR = .FALSE.

      IF (ILEN .LE. 0) RETURN
      IF (ILEN .GT. 80) THEN
         WRITE (LOUT,*)     ' Error in CKXNUM...line length > 80 '
         WRITE (LOUT,'(A)') LINE
         KERR = .TRUE.
         RETURN
      ENDIF

      ITEMP = LINE(:ILEN)
      IF (NEXP .LT. 0) THEN
         CALL IPPARR (ITEMP, -1, NEXP, RTEMP, NVAL, IERR, LOUT)
      ELSE
         CALL IPPARR (ITEMP, -1, -NEXP, RTEMP, NVAL, IERR, LOUT)
         IF (IERR .EQ. 1) THEN
            WRITE (LOUT, *)    ' Syntax errors in CKXNUM...'
            WRITE (LOUT,'(A)') LINE
            KERR = .TRUE.
         ELSEIF (NVAL .NE. NEXP) THEN
            WRITE (LOUT,*) ' Error in CKXNUM...'
            WRITE (LOUT,'(A)') LINE
            KERR = .TRUE.
            WRITE (LOUT,*) NEXP,' values expected, ', &
                           NVAL,' values found.'
         ENDIF
      ENDIF
      IF (NVAL .LE. ABS(NEXP)) THEN
         DO 20 N = 1, NVAL
            RVAL(N) = RTEMP(N)
   20    CONTINUE
      ENDIF

      RETURN
      END SUBROUTINE CKXNUM

!     *****************************************************************

      SUBROUTINE CKXTCP (P, T, X, ICKWRK, RCKWRK, C)
!     Returns the molar concentrations given the pressure,
!     temperature and mole fractions;  see Eq. (10).

!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature. [K]
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     C      - Molar concentrations of the species. [mol/cm3]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), X(*), C(*)

      PRUT = P/(RCKWRK(NcRU)*T)
      DO 150 K = 1, NKK
         C(K) = X(K)*PRUT
150   CONTINUE
      RETURN
      END SUBROUTINE CKXTCP

!     *****************************************************************

      SUBROUTINE CKXTCR (RHO, T, X, ICKWRK, RCKWRK, C)
!     Returns the molar concentrations given the mass density,
!     temperature and mole fractions;  see Eq. (11).

!     INPUT
!     RHO    - Mass density. [g/cm3]
!     T      - Temperature. [K]
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     C      - Molar concentrations of the species. [mol/cm3]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), X(*), C(*)
      SUM = 0.0
      DO 100 K = 1, NKK
         SUM = SUM + X(K)*RCKWRK(NcWT + K - 1)
  100 CONTINUE
      RHOW = RHO / SUM
      DO 200 K = 1, NKK
         C(K) = X(K)*RHOW
200   CONTINUE
      RETURN
      END SUBROUTINE CKXTCR

!     *****************************************************************

      SUBROUTINE CKXTY  (X, ICKWRK, RCKWRK, Y)
!     Returns the mass fractions given the mole fractions;
!     see Eq. (9).

!     INPUT
!     X      - Mole fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     Y      - Mass fractions of the species.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), X(*), Y(*)
      SUM = 0.0
      DO 100 K = 1, NKK
         SUM = SUM + X(K)*RCKWRK(NcWT + K - 1)
  100 CONTINUE

      DO 200 K = 1, NKK
         Y(K) = X(K)*RCKWRK(NcWT + K - 1)/SUM
200   CONTINUE
      RETURN
      END SUBROUTINE CKXTY

!     *****************************************************************

      SUBROUTINE CKYTCP (P, T, Y, ICKWRK, RCKWRK, C)
!     Returns the molar concentrations given the pressure,
!     temperature and mass fractions;  see Eq. (7).

!     INPUT
!     P      - Pressure. [dyn/cm2]
!     T      - Temperature. [K]
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     C      - Molar concentrations of the species. [mol/cm3]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), Y(*), C(*)
      SUMYOW = 0.0
      DO 150 K = 1, NKK
         SUMYOW = SUMYOW + Y(K)/RCKWRK(NcWT + K - 1)
150   CONTINUE
      SUMYOW = SUMYOW*T*RCKWRK(NcRU)
      DO 200 K = 1, NKK
         C(K) = P*Y(K)/(SUMYOW*RCKWRK(NcWT + K - 1))
200   CONTINUE
      RETURN
      END SUBROUTINE CKYTCP

!     *****************************************************************

      SUBROUTINE CKYTCR (RHO,T, Y, ICKWRK, RCKWRK, C)
!     Returns the molar concentrations given the mass density,
!     temperature and mass fractions;  see Eq. (8).

!     INPUT
!     RHO    - Mass density. [g/cm3]
!     T      - Temperature. [K]
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     C      - Molar concentrations of the species. [mol/cm3]
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), Y(*), C(*)

      DO 150 K = 1, NKK
         C(K) = RHO*Y(K)/RCKWRK(NcWT + K - 1)
150   CONTINUE
      RETURN
      END SUBROUTINE CKYTCR

!     *****************************************************************

      SUBROUTINE CKYTX  (Y, ICKWRK, RCKWRK, X)
!     Returns the mole fractions given the mass fractions;  see Eq. (6).

!     INPUT
!     Y      - Mass fractions of the species.
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     X      - Mole fractions of the species.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), Y(*), X(*)

      SUMYOW = 0.0
      DO 150 K = 1, NKK
         SUMYOW = SUMYOW + Y(K)/RCKWRK(NcWT + K - 1)
150   CONTINUE
      DO 200 K = 1, NKK
         X(K) = Y(K)/(SUMYOW*RCKWRK(NcWT + K - 1))
200   CONTINUE
      RETURN
      END SUBROUTINE CKYTX

!     *****************************************************************

      SUBROUTINE IPPARI(STRING, ICARD, NEXPEC, IVAL, NFOUND, IERR, LOUT)
!   BEGIN PROLOGUE  IPPARI
!   REFER TO  IPGETI
!   DATE WRITTEN  850625   (YYMMDD)
!   REVISION DATE 851725   (YYMMDD)
!   CATEGORY NO.  J3.,J4.,M2.
!   KEYWORDS  PARSE
!   AUTHOR  CLARK,G.L.,GROUP C-3 LOS ALAMOS NAT'L LAB
!   PURPOSE  Parses integer variables from a character variable.  Called
!            by IPGETI, the IOPAK routine used for interactive input.
!   DESCRIPTION
!
!-----------------------------------------------------------------------
!  IPPARI may be used for parsing an input record that contains integer
!  values, but was read into a character variable instead of directly
!  into integer variables.
!  The following benefits are gained by this approach:
!    - specification of only certain elements of the array is allowed,
!      thus letting the others retain default values
!    - variable numbers of values may be input in a record, up to a
!      specified maximum
!    - control remains with the calling program in case of an input
!      error
!    - diagnostics may be printed by IPPARI to indicate the nature
!      of input errors
!
!   The contents of STRING on input indicate which elements of IVAL
!   are to be changed from their entry values, and values to which
!   they should be changed on exit.  Commas and blanks serve as
!   delimiters, but multiple blanks are treated as a single delimeter.
!   Thus, an input record such as:
!     '   1,   2,,40000   , ,60'
!   is interpreted as the following set of instructions by IPGETR:
!
!     (1) set IVAL(1) = 1
!     (2) set IVAL(2) = 2
!     (3) leave IVAL(3) unchanged
!     (4) set IVAL(4) = 40000
!     (5) leave IVAL(5) unchanged
!     (6) set IVAL(6) = 60
!
!   IPPARI will print diagnostics on the default output device, if
!   desired.
!
!   IPPARI is part of IOPAK, and is written in ANSI FORTRAN 77
!
!   Examples:
!
!      Assume IVAL = (0, 0, 0) and NEXPEC = 3 on entry:
!
!   input string           IVAL on exit            IERR    NFOUND
!   -------------          ----------------------  ----    ------
!  '  2 ,   3 45 '         (2, 3, 45)                0       3
!  '2.15,,3'               (2, 0, 3)                 1       0
!  '3X, 25, 2'             (0, 0, 0)                 1       0
!  '10000'                 (10000, 0, 0)             2       1
!
!      Assume IVAL = (0, 0, 0, 0) and NEXPEC = -4 on entry:
!
!   input string           IVAL on exit            IERR    NFOUND
!   -------------          ----------------------  ----    ------
!  '1, 2'                  (1, 2)                    0       2
!  ',,37  400'             (0, 0, 37, 400)           0       4
!  ' 1,,-3,,5'             (1, 0, -3, 0)             3       4
!
!  arguments: (I=input,O=output)
!  -----------------------------
!  STRING (I) - the character string to be parsed.
!
!  ICARD  (I) - data statement number, and error processing flag
!         < 0 : no error messages printed
!         = 0 : print error messages, but not ICARD
!         > 0 : print error messages, and ICARD
!
!  NEXPEC (I) - number of real variables expected to be input.  If
!         < 0, the number is unknown, and any number of values
!         between 0 and abs(nexpec) may be input.  (see NFOUND)
!
!  PROMPT (I) - prompting string, character type.  A question
!         mark will be added to form the prompt at the screen.
!
!  IVAL (I,O) - the integer value or values to be modified.  On entry,
!       the values are printed as defaults.  The formal parameter
!       corresponding to IVAL must be dimensioned at least NEXPEC
!       in the calling program if NEXPEC > 1.
!
!  NFOUND (O) - the number of real values represented in STRING,
!         only in the case that there were as many or less than
!         NEXPEC.
!
!  IERR (O) - error flag:
!       = 0 if no errors found
!       = 1 syntax errors or illegal values found
!       = 2 for too few values found (NFOUND < NEXPEC)
!       = 3 for too many values found (NFOUND > NEXPEC)
!   REFERENCES  (NONE)
!   ROUTINES CALLED  IFIRCH,ILASCH
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
!
!
      CHARACTER STRING*(*), ITEMP*80
      DIMENSION IVAL(*)
      CHARACTER *8 FMT(14)
      LOGICAL OKINCR
!
!   FIRST EXECUTABLE STATEMENT  IPPARI
      IERR   = 0
      NFOUND = 0
      NEXP = IABS(NEXPEC)
      IE = ILASCH(STRING)
      IF (IE .EQ. 0) GO TO 500
      NC = 1
!
!--- OKINCR is a flag that indicates it's OK to increment
!--- NFOUND, the index of the array into which the value
!--- should be read.  It is set false when a space follows
!--- an integer value substring, to keep incrementing from
!--- occurring if a comma should be encountered before the
!--- next value.
!
      OKINCR = .TRUE.
!
!--- begin overall loop on characters in string
!
100   CONTINUE

      IF (STRING(NC:NC) .EQ. ',') THEN
         IF (OKINCR .OR. NC .EQ. IE) THEN
            NFOUND = NFOUND + 1
         ELSE
            OKINCR = .TRUE.
         ENDIF

         GO TO 450
      ENDIF
      IF (STRING(NC:NC) .EQ. ' ') GO TO 450

!--- first good character (non-delimeter) found - now find
!--- last good character

      IBS = NC
160   CONTINUE
      NC = NC + 1
      IF (NC .GT. IE) GO TO 180
      IF (STRING(NC:NC) .EQ. ' ')THEN
         OKINCR = .FALSE.
      ELSEIF (STRING(NC:NC) .EQ. ',')THEN
         OKINCR = .TRUE.
      ELSE
         GO TO 160
      ENDIF
!
!--- end of substring found - read value into integer array
!
180   CONTINUE
      NFOUND = NFOUND + 1
      IF (NFOUND .GT. NEXP) THEN
         IERR = 3
         GO TO 500
      ENDIF

      IES = NC - 1
      NCH = IES - IBS + 1
      DATA FMT/' (I1)', ' (I2)', ' (I3)', ' (I4)', ' (I5)', &
         ' (I6)', ' (I7)', ' (I8)', ' (I9)', '(I10)',       &
         '(I11)', '(I12)', '(I13)', '(I14)'/
      ITEMP = ' '
      ITEMP = STRING(IBS:IES)
      READ (ITEMP(1:NCH), FMT(NCH), ERR = 400) IVAL(NFOUND)
      GO TO 450
400   CONTINUE
      IERR = 1
      GO TO 510
450   CONTINUE
      NC = NC + 1
      IF (NC .LE. IE) GO TO 100

500   CONTINUE
      IF (NEXPEC .GT. 0 .AND. NFOUND .LT. NEXP) IERR = 2
510   CONTINUE

      IF (IERR .EQ. 0 .OR. ICARD .LT. 0)RETURN
      IF (ICARD .NE. 0) WRITE(LOUT,'(A,I3)') &
         '!! ERROR IN DATA STATEMENT NUMBER', ICARD
      IF (IERR .EQ. 1) WRITE(LOUT,'(A)')'SYNTAX ERROR, OR ILLEGAL VALUE'
      IF (IERR .EQ. 2) WRITE(LOUT,'(A,I2, A, I2)')          &
         ' TOO FEW DATA ITEMS.  NUMBER FOUND = ' , NFOUND,  &
         '  NUMBER EXPECTED = ', NEXPEC
      IF (IERR .EQ. 3) WRITE(LOUT,'(A,I2)') &
         ' TOO MANY DATA ITEMS.  NUMBER EXPECTED = ', NEXPEC
      END SUBROUTINE IPPARI

!     *****************************************************************

      SUBROUTINE IPPARR (STRING,ICARD,NEXPEC,RVAL,NFOUND,IERR,LOUT)
!   BEGIN PROLOGUE  IPPARR
!   REFER TO  IPGETR
!   DATE WRITTEN  850625   (YYMMDD)
!   REVISION DATE 851625   (YYMMDD)
!   CATEGORY NO.  J3.,J4.,M2.
!   KEYWORDS  PARSE
!   AUTHOR  CLARK,G.L.,GROUP C-3 LOS ALAMOS NAT'L LAB
!   PURPOSE  Parses real variables from a character variable.  Called
!            by IPGETR, the IOPAK routine used for interactive input.
!   DESCRIPTION

!-----------------------------------------------------------------------
!  IPPARR may be used for parsing an input record that contains real
!  values, but was read into a character variable instead of directly
!  into real variables.
!  The following benefits are gained by this approach:
!    - specification of only certain elements of the array is allowed,
!      thus letting the others retain default values
!    - variable numbers of values may be input in a record, up to a
!      specified maximum
!    - control remains with the calling program in case of an input
!      error
!    - diagnostics may be printed by IPPARR to indicate the nature
!      of input errors

!   The contents of STRING on input indicate which elements of RVAL
!   are to be changed from their entry values, and values to which
!   they should be changed on exit.  Commas and blanks serve as
!   delimiters, but multiple blanks are treated as a single delimeter.
!   Thus, an input record such as:
!     '   1.,   2,,4.e-5   , ,6.e-6'
!   is interpreted as the following set of instructions by IPGETR:

!     (1) set RVAL(1) = 1.0
!     (2) set RVAL(2) = 2.0
!     (3) leave RVAL(3) unchanged
!     (4) set RVAL(4) = 4.0E-05
!     (5) leave RVAL(5) unchanged
!     (6) set RVAL(6) = 6.0E-06

!   IPPARR will print diagnostics on the default output device, if
!   desired.

!   IPPARR is part of IOPAK, and is written in ANSI FORTRAN 77

!   Examples:

!      Assume RVAL = (0., 0., 0.) and NEXPEC = 3 on entry:

!   input string           RVAL on exit            IERR    NFOUND
!   -------------          ----------------------  ----    ------
!  '  2.34e-3,  3 45.1'    (2.34E-03, 3.0, 45.1)     0       3
!  '2,,3.-5'               (2.0, 0.0, 3.0E-05)       0       3
!  ',1.4,0.028E4'          (0.0, 1.4, 280.0)         0       3
!  '1.0, 2.a4, 3.0'        (1.0, 0.0, 0.0)           1       1
!  '1.0'                   (1.0, 0.0, 0.0)           2       1

!      Assume RVAL = (0.,0.,0.,0.) and NEXPEC = -4 on entry:

!   input string           RVAL on exit            IERR    NFOUND
!   -------------          ----------------------  ----    ------
!  '1.,2.'                 (1.0, 2.0)                0       2
!  ',,3  4.0'              (0.0, 0.0, 3.0, 4.0)      0       4
!  '1,,3,,5.0'             (0.0, 0.0, 3.0, 0.0)      3       4

!  arguments: (I=input,O=output)
!  -----------------------------
!  STRING (I) - the character string to be parsed.

!  ICARD  (I) - data statement number, and error processing flag
!         < 0 : no error messages printed
!         = 0 : print error messages, but not ICARD
!         > 0 : print error messages, and ICARD

!  NEXPEC (I) - number of real variables expected to be input.  If
!         < 0, the number is unknown, and any number of values
!         between 0 and abs(nexpec) may be input.  (see NFOUND)

!  PROMPT (I) - prompting string, character type.  A question
!         mark will be added to form the prompt at the screen.

!  RVAL (I,O) - the real value or values to be modified.  On entry,
!       the values are printed as defaults.  The formal parameter
!       corresponding to RVAL must be dimensioned at least NEXPEC
!       in the calling program if NEXPEC > 1.

!  NFOUND (O) - the number of real values represented in STRING,
!         only in the case that there were as many or less than
!         NEXPEC.

!  IERR (O) - error flag:
!       = 0 if no errors found
!       = 1 syntax errors or illegal values found
!       = 2 for too few values found (NFOUND < NEXPEC)
!       = 3 for too many values found (NFOUND > NEXPEC)
!-----------------------------------------------------------------------

!   REFERENCES  (NONE)
!   ROUTINES CALLED  IFIRCH,ILASCH
!   END PROLOGUE  IPPARR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

      CHARACTER STRING*(*), ITEMP*256, FMT_DYNAMIC*32
      DIMENSION RVAL(*)
      LOGICAL OKINCR

!   FIRST EXECUTABLE STATEMENT  IPPARR
      IERR   = 0
      NFOUND = 0
      NEXP = IABS(NEXPEC)
      IE = ILASCH(STRING)
      IF (IE .EQ. 0) GO TO 500
      NC = 1

!--- OKINCR is a flag that indicates it's OK to increment
!--- NFOUND, the index of the array into which the value
!--- should be read.  It is set negative when a space follows
!--- a real value substring, to keep incrementing from
!--- occurring if a comma should be encountered before the
!--- next value.

      OKINCR = .TRUE.

!--- begin overall loop on characters in string

100   CONTINUE

      IF (STRING(NC:NC) .EQ. ',') THEN
         IF (OKINCR) THEN
            NFOUND = NFOUND + 1
         ELSE
            OKINCR = .TRUE.
         ENDIF

         GO TO 450
      ENDIF
      IF (STRING(NC:NC) .EQ. ' ') GO TO 450

!--- first good character (non-delimeter) found - now find
!--- last good character

      IBS = NC
160   CONTINUE
      NC = NC + 1
      IF (NC .GT. IE) GO TO 180
      IF (STRING(NC:NC) .EQ. ' ')THEN
         OKINCR = .FALSE.
      ELSEIF (STRING(NC:NC) .EQ. ',')THEN
         OKINCR = .TRUE.
      ELSE
         GO TO 160
      ENDIF

!--- end of substring found - read value into real array

180   CONTINUE
      NFOUND = NFOUND + 1
      IF (NFOUND .GT. NEXP) THEN
         IERR = 3
         GO TO 500
      ENDIF

      IES = NC - 1
      NCH = IES - IBS + 1
      IF (NCH .GT. LEN(ITEMP)) THEN
         IERR = 1
         GO TO 510
      ENDIF
      ITEMP = ' '
      ITEMP = STRING(IBS:IES)
!     The original IOPAK table only defined E formats through E16.0.
!     yaml2ck legitimately emits round-trip decimal tokens such as
!     0.010000000000000004 (20 characters), which previously indexed
!     beyond that table and made valid PLOG lines fail as syntax errors.
!     Construct the identical Ew.0 descriptor for the actual width.
      WRITE (FMT_DYNAMIC,'("(E",I0,".0)")') NCH
      READ (ITEMP(:NCH), FMT_DYNAMIC, ERR = 400) RVAL(NFOUND)
      GO TO 450
400   CONTINUE
      IERR = 1
      GO TO 510
450   CONTINUE
      NC = NC + 1
      IF (NC .LE. IE) GO TO 100

500   CONTINUE
      IF (NEXPEC .GT. 0 .AND. NFOUND .LT. NEXP) IERR = 2
510   CONTINUE

      IF (IERR .EQ. 0 .OR. ICARD .LT. 0) RETURN
      IF (ICARD .NE. 0) WRITE(LOUT,'(A,I3)') &
         '!! ERROR IN DATA STATEMENT NUMBER', ICARD
      IF (IERR .EQ. 1) WRITE(LOUT,'(A)')'SYNTAX ERROR, OR ILLEGAL VALUE'
      IF (IERR .EQ. 2) WRITE(LOUT,'(A,I2, A, I2)')         &
         ' TOO FEW DATA ITEMS.  NUMBER FOUND = ' , NFOUND, &
         '  NUMBER EXPECTED = ', NEXPEC
      IF (IERR .EQ. 3) WRITE(LOUT,'(A,I2)') &
         ' TOO MANY DATA ITEMS.  NUMBER EXPECTED = ', NEXPEC
      END SUBROUTINE IPPARR

!     *****************************************************************

      INTEGER FUNCTION IFIRCH   (STRING)
!   BEGIN PROLOGUE  IFIRCH
!   DATE WRITTEN   850626
!   REVISION DATE  850626
!   CATEGORY NO.  M4.
!   KEYWORDS  CHARACTER STRINGS,SIGNIFICANT CHARACTERS
!   AUTHOR  CLARK,G.L.,GROUP C-3 LOS ALAMOS NAT'L LAB
!   PURPOSE  Determines first significant (non-blank) character
!            in character variable
!   DESCRIPTION

!-----------------------------------------------------------------------
!  IFIRCH locates the first non-blank character in a string of
!  arbitrary length.  If no characters are found, IFIRCH is set = 0.
!  When used with the companion routine ILASCH, the length of a string
!  can be determined, and/or a concatenated substring containing the
!  significant characters produced.
!-----------------------------------------------------------------------

!   REFERENCES  (NONE)
!   ROUTINES CALLED  (NONE)
!   END PROLOGUE IFIRCH
    IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

    CHARACTER* (*)STRING

!   FIRST EXECUTABLE STATEMENT IFIRCH
      NLOOP = LEN(STRING)

      IF (NLOOP.EQ.0 .OR. STRING.EQ.' ') THEN
         IFIRCH = 0
         RETURN
      ENDIF

      DO 100 I = 1, NLOOP
         IF (STRING(I:I) .NE. ' ') GO TO 120
100   CONTINUE

      IFIRCH = 0
      RETURN
120   CONTINUE
      IFIRCH = I
      END FUNCTION IFIRCH

!     *****************************************************************

      INTEGER FUNCTION ILASCH   (STRING)
!   BEGIN PROLOGUE  ILASCH
!   DATE WRITTEN   850626
!   REVISION DATE  850626
!   CATEGORY NO.  M4.
!   KEYWORDS  CHARACTER STRINGS,SIGNIFICANT CHARACTERS
!   AUTHOR  CLARK,G.L.,GROUP C-3 LOS ALAMOS NAT'L LAB
!   PURPOSE  Determines last significant (non-blank) character
!            in character variable
!   DESCRIPTION

!-----------------------------------------------------------------------
!  IFIRCH locates the last non-blank character in a string of
!  arbitrary length.  If no characters are found, ILASCH is set = 0.
!  When used with the companion routine IFIRCH, the length of a string
!  can be determined, and/or a concatenated substring containing the
!  significant characters produced.
!  Note that the FORTRAN intrinsic function LEN returns the length
!  of a character string as declared, rather than as filled.  The
!  declared length includes leading and trailing blanks, and thus is
!  not useful in generating 'significant' substrings.
!-----------------------------------------------------------------------

!   REFERENCES  (NONE)
!   ROUTINES CALLED  (NONE)
!   END PROLOGUE IFIRCH
    IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

    CHARACTER*(*) STRING

!   FIRST EXECUTABLE STATEMENT ILASCH
      NLOOP = LEN(STRING)
      IF (NLOOP.EQ.0 .OR. STRING.EQ.' ') THEN
         ILASCH = 0
         RETURN
      ENDIF

      DO 100 I = NLOOP, 1, -1
         ILASCH = I
         IF (STRING(I:I) .NE. ' ') RETURN
100   CONTINUE

      END FUNCTION ILASCH

!     *****************************************************************

      FUNCTION IPPLEN (LINE)
!     Returns the effective length of a character string, i.e.,
!     the index of the last character before an exclamation mark (!)
!     indicating a comment.
!     INPUT
!     LINE  - A character string.

!     OUTPUT
!     IPPLEN - The effective length of the character string.
!  END PROLOGUE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
      CHARACTER LINE*(*)

      IN = IFIRCH(LINE)
      IF (IN.EQ.0) THEN
         IPPLEN = 0
!     Do not join the two IF clauses to avoid array error (LINE(0:0))
      ELSEIF (LINE(IN:IN) .EQ. '!') THEN
         IPPLEN = 0
      ELSE
         IN = INDEX(LINE,'!')
         IF (IN .EQ. 0) THEN
            IPPLEN = ILASCH(LINE)
         ELSE
            IPPLEN = ILASCH(LINE(:IN-1))
         ENDIF
      ENDIF
      RETURN
      END FUNCTION IPPLEN

!     *****************************************************************

!     Added by Hai Wang, Penn State Univ. Oct.1, 1993
      SUBROUTINE CKFALX (I, LOUT, ICKWRK, RCKWRK, RD)
!     Put the perturbation factor of the ith reaction with fall-off
!     behavior.
!     INPUT
!     I      - Reaction number
!     LOUT   - a logical unit number on which to write diagnostic
!              messages
!     ICKWRK - Array of integer workspace.
!     RCKWRK - Array of real work space.

!     OUTPUT
!     RD     - perturbation factor
!  END PROLOGUE
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*)

      DO 200 N = 1 , NFAL
         IF (I .EQ. ICKWRK(IcFL+N-1)) GOTO 220
  200 CONTINUE
      WRITE (LOUT, 210) I
  210 FORMAT(/'WARNING: did not find reaction ',I3, &
              ' being defined as fall-off')

  220 CONTINUE
      NI = NcFL + (N-1) * NFAR
      RCKWRK(NI) = RD

      RETURN
      END SUBROUTINE CKFALX

!     *****************************************************************

!----------------------------------------------------------------------!
      FUNCTION UPCASE(ISTR, ILEN) result(upper)

      integer,          intent(in) :: ILEN
      CHARACTER(len=*), intent(in) :: ISTR
      character(len=len(istr))     :: upper
      INTEGER                      :: J, JJ, N
      character, dimension(26)     :: &
         LCASE = ['a','b','c','d','e','f','g','h','i','j','k','l','m', &
                  'n','o','p','q','r','s','t','u','v','w','x','y','z'],&
         UCASE = ['A','B','C','D','E','F','G','H','I','J','K','L','M', &
                  'N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
!
      upper = ' '
      upper = ISTR(:ILEN)
      JJ = MIN (LEN(upper), LEN(ISTR), ILEN)
      DO J = 1, JJ
         DO N = 1,26
            IF (ISTR(J:J) .EQ. LCASE(N)) upper(J:J) = UCASE(N)
         end do
      end do

      END FUNCTION UPCASE
!----------------------------------------------------------------------!


      end module chemkinII

!     *****************************************************************
!     **                                                             **
!     **   PLOG COLLECTION (cklink v2, stage 1 PLOG plumbing)        **
!     **                                                             **
!     **   Dynamic, growable parse-time storage for PLOG             **
!     **   (pressure-dependent Arrhenius) reactions. The classic     **
!     **   CHEMKIN interpreter threads every feature through the     **
!     **   CKAUXL/CKINTP argument lists and fixed-size arrays        **
!     **   (KORD(3000,9000) etc.); PLOG instead lives here as a      **
!     **   self-contained module so no giant argument chain or       **
!     **   static MAXPLOG x IDIM array is added (see plan.md         **
!     **   "推奨する内部データ構造").                                **
!     **                                                             **
!     **   CKAUXL calls plog_add_line() once per `PLOG / P A b E /`  **
!     **   line; CKINTP calls plog_finalize() before writing cklink  **
!     **   and reads the packed arrays out. Pressures must be         **
!     **   non-decreasing per reaction. Adjacent equal-pressure       **
!     **   entries are retained as the multiple Arrhenius terms that  **
!     **   CHEMKIN/Cantera sum at that pressure.                      **
!     *****************************************************************

      module plog_collect

      use working_precision, only: dp
      implicit none
      public

!     Dialect selector retained for on-disk/API compatibility. Both
!     accepted modes use standard grouped same-pressure PLOG semantics.
      integer, parameter :: PLOG_STRICT_CHEMKIN     = 0
      integer, parameter :: PLOG_PERMISSIVE_GROUPED = 1
      integer :: plog_dialect = PLOG_STRICT_CHEMKIN

!     Growable flat list of PLOG lines, in the order encountered.
!     Each entry i: pcl_reac(i)=reaction index, pcl_logP(i)=ln(P[Pa]),
!     pcl_A/pcl_b/pcl_EoverR(i)=Arrhenius params (A in input units,
!     b dimensionless, E/R in K). Grown geometrically to avoid a fixed
!     upper bound.
      integer :: plog_nlines = 0
      integer, dimension(:), allocatable :: pcl_reac
      real (dp), dimension(:), allocatable :: pcl_logP
      real (dp), dimension(:), allocatable :: pcl_A
      real (dp), dimension(:), allocatable :: pcl_b
      real (dp), dimension(:), allocatable :: pcl_EoverR

!     Packed output (filled by plog_finalize), mirrors the reacpar
!     packed layout so cklink v2 write is a direct dump:
!       plog_n_reactions              : # distinct PLOG reactions
!       plog_reaction(1:nr_plog)      : their global reaction indices (asc)
!       plog_node_ptr(0:nr_plog)      : node range per reaction (CSR)
!       plog_logP_out(1:n_nodes)      : ln(P[Pa]) per entry (non-decreasing)
!       plog_A_out/b_out/EoverR_out   : one Arrhenius set per entry
!     Adjacent entries with equal pressure form one pressure node and
!     are summed by the evaluator. The cklink v2 representation remains
!     backward compatible because unique-pressure files are unchanged.
      integer :: plog_n_reactions = 0
      integer :: plog_n_nodes     = 0
      integer, dimension(:), allocatable :: plog_reaction
      integer, dimension(:), allocatable :: plog_node_ptr
      real (dp), dimension(:), allocatable :: plog_logP_out
      real (dp), dimension(:), allocatable :: plog_A_out
      real (dp), dimension(:), allocatable :: plog_b_out
      real (dp), dimension(:), allocatable :: plog_EoverR_out

      contains

!     -------------------------------------------------------------
!     plog_reset: clear all state (call at start of a CKINTP run).
!     -------------------------------------------------------------
      subroutine plog_reset
      implicit none
      plog_nlines = 0
      plog_n_reactions = 0
      plog_n_nodes = 0
      if (allocated(pcl_reac))        deallocate(pcl_reac)
      if (allocated(pcl_logP))        deallocate(pcl_logP)
      if (allocated(pcl_A))           deallocate(pcl_A)
      if (allocated(pcl_b))           deallocate(pcl_b)
      if (allocated(pcl_EoverR))      deallocate(pcl_EoverR)
      if (allocated(plog_reaction))   deallocate(plog_reaction)
      if (allocated(plog_node_ptr))   deallocate(plog_node_ptr)
      if (allocated(plog_logP_out))   deallocate(plog_logP_out)
      if (allocated(plog_A_out))      deallocate(plog_A_out)
      if (allocated(plog_b_out))      deallocate(plog_b_out)
      if (allocated(plog_EoverR_out)) deallocate(plog_EoverR_out)
      end subroutine plog_reset

!     -------------------------------------------------------------
!     plog_grow: ensure the flat lists hold at least nreq entries,
!     preserving contents (geometric growth).
!     -------------------------------------------------------------
      subroutine plog_grow(nreq)
      implicit none
      integer, intent(in) :: nreq
      integer :: newcap, oldcap
      integer, dimension(:), allocatable :: itmp
      real (dp), dimension(:), allocatable :: rtmp
      if (allocated(pcl_reac)) then
         oldcap = size(pcl_reac)
      else
         oldcap = 0
      endif
      if (nreq <= oldcap) return
      newcap = max(16, oldcap*2)
      do while (newcap < nreq)
         newcap = newcap*2
      end do
!     integer field
      allocate(itmp(newcap)); itmp = 0
      if (oldcap > 0) itmp(1:oldcap) = pcl_reac(1:oldcap)
      call move_alloc(itmp, pcl_reac)
!     real fields
      allocate(rtmp(newcap)); rtmp = 0.0_dp
      if (oldcap > 0) rtmp(1:oldcap) = pcl_logP(1:oldcap)
      call move_alloc(rtmp, pcl_logP)
      allocate(rtmp(newcap)); rtmp = 0.0_dp
      if (oldcap > 0) rtmp(1:oldcap) = pcl_A(1:oldcap)
      call move_alloc(rtmp, pcl_A)
      allocate(rtmp(newcap)); rtmp = 0.0_dp
      if (oldcap > 0) rtmp(1:oldcap) = pcl_b(1:oldcap)
      call move_alloc(rtmp, pcl_b)
      allocate(rtmp(newcap)); rtmp = 0.0_dp
      if (oldcap > 0) rtmp(1:oldcap) = pcl_EoverR(1:oldcap)
      call move_alloc(rtmp, pcl_EoverR)
      end subroutine plog_grow

!     -------------------------------------------------------------
!     plog_add_line: record one `PLOG / P A b E /` line for reaction
!     `ireac`. P is in atm (CHEMKIN convention), A/b are Arrhenius
!     params in the mechanism's declared units, `eraw` is the activation
!     energy exactly as read (same raw units as PAR(3,*) before EFAC).
!     It is stored raw here and converted to E/R [K] later by
!     plog_apply_efac(), so PLOG E goes through the SAME unit handling
!     (CAL/KCAL/JOUL/KJOU/KELV) as every other Arrhenius E — no separate
!     conversion path that could drift (cf. the KJOU 4x bug).
!     Sets kerr=.true. on invalid input (non-positive pressure); the
!     pressure-ordering checks happen in plog_finalize
!     where all lines of a reaction are visible together.
!     -------------------------------------------------------------
      subroutine plog_add_line(ireac, pressure_atm, aval, bval, eraw, kerr, lout)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
      integer, intent(in) :: ireac, lout
      real (dp), intent(in) :: pressure_atm, aval, bval, eraw
      logical, intent(inout) :: kerr
      real (dp), parameter :: atm_to_pa = 101325.0_dp
      if (.not. (pressure_atm > 0.0_dp)) then
!        Fail-closed: non-positive PLOG pressure is a hard input error;
!        error stop (not STOP/return) so it is detectable by exit code.
         write(*   ,'(A,I6,A,1PE12.4,A)')                              &
            ' ERROR...PLOG pressure for reaction ', ireac,             &
            ' is not positive (', pressure_atm, ' atm)'
         write(lout,'(A,I6,A,1PE12.4,A)')                              &
            ' ERROR...PLOG pressure for reaction ', ireac,             &
            ' is not positive (', pressure_atm, ' atm)'
         kerr = .true.
         error stop 1
      endif
      if (.not. ieee_is_finite(aval) .or. .not. (aval > 0.0_dp)) then
         write(*   ,'(A,I6,A,1PE12.4)')                               &
            ' ERROR...PLOG A factor for reaction ', ireac,            &
            ' must be finite and positive; got ', aval
         write(lout,'(A,I6,A,1PE12.4)')                               &
            ' ERROR...PLOG A factor for reaction ', ireac,            &
            ' must be finite and positive; got ', aval
         kerr = .true.
         error stop 1
      endif
      call plog_grow(plog_nlines+1)
      plog_nlines = plog_nlines + 1
      pcl_reac(plog_nlines)   = ireac
      pcl_logP(plog_nlines)   = log(pressure_atm*atm_to_pa)
      pcl_A(plog_nlines)      = aval
      pcl_b(plog_nlines)      = bval
      pcl_EoverR(plog_nlines) = eraw   ! raw E; converted by plog_apply_efac
      end subroutine plog_add_line

!     -------------------------------------------------------------
!     plog_apply_efac: multiply the stored (raw) E of every PLOG line
!     belonging to reaction `ireac` by `efac`, converting E -> E/R [K].
!     Called from CPREAC once per reaction, right where PAR(3,II) is
!     converted, so PLOG activation energies use the identical EFAC.
!     Idempotency is the caller's responsibility (CPREAC runs once per
!     reaction); we only touch lines whose reaction index matches.
!     -------------------------------------------------------------
      subroutine plog_apply_efac(ireac, efac)
      implicit none
      integer, intent(in) :: ireac
      real (dp), intent(in) :: efac
      integer :: i
      if (plog_nlines <= 0) return
      do i = 1, plog_nlines
         if (pcl_reac(i) == ireac) pcl_EoverR(i) = pcl_EoverR(i) * efac
      end do
      end subroutine plog_apply_efac

!     -------------------------------------------------------------
!     plog_apply_afac: apply the same global MOLECULES -> MOLES
!     conversion used for the main Arrhenius A factor. PLOG replaces
!     the main rate expression, so every pressure-node A must receive
!     the identical reaction-order-dependent conversion.
!     -------------------------------------------------------------
      subroutine plog_apply_afac(ireac, afac)
      implicit none
      integer, intent(in) :: ireac
      real (dp), intent(in) :: afac
      integer :: i
      if (plog_nlines <= 0) return
      do i = 1, plog_nlines
         if (pcl_reac(i) == ireac) pcl_A(i) = pcl_A(i) * afac
      end do
      end subroutine plog_apply_afac

!     -------------------------------------------------------------
!     plog_finalize: turn the flat line list into packed, per-reaction
!     CSR arrays. Within each reaction, pressure entries must be
!     non-decreasing. Equal adjacent entries are legal grouped terms.
!     Lines are already grouped by reaction because CKAUXL processes one
!     reaction's aux lines consecutively; we still group defensively.
!     Sets kerr=.true. on a strict-dialect violation.
!     -------------------------------------------------------------
      subroutine plog_finalize(kerr, lout)
      implicit none
      logical, intent(inout) :: kerr
      integer, intent(in) :: lout
      integer :: i, j, r, nreac, node, ir
      integer, dimension(:), allocatable :: uniq
      real (dp), parameter :: logp_tol = 1.0e-9_dp

      plog_n_reactions = 0
      plog_n_nodes = 0
      if (plog_nlines <= 0) then
!        No PLOG reactions: leave empty packed arrays (size 0) so the
!        v2 writer emits a count of 0 and nothing else changes.
         allocate(plog_reaction(0), plog_node_ptr(0:0))
         plog_node_ptr(0) = 0
         allocate(plog_logP_out(0), plog_A_out(0), plog_b_out(0), &
                  plog_EoverR_out(0))
         return
      endif

!     Distinct reaction indices, in first-seen order (CKINTP numbers
!     reactions ascending, so this is ascending too).
      allocate(uniq(plog_nlines))
      nreac = 0
      do i = 1, plog_nlines
         ir = pcl_reac(i)
         if (nreac == 0) then
            nreac = 1
            uniq(1) = ir
         elseif (ir /= uniq(nreac)) then
!           Guard: a reaction's PLOG lines must be contiguous. If a
!           previously-seen reaction reappears, the mechanism interleaved
!           PLOG lines in a way we don't support -> refuse.
            do j = 1, nreac
               if (uniq(j) == ir) then
!                 Fail-closed: a plain STOP here is exit 0 (reads as
!                 success to a test harness); use error stop so an
!                 invalid PLOG mechanism is detectable by exit code.
                  write(*   ,'(A,I6,A)')                                &
                     ' ERROR...PLOG lines for reaction ', ir,           &
                     ' are not contiguous in the input'
                  write(lout,'(A,I6,A)')                                &
                     ' ERROR...PLOG lines for reaction ', ir,           &
                     ' are not contiguous in the input'
                  kerr = .true.
                  error stop 1
               endif
            end do
            nreac = nreac + 1
            uniq(nreac) = ir
         endif
      end do

      plog_n_reactions = nreac
      plog_n_nodes = plog_nlines
      allocate(plog_reaction(nreac), plog_node_ptr(0:nreac))
      allocate(plog_logP_out(plog_nlines), plog_A_out(plog_nlines), &
               plog_b_out(plog_nlines), plog_EoverR_out(plog_nlines))
      plog_node_ptr(0) = 0

      node = 0
      do r = 1, nreac
         plog_reaction(r) = uniq(r)
         do i = 1, plog_nlines
            if (pcl_reac(i) == uniq(r)) then
               node = node + 1
!              Standard PLOG: pressure entries are non-decreasing.
!              Equal adjacent pressures are multiple Arrhenius terms at
!              one pressure node and are summed during rate evaluation.
               if (node > plog_node_ptr(r-1)+1) then
                  if (pcl_logP(i) < plog_logP_out(node-1) - logp_tol) then
!                    Fail-closed (see note above): error stop, not STOP.
                     write(*   ,'(A,I6,A)')                             &
                        ' ERROR...PLOG pressures for reaction ',        &
                        uniq(r),                                        &
                        ' are out of order (pressures must be'//         &
                        ' non-decreasing)'
                     write(lout,'(A,I6,A)')                             &
                        ' ERROR...PLOG pressures for reaction ',        &
                        uniq(r),                                        &
                        ' are out of order (pressures must be'//         &
                        ' non-decreasing)'
                     kerr = .true.
                     error stop 1
                  endif
               endif
               plog_logP_out(node)    = pcl_logP(i)
               plog_A_out(node)       = pcl_A(i)
               plog_b_out(node)       = pcl_b(i)
               plog_EoverR_out(node)  = pcl_EoverR(i)
            endif
         end do
         plog_node_ptr(r) = node
      end do

      deallocate(uniq)
      end subroutine plog_finalize

      end module plog_collect

!     *****************************************************************
!     *****************************************************************
!     *****************************************************************












!     *****************************************************************
!     **                                                             **
!     **                     KIVA4 - CHEMISTRY                       **
!     **                                                             **
!     **   ChemkinII interpreter and runtime in FORTRAN 2003 format  **
!     **           Taken from version 4.2, 14/9/1993                 **
!     **                                                             **
!     **                                                             **
!     **   Modified by: Federico Perini                              **
!     **   Last update: thursday, 01/12/2011                         **
!     **                                                             **
!     *****************************************************************


      module chemkinII_interpreter

      implicit none
      public

!     *****************************************************************
      contains


      SUBROUTINE CKINTP
!
!----------------------------------------------------------------------!
!     VERSION 3.9
!
!=======================================================================
!
!     CKINTP interprets a formatted ASCII representation of a
!     chemical reaction mechanism and creates the binary file LINK
!     required by CHEMKIN.  CKINTP is dimensioned as follows:
!
!     MDIM = maximum number of elements in a problem;             (10)
!     KDIM = maximum number of species in a problem;             (100)
!     MAXTP= maximum number of temperatures used to fit            (3)
!            thermodynamic properties of species
!     NPC  = number of polynomial coefficients to fits             (5)
!     NPCP2= number of fit coefficients for a temperature range    (7)
!     IDIM = maximum number of reactions in a mechanism;         (500)
!     NPAR = number of Arrhenius parameters in a reaction;         (3)
!     NLAR = number of Landau-Teller parameters in a reaction;     (2)
!     NFAR = number of fall-off parameters in a reaction;          (8)
!     MAXSP= maximum number of species in a reaction               (6)
!     MAXTB= maximum number of third bodies for a reaction        (10)
!     LSYM = character string length of element and species names (18)
!
!     User input is read from LIN (Unit15), a thermodynamic database
!     is read from LTHRM (Unit17), printed output is assigned to LOUT
!     (Unit16), and binary data is written to LINC (Unit25).
!
!     REQUIRED ELEMENT INPUT: (Subroutine CKCHAR)          (DIMENSION)
!
!        The word 'ELEMENTS' followed by a list of element
!        names, terminated by the word 'END';
!
!        The resulting element data stored in LINK is:
!        MM       - integer number of elements found
!        ENAME(*) - CHARACTER*(*) array of element names        (MDIM)
!        AWT(*)   - real array of atomic weights;               (MDIM)
!                   default atomic weights are those on
!                   atomic weight charts; if an element
!                   is not on the periodic chart, or if
!                   it is desirable to alter its atomic
!                   weight, this value must be included
!                   after the element name, enclosed by
!                   slashed, i.e., D/2.014/
!
!     REQUIRED SPECIES INPUT: (Subroutine CKCHAR)
!
!        The word 'SPECIES' followed by a list of species
!        names, terminated by the word 'END';
!
!        The resulting species data stored in LINK is:
!        KK       - integer number of species found
!        KNAME(*) - CHARACTER*(*) array of species names        (KDIM)
!
!     OPTIONAL THERMODYNAMIC DATA: (Subroutine CKTHRM)
!     (If this feature is not used, thermodynamic properties are
!     obtained from a CHEMKIN database.)  The format for this option
!     is the word 'THERMO' followed by any number of 4-line data sets:
!
!     Line 1: species name, optional comments, elemental composition,
!             phase, T(low), T(high), T(mid), additional elemental
!             composition, card number (col. 80);
!             format(A10,A14,4(A2,I3),A1,E10.0,E10.0,E8.0,(A2,I3),I1)
!     Line 2: coefficients a(1--5) for upper temperature range,
!             card number (col. 80);
!             format(5(e15.0),I1)
!     Line 3: coefficients a(6--7) for upper temperature range,
!             coefficients a(1--3) for lower temperature range,
!             card number (col. 80);
!             format(5(e15.0),I1)
!     Line 4: coefficients a(4--7) for lower temperature range,
!             card number (col. 80);
!             format(4(e15.0),I1)
!
!     End of THERMO data is indicated by 'END' line or new keyword.
!
!        The resulting thermodynamic data stored in LINK are:
!        WTM(*)   - real array of molecular weights             (KDIM)
!        KNCF(*,*)- integer composition of species         (MDIM,KDIM)
!        KPHSE(*) - integer phase of a species;                 (KDIM)
!                   -1(solid), 0(gas), +1(liquid).
!        KCHRG(*) - ionic charge of a species;                  (KDIM)
!                   = 0 except in presence/absence of electrons
!                   = +n in absence of n electrons
!                   = -n in presence of n electons
!        NCHRG    - integer number of species with KCHRG<>0
!        NT(*)    - array of number of temperatures used        (KDIM)
!                   in fits
!        T(*,*)   - array of temperatures used in fits    (MAXTP,KDIM)
!        A(N,L,K) - Thermodynamic properties for      (NPC+2,NTR,KDIM)
!                   species K consists of polynomial
!                   coefficients for fits to
!                   CP/R = SUM (A(N,L,K)*Temperature**(N-1), N=1,NPC+2)
!                          where  T(L,K) <= Temperature < T(L+1,K),
!                   and,
!                   N=NPC+1 is formation enthalpy HO/R = A(NPC+1,L,K),
!                   N=NPC+2 is formation entropy  SO/R = A(NPC+2,L,K)
!
!     OPTIONAL REACTION INPUT:
!     Reaction data is input after all ELEMENT, SPECIES and THERMO
!     data in the following format:
!
!     1) (Subroutine CKREAC)
!        The first line contains the keyword 'REACTIONS' and an
!        optional description of units:
!
!           'MOLES' - (default), pre-exponential units are moles-sec-K;
!           'MOLECULES' - pre-exponential units are molecules and
!                         will be converted to moles.
!           'KELVINS' - activation energies are Kelvins, else the
!                       activation energies are converted to Kelvins;
!           'CAL/MOLE' - (default), activation energies are cal/mole;
!           'KCAL/MOLE' - activation energies are Kcal/mole;
!           'JOULES/MOLE' - activation energies are joules/mole;
!           'KJOULES/MOLE' - activation energies are Kjoules/mole.
!
!        A description of each reaction is expected to follow.
!        Required format for a reaction is a list of '+'-delimited
!        reactants, followed by a list of '+'-delimited reactants,
!        each preceded by its stoichiometric coefficient if greater
!        than 1;  separating the reactants from the products is a '='
!        if reversible reaction, else a '=>'.  Following the reaction
!        string on the same line are the space-delimited Arrhenius
!        coefficients.
!
!        If the reaction contains a third body, this is indicated by
!        by the presence of an 'M' as a reactant or product or both,
!        and enhancement factors for third-bodies may be defined on
!        additional lines as described in (2).
!
!        If the reaction contains a radiation wavelength, this is
!        indicated by the presence of an 'HV' either as a reactant
!        or as a product.  Unless otherwise defined on additional
!        lines as described in (2), the value of the wavelength is
!        -1.0 if a reactant or +1.0 if a product.
!
!        If the reaction is a fall-off reaction, this is indicated
!        either by a '(+M)' or a '(+KNAME(K))', and there must be
!        additional lines as described in (2) to define fall-off
!        parameters.
!
!    2)  (Subroutine CKAUXL)
!        Additional information for a reaction is given on lines
!        immediately following the reaction description; this data
!        will consist of a 'keyword' to denote the type of data,
!        followed by a '/', then the required parameters for the
!        keyword, followed by another '/'.  There may be more than
!        one keyword per line, and there may be any number of lines.
!        The keywords and required parameters are as follows:
!
!        KNAME(K)/efficiency value/ - species (K) is an enhanced
!                third body in the reaction
!        HV/wavelength/ - radiation wavelength parameter
!        LT/val1 val2/ - Landau-Teller coefficients
!        LOW/val1 val2 val3/ - low fall-off parameters
!        TROE/val1 val2 val3 val4/ - Troe fall-off parameters;
!                                    if val4 is omitted, a default
!                                    parameter will be used
!        SRI/val1 val2 val3 val4/ - SRI fall-off parameters;
!                                   if val4 is omitted, a default
!                                   parameter will be used
!           (it is an error to have both LT and Fall-off defined)
!        REV/par1 par2 par3/ - reverse parameters given
!        RLT/val1 val2/ - Landau-Teller coefficients for reverse
!           (it is an error if REV given and not RLT)
!        EIM/val1/ - Electon-impct reaction; val1 is the integer
!                    temperature dependence flag
!        JAN/val1...val9/ - coefficients for electron reactions in
!                    the form of Jannev, Langer & Post:
!                    k = SUM[an (lnT)^n]
!        FIT1/val1...val4/ - additional exponential terms for
!                    temperature powers > 1, e.g.,
!                    k = A T^B exp [ SUM (valn/T^n) ]
!        EXCI/val1/ - excitation reaction for energy los only;
!                     val1 is the energy loss per event in eV
!
!     The end of all reaction data is indicated by an 'END' card or
!     <eof>.
!
!     Resulting reaction data stored in LINC are:
!       II        - integer number of reactions found
!       PAR(*,*)  - array of real Arrhenius coefficients   (NPAR,IDIM)
!       NSPEC(*)  - total number of species in a reaction       (IDIM)
!                   if NSPEC < 0, reaction is irreversible
!       NREAC(*)  - number of reactants only                    (IDIM)
!       NUNK(*,*) - array of species indices for reaction (MAXSP,IDIM)
!       NU(*,*)   - array of stoichiometric coefficients  (MAXSP,IDIM)
!                   of species in a reaction, negative=reactant,
!                   positive=product
!
!       NWL       - number of reactions with radiation wavelength
!       IWL(*)    - the NWL reaction indices                    (IDIM)
!       WL(*)     - real radiation wavelengths                  (IDIM)
!
!       NTHB      - number of reactions with third bodies
!       ITHB      - the NTHB reaction indices                   (IDIM)
!       NTBS(*)   - total number of enhanced species for NTHB   (IDIM)
!       NKTB(*,*) - species indices of enhanced species   (MAXTB,IDIM)
!       AIK(*,*)  - enhancement factors                   (MAXTB,IDIM)
!
!       NFAL      - number of fall-off reactions
!       IFAL(*)   - the NFAL reaction indices                   (IDIM)
!       KFAL(*)   - integer species number for which
!                   concentrations are a factor in fall-off
!                   calculation
!       IFOP(*)   - integer fall-off type number                (IDIM)
!                   = 0 if fall-off reaction is found
!                   = 1 for Lindemann form
!                   = 2 for 6-parameter Troe form
!                   = 3 for 7-parameter Troe form
!                   = 4 for SRI form
!       PFAL(*,*) - fall-off parameters                    (NFAR,IDIM)
!
!       NLAN      - number of reactions with Landau-Teller
!       ILAN(*)   - the NLAN reaction indices                   (IDIM)
!       PLAN(*,*) - Landau-Teller parameters               (NLAR,IDIM)
!
!       NREV      - number of reactions with reverse parameters
!       IREV(*)   - the NREV reaction indices                   (IDIM)
!       RPAR(*,*) - parameters                             (NPAR,IDIM)
!
!       NRLT      - number of reactions with reverse parameters
!                   and Landau-Teller parameters
!       IRLT(*)   - the NRLT reaction indices                   (IDIM)
!       RLAN(*,*) - reverse Teller-Laudauer parameters     (NLAR,IDIM)
!       NEIM      - number of reactions with electron impact
!       IEIM(*)   - the NEIM reaction indices                   (IDIM)
!       ITDEP(*)  - the NEIM temperature dependence flags       (IDIM)
!
!       NJAN      - number of Jannev, Langer, Evans & Post reactions
!       IJAN(*)   - the NJAN reaction indices                   (IDIM)
!       PJAN(*,*) - coefficients for the NJAN reactions    (NJAR,IDIM)
!
!       NFT1      - number of reactions using fit #1
!       IFT1(*)   - the NFT1 reaction indices                   (IDIM)
!       PFT1(*,*) - additional exponential terms for fit#1 (NF1R,IDIM)
!
!       NEX!      - number of excitation reactions
!       IEXC(*)   - the NEXC reaction indices                   (IDIM)
!       PEXC(*)   - energy loss per event in units of eV        (IDIM)
!
!       NRNU      - number of reactions having real stoichiometry
!       IRNU(*)   - the NRNU reaction indices                   (IDIM)
!       RNU(*,*)  - matrix of real stoich. coefficients   (MAXSP,IDIM)
!
!       NORD      - number of reactions with modified species orders
!       IORD(*)   - the NORD reaction indices                   (IDIM)
!       KORD(*,*) - matrix of species indices whose order (MAXORD,IDIM)
!                   is modified
!       RORD(*,*) - matrix of species order values        (MAXORD,IDIM)
!----------------------------------------------------------------------!

      USE chemkinII,    only: IPPLEN, ILASCH, UPCASE, IFIRCH, IPPARR,&
                              chemdat
      USE chemistry_string_limits, only: species_name_len,            &
                                          mechanism_line_len
!ck2015
      USE chemistry_setup, only: mechdir
!     PLOG collection + cklink v2 packed arrays (stage 1 PLOG plumbing)
      USE plog_collect, only: plog_reset, plog_finalize,               &
                              plog_n_reactions, plog_n_nodes,          &
                              plog_reaction, plog_node_ptr,            &
                              plog_logP_out, plog_A_out, plog_b_out,   &
                              plog_EoverR_out

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
!
      PARAMETER (MDIM=10,KDIM=3000, MKDIM=MDIM*KDIM,IDIM=9000,        &
                 LSYM=species_name_len,                                &
                 NPAR=3, NPIDIM=IDIM*NPAR, NPC=5, NPCP2=NPC+2, MAXTP=3,&
                 NTR=MAXTP-1, NKTDIM=NTR*NPCP2*KDIM, MAXSP=8,MAXTB=10,&
                 NLAR=2, NSIDIM=MAXSP*IDIM, NTIDIM=MAXTB*IDIM,&
                 NLIDIM=NLAR*IDIM, NFAR=8, NFIDIM=NFAR*IDIM,&
                 NTDIM=KDIM*MAXTP, NIDIM=11*IDIM, LIN=40, LOUT=41,&
                 LTHRM=42, LINC=43, CKMIN=1.0E-3, MAXORD=KDIM,&
                 NOIDIM=MAXORD*IDIM)
!
!     Modern yaml2ck output is not restricted to the historical
!     80-column CHEMKIN card width (real mechanisms contain lines up to
!     at least 116 characters). Keep thermo-card readers at 80 columns,
!     but do not truncate mechanism/reaction and auxiliary-data lines.
      CHARACTER KNAME(KDIM)*(LSYM), ENAME(MDIM)*(LSYM),               &
                SUB(80)*mechanism_line_len, KEY(5)*4,                 &
                LINE*mechanism_line_len, IUNITS*80, AUNITS*4,         &
                EUNITS*4,&
!               FILETD holds trim(mechdir)//"therm.dat"; mechdir is
!               len=256 (chemistry_setup), so 80 silently truncated long
!               paths and therm.dat could not be opened. Match mechdir.
!     VERS/PREC are cklink header metadata with a legacy 16-character
!     record contract; they are independent of the species-name width.
                VERS*16, PREC*16, FILETD*265 !UPCASE*4
!     cklink v2 schema: a leading magic + integer schema version so the
!     reader can positively identify the format and refuse older/newer
!     files, replacing the fragile VERS-string check. Bump CK_SCHEMA on
!     any on-disk layout change.
      CHARACTER(len=8), PARAMETER :: CK_MAGIC = 'SCLKv2  '
      INTEGER, PARAMETER          :: CK_SCHEMA = 3
!
      DIMENSION AWT(MDIM), KNCF(MDIM,KDIM), WTM(KDIM), KPHSE(KDIM),&
                KCHRG(KDIM), A(NPCP2,NTR,KDIM), T(MAXTP,KDIM), NT(KDIM),&
                NSPEC(IDIM), NREAC(IDIM), NU(MAXSP,IDIM),&
                NUNK(MAXSP,IDIM), PAR(NPAR,IDIM), IDUP(IDIM),IREV(IDIM),&
                RPAR(NPAR,IDIM), ILAN(IDIM), PLAN(NLAR,IDIM),&
                IRLT(IDIM), RLAN(NLAR,IDIM), IWL(IDIM),  WL(IDIM),&
                IFAL(IDIM), IFOP(IDIM), KFAL(IDIM), PFAL(NFAR,IDIM),&
                ITHB(IDIM),NTBS(IDIM),AIK(MAXTB,IDIM),NKTB(MAXTB,IDIM),&
                IRNU(IDIM), RNU(MAXSP,IDIM), IORD(IDIM),&
                KORD(MAXORD,IDIM), RORD(MAXORD,IDIM)
      DIMENSION VALUE(5)
!
      LOGICAL KERR, THERMO, ITHRM(KDIM)
!
      PARAMETER (NJAR=9, NF1R=4, NJIDIM=NJAR*IDIM, NF1IDIM=NF1R*IDIM)
      DIMENSION IEIM(IDIM), ITDEP(IDIM), IJAN(IDIM), PJAN(NJAR,IDIM),&
                IFT1(IDIM), PFT1(NF1R,IDIM), IEXC(IDIM), PEXC(IDIM)
      DATA NEIM,NJAN,NFT1,NEXC/4*0/, IEIM/IDIM*0/, ITDEP/IDIM*0/,&
           IJAN/IDIM*0/, IFT1/IDIM*0/, PJAN/NJIDIM*0.0/,&
           PFT1/NF1IDIM*0.0/, PEXC/IDIM*0.0/
!
!     Initialize variables
!
      DATA KEY/'ELEM','SPEC','THER','REAC','END'/, KERR/.FALSE./,&
           ITASK,NCHRG,MM,KK,II,NLAN,NFAL,NTHB,NREV,NRLT,NWL,&
           NRNU,NORD/13*0/,&
           ENAME,AWT/MDIM*' ',MDIM*0.0/, THERMO/.TRUE./,&
           T/NTDIM*-1.0/, KNAME,WTM,NT,KPHSE,KCHRG,ITHRM&
           /KDIM*' ', KDIM*0.0, KDIM*3, KDIM*0, KDIM*0, KDIM*.FALSE./,&
           WL,IFOP,NTBS,IDUP /IDIM*0.0, IDIM*-1, IDIM*0, IDIM*0/,&
           NSPEC,NREAC,IREV,ILAN,IRLT,IWL,IFAL,KFAL,ITHB,IRNU,IORD&
           /NIDIM*0/
!
      DATA NUNK,NU/NSIDIM*0, NSIDIM*0/, NKTB,AIK/NTIDIM*0,NTIDIM*-1.0/
      DATA RNU/NSIDIM*0.0/, KORD/NOIDIM*0/, RORD/NOIDIM*0.0/
      DATA PAR,RPAR/NPIDIM*0.0, NPIDIM*0.0/
      DATA PLAN,RLAN/NLIDIM*0.0, NLIDIM*0.0/
      DATA PFAL/NFIDIM*0.0/, KNCF/MKDIM*0.0/, A/NKTDIM*0.0/


!ck2015      FILETD='therm.dat'
      FILETD=trim(mechdir)//"therm.dat"

!
!ck2015      OPEN (LOUT, FORM='FORMATTED', STATUS='UNKNOWN', FILE=chemdat)
      OPEN (LOUT, FORM='FORMATTED', STATUS='UNKNOWN', &
            FILE=trim(mechdir)//chemdat)
!
      VERS = '3.1'
!     Clear any PLOG state left from a previous CKINTP call in the same
!     process (cklink v2, stage 1 PLOG plumbing).
      CALL plog_reset
      WRITE  (LOUT, 15) VERS(:4)
   15 FORMAT (/ &
      ' CHEMKIN INTERPRETER OUTPUT: CHEMKIN-II Version ',A,' Aug. 1994' &
      /'                              DOUBLE PRECISION'/)
      PREC = 'DOUBLE'
!
!        START OF MECHANISM INTERPRETATION
!
!ck2015      OPEN (LIN, FORM='FORMATTED', STATUS='OLD', FILE='chem.inp',&
!ck2015            ERR=11111)
      OPEN (LIN, FORM='FORMATTED', STATUS='OLD', &
            FILE=trim(mechdir)//"chem.inp",ERR=11111)
      READ (LIN,'(A)',END=11111)
!
      REWIND (LIN)
  100 CONTINUE
      LINE = ' '
      READ (LIN,'(A)',END=5000) LINE
  105 CONTINUE
      ILEN = IPPLEN(LINE)
      IF (ILEN .EQ. 0) GO TO 100
!
      CALL CKISUB (LINE(:ILEN), SUB, NSUB)
!
!        IS THERE A KEYWORD?
!
      CALL CKCOMP ( UPCASE(SUB(1), 4) , KEY, 5, NKEY)
      IF (NKEY .GT. 0) ITASK = 0
!
      IF (NKEY.EQ.1 .OR. NKEY.EQ.2) THEN
!
!        ELEMENT OR SPECIES DATA
!
         ITASK = NKEY
         IF (NSUB .EQ. 1) GO TO 100
!
         DO 25 N = 2, NSUB
            SUB(N-1) = ' '
            SUB(N-1) = SUB(N)
   25    CONTINUE
         NSUB = NSUB-1
!
      ELSEIF (NKEY .EQ. 3) THEN
!
!        THERMODYNAMIC DATA
!
         IF (NSUB .GT. 1) THEN
            IF ( UPCASE(SUB(2), 3) .EQ. 'ALL') THEN
               THERMO = .FALSE.
               READ (LIN,'(A)') LINE
               CALL IPPARR (LINE, -1, 3, VALUE, NVAL, IER, LOUT)
               IF (NVAL .NE. 3 .OR. IER.NE.0) THEN
                  KERR = .TRUE.
                  WRITE (LOUT, 333)
               ELSE
                  TLO = VALUE(1)
                  TMID = VALUE(2)
                  THI = VALUE(3)
               ENDIF
            ENDIF
         ELSE
!
!           USE THERMODYNAMIC DATABASE FOR DEFAULT TLO,TMID,THI
            OPEN (LTHRM, FORM='FORMATTED', STATUS='OLD',&
                         FILE=FILETD, ERR=22222)
!
  311       CONTINUE
            READ (LTHRM,'(A)',END=22222) LINE
            IF (IPPLEN(LINE).LE.0 .OR. INDEX(LINE,'THERMO').GT.0&
                .OR. INDEX(LINE,'thermo').GT.0) GO TO 311
!
            CALL IPPARR (LINE, -1, 3, VALUE, NVAL, IER, LOUT)
            IF (NVAL .NE. 3 .OR. IER.NE.0) THEN
               KERR = .TRUE.
               WRITE (LOUT, 333)
            ELSE
               TLO = VALUE(1)
               TMID = VALUE(2)
               THI = VALUE(3)
            ENDIF
            CLOSE (LTHRM)
         ENDIF
!
         CALL CKTHRM (LIN, MDIM, ENAME, MM, AWT, KNAME, KK, KNCF,&
                      KPHSE, KCHRG, WTM, MAXTP, NT, NTR, TLO, TMID,&
                      THI, T, NPCP2, A, ITHRM, KERR, LOUT, LINE)
!
         IF (.NOT. THERMO)&
            CALL CKPRNT (MDIM, MAXTP, MM, ENAME, KK, KNAME, WTM, KPHSE,&
                         KCHRG, NT, T, TLO, TMID, THI, KNCF, ITHRM,&
                         LOUT, KERR)
         I1 = IFIRCH(LINE)
         IF (UPCASE(LINE(I1:), 4) .EQ. 'REAC') GO TO 105
!
      ELSEIF (NKEY .EQ. 4) THEN
!
         ITASK = 4
!        START OF REACTIONS; ARE UNITS SPECIFIED?
         CALL CKUNIT (LINE(:ILEN), AUNITS, EUNITS, IUNITS)
!
         IF (THERMO) THEN
!
!           THERMODYNAMIC DATA
            OPEN (LTHRM, FORM='FORMATTED', STATUS='OLD',&
                         FILE=FILETD, ERR=22222)
  312       CONTINUE
            READ (LTHRM,'(A)',END=22222) LINE
            IF (IPPLEN(LINE).LE.0 .OR. INDEX(LINE,'THERM').GT.0&
                .OR. INDEX(LINE,'therm').GT.0) GO TO 312
!
            CALL IPPARR (LINE, -1, 3, VALUE, NVAL, IER, LOUT)
            IF (NVAL .NE. 3 .OR. IER.NE.0) THEN
               KERR = .TRUE.
               WRITE (LOUT, 333)
            ELSE
               TLO = VALUE(1)
               TMID = VALUE(2)
               THI = VALUE(3)
            ENDIF
            CALL CKTHRM (LTHRM, MDIM, ENAME, MM, AWT, KNAME, KK, KNCF,&
                         KPHSE, KCHRG, WTM, MAXTP, NT, NTR, TLO, TMID,&
                         THI, T, NPCP2, A, ITHRM, KERR, LOUT, LINE)
            CALL CKPRNT (MDIM, MAXTP, MM, ENAME, KK, KNAME, WTM, KPHSE,&
                         KCHRG, NT, T, TLO, TMID, THI, KNCF, ITHRM,&
                         LOUT, KERR)
            THERMO = .FALSE.
            CLOSE (LTHRM)
         ENDIF
!
         WRITE (LOUT, 1800)
         GO TO 100
      ENDIF
!
      IF (ITASK .EQ. 1) THEN
!
!        ELEMENT DATA
!
         IF (MM .EQ. 0) THEN
            WRITE (LOUT, 200)
            WRITE (LOUT, 300)
            WRITE (LOUT, 200)
         ENDIF
!
         IF (NSUB .GT. 0) THEN
            M1 = MM +1
            CALL CKCHAR (SUB, NSUB, MDIM, ENAME, AWT, MM, KERR, LOUT)
            DO 110 M = M1, MM
               IF (AWT(M) .LE. 0) CALL CKAWTM (ENAME(M), AWT(M))
               WRITE (LOUT, 400) M,ENAME(M)(:4),AWT(M)
               IF (AWT(M) .LE. 0) THEN
                  KERR = .TRUE.
                  WRITE (LOUT, 1000) ENAME(M)
               ENDIF
  110       CONTINUE
         ENDIF
!
      ELSEIF (ITASK .EQ. 2) THEN
!
!        PROCESS SPECIES DATA
!
         IF (KK .EQ. 0) WRITE (LOUT, 200)
         IF (NSUB .GT. 0)&
         CALL CKCHAR (SUB, NSUB, KDIM, KNAME, WTM, KK, KERR, LOUT)
!
      ELSEIF (ITASK .EQ. 4) THEN
!
!        PROCESS REACTION DATA
!
         IND = 0
         DO 120 N = 1, NSUB
            IND = MAX(IND, INDEX(SUB(N),'/'))
            IF (UPCASE(SUB(N), 3) .EQ. 'DUP') IND = MAX(IND,1)
  120    CONTINUE
         IF (IND .GT. 0) THEN
!
!           AUXILIARY REACTION DATA
!
            CALL CKAUXL (SUB, NSUB, II, KK, KNAME, LOUT, MAXSP, NPAR,&
                         NSPEC, NTHB, ITHB, NTBS, MAXTB, NKTB, AIK,&
                         NFAL, IFAL, IDUP, NFAR, PFAL, IFOP, NLAN,&
                         ILAN, NLAR, PLAN, NREV, IREV, RPAR, NRLT, IRLT,&
                         RLAN, NWL, IWL, WL, KERR, NORD, IORD, MAXORD,&
                         KORD, RORD, NUNK, NU, NRNU, IRNU, RNU,&
                         NEIM, IEIM, ITDEP, NJAN, IJAN, NJAR, PJAN,&
                         NFT1, IFT1, NF1R, PFT1, NEXC, IEXC, PEXC)
!
         ELSE
!
!           THIS IS A REACTION STRING
!
            IF (II .LT. IDIM) THEN
!
               IF (II .GT. 0)&
!
!              CHECK PREVIOUS REACTION FOR COMPLETENESS
!
               CALL CPREAC (II, MAXSP, NSPEC, NPAR, PAR, RPAR,&
                            AUNITS, EUNITS, NREAC, NUNK, NU, KCHRG,&
                            MDIM, MM, KNCF, IDUP, NFAL, IFAL, KFAL,&
                            NFAR, PFAL, IFOP, NREV, IREV, NTHB, ITHB,&
                            NLAN, ILAN, NRLT, IRLT, KERR, LOUT, NRNU,&
                            IRNU, RNU, CKMIN)
!
!              NEW REACTION
!
               II = II+1
               CALL CKREAC (LINE(:ILEN), II, KK, KNAME, LOUT, MAXSP,&
                            NSPEC, NREAC, NUNK, NU, NPAR, PAR,&
                            NTHB, ITHB, NFAL, IFAL, KFAL, NWL,&
                            IWL, WL, NRNU, IRNU, RNU, KERR)
!
            ELSE
!              Fail immediately. Continuing with II pinned at IDIM makes
!              auxiliary records from later reactions attach to reaction
!              IDIM and can replace this primary error with a false PLOG
!              ordering diagnostic.
               WRITE (*,1070) IDIM
               WRITE (LOUT,1070) IDIM
               CLOSE (LIN)
               CLOSE (LOUT)
               ERROR STOP 1
            ENDIF
!
         ENDIF
      ENDIF
      GO TO 100
!
 5000 CONTINUE
!
!     END OF INPUT
!
      IF (II .GT. 0) THEN
!
!              CHECK FINAL REACTION FOR COMPLETENESS
!
          CALL CPREAC (II, MAXSP, NSPEC, NPAR, PAR, RPAR, AUNITS,&
                       EUNITS, NREAC, NUNK, NU, KCHRG, MDIM, MM,&
                       KNCF, IDUP, NFAL, IFAL, KFAL, NFAR, PFAL, IFOP,&
                       NREV, IREV, NTHB, ITHB, NLAN, ILAN, NRLT,&
                       IRLT, KERR, LOUT, NRNU, IRNU, RNU, CKMIN)
!
!              CHECK REACTIONS DECLARED AS DUPLICATES
!
         DO 500 I = 1, II
            IF (IDUP(I) .LT. 0) THEN
               KERR = .TRUE.
               WRITE (LOUT, 1095) I
            ENDIF
  500    CONTINUE
!
         WRITE (LOUT, '(/1X,A)') ' NOTE: '//IUNITS(:ILASCH(IUNITS))
!
      ELSEIF (THERMO) THEN
!
!        THERE WAS NO REACTION DATA, MAKE SURE SPECIES DATA IS COMPLETE
         OPEN (LTHRM, FORM='FORMATTED', STATUS='OLD',&
                      FILE=FILETD, ERR=22222)
!
  313    CONTINUE
         READ (LTHRM,'(A)',END=22222) LINE
         IF (IPPLEN(LINE).LE.0 .OR. INDEX(LINE,'THERM').GT.0&
             .OR. INDEX(LINE,'therm').GT.0) GO TO 313
!
         CALL IPPARR (LINE, -1, 3, VALUE, NVAL, IER, LOUT)
         IF (NVAL .NE. 3 .OR. IER.NE.0) THEN
            KERR = .TRUE.
            WRITE (LOUT, 333)
         ELSE
            TLO = VALUE(1)
            TMID = VALUE(2)
            THI = VALUE(3)
         ENDIF
         CALL CKTHRM (LTHRM, MDIM, ENAME, MM, AWT, KNAME, KK, KNCF,&
                      KPHSE, KCHRG, WTM, MAXTP, NT, NTR, TLO, TMID,&
                      THI, T, NPCP2, A, ITHRM, KERR, LOUT, LINE)
         CALL CKPRNT (MDIM, MAXTP, MM, ENAME, KK, KNAME, WTM, KPHSE,&
                      KCHRG, NT, T, TLO, TMID, THI, KNCF, ITHRM,&
                      LOUT, KERR)
         CLOSE  (LTHRM)
      ENDIF
!
      CLOSE (LIN)
!
!     Preserve the primary parser diagnostic. Do not run PLOG validation
!     and do not create a linking file when an earlier parse check failed.
      IF (KERR) THEN
         WRITE (LOUT, '(//A)')&
         ' ERROR...MECHANISM PARSING FAILED; NO LINKING FILE WAS WRITTEN'
         CLOSE (LOUT)
         ERROR STOP 1
      ENDIF
!
!     Validate and pack PLOG data before opening cklink. plog_finalize uses
!     ERROR STOP for malformed data, so no partial link exists on failure.
      CALL plog_finalize(KERR, LOUT)
      IF (KERR) THEN
         WRITE (LOUT, '(//A)')&
         ' ERROR...PLOG VALIDATION FAILED; NO LINKING FILE WAS WRITTEN'
         CLOSE (LOUT)
         ERROR STOP 1
      ENDIF
!
!     Create the link only after all parse-time validation succeeds.
      OPEN (LINC, FORM='UNFORMATTED', STATUS='REPLACE',&
                  FILE=trim(mechdir)//"cklink")
      REWIND LINC
!     v2 leading record: magic + integer schema version. Positively
!     identifies the format for the reader (SCcklink) — supersedes the
!     fragile VERS-string check.
      WRITE (LINC) CK_MAGIC, CK_SCHEMA
      WRITE (LINC) VERS, PREC, KERR

!
      DO 1150 K = 1, KK
         IF (KCHRG(K) .NE. 0) NCHRG = NCHRG+1
 1150 CONTINUE
!
      LENICK = 1 + (3 + MM)*KK + (2 + 2*MAXSP)*II + NLAN + NRLT&
                 + 3*NFAL + (2 + MAXTB)*NTHB + NREV + NWL + NRNU&
                 + NORD*(1 + MAXORD) + 2*NEIM + NJAN + NFT1&
                 + NEXC
!
      LENCCK = MM + KK
!
      LENRCK = 3 + MM + KK*(5 + MAXTP + NTR*NPCP2) + II*7 + NREV&
                 + NPAR*(II + NREV) + NLAR*(NLAN + NRLT)&
                 + NFAR*NFAL + MAXTB*NTHB + NWL + NRNU*MAXSP&
                 + NORD*MAXORD + NJAR*NJAN + NF1R*NFT1 + NEXC
!
      WRITE (LINC) LENICK, LENRCK, LENCCK, MM, KK, II, MAXSP,&
                   MAXTB, MAXTP, NPC, NPAR, NLAR, NFAR, NREV, NFAL,&
                   NTHB, NLAN, NRLT, NWL, NCHRG, NEIM, NJAR, NJAN,&
                   NF1R, NFT1, NEXC, NRNU, NORD, MAXORD, CKMIN
      WRITE (LINC) (ENAME(M), AWT(M), M = 1, MM)
      WRITE (LINC) (KNAME(K), (KNCF(M,K),M=1,MM), KPHSE(K),&
                    KCHRG(K), WTM(K), NT(K), (T(L,K),L=1,MAXTP),&
                    ((A(M,L,K), M=1,NPCP2), L=1,NTR), K = 1, KK)
!
      IF (II .GT. 0) THEN
!
         WRITE (LINC) (NSPEC(I), NREAC(I), (PAR(N,I), N = 1, NPAR),&
               (NU(M,I), NUNK(M,I), M = 1, MAXSP), I = 1, II)
!
         IF (NREV .GT. 0) WRITE (LINC)&
            (IREV(N),(RPAR(L,N),L=1,NPAR),N=1,NREV)
!
         IF (NFAL .GT. 0) WRITE (LINC)&
            (IFAL(N),IFOP(N),KFAL(N),(PFAL(L,N),L=1,NFAR), N = 1, NFAL)
!
         IF (NTHB .GT. 0) WRITE (LINC)&
            (ITHB(N),NTBS(N),(NKTB(M,N),AIK(M,N),M=1,MAXTB),N=1,NTHB)
!
         IF (NLAN .GT. 0) WRITE (LINC)&
            (ILAN(N), (PLAN(L,N), L = 1, NLAR), N = 1, NLAN)
!
         IF (NRLT .GT. 0) WRITE (LINC)&
            (IRLT(N), (RLAN(L,N), L = 1, NLAR), N=1,NRLT)
!
         IF (NWL .GT. 0) WRITE (LINC) (IWL(N), WL(N), N = 1, NWL)
!
         IF (NEIM .GT. 0) WRITE (LINC) (IEIM(N),ITDEP(N),N=1,NEIM)
!
         IF (NJAN .GT. 0) WRITE (LINC)&
            (IJAN(N), (PJAN(L,N), L = 1, NJAR), N = 1, NJAN)
!
         IF (NFT1 .GT. 0) WRITE (LINC)&
            (IFT1(N), (PFT1(L,N), L = 1, NF1R), N = 1, NFT1)
!
         IF (NEXC .GT. 0) WRITE (LINC)&
            (IEXC(N), PEXC(N), N=1, NEXC)

         IF (NRNU .GT. 0) WRITE (LINC)&
!
!            NRNU, total number of reactions with real stoichiometry
!
            (IRNU(N), (RNU(M,N), M = 1, MAXSP), N = 1, NRNU)
!
!            IRNU, reaction indices
!            RNU,  matrix of real stoichiometric coefficients
!
         IF (NORD .GT. 0) WRITE (LINC)&
!
!            NORD, total number of reactions which use "ORDER"
!
            (IORD(N), (KORD(L,N), RORD(L,N), L=1, MAXORD), N=1,NORD)
!
!            IORD, reaction indices
!            KORD, array of species indices with "ORDER" specified,
!                  -K for forward species, K for reverse species
!            RORD, array of order coefficients
!
      ELSE
         WRITE (LOUT, '(/A)')                       &
            ' WARNING...NO REACTION INPUT FOUND; ', &
            ' LINKING FILE HAS NO REACTION INFORMATION ON IT.'
      ENDIF
!
!     ---------------------------------------------------------------
!     cklink v2 PLOG SECTION (always written; counts are 0 when the
!     mechanism has no PLOG reactions, so the layout is fixed-position
!     and a no-PLOG mechanism's downstream numeric path is unchanged).
!     Records, in order:
!       1) counts        : n_plog_reactions, n_plog_nodes
!       2) reaction map  : plog_reaction(1:nr), plog_node_ptr(0:nr)
!       3) node data     : plog_logP, plog_A, plog_b, plog_EoverR
!                          (each length n_plog_nodes)
!       4) checksum      : integer sum, a cheap end-of-section sentinel
!     Units on disk: logP = ln(P[Pa]); A,b as declared; EoverR = E/R[K].
!     ---------------------------------------------------------------
      WRITE (LINC) plog_n_reactions, plog_n_nodes
      IF (plog_n_reactions .GT. 0) THEN
         WRITE (LINC) (plog_reaction(N), N = 1, plog_n_reactions)
         WRITE (LINC) (plog_node_ptr(N), N = 0, plog_n_reactions)
         WRITE (LINC) (plog_logP_out(N),   N = 1, plog_n_nodes)
         WRITE (LINC) (plog_A_out(N),      N = 1, plog_n_nodes)
         WRITE (LINC) (plog_b_out(N),      N = 1, plog_n_nodes)
         WRITE (LINC) (plog_EoverR_out(N), N = 1, plog_n_nodes)
      ENDIF
!     Section checksum: reactions + nodes + last node-ptr. Not crypto —
!     just catches a truncated/misaligned PLOG section on read.
      WRITE (LINC) plog_n_reactions + plog_n_nodes +                   &
                   MERGE(plog_node_ptr(plog_n_reactions), 0,           &
                         plog_n_reactions .GT. 0)
!
      WRITE (LOUT, '(///A)')&
         ' NO ERRORS FOUND ON INPUT...CHEMKIN LINKING FILE WRITTEN.'
!
      WRITE (LOUT, *) 'KERR = ',KERR

      WRITE (LOUT, '(/A,3(/A,I6))')&
            ' WORKING SPACE REQUIREMENTS ARE',&
            '    INTEGER:   ',LENICK,&
            '    REAL:      ',LENRCK,&
            '    CHARACTER: ',LENCCK
      CLOSE (LINC)
      CLOSE (LOUT)
!
!----------------------------------------------------------------------!
!
!     FORMATS
!
  200 FORMAT (26X,20('-'))
  300 FORMAT (26X,'ELEMENTS',5X,'ATOMIC',/26X,'CONSIDERED',3X,'WEIGHT')
  333 FORMAT (/6X,'Error...no TLO,TMID,THI given for THERMO ALL...'/)
  400 FORMAT (25X,I3,'. ',A4,G15.6)
!
 1000 FORMAT (6X,'Error...no atomic weight for element ',A)
 1070 FORMAT (6X,'ERROR...reaction count exceeds CKINTP capacity IDIM=',&
                   I0,'. Aborting before auxiliary/PLOG parsing.')
 1095 FORMAT (6X,'Error...no duplicate declared for reaction no.',I3)
 1800 FORMAT (///54X, '(k = A T**b exp(-E/RT))',/,&
              6X,'REACTIONS CONSIDERED',30X,'A',8X,'b',8X,'E',/)
!
      RETURN
11111 CONTINUE
      WRITE (LOUT,*) ' Error...cannot read chem.inp...'
      CLOSE (LIN)
      STOP 2
22222 CONTINUE
      WRITE (LOUT,*) ' Error...cannot read therm.dat...'
      CLOSE (LTHRM)
      STOP 2
      END SUBROUTINE CKINTP





!----------------------------------------------------------------------!
      SUBROUTINE CKCHAR (SUB, NSUB, NDIM, STRAY, RAY, NN, KERR, LOUT)
!
!     Extracts names and real values from an array of CHAR*(*)
!     substrings; stores names in STRAY array, real values in RAY;
!     i.e. can be used to store element and atomic weight data,
!     species names, etc.
!
!     Input:   SUB(N),N=1,NSUB  - array of CHAR*(*) substrings
!              NSUB             - number of substrings
!              NDIM             - size of STRAY,RAY arrays
!              NN               - actual number of STRAY found
!              STRAY(N),N=1,NN  - CHAR*(*) array
!              RAY(N),N=1,NN    - Real array
!              LOUT             - output unit for error messages
!     Output:  NN               - incremented if more STRAY found
!              STRAY(N),N=1,NN  - incremented array of STRAY
!              RAY(N),N=1,NN    - incremented array of reals
!              KERR             - logical, .TRUE. = error in data
!
!                                       F. Rupley, Div. 8245, 2/5/88
!----------------------------------------------------------------------!
      USE chemkinII, only: ILASCH, IPPARR, UPCASE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

      DIMENSION RAY(*), PAR(1)
      CHARACTER SUB(*)*(*), STRAY(*)*(*), ISTR*80!, UPCASE*4
      LOGICAL KERR

      ILEN = LEN(STRAY(1))

      DO 200 N = 1, NSUB
         IF ( UPCASE(SUB(N), 3) .EQ. 'END') RETURN
         ISTR = ' '
         I1 = INDEX(SUB(N),'/')
         IF (I1 .EQ. 1) THEN
            KERR = .TRUE.
            WRITE (LOUT, 130) SUB(N)(:ILASCH(SUB(N)))
         ELSE
            IF (I1 .LE. 0) THEN
               ISTR = SUB(N)
            ELSE
               ISTR = SUB(N)(:I1-1)
            ENDIF
            CALL CKCOMP (ISTR, STRAY, NN, INUM)

            IF (INUM .GT. 0) THEN
               WRITE (LOUT, 100) SUB(N)(:ILASCH(SUB(N)))
            ELSE
               IF (NN .LT. NDIM) THEN
                  IF (ISTR(ILEN+1:) .NE. ' ') THEN
                     WRITE (LOUT, 120) SUB(N)(:ILASCH(SUB(N)))
                     KERR = .TRUE.
                  ELSE
                     NN = NN + 1
                     STRAY(NN) = ' '
                     STRAY(NN) = ISTR(:ILEN)
                     IF (I1 .GT. 0) THEN
                        I2 = I1 + INDEX(SUB(N)(I1+1:),'/')
                        ISTR = ' '
                        ISTR = SUB(N)(I1+1:I2-1)
                        CALL IPPARR (ISTR, 1, 1, PAR, NVAL, IER, LOUT)
                        KERR = KERR .OR. (IER.NE.0)
                        RAY(NN) = PAR(1)
                     ENDIF
                  ENDIF
               ELSE
                  WRITE (LOUT, 110) SUB(N)(:ILASCH(SUB(N)))
                  KERR = .TRUE.
               ENDIF
            ENDIF
         ENDIF
  200 CONTINUE

  100 FORMAT (6X,'Warning...duplicate array element ignored...',A)
  110 FORMAT (6X,'Error...character array size too small for  ...',A)
  120 FORMAT (6X,'Error...character array element name too long...',A)
  130 FORMAT (6X,'Error...misplaced value...',A)
      END SUBROUTINE CKCHAR
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
!
      SUBROUTINE CKCOMP (IST, IRAY, II, I)
!
!  START PROLOGUE
!
!  SUBROUTINE CKCOMP (IST, IRAY, II, I)*
!     Returns the index of an element of a reference character
!     string array which corresponds to a character string;
!     leading and trailing blanks are ignored.
!
!
!  INPUT
!     IST   - A character string.
!                  Data type - CHARACTER*(*)
!     IRAY  - An array of character strings;
!             dimension IRAY(*) at least II
!                  Data type - CHARACTER*(*)
!     II    - The length of IRAY.
!                  Data type - integer scalar.
!
!  OUTPUT
!     I     - The first integer location in IRAY in which IST
!             corresponds to IRAY(I); if IST is not also an
!             entry in IRAY, I=0.
!
!  END PROLOGUE
!
      USE chemkinII, only: IFIRCH, ILASCH
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

!
      CHARACTER*(*) IST, IRAY(*)
!
      I = 0
      DO 10 N = II, 1, -1
         IS1 = IFIRCH(IST)
         IS2 = ILASCH(IST)
         IR1 = IFIRCH(IRAY(N))
         IR2 = ILASCH(IRAY(N))
         IF ( IS2.GE.IS1 .AND. IS2.GT.0 .AND. IS1>0 .AND. &
              IR2.GE.IR1 .AND. IR2.GT.0 .AND. IR1>0) THEN
              IF ( &
              IST(IS1:IS2).EQ.IRAY(N)(IR1:IR2) ) I=N
         ENDIF
   10 CONTINUE
      RETURN
      END SUBROUTINE CKCOMP
!



!----------------------------------------------------------------------!
      SUBROUTINE CKAWTM (ENAME, AWT)
!
!     Returns atomic weight of element ENAME.
!     Input:   ENAME - CHAR*(*) element name
!     Output:  AWT   - real atomic weight
!
!                                       F. Rupley, Div. 8245, 11/11/86
!----------------------------------------------------------------------!
      use chemkinII, only: UPCASE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

      PARAMETER (NATOM = 102)
      DIMENSION ATOM(NATOM)
      CHARACTER ENAME*(*), IATOM(NATOM)*2!, UPCASE*2
!
      DATA (IATOM(I),ATOM(I),I=1,40) /                               &
      'H ',  1.00797, 'HE',  4.00260, 'LI',  6.93900, 'BE',  9.01220,&
      'B ', 10.81100, 'C ', 12.01115, 'N ', 14.00670, 'O ', 15.99940,&
      'F ', 18.99840, 'NE', 20.18300, 'NA', 22.98980, 'MG', 24.31200,&
      'AL', 26.98150, 'SI', 28.08600, 'P ', 30.97380, 'S ', 32.06400,&
      'CL', 35.45300, 'AR', 39.94800, 'K ', 39.10200, 'CA', 40.08000,&
      'SC', 44.95600, 'TI', 47.90000, 'V ', 50.94200, 'CR', 51.99600,&
      'MN', 54.93800, 'FE', 55.84700, 'CO', 58.93320, 'NI', 58.71000,&
      'CU', 63.54000, 'ZN', 65.37000, 'GA', 69.72000, 'GE', 72.59000,&
      'AS', 74.92160, 'SE', 78.96000, 'BR', 79.90090, 'KR', 83.80000,&
      'RB', 85.47000, 'SR', 87.62000, 'Y ', 88.90500, 'ZR', 91.22000/
!
      DATA (IATOM(I),ATOM(I),I=41,80) /                              &
      'NB', 92.90600, 'MO', 95.94000, 'TC', 99.00000, 'RU',101.07000,&
      'RH',102.90500, 'PD',106.40000, 'AG',107.87000, 'CD',112.40000,&
      'IN',114.82000, 'SN',118.69000, 'SB',121.75000, 'TE',127.60000,&
      'I ',126.90440, 'XE',131.30000, 'CS',132.90500, 'BA',137.34000,&
      'LA',138.91000, 'CE',140.12000, 'PR',140.90700, 'ND',144.24000,&
      'PM',145.00000, 'SM',150.35000, 'EU',151.96000, 'GD',157.25000,&
      'TB',158.92400, 'DY',162.50000, 'HO',164.93000, 'ER',167.26000,&
      'TM',168.93400, 'YB',173.04000, 'LU',174.99700, 'HF',178.49000,&
      'TA',180.94800, 'W ',183.85000, 'RE',186.20000, 'OS',190.20000,&
      'IR',192.20000, 'PT',195.09000, 'AU',196.96700, 'HG',200.59000/
!
      DATA (IATOM(I),ATOM(I),I=81,NATOM) /                           &
      'TL',204.37000, 'PB',207.19000, 'BI',208.98000, 'PO',210.00000,&
      'AT',210.00000, 'RN',222.00000, 'FR',223.00000, 'RA',226.00000,&
      'AC',227.00000, 'TH',232.03800, 'PA',231.00000, 'U ',238.03000,&
      'NP',237.00000, 'PU',242.00000, 'AM',243.00000, 'CM',247.00000,&
      'BK',249.00000, 'CF',251.00000, 'ES',254.00000, 'FM',253.00000,&
      'D ',002.01410, 'E',5.45E-4/
!
      CALL CKCOMP ( UPCASE(ENAME, 2), IATOM, NATOM, L)
      IF (L .GT. 0) AWT = ATOM(L)
      RETURN
      END SUBROUTINE CKAWTM
!----------------------------------------------------------------------!


!----------------------------------------------------------------------!
      SUBROUTINE CKTHRM (LUNIT, MDIM, ENAME, MM, AWT, KNAME, KK, KNCF, &
                         KPHSE, KCHRG, WTM, MAXTP, NT, NTR, TLO, TMID, &
                         THI, T, NPCP2, A, ITHRM, KERR, LOUT, ISTR)
!
!     Finds thermodynamic data and elemental composition for species
!     Input:  LUNIT  - unit number for input of thermo properties
!             MDIM   - maximum number of elements allowed
!             ENAME(M),M=1,MM  - array of CHAR*(*) element names
!             MM     - total number of elements declared
!             AWT(M),M=1,MM    - array of atomic weights for elements
!             KNAME(K),K=1,KK  - array of CHAR*(*) species names
!             KK     - total number of species declared
!             LOUT   - output unit for messages
!             NT(K),K=1,KK - number of temperature values
!             NTR - number of temperature ranges
!     Output: KNCF(M,K) - elemental composition of species
!             KPHSE(K),K=1,KK - integer array, species phase
!             KCHRG(K),K=1,KK - integer array of species charge
!                      =0, if no electrons,
!                      =(-1)*number of electrons present
!             WTM(K),K=1,KK - array of molecular weights of species
!             A(M,L,K)- array of thermodynamic coefficients
!             T(N),N=1,NT - array of temperatures
!             KERR   - logical error flag
!----------------------------------------------------------------------!
      USE chemkinII, only: IPPLEN, IPPARR, UPCASE
      USE chemistry_string_limits, only: species_name_len
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
!
      DIMENSION WTM(*), NT(*), T(MAXTP,*), KPHSE(*), KNCF(MDIM,*), &
                KCHRG(*), A(NPCP2,NTR,*), AWT(*), VALUE(5)
      CHARACTER ENAME(*)*(*), KNAME(*)*(*), LINE(4)*80, ELEM*16
      CHARACTER ISTR*80, SUB(80)*80! UPCASE*4,
      LOGICAL KERR, ITHRM(*)
      INTEGER :: K
!
      IF (MM.LE.0 .OR. KK.LE.0) WRITE (LOUT, 80)
!
      GO TO 20
   10 CONTINUE
      ISTR = ' '
      READ (LUNIT,'(A)',END=40) ISTR
   20 CONTINUE
      ILEN = IPPLEN(ISTR)
      IF (ILEN .LE. 0) GO TO 10
!
      CALL CKISUB (ISTR(:ILEN), SUB, NSUB)
      IF (UPCASE(SUB(1), 3) .EQ. 'END' .OR. &
          UPCASE(SUB(1), 4) .EQ. 'REAC') RETURN
!
      IF (ILEN.LT.80 .OR. ISTR(80:80).NE.'1') GO TO 10
!     NASA-7 reserves columns 1:18 for the species identifier. Reading
!     that fixed field supports a full 18-character name even when the
!     date field starts immediately in column 19 with no separating blank.
      CALL CKCOMP (ISTR(:species_name_len), KNAME, KK, K)
!
      IF (K.LE.0)  GO TO 10
      IF (ITHRM(K)) GO TO 10
      ITHRM(K) = .TRUE.
      LINE(1) = ' '
      LINE(1) = ISTR
      L = 2
  111 CONTINUE
      READ (LUNIT,'(A)',END=40) LINE(L)
      IF (IPPLEN(LINE(L)) .GE. 80) THEN
         IF (LINE(L)(80:80) .EQ. '4') THEN
            GO TO 25
         ELSEIF (LINE(L)(80:80).EQ.'2' .OR. &
                 LINE(L)(80:80).EQ.'3') THEN
            L = L + 1
         ENDIF
      ENDIF
      GO TO 111
!
   25 CONTINUE
!



      ICOL = 20
      DO 60 I = 1, 5
         ICOL = ICOL + 5
         IF (I .EQ. 5) ICOL = 74
         ELEM  = LINE(1)(ICOL:ICOL+1)
         IELEM = 0
!
         IF (LINE(1)(ICOL+2:ICOL+4) .NE. ' ') THEN
            CALL IPPARR &
            (LINE(1)(ICOL+2:ICOL+4), 0, 1, VALUE, NVAL, IER, LOUT)
            IELEM = VALUE(1)
         ENDIF
!
         IF (ELEM.NE.' ' .AND. IELEM.NE.0) THEN
            IF (UPCASE(ELEM, 1) .EQ. 'E') &
                   KCHRG(K)=KCHRG(K)+IELEM*(-1)
            CALL CKCOMP (ELEM, ENAME, MM, M)
            IF (M .GT. 0) THEN
               KNCF(M,K) = IELEM
               WTM(K) = WTM(K) + AWT(M)*FLOAT(IELEM)
            ELSE
               WRITE (LOUT, 100) ELEM,KNAME(K)(:10)
               KERR = .TRUE.
            ENDIF
         ENDIF
   60 CONTINUE
!
      IF (UPCASE(LINE(1)(45:),1) .EQ. 'L') KPHSE(K)=1
      IF (UPCASE(LINE(1)(45:),1) .EQ. 'S') KPHSE(K)=-1
!
!-----Currently allows for three temperatures, two ranges;
!     in future, NT(K) may vary, NTR = NT(K)-1
!
      T(1,K) = TLO
      IF (LINE(1)(46:55) .NE. ' ') CALL IPPARR &
         (LINE(1)(46:55), 0, 1, T(1,K), NVAL, IER, LOUT)
!
      T(2,K) = TMID
      IF (LINE(1)(66:73) .NE. ' ') CALL IPPARR &
         (LINE(1)(66:73), 0, 1, T(2,K), NVAL, IER, LOUT)
!
      T(NT(K),K) = THI
      IF (LINE(1)(56:65) .NE. ' ') CALL IPPARR &
         (LINE(1)(56:65), 0, 1, T(NT(K),K), NVAL, IER, LOUT)
!
      READ (LINE(2)(:75),'(5E15.8)') (A(I,NTR,K),I=1,5)
      READ (LINE(3)(:75),'(5E15.8)') (A(I,NTR,K),I=6,7),(A(I,1,K),I=1,3)
      READ (LINE(4)(:60),'(4E15.8)') (A(I,1,K),I=4,7)
      GO TO 10
!
   40 RETURN
   80 FORMAT (6X,'Warning...THERMO cards misplaced will be ignored...')
  100 FORMAT (6X,'Error...element...',A,'not declared for...',A)
      END SUBROUTINE CKTHRM
!----------------------------------------------------------------------!


!----------------------------------------------------------------------!
      SUBROUTINE CKREAC (LINE, II, KK, KNAME, LOUT, MAXSP, NSPEC, NREAC,&
                         NUNK, NU, NPAR, PAR, NTHB, ITHB,               &
                         NFAL, IFAL, KFAL, NWL, IWL, WL,                &
                         NRNU, IRNU, RNU, KERR)
!
!     CKREAC parses the main CHAR*(*) line representing a gas-phase
!     reaction; first, the real Arrhenius parameters are located and
!     stored in PAR(N,I),N=1,NPAR, where I is the reaction number;
!     then a search is made over the reaction string:
!
!     '=','<=>': reaction I is reversible;
!     '=>'     : reaction I is irreversible;
!
!     '(+[n]KNAME(K))': reaction I is a fall-off reaction;
!                       NFAL is incremented, the total number of
!                       fall-off reactions;
!                       IFAL(NFAL)=I, KFAL(NFAL)=K;
!                       this species is eliminated from consideration
!                       as a reactant or product in this reaction.
!
!     '(+M)'   : reaction I is a fall-off reaction;
!                NFAL is incremented, IFAL(NFAL)=I, KFAL(NFAL)=0;
!
!     '+[n]KNAME(K)': NSPEC(I) is incremented, the total number of
!                     species for this reaction;
!                     n is an optional stoichiometric coefficient
!                     of KNAME(K), if omitted, n=1;
!                     if this string occurs before the =/-,
!                     NREAC(I) is incremented, the total number of
!                     reactants for this reaction, NUNK(N,I)=K, and
!                     NU(N,I) = -n, where N=1-3 is reserved for
!                     reactants;
!                     if this string occurs after the =/-,
!                     NUNK(N,I) = K, and NU(N,I) = n, where N=4-6
!                     is reserved for products;
!
!     '+M' : I is a third-body reaction; NTHB is incremented, the
!            total number of third-body reactions, and ITHB(NTHB)=I.
!
!     Input:  LINE  - a CHAR*(*) line (from data file)
!             II    - the index of this reaction, and the total number
!                     of reactions found so far.
!             KK    - actual integer number of species
!             KNAME(K),K=1,KK - array of CHAR*(*) species names
!             LOUT  - output unit for error messages
!             MAXSP - maximum number of species allowed in reaction
!             NPAR  - number of parameters expected
!     A '!' will comment out a line, or remainder of the line.
!
!     Output: NSPEC - total number of reactants+products in reaction
!             NREAC - number of reactants
!             NUNK  - the NSPEC species indices
!             NU    - the NSPEC stoichiometric coefficients
!             NFAL  - total number of fall-off reactions
!             IFAL  - the NFAL reaction indices
!             KFAL  - 3rd body species indices for the NFAL reactions
!             NTHB  - total number of 3rd-body reactions
!             ITHB  - the NTHB reaction indices
!             NWL   - number of radiation-enhanced reactions
!             IWL   - the NWL reaction indices
!             WL    - the NWL radiation wavelengths
!             KERR  - logical, .TRUE. = error in data file
!
!                                      F. Rupley, Div. 8245, 5/13/86
!----------------------------------------------------------------------!
      USE chemkinII, only: ILASCH, IPPARR, IPPARI, UPCASE
      USE chemistry_string_limits, only: mechanism_line_len
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

      DIMENSION NSPEC(*), NREAC(*), NUNK(MAXSP,*), NU(MAXSP,*),       &
                PAR(NPAR,*), IFAL(*), KFAL(*), ITHB(*), IWL(*), WL(*),&
                IRNU(*), RNU(MAXSP,*), IPLUS(20), RVAL(10), IVAL(10)
      CHARACTER KNAME(*)*(*), LINE*(*), CNUM(11)*1!, UPCASE*4
      CHARACTER(len=mechanism_line_len) :: ISTR, IREAC, IPROD, ISPEC, &
                                            INAME, ITEMP
      LOGICAL KERR, LTHB, LWL, LRSTO
      DATA CNUM/'.','0','1','2','3','4','5','6','7','8','9'/
!
      LTHB = .FALSE.
      LWL = .FALSE.
      NSPEC(II) = 0
      NREAC(II) = 0
!
!----------Find NPAR real parameters------------------------
!
      CALL IPNPAR (LINE, NPAR, ISTR, ISTART)
      CALL IPPARR (ISTR, 1, NPAR, PAR(1,II), NVAL, IER, LOUT)
      KERR = KERR .OR. (IER.NE.0)
!
!-----Remove blanks from reaction string
!
      INAME = ' '
      ILEN = 0
      DO 10 I = 1, ISTART-1
         IF (LINE(I:I) .NE. ' ') THEN
            ILEN = ILEN+1
            INAME(ILEN:ILEN) = LINE(I:I)
         ENDIF
   10 CONTINUE
!
!-----Find reaction string, product string
!
      I1 = 0
      I2 = 0
      DO 25 I = 1, ILEN
         IF (I1 .LE. 0) THEN
            IF (INAME(I:I+2) .EQ. '<=>') THEN
               I1 = I
               I2 = I+2
               IR = 1
            ELSEIF (INAME(I:I+1) .EQ. '=>') THEN
               I1 = I
               I2 = I+1
               IR = -1
            ELSEIF (I.GT.1) THEN
              IF (INAME(I:I    ).EQ.'='  .AND. &
                  INAME(I-1:I-1).NE.'=') THEN
               I1 = I
               I2 = I
               IR = 1
              ENDIF
            ENDIF
         ENDIF
   25 CONTINUE
!
      IF (ILASCH(INAME).GE.45 .AND. I1.GT.0) THEN
         WRITE (LOUT, 1900) II,INAME(:I1-1),(PAR(N,II),N=1,NPAR)
         WRITE (LOUT, 1920) INAME(I1:)
      ELSE
          WRITE (LOUT, 1900) II,INAME(:45),(PAR(N,II),N=1,NPAR)
      ENDIF
!
      IREAC = ' '
      IPROD = ' '
      IF (I1 .GT. 0) THEN
         IREAC = INAME(:I1-1)
         IPROD = INAME(I2+1:)
      ELSE
!
!-----did not find delimiter
!
         WRITE (LOUT, 660)
         KERR = .TRUE.
         RETURN
      ENDIF
!
      LRSTO = ((INDEX(IREAC,'.').GT.0) .OR. (INDEX(IPROD,'.').GT.0))
      IF (LRSTO) THEN
         NRNU = NRNU + 1
         IRNU(NRNU) = II
      ENDIF
!
      IF (INDEX(IREAC,'=>').GT.0 .OR. INDEX(IPROD,'=>').GT.0) THEN
!
!-----more than one '=>'
!
         WRITE (LOUT, 800)
         KERR = .TRUE.
         RETURN
      ENDIF
!
!-----Is this a fall-off reaction?
!
      IF (INDEX(IREAC,'(+').GT.0 .OR. INDEX(IPROD,'(+').GT.0) THEN
         KRTB = 0
         KPTB = 0
         DO 300 J = 1, 2
            ISTR = ' '
            KTB  = 0
            IF (J .EQ. 1) THEN
               ISTR = IREAC
            ELSE
               ISTR = IPROD
            ENDIF
!
            DO 35 N = 1, ILASCH(ISTR)-1
               IF (ISTR(N:N+1) .EQ. '(+') THEN
                  I1 = N+2
                  I2 = I1 + INDEX(ISTR(I1:),')')-1
                  IF (I2 .GT. I1) THEN
                     IF (ISTR(I1:I2-1).EQ.'M' .OR.  &
                         ISTR(I1:I2-1).EQ.'m') THEN
                         IF (KTB .NE. 0) THEN
                            WRITE (LOUT, 630)
                            KERR = .TRUE.
                            RETURN
                         ELSE
                            KTB = -1
                         ENDIF
                     ELSE
                        CALL CKCOMP (ISTR(I1:I2-1), KNAME, KK, KNUM)
                        IF (KNUM .GT. 0) THEN
                           IF (KTB .NE. 0) THEN
                              WRITE (LOUT, 630)
                              KERR = .TRUE.
                              RETURN
                           ELSE
                              KTB = KNUM
                           ENDIF
                        ENDIF
                     ENDIF
                     IF (KTB .NE. 0) THEN
                        ITEMP = ' '
                        IF (I1 .EQ. 1) THEN
                           ITEMP = ISTR(I2+1:)
                        ELSE
                           ITEMP = ISTR(:I1-3)//ISTR(I2+1:)
                        ENDIF
                        IF (J .EQ. 1) THEN
                           IREAC = ' '
                           IREAC = ITEMP
                           KRTB = KTB
                        ELSE
                           IPROD = ' '
                           IPROD = ITEMP
                           KPTB = KTB
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
   35       CONTINUE
  300    CONTINUE
!
         IF (KRTB.NE.0 .OR. KPTB.NE.0) THEN
!
!           does product third-body match reactant third-body
!
            IF (KRTB.LE.0 .AND. KPTB.LE.0) THEN
!
               NFAL = NFAL + 1
               IFAL(NFAL) = II
               KFAL(NFAL) = 0
!
               LTHB = .TRUE.
               NTHB = NTHB + 1
               ITHB(NTHB) = II
!
            ELSEIF (KRTB .EQ. KPTB) THEN
               NFAL = NFAL + 1
               IFAL(NFAL) = II
               KFAL(NFAL) = KRTB
!
            ELSE
!
               WRITE (LOUT, 640)
               KERR = .TRUE.
               RETURN
            ENDIF
         ENDIF
      ENDIF
!
!----------Find reactants, products-------------------------
!
      DO 600 J = 1, 2
         ISTR = ' '
         LTHB = .FALSE.
         IF (J .EQ. 1) THEN
            ISTR = IREAC
            NS = 0
         ELSE
            ISTR = IPROD
            NS = 3
         ENDIF
!
!-----------store pointers to '+'-signs
!
         NPLUS = 1
         IPLUS(NPLUS) = 0
         DO 500 L = 2, ILASCH(ISTR)-1
            IF (ISTR(L:L).EQ.'+') THEN
               NPLUS = NPLUS + 1
               IPLUS(NPLUS) = L
            ENDIF
  500    CONTINUE
         NPLUS = NPLUS + 1
         IPLUS(NPLUS) = ILASCH(ISTR)+1
!
         NSTART = 1
  505    CONTINUE
         N1 = NSTART
         DO 510 N = NPLUS, N1, -1
            ISPEC = ' '
            ISPEC = ISTR(IPLUS(N1)+1 : IPLUS(N)-1)
!
            IF (UPCASE(ISPEC, 1).EQ.'M' .AND. &
                     (ISPEC(2:2).EQ.' ' .OR. ISPEC(2:2).EQ.'+')) THEN
               IF (LTHB) THEN
                  WRITE (LOUT, 900)
                  KERR = .TRUE.
                  RETURN
               ELSEIF (NFAL.GT.0 .AND. IFAL(MAX(1,NFAL)).EQ.II) THEN
                  WRITE (LOUT, 640)
                  KERR = .TRUE.
                  RETURN
               ELSE
                  LTHB = .TRUE.
                  IF (NTHB.EQ.0 .OR. &
                     (NTHB.GT.0.AND.ITHB(MAX(1,NTHB)).NE.II)) THEN
                      NTHB = NTHB + 1
                      ITHB(NTHB) = II
                  ENDIF
                  IF (N .EQ. NPLUS) GO TO 600
                  NSTART = N
                  GO TO 505
               ENDIF
!
            ELSEIF (UPCASE(ISPEC, 2) .EQ. 'HV') THEN
               IF (LWL) THEN
                  WRITE (LOUT, 670)
                  KERR = .TRUE.
                  RETURN
               ELSE
                  LWL = .TRUE.
                  NWL = NWL + 1
                  IWL(NWL) = II
                  WL(NWL) = 1.0
                  IF (J .EQ. 1) WL(NWL) = -1.0
                  IF (N .EQ. NPLUS) GO TO 600
                  NSTART = N
                  GO TO 505
               ENDIF
            ENDIF
!
!-----------does this string start with a number?
!
            IND = 0
            DO 334 L = 1, LEN(ISPEC)
               NTEST = 0
               DO 333 M = 1, 11
                  IF (ISPEC(L:L) .EQ. CNUM(M)) THEN
                     NTEST=M
                     IND = L
                  ENDIF
  333          CONTINUE
               IF (NTEST .EQ. 0) GO TO 335
  334       CONTINUE
  335       CONTINUE
!
            RVAL = 1.0
            IVAL = 1
            IF (IND .GT. 0) THEN
               IF (LRSTO) THEN
                  CALL IPPARR (ISPEC(:IND), 1, 1, RVAL(1), NVAL, &
                               IER, LOUT)
               ELSE
                  CALL IPPARI (ISPEC(:IND), 1, 1, IVAL(1), NVAL, &
                              IER, LOUT)
               ENDIF
               IF (IER .EQ. 0) THEN
                  ITEMP = ' '
                  ITEMP = ISPEC(IND+1:)
                  ISPEC = ' '
                  ISPEC = ITEMP
               ELSE
                  KERR = .TRUE.
                  RETURN
               ENDIF
            ENDIF
!
            CALL CKCOMP (ISPEC, KNAME, KK, KNUM)
            IF (KNUM .EQ. 0) THEN
               IF ((N-N1) .GT. 1) GO TO 510
               WRITE (LOUT, 680) ISPEC(:ILASCH(ISPEC))
               KERR = .TRUE.
            ELSE
!
!--------------a species has been found
!
               IF (J .EQ. 1) THEN
                  IVAL(1) = -IVAL(1)
                  RVAL(1) = -RVAL(1)
               ENDIF
!
!--------------increment species coefficient count
!
               NNUM = 0
               IF (LRSTO) THEN
                  DO 110 K = 1, NS
                     IF (KNUM.EQ.NUNK(K,II) .AND. &
                         RNU(K,NRNU)/RVAL(1).GT.0) THEN
                         NNUM = K
                         RNU(NNUM,NRNU) = RNU(NNUM,NRNU) + RVAL(1)
                     ENDIF
  110             CONTINUE
               ELSE
                  DO 111 K = 1, NS
                     IF (KNUM.EQ.NUNK(K,II) .AND. &
                         NU(K,II)/IVAL(1).GT.0) THEN
                        NNUM=K
                        NU(NNUM,II) = NU(NNUM,II) + IVAL(1)
                     ENDIF
  111             CONTINUE
               ENDIF
!
               IF (NNUM .LE. 0) THEN
!
!-----------------are there too many species?
!
                  IF (J.EQ.1 .AND. NS.EQ.3) THEN
                     WRITE (LOUT, 690)
                     KERR = .TRUE.
                     RETURN
                  ELSEIF (J.EQ.2 .AND. NS.EQ.MAXSP) THEN
                     WRITE (LOUT, 700)
                     KERR = .TRUE.
                     RETURN
                  ELSE
!
!--------------------increment species count
!
                     NS = NS + 1
                     NSPEC(II) = NSPEC(II)+1
                     IF (J .EQ. 1) NREAC(II) = NS
                     NUNK(NS,II) = KNUM
                     IF (LRSTO) THEN
                        RNU(NS,NRNU) = RVAL(1)
                     ELSE
                        NU(NS,II)   = IVAL(1)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            IF (N .EQ. NPLUS) GO TO 600
            NSTART = N
            GO TO 505
!
  510    CONTINUE
  600 CONTINUE
!
      NSPEC(II) = IR*NSPEC(II)
!
  630 FORMAT (6X,'Error...more than one fall-off declaration...')
  640 FORMAT (6X,'Error in fall-off declaration...')
  650 FORMAT (6X,'Error...reaction string not found...')
  660 FORMAT (6X,'Error in reaction...')
  670 FORMAT (6X,'Error in HV declaration...')
  680 FORMAT (6X,'Error...undeclared species...',A)
  690 FORMAT (6X,'Error...more than 3 reactants...')
  700 FORMAT (6X,'Error...more than 3 products...')
  800 FORMAT (6X,'Error in reaction delimiter...')
  900 FORMAT (6X,'Error in third-body declaration...')
! 1900 FORMAT (I4,'. ',A,T51,E10.3,F7.3,F11.3)
 1900 FORMAT (I4,'. ', A, T53, 1PE8.2, 2X, 0PF5.1, 2X, F9.1)
 1920 FORMAT (6X,A)
      RETURN
      END SUBROUTINE CKREAC

!-----------------------------------------------------------------------

!----------------------------------------------------------------------!
      SUBROUTINE CKAUXL (SUB, NSUB, II, KK, KNAME, LOUT, MAXSP, NPAR,   &
                         NSPEC, NTHB, ITHB, NTBS, MAXTB, NKTB, AIK,     &
                         NFAL, IFAL, IDUP, NFAR, PFAL, IFOP, NLAN,      &
                         ILAN, NLAR, PLAN, NREV, IREV, RPAR, NRLT, IRLT,&
                         RLAN, NWL, IWL, WL, KERR, NORD, IORD, MAXORD,  &
                         KORD, RORD, NUNK, NU, NRNU, IRNU, RNU,         &
                         NEIM, IEIM, ITDEP, NJAN, IJAN, NJAR, PJAN,     &
                         NFT1, IFT1, NF1R, PFT1, NEXC, IEXC, PEXC)
!
!     CKAUXL parses the auxiliary CHAR*(*) lines representing
!     additional options for a gas-phase reaction; data is stored
!     based on finding a 'keyword' followed by its required
!     parameters:
!
!     KNAME(K)/val1/: this is an enhanced third-body;
!
!        if ITHB(NTHB) <> I, this is an error, reaction I is not a
!                            third-body reaction;
!        else NTBS(NTHB) is incremented,
!             AIK(NTBS(NTHB),NTHB) = K,
!             NKTB(NTBS(NTHB)),NTHB) = val1;
!
!     (LOW,TROE, and SRI define fall-off data):
!
!     LOW/val1 val2 val3/: PFAL(N,NFAL) = val(N),N=1,3;
!
!        if IFAL(NFAL)<>I, this is an error, reaction I is not a
!                          fall-off reaction;
!        if ILAN(NLAN)=I, this is an error, cannot have T-L numbers.
!        if IRLT(NRLT)=I, this is an error,         "
!        if IREV(NREV)=I, this is an error, cannot declare reverse
!                         parameters;
!        if IFOP(NFAL)>0, this is an error, LOW already declared;
!        else
!           IFOP(NFAL) = ABS(IFOP(NFAL))
!
!     TROE/val1 val2 val3 [val4]/: PFAL(N,NFAL) = val(N),N=4,7;
!
!        if IFAL(NFAL)<>I, this is an error, reaction I is not a
!                          fall-off reaction;
!        if ILAN(NLAN)=I, this is an error, cannot have T-L numbers.
!        if IRLT(NRLT)=I, this is an error,         "
!        if IREV(NREV)=I, this is an error, cannot declare reverse
!                         parameters;
!        if ABS(IFOP(NFAL)).GT.1, this is an error,
!        else
!        if 3 TROE values, IFOP(NFAL) = 3*IFOP(NFAL);
!        if 4 TROE values, IFOP(NFAL) = 4*IFOP(NFAL);
!
!     SRI/val1 val2 val3/: PFAL(N,NFAL) = val(N),N=4,6;
!
!        if IFAL(NFAL)<>I, this is an error, reaction I is not a
!                          fall-off reaction;
!        if ILAN(NLAN)=I, this is an error, cannot have T-L numbers.
!        if IRLT(NRLT)=I, this is an error,         "
!        if IREV(NREV)=I, this is an error, cannot declare reverse
!                         parameters;
!        if ABS(IFOP(NFAL))>1, this is an error;
!        else
!        if IFOP(NFAL)= 2*IFOP(NFAL);
!
!     LT/val1 val2/:
!        if IFAL(NFAL)=I, this is an error, cannot have fall-off and
!                         T-L numbers;
!        else increment NLAN, the number of T-L reactions,
!             ILAN(NLAN)=I, PLAN(N,NLAN)=val(N),N=1,2
!        if IREV(NREV)=I, need IRLT(NRLT)=I.
!
!     REV[ERSE]/val1 val2 val3/ :
!        if IFAL(NFAL)=I, this is an error;
!        if IREV(NREV)=I, this is an error, REV already declared;
!        if NSPEC(I)<0, this an error, as I is irreversible;
!        else increment NREV, the number of reactions with reverse
!             parameters given,
!             IREV(NREV)=I, RPAR(N,NREV)=val(N),N=1,3;
!             if ILAN(NLAN)=I, need IRLT(NRLT)=I;
!             if IRLT(NRLT)=I, need ILAN(NRLT)=I.
!
!     RLT/val1 val2/:
!       if IFAL(NFAL)=I, this is an error, cannot have fall-off and
!                        T-L numbers;
!       if IRLT(NRLT)=I, this is an error, RLT already declared;
!       else increment NRLT, the number of reactions with BOTH
!                      reverse parameters given, and T-L numbers;
!            IRLT(NRLT)=I, RLAN(N,NRLT)=val(N),N=1,2;
!            if IREV(NREV)<>I, need IREV(NREV)=I;
!            if ILAN(NREV)<>I, need ILAN(NLAN)=I;
!
!    DUP[LICATE]:
!       This reaction is allowed to be duplicated.
!
!    EIM/VAL1/:
!       if ITHB(NTHB)=I, this is an error, cannot have both
!       neutral 3rd-body dependence and e- impact
!
!     Input:  LINE - CHAR*(*) auxiliary information string
!             KK   - total number of species declared
!             KNAME- CHAR*(*) species names
!             LOUT - output unit for error messages
!             MAXSP- maximum third bodies allowed in a reaction
!     Output: NTHB - total number of reactions with third bodies
!             ITHB - the NTHB reaction indices
!             AIK  - non-zero third body enhancement factors
!             NKTB - array of species indices for the third body
!                         enchancement factors
!             NFAL - total number of fall-off reactions
!             IFAL - the NFAL reaction indices
!             IFOP - the NFAL fall-off types
!             PFAL - fall-off parameters
!             NLAN - total number of Landau-Teller reactions
!             ILAN - the NLAN reaction indices
!             NLAR - number of Landau-Teller numbers allowed
!             PLAN - array of Landau-Teller numbers
!             NRLT - total number of 'reverse' T-L reactions
!             IRLT - the NRLT reaction indices
!             RLAN - array of 'reverse' Landau-Teller numbers
!             NWL  - total number of radiation-enhanced reactions
!             IWL  - the NWL reaction indices
!             WL   - the NWL wavelengths
!             NEIM - total number of electron-impact reactions
!             IEIM - the NEIM reaction indices
!             ITDEP- the NEIM temperature dependence flags
!             NJAN - total number of Jannev, Langer, Evans & Post types
!             IJAN - the NJAN reaction indices
!             NJAR - number of coefficients required for J,L,E&P reacts
!             PJAN - array of coefficients of J,L,E&P reactions
!             NFT1 - number of fit#1 reactions
!             IFT1 - the NFT1 reacton indices
!             PFT1 - array of added exponential parameters for fit#1
!             NEXC - number of excitation reactions
!             IEXC - the NEXC reaction indices
!             PEXC - the NEXC energy loss (eV)
!             NRNU - number of real stoichiometry reactions
!             IRNU - the NRNU reaction indices
!             RNU  - matrix of coefficients for the NRNU reactions
!             NORD - number of change-of-order reactions
!             IORD - the IORD reaction indices
!             KORD - matrix of species indices for the NORD reactions
!             RORD - matrix of order values for the NORD reactions
!             KERR - logical, = .TRUE. if error found
!                                        F. Rupley, Div. 8245, 5/27/87
!----------------------------------------------------------------------!
      USE chemkinII, only: IPPARR, IPPARI, ILASCH, UPCASE
      USE chemistry_string_limits, only: mechanism_line_len
      USE plog_collect, only: plog_add_line   ! PLOG collection (cklink v2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
!
      DIMENSION NSPEC(*), ITHB(*), NTBS(*), NKTB(MAXTB,*), IDUP(*),   &
                AIK(MAXTB,*), IFAL(*), IFOP(*), PFAL(NFAR,*),         &
                ILAN(*), PLAN(NLAR,*), IREV(*), RPAR(NPAR,*), IRLT(*),&
                RLAN(NLAR,*), IWL(*), WL(*), VAL(1), IORD(*),         &
                KORD(MAXORD,*), RORD(MAXORD,*), NUNK(MAXSP,*),        &
                NU(MAXSP,*), IRNU(*), RNU(MAXSP,*)
!
      DIMENSION IEIM(*), ITDEP(*), IJAN(*), PJAN(NJAR,*), IFT1(*),    &
                PFT1(NF1R,*), IEXC(*), PEXC(*)
!
      DIMENSION PLOGV(4)   ! scratch for one `PLOG / P A b E /` line
!
      CHARACTER SUB(*)*(*), KNAME(*)*(*), KEY*80,                     & ! UPCASE*4,  &
                RSTR*mechanism_line_len, ISTR*mechanism_line_len
      LOGICAL KERR, LLAN, LRLT, LTHB, LFAL, LTRO, LSRI, LWL, LREV,    &
              LFORD, LRORD, LEIM, LJAN, LFT1, LEXC
!
      LTHB = (NTHB.GT.0 .AND. ITHB(MAX(1,NTHB)).EQ.II)
      LFAL = (NFAL.GT.0 .AND. IFAL(MAX(1,NFAL)).EQ.II)
      LWL  = (NWL .GT.0 .AND. IWL (MAX(1,NWL) ).EQ.II)
      LREV = (NREV.GT.0 .AND. IREV(MAX(1,NREV)).EQ.II)
      LLAN = (NLAN.GT.0 .AND. ILAN(MAX(1,NLAN)).EQ.II)
      LRLT = (NRLT.GT.0 .AND. IRLT(MAX(1,NRLT)).EQ.II)
      LTRO = (NFAL.GT.0 .AND. IFAL(MAX(1,NFAL)).EQ.II .AND. IFOP(MAX(1,NFAL)).GT.2)
      LSRI = (NFAL.GT.0 .AND. IFAL(MAX(1,NFAL)).EQ.II .AND. IFOP(MAX(1,NFAL)).EQ.2)
      LEIM = (NEIM.GT.0 .AND. IEIM(MAX(1,NEIM)).EQ.II)
      LJAN = (NJAN.GT.0 .AND. IJAN(MAX(1,NJAN)).EQ.II)
      LFT1 = (NFT1.GT.0 .AND. IFT1(MAX(1,NFT1)).EQ.II)
      LEXC = (NEXC.GT.0 .AND. IEXC(MAX(1,NEXC)).EQ.II)
!
      DO 500 N = 1, NSUB
         ILEN = ILASCH(SUB(N))
         KEY = ' '
!
         IF ( UPCASE(SUB(N), 3) .EQ. 'DUP') THEN
            IDUP(II) = -1
            WRITE (LOUT, 4000)
            GO TO 500
         ELSE
            I1 = INDEX(SUB(N),'/')
            I2 = INDEX(SUB(N)(I1+1:),'/')
            IF (I1.LE.0 .OR. I2.LE.0) THEN
               KERR = .TRUE.
               WRITE (LOUT, 2090) SUB(N)(:ILEN)
               GO TO 500
            ENDIF
            KEY = SUB(N)(:I1-1)
            RSTR = ' '
            RSTR = SUB(N)(I1+1:I1+I2-1)
         ENDIF
!
         IF (UPCASE(KEY, 3).EQ.'LOW' .OR. &
             UPCASE(KEY, 4).EQ.'TROE'.OR. &
             UPCASE(KEY, 3).EQ.'SRI') THEN
!
!        FALL-OFF DATA
!
            IF ((.NOT.LFAL) .OR. LLAN .OR. LRLT .OR. LREV) THEN
               KERR = .TRUE.
               IF (.NOT. LFAL) WRITE (LOUT, 1050) SUB(N)(:ILEN)
               IF (LLAN)       WRITE (LOUT, 1060) SUB(N)(:ILEN)
               IF (LRLT)       WRITE (LOUT, 1070) SUB(N)(:ILEN)
               IF (LREV)       WRITE (LOUT, 1090) SUB(N)(:ILEN)
            ELSE
!
               IF (UPCASE(KEY, 3) .EQ. 'LOW') THEN
                  IF (IFOP(NFAL) .GT. 0) THEN
                     WRITE (LOUT, 2000) SUB(N)(:ILEN)
                     KERR = .TRUE.
                  ELSE
                     IFOP(NFAL) = ABS(IFOP(NFAL))
                     CALL IPPARR (RSTR,1,3,PFAL(1,NFAL),NVAL,IER,LOUT)
                     KERR = KERR .OR. (IER.NE.0)
                     WRITE (LOUT, 3050) (PFAL(L,NFAL),L=1,3)
                  ENDIF
!
               ELSEIF (UPCASE(KEY, 4) .EQ. 'TROE') THEN
                  IF (LTRO .OR. LSRI) THEN
                     KERR = .TRUE.
                     IF (LTRO) WRITE (LOUT, 2010) SUB(N)(:ILEN)
                     IF (LSRI) WRITE (LOUT, 2030) SUB(N)(:ILEN)
                  ELSE
                     LTRO = .TRUE.
                     CALL IPPARR (RSTR,1,-4,PFAL(4,NFAL),NVAL,IER,LOUT)
                     IF (NVAL .EQ. 3) THEN
                        IFOP(NFAL) = 3*IFOP(NFAL)
                        WRITE (LOUT, 3080) (PFAL(L,NFAL),L=4,6)
                     ELSEIF (NVAL .EQ. 4) THEN
                        IFOP(NFAL) = 4*IFOP(NFAL)
                        WRITE (LOUT, 3090) (PFAL(L,NFAL),L=4,7)
                     ELSE
                        WRITE (LOUT, 2020) SUB(N)(:ILEN)
                        KERR = .TRUE.
                     ENDIF
                  ENDIF
!
               ELSEIF (UPCASE(KEY, 3) .EQ. 'SRI') THEN
                  IF (LTRO .OR. LSRI) THEN
                     KERR = .TRUE.
                     IF (LTRO) WRITE (LOUT, 2030) SUB(N)(:ILEN)
                     IF (LSRI) WRITE (LOUT, 2040) SUB(N)(:ILEN)
                  ELSE
                     LSRI = .TRUE.
                     IFOP(NFAL) = 2*IFOP(NFAL)
                     CALL IPPARR (RSTR,1,-5,PFAL(4,NFAL),NVAL,IER,LOUT)
                     IF (NVAL .EQ. 3) THEN
                        PFAL(7,NFAL) = 1.0
                        PFAL(8,NFAL) = 0.0
                        WRITE (LOUT, 3060) (PFAL(L,NFAL),L=4,6)
                     ELSEIF (NVAL .EQ. 5) THEN
                        WRITE (LOUT, 3070) (PFAL(L,NFAL),L=4,8)
                     ELSE
                        WRITE (LOUT, 2020) SUB(N)(:ILEN)
                        KERR = .TRUE.
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
!
         ELSEIF (UPCASE(KEY, 3) .EQ. 'REV') THEN
!
!        REVERSE ARRHENIUS PARAMETERS
!
            IF (LFAL .OR. LREV .OR. NSPEC(II).LT.0) THEN
               KERR = .TRUE.
               IF (LFAL) WRITE (LOUT, 1090) SUB(N)(:ILEN)
               IF (LREV) WRITE (LOUT, 2050) SUB(N)(:ILEN)
               IF (NSPEC(II) .LT. 0) WRITE (LOUT, 2060) SUB(N)(:ILEN)
            ELSE
               LREV = .TRUE.
               NREV = NREV+1
               IREV(NREV) = II
               CALL IPPARR (RSTR,1,NPAR,RPAR(1,NREV),NVAL,IER,LOUT)
               KERR = KERR .OR. (IER.NE.0)
               WRITE (LOUT, 1900) '   Reverse Arrhenius coefficients:', &
                                 (RPAR(L,NREV),L=1,3)
            ENDIF
!
         ELSEIF (UPCASE(KEY, 3) .EQ. 'RLT') THEN
!
!        REVERSE LANDAU-TELLER PARAMETERS
!
            IF (LFAL .OR. LRLT .OR. NSPEC(II).LT.0) THEN
               KERR = .TRUE.
               IF (LFAL) WRITE (LOUT, 1070) SUB(N)(:ILEN)
               IF (LRLT) WRITE (LOUT, 2080) SUB(N)(:ILEN)
               IF (NSPEC(II) .LT. 0) WRITE (LOUT, 1080) SUB(N)(:ILEN)
            ELSE
               LRLT = .TRUE.
               NRLT = NRLT + 1
               IRLT(NRLT) = II
               CALL IPPARR (RSTR,1,NLAR,RLAN(1,NRLT),NVAL,IER,LOUT)
               KERR = KERR .OR. (IER.NE.0)
               WRITE (LOUT, 3040) (RLAN(L,NRLT),L=1,2)
            ENDIF
!
         ELSEIF (UPCASE(KEY, 2) .EQ. 'HV') THEN
!
!        RADIATION WAVELENGTH ENHANCEMENT FACTOR
!
            IF (.NOT.LWL) THEN
               WRITE (LOUT, 1000) SUB(N)(:ILEN)
               KERR = .TRUE.
            ELSE
               CALL IPPARR (RSTR,1,1,VAL,NVAL,IER,LOUT)
               IF (IER .EQ. 0) THEN
                  WL(NWL) = WL(NWL)*VAL(1)
                  WRITE (LOUT, 3020) ABS(WL(NWL))
               ELSE
                  WRITE (LOUT, 1000) SUB(N)(:ILEN)
                  KERR = .TRUE.
               ENDIF
            ENDIF
!
         ELSEIF (UPCASE(KEY, 2) .EQ. 'LT') THEN
!
!        LANDAU-TELLER PARAMETERS
!
            IF (LFAL .OR. LLAN) THEN
               KERR = .TRUE.
               IF (LFAL) WRITE (LOUT, 1060) SUB(N)(:ILEN)
               IF (LLAN) WRITE (LOUT, 2070) SUB(N)(:ILEN)
            ELSE
               LLAN = .TRUE.
               NLAN = NLAN + 1
               ILAN(NLAN) = II
               CALL IPPARR (RSTR,1,NLAR,PLAN(1,NLAN),NVAL,IER,LOUT)
               IF (IER .NE. 0) THEN
                  WRITE (LOUT, 1010) SUB(N)(:ILEN)
                  KERR = .TRUE.
               ENDIF
               WRITE (LOUT, 3000) (PLAN(L,NLAN),L=1,2)
            ENDIF
!
         ELSEIF (UPCASE(KEY,4).EQ.'FORD' .OR. &
                 UPCASE(KEY,4).EQ.'RORD') THEN
             LFORD = (UPCASE(KEY,4) .EQ. 'FORD')
             LRORD = (UPCASE(KEY,4) .EQ. 'RORD')
             IF (LRORD .AND. NSPEC(II).LT.0) THEN
                KERR = .TRUE.
                WRITE (LOUT, 2065)
             ELSE
             IF (NORD.EQ.0 .OR.(NORD.GT.0 .AND. IORD(NORD).NE.II)) THEN
                NORD = NORD + 1
                IORD(NORD) = II
                NKORD = 0
!
                IF (NRNU.GT.0 .AND. IRNU(NRNU).EQ.II) THEN
                   DO 111 L = 1, 6
                      IF (NUNK(L,II) .NE. 0) THEN
                         NKORD = NKORD + 1
                         IF (RNU(L,NRNU) .LT. 0.0) THEN
                            KORD(NKORD,NORD) = -NUNK(L,II)
                            RORD(NKORD,NORD) = ABS(RNU(L,NRNU))
                         ELSE
                            KORD(NKORD,NORD) = NUNK(L,II)
                            RORD(NKORD,NORD) = RNU(L,NRNU)
                         ENDIF
                      ENDIF
  111              CONTINUE
               ELSE
                   DO 113 L = 1, 6
                      IF (NUNK(L,II) .NE. 0) THEN
                         NKORD = NKORD + 1
                         IF (NU(L,II) .LT. 0) THEN
                            KORD(NKORD,NORD) = -NUNK(L,II)
                            RORD(NKORD,NORD) =  IABS(NU(L,II))
                         ELSE
                            KORD(NKORD,NORD) = NUNK(L,II)
                            RORD(NKORD,NORD) = NU(L,II)
                         ENDIF
                      ENDIF
  113              CONTINUE
                ENDIF
             ENDIF
             ENDIF
!
             CALL IPNPAR (RSTR, 1, ISTR, ISTART)
             IF (ISTART .GE. 1) THEN
                CALL IPPARR (ISTR, 1, 1, VAL, NVAL, IER, LOUT)
                CALL CKCOMP (RSTR(:ISTART-1), KNAME, KK, K)
                IF (LFORD) K = -K
                NK = 0
                DO 121 L = 1, MAXORD
!
                   IF (KORD(L,NORD).EQ.0) THEN
                      NK = L
                      GO TO 122
                   ELSEIF (KORD(L,NORD).EQ.K) THEN
                      IF (LFORD) THEN
                         WRITE (LOUT,*)                  &
      '       Warning...changing order for reactant...', &
                         KNAME(IABS(K))
                      ELSE
                         WRITE (LOUT,*)                  &
      '       Warning...changing order for product...',  &
                         KNAME(K)
                      ENDIF
                      NK = L
                      GO TO 122
                   ENDIF
  121           CONTINUE
  122           CONTINUE
                KORD(NK,NORD) = K
                RORD(NK,NORD) = VAL(1)
                IF (LFORD) THEN
                   WRITE (LOUT, 3015) KNAME(IABS(K)),VAL(1)
                ELSE
                   WRITE (LOUT, 3016) KNAME(K),VAL(1)
                ENDIF
            ENDIF
!
         ELSEIF (UPCASE(KEY, 4) .EQ. 'PLOG') THEN
!
!        PRESSURE-DEPENDENT ARRHENIUS (PLOG) — cklink v2, stage 1.
!        One line per pressure node: PLOG / P[atm] A b E /.
!        Parsed here and accumulated in the plog_collect module (NOT
!        threaded through the CKAUXL/CKINTP argument lists). The E value
!        is stored RAW and converted to E/R[K] later in CPREAC via the
!        same EFAC as PAR(3,II). Pressure-ordering checks run in
!        plog_finalize after all entries are visible; equal-pressure
!        entries are grouped and summed by the evaluator.
!
            IF (LTHB .OR. LFAL .OR. LTRO .OR. LSRI .OR. LREV) THEN
               WRITE(*   ,'(A,I0,A)')                                  &
                  ' ERROR...PLOG reaction ', II,                       &
                  ' also uses REV/third-body/falloff syntax'
               WRITE(LOUT,'(A,I0,A)')                                  &
                  ' ERROR...PLOG reaction ', II,                       &
                  ' also uses REV/third-body/falloff syntax'
               KERR = .TRUE.
               ERROR STOP 1
            ENDIF
            CALL IPPARR (RSTR, 1, 4, PLOGV, NVAL, IER, LOUT)
            IF (IER .NE. 0 .OR. NVAL .NE. 4) THEN
               WRITE (LOUT,'(A,A)')                                     &
                  ' ERROR...PLOG requires exactly 4 values',           &
                  ' (P[atm] A b E): '//SUB(N)(:ILEN)
               KERR = .TRUE.
            ELSE
!              PLOGV = [P_atm, A, b, E_raw]; add_line validates P>0 and
!              converts P(atm)->ln(P[Pa]).
               CALL plog_add_line(II, PLOGV(1), PLOGV(2), PLOGV(3),     &
                                  PLOGV(4), KERR, LOUT)
               WRITE (LOUT,'(A,1PE10.3,A,3(1PE12.4))')                  &
                  '        PLOG: P=', PLOGV(1), ' atm  A,b,E=',         &
                  PLOGV(2), PLOGV(3), PLOGV(4)
            ENDIF
!
         ELSEIF (UPCASE(KEY, 3) .EQ. 'EIM') THEN
!
!        ELECTRON IMPACT OR THIRD-BODY REACTIONS
!
            NEIM = NEIM + 1
            IEIM(NEIM) = II
            IF (LTHB) THEN
               WRITE (LOUT, 1100) SUB(N)(:ILEN)
               KERR = .TRUE.
            ENDIF
            CALL IPPARI (RSTR, 1, 1, ITDEP(NEIM), NVAL, IER, LOUT)
            KERR = KERR .OR. (IER.NE.0) .OR. (NVAL.NE.1)
            WRITE (LOUT, 3100) ITDEP(NEIM)
!
         ELSEIF (UPCASE(KEY, 3) .EQ. 'JAN') THEN
!
!        JANNEV, LANGER, EVANS & POST TYPE REACTIONS
!
            NJAN = NJAN + 1
            IJAN(NJAN) = II
            CALL IPPARR (RSTR,1,NJAR,PJAN(1,NJAN),NVAL,IER,LOUT)
            IF (IER .NE. 0) THEN
               WRITE (LOUT, 1110) SUB(N)(:ILEN)
               KERR = .TRUE.
            ENDIF
            WRITE (LOUT, 3110) (PJAN(L,NJAN), L = 1, NJAR)
!
         ELSEIF (UPCASE(KEY, 4) .EQ. 'FIT1') THEN
!
!        MISCELLANEOUS FIT #1: k = A * T^B * exp [SUM(Vn/T^n)]
!
            NFT1 = NFT1 + 1
            IFT1(NFT1) = II
            CALL IPPARR (RSTR,1,NF1R,PFT1(1,NFT1),NVAL,IER,LOUT)
            IF (IER .NE. 0) THEN
               WRITE (LOUT, 1112) SUB(N)(:ILEN)
               KERR = .TRUE.
            ENDIF
            WRITE (LOUT, 3112) (PFT1(L,NFT1), L = 1, NF1R)
!
         ELSEIF (UPCASE(KEY, 4) .EQ. 'EXCI') THEN
!
!        EXCITATION-ONLY REACTION DESCRIPTION (FOR ENERGY LOSS)
!
            NEXC = NEXC + 1
            IEXC(NEXC) = II
            CALL IPPARR (RSTR,1,1,PEXC(NEXC),NVAL,IER,LOUT)
            KERR = KERR .OR. (IER.NE.0) .OR. (NVAL.NE.1)
            WRITE (LOUT, 3114) PEXC(NEXC)
!
         ELSE
!
!        ENHANCED THIRD BODIES
!
            CALL CKCOMP (KEY, KNAME, KK, K)
            IF (K .EQ. 0) THEN
               WRITE (LOUT, 1040) KEY(:ILASCH(KEY))
               KERR = .TRUE.
            ELSE
               IF (.NOT.LTHB) THEN
                  KERR = .TRUE.
                  WRITE (LOUT, 1020) SUB(N)(:ILEN)
               ELSE
                  IF (NTBS(NTHB) .EQ. MAXTB) THEN
                     KERR = .TRUE.
                     WRITE (LOUT, 1030) SUB(N)(:ILEN)
                  ELSE
                     CALL IPPARR (RSTR, 1, 1, VAL, NVAL, IER, LOUT)
                     IF (IER .EQ. 0) THEN
                        WRITE (LOUT, 3010) KNAME(K),VAL(1)
                        NTBS(NTHB) = NTBS(NTHB) + 1
                        NKTB(NTBS(NTHB),NTHB) = K
                        AIK(NTBS(NTHB),NTHB) = VAL(1)
                     ELSE
                        WRITE (LOUT, 1020) SUB(N)(:ILEN)
                        KERR = .TRUE.
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
  500 CONTINUE
!
!     FORMATS
!
 1000 FORMAT (6X,'Error in HV declaration...',A)
 1010 FORMAT (6X,'Error in LT declaration..',A)
 1020 FORMAT (6X,'Error in 3rd-body declaration...',A)
 1030 FORMAT (6X,'Error...more than MAXTB 3rd bodies...',A)
 1040 FORMAT (6X,'Error...undeclared species...',A)
 1050 FORMAT (6X,'Error...this is not a fall-off reaction...',A)
 1060 FORMAT (6X,'Error...LT declared in fall-off reaction...',A)
 1070 FORMAT (6X,'Error...RLT declared in fall-off reaction...',A)
 1080 FORMAT (6X,'Error...RLT declared in irreversible reaction...',A)
 1090 FORMAT (6X,'Error...REV declared in fall-off reaction...',A)
 1100 FORMAT (6X,'Error...EIM declared in heavy 3rd-body reaction...',A)
 1110 FORMAT (6X,'Error in JAN declaration...',A)
 1112 FORMAT (6X,'Error in FIT1 declaration...',A)
 2000 FORMAT (6X,'Error...LOW declared more than once...',A)
 2010 FORMAT (6X,'Error...TROE declared more than once...',A)
 2020 FORMAT (6X,'Error in fall-off parameters...',A)
 2030 FORMAT (6X,'Error...cannot use both TROE and SRI...',A)
 2040 FORMAT (6X,'Error...SRI declared more than once...',A)
 2050 FORMAT (6X,'Error...REV declared more than once...',A)
 2060 FORMAT (6X,'Error...REV declared for irreversible reaction...',A)
 2065 FORMAT (6X,'Error...RORD declared for irreversible reaction...')
 2070 FORMAT (6X,'Error...LT declared more than once...',A)
 2080 FORMAT (6X,'Error...RLT declared more than once...',A)
 2090 FORMAT (6X,'Error in auxiliary data...',A)
 3000 FORMAT (9X,'Landau-Teller parameters: B=',E12.5,', !=',E12.5)
 3010 FORMAT (9X,A16,' Enhanced by ',1PE12.3)
 3015 FORMAT (7X,A16,' Forward order ',1PE12.3)
 3016 FORMAT (7X,A16,' Reverse order ',1PE12.3)
 3020 FORMAT (9X,'Radiation wavelength (A): ',F10.2)
! 1900 FORMAT (6X,A,T51,E10.3,F7.3,F11.3)
 1900 FORMAT (6X, A, T53, 1PE8.2, 2X, 0PF5.1, 2X, F9.1)
 3040 FORMAT (9X,'Reverse Landau-Teller parameters: B=',E12.5,', !=',E12.5)
 3050 FORMAT (6X,'Low pressure limit:',3E13.5)
 3060 FORMAT (6X,'SRI centering:     ',3E13.5)
 3070 FORMAT (6X,'SRI centering:     ',5E13.5)
 3080 FORMAT (6X,'TROE centering:    ',3E13.5)
 3090 FORMAT (6X,'TROE centering:    ',4E13.5)
 3100 FORMAT (6X,'Electron 3rd-body reaction; Temp. Dependence =',I5)
 3110 FORMAT (6X,'Jannev, Langer, Evans & Post type reaction:' &
             /9X,'Coefficients: ',5(E10.3,1X)                  &
            /23X,5(E10.3,1X))
 3112 FORMAT (6X,'Modified fit#1:  k= A * T^B * exp [SUM(Vn/T^n)]...', &
             /9X,'Added parameters: ',E10.3,3(/27X,E10.3))
 3114 FORMAT (6X,'Excitation reaction only; Energy loss =',e10.3,' eV')
 4000 FORMAT (6X,'Declared duplicate reaction...')
      END SUBROUTINE CKAUXL

!-------------------------------------------------------------------------

!----------------------------------------------------------------------!
      SUBROUTINE CKPRNT (MDIM, MAXTP, MM, ENAME, KK, KNAME, WTM,    &
                         KPHSE, KCHRG, NT, T, TLO, TMID, THI, KNCF, &
                         ITHRM, LOUT, KERR)
!
!     Prints species interpreter output and checks for completeness.
!----------------------------------------------------------------------!
      USE chemkinII, only: ILASCH
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
!
      DIMENSION WTM(*), KPHSE(*), KCHRG(*), T(MAXTP,*), &
                NT(*), KNCF(MDIM,*), IPLUS(10)
      LOGICAL KERR, ITHRM(*)
      CHARACTER ENAME(*)*(*), KNAME(*)*(*), IPHSE(3)*1, INUM(10)*1
      DATA IPHSE/'S','G','L'/
      DATA INUM/'0','1','2','3','4','5','6','7','8','9'/
!
      WRITE (LOUT, 400) (ENAME(M), M = 1, MM)
      WRITE (LOUT, 300)
!
      DO 100 K = 1, KK
!
         IF (T(1,K) .LT. 0.0) T(1,K) = TLO
         IF (T(2,K) .LT. 0.0) T(2,K) = TMID
         IF (T(3,K) .LT. 0.0) T(NT(K),K) = THI
         WRITE (LOUT, 500) K, KNAME(K), IPHSE(KPHSE(K)+2), KCHRG(K), &
                          WTM(K), INT(T(1,K)), INT(T(NT(K),K)),      &
                         (KNCF(M,K),M=1,MM)
         IF (T(1,K) .GE. T(NT(K),K)) THEN
            KERR = .TRUE.
            WRITE (LOUT, 240)
         ENDIF
         IF (T(1,K) .GT. T(2,K)) THEN
            WRITE (LOUT, 250)
            KERR = .TRUE.
         ENDIF
         IF (T(NT(K),K) .LT. T(2,K)) THEN
            WRITE (LOUT, 260)
            KERR = .TRUE.
         ENDIF
!
!        each species must have thermodynamic data
!
         IF (.NOT. ITHRM(K)) THEN
            KERR = .TRUE.
            WRITE (LOUT, 200)
         ENDIF
!
!        a species cannot start with a number
!
         CALL CKCOMP (KNAME(K)(:1), INUM, 10, I)
         IF (I .GT. 0) THEN
            KERR = .TRUE.
            WRITE (LOUT, 210)
         ENDIF
!
!        if '+' sign is used in a species name,
!           examples of legal species symbols with + are:
!           OH(+)2, OH(+2), OH+, OH++, OH+++, OH(+), OH(++),
!           OH[+OH], OH2+, OH+2
!
!           examples of illegal species symbols with + are:
!           +OH        (symbol starts with a +, this will cause
!                       confusion in a reaction)
!           OH(+OH)    (symbol in parentheses is another species-
!                       this arrangement is reserved for a fall-off
!                       reaction)
!           OH+OH      (plus delimits other species names, this
!                       will cause confusion in a reaction)
!
         NPLUS = 0
         DO 50 N = 1, ILASCH(KNAME(K))
            IF (KNAME(K)(N:N) .EQ. '+') THEN
               NPLUS = NPLUS + 1
               IPLUS(NPLUS) = N
            ENDIF
   50    CONTINUE
         DO 60 N = 1, NPLUS
            I1 = IPLUS(N)
            IF (I1 .EQ. 1) THEN
               WRITE (LOUT, 220)
               KERR = .TRUE.
            ELSE
!
!              is there another species name in parentheses
!
               IF (KNAME(K)(I1-1:I1-1) .EQ. '(') THEN
                  I1 = I1 + 1
                  I2 = I1 + INDEX(KNAME(K)(I1:),')')-1
                  IF (I2 .GT. I1) THEN
                     CALL CKCOMP (KNAME(K)(I1:I2-1), KNAME, KK, KNUM)
                     IF (KNUM .GT. 0) THEN
                        WRITE (LOUT, 230)
                        KERR = .TRUE.
                     ENDIF
                  ENDIF
               ENDIF
!
!              is there another species name after a +
!
               I1 = I1 + 1
               IF (N .LT. NPLUS) THEN
                  DO 55 L = N+1, NPLUS
                     I2 = IPLUS(L)
                     IF (I2 .GT. I1) THEN
                        CALL CKCOMP (KNAME(K)(I1:I2-1),KNAME,KK,KNUM)
                        IF (KNUM .GT. 0) THEN
                           WRITE (LOUT, 230)
                           KERR = .TRUE.
                        ENDIF
                     ENDIF
   55             CONTINUE
               ENDIF
!
               I2 = ILASCH(KNAME(K))
               IF (I2 .GE. I1) THEN
                  CALL CKCOMP (KNAME(K)(I1:I2), KNAME, KK, KNUM)
                  IF (KNUM .GT. 0) THEN
                     WRITE (LOUT, 230)
                     KERR = .TRUE.
                  ENDIF
               ENDIF
            ENDIF
   60    CONTINUE
!
  100 CONTINUE
      WRITE (LOUT, 300)
      RETURN
!
  200 FORMAT (6X,'Error...no thermodynamic properties for species')
  210 FORMAT (6X,'Error...species starts with a number')
  220 FORMAT (6X,'Error...species starts with a plus')
  230 FORMAT (6X,'Error...illegal + in species name')
  240 FORMAT (6X,'Error...High temperature must be < Low temperature')
  250 FORMAT (6X,'Error...Low temperature must be <= Mid temperature')
  260 FORMAT (6X,'Error...High temperature must be => Mid temperature')
  300 FORMAT (1X,79('-'))
!
  400 FORMAT (1X,79('-'),/T27,'C',/T24,'P  H',/T24,'H  A',/T24,'A  R',  &
             /1X,'SPECIES',T24,'S  G',T30,'MOLECULAR',T41,'TEMPERATURE',&
             T54,'ELEMENT COUNT',                                       &
             /1X,'CONSIDERED',T24,'E  E',T30,'WEIGHT',T41,'LOW',        &
             T48,'HIGH',T54,15(A3))
  500 FORMAT (1X,I3,'. ',A16,T24,A1,T26,I2,T29,F10.5,T39,I6,T46,I6, &
             T53,15(I3))
      END SUBROUTINE CKPRNT
!----------------------------------------------------------------------!

      SUBROUTINE CPREAC (II, MAXSP, NSPEC, NPAR, PAR, RPAR, AUNITS,     &
                         EUNITS, NREAC, NUNK, NU, KCHRG, MDIM, MM, KNCF,&
                         IDUP, NFAL, IFAL, KFAL, NFAR, PFAL, IFOP, NREV,&
                         IREV, NTHB, ITHB, NLAN, ILAN, NRLT, IRLT, KERR,&
                         LOUT, NRNU, IRNU, RNU, CKMIN)
!
!     Prints reaction interpreter output and checks for reaction
!     balance, duplication, and missing data in 'REV' reactions;
!     correct units of Arrhenius parameters
!
!     Input: II     - the index number of the reaction
!            MAXSP  - maximum number of species allowed in a reaction
!            NSPEC  - array of the number of species in the reactions
!            NPAR   - the number of Arrhenius parameters required
!            PAR    - matrix of Arrhenius parameters for the reactions
!            RPAR   - matrix of reverse Arrhenius parameters for the
!                     reactions which declared them
!            AUNITS - character string which describes the input units
!                     of A, the pre-exponential factor PAR(1,I)
!            EUNITS - character string which describes the input units
!                     of E, the activation energy PAR(3,I)
!            NREAC  - array of the number of reactants in the reactions
!            NUNK   - matrix of the species indices of the reactants
!                     and products in the reactions
!            NU     - matrix of the stoichiometric coefficients of the
!                     reactants and products in the reactions
!            KCHRG  - array of the electronic charges of the species
!            MDIM   - the maximum number of elements allowed
!            MM     - the actual number of elements declared
!            KNCF   - matrix of elemental composition of the species
!            IDUP   - array of integer flags to indicate duplicate
!                     reactions
!            NFAL   - total number of reactions with fall-off
!            IFAL   - the NFAL reaction indices
!            NFAR   - maximum number of fall-off parameters allowed
!            PFAL   - matrix of fall-off parameters for the NFAL
!                     reactions
!            IFOP   - the NFAL fall-off types
!            NREV   - total number of reactions with reverse parameters
!            IREV   - the NREV reaction indices
!            NTHB   - total number of reactions with third-bodies
!            ITHB   - the NTHB reaction indices
!            NLAN   - total number of reactions with Landauer-Teller
!                     parameters
!            ILAN   - the NLAN reaction indices
!            NRLT   - total number of reactions with reverse
!                     Landauer-Teller parameters
!            IRLT   - the NRLT reaction indices
!            KERR   - logical error flag
!            LOUT   - unit number for output messages
!
!----------------------------------------------------------------------!
!     (Value of Avrogadro's Constant from 1986 CODATA
!      recommended values (1993 CRC)
!      J. Research National Bureal of Standards, 92, 95, 1987
!      6.0221367(39) mol-1 )
!
      USE plog_collect, only: plog_apply_efac, plog_apply_afac
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
      DOUBLE PRECISION RU_JOUL,        AVAG,                ONE
      PARAMETER (RU_JOUL = 8.3140D0, AVAG = 6.0221367D23, ONE=1.0D0)

      DIMENSION NSPEC(*), PAR(NPAR,*), RPAR(NPAR,*), NREAC(*),      &
                NUNK(MAXSP,*), NU(MAXSP,*), KCHRG(*), KNCF(MDIM,*), &
                IDUP(*), IFAL(*), KFAL(*), PFAL(NFAR,*), IFOP(*),   &
                IREV(*), ITHB(*), ILAN(*), IRLT(*), IRNU(*),        &
                RNU(MAXSP,*)
      CHARACTER*(*) AUNITS, EUNITS
      LOGICAL IERR,KERR,LREV,LLAN,LRLT,LREAL
!
!     Fortran does not guarantee short-circuit evaluation of .AND.;
!     do not form IRNU(0) when no real-stoichiometry reaction exists.
      LREAL = .FALSE.
      IF (NRNU.GT.0) LREAL = II.EQ.IRNU(NRNU)
      IF (LREAL) THEN
         CALL CKRBAL (MAXSP, NUNK(1,II), RNU(1,NRNU), MDIM, MM, KCHRG,&
                      KNCF, CKMIN, IERR)
      ELSE
         CALL CKBAL (MAXSP, NUNK(1,II), NU(1,II), MDIM, MM, KCHRG, KNCF,&
                     IERR)
      ENDIF
!
      IF (IERR) THEN
         KERR = .TRUE.
         WRITE (LOUT, 1060)
      ENDIF
!
      CALL CKDUP (II, MAXSP, NSPEC, NREAC, NU, NUNK, NFAL, IFAL, KFAL,&
                  ISAME)
!
      IF (ISAME .GT. 0) THEN
         IF (IDUP(ISAME).NE.0 .AND. IDUP(II).NE.0) THEN
            IDUP(ISAME) = ABS(IDUP(ISAME))
            IDUP(II)    = ABS(IDUP(II))
         ELSE
            N1 = 0
            N2 = 0
            IF (NTHB .GT. 1) THEN
               DO 150 N = 1, NTHB
                  IF (ITHB(N) .EQ. ISAME) N1 = 1
                  IF (ITHB(N) .EQ. II)    N2 = 1
  150          CONTINUE
            ENDIF
            IF (N1 .EQ. N2) THEN
               KERR = .TRUE.
               WRITE (LOUT, 1050) ISAME
            ENDIF
         ENDIF
      ENDIF
!

      IF (NFAL.GT.0) THEN
         IF (IFAL(NFAL).EQ.II .AND. IFOP(NFAL).LT.0) THEN
           KERR = .TRUE.
           WRITE (LOUT, 1020)
         ENDIF
      ENDIF
!
      LREV = (NREV.GT.0 .AND. IREV(MAX(1,NREV)).EQ.II)
      LLAN = (NLAN.GT.0 .AND. ILAN(MAX(1,NLAN)).EQ.II)
      LRLT = (NRLT.GT.0 .AND. IRLT(MAX(1,NRLT)).EQ.II)
      IF (LREV .AND. LLAN .AND. (.NOT.LRLT)) THEN
         KERR = .TRUE.
         WRITE (LOUT, 1030)
      ENDIF
      IF (LRLT .AND. (.NOT.LLAN)) THEN
         KERR = .TRUE.
         WRITE (LOUT, 1040)
      ENDIF
      IF (LRLT .AND. (.NOT.LREV)) THEN
         KERR = .TRUE.
         WRITE (LOUT, 1045)
      ENDIF
!
      IF (EUNITS .EQ. 'KELV') THEN
         EFAC = 1.0
      ELSEIF (EUNITS .EQ. 'CAL/') THEN
!        convert E from cal/mole to Kelvin
         EFAC = 4.184  / RU_JOUL
      ELSEIF (EUNITS .EQ. 'KCAL') THEN
!        convert E from kcal/mole to Kelvin
         EFAC = 4184.0 / RU_JOUL
      ELSEIF (EUNITS .EQ. 'JOUL') THEN
!        convert E from Joules/mole to Kelvin
         EFAC = 1.00  / RU_JOUL
      ELSEIF (EUNITS .EQ. 'KJOU') THEN
!        convert E from Kjoules/mole to Kelvin
!        1 kJ/mol = 1000 J/mol (was erroneously 4000.0, a 4x error in E/R)
         EFAC = 1000.0 / RU_JOUL
      ENDIF
      PAR(3,II) = PAR(3,II) * EFAC
!
!     PLOG activation energies use the identical EFAC (E -> E/R [K]).
      CALL plog_apply_efac(II, EFAC)
!
!      IF (NREV.GT.0 .AND. IREV(NREV).EQ.II) RPAR(3,II)=RPAR(3,II)*EFAC
!      IF (NFAL.GT.0 .AND. IFAL(NFAL).EQ.II) PFAL(3,II)=PFAL(3,II)*EFAC
!
      IF (NREV.GT.0 .AND. IREV(MAX(1,NREV)).EQ.II) &
          RPAR(3,NREV) = RPAR(3,NREV) * EFAC
      IF (NFAL.GT.0 .AND. IFAL(MAX(1,NFAL)).EQ.II) &
          PFAL(3,NFAL) = PFAL(3,NFAL) * EFAC
!
      IF (AUNITS .EQ. 'MOLC') THEN
         NSTOR = 0
         NSTOP = 0
         DO 50 N = 1, MAXSP
            IF (NU(N,II) .LT. 0) THEN
!              sum of stoichiometric coefficients of reactants
               NSTOR = NSTOR + ABS(NU(N,II))
            ELSEIF (NU(N,II) .GT. 0) THEN
!              sum of stoichiometric coefficients of products
               NSTOP = NSTOP + NU(N,II)
            ENDIF
   50    CONTINUE
!        PLOG is valid only for an ordinary gas-phase reaction (the
!        combination guard runs when cklink is loaded). Its A-factor
!        conversion therefore uses the unmodified forward molecularity.
         IF (NSTOR.GT.0) CALL plog_apply_afac(II, AVAG**(NSTOR-1))
!
         IF (NFAL.GT.0 .AND. IFAL(NFAL).EQ.II) THEN
!
!           fall-off reaction, "(+M)" or "(+species name)" does not
!           count except in "LOW" A-factor;
!           reverse-rate declarations are not allowed
!
            IF (NSTOR.GT.0) PAR(1,II) = PAR(1,II) * AVAG**(NSTOR-1)
            NSTOR = NSTOR + 1
            IF (NSTOR.GT.0) PFAL(1,NFAL) = PFAL(1,NFAL)*AVAG**(NSTOR-1)
!
         ELSEIF (NTHB.GT.0 .AND. ITHB(NTHB).EQ.II) THEN
!
!           third body reaction, "+M" counts as species in
!           forward and reverse A-factor conversion
!
            NSTOR = NSTOR + 1
            NSTOP = NSTOP + 1
            IF (NSTOR.GT.0) PAR(1,II) = PAR(1,II) * AVAG**(NSTOR-1)
            IF (NREV.GT.0 .AND. IREV(NREV).EQ.II .AND. NSTOP.GT.0) &
                RPAR(1,NREV) = RPAR(1,NREV) * AVAG**(NSTOP-1)
!
         ELSE
!
!           not third-body or fall-off reaction, but may have
!           reverse rates.
!
            IF (NSTOR .GT. 0) PAR(1,II) = PAR(1,II) * AVAG**(NSTOR-1)
            IF (NREV.GT.0 .AND. IREV(NREV).EQ.II .AND. NSTOP.GT.0) &
                RPAR(1,NREV) = RPAR(1,NREV) * AVAG**(NSTOP-1)
         ENDIF
      ENDIF
!
 1020 FORMAT (6X,'Error...no LOW parameters given for fall-off...')
 1030 FORMAT (6X,'Error...reverse T-L required...')
 1040 FORMAT (6X,'Error...forward T-L required...')
 1045 FORMAT (6X,'Error...REV parameters must be given with RTL...')
 1050 FORMAT (6X,'Error...undeclared duplicate to reaction number ',I3)
 1060 FORMAT (6X,'Error...reaction does not balance...')
      RETURN
      END SUBROUTINE CPREAC

!-------------------------------------------------------------------------------

!----------------------------------------------------------------------!
      SUBROUTINE CKBAL (MXSPEC, KSPEC, KCOEF, MDIM, MM, KCHRG, KNCF,   &
                        IERR)
!
!     Checks elemental balance of reactants vs. products.
!     Checks charge balance of reaction.
!
!     Input:  MXSPEC - number of species allowed in a reaction
!             KSPEC(N),N=1,MXSPEC- array of species indices in reaction
!             KCOEF(N) - stoichiometric coefficients of the species
!             MDIM  - maximum number of elements allowed
!             MM    - actual integer number of elements
!             KCHRG(K) - ionic charge Kth species
!             KNCF(M,K)- integer elemental composition of Kth species
!     Output: KERR  - logical, =.TRUE. if reaction does not balance
!                                      F. Rupley, Div. 8245, 5/13/86
!----------------------------------------------------------------------!

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

!
      DIMENSION KSPEC(*), KCOEF(*), KNCF(MDIM,*), KCHRG(*)
      LOGICAL IERR
!
      IERR = .FALSE.
!
!     charge balance
!
      KBAL = 0
      DO 50 N = 1, MXSPEC
         IF (KSPEC(N) .NE. 0) &
         KBAL = KBAL + KCOEF(N)*KCHRG(KSPEC(N))
   50 CONTINUE
      IF (KBAL .NE. 0) IERR = .TRUE.
!
!     element balance
!
      DO 100 M = 1, MM
         MBAL = 0
         DO 80 N = 1, MXSPEC
            IF (KSPEC(N) .NE. 0) &
            MBAL = MBAL + KCOEF(N)*KNCF(M,KSPEC(N))
   80    CONTINUE
         IF (MBAL .NE. 0) IERR = .TRUE.
  100 CONTINUE
      RETURN
      END SUBROUTINE CKBAL

!-----------------------------------------------------------------------!

!----------------------------------------------------------------------!
      SUBROUTINE CKRBAL (MXSPEC, KSPEC, RCOEF, MDIM, MM, KCHRG, KNCF, &
                         CKMIN, IERR)
!
!     Checks elemental balance of reactants vs. products.
!     Checks charge balance of reaction.
!
!     Input:  MXSPEC - number of species allowed in a reaction
!             KSPEC(N),N=1,MXSPEC- array of species indices in reaction
!             RCOEF(N) - stoichiometric coefficients of the species
!             MDIM  - maximum number of elements allowed
!             MM    - actual integer number of elements
!             KCHRG(K) - ionic charge Kth species
!             KNCF(M,K)- integer elemental composition of Kth species
!     Output: KERR  - logical, =.TRUE. if reaction does not balance
!                                      F. Rupley, Div. 8245, 5/13/86
!----------------------------------------------------------------------!

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

      DIMENSION KSPEC(*), RCOEF(*), KNCF(MDIM,*), KCHRG(*)
      LOGICAL IERR
!
      IERR = .FALSE.
!
!     charge balance
!
      SBAL = 0
      DO 50 N = 1, MXSPEC
         IF (KSPEC(N) .NE. 0)&
         SBAL = SBAL + RCOEF(N)*KCHRG(KSPEC(N))
   50 CONTINUE
      IF (ABS(SBAL) .GT. CKMIN) IERR = .TRUE.
!
!     element balance
!
      DO 100 M = 1, MM
         SMBAL = 0
         DO 80 N = 1, MXSPEC
            IF (KSPEC(N) .NE. 0) &
            SMBAL = SMBAL + RCOEF(N)*KNCF(M,KSPEC(N))
   80    CONTINUE
         IF (ABS(SMBAL) .GT. CKMIN) IERR = .TRUE.
  100 CONTINUE
      RETURN
      END SUBROUTINE CKRBAL

!----------------------------------------------------------------------!
      SUBROUTINE CKDUP (I, MAXSP, NS, NR, NU, NUNK, NFAL, IFAL, KFAL, &
                        ISAME)
!
!     Checks reaction I against the (I-1) reactions for duplication
!----------------------------------------------------------------------!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

      DIMENSION NS(*), NR(*), NU(MAXSP,*), NUNK(MAXSP,*), IFAL(*), &
                KFAL(*)
!
      ISAME = 0
      NRI = NR(I)
      NPI = ABS(NS(I)) - NR(I)
!
      DO 500 J = 1, I-1
!
         NRJ = NR(J)
         NPJ = ABS(NS(J)) - NR(J)
!
         IF (NRJ.EQ.NRI .AND. NPJ.EQ.NPI) THEN
!
            NSAME = 0
            DO 20 N = 1, MAXSP
               KI = NUNK(N,I)
               NI = NU(N,I)
!
               DO 15 L = 1, MAXSP
                  KJ = NUNK(L,J)
                  NJ = NU(L,J)
                  IF (NJ.NE.0 .AND. KJ.EQ.KI .AND. NJ.EQ.NI) &
                  NSAME = NSAME + 1
   15          CONTINUE
   20       CONTINUE
!
            IF (NSAME .EQ. ABS(NS(J))) THEN
!
!           same products, reactants, coefficients, check fall-off
!           third body
!
               IF (NFAL.GT.0 .AND. IFAL(NFAL).EQ.I) THEN
                  DO 22 N = 1, NFAL-1
                     IF (J.EQ.IFAL(N) .AND. KFAL(N).EQ.KFAL(NFAL)) THEN
                        ISAME = J
                        RETURN
                     ENDIF
   22             CONTINUE
                  RETURN
               ENDIF
!
               ISAME = J
               RETURN
            ENDIF
         ENDIF
!
         IF (NPI.EQ.NRJ .AND. NPJ.EQ.NRI) THEN
!
            NSAME = 0
            DO 30 N = 1, MAXSP
               KI = NUNK(N,I)
               NI = NU(N,I)
!
               DO 25 L = 1, MAXSP
                  KJ = NUNK(L,J)
                  NJ = NU(L,J)
                  IF (NJ.NE.0 .AND. KJ.EQ.KI .AND. -NJ.EQ.NI) &
                  NSAME = NSAME + 1
   25          CONTINUE
   30       CONTINUE
!
            IF (NSAME.EQ.ABS(NS(J)) .AND.   &
                (NS(J).GT.0 .OR. NS(I).GT.0)) THEN
!
!           same products as J reactants, and vice-versa
!
               IF (NFAL.GT.0 .AND. IFAL(NFAL).EQ.I) THEN
                  DO 32 N = 1, NFAL-1
                     IF (J.EQ.IFAL(N) .AND. KFAL(N).EQ.KFAL(NFAL)) THEN
                        ISAME = J
                        RETURN
                     ENDIF
   32             CONTINUE
                  RETURN
               ENDIF
!
               ISAME = J
               RETURN
            ENDIF
         ENDIF
!
  500 CONTINUE
      RETURN
      END SUBROUTINE CKDUP


!----------------------------------------------------------------------!
      SUBROUTINE CKISUB (LINE, SUB, NSUB)
!
!     Generates an array of CHAR*(*) substrings from a CHAR*(*) string,
!     using blanks or tabs as delimiters
!
!     Input:  LINE  - a CHAR*(*) line
!     Output: SUB   - a CHAR*(*) array of substrings
!             NSUB  - number of substrings found
!     A '!' will comment out a line, or remainder of the line.
!                                      F. Rupley, Div. 8245, 5/15/86
!----------------------------------------------------------------------!
      USE chemkinII, only: IPPLEN, ILASCH, IFIRCH
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
!
      CHARACTER*(*) SUB(*), LINE
      NSUB = 0
!
      DO 5 N = 1, LEN(LINE)
        IF (ICHAR(LINE(N:N)) .EQ. 9) LINE(N:N) = ' '
    5 CONTINUE
!
      IF (IPPLEN(LINE) .LE. 0) RETURN
!
      ILEN = ILASCH(LINE)
!
      NSTART = IFIRCH(LINE)
   10 CONTINUE
      ISTART = NSTART
      NSUB = NSUB + 1
      SUB(NSUB) = ' '
!
      DO 100 I = ISTART, ILEN
         ILAST = INDEX(LINE(ISTART:),' ') - 1
         IF (ILAST .GT. 0) THEN
            ILAST = ISTART + ILAST - 1
         ELSE
            ILAST = ILEN
         ENDIF
         SUB(NSUB) = LINE(ISTART:ILAST)
         IF (ILAST .EQ. ILEN) RETURN
!
         NSTART = ILAST + IFIRCH(LINE(ILAST+1:))
!
!        Does SUB have any slashes?
!
         I1 = INDEX(SUB(NSUB),'/')
         IF (I1 .LE. 0) THEN
            IF (LINE(NSTART:NSTART) .NE. '/') GO TO 10
            NEND = NSTART + INDEX(LINE(NSTART+1:),'/')
            IND = INDEX(SUB(NSUB),' ')
            SUB(NSUB)(IND:) = LINE(NSTART:NEND)
            IF (NEND .EQ. ILEN) RETURN
            NSTART = NEND + IFIRCH(LINE(NEND+1:))
            GO TO 10
         ENDIF
!
!        Does SUB have 2 slashes?
!
         I2 = INDEX(SUB(NSUB)(I1+1:),'/')
         IF (I2 .GT. 0) GO TO 10
!
         NEND = NSTART + INDEX(LINE(NSTART+1:),'/')
         IND = INDEX(SUB(NSUB),' ') + 1
         SUB(NSUB)(IND:) = LINE(NSTART:NEND)
         IF (NEND .EQ. ILEN) RETURN
         NSTART = NEND + IFIRCH(LINE(NEND+1:))
         GO TO 10
  100 CONTINUE
      RETURN
      END SUBROUTINE CKISUB
!----------------------------------------------------------------------!
      SUBROUTINE IPNPAR (LINE, NPAR, IPAR, ISTART)
!
!     Returns CHAR*(*) IPAR substring of CHAR*(*) string LINE which
!     contains NPAR real parameters
!
!     Input:     LINE - a CHAR*(*) line
!                NPAR - number of parameters expected
!     Output:    IPAR - the substring of parameters only
!                ISTART - the starting location of IPAR substring
!     A '!' will comment out a line, or remainder of the line.
!                                      F. Rupley, Div. 8245, 5/14/86
!----------------------------------------------------------------------!
!*****precision > double
      USE chemkinII, only : IPPLEN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
!*****END precision > double
!*****precision > single
!      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
!*****END precision > single
!
      CHARACTER*(*) LINE,IPAR
!
!----------Find Comment String (! signifies comment)
!
      ILEN = IPPLEN(LINE)
      ISTART = 0
      N = 0
      IF (ILEN.GT.0) THEN
         DO 40 I = ILEN, 1, -1
            ISTART = I
            IPAR = ' '
            IPAR = LINE(ISTART:ILEN)
            IF (LINE(I:I).NE.' ') THEN
               IF (I .EQ. 1) RETURN
               IF (LINE(I-1:I-1) .EQ. ' ') THEN
                  N = N + 1
                  IF (N .EQ. NPAR) RETURN
               ENDIF
            ENDIF
   40    CONTINUE
      ENDIF
      RETURN
      END SUBROUTINE IPNPAR

!----------------------------------------------------------------------!
      SUBROUTINE CKUNIT (LINE, AUNITS, EUNITS, IUNITS)
      USE chemkinII, only: ILASCH, UPCASE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

      CHARACTER*(*) LINE, IUNITS, AUNITS, EUNITS
!      CHARACTER*4 UPCASE
!      DIMENSION UPCASE(4)
!
      AUNITS = ' '
      EUNITS = ' '
      IUNITS = ' '
      LCHAR = ILASCH(LINE)
      DO 85 N = 1, ILASCH(LINE)-3
         IND = ILASCH(IUNITS)
         IF (EUNITS .EQ. ' ') THEN
            IF (UPCASE(LINE(N:), 4)     .EQ. 'CAL/') THEN
               EUNITS = 'CAL/'
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units cal/mole'
               ELSE
                  IUNITS(IND:) = ', E units cal/mole'
               ENDIF
            ELSEIF (UPCASE(LINE(N:), 4) .EQ. 'KCAL') THEN
               EUNITS = 'KCAL'
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units Kcal/mole'
               ELSE
                  IUNITS(IND:) = ', E units Kcal/mole'
               ENDIF
            ELSEIF (UPCASE(LINE(N:), 4) .EQ. 'JOUL') THEN
               EUNITS = 'JOUL'
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units Joules/mole'
               ELSE
                  IUNITS(IND:) = ', E units Joules/mole'
               ENDIF
            ELSEIF (UPCASE(LINE(N:), 4) .EQ. 'KJOU') THEN
               EUNITS = 'KJOU'
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units Kjoule/mole'
               ELSE
                  IUNITS(IND:) = ', E units Kjoule/mole'
               ENDIF
            ELSEIF (UPCASE(LINE(N:), 4) .EQ. 'KELV') THEN
               EUNITS = 'KELV'
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units Kelvins'
               ELSE
                  IUNITS(IND:) = ', E units Kelvins'
               ENDIF
            ENDIF
         ENDIF
         IF (AUNITS .EQ. ' ') THEN
            IF (UPCASE(LINE(N:), 4) .EQ. 'MOLE') THEN
               IF (N+4.LE.ILASCH(LINE) .AND. &
                          UPCASE(LINE(N+4:),1).EQ.'C') THEN
!
                  AUNITS = 'MOLC'
                  IF (IUNITS .EQ. ' ') THEN
                     IUNITS = 'A units molecules'
                  ELSE
                      IUNITS(IND:) = ', A units molecules'
                  ENDIF
               ELSE
                  AUNITS = 'MOLE'
                  IF (IUNITS .EQ. ' ') THEN
                     IUNITS = 'A units mole-cm-sec-K'
                  ELSE
                     IUNITS(IND:) = ', A units mole-cm-sec-K'
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
   85 CONTINUE
!
      IF (AUNITS .EQ. ' ') THEN
         AUNITS = 'MOLE'
         IND = ILASCH(IUNITS) + 1
         IF (IND .GT. 1) THEN
            IUNITS(IND:) = ', A units mole-cm-sec-K'
         ELSE
            IUNITS(IND:) = ' A units mole-cm-sec-K'
         ENDIF
      ENDIF
!
      IF (EUNITS .EQ. ' ') THEN
         EUNITS = 'CAL/'
         IND = ILASCH(IUNITS) + 1
         IF (IND .GT. 1) THEN
            IUNITS(IND:) = ', E units cal/mole'
         ELSE
            IUNITS(IND:) = ' E units cal/mole'
         ENDIF
      ENDIF
!
      RETURN
      END SUBROUTINE CKUNIT
!

      end module chemkinII_interpreter

!     *****************************************************************
!     **                                                             **
!     **                     KIVA4 - CHEMISTRY                       **
!     **                                                             **
!     **               CHEMKIN-II - KIVA4 INTERFACE                  **
!     **                                                             **
!     **                                                             **
!     **   Modified by: Federico Perini                              **
!     **   Last update: wedesday, 30/11/2011                         **
!     **                                                             **
!     *****************************************************************

      module chemkin_kiva

      use working_precision, only: dp
      use chemistry_string_limits, only: species_name_len

      implicit none
      public


!        chemkin problem variables
         real (dp) :: CKp, CKrho
!$OMP THREADPRIVATE(CKp, CKrho)

!        chemkinII file names
         character(len=*), parameter :: chembin = 'chem.bin', &
                                        chemdat = 'dat.chemkin'

!        Units for chemkinII link
         integer, parameter :: lin = 40, linc = 41, lout = 42

!        Working arrays lengths
         integer :: leniwk, lenrwk, lencwk

!        Chemkin Working arrays (MUST BE DOUBLE PRECISION!)
         integer,           dimension(:,:), allocatable, target :: intwork
         double precision,  dimension(:,:), allocatable, target :: reawork
         character(len=species_name_len), dimension(:,:), allocatable, target :: chawork

         integer,           dimension(:), allocatable :: ICK
         double precision,  dimension(:), allocatable :: RCK
         character(len=species_name_len), dimension(:), allocatable :: CCK

!        Chemkin problem arrays (elements, species, reactions, eqs)
         integer :: ckne, ckns, cknr, ckneq

!        Chemkin working arrays subparts (pointers) *******************

!        Species molar weights [g/mol]
!         real (dp)       , dimension(:), pointer :: CKmw
         real (dp)       , dimension(:), allocatable :: CKmw

!        VODE working arrays
         integer :: cklrw, ckliw
         DOUBLE PRECISION, dimension(:,:), allocatable :: ckrwork
         integer,          dimension(:,:), allocatable :: ckiwork

!        Integration monitoring arrays *********************************

!        Calls to the constant volume ODE and to the jacobian routine
         integer :: CKncJAC   = 0
         integer :: CKncCONV  = 0

!        Number of integration steps
         integer :: CKnsteps  = 0
!        Number of LU decompositions
         integer :: CKnLUdec  = 0
!        Number of Newton iterations
         integer :: CKnNewton = 0
!$OMP THREADPRIVATE(CKncJAC,CKncCONV,CKnsteps,CKnLUdec,CKnNewton)

      contains


!         *************************************************************
!         **              KIVA4 - CHEMKIN II interface               **
!         **     Initializaion of chemkin arrays from chem.bin       **
!         **                                                         **
!         **   Author:      Federico Perini                          **
!         **   Last update: tuesday, 22/11/2011                      **
!         **                                                         **
!         *************************************************************
          subroutine chemkin_initialize
          use chemkinII, only: cklen, ckinit, ckindx, NcWT
!ck2015
          use chemistry_setup, only: mechdir
!          use openmpvars, only: nthreads
!          use omp_lib
          implicit none

          integer :: nfit, i, nthreads
          logical :: present1, present2
          character(len=20) :: filename
          character(len=*), parameter :: &
             fmt_errfile = "(' missing file error in chemkinII'," &
                           //"' - KIVA init:',A20)",              &
             fmt_ckhead  = "(1x,' SpeedCHEM - chemkinII interface '/" &
                           //"'  ------------------------------- ')"

!            Opening units:
!            1) linc = input file generated by chemkin interp
             inquire(file = chembin,  exist = present1)
!ck2015             inquire(file = 'cklink', exist = present2)
             inquire(file = trim(mechdir)//"cklink", exist = present2)
             if ((.not.present1).and.(.not.present2)) then
                write(*,fmt_errfile)chembin
                stop
             endif

             if (present1) filename = trim(chembin)
!ck2015             if (present2) filename = 'cklink'
             if (present2) filename = trim(mechdir)//"cklink"

             open(unit=linc, file = filename, form   = 'unformatted',&
                                              status = 'unknown'     )

!            2) lout = output file (dat.chem, usually)
!ck2015             open(unit=lout, file = chemdat, action = 'write')
             open(unit=lout, file = trim(mechdir)//chemdat, action = 'write')

             write(lout, *         )
             write(lout, fmt_ckhead)
             write(lout, *         )

!            **********************************************************

!            Retrieve set number of threads
             nthreads = 1
!             if (omp_in_parallel()) then
!	            nthreads = omp_get_num_threads()
!	         else
!                nthreads = omp_get_num_procs()
!	         endif


!            Retrieve working arrays lengths and allocate working space
             call cklen (linc, lout, leniwk, lenrwk, lencwk)

             allocate ( intwork(leniwk,nthreads), ICK(leniwk), &
                        reawork(lenrwk,nthreads), RCK(lenrwk), &
                        chawork(lencwk,nthreads), CCK(lencwk)  )

             ICK = 0
             RCK = 0.d0
             CCK = ' '

!            Initialize chemistry problem and fill in working arrays
             call ckinit(leniwk,lenrwk,lencwk, linc, lout,    &
                         ICK, RCK, CCK)

!            Retrieve reaction mechanism dimensions
             call ckindx(ICK, RCK, ckne, ckns, cknr, nfit)

             do i = 1, nthreads
               intwork(1:leniwk,i) = ICK
               reawork(1:lenrwk,i) = RCK
               chawork(1:lencwk,i) = CCK
             end do

!            Initialize and allocate VODE integrator working arrays

             ckneq = ckns + 1
             ckliw = 30 + ckneq
             cklrw = 22 + 9*ckneq + 2 * ckneq**2

             allocate ( ckiwork(ckliw,nthreads),                   &
                        ckrwork(cklrw,nthreads) )

             ckiwork = 0
             ckrwork = 0.d0

!            Assign working array portions with variable names

!            Species molecular weights [g/mol]
!             if (associated(CKmw)) nullify(CKmw)
!                                   CKmw => RCKWRK(NcWT:NcWT+ckns-1, 1)

             allocate(CKmw(ckns))
             CKmw = reawork(NcWT:NcWT+ckns-1,1)

!            Closing files
             close(linc)
             close(lout)

             write(*,*)'end chemkin initialize'

          end subroutine chemkin_initialize


!         *************************************************************
!         **              KIVA4 - CHEMKIN II interface               **
!         **     Integration of the chemistry ODE system (conV)      **
!         **                                                         **
!         **   Author:      Federico Perini                          **
!         **   Last update: friday, 02/12/2011                       **
!         **                                                         **
!         *************************************************************
          subroutine ckII_integrate(neq,fun,rtol,atol,t0,tf,yin)

          use omp_lib

          implicit none

!   	  ** Problem - related variables
          integer,          intent(in)    :: neq
          real (dp)       , intent(in)    :: rtol
          real (dp)       , intent(in)    :: t0, tf
          real (dp)       , dimension(neq), intent(in)    :: atol
          real (dp)       , dimension(neq), intent(inout) :: yin
          real (dp)       , dimension(neq)                :: yin0
!          real (dp)       , dimension(:),   pointer       :: rwrk
!          integer,          dimension(:),   pointer       :: iwrk
          external           :: fun, dummy
          integer            :: istate, ipar, nthread, it
          real (dp)          :: rpar, rto2

          integer            :: itask, iopt, itol, mtype
          integer, parameter :: nit = 5

          character(len=*), parameter :: &
            fmt_errvode = "(' Error in CHEMKIN-II VODE integration ')"

!         T = yin(1)     = temperature [K]
!         Y = yin(2:neq) = species mass fractions [-]

!         ** Main integration *****************************************
          itask  = 1
          istate = 1
          iopt   = 1
          itol   = 2
          mtype  = 22
          it     = 0

          nthread = 1!omp_get_thread_num()+1

          ckiwork(1:ckliw, nthread) = 0
          ckrwork(1:cklrw, nthread) = 0.d0

          ckiwork(6, nthread) = 10000 ! Arbitrary number of steps

          intwork(1:leniwk, nthread) = ICK(1:leniwk)
          reawork(1:lenrwk, nthread) = RCK(1:lenrwk)

          yin0 = yin

!         Try to integrate up to nit times, at increasing relative
!         tolerance values
          integration_attempts: do while (it<nit .and. istate<2)

             yin    = yin0
             it     = it + 1
             istate = 1
             rto2   = rtol * 10.d0**(it-1)

             call DVODE(fun, neq, yin, t0, tf, itol, rto2, atol, itask,  &
                        istate, iopt, ckrwork(1:cklrw, nthread), cklrw,  &
                        ckiwork(1:ckliw,nthread), ckliw, dummy, mtype)

!            ** Integration statistics ***********************************
             CKncJAC   = CKncJAC   + ckiwork(13,nthread)
             CKncCONV  = CKncCONV  + ckiwork(12,nthread)
             CKnsteps  = CKnsteps  + ckiwork(11,nthread)
             CKnLUDEC  = CKnLUDEC  + ckiwork(19,nthread)
             CKnNewton = CKnNewton + ckiwork(20,nthread)

          end do integration_attempts

          if (it>1)write(*,*)'Completed in ',it,' iterations'


!         Check for errors
          if (istate < 2) then
             write(*,fmt_errvode)
             stop
          endif

          end subroutine ckII_integrate

!         *************************************************************
!         **              KIVA4 - CHEMKIN II interface               **
!         **     constant-volume reactor ODE system computation      **
!         **                                                         **
!         **   Author:      Federico Perini                          **
!         **   Last update: friday, 02/12/2011                       **
!         **                                                         **
!         *************************************************************
          subroutine ckII_conV(neq,time,yin,dyindt)

          use chemkinII, only: ckcvbs, ckwyr, ckums, ckuml, ckpy
          use kinetics_mod, only: SCdwdt
!          use omp_lib
          implicit none

          integer,          intent(in)          :: neq
          real (dp)       , intent(in)          :: time
          real (dp)       , intent(in),  target :: yin   (neq)
          real (dp)       , intent(out), target :: dyindt(neq)

          real (dp)                             :: cv
          real (dp)                             :: dwdt(neq-1), &
                                                   umol(neq-1)
          double precision                      :: dcv, dCKp
          double precision                      :: ddwdt(neq-1), &
                                                   dumol(neq-1)
          real (dp)       ,                 pointer  :: T, dTdt
          real (dp)       , dimension(:),   pointer  :: Y, dYdt

          integer                                    :: n, j
          integer,          dimension(:),   pointer  :: iw
          double precision, dimension(:),   pointer  :: rw


!         Pointer associations
          T => yin(1)
          Y => yin(2:neq)

          dTdt => dyindt(1)
          dYdt => dyindt(2:neq)

!         Thread-related stuff
          n       = 1!omp_get_thread_num() + 1

!         OpenMP Parallel-safe working array assignment
          iw  => intwork(1:leniwk,n)
          rw  => reawork(1:lenrwk,n)

!         NB: problem constant density CKrho [g/cm3]
!             has to be set prior to calling this subroutine!

!         Compute current system pressure
          call ckpy (dble(CKrho), dble(T), dble(Y), iw, rw, dCKp)
          CKp = real (dCKp, dp)

!         Compute mixture constant volume specific heat [erg/g/K]
          call ckcvbs(dble(T), dble(Y), iw, rw, dcv)
          cv = real( dcv, dp )

!         Compute species production rates dwdt [mol/cm3/s] at
!         constant volume
          call ckwyr (dble(CKrho), dble(T), dble(Y), iw, rw, ddwdt)
          dwdt = real(ddwdt, dp)

!         Compute internal energies of the species in molar units [erg/mol]
          call ckuml (dble(T), iw, rw, dumol)
          umol = real(dumol, dp)

!         Compute species mass fraction change for constant vol [1/s]
          dYdt = dwdt * CKmw / CKrho

!         Compute temperature change for constant volume reactor [K/s]
          dTdt = - sum( umol * dwdt ) / (CKrho * cv)

          end subroutine ckII_conV


      end module chemkin_kiva
