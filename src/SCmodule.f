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

c     *****************************************************************
c     **                                                             **
c     **                     SpeedCHEM FORTRAN                       **
c     **    A package for the integration and solution of chemical   **
c     **       kinetics problems                                     **
c     **                                                             **
c     **                        "Sempre fortissimo, e con strepito"  **
c     **                                                   F. Liszt  **
c     **                                                             **
c     **   Author:      (C) Federico Perini                          **
c     **   Last update: monday, 23/04/2012                           **
c     **                                                             **
c     *****************************************************************


c     ** GENERAL SETUP OF CHEMISTRY INTEGRATION ***********************
      module chemistry_setup

        use working_precision, only: dp
        implicit none
        public

c       Mechanism label
        character(len=80) :: mechanism

ck2015       Mechanism file path
        character(len=256) :: mechdir
c       Label for chemistry ODE solver
        character(len=15) :: solver

c       Label for current host
        character(len=30) :: host_name

c       Switch for choosing solution program (SpeedCHEM vs CHEMKIN-II)
c       SpeedCHEM is used by default
        logical :: use_speedchem = .true.

c       Switch for accurate vs. tabulated temperature-dependent
c       parameters
        logical :: accurate_scthermo

c       Switch for using analytical jacobian
        logical :: analytical_jac

c       Switch for simplified third-body molecularity treatment
c       Complete:  Meff = sum( alpha_i * C_i )
C       Sparse  :  Meff = C_tot + sum ( (1-alpha_i) * C_i )
        logical :: simplified_for_sparsity

c       Switch for saving temperature-dependent parameters between
c       evaluation of ODE function and of its jacobian
        logical :: save_thermal_parameters

c       Switch to check for the reaction mechanism consistency after
c       that the reaction data have been imported
        logical :: check_reaction_mechanism

c       Switch to choose whether to print out dat.jacobian with sparsity
c       pattern
        logical :: print_out_jacobian

c       Switch for using separate tolerances for temperature and species
        logical :: separate_tols

c       Switch for permutating species indices for minimum LU decompo-
c       sition sparse array storage
        logical :: permutate_species = .false.

c       Integration tolerances: relative, absolute for species,
c       absolute for temperature
        real (dp)        :: TOLR, YTOLA, TTOLA

c       Number of intermediate watchiing steps
        integer :: nwatch

c       Accuracy and bounds for temperature-dependent data tabulation
        real (dp)        :: Temp_table_accuracy
        real (dp)        :: rec_Ttable_accuracy
        real (dp)        :: Temp_HIlim
        real (dp)        :: Temp_LOlim
        integer          :: tab_nsteps


        contains


c     *****************************************************************
c     **                     SpeedCHEM FORTRAN                       **
c     **     Import solver and solution parameters from itapeChem    **
c     **                                                             **
c     **   Author:      Federico Perini                              **
c     **   Last update: tuesday, 22/11/2011                          **
c     **                                                             **
c     *****************************************************************
        subroutine solver_setup
        implicit none

        logical :: present
        integer :: idummy
        character(len=15) :: tmpsolver, programme
        character(len=*), parameter :: itapeTSS = "itapeTSS",
     &                                 itapeCH  = "itapeChem"

c       Output labels
        character(len=*), parameter ::
     &    fmt_nofile = "(' itapeChem not found, assuming default " //
     &                     "values ')",
     &    fmt_bar = "(' --------------------------'"//
     &              "'----------------------------')",
     &    fmt_prog   = "(' chemistry program     : ',A15)",
     &    fmt_solver = "(' integration solver    : ',A15)",
     &    fmt_jac    = "(' jacobian formulation  : ',A10)",
     &    fmt_thermo = "(' thermal data accuracy : ',A9 )",
     &    fmt_save   = "(' thermal data behavior : ',A9 )",
     &    fmt_spars  = "(' third-body formulation: ',A8 )",
     &    fmt_tols   = "(' integration tolerances  ')",
     &    fmt_rtol   = "('    relative           : ',1P,E8.2)",
     &    fmt_atol   = "('    absolute           : ',1P,E8.2)",
     &    fmt_ytol   = "('    absolute (species) : ',1P,E8.2)",
     &    fmt_ttol   = "('    absolute (temp)    : ',1P,E8.2)",
     &    fmt_erprog = "(' ERROR - program not recognized: ',A15)",
     &    fmt_ttab   = "(' temperature tabulation  ')",
     &    fmt_tstp   = "('    sampling step [K]  : ',1PE9.3)",
     &    fmt_lowt   = "('    low   limit   [K]  : ',1PE9.3)",
     &    fmt_hit    = "('    high  limit   [K]  : ',1PE9.3)",
     &    fmt_host   = "(' running host          : ',A15)"

c       Import multi timescale integration settings
        inquire(file=itapeCH,exist=present)

        if (present) then

         open(unit = 110,file=itapeCH)
         read(110,*)
         read(110,*) programme
         read(110,*) tmpsolver
         read(110,*) idummy

         if (idummy == 1) analytical_jac = .true.
         if (idummy == 0) analytical_jac = .false.

c        Convert string to uppercase
         call s_cap(programme)
         call s_cap(tmpsolver)
c        Compute solver string including Jacobian
         solver = trim(adjustl(tmpsolver))
         if (analytical_jac) solver = trim(solver)//'JAC'

         if (trim(adjustl(programme)) == 'CHEMKIN') then
            use_speedchem = .false.
         elseif (trim(adjustl(programme)) == 'SPEEDCHEM') then
            use_speedchem = .true.
         else
            write(*,fmt_erprog)trim(adjustl(programme))
            stop
         endif

         read(110,*) idummy

         if (idummy == 1) accurate_scthermo = .false.
         if (idummy == 0) accurate_scthermo = .true.

         read(110,*) idummy

         if (idummy == 1) save_thermal_parameters = .true.
         if (idummy == 0) save_thermal_parameters = .false.

         read(110,*) idummy

         if (idummy == 0) simplified_for_sparsity = .false.
         if (idummy == 1) simplified_for_sparsity = .true.

         if (.not.analytical_jac) simplified_for_sparsity = .false.

         read(110,*) idummy

         if (idummy == 0) print_out_jacobian = .false.
         if (idummy == 1) print_out_jacobian = .true.

         read(110,*) idummy

         if (idummy == 0) check_reaction_mechanism = .false.
         if (idummy == 1) check_reaction_mechanism = .true.

         read(110,*) idummy

         if (idummy == 0) separate_tols = .false.
         if (idummy == 1) separate_tols = .true.

         read(110,*) TOLR
         read(110,*) YTOLA
         read(110,*) TTOLA

c        Reading tabulation parameters
         read(110,*) Temp_table_accuracy
         read(110,*) Temp_LOlim
         read(110,*) Temp_HIlim

c        Number of intermediate output points per simulation
         read(110,*) nwatch



         close(110)

        else

c          Apply default setup if 'itapeChem' is not present
           programme                = "SPEEDCHEM"
           tmpsolver                = "LSODES"
           solver                   = "LSODESJAC"
           analytical_jac           = .true.
           accurate_scthermo        = .false.
           save_thermal_parameters  = .false.
           separate_tols            = .false.
           simplified_for_sparsity  = .true.
ck2015           check_reaction_mechanism = .true.
           check_reaction_mechanism = .false.
           print_out_jacobian       = .false.
           TOLR                     = 1.0e-04_dp
           YTOLA                    = 1.0e-15_dp
           TTOLA                    = 1.0e-02_dp
           nwatch                   = 2
           Temp_table_accuracy      = 1.e1_dp
           Temp_LOlim               = 300.e0_dp
           Temp_HIlim               = 3500.e0_dp

           write(*,fmt_nofile)
        endif

c       ** Output solver setup data to screen **************************

        write(*,fmt_prog  )programme

        write(*,fmt_solver)tmpsolver

        if (analytical_jac) then
           write(*,fmt_jac)"ANALYTICAL"
        else
           write(*,fmt_jac)"NUMERICAL "
        endif

        if (accurate_scthermo) then
           write(*,fmt_thermo)"ALGEBRAIC"
        else
           write(*,fmt_thermo)"TABULATED"
           write(*,fmt_ttab  )
           write(*,fmt_tstp  )Temp_table_accuracy
           write(*,fmt_lowt  )Temp_LOlim
           write(*,fmt_hit   )Temp_HIlim
        endif

        if (save_thermal_parameters) then
           write(*,fmt_save  )"SAVE     "
        else
           write(*,fmt_save  )"CALCULATE"
        endif

c       No need to save data if there is no Jacobian
c        if (.not.analytical_jac) save_thermal_parameters = .false.
c       [fp] function deactivated
c        save_thermal_parameters = .false.

        if (simplified_for_sparsity) then
           write(*,fmt_spars)"SPARSE  "
        else
           write(*,fmt_spars)"COMPLETE"
        endif

        write(*,fmt_tols)
        if (separate_tols) then
           write(*,fmt_rtol)TOLR
           write(*,fmt_ytol)YTOLA
           write(*,fmt_ttol)TTOLA
        else
           write(*,fmt_rtol)TOLR
           write(*,fmt_atol)YTOLA
        endif

        write(*,*)

        ! Output of running host machine
        call hostnm(host_name)
        write(*,fmt_host)adjustl(host_name)

        write(*,*)

        end subroutine solver_setup




      end module chemistry_setup


c     *****************************************************************
c     **                                                             **
c     **                     SpeedCHEM FORTRAN                       **
c     **                                                             **
c     **         Module for sparse chemistry computations            **
c     **                                                             **
c     **                                                             **
c     **   Author:      Federico Perini                              **
c     **   Last update: wednesday, 23/11/2011                        **
c     **                                                             **
c     *****************************************************************
      module sparse_chemistry

         use working_precision, only: dp
         use sparse_definitions
         use sparse_algebra

         implicit none
         public

!        Sparse matrix storage for:
!        --------------------------

!        nudiff = stoich_products - stoich_reactants
         type(sparse)    :: nudiff_sparse
         type(sparse)    :: nudiffT_sparse
         type(sparse)    :: nudiff_EQREV_sp
         type(sparse)    :: stoich_r_sp, stoich_r_eff_sp
         type(sparse)    :: stoich_p_sp, stoich_p_eff_sp

         type(sparseint) :: inudiff_sparse
         type(sparseint) :: inudiffT_sparse
         type(sparseint) :: istoich_r_sp, istoich_r_eff_sp
         type(sparseint) :: istoich_p_sp, istoich_p_eff_sp

!        Element matrix
         type(sparse)    :: EM_sp

!        Jacobian sparsity pattern and jacobian values
!         type(sparse)            :: ijac_sparse
         type(sparse), target    :: JAC_sparse
         type(sparse)            :: JACT_sparse
!$OMP THREADPRIVATE(JAC_SPARSE,JACT_sparse)

!        Species vs species subpart of the jacobian matrix
!        JACYY(:,:) = JAC(2:neq, 2:neq); JACYYT = transpose(JACYY)
         type(sparse)    :: JACYY_sparse, JACYYT_sparse
!$OMP THREADPRIVATE(JACYY_sparse,JACYYT_sparse)

!        Premultiplied nudiff * molar volume = nudiff * MW / rho
         type(sparse)    :: nudiffT_molarv_sparse
!$OMP THREADPRIVATE(nudiffT_molarv_sparse)

!        Derivative of reaction progress variable with respect to
!        species mass fractions dq_dY [nr x ns]
         type(sparse)    :: dq_dY_sparse, dq_dY_T_sparse
!$OMP THREADPRIVATE(dq_dY_sparse, dq_dY_T_sparse)

!        Third-body enhanced molecularity coefficients [ntbALL x ns]
         type(sparse)    :: third_body_sp
         type(sparse)    :: tb_beta_sp, tbb_uMW_sp



         contains

           !   *********************************************************
           !   **           Chemistry sparse matrix setup             **
           !   **                                                     **
           !   **   Author:      Federico Perini                      **
           !   **   Last update: wednesday, 23/11/2011                **
           !   *********************************************************

           subroutine sparse_chemistry_setup(inotrev,nnotrev,Ainf)

              implicit none

              integer, dimension(:), intent(in) :: inotrev
              integer, intent(in) :: nnotrev
              real (dp)       , dimension(:), intent(in) :: Ainf

              integer :: i

              ! Convert reaction mechanism matrices into sparse format
              ! Convert stoichiometric coefficients into sparse format
              call sparse_compress(nudiff_sparse)
              nudiffT_sparse = sparse_transpose(nudiff_sparse)

              ! Convert stoichiometric coefficients into sparse int
              istoich_r_sp    = stoich_r_sp
              istoich_p_sp    = stoich_p_sp
              inudiff_sparse  = nudiff_sparse
              inudiffT_sparse = nudiffT_sparse

              ! Subset for effective reversible reactions only
              stoich_p_eff_sp = stoich_p_sp
              if (nnotrev < stoich_p_sp%nr) then
                 do i = 1, nnotrev
                    call remove_line(stoich_p_eff_sp,inotrev(i))
                 end do
              endif

              ! Subset for effecive reactions with non-zero forward
              ! reaction rate constant
              stoich_r_eff_sp = stoich_r_sp
              do i = 1, stoich_r_sp%nr
               if (Ainf(i)==0.e0_dp)call remove_line(stoich_r_eff_sp,i)
              end do

              istoich_r_eff_sp = stoich_r_eff_sp
              istoich_p_eff_sp = stoich_p_eff_sp

           end subroutine sparse_chemistry_setup


c        **************************************************************
c        **  Permutate mechanism definitions in sparse matrix form   **
c        **************************************************************

         subroutine permutate_sparse_matrices(new_order)
         implicit none

         integer, dimension(:), intent(in) :: new_order


c        Permutate matrices related to stoichiometric coefficients
         nudiff_sparse   = sparse_col_permutation(nudiff_sparse,
     &                                            new_order)
         nudiffT_sparse  = sparse_row_permutation(nudiffT_sparse,
     &                                            new_order)
         if (nudiff_EQREV_sp%n > 0)
     &   nudiff_EQREV_sp = sparse_col_permutation(nudiff_EQREV_sp,
     &                                            new_order)
         stoich_r_sp     = sparse_col_permutation(stoich_r_sp,
     &                                            new_order)
         stoich_r_eff_sp = sparse_col_permutation(stoich_r_eff_sp,
     &                                            new_order)
         stoich_p_sp     = sparse_col_permutation(stoich_p_sp,
     &                                            new_order)
         stoich_p_eff_sp = sparse_col_permutation(stoich_p_eff_sp,
     &                                            new_order)
         inudiff_sparse  = sparseint_col_permutation(inudiff_sparse,
     &                                               new_order)
         inudiffT_sparse = sparseint_row_permutation(inudiffT_sparse,
     &                                               new_order)
         istoich_r_sp    = sparseint_col_permutation(istoich_r_sp,
     &                                               new_order)
         istoich_p_sp    = sparseint_col_permutation(istoich_p_sp,
     &                                               new_order)

!        Element matrix
         EM_sp      = sparse_row_permutation(EM_sp, new_order)

!        Permutate third-body enhanced molecularity coefficients
         tb_beta_sp = sparse_col_permutation(tb_beta_sp,new_order)
         tbb_uMW_sp = sparse_col_permutation(tbb_uMW_sp,new_order)

         end subroutine permutate_sparse_matrices



      end module sparse_chemistry





c     *****************************************************************
c     **                                                             **
c     **                     SpeedCHEM MAIN MODULE                   **
c     **          Contains reaction mechanism information            **
c     **                                                             **
c     **                                                             **
c     **   Author:      Federico Perini                              **
c     **   Last update: monday, 23/04/2012                           **
c     **                                                             **
c     *****************************************************************

      module speedchem

      use working_precision, only: dp
      implicit none
      public

c     ** number of elements, species and reactions
      integer :: nel,ns,nr


      integer, dimension(:), allocatable :: species, reactions
      integer, dimension(:), allocatable :: species_permutations,
     &                                      species_inverse_permutations

c     ** number of chemistry problem equations
      integer          :: neq
      real (dp)        :: uneq

c     ** arrhenius parameters
      real (dp)       , dimension(:), allocatable :: A0,AREV,AINF,ARi,
     &                                               b0,bREV,bINF,bRi,
     &                                               E0,EREV,EINF,ERi

c     ** species molecular weights and reciprocal
      real (dp)       , dimension(:), allocatable :: SCMW, uMW

c     ** elements molecular weights and reciprocal
      real (dp)       , dimension(:), allocatable :: elMW, uelMW

c     ** reaction reversibility and declared reverse data
      logical, dimension(:), allocatable :: reversibile, REV
      integer, dimension(:), allocatable :: inotrev
      integer                            :: nnotrev


c     ** declared duplicate reactions
      logical, dimension(:), allocatable :: duplicate


c     ** reaction coefficients
      real (dp)       , dimension(:,:), allocatable :: stoich_r,
     &                                                 stoich_p,
     &                                                 nudiff,
     &                                                 nudiffT,
     &                                                 third_body,
     &                                                 third_body_beta

      real (dp)       , dimension(:,:), allocatable :: stoich_r_pack,
     &                                                 stoich_p_pack,
     &                                                 nudiff_pack

      integer,          dimension(:,:), allocatable :: i_stoich_r,
     &                                                 i_stoich_p,
     &                                                 i_nudiff

      integer,          dimension(:),   allocatable :: n_stoich_r,
     &                                                 n_stoich_p,
     &                                                 n_nudiff

      real (dp)       , dimension(:), allocatable :: sumnudiff
      integer,          dimension(:), allocatable :: isumnudiff

c     ** pressure-dependent reactions
      logical, dimension(:), allocatable :: HIGH, LOW, TROE, THREE
      real (dp)       , dimension(:), allocatable ::   aTROE, T1TROE,
     &                                                 T2TROE, T3TROE

c     ** Three-body reaction indices
      integer, dimension(:), allocatable :: iTHREE ! All third body
      integer, dimension(:), allocatable :: iTB ! simple third body only

      integer                            :: nTHREE, nTB

c     ** character strings for elements and species
      character(len=18), dimension(:), allocatable :: specie
      character(len=2) , dimension(:), allocatable :: elementi

c     ** Storage of chemistry jacobian sparsity ************************
c     ijac contains sparsity information,
c     rowjac the row index, coljac the column index of the nonzero
c     elements
      logical                              :: sparse_jac = .false.
      integer                              :: njac
c      integer, dimension(:,:), allocatable :: ijac, ijacYY
      logical, dimension(:,:), allocatable :: ljac
      integer, dimension(:),   allocatable :: rowjac, coljac
!$OMP THREADPRIVATE(sparse_jac,njac, rowjac, coljac, ljac)

c     ** Index of the oxidizer within the species array
      integer :: iO2




c     ******************************************************************
      contains

c     Exporting mechanism information to file
      subroutine SCmechanism_to_file
c     *****************************************************************
c     **                                                             **
c     **                     SpeedCHEM FORTRAN                       **
c     **                                                             **
c     **     Export SpeedCHEM Chemistry Mechanism to file            **
c     **                                                             **
c     **                                                             **
c     **   Author:      Federico Perini                              **
c     **   Last update: sunday, 23/10/2011                           **
c     **                                                             **
c     *****************************************************************

      use sparse_chemistry, only:   stoich_r_sp, stoich_p_sp,
     &                              third_body_sp, tb_beta_sp
      use sparse_algebra,   only:   sparse_value, sparse_internal_count
ck2015 mechdir 1 line
      use chemistry_setup,  only:   mechdir


      implicit none
      integer                     :: i,j,countns
      integer                     :: count_tb(nr), count_tb_beta(nTHREE)
      real (dp)                   :: current_value
      character(len=255)          :: reaction_line, tmp, tmpval

c     String formats
      character(len=*), parameter ::
     &
     & fmt_head = "(' Reaction mechanism output information ' /
     &              ' ------------------------------------- '/)",
     & fmt_reac = "(' Reaction information [cal,mol,K]      ', 20x,
     &                '  A',15x,' b ',8x,' E  ' /
     &              ' ------------------------------------- ', 20x,
     &                (46('-'))/)",
     & fmt_el  = "(' Element information                   ' /
     &              ' ------------------------------------- '/)",
     & fmt_sp  = "(' Species information                   ' /
     &              ' ------------------------------------- '/)",
     & fmt_high = "(7x,' HIGH / ',1pE13.6,' , ',E10.3,' , ',
     &                            1pE13.6,' / ')",
     & fmt_low  = "(7x,' LOW  / ',1pE13.6,' , ',E10.3,' , ',
     &                            1pE13.6,' / ')",
     & fmt_rev  = "(52x,' REV   ',1pE13.6,' ,  ',F7.2 ,1x,' ,  ',
     &                            1pE13.6,'   ')",
     & fmt_troe = "(7x,' TROE / ',(3(1pE13.6,' , ')),1pE13.6,' / ')",
     & fmt_nel  = "(' Number of elements : ',i5)",
     & fmt_ns   = "(' Number of species  : ',i5)",
     & fmt_nr   = "(' Number of reactions: ',i5)"

c     ******************************************************************

ck2015      open(unit = 13, file = 'dat.mechanism')
      open(unit = 13, file = trim(mechdir)//"dat.mechanism")
      write(13,fmt_head)

      write(13,fmt_nel )nel
      write(13,fmt_ns  )ns
      write(13,fmt_nr  )nr

      write(13,*)


c     Storing enhanced three body reactions coefficients
      call sparse_internal_count(third_body_sp, count_tb,      dim=2)
      call sparse_internal_count(tb_beta_sp,    count_tb_beta, dim=2)

c     Writing reaction information
      write(13,fmt_reac)
      do i = 1, nr

c         Reaction number
          write(reaction_line,"(I5)")i

          reaction_line = trim(reaction_line)//'.  '

c         Add reactant information to reaction line string
          countns = 0
          do j = 1, ns
               current_value = sparse_value(stoich_r_sp,i,j)
               if (current_value/=0.e0_dp) then
                  countns = countns + 1
                  tmp = ''
ck2015s for real stoichometry coefficients
c                  if(current_value>1.e0_dp)
c    &                write(tmp,"(I2)") int(current_value+tiny(0.e0_dp))
                    if(current_value-int(current_value) > 0.e0_dp) then
                     write(tmp,"(f3.1)") current_value
                    else
                     if(current_value>1.e0_dp)
     &               write(tmp,"(I2)") int(current_value+tiny(0.e0_dp))
                    endif
ck2015e

                  if (countns>1) tmp = ' + '//trim(tmp)

                     reaction_line = trim(reaction_line) // trim(tmp) //
     &                               ' ' // specie(j)

               endif
          end do

          if (count_tb(i)>0) then
             if (HIGH(i).or.LOW(i)) then
                reaction_line = trim(reaction_line) // ' (+M)'
             else
                reaction_line = trim(reaction_line) // ' +M '
             endif
          endif

c         Add reversible / irreversible reaction equals sign
          if (reversibile(i)) then
             reaction_line = trim(reaction_line) // ' <=> '
          else
             reaction_line = trim(reaction_line) // ' => '
          endif

c         Add reactant information to reaction line string
          countns = 0
          do j = 1, ns
               current_value = sparse_value(stoich_p_sp,i,j)
               if (current_value/=0.e0_dp) then
                  countns = countns + 1
                  tmp = ''
                  if(current_value>1.e0_dp)
     &                write(tmp,"(I2)") int(current_value+tiny(0.e0_dp))

                  if (countns>1) tmp = ' + '//trim(tmp)

                     reaction_line = trim(reaction_line) // trim(tmp) //
     &                               ' ' // specie(j)

               endif
          end do

          if (count_tb(i)>0) then
             if (HIGH(i).or.LOW(i)) then
                reaction_line = trim(reaction_line) // ' (+M)'
             else
                reaction_line = trim(reaction_line) // ' +M'
             endif
          endif

c         Adding Arrhenius parameters to reaction string
          if (HIGH(i)) then
             write(reaction_line(59:71),"(1pE13.6)")A0(i)
             write(reaction_line(76:82),"(F7.2   )")b0(i)
             write(reaction_line(87:100),"(1pE14.6)")E0(i)
          else
             write(reaction_line(59:71),"(1pE13.6)")Ainf(i)
             write(reaction_line(76:82),"(F7.2   )")binf(i)
             write(reaction_line(87:100),"(1pE14.6)")Einf(i)
          endif

c         Write reaction line to file
          write(13,"(A101)") trim(reaction_line)

c         Check if reaction contains additional data (e.g. HIGH or LOW
c         Arrhenius parameter specifications, or third-body effective
c         molecularities
          if (HIGH(i)) write(13,fmt_high)Ainf(i),binf(i),Einf(i)
          if (LOW(i))  write(13,fmt_low )A0  (i),b0  (i),E0  (i)
c         NB: Troe is ordered as can be read by Chemkin!!
          if (TROE(i)) write(13,fmt_troe)aTROE(i),T3TROE(i),
     &                                  T1TROE(i),T2TROE(i)

          if (REV(i))  write(13,fmt_rev )AREV(i),bREV(i),EREV(i)

          if (count_tb(i)>0) then

             tmp = ''
             tmpval = ''
             countns = 0
             do j = 1, ns

                if (sparse_value(third_body_sp,i,j)/=1.e0_dp) then
                    countns = countns + 1
                    write(tmpval,"(F7.2)")
     &              sparse_value(third_body_sp,i,j)
                    tmp = trim(tmp)//'  '//trim(specie(j)) // ' / ' //
     &                    trim(tmpval) // ' / '

                endif

                if (countns >= 5) then
                   countns = 0
c                   write(13,*)tmp
C                  DA SISTEMARE!
!                   if (count_tb(i)/=count_tb_beta(i))
                   write(13,*)'       ',trim(tmp)
                   tmp = ''
                endif
             end do
             write(13,*)'       ',trim(tmp)
          endif ! Reaction is pressure dependent

      end do ! Reaction information cycle

      write(13,*)

c     Write element information
      write(13,fmt_el)
      tmp = elementi(1)
      do i = 2, nel
         tmp = trim(tmp) // ' ' // elementi(i)
      end do
      write(13,*)trim(tmp)

      write(13,*)
      write(13,*)

c     Write species information
      write(13,fmt_sp)
      tmp = ' ' // specie(1)
      countns = 1
      do i = 2, ns
         tmp = trim(tmp) // ' ' // specie(i)
         countns = countns + 1
         if (countns >= 5) then
            countns = 0
            write(13,*)trim(tmp)
            tmp = ''
         endif

      end do
      if (countns>0)write(13,*)trim(tmp)

      write(13,*)

      close(13)

      end subroutine SCmechanism_to_file


      subroutine SCmechanism_to_chemkin
c     *****************************************************************
c     **     Export SpeedCHEM Chemistry in CHEMKIN ascii format      **
c     *****************************************************************
ck2015      use chemistry_setup,  only: mechanism
      use chemistry_setup,  only: mechanism, mechdir
      use sparse_chemistry, only: stoich_p_sp, stoich_r_sp,
     &                            third_body_sp
      use sparse_algebra,   only: sparse_value, sparse_internal_count

      implicit none
      logical                :: orderfile_present
      integer                :: i,j,lr,countns
      real (dp)              :: current_value
      integer, dimension(ns) :: tmpip
      integer, dimension(nr) :: count_tb

      character(len=255)          :: reaction_line, tmp, tmpval

c     String formats
      character(len=*), parameter ::
     &
     & fmt_mechhd = "('!',A)",
     & fmt_elems  = "('ELEMENTS')",
     & fmt_specs  = "('SPECIES')",
     & fmt_reacs  = "('REACTIONS')",
     & fmt_end    = "('END')",
     & fmt_reac   = "(' Reaction information [cal,mol,K]      ', 20x,
     &                '  A',15x,' b ',8x,' E  ' /
     &                ' ------------------------------------- ', 20x,
     &                (46('-'))/)",
     & fmt_el     = "(' Element information                   ' /
     &                ' ------------------------------------- ' /)",
     & fmt_sp     = "(' Species information                   ' /
     &              ' ------------------------------------- '/)",
     & fmt_high   = "(7x,' HIGH / ',1pE13.6,' ',E10.3,'  ',
     &                            1pE13.6,' / ')",
     & fmt_low    = "(7x,' LOW  / ',1pE13.6,' ',E10.3,'  ',
     &                            1pE13.6,' / ')",
     & fmt_rev    = "(7x,' REV  / ',1pE13.6,' ',E10.3,'  ',
     &                            1pE13.6,' / ')",
     & fmt_troe   = "(7x,' TROE / ',(3(1pE13.6,' '))' / ')",
     & fmt_troe4  = "(7x,' TROE / ',(3(1pE13.6,' ')),1pE13.6,' / ')",
     & fmt_nel    = "(' Number of elements : ',i5)",
     & fmt_ns     = "(' Number of species  : ',i5)",
     & fmt_nr     = "(' Number of reactions: ',i5)",
     & fmt_dup    = "(7X,' DUPLICATE')"

c     ******************************************************************

ck2015      open(unit = 13, file = 'dat.cheminp')
      open(unit = 13, file = trim(mechdir)//"dat.cheminp")
      write(13,fmt_mechhd)mechanism

c     Element data
      write(13,fmt_elems)
         tmp = elementi(1)
         do i = 2, nel
            tmp = trim(tmp) // ' ' // elementi(i)
         end do
         write(13,*)trim(tmp)
      write(13,fmt_end)


c     Species data
      write(13,fmt_specs)
         tmp = ' ' // specie(1)
         countns = 1
         do i = 2, ns
            tmp = trim(tmp) // ' ' // specie(i)
            countns = countns + 1
            if (countns >= 5 .or. len(trim(tmp))>=70) then
                countns = 0
               write(13,*)trim(tmp)
               tmp = ''
            endif
         end do
         if (countns>0)write(13,*)trim(tmp)
      write(13,fmt_end)

c     Three body reactions
      call sparse_internal_count(third_body_sp, count_tb, dim = 2)

c     Writing reaction information
      write(13,fmt_reacs)
      reaction_loop: do i = 1, nr

          reaction_line = ' '

c         Add reactant information to reaction line string
          countns = 0
          do j = 1, ns
               current_value = sparse_value(stoich_r_sp,i,j)
               if (current_value/=0.e0_dp) then
                  countns = countns + 1
                  tmp = ''
ck2015s for real stoichometric coefficients
c                 if(current_value>1.e0_dp)
c    &                write(tmp,"(I2)") int(current_value+tiny(0.e0_dp))
                    if(current_value-int(current_value) > 0.e0_dp) then
                     write(tmp,"(f3.1)") current_value
                    else
                     if(current_value>1.e0_dp)
     &               write(tmp,"(I2)") int(current_value+tiny(0.e0_dp))
                    endif
ck2015e
                  if (countns>1) tmp = '+'//trim(tmp)

                     reaction_line = trim(reaction_line) // trim(tmp)
     &                               // specie(j)

               endif
          end do

          if (count_tb(i)>0) then
             if (HIGH(i).or.LOW(i)) then
                reaction_line = trim(reaction_line) // '(+M)'
             else
                reaction_line = trim(reaction_line) // '+M'
             endif
          endif

c         Add reversible / irreversible reaction equals sign
          if (reversibile(i)) then
             reaction_line = trim(reaction_line) // '<=>'
          else
             reaction_line = trim(reaction_line) // '=>'
          endif

c         Add reactant information to reaction line string
          countns = 0
          do j = 1, ns
               current_value = sparse_value(stoich_p_sp,i,j)
               if (current_value/=0.e0_dp) then
                  countns = countns + 1
                  tmp = ''
                  if(current_value>1.e0_dp)
     &                write(tmp,"(I2)") int(current_value+tiny(0.e0_dp))

                  if (countns>1) tmp = '+'//trim(tmp)

                     reaction_line = trim(reaction_line) // trim(tmp)
     &                               // specie(j)

               endif
          end do

          if (count_tb(i)>0) then
             if (HIGH(i).or.LOW(i)) then
                reaction_line = trim(reaction_line) // '(+M)'
             else
                reaction_line = trim(reaction_line) // '+M'
             endif
          endif

c         Reaction line length
          lr = len(trim(adjustl(reaction_line)))

          if (lr<38) then
c         Adding Arrhenius parameters to reaction string
          if (HIGH(i)) then
             write(reaction_line(39:51),"(1pE13.6)")A0(i)
             write(reaction_line(56:62),"(F7.2   )")b0(i)
             write(reaction_line(67:80),"(1pE14.6)")E0(i)
          else
             write(reaction_line(39:51),"(1pE13.6)")Ainf(i)
             write(reaction_line(56:62),"(F7.2   )")binf(i)
             write(reaction_line(67:80),"(1pE14.6)")Einf(i)
          endif

c         Write reaction line to file
          write(13,"(A80)") trim(reaction_line)

          else ! (if lr < 38)

c         Adding Arrhenius parameters to reaction string
          if (HIGH(i)) then
             write(reaction_line(lr+2:lr+14),"(1pE13.6)")A0(i)
             write(reaction_line(lr+16:lr+22),"(F7.2   )")b0(i)
             write(reaction_line(lr+24:lr+37),"(1pE14.6)")E0(i)
          else
             write(reaction_line(lr+2:lr+14),"(1pE13.6)")Ainf(i)
             write(reaction_line(lr+16:lr+22),"(F7.2   )")binf(i)
             write(reaction_line(lr+24:lr+37),"(1pE14.6)")Einf(i)
          endif

c         Write reaction line to file

          write(13,"(2(A80))") trim(reaction_line(1:80)),
     &                         trim(reaction_line(81:lr+37))

          endif ! (lr < 38)

c         Check if reaction contains additional data (e.g. HIGH or LOW
c         Arrhenius parameter specifications, or third-body effective
c         molecularities
          if (duplicate(i)) write(13,fmt_dup)
          if (HIGH(i))      write(13,fmt_high)Ainf(i),binf(i),Einf(i)
          if (LOW(i))       write(13,fmt_low )A0  (i),b0  (i),E0  (i)

c         NB: Troe is ordered as can be read by Chemkin!!
          if (TROE(i) .and. T2TROE(i)==0.e0_dp)
     &        write(13,fmt_troe )aTROE(i),T3TROE(i),T1TROE(i)
          if (TROE(i) .and. T2TROE(i)/=0.e0_dp)
     &        write(13,fmt_troe4)aTROE(i),T3TROE(i),T1TROE(i),T2TROE(i)

          if (REV(i))  write(13,fmt_rev )AREV(i),bREV(i),EREV(i)

          if (count_tb(i)>0) then

c            Insert blank space
c             if (count(third_body(i,:)==1.e0_dp)/=ns) write(13,*)

             tmp = ''
             tmpval = ''
             countns = 0
             do j = 1, ns

                if (sparse_value(third_body_sp,i,j)/=1.e0_dp) then
                    countns = countns + 1
                    write(tmpval,"(F7.2)")
     &              sparse_value(third_body_sp,i,j)
                    tmp = trim(tmp)//trim(specie(j)) // ' /' //
     &                    trim(tmpval) // '/ '

                endif

                if (countns >= 5) then
                   countns = 0
c                   write(13,*)tmp
c                   if (count(third_body(i,:)==1.e0_dp)/=ns)
                   write(13,*)'       ',trim(tmp)
                   tmp = ''
                endif
             end do
             write(13,*)'       ',trim(tmp)
          endif ! Reaction is pressure dependent

      end do reaction_loop! Reaction information cycle

      if (nr>0) write(13,fmt_end)

      write(13,*)

      close(13)

      end subroutine SCmechanism_to_chemkin


c        **************************************************************
c        **     Initialise stoichiometric species  indices           **
c        **************************************************************
         subroutine init_stoich_indices

          implicit none

          integer :: i, ir, is

          allocate(n_stoich_r(ns), n_stoich_p(ns), n_nudiff(ns))

          do is = 1, ns
             n_stoich_r(is) = count(stoich_r(:,is)/=0.e0_dp)
             n_stoich_p(is) = count(stoich_p(:,is)/=0.e0_dp)
             n_nudiff  (is) = count(nudiff  (:,is)/=0.e0_dp)
          end do

          allocate(i_stoich_r(maxval(n_stoich_r,1), ns))
          allocate(i_stoich_p(maxval(n_stoich_p,1), ns))
          allocate(i_nudiff  (maxval(n_nudiff  ,1), ns))
          allocate(stoich_r_pack(maxval(n_stoich_r,1), ns))
          allocate(stoich_p_pack(maxval(n_stoich_p,1), ns))
          allocate(nudiff_pack  (maxval(n_nudiff  ,1), ns))

          do is = 1, ns

            i_stoich_r(1:n_stoich_r(is),is) =
     &      pack(reactions,stoich_r(1:nr,is)/=0.e0_dp)

            i_stoich_p(1:n_stoich_p(is),is) =
     &      pack(reactions,stoich_p(1:nr,is)/=0.e0_dp)

            i_nudiff  (1:n_nudiff  (is),is) =
     &      pack(reactions,nudiff  (1:nr,is)/=0.e0_dp)

            stoich_r_pack(1:n_stoich_r(is),is) =
     &      pack(stoich_r(1:nr,is),stoich_r(1:nr,is)/=0.e0_dp)

            stoich_p_pack(1:n_stoich_p(is),is) =
     &      pack(stoich_p(1:nr,is),stoich_p(1:nr,is)/=0.e0_dp)

            nudiff_pack  (1:n_nudiff  (is),is) =
     &      pack(nudiff  (1:nr,is),nudiff  (1:nr,is)/=0.e0_dp)

          end do


         end subroutine init_stoich_indices





c        ***************************************************************
c        **                                                           **
c        **     Computing mole fractions based on mass fractions      **
c        **                                                           **
c        ***************************************************************

         function massfr_to_molefr(Y) result(X)
         implicit none

         real (dp)       , dimension(ns), intent(in)  :: Y
         real (dp)       , dimension(ns)              :: X
         real (dp)        :: uavgMW

c        **************************************************************

c        Computing reciprocal of average mixture molecular weight [g/mol]
         uavgMW = 1.e0_dp/sum(Y*uMW)

c        Computing mole fractions
         X = uMW * Y * uavgMW

         end function massfr_to_molefr


c        ***************************************************************
c        **                                                           **
c        **     Computing mass fractions based on mole fractions      **
c        **                                                           **
c        ***************************************************************
         subroutine molefr_to_massfr(X,Y)
         implicit none

         real (dp)       , dimension(size(SCMW)), intent(in)  :: X
         real (dp)       , dimension(size(SCMW)), intent(out) :: Y
         real (dp)        :: avgMW

c        ***************************************************************

c        Computing average mixture molecular weight [g/mol]
         avgMW = sum(X*SCMW)

c        Computing mass fractions
         Y = X * SCMW / avgMW

         end subroutine molefr_to_massfr

c        **************************************************************
c          Permutates all the arrays depending on species into a new
c          species ordering previously defined
c        **************************************************************
         subroutine permutate_mechanism_species(new_order)
         implicit none

         integer, dimension(:), intent(in) :: new_order
         integer                           :: i, j


c        Permutate species molecular weights
         SCMW = SCMW(new_order)
         uMW  = uMW (new_order)

c        Permutate stoichiometric coefficients arrays
         stoich_r_pack(:,1:ns) = stoich_r_pack(:,new_order)
         stoich_p_pack(:,1:ns) = stoich_p_pack(:,new_order)
         nudiff_pack  (:,1:ns) = nudiff_pack  (:,new_order)

         i_stoich_r(:,1:ns) = i_stoich_r(:,new_order)
         i_stoich_p(:,1:ns) = i_stoich_p(:,new_order)
         i_nudiff  (:,1:ns) = i_nudiff  (:,new_order)

c        Species array is not permutated as that is the official
c        species order
!        character(len=18), dimension(:), allocatable :: specie

         end subroutine permutate_mechanism_species

           !   *********************************************************
           !   **   Compute optimal species ordering                  **
           !   **                                                     **
           !   **   Author:      Federico Perini                      **
           !   **   Last update: wednesday, 23/11/2011                **
           !   *********************************************************

           subroutine compute_optimal_ordering

              use utilities,        only: force_allocate
              use sparse_definitions
              use sparse_algebra,   only: identity
              use sparse_chemistry, only: JAC_sparse
              use chemistry_setup,  only: permutate_species

              implicit none

              type(sparse) :: tmp_mat
              integer      :: j

              if (.not.allocated(JAC_sparse)) then
                 write(*,*)'Error: cannot compute permutations'
                 write(*,*)'Jacobian has to be initialised'
                 stop
              endif

              ! Diagonal elements must be nonzero!
              tmp_mat = JAC_sparse + identity(ns+1,tiny(0.e0_dp))

              call force_allocate(species_permutations,         ns+1)
              call force_allocate(species_inverse_permutations, ns+1)

              species_permutations(1:tmp_mat%nr)
     &                                     = optimal_ordering(tmp_mat)
              species_inverse_permutations
     &                            (species_permutations(1:tmp_mat%nr))
     &                                     = [(j, j=1,ns+1)]

              ! Deallocate temporary storage
              call deallocate(tmp_mat)


           end subroutine compute_optimal_ordering

c        **************************************************************
c        **     Check for duplicate reactions                        **
c        **************************************************************

         subroutine check_duplicate_reactions

         use sparse_chemistry, only: stoich_r_sp, stoich_p_sp

         implicit none

         logical, dimension(nr) :: tmp_dup
         integer :: i, j
         logical :: check_reac, check_prod
         logical, parameter     :: verbose = .false.

         character(len=*), parameter ::
     &     fmt_strt = "(' Checking for duplicate reactions... ')",
     &     fmt_dup  = "(' Reactions no. ',I5,' and ',I5,' are ',"//
     &                "'duplicate')"


c        Initialise temporary duplicate reaction array
         tmp_dup = .false.

         if (verbose) write(*,fmt_strt)

c        Start check
         dup_check: do i = 1, nr
            inner_reac: do j = i + 1, nr

c             Reactants are equal if both the participating species
c             and the stoichiometric coefficients are the same
              check_reac = (stoich_r_sp%IA(i+1)-stoich_r_sp%IA(i) ==
     &                      stoich_r_sp%IA(j+1)-stoich_r_sp%IA(j) )
              if (check_reac)
     &        check_reac = all( stoich_r_sp%A(
     &                       stoich_r_sp%IA(i):stoich_r_sp%IA(i+1)-1)
     &                   ==  stoich_r_sp%A(
     &                       stoich_r_sp%IA(j):stoich_r_sp%IA(j+1)-1)
     &                 .and. stoich_r_sp%JA(
     &                       stoich_r_sp%IA(i):stoich_r_sp%IA(i+1)-1)
     &                   ==  stoich_r_sp%JA(
     &                       stoich_r_sp%IA(j):stoich_r_sp%IA(j+1)-1) )

c             Products are equal if both the participating species
c             and the stoichiometric coefficients are the same
              check_prod = (stoich_p_sp%IA(i+1)-stoich_p_sp%IA(i) ==
     &                      stoich_p_sp%IA(j+1)-stoich_p_sp%IA(j) )
              if (check_prod)
     &        check_prod = all( stoich_p_sp%A(
     &                       stoich_p_sp%IA(i):stoich_p_sp%IA(i+1)-1)
     &                   ==  stoich_p_sp%A(
     &                       stoich_p_sp%IA(j):stoich_p_sp%IA(j+1)-1)
     &                 .and. stoich_p_sp%JA(
     &                       stoich_p_sp%IA(i):stoich_p_sp%IA(i+1)-1)
     &                   ==  stoich_p_sp%JA(
     &                       stoich_p_sp%IA(j):stoich_p_sp%IA(j+1)-1) )

              found_dup: if (check_reac .and. check_prod) then

                  if (verbose) write(*,fmt_dup)i,j
                  tmp_dup(i) = .true.
                  tmp_dup(j) = .true.
              endif found_dup

            end do inner_reac

         end do dup_check

c        Finalise results into public variable
         if (verbose) write(*,*)
         if (allocated(duplicate)) deallocate(duplicate)
         allocate(duplicate(nr))

         duplicate = tmp_dup

         end subroutine check_duplicate_reactions

      end module speedchem

      module find_mod

      use working_precision, only: dp
      integer, allocatable, public :: indices(:)
      integer, allocatable, public :: i2D1(:),i2D2(:)


!$OMP THREADPRIVATE(indices,i2D1,i2D2)
      contains

      subroutine find_indices(mask)
      logical, intent(in) :: mask(:)
      integer :: i

      if (allocated(indices)) then
         deallocate(indices)
      end if

      allocate(indices(count(mask)))

      indices = pack( [(i,i=1,size(mask))], mask)

      end subroutine find_indices

c     ** Finding indices for 2D matrix
      subroutine find_indices2D(mask)
      logical, intent(in) :: mask(:,:)
      integer :: i,j
      integer :: icol(size(mask,2)),irow(size(mask,1))
      integer :: mcol(size(mask,1),size(mask,2))
      integer :: mrow(size(mask,1),size(mask,2))


      irow(:) = [(i,i=1,size(mask,1))]
      icol(:) = [(i,i=1,size(mask,2))]

      mrow = spread(irow, 2, size(mask,2))
      mcol = spread(icol, 1, size(mask,1))


      if (allocated(i2D1)) then
         deallocate(i2D1)
         deallocate(i2D2)
      end if

      allocate(i2D1(count(mask)))
      allocate(i2D2(count(mask)))

      i2D1 = pack( mrow, mask)
      i2D2 = pack( mcol, mask)

      end subroutine find_indices2D



      end module find_mod





c     *****************************************************************
c     **                                                             **
c     **                     SpeedCHEM FORTRAN                       **
c     **                                                             **
c     **               Physical universal constants                  **
c     **                                                             **
c     **   Author:      Federico Perini                              **
c     **   Last update: tuesday, 30/11/2011                          **
c     **                                                             **
c     *****************************************************************

      module universal_constants

      use working_precision, only: dp

      implicit none
      public

c     Real number constants
      real (dp), parameter :: zero  = 0.0_dp,
     &                        one   = 1.0_dp,
     &                        two   = 2.0_dp,
     &                        three = 3.0_dp,
     &                        four  = 4.0_dp,
     &                        five  = 5.0_dp,
     &                        six   = 6.0_dp,
     &                        seven = 7.0_dp,
     &                        eight = 8.0_dp,
     &                        nine  = 9.0_dp,
     &                        ten   = 10.0_dp,
     &                        twlve = 12.0_dp


c     Real number fractions
      real (dp), parameter :: u2 = one/two,
     &                        u3 = one/three,
     &                        u4 = one/four,
     &                        u5 = one/five,
     &                        u6 = one/six,
     &                        u10= one/ten,
     &                        u12= one/twlve,
     &                        u20= u2/ten

c     Orders of magnitude
      real (dp), parameter :: nano    = 1.0e-09_dp,
     &                        micro   = 1.0e-06_dp,
     &                        milli   = 1.0e-03_dp,
     &                        kilo    = 1.0e+03_dp,
     &                        mega    = 1.0e+06_dp,
     &                        tera    = 1.0e+09_dp

c     Natural logarithms
      real (dp)       , parameter ::  ln10 = log(ten)
      real (dp)       , parameter :: uln10 = one/log(ten)


c     ** Unit Conversion constants *************************************

c     Pressure
      real (dp)       , parameter :: bar_to_dyncm2 = 1.0e06_dp
      real (dp)       , parameter :: dyncm2_to_Pa  = 1.0e-01_dp
      real (dp)       , parameter :: bar_to_Pa     = 1.0e05_dp
      real (dp)       , parameter :: atm_to_Pa     = 101325.0_dp

c     Energy
      real (dp)       , parameter :: cal_to_joule  = 4.184e0_dp
      real (dp)       , parameter :: joule_to_cal  = one/cal_to_joule

      real (dp)       , parameter :: joule_to_kcal = milli*joule_to_cal
      real (dp)       , parameter :: kcal_to_joule = kilo*cal_to_joule

      real (dp)       , parameter :: joule_to_erg  = ten*mega
      real (dp)       , parameter :: erg_to_joule  = u10*micro

c     ** Physical constants ********************************************

c     Universal gas constant [J / mol K]
      real (dp)       , parameter :: R         = 8.31446210000000_dp
      real (dp)       , parameter :: uR        = one/R

c     Universal gas constant, calorie units [cal / mol K]
      real (dp)       , parameter :: Rcal      = R * joule_to_cal
      real (dp)       , parameter :: uRcal     = one/Rcal

c     Universal gas constant, CGS units [erg / mol K]
      real (dp)       , parameter :: Rerg      = R * joule_to_erg
      real (dp)       , parameter :: uRerg     = one/Rerg

c     Standard pressure [Pa]
      real (dp)       , parameter :: Patm      = atm_to_Pa
      real (dp)       , parameter :: uPatm     = one/Patm


c       ** Definition of a type for physical constants ****************
c          (each unit system will correspond to a variable containing
c           the constant values)
        type physchem_const

c          Atomic mass unit
           real (dp)        :: mu

c          Avogadro's number
           real (dp)        :: NA, uNA

c          Boltzmann constant
           real (dp)        :: kB, ukB

c          Gas constant
           real (dp)        :: R, uR

        end type physchem_const

c       ***************************************************************
        type(physchem_const) :: SIpc
        type(physchem_const) :: CGSpc

        contains

c       Setting universal constants
        subroutine set_universal_const
        implicit none

           SIpc%mu = 1.66053892173e-27_dp ! Atomic mass unit  [kg]
           SIpc%NA = 6.0221412927e+23_dp  ! Avogadro's number [1/mol]
           SIpc%kB = 1.380648813e-23_dp   ! Boltzmann's const [J/K]
           SIpc%R  = SIpc%NA*SIpc%kB      ! Gas constant      [J/mol/K]


c          Reciprocal values
           SIpc%uNA = one/SIpc%NA
           SIpc%ukB = one/SIpc%kB
           SIpc%uR  = one/SIpc%R


        end subroutine set_universal_const


      end module universal_constants



c     ** Module for thermodynamic database ****************************
      module SCthermodata

      use working_precision, only: dp
      implicit none
      public

c     all coefficient in thermo database have number of elements
c     equal to the total number of species
      real (dp)       , dimension(:), allocatable , target ::
     &                                               tsw, aL, bL, cL,
     &                                               dL, eL, fL, gL,
     &                                               aH, bH, cH, dH,
     &                                               eH, fH, gH

c     Computationally efficient coefficient storage
      real (dp)       , dimension(:,:), allocatable :: gibbsL, gibbsH


c     Arrays for tabulated data (faster computation)
      real (dp)       , dimension(:,:), allocatable :: tab_CpsuR,
     &                                                 tab_HsuRT,
     &                                                 tab_SsuR,
     &                                                 tab_dHdt,
     &                                                 tab_dUdt,
     &                                                 tab_dGdt,
     &                                                 tab_dDG0dt,
     &                                                 tab_dCvdT


c     Arrays for tabulated equilibrium data (faster computation)
c     uKp = reciprocal of equilibrium constant (pressure formulation)
c     uKc = reciprocal of equilibrium constant (concentration formulat.)
      real (dp)       , dimension(:,:), allocatable :: tab_uKp, tab_uKc,
     &                                                 tab_dKcdT

c     Element labels and number of elements within the species
      integer, dimension(:), allocatable :: nels,nel1,nel2,nel3,nel4
      character*2, dimension(:), allocatable :: el1,el2,el3,el4

c     Elements matrix has dimension (EM(ns,nel))
      integer, dimension(:,:), allocatable :: EM
      real (dp)       , dimension(:,:), allocatable :: Emfr

c     Indices of most important elements
      integer :: iC, iO, iH, iN, iAr

      contains
c        **************************************************************

c        **************************************************************
c           Get element indices
c        **************************************************************
         function get_el(el) result(iel)

         use speedchem, only: elementi, nel
         use utilities, only: find_stringi_in_array

         implicit none

         character(len=*), intent(in) :: el
         integer                      :: iel, i

         character(len=*), parameter  ::
     &     fmt_er = "(' Element not supported in get_el: ',A2)",
     &     fmt_e2 = "(' Element not found in mechanism : ',A2)",
     &     fmt_le = "(' Available elements:',100(1x,A2))"

         select case (trim(adjustl(el)))

            case ("C","c")

             if (.not.iC>0) iC = find_stringi_in_array(el,elementi)
             iel = iC

            case ("H","h")

             if (.not.iH>0) iH = find_stringi_in_array(el,elementi)
             iel = iH

            case ("O","o")

             if (.not.iO>0) iO = find_stringi_in_array(el,elementi)
             iel = iO


            case ("N","n")

             if (.not.iN>0) iN = find_stringi_in_array(el,elementi)
             iel = iN

            case ("AR","ar","Ar","aR")

             if (.not.iAr>0) iAr = find_stringi_in_array(el,elementi)
             iel = iAr

            case default
              write(*,fmt_er)trim(adjustl(el))
              stop
         end select

         if (iel == 0) then
            write(*,fmt_e2)trim(adjustl(el))
            write(*,fmt_le)(elementi(i),i=1,nel)
            stop
         endif

         return
         end function get_el

c        **************************************************************

c        **************************************************************
c           Allocate computationally efficient storage
c        **************************************************************
         subroutine store_thermo_coeffs

         use speedchem, only: ns
         use universal_constants, only: u2, u6, u12, u20
         implicit none

         integer :: is

c        Coefficients for non-dimensional gibbs free energy ***********
         allocate( gibbsL(7, ns), gibbsH(7, ns) )

         do is = 1, ns

            gibbsL(1,is) = aL(is) - gL(is)
            gibbsL(2,is) = - u2 * bL(is)
            gibbsL(3,is) = - u6 * cL(is)
            gibbsL(4,is) = - u12 * dL(is)
            gibbsL(5,is) = - u20 * eL(is)
            gibbsL(6,is) = fL(is)
            gibbsL(7,is) = - aL(is)

            gibbsH(1,is) = aH(is) - gH(is)
            gibbsH(2,is) = - u2 * bH(is)
            gibbsH(3,is) = - u6 * cH(is)
            gibbsH(4,is) = - u12 * dH(is)
            gibbsH(5,is) = - u20 * eH(is)
            gibbsH(6,is) = fH(is)
            gibbsH(7,is) = - aH(is)

         end do


         end subroutine store_thermo_coeffs

c        **************************************************************
c          Permutates all the arrays containing species thermodynamic
c          properties into a new species ordering previously defined
c        **************************************************************
         subroutine permutate_thermo(new_order)
         implicit none

         integer, dimension(:), intent(in) :: new_order
         integer                           :: j

c        Preliminary check
         if (size(new_order) /= size(tsw)) then
           write(*,*)'Wrong indexing array in permutate_thermo'
           write(*,*)size(new_order),size(tsw)
           stop
         endif

c        Permutate JANAF tables coefficients
         tsw = tsw(new_order)
         aL  = aL (new_order)
         bL  = bL (new_order)
         cL  = cL (new_order)
         dL  = dL (new_order)
         eL  = eL (new_order)
         fL  = fL (new_order)
         gL  = gL (new_order)
         aH  = aH (new_order)
         bH  = bH (new_order)
         cH  = cH (new_order)
         dH  = dH (new_order)
         eH  = eH (new_order)
         fH  = fH (new_order)
         gH  = gH (new_order)


c        Permutate element data
         do j = 1, size(EM, 1)
           EM  (j, :) = EM(new_order(j), :)
           Emfr(j, :) = Emfr(new_order(j), :)
         end do


         end subroutine permutate_thermo


c        **************************************************************
c        **                                                          **
c        **                ELEMENT MATRIX COMPUTATION                **
c        **                                                          **
c        **      Costruzione della matrice degli elementi E per il   **
c        **             problema dell equilibrio chimico             **
c        **                                                          **
c        **************************************************************
         subroutine element_matrix

         use speedchem,    only : specie, elementi, nel, ns, nr
!      use SCthermodata, only : EM, nels,nel1,nel2,nel3,nel4,
!     &                             el1,el2,el3,el4
         use sparse_chemistry, only: EM_sp
         use sparse_algebra,   only: dense_to_sparse

         implicit none

c        Local variables **********************************************
         integer :: isp, ie, iee, ielassign
         integer, dimension(4) :: ndummy
         character*2, allocatable, dimension(:) :: dummy

c        **************************************************************

         allocate(dummy(4))

         if (.not.allocated(EM)) then
           allocate(EM(ns,nel))
           EM(:,:) = 0
         else
           call dense_to_sparse(real(EM, dp),EM_sp)
           return
         endif

c        Reducing mechanism element labels into lowercase characters
         do ie = 1,nel
            call lower_case(elementi(ie))
         end do

c        Species main loop
         do isp = 1,ns

c           Sono memorizzati 4 elementi per ogni specie
            dummy(:) = '  '
            call lower_case(el1(isp))
            call lower_case(el2(isp))
            call lower_case(el3(isp))
            call lower_case(el4(isp))
            if (nels(isp).gt.0)read(el1(isp),*)dummy(1)
            if (nels(isp).gt.1)read(el2(isp),*)dummy(2)
            if (nels(isp).gt.2)read(el3(isp),*)dummy(3)
            if (nels(isp).gt.3)read(el4(isp),*)dummy(4)

            ndummy(1) = nel1(isp)
            ndummy(2) = nel2(isp)
            ndummy(3) = nel3(isp)
            ndummy(4) = nel4(isp)

c           Every species can have up to four elements
            do ie = 1,nels(isp)

c              Look for the element name
               ielassign = 0
               do iee = 1,nel
                  if(dummy(ie).eq.elementi(iee)) ielassign = iee
               end do

c              Storing element number
               EM(isp,ielassign) = ndummy(ie)

            end do

            if (sum(EM(isp,:)).eq.0) then
              write(*,*)'Error: species ',specie(isp),' has 0 elements.'
              stop
            endif



         end do

c        Store element matrix in sparse form
         call dense_to_sparse(real(EM, dp),EM_sp)

         deallocate(dummy)

      end subroutine element_matrix

c        **************************************************************
c        **             ELEMENT MASS FRACTIONS COMPUTATION           **
c        **   Computing mass fractions of elements in each species   **
c        **************************************************************
         subroutine element_mass_fracs
         use speedchem,    only : specie, elementi, nel, ns, nr

         implicit none

         integer :: i, isp, iel
         real (dp)       , dimension(nel) :: amatom
         integer, parameter :: NATOM = 102
         real (dp)       , dimension(NATOM) :: ATOM
         CHARACTER :: IATOM(NATOM)*2, dummy*2, tempel*2
C
         DATA (IATOM(I),ATOM(I),I=1,40) /
     *   'H ',  1.00797, 'HE',  4.00260, 'LI',  6.93900, 'BE',  9.01220,
     *   'B ', 10.81100, 'C ', 12.01115, 'N ', 14.00670, 'O ', 15.99940,
     *   'F ', 18.99840, 'NE', 20.18300, 'NA', 22.98980, 'MG', 24.31200,
     *   'AL', 26.98150, 'SI', 28.08600, 'P ', 30.97380, 'S ', 32.06400,
     *   'CL', 35.45300, 'AR', 39.94800, 'K ', 39.10200, 'CA', 40.08000,
     *   'SC', 44.95600, 'TI', 47.90000, 'V ', 50.94200, 'CR', 51.99600,
     *   'MN', 54.93800, 'FE', 55.84700, 'CO', 58.93320, 'NI', 58.71000,
     *   'CU', 63.54000, 'ZN', 65.37000, 'GA', 69.72000, 'GE', 72.59000,
     *   'AS', 74.92160, 'SE', 78.96000, 'BR', 79.90090, 'KR', 83.80000,
     *   'RB', 85.47000, 'SR', 87.62000, 'Y ', 88.90500, 'ZR', 91.22000/
C
         DATA (IATOM(I),ATOM(I),I=41,80) /
     *   'NB', 92.90600, 'MO', 95.94000, 'TC', 99.00000, 'RU',101.07000,
     *   'RH',102.90500, 'PD',106.40000, 'AG',107.87000, 'CD',112.40000,
     *   'IN',114.82000, 'SN',118.69000, 'SB',121.75000, 'TE',127.60000,
     *   'I ',126.90440, 'XE',131.30000, 'CS',132.90500, 'BA',137.34000,
     *   'LA',138.91000, 'CE',140.12000, 'PR',140.90700, 'ND',144.24000,
     *   'PM',145.00000, 'SM',150.35000, 'EU',151.96000, 'GD',157.25000,
     *   'TB',158.92400, 'DY',162.50000, 'HO',164.93000, 'ER',167.26000,
     *   'TM',168.93400, 'YB',173.04000, 'LU',174.99700, 'HF',178.49000,
     *   'TA',180.94800, 'W ',183.85000, 'RE',186.20000, 'OS',190.20000,
     *   'IR',192.20000, 'PT',195.09000, 'AU',196.96700, 'HG',200.59000/
C
         DATA (IATOM(I),ATOM(I),I=81,NATOM) /
     *   'TL',204.37000, 'PB',207.19000, 'BI',208.98000, 'PO',210.00000,
     *   'AT',210.00000, 'RN',222.00000, 'FR',223.00000, 'RA',226.00000,
     *   'AC',227.00000, 'TH',232.03800, 'PA',231.00000, 'U ',238.03000,
     *   'NP',237.00000, 'PU',242.00000, 'AM',243.00000, 'CM',247.00000,
     *   'BK',249.00000, 'CF',251.00000, 'ES',254.00000, 'FM',253.00000,
     *   'D ',002.01410, 'E',5.45D-4/

c        **************************************************************

c        Initialise matrix of elemental mass fractions
         if (.not.allocated(Emfr))allocate(Emfr(ns,nel))
         Emfr = 0.e0_dp


c        Find out the atomic mass numbers of the elements considered in
c        current problem
         mechanism_elements: do iel = 1,nel

c           Read element in lower case
            call lower_case(elementi(iel))
            read(elementi(iel),*)tempel

            periodic_table_elements: do i = 1,natom

               call lower_case(IATOM(i))
               read(IATOM(i),*)dummy


               if (trim(dummy).eq.trim(tempel)) then
                  amatom(iel) = ATOM(i)
                  exit periodic_table_elements
               endif

            end do periodic_table_elements
         end do mechanism_elements

c        Normalizing atom masses over the mass of each specie
         mechanism_species: do isp = 1,ns
            Emfr(isp,:) = amatom * EM(isp,:)
            Emfr(isp,:) = Emfr(isp,:)/sum(Emfr(isp,:))
         end do mechanism_species


      end subroutine element_mass_fracs


c        ***************************************************************
c        **                                                           **
c        **     Check atom conservation in the mechanism              **
c        **                                                           **
c        ***************************************************************

         subroutine check_atom_conservation
         use speedchem, only: nr, ns, nel, stoich_r, stoich_p,
     &                        species, elementi
         use sparse_chemistry, only: stoich_r_sp, stoich_p_sp, EM_sp
         use sparse_algebra


         implicit none

         integer :: atom_number_products (nel),
     &              atom_number_reactants(nel),
     &              ispp(ns), ispr(ns)
         integer :: i, j, k, l, nspp, nspr

         type(sparse) :: atoms_prod, atoms_reac ! dimension: (nr x nel)

         logical :: wrong_found = .false.

         character(len=*), parameter ::
     &     fmt_er = "(' Missing atom conservation in reaction ',I5)",
     &     fmt_r  = "(' Atoms in reactants: ',(100(A2,'=',1X,I2,1x)))",
     &     fmt_p  = "(' Atoms in products : ',(100(A2,'=',1X,I2,1x)))",
     &     fmt_ok = "(' Element conservation check passed.')"


c          Compute element counts for each reaction, in reactants and
c          products
           call sparse_symbolic_mm(stoich_p_sp, EM_sp, atoms_prod)
           call sparse_2_matmul   (stoich_p_sp, EM_sp, atoms_prod)

           call sparse_symbolic_mm(stoich_r_sp, EM_sp, atoms_reac)
           call sparse_2_matmul   (stoich_r_sp, EM_sp, atoms_reac)

           reactions_loop: do i = 1, nr

               atom_number_reactants = 0
               atom_number_products  = 0
               do j = atoms_prod%IA(i), atoms_prod%IA(i+1)-1
               atom_number_products(atoms_prod%JA(j)) =
     &                                            int(atoms_prod%A(j))
               end do

               do k = atoms_reac%IA(i), atoms_reac%IA(i+1)-1
               atom_number_reactants(atoms_reac%JA(k)) =
     &                                            int(atoms_reac%A(k))
               end do

              if (any(atom_number_reactants /=
     &                atom_number_products     )) then

                 wrong_found = .true.

                 write( *, fmt_er) i

                 write( *, fmt_r )(elementi(k),
     &                             atom_number_reactants(k),k=1,nel)

                 write( *, fmt_p )(elementi(k),
     &                             atom_number_products (k),k=1,nel)

              endif

           end do reactions_loop

           if (wrong_found) stop

           write(*,fmt_ok)

         end subroutine check_atom_conservation


c        ***************************************************************
c        **  Check polynomial coefficiens at the interface between    **
c        **  low- and high-temperature ranges, for accuracy           **
c        **                                                           **
c        **  Author: Federico Perini                                  **
c        **  Date  : 26/07/2012                                       **
c        **                                                           **
c        ***************************************************************
         subroutine check_janaf_polynomials
         use speedchem,           only: ns, specie
ck2015 mechdir 1 line
         use chemistry_setup,     only: mechdir
         use universal_constants, only: u2,u3,u4,u5
         use utilities,           only: relative_error
         implicit none

         logical,          parameter     :: verbose = .true.
         real (dp)       , parameter     :: tolerance = 1.e-02_dp
         integer,          parameter     :: lowT  = 300,
     &                                      highT = 3000,
     %                                      stepT = 10


         integer                         :: n_error, j, i
         logical,          dimension(ns) :: unaccurate
         real (dp)                       :: T, Ta(6)
         real (dp)       , dimension(ns) :: LCp, LU, LS, LH,
     &                                      HCp, HU, HS, HH
         real (dp)                       :: tCp, tU, tS, tH
         real (dp)       , dimension(ns) :: Cp_error, U_error, S_error,
     &                                      H_error

         real (dp)       ,   pointer     :: a,b,c,d,e,f,g

         character(len=28)           :: nfile
         character(len=*), parameter ::
     &     fmt_scn = "(' There are warnings on JANAF polynomials. ')",
     &     fmt_wrn = "(' JANAF polynomial warnings: ',I3,' found. ')",
     &     fmt_spc = "(' Species ',I4,': ',A18)",
     &     fmt_pct = "(' ',A15,' values differ by ',F5.2,'%')",
     &     fmt_tmp = "(' at low-high temperature switch: ',F7.2,' K')",
     &     fmt_nfi = "('dat.janaf_',A18)",
     &     fmt_hed = "('  Temp [K]    Cp/R   U/RT    S/R    H/RT ')",
     &     fmt_dat = "(5(1x,1pE12.5,1x))"


c        Compute SPECIFIC HEAT values (Cp/R)
         LCp  = aL +bL*tsw +cL*tsw**2 +dL*tsw**3 +eL*tsw**4
         HCp  = aH +bH*tsw +cH*tsw**2 +dH*tsw**3 +eH*tsw**4

c        Compute internal energy values (U/RT)
         LU   =(aL - 1.e0_dp) + u2*bL*tsw + cL*u3*tsw**2
     &         + dL*u4*tsw**3 + eL*u5*tsw**4 + fL/tsw
         HU   =(aH - 1.e0_dp) + u2*bH*tsw + cH*u3*tsw**2
     &         + dH*u4*tsw**3 + eH*u5*tsw**4 + fH/tsw

c        Compute entropy values (S/R)
         LS =  aL*log(tsw) + bL*tsw + u2*cL*tsw**2 + dL*u3*tsw**3
     &      +  u4*eL*tsw**4 + gL
         HS =  aH*log(tsw) + bH*tsw + u2*cH*tsw**2 + dH*u3*tsw**3
     &      +  u4*eH*tsw**4 + gH

c        Compute enthalpy values (H/RT)
         LH = aL + u2*bL*tsw + cL*u3*tsw**2 + dL*u4*tsw**3
     &      + u5*eL*tsw**4 + fL/tsw
         HH = aH + u2*bH*tsw + cH*u3*tsw**2 + dH*u4*tsw**3
     &      + u5*eH*tsw**4 + fH/tsw
ck2015
         return

c        Compute relative differences at switch temperature
         Cp_error = relative_error(LCp,HCp)
         U_error  = relative_error(LU, HU )
         S_error  = relative_error(LS, HS )
         H_error  = relative_error(LH, HH )

c        Evaluate unaccurate polynomials
         unaccurate = (Cp_error > tolerance) .or.
     &                (U_error  > tolerance) .or.
     &                (H_error  > tolerance) .or.
     &                (S_error  > tolerance)
         n_error    = count(unaccurate)

c        Unaccurate values found
         warnings_found: if (n_error > 0) then

            write(*,fmt_scn)

ck2015            open(unit = 13, file = 'dat.janafwarn')
            open(unit = 13, file = trim(mechdir)//"dat.janafwarn")

            write(13,      *)
            write(13,fmt_wrn)n_error


            do j = 1, ns
               if (unaccurate(j)) then
                  write(13,      *)
                  write(13, fmt_spc)j, specie(j)
                  write(13, fmt_pct)'specific heat  ',1.d2*Cp_error(j)
                  write(13, fmt_pct)'internal energy',1.d2* U_error(j)
                  write(13, fmt_pct)'entropy        ',1.d2* S_error(j)
                  write(13, fmt_pct)'enthalpy       ',1.d2* H_error(j)
                  write(13, fmt_tmp)tsw(j)
               endif
            end do

            close(13)

C           If requested, print temperature polynomials to files
            prtdetails: if (verbose) then

               do j = 1, ns
                  if (unaccurate(j)) then
                     write(nfile,fmt_nfi)adjustl(specie(j))
ck2015                     open (unit = 13, file = nfile)
                     open (unit = 13, file = trim(mechdir)//nfile)

                     write(13,fmt_hed)
                     temploop2: do i = lowT, highT, stepT
                       T  = real(i, dp)

                     if (T>=tsw(j)) then
                        a => aH(j)
                        b => bH(j)
                        c => cH(j)
                        d => dH(j)
                        e => eH(j)
                        f => fH(j)
                        g => gH(j)
                     else
                        a => aL(j)
                        b => bL(j)
                        c => cL(j)
                        d => dL(j)
                        e => eL(j)
                        f => fL(j)
                        g => gL(j)
                     endif

                     tU  = (a-1.e0_dp) + u2*b*T + c*u3*T**2
     &                   + d*u4*T**3 + e*u5*T**4 + f/T

                     tCp = a + b*T +c *T**2 +d*T**3 +e*T**4

                     tS  = a*log(T) + b*T + u2*c*T**2
     &                   + d*u3*T**3 +  u4*e*T**4 + g

                     tH  = a + u2*b*T + c*u3*T**2 + d*u4*T**3
     &                   + u5*e*T**4 + f/T
!
!              elsewhere
!
!                tU = (aL-1.e0_dp)*T + u2*bL*T**2 + cL*u3*T**3
!     &               + dL*u4*T**4 + eL*u5*T**5 + fL
!
!                tCp = aL +bL*T +cL*T**2 +dL*T**3 +eL*T**4
!
!                tS  = aL*log(T) + bL*T + u2*cL*T**2
!     &              + dL*u3*T**3 +  u4*eL*T**4 + gL
!
!                tH  = aL + u2*bL*T + cL*u3*T**2 + dL*u4*T**3
!     &              + u5*eL*T**4 + fL/T


!              end where


                       write(13,fmt_dat) T, tCp, tU, tS, tH
                     end do temploop2



                     close(13)

                  endif
               end do




            endif prtdetails


         end if warnings_found
         end subroutine check_janaf_polynomials



c        *****************************************************************
c        Writes JANAF thermodynamic coefficients card to file
c        isp   = number of the species to print out
c        ifile = logical unit number of the file
c
c        NB every species here is considered as gaseous
c
c        Example
!ar                120186ar  1               g  0300.00   5000.00  1000.00      1
! 0.02500000e+02 0.00000000e+00 0.00000000e+00 0.00000000e+00 0.00000000e+00    2
!-0.07453750e+04 0.04366001e+02 0.02500000e+02 0.00000000e+00 0.00000000e+00    3
! 0.00000000e+00 0.00000000e+00-0.07453750e+04 0.04366001e+02                   4
c        *****************************************************************

         subroutine write_janaf_card(isp,ifile)

         use speedchem, only: elementi, specie, nel

         implicit none

         integer, intent(in) :: isp, ifile

         character(len=24) :: spname
         character(len=2)  :: elname(4)
         character(len=3)  :: elcount(4)

         integer           :: i, j, nelms, elemsn(4), ielems(4)

         character(len=*), parameter ::
     &     fmt_line1 = "(A24,4(A2,A3),A1,3(2X,F8.3),4x,I1)",
     &     fmt_line2 = "(5(1PE15.8),4x,I1)",
     &     fmt_line4 = "(4(1PE15.8),19x,I1)"

c          Initialise arrays
           spname(1:24) = ' '
           elname(1:4)(1:2) = ' '

c          Full species name
           spname = trim(adjustl(specie(isp)))
           spname(len(trim(spname))+1:24) = ' '

c          Element numbers
           nelms = count(EM(isp,:)/=0)
           ielems(1:nelms) = pack([(j,j=1,nel)],mask=EM(isp,:)/=0)

           do i = 1, nelms

              elemsn(i) = EM(isp, ielems(i))
              elname(i)(1:2) = elementi(ielems(i))

              if (elemsn(i) > 0) then
                 write(elcount(i)(1:3),"(I3)")elemsn(i)
              else
                 elcount(i)(1:3) = ' '
              endif

           end do

           do i = nelms+1,4
              elcount(i)(1:3) = ' '
              elname(i)(1:2)  = ' '
           end do

           write(ifile,fmt_line1)spname,elname(1)(1:2), elcount(1)(1:3),
     &                                  elname(2)(1:2), elcount(2)(1:3),
     &                                  elname(3)(1:2), elcount(3)(1:3),
     &                                  elname(4)(1:2), elcount(4)(1:3),
     &                                  'g',300.e0_dp, 5000.e0_dp,
     &                                  tsw(isp), 1

           j = isp

           write(ifile,fmt_line2) aH(j), bH(j), cH(j), dH(j), eH(j), 2
           write(ifile,fmt_line2) fH(j), gH(j), aL(j), bL(j), cL(j), 3
           write(ifile,fmt_line4) dL(j), eL(j), fL(j), gL(j),        4

         end subroutine write_janaf_card

c        ***************************************************************
c        Writes JANAF thermodynamic coefficients for the whole mechanism
c        ***************************************************************

         subroutine scthermo_to_janaf

ck2015         use chemistry_setup, only: mechanism
         use chemistry_setup, only: mechanism, mechdir
         use speedchem, only: ns

         implicit none

         character(len=*), parameter ::
     &     fmt_thermo = "('THERMO')",
     &     fmt_end    = "('END')",
     &     fmt_name   = "('!',A80)",
     &     fmt_trange = "(1x,3(F10.3))"

         integer :: ifile, j

c        Open file unit
         ifile = 30
ck2015         open(unit=ifile, file = 'dat.cktherm')
         open(unit=ifile, file = trim(mechdir)//"dat.cktherm")

c           Write mechanism header
            write(ifile,fmt_name)mechanism

c           Write temperature ranges for THERMO ALL data
            write(ifile,fmt_thermo)
            write(ifile,fmt_trange)3.0d2,1.0d3,5.0d3

c           Write JANAF card for each species
            do j = 1, ns
               call write_janaf_card(j,ifile)
            end do

            write(ifile,fmt_end)

         close(ifile)

         end subroutine scthermo_to_janaf


c     *****************************************************************
c        Compute indices for tabulated data structures
c     *****************************************************************
      subroutine table_indices(T, iT, frac)

      use chemistry_setup, only: Temp_LOlim, rec_Ttable_accuracy,
     &                           tab_nsteps

      implicit none

      real (dp)       , intent(in)   :: T
      integer,          intent(out)  :: iT
      real (dp)       , intent(out)  :: frac
      real (dp)                      :: tenthT

         tenthT = rec_Ttable_accuracy * T
         iT     = int(tenthT-Temp_LOlim*rec_Ttable_accuracy)+1
         frac   =    (tenthT-int(tenthT))

!         write(*,*)T,iT,frac,rec_Ttable_accuracy
!         pause

!         if (.not.(iT>=1.and.iT<=tab_nsteps)) then
!           write(* ,900)iT,T
!           stop
!         endif

 900  format(' error in thermodynamic data calculation: T(',I3,') = ',
     &         F7.2,' K')
      end subroutine table_indices

c     *****************************************************************
c        Compute fractions for degree 3 interpolation of tables
c        (4 tabulated points needed)
c     *****************************************************************
      function interp3_coefs(frac) result(coefs)

      use universal_constants, only: u2, u3, u6

      implicit none

      real (dp)       , intent(in)   :: frac
      real (dp)       , dimension(4) :: coefs
      real (dp)                      :: frac2

         frac2 = frac * frac

         coefs(1) = frac*(u2*frac-u6*frac2-u3)
         coefs(2) = (u2*frac-1.e0_dp)*(frac2-1.e0_dp)
         coefs(3) = frac*(u2*frac-u2*frac2+1.e0_dp)
         coefs(4) = frac*u6*(frac2-1.e0_dp)

      end function interp3_coefs

c     *****************************************************************
c        Compute fractions for degree 4 interpolation of tables
c        (5 tabulated points needed)
c     *****************************************************************
      function interp4_coefs(frac) result(coefs)

      use universal_constants, only: u2, u3, u4, u12, one, two

      implicit none

      real (dp)       , intent(in)   :: frac
      real (dp)       , dimension(5) :: coefs
      real (dp)                      :: frac2, frac2m1

         frac2   = frac * frac
         frac2m1 = frac2 - 1.e0_dp

         coefs(1) = u12 * frac * frac2m1 * (u2*frac - one)
         coefs(2) = u3  * frac * (frac - one) * (two - u2 * frac2)
         coefs(3) = (u4 * frac2 - one) * frac2m1
         coefs(4) = u3  * frac * (frac + one) * (two - u2 * frac2)
         coefs(5) = u12 * frac * frac2m1 * (u2 * frac + one)

      end function interp4_coefs

c     *****************************************************************
c        Compute temperature array containing temperature functions
c                   Ta = [T, T^2, T^3, T^4, 1/T, ln(T)]
c     *****************************************************************
      function temperature_array(T) result(Ta)
      implicit none

      real (dp)       , intent(in)   :: T
      real (dp)       , dimension(6) :: Ta

         Ta(1) = T           ! T
         Ta(2) = Ta(1) * T   ! T**2
         Ta(3) = Ta(2) * T   ! T**3
         Ta(4) = Ta(3) * T   ! T**4
         Ta(5) = 1.e0_dp  / T   ! 1/T
         Ta(6) = log(T)      ! ln(T)

      end function temperature_array

c     *****************************************************************
c        Logical swithch for decision about accurate or tabulated
c        computation of temperature-dependent properties
c     *****************************************************************
      function use_table(T) result(usetable)
      use chemistry_setup, only: accurate_scthermo, Temp_LOlim,
     &                           Temp_HIlim, Temp_table_accuracy
      implicit none

      real (dp)       , intent(in)   :: T
      logical                        :: usetable


         usetable = .false.
         if (accurate_scthermo) return

         ! Check for temperature bounds
         if (T >= Temp_LOlim + 3*Temp_table_accuracy .and.
     &       T <= Temp_HIlim - 2*Temp_table_accuracy )
     &   usetable = .true.

      end function use_table




c     *****************************************************************
c     **                                                             **
c     **                     SpeedCHEM FORTRAN                       **
c     **                                                             **
c     **    Tabulation of thermodynamic polynomials at 10K steps     **
c     **   also including derivatives with respect to temperature    **
c     **                                                             **
c     **                                                             **
c     **   Author:      Federico Perini                              **
c     **   Last update: tuesday, 28/06/2011                          **
c     **                                                             **
c     *****************************************************************


      subroutine tabulate_SCthermo

      use chemistry_setup, only: tab_nsteps,
     &                      Temp_table_accuracy, Temp_LOlim, Temp_HIlim,
     &                      rec_Ttable_accuracy
      use speedchem,       only: ns, nr

      use sparse_chemistry, only: nudiff_sparse
      use sparse_algebra

      use universal_constants, only: u2, u3, u4, u5, R

      implicit none
      integer :: i, k, nsteps
      real (dp)        :: T, T2, T3, T4, unosuT, lnT
      real (dp)       , dimension(ns) :: C,H,S,dHdT,dUdT,dGdT,dCvdT

      tab_nsteps = int( ( Temp_HIlim - Temp_LOlim )
     &                    / Temp_table_accuracy)

      rec_Ttable_accuracy = 1.e0_dp/Temp_table_accuracy

      allocate(tab_CpsuR (ns,-1:tab_nsteps),tab_HsuRT(ns,-1:tab_nsteps),
     &         tab_SsuR  (ns,-1:tab_nsteps),tab_dGdt (ns,-1:tab_nsteps),
     &         tab_dDG0dt(nr,-1:tab_nsteps),tab_dCvdT(ns,-1:tab_nsteps))


      do i=1,tab_nsteps + 2


         T = Temp_LOlim + real(i-2, dp) * Temp_table_accuracy

         T2     = T  * T
         T3     = T2 * T
         T4     = T3 * T
         unosuT = 1/T
         lnT    = log(T)

         where (tsw.gt.T)
           C = aL + bL*T + cL * T2 + dL * T3 + eL * T4
           H = aL + u2*bL*T + cL*T2*u3 + dL*T3*u4 + u5*eL*T4 +
     &           fL*unosuT
           S =  aL*lnT + bL*T + u2*cL*T2 + dL*u3*T3 + u4*eL*T4
     &           + gL

c          Constant volume specific heat derivative in dT [J/mol/K]
           dCvdT = R * (  bL + 2 * cL * T + 3 * dL * T2 +
     &                    4 * eL * T3 )

c          Enthaply derivative in dT [J/mol/K^2]
!           dHdT = R * (aL + bL * T + cL * T2 + dL * T3 + eL * T4)

c          Gibbs potential derivative in dT [J/mol/K^2]
           dGdT = - R * ( gL + aL * lnT + bL * T + u2 * cL * T2 +
     &                      u3 * dL * T3 + u4 * eL * T4 )

         elsewhere
           C = aH + bH*T + cH * T2 + dH * T3 + eH * T4
           H = aH + u2*bH*T + cH*T2*u3 + dH*T3*u4 + u5*eH*T4 +
     &            fH*unosuT
           S = aH*lnT + bH*T + u2*cH*T2 + dH*u3*T3 + u4*eH*T4
     &          + gH

c          Constant volume specific heat derivative in dT [J/mol/K]
           dCvdT = R * (  bH + 2 * cH * T + 3 * dH * T2 +
     &                    4 * eH * T3 )

c          Enthaply derivative in dT [J/mol/K^2]
!           dHdT = R * (aH + bH * T + cH * T2 + dH * T3 + eH * T4)

c          Gibbs potential derivative in dT [J/mol/K^2]
           dGdT = - R * ( gH + aH * lnT + bH * T + u2 * cH * T2 +
     &                      u3 * dH * T3 + u4 * eH * T4 )

         end where

c        Internal energy derivative in dT [J/mol/K^2]
!         dUdt = dHdT - R

c        Gibbs reaction free energy in dT [J/mol/K^2]

         tab_dDG0dt(1:nr,i-2) = nudiff_sparse * dGdt

         tab_CpsuR(:,i-2) = C
         tab_HsuRT(:,i-2) = H
         tab_SsuR (:,i-2) = S

!         tab_dHdT (:,i-2) = dHdT
!         tab_dUdT (:,i-2) = dUdT
         tab_dGdT (:,i-2) = dGdT

         tab_dCvdT(:,i-2)= dCvdT

      end do


      write(*,900)

 900  format(' thermodynamic data tabulated.')

      end subroutine tabulate_SCthermo

      end module SCthermodata


c     ** Module for thermodynamic properties **************************
      module SCspeciesthermo

      use working_precision, only: dp
      implicit none
      public

c     Species molar quantities
      real (dp)       , dimension(:), allocatable :: Cpmol,Cvmol,Hmol,
     &                                               Smol

c     Species mass specific quantities
      real (dp)       , dimension(:), allocatable :: cpm,cvm,hm


      contains

c        ***************************************************************
c        ** Species nondimensional constant volume specific heats     **
c        ***************************************************************
         subroutine CvuRmol(T,CvuR,iT,frac,mask)

         use SCthermodata,    only:tsw, aL, bL, cL, dL, eL,
     &                                  aH, bH, cH, dH, eH,
     &                                  tab_CpsuR, interp4_coefs
         use chemistry_setup, only: accurate_scthermo
         use universal_constants, only: u2, u3, u6
         use SCthermodata,        only: use_table

         implicit none

         real (dp)       , dimension(6), intent(in)  :: T
         integer, optional, intent(in)               :: iT
         real (dp)       , optional, intent(in)      :: frac
         logical, dimension(size(aL)), intent(in), optional :: mask
         real (dp)       , dimension(size(aL)), intent(out) :: CvuR
         real (dp)        :: frac2, coefs(5)
         logical, parameter :: parabolic = .true.

         accuracy: if (.not.(use_table(T(1))).or.(.not.present(iT)))
     &             then

         if (present(mask)) then

            where (mask)
               where (tsw.gt.T(1) )
                 CvuR = (aL +bL*T(1) +cL*T(2) +dL*T(3) +eL*T(4))-1.e0_dp
               elsewhere
                 CvuR = (aH +bH*T(1) +cH*T(2) +dH*T(3) +eH*T(4))-1.e0_dp
               end where
            elsewhere
              CvuR = 0.e0_dp
            end where

         else ! .not.present(mask)

            where (tsw.gt.T(1) )
              CvuR  = (aL +bL*T(1) +cL*T(2) +dL*T(3) +eL*T(4))-1.e0_dp
            elsewhere
              CvuR  = (aH +bH*T(1) +cH*T(2) +dH*T(3) +eH*T(4))-1.e0_dp
            end where

         endif

         else

         if (present(mask)) then

          where (mask)
            CvuR = ( tab_CpsuR(:,iT-1)*(1.e0_dp-frac)
     &             + tab_CpsuR(:,iT  )*frac )
     &             - 1.e0_dp
          elsewhere
            CvuR = 0.e0_dp
          end where

         else

            tabulation_order: if (parabolic) then

                coefs = interp4_coefs(frac)
                CvuR  = tab_CpsuR(:,iT-3)*coefs(1) +
     &                  tab_CpsuR(:,iT-2)*coefs(2) +
     &                  tab_CpsuR(:,iT-1)*coefs(3) +
     &                  tab_CpsuR(:,iT  )*coefs(4) +
     &                  tab_CpsuR(:,iT+1)*coefs(5)

            else

               CvuR = tab_CpsuR(:,iT-1)*(1.e0_dp-frac) +
     &                tab_CpsuR(:,iT)*frac

            endif tabulation_order

            CvuR = CvuR - 1.e0_dp

         endif


         endif accuracy

         end subroutine CvuRmol

c        ***************************************************************
c        ** Species nondimensional constant volume specific heats     **
c        ***************************************************************
         subroutine CpuRmol(T,CpuR,iT,frac,mask)

         use SCthermodata,    only:tsw, aL, bL, cL, dL, eL,
     &                                  aH, bH, cH, dH, eH,
     &                                  tab_CpsuR, interp4_coefs,
     &                                  use_table
         use chemistry_setup, only: accurate_scthermo
         use universal_constants, only: u2, u3, u6

         implicit none

         real (dp)       , dimension(6), intent(in)  :: T
         integer, optional, intent(in)               :: iT
         real (dp)       , optional, intent(in)      :: frac
         logical, dimension(size(aL)), intent(in), optional :: mask
         real (dp)       , dimension(size(aL)), intent(out) :: CpuR
         real (dp)        :: frac2, coefs(5)
         logical, parameter :: parabolic = .true.

         accuracy: if ((.not.use_table(T(1))).or.(.not.present(iT)))
     &                then

         if (present(mask)) then

            where (mask)
               where (tsw.gt.T(1) )
                 CpuR = (aL +bL*T(1) +cL*T(2) +dL*T(3) +eL*T(4))
               elsewhere
                 CpuR = (aH +bH*T(1) +cH*T(2) +dH*T(3) +eH*T(4))
               end where
            elsewhere
              CpuR = 0.e0_dp
            end where

         else ! .not.present(mask)

            where (tsw.gt.T(1) )
              CpuR  = (aL +bL*T(1) +cL*T(2) +dL*T(3) +eL*T(4))
            elsewhere
              CpuR  = (aH +bH*T(1) +cH*T(2) +dH*T(3) +eH*T(4))
            end where

         endif

         else

         if (present(mask)) then

          where (mask)
            CpuR = ( tab_CpsuR(:,iT-1)*(1.e0_dp-frac)
     &             + tab_CpsuR(:,iT  )*frac )
          elsewhere
            CpuR = 0.e0_dp
          end where

         else

            tabulation_order: if (parabolic) then

                coefs = interp4_coefs(frac)
                CpuR  = tab_CpsuR(:,iT-3)*coefs(1) +
     &                  tab_CpsuR(:,iT-2)*coefs(2) +
     &                  tab_CpsuR(:,iT-1)*coefs(3) +
     &                  tab_CpsuR(:,iT  )*coefs(4) +
     &                  tab_CpsuR(:,iT+1)*coefs(5)

            else

               CpuR = tab_CpsuR(:,iT-1)*(1.e0_dp-frac) +
     &                tab_CpsuR(:,iT)*frac

            endif tabulation_order


         endif


         endif accuracy

         end subroutine CpuRmol


c        ***************************************************************
c        ** Species internal energies in moles [J/mol]                **
c        ***************************************************************
         function int_energy(T,iT,frac) result(Umol)

         use speedchem,       only: ns
         use chemistry_setup, only: accurate_scthermo
         use SCthermodata,    only: tsw, aL, bL, cL, dL, eL, fL,
     &                                   aH, bH, cH, dH, eH, fH,
     &                                   tab_HsuRT, tab_CpsuR,
     &                                   interp4_coefs, use_table
         use universal_constants, only: R, u2, u3, u4, u5, u6

         implicit none

         real (dp)       , dimension(6), intent(in)  :: T
         integer, optional, intent(in)               :: iT
         real (dp)       , optional, intent(in)      :: frac
         real (dp)       , dimension(ns)             :: Umol
         real (dp)                                   :: coefs(5)
         logical, parameter :: parabolic = .true.


         if ((.not.use_table(T(1))).or.(.not.present(iT))) then

            where (tsw.gt.T(1) )
              Umol  = (aL - 1.e0_dp)*T(1) + u2*bL*T(2) + cL*T(3)*u3
     &                + dL*T(4)*u4 + u5*eL*T(4)*T(1) + fL
            elsewhere
              Umol  = (aH - 1.e0_dp)*T(1) + u2*bH*T(2) + cH*T(3)*u3
     &                + dH*T(4)*u4 + u5*eH*T(4)*T(1) + fH
            end where

            Umol = Umol * R

         else

            tabulation_order: if (parabolic) then

                coefs = interp4_coefs(frac)
                Umol  = tab_HsuRT(:,iT-3)*coefs(1) +
     &                  tab_HsuRT(:,iT-2)*coefs(2) +
     &                  tab_HsuRT(:,iT-1)*coefs(3) +
     &                  tab_HsuRT(:,iT  )*coefs(4) +
     &                  tab_HsuRT(:,iT+1)*coefs(5)


               Umol = (Umol - 1.e0_dp) * R * T(1)

            else

               Umol = tab_HsuRT(:,iT-1)*(1.e0_dp-frac) +
     &                tab_HsuRT(:,iT)*frac

               Umol = (Umol - 1.e0_dp) * R * T(1)

            endif tabulation_order



         endif


         end function int_energy

c        ***************************************************************
c        ** Species enthalpies in moles [J/mol]                       **
c        ***************************************************************
         function enthalpy(T,iT,frac) result(Hmol)

         use speedchem,       only: ns
         use chemistry_setup, only: accurate_scthermo
         use SCthermodata,    only: tsw, aL, bL, cL, dL, eL, fL,
     &                                   aH, bH, cH, dH, eH, fH,
     &                                   tab_HsuRT, tab_CpsuR,
     &                                   interp4_coefs, use_table
         use universal_constants, only: R, u2, u3, u4, u5, u6

         implicit none

         real (dp)       , dimension(6), intent(in)  :: T
         integer, optional, intent(in)               :: iT
         real (dp)       , optional, intent(in)      :: frac
         real (dp)       , dimension(ns)             :: Hmol
         real (dp)                                   :: coefs(5)
         logical, parameter :: parabolic = .true.


         if ((.not.use_table(T(1))).or.(.not.present(iT))) then

            where (tsw.gt.T(1) )
              Hmol  = aL*T(1) + u2*bL*T(2) + cL*T(3)*u3
     &                + dL*T(4)*u4 + u5*eL*T(4)*T(1) + fL
            elsewhere
              Hmol  = aH*T(1) + u2*bH*T(2) + cH*T(3)*u3
     &                + dH*T(4)*u4 + u5*eH*T(4)*T(1) + fH
            end where

            Hmol = Hmol * R

         else

            tabulation_order: if (parabolic) then

                coefs = interp4_coefs(frac)
                Hmol  = tab_HsuRT(:,iT-3)*coefs(1) +
     &                  tab_HsuRT(:,iT-2)*coefs(2) +
     &                  tab_HsuRT(:,iT-1)*coefs(3) +
     &                  tab_HsuRT(:,iT  )*coefs(4) +
     &                  tab_HsuRT(:,iT+1)*coefs(5)


               Hmol = Hmol * R * T(1)

            else

               Hmol = tab_HsuRT(:,iT-1)*(1.e0_dp-frac) +
     &                tab_HsuRT(:,iT)*frac

               Hmol = Hmol * R * T(1)

            endif tabulation_order



         endif


         end function enthalpy



c        ***************************************************************
c        ** Derivatives of specific heats at conV in moles [J/mol/K]  **
c        ***************************************************************
         function dCv_dT(T,iT,frac) result(dCvdT)

         use speedchem,       only: ns
         use chemistry_setup, only: accurate_scthermo
         use SCthermodata,    only: tsw, bL, cL, dL, eL,
     &                                   bH, cH, dH, eH,
     &                                   tab_dCvdT, interp4_coefs,
     &                                   use_table
         use universal_constants, only: R

         implicit none

         real (dp)       , dimension(6), intent(in)  :: T
         integer,          optional,     intent(in)  :: iT
         real (dp)       , optional,     intent(in)  :: frac
         real (dp)       , dimension(ns)             :: dCvdT
         real (dp)                                   :: coefs(5)

         logical, parameter :: parabolic = .true.


         if ((.not.use_table(T(1))).or.(.not.present(iT))) then

            where (tsw.gt.T(1) )
              dCvdT = R* (bL + 2*cL*T(1) + 3*dL*T(2) + 4*eL*T(3))
            elsewhere
              dCvdT = R* (bH + 2*cH*T(1) + 3*dH*T(2) + 4*eH*T(3))
            end where

!            dCvdT = dCvdT * R

         else

            coefs = interp4_coefs(frac)

            tabulation_order: if (parabolic) then

               dCvdT  = tab_dCvdT(:,iT-3)*coefs(1) +
     &                  tab_dCvdT(:,iT-2)*coefs(2) +
     &                  tab_dCvdT(:,iT-1)*coefs(3) +
     &                  tab_dCvdT(:,iT  )*coefs(4) +
     &                  tab_dCvdT(:,iT+1)*coefs(5)

            else

               dCvdT  = tab_dCvdT(:,iT-1)*(1.e0_dp-frac) +
     &                  tab_dCvdT(:,iT  )*frac

            endif tabulation_order

         endif


         end function dCv_dT



c        ***************************************************************
c        ** Species nondimensional gibbs potentials                   **
c        ***************************************************************
         function gibbs(T) result(G)

         use speedchem,       only: ns
         use SCthermodata,    only: tsw, aL, bL, cL, dL, eL, fL, gL,
     &                                   aH, bH, cH, dH, eH, fH, gH,
     &                                   gibbsL, gibbsH

         use universal_constants, only: u2, u6, u12, u20

         implicit none

         real (dp)       , dimension(6), intent(in)  :: T
         real (dp)       , dimension(ns)             :: G
         integer :: i

         if (allocated(gibbsL)) then

         do i = 1, ns
            if (tsw(i) > T(1)) then
              G(i) = gibbsL(1,i) + dot_product(gibbsL(2:7,i),T)
            else
              G(i) = gibbsH(1,i) + dot_product(gibbsH(2:7,i),T)
            endif
         end do

         else

            where (tsw.gt.T(1) )
              G = aL*(1.e0_dp-T(6)) -u2*bL*T(1) -u6*cL*T(2) -u12*dL*T(3)
     &            - u20*eL*T(4) + fL*T(5) - gL

            elsewhere
              G = aH*(1.e0_dp-T(6)) -u2*bH*T(1) -u6*cH*T(2) -u12*dH*T(3)
     &            - u20*eH*T(4) + fH*T(5) - gH
            end where

         endif

         end function gibbs

c        ***************************************************************
c        ** Species nondimensional gibbs potential derivatives with   **
c        ** respect to temperature [1/K]                              **
c        ***************************************************************
         function dgibbs_dT(T) result(dGdT)

         use speedchem,       only: ns
         use SCthermodata,    only: tsw, aL, bL, cL, dL, eL, fL, gL,
     &                                   aH, bH, cH, dH, eH, fH, gH

         use universal_constants, only: u2, u3, u4, u5

         implicit none

         real (dp)       , dimension(6), intent(in)  :: T
         real (dp)       , dimension(ns)             :: dGdT

            where (tsw.gt.T(1) )
c              dGdT = - ( gL + aL * T(6) + bL * T(1) + u2 * cL * T(2)
c     &               +   u3 * dL * T(3) + u4 * eL * T(4) )
              dGdt = - ( aL * T(5) + u2 * bL + u3 * cL * T(1)
     &               +   u4 * dL * T(2) + u5 * eL * T(3) + fL * T(5)**2)


            elsewhere
c              dGdT = - ( gH + aH * T(6) + bH * T(1) + u2 * cH * T(2)
c     &               +   u3 * dH * T(3) + u4 * eH * T(4) )
              dGdt = - ( aH * T(5) + u2 * bH + u3 * cH * T(1)
     &               +   u4 * dH * T(2) + u5 * eH * T(3) + fH * T(5)**2)


            end where


         end function dgibbs_dT

      end module SCspeciesthermo

c     ** Module for mixture-averaged thermodynamic properties *********

      module SCmixturethermo

      use working_precision, only: dp
      implicit none
      public

      real (dp)        :: MWm, Rm, cp, cv, h, e
!$OMP THREADPRIVATE(cp,cv)
      real (dp)       , target :: SCP, SCrho
!$OMP THREADPRIVATE(SCP,SCrho)


      contains

c        ***************************************************************
c        ** Species molar volumes multiplied by stoich coefficients   **
c        ***************************************************************

         subroutine molar_volumes

         use speedchem,        only: nr, ns, SCMW,
     &                               nudiff_pack,n_nudiff,i_nudiff
         use sparse_chemistry, only: nudiffT_molarv_sparse,
     &                               nudiffT_sparse
         use sparse_algebra,   only: sparse_allocated,
     &                               sparse_col_prod

         implicit none

         real (dp)       , dimension(:,:), allocatable :: nudiffT_molarv
         real (dp)                                     :: thouurho
         integer :: j, k

         thouurho = 1.0d+03/SCrho

c        Check for need to allocate nudiffT_molarv_sparse matrix
         if (.not.sparse_allocated(nudiffT_molarv_sparse)) then
            nudiffT_molarv_sparse = nudiffT_sparse
         endif

c        If sparse molar volumes array has already been initialised
           do k = 1, ns
           if (nudiffT_molarv_sparse%IA(k+1) >=
     &         nudiffT_molarv_sparse%IA(k)       )
     &         nudiffT_molarv_sparse%A( nudiffT_molarv_sparse%IA(k) :
     &                              nudiffT_molarv_sparse%IA(k+1)-1 )
     &          = nudiffT_sparse%A( nudiffT_molarv_sparse%IA(k) :
     &                              nudiffT_molarv_sparse%IA(k+1)-1 )
     &          * thouurho * SCMW(k)
           end do


         end subroutine molar_volumes

c        ***************************************************************
c        ** Mixture average constant volume specific heat [J/kg/K]    **
c        ***************************************************************
         function cvmas(Ta,Y,iT,frac,nonzeroes) result(cvmixture)

         use speedchem,           only: ns, uMW
         use SCthermodata,        only: use_table
         use SCspeciesthermo,     only: CvuRmol
         use universal_constants, only: R, kilo

         implicit none

         real (dp)       , dimension(6), intent(in)  :: Ta
         real (dp)       , dimension(ns),intent(in)  :: Y
         integer, optional, intent(in)               :: iT
         real (dp)       , optional, intent(in)      :: frac
         logical,          dimension(ns), optional   :: nonzeroes

         real (dp)       , dimension(ns)             :: CvuR
         real (dp)                                   :: cvmixture


c        Retrieve non dimensional constant volume specific heats [-]
         if (use_table(Ta(1)).and.present(iT)) then
             call CvuRmol(Ta,CvuR,iT,frac)!,nonzeroes)
         else
             call CvuRmol(Ta,CvuR)!,nonzeroes)
         endif

c        Mixture averaged mass value [J/kg/K]
         if (present(nonzeroes)) then
           cvmixture = kilo * R * sum(Y * CvuR * uMW, nonzeroes)
         else
           cvmixture = kilo * R * sum(Y * CvuR * uMW)
         endif

         end function cvmas

c        ***************************************************************
c        ** Mixture average constant pressure specific heat [J/kg/K]  **
c        ***************************************************************
         function cpmas(Ta,Y,iT,frac,nonzeroes) result(cpmixture)

         use speedchem,           only: ns, uMW
         use SCthermodata,        only: use_table
         use SCspeciesthermo,     only: CpuRmol
         use universal_constants, only: R, kilo

         implicit none

         real (dp)       , dimension(6), intent(in)  :: Ta
         real (dp)       , dimension(ns),intent(in)  :: Y
         integer, optional, intent(in)               :: iT
         real (dp)       , optional, intent(in)      :: frac
         logical,          dimension(ns), optional   :: nonzeroes

         real (dp)       , dimension(ns)             :: CpuR
         real (dp)                                   :: cpmixture


c        Retrieve non dimensional constant volume specific heats [-]
         if (use_table(Ta(1)).and.present(iT)) then
             call CpuRmol(Ta,CpuR,iT,frac)!,nonzeroes)
         else
             call CpuRmol(Ta,CpuR)!,nonzeroes)
         endif

c        Mixture averaged mass value [J/kg/K]
         if (present(nonzeroes)) then
           cpmixture = kilo * R * sum(Y * CpuR * uMW, nonzeroes)
         else
           cpmixture = kilo * R * sum(Y * CpuR * uMW)
         endif

         end function cpmas

c        ***************************************************************
c        ** Mixture average specific internal energy [J/kg]           **
c        ***************************************************************
         function umas(Ta,Y,iT,frac,nonzeroes) result(umixture)

         use speedchem,           only: ns, uMW
         use SCthermodata,        only: use_table
         use SCspeciesthermo,     only: int_energy
         use universal_constants, only: kilo

         implicit none

         real (dp)       , dimension(6), intent(in)  :: Ta
         real (dp)       , dimension(ns),intent(in)  :: Y
         integer,          optional,     intent(in)  :: iT
         real (dp)       , optional,     intent(in)  :: frac
         logical,          dimension(ns), optional   :: nonzeroes

         real (dp)       , dimension(ns)             :: u
         real (dp)                                   :: umixture


c        Retrieve molar internal energies of the species [J/mol]
         if (use_table(Ta(1)).and.present(iT)) then
             u = int_energy(Ta,iT,frac)
         else
             u = int_energy(Ta)
         endif

c        Mixture averaged mass value [J/kg/K]
         if (present(nonzeroes)) then
           umixture = kilo * sum(Y * u * uMW, nonzeroes)
         else
           umixture = kilo * sum(Y * u * uMW)
         endif

         end function umas

c        ***************************************************************
c        ** Mixture average gamma = cp/cv                             **
c        ***************************************************************
         function gammamix(Ta,Y,iT,frac) result(gammixture)

         use speedchem,           only: ns, uMW
         use SCthermodata,        only: use_table
         use SCspeciesthermo,     only: CpuRmol
         use universal_constants, only: R

         implicit none

         real (dp)       , dimension(6), intent(in)  :: Ta
         real (dp)       , dimension(ns),intent(in)  :: Y
         integer, optional, intent(in)               :: iT
         real (dp)       , optional, intent(in)      :: frac

         real (dp)       , dimension(ns)             :: CpuR
         real (dp)                                   :: gammixture,
     &                                                  cpurmixture


c        Retrieve non dimensional constant volume specific heats [-]
         if (use_table(Ta(1)).and.present(iT)) then
             call CpuRmol(Ta,CpuR,iT,frac)
         else
             call CpuRmol(Ta,CpuR)
         endif

c        average mixture Cp/R [-]
         cpurmixture = sum(Y * CpuR * uMW)

c        Mixture averaged gamma value [-]
         gammixture = cpurmixture / ( cpurmixture - sum(Y * uMW) )

         end function gammamix

c        ***************************************************************
c        **                                                           **
c        **     Computing concentration based on mass fractions       **
c        **                                                           **
c        ***************************************************************
         subroutine massfr_to_conc (Y,T,C)
         use speedchem,           only: ns, uMW
         use universal_constants, only: uR, micro, one
         implicit none

         real (dp)       , intent(in) :: T
         real (dp)       , dimension(ns), intent(in)  :: Y
         real (dp)       , dimension(ns), intent(out) :: C
         real (dp)                                    :: uT

c        **************************************************************

c        Reciprocal of temperature [1/K]
         uT = one/T

c        Concentration in mole units [mol/cm3]
         C = micro * SCP * uR * uT * (Y*uMW) / sum(Y*uMW)

         end subroutine massfr_to_conc



c        **************************************************************
c        **         Computing average mixture density [kg/m^3]       **
c        **************************************************************
         subroutine rhoY (Y, T)
         use speedchem,           only : nr, ns, uMW
         use universal_constants, only : uR, one, kilo, milli

         implicit none
         real (dp)       , dimension(ns), intent(in) :: Y
         real (dp)       , dimension(ns)             :: C
         real (dp)       , intent(in)                :: T
         real (dp)                                   :: uT

!        Computing species concentrations
!         call massfr_to_conc (Y,T,C)
!         SCrho = kilo * sum(SCMW*C)

c        Reciprocal of temperature [1/K]
         uT = one/T

c        Density value [kg/m3] !NB: pressure SCP must be initialised!
         SCrho = milli * SCP * uR * uT * sum(Y / sum(Y*uMW))

         end subroutine rhoY

c        **************************************************************
c        **          Computing average reactor pressure [Pa]         **
c        **************************************************************
         function pressurerhoT(T,Y) result(pressure)
         use speedchem,           only: ns,nr, SCMW, uMW
         use universal_constants, only: R, kilo


         implicit none
         real (dp)       , intent(in) :: T
         real (dp)       , dimension(ns), intent(in) :: Y
         real (dp)        :: pressure

         pressure = kilo * SCrho * R * T * sum(Y*uMW)

         end function pressurerhoT

c        ***************************************************************
c        ** Computing mixture EQUIVALENCE RATIO                       **
c        **                                                           **
c        ** INPUT DATA                                                **
c        ** X = array of molar fractions [-]                          **
c        ***************************************************************
         function equivalence_ratio(X) result(phi)
         use speedchem,           only: ns, nel
         use sparse_algebra,      only: sparse_matmulT
         use sparse_chemistry,    only: EM_sp
         use scthermodata,        only: get_el
         use universal_constants, only: two, u2

         real (dp), dimension(ns) , intent(in)  :: X
         real (dp)                              :: phi
         real (dp), dimension(nel)              :: Xel

         ! Sum total molar fractions of elements in the
         ! mixture
         Xel = sparse_matmulT(EM_sp,X)

         ! Compute global equivalence ratio
         phi = (two*Xel(get_el('C')) + u2*Xel(get_el('H')))
     &       / Xel(get_el('O'))


         end function equivalence_ratio


      end module SCmixturethermo


c     ** Module for Troe Parameters computation ***********************

      module troepar

      use working_precision, only: dp
      implicit none
      public

c        THIRD-BODY REACTIONS *****************************************

c        INDICES in main reaction indexing: (j = 1, nr) ***************

c        nTBALL  = total number of all three body reactions
c        iTBALL  = indices of all three body reactions
         integer                                       :: ntbALL
         integer, dimension(:), allocatable            :: itbALL
         logical, dimension(:), allocatable            :: ltbALL

c        nTBTROE = number of fall-off reactions (ie TROE or LINDEMANN)
c        iTBTROE = indices of fall-off reactions
         integer                                       :: ntbFALL
         integer, dimension(:), allocatable            :: itbFALL
         logical, dimension(:), allocatable            :: ltbFALL

c        nTBTROE = total number of Troe reactions
c        iTBTROE = indices of Troe reactions
         integer                                       :: ntbTROE
         integer, dimension(:), allocatable            :: itbTROE
         logical, dimension(:), allocatable            :: ltbTROE

c        Sub-Indexes of reactions with all 4 troe parameters present
c        (This indexes are positions in the itbTROE array, not in the
c        whole reactions array!)
         integer                                       :: nTROE4
         integer, dimension(:), allocatable            :: iTROE4
         logical, dimension(:), allocatable            :: lTROE4

c        nTBLIND = total number of Lindemann reactions
c        iTBLIND = indices of Lindemann reactions
         integer                                       :: ntbLIND
         integer, dimension(:), allocatable            :: itbLIND
         logical, dimension(:), allocatable            :: ltbLIND

c        nTBSIMP = total number of simple three body reactions
c        iTBSIMP = indices of simple three body reactions
         integer                                       :: ntbSIMP
         integer, dimension(:), allocatable            :: itbSIMP
         logical, dimension(:), allocatable            :: ltbSIMP

c        SUB-INDICES in TBALL reaction indexing: (j = 1, nTBALL) ******

         integer, dimension(:), allocatable :: iTROEitbALL,
     &                                         iLINDitbALL,
     &                                         iSIMPitbALL,
     &                                         iFALLitbALL,
     &                                         iTROEitbFALL,
     &                                         iLINDitbFALL


      real (dp)        :: a1,a2
      integer, dimension(:), allocatable :: todotroe,zeroT2
      real (dp)       , dimension(:), allocatable :: aT2, uT1T2, T2T2,
     &                                               uT3T2

!      real (dp)       , dimension(:), allocatable :: expuT3T2,
!     &                                               expuT1T2

c     Variables for tabulated Troe parameters
      real (dp)       , dimension(:,:), allocatable :: tab_troefactor,
     &                                                 tab_log10Fcent

c     Table accuracies
c     1] Pr values are in the [0,1] range
      real (dp)        :: Pr_table_accuracy = 1.e-2_dp


      integer :: Pr_intervals, troeTemp_intervals

      contains
c     *****************************************************************

c        **************************************************************
c        **         Initialise third-body reaction indices           **
c        **************************************************************
         subroutine init_thirdbody_indices

          use speedchem, only: ns, nr, species, reactions,
     &                         LOW, HIGH, TROE, third_body, T2TROE
          use sparse_chemistry, only: third_body_sp,
     &                                sparse_internal_sum
          use find_mod
          implicit none

          integer :: i, ir
          real (dp)       , dimension(nr) :: sumtb2
          character(len=*), parameter :: fmt = "(A10,100(1x,i4))",
     &       fmt_hed = "(' Non-Arrhenius reactions stats ',/,"//
     &                 "1x,41('-'))",
     &       fmt_all = "(' Three-body and pressure-dependent : ',I5)",
     &       fmt_tro = "(' Troe pressure-dependent           : ',I5)",
     &       fmt_thr = "(' simple three-body                 : ',I5)",
     &       fmt_lin = "(' Lindemann pressure-dependent      : ',I5)",
     &       fmt_clo = "(1x,41('-'))"
          integer, dimension(:), allocatable :: iALLtmp



          allocate( ltbLIND(nr), ltbTROE(nr), ltbSIMP(nr),
     &              ltbALL(nr),  ltbFALL(nr))

c         All fall-off (Lindemann or Troe) reactions
          ltbFALL  = (LOW .or. HIGH)
          ntbFALL  = count (ltbFALL)
          allocate (itbFALL(ntbFALL))
          itbFALL  = pack  (reactions, ltbFALL)


c         Lindemann reactions
          ltbLIND  = (LOW .or. HIGH) .and. (.not.TROE)
          ntbLIND  = count (ltbLIND)
          allocate (itbLIND(ntbLIND))
          itbLIND  = pack  (reactions, ltbLIND)


c         Troe reactions
          ltbTROE  = (LOW .or. HIGH) .and. TROE
          ntbTROE  = count (ltbTROE)
          allocate (itbTROE(ntbTROE))
          itbTROE  = pack  (reactions, ltbTROE)


c         Troe reactions with three parameters
          allocate (lTROE4(ntbTROE))
          lTROE4   = T2TROE(itbTROE) /= 0.e0_dp
          nTROE4   = count (lTROE4)
          allocate (iTROE4(nTROE4))
          iTROE4   = pack  ([(i,i=1,ntbTROE)], lTROE4)


c         Simple third-body reactions
          call sparse_internal_sum(third_body_sp,sumtb2, dim=2)
          ltbSIMP  = (.not.(HIGH.or.LOW)) .and. sumtb2>0.e0_dp
          ntbSIMP  = count (ltbSIMP)
          allocate (itbSIMP(ntbSIMP))
          itbSIMP  = pack  (reactions, ltbSIMP)


c         All three body reactions
          ltbALL   = sumtb2>0.e0_dp
          ntbALL   = count (ltbALL)
          allocate (itbALL(ntbALL))
          itbALL   = pack  (reactions, ltbALL)


c         SUB-INDICES referring to position in itbALL array
          allocate (iSIMPitbALL (ntbSIMP),
     &              iLINDitbALL (ntbLIND),
     &              iTROEitbALL (ntbTROE),
     &              iTROEitbFALL(ntbTROE),
     &              iLINDitbFALL(ntbLIND),
     &              iFALLitbALL (ntbFALL))

          do i = 1, ntbSIMP
             call find_indices(itbALL == itbSIMP(i))
             iSIMPitbALL(i) = indices(1)
          end do

          do i = 1, ntbLIND
             call find_indices(itbALL == itbLIND(i))
             iLINDitbALL(i) = indices(1)

             call find_indices(itbFALL == itbLIND(i))
             iLINDitbFALL(i) = indices(1)
          end do

          do i = 1, ntbTROE
             call find_indices(itbALL == itbTROE(i))
             iTROEitbALL(i) = indices(1)

             call find_indices(itbFALL == itbTROE(i))
             iTROEitbFALL(i) = indices(1)

          end do

          do i = 1, ntbFALL
             call find_indices(itbALL == itbFALL(i))
             iFALLitbALL(i)  = indices(1)
          end do

!         printout: if (verbose) then

            write(*,       *)
            write(*, fmt_hed)
            write(*, fmt_all)ntbALL
            write(*, fmt_thr)ntbSIMP
            write(*, fmt_tro)ntbTROE
            write(*, fmt_lin)ntbLIND
            write(*, fmt_clo)
            write(*,       *)


         end subroutine init_thirdbody_indices




      subroutine tabulate_troepars
c     *****************************************************************
c     **                                                             **
c     **                     SpeedCHEM FORTRAN                       **
c     **                                                             **
c     **     Tabulate functions for computing TROE molecularity      **
c     **                                                             **
c     **                                                             **
c     **   Author:      Federico Perini                              **
c     **   Last update: sunday, 30/10/2011                           **
c     **                                                             **
c     *****************************************************************
      use chemistry_setup,     only: Temp_table_accuracy, Temp_LOlim,
     &                               Temp_HIlim,tab_nsteps
      use universal_constants, only: one

      implicit none

      integer                       :: i, j
      real (dp)                     :: T, unosuT
      real (dp), dimension(ntbTROE) :: Fcent

c     String formats
      character(len=*), parameter   ::
     &      fmt_end = "(' Troe parameters tabulated. ')"

c     Compute number of subdivision in temperature and Pressure ratio
c     ranges
      Pr_intervals  = int(one/Pr_table_accuracy)

      tab_nsteps    = int( ( Temp_HIlim - Temp_LOlim )
     &                    / Temp_table_accuracy)

      allocate(tab_log10Fcent(ntbTROE,-1:tab_nsteps))

c     ** TROE *********************************************************
c     Computing reaction rate constants according to Troe's formulation
      do i=1,tab_nsteps + 2

         T = Temp_LOlim + real(i-2, dp) * Temp_table_accuracy
         unosuT = one/T

c     ** Computing troe centering factor
         Fcent = (one - aT2) * exp(-T    * uT3T2 )
     &         +        aT2  * exp(-T    * uT1T2 )
     &         +               exp(-T2T2 * unosuT)

         Fcent(zeroT2) = Fcent(zeroT2) - one

         tab_log10Fcent(:,i-2) = log10(Fcent)

      end do

      write(*,fmt_end)

      end subroutine tabulate_troepars

c        ***************************************************************
c        ** Compute enhancement factor for TROE reactions             **
c        ***************************************************************
         function troe_logfac(Ta,iT,frac) result(log10Fcent)

         use speedchem,       only: ns
         use chemistry_setup, only: accurate_scthermo

         use SCthermodata,    only: tsw, bL, cL, dL, eL,
     &                                   bH, cH, dH, eH,
     &                                   tab_dCvdT, interp4_coefs,
     &                                   use_table

         implicit none

         real (dp)       , dimension(6), intent(in)  :: Ta
         integer,          optional,     intent(in)  :: iT
         real (dp)       , optional,     intent(in)  :: frac
         real (dp)       , dimension(ntbTROE)        :: log10Fcent,
     &                                                  Fcent
         real (dp)                                   :: coefs(5)

         integer :: j
         logical, parameter :: parabolic = .true.

         if ((.not.use_table(Ta(1))).or.(.not.present(iT))) then

             Fcent = (1.e0_dp - aT2) * exp(-uT3T2 * Ta(1))
     &                     + aT2  * exp(-uT1T2 * Ta(1))

             Fcent(iTROE4) = Fcent(iTROE4) + exp(-T2T2(iTROE4) * Ta(5))

             log10Fcent    = log10(Fcent)

         else

            coefs = interp4_coefs(frac)

            tabulation_order: if (parabolic) then

                log10Fcent  = tab_log10Fcent(:,iT-3)*coefs(1) +
     &                        tab_log10Fcent(:,iT-2)*coefs(2) +
     &                        tab_log10Fcent(:,iT-1)*coefs(3) +
     &                        tab_log10Fcent(:,iT+0)*coefs(4) +
     &                        tab_log10Fcent(:,iT+1)*coefs(5)


             else

                log10Fcent = tab_log10Fcent (:,iT-1)*(1.e0_dp-frac) +
     &                       tab_log10Fcent (:,iT  )*frac

             endif tabulation_order
         endif


         end function troe_logfac

      end module troepar

c     ** Module for reaction rates computation ************************
      module reacpar

      use working_precision, only: dp
      implicit none
      public


c        REVERSIBLE REACTIONS *****************************************
c        nTREV = total number of reversible reactions
c        iTREV = indices of reversible reactions
         integer                                       :: nTREV
         integer, dimension(:), allocatable            :: iTREV

c        nXREV = total number of reversible reaction with explicit
c                Arrhenius coefficients for reverse reaction rate
c        iXREV = corresponding reaction indices
         integer                                       :: nXREV
         integer, dimension(:), allocatable            :: iXREV

c        nEQREV = total number of reversible reaction with equilibrium-
c                based reverse reaction rates
c        iEQREV = corresponding reaction indices
         integer                                       :: nEQREV
         integer, dimension(:), allocatable            :: iEQREV


c        THIRD-BODY REACTIONS *****************************************
c        Third body molecularity beta reaction coefficients
c        in packed form
         real (dp)       , dimension(:,:), allocatable :: tb_beta_pack
         integer,          dimension(:,:), allocatable :: is_beta_pack
         integer,          dimension(:),   allocatable :: n_tb_beta


         integer, dimension(:), allocatable :: Arrhreac,
     &                                         Lindreac,
     &                                         Troereac,
     &                                         Revreac

c        Storage of equilibrium constants
         real (dp)       , parameter     :: uKc_RTOL=1e-15_dp
         real (dp)       , dimension(:), allocatable   :: store_uKc
         real (dp)                                     :: store_uKc_T


      contains

c        ***************************************************************
c        ** Allocate and compute packed reaction indices              **
c        ***************************************************************

         subroutine init_reaction_indices

         use speedchem, only : ns, nr, nTHREE, iTHREE, third_body_beta,
     &                         AREV,bREV,EREV,ARi,bRi,ERi,reversibile,
     &                         REV, reactions, species, uMW

         use sparse_chemistry, only: nudiff_EQREV_sp, tb_beta_sp,
     &                               tbb_uMW_sp, nudiff_sparse,
     &                               third_body_sp
         use sparse_algebra, only:sparse, dense_to_sparse, sparse_value,
     &                            sparse_row_prod, sparse_internal_count

         implicit none

         real (dp)       , dimension(:,:), allocatable :: nudiff_EQREV,
     &                                                    tbb_T_pack
         integer :: i,j, ir, nmax, itmp
         integer, dimension(nTHREE) :: count_tb


!        Count number of species involved in three body reactions
         call sparse_internal_count(tb_beta_sp, count_tb, dim=2)

c        Allocate number of species with beta /= 0 in each tb reaction
         nmax = maxval(count_tb)
         allocate(n_tb_beta   (nTHREE),
     &            tb_beta_pack(nmax, nTHREE), is_beta_pack(nmax,nTHREE))



         do i = 1, nTHREE

            ir = iTHREE(i)

            n_tb_beta(i) = count_tb(i)

            if (n_tb_beta(i) > 0) then

               tb_beta_pack(1:n_tb_beta(i),i) =
     &         tb_beta_sp% A(tb_beta_sp%IA(i):tb_beta_sp%IA(i+1)-1)

               is_beta_pack(1:n_tb_beta(i),i) =
     &         tb_beta_sp%JA(tb_beta_sp%IA(i):tb_beta_sp%IA(i+1)-1)

            endif

         end do

c        Initialise third body coefficients divided by molecular weights
         tbb_uMW_sp =  sparse_row_prod(tb_beta_sp, uMW)



c        ******************* REVERSIBLE REACTIONS *********************

c        General reversible reactions
         nTREV = count(reversibile)

         reverse_present: if (nTREV > 0) then
            allocate(iTREV(nTREV))
            iTREV = pack(reactions, mask = reversibile)
         endif reverse_present


c        Explicit reverse reaction rates
         nXREV = count(REV)

         explicit_reverse_present: if (nXREV > 0) then
            allocate(iXREV(nXREV),ARi(nXREV),bRi(nXREV),ERi(nXREV))
            iXREV = pack(reactions, REV)
            ARi   = pack(AREV,      REV)
            bRi   = pack(bREV,      REV)
            ERi   = pack(EREV,      REV)
         endif explicit_reverse_present

c        Equilibrium-based reverse reaction rates
         nEQREV = count(reversibile.and.(.not.REV))

         equilibrium_reverse_present: if (nEQREV > 0) then
            allocate(iEQREV(nEQREV))

c           Reactions for which an equilibrium-derived reverse
c           reaction rate constant is present are the reversible ones
c           for which no explicit reverse coefficients have been
c           provided
            iEQREV = pack(reactions, reversibile.and.(.not.REV))

c           Equilibrium reactions coefficients in sparse format
            if (nEQREV > 0) then
            allocate(nudiff_EQREV(nEQREV,ns))
            do j = 1, ns
               do i = 1, nEQREV
                  nudiff_EQREV(i,j) = sparse_value(nudiff_sparse,
     &                                             iEQREV(i),j)
               end do
            end do
            call dense_to_sparse(nudiff_EQREV,nudiff_EQREV_sp)
            endif

         endif equilibrium_reverse_present


         end subroutine init_reaction_indices



c        ***************************************************************
c        ** Reciprocal of equilibrium constants in concentration units**
c        ***************************************************************
         function uequilC(T,iT,frac) result(uKc)

         use chemistry_setup,  only: accurate_scthermo
         use speedchem,        only: nr, ns, isumnudiff
         use SCspeciesthermo,  only: gibbs
         use SCthermodata,     only: tab_uKc, interp4_coefs, use_table
         use sparse_chemistry, only: nudiff_EQREV_sp
         use universal_constants, only: R, uPatm, u2, u3, u6
         use sparse_definitions
         implicit none

         real (dp)       , dimension(6), intent(in)  :: T
         integer, optional, intent(in)               :: iT
         real (dp)       , optional, intent(in)      :: frac
         real (dp)       , dimension(nEQREV)         :: dG0, uKp, uKc
         real (dp)       , dimension(ns)             :: G
         real (dp)                                   :: ufaceq, coefs(5)
         integer                                     :: i,ir

         logical, parameter :: parabolic = .true.

         if ((.not.use_table(T(1))).or.(.not.present(iT))) then

c          Nondimensional gibbs free energy [-]
           G = gibbs(T)

c          Sparse multiplication      dG0 = matmul(nudiff,gibbs)
           dG0 = nudiff_EQREV_sp * G

c          ************************************************************

c          Computing concentration-based equilibrium constant
           ufaceq = 1.0d+6 *  R *  T(1) * uPatm

c          Computing reciprocal of pressure-based equilibrium constant

            uKc = exp(+dG0) * ufaceq**isumnudiff(iEQREV)

         else ! .not.accurate_scthermo

            tabulation_order: if (parabolic) then

               coefs = interp4_coefs(frac)
               uKc  = tab_uKc(:,iT-3)*coefs(1) +
     &                tab_uKc(:,iT-2)*coefs(2) +
     &                tab_uKc(:,iT-1)*coefs(3) +
     &                tab_uKc(:,iT  )*coefs(4) +
     &                tab_uKc(:,iT+1)*coefs(5)

            else

               uKc = tab_uKc(:,iT-1)*(1.e0_dp-frac) + tab_uKc(:,iT)*frac

            endif tabulation_order

         endif ! accurate_scthermo for equilibrium constant


         end function uequilC

c        ***************************************************************
c        ** Reciprocal of equilibrium constants in concentration units
c        ** and derivative of the equilibrium constants with respect
c        ** to temperature [(Kc_units)/K]
c        ***************************************************************
         function uequilC_and_derivative(T,iT,frac) result (uKc_dKcdT)

         use chemistry_setup,  only: accurate_scthermo
         use speedchem,        only: nr, ns, isumnudiff
         use SCspeciesthermo,  only: gibbs, dgibbs_dT
         use SCthermodata,     only: tab_uKc, tab_dKcdT, interp4_coefs,
     &                               use_table
         use sparse_chemistry, only: nudiff_EQREV_sp
         use universal_constants, only: R, uR, uPatm, u2, u3, u6
         use sparse_definitions
         implicit none

         real (dp)       , dimension(6), intent(in)  :: T
         integer, optional, intent(in)               :: iT
         real (dp)       , optional, intent(in)      :: frac
         real (dp)       , dimension(nEQREV,2)       :: uKc_dKcdT
         real (dp)       , dimension(nEQREV)         :: dG0, dDG0dT,
     &                                                  dKcdT
         real (dp)       , dimension(ns)             :: G, dGdT
         real (dp)                                   :: ufaceq, coefs(5)

         logical, parameter :: parabolic = .true.

         if ((.not.use_table(T(1))).or.(.not.present(iT))) then

c           Nondimensional gibbs free energy [-] and derivative [1/K]
            G    = gibbs(T)
            dGdT = dgibbs_dT(T)

c           Sparse multiplication      dG0 = matmul(nudiff,gibbs)
            dG0 = nudiff_EQREV_sp * G

c           Sparse multiplication   dDG0dT = matmul(nudiff,dgibbs_dT)
            dDG0dT = nudiff_EQREV_sp * dGdT


c           ***********************************************************

c           Computing concentration-based equilibrium constant
            ufaceq = 1.0d+6 *  R *  T(1) * uPatm

c           Computing reciprocal of pressure-based equilibrium constant
            uKc_dKcdT(:,1) = exp(+dG0) * ufaceq**isumnudiff(iEQREV)

c           Computing derivative of (not reciprocal) equilibrium const
            uKc_dKcdT(:,2) = 1.e0_dp/uKc_dKcdT(:,1)
     &                     * ( -isumnudiff(iEQREV)*T(5) - dDG0dT )


         else ! .not.accurate_scthermo

          tabulation_order: if (parabolic) then

                coefs = interp4_coefs(frac)

                uKc_dKcdT(:,1)  = tab_uKc(:,iT-3)*coefs(1) +
     &                            tab_uKc(:,iT-2)*coefs(2) +
     &                            tab_uKc(:,iT-1)*coefs(3) +
     &                            tab_uKc(:,iT+0)*coefs(4) +
     &                            tab_uKc(:,iT+1)*coefs(5)

                uKc_dKcdT(:,2)  = tab_dKcdT(:,iT-3)*coefs(1) +
     &                            tab_dKcdT(:,iT-2)*coefs(2) +
     &                            tab_dKcdT(:,iT-1)*coefs(3) +
     &                            tab_dKcdT(:,iT+0)*coefs(4) +
     &                            tab_dKcdT(:,iT+1)*coefs(5)

          else

            uKc_dKcdT(:,1) = tab_uKc(:,iT-1)*(1.e0_dp-frac)
     &                     + tab_uKc(:,iT)  *      frac

            uKc_dKcdT(:,2) = tab_dKcdT(:,iT-1)*(1.e0_dp-frac)
     &                     + tab_dKcdT(:,iT)  *      frac

          endif tabulation_order


         endif ! accurate_scthermo for equilibrium constant


         end function uequilC_and_derivative

c     *****************************************************************
c     **                                                             **
c     **                     SpeedCHEM FORTRAN                       **
c     **                                                             **
c     **    Tabulation of reaction equilibrium constants with T      **
c     **                                                             **
c     **                                                             **
c     **   Author:      Federico Perini                              **
c     **   Last update: tuesday, 15/11/2011                          **
c     **                                                             **
c     *****************************************************************


      subroutine tabulate_equilibrium

      use chemistry_setup, only: Temp_table_accuracy, Temp_LOlim,
     &                       Temp_HIlim, tab_nsteps

      use speedchem, only : ns,nr,isumnudiff

      use scthermodata, only: tab_uKC, tab_uKp, tab_dKcdT,
     &                        tsw, aL, bL, cL, dL, eL, fL, gL,
     &                             aH, bH, cH, dH, eH, fH, gH,
     &                        temperature_array, use_table
      use sparse_definitions
      use sparse_chemistry, only: nudiff_EQREV_sp
      use SCspeciesthermo,  only: gibbs, dgibbs_dT
      use universal_constants, only: u2,u3,u4,u5, R, uR, Patm, uPatm

      implicit none
      integer :: i, j, k, nsteps
      real (dp)        :: T, T2, T3, T4, unosuT, lnT
      real (dp)       , dimension(ns)    :: G, dGdT
      real (dp)       , dimension(nEQREV)::dG0, uKp, uKc, dDG0dT, dKcdT
      real (dp)                   :: faceq, ufaceq
      real (dp)       , dimension(6)     :: Ta

c     Only perform tabulation if there are equilibrium-based reversible
c     reactions
      tabulate: if (nEQREV > 0) then

      tab_nsteps = int( ( Temp_HIlim - Temp_LOlim )
     &                   / Temp_table_accuracy)

      allocate(tab_uKc  (nEQREV,-1:tab_nsteps),
     &         tab_uKp  (nEQREV,-1:tab_nsteps),
     &         tab_dKcdT(nEQREV,-1:tab_nsteps) )

      do i=1,tab_nsteps + 2



         T = Temp_LOlim + real(i-2, dp) * Temp_table_accuracy

         Ta = temperature_array(T)

         T2     = Ta(2)
         T3     = Ta(3)
         T4     = Ta(4)
         unosuT = Ta(5)
         lnT    = Ta(6)

c        Compute nondimensional gibbs potential T [-] and
c        its temperature derivative [1/K]
         G    = gibbs(Ta)
         dGdt = dgibbs_dT(Ta)

         dG0    = (nudiff_EQREV_sp * G   )
         dDG0dT = (nudiff_EQREV_sp * dGdT)

c        **************************************************************
c        Computing reciprocal of pressure-based equilibrium constant
         uKp = exp(+dG0)! * uR * unosuT)

c        Computing concentration-based equilibrium constant
         ufaceq = 1.0d+6 *  R *  T * uPatm

         uKc = uKp * ufaceq**isumnudiff(iEQREV)
         dKcdT = 1.e0_dp/uKc * ( -isumnudiff(iEQREV)*Ta(5) - dDG0dT )

         tab_uKp(:,i-2) = uKp
         tab_uKc(:,i-2) = uKc
         tab_dKcdT(:,i-2) = dKcdT

      end do

      write(*,900)

      else

      write(*,905)

      end if tabulate

 900  format(' equilibrium data tabulated.')
 905  format(' no equilibrium parameters needed.')

      end subroutine tabulate_equilibrium






      end module reacpar



c     ** Module for computing reaction progress rates *****************
      module kinetics_mod

      use working_precision, only: dp
      implicit none
      public

      integer, dimension(:), allocatable :: ir,jr,ijr,ip,jp,ijp
      integer, dimension(:,:), allocatable :: i1r, i1p
      integer, dimension(:)  , allocatable :: i2r, i2p
      integer, dimension(:), allocatable :: i2D1r, i2D2r, i2D1p, i2D2p
      real (dp)       , dimension(:,:), allocatable :: unitr, unitp
      integer, dimension(:,:), allocatable :: indice_r, indice_p
      real (dp)       , dimension(:), allocatable :: v_stoich_r,
     &                                               v_stoich_p
      integer,          dimension(:), allocatable :: iv_stoich_r,
     &                                               iv_stoich_p
      integer, dimension (:), allocatable ::         v_stoich_r1,
     &                                               v_stoich_r2,
     &                                               v_stoich_ro,
     &                                               v_stoich_p1,
     &                                               v_stoich_p2,
     &                                               v_stoich_po
      integer :: num_vro, num_vpo, num_vr2, num_vp2

c     Variables for tabulated kinetics parameters
      real (dp)       , dimension(:,:), allocatable, target ::
     &                                               tab_k0, tab_kinf,
     &                                               tab_dkinfdt,
     &                                               tab_dk0dt,
     &                                               tab_Xkb, tab_Xdkbdt

      real (dp)       , dimension(:,:), allocatable :: tab_kinfT

      integer, dimension(:), allocatable :: iA0

c     Arrays for saved data for the Jacobian matrix
      real (dp)       , dimension(:), allocatable :: save_k0, save_kinf
      logical                                     :: lsavek = .false.
!$OMP THREADPRIVATE(save_k0, save_kinf, lsavek)


c     Array for storing species rates of change
      real (dp)       ,               parameter   :: T_RTOL = 1.e-15_dp
      real (dp)       , dimension(:), allocatable :: store_dwdt
      real (dp)         :: store_dwdt_T = 0.e0_dp

c     Array for storing reaction rate constants
      real (dp)       , dimension(:), allocatable :: store_kinf
      real (dp)                                   :: store_kinf_T

c     Array for storing mass action productories
      real (dp)       , dimension(:), allocatable, target :: store_qf,
     &                                                       store_qb
      real (dp)       , dimension(:), allocatable :: store_C


      contains

c     *****************************************************************
c     **                                                             **
c     **                     SpeedCHEM FORTRAN                       **
c     **                                                             **
c     **         Compute ordinary forward reaction rates             **
c     **                                                             **
c     **                                                             **
c     **   Author:      Federico Perini                              **
c     **   Last update: sunday, 03/12/2011                           **
c     **                                                             **
c     *****************************************************************
      subroutine reaction_rates(Ta, k0, kinf, iT, frac)

      use chemistry_setup, only: accurate_scthermo,
     &                           Temp_HIlim, Temp_LOlim
      use speedchem,       only: ns, nr, A0, b0, E0, Ainf, binf, Einf
      use SCthermodata,    only: interp4_coefs, use_table
      use troepar,         only: ntbFALL, itbFALL
      use universal_constants, only: uRcal, u2, u3, u6, one

      implicit none

      real (dp)       , dimension(6),       intent(in)  :: Ta
      real (dp)       , dimension(ntbFALL), intent(out) :: k0
      real (dp)       , dimension(nr),      intent(out) :: kinf
      integer,          optional,           intent(in)  :: iT
      real (dp)       , optional,           intent(in)  :: frac
      real (dp)                                         :: uRcT,coefs(5)
      integer :: i, ii, j
      logical, parameter :: parabolic = .true.
      real (dp)       , dimension(:,:), pointer    :: kpoint

      if ( (.not.use_table(Ta(1)))   .or. (.not.present(iT)) .or.
     &    (.not.present(frac))  )   then

          uRcT = uRcal * Ta(5)

          if (ntbFALL>0)
     &    k0 = A0(itbFALL) * exp(b0(itbFALL)* Ta(6) - E0(itbFALL)* uRcT)

!          kinf = retrieve_kinf(Ta)
          kinf    = Ainf   * exp(binf   * Ta(6) - Einf   * uRcT)


      else !tabulated data

          tabulation_order: if (parabolic) then

          coefs = interp4_coefs(frac)

!		  kinf = retrieve_kinf(Ta, iT, frac, coefs)

           if (ntbFALL>0)
     &             k0     = tab_k0  (:,iT-3)*coefs(1) +
     &                      tab_k0  (:,iT-2)*coefs(2) +
     &                      tab_k0  (:,iT-1)*coefs(3) +
     &                      tab_k0  (:,iT  )*coefs(4) +
     &                      tab_k0  (:,iT+1)*coefs(5)

                   kinf   = tab_kinf(:,iT-3)*coefs(1) +
     &                      tab_kinf(:,iT-2)*coefs(2) +
     &                      tab_kinf(:,iT-1)*coefs(3) +
     &                      tab_kinf(:,iT  )*coefs(4) +
     &                      tab_kinf(:,iT+1)*coefs(5)


          else

          if (ntbFALL>0)
     &    k0    = tab_k0  (:,iT-1)*(one-frac) + tab_k0  (:,iT)*frac
          kinf  = tab_kinf(:,iT-1)*(one-frac) + tab_kinf(:,iT)*frac

          endif tabulation_order

      endif

      end subroutine reaction_rates


      function retrieve_kinf(Ta, iT, frac, coefs) result(kinf)

      use chemistry_setup, only: accurate_scthermo,Temp_HIlim,Temp_LOlim
      use speedchem,           only: nr, Ainf, binf, Einf
      use scthermodata,        only: use_table
      use universal_constants, only: uRcal, one

      implicit none

      real (dp)       , dimension(6),  intent(in)  :: Ta

      integer,          optional,      intent(in)  :: iT
      real (dp)       , optional,      intent(in)  :: frac
      real (dp)       , optional,      intent(in)  :: coefs(5)

      real (dp)       , dimension(nr)              :: kinf
      real (dp)                                    :: uRcT
      integer                                      :: i, j
      logical,          parameter                  :: parabolic = .true.

      if (.not.allocated(store_kinf))allocate(store_kinf(nr))

      stored: if (abs(Ta(1) - store_kinf_T)/Ta(1) <= T_RTOL) then
         kinf = store_kinf
      else

         if ( (.not.use_table(Ta(1))) .or. (.not.present(iT)) .or.
     &       (.not.present(frac)) .or. (.not.present(coefs)) )   then
                uRcT = uRcal * Ta(5)
                kinf    = Ainf   * exp(binf   * Ta(6) - Einf   * uRcT)
            else !tabulated data

             tabulation_order: if (parabolic) then

                kinf   = tab_kinf(:,iT-3)*coefs(1) +
     &                   tab_kinf(:,iT-2)*coefs(2) +
     &                   tab_kinf(:,iT-1)*coefs(3) +
     &                   tab_kinf(:,iT  )*coefs(4) +
     &                   tab_kinf(:,iT+1)*coefs(5)

             else

                kinf   = tab_kinf(:,iT-1)*(one-frac)
     &                 + tab_kinf(:,iT  )*frac

             endif tabulation_order

         endif

         store_kinf   = kinf
         store_kinf_T = Ta(1)

      endif stored



      end function retrieve_kinf


c     *****************************************************************
c     **                                                             **
c     **         Compute Explicit backward reaction rates            **
c     **                                                             **
c     *****************************************************************
      subroutine explicit_rev_reaction_rates(Ta, kbXrev, iT, frac)

      use chemistry_setup, only: accurate_scthermo,
     &                           Temp_HIlim, Temp_LOlim
      use speedchem,       only: ns, nr, ARi, bRi, ERi
      use reacpar,         only: nXREV
      use SCthermodata,    only: interp4_coefs, use_table
      use troepar,         only: ntbFALL, itbFALL
      use universal_constants, only: uRcal, u2, u3, u6, one

      implicit none

      real (dp)       , dimension(6),     intent(in)  :: Ta
      real (dp)       , dimension(nXREV), intent(out) :: kbXrev
      integer,          optional,         intent(in)  :: iT
      real (dp)       , optional,         intent(in)  :: frac
      real (dp)                                       :: uRcT, coefs(5)
      integer                                         :: ii, j
      logical, parameter :: parabolic = .true.

      if ( (.not.use_table(Ta(1))) .or. (.not.present(iT)) .or.
     &    (.not.present(frac)) .or.
     &    (Ta(1)>Temp_HIlim .or. Ta(1)<Temp_LOlim) )   then

          uRcT = uRcal * Ta(5)
          kbXrev = ARi*exp(bRi*Ta(6)-ERi*uRcT)


      else !tabulated data

          tabulation_order: if (parabolic) then

            coefs = interp4_coefs(frac)

            kbXrev  = tab_Xkb(:,iT-3)*coefs(1) +
     &                tab_Xkb(:,iT-2)*coefs(2) +
     &                tab_Xkb(:,iT-1)*coefs(3) +
     &                tab_Xkb(:,iT  )*coefs(4) +
     &                tab_Xkb(:,iT+1)*coefs(5)

          else

            kbXrev  = tab_Xkb(:,iT-1)*(one-frac) + tab_Xkb(:,iT)*frac

          endif tabulation_order

      endif

      end subroutine explicit_rev_reaction_rates

c     *****************************************************************
c     **                                                             **
c     **   Compute Explicit backward reaction rates and derivatives  **
c     **                                                             **
c     *****************************************************************
      subroutine explicit_rev_reaction_rates_deriv
     &                    (Ta, kbXrev, dkbXrevdT, iT, frac)

      use chemistry_setup, only: accurate_scthermo,
     &                           Temp_HIlim, Temp_LOlim
      use speedchem,       only: ns, nr, ARi, bRi, ERi
      use reacpar,         only: nXREV
      use SCthermodata,    only: interp4_coefs, use_table
      use troepar,         only: ntbFALL, itbFALL
      use universal_constants, only: uRcal, u2, u3, u6, one

      implicit none

      real (dp)       , dimension(6),     intent(in)  :: Ta
      real (dp)       , dimension(nXREV), intent(out) :: kbXrev
      real (dp)       , dimension(nXREV), intent(out) :: dkbXrevdT
      integer,          optional,         intent(in)  :: iT
      real (dp)       , optional,         intent(in)  :: frac
      real (dp)                                       :: uRcT, coefs(5)
      integer                                         :: ii, j
      logical, parameter :: parabolic = .true.

      if ( (.not.use_table(Ta(1)))   .or. (.not.present(iT)) .or.
     &    (.not.present(frac))  )   then

          uRcT = uRcal * Ta(5)
          kbXrev    = ARi*exp(bRi*Ta(6)-ERi*uRcT)
          dkbXrevdT = kbXrev * Ta(5) * (bRi + Eri * uRcT)

      else !tabulated data

          tabulation_order: if (parabolic) then

            coefs = interp4_coefs(frac)

            kbXrev  = tab_Xkb(:,iT-3)*coefs(1) +
     &                tab_Xkb(:,iT-2)*coefs(2) +
     &                tab_Xkb(:,iT-1)*coefs(3) +
     &                tab_Xkb(:,iT+0)*coefs(4) +
     &                tab_Xkb(:,iT+1)*coefs(5)

            dkbXrevdT  = tab_Xdkbdt(:,iT-3)*coefs(1) +
     &                   tab_Xdkbdt(:,iT-2)*coefs(2) +
     &                   tab_Xdkbdt(:,iT-1)*coefs(3) +
     &                   tab_Xdkbdt(:,iT+0)*coefs(4) +
     &                   tab_Xdkbdt(:,iT+1)*coefs(5)

          else

            kbXrev  = tab_Xkb(:,iT-1)*(one-frac) + tab_Xkb(:,iT)*frac

            dkbXrevdT = tab_Xdkbdt(:,iT-1)*(one-frac)
     &                + tab_Xdkbdt(:,iT  )*frac

          endif tabulation_order

      endif

      end subroutine explicit_rev_reaction_rates_deriv



c     *****************************************************************
c     **   Compute ordinary forward reaction rates and derivatives   **
c     *****************************************************************
      subroutine reaction_rates_and_derivative
     &           (Ta, k0, dk0dT, kinf, dkinfdT, iT, frac)

      use chemistry_setup, only: accurate_scthermo,
     &                           Temp_HIlim, Temp_LOlim
      use speedchem,       only: ns,nr, A0, b0, E0, Ainf, binf, Einf
      use troepar,         only: itbFALL,ntbFALL
      use SCthermodata,    only: interp4_coefs, use_table
      use universal_constants, only: uRcal, u2, one

      implicit none

      real (dp)       , dimension(6),       intent(in)  :: Ta
      real (dp)       , dimension(nr),      intent(out) :: kinf,dkinfdt
      real (dp)       , dimension(ntbFALL), intent(out) :: k0,dk0dT
      integer,          optional,           intent(in)  :: iT
      real (dp)       , optional,           intent(in)  :: frac
      real (dp)                                         :: uRcT,coefs(5)
      integer                                           :: ii, j
      logical, parameter :: parabolic = .true.

      uRcT = uRcal * Ta(5)

      if ( (.not.use_table(Ta(1)))   .or. (.not.present(iT)) .or.
     &    (.not.present(frac)) .or.
     &    (Ta(1)>Temp_HIlim .or. Ta(1)<Temp_LOlim) )   then

          if (ntbFALL>0) then
            k0             = A0(itbFALL)
     &                     * exp(b0(itbFALL)* Ta(6) - E0(itbFALL)* uRcT)
            dk0dT          = k0
     &                     * Ta(5) * (b0(itbFALL) + E0(itbFALL) * uRcT)
          endif
          kinf     = Ainf   * exp(binf   * Ta(6) - Einf   * uRcT)

!          kinf = retrieve_kinf(Ta)

          dkinfdT  = kinf   * Ta(5) * (binf + Einf * uRcT)



      else !tabulated data

          tabulation_order: if (parabolic) then

             coefs = interp4_coefs(frac)

!             kinf = retrieve_kinf(Ta, iT, frac, coefs)

             if (ntbFALL>0) then

              k0          =       tab_k0(:,iT-3)*coefs(1) +
     &                            tab_k0(:,iT-2)*coefs(2) +
     &                            tab_k0(:,iT-1)*coefs(3) +
     &                            tab_k0(:,iT  )*coefs(4) +
     &                            tab_k0(:,iT+1)*coefs(5)

             dk0dT        =       tab_dk0dT(:,iT-3)*coefs(1) +
     &                            tab_dk0dT(:,iT-2)*coefs(2) +
     &                            tab_dk0dT(:,iT-1)*coefs(3) +
     &                            tab_dk0dT(:,iT  )*coefs(4) +
     &                            tab_dk0dT(:,iT+1)*coefs(5)
             endif

           kinf           =       tab_kinf(1:nr,iT-3)*coefs(1) +
     &                            tab_kinf(1:nr,iT-2)*coefs(2) +
     &                            tab_kinf(1:nr,iT-1)*coefs(3) +
     &                            tab_kinf(1:nr,iT  )*coefs(4) +
     &                            tab_kinf(1:nr,iT+1)*coefs(5)



           dkinfdT        =       tab_dkinfdT(1:nr,iT-3)*coefs(1) +
     &                            tab_dkinfdT(1:nr,iT-2)*coefs(2) +
     &                            tab_dkinfdT(1:nr,iT-1)*coefs(3) +
     &                            tab_dkinfdT(1:nr,iT  )*coefs(4) +
     &                            tab_dkinfdT(1:nr,iT+1)*coefs(5)

          else

      if (ntbFALL>0) then
       k0    = tab_k0   (:,iT-1)*(one-frac)+tab_k0   (:,iT)*frac
       dk0dT = tab_dk0dT(:,iT-1)*(one-frac)+tab_dk0dT(:,iT)*frac
      endif
       kinf    = tab_kinf   (:,iT-1)*(one-frac) +tab_kinf   (:,iT)*frac
       dkinfdT = tab_dkinfdT(:,iT-1)*(one-frac) +tab_dkinfdT(:,iT)*frac

          endif tabulation_order

      endif

      end subroutine reaction_rates_and_derivative

c     *****************************************************************
c     **   Compute ordinary forward reaction rates derivatives only  **
c     *****************************************************************
      subroutine reaction_rates_derivatives
     &           (Ta, dk0dT, dkinfdT, iT, frac)

      use chemistry_setup, only: accurate_scthermo,
     &                           Temp_HIlim, Temp_LOlim
      use speedchem,       only: ns,nr, A0, b0, E0, Ainf, binf, Einf
      use troepar,         only: itbFALL,ntbFALL
      use SCthermodata,    only: interp4_coefs, use_table
      use universal_constants, only: uRcal, u2, one

      implicit none

      real (dp)       , dimension(6),  intent(in)  :: Ta
      real (dp)       , dimension(nr), intent(out) :: dkinfdT
      real (dp)       , dimension(ntbFALL), intent(out) :: dk0dT
      integer,          optional,      intent(in)  :: iT
      real (dp)       , optional,      intent(in)  :: frac
      real (dp)                                    :: uRcT, coefs(5)
      integer :: ii, j
      logical, parameter :: parabolic = .true.

      uRcT = uRcal * Ta(5)

      if ( (.not.use_table(Ta(1))) .or. (.not.present(iT)) .or.
     &    (.not.present(frac)) .or.
     &    (Ta(1)>Temp_HIlim .or. Ta(1)<Temp_LOlim) )   then

          if (ntbFALL>0)
     &    dk0dT(itbFALL) = A0(itbFALL)
     &                   * exp(b0(itbFALL)* Ta(6) - E0(itbFALL)* uRcT)
     &                   * Ta(5) * (b0(itbFALL) + E0(itbFALL) * uRcT)

          dkinfdT  = Ainf   * exp(binf   * Ta(6) - Einf   * uRcT)
     &             * Ta(5)  * (binf + Einf * uRcT)

      else !tabulated data

          tabulation_order: if (parabolic) then

             coefs = interp4_coefs(frac)
             if (ntbFALL>0) then

           dk0dT          = tab_dk0dT(:,iT-3)*coefs(1) +
     &                            tab_dk0dT(:,iT-2)*coefs(2) +
     &                            tab_dk0dT(:,iT-1)*coefs(3) +
     &                            tab_dk0dT(:,iT  )*coefs(4) +
     &                            tab_dk0dT(:,iT+1)*coefs(5)
             endif

           dkinfdT        = tab_dkinfdT(:,iT-3)*coefs(1) +
     &                            tab_dkinfdT(:,iT-2)*coefs(2) +
     &                            tab_dkinfdT(:,iT-1)*coefs(3) +
     &                            tab_dkinfdT(:,iT  )*coefs(4) +
     &                            tab_dkinfdT(:,iT+1)*coefs(5)

          else

            if (ntbFALL>0)
     &      dk0dT   = tab_dk0dT(:,iT-1)*(one-frac)+tab_dk0dT(:,iT)*frac

            dkinfdT = tab_dkinfdT(:,iT-1)*(one-frac)
     &              + tab_dkinfdT(:,iT)*frac

          endif tabulation_order

      endif

      end subroutine reaction_rates_derivatives


c     *****************************************************************
c     **                                                             **
c     **                     SpeedCHEM FORTRAN                       **
c     **                                                             **
c     **       Compute species molar change rates [mol/cm3/s]        **
c     **                                                             **
c     **                                                             **
c     **   Author:      Federico Perini                              **
c     **   Last update: sunday, 03/12/2011                           **
c     **                                                             **
c     *****************************************************************
      function SCdwdt(Ta,Y,iT,frac,nonzeroes) result(dwdt)


      use chemistry_setup,     only: accurate_scthermo,
     &                               simplified_for_sparsity,
     &                               Temp_HIlim, Temp_LOlim,
     &                               save_thermal_parameters
      use speedchem,           only: ns, nr, iTB, nTB,
     &                               ARi, bri, ERi, reversibile,inotrev,
     &                               A0 , b0 , E0 , Ainf, binf, Einf,
     &                               third_body_beta, nTHREE, iTHREE,
     &                               LOW, HIGH, uMW, nnotrev
      use troepar,             only: itbALL, itbTROE, iTROE4, itbLIND,
     &                               ntbALL, ntbTROE, nTROE4, ntbLIND,
     &                               ntbSIMP, itbSIMP,iSIMPitbALL,
     &                               iLINDitbALL, iTROEitbALL, ntbFALL,
     &                               todotroe, aT2, uT1T2, T2T2, uT3T2,
     &                               zeroT2, troe_logfac, iFALLitbALL,
     &                               itbFALL, iTROEitbFALL, iLINDitbFALL
      use reacpar,             only: troereac, Lindreac, uequilC,
     &                               n_tb_beta,tb_beta_pack,
     &                               nEQREV,iEQREV,nXREV,iXREV,
     &                               is_beta_pack
      use SCmixturethermo,     only: SCP, SCrho
      use SCthermodata,        only: interp4_coefs, use_table
      use sparse_chemistry,    only: nudiffT_sparse, tb_beta_sp,
     &                               istoich_r_sp, istoich_p_sp,
     &                               inudiffT_sparse
      use universal_constants, only: uRcal, uR, u2, u3, u6, milli, one,
     &                               ten
      use sparse_definitions

      implicit none

      real (dp)       , dimension(6),      intent(in) :: Ta
      real (dp)       , dimension(ns),     intent(in) :: Y
      logical,          dimension(ns),optional, intent(in) :: nonzeroes
      real (dp)       , optional,          intent(in) :: frac
      integer,          optional,          intent(in) :: iT

      real (dp)       , dimension(ns)      :: YY, dwdt

      real (dp)       , dimension(nr)      :: q, qf, qb, kf, kb,
     &                                        uKc, kinf

      real (dp)       , dimension(ntbFALL) :: k0, Pr, k0M
      real (dp)       , dimension(ntbALL)  :: M


      real (dp)       , dimension(ns)      :: C, dCtot_dY, X


      real (dp)       , dimension(ntbTROE) :: Pr2,Fcent,troecor,Pr3,
     &                                        log10Fcent,
     &                                        log10Pr,fattore,
     &                                        ctroe,ntroe,logF,
     &                                        ten_pow_logF

      real (dp)       , dimension(nEQREV)  :: uKceq
      real (dp)       , dimension(nXREV)   :: kbXrev

      real (dp)                            :: uRcT, uRT, C_pow_nu, Ctot,
     &                                        uPr0, coefs(5),
     &                                        facf, facb

      real (dp)       , parameter          :: dtroe = 0.14e0_dp
      real (dp)       , parameter          :: zero  = 0.0_dp

      integer :: i, j, jf, jb, if, ib, ii, iii, ire, ivsp, ivsr, ifl

      logical :: use_tabulated_thermo

      logical, parameter :: parabolic = .true.


c     *****************************************************************

c     Species concentrations [mol/cm3]
      dCtot_dY = milli * SCrho * uMW
      C        = dCtot_dY * Y
      Ctot     = sum(C)

c     ** Compute forward reaction rates *******************************
      use_tabulated_thermo = use_table(Ta(1))

c     computing reaction rate constants [cm^3, mol, s, K]
      if (.not.use_tabulated_thermo) then
         call reaction_rates(Ta, k0, kinf)
         if (ntbTROE>0)log10Fcent = troe_logfac(Ta)
      else
         call reaction_rates(Ta, k0, kinf, iT, frac)
         if (ntbTROE>0)log10Fcent = troe_logfac(Ta,iT,frac)
      endif

c     ** Compute forward reaction rates *******************************
      kf = kinf

c     ** Compute backward reaction rates from either equilibrium, *****
c     ** specified Arrhenius parameters or irreversible behaviour *****

c     1) Non-reversible reactions
      if(nnotrev>0)kb(inotrev) = zero


c     2) Equilibrium-computed backward reaction rates
      if (nEQREV > 0) then
        if (.not.use_tabulated_thermo) then
           uKceq = uequilC(Ta)
        else
           uKceq = uequilC(Ta,iT,frac)
        endif
        kb(iEQREV) = kf(iEQREV) * uKceq
      endif

c     2) Explicit backward reaction rates
      if (nXREV > 0) then
        if (.not.use_tabulated_thermo) then
           kbXrev = ARi*exp(bRi*Ta(6)-ERi*uRcal * Ta(5))
        else
           call explicit_rev_reaction_rates(Ta, kbXrev, iT, frac)
        endif
        kb(iXREV) = kbXrev
      endif

c     Computing effective molecularity of the reactions
      if (ntbALL > 0) M = Ctot - tb_beta_sp * C

c     Computing reduced pressure values for falloff reactions
      reduced_pressures: if (ntbFALL > 0) then
         k0M = k0 * M(iFALLitbALL)
         Pr  = k0M/(kinf(itbFALL) + k0M)
      end if reduced_pressures


c     Updating forward reaction rate constant for Lindemann reactions
c      kf(itbLIND) = kf(itbLIND) * Pr(iLINDitbALL)
      if (ntbLIND > 0) then
        kf(itbLIND) = kf(itbLIND) * Pr(iLINDitbFALL)
        kb(itbLIND) = kb(itbLIND) * Pr(iLINDitbFALL)
      endif

c      lindemann: do j = 1, ntbLIND
c         ire = itbLIND(j)
c         ifl = iLINDitbFALL(j)
c         kf(ire) = kf(ire) * Pr(ifl)
c         kb(ire) = kb(ire) * Pr(ifl)
cc      end do lindemann
c      endif

c     ** TROE *********************************************************
c     Computing reaction rate constants according to Troe's formulation
      troefactors: if (ntbTROE > 0) then

         Pr2 = Pr(iTROEitbFALL)

c        ** Computing troe centering factor
         log10Pr    = log10(Pr2/(one-Pr2))

c        ** Troe model parameters
         ctroe = -0.40_dp - 0.67_dp * log10Fcent
         ntroe =  0.75_dp - 1.27_dp * log10Fcent

         fattore = ((log10Pr+ctroe)/(ntroe-dtroe*(log10Pr+ctroe)))
         logF = log10Fcent/(one + fattore*fattore)

         ten_pow_logF = ten ** logF
         troecor      = Pr2 * ten_pow_logF

         kf(itbTROE)  = kf(itbTROE) * troecor
         kb(itbTROE)  = kb(itbTROE) * troecor

      endif troefactors

c     ** Law of mass action productories
      call mass_action_productories(C,qf,qb)

      qf = qf * kf
      qb = qb * kb

c     ** Complete reaction progress variable formulation **************
      q = qf - qb

c     ** Apply molecularities of simple thirdbody reactions ***********
      if (ntbSIMP > 0) q(itbSIMP) = q(itbSIMP) * M(iSIMPitbALL)

c     ** Final, sparse multiplication: dwdt = matmul(nudiffT,q) *******
      dwdt = nudiffT_sparse * q

      end function SCdwdt

c     *****************************************************************
c     **                                                             **
c     ** Compute (or retrieve, if already available) reaction rates  **
c     **                                                             **
c     *****************************************************************
      function retrieve_dwdt(q,T) result(dwdt)
      use speedchem,        only: ns, nr
      use sparse_chemistry, only: nudiffT_sparse
      use sparse_definitions

      implicit none

      real (dp)       , dimension(nr), intent(in) :: q
      real (dp)       ,                intent(in) :: T
      real (dp)       , dimension(ns)             :: dwdt

      integer :: j
      if (.not.allocated(store_dwdt))allocate(store_dwdt(ns))

      stored: if (abs(T - store_dwdt_T)/T <= T_RTOL) then

         dwdt = store_dwdt

      else
c        ** sparse multiplication: dwdt = matmul(nudiffT,q) *******
         dwdt         = nudiffT_sparse * q
         store_dwdt_T = T
         store_dwdt   = dwdt
      endif stored


      end function retrieve_dwdt


c     *****************************************************************
c     **                                                             **
c     ** Compute (or retrieve, if already available) mass action     **
c     ** productories for the reactions:                             **
c     **   qf(k) = product(C(i)**stoich_r(k,i)), i=1,...,ns          **
c     **   qb(k) = product(C(i)**stoich_p(k,i)), i=1,...,ns          **
c     **   for all k = 1, ..., nr                                    **
c     **                                                             **
c     *****************************************************************

      subroutine mass_action_productories(C,qf,qb)
      use speedchem,      only: ns, nr
      use chemistry_setup, only: save_thermal_parameters
      use sparse_chemistry, only: istoich_r_eff_sp,istoich_p_eff_sp,
     &                            row_power_to_sparse_product
      use universal_constants, only: micro, one
      implicit none

      real (dp)       , dimension(ns), intent(in)  :: C
      real (dp)       , dimension(nr), intent(out) :: qf, qb


      logical          :: out_of_tolerance
      real (dp)        :: C_pow_nu
      integer          :: i, ii, iii, j, ivsr, ivsp


      if (save_thermal_parameters) then

c        ** Check for allocation of the storage arrays ****************
         if (.not.allocated(store_C ))allocate(store_C (ns))
         if (.not.allocated(store_qf))allocate(store_qf(nr))
         if (.not.allocated(store_qb))allocate(store_qb(nr))

         out_of_tolerance = .false.

         check_tolerance: do i = 1, ns
           if (abs(1.e0_dp - C(i)*store_C(i)) > micro) then
              out_of_tolerance = .true.
              exit check_tolerance
           endif
         end do check_tolerance

      else

         out_of_tolerance = .true.

      endif

      evaluation_needed: if (out_of_tolerance) then

c        ** Forward progress variable formulation *********************
         qf = 1.e0_dp
         forward_progress_variable: do j=1,size(i2r)
           iii     = i2D1r(j)
           ii      = ijr(j)
           ivsr    = iv_stoich_r(j)

           C_pow_nu = C(ii)
           do i = 1, ivsr-1
              C_pow_nu = C_pow_nu * C(ii)
           end do
           qf(iii) = qf(iii) * C_pow_nu

          end do forward_progress_variable
c          qf = row_power_to_sparse_product(C,istoich_r_eff_sp)

c         ** Backward progress variable formulation *******************

!          qb = row_power_to_sparse_product(C,istoich_p_eff_sp)
          qb = 1.e0_dp
          if (allocated(i2p)) then
             backward_progress_variable: do j=1,size(i2p)
               iii     = i2D1p(j)
               ii      = ijp(j)
               ivsp    = iv_stoich_p(j)
               C_pow_nu = C(ii)

               do i = 1, ivsp-1
                  C_pow_nu = C_pow_nu * C(ii)
               end do
               qb(iii) = qb(iii) * C_pow_nu

            end do backward_progress_variable
         endif

c        Store values for retrieval
         if (save_thermal_parameters) then
            store_C  = one/C
            store_qf = qf
            store_qb = qb
         endif

      else

         qf = store_qf
         qb = store_qb

      endif evaluation_needed

      end subroutine mass_action_productories


      subroutine tabulate_kinetics
c     *****************************************************************
c     **                                                             **
c     **                     SpeedCHEM FORTRAN                       **
c     **                                                             **
c     **     Tabulate functions for computing reaction rates         **
c     **                                                             **
c     **                                                             **
c     **   Author:      Federico Perini                              **
c     **   Last update: sunday, 30/10/2011                           **
c     **                                                             **
c     *****************************************************************

      use speedchem, only: nr,A0,b0,E0,Ainf,binf,Einf,ARi,bRi,ERi
      use reacpar,   only: nXREV
      use chemistry_setup,only: Temp_table_accuracy, Temp_LOlim,
     &                     Temp_HIlim, tab_nsteps
      use troepar,   only: ntbFALL, itbFALL, ltbFALL
      use universal_constants, only: Rcal, uRcal, zero, one
      implicit none

      integer :: i, nA0
      real (dp)        :: T, uT
      real (dp)       , dimension(nr) :: k0, kinf
      real (dp)       , dimension(nXREV) :: kb

c     String formats
      character(len=*), parameter ::
     &      fmt_end = "(' kinetics parameters tabulated. ')"

c     Compute number of subdivision in temperature and Pressure ratio
c     ranges
      tab_nsteps = int( ( Temp_HIlim - Temp_LOlim )
     &                    / Temp_table_accuracy)

c     ** Allocate data
      nA0 = count(A0 /= zero)

      if (ntbFALL>0) then
         allocate(   tab_k0   (ntbFALL,-1:tab_nsteps),
     &               tab_dk0dT(ntbFALL,-1:tab_nsteps)     )

         tab_k0    = zero
         tab_dk0dT = zero

c         iA0 = pack([(i,i=1,nr)],mask = A0/=0.e0_dp)

      endif

      if (nXREV>0) then

         allocate(   tab_Xkb   (nXREV  ,-1:tab_nsteps) ,
     &               tab_Xdkbdt(nXREV  ,-1:tab_nsteps) )
         tab_Xkb    = zero
         tab_Xdkbdt = zero

      endif

      allocate(tab_kinf   (nr,-1:tab_nsteps),
     &         tab_dkinfdt(nr,-1:tab_nsteps))
      tab_kinf    = zero
      tab_dkinfdT = zero

      do i=1,tab_nsteps + 2

         T  = Temp_LOlim + real(i-2, dp) * Temp_table_accuracy
         uT = one/T

c     ** Computing forward reaction rate constants ********************
         k0   = A0*T**b0*exp(-E0*uRcal*uT)
         kinf = Ainf*T**binf*exp(-Einf*uRcal*uT)

         if (ntbFALL>0) then
           tab_k0   (:,i-2) = pack(k0,ltbFALL)
           tab_dk0dt(:,i-2) = pack(k0*uT*(b0 + E0*uRcal*uT),ltbFALL)
         endif

         tab_kinf   (:,i-2) = kinf
         tab_dkinfdt(:,i-2) = kinf*uT*(binf + Einf*uRcal*uT)

c     ** Computing backward reaction rate constants *******************
         if (nXREV>0) then
          kb                = ARi * T**bRi*exp(-ERi * uRcal * uT)
          tab_Xkb   (:,i-2) = kb
          tab_Xdkbdt(:,i-2) = kb * uT * (bRi + ERi * uRcal * uT)
         endif



      end do

      write(*,fmt_end)

      end subroutine tabulate_kinetics


c        **************************************************************
c          Permutates all the arrays containing species index
c          into a new species ordering previously defined
c        **************************************************************
         subroutine permutate_kinetics(new_order)
         implicit none

         integer, dimension(:), intent(in)   :: new_order
         integer                             :: i, j, tmp
         integer, dimension(size(new_order)) :: inverse_index
         integer, dimension(size(ijr))       :: sorted_indexr
         integer, dimension(size(ijp))       :: sorted_indexp

         inverse_index(new_order) = [(j,j=1,size(new_order))]

         ijr         = inverse_index(ijr)
         iv_stoich_r = inverse_index(iv_stoich_r)

c        Reorder rate of progress variable arrays for more efficient
c        evaluation of the productories q = prod(C)^stoich_coefs
c        i2D1r, ijr, iv_stoich_r
         sorted_indexr = [(j,j=1,size(ijr))]
         do i = 1, size(ijr)
           do j = i+1, size(ijr)
              swap_ij: if (ijr(i) > ijr(j)) then
                 tmp = sorted_indexr(j)
                 sorted_indexr(j) = sorted_indexr(i)
                 sorted_indexr(i) = tmp
              endif swap_ij
           end do
         end do

         ijr   = ijr  (sorted_indexr)
         i2D1r = i2D1r(sorted_indexr)
         iv_stoich_r = iv_stoich_r(sorted_indexr)

         if (allocated(ijp)) then
            ijp         = inverse_index(ijp)
            iv_stoich_p = inverse_index(iv_stoich_p)

c           Reorder rate of progress variable arrays for more efficient
c           evaluation of the productories q = prod(C)^stoich_coefs
c           i2D1p, ijp, iv_stoich_p
            sorted_indexp = [(j,j=1,size(ijp))]
            do i = 1, size(ijp)
              do j = i+1, size(ijp)
                 swap_ij2: if (ijp(i) > ijp(j)) then
                    tmp = sorted_indexp(j)
                    sorted_indexp(j) = sorted_indexp(i)
                    sorted_indexp(i) = tmp
                 endif swap_ij2
              end do
            end do

            ijp   = ijp  (sorted_indexp)
            i2D1p = i2D1p(sorted_indexp)
            iv_stoich_p = iv_stoich_p(sorted_indexp)



         endif




         end subroutine permutate_kinetics






      end module kinetics_mod

!     *****************************************************************
!     **                                                             **
!     **   ODE solver parameters and arrays                          **
!     **                                                             **
!     **   Author:      Federico Perini                              **
!     **   Last update: monday, 23/04/2012                           **
!     **                                                             **
!     *****************************************************************
      module ode_solver

      use working_precision, only: dp
      use sparse_definitions
      use dvode_f90_m, only: vode_opts

      implicit none
      public

!       Array length parameters
        integer :: lrw, liw
!OMP    THREADPRIVATE(lrw,liw)

!       Integration method parameters and option switches
        integer :: method, itol, iopt, itask, istate, ijac, imas, iout,
     &             ifcn, idfx, mumas, mlmas, mujac, mljac, maxk
!$OMP   THREADPRIVATE(method,itol,iopt,itask,istate,ijac,imas,iout,ifcn,
!$OMP&  idfx,mumas,mlmas,mujac,mljac,maxk)

!       Runtime parameters
        integer :: maxnsteps
        real (dp)                                   :: hs
        real (dp)                                   :: rtol
        real (dp)       , dimension(:), allocatable :: atol
!$OMP   THREADPRIVATE(maxnsteps, hs, rtol, atol)


!       Fortran 90 implementation of VODE related data
        type(vode_opts)                 :: vf90_opts
        integer,          dimension(31) :: istats
        double precision, dimension(22) :: rstats
        integer                         :: nIA, nJA, nPD, nVF90JAC
!$OMP   THREADPRIVATE(vf90_opts,istats,rstats,nIA,nJA,nPD,nVF90JAC)

!       Working arrays
        real (dp)       , dimension(:), allocatable :: rwork
        integer,          dimension(:), allocatable :: iwork
!$OMP   THREADPRIVATE(rwork,iwork)

!       Working arrays for RADAU5 sparse matrices storage
        integer,          dimension(:), allocatable :: iper1, iiper1
        integer,          dimension(:), allocatable :: iper2, iiper2
        integer,          dimension(:), allocatable :: iLU1,  iLU2
        real (dp)       , dimension(:), allocatable :: rLU1,  rLU2
        integer                                     :: liLU1, lrLU1,
     &                                                 liLU2, lrLU2,
     &                                                 strg1, strg2


        type(sparse)         :: JACT_VF90
        type(sparse_ordered) :: R5_sys1, R5_sys2, DASPK_sys

!$OMP   THREADPRIVATE(iper1,iiper1,iper2,iiper2,iLU1,iLU2,rLU1,rLU2,
!$OMP&                liLU1,lrLU1,liLU2,lrLU2,strg1,strg2)

!       Option array for DASPK ode solver
        integer,          dimension(20)             :: daspkinfo
!$OMP   THREADPRIVATE(daspkinfo)

!       --- Integration monitoring arrays ---
!       Calls to the ODE and to the jacobian routine
        integer(8) :: ncJAC   = 0
        integer(8) :: ncCONV  = 0
!       Number of integration steps
        integer(8) :: nsteps  = 0
!       Number of LU decompositions
        integer(8) :: nLUdec  = 0
!       Number of Newton iterations
        integer(8) :: nNewton = 0

!       OpenMP-related stuff
        integer :: maxnthreads, actnthreads


        contains

!          ************************************************************
!          ** Parallel initialisation                                **
!          ************************************************************

           subroutine set_parallel(nt)
!           use omp_lib
           implicit none

           integer, intent(in) :: nt


!          Get maximum number of threads
           maxnthreads = 1!omp_get_num_procs()

!          Set used number of threads
           if (nt>0) then
              actnthreads = min(nt, maxnthreads)
           else
              actnthreads = 1
           endif

!           call omp_set_num_threads(actnthreads)

           end subroutine set_parallel


!          ************************************************************
!          ** Dummy, empty subroutine                                **
!          ************************************************************
           subroutine dummy
           implicit none
           end subroutine dummy

      end module ode_solver