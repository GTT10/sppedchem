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

!     ******************************************************************
!     **                                                              **
!     **                    Sparse matrix library                     **
!     **                                                              **
!     **   A Fortran 2008 library for operations with sparse matrices **
!     **                                                              **
!     **   Author:      (C) Federico Perini                           **
!     **                                                              **
!     **   Date created: sunday, 29/07/2012                           **
!     **   Last update : sunday, 29/07/2012                           **
!     **                                                              **
!     **                                                              **
!     ******************************************************************

      module sparse_definitions

!     ** Working precision has to be set in the working_precision ******
!     ** module, by setting dp for real numbers                       **
      use working_precision

      implicit none
      private

!     ** Procedures that can be called from outside the module *********
      public :: allocate, deallocate, assignment(=), operator(+), &
                operator(-), operator(*), operator(.backslash.),  &
                operator(==), allocated

!     ** Public, callable subroutines and functions
      public :: sparse_compress, optimal_ordering, sparseLU

      ! Print matrix details to screen or to output file
      public :: matrix_details


!     ** SPARSE MATRIX TYPES *******************************************
!     Define sparse matrix formats: double precision, integer, logical
!     NB: Universal format: CSR and Yale formats together

!     ******************************************************************
!     ** DOUBLE PRECISION sparse matrix format                        **
!     ******************************************************************
      type, public :: sparse
           integer   :: n  ! Number of sparse elements
           integer   :: nr ! Number of dense rows
           integer   :: nc ! Number of dense columns
           real (dp), dimension(:), allocatable :: A
           integer,          dimension(:), allocatable :: IA
           integer,          dimension(:), allocatable :: JA
           real (dp), dimension(:), pointer     :: val
           integer,          dimension(:), pointer     :: ptrB
           integer,          dimension(:), pointer     :: ptrE
           integer,          dimension(:), pointer     :: col
      end type sparse

!     ******************************************************************
!     ** INTEGER          sparse matrix format                        **
!     ******************************************************************
      type, public :: sparseint
           integer   :: n
           integer   :: nr ! Number of dense rows
           integer   :: nc ! Number of dense columns
           integer, dimension(:), allocatable :: A
           integer, dimension(:), allocatable :: IA
           integer, dimension(:), allocatable :: JA
           integer, dimension(:), pointer     :: val
           integer, dimension(:), pointer     :: ptrB
           integer, dimension(:), pointer     :: ptrE
           integer, dimension(:), pointer     :: col
      end type sparseint

!     ******************************************************************
!     ** LOGICAL          sparse matrix format                        **
!     ******************************************************************
      type, public :: sparselog
           integer   :: n
           integer   :: nr ! Number of dense rows
           integer   :: nc ! Number of dense columns
           logical, dimension(:), allocatable :: A
           integer, dimension(:), allocatable :: IA
           integer, dimension(:), allocatable :: JA
           logical, dimension(:), pointer     :: val
           integer, dimension(:), pointer     :: ptrB
           integer, dimension(:), pointer     :: ptrE
           integer, dimension(:), pointer     :: col
      end type sparselog


!     ******************************************************************
!     ** DOUBLE PRECISION square sparse matrix format WITH REORDERING **
!     ** (for the solution of sparse linear systems)                  **
!     **                                                              **
!     ** More than the traditional arrays, this format contains:      **
!     ** - permutations = the permutation indices for rows/columns    **
!     **                  (to reduce fillin during LU factorization)  **
!     ** - inverse_permutations = the inverse indices                 **
!     **        permutations(i) = j <=> inverse_permutations(j) = i   **
!     ** - real_space   = double precision space for LU factoriz.     **
!     **                  using Yale sparse matrix routines           **
!     ** - int_space    = integer space for LU factorization          **
!     **                                                              **
!     ** NB: type(sparse_ordered) matrices cannot be allocated by     **
!     **     zeroes, they have to be passed by equivalence or         **
!     **     assignment (=) from a type(sparse) matrix.               **
!     ** NB: real_space and int_space contain information including   **
!     **     the positions from symbolic factorization                **
!     ******************************************************************

!     NB if compiler supports Fortran 0X, the inherited type formulation
!        is better because all the subroutine involving type(sparse)
!        will also work with type(sparse_ordered)
!      type, public, extends(sparse) :: sparse_ordered
      type, public :: sparse_ordered

           integer   :: n  ! Number of sparse elements
           integer   :: nr ! Number of dense rows
           integer   :: nc ! Number of dense columns
           real (dp)       , dimension(:), allocatable :: A
           integer,          dimension(:), allocatable :: IA
           integer,          dimension(:), allocatable :: JA
           real (dp)       , dimension(:), pointer     :: val
           integer,          dimension(:), pointer     :: ptrB
           integer,          dimension(:), pointer     :: ptrE
           integer,          dimension(:), pointer     :: col

           logical                             :: is_ordered! = .false.
           integer, dimension(:), allocatable  :: perm
           integer, dimension(:), allocatable  :: inv_perm

           logical :: symbolically_factorized 
           logical :: numerically_factorized  
           integer                                     :: lR, lI
           integer,          dimension(:), allocatable :: int_space
           real (dp)       , dimension(:), allocatable :: real_space

      end type sparse_ordered

!     ******************************************************************
!     ** OVERLOADED OPERATORS                                         **
!     ******************************************************************

      ! Initialise an empty matrix by dimensions (nrows, ncols, nelems)
      interface allocate
         module procedure sparse_allocate
         module procedure sparse_ordered_allocate
         module procedure sparseint_allocate
      end interface allocate

      ! Initialise an empty matrix or deallocate a matrix
      ! [usage] call deallocate(matrix)
      interface deallocate
         module procedure sparse_init
         module procedure sparseint_init
         module procedure sparselog_init
         module procedure sparse_ordered_init
      end interface deallocate

      ! Assignment operator
      interface assignment (=)

         ! Assign dense format matrices to sparse format
         module procedure double_to_sparse
         module procedure double_to_sparseint
         module procedure sparse_to_double
         module procedure sparse_ordered_to_double

         ! Equivalence between two sparse matrix kinds
         module procedure sparse_equivalence
         module procedure sparseint_equivalence
         module procedure sparse_to_sparseint_equivalence
         module procedure sparseint_to_sparse_equivalence
         module procedure sparse_to_sparse_ordered_equivalence
         module procedure sparse_ordered_to_sparse_equivalence

      end interface !assignment (=)

      ! Sum and difference operators
      interface operator (+)
         module procedure sparse_plus_sparse
      end interface

      interface operator (-)
         module procedure sparse_minus_sparse
         module procedure sparse_negative
      end interface

      ! Product
      interface operator (*)
         module procedure sparse_dot_double
         module procedure double_dot_sparse
         module procedure sparse_matmul
         module procedure sparseint_matmul
      end interface

      ! Linear system operator!
      interface operator (.backslash.)
         module procedure sparse_ordered_linear_system
!         module procedure sparse_ordered_linear_system_double
         module procedure sparse_linear_system
!         module procedure sparse_linear_system_double
      end interface

      ! Comparison between sparse matrices
      interface operator (==)
         module procedure sparse_comparison_sparse
         module procedure sparse_ordered_comparison_sparse
         module procedure sparse_comparison_sparse_ordered
      end interface

      interface same_structure
         module procedure sparse_structurecomp_sparse
         module procedure sparse_ordered_structurecomp_sparse
         module procedure sparse_structurecomp_sparse_ordered
      end interface same_structure

      interface allocated
         module procedure sparse_ordered_allocated
         module procedure sparse_allocated
         module procedure sparseint_allocated
      end interface allocated

      interface sparseLU
         module procedure sparse_ordered_numerical_factorization
      end interface sparseLU

      interface matrix_details
         module procedure print_sparse_details
         module procedure print_sparseint_details
         module procedure print_sparse_ordered_details
      end interface matrix_details

      interface zeros
         module procedure zeros_sparse
         module procedure zeros_sparseint
         module procedure zeros_sparse_ordered
      end interface zeros

!     ******************************************************************
      contains

!     ******************************************************************
!     ** DEALLOCATION and INITIALISATION -related routines            **
!     ******************************************************************

!     ** sparse initialisation *****************************************
      subroutine sparse_init(spmatrix)
      implicit none

      type(sparse), intent(inout) :: spmatrix

      nullify(spmatrix%val )
      nullify(spmatrix%col )
      nullify(spmatrix%ptrB)
      nullify(spmatrix%ptrE)

      if (allocated(spmatrix%JA)  ) deallocate(spmatrix%JA)
      if (allocated(spmatrix%IA)  ) deallocate(spmatrix%IA)
      if (allocated(spmatrix%A )  ) deallocate(spmatrix%A)

      spmatrix%n  = -1
      spmatrix%nr = -1
      spmatrix%nc = -1

      end subroutine sparse_init

!     ** sparseint initialisation **************************************
      subroutine sparseint_init(spmatrix)
      implicit none

      type(sparseint), intent(inout) :: spmatrix

      nullify(spmatrix%val )
      nullify(spmatrix%col )
      nullify(spmatrix%ptrB)
      nullify(spmatrix%ptrE)

      if (allocated(spmatrix%JA)  ) deallocate(spmatrix%JA)
      if (allocated(spmatrix%IA)  ) deallocate(spmatrix%IA)
      if (allocated(spmatrix%A )  ) deallocate(spmatrix%A)

      spmatrix%n  = -1
      spmatrix%nr = -1
      spmatrix%nc = -1

      end subroutine sparseint_init

!     ** sparseint initialisation **************************************
      subroutine sparselog_init(spmatrix)
      implicit none

      type(sparselog), intent(inout) :: spmatrix

      nullify(spmatrix%val )
      nullify(spmatrix%col )
      nullify(spmatrix%ptrB)
      nullify(spmatrix%ptrE)

      if (allocated(spmatrix%JA)  ) deallocate(spmatrix%JA)
      if (allocated(spmatrix%IA)  ) deallocate(spmatrix%IA)
      if (allocated(spmatrix%A )  ) deallocate(spmatrix%A)

      spmatrix%n  = -1
      spmatrix%nr = -1
      spmatrix%nc = -1

      end subroutine sparselog_init

!     ** sparse initialisation *****************************************
      subroutine sparse_ordered_init(spmatrix)
      implicit none

      type(sparse_ordered), intent(inout) :: spmatrix

      nullify(spmatrix%val )
      nullify(spmatrix%col )
      nullify(spmatrix%ptrB)
      nullify(spmatrix%ptrE)

      if (allocated(spmatrix%JA)  ) deallocate(spmatrix%JA)
      if (allocated(spmatrix%IA)  ) deallocate(spmatrix%IA)
      if (allocated(spmatrix%A )  ) deallocate(spmatrix%A)
      if (allocated(spmatrix%perm )  ) deallocate(spmatrix%perm)
      if (allocated(spmatrix%inv_perm)  ) deallocate(spmatrix%inv_perm)
      if (allocated(spmatrix%real_space)) deallocate(spmatrix%real_space)
      if (allocated(spmatrix%int_space )) deallocate(spmatrix%int_space )

      spmatrix%n  = -1
      spmatrix%nr = -1
      spmatrix%nc = -1
      spmatrix%lR = -1
      spmatrix%lI = -1

      end subroutine sparse_ordered_init


!     ******************************************************************
!     ** ALLOCATION and DEALLOCATION -related routines                **
!     ******************************************************************

      ! Allocate sparse matrix with correct dimensions - but leave
      ! content blank
      subroutine sparse_allocate(nrows, ncols, nels, spmatrix)
      use utilities, only: extend
      implicit none

      integer,      intent(in)                :: nrows, ncols, nels
      type(sparse), intent(inout), target     :: spmatrix

      character(len=*), parameter ::                                   &
         fmt_ersize = "(' Error: non-2D matrix made sparse ',          &
                        ' in sparse_allocate_fun:' )"

!     A zero-row matrix is a valid representation for an optional
!     chemistry feature that is absent from a mechanism.  It still has
!     a well-defined column count and one CSR row-pointer sentinel.
      if (nrows<0.or.ncols<1.or.nels<0) then
         write(*,fmt_ersize)
         write(*,*)'nrow=',nrows,'nels=',nels,'ncols=',ncols
         stop
      endif

!     Initialize sparse matrix array
!      call deallocate(spmatrix)

!     Assign data
      spmatrix%n  = nels
      spmatrix%nr = nrows
      spmatrix%nc = ncols

!     Allocate A, IA, JA arrays
      call extend(spmatrix%A , nels)
      call extend(spmatrix%JA, nels)
      call extend(spmatrix%IA, nrows+1)

!     Associate pointers for CSR format
      spmatrix%val  => spmatrix%A  (1:nels)
      spmatrix%col  => spmatrix%JA (1:nels)
      spmatrix%ptrB => spmatrix%IA (1:nrows  )
      spmatrix%ptrE => spmatrix%IA (2:nrows+1)

      end subroutine sparse_allocate

      ! Allocate sparse matrix with correct dimensions - but leave
      ! content blank
      subroutine sparse_ordered_allocate(nrows, ncols, nels, spmatrix)
      use utilities, only: extend
      implicit none

      integer,              intent(in)                :: nrows, ncols, &
                                                         nels
      type(sparse_ordered), intent(inout), target     :: spmatrix

      character(len=*), parameter ::                                   &
         fmt_ersize = "(' Error: non-2D matrix made sparse ',          &
                        ' in sparse_allocate_fun:' )"

!     Halt on non-2D matrix
      if (nrows<1.or.ncols<1) then
         write(*,fmt_ersize)
         write(*,*)'nrow=',nrows,'nels=',nels,'ncols=',ncols
         stop
      endif

!     Initialize sparse matrix array
!      call deallocate(spmatrix)

!     Assign data
      spmatrix%n  = nels
      spmatrix%nr = nrows
      spmatrix%nc = ncols

!     Allocate A, IA, JA arrays
      call extend(spmatrix%A , nels)
      call extend(spmatrix%JA, nels)
      call extend(spmatrix%IA, nrows+1)

!     Associate pointers for CSR format
      spmatrix%val  => spmatrix%A  (1:nels)
      spmatrix%col  => spmatrix%JA (1:nels)
      spmatrix%ptrB => spmatrix%IA (1:nrows  )
      spmatrix%ptrE => spmatrix%IA (2:nrows+1)

!     Prepare room for ordering and factorization
      spmatrix%is_ordered = .false.
      spmatrix%numerically_factorized = .false.
      spmatrix%symbolically_factorized = .false.

      if (allocated(spmatrix%perm    ))deallocate(spmatrix%perm)
      if (allocated(spmatrix%inv_perm))deallocate(spmatrix%inv_perm)

      spmatrix %lR = 0
      spmatrix %lI = 0
      if (allocated(spmatrix%int_space ))deallocate(spmatrix%int_space)
      if (allocated(spmatrix%real_space))deallocate(spmatrix%real_space)

      end subroutine sparse_ordered_allocate


      ! Allocate sparse matrix with correct dimensions - but leave
      ! content blank
      subroutine sparseint_allocate(nrows, ncols, nels, spmatrix)
      use utilities, only: extend
      implicit none

      integer,         intent(in)                :: nrows, ncols, nels
      type(sparseint), intent(inout), target     :: spmatrix

      character(len=*), parameter ::                                   &
         fmt_ersize = "(' Error: non-2D matrix made sparse ',          &
                        ' in sparseint_allocate_fun:' )"

!     Halt on non-2D matrix
      if (nrows<1.or.ncols<1) then
         write(*,fmt_ersize)
         write(*,*)'nrow=',nrows,'nels=',nels,'ncols=',ncols
         stop
      endif

!     Assign data
      spmatrix%n  = nels
      spmatrix%nr = nrows
      spmatrix%nc = ncols

!     Allocate A, IA, JA arrays
      call extend(spmatrix%A , nels)
      call extend(spmatrix%JA, nels)
      call extend(spmatrix%IA, nrows+1)

!     Associate pointers for CSR format
      spmatrix%val  => spmatrix%A  (1:nels)
      spmatrix%col  => spmatrix%JA (1:nels)
      spmatrix%ptrB => spmatrix%IA (1:nrows  )
      spmatrix%ptrE => spmatrix%IA (2:nrows+1)

      end subroutine sparseint_allocate

!     ******************************************************************
!     ** CONVERSION of double precision matrix into SPARSE DOUBLE     **
!     ******************************************************************
      subroutine double_to_sparse(spmatrix,dense)
      implicit none

      real (dp)       , dimension(:,:), intent(in)    :: dense
      type(sparse)    , target,         intent(inout) :: spmatrix
      integer         , dimension(size(dense,1))      :: nonzero_row
      real (dp)       , parameter                     :: zero = 0.0_dp
      integer :: m,n,i,j,iia,nA, nIA, ntot,c, nr, nc, nemptyprev
      logical :: isfirst, updateprevious

      character(len=*), parameter ::                                   &
         fmt_ersize = "(' Error: non-2D matrix made sparse ',          &
                        ' in double_to_sparse:')"

      ! Gather matrix dimensions
      m = size(dense,1)
      n = size(dense,2)

!     Halt on non-2D matrix
      if (m*n/=size(dense)) then
         write(*,fmt_ersize)
         write(*,*)'nrow = ',m,'ncol = ',n
         stop
      endif

!     Store number of nonzero elements in each row
      nonzero_row = count(dense /= zero, dim=2)

!     Initialize sparse matrix array
      nA = sum(nonzero_row)

      nIA = size(dense,1)
      nr  = size(dense,1)
      nc  = size(dense,2)

!     Initialise sparse matrix
      call sparse_allocate(nr, nc, nA, spmatrix)

!     Store data into sparse format
      ntot = 0
      iia  = 0
      nemptyprev = 0

      row: do i = 1, spmatrix%nr
           isfirst = .true.
           c = nonzero_row(i)

           emptyrow: if (c==0) then ! Special treatment for empty row!

              updateprevious = .true.
              nemptyprev = nemptyprev + 1

           else

              column: do j = 1, spmatrix%nc

                 if (dense(i,j) /= zero) then
                    ntot = ntot + 1
                    spmatrix%A (ntot) = dense(i,j)
                    spmatrix%JA(ntot) = j

                    if (isfirst) then
                       iia     = i
                       isfirst = .false.
                       spmatrix%IA(iia:iia+1) = ntot
                       if (updateprevious) then
                           updateprevious = .false.
                           spmatrix%IA(iia-nemptyprev:iia-1) = ntot
                           nemptyprev         = 0
                       endif
                    endif
                 endif
              end do column
           endif emptyrow
      end do row

!     Storing last column length
      if (.not.updateprevious) then
         spmatrix%IA(m+1) = spmatrix%IA(m) + nonzero_row(m)
      else
         spmatrix%IA(m+1-nemptyprev:m+1) = ntot+1
      endif

!     Correct storage check
      if (nA/=ntot.or.ntot/=spmatrix%IA(m+1)-1) then
         write(*,*)'sparse storage error: ntot     = ',ntot
         write(*,*)'                      nA       = ',nA
         write(*,*)'                      IA(nr+1) = ',spmatrix%IA(m+1)
         write(*,*)'matrix rows:   nr = ',spmatrix%nr
         write(*,*)'matrix cols:   nc = ',spmatrix%nc
         stop
      endif

      end subroutine double_to_sparse


      subroutine sparse_to_double(dense, matrix)
      implicit none

      real (dp)       , dimension(:,:), allocatable, intent(inout) :: dense
      type (sparse)   , target,                      intent(in)    :: matrix
      real (dp)       , parameter :: zero = 0.0_dp
      integer                     :: i, j, jj


      if (allocated(dense)) deallocate(dense)

      if (matrix%nr == 0 .or. matrix%nc == 0 .or. matrix%n == 0) return

      ! Allocate matrix
      allocate(dense(matrix%nr, matrix%nc))

      ! Initialise matrix
      dense = zero

      do i = 1, matrix%nr
         do jj = matrix%IA(i), matrix%IA(i+1) - 1
            j  = matrix%JA(jj)
            dense(i,j) = matrix%A(jj)
         end do
      end do


      end subroutine sparse_to_double

      subroutine sparse_ordered_to_double(dense, matrix)
      implicit none

      real (dp)       , dimension(:,:), allocatable, intent(inout) :: dense
      type(sparse_ordered), target,                  intent(in)    :: matrix
      integer :: i, j, jj


      if (allocated(dense)) deallocate(dense)

      if (matrix%nr == 0 .or. matrix%nc == 0 .or. matrix%n == 0) return

      ! Allocate matrix
      allocate(dense(matrix%nr, matrix%nc))

      ! Initialise matrix
      dense = 0.0_dp

      do i = 1, matrix%nr
         do jj = matrix%IA(i), matrix%IA(i+1) - 1
            j  = matrix%JA(jj)
            dense(i,j) = matrix%A(jj)
         end do
      end do

      end subroutine sparse_ordered_to_double

!     ******************************************************************
!     ** CONVERSION of double precision matrix into SPARSE INTEGER    **
!     ******************************************************************
      subroutine double_to_sparseint(spmatrix,dense)
      implicit none

      real (dp)       , dimension(:,:), intent(in)    :: dense
      type(sparseint) , target,         intent(inout) :: spmatrix
      integer         , dimension(size(dense,1))      :: nonzero_row
      integer :: m,n,i,j,iia,nA, nIA, ntot,c, nr, nc, nemptyprev
      logical :: isfirst, updateprevious

      character(len=*), parameter ::                                   &
         fmt_ersize = "(' Error: non-2D matrix made sparse ',          &
                        ' in double_to_sparse:')"

      ! Gather matrix dimensions
      m = size(dense,1)
      n = size(dense,2)

!     Halt on non-2D matrix
      if (m*n/=size(dense)) then
         write(*,fmt_ersize)
         write(*,*)'nrow = ',m,'ncol = ',n
         stop
      endif

!     Store number of nonzero elements in each row
      nonzero_row = count(dense /= 0.0_dp, dim=2)

!     Initialize sparse matrix array
      nA = sum(nonzero_row)

      nIA = size(dense,1)
      nr  = size(dense,1)
      nc  = size(dense,2)

!     Initialise sparse matrix
      call allocate(nr, nc, nA, spmatrix)

!     Store data into sparse format
      ntot = 0
      iia  = 0
      nemptyprev = 0

      row: do i = 1, spmatrix%nr
           isfirst = .true.
           c = nonzero_row(i)

           emptyrow: if (c==0) then ! Special treatment for empty row!

              updateprevious = .true.
              nemptyprev = nemptyprev + 1

           else

              column: do j = 1, spmatrix%nc

                 if (dense(i,j) /= 0.e0_dp) then! .or. i==j) then
                    ntot = ntot + 1
                    spmatrix%A (ntot) = int(dense(i,j)+tiny(0.e0_dp))
                    spmatrix%JA(ntot) = j

                    if (isfirst) then
                       iia     = i
                       isfirst = .false.
                       spmatrix%IA(iia:iia+1) = ntot
                       if (updateprevious) then
                           updateprevious = .false.
                           spmatrix%IA(iia-nemptyprev:iia-1) = ntot
                           nemptyprev         = 0
                       endif
                    endif
                 endif
              end do column
           endif emptyrow
      end do row

!     Storing last column length
      if (.not.updateprevious) then
         spmatrix%IA(m+1) = spmatrix%IA(m) + nonzero_row(m)
      else
         spmatrix%IA(m+1-nemptyprev:m+1) = ntot+1
      endif

!     Correct storage check
      if (nA/=ntot.or.ntot/=spmatrix%IA(m+1)-1) then
         write(*,*)'sparse storage error: ntot     = ',ntot
         write(*,*)'                      nA       = ',nA
         write(*,*)'                      IA(nr+1) = ',spmatrix%IA(m+1)
         write(*,*)'matrix rows:   nr = ',spmatrix%nr
         write(*,*)'matrix cols:   nc = ',spmatrix%nc
         stop
      endif

      end subroutine double_to_sparseint


!     ******************************************************************
!     ** ZEROS = subroutine for filling existing matrices with zeros  **
!     **         NB sparsity structure is destroyed!                  **
!     ******************************************************************

      subroutine zeros_sparse(A)
      implicit none

      type (sparse), intent(inout) :: A

      ! Do nothing if matrix is not allocated
      if (.not.allocated(A)) return

      A%A (1:size(A%A) ) = 0.0_dp
      A%IA(1:size(A%IA)) = 0
      A%JA(1:size(A%JA)) = 0

      end subroutine zeros_sparse

      subroutine zeros_sparseint(A)
      implicit none

      type (sparseint), intent(inout) :: A

      ! Do nothing if matrix is not allocated
      if (.not.allocated(A)) return

      A%A (1:size(A%A) ) = 0
      A%IA(1:size(A%IA)) = 0
      A%JA(1:size(A%JA)) = 0

      end subroutine zeros_sparseint

      subroutine zeros_sparse_ordered(A)
      implicit none

      type (sparse_ordered), intent(inout) :: A

      ! Do nothing if matrix is not allocated
      if (.not.allocated(A)) return

      A%A (1:size(A%A) ) = 0
      A%IA(1:size(A%IA)) = 0
      A%JA(1:size(A%JA)) = 0

      end subroutine zeros_sparse_ordered


!     ******************************************************************
!     ** MATRIX MULTIPLICATION                                        **
!     **                                                              **
!     ******************************************************************

!     ************************************************************
!     ** Sparse matrix - dense vector multiplication, C = A*b   **
!     ************************************************************
      function sparse_matmul(A,b) result(C)
      implicit none

      type(sparse),                      intent(in) :: A
      real (dp)       , dimension(A%nc), intent(in) :: b
      real (dp)       , dimension(A%nr)             :: C

      integer, parameter     :: maxcache = 32000, cacheline = 4
      integer                :: i, ii, j, k, &
                                brow, blkrows, blkrowb, blkrowe

      if (.not.allocated(A)) then
         write(*,*)' Unassociated sparse array in sparse_matmul '
         call print_sparse_details(A,'input matrix')
         write(*,*)' dense  matrix: ',size(b)
         stop
      endif

      ! Calculate number of block rows, based on CACHELINE:
      brow = cacheline
      blkrows = A%nr/brow ! division among integers!
      if ( blkrows*brow /= A%nr ) blkrows = blkrows + 1

      ! Loop through blkrows block rows:
      loop_block_rows: do i = 1 , blkrows
         blkrowb = (i-1)*brow + 1
         blkrowe = blkrowb + brow - 1
         if (blkrowe >= A%nr ) blkrowe = A%nr

         !Now, loop through the only column rhscols block of c & b:
         !Loop through the brow rows in this block:
         loop_block_columns: do j = blkrowb,blkrowe


             ! Sparse array integer dot product
             c(j) = 0.e0_dp
             ddoti: do ii = A%IA(j), A%IA(j+1)-1
                 c(j) = c(j) + A%A(ii) * b(A%JA(ii))
             end do ddoti

         end do loop_block_columns

      end do loop_block_rows

      end function sparse_matmul

!     ************************************************************
!     ** Sparse matrix - dense vector multiplication, C = A*b   **
!     ************************************************************
      function sparseint_matmul(A,b) result(C)
      implicit none

      type(sparseint),                   intent(in) :: A
      real (dp)       , dimension(A%nc), intent(in) :: b
      real (dp)       , dimension(A%nr)             :: C

      integer, parameter     :: maxcache = 32000, cacheline = 8
      integer                :: i, ii, j, k, &
                                brow, blkrows, blkrowb, blkrowe

      if (.not.allocated(A)) then
         write(*,*)' Unassociated sparse array in sparseint_matmul '
         call matrix_details(A,'integer input matrix')
         write(*,*)' dense  matrix: ',size(b)
         stop
      endif

      ! Calculate number of block rows, based on CACHELINE:
      brow = cacheline
      blkrows = A%nr/brow ! division among integers!
      if ( blkrows*brow /= A%nr ) blkrows = blkrows + 1

      ! Loop through blkrows block rows:
      loop_block_rows: do i = 1 , blkrows
         blkrowb = (i-1)*brow + 1
         blkrowe = blkrowb + brow - 1
         if (blkrowe >= A%nr ) blkrowe = A%nr

         !Now, loop through the only column rhscols block of c & b:
         !Loop through the brow rows in this block:
         loop_block_columns: do j = blkrowb,blkrowe


             ! Sparse array integer dot product
             c(j) = 0.e0_dp
             ddoti: do ii = A%IA(j), A%IA(j+1)-1
                 c(j) = c(j) + A%A(ii) * b(A%JA(ii))
             end do ddoti

         end do loop_block_columns

      end do loop_block_rows

      end function sparseint_matmul


!     ******************************************************************
!     ** EQUIVALENCE -related routines: for assigning one sparse      **
!     ** matrix to another                                            **
!     ******************************************************************

!     ************************************************************
!     ** Sparse Equivalence                                     **
!     ** This routine assign one sparse matrix to another       **
!     ** sparse matrix object, initialising and assigning all   **
!     ** the arrays                                             **
!     ************************************************************

      subroutine sparse_equivalence(B,A)
      implicit none

      type(sparse), intent(in)      :: A
      type(sparse), intent(inout)   :: B

!     Preliminary check.  A matrix with valid dimensions and no stored
!     entries is still an initialized sparse matrix.
      if (A%nr<0.or.A%nc<1.or.A%n<0) then
         call matrix_details(A,'input matrix, A')
         write(*,*)'sparse_equivalence'
         write(*,*)'Error in sparse_equivalence: A not init'
         write(*,*)'A%nr=',A%nr,' A%nc=',A%nc,' A%n=',A%n
         stop
      endif

!     Initialize matrix B with A dimensions
      call allocate(A%nr, A%nc, A%n, B)

!     Assign values
      B%A  (1:A%n) = A%A(1:A%n)
      B%IA = A%IA
      B%JA (1:A%n) = A%JA(1:A%n)

      end subroutine sparse_equivalence

      subroutine sparse_ordered_to_sparse_equivalence(B,A)
      implicit none

      type(sparse_ordered), intent(in)    :: A
      type(sparse),         intent(out)   :: B

!     Preliminary check
      if (.not.(A%nr>0).or.(.not.(A%nc>0)).or.(.not.(A%n>0))) then
         write(*,*)'sparse_ordered_to_sparse_equivalence'
         write(*,*)'Error in sparse_equivalence: A not init'
         write(*,*)'A%nr=',A%nr,' A%nc=',A%nc,' A%n=',A%n
         stop
      endif

!     Simplify if the structure is the same
!      if (same_structure(A,B)) then
!         B%A = A%A
!         return
!      endif

!     Initialize matrix B with A dimensions
      call allocate(A%nr, A%nc, A%n, B)

!     Assign values
      B%A  (1:A%n) = A%A(1:A%n)
      B%IA = A%IA
      B%JA (1:A%n) = A%JA(1:A%n)

      end subroutine sparse_ordered_to_sparse_equivalence


!     ************************************************************
!     ** Sparse Equivalence                                     **
!     ** This routine assign one sparse matrix to another       **
!     ** sparse matrix object, initialising and assigning all   **
!     ** the arrays                                             **
!     ************************************************************

      subroutine sparseint_equivalence(B,A)
      implicit none

      type(sparseint), intent(in)  :: A
      type(sparseint), intent(out) :: B

!     Preliminary check
      if (.not.(A%nr>0).or.(.not.(A%nc>0)).or.(.not.(A%n>0))) then
         write(*,*)'sparseint_equivalence'
         write(*,*)'Error in sparse_equivalence: A not init'
         write(*,*)'A%nr=',A%nr,' A%nc=',A%nc,' A%n=',A%n
         stop
      endif

!     Initialize matrix B with A dimensions
      call allocate(A%nr, A%nc, A%n, B)

!     Assign values
      B%A  (1:A%n) = A%A(1:A%n)
      B%IA = A%IA
      B%JA (1:A%n) = A%JA(1:A%n)

      end subroutine sparseint_equivalence

      subroutine sparse_to_sparseint_equivalence(B,A)
      implicit none

      type(sparse   ), intent(in)  :: A
      type(sparseint), intent(out) :: B

!     Preliminary check
      if (.not.(A%nr>0).or.(.not.(A%nc>0)).or.(.not.(A%n>0))) then
         write(*,*)'sparse_to_sparseint_equivalence'
         write(*,*)'Error in sparse_equivalence: A not init'
         write(*,*)'A%nr=',A%nr,' A%nc=',A%nc,' A%n=',A%n
         stop
      endif

!     Initialize matrix B with A dimensions
      call allocate(A%nr, A%nc, A%n, B)

!     Assign values
      B%A  (1:A%n) = int(A%A(1:A%n))
      B%IA = A%IA
      B%JA (1:A%n) = A%JA(1:A%n)

      end subroutine sparse_to_sparseint_equivalence

      subroutine sparseint_to_sparse_equivalence(B,A)
      implicit none

      type(sparseint), intent(in)  :: A
      type(sparse   ), intent(out) :: B

!     Preliminary check
      if (.not.(A%nr>0).or.(.not.(A%nc>0)).or.(.not.(A%n>0))) then
         write(*,*)'sparseint_to_sparse_equivalence'
         write(*,*)'Error in sparse_equivalence: A not init'
         write(*,*)'A%nr=',A%nr,' A%nc=',A%nc,' A%n=',A%n
         stop
      endif

!     Initialize matrix B with A dimensions
      call allocate(A%nr, A%nc, A%n, B)

!     Assign values
      B%A  (1:B%n) = real(A%A(1:B%n), dp)
      B%IA = A%IA
      B%JA (1:B%n) = A%JA(1:B%n)

      end subroutine sparseint_to_sparse_equivalence

!     ************************************************************
!     ** Sparse Ordered matrix assignment                       **
!     ** This routine assign one sparse matrix to the ordered   **
!     ** form, i.e. adds the optimal ordering information and   **
!     ** initialises working arrays for the solution of linear  **
!     ** systems                                                **
!     ************************************************************

      subroutine sparse_to_sparse_ordered_equivalence(B,A)
      use utilities, only: extend
      implicit none

      type(sparse),         intent(in)            :: A
      type(sparse_ordered), intent(inout), target :: B

      type(sparse)                                :: tmpA
      integer                                     :: j
      logical                                     :: full_assignment


!     Preliminary check
      if (.not.(A%nr>0).or.(.not.(A%nc>0)).or.(.not.(A%n>0))) then
         write(*,*)'sparse_to_sparse_ordered_equivalence'
         write(*,*)'Error in sparse_equivalence: A not init'
         write(*,*)'A%nr=',A%nr,' A%nc=',A%nc,' A%n=',A%n
         stop
      endif

      if (A%nr /= A%nc) then
         write(*,*)'sparse_ordered matrices have to be square'
         write(*,*)'A%nr= ',A%nr, ', A%nc=',A%nc
         stop
      endif

      ! Check if a full assignment is needed
      full_assignment = .false.

      if (.not.allocated(B)) full_assignment = .true.
      if (.not.allocated(A)) full_assignment = .true.

      if (allocated(B) .and. allocated(A)) then

         if (.not.same_structure(B,A)) full_assignment = .true.

      endif


      equal_structure: if (.not.full_assignment) then

         B%A  (1:B%n) = A%A (1:A%n)
         B%JA (1:B%n) = A%JA(1:A%n)
         B%numerically_factorized = .false.

         return

      else

!        Allocate matrix B with A dimensions
         call extend(B%A,  A%n       )
         call extend(B%IA, A%nr+1    )
         call extend(B%JA, A%n       )

!        Assign values
         B%n  = A%n
         B%nr = A%nr
         B%nc = A%nc
         B%A  (1:B%n) = A%A (1:B%n)
         B%IA         = A%IA
         B%JA (1:B%n) = A%JA(1:B%n)

!        Allocate pointers
!        Associate pointers for CSR format
         B%val  => B%A  (1:B%n)
         B%col  => B%JA (1:B%n)
         B%ptrB => B%IA (1:B%nr  )
         B%ptrE => B%IA (2:B%nr+1)

!        Allocate linear systems-related data: permutations
         call extend(B%perm    , B%nr)
         call extend(B%inv_perm, B%nr)

!        Compute optimal permutations and inverse array
         B%perm             = optimal_ordering(A)
         B%inv_perm(B%perm) = [(j,j=1,B%nr)]
         B%is_ordered       = .true.

         B%lR = 0
         B%lI = 0

         ! Compute symbolic LU factorization
         call sparse_ordered_symbolic_factorization(B)

         ! Deallocate temporary matrix
!         call deallocate(tmpA)

      endif equal_structure


      end subroutine sparse_to_sparse_ordered_equivalence

!     ******************************************************************
!     ** COMPARISON-RELATED ROUTINES                                  **
!     ******************************************************************
      function sparse_ordered_comparison_sparse(A,B) result(equal)

      use utilities, only: sort

      implicit none

      type(sparse_ordered), intent(in)     :: A
      type(sparse),         intent(in)     :: B
      integer,              dimension(A%n), target  :: AtmpJA
      integer,              dimension(B%n), target  :: BtmpJA
      logical                              :: equal
      integer                              :: i, j, jj
      integer,              dimension(:),   pointer :: Acol, Bcol

      ! Initialize as not equal
      equal = .false.

      if (A%nr /= B%nr) return
      if (A%nc /= B%nc) return
      if (A%n  /= B%n ) return
      if (size(A%IA)/=size(B%IA)) return
      if (any(A%A(1:A%n)  /= B%A(1:A%n)) ) return
      if (size(A%IA) == size(B%IA)) then
         if (any(A%IA /= B%IA)) return
      endif

      AtmpJA = A%JA(1:A%n)
      BtmpJA = B%JA(1:B%n)

      ! Column indices could have the same pattern, but
      ! different ordering in the JA array: so,
      ! sort both arrays and compare the sorted ones
      if (A%n == B%n) then
      if (any(A%JA(1:A%n)/=B%JA(1:B%n))) then
      rows: do i = 1, A%nr
         nonzero: if (A%IA(i+1) > A%IA(i)) then
            Acol => AtmpJA(A%IA(i):A%IA(i+1)-1)
            Bcol => BtmpJA(B%IA(i):B%IA(i+1)-1)

            ! Sort column indices
            call sort(Bcol)
            call sort(Acol)

            ! Exit on different value
            if (any(Bcol /= Acol))return

         endif nonzero
      end do rows
      endif
      endif

      ! All the checks passed; matrices are equal
      equal = .true.

      end function sparse_ordered_comparison_sparse

      function sparse_comparison_sparse_ordered(A,B) result(equal)
      implicit none

      type(sparse),                 intent(in) :: A
      type(sparse_ordered),         intent(in) :: B
      logical                                  :: equal

      ! Initialize as not equal
      equal = .false.

      if (A%nr /= B%nr) return
      if (A%nc /= B%nc) return
      if (A%n  /= B%n ) return

      if (any(A%A (1:A%n) /= B%A (1:A%n))) return
      if (any(A%IA        /= B%IA       )) return
      if (any(A%JA(1:A%n) /= B%JA(1:A%n))) return

      ! All the checks passed; matrices are equal
      equal = .true.

      return
      end function sparse_comparison_sparse_ordered


      function sparse_comparison_sparse(A,B) result(equal)
      implicit none

      type(sparse), intent(in) :: A, B
      logical                  :: equal


      ! Initialize as not equal
      equal = .false.

      if (A%nr /= B%nr) return
      if (A%nc /= B%nc) return
      if (A%n  /= B%n ) return

!      if (size(A%IA)/=size(B%IA)) return
!      if (size(A%JA)/=size(B%JA)) return
!      if (size(A%A)/=size(B%A)) return

      if (any(A%A (1:A%n) /= B%A (1:A%n))) return
      if (any(A%IA        /= B%IA       )) return
      if (any(A%JA(1:A%n) /= B%JA(1:A%n))) return

      ! All the checks passed; matrices are equal
      equal = .true.

      return
      end function sparse_comparison_sparse

      function sparse_structurecomp_sparse(A,B) result(equal)

      use utilities, only: sort

      implicit none

      type(sparse),         intent(in)     :: A, B
      integer,              dimension(A%n), target  :: AtmpJA
      integer,              dimension(B%n), target  :: BtmpJA
      logical                              :: equal
      integer                              :: i, j, jj
      integer,              dimension(:),   pointer :: Acol, Bcol

      ! Initialize as not equal
      equal = .false.

      if (.not. allocated(A)) return
      if (.not. allocated(B)) return
      if (A%nr /= B%nr) return
      if (A%nc /= B%nc) return
      if (A%n  /= B%n ) return
      if (size(A%IA)/=size(B%IA)) return
!      if (size(A%JA)/=size(B%JA)) return
!      if (size(A%A)/=size(B%A)) return
      if (size(A%IA) == size(B%IA)) then
         if (any(A%IA /= B%IA)) return
      endif

      AtmpJA = A%JA(1:A%n)
      BtmpJA = B%JA(1:B%n)

      ! Column indices could have the same pattern, but
      ! different ordering in the JA array: so,
      ! sort both arrays and compare the sorted ones
      if (A%n == B%n) then
      if (any(A%JA(1:A%n)/=B%JA(1:B%n))) then
      rows: do i = 1, A%nr
         nonzero: if (A%IA(i+1) > A%IA(i)) then
            Acol => AtmpJA(A%IA(i):A%IA(i+1)-1)
            Bcol => BtmpJA(B%IA(i):B%IA(i+1)-1)

            ! Sort column indices
            call sort(Bcol)
            call sort(Acol)

            ! Exit on different value
            if (any(Bcol /= Acol)) return

         endif nonzero
      end do rows
      endif
      endif

      ! All the checks passed; matrices are equal
      equal = .true.

      return

      end function sparse_structurecomp_sparse

      function sparse_ordered_structurecomp_sparse(A,B) result(equal)

      use utilities, only: sort

      implicit none

      type(sparse_ordered), intent(in)     :: A
      type(sparse),         intent(in)     :: B
      integer,              dimension(A%n), target  :: AtmpJA
      integer,              dimension(B%n), target  :: BtmpJA
      logical                              :: equal
      integer                              :: i, j, jj
      integer,              dimension(:),   pointer :: Acol, Bcol

      ! Initialize as not equal
      equal = .false.

      if (.not. allocated(A)) return
      if (.not. allocated(B)) return
      if (A%nr /= B%nr) return
      if (A%nc /= B%nc) return
      if (A%n  /= B%n ) return
      if (size(A%IA)/=size(B%IA)) return
!      if (size(A%JA)/=size(B%JA)) return
!      if (size(A%A)/=size(B%A)) return
      if (size(A%IA) == size(B%IA)) then
         if (any(A%IA /= B%IA)) return
      endif

      AtmpJA = A%JA(1:A%n)
      BtmpJA = B%JA(1:B%n)

      ! Column indices could have the same pattern, but
      ! different ordering in the JA array: so,
      ! sort both arrays and compare the sorted ones
      if (A%n == B%n) then
      if (any(A%JA(1:A%n)/=B%JA(1:B%n))) then
      rows: do i = 1, A%nr
         nonzero: if (A%IA(i+1) > A%IA(i)) then
            Acol => AtmpJA(A%IA(i):A%IA(i+1)-1)
            Bcol => BtmpJA(B%IA(i):B%IA(i+1)-1)

            ! Sort column indices
            call sort(Bcol)
            call sort(Acol)

            ! Exit on different value
            if (any(Bcol /= Acol)) return

         endif nonzero
      end do rows
      endif
      endif

      ! All the checks passed; matrices are equal
      equal = .true.

      return

      end function sparse_ordered_structurecomp_sparse

      function sparse_structurecomp_sparse_ordered(A,B) result(equal)

      use utilities, only: sort

      implicit none

      type(sparse),         intent(in)     :: A
      type(sparse_ordered), intent(in)     :: B
      integer,              dimension(A%n), target  :: AtmpJA
      integer,              dimension(B%n), target  :: BtmpJA
      logical                              :: equal
      integer                              :: i, j, jj
      integer,              dimension(:),   pointer :: Acol, Bcol

      ! Initialize as not equal
      equal = .false.

      if (.not. allocated(A)) return
      if (.not. allocated(B)) return
      if (A%nr /= B%nr) return
      if (A%nc /= B%nc) return
      if (A%n  /= B%n ) return
      if (size(A%IA)/=size(B%IA)) return
!      if (size(A%JA)/=size(B%JA)) return
!      if (size(A%A)/=size(B%A)) return
      if (size(A%IA) == size(B%IA)) then
         if (any(A%IA /= B%IA)) return
      endif

      AtmpJA = A%JA(1:A%n)
      BtmpJA = B%JA(1:B%n)

      ! Column indices could have the same pattern, but
      ! different ordering in the JA array: so,
      ! sort both arrays and compare the sorted ones
      if (A%n == B%n) then
      if (any(A%JA(1:A%n)/=B%JA(1:B%n))) then
      rows: do i = 1, A%nr
         nonzero: if (A%IA(i+1) > A%IA(i)) then
            Acol => AtmpJA(A%IA(i):A%IA(i+1)-1)
            Bcol => BtmpJA(B%IA(i):B%IA(i+1)-1)

            ! Sort column indices
            call sort(Bcol)
            call sort(Acol)

            ! Exit on different value
            if (any(Bcol /= Acol)) return

         endif nonzero
      end do rows
      endif
      endif

      ! All the checks passed; matrices are equal
      equal = .true.

      return
      end function sparse_structurecomp_sparse_ordered

!     ******************************************************************
!     ** SYMBOLIC  LU FACTORIZATION of a sparse_ordered matrix        **
!     ******************************************************************
      subroutine sparse_ordered_symbolic_factorization(B)
      use utilities, only: extend
      implicit none

      type(sparse_ordered), intent(inout) :: B
      integer, parameter                  :: safety = 1
      integer                             :: j, error_flag, free_space
      real (dp)       , dimension(B%nr)   :: tmp

      if (B%n == 0 .or. B%nr == 0) return

!     First hypothesis on required workspace
      if (B%lI==0) B%lR = safety * ( 9*(B%nr+2) + 2 * B%n   )
      if (B%lI==0) B%lI = real_to_int_length * B%lR
      j    = 0

!     Allocate WORKSPACE for linear systems computations
      find_space: do while (j < 100)

        j = j + 1

        ! Tentative allocation
        B%numerically_factorized = .false.
        call extend(B%real_space, B%lR)
        call extend(B% int_space, B%lI)

        tmp          = 0.e0_dp
        B%real_space = 0.e0_dp
        B%int_space  = 0

        ! Run symbolic factorization routine
        error_flag = 0
        free_space = 0

        call cdrv(B%nr, B%perm, B%perm, B%inv_perm, B%IA, B%col, B%val,&
                  tmp, tmp, B%lR, B%int_space, B%real_space,           &
                  free_space, 5, error_flag)

        if ( error_flag == 0 .or. j>10 ) exit find_space

        ! Update dimensions for new attempt
        B%lR = B%lR - free_space
        B%lI = B%lR * real_to_int_length

        if (error_flag == 10*B%nr + 1) then
           write(*,*)'insufficient space in CDRV',free_space,'storage',B%lR
           write(*,*)'possibly meaning singular matrix'
           write(*,*)'error_flag',error_flag
        endif

      end do find_space

      ! Logical flag
      B%symbolically_factorized = .true.

      ! Halt on error after 100 iterations
      if (error_flag /= 0) then
         write(*,*)'Persisting error in sparse_ordered initialisation'
         write(*,*)'CDRV error flag   = ',error_flag
         write(*,*)'Missing space     = ',free_space
         write(*,*)'Matrix dimensions = ',B%nr
         stop
      endif

      return
      end subroutine sparse_ordered_symbolic_factorization


!     ******************************************************************
!     ** NUMERICAL LU FACTORIZATION of a sparse_ordered matrix        **
!     ******************************************************************
      subroutine sparse_ordered_numerical_factorization(A)
      use utilities, only: change_size

      implicit none
      type(sparse_ordered), intent(inout) :: A

      real (dp)       , dimension(A%nr)   :: tmpx
      real (dp)       , dimension(:), allocatable :: rwork
      integer,          dimension(:), allocatable :: iwork

      integer :: free_space, error_flag, j

      j = 0

      ! Pursue factorization again if space is not enough
      attempt_loop: do while (j < 100)

        j = j + 1

        ! Initialise free space and error_flag arrays
        free_space = 0
        error_flag = 0
        tmpx       = 0.e0_dp

        call cdrv(A%nr, A%perm, A%perm, A%inv_perm, A%IA, A%JA(1:A%n),     &
                  A%A(1:A%n), tmpx, tmpx, A%lR, A%int_space, A%real_space, &
                  free_space, 2, error_flag)

!        write(*,*)j,free_space,error_flag

        if (error_flag == 0 .and. free_space>=0) exit attempt_loop

        ! Error flags for insufficient space
        if (error_flag ==  7 * A%nr + 1 .or. &
            error_flag == 10 * A%nr + 1 .or. &
            error_flag ==  4 * A%nr + 1 .or. &
            free_space < 0                   )     then

           A%lR = A%lR - free_space
           A%lI = A%lR * real_to_int_length

           call change_size(A%real_space,A%lR)
           call change_size(A% int_space,A%lI)

           ! Symbolic factorization is needed again if the array
           ! dimensions are changed
           call sparse_ordered_symbolic_factorization(A)

        endif

      end do attempt_loop

!      write(*,*)'stop',j,error_flag,free_space,A%lr,A%li,A%nr

      if (error_flag/=0 .or. free_space < 0) then
         call matrix_details(A,'matrix to be factorized')
         call print_sparse_ordered_to_file(A,'LU_error_matrix.dat')
         write(*,*)'Sparse numerical factorization:'
         write(*,*)'cdrv exited with error flag ',error_flag
         write(*,*)'working storage free space: ',free_space
         write(*,*)'working storage provided  : ',A%lR
         stop
      endif

      ! Set numerical factorization flag
      A%numerically_factorized = .true.

      end subroutine sparse_ordered_numerical_factorization

!     ******************************************************************
!     ** SOLUTION OF A SPARSE LINEAR SYSTEM                           **
!     ******************************************************************
      function sparse_ordered_linear_system(A,b) result(x)
      implicit none

      type(sparse_ordered),                intent(in)    :: A
      real (dp)       , dimension(:),      intent(in)       :: b
      real (dp)       , dimension(size(b))                  :: x

      integer :: free_space, error_flag, j

      ! Preliminary checks on dimensions
      if (size(b) /= A%nr) then
         write(*,*)'Linear system: rank mismatch '
         write(*,*)'Coefficients array: ',size(b)
         write(*,*)'System matrix     : ',A%nr
         stop
      end if

      if (.not.A%is_ordered) then
        write(*,*)'Unordered matrix in linear system'
        stop
      endif

      if (.not.A%symbolically_factorized) then
        write(*,*)'Unfactorized matrix in linear system'
        stop
      endif

      if (.not.A%numerically_factorized) then
          write(*,*)' numerically Unfactorized matrix'
          stop
      endif


      ! Initialise free space and error_flag arrays
      free_space = 0
      error_flag = 0
      x          = 0.e0_dp

!      call matrix_details(A)
!      pause

      call cdrv(A%nr, A%perm, A%perm, A%inv_perm, A%IA, A%JA(1:A%n), &
                A%A(1:A%n), b, x, A%lR, A%int_space, A%real_space,   &
                free_space, 3, error_flag)

!      write(*,*)free_space, error_flag


      if (error_flag/=0 .or. free_space<0) then
         write(*,*)'Sparse linear system solution error.'
         write(*,*)'cdrv exited with error flag ',error_flag
         write(*,*)'working storage free space: ',free_space
         write(*,*)'working storage provided  : ',A%lR
         stop
      endif


      end function sparse_ordered_linear_system

!     ******************************************************************
      function sparse_ordered_linear_system_double(A,b) result(x)
      implicit none

      type(sparse_ordered),                intent(in)    :: A
      double precision, dimension(:),      intent(in)    :: b
      double precision, dimension(size(b))               :: x

      real (dp),        dimension(size(b))               :: tmpx, tmpb

      integer :: free_space, error_flag, j

      tmpb = real(b, dp)
      tmpx = sparse_ordered_linear_system(A,tmpb)
      x    = dble(tmpx )


      end function sparse_ordered_linear_system_double


!     ******************************************************************
      function sparse_linear_system(A,b) result(x)
      implicit none

      type(sparse),                        intent(in) :: A
      real (dp)       , dimension(:),      intent(in) :: b
      real (dp)       , dimension(size(b))            :: x

      type(sparse_ordered)                            :: tmp

      ! Assign and reorder matrix
      tmp = A

      ! Numerically factorize tmp
      call sparse_ordered_numerical_factorization(tmp)

      ! Solve linear system
      x = sparse_ordered_linear_system(tmp,b)

      ! Delete temporary variable
      call deallocate(tmp)

      end function sparse_linear_system

!     ******************************************************************
      function sparse_linear_system_double(A,b) result(x)
      implicit none

      type(sparse),                        intent(in) :: A
      double precision, dimension(:),      intent(in) :: b
      double precision, dimension(size(b))            :: x

      type(sparse_ordered)                            :: tmp
      real(dp)                                        :: tmpx

      ! Assign and reorder matrix
      tmp = A

      ! Numerically factorize tmp
      call sparse_ordered_numerical_factorization(tmp)

      ! Solve linear system
      x = sparse_ordered_linear_system(tmp,real(b, kind=dp))
!      x    = dble(tmpx)

      ! Delete temporary variable
      call deallocate(tmp)

      end function sparse_linear_system_double

!     ******************************************************************
!     ** ELEMENTARY OPERATIONS                                        **
!     ******************************************************************

!     ******************************************************************
!     ** Performs the sum of two sparse matrices                      **
!     ** C = A + B                                                    **
!     ******************************************************************
      function sparse_plus_sparse(A,B) result(C)
      implicit none

      type(sparse)    , intent(in) :: A, B
      type(sparse)                 :: C

      integer, dimension(A%nc)     :: itmp
      integer                      :: cur_element, i, k, ka, kb, jpos
      logical                      :: noA, noB

!     Preliminary check on matrix dimensions
      if (A%nr /= B%nr) then
         write(*,*)'Sparse sum: wrong number of rows'
         call matrix_details(A,'left in sum' )
         call matrix_details(B,'right in sum')
         stop
      endif

      if (A%nc /= B%nc) then
         write(*,*)'Sparse addition: wrong number of columns'
         write(*,*)'A=',A%nc,' B=',B%nc
         error stop 1
      endif

!      noA = .false.!.not.allocated(A)
!      noB = .false.!.not.allocated(B)
!
!      full_sum: if (.not.(noA.or.noB)) then

!     Initialise temporary array
      cur_element = A%n + B%n! - overlaps(A,B)
      call allocate(A%nr, A%nc, cur_element, C)

!     Initialise number of nonzero elements in matrix C
      cur_element  = 0

!     Initialise first element
      C%IA(1) = 1

!     Initialise column working array
      itmp(1:A%nc) = 0

      rows: do i = 1, A%nr

             
            thisrowA: do ka = A%IA(i), A%IA(i+1) - 1
                        if (A%A(ka) == 0.e0_dp) cycle thisrowA
                        cur_element = cur_element + 1
                        C%JA(cur_element) = A%JA(ka)
                        C%A (cur_element) = A%A (ka)
                        itmp(A%JA(ka)) = cur_element
                      end do thisrowA

            thisrowB: do kb = B%IA(i), B%IA(i+1) - 1
                        if (B%A(kb)==0.e0_dp) cycle thisrowB
                        jpos = itmp(B%JA(kb))
                        newposition: if (jpos == 0) then
                          cur_element = cur_element + 1
                          C%JA(cur_element) = B%JA(kb)
                          C%A(cur_element)  = B%A(kb)
                          itmp(B%JA(kb))      = cur_element
                        else
                          C%A(jpos)  = C%A(jpos) + B%A(kb)

                          delete_empty: if (C%A(jpos)==0.e0_dp) then
                            cur_element = cur_element - 1
                            C% A(jpos:cur_element) = C% A(jpos+1:cur_element+1)
                            C%JA(jpos:cur_element) = C%JA(jpos+1:cur_element+1)
                            itmp(B%JA(kb)) = 0
                            itmp(C%JA(jpos:cur_element)) = itmp(C%JA(jpos:cur_element)) -1
                          endif delete_empty

                        endif newposition
                      end do thisrowB

            ! Clean working array, prepare for next row
            clean_itmp: do k = C%IA(i), cur_element
                           itmp(C%JA(k)) = 0
                        end do clean_itmp

            ! Update number of elements in current row
            C%IA(i+1) = cur_element+1

      end do rows

      C%n = cur_element
!      call sparse_compress(C)

!      else
!
!        if (noA.and.noB) then
!           call deallocate(C)
!        elseif (noA) then
!           call sparse_equivalence(C,B)
!        elseif (noB) then
!           call sparse_equivalence(C,A)
!        endif
!
!      endif full_sum

      end function sparse_plus_sparse

!     ******************************************************************
!     ** Generate a identity matrix of order n in sparse form         **
!     ** beta is a real factor that multiplies the matrix             **
!     ******************************************************************
      function identity(n,beta) result(I)
      implicit none

      integer, intent(in)          :: n
      real (dp)       , optional, intent(in) :: beta
      type(sparse)                 :: I
      integer                      :: j

      if (n<2) then
         write(*,*)'An identity matrix must have n>1; n =',n
         stop
      endif

!     Allocate space for the sparse matrix
      call allocate(n, n, n, I)

      if (present(beta)) then
         I%A(1:n)  = beta
      else
         I%A(1:n)  = 1.e0_dp
      endif
      I%IA      = [(j,j=1,n+1)]
      I%JA(1:n) = [(j,j=1,n)]

      end function identity



!     ******************************************************************
!     ** Negative of a sparse matrix                                  **
!     ** B = -A                                                       **
!     ******************************************************************
      function sparse_negative(A) result(B)
      implicit none

      type(sparse)    , intent(in) :: A
      type(sparse)                 :: B

      call allocate(A%nr, A%nc, A%n, B)

      B%A  (1:B%n) = -A%A(1:B%n)
      B%IA = A%IA
      B%JA (1:B%n) = A%JA(1:B%n)

      end function sparse_negative



!     ******************************************************************
!     ** Performs the product of every element in a sparse matrix     **
!     ** by a real number: C = a*B                                    **
!     ******************************************************************
      function double_dot_sparse(a,B) result(C)
      implicit none

      real (dp)       , intent(in) :: a
      type(sparse)    , intent(in) :: B
      type(sparse)                 :: C


!     Preserve dimensions even when the numerical result is zero.
!     Returning an uninitialized object loses its column count and makes
!     a later sparse sum fail for absent optional chemistry terms.
      call allocate(B%nr, B%nc, B%n, C)
      C%IA = B%IA
      C%JA (1:C%n) = B%JA(1:C%n)
      C%A  (1:C%n) = B%A (1:C%n) * a

      end function double_dot_sparse

      function sparse_dot_double(B,a) result(C)
      implicit none

      real (dp)       , intent(in) :: a
      type(sparse)    , intent(in) :: B
      type(sparse)                 :: C

      if (a == 0.e0_dp) then
         call deallocate(C)
      else
         ! Initialize matrix
         call allocate(B%nr, B%nc, B%n, C)

         ! Perform multiplication
         C%IA = B%IA
         C%JA (1:C%n) = B%JA(1:C%n)
         C%A  (1:C%n) = B%A (1:C%n) * a
      endif

      end function sparse_dot_double


!     ******************************************************************
!     ** Performs the difference of two sparse matrices               **
!     ** C = A - B                                                    **
!     ******************************************************************
      function sparse_minus_sparse(A,B) result(C)
      implicit none

      type(sparse)    , intent(in) :: A, B
      type(sparse)                 :: C

      integer, dimension(A%nc)     :: itmp
      integer                      :: cur_element, i, k, ka, kb, jpos
      logical                      :: noA, noB

!     Preliminary check on matrix dimensions
      if (A%nr /= B%nr) then
         write(*,*)'Sparse difference: wrong number of rows'
         call matrix_details(A,'left in difference' )
         call matrix_details(B,'right in difference')
         stop
      endif

      if (A%nc /= B%nc) then
         write(*,*)'Sparse difference: wrong number of columns'
         write(*,*)'A=',A%nc,' B=',B%nc
         error stop 1
      endif

!      noA = .false.!.not.allocated(A)
!      noB = .false.!.not.allocated(B)
!
!      full_sum: if (.not.(noA.or.noB)) then

!     Initialise temporary array
      cur_element = A%n + B%n! - overlaps(A,B)
      call allocate(A%nr, A%nc, cur_element, C)

!     Initialise number of nonzero elements in matrix C
      cur_element  = 0

!     Initialise first element
      C%IA(1) = 1

!     Initialise column working array
      itmp(1:A%nc) = 0

      rows: do i = 1, A%nr
            thisrowA: do ka = A%IA(i), A%IA(i+1) - 1
!                        if (A%A(ka)==0.e0_dp) cycle thisrowA
                        cur_element = cur_element + 1
                        C%JA(cur_element) = A%JA(ka)
                        C%A (cur_element) = A%A (ka)
                        itmp(A%JA(ka)) = cur_element
                      end do thisrowA
            thisrowB: do kb = B%IA(i), B%IA(i+1) - 1
!                        if (B%A(kb)==0.e0_dp) cycle thisrowB
                        jpos = itmp(B%JA(kb))
                        newposition: if (jpos == 0) then
                          cur_element = cur_element + 1
                          C%JA(cur_element) = B%JA(kb)
                          C%A(cur_element)  = - B%A(kb)
                          itmp(B%JA(kb))    = cur_element
                        else
                          C%A(jpos)  = C%A(jpos) - B%A(kb)

                          delete_empty: if (C%A(jpos)==0.e0_dp) then
                            cur_element = cur_element - 1
                            C% A(jpos:cur_element) = C% A(jpos+1:cur_element+1)
                            C%JA(jpos:cur_element) = C%JA(jpos+1:cur_element+1)
                            itmp(B%JA(kb)) = 0
                            itmp(C%JA(jpos:cur_element)) = itmp(C%JA(jpos:cur_element)) -1
                          endif delete_empty


                        endif newposition
                      end do thisrowB
            ! Clean working array, prepare for next row
            clean_itmp: do k = C%IA(i), cur_element
                           itmp(C%JA(k)) = 0
                        end do clean_itmp

            ! Update number of elements in current row
            C%IA(i+1) = cur_element+1

      end do rows


      C%n = cur_element
!     Clean zero elements
!      call sparse_compress(C)

!      else
!
!        if (noA.and.noB) then
!           call deallocate(C)
!        elseif (noA) then
!           call sparse_equivalence(C,B)
!        elseif (noB) then
!           call sparse_equivalence(C,A)
!        endif
!
!      endif full_sum

!     ******************************************************************

      end function sparse_minus_sparse

!     ************************************************************
!     ** Compress sparse matrix by removing zero elements       **
!     ************************************************************

      subroutine sparse_compress(A)
      use utilities, only: change_size
      implicit none

      type(sparse), intent(inout), target        :: A
      integer :: i, j, r0, rf, removed, n_null
      logical,          dimension(A%n)           :: is_null

      if (.not.(A%nr>0.and.A%nc>0.and.A%n>0)) then
         write(*,*)'sparse_compress: wrong matrix dim:'
         write(*,*)A%nr, A%nc, A%n
         stop
      endif

      ! Logical positions of empty values
      is_null = (A%A(1:A%n) == 0.e0_dp)

      ! Number of null elements
      n_null  = count(is_null)

      ! If nothing to do
      if (n_null == 0) return

      ! Check the lines where elements have to be removed
      ! rem_per_line = 0
      removed = 0
      do i = 1, A%nr
         r0 = A%IA(i)
         rf = A%IA(i+1)-1
!         rem_per_line(i) = count(is_null(r0:rf))
         A%IA(i) = A%IA(i) - removed
         removed = removed + count(is_null(r0:rf))!rem_per_line(i)
      end do
      A%IA(A%nr+1) = A%IA(A%nr+1) - removed


!     Update number of nonzero elements
      A%A (1:) = pack(A%A (1:A%n), .not.is_null)
      A%JA(1:) = pack(A%JA(1:A%n), .not.is_null)


!     Update sparse matrix arrays
      A%n  = A%n - n_null


!     Attempt not to reallocate, but leave blank trailing space
!      call change_size(A%A , A%n)
!      call change_size(A%JA, A%n)

!      deallocate(A%A)
!      deallocate(A%JA)
!
!      allocate(A%A(A%n), A%JA(A%n))
!      A%A      = tmpA(1:A%n)
!      A%JA     = tmpJA(1:A%n)

!      removed = 0
!      do j = 1, A%nr
!         removed = removed + rem_per_line(j)
!         A%IA(j+1) = A%IA(j+1) - removed
!      end do

!     Associate pointers for CSR format
      A%val  => A%A  (1:A%n)
      A%col  => A%JA (1:A%n)

      end subroutine sparse_compress


!     ******************************************************************
!     ** SOLUTION OF LINEAR SYSTEMS - related subroutines             **
!     ******************************************************************

!     ******************************************************************
!     ** Compute optimal ordering for the sparse matrix rows/columns  **
!     ** that reduces fill-in during LU factorization                 **
!     ******************************************************************
      function optimal_ordering(matrix) result(permutations)
      implicit none

      type(sparse),                  intent(in)  :: matrix
      integer, dimension(matrix%nr)              :: permutations
      integer, dimension(matrix%nr)              :: inv_permutations

      integer                            :: lwork, error_flag
      integer, dimension(:), allocatable :: iwork

      ! Allocate temporary storage for the ordering routine
      ! (minimum is 3*nrows + 4 * number of elements in the upper
      ! triangular part of the matrix. Put all size for safety)
      lwork = 3 * matrix%nr + 4 * matrix%n

      allocate(iwork(lwork))
      iwork = 0

      ! Initialise permutations arrays
      permutations     = 0
      inv_permutations = 0

      ! Call driver to the Yale optimal ordering routine
      call odrv(matrix%nr, matrix%IA, matrix%JA(1:matrix%n), &
                matrix%A(1:matrix%n), permutations,          &
                inv_permutations, lwork, iwork, 1, error_flag)

      if (error_flag /= 0) then
         write(*,*)'Sparse matrix optimal ordering algorithm FAILED.'
         write(*,*)'Exited with error code = ',error_flag
         write(*,*)'Sparse matrix size     = ',matrix%nr
         stop
      endif

      deallocate(iwork)

      end function optimal_ordering

!          ************************************************************
!          ** overlapping_elements:                                  **
!          ** returns the number of overlapping elements in two      **
!          ** sparse matrix structures                               **
!          ************************************************************

           function overlaps(A,B) result(n_overlaps)
           implicit none

           type(sparse), intent(in) :: A, B
           integer                  :: i,ia,ib,ja,jb, n_overlaps

!          Initialise null number of overlaps
           n_overlaps = 0

!          Check that the two matrices have the same dimensions
           if ( A%nr /= B%nr .or. A%nc /= B%nc) then
             write(*,*)'Cannot compare matrices with different sizes'
             write(*,*)'in overlaps, size(A) = ',A%nr,'x',A%nc
             write(*,*)'             size(B) = ',B%nr,'x',B%nc
             stop
           endif

!          Loop to count
           rows: do i = 1, A%nr
              Acols: do ia = A%IA(i), A%IA(i+1) - 1
                   ja = A%JA(ia)
                 Bcols: do ib = B%IA(i), B%IA(i+1) - 1
                   jb = B%JA(ib)
                   if (jb == ja) then
                      n_overlaps = n_overlaps + 1
                      cycle Acols
                   endif
                 end do Bcols
              end do Acols
           end do rows

           end function overlaps

           !   *********************************************************
           !   **  Print sparse matrix contents to file               **
           !   **                                                     **
           !   **   Author:      Federico Perini                      **
           !   **   Last update: wednesday, 20/05/2012                **
           !   *********************************************************

           subroutine print_sparse_to_file(spmat,filename)
           implicit none

           type(sparse),     intent(in)                   :: spmat
           character(len=*), intent(in)                   :: filename
           real (dp)       , dimension(:,:), allocatable  :: dense
           integer                                        :: i, j
           character(len=*), parameter                    :: &
     &       fmt_er  = "(' Unable to print unallocated matrix ',&
     &                    'in print_sparse_to_file')",&
     &       fmt_out = "(1000(1x,E26.19))"

           ! Preliminary check
           if (.not.(spmat%nr/=0 .and. spmat%nc/=0)) then
              write(*,fmt_er)
              stop
           endif

           ! Allocate dense matrix
           allocate(dense(spmat%nr,spmat%nc))
           dense = spmat

           open(unit = 1000, file = filename)

              rows: do i = 1, spmat%nr
                 write(1000,fmt_out)(dense(i,j),j = 1,spmat%nc)
              end do rows

           close(1000)
           deallocate(dense)

           end subroutine print_sparse_to_file
!          ************************************************************

           subroutine print_sparse_ordered_to_file(spmat,filename)
           implicit none

           type(sparse_ordered),     intent(in)           :: spmat
           character(len=*), intent(in)                   :: filename
           real (dp)       , dimension(:,:), allocatable  :: dense
           integer                                        :: i, j
           character(len=*), parameter                    :: &
     &       fmt_er  = "(' Unable to print unallocated matrix ',&
     &                    'in print_sparse_to_file')",&
     &       fmt_out = "(1000(1x,E26.19))"

           ! Preliminary check
           if (.not.(spmat%nr/=0 .and. spmat%nc/=0)) then
              write(*,fmt_er)
              stop
           endif

           ! Allocate dense matrix
           allocate(dense(spmat%nr,spmat%nc))
           dense = spmat

           open(unit = 1000, file = filename)

              rows: do i = 1, spmat%nr
                 write(1000,fmt_out)(dense(i,j),j = 1,spmat%nc)
              end do rows

           close(1000)
           deallocate(dense)

           end subroutine print_sparse_ordered_to_file
!          ************************************************************


           !   *********************************************************
           !   **  Check if sparse matrix is correctly allocated      **
           !   *********************************************************

           function sparse_allocated(spmatrix) result(isallocated)
           implicit none

           type(sparse), intent(in) :: spmatrix
           logical                  :: isallocated

             ! Initialise as not correctly allocated
             isallocated = .false.

             ! Check that all the allocation parameters are satisfied
             if ( associated(spmatrix%val ) .and. &
                  associated(spmatrix%col ) .and. &
                  associated(spmatrix%ptrB) .and. &
                  associated(spmatrix%ptrE) .and. &
                  allocated (spmatrix%A   ) .and. &
                  allocated (spmatrix%JA  ) .and. &
                  allocated (spmatrix%IA  ) .and. &
                  spmatrix%n  <= size(spmatrix%A ) .and. &
                  spmatrix%n  <= size(spmatrix%JA) .and. &
                  spmatrix%nr == size(spmatrix%IA) - 1 .and. &
                  spmatrix%nc >  0  .and. &
                  spmatrix%nr >= 0 ) &
             isallocated = .true.


           end function sparse_allocated

           function sparseint_allocated(spmatrix) result(isallocated)
           implicit none

           type(sparseint), intent(in) :: spmatrix
           logical                     :: isallocated

             ! Initialise as not correctly allocated
             isallocated = .false.

             ! Check that all the allocation parameters are satisfied
             if ( associated(spmatrix%val ) .and. &
                  associated(spmatrix%col ) .and. &
                  associated(spmatrix%ptrB) .and. &
                  associated(spmatrix%ptrE) .and. &
                  allocated (spmatrix%A   ) .and. &
                  allocated (spmatrix%JA  ) .and. &
                  allocated (spmatrix%IA  ) .and. &
                  spmatrix%n  <= size(spmatrix%A ) .and. &
                  spmatrix%n  <= size(spmatrix%JA) .and. &
                  spmatrix%nr == size(spmatrix%IA) - 1 .and. &
                  spmatrix%nc >  0  .and. &
                  spmatrix%nr >  0 ) &
             isallocated = .true.


           end function sparseint_allocated


           !   *********************************************************
           !   **  Check if sparse matrix is correctly allocated      **
           !   *********************************************************

           function sparse_ordered_allocated(spmatrix) result(isallocated)
           implicit none

           type(sparse_ordered), intent(in) :: spmatrix
           logical                          :: isallocated

             ! Initialise as not correctly allocated
             isallocated = .false.

             ! Check that all the allocation parameters are satisfied
             if ( .not.associated(spmatrix%val ) ) return
             if ( .not.associated(spmatrix%col ) )return
             if ( .not.associated(spmatrix%ptrB) )return
             if ( .not.associated(spmatrix%ptrE) )return
             if ( .not.allocated (spmatrix%A   ) )return
             if ( .not.allocated (spmatrix%JA  ) ) return
             if ( .not.allocated (spmatrix%IA  ) ) return
             if ( spmatrix%n  > size(spmatrix%A ) ) return
             if ( spmatrix%n  > size(spmatrix%JA) ) return
             if ( spmatrix%nr /= size(spmatrix%IA) - 1 )return
             if ( spmatrix%nc ==  0  ) return
             if ( spmatrix%nr ==  0  ) return

             isallocated = .true.


           end function sparse_ordered_allocated


!     ******************************************************************
!     ** PRINT MATRIX PROPERTIES                                      **
!     ** Print to screen or other fortran unit the detailed matrix    **
!     ** properties                                                   **
!     ******************************************************************
      subroutine print_sparse_details(A, label, unit)
      implicit none

      type(sparse),      intent(in)           :: A
      character(len=*),  intent(in), optional :: label
      integer,           intent(in), optional :: unit

      character(len=20)  :: A_name
      integer            :: unit_number

      character(len=*), parameter ::                                   &
        fmt_head  = "(1x,71('-'))",                                    &
        fmt_head2 = "(2x,'Properties of sparse matrix object: ',A20)", &
        fmt_head3 = "(1x,'Designated integer dimension values: ')",    &
        fmt_nr    = "(4x,'- matrix%nr = ',I5)",                        &
        fmt_nc    = "(4x,'- matrix%nc = ',I5)",                        &
        fmt_n     = "(4x,'- matrix%n  = ',I5)",                        &
        fmt_head4 = "(1x,'Allocatable array components:  ')",          &
        fmt_A     = "(4x,'- matrix%A  = ',A3,' allocated, size=',I5)", &
        fmt_IA    = "(4x,'- matrix%IA = ',A3,' allocated, size=',I5)", &
        fmt_JA    = "(4x,'- matrix%JA = ',A3,' allocated, size=',I5)", &
        fmt_head5 = "(1x,'Pointer definitions: ')",                    &
        fmt_val   = "(4x,'- matrix%val = ',A3,' associated, size',I5)",&
        fmt_col   = "(4x,'- matrix%col = ',A3,' associated, size',I5)",&
        fmt_ptrB  = "(4x,'- matrix%ptrB= ',A3,' associated, size',I5)",&
        fmt_ptrE  = "(4x,'- matrix%ptrE= ',A3,' associated, size',I5)",&
        fmt_head6 = "(1x,'Check for null components:')",               &
        fmt_nA    = "(4x,'- in matrix%A : ',I5,' zeros found ')",      &
        fmt_nIA   = "(4x,'- in matrix%IA: ',I5,' zeros found ')",      &
        fmt_nJA   = "(4x,'- in matrix%JA: ',I5,' zeros found ')"


      if (     present(label)) A_name = label
      if (.not.present(label)) A_name = adjustl('(no_name)')

      if (present(unit)) then
         unit_number = unit
      else
         unit_number = 6 ! print to screen
      endif

      ! Print header
      write(unit_number,fmt_head )
      write(unit_number,fmt_head2)A_name
      write(unit_number,fmt_head )
      write(unit_number,*        )

      ! Print static matrix integers
      write(unit_number,fmt_head3)
      write(unit_number,fmt_nr   )A%nr
      write(unit_number,fmt_nc   )A%nc
      write(unit_number,fmt_n    )A%n
      write(unit_number,*        )

      ! Print allocatable components details
      write(unit_number,fmt_head4)
      write(unit_number,fmt_A    )logprt(allocated(A%A)),size(A%A)
      write(unit_number,fmt_IA   )logprt(allocated(A%IA)),size(A%IA)
      write(unit_number,fmt_JA   )logprt(allocated(A%JA)),size(A%JA)
      write(unit_number,*        )

      ! Print pointer components details
      write(unit_number,fmt_head5)
      write(unit_number,fmt_val  )logprt(associated(A%val )),size(A%val )
      write(unit_number,fmt_col  )logprt(associated(A%col )),size(A%col )
      write(unit_number,fmt_ptrB )logprt(associated(A%ptrB)),size(A%ptrB)
      write(unit_number,fmt_ptrE )logprt(associated(A%ptrE)),size(A%ptrE)
      write(unit_number,*        )

      ! Check for zeroes
      write(unit_number,fmt_head6)
      write(unit_number,fmt_nA   )count(A%A(1:A%n)==0.e0_dp)
      write(unit_number,fmt_nIA  )count(A%IA(1:A%nr+1)==0)
      write(unit_number,fmt_nJA  )count(A%JA(1:A%n)==0)
      write(unit_number,*        )

      write(unit_number,fmt_head )




      end subroutine print_sparse_details

!     ******************************************************************
!     ** PRINT MATRIX PROPERTIES                                      **
!     ** Print to screen or other fortran unit the detailed matrix    **
!     ** properties                                                   **
!     ******************************************************************
      subroutine print_sparseint_details(A, label, unit)
      implicit none

      type(sparseint),   intent(in)           :: A
      character(len=*),  intent(in), optional :: label
      integer,           intent(in), optional :: unit

      character(len=20)  :: A_name
      integer            :: unit_number

      character(len=*), parameter ::                                   &
        fmt_head  = "(1x,71('-'))",                                    &
        fmt_head2 = "(2x,'Properties of sparse matrix object: ',A20)", &
        fmt_head3 = "(1x,'Designated integer dimension values: ')",    &
        fmt_nr    = "(4x,'- matrix%nr = ',I5)",                        &
        fmt_nc    = "(4x,'- matrix%nc = ',I5)",                        &
        fmt_n     = "(4x,'- matrix%n  = ',I5)",                        &
        fmt_head4 = "(1x,'Allocatable array components:  ')",          &
        fmt_A     = "(4x,'- matrix%A  = ',A3,' allocated, size=',I5)", &
        fmt_IA    = "(4x,'- matrix%IA = ',A3,' allocated, size=',I5)", &
        fmt_JA    = "(4x,'- matrix%JA = ',A3,' allocated, size=',I5)", &
        fmt_head5 = "(1x,'Pointer definitions: ')",                    &
        fmt_val   = "(4x,'- matrix%val = ',A3,' associated, size',I5)",&
        fmt_col   = "(4x,'- matrix%col = ',A3,' associated, size',I5)",&
        fmt_ptrB  = "(4x,'- matrix%ptrB= ',A3,' associated, size',I5)",&
        fmt_ptrE  = "(4x,'- matrix%ptrE= ',A3,' associated, size',I5)",&
        fmt_head6 = "(1x,'Check for null components:')",               &
        fmt_nA    = "(4x,'- in matrix%A : ',I5,' zeros found ')",      &
        fmt_nIA   = "(4x,'- in matrix%IA: ',I5,' zeros found ')",      &
        fmt_nJA   = "(4x,'- in matrix%JA: ',I5,' zeros found ')"


      if (     present(label)) A_name = label
      if (.not.present(label)) A_name = adjustl('(no_name)')

      if (present(unit)) then
         unit_number = unit
      else
         unit_number = 6 ! print to screen
      endif

      ! Print header
      write(unit_number,fmt_head )
      write(unit_number,fmt_head2)A_name
      write(unit_number,fmt_head )
      write(unit_number,*        )

      ! Print static matrix integers
      write(unit_number,fmt_head3)
      write(unit_number,fmt_nr   )A%nr
      write(unit_number,fmt_nc   )A%nc
      write(unit_number,fmt_n    )A%n
      write(unit_number,*        )

      ! Print allocatable components details
      write(unit_number,fmt_head4)
      write(unit_number,fmt_A    )logprt(allocated(A%A)),size(A%A)
      write(unit_number,fmt_IA   )logprt(allocated(A%IA)),size(A%IA)
      write(unit_number,fmt_JA   )logprt(allocated(A%JA)),size(A%JA)
      write(unit_number,*        )

      ! Print pointer components details
      write(unit_number,fmt_head5)
      write(unit_number,fmt_val  )logprt(associated(A%val )),size(A%val )
      write(unit_number,fmt_col  )logprt(associated(A%col )),size(A%col )
      write(unit_number,fmt_ptrB )logprt(associated(A%ptrB)),size(A%ptrB)
      write(unit_number,fmt_ptrE )logprt(associated(A%ptrE)),size(A%ptrE)
      write(unit_number,*        )

      ! Check for zeroes
      write(unit_number,fmt_head6)
      write(unit_number,fmt_nA   )count(A%A(1:A%n)==0.e0_dp)
      write(unit_number,fmt_nIA  )count(A%IA(1:A%nr+1)==0)
      write(unit_number,fmt_nJA  )count(A%JA(1:A%n)==0)
      write(unit_number,*        )

      write(unit_number,fmt_head )




      end subroutine print_sparseint_details

      ! Print details of a sparse_ordered matrix
      subroutine print_sparse_ordered_details(A, label, unit)
      implicit none

      type(sparse_ordered),  intent(in)           :: A
      character(len=*),      intent(in), optional :: label
      integer,               intent(in), optional :: unit

      character(len=20)  :: A_name
      integer            :: unit_number

      character(len=*), parameter ::                                   &
        fmt_head  = "(1x,71('-'))",                                    &
        fmt_head2 = "(2x,'Properties of sparse matrix object: ',A20)", &
        fmt_head3 = "(1x,'Designated integer dimension values: ')",    &
        fmt_nr    = "(4x,'- matrix%nr = ',I5)",                        &
        fmt_nc    = "(4x,'- matrix%nc = ',I5)",                        &
        fmt_n     = "(4x,'- matrix%n  = ',I5)",                        &
        fmt_head4 = "(1x,'Allocatable array components:  ')",          &
        fmt_A     = "(4x,'- matrix%A  = ',A3,' allocated, size=',I5)", &
        fmt_IA    = "(4x,'- matrix%IA = ',A3,' allocated, size=',I5)", &
        fmt_JA    = "(4x,'- matrix%JA = ',A3,' allocated, size=',I5)", &
        fmt_head5 = "(1x,'Pointer definitions: ')",                    &
        fmt_val   = "(4x,'- matrix%val = ',A3,' associated, size',I5)",&
        fmt_col   = "(4x,'- matrix%col = ',A3,' associated, size',I5)",&
        fmt_ptrB  = "(4x,'- matrix%ptrB= ',A3,' associated, size',I5)",&
        fmt_ptrE  = "(4x,'- matrix%ptrE= ',A3,' associated, size',I5)",&
        fmt_head6 = "(1x,'Linear systems related arrays: ')",          &
        fmt_ord   = "(4x,'- rows/columns orderind  ',A3,' complete')", &
        fmt_sym   = "(4x,'- symbolic factorization ',A3,' complete')", &
        fmt_num   = "(4x,'- numerical factorization',A3,' complete')", &
        fmt_head7 = "(1x,'Check for null components:')",               &
        fmt_nA    = "(4x,'- in matrix%A : ',I5,' zeros found ')",      &
        fmt_nIA   = "(4x,'- in matrix%IA: ',I5,' zeros found ')",      &
        fmt_nJA   = "(4x,'- in matrix%JA: ',I5,' zeros found ')"



      !!!!!!! TO DO: ORDERING AND LINEAR-SYSTEM RELATED PART!

      if (     present(label)) A_name = label
      if (.not.present(label)) A_name = adjustl('(no_name)')

      if (present(unit)) then
         unit_number = unit
      else
         unit_number = 6 ! print to screen
      endif

      ! Print header
      write(unit_number,fmt_head )
      write(unit_number,fmt_head2)A_name
      write(unit_number,fmt_head )
      write(unit_number,*        )

      ! Print static matrix integers
      write(unit_number,fmt_head3)
      write(unit_number,fmt_nr   )A%nr
      write(unit_number,fmt_nc   )A%nc
      write(unit_number,fmt_n    )A%n
      write(unit_number,*        )

      ! Print allocatable components details
      write(unit_number,fmt_head4)
      write(unit_number,fmt_A    )logprt(allocated(A%A)),size(A%A)
      write(unit_number,fmt_IA   )logprt(allocated(A%IA)),size(A%IA)
      write(unit_number,fmt_JA   )logprt(allocated(A%JA)),size(A%JA)
      write(unit_number,*        )

      ! Print pointer components details
      write(unit_number,fmt_head5)
      write(unit_number,fmt_val  )logprt(associated(A%val )),size(A%val )
      write(unit_number,fmt_col  )logprt(associated(A%col )),size(A%col )
      write(unit_number,fmt_ptrB )logprt(associated(A%ptrB)),size(A%ptrB)
      write(unit_number,fmt_ptrE )logprt(associated(A%ptrE)),size(A%ptrE)
      write(unit_number,*        )

      ! Linear systems-related details
      write(unit_number,fmt_head6)
      write(unit_number,fmt_ord  )logprt(A%is_ordered)
      write(unit_number,fmt_sym  )logprt(A%symbolically_factorized)
      write(unit_number,fmt_num  )logprt(A%numerically_factorized)
      write(unit_number,*        )

      ! Check for zeroes
      write(unit_number,fmt_head7)
      write(unit_number,fmt_nA   )count(A%A(1:A%n)==0.e0_dp)
      write(unit_number,fmt_nIA  )count(A%IA(1:A%nr+1)==0)
      write(unit_number,fmt_nJA  )count(A%JA(1:A%n)==0)
      write(unit_number,*        )


      write(unit_number,fmt_head )




      end subroutine print_sparse_ordered_details


      ! ****************************************************************
      ! Function that gives a three-character flag for true/false items
      ! ****************************************************************
      function logprt(log_value) result(label)
      implicit none

      logical, intent(in) :: log_value
      character(len=3)    :: label

      if (     log_value) label = ' IS'
      if (.not.log_value) label = 'NOT'

      end function logprt

      end module sparse_definitions
