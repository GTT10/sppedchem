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
!     **                  Functions and utilities                    **
!     **                                                             **
!     **   Author:      Federico Perini                              **
!     **   Date created: tuesday, 25/05/2010                         **
!     **   Last update : tuesday, 31/07/2012                         **
!     **                                                             **
!     *****************************************************************

module utilities

   use working_precision

   implicit none
   private

   ! ****************************************************************
   ! ** ROUTINES CALLABLE FROM PRESENT MODULE                      **
   ! ****************************************************************

   ! Compute relative error between two values/array/matrices
   ! function: relative_error(A,B)
   public :: relative_error

   ! Sort arrays in ascending order using quicksort algorithm
   ! subroutine: call sort(A)
   public :: sort

   ! Allocate/Rellocate an array with new dimensions
   ! subroutine: call force_allocate(A, dim1)
   !             call force_allocate(A, dim1, dim2)
   public :: force_allocate
   public :: extend

   ! Extend/shrink an array's dimensions keeping the contained
   ! data possible in the new shape
   ! subroutine: call change_size(A, dim1)
   public :: change_size

   ! Detect NaN or Inf values in arrays
!      public :: isnan

   ! Convert integer into string of given length, with
   ! leading zeros
   public :: int_to_string
   public :: find_string_in_array, find_stringi_in_array

   ! ****************************************************************
   ! ** RELATIVE ERROR                                             **
   ! ****************************************************************
   interface relative_error
      module procedure dble_relative_error
      module procedure int_relative_error
   end interface relative_error

   ! ****************************************************************
   ! ** QUICK SORT                                                 **
   ! ****************************************************************
   interface sort
      module procedure double_quick_sort
      module procedure int_quick_sort
   end interface sort

   ! ****************************************************************
   ! ** ALLOCATION AND REALLOCATION OF ARRAYS                      **
   ! ****************************************************************
   interface force_allocate
      module procedure reallocate_double_array
      module procedure reallocate_double_matrix
      module procedure reallocate_int_array
      module procedure reallocate_int_matrix
   end interface force_allocate

   interface force_deallocate
      module procedure deallocate_double_array
      module procedure deallocate_double_matrix
      module procedure deallocate_int_array
      module procedure deallocate_int_matrix
   end interface force_deallocate

   interface extend
      module procedure extend_double_array
      module procedure extend_char_array
      module procedure extend_double_matrix
      module procedure extend_int_array
   end interface extend

   interface change_size
      module procedure change_size_double_array
      module procedure change_size_int_array
   end interface change_size


   ! ****************************************************************
   ! ****************************************************************
contains


   ! ****************************************************************
   ! ** Relative and absolute error                                **
   ! ****************************************************************
   elemental function dble_relative_error(a,b) result(epsi)
      implicit none

      real (dp)       , intent(in) :: a, b
      real (dp)                    :: epsi

      denominator_case: if (b == 0.e0_dp) then

         if (a==0.e0_dp) then
            epsi = 0.e0_dp
         else
            epsi = huge(0.e0_dp)
         endif

      else
         epsi = abs(a-b)/b
      endif denominator_case

   end function dble_relative_error

   elemental function int_relative_error(a,b) result(epsi)
      implicit none

      integer,          intent(in) :: a, b
      real (dp)                    :: epsi

      denominator_case: if (b == 0) then

         if (a==0) then
            epsi = 0.e0_dp
         else
            epsi = huge(0.e0_dp)
         endif

      else
         epsi = abs(real(a-b, dp))/real(b, dp)
      endif denominator_case

   end function int_relative_error

!     ******************************************************************
!     ** ALLOCATION and REALLOCATION                                  **
!     ******************************************************************
   subroutine reallocate_double_array(A,dim)
      implicit none

      real (dp)       , dimension(:), allocatable, intent(inout) :: A
      integer,                                     intent(in)    :: dim

      if (dim < 0) then
         write(*,*)'Cannot allocate with negative size:',dim
         stop
      endif

      ! Check for current allocation
      allocation_state: if (allocated(A)) then

         ! Check for current array size
         change_needed: if (size(A) /= dim) then
            deallocate(A)
            allocate(A(dim))
         endif change_needed

      else

         ! Allocate array for sure
         allocate(A(dim))

      endif allocation_state

   end subroutine reallocate_double_array

   subroutine reallocate_int_array(A,dim)
      implicit none

      integer, dimension(:), allocatable, intent(inout) :: A
      integer,                            intent(in)    :: dim

      if (dim < 0) then
         write(*,*)'Cannot allocate with negative size:',dim
         stop
      endif

      ! Check for current allocation
      allocation_state: if (allocated(A)) then

         ! Check for current array size
         change_needed: if (size(A) /= dim) then
            deallocate(A)
            allocate(A(dim))
         endif change_needed

      else

         ! Allocate array for sure
         allocate(A(dim))

      endif allocation_state

   end subroutine reallocate_int_array

   subroutine reallocate_double_matrix(A,dim1,dim2)
      implicit none

      real (dp)       , dimension(:,:), allocatable, intent(inout) :: A
      integer,                              intent(in)    :: dim1,dim2

      if (dim1 < 0 .or. dim2<0) then
         write(*,*)'Cannot allocate with negative size:',dim1,dim2
         stop
      endif

      ! Check for current allocation
      allocation_state: if (allocated(A)) then

         ! Check for current array size
         change_needed: if (size(A,1) /= dim1 .or. size(A,2)/=dim2) then
            deallocate(A)
            allocate(A(dim1,dim2))
         endif change_needed

      else

         ! Allocate array for sure
         allocate(A(dim1,dim2))

      endif allocation_state

   end subroutine reallocate_double_matrix

   subroutine reallocate_int_matrix(A,dim1,dim2)
      implicit none

      integer, dimension(:,:), allocatable, intent(inout) :: A
      integer,                              intent(in)    :: dim1,dim2

      if (dim1 < 0 .or. dim2<0) then
         write(*,*)'Cannot allocate with negative size:',dim1,dim2
         stop
      endif

      ! Check for current allocation
      allocation_state: if (allocated(A)) then

         ! Check for current array size
         change_needed: if (size(A,1) /= dim1 .or. size(A,2)/=dim2) then
            deallocate(A)
            allocate(A(dim1,dim2))
         endif change_needed

      else

         ! Allocate array for sure
         allocate(A(dim1,dim2))

      endif allocation_state

   end subroutine reallocate_int_matrix


!     ******************************************************************
!     ** DEALLOCATE ARRAYS                                            **
!     ******************************************************************
   subroutine deallocate_double_array(A)
      implicit none

      real (dp)       , dimension(:), allocatable, intent(inout) :: A

      ! Complete deallocation
      if (allocated(A)) deallocate(A)

   end subroutine deallocate_double_array

   subroutine deallocate_double_matrix(A)
      implicit none

      real (dp)       , dimension(:,:), allocatable, intent(inout) :: A

      ! Complete deallocation
      if (allocated(A)) deallocate(A)

   end subroutine deallocate_double_matrix

   subroutine deallocate_int_array(A)
      implicit none

      integer, dimension(:), allocatable, intent(inout) :: A

      ! Complete deallocation
      if (allocated(A)) deallocate(A)

   end subroutine deallocate_int_array

   subroutine deallocate_int_matrix(A)
      implicit none

      integer, dimension(:,:), allocatable, intent(inout) :: A

      ! Complete deallocation
      if (allocated(A)) deallocate(A)

   end subroutine deallocate_int_matrix

!     ******************************************************************
!     ** CHANGE SIZE OF ARRAYS                                        **
!     ******************************************************************

   subroutine change_size_double_array(A,dim)
      implicit none

      real (dp)       , dimension(:), allocatable, intent(inout) :: A
      real (dp)       , dimension(:), allocatable                :: tmp
      integer,                                     intent(in)    :: dim


      if (dim <= 0) then
         write(*,*)'Cannot change size with negative size:',dim
         stop
      endif

      ! Check for current allocation
      allocation_state: if (allocated(A)) then

         ! Check for current array size
         change_needed: if (size(A) /= dim) then

            allocate(tmp(size(A)))
            tmp = A

            deallocate(A)
            allocate(A(dim))

            A(1:min(dim,size(tmp))) = tmp(1:min(dim,size(tmp)))
            A(min(dim,size(tmp))+1:dim) = 0.e0_dp

            deallocate(tmp)
         endif change_needed

      else

         ! Allocate array for sure
         allocate(A(dim))

      endif allocation_state

   end subroutine change_size_double_array

   subroutine change_size_int_array(A,dim)
      implicit none

      integer, dimension(:), allocatable, intent(inout) :: A
      integer, dimension(:), allocatable                :: tmp
      integer,                            intent(in)    :: dim


      if (dim <= 0) then
         write(*,*)'Cannot change size with negative size:',dim
         stop
      endif

      ! Check for current allocation
      allocation_state: if (allocated(A)) then

         ! Check for current array size
         change_needed: if (size(A) /= dim) then

            allocate(tmp(size(A)))
            tmp = A

            deallocate(A)
            allocate(A(dim))

            A(1:min(dim,size(tmp))) = tmp(1:min(dim,size(tmp)))
            A(min(dim,size(tmp))+1:dim) = 0

            deallocate(tmp)
         endif change_needed

      else

         ! Allocate array for sure
         allocate(A(dim))

      endif allocation_state

   end subroutine change_size_int_array




!     ******************************************************************
!     ** QUICKSORT algorithm                                          **
!     ******************************************************************

   recursive subroutine double_quick_sort(list)

      real (dp)       , dimension(:), intent(inout) :: list

      integer            :: i, j, n
      integer, parameter :: max_simple_sort_size = 6
      real (dp)          :: chosen, temp

      n = size(list)

      choose_sorting_algorithm: if (n <= max_simple_sort_size) then
         ! Use interchange sort for small lists
         call double_interchange_sort(list)
      else
         ! Use partition (“quick”) sort if the list is big
         chosen = list(int(n/2))
         i = 0
         j = n + 1

         scan_lists: do
            ! Scan list from left end
            ! until element >= chosen is found
            scan_from_left: do
               i = i + 1
               if (list(i) >= chosen) exit scan_from_left
            end do scan_from_left

            ! Scan list from right end
            ! until element <= chosen is found

            scan_from_right: do
               j = j - 1
               if (list(j) <= chosen) exit scan_from_right
            end do scan_from_right

            swap: if (i < j) then
               ! Swap two out of place elements
               temp = list(i)
               list(i) = list(j)
               list(j) = temp
            else if (i == j) then
               i = i + 1
               exit
            else
               exit
            endif swap

         end do scan_lists

         if (1 < j) call double_quick_sort(list(:j))
         if (i < n) call double_quick_sort(list(i:))

      end if choose_sorting_algorithm ! test for small array

   end subroutine double_quick_sort

   recursive subroutine int_quick_sort(list)

      integer, dimension(:), intent(inout) :: list

      integer            :: i, j, n, chosen, temp
      integer, parameter :: max_simple_sort_size = 6

      n = size(list)

      choose_sorting_algorithm: if (n <= max_simple_sort_size) then
         ! Use interchange sort for small lists
         call int_interchange_sort(list)
      else
         ! Use partition (“quick”) sort if the list is big
         chosen = list(int(n/2))
         i = 0
         j = n + 1

         scan_lists: do
            ! Scan list from left end
            ! until element >= chosen is found
            scan_from_left: do
               i = i + 1
               if (list(i) >= chosen) exit scan_from_left
            end do scan_from_left

            ! Scan list from right end
            ! until element <= chosen is found

            scan_from_right: do
               j = j - 1
               if (list(j) <= chosen) exit scan_from_right
            end do scan_from_right

            swap: if (i < j) then
               ! Swap two out of place elements
               temp = list(i)
               list(i) = list(j)
               list(j) = temp
            else if (i == j) then
               i = i + 1
               exit
            else
               exit
            endif swap

         end do scan_lists

         if (1 < j) call int_quick_sort(list(:j))
         if (i < n) call int_quick_sort(list(i:))

      end if choose_sorting_algorithm ! test for small array

   end subroutine int_quick_sort


   subroutine double_interchange_sort(list)

      real (dp)       , dimension(:), intent(inout) :: list

      integer :: i, j
      real (dp)        :: temp

      do i = 1, size(list) - 1
         do j = i + 1, size(list)
            if (list(i) > list(j)) then
               temp = list(i)
               list(i) = list(j)
               list(j) = temp
            end if
         end do
      end do

   end subroutine double_interchange_sort

   subroutine int_interchange_sort(list)

      integer, dimension(:), intent(inout) :: list

      integer :: i, j, temp

      do i = 1, size(list) - 1
         do j = i + 1, size(list)
            if (list(i) > list(j)) then
               temp = list(i)
               list(i) = list(j)
               list(j) = temp
            end if
         end do
      end do

   end subroutine int_interchange_sort


   subroutine extend_double_array(A,dim)
      implicit none

      real (dp)       , dimension(:), allocatable, intent(inout) :: A
      integer,                                     intent(in)    :: dim
      integer :: state
      character(len=*), parameter ::&
      &fmt_er = "(' Error: allocation failed. Size= ',I10,' state=',I9)"

      if (dim < 0) then
         write(*,*)'Cannot allocate with negative size:',dim
         stop
      endif
      state = 0
      ! Check for current allocation
      allocation_state: if (allocated(A)) then

         ! Check for current array size
         change_needed: if (size(A) < dim) then
            deallocate(A)
            allocate(A(dim), stat=state)
         endif change_needed

      else

         ! Allocate array for sure
         allocate(A(dim), stat=state)

      endif allocation_state

      if (state /= 0) then
         write(*,fmt_er)dim,state
         stop
      endif

   end subroutine extend_double_array

   subroutine extend_int_array(A,dim)
      implicit none

      integer, dimension(:), allocatable, intent(inout) :: A
      integer,                            intent(in)    :: dim
      integer :: state
      character(len=*), parameter ::&
      &fmt_er = "(' Error: allocation failed. Size= ',I10,' state=',I9)"

      if (dim < 0) then
         write(*,*)'Cannot allocate with negative size:',dim
         stop
      endif
      state = 0
      ! Check for current allocation
      allocation_state: if (allocated(A)) then

         ! Check for current array size
         change_needed: if (size(A) < dim) then
            deallocate(A)
            allocate(A(dim), stat = state)
         endif change_needed

      else

         ! Allocate array for sure
         allocate(A(dim), stat = state)

      endif allocation_state

      if (state /= 0) then
         write(*,fmt_er)dim,state
         stop
      endif

   end subroutine extend_int_array

   subroutine extend_char_array(A,dim)
      implicit none

      character(len=*), dimension(:), allocatable, intent(inout) :: A
      integer,                                     intent(in)    :: dim
      integer :: lgth
      integer :: state
      character(len=*), parameter ::&
      &fmt_er = "(' Error: allocation failed. Size= ',I10,' state=',I9)"

      if (dim < 0) then
         write(*,*)'Cannot allocate with negative size:',dim
         stop
      endif

      state = 0
      ! Check for current allocation
      allocation_state: if (allocated(A)) then

         ! Check for current array size
         change_needed: if (size(A) < dim) then
            deallocate(A)
            allocate(A(dim), stat=state)
         endif change_needed

      else

         ! Allocate array for sure
         allocate(A(dim), stat=state)

      endif allocation_state

      if (state /= 0) then
         write(*,fmt_er)dim,state
         stop
      endif

   end subroutine extend_char_array



   subroutine extend_double_matrix(A,dim1, dim2)
      implicit none

      real (dp)       , dimension(:,:), allocatable, intent(inout) :: A
      integer,                                  intent(in)    :: dim1
      integer,                                  intent(in)    :: dim2

      if (dim1 < 0 .or. dim2<0) then
         write(*,*)'Cannot allocate with negative size:',dim1,dim2
         stop
      endif

      ! Check for current allocation
      allocation_state: if (allocated(A)) then

         ! Check for current array size
         change_needed: if (size(A,1) < dim1 .or. size(A,2)<dim2) then
            deallocate(A)
            allocate(A(dim1,dim2))
         endif change_needed

      else

         ! Allocate array for sure
         allocate(A(dim1,dim2))

      endif allocation_state

   end subroutine extend_double_matrix


   ! ****************************************************************
   ! ** Convert integer into string with leading zeros             **
   ! ****************************************************************
   function int_to_string(num,length) result(str)
      implicit none

      integer, intent(in)   :: num, length
      character(len=length) :: str

      integer               :: i, numtmp, numprt, ipos

      ! Prepare heading zeros
      do i = 1, length
         str (i:i) = '0'
      end do

      if (num < 0) return

      ! Initialise number
      numtmp = num

      ipos = length
      do while (numtmp > 0)

         ! Check for number too long
         if (ipos < 1) then
            do i = 1, length
               str(i:i) = '*'
            end do
            exit
         endif

         numprt = numtmp - int(numtmp/10)*10
         numtmp = numtmp / 10
         write(str(ipos:ipos),'(I1)')numprt
         ipos = ipos - 1

      end do

   end function int_to_string

   ! ****************************************************************
   ! ** Find an equal string in an array of strings                **
   ! ****************************************************************
   function find_string_in_array(string,strarray) result(ifind)
      implicit none

      character(len=*),               intent(in) :: string
      character(len=*), dimension(:), intent(in) :: strarray
      integer                                    :: ifind

      integer                                    :: i, n_list
      integer                                    :: lenl, lens

      ! Length of string and list
      lenl = len(strarray)
      lens = len(string)

      n_list = size(strarray,1)

      ! Initialise ifind as not found
      ifind = 0

      do i = 1, n_list

         if (trim(adjustl(string)) == trim(adjustl(strarray(i)))) then

            ifind = i
            return

         endif

      end do

   end function find_string_in_array


   ! ****************************************************************
   ! ** Find an equal string in an array of strings                **
   ! ** case insensitive version                                   **
   ! ****************************************************************
   function find_stringi_in_array(string,strarray) result(ifind)
      implicit none

      interface
         function fs_cap(s)
            character(len=*), intent(in) :: s
            character(len=len(s))        :: fs_cap
         end function fs_cap
      end interface


      character(len=*),               intent(in) :: string
      character(len=*), dimension(:), intent(in) :: strarray
      integer                                    :: ifind

      integer                                    :: i, n_list
      integer                                    :: lenl, lens

      ! Length of string and list
      lenl = len(strarray)
      lens = len(string)

      n_list = size(strarray,1)

      ! Initialise ifind as not found
      ifind = 0

      do i = 1, n_list

         if ( trim(adjustl(fs_cap(string)))   ==&
         &trim(adjustl(fs_cap(strarray(i))))) then

            ifind = i
            return

         endif

      end do

   end function find_stringi_in_array



end module utilities















!c     *****************************************************************
!c     **                                                             **
!c     **                     SpeedCHEM FORTRAN                       **
!c     **                                                             **
!c     **     Computing mass fraction conservation for a reacting     **
!c     **                        environment                          **
!c     **                                                             **
!c     **                                                             **
!c     **   Author:      Federico Perini                              **
!c     **   Last update: wedesday, 26/05/2010                         **
!c     **                                                             **
!c     *****************************************************************
!
!      subroutine massfr_conservation (domegadt,rho,dYdt)
!      use speedchem, only: ns,SCMW
!      implicit none
!
!      real (dp)       , dimension(ns), intent(in) :: domegadt
!      real (dp)       , intent(in) :: rho
!      real (dp)       , dimension(ns), intent(out) :: dYdt
!
!c     *****************************************************************
!
!      dYdt = 1000.d0 * domegadt * SCMW / rho
!
!      end subroutine






!c     *****************************************************************
!c     **                                                             **
!c     **                     SpeedCHEM FORTRAN                       **
!c     **                                                             **
!c     **   Computing reduced pressure value for pressure-dependent   **
!c     **                  reaction rate constants                    **
!c     **                                                             **
!c     **                                                             **
!c     **   Author:      Federico Perini                              **
!c     **   Last update: wedesday, 26/05/2010                         **
!c     **                                                             **
!c     *****************************************************************
!
!      subroutine redP (k0,kinf,M,Pr)
!      use speedchem, only: nr
!      implicit none
!
!      real (dp)       , dimension(nr), intent(in) :: k0,kinf,M
!      real (dp)       , dimension(nr), intent(out) :: Pr
!
!c     *****************************************************************
!
!      Pr = k0 * M / kinf
!
!      end subroutine







!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **      Converting dense matrix into sparse CSR format         **
!     **                                                             **
!     **                                                             **
!     **   Author:      Federico Perini                              **
!     **   Last update: thursday, 10/06/2010                         **
!     **                                                             **
!     *****************************************************************

!      subroutine CSR(nrows,ncols,matrice,vals,col_ind,row_ptrB,
!     &               row_ptrE,nvals,r)
!
!c      use find_mod
!
!      implicit none
!c     Subroutine input and output
!      integer, intent(in) :: nrows,ncols
!      real (dp)       , dimension(nrows,ncols), intent(in) :: matrice
!      real (dp)       , dimension(nrows*ncols), intent(out) :: vals
!
!      integer, dimension(nrows*ncols), intent(out) :: col_ind
!      integer, dimension(nrows), intent(out) :: row_ptrB, row_ptrE
!      integer, intent(out) :: nvals, r
!
!c     Other variables
!      integer :: i,j,k
!
!
!c     *****************************************************************
!
!      nvals = 0
!      r     = 0
!      k     = 0
!      do i=1,nrows
!        do j=1,ncols
!
!	    if (matrice(i,j).ne.0.d0) then
!	    nvals = nvals + 1
!	    vals(nvals)    = matrice(i,j)
!	    col_ind(nvals) = j
!
!	    if (i.ne.k) then
!	      r = r + 1
!	      row_ptrB(r) = nvals
!          if (r.gt.1) row_ptrE(r-1) = nvals
!	      k = i
!	    end if
!
!	  end if
!
!      end do
!      end do
!
!	  row_ptrE(r) = nvals+1

!      write(*,*)nvals,r
!     Counting number of matrix elements and initialising sparse
!      call find_indices2D(matrice.ne.0.d0)

!      write(*,*)nrows,ncols

!      nvals = size(i2D1)

!      vals(1:nvals)    = 0.d0
!      col_ind(1:nvals) = 0


!     Reordering i2D1 and i2D2 in row-wise format
!      do i=1,nvals
!        do j=i+1,nvals

!	 if (i2D1(j).le.i2D1(i)) then
!
!	  k = i2D1(j)
!	  i2D1(j) = i2D1(i)
!	  i2D1(i) = k
!
!	  k = i2D2(j)
!	  i2D2(j) = i2D2(i)
!	  i2D2(i) = k
!
!	 endif
!
!	end do
!
!	write(*,"(2I5)")i2D1(i),i2D2(i)
!     end do
!  c    pause

!     Counting number of rows used
!      k = 0
! c     r = 0

!      do i=1,nvals
!        if (i2D1(i).ne.k) then
!	k = i2D1(i)
!	r = r+1
!	end if
!      end do
!
!
!      k = 0
!      r = 0
!      do i = 1,nvals
!
!        vals(i)    = matrice(i2D1(i),i2D2(i))
!	col_ind(i) = i2D2(i)
!
!	if (i2D1(i).gt.k) then
!	  k = k+1
!	  row_ptr(k) = i2D1(i)
!	  r = row_ptr(k)
!	endif
!
!      end do


!      end subroutine CSR



elemental subroutine lower_case(word)
!     ! convert a word to lower case
   character (len=*) , intent(in out) :: word
   integer :: i,ic,nlen
   nlen = len(word)
   do i=1,nlen
      ic = ichar(word(i:i))
      if (ic >= 65 .and. ic < 90) word(i:i) = char(ic+32)
   end do
end subroutine lower_case










subroutine SC_ch_cap ( ch )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Discussion:
!
!    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions,
!    which guarantee the ASCII collating sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character CH, the character to capitalize.
!
   implicit none

   character              ch
   integer   ( kind = 4 ) itemp

   itemp = iachar ( ch )

   if ( 97 <= itemp .and. itemp <= 122 ) then
      ch = achar ( itemp - 32 )
   end if

   return

end


function fs_cap ( s )

!*****************************************************************************80
!
!! S_CAP replaces any lowercase letters by uppercase ones in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
   implicit none

   character              ch
   integer   ( kind = 4 ) i
   character ( len = * ), intent(in) ::  s
   character ( len = len(s) )        ::  fs_cap
   integer   ( kind = 4 ) s_length

   s_length = len_trim ( s )
   fs_cap    = s


   do i = 1, s_length

      ch = s(i:i)
      call SC_ch_cap ( ch )
      fs_cap(i:i) = ch

   end do

   return
end function fs_cap



subroutine s_cap ( s )

!*****************************************************************************80
!
!! S_CAP replaces any lowercase letters by uppercase ones in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
   implicit none

   character              ch
   integer   ( kind = 4 ) i
   character ( len = * )  s
   integer   ( kind = 4 ) s_length

   s_length = len_trim ( s )


   do i = 1, s_length

      ch = s(i:i)
      call SC_ch_cap ( ch )
      s(i:i) = ch

   end do

   return
end subroutine s_cap



