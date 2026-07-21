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
!     **         Module for sparse chemistry computations            **
!     **                                                             **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: wednesday, 23/11/2011                        **
!     **                                                             **
!     *****************************************************************
module sparse_algebra

   ! ** IMPORTANT: Working Precision *****************************
   ! ** The current working precision has to be set in the      **
   ! ** working_precision module, and is common to all the      **
   ! ** code. Note that MPI doesn't support quadruple precision **
   ! ** in Fortran.                                             **
   ! *************************************************************
   use working_precision

!        The type components and basic operations with sparse matrices
!        are defined in the sparse_definitions module
   use sparse_definitions

   implicit none
   public

   ! *************************************************************
   ! ** Interfaces for operations in common between             **
   ! ** sparse integer and sparse double matrices               **
   ! *************************************************************


   ! Deallocation operator
   interface spdeallocate
      module procedure sparse_nullify_general
      module procedure sparse_nullify_int
   end interface

   ! Allocation operator
   interface spallocate
      module procedure sparse_allocate_det
      module procedure sparseint_allocate_det
   end interface


   ! Greater-than (.gt., >) sign
   interface operator (>)
      module procedure sparse_gt_dble
      module procedure sparse_gt_array
   end interface !operator (>)







!        Arrays for sparse matrix computations
   integer :: lspwork,ls2work
   real (dp)       , dimension(:), allocatable :: spwork,s2work
!$ OMP    THREADPRIVATE(lspwork,ls2work,spwork,s2work)

contains


   !   *********************************************************
   !   **  Extract indices of the dense matrix from the       **
   !   **  sparse memory allocation, as a function of the     **
   !   **  array position in sparse form                      **
   !   **                                                     **
   !   **   Author:      Federico Perini                      **
   !   **   Last update: wednesday, 29/03/2012                **
   !   *********************************************************

   subroutine extract_dense_indices(spmatrix,Iidx,Jidx)
      implicit none

      type(sparse),                   intent(in)  :: spmatrix
      integer, dimension(spmatrix%n), intent(out) :: Iidx, Jidx
      integer :: i


!          Assign data
      do i = 1, spmatrix%nr

!            Column index
         Jidx(spmatrix%IA(i):spmatrix%IA(i+1)-1) =&
         &spmatrix%JA(spmatrix%IA(i):spmatrix%IA(i+1)-1)

!            Row index
         Iidx(spmatrix%IA(i):spmatrix%IA(i+1)-1) = i

      end do


   end subroutine extract_dense_indices
!          ************************************************************

   !   *********************************************************
   !   **  Extract position of element sparse_matrix%A in the **
   !   **  sparse matrix storage, corresponding to (i,j)      **
   !   **  element in dense form                              **
   !   **  NB: output is -1 if element isn't found            **
   !   **                                                     **
   !   **   Author:      Federico Perini                      **
   !   **   Last update: wednesday, 29/03/2012                **
   !   *********************************************************

   function extract_sparse_index(spmatrix,i,j) result(k)
      implicit none

      type(sparse),                   intent(in)  :: spmatrix
      integer,                        intent(in)  :: i, j
      integer                                     :: k

      integer :: ispcol

!          Preliminary checks
!           if (.not. (i > 0 .and. i <= spmatrix%nr)) then
!             write(*,*)' wrong row index in extract_sparse_index '
!             write(*,*)' i = ',i,'; spmatrix%nr = ',spmatrix%nr
!             stop
!           endif
!           if (.not. (j > 0 .and. j <= spmatrix%nc)) then
!             write(*,*)' wrong column index in extract_sparse_index '
!             write(*,*)' j = ',j,'; spmatrix%nc = ',spmatrix%nc
!             stop
!           endif

!          Initialize sparse index
      k = -1

!          Find sparse index
      if (i > 0 .and. i <= spmatrix%nr) then
         column_find: do ispcol = spmatrix%IA(i), spmatrix%IA(i+1)-1
            if (spmatrix%JA(ispcol) == j) then
               k = ispcol
               exit column_find
            endif
         end do column_find
      endif

!           if (k == -1) then
!              write(*,*)' Element (i,j) = (',i,',',j,') not found'
!              stop
!           endif


   end function extract_sparse_index
!          ************************************************************

   !   *********************************************************
   !   **  Adds one value to the sparse matrix representation **
   !   **  given the row and column indices in the full       **
   !   **  format                                             **
   !   **                                                     **
   !   **   Author:      Federico Perini                      **
   !   **   Last update: sunday, 22/07/2012                   **
   !   *********************************************************

   subroutine add_value(matrix,irow,jcol,val)
      use working_precision, only: dp
      use utilities, only: change_size, force_allocate
      implicit none

      type(sparse),     intent(inout), target :: matrix
      integer,          intent(in)            :: irow,jcol
      real (dp)       , intent(in)            :: val

      integer                         :: num_add_rows, j, japos,&
      &jaorig
      logical                         :: extend_rows,&
      &extend_columns,&
      &already_present
      integer,          dimension(:), allocatable :: JAtmp, IAtmp
      real (dp)       , dimension(:), allocatable :: Atmp




      ! Do nothing if value is zero
      if (val == 0.e0_dp) return


      ! If matrix was empty, it has to be initialised to contain
      ! just that value
      already_init: if ((.not.sparse_allocated(matrix))&
      &.or.matrix%n==0) then

         call allocate(max(irow, matrix%nr),&
         &max(jcol, matrix%nc), 1, matrix)

         matrix%A(1)  = val
         matrix%IA = [(1,j=1,irow),(2,j=irow+1,matrix%nr+1)]
         matrix%JA(1) = jcol

      else


         ! Value already exists?
         jaorig = extract_sparse_index(matrix,irow,jcol)
         already_present = jaorig /= -1


         no_mods: if (already_present) then

            japos = jaorig

            matrix%A(japos) = val

         else ! .not.already_present

            ! Preliminary checks on the total matrix dimensions
            extend_rows     = irow > matrix%nr
            extend_columns  = jcol > matrix%nc

            ! Prepare temporary arrays of column indices and values
!                 allocate(JAtmp(size(matrix%JA)+1))
!                 allocate(Atmp (size(matrix%A )+1))
            call change_size(matrix%JA,matrix%n+1)
            call change_size(matrix%A ,matrix%n+1)


            ! Row extension means that matrix%IA array has to be u
            ! pdated
            updated_IA: if (extend_rows) then

               ! Number of additional rows
               num_add_rows = irow - matrix%nr

               ! Reassign the IA array of indices

               allocate(IAtmp(matrix%nr+num_add_rows+1))

               ! Old matrix structure is preserved
               IAtmp(1:matrix%nr+1) = matrix%IA

               ! Add empty lines
               IAtmp(matrix%nr+2:matrix%nr+num_add_rows) =&
               &IAtmp(matrix%nr+1)

               ! Last line has just one entry (the currently added
               ! value)
               IAtmp(matrix%nr+num_add_rows+1) = IAtmp(matrix%nr+1)&
               &+ 1

               ! Update sparse matrix number of rows
               matrix%nr = irow

               ! Update matrix arrays using temporary storage
               matrix%JA(matrix%n+1) = jcol
               matrix%A(matrix%n+1)  = val

               call force_allocate(matrix%IA, matrix%nr+1)
               matrix%IA = IAtmp

               ! Deallocate temporary arrays
               deallocate(IAtmp)

            else

               ! If cannot find better values, put in the end of row
               japos = matrix%IA(irow+1)

               ! Find the correct position to insert the column
               ! index and the value
               do j = matrix%IA(irow+1)-1, matrix%IA(irow), -1
                  if (matrix%JA(j) > jcol) then
                     japos = j
                  else
                     exit
                  endif
               end do

               ! Update matrix arrays using temporary storage
               matrix%JA(1:matrix%n+1) = [ matrix%JA(1:japos-1),&
               &jcol,&
               &matrix%JA(japos:matrix%n)]

               matrix%A(1:matrix%n+1)  = [ matrix%A(1:japos-1),&
               &val,&
               &matrix%A(japos:matrix%n)]


               matrix%IA(irow+1:matrix%nr+1) =&
               &matrix%IA(irow+1:matrix%nr+1) + 1


            endif updated_IA

            ! Reallocate and reassign A and JA into matrix

            matrix%n = matrix%n + 1

            ! Reassociate pointers for CSR format
            matrix%val  => matrix%A  (1:matrix%n)
            matrix%col  => matrix%JA (1:matrix%n)
            matrix%ptrB => matrix%IA (1:matrix%nr  )
            matrix%ptrE => matrix%IA (2:matrix%nr+1)

            ! Column extension just needs to update the matrix%nc
            ! val
            if (extend_columns) matrix%nc = jcol

         endif no_mods

      endif already_init

   end subroutine add_value

   !   *********************************************************
   !   **  Adds one whole line to the sparse matrix represen- **
   !   **  tation given the row and column indices in the     **
   !   **  full format. Values possibly already present in the**
   !   **  row are overwritten.                               **
   !   **                                                     **
   !   **   Author:      Federico Perini                      **
   !   **   Last update: thursday, 26/07/2012                 **
   !   *********************************************************

   subroutine add_line(matrix,irow,val)
      implicit none

      type(sparse),     intent(inout), target      :: matrix
      integer,          intent(in)                 :: irow
      real (dp)       , dimension(:),  intent(in)  :: val

      integer                         :: j, old_nr, irow_last_new,&
      &irow_first, irow_last


      integer                                     :: eff_n
      integer,          dimension(:), allocatable :: eff_jcol
      real (dp)       , dimension(:), allocatable :: eff_val
      integer,          dimension(:), allocatable :: IAtmp, JAtmp
      real (dp)       , dimension(:), allocatable :: Atmp
      real (dp)       , parameter                 :: zero = 0.0_dp

      ! Gather number of elements and their effective values and
      ! column indices
      eff_n = count(val /= zero)

      ! Do nothing if values are all zero
      if (eff_n == 0) return

      ! Gather effective elements and their column indices
      allocate(eff_val(eff_n), eff_jcol(eff_n))
      eff_val  = pack(val,                 val /= zero)
      eff_jcol = pack([(j,j=1,size(val))], val /= zero)

      ! If matrix was empty, it has to be initialised to contain
      ! just that line
      already_init: if (.not.sparse_allocated(matrix) ) then

         call allocate(irow,  maxval(eff_jcol),&
         &eff_n, matrix          )

         matrix%A(1:eff_n)  = eff_val
         matrix%IA = [(1, j=1,irow), 1+eff_n]
         matrix%JA(1:eff_n) = eff_jcol

      else

         ! Store current matrix data
         allocate(JAtmp(matrix%n), Atmp(matrix%n))
         allocate(IAtmp(matrix%nr+1))
         IAtmp = matrix%IA
         JAtmp = matrix%JA(1:matrix%n)
         Atmp  = matrix%A (1:matrix%n)

         ! Dimensions of this matrix already contain row irow?
         extend_matrix: if (matrix%n == 0) then

            call allocate(matrix%nr,&
            &max(maxval(eff_jcol),matrix%nc),&
            &eff_n, matrix          )

            matrix%A(1:eff_n)  = eff_val
            matrix%IA = [ (1, j=1,irow),&
            &(1+eff_n, j=irow+1,matrix%nr+1) ]
            matrix%JA(1:eff_n) = eff_jcol

         elseif (irow > matrix%nr) then

            ! Old number of rows
            old_nr = matrix%nr

            ! Add a row to that matrix(updates all matrix% values)
            call allocate(irow,&
            &max(matrix%nc,maxval(eff_jcol)),&
            &matrix%n + eff_n,&
            &matrix)


            matrix%A    = [Atmp , eff_val]
            matrix%JA   = [JAtmp, eff_jcol]
            matrix%IA(1:old_nr+1)    = IAtmp
            matrix%IA(old_nr+1:irow) = IAtmp(old_nr+1)
            matrix%IA(irow+1)        = IAtmp(old_nr+1)+eff_n


         else

            ! Matrix already contains that row

            ! Gather current number of row elements
            irow_first = matrix%IA(irow)
            irow_last  = matrix%IA(irow+1)-1

            ! Add a row to that matrix(updates all matrix% values)
            call allocate(matrix%nr,&
            &max(matrix%nc,maxval(eff_jcol)),&
            &matrix%n + eff_n, matrix)

            ! 1) Previous lines
            matrix%A (1:irow_first-1) =  Atmp(1:irow_first-1)
            matrix%JA(1:irow_first-1) = JAtmp(1:irow_first-1)
            matrix%IA(1:irow)         = IAtmp(1:irow)

            ! 2) Inserted line
            irow_last_new = irow_first + eff_n - 1
            matrix%A (irow_first:irow_last_new) = eff_val
            matrix%JA(irow_first:irow_last_new) = eff_jcol
            matrix%IA(irow+1) = matrix%IA(irow) + eff_n

            ! 3) Remaining lines
            matrix%A(irow_last_new+1:matrix%n )= Atmp(irow_first:)
            matrix%JA(irow_last_new+1:matrix%n)=JAtmp(irow_first:)

            matrix%IA(irow+2:matrix%nr+1) =&
            &IAtmp(irow+2:matrix%nr+1) + eff_n&
            &- (irow_last - irow_first + 1)


         endif extend_matrix

         deallocate(IAtmp, JAtmp, Atmp)

      endif already_init

      deallocate(eff_val, eff_jcol)

   end subroutine add_line


   !   *********************************************************
   !   **  Extracts all the (i,j) indices of the nonzero      **
   !   **  elements from the sparse matrix representation,    **
   !   **  in columnwise order                                **
   !   **  Optional: also the values are exported             **
   !   **   Author:      Federico Perini                      **
   !   **   Last update: wednesday, 25/07/2012                **
   !   *********************************************************

   subroutine extract_rowcol_indices_columnwise(mat,irow,jcol,&
   &vals)
      implicit none

      type(sparse),                     intent(in)    :: mat
      integer,          dimension(:),   intent(inout) :: irow,jcol
      real (dp)       , dimension(:),   intent(inout),&
      &optional      :: vals

      type(sparse)                    :: tmp

      integer                         :: i, j, k, icount

      if (size(irow)<mat%n.or.size(jcol)<mat%n) then
         write(*,*)'Insufficient space provided in extract_rowcol'
         write(*,*)'size(irow)',size(irow)
         write(*,*)'size(jcol)',size(jcol)
         stop
      endif
      if (present(vals)) then
         if (size(vals)<mat%n) then
            write(*,*)'Insufficient space provided in extract_rowcol'
            write(*,*)'size(vals)',size(vals)
            stop
         endif
      endif


      if ((.not.mat%n>0).or.(.not.sparse_allocated(mat))) return

      !   Compute the transposed matrix (CSR indices of the
      !   transposed are equal to the columnwise representation
      !   of the non-transposed matrix)
      tmp = sparse_transpose(mat)


      !   Row index of the transposed is column index of mat
      icount = 0
      matj: do j = 1, tmp%nr
         mati: do k = tmp%IA(j), tmp%IA(j+1)-1
            i           = tmp%JA(k)
            icount      = icount + 1

            if (icount > mat%n) exit matj

            irow(icount) = i
            jcol(icount) = j

            if (present(vals)) vals(icount) = tmp%A(k)

         end do mati
      end do matj

      if (icount /= tmp%n) then
         write(*,*)'allocation error in extract_rowcol_indices'
         write(*,*)'spmatrix has ',icount,' elements'
         write(*,*)'should be: ',tmp%n
         stop
      endif

      call spdeallocate(tmp)

   end subroutine extract_rowcol_indices_columnwise




   !   *********************************************************
   !   **  Remove one value from the sparse matrix            **
   !   **  given the row and column indices in the full       **
   !   **  format                                             **
   !   **                                                     **
   !   **   Author:      Federico Perini                      **
   !   **   Last update: monday, 23/07/2012                   **
   !   *********************************************************

   subroutine remove_value(matrix,irow,jcol)
      use utilities, only: change_size
      implicit none

      type(sparse),     intent(inout), target :: matrix
      integer,          intent(in)            :: irow,jcol

      integer                         :: jaorig
      logical                         :: already_present
      integer,          dimension(:), allocatable :: JAtmp
      real (dp)       , dimension(:), allocatable :: Atmp


      ! If matrix was empty, it remains empty
      if (.not.sparse_allocated(matrix)) return

      ! Value already exists?
      jaorig = extract_sparse_index(matrix,irow,jcol)
      already_present = jaorig /= -1

      ! If value wasn't present, return
      if (.not.already_present) return

      ! Prepare temporary arrays of column indices and values
!           allocate(JAtmp(matrix%n-1))
!           allocate(Atmp (matrix%n-1))

      matrix%JA(1:matrix%n-1) =&
      &[matrix%JA(1:jaorig-1),matrix%JA(jaorig+1:matrix%n)]
      matrix%A (1:matrix%n-1) =&
      &[matrix% A(1:jaorig-1),matrix% A(jaorig+1:matrix%n)]

      ! Fix the row array by removing the current element to the
      ! count
      matrix%IA(irow+1:matrix%nr+1) =&
      &matrix%IA(irow+1:matrix%nr+1) - 1

      ! Update number of elements
      matrix%n = matrix%n - 1

      ! Safety, but should never be needed
!           if (size(matrix%A)<matrix%n)
!     &        call change_size(matrix%A,matrix%n)

!           if (size(matrix%JA)<matrix%n)
!     &        call change_size(matrix%JA,matrix%n)

!           deallocate(matrix%A, matrix%JA)
!           allocate(matrix%A(matrix%n), matrix%JA(matrix%n))
!           matrix%A  = Atmp
!           matrix%JA = JAtmp

      ! Reassociate pointers for CSR format
      matrix%val  => matrix%A  (1:matrix%n)
      matrix%col  => matrix%JA (1:matrix%n)
      matrix%ptrB => matrix%IA (1:matrix%nr  )
      matrix%ptrE => matrix%IA (2:matrix%nr+1)

      ! Matrix is empty now?
      if (matrix%n==0) call sparse_nullify_general(matrix)

!           deallocate(Atmp, JAtmp)

   end subroutine remove_value

   !   *********************************************************
   !   **  Remove one entire line from the sparse matrix      **
   !   **  given the row and column indices in the full       **
   !   **  format                                             **
   !   **                                                     **
   !   **   Author:      Federico Perini                      **
   !   **   Last update: monday, 23/07/2012                   **
   !   *********************************************************

   subroutine remove_line(matrix,irow)
      implicit none

      type(sparse),     intent(inout), target :: matrix
      integer,          intent(in)            :: irow

      integer                         :: nelems, ifirst, ilast

      integer,          dimension(:), allocatable :: JAtmp
      real (dp)       , dimension(:), allocatable :: Atmp

      character(len=*), parameter ::&
      &fmt_er = "(' Negative number of element in spmatrix: ',I5)",&
      &fmt_wrn= "(' Warning: deleted all intems in spmatrix ')"


      ! If matrix was empty, it remains empty
      if (.not.sparse_allocated(matrix)) return

      ! How many elements are there in the line?
      nelems = matrix%IA(irow+1) - matrix%IA(irow)

      ! If line was already empty, return
      if (nelems == 0) return

      ! Prepare temporary arrays of column indices and values
      allocate(JAtmp(matrix%n-nelems))
      allocate(Atmp (matrix%n-nelems))

      ! Compute first and last element of the row
      ifirst = matrix%IA(irow)
      ilast  = matrix%IA(irow+1)-1

      JAtmp = [matrix%JA(1:ifirst-1),matrix%JA(ilast+1:matrix%n)]
      Atmp  = [matrix% A(1:ifirst-1),matrix% A(ilast+1:matrix%n)]

      deallocate(matrix%A, matrix%JA)

      ! Fix the row array by removing the number of elements in
      ! row irow to the count
      matrix%IA(irow+1:matrix%nr+1) =&
      &matrix%IA(irow+1:matrix%nr+1) - nelems

      ! Update number of elements
      matrix%n = matrix%n - nelems

      if (matrix%n<0) then
         write(*,fmt_er)matrix%n
         stop
      endif

      if (matrix%n==0) then
         write(*,fmt_wrn)
         stop
      endif

      allocate(matrix%A(matrix%n), matrix%JA(matrix%n))
      matrix%A(1:matrix%n)  = Atmp
      matrix%JA(1:matrix%n) = JAtmp

      ! Reassociate pointers for CSR format
      matrix%val  => matrix%A  (1:matrix%n)
      matrix%col  => matrix%JA (1:matrix%n)
      matrix%ptrB => matrix%IA (1:matrix%nr  )
      matrix%ptrE => matrix%IA (2:matrix%nr+1)

      ! Matrix is empty now?
      if (matrix%n==0) call sparse_nullify_general(matrix)

      deallocate(Atmp, JAtmp)

   end subroutine remove_line


   !   *********************************************************
   !   **  Print sparse matrix pattern to file                **
   !   **                                                     **
   !   **   Author:      Federico Perini                      **
   !   **   Last update: wednesday, 28/03/2012                **
   !   *********************************************************

   subroutine print_sparsity_to_file(spmat,filename)
      implicit none

      type(sparse),     intent(in)                   :: spmat
      character(len=*), intent(in)                   :: filename
      real (dp)       , dimension(:,:), allocatable  :: dense
      integer                                        :: i, j
      character(len=*), parameter                    ::&
      &fmt_er  = "(' Unable to print unallocated matrix ',&
      &                    'in print_sparsity_to_file')",&
      &fmt_out = "(10000(1x,I2))"

      ! Preliminary check
      if (.not.(spmat%nr/=0 .and. spmat%nc/=0)) then
         write(*,fmt_er)
         stop
      endif

      ! Allocate dense matrix
      allocate(dense(spmat%nr,spmat%nc))
      call sparse_to_dense(spmat,dense)

      ! Normalize to 1
      where (dense/=0.e0_dp)
         dense = 1.e0_dp
      elsewhere
         dense = 0.e0_dp
      end where

      open(unit = 1000, file = filename)

      rows: do i = 1, spmat%nr
         write(1000,fmt_out)(int(dense(i,j)),j=1,spmat%nc)
      end do rows

      close(1000)
      deallocate(dense)

   end subroutine print_sparsity_to_file
!          ************************************************************

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
      character(len=*), parameter                    ::&
      &fmt_er  = "(' Unable to print unallocated matrix ',&
      &                    'in print_sparse_to_file')",&
      &fmt_out = "(10000(1x,E15.8))"

      ! Preliminary check
      if (.not.(spmat%nr/=0 .and. spmat%nc/=0)) then
         write(*,fmt_er)
         stop
      endif

      ! Allocate dense matrix
      allocate(dense(spmat%nr,spmat%nc))
      call sparse_to_dense(spmat,dense)

      open(unit = 1000, file = filename)

      rows: do i = 1, spmat%nr
         write(1000,fmt_out)(dense(i,j),j = 1,spmat%nc)
      end do rows

      close(1000)
      deallocate(dense)

   end subroutine print_sparse_to_file
!          ************************************************************

!     ******************************************************************
!     ** DEALLOCATION and INITIALISATION -related routines            **
!     ******************************************************************
   subroutine sparse_nullify_general(spmatrix)
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

   end subroutine sparse_nullify_general


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
      if ( associated(spmatrix%val ) .and.&
      &associated(spmatrix%col ) .and.&
      &associated(spmatrix%ptrB) .and.&
      &associated(spmatrix%ptrE) .and.&
      &allocated (spmatrix%A   ) .and.&
      &allocated (spmatrix%JA  ) .and.&
      &allocated (spmatrix%IA  ) .and.&
      &spmatrix%n  <= size(spmatrix%A ) .and.&
      &spmatrix%n  <= size(spmatrix%JA) .and.&
      &spmatrix%nr == size(spmatrix%IA) - 1 .and.&
      &spmatrix%nc >  0  )&
      &isallocated = .true.


   end function sparse_allocated



   subroutine compute_sparse_span(dense, nels)
      implicit none

      real (dp)       , dimension(:,:), intent(in)  :: dense
      integer,                          intent(out) :: nels

!          n = number of matrix rows plus 1
      if (size(dense,1)<1) then
         write(*,*)'Error: 0-row matrix in sparse format, n=',&
         &size(dense,1)
         stop
      endif

!          Sparse elements are all nonzero elements in the matrix,
!          plus diagonal elements even if zero

!           do i = 1, min(size(dense,1), size(dense,2))
!             if (dense(i,i)==0.e0_dp) nels = nels + 1
!           end do

      nels = count(dense/=0.e0_dp)
!           do i = 1, min(size(dense,1),size(dense,2))
!             if (dense(i,i)==0.e0_dp) nels = nels + 1
!           end do


   end subroutine compute_sparse_span
!          ************************************************************

   subroutine sparse_nullify_int(spmatrix)
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

   end subroutine sparse_nullify_int



   ! Convert dense matrix into sparse Yale format
   subroutine dense_to_sparse(dense,spmatrix)
      implicit none

      real (dp)       , dimension(:,:), intent(in)    :: dense
      type(sparse), intent(out), target               :: spmatrix
      integer :: m,n,i,j,iia,nA, nIA, ntot,c, nr, nc, nemptyprev
      logical :: isfirst, updateprevious

      character(len=*), parameter ::&
      &fmt_ersize = "(' Error: non-2D matrix made sparse ',&
      &                      ' in dense_to_sparse:')"
      m = size(dense,1)
      n = size(dense,2)

!            Halt on non-2D matrix
      if (m*n/=size(dense)) then
         write(*,fmt_ersize)
         write(*,*)'nrow = ',m,'ncol = ',n
         stop
      endif

!            Initialize sparse matrix array
      call compute_sparse_span(dense, nA)
      nIA = size(dense,1)
      nr  = size(dense,1)
      nc  = size(dense,2)
      call allocate(nr, nc, nA, spmatrix)

!            Store data into sparse format
      ntot = 0
      iia  = 0
      nemptyprev = 0
      row: do i = 1, spmatrix%nr
         isfirst = .true.
         c = count(dense(i,:)/=0.e0_dp)

         emptyrow: if (c==0) then ! Special treatment for empty row!

            updateprevious = .true.
            nemptyprev = nemptyprev + 1

         else


            column: do j = 1, spmatrix%nc

               if (dense(i,j) /= 0.e0_dp) then

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

!            Storing last column length
      if (.not.updateprevious) then
         spmatrix%IA(m+1) = spmatrix%IA(m)&
         &+ count(dense(m,:) /= 0.e0_dp)
!                 if (m <= n) then
!                    if (dense(m,m)==0.e0_dp) then
!                       spmatrix%IA(m+1) = spmatrix%IA(m+1)+1
!                    endif
!                 endif
      else
         spmatrix%IA(m+1-nemptyprev:m+1) = ntot+1
      endif

!             if (m<=n) then
!               if (dense(m,m)==0.e0_dp)
!     &            spmatrix%IA(m+1) = spmatrix%IA(m+1)+1
!             endif
!            Correct storage check

      if (nA/=ntot.or.ntot/=spmatrix%IA(m+1)-1) then

         write(*,*)'sparse storage error: ntot     = ',ntot
         write(*,*)'                      nA       = ',nA
         write(*,*)'                      IA(nr+1) = ',&
         &spmatrix%IA(m+1)
         write(*,*)'matrix rows:   nr = ',spmatrix%nr
         write(*,*)'matrix cols:   nc = ',spmatrix%nc

         stop

      endif


   end subroutine dense_to_sparse

!     *****************************************************************
   subroutine dense_to_sparseint(dense,spmatrix)
      implicit none

      real (dp)       , dimension(:,:), intent(in) :: dense
      type(sparseint), intent(out), target         :: spmatrix
      integer :: m,n,i,j,iia,nA, nIA, ntot,c, nr, nc, nemptyprev
      logical :: isfirst, updateprevious

      character(len=*), parameter ::&
      &fmt_ersize = "(' Error: non-2D matrix made sparse ')"
      m = size(dense,1)
      n = size(dense,2)

!            Halt on non-2D matrix
      if (m*n/=size(dense)) then
         write(*,fmt_ersize)
         write(*,*)'nrow = ',m,'ncol = ',n
         stop
      endif

!            Initialize sparse matrix array
      call compute_sparse_span(dense, nA)
      nIA = size(dense,1)
      nr  = size(dense,1)
      nc  = size(dense,2)
      call sparseint_allocate_det(nr, nc, nA, spmatrix)


!            Store data into sparse format
      ntot = 0
      iia  = 0
      nemptyprev = 0
      row: do i = 1, m
         isfirst = .true.
         c = count(dense(i,:)/=0.e0_dp)

         if (c==0) then ! Special treatment for empty row!

            updateprevious = .true.
            nemptyprev = nemptyprev + 1

         else



            column: do j = 1, n

               if (dense(i,j) /= 0.e0_dp) then

                  ntot = ntot + 1

                  spmatrix%A(ntot)  = int(dense(i,j))
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

         endif

      end do row

!            Storing last column length
      if (.not.updateprevious) then
         spmatrix%IA(m+1) = spmatrix%IA(m)&
         &+ count(dense(m,:) /= 0.e0_dp)
!                 if (m <= n) then
!                    if (dense(m,m)==0.e0_dp) then
!                       spmatrix%IA(m+1) = spmatrix%IA(m+1)+1
!                    endif
!                 endif
      else
         spmatrix%IA(m+1-nemptyprev:m+1) = ntot+1
      endif

!            Correct storage check

      if (nA/=ntot.or.ntot/=spmatrix%IA(m+1)-1) then

         write(*,*)'sparse storage error: ntot     = ',ntot
         write(*,*)'                      nA       = ',nA
         write(*,*)'                      IA(nr+1) = ',&
         &spmatrix%IA(m+1)
         write(*,*)'matrix rows:   nr = ',spmatrix%nr
         write(*,*)'matrix cols:   nc = ',spmatrix%nc

         stop

      endif

   end subroutine dense_to_sparseint

   ! Allocate sparse matrix with correct dimensions - but leave
   ! content blank
   subroutine sparseint_allocate_det(nrows,ncols,nels,spmatrix)
      implicit none

      integer, intent(in) :: nrows, ncols, nels
      type(sparseint), intent(inout), target :: spmatrix

      character(len=*), parameter ::&
      &fmt_ersize = "(' Error: non-2D matrix made sparse ',&
      &                      ' in sparse_allocate_det:' )"

!            Halt on non-2D matrix
      if (nrows<1.or.nels<1) then
         write(*,fmt_ersize)
         write(*,*)'nrow=',nrows,'nels=',nels
         stop
      endif

!            NB DO NOT CHANGE ARRAY ALLOCATION OUT OF THE ALLOCATION
!            ROUTINES!

!            Initialize sparse matrix array
      call sparse_nullify_int(spmatrix)

!            Assign data
      spmatrix%n  = nels
      spmatrix%nr = nrows
      spmatrix%nc = ncols

      allocate(spmatrix%A (nels)     )
      allocate(spmatrix%IA(nrows + 1))
      allocate(spmatrix%JA(nels)     )

      spmatrix%A  = 0
      spmatrix%IA = 0
      spmatrix%JA = 0

!            Associate pointers for CSR format
      spmatrix%val  => spmatrix%A  (1:nels)
      spmatrix%col  => spmatrix%JA (1:nels)
      spmatrix%ptrB => spmatrix%IA (1:nrows  )
      spmatrix%ptrE => spmatrix%IA (2:nrows+1)

   end subroutine sparseint_allocate_det

   ! **********************************************************
   ! **                                                      **
   ! ** Convert sparse Yale format into dense matrix         **
   ! **                                                      **
   ! **********************************************************
   subroutine sparse_to_dense(spmatrix,dense)
      implicit none
      type(sparse),                     intent(in)      :: spmatrix
      real (dp)       , dimension(spmatrix%nr,spmatrix%nc),&
      &intent(out)     :: dense

      integer :: i

      dense = 0.e0_dp

      if ( size(dense,1) > size(spmatrix%IA) ) then
         write(*,*) 'Matrix size error in sparse_to_dense'
         write(*,*) 'nrow dense : ',size(dense,1)
         write(*,*) 'nrow sparse: ',size(spmatrix%IA)
         stop
      endif

!          Assign data
      do i = 1, spmatrix%nr
         dense(i,spmatrix%JA(spmatrix%IA(i):spmatrix%IA(i+1)-1)) =&
         &spmatrix%A (spmatrix%IA(i):spmatrix%IA(i+1)-1)
      end do

   end subroutine sparse_to_dense

   ! **********************************************************
   ! **                                                      **
   ! ** Convert sparse Yale format into dense matrix         **
   ! **                                                      **
   ! **********************************************************

   ! **********************************************************
   ! The dense matrix is filled in columnwise order exploiting
   ! transposition properties
   ! tr = logical value. If .true., spmatrix already contains
   !      the transposed values for direct column designation
   ! **********************************************************

   subroutine sparse_to_dense_columnwise(tr,spmatrix,dense)
      implicit none
      logical,                          intent(in)    :: tr
      type(sparse),                     intent(in)    :: spmatrix
      real (dp)       , dimension(spmatrix%nr,spmatrix%nc),&
      &intent(out)   :: dense
      type(sparse)                                    :: sptransp
      integer :: i

      dense = 0.e0_dp

      if ( size(dense,1) > size(spmatrix%IA) ) then
         write(*,*) 'Matrix size error in sparse_to_dense'
         write(*,*) 'nrow dense : ',size(dense,1)
         write(*,*) 'nrow sparse: ',size(spmatrix%IA)
         stop
      endif

      if (.not.tr) then


         sptransp = sparse_transpose(spmatrix)

         do i = 1, spmatrix%nc
            dense(sptransp%JA(sptransp%IA(i):sptransp%IA(i+1)-1),i) =&
            &sptransp%A (sptransp%IA(i):sptransp%IA(i+1)-1)
         end do

         call sparse_nullify_general(sptransp)

      else

         do i = 1, spmatrix%nr
            dense(spmatrix%JA(spmatrix%IA(i):spmatrix%IA(i+1)-1),i) =&
            &spmatrix%A (spmatrix%IA(i):spmatrix%IA(i+1)-1)
         end do


      endif

   end subroutine sparse_to_dense_columnwise



   ! Allocate sparse matrix with correct dimensions - but leave
   ! content blank
   subroutine sparse_allocate(dense,spmatrix)
      implicit none

      real (dp)       , dimension(:,:), intent(in) :: dense
      type(sparse), intent(out), target :: spmatrix
      integer :: m, n, j

      character(len=*), parameter ::&
      &fmt_ersize = "(' Error: non-2D matrix made sparse ',&
      &                       ' in sparse-allocate: ')"
      m = size(dense,1)
      n = size(dense,2)

!            Halt on non-2D matrix
      if (m*n/=size(dense)) then
         write(*,fmt_ersize)
         write(*,*)'nrow = ',m,' ncol = ',n
         stop
      endif

!            NB DO NOT CHANGE ARRAY ALLOCATION OUT OF THE ALLOCATION
!            ROUTINES!

!            Initialize sparse matrix array
      call sparse_nullify_general(spmatrix)

!            Assign data
      if (count(dense/=0.e0_dp) == 0) then
         spmatrix%n = m*n
      else
         spmatrix%n = count(dense/=0.e0_dp)

         do j = 1, min(m,n)
            if (dense(j,j)==0.e0_dp) spmatrix%n = spmatrix%n + 1
         end do

      endif

      spmatrix%nr = m
      spmatrix%nc = n

      allocate(spmatrix%A (spmatrix%n))
      allocate(spmatrix%IA(m + 1)     )
      allocate(spmatrix%JA(spmatrix%n))

      spmatrix%A  = 0.e0_dp
      spmatrix%IA = 0
      spmatrix%JA = 0

!            Associate pointers for CSR format
      spmatrix%val  => spmatrix%A  (1:spmatrix%n)
      spmatrix%col  => spmatrix%JA (1:spmatrix%n)
      spmatrix%ptrB => spmatrix%IA (1:spmatrix%nr  )
      spmatrix%ptrE => spmatrix%IA (2:spmatrix%nr+1)


   end subroutine sparse_allocate

   ! ***********************************************************
   ! ** Multiply whole sparse matrix by a real number         **
   ! ***********************************************************

   function sparse_by_real(val,mat) result(mat2)
      implicit none

      real (dp)       , intent(in) :: val
      type(sparse),     intent(in) :: mat
      type(sparse)                 :: mat2

      if (val==0.e0_dp) then
         ! Multiplying by zero leads to empty matrix
         call sparse_nullify_general(mat2)
         return
      else

         ! Allocate matrix and setup data
         call allocate(mat%nr,mat%nc,mat%n,mat2)

         mat2%A  = mat%A * val
         mat2%IA = mat%IA
         mat2%JA = mat%JA

      endif

   end function sparse_by_real

   ! ***********************************************************
   ! ** Compute negative of a sparse matrix                   **
   ! ***********************************************************

   function sparse_neg(mat) result(mat2)
      implicit none

      type(sparse),     intent(in) :: mat
      type(sparse)                 :: mat2

      ! Allocate matrix and setup data
      call allocate(mat%nr,mat%nc,mat%n,mat2)

      mat2%A  = - mat%A
      mat2%IA = mat%IA
      mat2%JA = mat%JA

   end function sparse_neg

   ! Allocate sparse matrix with correct dimensions - but leave
   ! content blank
   subroutine sparse_allocate_det(nrows, ncols, nels, spmatrix)
      implicit none

      integer, intent(in) :: nrows, ncols, nels
      type(sparse), intent(inout), target :: spmatrix
      logical :: do_A, do_IA, do_JA

      character(len=*), parameter ::&
      &fmt_ersize = "(' Error: non-2D matrix made sparse ',&
      &                      ' in sparse_allocate_det:' )"

!            Halt on non-2D matrix
      if (nrows<1.or.ncols<1) then
         write(*,fmt_ersize)
         write(*,*)'nrow=',nrows,'nels=',nels,'ncols=',ncols
         stop
      endif

!            NB DO NOT CHANGE ARRAY ALLOCATION OUT OF THE ALLOCATION
!            ROUTINES!
!
!             write(*,*)'nels',nels,size(spmatrix%A)
!             write(*,*)'nrows',nrows,size(spmatrix%IA)
!             write(*,*)'ncols',ncols,size(spmatrix%JA)

!            Check if matrix is already allocated
      do_A  = size(spmatrix%A )<nels
      do_IA = size(spmatrix%IA)/=nrows+1
      do_JA = size(spmatrix%JA)<nels
      do_A  = do_A  .or..not.allocated(spmatrix%A )
      do_IA = do_IA .or..not.allocated(spmatrix%IA)
      do_JA = do_JA .or..not.allocated(spmatrix%JA)


!            Initialize sparse matrix array
      if (do_A.and.do_IA.and.do_JA)&
      &call sparse_nullify_general(spmatrix)

!            Assign data
      spmatrix%n  = nels
      spmatrix%nr = nrows
      spmatrix%nc = ncols

      if (do_A ) then
         if (allocated(spmatrix%A))deallocate(spmatrix%A)
         allocate(spmatrix%A (nels)     )
      endif
      if (do_IA) then
         if (allocated(spmatrix%IA))deallocate(spmatrix%IA)
         allocate(spmatrix%IA(nrows + 1))
      endif
      if (do_JA) then
         if (allocated(spmatrix%JA))deallocate(spmatrix%JA)
         allocate(spmatrix%JA(nels)     )
      endif

      spmatrix%A  = 0.e0_dp
      spmatrix%IA = 0
      spmatrix%JA = 0

!            Associate pointers for CSR format
      spmatrix%val  => spmatrix%A  (1:nels)
      spmatrix%col  => spmatrix%JA (1:nels)
      spmatrix%ptrB => spmatrix%IA (1:nrows  )
      spmatrix%ptrE => spmatrix%IA (2:nrows+1)

   end subroutine sparse_allocate_det


!          ************************************************************
!          ** Sparse matrix - dense vector multiplication, C = A*b   **
!          ** operating with transposed matrix, AT
!          ************************************************************
   function sparse_matmulT(AT,b) result(C)
      implicit none

      type(sparse),                      intent(in) :: AT
      real (dp)       , dimension(AT%nr), intent(in) :: b
      real (dp)       , dimension(AT%nc)             :: C

      integer, parameter     :: maxcache = 32000, cacheline = 4
      integer                :: bcol, blkcols, blkcolb, blkcole,&
      &blknnz, i, ii, j, je, jb

      if((.not.associated(AT%val))  .or.&
      &(.not.associated(AT%col))  .or.&
      &(.not.associated(AT%ptrB)) .or.&
      &(.not.associated(AT%ptrE))) then
         write(*,*)' Unassociated sparse array in sparse_matmulT '
         write(*,*)' sparse matrix: ',AT%nr,AT%nc,AT%n
         write(*,*)' dense  matrix: ',size(b)
         stop
      endif

!          Check array dimensions
      if (AT%nr /= size(AT%IA) -1) then
         write(*,*)' Wrong array dimensions in sparse_matmulT, ',&
         &' spA = ',size(AT%IA)-1,' C ',size(C)
         stop
      endif

!          Initialize C
      C(:) = 0.e0_dp

!          Calculate number of block columns of A (not AT),
!          based on CACHELINE:
      bcol = cacheline
      blkcols = AT%nr/bcol ! division among integers!
      if ( blkcols*bcol /= AT%nr ) blkcols = blkcols + 1

!          Loop through blkrows block columns:
      loop_block_cols: do i = 1 , blkcols

         blkcolb = (i-1)*bcol + 1
         blkcole = blkcolb + bcol - 1
         if (blkcole >= AT%nr ) blkcole = AT%nr

!             Count the number of nonzeros in this block column:
!              blknnz = 0
!              do j = blkcolb,blkcole
!                 blknnz = blknnz + AT%IA(j+1) - AT%IA(j)
!              end do


!             Now, loop through the only column rhscols block of c & b:
!             Loop through the brow columns in this block:
         loop_block_columns: do j = blkcolb,blkcole
            jb = AT%IA(j)
            je = AT%IA(j+1)

            if (b(j)==0.e0_dp) cycle loop_block_columns

            daxpyi: do ii = jb, je-1
               C(AT%JA(ii)) = C(AT%JA(ii)) + b(j) * AT%A(ii)
            end do daxpyi

         end do loop_block_columns

      end do loop_block_cols



   end function sparse_matmulT



!          ************************************************************
!          ** Elevate all the elements in an array to the powers     **
!          ** contained in a sparse matrix, whose number of rows     **
!          ** is the same as the dimensions of the array:            **
!          ** sparse_result(i,j) = a(i)**sparse_matrix(i,j)          **
!          ************************************************************
   function column_sparse_power(A,b) result(C)
      implicit none

      type(sparse),                      intent(in) :: A
      real (dp)       , dimension(A%nr), intent(in) :: b
      type(sparse)                                  :: C

      integer, parameter     :: maxcache = 32000, cacheline = 4
      integer                :: brow, blkrows, blkrowb, blkrowe,&
      &blknnz, i, ii, j

      if((.not.associated(A%val))  .or.&
      &(.not.associated(A%col))  .or.&
      &(.not.associated(A%ptrB)) .or.&
      &(.not.associated(A%ptrE))) then
         write(*,*)' Unassociated sparse array in sparse_matmul '
         write(*,*)' sparse matrix: ',A%nr,A%nc,A%n
         write(*,*)' dense  array : ',size(b)
         stop
      endif

!          Check array dimensions
      if (A%nr /= size(b) ) then
         write(*,*)' Wrong array dimensions in col_sparse_power ',&
         &' spA = ',A%nr,' b ',size(b)
         stop
      endif

!          Allocate result matrix
      call allocate(A%nr, A%nc, A%n, C)

!          Calculate number of block rows, based on CACHELINE:
      brow = cacheline
      blkrows = A%nr/brow ! division among integers!
      if ( blkrows*brow /= A%nr ) blkrows = blkrows + 1

!          Loop through blkrows block rows:
!
      loop_block_rows: do i = 1 , blkrows

         blkrowb = (i-1)*brow + 1
         blkrowe = blkrowb + brow - 1
         if (blkrowe >= A%nr ) blkrowe = A%nr

!             Count the number of nonzeros in this block of rows:
!             Now, loop through the only column rhscols block of c & b:
!             Loop through the brow rows in this block:
         loop_block_columns: do j = blkrowb,blkrowe
            rowvals: do ii = A%IA(j), A%IA(j+1)-1
               C%A(ii) = b(j)**A%A(ii)
            end do rowvals
         end do loop_block_columns

      end do loop_block_rows



   end function column_sparse_power



!          ************************************************************
!          ** Sparse matrix - dense vector multiplication, C = A*b   **
!          ************************************************************
   function sparseint_matmul(A,b) result(C)
      implicit none

      type(sparseint),                   intent(in) :: A
      real (dp)       , dimension(A%nc), intent(in) :: b
      real (dp)       , dimension(A%nr)             :: C

      integer, parameter     :: maxcache = 32000, cacheline = 8
      integer                :: brow, blkrows, blkrowb, blkrowe,&
      &blknnz, i, ii, j


!          Check array dimensions
!           lspwork = max(A%nr,A%nc)

      if((.not.associated(A%val))  .or.&
      &(.not.associated(A%col))  .or.&
      &(.not.associated(A%ptrB)) .or.&
      &(.not.associated(A%ptrE))) then
         write(*,*)' Unassociated sparse array in sparse_matmul '
         write(*,*)' sparse matrix: ',A%nr,A%nc,A%n
         write(*,*)' dense  matrix: ',size(b)
         stop
      endif

!          Check array dimensions
      if (A%nr /= size(A%IA) -1) then
         write(*,*)' Wrong array dimensions in sparse_matmul, ',&
         &' spA = ',size(A%IA)-1,' C ',size(C)
         stop
      endif

!          Calculate number of block rows, based on CACHELINE:
      brow = cacheline
      blkrows = A%nr/brow ! division among integers!
      if ( blkrows*brow /= A%nr ) blkrows = blkrows + 1

!          Loop through blkrows block rows:
!
      loop_block_rows: do i = 1 , blkrows

         blkrowb = (i-1)*brow + 1
         blkrowe = blkrowb + brow - 1
         if (blkrowe >= A%nr ) blkrowe = A%nr

!             Count the number of nonzeros in this block of rows:
!              blknnz = 0
!              do j = blkrowb,blkrowe
!                 blknnz = blknnz + A%IA(j+1) - A%IA(j)
!              end do

!             Now, loop through the only column rhscols block of c & b:
!             Loop through the brow rows in this block:
         loop_block_columns: do j = blkrowb,blkrowe
            c(j) = 0.e0_dp
            ddoti: do ii = A%IA(j), A%IA(j+1)-1
               c(j) = c(j) + A%A(ii) * b(A%JA(ii))
            end do ddoti


         end do loop_block_columns

      end do loop_block_rows



   end function sparseint_matmul



!          ************************************************************
!          ** Sparse matrix transposition                            **
!          ************************************************************
   function sparse_transpose(A) result(AT)
      implicit none

      type(sparse), intent(in)    :: A
      type(sparse)                :: AT

      integer                    :: nels, i, j, idx

      if((.not.associated(A%val))  .or.&
      &(.not.associated(A%col))  .or.&
      &(.not.associated(A%ptrB)) .or.&
      &(.not.associated(A%ptrE))) then
         write(*,*)'Unassociated sparse array in sparse_transpose'
         stop
      endif

      nels = A%n!size(A%A)

      if (.not.(nels>0)) then
         write(*,*)'Error in sparse matrix allocation for transp'
         write(*,*)'elements',nels
         write(*,*)'Values',A%A
         stop
      endif

!          Reset matrix allocation
      call allocate(A%nc, A%nr, A%n, AT)

      AT%IA(1:A%nc+1) = 0
      AT%A (1:A%n   ) = 0.e0_dp

      AT%IA(1) = 1

!          Count indices for each column
      do i = 1, A%nr
         do j = A%IA(i), A%IA(i+1) - 1
            AT%IA(A%JA(j) + 1) = AT%IA(A%JA(j)+1)+1
         end do
      end do

      do i = 1, A%nc
         AT%IA(i + 1) = AT%IA(i) + AT%IA(i+1)
      end do

!          Build AT%JA array
      do i = 1, A%nr
         do j = A%IA(i), A%IA(i+1) - 1
            idx = A%JA(j)
            AT%JA(AT%IA(idx)) = i
            AT%A (AT%IA(idx)) = A%A(j)
            AT%IA(idx)        = AT%IA(idx)+1
         end do
      end do

!          Fix AT%IA array
      do i = A%nc, 2, -1
         AT%IA(i) = AT%IA(i-1)
      end do
      AT%IA(1) = 1

   end function sparse_transpose

!          ************************************************************
!          ** Sparse matrix transposition                            **
!          ************************************************************
   subroutine sparseint_transpose(A,AT)
      implicit none

      type(sparseint),                   intent(in)    :: A
      type(sparseint),                   intent(inout) :: AT

      integer                :: nels, i, j, idx
!           integer, intent(in)    :: nrow, ncol


      if((.not.associated(A%val))  .or.&
      &(.not.associated(A%col))  .or.&
      &(.not.associated(A%ptrB)) .or.&
      &(.not.associated(A%ptrE))) then
         write(*,*)'Unassociated sparse array in sparse_transpose'
         stop
      endif

      nels = A%n

      if (.not.(nels>0)) then
         write(*,*)'Error in sparse matrix allocation for transp'
         write(*,*)'elements',nels
         write(*,*)'Values',A%A
         stop
      endif

!          Reset matrix allocation
      call sparseint_allocate_det(A%nc, A%nr, A%n, AT)

      AT%IA(1:A%nc+1) = 0
      AT%A (1:A%n   ) = 0

      AT%IA(1) = 1

!          Count indices for each column
      do i = 1, A%nr
         do j = A%IA(i), A%IA(i+1) - 1
            AT%IA(A%JA(j) + 1) = AT%IA(A%JA(j)+1)+1
         end do
      end do

      do i = 1, A%nc
         AT%IA(i + 1) = AT%IA(i) + AT%IA(i+1)
      end do

!          Build AT%JA array
      do i = 1, A%nr
         do j = A%IA(i), A%IA(i+1) - 1
            idx = A%JA(j)
            AT%JA(AT%IA(idx)) = i
            AT%A (AT%IA(idx)) = A%A(j)
            AT%IA(idx)        = AT%IA(idx)+1
         end do
      end do

!          Fix AT%IA array
      do i = A%nc, 2, -1
         AT%IA(i) = AT%IA(i-1)
      end do
      AT%IA(1) = 1

   end subroutine sparseint_transpose



!          ************************************************************
!          ** Sparse matrix transposition                            **
!          ** Output: not a full sparse matrix but only the array    **
!          ** of the real values, not the sparsity structure         **
!          ************************************************************
   subroutine sparse_transpose_valuesonly(A,AT_values)
      implicit none

      type(sparse),                      intent(in)    :: A
      real (dp)       , dimension(A%n),  intent(out)   :: AT_values

      integer                    :: nels, i, j, idx
      integer, dimension(A%nc+1) :: IA

      if((.not.associated(A%val))  .or.&
      &(.not.associated(A%col))  .or.&
      &(.not.associated(A%ptrB)) .or.&
      &(.not.associated(A%ptrE))) then
         write(*,*)'Unassociated sparse array in sparse_transpose'
         stop
      endif

      nels = A%n

      if (.not.(nels>0)) then
         write(*,*)'Error in sparse matrix allocation for transp'
         write(*,*)'elements',nels
         write(*,*)'Values',A%A
         stop
      endif

      IA(1)            = 1
      IA(2:A%nc+1)     = 0
      AT_values(1:A%n) = 0.e0_dp

!          Count indices for each column
      do i = 1, A%nr
         do j = A%IA(i), A%IA(i+1) - 1
            IA(A%JA(j) + 1) = IA(A%JA(j)+1)+1
         end do
      end do

      do i = 1, A%nc
         IA(i + 1) = IA(i) + IA(i+1)
      end do

!          Build AT%JA array
      do i = 1, A%nr
         do j = A%IA(i), A%IA(i+1) - 1
            idx = A%JA(j)
            AT_values (IA(idx)) = A%A(j)
            IA(idx)        = IA(idx)+1
         end do
      end do

!          Fix AT%IA array
      do i = A%nc, 2, -1
         IA(i) = IA(i-1)
      end do
      IA(1) = 1

   end subroutine sparse_transpose_valuesonly


!          ************************************************************
!          ** Sparse matrix transposition                            **
!          ** Output: the array of permutations that the values A%A  **
!          ** need to have in order to become values of AT%A         **
!          ************************************************************
   subroutine sparse_transpose_permutations(A,I_ATval)
      implicit none

      type(sparse),                      intent(in)    :: A
      integer, dimension(A%n),           intent(out)   :: I_ATval

      integer                    :: nels, i, j, idx
      integer, dimension(A%nc+1) :: IA

      if((.not.associated(A%val))  .or.&
      &(.not.associated(A%col))  .or.&
      &(.not.associated(A%ptrB)) .or.&
      &(.not.associated(A%ptrE))) then
         write(*,*)'Unassociated sparse array in sparse_transpose'
         stop
      endif

      nels = A%n

      if (.not.(nels>0)) then
         write(*,*)'Error in sparse matrix allocation for transp'
         write(*,*)'elements',nels
         write(*,*)'Values',A%A
         stop
      endif

      IA(1)            = 1
      IA(2:A%nc+1)     = 0
      I_ATval          = 0

!          Count indices for each column
      do i = 1, A%nr
         do j = A%IA(i), A%IA(i+1) - 1
            IA(A%JA(j) + 1) = IA(A%JA(j)+1)+1
         end do
      end do

      do i = 1, A%nc
         IA(i + 1) = IA(i) + IA(i+1)
      end do

!          Build AT%JA array
      do i = 1, A%nr
         do j = A%IA(i), A%IA(i+1) - 1
            idx = A%JA(j)
            I_ATval(j) = IA(idx)
            IA(idx)        = IA(idx)+1
         end do
      end do

   end subroutine sparse_transpose_permutations



!          ************************************************************
!          ** Symbolic sparse matrix multiplication: indexing C for  **
!          ** sparse(A) * sparse(B) = sparse(C)                      **
!          ************************************************************
   subroutine sparse_symbolic_mm(spA, spB, spC)
      implicit none

      type(sparse), intent(in)  :: spA, spB
      type(sparse), intent(out) :: spC
      type(sparse)              :: tmp

      character(len=*), parameter :: er1= "(' Unallocated matrix ',&
      &      'in symbolic sparse matmul: ',A1)",&
      &er2 = "(' wrong matrix dimension for sparse matmul: ( ',I4,&
      &            ' x ',I4,' ) * ( ',I4,' x ',I4,' ) ')"

      integer :: nrA, ncA, nrB, ncB, nrC, ncC, nelC,&
      &minmn, minlm, length, istart, i, j, k, jj
      integer, parameter :: diagA = 0, diagB = 0, diagC = 0
      integer, dimension(:), allocatable :: idx

!          Check on matrix dimensions!
      if (.not.allocated(spA%A))write(*,er1)'A'
      if (.not.allocated(spB%A))write(*,er1)'B'

      nrA = spA%nr
      nrB = spB%nr
      ncA = spA%nc
      ncB = spB%nc

      if (ncA/=nrB) write(*,er2)nrA,ncA,nrB,ncB

      nrC = nrA
      ncC = ncB

!          Prepare spC matrix for intent(out)
      call allocate(nrC,ncC,nrC*ncC,tmp)

!          Start the smmp symbolic matrix multiply routine
      allocate(idx(max(nrA,ncA,ncB)))
      idx(:) = 0

      tmp%IA(1) = 1

!          n = nrA
!          m = ncA
!          l = ncB
      minlm = min(ncB,ncA)
      minmn = min(ncA,nrA)

      main_loop: do i = 1, nrA

         istart = -1
         length = 0

         merge_row_lists: do jj = spA%IA(i), spA%IA(i+1)

!               a = d + ...
            if (jj == spA%IA(i+1)) then
               cycle merge_row_lists
            else
               j=spA%JA(jj)
            endif
!               b = d + ...
            b_length: do k = spB%IA(j),spB%IA(j+1)-1
               if(idx(spB%JA(k)) == 0) then
                  idx(spB%JA(k)) = istart
                  istart = spB%JA(k)
                  length = length + 1
               endif
            end do b_length


         end do merge_row_lists ! 30


!             row i of spC%JA
         tmp%IA(i+1) = tmp%IA(i) + length

         spC_rowi: do j = tmp%IA(i), tmp%IA(i+1)-1
            tmp%JA(j)      = istart
            istart         = idx(istart)
            idx(tmp%JA(j)) = 0
         end do spC_rowi
         idx(i)           = 0

      end do main_loop ! 50

      nelC = count(tmp%JA/=0)

      call allocate(nrC,ncC,nelC,spC)

      spC%A  = tmp%A(1:nelC)
      spC%IA = tmp%IA
      spc%JA = tmp%JA(1:nelC)

      call sparse_nullify_general(tmp)


      deallocate(idx)

   end subroutine sparse_symbolic_mm

!          ************************************************************
!          ** Sparse matrix * sparse matrix multiplication:          **
!          ** sparse(A) * sparse(B) = sparse(C)                      **
!          ** NB sparse(C) has prior to be indexed using the         **
!          **    symbolic routine sparse_symbolic_mm                 **
!          ************************************************************
   subroutine sparse_2_matmul(spA, spB, spC)
      use utilities, only: extend
      implicit none

      type(sparse), intent(in)    :: spA, spB
      type(sparse), intent(inout) :: spC

      integer, parameter :: cacheline = 4
      integer            :: maxdim, i, j, ii, jj,&
      &k, blkrows, blkrowb, blkrowe
      real (dp)        :: ajj

!          Check initialization
      if  ( (.not.(spC%n>0)).or.(.not.(spC%nr>0)) .or.&
      &(.not.(spC%nc>0)) ) then

         write(*,*) 'Sparse_2_matmul: need to initialize result ',&
         &' matrix sparsity pattern '
         stop

      endif

!          Check matrix dimensions
      if (spA%nr/=spC%nr .or. spA%nc/=spB%nr .or. spB%nc/=spC%nc)&
      &then

         write(*,*) 'Wrong matrix dimensions in sparse_2_matmul'
         write(*,*) 'Matrix A: ',spA%nr,' x ',spA%nc
         write(*,*) 'Matrix B: ',spB%nr,' x ',spB%nc
         write(*,*) 'Matrix C: ',spC%nr,' x ',spC%nc

         stop

      endif

!        Initialize matrix dimensions invariants
      maxdim = max(spA%nr,spA%nc,spB%nc)
!         minA   = min(spA%nr,spA%nc)
!         minB   = min(spB%nr,spB%nc)
!         minC   = min(spC%nr,spC%nc)


!        Initialize and check working array
      call extend(s2work, maxdim)
!           if (.not.allocated(s2work)) then
!              ls2work = max(spA%nr,spA%nc,spB%nc)
!              allocate(s2work(ls2work))
!              s2work = 0.e0_dp
!           endif
!           if (size(s2work)<max(spA%nr,spA%nc,spB%nc)) then
!              if(allocated(s2work))deallocate(s2work)
!              ls2work = max(spA%nr,spA%nc,spB%nc)
!              allocate(s2work(ls2work))
!              s2work = 0.e0_dp
!           endif



!        Initialize working array
      s2work(1:maxdim) = 0.e0_dp

!        Calculate number of block rows, based on CACHELINE:
      blkrows = spA%nr/cacheline ! division among integers!
      if ( blkrows*cacheline /= spA%nr ) blkrows = blkrows + 1

      loop_A_blocks: do ii = 1, blkrows

         blkrowb = (ii-1)*cacheline + 1
         blkrowe = blkrowb + cacheline - 1
         if (blkrowe >= spA%nr ) blkrowe = spA%nr

         loop_A_rows: do i = blkrowb,blkrowe!1, spA%nr

            loop_A_cols: do jj = spA%IA(i), spA%IA(i+1)

               ! a = d + ...
               if (jj == spA%IA(i+1)) cycle loop_A_cols
               j   = spA%JA(jj)
               ajj = spA%A (jj)

               ! b = d + ...
               loop_B_cols: do k = spB%IA(j),spB%IA(j+1)-1
!                  kk = spB%JA(k)
                  s2work(spB%JA(k)) = s2work(spB%JA(k)) + ajj*spB%A(k)
!                  s2work(kk) = s2work(kk) + ajj * spB%A(k)
               end do loop_B_cols

            end do loop_A_cols

            ! c = d + ...
            loop_C_elems: do j = spC%IA(i),spC%IA(i+1)-1

               spC%A(j) = s2work(spC%JA(j))
               s2work(spC%JA(j)) = 0.e0_dp

            end do loop_C_elems

         end do loop_A_rows

      end do loop_A_blocks

      return
   end subroutine sparse_2_matmul


!          ************************************************************
!          ** Sparse value gather                                    **
!          ** Extract a single value from a sparse matrix, given     **
!          ** row and column indices (i,j) of the dense form         **
!          ************************************************************

   function sparse_value(matrix,i,j) result(x)
      implicit none

      type(sparse), intent(in) :: matrix
      integer,      intent(in) :: i, j
      real (dp)                :: x
      integer                  :: i0, i1, k


!          Check for stupid errors
      if (i > matrix%nr .or. j > matrix%nc) then
         write(*,*)' (i,j) couple exceeds matrix dimensions '
         write(*,*)' in sparse_value; i,j=',i,j
         write(*,*)' matrix dims=',matrix%nr, matrix%nc
         stop
      endif

      x = 0.e0_dp

      i0 = matrix%IA(i)
      i1 = matrix%IA(i+1)-1

      set_value: do k = i0,i1
         if (matrix%JA(k) == j) then
            x = matrix%A(k)
            return
         endif
      end do set_value


   end function sparse_value




!          ************************************************************
!          ** Sparse square permutation                              **
!          ** Permutates matrix rows and colum indexes according to  **
!          ** the integer permutation array IP                       **
!          ** - Works only on square matrices                        **
!          ************************************************************

   function sparse_square_permutation(A,IP) result(B)
      implicit none

      type(sparse),             intent(in) :: A
      integer, dimension(:),    intent(in) :: IP
      type(sparse)                         :: B

      integer, dimension(A%nr)             :: IIP
      integer                              :: i, j, ir, ic, nel_B

!          Preliminary check
      if (A%nr /= A%nc) then
         write(*,*)'Cannot permutate non-square matrix'
         stop
      endif

      if (size(IP) /= A%nr) then
         write(*,*)'Permutation index array length is wrong'
         write(*,*)'length = ',size(IP), ' should be ',A%nr
         stop
      endif

!          Generate inverted pointer array:
!          If IP(i) = j,  then j = row index of A and i = row index of B
!            IIP(j) = i
      IIP(IP(1:A%nr)) = [(i,i=1,A%nr)]

!          Allocate B with the same number of elements as A
      call allocate(A%nr, A%nc, A%n, B)

!          Sparse matrices follow row-based indexing
      nel_b   = 0
      B%IA(1) = 1
      rows: do i = 1, B%nr

!                   Permutated row index
         ir = IP(i)

         cols: do j = A%IA(ir), A%IA(ir+1)-1



!                            Permutated column index
            ic = IIP(A%JA(j))

!                            Update counter of elements in matrix B
            nel_B = nel_B + 1

!                            Store element in B, and its column index
            B%A (nel_B) = A%A(j)
            B%JA(nel_B) = ic

         end do cols

!                    Update number of elements in current row
         B%IA(i+1) = nel_B + 1

      end do rows

   end function sparse_square_permutation

!          ************************************************************
!          ** Sparse row permutation                                 **
!          ** Permutates matrix rows indexes according to the integer**
!          ** permutation array IP                                   **
!          ************************************************************

   function sparse_row_permutation(A,IP) result(B)
      implicit none

      type(sparse),             intent(in) :: A
      integer, dimension(:),    intent(in) :: IP
      type(sparse)                         :: B

      integer, dimension(A%nr)             :: IIP
      integer                              :: i, j, ir, nel_B

      if (size(IP) /= A%nr) then
         write(*,*)'Permutation index array length is wrong'
         write(*,*)'in sparse_row_permutation'
         write(*,*)'length = ',size(IP), ' should be ',A%nr
         stop
      endif

!          Generate inverted pointer array:
!          If IP(i) = j,  then j = row index of A and i = row index of B
!            IIP(j) = i
      IIP(IP(1:A%nr)) = [(i,i=1,A%nr)]

!          Allocate B with the same number of elements as A
      call allocate(A%nr, A%nc, A%n, B)

!          Sparse matrices follow row-based indexing
      nel_b   = 0
      B%IA(1) = 1
      rows: do i = 1, B%nr

!                   Permutated row index
         ir = IP(i)

         cols: do j = A%IA(ir), A%IA(ir+1)-1

!                            Update counter of elements in matrix B
            nel_B = nel_B + 1

!                            Column data is not permutated
            B%A (nel_B) = A%A(j)
            B%JA(nel_B) = A%JA(j)

         end do cols

!                    Update number of elements in current row
         B%IA(i+1) = nel_B + 1

      end do rows

   end function sparse_row_permutation

   function sparseint_row_permutation(A,IP) result(B)
      implicit none

      type(sparseint),          intent(in) :: A
      integer, dimension(:),    intent(in) :: IP
      type(sparseint)                      :: B

      integer, dimension(A%nr)             :: IIP
      integer                              :: i, j, ir, nel_B

      if (size(IP) /= A%nr) then
         write(*,*)'Permutation index array length is wrong'
         write(*,*)'in sparseint_row_permutation'
         write(*,*)'length = ',size(IP), ' should be ',A%nr
         stop
      endif

!          Generate inverted pointer array:
!          If IP(i) = j,  then j = row index of A and i = row index of B
!            IIP(j) = i
      IIP(IP(1:A%nr)) = [(i,i=1,A%nr)]

!          Allocate B with the same number of elements as A
      call sparseint_allocate_det(A%nr, A%nc, A%n, B)

!          Sparse matrices follow row-based indexing
      nel_b   = 0
      B%IA(1) = 1
      rows: do i = 1, B%nr

!                   Permutated row index
         ir = IP(i)

         cols: do j = A%IA(ir), A%IA(ir+1)-1

!                            Update counter of elements in matrix B
            nel_B = nel_B + 1

!                            Column data is not permutated
            B%A (nel_B) = A%A(j)
            B%JA(nel_B) = A%JA(j)

         end do cols

!                    Update number of elements in current row
         B%IA(i+1) = nel_B + 1

      end do rows

   end function sparseint_row_permutation

!          ************************************************************
!          ** Sparse column permutation                              **
!          ** Permutates matrix column indexes according to the      **
!          ** integer permutation array IP                           **
!          ************************************************************

   function sparse_col_permutation(A,IP) result(B)
      implicit none

      type(sparse),             intent(in) :: A
      integer, dimension(:),    intent(in) :: IP
      type(sparse)                         :: B

      integer, dimension(A%nc)             :: IIP
      integer                              :: i, j, ir, ic

!          Preliminary check
      if (size(IP) /= A%nc) then
         write(*,*)'Permutation index array length is wrong'
         write(*,*)'in sparse_column_permutation'
         write(*,*)'length = ',size(IP), ' should be ',A%nc
         stop
      endif

!          Generate inverted pointer array:
!          If IP(i) = j,  then j = col index of A and i = col index of B
!            IIP(j) = i
      IIP(IP(1:A%nc)) = [(i,i=1,A%nc)]

!          Allocate B with the same number of elements as A
      call allocate(A%nr, A%nc, A%n, B)

!          Sparse matrices follow row-based indexing

!          The number of values per each row is unchanged by permutation
      B%IA = A%IA

      rows: do ir = 1, B%nr

         cols: do j = A%IA(ir), A%IA(ir+1)-1

!                            Permutated column index
            ic = IIP(A%JA(j))

!                            Store element in B, and its column index
            B%A (j) = A%A(j)
            B%JA(j) = ic

         end do cols

      end do rows

   end function sparse_col_permutation

   function sparseint_col_permutation(A,IP) result(B)
      implicit none

      type(sparseint),          intent(in) :: A
      integer, dimension(:),    intent(in) :: IP
      type(sparseint)                      :: B

      integer, dimension(A%nc)             :: IIP
      integer                              :: i, j, ir, ic

!          Preliminary check
      if (size(IP) /= A%nc) then
         write(*,*)'Permutation index array length is wrong'
         write(*,*)'in sparseint_column_permutation'
         write(*,*)'length = ',size(IP), ' should be ',A%nr
         stop
      endif

!          Generate inverted pointer array:
!          If IP(i) = j,  then j = col index of A and i = col index of B
!            IIP(j) = i
      IIP(IP(1:A%nc)) = [(i,i=1,A%nc)]

!          Allocate B with the same number of elements as A
      call sparseint_allocate_det(A%nr, A%nc, A%n, B)

!          Sparse matrices follow row-based indexing

!          The number of values per each row is unchanged by permutation
      B%IA = A%IA

      rows: do ir = 1, B%nr

         cols: do j = A%IA(ir), A%IA(ir+1)-1

!                            Permutated column index
            ic = IIP(A%JA(j))

!                            Store element in B, and its column index
            B%A (j) = A%A(j)
            B%JA(j) = ic

         end do cols

      end do rows

   end function sparseint_col_permutation



!          ************************************************************
!          ** Generate a identity matrix of order n in sparse form   **
!          ** beta is a real factor that multiplies the matrix       **
!          ************************************************************
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

!          Allocate space for the sparse matrix
      call allocate(n, n, n, I)

      if (present(beta)) then
         I%A(1:n)  = beta
      else
         I%A(1:n)  = 1.e0_dp
      endif
      I%IA      = [(j,j=1,n+1)]
      I%JA(1:n) = [(j,j=1,n)]


   end function identity

!          ************************************************************
!          ** Performs the sum of two sparse matrices                **
!          ** C = (alpha * A) + (beta * B)                           **
!          ** alpha, beta real coefficients                          **
!          ************************************************************
   function sparse_sum(alpha,A,beta,B) result(C)
      implicit none

      real (dp)       , intent(in) :: alpha, beta
      type(sparse)    , intent(in) :: A, B
      type(sparse)                 :: C

      integer, dimension(A%nc)     :: itmp
      integer                      :: cur_element, i, k, ka, kb,&
      &jpos

!          Preliminary check on matrix dimensions
      if (A%nr /= B%nr) then
         write(*,*)'Sparse sum: wrong number of rows'
         write(*,*)'A=',A%nr,' B=',B%nr
         stop
      endif

      if (A%nc /= B%nc) then
         write(*,*)'Sparse sum: wrong number of columns'
         write(*,*)'A=',A%nc,' B=',B%nc
         stop
      endif

      if (alpha==0.e0_dp .or. beta==0.e0_dp) then
         write(*,*)'Error: coefficients in sparse sum must be /=0'
         stop
      endif

!          Initialise temporary array
      cur_element = A%n + B%n! - overlaps(A,B)
      call allocate(A%nr, A%nc, cur_element, C)

!          Initialise number of nonzero elements in matrix C
      cur_element  = 0

!          Initialise first element
      C%IA(1) = 1

!          Initialise column working array
      itmp(1:A%nc) = 0

      rows: do i = 1, A%nr

         thisrowA: do ka = A%IA(i), A%IA(i+1) - 1

            if (A%A(ka)==0.e0_dp) cycle thisrowA
            cur_element = cur_element + 1
            C%JA(cur_element) = A%JA(ka)
            C%A (cur_element) = alpha * A%A (ka)
            itmp(A%JA(ka)) = cur_element

         end do thisrowA

         thisrowB: do kb = B%IA(i), B%IA(i+1) - 1

            if (B%A(kb)==0.e0_dp) cycle thisrowB

            jpos = itmp(B%JA(kb))

            newposition: if (jpos == 0) then
               cur_element = cur_element + 1
               C%JA(cur_element) = B%JA(kb)
               C%A(cur_element)  = beta * B%A(kb)
               itmp(B%JA(kb))      = cur_element
            else
               C%A(jpos)         = C%A(jpos)&
               &+ beta * B%A(kb)

               delete_empty: if (C%A(jpos)==0.e0_dp) then
                  cur_element = cur_element - 1
                  C% A(jpos:cur_element) =&
                  &C% A(jpos+1:cur_element+1)
                  C%JA(jpos:cur_element) =&
                  &C%JA(jpos+1:cur_element+1)
                  itmp(B%JA(kb)) = 0
                  itmp(C%JA(jpos:cur_element)) =&
                  &itmp(C%JA(jpos:cur_element)) -1
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
!           call sparse_compress(C)

   end function sparse_sum

!          ************************************************************
!          ** Performs the sum of two sparse matrices                **
!          ** C = (A) + (B)                                          **
!          ** where matrix B has only a partial number of the rows   **
!          ** that matrix A has, and the array iindexB gives the     **
!          ** correspondence indices for the rows of A the rows of   **
!          ** B have to be summed to                                 **
!          ************************************************************
   function sparse_partial_sum(A,B,iB) result(C)
      implicit none

      type(sparse)    , intent(in) :: A, B
      integer         , intent(in) :: iB(:)
      type(sparse)                 :: C


      real (dp), parameter         :: zero = 0.0_dp
      integer, dimension(A%nc)     :: itmp
      integer                      :: cur_element, i, k, ka, kb,&
      &jpos, countrowB,&
      &irowAB
      integer                      :: count_new_pos
      logical                      :: noA, noB
!          Preliminary check on matrix dimensions
      if (size(iB) /= B%nr) then
         write(*,*)'Sparse partial sum: wrong index array'
         write(*,*)'i=',size(iB),' B=',B%nr
         stop
      endif

      if (A%nc /= B%nc) then
         write(*,*)'Sparse sum: wrong number of columns'
         write(*,*)'A=',A%nc,' B=',B%nc
         stop
      endif

!          Initialise temporary array
      call allocate(A%nr, A%nc, A%n+B%n, C)

!          Initialise number of nonzero elements in matrix C
      cur_element   = 0

!          Initialise first element
      C%IA(1) = 1

!          Initialise column working array
      itmp(1:A%nc) = 0

!          Initialise first row of matrix B to be matched with A
      countrowB = 1
      irowAB    = iB(1)

      rows: do i = 1, A%nr

         thisrowA: do ka = A%IA(i), A%IA(i+1) - 1

            if (A%A(ka)==0.e0_dp) cycle thisrowA

            cur_element = cur_element + 1
            C%JA(cur_element) = A%JA(ka)
            C%A (cur_element) = A%A (ka)
            itmp(A%JA(ka)) = cur_element

         end do thisrowA

         matchABrow: if (i == irowAB) then

            thisrowB: do kb = B%IA(countrowB),&
            &B%IA(countrowB+1) - 1

               if (B%A(kb)==0.e0_dp) cycle thisrowB

               jpos = itmp(B%JA(kb))

               newposition: if (jpos == 0) then
                  cur_element = cur_element + 1
                  C%JA(cur_element) = B%JA(kb)
                  C%A(cur_element)  = B%A(kb)
                  itmp(B%JA(kb))      = cur_element
               else
                  C%A(jpos)         = C%A(jpos)&
                  &+ B%A(kb)

                  delete_empty: if (C%A(jpos)==zero) then

                     cur_element = cur_element - 1
                     C% A(jpos:cur_element) =&
                     &C% A(jpos+1:cur_element+1)
                     C%JA(jpos:cur_element) =&
                     &C%JA(jpos+1:cur_element+1)
                     itmp(B%JA(kb)) = 0
                     itmp(C%JA(jpos:cur_element)) =&
                     &itmp(C%JA(jpos:cur_element)) -1
                  endif delete_empty



               endif newposition

            end do thisrowB

!                   Check which is the next line to be summed
            countrowB = countrowB + 1
            if (countrowB<=size(iB)) irowAB = iB(countrowB)

         endif matchABrow

         ! Clean working array, prepare for next row
         clean_itmp: do k = C%IA(i), cur_element
            itmp(C%JA(k)) = 0
         end do clean_itmp

         ! Update number of elements in current row
         C%IA(i+1) = cur_element+1

      end do rows

      C%n = cur_element




   end function sparse_partial_sum

!          ************************************************************
!          ** Performs the internal sum of elements over rows or     **
!          ** columns (equal to fortran's intrinsic sum)             **
!          ************************************************************

   subroutine sparse_internal_sum(A,sm,dim,mask)
      implicit none

      type(sparse),          intent(in)             :: A
      real (dp)       ,      dimension(:), intent(inout) :: sm
      integer,               intent(in),   optional :: dim
      logical, dimension(size(sm)), intent(in),   optional :: mask

      integer :: d, i, j

!          Dimensional checks
      if (present(dim)) then
         d = dim
      else
         d = 1
      endif

      if (.not.(d == 1 .or. d == 2)) then
         write(*,*)'sparse_internal_sum: wrong dim=',d
         stop
      endif

      if (d==1 .and. size(sm) /= A%nc) then
         write(*,*)'sparse_internal_sum: wrong array dim:',size(sm)
         stop
      endif

      if (d==2 .and. size(sm) /= A%nr) then
         write(*,*)'sparse_internal_sum: wrong array dim:',size(sm)
         stop
      endif

      mask_check: if (present(mask)) then
         if ( (d==1 .and. size(mask)/= A%nc)   .or.&
         &(d==2 .and. size(mask)/= A%nr) ) then
            write(*,*)'sparse_internal_sum: wrong mask dimension'
            stop
         endif
      endif mask_check


!          ** dim = 1: sum rows ***************************************
      dimens: if (d == 1) then

         dim_1_sum: if (present(mask)) then

            sm(1:A%nc) = 0.e0_dp
            sparse_els_mask: do i = 1, A%n
               j = A%JA(i)
               if (mask(j)) sm(j) = sm(j) + A%A(i)
            end do sparse_els_mask



         else

            sm(1:A%nc) = 0.e0_dp
            sparse_els: do i = 1, A%n
               j = A%JA(i)
               sm(j) = sm(j) + A%A(i)
            end do sparse_els
         endif dim_1_sum



!          ** dim = 2: sum columns ************************************
      elseif (d == 2) then


         dim_2_sum: if (.not.present(mask)) then
            rowsmask: do i = 1, A%nr
               sm(i) = 0.e0_dp
               row_elems_mask: do j = A%IA(i), A%IA(i+1)-1
                  sm(i) = sm(i) + A%A(j)
               end do row_elems_mask
            end do rowsmask

         else

            rows: do i = 1, A%nr
               sm(i) = 0.e0_dp
               if (mask(i)) then
                  row_elems: do j = A%IA(i), A%IA(i+1)-1
                     sm(i) = sm(i) + A%A(j)
                  end do row_elems
               endif
            end do rows

         endif dim_2_sum

      endif dimens

   end subroutine sparse_internal_sum

   subroutine sparseint_internal_sum(A,sm,dim,mask)
      implicit none

      type(sparseint),       intent(in)             :: A
      integer, dimension(:), intent(inout)          :: sm
      integer,               intent(in),   optional :: dim
      logical, dimension(size(sm)), intent(in),   optional :: mask

      integer :: d, i, j

!          Dimensional checks
      if (present(dim)) then
         d = dim
      else
         d = 1
      endif

      if (.not.(d == 1 .or. d == 2)) then
         write(*,*)'sparse_internal_sum: wrong dim=',d
         stop
      endif

      if (d==1 .and. size(sm) /= A%nc) then
         write(*,*)'sparse_internal_sum: wrong array dim:',size(sm)
         stop
      endif

      if (d==2 .and. size(sm) /= A%nr) then
         write(*,*)'sparse_internal_sum: wrong array dim:',size(sm)
         stop
      endif

      mask_check: if (present(mask)) then
         if ( (d==1 .and. size(mask)/= A%nc)   .or.&
         &(d==2 .and. size(mask)/= A%nr) ) then
            write(*,*)'sparse_internal_sum: wrong mask dimension'
            stop
         endif
      endif mask_check


!          ** dim = 1: sum rows ***************************************
      dimens: if (d == 1) then

         dim_1_sum: if (present(mask)) then

            sm(1:A%nc) = 0
            sparse_els_mask: do i = 1, A%n
               j = A%JA(i)
               if (mask(j)) sm(j) = sm(j) + A%A(i)
            end do sparse_els_mask



         else

            sm(1:A%nc) = 0
            sparse_els: do i = 1, A%n
               j = A%JA(i)
               sm(j) = sm(j) + A%A(i)
            end do sparse_els
         endif dim_1_sum



!          ** dim = 2: sum columns ************************************
      elseif (d == 2) then


         dim_2_sum: if (.not.present(mask)) then
            rowsmask: do i = 1, A%nr
               sm(i) = 0
               row_elems_mask: do j = A%IA(i), A%IA(i+1)-1
                  sm(i) = sm(i) + A%A(j)
               end do row_elems_mask
            end do rowsmask

         else

            rows: do i = 1, A%nr
               sm(i) = 0
               if (mask(i)) then
                  row_elems: do j = A%IA(i), A%IA(i+1)-1
                     sm(i) = sm(i) + A%A(j)
                  end do row_elems
               endif
            end do rows

         endif dim_2_sum

      endif dimens

   end subroutine sparseint_internal_sum

   ! ***********************************************************
   ! ** Counts the number of elements in rows or columns      **
   ! ** in a sparse matrix, eventually given a logical array  **
   ! ** mask                                                  **
   ! **                                                       **
   ! ** Federico Perini, 23/07/2012 (C)                       **
   ! **                                                       **
   ! ***********************************************************

   subroutine sparse_internal_count(A,cnt,dim,mask)
      implicit none

      type(sparse),          intent(in)                  :: A
      integer,               dimension(:), intent(inout) :: cnt
      integer,               intent(in),   optional      :: dim
      logical, dimension(size(cnt)), intent(in), optional:: mask

      integer :: d, i, j

!          Dimensional checks
      if (present(dim)) then
         d = dim
      else
         d = 1
      endif

      if (.not.(d == 1 .or. d == 2)) then
         write(*,*)'sparse_internal_count: wrong dim=',d
         stop
      endif

      if (d==1 .and. size(cnt) /= A%nc) then
         write(*,*)'sparse_internal_count: wrong array dim:',&
         &size(cnt),', A%nc=',A%nc
         stop
      endif

      if (d==2 .and. size(cnt) /= A%nr) then
         write(*,*)'sparse_internal_count: wrong array dim:',&
         &size(cnt),', A%nr=',A%nr
         stop
      endif

      mask_check: if (present(mask)) then
         if ( (d==1 .and. size(mask)/= A%nc)   .or.&
         &(d==2 .and. size(mask)/= A%nr) ) then
            write(*,*)'sparse_internal_count: wrong mask dimension'
            stop
         endif
      endif mask_check


!          ** dim = 1: count rows **************************************
      dimens: if (d == 1) then

         dim_1_sum: if (present(mask)) then


            cnt(1:A%nc) = 0
            sparse_els_mask: do i = 1, A%n
               j = A%JA(i)
               if (mask(j)) cnt(j) = cnt(j) + 1
            end do sparse_els_mask

         else

            cnt(1:A%nc) = 0
            sparse_els: do i = 1, A%n
               j = A%JA(i)
               cnt(j) = cnt(j) + 1
            end do sparse_els
         endif dim_1_sum

!          ** dim = 2: count over columns *****************************
      elseif (d == 2) then

         dim_2_sum: if (.not.present(mask)) then
            rowsmask: do i = 1, A%nr
               cnt(i) = A%IA(i+1) - A%IA(i)
            end do rowsmask

         else

            rows: do i = 1, A%nr
               cnt(i) = 0
               if (mask(i)) cnt(i) = A%IA(i+1) - A%IA(i)
            end do rows

         endif dim_2_sum

      endif dimens

   end subroutine sparse_internal_count


!          ************************************************************
!          ** Logical operator for ">" (.gt.) sign. Output is a      **
!          ** logical mask whose size is size(spmatrix%A)            **
!          ************************************************************
   function sparse_gt_array(spmatrix,val) result(gtmask)

      type(sparse),             intent(in) :: spmatrix
      real (dp)       , dimension(size(spmatrix%A)),&
      &intent(in) :: val
      logical, dimension(size(spmatrix%A)) :: gtmask

      ! Do the logical ">" operator
      where (spmatrix%A > val)
         gtmask = .true.
      elsewhere
         gtmask = .false.
      end where

   end function sparse_gt_array


   function sparse_gt_dble(spmatrix,val) result(gtmask)

      type(sparse),             intent(in) :: spmatrix
      real (dp)       ,         intent(in) :: val
      logical, dimension(size(spmatrix%A)) :: gtmask

      ! Do the logical ">" operator
      where (spmatrix%A > val)
         gtmask = .true.
      elsewhere
         gtmask = .false.
      end where

   end function sparse_gt_dble



!          ************************************************************
!          ** Sort sparse matrix columns in ascending order          **
!          ************************************************************

   subroutine sparse_sort(A)
      implicit none

      type(sparse), intent(inout), target         :: A
      integer            :: i, ji, jo, r0, rf
      integer            :: tmpJA
      real (dp)          :: tmpA

      if (.not.(A%nr>0.and.A%nc>0.and.A%n>0)) then
         write(*,*)'sparse_compress: wrong matrix dim:'
         write(*,*)A%nr, A%nc, A%n
         stop
      endif



      row: do i = 1, A%nr

         r0 = A%IA(i)
         rf = A%IA(i+1)-1

         outer_cols: do jo = r0, rf

            loop_inner: do ji = r0+1, rf

               sort: if (A%JA(ji)>A%JA(jo)) then

                  tmpJA = A%JA(ji)
                  tmpA  = A%A (ji)
                  A%JA(ji) = A%JA(jo)
                  A%A (ji) = A%A (jo)
                  A%JA(jo) = tmpJA
                  A%A (jo) = tmpA

               endif sort

            end do loop_inner

         end do outer_cols

      end do row

   end subroutine sparse_sort

!          ************************************************************
!          ** Performs the internal product of elements over rows or **
!          ** columns (equal to fortran's intrinsic product)         **
!          ************************************************************

   subroutine sparse_internal_prod(A,sm,dim,mask)
      implicit none

      type(sparse),          intent(in)             :: A
      real (dp)       ,      intent(inout),dimension(:) :: sm
      integer,               intent(in),   optional :: dim
      logical, dimension(:), intent(in),   optional :: mask

      integer :: d, i, j

!          Dimensional checks
      if (present(dim)) then
         d = dim
      else
         d = 1
      endif

      if (.not.(d == 1 .or. d == 2)) then
         write(*,*)'sparse_internal_sum: wrong dim=',d
         stop
      endif

      if (d==1 .and. size(sm) /= A%nc) then
         write(*,*)'sparse_internal_sum: wrong array dim:',size(sm)
         stop
      endif

      if (d==2 .and. size(sm) /= A%nr) then
         write(*,*)'sparse_internal_sum: wrong array dim:',size(sm)
         stop
      endif

      mask_check: if (present(mask)) then
         if ( (d==1 .and. size(mask)/= A%nc)   .or.&
         &(d==2 .and. size(mask)/= A%nr) ) then
            write(*,*)'sparse_internal_sum: wrong mask dimension'
            stop
         endif
      endif mask_check


!          ** dim = 1: sum rows ***************************************
      dimens: if (d == 1) then

         dim_1_sum: if (present(mask)) then

            sm(1:A%nc) = 1.e0_dp
            sparse_els_mask: do i = 1, A%n
               j = A%JA(i)
               if (mask(j)) sm(j) = sm(j) * A%A(i)
            end do sparse_els_mask



         else

            sm(1:A%nc) = 1.e0_dp
            sparse_els: do i = 1, A%n
               j = A%JA(i)
               sm(j) = sm(j) * A%A(i)
            end do sparse_els
         endif dim_1_sum



!          ** dim = 2: sum columns ************************************
      elseif (d == 2) then


         dim_2_sum: if (.not.present(mask)) then
            rowsmask: do i = 1, A%nr
               sm(i) = 1.e0_dp
               row_elems_mask: do j = A%IA(i), A%IA(i+1)-1
                  sm(i) = sm(i) * A%A(j)
               end do row_elems_mask
            end do rowsmask

         else

            rows: do i = 1, A%nr
               sm(i) = 1.e0_dp
               if (mask(i)) then
                  row_elems: do j = A%IA(i), A%IA(i+1)-1
                     sm(i) = sm(i) * A%A(j)
                  end do row_elems
               endif
            end do rows

         endif dim_2_sum

      endif dimens

   end subroutine sparse_internal_prod

!          ************************************************************
!          ** Performs the internal product of elements over rows or **
!          ** columns (equal to fortran's intrinsic product)         **
!          ************************************************************

   function sparse_internal_prod1(A) result(pd)
      implicit none

      type(sparse),          intent(in)             :: A
      real (dp)       ,      dimension(A%nc)        :: pd


      integer :: i, j


      if (size(pd) /= A%nc) then
         write(*,*)'sparse_int_prod: wrong array dim:',size(pd)
         stop
      endif


!          ** dim = 1: sum rows ***************************************
      pd(1:A%nc) = 1.e0_dp
      sparse_els: do i = 1, A%n
         j = A%JA(i)
         pd(j) = pd(j) * A%A(i)
      end do sparse_els

   end function sparse_internal_prod1

!          ************************************************************
!          ** Multiplies an entire sparse matrix by a real           **
!          ** B = alpha * A                                          **
!          ************************************************************
   function sparse_real_prod(alpha,A) result(B)
      implicit none

      real (dp)       , intent(in) :: alpha
      type(sparse)    , intent(in) :: A
      type(sparse)                 :: B


!          Initialise B matrix
      call allocate(A%nr, A%nc, A%n, B)

!          Copying position indices
      B%IA = A%IA
      B%JA = A%JA

!          Performing multiplication
      B%A = alpha * A%A

   end function sparse_real_prod

!          ************************************************************
!          ** Multiplies each row of the sparse matrix by real values**
!          ** defined in the array b (1 x A%ncols)                   **
!          ************************************************************
   function sparse_row_prod(A,b) result(Amult)
      implicit none

      type(sparse)    ,               intent(in) :: A
      real (dp)       , dimension(:), intent(in) :: b
      type(sparse)                               :: Amult

      integer                                    :: i, j, jj

!          Check on the dimensions of b
      if (size(b) /= A%nc) then
         write(*,*)'Error in sparse_row_prod'
         write(*,*)'A%ncols = ',A%nc
         write(*,*)'array b = ',size(b)
         stop
      endif

!          Initialise Amult matrix
      call allocate(A%nr, A%nc, A%n, Amult)

!          Copying position indices
      Amult%IA = A%IA
      Amult%JA = A%JA

!          Performing multiplication
      rows: do i = 1, A%nr
         cols: do j = A%IA(i), A%IA(i+1)-1
            jj         = A%JA(j)
            Amult%A(j) = A%A (j) * b(jj)
         end do cols
      end do rows

   end function sparse_row_prod

!          ************************************************************
!          ** Multiplies each row of the sparse matrix by real values**
!          ** defined in the array b (1 x A%ncols)                   **
!          ************************************************************
   subroutine sparse_row_prod_valonly(A,b)
      implicit none

      type(sparse)    ,               intent(inout) :: A
      real (dp)       , dimension(:), intent(in)    :: b

      integer                                    :: i, j, jj

!          Check on the dimensions of b
      if (size(b) /= A%nc) then
         write(*,*)'Error in sparse_row_prod'
         write(*,*)'A%ncols = ',A%nc
         write(*,*)'array b = ',size(b)
         stop
      endif

!          Performing multiplication
      rows: do i = 1, A%nr
         cols: do j = A%IA(i), A%IA(i+1)-1
            jj     = A%JA(j)
            A%A(j) = A%A (j) * b(jj)
         end do cols
      end do rows

   end subroutine sparse_row_prod_valonly

!          ************************************************************
!          ** Computes power of an array to a sparse matrix, where   **
!          ** the array size is equal to the number of columns in the**
!          ** matrix (i.e. the array is superposed to every line)    **
!          ** and then its inner product along the columns is comp.  **
!          ** prod(j) = product(b(i)^A(j,i))
!          ************************************************************
   function row_power_to_sparse_product(b,A) result(prod)
      implicit none

      type(sparseint) ,                  intent(in)    :: A
      real (dp)       , dimension(A%nc), intent(in)    :: b
      real (dp)       , dimension(A%nr)                :: prod
      real (dp)                                        :: rval

      integer                                    :: i, j, jj

      empty_matrix: if (A%n == 0) then
         prod = 0.e0_dp
         return
      endif empty_matrix

!          Performing multiplication
      rows: do i = 1, A%nr
         if (A%IA(i+1)==A%IA(i)) then
            prod(i) = 0.e0_dp
         else

            prod(i) = 1.e0_dp
            cols: do j = A%IA(i), A%IA(i+1)-1
               jj      = A%JA(j)
               prod(i) = prod(i) * b(jj)
               if (A%A(j)> 1) then
                  prod(i)=prod(i)*b(jj)
                  if (A%A(j)> 2)&
                  &prod(i)=prod(i)*b(jj)**(A%A(j)-2)
               elseif (A%A(j)<0) then
                  rval = 1.e0_dp/b(jj)
                  prod(i) = prod(i) * rval**abs(A%A(j)-1)
               endif
            end do cols
         endif
      end do rows

   end function row_power_to_sparse_product

!          ************************************************************
!          ** Multiplies each column of the sparse matrix by real    **
!          ** values defined in the array b (1 x A%nrows)            **
!          ************************************************************
   function sparse_col_prod(A,b) result(Amult)
      implicit none

      type(sparse)    ,               intent(in) :: A
      real (dp)       , dimension(:), intent(in) :: b
      type(sparse)                               :: Amult

      integer                                    :: i, j

!          Check on the dimensions of b
      if (size(b) /= A%nr) then
         write(*,*)'Error in sparse_col_prod'
         write(*,*)'A%nr    = ',A%nr
         write(*,*)'array b = ',size(b)
         stop
      endif


!          Initialise Amult matrix
      call allocate(A%nr, A%nc, A%n, Amult)

!          Copying position indices
      Amult%IA = A%IA
      Amult%JA(1:A%n) = A%JA(1:A%n)

!          Performing multiplication
      rows: do i = 1, A%nr
         cols: do j = A%IA(i), A%IA(i+1)-1
            Amult%A(j) = A%A(j) * b(i)
         end do cols
      end do rows

   end function sparse_col_prod

!          ************************************************************
!          ** Multiplies each column of the sparse matrix by real    **
!          ** values defined in the array b (1 x A%nrows)            **
!          ************************************************************
   subroutine sparse_col_prod_valonly(A,b)
      implicit none

      type(sparse)    ,               intent(inout) :: A
      real (dp)       , dimension(:), intent(in)    :: b
      integer                                    :: i, j

!          Check on the dimensions of b
      if (size(b) /= A%nr) then
         write(*,*)'Error in sparse_col_prod'
         write(*,*)'A%nr    = ',A%nr
         write(*,*)'array b = ',size(b)
         stop
      endif

!          Performing multiplication
      rows: do i = 1, A%nr
         cols: do j = A%IA(i), A%IA(i+1)-1
            A%A(j) = A%A(j) * b(i)
         end do cols
      end do rows

   end subroutine sparse_col_prod_valonly

!          ************************************************************
!          ** Builds a new sparse matrix by blocking horizontally    **
!          ** two sparse matrices                                    **
!          **   ----------         ----- -----                       **
!          **   |   C    |    =    | A | | B |                       **
!          **   ----------         ----- -----                       **
!          ************************************************************
   function sparse_block_horizontal(A,B) result(C)
      implicit none

      type(sparse)    , intent(in) :: A, B
      type(sparse)                 :: C

      integer :: i, ia, ib, ic

!          Check on the number of rows
      if (A%nr /= B%nr) then
         write(*,*)'Cannot block with different number of rows'
         write(*,*)'A%nr = ',A%nr,' B%nr = ',B%nr
         stop
      endif

!          Initialize matrix C
      call allocate(A%nr,A%nc+B%nc,A%n+B%n,C)

!          Assign matrix
      ic = 0
      C%IA(1) = 1
      rows: do i = 1, A%nr

         colsA: do ia = A%IA(i), A%IA(i+1)-1

            ic = ic + 1
            C%A (ic) = A%A (ia)
            C%JA(ic) = A%JA(ia)

         end do colsA

         colsB: do ib = B%IA(i), B%IA(i+1)-1

            ic = ic + 1
            C%A (ic) = B%A (ib)
            C%JA(ic) = B%JA(ib) + A%nc

         end do colsB

         C%IA(i+1) = ic+1

      end do rows


   end function sparse_block_horizontal

!          ************************************************************
!          ** Builds a new sparse matrix by blocking vertically      **
!          ** two sparse matrices                                    **
!          **   -----         -----                                  **
!          **   |   |         | A |                                  **
!          **   | C |    =    -----                                  **
!          **   |   |         -----                                  **
!          **   |   |         | B |                                  **
!          **   -----         -----                                  **
!          ************************************************************
   function sparse_block_vertical(A,B) result(C)
      implicit none

      type(sparse)    , intent(in) :: A, B
      type(sparse)                 :: C
      integer :: j


!          Check on the number of rows
      if (A%nc /= B%nc) then
         write(*,*)'Cannot block with different number of columns'
         write(*,*)'A%nc = ',A%nc,' B%nc = ',B%nc
         stop
      endif

!          Initialize matrix C
      call allocate(A%nr+B%nr,A%nc,A%n+B%n,C)



!          The element and column indexes do not change; they only
!          have to be merged between A and B
      C%A (1:A%n)     = A%A(1:A%n)
      C%JA(1:A%n)     = A%JA(1:A%n)

      C%A (A%n+1:C%n) = B%A(1:B%n)
      C%JA(A%n+1:C%n) = B%JA(1:B%n)

!          Row indexes of the first part are the same as the A matrix;
!          those of the B part just need be added the total amount of
!          elements of A (check empty lines!!)
      C%IA(1:A%nr)        = A%IA(1:A%nr)
      C%IA(A%nr+1:C%nr+1) = B%IA(1:B%nr+1) + A%IA(A%nr+1) - 1

!
!           do j = A%n+1, C%n
!             write(*,*)C%A(j),B%A(j-A%n),C%JA(j),B%JA(j-A%n)
!           end do

!           call matrix_details(A,'A')
!           call matrix_details(B,'B')
!           call matrix_details(C,'C')

   end function sparse_block_vertical

!          ************************************************************
!          ** Builds a new sparse matrix by blocking diagonally      **
!          ** two sparse matrices. Input matrix dimensions are       **
!          ** arbitrary                                              **
!          **   -----------         ------------                     **
!          **   |         |         | A |      |                     **
!          **   |    C    |    =    ------------                     **
!          **   |         |         |   |   B  |                     **
!          **   |         |         |   |      |                     **
!          **   -----------         ------------                     **
!          ************************************************************
   function sparse_block_diagonal(A,B) result(C)
      implicit none

      type(sparse)    , intent(in) :: A, B
      type(sparse)                 :: C

!          Initialize matrix C
      call allocate(A%nr+B%nr,A%nc+B%nc,A%n+B%n,C)

!          The element and column indexes do not change; they only
!          have to be merged between A and B
      C%A (1:A%n)     = A%A(1:A%n)
      C%JA(1:A%n)     = A%JA(1:A%n)

      C%A (A%n+1:C%n) = B%A(1:B%n)
      C%JA(A%n+1:C%n) = B%JA(1:B%n) + A%nc

!          Row indexes of the first part are the same as the A matrix;
!          those of the B part just need be added the total amount of
!          elements of A (check empty lines!!)
      C%IA(1:A%nr)        = A%IA(1:A%nr)
      C%IA(A%nr+1:C%nr+1) = B%IA(1:B%nr+1) + A%IA(A%nr+1) - 1

   end function sparse_block_diagonal

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


end module sparse_algebra

