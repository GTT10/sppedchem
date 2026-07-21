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
c     **                                                             **
c     **    Module for MPI communications of sparse matrix algebra   **
c     **                                                             **
c     **                                                             **
c     **   Author:      Federico Perini                              **
c     **   Last update: friday, 01/06/2012                           **
c     **                                                             **
c     *****************************************************************
      module sparse_MPI

         use working_precision, only: dp

         implicit none
         private

c        ** Procedures accessible from the outside
         public :: mpi_broadcast
         public :: mpi_initialise
         public :: mpi_pointer_broadcast
         public :: current_mpi_cpuid
         public :: total_mpi_cpus

c        ** The "verbose" routine sets output options for the
c        ** broadcasting routines
         logical, parameter :: verbose = .false.


         interface mpi_broadcast
            module procedure sparse_broadcast_all
            module procedure sparseint_broadcast_all
            module procedure sparse_ordered_broadcast
            module procedure double_matrix_broadcast
            module procedure integer_matrix_broadcast
            module procedure double_array_broadcast
            module procedure integer_array_broadcast
            module procedure integer_broadcast
            module procedure arbitrary_broadcast
            module procedure logical_broadcast
            module procedure character_broadcast
            module procedure logical_array_broadcast
            module procedure character_array_broadcast
            module procedure fixed_integer_array_broadcast
            module procedure fixed_double_array_broadcast
            module procedure fixed_character_array_broadcast
         end interface mpi_broadcast

         interface mpi_pointer_broadcast
!            module procedure double_pointer_broadcast
            module procedure arbitrary_pointer_broadcast
         end interface mpi_pointer_broadcast


         contains

c        **************************************************************
c        **                                                          **
c        **  Subroutine that initializes the MPI space and returns   **
c        **  the number of processor and the current processor ID    **
c        **                                                          **
c        **************************************************************
         subroutine mpi_initialise(tot_cpu_number, cur_cpu_id)
            use mpi
         implicit none

         integer, intent(out) :: tot_cpu_number, cur_cpu_id

         integer :: status(MPI_STATUS_SIZE)
         integer :: i_error


         call MPI_Init(i_error)

         if (i_error/=0)
     &    write(*,*)' Warning: MPI init exited with error ',i_error

         call MPI_Comm_size(MPI_COMM_WORLD, tot_cpu_number, i_error)
         call MPI_Comm_rank(MPI_COMM_WORLD, cur_cpu_id,     i_error)

         end subroutine mpi_initialise





c        **************************************************************
c        **                                                          **
c        **  Function that computes the flag of the MPI real number  **
c        **  datatype corresponding to the current machine precision **
c        **                                                          **
c        **************************************************************

         function real_precision_flag(x) result(MPIflag)
            use mpi
         implicit none

         real (dp), intent(in) :: x
         integer               :: MPIflag

         integer :: status(MPI_STATUS_SIZE), xtype

         character(len=*), parameter ::
     &     fmt_er = "(1x,'MPI does not support reals of type',1x,I2)"


c        Gather Fortran type
         xtype = kind(x)

         if (xtype == 4) then
           MPIflag = MPI_REAL4
         elseif (xtype == 8) then
           MPIflag = MPI_REAL8
         elseif (xtype == 16) then
           MPIflag = MPI_REAL16
         else
           MPIflag = -1
           write(*,fmt_er)xtype
           stop
         endif

         end function real_precision_flag

c        **************************************************************
c        **                                                          **
c        **  Broadcast sparse matrix from node 0 to all the other    **
c        **  nodes                                                   **
c        **                                                          **
c        **************************************************************
         subroutine MSR_broadcast_all(matrix)
            use sparse_algebra
            use mpi

         implicit none

         integer :: status(MPI_STATUS_SIZE)
         type(sparse), intent(inout) :: matrix
         integer  :: root, node, nrows, ncols, nelems, ierror, nnodes
         integer  :: sparse_size(3), j, rp, tag, cw, inode, mi
         external :: exita

         root = 0
         mi   = MPI_INTEGER
         rp   = real_precision_flag(0.0_dp)
         cw   = MPI_COMM_WORLD

c        Get current processor id
         call MPI_COMM_RANK(MPI_COMM_WORLD, node, ierror)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, nnodes, ierror)

c        Matrix doesn't need to be initialised
         masternode: if (node == 0) then

            if (.not.matrix%n>0) then


               if (verbose) write(*,*)'Warning: passed unallocated '//
     &                                'matrix in sparse_broadcast'
               sparse_size = [0,0,0]

            else

            nrows  = matrix%nr
            ncols  = matrix%nc
            nelems = matrix%n

            sparse_size = [nrows, ncols, nelems]

            endif

         endif masternode

c        Send sparse matrix dimensions to other nodes
         call MPI_BCAST(sparse_size,3,mi,root,cw,ierror)

         if (sum(sparse_size)>0) then

         if (node /= 0) then

c        Other nodes receive matrix dimensions to initialise it
            nrows = sparse_size(1)
            ncols = sparse_size(2)
            nelems = sparse_size(3)

            call sparse_allocate_det(nrows,ncols,nelems,matrix)

         endif

c        Now that all the matrices are initialised, broadcast data
         call MPI_BCAST(matrix%A, nelems, rp,root,cw,ierror)
         call MPI_BCAST(matrix%IA,nrows+1,mi,root,cw,ierror)
         call MPI_BCAST(matrix%JA,nelems, mi,root,cw,ierror)

         else
            if (node /= 0) then
              call sparse_nullify_general(matrix)
            endif
         endif

         end subroutine MSR_broadcast_all


c        **************************************************************
c        **                                                          **
c        **  Broadcast sparse matrix from node 0 to all the other    **
c        **  nodes                                                   **
c        **                                                          **
c        **************************************************************
         subroutine sparse_broadcast_all(matrix)
         use sparse_algebra
         use mpi

         implicit none

         integer :: status(MPI_STATUS_SIZE)
         type(sparse), intent(inout) :: matrix
         integer  :: root, node, nrows, ncols, nelems, ierror, nnodes
         integer  :: sparse_size(3), j, rp, tag, cw, inode, mi
         external :: exita

         root = 0
         mi   = MPI_INTEGER
         rp   = real_precision_flag(0.0_dp)
         cw   = MPI_COMM_WORLD

c        Get current processor id
         call MPI_COMM_RANK(MPI_COMM_WORLD, node, ierror)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, nnodes, ierror)

c        Matrix doesn't need to be initialised
         masternode: if (node == 0) then

            if (.not.matrix%n>0) then


               if (verbose) write(*,*)'Warning: passed unallocated '//
     &                                'matrix in sparse_broadcast'
               sparse_size = [0,0,0]

            else

            nrows  = matrix%nr
            ncols  = matrix%nc
            nelems = matrix%n

            sparse_size = [nrows, ncols, nelems]

            endif

         endif masternode

c        Send sparse matrix dimensions to other nodes
         call MPI_BCAST(sparse_size,3,mi,root,cw,ierror)

         if (sum(sparse_size)>0) then

         if (node /= 0) then

c        Other nodes receive matrix dimensions to initialise it
            nrows = sparse_size(1)
            ncols = sparse_size(2)
            nelems = sparse_size(3)

            call sparse_allocate_det(nrows,ncols,nelems,matrix)

         endif

c        Now that all the matrices are initialised, broadcast data
         call MPI_BCAST(matrix%A, nelems, rp,root,cw,ierror)
         call MPI_BCAST(matrix%IA,nrows+1,mi,root,cw,ierror)
         call MPI_BCAST(matrix%JA,nelems, mi,root,cw,ierror)

         else
            if (node /= 0) then
              call sparse_nullify_general(matrix)
            endif
         endif

         end subroutine sparse_broadcast_all


c        **************************************************************
c        **                                                          **
c        **  Broadcast sparse_ordered  matrix from node 0 to all     **
c        **  the other nodes, including (if present) information     **
c        **  regarding ordering and factorization                    **
c        **                                                          **
c        **************************************************************
         subroutine sparse_ordered_broadcast(matrix)
         use sparse_algebra
         use sparse_definitions
         use mpi

         implicit none

         integer :: status(MPI_STATUS_SIZE)
         type(sparse_ordered), intent(inout) :: matrix
         integer :: root, node, nrows, ncols, nelems, ierror, nnodes
         integer :: sparse_size(3), j, rp, tag, cw, inode, mi, ml
         logical :: ordering_flags(3)
         integer :: factorization_dims(2)
         external :: exita

         root = 0
         mi   = MPI_INTEGER
         ml   = MPI_LOGICAL
         rp   = real_precision_flag(0.0_dp)
         cw   = MPI_COMM_WORLD

c        Get current processor id
         call MPI_COMM_RANK(MPI_COMM_WORLD, node, ierror)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, nnodes, ierror)

c        Matrix doesn't need to be initialised
         masternode: if (node == 0) then

            if (.not.matrix%n>0) then


               if (verbose) write(*,*)'Warning: passed unallocated '//
     &                                'matrix in sparse_broadcast'

               sparse_size    = [0,0,0]
               ordering_flags = [.false.,.false.,.false.]
               factorization_dims = [0,0]

            else

            nrows  = matrix%nr
            ncols  = matrix%nc
            nelems = matrix%n

            sparse_size = [nrows, ncols, nelems]

            ordering_flags(1) = matrix%is_ordered
            ordering_flags(2) = matrix%symbolically_factorized
            ordering_flags(3) = matrix%numerically_factorized

            factorization_dims = [matrix%lR, matrix%lI]

            endif

         endif masternode

c        Send sparse matrix dimensions to other nodes
         call MPI_BCAST(sparse_size       ,3,mi,root,cw,ierror)
         call MPI_BCAST(ordering_flags    ,3,ml,root,cw,ierror)
         call MPI_BCAST(factorization_dims,2,mi,root,cw,ierror)

         if (sum(sparse_size)>0) then

         if (node /= 0) then

c        Other nodes receive matrix dimensions to initialise it
            nrows = sparse_size(1)
            ncols = sparse_size(2)
            nelems = sparse_size(3)

            call allocate(nrows,ncols,nelems,matrix)

            matrix%is_ordered              = ordering_flags(1)
            matrix%symbolically_factorized = ordering_flags(2)
            matrix%numerically_factorized  = ordering_flags(3)

            matrix%lR = factorization_dims(1)
            matrix%lI = factorization_dims(2)

            if (matrix%is_ordered) then
                allocate(matrix%perm(nrows),matrix%inv_perm(nrows))
            endif

         if (matrix%lR>0) allocate(matrix%real_space(matrix%lR))
         if (matrix%lI>0) allocate(matrix%int_space (matrix%lI))

         endif

c        Now that all the matrices are initialised, broadcast data
         call MPI_BCAST(matrix%A, nelems, rp,root,cw,ierror)
         call MPI_BCAST(matrix%IA,nrows+1,mi,root,cw,ierror)
         call MPI_BCAST(matrix%JA,nelems, mi,root,cw,ierror)

         ordering: if (matrix%is_ordered) then
           call MPI_BCAST(matrix%perm    ,nrows, mi,root,cw,ierror)
           call MPI_BCAST(matrix%inv_perm,nrows, mi,root,cw,ierror)
         endif ordering

         if (matrix%lR>0)
     &    call MPI_BCAST(matrix%real_space,matrix%lR, rp,root,cw,ierror)
         if (matrix%lI>0)
     &    call MPI_BCAST(matrix% int_space,matrix%lI, mi,root,cw,ierror)

         else

            if (node /= 0) call deallocate(matrix)

         endif

         end subroutine sparse_ordered_broadcast


c        **************************************************************
c        **                                                          **
c        **  Broadcast sparse matrix from node 0 to all the other    **
c        **  nodes (integer matrix type)                             **
c        **                                                          **
c        **************************************************************
         subroutine sparseint_broadcast_all(matrix)
         use sparse_algebra
         use mpi

         implicit none

         integer :: status(MPI_STATUS_SIZE)
         type(sparseint), intent(inout) :: matrix
         integer :: root, node, nrows, ncols, nelems, ierror, nnodes
         integer :: sparse_size(3), j, rp, tag, cw, inode, mi
         external :: exita

         root = 0
         mi   = MPI_INTEGER
         rp   = real_precision_flag(0.0_dp)
         cw   = MPI_COMM_WORLD

c        Get current processor id
         call MPI_COMM_RANK(MPI_COMM_WORLD, node, ierror)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, nnodes, ierror)

c        Matrix doesn't need to be initialised
         masternode: if (node == 0) then

            if (.not.matrix%n>0) then


               if (verbose) write(*,*)'Warning: passed unallocated '//
     &                                'matrix in sparse_broadcast'
               sparse_size = [0,0,0]

            else

            nrows  = matrix%nr
            ncols  = matrix%nc
            nelems = matrix%n

            sparse_size = [nrows, ncols, nelems]

            endif

         endif masternode

c        Send sparse matrix dimensions to other nodes
         call MPI_BCAST(sparse_size,3,mi,root,cw,ierror)

         if (sum(sparse_size)>0) then

         if (node /= 0) then

c        Other nodes receive matrix dimensions to initialise it
            nrows = sparse_size(1)
            ncols = sparse_size(2)
            nelems = sparse_size(3)

            call sparseint_allocate_det(nrows,ncols,nelems,matrix)

         endif

c        Now that all the matrices are initialised, broadcast data
         call MPI_BCAST(matrix%A, nelems, mi,root,cw,ierror)
         call MPI_BCAST(matrix%IA,nrows+1,mi,root,cw,ierror)
         call MPI_BCAST(matrix%JA,nelems, mi,root,cw,ierror)

         else
            if (node /= 0) then
              call sparse_nullify_int(matrix)
            endif
         endif

         end subroutine sparseint_broadcast_all


c        **************************************************************
c        **                                                          **
c        **  Broadcast double prec matrix from node 0 to all the     **
c        **  other nodes                                             **
c        **                                                          **
c        **************************************************************
         subroutine double_matrix_broadcast(matrix)
            use mpi

         implicit none

         integer :: status(MPI_STATUS_SIZE)
         real (dp)       , dimension(:,:), allocatable,
     &                                     intent(inout) :: matrix
         integer :: root, node, nrows, ncols, nelems, ierror, nnodes
         integer :: matrix_size(4), j, rp, tag, cw, inode, mi
         integer :: l1, u1, l2, u2
         external :: exita

         root = 0
         mi   = MPI_INTEGER
         rp   = real_precision_flag(0.0_dp)
         cw   = MPI_COMM_WORLD

c        Get current processor id
         call MPI_COMM_RANK(MPI_COMM_WORLD, node, ierror)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, nnodes, ierror)

c        Matrix doesn't need to be initialised
         masternode: if (node == 0) then

            if (.not.allocated(matrix)) then

               if (verbose) write(*,*) ' Warning: passed unallocated '//
     &                             'matrix in double_matrix_broadcast'

               matrix_size = [0,0,0,0]

            else

               nrows  = size(matrix,1)
               ncols  = size(matrix,2)
               nelems = nrows*ncols
               l1     = lbound(matrix,1)
               u1     = ubound(matrix,1)
               l2     = lbound(matrix,2)
               u2     = ubound(matrix,2)

               matrix_size = [l1, u1, l2, u2]

            endif

         endif masternode

c        Send sparse matrix dimensions to other nodes
         call MPI_BCAST(matrix_size,4,mi,root,cw,ierror)

         if (sum(matrix_size)>0) then

            if (node /= 0) then

c           Other nodes receive matrix dimensions to initialise it
               l1 = matrix_size(1)
               u1 = matrix_size(2)
               l2 = matrix_size(3)
               u2 = matrix_size(4)

               nrows = u1 - l1 + 1

               if (allocated(matrix)) deallocate(matrix)
               allocate(matrix(l1:u1,l2:u2))

            endif

c           Now that all the matrices are initialised, broadcast data
            do j = l2, u2
               call MPI_BCAST(matrix(l1:u1,j),nrows,rp,root,cw,ierror)
           end do

         else

            if (allocated(matrix)) deallocate(matrix)

         endif

         end subroutine double_matrix_broadcast


c        **************************************************************
c        **                                                          **
c        **  Broadcast integer matrix from node 0 to all the         **
c        **  other nodes                                             **
c        **                                                          **
c        **************************************************************
         subroutine integer_matrix_broadcast(matrix)

         use mpi

         implicit none

         integer :: status(MPI_STATUS_SIZE)
         integer, dimension(:,:), allocatable, intent(inout) :: matrix
         integer :: root, node, nrows, ncols, nelems, ierror, nnodes
         integer :: matrix_size(4), j, rp, tag, cw, inode, mi
         integer :: l1, u1, l2, u2

         root = 0
         mi   = MPI_INTEGER
         rp   = real_precision_flag(0.0_dp)
         cw   = MPI_COMM_WORLD

c        Get current processor id
         call MPI_COMM_RANK(MPI_COMM_WORLD, node, ierror)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, nnodes, ierror)

c        Matrix doesn't need to be initialised
         masternode: if (node == 0) then

            if (.not.allocated(matrix)) then

               if (verbose) write(*,*)'Warning: broadcast unallocated'//
     &                           ' matrix in integer_matrix_broadcast'

               matrix_size = [0,0,0,0]

            else

               nrows  = size(matrix,1)
               ncols  = size(matrix,2)
               nelems = nrows*ncols

               l1     = lbound(matrix,1)
               u1     = ubound(matrix,1)
               l2     = lbound(matrix,2)
               u2     = ubound(matrix,2)

               matrix_size = [l1, u1, l2, u2]

            endif

         endif masternode

c        Send sparse matrix dimensions to other nodes
         call MPI_BCAST(matrix_size,4,mi,root,cw,ierror)

         if (sum(matrix_size)>0) then

            if (node /= 0) then

c              Other nodes receive matrix dimensions to initialise it
               l1 = matrix_size(1)
               u1 = matrix_size(2)
               l2 = matrix_size(3)
               u2 = matrix_size(4)

               nrows = u1 - l1 + 1

               if (allocated(matrix)) deallocate(matrix)
               allocate(matrix(l1:u1,l2:u2))

            endif

c           Now that all the matrices are initialised, broadcast data
            do j = l2, u2
               call MPI_BCAST(matrix(l1:u1,j),nrows,mi,root,cw,ierror)
            end do

         else

            if (allocated(matrix)) deallocate(matrix)

         endif

         end subroutine integer_matrix_broadcast


c        **************************************************************
c        **                                                          **
c        **  Broadcast double prec array from node 0 to all the      **
c        **  other nodes                                             **
c        **                                                          **
c        **************************************************************
         subroutine double_array_broadcast(array)

         use mpi

         implicit none

         integer :: status(MPI_STATUS_SIZE)
         real (dp)       , dimension(:), allocatable,
     &                                     intent(inout) :: array
         integer :: root, node, ierror, nnodes
         integer :: array_size(2), j, rp, tag, cw, inode, mi, l1, u1
         external :: exita

         root = 0
         mi   = MPI_INTEGER
         rp   = real_precision_flag(0.0_dp)
         cw   = MPI_COMM_WORLD

c        Get current processor id
         call MPI_COMM_RANK(MPI_COMM_WORLD, node, ierror)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, nnodes, ierror)

c        Matrix doesn't need to be initialised
         masternode: if (node == 0) then

            if (.not.allocated(array)) then
               if (verbose) write(*,*)'Warning: broadcast unallocated '
     &                             // 'array in double_array_broadcast'

               array_size = [0,0]

            else

               array_size(1)  = lbound(array,1)
               array_size(2)  = ubound(array,1)

            endif

         endif masternode

c        Broadcast array dimensions to other nodes
         call MPI_BCAST(array_size,2,mi,root,cw,ierror)

         if (sum(array_size)>0) then

            if (node /= 0) then

c           Other nodes receive matrix dimensions to initialise it
               if (allocated(array)) deallocate(array)
               allocate(array(array_size(1):array_size(2)))
            endif

c           Now that all the arrays are initialised, broadcast data
            call MPI_BCAST(array,array_size(2)-array_size(1)+1,
     &                     rp,root,cw,ierror)

         else

            if (allocated(array)) deallocate(array)

         endif

         end subroutine double_array_broadcast

c        **************************************************************
c        **                                                          **
c        **  Broadcast double prec array from node 0 to all the      **
c        **  other nodes                                             **
c        **                                                          **
c        **************************************************************
         subroutine double_pointer_broadcast(array)

         use mpi

         implicit none

         integer :: status(MPI_STATUS_SIZE)
         double precision, dimension(:), pointer,
     &                                     intent(inout) :: array
         integer :: root, node, ierror, nnodes
         integer :: array_size(2), j, rp, tag, cw, inode, mi, l1, u1
         external :: exita

         root = 0
         mi   = MPI_INTEGER
         rp   = real_precision_flag(0.0_dp)
         cw   = MPI_COMM_WORLD

c        Get current processor id
         call MPI_COMM_RANK(MPI_COMM_WORLD, node, ierror)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, nnodes, ierror)

c        Matrix doesn't need to be initialised
         masternode: if (node == 0) then

            if (.not.associated(array)) then
               if (verbose) write(*,*)'Warning: broadcast unallocated '
     &                             // 'array in double_pointe_broadcast'

               array_size = [0,0]

            else

               array_size(1)  = lbound(array,1)
               array_size(2)  = ubound(array,1)

            endif

         endif masternode

c        Broadcast array dimensions to other nodes
         call MPI_BCAST(array_size,2,mi,root,cw,ierror)

         if (sum(array_size)>0) then

            if (node /= 0) then

c           Other nodes receive matrix dimensions to initialise it
               array => null()
               allocate(array(array_size(1):array_size(2)))
            endif

c           Now that all the arrays are initialised, broadcast data
            call MPI_BCAST(array,array_size(2)-array_size(1)+1,
     &                     rp,root,cw,ierror)

         else

            array => null()

         endif

         end subroutine double_pointer_broadcast

         subroutine arbitrary_pointer_broadcast(array)

         use mpi

         implicit none

         integer :: status(MPI_STATUS_SIZE)
         real (dp)       , dimension(:), pointer,
     &                                     intent(inout) :: array
         integer :: root, node, ierror, nnodes
         integer :: array_size(2), j, rp, tag, cw, inode, mi, l1, u1
         external :: exita

         root = 0
         mi   = MPI_INTEGER
         rp   = real_precision_flag(0.0_dp)
         cw   = MPI_COMM_WORLD

c        Get current processor id
         call MPI_COMM_RANK(MPI_COMM_WORLD, node, ierror)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, nnodes, ierror)

c        Matrix doesn't need to be initialised
         masternode: if (node == 0) then

            if (.not.associated(array)) then
               if (verbose) write(*,*)'Warning: broadcast unallocated '
     &                             // 'array in double_pointe_broadcast'

               array_size = [0,0]

            else

               array_size(1)  = lbound(array,1)
               array_size(2)  = ubound(array,1)

            endif

         endif masternode

c        Broadcast array dimensions to other nodes
         call MPI_BCAST(array_size,2,mi,root,cw,ierror)

         if (sum(array_size)>0) then

            if (node /= 0) then

c           Other nodes receive matrix dimensions to initialise it
               array => null()
               allocate(array(array_size(1):array_size(2)))
            endif

c           Now that all the arrays are initialised, broadcast data
            call MPI_BCAST(array,array_size(2)-array_size(1)+1,
     &                     rp,root,cw,ierror)

         else

            array => null()

         endif

         end subroutine arbitrary_pointer_broadcast


c        **************************************************************
c        **                                                          **
c        **  Broadcast integer array from node 0 to all the          **
c        **  other nodes                                             **
c        **                                                          **
c        **************************************************************
         subroutine integer_array_broadcast(array)

            use mpi
         implicit none

         integer :: status(MPI_STATUS_SIZE)
         integer, dimension(:), allocatable, intent(inout) :: array
         integer :: root, node, ierror, nnodes
         integer :: array_size(2), j, rp, tag, cw, inode, mi
         external :: exita

         root = 0
         mi   = MPI_INTEGER
         rp   = real_precision_flag(0.0_dp)
         cw   = MPI_COMM_WORLD

c        Get current processor id
         call MPI_COMM_RANK(MPI_COMM_WORLD, node, ierror)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, nnodes, ierror)

c        Matrix doesn't need to be initialised
         masternode: if (node == 0) then

            if (.not.allocated(array)) then

               if (verbose) write(*,*)'Warning: broadcast unallocated '
     &                             //'array in integer_array_broadcast'
               array_size = [0,0]
            else
               array_size(1) = lbound(array,1)
               array_size(2) = ubound(array,1)
            endif

         endif masternode

c           Broadcast array dimensions to other nodes
            call MPI_BCAST(array_size,2,mi,root,cw,ierror)

            if (sum(array_size)>0) then

            if (node /= 0) then

c           Other nodes receive matrix dimensions to initialise it
               if (allocated(array)) deallocate(array)
               allocate(array(array_size(1):array_size(2)))
            endif

c           Now that all the arrays are initialised, broadcast data
            call MPI_BCAST(array,array_size(2)-array_size(1)+1,
     &                     mi,root,cw,ierror)

         else
            if (allocated(array)) deallocate(array)
         endif

         end subroutine integer_array_broadcast



c        **************************************************************
c        **                                                          **
c        **  Broadcast integer value from node 0 to all the          **
c        **  other nodes                                             **
c        **                                                          **
c        **************************************************************
         subroutine integer_broadcast(variable)

         use mpi

         implicit none

         integer :: status(MPI_STATUS_SIZE)
         integer, intent(inout) :: variable
         integer :: root, ierror, rp, cw, mi

         root = 0
         mi   = MPI_INTEGER
         rp   = real_precision_flag(0.0_dp)
         cw   = MPI_COMM_WORLD

c        Broadcast datum
         call MPI_BCAST(variable,1,mi,root,cw,ierror)

         end subroutine integer_broadcast


c        **************************************************************
c        **                                                          **
c        **  Broadcast double precision value from node 0 to all the **
c        **  other nodes                                             **
c        **                                                          **
c        **************************************************************
         subroutine double_broadcast(variable)

         use mpi

         implicit none

         integer :: status(MPI_STATUS_SIZE)
         double precision, intent(inout) :: variable
         integer :: root, ierror, rp, cw

         root = 0
         rp   = MPI_DOUBLE_PRECISION
         cw   = MPI_COMM_WORLD

c        Broadcast datum
         call MPI_BCAST(variable,1,rp,root,cw,ierror)

         end subroutine double_broadcast

         subroutine arbitrary_broadcast(variable)

         use mpi

         implicit none

         integer :: status(MPI_STATUS_SIZE)
         real (dp)       , intent(inout) :: variable
         integer :: root, ierror, rp, cw

         root = 0
         rp   = real_precision_flag(0.0_dp)
         cw   = MPI_COMM_WORLD

c        Broadcast datum
         call MPI_BCAST(variable,1,rp,root,cw,ierror)

         end subroutine arbitrary_broadcast

c        **************************************************************
c        **                                                          **
c        **  Broadcast logical value from node 0 to all the nodes    **
c        **                                                          **
c        **************************************************************
         subroutine logical_broadcast(variable)

         use mpi

         implicit none

         integer :: status(MPI_STATUS_SIZE)
         logical, intent(inout) :: variable
         integer :: root, ierror, ml, cw

         root = 0
         ml   = MPI_LOGICAL
         cw   = MPI_COMM_WORLD

c        Broadcast datum
         call MPI_BCAST(variable,1,ml,root,cw,ierror)

         end subroutine logical_broadcast


c        **************************************************************
c        **                                                          **
c        **  Broadcast character string from node 0 to all the nodes **
c        **                                                          **
c        **************************************************************
         subroutine character_broadcast(chlen,variable)

         use mpi

         implicit none

         integer :: status(MPI_STATUS_SIZE)
         integer, intent(in) :: chlen
         character(len=chlen), intent(inout) :: variable
         integer :: root, ierror, mc, cw

         root = 0
         mc   = MPI_CHARACTER
         cw   = MPI_COMM_WORLD

c        Broadcast datum
         call MPI_BCAST(variable,chlen,mc,root,cw,ierror)

         end subroutine character_broadcast


c        **************************************************************
c        **                                                          **
c        **  Broadcast logical array from node 0 to all the          **
c        **  other nodes                                             **
c        **                                                          **
c        **************************************************************
         subroutine logical_array_broadcast(array)

         use mpi

         implicit none

         integer :: status(MPI_STATUS_SIZE)
         logical, dimension(:), allocatable, intent(inout) :: array
         integer :: root, node, ierror, nnodes
         integer :: array_size(2), j, rp, tag, cw, inode, mi
         external :: exita

         root = 0
         mi   = MPI_LOGICAL
         cw   = MPI_COMM_WORLD

c        Get current processor id
         call MPI_COMM_RANK(MPI_COMM_WORLD, node, ierror)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, nnodes, ierror)

c        Matrix doesn't need to be initialised
         masternode: if (node == 0) then

            if (.not.allocated(array)) then
               if (verbose) write(*,*)'Warning: broadcast unallocated '
     &                             //'array in logical_array_broadcast'
               array_size = [0,0]
            else
               array_size(1)  = lbound(array,1)
               array_size(2)  = ubound(array,1)
            endif

         endif masternode


c        Broadcast array dimensions to other nodes
         call MPI_BCAST(array_size,2,mi,root,cw,ierror)

         if (sum(array_size)>0) then

            if (node /= 0) then

c           Other nodes receive matrix dimensions to initialise it
               if (allocated(array)) deallocate(array)
               allocate(array(array_size(1):array_size(2)))
            endif

c           Now that all the arrays are initialised, broadcast data
            call MPI_BCAST(array,array_size(2)-array_size(1)+1,
     &                     mi,root,cw,ierror)

         else

            if (allocated(array)) deallocate(array)

         endif

         end subroutine logical_array_broadcast


c        **************************************************************
c        **                                                          **
c        **  Broadcast character array from node 0 to all the        **
c        **  other nodes                                             **
c        **                                                          **
c        **************************************************************
         subroutine character_array_broadcast(chlen,array)

         use mpi

         implicit none

         integer :: status(MPI_STATUS_SIZE)
         integer, intent(in) :: chlen
         character(len=chlen), dimension(:),
     &                         allocatable , intent(inout) :: array
         integer :: root, node, ierror, nnodes
         integer :: array_size(2), j, rp, tag, cw, inode, mi, mc
         external :: exita

         root = 0
         mi   = MPI_INTEGER
         mc   = MPI_CHARACTER
         cw   = MPI_COMM_WORLD

c        Get current processor id
         call MPI_COMM_RANK(MPI_COMM_WORLD, node, ierror)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, nnodes, ierror)

c        Matrix doesn't need to be initialised
         masternode: if (node == 0) then

            if (.not.allocated(array)) then

               if (verbose) write(*,*)'Warning: broadcast unallocated '
     &                          // 'array in character_array_broadcast'
               array_size = [0,0]

            else

               array_size(1)  = lbound(array,1)
               array_size(2)  = ubound(array,1)

            endif

         endif masternode

c        Broadcast array dimensions to other nodes
         call MPI_BCAST(array_size,2,mi,root,cw,ierror)

         if (sum(array_size)>0) then

         if (node /= 0) then

c        Other nodes receive matrix dimensions to initialise it
            if (allocated(array)) deallocate(array)
            allocate(array(array_size(1):array_size(2)))
         endif

c        Now that all the arrays are initialised, broadcast data
         do j = array_size(1), array_size(2)
           call MPI_BCAST(array(j)(1:chlen),chlen,mc,root,cw,ierror)
         end do

         else

           if (allocated(array)) deallocate(array)

         endif



         end subroutine character_array_broadcast


c        **************************************************************
c        **                                                          **
c        **  Broadcast fixed dimension integer array from node 0 to  **
c        **  all the other nodes                                     **
c        **                                                          **
c        **************************************************************
         subroutine fixed_integer_array_broadcast(width,array)

         use mpi

         implicit none

         integer :: status(MPI_STATUS_SIZE)
         integer,                   intent(in)    :: width
         integer, dimension(width), intent(inout) :: array
         integer :: root, node, ierror, nnodes
         integer :: array_size(2), j, rp, tag, cw, inode, mi
         external :: exita

         root = 0
         mi   = MPI_INTEGER
         rp   = real_precision_flag(0.0_dp)
         cw   = MPI_COMM_WORLD

c        Get current processor id
         call MPI_COMM_RANK(MPI_COMM_WORLD, node,   ierror)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, nnodes, ierror)

c        Matrix doesn't need to be initialised
         masternode: if (node == 0) then

            array_size(1) = lbound(array,1)
            array_size(2) = ubound(array,1)

         endif masternode

c        Broadcast array dimensions to other nodes
         call MPI_BCAST(array_size,2,mi,root,cw,ierror)

c        Now that all the arrays are initialised, broadcast data
         call MPI_BCAST(array,array_size(2)-array_size(1)+1,
     &                  mi,root,cw,ierror)

         end subroutine fixed_integer_array_broadcast


c        **************************************************************
c        **                                                          **
c        **  Broadcast character array from node 0 to all the        **
c        **  other nodes                                             **
c        **                                                          **
c        **************************************************************
         subroutine fixed_character_array_broadcast(chlen,width,array)

         use mpi

         implicit none

         integer :: status(MPI_STATUS_SIZE)
         integer, intent(in) :: chlen, width
         character(len=chlen), dimension(width), intent(inout) :: array
         integer :: root, node, ierror, nnodes
         integer :: array_size(2), j, rp, tag, cw, inode, mi, mc
         external :: exita

         root = 0
         mi   = MPI_INTEGER
         mc   = MPI_CHARACTER
         cw   = MPI_COMM_WORLD

c        Get current processor id
         call MPI_COMM_RANK(MPI_COMM_WORLD, node, ierror)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, nnodes, ierror)

c        Matrix doesn't need to be initialised
         masternode: if (node == 0) then

            array_size(1)  = lbound(array,1)
            array_size(2)  = ubound(array,1)

         endif masternode

c        Broadcast array dimensions to other nodes
         call MPI_BCAST(array_size,2,mi,root,cw,ierror)

c        Now that all the arrays are initialised, broadcast data
         do j = array_size(1), array_size(2)
           call MPI_BCAST(array(j)(1:chlen),chlen,mc,root,cw,ierror)
         end do

         end subroutine fixed_character_array_broadcast

c        **************************************************************
c        **                                                          **
c        **  Broadcast double prec array from node 0 to all the      **
c        **  other nodes                                             **
c        **                                                          **
c        **************************************************************
         subroutine fixed_double_array_broadcast(width,array)

         use mpi

         implicit none

         integer, intent(in) :: width
         integer :: status(MPI_STATUS_SIZE)
         real (dp), dimension(width), intent(inout) :: array
         integer :: root, node, ierror, nnodes
         integer :: array_size(2), j, rp, tag, cw, inode, mi, l1, u1
         external :: exita

         root = 0
         mi   = MPI_INTEGER
         rp   = real_precision_flag(0.0_dp)
         cw   = MPI_COMM_WORLD

c        Get current processor id
         call MPI_COMM_RANK(MPI_COMM_WORLD, node, ierror)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, nnodes, ierror)

c        Matrix doesn't need to be initialised
         masternode: if (node == 0) then

            array_size(1)  = lbound(array,1)
            array_size(2)  = ubound(array,1)

         endif masternode

c        Broadcast array dimensions to other nodes
         call MPI_BCAST(array_size,2,mi,root,cw,ierror)

c        Now that all the arrays are initialised, broadcast data
         call MPI_BCAST(array,array_size(2)-array_size(1)+1,
     &                  rp,root,cw,ierror)

         end subroutine fixed_double_array_broadcast

c        **************************************************************
c        **                                                          **
c        **  Returns the current node number (from 0 to nprocs-1)    **
c        **                                                          **
c        **************************************************************
         function current_MPI_CPUID() result(nidproc)

         use mpi

         implicit none

         integer :: status(MPI_STATUS_SIZE)
         integer :: nidproc
         integer :: ierror

c        Get current processor id
         call MPI_COMM_RANK(MPI_COMM_WORLD, nidproc, ierror)

         if (ierror/=0) then
            write(*,*)'Error in current_MPI_CPUID: ierror = ',ierror
            stop
         endif

         end function current_MPI_CPUID

c        **************************************************************
c        **                                                          **
c        **  Returns the total number of cpus                        **
c        **                                                          **
c        **************************************************************
         function total_MPI_CPUs() result(noprocs)

         use mpi

         implicit none

         integer :: status(MPI_STATUS_SIZE)
         integer :: noprocs
         integer :: ierror

c        Get total number of processors
         call MPI_COMM_SIZE(MPI_COMM_WORLD, noprocs, ierror)

         if (ierror/=0) then
            write(*,*)'Error in total_MPI_CPUs: ierror = ',ierror
            stop
         endif

         end function total_MPI_CPUs




      end module sparse_MPI










