program driver_long_strings
   use mpi
   use chemistry_string_limits, only: species_name_len
   use chemistry_setup, only: mechdir, use_speedchem
   use speedchem, only: ns, nr, specie
   use sparse_chemistry, only: stoich_r_sp, stoich_p_sp
   implicit none

   integer, parameter :: nlong = 6
   character(len=species_name_len), parameter :: expected(nlong) = [ &
      'LONGSPECIES0000001', 'LONGSPECIES0000002',                    &
      'LONGSPECIES0000003', 'LONGSPECIES0000004',                    &
      'LONGSPECIES0000005', 'LONGSPECIES0000006' ]
   character(len=256) :: mech
   integer :: env_len, env_stat, ierr, rank, nproc
   integer :: i, j, local_fail, global_fail
   logical :: found

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

   if (rank == 0) then
      call get_environment_variable('SC_MECHDIR', mech, env_len, env_stat)
      if (env_stat /= 0 .or. env_len == 0) mech = 'test/data_long_strings/'
      if (mech(len_trim(mech):len_trim(mech)) /= '/') mech = trim(mech)//'/'
      mechdir = mech
      use_speedchem = .true.
      call chemistry_input
   endif

   call SCbroadcast

   local_fail = 0
   if (species_name_len /= 18 .or. len(specie) /= species_name_len) local_fail = 1
   if (ns /= nlong + 1 .or. nr /= 1) local_fail = 1
   if (stoich_r_sp%n /= 3 .or. stoich_p_sp%n /= 3) local_fail = 1

   do i = 1, nlong
      found = .false.
      do j = 1, ns
         if (specie(j) == expected(i)) then
            found = .true.
            exit
         endif
      enddo
      if (.not. found) local_fail = 1
   enddo

   call MPI_Allreduce(local_fail, global_fail, 1, MPI_INTEGER, MPI_MAX, &
                      MPI_COMM_WORLD, ierr)
   if (global_fail /= 0) then
      if (rank == 0) write(*,'(a)') &
         'RESULT: FAIL - long chemistry strings did not round-trip'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
   endif

   if (rank == 0) write(*,'(a,i0,a)') &
      'RESULT: PASS - 18-character species and >80-character records on ', &
      nproc, ' MPI rank(s)'
   call MPI_Finalize(ierr)
end program driver_long_strings
