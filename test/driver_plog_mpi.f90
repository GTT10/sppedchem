program driver_plog_mpi
   use mpi
   use working_precision, only: dp
   use chemistry_setup, only: mechdir, use_speedchem
   use reacpar, only: n_plog_reactions, n_plog_nodes, plog_reaction, &
                       plog_node_ptr, plog_logP, plog_A, plog_b,     &
                       plog_EoverR
   implicit none
   integer :: ierr, rank, nproc
   real(dp) :: local(7), vmin(7), vmax(7)

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
   call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
   if (rank == 0) then
      mechdir = 'test/data_plog/'
      use_speedchem = .true.
      call chemistry_input
   endif
   call SCbroadcast

   local = [real(n_plog_reactions,dp), real(n_plog_nodes,dp),         &
            real(sum(plog_reaction),dp), real(sum(plog_node_ptr),dp), &
            sum(plog_logP), sum(log(plog_A)), sum(plog_b+plog_EoverR)]
   call MPI_Allreduce(local,vmin,size(local),MPI_DOUBLE_PRECISION,     &
                      MPI_MIN,MPI_COMM_WORLD,ierr)
   call MPI_Allreduce(local,vmax,size(local),MPI_DOUBLE_PRECISION,     &
                      MPI_MAX,MPI_COMM_WORLD,ierr)
   if (nproc < 2 .or. any(vmin /= vmax) .or.                         &
       n_plog_reactions /= 1 .or. n_plog_nodes /= 3) then
      if (rank == 0) write(*,'(a)') 'RESULT: FAIL - MPI PLOG broadcast mismatch'
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
   endif
   if (rank == 0) write(*,'(a,i0,a)')                                &
      'RESULT: PASS - PLOG data identical on ',nproc,' MPI ranks'
   call MPI_Finalize(ierr)
end program driver_plog_mpi
