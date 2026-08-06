program driver_plog_mpi_real
   use mpi
   use working_precision, only: dp
   use chemistry_setup, only: mechdir, use_speedchem
   use speedchem, only: ns, nr, neq, specie, molefr_to_massfr
   use reacpar, only: n_plog_reactions, n_plog_nodes, plog_reaction, &
                       plog_kinf_eval
   use speedchem_conV, only: SC_conV, constV_jac_sparse
   use sparse_chemistry, only: JAC_sparse
   use SCmixturethermo, only: SCrho, SCP, rhoY, molar_volumes,        &
                              pressurerhoT
   implicit none

   integer :: ierr, rank, nproc, i, ih2, io2, iar, ih, iho2
   integer :: env_len, env_stat
   character(len=256) :: env
   real(dp), allocatable :: y(:), x(:), rhs(:), kinf(:)
   real(dp) :: ta(6), pressure
   real(dp) :: local(10), vmin(10), vmax(10)

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
   call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
   if (rank == 0) then
      call get_environment_variable('SC_MECHDIR',env,env_len,env_stat)
      if (env_stat /= 0 .or. env_len == 0) then
         write(*,'(a)') 'ERROR: SC_MECHDIR is required'
         call MPI_Abort(MPI_COMM_WORLD,1,ierr)
      endif
      mechdir = trim(env)
      if (mechdir(len_trim(mechdir):len_trim(mechdir)) /= '/')        &
         mechdir = trim(mechdir)//'/'
      use_speedchem = .true.
      call chemistry_input
   endif
   call SCbroadcast

   allocate(y(neq),x(ns),rhs(neq),kinf(nr))
   ih2=0; io2=0; iar=0; ih=0; iho2=0
   do i=1,ns
      if (trim(specie(i)) == 'H2')  ih2=i
      if (trim(specie(i)) == 'O2')  io2=i
      if (trim(specie(i)) == 'AR')  iar=i
      if (trim(specie(i)) == 'H')   ih=i
      if (trim(specie(i)) == 'HO2') iho2=i
   enddo
   if (min(ih2,io2,iar,ih,iho2) == 0) then
      if (rank == 0) write(*,'(a)') 'ERROR: required C3Mech species missing'
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
   endif

   x=0.0_dp
   x(ih2)=2.0_dp
   x(io2)=1.0_dp
   x(iar)=7.0_dp
   x(ih)=1.0e-5_dp
   x(iho2)=1.0e-8_dp
   x=x/sum(x)
   y=0.0_dp
   y(1)=1000.0_dp
   call molefr_to_massfr(x,y(2:neq))
   SCP=7.071067811865476_dp*101325.0_dp
   call rhoY(y(2:neq),y(1))
   call molar_volumes
   pressure=pressurerhoT(y(1),y(2:neq))

   ta(1)=y(1)
   ta(2)=ta(1)*ta(1)
   ta(3)=ta(2)*ta(1)
   ta(4)=ta(3)*ta(1)
   ta(5)=1.0_dp/ta(1)
   ta(6)=log(ta(1))
   kinf=0.0_dp
   call plog_kinf_eval(ta,pressure,kinf)
   call SC_conV(neq,0.0_dp,y,rhs)
   call constV_jac_sparse(neq,0.0_dp,y)

   local = [real(n_plog_reactions,dp), real(n_plog_nodes,dp),         &
            sum(kinf(plog_reaction)), sum(abs(kinf(plog_reaction))),  &
            sum(rhs), sum(abs(rhs)), sum(JAC_sparse%A),               &
            sum(abs(JAC_sparse%A)), rhs(1),                           &
            real(JAC_sparse%n,dp)]
   call MPI_Allreduce(local,vmin,size(local),MPI_DOUBLE_PRECISION,     &
                      MPI_MIN,MPI_COMM_WORLD,ierr)
   call MPI_Allreduce(local,vmax,size(local),MPI_DOUBLE_PRECISION,     &
                      MPI_MAX,MPI_COMM_WORLD,ierr)
   if (nproc < 2 .or. any(vmin /= vmax) .or.                         &
       n_plog_reactions /= 1 .or. n_plog_nodes /= 11) then
      if (rank == 0) write(*,'(a)')                                   &
         'RESULT: FAIL - public real-PLOG MPI rate/Jacobian mismatch'
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
   endif
   if (rank == 0) write(*,'(a,i0,a,es12.4,a,i0)')                     &
      'RESULT: PASS - public real-PLOG rates/Jacobian identical on ', &
      nproc,' ranks; dTdt=',rhs(1),' jac_nnz=',JAC_sparse%n
   call MPI_Finalize(ierr)
end program driver_plog_mpi_real
