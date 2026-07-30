program driver_plog_integrate
   use working_precision, only: dp
   use chemistry_setup, only: mechdir, use_speedchem, solver
   use speedchem, only: ns, neq, specie
   use SCmixturethermo, only: SCrho, molar_volumes
   use ieee_arithmetic, only: ieee_is_finite
   implicit none
   real(dp), allocatable :: y(:), atol(:)
   real(dp) :: t0, tf
   integer :: i, ih, io2, in2

   mechdir = 'test/data_plog/'
   use_speedchem = .true.
   call chemistry_input
   if (trim(solver) /= 'LSODESJAC') error stop 1

   ih=0; io2=0; in2=0
   do i=1,ns
      if (trim(specie(i)) == 'H') ih=i
      if (trim(specie(i)) == 'O2') io2=i
      if (trim(specie(i)) == 'N2') in2=i
   enddo
   if (min(ih,io2,in2) == 0) error stop 1
   allocate(y(neq),atol(neq))
   y=0.0_dp; y(1)=1100.0_dp
   y(ih+1)=1.0e-3_dp; y(io2+1)=0.23_dp; y(in2+1)=0.769_dp
   atol=1.0e-15_dp
   SCrho=1.0_dp
   call molar_volumes
   t0=0.0_dp; tf=1.0e-8_dp
   call chemistry_ODE_integrate(neq,1.0e-7_dp,atol,t0,tf,y)
   if (.not.all(ieee_is_finite(y)) .or. t0 < tf) then
      write(*,'(a)') 'RESULT: FAIL - PLOG LSODESJAC integration'
      error stop 1
   endif
   write(*,'(a,es12.4)') 'RESULT: PASS - PLOG LSODESJAC integration T=',y(1)
end program driver_plog_integrate
