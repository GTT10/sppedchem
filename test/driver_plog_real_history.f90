program driver_plog_real_history
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use working_precision, only: dp
   use chemistry_setup, only: mechdir, use_speedchem, solver
   use speedchem, only: ns, nr, neq, specie, molefr_to_massfr,       &
                        massfr_to_molefr
   use SCmixturethermo, only: SCrho, SCP, rhoY, molar_volumes,       &
                              pressurerhoT
   implicit none

   real(dp), allocatable :: y(:), atol(:), x(:)
   real(dp) :: t0, target, tend, temperature0, pressure0
   real(dp) :: prev_time, prev_temp, max_dTdt, dTdt, idt50, frac
   real(dp) :: wall0, wall1, mass_sum, max_mass_error
   integer :: i, nout, ih2, io2, iar, ih, ioh, iho2, ih2o
   character(len=256) :: env
   integer :: env_len, env_stat

   call get_environment_variable('SC_MECHDIR', env, env_len, env_stat)
   if (env_stat /= 0 .or. env_len == 0) then
      write(*,'(a)') 'ERROR: SC_MECHDIR is required'
      error stop 1
   endif
   mechdir = trim(env)
   if (mechdir(len_trim(mechdir):len_trim(mechdir)) /= '/')           &
      mechdir = trim(mechdir)//'/'
   use_speedchem = .true.

   temperature0 = env_real('SC_T0', 1200.0_dp)
   pressure0    = env_real('SC_P0', 10.0_dp*101325.0_dp)
   tend         = env_real('SC_TEND', 2.0e-3_dp)
   nout         = env_integer('SC_NOUT', 400)
   if (temperature0 <= 0.0_dp .or. pressure0 <= 0.0_dp .or.          &
       tend <= 0.0_dp .or. nout < 2) then
      write(*,'(a)') 'ERROR: invalid real-mechanism history settings'
      error stop 1
   endif

   call chemistry_input

   ih2  = species_index('H2')
   io2  = species_index('O2')
   iar  = species_index('AR')
   ih   = species_index('H')
   ioh  = species_index('OH')
   iho2 = species_index('HO2')
   ih2o = species_index('H2O')
   if (min(ih2, io2, iar, ih, ioh, iho2, ih2o) == 0) then
      write(*,'(a)') 'ERROR: required C3Mech history species are missing'
      error stop 1
   endif

   allocate(y(neq), atol(neq), x(ns))
   x = 0.0_dp
   x(ih2) = 2.0_dp
   x(io2) = 1.0_dp
   x(iar) = 7.0_dp
   x = x/sum(x)

   y = 0.0_dp
   y(1) = temperature0
   call molefr_to_massfr(x, y(2:neq))
   atol = 1.0e-18_dp
   atol(1) = 1.0e-8_dp

!  Constant-volume state: establish density from the requested initial
!  T/P/composition, then keep SCrho fixed throughout integration.
   SCP = pressure0
   call rhoY(y(2:neq), y(1))
   call molar_volumes

   write(*,'(a,a)') '# REAL_PLOG solver=',trim(solver)
   write(*,'(a,i0,a,i0,a,es24.16,a,es24.16,a,es24.16)')               &
      '# REAL_PLOG ns=',ns,' nr=',nr,' rho=',SCrho,' T0=',temperature0,&
      ' P0=',pressure0
   write(*,'(a)')                                                      &
      'kind,time_s,T_K,P_Pa,X_H2,X_O2,X_H,X_OH,X_HO2,X_H2O,X_AR,sumY'

   t0 = 0.0_dp
   prev_time = 0.0_dp
   prev_temp = y(1)
   max_dTdt = -huge(1.0_dp)
   idt50 = -1.0_dp
   max_mass_error = 0.0_dp
   call cpu_time(wall0)
   call emit_history(t0)

   do i = 1, nout
      target = tend*real(i,dp)/real(nout,dp)
      call chemistry_ODE_integrate(neq, 1.0e-8_dp, atol, t0, target, y)
      if (.not. all(ieee_is_finite(y))) then
         write(*,'(a)') 'ERROR: non-finite state in real PLOG history'
         error stop 1
      endif
      if (t0 <= prev_time) then
         write(*,'(a)') 'ERROR: integrator did not advance time'
         error stop 1
      endif
      dTdt = (y(1)-prev_temp)/(t0-prev_time)
      max_dTdt = max(max_dTdt,dTdt)
      if (idt50 < 0.0_dp .and. y(1) >= temperature0+50.0_dp) then
         frac = (temperature0+50.0_dp-prev_temp)/(y(1)-prev_temp)
         idt50 = prev_time + frac*(t0-prev_time)
      endif
      prev_time = t0
      prev_temp = y(1)
      call emit_history(t0)
   enddo
   call cpu_time(wall1)

   write(*,'(a,",",a,4(",",a,",",es24.16))') 'SUMMARY',           &
      trim(solver), 'idt50_s', idt50, 'max_dTdt_Kps', max_dTdt,       &
      'cpu_s', wall1-wall0, 'max_sumY_error', max_mass_error

contains

   subroutine emit_history(time)
      real(dp), intent(in) :: time
      real(dp) :: pressure
      x = massfr_to_molefr(y(2:neq))
      pressure = pressurerhoT(y(1),y(2:neq))
      mass_sum = sum(y(2:neq))
      max_mass_error = max(max_mass_error,abs(mass_sum-1.0_dp))
      write(*,'(a,11(",",es24.16))') 'HIST',time,y(1),pressure,       &
         x(ih2),x(io2),x(ih),x(ioh),x(iho2),x(ih2o),x(iar),mass_sum
   end subroutine emit_history

   integer function species_index(name) result(idx)
      character(len=*), intent(in) :: name
      integer :: k
      idx = 0
      do k=1,ns
         if (trim(specie(k)) == name) then
            idx = k
            return
         endif
      enddo
   end function species_index

   real(dp) function env_real(name, default_value) result(value)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: default_value
      character(len=64) :: raw
      integer :: length, status, ios
      value = default_value
      call get_environment_variable(name,raw,length,status)
      if (status == 0 .and. length > 0) then
         read(raw,*,iostat=ios) value
         if (ios /= 0) error stop 1
      endif
   end function env_real

   integer function env_integer(name, default_value) result(value)
      character(len=*), intent(in) :: name
      integer, intent(in) :: default_value
      character(len=64) :: raw
      integer :: length, status, ios
      value = default_value
      call get_environment_variable(name,raw,length,status)
      if (status == 0 .and. length > 0) then
         read(raw,*,iostat=ios) value
         if (ios /= 0) error stop 1
      endif
   end function env_integer

end program driver_plog_real_history
