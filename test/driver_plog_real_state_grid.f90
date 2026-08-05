program driver_plog_real_state_grid
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use working_precision, only: dp
   use chemistry_setup, only: mechdir, use_speedchem
   use speedchem, only: ns, neq, specie, molefr_to_massfr
   use speedchem_conV, only: SC_conV
   use SCmixturethermo, only: SCrho, SCP, rhoY, molar_volumes,        &
                              pressurerhoT
   implicit none

   real(dp), parameter :: atm_to_pa = 101325.0_dp
   real(dp), parameter :: temperatures(3) = [ &
      700.0_dp, 1000.0_dp, 1500.0_dp ]
   real(dp), parameter :: pressures_atm(11) = [ &
      0.005_dp, 0.01_dp, 0.03162277660168379_dp, 0.1_dp, &
      1.0_dp, 5.0_dp, 7.071067811865476_dp, 10.0_dp, &
      30.0_dp, 100.0_dp, 300.0_dp ]

   character(len=256) :: env
   integer :: env_len, env_stat
   integer :: i, it, ip, icase, ih2, io2, iar, ih, iho2
   integer :: nfinite, nnonfinite
   real(dp) :: requested_pressure, reconstructed_pressure
   real(dp), allocatable :: x(:), ymass(:), state(:), rhs(:)

   call get_environment_variable('SC_MECHDIR', env, env_len, env_stat)
   if (env_stat /= 0 .or. env_len == 0) then
      write(*,'(a)') 'ERROR: SC_MECHDIR is required'
      error stop 1
   endif
   mechdir = trim(env)
   if (mechdir(len_trim(mechdir):len_trim(mechdir)) /= '/')           &
      mechdir = trim(mechdir)//'/'
   use_speedchem = .true.
   call chemistry_input

   ih2  = species_index('H2')
   io2  = species_index('O2')
   iar  = species_index('AR')
   ih   = species_index('H')
   iho2 = species_index('HO2')
   if (min(ih2, io2, iar, ih, iho2) == 0) then
      write(*,'(a)') 'ERROR: C3Mech H2/O2/AR/H/HO2 species are required'
      error stop 1
   endif

   allocate(x(ns), ymass(ns), state(neq), rhs(neq))
   x = 0.0_dp
   x(ih2)  = 2.0_dp
   x(io2)  = 1.0_dp
   x(iar)  = 7.0_dp
   x(ih)   = 1.0e-5_dp
   x(iho2) = 1.0e-8_dp
   x = x/sum(x)
   call molefr_to_massfr(x, ymass)

   if (abs(sum(ymass)-1.0_dp) > 1.0e-12_dp) then
      write(*,'(a,es24.16)') 'ERROR: invalid mass-fraction sum: ',     &
         sum(ymass)
      error stop 1
   endif

   write(*,'(a)',advance='no') 'SPECIES'
   do i = 1, ns
      write(*,'(",",a)',advance='no') trim(specie(i))
   enddo
   write(*,*)
   write(*,'(a,5(",",es24.16))') 'MIXTURE_X_H2_O2_AR_H_HO2',         &
      x(ih2), x(io2), x(iar), x(ih), x(iho2)
   write(*,'(a)') 'kind,case,T_K,P_Pa,rho_kgpm3,dTdt_Kps,dYdt_1,...'

   icase = 0
   nfinite = 0
   nnonfinite = 0
   do it = 1, size(temperatures)
      do ip = 1, size(pressures_atm)
         icase = icase + 1
         requested_pressure = pressures_atm(ip)*atm_to_pa

         state(1) = temperatures(it)
         state(2:neq) = ymass
         SCP = requested_pressure
         call rhoY(state(2:neq), state(1))
         call molar_volumes
         reconstructed_pressure = pressurerhoT(state(1), state(2:neq))

         if (abs(reconstructed_pressure/requested_pressure-1.0_dp) >  &
             1.0e-10_dp) then
            write(*,'(a,i0,2(a,es24.16))')                            &
               'ERROR: pressure reconstruction failed at case ',      &
               icase, ' requested=', requested_pressure,              &
               ' reconstructed=', reconstructed_pressure
            error stop 1
         endif

         call SC_conV(neq, 0.0_dp, state, rhs)
         if (.not. all(ieee_is_finite(rhs))) then
            nnonfinite = nnonfinite + 1
            write(*,'(a,",",i0,3(",",es24.16))',advance='no')       &
               'NONFINITE', icase, state(1), reconstructed_pressure, SCrho
            do i = 1, neq
               write(*,'(",",l1)',advance='no') ieee_is_finite(rhs(i))
            enddo
            write(*,*)
            cycle
         endif

         nfinite = nfinite + 1
         write(*,'(a,",",i0,4(",",es24.16))',advance='no') 'STATE', &
            icase, state(1), reconstructed_pressure, SCrho, rhs(1)
         do i = 2, neq
            write(*,'(",",es24.16)',advance='no') rhs(i)
         enddo
         write(*,*)
      enddo
   enddo

   write(*,'(a,i0,a,i0)') 'STATE_GRID_SUMMARY finite=', nfinite,      &
      ' nonfinite=', nnonfinite
   if (nfinite == 0) error stop 1

contains

   integer function species_index(name) result(idx)
      character(len=*), intent(in) :: name
      integer :: k
      idx = 0
      do k = 1, ns
         if (trim(specie(k)) == name) then
            idx = k
            return
         endif
      enddo
   end function species_index

end program driver_plog_real_state_grid
