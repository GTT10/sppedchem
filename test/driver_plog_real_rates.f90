program driver_plog_real_rates
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use working_precision, only: dp
   use chemistry_setup, only: mechdir, use_speedchem
   use speedchem, only: nr
   use reacpar, only: n_plog_reactions, n_plog_nodes, plog_reaction, &
                      plog_kinf_eval
   implicit none

   real(dp), parameter :: atm_to_pa = 101325.0_dp
   real(dp), parameter :: temperatures(5) = [ &
      500.0_dp, 700.0_dp, 1000.0_dp, 1500.0_dp, 2500.0_dp ]
   real(dp), parameter :: pressures_atm(17) = [ &
      0.005_dp, 0.01_dp, 0.03162277660168379_dp, 0.1_dp, &
      0.31622776601683794_dp, 1.0_dp, 2.23606797749979_dp, &
      5.0_dp, 7.071067811865476_dp, 10.0_dp, 20.0_dp, &
      30.0_dp, 40.0_dp, 50.0_dp, 100.0_dp, 300.0_dp, 500.0_dp ]

   character(len=256) :: env
   integer :: env_len, env_stat, it, ip, ipl, ir
   real(dp) :: ta(6), pressure_pa
   real(dp), allocatable :: kinf(:)

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

   if (n_plog_reactions <= 0 .or. n_plog_nodes <= 0) then
      write(*,'(a)') 'ERROR: loaded mechanism contains no PLOG reactions'
      error stop 1
   endif

   allocate(kinf(nr))
   write(*,'(a,i0,a,i0)') '# REAL_PLOG_RATES reactions=',             &
      n_plog_reactions, ' terms=', n_plog_nodes
   write(*,'(a)') 'kind,packed_index,reaction_index,T_K,P_Pa,k_ck'

   do it = 1, size(temperatures)
      ta(1) = temperatures(it)
      ta(2) = ta(1)*ta(1)
      ta(3) = ta(2)*ta(1)
      ta(4) = ta(3)*ta(1)
      ta(5) = 1.0_dp/ta(1)
      ta(6) = log(ta(1))

      do ip = 1, size(pressures_atm)
         pressure_pa = pressures_atm(ip)*atm_to_pa
         kinf = -1.0_dp
         call plog_kinf_eval(ta, pressure_pa, kinf)

         do ipl = 1, n_plog_reactions
            ir = plog_reaction(ipl)
            if (.not. ieee_is_finite(kinf(ir)) .or. kinf(ir) <= 0.0_dp) then
               write(*,'(a,i0,a,es24.16,a,es24.16)')                  &
                  'ERROR: non-finite/non-positive PLOG k at reaction ', &
                  ir, ' T=', ta(1), ' P=', pressure_pa
               error stop 1
            endif
            write(*,'(a,",",i0,",",i0,3(",",es24.16))') 'RATE',    &
               ipl, ir, ta(1), pressure_pa, kinf(ir)
         enddo
      enddo
   enddo

end program driver_plog_real_rates
