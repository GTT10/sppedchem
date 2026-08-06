program driver_plog_extreme

   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use working_precision, only: dp
   use reacpar, only: n_plog_reactions, n_plog_nodes, n_plog_terms, &
                      plog_reaction, plog_node_ptr, plog_term_ptr,   &
                      plog_logP, plog_A, plog_b, plog_EoverR,       &
                      plog_kinf_eval

   implicit none

   integer, parameter :: nr_test = 3
   integer, parameter :: irx = 2
   real(dp), parameter :: temperature = 1000.0_dp
   real(dp), parameter :: pressure = 101325.0_dp
   real(dp), parameter :: target1 = 1.0e160_dp
   real(dp), parameter :: target2 = 5.0e159_dp
   real(dp), parameter :: rtol = 2.0e-13_dp

   real(dp) :: ta(6), kinf(nr_test), dkinfdt(nr_test), slope(nr_test)
   real(dp) :: expected_k, expected_g, expected_dkdt
   real(dp) :: g1, g2, rel_k, rel_dkdt

   ! Each raw product A*T**b overflows in IEEE double precision, while the
   ! complete Arrhenius expression remains finite after exp(-E/RT). This
   ! catches regressions from logarithmic evaluation back to direct products.
   n_plog_reactions = 1
   n_plog_nodes = 2
   n_plog_terms = 2
   allocate(plog_reaction(1), plog_node_ptr(0:1), plog_term_ptr(0:2))
   allocate(plog_logP(2), plog_A(2), plog_b(2), plog_EoverR(2))

   plog_reaction = [irx]
   plog_node_ptr = [0, 2]
   plog_term_ptr = [0, 1, 2]
   plog_logP = log(pressure)
   plog_A = [1.0e300_dp, 5.0e299_dp]
   plog_b = [20.0_dp, 24.0_dp]
   plog_EoverR(1) = temperature *                              &
      (log(plog_A(1)) + plog_b(1)*log(temperature) - log(target1))
   plog_EoverR(2) = temperature *                              &
      (log(plog_A(2)) + plog_b(2)*log(temperature) - log(target2))

   ta(1) = temperature
   ta(2) = temperature*temperature
   ta(3) = ta(2)*temperature
   ta(4) = ta(3)*temperature
   ta(5) = 1.0_dp/temperature
   ta(6) = log(temperature)

   kinf = -7.0_dp
   dkinfdt = -9.0_dp
   slope = -11.0_dp
   call plog_kinf_eval(ta, pressure, kinf, dkinfdt, slope)

   expected_k = target1 + target2
   g1 = plog_b(1)/temperature + plog_EoverR(1)/(temperature*temperature)
   g2 = plog_b(2)/temperature + plog_EoverR(2)/(temperature*temperature)
   expected_g = (target1*g1 + target2*g2)/expected_k
   expected_dkdt = expected_k*expected_g

   rel_k = abs(kinf(irx)-expected_k)/expected_k
   rel_dkdt = abs(dkinfdt(irx)-expected_dkdt)/abs(expected_dkdt)

   if (.not. ieee_is_finite(kinf(irx)) .or.                         &
       .not. ieee_is_finite(dkinfdt(irx)) .or.                      &
       rel_k > rtol .or. rel_dkdt > rtol .or.                       &
       slope(irx) /= 0.0_dp .or.                                   &
       kinf(1) /= -7.0_dp .or. kinf(3) /= -7.0_dp .or.              &
       dkinfdt(1) /= -9.0_dp .or. dkinfdt(3) /= -9.0_dp .or.        &
       slope(1) /= -11.0_dp .or. slope(3) /= -11.0_dp) then
      write(*,'(a,4(1x,es23.15e3))')                                &
         'RESULT: FAIL - extreme grouped PLOG evaluation',          &
         kinf(irx), expected_k, rel_k, rel_dkdt
      error stop 1
   endif

   write(*,'(a,3(1x,es13.5))')                                      &
      'RESULT: PASS - extreme grouped PLOG remains finite',          &
      kinf(irx), rel_k, rel_dkdt

   deallocate(plog_reaction, plog_node_ptr, plog_term_ptr)
   deallocate(plog_logP, plog_A, plog_b, plog_EoverR)
   n_plog_reactions = 0
   n_plog_nodes = 0
   n_plog_terms = 0

end program driver_plog_extreme
