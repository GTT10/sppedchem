! Verify the complete constant-volume analytic Jacobian for a genuinely
! pressure-dependent PLOG rate against an independent central difference
! of the production RHS. A compact four-reaction mechanism keeps the PLOG
! contribution well-conditioned relative to the rest of the chemistry.
program driver_plog_jacobian

   use working_precision, only: dp
   use chemistry_setup,   only: mechdir, use_speedchem, permutate_species
   use speedchem,         only: ns, nr, neq, uMW, Ainf
   use kinetics_mod,      only: tab_kinf, tab_dkinfdT
   use reacpar,           only: n_plog_reactions, n_plog_nodes, plog_A
   use speedchem_conV,    only: SC_conV, constV_jac_sparse
   use sparse_chemistry,  only: JAC_sparse
   use sparse_algebra,    only: sparse_to_dense
   use SCmixturethermo,   only: SCrho, SCP, pressurerhoT, molar_volumes
   use universal_constants, only: R
   use ieee_arithmetic,   only: ieee_is_finite

   implicit none

   real(dp), allocatable :: y(:), yp(:), ym(:), fp(:), fm(:)
   real(dp), allocatable :: jac(:,:), jac_fd(:,:), jac_base(:,:),      &
                            jac_fd_base(:,:), plog_A_saved(:)
   real(dp) :: target_pressure, h, scale_floor, rel_l2, max_scaled
   real(dp) :: norm_ref, norm_diff
   integer :: i

   mechdir = 'test/data_plog/'
   use_speedchem = .true.
   call chemistry_input

   if (permutate_species) then
      write(*,'(a)') 'RESULT: FAIL - Jacobian test requires the default unpermuted state'
      error stop 1
   endif
   if (n_plog_reactions /= 1 .or. n_plog_nodes /= 3) then
      write(*,'(a,2(i0,1x))') 'RESULT: FAIL - unexpected PLOG topology ', &
                               n_plog_reactions, n_plog_nodes
      error stop 1
   endif

!  Isolate the PLOG reaction numerically. The mechanism deliberately also
!  contains ordinary Arrhenius rows to test packing, but those rows add
!  cancellation noise without testing the PLOG Jacobian.
   Ainf = 0.0_dp
   tab_kinf = 0.0_dp
   tab_dkinfdT = 0.0_dp

   allocate(y(neq), yp(neq), ym(neq), fp(neq), fm(neq))
   allocate(jac(neq,neq), jac_fd(neq,neq), jac_base(neq,neq),          &
            jac_fd_base(neq,neq), plog_A_saved(n_plog_nodes))

!  A smooth, strictly-positive composition exercises every independent
!  species column without crossing the Y=0 regularisation branch.
   y(1) = 1100.0_dp
   do i = 1, ns
      y(i+1) = real(1 + mod(17*i,23),dp)
   enddo
   y(2:neq) = y(2:neq) / sum(y(2:neq))

!  Put the state at 3 atm, strictly inside the 1--100 atm PLOG interval.
!  pressurerhoT = 1000*rho*R*T*sum(Y/W).
   target_pressure = 3.0_dp*101325.0_dp
   SCrho = target_pressure / (1000.0_dp*R*y(1)*sum(y(2:neq)*uMW))
   call molar_volumes

!  Capture a baseline first. Subtracting this pair isolates the PLOG
!  contribution from pre-existing approximations in the legacy Jacobian.
   plog_A_saved = plog_A
   call evaluate_pair(jac_base,jac_fd_base)

!  Change the PLOG pressure slope while preserving positive rates.  The
!  common scale makes this reaction dominate finite-difference roundoff;
!  the baseline subtraction still removes all non-PLOG contributions.
   plog_A = plog_A_saved * [0.4_dp, 1.0_dp, 4.0_dp]
   call evaluate_pair(jac,jac_fd)

   jac    = jac    - jac_base
   jac_fd = jac_fd - jac_fd_base

   if (.not. all(ieee_is_finite(jac)) .or.                         &
       .not. all(ieee_is_finite(jac_fd))) then
      write(*,'(a)') 'RESULT: FAIL - non-finite analytic or finite-difference Jacobian'
      error stop 1
   endif

   norm_ref  = sqrt(sum(jac_fd*jac_fd))
   norm_diff = sqrt(sum((jac-jac_fd)*(jac-jac_fd)))
   rel_l2 = norm_diff/max(norm_ref,tiny(1.0_dp))
   scale_floor = maxval(abs(jac_fd))*1.0e-10_dp
   max_scaled = maxval(abs(jac-jac_fd) /                              &
                       max(max(abs(jac),abs(jac_fd)),scale_floor))

   SCP = pressurerhoT(y(1),y(2:neq))
   write(*,'(a,es14.6)') 'PLOG_JAC pressure_Pa=', SCP
   write(*,'(a,es14.6)') 'PLOG_JAC rel_l2=', rel_l2
   write(*,'(a,es14.6)') 'PLOG_JAC max_scaled=', max_scaled
   write(*,'(a,es14.6)') 'PLOG_JAC species_rel_l2=',                 &
      sqrt(sum((jac(2:neq,:)-jac_fd(2:neq,:))**2))/                  &
      max(sqrt(sum(jac_fd(2:neq,:)**2)),tiny(1.0_dp))

!  The legacy analytic Jacobian and tabulated thermo are not expected to
!  match finite differences to machine epsilon. These gates are tight
!  enough to catch a missing s/T term or missing dense Y coupling.
   if (rel_l2 > 2.0e-3_dp .or. max_scaled > 5.0e-2_dp) then
      write(*,'(a)') 'RESULT: FAIL - analytic PLOG Jacobian disagrees with central difference'
      error stop 1
   endif

   plog_A = plog_A_saved
   write(*,'(a)') 'RESULT: PASS - analytic PLOG Jacobian matches central difference'

contains

   subroutine evaluate_pair(jac_analytic,jac_numeric)
      real(dp), intent(out) :: jac_analytic(neq,neq), jac_numeric(neq,neq)
      integer :: col

!     Analytic sparse Jacobian -> dense reference layout.
      call SC_conV(neq,0.0_dp,y,fp)
      call constV_jac_sparse(neq,0.0_dp,y)
      call sparse_to_dense(JAC_sparse,jac_analytic)

!     Independent central differences. Density is held fixed; each RHS
!     call recomputes pressure from its own perturbed T,Y state.
      do col = 1, neq
         yp = y
         ym = y
         if (col == 1) then
            h = max(1.0e-4_dp, abs(y(col))*1.0e-6_dp)
         else
            h = max(1.0e-9_dp, abs(y(col))*1.0e-5_dp)
         endif
         yp(col) = yp(col) + h
         ym(col) = ym(col) - h
         call SC_conV(neq,0.0_dp,yp,fp)
         call SC_conV(neq,0.0_dp,ym,fm)
         jac_numeric(:,col) = (fp-fm)/(2.0_dp*h)
      enddo
   end subroutine evaluate_pair

end program driver_plog_jacobian
