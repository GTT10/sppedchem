program driver_mechanism_reload

   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use working_precision, only: dp
   use chemistry_setup, only: mechdir, use_speedchem
   use speedchem, only: ns, nr, neq, A0
   use reacpar, only: n_plog_reactions, plog_reaction
   use sparse_chemistry, only: JAC_sparse
   use sparse_definitions, only: sparse_is_allocated => allocated
   use speedchem_conV, only: SC_conV

   implicit none

   character(len=256) :: mechanism_a, mechanism_b
   real(dp), allocatable :: rhs(:), rhs_a(:), integrated(:), integrated_a(:)
   real(dp) :: fingerprint_a, fingerprint_b

   call get_command_argument(1, mechanism_a)
   call get_command_argument(2, mechanism_b)
   if (len_trim(mechanism_a) == 0 .or. len_trim(mechanism_b) == 0) then
      write(*,'(a)') 'ERROR: usage: driver_mechanism_reload MECH_A MECH_B'
      error stop 1
   endif

   call normalize_dir(mechanism_a)
   call normalize_dir(mechanism_b)
   use_speedchem = .true.

!  A: bundled non-PLOG mechanism.
   mechdir = mechanism_a
   call chemistry_input
   call require(ns == 51 .and. nr == 153, 'unexpected mechanism A size')
   call require(n_plog_reactions == 0, 'mechanism A unexpectedly has PLOG')
   call evaluate_rhs(rhs)
   allocate(rhs_a(size(rhs)))
   rhs_a = rhs
   fingerprint_a = sum(abs(rhs))
   deallocate(rhs)
   call integrate_probe(integrated_a)
   call chemistry_finalize
   call require_finalized('after mechanism A')

!  B: compact PLOG mechanism with different dimensions and sparse pattern.
   mechdir = mechanism_b
   call chemistry_input
   call require(ns == 9 .and. nr == 4, 'unexpected mechanism B size')
   call require(n_plog_reactions == 1, 'mechanism B lost its PLOG data')
   call evaluate_rhs(rhs)
   fingerprint_b = sum(abs(rhs))
   call require(fingerprint_b > 0.0_dp, 'mechanism B RHS is empty')
   deallocate(rhs)
   call integrate_probe(integrated)
   deallocate(integrated)
   call chemistry_finalize
   call require_finalized('after mechanism B')

!  A again: this catches stale dimensions, PLOG tags, caches, and sparse data.
   mechdir = mechanism_a
   call chemistry_input
   call require(ns == 51 .and. nr == 153, 'mechanism A reload size mismatch')
   call require(n_plog_reactions == 0, 'PLOG state leaked into mechanism A')
   call evaluate_rhs(rhs)
   call require(size(rhs) == size(rhs_a), 'mechanism A RHS size changed')
   call require(all(rhs == rhs_a), 'mechanism A RHS changed after A-B-A reload')
   call require(sum(abs(rhs)) == fingerprint_a,                         &
                'mechanism A fingerprint changed after reload')
   call integrate_probe(integrated)
   call require(all(integrated == integrated_a),                        &
                'mechanism A integration changed after A-B-A reload')
   deallocate(integrated, integrated_a)
   deallocate(rhs, rhs_a)
   call chemistry_finalize
   call require_finalized('after final mechanism A')
   call chemistry_finalize
   call require_finalized('after repeated finalization')

   write(*,'(a,1pe13.5,a,1pe13.5)')                                    &
      'RESULT: PASS - A-B-A reload is exact; fingerprints A=',         &
      fingerprint_a, ' B=', fingerprint_b

contains

   subroutine normalize_dir(path)
      character(len=*), intent(inout) :: path
      if (path(len_trim(path):len_trim(path)) /= '/') path = trim(path)//'/'
   end subroutine normalize_dir

   subroutine evaluate_rhs(values)
      real(dp), allocatable, intent(out) :: values(:)
      real(dp), allocatable :: state(:)

      allocate(state(neq), values(neq))
      state(1) = 1000.0_dp
      state(2:neq) = 1.0_dp/real(ns,dp)
      call SC_conV(neq, 0.0_dp, state, values)
      call require(all(ieee_is_finite(values)), 'RHS contains non-finite data')
      deallocate(state)
   end subroutine evaluate_rhs

   subroutine integrate_probe(state)
      real(dp), allocatable, intent(out) :: state(:)
      real(dp), allocatable :: tolerance(:)
      real(dp) :: t0, tf

      allocate(state(neq), tolerance(neq))
      state(1) = 1000.0_dp
      state(2:neq) = 1.0_dp/real(ns,dp)
      tolerance = 1.0e-14_dp
      t0 = 0.0_dp
      tf = 1.0e-10_dp
      call chemistry_ODE_integrate(neq, 1.0e-6_dp, tolerance, t0, tf, state)
      call require(all(ieee_is_finite(state)),                          &
                   'integration contains non-finite data')
      call require(abs(t0-tf) <= epsilon(tf), 'integration did not reach tf')
      deallocate(tolerance)
   end subroutine integrate_probe

   subroutine require_finalized(where)
      character(len=*), intent(in) :: where
      call require(ns == 0 .and. nr == 0 .and. neq == 0,                &
                   trim(where)//': dimensions were not reset')
      call require(.not.allocated(A0),                                  &
                   trim(where)//': mechanism arrays remain allocated')
      call require(.not.allocated(plog_reaction),                       &
                   trim(where)//': PLOG arrays remain allocated')
      call require(.not.sparse_is_allocated(JAC_sparse),                &
                   trim(where)//': sparse Jacobian remains allocated')
   end subroutine require_finalized

   subroutine require(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not.condition) then
         write(*,'(a)') 'ERROR: '//trim(message)
         error stop 1
      endif
   end subroutine require

end program driver_mechanism_reload
