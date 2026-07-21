!   ***********************************************************************************************
!   **  SpeedCHEM smoke-test driver                                                              **
!   **                                                                                           **
!   **  Exercises the full library pipeline against a real CHEMKIN mechanism:                    **
!   **                                                                                           **
!   **     chem.inp + therm.dat  --[CKINTP]-->  cklink  --[SCcklink/setup]-->  internal state    **
!   **                           --[chemistry_ODE_integrate]-->  time-integrated solution        **
!   **                                                                                           **
!   **  It is deliberately minimal: it sets `mechdir`, runs the standard `chemistry_input`       **
!   **  orchestrator (which interprets the mechanism and builds the sparse chemistry), seeds a   **
!   **  plausible constant-volume initial state (temperature + species molar concentrations),    **
!   **  integrates a short interval, and asserts the result is finite.                           **
!   **                                                                                           **
!   **  Exit status: 0 = pass, non-zero = failure (so CI / run_tests.sh can gate on it).         **
!   **                                                                                           **
!   **  The mechanism directory is taken from the SC_MECHDIR environment variable (must end      **
!   **  with a path separator, e.g. "/mnt/.../gttdata/"), falling back to a compiled-in default. **
!   ***********************************************************************************************

program driver_smoke

    use working_precision,   only: dp
    use chemistry_setup,     only: mechdir, solver, use_speedchem
    use speedchem,           only: ns, nr, neq
    use ieee_arithmetic,     only: ieee_is_finite

    implicit none

    ! --- Configuration -----------------------------------------------------------------------
    character(len=256) :: mech
    integer            :: env_len, env_stat

    ! --- Problem state -----------------------------------------------------------------------
    real (dp), allocatable :: yin(:), atol(:)
    real (dp)              :: rtol, t0, tf
    real (dp)              :: T_init, conc_init
    integer                :: i, n_nonfinite

    ! --- Exit bookkeeping --------------------------------------------------------------------
    integer, parameter :: EXIT_OK = 0, EXIT_FAIL = 1
    character(len=*), parameter :: bar = &
        '----------------------------------------------------------------------'

    write(*,'(a)') bar
    write(*,'(a)') ' SpeedCHEM smoke test'
    write(*,'(a)') bar

    ! -----------------------------------------------------------------------------------------
    ! 1. Resolve the mechanism directory.
    !    `mechdir` is concatenated with bare filenames ("chem.inp", "therm.dat", "cklink"),
    !    so it MUST end with a path separator.
    ! -----------------------------------------------------------------------------------------
    call get_environment_variable('SC_MECHDIR', mech, env_len, env_stat)
    if (env_stat /= 0 .or. env_len == 0) then
        mech = '/mnt/sn850/Projects/gtt/input/gttdata/'
        write(*,'(a)') ' SC_MECHDIR not set; using default:'
    else
        write(*,'(a)') ' Using SC_MECHDIR:'
    end if
    if (mech(len_trim(mech):len_trim(mech)) /= '/') mech = trim(mech)//'/'
    write(*,'(3a)') '   ', trim(mech)

    mechdir       = mech
    use_speedchem = .true.

    ! -----------------------------------------------------------------------------------------
    ! 2. Run the standard setup orchestrator. This interprets chem.inp/therm.dat into the
    !    intermediate `cklink` (via CKINTP), links it into SpeedCHEM's internal structures,
    !    builds the sparse Jacobian sparsity pattern, and initialises the ODE solver workspace.
    !    On return, ns/nr/neq describe the loaded mechanism.
    ! -----------------------------------------------------------------------------------------
    write(*,'(a)') bar
    write(*,'(a)') ' Interpreting mechanism and building sparse chemistry...'
    call chemistry_input

    write(*,'(a,i0)')  '   species  (ns)  = ', ns
    write(*,'(a,i0)')  '   reactions(nr)  = ', nr
    write(*,'(a,i0)')  '   ODE size (neq) = ', neq
    write(*,'(2a)')    '   solver         = ', trim(solver)

    if (ns <= 0 .or. neq /= ns + 1) then
        write(*,'(a)') ' FAIL: mechanism setup produced an implausible size (expected neq = ns+1).'
        call finish(EXIT_FAIL)
    end if

    ! -----------------------------------------------------------------------------------------
    ! 3. Seed a constant-volume initial state.
    !    State-vector layout for this build: yin(1) = temperature [K],
    !    yin(2:neq) = species molar concentrations [mol/m3].
    !    We use a mild, uniform seed: the production mechanism is a Zero-RK metadata shell whose
    !    reaction rates are ~1e-30, so nothing should actually react. A converging, finite
    !    integration is exactly the pass condition we want for a smoke test.
    ! -----------------------------------------------------------------------------------------
    allocate(yin(neq), atol(neq))

    T_init    = 1000.0_dp          ! K
    conc_init = 1.0e-3_dp          ! mol/m3, small uniform seed

    yin(1)      = T_init
    yin(2:neq)  = conc_init

    rtol      = 1.0e-6_dp
    atol      = 1.0e-12_dp
    t0        = 0.0_dp
    tf        = 1.0e-6_dp          ! s

    write(*,'(a)') bar
    write(*,'(a)')        ' Integrating constant-volume chemistry...'
    write(*,'(a,1pe10.3,a,1pe10.3,a)') '   t: ', t0, ' -> ', tf, ' s'
    write(*,'(a,0pf8.2,a)')            '   T_init = ', T_init, ' K'

    call chemistry_ODE_integrate(neq, rtol, atol, t0, tf, yin)

    ! -----------------------------------------------------------------------------------------
    ! 4. Assert the solution is finite (no NaN/Inf from a diverged or corrupted solve).
    ! -----------------------------------------------------------------------------------------
    n_nonfinite = 0
    do i = 1, neq
        if (.not. ieee_is_finite(yin(i))) n_nonfinite = n_nonfinite + 1
    end do

    write(*,'(a)') bar
    write(*,'(a,0pf10.3,a)') '   T_final = ', yin(1), ' K'
    write(*,'(a,i0,a,i0)')   '   non-finite entries: ', n_nonfinite, ' / ', neq

    if (n_nonfinite > 0) then
        write(*,'(a)') ' FAIL: integration produced non-finite values.'
        call finish(EXIT_FAIL)
    end if

    if (yin(1) <= 0.0_dp) then
        write(*,'(a)') ' FAIL: final temperature is non-physical (<= 0 K).'
        call finish(EXIT_FAIL)
    end if

    write(*,'(a)') bar
    write(*,'(a)') ' PASS: pipeline ran and produced a finite constant-volume solution.'
    write(*,'(a)') bar
    call finish(EXIT_OK)

contains

    subroutine finish(code)
        integer, intent(in) :: code
        if (allocated(yin))  deallocate(yin)
        if (allocated(atol)) deallocate(atol)
        call exit(code)
    end subroutine finish

end program driver_smoke
