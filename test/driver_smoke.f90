!   ***********************************************************************************************
!   **  SpeedCHEM ignition smoke-test driver                                                     **
!   **                                                                                           **
!   **  Exercises the full library pipeline against a real CHEMKIN mechanism and verifies that   **
!   **  a hot, near-stoichiometric fuel/air mixture actually ignites:                            **
!   **                                                                                           **
!   **     chem.inp + therm.dat  --[CKINTP]-->  cklink  --[SCcklink/setup]-->  internal state    **
!   **                           --[chemistry_ODE_integrate]-->  time-integrated solution        **
!   **                                                                                           **
!   **  Steps:                                                                                    **
!   **    1. Set `mechdir` (from SC_MECHDIR env, default = bundled test/data/).                   **
!   **    2. Run `chemistry_input` to interpret+link the mechanism and build sparse chemistry.    **
!   **    3. Seed a constant-volume state: high temperature + stoichiometric fuel/O2/N2, looked   **
!   **       up by species name so it is independent of internal ordering.                        **
!   **    4. Integrate long enough for autoignition.                                              **
!   **    5. Assert the solution is finite AND the temperature rose significantly (ignition).     **
!   **                                                                                           **
!   **  Exit status: 0 = pass, non-zero = failure (so CI / run_tests.sh can gate on it).         **
!   ***********************************************************************************************

program driver_smoke

    use working_precision,   only: dp
    use chemistry_setup,     only: mechdir, solver, use_speedchem
    use speedchem,           only: ns, nr, neq, specie
    use ieee_arithmetic,     only: ieee_is_finite

    implicit none

    ! --- Configuration -----------------------------------------------------------------------
    character(len=256) :: mech
    integer            :: env_len, env_stat
    ! Optional solver override (SC_SOLVER): pick a numeric-Jacobian solver for
    ! Optional solver override for targeted solver regression tests.
    character(len=15)  :: solver_override
    integer            :: sv_len, sv_stat

    ! --- Problem state -----------------------------------------------------------------------
    real (dp), allocatable :: yin(:), atol(:)
    real (dp)              :: rtol, t0, tf
    real (dp)              :: T_init, T_final, T_rise
    integer                :: i, n_nonfinite

    ! --- Mixture setup -----------------------------------------------------------------------
    ! Constant-volume autoignition of an n-hexadecane / air mixture.
    ! Global reaction (stoichiometric): C16H34 + 24.5 O2 -> 16 CO2 + 17 H2O.
    integer            :: i_fuel, i_o2, i_n2
    real (dp)          :: c_fuel, c_o2, c_n2
    real (dp), parameter :: phi        = 1.0_dp        ! equivalence ratio
    real (dp), parameter :: o2_stoich  = 24.5_dp       ! moles O2 per mole C16H34
    real (dp), parameter :: n2_per_o2  = 3.76_dp       ! air composition
    real (dp), parameter :: c_fuel_ref = 1.0_dp        ! mol/m3 fuel (scale is arbitrary here)
    real (dp), parameter :: T_ignite_threshold = 100.0_dp  ! K rise counted as "ignition"

    ! --- Exit bookkeeping --------------------------------------------------------------------
    integer, parameter :: EXIT_OK = 0, EXIT_FAIL = 1
    character(len=*), parameter :: bar = &
        '----------------------------------------------------------------------'

    write(*,'(a)') bar
    write(*,'(a)') ' SpeedCHEM ignition smoke test'
    write(*,'(a)') bar

    ! -----------------------------------------------------------------------------------------
    ! 1. Resolve the mechanism directory (must end with a path separator).
    ! -----------------------------------------------------------------------------------------
    call get_environment_variable('SC_MECHDIR', mech, env_len, env_stat)
    if (env_stat /= 0 .or. env_len == 0) then
        mech = 'test/data/'
        write(*,'(a)') ' SC_MECHDIR not set; using default:'
    else
        write(*,'(a)') ' Using SC_MECHDIR:'
    end if
    if (mech(len_trim(mech):len_trim(mech)) /= '/') mech = trim(mech)//'/'
    write(*,'(3a)') '   ', trim(mech)

    mechdir       = mech
    use_speedchem = .true.

    ! -----------------------------------------------------------------------------------------
    ! 2. Interpret + link the mechanism and build the sparse chemistry / solver workspace.
    ! -----------------------------------------------------------------------------------------
    write(*,'(a)') bar
    write(*,'(a)') ' Interpreting mechanism and building sparse chemistry...'
    call chemistry_input

    ! Optional: override the solver.
    call get_environment_variable('SC_SOLVER', solver_override, sv_len, sv_stat)
    if (sv_stat == 0 .and. sv_len > 0) then
        solver = trim(adjustl(solver_override))
        write(*,'(2a)') '   solver override: ', trim(solver)
    end if

    write(*,'(a,i0)')  '   species  (ns)  = ', ns
    write(*,'(a,i0)')  '   reactions(nr)  = ', nr
    write(*,'(a,i0)')  '   ODE size (neq) = ', neq
    write(*,'(2a)')    '   solver         = ', trim(solver)

    if (ns <= 0 .or. neq /= ns + 1) then
        write(*,'(a)') ' FAIL: mechanism setup produced an implausible size (expected neq = ns+1).'
        call finish(EXIT_FAIL)
    end if

    ! -----------------------------------------------------------------------------------------
    ! 3. Seed a hot, stoichiometric fuel/air state.
    !    State-vector layout for this build: yin(1) = temperature [K],
    !    yin(2:neq) = species molar concentrations [mol/m3]. Species indices are resolved by
    !    NAME so the seed is robust to internal reordering (permutation for sparse fill-in).
    ! -----------------------------------------------------------------------------------------
    i_fuel = species_index('C16H34')
    if (i_fuel == 0) i_fuel = species_index('IC16H34')   ! fall back to iso-cetane
    i_o2   = species_index('O2')
    i_n2   = species_index('N2')

    if (i_fuel == 0 .or. i_o2 == 0 .or. i_n2 == 0) then
        write(*,'(a)') ' FAIL: could not locate fuel / O2 / N2 in the mechanism species list.'
        write(*,'(a,i0,a,i0,a,i0)') '   (fuel=', i_fuel, ' O2=', i_o2, ' N2=', i_n2, ')'
        call finish(EXIT_FAIL)
    end if

    c_fuel = c_fuel_ref
    c_o2   = o2_stoich * c_fuel / phi
    c_n2   = n2_per_o2 * c_o2

    allocate(yin(neq), atol(neq))
    yin        = 0.0_dp
    T_init     = 1400.0_dp                 ! K — hot enough for prompt autoignition
    yin(1)     = T_init
    yin(1 + i_fuel) = c_fuel               ! +1: index 1 is temperature
    yin(1 + i_o2)   = c_o2
    yin(1 + i_n2)   = c_n2

    rtol      = 1.0e-6_dp
    atol      = 1.0e-15_dp
    t0        = 0.0_dp
    tf        = 1.0e-2_dp                   ! s — long enough to capture ignition

    write(*,'(a)') bar
    write(*,'(a)')            ' Integrating constant-volume autoignition...'
    write(*,'(a,a,a,i0,a)')   '   fuel  : ', trim(specie(i_fuel)), ' (index ', i_fuel, ')'
    write(*,'(a,1pe10.3,a,1pe10.3,a)') '   t     : ', t0, ' -> ', tf, ' s'
    write(*,'(a,0pf8.2,a)')            '   T_init: ', T_init, ' K'

    call chemistry_ODE_integrate(neq, rtol, atol, t0, tf, yin)

    ! -----------------------------------------------------------------------------------------
    ! 4. Checks: solution finite, temperature physical, and ignition occurred.
    ! -----------------------------------------------------------------------------------------
    n_nonfinite = 0
    do i = 1, neq
        if (.not. ieee_is_finite(yin(i))) n_nonfinite = n_nonfinite + 1
    end do

    T_final = yin(1)
    T_rise  = T_final - T_init

    write(*,'(a)') bar
    write(*,'(a,0pf10.3,a)') '   T_final = ', T_final, ' K'
    write(*,'(a,0pf10.3,a)') '   T_rise  = ', T_rise,  ' K'
    write(*,'(a,i0,a,i0)')   '   non-finite entries: ', n_nonfinite, ' / ', neq

    if (n_nonfinite > 0) then
        write(*,'(a)') ' FAIL: integration produced non-finite values.'
        call finish(EXIT_FAIL)
    end if

    if (T_final <= 0.0_dp) then
        write(*,'(a)') ' FAIL: final temperature is non-physical (<= 0 K).'
        call finish(EXIT_FAIL)
    end if

    if (T_rise < T_ignite_threshold) then
        write(*,'(a,0pf6.1,a)') ' FAIL: no ignition — temperature rise below threshold (', &
                                 T_ignite_threshold, ' K).'
        call finish(EXIT_FAIL)
    end if

    write(*,'(a)') bar
    write(*,'(a,0pf7.1,a)') ' PASS: mixture ignited (T rose ', T_rise, ' K) with a finite solution.'
    write(*,'(a)') bar
    call finish(EXIT_OK)

contains

    ! Case-sensitive lookup of a species by name; returns its 1-based index in the species
    ! list, or 0 if not present. `specie` is populated by SCcklink after linking.
    integer function species_index(name) result(idx)
        character(len=*), intent(in) :: name
        integer :: k
        idx = 0
        do k = 1, ns
            if (trim(adjustl(specie(k))) == name) then
                idx = k
                return
            end if
        end do
    end function species_index

    subroutine finish(code)
        integer, intent(in) :: code
        if (allocated(yin))  deallocate(yin)
        if (allocated(atol)) deallocate(atol)
        call exit(code)
    end subroutine finish

end program driver_smoke
