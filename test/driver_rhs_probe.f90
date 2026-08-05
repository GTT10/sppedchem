!   ***********************************************************************************************
!   **  SpeedCHEM RHS probe driver  (stage 2 PLOG integration check)                              **
!   **                                                                                           **
!   **  Links a mechanism, seeds ONE physical constant-volume state (hot stoichiometric          **
!   **  n-hexadecane/air), evaluates d(state)/dt exactly once via SC_conV, and prints the          **
!   **  derivative as canonical text. It does not integrate, so it remains independent of         **
!   **  the stiff-solver machinery.                                                               **
!   **                                                                                           **
!   **  Used to prove that a PLOG mechanism whose nodes reproduce a reaction's original           **
!   **  Arrhenius rate yields the SAME RHS as the plain-Arrhenius mechanism. The input mole       **
!   **  fractions are converted to the mass-fraction state required by SC_conV, and density is    **
!   **  initialized from the requested T/P before the RHS call.                                  **
!   ***********************************************************************************************

program driver_rhs_probe

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use working_precision, only: dp
    use chemistry_setup,   only: mechdir, use_speedchem
    use speedchem,         only: ns, nr, neq, specie, molefr_to_massfr
    use speedchem_conV,    only: SC_conV
    use SCmixturethermo,   only: SCP, SCrho, rhoY, molar_volumes,       &
                                 pressurerhoT

    implicit none

    character(len=256) :: mech
    integer            :: env_len, env_stat, i, i_fuel, i_o2, i_n2
    real (dp), allocatable :: yin(:), dyindt(:), x(:)
    real (dp), parameter :: temperature0 = 1400.0_dp
    real (dp), parameter :: pressure0    = 2.0e6_dp
    real (dp), parameter :: o2_stoich    = 24.5_dp
    real (dp), parameter :: n2_per_o2    = 3.76_dp
    real (dp) :: reconstructed_pressure

    call get_environment_variable('SC_MECHDIR', mech, env_len, env_stat)
    if (env_stat /= 0 .or. env_len == 0) mech = 'test/data/'
    if (mech(len_trim(mech):len_trim(mech)) /= '/') mech = trim(mech)//'/'
    mechdir       = mech
    use_speedchem = .true.

    call chemistry_input

    i_fuel = species_index('C16H34')
    if (i_fuel == 0) i_fuel = species_index('IC16H34')
    i_o2   = species_index('O2')
    i_n2   = species_index('N2')
    if (i_fuel == 0 .or. i_o2 == 0 .or. i_n2 == 0) then
        write(*,'(a)') '# FAIL: fuel/O2/N2 not found'
        error stop 2
    end if

    allocate(yin(neq), dyindt(neq), x(ns))
    x = 0.0_dp
    x(i_fuel) = 1.0_dp
    x(i_o2)   = o2_stoich
    x(i_n2)   = n2_per_o2 * o2_stoich
    x = x / sum(x)

    yin = 0.0_dp
    yin(1) = temperature0
    call molefr_to_massfr(x, yin(2:neq))

    if (abs(sum(yin(2:neq)) - 1.0_dp) > 1.0e-13_dp) then
        write(*,'(a,es23.15e3)') '# FAIL: mass fractions do not sum to one: ', &
                                  sum(yin(2:neq))
        error stop 2
    endif

    ! Establish the constant-volume density from the requested initial
    ! pressure, temperature, and composition. SC_conV then recomputes pressure
    ! from this fixed density and its input state on every evaluation.
    SCP = pressure0
    call rhoY(yin(2:neq), yin(1))
    call molar_volumes
    reconstructed_pressure = pressurerhoT(yin(1), yin(2:neq))
    if (.not. ieee_is_finite(SCrho) .or. SCrho <= 0.0_dp .or.          &
        abs(reconstructed_pressure/pressure0 - 1.0_dp) > 1.0e-12_dp) then
        write(*,'(a,2(1x,es23.15e3))')                                 &
            '# FAIL: invalid initialized rho/P:', SCrho, reconstructed_pressure
        error stop 2
    endif

    call SC_conV(neq, 0.0_dp, yin, dyindt)
    if (.not. all(ieee_is_finite(dyindt))) then
        write(*,'(a)') '# FAIL: RHS contains non-finite values'
        error stop 2
    endif

    ! Canonical dump: ns, nr, then every derivative component. Fixed
    ! formatting keeps the plain and PLOG-equivalent runs byte-comparable.
    write(*,'(a,i0,a,i0)') '# RHS probe ns=', ns, ' nr=', nr
    do i = 1, neq
        write(*,'(i0,1x,es23.15e3)') i, dyindt(i)
    end do

contains

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

end program driver_rhs_probe
