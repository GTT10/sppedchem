!   ***********************************************************************************************
!   **  SpeedCHEM RHS probe driver  (stage 2 PLOG integration check)                              **
!   **                                                                                           **
!   **  Links a mechanism, seeds ONE constant-volume state (hot stoichiometric n-hexadecane/air) **
!   **  and evaluates the RHS derivative vector d(state)/dt exactly once via SC_conV, then prints **
!   **  it as canonical text. It does NOT integrate, so it is independent of the stiff-solver     **
!   **  machinery.                                                                                **
!   **                                                                                           **
!   **  Used to prove that a PLOG mechanism whose PLOG nodes reproduce a reaction's original      **
!   **  Arrhenius rate yields the SAME RHS as the plain-Arrhenius mechanism: run this driver on   **
!   **  test/data/ (plain PRF) and test/data_plog_prf/ (one reaction PLOG-ified with identical    **
!   **  nodes) and diff the two dumps -- they must match, which exercises plog_kinf_eval inside   **
!   **  the real RHS (mass_action) that every solver calls.                                       **
!   ***********************************************************************************************

program driver_rhs_probe

    use working_precision, only: dp
    use chemistry_setup,   only: mechdir, use_speedchem
    use speedchem,         only: ns, nr, neq, specie
    use speedchem_conV,    only: SC_conV

    implicit none

    character(len=256) :: mech
    integer            :: env_len, env_stat, i, i_fuel, i_o2, i_n2
    real (dp), allocatable :: yin(:), dyindt(:)
    real (dp), parameter :: o2_stoich = 24.5_dp, n2_per_o2 = 3.76_dp

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
        call exit(2)
    end if

    allocate(yin(neq), dyindt(neq))
    yin = 0.0_dp
    yin(1)          = 1400.0_dp
    yin(1 + i_fuel) = 1.0_dp
    yin(1 + i_o2)   = o2_stoich
    yin(1 + i_n2)   = n2_per_o2 * o2_stoich

    ! Single RHS evaluation.
    call SC_conV(neq, 0.0_dp, yin, dyindt)

    ! Canonical dump: ns, nr, then every derivative component. Fixed
    ! formatting so two runs are byte-comparable.
    write(*,'(a,i0,a,i0)') '# RHS probe ns=', ns, ' nr=', nr
    do i = 1, neq
        write(*,'(i0,1x,es23.15e3)') i, dyindt(i)
    end do

    call exit(0)

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
