!   ***********************************************************************************************
!   **  SpeedCHEM PLOG parse-only driver  (cklink v2, stage 1)                                    **
!   **                                                                                           **
!   **  Exercises the PLOG data-plumbing path WITHOUT integrating:                               **
!   **                                                                                           **
!   **     chem.inp + therm.dat --[CKINTP]--> cklink(v2) --[SCcklink]--> reacpar PLOG arrays      **
!   **                          --[plog_dump_canonical]--> stable text on stdout                 **
!   **                                                                                           **
!   **  It stops right after linking (never builds the ODE workspace or integrates), so a PLOG   **
!   **  mechanism can be parsed and round-tripped even though PLOG rate evaluation is not yet     **
!   **  implemented (stage 2). SCcklink only skips its PLOG fail-closed guard when SC_PLOG_DUMP   **
!   **  is set in the environment; this driver sets it defensively too.                          **
!   **                                                                                           **
!   **  Usage:  SC_MECHDIR=test/data_plog/ SC_PLOG_DUMP=1 driver_plog                            **
!   **  Output: the canonical PLOG dump (see reacpar::plog_dump_canonical); exit 0 on success.   **
!   ***********************************************************************************************

program driver_plog

    use working_precision,     only: dp
    use chemistry_setup,       only: mechdir, use_speedchem
    use chemkinii_interpreter                       ! ckintp
    use reacpar,               only: n_plog_reactions, n_plog_nodes, &
                                     plog_dump_canonical
    use speedchem,             only: ns, nr

    implicit none

    character(len=256) :: mech
    integer            :: env_len, env_stat

    ! ------------------------------------------------------------------
    ! 1. Mechanism directory (must end with a path separator).
    ! ------------------------------------------------------------------
    call get_environment_variable('SC_MECHDIR', mech, env_len, env_stat)
    if (env_stat /= 0 .or. env_len == 0) then
        mech = 'test/data_plog/'
    end if
    if (mech(len_trim(mech):len_trim(mech)) /= '/') mech = trim(mech)//'/'

    mechdir       = mech
    use_speedchem = .true.

    ! ------------------------------------------------------------------
    ! 2. Interpret + link only (no ODE setup, no integration).
    !    CKINTP writes cklink v2 (with the PLOG section); SCcklink reads
    !    it back into the reacpar packed arrays. With SC_PLOG_DUMP set,
    !    SCcklink does not fail-closed on PLOG presence.
    ! ------------------------------------------------------------------
    call ckintp
    call SCcklink

    write(*,'(a,i0,a,i0,a)') '# linked: ', ns, ' species, ', nr, ' reactions'

    ! ------------------------------------------------------------------
    ! 3. Canonical PLOG dump to stdout (stable, diff-able).
    ! ------------------------------------------------------------------
    call plog_dump_canonical(6)

    call exit(0)

end program driver_plog
