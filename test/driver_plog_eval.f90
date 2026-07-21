!   ***********************************************************************************************
!   **  SpeedCHEM PLOG rate-evaluation unit test  (stage 2)                                       **
!   **                                                                                           **
!   **  Exercises reacpar::plog_kinf_eval in isolation -- no mechanism file, no sparse setup,     **
!   **  no integration. It populates the packed PLOG arrays by hand for a single 3-node PLOG      **
!   **  reaction, then checks the interpolated forward rate constant against independently        **
!   **  hand-computed reference values at:                                                        **
!   **    - each pressure node exactly                                                            **
!   **    - the geometric mean of two adjacent nodes (mid-interval in ln P)                       **
!   **    - below the lowest node and above the highest node (nearest-endpoint clamp)             **
!   **    - several temperatures                                                                  **
!   **                                                                                           **
!   **  The reference k is computed here with the very same closed-form definition the plan       **
!   **  specifies (log-linear interpolation of ln k in ln P; per-node Arrhenius), but written     **
!   **  independently of the library routine, so a coding error in plog_kinf_eval that diverged   **
!   **  from the definition would be caught. Exit 0 = all checks within tolerance.                **
!   ***********************************************************************************************

program driver_plog_eval

    use working_precision, only: dp
    use reacpar,           only: n_plog_reactions, n_plog_nodes, n_plog_terms, &
                                 plog_reaction, plog_node_ptr,                 &
                                 plog_logP, plog_A, plog_b, plog_EoverR,       &
                                 plog_kinf_eval

    implicit none

    integer, parameter :: NR_TEST = 4          ! pretend there are 4 reactions
    integer, parameter :: IRX     = 3          ! the PLOG reaction is #3
    real (dp), parameter :: atm = 101325.0_dp

    ! Three PLOG nodes: P[atm], A, b, E/R[K]  (E/R chosen directly in K here)
    real (dp), parameter :: P1 = 0.1_dp,  A1 = 1.0e11_dp, B1 = 0.10_dp, ER1 = 150.0_dp
    real (dp), parameter :: P2 = 1.0_dp,  A2 = 1.0e12_dp, B2 = 0.00_dp, ER2 = 250.0_dp
    real (dp), parameter :: P3 = 10.0_dp, A3 = 1.0e13_dp, B3 =-0.10_dp, ER3 = 350.0_dp

    real (dp) :: kinf(NR_TEST)
    real (dp) :: Ta(6)
    integer   :: nfail
    real (dp), parameter :: rtol = 1.0e-12_dp

    nfail = 0

    ! ---- Populate the packed PLOG arrays for one 3-node reaction --------------
    n_plog_reactions = 1
    n_plog_nodes     = 3
    n_plog_terms     = 3
    allocate(plog_reaction(1), plog_node_ptr(0:1))
    allocate(plog_logP(3), plog_A(3), plog_b(3), plog_EoverR(3))
    plog_reaction(1) = IRX
    plog_node_ptr(0) = 0
    plog_node_ptr(1) = 3
    plog_logP   = [ log(P1*atm), log(P2*atm), log(P3*atm) ]
    plog_A      = [ A1, A2, A3 ]
    plog_b      = [ B1, B2, B3 ]
    plog_EoverR = [ ER1, ER2, ER3 ]

    ! ---- Checks ---------------------------------------------------------------
    ! At each node, k must equal that node's plain Arrhenius value.
    call check('node1 @ 900K',  900.0_dp,  P1*atm, ref_node(1, 900.0_dp))
    call check('node2 @ 900K',  900.0_dp,  P2*atm, ref_node(2, 900.0_dp))
    call check('node3 @ 900K',  900.0_dp,  P3*atm, ref_node(3, 900.0_dp))
    call check('node2 @ 1500K', 1500.0_dp, P2*atm, ref_node(2, 1500.0_dp))

    ! Geometric-mean pressure between nodes 1-2 and 2-3 (theta = 0.5 in ln P).
    call check('geo(1,2) @ 1200K', 1200.0_dp, sqrt(P1*P2)*atm, &
               ref_interp(1, 1200.0_dp, sqrt(P1*P2)*atm))
    call check('geo(2,3) @ 1200K', 1200.0_dp, sqrt(P2*P3)*atm, &
               ref_interp(2, 1200.0_dp, sqrt(P2*P3)*atm))

    ! An arbitrary intermediate pressure (not the midpoint) between nodes 2-3.
    call check('P=3atm @ 1000K', 1000.0_dp, 3.0_dp*atm, &
               ref_interp(2, 1000.0_dp, 3.0_dp*atm))

    ! Below lowest node -> node 1; above highest node -> node 3 (clamp).
    call check('below-min @ 800K', 800.0_dp, 0.01_dp*atm, ref_node(1, 800.0_dp))
    call check('above-max @ 800K', 800.0_dp, 100.0_dp*atm, ref_node(3, 800.0_dp))

    ! Just below / just above an interior node (node 2 = 1 atm).
    call check('just<node2 @ 1100K', 1100.0_dp, 0.999_dp*atm, &
               ref_interp(1, 1100.0_dp, 0.999_dp*atm))
    call check('just>node2 @ 1100K', 1100.0_dp, 1.001_dp*atm, &
               ref_interp(2, 1100.0_dp, 1.001_dp*atm))

    if (nfail == 0) then
        write(*,'(a)') 'RESULT: PASS - all PLOG rate checks within tolerance'
        call exit(0)
    else
        write(*,'(a,i0,a)') 'RESULT: FAIL - ', nfail, ' PLOG rate check(s) off'
        call exit(1)
    end if

contains

    ! Set Ta from T and call the library evaluator; return kinf(IRX).
    real (dp) function lib_k(T, P_pa) result(k)
        real (dp), intent(in) :: T, P_pa
        Ta(1) = T
        Ta(2) = T*T
        Ta(3) = Ta(2)*T
        Ta(4) = Ta(3)*T
        Ta(5) = 1.0_dp / T
        Ta(6) = log(T)
        kinf = -1.0_dp
        call plog_kinf_eval(Ta, P_pa, kinf)
        k = kinf(IRX)
    end function lib_k

    ! Reference single-node Arrhenius rate: A * exp(b*lnT - (E/R)/T).
    real (dp) function ref_node(node, T) result(k)
        integer,   intent(in) :: node
        real (dp), intent(in) :: T
        k = plog_A(node) * exp(plog_b(node)*log(T) - plog_EoverR(node)/T)
    end function ref_node

    ! Reference log-linear interpolation of ln k in ln P, between node j
    ! and j+1, at pressure P_pa. Independent restatement of the definition.
    real (dp) function ref_interp(j, T, P_pa) result(k)
        integer,   intent(in) :: j
        real (dp), intent(in) :: T, P_pa
        real (dp) :: theta, lnk
        theta = (log(P_pa) - plog_logP(j)) / (plog_logP(j+1) - plog_logP(j))
        lnk   = (1.0_dp - theta)*log(ref_node(j, T)) + theta*log(ref_node(j+1, T))
        k = exp(lnk)
    end function ref_interp

    subroutine check(label, T, P_pa, kref)
        character(len=*), intent(in) :: label
        real (dp),        intent(in) :: T, P_pa, kref
        real (dp) :: kgot, relerr
        kgot   = lib_k(T, P_pa)
        relerr = abs(kgot - kref) / max(abs(kref), tiny(1.0_dp))
        if (relerr <= rtol) then
            write(*,'(a,a,a,es20.12,a,es9.2)') '  ok   ', label, ' k=', kgot, &
                                               ' relerr=', relerr
        else
            write(*,'(a,a,a,es20.12,a,es20.12,a,es9.2)') '  FAIL ', label, &
                ' got=', kgot, ' ref=', kref, ' relerr=', relerr
            nfail = nfail + 1
        end if
    end subroutine check

end program driver_plog_eval
