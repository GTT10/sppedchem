!   ***********************************************************************************************
!   **                                                                                           **
!   **         |\   -  -.   ./                                                                   **
!   **         | \./ \/ | ./ /     _____                     __________  __________  ___         **
!   **       __|         /  /     / ___/____  ___  ___  ____/ / ____/ / / / ____/  |/  /         **
!   **       \    .        /.-/   \__ \/ __ \/ _ \/ _ \/ __  / /   / /_/ / __/ / /|_/ /          **
!   **        \   |\.|\/|    /   ___/ / /_/ /  __/  __/ /_/ / /___/ __  / /___/ /  / /           **
!   **         \__\     /___/   /____/ .___/\___/\___/\__,_/\____/_/ /_/_____/_/  /_/            **
!   **                              /_/                                                          **
!   **                                                                                           **
!   **       SpeedCHEM - A fast, portable Fortran library for Chemical Kinetics problems         **
!   **                          http://www.federicoperini.info/speedchem                         **
!   **                 Copyright (C) 2010-2013 Federico Perini. Version 1.0                      **
!   **                                                                                           **
!   ***********************************************************************************************
!   **                                                                                           **
!   ** License                                                                                   **
!   ** This file is part of SpeedCHEM.                                                           **
!   **                                                                                           **
!   ** SpeedCHEM is free software: you can redistribute it and/or modify                         **
!   ** it under the terms of the GNU General Public License as                                   **
!   ** published by the Free Software Foundation, either version 3 of the                        **
!   ** License, or any later version.                                                            **
!   **                                                                                           **
!   ** SpeedCHEM is distributed in the hope that it will be useful,                              **
!   ** but WITHOUT ANY WARRANTY; without even the implied warranty of                            **
!   ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                             **
!   ** GNU General Public License for more details.                                              **
!   **                                                                                           **
!   ** You should have received a copy of the GNU General Public License                         **
!   ** along with SpeedCHEM. If not, see <http://www.gnu.org/licenses/>.                         **
!   **                                                                                           **
!   ** Please send comments, requests, bug reports to federico.perini@gmail.com                  **
!   **                                                                                           **
!   ***********************************************************************************************

!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **                    Variable Allocation                      **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: tuesday, 25/05/2010                          **
!     **                                                             **
!     *****************************************************************

subroutine SCallocate

   use speedchem
   use SCthermodata
   use troepar
   use reacpar
   use kinetics_mod

!      use cellcluster
   implicit none
   integer i,j

!     Variables depending on number of reactions only
   allocate(A0(nr),AREV(nr),AINF(nr),b0(nr),bREV(nr),bINF(nr),E0(nr),&
   &EREV(nr),EINF(nr))
   allocate(reversibile(nr),REV(nr))
   allocate(HIGH(nr),LOW(nr),TROE(nr),THREE(nr))
   allocate(aTROE(nr),T1TROE(nr),T2TROE(nr),T3TROE(nr))


!     Variables depending on number of species only
   allocate(SCMW(ns),uMW(ns))
   allocate(specie(ns))

!     Variables depending on number of elements only
   allocate(elementi(nel), elMW(nel), uelMW(nel))

!     Variables depending on number of reactions and species
!      allocate(!nudiff(nr,ns),
!     &         third_body(nr,ns),!reagenti(nr,ns),prodotti(nr,ns),
!     &         third_body_beta(ns,nr))!,
!     &                            third_body_betaT(nr,ns))
!      allocate(nudiffT(ns,nr))
!      allocate(nudiffmask(nr,ns))

!     Thermodynamic database
   allocate(tsw(ns),aL(ns),bL(ns),cL(ns),dL(ns),eL(ns),fL(ns),gL(ns))
   allocate(aH(ns),bH(ns),cH(ns),dH(ns),eH(ns),fH(ns),gH(ns))

!     Allocate arrays for temperature-dependent data in case
!     the "SAVE" attribute has been chosen for the Jacobian evaluation


!      if (analytical_jac .and. save_thermal_parameters) then

!         !$OMP PARALLEL
!         !$OMP CRITICAL

!        Reaction rate constants
!         allocate(save_k0(nr), save_kinf(nr))


!         save_k0   = 0.d0
!         save_kinf = 0.d0

!         !$OMP END CRITICAL
!         !$OMP END PARALLEL

!      endif



end subroutine SCallocate



!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **              Deallocate unuseful dense data                 **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: tuesday, 15/01/2012                          **
!     **   Usage, copying, distribution of the proprietary routines  **
!     **   or any modified versions of them are not allowed without  **
!     **   approval of the author.                                   **
!     **                                                             **
!     *****************************************************************

subroutine SCdeallocate

   use speedchem, only : third_body, third_body_beta

   implicit none


!     This statement has to follow the first call to the chemistry
!     jacobian routine, to allow the sparse jacobian structure be
!     evaluated and stored in sparse format
!      deallocate(nudiff)

!      deallocate(third_body, third_body_beta)

end subroutine SCdeallocate



!     *****************************************************************
!     **                                                             **
!     **              Finalize a loaded SpeedCHEM mechanism          **
!     **                                                             **
!     ** chemistry_input owns module-global mechanism, thermodynamic,**
!     ** kinetics, sparse-Jacobian, and solver work arrays.  This     **
!     ** routine releases that complete state so chemistry_input can  **
!     ** be called again for a different mechanism in the same       **
!     ** process.  SCdeallocate above has a different legacy purpose: **
!     ** it is called in the middle of chemistry_input and therefore  **
!     ** must remain limited to disposable import storage.            **
!     **                                                             **
!     *****************************************************************

subroutine chemistry_finalize

   use speedchem, only: nel, ns, nr, neq, uneq, species, reactions,    &
      species_permutations, species_inverse_permutations, A0, AREV,    &
      AINF, ARi, b0, bREV, bINF, bRi, E0, EREV, EINF, ERi, SCMW, uMW, &
      elMW, uelMW, reversibile, REV, inotrev, nnotrev, duplicate,      &
      stoich_r, stoich_p, nudiff, nudiffT, third_body,                 &
      third_body_beta, stoich_r_pack, stoich_p_pack, nudiff_pack,     &
      i_stoich_r, i_stoich_p, i_nudiff, n_stoich_r, n_stoich_p,       &
      n_nudiff, sumnudiff, isumnudiff, HIGH, LOW, TROE, THREE,         &
      aTROE, T1TROE, T2TROE, T3TROE, iTHREE, iTB, nTHREE, nTB,        &
      specie, elementi, sparse_jac, njac, ljac, rowjac, coljac, iO2
   use sparse_chemistry, only: nudiff_sparse, nudiffT_sparse,          &
      nudiff_EQREV_sp, stoich_r_sp, stoich_r_eff_sp, stoich_p_sp,     &
      stoich_p_eff_sp, inudiff_sparse, inudiffT_sparse,               &
      istoich_r_sp, istoich_r_eff_sp, istoich_p_sp,                   &
      istoich_p_eff_sp, EM_sp, JAC_sparse, JACT_sparse, JACYY_sparse, &
      JACYYT_sparse, nudiffT_molarv_sparse, dq_dY_sparse,              &
      dq_dY_T_sparse, third_body_sp, tb_beta_sp, tbb_uMW_sp
   use sparse_definitions, only: sparse_reset => deallocate
   use find_mod, only: indices, i2D1, i2D2
   use SCthermodata, only: tsw, aL, bL, cL, dL, eL, fL, gL, aH, bH,   &
      cH, dH, eH, fH, gH, gibbsL, gibbsH, tab_CpsuR, tab_HsuRT,       &
      tab_SsuR, tab_dHdt, tab_dUdt, tab_dGdt, tab_dDG0dt, tab_dCvdT,  &
      tab_uKp, tab_uKc, tab_dKcdT, nels, nel1, nel2, nel3, nel4,      &
      el1, el2, el3, el4, EM, Emfr, iC, iO, iH, iN, iAr
   use SCspeciesthermo, only: Cpmol, Cvmol, Hmol, Smol, cpm, cvm, hm
   use troepar, only: ntbALL, itbALL, ltbALL, ntbFALL, itbFALL,       &
      ltbFALL, ntbTROE, itbTROE, ltbTROE, nTROE4, iTROE4, lTROE4,    &
      ntbLIND, itbLIND, ltbLIND, ntbSIMP, itbSIMP, ltbSIMP,           &
      iTROEitbALL, iLINDitbALL, iSIMPitbALL, iFALLitbALL,             &
      iTROEitbFALL, iLINDitbFALL, todotroe, zeroT2, aT2, uT1T2,      &
      T2T2, uT3T2, tab_troefactor, tab_log10Fcent
   use reacpar, only: nTREV, iTREV, nXREV, iXREV, nEQREV, iEQREV,     &
      tb_beta_pack, is_beta_pack, n_tb_beta, Arrhreac, Lindreac,      &
      Troereac, Revreac, rate_form, n_plog_reactions, n_plog_nodes,   &
      n_plog_terms, plog_reaction, plog_node_ptr, plog_logP,          &
      plog_term_ptr, plog_A, plog_b, plog_EoverR, store_uKc,          &
      store_uKc_T
   use kinetics_mod, only: ir, jr, ijr, ip, jp, ijp, i1r, i1p, i2r,  &
      i2p, i2D1r, i2D2r, i2D1p, i2D2p, unitr, unitp, indice_r,       &
      indice_p, v_stoich_r, v_stoich_p, iv_stoich_r, iv_stoich_p,    &
      v_stoich_r1, v_stoich_r2, v_stoich_ro, v_stoich_p1,            &
      v_stoich_p2, v_stoich_po, num_vro, num_vpo, num_vr2, num_vp2, &
      tab_k0, tab_kinf, tab_dkinfdt, tab_dk0dt, tab_Xkb,              &
      tab_Xdkbdt, tab_kinfT, iA0, save_k0, save_kinf, lsavek,         &
      store_dwdt, store_dwdt_T, store_kinf, store_kinf_T, store_qf,  &
      store_qb, store_C, store_valid
   use ode_solver, only: atol, rwork, iwork, iper1, iiper1, iper2,    &
      iiper2, iLU1, iLU2, rLU1, rLU2, JACT_VF90, R5_sys1, R5_sys2,   &
      DASPK_sys, vf90_opts, ncJAC, ncCONV, nsteps, nLUdec, nNewton,   &
      lrw, liw
   use plog_collect, only: plog_reset
   use chemkin_kiva, only: intwork, reawork, chawork, ICK, RCK, CCK,  &
      CKmw, ckrwork, ckiwork, ckne, ckns, cknr, ckneq

   implicit none

!  Shared mechanism arrays.
   if (allocated(species)) deallocate(species)
   if (allocated(reactions)) deallocate(reactions)
   if (allocated(species_permutations)) deallocate(species_permutations)
   if (allocated(species_inverse_permutations))                        &
      deallocate(species_inverse_permutations)

   if (allocated(A0)) deallocate(A0, AREV, AINF, b0, bREV, bINF, E0,  &
                                  EREV, EINF)
   if (allocated(ARi)) deallocate(ARi)
   if (allocated(bRi)) deallocate(bRi)
   if (allocated(ERi)) deallocate(ERi)
   if (allocated(SCMW)) deallocate(SCMW, uMW)
   if (allocated(elMW)) deallocate(elMW, uelMW)
   if (allocated(reversibile)) deallocate(reversibile, REV)
   if (allocated(inotrev)) deallocate(inotrev)
   if (allocated(duplicate)) deallocate(duplicate)

   if (allocated(stoich_r)) deallocate(stoich_r)
   if (allocated(stoich_p)) deallocate(stoich_p)
   if (allocated(nudiff)) deallocate(nudiff)
   if (allocated(nudiffT)) deallocate(nudiffT)
   if (allocated(third_body)) deallocate(third_body)
   if (allocated(third_body_beta)) deallocate(third_body_beta)
   if (allocated(stoich_r_pack)) deallocate(stoich_r_pack)
   if (allocated(stoich_p_pack)) deallocate(stoich_p_pack)
   if (allocated(nudiff_pack)) deallocate(nudiff_pack)
   if (allocated(i_stoich_r)) deallocate(i_stoich_r)
   if (allocated(i_stoich_p)) deallocate(i_stoich_p)
   if (allocated(i_nudiff)) deallocate(i_nudiff)
   if (allocated(n_stoich_r)) deallocate(n_stoich_r)
   if (allocated(n_stoich_p)) deallocate(n_stoich_p)
   if (allocated(n_nudiff)) deallocate(n_nudiff)
   if (allocated(sumnudiff)) deallocate(sumnudiff)
   if (allocated(isumnudiff)) deallocate(isumnudiff)

   if (allocated(HIGH)) deallocate(HIGH, LOW, TROE, THREE)
   if (allocated(aTROE)) deallocate(aTROE, T1TROE, T2TROE, T3TROE)
   if (allocated(iTHREE)) deallocate(iTHREE)
   if (allocated(iTB)) deallocate(iTB)
   if (allocated(specie)) deallocate(specie)
   if (allocated(elementi)) deallocate(elementi)

!  Sparse shared storage.
   call sparse_reset(nudiff_sparse)
   call sparse_reset(nudiffT_sparse)
   call sparse_reset(nudiff_EQREV_sp)
   call sparse_reset(stoich_r_sp)
   call sparse_reset(stoich_r_eff_sp)
   call sparse_reset(stoich_p_sp)
   call sparse_reset(stoich_p_eff_sp)
   call sparse_reset(inudiff_sparse)
   call sparse_reset(inudiffT_sparse)
   call sparse_reset(istoich_r_sp)
   call sparse_reset(istoich_r_eff_sp)
   call sparse_reset(istoich_p_sp)
   call sparse_reset(istoich_p_eff_sp)
   call sparse_reset(EM_sp)
   call sparse_reset(third_body_sp)
   call sparse_reset(tb_beta_sp)
   call sparse_reset(tbb_uMW_sp)

!  Thermodynamic data.
   if (allocated(tsw)) deallocate(tsw, aL, bL, cL, dL, eL, fL, gL,  &
                                  aH, bH, cH, dH, eH, fH, gH)
   if (allocated(gibbsL)) deallocate(gibbsL, gibbsH)
   if (allocated(tab_CpsuR)) deallocate(tab_CpsuR)
   if (allocated(tab_HsuRT)) deallocate(tab_HsuRT)
   if (allocated(tab_SsuR)) deallocate(tab_SsuR)
   if (allocated(tab_dHdt)) deallocate(tab_dHdt)
   if (allocated(tab_dUdt)) deallocate(tab_dUdt)
   if (allocated(tab_dGdt)) deallocate(tab_dGdt)
   if (allocated(tab_dDG0dt)) deallocate(tab_dDG0dt)
   if (allocated(tab_dCvdT)) deallocate(tab_dCvdT)
   if (allocated(tab_uKp)) deallocate(tab_uKp)
   if (allocated(tab_uKc)) deallocate(tab_uKc)
   if (allocated(tab_dKcdT)) deallocate(tab_dKcdT)
   if (allocated(nels)) deallocate(nels)
   if (allocated(nel1)) deallocate(nel1)
   if (allocated(nel2)) deallocate(nel2)
   if (allocated(nel3)) deallocate(nel3)
   if (allocated(nel4)) deallocate(nel4)
   if (allocated(el1)) deallocate(el1)
   if (allocated(el2)) deallocate(el2)
   if (allocated(el3)) deallocate(el3)
   if (allocated(el4)) deallocate(el4)
   if (allocated(EM)) deallocate(EM)
   if (allocated(Emfr)) deallocate(Emfr)
   if (allocated(Cpmol)) deallocate(Cpmol)
   if (allocated(Cvmol)) deallocate(Cvmol)
   if (allocated(Hmol)) deallocate(Hmol)
   if (allocated(Smol)) deallocate(Smol)
   if (allocated(cpm)) deallocate(cpm)
   if (allocated(cvm)) deallocate(cvm)
   if (allocated(hm)) deallocate(hm)

!  Reaction-index and pressure-dependent storage.
   if (allocated(itbALL)) deallocate(itbALL, ltbALL)
   if (allocated(itbFALL)) deallocate(itbFALL, ltbFALL)
   if (allocated(itbTROE)) deallocate(itbTROE, ltbTROE)
   if (allocated(iTROE4)) deallocate(iTROE4, lTROE4)
   if (allocated(itbLIND)) deallocate(itbLIND, ltbLIND)
   if (allocated(itbSIMP)) deallocate(itbSIMP, ltbSIMP)
   if (allocated(iTROEitbALL)) deallocate(iTROEitbALL)
   if (allocated(iLINDitbALL)) deallocate(iLINDitbALL)
   if (allocated(iSIMPitbALL)) deallocate(iSIMPitbALL)
   if (allocated(iFALLitbALL)) deallocate(iFALLitbALL)
   if (allocated(iTROEitbFALL)) deallocate(iTROEitbFALL)
   if (allocated(iLINDitbFALL)) deallocate(iLINDitbFALL)
   if (allocated(todotroe)) deallocate(todotroe)
   if (allocated(zeroT2)) deallocate(zeroT2)
   if (allocated(aT2)) deallocate(aT2)
   if (allocated(uT1T2)) deallocate(uT1T2)
   if (allocated(T2T2)) deallocate(T2T2)
   if (allocated(uT3T2)) deallocate(uT3T2)
   if (allocated(tab_troefactor)) deallocate(tab_troefactor)
   if (allocated(tab_log10Fcent)) deallocate(tab_log10Fcent)

   if (allocated(iTREV)) deallocate(iTREV)
   if (allocated(iXREV)) deallocate(iXREV)
   if (allocated(iEQREV)) deallocate(iEQREV)
   if (allocated(tb_beta_pack)) deallocate(tb_beta_pack)
   if (allocated(is_beta_pack)) deallocate(is_beta_pack)
   if (allocated(n_tb_beta)) deallocate(n_tb_beta)
   if (allocated(Arrhreac)) deallocate(Arrhreac)
   if (allocated(Lindreac)) deallocate(Lindreac)
   if (allocated(Troereac)) deallocate(Troereac)
   if (allocated(Revreac)) deallocate(Revreac)
   if (allocated(rate_form)) deallocate(rate_form)
   if (allocated(plog_reaction)) deallocate(plog_reaction)
   if (allocated(plog_node_ptr)) deallocate(plog_node_ptr)
   if (allocated(plog_logP)) deallocate(plog_logP)
   if (allocated(plog_term_ptr)) deallocate(plog_term_ptr)
   if (allocated(plog_A)) deallocate(plog_A)
   if (allocated(plog_b)) deallocate(plog_b)
   if (allocated(plog_EoverR)) deallocate(plog_EoverR)
   if (allocated(store_uKc)) deallocate(store_uKc)

!  Kinetics maps, tabulations, and caches.
   if (allocated(ir)) deallocate(ir)
   if (allocated(jr)) deallocate(jr)
   if (allocated(ijr)) deallocate(ijr)
   if (allocated(ip)) deallocate(ip)
   if (allocated(jp)) deallocate(jp)
   if (allocated(ijp)) deallocate(ijp)
   if (allocated(i1r)) deallocate(i1r)
   if (allocated(i1p)) deallocate(i1p)
   if (allocated(i2r)) deallocate(i2r)
   if (allocated(i2p)) deallocate(i2p)
   if (allocated(i2D1r)) deallocate(i2D1r)
   if (allocated(i2D2r)) deallocate(i2D2r)
   if (allocated(i2D1p)) deallocate(i2D1p)
   if (allocated(i2D2p)) deallocate(i2D2p)
   if (allocated(unitr)) deallocate(unitr)
   if (allocated(unitp)) deallocate(unitp)
   if (allocated(indice_r)) deallocate(indice_r)
   if (allocated(indice_p)) deallocate(indice_p)
   if (allocated(v_stoich_r)) deallocate(v_stoich_r)
   if (allocated(v_stoich_p)) deallocate(v_stoich_p)
   if (allocated(iv_stoich_r)) deallocate(iv_stoich_r)
   if (allocated(iv_stoich_p)) deallocate(iv_stoich_p)
   if (allocated(v_stoich_r1)) deallocate(v_stoich_r1)
   if (allocated(v_stoich_r2)) deallocate(v_stoich_r2)
   if (allocated(v_stoich_ro)) deallocate(v_stoich_ro)
   if (allocated(v_stoich_p1)) deallocate(v_stoich_p1)
   if (allocated(v_stoich_p2)) deallocate(v_stoich_p2)
   if (allocated(v_stoich_po)) deallocate(v_stoich_po)
   if (allocated(tab_k0)) deallocate(tab_k0)
   if (allocated(tab_kinf)) deallocate(tab_kinf)
   if (allocated(tab_dkinfdt)) deallocate(tab_dkinfdt)
   if (allocated(tab_dk0dt)) deallocate(tab_dk0dt)
   if (allocated(tab_Xkb)) deallocate(tab_Xkb)
   if (allocated(tab_Xdkbdt)) deallocate(tab_Xdkbdt)
   if (allocated(tab_kinfT)) deallocate(tab_kinfT)
   if (allocated(iA0)) deallocate(iA0)
   if (allocated(store_dwdt)) deallocate(store_dwdt)
   if (allocated(store_kinf)) deallocate(store_kinf)
   if (allocated(store_qf)) deallocate(store_qf)
   if (allocated(store_qb)) deallocate(store_qb)
   if (allocated(store_C)) deallocate(store_C)

!  Sparse and allocatable data created separately for each OpenMP thread.
!$omp parallel
   if (allocated(ljac)) deallocate(ljac)
   if (allocated(rowjac)) deallocate(rowjac)
   if (allocated(coljac)) deallocate(coljac)
   if (allocated(indices)) deallocate(indices)
   if (allocated(i2D1)) deallocate(i2D1)
   if (allocated(i2D2)) deallocate(i2D2)
   if (allocated(save_k0)) deallocate(save_k0)
   if (allocated(save_kinf)) deallocate(save_kinf)
   if (allocated(atol)) deallocate(atol)
   if (allocated(rwork)) deallocate(rwork)
   if (allocated(iwork)) deallocate(iwork)
   if (allocated(iper1)) deallocate(iper1)
   if (allocated(iiper1)) deallocate(iiper1)
   if (allocated(iper2)) deallocate(iper2)
   if (allocated(iiper2)) deallocate(iiper2)
   if (allocated(iLU1)) deallocate(iLU1)
   if (allocated(iLU2)) deallocate(iLU2)
   if (allocated(rLU1)) deallocate(rLU1)
   if (allocated(rLU2)) deallocate(rLU2)
   if (associated(vf90_opts%ATOL)) deallocate(vf90_opts%ATOL)
   if (associated(vf90_opts%RTOL)) deallocate(vf90_opts%RTOL)
   nullify(vf90_opts%ATOL, vf90_opts%RTOL)
   call sparse_reset(JAC_sparse)
   call sparse_reset(JACT_sparse)
   call sparse_reset(JACYY_sparse)
   call sparse_reset(JACYYT_sparse)
   call sparse_reset(nudiffT_molarv_sparse)
   call sparse_reset(dq_dY_sparse)
   call sparse_reset(dq_dY_T_sparse)
   sparse_jac = .false.
   njac = 0
   lsavek = .false.
   lrw = 0
   liw = 0
!$omp end parallel

!  Solver-owned sparse workspaces are shared.
   call sparse_reset(JACT_VF90)
   call sparse_reset(R5_sys1)
   call sparse_reset(R5_sys2)
   call sparse_reset(DASPK_sys)

!  Optional CHEMKIN runtime storage.
   if (allocated(intwork)) deallocate(intwork)
   if (allocated(reawork)) deallocate(reawork)
   if (allocated(chawork)) deallocate(chawork)
   if (allocated(ICK)) deallocate(ICK)
   if (allocated(RCK)) deallocate(RCK)
   if (allocated(CCK)) deallocate(CCK)
   if (allocated(CKmw)) deallocate(CKmw)
   if (allocated(ckrwork)) deallocate(ckrwork)
   if (allocated(ckiwork)) deallocate(ckiwork)

!  Reset scalar state and the interpreter-side PLOG collector.
   call plog_reset
   nel = 0
   ns = 0
   nr = 0
   neq = 0
   uneq = 0.0
   nnotrev = 0
   nTHREE = 0
   nTB = 0
   iO2 = 0
   ntbALL = 0
   ntbFALL = 0
   ntbTROE = 0
   nTROE4 = 0
   ntbLIND = 0
   ntbSIMP = 0
   nTREV = 0
   nXREV = 0
   nEQREV = 0
   n_plog_reactions = 0
   n_plog_nodes = 0
   n_plog_terms = 0
   store_uKc_T = 0.0
   store_dwdt_T = 0.0
   store_kinf_T = -1.0
   store_valid = .false.
   num_vro = 0
   num_vpo = 0
   num_vr2 = 0
   num_vp2 = 0
   iC = 0
   iO = 0
   iH = 0
   iN = 0
   iAr = 0
   ncJAC = 0
   ncCONV = 0
   nsteps = 0
   nLUdec = 0
   nNewton = 0
   ckne = 0
   ckns = 0
   cknr = 0
   ckneq = 0

end subroutine chemistry_finalize


