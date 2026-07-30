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
!     **   SpeedCHEM - MPI parallel detailed chemistry link          **
!     **                                                             **
!     **   This subroutine broadcasts the reaction mechanism data    **
!     **   of chemistry modules from node 0 to all the other nodes.  **
!     **   The routine contains MPI_BCAST commands and should be     **
!     **   executed by all the instances of the program              **
!     **                                                             **
!     **                                                             **
!     **   Author:      Federico Perini                              **
!     **   Last update: friday, 01/06/2012                           **
!     **                                                             **
!     *****************************************************************

subroutine SCbroadcast

   use sparse_mpi
   use chemistry_setup
   use speedchem
   use sparse_chemistry
   use SCthermodata
   use troepar
   use reacpar
   use kinetics_mod
   use ode_solver

   implicit none

   integer :: j

!     ** Broadcast mechanism setup data ** (from chemistry_setup) *****

   call mpi_broadcast(80,mechanism)
   call mpi_broadcast(15,solver   )

   call mpi_broadcast(use_speedchem)
   call mpi_broadcast(accurate_scthermo)
   call mpi_broadcast(analytical_jac)
   call mpi_broadcast(simplified_for_sparsity)
   call mpi_broadcast(save_thermal_parameters)
   call mpi_broadcast(check_reaction_mechanism)
   call mpi_broadcast(print_out_jacobian)
   call mpi_broadcast(separate_tols)
   call mpi_broadcast(permutate_species)

   call mpi_broadcast(TOLR)
   call mpi_broadcast(YTOLA)
   call mpi_broadcast(TTOLA)
   call mpi_broadcast(nwatch)

   call mpi_broadcast(Temp_table_accuracy)
   call mpi_broadcast(rec_Ttable_accuracy)
   call mpi_broadcast(Temp_HIlim)
   call mpi_broadcast(Temp_LOlim)
   call mpi_broadcast(tab_nsteps)


!     ** Broadcast mechanism information ** (from speedchem) **********
   call mpi_broadcast(nel)
   call mpi_broadcast(ns)
   call mpi_broadcast(nr)

   call mpi_broadcast(species)
   call mpi_broadcast(reactions)
   call mpi_broadcast(species_permutations)
   call mpi_broadcast(species_inverse_permutations)

   call mpi_broadcast(neq)
   call mpi_broadcast(uneq)

   call mpi_broadcast(A0)
   call mpi_broadcast(AREV)
   call mpi_broadcast(AINF)
   call mpi_broadcast(ARi)
   call mpi_broadcast(b0)
   call mpi_broadcast(bREV)
   call mpi_broadcast(bINF)
   call mpi_broadcast(bRi)
   call mpi_broadcast(E0)
   call mpi_broadcast(EREV)
   call mpi_broadcast(EINF)
   call mpi_broadcast(ERi)
   call mpi_broadcast(SCMW)
   call mpi_broadcast(uMW)
   call mpi_broadcast(elMW)
   call mpi_broadcast(uelMW)

   call mpi_broadcast(reversibile)
   call mpi_broadcast(REV)
   call mpi_broadcast(duplicate)
   call mpi_broadcast(inotrev)
   call mpi_broadcast(nnotrev)

   call mpi_broadcast(stoich_r_pack)
   call mpi_broadcast(stoich_p_pack)
   call mpi_broadcast(nudiff_pack)

   call mpi_broadcast(i_stoich_r)
   call mpi_broadcast(i_stoich_p)
   call mpi_broadcast(i_nudiff)

   call mpi_broadcast(n_stoich_r)
   call mpi_broadcast(n_stoich_p)
   call mpi_broadcast(n_nudiff)

   call mpi_broadcast(isumnudiff)
   call mpi_broadcast(HIGH)
   call mpi_broadcast(LOW)
   call mpi_broadcast(TROE)
   call mpi_broadcast(THREE)

   call mpi_broadcast(aTROE)
   call mpi_broadcast(T1TROE)
   call mpi_broadcast(T2TROE)
   call mpi_broadcast(T3TROE)

   call mpi_broadcast(iTHREE)
   call mpi_broadcast(iTB)

   call mpi_broadcast(nTHREE)
   call mpi_broadcast(nTB)

   call mpi_broadcast(18,specie)
   call mpi_broadcast(2 ,elementi)

   call mpi_broadcast(sparse_jac)
   call mpi_broadcast(njac)
   call mpi_broadcast(iO2)

!     ** Broadcast sparse chemistry information (sparse_chemistry) ****
   call mpi_broadcast(nudiff_sparse)
   call mpi_broadcast(nudiffT_sparse)
   call mpi_broadcast(nudiff_EQREV_sp)
   call mpi_broadcast(stoich_r_sp)
   call mpi_broadcast(stoich_p_sp)
   call mpi_broadcast(stoich_r_eff_sp)
   call mpi_broadcast(stoich_p_eff_sp)
   call mpi_broadcast(EM_sp)

   call mpi_broadcast(JAC_sparse)
   call mpi_broadcast(JACT_sparse)
   call mpi_broadcast(JACYY_sparse)
   call mpi_broadcast(JACYYT_sparse)
   call mpi_broadcast(dq_dY_sparse)
   call mpi_broadcast(dq_dY_T_sparse)

   call mpi_broadcast(third_body_sp)
   call mpi_broadcast(tb_beta_sp)
   call mpi_broadcast(tbb_uMW_sp)

   call mpi_broadcast(inudiff_sparse)
   call mpi_broadcast(inudiffT_sparse)
   call mpi_broadcast(istoich_r_sp)
   call mpi_broadcast(istoich_p_sp)
   call mpi_broadcast(istoich_r_eff_sp)
   call mpi_broadcast(istoich_p_eff_sp)

!     ** Broadcast thermodynamic data information (SCthermodata) ******
   call mpi_broadcast(tsw)
   call mpi_broadcast(aL)
   call mpi_broadcast(bL)
   call mpi_broadcast(cL)
   call mpi_broadcast(dL)
   call mpi_broadcast(eL)
   call mpi_broadcast(fL)
   call mpi_broadcast(gL)
   call mpi_broadcast(aH)
   call mpi_broadcast(bH)
   call mpi_broadcast(cH)
   call mpi_broadcast(dH)
   call mpi_broadcast(eH)
   call mpi_broadcast(fH)
   call mpi_broadcast(gH)

   call mpi_broadcast(tab_CpsuR)
   call mpi_broadcast(tab_HsuRT)
   call mpi_broadcast(tab_SsuR)
   call mpi_broadcast(tab_dGdt)
   call mpi_broadcast(tab_dDG0dt)
   call mpi_broadcast(tab_dCvdT)
   call mpi_broadcast(tab_uKp)
   call mpi_broadcast(tab_uKc)
   call mpi_broadcast(tab_dKcdT)
   call mpi_broadcast(EM)
   call mpi_broadcast(Emfr)

!     ** Broadcast TROE reactions parameters (module troepar) *********
   call mpi_broadcast(ntbALL)
   call mpi_broadcast(itbALL)
   call mpi_broadcast(ltbALL)

   call mpi_broadcast(ntbFALL)
   call mpi_broadcast(itbFALL)
   call mpi_broadcast(ltbFALL)

   call mpi_broadcast(ntbTROE)
   call mpi_broadcast(itbTROE)
   call mpi_broadcast(ltbTROE)

   call mpi_broadcast(nTROE4)
   call mpi_broadcast(iTROE4)
   call mpi_broadcast(lTROE4)

   call mpi_broadcast(ntbLIND)
   call mpi_broadcast(itbLIND)
   call mpi_broadcast(ltbLIND)

   call mpi_broadcast(ntbSIMP)
   call mpi_broadcast(itbSIMP)
   call mpi_broadcast(ltbSIMP)

   call mpi_broadcast(iTROEitbALL)
   call mpi_broadcast(iLINDitbALL)
   call mpi_broadcast(iSIMPitbALL)
   call mpi_broadcast(iFALLitbALL)
   call mpi_broadcast(iTROEitbFALL)
   call mpi_broadcast(iLINDitbFALL)

   call mpi_broadcast(todotroe)
   call mpi_broadcast(zeroT2)

   call mpi_broadcast(aT2)
   call mpi_broadcast(uT1T2)
   call mpi_broadcast(T2T2)
   call mpi_broadcast(uT3T2)

   call mpi_broadcast(tab_log10Fcent)

!     ** Broadcast reaction rate parameters ***************************
   call mpi_broadcast(nTREV)
   call mpi_broadcast(iTREV)
   call mpi_broadcast(nXREV)
   call mpi_broadcast(iXREV)
   call mpi_broadcast(nEQREV)
   call mpi_broadcast(iEQREV)

   call mpi_broadcast(tb_beta_pack)
   call mpi_broadcast(is_beta_pack)
   call mpi_broadcast(n_tb_beta)

   call mpi_broadcast(Arrhreac)
   call mpi_broadcast(Lindreac)
   call mpi_broadcast(Troereac)
   call mpi_broadcast(Revreac )

!  PLOG (cklink v2): per-reaction rate-form tag plus
!  the packed pressure-dependent-Arrhenius arrays. Scalars first so
!  workers know the sizes; the array broadcasts preserve lbounds, so the
!  0-based *_ptr arrays keep their CSR indexing. A no-PLOG mechanism
!  broadcasts empty/zero data and never uses it.
   call mpi_broadcast(rate_form)
   call mpi_broadcast(n_plog_reactions)
   call mpi_broadcast(n_plog_nodes)
   call mpi_broadcast(n_plog_terms)
   call mpi_broadcast(plog_reaction)
   call mpi_broadcast(plog_node_ptr)
   call mpi_broadcast(plog_term_ptr)
   call mpi_broadcast(plog_logP)
   call mpi_broadcast(plog_A)
   call mpi_broadcast(plog_b)
   call mpi_broadcast(plog_EoverR)


   call mpi_broadcast(ip)
   call mpi_broadcast(jp)
   call mpi_broadcast(ir)
   call mpi_broadcast(jr)
   call mpi_broadcast(v_stoich_p)
   call mpi_broadcast(v_stoich_r)
   call mpi_broadcast(indice_r)
   call mpi_broadcast(indice_p)
   call mpi_broadcast(ijp)
   call mpi_broadcast(ijr)
   call mpi_broadcast(iv_stoich_p)
   call mpi_broadcast(iv_stoich_r)
   call mpi_broadcast(i2D1p)
   call mpi_broadcast(i2D1r)
   call mpi_broadcast(i2D2p)
   call mpi_broadcast(i2D2r)
   call mpi_broadcast(i1p)
   call mpi_broadcast(i1r)
   call mpi_broadcast(i2p)
   call mpi_broadcast(i2r)

   call mpi_broadcast(tab_k0)
   call mpi_broadcast(tab_kinf)
   call mpi_broadcast(tab_dkinfdt)
   call mpi_broadcast(tab_dk0dt)
   call mpi_broadcast(tab_Xkb)
   call mpi_broadcast(tab_Xdkbdt)

   call mpi_broadcast(iA0)


!     ** Broadcast ODE solver variables *******************************

   call mpi_broadcast(lrw)
   call mpi_broadcast(liw)
   call mpi_broadcast(method)
   call mpi_broadcast(itol)
   call mpi_broadcast(iopt)
   call mpi_broadcast(itask)
   call mpi_broadcast(istate)
   call mpi_broadcast(ijac)
   call mpi_broadcast(imas)
   call mpi_broadcast(iout)
   call mpi_broadcast(ifcn)
   call mpi_broadcast(idfx)
   call mpi_broadcast(mumas)
   call mpi_broadcast(mlmas)
   call mpi_broadcast(mujac)
   call mpi_broadcast(mljac)
   call mpi_broadcast(maxk)
   call mpi_broadcast(maxnsteps)

   call mpi_broadcast(hs)
   call mpi_broadcast(rtol)
   call mpi_broadcast(atol)

   call mpi_broadcast(rwork)
   call mpi_broadcast(iwork)

   call mpi_broadcast(iper1)
   call mpi_broadcast(iper2)
   call mpi_broadcast(iiper1)
   call mpi_broadcast(iiper2)
   call mpi_broadcast(iLU1)
   call mpi_broadcast(iLU2)

   call mpi_broadcast(rLU1)
   call mpi_broadcast(rLU2)

   call mpi_broadcast(liLU1)
   call mpi_broadcast(liLU2)
   call mpi_broadcast(lrLU1)
   call mpi_broadcast(lrLU2)
   call mpi_broadcast(strg1)
   call mpi_broadcast(strg2)
   call mpi_broadcast(maxnthreads)
   call mpi_broadcast(actnthreads)


!     ** Broadcast VODE-F90 related data *******************************
   call mpi_broadcast(nIA)
   call mpi_broadcast(nJA)
   call mpi_broadcast(nPD)
   call mpi_broadcast(nVF90JAC)
   call mpi_broadcast(vf90_opts%MF)
   call mpi_broadcast(vf90_opts%METH)
   call mpi_broadcast(vf90_opts%MITER)
   call mpi_broadcast(vf90_opts%MOSS)
   call mpi_broadcast(vf90_opts%ITOL)
   call mpi_broadcast(vf90_opts%IOPT)
   call mpi_broadcast(vf90_opts%NG)
   call mpi_broadcast(vf90_opts%DENSE)
   call mpi_broadcast(vf90_opts%BANDED)
   call mpi_broadcast(vf90_opts%SPARSE)
   call mpi_pointer_broadcast(vf90_opts%RTOL)
   call mpi_pointer_broadcast(vf90_opts%ATOL)

   if (current_MPI_CPUID() > 0) call ODE_solver_speedchem_init


!     ** Broadcast DASPK related data *********************************
   DASPK: do j = 1, 20
      call mpi_broadcast(daspkinfo(j))
   end do DASPK

!     *****************************************************************
end subroutine SCbroadcast

