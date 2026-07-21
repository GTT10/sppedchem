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



