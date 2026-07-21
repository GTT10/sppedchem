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

c     *****************************************************************
c     **                                                             **
c     **                     SpeedCHEM FORTRAN                       **
c     **                                                             **
c     **                    Variable Allocation                      **
c     **                                                             **
c     **   Author:      (C) Federico Perini                          **
c     **   Last update: tuesday, 25/05/2010                          **
c     **                                                             **
c     *****************************************************************

      subroutine SCallocate

      use speedchem
      use SCthermodata
      use troepar
      use reacpar
      use kinetics_mod

!      use cellcluster
      implicit none  
      integer i,j

c     Variables depending on number of reactions only
      allocate(A0(nr),AREV(nr),AINF(nr),b0(nr),bREV(nr),bINF(nr),E0(nr),
     &         EREV(nr),EINF(nr))      
      allocate(reversibile(nr),REV(nr))
      allocate(HIGH(nr),LOW(nr),TROE(nr),THREE(nr))
      allocate(aTROE(nr),T1TROE(nr),T2TROE(nr),T3TROE(nr))
  

c     Variables depending on number of species only
      allocate(SCMW(ns),uMW(ns))  
      allocate(specie(ns))

c     Variables depending on number of elements only
      allocate(elementi(nel), elMW(nel), uelMW(nel))

c     Variables depending on number of reactions and species
c      allocate(!nudiff(nr,ns),
c     &         third_body(nr,ns),!reagenti(nr,ns),prodotti(nr,ns),
c     &         third_body_beta(ns,nr))!,
c     &                            third_body_betaT(nr,ns))
!      allocate(nudiffT(ns,nr))
c      allocate(nudiffmask(nr,ns))

c     Thermodynamic database
      allocate(tsw(ns),aL(ns),bL(ns),cL(ns),dL(ns),eL(ns),fL(ns),gL(ns))
      allocate(aH(ns),bH(ns),cH(ns),dH(ns),eH(ns),fH(ns),gH(ns))

c     Allocate arrays for temperature-dependent data in case
c     the "SAVE" attribute has been chosen for the Jacobian evaluation


c      if (analytical_jac .and. save_thermal_parameters) then

c         !$OMP PARALLEL
c         !$OMP CRITICAL

c        Reaction rate constants
c         allocate(save_k0(nr), save_kinf(nr))


c         save_k0   = 0.d0
c         save_kinf = 0.d0

c         !$OMP END CRITICAL
c         !$OMP END PARALLEL

c      endif



      end subroutine SCallocate



c     *****************************************************************
c     **                                                             **
c     **                     SpeedCHEM FORTRAN                       **
c     **                                                             **
c     **              Deallocate unuseful dense data                 **
c     **                                                             **
c     **   Author:      (C) Federico Perini                          **
c     **   Last update: tuesday, 15/01/2012                          **
c     **   Usage, copying, distribution of the proprietary routines  **
c     **   or any modified versions of them are not allowed without  **
c     **   approval of the author.                                   **
c     **                                                             **
c     *****************************************************************

      subroutine SCdeallocate

      use speedchem, only : third_body, third_body_beta

      implicit none


c     This statement has to follow the first call to the chemistry
c     jacobian routine, to allow the sparse jacobian structure be
c     evaluated and stored in sparse format
c      deallocate(nudiff)

c      deallocate(third_body, third_body_beta)

      end subroutine SCdeallocate


  
