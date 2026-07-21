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

!   ********************************************************************
!   **                                                                **
!   **                Sparse ODE - Working precision                  **
!   **   This module contains the constants that define the desired   **
!   **   working precision, and should be linked to from every other  **
!   **   module.                                                      **
!   **                                                                **
!   **                                                                **
!   **   Author:      (C) Federico Perini                             **
!   **   Last update: thursday, 16/08/2012                            **
!   **                                                                **
!   **                                                                **
!   ********************************************************************

    module working_precision
        implicit none
        private

        ! Parameters setting the integer and real precision.
        ! Now working in quadruple precision
        integer, parameter, public :: int_p = kind(12345)
        integer, parameter, public :: dp    = 8

        ! Integer parameter representing the ratio in length
        ! between an real and an integer number storages
        integer, parameter, public :: real_to_int_length = dp/int_p

        ! ** CURRENT MACHINE NUMBERS CONSTANTS *************************
        real (dp), parameter, public ::                                &
           small  = tiny     (0.0_dp),                                 &
           big    = huge     (0.0_dp),                                 &
           uround = epsilon  (0.0_dp)

        integer (int_p), parameter, public ::                          &
           rdigit = digits   (0.0_dp),                                 &
           rprec  = precision(0.0_dp)

        ! ** PUBLIC PROCEDURES *****************************************
        public :: print_working_precision_properties

        ! ** CONVERSION TO AND FROM DOUBLE PRECISION FORMAT ************
!        interface assignment(=)
!!            module procedure double_to_arbitrary
!            module procedure arbitrary_to_double
!        end interface

        contains


        ! **************************************************************
        ! ** Conversion of arbitrary precision number to double       **
        ! ** and viceversa                                            **
        ! **************************************************************

!        elemental subroutine double_to_arbitrary(arbitrary,doble)
!        implicit none
!
!        double precision, intent(in)  :: doble
!        real (kind = dp), intent(out) :: arbitrary
!
!        arbitrary = real(doble, kind = dp)
!
!        end subroutine double_to_arbitrary


!        elemental subroutine arbitrary_to_double(doble,arbitrary)
!        implicit none

!        real (kind = dp), intent(in)   :: arbitrary
!        double precision, intent(out)  :: doble

!        doble = dble(arbitrary)

!        end subroutine arbitrary_to_double


        ! **************************************************************
        ! ** Print to screen the current working precision            **
        ! ** details                                                  **
        ! **************************************************************
        subroutine print_working_precision_properties(desired_unit)
            implicit none

            integer (int_p), intent(in), optional :: desired_unit
            integer (int_p)                       :: this_unit


            character (len=*), parameter ::                            &
               fmt_head = "(72('-'))",                                 &
               fmt_hed2 = "(19x,'Current Working Precision Details')", &
               fmt_dig  = "(' Digit storage bits        :',I15)",      &
               fmt_eps  = "(' Machine round off value   :',1pe15.5e4)",&
               fmt_tiny = "(' Smallest positive number  :',1pe15.5e4)",&
               fmt_huge = "(' Biggest positive number   :',1pe15.5e4)",&
               fmt_prec = "(' Decimal precision         :',I15)"      ,&
               fmt_kind = "(' Fortran real KIND value   :',I15)"      ,&
               fmt_lrat = "(' Real to integer word ratio:',I15)"

            ! Choose unit for output (default = screen, unit=6)
            this_unit = 6
            if (present(desired_unit)) this_unit = desired_unit

            write(this_unit, fmt_head)
            write(this_unit, fmt_hed2)
            write(this_unit, fmt_head)
            write(this_unit,        *)

            write(this_unit, fmt_kind)dp
            write(this_unit, fmt_dig )rdigit
            write(this_unit, fmt_eps )uround
            write(this_unit, fmt_tiny)small
            write(this_unit, fmt_huge)big
            write(this_unit, fmt_prec)rprec
            write(this_unit, fmt_lrat)real_to_int_length

            write(this_unit,        *)
            write(this_unit, fmt_head)
            write(this_unit,        *)

        end subroutine print_working_precision_properties



    end module working_precision
