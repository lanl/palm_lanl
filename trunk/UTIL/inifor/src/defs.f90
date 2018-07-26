!> @file src/defs.f90
!------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 2017-2018 Leibniz Universitaet Hannover
! Copyright 2017-2018 Deutscher Wetterdienst Offenbach
!------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: defs.f90 2718 2018-01-02 08:49:38Z maronga $
! Initial revision
!
! 
!
! Authors:
! --------
! @author Eckhard Kadasch
!
! Description:
! ------------
!> The defs module provides global constants used in INIFOR.
!------------------------------------------------------------------------------!
 MODULE defs
 
 IMPLICIT NONE
 
 ! Parameters for type definitions
 INTEGER, PARAMETER  ::  dp    = 8   !< double precision (8 bytes = 64 bits)
 INTEGER, PARAMETER  ::  sp    = 4   !< single precision (4 bytes = 32 bits)
 INTEGER, PARAMETER  ::  hp    = 2   !< half precision (2 bytes = 16 bits)
 INTEGER, PARAMETER  ::  PATH  = 140 !< length of file path strings
 INTEGER, PARAMETER  ::  LNAME = 150 !< length of long name strings
 INTEGER, PARAMETER  ::  SNAME = 40  !< length of short name strings
 INTEGER, PARAMETER  ::  DATE  = 10  !< length of date strings

 ! Trigonomentry
 REAL(dp), PARAMETER ::  PI = 3.14159265358979323846264338_dp !< Ratio of a circle's circumference to its diamter [-]
 REAL(dp), PARAMETER ::  TO_RADIANS = PI / 180.0_dp           !< Conversion factor from degrees to radiant [-]
 REAL(dp), PARAMETER ::  TO_DEGREES = 180.0_dp / PI           !< Conversion factor from radians to degrees [-]

 ! COSMO-DE parameters
 INTEGER, PARAMETER  ::  WATER_ID = 9                         !< Integer corresponding to the water soil type in COSMO-DE [-]
 REAL(dp), PARAMETER ::  EARTH_RADIUS = 6371229.0_dp          !< Earth radius used in COSMO-DE [m]
 REAL(dp), PARAMETER ::  P_SL = 1e5_dp                        !< Reference pressure for computation of COSMO-DE's basic state pressure [Pa]
 REAL(dp), PARAMETER ::  T_SL = 288.15_dp                     !< Reference temperature for computation of COSMO-DE's basic state pressure [K]
 REAL(dp), PARAMETER ::  BETA = 42.0_dp                       !< logarithmic lapse rate, dT / d ln(p), for computation of COSMO-DE's basic state pressure [K]
 REAL(dp), PARAMETER ::  RD   = 287.05_dp                     !< specific gar constant of dry air, used in computation of COSMO-DE's basic state [J/kg/K]
 REAL(dp), PARAMETER ::  G    = 9.80665_dp                    !< acceleration of Earth's gravity, used in computation of COSMO-DE's basic state [m/s/s]
 REAL(dp), PARAMETER ::  RHO_L = 1e3_dp                       !< density of liquid water, used to convert W_SO from [kg/m^2] to [m^3/m^3], in [kg/m^3]

 ! PALM-4U parameters
 REAL(dp), PARAMETER ::  P_REF   = 1e5_dp                     !< Reference pressure for potential temperature [Pa]
 REAL(dp), PARAMETER ::  RD_PALM = 287.0_dp                   !< specific gas constant of dry air, used in computation of PALM-4U's potential temperature [J/kg/K]
 REAL(dp), PARAMETER ::  CP_PALM = 1005.0_dp                  !< heat capacity of dry air at constant pressure, used in computation of PALM-4U's potential temperature [J/kg/K]

 ! INIFOR parameters
 INTEGER, PARAMETER          ::  FILL_ITERATIONS = 10         !< Number of iterations for extrapolating soil data into COSMO-DE water cells [-]
 REAL(dp), PARAMETER         ::  FORCING_FREQ = 3600.0_dp     !< Reference pressure for potential temperature [Pa]
 CHARACTER(LEN=*), PARAMETER ::  VERSION = '1.1.4'            !< path to script for generating input file names

 END MODULE defs
