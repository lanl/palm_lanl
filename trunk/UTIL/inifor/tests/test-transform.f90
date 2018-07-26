!> @file tests/test-transform.f90
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
! $Id: test-transform.f90 2718 2018-01-02 08:49:38Z maronga $
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
!> This program tests INIFOR's rotated-pole coordinate transforms.
!------------------------------------------------------------------------------!
 PROGRAM test_transform

    USE defs, ONLY :  TO_RADIANS, TO_DEGREES
    USE grid, ONLY :  grid_definition, init_grid_definition
    USE transform, ONLY :  phi2phirot, rla2rlarot, phirot2phi, rlarot2rla
    USE test_utils
    
    IMPLICIT NONE

    CHARACTER(LEN=30) ::  title
    LOGICAL           ::  res

    REAL, PARAMETER ::  lx = 100., ly = 200., lz = 300.

    ! Angels in degrees
    REAL ::  phi, lambda, phi_c, lambda_c, phi2, lambda2, phi_n, lambda_n

    title = "rotation north-east"
    CALL begin_test(title, res)
    ! Arange
    phi    = 52.5166670000000000
    lambda = 13.3833330000000000

    phi_n    =   40.
    lambda_n = -170.

    ! Act
    phi_c     = phi2phirot(phi,   lambda, phi_n, lambda_n)
    lambda_c  = rla2rlarot(phi,   lambda, phi_n, lambda_n, 0.)

    phi2      = phirot2phi(phi_c, lambda_c, phi_n, lambda_n, 0.)
    lambda2   = rlarot2rla(phi_c, lambda_c, phi_n, lambda_n, 0.)

    ! Assert
    res = assert_equal( (/phi, lambda/), (/phi2, lambda2/),  "rotated grid transformations" )
    PRINT *, " Angles before transformation:  ", phi, lambda
    PRINT *, " and after back transformation: ", phi2, lambda2

    CALL end_test(title, res)

    title = "rotation south-west"
    CALL begin_test(title, res)
    ! Arange
    phi    = 49.
    lambda =  9.

    phi_n    =   40.
    lambda_n = -170.

    ! Act
    phi_c     = phi2phirot(phi,   lambda, phi_n, lambda_n)
    lambda_c  = rla2rlarot(phi,   lambda, phi_n, lambda_n, 0.)

    phi2      = phirot2phi(phi_c, lambda_c, phi_n, lambda_n, 0.)
    lambda2   = rlarot2rla(phi_c, lambda_c, phi_n, lambda_n, 0.)

    ! Assert
    res = assert_equal( (/phi, lambda/), (/phi2, lambda2/),  "rotated grid transformations" )
    PRINT *, " Angles before transformation:  ", phi, lambda
    PRINT *, " and after back transformation: ", phi2, lambda2

    CALL end_test(title, res)
 END PROGRAM test_transform
