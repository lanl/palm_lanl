!> @file tests/test-centre-velocites.f90
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
! $Id: test-centre-velocities.f90 2718 2018-01-02 08:49:38Z maronga $
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
!> This program tests INIFOR's central velocity interpolation.
!------------------------------------------------------------------------------!
 PROGRAM test_centre_velocities

    USE test_utils
    USE transform, &
        ONLY :  centre_velocities
    
    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER ::  title = "centre velocities"
    LOGICAL                     ::  res
    REAL, DIMENSION(3,3,1)      ::  u_face, u_centre, u_ref
    REAL, DIMENSION(3,3,1)      ::  v_face, v_centre, v_ref
    INTEGER                     ::  i

    CALL begin_test(title, res)

    ! Arange
    u_face = RESHAPE( (/1, 2, 3, 1, 2, 3, 1, 2, 3/), SHAPE(u_face))
    v_face = RESHAPE( (/1, 1, 1, 2, 2, 2, 3, 3, 3/), SHAPE(v_face))

    u_ref = RESHAPE( (/0.0, 1.5, 2.5, 0.0, 1.5, 2.5, 0.0, 1.5, 2.5/), SHAPE(u_ref))
    v_ref = RESHAPE( (/0.0, 0.0, 0.0, 1.5, 1.5, 1.5, 2.5, 2.5, 2.5/), SHAPE(v_ref))

    u_centre = 0.0
    v_centre = 0.0

    ! Act
    CALL centre_velocities(u_face, v_face, u_centre, v_centre)

    ! Assert that the correct central velocities u_centre and v_centre are computed.
    DO i = 1, 3
       res = res .AND. assert_equal(u_centre(:,i,1), u_ref(:,i,1), 'centering u')
       res = res .AND. assert_equal(v_centre(i,:,1), v_ref(i,:,1), 'centering v')
    END DO

    CALL end_test(title, res)

 END PROGRAM test_centre_velocities
