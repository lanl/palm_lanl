!> @file tests/test-prototype.f90
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
! $Id: test-prototype.f90 2718 2018-01-02 08:49:38Z maronga $
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
!> This prototype test does nothing, but serves as a template for programming
!> new INIFOR tests.
!------------------------------------------------------------------------------!
 PROGRAM test_prototype

    USE test_utils
    
    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER ::  title = "Prototype"
    LOGICAL                     ::  res

    CALL begin_test(title, res)

    ! Arange

    ! Act

    ! Assert

    CALL end_test(title, res)

 END PROGRAM test_prototype
