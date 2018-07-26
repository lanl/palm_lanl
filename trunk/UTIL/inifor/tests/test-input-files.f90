!> @file tests/test-input-files.f90
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
! $Id: test-input-files.f90 2718 2018-01-02 08:49:38Z maronga $
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
!> This program tests INIFOR's timestamping used for generating input file
!> names.
!------------------------------------------------------------------------------!
 PROGRAM test_input_files

    USE defs,                                                                  &
        ONLY :  PATH
    USE grid,                                                                  & 
        ONLY :  input_file_list
    USE test_utils
    
    IMPLICIT NONE

    CHARACTER(LEN=50)                              ::  title
    CHARACTER(LEN=PATH), ALLOCATABLE, DIMENSION(:) ::  file_list, ref_list
    LOGICAL                                        ::  res
    INTEGER                                        ::  i     

    title = "input files - daylight saving to standard time"
    CALL begin_test(title, res)

    ! Arange
    ! ...a date range that inlcudes a shift from daylight saving time to
    ! standard time (29.10.2017). Since all time stamps in COSMO-DE input files
    ! are in UTC, this should not the naming cadence.
    ALLOCATE( ref_list(6) )
    ref_list(1)  = './laf2017102823-test.nc'
    ref_list(2)  = './laf2017102900-test.nc'
    ref_list(3)  = './laf2017102901-test.nc'
    ref_list(4)  = './laf2017102902-test.nc'
    ref_list(5)  = './laf2017102903-test.nc'
    ref_list(6)  = './laf2017102904-test.nc'

    ! Act
    CALL input_file_list(start_date_string='2017102823',                       &
                         start_hour=0, end_hour=5, step_hour=1,                &
                         path='./', prefix="laf", suffix='-test',              &
                         file_list=file_list)

    ! Assert
    DO i = 1, 6
       res = res .AND. (TRIM(ref_list(i)) .EQ. TRIM(file_list(i)))
    END DO

    DEALLOCATE( ref_list, file_list )
    CALL end_test(title, res)


    title = "input files - leap day"
    CALL begin_test(title, res)

    ! Arange
    ! ...a date range that inlcudes a leap day (29. Feb. 2016) which should be
    ! inlcuded in UTC time stamps.
    ALLOCATE( ref_list(2) )
    ref_list(1)  = './laf2016022823-test.nc'
    ref_list(2)  = './laf2016022900-test.nc'

    ! Act
    CALL input_file_list(start_date_string='2016022823',                       &
                         start_hour=0, end_hour=1, step_hour=1,                &
                         path='./', prefix="laf", suffix='-test',              &
                         file_list=file_list)

    ! Assert
    DO i = 1, 2
       res = res .AND. (TRIM(ref_list(i)) .EQ. TRIM(file_list(i)))
    END DO

    DEALLOCATE( ref_list, file_list )
    CALL end_test(title, res)

 END PROGRAM test_input_files
