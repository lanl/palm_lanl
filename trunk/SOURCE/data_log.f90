!> @file data_log.f90
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
! Copyright 1997-2018 Leibniz Universitaet Hannover
!------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: data_log.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2101 2017-01-05 16:42:31Z suehring
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! RCS Log replace by Id keyword, revision history cleaned up
!
! Revision 1.1  2006/02/23 10:09:29  raasch
! Initial revision
!
!
! Description:
! ------------
!> Complete logging of data
!------------------------------------------------------------------------------!
 SUBROUTINE data_log( array, i1, i2, j1, j2, k1, k2 )
 
#if defined( __logging )

    USE control_parameters,                                                    &
        ONLY:  log_message, simulated_time
        
    USE kinds
        
    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i1  !< 
    INTEGER(iwp) ::  i2  !< 
    INTEGER(iwp) ::  j1  !< 
    INTEGER(iwp) ::  j2  !< 
    INTEGER(iwp) ::  k1  !< 
    INTEGER(iwp) ::  k2  !< 

    REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2) ::  array  !< 


!
!-- Open the file for data logging
    CALL check_open( 20 )

!
!-- Write the message string
    WRITE ( 20 )  log_message

!
!-- Write the simulated time and the array indices
    WRITE ( 20 )  simulated_time, i1, i2, j1, j2, k1, k2

!
!-- Write the array
    WRITE ( 20 )  array

#endif
 END SUBROUTINE data_log



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Complete logging of data for 2d arrays
!------------------------------------------------------------------------------!
 
 SUBROUTINE data_log_2d( array, i1, i2, j1, j2)

#if defined( __logging )

    USE control_parameters,                                                    &
        ONLY:  log_message, simulated_time

    USE kinds
            
    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i1  !< 
    INTEGER(iwp) ::  i2  !< 
    INTEGER(iwp) ::  j1  !< 
    INTEGER(iwp) ::  j2  !< 

    REAL(wp), DIMENSION(i1:i2,j1:j2) ::  array  !< 


!
!-- Open the file for data logging
    CALL check_open( 20 )

!
!-- Write the message string
    WRITE ( 20 )  log_message

!
!-- Write the simulated time and the array indices
    WRITE ( 20 )  simulated_time, i1, i2, j1, j2

!
!-- Write the array
    WRITE ( 20 )  array

#endif
 END SUBROUTINE data_log_2d



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Complete logging of data for 2d integer arrays
!------------------------------------------------------------------------------!
 
 SUBROUTINE data_log_2d_int( array, i1, i2, j1, j2)

#if defined( __logging )

    USE control_parameters,                                                    &
        ONLY:  log_message, simulated_time

    USE kinds
            
    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i1  !< 
    INTEGER(iwp) ::  i2  !< 
    INTEGER(iwp) ::  j1  !< 
    INTEGER(iwp) ::  j2  !< 

    INTEGER(iwp), DIMENSION(i1:i2,j1:j2) ::  array  !< 


!
!-- Open the file for data logging
    CALL check_open( 20 )

!
!-- Write the message string
    WRITE ( 20 )  log_message

!
!-- Write the simulated time and the array indices
    WRITE ( 20 )  simulated_time, i1, i2, j1, j2

!
!-- Write the array
    WRITE ( 20 )  array

#endif
 END SUBROUTINE data_log_2d_int
