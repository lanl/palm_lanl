!> @file timestep_scheme_steering.f90
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
! ------------------
! 
! 
! Former revisions:
! -----------------
! $Id: timestep_scheme_steering.f90 2718 2018-01-02 08:49:38Z maronga $
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
! 1322 2014-03-20 16:38:49Z raasch
! REAL constants defined as wp-kind
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1001 2012-09-13 14:08:46Z raasch
! all actions concerning leapfrog scheme removed
!
! Revision 1.1  2004/01/28 15:34:47  raasch
! Initial revision
!
!
! Description:
! ------------
!> Depending on the timestep scheme set the steering factors for the prognostic
!> equations.
!------------------------------------------------------------------------------!
 SUBROUTINE timestep_scheme_steering
 

    USE control_parameters,                                                    &
        ONLY:  intermediate_timestep_count, timestep_scheme, tsc

    USE kinds

    IMPLICIT NONE


    IF ( timestep_scheme(1:5) == 'runge' )  THEN
!
!--    Runge-Kutta schemes (here the factors depend on the respective
!--    intermediate step)
       IF ( timestep_scheme == 'runge-kutta-2' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             tsc(1:5) = (/ 1.0_wp, 1.0_wp,  0.0_wp, 0.0_wp, 0.0_wp /)
          ELSE
             tsc(1:5) = (/ 1.0_wp, 0.5_wp, -0.5_wp, 0.0_wp, 1.0_wp /)
          ENDIF
       ELSE
          IF ( intermediate_timestep_count == 1 )  THEN
             tsc(1:5) = (/ 1.0_wp,  1.0_wp /  3.0_wp,           0.0_wp, 0.0_wp, 0.0_wp /)
          ELSEIF ( intermediate_timestep_count == 2 )  THEN
             tsc(1:5) = (/ 1.0_wp, 15.0_wp / 16.0_wp, -25.0_wp/48.0_wp, 0.0_wp, 0.0_wp /)
          ELSE
             tsc(1:5) = (/ 1.0_wp,  8.0_wp / 15.0_wp,   1.0_wp/15.0_wp, 0.0_wp, 1.0_wp /)
          ENDIF          
       ENDIF

    ELSEIF ( timestep_scheme == 'euler' )  THEN
!
!--    Euler scheme
       tsc(1:5) = (/ 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp /)

    ENDIF


 END SUBROUTINE timestep_scheme_steering
