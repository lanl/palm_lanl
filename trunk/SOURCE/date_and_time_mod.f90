!> @file date_and_time_mod.f90
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
! $Id: date_and_time_mod.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2701 2017-12-15 15:40:50Z suehring
! Changes from last commit documented
!
! 2698 2017-12-14 18:46:24Z suehring
! Bugfix in definition of d_seconds_year.
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! 
! 2544 2017-10-13 18:09:32Z maronga
! Initial revision
! 
! 

!
! Description:
! ------------
!> This routine calculates all needed information on date and time used by
!> other modules
!------------------------------------------------------------------------------!
 MODULE date_and_time_mod
 
    USE control_parameters,                                                    &
        ONLY: time_since_reference_point 
 
    USE kinds

    IMPLICIT NONE

    PRIVATE
    
    PUBLIC   calc_date_and_time, d_hours_day, d_seconds_hour, d_seconds_year,  &
             day_of_year, day_of_year_init, time_utc, time_utc_init


    INTEGER(iwp) ::  day_of_year              !< day of the year (DOY)
    INTEGER(iwp) ::  day_of_year_init = 172   !< DOY at model start (default: 21 June)

    REAL(wp) ::  time_utc                     !< current model time in UTC
    REAL(wp) ::  time_utc_init = 43200.0_wp   !< UTC time at model start

    REAL(wp), PARAMETER ::  d_hours_day    = 1.0_wp / 24.0_wp       !< inverse of hours per day (1/24)
    REAL(wp), PARAMETER ::  d_seconds_hour = 1.0_wp / 3600.0_wp     !< inverse of seconds per hour (1/3600)
    REAL(wp), PARAMETER ::  d_seconds_year = 1.0_wp / 31536000.0_wp !< inverse of the seconds per year (1/(365*86400))
    
    SAVE

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate current day of the year and time in UTC
!------------------------------------------------------------------------------!
 
    SUBROUTINE calc_date_and_time

       IMPLICIT NONE

!
!--    Calculate current day of the year 
       day_of_year = day_of_year_init + INT(FLOOR( (time_utc_init + time_since_reference_point)&
                               / 86400.0_wp ), KIND=iwp)
                        
       time_utc = MOD((time_utc_init + time_since_reference_point), 86400.0_wp)
       
       

    END SUBROUTINE calc_date_and_time



 END MODULE date_and_time_mod
