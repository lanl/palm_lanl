!> @file user_init_flight.f90
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
! $Id: user_init_flight.f90 2718 2018-01-02 08:49:38Z maronga $
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
! 1957 2016-07-07 10:43:48Z suehring
! Initial revision
!
! Description:
! ------------
!> Execution of user-defined initialization for flight measurements.
!------------------------------------------------------------------------------!
 SUBROUTINE user_init_flight( init, k, id, label_leg )
 

    USE control_parameters
    
    USE indices
    
    USE kinds
    
    USE netcdf_interface,                                                      &
        ONLY: dofl_label, dofl_unit
    
    USE user

    IMPLICIT NONE
    
    CHARACTER(LEN=10), OPTIONAL ::  label_leg     !< label of the respective leg
   
    INTEGER(iwp), OPTIONAL                ::  id  !< variable index
    INTEGER(iwp), OPTIONAL, INTENT(INOUT) ::  k   !< index for respective variable and leg
    
    LOGICAL ::  init  !< variable to recognize initial call 
!
!-- Sample for user-defined flight-time series.
!-- For each quantity you have to give a label and a unit, which will be used
!-- for the output into NetCDF file. They must not contain more than
!-- twenty characters. 


    IF ( init )  THEN
!
!--    The number of user-defined quantity has to be increased appropriately.
!--    In the following example, 2 user-defined quantities are added.
!        num_var_fl_user = num_var_fl_user + 2 

       init = .FALSE.
    
    ELSE
    
!
!--    Please add the respective number of new variables as following:
      
!        SELECT CASE ( id )
!        
!           CASE ( 1 )
!              dofl_label(k)   = TRIM(label_leg) // '_' // 'abs_u'
!              dofl_unit(k)    = 'm/s'
!              k               = k + 1
!              
!           CASE ( 2 )
!     
!              dofl_label(k)   = TRIM(label_leg) // '_' // 'abs_v'
!              dofl_unit(k)    = 'm/s'
!              k               = k + 1
!              
!        END SELECT

    ENDIF
       
 END SUBROUTINE user_init_flight

