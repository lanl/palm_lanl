!> @file user_write_restart_data_mod.f90
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
! $Id: user_write_restart_data_mod.f90 2894 2018-03-15 09:17:58Z Giersch $
! Initial revision
! 
!
! Description:
! ------------
!> Writes user specific restart data into binary file(s) for restart runs.
!------------------------------------------------------------------------------!
 MODULE user_write_restart_data_mod


    USE user
     

    IMPLICIT NONE


    INTERFACE user_wrd_global
       MODULE PROCEDURE user_wrd_global
    END INTERFACE user_wrd_global

    INTERFACE user_wrd_local
       MODULE PROCEDURE user_wrd_local
    END INTERFACE user_wrd_local


    PUBLIC user_wrd_global, user_wrd_local


 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Writes global and user-defined restart data into binary file(s) for restart
!> runs.
!------------------------------------------------------------------------------!
    SUBROUTINE user_wrd_global  


       IMPLICIT NONE

       
!       CALL wrd_write_string( 'global_parameter' )
!       WRITE ( 14 )  global_parameter

       
    END SUBROUTINE user_wrd_global    


! Description:
! ------------
!> Writes processor specific and user-defined restart data into binary file(s) 
!> for restart runs.
!------------------------------------------------------------------------------!
    SUBROUTINE user_wrd_local


       IMPLICIT NONE


!
!-- Here the user-defined actions at the end of a job follow.
!-- Sample for user-defined output:
!          IF ( ALLOCATED( u2_av ) )  THEN
!             CALL wrd_write_string( 'u2_av' )  
!             WRITE ( 14 )  u2_av
!          ENDIF



    END SUBROUTINE user_wrd_local


 END MODULE user_write_restart_data_mod