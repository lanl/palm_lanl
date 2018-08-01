!> @file user_dvrp_coltab.f90
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
! $Id: user_dvrp_coltab.f90 2718 2018-01-02 08:49:38Z maronga $
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
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 211 2008-11-11 04:46:24Z raasch
! Former file user_interface.f90 split into one file per subroutine
!
! Description:
! ------------
!> Definition of the colour table to be used by the dvrp software.
!------------------------------------------------------------------------------!
 SUBROUTINE user_dvrp_coltab( mode, variable )
 
    
    USE control_parameters
        
    USE dvrp_variables
    
    USE kinds
    
    USE pegrid
    
    USE user

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  mode       !< 
    CHARACTER (LEN=*) ::  variable   !< 


!
!-- Here the user-defined actions follow
    SELECT CASE ( mode )

       CASE ( 'particles' )

       CASE ( 'slicer' )

       CASE DEFAULT
          message_string = 'unknown mode "' // mode // '"'
          CALL message( 'user_dvrp_coltab', 'UI0004', 1, 2, 0, 6, 0 )


    END SELECT

 END SUBROUTINE user_dvrp_coltab

