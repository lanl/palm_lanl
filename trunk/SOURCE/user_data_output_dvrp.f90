!> @file user_data_output_dvrp.f90
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
! $Id: user_data_output_dvrp.f90 3049 2018-05-29 13:52:36Z Giersch $
! Error messages revised
! 
! 2718 2018-01-02 08:49:38Z Giersch
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
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! old module precision_kind is removed,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 211 2008-11-11 04:46:24Z raasch
! Former file user_interface.f90 split into one file per subroutine
!
! Description:
! ------------
!> Execution of user-defined dvrp output
!------------------------------------------------------------------------------!
 SUBROUTINE user_data_output_dvrp( output_variable, local_pf )
 

    USE control_parameters

    USE dvrp_variables

    USE indices

    USE kinds

    USE pegrid

    USE user

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  output_variable !< 

    INTEGER(iwp) ::  i !< 
    INTEGER(iwp) ::  j !< 
    INTEGER(iwp) ::  k !< 

    REAL(wp), DIMENSION(nxl_dvrp:nxr_dvrp+1,nys_dvrp:nyn_dvrp+1,nzb:nz_do3d) :: &
              local_pf !< 

!
!-- Here the user-defined DVRP output follows:

!
!-- Move original array to intermediate array
    SELECT CASE ( output_variable )

!       CASE ( 'u2', 'u2_xy', 'u2_xz', 'u2_yz'  )
!!
!!--       Here the user can add user_defined output quantities. 
!!--       Uncomment and extend the following lines, if necessary.
!          DO  i = nxl_dvrp, nxr_dvrp+1
!             DO  j = nys_dvrp, nyn_dvrp+1
!                DO  k = nzb, nz_do3d
!                   local_pf(i,j,k) = u2(k,j,i)
!                ENDDO
!             ENDDO
!          ENDDO


       CASE DEFAULT
!
!--       The DEFAULT case is reached if output_variable contains a
!--       wrong character string that is neither recognized in data_output_dvrp
!--       nor here in user_data_output_dvrp.
          WRITE( message_string, * ) 'no output possible for: ',               &
                                     output_variable
          CALL message( 'user_data_output_dvrp', 'UI0003', 0, 1, 0, 6, 0 )
          

    END SELECT


 END SUBROUTINE user_data_output_dvrp

