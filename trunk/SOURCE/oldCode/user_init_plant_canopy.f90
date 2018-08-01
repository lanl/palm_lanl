!> @file user_init_plant_canopy.f90
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
! $Id: user_init_plant_canopy.f90 2718 2018-01-02 08:49:38Z maronga $
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
! 1951 2016-06-20 10:05:12Z suehring
! Bugfix in example initialization
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1484 2014-10-21 10:53:05Z kanani
! Changes in the course of the canopy-model modularization:
!   module plant_canopy_model_mod added,
!   definition of array cdc (canopy drag coefficient) removed, since it is now
!   defined purely as a single constant value (see module plant_canopy_model)
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute 
!
! 1320 2014-03-20 08:40:49Z raasch
! kind-parameters added to all INTEGER and REAL declaration statements, 
! kinds are defined in new module kinds, 
! old module precision_kind is removed,
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
!> Initialisation of the leaf area density array (for scalar grid points) and 
!> the array of the canopy drag coefficient, if the user has not chosen any 
!> of the default cases 
!------------------------------------------------------------------------------!
 SUBROUTINE user_init_plant_canopy
 

    USE arrays_3d
    
    USE control_parameters
    
    USE indices
    
    USE kinds

    USE plant_canopy_model_mod
    
    USE user

    IMPLICIT NONE

    INTEGER(iwp) :: i   !< running index
    INTEGER(iwp) :: j   !< running index

!
!-- Here the user-defined grid initializing actions follow:

!
!-- Set the 3D-array lad_s for user defined canopies
    SELECT CASE ( TRIM( canopy_mode ) )

       CASE ( 'block' )
!
!--       Not allowed here since this is the standard case used in init_3d_model.

       CASE ( 'user_defined_canopy_1' )
!
!--       Here the user can define his own forest topography. 
!--       The following lines contain an example, where the plant canopy extends
!--       only over the second half of the model domain along x.
!--       Attention: DO-loops have to include the ghost points (nxlg-nxrg, 
!--       nysg-nyng), because no exchange of ghost point information is intended,
!--       in order to minimize communication between CPUs.
!          DO  i = nxlg, nxrg
!             IF ( i >= INT( ( nx+1 ) / 2 ) ) THEN
!                DO  j = nysg, nyng
!                   lad_s(:,j,i) = lad(:)
!                ENDDO
!             ELSE
!                lad_s(:,:,i) = 0.0_wp
!             ENDIF
!          ENDDO 
!
!--       After definition, please
!--       remove the following three lines!
          message_string = 'canopy_mode "' // canopy_mode // &
                           '" not available yet'
          CALL message( 'user_init_plant_canopy', 'UI0007', 0, 1, 0, 6, 0 )
          
       CASE DEFAULT
!
!--       The DEFAULT case is reached if the parameter canopy_mode contains a
!--       wrong character string that is neither recognized in init_3d_model nor
!--       here in user_init_plant_canopy.
          message_string = 'unknown canopy_mode "' // canopy_mode // '"'
          CALL message( 'user_init_plant_canopy', 'UI0008', 1, 2, 0, 6, 0 )

    END SELECT


 END SUBROUTINE user_init_plant_canopy

