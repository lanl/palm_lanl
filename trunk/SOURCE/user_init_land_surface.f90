!> @file user_init_land_surface.f90
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
! $Id: user_init_land_surface.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Load surface_mod
! Add simple example how to access surface data type
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1972 2016-07-26 07:52:02Z maronga
! Update of use statements
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1585 2015-04-30 07:05:52Z maronga
! Changed description text
! 
! 1496 2014-12-02 17:25:50Z maronga
! Initial revision
! 
! Description:
! ------------
!> Execution of user-defined actions to initiate the land surface model
!------------------------------------------------------------------------------!
 SUBROUTINE user_init_land_surface
 

    USE control_parameters
    
    USE indices
    
    USE kinds
    
    USE land_surface_model_mod

    USE netcdf_interface,                                                      &
        ONLY: dots_label, dots_unit, dots_num
    
    USE pegrid

    USE surface_mod    

    USE user

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< grid index
    INTEGER(iwp) ::  j  !< grid index
    INTEGER(iwp) ::  m  !< running index on 1D wall-type grid

!
!-- Here the user-defined land surface initialization actions follow.
!-- Example: set roughness length at natural land-surface
!     DO  m = 1, surf_lsm_h%ns
!        surf_lsm_h%z0(m) = 0.1_wp
!     ENDDO


 END SUBROUTINE user_init_land_surface

