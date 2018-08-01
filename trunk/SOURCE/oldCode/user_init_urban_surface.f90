!> @file user_init_urban_surface.f90
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
! Copyright 2015 Czech Technical University in Prague
! Copyright 1997-2018 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: user_init_urban_surface.f90 2718 2018-01-02 08:49:38Z maronga $
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
! 2011 2016-09-19 17:29:57Z kanani
! Flag urban_surface is now defined in module control_parameters.
! 
! 2007 2016-08-24 15:47:17Z kanani
! Initial revision
!
!
! Description:
! ------------
!> Execution of user-defined actions to initiate the urban surface model
!------------------------------------------------------------------------------!
 SUBROUTINE user_init_urban_surface

    USE arrays_3d
    
    USE control_parameters,                                                    &
        ONLY:  urban_surface
    
    USE indices
    
    USE kinds

    USE urban_surface_mod

    USE surface_mod    
    
    USE user

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< grid index
    INTEGER(iwp) ::  j  !< grid index
    INTEGER(iwp) ::  m  !< running index on 1D wall-type grid

!
!-- Here the user-defined urban surface initialization actions follow.
!-- Example: set roughness length at urban surface
!     DO  m = 1, surf_usm_h%ns
!        surf_usm_h%z0(m) = 0.1_wp
!     ENDDO



 END SUBROUTINE user_init_urban_surface

