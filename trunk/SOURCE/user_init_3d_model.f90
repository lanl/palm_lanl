!> @file user_init_3d_model.f90
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
! $Id: user_init_3d_model.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2618 2017-11-16 15:37:30Z suehring
! Provide example for user-defined initialization of surface-related quantities
! 
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! +surface_mod
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
!
! 1320 2014-03-20 08:40:49Z raasch
! small changes in layout
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 211 2008-11-11 04:46:24Z raasch
! Former file user_interface.f90 split into one file per subroutine
!
! Description:
! ------------
!> Allows the complete initialization of the 3d model.
!>
!> @attention The user is responsible to set at least all those quantities which
!>            are normally set within init_3d_model!
!------------------------------------------------------------------------------!
 SUBROUTINE user_init_3d_model
 

    USE arrays_3d
    
    USE control_parameters
    
    USE indices
    
    USE kinds

    USE surface_mod
    
    USE user

    IMPLICIT NONE

    INTEGER(iwp) ::  l !< running index surface orientation
    INTEGER(iwp) ::  m !< running index surface elements

!
!-- Initialization of surface-related quantities.
!-- The following example shows required initialization of surface quantitites
!-- at default-type upward-facing surfaces.  
!   DO  m = 1, surf_def_h(0)%ns
!      surf_def_h(0)%ol(m)   = ...    ! Obukhov length
!      surf_def_h(0)%us(m  ) = ...    ! friction velocity
!      surf_def_h(0)%usws(m) = ...    ! vertical momentum flux, u-component
!      surf_def_h(0)%vsws(m) = ...    ! vertical momentum flux, v-component
!      surf_def_h(0)%z0(m)   = ...    ! roughness length for momentum
!      IF ( .NOT. neutral )  THEN
!         surf_def_h(0)%ts(m)   = ... ! scaling parameter
!         surf_def_h(0)%shf(m)  = ... ! surface sensible heat flux
!         surf_def_h(0)%z0h(m)  = ... ! roughness length for heat
!      ENDIF
!      IF ( humditiy )  THEN
!         surf_def_h(0)%qs(m)   = ... ! scaling parameter
!         surf_def_h(0)%qsws(m) = ... ! surface latent heat flux
!         surf_def_h(0)%z0q(m)  = ... ! roughness length for moisture
!      ENDIF
!      IF ( passive_scalar )  THEN
!         surf_def_h(0)%ss(m)   = ... ! scaling parameter
!         surf_def_h(0)%ssws(m) = ... ! surface latent heat flux
!      ENDIF
!   ENDDO 
!
!-- Same for natural and urban type surfaces
!   DO  m = 1, surf_lsm_h%ns
!      ...
!   ENDDO 
!   DO  m = 1, surf_usm_h%ns
!      ...
!   ENDDO
!
!-- Also care for vertically aligned surfaces (default-, natural-, and 
!-- urban-type).
!   DO  l = 0, 3
!      DO  m = 1, surf_def_v(l)%ns
!         ...
!      ENDDO
!      DO  m = 1, surf_lsm_v(l)%ns
!         ...
!      ENDDO
!      DO  m = 1, surf_usm_v(l)%ns
!         ...
!      ENDDO
!   ENDDO
!
!
!-- In the following, initialize 3D quantities, e.g. u, v, w, pt, etc..

 END SUBROUTINE user_init_3d_model

