!> @file init_ocean.f90
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
! ------------------
! $Id: init_ocean.f90 2846 2018-03-01 08:48:47Z raasch $
! FORTRAN bugfix for r2845
!
! 2845 2018-03-01 08:32:34Z raasch
! bugfix: set kinematic viscosity for sea water
!
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2502 2017-09-26 15:51:59Z gronemeier
! Bugfix: use equation of state for seawater to calculate rho_ocean_init at nzt
!
! 2195 2017-03-23 08:15:17Z raasch
!
! 2194 2017-03-23 08:03:04Z raasch
! bugfix: density is now used in single-value reference state instead of
! potential density
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho_init to rho_ocean_init
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL constants defined as wp_kind
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1179 2013-06-14 05:57:58Z raasch
! Initial density profile is stored in array hom
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 97 2007-06-21 08:23:15Z raasch
! Initial revision
!
! Description:
! ------------
!> Initialization of quantities needed for the ocean version
!------------------------------------------------------------------------------!
 SUBROUTINE init_ocean


    USE arrays_3d,                                                             &
        ONLY:  dzu, dzw, hyp, pt_init, ref_state, rho_ambient, rho_ref_zu,     &
               rho_ref_zw, sa_init, u_init, v_init, zu, zw

    USE cloud_parameters,                                                      &
        ONLY:  cp

    USE constants,                                                             &
        ONLY:  cpw

    USE control_parameters,                                                    &
        ONLY:  ambient_density_for_buoyancy, cos_alpha_surface,                &
               dpdxy, dpdxy_loc, dpdx, dpdy, dpdx_phase, dpdy_phase,           &
               f, g, initialize_to_geostrophic,                                &
               message_string, molecular_viscosity, prandtl_number,            &
               prho_reference, rho_reference, rho_surface, surface_pressure,   &
               use_single_reference_value, stokes_force

    USE eqn_state_seawater_mod,                                                &
        ONLY:  eqn_state_seawater, eqn_state_seawater_func

    USE indices,                                                               &
        ONLY:  nzb, nzt

    USE kinds

    USE pegrid

    USE statistics,                                                            &
        ONLY:  hom, statistic_regions

    USE stokes_drift_mod,                                                      &
        ONLY:  init_stokes_drift

    IMPLICIT NONE

    INTEGER(iwp) ::  k !<
    INTEGER(iwp) ::  n !<

    REAL(wp)     ::  pt_l !<
    REAL(wp)     ::  sa_l !<

    ALLOCATE( hyp(nzb:nzt+1) )

!
!-- Set kinematic viscosity to sea water at 20C
!-- WARNING: This value is especially used for calculating terminal fall
!--          velocities in the LPM, so lpm_init should always be called after
!--          init_ocean!
    molecular_viscosity = 1.05E-6_wp

!
!-- Set specific heat to specific heat of water 
    cp = cpw

!
!-- Calculate initial vertical profile of hydrostatic pressure (in Pa)
!-- and the reference density (used later in buoyancy term)
!-- First step: Calculate pressure using reference density
    hyp(nzt+1) = surface_pressure * 100.0_wp

    rho_surface = eqn_state_seawater_func( hyp(nzt+1), pt_init(nzt+1), sa_init(nzt+1) )
    rho_ref_zw(nzt+1) = rho_surface
    
    hyp(nzt)      = hyp(nzt+1) + rho_surface * g * 0.5_wp * dzu(nzt+1)

    DO  k = nzt-1, nzb+1, -1
       hyp(k) = hyp(k+1) + rho_surface * g * cos_alpha_surface * dzu(k)
    ENDDO
    hyp(nzb) = hyp(nzb+1) + rho_surface * g * cos_alpha_surface * dzu(nzb+1)

!
!-- Second step: Iteratively calculate in situ density (based on presssure)
!-- and pressure (based on in situ density)
    DO  n = 1, 5

       DO  k = nzt, nzb, -1

!--       Calculate initial profiles on the zw grid
          sa_l = 0.5_wp * ( sa_init(k) + sa_init(k+1) ) 
          pt_l = 0.5_wp * ( pt_init(k) + pt_init(k+1) )

          rho_ref_zw(k) = eqn_state_seawater_func(                             &
                                0.5_wp * ( hyp(k) + hyp(k+1) ), pt_l, sa_l )

       ENDDO

       DO  k = nzt, nzb, -1
          hyp(k) = hyp(k+1) + g * cos_alpha_surface *                          &
                         0.5_wp * ( rho_ref_zw(k) + rho_ref_zw(k+1) ) * dzw(k+1)
       ENDDO

    ENDDO

!
!-- Define reference in situ density on the zw grid
    rho_reference = 0.0_wp
    DO  k = nzt, nzb, -1
       sa_l = 0.5_wp * ( sa_init(k) + sa_init(k+1) )
       pt_l = 0.5_wp * ( pt_init(k) + pt_init(k+1) )
       rho_ref_zw(k) = eqn_state_seawater_func( 0.5_wp * ( hyp(k) + hyp(k+1) ),&
                                                pt_l, sa_l )
       rho_reference = rho_reference + rho_ref_zw(k) * dzu(k+1)
       rho_ref_zu(k) = eqn_state_seawater_func( hyp(k), pt_init(k), sa_init(k) )
    ENDDO
    rho_reference = rho_reference / ( zu(nzt+1) - zu(nzb) )

!
!-- Define ambient density profile on zw grid
!-- This ensures that buoyancy is computed relative to the far-field (bottom of 
!-- the domain), hydrostatic, initial density
    IF ( ambient_density_for_buoyancy ) THEN
       DO  k = nzt, nzb, -1
          rho_ambient(k) = eqn_state_seawater_func( 0.5_wp * ( hyp(k) + hyp(k+1) ),&
                                                   pt_init(0), sa_init(0) )
       ENDDO   
    ENDIF
!
!-- Calculate the reference potential density on the zw grid
    prho_reference = 0.0_wp
    DO  k = nzt, nzb, -1
       sa_l = 0.5_wp * ( sa_init(k) + sa_init(k+1) )
       pt_l = 0.5_wp * ( pt_init(k) + pt_init(k+1) )

       prho_reference = prho_reference + dzu(k+1) * &
                        eqn_state_seawater_func( 0.0_wp, pt_l, sa_l )
    ENDDO
    prho_reference = prho_reference / ( zu(nzt+1) - zu(nzb) )

!
!-- Store initial density profile
    hom(:,1,77,:)  = SPREAD( rho_ref_zw(:), 2, statistic_regions+1 )

!
!-- Set the reference state to be used in the buoyancy terms
    IF ( use_single_reference_value )  THEN
       ref_state(:) = rho_reference
    ELSE
       ref_state(:) = rho_ref_zw(:)
    ENDIF

    IF (initialize_to_geostrophic) THEN
       !-- For sinusoidally varying pressure gradient, solve for 
       !-- instantaneous pressure gradient
       dpdxy_loc = dpdxy
       IF ( ANY( dpdx /= 0.0_wp ) .OR. ANY( dpdy /= 0.0_wp) ) THEN
          DO k = 1, 30
             dpdxy_loc(1) = dpdxy_loc(1) + dpdx(k)*cos(-1.0_wp*dpdx_phase(k))
             dpdxy_loc(2) = dpdxy_loc(2) + dpdy(k)*cos(-1.0_wp*dpdy_phase(k))
          ENDDO
       ENDIF
          
       !-- Update u_init and v_init used in rayleigh damping scheme
       u_init(:) = -1.0_wp*dpdxy_loc(2)/(rho_reference*f)
       v_init(:) =         dpdxy_loc(1)/(rho_reference*f)
          
    ENDIF
    
!
!-- Initialize Stokes drift, if required
    IF ( stokes_force ) THEN
       CALL init_stokes_drift
    ENDIF

 END SUBROUTINE init_ocean
