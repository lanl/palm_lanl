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
        ONLY:  dzu, hyp, pt_init, ref_state, sa_init, zu, zw

    USE control_parameters,                                                    &
        ONLY:  g, molecular_viscosity, prho_reference, rho_surface,            &
               rho_reference, surface_pressure, use_single_reference_value,    &
               stokes_force

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

    REAL(wp), DIMENSION(nzb:nzt+1) ::  rho_ocean_init !<

    ALLOCATE( hyp(nzb:nzt+1) )

!
!-- Set water density near the ocean surface
    rho_surface = 1027.62_wp

!
!-- Set kinematic viscosity to sea water at 20C
!-- WARNING: This value is especially used for calculating terminal fall
!--          velocities in the LPM, so lpm_init should always be called after
!--          init_ocean!
    molecular_viscosity = 1.05E-6_wp

!
!-- Calculate initial vertical profile of hydrostatic pressure (in Pa)
!-- and the reference density (used later in buoyancy term)
!-- First step: Calculate pressure using reference density
    hyp(nzt+1) = surface_pressure * 100.0_wp

    hyp(nzt)      = hyp(nzt+1) + rho_surface * g * 0.5_wp * dzu(nzt+1)
    rho_ocean_init(nzt+1) = rho_surface

    DO  k = nzt-1, 1, -1
       hyp(k) = hyp(k+1) + rho_surface * g * dzu(k)
    ENDDO
    hyp(0) = hyp(1) + rho_surface * g * dzu(1)

!
!-- Second step: Iteratively calculate in situ density (based on presssure)
!-- and pressure (based on in situ density)
    DO  n = 1, 5

       rho_reference = rho_surface * 0.5_wp * dzu(nzt+1)

       DO  k = nzt, 0, -1

          sa_l = 0.5_wp * ( sa_init(k) + sa_init(k+1) )
          pt_l = 0.5_wp * ( pt_init(k) + pt_init(k+1) )

          rho_ocean_init(k) = eqn_state_seawater_func( hyp(k), pt_l, sa_l )

          rho_reference = rho_reference + rho_ocean_init(k) * dzu(k+1)

       ENDDO

       rho_reference = rho_reference / ( zw(nzt) - zu(nzb) )

    
       DO  k = nzt, 0, -1
          hyp(k) = hyp(k+1) + g * 0.5_wp * ( rho_ocean_init(k)                 &
                                           + rho_ocean_init(k+1) ) * dzu(k+1)
       ENDDO

    ENDDO

!
!-- Calculate the reference potential density
    prho_reference = 0.0_wp
    DO  k = 0, nzt

       sa_l = 0.5_wp * ( sa_init(k) + sa_init(k+1) )
       pt_l = 0.5_wp * ( pt_init(k) + pt_init(k+1) )

       prho_reference = prho_reference + dzu(k+1) * &
                        eqn_state_seawater_func( 0.0_wp, pt_l, sa_l )

    ENDDO

    prho_reference = prho_reference / ( zu(nzt) - zu(nzb) )


!
!-- Calculate the 3d array of initial in situ and potential density,
!-- based on the initial temperature and salinity profile
    CALL eqn_state_seawater

!
!-- Store initial density profile
    hom(:,1,77,:)  = SPREAD( rho_ocean_init(:), 2, statistic_regions+1 )

!
!-- Set the reference state to be used in the buoyancy terms
    IF ( use_single_reference_value )  THEN
       ref_state(:) = rho_reference
    ELSE
       ref_state(:) = rho_ocean_init(:)
    ENDIF

!
!-- Initialize Stokes drift, if required
    IF ( stokes_force ) THEN
       CALL init_stokes_drift
    ENDIF

 END SUBROUTINE init_ocean
