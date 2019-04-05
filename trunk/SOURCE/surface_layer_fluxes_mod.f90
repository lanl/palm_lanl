!> @file surface_layer_fluxes_mod.f90
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
!
!------------------------------------------------------------------------------!
!
! Current revisions:
! ------------------
! 
! 2018-11-21 cbegeman
! Add McPhee MOST parameterization of surface layer fluxes for ice-ocean 
! interface
! 
! Former revisions:
! -----------------
! $Id: surface_layer_fluxes_mod.f90 3045 2018-05-28 07:55:41Z Giersch $
! Error message revised
! 
! 2766 2018-01-22 17:17:47Z kanani
! Removed preprocessor directive __chem
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! - Change in file header (GPL part)
! - Implementation of chemistry module (FK)
! - Added calculation of pt1 and qv1 for all surface types. Added calculation of
!   pt_surface for default-type surfaces (BM)
! - Add flag to disable computation of qsws in case of urban surface (MS)
! 
! 2547 2017-10-16 12:41:56Z schwenkel
! extended by cloud_droplets option
! 
! 2321 2017-07-24 15:57:07Z schwenkel
! Bugfix: Correct index in lookup table for Obukhov length 
! 
! 2299 2017-06-29 10:14:38Z suehring
! Adjusted for allow separate spinups of LSM and atmosphere code
! 
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison' 
! includes two more prognostic equations for cloud drop concentration (nc)  
! and cloud water content (qc). 
! 
! 2281 2017-06-13 11:34:50Z suehring
! Clean-up unnecessary index access to surface type
! 
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new surface concept
! OpenMP bugfix
! 
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC directives and related code removed
! 
! 2091 2016-12-21 16:38:18Z suehring
! Bugfix in calculation of vsws ( incorrect linear interpolation of us ) 
! 
! 2076 2016-12-02 13:54:20Z raasch
! further openmp bugfix for lookup method
! 
! 2073 2016-11-30 14:34:05Z raasch
! openmp bugfix for lookup method
! 
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
! 
! 2011 2016-09-19 17:29:57Z kanani
! Flag urban_surface is now defined in module control_parameters.
! 
! 2007 2016-08-24 15:47:17Z kanani
! Account for urban surface model in computation of vertical kinematic heatflux
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1992 2016-08-12 15:14:59Z suehring
! Minor bug, declaration of look-up index as INTEGER
! 
! 1960 2016-07-12 16:34:24Z suehring
! Treat humidity and passive scalar separately
! 
! 1929 2016-06-09 16:25:25Z suehring
! Bugfix: avoid segmentation fault in case one grid point is horizontally 
! completely surrounded by topography
! 
! 1920 2016-05-30 10:50:15Z suehring
! Avoid segmentation fault (see change in 1915) by different initialization of 
! us instead of adding a very small number in the denominator
!
! 1915 2016-05-27 11:05:02Z suehring
! Bugfix: avoid segmentation fault in case of most_method = 'circular' at first
! timestep
!
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
!
! 
! 1822 2016-04-07 07:49:42Z hoffmann
! icloud_scheme replaced by microphysics_*
!
! 1788 2016-03-10 11:01:04Z maronga
! Added parameter z0q which replaces z0h in the similarity functions for
! humidity.
! Syntax layout improved.
! 
! 1757 2016-02-22 15:49:32Z maronga
! Minor fixes.
!
! 1749 2016-02-09 12:19:56Z raasch
! further OpenACC adjustments
!
! 1747 2016-02-08 12:25:53Z raasch
! adjustments for OpenACC usage
!
! 1709 2015-11-04 14:47:01Z maronga
! Bugfix: division by zero could occur when calculating rib at low wind speeds
! Bugfix: calculation of uv_total for neutral = .T., initial value for ol for
! neutral = .T.
! 
! 1705 2015-11-02 14:28:56Z maronga
! Typo removed
!
! 1697 2015-10-28 17:14:10Z raasch
! FORTRAN and OpenMP errors removed
!
! 1696 2015-10-27 10:03:34Z maronga
! Modularized and completely re-written version of prandtl_fluxes.f90. In the
! course of the re-writing two additional methods have been implemented. See
! updated description.
!
! 1551 2015-03-03 14:18:16Z maronga
! Removed land surface model part. The surface fluxes are now always calculated
! within prandtl_fluxes, based on the given surface temperature/humidity (which 
! is either provided by the land surface model, by large scale forcing data, or
! directly prescribed by the user.
! 
! 1496 2014-12-02 17:25:50Z maronga
! Adapted for land surface model
! 
! 1494 2014-11-21 17:14:03Z maronga
! Bugfixes: qs is now calculated before calculation of Rif. calculation of
! buoyancy flux in Rif corrected (added missing humidity term), allow use of 
! topography for coupled runs (not tested)
! 
! 1361 2014-04-16 15:17:48Z hoffmann
! Bugfix: calculation of turbulent fluxes of rain water content (qrsws) and rain 
! drop concentration (nrsws) added
! 
! 1340 2014-03-25 19:45:13Z kanani
! REAL constants defined as wp-kind
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! old module precision_kind is removed,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1276 2014-01-15 13:40:41Z heinze
! Use LSF_DATA also in case of Dirichlet bottom boundary condition for scalars
!
! 1257 2013-11-08 15:18:40Z raasch
! openACC "kernels do" replaced by "kernels loop", "loop independent" added
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1015 2012-09-27 09:23:24Z raasch
! OpenACC statements added
!
! 978 2012-08-09 08:28:32Z fricke
! roughness length for scalar quantities z0h added
!
! Revision 1.1  1998/01/23 10:06:06  raasch
! Initial revision
!
!
! Description:
! ------------
!> Diagnostic computation of vertical fluxes in the constant flux layer from the
!> values of the variables at grid point k=1. Three different methods are
!> available:
!> 1) the "old" version (most_method = 'circular') which is fast, but inaccurate
!> 2) a Newton iteration method (most_method = 'newton'), which is accurate, but
!>    slower
!> 3) a method using a lookup table which is fast and accurate. Note, however,
!>    that this method cannot be used in case of roughness heterogeneity
!>
!> @todo (re)move large_scale_forcing actions
!> @todo check/optimize OpenMP directives
!> @todo simplify if conditions (which flux need to be computed in which case)
!------------------------------------------------------------------------------!
 MODULE surface_layer_fluxes_mod

    USE arrays_3d,                                                             &
        ONLY:  alpha_T, beta_S, e, hyp, kh, nc, nr, pt, q, ql, qc, qr,         &
               rho_ocean, s, sa, u, v, vpt, w, zu, zw,                         &
               drho_ref_zw, rho_ref_zw

    USE chem_modules,                                                          &
        ONLY:  constant_csflux, nvar

    USE cloud_parameters,                                                      &
        ONLY:  cp, l_d_cp, pt_d_t

    USE constants,                                                             &
        ONLY:  pi, cpw

    USE cpulog

    USE control_parameters,                                                    &
        ONLY:  air_chemistry, c1, c2, c3, cloud_droplets, cloud_physics,       &
               constant_flux_layer, constant_heatflux, constant_scalarflux,    &     
               constant_waterflux, coupling_mode, drag_coeff, f, g, humidity,  &
               ibc_e_b, ibc_e_t, ibc_pt_b, ibc_pt_t, initializing_actions,     &
               intermediate_timestep_count, intermediate_timestep_count_max,   &
               k_offset_mcphee, kappa, l_m, land_surface, large_scale_forcing, &
               lsf_surf,                                                       &
               message_string, microphysics_morrison, microphysics_seifert,    &
               molecular_viscosity, most_method, neutral, ocean,               &
               passive_scalar, prandtl_number, pt_surface, q_surface,          &
               run_coupled, schmidt_number, surface_pressure, simulated_time,  &
               terminate_run, time_since_reference_point, urban_surface,       &
               z_offset_mcphee, zeta_max, zeta_min

    USE eqn_state_seawater_mod,                                                &
        ONLY: pt_freezing, pt_freezing_SA

    USE grid_variables,                                                        &
        ONLY:  dx, dy  

    USE indices,                                                               &
        ONLY:  nxl, nxr, nys, nyn, nzb, nzt

    USE kinds

    USE pegrid

    USE land_surface_model_mod,                                                &
        ONLY:  aero_resist_kray, skip_time_do_lsm

    USE surface_mod,                                                           &
        ONLY :  surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v, surf_type,     &
                surf_usm_h, surf_usm_v
        

    IMPLICIT NONE

    INTEGER(iwp) ::  i              !< loop index x direction
    INTEGER(iwp) ::  j              !< loop index y direction
    INTEGER(iwp) ::  k              !< loop index z direction
    INTEGER(iwp) ::  l              !< loop index for surf type
    INTEGER(iwp) ::  li_bnd  = 7500 !< Lookup table index of the last time step
    INTEGER(iwp) ::  m              !< loop index surface

    INTEGER(iwp), PARAMETER ::  num_steps = 15000  !< number of steps in the lookup table

    LOGICAL      ::  coupled_run       !< Flag for coupled atmosphere-ocean runs
    LOGICAL      ::  downward = .FALSE.!< Flag indicating downward-facing horizontal surface
    LOGICAL      ::  mom_uv  = .FALSE. !< Flag indicating calculation of usvs and vsus at vertical surfaces
    LOGICAL      ::  mom_w   = .FALSE. !< Flag indicating calculation of wsus and wsvs at vertical surfaces
    LOGICAL      ::  mom_tke = .FALSE. !< Flag indicating calculation of momentum fluxes at vertical surfaces used for TKE production 
    LOGICAL      ::  surf_vertical     !< Flag indicating vertical surfaces

    REAL(wp), DIMENSION(0:num_steps-1) :: rib_tab,  & !< Lookup table bulk Richardson number
                                          ol_tab      !< Lookup table values of L

    REAL(wp)     ::  dheatflux,                & !< Change in heat flux over iteration
                     dheatflux_tol = 1.0E-6_wp,& !< Maximum acceptable change in heat flux over iteration !CB
                     eta_star,                 & !< Stability factor 
                     e_s,                      & !< Saturation water vapor pressure
                     ol_max = 1.0E6_wp,        & !< Maximum Obukhov length
                     rib_max,                  & !< Maximum Richardson number in lookup table
                     rib_min,                  & !< Minimum Richardson number in lookup table
                     z_mo                        !< Height of the constant flux layer where MOST is assumed
    REAL(wp)     ::  ri_crit = 0.2_wp            !< critical flux Richardson number
    REAL(wp)     ::  xi_N = 0.0052_wp            !< non-dimensional surface layer extent

    TYPE(surf_type), POINTER ::  surf     !< surf-type array, used to generalize subroutines


    SAVE

    PRIVATE

    PUBLIC init_surface_layer_fluxes, surface_layer_fluxes

    INTERFACE init_surface_layer_fluxes
       MODULE PROCEDURE init_surface_layer_fluxes
    END INTERFACE init_surface_layer_fluxes

    INTERFACE surface_layer_fluxes
       MODULE PROCEDURE surface_layer_fluxes
    END INTERFACE surface_layer_fluxes


 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Main routine to compute the surface fluxes
!------------------------------------------------------------------------------!
    SUBROUTINE surface_layer_fluxes

       IMPLICIT NONE

       !CB INTEGER(iwp) :: m
       REAL(wp), DIMENSION(:), ALLOCATABLE :: shf_p
 
       surf_vertical = .FALSE.
       downward      = .FALSE.

!
!--    Derive potential temperature and specific humidity at first grid level 
!--    from the fields pt and q
       IF ( TRIM(constant_flux_layer) == 'bottom' )  THEN
!
!--       First call for horizontal default-type surfaces (l=0 - upward facing, 
!--       l=1 - downward facing)
          DO  l = 0, 1
             IF ( surf_def_h(l)%ns >= 1 )  THEN
                surf => surf_def_h(l)
                CALL calc_pt_q
                IF ( .NOT. neutral )  CALL calc_pt_surface
             ENDIF
          ENDDO
!
!--       Call for natural-type horizontal surfaces
          IF ( surf_lsm_h%ns >= 1 )  THEN
             surf => surf_lsm_h
             CALL calc_pt_q
          ENDIF

!
!--       Call for urban-type horizontal surfaces
          IF ( surf_usm_h%ns >= 1 )  THEN
             surf => surf_usm_h
             CALL calc_pt_q
          ENDIF

!
!--       Call for natural-type vertical surfaces
          DO  l = 0, 3
             IF ( surf_lsm_v(l)%ns >= 1 )  THEN
                surf => surf_lsm_v(l)
                CALL calc_pt_q
             ENDIF

!
!--       Call for urban-type vertical surfaces
             IF ( surf_usm_v(l)%ns >= 1 )  THEN
                surf => surf_usm_v(l)
                CALL calc_pt_q
             ENDIF
          ENDDO
       
       ENDIF

       IF ( TRIM(constant_flux_layer) == 'top' )  THEN
!
!--       First call for horizontal default-type surfaces 
          IF ( surf_def_h(2)%ns >= 1 )  THEN
             surf => surf_def_h(2)
!
!--          Derive scalar properties at the nth grid level 
!--          from the surface from their 3-d arrays
             IF ( ocean .AND. trim(most_method) == 'mcphee' ) THEN
                surf => surf_def_h(2)
                CALL calc_pt_sa
             ELSE
                CALL calc_pt_q
                IF ( .NOT. neutral )  CALL calc_pt_surface
             ENDIF
          ENDIF
       ENDIF

!
!--    First, calculate the new Obukhov length, then new friction velocity,
!--    followed by the new scaling parameters (th*, q*, etc.), and the new
!--    surface fluxes if required. The old routine ("circular") requires a 
!--    different order of calls as the scaling parameters from the previous time
!--    steps are used to calculate the Obukhov length

!
!--    Depending on setting of most_method use the "old" routine
!--    Note, each routine is called for different surface types. 
!--    First call for default-type horizontal surfaces, for natural- and 
!--    urban-type surfaces. Note, at this place only upward-facing horizontal
!--    surfaces are treted.
       IF ( TRIM(constant_flux_layer) == 'bottom' ) THEN 
          
          IF ( most_method == 'circular' )  THEN
!
!--          Default-type upward-facing horizontal surfaces
             IF ( surf_def_h(0)%ns >= 1 )  THEN
                surf => surf_def_h(0)
                CALL calc_scaling_parameters
                CALL calc_uvw_abs
                IF ( .NOT. neutral )  CALL calc_ol
                CALL calc_us
                CALL calc_surface_fluxes
             ENDIF
!
!--          Natural-type horizontal surfaces
             IF ( surf_lsm_h%ns >= 1 )  THEN
                surf => surf_lsm_h
                CALL calc_scaling_parameters
                CALL calc_uvw_abs
                IF ( .NOT. neutral )  CALL calc_ol
                CALL calc_us
                CALL calc_surface_fluxes
             ENDIF
!
!--          Urban-type horizontal surfaces
             IF ( surf_usm_h%ns >= 1 )  THEN
                surf => surf_usm_h
                CALL calc_scaling_parameters
                CALL calc_uvw_abs
                IF ( .NOT. neutral )  CALL calc_ol
                CALL calc_us
                CALL calc_surface_fluxes
             ENDIF
!
!--       Use either Newton iteration or a lookup table for the bulk Richardson
!--       number to calculate the Obukhov length 
          ELSEIF ( most_method == 'newton'  .OR.  most_method == 'lookup' )  THEN
!
!--          Default-type upward-facing horizontal surfaces
             IF ( surf_def_h(0)%ns >= 1 )  THEN
                surf => surf_def_h(0)
                CALL calc_uvw_abs
                IF ( .NOT. neutral )  CALL calc_ol
                CALL calc_us
                CALL calc_scaling_parameters
                CALL calc_surface_fluxes
             ENDIF
!
!--          Natural-type horizontal surfaces
             IF ( surf_lsm_h%ns >= 1 )  THEN
                surf => surf_lsm_h
                CALL calc_uvw_abs
                IF ( .NOT. neutral )  CALL calc_ol
                CALL calc_us
                CALL calc_scaling_parameters
                CALL calc_surface_fluxes
             ENDIF
!
!--          Urban-type horizontal surfaces
             IF ( surf_usm_h%ns >= 1 )  THEN
                surf => surf_usm_h
                CALL calc_uvw_abs
                IF ( .NOT. neutral )  CALL calc_ol
                CALL calc_us
                CALL calc_scaling_parameters
                CALL calc_surface_fluxes
             ENDIF

          ENDIF
!
!--       Treat downward-facing horizontal surfaces. Note, so far, these are 
!--       always default type. Stratification is not considered
!--       in this case, hence, no further distinction between different 
!--       most_method is required.  
          IF ( surf_def_h(1)%ns >= 1 )  THEN
             downward = .TRUE.
             surf => surf_def_h(1)
             CALL calc_uvw_abs
             CALL calc_us
             CALL calc_surface_fluxes
             downward = .FALSE.
          ENDIF
       ENDIF
      
!
!--    Treat downward-facing horizontal surface at the top.
       IF ( surf_def_h(2)%ns >= 1 .AND. TRIM(constant_flux_layer) == 'top' )  THEN
          
          downward = .TRUE.
          surf => surf_def_h(2)

          IF ( trim(most_method) == 'mcphee' ) THEN
             CALL calc_uvw_abs 
             CALL calc_us
             CALL calc_3eqn
          ELSE
             CALL calc_scaling_parameters
             CALL calc_uvw_abs
             IF ( .NOT. neutral )  CALL calc_ol
             CALL calc_us
             CALL calc_surface_fluxes
          ENDIF
          
          downward = .FALSE.
       
       ENDIF
!
!--    Calculate surfaces fluxes at vertical surfaces for momentum 
!--    and subgrid-scale TKE.
!--    No stability is considered. Therefore, scaling parameters and Obukhov-
!--    length do not need to be calculated and no distinction in 'circular',
!--    'Newton' or 'lookup' is necessary so far. 
!--    Note, this will change if stability is once considered.
       surf_vertical = .TRUE.
!
!--    Calculate horizontal momentum fluxes at north- and south-facing 
!--    surfaces(usvs).
!--    For default-type surfaces
       mom_uv = .TRUE. 
       DO  l = 0, 1
          IF ( surf_def_v(l)%ns >= 1 .AND. TRIM(constant_flux_layer) == 'bottom' )  THEN
             surf => surf_def_v(l)
!
!--          Compute surface-parallel velocity
             CALL calc_uvw_abs_v_ugrid
!
!--          Compute respective friction velocity on staggered grid
             CALL calc_us
!
!--          Compute respective surface fluxes for momentum and TKE
             CALL calc_surface_fluxes
          ENDIF
       ENDDO
!
!--    For natural-type surfaces. Please note, even though stability is not
!--    considered for the calculation of momentum fluxes at vertical surfaces,
!--    scaling parameters and Obukov length are calculated nevertheless in this
!--    case. This is due to the requirement of ts in parameterization of heat
!--    flux in land-surface model in case of aero_resist_kray is not true.
       IF ( TRIM(constant_flux_layer) == 'bottom' ) THEN
          IF ( .NOT. aero_resist_kray )  THEN
             IF ( most_method == 'circular' )  THEN
                DO  l = 0, 1
                   IF ( surf_lsm_v(l)%ns >= 1 )  THEN
                      surf => surf_lsm_v(l)
!
!--                   Compute scaling parameters
                      CALL calc_scaling_parameters
!
!--                   Compute surface-parallel velocity
                      CALL calc_uvw_abs_v_ugrid
!
!--                   Compute Obukhov length
                      IF ( .NOT. neutral )  CALL calc_ol
!
!--                   Compute respective friction velocity on staggered grid
                      CALL calc_us
!
!--                   Compute respective surface fluxes for momentum and TKE
                      CALL calc_surface_fluxes
                   ENDIF
                ENDDO
             ELSE
                DO  l = 0, 1
                   IF ( surf_lsm_v(l)%ns >= 1 )  THEN
                      surf => surf_lsm_v(l)
!
!--                   Compute surface-parallel velocity
                      CALL calc_uvw_abs_v_ugrid
!
!--                   Compute Obukhov length
                      IF ( .NOT. neutral )  CALL calc_ol
!
!--                   Compute respective friction velocity on staggered grid
                      CALL calc_us
!
!--                   Compute scaling parameters
                      CALL calc_scaling_parameters
!
!--                   Compute respective surface fluxes for momentum and TKE
                      CALL calc_surface_fluxes
                   ENDIF
                ENDDO
             ENDIF
!
!--       No ts is required, so scaling parameters and Obukhov length do not need
!--       to be computed.
          ELSE
             DO  l = 0, 1
                IF ( surf_lsm_v(l)%ns >= 1 )  THEN
                   surf => surf_lsm_v(l)
!   
!--                Compute surface-parallel velocity
                   CALL calc_uvw_abs_v_ugrid
!
!--                Compute respective friction velocity on staggered grid
                   CALL calc_us
!
!--                Compute respective surface fluxes for momentum and TKE
                   CALL calc_surface_fluxes
                ENDIF
             ENDDO
          ENDIF
!
!--       For urban-type surfaces
          DO  l = 0, 1
             IF ( surf_usm_v(l)%ns >= 1 )  THEN
                surf => surf_usm_v(l)
!
!--             Compute surface-parallel velocity
                CALL calc_uvw_abs_v_ugrid
!
!--             Compute respective friction velocity on staggered grid
                CALL calc_us
!
!--             Compute respective surface fluxes for momentum and TKE
                CALL calc_surface_fluxes
             ENDIF
          ENDDO
!
!--       Calculate horizontal momentum fluxes at east- and west-facing 
!--       surfaces (vsus).
!--       For default-type surfaces
          DO  l = 2, 3
             IF ( surf_def_v(l)%ns >= 1 )  THEN
                surf => surf_def_v(l)
!
!--             Compute surface-parallel velocity
                CALL calc_uvw_abs_v_vgrid
!
!--             Compute respective friction velocity on staggered grid
                CALL calc_us
!
!--             Compute respective surface fluxes for momentum and TKE
                CALL calc_surface_fluxes
             ENDIF
          ENDDO
!
!--       For natural-type surfaces. Please note, even though stability is not
!--       considered for the calculation of momentum fluxes at vertical surfaces,
!--       scaling parameters and Obukov length are calculated nevertheless in this
!--       case. This is due to the requirement of ts in parameterization of heat
!--       flux in land-surface model in case of aero_resist_kray is not true.
          IF ( .NOT. aero_resist_kray )  THEN
             IF ( most_method == 'circular' )  THEN
                DO  l = 2, 3
                   IF ( surf_lsm_v(l)%ns >= 1 )  THEN
                      surf => surf_lsm_v(l)
!
!--                   Compute scaling parameters
                      CALL calc_scaling_parameters
!
!--                   Compute surface-parallel velocity
                      CALL calc_uvw_abs_v_vgrid
!
!--                   Compute Obukhov length
                      IF ( .NOT. neutral )  CALL calc_ol
!
!--                   Compute respective friction velocity on staggered grid
                      CALL calc_us
!
!--                   Compute respective surface fluxes for momentum and TKE
                      CALL calc_surface_fluxes
                   ENDIF
                ENDDO
             ELSE
                DO  l = 2, 3
                   IF ( surf_lsm_v(l)%ns >= 1 )  THEN
                      surf => surf_lsm_v(l)
!
!--                   Compute surface-parallel velocity
                      CALL calc_uvw_abs_v_vgrid
!
!--                   Compute Obukhov length
                      IF ( .NOT. neutral )  CALL calc_ol
!
!--                   Compute respective friction velocity on staggered grid
                      CALL calc_us
!
!--                   Compute scaling parameters
                      CALL calc_scaling_parameters
!
!--                   Compute respective surface fluxes for momentum and TKE
                      CALL calc_surface_fluxes
                   ENDIF
                ENDDO
             ENDIF
          ELSE
             DO  l = 2, 3
                IF ( surf_lsm_v(l)%ns >= 1 )  THEN
                   surf => surf_lsm_v(l)
!   
!--                Compute surface-parallel velocity
                   CALL calc_uvw_abs_v_vgrid
!
!--                Compute respective friction velocity on staggered grid
                   CALL calc_us
!
!--                Compute respective surface fluxes for momentum and TKE
                   CALL calc_surface_fluxes
                ENDIF
             ENDDO
          ENDIF
!
!--       For urban-type surfaces
          DO  l = 2, 3
             IF ( surf_usm_v(l)%ns >= 1 )  THEN
                surf => surf_usm_v(l)
!
!--             Compute surface-parallel velocity
                CALL calc_uvw_abs_v_vgrid
!
!--             Compute respective friction velocity on staggered grid
                CALL calc_us
!
!--             Compute respective surface fluxes for momentum and TKE
                CALL calc_surface_fluxes
             ENDIF
          ENDDO
          mom_uv = .FALSE.
!
!--       Calculate horizontal momentum fluxes of w (wsus and wsvs) at vertical 
!--       surfaces.
          mom_w = .TRUE.
!
!--       Default-type surfaces
          DO  l = 0, 3
             IF ( surf_def_v(l)%ns >= 1 )  THEN
                surf => surf_def_v(l)
                CALL calc_uvw_abs_v_wgrid
                CALL calc_us
                CALL calc_surface_fluxes
             ENDIF
          ENDDO 
!
!--       Natural-type surfaces
          DO  l = 0, 3
             IF ( surf_lsm_v(l)%ns >= 1 )  THEN
                surf => surf_lsm_v(l)
                CALL calc_uvw_abs_v_wgrid
                CALL calc_us
                CALL calc_surface_fluxes
             ENDIF
          ENDDO 
!
!--       Urban-type surfaces
          DO  l = 0, 3
             IF ( surf_usm_v(l)%ns >= 1 )  THEN
                surf => surf_usm_v(l)
                CALL calc_uvw_abs_v_wgrid
                CALL calc_us
                CALL calc_surface_fluxes
             ENDIF
          ENDDO 
          mom_w = .FALSE.   
!
!--       Calculate momentum fluxes usvs, vsus, wsus and wsvs at vertical 
!--       surfaces for TKE production. Note, here, momentum fluxes are defined 
!--       at grid center and are not staggered as before.
          mom_tke = .TRUE.
!
!--       Default-type surfaces
          DO  l = 0, 3
             IF ( surf_def_v(l)%ns >= 1 )  THEN
                surf => surf_def_v(l)
                CALL calc_uvw_abs_v_sgrid
                CALL calc_us
                CALL calc_surface_fluxes
             ENDIF
          ENDDO 
!
!--       Natural-type surfaces
          DO  l = 0, 3
             IF ( surf_lsm_v(l)%ns >= 1 )  THEN
                surf => surf_lsm_v(l)
                CALL calc_uvw_abs_v_sgrid
                CALL calc_us
                CALL calc_surface_fluxes
             ENDIF
          ENDDO 
!
!--       Urban-type surfaces
          DO  l = 0, 3
             IF ( surf_usm_v(l)%ns >= 1 )  THEN
                surf => surf_usm_v(l)
                CALL calc_uvw_abs_v_sgrid
                CALL calc_us
                CALL calc_surface_fluxes
             ENDIF
          ENDDO 
          mom_tke = .FALSE.
  
       ENDIF

    END SUBROUTINE surface_layer_fluxes


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializing actions for the surface layer routine. Basically, this involves
!> the preparation of a lookup table for the the bulk Richardson number vs 
!> Obukhov length L when using the lookup table method.
!------------------------------------------------------------------------------!
    SUBROUTINE init_surface_layer_fluxes

       IMPLICIT NONE

       LOGICAL :: terminate_run_l = .FALSE.    !< Flag to terminate run (global)

       IF ( most_method == 'lookup' .AND. TRIM(constant_flux_layer) == 'bottom' )  THEN
!
!--       Check roughness homogeneity between different surface types. 
          IF ( surf_lsm_h%ns >= 1  .AND.  surf_def_h(0)%ns >= 1 )  THEN
             IF ( MINVAL( surf_lsm_h%z0h ) /= MAXVAL( surf_def_h(0)%z0h )  .OR.    &
                  MINVAL( surf_lsm_h%z0  ) /= MAXVAL( surf_def_h(0)%z0  ) )  THEN
                terminate_run_l = .TRUE.
             ENDIF
          ENDIF
          IF ( surf_usm_h%ns >= 1  .AND.  surf_def_h(0)%ns >= 1 )  THEN
             IF ( MINVAL( surf_usm_h%z0h ) /= MAXVAL( surf_def_h(0)%z0h )  .OR.    &
                  MINVAL( surf_usm_h%z0  ) /= MAXVAL( surf_def_h(0)%z0  ) )  THEN
                terminate_run_l = .TRUE.
             ENDIF
          ENDIF
          IF ( surf_usm_h%ns >= 1  .AND.  surf_lsm_h%ns >= 1 )  THEN
             IF ( MINVAL( surf_usm_h%z0h ) /= MAXVAL( surf_lsm_h%z0h )  .OR.       &
                  MINVAL( surf_usm_h%z0  ) /= MAXVAL( surf_lsm_h%z0  ) )  THEN
                terminate_run_l = .TRUE.
             ENDIF
          ENDIF


#if defined( __parallel )
!
!--       Make a logical OR for all processes. Force termination of model if result
!--       is TRUE
          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_ALLREDUCE( terminate_run_l, terminate_run, 1, MPI_LOGICAL,  &
                              MPI_LOR, comm2d, ierr )
#else
          terminate_run = terminate_run_l
#endif

          IF ( terminate_run )  THEN
             message_string = 'most_method = "lookup" cannot be used in ' //   &
                              'combination with a prescribed roughness '  //   &
                              'heterogeneity'
             CALL message( 'surface_layer_fluxes', 'PA0116', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       IF ( surf_def_h(0)%ns >= 1 .AND. TRIM(constant_flux_layer) == 'bottom' )&
          THEN
          surf => surf_def_h(0)
          CALL init_surface_layer_fluxes_m
       ENDIF
       IF ( surf_lsm_h%ns    >= 1 .AND. TRIM(constant_flux_layer) == 'bottom' ) THEN
          surf => surf_lsm_h
          CALL init_surface_layer_fluxes_m
       ENDIF
       IF ( surf_usm_h%ns    >= 1 .AND. TRIM(constant_flux_layer) == 'bottom' ) THEN
          surf => surf_usm_h
          CALL init_surface_layer_fluxes_m
       ENDIF
       IF ( surf_def_h(2)%ns >= 1 .AND. TRIM(constant_flux_layer) == 'top' )   &
          THEN
          surf => surf_def_h(2)
          CALL init_surface_layer_fluxes_m
       ENDIF

    END SUBROUTINE init_surface_layer_fluxes

!-- Subroutine to initialize each surface
    SUBROUTINE init_surface_layer_fluxes_m

       IMPLICIT NONE

       INTEGER(iwp) :: li,         & !< Index for loop to create lookup table
                       num_steps_n   !< Number of non-stretched zeta steps

       LOGICAL :: terminate_run_l = .FALSE.    !< Flag to terminate run (global)

       REAL(wp), PARAMETER ::  zeta_stretch = -10.0_wp !< Start of stretching in the free convection limit
                               
       REAL(wp), DIMENSION(:), ALLOCATABLE :: zeta_tmp

       REAL(wp) :: zeta_step,            & !< Increment of zeta
                   regr      = 1.01_wp,  & !< Stretching factor of zeta_step in the free convection limit
                   regr_old  = 1.0E9_wp, & !< Stretching factor of last iteration step
                   z0h_min   = 0.0_wp,   & !< Minimum value of z0h to create table
                   z0_min    = 0.0_wp      !< Minimum value of z0 to create table


!
!--    In case of runs with neutral statification, set Obukhov length to a
!--    large value
       IF ( neutral .AND. surf%ns >= 1 )  surf%ol = 1.0E10_wp

       IF ( most_method == 'lookup' )  THEN

!
!--       Check for roughness heterogeneity. In that case terminate run and
!--       inform user. Check for both, natural and non-natural walls.
          IF ( surf%ns >= 1 )  THEN
             IF ( MINVAL( surf%z0h ) /= MAXVAL( surf%z0h )  .OR. &
                  MINVAL( surf%z0  ) /= MAXVAL( surf%z0  ) )  THEN
                terminate_run_l = .TRUE.
             ENDIF
          ENDIF

#if defined( __parallel )
!
!--       Make a logical OR for all processes. Force termiation of model if result
!--       is TRUE
          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_ALLREDUCE( terminate_run_l, terminate_run, 1, MPI_LOGICAL,  &
                              MPI_LOR, comm2d, ierr )
#else
          terminate_run = terminate_run_l
#endif

          IF ( terminate_run )  THEN
             message_string = 'most_method = "lookup" cannot be used in ' //   &
                              'combination with a prescribed roughness '  //   &
                              'heterogeneity'
             CALL message( 'surface_layer_fluxes', 'PA0116', 1, 2, 0, 6, 0 )
          ENDIF

          ALLOCATE(  zeta_tmp(0:num_steps-1) )

!
!--       Use the lowest possible value for z_mo
          IF ( TRIM(constant_flux_layer) == 'top'    ) z_mo = zu(nzt+1) - zw(nzt)
          IF ( TRIM(constant_flux_layer) == 'bottom' ) z_mo = zu(nzb+1) - zw(nzb)

!
!--       Calculate z/L range from zeta_stretch to zeta_max using 90% of the
!--       available steps (num_steps). The calculation is done with negative
!--       values of zeta in order to simplify the stretching in the free
!--       convection limit for the remaining 10% of steps.
          zeta_tmp(0) = - zeta_max
          num_steps_n = ( num_steps * 9 / 10 ) - 1
          zeta_step   = (zeta_max - zeta_stretch) / REAL(num_steps_n)

          DO li = 1, num_steps_n
             zeta_tmp(li) = zeta_tmp(li-1) + zeta_step
          ENDDO

!
!--       Calculate stretching factor for the free convection range
          DO  WHILE ( ABS( (regr-regr_old) / regr_old ) > 1.0E-10_wp )
             regr_old = regr
             regr = ( 1.0_wp - ( -zeta_min / zeta_step ) * ( 1.0_wp - regr )   &
                    )**( 10.0_wp / REAL(num_steps) )
          ENDDO

!
!--       Calculate z/L range from zeta_min to zeta_stretch
          DO li = num_steps_n+1, num_steps-1
             zeta_tmp(li) = zeta_tmp(li-1) + zeta_step
             zeta_step = zeta_step * regr
          ENDDO

!
!--       Save roughness lengths to temporary variables
          z0h_min = surf%z0h(1)
          z0_min  = surf%z0(1)
!
!--       Calculate lookup table for the Richardson number versus Obukhov length
!--       The Richardson number (rib) is defined depending on the choice of 
!--       boundary conditions for temperature
          IF ( ibc_pt_b == 1 .OR. ibc_pt_t == 1 )  THEN
             DO li = 0, num_steps-1
                ol_tab(li)  = - z_mo / zeta_tmp(num_steps-1-li)
                rib_tab(li) = z_mo / ol_tab(li)  / ( LOG( z_mo / z0_min )       &
                                                - psi_m( z_mo / ol_tab(li) )    &
                                                + psi_m( z0_min / ol_tab(li) )  &
                                                  )**3
             ENDDO  
          ELSE
             DO li = 0, num_steps-1
                ol_tab(li)  = - z_mo / zeta_tmp(num_steps-1-li)
                rib_tab(li) = z_mo / ol_tab(li)  * ( LOG( z_mo / z0h_min )     &
                                              - psi_h( z_mo / ol_tab(li) )     &
                                              + psi_h( z0h_min / ol_tab(li) )  &
                                            )                                  &
                                          / ( LOG( z_mo / z0_min )             &
                                              - psi_m( z_mo / ol_tab(li) )     &
                                              + psi_m( z0_min / ol_tab(li) )   &
                                            )**2
             ENDDO
          ENDIF

!
!--       Determine minimum values of rib in the lookup table. Set upper limit 
!--       to critical Richardson number (0.25)
          rib_min  = MINVAL(rib_tab)
          rib_max  = 0.25 !MAXVAL(rib_tab)

          DEALLOCATE( zeta_tmp )
       ENDIF

    END SUBROUTINE init_surface_layer_fluxes_m


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the absolute value of the horizontal velocity (relative to the
!> surface) for horizontal surface elements. This is required by all methods.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_uvw_abs

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  ibit          !< k-index away from surface for velocity
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  m             !< running index surface elements

       ibit = MERGE( -1 * k_offset_mcphee, 0, TRIM(most_method) == 'mcphee' )

       DO  m = 1, surf%ns

          i   = surf%i(m)            
          j   = surf%j(m)
          k   = surf%k(m)
!
!--       Compute the absolute value of the horizontal velocity.
!--       (relative to the surface in case the lower surface is the ocean).
!--       Please note, in new surface modelling concept the index values changed,
!--       i.e. the reference grid point is not the surface-grid point itself but
!--       the first grid point outside of the topography. 
!--       Note, in case of coupled ocean-atmosphere simulations relative velocity
!--       with respect to the ocean surface is used, hence, (k-1,j,i) values
!--       are used to calculate the absolute velocity. 
!--       However, this do not apply for downward-facing walls. To mask this, 
!--       use ibit, which checks for upward/downward-facing surfaces.
!--       Note: a previous version of the code computed the average relative 
!--       velocity over two grid cells: i and i+1 for u and j and j+1 for v
          IF ( TRIM(most_method) == 'mcphee' ) THEN
             surf%uvw_abs(m) = SQRT( ( u(k-surf%koff,j,i) - u(k+ibit,j,i) ) **2 +&
                                     ( v(k-surf%koff,j,i) - v(k+ibit,j,i) ) **2  &
                                   ) 
          ELSE
             surf%uvw_abs(m) = SQRT( ( u(k,j,i) - u(k+ibit,j,i) ) **2 +        &
                                     ( v(k,j,i) - v(k+ibit,j,i) ) **2          &
                                   )
          ENDIF
       ENDDO

    END SUBROUTINE calc_uvw_abs


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the absolute value of the horizontal velocity (relative to the
!> surface) for horizontal surface elements. This is required by all methods.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_uvw_abs_v_ugrid

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  m             !< running index surface elements

       REAL(wp)     ::  u_i
       REAL(wp)     ::  w_i


       DO  m = 1, surf%ns
          i   = surf%i(m)            
          j   = surf%j(m)
          k   = surf%k(m)
!
!--       Compute the absolute value of the surface parallel velocity on u-grid.
          u_i  = u(k,j,i)
          w_i  = 0.25_wp * ( w(k-1,j,i-1) + w(k-1,j,i) +                       &
                             w(k,j,i-1)   + w(k,j,i) ) 

          surf%uvw_abs(m) = SQRT( u_i**2 + w_i**2 ) 

       ENDDO

    END SUBROUTINE calc_uvw_abs_v_ugrid

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the absolute value of the horizontal velocity (relative to the
!> surface) for horizontal surface elements. This is required by all methods.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_uvw_abs_v_vgrid

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  m             !< running index surface elements

       REAL(wp)     ::  v_i
       REAL(wp)     ::  w_i


       DO  m = 1, surf%ns
          i   = surf%i(m)            
          j   = surf%j(m)
          k   = surf%k(m)

          v_i  = u(k,j,i)
          w_i  = 0.25_wp * ( w(k-1,j-1,i) + w(k-1,j,i) +                       &
                             w(k,j-1,i)   + w(k,j,i) ) 

          surf%uvw_abs(m) = SQRT( v_i**2 + w_i**2 ) 

       ENDDO

    END SUBROUTINE calc_uvw_abs_v_vgrid

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the absolute value of the horizontal velocity (relative to the
!> surface) for horizontal surface elements. This is required by all methods.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_uvw_abs_v_wgrid

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  m             !< running index surface elements

       REAL(wp)     ::  u_i
       REAL(wp)     ::  v_i
       REAL(wp)     ::  w_i
!
!--    North- (l=0) and south-facing (l=1) surfaces 
       IF ( l == 0  .OR.  l == 1 )  THEN
          DO  m = 1, surf%ns
             i   = surf%i(m)            
             j   = surf%j(m)
             k   = surf%k(m)

             u_i  = 0.25_wp * ( u(k+1,j,i+1) + u(k+1,j,i) +                    &
                                u(k,j,i+1)   + u(k,j,i) )
             v_i  = 0.0_wp
             w_i  = w(k,j,i)

             surf%uvw_abs(m) = SQRT( u_i**2 + v_i**2 + w_i**2 ) 
          ENDDO
!
!--    East- (l=2) and west-facing (l=3) surfaces
       ELSE
          DO  m = 1, surf%ns
             i   = surf%i(m)            
             j   = surf%j(m)
             k   = surf%k(m)

             u_i  = 0.0_wp
             v_i  = 0.25_wp * ( v(k+1,j+1,i) + v(k+1,j,i) +                    &
                                v(k,j+1,i)   + v(k,j,i) )
             w_i  = w(k,j,i)

             surf%uvw_abs(m) = SQRT( u_i**2 + v_i**2 + w_i**2 ) 
          ENDDO
       ENDIF            

    END SUBROUTINE calc_uvw_abs_v_wgrid

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the absolute value of the horizontal velocity (relative to the
!> surface) for horizontal surface elements. This is required by all methods.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_uvw_abs_v_sgrid

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  m             !< running index surface elements

       REAL(wp)     ::  u_i
       REAL(wp)     ::  v_i
       REAL(wp)     ::  w_i

!
!--    North- (l=0) and south-facing (l=1) walls 
       IF ( l == 0  .OR.  l == 1 )  THEN
          DO  m = 1, surf%ns
             i   = surf%i(m)            
             j   = surf%j(m)
             k   = surf%k(m)

             u_i = 0.5_wp * ( u(k,j,i) + u(k,j,i+1) )
             v_i = 0.0_wp
             w_i = 0.5_wp * ( w(k,j,i) + w(k-1,j,i) )

             surf%uvw_abs(m) = SQRT( u_i**2 + v_i**2 + w_i**2 ) 
          ENDDO
!
!--    East- (l=2) and west-facing (l=3) walls 
       ELSE
          DO  m = 1, surf%ns
             i   = surf%i(m)            
             j   = surf%j(m)
             k   = surf%k(m)

             u_i = 0.0_wp
             v_i = 0.5_wp * ( v(k,j,i) + v(k,j+1,i) )
             w_i = 0.5_wp * ( w(k,j,i) + w(k-1,j,i) )

             surf%uvw_abs(m) = SQRT( u_i**2 + v_i**2 + w_i**2 ) 
          ENDDO
       ENDIF  

    END SUBROUTINE calc_uvw_abs_v_sgrid


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the Obukhov length (L) and Richardson flux number (z/L)
!------------------------------------------------------------------------------!
    SUBROUTINE calc_ol

       IMPLICIT NONE

       INTEGER(iwp) ::  ibit          !< flag to mask computation of relative velocity in case of downward-facing surfaces
       INTEGER(iwp) ::  iter    !< Newton iteration step
       INTEGER(iwp) ::  li      !< look index
       INTEGER(iwp) ::  m       !< loop variable over all horizontal wall elements 

       REAL(wp)     :: f,      & !< Function for Newton iteration: f = Ri - [...]/[...]^2 = 0
                       f_d_ol, & !< Derivative of f
                       ol_l,   & !< Lower bound of L for Newton iteration
                       ol_m,   & !< Previous value of L for Newton iteration
                       ol_old, & !< Previous time step value of L
                       ol_u      !< Upper bound of L for Newton iteration

!
!--    Reference velocity is at k+1 for downward-facing surfaces, 
!--    k-1 for upward-facing surfaces.
       ibit = MERGE( 1, -1, downward )

       IF ( TRIM( most_method ) /= 'mcphee' .AND. TRIM( most_method ) /= 'circular' )  THEN
!
!--       Evaluate bulk Richardson number (calculation depends on
!--       definition based on setting of boundary conditions
          IF ( ibc_pt_b /= 1 .OR. ibc_pt_t /= 1 )  THEN
             IF ( humidity )  THEN
                !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
                DO  m = 1, surf%ns

                   i   = surf%i(m)            
                   j   = surf%j(m)
                   k   = surf%k(m)

                   z_mo = surf%z_mo(m)

                   surf%rib(m) = g * z_mo *                                    &
                                         ( vpt(k,j,i) - vpt(k+ibit,j,i) ) /    &
                      ( surf%uvw_abs(m)**2 * vpt(k,j,i) + 1.0E-20_wp )
                ENDDO
             ELSE
                !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
                DO  m = 1, surf%ns

                   i   = surf%i(m)            
                   j   = surf%j(m)
                   k   = surf%k(m)

                   z_mo = surf%z_mo(m)

                   surf%rib(m) = g * z_mo *                                    &
                                         ( pt(k,j,i) - pt(k+ibit,j,i)   ) /    &
                      ( surf%uvw_abs(m)**2 * pt(k,j,i)  + 1.0E-20_wp )
                ENDDO
             ENDIF
          ELSE ! not circular
             IF ( humidity )  THEN
                !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
                DO  m = 1, surf%ns

                   i   = surf%i(m)            
                   j   = surf%j(m)
                   k   = surf%k(m)

                   z_mo = surf%z_mo(m)

                   surf%rib(m) = - g * z_mo * ( ( 1.0_wp + 0.61_wp             &
                           * q(k,j,i) ) * surf%shf(m) + 0.61_wp                &
                           * pt(k,j,i) * surf%qsws(m) ) *                      &
                             drho_ref_zw(k+ibit)                /              &
                         ( surf%uvw_abs(m)**3 * vpt(k,j,i) * kappa**2          &
                           + 1.0E-20_wp )
                ENDDO
             ELSE ! not humidity
                !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
                DO  m = 1, surf%ns

                   i   = surf%i(m)            
                   j   = surf%j(m)
                   k   = surf%k(m)

                   z_mo = surf%z_mo(m)

                   surf%rib(m) = - g * z_mo * surf%shf(m) *                    &
                                        drho_ref_zw(k+ibit)            /       &
                        ( surf%uvw_abs(m)**3 * pt(k,j,i) * kappa**2            &
                           + 1.0E-20_wp )
                ENDDO
             ENDIF
          ENDIF

       ENDIF


!
!--    Calculate the Obukhov length either using a Newton iteration
!--    method, via a lookup table, or using the old circular way
       IF ( TRIM( most_method ) == 'newton' )  THEN

          DO  m = 1, surf%ns

             i   = surf%i(m)            
             j   = surf%j(m)

             z_mo = surf%z_mo(m)

!
!--          Store current value in case the Newton iteration fails
             ol_old = surf%ol(m)

!
!--          Ensure that the bulk Richardson number and the Obukhov 
!--          length have the same sign
             IF ( surf%rib(m) * surf%ol(m) < 0.0_wp  .OR.                      &
                  ABS( surf%ol(m) ) == ol_max )  THEN
                IF ( surf%rib(m) > 1.0_wp ) surf%ol(m) =  0.01_wp
                IF ( surf%rib(m) < 0.0_wp ) surf%ol(m) = -0.01_wp
             ENDIF
!
!--          Iteration to find Obukhov length
             iter = 0
             DO
                iter = iter + 1
!
!--             In case of divergence, use the value of the previous time step
                IF ( iter > 1000 )  THEN
                   surf%ol(m) = ol_old
                   EXIT
                ENDIF

                ol_m = surf%ol(m)
                ol_l = ol_m - 0.001_wp * ol_m
                ol_u = ol_m + 0.001_wp * ol_m


                IF ( ibc_pt_b /= 1 .OR. ibc_pt_t /= 1 )  THEN
!
!--                Calculate f = Ri - [...]/[...]^2 = 0
                   f = surf%rib(m) - ( z_mo / ol_m ) * (                       &
                                                 LOG( z_mo / surf%z0h(m) )     &
                                                 - psi_h( z_mo / ol_m )        &
                                                 + psi_h( surf%z0h(m) /        &
                                                          ol_m )               &
                                                               )               &
                                              / ( LOG( z_mo / surf%z0(m) )     &
                                                 - psi_m( z_mo / ol_m )        &
                                                 + psi_m( surf%z0(m) /         &
                                                          ol_m )               &
                                                 )**2

!
!--                 Calculate df/dL
                    f_d_ol = ( - ( z_mo / ol_u ) * ( LOG( z_mo /               &
                                                             surf%z0h(m) )     &
                                            - psi_h( z_mo / ol_u )             &
                                            + psi_h( surf%z0h(m) / ol_u )      &
                                              )                                &
                                            / ( LOG( z_mo / surf%z0(m) )       &
                                            - psi_m( z_mo / ol_u )             &
                                            + psi_m( surf%z0(m) / ol_u )       &
                                              )**2                             &
                           + ( z_mo / ol_l ) * ( LOG( z_mo / surf%z0h(m) )     &
                                            - psi_h( z_mo / ol_l )             &
                                            + psi_h( surf%z0h(m) / ol_l )      &
                                              )                                &
                                            / ( LOG( z_mo / surf%z0(m) )       &
                                            - psi_m( z_mo / ol_l )             &
                                            + psi_m( surf%z0(m) / ol_l )       &
                                              )**2                             &
                             ) / ( ol_u - ol_l )
                ELSE
!
!--                Calculate f = Ri - 1 /[...]^3 = 0
                   f = surf%rib(m) - ( z_mo / ol_m ) /                         &
                                                ( LOG( z_mo / surf%z0(m) )     &
                                            - psi_m( z_mo / ol_m )             &
                                            + psi_m( surf%z0(m) / ol_m )       &
                                                )**3

!
!--                Calculate df/dL
                   f_d_ol = ( - ( z_mo / ol_u ) / ( LOG( z_mo /                &
                                                            surf%z0(m) )       &
                                            - psi_m( z_mo / ol_u )             &
                                            + psi_m( surf%z0(m) / ol_u )       &
                                                     )**3                      &
                           + ( z_mo / ol_l ) / ( LOG( z_mo / surf%z0(m) )      &
                                            - psi_m( z_mo / ol_l )             &
                                            + psi_m( surf%z0(m) / ol_l )       &
                                               )**3                            &
                                  ) / ( ol_u - ol_l )
                ENDIF
!
!--             Calculate new L
                surf%ol(m) = ol_m - f / f_d_ol

!
!--             Ensure that the bulk Richardson number and the Obukhov 
!--             length have the same sign and ensure convergence.
                IF ( surf%ol(m) * ol_m < 0.0_wp )  surf%ol(m) = ol_m * 0.5_wp
!
!--             If unrealistic value occurs, set L to the maximum
!--             value that is allowed
                IF ( ABS( surf%ol(m) ) > ol_max )  THEN
                   surf%ol(m) = ol_max
                   EXIT
                ENDIF
!
!--             Check for convergence
                IF ( ABS( ( surf%ol(m) - ol_m ) /                              &
                     surf%ol(m) ) < 1.0E-4_wp )  THEN
                   EXIT
                ELSE
                   CYCLE
                ENDIF

             ENDDO
          ENDDO

       ELSEIF ( TRIM( most_method ) == 'lookup' )  THEN

          !$OMP PARALLEL DO PRIVATE( i, j, z_mo, li ) FIRSTPRIVATE( li_bnd ) LASTPRIVATE( li_bnd )
          DO  m = 1, surf%ns

             i   = surf%i(m)            
             j   = surf%j(m)
!
!--          If the bulk Richardson number is outside the range of the lookup
!--          table, set it to the exceeding threshold value
             IF ( surf%rib(m) < rib_min )  surf%rib(m) = rib_min
             IF ( surf%rib(m) > rib_max )  surf%rib(m) = rib_max

!
!--          Find the correct index bounds for linear interpolation. As the
!--          Richardson number will not differ very much from time step to
!--          time step , use the index from the last step and search in the 
!--          correct direction
             li = li_bnd
             IF ( rib_tab(li) - surf%rib(m) > 0.0_wp )  THEN
                DO  WHILE ( rib_tab(li-1) - surf%rib(m) > 0.0_wp  .AND.  li > 0 )
                   li = li-1
                ENDDO
             ELSE
                DO  WHILE ( rib_tab(li) - surf%rib(m) < 0.0_wp                 &
                           .AND.  li < num_steps-1 )
                   li = li+1
                ENDDO
             ENDIF
             li_bnd = li

!
!--          Linear interpolation to find the correct value of z/L
             surf%ol(m) = ( ol_tab(li-1) + ( ol_tab(li) - ol_tab(li-1) )       &
                         / (  rib_tab(li) - rib_tab(li-1) )                    &
                         * ( surf%rib(m) - rib_tab(li-1) ) )

          ENDDO

       ELSEIF ( TRIM( most_method ) == 'circular' )  THEN

          IF ( .NOT. humidity )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
             DO  m = 1, surf%ns

                i   = surf%i(m)            
                j   = surf%j(m)
                k   = surf%k(m)

                z_mo = surf%z_mo(m)

                surf%ol(m) =  ( pt(k,j,i) *  surf%us(m)**2 ) /                 &
                                       ( kappa * g *                           &
                                         surf%ts(m) + 1E-30_wp )
!
!--             Limit the value range of the Obukhov length.
!--             This is necessary for very small velocities (u,v --> 1), because
!--             the absolute value of ol can then become very small, which in
!--             consequence would result in very large shear stresses and very
!--             small momentum fluxes (both are generally unrealistic).
                IF ( ( z_mo / ( surf%ol(m) + 1E-30_wp ) ) < zeta_min )         &
                   surf%ol(m) = z_mo / zeta_min
                IF ( ( z_mo / ( surf%ol(m) + 1E-30_wp ) ) > zeta_max )         &
                   surf%ol(m) = z_mo / zeta_max

             ENDDO
          ELSEIF ( cloud_physics  .OR.  cloud_droplets )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
             DO  m = 1, surf%ns

                i   = surf%i(m)            
                j   = surf%j(m)
                k   = surf%k(m)

                z_mo = surf%z_mo(m)


                surf%ol(m) =  ( vpt(k,j,i) * surf%us(m)**2 ) /                 &
                    ( kappa * g * ( surf%ts(m) +                               &
                         0.61_wp * surf%pt1(m) * surf%us(m)                    &
                       + 0.61_wp * surf%qv1(m) * surf%ts(m) -                  &
                                   surf%ts(m)  * ql(k,j,i) ) + 1E-30_wp )
!
!--             Limit the value range of the Obukhov length.
!--             This is necessary for very small velocities (u,v --> 1), because
!--             the absolute value of ol can then become very small, which in
!--             consequence would result in very large shear stresses and very
!--             small momentum fluxes (both are generally unrealistic).
                IF ( ( z_mo / ( surf%ol(m) + 1E-30_wp ) ) < zeta_min )         &
                   surf%ol(m) = z_mo / zeta_min
                IF ( ( z_mo / ( surf%ol(m) + 1E-30_wp ) ) > zeta_max )         &
                   surf%ol(m) = z_mo / zeta_max

             ENDDO
          ELSE

             !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
             DO  m = 1, surf%ns

                i   = surf%i(m)            
                j   = surf%j(m)
                k   = surf%k(m)

                z_mo = surf%z_mo(m)

                surf%ol(m) =  ( vpt(k,j,i) *  surf%us(m)**2 ) /                &
                    ( kappa * g * ( surf%ts(m) + 0.61_wp * pt(k,j,i) *         &
                                    surf%qs(m) + 0.61_wp * q(k,j,i)  *         &
                                    surf%ts(m) ) + 1E-30_wp )

!
!--             Limit the value range of the Obukhov length.
!--             This is necessary for very small velocities (u,v --> 1), because
!--             the absolute value of ol can then become very small, which in
!--             consequence would result in very large shear stresses and very
!--             small momentum fluxes (both are generally unrealistic).
                IF ( ( z_mo / ( surf%ol(m) + 1E-30_wp ) ) < zeta_min )         &
                   surf%ol(m) = z_mo / zeta_min
                IF ( ( z_mo / ( surf%ol(m) + 1E-30_wp ) ) > zeta_max )         &
                   surf%ol(m) = z_mo / zeta_max

             ENDDO

          ENDIF

       ELSEIF ( trim(most_method) == 'mcphee' ) THEN

          !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
          DO  m = 1, surf%ns

             i   = surf%i(m)            
             j   = surf%j(m)
             k   = surf%k(m)

             surf%ol(m) = rho_ocean(k,j,i) * surf%us(m)**3 /                   &
                          (g * kappa *                                         &
                            ( beta_S(k,j,i)*surf%sasws(m) -                    &
                              alpha_T(k,j,i)*surf%shf(m)    ) + 1E-30_wp )

          ENDDO

       ENDIF

    END SUBROUTINE calc_ol

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the Obukhov length for an individual horizontal surface element. 
!> This is required by mcphee method only.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_ol_m(m)

       IMPLICIT NONE

       INTEGER(iwp) ::  i,j,k   !< index variables 
       INTEGER(iwp) ::  m       !< loop variable over all horizontal surf elements 

       i   = surf%i(m)            
       j   = surf%j(m)
       k   = surf%k(m)

       surf%ol(m) = rho_ocean(k,j,i) * surf%us(m)**3 /                   &
                    (g * kappa *                                         &
                      ( beta_S(k,j,i)*surf%sasws(m) -                    &
                        alpha_T(k,j,i)*surf%shf(m)    ) + 1E-30_wp )

    END SUBROUTINE calc_ol_m

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the friction velocity u*. This is required by all methods.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_us

       IMPLICIT NONE

       INTEGER(iwp) ::  m       !< loop variable over all horizontal surf elements 

!
!--    Compute u* at horizontal surfaces at the scalars' grid points
       IF ( .NOT. surf_vertical )  THEN
!
!--       Compute u* at upward-facing surfaces
          IF ( .NOT. downward )  THEN
             !$OMP PARALLEL  DO PRIVATE( z_mo )
             DO  m = 1, surf%ns

                z_mo = surf%z_mo(m)
!
!--             Compute u* at the scalars' grid points
                surf%us(m) = kappa * surf%uvw_abs(m) /                         &
                            ( LOG( z_mo / surf%z0(m) )                         &
                           - psi_m( z_mo / surf%ol(m) )                        &
                           + psi_m( surf%z0(m) / surf%ol(m) ) )
   
             ENDDO
!
!--       Compute u* at downward-facing surfaces. This case, do not consider
!--       any stability. 
          ELSE
             IF ( trim(most_method) == 'mcphee' ) THEN
                DO  m = 1, surf%ns
                   ! Equivalent to: surf%us(m) = SQRT(drag_coeff) * surf%uvw_abs(m)
                   surf%us(m) = kappa * surf%uvw_abs(m) / LOG( z_offset_mcphee / surf%z0(m) ) 
                                
                ENDDO
             ELSE
                !$OMP PARALLEL  DO PRIVATE( z_mo )
                DO  m = 1, surf%ns

                   z_mo = surf%z_mo(m)
!
!--                Compute u* at the scalars' grid points
                   surf%us(m) = kappa * surf%uvw_abs(m) / LOG( z_mo / surf%z0(m) )
             
                ENDDO
             ENDIF
          ENDIF
!
!--    Compute u* at vertical surfaces at the u/v/v grid, respectively. 
!--    No stability is considered in this case.
       ELSE
          !$OMP PARALLEL DO PRIVATE( z_mo )
          DO  m = 1, surf%ns
             z_mo = surf%z_mo(m)

             surf%us(m) = kappa * surf%uvw_abs(m) / LOG( z_mo / surf%z0(m) )

          ENDDO
       ENDIF

    END SUBROUTINE calc_us

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate potential temperature and specific humidity at first grid level
!> ( only for upward-facing surfs ).
!------------------------------------------------------------------------------!
    SUBROUTINE calc_pt_q

       IMPLICIT NONE

       INTEGER(iwp) ::  m       !< loop variable over all horizontal surf elements 

       !$OMP PARALLEL DO PRIVATE( i, j, k )
       DO  m = 1, surf%ns 

          i   = surf%i(m)            
          j   = surf%j(m)
          k   = surf%k(m)

          IF ( cloud_physics ) THEN
             surf%pt1(m) = pt(k,j,i) + l_d_cp * pt_d_t(k) * ql(k,j,i)
             surf%qv1(m) = q(k,j,i) - ql(k,j,i)
          ELSEIF( cloud_droplets ) THEN
             surf%pt1(m) = pt(k,j,i) + l_d_cp * pt_d_t(k) * ql(k,j,i)
             surf%qv1(m) = q(k,j,i) 
          ELSE
             surf%pt1(m) = pt(k,j,i)
             IF ( humidity )  THEN
                surf%qv1(m) = q(k,j,i)
             ELSE
                surf%qv1(m) = 0.0_wp
             ENDIF
          ENDIF 

       ENDDO

    END SUBROUTINE calc_pt_q

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate potential temperature and salinity at nth level below the top surface
!------------------------------------------------------------------------------!
    SUBROUTINE calc_pt_sa

       IMPLICIT NONE

       INTEGER(iwp) ::  m       !< loop variable over all horizontal surf elements 
       INTEGER(iwp) ::  ibit    !< k offset for mcphee case 

       ibit = MERGE( -1 * k_offset_mcphee, 0, TRIM(most_method) == 'mcphee' )
       
       !$OMP PARALLEL DO PRIVATE( i, j, k )
       DO  m = 1, surf%ns

          i   = surf%i(m)            
          j   = surf%j(m)
          k   = surf%k(m)

          surf%pt1(m) = pt(k+ibit,j,i)
          surf%sa1(m) = sa(k+ibit,j,i)

       ENDDO

    END SUBROUTINE calc_pt_sa

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate potential temperature and specific humidity at first grid level
!> ( only for upward-facing surfs )
!------------------------------------------------------------------------------!
    SUBROUTINE calc_pt_surface

       IMPLICIT NONE

       INTEGER(iwp) ::  koff    !< index offset between surface and atmosphere grid point (-1 for upward-, +1 for downward-facing walls) 
       INTEGER(iwp) ::  m       !< loop variable over all horizontal surf elements 
       
       koff = surf%koff
       !$OMP PARALLEL DO PRIVATE( i, j, k )
       DO  m = 1, surf%ns 

          i   = surf%i(m)            
          j   = surf%j(m)
          k   = surf%k(m)

          surf%pt_surface(m) = pt(k+koff,j,i)

       ENDDO

    END SUBROUTINE calc_pt_surface

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the other MOST scaling parameters theta*, q*, (qc*, qr*, nc*, nr*)
!------------------------------------------------------------------------------!
    SUBROUTINE calc_scaling_parameters

       IMPLICIT NONE


       INTEGER(iwp) ::  ibit          !< flag to mask computation of relative velocity in case of downward-facing surfaces
       INTEGER(iwp)  ::  m       !< loop variable over all horizontal surf elements 
       INTEGER(iwp)  ::  lsp     !< running index for chemical species

!
!--    Reference velocity is at k+1 for downward-facing surfaces, 
!--    k-1 for upward-facing surfaces.
       ibit = MERGE( 1, -1, downward )

! 
!--    Compute theta* at horizontal surfaces
       IF ( constant_heatflux  .AND.  .NOT. surf_vertical )  THEN
!
!--       For a given heat flux in the surface layer:

          !$OMP PARALLEL DO PRIVATE( i, j, k )
          DO  m = 1, surf%ns 

             i   = surf%i(m)            
             j   = surf%j(m)
             k   = surf%k(m)

             surf%ts(m) = -surf%shf(m) * drho_ref_zw(k+ibit) /                 &
                                  ( surf%us(m) + 1E-30_wp ) 

!
!--          ts must be limited, because otherwise overflow may occur in case
!--          of us=0 when computing ol further below
             IF ( surf%ts(m) < -1.05E5_wp )  surf%ts(m) = -1.0E5_wp
             IF ( surf%ts(m) >  1.0E5_wp  )  surf%ts(m) =  1.0E5_wp

          ENDDO

       ELSEIF ( .NOT. surf_vertical ) THEN
!
!--       For a given surface temperature:
          IF ( large_scale_forcing  .AND.  lsf_surf )  THEN

             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, surf%ns 
                i   = surf%i(m)            
                j   = surf%j(m)
                k   = surf%k(m)

                pt(k+ibit,j,i) = pt_surface
             ENDDO
          ENDIF

          !$OMP PARALLEL DO PRIVATE( z_mo )
          DO  m = 1, surf%ns   

             z_mo = surf%z_mo(m)

             surf%ts(m) = kappa * ( surf%pt1(m) - surf%pt_surface(m) )      &
                                  / ( LOG( z_mo / surf%z0h(m) )             &
                                      - psi_h( z_mo / surf%ol(m) )          &
                                      + psi_h( surf%z0h(m) / surf%ol(m) ) )

          ENDDO

       ENDIF
! 
!--    Compute theta* at vertical surfaces. This is only required in case of 
!--    land-surface model, in order to compute aerodynamical resistance.
       IF ( surf_vertical )  THEN
          !$OMP PARALLEL DO PRIVATE( i, j )
          DO  m = 1, surf%ns 

             i          =  surf%i(m)            
             j          =  surf%j(m)
             surf%ts(m) = -surf%shf(m) / ( surf%us(m) + 1E-30_wp )
!
!--          ts must be limited, because otherwise overflow may occur in case
!--          of us=0 when computing ol further below
             IF ( surf%ts(m) < -1.05E5_wp )  surf%ts(m) = -1.0E5_wp
             IF ( surf%ts(m) >  1.0E5_wp  )  surf%ts(m) =  1.0E5_wp

          ENDDO
       ENDIF

!
!--    If required compute q* at horizontal surfaces
       IF ( humidity )  THEN
          IF ( constant_waterflux  .AND.  .NOT. surf_vertical )  THEN
!
!--          For a given water flux in the surface layer
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, surf%ns  

                i   = surf%i(m)            
                j   = surf%j(m)
                k   = surf%k(m)
                surf%qs(m) = -surf%qsws(m) * drho_ref_zw(k+ibit) /             &
                                               ( surf%us(m) + 1E-30_wp )

             ENDDO

          ELSEIF ( .NOT. surf_vertical )  THEN 
             coupled_run = ( coupling_mode == 'atmosphere_to_ocean'  .AND.     &
                             run_coupled )

             IF ( large_scale_forcing  .AND.  lsf_surf )  THEN
                !$OMP PARALLEL DO PRIVATE( i, j, k )
                DO  m = 1, surf%ns  

                   i   = surf%i(m)           
                   j   = surf%j(m)
                   k   = surf%k(m)
                   q(k+ibit,j,i) = q_surface
                   
                ENDDO
             ENDIF

!
!--          Assume saturation for atmosphere coupled to ocean (but not
!--          in case of precursor runs)
             IF ( coupled_run )  THEN
                !$OMP PARALLEL DO PRIVATE( i, j, k, e_s )
                DO  m = 1, surf%ns   
                   i   = surf%i(m)            
                   j   = surf%j(m)
                   k   = surf%k(m)
                   e_s = 6.1_wp * &
                              EXP( 0.07_wp * ( MIN(pt(k+ibit,j,i),pt(k,j,i))   &
                                               - 273.15_wp ) )
                   q(k+ibit,j,i) = 0.622_wp * e_s / ( surface_pressure - e_s )
                ENDDO
             ENDIF

             IF ( cloud_physics  .OR.  cloud_droplets )  THEN
               !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
                DO  m = 1, surf%ns   

                   i   = surf%i(m)            
                   j   = surf%j(m)
                   k   = surf%k(m)
   
                   z_mo = surf%z_mo(m)

                   surf%qs(m) = kappa * ( surf%qv1(m) - q(k+ibit,j,i) )        &
                                        / ( LOG( z_mo / surf%z0q(m) )          &
                                            - psi_h( z_mo / surf%ol(m) )       &
                                            + psi_h( surf%z0q(m) /             &
                                                     surf%ol(m) ) )
                ENDDO
             ELSE
                !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
                DO  m = 1, surf%ns   

                   i   = surf%i(m)            
                   j   = surf%j(m)
                   k   = surf%k(m)
   
                   z_mo = surf%z_mo(m)

                   surf%qs(m) = kappa * ( q(k,j,i) - q(k-1,j,i) )              &
                                        / ( LOG( z_mo / surf%z0q(m) )          &
                                            - psi_h( z_mo / surf%ol(m) )       &
                                            + psi_h( surf%z0q(m) /             &
                                                     surf%ol(m) ) )
                ENDDO
             ENDIF
          ENDIF
! 
!--       Compute q* at vertical surfaces
          IF ( surf_vertical )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j )
             DO  m = 1, surf%ns  

                i          =  surf%i(m)            
                j          =  surf%j(m)
                surf%qs(m) = -surf%qsws(m) / ( surf%us(m) + 1E-30_wp )

             ENDDO
          ENDIF
       ENDIF
       
!
!--    If required compute s*
       IF ( passive_scalar )  THEN
!
!--       At horizontal surfaces
          IF ( constant_scalarflux  .AND.  .NOT. surf_vertical )  THEN
!
!--          For a given scalar flux in the surface layer
             !$OMP PARALLEL DO PRIVATE( i, j )
             DO  m = 1, surf%ns  
                i   = surf%i(m)            
                j   = surf%j(m)
                surf%ss(m) = -surf%ssws(m) / ( surf%us(m) + 1E-30_wp )
             ENDDO
          ENDIF
!
!--       At vertical surfaces
          IF ( surf_vertical )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j )
             DO  m = 1, surf%ns  
                i          =  surf%i(m)            
                j          =  surf%j(m)
                surf%ss(m) = -surf%ssws(m) / ( surf%us(m) + 1E-30_wp )
             ENDDO
          ENDIF
       ENDIF

!
!--    If required compute cs* (chemical species)
       IF ( air_chemistry  )  THEN  
!
!--       At horizontal surfaces                             
          DO  lsp = 1, nvar
             IF ( constant_csflux(lsp)  .AND.  .NOT.  surf_vertical )  THEN
!--             For a given chemical species' flux in the surface layer
                !$OMP PARALLEL DO PRIVATE( i, j )
                DO  m = 1, surf%ns  
                   i   = surf%i(m)            
                   j   = surf%j(m)
                   surf%css(lsp,m) = -surf%cssws(lsp,m) / ( surf%us(m) + 1E-30_wp )
                ENDDO
             ENDIF
          ENDDO
!
!--       At vertical surfaces
          IF ( surf_vertical )  THEN
             DO  lsp = 1, nvar
                !$OMP PARALLEL DO PRIVATE( i, j )
                DO  m = 1, surf%ns  
                   i   = surf%i(m)            
                   j   = surf%j(m)
                   surf%css(lsp,m) = -surf%cssws(lsp,m) / ( surf%us(m) + 1E-30_wp )
                ENDDO
             ENDDO
          ENDIF
       ENDIF

!
!--    If required compute qc* and nc*
       IF ( cloud_physics  .AND.  microphysics_morrison  .AND.                 &
            .NOT. surf_vertical )  THEN
          !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
          DO  m = 1, surf%ns    

             i    = surf%i(m)            
             j    = surf%j(m)
             k    = surf%k(m)

             z_mo = surf%z_mo(m)

             surf%qcs(m) = kappa * ( qc(k,j,i) - qc(k-1,j,i) )                 &
                                 / ( LOG( z_mo / surf%z0q(m) )                 &
                                 - psi_h( z_mo / surf%ol(m) )                  &
                                 + psi_h( surf%z0q(m) / surf%ol(m) ) )

             surf%ncs(m) = kappa * ( nc(k,j,i) - nc(k-1,j,i) )                 &
                                 / ( LOG( z_mo / surf%z0q(m) )                 &
                                 - psi_h( z_mo / surf%ol(m) )                  &
                                 + psi_h( surf%z0q(m) / surf%ol(m) ) )
          ENDDO

       ENDIF

!
!--    If required compute qr* and nr*
       IF ( cloud_physics  .AND.  microphysics_seifert  .AND.                  &
            .NOT. surf_vertical )  THEN
          !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
          DO  m = 1, surf%ns    

             i    = surf%i(m)            
             j    = surf%j(m)
             k    = surf%k(m)

             z_mo = surf%z_mo(m)

             surf%qrs(m) = kappa * ( qr(k,j,i) - qr(k-1,j,i) )                 &
                                 / ( LOG( z_mo / surf%z0q(m) )                 &
                                 - psi_h( z_mo / surf%ol(m) )                  &
                                 + psi_h( surf%z0q(m) / surf%ol(m) ) )

             surf%nrs(m) = kappa * ( nr(k,j,i) - nr(k-1,j,i) )                 &
                                 / ( LOG( z_mo / surf%z0q(m) )                 &
                                 - psi_h( z_mo / surf%ol(m) )                  &
                                 + psi_h( surf%z0q(m) / surf%ol(m) ) )
          ENDDO

       ENDIF

    END SUBROUTINE calc_scaling_parameters



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate surface fluxes usws, vsws, shf, qsws, (sasws, qcsws, qrsws, ncsws, 
!> nrsws).
!------------------------------------------------------------------------------!
    SUBROUTINE calc_surface_fluxes

       IMPLICIT NONE

       COMPLEX       ::  im = (0,-1), tau

       INTEGER(iwp)  ::  ibit    !< flag to mask computation of relative velocity in case of downward-facing surfaces
       INTEGER(iwp)  ::  m       !< loop variable over all horizontal surf elements
       INTEGER(iwp)  ::  lsp     !< running index for chemical species

       REAL(wp)                            ::  dum         !< dummy to precalculate logarithm
       REAL(wp)                            ::  flag_u      !< flag indicating u-grid, used for calculation of horizontal momentum fluxes at vertical surfaces
       REAL(wp)                            ::  flag_v      !< flag indicating v-grid, used for calculation of horizontal momentum fluxes at vertical surfaces
       REAL(wp)                            ::  melt_min = 0.0_wp,melt_max = 1e3/3e7
       REAL(wp)                            ::  s_factor    !< term in mcphee surface flux calculation
       REAL(wp)                            ::  us_x,us_y,us_mag
       REAL(wp)                            ::  zeta
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  u_i         !< u-component interpolated onto scalar grid point, required for momentum fluxes at vertical surfaces 
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  v_i         !< v-component interpolated onto scalar grid point, required for momentum fluxes at vertical surfaces 
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  w_i         !< w-component interpolated onto scalar grid point, required for momentum fluxes at vertical surfaces 

!
!--    Reference velocity is at k+1 for downward-facing surfaces, 
!--    k-1 for upward-facing surfaces.
       ibit = MERGE( 1, -1, downward )

!
!--    Calcuate surface fluxes at horizontal walls
       IF ( .NOT. surf_vertical )  THEN
!
!--       Compute u'w' for the total model domain at upward-facing surfaces.
!--       First compute the corresponding component of u* and square it.
          IF ( .NOT. downward )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
             DO  m = 1, surf%ns  
   
                i = surf%i(m)            
                j = surf%j(m)
                k = surf%k(m)


                surf%usws(m) = kappa * ( u(k,j,i) - u(k-1,j,i) )               &
                              / ( LOG( surf%z_mo(m) / surf%z0(m) )             &
                                  - psi_m( surf%z_mo(m) / surf%ol(m) )         &
                                  + psi_m( surf%z0(m) / surf%ol(m) ) )
                surf%vsws(m) = kappa * ( v(k,j,i) - v(k-1,j,i) )               &
                           / ( LOG( surf%z_mo(m) / surf%z0(m) )                &
                               - psi_m( surf%z_mo(m) / surf%ol(m) )            &
                               + psi_m( surf%z0(m) / surf%ol(m) ) )
!
!--             Please note, the computation of usws is not fully accurate. Actually 
!--             a further interpolation of us onto the u-grid, where usws is defined, 
!--             is required. However, this is not done as this would require several
!--             data transfers between 2D-grid and the surf-type. 
!--             The impact of the missing interpolation is negligible as several 
!--             tests had shown. 
!--             Same also for ol.  
                surf%usws(m) = -surf%usws(m) * surf%us(m) * rho_ref_zw(k-1)
                surf%vsws(m) = -surf%vsws(m) * surf%us(m) * rho_ref_zw(k-1)
             
             ENDDO
!
!--       At downward-facing surfaces
          ELSE

             !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
             DO  m = 1, surf%ns  
   
                i = surf%i(m)            
                j = surf%j(m)
                k = surf%k(m)

                IF ( trim(most_method) == 'mcphee' ) THEN

                   us_x = kappa * ( u(k,j,i) - u(k+1,j,i) )                    &
                          / LOG( z_offset_mcphee / surf%z0(m) )
                   us_y = kappa * ( v(k,j,i) - v(k+1,j,i) )                    &
                          / LOG( z_offset_mcphee / surf%z0(m) )
                   eta_star = ( 1.0_wp + ( xi_N * surf%us(m) ) /               &
                              ( ABS(f) * surf%ol(m) * ri_crit ) )**-0.5
                   zeta = ABS(f) * -1.0_wp * z_offset_mcphee / ( eta_star * surf%us(m) )
                   tau = EXP( SQRT( im / ( kappa * xi_N ) ) * zeta )
                   surf%usws(m) = rho_ref_zw(k+1) * surf%us(m) *               &
                                  ( us_x * REAL(tau) - us_y * IMAG(tau) )
                   surf%vsws(m) = rho_ref_zw(k+1) * surf%us(m) *               &
                                  ( us_x * IMAG(tau) + us_y * REAL(tau) )
                
                ELSE

                   surf%usws(m) = surf%us(m) * rho_ref_zw(k+1) * kappa *       &
                                  ( u(k,j,i) - u(k+1,j,i) )                    &
                                  / LOG( surf%z_mo(m) / surf%z0(m) )
                   surf%vsws(m) = surf%us(m) * rho_ref_zw(k+1) * kappa *       &
                                  ( v(k,j,i) - v(k+1,j,i) )                    &
                                  / LOG( surf%z_mo(m) / surf%z0(m) )

                ENDIF
             ENDDO     
          ENDIF

!
!--       Compute the vertical kinematic heat flux
          IF ( trim(most_method) == 'mcphee' .AND. downward ) THEN
             
             !$OMP PARALLEL DO PRIVATE( i, j, k, s_factor )
             DO  m = 1, surf%ns  
   
                i = surf%i(m)            
                j = surf%j(m)
                k = surf%k(m)
                
                s_factor = -1.0_wp * surf%gamma_S(m) *                           &
                           ( surf%sa_io(m) - surf%sa1(m) ) / surf%sa_io(m)

                surf%shf(m)   = -1.0_wp * rho_ocean(k,j,i) *                     &
                                ( surf%gamma_T(m) + s_factor ) *                 &
                                ( surf%pt_io(m) - surf%pt1(m) )
                surf%sasws(m) = -1.0_wp * rho_ocean(k,j,i) *                     &
                                ( surf%gamma_S(m) + s_factor ) *                 &
                                ( surf%sa_io(m) - surf%sa1(m) )
                surf%melt(m)  = MAX( 0.0_wp,                                     &
                                     s_factor * rho_ocean(k,j,i)/1e3 )

                IF ( surf%melt(m) < melt_min ) CALL location_message('Melt rate < minimum',.TRUE.)
                IF ( surf%melt(m) > melt_max ) CALL location_message('Melt rate > maximum',.TRUE.)
             
             ENDDO
             
          ELSEIF (  .NOT.  constant_heatflux  .AND.  ( ( time_since_reference_point&
               <=  skip_time_do_lsm  .AND. simulated_time > 0.0_wp ) .OR.      &
               .NOT.  land_surface )  .AND.  .NOT. urban_surface )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, surf%ns 
                i    = surf%i(m)            
                j    = surf%j(m)
                k    = surf%k(m)
                surf%shf(m) = -surf%ts(m) * surf%us(m) * rho_ref_zw(k+ibit)
             ENDDO

          ENDIF
!
!--       Compute the vertical water flux
          IF (  .NOT.  constant_waterflux  .AND.  humidity  .AND.              &
               ( ( time_since_reference_point <= skip_time_do_lsm  .AND.       &
               simulated_time > 0.0_wp ) .OR.  .NOT.  land_surface  )  .AND.   &
               .NOT.  urban_surface  .AND.  .NOT. downward )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, surf%ns 
                i    = surf%i(m)            
                j    = surf%j(m)
                k    = surf%k(m)
                surf%qsws(m) = -surf%qs(m) * surf%us(m) * rho_ref_zw(k+ibit)
             ENDDO
          ENDIF
!
!--       Compute the vertical scalar flux
          IF (  .NOT.  constant_scalarflux  .AND.  passive_scalar  .AND.       &
                .NOT.  downward )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j )
             DO  m = 1, surf%ns   

                i    = surf%i(m)            
                j    = surf%j(m)
                surf%ssws(m) = -surf%ss(m) * surf%us(m)

             ENDDO
          ENDIF   
!
!--       Compute the vertical chemical species' flux
          DO  lsp = 1, nvar
             IF (  .NOT.  constant_csflux(lsp)  .AND.  air_chemistry  .AND.    &
                   .NOT.  downward )  THEN
                !$OMP PARALLEL DO PRIVATE( i, j )
                DO  m = 1, surf%ns   
                   i     = surf%i(m)            
                   j     = surf%j(m)
                   surf%cssws(lsp,m) = -surf%css(lsp,m) * surf%us(m)
                ENDDO
             ENDIF   
          ENDDO

!
!--       Compute (turbulent) fluxes of cloud water content and cloud drop conc.
          IF ( cloud_physics  .AND.  microphysics_morrison  .AND.              &
               .NOT. downward)  THEN
             !$OMP PARALLEL DO PRIVATE( i, j )
             DO  m = 1, surf%ns   

                i    = surf%i(m)            
                j    = surf%j(m)

                surf%qcsws(m) = -surf%qcs(m) * surf%us(m)
                surf%ncsws(m) = -surf%ncs(m) * surf%us(m)
             ENDDO
          ENDIF    
!
!--       Compute (turbulent) fluxes of rain water content and rain drop conc.
          IF ( cloud_physics  .AND.  microphysics_seifert  .AND.               &
               .NOT. downward)  THEN
             !$OMP PARALLEL DO PRIVATE( i, j )
             DO  m = 1, surf%ns   

                i    = surf%i(m)            
                j    = surf%j(m)

                surf%qrsws(m) = -surf%qrs(m) * surf%us(m)
                surf%nrsws(m) = -surf%nrs(m) * surf%us(m)
             ENDDO
          ENDIF

!
!--       Bottom boundary condition for the TKE. 
          IF ( ibc_e_b == 2 .AND. .NOT. downward )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, surf%ns   

                i    = surf%i(m)            
                j    = surf%j(m)
                k    = surf%k(m)

                e(k,j,i) = ( surf%us(m) / 0.1_wp )**2
!
!--             As a test: cm = 0.4
!               e(k,j,i) = ( us(j,i) / 0.4_wp )**2
                e(k-1,j,i)   = e(k,j,i)

             ENDDO
          ENDIF
!
!--       Top boundary condition for the TKE. 
          IF ( ibc_e_t == 2 .AND. downward )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, surf%ns   

                i    = surf%i(m)            
                j    = surf%j(m)
                k    = surf%k(m)

                e(k,j,i) = ( surf%us(m) / 0.1_wp )**2
                e(k+1,j,i)   = e(k,j,i)

             ENDDO
          ENDIF
!
!--    Calcuate surface fluxes at vertical surfaces. No stability is considered. 
       ELSE
!
!--       Compute usvs l={0,1} and vsus l={2,3}
          IF ( mom_uv )  THEN
!
!--          Generalize computation by introducing flags. At north- and south-
!--          facing surfaces u-component is used, at east- and west-facing
!--          surfaces v-component is used.
             flag_u = MERGE( 1.0_wp, 0.0_wp, l == 0  .OR.  l == 1 )   
             flag_v = MERGE( 1.0_wp, 0.0_wp, l == 2  .OR.  l == 3 )   
             !$OMP PARALLEL  DO PRIVATE( i, j, k, z_mo )
             DO  m = 1, surf%ns  
                i = surf%i(m)            
                j = surf%j(m)
                k = surf%k(m)

                z_mo = surf%z_mo(m)

                surf%mom_flux_uv(m) = kappa *                                  &
                                ( flag_u * u(k,j,i) + flag_v * v(k,j,i) )  /   &
                                                        LOG( z_mo / surf%z0(m) )

                surf%mom_flux_uv(m) = - surf%mom_flux_uv(m) * surf%us(m)

             ENDDO
          ENDIF
!
!--       Compute wsus l={0,1} and wsvs l={2,3}
          IF ( mom_w )  THEN
             !$OMP PARALLEL  DO PRIVATE( i, j, k, z_mo )
             DO  m = 1, surf%ns  
                i = surf%i(m)            
                j = surf%j(m)
                k = surf%k(m)

                z_mo = surf%z_mo(m)

                surf%mom_flux_w(m) = kappa * w(k,j,i) / LOG( z_mo / surf%z0(m) )

                surf%mom_flux_w(m) = - surf%mom_flux_w(m) * surf%us(m)

             ENDDO
          ENDIF
!
!--       Compute momentum fluxes used for subgrid-scale TKE production at 
!--       vertical surfaces. In constrast to the calculated momentum fluxes at 
!--       vertical surfaces before, which are defined on the u/v/w-grid, 
!--       respectively), the TKE fluxes are defined at the scalar grid. 
!--       
          IF ( mom_tke )  THEN
!
!--          Precalculate velocity components at scalar grid point. 
             ALLOCATE( u_i(1:surf%ns) )
             ALLOCATE( v_i(1:surf%ns) )
             ALLOCATE( w_i(1:surf%ns) )

             IF ( l == 0  .OR.  l == 1 )  THEN
                !$OMP PARALLEL DO PRIVATE( i, j, k )
                DO  m = 1, surf%ns  
                   i = surf%i(m)            
                   j = surf%j(m)
                   k = surf%k(m)

                   u_i(m) = 0.5_wp * ( u(k,j,i) + u(k,j,i+1) )
                   v_i(m) = 0.0_wp
                   w_i(m) = 0.5_wp * ( w(k,j,i) + w(k-1,j,i) )
                ENDDO
             ELSE
                !$OMP PARALLEL DO PRIVATE( i, j, k )
                DO  m = 1, surf%ns  
                   i = surf%i(m)            
                   j = surf%j(m)
                   k = surf%k(m)

                   u_i(m) = 0.0_wp
                   v_i(m) = 0.5_wp * ( v(k,j,i) + v(k,j+1,i) )
                   w_i(m) = 0.5_wp * ( w(k,j,i) + w(k-1,j,i) )
                ENDDO
             ENDIF

             !$OMP PARALLEL DO PRIVATE( i, j, dum, z_mo )
             DO  m = 1, surf%ns  
                i = surf%i(m)            
                j = surf%j(m)

                z_mo = surf%z_mo(m)

                dum = kappa / LOG( z_mo / surf%z0(m) )
!
!--             usvs (l=0,1) and vsus (l=2,3)
                surf%mom_flux_tke(0,m) = dum * ( u_i(m) + v_i(m) )
!
!--             wsvs (l=0,1) and wsus (l=2,3)
                surf%mom_flux_tke(1,m) = dum * w_i(m)

                surf%mom_flux_tke(0:1,m) =                                     &
                               - surf%mom_flux_tke(0:1,m) * surf%us(m)
             ENDDO
!
!--          Deallocate temporary arrays
             DEALLOCATE( u_i )             
             DEALLOCATE( v_i )  
             DEALLOCATE( w_i )  
          ENDIF
       ENDIF

    END SUBROUTINE calc_surface_fluxes


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate surface fluxes (shf, sasws) and ice-ocean interface pt, sa using
!> 3 equations as presented in Asay-Davis et al. (2016) Equations (24)-(26).
!> A variant on Equation (25), the freezing point equation, is used:
!> pt_io = pt_io_p + dptf_dsa*(sa_io - sa_io_p) 
!> Heat and salt transport across the surface layer are parameterized following
!> McPhee et al. (1987) and McPhee et al. (2008)
!> Note: technically, temperature gradients should be in terms of in situ
!> temperature rather than the present method, which uses potential temperature
!> below the boundary and conservative temperature freezing point at the 
!> boundary. For the range of P,SA conditions on the continental shelf, 
!> this difference in temperature is on the order of 0.007 K.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_3eqn

       IMPLICIT NONE

       INTEGER(iwp)  ::  m                   !< loop variable over all horizontal surf elements
       INTEGER(iwp)  ::  n                   !< iteration counter
       INTEGER(iwp)  ::  nn = 0              !< iteration counter
       INTEGER(iwp)  ::  max_iter = 10       !< maximum number of iterations

       REAL(wp) ::  a, b, c                  !< coefficients of quadratic equation
       REAL(wp) ::  dptf_dsa                 !< derivative of conservative temperature 
                                             !< at freezing point with respect to salinity
       REAL(wp) ::  sa_io_k                  !< change in salinity over kth iteration
       REAL(wp) ::  sa_io_p, pt_io_p         !< interface properties at previous timestep
       REAL(wp) ::  sa_io_min, sa_io_max     !< acceptable range of values
       REAL(wp) ::  dsa_tol = 0.001          !< minimum acceptable change in salinity over kth iteration
       REAL(wp) ::  h_nu                     !< thickness of viscous sublayer
       REAL(wp) ::  term1                    !< term in 3 equations
       REAL(wp) ::  Gamma_turb               !< turbulent contribution to exchange velocity
       REAL(wp) ::  Gamma_mol_T, Gamma_mol_S !< molecular contribution to exchange velocity

       sa_io_min = 0.0_wp
       sa_io_max = 40.0_wp

       DO m = 1, surf%ns

          !$OMP PARALLEL DO PRIVATE( i, j, k )
          i = surf%i(m)            
          j = surf%j(m)
          k = surf%k(m)

!--       Store sa, pt at the boundary at the previous time step
          sa_io_p = surf%sa_io(m)

          dptf_dsa = pt_freezing_SA(surface_pressure/1e2,sa_io_p)
          pt_io_p  = pt_freezing(surface_pressure/1e2,sa_io_p)

!--       Loop until a grid cell's interface salinity changes by less than the
!--       tolerance or until maximum number of iterations is reached
          DO n = 1, max_iter

             sa_io_k = surf%sa_io(m)
             
!--          Update Obukhov length each iteration because it depends on the
!--          buoyancy flux  
             CALL calc_ol_m(m)

!--          If buoyancy flux is 0, set stability parameter at neutral value
             IF (surf%ol(m) <= 0.0_wp) THEN
                eta_star = 1.0_wp
!
!--             Buoyancy flux is destabilizing, so assume double diffusion is not relevant
!--             following the recommendation of McPhee et al. (2008)
                surf%gamma_T(m) = 0.0057_wp * surf%us(m)
                surf%gamma_S(m) = 0.0057_wp * surf%us(m)

!--          If Obukhov length is positive, then buoyancy flux is stabilizing
!--          and heat and salt transfer are governed by McPhee et al. (1987) 
!--          parameterization
             ELSE

                ! viscous sublayer thickness, Tennekes and Lumley (1972) p. 160
                h_nu = 5.0_wp * molecular_viscosity / surf%us(m)
                
                eta_star = ( 1.0_wp + ( xi_N * surf%us(m) ) /                  &
                           ( ABS(f) * surf%ol(m) * ri_crit ) )**-0.5
                Gamma_turb = ( 1 / kappa ) *                                   &
                                    LOG( ( surf%us(m) * xi_N * eta_star**2 ) / &
                                         ( ABS(f) * h_nu ) )                   &
                             + ( 1 / ( 2 * xi_N * eta_star ) )                 &
                             - ( 1 / kappa )
                Gamma_mol_T = 12.5_wp * prandtl_number**(2/3) - 6
                Gamma_mol_S = 12.5_wp * schmidt_number**(2/3) - 6

                surf%gamma_T(m) = surf%us(m) / (Gamma_turb + Gamma_mol_T)
                surf%gamma_S(m) = surf%us(m) / (Gamma_turb + Gamma_mol_S)

             ENDIF

!
!--          Solve 3 equations for salinity at ice-ocean interface
             term1 = pt_io_p - ( dptf_dsa * sa_io_k ) - surf%pt1(m)
             a = dptf_dsa * surf%gamma_T(m)
             b = ( surf%gamma_T(m) * term1 ) - ( surf%gamma_S(m) * ( l_m / cpw ) ) 
             c = surf%gamma_S(m) * surf%sa1(m) * ( l_m / cpw )

             surf%sa_io(m) = 2.0_wp*c / (-1.0_wp*b + SQRT(b**2 - 4.0_wp*a*c))
             IF ( surf%sa_io(m) < sa_io_min )                                  &
                CALL location_message('Interface salinity is below minimum value',.TRUE.)
             IF ( surf%sa_io(m) > sa_io_max )                                  &
                CALL location_message('Interface salinity exceeds minimum value',.TRUE.)

!
!--          Calculate approximate freezing temperature at ice-ocean interface
             surf%pt_io(m) = pt_io_p + ( dptf_dsa * ( surf%sa_io(m) - sa_io_p ) )

             nn = nn + 1

!
!--          Stop iterations if interface salinity change is small enough
             IF ( ABS(surf%sa_io(m) - sa_io_k) < dsa_tol ) EXIT

          ENDDO ! end iterations

       ENDDO ! end loop over surface indices

!
!--    Update surface fluxes
       downward = .TRUE.
       CALL calc_surface_fluxes
       downward = .FALSE.

       write(message_string,*) 'us = ',surf%us(1)
       CALL location_message(message_string,.TRUE.)
       write(message_string,*) 'melt = ',surf%melt(1)
       CALL location_message(message_string,.TRUE.)
       IF ( isnan(surf%melt(1)) ) THEN
          CALL message( 'surface_layer_fluxes', 'PA0655', 1, 2, 0, 6, 0 )
       ENDIF
       nn = nn/surf%ns
       WRITE(message_string,*) 'Average 3 eqn iterations = ', nn
       CALL location_message(message_string,.TRUE.)
       
    END SUBROUTINE calc_3eqn

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Integrated stability function for momentum
!------------------------------------------------------------------------------!
    FUNCTION psi_m( zeta ) 
       
       USE kinds

       IMPLICIT NONE 

       REAL(wp)            :: psi_m !< Integrated similarity function result
       REAL(wp)            :: zeta  !< Stability parameter z/L
       REAL(wp)            :: x     !< dummy variable

       REAL(wp), PARAMETER :: a = 1.0_wp            !< constant
       REAL(wp), PARAMETER :: b = 0.66666666666_wp  !< constant
       REAL(wp), PARAMETER :: c = 5.0_wp            !< constant
       REAL(wp), PARAMETER :: d = 0.35_wp           !< constant
       REAL(wp), PARAMETER :: c_d_d = c / d         !< constant
       REAL(wp), PARAMETER :: bc_d_d = b * c / d    !< constant


       IF ( zeta < 0.0_wp )  THEN
          x = SQRT( SQRT( 1.0_wp  - 16.0_wp * zeta ) )
          psi_m = pi * 0.5_wp - 2.0_wp * ATAN( x ) + LOG( ( 1.0_wp + x )**2    &
                  * ( 1.0_wp + x**2 ) * 0.125_wp )
       ELSE

          psi_m = - b * ( zeta - c_d_d ) * EXP( -d * zeta ) - a * zeta         &
                   - bc_d_d
!
!--       Old version for stable conditions (only valid for z/L < 0.5)
!--       psi_m = - 5.0_wp * zeta

       ENDIF

    END FUNCTION psi_m

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Integrated stability function for heat and moisture
!------------------------------------------------------------------------------!
    FUNCTION psi_h( zeta ) 
       
       USE kinds

       IMPLICIT NONE 

       REAL(wp)            :: psi_h !< Integrated similarity function result
       REAL(wp)            :: zeta  !< Stability parameter z/L
       REAL(wp)            :: x     !< dummy variable

       REAL(wp), PARAMETER :: a = 1.0_wp            !< constant
       REAL(wp), PARAMETER :: b = 0.66666666666_wp  !< constant
       REAL(wp), PARAMETER :: c = 5.0_wp            !< constant
       REAL(wp), PARAMETER :: d = 0.35_wp           !< constant
       REAL(wp), PARAMETER :: c_d_d = c / d         !< constant
       REAL(wp), PARAMETER :: bc_d_d = b * c / d    !< constant


       IF ( zeta < 0.0_wp )  THEN
          x = SQRT( 1.0_wp  - 16.0_wp * zeta )
          psi_h = 2.0_wp * LOG( (1.0_wp + x ) / 2.0_wp )
       ELSE
          psi_h = - b * ( zeta - c_d_d ) * EXP( -d * zeta ) - (1.0_wp          &
                  + 0.66666666666_wp * a * zeta )**1.5_wp - bc_d_d             &
                  + 1.0_wp
!
!--       Old version for stable conditions (only valid for z/L < 0.5)
!--       psi_h = - 5.0_wp * zeta
       ENDIF

    END FUNCTION psi_h

 END MODULE surface_layer_fluxes_mod
