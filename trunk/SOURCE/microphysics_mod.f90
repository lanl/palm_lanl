!> @file microphysics_mod.f90
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
! ------------------
! 
! 
! Former revisions:
! -----------------
! $Id: microphysics_mod.f90 3026 2018-05-22 10:30:53Z schwenkel $
! Changed the name specific humidity to mixing ratio, since we are computing
! mixing ratios.
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2701 2017-12-15 15:40:50Z suehring
! Changes from last commit documented
! 
! 2698 2017-12-14 18:46:24Z suehring
! Bugfix in get_topography_top_index
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2608 2017-11-13 14:04:26Z schwenkel
! Calculation of supersaturation in external module (diagnostic_quantities_mod).
! Change: correct calculation of saturation specific humidity to saturation
! mixing ratio (the factor of 0.378 vanishes).
! 
! 2522 2017-10-05 14:20:37Z schwenkel
! Minor bugfix
! 
! 2375 2017-08-29 14:10:28Z schwenkel
! Improved aerosol initilization and some minor bugfixes 
! for droplet sedimenation
! 
! 2318 2017-07-20 17:27:44Z suehring
! Get topography top index via Function call 
!
! 2317 2017-07-20 17:27:19Z suehring
! s1 changed to log_sigma
!
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison'
! includes two more prognostic equations for cloud drop concentration (nc)
! and cloud water content (qc).
! - The process of activation is parameterized with a simple Twomey
!   activion scheme or with considering solution and curvature
!   effects (Khvorostyanov and Curry ,2006).
! - The saturation adjustment scheme is replaced by the parameterization
!   of condensation rates (Khairoutdinov and Kogan, 2000, Mon. Wea. Rev.,128).
! - All other microphysical processes of Seifert and Beheng are used.
!   Additionally, in those processes the reduction of cloud number concentration
!   is considered.
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography and surface concept
!
! 2155 2017-02-21 09:57:40Z hoffmann
! Bugfix in the calculation of microphysical quantities on ghost points.
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
! Adapted for modularization of microphysics.
!
! 1845 2016-04-08 08:29:13Z raasch
! nzb_2d replaced by nzb_s_inner, Kessler precipitation is stored at surface
! point (instead of one point above surface)
!
! 1831 2016-04-07 13:15:51Z hoffmann
! turbulence renamed collision_turbulence,
! drizzle renamed cloud_water_sedimentation. cloud_water_sedimentation also
! avaialble for microphysics_kessler.
!
! 1822 2016-04-07 07:49:42Z hoffmann
! Unused variables removed.
! Kessler scheme integrated.
!
! 1691 2015-10-26 16:17:44Z maronga
! Added new routine calc_precipitation_amount. The routine now allows to account
! for precipitation due to sedimenation of cloud (fog) droplets
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1646 2015-09-02 16:00:10Z hoffmann
! Bugfix: Wrong computation of d_mean.
!
! 1361 2014-04-16 15:17:48Z hoffmann
! Bugfix in sedimentation_rain: Index corrected.
! Vectorized version of adjust_cloud added.
! Little reformatting of the code.
!
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1346 2014-03-27 13:18:20Z heinze
! Bugfix: REAL constants provided with KIND-attribute especially in call of
! intrinsic function like MAX, MIN, SIGN
!
! 1334 2014-03-25 12:21:40Z heinze
! Bugfix: REAL constants provided with KIND-attribute
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL constants defined as wp-kind
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1241 2013-10-30 11:36:58Z heinze
! hyp and rho_ocean have to be calculated at each time step if data from external
! file LSF_DATA are used
!
! 1115 2013-03-26 18:16:16Z hoffmann
! microphyical tendencies are calculated in microphysics_control in an optimized
! way; unrealistic values are prevented; bugfix in evaporation; some reformatting
!
! 1106 2013-03-04 05:31:38Z raasch
! small changes in code formatting
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
! file put under GPL
!
! 1065 2012-11-22 17:42:36Z hoffmann
! Sedimentation process implemented according to Stevens and Seifert (2008).
! Turbulence effects on autoconversion and accretion added (Seifert, Nuijens
! and Stevens, 2010).
!
! 1053 2012-11-13 17:11:03Z hoffmann
! initial revision
!
! Description:
! ------------
!> Calculate bilk cloud microphysics.
!------------------------------------------------------------------------------!
 MODULE microphysics_mod

    USE  kinds

    IMPLICIT NONE

    LOGICAL ::  cloud_water_sedimentation = .FALSE.       !< cloud water sedimentation
    LOGICAL ::  curvature_solution_effects_bulk = .FALSE. !< flag for considering koehler theory
    LOGICAL ::  limiter_sedimentation = .TRUE.            !< sedimentation limiter
    LOGICAL ::  collision_turbulence = .FALSE.            !< turbulence effects
    LOGICAL ::  ventilation_effect = .TRUE.               !< ventilation effect

    REAL(wp) ::  a_1 = 8.69E-4_wp          !< coef. in turb. parametrization (cm-2 s3)
    REAL(wp) ::  a_2 = -7.38E-5_wp         !< coef. in turb. parametrization (cm-2 s3)
    REAL(wp) ::  a_3 = -1.40E-2_wp         !< coef. in turb. parametrization
    REAL(wp) ::  a_term = 9.65_wp          !< coef. for terminal velocity (m s-1)
    REAL(wp) ::  a_vent = 0.78_wp          !< coef. for ventilation effect
    REAL(wp) ::  b_1 = 11.45E-6_wp         !< coef. in turb. parametrization (m)
    REAL(wp) ::  b_2 = 9.68E-6_wp          !< coef. in turb. parametrization (m)
    REAL(wp) ::  b_3 = 0.62_wp             !< coef. in turb. parametrization
    REAL(wp) ::  b_term = 9.8_wp           !< coef. for terminal velocity (m s-1)
    REAL(wp) ::  b_vent = 0.308_wp         !< coef. for ventilation effect
    REAL(wp) ::  beta_cc = 3.09E-4_wp      !< coef. in turb. parametrization (cm-2 s3)
    REAL(wp) ::  c_1 = 4.82E-6_wp          !< coef. in turb. parametrization (m)
    REAL(wp) ::  c_2 = 4.8E-6_wp           !< coef. in turb. parametrization (m)
    REAL(wp) ::  c_3 = 0.76_wp             !< coef. in turb. parametrization
    REAL(wp) ::  c_const = 0.93_wp         !< const. in Taylor-microscale Reynolds number
    REAL(wp) ::  c_evap = 0.7_wp           !< constant in evaporation
    REAL(wp) ::  c_term = 600.0_wp         !< coef. for terminal velocity (m-1)
    REAL(wp) ::  diff_coeff_l = 0.23E-4_wp  !< diffusivity of water vapor (m2 s-1)
    REAL(wp) ::  eps_sb = 1.0E-10_wp       !< threshold in two-moments scheme
    REAL(wp) ::  eps_mr = 0.0_wp           !< threshold for morrison scheme
    REAL(wp) ::  k_cc = 9.44E09_wp         !< const. cloud-cloud kernel (m3 kg-2 s-1)
    REAL(wp) ::  k_cr0 = 4.33_wp           !< const. cloud-rain kernel (m3 kg-1 s-1)
    REAL(wp) ::  k_rr = 7.12_wp            !< const. rain-rain kernel (m3 kg-1 s-1)
    REAL(wp) ::  k_br = 1000.0_wp          !< const. in breakup parametrization (m-1)
    REAL(wp) ::  k_st = 1.2E8_wp           !< const. in drizzle parametrization (m-1 s-1)
    REAL(wp) ::  kappa_rr = 60.7_wp        !< const. in collision kernel (kg-1/3)
    REAL(wp) ::  kin_vis_air = 1.4086E-5_wp  !< kin. viscosity of air (m2 s-1)
    REAL(wp) ::  prec_time_const = 0.001_wp  !< coef. in Kessler scheme (s-1)
    REAL(wp) ::  ql_crit = 0.0005_wp       !< coef. in Kessler scheme (kg kg-1)
    REAL(wp) ::  schmidt_p_1d3=0.8921121_wp  !< Schmidt number**0.33333, 0.71**0.33333
    REAL(wp) ::  sigma_gc = 1.3_wp         !< geometric standard deviation cloud droplets
    REAL(wp) ::  thermal_conductivity_l = 2.43E-2_wp  !< therm. cond. air (J m-1 s-1 K-1)
    REAL(wp) ::  w_precipitation = 9.65_wp  !< maximum terminal velocity (m s-1)
    REAL(wp) ::  x0 = 2.6E-10_wp           !< separating drop mass (kg)
    REAL(wp) ::  xamin = 5.24E-19_wp       !< average aerosol mass (kg) (~ 0.05µm)
    REAL(wp) ::  xcmin = 4.18E-15_wp       !< minimum cloud drop size (kg) (~ 1µm)
    REAL(wp) ::  xcmax = 2.6E-10_wp        !< maximum cloud drop size (kg) (~ 40µm)
    REAL(wp) ::  xrmin = 2.6E-10_wp        !< minimum rain drop size (kg)
    REAL(wp) ::  xrmax = 5.0E-6_wp         !< maximum rain drop site (kg)

    REAL(wp) ::  c_sedimentation = 2.0_wp        !< Courant number of sedimentation process
    REAL(wp) ::  dpirho_l                        !< 6.0 / ( pi * rho_l )
    REAL(wp) ::  dry_aerosol_radius = 0.05E-6_wp !< dry aerosol radius
    REAL(wp) ::  dt_micro                        !< microphysics time step
    REAL(wp) ::  sigma_bulk = 2.0_wp             !< width of aerosol spectrum
    REAL(wp) ::  na_init = 100.0E6_wp            !< Total particle/aerosol concentration (cm-3)
    REAL(wp) ::  nc_const = 70.0E6_wp            !< cloud droplet concentration
    REAL(wp) ::  dt_precipitation = 100.0_wp     !< timestep precipitation (s)
    REAL(wp) ::  sed_qc_const                    !< const. for sedimentation of cloud water
    REAL(wp) ::  pirho_l                         !< pi * rho_l / 6.0;

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  nc_1d  !<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  nr_1d  !<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_1d  !<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  qc_1d  !<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  qr_1d  !<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  q_1d   !<

    SAVE

    PRIVATE
    PUBLIC microphysics_control, microphysics_init

    PUBLIC cloud_water_sedimentation, collision_turbulence,                    &
           curvature_solution_effects_bulk, c_sedimentation,                   &
           dry_aerosol_radius, dt_precipitation,                               &
           limiter_sedimentation, na_init, nc_const, sigma_gc, sigma_bulk,     &
           ventilation_effect


    INTERFACE microphysics_control
       MODULE PROCEDURE microphysics_control
       MODULE PROCEDURE microphysics_control_ij
    END INTERFACE microphysics_control

    INTERFACE adjust_cloud
       MODULE PROCEDURE adjust_cloud
       MODULE PROCEDURE adjust_cloud_ij
    END INTERFACE adjust_cloud

    INTERFACE activation
       MODULE PROCEDURE activation
       MODULE PROCEDURE activation_ij
    END INTERFACE activation

    INTERFACE condensation
       MODULE PROCEDURE condensation
       MODULE PROCEDURE condensation_ij
    END INTERFACE condensation

    INTERFACE autoconversion
       MODULE PROCEDURE autoconversion
       MODULE PROCEDURE autoconversion_ij
    END INTERFACE autoconversion

    INTERFACE autoconversion_kessler
       MODULE PROCEDURE autoconversion_kessler
       MODULE PROCEDURE autoconversion_kessler_ij
    END INTERFACE autoconversion_kessler

    INTERFACE accretion
       MODULE PROCEDURE accretion
       MODULE PROCEDURE accretion_ij
    END INTERFACE accretion

    INTERFACE selfcollection_breakup
       MODULE PROCEDURE selfcollection_breakup
       MODULE PROCEDURE selfcollection_breakup_ij
    END INTERFACE selfcollection_breakup

    INTERFACE evaporation_rain
       MODULE PROCEDURE evaporation_rain
       MODULE PROCEDURE evaporation_rain_ij
    END INTERFACE evaporation_rain

    INTERFACE sedimentation_cloud
       MODULE PROCEDURE sedimentation_cloud
       MODULE PROCEDURE sedimentation_cloud_ij
    END INTERFACE sedimentation_cloud

    INTERFACE sedimentation_rain
       MODULE PROCEDURE sedimentation_rain
       MODULE PROCEDURE sedimentation_rain_ij
    END INTERFACE sedimentation_rain

    INTERFACE calc_precipitation_amount
       MODULE PROCEDURE calc_precipitation_amount
       MODULE PROCEDURE calc_precipitation_amount_ij
    END INTERFACE calc_precipitation_amount

 CONTAINS
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of bulk microphysics
!------------------------------------------------------------------------------!
    SUBROUTINE microphysics_init

       USE arrays_3d,                                                          &
           ONLY:  dzu

       USE constants,                                                          &
           ONLY:  pi

       USE cloud_parameters,                                                   &
           ONLY:  molecular_weight_of_solute, rho_l, rho_s, vanthoff

       USE control_parameters,                                                 &
           ONLY:  aerosol_nacl, aerosol_c3h4o4 , aerosol_nh4no3,               &
                  microphysics_morrison, microphysics_seifert

       USE indices,                                                            &
           ONLY:  nzb, nzt

       IMPLICIT NONE

!
!--    constant for the sedimentation of cloud water (2-moment cloud physics)
       sed_qc_const = k_st * ( 3.0_wp / ( 4.0_wp * pi * rho_l )                &
                          )**( 2.0_wp / 3.0_wp ) *                             &
                      EXP( 5.0_wp * LOG( sigma_gc )**2 )

!
!--    Calculate timestep according to precipitation
       IF ( microphysics_seifert )  THEN
          dt_precipitation = c_sedimentation * MINVAL( dzu(nzb+2:nzt) ) /      &
                             w_precipitation
       ENDIF

!
!--    Set constants for certain aerosol type
       IF ( microphysics_morrison )  THEN
          IF ( aerosol_nacl ) THEN
             molecular_weight_of_solute = 0.05844_wp 
             rho_s                      = 2165.0_wp
             vanthoff                   = 2.0_wp
          ELSEIF ( aerosol_c3h4o4 ) THEN
             molecular_weight_of_solute = 0.10406_wp 
             rho_s                      = 1600.0_wp
             vanthoff                   = 1.37_wp
          ELSEIF ( aerosol_nh4no3 ) THEN
             molecular_weight_of_solute = 0.08004_wp 
             rho_s                      = 1720.0_wp
             vanthoff                   = 2.31_wp
          ENDIF
       ENDIF
 
!
!--    Pre-calculate frequently calculated fractions of pi and rho_l
       pirho_l  = pi * rho_l / 6.0_wp
       dpirho_l = 1.0_wp / pirho_l

!
!--    Allocate 1D microphysics arrays
       ALLOCATE ( pt_1d(nzb:nzt+1), q_1d(nzb:nzt+1),         &
                  qc_1d(nzb:nzt+1) )

       IF ( microphysics_morrison .OR. microphysics_seifert ) THEN
          ALLOCATE ( nc_1d(nzb:nzt+1) )
       ENDIF

       IF ( microphysics_seifert )  THEN
          ALLOCATE ( nr_1d(nzb:nzt+1), qr_1d(nzb:nzt+1) )
       ENDIF

    END SUBROUTINE microphysics_init


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Control of microphysics for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE microphysics_control

       USE arrays_3d,                                                          &
           ONLY:  hyp, pt_init, prr, zu

       USE cloud_parameters,                                                   &
           ONLY:  cp, hyrho, pt_d_t, r_d, t_d_pt

       USE control_parameters,                                                 &
           ONLY:  call_microphysics_at_all_substeps, dt_3d, g,                 &
                  intermediate_timestep_count, large_scale_forcing,            &
                  lsf_surf, microphysics_kessler, microphysics_morrison,       &
                  microphysics_seifert, pt_surface,                            &
                  rho_surface,surface_pressure

       USE indices,                                                            &
           ONLY:  nzb, nzt

       USE kinds

       USE statistics,                                                         &
           ONLY:  weight_pres

       IMPLICIT NONE

       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  t_surface         !<

       IF ( large_scale_forcing  .AND.  lsf_surf ) THEN
!
!--       Calculate:
!--       pt / t : ratio of potential and actual temperature (pt_d_t)
!--       t / pt : ratio of actual and potential temperature (t_d_pt)
!--       p_0(z) : vertical profile of the hydrostatic pressure (hyp)
          t_surface = pt_surface * ( surface_pressure / 1000.0_wp )**0.286_wp
          DO  k = nzb, nzt+1
             hyp(k)    = surface_pressure * 100.0_wp * &
                         ( ( t_surface - g / cp * zu(k) ) /                    &
                         t_surface )**(1.0_wp / 0.286_wp)
             pt_d_t(k) = ( 100000.0_wp / hyp(k) )**0.286_wp
             t_d_pt(k) = 1.0_wp / pt_d_t(k)
             hyrho(k)  = hyp(k) / ( r_d * t_d_pt(k) * pt_init(k) )
          ENDDO

!
!--       Compute reference density
          rho_surface = surface_pressure * 100.0_wp / ( r_d * t_surface )
       ENDIF

!
!--    Compute length of time step
       IF ( call_microphysics_at_all_substeps )  THEN
          dt_micro = dt_3d * weight_pres(intermediate_timestep_count)
       ELSE
          dt_micro = dt_3d
       ENDIF

!
!--    Reset precipitation rate
       IF ( intermediate_timestep_count == 1 )  prr = 0.0_wp

!
!--    Compute cloud physics
       IF ( microphysics_kessler )  THEN

          CALL autoconversion_kessler
          IF ( cloud_water_sedimentation )  CALL sedimentation_cloud

       ELSEIF ( microphysics_seifert )  THEN

          CALL adjust_cloud
          IF ( microphysics_morrison )  CALL activation
          IF ( microphysics_morrison )  CALL condensation
          CALL autoconversion
          CALL accretion
          CALL selfcollection_breakup
          CALL evaporation_rain
          CALL sedimentation_rain
          IF ( cloud_water_sedimentation )  CALL sedimentation_cloud

       ENDIF

       CALL calc_precipitation_amount

    END SUBROUTINE microphysics_control

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Adjust number of raindrops to avoid nonlinear effects in sedimentation and
!> evaporation of rain drops due to too small or too big weights
!> of rain drops (Stevens and Seifert, 2008).
!------------------------------------------------------------------------------!
    SUBROUTINE adjust_cloud

       USE arrays_3d,                                                          &
           ONLY:  qc, nc, qr, nr

       USE cloud_parameters,                                                   &
           ONLY:  hyrho

       USE control_parameters,                                                 &
           ONLY:  microphysics_morrison

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nyng, nysg, nzb, nzt, wall_flags_0

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       CALL cpu_log( log_point_s(54), 'adjust_cloud', 'start' )

       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb+1, nzt
                IF ( qr(k,j,i) <= eps_sb )  THEN
                   qr(k,j,i) = 0.0_wp
                   nr(k,j,i) = 0.0_wp
                ELSE
                   IF ( nr(k,j,i) * xrmin > qr(k,j,i) * hyrho(k) )  THEN
                      nr(k,j,i) = qr(k,j,i) * hyrho(k) / xrmin *               &
                                       MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,i), 0 ) )
                   ELSEIF ( nr(k,j,i) * xrmax < qr(k,j,i) * hyrho(k) )  THEN
                      nr(k,j,i) = qr(k,j,i) * hyrho(k) / xrmax *               &
                                       MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,i), 0 ) )
                   ENDIF
                ENDIF
                IF ( microphysics_morrison ) THEN
                   IF ( qc(k,j,i) <= eps_sb )  THEN
                      qc(k,j,i) = 0.0_wp
                      nc(k,j,i) = 0.0_wp
                   ELSE
                      IF ( nc(k,j,i) * xcmin > qc(k,j,i) * hyrho(k) )  THEN
                         nc(k,j,i) = qc(k,j,i) * hyrho(k) / xcmin *            &
                                          MERGE( 1.0_wp, 0.0_wp,               &
                                              BTEST( wall_flags_0(k,j,i), 0 ) )
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point_s(54), 'adjust_cloud', 'stop' )

    END SUBROUTINE adjust_cloud

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate number of activated condensation nucleii after simple activation
!> scheme of Twomey, 1959.
!------------------------------------------------------------------------------!
    SUBROUTINE activation

       USE arrays_3d,                                                          &
           ONLY:  hyp, nc, nr, pt, q,  qc, qr

       USE cloud_parameters,                                                   &
           ONLY:  hyrho, l_d_cp, l_d_r, l_v, molecular_weight_of_solute,       &
                  molecular_weight_of_water, rho_l, rho_s, r_v, t_d_pt,        &
                  vanthoff

       USE constants,                                                          &
           ONLY:  pi

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE diagnostic_quantities_mod,                                          &
           ONLY: e_s, magnus, q_s, sat, supersaturation, t_l

       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nysg, nyng, nzb, nzt

       USE kinds

       USE control_parameters,                                                 &
          ONLY: simulated_time

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  activ             !<
       REAL(wp)     ::  afactor           !<
       REAL(wp)     ::  beta_act          !<
       REAL(wp)     ::  bfactor           !<
       REAL(wp)     ::  k_act             !<
       REAL(wp)     ::  n_act             !<
       REAL(wp)     ::  n_ccn             !<
       REAL(wp)     ::  s_0               !<
       REAL(wp)     ::  sat_max           !<
       REAL(wp)     ::  sigma             !<
       REAL(wp)     ::  sigma_act         !<

       CALL cpu_log( log_point_s(65), 'activation', 'start' )

       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb+1, nzt

!
!--             Call calculation of supersaturation located 
!--             in diagnostic_quantities_mod
                CALL supersaturation ( i, j, k )
!
!--             Prescribe parameters for activation
!--             (see: Bott + Trautmann, 2002, Atm. Res., 64)
                k_act  = 0.7_wp
                activ  = 0.0_wp


                IF ( sat > 0.0 .AND. .NOT. curvature_solution_effects_bulk ) THEN
!
!--                Compute the number of activated Aerosols
!--                (see: Twomey, 1959, Pure and applied Geophysics, 43)
                   n_act     = na_init * sat**k_act
!
!--                Compute the number of cloud droplets
!--                (see: Morrison + Grabowski, 2007, JAS, 64)
!                  activ = MAX( n_act - nc(k,j,i), 0.0_wp) / dt_micro

!
!--                Compute activation rate after Khairoutdinov and Kogan
!--                (see: Khairoutdinov + Kogan, 2000, Mon. Wea. Rev., 128)
                   sat_max = 1.0_wp / 100.0_wp
                   activ   = MAX( 0.0_wp, ( (na_init + nc(k,j,i) ) * MIN       &
                      ( 1.0_wp, ( sat / sat_max )**k_act) - nc(k,j,i) ) ) /    &
                       dt_micro
                ELSEIF ( sat > 0.0 .AND. curvature_solution_effects_bulk ) THEN
!
!--                Curvature effect (afactor) with surface tension
!--                parameterization by Straka (2009)
                   sigma = 0.0761_wp - 0.000155_wp * ( t_l - 273.15_wp )
                   afactor = 2.0_wp * sigma / ( rho_l * r_v * t_l )
!
!--                Solute effect (bfactor)
                   bfactor = vanthoff * molecular_weight_of_water *            &
                       rho_s / ( molecular_weight_of_solute * rho_l )

!
!--                Prescribe power index that describes the soluble fraction
!--                of an aerosol particle (beta) 
!--                (see: Morrison + Grabowski, 2007, JAS, 64)
                   beta_act  = 0.5_wp
                   sigma_act = sigma_bulk**( 1.0_wp + beta_act )
!
!--                Calculate mean geometric supersaturation (s_0) with
!--                parameterization by Khvorostyanov and Curry (2006)
                   s_0   = dry_aerosol_radius **(- ( 1.0_wp + beta_act ) ) *    &
                       ( 4.0_wp * afactor**3 / ( 27.0_wp * bfactor ) )**0.5_wp

!
!--                Calculate number of activated CCN as a function of
!--                supersaturation and taking Koehler theory into account
!--                (see: Khvorostyanov + Curry, 2006, J. Geo. Res., 111)
                   n_ccn = ( na_init / 2.0_wp ) * ( 1.0_wp - ERF(              &
                      LOG( s_0 / sat ) / ( SQRT(2.0_wp) * LOG(sigma_act) ) ) )
                   activ = MAX( ( n_ccn - nc(k,j,i) ) / dt_micro, 0.0_wp )
                ENDIF

                nc(k,j,i) = MIN( (nc(k,j,i) + activ * dt_micro), na_init)

             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point_s(65), 'activation', 'stop' )

    END SUBROUTINE activation


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate condensation rate for cloud water content (after Khairoutdinov and
!> Kogan, 2000).
!------------------------------------------------------------------------------!
    SUBROUTINE condensation

       USE arrays_3d,                                                          &
           ONLY:  hyp, nr, pt, q,  qc, qr, nc

       USE cloud_parameters,                                                   &
           ONLY:  hyrho, l_d_cp, l_d_r, l_v, r_v, t_d_pt

       USE constants,                                                          &
           ONLY:  pi

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE diagnostic_quantities_mod,                                          &
           ONLY: e_s, magnus, q_s, sat, supersaturation, t_l

       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nysg, nyng, nzb, nzt

       USE kinds

       USE control_parameters,                                                 &
          ONLY: simulated_time

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  cond              !<
       REAL(wp)     ::  cond_max          !<
       REAL(wp)     ::  dc                !<
       REAL(wp)     ::  evap              !<
       REAL(wp)     ::  evap_nc           !<
       REAL(wp)     ::  g_fac             !<
       REAL(wp)     ::  nc_0              !<
       REAL(wp)     ::  temp              !<
       REAL(wp)     ::  xc                !<

       CALL cpu_log( log_point_s(66), 'condensation', 'start' )

       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb+1, nzt
!
!--             Call calculation of supersaturation located 
!--             in diagnostic_quantities_mod
                CALL supersaturation ( i, j, k )
!
!--             Actual temperature:
                temp = t_l + l_d_cp * ( qc(k,j,i) + qr(k,j,i) )

                g_fac  = 1.0_wp / ( ( l_v / ( r_v * temp ) - 1.0_wp ) *        &
                                    l_v / ( thermal_conductivity_l * temp )    &
                                    + r_v * temp / ( diff_coeff_l * e_s )      &
                                    )
!
!--             Mean weight of cloud drops
                IF ( nc(k,j,i) <= 0.0_wp) CYCLE
                xc = MAX( (hyrho(k) * qc(k,j,i) / nc(k,j,i)), xcmin)
!
!--             Weight averaged diameter of cloud drops:
                dc   = ( xc * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--             Integral diameter of cloud drops
                nc_0 = nc(k,j,i) * dc
!
!--             Condensation needs only to be calculated in supersaturated regions
                IF ( sat > 0.0_wp )  THEN
!
!--                Condensation rate of cloud water content
!--                after KK scheme.
!--                (see: Khairoutdinov + Kogan, 2000, Mon. Wea. Rev.,128)
                   cond      = 2.0_wp * pi * nc_0 * g_fac * sat / hyrho(k)
                   cond_max  = q(k,j,i) - q_s - qc(k,j,i) - qr(k,j,i)
                   cond      = MIN( cond, cond_max / dt_micro )

                   qc(k,j,i) = qc(k,j,i) + cond * dt_micro
                ELSEIF ( sat < 0.0_wp ) THEN
                   evap      = 2.0_wp * pi * nc_0 * g_fac * sat / hyrho(k)
                   evap      = MAX( evap, -qc(k,j,i) / dt_micro )

                   qc(k,j,i) = qc(k,j,i) + evap * dt_micro
                ENDIF
                IF ( nc(k,j,i) * xcmin > qc(k,j,i) * hyrho(k) )  THEN
                   nc(k,j,i) = qc(k,j,i) * hyrho(k) / xcmin
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point_s(66), 'condensation', 'stop' )

    END SUBROUTINE condensation


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Autoconversion rate (Seifert and Beheng, 2006).
!------------------------------------------------------------------------------!
    SUBROUTINE autoconversion

       USE arrays_3d,                                                          &
           ONLY:  diss, dzu, nc, nr, qc, qr

       USE cloud_parameters,                                                   &
           ONLY:  hyrho

       USE control_parameters,                                                 &
           ONLY:  microphysics_morrison, rho_surface

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE grid_variables,                                                     &
           ONLY:  dx, dy

       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nyng, nysg, nzb, nzt, wall_flags_0

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  alpha_cc          !<
       REAL(wp)     ::  autocon           !<
       REAL(wp)     ::  dissipation       !<
       REAL(wp)     ::  flag              !< flag to mask topography grid points
       REAL(wp)     ::  k_au              !<
       REAL(wp)     ::  l_mix             !<
       REAL(wp)     ::  nc_auto           !<
       REAL(wp)     ::  nu_c              !<
       REAL(wp)     ::  phi_au            !<
       REAL(wp)     ::  r_cc              !<
       REAL(wp)     ::  rc                !<
       REAL(wp)     ::  re_lambda         !<
       REAL(wp)     ::  sigma_cc          !<
       REAL(wp)     ::  tau_cloud         !<
       REAL(wp)     ::  xc                !<

       CALL cpu_log( log_point_s(55), 'autoconversion', 'start' )

       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

                nc_auto = MERGE( nc(k,j,i), nc_const, microphysics_morrison )

                IF ( qc(k,j,i) > eps_sb  .AND.  nc_auto > eps_mr )  THEN

                   k_au = k_cc / ( 20.0_wp * x0 )
!
!--                Intern time scale of coagulation (Seifert and Beheng, 2006):
!--                (1.0_wp - qc(k,j,i) / ( qc(k,j,i) + qr(k,j,i) ))
                   tau_cloud = MAX( 1.0_wp - qc(k,j,i) / ( qr(k,j,i) +         &
                                    qc(k,j,i) ), 0.0_wp )
!
!--                Universal function for autoconversion process
!--                (Seifert and Beheng, 2006):
                   phi_au = 600.0_wp * tau_cloud**0.68_wp *                    &
                            ( 1.0_wp - tau_cloud**0.68_wp )**3
!
!--                Shape parameter of gamma distribution (Geoffroy et al., 2010):
!--                (Use constant nu_c = 1.0_wp instead?)
                   nu_c = 1.0_wp !MAX( 0.0_wp, 1580.0_wp * hyrho(k) * qc(k,j,i) - 0.28_wp )
!
!--                Mean weight of cloud droplets:
                   xc = MAX( hyrho(k) * qc(k,j,i) / nc_auto, xcmin) 
!
!--                Parameterized turbulence effects on autoconversion (Seifert,
!--                Nuijens and Stevens, 2010)
                   IF ( collision_turbulence )  THEN
!
!--                   Weight averaged radius of cloud droplets:
                      rc = 0.5_wp * ( xc * dpirho_l )**( 1.0_wp / 3.0_wp )

                      alpha_cc = ( a_1 + a_2 * nu_c ) / ( 1.0_wp + a_3 * nu_c )
                      r_cc     = ( b_1 + b_2 * nu_c ) / ( 1.0_wp + b_3 * nu_c )
                      sigma_cc = ( c_1 + c_2 * nu_c ) / ( 1.0_wp + c_3 * nu_c )
!
!--                   Mixing length (neglecting distance to ground and
!--                   stratification)
                      l_mix = ( dx * dy * dzu(k) )**( 1.0_wp / 3.0_wp )
!
!--                   Limit dissipation rate according to Seifert, Nuijens and
!--                   Stevens (2010)
                      dissipation = MIN( 0.06_wp, diss(k,j,i) )
!
!--                   Compute Taylor-microscale Reynolds number:
                      re_lambda = 6.0_wp / 11.0_wp *                           &
                                  ( l_mix / c_const )**( 2.0_wp / 3.0_wp ) *   &
                                  SQRT( 15.0_wp / kin_vis_air ) *              &
                                  dissipation**( 1.0_wp / 6.0_wp )
!
!--                   The factor of 1.0E4 is needed to convert the dissipation
!--                   rate from m2 s-3 to cm2 s-3.
                      k_au = k_au * ( 1.0_wp +                                 &
                                      dissipation * 1.0E4_wp *                 &
                                      ( re_lambda * 1.0E-3_wp )**0.25_wp *     &
                                      ( alpha_cc * EXP( -1.0_wp * ( ( rc -     &
                                                                      r_cc ) / &
                                                        sigma_cc )**2          &
                                                      ) + beta_cc              &
                                      )                                        &
                                    )
                   ENDIF
!
!--                Autoconversion rate (Seifert and Beheng, 2006):
                   autocon = k_au * ( nu_c + 2.0_wp ) * ( nu_c + 4.0_wp ) /    &
                             ( nu_c + 1.0_wp )**2 * qc(k,j,i)**2 * xc**2 *     &
                             ( 1.0_wp + phi_au / ( 1.0_wp - tau_cloud )**2 ) * &
                             rho_surface
                   autocon = MIN( autocon, qc(k,j,i) / dt_micro )

                   qr(k,j,i) = qr(k,j,i) + autocon * dt_micro * flag
                   qc(k,j,i) = qc(k,j,i) - autocon * dt_micro * flag
                   nr(k,j,i) = nr(k,j,i) + autocon / x0 * hyrho(k) * dt_micro  &
                                                              * flag
                   IF ( microphysics_morrison ) THEN
                      nc(k,j,i) = nc(k,j,i) - MIN( nc(k,j,i), 2.0_wp *         &
                                    autocon / x0 * hyrho(k) * dt_micro * flag )
                   ENDIF

                ENDIF

             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point_s(55), 'autoconversion', 'stop' )

    END SUBROUTINE autoconversion


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Autoconversion process (Kessler, 1969).
!------------------------------------------------------------------------------!
    SUBROUTINE autoconversion_kessler

       USE arrays_3d,                                                          &
           ONLY:  dzw, pt, prr, q, qc

       USE cloud_parameters,                                                   &
           ONLY:  l_d_cp, pt_d_t

       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nyng, nysg, nzb, nzt, wall_flags_0

       USE kinds

       USE surface_mod,                                                        &
           ONLY:  get_topography_top_index_ji


       IMPLICIT NONE

       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  k      !<
       INTEGER(iwp) ::  k_wall !< topgraphy top index

       REAL(wp)    ::  dqdt_precip !<
       REAL(wp)    ::  flag        !< flag to mask topography grid points

       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
!
!--          Determine vertical index of topography top
             k_wall = get_topography_top_index_ji( j, i, 's' )
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

                IF ( qc(k,j,i) > ql_crit )  THEN
                   dqdt_precip = prec_time_const * ( qc(k,j,i) - ql_crit )
                ELSE
                   dqdt_precip = 0.0_wp
                ENDIF

                qc(k,j,i) = qc(k,j,i) - dqdt_precip * dt_micro * flag
                q(k,j,i)  = q(k,j,i)  - dqdt_precip * dt_micro * flag
                pt(k,j,i) = pt(k,j,i) + dqdt_precip * dt_micro * l_d_cp *      &
                                        pt_d_t(k)              * flag

!
!--             Compute the rain rate (stored on surface grid point)
                prr(k_wall,j,i) = prr(k_wall,j,i) + dqdt_precip * dzw(k) * flag

             ENDDO
          ENDDO
       ENDDO

   END SUBROUTINE autoconversion_kessler


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Accretion rate (Seifert and Beheng, 2006).
!------------------------------------------------------------------------------!
    SUBROUTINE accretion

       USE arrays_3d,                                                          &
           ONLY:  diss, qc, qr, nc

       USE cloud_parameters,                                                   &
           ONLY:  hyrho

       USE control_parameters,                                                 &
           ONLY:  microphysics_morrison, rho_surface

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nyng, nysg, nzb, nzt, wall_flags_0

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  accr              !<
       REAL(wp)     ::  flag              !< flag to mask topography grid points
       REAL(wp)     ::  k_cr              !<
       REAL(wp)     ::  nc_accr           !<
       REAL(wp)     ::  phi_ac            !<
       REAL(wp)     ::  tau_cloud         !<
       REAL(wp)     ::  xc                !<


       CALL cpu_log( log_point_s(56), 'accretion', 'start' )

       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

                nc_accr = MERGE( nc(k,j,i), nc_const, microphysics_morrison )

                IF ( ( qc(k,j,i) > eps_sb )  .AND.  ( qr(k,j,i) > eps_sb )    &
                                             .AND.  ( nc_accr > eps_mr ) ) THEN
!
!--                Intern time scale of coagulation (Seifert and Beheng, 2006):
                   tau_cloud = 1.0_wp - qc(k,j,i) / ( qc(k,j,i) + qr(k,j,i) )
!
!--                Universal function for accretion process (Seifert and
!--                Beheng, 2001):
                   phi_ac = ( tau_cloud / ( tau_cloud + 5.0E-5_wp ) )**4

!
!--                Mean weight of cloud drops
                   xc = MAX( (hyrho(k) * qc(k,j,i) / nc_accr), xcmin)
!
!--                Parameterized turbulence effects on autoconversion (Seifert,
!--                Nuijens and Stevens, 2010). The factor of 1.0E4 is needed to
!--                convert the dissipation rate (diss) from m2 s-3 to cm2 s-3.
                   IF ( collision_turbulence )  THEN
                      k_cr = k_cr0 * ( 1.0_wp + 0.05_wp *                      &
                                       MIN( 600.0_wp,                          &
                                            diss(k,j,i) * 1.0E4_wp )**0.25_wp  &
                                     )
                   ELSE
                      k_cr = k_cr0
                   ENDIF
!
!--                Accretion rate (Seifert and Beheng, 2006):
                   accr = k_cr * qc(k,j,i) * qr(k,j,i) * phi_ac *              &
                          SQRT( rho_surface * hyrho(k) )
                   accr = MIN( accr, qc(k,j,i) / dt_micro )

                   qr(k,j,i) = qr(k,j,i) + accr * dt_micro * flag
                   qc(k,j,i) = qc(k,j,i) - accr * dt_micro * flag
                   IF ( microphysics_morrison )  THEN
                      nc(k,j,i) = nc(k,j,i) - MIN( nc(k,j,i),                  &
                                         accr / xc * hyrho(k) * dt_micro * flag)
                   ENDIF

                ENDIF

             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point_s(56), 'accretion', 'stop' )

    END SUBROUTINE accretion


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Collisional breakup rate (Seifert, 2008).
!------------------------------------------------------------------------------!
    SUBROUTINE selfcollection_breakup

       USE arrays_3d,                                                          &
           ONLY:  nr, qr

       USE cloud_parameters,                                                   &
           ONLY:  hyrho

       USE control_parameters,                                                 &
           ONLY:  rho_surface

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nyng, nysg, nzb, nzt, wall_flags_0

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  breakup           !<
       REAL(wp)     ::  dr                !<
       REAL(wp)     ::  flag              !< flag to mask topography grid points
       REAL(wp)     ::  phi_br            !<
       REAL(wp)     ::  selfcoll          !<

       CALL cpu_log( log_point_s(57), 'selfcollection', 'start' )

       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

                IF ( qr(k,j,i) > eps_sb )  THEN
!
!--                Selfcollection rate (Seifert and Beheng, 2001):
                   selfcoll = k_rr * nr(k,j,i) * qr(k,j,i) *                   &
                              SQRT( hyrho(k) * rho_surface )
!
!--                Weight averaged diameter of rain drops:
                   dr = ( hyrho(k) * qr(k,j,i) /                               &
                          nr(k,j,i) * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--                Collisional breakup rate (Seifert, 2008):
                   IF ( dr >= 0.3E-3_wp )  THEN
                      phi_br  = k_br * ( dr - 1.1E-3_wp )
                      breakup = selfcoll * ( phi_br + 1.0_wp )
                   ELSE
                      breakup = 0.0_wp
                   ENDIF

                   selfcoll = MAX( breakup - selfcoll, -nr(k,j,i) / dt_micro )
                   nr(k,j,i) = nr(k,j,i) + selfcoll * dt_micro * flag

                ENDIF
             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point_s(57), 'selfcollection', 'stop' )

    END SUBROUTINE selfcollection_breakup


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Evaporation of precipitable water. Condensation is neglected for
!> precipitable water.
!------------------------------------------------------------------------------!
    SUBROUTINE evaporation_rain

       USE arrays_3d,                                                          &
           ONLY:  hyp, nr, pt, q,  qc, qr

       USE cloud_parameters,                                                   &
           ONLY:  hyrho, l_d_cp, l_d_r, l_v, r_v, t_d_pt

       USE constants,                                                          &
           ONLY:  pi

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE diagnostic_quantities_mod,                                          &
           ONLY: e_s, magnus, q_s, sat, supersaturation, t_l

       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nyng, nysg, nzb, nzt, wall_flags_0

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  dr                !<
       REAL(wp)     ::  evap              !<
       REAL(wp)     ::  evap_nr           !<
       REAL(wp)     ::  f_vent            !<
       REAL(wp)     ::  flag              !< flag to mask topography grid points
       REAL(wp)     ::  g_evap            !<
       REAL(wp)     ::  lambda_r          !<
       REAL(wp)     ::  mu_r              !<
       REAL(wp)     ::  mu_r_2            !<
       REAL(wp)     ::  mu_r_5d2          !<
       REAL(wp)     ::  nr_0              !<
       REAL(wp)     ::  temp              !<
       REAL(wp)     ::  xr                !<

       CALL cpu_log( log_point_s(58), 'evaporation', 'start' )

       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

                IF ( qr(k,j,i) > eps_sb )  THEN

!
!--                Call calculation of supersaturation located 
!--                in diagnostic_quantities_mod
                   CALL supersaturation ( i, j, k )
!
!--                Evaporation needs only to be calculated in subsaturated regions
                   IF ( sat < 0.0_wp )  THEN
!
!--                   Actual temperature:
                      temp = t_l + l_d_cp * ( qc(k,j,i) + qr(k,j,i) )

                      g_evap = 1.0_wp / ( ( l_v / ( r_v * temp ) - 1.0_wp ) *     &
                                          l_v / ( thermal_conductivity_l * temp ) &
                                          + r_v * temp / ( diff_coeff_l * e_s )   &
                                        )
!
!--                   Mean weight of rain drops
                      xr = hyrho(k) * qr(k,j,i) / nr(k,j,i)
!
!--                   Weight averaged diameter of rain drops:
                      dr = ( xr * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--                   Compute ventilation factor and intercept parameter
!--                   (Seifert and Beheng, 2006; Seifert, 2008):
                      IF ( ventilation_effect )  THEN
!
!--                      Shape parameter of gamma distribution (Milbrandt and Yau,
!--                      2005; Stevens and Seifert, 2008):
                         mu_r = 10.0_wp * ( 1.0_wp + TANH( 1.2E3_wp *          &
                                                          ( dr - 1.4E-3_wp ) ) )
!
!--                      Slope parameter of gamma distribution (Seifert, 2008):
                         lambda_r = ( ( mu_r + 3.0_wp ) * ( mu_r + 2.0_wp ) *  &
                                      ( mu_r + 1.0_wp )                        &
                                    )**( 1.0_wp / 3.0_wp ) / dr

                         mu_r_2   = mu_r + 2.0_wp
                         mu_r_5d2 = mu_r + 2.5_wp

                         f_vent = a_vent * gamm( mu_r_2 ) *                     &
                                  lambda_r**( -mu_r_2 ) + b_vent *              &
                                  schmidt_p_1d3 * SQRT( a_term / kin_vis_air ) *&
                                  gamm( mu_r_5d2 ) * lambda_r**( -mu_r_5d2 ) *  &
                                  ( 1.0_wp -                                    &
                                    0.5_wp * ( b_term / a_term ) *              &
                                    ( lambda_r / ( c_term + lambda_r )          &
                                    )**mu_r_5d2 -                               &
                                    0.125_wp * ( b_term / a_term )**2 *         &
                                    ( lambda_r / ( 2.0_wp * c_term + lambda_r ) &
                                    )**mu_r_5d2 -                               &
                                    0.0625_wp * ( b_term / a_term )**3 *        &
                                    ( lambda_r / ( 3.0_wp * c_term + lambda_r ) &
                                    )**mu_r_5d2 -                               &
                                    0.0390625_wp * ( b_term / a_term )**4 *     &
                                    ( lambda_r / ( 4.0_wp * c_term + lambda_r ) &
                                    )**mu_r_5d2                                 &
                                  )

                         nr_0   = nr(k,j,i) * lambda_r**( mu_r + 1.0_wp ) /    &
                                  gamm( mu_r + 1.0_wp )
                      ELSE
                         f_vent = 1.0_wp
                         nr_0   = nr(k,j,i) * dr
                      ENDIF
   !
   !--                Evaporation rate of rain water content (Seifert and
   !--                Beheng, 2006):
                      evap    = 2.0_wp * pi * nr_0 * g_evap * f_vent * sat /   &
                                hyrho(k)
                      evap    = MAX( evap, -qr(k,j,i) / dt_micro )
                      evap_nr = MAX( c_evap * evap / xr * hyrho(k),            &
                                     -nr(k,j,i) / dt_micro )

                      qr(k,j,i) = qr(k,j,i) + evap    * dt_micro * flag
                      nr(k,j,i) = nr(k,j,i) + evap_nr * dt_micro * flag

                   ENDIF
                ENDIF

             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point_s(58), 'evaporation', 'stop' )

    END SUBROUTINE evaporation_rain


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sedimentation of cloud droplets (Ackermann et al., 2009, MWR).
!------------------------------------------------------------------------------!
    SUBROUTINE sedimentation_cloud

       USE arrays_3d,                                                          &
           ONLY:  ddzu, dzu, nc, pt, prr, q, qc

       USE cloud_parameters,                                                   &
           ONLY:  hyrho, l_d_cp, pt_d_t

       USE control_parameters,                                                 &
           ONLY:  call_microphysics_at_all_substeps,                           &
                  intermediate_timestep_count, microphysics_morrison

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nyng, nysg, nzb, nzt, wall_flags_0

       USE kinds

       USE statistics,                                                         &
           ONLY:  weight_substep


       IMPLICIT NONE

       INTEGER(iwp) ::  i             !<
       INTEGER(iwp) ::  j             !<
       INTEGER(iwp) ::  k             !<

       REAL(wp) ::  flag              !< flag to mask topography grid points
       REAL(wp) ::  nc_sedi           !<

       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_qc !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_nc !<


       CALL cpu_log( log_point_s(59), 'sed_cloud', 'start' )

       sed_qc(nzt+1) = 0.0_wp
       sed_nc(nzt+1) = 0.0_wp

       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzt, nzb+1, -1
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )
                nc_sedi = MERGE ( nc(k,j,i), nc_const, microphysics_morrison )

!
!--             Sedimentation fluxes for number concentration are only calculated
!--             for cloud_scheme = 'morrison'
                IF ( microphysics_morrison ) THEN
                   IF ( qc(k,j,i) > eps_sb  .AND.  nc(k,j,i) > eps_mr )  THEN
                      sed_nc(k) = sed_qc_const *                               &
                              ( qc(k,j,i) * hyrho(k) )**( 2.0_wp / 3.0_wp ) *  &
                              ( nc(k,j,i) )**( 1.0_wp / 3.0_wp )
                   ELSE
                      sed_nc(k) = 0.0_wp
                   ENDIF

                   sed_nc(k) = MIN( sed_nc(k), hyrho(k) * dzu(k+1) *           &
                                    nc(k,j,i) / dt_micro + sed_nc(k+1)         &
                                  ) * flag

                   nc(k,j,i) = nc(k,j,i) + ( sed_nc(k+1) - sed_nc(k) ) *       &
                                         ddzu(k+1) / hyrho(k) * dt_micro * flag
                ENDIF

                IF ( qc(k,j,i) > eps_sb .AND.  nc_sedi > eps_mr )  THEN
                   sed_qc(k) = sed_qc_const * nc_sedi**( -2.0_wp / 3.0_wp ) *  &
                               ( qc(k,j,i) * hyrho(k) )**( 5.0_wp / 3.0_wp ) * &
                                                                           flag
                ELSE
                   sed_qc(k) = 0.0_wp
                ENDIF

                sed_qc(k) = MIN( sed_qc(k), hyrho(k) * dzu(k+1) * q(k,j,i) /   &
                                            dt_micro + sed_qc(k+1)             &
                               ) * flag

                q(k,j,i)  = q(k,j,i)  + ( sed_qc(k+1) - sed_qc(k) ) *          &
                                        ddzu(k+1) / hyrho(k) * dt_micro * flag
                qc(k,j,i) = qc(k,j,i) + ( sed_qc(k+1) - sed_qc(k) ) *          &
                                        ddzu(k+1) / hyrho(k) * dt_micro * flag
                pt(k,j,i) = pt(k,j,i) - ( sed_qc(k+1) - sed_qc(k) ) *          &
                                        ddzu(k+1) / hyrho(k) * l_d_cp *        &
                                        pt_d_t(k) * dt_micro            * flag

!
!--             Compute the precipitation rate due to cloud (fog) droplets
                IF ( call_microphysics_at_all_substeps )  THEN
                   prr(k,j,i) = prr(k,j,i) +  sed_qc(k) / hyrho(k)             &
                                * weight_substep(intermediate_timestep_count)  &
                                * flag
                ELSE
                   prr(k,j,i) = prr(k,j,i) + sed_qc(k) / hyrho(k) * flag
                ENDIF

             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point_s(59), 'sed_cloud', 'stop' )

    END SUBROUTINE sedimentation_cloud


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computation of sedimentation flux. Implementation according to Stevens
!> and Seifert (2008). Code is based on UCLA-LES.
!------------------------------------------------------------------------------!
    SUBROUTINE sedimentation_rain

       USE arrays_3d,                                                          &
           ONLY:  ddzu, dzu, nr, pt, prr, q, qr

       USE cloud_parameters,                                                   &
           ONLY:  hyrho, l_d_cp, pt_d_t

       USE control_parameters,                                                 &
           ONLY:  call_microphysics_at_all_substeps, intermediate_timestep_count
       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nyng, nysg, nzb, nzt, wall_flags_0

       USE kinds

       USE statistics,                                                         &
           ONLY:  weight_substep

       USE surface_mod,                                                        &
           ONLY :  bc_h

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  k_run         !<
       INTEGER(iwp) ::  l             !< running index of surface type
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp)     ::  c_run                      !<
       REAL(wp)     ::  d_max                      !<
       REAL(wp)     ::  d_mean                     !<
       REAL(wp)     ::  d_min                      !<
       REAL(wp)     ::  dr                         !<
       REAL(wp)     ::  flux                       !<
       REAL(wp)     ::  flag                       !< flag to mask topography grid points
       REAL(wp)     ::  lambda_r                   !<
       REAL(wp)     ::  mu_r                       !<
       REAL(wp)     ::  z_run                      !<

       REAL(wp), DIMENSION(nzb:nzt+1) ::  c_nr     !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  c_qr     !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  nr_slope !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  qr_slope !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_nr   !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_qr   !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  w_nr     !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  w_qr     !<

       CALL cpu_log( log_point_s(60), 'sed_rain', 'start' )

!
!--    Compute velocities
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

                IF ( qr(k,j,i) > eps_sb )  THEN
!
!--                Weight averaged diameter of rain drops:
                   dr = ( hyrho(k) * qr(k,j,i) /                               &
                          nr(k,j,i) * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--                Shape parameter of gamma distribution (Milbrandt and Yau, 2005;
!--                Stevens and Seifert, 2008):
                   mu_r = 10.0_wp * ( 1.0_wp + TANH( 1.2E3_wp *                &
                                                     ( dr - 1.4E-3_wp ) ) )
!
!--                Slope parameter of gamma distribution (Seifert, 2008):
                   lambda_r = ( ( mu_r + 3.0_wp ) * ( mu_r + 2.0_wp ) *        &
                                ( mu_r + 1.0_wp ) )**( 1.0_wp / 3.0_wp ) / dr

                   w_nr(k) = MAX( 0.1_wp, MIN( 20.0_wp,                        &
                                               a_term - b_term * ( 1.0_wp +    &
                                                  c_term /                     &
                                                  lambda_r )**( -1.0_wp *      &
                                                  ( mu_r + 1.0_wp ) )          &
                                              )                                &
                                ) * flag

                   w_qr(k) = MAX( 0.1_wp, MIN( 20.0_wp,                        &
                                               a_term - b_term * ( 1.0_wp +    &
                                                  c_term /                     &
                                                  lambda_r )**( -1.0_wp *      &
                                                  ( mu_r + 4.0_wp ) )          &
                                             )                                 &
                                ) * flag
                ELSE
                   w_nr(k) = 0.0_wp
                   w_qr(k) = 0.0_wp
                ENDIF
             ENDDO
!
!--          Adjust boundary values using surface data type.
!--          Upward-facing
             surf_s = bc_h(0)%start_index(j,i)
             surf_e = bc_h(0)%end_index(j,i)
             DO  m = surf_s, surf_e
                k         = bc_h(0)%k(m)
                w_nr(k-1) = w_nr(k)
                w_qr(k-1) = w_qr(k)
             ENDDO
!
!--          Downward-facing
             surf_s = bc_h(1)%start_index(j,i)
             surf_e = bc_h(1)%end_index(j,i)
             DO  m = surf_s, surf_e
                k         = bc_h(1)%k(m)
                w_nr(k+1) = w_nr(k)
                w_qr(k+1) = w_qr(k)
             ENDDO
!
!--          Model top boundary value
             w_nr(nzt+1) = 0.0_wp
             w_qr(nzt+1) = 0.0_wp
!
!--          Compute Courant number
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

                c_nr(k) = 0.25_wp * ( w_nr(k-1) +                              &
                                      2.0_wp * w_nr(k) + w_nr(k+1) ) *         &
                          dt_micro * ddzu(k) * flag
                c_qr(k) = 0.25_wp * ( w_qr(k-1) +                              &
                                      2.0_wp * w_qr(k) + w_qr(k+1) ) *         &
                          dt_micro * ddzu(k) * flag
             ENDDO
!
!--          Limit slopes with monotonized centered (MC) limiter (van Leer, 1977):
             IF ( limiter_sedimentation )  THEN

                DO k = nzb+1, nzt
!
!--                Predetermine flag to mask topography
                   flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

                   d_mean = 0.5_wp * ( qr(k+1,j,i) - qr(k-1,j,i) )
                   d_min  = qr(k,j,i) - MIN( qr(k+1,j,i), qr(k,j,i), qr(k-1,j,i) )
                   d_max  = MAX( qr(k+1,j,i), qr(k,j,i), qr(k-1,j,i) ) - qr(k,j,i)

                   qr_slope(k) = SIGN(1.0_wp, d_mean) * MIN ( 2.0_wp * d_min,  &
                                                              2.0_wp * d_max,  &
                                                              ABS( d_mean ) )  &
                                                      * flag

                   d_mean = 0.5_wp * ( nr(k+1,j,i) - nr(k-1,j,i) )
                   d_min  = nr(k,j,i) - MIN( nr(k+1,j,i), nr(k,j,i), nr(k-1,j,i) )
                   d_max  = MAX( nr(k+1,j,i), nr(k,j,i), nr(k-1,j,i) ) - nr(k,j,i)

                   nr_slope(k) = SIGN(1.0_wp, d_mean) * MIN ( 2.0_wp * d_min,  &
                                                              2.0_wp * d_max,  &
                                                              ABS( d_mean ) )
                ENDDO

             ELSE

                nr_slope = 0.0_wp
                qr_slope = 0.0_wp

             ENDIF

             sed_nr(nzt+1) = 0.0_wp
             sed_qr(nzt+1) = 0.0_wp
!
!--          Compute sedimentation flux
             DO  k = nzt, nzb+1, -1
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )
!
!--             Sum up all rain drop number densities which contribute to the flux
!--             through k-1/2
                flux  = 0.0_wp
                z_run = 0.0_wp ! height above z(k)
                k_run = k
                c_run = MIN( 1.0_wp, c_nr(k) )
                DO WHILE ( c_run > 0.0_wp  .AND.  k_run <= nzt )
                   flux  = flux + hyrho(k_run) *                               &
                           ( nr(k_run,j,i) + nr_slope(k_run) *                 &
                           ( 1.0_wp - c_run ) * 0.5_wp ) * c_run * dzu(k_run)  &
                                              * flag
                   z_run = z_run + dzu(k_run) * flag
                   k_run = k_run + 1          * flag
                   c_run = MIN( 1.0_wp, c_nr(k_run) - z_run * ddzu(k_run) )    &
                                              * flag
                ENDDO
!
!--             It is not allowed to sediment more rain drop number density than
!--             available
                flux = MIN( flux,                                              &
                            hyrho(k) * dzu(k+1) * nr(k,j,i) + sed_nr(k+1) *    &
                            dt_micro                                           &
                          )

                sed_nr(k) = flux / dt_micro * flag
                nr(k,j,i) = nr(k,j,i) + ( sed_nr(k+1) - sed_nr(k) ) *          &
                                        ddzu(k+1) / hyrho(k) * dt_micro * flag
!
!--             Sum up all rain water content which contributes to the flux
!--             through k-1/2
                flux  = 0.0_wp
                z_run = 0.0_wp ! height above z(k)
                k_run = k
                c_run = MIN( 1.0_wp, c_qr(k) )

                DO WHILE ( c_run > 0.0_wp  .AND.  k_run <= nzt )

                   flux  = flux + hyrho(k_run) * ( qr(k_run,j,i) +             &
                                  qr_slope(k_run) * ( 1.0_wp - c_run ) *       &
                                  0.5_wp ) * c_run * dzu(k_run) * flag
                   z_run = z_run + dzu(k_run)                   * flag
                   k_run = k_run + 1                            * flag
                   c_run = MIN( 1.0_wp, c_qr(k_run) - z_run * ddzu(k_run) )    &
                                                                * flag

                ENDDO
!
!--             It is not allowed to sediment more rain water content than
!--             available
                flux = MIN( flux,                                              &
                            hyrho(k) * dzu(k) * qr(k,j,i) + sed_qr(k+1) *      &
                            dt_micro                                           &
                          )

                sed_qr(k) = flux / dt_micro * flag

                qr(k,j,i) = qr(k,j,i) + ( sed_qr(k+1) - sed_qr(k) ) *          &
                                        ddzu(k+1) / hyrho(k) * dt_micro * flag
                q(k,j,i)  = q(k,j,i)  + ( sed_qr(k+1) - sed_qr(k) ) *          &
                                        ddzu(k+1) / hyrho(k) * dt_micro * flag
                pt(k,j,i) = pt(k,j,i) - ( sed_qr(k+1) - sed_qr(k) ) *          &
                                        ddzu(k+1) / hyrho(k) * l_d_cp *        &
                                        pt_d_t(k) * dt_micro            * flag
!
!--             Compute the rain rate
                IF ( call_microphysics_at_all_substeps )  THEN
                   prr(k,j,i) = prr(k,j,i) + sed_qr(k) / hyrho(k)             &
                                * weight_substep(intermediate_timestep_count) &
                                * flag
                ELSE
                   prr(k,j,i) = prr(k,j,i) + sed_qr(k) / hyrho(k) * flag
                ENDIF

             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point_s(60), 'sed_rain', 'stop' )

    END SUBROUTINE sedimentation_rain


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computation of the precipitation amount due to gravitational settling of
!> rain and cloud (fog) droplets
!------------------------------------------------------------------------------!
    SUBROUTINE calc_precipitation_amount

       USE arrays_3d,                                                          &
           ONLY:  precipitation_amount, prr

       USE cloud_parameters,                                                   &
           ONLY:  hyrho

       USE control_parameters,                                                 &
           ONLY:  call_microphysics_at_all_substeps, dt_do2d_xy, dt_3d,        &
                  intermediate_timestep_count, intermediate_timestep_count_max,&
                  precipitation_amount_interval, time_do2d_xy

       USE indices,                                                            &
           ONLY:  nxl, nxr, nys, nyn, nzb, nzt, wall_flags_0

       USE kinds

       USE surface_mod,                                                        &
           ONLY :  bc_h

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index y direction
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       IF ( ( dt_do2d_xy - time_do2d_xy ) < precipitation_amount_interval .AND.&
            ( .NOT. call_microphysics_at_all_substeps .OR.                     &
            intermediate_timestep_count == intermediate_timestep_count_max ) ) &
       THEN
!
!--       Run over all upward-facing surface elements, i.e. non-natural,
!--       natural and urban
          DO  m = 1, bc_h(0)%ns
             i = bc_h(0)%i(m)
             j = bc_h(0)%j(m)
             k = bc_h(0)%k(m)
             precipitation_amount(j,i) = precipitation_amount(j,i) +           &
                                               prr(k,j,i) * hyrho(k) * dt_3d
          ENDDO

       ENDIF

    END SUBROUTINE calc_precipitation_amount


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Control of microphysics for grid points i,j
!------------------------------------------------------------------------------!

    SUBROUTINE microphysics_control_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  hyp, nc, nr, pt, pt_init, prr, q, qc, qr, zu

       USE cloud_parameters,                                                   &
           ONLY:  cp, hyrho, pt_d_t, r_d, t_d_pt

       USE control_parameters,                                                 &
           ONLY:  call_microphysics_at_all_substeps, dt_3d, g,                 &
                  intermediate_timestep_count, large_scale_forcing,            &
                  lsf_surf, microphysics_morrison, microphysics_seifert,       &
                  microphysics_kessler, pt_surface, rho_surface,               &
                  surface_pressure

       USE indices,                                                            &
           ONLY:  nzb, nzt

       USE kinds

       USE statistics,                                                         &
           ONLY:  weight_pres

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  t_surface         !<

       IF ( large_scale_forcing  .AND.  lsf_surf ) THEN
!
!--       Calculate:
!--       pt / t : ratio of potential and actual temperature (pt_d_t)
!--       t / pt : ratio of actual and potential temperature (t_d_pt)
!--       p_0(z) : vertical profile of the hydrostatic pressure (hyp)
          t_surface = pt_surface * ( surface_pressure / 1000.0_wp )**0.286_wp
          DO  k = nzb, nzt+1
             hyp(k)    = surface_pressure * 100.0_wp * &
                         ( ( t_surface - g / cp * zu(k) ) / t_surface )**(1.0_wp / 0.286_wp)
             pt_d_t(k) = ( 100000.0_wp / hyp(k) )**0.286_wp
             t_d_pt(k) = 1.0_wp / pt_d_t(k)
             hyrho(k)  = hyp(k) / ( r_d * t_d_pt(k) * pt_init(k) )
          ENDDO
!
!--       Compute reference density
          rho_surface = surface_pressure * 100.0_wp / ( r_d * t_surface )
       ENDIF

!
!--    Compute length of time step
       IF ( call_microphysics_at_all_substeps )  THEN
          dt_micro = dt_3d * weight_pres(intermediate_timestep_count)
       ELSE
          dt_micro = dt_3d
       ENDIF

!
!--    Use 1d arrays
       q_1d(:)  = q(:,j,i)
       pt_1d(:) = pt(:,j,i)
       qc_1d(:) = qc(:,j,i)
       IF ( microphysics_seifert )  THEN
          qr_1d(:) = qr(:,j,i)
          nr_1d(:) = nr(:,j,i)
       ENDIF
       IF ( microphysics_morrison )  THEN
          nc_1d(:) = nc(:,j,i)
       ENDIF


!
!--    Reset precipitation rate
       IF ( intermediate_timestep_count == 1 )  prr(:,j,i) = 0.0_wp

!
!--    Compute cloud physics
       IF( microphysics_kessler )  THEN

          CALL autoconversion_kessler( i,j )
          IF ( cloud_water_sedimentation )  CALL sedimentation_cloud( i,j )

       ELSEIF ( microphysics_seifert )  THEN

          CALL adjust_cloud( i,j )
          IF ( microphysics_morrison ) CALL activation( i,j )
          IF ( microphysics_morrison ) CALL condensation( i,j )
          CALL autoconversion( i,j )
          CALL accretion( i,j )
          CALL selfcollection_breakup( i,j )
          CALL evaporation_rain( i,j )
          CALL sedimentation_rain( i,j )
          IF ( cloud_water_sedimentation )  CALL sedimentation_cloud( i,j )

       ENDIF

       CALL calc_precipitation_amount( i,j )

!
!--    Store results on the 3d arrays
       q(:,j,i)  = q_1d(:)
       pt(:,j,i) = pt_1d(:)
       IF ( microphysics_morrison ) THEN
          qc(:,j,i) = qc_1d(:)
          nc(:,j,i) = nc_1d(:)
       ENDIF
       IF ( microphysics_seifert )  THEN
          qr(:,j,i) = qr_1d(:)
          nr(:,j,i) = nr_1d(:)
       ENDIF

    END SUBROUTINE microphysics_control_ij

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Adjust number of raindrops to avoid nonlinear effects in
!> sedimentation and evaporation of rain drops due to too small or
!> too big weights of rain drops (Stevens and Seifert, 2008).
!> The same procedure is applied to cloud droplets if they are determined
!> prognostically. Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE adjust_cloud_ij( i, j )

       USE cloud_parameters,                                                   &
           ONLY:  hyrho

       USE control_parameters,                                                 &
           ONLY: microphysics_morrison

       USE indices,                                                            &
           ONLY:  nzb, nzt, wall_flags_0

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp) ::  flag                  !< flag to indicate first grid level above surface

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

          IF ( qr_1d(k) <= eps_sb )  THEN
             qr_1d(k) = 0.0_wp
             nr_1d(k) = 0.0_wp
          ELSE
!
!--          Adjust number of raindrops to avoid nonlinear effects in
!--          sedimentation and evaporation of rain drops due to too small or
!--          too big weights of rain drops (Stevens and Seifert, 2008).
             IF ( nr_1d(k) * xrmin > qr_1d(k) * hyrho(k) )  THEN
                nr_1d(k) = qr_1d(k) * hyrho(k) / xrmin * flag
             ELSEIF ( nr_1d(k) * xrmax < qr_1d(k) * hyrho(k) )  THEN
                nr_1d(k) = qr_1d(k) * hyrho(k) / xrmax * flag
             ENDIF

          ENDIF

          IF ( microphysics_morrison ) THEN
             IF ( qc_1d(k) <= eps_sb )  THEN
                qc_1d(k) = 0.0_wp
                nc_1d(k) = 0.0_wp
             ELSE
                IF ( nc_1d(k) * xcmin > qc_1d(k) * hyrho(k) )  THEN
                   nc_1d(k) = qc_1d(k) * hyrho(k) / xamin * flag
                ENDIF
             ENDIF
          ENDIF

       ENDDO

    END SUBROUTINE adjust_cloud_ij

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate number of activated condensation nucleii after simple activation
!> scheme of Twomey, 1959.
!------------------------------------------------------------------------------!
    SUBROUTINE activation_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  hyp, nr, pt, q,  qc, qr, nc

       USE cloud_parameters,                                                   &
           ONLY:  hyrho, l_d_cp, l_d_r, l_v, molecular_weight_of_solute,       &
                  molecular_weight_of_water, rho_l, rho_s, r_v, t_d_pt,        &
                  vanthoff

       USE constants,                                                          &
           ONLY:  pi

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nysg, nyng, nzb, nzt

       USE kinds

       USE control_parameters,                                                 &
          ONLY: simulated_time

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  activ             !<
       REAL(wp)     ::  afactor           !<
       REAL(wp)     ::  alpha             !<
       REAL(wp)     ::  beta_act          !<
       REAL(wp)     ::  bfactor           !<
       REAL(wp)     ::  e_s               !<
       REAL(wp)     ::  k_act             !<
       REAL(wp)     ::  n_act             !<
       REAL(wp)     ::  n_ccn             !<
       REAL(wp)     ::  q_s               !<
       REAL(wp)     ::  s_0               !<
       REAL(wp)     ::  sat               !<
       REAL(wp)     ::  sat_max           !<
       REAL(wp)     ::  sigma             !<
       REAL(wp)     ::  sigma_act         !<
       REAL(wp)     ::  t_int             !<
       REAL(wp)     ::  t_l               !<

       DO  k = nzb+1, nzt
!
!--       Actual liquid water temperature:
          t_l = t_d_pt(k) * pt_1d(k)

!
!--       Calculate actual temperature
          t_int = pt_1d(k) * ( hyp(k) / 100000.0_wp )**0.286_wp
!
!--             Saturation vapor pressure at t_l:
          e_s = 610.78_wp * EXP( 17.269_wp * ( t_l - 273.16_wp ) /             &
                                 ( t_l - 35.86_wp )                            &
                                 )
!
!--       Computation of saturation mixing ratio:
          q_s   = 0.622_wp * e_s / ( hyp(k) -  e_s )
          alpha = 0.622_wp * l_d_r * l_d_cp / ( t_l * t_l )
          q_s   = q_s * ( 1.0_wp + alpha * q_1d(k) ) /                         &
                        ( 1.0_wp + alpha * q_s )

!--       Supersaturation:
          sat   = ( q_1d(k) - qr_1d(k) - qc_1d(k) ) / q_s - 1.0_wp

!
!--       Prescribe parameters for activation
!--       (see: Bott + Trautmann, 2002, Atm. Res., 64)
          k_act  = 0.7_wp
          activ  = 0.0_wp

          IF ( sat >= 0.0 .AND. .NOT. curvature_solution_effects_bulk )  THEN
!
!--          Compute the number of activated Aerosols
!--          (see: Twomey, 1959, Pure and applied Geophysics, 43)
             n_act     = na_init * sat**k_act
!
!--          Compute the number of cloud droplets
!--          (see: Morrison + Grabowski, 2007, JAS, 64)
!            activ = MAX( n_act - nc_d1(k), 0.0_wp) / dt_micro

!
!--          Compute activation rate after Khairoutdinov and Kogan
!--          (see: Khairoutdinov + Kogan, 2000, Mon. Wea. Rev., 128)
             sat_max = 0.8_wp / 100.0_wp
             activ   = MAX( 0.0_wp, ( (na_init + nc_1d(k) ) * MIN              &
                 ( 1.0_wp, ( sat / sat_max )**k_act) - nc_1d(k) ) ) /          &
                  dt_micro

             nc_1d(k) = MIN( (nc_1d(k) + activ * dt_micro), na_init)
          ELSEIF ( sat >= 0.0 .AND. curvature_solution_effects_bulk )  THEN
!
!--          Curvature effect (afactor) with surface tension
!--          parameterization by Straka (2009)
             sigma = 0.0761_wp - 0.000155_wp * ( t_int - 273.15_wp )
             afactor = 2.0_wp * sigma / ( rho_l * r_v * t_int )
!
!--          Solute effect (bfactor)
             bfactor = vanthoff * molecular_weight_of_water *                  &
                  rho_s / ( molecular_weight_of_solute * rho_l )

!
!--          Prescribe power index that describes the soluble fraction
!--          of an aerosol particle (beta).
!--          (see: Morrison + Grabowski, 2007, JAS, 64)
             beta_act  = 0.5_wp
             sigma_act = sigma_bulk**( 1.0_wp + beta_act )
!
!--          Calculate mean geometric supersaturation (s_0) with
!--          parameterization by Khvorostyanov and Curry (2006)
             s_0   = dry_aerosol_radius **(- ( 1.0_wp + beta_act ) ) *         &
               ( 4.0_wp * afactor**3 / ( 27.0_wp * bfactor ) )**0.5_wp

!
!--          Calculate number of activated CCN as a function of
!--          supersaturation and taking Koehler theory into account
!--          (see: Khvorostyanov + Curry, 2006, J. Geo. Res., 111)
             n_ccn = ( na_init / 2.0_wp ) * ( 1.0_wp - ERF(                    &
                LOG( s_0 / sat ) / ( SQRT(2.0_wp) * LOG(sigma_act) ) ) )
             activ = MAX( ( n_ccn ) / dt_micro, 0.0_wp )

             nc_1d(k) = MIN( (nc_1d(k) + activ * dt_micro), na_init)
          ENDIF

       ENDDO

    END SUBROUTINE activation_ij

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate condensation rate for cloud water content (after Khairoutdinov and
!> Kogan, 2000).
!------------------------------------------------------------------------------!
    SUBROUTINE condensation_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  hyp, nr, pt, q,  qc, qr, nc

       USE cloud_parameters,                                                   &
           ONLY:  hyrho, l_d_cp, l_d_r, l_v, r_v, t_d_pt

       USE constants,                                                          &
           ONLY:  pi

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nysg, nyng, nzb, nzt

       USE kinds

       USE control_parameters,                                                 &
          ONLY: simulated_time

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  alpha             !<
       REAL(wp)     ::  cond              !<
       REAL(wp)     ::  cond_max          !<
       REAL(wp)     ::  dc                !<
       REAL(wp)     ::  e_s               !<
       REAL(wp)     ::  evap              !<
       REAL(wp)     ::  evap_nc           !<
       REAL(wp)     ::  g_fac             !<
       REAL(wp)     ::  nc_0              !<
       REAL(wp)     ::  q_s               !<
       REAL(wp)     ::  sat               !<
       REAL(wp)     ::  t_l               !<
       REAL(wp)     ::  temp              !<
       REAL(wp)     ::  xc                !<


       DO  k = nzb+1, nzt
!
!--       Actual liquid water temperature:
          t_l = t_d_pt(k) * pt_1d(k)
!
!--             Saturation vapor pressure at t_l:
          e_s = 610.78_wp * EXP( 17.269_wp * ( t_l - 273.16_wp ) /    &
                                 ( t_l - 35.86_wp )                   &
                                 )
!
!--       Computation of saturation mixing ratio:
          q_s   = 0.622_wp * e_s / ( hyp(k) - e_s )
          alpha = 0.622_wp * l_d_r * l_d_cp / ( t_l * t_l )
          q_s   = q_s * ( 1.0_wp + alpha * q_1d(k) ) /               &
                        ( 1.0_wp + alpha * q_s )

!--       Supersaturation:
          sat   = ( q_1d(k) - qr_1d(k) - qc_1d(k) ) / q_s - 1.0_wp


!
!--       Actual temperature:
          temp = t_l + l_d_cp * ( qc_1d(k) + qr_1d(k) )

          g_fac  = 1.0_wp / ( ( l_v / ( r_v * temp ) - 1.0_wp ) *        &
                              l_v / ( thermal_conductivity_l * temp )    &
                              + r_v * temp / ( diff_coeff_l * e_s )      &
                            )
!
!--       Mean weight of cloud drops
          IF ( nc_1d(k) <= 0.0_wp) CYCLE
          xc = MAX( (hyrho(k) * qc_1d(k) / nc_1d(k)), xcmin)
!
!--       Weight averaged diameter of cloud drops:
          dc   = ( xc * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--       Integral diameter of cloud drops
          nc_0 = nc_1d(k) * dc
!
!--       Condensation needs only to be calculated in supersaturated regions
          IF ( sat > 0.0_wp )  THEN
!
!--          Condensation rate of cloud water content
!--          after KK scheme.
!--          (see: Khairoutdinov + Kogan, 2000, Mon. Wea. Rev.,128)
             cond      = 2.0_wp * pi * nc_0 * g_fac * sat / hyrho(k)
             cond_max  = q_1d(k) - q_s - qc_1d(k) - qr_1d(k)
             cond      = MIN( cond, cond_max / dt_micro )

             qc_1d(k) = qc_1d(k) + cond * dt_micro
          ELSEIF ( sat < 0.0_wp ) THEN
             evap      = 2.0_wp * pi * nc_0 * g_fac * sat / hyrho(k)
             evap      = MAX( evap, -qc_1d(k) / dt_micro )

             qc_1d(k) = qc_1d(k) + evap * dt_micro
          ENDIF
       ENDDO

    END SUBROUTINE condensation_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Autoconversion rate (Seifert and Beheng, 2006). Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE autoconversion_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  diss, dzu

       USE cloud_parameters,                                                   &
           ONLY:  hyrho

       USE control_parameters,                                                 &
           ONLY:  microphysics_morrison, rho_surface

       USE grid_variables,                                                     &
           ONLY:  dx, dy

       USE indices,                                                            &
           ONLY:  nzb, nzt, wall_flags_0

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  alpha_cc          !<
       REAL(wp)     ::  autocon           !<
       REAL(wp)     ::  dissipation       !<
       REAL(wp)     ::  flag              !< flag to indicate first grid level above surface
       REAL(wp)     ::  k_au              !<
       REAL(wp)     ::  l_mix             !<
       REAL(wp)     ::  nc_auto           !<
       REAL(wp)     ::  nu_c              !<
       REAL(wp)     ::  phi_au            !<
       REAL(wp)     ::  r_cc              !<
       REAL(wp)     ::  rc                !<
       REAL(wp)     ::  re_lambda         !<
       REAL(wp)     ::  sigma_cc          !<
       REAL(wp)     ::  tau_cloud         !<
       REAL(wp)     ::  xc                !<

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )
          nc_auto = MERGE ( nc_1d(k), nc_const, microphysics_morrison )

          IF ( qc_1d(k) > eps_sb  .AND.  nc_auto > eps_mr )  THEN

             k_au = k_cc / ( 20.0_wp * x0 )
!
!--          Intern time scale of coagulation (Seifert and Beheng, 2006):
!--          (1.0_wp - qc(k,j,i) / ( qc(k,j,i) + qr_1d(k) ))
             tau_cloud = MAX( 1.0_wp - qc_1d(k) / ( qr_1d(k) + qc_1d(k) ),     & 
                              0.0_wp )
!
!--          Universal function for autoconversion process
!--          (Seifert and Beheng, 2006):
             phi_au = 600.0_wp * tau_cloud**0.68_wp * ( 1.0_wp - tau_cloud**0.68_wp )**3
!
!--          Shape parameter of gamma distribution (Geoffroy et al., 2010):
!--          (Use constant nu_c = 1.0_wp instead?)
             nu_c = 1.0_wp !MAX( 0.0_wp, 1580.0_wp * hyrho(k) * qc_1d(k) - 0.28_wp )
!
!--          Mean weight of cloud droplets:
             xc = hyrho(k) * qc_1d(k) / nc_auto
!
!--          Parameterized turbulence effects on autoconversion (Seifert,
!--          Nuijens and Stevens, 2010)
             IF ( collision_turbulence )  THEN
!
!--             Weight averaged radius of cloud droplets:
                rc = 0.5_wp * ( xc * dpirho_l )**( 1.0_wp / 3.0_wp )

                alpha_cc = ( a_1 + a_2 * nu_c ) / ( 1.0_wp + a_3 * nu_c )
                r_cc     = ( b_1 + b_2 * nu_c ) / ( 1.0_wp + b_3 * nu_c )
                sigma_cc = ( c_1 + c_2 * nu_c ) / ( 1.0_wp + c_3 * nu_c )
!
!--             Mixing length (neglecting distance to ground and stratification)
                l_mix = ( dx * dy * dzu(k) )**( 1.0_wp / 3.0_wp )
!
!--             Limit dissipation rate according to Seifert, Nuijens and
!--             Stevens (2010)
                dissipation = MIN( 0.06_wp, diss(k,j,i) )
!
!--             Compute Taylor-microscale Reynolds number:
                re_lambda = 6.0_wp / 11.0_wp *                                 &
                            ( l_mix / c_const )**( 2.0_wp / 3.0_wp ) *         &
                            SQRT( 15.0_wp / kin_vis_air ) *                    &
                            dissipation**( 1.0_wp / 6.0_wp )
!
!--             The factor of 1.0E4 is needed to convert the dissipation rate
!--             from m2 s-3 to cm2 s-3.
                k_au = k_au * ( 1.0_wp +                                       &
                                dissipation * 1.0E4_wp *                       &
                                ( re_lambda * 1.0E-3_wp )**0.25_wp *           &
                                ( alpha_cc * EXP( -1.0_wp * ( ( rc - r_cc ) /  &
                                                  sigma_cc )**2                &
                                                ) + beta_cc                    &
                                )                                              &
                              )
             ENDIF
!
!--          Autoconversion rate (Seifert and Beheng, 2006):
             autocon = k_au * ( nu_c + 2.0_wp ) * ( nu_c + 4.0_wp ) /          &
                       ( nu_c + 1.0_wp )**2 * qc_1d(k)**2 * xc**2 *            &
                       ( 1.0_wp + phi_au / ( 1.0_wp - tau_cloud )**2 ) *       &
                       rho_surface
             autocon = MIN( autocon, qc_1d(k) / dt_micro )

             qr_1d(k) = qr_1d(k) + autocon * dt_micro                 * flag
             qc_1d(k) = qc_1d(k) - autocon * dt_micro                 * flag
             nr_1d(k) = nr_1d(k) + autocon / x0 * hyrho(k) * dt_micro * flag
             IF ( microphysics_morrison ) THEN
                nc_1d(k) = nc_1d(k) - MIN( nc_1d(k), 2.0_wp *                &
                                    autocon / x0 * hyrho(k) * dt_micro * flag )
             ENDIF

          ENDIF

       ENDDO

    END SUBROUTINE autoconversion_ij

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Autoconversion process (Kessler, 1969).
!------------------------------------------------------------------------------!
    SUBROUTINE autoconversion_kessler_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  dzw, prr

       USE cloud_parameters,                                                   &
           ONLY:  l_d_cp, pt_d_t

       USE indices,                                                            &
           ONLY:  nzb, nzt, wall_flags_0

       USE kinds

       USE surface_mod,                                                        &
           ONLY:  get_topography_top_index_ji


       IMPLICIT NONE

       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  k      !<
       INTEGER(iwp) ::  k_wall !< topography top index

       REAL(wp)    ::  dqdt_precip !<
       REAL(wp)    ::  flag              !< flag to indicate first grid level above surface

!
!--    Determine vertical index of topography top
       k_wall = get_topography_top_index_ji( j, i, 's' )
       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

          IF ( qc_1d(k) > ql_crit )  THEN
             dqdt_precip = prec_time_const * ( qc_1d(k) - ql_crit )
          ELSE
             dqdt_precip = 0.0_wp
          ENDIF

          qc_1d(k) = qc_1d(k) - dqdt_precip * dt_micro * flag
          q_1d(k)  = q_1d(k)  - dqdt_precip * dt_micro * flag
          pt_1d(k) = pt_1d(k) + dqdt_precip * dt_micro * l_d_cp * pt_d_t(k) * flag

!
!--       Compute the rain rate (stored on surface grid point)
          prr(k_wall,j,i) = prr(k_wall,j,i) + dqdt_precip * dzw(k) * flag

       ENDDO

    END SUBROUTINE autoconversion_kessler_ij

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Accretion rate (Seifert and Beheng, 2006). Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE accretion_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  diss

       USE cloud_parameters,                                                   &
           ONLY:  hyrho

       USE control_parameters,                                                 &
           ONLY:  microphysics_morrison, rho_surface

       USE indices,                                                            &
           ONLY:  nzb, nzt, wall_flags_0

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  accr              !<
       REAL(wp)     ::  flag              !< flag to indicate first grid level above surface
       REAL(wp)     ::  k_cr              !<
       REAL(wp)     ::  nc_accr           !<
       REAL(wp)     ::  phi_ac            !<
       REAL(wp)     ::  tau_cloud         !<
       REAL(wp)     ::  xc                !<


       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )
          nc_accr = MERGE ( nc_1d(k), nc_const, microphysics_morrison )

          IF ( ( qc_1d(k) > eps_sb )  .AND.  ( qr_1d(k) > eps_sb )  .AND.  &
               ( nc_accr > eps_mr ) )  THEN
!
!--          Intern time scale of coagulation (Seifert and Beheng, 2006):
             tau_cloud = 1.0_wp - qc_1d(k) / ( qc_1d(k) + qr_1d(k) )
!
!--          Universal function for accretion process
!--          (Seifert and Beheng, 2001):
             phi_ac = ( tau_cloud / ( tau_cloud + 5.0E-5_wp ) )**4

!
!--          Mean weight of cloud drops
             xc = MAX( (hyrho(k) * qc_1d(k) / nc_accr), xcmin)
!
!--          Parameterized turbulence effects on autoconversion (Seifert,
!--          Nuijens and Stevens, 2010). The factor of 1.0E4 is needed to
!--          convert the dissipation rate (diss) from m2 s-3 to cm2 s-3.
             IF ( collision_turbulence )  THEN
                k_cr = k_cr0 * ( 1.0_wp + 0.05_wp *                            &
                                 MIN( 600.0_wp,                                &
                                      diss(k,j,i) * 1.0E4_wp )**0.25_wp        &
                               )
             ELSE
                k_cr = k_cr0
             ENDIF
!
!--          Accretion rate (Seifert and Beheng, 2006):
             accr = k_cr * qc_1d(k) * qr_1d(k) * phi_ac *                      &
                    SQRT( rho_surface * hyrho(k) )
             accr = MIN( accr, qc_1d(k) / dt_micro )

             qr_1d(k) = qr_1d(k) + accr * dt_micro * flag
             qc_1d(k) = qc_1d(k) - accr * dt_micro * flag
             IF ( microphysics_morrison )  THEN
                nc_1d(k) = nc_1d(k) - MIN( nc_1d(k), accr / xc *               &
                                             hyrho(k) * dt_micro * flag        &
                                         )
             ENDIF


          ENDIF

       ENDDO

    END SUBROUTINE accretion_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Collisional breakup rate (Seifert, 2008). Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE selfcollection_breakup_ij( i, j )

       USE cloud_parameters,                                                   &
           ONLY:  hyrho

       USE control_parameters,                                                 &
           ONLY:  rho_surface

       USE indices,                                                            &
           ONLY:  nzb, nzt, wall_flags_0

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  breakup           !<
       REAL(wp)     ::  dr                !<
       REAL(wp)     ::  flag              !< flag to indicate first grid level above surface
       REAL(wp)     ::  phi_br            !<
       REAL(wp)     ::  selfcoll          !<

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

          IF ( qr_1d(k) > eps_sb )  THEN
!
!--          Selfcollection rate (Seifert and Beheng, 2001):
             selfcoll = k_rr * nr_1d(k) * qr_1d(k) * SQRT( hyrho(k) * rho_surface )
!
!--          Weight averaged diameter of rain drops:
             dr = ( hyrho(k) * qr_1d(k) / nr_1d(k) * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--          Collisional breakup rate (Seifert, 2008):
             IF ( dr >= 0.3E-3_wp )  THEN
                phi_br  = k_br * ( dr - 1.1E-3_wp )
                breakup = selfcoll * ( phi_br + 1.0_wp )
             ELSE
                breakup = 0.0_wp
             ENDIF

             selfcoll = MAX( breakup - selfcoll, -nr_1d(k) / dt_micro )
             nr_1d(k) = nr_1d(k) + selfcoll * dt_micro * flag

          ENDIF
       ENDDO

    END SUBROUTINE selfcollection_breakup_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Evaporation of precipitable water. Condensation is neglected for
!> precipitable water. Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE evaporation_rain_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  hyp

       USE cloud_parameters,                                                   &
           ONLY:  hyrho, l_d_cp, l_d_r, l_v, r_v, t_d_pt

       USE constants,                                                          &
           ONLY:  pi

       USE indices,                                                            &
           ONLY:  nzb, nzt, wall_flags_0

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  alpha             !<
       REAL(wp)     ::  dr                !<
       REAL(wp)     ::  e_s               !<
       REAL(wp)     ::  evap              !<
       REAL(wp)     ::  evap_nr           !<
       REAL(wp)     ::  f_vent            !<
       REAL(wp)     ::  flag              !< flag to indicate first grid level above surface
       REAL(wp)     ::  g_evap            !<
       REAL(wp)     ::  lambda_r          !<
       REAL(wp)     ::  mu_r              !<
       REAL(wp)     ::  mu_r_2            !<
       REAL(wp)     ::  mu_r_5d2          !<
       REAL(wp)     ::  nr_0              !<
       REAL(wp)     ::  q_s               !<
       REAL(wp)     ::  sat               !<
       REAL(wp)     ::  t_l               !<
       REAL(wp)     ::  temp              !<
       REAL(wp)     ::  xr                !<

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

          IF ( qr_1d(k) > eps_sb )  THEN
!
!--          Actual liquid water temperature:
             t_l = t_d_pt(k) * pt_1d(k)
!
!--          Saturation vapor pressure at t_l:
             e_s = 610.78_wp * EXP( 17.269_wp * ( t_l - 273.16_wp ) /          &
                                    ( t_l - 35.86_wp )                         &
                                  )
!
!--          Computation of saturation mixing ratio:
             q_s   = 0.622_wp * e_s / ( hyp(k) - e_s )
             alpha = 0.622_wp * l_d_r * l_d_cp / ( t_l * t_l )
             q_s   = q_s * ( 1.0_wp + alpha * q_1d(k) ) / ( 1.0_wp + alpha * q_s )
!
!--          Supersaturation:
             sat   = ( q_1d(k) - qr_1d(k) - qc_1d(k) ) / q_s - 1.0_wp
!
!--          Evaporation needs only to be calculated in subsaturated regions
             IF ( sat < 0.0_wp )  THEN
!
!--             Actual temperature:
                temp = t_l + l_d_cp * ( qc_1d(k) + qr_1d(k) )

                g_evap = 1.0_wp / ( ( l_v / ( r_v * temp ) - 1.0_wp ) * l_v /  &
                                    ( thermal_conductivity_l * temp ) +        &
                                    r_v * temp / ( diff_coeff_l * e_s )        &
                                  )
!
!--             Mean weight of rain drops
                xr = hyrho(k) * qr_1d(k) / nr_1d(k)
!
!--             Weight averaged diameter of rain drops:
                dr = ( xr * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--             Compute ventilation factor and intercept parameter
!--             (Seifert and Beheng, 2006; Seifert, 2008):
                IF ( ventilation_effect )  THEN
!
!--                Shape parameter of gamma distribution (Milbrandt and Yau, 2005;
!--                Stevens and Seifert, 2008):
                   mu_r = 10.0_wp * ( 1.0_wp + TANH( 1.2E3_wp * ( dr - 1.4E-3_wp ) ) )
!
!--                Slope parameter of gamma distribution (Seifert, 2008):
                   lambda_r = ( ( mu_r + 3.0_wp ) * ( mu_r + 2.0_wp ) *        &
                                ( mu_r + 1.0_wp )                              &
                              )**( 1.0_wp / 3.0_wp ) / dr

                   mu_r_2   = mu_r + 2.0_wp
                   mu_r_5d2 = mu_r + 2.5_wp

                   f_vent = a_vent * gamm( mu_r_2 ) * lambda_r**( -mu_r_2 ) +  &
                            b_vent * schmidt_p_1d3 *                           &
                            SQRT( a_term / kin_vis_air ) * gamm( mu_r_5d2 ) *  &
                            lambda_r**( -mu_r_5d2 ) *                          &
                            ( 1.0_wp -                                         &
                              0.5_wp * ( b_term / a_term ) *                   &
                              ( lambda_r / ( c_term + lambda_r )               &
                              )**mu_r_5d2 -                                    &
                              0.125_wp * ( b_term / a_term )**2 *              &
                              ( lambda_r / ( 2.0_wp * c_term + lambda_r )      &
                              )**mu_r_5d2 -                                    &
                              0.0625_wp * ( b_term / a_term )**3 *             &
                              ( lambda_r / ( 3.0_wp * c_term + lambda_r )      &
                              )**mu_r_5d2 -                                    &
                              0.0390625_wp * ( b_term / a_term )**4 *          &
                              ( lambda_r / ( 4.0_wp * c_term + lambda_r )      &
                              )**mu_r_5d2                                      &
                            )

                   nr_0   = nr_1d(k) * lambda_r**( mu_r + 1.0_wp ) /           &
                            gamm( mu_r + 1.0_wp )
                ELSE
                   f_vent = 1.0_wp
                   nr_0   = nr_1d(k) * dr
                ENDIF
!
!--             Evaporation rate of rain water content (Seifert and Beheng, 2006):
                evap    = 2.0_wp * pi * nr_0 * g_evap * f_vent * sat / hyrho(k)
                evap    = MAX( evap, -qr_1d(k) / dt_micro )
                evap_nr = MAX( c_evap * evap / xr * hyrho(k),                  &
                               -nr_1d(k) / dt_micro )

                qr_1d(k) = qr_1d(k) + evap    * dt_micro * flag
                nr_1d(k) = nr_1d(k) + evap_nr * dt_micro * flag

             ENDIF
          ENDIF

       ENDDO

    END SUBROUTINE evaporation_rain_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sedimentation of cloud droplets (Ackermann et al., 2009, MWR).
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE sedimentation_cloud_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  ddzu, dzu, prr

       USE cloud_parameters,                                                   &
           ONLY:  hyrho, l_d_cp, pt_d_t

       USE control_parameters,                                                 &
           ONLY:  call_microphysics_at_all_substeps,                           &
                  intermediate_timestep_count, microphysics_morrison

       USE indices,                                                            &
           ONLY:  nzb, nzb, nzt, wall_flags_0

       USE kinds

       USE statistics,                                                         &
           ONLY:  weight_substep

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !<
       INTEGER(iwp) ::  j             !<
       INTEGER(iwp) ::  k             !<

       REAL(wp)     ::  flag    !< flag to indicate first grid level above surface
       REAL(wp)     ::  nc_sedi !<

       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_nc  !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_qc  !<

       sed_qc(nzt+1) = 0.0_wp
       sed_nc(nzt+1) = 0.0_wp


       DO  k = nzt, nzb+1, -1
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )
          nc_sedi = MERGE( nc_1d(k), nc_const, microphysics_morrison )
!
!--       Sedimentation fluxes for number concentration are only calculated
!--       for cloud_scheme = 'morrison'
          IF ( microphysics_morrison ) THEN
             IF ( qc_1d(k) > eps_sb  .AND.  nc_1d(k) > eps_mr )  THEN
                sed_nc(k) = sed_qc_const *                                     &
                            ( qc_1d(k) * hyrho(k) )**( 2.0_wp / 3.0_wp ) *     &
                            ( nc_1d(k) )**( 1.0_wp / 3.0_wp )
             ELSE
                sed_nc(k) = 0.0_wp
             ENDIF

             sed_nc(k) = MIN( sed_nc(k), hyrho(k) * dzu(k+1) *                 &
                              nc_1d(k) / dt_micro + sed_nc(k+1)                &
                            ) * flag

             nc_1d(k) = nc_1d(k) + ( sed_nc(k+1) - sed_nc(k) ) *               &
                                         ddzu(k+1) / hyrho(k) * dt_micro * flag
          ENDIF

          IF ( qc_1d(k) > eps_sb  .AND.  nc_sedi > eps_mr )  THEN
             sed_qc(k) = sed_qc_const * nc_sedi**( -2.0_wp / 3.0_wp ) *        &
                         ( qc_1d(k) * hyrho(k) )**( 5.0_wp / 3.0_wp )  * flag
          ELSE
             sed_qc(k) = 0.0_wp
          ENDIF

          sed_qc(k) = MIN( sed_qc(k), hyrho(k) * dzu(k+1) * q_1d(k) /          &
                                      dt_micro + sed_qc(k+1)                   &
                         ) * flag

          q_1d(k)  = q_1d(k)  + ( sed_qc(k+1) - sed_qc(k) ) * ddzu(k+1) /      &
                                hyrho(k) * dt_micro * flag
          qc_1d(k) = qc_1d(k) + ( sed_qc(k+1) - sed_qc(k) ) * ddzu(k+1) /      &
                                hyrho(k) * dt_micro * flag
          pt_1d(k) = pt_1d(k) - ( sed_qc(k+1) - sed_qc(k) ) * ddzu(k+1) /      &
                                hyrho(k) * l_d_cp * pt_d_t(k) * dt_micro * flag

!
!--       Compute the precipitation rate of cloud (fog) droplets
          IF ( call_microphysics_at_all_substeps )  THEN
             prr(k,j,i) = prr(k,j,i) + sed_qc(k) / hyrho(k) *                  &
                              weight_substep(intermediate_timestep_count) * flag
          ELSE
             prr(k,j,i) = prr(k,j,i) + sed_qc(k) / hyrho(k) * flag
          ENDIF

       ENDDO

    END SUBROUTINE sedimentation_cloud_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computation of sedimentation flux. Implementation according to Stevens
!> and Seifert (2008). Code is based on UCLA-LES. Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE sedimentation_rain_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  ddzu, dzu, prr

       USE cloud_parameters,                                                   &
           ONLY:  hyrho, l_d_cp, pt_d_t

       USE control_parameters,                                                 &
           ONLY:  call_microphysics_at_all_substeps, intermediate_timestep_count

       USE indices,                                                            &
           ONLY:  nzb, nzb, nzt, wall_flags_0

       USE kinds

       USE statistics,                                                         &
           ONLY:  weight_substep

       USE surface_mod,                                                        &
           ONLY :  bc_h

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  k_run         !<
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp)     ::  c_run                      !<
       REAL(wp)     ::  d_max                      !<
       REAL(wp)     ::  d_mean                     !<
       REAL(wp)     ::  d_min                      !<
       REAL(wp)     ::  dr                         !<
       REAL(wp)     ::  flux                       !<
       REAL(wp)     ::  flag                       !< flag to indicate first grid level above surface
       REAL(wp)     ::  lambda_r                   !<
       REAL(wp)     ::  mu_r                       !<
       REAL(wp)     ::  z_run                      !<

       REAL(wp), DIMENSION(nzb:nzt+1) ::  c_nr     !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  c_qr     !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  nr_slope !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  qr_slope !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_nr   !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_qr   !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  w_nr     !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  w_qr     !<

!
!--    Compute velocities
       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

          IF ( qr_1d(k) > eps_sb )  THEN
!
!--          Weight averaged diameter of rain drops:
             dr = ( hyrho(k) * qr_1d(k) / nr_1d(k) * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--          Shape parameter of gamma distribution (Milbrandt and Yau, 2005;
!--          Stevens and Seifert, 2008):
             mu_r = 10.0_wp * ( 1.0_wp + TANH( 1.2E3_wp * ( dr - 1.4E-3_wp ) ) )
!
!--          Slope parameter of gamma distribution (Seifert, 2008):
             lambda_r = ( ( mu_r + 3.0_wp ) * ( mu_r + 2.0_wp ) *              &
                          ( mu_r + 1.0_wp ) )**( 1.0_wp / 3.0_wp ) / dr

             w_nr(k) = MAX( 0.1_wp, MIN( 20.0_wp,                              &
                                         a_term - b_term * ( 1.0_wp +          &
                                            c_term / lambda_r )**( -1.0_wp *   &
                                            ( mu_r + 1.0_wp ) )                &
                                        )                                      &
                          ) * flag
             w_qr(k) = MAX( 0.1_wp, MIN( 20.0_wp,                              &
                                         a_term - b_term * ( 1.0_wp +          &
                                            c_term / lambda_r )**( -1.0_wp *   &
                                            ( mu_r + 4.0_wp ) )                &
                                       )                                       &
                          ) * flag
          ELSE
             w_nr(k) = 0.0_wp
             w_qr(k) = 0.0_wp
          ENDIF
       ENDDO
!
!--    Adjust boundary values using surface data type.
!--    Upward facing non-natural
       surf_s = bc_h(0)%start_index(j,i)
       surf_e = bc_h(0)%end_index(j,i)
       DO  m = surf_s, surf_e
          k         = bc_h(0)%k(m)
          w_nr(k-1) = w_nr(k)
          w_qr(k-1) = w_qr(k)
       ENDDO
!
!--    Downward facing non-natural
       surf_s = bc_h(1)%start_index(j,i)
       surf_e = bc_h(1)%end_index(j,i)
       DO  m = surf_s, surf_e
          k         = bc_h(1)%k(m)
          w_nr(k+1) = w_nr(k)
          w_qr(k+1) = w_qr(k)
       ENDDO
!
!--    Neumann boundary condition at model top
       w_nr(nzt+1) = 0.0_wp
       w_qr(nzt+1) = 0.0_wp
!
!--    Compute Courant number
       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

          c_nr(k) = 0.25_wp * ( w_nr(k-1) + 2.0_wp * w_nr(k) + w_nr(k+1) ) *   &
                    dt_micro * ddzu(k) * flag
          c_qr(k) = 0.25_wp * ( w_qr(k-1) + 2.0_wp * w_qr(k) + w_qr(k+1) ) *   &
                    dt_micro * ddzu(k) * flag
       ENDDO
!
!--    Limit slopes with monotonized centered (MC) limiter (van Leer, 1977):
       IF ( limiter_sedimentation )  THEN

          DO k = nzb+1, nzt
!
!--          Predetermine flag to mask topography
             flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

             d_mean = 0.5_wp * ( qr_1d(k+1) - qr_1d(k-1) )
             d_min  = qr_1d(k) - MIN( qr_1d(k+1), qr_1d(k), qr_1d(k-1) )
             d_max  = MAX( qr_1d(k+1), qr_1d(k), qr_1d(k-1) ) - qr_1d(k)

             qr_slope(k) = SIGN(1.0_wp, d_mean) * MIN ( 2.0_wp * d_min,        &
                                                        2.0_wp * d_max,        &
                                                        ABS( d_mean ) ) * flag

             d_mean = 0.5_wp * ( nr_1d(k+1) - nr_1d(k-1) )
             d_min  = nr_1d(k) - MIN( nr_1d(k+1), nr_1d(k), nr_1d(k-1) )
             d_max  = MAX( nr_1d(k+1), nr_1d(k), nr_1d(k-1) ) - nr_1d(k)

             nr_slope(k) = SIGN(1.0_wp, d_mean) * MIN ( 2.0_wp * d_min,        &
                                                        2.0_wp * d_max,        &
                                                        ABS( d_mean ) ) * flag
          ENDDO

       ELSE

          nr_slope = 0.0_wp
          qr_slope = 0.0_wp

       ENDIF

       sed_nr(nzt+1) = 0.0_wp
       sed_qr(nzt+1) = 0.0_wp
!
!--    Compute sedimentation flux
       DO  k = nzt, nzb+1, -1
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )
!
!--       Sum up all rain drop number densities which contribute to the flux
!--       through k-1/2
          flux  = 0.0_wp
          z_run = 0.0_wp ! height above z(k)
          k_run = k
          c_run = MIN( 1.0_wp, c_nr(k) )
          DO WHILE ( c_run > 0.0_wp  .AND.  k_run <= nzt )
             flux  = flux + hyrho(k_run) *                                     &
                     ( nr_1d(k_run) + nr_slope(k_run) * ( 1.0_wp - c_run ) *   &
                     0.5_wp ) * c_run * dzu(k_run) * flag
             z_run = z_run + dzu(k_run)            * flag
             k_run = k_run + 1                     * flag
             c_run = MIN( 1.0_wp, c_nr(k_run) - z_run * ddzu(k_run) ) * flag
          ENDDO
!
!--       It is not allowed to sediment more rain drop number density than
!--       available
          flux = MIN( flux,                                                    &
                      hyrho(k) * dzu(k+1) * nr_1d(k) + sed_nr(k+1) * dt_micro )

          sed_nr(k) = flux / dt_micro * flag
          nr_1d(k)  = nr_1d(k) + ( sed_nr(k+1) - sed_nr(k) ) * ddzu(k+1) /     &
                                    hyrho(k) * dt_micro * flag
!
!--       Sum up all rain water content which contributes to the flux
!--       through k-1/2
          flux  = 0.0_wp
          z_run = 0.0_wp ! height above z(k)
          k_run = k
          c_run = MIN( 1.0_wp, c_qr(k) )

          DO WHILE ( c_run > 0.0_wp  .AND.  k_run <= nzt )

             flux  = flux + hyrho(k_run) *                                     &
                     ( qr_1d(k_run) + qr_slope(k_run) * ( 1.0_wp - c_run ) *   &
                     0.5_wp ) * c_run * dzu(k_run) * flag
             z_run = z_run + dzu(k_run)            * flag
             k_run = k_run + 1                     * flag
             c_run = MIN( 1.0_wp, c_qr(k_run) - z_run * ddzu(k_run) ) * flag

          ENDDO
!
!--       It is not allowed to sediment more rain water content than available
          flux = MIN( flux,                                                    &
                      hyrho(k) * dzu(k) * qr_1d(k) + sed_qr(k+1) * dt_micro )

          sed_qr(k) = flux / dt_micro * flag

          qr_1d(k) = qr_1d(k) + ( sed_qr(k+1) - sed_qr(k) ) * ddzu(k+1) /      &
                                hyrho(k) * dt_micro * flag
          q_1d(k)  = q_1d(k)  + ( sed_qr(k+1) - sed_qr(k) ) * ddzu(k+1) /      &
                                hyrho(k) * dt_micro * flag
          pt_1d(k) = pt_1d(k) - ( sed_qr(k+1) - sed_qr(k) ) * ddzu(k+1) /      &
                                hyrho(k) * l_d_cp * pt_d_t(k) * dt_micro * flag
!
!--       Compute the rain rate
          IF ( call_microphysics_at_all_substeps )  THEN
             prr(k,j,i) = prr(k,j,i) + sed_qr(k) / hyrho(k)                    &
                          * weight_substep(intermediate_timestep_count) * flag
          ELSE
             prr(k,j,i) = prr(k,j,i) + sed_qr(k) / hyrho(k) * flag
          ENDIF

       ENDDO

    END SUBROUTINE sedimentation_rain_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine computes the precipitation amount due to gravitational
!> settling of rain and cloud (fog) droplets
!------------------------------------------------------------------------------!
    SUBROUTINE calc_precipitation_amount_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  precipitation_amount, prr

       USE cloud_parameters,                                                   &
           ONLY:  hyrho

       USE control_parameters,                                                 &
           ONLY:  call_microphysics_at_all_substeps, dt_do2d_xy, dt_3d,        &
                  intermediate_timestep_count, intermediate_timestep_count_max,&
                  precipitation_amount_interval, time_do2d_xy

       USE indices,                                                            &
           ONLY:  nzb, nzt, wall_flags_0

       USE kinds

       USE surface_mod,                                                        &
           ONLY :  bc_h

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       IF ( ( dt_do2d_xy - time_do2d_xy ) < precipitation_amount_interval .AND.&
            ( .NOT. call_microphysics_at_all_substeps .OR.                     &
            intermediate_timestep_count == intermediate_timestep_count_max ) ) &
       THEN

          surf_s = bc_h(0)%start_index(j,i)
          surf_e = bc_h(0)%end_index(j,i)
          DO  m = surf_s, surf_e
             k                         = bc_h(0)%k(m)
             precipitation_amount(j,i) = precipitation_amount(j,i) +           &
                                               prr(k,j,i) * hyrho(k) * dt_3d
          ENDDO

       ENDIF

    END SUBROUTINE calc_precipitation_amount_ij

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This function computes the gamma function (Press et al., 1992).
!> The gamma function is needed for the calculation of the evaporation
!> of rain drops.
!------------------------------------------------------------------------------!
    FUNCTION gamm( xx )

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  j            !<

       REAL(wp)     ::  gamm         !<
       REAL(wp)     ::  ser          !<
       REAL(wp)     ::  tmp          !<
       REAL(wp)     ::  x_gamm       !<
       REAL(wp)     ::  xx           !<
       REAL(wp)     ::  y_gamm       !<


       REAL(wp), PARAMETER  ::  stp = 2.5066282746310005_wp               !<
       REAL(wp), PARAMETER  ::  cof(6) = (/ 76.18009172947146_wp,      &
                                           -86.50532032941677_wp,      &
                                            24.01409824083091_wp,      &
                                            -1.231739572450155_wp,     &
                                             0.1208650973866179E-2_wp, &
                                            -0.5395239384953E-5_wp /)     !<

       x_gamm = xx
       y_gamm = x_gamm
       tmp = x_gamm + 5.5_wp
       tmp = ( x_gamm + 0.5_wp ) * LOG( tmp ) - tmp
       ser = 1.000000000190015_wp

       DO  j = 1, 6
          y_gamm = y_gamm + 1.0_wp
          ser    = ser + cof( j ) / y_gamm
       ENDDO

!
!--    Until this point the algorithm computes the logarithm of the gamma
!--    function. Hence, the exponential function is used.
!       gamm = EXP( tmp + LOG( stp * ser / x_gamm ) )
       gamm = EXP( tmp ) * stp * ser / x_gamm

       RETURN

    END FUNCTION gamm

 END MODULE microphysics_mod
