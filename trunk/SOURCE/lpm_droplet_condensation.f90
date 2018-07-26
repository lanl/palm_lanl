!> @file lpm_droplet_condensation.f90
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
! $Id: lpm_droplet_condensation.f90 3049 2018-05-29 13:52:36Z Giersch $
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
! 
! 3039 2018-05-24 13:13:11Z schwenkel
! bugfix for lcm with grid stretching
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2608 2017-11-13 14:04:26Z schwenkel
! Calculation of magnus equation in external module (diagnostic_quantities_mod).
! 
! 2375 2017-08-29 14:10:28Z schwenkel
! Changed ONLY-dependencies
!
! 2312 2017-07-14 20:26:51Z hoffmann
! Rosenbrock scheme improved. Gas-kinetic effect added.
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1890 2016-04-22 08:52:11Z hoffmann
! Some improvements of the Rosenbrock method. If the Rosenbrock method needs more
! than 40 iterations to find a sufficient time setp, the model is not aborted.
! This might lead to small erros. But they can be assumend as negligible, since
! the maximum timesetp of the Rosenbrock method after 40 iterations will be
! smaller than 10^-16 s.
!
! 1871 2016-04-15 11:46:09Z hoffmann
! Initialization of aerosols added.
!
! 1849 2016-04-08 11:33:18Z hoffmann
! Interpolation of supersaturation has been removed because it is not in
! accordance with the release/depletion of latent heat/water vapor in
! interaction_droplets_ptq.
! Calculation of particle Reynolds number has been corrected.
! eps_ros added from modules.
!
! 1831 2016-04-07 13:15:51Z hoffmann
! curvature_solution_effects moved to particle_attributes
!
! 1822 2016-04-07 07:49:42Z hoffmann
! Unused variables removed.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1359 2014-04-11 17:15:14Z hoffmann
! New particle structure integrated.
! Kind definition added to all floating point numbers.
!
! 1346 2014-03-27 13:18:20Z heinze
! Bugfix: REAL constants provided with KIND-attribute especially in call of
! intrinsic function like MAX, MIN, SIGN
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
! 1318 2014-03-17 13:35:16Z raasch
! module interfaces removed
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1071 2012-11-29 16:54:55Z franke
! Ventilation effect for evaporation of large droplets included
! Check for unreasonable results included in calculation of Rosenbrock method
! since physically unlikely results were observed and for the same
! reason the first internal time step in Rosenbrock method should be < 1.0E02 in
! case of evaporation
! Unnecessary calculation of ql_int removed
! Unnecessary calculations in Rosenbrock method (d2rdt2, drdt_m, dt_ros_last)
! removed
! Bugfix: factor in calculation of surface tension changed from 0.00155 to
! 0.000155
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 849 2012-03-15 10:35:09Z raasch
! initial revision (former part of advec_particles)
!
!
! Description:
! ------------
!> Calculates change in droplet radius by condensation/evaporation, using
!> either an analytic formula or by numerically integrating the radius growth
!> equation including curvature and solution effects using Rosenbrocks method
!> (see Numerical recipes in FORTRAN, 2nd edition, p. 731).
!> The analytical formula and growth equation follow those given in
!> Rogers and Yau (A short course in cloud physics, 3rd edition, p. 102/103).
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_droplet_condensation (ip,jp,kp)


    USE arrays_3d,                                                             &
        ONLY:  dzw, hyp, pt, q, ql_c, ql_v

    USE cloud_parameters,                                                      &
        ONLY:  l_d_rv, l_v, molecular_weight_of_solute,                        &
               molecular_weight_of_water, rho_l, rho_s, r_v, vanthoff

    USE constants,                                                             &
        ONLY:  pi

    USE control_parameters,                                                    &
        ONLY:  dt_3d, dz, message_string, molecular_viscosity, rho_surface

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE diagnostic_quantities_mod,                                             &
        ONLY:  magnus

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE lpm_collision_kernels_mod,                                             &
        ONLY:  rclass_lbound, rclass_ubound

    USE kinds

    USE particle_attributes,                                                   &
        ONLY:  curvature_solution_effects, hall_kernel, number_of_particles,   &
               particles, radius_classes, use_kernel_tables, wang_kernel

    IMPLICIT NONE

    INTEGER(iwp) :: ip                         !<
    INTEGER(iwp) :: jp                         !<
    INTEGER(iwp) :: kp                         !<
    INTEGER(iwp) :: n                          !<

    REAL(wp) ::  afactor                       !< curvature effects
    REAL(wp) ::  arg                           !<
    REAL(wp) ::  bfactor                       !< solute effects
    REAL(wp) ::  ddenom                        !<
    REAL(wp) ::  delta_r                       !<
    REAL(wp) ::  diameter                      !< diameter of cloud droplets
    REAL(wp) ::  diff_coeff                    !< diffusivity for water vapor
    REAL(wp) ::  drdt                          !<
    REAL(wp) ::  dt_ros                        !<
    REAL(wp) ::  dt_ros_sum                    !<
    REAL(wp) ::  d2rdtdr                       !<
    REAL(wp) ::  e_a                           !< current vapor pressure
    REAL(wp) ::  e_s                           !< current saturation vapor pressure
    REAL(wp) ::  error                         !< local truncation error in Rosenbrock
    REAL(wp) ::  k1                            !<
    REAL(wp) ::  k2                            !<
    REAL(wp) ::  r_err                         !< First order estimate of Rosenbrock radius
    REAL(wp) ::  r_ros                         !< Rosenbrock radius
    REAL(wp) ::  r_ros_ini                     !< initial Rosenbrock radius
    REAL(wp) ::  r0                            !< gas-kinetic lengthscale
    REAL(wp) ::  sigma                         !< surface tension of water
    REAL(wp) ::  thermal_conductivity          !< thermal conductivity for water
    REAL(wp) ::  t_int                         !< temperature
    REAL(wp) ::  w_s                           !< terminal velocity of droplets
    REAL(wp) ::  re_p                          !< particle Reynolds number
!
!-- Parameters for Rosenbrock method (see Verwer et al., 1999)
    REAL(wp), PARAMETER :: prec = 1.0E-3_wp     !< precision of Rosenbrock solution
    REAL(wp), PARAMETER :: q_increase = 1.5_wp  !< increase factor in timestep
    REAL(wp), PARAMETER :: q_decrease = 0.9_wp  !< decrease factor in timestep
    REAL(wp), PARAMETER :: gamma = 0.292893218814_wp !< = 1.0 - 1.0 / SQRT(2.0)
!
!-- Parameters for terminal velocity
    REAL(wp), PARAMETER ::  a_rog = 9.65_wp      !< parameter for fall velocity
    REAL(wp), PARAMETER ::  b_rog = 10.43_wp     !< parameter for fall velocity
    REAL(wp), PARAMETER ::  c_rog = 0.6_wp       !< parameter for fall velocity
    REAL(wp), PARAMETER ::  k_cap_rog = 4.0_wp   !< parameter for fall velocity
    REAL(wp), PARAMETER ::  k_low_rog = 12.0_wp  !< parameter for fall velocity
    REAL(wp), PARAMETER ::  d0_rog = 0.745_wp    !< separation diameter

    REAL(wp), DIMENSION(number_of_particles) ::  ventilation_effect     !<
    REAL(wp), DIMENSION(number_of_particles) ::  new_r                  !<

    CALL cpu_log( log_point_s(42), 'lpm_droplet_condens', 'start' )

!
!-- Absolute temperature
    t_int = pt(kp,jp,ip) * ( hyp(kp) / 100000.0_wp )**0.286_wp
!
!-- Saturation vapor pressure (Eq. 10 in Bolton, 1980)
    e_s = magnus( t_int )
!
!-- Current vapor pressure
    e_a = q(kp,jp,ip) * hyp(kp) / ( q(kp,jp,ip) + 0.622_wp )
!
!-- Thermal conductivity for water (from Rogers and Yau, Table 7.1)
    thermal_conductivity = 7.94048E-05_wp * t_int + 0.00227011_wp
!
!-- Moldecular diffusivity of water vapor in air (Hall und Pruppacher, 1976)
    diff_coeff           = 0.211E-4_wp * ( t_int / 273.15_wp )**1.94_wp * &
                           ( 101325.0_wp / hyp(kp) )
!
!-- Lengthscale for gas-kinetic effects (from Mordy, 1959, p. 23):
    r0 = diff_coeff / 0.036_wp * SQRT( 2.0_wp * pi / ( r_v * t_int ) )
!
!-- Calculate effects of heat conductivity and diffusion of water vapor on the
!-- diffusional growth process (usually known as 1.0 / (F_k + F_d) )
    ddenom  = 1.0_wp / ( rho_l * r_v * t_int / ( e_s * diff_coeff ) +          &
                         ( l_v / ( r_v * t_int ) - 1.0_wp ) * rho_l *          &
                         l_v / ( thermal_conductivity * t_int )                &
                       )
    new_r = 0.0_wp
!
!-- Determine ventilation effect on evaporation of large drops
    DO  n = 1, number_of_particles

       IF ( particles(n)%radius >= 4.0E-5_wp  .AND.  e_a / e_s < 1.0_wp )  THEN
!
!--       Terminal velocity is computed for vertical direction (Rogers et al.,
!--       1993, J. Appl. Meteorol.)
          diameter = particles(n)%radius * 2000.0_wp !diameter in mm
          IF ( diameter <= d0_rog )  THEN
             w_s = k_cap_rog * diameter * ( 1.0_wp - EXP( -k_low_rog * diameter ) )
          ELSE
             w_s = a_rog - b_rog * EXP( -c_rog * diameter )
          ENDIF
!
!--       Calculate droplet's Reynolds number
          re_p = 2.0_wp * particles(n)%radius * w_s / molecular_viscosity
!
!--       Ventilation coefficient (Rogers and Yau, 1989):
          IF ( re_p > 2.5_wp )  THEN
             ventilation_effect(n) = 0.78_wp + 0.28_wp * SQRT( re_p )
          ELSE
             ventilation_effect(n) = 1.0_wp + 0.09_wp * re_p
          ENDIF
       ELSE
!
!--       For small droplets or in supersaturated environments, the ventilation
!--       effect does not play a role
          ventilation_effect(n) = 1.0_wp
       ENDIF
    ENDDO

    IF( .NOT. curvature_solution_effects ) then
!
!--    Use analytic model for diffusional growth including gas-kinetic
!--    effects (Mordy, 1959) but without the impact of aerosols.
       DO  n = 1, number_of_particles
          arg      = ( particles(n)%radius + r0 )**2 + 2.0_wp * dt_3d * ddenom * &
                                                       ventilation_effect(n) *   &
                                                       ( e_a / e_s - 1.0_wp )
          arg      = MAX( arg, ( 0.01E-6 + r0 )**2 )
          new_r(n) = SQRT( arg ) - r0
       ENDDO

    ELSE
!
!--    Integrate the diffusional growth including gas-kinetic (Mordy, 1959),
!--    as well as curvature and solute effects (e.g., Köhler, 1936).
!
!--    Curvature effect (afactor) with surface tension (sigma) by Straka (2009)
       sigma = 0.0761_wp - 0.000155_wp * ( t_int - 273.15_wp )
!
!--    Solute effect (afactor)
       afactor = 2.0_wp * sigma / ( rho_l * r_v * t_int )

       DO  n = 1, number_of_particles
!
!--       Solute effect (bfactor)
          bfactor = vanthoff * rho_s * particles(n)%aux1**3 *                    &
                    molecular_weight_of_water / ( rho_l * molecular_weight_of_solute )

          dt_ros     = particles(n)%aux2  ! use previously stored Rosenbrock timestep
          dt_ros_sum = 0.0_wp

          r_ros     = particles(n)%radius  ! initialize Rosenbrock particle radius
          r_ros_ini = r_ros
!
!--       Integrate growth equation using a 2nd-order Rosenbrock method
!--       (see Verwer et al., 1999, Eq. (3.2)). The Rosenbrock method adjusts
!--       its with internal timestep to minimize the local truncation error.
          DO WHILE ( dt_ros_sum < dt_3d )

             dt_ros = MIN( dt_ros, dt_3d - dt_ros_sum )

             DO

                drdt = ddenom * ventilation_effect(n) * ( e_a / e_s - 1.0 -    &
                                                          afactor / r_ros +    &
                                                          bfactor / r_ros**3   &
                                                        ) / ( r_ros + r0 )

                d2rdtdr = -ddenom * ventilation_effect(n) * (                  &
                                                (e_a / e_s - 1.0) * r_ros**4 - &
                                                afactor * r0 * r_ros**2 -      &
                                                2.0 * afactor * r_ros**3 +     &
                                                3.0 * bfactor * r0 +           &
                                                4.0 * bfactor * r_ros          &
                                                            )                  &
                          / ( r_ros**4 * ( r_ros + r0 )**2 )

                k1    = drdt / ( 1.0 - gamma * dt_ros * d2rdtdr )

                r_ros = MAX(r_ros_ini + k1 * dt_ros, particles(n)%aux1)
                r_err = r_ros

                drdt  = ddenom * ventilation_effect(n) * ( e_a / e_s - 1.0 -   &
                                                           afactor / r_ros +   &
                                                           bfactor / r_ros**3  &
                                                         ) / ( r_ros + r0 )

                k2 = ( drdt - dt_ros * 2.0 * gamma * d2rdtdr * k1 ) / &
                     ( 1.0 - dt_ros * gamma * d2rdtdr )

                r_ros = MAX(r_ros_ini + dt_ros * ( 1.5 * k1 + 0.5 * k2), particles(n)%aux1)
   !
   !--          Check error of the solution, and reduce dt_ros if necessary.
                error = ABS(r_err - r_ros) / r_ros
                IF ( error .GT. prec )  THEN
                   dt_ros = SQRT( q_decrease * prec / error ) * dt_ros
                   r_ros  = r_ros_ini
                ELSE
                   dt_ros_sum = dt_ros_sum + dt_ros
                   dt_ros     = q_increase * dt_ros
                   r_ros_ini  = r_ros
                   EXIT
                ENDIF

             END DO

          END DO !Rosenbrock loop
!
!--       Store new particle radius
          new_r(n) = r_ros
!
!--       Store internal time step value for next PALM step
          particles(n)%aux2 = dt_ros

       ENDDO !Particle loop

    ENDIF

    DO  n = 1, number_of_particles
!
!--    Sum up the change in liquid water for the respective grid
!--    box for the computation of the release/depletion of water vapor
!--    and heat.
       ql_c(kp,jp,ip) = ql_c(kp,jp,ip) + particles(n)%weight_factor *          &
                                   rho_l * 1.33333333_wp * pi *                &
                                   ( new_r(n)**3 - particles(n)%radius**3 ) /  &
                                   ( rho_surface * dx * dy * dzw(kp) )
!
!--    Check if the increase in liqid water is not too big. If this is the case,
!--    the model timestep might be too long.
       IF ( ql_c(kp,jp,ip) > 100.0_wp )  THEN
          WRITE( message_string, * ) 'k=',kp,' j=',jp,' i=',ip,                &
                       ' ql_c=',ql_c(kp,jp,ip), '&part(',n,')%wf=',            &
                       particles(n)%weight_factor,' delta_r=',delta_r
          CALL message( 'lpm_droplet_condensation', 'PA0143', 2, 2, -1, 6, 1 )
       ENDIF
!
!--    Check if the change in the droplet radius is not too big. If this is the
!--    case, the model timestep might be too long.
       delta_r = new_r(n) - particles(n)%radius
       IF ( delta_r < 0.0_wp  .AND. new_r(n) < 0.0_wp )  THEN
          WRITE( message_string, * ) '#1 k=',kp,' j=',jp,' i=',ip,             &
                       ' e_s=',e_s, ' e_a=',e_a,' t_int=',t_int,               &
                       '&delta_r=',delta_r,                                    &
                       ' particle_radius=',particles(n)%radius
          CALL message( 'lpm_droplet_condensation', 'PA0144', 2, 2, -1, 6, 1 )
       ENDIF
!
!--    Sum up the total volume of liquid water (needed below for
!--    re-calculating the weighting factors)
       ql_v(kp,jp,ip) = ql_v(kp,jp,ip) + particles(n)%weight_factor * new_r(n)**3
!
!--    Determine radius class of the particle needed for collision
       IF ( use_kernel_tables )  THEN
          particles(n)%class = ( LOG( new_r(n) ) - rclass_lbound ) /           &
                               ( rclass_ubound - rclass_lbound ) *             &
                               radius_classes
          particles(n)%class = MIN( particles(n)%class, radius_classes )
          particles(n)%class = MAX( particles(n)%class, 1 )
       ENDIF
 !
 !--   Store new radius to particle features
       particles(n)%radius = new_r(n)

    ENDDO

    CALL cpu_log( log_point_s(42), 'lpm_droplet_condens', 'stop' )


 END SUBROUTINE lpm_droplet_condensation
