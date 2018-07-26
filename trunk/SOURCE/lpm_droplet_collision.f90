!> @file lpm_droplet_collision.f90
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
! $Id: lpm_droplet_collision.f90 3039 2018-05-24 13:13:11Z schwenkel $
! bugfix for lcm with grid stretching
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2375 2017-08-29 14:10:28Z schwenkel
! Changed ONLY-dependencies
!
! 2312 2017-07-14 20:26:51Z hoffmann
! Consideration of aerosol mass during collision. Average impact algorithm has
! been removed.
!
! 2274 2017-06-09 13:27:48Z Giersch
! Changed error messages
!
! 2123 2017-01-18 12:34:59Z hoffmann
!
! 2122 2017-01-18 12:22:54Z hoffmann
! Some reformatting of the code.
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1884 2016-04-21 11:11:40Z hoffmann
! Conservation of mass should only be checked if collisions took place.
!
! 1860 2016-04-13 13:21:28Z hoffmann
! Interpolation of dissipation rate adjusted to more reasonable values.
!
! 1822 2016-04-07 07:49:42Z hoffmann
! Integration of a new collision algortithm based on Shima et al. (2009) and
! Soelch and Kaercher (2010) called all_or_nothing. The previous implemented
! collision algorithm is called average_impact. Moreover, both algorithms are
! now positive definit due to their construction, i.e., no negative weighting
! factors should occur.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1359 2014-04-11 17:15:14Z hoffmann
! New particle structure integrated.
! Kind definition added to all floating point numbers.
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
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1071 2012-11-29 16:54:55Z franke
! Calculation of Hall and Wang kernel now uses collision-coalescence formulation
! proposed by Wang instead of the continuous collection equation (for more
! information about new method see PALM documentation)
! Bugfix: message identifiers added
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
!> Calculates change in droplet radius by collision. Droplet collision is
!> calculated for each grid box seperately. Collision is parameterized by
!> using collision kernels. Two different kernels are available:
!> Hall kernel: Kernel from Hall (1980, J. Atmos. Sci., 2486-2507), which
!>              considers collision due to pure gravitational effects.
!> Wang kernel: Beside gravitational effects (treated with the Hall-kernel) also
!>              the effects of turbulence on the collision are considered using
!>              parameterizations of Ayala et al. (2008, New J. Phys., 10,
!>              075015) and Wang and Grabowski (2009, Atmos. Sci. Lett., 10,
!>              1-8). This kernel includes three possible effects of turbulence:
!>              the modification of the relative velocity between the droplets,
!>              the effect of preferential concentration, and the enhancement of
!>              collision efficiencies.
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_droplet_collision (i,j,k)

    USE arrays_3d,                                                             &
        ONLY:  diss, dzw, ql_v, ql_vp

    USE cloud_parameters,                                                      &
        ONLY:  rho_l, rho_s

    USE constants,                                                             &
        ONLY:  pi

    USE control_parameters,                                                    &
        ONLY:  dt_3d, message_string, dz

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE kinds

    USE lpm_collision_kernels_mod,                                             &
        ONLY:  ckernel, recalculate_kernel

    USE particle_attributes,                                                   &
        ONLY:  curvature_solution_effects, dissipation_classes, hall_kernel,   &
               iran_part, number_of_particles, particles, particle_type,       &
               prt_count, use_kernel_tables, wang_kernel

    USE random_function_mod,                                                   &
        ONLY:  random_function

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  eclass   !<
    INTEGER(iwp) ::  i        !<
    INTEGER(iwp) ::  j        !<
    INTEGER(iwp) ::  k        !<
    INTEGER(iwp) ::  n        !<
    INTEGER(iwp) ::  m        !<
    INTEGER(iwp) ::  rclass_l !<
    INTEGER(iwp) ::  rclass_s !<

    REAL(wp) ::  collection_probability  !< probability for collection
    REAL(wp) ::  ddV                     !< inverse grid box volume
    REAL(wp) ::  epsilon                 !< dissipation rate
    REAL(wp) ::  factor_volume_to_mass   !< 4.0 / 3.0 * pi * rho_l
    REAL(wp) ::  xm                      !< droplet mass of super-droplet m
    REAL(wp) ::  xn                      !< droplet mass of super-droplet n
    REAL(wp) ::  xsm                     !< aerosol mass of super-droplet m
    REAL(wp) ::  xsn                     !< aerosol mass of super-droplet n

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  weight    !< weighting factor
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  mass      !< total mass of super droplet
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  aero_mass !< total aerosol mass of super droplet

    CALL cpu_log( log_point_s(43), 'lpm_droplet_coll', 'start' )

    number_of_particles   = prt_count(k,j,i)
    factor_volume_to_mass = 4.0_wp / 3.0_wp * pi * rho_l
    ddV                   = 1.0_wp / ( dx * dy * dzw(k) )
!
!-- Collision requires at least one super droplet inside the box
    IF ( number_of_particles > 0 )  THEN

       IF ( use_kernel_tables )  THEN
!
!--       Fast method with pre-calculated collection kernels for
!--       discrete radius- and dissipation-classes.
          IF ( wang_kernel )  THEN
             eclass = INT( diss(k,j,i) * 1.0E4_wp / 600.0_wp * &
                           dissipation_classes ) + 1
             epsilon = diss(k,j,i)
          ELSE
             epsilon = 0.0_wp
          ENDIF

          IF ( hall_kernel  .OR.  epsilon * 1.0E4_wp < 0.001_wp )  THEN
             eclass = 0   ! Hall kernel is used
          ELSE
             eclass = MIN( dissipation_classes, eclass )
          ENDIF

       ELSE
!
!--       Collection kernels are re-calculated for every new
!--       grid box. First, allocate memory for kernel table.
!--       Third dimension is 1, because table is re-calculated for
!--       every new dissipation value.
          ALLOCATE( ckernel(1:number_of_particles,1:number_of_particles,1:1) )
!
!--       Now calculate collection kernel for this box. Note that
!--       the kernel is based on the previous time step
          CALL recalculate_kernel( i, j, k )

       ENDIF
!
!--    Temporary fields for total mass of super-droplet, aerosol mass, and
!--    weighting factor are allocated.
       ALLOCATE(mass(1:number_of_particles), weight(1:number_of_particles))
       IF ( curvature_solution_effects )  ALLOCATE(aero_mass(1:number_of_particles))

       mass(1:number_of_particles)   = particles(1:number_of_particles)%weight_factor * &
                                       particles(1:number_of_particles)%radius**3     * &
                                       factor_volume_to_mass

       weight(1:number_of_particles) = particles(1:number_of_particles)%weight_factor

       IF ( curvature_solution_effects )  THEN
          aero_mass(1:number_of_particles) = particles(1:number_of_particles)%weight_factor * &
                                             particles(1:number_of_particles)%aux1**3       * &
                                             4.0 / 3.0 * pi * rho_s
       ENDIF
!
!--    Calculate collision/coalescence
       DO  n = 1, number_of_particles

          DO  m = n, number_of_particles
!
!--          For collisions, the weighting factor of at least one super-droplet
!--          needs to be larger or equal to one.
             IF ( MIN( weight(n), weight(m) ) .LT. 1.0 )  CYCLE
!
!--          Get mass of individual droplets (aerosols)
             xn = mass(n) / weight(n)
             xm = mass(m) / weight(m)
             IF ( curvature_solution_effects )  THEN
                xsn = aero_mass(n) / weight(n)
                xsm = aero_mass(m) / weight(m)
             ENDIF
!
!--          Probability that the necessary collisions take place
             IF ( use_kernel_tables )  THEN
                rclass_l = particles(n)%class
                rclass_s = particles(m)%class

                collection_probability  = MAX( weight(n), weight(m) ) *     &
                                          ckernel(rclass_l,rclass_s,eclass) * ddV * dt_3d
             ELSE
                collection_probability  = MAX( weight(n), weight(m) ) *     &
                                          ckernel(n,m,1) * ddV * dt_3d
             ENDIF
!
!--          Calculate the number of collections and consider multiple collections.
!--          (Accordingly, p_crit will be 0.0, 1.0, 2.0, ...)
             IF ( collection_probability - FLOOR(collection_probability)    &
                  .GT. random_function( iran_part ) )  THEN
                collection_probability = FLOOR(collection_probability) + 1.0_wp
             ELSE
                collection_probability = FLOOR(collection_probability)
             ENDIF

             IF ( collection_probability .GT. 0.0 )  THEN
!
!--             Super-droplet n collects droplets of super-droplet m
                IF ( weight(n) .LT. weight(m) )  THEN

                   mass(n)   = mass(n)   + weight(n) * xm * collection_probability
                   weight(m) = weight(m) - weight(n)      * collection_probability
                   mass(m)   = mass(m)   - weight(n) * xm * collection_probability
                   IF ( curvature_solution_effects )  THEN
                      aero_mass(n) = aero_mass(n) + weight(n) * xsm * collection_probability
                      aero_mass(m) = aero_mass(m) - weight(n) * xsm * collection_probability
                   ENDIF

                ELSEIF ( weight(m) .LT. weight(n) )  THEN

                   mass(m)   = mass(m)   + weight(m) * xn * collection_probability
                   weight(n) = weight(n) - weight(m)      * collection_probability
                   mass(n)   = mass(n)   - weight(m) * xn * collection_probability
                   IF ( curvature_solution_effects )  THEN
                      aero_mass(m) = aero_mass(m) + weight(m) * xsn * collection_probability
                      aero_mass(n) = aero_mass(n) - weight(m) * xsn * collection_probability
                   ENDIF

                ELSE
!
!--                Collisions of particles of the same weighting factor.
!--                Particle n collects 1/2 weight(n) droplets of particle m,
!--                particle m collects 1/2 weight(m) droplets of particle n.
!--                The total mass mass changes accordingly.
!--                If n = m, the first half of the droplets coalesces with the
!--                second half of the droplets; mass is unchanged because
!--                xm = xn for n = m.
!--
!--                Note: For m = n this equation is an approximation only
!--                valid for weight >> 1 (which is usually the case). The
!--                approximation is weight(n)-1 = weight(n).
                   mass(n)   = mass(n)   + 0.5_wp * weight(n) * ( xm - xn )
                   mass(m)   = mass(m)   + 0.5_wp * weight(m) * ( xn - xm )
                   IF ( curvature_solution_effects )  THEN
                      aero_mass(n) = aero_mass(n) + 0.5_wp * weight(n) * ( xsm - xsn )
                      aero_mass(m) = aero_mass(m) + 0.5_wp * weight(m) * ( xsn - xsm )
                   ENDIF
                   weight(n) = weight(n) - 0.5_wp * weight(m)
                   weight(m) = weight(n)

                ENDIF

             ENDIF

          ENDDO

          ql_vp(k,j,i) = ql_vp(k,j,i) + mass(n) / factor_volume_to_mass

       ENDDO

       IF ( ANY(weight < 0.0_wp) )  THEN
             WRITE( message_string, * ) 'negative weighting factor'
             CALL message( 'lpm_droplet_collision', 'PA0028',      &
                            2, 2, -1, 6, 1 )
       ENDIF

       particles(1:number_of_particles)%radius = ( mass(1:number_of_particles) /   &
                                                   ( weight(1:number_of_particles) &
                                                     * factor_volume_to_mass       &
                                                   )                               &
                                                 )**0.33333333333333_wp

       IF ( curvature_solution_effects )  THEN
          particles(1:number_of_particles)%aux1 = ( aero_mass(1:number_of_particles) / &
                                                    ( weight(1:number_of_particles)    &
                                                      * 4.0_wp / 3.0_wp * pi * rho_s   &
                                                    )                                  &
                                                  )**0.33333333333333_wp
       ENDIF

       particles(1:number_of_particles)%weight_factor = weight(1:number_of_particles)

       DEALLOCATE( weight, mass )
       IF ( curvature_solution_effects )  DEALLOCATE( aero_mass )
       IF ( .NOT. use_kernel_tables )  DEALLOCATE( ckernel )

!
!--    Check if LWC is conserved during collision process
       IF ( ql_v(k,j,i) /= 0.0_wp )  THEN
          IF ( ql_vp(k,j,i) / ql_v(k,j,i) >= 1.0001_wp  .OR.                      &
               ql_vp(k,j,i) / ql_v(k,j,i) <= 0.9999_wp )  THEN
             WRITE( message_string, * ) ' LWC is not conserved during',           &
                                        ' collision! ',                           &
                                        ' LWC after condensation: ', ql_v(k,j,i), &
                                        ' LWC after collision: ', ql_vp(k,j,i)
             CALL message( 'lpm_droplet_collision', 'PA0040', 2, 2, -1, 6, 1 )
          ENDIF
       ENDIF

    ENDIF

    CALL cpu_log( log_point_s(43), 'lpm_droplet_coll', 'stop' )

 END SUBROUTINE lpm_droplet_collision
