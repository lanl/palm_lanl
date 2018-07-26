!> @file lpm_boundary_conds.f90
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
! $Id: lpm_boundary_conds.f90 3067 2018-06-12 14:04:34Z suehring $
! Remove write statements
! 
! 2801 2018-02-14 16:01:55Z thiele
! Introduce particle transfer in nested models.
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2701 2017-12-15 15:40:50Z suehring
! Changes from last commit documented
! 
! 2698 2017-12-14 18:46:24Z suehring
! Particle reflections at downward-facing walls implemented. Moreover, 
! reflections are adjusted to revised particle grid box location. 
! (responsible Philipp Thiele)
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2606 2017-11-10 10:36:31Z schwenkel
! Changed particle box locations: center of particle box now coincides 
! with scalar grid point of same index.
! Renamed module and subroutines: lpm_pack_arrays_mod -> lpm_pack_and_sort_mod
! lpm_pack_all_arrays -> lpm_sort_in_subboxes, lpm_pack_arrays -> lpm_pack
! lpm_sort -> lpm_sort_timeloop_done
! 
! 2318 2017-07-20 17:27:44Z suehring
! Get topography top index via Function call 
! 
! 2317 2017-07-20 17:27:19Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography and surface concept
! Rename character range into location, as range is an intrinsic. 
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1929 2016-06-09 16:25:25Z suehring
! Rewritten wall reflection 
!
! 1822 2016-04-07 07:49:42Z hoffmann
! Tails removed. Unused variables removed.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1359 2014-04-11 17:15:14Z hoffmann
! New particle structure integrated. 
! Kind definition added to all floating point numbers.
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 849 2012-03-15 10:35:09Z raasch
! routine renamed lpm_boundary_conds, bottom and top boundary conditions
! included (former part of advec_particles)
!
! 824 2012-02-17 09:09:57Z raasch
! particle attributes speed_x|y|z_sgs renamed rvar1|2|3
!
! Initial version (2007/03/09)
!
! Description:
! ------------
!> Boundary conditions for the Lagrangian particles.
!> The routine consists of two different parts. One handles the bottom (flat)
!> and top boundary. In this part, also particles which exceeded their lifetime
!> are deleted.
!> The other part handles the reflection of particles from vertical walls.
!> This part was developed by Jin Zhang during 2006-2007.
!>
!> To do: Code structure for finding the t_index values and for checking the
!> -----  reflection conditions is basically the same for all four cases, so it
!>        should be possible to further simplify/shorten it.
!>
!> THE WALLS PART OF THIS ROUTINE HAS NOT BEEN TESTED FOR OCEAN RUNS SO FAR!!!!
!> (see offset_ocean_*)
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_boundary_conds( location , i, j, k )
 

    USE arrays_3d,                                                             &
        ONLY:  zu, zw

    USE control_parameters,                                                    &
        ONLY:  dz, message_string, particle_maximum_age

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE grid_variables,                                                        &
        ONLY:  ddx, dx, ddy, dy

    USE indices,                                                               &
        ONLY:  nxl, nxr, nyn, nys, nz, nzb, wall_flags_0,nyng,nysg

    USE kinds

    USE particle_attributes,                                                   &
        ONLY:  deleted_particles, ibc_par_b, ibc_par_t, number_of_particles,   &
               particles, particle_type, offset_ocean_nzt_m1,                  &
               use_sgs_for_particles

    USE pegrid

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  location     !<
    
    INTEGER(iwp), INTENT(IN) ::  i !< 
    INTEGER(iwp), INTENT(IN) ::  j !<
    INTEGER(iwp), INTENT(IN) ::  k !<
    
    INTEGER(iwp) ::  inc            !< dummy for sorting algorithmus
    INTEGER(iwp) ::  ir             !< dummy for sorting algorithmus
    INTEGER(iwp) ::  i1             !< grid index (x) of old particle position
    INTEGER(iwp) ::  i2             !< grid index (x) of current particle position
    INTEGER(iwp) ::  i3             !< grid index (x) of intermediate particle position
    INTEGER(iwp) ::  jr             !< dummy for sorting algorithmus
    INTEGER(iwp) ::  j1             !< grid index (y) of old particle position
    INTEGER(iwp) ::  j2             !< grid index (y) of current particle position
    INTEGER(iwp) ::  j3             !< grid index (y) of intermediate particle position
    INTEGER(iwp) ::  k1             !< grid index (z) of old particle position
    INTEGER(iwp) ::  k2             !< grid index (z) of current particle position
    INTEGER(iwp) ::  k3             !< grid index (z) of intermediate particle position
    INTEGER(iwp) ::  n              !< particle number
    INTEGER(iwp) ::  t_index        !< running index for intermediate particle timesteps in reflection algorithmus
    INTEGER(iwp) ::  t_index_number !< number of intermediate particle timesteps in reflection algorithmus
    INTEGER(iwp) ::  tmp_x          !< dummy for sorting algorithm
    INTEGER(iwp) ::  tmp_y          !< dummy for sorting algorithm
    INTEGER(iwp) ::  tmp_z          !< dummy for sorting algorithm

    INTEGER(iwp), DIMENSION(0:10) :: x_ind(0:10) = 0 !< index array (x) of intermediate particle positions 
    INTEGER(iwp), DIMENSION(0:10) :: y_ind(0:10) = 0 !< index array (y) of intermediate particle positions
    INTEGER(iwp), DIMENSION(0:10) :: z_ind(0:10) = 0 !< index array (z) of intermediate particle positions
    
    LOGICAL  ::  cross_wall_x    !< flag to check if particle reflection along x is necessary
    LOGICAL  ::  cross_wall_y    !< flag to check if particle reflection along y is necessary
    LOGICAL  ::  cross_wall_z    !< flag to check if particle reflection along z is necessary
    LOGICAL  ::  reflect_x       !< flag to check if particle is already reflected along x
    LOGICAL  ::  reflect_y       !< flag to check if particle is already reflected along y
    LOGICAL  ::  reflect_z       !< flag to check if particle is already reflected along z
    LOGICAL  ::  tmp_reach_x     !< dummy for sorting algorithmus
    LOGICAL  ::  tmp_reach_y     !< dummy for sorting algorithmus
    LOGICAL  ::  tmp_reach_z     !< dummy for sorting algorithmus
    LOGICAL  ::  x_wall_reached  !< flag to check if particle has already reached wall
    LOGICAL  ::  y_wall_reached  !< flag to check if particle has already reached wall
    LOGICAL  ::  z_wall_reached  !< flag to check if particle has already reached wall

    LOGICAL, DIMENSION(0:10) ::  reach_x  !< flag to check if particle is at a yz-wall
    LOGICAL, DIMENSION(0:10) ::  reach_y  !< flag to check if particle is at a xz-wall
    LOGICAL, DIMENSION(0:10) ::  reach_z  !< flag to check if particle is at a xy-wall

    REAL(wp) ::  dt_particle    !< particle timestep
    REAL(wp) ::  dum            !< dummy argument
    REAL(wp) ::  eps = 1E-10_wp !< security number to check if particle has reached a wall
    REAL(wp) ::  pos_x          !< intermediate particle position (x)
    REAL(wp) ::  pos_x_old      !< particle position (x) at previous particle timestep
    REAL(wp) ::  pos_y          !< intermediate particle position (y)
    REAL(wp) ::  pos_y_old      !< particle position (y) at previous particle timestep
    REAL(wp) ::  pos_z          !< intermediate particle position (z)
    REAL(wp) ::  pos_z_old      !< particle position (z) at previous particle timestep
    REAL(wp) ::  prt_x          !< current particle position (x)
    REAL(wp) ::  prt_y          !< current particle position (y)
    REAL(wp) ::  prt_z          !< current particle position (z)
    REAL(wp) ::  t_old          !< previous reflection time
    REAL(wp) ::  tmp_t          !< dummy for sorting algorithmus
    REAL(wp) ::  xwall          !< location of wall in x
    REAL(wp) ::  ywall          !< location of wall in y
    REAL(wp) ::  zwall          !< location of wall in z

    REAL(wp), DIMENSION(0:10) ::  t  !< reflection time


    IF ( location == 'bottom/top' )  THEN

!
!--    Apply boundary conditions to those particles that have crossed the top or
!--    bottom boundary and delete those particles, which are older than allowed
       DO  n = 1, number_of_particles

!
!--       Stop if particles have moved further than the length of one 
!--       PE subdomain (newly released particles have age = age_m!)
          IF ( particles(n)%age /= particles(n)%age_m )  THEN
             IF ( ABS(particles(n)%speed_x) >                                  &
                  ((nxr-nxl+2)*dx)/(particles(n)%age-particles(n)%age_m)  .OR. &
                  ABS(particles(n)%speed_y) >                                  &
                  ((nyn-nys+2)*dy)/(particles(n)%age-particles(n)%age_m) )  THEN

                  WRITE( message_string, * )  'particle too fast.  n = ',  n 
                  CALL message( 'lpm_boundary_conds', 'PA0148', 2, 2, -1, 6, 1 )
             ENDIF
          ENDIF

          IF ( particles(n)%age > particle_maximum_age  .AND.  &
               particles(n)%particle_mask )                              &
          THEN
             particles(n)%particle_mask  = .FALSE.
             deleted_particles = deleted_particles + 1
          ENDIF

          IF ( particles(n)%z >= zw(nz)  .AND.  particles(n)%particle_mask )  THEN
             IF ( ibc_par_t == 1 )  THEN
!
!--             Particle absorption
                particles(n)%particle_mask  = .FALSE.
                deleted_particles = deleted_particles + 1
             ELSEIF ( ibc_par_t == 2 )  THEN
!
!--             Particle reflection
                particles(n)%z       = 2.0_wp * zw(nz) - particles(n)%z
                particles(n)%speed_z = -particles(n)%speed_z
                IF ( use_sgs_for_particles  .AND. &
                     particles(n)%rvar3 > 0.0_wp )  THEN
                   particles(n)%rvar3 = -particles(n)%rvar3
                ENDIF
             ENDIF
          ENDIF
          
          IF ( particles(n)%z < zw(0)  .AND.  particles(n)%particle_mask )  THEN
             IF ( ibc_par_b == 1 )  THEN
!
!--             Particle absorption
                particles(n)%particle_mask  = .FALSE.
                deleted_particles = deleted_particles + 1
             ELSEIF ( ibc_par_b == 2 )  THEN
!
!--             Particle reflection
                particles(n)%z       = 2.0_wp * zw(0) - particles(n)%z
                particles(n)%speed_z = -particles(n)%speed_z
                IF ( use_sgs_for_particles  .AND. &
                     particles(n)%rvar3 < 0.0_wp )  THEN
                   particles(n)%rvar3 = -particles(n)%rvar3
                ENDIF
             ENDIF
          ENDIF
       ENDDO

    ELSEIF ( location == 'walls' )  THEN


       CALL cpu_log( log_point_s(48), 'lpm_wall_reflect', 'start' )

       DO  n = 1, number_of_particles
!
!--       Recalculate particle timestep
          dt_particle = particles(n)%age - particles(n)%age_m
!
!--       Obtain x/y indices for current particle position
          i2 = particles(n)%x * ddx
          j2 = particles(n)%y * ddy
          IF (zw(k)   < particles(n)%z ) k2 = k + 1
          IF (zw(k-1) > particles(n)%z ) k2 = k - 1 
!
!--       Save current particle positions
          prt_x = particles(n)%x
          prt_y = particles(n)%y
          prt_z = particles(n)%z
!
!--       Recalculate old particle positions
          pos_x_old = particles(n)%x - particles(n)%speed_x * dt_particle
          pos_y_old = particles(n)%y - particles(n)%speed_y * dt_particle
          pos_z_old = particles(n)%z - particles(n)%speed_z * dt_particle
!
!--       Obtain x/y indices for old particle positions
          i1 = i
          j1 = j
          k1 = k
!
!--       Determine horizontal as well as vertical walls at which particle can 
!--       be potentially reflected. 
!--       Start with walls aligned in yz layer.
!--       Wall to the right 
          IF ( prt_x > pos_x_old )  THEN
             xwall = ( i1 + 1 ) * dx
!
!--       Wall to the left
          ELSE
             xwall = i1 * dx
          ENDIF
!
!--       Walls aligned in xz layer
!--       Wall to the north
          IF ( prt_y > pos_y_old )  THEN
             ywall = ( j1 +1 ) * dy
!--       Wall to the south
          ELSE
             ywall = j1 * dy
          ENDIF

          IF ( prt_z > pos_z_old ) THEN
             zwall = zw(k)
          ELSE
             zwall = zw(k-1)
          ENDIF     
!
!--       Initialize flags to check if particle reflection is necessary
          cross_wall_x = .FALSE.
          cross_wall_y = .FALSE.
          cross_wall_z = .FALSE.
!
!--       Initialize flags to check if a wall is reached
          reach_x      = .FALSE.
          reach_y      = .FALSE.
          reach_z      = .FALSE.
!
!--       Initialize flags to check if a particle was already reflected
          reflect_x    = .FALSE.
          reflect_y    = .FALSE.
          reflect_z    = .FALSE.
!
!--       Initialize flags to check if a wall is already crossed.
!--       ( Required to obtain correct indices. )
          x_wall_reached = .FALSE.
          y_wall_reached = .FALSE.
          z_wall_reached = .FALSE.
!
!--       Initialize time array 
          t     = 0.0_wp
!
!--       Check if particle can reach any wall. This case, calculate the 
!--       fractional time needed to reach this wall. Store this fractional
!--       timestep in array t. Moreover, store indices for these grid
!--       boxes where the respective wall belongs to.  
!--       Start with x-direction.
          t_index    = 1
          t(t_index) = ( xwall - pos_x_old )                                   &
                     / MERGE( MAX( prt_x - pos_x_old,  1E-30_wp ),             &
                              MIN( prt_x - pos_x_old, -1E-30_wp ),             &
                              prt_x > pos_x_old )
          x_ind(t_index)   = i2
          y_ind(t_index)   = j1
          z_ind(t_index)   = k1
          reach_x(t_index) = .TRUE.
          reach_y(t_index) = .FALSE.
          reach_z(t_index) = .FALSE.
!
!--       Store these values only if particle really reaches any wall. t must 
!--       be in a interval between [0:1].
          IF ( t(t_index) <= 1.0_wp .AND. t(t_index) >= 0.0_wp )  THEN
             t_index      = t_index + 1
             cross_wall_x = .TRUE.
          ENDIF
!
!--       y-direction
          t(t_index) = ( ywall - pos_y_old )                                   &
                     / MERGE( MAX( prt_y - pos_y_old,  1E-30_wp ),             &
                              MIN( prt_y - pos_y_old, -1E-30_wp ),             &
                              prt_y > pos_y_old )
          x_ind(t_index)   = i1
          y_ind(t_index)   = j2
          z_ind(t_index)   = k1
          reach_x(t_index) = .FALSE.
          reach_y(t_index) = .TRUE.
          reach_z(t_index) = .FALSE.
          IF ( t(t_index) <= 1.0_wp .AND. t(t_index) >= 0.0_wp )  THEN
             t_index      = t_index + 1
             cross_wall_y = .TRUE.
          ENDIF
!
!--       z-direction
          t(t_index) = (zwall - pos_z_old )                                    &
                     / MERGE( MAX( prt_z - pos_z_old,  1E-30_wp ),             &
                              MIN( prt_z - pos_z_old, -1E-30_wp ),             &
                              prt_z > pos_z_old )
                     
          x_ind(t_index)   = i1
          y_ind(t_index)   = j1
          z_ind(t_index)   = k2
          reach_x(t_index) = .FALSE.
          reach_y(t_index) = .FALSE.
          reach_z(t_index) = .TRUE.
          IF( t(t_index) <= 1.0_wp .AND. t(t_index) >= 0.0_wp) THEN
             t_index      = t_index + 1
             cross_wall_z = .TRUE.
          ENDIF
          
          t_index_number = t_index - 1
!
!--       Carry out reflection only if particle reaches any wall
          IF ( cross_wall_x .OR. cross_wall_y .OR. cross_wall_z )  THEN
!
!--          Sort fractional timesteps in ascending order. Also sort the 
!--          corresponding indices and flag according to the time interval a  
!--          particle reaches the respective wall.
             inc = 1
             jr  = 1
             DO WHILE ( inc <= t_index_number )
                inc = 3 * inc + 1
             ENDDO

             DO WHILE ( inc > 1 )
                inc = inc / 3
                DO  ir = inc+1, t_index_number
                   tmp_t       = t(ir)
                   tmp_x       = x_ind(ir)
                   tmp_y       = y_ind(ir)
                   tmp_z       = z_ind(ir)
                   tmp_reach_x = reach_x(ir)
                   tmp_reach_y = reach_y(ir)
                   tmp_reach_z = reach_z(ir)
                   jr    = ir
                   DO WHILE ( t(jr-inc) > tmp_t )
                      t(jr)       = t(jr-inc)
                      x_ind(jr)   = x_ind(jr-inc)
                      y_ind(jr)   = y_ind(jr-inc)
                      z_ind(jr)   = z_ind(jr-inc)
                      reach_x(jr) = reach_x(jr-inc)
                      reach_y(jr) = reach_y(jr-inc)
                      reach_z(jr) = reach_z(jr-inc)
                      jr    = jr - inc
                      IF ( jr <= inc )  EXIT
                   ENDDO
                   t(jr)       = tmp_t
                   x_ind(jr)   = tmp_x
                   y_ind(jr)   = tmp_y
                   z_ind(jr)   = tmp_z
                   reach_x(jr) = tmp_reach_x
                   reach_y(jr) = tmp_reach_y
                   reach_z(jr) = tmp_reach_z
                ENDDO
             ENDDO
!
!--          Initialize temporary particle positions
             pos_x = pos_x_old
             pos_y = pos_y_old
             pos_z = pos_z_old
!
!--          Loop over all times a particle possibly moves into a new grid box
             t_old = 0.0_wp
             DO t_index = 1, t_index_number 
!           
!--             Calculate intermediate particle position according to the 
!--             timesteps a particle reaches any wall.
                pos_x = pos_x + ( t(t_index) - t_old ) * dt_particle           &
                                                       * particles(n)%speed_x
                pos_y = pos_y + ( t(t_index) - t_old ) * dt_particle           &
                                                       * particles(n)%speed_y
                pos_z = pos_z + ( t(t_index) - t_old ) * dt_particle           &
                                                       * particles(n)%speed_z
!
!--             Obtain x/y grid indices for intermediate particle position from
!--             sorted index array
                i3 = x_ind(t_index)
                j3 = y_ind(t_index)
                k3 = z_ind(t_index)
!
!--             Check which wall is already reached
                IF ( .NOT. x_wall_reached )  x_wall_reached = reach_x(t_index) 
                IF ( .NOT. y_wall_reached )  y_wall_reached = reach_y(t_index)
                IF ( .NOT. z_wall_reached )  z_wall_reached = reach_z(t_index)
!
!--             Check if a particle needs to be reflected at any yz-wall. If 
!--             necessary, carry out reflection. Please note, a security 
!--             constant is required, as the particle position does not 
!--             necessarily exactly match the wall location due to rounding 
!--             errors.
                IF ( reach_x(t_index)                      .AND.               & 
                     ABS( pos_x - xwall ) < eps            .AND.               &
                     .NOT. BTEST(wall_flags_0(k3,j3,i3),0) .AND.               &
                     .NOT. reflect_x )  THEN
! 
! 
!--                Reflection in x-direction. 
!--                Ensure correct reflection by MIN/MAX functions, depending on
!--                direction of particle transport. 
!--                Due to rounding errors pos_x does not exactly match the wall
!--                location, leading to erroneous reflection.             
                   pos_x = MERGE( MIN( 2.0_wp * xwall - pos_x, xwall ),        &
                                  MAX( 2.0_wp * xwall - pos_x, xwall ),        &
                                  particles(n)%x > xwall )
!
!--                Change sign of particle speed                     
                   particles(n)%speed_x = - particles(n)%speed_x
!
!--                Also change sign of subgrid-scale particle speed
                   particles(n)%rvar1 = - particles(n)%rvar1
!
!--                Set flag that reflection along x is already done
                   reflect_x          = .TRUE.
!
!--                As the particle does not cross any further yz-wall during 
!--                this timestep, set further x-indices to the current one.
                   x_ind(t_index:t_index_number) = i1
!
!--             If particle already reached the wall but was not reflected, 
!--             set further x-indices to the new one.
                ELSEIF ( x_wall_reached .AND. .NOT. reflect_x )  THEN
                    x_ind(t_index:t_index_number) = i2
                ENDIF !particle reflection in x direction done

!
!--             Check if a particle needs to be reflected at any xz-wall. If 
!--             necessary, carry out reflection. Please note, a security 
!--             constant is required, as the particle position does not 
!--             necessarily exactly match the wall location due to rounding 
!--             errors. 
                IF ( reach_y(t_index)                      .AND.               & 
                     ABS( pos_y - ywall ) < eps            .AND.               &
                     .NOT. BTEST(wall_flags_0(k3,j3,i3),0) .AND.               &
                     .NOT. reflect_y )  THEN
! 
! 
!--                Reflection in y-direction. 
!--                Ensure correct reflection by MIN/MAX functions, depending on
!--                direction of particle transport. 
!--                Due to rounding errors pos_y does not exactly match the wall
!--                location, leading to erroneous reflection.             
                   pos_y = MERGE( MIN( 2.0_wp * ywall - pos_y, ywall ),        &
                                  MAX( 2.0_wp * ywall - pos_y, ywall ),        &
                                  particles(n)%y > ywall )
!
!--                Change sign of particle speed                     
                   particles(n)%speed_y = - particles(n)%speed_y
!
!--                Also change sign of subgrid-scale particle speed
                   particles(n)%rvar2 = - particles(n)%rvar2
!
!--                Set flag that reflection along y is already done
                   reflect_y          = .TRUE.
!
!--                As the particle does not cross any further xz-wall during 
!--                this timestep, set further y-indices to the current one.
                   y_ind(t_index:t_index_number) = j1
!
!--             If particle already reached the wall but was not reflected, 
!--             set further y-indices to the new one.
                ELSEIF ( y_wall_reached .AND. .NOT. reflect_y )  THEN
                    y_ind(t_index:t_index_number) = j2
                ENDIF !particle reflection in y direction done
                
!
!--             Check if a particle needs to be reflected at any xy-wall. If 
!--             necessary, carry out reflection. Please note, a security 
!--             constant is required, as the particle position does not 
!--             necessarily exactly match the wall location due to rounding 
!--             errors. 
                IF ( reach_z(t_index)                      .AND.               & 
                     ABS( pos_z - zwall ) < eps            .AND.               &
                     .NOT. BTEST(wall_flags_0(k3,j3,i3),0) .AND.               &
                     .NOT. reflect_z )  THEN
! 
! 
!--                Reflection in z-direction. 
!--                Ensure correct reflection by MIN/MAX functions, depending on
!--                direction of particle transport. 
!--                Due to rounding errors pos_z does not exactly match the wall
!--                location, leading to erroneous reflection.             
                   pos_z = MERGE( MIN( 2.0_wp * zwall - pos_z, zwall ),        &
                                  MAX( 2.0_wp * zwall - pos_z, zwall ),        &
                                  particles(n)%z > zwall )
!
!--                Change sign of particle speed                     
                   particles(n)%speed_z = - particles(n)%speed_z
!
!--                Also change sign of subgrid-scale particle speed
                   particles(n)%rvar3 = - particles(n)%rvar3
!
!--                Set flag that reflection along z is already done
                   reflect_z          = .TRUE.
!
!--                As the particle does not cross any further xy-wall during 
!--                this timestep, set further z-indices to the current one.
                   z_ind(t_index:t_index_number) = k1
!
!--             If particle already reached the wall but was not reflected, 
!--             set further z-indices to the new one.
                ELSEIF ( z_wall_reached .AND. .NOT. reflect_z )  THEN
                    z_ind(t_index:t_index_number) = k2
                ENDIF !particle reflection in z direction done                
                
!
!--             Swap time
                t_old = t(t_index)

             ENDDO
!
!--          If a particle was reflected, calculate final position from last
!--          intermediate position.
             IF ( reflect_x .OR. reflect_y .OR. reflect_z )  THEN

                particles(n)%x = pos_x + ( 1.0_wp - t_old ) * dt_particle      &
                                                         * particles(n)%speed_x
                particles(n)%y = pos_y + ( 1.0_wp - t_old ) * dt_particle      &
                                                         * particles(n)%speed_y
                particles(n)%z = pos_z + ( 1.0_wp - t_old ) * dt_particle      &
                                                         * particles(n)%speed_z

             ENDIF

          ENDIF

       ENDDO

       CALL cpu_log( log_point_s(48), 'lpm_wall_reflect', 'stop' )

    ENDIF

 END SUBROUTINE lpm_boundary_conds
