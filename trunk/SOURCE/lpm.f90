MODULE lpm_mod

!> @file lpm.f90
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
! $Id: lpm.f90 2801 2018-02-14 16:01:55Z thiele $
! Changed lpm from subroutine to module.
! Introduce particle transfer in nested models.
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2701 2017-12-15 15:40:50Z suehring
! Changes from last commit documented
! 
! 2698 2017-12-14 18:46:24Z suehring
! Grid indices passed to lpm_boundary_conds. (responsible Philipp Thiele)
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
! 2418 2017-09-06 15:24:24Z suehring
! Major bugfixes in modeling SGS particle speeds (since revision 1359).
! Particle sorting added to distinguish between already completed and 
! non-completed particles.
! 
! 2263 2017-06-08 14:59:01Z schwenkel
! Implemented splitting and merging algorithm
! 
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography concept
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1936 2016-06-13 13:37:44Z suehring
! Call routine for deallocation of unused memory.
! Formatting adjustments
! 
! 1929 2016-06-09 16:25:25Z suehring
! Call wall boundary conditions only if particles are in the vertical range of 
! topography.
!
! 1822 2016-04-07 07:49:42Z hoffmann
! Tails removed.
!
! Initialization of sgs model not necessary for the use of cloud_droplets and
! use_sgs_for_particles.
!
! lpm_release_set integrated.
!
! Unused variabled removed.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1416 2014-06-04 16:04:03Z suehring
! user_lpm_advec is called for each gridpoint.
! Bugfix: in order to prevent an infinite loop, time_loop_done is set .TRUE. 
! at the head of the do-loop.  
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
! 1318 2014-03-17 13:35:16Z raasch
! module interfaces removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 851 2012-03-15 14:32:58Z raasch
! Bugfix: resetting of particle_mask and tail mask moved from routine
! lpm_exchange_horiz to here (end of sub-timestep loop)
!
! 849 2012-03-15 10:35:09Z raasch
! original routine advec_particles split into several subroutines and renamed
! lpm
!
! 831 2012-02-22 00:29:39Z raasch
! thermal_conductivity_l and diff_coeff_l now depend on temperature and
! pressure
!
! 828 2012-02-21 12:00:36Z raasch
! fast hall/wang kernels with fixed radius/dissipation classes added,
! particle feature color renamed class, routine colker renamed
! recalculate_kernel,
! lower limit for droplet radius changed from 1E-7 to 1E-8
!
! Bugfix: transformation factor for dissipation changed from 1E5 to 1E4
!
! 825 2012-02-19 03:03:44Z raasch
! droplet growth by condensation may include curvature and solution effects,
! initialisation of temporary particle array for resorting removed,
! particle attributes speed_x|y|z_sgs renamed rvar1|2|3,
! module wang_kernel_mod renamed lpm_collision_kernels_mod,
! wang_collision_kernel renamed wang_kernel
!
!
! Revision 1.1  1999/11/25 16:16:06  raasch
! Initial revision
!
!
! Description:
! ------------
!> Particle advection
!------------------------------------------------------------------------------!
 

    USE arrays_3d,                                                             &
        ONLY:  ql_c, ql_v, ql_vp

    USE control_parameters,                                                    &
        ONLY:  cloud_droplets, dt_3d, dt_3d_reached, dt_3d_reached_l,          &
               molecular_viscosity, simulated_time, topography

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s

    USE indices,                                                               &
        ONLY: nxl, nxr, nys, nyn, nzb, nzb_max, nzt

    USE kinds

    USE lpm_exchange_horiz_mod,                                                &
        ONLY:  dealloc_particles_array, lpm_exchange_horiz, lpm_move_particle

    USE lpm_init_mod,                                                          &
        ONLY: lpm_create_particle, PHASE_RELEASE

    USE lpm_pack_and_sort_mod

    USE particle_attributes,                                                   &
        ONLY:  collision_kernel, deleted_particles, deallocate_memory,         &
               dt_write_particle_data, dt_prel, end_time_prel,                 &
               grid_particles, merging,  number_of_particles,                  &
               number_of_particle_groups, particles, particle_groups,          &
               prt_count, splitting, step_dealloc, time_prel,                  &
               time_write_particle_data, trlp_count_sum, trlp_count_recv_sum,  &
               trnp_count_sum, trnp_count_recv_sum, trrp_count_sum,            & 
               trrp_count_recv_sum, trsp_count_sum, trsp_count_recv_sum,       &
               use_sgs_for_particles, write_particle_statistics
                              
    USE pegrid

 USE pmc_particle_interface,                                                &
     ONLY: pmcp_c_get_particle_from_parent, pmcp_p_fill_particle_win,       &
           pmcp_c_send_particle_to_parent, pmcp_p_empty_particle_win,       &
           pmcp_p_delete_particles_in_fine_grid_area

    USE pmc_interface,                                                      &
       ONLY: nested_run

 IMPLICIT NONE
 PRIVATE
 SAVE

 INTERFACE lpm
    MODULE PROCEDURE lpm
 END INTERFACE lpm

 PUBLIC lpm

CONTAINS
 SUBROUTINE lpm
    IMPLICIT NONE

    INTEGER(iwp)       ::  i                  !<
    INTEGER(iwp)       ::  ie                 !<
    INTEGER(iwp)       ::  is                 !<
    INTEGER(iwp)       ::  j                  !<
    INTEGER(iwp)       ::  je                 !<
    INTEGER(iwp)       ::  js                 !<
    INTEGER(iwp), SAVE ::  lpm_count = 0      !<
    INTEGER(iwp)       ::  k                  !<
    INTEGER(iwp)       ::  ke                 !<
    INTEGER(iwp)       ::  ks                 !<
    INTEGER(iwp)       ::  m                  !<
    INTEGER(iwp), SAVE ::  steps = 0          !<

    LOGICAL            ::  first_loop_stride  !<

    CALL cpu_log( log_point(25), 'lpm', 'start' )

!
!-- Write particle data at current time on file.
!-- This has to be done here, before particles are further processed,
!-- because they may be deleted within this timestep (in case that
!-- dt_write_particle_data = dt_prel = particle_maximum_age).
    time_write_particle_data = time_write_particle_data + dt_3d
    IF ( time_write_particle_data >= dt_write_particle_data )  THEN

       CALL lpm_data_output_particles
!
!--    The MOD function allows for changes in the output interval with restart
!--    runs.
       time_write_particle_data = MOD( time_write_particle_data, &
                                  MAX( dt_write_particle_data, dt_3d ) )
    ENDIF

!
!-- Initialize arrays for marking those particles to be deleted after the
!-- (sub-) timestep
    deleted_particles = 0

!
!-- Initialize variables used for accumulating the number of particles
!-- exchanged between the subdomains during all sub-timesteps (if sgs
!-- velocities are included). These data are output further below on the
!-- particle statistics file.
    trlp_count_sum      = 0
    trlp_count_recv_sum = 0
    trrp_count_sum      = 0
    trrp_count_recv_sum = 0
    trsp_count_sum      = 0
    trsp_count_recv_sum = 0
    trnp_count_sum      = 0
    trnp_count_recv_sum = 0


!
!-- Calculate exponential term used in case of particle inertia for each
!-- of the particle groups
    DO  m = 1, number_of_particle_groups
       IF ( particle_groups(m)%density_ratio /= 0.0_wp )  THEN
          particle_groups(m)%exp_arg  =                                        &
                    4.5_wp * particle_groups(m)%density_ratio *                &
                    molecular_viscosity / ( particle_groups(m)%radius )**2

          particle_groups(m)%exp_term = EXP( -particle_groups(m)%exp_arg *     &
                    dt_3d )
       ENDIF
    ENDDO

!
!-- If necessary, release new set of particles
    IF ( time_prel >= dt_prel  .AND.  end_time_prel > simulated_time )  THEN

       CALL lpm_create_particle(PHASE_RELEASE)
!
!--    The MOD function allows for changes in the output interval with
!--    restart runs.
       time_prel = MOD( time_prel, MAX( dt_prel, dt_3d ) )

    ENDIF
!
!-- Reset summation arrays
    IF ( cloud_droplets)  THEN
       ql_c  = 0.0_wp
       ql_v  = 0.0_wp
       ql_vp = 0.0_wp
    ENDIF

    first_loop_stride = .TRUE.
    grid_particles(:,:,:)%time_loop_done = .TRUE.
!
!-- Timestep loop for particle advection.
!-- This loop has to be repeated until the advection time of every particle
!-- (within the total domain!) has reached the LES timestep (dt_3d).
!-- In case of including the SGS velocities, the particle timestep may be
!-- smaller than the LES timestep (because of the Lagrangian timescale 
!-- restriction) and particles may require to undergo several particle 
!-- timesteps, before the LES timestep is reached. Because the number of these 
!-- particle timesteps to be carried out is unknown at first, these steps are 
!-- carried out in the following infinite loop with exit condition.
    DO
       CALL cpu_log( log_point_s(44), 'lpm_advec', 'start' )
       CALL cpu_log( log_point_s(44), 'lpm_advec', 'pause' )
       
!
!--    If particle advection includes SGS velocity components, calculate the
!--    required SGS quantities (i.e. gradients of the TKE, as well as 
!--    horizontally averaged profiles of the SGS TKE and the resolved-scale 
!--    velocity variances)
       IF ( use_sgs_for_particles  .AND.  .NOT. cloud_droplets )  THEN
          CALL lpm_init_sgs_tke
       ENDIF

!
!--    In case SGS-particle speed is considered, particles may carry out 
!--    several particle timesteps. In order to prevent unnecessary 
!--    treatment of particles that already reached the final time level, 
!--    particles are sorted into contiguous blocks of finished and 
!--    not-finished particles, in addition to their already sorting
!--    according to their sub-boxes. 
       IF ( .NOT. first_loop_stride  .AND.  use_sgs_for_particles )            &
          CALL lpm_sort_timeloop_done

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt

                number_of_particles = prt_count(k,j,i)
!
!--             If grid cell gets empty, flag must be true
                IF ( number_of_particles <= 0 )  THEN
                   grid_particles(k,j,i)%time_loop_done = .TRUE.
                   CYCLE
                ENDIF

                IF ( .NOT. first_loop_stride  .AND.  &
                     grid_particles(k,j,i)%time_loop_done ) CYCLE

                particles => grid_particles(k,j,i)%particles(1:number_of_particles)

                particles(1:number_of_particles)%particle_mask = .TRUE.
!
!--             Initialize the variable storing the total time that a particle 
!--             has advanced within the timestep procedure
                IF ( first_loop_stride )  THEN
                   particles(1:number_of_particles)%dt_sum = 0.0_wp
                ENDIF
!
!--             Particle (droplet) growth by condensation/evaporation and 
!--             collision
                IF ( cloud_droplets  .AND.  first_loop_stride)  THEN
!
!--                Droplet growth by condensation / evaporation
                   CALL lpm_droplet_condensation(i,j,k)
!
!--                Particle growth by collision
                   IF ( collision_kernel /= 'none' )  THEN
                      CALL lpm_droplet_collision(i,j,k)
                   ENDIF

                ENDIF
!
!--             Initialize the switch used for the loop exit condition checked
!--             at the end of this loop. If at least one particle has failed to
!--             reach the LES timestep, this switch will be set false in 
!--             lpm_advec.
                dt_3d_reached_l = .TRUE.

!
!--             Particle advection
                CALL lpm_advec(i,j,k)
!
!--             Particle reflection from walls. Only applied if the particles 
!--             are in the vertical range of the topography. (Here, some
!--             optimization is still possible.)
                IF ( topography /= 'flat' .AND. k < nzb_max + 2 )  THEN
                   CALL  lpm_boundary_conds( 'walls', i, j, k )
                ENDIF
!
!--             User-defined actions after the calculation of the new particle 
!--             position
                CALL user_lpm_advec(i,j,k)
!
!--             Apply boundary conditions to those particles that have crossed 
!--             the top or bottom boundary and delete those particles, which are
!--             older than allowed
                CALL lpm_boundary_conds( 'bottom/top', i, j, k )
!
!---            If not all particles of the actual grid cell have reached the 
!--             LES timestep, this cell has to do another loop iteration. Due to
!--             the fact that particles can move into neighboring grid cells, 
!--             these neighbor cells also have to perform another loop iteration.
!--             Please note, this realization does not work properly if 
!--             particles move into another subdomain. 
                IF ( .NOT. dt_3d_reached_l )  THEN
                   ks = MAX(nzb+1,k-1)
                   ke = MIN(nzt,k+1)
                   js = MAX(nys,j-1)
                   je = MIN(nyn,j+1)
                   is = MAX(nxl,i-1)
                   ie = MIN(nxr,i+1)
                   grid_particles(ks:ke,js:je,is:ie)%time_loop_done = .FALSE.
                ELSE
                   grid_particles(k,j,i)%time_loop_done = .TRUE.
                ENDIF

             ENDDO
          ENDDO
       ENDDO

       steps = steps + 1
       dt_3d_reached_l = ALL(grid_particles(:,:,:)%time_loop_done)
!
!--    Find out, if all particles on every PE have completed the LES timestep
!--    and set the switch corespondingly
#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( dt_3d_reached_l, dt_3d_reached, 1, MPI_LOGICAL, &
                           MPI_LAND, comm2d, ierr )
#else
       dt_3d_reached = dt_3d_reached_l
#endif

       CALL cpu_log( log_point_s(44), 'lpm_advec', 'stop' )

!
!--    Increment time since last release
       IF ( dt_3d_reached )  time_prel = time_prel + dt_3d
!
!--    Apply splitting and merging algorithm
       IF ( cloud_droplets )  THEN
          IF ( splitting ) THEN
             CALL lpm_splitting
          ENDIF
          IF ( merging ) THEN
             CALL lpm_merging
          ENDIF
       ENDIF
!
!--    Move Particles local to PE to a different grid cell
       CALL lpm_move_particle

!
!--    Horizontal boundary conditions including exchange between subdmains
       CALL lpm_exchange_horiz

       IF ( .NOT. dt_3d_reached .OR. .NOT. nested_run )   THEN    ! IF .FALSE., lpm_sort_in_subboxes is done inside pcmp
!
!--       Pack particles (eliminate those marked for deletion),
!--       determine new number of particles
          CALL lpm_sort_in_subboxes
!
!--       Initialize variables for the next (sub-) timestep, i.e., for marking
!--       those particles to be deleted after the timestep
          deleted_particles = 0
       ENDIF

       IF ( dt_3d_reached )  EXIT

       first_loop_stride = .FALSE.
    ENDDO   ! timestep loop
!    
!-- in case of nested runs do the transfer of particles after every full model time step
    IF ( nested_run )   THEN
       CALL particles_from_parent_to_child
       CALL particles_from_child_to_parent
       CALL pmcp_p_delete_particles_in_fine_grid_area

       CALL lpm_sort_in_subboxes

       deleted_particles = 0
    ENDIF

!
!-- Calculate the new liquid water content for each grid box
    IF ( cloud_droplets )  CALL lpm_calc_liquid_water_content
!
!-- Deallocate unused memory
    IF ( deallocate_memory  .AND.  lpm_count == step_dealloc )  THEN
       CALL dealloc_particles_array
       lpm_count = 0
    ELSEIF ( deallocate_memory )  THEN
       lpm_count = lpm_count + 1
    ENDIF

!
!-- Set particle attributes.
!-- Feature is not available if collision is activated, because the respective
!-- particle attribute (class) is then used for storing the particle radius
!-- class.
    IF ( collision_kernel == 'none' )  CALL lpm_set_attributes

!
!-- Set particle attributes defined by the user
    CALL user_lpm_set_attributes

!
!-- Write particle statistics (in particular the number of particles
!-- exchanged between the subdomains) on file
    IF ( write_particle_statistics )  CALL lpm_write_exchange_statistics

    CALL cpu_log( log_point(25), 'lpm', 'stop' )

 END SUBROUTINE lpm

 SUBROUTINE particles_from_parent_to_child
    IMPLICIT NONE

    CALL pmcp_c_get_particle_from_parent                         ! Child actions
    CALL pmcp_p_fill_particle_win                                ! Parent actions

    RETURN
 END SUBROUTINE particles_from_parent_to_child

 SUBROUTINE particles_from_child_to_parent
    IMPLICIT NONE

    CALL pmcp_c_send_particle_to_parent                         ! Child actions
    CALL pmcp_p_empty_particle_win                              ! Parent actions

    RETURN
 END SUBROUTINE particles_from_child_to_parent


END MODULE lpm_mod
