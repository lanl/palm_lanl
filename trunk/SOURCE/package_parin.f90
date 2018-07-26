!> @file package_parin.f90
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
! $Id: package_parin.f90 3049 2018-05-29 13:52:36Z Giersch $
! Error messages revised
! 
! 2932 2018-03-26 09:39:22Z Giersch
! renamed particles_par to particle_parameters
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2375 2017-08-29 14:10:28Z schwenkel
! Added aerosol_species
!
! 2312 2017-07-14 20:26:51Z hoffmann
! Aerosol initialization improved.
!
! 2263 2017-06-08 14:59:01Z schwenkel
! Implemented splitting and merging algorithm
!
! 2183 2017-03-17 14:29:15Z schwenkel
!
! 2182 2017-03-17 14:27:40Z schwenkel
! Added parameters for simplified particle initialization.
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1936 2016-06-13 13:37:44Z suehring
! +deallocate_memory, step_dealloc
!
! 1871 2016-04-15 11:46:09Z hoffmann
! Initialization of aerosols added.
!
! 1833 2016-04-07 14:23:03Z raasch
! reading of spectra_par moved to spectra_mod
!
! 1831 2016-04-07 13:15:51Z hoffmann
! curvature_solution_effects added
!
! 1826 2016-04-07 12:01:39Z maronga
! Reading of &radiation_par moved to radiation_model_mod.
! Reading of &canopy_par moved to plant_canopy_model_mod.
!
! 822 2016-04-07 07:49:42Z hoffmann
! +collision_algorithm
! Tails removed.
!
! 1817 2016-04-06 15:44:20Z maronga
! Reading of &lsm_par moved to land_surface_model_mod.
!
! 1788 2016-03-10 11:01:04Z maronga
! Parameter dewfall removed.
!
! 1786 2016-03-08 05:49:27Z raasch
! cpp-direktives for spectra removed
!
! 1757 2016-02-22 15:49:32Z maronga
! Added parameter unscheduled_radiation_calls
!
! 1691 2015-10-26 16:17:44Z maronga
! Added skip_time_do_lsm, skip_time_do_radiation, and emissivity
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1585 2015-04-30 07:05:52Z maronga
! Added several radiation_par parameters
!
! 1575 2015-03-27 09:56:27Z raasch
! +seed_follows_topography in particles_par
!
! 1553 2015-03-03 17:33:54Z maronga
! Resorting of lsm_par
!
! 1551 2015-03-03 14:18:16Z maronga
! Several changes in the liste for land surface model and radiation model
!
! 1496 2014-12-02 17:25:50Z maronga
! Added support for the land surface model and radiation scheme
!
! 1484 2014-10-21 10:53:05Z kanani
! Changes due to new module structure of the plant canopy model:
!   module plant_canopy_model_mod added,
!   new package/namelist canopy_par added, i.e. the canopy model is no longer
!   steered over the inipar-namelist,
!   drag_coefficient, leaf_surface_concentration and scalar_exchange_coefficient
!   renamed to canopy_drag_coeff, leaf_surface_conc and leaf_scalar_exch_coeff.
! Changed statement tags in CONTINUE-statement
!
! 1367 2014-04-23 15:18:30Z witha
! Bugfix: module kinds must be used
!
! 1359 2014-04-11 17:15:14Z hoffmann
! +alloc_factor, + min_nr_particle
! -dt_sort_particles, -maximum_number_of_particles
!
! 1340 2014-03-25 19:45:13Z kanani
! REAL constants defined as wp-kinds
!
! 1324 2014-03-21 09:13:16Z suehring
! Bugfix: Missing variable dt_data_output output added to ONLY statement
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
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 828 2012-02-21 12:00:36Z raasch
! +dissipation_classes, radius_classes in parpar
!
! 825 2012-02-19 03:03:44Z raasch
! wang_collision_kernel and turbulence_effects_on_collision in particles_par
! replaced by collision_kernel
!
! Revision 1.1  2000/12/28 13:21:57  raasch
! Initial revision
!
!
! Description:
! ------------
!> This subroutine reads from the NAMELIST file variables controling model
!> software packages which are used optionally in the run.
!>
!> @todo Perform all actions in the respective submodules and remove
!>       package_parin
!------------------------------------------------------------------------------!
 SUBROUTINE package_parin


    USE control_parameters,                                                    &
        ONLY:  dt_data_output, dt_dopts, dt_dvrp, message_string,              &
               particle_maximum_age, threshold

    USE dvrp_variables,                                                        &
        ONLY:  clip_dvrp_l, clip_dvrp_n, clip_dvrp_r, clip_dvrp_s,             &
               cluster_size, color_interval, dvrpsize_interval,                &
               dvrp_directory, dvrp_file, dvrp_host, dvrp_output,              &
               dvrp_password, dvrp_username, groundplate_color,                &
               isosurface_color, mode_dvrp, particle_color,                    &
               particle_dvrpsize, pathlines_fadeintime,                        &
               pathlines_fadeouttime, pathlines_linecount,                     &
               pathlines_maxhistory, pathlines_wavecount,                      &
               pathlines_wavetime, slicer_range_limits_dvrp, superelevation,   &
               superelevation_x, superelevation_y, topography_color,           &
               vc_alpha, vc_gradient_normals, vc_mode, vc_size_x, vc_size_y,   &
               vc_size_z

    USE kinds

    USE particle_attributes,                                                   &
        ONLY:  aero_species, aero_type, aero_weight, alloc_factor, bc_par_b,   &
               bc_par_lr, bc_par_ns, bc_par_t, collision_kernel,               &
               curvature_solution_effects, deallocate_memory, density_ratio,   &
               dissipation_classes, dt_min_part, dt_prel,                      &
               dt_write_particle_data, end_time_prel, initial_weighting_factor,&
               log_sigma, max_number_particles_per_gridbox,                    &
               merging, min_nr_particle, na,                                   &
               number_concentration, number_of_particle_groups,                &
               number_particles_per_gridbox, particles_per_point,              &
               particle_advection, particle_advection_start, pdx, pdy, pdz,    &
               psb, psl, psn, psr, pss, pst, radius, radius_classes,           &
               radius_merge, radius_split, random_start_position,              &
               read_particles_from_restartfile, rm,                            &
               seed_follows_topography, splitting, splitting_factor,           &
               splitting_factor_max, splitting_function, splitting_mode,       &
               step_dealloc, use_sgs_for_particles,                            &
               vertical_particle_advection, weight_factor_merge,               &
               weight_factor_split, write_particle_statistics

    IMPLICIT NONE

    CHARACTER (LEN=80) ::  line  !<



    NAMELIST /dvrp_graphics_par/  clip_dvrp_l, clip_dvrp_n, clip_dvrp_r,       &
                                  clip_dvrp_s, cluster_size, color_interval,   &
                                  dt_dvrp, dvrpsize_interval, dvrp_directory,  &
                                  dvrp_file, dvrp_host, dvrp_output,           &
                                  dvrp_password, dvrp_username,                &
                                  groundplate_color, isosurface_color,         &
                                  mode_dvrp, particle_color, particle_dvrpsize,&
                                  pathlines_fadeintime, pathlines_fadeouttime, &
                                  pathlines_linecount, pathlines_maxhistory,   &
                                  pathlines_wavecount, pathlines_wavetime,     &
                                  slicer_range_limits_dvrp, superelevation,    &
                                  superelevation_x, superelevation_y,          &
                                  threshold, topography_color, vc_alpha,       &
                                  vc_gradient_normals, vc_mode, vc_size_x,     &
                                  vc_size_y, vc_size_z

    NAMELIST /particles_par/      aero_species, aero_type, aero_weight,        &
                                  alloc_factor, bc_par_b, bc_par_lr,           &
                                  bc_par_ns, bc_par_t,                         &
                                  collision_kernel, curvature_solution_effects,&
                                  deallocate_memory, density_ratio,            &
                                  dissipation_classes, dt_dopts,               &
                                  dt_min_part, dt_prel,                        &
                                  dt_write_particle_data,                      &
                                  end_time_prel, initial_weighting_factor,     &
                                  log_sigma,                                   &
                                  max_number_particles_per_gridbox, merging,   &
                                  min_nr_particle,                             &
                                  na, number_concentration,                    &
                                  number_of_particle_groups,                   &
                                  number_particles_per_gridbox,                &
                                  particles_per_point,                         &
                                  particle_advection_start,                    &
                                  particle_maximum_age, pdx, pdy, pdz, psb,    &
                                  psl, psn, psr, pss, pst, radius,             &
                                  radius_classes, radius_merge, radius_split,  &
                                  random_start_position,                       &
                                  read_particles_from_restartfile, rm,         &
                                  seed_follows_topography,                     &
                                  splitting, splitting_factor,                 &
                                  splitting_factor_max, splitting_function,    &
                                  splitting_mode, step_dealloc,                &
                                  use_sgs_for_particles,                       &
                                  vertical_particle_advection,                 &
                                  weight_factor_merge, weight_factor_split,    &
                                  write_particle_statistics

                                  
    NAMELIST /particle_parameters/                                             &
                                  aero_species, aero_type, aero_weight,        &
                                  alloc_factor, bc_par_b, bc_par_lr,           &
                                  bc_par_ns, bc_par_t,                         &
                                  collision_kernel, curvature_solution_effects,&
                                  deallocate_memory, density_ratio,            &
                                  dissipation_classes, dt_dopts,               &
                                  dt_min_part, dt_prel,                        &
                                  dt_write_particle_data,                      &
                                  end_time_prel, initial_weighting_factor,     &
                                  log_sigma,                                   &
                                  max_number_particles_per_gridbox, merging,   &
                                  min_nr_particle,                             &
                                  na, number_concentration,                    &
                                  number_of_particle_groups,                   &
                                  number_particles_per_gridbox,                &
                                  particles_per_point,                         &
                                  particle_advection_start,                    &
                                  particle_maximum_age, pdx, pdy, pdz, psb,    &
                                  psl, psn, psr, pss, pst, radius,             &
                                  radius_classes, radius_merge, radius_split,  &
                                  random_start_position,                       &
                                  read_particles_from_restartfile, rm,         &
                                  seed_follows_topography,                     &
                                  splitting, splitting_factor,                 &
                                  splitting_factor_max, splitting_function,    &
                                  splitting_mode, step_dealloc,                &
                                  use_sgs_for_particles,                       &
                                  vertical_particle_advection,                 &
                                  weight_factor_merge, weight_factor_split,    &
                                  write_particle_statistics
!
!-- Position the namelist-file at the beginning (it was already opened in
!-- parin), search for the namelist-group of the package and position the
!-- file at this line. Do the same for each optionally used package.
    line = ' '


#if defined( __dvrp_graphics )
    REWIND ( 11 )
    line = ' '
    DO   WHILE ( INDEX( line, '&dvrp_graphics_par' ) == 0 )
       READ ( 11, '(A)', END=20 )  line
    ENDDO
    BACKSPACE ( 11 )

!
!-- Read user-defined namelist
    READ ( 11, dvrp_graphics_par )

 20 CONTINUE
#endif

!
!-- Try to find particles package
    REWIND ( 11 )
    line = ' '
    DO   WHILE ( INDEX( line, '&particle_parameters' ) == 0 )
       READ ( 11, '(A)', END=30 )  line
    ENDDO
    BACKSPACE ( 11 )

!
!-- Read user-defined namelist
    READ ( 11, particle_parameters )

!
!-- Set flag that indicates that particles are switched on
    particle_advection = .TRUE.
    
    GOTO 31

!
!-- Try to find particles package (old namelist)
30  REWIND ( 11 )
    line = ' '
    DO   WHILE ( INDEX( line, '&particles_par' ) == 0 )
       READ ( 11, '(A)', END=31 )  line
    ENDDO
    BACKSPACE ( 11 )

!
!-- Read user-defined namelist
    READ ( 11, particles_par )
    
    
    message_string = 'namelist particles_par is deprecated and will be ' //    &
                     'removed in near future. Please use namelist ' //         &
                     'particle_parameters instead'
    CALL message( 'package_parin', 'PA0487', 0, 1, 0, 6, 0 )

!
!-- Set flag that indicates that particles are switched on
    particle_advection = .TRUE.

 31 CONTINUE

 END SUBROUTINE package_parin
