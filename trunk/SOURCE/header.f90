!> @file header.f90
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
! 2018-10-25 cbegeman
! Add dirichlet bottom boundary conditions for salinity
! 
! Former revisions:
! -----------------
! $Id: header.f90 3083 2018-06-19 14:03:12Z gronemeier $
! Print RANS-mode constants
! 
! 3065 2018-06-12 07:03:02Z Giersch
! Header output concerning stretching revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
! 
! 2967 2018-04-13 11:22:08Z raasch
! bugfix: missing parallel cpp-directives added
! 
! 2883 2018-03-14 08:29:10Z Giersch
! Format of the output of dt_dopr_listing (325) has been changed
! 
! 2817 2018-02-19 16:32:21Z knoop
! Preliminary gust module interface implemented
! 
! 2776 2018-01-31 10:44:42Z Giersch
! Variable synthetic_turbulence_generator has been abbreviated
! 
! 2746 2018-01-15 12:06:04Z suehring
! Move flag plant canopy to modules
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
! Print information about turbulence closure (TG)
! Print information about inifor initialization (MS)
! 
! 2575 2017-10-24 09:57:58Z maronga
! Added output for complex terrain simulations
! 
! 2544 2017-10-13 18:09:32Z maronga
! Moved initial day of year and time to inipar. 
! 
! 2339 2017-08-07 13:55:26Z gronemeier
! corrected timestamp in header
! 
! 2338 2017-08-07 12:15:38Z gronemeier
! Modularize 1D model
! 
! 2320 2017-07-21 12:47:43Z suehring
! Modularize large-scale forcing and nudging
! 
! 2300 2017-06-29 13:31:14Z raasch
! host-specific code removed
! 
! 2299 2017-06-29 10:14:38Z maronga
! Modified output for spinups
! 
! 2298 2017-06-29 09:28:18Z raasch
! MPI2 related parts removed
! 
! 2270 2017-06-09 12:18:47Z maronga
! Renamed Prandtl layer to constant flux layer
! 
! 2259 2017-06-08 09:09:11Z gronemeier
! Implemented synthetic turbulence generator
!
! 2258 2017-06-08 07:55:13Z suehring
! Bugfix, add pre-preprocessor directives to enable non-parrallel mode
! 
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography and surface concept
! Generic tunnel setup added
! 
! 2200 2017-04-11 11:37:51Z suehring
! monotonic_adjustment removed
! 
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC relatec code removed
! 
! 2073 2016-11-30 14:34:05Z raasch
! small bugfix concerning output of scalar profiles
! 
! 2050 2016-11-08 15:00:55Z gronemeier
! Implement turbulent outflow condition
! 
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1992 2016-08-12 15:14:59Z suehring
! Adapted for top_scalarflux
! 
! 1960 2016-07-12 16:34:24Z suehring
! Treat humidity and passive scalar separately.
! Modify misleading information concerning humidity.
! Bugfix, change unit for humidity flux.
! 
! 1957 2016-07-07 10:43:48Z suehring
! flight module added
!
! 1931 2016-06-10 12:06:59Z suehring
! Rename multigrid into multigrid_noopt
!
! 1902 2016-05-09 11:18:56Z suehring
! Write information about masking_method only for multigrid solver
!
! 1849 2016-04-08 11:33:18Z hoffmann
! Adapted for modularization of microphysics
!
! 1833 2016-04-07 14:23:03Z raasch
! spectrum renamed spectra_mod, output of spectra related quantities moved to
! spectra_mod
!
! 1831 2016-04-07 13:15:51Z hoffmann
! turbulence renamed collision_turbulence,
! drizzle renamed cloud_water_sedimentation
!
! 1826 2016-04-07 12:01:39Z maronga
! Moved radiation model header output to the respective module.
! Moved canopy model header output to the respective module.
!
! 1822 2016-04-07 07:49:42Z hoffmann
! Tails removed. icloud_scheme replaced by microphysics_*
!
! 1817 2016-04-06 15:44:20Z maronga
! Moved land_surface_model header output to the respective module.
!
! 1808 2016-04-05 19:44:00Z raasch
! routine local_flush replaced by FORTRAN statement
!
! 1797 2016-03-21 16:50:28Z raasch
! output of nesting datatransfer mode
!
! 1791 2016-03-11 10:41:25Z raasch
! output of nesting informations of all domains
!
! 1788 2016-03-10 11:01:04Z maronga
! Parameter dewfall removed
!
! 1786 2016-03-08 05:49:27Z raasch
! cpp-direktives for spectra removed
!
! 1783 2016-03-06 18:36:17Z raasch
! netcdf module and variable names changed, output of netcdf_deflate
!
! 1764 2016-02-28 12:45:19Z raasch
! output of nesting informations
!
! 1697 2015-10-28 17:14:10Z raasch
! small E- and F-FORMAT changes to avoid informative compiler messages about
! insufficient field width
!
! 1691 2015-10-26 16:17:44Z maronga
! Renamed prandtl_layer to constant_flux_layer, renames rif_min/rif_max to
! zeta_min/zeta_max.
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1675 2015-10-02 08:28:59Z gronemeier
! Bugfix: Definition of topography grid levels
! 
! 1660 2015-09-21 08:15:16Z gronemeier
! Bugfix: Definition of building/street canyon height if vertical grid stretching
!         starts below the maximum topography height.
! 
! 1590 2015-05-08 13:56:27Z maronga
! Bugfix: Added TRIM statements for character strings for LSM and radiation code
! 
! 1585 2015-04-30 07:05:52Z maronga
! Further output for radiation model(s).
!
! 1575 2015-03-27 09:56:27Z raasch
! adjustments for psolver-queries, output of seed_follows_topography
!
! 1560 2015-03-06 10:48:54Z keck
! output for recycling y shift
! 
! 1557 2015-03-05 16:43:04Z suehring
! output for monotonic limiter
!
! 1551 2015-03-03 14:18:16Z maronga
! Added informal output for land surface model and radiation model. Removed typo.
! 
! 1496 2014-12-02 17:25:50Z maronga
! Renamed: "radiation -> "cloud_top_radiation"
! 
! 1484 2014-10-21 10:53:05Z kanani
! Changes due to new module structure of the plant canopy model:
!   module plant_canopy_model_mod and output for new canopy model parameters 
!   (alpha_lad, beta_lad, lai_beta,...) added,
!   drag_coefficient, leaf_surface_concentration and scalar_exchange_coefficient
!   renamed to canopy_drag_coeff, leaf_surface_conc and leaf_scalar_exch_coeff,
!   learde renamed leaf_area_density.
! Bugfix: DO-WHILE-loop for lad header information additionally restricted 
! by maximum number of gradient levels (currently 10)
!
! 1482 2014-10-18 12:34:45Z raasch
! information about calculated or predefined virtual processor topology adjusted
!
! 1468 2014-09-24 14:06:57Z maronga
! Adapted for use on up to 6-digit processor cores
! 
! 1429 2014-07-15 12:53:45Z knoop
! header exended to provide ensemble_member_nr if specified
! 
! 1376 2014-04-26 11:21:22Z boeske
! Correction of typos
! 
! 1365 2014-04-22 15:03:56Z boeske
! New section 'Large scale forcing and nudging':
! output of large scale forcing and nudging information,
! new section for initial profiles created
! 
! 1359 2014-04-11 17:15:14Z hoffmann
! dt_sort_particles removed
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute 
! 
! 1327 2014-03-21 11:00:16Z raasch
! parts concerning iso2d and avs output removed,
! -netcdf output queries
!
! 1324 2014-03-21 09:13:16Z suehring
! Bugfix: module spectrum added
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL functions provided with KIND-attribute,
! some REAL constants defined as wp-kind
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements, 
! kinds are defined in new module kinds, 
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
! 
! 1308 2014-03-13 14:58:42Z fricke
! output of the fixed number of output time levels
! output_format adjusted for masked data if netcdf_data_format > 5
!
! 1299 2014-03-06 13:15:21Z heinze
! output for using large_scale subsidence in combination 
! with large_scale_forcing
! reformatting, more detailed explanations
!
! 1241 2013-10-30 11:36:58Z heinze
! output for nudging + large scale forcing from external file
!
! 1216 2013-08-26 09:31:42Z raasch
! output for transpose_compute_overlap
!
! 1212 2013-08-15 08:46:27Z raasch
! output for poisfft_hybrid removed
!
! 1179 2013-06-14 05:57:58Z raasch
! output of reference_state, use_reference renamed use_single_reference_value
!
! 1159 2013-05-21 11:58:22Z fricke
! +use_cmax
!
! 1115 2013-03-26 18:16:16Z hoffmann
! descriptions for Seifert-Beheng-cloud-physics-scheme added
!
! 1111 2013-03-08 23:54:10Z raasch
! output of accelerator board information
! ibc_p_b = 2 removed
!
! 1108 2013-03-05 07:03:32Z raasch
! bugfix for r1106
!
! 1106 2013-03-04 05:31:38Z raasch
! some format changes for coupled runs
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1031 2012-10-19 14:35:30Z raasch
! output of netCDF data format modified
!
! 1015 2012-09-27 09:23:24Z raasch
! output of Adjustment of mixing length to the Prandtl mixing length at first
! grid point above ground removed
!
! 1003 2012-09-14 14:35:53Z raasch
! output of information about equal/unequal subdomain size removed
!
! 1001 2012-09-13 14:08:46Z raasch
! all actions concerning leapfrog- and upstream-spline-scheme removed
!
! 978 2012-08-09 08:28:32Z fricke
! -km_damp_max, outflow_damping_width
! +pt_damping_factor, pt_damping_width
! +z0h
!
! 964 2012-07-26 09:14:24Z raasch
! output of profil-related quantities removed
!
! 940 2012-07-09 14:31:00Z raasch
! Output in case of simulations for pure neutral stratification (no pt-equation
! solved)
!
! 927 2012-06-06 19:15:04Z raasch
! output of masking_method for mg-solver
!
! 868 2012-03-28 12:21:07Z raasch
! translation velocity in Galilean transformation changed to 0.6 * ug
!
! 833 2012-02-22 08:55:55Z maronga
! Adjusted format for leaf area density
!
! 828 2012-02-21 12:00:36Z raasch
! output of dissipation_classes + radius_classes
!
! 825 2012-02-19 03:03:44Z raasch
! Output of cloud physics parameters/quantities complemented and restructured
!
! Revision 1.1  1997/08/11 06:17:20  raasch
! Initial revision
!
!
! Description:
! ------------
!> Writing a header with all important information about the current run.
!> This subroutine is called three times, two times at the beginning 
!> (writing information on files RUN_CONTROL and HEADER) and one time at the
!> end of the run, then writing additional information about CPU-usage on file
!> header.
!-----------------------------------------------------------------------------!
 SUBROUTINE header
 

    USE arrays_3d,                                                             &
        ONLY:  pt_init, q_init, s_init, sa_init, ug, vg, w_subs, zu, zw
        
    USE control_parameters
        
    USE cloud_parameters,                                                      &
        ONLY:  cp, l_v, r_d

    USE cpulog,                                                                &
        ONLY:  log_point_s
            
    USE date_and_time_mod,                                                     &
        ONLY:  day_of_year_init, time_utc_init

    USE dvrp_variables,                                                        &
        ONLY:  use_seperate_pe_for_dvrp_output
        
    USE flight_mod,                                                            &
        ONLY:  flight_header
        
    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE gust_mod,                                                              &
        ONLY: gust_header, gust_module_enabled
        
    USE indices,                                                               &
        ONLY:  mg_loc_ind, nnx, nny, nnz, nx, ny, nxl_mg, nxr_mg, nyn_mg,      &
               nys_mg, nzt, nzt_mg
        
    USE kinds
 
    USE land_surface_model_mod,                                                &
        ONLY: lsm_header

    USE lsf_nudging_mod,                                                       &
        ONLY:  lsf_nudging_header

    USE microphysics_mod,                                                      &
        ONLY:  cloud_water_sedimentation, collision_turbulence,                &
               c_sedimentation, limiter_sedimentation, nc_const,               &
               ventilation_effect

    USE model_1d_mod,                                                          &
        ONLY:  damp_level_ind_1d, dt_pr_1d, dt_run_control_1d, end_time_1d

    USE netcdf_interface,                                                      &
        ONLY:  netcdf_data_format, netcdf_data_format_string, netcdf_deflate

    USE particle_attributes,                                                   &
        ONLY:  bc_par_b, bc_par_lr, bc_par_ns, bc_par_t, collision_kernel,     &
               curvature_solution_effects,                                     &
               density_ratio, dissipation_classes, dt_min_part, dt_prel,       &
               dt_write_particle_data, end_time_prel,                          &
               number_of_particle_groups, particle_advection,                  &
               particle_advection_start,                                       &
               particles_per_point, pdx, pdy, pdz,  psb, psl, psn, psr, pss,   &
               pst, radius, radius_classes, random_start_position,             &
               seed_follows_topography,                                        &
               total_number_of_particles, use_sgs_for_particles,               &
               vertical_particle_advection, write_particle_statistics
        
    USE pegrid

    USE plant_canopy_model_mod,                                                &
        ONLY:  pcm_header

#if defined( __parallel )
    USE pmc_handle_communicator,                                               &
        ONLY:  pmc_get_model_info
#endif

    USE pmc_interface,                                                         &
        ONLY:  nested_run, nesting_datatransfer_mode, nesting_mode

    USE radiation_model_mod,                                                   &
        ONLY:  radiation, radiation_header
    
    USE spectra_mod,                                                           &
        ONLY:  calculate_spectra, spectra_header

    USE surface_mod,                                                           &
        ONLY:  surf_def_h, get_topography_top_index_ji

    USE synthetic_turbulence_generator_mod,                                    &
        ONLY:  stg_header

    USE turbulence_closure_mod,                                                &
        ONLY:  rans_const_c, rans_const_sigma

    IMPLICIT NONE

    CHARACTER (LEN=1)  ::  prec                !<
    
    CHARACTER (LEN=2)  ::  do2d_mode           !<
    
    CHARACTER (LEN=5)  ::  section_chr         !<
    
    CHARACTER (LEN=10) ::  coor_chr            !<
    CHARACTER (LEN=10) ::  host_chr            !<
    
    CHARACTER (LEN=16) ::  begin_chr           !<
    
    CHARACTER (LEN=26) ::  ver_rev             !<

    CHARACTER (LEN=32) ::  cpl_name            !<
    
    CHARACTER (LEN=40) ::  output_format       !<
    
    CHARACTER (LEN=70) ::  char1               !<
    CHARACTER (LEN=70) ::  char2               !<
    CHARACTER (LEN=70) ::  dopr_chr            !<
    CHARACTER (LEN=70) ::  do2d_xy             !<
    CHARACTER (LEN=70) ::  do2d_xz             !<
    CHARACTER (LEN=70) ::  do2d_yz             !<
    CHARACTER (LEN=70) ::  do3d_chr            !<
    CHARACTER (LEN=70) ::  domask_chr          !<
    CHARACTER (LEN=70) ::  run_classification  !<
    
    CHARACTER (LEN=85) ::  r_upper             !<
    CHARACTER (LEN=85) ::  r_lower             !<
    
    CHARACTER (LEN=86) ::  coordinates         !<
    CHARACTER (LEN=86) ::  gradients           !<
    CHARACTER (LEN=86) ::  slices              !<
    CHARACTER (LEN=86) ::  temperatures        !<
    CHARACTER (LEN=86) ::  ugcomponent         !<
    CHARACTER (LEN=86) ::  vgcomponent         !<

    CHARACTER (LEN=1), DIMENSION(1:3) ::  dir = (/ 'x', 'y', 'z' /)  !<

    INTEGER(iwp) ::  av             !<
    INTEGER(iwp) ::  bh             !<
    INTEGER(iwp) ::  blx            !<
    INTEGER(iwp) ::  bly            !<
    INTEGER(iwp) ::  bxl            !<
    INTEGER(iwp) ::  bxr            !<
    INTEGER(iwp) ::  byn            !<
    INTEGER(iwp) ::  bys            !<
    INTEGER(iwp) ::  ch             !<
    INTEGER(iwp) ::  count          !<
    INTEGER(iwp) ::  cpl_parent_id  !<
    INTEGER(iwp) ::  cwx            !<
    INTEGER(iwp) ::  cwy            !<
    INTEGER(iwp) ::  cxl            !<
    INTEGER(iwp) ::  cxr            !<
    INTEGER(iwp) ::  cyn            !<
    INTEGER(iwp) ::  cys            !<
    INTEGER(iwp) ::  dim            !<
    INTEGER(iwp) ::  i              !<
    INTEGER(iwp) ::  io             !<
    INTEGER(iwp) ::  j              !<
    INTEGER(iwp) ::  k              !<
    INTEGER(iwp) ::  l              !<
    INTEGER(iwp) ::  ll             !<
    INTEGER(iwp) ::  my_cpl_id      !<
    INTEGER(iwp) ::  n              !<
    INTEGER(iwp) ::  ncpl           !<
    INTEGER(iwp) ::  npe_total      !<
    

    REAL(wp) ::  cpuseconds_per_simulated_second  !<
    REAL(wp) ::  lower_left_coord_x               !< x-coordinate of nest domain
    REAL(wp) ::  lower_left_coord_y               !< y-coordinate of nest domain

!
!-- Open the output file. At the end of the simulation, output is directed
!-- to unit 19.
    IF ( ( runnr == 0 .OR. force_print_header )  .AND. &
         .NOT. simulated_time_at_begin /= simulated_time )  THEN
       io = 15   !  header output on file RUN_CONTROL
    ELSE
       io = 19   !  header output on file HEADER
    ENDIF
    CALL check_open( io )

!
!-- At the end of the run, output file (HEADER) will be rewritten with
!-- new information
    IF ( io == 19 .AND. simulated_time_at_begin /= simulated_time ) REWIND( 19 )

!
!-- Determine kind of model run
    IF ( TRIM( initializing_actions ) == 'read_restart_data' )  THEN
       run_classification = 'restart run'
    ELSEIF ( TRIM( initializing_actions ) == 'cyclic_fill' )  THEN
       run_classification = 'run with cyclic fill of 3D - prerun data'
    ELSEIF ( INDEX( initializing_actions, 'set_constant_profiles' ) /= 0 )  THEN
       run_classification = 'run without 1D - prerun'
    ELSEIF ( INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN
       run_classification = 'run with 1D - prerun'
    ELSEIF ( INDEX( initializing_actions, 'inifor' ) /= 0 )  THEN
       run_classification = 'run initialized with COSMO data'
    ELSEIF ( INDEX( initializing_actions, 'by_user' ) /=0 )  THEN
       run_classification = 'run initialized by user'
    ELSEIF ( INDEX( initializing_actions, 'initialize_vortex' ) /=0 )  THEN
       run_classification = 'run additionally initialized by a Rankine-vortex'
    ELSEIF ( INDEX( initializing_actions, 'initialize_ptanom' ) /=0 )  THEN
       run_classification = 'run additionally initialized by temperature anomaly'
    ELSEIF ( INDEX( initializing_actions, 'initialize_2D_bubble' ) /=0 )  THEN
       run_classification = 'run additionally initialized by 2D bubble'
    ELSEIF ( INDEX( initializing_actions, 'initialize_3D_bubble' ) /=0 )  THEN
       run_classification = 'run additionally initialized by 3D bubble'
    ELSE
       message_string = ' unknown action(s): ' // TRIM( initializing_actions )
       CALL message( 'header', 'PA0191', 0, 0, 0, 6, 0 )
    ENDIF
    IF ( nested_run )  run_classification = 'nested ' // run_classification
    IF ( ocean )  THEN
       run_classification = 'ocean - ' // run_classification
    ELSE
       run_classification = 'atmosphere - ' // run_classification
    ENDIF

!
!-- Run-identification, date, time, host
    host_chr = host(1:10)
    ver_rev = TRIM( version ) // '  ' // TRIM( revision )
    WRITE ( io, 100 )  ver_rev, TRIM( run_classification )
    IF ( TRIM( coupling_mode ) /= 'uncoupled' )  THEN
       WRITE ( io, 101 )  coupling_mode
    ENDIF
#if defined( __parallel )
    IF ( coupling_start_time /= 0.0_wp  .AND. .NOT. spinup )  THEN
       IF ( coupling_start_time > simulated_time_at_begin )  THEN
          WRITE ( io, 109 )
       ELSE
          WRITE ( io, 114 )
       ENDIF
    ENDIF
#endif
    IF ( ensemble_member_nr /= 0 )  THEN
       WRITE ( io, 512 )  run_date, run_identifier, run_time, runnr,           &
                       ADJUSTR( host_chr ), ensemble_member_nr
    ELSE
       WRITE ( io, 102 )  run_date, run_identifier, run_time, runnr,           &
                       ADJUSTR( host_chr )
    ENDIF
#if defined( __parallel )
    IF ( npex == -1  .AND.  npey == -1 )  THEN
       char1 = 'calculated'
    ELSE
       char1 = 'predefined'
    ENDIF
    IF ( threads_per_task == 1 )  THEN
       WRITE ( io, 103 )  numprocs, pdims(1), pdims(2), TRIM( char1 )
    ELSE
       WRITE ( io, 104 )  numprocs*threads_per_task, numprocs, &
                          threads_per_task, pdims(1), pdims(2), TRIM( char1 )
    ENDIF

    IF ( pdims(2) == 1 )  THEN
       WRITE ( io, 107 )  'x'
    ELSEIF ( pdims(1) == 1 )  THEN
       WRITE ( io, 107 )  'y'
    ENDIF
    IF ( use_seperate_pe_for_dvrp_output )  WRITE ( io, 105 )
    IF ( numprocs /= maximum_parallel_io_streams )  THEN
       WRITE ( io, 108 )  maximum_parallel_io_streams
    ENDIF
#endif

!
!-- Nesting informations
    IF ( nested_run )  THEN

#if defined( __parallel )
       WRITE ( io, 600 )  TRIM( nesting_mode ),                                &
                          TRIM( nesting_datatransfer_mode )
       CALL pmc_get_model_info( ncpl = ncpl, cpl_id = my_cpl_id )

       DO  n = 1, ncpl
          CALL pmc_get_model_info( request_for_cpl_id = n, cpl_name = cpl_name,&
                                   cpl_parent_id = cpl_parent_id,              &
                                   lower_left_x = lower_left_coord_x,          &
                                   lower_left_y = lower_left_coord_y,          &
                                   npe_total = npe_total )
          IF ( n == my_cpl_id )  THEN
             char1 = '*'
          ELSE
             char1 = ' '
          ENDIF
          WRITE ( io, 601 )  TRIM( char1 ), n, cpl_parent_id, npe_total,       &
                             lower_left_coord_x, lower_left_coord_y,           &
                             TRIM( cpl_name )
       ENDDO
#endif

    ENDIF
    WRITE ( io, 99 )

!
!-- Numerical schemes
    WRITE ( io, 110 )
    IF ( rans_mode )  THEN
       WRITE ( io, 124 )  TRIM( turbulence_closure ), 'RANS'
    ELSE
       WRITE ( io, 124 )  TRIM( turbulence_closure ), 'LES'
    ENDIF
    WRITE ( io, 121 )  TRIM( approximation )
    IF ( psolver(1:7) == 'poisfft' )  THEN
       WRITE ( io, 111 )  TRIM( fft_method )
       IF ( transpose_compute_overlap )  WRITE( io, 115 )
    ELSEIF ( psolver == 'sor' )  THEN
       WRITE ( io, 112 )  nsor_ini, nsor, omega_sor
    ELSEIF ( psolver(1:9) == 'multigrid' )  THEN
       WRITE ( io, 135 )  TRIM(psolver), cycle_mg, maximum_grid_level, ngsrb
       IF ( mg_cycles == -1 )  THEN
          WRITE ( io, 140 )  residual_limit
       ELSE
          WRITE ( io, 141 )  mg_cycles
       ENDIF
       IF ( mg_switch_to_pe0_level == 0 )  THEN
          WRITE ( io, 136 )  nxr_mg(1)-nxl_mg(1)+1, nyn_mg(1)-nys_mg(1)+1, &
                             nzt_mg(1)
       ELSEIF (  mg_switch_to_pe0_level /= -1 )  THEN
          WRITE ( io, 137 )  mg_switch_to_pe0_level,            &
                             mg_loc_ind(2,0)-mg_loc_ind(1,0)+1, &
                             mg_loc_ind(4,0)-mg_loc_ind(3,0)+1, &
                             nzt_mg(mg_switch_to_pe0_level),    &
                             nxr_mg(1)-nxl_mg(1)+1, nyn_mg(1)-nys_mg(1)+1, &
                             nzt_mg(1)
       ENDIF
       IF ( psolver == 'multigrid_noopt' .AND. masking_method )  WRITE ( io, 144 )
    ENDIF
    IF ( call_psolver_at_all_substeps  .AND. timestep_scheme(1:5) == 'runge' ) &
    THEN
       WRITE ( io, 142 )
    ENDIF

    IF ( momentum_advec == 'pw-scheme' )  THEN
       WRITE ( io, 113 )
    ELSEIF (momentum_advec == 'ws-scheme' )  THEN
       WRITE ( io, 503 )
    ENDIF
    IF ( scalar_advec == 'pw-scheme' )  THEN
       WRITE ( io, 116 )
    ELSEIF ( scalar_advec == 'ws-scheme' )  THEN
       WRITE ( io, 504 )
    ELSE
       WRITE ( io, 118 )
    ENDIF

    WRITE ( io, 139 )  TRIM( loop_optimization )

    IF ( galilei_transformation )  THEN
       IF ( use_ug_for_galilei_tr )  THEN
          char1 = '0.6 * geostrophic wind'
       ELSE
          char1 = 'mean wind in model domain'
       ENDIF
       IF ( simulated_time_at_begin == simulated_time )  THEN
          char2 = 'at the start of the run'
       ELSE
          char2 = 'at the end of the run'
       ENDIF
       WRITE ( io, 119 )  TRIM( char1 ), TRIM( char2 ),                        &
                          advected_distance_x/1000.0_wp,                       &
                          advected_distance_y/1000.0_wp
    ENDIF
    WRITE ( io, 122 )  timestep_scheme
    IF ( use_upstream_for_tke )  WRITE ( io, 143 )
    IF ( rayleigh_damping_factor /= 0.0_wp )  THEN
       IF ( .NOT. ocean )  THEN
          WRITE ( io, 123 )  'above', rayleigh_damping_height, &
               rayleigh_damping_factor
       ELSE
          WRITE ( io, 123 )  'below', rayleigh_damping_height, &
               rayleigh_damping_factor
       ENDIF
    ENDIF
    IF ( neutral )  WRITE ( io, 131 )  pt_surface
    IF ( humidity )  THEN
       IF ( .NOT. cloud_physics )  THEN
          WRITE ( io, 129 )
       ELSE
          WRITE ( io, 130 )
       ENDIF
    ENDIF
    IF ( passive_scalar )  WRITE ( io, 134 )
    IF ( conserve_volume_flow )  THEN
       WRITE ( io, 150 )  conserve_volume_flow_mode
       IF ( TRIM( conserve_volume_flow_mode ) == 'bulk_velocity' )  THEN
          WRITE ( io, 151 )  u_bulk, v_bulk
       ENDIF
    ELSEIF ( dp_external )  THEN
       IF ( dp_smooth )  THEN
          WRITE ( io, 152 )  dpdxy, dp_level_b, ', vertically smoothed.'
       ELSE
          WRITE ( io, 152 )  dpdxy, dp_level_b, '.'
       ENDIF
    ENDIF
    WRITE ( io, 99 )

!
!-- Runtime and timestep information
    WRITE ( io, 200 )
    IF ( .NOT. dt_fixed )  THEN
       WRITE ( io, 201 )  dt_max, cfl_factor
    ELSE
       WRITE ( io, 202 )  dt
    ENDIF
    WRITE ( io, 203 )  simulated_time_at_begin, end_time

    IF ( time_restart /= 9999999.9_wp  .AND. &
         simulated_time_at_begin == simulated_time )  THEN
       IF ( dt_restart == 9999999.9_wp )  THEN
          WRITE ( io, 204 )  ' Restart at:       ',time_restart
       ELSE
          WRITE ( io, 205 )  ' Restart at:       ',time_restart, dt_restart
       ENDIF
    ENDIF

    IF ( simulated_time_at_begin /= simulated_time )  THEN
       i = MAX ( log_point_s(10)%counts, 1 )
       IF ( ( simulated_time - simulated_time_at_begin ) == 0.0_wp )  THEN
          cpuseconds_per_simulated_second = 0.0_wp
       ELSE
          cpuseconds_per_simulated_second = log_point_s(10)%sum / &
                                            ( simulated_time -    &
                                              simulated_time_at_begin )
       ENDIF
       WRITE ( io, 206 )  simulated_time, log_point_s(10)%sum,      &
                          log_point_s(10)%sum / REAL( i, KIND=wp ), &
                          cpuseconds_per_simulated_second
       IF ( time_restart /= 9999999.9_wp  .AND.  time_restart < end_time )  THEN
          IF ( dt_restart == 9999999.9_wp )  THEN
             WRITE ( io, 204 )  ' Next restart at:     ',time_restart
          ELSE
             WRITE ( io, 205 )  ' Next restart at:     ',time_restart, dt_restart
          ENDIF
       ENDIF
    ENDIF


!
!-- Start time for coupled runs, if independent precursor runs for atmosphere
!-- and ocean are used or have been used. In this case, coupling_start_time
!-- defines the time when the coupling is switched on.
    IF ( coupling_start_time /= 0.0_wp )  THEN
       WRITE ( io, 207 )  coupling_start_time
    ENDIF

!
!-- Computational grid
    IF ( .NOT. ocean )  THEN
       WRITE ( io, 250 )  dx, dy
       
       DO i = 1, number_stretch_level_start+1
          WRITE ( io, 253 )  i, dz(i)
       ENDDO
       
       WRITE( io, 251 ) (nx+1)*dx, (ny+1)*dy, zu(nzt+1)
       
       IF ( ANY( dz_stretch_level_start_index < nzt+1 ) )  THEN
          WRITE( io, '(A)', advance='no') ' Vertical stretching starts at height:'
          DO i = 1, number_stretch_level_start
             WRITE ( io, '(F10.1,A3)', advance='no' )  dz_stretch_level_start(i), ' m,'
          ENDDO
          WRITE( io, '(/,A)', advance='no') ' Vertical stretching starts at index: '
          DO i = 1, number_stretch_level_start
             WRITE ( io, '(I12,A1)', advance='no' )  dz_stretch_level_start_index(i), ','
          ENDDO
          WRITE( io, '(/,A)', advance='no') ' Vertical stretching ends at height:  '
          DO i = 1, number_stretch_level_start
             WRITE ( io, '(F10.1,A3)', advance='no' )  dz_stretch_level_end(i), ' m,'
          ENDDO
          WRITE( io, '(/,A)', advance='no') ' Vertical stretching ends at index:   '
          DO i = 1, number_stretch_level_start
             WRITE ( io, '(I12,A1)', advance='no' )  dz_stretch_level_end_index(i), ','
          ENDDO
          WRITE( io, '(/,A)', advance='no') ' Factor used for stretching:          '
          DO i = 1, number_stretch_level_start
             WRITE ( io, '(F12.3,A1)', advance='no' )  dz_stretch_factor_array(i), ','
          ENDDO
       ENDIF
       
    ELSE
       WRITE ( io, 250 )  dx, dy
       DO i = 1, number_stretch_level_start+1
          WRITE ( io, 253 )  i, dz(i)
       ENDDO
       
       WRITE ( io, 251 ) (nx+1)*dx, (ny+1)*dy, zu(0)
       
       IF ( ANY( dz_stretch_level_start_index > 0 ) )  THEN
          WRITE( io, '(A)', advance='no') ' Vertical stretching starts at height:'
          DO i = 1, number_stretch_level_start
             WRITE ( io, '(F10.1,A3)', advance='no' )  dz_stretch_level_start(i), ' m,'
          ENDDO
          WRITE( io, '(/,A)', advance='no') ' Vertical stretching starts at index: '
          DO i = 1, number_stretch_level_start
             WRITE ( io, '(I12,A1)', advance='no' )  dz_stretch_level_start_index(i), ','
          ENDDO
          WRITE( io, '(/,A)', advance='no') ' Vertical stretching ends at height:  '
          DO i = 1, number_stretch_level_start
             WRITE ( io, '(F10.1,A3)', advance='no' )  dz_stretch_level_end(i), ' m,'
          ENDDO
          WRITE( io, '(/,A)', advance='no') ' Vertical stretching ends at index:   '
          DO i = 1, number_stretch_level_start
             WRITE ( io, '(I12,A1)', advance='no' )  dz_stretch_level_end_index(i), ','
          ENDDO
          WRITE( io, '(/,A)', advance='no') ' Factor used for stretching:          '
          DO i = 1, number_stretch_level_start
             WRITE ( io, '(F12.3,A1)', advance='no' )  dz_stretch_factor_array(i), ','
          ENDDO
       ENDIF
    ENDIF
    WRITE ( io, 254 )  nx, ny, nzt+1, MIN( nnx, nx+1 ), MIN( nny, ny+1 ),      &
                       MIN( nnz+2, nzt+2 )
    IF ( sloping_surface )  WRITE ( io, 260 )  alpha_surface

!
!-- Large scale forcing and nudging
    IF ( large_scale_forcing )  CALL lsf_nudging_header( io )

!
!-- Profile for the large scale vertial velocity
!-- Building output strings, starting with surface value
    IF ( large_scale_subsidence )  THEN
       temperatures = '   0.0'
       gradients = '------'
       slices = '     0'
       coordinates = '   0.0'
       i = 1
       DO  WHILE ( subs_vertical_gradient_level_i(i) /= -9999 )

          WRITE (coor_chr,'(E10.2,7X)')  &
                                w_subs(subs_vertical_gradient_level_i(i))
          temperatures = TRIM( temperatures ) // ' ' // TRIM( coor_chr )

          WRITE (coor_chr,'(E10.2,7X)')  subs_vertical_gradient(i)
          gradients = TRIM( gradients ) // ' ' // TRIM( coor_chr )

          WRITE (coor_chr,'(I10,7X)')  subs_vertical_gradient_level_i(i)
          slices = TRIM( slices ) // ' ' // TRIM( coor_chr )

          WRITE (coor_chr,'(F10.2,7X)')  subs_vertical_gradient_level(i)
          coordinates = TRIM( coordinates ) // ' '  // TRIM( coor_chr )

          IF ( i == 10 )  THEN
             EXIT
          ELSE
             i = i + 1
          ENDIF

       ENDDO

 
       IF ( .NOT. large_scale_forcing )  THEN
          WRITE ( io, 426 )  TRIM( coordinates ), TRIM( temperatures ), &
                             TRIM( gradients ), TRIM( slices )
       ENDIF


    ENDIF

!-- Profile of the geostrophic wind (component ug)
!-- Building output strings
    WRITE ( ugcomponent, '(F6.2)' )  ug_surface
    gradients = '------'
    slices = '     0'
    coordinates = '   0.0'
    i = 1
    DO  WHILE ( ug_vertical_gradient_level_ind(i) /= -9999 )
      
       WRITE (coor_chr,'(F6.2,1X)')  ug(ug_vertical_gradient_level_ind(i))
       ugcomponent = TRIM( ugcomponent ) // '  ' // TRIM( coor_chr )

       WRITE (coor_chr,'(F6.2,1X)')  ug_vertical_gradient(i)
       gradients = TRIM( gradients ) // '  ' // TRIM( coor_chr )

       WRITE (coor_chr,'(I6,1X)')  ug_vertical_gradient_level_ind(i)
       slices = TRIM( slices ) // '  ' // TRIM( coor_chr )

       WRITE (coor_chr,'(F6.1,1X)')  ug_vertical_gradient_level(i)
       coordinates = TRIM( coordinates ) // '  ' // TRIM( coor_chr )

       IF ( i == 10 )  THEN
          EXIT
       ELSE
          i = i + 1
       ENDIF

    ENDDO

    IF ( .NOT. large_scale_forcing )  THEN
       WRITE ( io, 423 )  TRIM( coordinates ), TRIM( ugcomponent ), &
                          TRIM( gradients ), TRIM( slices )
    ENDIF

!-- Profile of the geostrophic wind (component vg)
!-- Building output strings
    WRITE ( vgcomponent, '(F6.2)' )  vg_surface
    gradients = '------'
    slices = '     0'
    coordinates = '   0.0'
    i = 1
    DO  WHILE ( vg_vertical_gradient_level_ind(i) /= -9999 )

       WRITE (coor_chr,'(F6.2,1X)')  vg(vg_vertical_gradient_level_ind(i))
       vgcomponent = TRIM( vgcomponent ) // '  ' // TRIM( coor_chr )

       WRITE (coor_chr,'(F6.2,1X)')  vg_vertical_gradient(i)
       gradients = TRIM( gradients ) // '  ' // TRIM( coor_chr )

       WRITE (coor_chr,'(I6,1X)')  vg_vertical_gradient_level_ind(i)
       slices = TRIM( slices ) // '  ' // TRIM( coor_chr )

       WRITE (coor_chr,'(F6.1,1X)')  vg_vertical_gradient_level(i)
       coordinates = TRIM( coordinates ) // '  ' // TRIM( coor_chr )

       IF ( i == 10 )  THEN
          EXIT
       ELSE
          i = i + 1
       ENDIF
 
    ENDDO

    IF ( .NOT. large_scale_forcing )  THEN
       WRITE ( io, 424 )  TRIM( coordinates ), TRIM( vgcomponent ), &
                          TRIM( gradients ), TRIM( slices )
    ENDIF

!
!-- Topography
    WRITE ( io, 270 )  topography
    SELECT CASE ( TRIM( topography ) )

       CASE ( 'flat' )
          ! no actions necessary

       CASE ( 'single_building' )
          blx = INT( building_length_x / dx )
          bly = INT( building_length_y / dy )
          bh  = MINLOC( ABS( zw - building_height ), 1 ) - 1
          IF ( ABS( zw(bh  ) - building_height ) == &
               ABS( zw(bh+1) - building_height )    )  bh = bh + 1

          IF ( building_wall_left == 9999999.9_wp )  THEN
             building_wall_left = ( nx + 1 - blx ) / 2 * dx
          ENDIF
          bxl = INT ( building_wall_left / dx + 0.5_wp )
          bxr = bxl + blx

          IF ( building_wall_south == 9999999.9_wp )  THEN
             building_wall_south = ( ny + 1 - bly ) / 2 * dy
          ENDIF
          bys = INT ( building_wall_south / dy + 0.5_wp )
          byn = bys + bly

          WRITE ( io, 271 )  building_length_x, building_length_y, &
                             building_height, bxl, bxr, bys, byn

       CASE ( 'single_street_canyon' )
          ch  = MINLOC( ABS( zw - canyon_height ), 1 ) - 1
          IF ( ABS( zw(ch  ) - canyon_height ) == &
               ABS( zw(ch+1) - canyon_height )    )  ch = ch + 1
          IF ( canyon_width_x /= 9999999.9_wp )  THEN
!
!--          Street canyon in y direction
             cwx = NINT( canyon_width_x / dx )
             IF ( canyon_wall_left == 9999999.9_wp )  THEN
                canyon_wall_left = ( nx + 1 - cwx ) / 2 * dx
             ENDIF
             cxl = NINT( canyon_wall_left / dx )
             cxr = cxl + cwx
             WRITE ( io, 272 )  'y', canyon_height, ch, 'u', cxl, cxr

          ELSEIF ( canyon_width_y /= 9999999.9_wp )  THEN
!
!--          Street canyon in x direction
             cwy = NINT( canyon_width_y / dy )
             IF ( canyon_wall_south == 9999999.9_wp )  THEN
                canyon_wall_south = ( ny + 1 - cwy ) / 2 * dy
             ENDIF
             cys = NINT( canyon_wall_south / dy )
             cyn = cys + cwy
             WRITE ( io, 272 )  'x', canyon_height, ch, 'v', cys, cyn
          ENDIF

       CASE ( 'tunnel' )
          IF ( tunnel_width_x /= 9999999.9_wp )  THEN
!
!--          Tunnel axis in y direction
             IF ( tunnel_length == 9999999.9_wp  .OR.                          &
                  tunnel_length >= ( nx + 1 ) * dx )  THEN
                WRITE ( io, 273 )  'y', tunnel_height, tunnel_wall_depth,      &
                                        tunnel_width_x
             ELSE
                WRITE ( io, 274 )  'y', tunnel_height, tunnel_wall_depth,      &
                                        tunnel_width_x, tunnel_length
             ENDIF

          ELSEIF ( tunnel_width_y /= 9999999.9_wp )  THEN
!
!--          Tunnel axis in x direction
             IF ( tunnel_length == 9999999.9_wp  .OR.                          &
                  tunnel_length >= ( ny + 1 ) * dy )  THEN
                WRITE ( io, 273 )  'x', tunnel_height, tunnel_wall_depth,      &
                                        tunnel_width_y
             ELSE
                WRITE ( io, 274 )  'x', tunnel_height, tunnel_wall_depth,      &
                                        tunnel_width_y, tunnel_length
             ENDIF
          ENDIF

    END SELECT

    IF ( TRIM( topography ) /= 'flat' )  THEN
       IF ( TRIM( topography_grid_convention ) == ' ' )  THEN
          IF ( TRIM( topography ) == 'single_building' .OR.  &
               TRIM( topography ) == 'single_street_canyon' )  THEN
             WRITE ( io, 278 )
          ELSEIF ( TRIM( topography ) == 'read_from_file' )  THEN
             WRITE ( io, 279 )
          ENDIF
       ELSEIF ( TRIM( topography_grid_convention ) == 'cell_edge' )  THEN
          WRITE ( io, 278 )
       ELSEIF ( TRIM( topography_grid_convention ) == 'cell_center' )  THEN
          WRITE ( io, 279 )
       ENDIF
    ENDIF

!-- Complex terrain
    IF ( complex_terrain )  THEN
       WRITE( io, 280 ) 
       IF ( turbulent_inflow )  THEN
          WRITE( io, 281 )  zu( get_topography_top_index_ji( 0, 0, 's' ) )
       ENDIF
       IF ( TRIM( initializing_actions ) == 'cyclic_fill' )  THEN
          WRITE( io, 282 )
       ENDIF
    ENDIF

    IF ( syn_turb_gen )  CALL stg_header ( io )

    IF ( plant_canopy )  CALL pcm_header ( io )

    IF ( land_surface )  CALL lsm_header ( io )

    IF ( radiation )  CALL radiation_header ( io )

    IF ( gust_module_enabled )  CALL gust_header ( io )

!
!-- Boundary conditions
    IF ( ibc_p_b == 0 )  THEN
       r_lower = 'p(0)     = 0      |'
    ELSEIF ( ibc_p_b == 1 )  THEN
       r_lower = 'p(0)     = p(1)   |'
    ENDIF
    IF ( ibc_p_t == 0 )  THEN
       r_upper  = 'p(nzt+1) = 0      |'
    ELSE
       r_upper  = 'p(nzt+1) = p(nzt) |'
    ENDIF

    IF ( ibc_uv_b == 0 )  THEN
       r_lower = TRIM( r_lower ) // ' uv(0)     = -uv(1)                |'
    ELSE
       r_lower = TRIM( r_lower ) // ' uv(0)     = uv(1)                 |'
    ENDIF
    IF ( TRIM( bc_uv_t ) == 'dirichlet_0' )  THEN
       r_upper  = TRIM( r_upper  ) // ' uv(nzt+1) = 0                     |'
    ELSEIF ( ibc_uv_t == 0 )  THEN
       r_upper  = TRIM( r_upper  ) // ' uv(nzt+1) = ug(nzt+1), vg(nzt+1)  |'
    ELSE
       r_upper  = TRIM( r_upper  ) // ' uv(nzt+1) = uv(nzt)               |'
    ENDIF

    IF ( ibc_pt_b == 0 )  THEN
       IF ( land_surface )  THEN
          r_lower = TRIM( r_lower ) // ' pt(0)     = from soil model'
       ELSE
          r_lower = TRIM( r_lower ) // ' pt(0)     = pt_surface'
       ENDIF
    ELSEIF ( ibc_pt_b == 1 )  THEN
       r_lower = TRIM( r_lower ) // ' pt(0)     = pt(1)'
    ELSEIF ( ibc_pt_b == 2 )  THEN
       r_lower = TRIM( r_lower ) // ' pt(0)     = from coupled model'
    ENDIF
    IF ( ibc_pt_t == 0 )  THEN
       r_upper  = TRIM( r_upper  ) // ' pt(nzt+1) = pt_top'
    ELSEIF( ibc_pt_t == 1 )  THEN
       r_upper  = TRIM( r_upper  ) // ' pt(nzt+1) = pt(nzt)'
    ELSEIF( ibc_pt_t == 2 )  THEN
       r_upper  = TRIM( r_upper  ) // ' pt(nzt+1) = pt(nzt) + dpt/dz_ini'

    ENDIF

    WRITE ( io, 300 )  r_lower, r_upper

    IF ( .NOT. constant_diffusion )  THEN
       IF ( ibc_e_b == 1 )  THEN
          r_lower = 'e(0)     = e(1)'
       ELSE
          r_lower = 'e(0)     = e(1) = (u*/0.1)**2'
       ENDIF
       r_upper = 'e(nzt+1) = e(nzt) = e(nzt-1)'

       WRITE ( io, 301 )  'e', r_lower, r_upper       

    ENDIF

    IF ( ocean )  THEN
       IF ( ibc_sa_b == 0 ) THEN
          r_lower = 'sa(0)    = sa_surface'
       ELSE
          r_lower = 'sa(0)    = sa(1)'
       ENDIF
       IF ( ibc_sa_t == 0 )  THEN
          r_upper =  'sa(nzt+1) = sa_top'
       ELSE
          r_upper =  'sa(nzt+1) = sa(nzt)'
       ENDIF
       WRITE ( io, 301 ) 'sa', r_lower, r_upper
    ENDIF

    IF ( humidity )  THEN
       IF ( ibc_q_b == 0 )  THEN
          IF ( land_surface )  THEN
             r_lower = 'q(0)     = from soil model'
          ELSE
             r_lower = 'q(0)     = q_surface'
          ENDIF

       ELSE
          r_lower = 'q(0)      = q(1)'
       ENDIF
       IF ( ibc_q_t == 0 )  THEN
          r_upper =  'q(nzt+1) = q_top'
       ELSE
          r_upper =  'q(nzt+1) = q(nzt) + dq/dz'
       ENDIF
       WRITE ( io, 301 ) 'q', r_lower, r_upper
    ENDIF

    IF ( passive_scalar )  THEN
       IF ( ibc_s_b == 0 )  THEN
          r_lower = 's(0)      = s_surface'
       ELSE
          r_lower = 's(0)      = s(1)'
       ENDIF
       IF ( ibc_s_t == 0 )  THEN
          r_upper =  's(nzt+1) = s_top'
       ELSEIF ( ibc_s_t == 1 )  THEN
          r_upper =  's(nzt+1) = s(nzt)'
       ELSEIF ( ibc_s_t == 2 )  THEN
          r_upper =  's(nzt+1) = s(nzt) + ds/dz'
       ENDIF
       WRITE ( io, 301 ) 's', r_lower, r_upper
    ENDIF

    IF ( use_surface_fluxes )  THEN
       WRITE ( io, 303 )
       IF ( constant_heatflux )  THEN
          IF ( large_scale_forcing .AND. lsf_surf )  THEN
             IF ( surf_def_h(0)%ns >= 1 )  WRITE ( io, 306 )  surf_def_h(0)%shf(1)
          ELSE
             WRITE ( io, 306 )  surface_heatflux
          ENDIF
          IF ( random_heatflux )  WRITE ( io, 307 )
       ENDIF
       IF ( humidity  .AND.  constant_waterflux )  THEN
          IF ( large_scale_forcing .AND. lsf_surf )  THEN
             WRITE ( io, 311 ) surf_def_h(0)%qsws(1)
          ELSE
             WRITE ( io, 311 ) surface_waterflux
          ENDIF
       ENDIF
       IF ( passive_scalar  .AND.  constant_scalarflux )  THEN
          WRITE ( io, 313 ) surface_scalarflux
       ENDIF
    ENDIF

    IF ( use_top_fluxes )  THEN
       WRITE ( io, 304 )
       IF ( coupling_mode == 'uncoupled' )  THEN
          WRITE ( io, 320 )  top_momentumflux_u, top_momentumflux_v
          IF ( constant_top_heatflux )  THEN
             WRITE ( io, 306 )  top_heatflux
          ENDIF
       ELSEIF ( coupling_mode == 'ocean_to_atmosphere' )  THEN
          WRITE ( io, 316 )
       ENDIF
       IF ( ocean  .AND.  constant_top_salinityflux )                          &
          WRITE ( io, 309 )  top_salinityflux
       IF ( humidity       )  WRITE ( io, 315 )
       IF ( passive_scalar .AND.  constant_top_scalarflux )                    &
          WRITE ( io, 302 ) top_scalarflux
    ENDIF

    IF ( .NOT. TRIM(constant_flux_layer) == 'none' )  THEN
       WRITE ( io, 305 )  (zu(1)-zu(0)), roughness_length,                     &
                          z0h_factor*roughness_length, kappa,                  &
                          zeta_min, zeta_max
       IF ( .NOT. constant_heatflux )  WRITE ( io, 308 )
       IF ( humidity  .AND.  .NOT. constant_waterflux )  THEN
          WRITE ( io, 312 )
       ENDIF
       IF ( passive_scalar  .AND.  .NOT. constant_scalarflux )  THEN
          WRITE ( io, 314 )
       ENDIF
    ELSE
       IF ( INDEX(initializing_actions, 'set_1d-model_profiles') /= 0 )  THEN
          WRITE ( io, 310 )  zeta_min, zeta_max
       ENDIF
    ENDIF

    WRITE ( io, 317 )  bc_lr, bc_ns
    IF ( .NOT. bc_lr_cyc  .OR.  .NOT. bc_ns_cyc )  THEN
       WRITE ( io, 318 )  use_cmax, pt_damping_width, pt_damping_factor       
       IF ( turbulent_inflow )  THEN
          IF ( .NOT. recycling_yshift ) THEN
             WRITE ( io, 319 )  recycling_width, recycling_plane, &
                                inflow_damping_height, inflow_damping_width
          ELSE
             WRITE ( io, 322 )  recycling_width, recycling_plane, &
                                inflow_damping_height, inflow_damping_width
          END IF
       ENDIF
       IF ( turbulent_outflow )  THEN
          WRITE ( io, 323 )  outflow_source_plane, INT(outflow_source_plane/dx)
       ENDIF
    ENDIF

!
!-- Initial Profiles
    WRITE ( io, 321 )
!
!-- Initial wind profiles
    IF ( u_profile(1) /= 9999999.9_wp )  WRITE ( io, 427 )

!
!-- Initial temperature profile
!-- Building output strings, starting with surface temperature
    WRITE ( temperatures, '(F6.2)' )  pt_surface
    gradients = '------'
    slices = '     0'
    coordinates = '   0.0'
    i = 1
    DO  WHILE ( pt_vertical_gradient_level_ind(i) /= -9999 )

       WRITE (coor_chr,'(F7.2)')  pt_init(pt_vertical_gradient_level_ind(i))
       temperatures = TRIM( temperatures ) // ' ' // TRIM( coor_chr )

       WRITE (coor_chr,'(F7.2)')  pt_vertical_gradient(i)
       gradients = TRIM( gradients ) // ' ' // TRIM( coor_chr )

       WRITE (coor_chr,'(I7)')  pt_vertical_gradient_level_ind(i)
       slices = TRIM( slices ) // ' ' // TRIM( coor_chr )

       WRITE (coor_chr,'(F7.1)')  pt_vertical_gradient_level(i)
       coordinates = TRIM( coordinates ) // ' '  // TRIM( coor_chr )

       IF ( i == 10 )  THEN
          EXIT
       ELSE
          i = i + 1
       ENDIF

    ENDDO

    IF ( .NOT. nudging )  THEN
       WRITE ( io, 420 )  TRIM( coordinates ), TRIM( temperatures ), &
                          TRIM( gradients ), TRIM( slices )
    ELSE
       WRITE ( io, 428 ) 
    ENDIF

!
!-- Initial humidity profile
!-- Building output strings, starting with surface humidity
    IF ( humidity )  THEN
       WRITE ( temperatures, '(E8.1)' )  q_surface
       gradients = '--------'
       slices = '       0'
       coordinates = '     0.0'
       i = 1
       DO  WHILE ( q_vertical_gradient_level_ind(i) /= -9999 )
          
          WRITE (coor_chr,'(E8.1,4X)')  q_init(q_vertical_gradient_level_ind(i))
          temperatures = TRIM( temperatures ) // '  ' // TRIM( coor_chr )

          WRITE (coor_chr,'(E8.1,4X)')  q_vertical_gradient(i)
          gradients = TRIM( gradients ) // '  ' // TRIM( coor_chr )
          
          WRITE (coor_chr,'(I8,4X)')  q_vertical_gradient_level_ind(i)
          slices = TRIM( slices ) // '  ' // TRIM( coor_chr )
          
          WRITE (coor_chr,'(F8.1,4X)')  q_vertical_gradient_level(i)
          coordinates = TRIM( coordinates ) // '  '  // TRIM( coor_chr )

          IF ( i == 10 )  THEN
             EXIT
          ELSE
             i = i + 1
          ENDIF

       ENDDO

       IF ( .NOT. nudging )  THEN
          WRITE ( io, 421 )  TRIM( coordinates ), TRIM( temperatures ),        &
                             TRIM( gradients ), TRIM( slices )
       ENDIF
    ENDIF
!
!-- Initial scalar profile
!-- Building output strings, starting with surface humidity
    IF ( passive_scalar )  THEN
       WRITE ( temperatures, '(E8.1)' )  s_surface
       gradients = '--------'
       slices = '       0'
       coordinates = '     0.0'
       i = 1
       DO  WHILE ( s_vertical_gradient_level_ind(i) /= -9999 )
          
          WRITE (coor_chr,'(E8.1,4X)')  s_init(s_vertical_gradient_level_ind(i))
          temperatures = TRIM( temperatures ) // '  ' // TRIM( coor_chr )

          WRITE (coor_chr,'(E8.1,4X)')  s_vertical_gradient(i)
          gradients = TRIM( gradients ) // '  ' // TRIM( coor_chr )
          
          WRITE (coor_chr,'(I8,4X)')  s_vertical_gradient_level_ind(i)
          slices = TRIM( slices ) // '  ' // TRIM( coor_chr )
          
          WRITE (coor_chr,'(F8.1,4X)')  s_vertical_gradient_level(i)
          coordinates = TRIM( coordinates ) // '  '  // TRIM( coor_chr )

          IF ( i == 10 )  THEN
             EXIT
          ELSE
             i = i + 1
          ENDIF

       ENDDO

       WRITE ( io, 422 )  TRIM( coordinates ), TRIM( temperatures ),           &
                          TRIM( gradients ), TRIM( slices )
    ENDIF    

!
!-- Initial salinity profile
!-- Building output strings, starting with surface salinity
    IF ( ocean )  THEN
       WRITE ( temperatures, '(F6.2)' )  sa_surface
       gradients = '------'
       slices = '     0'
       coordinates = '   0.0'
       i = 1
       DO  WHILE ( sa_vertical_gradient_level_ind(i) /= -9999 )

          WRITE (coor_chr,'(F7.2)')  sa_init(sa_vertical_gradient_level_ind(i))
          temperatures = TRIM( temperatures ) // ' ' // TRIM( coor_chr )

          WRITE (coor_chr,'(F7.2)')  sa_vertical_gradient(i)
          gradients = TRIM( gradients ) // ' ' // TRIM( coor_chr )

          WRITE (coor_chr,'(I7)')  sa_vertical_gradient_level_ind(i)
          slices = TRIM( slices ) // ' ' // TRIM( coor_chr )

          WRITE (coor_chr,'(F7.1)')  sa_vertical_gradient_level(i)
          coordinates = TRIM( coordinates ) // ' '  // TRIM( coor_chr )

          IF ( i == 10 )  THEN
             EXIT
          ELSE
             i = i + 1
          ENDIF

       ENDDO

       WRITE ( io, 425 )  TRIM( coordinates ), TRIM( temperatures ), &
                          TRIM( gradients ), TRIM( slices )
    ENDIF


!
!-- Listing of 1D-profiles
    WRITE ( io, 325 )  dt_dopr_listing
    IF ( averaging_interval_pr /= 0.0_wp )  THEN
       WRITE ( io, 326 )  averaging_interval_pr, dt_averaging_input_pr
    ENDIF

!
!-- DATA output
    WRITE ( io, 330 )
    IF ( averaging_interval_pr /= 0.0_wp )  THEN
       WRITE ( io, 326 )  averaging_interval_pr, dt_averaging_input_pr
    ENDIF

!
!-- 1D-profiles
    dopr_chr = 'Profile:'
    IF ( dopr_n /= 0 )  THEN
       WRITE ( io, 331 )

       output_format = ''
       output_format = netcdf_data_format_string
       IF ( netcdf_deflate == 0 )  THEN
          WRITE ( io, 344 )  output_format
       ELSE
          WRITE ( io, 354 )  TRIM( output_format ), netcdf_deflate
       ENDIF

       DO  i = 1, dopr_n
          dopr_chr = TRIM( dopr_chr ) // ' ' // TRIM( data_output_pr(i) ) // ','
          IF ( LEN_TRIM( dopr_chr ) >= 60 )  THEN
             WRITE ( io, 332 )  dopr_chr
             dopr_chr = '       :'
          ENDIF
       ENDDO

       IF ( dopr_chr /= '' )  THEN
          WRITE ( io, 332 )  dopr_chr
       ENDIF
       WRITE ( io, 333 )  dt_dopr, averaging_interval_pr, dt_averaging_input_pr
       IF ( skip_time_dopr /= 0.0_wp )  WRITE ( io, 339 )  skip_time_dopr
    ENDIF

!
!-- 2D-arrays
    DO  av = 0, 1

       i = 1
       do2d_xy = ''
       do2d_xz = ''
       do2d_yz = ''
       DO  WHILE ( do2d(av,i) /= ' ' )

          l = MAX( 2, LEN_TRIM( do2d(av,i) ) )
          do2d_mode = do2d(av,i)(l-1:l)

          SELECT CASE ( do2d_mode )
             CASE ( 'xy' )
                ll = LEN_TRIM( do2d_xy )
                do2d_xy = do2d_xy(1:ll) // ' ' // do2d(av,i)(1:l-3) // ','
             CASE ( 'xz' )
                ll = LEN_TRIM( do2d_xz )
                do2d_xz = do2d_xz(1:ll) // ' ' // do2d(av,i)(1:l-3) // ','
             CASE ( 'yz' )
                ll = LEN_TRIM( do2d_yz )
                do2d_yz = do2d_yz(1:ll) // ' ' // do2d(av,i)(1:l-3) // ','
          END SELECT

          i = i + 1

       ENDDO

       IF ( ( ( do2d_xy /= ''  .AND.  section(1,1) /= -9999 )  .OR.    &
              ( do2d_xz /= ''  .AND.  section(1,2) /= -9999 )  .OR.    &
              ( do2d_yz /= ''  .AND.  section(1,3) /= -9999 ) ) )  THEN

          IF (  av == 0 )  THEN
             WRITE ( io, 334 )  ''
          ELSE
             WRITE ( io, 334 )  '(time-averaged)'
          ENDIF

          IF ( do2d_at_begin )  THEN
             begin_chr = 'and at the start'
          ELSE
             begin_chr = ''
          ENDIF

          output_format = ''
          output_format = netcdf_data_format_string
          IF ( netcdf_deflate == 0 )  THEN
             WRITE ( io, 344 )  output_format
          ELSE
             WRITE ( io, 354 )  TRIM( output_format ), netcdf_deflate
          ENDIF

          IF ( do2d_xy /= ''  .AND.  section(1,1) /= -9999 )  THEN
             i = 1
             slices = '/'
             coordinates = '/'
!
!--          Building strings with index and coordinate information of the
!--          slices
             DO  WHILE ( section(i,1) /= -9999 )

                WRITE (section_chr,'(I5)')  section(i,1)
                section_chr = ADJUSTL( section_chr )
                slices = TRIM( slices ) // TRIM( section_chr ) // '/'

                IF ( section(i,1) == -1 )  THEN
                   WRITE (coor_chr,'(F10.1)')  -1.0_wp
                ELSE
                   WRITE (coor_chr,'(F10.1)')  zu(section(i,1))
                ENDIF
                coor_chr = ADJUSTL( coor_chr )
                coordinates = TRIM( coordinates ) // TRIM( coor_chr ) // '/'

                i = i + 1
             ENDDO
             IF ( av == 0 )  THEN
                WRITE ( io, 335 )  'XY', do2d_xy, dt_do2d_xy, &
                                   TRIM( begin_chr ), 'k', TRIM( slices ), &
                                   TRIM( coordinates )
                IF ( skip_time_do2d_xy /= 0.0_wp )  THEN
                   WRITE ( io, 339 )  skip_time_do2d_xy
                ENDIF
             ELSE
                WRITE ( io, 342 )  'XY', do2d_xy, dt_data_output_av, &
                                   TRIM( begin_chr ), averaging_interval, &
                                   dt_averaging_input, 'k', TRIM( slices ), &
                                   TRIM( coordinates )
                IF ( skip_time_data_output_av /= 0.0_wp )  THEN
                   WRITE ( io, 339 )  skip_time_data_output_av
                ENDIF
             ENDIF
             IF ( netcdf_data_format > 4 )  THEN
                WRITE ( io, 352 )  ntdim_2d_xy(av)
             ELSE
                WRITE ( io, 353 )
             ENDIF
          ENDIF

          IF ( do2d_xz /= ''  .AND.  section(1,2) /= -9999 )  THEN
             i = 1
             slices = '/'
             coordinates = '/'
!
!--          Building strings with index and coordinate information of the
!--          slices
             DO  WHILE ( section(i,2) /= -9999 )

                WRITE (section_chr,'(I5)')  section(i,2)
                section_chr = ADJUSTL( section_chr )
                slices = TRIM( slices ) // TRIM( section_chr ) // '/'

                WRITE (coor_chr,'(F10.1)')  section(i,2) * dy
                coor_chr = ADJUSTL( coor_chr )
                coordinates = TRIM( coordinates ) // TRIM( coor_chr ) // '/'

                i = i + 1
             ENDDO
             IF ( av == 0 )  THEN
                WRITE ( io, 335 )  'XZ', do2d_xz, dt_do2d_xz, &
                                   TRIM( begin_chr ), 'j', TRIM( slices ), &
                                   TRIM( coordinates )
                IF ( skip_time_do2d_xz /= 0.0_wp )  THEN
                   WRITE ( io, 339 )  skip_time_do2d_xz
                ENDIF
             ELSE
                WRITE ( io, 342 )  'XZ', do2d_xz, dt_data_output_av, &
                                   TRIM( begin_chr ), averaging_interval, &
                                   dt_averaging_input, 'j', TRIM( slices ), &
                                   TRIM( coordinates )
                IF ( skip_time_data_output_av /= 0.0_wp )  THEN
                   WRITE ( io, 339 )  skip_time_data_output_av
                ENDIF
             ENDIF
             IF ( netcdf_data_format > 4 )  THEN
                WRITE ( io, 352 )  ntdim_2d_xz(av)
             ELSE
                WRITE ( io, 353 )
             ENDIF
          ENDIF

          IF ( do2d_yz /= ''  .AND.  section(1,3) /= -9999 )  THEN
             i = 1
             slices = '/'
             coordinates = '/'
!
!--          Building strings with index and coordinate information of the
!--          slices
             DO  WHILE ( section(i,3) /= -9999 )

                WRITE (section_chr,'(I5)')  section(i,3)
                section_chr = ADJUSTL( section_chr )
                slices = TRIM( slices ) // TRIM( section_chr ) // '/'

                WRITE (coor_chr,'(F10.1)')  section(i,3) * dx
                coor_chr = ADJUSTL( coor_chr )
                coordinates = TRIM( coordinates ) // TRIM( coor_chr ) // '/'

                i = i + 1
             ENDDO
             IF ( av == 0 )  THEN
                WRITE ( io, 335 )  'YZ', do2d_yz, dt_do2d_yz, &
                                   TRIM( begin_chr ), 'i', TRIM( slices ), &
                                   TRIM( coordinates )
                IF ( skip_time_do2d_yz /= 0.0_wp )  THEN
                   WRITE ( io, 339 )  skip_time_do2d_yz
                ENDIF
             ELSE
                WRITE ( io, 342 )  'YZ', do2d_yz, dt_data_output_av, &
                                   TRIM( begin_chr ), averaging_interval, &
                                   dt_averaging_input, 'i', TRIM( slices ), &
                                   TRIM( coordinates )
                IF ( skip_time_data_output_av /= 0.0_wp )  THEN
                   WRITE ( io, 339 )  skip_time_data_output_av
                ENDIF
             ENDIF
             IF ( netcdf_data_format > 4 )  THEN
                WRITE ( io, 352 )  ntdim_2d_yz(av)
             ELSE
                WRITE ( io, 353 )
             ENDIF
          ENDIF

       ENDIF

    ENDDO

!
!-- 3d-arrays
    DO  av = 0, 1

       i = 1
       do3d_chr = ''
       DO  WHILE ( do3d(av,i) /= ' ' )

          do3d_chr = TRIM( do3d_chr ) // ' ' // TRIM( do3d(av,i) ) // ','
          i = i + 1

       ENDDO

       IF ( do3d_chr /= '' )  THEN
          IF ( av == 0 )  THEN
             WRITE ( io, 336 )  ''
          ELSE
             WRITE ( io, 336 )  '(time-averaged)'
          ENDIF

          output_format = netcdf_data_format_string
          IF ( netcdf_deflate == 0 )  THEN
             WRITE ( io, 344 )  output_format
          ELSE
             WRITE ( io, 354 )  TRIM( output_format ), netcdf_deflate
          ENDIF

          IF ( do3d_at_begin )  THEN
             begin_chr = 'and at the start'
          ELSE
             begin_chr = ''
          ENDIF
          IF ( av == 0 )  THEN
             WRITE ( io, 337 )  do3d_chr, dt_do3d, TRIM( begin_chr ), &
                                zu(nz_do3d), nz_do3d
          ELSE
             WRITE ( io, 343 )  do3d_chr, dt_data_output_av,           &
                                TRIM( begin_chr ), averaging_interval, &
                                dt_averaging_input, zu(nz_do3d), nz_do3d
          ENDIF

          IF ( netcdf_data_format > 4 )  THEN
             WRITE ( io, 352 )  ntdim_3d(av)
          ELSE
             WRITE ( io, 353 )
          ENDIF

          IF ( av == 0 )  THEN
             IF ( skip_time_do3d /= 0.0_wp )  THEN
                WRITE ( io, 339 )  skip_time_do3d
             ENDIF
          ELSE
             IF ( skip_time_data_output_av /= 0.0_wp )  THEN
                WRITE ( io, 339 )  skip_time_data_output_av
             ENDIF
          ENDIF

       ENDIF

    ENDDO

!
!-- masked arrays
    IF ( masks > 0 )  WRITE ( io, 345 )  &
         mask_scale_x, mask_scale_y, mask_scale_z
    DO  mid = 1, masks
       DO  av = 0, 1

          i = 1
          domask_chr = ''
          DO  WHILE ( domask(mid,av,i) /= ' ' )
             domask_chr = TRIM( domask_chr ) // ' ' //  &
                          TRIM( domask(mid,av,i) ) // ','
             i = i + 1
          ENDDO

          IF ( domask_chr /= '' )  THEN
             IF ( av == 0 )  THEN
                WRITE ( io, 346 )  '', mid
             ELSE
                WRITE ( io, 346 )  ' (time-averaged)', mid
             ENDIF

             output_format = netcdf_data_format_string
!--          Parallel output not implemented for mask data, hence
!--          output_format must be adjusted.
             IF ( netcdf_data_format == 5 ) output_format = 'netCDF4/HDF5'
             IF ( netcdf_data_format == 6 ) output_format = 'netCDF4/HDF5 classic'
             IF ( netcdf_deflate == 0 )  THEN
                WRITE ( io, 344 )  output_format
             ELSE
                WRITE ( io, 354 )  TRIM( output_format ), netcdf_deflate
             ENDIF

             IF ( av == 0 )  THEN
                WRITE ( io, 347 )  domask_chr, dt_domask(mid)
             ELSE
                WRITE ( io, 348 )  domask_chr, dt_data_output_av, &
                                   averaging_interval, dt_averaging_input
             ENDIF

             IF ( av == 0 )  THEN
                IF ( skip_time_domask(mid) /= 0.0_wp )  THEN
                   WRITE ( io, 339 )  skip_time_domask(mid)
                ENDIF
             ELSE
                IF ( skip_time_data_output_av /= 0.0_wp )  THEN
                   WRITE ( io, 339 )  skip_time_data_output_av
                ENDIF
             ENDIF
!
!--          output locations
             DO  dim = 1, 3
                IF ( mask(mid,dim,1) >= 0.0_wp )  THEN
                   count = 0
                   DO  WHILE ( mask(mid,dim,count+1) >= 0.0_wp )
                      count = count + 1
                   ENDDO
                   WRITE ( io, 349 )  dir(dim), dir(dim), mid, dir(dim), &
                                      mask(mid,dim,:count)
                ELSEIF ( mask_loop(mid,dim,1) < 0.0_wp .AND.  &
                         mask_loop(mid,dim,2) < 0.0_wp .AND.  &
                         mask_loop(mid,dim,3) == 0.0_wp )  THEN
                   WRITE ( io, 350 )  dir(dim), dir(dim)
                ELSEIF ( mask_loop(mid,dim,3) == 0.0_wp )  THEN
                   WRITE ( io, 351 )  dir(dim), dir(dim), mid, dir(dim), &
                                      mask_loop(mid,dim,1:2)
                ELSE
                   WRITE ( io, 351 )  dir(dim), dir(dim), mid, dir(dim), &
                                      mask_loop(mid,dim,1:3)
                ENDIF
             ENDDO
          ENDIF

       ENDDO
    ENDDO

!
!-- Timeseries
    IF ( dt_dots /= 9999999.9_wp )  THEN
       WRITE ( io, 340 )

       output_format = netcdf_data_format_string
       IF ( netcdf_deflate == 0 )  THEN
          WRITE ( io, 344 )  output_format
       ELSE
          WRITE ( io, 354 )  TRIM( output_format ), netcdf_deflate
       ENDIF
       WRITE ( io, 341 )  dt_dots
    ENDIF

#if defined( __dvrp_graphics )
!
!-- Dvrp-output
    IF ( dt_dvrp /= 9999999.9_wp )  THEN
       WRITE ( io, 360 )  dt_dvrp, TRIM( dvrp_output ), TRIM( dvrp_host ), &
                          TRIM( dvrp_username ), TRIM( dvrp_directory )
       i = 1
       l = 0
       m = 0
       DO WHILE ( mode_dvrp(i) /= ' ' )
          IF ( mode_dvrp(i)(1:10) == 'isosurface' )  THEN
             READ ( mode_dvrp(i), '(10X,I2)' )  j
             l = l + 1
             IF ( do3d(0,j) /= ' ' )  THEN
                WRITE ( io, 361 )  TRIM( do3d(0,j) ), threshold(l), &
                                   isosurface_color(:,l)
             ENDIF
          ELSEIF ( mode_dvrp(i)(1:6) == 'slicer' )  THEN
             READ ( mode_dvrp(i), '(6X,I2)' )  j
             m = m + 1
             IF ( do2d(0,j) /= ' ' )  THEN
                WRITE ( io, 362 )  TRIM( do2d(0,j) ), &
                                   slicer_range_limits_dvrp(:,m)
             ENDIF
          ENDIF
          i = i + 1
       ENDDO

       WRITE ( io, 365 )  groundplate_color, superelevation_x, &
                          superelevation_y, superelevation, clip_dvrp_l, &
                          clip_dvrp_r, clip_dvrp_s, clip_dvrp_n

       IF ( TRIM( topography ) /= 'flat' )  THEN
          WRITE ( io, 366 )  topography_color
          IF ( cluster_size > 1 )  THEN
             WRITE ( io, 367 )  cluster_size
          ENDIF
       ENDIF

    ENDIF
#endif
!
!-- Output of virtual flight information
    IF ( virtual_flight )  CALL flight_header( io )

!
!-- Output of spectra related quantities
    IF ( calculate_spectra )  CALL spectra_header( io )

    WRITE ( io, 99 )

!
!-- Physical quantities
    WRITE ( io, 400 )

!
!-- Geostrophic parameters
    WRITE ( io, 410 )  latitude, longitude, omega, f, fs

 !
!-- Geostrophic parameters
    WRITE ( io, 456 )  day_of_year_init, time_utc_init
    
!
!-- Other quantities
    WRITE ( io, 411 )  g

    WRITE ( io, 412 )  TRIM( reference_state )
    IF ( use_single_reference_value )  THEN
       IF ( ocean )  THEN
          WRITE ( io, 413 )  prho_reference
       ELSE
          WRITE ( io, 414 )  pt_reference
       ENDIF
    ENDIF

!
!-- Cloud physics parameters
    IF ( cloud_physics )  THEN
       WRITE ( io, 415 )
       WRITE ( io, 416 ) surface_pressure, r_d, rho_surface, cp, l_v
       IF ( microphysics_seifert )  THEN
          WRITE ( io, 510 ) 1.0E-6_wp * nc_const
          WRITE ( io, 511 ) c_sedimentation
       ENDIF
    ENDIF

!
!-- Cloud physcis parameters / quantities / numerical methods
    WRITE ( io, 430 )
    IF ( humidity .AND. .NOT. cloud_physics .AND. .NOT. cloud_droplets)  THEN
       WRITE ( io, 431 )
    ELSEIF ( humidity  .AND.  cloud_physics )  THEN
       WRITE ( io, 432 )
       IF ( cloud_top_radiation )  WRITE ( io, 132 )
       IF ( microphysics_kessler )  THEN
          WRITE ( io, 133 )
       ELSEIF ( microphysics_seifert )  THEN
          IF ( cloud_water_sedimentation )  WRITE ( io, 506 )
          WRITE ( io, 505 )
          IF ( collision_turbulence )  WRITE ( io, 507 )
          IF ( ventilation_effect )  WRITE ( io, 508 )
          IF ( limiter_sedimentation )  WRITE ( io, 509 )
       ENDIF
    ELSEIF ( humidity  .AND.  cloud_droplets )  THEN
       WRITE ( io, 433 )
       IF ( curvature_solution_effects )  WRITE ( io, 434 )
       IF ( collision_kernel /= 'none' )  THEN
          WRITE ( io, 435 )  TRIM( collision_kernel )
          IF ( collision_kernel(6:9) == 'fast' )  THEN
             WRITE ( io, 436 )  radius_classes, dissipation_classes
          ENDIF
       ELSE
          WRITE ( io, 437 )
       ENDIF
    ENDIF

!
!-- LES / turbulence parameters
    WRITE ( io, 450 )

!--
! ... LES-constants used must still be added here
!--
    IF ( constant_diffusion )  THEN
       WRITE ( io, 451 )  km_constant, km_constant/prandtl_number, &
                          prandtl_number
    ENDIF
    IF ( .NOT. constant_diffusion)  THEN
       IF ( e_init > 0.0_wp )  WRITE ( io, 455 )  e_init
       IF ( e_min > 0.0_wp )  WRITE ( io, 454 )  e_min
       IF ( wall_adjustment )  WRITE ( io, 453 )  wall_adjustment_factor
    ENDIF
    IF ( rans_mode )  THEN
       WRITE ( io, 457 )  rans_const_c, rans_const_sigma
    ENDIF
!
!-- Special actions during the run
    WRITE ( io, 470 )
    IF ( create_disturbances )  THEN
       WRITE ( io, 471 )  dt_disturb, disturbance_amplitude,                   &
                          zu(disturbance_level_ind_b), disturbance_level_ind_b,&
                          zu(disturbance_level_ind_t), disturbance_level_ind_t
       IF ( .NOT. bc_lr_cyc  .OR.  .NOT. bc_ns_cyc )  THEN
          WRITE ( io, 472 )  inflow_disturbance_begin, inflow_disturbance_end
       ELSE
          WRITE ( io, 473 )  disturbance_energy_limit
       ENDIF
       WRITE ( io, 474 )  TRIM( random_generator )
    ENDIF
    IF ( pt_surface_initial_change /= 0.0_wp )  THEN
       WRITE ( io, 475 )  pt_surface_initial_change
    ENDIF
    IF ( humidity  .AND.  q_surface_initial_change /= 0.0_wp )  THEN
       WRITE ( io, 476 )  q_surface_initial_change       
    ENDIF
    IF ( passive_scalar  .AND.  q_surface_initial_change /= 0.0_wp )  THEN
       WRITE ( io, 477 )  q_surface_initial_change       
    ENDIF

    IF ( particle_advection )  THEN
!
!--    Particle attributes
       WRITE ( io, 480 )  particle_advection_start, dt_prel, bc_par_lr, &
                          bc_par_ns, bc_par_b, bc_par_t, particle_maximum_age, &
                          end_time_prel
       IF ( use_sgs_for_particles )  WRITE ( io, 488 )  dt_min_part
       IF ( random_start_position )  WRITE ( io, 481 )
       IF ( seed_follows_topography )  WRITE ( io, 496 )
       IF ( particles_per_point > 1 )  WRITE ( io, 489 )  particles_per_point
       WRITE ( io, 495 )  total_number_of_particles
       IF ( dt_write_particle_data /= 9999999.9_wp )  THEN
          WRITE ( io, 485 )  dt_write_particle_data
          IF ( netcdf_data_format > 1 )  THEN
             output_format = 'netcdf (64 bit offset) and binary'
          ELSE
             output_format = 'netcdf and binary'
          ENDIF
          IF ( netcdf_deflate == 0 )  THEN
             WRITE ( io, 344 )  output_format
          ELSE
             WRITE ( io, 354 )  TRIM( output_format ), netcdf_deflate
          ENDIF
       ENDIF
       IF ( dt_dopts /= 9999999.9_wp )  WRITE ( io, 494 )  dt_dopts
       IF ( write_particle_statistics )  WRITE ( io, 486 )

       WRITE ( io, 487 )  number_of_particle_groups

       DO  i = 1, number_of_particle_groups
          IF ( i == 1  .AND.  density_ratio(i) == 9999999.9_wp )  THEN
             WRITE ( io, 490 )  i, 0.0_wp
             WRITE ( io, 492 )
          ELSE
             WRITE ( io, 490 )  i, radius(i)
             IF ( density_ratio(i) /= 0.0_wp )  THEN
                WRITE ( io, 491 )  density_ratio(i)
             ELSE
                WRITE ( io, 492 )
             ENDIF
          ENDIF
          WRITE ( io, 493 )  psl(i), psr(i), pss(i), psn(i), psb(i), pst(i), &
                             pdx(i), pdy(i), pdz(i)
          IF ( .NOT. vertical_particle_advection(i) )  WRITE ( io, 482 )
       ENDDO

    ENDIF


!
!-- Parameters of 1D-model
    IF ( INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN
       WRITE ( io, 500 )  end_time_1d, dt_run_control_1d, dt_pr_1d, &
                          mixing_length_1d, dissipation_1d
       IF ( damp_level_ind_1d /= nzt+1 )  THEN
          WRITE ( io, 502 )  zu(damp_level_ind_1d), damp_level_ind_1d
       ENDIF
    ENDIF

!
!-- User-defined information
    CALL user_header( io )

    WRITE ( io, 99 )

!
!-- Write buffer contents to disc immediately
    FLUSH( io )

!
!-- Here the FORMATs start

 99 FORMAT (1X,78('-'))
100 FORMAT (/1X,'******************************',4X,44('-')/        &
            1X,'* ',A,' *',4X,A/                               &
            1X,'******************************',4X,44('-'))
101 FORMAT (35X,'coupled run: ',A/ &
            35X,42('-'))
102 FORMAT (/' Date:                 ',A8,4X,'Run:       ',A20/      &
            ' Time:                 ',A8,4X,'Run-No.:   ',I2.2/     &
            ' Run on host:        ',A10)
#if defined( __parallel )
103 FORMAT (' Number of PEs:',10X,I6,4X,'Processor grid (x,y): (',I4,',',I4, &
              ')',1X,A)
104 FORMAT (' Number of PEs:',10X,I6,4X,'Tasks:',I4,'   threads per task:',I4/ &
              35X,'Processor grid (x,y): (',I4,',',I4,')',1X,A)
105 FORMAT (35X,'One additional PE is used to handle'/37X,'the dvrp output!')
107 FORMAT (35X,'A 1d-decomposition along ',A,' is used')
108 FORMAT (35X,'Max. # of parallel I/O streams is ',I5)
109 FORMAT (35X,'Precursor run for coupled atmos-ocean run'/ &
            35X,42('-'))
114 FORMAT (35X,'Coupled atmosphere-ocean run following'/ &
            35X,'independent precursor runs'/             &
            35X,42('-'))
#endif
110 FORMAT (/' Numerical Schemes:'/ &
             ' -----------------'/)
124 FORMAT (' --> Use the ',A,' turbulence closure (',A,' mode).')
121 FORMAT (' --> Use the ',A,' approximation for the model equations.')
111 FORMAT (' --> Solve perturbation pressure via FFT using ',A,' routines')
112 FORMAT (' --> Solve perturbation pressure via SOR-Red/Black-Schema'/ &
            '     Iterations (initial/other): ',I3,'/',I3,'  omega =',F6.3)
113 FORMAT (' --> Momentum advection via Piascek-Williams-Scheme (Form C3)', &
                  ' or Upstream')
115 FORMAT ('     FFT and transpositions are overlapping')
116 FORMAT (' --> Scalar advection via Piascek-Williams-Scheme (Form C3)', &
                  ' or Upstream')
118 FORMAT (' --> Scalar advection via Bott-Chlond-Scheme')
119 FORMAT (' --> Galilei-Transform applied to horizontal advection:'/ &
            '     translation velocity = ',A/ &
            '     distance advected ',A,':  ',F8.3,' km(x)  ',F8.3,' km(y)')
122 FORMAT (' --> Time differencing scheme: ',A)
123 FORMAT (' --> Rayleigh-Damping active, starts ',A,' z = ',F8.2,' m'/ &
            '     maximum damping coefficient:',F6.3, ' 1/s')
129 FORMAT (' --> Additional prognostic equation for the specific humidity')
130 FORMAT (' --> Additional prognostic equation for the total water content')
131 FORMAT (' --> No pt-equation solved. Neutral stratification with pt = ', &
                  F6.2, ' K assumed')
132 FORMAT ('     Parameterization of long-wave radiation processes via'/ &
            '     effective emissivity scheme')
133 FORMAT ('     Precipitation parameterization via Kessler-Scheme')
134 FORMAT (' --> Additional prognostic equation for a passive scalar')
135 FORMAT (' --> Solve perturbation pressure via ',A,' method (', &
                  A,'-cycle)'/ &
            '     number of grid levels:                   ',I2/ &
            '     Gauss-Seidel red/black iterations:       ',I2)
136 FORMAT ('     gridpoints of coarsest subdomain (x,y,z): (',I3,',',I3,',', &
                  I3,')')
137 FORMAT ('     level data gathered on PE0 at level:     ',I2/ &
            '     gridpoints of coarsest subdomain (x,y,z): (',I3,',',I3,',', &
                  I3,')'/ &
            '     gridpoints of coarsest domain (x,y,z):    (',I3,',',I3,',', &
                  I3,')')
139 FORMAT (' --> Loop optimization method: ',A)
140 FORMAT ('     maximum residual allowed:                ',E10.3)
141 FORMAT ('     fixed number of multigrid cycles:        ',I4)
142 FORMAT ('     perturbation pressure is calculated at every Runge-Kutta ', &
                  'step')
143 FORMAT ('     Euler/upstream scheme is used for the SGS turbulent ', &
                  'kinetic energy')
144 FORMAT ('     masking method is used')
150 FORMAT (' --> Volume flow at the right and north boundary will be ', &
                  'conserved'/ &
            '     using the ',A,' mode')
151 FORMAT ('     with u_bulk = ',F7.3,' m/s and v_bulk = ',F7.3,' m/s')
152 FORMAT (' --> External pressure gradient directly prescribed by the user:',&
           /'     ',2(1X,E12.5),'Pa/m in x/y direction', &
           /'     starting from dp_level_b =', F8.3, 'm', A /)
200 FORMAT (//' Run time and time step information:'/ &
             ' ----------------------------------'/)
201 FORMAT ( ' Timestep:             variable     maximum value: ',F6.3,' s', &
             '    CFL-factor:',F5.2)
202 FORMAT ( ' Timestep:          dt = ',F6.3,' s'/)
203 FORMAT ( ' Start time:          ',F9.3,' s'/ &
             ' End time:            ',F9.3,' s')
204 FORMAT ( A,F9.3,' s')
205 FORMAT ( A,F9.3,' s',5X,'restart every',17X,F9.3,' s')
206 FORMAT (/' Time reached:        ',F9.3,' s'/ &
             ' CPU-time used:       ',F9.3,' s     per timestep:               ', &
               '  ',F9.3,' s'/                                                    &
             '                                      per second of simulated tim', &
               'e: ',F9.3,' s')
207 FORMAT ( ' Coupling start time: ',F9.3,' s')
250 FORMAT (//' Computational grid and domain size:'/ &
              ' ----------------------------------'// &
              ' Grid length:      dx =    ',F8.3,' m    dy =    ',F8.3, ' m')
251 FORMAT (  /' Domain size:       x = ',F10.3,' m     y = ',F10.3, &
              ' m  z(u) = ',F10.3,' m'/)
253 FORMAT ( '                dz(',I1,') =    ', F8.3, ' m')
254 FORMAT (//' Number of gridpoints (x,y,z):  (0:',I4,', 0:',I4,', 0:',I4,')'/ &
            ' Subdomain size (x,y,z):        (  ',I4,',   ',I4,',   ',I4,')'/)
260 FORMAT (/' The model has a slope in x-direction. Inclination angle: ',F6.2,&
             ' degrees')
270 FORMAT (//' Topography information:'/ &
              ' ----------------------'// &
              1X,'Topography: ',A)
271 FORMAT (  ' Building size (x/y/z) in m: ',F5.1,' / ',F5.1,' / ',F5.1/ &
              ' Horizontal index bounds (l/r/s/n): ',I4,' / ',I4,' / ',I4, &
                ' / ',I4)
272 FORMAT (  ' Single quasi-2D street canyon of infinite length in ',A, &
              ' direction' / &
              ' Canyon height: ', F6.2, 'm, ch = ', I4, '.'      / &
              ' Canyon position (',A,'-walls): cxl = ', I4,', cxr = ', I4, '.')
273 FORMAT (  ' Tunnel of infinite length in ',A, &
              ' direction' / &
              ' Tunnel height: ', F6.2, / &
              ' Tunnel-wall depth: ', F6.2      / &
              ' Tunnel width: ', F6.2 )
274 FORMAT (  ' Tunnel in ', A, ' direction.' / &
              ' Tunnel height: ', F6.2, / &    
              ' Tunnel-wall depth: ', F6.2      / &
              ' Tunnel width: ', F6.2, / &
              ' Tunnel length: ', F6.2 )
278 FORMAT (' Topography grid definition convention:'/ &
            ' cell edge (staggered grid points'/  &
            ' (u in x-direction, v in y-direction))' /)
279 FORMAT (' Topography grid definition convention:'/ &
            ' cell center (scalar grid points)' /)
280 FORMAT (' Complex terrain simulation is activated.')
281 FORMAT ('    --> Mean inflow profiles are adjusted.' / &
            '    --> Elevation of inflow boundary: ', F7.1, ' m' )
282 FORMAT ('    --> Initial data from 3D-precursor run is shifted' / &
            '        vertically depending on local surface height.')
300 FORMAT (//' Boundary conditions:'/ &
             ' -------------------'// &
             '                     p                    uv             ', &
             '                     pt'// &
             ' B. bound.: ',A/ &
             ' T. bound.: ',A)
301 FORMAT (/'                     ',A// &
             ' B. bound.: ',A/ &
             ' T. bound.: ',A)
303 FORMAT (/' Bottom surface fluxes are used in diffusion terms at k=1')
304 FORMAT (/' Top surface fluxes are used in diffusion terms at k=nzt')
305 FORMAT (//'    Constant flux layer between bottom surface and first ',     &
              'computational u,v-level:'// &
             '       z_mo = ',F6.2,' m   z0 =',F7.4,' m   z0h =',F8.5,&
             ' m   kappa =',F5.2/ &
             '       Rif value range:   ',F8.2,' <= rif <=',F6.2)
306 FORMAT ('       Predefined constant heatflux:   ',F9.6,' K m/s')
307 FORMAT ('       Heatflux has a random normal distribution')
308 FORMAT ('       Predefined surface temperature')
309 FORMAT ('       Predefined constant salinityflux:   ',F9.6,' kg/m3 psu m/s')
310 FORMAT (//'    1D-Model:'// &
             '       Rif value range:   ',F6.2,' <= rif <=',F6.2)
311 FORMAT ('       Predefined constant humidity flux: ',E10.3,' kg/kg m/s')
312 FORMAT ('       Predefined surface humidity')
313 FORMAT ('       Predefined constant scalar flux: ',E10.3,' kg/(m**2 s)')
314 FORMAT ('       Predefined scalar value at the surface')
302 FORMAT ('       Predefined constant scalarflux:   ',F9.6,' kg/(m**2 s)')
315 FORMAT ('       Humidity flux at top surface is 0.0')
316 FORMAT ('       Sensible heatflux and momentum flux from coupled ', &
                    'atmosphere model')
317 FORMAT (//' Lateral boundaries:'/ &
            '       left/right:  ',A/    &
            '       north/south: ',A)
318 FORMAT (/'       use_cmax: ',L1 / &
            '       pt damping layer width = ',F8.2,' m, pt ', &
                    'damping factor =',F7.4)
319 FORMAT ('       turbulence recycling at inflow switched on'/ &
            '       width of recycling domain: ',F7.1,' m   grid index: ',I4/ &
            '       inflow damping height: ',F6.1,' m   width: ',F6.1,' m')
320 FORMAT ('       Predefined constant momentumflux:  u: ',F9.6,' m**2/s**2'/ &
            '                                          v: ',F9.6,' m**2/s**2')
321 FORMAT (//' Initial profiles:'/ &
              ' ----------------')
322 FORMAT ('       turbulence recycling at inflow switched on'/ &
            '       y shift of the recycled inflow turbulence switched on'/ &
            '       width of recycling domain: ',F7.1,' m   grid index: ',I4/ &
            '       inflow damping height: ',F6.1,' m   width: ',F6.1,' m'/)
323 FORMAT ('       turbulent outflow conditon switched on'/ &
            '       position of outflow source plane: ',F7.1,' m   ', &
                    'grid index: ', I4)
325 FORMAT (//' List output:'/ &
             ' -----------'//  &
            '    1D-Profiles:'/    &
            '       Output every             ',F10.2,' s')
326 FORMAT ('       Time averaged over       ',F8.2,' s'/ &
            '       Averaging input every    ',F8.2,' s')
330 FORMAT (//' Data output:'/ &
             ' -----------'/)
331 FORMAT (/'    1D-Profiles:')
332 FORMAT (/'       ',A)
333 FORMAT ('       Output every             ',F8.2,' s',/ &
            '       Time averaged over       ',F8.2,' s'/ &
            '       Averaging input every    ',F8.2,' s')
334 FORMAT (/'    2D-Arrays',A,':')
335 FORMAT (/'       ',A2,'-cross-section  Arrays: ',A/ &
            '       Output every             ',F8.2,' s  ',A/ &
            '       Cross sections at ',A1,' = ',A/ &
            '       scalar-coordinates:   ',A,' m'/)
336 FORMAT (/'    3D-Arrays',A,':')
337 FORMAT (/'       Arrays: ',A/ &
            '       Output every             ',F8.2,' s  ',A/ &
            '       Upper output limit at    ',F8.2,' m  (GP ',I4,')'/)
339 FORMAT ('       No output during initial ',F8.2,' s')
340 FORMAT (/'    Time series:')
341 FORMAT ('       Output every             ',F8.2,' s'/)
342 FORMAT (/'       ',A2,'-cross-section  Arrays: ',A/ &
            '       Output every             ',F8.2,' s  ',A/ &
            '       Time averaged over       ',F8.2,' s'/ &
            '       Averaging input every    ',F8.2,' s'/ &
            '       Cross sections at ',A1,' = ',A/ &
            '       scalar-coordinates:   ',A,' m'/)
343 FORMAT (/'       Arrays: ',A/ &
            '       Output every             ',F8.2,' s  ',A/ &
            '       Time averaged over       ',F8.2,' s'/ &
            '       Averaging input every    ',F8.2,' s'/ &
            '       Upper output limit at    ',F8.2,' m  (GP ',I4,')'/)
344 FORMAT ('       Output format: ',A/)
345 FORMAT (/'    Scaling lengths for output locations of all subsequent mask IDs:',/ &
            '       mask_scale_x (in x-direction): ',F9.3, ' m',/ &
            '       mask_scale_y (in y-direction): ',F9.3, ' m',/ &
            '       mask_scale_z (in z-direction): ',F9.3, ' m' )
346 FORMAT (/'    Masked data output',A,' for mask ID ',I2, ':')
347 FORMAT ('       Variables: ',A/ &
            '       Output every             ',F8.2,' s')
348 FORMAT ('       Variables: ',A/ &
            '       Output every             ',F8.2,' s'/ &
            '       Time averaged over       ',F8.2,' s'/ &
            '       Averaging input every    ',F8.2,' s')
349 FORMAT (/'       Output locations in ',A,'-direction in multiples of ', &
            'mask_scale_',A,' predefined by array mask_',I2.2,'_',A,':'/ &
            13('       ',8(F8.2,',')/) )
350 FORMAT (/'       Output locations in ',A,'-direction: ', &
            'all gridpoints along ',A,'-direction (default).' )
351 FORMAT (/'       Output locations in ',A,'-direction in multiples of ', &
            'mask_scale_',A,' constructed from array mask_',I2.2,'_',A,'_loop:'/ &
            '          loop begin:',F8.2,', end:',F8.2,', stride:',F8.2 )
352 FORMAT  (/'       Number of output time levels allowed: ',I3 /)
353 FORMAT  (/'       Number of output time levels allowed: unlimited' /)
354 FORMAT ('       Output format: ',A, '   compressed with level: ',I1/)
#if defined( __dvrp_graphics )
360 FORMAT ('    Plot-Sequence with dvrp-software:'/ &
            '       Output every      ',F7.1,' s'/ &
            '       Output mode:      ',A/ &
            '       Host / User:      ',A,' / ',A/ &
            '       Directory:        ',A// &
            '       The sequence contains:')
361 FORMAT (/'       Isosurface of "',A,'"    Threshold value: ', E12.3/ &
            '          Isosurface color: (',F4.2,',',F4.2,',',F4.2,') (R,G,B)')
362 FORMAT (/'       Slicer plane ',A/ &
            '       Slicer limits: [',F6.2,',',F6.2,']')
365 FORMAT (/'       Groundplate color: (',F4.2,',',F4.2,',',F4.2,') (R,G,B)'/ &
            '       Superelevation along (x,y,z): (',F4.1,',',F4.1,',',F4.1, &
                     ')'/ &
            '       Clipping limits: from x = ',F9.1,' m to x = ',F9.1,' m'/ &
            '                        from y = ',F9.1,' m to y = ',F9.1,' m')
366 FORMAT (/'       Topography color: (',F4.2,',',F4.2,',',F4.2,') (R,G,B)')
367 FORMAT ('       Polygon reduction for topography: cluster_size = ', I1)
#endif
400 FORMAT (//' Physical quantities:'/ &
              ' -------------------'/)
410 FORMAT ('    Geograph. latitude  :   latitude  = ',F4.1,' degr'/   &
            '    Geograph. longitude :   longitude = ',F4.1,' degr'/   &
            '    Angular velocity    :   omega  =',E10.3,' rad/s'/  &
            '    Coriolis parameter  :   f      = ',F9.6,' 1/s'/    &
            '                            f*     = ',F9.6,' 1/s')
411 FORMAT (/'    Gravity             :   g      = ',F4.1,' m/s**2')
412 FORMAT (/'    Reference state used in buoyancy terms: ',A)
413 FORMAT ('       Reference density in buoyancy terms: ',F8.3,' kg/m**3')
414 FORMAT ('       Reference temperature in buoyancy terms: ',F8.4,' K')
415 FORMAT (/' Cloud physics parameters:'/ &
             ' ------------------------'/)
416 FORMAT ('    Surface pressure   :   p_0   = ',F7.2,' hPa'/      &
            '    Gas constant       :   R     = ',F5.1,' J/(kg K)'/ &
            '    Density of air     :   rho_0 =',F6.3,' kg/m**3'/  &
            '    Specific heat cap. :   c_p   = ',F6.1,' J/(kg K)'/ &
            '    Vapourization heat :   L_v   =',E9.2,' J/kg')
418 FORMAT (/'    Day of the year at model start :   day_init      =     ',I3 &
            /'    UTC time at model start        :   time_utc_init = ',F7.1' s')
420 FORMAT (/'    Characteristic levels of the initial temperature profile:'// &
            '       Height:        ',A,'  m'/ &
            '       Temperature:   ',A,'  K'/ &
            '       Gradient:      ',A,'  K/100m'/ &
            '       Gridpoint:     ',A)
421 FORMAT (/'    Characteristic levels of the initial humidity profile:'// &
            '       Height:      ',A,'  m'/ &
            '       Humidity:    ',A,'  kg/kg'/ &
            '       Gradient:    ',A,'  (kg/kg)/100m'/ &
            '       Gridpoint:   ',A)
422 FORMAT (/'    Characteristic levels of the initial scalar profile:'// &
            '       Height:                  ',A,'  m'/ &
            '       Scalar concentration:    ',A,'  kg/m**3'/ &
            '       Gradient:                ',A,'  (kg/m**3)/100m'/ &
            '       Gridpoint:               ',A)
423 FORMAT (/'    Characteristic levels of the geo. wind component ug:'// &
            '       Height:      ',A,'  m'/ &
            '       ug:          ',A,'  m/s'/ &
            '       Gradient:    ',A,'  1/100s'/ &
            '       Gridpoint:   ',A)
424 FORMAT (/'    Characteristic levels of the geo. wind component vg:'// &
            '       Height:      ',A,'  m'/ &
            '       vg:          ',A,'  m/s'/ &
            '       Gradient:    ',A,'  1/100s'/ &
            '       Gridpoint:   ',A)
425 FORMAT (/'    Characteristic levels of the initial salinity profile:'// &
            '       Height:     ',A,'  m'/ &
            '       Salinity:   ',A,'  psu'/ &
            '       Gradient:   ',A,'  psu/100m'/ &
            '       Gridpoint:  ',A)
426 FORMAT (/'    Characteristic levels of the subsidence/ascent profile:'// &
            '       Height:      ',A,'  m'/ &
            '       w_subs:      ',A,'  m/s'/ &
            '       Gradient:    ',A,'  (m/s)/100m'/ &
            '       Gridpoint:   ',A)
427 FORMAT (/'    Initial wind profiles (u,v) are interpolated from given'// &
                  ' profiles')
428 FORMAT (/'    Initial profiles (u, v, pt, q) are taken from file '/ &
             '    NUDGING_DATA')
430 FORMAT (//' Cloud physics quantities / methods:'/ &
              ' ----------------------------------'/)
431 FORMAT ('    Humidity is considered, bu no condensation')
432 FORMAT ('    Bulk scheme with liquid water potential temperature and'/ &
            '    total water content is used.'/ &
            '    Condensation is parameterized via 0% - or 100% scheme.')
433 FORMAT ('    Cloud droplets treated explicitly using the Lagrangian part', &
                 'icle model')
434 FORMAT ('    Curvature and solution effecs are considered for growth of', &
                 ' droplets < 1.0E-6 m')
435 FORMAT ('    Droplet collision is handled by ',A,'-kernel')
436 FORMAT ('       Fast kernel with fixed radius- and dissipation classes ', &
                    'are used'/ &
            '          number of radius classes:       ',I3,'    interval ', &
                       '[1.0E-6,2.0E-4] m'/ &
            '          number of dissipation classes:   ',I2,'    interval ', &
                       '[0,1000] cm**2/s**3')
437 FORMAT ('    Droplet collision is switched off')
450 FORMAT (//' LES / Turbulence quantities:'/ &
              ' ---------------------------'/)
451 FORMAT ('    Diffusion coefficients are constant:'/ &
            '    Km = ',F6.2,' m**2/s   Kh = ',F6.2,' m**2/s   Pr = ',F5.2)
453 FORMAT ('    Mixing length is limited to',F5.2,' * z')
454 FORMAT ('    TKE is not allowed to fall below ',E9.2,' (m/s)**2')
455 FORMAT ('    initial TKE is prescribed as ',E9.2,' (m/s)**2')
456 FORMAT ('    Day of the year at model start :   day_init = ',I3             &
            /'    UTC time at model start        :   time_utc_init = ',F7.1' s')
457 FORMAT ('    RANS-mode constants: c_0 = ',F9.5/         &
            '                         c_1 = ',F9.5/         &
            '                         c_2 = ',F9.5/         &
            '                         c_3 = ',F9.5/         &
            '                         c_4 = ',F9.5/         &
            '                         sigma_e    = ',F9.5/  &
            '                         sigma_diss = ',F9.5)
470 FORMAT (//' Actions during the simulation:'/ &
              ' -----------------------------'/)
471 FORMAT ('    Disturbance impulse (u,v) every :   ',F6.2,' s'/            &
            '    Disturbance amplitude           :    ',F5.2, ' m/s'/       &
            '    Lower disturbance level         : ',F8.2,' m (GP ',I4,')'/  &
            '    Upper disturbance level         : ',F8.2,' m (GP ',I4,')')
472 FORMAT ('    Disturbances continued during the run from i/j =',I4, &
                 ' to i/j =',I4)
473 FORMAT ('    Disturbances cease as soon as the disturbance energy exceeds',&
                 F6.3, ' m**2/s**2')
474 FORMAT ('    Random number generator used    : ',A/)
475 FORMAT ('    The surface temperature is increased (or decreased, ', &
                 'respectively, if'/ &
            '    the value is negative) by ',F5.2,' K at the beginning of the',&
                 ' 3D-simulation'/)
476 FORMAT ('    The surface humidity is increased (or decreased, ',&
                 'respectively, if the'/ &
            '    value is negative) by ',E8.1,' kg/kg at the beginning of', &
                 ' the 3D-simulation'/)
477 FORMAT ('    The scalar value is increased at the surface (or decreased, ',&
                 'respectively, if the'/ &
            '    value is negative) by ',E8.1,' kg/m**3 at the beginning of', &
                 ' the 3D-simulation'/)
480 FORMAT ('    Particles:'/ &
            '    ---------'// &
            '       Particle advection is active (switched on at t = ', F7.1, &
                    ' s)'/ &
            '       Start of new particle generations every  ',F6.1,' s'/ &
            '       Boundary conditions: left/right: ', A, ' north/south: ', A/&
            '                            bottom:     ', A, ' top:         ', A/&
            '       Maximum particle age:                 ',F9.1,' s'/ &
            '       Advection stopped at t = ',F9.1,' s'/)
481 FORMAT ('       Particles have random start positions'/)
482 FORMAT ('          Particles are advected only horizontally'/)
485 FORMAT ('       Particle data are written on file every ', F9.1, ' s')
486 FORMAT ('       Particle statistics are written on file'/)
487 FORMAT ('       Number of particle groups: ',I2/)
488 FORMAT ('       SGS velocity components are used for particle advection'/ &
            '          minimum timestep for advection:', F8.5/)
489 FORMAT ('       Number of particles simultaneously released at each ', &
                    'point: ', I5/)
490 FORMAT ('       Particle group ',I2,':'/ &
            '          Particle radius: ',E10.3, 'm')
491 FORMAT ('          Particle inertia is activated'/ &
            '             density_ratio (rho_fluid/rho_particle) =',F6.3/)
492 FORMAT ('          Particles are advected only passively (no inertia)'/)
493 FORMAT ('          Boundaries of particle source: x:',F8.1,' - ',F8.1,' m'/&
            '                                         y:',F8.1,' - ',F8.1,' m'/&
            '                                         z:',F8.1,' - ',F8.1,' m'/&
            '          Particle distances:  dx = ',F8.1,' m  dy = ',F8.1, &
                       ' m  dz = ',F8.1,' m'/)
494 FORMAT ('       Output of particle time series in NetCDF format every ', &
                    F8.2,' s'/)
495 FORMAT ('       Number of particles in total domain: ',I10/)
496 FORMAT ('       Initial vertical particle positions are interpreted ', &
                    'as relative to the given topography')
500 FORMAT (//' 1D-Model parameters:'/                           &
              ' -------------------'//                           &
            '    Simulation time:                   ',F8.1,' s'/ &
            '    Run-controll output every:         ',F8.1,' s'/ &
            '    Vertical profile output every:     ',F8.1,' s'/ &
            '    Mixing length calculation:         ',A/         &
            '    Dissipation calculation:           ',A/)
502 FORMAT ('    Damping layer starts from ',F7.1,' m (GP ',I4,')'/)
503 FORMAT (' --> Momentum advection via Wicker-Skamarock-Scheme 5th order')
504 FORMAT (' --> Scalar advection via Wicker-Skamarock-Scheme 5th order')
505 FORMAT ('    Precipitation parameterization via Seifert-Beheng-Scheme')
506 FORMAT ('    Cloud water sedimentation parameterization via Stokes law')
507 FORMAT ('    Turbulence effects on precipitation process')
508 FORMAT ('    Ventilation effects on evaporation of rain drops')
509 FORMAT ('    Slope limiter used for sedimentation process')
510 FORMAT ('    Droplet density    :   N_c   = ',F6.1,' 1/cm**3')
511 FORMAT ('    Sedimentation Courant number:                  '/&
            '                               C_s   =',F4.1,'        ')
512 FORMAT (/' Date:                 ',A8,6X,'Run:       ',A20/      &
            ' Time:                 ',A8,6X,'Run-No.:   ',I2.2/     &
            ' Run on host:        ',A10,6X,'En-No.:    ',I2.2)
600 FORMAT (/' Nesting informations:'/ &
            ' --------------------'/ &
            ' Nesting mode:                     ',A/ &
            ' Nesting-datatransfer mode:        ',A// &
            ' Nest id  parent  number   lower left coordinates   name'/ &
            ' (*=me)     id    of PEs      x (m)     y (m)' )
601 FORMAT (2X,A1,1X,I2.2,6X,I2.2,5X,I5,5X,F8.2,2X,F8.2,5X,A)

 END SUBROUTINE header
