!> @file parin.f90
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
! $Id: parin.f90 3083 2018-06-19 14:03:12Z gronemeier $
! Added rans_const_c and rans_const_sigma as input parameters (TG)
! 
! 3065 2018-06-12 07:03:02Z Giersch
! New initialization parameters added
! 
! 3049 2018-05-29 13:52:36Z Giersch
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! z_max_do2d removed, error messages revised
! 
! 2995 2018-04-19 12:13:16Z Giersch
! time_since_reference_point must be calculated/initialized before the first  
! call of functions related to the radiation model which occur in 
! time_integration_spinup or time_integration
! 
! 2980 2018-04-17 15:19:27Z suehring
! Revise message call
! 
! 2975 2018-04-16 15:22:20Z suehring
! - Informative message when initializing_actions has been changed 
!   to set_constant_profile in child domain 
! - Change location in message call
! 
! 2967 2018-04-13 11:22:08Z raasch
! bugfix: missing parallel cpp-directives added
! 
! 2941 2018-04-03 11:54:58Z kanani
! Fix for spinup in case of restart run
! 
! 2938 2018-03-27 15:52:42Z suehring
! Change initialization in case child domain should be initialized with Inifor.
! 
! 2936 2018-03-27 14:49:27Z suehring
! inipar renamed to initialization_parameters.
! d3par renamed to runtime_parameters.
! 
! 2921 2018-03-22 15:05:23Z Giersch
! Activation of spinup has been moved from lsm/usm_parin to parin itself
! 
! 2906 2018-03-19 08:56:40Z Giersch
! ENVIRONMENT variables read/write_svf has been added 
! 
! 2894 2018-03-15 09:17:58Z Giersch
! read_var_list has been renamed to rrd_global, all module related _parin 
! routines are called before reading the global restart data to overwrite them
! in case of restart runs
! 
! 2881 2018-03-13 16:24:40Z suehring
! Added flag for switching on/off calculation of soil moisture
! 
! 2849 2018-03-05 10:49:33Z Giersch
! Position of d3par namelist in parameter file is unimportant now 
! 
! 2826 2018-02-21 12:39:28Z Giersch
! Bugfix in setting the default boundary conditions for nest domains
!
! 2817 2018-02-19 16:32:21Z knoop
! Preliminary gust module interface implemented
!
! 2773 2018-01-30 14:12:54Z suehring
! Nesting for chemical species implemented
! 
! 2766 2018-01-22 17:17:47Z kanani
! Removed preprocessor directive __chem
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Implementation of uv exposure model (FK)
! Added rans_mode and turbulence_closure to inipar (TG)
! Implementation of chemistry module
! Sorting of USE list (FK)
! Forcing implemented, and initialization with inifor (MS)
! 
! 2600 2017-11-01 14:11:20Z raasch
! some comments added and variables renamed concerning r2599
! 
! 2599 2017-11-01 13:18:45Z hellstea
! The i/o grouping is updated to work correctly also in nested runs.
! 
! 2575 2017-10-24 09:57:58Z maronga
! Renamed phi -> latitude, added longitude
! 
! 2563 2017-10-19 15:36:10Z Giersch
! Changed position where restart files are closed.
! 
! 2550 2017-10-16 17:12:01Z boeske
! Added complex_terrain
! 
! 2544 2017-10-13 18:09:32Z maronga
! Moved day_of_year_init and time_utc_init to inipar.
! 
! 2397 2017-09-04 16:22:48Z suehring
! Enable initialization of 3d model by user in the child domain. 
! 
! 2375 2017-08-29 14:10:28Z schwenkel
! Added aerosol initialization for bulk microphysics
! 
! 2372 2017-08-25 12:37:32Z sward
! y_shift added to namelist
! 
! 2365 2017-08-21 14:59:59Z kanani
! Vertical grid nesting: add vnest_start_time to d3par (SadiqHuq)
! 
! 2339 2017-08-07 13:55:26Z gronemeier
! corrected timestamp in header
! 
! 2338 2017-08-07 12:15:38Z gronemeier
! Modularize 1D model
! 
! 2310 2017-07-11 09:37:02Z gronemeier
! Bugfix: re-arranged call for error messages for ENVPAR file
! 
! 2304 2017-07-04 14:35:55Z suehring
! Bugfix, enable restarts for child domain.
! 
! 2298 2017-06-29 09:28:18Z raasch
! -return_addres, return_username in ENVPAR, -cross_ts_uymax, cross_ts_uymin in
! d3par
! 
! 2296 2017-06-28 07:53:56Z maronga
! Added parameters for model spinup
! 
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison' 
! includes two more prognostic equations for cloud drop concentration (nc)  
! and cloud water content (qc). 
! 
! 2267 2017-06-09 09:33:25Z gronemeier
! Bugfix: removed skipping of reading namelists in case of omitted d3par
! 
! 2259 2017-06-08 09:09:11Z gronemeier
! Implemented synthetic turbulence generator
!
! 2233 2017-05-30 18:08:54Z suehring
! 
! 2232 2017-05-30 17:47:52Z suehring
! typo corrected
! +wall_salinityflux
! +tunnel_height, tunnel_lenght, tunnel_width_x, tunnel_width_y, 
!  tunnel_wall_depth
! 
! 2118 2017-01-17 16:38:49Z raasch
! -background_communication from inipar
! 
! 2050 2016-11-08 15:00:55Z gronemeier
! Implement turbulent outflow condition
! 
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
! 
! 2035 2016-10-24 15:06:17Z suehring
! Remove check for npex and npey in nesting case
! 
! 2011 2016-09-19 17:29:57Z kanani
! Added flag lsf_exception to allow explicit enabling of large scale forcing 
! together with buildings on flat terrain.
! 
! 2007 2016-08-24 15:47:17Z kanani
! Added call to urban surface model for reading of &urban_surface_par
! 
! 2004 2016-08-24 10:25:59Z suehring
! Humidity and passive scalar treated separately in nesting mode
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1992 2016-08-12 15:14:59Z suehring
! +top_scalarflux
! 
! 1960 2016-07-12 16:34:24Z suehring
! Allocation of s_init
! 
! 1957 2016-07-07 10:43:48Z suehring
! flight module added
! 
! 1955 2016-07-01 12:38:59Z hellstea
! The parameter intializating_actions is set to 'set_constant_profiles for 
! all nest domains in order to make sure that diagnostic variables are properly 
! initialized for nest domains. Prognostic variables are later initialized by 
! interpolation from the parent domain. 
!
! 1917 2016-05-27 14:28:12Z witha
! Initial version of purely vertical nesting introduced.
! 
! 1914 2016-05-26 14:44:07Z witha
! Added call to wind turbine model for reading of &wind_turbine_par
!
! 1849 2016-04-08 11:33:18Z hoffmann
! Adapted for modularization of microphysics
!
! 1833 2016-04-07 14:23:03Z raasch
! call of spectra_parin
! 
! 1831 2016-04-07 13:15:51Z hoffmann
! turbulence renamed collision_turbulence, drizzle renamed
! cloud_water_sedimentation
! curvature_solution_effects removed
!
! 1826 2016-04-07 12:01:39Z maronga
! Added call to radiation model for reading of &radiation_par.
! Added call to plant canopy model for reading of &canopy_par.
! 
! 1817 2016-04-06 15:44:20Z maronga
! Added call to land surface model for reading of &lsm_par
! 
! 1804 2016-04-05 16:30:18Z maronga
! Removed code for parameter file check (__check)
!
! 1783 2016-03-06 18:36:17Z raasch
! +netcdf_deflate in d3par, netcdf module and variable names changed
!
! 1764 2016-02-28 12:45:19Z raasch
! cpp-statements for nesting removed, explicit settings of boundary conditions
! in nest domains,
! bugfix: npex/npey message moved from inipar to d3par
! bugfix: check of lateral boundary conditions from check_parameters to here,
! because they will be already used in init_pegrid and init_grid
!
! 1762 2016-02-25 12:31:13Z hellstea
! Introduction of nested domain feature
!
! 1691 2015-10-26 16:17:44Z maronga
! Added parameter most_method. Renamed prandtl_layer to constant_flux_layer.
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1560 2015-03-06 10:48:54Z keck
! +recycling_yshift
! 
! 1496 2014-12-02 17:25:50Z maronga
! Renamed: "radiation -> "cloud_top_radiation"
! 
! 1484 2014-10-21 10:53:05Z kanani
! Changes due to new module structure of the plant canopy model:
!   canopy-model related parameters moved to new package canopy_par in 
!   subroutine package_parin
! 
! 1429 2014-07-15 12:53:45Z knoop
! +ensemble_member_nr to prepare the random_generator for ensemble runs
! 
! 1402 2014-05-09 14:25:13Z raasch
! location messages modified, batch_job included in envpar-NAMELIST
! 
! 1384 2014-05-02 14:31:06Z raasch
! location messages added
! 
! 1365 2014-04-22 15:03:56Z boeske
! Usage of large scale forcing enabled:
! +use_subsidence_tendencies
! 
! 1361 2014-04-16 15:17:48Z hoffmann
! +call_microphysics_at_all_substeps
! 
! 1359 2014-04-11 17:15:14Z hoffmann
! REAL constants provided with KIND-attribute 
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute 
! 
! 1327 2014-03-21 11:00:16Z raasch
! -data_output_format, do3d_compress, do3d_comp_prec
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
! 1318 2014-03-17 13:35:16Z raasch
! +cpu_log_barrierwait in d3par
!
! 1301 2014-03-06 13:29:46Z heinze
! +large_scale_subsidence
!
! 1241 2013-10-30 11:36:58Z heinze
! +nudging
! +large_scale_forcing
!
! 1216 2013-08-26 09:31:42Z raasch
! +transpose_compute_overlap in inipar
!
! 1195 2013-07-01 12:27:57Z heinze
! Bugfix: allocate ref_state
!
! 1179 2013-06-14 05:57:58Z raasch
! +reference_state in inipar
!
! 1159 2013-05-21 11:58:22Z fricke
! +use_cmax
!
! 1128 2013-04-12 06:19:32Z raasch
! +background_communication in inipar
!
! 1115 2013-03-26 18:16:16Z hoffmann
! unused variables removed
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1065 2012-11-22 17:42:36Z hoffmann
! +nc, c_sedimentation, limiter_sedimentation, turbulence
! -mu_constant, mu_constant_value
!
! 1053 2012-11-13 17:11:03Z hoffmann
! necessary expansions according to the two new prognostic equations (nr, qr) 
! of the two-moment cloud physics scheme and steering parameters:
! +*_init, *_surface, *_surface_initial_change, *_vertical_gradient, 
! +*_vertical_gradient_level, surface_waterflux_*, 
! +cloud_scheme, drizzle, mu_constant, mu_constant_value, ventilation_effect
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1015 2012-09-27 09:23:24Z raasch
! -adjust_mixing_length
!
! 1003 2012-09-14 14:35:53Z raasch
! -grid_matching
!
! 1001 2012-09-13 14:08:46Z raasch
! -cut_spline_overshoot, long_filter_factor, overshoot_limit_*, ups_limit_*
!
! 996 2012-09-07 10:41:47Z raasch
! -use_prior_plot1d_parameters
!
! 978 2012-08-09 08:28:32Z fricke
! -km_damp_max, outflow_damping_width
! +pt_damping_factor, pt_damping_width
! +z0h_factor
!
! 964 2012-07-26 09:14:24Z raasch
! -cross_normalized_x, cross_normalized_y, cross_xtext, z_max_do1d,
! z_max_do1d_normalized
!
! 940 2012-07-09 14:31:00Z raasch
! +neutral in inipar
!
! 927 2012-06-06 19:15:04Z raasch
! +masking_method in inipar
!
! 824 2012-02-17 09:09:57Z raasch
! +curvature_solution_effects in inipar
!
! 809 2012-01-30 13:32:58Z maronga
! Bugfix: replaced .AND. and .NOT. with && and ! in the preprocessor directives
!
! 807 2012-01-25 11:53:51Z maronga
! New cpp directive "__check" implemented which is used by check_namelist_files
!
! Revision 1.1  1997/07/24 11:22:50  raasch
! Initial revision
!
!
! Description:
! ------------
!> This subroutine reads variables controling the run from the NAMELIST files
!------------------------------------------------------------------------------!
 SUBROUTINE parin
 

    USE arrays_3d,                                                             &
        ONLY:  pt_init, q_init, ref_state, s_init, sa_init,                    &     
               ug, u_init, v_init, vg

    USE chemistry_model_mod,                                                   &
        ONLY:  chem_parin
        
    USE chem_modules

    USE control_parameters

    USE cpulog,                                                                &
        ONLY:  cpu_log_barrierwait

    USE date_and_time_mod,                                                     &
        ONLY:  day_of_year_init, time_utc_init

    USE dvrp_variables,                                                        &
        ONLY:  local_dvrserver_running

    USE flight_mod,                                                            &
        ONLY:  flight_parin

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE gust_mod,                                                              &
        ONLY: gust_parin

    USE indices,                                                               &
        ONLY:  nx, ny, nz

    USE kinds

    USE land_surface_model_mod,                                                &
        ONLY: lsm_parin

    USE microphysics_mod,                                                      &
        ONLY:  c_sedimentation, cloud_water_sedimentation,                     &
               collision_turbulence, curvature_solution_effects_bulk,          &
               dry_aerosol_radius, limiter_sedimentation, na_init, nc_const,   &
               sigma_bulk, ventilation_effect

    USE model_1d_mod,                                                          &
        ONLY:  damp_level_1d, dt_pr_1d, dt_run_control_1d, end_time_1d

    USE netcdf_interface,                                                      &
        ONLY:  netcdf_data_format, netcdf_deflate, netcdf_precision

    USE pegrid
               
    USE plant_canopy_model_mod,                                                &
         ONLY: pcm_parin

    USE pmc_interface,                                                         &
        ONLY:  nested_run, nesting_mode

    USE profil_parameter,                                                      &
        ONLY:  cross_profiles, profile_columns, profile_rows

    USE progress_bar,                                                          &
        ONLY :  batch_job

    USE radiation_model_mod,                                                   &
        ONLY: radiation_parin

    USE read_restart_data_mod,                                                 &
        ONLY:  rrd_global     

    USE spectra_mod,                                                           &
        ONLY :  spectra_parin

    USE statistics,                                                            &
        ONLY:  hom, hom_sum, pr_palm, region, statistic_regions

    USE synthetic_turbulence_generator_mod,                                    &
        ONLY:  stg_parin

    USE turbulence_closure_mod,                                                &
        ONLY:  rans_const_c, rans_const_sigma

    USE urban_surface_mod,                                                     &
        ONLY: usm_parin

    USE uv_exposure_model_mod,                                                 &
        ONLY:  uvem_parin

    USE vertical_nesting_mod,                                                  &
        ONLY:  vnest_start_time

    USE wind_turbine_model_mod,                                                &
        ONLY:  wtm_parin


    IMPLICIT NONE

    CHARACTER (LEN=80) ::  line  !< dummy string that contains the current line of the parameter file 

    INTEGER(iwp) ::  global_id      !< process id with respect to MPI_COMM_WORLD
    INTEGER(iwp) ::  global_procs   !< # of procs with respect to MPI_COMM_WORLD
    INTEGER(iwp) ::  i              !<
    INTEGER(iwp) ::  ioerr          !< error flag for open/read/write

    NAMELIST /inipar/  aerosol_bulk, alpha_surface, approximation, bc_e_b,     &
                       bc_lr, bc_ns, bc_p_b, bc_p_t, bc_pt_b, bc_pt_t, bc_q_b, &
             bc_q_t,bc_s_b, bc_s_t, bc_sa_t, bc_uv_b, bc_uv_t,                 &
             bottom_salinityflux, building_height, building_length_x,          &
             building_length_y, building_wall_left, building_wall_south,       &
             calc_soil_moisture_during_spinup,                                 &
             call_psolver_at_all_substeps, call_microphysics_at_all_substeps,  &
             canyon_height,                                                    &
             canyon_width_x, canyon_width_y, canyon_wall_left,                 &
             canyon_wall_south, c_sedimentation, cfl_factor, cloud_droplets,   &
             cloud_physics, cloud_scheme, cloud_top_radiation,                 &
             cloud_water_sedimentation,                                        &
             collective_wait, collision_turbulence, complex_terrain,           &
             conserve_volume_flow,                                             &
             conserve_volume_flow_mode, constant_flux_layer,                   &
             coupling_start_time, curvature_solution_effects_bulk,             &
             cycle_mg, damp_level_1d,                                          &
             data_output_during_spinup,                                        &
             day_of_year_init,                                                 &
             dissipation_1d,                                                   &
             dp_external, dp_level_b, dp_smooth, dpdxy, dry_aerosol_radius,    &
             dt, dt_pr_1d, dt_run_control_1d, dt_spinup, dx, dy, dz, dz_max,   &
             dz_stretch_factor, dz_stretch_level, dz_stretch_level_start,      &
             dz_stretch_level_end, end_time_1d, ensemble_member_nr, e_init,    &
             e_min, fft_method, flux_input_mode, flux_output_mode, forcing,    &
             galilei_transformation, humidity,                                 &
             inflow_damping_height, inflow_damping_width,                      &
             inflow_disturbance_begin, inflow_disturbance_end,                 &
             initializing_actions, km_constant,                                &
             large_scale_forcing, large_scale_subsidence, latitude,            &
             limiter_sedimentation, longitude,                                 &
             loop_optimization, lsf_exception, masking_method, mg_cycles,      &
             mg_switch_to_pe0_level, mixing_length_1d, momentum_advec,         &
             most_method, na_init, nc_const, netcdf_precision, neutral, ngsrb, &
             nsor, nsor_ini, nudging, nx, ny, nz, ocean, idealized_diurnal,    &
             linear_eqnOfState, rho_ref, fixed_alpha, alpha_const,             &
             beta_const, pt_ref, sa_ref,                                       & 
             omega, omega_sor, wb_solar,     &
             outflow_source_plane, passive_scalar,                             &
             prandtl_number, precipitation, psolver, pt_damping_factor,        &
             pt_damping_width, pt_reference, pt_surface,                       &
             pt_surface_initial_change, pt_vertical_gradient,                  &
             pt_vertical_gradient_level, q_surface, q_surface_initial_change,  &
             q_vertical_gradient, q_vertical_gradient_level,                   &
             random_generator, random_heatflux, rans_const_c, rans_const_sigma,&
             rans_mode,                                                        &
             rayleigh_damping_factor, rayleigh_damping_height,                 &
             recycling_width, recycling_yshift,                                &
             reference_state, residual_limit,                                  &
             roughness_length, sa_surface,                                     &
             sa_vertical_gradient, sa_vertical_gradient_level, scalar_advec,   &
             scalar_rayleigh_damping, sigma_bulk,                              &
             spinup_time, spinup_pt_amplitude, spinup_pt_mean,                 &
             statistic_regions, subs_vertical_gradient,                        &
             subs_vertical_gradient_level, surface_heatflux, surface_pressure, &
             surface_scalarflux, surface_waterflux,                            &
             s_surface, s_surface_initial_change, s_vertical_gradient,         &
             s_vertical_gradient_level, time_utc_init, timestep_scheme,        &
             topography, topography_grid_convention, top_heatflux,             &
             top_momentumflux_u, top_momentumflux_v, top_salinityflux,         &
             top_scalarflux, transpose_compute_overlap,                        &
             tunnel_height, tunnel_length, tunnel_width_x, tunnel_width_y,     &
             tunnel_wall_depth, turbulence_closure,                            &
             turbulent_inflow, turbulent_outflow,                              &
             use_subsidence_tendencies, ug_surface, ug_vertical_gradient,      &
             ug_vertical_gradient_level, use_surface_fluxes, use_cmax,         &
             use_top_fluxes, use_ug_for_galilei_tr, use_upstream_for_tke,      &
             uv_heights, u_bulk, u_profile, vg_surface, vg_vertical_gradient,  &
             vg_vertical_gradient_level, v_bulk, v_profile, ventilation_effect,&
             wall_adjustment, wall_heatflux, wall_humidityflux,                &
             wall_salinityflux, wall_scalarflux, y_shift, zeta_max, zeta_min,  &
             z0h_factor

    NAMELIST /initialization_parameters/  aerosol_bulk, alpha_surface,         &
             approximation, bc_e_b,                                            &
             bc_lr, bc_ns, bc_p_b, bc_p_t, bc_pt_b, bc_pt_t, bc_q_b,           &
             bc_q_t,bc_s_b, bc_s_t, bc_sa_t, bc_uv_b, bc_uv_t,                 &
             bottom_salinityflux, building_height, building_length_x,          &
             building_length_y, building_wall_left, building_wall_south,       &
             calc_soil_moisture_during_spinup,                                 &
             call_psolver_at_all_substeps, call_microphysics_at_all_substeps,  &
             canyon_height,                                                    &
             canyon_width_x, canyon_width_y, canyon_wall_left,                 &
             canyon_wall_south, c_sedimentation, cfl_factor, cloud_droplets,   &
             cloud_physics, cloud_scheme, cloud_top_radiation,                 &
             cloud_water_sedimentation,                                        &
             collective_wait, collision_turbulence, complex_terrain,           &
             conserve_volume_flow,                                             &
             conserve_volume_flow_mode, constant_flux_layer,                   &
             coupling_start_time, curvature_solution_effects_bulk,             &
             cycle_mg, damp_level_1d,                                          &
             data_output_during_spinup,                                        &
             day_of_year_init,                                                 &
             dissipation_1d,                                                   &
             dp_external, dp_level_b, dp_smooth, dpdxy, dry_aerosol_radius,    &
             dt, dt_pr_1d, dt_run_control_1d, dt_spinup, dx, dy, dz, dz_max,   &
             dz_stretch_factor, dz_stretch_level, dz_stretch_level_start,      &
             dz_stretch_level_end, end_time_1d, ensemble_member_nr, e_init,    &
             e_min, fft_method, flux_input_mode, flux_output_mode, forcing,    &
             galilei_transformation, humidity,                                 &
             inflow_damping_height, inflow_damping_width,                      &
             inflow_disturbance_begin, inflow_disturbance_end,                 &
             initializing_actions, km_constant,                                &
             large_scale_forcing, large_scale_subsidence, latitude,            &
             limiter_sedimentation, longitude,                                 &
             loop_optimization, lsf_exception, masking_method, mg_cycles,      &
             mg_switch_to_pe0_level, mixing_length_1d, momentum_advec,         &
             most_method, na_init, nc_const, netcdf_precision, neutral, ngsrb, &
             nsor, nsor_ini, nudging, nx, ny, nz, ocean, idealized_diurnal,    &
             linear_eqnOfState, rho_ref, fixed_alpha, alpha_const,             &
             beta_const, pt_ref, sa_ref,                                       &
             omega, omega_sor, wb_solar,    &
             outflow_source_plane, passive_scalar,                             &
             prandtl_number, precipitation, psolver, pt_damping_factor,        &
             pt_damping_width, pt_reference, pt_surface,                       &
             pt_surface_initial_change, pt_vertical_gradient,                  &
             pt_vertical_gradient_level, q_surface, q_surface_initial_change,  &
             q_vertical_gradient, q_vertical_gradient_level,                   &
             random_generator, random_heatflux, rans_const_c, rans_const_sigma,&
             rans_mode,                                                        &
             rayleigh_damping_factor, rayleigh_damping_height,                 &
             recycling_width, recycling_yshift,                                &
             reference_state, residual_limit,                                  &
             roughness_length, sa_surface,                                     &
             sa_vertical_gradient, sa_vertical_gradient_level, scalar_advec,   &
             scalar_rayleigh_damping, sigma_bulk,                              &
             spinup_time, spinup_pt_amplitude, spinup_pt_mean,                 &
             statistic_regions, subs_vertical_gradient,                        &
             subs_vertical_gradient_level, surface_heatflux, surface_pressure, &
             surface_scalarflux, surface_waterflux,                            &
             s_surface, s_surface_initial_change, s_vertical_gradient,         &
             s_vertical_gradient_level, time_utc_init, timestep_scheme,        &
             topography, topography_grid_convention, top_heatflux,             &
             top_momentumflux_u, top_momentumflux_v, top_salinityflux,         &
             top_scalarflux, transpose_compute_overlap,                        &
             tunnel_height, tunnel_length, tunnel_width_x, tunnel_width_y,     &
             tunnel_wall_depth, turbulence_closure,                            &
             turbulent_inflow, turbulent_outflow,                              &
             use_subsidence_tendencies, ug_surface, ug_vertical_gradient,      &
             ug_vertical_gradient_level, use_surface_fluxes, use_cmax,         &
             use_top_fluxes, use_ug_for_galilei_tr, use_upstream_for_tke,      &
             uv_heights, u_bulk, u_profile, vg_surface, vg_vertical_gradient,  &
             vg_vertical_gradient_level, v_bulk, v_profile, ventilation_effect,&
             wall_adjustment, wall_heatflux, wall_humidityflux,                &
             wall_salinityflux, wall_scalarflux, y_shift, zeta_max, zeta_min,  &
             z0h_factor
             
    NAMELIST /d3par/  averaging_interval, averaging_interval_pr,               &
             cpu_log_barrierwait, create_disturbances,                         &
             cross_profiles, data_output, data_output_masks,                   &
             data_output_pr, data_output_2d_on_each_pe, disturbance_amplitude, &
             disturbance_energy_limit, disturbance_level_b,                    &
             disturbance_level_t, do2d_at_begin, do3d_at_begin,                &
             dt, dt_averaging_input, dt_averaging_input_pr,                    &
             dt_coupling, dt_data_output, dt_data_output_av, dt_disturb,       &
             dt_domask, dt_dopr, dt_dopr_listing, dt_dots, dt_do2d_xy,         &
             dt_do2d_xz, dt_do2d_yz, dt_do3d, dt_max, dt_restart,              &
             dt_run_control,end_time, force_print_header, mask_scale_x,        &
             mask_scale_y, mask_scale_z, mask_x, mask_y, mask_z, mask_x_loop,  &
             mask_y_loop, mask_z_loop, netcdf_data_format, netcdf_deflate,     &
             normalizing_region, npex, npey, nz_do3d,                          &
             precipitation_amount_interval, profile_columns, profile_rows,     &
             restart_time, section_xy, section_xz, section_yz,                 &
             skip_time_data_output, skip_time_data_output_av, skip_time_dopr,  &
             skip_time_do2d_xy, skip_time_do2d_xz, skip_time_do2d_yz,          &
             skip_time_do3d, skip_time_domask, synchronous_exchange,           &
             termination_time_needed, vnest_start_time

    NAMELIST /runtime_parameters/  averaging_interval, averaging_interval_pr,  &
             cpu_log_barrierwait, create_disturbances,                         &
             cross_profiles, data_output, data_output_masks,                   &
             data_output_pr, data_output_2d_on_each_pe, disturbance_amplitude, &
             disturbance_energy_limit, disturbance_level_b,                    &
             disturbance_level_t, do2d_at_begin, do3d_at_begin,                &
             dt, dt_averaging_input, dt_averaging_input_pr,                    &
             dt_coupling, dt_data_output, dt_data_output_av, dt_disturb,       &
             dt_domask, dt_dopr, dt_dopr_listing, dt_dots, dt_do2d_xy,         &
             dt_do2d_xz, dt_do2d_yz, dt_do3d, dt_max, dt_restart,              &
             dt_run_control,end_time, force_print_header, mask_scale_x,        &
             mask_scale_y, mask_scale_z, mask_x, mask_y, mask_z, mask_x_loop,  &
             mask_y_loop, mask_z_loop, netcdf_data_format, netcdf_deflate,     &
             normalizing_region, npex, npey, nz_do3d,                          &
             precipitation_amount_interval, profile_columns, profile_rows,     &
             restart_time, section_xy, section_xz, section_yz,                 &
             skip_time_data_output, skip_time_data_output_av, skip_time_dopr,  &
             skip_time_do2d_xy, skip_time_do2d_xz, skip_time_do2d_yz,          &
             skip_time_do3d, skip_time_domask, synchronous_exchange,           &
             termination_time_needed, vnest_start_time

    NAMELIST /envpar/  batch_job, host, local_dvrserver_running,               &
                       maximum_cpu_time_allowed, maximum_parallel_io_streams,  &
                       read_svf, revision, run_identifier, tasks_per_node,     &
                       write_binary, write_svf

!
!-- First read values of environment variables (this NAMELIST file is
!-- generated by palmrun)
    CALL location_message( 'reading environment parameters from ENVPAR', .FALSE. )

    OPEN ( 90, FILE='ENVPAR', STATUS='OLD', FORM='FORMATTED', IOSTAT=ioerr )

    IF ( ioerr /= 0 )  THEN
       message_string = 'local file ENVPAR not found' //                       &
                        '&some variables for steering may not be properly set'
       CALL message( 'parin', 'PA0276', 0, 1, 0, 6, 0 )
    ELSE
       READ ( 90, envpar, IOSTAT=ioerr )
       IF ( ioerr < 0 )  THEN
          message_string = 'no envpar-NAMELIST found in local file '  //       &
                           'ENVPAR& or some variables for steering may '  //   &
                           'not be properly set'
          CALL message( 'parin', 'PA0278', 0, 1, 0, 6, 0 )
       ELSEIF ( ioerr > 0 )  THEN
          message_string = 'errors in local file ENVPAR' //                    &
                           '&some variables for steering may not be properly set'
          CALL message( 'parin', 'PA0277', 0, 1, 0, 6, 0 )
       ENDIF
       CLOSE ( 90 )
    ENDIF

    CALL location_message( 'finished', .TRUE. )
!
!-- Calculate the number of groups into which parallel I/O is split.
!-- The default for files which are opened by all PEs (or where each
!-- PE opens his own independent file) is, that all PEs are doing input/output
!-- in parallel at the same time. This might cause performance or even more
!-- severe problems depending on the configuration of the underlying file
!-- system.
!-- Calculation of the number of blocks and the I/O group must be based on all
!-- PEs involved in this run. Since myid and numprocs are related to the
!-- comm2d communicator, which gives only a subset of all PEs in case of
!-- nested runs, that information must be inquired again from the global
!-- communicator.
!-- First, set the default:
#if defined( __parallel )
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, global_id, ierr )
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, global_procs, ierr )
#else
    global_id    = 0
    global_procs = 1
#endif
    IF ( maximum_parallel_io_streams == -1  .OR.                               &
         maximum_parallel_io_streams > global_procs )  THEN
       maximum_parallel_io_streams = global_procs
    ENDIF
!
!-- Now calculate the number of io_blocks and the io_group to which the
!-- respective PE belongs. I/O of the groups is done in serial, but in parallel
!-- for all PEs belonging to the same group.
    io_blocks = global_procs / maximum_parallel_io_streams
    io_group  = MOD( global_id+1, io_blocks )
    
    CALL location_message( 'reading NAMELIST parameters from PARIN', .FALSE. )
!
!-- Data is read in parallel by groups of PEs
    DO  i = 0, io_blocks-1
       IF ( i == io_group )  THEN

!
!--       Open the NAMELIST-file which is send with this job
          CALL check_open( 11 )

!
!--       Read the control parameters for initialization.
!--       The namelist "inipar" must be provided in the NAMELIST-file.
          READ ( 11, initialization_parameters, ERR=10, END=11 )

          GOTO 12

 10       message_string = 'errors in initialization_parameters & or no ' //  &
                           'initialization_parameters-namelist ' //           &
                           'found (CRAY-machines only)'
          CALL message( 'parin', 'PA0271', 1, 2, 0, 6, 0 )

 11       REWIND ( 11 )
          READ ( 11, inipar, ERR=13, END=14 )
 
          message_string = 'namelist inipar is deprecated and will be ' //    &
                          'removed in near future. & Please use namelist ' // &
                          'initialization_parameters instead'
          CALL message( 'parin', 'PA0017', 0, 1, 0, 6, 0 )
 
          GOTO 12
 
 13       message_string = 'errors in inipar & or no inipar-namelist ' //      &
                           'found (CRAY-machines only)'
          CALL message( 'parin', 'PA0271', 1, 2, 0, 6, 0 )
          
 14       message_string = 'no initialization_parameters-namelist found'
          CALL message( 'parin', 'PA0272', 1, 2, 0, 6, 0 )

!
!--       Try to read runtime parameters given by the user for this run
!--       (namelist "runtime_parameters"). The namelist "runtime_parmeters"    
!--       can be omitted. In that case default values are used for the         
!--       parameters.
 12       line = ' '

          REWIND ( 11 )
          line = ' '
          DO   WHILE ( INDEX( line, '&runtime_parameters' ) == 0 )
             READ ( 11, '(A)', END=20 )  line
          ENDDO
          BACKSPACE ( 11 )

!
!--       Read namelist
          READ ( 11, runtime_parameters )

          GOTO 21
          
 20       REWIND ( 11 )
          line = ' '
          DO   WHILE ( INDEX( line, '&d3par' ) == 0 )
             READ ( 11, '(A)', END=21 )  line
          ENDDO
          BACKSPACE ( 11 )
 
 !
!--       Read namelist
          READ ( 11, d3par )
 
          message_string = 'namelist d3par is deprecated and will be ' //      &
                          'removed in near future. &Please use namelist ' //   &
                          'runtime_parameters instead'
          CALL message( 'parin', 'PA0487', 0, 1, 0, 6, 0 )
          
 21       CONTINUE

!
!--       Check if land surface model is used and read &lsm_par if required
          CALL lsm_parin

!
!--       Check if urban surface model is used and read &urban_surface_par if required
          CALL usm_parin

!
!--       Check if spectra shall be calculated and read spectra_par if required
          CALL spectra_parin

!
!--       Check if radiation model is used and read &radiation_par if required
          CALL radiation_parin

!
!--       Check if gust module is used and read &gust_par if required
          CALL gust_parin
 
 
!--       Check if plant canopy model is used and read &canopy_par if required
          CALL pcm_parin
 
!
!--       Read control parameters for optionally used model software packages
          CALL package_parin

!
!--       Check if wind turbine model is used and read &wind_turbine_par if
!--       required
          CALL wtm_parin
!
!--       Check if virtual flights should be carried out and read &flight_par
!--       if required
          CALL flight_parin
!
!--       Check if synthetic turbulence generator is used and read stg_par if
!--       required
          CALL stg_parin
!
!--       Read chemistry variables
          CALL chem_parin
!
!--       Check if uv exposure model is used and read &uvexposure_par
          CALL uvem_parin
!
!--       Read user-defined variables
          CALL user_parin

!
!--       If required, read control parameters from restart file (produced by
!--       a prior run). All PEs are reading from file created by PE0 (see
!--       check_open)
          IF ( TRIM( initializing_actions ) == 'read_restart_data' )  THEN

             CALL rrd_global
!
!--          Increment the run count
             runnr = runnr + 1
          ENDIF

!
!--       Activate spinup
          IF ( land_surface .OR. urban_surface )  THEN
             IF ( spinup_time > 0.0_wp )  THEN
                coupling_start_time = spinup_time
                time_since_reference_point = simulated_time - coupling_start_time
                IF ( spinup_pt_mean == 9999999.9_wp )  THEN
                   spinup_pt_mean = pt_surface
                ENDIF
                end_time = end_time + spinup_time
                IF ( TRIM( initializing_actions ) /= 'read_restart_data' )     &
                   spinup = .TRUE.
             ENDIF
          ENDIF

!
!--       In case of nested runs, explicitly set nesting boundary conditions. 
!--       This will overwrite the user settings and basic defaults.
!--       bc_lr and bc_ns always need to be cyclic for vertical nesting.
          IF ( nested_run )  THEN
             IF ( nesting_mode == 'vertical' )  THEN
                IF (bc_lr /= 'cyclic' .OR. bc_ns /= 'cyclic' )  THEN
                   WRITE ( message_string, *) 'bc_lr and bc_ns were set to ,', &
                        'cyclic for vertical nesting'
                   CALL message( 'parin', 'PA0428', 0, 0, 0, 6, 0 )
                   bc_lr   = 'cyclic'
                   bc_ns   = 'cyclic'
                ENDIF
                IF ( nest_domain )  THEN
                   bc_uv_t  = 'nested'
                   bc_pt_t  = 'nested'
                   bc_q_t   = 'nested'
                   bc_s_t   = 'nested'
                   bc_cs_t  = 'nested'
                   bc_p_t   = 'neumann'  
                ENDIF
!
!--          For other nesting modes only set boundary conditions for 
!--          nested domains.
             ELSE 
                IF ( nest_domain )  THEN
                   bc_lr    = 'nested'
                   bc_ns    = 'nested'
                   bc_uv_t  = 'nested'
                   bc_pt_t  = 'nested'
                   bc_q_t   = 'nested'
                   bc_s_t   = 'nested'
                   bc_cs_t  = 'nested'
                   bc_p_t   = 'neumann'
                ENDIF
             ENDIF
          ENDIF

          IF ( forcing )  THEN
             bc_lr    = 'forcing'
             bc_ns    = 'forcing'
             bc_uv_t  = 'forcing'
             bc_pt_t  = 'forcing'
             bc_q_t   = 'forcing'
             bc_s_t   = 'forcing'  ! scalar boundary condition is not clear
             bc_cs_t  = 'forcing'  ! same for chemical species
             bc_p_t   = 'neumann'
          ENDIF

!         
!--       In case of nested runs, make sure that initializing_actions =
!--       'set_constant_profiles' even though the constant-profiles 
!--       initializations for the prognostic variables will be overwritten 
!--       by pmci_child_initialize and pmci_parent_initialize. This is, 
!--       however, important e.g. to make sure that diagnostic variables 
!--       are set properly. An exception is made in case of restart runs and
!--       if user decides to do everything by its own.
          IF ( nest_domain  .AND.  .NOT. (                                     &
               TRIM( initializing_actions ) == 'read_restart_data'      .OR.   &
               TRIM( initializing_actions ) == 'set_constant_profiles'  .OR.   &
               TRIM( initializing_actions ) == 'by_user' ) )  THEN
             message_string = 'initializing_actions = ' //                     &
                              TRIM( initializing_actions ) // ' has been ' //  &
                              'changed to set_constant_profiles in child ' //  &
                              'domain.' 
             CALL message( 'parin', 'PA0492', 0, 0, 0, 6, 0 )

             initializing_actions = 'set_constant_profiles'
          ENDIF
            
!
!--       Check validity of lateral boundary conditions. This has to be done
!--       here because they are already used in init_pegrid and init_grid and
!--       therefore cannot be check in check_parameters
          IF ( bc_lr /= 'cyclic'  .AND.  bc_lr /= 'dirichlet/radiation'  .AND. &
               bc_lr /= 'radiation/dirichlet'  .AND.  bc_lr /= 'nested'  .AND. &
               bc_lr /= 'forcing' )  THEN
             message_string = 'unknown boundary condition: bc_lr = "' // &
                              TRIM( bc_lr ) // '"'
             CALL message( 'parin', 'PA0049', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( bc_ns /= 'cyclic'  .AND.  bc_ns /= 'dirichlet/radiation'  .AND. &
               bc_ns /= 'radiation/dirichlet'  .AND.  bc_ns /= 'nested'  .AND. &
               bc_ns /= 'forcing' )  THEN
             message_string = 'unknown boundary condition: bc_ns = "' // &
                              TRIM( bc_ns ) // '"'
             CALL message( 'parin', 'PA0050', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       Set internal variables used for speed optimization in if clauses
          IF ( bc_lr /= 'cyclic' )               bc_lr_cyc    = .FALSE.
          IF ( bc_lr == 'dirichlet/radiation' )  bc_lr_dirrad = .TRUE.
          IF ( bc_lr == 'radiation/dirichlet' )  bc_lr_raddir = .TRUE.
          IF ( bc_ns /= 'cyclic' )               bc_ns_cyc    = .FALSE.
          IF ( bc_ns == 'dirichlet/radiation' )  bc_ns_dirrad = .TRUE.
          IF ( bc_ns == 'radiation/dirichlet' )  bc_ns_raddir = .TRUE.

!
!--       Definition of names of areas used for computing statistics. They must
!--       be defined at this place, because they are allowed to be redefined by
!--       the user in user_parin.
          region = 'total domain'

!
!--       Check in case of initial run, if the grid point numbers are well
!--       defined and allocate some arrays which are already needed in
!--       init_pegrid or check_parameters. During restart jobs, these arrays
!--       will be allocated in rrd_global. All other arrays are allocated
!--       in init_3d_model.
          IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN

             IF ( nx <= 0 )  THEN
                WRITE( message_string, * ) 'no value or wrong value given',    &
                                           ' for nx: nx=', nx
                CALL message( 'parin', 'PA0273', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( ny <= 0 )  THEN
                WRITE( message_string, * ) 'no value or wrong value given',    &
                                           ' for ny: ny=', ny
                CALL message( 'parin', 'PA0274', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( nz <= 0 )  THEN
                WRITE( message_string, * ) 'no value or wrong value given',    &
                                           ' for nz: nz=', nz
                CALL message( 'parin', 'PA0275', 1, 2, 0, 6, 0 )
             ENDIF
!
!--          ATTENTION: in case of changes to the following statement please
!--                  also check the allocate statement in routine rrd_global
             ALLOCATE( pt_init(0:nz+1), q_init(0:nz+1), s_init(0:nz+1),        &
                       ref_state(0:nz+1), sa_init(0:nz+1), ug(0:nz+1),         &
                       u_init(0:nz+1), v_init(0:nz+1), vg(0:nz+1),             &
                       hom(0:nz+1,2,pr_palm+max_pr_user,0:statistic_regions),  &
                       hom_sum(0:nz+1,pr_palm+max_pr_user,0:statistic_regions) )

             hom = 0.0_wp

          ENDIF

!
!--       NAMELIST-file is not needed anymore
          CALL close_file( 11 )

       ENDIF
#if defined( __parallel )
       CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
#endif
    ENDDO

    CALL location_message( 'finished', .TRUE. )

 END SUBROUTINE parin
