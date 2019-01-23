!> @file read_restart_data_mod.f90
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
! $Id: read_restart_data_mod.f90 3065 2018-06-12 07:03:02Z Giersch $
! New parameters concerning vertical grid stretching have been added
! 
! 3056 2018-06-04 07:49:35Z Giersch
! found variable has to be set to false inside overlap loop
! 
! 3049 2018-05-29 13:52:36Z Giersch
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
! 
! 3004 2018-04-27 12:33:25Z Giersch
! precipitation_rate_av removed
! 
! 3003 2018-04-23 10:22:58Z Giersch
! z_i is also read to use the last known inversion height from the 
! initial run as the first inversion height which is written into the
! run control file 
! 
! 2956 2018-04-10 11:01:03Z Giersch
! spectrum_x and spectrum_y have been moved to global data
! 
! 2921 2018-03-22 15:05:23Z Giersch
! spinup_time, day_of_year_init and time_utc_init are also read now
! 
! 2912 2018-03-20 13:00:05Z knoop
! Added gust module interface calls
! 
! 2894 2018-03-15 09:17:58Z Giersch
! Initial revision
! 
!
! Description:
! ------------
!> Reads restart data from restart-file(s) (binary format).
!------------------------------------------------------------------------------!
 MODULE read_restart_data_mod


    USE control_parameters
     

    IMPLICIT NONE


    INTERFACE rrd_global
       MODULE PROCEDURE rrd_global
    END INTERFACE rrd_global

    INTERFACE rrd_read_parts_of_global
       MODULE PROCEDURE rrd_read_parts_of_global
    END INTERFACE rrd_read_parts_of_global

    INTERFACE rrd_local
       MODULE PROCEDURE rrd_local
    END INTERFACE rrd_local

    INTERFACE rrd_skip_global
       MODULE PROCEDURE rrd_skip_global
    END INTERFACE rrd_skip_global


    PUBLIC rrd_global, rrd_read_parts_of_global, rrd_local,      & 
           rrd_skip_global


 CONTAINS

! Description:
! ------------
!> Reads values of global control variables from restart-file (binary format) 
!> created by PE0 of the previous run
!------------------------------------------------------------------------------!
    SUBROUTINE rrd_global


       USE arrays_3d,                                                          &
           ONLY:  inflow_damping_factor, mean_inflow_profiles, pt_init,        &
                  q_init, ref_state, s_init, sa_init, u_init, ug, v_init, vg

       USE date_and_time_mod,                                                  &
           ONLY:  day_of_year_init, time_utc_init

           USE grid_variables,                                                     &
           ONLY:  dx, dy

           USE indices,                                                            &
           ONLY:  nz, nx, nx_on_file, ny, ny_on_file

           USE netcdf_interface,                                                   &
           ONLY:  netcdf_precision, output_for_t0

       USE pegrid

       USE statistics,                                                         &
           ONLY:  statistic_regions, hom, hom_sum, pr_palm, u_max, u_max_ijk,  &
                  v_max, v_max_ijk, w_max, w_max_ijk, z_i

       IMPLICIT NONE

       CHARACTER (LEN=10) ::  binary_version_global, version_on_file

       LOGICAL ::  found  


       CALL check_open( 13 )

!
!-- Make version number check first
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)
       READ ( 13 )  version_on_file
!
!-- Read number of PEs and horizontal index bounds of all PEs used in previous
!-- run
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( TRIM( restart_string(1:length) ) /= 'numprocs' )  THEN
          WRITE( message_string, * ) 'numprocs not found in data from prior ', &
                                     'run on PE ', myid
          CALL message( 'rrd_global', 'PA0297', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  numprocs_previous_run

       IF ( .NOT. ALLOCATED( hor_index_bounds_previous_run ) )  THEN
          ALLOCATE( hor_index_bounds_previous_run(4,0:numprocs_previous_run-1) )
       ENDIF

       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'hor_index_bounds' )  THEN
          WRITE( message_string, * ) 'hor_index_bounds not found in data ',    &
                                     'from prior run on PE ', myid
          CALL message( 'rrd_global', 'PA0298', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  hor_index_bounds_previous_run

!
!-- Read vertical number of gridpoints and number of different areas used
!-- for computing statistics. Allocate arrays depending on these values,
!-- which are needed for the following read instructions.
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'nz' )  THEN
          WRITE( message_string, * ) 'nz not found in data from prior run ',   &
                                     'on PE ', myid
          CALL message( 'rrd_global', 'PA0299', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  nz

       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'max_pr_user' )  THEN
          WRITE( message_string, * ) 'max_pr_user not found in data from ',    &
                                     'prior run on PE ', myid
          CALL message( 'rrd_global', 'PA0300', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  max_pr_user    ! This value is checked against the number of
                                   ! user profiles given for the current run
                                   ! in routine user_parin (it has to match)

       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'statistic_regions' )  THEN
          WRITE( message_string, * ) 'statistic_regions not found in data ',   &
                                     'from prior run on PE ', myid
          CALL message( 'rrd_global', 'PA0301', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  statistic_regions
       IF ( .NOT. ALLOCATED( ug ) )  THEN
          ALLOCATE( ug(0:nz+1), u_init(0:nz+1), vg(0:nz+1),                    &
                    v_init(0:nz+1), pt_init(0:nz+1), q_init(0:nz+1),           &
                    ref_state(0:nz+1), s_init(0:nz+1), sa_init(0:nz+1),        &
                    hom(0:nz+1,2,pr_palm+max_pr_user,0:statistic_regions),     &
                    hom_sum(0:nz+1,pr_palm+max_pr_user,0:statistic_regions) )
       ENDIF

!
!-- Now read all control parameters:
!-- Caution: When the following read instructions have been changed, the
!-- -------  version number stored in the variable binary_version_global has to
!--          be increased. The same changes must also be done in 
!--          wrd_write_global.
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       DO WHILE ( restart_string(1:length) /= 'binary_version_local' )

          found = .FALSE.

          SELECT CASE ( restart_string(1:length) )

             CASE ( 'advected_distance_x' )
                READ ( 13 )  advected_distance_x
             CASE ( 'advected_distance_y' )
                READ ( 13 )  advected_distance_y
             CASE ( 'alpha_surface' )
                READ ( 13 )  alpha_surface
             CASE ( 'average_count_pr' )
                READ ( 13 )  average_count_pr
             CASE ( 'average_count_3d' )
                READ ( 13 )  average_count_3d
             CASE ( 'bc_e_b' )
                READ ( 13 )  bc_e_b
             CASE ( 'bc_lr' )
                READ ( 13 )  bc_lr
             CASE ( 'bc_ns' )
                READ ( 13 )  bc_ns
             CASE ( 'bc_p_b' )
                READ ( 13 )  bc_p_b
             CASE ( 'bc_p_t' )
                READ ( 13 )  bc_p_t
             CASE ( 'bc_pt_b' )
                READ ( 13 )  bc_pt_b
             CASE ( 'bc_pt_t' )
                READ ( 13 )  bc_pt_t
             CASE ( 'bc_pt_t_val' )
                READ ( 13 )  bc_pt_t_val
             CASE ( 'bc_q_b' )
                READ ( 13 )  bc_q_b
             CASE ( 'bc_q_t' )
                READ ( 13 )  bc_q_t
             CASE ( 'bc_q_t_val' )
                READ ( 13 )  bc_q_t_val
             CASE ( 'bc_s_b' )
                READ ( 13 )  bc_s_b
             CASE ( 'bc_s_t' )
                READ ( 13 )  bc_s_t
             CASE ( 'bc_sa_b' )
                READ ( 13 )  bc_sa_b
             CASE ( 'bc_sa_t' )
                READ ( 13 )  bc_sa_t
             CASE ( 'bc_uv_b' )
                READ ( 13 )  bc_uv_b
             CASE ( 'bc_uv_t' )
                READ ( 13 )  bc_uv_t
             CASE ( 'bottom_salinityflux' )
                READ ( 13 )  bottom_salinityflux
             CASE ( 'building_height' )
                READ ( 13 )  building_height
             CASE ( 'building_length_x' )
                READ ( 13 )  building_length_x
             CASE ( 'building_length_y' )
                READ ( 13 )  building_length_y
             CASE ( 'building_wall_left' )
                READ ( 13 )  building_wall_left
             CASE ( 'building_wall_south' )
                READ ( 13 )  building_wall_south
             CASE ( 'call_psolver_at_all_substeps' )
                READ ( 13 )  call_psolver_at_all_substeps
             CASE ( 'canyon_height' )
                READ ( 13 )  canyon_height
             CASE ( 'canyon_wall_left' )
                READ ( 13 )  canyon_wall_left
             CASE ( 'canyon_wall_south' )
                READ ( 13 )  canyon_wall_south
             CASE ( 'canyon_width_x' )
                READ ( 13 )  canyon_width_x
             CASE ( 'canyon_width_y' )
                READ ( 13 )  canyon_width_y
             CASE ( 'cfl_factor' )
                READ ( 13 )  cfl_factor
             CASE ( 'cloud_droplets' )
                READ ( 13 )  cloud_droplets
             CASE ( 'cloud_physics' )
                READ ( 13 )  cloud_physics
             CASE ( 'cloud_scheme' )
                READ ( 13 )  cloud_scheme
             CASE ( 'cloud_top_radiation' )
                READ ( 13 )  cloud_top_radiation
             CASE ( 'collective_wait' )
                READ ( 13 )  collective_wait
             CASE ( 'conserve_volume_flow' )
                READ ( 13 )  conserve_volume_flow
             CASE ( 'conserve_volume_flow_mode' )
                READ ( 13 )  conserve_volume_flow_mode
             CASE ( 'constant_flux_layer' )
                READ ( 13 )  constant_flux_layer
             CASE ( 'coupling_start_time' )
                READ ( 13 )  coupling_start_time
             CASE ( 'current_timestep_number' )
                READ ( 13 )  current_timestep_number
             CASE ( 'cycle_mg' )
                READ ( 13 )  cycle_mg
             CASE ( 'day_of_year_init' )
                READ ( 13 )  day_of_year_init
             CASE ( 'do2d_xy_time_count' )
                READ ( 13 )  do2d_xy_time_count
             CASE ( 'do2d_xz_time_count' )
                READ ( 13 )  do2d_xz_time_count
             CASE ( 'do2d_yz_time_count' )
                READ ( 13 )  do2d_yz_time_count
             CASE ( 'do3d_time_count' )
                READ ( 13 )  do3d_time_count
             CASE ( 'dp_external' )
                READ ( 13 )  dp_external
             CASE ( 'dp_level_b' )
                READ ( 13 )  dp_level_b
             CASE ( 'dp_smooth' )
                READ ( 13 )  dp_smooth
             CASE ( 'dpdxy' )
                READ ( 13 )  dpdxy
             CASE ( 'dt_3d' )
                READ ( 13 )  dt_3d
             CASE ( 'dx' )
                READ ( 13 )  dx
             CASE ( 'dy' )
                READ ( 13 )  dy
             CASE ( 'dz' )
                READ ( 13 )  dz
             CASE ( 'dz_max' )
                READ ( 13 )  dz_max
             CASE ( 'dz_stretch_factor' )
                READ ( 13 )  dz_stretch_factor
             CASE ( 'dz_stretch_factor_array' )
                READ ( 13 )  dz_stretch_factor_array
             CASE ( 'dz_stretch_level' )
                READ ( 13 )  dz_stretch_level
             CASE ( 'dz_stretch_level_end' )
                READ ( 13 )  dz_stretch_level_end
             CASE ( 'dz_stretch_level_start' )
                READ ( 13 )  dz_stretch_level_start
             CASE ( 'e_min' )
                READ ( 13 )  e_min
             CASE ( 'fft_method' )
                READ ( 13 )  fft_method
             CASE ( 'first_call_lpm' )
                READ ( 13 )  first_call_lpm
             CASE ( 'galilei_transformation' )
                READ ( 13 )  galilei_transformation
             CASE ( 'hom' )
                READ ( 13 )  hom
             CASE ( 'hom_sum' )
                READ ( 13 )  hom_sum
             CASE ( 'humidity' )
                READ ( 13 )  humidity
             CASE ( 'inflow_damping_factor' )
                IF ( .NOT. ALLOCATED( inflow_damping_factor ) )  THEN
                   ALLOCATE( inflow_damping_factor(0:nz+1) )
                ENDIF
                READ ( 13 )  inflow_damping_factor
             CASE ( 'inflow_damping_height' )
                READ ( 13 )  inflow_damping_height
             CASE ( 'inflow_damping_width' )
                READ ( 13 )  inflow_damping_width
             CASE ( 'inflow_disturbance_begin' )
                READ ( 13 )  inflow_disturbance_begin
             CASE ( 'inflow_disturbance_end' )
                READ ( 13 )  inflow_disturbance_end
             CASE ( 'km_constant' )
                READ ( 13 )  km_constant
             CASE ( 'large_scale_forcing' )
                READ ( 13 )  large_scale_forcing
             CASE ( 'large_scale_subsidence' )
                READ ( 13 )  large_scale_subsidence
             CASE ( 'latitude' )
                READ ( 13 )  latitude
             CASE ( 'longitude' )
                READ ( 13 )  longitude
             CASE ( 'loop_optimization' )
                READ ( 13 )  loop_optimization
             CASE ( 'masking_method' )
                READ ( 13 )  masking_method
             CASE ( 'mean_inflow_profiles' )
                IF ( .NOT. ALLOCATED( mean_inflow_profiles ) )  THEN
                   ALLOCATE( mean_inflow_profiles(0:nz+1,7) )
                ENDIF
                READ ( 13 )  mean_inflow_profiles
             CASE ( 'mg_cycles' )
                READ ( 13 )  mg_cycles
             CASE ( 'mg_switch_to_pe0_level' )
                READ ( 13 )  mg_switch_to_pe0_level
             CASE ( 'momentum_advec' )
                READ ( 13 )  momentum_advec
             CASE ( 'most_method' )
                READ ( 13 )  most_method
             CASE ( 'netcdf_precision' )
                READ ( 13 )  netcdf_precision
             CASE ( 'neutral' )
                READ ( 13 )  neutral
             CASE ( 'ngsrb' )
                READ ( 13 )  ngsrb
             CASE ( 'nsor' )
                READ ( 13 )  nsor
             CASE ( 'nsor_ini' )
                READ ( 13 )  nsor_ini
             CASE ( 'nudging' )
                READ ( 13 )  nudging
             CASE ( 'num_leg' )
                READ ( 13 )  num_leg
             CASE ( 'nx' )
                READ ( 13 )  nx
                nx_on_file = nx
             CASE ( 'ny' )
                READ ( 13 )  ny
                ny_on_file = ny
             CASE ( 'ocean' )
                READ ( 13 )  ocean
             CASE ( 'old_dt' )
                READ ( 13 )  old_dt
             CASE ( 'omega' )
                READ ( 13 )  omega
             CASE ( 'omega_sor' )
                READ ( 13 )  omega_sor
             CASE ( 'output_for_t0' )
                READ (13)    output_for_t0
             CASE ( 'passive_scalar' )
                READ ( 13 )  passive_scalar
             CASE ( 'prandtl_number' )
                READ ( 13 )  prandtl_number
             CASE ( 'precipitation' )
                READ ( 13 ) precipitation
             CASE ( 'psolver' )
                READ ( 13 )  psolver
             CASE ( 'pt_damping_factor' )
                READ ( 13 )  pt_damping_factor
             CASE ( 'pt_damping_width' )
                READ ( 13 )  pt_damping_width
             CASE ( 'pt_init' )
                READ ( 13 )  pt_init
             CASE ( 'pt_reference' )
                READ ( 13 )  pt_reference
             CASE ( 'pt_surface' )
                READ ( 13 )  pt_surface
             CASE ( 'pt_surface_initial_change' )
                READ ( 13 )  pt_surface_initial_change
             CASE ( 'pt_vertical_gradient' )
                READ ( 13 )  pt_vertical_gradient
             CASE ( 'pt_vertical_gradient_level' )
                READ ( 13 )  pt_vertical_gradient_level
             CASE ( 'pt_vertical_gradient_level_ind' )
                READ ( 13 )  pt_vertical_gradient_level_ind
             CASE ( 'q_init' )
                READ ( 13 )  q_init
             CASE ( 'q_surface' )
                READ ( 13 )  q_surface
             CASE ( 'q_surface_initial_change' )
                READ ( 13 )  q_surface_initial_change
             CASE ( 'q_vertical_gradient' )
                READ ( 13 )  q_vertical_gradient
             CASE ( 'q_vertical_gradient_level' )
                READ ( 13 )  q_vertical_gradient_level
             CASE ( 'q_vertical_gradient_level_ind' )
                READ ( 13 )  q_vertical_gradient_level_ind
             CASE ( 'random_generator' )
                READ ( 13 )  random_generator
             CASE ( 'random_heatflux' )
                READ ( 13 )  random_heatflux
             CASE ( 'rans_mode' )
                READ ( 13 )  rans_mode
             CASE ( 'rayleigh_damping_factor' )
                READ ( 13 )  rayleigh_damping_factor
             CASE ( 'rayleigh_damping_height' )
                READ ( 13 )  rayleigh_damping_height
             CASE ( 'recycling_width' )
                READ ( 13 )  recycling_width
             CASE ( 'recycling_yshift' )
                READ ( 13 ) recycling_yshift
             CASE ( 'ref_state' )
                READ ( 13 )  ref_state
             CASE ( 'reference_state' )
                READ ( 13 )  reference_state
             CASE ( 'residual_limit' )
                READ ( 13 )  residual_limit
             CASE ( 'roughness_length' )
                READ ( 13 )  roughness_length
             CASE ( 'run_coupled' )
                READ ( 13 )  run_coupled
             CASE ( 'runnr' )
                READ ( 13 )  runnr
             CASE ( 's_init' )
                READ ( 13 )  s_init
             CASE ( 's_surface' )
                READ ( 13 )  s_surface
             CASE ( 's_surface_initial_change' )
                READ ( 13 )  s_surface_initial_change
             CASE ( 's_vertical_gradient' )
                READ ( 13 )  s_vertical_gradient
             CASE ( 's_vertical_gradient_level' )
                READ ( 13 )  s_vertical_gradient_level
             CASE ( 's_vertical_gradient_level_ind' )
                READ ( 13 )  s_vertical_gradient_level_ind
             CASE ( 'sa_init' )
                READ ( 13 )  sa_init
             CASE ( 'sa_surface' )
                READ ( 13 )  sa_surface
             CASE ( 'sa_vertical_gradient' )
                READ ( 13 )  sa_vertical_gradient
             CASE ( 'sa_vertical_gradient_level' )
                READ ( 13 )  sa_vertical_gradient_level
             CASE ( 'scalar_advec' )
                READ ( 13 )  scalar_advec
             CASE ( 'simulated_time' )
                READ ( 13 )  simulated_time
            CASE ( 'spinup_time' )
                READ ( 13 )  spinup_time
             CASE ( 'surface_heatflux' )
                READ ( 13 )  surface_heatflux
             CASE ( 'surface_pressure' )
                READ ( 13 )  surface_pressure
             CASE ( 'surface_scalarflux' )
                READ ( 13 )  surface_scalarflux
             CASE ( 'surface_waterflux' )
                READ ( 13 )  surface_waterflux
             CASE ( 'time_coupling' )
                READ ( 13 )  time_coupling
             CASE ( 'time_disturb' )
                READ ( 13 )  time_disturb
             CASE ( 'time_do2d_xy' )
                READ ( 13 )  time_do2d_xy
             CASE ( 'time_do2d_xz' )
                READ ( 13 )  time_do2d_xz
             CASE ( 'time_do2d_yz' )
                READ ( 13 )  time_do2d_yz
             CASE ( 'time_do3d' )
                READ ( 13 )  time_do3d
             CASE ( 'time_do_av' )
                READ ( 13 )  time_do_av
             CASE ( 'time_do_sla' )
                READ ( 13 )  time_do_sla
             CASE ( 'time_domask' )
                READ ( 13 )  time_domask
             CASE ( 'time_dopr' )
                READ ( 13 )  time_dopr
             CASE ( 'time_dopr_av' )
                READ ( 13 )  time_dopr_av
             CASE ( 'time_dopr_listing' )
                READ ( 13 )  time_dopr_listing
             CASE ( 'time_dopts' )
                READ ( 13 )  time_dopts
             CASE ( 'time_dosp' )
                READ ( 13 )  time_dosp
             CASE ( 'time_dots' )
                READ ( 13 )  time_dots
             CASE ( 'time_restart' )
                READ ( 13 )  time_restart
             CASE ( 'time_run_control' )
                READ ( 13 )  time_run_control
             CASE ( 'time_since_reference_point' )
                READ ( 13 )  time_since_reference_point
             CASE ( 'time_utc_init' )
                READ ( 13 )  time_utc_init
             CASE ( 'timestep_scheme' )
                READ ( 13 )  timestep_scheme
             CASE ( 'top_heatflux' )
                READ ( 13 )  top_heatflux
             CASE ( 'top_momentumflux_u' )
                READ ( 13 )  top_momentumflux_u
             CASE ( 'top_momentumflux_v' )
                READ ( 13 )  top_momentumflux_v
             CASE ( 'top_salinityflux' )
                READ ( 13 )  top_salinityflux
             CASE ( 'top_scalarflux' )
                READ ( 13 )  top_scalarflux
             CASE ( 'topography' )
                READ ( 13 )  topography
             CASE ( 'topography_grid_convention' )
                READ ( 13 )  topography_grid_convention
             CASE ( 'tsc' )
                READ ( 13 )  tsc
             CASE ( 'tunnel_height' )
                READ ( 13 )  tunnel_height
             CASE ( 'tunnel_length' )
                READ ( 13 )  tunnel_length
             CASE ( 'tunnel_wall_depth' )
                READ ( 13 )  tunnel_wall_depth
             CASE ( 'tunnel_width_x' )
                READ ( 13 )  tunnel_width_x
             CASE ( 'tunnel_width_y' )
                READ ( 13 )  tunnel_width_y
             CASE ( 'turbulence_closure' )
                READ ( 13 )  turbulence_closure
             CASE ( 'turbulent_inflow' )
                READ ( 13 )  turbulent_inflow
             CASE ( 'u_bulk' )
                READ ( 13 )  u_bulk
             CASE ( 'u_init' )
                READ ( 13 )  u_init
             CASE ( 'u_max' )
                READ ( 13 )  u_max
             CASE ( 'u_max_ijk' )
                READ ( 13 )  u_max_ijk
             CASE ( 'ug' )
                READ ( 13 )  ug
             CASE ( 'ug_surface' )
                READ ( 13 )  ug_surface
             CASE ( 'ug_vertical_gradient' )
                READ ( 13 )  ug_vertical_gradient
             CASE ( 'ug_vertical_gradient_level' )
                READ ( 13 )  ug_vertical_gradient_level
             CASE ( 'ug_vertical_gradient_level_ind' )
                READ ( 13 )  ug_vertical_gradient_level_ind
             CASE ( 'use_surface_fluxes' )
                READ ( 13 )  use_surface_fluxes
             CASE ( 'use_top_fluxes' )
                READ ( 13 )  use_top_fluxes
             CASE ( 'use_ug_for_galilei_tr' )
                READ ( 13 )  use_ug_for_galilei_tr
             CASE ( 'use_upstream_for_tke' )
                READ ( 13 )  use_upstream_for_tke
             CASE ( 'v_bulk' )
                READ ( 13 )  v_bulk
             CASE ( 'v_init' )
                READ ( 13 )  v_init
             CASE ( 'v_max' )
                READ ( 13 )  v_max
             CASE ( 'v_max_ijk' )
                READ ( 13 )  v_max_ijk
             CASE ( 'vg' )
                READ ( 13 )  vg
             CASE ( 'vg_surface' )
                READ ( 13 )  vg_surface
             CASE ( 'vg_vertical_gradient' )
                READ ( 13 )  vg_vertical_gradient
             CASE ( 'vg_vertical_gradient_level' )
                READ ( 13 )  vg_vertical_gradient_level
             CASE ( 'vg_vertical_gradient_level_ind' )
                READ ( 13 )  vg_vertical_gradient_level_ind
             CASE ( 'volume_flow_area' )
                READ ( 13 )  volume_flow_area
             CASE ( 'volume_flow_initial' )
                READ ( 13 )  volume_flow_initial
             CASE ( 'subs_vertical_gradient' )
                READ ( 13 )  subs_vertical_gradient
             CASE ( 'subs_vertical_gradient_level' )
                READ ( 13 )  subs_vertical_gradient_level
             CASE ( 'subs_vertical_gradient_level_i' )
                READ ( 13 )  subs_vertical_gradient_level_i
             CASE ( 'w_max' )
                READ ( 13 )  w_max
             CASE ( 'w_max_ijk' )
                READ ( 13 )  w_max_ijk
             CASE ( 'wall_adjustment' )
                READ ( 13 )  wall_adjustment
             CASE ( 'wall_heatflux' )
                READ ( 13 )  wall_heatflux
             CASE ( 'wall_humidityflux' )
                READ ( 13 )  wall_humidityflux
             CASE ( 'wall_salinityflux' )
                READ ( 13 )  wall_salinityflux
             CASE ( 'wall_scalarflux' )
                READ ( 13 )  wall_scalarflux
             CASE ( 'y_shift' )
                READ ( 13 )  y_shift
             CASE ( 'z0h_factor' )
                READ ( 13 )  z0h_factor
             CASE ( 'zeta_max' )
                READ ( 13 )  zeta_max
             CASE ( 'zeta_min' )
                READ ( 13 )  zeta_min
             CASE ( 'z_i' )
                READ ( 13 )  z_i


             CASE DEFAULT

                IF ( .NOT. found )  THEN
                   WRITE( message_string, * ) 'unknown variable named "',      &
                                           restart_string(1:length),           &
                                          '" found in global data from ',      & 
                                          'prior run on PE ', myid
                CALL message( 'rrd_global', 'PA0302', 1, 2, 0, 6, 0 )
 
                ENDIF

          END SELECT
!
!--       Read next string
          READ ( 13 )  length
          READ ( 13 )  restart_string(1:length)    

       ENDDO
 

    CALL close_file( 13 )

    
    END SUBROUTINE rrd_global



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Skipping the global control variables from restart-file (binary format)
!> except some information needed when reading restart data from a previous
!> run which used a smaller total domain or/and a different domain decomposition
!> (initializing_actions  == 'cyclic_fill').
!------------------------------------------------------------------------------!

    SUBROUTINE rrd_read_parts_of_global


       USE arrays_3d,                                                             &
           ONLY:  inflow_damping_factor, mean_inflow_profiles, ref_state, ug, vg

       USE indices,                                                               &
           ONLY:  nz, nx, nx_on_file, ny, ny_on_file

       USE kinds

       USE pegrid

       USE statistics,                                                            &
           ONLY:  statistic_regions, hom, hom_sum, pr_palm, u_max, u_max_ijk,     &
                  v_max, v_max_ijk, w_max, w_max_ijk

       IMPLICIT NONE

       CHARACTER (LEN=10) ::  version_on_file
       CHARACTER (LEN=20) ::  momentum_advec_check
       CHARACTER (LEN=20) ::  scalar_advec_check
       CHARACTER (LEN=1)  ::  cdum

       INTEGER(iwp) ::  max_pr_user_on_file
       INTEGER(iwp) ::  nz_on_file
       INTEGER(iwp) ::  statistic_regions_on_file
       INTEGER(iwp) ::  tmp_mpru
       INTEGER(iwp) ::  tmp_sr

       REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  hom_sum_on_file
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  hom_on_file


       CALL check_open( 13 )

       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)
       READ ( 13 )  version_on_file

!
!-- Read number of PEs and horizontal index bounds of all PEs used in previous
!-- run
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'numprocs' )  THEN
          WRITE( message_string, * ) 'numprocs not found in data from prior ', &
                                     'run on PE ', myid
          CALL message( 'rrd_read_parts_of_global', 'PA0297', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  numprocs_previous_run

       IF ( .NOT. ALLOCATED( hor_index_bounds_previous_run ) )  THEN
          ALLOCATE( hor_index_bounds_previous_run(4,0:numprocs_previous_run-1) )
       ENDIF

       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'hor_index_bounds' )  THEN
          WRITE( message_string, * ) 'hor_index_bounds not found in data ',    &
                                     'from prior run on PE ', myid
          CALL message( 'rrd_read_parts_of_global', 'PA0298', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  hor_index_bounds_previous_run

!
!-- Read vertical number of gridpoints and number of different areas used
!-- for computing statistics. Allocate arrays depending on these values,
!-- which are needed for the following read instructions.
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'nz' )  THEN
          message_string = 'nz not found in restart data file'
          CALL message( 'rrd_read_parts_of_global', 'PA0303', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  nz_on_file
       IF ( nz_on_file /= nz )  THEN
          WRITE( message_string, * ) 'mismatch concerning number of ',         &
                                     'gridpoints along z:',                    &
                                     '&nz on file    = "', nz_on_file, '"',    &
                                     '&nz from run   = "', nz, '"'
          CALL message( 'rrd_read_parts_of_global', 'PA0304', 1, 2, 0, 6, 0 )
       ENDIF

       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'max_pr_user' )  THEN
          message_string = 'max_pr_user not found in restart data file'
          CALL message( 'rrd_read_parts_of_global', 'PA0305', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  max_pr_user_on_file
       IF ( max_pr_user_on_file /= max_pr_user )  THEN
          WRITE( message_string, * ) 'number of user profiles on res',         &
                                     'tart data file differs from the ',       &
                                     'current run:&max_pr_user on file    = "',&
                                     max_pr_user_on_file, '"',                 &
                                     '&max_pr_user from run   = "',            &
                                     max_pr_user, '"'
          CALL message( 'rrd_read_parts_of_global', 'PA0306', 0, 0, 0, 6, 0 )
          tmp_mpru = MIN( max_pr_user_on_file, max_pr_user )
       ELSE
          tmp_mpru = max_pr_user
       ENDIF

       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'statistic_regions' )  THEN
          message_string = 'statistic_regions not found in restart data file'
          CALL message( 'rrd_read_parts_of_global', 'PA0307', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  statistic_regions_on_file
       IF ( statistic_regions_on_file /= statistic_regions )  THEN
          WRITE( message_string, * ) 'statistic regions on restart data file ',&
                                     'differ from the current run:',           &
                                     '&statistic regions on file    = "',      &
                                     statistic_regions_on_file, '"',           &
                                     '&statistic regions from run   = "',      &
                                      statistic_regions, '"',                  &
                                     '&statistic data may be lost!'
          CALL message( 'rrd_read_parts_of_global', 'PA0308', 0, 1, 0, 6, 0 )
          tmp_sr = MIN( statistic_regions_on_file, statistic_regions )
       ELSE
          tmp_sr = statistic_regions
       ENDIF

!
!-- Now read and check some control parameters and skip the rest
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       DO  WHILE ( restart_string(1:length) /= 'binary_version_local' )

          SELECT CASE ( restart_string(1:length) )

             CASE ( 'average_count_pr' )
                READ ( 13 )  average_count_pr
                IF ( average_count_pr /= 0 )  THEN
                   WRITE( message_string, * ) 'inflow profiles not ',          &
                                  'temporally averaged. &Averaging will be ',  &
                                  'done now using', average_count_pr,          & 
                                  ' samples.'
                   CALL message( 'rrd_read_parts_of_global', 'PA0309',         &
                                 0, 1, 0, 6, 0 )
                ENDIF

             CASE ( 'hom' )
                ALLOCATE( hom_on_file(0:nz+1,2,pr_palm+max_pr_user_on_file,    &
                          0:statistic_regions_on_file) )
                READ ( 13 )  hom_on_file
                hom(:,:,1:pr_palm+tmp_mpru,0:tmp_sr) =                         &
                             hom_on_file(:,:,1:pr_palm+tmp_mpru,0:tmp_sr)
                DEALLOCATE( hom_on_file )

             CASE ( 'hom_sum' )
                ALLOCATE( hom_sum_on_file(0:nz+1,pr_palm+max_pr_user_on_file,  &
                          0:statistic_regions_on_file) )
                READ ( 13 )  hom_sum_on_file
                hom_sum(:,1:pr_palm+tmp_mpru,0:tmp_sr) =                       &
                             hom_sum_on_file(:,1:pr_palm+tmp_mpru,0:tmp_sr)
                DEALLOCATE( hom_sum_on_file )

             CASE ( 'momentum_advec' )
                momentum_advec_check = momentum_advec
                READ ( 13 )  momentum_advec
                IF ( TRIM( momentum_advec_check ) /= TRIM( momentum_advec ) )  &
                THEN
                   WRITE( message_string, * ) 'momentum_advec of the restart ',&
                                  'run differs from momentum_advec of the ',   &
                                  'initial run.'
                   CALL message( 'rrd_read_parts_of_global', 'PA0100',         &
                                 1, 2, 0, 6, 0 )
                ENDIF

             CASE ( 'nx' )
                READ ( 13 )  nx_on_file

             CASE ( 'ny' )
                READ ( 13 )  ny_on_file

             CASE ( 'ref_state' )
                READ ( 13 )  ref_state

             CASE ( 'scalar_advec' )
                scalar_advec_check = scalar_advec
                READ ( 13 )  scalar_advec
                IF ( TRIM( scalar_advec_check ) /= TRIM( scalar_advec ) )      &
                THEN
                   WRITE( message_string, * ) 'scalar_advec of the restart ',  &
                                  'run differs from scalar_advec of the ',     &
                                  'initial run.'
                   CALL message( 'rrd_read_parts_of_global', 'PA0101',         &
                                 1, 2, 0, 6, 0 )
                ENDIF

             CASE DEFAULT

                READ ( 13 )  cdum

          END SELECT

          READ ( 13 )  length
          READ ( 13 )  restart_string(1:length)

       ENDDO

!
!-- Calculate the temporal average of vertical profiles, if neccessary
    IF ( average_count_pr /= 0 )  THEN
       hom_sum = hom_sum / REAL( average_count_pr, KIND=wp )
    ENDIF


    CALL close_file( 13 )


    END SUBROUTINE rrd_read_parts_of_global


! Description:
! ------------
!> Reads processor specific data of variables and arrays from restart file 
!> (binary format).
!------------------------------------------------------------------------------!
 SUBROUTINE rrd_local 
 

    USE arrays_3d,                                                             &
        ONLY:  e, kh, km, p, pt, q, ql, qc, nc, nr, prr, precipitation_amount, &
               qr, s, sa, u, u_m_l, u_m_n, u_m_r, u_m_s, v, v_m_l, v_m_n,      &
               v_m_r, v_m_s, vpt, w, w_m_l, w_m_n, w_m_r, w_m_s

    USE averaging

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, nx_on_file, ny, nys, nysg, nyn, &
               nyng, ny_on_file, nzb, nzt

    USE kinds

    USE pegrid

    USE random_function_mod,                                                   &
        ONLY:  random_iv, random_iy

    USE random_generator_parallel,                                             &
        ONLY:  id_random_array, seq_random_array

    USE surface_mod,                                                           &
        ONLY :  surface_rrd_local

    IMPLICIT NONE

    CHARACTER (LEN=7)  ::  myid_char_save
    CHARACTER (LEN=10) ::  binary_version_local
    CHARACTER (LEN=10) ::  version_on_file

    INTEGER(iwp) ::  files_to_be_opened  !<
    INTEGER(iwp) ::  i                   !<
    INTEGER(iwp) ::  j                   !<
    INTEGER(iwp) ::  k                   !<
    INTEGER(iwp) ::  myid_on_file        !<
    INTEGER(iwp) ::  numprocs_on_file    !<
    INTEGER(iwp) ::  nxlc                !<
    INTEGER(iwp) ::  nxlf                !<
    INTEGER(iwp) ::  nxlpr               !<
    INTEGER(iwp) ::  nxl_on_file         !<
    INTEGER(iwp) ::  nxrc                !<
    INTEGER(iwp) ::  nxrf                !<
    INTEGER(iwp) ::  nxrpr               !<
    INTEGER(iwp) ::  nxr_on_file         !<
    INTEGER(iwp) ::  nync                !<
    INTEGER(iwp) ::  nynf                !<
    INTEGER(iwp) ::  nynpr               !<
    INTEGER(iwp) ::  nyn_on_file         !<
    INTEGER(iwp) ::  nysc                !<
    INTEGER(iwp) ::  nysf                !<
    INTEGER(iwp) ::  nyspr               !<
    INTEGER(iwp) ::  nys_on_file         !<
    INTEGER(iwp) ::  nzb_on_file         !<
    INTEGER(iwp) ::  nzt_on_file         !<
    INTEGER(iwp) ::  offset_x            !<
    INTEGER(iwp) ::  offset_y            !<
    INTEGER(iwp) ::  shift_x             !<
    INTEGER(iwp) ::  shift_y             !<

    INTEGER(iwp), DIMENSION(numprocs_previous_run) ::  file_list       !<
    INTEGER(iwp), DIMENSION(numprocs_previous_run) ::  overlap_count   !<

    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  nxlfa      !<
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  nxrfa      !<
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  nynfa      !<
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  nysfa      !<
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  offset_xa  !<
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  offset_ya  !<

    INTEGER(isp), DIMENSION(:,:),   ALLOCATABLE ::  tmp_2d_id_random   !< temporary array for storing random generator data
    INTEGER(isp), DIMENSION(:,:,:), ALLOCATABLE ::  tmp_2d_seq_random  !< temporary array for storing random generator data

    LOGICAL ::  found

    REAL(wp) ::  rdummy

    REAL(wp), DIMENSION(:,:),   ALLOCATABLE   ::  tmp_2d      !< temporary array for storing 2D data
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3d      !< temporary array for storing 3D data
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwul   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwun   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwur   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwus   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwvl   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwvn   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwvr   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwvs   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwwl   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwwn   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwwr   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3dwws   !<


!
!-- Read data from previous model run.
    CALL cpu_log( log_point_s(14), 'rrd_local', 'start' )

!
!-- Check which of the restart files contain data needed for the subdomain
!-- of this PE
    files_to_be_opened = 0

    DO  i = 1, numprocs_previous_run
!
!--    Store array bounds of the previous run ("pr") in temporary scalars
       nxlpr = hor_index_bounds_previous_run(1,i-1)
       nxrpr = hor_index_bounds_previous_run(2,i-1)
       nyspr = hor_index_bounds_previous_run(3,i-1)
       nynpr = hor_index_bounds_previous_run(4,i-1)

!
!--    Determine the offsets. They may be non-zero in case that the total domain
!--    on file is smaller than the current total domain.
       offset_x = ( nxl / ( nx_on_file + 1 ) ) * ( nx_on_file + 1 )
       offset_y = ( nys / ( ny_on_file + 1 ) ) * ( ny_on_file + 1 )

!
!--    Start with this offset and then check, if the subdomain on file
!--    matches another time(s) in the current subdomain by shifting it
!--    for nx_on_file+1, ny_on_file+1 respectively
    
       shift_y = 0
       j       = 0
       DO WHILE (  nyspr+shift_y <= nyn-offset_y )
          
          IF ( nynpr+shift_y >= nys-offset_y ) THEN 

             shift_x = 0
             DO WHILE ( nxlpr+shift_x <= nxr-offset_x )
                
                IF ( nxrpr+shift_x >= nxl-offset_x ) THEN
                   j = j +1
                   IF ( j > 1000 )  THEN
!
!--                   Array bound exceeded
                      message_string = 'data from subdomain of previous' //    &
                                       ' run mapped more than 1000 times'
                      CALL message( 'rrd_local', 'PA0284', 2, 2, -1,           &
                                       6, 1 )
                   ENDIF

                   IF ( j == 1 )  THEN
                      files_to_be_opened = files_to_be_opened + 1
                      file_list(files_to_be_opened) = i-1
                   ENDIF
                      
                   offset_xa(files_to_be_opened,j) = offset_x + shift_x
                   offset_ya(files_to_be_opened,j) = offset_y + shift_y
!
!--                Index bounds of overlapping data
                   nxlfa(files_to_be_opened,j) = MAX( nxl-offset_x-shift_x,    &
                                                      nxlpr )
                   nxrfa(files_to_be_opened,j) = MIN( nxr-offset_x-shift_x,    &
                                                      nxrpr )
                   nysfa(files_to_be_opened,j) = MAX( nys-offset_y-shift_y,    & 
                                                      nyspr )
                   nynfa(files_to_be_opened,j) = MIN( nyn-offset_y-shift_y,    & 
                                                      nynpr )

                ENDIF

                shift_x = shift_x + ( nx_on_file + 1 )
             ENDDO
       
          ENDIF
             
          shift_y = shift_y + ( ny_on_file + 1 )             
       ENDDO
          
       IF ( j > 0 )  overlap_count(files_to_be_opened) = j
          
    ENDDO
   
!
!-- Save the id-string of the current process, since myid_char may now be used
!-- to open files created by PEs with other id.
    myid_char_save = myid_char

    IF ( files_to_be_opened /= 1  .OR.  numprocs /= numprocs_previous_run )    &
    THEN
       WRITE( message_string, * ) 'number of PEs or virtual PE-grid changed ', &
                        'in restart run & PE', myid, ' will read from files ', &
                         file_list(1:files_to_be_opened)
       CALL message( 'rrd_local', 'PA0285', 0, 0, 0, 6, 0 )
    ENDIF

!
!-- Read data from all restart files determined above
    DO  i = 1, files_to_be_opened
 
       j = file_list(i)
!
!--    Set the filename (underscore followed by four digit processor id)
       WRITE (myid_char,'(''_'',I6.6)')  j

!
!--    Open the restart file. If this file has been created by PE0 (_000000),
!--    the global variables at the beginning of the file have to be skipped
!--    first.
       CALL check_open( 13 )
       IF ( j == 0 )  CALL rrd_skip_global

!
!--    First compare the version numbers
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)
       READ ( 13 )  version_on_file
!
!--    Read number of processors, processor-id, and array ranges.
!--    Compare the array ranges with those stored in the index bound array.
       READ ( 13 )  numprocs_on_file, myid_on_file, nxl_on_file, nxr_on_file,  &
                    nys_on_file, nyn_on_file, nzb_on_file, nzt_on_file

       IF ( nxl_on_file /= hor_index_bounds_previous_run(1,j) )  THEN
          WRITE( message_string, * ) 'problem with index bound nxl on ',       &
                            'restart file "', myid_char, '"',                  &
                            '&nxl = ', nxl_on_file, ' but it should be',       &
                            '&= ', hor_index_bounds_previous_run(1,j),         &
                            '&from the index bound information array'
          CALL message( 'rrd_local', 'PA0287', 2, 2, -1, 6, 1 )
       ENDIF

       IF ( nxr_on_file /= hor_index_bounds_previous_run(2,j) )  THEN
           WRITE( message_string, * ) 'problem with index bound nxr on ',      &
                               'restart file "', myid_char, '"'  ,             &
                               ' nxr = ', nxr_on_file, ' but it should be',    &
                               ' = ', hor_index_bounds_previous_run(2,j),      &
                               ' from the index bound information array' 
          CALL message( 'rrd_local', 'PA0288', 2, 2, -1, 6, 1 )

       ENDIF

       IF ( nys_on_file /= hor_index_bounds_previous_run(3,j) )  THEN
          WRITE( message_string, * ) 'problem with index bound nys on ',       &
                                 'restart file "', myid_char, '"',             &
                                 '&nys = ', nys_on_file, ' but it should be',  &
                                 '&= ', hor_index_bounds_previous_run(3,j),    &
                                 '&from the index bound information array'
          CALL message( 'rrd_local', 'PA0289', 2, 2, -1, 6, 1 ) 
       ENDIF

       IF ( nyn_on_file /= hor_index_bounds_previous_run(4,j) )  THEN
          WRITE( message_string, * ) 'problem with index bound nyn on ',       &
                               'restart file "', myid_char, '"',               &
                               '&nyn = ', nyn_on_file, ' but it should be',    &
                               '&= ', hor_index_bounds_previous_run(4,j),      &
                               '&from the index bound information array'
          CALL message( 'rrd_local', 'PA0290', 2, 2, -1, 6, 1 ) 
       ENDIF

       IF ( nzb_on_file /= nzb )  THEN
          WRITE( message_string, * ) 'mismatch between actual data and data ', &
                                     'from prior run on PE ', myid,            &
                                     '&nzb on file = ', nzb_on_file,           &
                                     '&nzb         = ', nzb
          CALL message( 'rrd_local', 'PA0291', 1, 2, 0, 6, 0 )  
       ENDIF

       IF ( nzt_on_file /= nzt )  THEN
          WRITE( message_string, * ) 'mismatch between actual data and data ', &
                                     'from prior run on PE ', myid,            &
                                     '&nzt on file = ', nzt_on_file,           &
                                     '&nzt         = ', nzt
          CALL message( 'rrd_local', 'PA0292', 1, 2, 0, 6, 0 )  
       ENDIF

!
!--    Allocate temporary arrays sized as the arrays on the restart file
       ALLOCATE( tmp_2d(nys_on_file-nbgp:nyn_on_file+nbgp,                     &
                        nxl_on_file-nbgp:nxr_on_file+nbgp),                    &
                 tmp_3d(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,           &
                        nxl_on_file-nbgp:nxr_on_file+nbgp) )

!
!--    Read arrays
!--    ATTENTION: If the following read commands have been altered, the 
!--    ---------- version number of the variable binary_version_local must 
!--               be altered, too. Furthermore, the output list of arrays in 
!--               wrd_write_local must also be altered 
!--               accordingly.
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)
       

!
!--    Loop over processor specific field data
       DO  WHILE ( restart_string(1:length) /= '*** end ***' )

!
!--       Map data on file as often as needed (data are read only for k=1)
          DO  k = 1, overlap_count(i)

             found = .FALSE.
             
!
!--          Get the index range of the subdomain on file which overlap with
!--          the current subdomain
             nxlf = nxlfa(i,k)
             nxlc = nxlfa(i,k) + offset_xa(i,k)
             nxrf = nxrfa(i,k)
             nxrc = nxrfa(i,k) + offset_xa(i,k)
             nysf = nysfa(i,k)
             nysc = nysfa(i,k) + offset_ya(i,k)
             nynf = nynfa(i,k)
             nync = nynfa(i,k) + offset_ya(i,k)


             SELECT CASE ( restart_string(1:length) )

                CASE ( 'ghf_av' )
                   IF ( .NOT. ALLOCATED( ghf_av ) )  THEN
                      ALLOCATE( ghf_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   ghf_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                      tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'e' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   e(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =              &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'e_av' )
                   IF ( .NOT. ALLOCATED( e_av ) )  THEN
                      ALLOCATE( e_av(nzb:nzt+1,nys-nbgp:nyn+nbgp,              &
                                     nxl-nbgp:nxr+nbgp) )    
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   e_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'iran' ) ! matching random numbers is still unresolved
                                ! issue
                   IF ( k == 1 )  READ ( 13 )  iran

                CASE ( 'kh' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   kh(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =             &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'kh_av' )
                   IF ( .NOT. ALLOCATED( kh_av ) )  THEN
                      ALLOCATE( kh_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg ))
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   kh_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'km' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   km(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =             &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'km_av' )
                   IF ( .NOT. ALLOCATED( km_av ) )  THEN
                      ALLOCATE( km_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg ))
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   km_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'lpt_av' )
                   IF ( .NOT. ALLOCATED( lpt_av ) )  THEN
                      ALLOCATE( lpt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg ))
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   lpt_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =         &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'lwp_av' )
                   IF ( .NOT. ALLOCATED( lwp_av ) )  THEN
                      ALLOCATE( lwp_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   lwp_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                      tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'nc' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   nc(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =             &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'nc_av' )
                   IF ( .NOT. ALLOCATED( nc_av ) )  THEN
                      ALLOCATE( nc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   nc_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)


                CASE ( 'nr' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   nr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =             &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'nr_av' )
                   IF ( .NOT. ALLOCATED( nr_av ) )  THEN
                      ALLOCATE( nr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   nr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'p' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   p(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =              &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'p_av' )
                   IF ( .NOT. ALLOCATED( p_av ) )  THEN
                      ALLOCATE( p_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   p_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'pc_av' )
                   IF ( .NOT. ALLOCATED( pc_av ) )  THEN
                      ALLOCATE( pc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   pc_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'pr_av' )
                   IF ( .NOT. ALLOCATED( pr_av ) )  THEN
                      ALLOCATE( pr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   pr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'prr' )
                   IF ( .NOT. ALLOCATED( prr ) )  THEN
                      ALLOCATE( prr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   prr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =            &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'prr_av' )
                   IF ( .NOT. ALLOCATED( prr_av ) )  THEN
                      ALLOCATE( prr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   prr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =         &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'precipitation_amount' )
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   precipitation_amount(nysc-nbgp:nync+nbgp,                   &
                                        nxlc-nbgp:nxrc+nbgp)  =                &
                         tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'pt' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   pt(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =             &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'pt_av' )
                   IF ( .NOT. ALLOCATED( pt_av ) )  THEN
                      ALLOCATE( pt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   pt_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'q' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   q(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =              &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'q_av' )
                   IF ( .NOT. ALLOCATED( q_av ) )  THEN
                      ALLOCATE( q_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg ))
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   q_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qc' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   qc(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =             &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qc_av' )
                   IF ( .NOT. ALLOCATED( qc_av ) )  THEN
                      ALLOCATE( qc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   qc_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'ql' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   ql(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =             &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'ql_av' )
                   IF ( .NOT. ALLOCATED( ql_av ) )  THEN
                      ALLOCATE( ql_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   ql_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qr' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   qr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =             &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qr_av' )
                   IF ( .NOT. ALLOCATED( qr_av ) )  THEN
                      ALLOCATE( qr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   qr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'ql_c_av' )
                   IF ( .NOT. ALLOCATED( ql_c_av ) )  THEN
                      ALLOCATE( ql_c_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   ql_c_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =        &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'ql_v_av' )
                   IF ( .NOT. ALLOCATED( ql_v_av ) )  THEN
                      ALLOCATE( ql_v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   ql_v_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =        &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'ql_vp_av' )
                   IF ( .NOT. ALLOCATED( ql_vp_av ) )  THEN
                      ALLOCATE( ql_vp_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   ql_vp_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =       &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qsws_av' )
                   IF ( .NOT. ALLOCATED( qsws_av ) )  THEN
                      ALLOCATE( qsws_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF  
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   qsws_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =         &
                      tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'qv_av' )
                   IF ( .NOT. ALLOCATED( qv_av ) )  THEN
                      ALLOCATE( qv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   qv_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'r_a_av' )
                   IF ( .NOT. ALLOCATED( r_a_av ) )  THEN
                      ALLOCATE( r_a_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   r_a_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                      tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'random_iv' )  ! still unresolved issue
                   IF ( k == 1 )  READ ( 13 )  random_iv
                   IF ( k == 1 )  READ ( 13 )  random_iy

                CASE ( 'seq_random_array' )
                   ALLOCATE( tmp_2d_id_random(nys_on_file:nyn_on_file,         &
                                              nxl_on_file:nxr_on_file) )
                   ALLOCATE( tmp_2d_seq_random(5,nys_on_file:nyn_on_file,      &
                                                 nxl_on_file:nxr_on_file) )
                   IF ( .NOT. ALLOCATED( id_random_array ) )  THEN
                      ALLOCATE( id_random_array(nys:nyn,nxl:nxr) )
                   ENDIF
                   IF ( .NOT. ALLOCATED( seq_random_array ) )  THEN
                      ALLOCATE( seq_random_array(5,nys:nyn,nxl:nxr) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d_id_random
                   IF ( k == 1 )  READ ( 13 )  tmp_2d_seq_random
                   id_random_array(nysc:nync,nxlc:nxrc) =                      &
                      tmp_2d_id_random(nysf:nynf,nxlf:nxrf)
                   seq_random_array(:,nysc:nync,nxlc:nxrc) =                   &
                      tmp_2d_seq_random(:,nysf:nynf,nxlf:nxrf)
                   DEALLOCATE( tmp_2d_id_random, tmp_2d_seq_random )

                CASE ( 'solar3d_av' )
                   IF ( .NOT. ALLOCATED( solar3d_av ) )  THEN
                      ALLOCATE( solar3d_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   solar3d_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =   &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)


                CASE ( 'rho_ocean_av' )
                   IF ( .NOT. ALLOCATED( rho_ocean_av ) )  THEN
                      ALLOCATE( rho_ocean_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   rho_ocean_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =   &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'alpha_T_av' )
                   IF ( .NOT. ALLOCATED( alpha_T_av ) )  THEN
                      ALLOCATE( alpha_T_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   alpha_T_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =   &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'beta_S_av' )
                   IF ( .NOT. ALLOCATED( beta_S_av ) )  THEN
                      ALLOCATE( beta_S_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   beta_S_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =   &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 's' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   s(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =              &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 's_av' )
                   IF ( .NOT. ALLOCATED( s_av ) )  THEN
                      ALLOCATE( s_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg))
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   s_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'sa' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   sa(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =             &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'sa_av' )
                   IF ( .NOT. ALLOCATED( sa_av ) )  THEN
                      ALLOCATE( sa_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   sa_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'shf_av' )
                   IF ( .NOT. ALLOCATED( shf_av ) )  THEN
                      ALLOCATE( shf_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   shf_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =          &
                      tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'shf_sol_av' )
                   IF ( .NOT. ALLOCATED( shf_sol_av ) )  THEN
                      ALLOCATE( shf_sol_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   shf_sol_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =          &
                      tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'ssws_av' )
                   IF ( .NOT. ALLOCATED( ssws_av ) )  THEN
                      ALLOCATE( ssws_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF  
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   ssws_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =         &
                      tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
                 
                CASE ( 'ts_av' )
                   IF ( .NOT. ALLOCATED( ts_av ) )  THEN
                      ALLOCATE( ts_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   ts_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =           &
                        tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'tsurf_av' )
                   IF ( .NOT. ALLOCATED( tsurf_av ) )  THEN
                      ALLOCATE( tsurf_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   tsurf_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =         &
                      tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'u' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   u(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =              &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'u_av' )
                   IF ( .NOT. ALLOCATED( u_av ) )  THEN
                      ALLOCATE( u_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   u_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'u_m_l' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwul(nzb:nzt+1,                           &
                                nys_on_file-nbgp:nyn_on_file+nbgp,1:2) )
                      READ ( 13 )  tmp_3dwul
                   ENDIF
                   IF ( outflow_l )  THEN
                      u_m_l(:,nysc-nbgp:nync+nbgp,:) =                         & 
                         tmp_3dwul(:,nysf-nbgp:nynf+nbgp,:)
                   ENDIF

                CASE ( 'u_m_n' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwun(nzb:nzt+1,ny-1:ny,                   &
                                          nxl_on_file-nbgp:nxr_on_file+nbgp) )
                      READ ( 13 )  tmp_3dwun
                   ENDIF
                   IF ( outflow_n )  THEN
                      u_m_n(:,:,nxlc-nbgp:nxrc+nbgp) =                         & 
                         tmp_3dwun(:,:,nxlf-nbgp:nxrf+nbgp)
                   ENDIF

                CASE ( 'u_m_r' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwur(nzb:nzt+1,                           &
                                nys_on_file-nbgp:nyn_on_file+nbgp,nx-1:nx) )
                      READ ( 13 )  tmp_3dwur
                   ENDIF
                   IF ( outflow_r )  THEN
                      u_m_r(:,nysc-nbgp:nync+nbgp,:) =                         & 
                         tmp_3dwur(:,nysf-nbgp:nynf+nbgp,:)
                   ENDIF

                CASE ( 'u_m_s' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwus(nzb:nzt+1,0:1,                       &
                                          nxl_on_file-nbgp:nxr_on_file+nbgp) )
                      READ ( 13 )  tmp_3dwus
                   ENDIF
                   IF ( outflow_s )  THEN
                      u_m_s(:,:,nxlc-nbgp:nxrc+nbgp) =                         &
                         tmp_3dwus(:,:,nxlf-nbgp:nxrf+nbgp)
                   ENDIF

                CASE ( 'us_av' )
                   IF ( .NOT. ALLOCATED( us_av ) )  THEN
                      ALLOCATE( us_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   us_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =           &
                      tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'v' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   v(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =              &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'v_av' )
                   IF ( .NOT. ALLOCATED( v_av ) )  THEN
                      ALLOCATE( v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   v_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'v_m_l' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwvl(nzb:nzt+1,                           &
                                nys_on_file-nbgp:nyn_on_file+nbgp,0:1) )
                      READ ( 13 )  tmp_3dwvl
                   ENDIF
                   IF ( outflow_l )  THEN
                      v_m_l(:,nysc-nbgp:nync+nbgp,:) =                         & 
                         tmp_3dwvl(:,nysf-nbgp:nynf+nbgp,:)
                   ENDIF

                CASE ( 'v_m_n' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwvn(nzb:nzt+1,ny-1:ny,                   &
                                          nxl_on_file-nbgp:nxr_on_file+nbgp) )
                      READ ( 13 )  tmp_3dwvn
                   ENDIF
                   IF ( outflow_n )  THEN
                      v_m_n(:,:,nxlc-nbgp:nxrc+nbgp) =                         &
                         tmp_3dwvn(:,:,nxlf-nbgp:nxrf+nbgp)
                   ENDIF

                CASE ( 'v_m_r' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwvr(nzb:nzt+1,                           &
                                nys_on_file-nbgp:nyn_on_file+nbgp,nx-1:nx) )
                      READ ( 13 )  tmp_3dwvr
                   ENDIF
                   IF ( outflow_r )  THEN
                      v_m_r(:,nysc-nbgp:nync+nbgp,:) =                         &
                         tmp_3dwvr(:,nysf-nbgp:nynf+nbgp,:)
                   ENDIF

                CASE ( 'v_m_s' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwvs(nzb:nzt+1,1:2,                       &
                                          nxl_on_file-nbgp:nxr_on_file+nbgp) )
                      READ ( 13 )  tmp_3dwvs
                   ENDIF
                   IF ( outflow_s )  THEN
                      v_m_s(:,:,nxlc-nbgp:nxrc+nbgp) =                         &
                         tmp_3dwvs(:,:,nxlf-nbgp:nxrf+nbgp)
                   ENDIF

                CASE ( 'vpt' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   vpt(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =            &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'vpt_av' )
                   IF ( .NOT. ALLOCATED( vpt_av ) )  THEN
                      ALLOCATE( vpt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   vpt_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =         &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'w' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   w(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =              &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'w_av' )
                   IF ( .NOT. ALLOCATED( w_av ) )  THEN
                      ALLOCATE( w_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   w_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'w_m_l' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwwl(nzb:nzt+1,&
                                nys_on_file-nbgp:nyn_on_file+nbgp,0:1) )
                      READ ( 13 )  tmp_3dwwl
                   ENDIF
                   IF ( outflow_l )  THEN
                      w_m_l(:,nysc-nbgp:nync+nbgp,:) =                         & 
                         tmp_3dwwl(:,nysf-nbgp:nynf+nbgp,:)
                   ENDIF

                CASE ( 'w_m_n' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwwn(nzb:nzt+1,ny-1:ny, &
                                          nxl_on_file-nbgp:nxr_on_file+nbgp) )
                      READ ( 13 )  tmp_3dwwn
                   ENDIF
                   IF ( outflow_n )  THEN
                      w_m_n(:,:,nxlc-nbgp:nxrc+nbgp) =                         &
                         tmp_3dwwn(:,:,nxlf-nbgp:nxrf+nbgp)
                   ENDIF

                CASE ( 'w_m_r' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwwr(nzb:nzt+1,&
                                nys_on_file-nbgp:nyn_on_file+nbgp,nx-1:nx) )
                      READ ( 13 )  tmp_3dwwr
                   ENDIF
                   IF ( outflow_r )  THEN
                      w_m_r(:,nysc-nbgp:nync+nbgp,:) =                         & 
                         tmp_3dwwr(:,nysf-nbgp:nynf+nbgp,:)
                   ENDIF

                CASE ( 'w_m_s' )
                   IF ( k == 1 )  THEN
                      ALLOCATE( tmp_3dwws(nzb:nzt+1,0:1, &
                                          nxl_on_file-nbgp:nxr_on_file+nbgp) )
                      READ ( 13 )  tmp_3dwws
                   ENDIF
                   IF ( outflow_s )  THEN
                      w_m_s(:,:,nxlc-nbgp:nxrc+nbgp) =                         & 
                         tmp_3dwws(:,:,nxlf-nbgp:nxrf+nbgp)
                   ENDIF

                CASE ( 'z0_av' )
                   IF ( .NOT. ALLOCATED( z0_av ) )  THEN
                      ALLOCATE( z0_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   z0_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =           &
                      tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'z0h_av' )
                   IF ( .NOT. ALLOCATED( z0h_av ) )  THEN
                      ALLOCATE( z0h_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   z0h_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =          &
                       tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'z0q_av' )
                   IF ( .NOT. ALLOCATED( z0q_av ) )  THEN
                      ALLOCATE( z0q_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   z0q_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =          &
                       tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE DEFAULT

!
!--                Read surface related variables
                   IF ( .NOT. found ) CALL surface_rrd_local( i, k, nxlf,      &
                                           nxlc, nxl_on_file, nxrf, nxrc,      & 
                                           nxr_on_file, nynf, nync,            & 
                                           nyn_on_file, nysf, nysc,            &
                                           nys_on_file, found )

                   IF ( .NOT. found )  THEN
                      WRITE( message_string, * ) 'unknown variable named "',   &
                                                 restart_string(1:length),     &
                                                '" found in subdomain data ',  & 
                                                'from prior run on PE ', myid
                      CALL message( 'rrd_local', 'PA0302', 1, 2, 0, 6, 0 )
 
                   ENDIF

             END SELECT

          ENDDO ! overlaploop 

!
!--       Deallocate arrays needed for specific variables only
          IF ( ALLOCATED( tmp_3dwul ) )  DEALLOCATE( tmp_3dwul )
          IF ( ALLOCATED( tmp_3dwun ) )  DEALLOCATE( tmp_3dwun )
          IF ( ALLOCATED( tmp_3dwur ) )  DEALLOCATE( tmp_3dwur )
          IF ( ALLOCATED( tmp_3dwus ) )  DEALLOCATE( tmp_3dwus )
          IF ( ALLOCATED( tmp_3dwvl ) )  DEALLOCATE( tmp_3dwvl )
          IF ( ALLOCATED( tmp_3dwvn ) )  DEALLOCATE( tmp_3dwvn )
          IF ( ALLOCATED( tmp_3dwvr ) )  DEALLOCATE( tmp_3dwvr )
          IF ( ALLOCATED( tmp_3dwvs ) )  DEALLOCATE( tmp_3dwvs )
          IF ( ALLOCATED( tmp_3dwwl ) )  DEALLOCATE( tmp_3dwwl )
          IF ( ALLOCATED( tmp_3dwwn ) )  DEALLOCATE( tmp_3dwwn )
          IF ( ALLOCATED( tmp_3dwwr ) )  DEALLOCATE( tmp_3dwwr )
          IF ( ALLOCATED( tmp_3dwws ) )  DEALLOCATE( tmp_3dwws )

!
!--       Read next character string
          READ ( 13 )  length
          READ ( 13 )  restart_string(1:length)

       ENDDO ! dataloop

!
!--    Close the restart file
       CALL close_file( 13 )

       DEALLOCATE( tmp_2d, tmp_3d )

    ENDDO  ! loop over restart files


!
!-- Restore the original filename for the restart file to be written
    myid_char = myid_char_save

!
!-- End of time measuring for reading binary data
    CALL cpu_log( log_point_s(14), 'rrd_local', 'stop' )

 END SUBROUTINE rrd_local


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Skipping the global control variables from restart-file (binary format)
!------------------------------------------------------------------------------!

    SUBROUTINE rrd_skip_global


       IMPLICIT NONE

       CHARACTER (LEN=10) ::  version_on_file

       CHARACTER (LEN=1) ::  cdum


       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       DO  WHILE ( restart_string(1:length) /= 'binary_version_local' )

          READ ( 13 )  cdum
          READ ( 13 )  length 
          READ ( 13 )  restart_string(1:length)

       ENDDO

       BACKSPACE ( 13 )
       BACKSPACE ( 13 )


    END SUBROUTINE rrd_skip_global


 END MODULE read_restart_data_mod
