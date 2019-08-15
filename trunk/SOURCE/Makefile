#--------------------------------------------------------------------------------#
# This file is part of the PALM model system.
#
# PALM is free software: you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# PALM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PALM. If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 1997-2018  Leibniz Universitaet Hannover
#--------------------------------------------------------------------------------#
#
# WARNING: don't write filenames with extension .f90 in this header!!!!
#
# Current revisions:
# ------------------
#
#
# Former revisions:
# -----------------
# $Id: Makefile 3083 2018-06-19 14:03:12Z gronemeier $
# add turbulence_closure_mod to parin
#
# 2963 2018-04-12 14:47:44Z suehring
# Introduce index for vegetation/wall, pavement/green-wall and water/window
# surfaces, for clearer access of surface fraction, albedo, emissivity, etc. .
#
# 2955 2018-04-09 15:14:01Z suehring
# Add log-points to measure CPU time of NetCDF data input.
#
# 2938 2018-03-27 15:52:42Z suehring
# No initialization of child domains via dynamic input file, except for soil
# moisture and temperature
# Apply turbulence generator at non-cyclic lateral boundary in nesting case
#
# 2936 2018-03-27 14:49:27Z suehring
# Added dependencies for parent and child synchronization
#
# 2921 2018-03-22 15:05:23Z Giersch
# date_and_time_mod dependency has been added to read/write_restart_data_mod
#
# 2918 2018-03-21 15:52:14Z gronemeier
# read/write_3d_binary and read/write_var_list has been removed,
# read/write_restart_data_mod, wrd_write_string and
# user_read/write_restart_data_mod has been added, dependencies with respect to
# the aforementioned routines have been added/removed
#
# 2847 2018-03-02 21:45:58Z suehring
# Changed format and enforced sorting
#
# 2817 2018-02-19 16:32:21Z knoop
# Preliminary gust module interface implemented
#
# 2802 2018-02-14 16:21:39Z thiele
# Changed lpm from subroutine to module.
# Introduce particle transfer in nested models.
#
# 2773 2018-01-30 14:12:54Z suehring
# Nesting of chemical species
#
# 2718 2018-01-02 08:49:38Z maronga
# Corrected "Former revisions" section
#
# 2697 2017-12-14 17:57:20Z kanani
# Bugfix, missing dependencies
#
# 2696 2017-12-14 17:12:51Z kanani
# Change in file header (GPL part)
# Implementation of uv exposure model (FK)
# Bugfix, removed loop dependcy for vertical_nesting_mod and
# turbulence_closure_mod, added depencies for vertical_nesting_mod (TG)
# implemented turbulence_closure_mod (TG)
# Additional dependencies for chemistry module (FK)
# Additional dependencies for chemistry module (FK)
# Implementation of chemistry module (FK)
# For LSM, add dependency on calc_mean_profile (??)
# poismg_noopt modularized and renamed into poismg_noopt_mod
# add dependencies for netcdf_data_input_mod, calc_mean_profile,
# radiation_module_mod, land_surface_model_mod (MS)
#
# 2608 2017-11-13 14:04:26Z schwenkel
# Added diagnostic_quantities_mod
#
# 2600 2017-11-01 14:11:20Z raasch
# comment line concerning bound checks removed
#
# 2599 2017-11-01 13:18:45Z hellstea
# virtual_flight_mod, synthetic_turbulence_generator_mod and
# wind_turbine_model_mod were added to read_var_list
#
# 2563 2017-10-19 15:36:10Z Giersch
# wind_turbine_model_mod and synthetic_turbulence_generator_mod were added to
# write_var_list and virtual_flight_mod was deleted from read_var_list
#
# 2544 2017-10-13 18:09:32Z maronga
# Added date_and_time_mod
#
# 2371 2017-08-24 13:01:17Z kanani
# Corrected dependencies for vertical_nesting_mod
#
# 2370 2017-08-23 06:11:43Z raasch
# dependency bugfix for synthetic_turbulence_generator
#
# 2365 2017-08-21 14:59:59Z kanani
# Added dependencies for vertical_nesting_mod
#
# 2339 2017-08-07 13:55:26Z gronemeier
# corrected timestamp in header
#
# 2338 2017-08-07 12:15:38Z gronemeier
# Modularize 1D model
#
# 2320 2017-07-21 12:47:43Z suehring
# -ls_forcing nudging
# +large_scale_forcing_nudging
#
# 2318 2017-07-20 17:27:44Z suehring
# Add further dependencies on surface_mod
#
# 2317 2017-07-20 17:27:19Z suehring
# Added time_integration_spinup
#
# 2269 2017-06-09 11:57:32Z suehring
# Add dependency in read_3d_binary
#
# 2263 2017-06-08 14:59:01Z schwenkel
# Implemented splitting and merging algorithm
#
# 2259 2017-06-08 09:09:11Z gronemeier
# Implemented synthetic turbulence generator
#
# 2256 2017-06-07 13:58:08Z suehring
# Remove ring dependency in init_pegrid
#
# 2238 2017-05-31 16:49:16Z suehring
# Bugfix, further missing dependency on surface_mod
#
# 2237 2017-05-31 10:34:53Z suehring
# Bugfix, add dependencies on surface_mod for surface_coupler,
# plant_canopy_model_mod and ls_forcing_mod
#
# 2233 2017-05-30 18:08:54Z suehring
#
# 2232 2017-05-30 17:47:52Z suehring
# +dependencies for surface_mod
# -wall_fluxes
#
# 2130 2017-01-24 16:25:39Z raasch
# dependency for timestep updated
#
# 2118 2017-01-17 16:38:49Z raasch
# -cuda_fft_interfaces_mod
#
# 2050 2016-11-08 15:00:55Z gronemeier
# Implement turbulent outflow method
#
# 2007 2016-08-24 15:47:17Z kanani
# urban surface module added,
# cleaned up some lines (compiler flags/options), which were accidentally
# added in rev1938
#
# 1998 2016-08-20 18:45:34Z knoop
# Bugfix: added netcdf_interface to dependency list for user_init_land_surface
#
# 1986 2016-08-10 14:07:17Z gronemeier
# POSIX-calls module added
#
# 1972 2016-07-26 07:52:02Z maronga
# Removed some dependencies due to further modularization of land surface model
#
# 1957 2016-07-07 10:43:48Z suehring
# flight module added
#
# 1946 2016-06-15 11:38:13Z suehring
# Bugfix dependency init_grid on advec_ws
#
# 1938 2016-06-13 15:26:05Z hellstea
# Some dependency errors corrected
#
# 1934 2016-06-13 09:46:57Z hellstea
# poismg renamed poismg_noopt, poismg_fast_mod renamed poismg_mod
#
# 1914 2016-05-26 14:44:07Z witha
# Added wind_turbine_model_mod
#
# 1900 2016-05-04 15:27:53Z raasch
# pmc_general renamed pmc_general_mod
#
# 1871 2016-04-15 11:46:09Z hoffmann
# dependencies for package_parin corrected
#
# 1850 2016-04-08 13:29:27Z maronga
# Adapted for modularization of microphysics
# Several files renamed --> _mod
# Bugfix for previous commit
#
# 1849 2016-04-08 11:33:18Z hoffmann
# Adapted for modularization of microphysics
#
# 1837 2016-04-07 18:21:57Z raasch
# further bug removed from r1835
#
# 1835 2016-04-07 16:31:50Z raasch
# bugfix: merge comments removed
#
# 1833 2016-04-07 14:23:03Z raasch
# spectrum renamed spectra_mod, depencies for spectra changed
#
# 1826 2016-04-07 12:01:39Z maronga
# Renamed radiation_model to radiation_model_mod.
# Renamed plant_canopy_model to plant_canopy_model_mod.
#
# 1822 2016-04-07 07:49:42Z hoffmann
# Tails removed. lpm_release_set removed. calc_precipitation, impact_of_latent_heat
# removed.
#
# 1817 2016-04-06 15:44:20Z maronga
# Renamed land_surface_model to land_surface_model_mod.
#
# 1808 2016-04-05 19:44:00Z raasch
# -local_flush, -local_getenv
#
# 1797 2016-03-21 16:50:28Z raasch
# dependcy for check_for_restart updated
#
# 1791 2016-03-11 10:41:25Z raasch
# dependencies of header changed
#
# 1786 2016-03-08 05:49:27Z raasch
# rename calc_spectra to spectrum + modularization of spectrum
#
# 2016-03-06 18:36:17Z raasch
# name change of netcdf routines and respective dependencies
#
# 1779 2016-03-03 08:01:28Z raasch
# dependencies changed for init_peprid
#
# 1766 2016-02-29 08:37:15Z raasch
# update dependency for swap_timelevel
#
# 1764 2016-02-28 12:45:19Z raasch
# several updates of dependencies
#
# 1762 2016-02-25 12:31:13Z hellstea
# +pmc_interface, +pmc routines
#
# 1747 2016-02-08 12:25:53Z raasch
# dependencies for routine palm updated
#
# 1691 2015-10-26 16:17:44Z maronga
# Replaced prandtl_fluxes with surface_layer_fluxes. Added radiation_model to
# prognostic_equations
#
# 1585 2015-04-30 07:05:52Z maronga
# Added user_init_radiation.f90
#
# 1575 2015-03-27 09:56:27Z raasch
# +poismg_fast
#
# 1551 2015-03-03 14:18:16Z maronga
# Bugfix: further adjustments for the land surface model and radiation model
#
# 1517 2015-01-07 19:12:25Z hoffmann
# advec_s_bc added to prognostic_equations
#
# 1500 2014-12-03 17:42:41Z maronga
# Bugfix: missing adjustments for land surface model and radiation model
#
# 1496 2014-12-02 17:25:50Z maronga
# Added land surface model and radiation model files: land_surface_model,
# radiation_model, user_init_land_surface
#
# 1484 2014-10-21 10:53:05Z kanani
# plant_canopy_model-dependency added for check_parameters, header, init_3d_model,
# package_parin, read_var_list, user_init_plant_canopy, write_var_list
#
# 1444 2014-08-02 20:10:32Z letzel
# bugfix: cpulog added to lpm_advec
#
# 1404 2014-05-14 09:01:39Z keck
# bugfix: dependencies added for progress_bar
#
# 1402 2014-05-09 14:25:13Z raasch
# progress_bar added
#
# 1400 2014-05-09 14:03:54Z knoop
# Added new module random_generator_parallel
#
# 1380 2014-04-28 12:40:45Z heinze
# bugfix: mod_particle_attributes  added to check_open
# nudging added to time_integration
#
# 1374 2014-04-25 12:55:07Z raasch
# bugfix: missing dependency added for check_open
#
# 1365 2014-04-22 15:03:56Z boeske
# Added new module calc_mean_profile, previously in module buoyancy,
# removed buoyancy dependency from nudging
#
# 1363 2014-04-17 12:28:49Z keck
# bugfix: cpulog added to lpm_pack_arrays
#
# 1361 2014-04-16 15:17:48Z hoffmann
# cpulog added to microphysics
#
# 1359 2014-04-11 17:15:14Z hoffmann
# mod_particle_attributes added, lpm_sort_arrays removed,
# lpm_extend_particle_array removed
#
# 1336 2014-03-25 15:07:20Z raasch
# -write_compressed
#
# 2014-03-20 16:38:49Z raasch
# cpulog-dependency added for disturb_heatflux
#
# 1320 2014-03-20 08:40:49Z raasch
# +mod_kinds in dependency list
#
# 1318 2014-03-17 13:35:16Z raasch
# +cpulog, dependencies for cpulog added
# -cpu_log, cpu_statistics
#
# 1311 2014-03-14 12:13:39Z heinze
# bugfix: buoyancy dependency added to nudging
#
# 1306 2014-03-13 14:30:59Z raasch
# +mod_kinds
#
# 1241 2013-10-30 11:36:58Z heinze
# +nudging, ls_forcing
#
# 1212 2013-08-15 08:46:27Z raasch
# +tridia_solver, -poisfft_hybrid
#
# 1210 2013-08-14 10:58:20Z raasch
# object file list replaced by one line statement
#
# 1128 2013-04-12 06:19:32Z raasch
# dependencies for time_integration modified
#
# 1111 2013-03-08 23:54:10Z raasch
# dependencies removed from init_pegrid
# bugfix: dependency added for cuda_fft_interfaces
#
# 1106 2013-03-04 05:31:38Z raasch
# +cuda_fft_interfaces
#
# 1053 2012-11-13 17:11:03Z hoffmann
# +microphysics for two-moment cloud physics
#
# 1036 2012-10-22 13:43:42Z raasch
# code put under GPL (PALM 3.9)
#
# 1001 2012-09-13 14:08:46Z raasch
# -asselin_filter, advec_s|u|v|w_ups, spline_x|y|z
#
# 872 2012-03-29 15:57:54Z franke
# lpm depends on modules
#
# 849 2012-03-15 10:35:09Z raasch
# +lpm_advec, lpm_data_output_particles, lpm_droplet_collision,
# lpm_droplet_condensation,
# lpm_exchange_horiz, lpm_extend_particle_array, lpm_extend_tails,
# lpm_extend_tail_array, lpm_init_sgs_tke, lpm_pack_arrays,
# lpm_read_restart_file, lpm_release_set, lpm_sort_arrays,
# lpm_write_exchange_statistics, lpm_write_restart_file
# advec_particles renamed lpm,
# init_particles renamed lpm_init,
# user_advec_particles renamed user_lpm_advec,
# particle_boundary_conds renamed lpm_boundary_conds,
# set_particle_attributes renamed lpm_set_attributes
# user_init_particles renamed user_lpm_init,
# user_particle_attributes renamed user_lpm_set_attributes
#
# 828 2012-02-21 12:00:36Z raasch
# init_particles depends on lpm_collision_kernels
#
# 825 2012-02-19 03:03:44Z raasch
# wang_kernel renamed lpm_collision_kernels
#
# 790 2011-11-29 03:11:20Z raasch
# +wang_kernel
#
# 667 2010-12-23 12:06:00Z suehring/gryschka
# +advec_ws
#
# 421 2010-01-25 15:50:10Z raasch
# comment line for simple install method added, which will be replaced
# by an include command for preprocessor and compiler options when the simple
# install process is used instead of mbuild/mrun,
# +set_particle_attributes, set_slicer_attributes_dvrp
# +subsidence
#
# 410 2009-12-04 17:05:40Z letzel
# +data_output_mask, init_masks, user_data_output_mask
#
# 217 2008-12-09 18:00:48Z letzel
# former file user_interface split into one file per subroutine
# +init_coupling, message
#
# 151 2008-03-07 13:42:18Z raasch
# +plant_canopy_model, inflow_turbulence
# +surface_coupler
#
# 96 2007-06-04 08:07:41Z raasch
# +eqn_state_seawater, init_ocean
#
# 82 2007-04-16 15:40:52Z raasch
# +local_flush
#
# 58 2007-03-09 14:27:38Z raasch
# default suffixes removed from the suffix list to avoid calling of m2c in
# case of .mod files
# +wall_fluxes, +particle_boundary_conds
#
# RCS Log replace by Id keyword, revision history cleaned up
#
# Revision 1.13  2006/08/04 14:45:31  raasch
# +data_output_ptseries
#
# Revision 1.1  2002/06/11 13:42:20  raasch
# Initial revision
#
#
# Description:
# ------------
# Makefile for the parallelized LES model (PALM) to be used by the
# shellscript "mrun"
#------------------------------------------------------------------------------!

PROG =  palm

SOURCES = \
	advec_s_bc.f90 \
	advec_s_pw.f90 \
	advec_s_up.f90 \
	advec_u_pw.f90 \
	advec_u_up.f90 \
	advec_v_pw.f90 \
	advec_v_up.f90 \
	advec_w_pw.f90 \
	advec_ws.f90 \
	advec_w_up.f90 \
	average_3d_data.f90 \
	boundary_conds.f90 \
	buoyancy.f90 \
	calc_liquid_water_content.f90 \
	calc_mean_profile.f90 \
	calc_radiation.f90 \
	check_for_restart.f90 \
	check_open.f90 \
	check_parameters.f90 \
	chem_gasphase_mod.f90 \
	chemistry_model_mod.f90 \
	chem_modules.f90 \
	chem_photolysis_mod.f90 \
	close_file.f90 \
	compute_vpt.f90 \
	coriolis.f90 \
	cpulog_mod.f90 \
	data_log.f90 \
	data_output_2d.f90 \
	data_output_3d.f90 \
	data_output_dvrp.f90 \
	data_output_flight.f90\
	data_output_mask.f90 \
	data_output_profiles.f90 \
	data_output_ptseries.f90 \
	data_output_spectra.f90 \
	data_output_tseries.f90 \
	date_and_time_mod.f90 \
	diagnostic_quantities_mod.f90 \
	diffusion_s.f90 \
	diffusion_u.f90 \
	diffusion_v.f90 \
	diffusion_w.f90 \
	disturb_field.f90 \
	disturb_heatflux.f90 \
	eqn_state_seawater.f90 \
	exchange_horiz_2d.f90 \
	exchange_horiz.f90 \
	fft_xy_mod.f90 \
	flow_statistics.f90 \
	global_min_max.f90 \
	gust_mod.f90 \
	header.f90 \
	inflow_turbulence.f90 \
	init_3d_model.f90 \
	init_advec.f90 \
	init_cloud_physics.f90 \
	init_coupling.f90 \
	init_dvrp.f90 \
	init_grid.f90 \
	init_masks.f90 \
	init_ocean.f90 \
	init_pegrid.f90 \
	init_pt_anomaly.f90 \
	init_rankine.f90 \
	init_slope.f90 \
	interaction_droplets_ptq.f90 \
	land_surface_model_mod.f90 \
	large_scale_forcing_nudging_mod.f90 \
	local_stop.f90 \
	local_system.f90 \
	local_tremain.f90 \
	local_tremain_ini.f90 \
	lpm_advec.f90 \
	lpm_boundary_conds.f90 \
	lpm_calc_liquid_water_content.f90 \
	lpm_collision_kernels.f90 \
	lpm_data_output_particles.f90 \
	lpm_droplet_collision.f90 \
	lpm_droplet_condensation.f90 \
	lpm_exchange_horiz.f90 \
	lpm.f90 \
	lpm_init.f90 \
	lpm_init_sgs_tke.f90 \
	lpm_merging.f90 \
	lpm_pack_arrays.f90 \
	lpm_read_restart_file.f90 \
	lpm_set_attributes.f90 \
	lpm_splitting.f90 \
	lpm_write_exchange_statistics.f90 \
	lpm_write_restart_file.f90 \
	message.f90 \
	microphysics_mod.f90 \
	model_1d_mod.f90 \
	mod_kinds.f90 \
	mod_particle_attributes.f90 \
	modules.f90 \
	netcdf_data_input_mod.f90 \
	netcdf_interface_mod.f90 \
	outflow_turbulence.f90 \
	package_parin.f90 \
	palm.f90 \
	parin.f90 \
	plant_canopy_model_mod.f90 \
	pmc_child_mod.f90 \
	pmc_general_mod.f90 \
	pmc_handle_communicator_mod.f90 \
	pmc_interface_mod.f90 \
	pmc_mpi_wrapper_mod.f90 \
	pmc_parent_mod.f90 \
	pmc_particle_interface.f90 \
	poisfft_mod.f90 \
	poismg_mod.f90 \
	poismg_noopt_mod.f90 \
	posix_calls_from_fortran.f90 \
	pres.f90 \
	print_1d.f90 \
	prognostic_equations.f90 \
	progress_bar_mod.f90 \
	radiation_model_mod.f90 \
	random_function_mod.f90 \
	random_gauss.f90 \
	random_generator_parallel_mod.f90 \
        read_restart_data_mod.f90 \
	run_control.f90 \
	set_slicer_attributes_dvrp.f90 \
	singleton_mod.f90 \
	sor.f90 \
	spectra_mod.f90 \
	stokes_drift_mod.f90 \
	stokes_force_mod.f90 \
	subsidence_mod.f90 \
	sum_up_3d_data.f90 \
	surface_coupler.f90 \
	surface_layer_fluxes_mod.f90 \
	surface_mod.f90 \
	swap_timelevel.f90 \
	synthetic_turbulence_generator_mod.f90 \
	temperton_fft_mod.f90 \
	time_integration.f90 \
	time_integration_spinup.f90 \
	timestep.f90 \
	timestep_scheme_steering.f90 \
	time_to_string.f90 \
	transpose.f90 \
	tridia_solver_mod.f90 \
	turbulence_closure_mod.f90 \
	urban_surface_mod.f90 \
	user_3d_data_averaging.f90 \
	user_actions.f90 \
	user_additional_routines.f90 \
	user_check_data_output.f90 \
	user_check_data_output_pr.f90 \
	user_check_parameters.f90 \
	user_data_output_2d.f90 \
	user_data_output_3d.f90 \
	user_data_output_dvrp.f90 \
	user_data_output_mask.f90 \
	user_define_netcdf_grid.f90 \
	user_dvrp_coltab.f90 \
	user_flight.f90\
	user_header.f90 \
	user_init_3d_model.f90 \
	user_init.f90 \
	user_init_flight.f90 \
	user_init_grid.f90 \
	user_init_land_surface.f90 \
	user_init_plant_canopy.f90 \
	user_init_radiation.f90 \
	user_init_urban_surface.f90 \
	user_last_actions.f90 \
	user_lpm_advec.f90 \
	user_lpm_init.f90 \
	user_lpm_set_attributes.f90 \
	user_module.f90 \
	user_parin.f90 \
	user_read_restart_data_mod.f90 \
	user_spectra.f90 \
	user_statistics.f90 \
	user_write_restart_data_mod.f90 \
	uv_exposure_model_mod.f90 \
	vertical_nesting_mod.f90 \
	virtual_flight_mod.f90 \
	wind_turbine_model_mod.f90 \
	wrd_write_string.f90 \
        write_restart_data_mod.f90


OBJS=$(SOURCES:.f90=.o)

CC = cc
CFLAGS = -O

F90 =
COPT =
F90FLAGS =
LDFLAGS =


# The following line is needed for palm_simple_install, don't remove it!
#to_be_replaced_by_include


.SUFFIXES:
.SUFFIXES: .o .f90


all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) -o $(PROG) $(OBJS) $(LDFLAGS)

clean:
	rm -f $(PROG) $(OBJS) *.mod *.i *.lst

.f90.o:
	$(F90) $(F90FLAGS) $(COPT) -c $<


advec_s_bc.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
advec_s_pw.o: \
	mod_kinds.o \
	modules.o
advec_s_up.o: \
	mod_kinds.o \
	modules.o
advec_u_pw.o: \
	mod_kinds.o \
	modules.o
advec_u_up.o: \
	mod_kinds.o \
	modules.o
advec_v_pw.o: \
	mod_kinds.o \
	modules.o
advec_v_up.o: \
	mod_kinds.o \
	modules.o
advec_ws.o: \
	mod_kinds.o \
	modules.o
advec_w_pw.o: \
	mod_kinds.o \
	modules.o
advec_w_up.o: \
	mod_kinds.o \
	modules.o
average_3d_data.o: \
	chemistry_model_mod.o \
	cpulog_mod.o \
	exchange_horiz_2d.o \
	gust_mod.o \
	land_surface_model_mod.o \
	mod_kinds.o \
	modules.o \
	radiation_model_mod.o \
	turbulence_closure_mod.o \
	urban_surface_mod.o
boundary_conds.o: \
	chemistry_model_mod.o \
	mod_kinds.o \
	modules.o \
	pmc_interface_mod.o \
	surface_mod.o \
	turbulence_closure_mod.o
buoyancy.o: \
	mod_kinds.o \
	modules.o
calc_mean_profile.o: \
	mod_kinds.o \
	modules.o
calc_liquid_water_content.o: \
	diagnostic_quantities_mod.o \
	mod_kinds.o \
	modules.o
calc_radiation.o: \
	mod_kinds.o \
	modules.o
check_for_restart.o: \
	mod_kinds.o \
	modules.o \
	pmc_interface_mod.o
check_open.o: \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	netcdf_interface_mod.o \
	posix_calls_from_fortran.o
check_parameters.o: \
	chemistry_model_mod.o \
	gust_mod.o \
	land_surface_model_mod.o \
	large_scale_forcing_nudging_mod.o \
	microphysics_mod.o \
	model_1d_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	netcdf_interface_mod.o \
	plant_canopy_model_mod.o \
	pmc_interface_mod.o \
	radiation_model_mod.o \
	spectra_mod.o \
	subsidence_mod.o \
	synthetic_turbulence_generator_mod.o \
	turbulence_closure_mod.o \
	urban_surface_mod.o \
	uv_exposure_model_mod.o \
	vertical_nesting_mod.o \
	wind_turbine_model_mod.o
chemistry_model_mod.o: \
	chem_gasphase_mod.o \
	chem_modules.o \
	diffusion_s.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	surface_mod.o
chem_gasphase_mod.o: \
	mod_kinds.o \
	modules.o
chem_modules.o: \
	chem_gasphase_mod.o \
	mod_kinds.o \
	modules.o
chem_photolysis_mod.o: \
	chem_gasphase_mod.o \
	chemistry_model_mod.o \
	mod_kinds.o \
	modules.o \
	radiation_model_mod.o
close_file.o: \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o
compute_vpt.o: \
	mod_kinds.o \
	modules.o
coriolis.o: \
	mod_kinds.o \
	modules.o
cpulog_mod.o: \
	mod_kinds.o \
	modules.o
cpu_statistics.o: \
	mod_kinds.o \
	modules.o
data_log.o: \
	mod_kinds.o \
	modules.o
data_output_dvrp.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
data_output_mask.o: \
	cpulog_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	netcdf_interface_mod.o
data_output_profiles.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o
data_output_ptseries.o: \
	cpulog_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	netcdf_interface_mod.o
data_output_spectra.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	spectra_mod.o
data_output_tseries.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o
data_output_flight.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	virtual_flight_mod.o
data_output_2d.o: \
	cpulog_mod.o \
	gust_mod.o \
	land_surface_model_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	netcdf_interface_mod.o \
	radiation_model_mod.o \
	surface_mod.o \
	urban_surface_mod.o \
	uv_exposure_model_mod.o
data_output_3d.o: \
	chemistry_model_mod.o \
	cpulog_mod.o \
	gust_mod.o \
	land_surface_model_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	netcdf_interface_mod.o \
	turbulence_closure_mod.o \
	urban_surface_mod.o
date_and_time_mod.o: \
	mod_kinds.o \
	modules.o
diagnostic_quantities_mod.o: \
	mod_kinds.o \
	modules.o
diffusion_s.o: \
	mod_kinds.o \
	modules.o \
	surface_mod.o
diffusion_u.o: \
	mod_kinds.o \
	modules.o \
	surface_mod.o
diffusion_v.o: \
	mod_kinds.o \
	modules.o \
	surface_mod.o
diffusion_w.o: \
	mod_kinds.o \
	modules.o \
	surface_mod.o
disturb_field.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	random_function_mod.o \
	random_generator_parallel_mod.o
disturb_heatflux.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	random_generator_parallel_mod.o \
	surface_mod.o
eqn_state_seawater.o: \
	mod_kinds.o \
	modules.o \
	surface_mod.o
exchange_horiz.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
exchange_horiz_2d.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	pmc_interface_mod.o
fft_xy_mod.o: \
	mod_kinds.o \
	modules.o \
	singleton_mod.o \
	temperton_fft_mod.o
flow_statistics.o: \
	cpulog_mod.o \
	gust_mod.o \
	land_surface_model_mod.o \
	large_scale_forcing_nudging_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	radiation_model_mod.o \
	surface_mod.o
global_min_max.o: \
	mod_kinds.o \
	modules.o
gust_mod.o: \
	mod_kinds.o \
	modules.o
header.o: \
	cpulog_mod.o \
	date_and_time_mod.o \
	gust_mod.o \
	land_surface_model_mod.o \
	large_scale_forcing_nudging_mod.o \
	microphysics_mod.o \
	model_1d_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	plant_canopy_model_mod.o \
	pmc_handle_communicator_mod.o \
	pmc_interface_mod.o \
	radiation_model_mod.o \
	spectra_mod.o \
	subsidence_mod.o \
	surface_mod.o \
	virtual_flight_mod.o
inflow_turbulence.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
init_3d_model.o: \
	advec_ws.o \
	chemistry_model_mod.o \
	cpulog_mod.o \
	disturb_heatflux.o \
	gust_mod.o \
	land_surface_model_mod.o \
	large_scale_forcing_nudging_mod.o \
	lpm_init.o \
	microphysics_mod.o \
	model_1d_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	netcdf_data_input_mod.o \
	netcdf_interface_mod.o \
	plant_canopy_model_mod.o \
	pmc_interface_mod.o \
	radiation_model_mod.o \
	random_function_mod.o \
	random_generator_parallel_mod.o \
	read_restart_data_mod.o \
	surface_layer_fluxes_mod.o \
	surface_mod.o \
	turbulence_closure_mod.o \
	urban_surface_mod.o \
	uv_exposure_model_mod.o \
	virtual_flight_mod.o \
	wind_turbine_model_mod.o
init_advec.o: \
	mod_kinds.o \
	modules.o
init_cloud_physics.o: \
	mod_kinds.o \
	modules.o
init_coupling.o: \
	mod_kinds.o \
	modules.o \
	vertical_nesting_mod.o
init_dvrp.o: \
	mod_kinds.o \
	modules.o
init_grid.o: \
	advec_ws.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	netcdf_interface_mod.o \
	poismg_noopt_mod.o \
	surface_mod.o \
	vertical_nesting_mod.o
init_masks.o: \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o
init_ocean.o: \
	eqn_state_seawater.o \
	mod_kinds.o \
	modules.o \
	stokes_drift_mod.o
init_pegrid.o: \
	mod_kinds.o \
	modules.o \
	spectra_mod.o \
	synthetic_turbulence_generator_mod.o \
	vertical_nesting_mod.o
init_pt_anomaly.o: \
	mod_kinds.o \
	modules.o
init_rankine.o: \
	mod_kinds.o \
	modules.o
init_slope.o: \
	eqn_state_seawater.o \
	mod_kinds.o \
	modules.o
interaction_droplets_ptq.o: \
	mod_kinds.o \
	modules.o
land_surface_model_mod.o: \
	calc_mean_profile.o \
	diagnostic_quantities_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	pmc_interface_mod.o \
	radiation_model_mod.o \
	surface_mod.o
large_scale_forcing_nudging_mod.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	surface_mod.o
local_stop.o: \
	mod_kinds.o \
	modules.o \
	pmc_interface_mod.o \
	vertical_nesting_mod.o
local_tremain.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
local_tremain_ini.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
lpm.o: \
	cpulog_mod.o \
	lpm_exchange_horiz.o \
	lpm_init.o \
	lpm_merging.o \
	lpm_pack_arrays.o \
	lpm_splitting.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	pmc_interface_mod.o \
	pmc_particle_interface.o
lpm_advec.o: \
	cpulog_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	surface_mod.o
lpm_boundary_conds.o: \
	cpulog_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	surface_mod.o
lpm_calc_liquid_water_content.o: \
	cpulog_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o
lpm_collision_kernels.o: \
	cpulog_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	user_module.o
lpm_data_output_particles.o: \
	cpulog_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	netcdf_interface_mod.o
lpm_droplet_collision.o: \
	cpulog_mod.o \
	lpm_collision_kernels.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o
lpm_droplet_condensation.o: \
	cpulog_mod.o \
	diagnostic_quantities_mod.o \
	lpm_collision_kernels.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o
lpm_exchange_horiz.o: \
	cpulog_mod.o \
	lpm_pack_arrays.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	netcdf_interface_mod.o
lpm_init.o: \
	diagnostic_quantities_mod.o \
	lpm_collision_kernels.o \
	lpm_exchange_horiz.o \
	lpm_pack_arrays.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	netcdf_interface_mod.o \
	pmc_particle_interface.o \
	random_function_mod.o \
	random_generator_parallel_mod.o \
	surface_mod.o
lpm_init_sgs_tke.o: \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	surface_mod.o
lpm_merging.o: \
	cpulog_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o
lpm_pack_arrays.o: \
	cpulog_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o
lpm_read_restart_file.o: \
	lpm_pack_arrays.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o
lpm_set_attributes.o: \
	cpulog_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o
lpm_splitting.o: \
	cpulog_mod.o \
	lpm_exchange_horiz.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o
lpm_write_exchange_statistics.o: \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	pmc_particle_interface.o
lpm_write_restart_file.o: \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o
message.o: \
	mod_kinds.o \
	modules.o \
	pmc_interface_mod.o
microphysics_mod.o: \
	cpulog_mod.o \
	diagnostic_quantities_mod.o \
	mod_kinds.o \
	modules.o \
	surface_mod.o
model_1d_mod.o: \
	mod_kinds.o \
	modules.o
modules.o: \
	mod_kinds.o
mod_kinds.o: \

mod_particle_attributes.o: \
	mod_kinds.o
netcdf_data_input_mod.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o   \
	surface_mod.o
netcdf_interface_mod.o: \
	chemistry_model_mod.o \
	gust_mod.o \
	land_surface_model_mod.o \
	mod_kinds.o \
	modules.o \
	radiation_model_mod.o \
	spectra_mod.o \
	turbulence_closure_mod.o \
	urban_surface_mod.o \
	uv_exposure_model_mod.o
outflow_turbulence.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
package_parin.o: \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o
palm.o: \
	chemistry_model_mod.o \
	chem_photolysis_mod.o \
	cpulog_mod.o \
	land_surface_model_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	pmc_interface_mod.o \
	pmc_particle_interface.o \
	surface_layer_fluxes_mod.o \
	write_restart_data_mod.o
parin.o: \
	chemistry_model_mod.o \
	chem_modules.o \
	cpulog_mod.o \
	date_and_time_mod.o \
	gust_mod.o \
	land_surface_model_mod.o \
	microphysics_mod.o \
	model_1d_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	plant_canopy_model_mod.o \
	pmc_interface_mod.o \
	progress_bar_mod.o \
	radiation_model_mod.o \
	read_restart_data_mod.o \
	spectra_mod.o \
	synthetic_turbulence_generator_mod.o \
	turbulence_closure_mod.o \
	uv_exposure_model_mod.o \
	vertical_nesting_mod.o \
	virtual_flight_mod.o \
	wind_turbine_model_mod.o
plant_canopy_model_mod.o: \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	surface_mod.o
pmc_particle_interface.o: \
	lpm_exchange_horiz.o \
	lpm_pack_arrays.o \
	mod_particle_attributes.o \
	modules.o \
	pmc_child_mod.o \
	pmc_interface_mod.o \
	pmc_parent_mod.o
pmc_interface_mod.o: \
	chemistry_model_mod.o \
	chem_modules.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	pmc_child_mod.o \
	pmc_general_mod.o \
	pmc_handle_communicator_mod.o \
	pmc_mpi_wrapper_mod.o \
	pmc_parent_mod.o \
	surface_mod.o
pmc_child_mod.o: \
	mod_kinds.o \
	pmc_general_mod.o \
	pmc_handle_communicator_mod.o \
	pmc_mpi_wrapper_mod.o
pmc_general_mod.o: \
	mod_kinds.o
pmc_handle_communicator_mod.o: \
	mod_kinds.o \
	modules.o \
	pmc_general_mod.o
pmc_mpi_wrapper_mod.o: \
	pmc_handle_communicator_mod.o
pmc_parent_mod.o: \
	pmc_general_mod.o \
	pmc_handle_communicator_mod.o \
	pmc_mpi_wrapper_mod.o
poisfft_mod.o: \
	cpulog_mod.o \
	fft_xy_mod.o \
	mod_kinds.o \
	modules.o \
	tridia_solver_mod.o
poismg_mod.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
poismg_noopt_mod.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
posix_calls_from_fortran.o: \

pres.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	pmc_interface_mod.o \
	poisfft_mod.o \
	poismg_mod.o \
	poismg_noopt_mod.o \
	surface_mod.o
print_1d.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
prognostic_equations.o: \
	advec_s_bc.o \
	advec_s_pw.o \
	advec_s_up.o \
	advec_u_pw.o \
	advec_u_up.o \
	advec_v_pw.o \
	advec_v_up.o \
	advec_w_pw.o \
	advec_ws.o \
	advec_w_up.o \
	buoyancy.o \
	calc_radiation.o \
	chemistry_model_mod.o \
	chem_modules.o \
	chem_photolysis_mod.o \
	coriolis.o \
	cpulog_mod.o \
	cpulog_mod.o \
	diffusion_s.o \
	diffusion_u.o \
	diffusion_v.o \
	diffusion_w.o \
	eqn_state_seawater.o \
	large_scale_forcing_nudging_mod.o \
	microphysics_mod.o \
	mod_kinds.o \
	mod_kinds.o \
	modules.o \
	modules.o \
	plant_canopy_model_mod.o \
	radiation_model_mod.o \
	stokes_force_mod.o \
	subsidence_mod.o \
	surface_mod.o \
	surface_mod.o \
	turbulence_closure_mod.o \
	user_actions.o \
	wind_turbine_model_mod.o
progress_bar_mod.o: \
	mod_kinds.o \
	modules.o
radiation_model_mod.o: \
	calc_mean_profile.o \
	date_and_time_mod.o \
	microphysics_mod.o \
	mod_particle_attributes.o \
	modules.o \
	netcdf_data_input_mod.o \
	plant_canopy_model_mod.o \
	surface_mod.o
random_function_mod.o: \
	mod_kinds.o
random_gauss.o: \
	mod_kinds.o \
	random_function_mod.o \
	random_generator_parallel_mod.o
random_generator_parallel_mod.o: \
	mod_kinds.o \
	modules.o
read_restart_data_mod.o: \
	chemistry_model_mod.o \
	cpulog_mod.o \
	date_and_time_mod.o \
	gust_mod.o \
	land_surface_model_mod.o \
	microphysics_mod.o \
	modules.o \
	mod_kinds.o \
	model_1d_mod.o \
	netcdf_interface_mod.o \
	radiation_model_mod.o \
	random_function_mod.o\
	random_generator_parallel_mod.o \
	spectra_mod.o \
	surface_mod.o \
	synthetic_turbulence_generator_mod.o \
	urban_surface_mod.o \
	user_read_restart_data_mod.o \
	vertical_nesting_mod.o \
	virtual_flight_mod.o \
	wind_turbine_model_mod.o
run_control.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
set_slicer_attributes_dvrp.o: \
	mod_kinds.o \
	modules.o
singleton_mod.o: \
	mod_kinds.o
sor.o: \
	mod_kinds.o \
	modules.o
spectra_mod.o: \
	cpulog_mod.o \
	fft_xy_mod.o \
	mod_kinds.o \
	modules.o
stokes_drift_mod.o: \
	mod_kinds.o \
	modules.o
stokes_force_mod.o: \
	mod_kinds.o \
	modules.o
subsidence_mod.o: \
	mod_kinds.o \
	modules.o
sum_up_3d_data.o: \
	chemistry_model_mod.o \
	cpulog_mod.o \
	gust_mod.o \
	land_surface_model_mod.o \
	mod_kinds.o \
	modules.o \
	radiation_model_mod.o \
	surface_mod.o \
	turbulence_closure_mod.o \
	urban_surface_mod.o \
	uv_exposure_model_mod.o
surface_coupler.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	surface_mod.o
surface_layer_fluxes_mod.o: \
	chem_modules.o \
    eqn_state_seawater.o \
	exchange_horiz_2d.o \
	land_surface_model_mod.o \
	mod_kinds.o \
	modules.o \
	surface_mod.o \
	urban_surface_mod.o
surface_mod.o: \
	chem_modules.o \
	model_1d_mod.o \
	mod_kinds.o \
	modules.o
swap_timelevel.o: \
	chemistry_model_mod.o \
	cpulog_mod.o \
	gust_mod.o \
	land_surface_model_mod.o \
	mod_kinds.o \
	modules.o \
	pmc_interface_mod.o \
	turbulence_closure_mod.o \
	urban_surface_mod.o
synthetic_turbulence_generator_mod.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	pmc_interface_mod.o
temperton_fft_mod.o: \
	mod_kinds.o \
	modules.o
time_integration.o: \
	advec_ws.o \
	buoyancy.o \
	calc_mean_profile.o \
	chemistry_model_mod.o \
	chem_modules.o \
	cpulog_mod.o \
	data_output_flight.o \
	disturb_heatflux.o \
	gust_mod.o \
	interaction_droplets_ptq.o \
	land_surface_model_mod.o \
	large_scale_forcing_nudging_mod.o \
	lpm.o \
	microphysics_mod.o \
	mod_kinds.o \
	modules.o \
	pmc_interface_mod.o \
	prognostic_equations.o \
	progress_bar_mod.o \
	radiation_model_mod.o \
	spectra_mod.o \
	stokes_force_mod.o \
	surface_layer_fluxes_mod.o \
	surface_mod.o \
	synthetic_turbulence_generator_mod.o \
	turbulence_closure_mod.o\
	urban_surface_mod.o \
	user_actions.o \
	uv_exposure_model_mod.o \
	vertical_nesting_mod.o \
	virtual_flight_mod.o \
	wind_turbine_model_mod.o
time_integration_spinup.o: \
	cpulog_mod.o \
	date_and_time_mod.o \
	land_surface_model_mod.o \
	mod_kinds.o \
	modules.o \
	pmc_interface_mod.o \
	radiation_model_mod.o \
	surface_layer_fluxes_mod.o \
	urban_surface_mod.o
time_to_string.o: \
	mod_kinds.o
timestep.o: \
	cpulog_mod.o \
	microphysics_mod.o \
	mod_kinds.o \
	modules.o \
	pmc_interface_mod.o \
	vertical_nesting_mod.o
timestep_scheme_steering.o: \
	mod_kinds.o \
	modules.o
transpose.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
tridia_solver_mod.o: \
	mod_kinds.o \
	modules.o
turbulence_closure_mod.o: \
	advec_s_bc.o \
	advec_s_pw.o \
	advec_s_up.o \
	advec_ws.o \
	cpulog_mod.o \
	microphysics_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	plant_canopy_model_mod.o \
	pmc_interface_mod.o \
	stokes_drift_mod.o \
	surface_mod.o \
	user_actions.o
urban_surface_mod.o: \
	date_and_time_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	plant_canopy_model_mod.o \
	radiation_model_mod.o \
	surface_mod.o
user_3d_data_averaging.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_actions.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	user_module.o
user_additional_routines.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_check_data_output.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_check_data_output_pr.o: \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	user_module.o
user_check_parameters.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_data_output_2d.o: \
	mod_kinds.o \
	modules.o \
	surface_mod.o \
	user_module.o
user_data_output_3d.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_data_output_mask.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_data_output_dvrp.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_define_netcdf_grid.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_dvrp_coltab.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_flight.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_header.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_init.o: \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	user_module.o
user_init_3d_model.o: \
	mod_kinds.o \
	modules.o \
	surface_mod.o \
	user_module.o
user_init_flight.o: \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	user_module.o
user_init_grid.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_init_land_surface.o: \
	land_surface_model_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	surface_mod.o \
	user_module.o
user_init_plant_canopy.o: \
	mod_kinds.o \
	modules.o \
	plant_canopy_model_mod.o \
	user_module.o
user_init_radiation.o: \
	mod_kinds.o \
	modules.o \
	radiation_model_mod.o \
	user_module.o
user_init_urban_surface.o: \
	mod_kinds.o \
	modules.o \
	surface_mod.o \
	urban_surface_mod.o \
	user_module.o
user_last_actions.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_lpm_advec.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_lpm_init.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_lpm_set_attributes.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_module.o: \
	mod_kinds.o
user_parin.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_read_restart_data_mod.o: \
	modules.o \
	user_module.o
user_spectra.o: \
	mod_kinds.o \
	modules.o \
	spectra_mod.o \
	user_module.o
user_statistics.o: \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	user_module.o
user_write_restart_data_mod.o: \
	user_module.o
uv_exposure_model_mod.o: \
	date_and_time_mod.o \
	mod_kinds.o \
	modules.o \
	radiation_model_mod.o
vertical_nesting_mod.o: \
	mod_kinds.o \
	modules.o \
	surface_mod.o \
	turbulence_closure_mod.o
virtual_flight_mod.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	user_flight.o \
	user_init_flight.o
wind_turbine_model_mod.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
wrd_write_string.o: \
	mod_kinds.o
write_restart_data_mod.o: \
	chemistry_model_mod.o \
	date_and_time_mod.o \
	gust_mod.o \
	land_surface_model_mod.o \
	microphysics_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	model_1d_mod.o \
	modules.o \
	netcdf_interface_mod.o \
	radiation_model_mod.o \
	random_function_mod.o \
	random_generator_parallel_mod.o \
	spectra_mod.o surface_mod.o \
	synthetic_turbulence_generator_mod.o \
	urban_surface_mod.o \
	user_write_restart_data_mod.o \
	vertical_nesting_mod.o \
	virtual_flight_mod.o \
	wind_turbine_model_mod.o
