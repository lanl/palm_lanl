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
#
# Description:
# ------------
# Makefile for the parallelized LES model (PALM) to be used by the
# shellscript "mrun"
#------------------------------------------------------------------------------!

OS := $(shell uname -s)

PROG =  palm

SOURCES = \
	advec_ws.f90 \
	average_3d_data.f90 \
	boundary_conds.f90 \
	buoyancy.f90 \
	calc_mean_profile.f90 \
	check_for_restart.f90 \
	check_open.f90 \
	check_parameters.f90 \
	close_file.f90 \
	compute_vpt.f90 \
	coriolis.f90 \
	cpulog_mod.f90 \
	data_log.f90 \
	data_output_3d.f90 \
	data_output_profiles.f90 \
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
	header.f90 \
	inflow_turbulence.f90 \
	init_3d_model.f90 \
	init_advec.f90 \
	init_grid.f90 \
	init_ocean.f90 \
	init_pegrid.f90 \
	init_pt_anomaly.f90 \
	init_rankine.f90 \
	init_slope.f90 \
	local_stop.f90 \
	local_system.f90 \
	local_tremain.f90 \
	local_tremain_ini.f90 \
	message.f90 \
	mod_kinds.f90 \
	modules.f90 \
	netcdf_data_input_mod.f90 \
	netcdf_interface_mod.f90 \
	outflow_turbulence.f90 \
	palm.f90 \
	parin.f90 \
	poisfft_mod.f90 \
	poismg_mod.f90 \
	poismg_noopt_mod.f90 \
	posix_calls_from_fortran.f90 \
	pres.f90 \
	prognostic_equations.f90 \
	progress_bar_mod.f90 \
	random_function_mod.f90 \
	random_gauss.f90 \
	random_generator_parallel_mod.f90 \
        read_restart_data_mod.f90 \
	run_control.f90 \
	singleton_mod.f90 \
	sor.f90 \
	stokes_drift_mod.f90 \
	stokes_force_mod.f90 \
	sum_up_3d_data.f90 \
	surface_coupler.f90 \
	surface_layer_fluxes_mod.f90 \
	surface_mod.f90 \
	swap_timelevel.f90 \
	temperton_fft_mod.f90 \
	time_integration.f90 \
	time_integration_spinup.f90 \
	timestep.f90 \
	timestep_scheme_steering.f90 \
	time_to_string.f90 \
	transpose.f90 \
	tridia_solver_mod.f90 \
	turbulence_closure_mod.f90 \
	wrd_write_string.f90 \
        write_restart_data_mod.f90


OBJS=$(SOURCES:.f90=.o)

CC = cc
CFLAGS = -O

include MAKE.inc



.SUFFIXES:
.SUFFIXES: .o .f90


all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) -o $(PROG) $(OBJS) $(LDFLAGS)

libs: libLES.a

libLES.a: $(OBJS)
ifeq ($(OS),Linux)
	ar cr -o libLES.a *.o 
endif
ifeq ($(OS),Darwin)
	libtool -static -o libLES.a *.o
endif
	ranlib libLES.a

clean:
	rm -f libLES.a $(PROG) $(OBJS) *.mod *.i *.lst

.f90.o:
	$(F90) $(F90FLAGS) $(COPT) -c $<


advec_ws.o: \
	mod_kinds.o \
	modules.o
average_3d_data.o: \
	cpulog_mod.o \
	exchange_horiz_2d.o \
	mod_kinds.o \
	modules.o \
	turbulence_closure_mod.o
boundary_conds.o: \
	mod_kinds.o \
	modules.o \
	surface_mod.o \
	turbulence_closure_mod.o
buoyancy.o: \
	mod_kinds.o \
	modules.o
calc_mean_profile.o: \
	mod_kinds.o \
	modules.o
check_for_restart.o: \
	mod_kinds.o \
	modules.o
check_open.o: \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	posix_calls_from_fortran.o
check_parameters.o: \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	netcdf_interface_mod.o \
	turbulence_closure_mod.o
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
data_output_profiles.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o
data_output_tseries.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o
data_output_3d.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	turbulence_closure_mod.o
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
	modules.o
fft_xy_mod.o: \
	mod_kinds.o \
	modules.o \
	singleton_mod.o \
	temperton_fft_mod.o
flow_statistics.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	surface_mod.o
global_min_max.o: \
	mod_kinds.o \
	modules.o
header.o: \
	cpulog_mod.o \
	date_and_time_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	surface_mod.o
inflow_turbulence.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
init_3d_model.o: \
	advec_ws.o \
	cpulog_mod.o \
	disturb_heatflux.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	netcdf_interface_mod.o \
	random_function_mod.o \
	random_generator_parallel_mod.o \
	read_restart_data_mod.o \
	surface_layer_fluxes_mod.o \
	surface_mod.o \
	turbulence_closure_mod.o
init_advec.o: \
	mod_kinds.o \
	modules.o
init_coupling.o: \
	mod_kinds.o \
	modules.o
init_grid.o: \
	advec_ws.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	netcdf_interface_mod.o \
	poismg_noopt_mod.o \
	surface_mod.o
init_ocean.o: \
	eqn_state_seawater.o \
	mod_kinds.o \
	modules.o \
	stokes_drift_mod.o
init_pegrid.o: \
	mod_kinds.o \
	modules.o
init_pt_anomaly.o: \
	mod_kinds.o \
	modules.o
init_rankine.o: \
	mod_kinds.o \
	modules.o
init_slope.o: \
	mod_kinds.o \
	modules.o
local_stop.o: \
	mod_kinds.o \
	modules.o
local_tremain.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
local_tremain_ini.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
message.o: \
	mod_kinds.o \
	modules.o
modules.o: \
	mod_kinds.o
mod_kinds.o: \

netcdf_data_input_mod.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o   \
	surface_mod.o
netcdf_interface_mod.o: \
	mod_kinds.o \
	turbulence_closure_mod.o \
	modules.o
outflow_turbulence.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
palm.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	surface_layer_fluxes_mod.o \
	write_restart_data_mod.o \
	init_3d_model.o \
	tridia_solver_mod.o \
	turbulence_closure_mod.o
parin.o: \
	cpulog_mod.o \
	date_and_time_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	progress_bar_mod.o \
	read_restart_data_mod.o \
	turbulence_closure_mod.o
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
	poisfft_mod.o \
	poismg_mod.o \
	poismg_noopt_mod.o \
	surface_mod.o
prognostic_equations.o: \
	advec_ws.o \
	buoyancy.o \
	coriolis.o \
	cpulog_mod.o \
	cpulog_mod.o \
	diffusion_s.o \
	diffusion_u.o \
	diffusion_v.o \
	diffusion_w.o \
	eqn_state_seawater.o \
	mod_kinds.o \
	modules.o \
	stokes_force_mod.o \
	surface_mod.o \
	turbulence_closure_mod.o
progress_bar_mod.o: \
	mod_kinds.o \
	modules.o
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
	cpulog_mod.o \
	date_and_time_mod.o \
	modules.o \
	mod_kinds.o \
	netcdf_interface_mod.o \
	random_function_mod.o\
	random_generator_parallel_mod.o \
	surface_mod.o
run_control.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
singleton_mod.o: \
	mod_kinds.o
sor.o: \
	mod_kinds.o \
	modules.o
stokes_drift_mod.o: \
	mod_kinds.o \
	modules.o
stokes_force_mod.o: \
	mod_kinds.o \
	modules.o
sum_up_3d_data.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	turbulence_closure_mod.o \
	surface_mod.o
surface_coupler.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	surface_mod.o
surface_layer_fluxes_mod.o: \
	exchange_horiz_2d.o \
	mod_kinds.o \
	modules.o \
	surface_mod.o
surface_mod.o: \
	mod_kinds.o \
	modules.o
swap_timelevel.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	turbulence_closure_mod.o
temperton_fft_mod.o: \
	mod_kinds.o \
	modules.o
time_integration.o: \
	advec_ws.o \
	buoyancy.o \
	calc_mean_profile.o \
	cpulog_mod.o \
	disturb_heatflux.o \
	mod_kinds.o \
	modules.o \
	prognostic_equations.o \
	progress_bar_mod.o \
	stokes_force_mod.o \
	surface_layer_fluxes_mod.o \
	surface_mod.o \
	turbulence_closure_mod.o
time_integration_spinup.o: \
	cpulog_mod.o \
	date_and_time_mod.o \
	mod_kinds.o \
	modules.o \
	surface_layer_fluxes_mod.o
time_to_string.o: \
	mod_kinds.o
timestep.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
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
	advec_ws.o \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	stokes_drift_mod.o \
	stokes_force_mod.o \
	surface_mod.o
wrd_write_string.o: \
	mod_kinds.o
write_restart_data_mod.o: \
	date_and_time_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	random_function_mod.o \
	random_generator_parallel_mod.o
