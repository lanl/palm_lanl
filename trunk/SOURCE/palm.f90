!> @file palm.f90
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
! Description:
! ------------
!> Large-Eddy Simulation (LES) model for the convective boundary layer,
!> optimized for use on parallel machines (implementation realized using the
!> Message Passing Interface (MPI)). The model can also be run on vector machines
!> (less well optimized) and workstations. Versions for the different types of
!> machines are controlled via cpp-directives.
!> Model runs are only feasible using the ksh-script mrun.
!>
!> @todo create routine last_actions instead of calling lsm_last_actions etc.
!> @todo move chem_init call to init_3d_model or to check_parameters
!------------------------------------------------------------------------------!
 PROGRAM palm

    USE arrays_3d

    USE control_parameters

    USE configure_3D_MODEL

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s, cpu_statistics

    USE indices

    USE netcdf_data_input_mod,                                                 &
        ONLY:  netcdf_data_input_inquire_file, netcdf_data_input_init,         &
               netcdf_data_input_surface_data, netcdf_data_input_topo

    USE fft_xy,                                                                &
        ONLY: fft_finalize

    USE kinds

    USE pegrid

    USE write_restart_data_mod,                                                &
        ONLY:  wrd_global, wrd_local

#if defined( __cudaProfiler )
    USE cudafor
#endif

    IMPLICIT NONE

!
!-- Local variables
   CHARACTER(LEN=9)  ::  time_to_string  !<
   CHARACTER(LEN=10) ::  env_string      !< to store string of environment var
   INTEGER(iwp)      ::  env_stat        !< to hold status of GET_ENV
   INTEGER(iwp)      ::  myid_openmpi    !< OpenMPI local rank for CUDA aware MPI

   ! will need to interpolate profiles to an LES grid, or maybe assume the same???
   ! at end assign fluxes tot mpas variables and end routine.

#if defined( __parallel )
!
!-- MPI initialisation. comm2d is preliminary set, because
!-- it will be defined in init_pegrid but is used before in cpu_log.
    CALL MPI_INIT( ierr )

    comm_palm = MPI_COMM_WORLD
    comm2d = comm_palm
!
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
    !
#endif

!
!-- Initialize measuring of the CPU-time remaining to the run
    CALL local_tremain_ini
!
!-- Start of total CPU time measuring.
    CALL cpu_log( log_point(1), 'total', 'start' )
    CALL cpu_log( log_point(2), 'initialisation', 'start' )
!
!-- Read control parameters from NAMELIST files and read environment-variables
    CALL parin
!
!-- Determine processor topology and local array indices
    CALL init_pegrid
!
!-- Check if input file according to input-data standard exists
    CALL netcdf_data_input_inquire_file
!
!-- Generate grid parameters, initialize generic topography and further process
!-- topography information if required
    CALL init_grid
!
!-- Read global attributes if available.
    CALL netcdf_data_input_init
!
!-- Read surface classification data, e.g. vegetation and soil types, water
!-- surfaces, etc., if available. Some of these data is required before
!-- check parameters is invoked.
    CALL netcdf_data_input_surface_data
!
!-- Check control parameters and deduce further quantities
    CALL check_parameters
!
!-- Initialize all necessary variables
    CALL init_3d_model
!
!-- Output of program header
    IF ( myid == 0 )  CALL header

    CALL cpu_log( log_point(2), 'initialisation', 'stop' )
!
!-- Set start time in format hh:mm:ss
    simulated_time_chr = time_to_string( time_since_reference_point )

    IF ( do3d_at_begin )  THEN
       CALL data_output_3d( 0 )
    ENDIF

#if defined( __cudaProfiler )
!-- Only profile time_integration
    CALL cudaProfilerStart()
#endif
!
!-- Integration of the model equations using timestep-scheme
    CALL time_integration

#if defined( __cudaProfiler )
!-- Only profile time_integration
    CALL cudaProfilerStop()
#endif

!
!-- If required, repeat output of header including the required CPU-time
    IF ( myid == 0 )  CALL header
!
!-- If required, final  user-defined actions, and
!-- last actions on the open files and close files. Unit 14 was opened
!-- in wrd_local but it is closed here, to allow writing on this
!-- unit in routine user_last_actions.
    CALL cpu_log( log_point(4), 'last actions', 'start' )

    CALL close_file( 0 )

    CALL cpu_log( log_point(4), 'last actions', 'stop' )

!
!-- Take final CPU-time for CPU-time analysis
    CALL cpu_log( log_point(1), 'total', 'stop' )
    CALL cpu_statistics

    CALL fft_finalize
#if defined( __parallel )
    CALL MPI_FINALIZE( ierr )
#endif

END PROGRAM palm
