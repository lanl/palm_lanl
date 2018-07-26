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
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: palm.f90 2977 2018-04-17 10:27:57Z kanani $
! Deduct spinup_time from RUN_CONTROL output of main 3d run 
! (use time_since_reference_point instead of simulated_time)
! 
! 2951 2018-04-06 09:05:08Z kanani
! Add log_point_s for pmci_init
! 
! 2903 2018-03-16 08:17:06Z hellstea
! Nesting-related calls to pmci_ensure_nest_mass_conservation and pres after
! the nest initialization are removed as they may create unwanted initial
! perturbation in some cases.
! 
! 2894 2018-03-15 09:17:58Z Giersch
! Modified todo list, _write_restart_data/_last_actions has been renamed to 
! _wrd_local, unit 14 will be opened now for each io_group
! write_3d_binary is called wrd_local now, wrd_global moved from wrd_local to 
! palm.f90, unit 14 is closed directly after the wrd_local call, Module related
! routines for writing restart data have been moved to wrd_local
! 
! 2801 2018-02-14 16:01:55Z suehring
! Changed lpm from subroutine to module.
! Introduce particle transfer in nested models.
! 
! 2766 2018-01-22 17:17:47Z kanani
! Removed preprocessor directive __chem
! 
! 2720 2018-01-02 16:27:15Z kanani
! Version update to 5.0
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Implementation of chemistry module (FK)
! Introduce input-data standard
! Rename lsm_last_actions into lsm_write_restart_data
! Move usm_write_restart_data into io_blocks loop (MS)
! 
! 2512 2017-10-04 08:26:59Z raasch
! user interface required revision updated
! 
! 2320 2017-07-21 12:47:43Z suehring
! Modularize large-scale forcing and nudging
! 
! 2298 2017-06-29 09:28:18Z raasch
! type of write_binary changed from CHARACTER to LOGICAL,
! user interface required revision updated, MPI2 related part removed
! 
! 2296 2017-06-28 07:53:56Z maronga
! Added call to new spinup routine
! 
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison' 
! includes two more prognostic equations for cloud drop concentration (nc)  
! and cloud water content (qc). 
! 
! 2261 2017-06-08 14:25:57Z raasch
! output of run number for mrun to create unified cycle numbers
! 
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Renamed wall_flags_0 and wall_flags_00 into advc_flags_1 and advc_flags_2, 
! respectively, within copyin statement. Moreover, introduced further flag 
! array wall_flags_0. 
! Remove unused variables from ONLY list.
! 
! 2178 2017-03-17 11:07:39Z hellstea
! Calls for pmci_ensure_nest_mass_conservation and pres are added after
! the nest initialization
!
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC directives and related code removed
! 
! 2011 2016-09-19 17:29:57Z kanani
! Flag urban_surface is now defined in module control_parameters.
! 
! 2007 2016-08-24 15:47:17Z kanani
! Temporarily added CALL for writing of restart data for urban surface model
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1976 2016-07-27 13:28:04Z maronga
! Added call to radiation_last_actions for binary output of land surface model 
! data
! 
! 1972 2016-07-26 07:52:02Z maronga
! Added call to lsm_last_actions for binary output of land surface model data
! 
! 1960 2016-07-12 16:34:24Z suehring
! Separate humidity and passive scalar
! 
! 1834 2016-04-07 14:34:20Z raasch
! Initial version of purely vertical nesting introduced. 
! 
! 1833 2016-04-07 14:23:03Z raasch
! required user interface version changed
!
! 1808 2016-04-05 19:44:00Z raasch
! routine local_flush replaced by FORTRAN statement
!
! 1783 2016-03-06 18:36:17Z raasch
! required user interface version changed
!
! 1781 2016-03-03 15:12:23Z raasch
! pmc initialization moved from time_integration to here
!
! 1779 2016-03-03 08:01:28Z raasch
! setting of nest_domain and coupling_char moved to the pmci
!
! 1764 2016-02-28 12:45:19Z raasch
! cpp-statements for nesting removed, communicator settings cleaned up
!
! 1762 2016-02-25 12:31:13Z hellstea
! Introduction of nested domain feature
!
! 1747 2016-02-08 12:25:53Z raasch
! OpenACC-adjustment for new surface layer parameterization
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
!
! 1668 2015-09-23 13:45:36Z raasch
! warning replaced by abort in case of failed user interface check
!
! 1666 2015-09-23 07:31:10Z raasch
! check for user's interface version added
!
! 1482 2014-10-18 12:34:45Z raasch
! adjustments for using CUDA-aware OpenMPI
!
! 1468 2014-09-24 14:06:57Z maronga
! Adapted for use on up to 6-digit processor cores
! 
! 1402 2014-05-09 14:25:13Z raasch
! location messages added
! 
! 1374 2014-04-25 12:55:07Z raasch
! bugfix: various modules added
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
! module interfaces removed
!
! 1241 2013-10-30 11:36:58Z heinze
! initialization of nuding and large scale forcing from external file
!
! 1221 2013-09-10 08:59:13Z raasch
! +wall_flags_00, rflags_invers, rflags_s_inner in copyin statement
!
! 1212 2013-08-15 08:46:27Z raasch
! +tri in copyin statement
!
! 1179 2013-06-14 05:57:58Z raasch
! ref_state added to copyin-list
!
! 1113 2013-03-10 02:48:14Z raasch
! openACC statements modified
!
! 1111 2013-03-08 23:54:10Z raasch
! openACC statements updated
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1015 2012-09-27 09:23:24Z raasch
! Version number changed from 3.8 to 3.8a.
! OpenACC statements added + code changes required for GPU optimization
!
! 849 2012-03-15 10:35:09Z raasch
! write_particles renamed lpm_write_restart_file
!
! Revision 1.1  1997/07/24 11:23:35  raasch
! Initial revision
!
!
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

    USE chemistry_model_mod,                                                   &
        ONLY:  chem_init

    USE chem_photolysis_mod,                                                   &
        ONLY:  photolysis_init

    USE control_parameters,                                                    &
        ONLY:  air_chemistry,                                                  &
               cloud_physics, constant_diffusion, coupling_char, coupling_mode,&
               do2d_at_begin, do3d_at_begin, humidity, initializing_actions,   &
               io_blocks, io_group, large_scale_forcing,                       &
               message_string, microphysics_morrison, microphysics_seifert,    &
               nest_domain, neutral, nudging, passive_scalar, runnr,           &
               simulated_time, simulated_time_chr, spinup,                     &
               time_since_reference_point,                                     &
               user_interface_current_revision,                                &
               user_interface_required_revision, version, wall_heatflux,       &
               write_binary

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s, cpu_statistics

    USE indices,                                                               &
        ONLY:  nbgp

    USE netcdf_data_input_mod,                                                 &
        ONLY:  netcdf_data_input_inquire_file, netcdf_data_input_init,         &
               netcdf_data_input_surface_data, netcdf_data_input_topo

    USE kinds

    USE particle_attributes,                                                   &
        ONLY:  particle_advection

    USE pegrid

    USE pmc_particle_interface,                                                &
        ONLY: pmcp_g_alloc_win

    USE pmc_interface,                                                         &
        ONLY:  cpl_id, nested_run, pmci_child_initialize, pmci_init,           &
               pmci_modelconfiguration, pmci_parent_initialize,                &
               pmci_ensure_nest_mass_conservation

    USE write_restart_data_mod,                                                &
        ONLY:  wrd_global, wrd_local


    IMPLICIT NONE

!
!-- Local variables
    CHARACTER(LEN=9)  ::  time_to_string  !<
    CHARACTER(LEN=10) ::  env_string      !< to store string of environment var
    INTEGER(iwp)      ::  env_stat        !< to hold status of GET_ENV
    INTEGER(iwp)      ::  i               !<
    INTEGER(iwp)      ::  myid_openmpi    !< OpenMPI local rank for CUDA aware MPI

    version = 'PALM 5.0'
    user_interface_required_revision = 'r2512'

#if defined( __parallel )
!
!-- MPI initialisation. comm2d is preliminary set, because
!-- it will be defined in init_pegrid but is used before in cpu_log.
    CALL MPI_INIT( ierr )

!
!-- Initialize the coupling for nested-domain runs
!-- comm_palm is the communicator which includes all PEs (MPI processes)
!-- available for this (nested) model. If it is not a nested run, comm_palm
!-- is returned as MPI_COMM_WORLD
    CALL cpu_log( log_point_s(70), 'pmci_init', 'start' )
    CALL pmci_init( comm_palm )
    CALL cpu_log( log_point_s(70), 'pmci_init', 'stop' )
    comm2d = comm_palm
!
!-- Get the (preliminary) number of MPI processes and the local PE-id (in case
!-- of a further communicator splitting in init_coupling, these numbers will
!-- be changed in init_pegrid).
    IF ( nested_run )  THEN

       CALL MPI_COMM_SIZE( comm_palm, numprocs, ierr )
       CALL MPI_COMM_RANK( comm_palm, myid, ierr )

    ELSE

       CALL MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
       CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
!
!--    Initialize PE topology in case of coupled atmosphere-ocean runs (comm_palm
!--    will be splitted in init_coupling)
       CALL init_coupling
    ENDIF
#endif

!
!-- Initialize measuring of the CPU-time remaining to the run
    CALL local_tremain_ini

!
!-- Start of total CPU time measuring.
    CALL cpu_log( log_point(1), 'total', 'start' )
    CALL cpu_log( log_point(2), 'initialisation', 'start' )

!
!-- Open a file for debug output
    WRITE (myid_char,'(''_'',I6.6)')  myid
    OPEN( 9, FILE='DEBUG'//TRIM( coupling_char )//myid_char, FORM='FORMATTED' )

!
!-- Initialize dvrp logging. Also, one PE maybe split from the global
!-- communicator for doing the dvrp output. In that case, the number of
!-- PEs available for PALM is reduced by one and communicator comm_palm
!-- is changed respectively.
#if defined( __parallel )
    CALL MPI_COMM_RANK( comm_palm, myid, ierr )
!
!-- TEST OUTPUT (TO BE REMOVED)
    WRITE(9,*) '*** coupling_mode = "', TRIM( coupling_mode ), '"'
    FLUSH( 9 )
    IF ( TRIM( coupling_mode ) /= 'uncoupled' )  THEN
       PRINT*, '*** PE', myid, ' Global target PE:', target_id, &
               TRIM( coupling_mode )
    ENDIF
#endif

    CALL init_dvrp_logging

!
!-- Read control parameters from NAMELIST files and read environment-variables
    CALL parin

!
!-- Check for the user's interface version
    IF ( user_interface_current_revision /= user_interface_required_revision )  &
    THEN
       message_string = 'current user-interface revision "' //                  &
                        TRIM( user_interface_current_revision ) // '" does ' // &
                        'not match the required revision ' //                   &
                        TRIM( user_interface_required_revision )
        CALL message( 'palm', 'PA0169', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Determine processor topology and local array indices
    CALL init_pegrid
!
!-- Check if input file according to input-data standard exists
    CALL netcdf_data_input_inquire_file
!
!-- Read topography input data if required. This is required before the 
!-- numerical grid is finally created in init_grid
    CALL netcdf_data_input_topo  
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
!-- Initialize chemistry (called before check_parameters due to dependencies)
!-- --> Needs to be moved!! What is the dependency about?
! IF (  TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
    IF ( air_chemistry )  THEN
       CALL chem_init
       CALL photolysis_init   ! probably also required for restart
    ENDIF
! END IF
!
!-- Check control parameters and deduce further quantities
    CALL check_parameters

!
!-- Initialize all necessary variables
    CALL init_3d_model

!
!-- Coupling protocol setup for nested-domain runs
    IF ( nested_run )  THEN
       CALL pmci_modelconfiguration
!
!--    Receive and interpolate initial data on children.
!--    Child initialization must be made first if the model is both child and
!--    parent if necessary
       IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
          CALL pmci_child_initialize
!
!--       Send initial condition data from parent to children
          CALL pmci_parent_initialize
!
!--    Exchange_horiz is needed after the nest initialization
          IF ( nest_domain )  THEN
             CALL exchange_horiz( u, nbgp )
             CALL exchange_horiz( v, nbgp )
             CALL exchange_horiz( w, nbgp )
             IF ( .NOT. neutral )  THEN
                CALL exchange_horiz( pt, nbgp )
             ENDIF
             IF ( .NOT. constant_diffusion )  CALL exchange_horiz( e, nbgp )
             IF ( humidity )  THEN
                CALL exchange_horiz( q, nbgp )
                IF ( cloud_physics  .AND.  microphysics_morrison )  THEN
                  CALL exchange_horiz( qc, nbgp )
                  CALL exchange_horiz( nc, nbgp )
                ENDIF
                IF ( cloud_physics  .AND.  microphysics_seifert )  THEN
                   CALL exchange_horiz( qr, nbgp ) 
                   CALL exchange_horiz( nr, nbgp )
                ENDIF
             ENDIF
             IF ( passive_scalar )  CALL exchange_horiz( s, nbgp )
          ENDIF
       ENDIF

       CALL pmcp_g_alloc_win                    ! Must be called after pmci_child_initialize and pmci_parent_initialize
    ENDIF

!
!-- Output of program header
    IF ( myid == 0 )  CALL header

    CALL cpu_log( log_point(2), 'initialisation', 'stop' )

!
!-- Integration of the non-atmospheric equations (land surface model, urban
!-- surface model)
    IF ( spinup )  THEN
       CALL time_integration_spinup
    ENDIF

!
!-- Set start time in format hh:mm:ss
    simulated_time_chr = time_to_string( time_since_reference_point )

!
!-- If required, output of initial arrays
    IF ( do2d_at_begin )  THEN
       CALL data_output_2d( 'xy', 0 )
       CALL data_output_2d( 'xz', 0 )
       CALL data_output_2d( 'yz', 0 )
    ENDIF

    IF ( do3d_at_begin )  THEN
       CALL data_output_3d( 0 )
    ENDIF

!
!-- Integration of the model equations using timestep-scheme
    CALL time_integration

!
!-- If required, write binary data for restart runs
    IF ( write_binary )  THEN

       CALL cpu_log( log_point(22), 'wrd_local', 'start' )

       CALL location_message( 'writing restart data', .FALSE. )

       DO  i = 0, io_blocks-1
          IF ( i == io_group )  THEN

!
!--          Open binary file 
             CALL check_open( 14 )
!
!--          Write control parameters and other global variables for restart.
             IF ( myid == 0 )  CALL wrd_global
!
!--          Write processor specific flow field data for restart runs
             CALL wrd_local
!
!--          Close binary file 
             CALL close_file( 14 )

          ENDIF
#if defined( __parallel )
          CALL MPI_BARRIER( comm2d, ierr )
#endif
       ENDDO

       CALL location_message( 'finished', .TRUE. )

       CALL cpu_log( log_point(22), 'wrd_local', 'stop' )

!
!--    If required, write particle data in own restart files
       IF ( particle_advection )  CALL lpm_write_restart_file
       
    ENDIF

!
!-- If required, repeat output of header including the required CPU-time
    IF ( myid == 0 )  CALL header
!
!-- If required, final  user-defined actions, and
!-- last actions on the open files and close files. Unit 14 was opened
!-- in wrd_local but it is closed here, to allow writing on this
!-- unit in routine user_last_actions.
    CALL cpu_log( log_point(4), 'last actions', 'start' )
          
    CALL user_last_actions
    CALL close_file( 0 )
    CALL close_dvrp

    CALL cpu_log( log_point(4), 'last actions', 'stop' )

!
!-- Write run number to file (used by mrun to create unified cycle numbers for
!-- output files
    IF ( myid == 0  .AND.  runnr > 0 )  THEN
       OPEN( 90, FILE='RUN_NUMBER', FORM='FORMATTED' )
       WRITE( 90, '(I4)' )  runnr
       CLOSE( 90 )
    ENDIF

!
!-- Take final CPU-time for CPU-time analysis
    CALL cpu_log( log_point(1), 'total', 'stop' )
    CALL cpu_statistics

#if defined( __parallel )
    CALL MPI_FINALIZE( ierr )
#endif

 END PROGRAM palm
