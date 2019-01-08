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
 SUBROUTINE palm(T_mpas,S_mpas,U_mpas,V_mpas,wt,ws,uw,vw,zmid_mpas,zedge_mpas, &
             wtflux,wsflux,uwflux,vwflux,zuLES)

    USE arrays_3d

    USE control_parameters,                                                    &
        ONLY:  data_output, data_output_pr,constant_diffusion, do3d_at_begin,               &
               section_xy, initializing_actions, io_blocks, io_group,                      &
               message_string, runnr, simulated_time, simulated_time_chr,      &
               time_since_reference_point, write_binary, top_heatflux,     &
               top_momentumflux_u, top_momentumflux_v, top_salinityflux

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s, cpu_statistics

    USE indices,                                                               &
        ONLY:  nbgp, nz, nzt, nzb

    USE netcdf_data_input_mod,                                                 &
        ONLY:  netcdf_data_input_inquire_file, netcdf_data_input_init,         &
               netcdf_data_input_surface_data, netcdf_data_input_topo

    USE kinds

    USE pegrid

    USE write_restart_data_mod,                                                &
        ONLY:  wrd_global, wrd_local

    IMPLICIT NONE

!
! -- Variables from MPAS
   integer(iwp) :: nVertLevels, i, j, k
   Real(wp),allocatable,dimension(:)   :: T_mpas, S_mpas, U_mpas, V_mpas
   Real(wp),allocatable,dimension(:)   :: wt, ws, uw, vw, zmid_mpas, zedge_mpas
   Real(wp) :: wtflux, wsflux, uwflux, vwflux, zuLES
!
!-- Local variables
    CHARACTER(LEN=9)  ::  time_to_string  !<
    CHARACTER(LEN=10) ::  env_string      !< to store string of environment var
    INTEGER(iwp)      ::  env_stat        !< to hold status of GET_ENV
    INTEGER(iwp)      ::  myid_openmpi    !< OpenMPI local rank for CUDA aware MPI
   Real(wp) :: coeff1, coeff2

   nVertLevels = 20
   allocate(T_mpas(nVertLevels),S_mpas(nVertLevels),U_mpas(nVertLevels),V_mpas(nVertLevels))
   allocate(zmid_mpas(nVertLevels),zedge_mpas(nVertLevels+1))
   allocate(wt(nVertLevels),ws(nVertLevels),uw(nVertLevels),vw(nVertLevels))

   wtflux = 1.78e-5
   wsflux = 1e-4 
   uwflux = 0.0
   vwflux = 0.0


   zedge_mpas(1) = 0.0_wp
   zmid_mpas(1) = -5.0_wp

   do i=2,nVertLevels
      zmid_mpas(i) = zmid_mpas(i-1) - 10.0_wp
      zedge_mpas(i) = zedge_mpas(i-1) - 10.0_wp
   enddo


   zedge_mpas(nVertLevels+1) = zedge_mpas(nVertLevels) - 10.0_wp

   U_mpas(:) = 0.0_wp
   V_mpas(:) = 0.0_wp
   S_mpas(:) = 34.0_wp
   do i=1,nVertLevels
      T_mpas(i) = 293.15 + 0.005*zmid_mpas(i)
   enddo

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

!TODO add check for right / acceptable range. 
    top_momentumflux_u = uwflux
    top_momentumflux_v = vwflux
    top_heatflux = wtflux
    top_salinityflux = wsflux 
   
    !TODO ooverride the LES setting from a namelist
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
!-- Set netcdf output fields -- not sure how or if to do this once hooked to MPAS
!
    section_xy(1) = 8

    data_output(1) = 'shf*_xy'
    data_output(2) = 'e'
    data_output(3) = 'pt'
    data_output(4) = 'sa'
    data_output(5) = 'u'
    data_output(6) = 'v'
    data_output(7) = 'w'
    data_output(8) = 'rho_ocean'
    data_output(9) = 'alpha_T'

    data_output_pr(1) = 'e'
    data_output_pr(2) = 'e*'
    data_output_pr(3) = '#pt'
    data_output_pr(4) = '#sa'
    data_output_pr(5) = 'p'
    data_output_pr(6) = 'hyp'
    data_output_pr(7) = 'km'
    data_output_pr(8) = 'kh'
    data_output_pr(9) = 'l'
    data_output_pr(10) = '#u'
    data_output_pr(11) = '#v'
    data_output_pr(12) = 'w'
    data_output_pr(13) = 'prho'
    data_output_pr(14) = 'w"u"'
    data_output_pr(15) = 'w*u*'
    data_output_pr(16) = 'w"v"'
    data_output_pr(17) = 'w*v*'
    data_output_pr(18) = 'w"pt"'
    data_output_pr(19) = 'w*pt*'
    data_output_pr(20) = 'w"sa"'
    data_output_pr(21) = 'w*sa*'
    data_output_pr(22) = 'w*e*'
    data_output_pr(23) = 'u*2'
    data_output_pr(24) = 'v*2'
    data_output_pr(25) = 'w*2'
    data_output_pr(26) = 'pt*2'
    data_output_pr(27) = 'w*3'
    data_output_pr(28) = 'Sw'
    data_output_pr(29) = 'w*2pt*'
    data_output_pr(30) = 'w*pt*2'
    data_output_pr(31) = 'w*u*u*:dz'
    data_output_pr(32) = 'w*p*:dz'
    data_output_pr(33) = 'rho_ocean'
    data_output_pr(34) = 'alpha_T'

!-- Determine processor topology and local array indices
    CALL init_pegrid

    allocate(zu(nzb:nzt+1))

    zu(nzt+1) = 0.5
    do i = nzt,nzb,-1
      zu(i) = zu(i+1) - 1.0
    enddo

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

!-- Check control parameters and deduce further quantities
    CALL check_parameters
! interpolate mpas data to les and send to init variable
!
!-- Initialize all necessary variables
    CALL init_3d_model
    i = 1
    do while (i < nVertLevels-1)
      do j=nzt,nzb,-1
         if(zu(j) < zmid_mpas(i+1)) then
           i = i+1
           i = min(i,nVertLevels-1) 
         endif

         coeff2 = (T_mpas(i) - T_mpas(i+1)) / (zmid_mpas(i) - zmid_mpas(i+1))
         coeff1 = T_mpas(i+1) - coeff2*zmid_mpas(i+1)
         pt(j,:,:) = coeff2*zu(j) + coeff1

         coeff2 = (S_mpas(i) - S_mpas(i+1)) / (zmid_mpas(i) - zmid_mpas(i+1))
         coeff1 = S_mpas(i+1) - coeff2*zmid_mpas(i+1)
         sa(j,:,:) = coeff2*zu(j) + coeff1

         coeff2 = (U_mpas(i) - U_mpas(i+1)) / (zmid_mpas(i) - zmid_mpas(i+1))
         coeff1 = U_mpas(i+1) - coeff2*zmid_mpas(i+1)
         u(j,:,:) = coeff2*zu(j) + coeff1

         coeff2 = (V_mpas(i) - V_mpas(i+1)) / (zmid_mpas(i) - zmid_mpas(i+1))
         coeff1 = V_mpas(i+1) - coeff2*zmid_mpas(i+1)
         v(j,:,:) = coeff2*zu(j) + coeff1

      enddo
    enddo

    pt(nzt+1,:,:) = pt(nzt,:,:)
    sa(nzt+1,:,:) = sa(nzt,:,:)
    u(nzt+1,:,:) = u(nzt,:,:)
    v(nzt+1,:,:) = v(nzt,:,:)

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

!
!-- Integration of the model equations using timestep-scheme
    CALL time_integration

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

#if defined( __parallel )
    CALL MPI_FINALIZE( ierr )
#endif

 END SUBROUTINE palm
