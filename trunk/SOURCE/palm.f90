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
 subroutine palm(T_mpas,S_mpas,U_mpas,V_mpas,lt_mpas, &
             f_mpas,nVertLevels,wtflux,wsflux,uwflux, &
             vwflux,fac,dep1,dep2,dzLES,nzLES,        &
             dtDataOutput, dtDisturb, dtDataOutputAv, &
             dtDopr, endTime, dtDots,                &
             tIncrementLES,sIncrementLES,             &
             uIncrementLES,vIncrementLES)


    USE arrays_3d

    USE control_parameters,                                                     &
        ONLY:  data_output, data_output_pr,constant_diffusion, do3d_at_begin,   &
               section_xy, initializing_actions, io_blocks, io_group,           &
               message_string, runnr, simulated_time, simulated_time_chr,       &
               time_since_reference_point, write_binary, top_heatflux,          &
               top_momentumflux_u, top_momentumflux_v, top_salinityflux, f,     &
               end_time, dt_dopr, dt_data_output, dt_data_output_av, dt_disturb, &
               dt_dots, ideal_solar_division, ideal_solar_efolding1,            &
               ideal_solar_efolding2

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

    use statistics, ONLY: hom, statistic_regions

    IMPLICIT NONE

!
! -- Variables from MPAS
   integer(iwp) :: nVertLevels, i, j, k, knt, nzLES, iz
   Real(wp),allocatable,dimension(:),intent(inout)   :: T_mpas, S_mpas, U_mpas, V_mpas
   Real(wp),allocatable,dimension(:),intent(inout)   :: tIncrementLES, sIncrementLES, &
                                                        uIncrementLES, vIncrementLES
   Real(wp),allocatable,dimension(:)   :: lt_mpas, T_mpas2, S_mpas2, U_mpas2, V_mpas2
   Real(wp),allocatable,dimension(:)   :: Tles, Sles, Ules, Vles, zmid, zedge
   real(wp),allocatable,dimension(:)   :: zeLES, wtLES, wsLES, wuLES, wvLES
   Real(wp) :: wtflux, wsflux, uwflux, vwflux, dzLES, z_fac, z_frst, z_cntr
   real(wp) :: z_fac1, z_fac2, z_facn, tol, test, f_mpas, fac, dep1, dep2
   real(wp) :: dtDataOutput, dtDisturb, dtDataOutputAv, dtDopr, endTime, dtDots
!
!-- Local variables
   CHARACTER(LEN=9)  ::  time_to_string  !<
   CHARACTER(LEN=10) ::  env_string      !< to store string of environment var
   INTEGER(iwp)      ::  env_stat        !< to hold status of GET_ENV
   INTEGER(iwp)      ::  myid_openmpi    !< OpenMPI local rank for CUDA aware MPI
   Real(wp) :: coeff1, coeff2

!more arguments to send
! dt_data_output, dt_disturb, dt_data_output_av, dt_dopr
! end_time

   dt_data_output = dtDataOutput
   dt_disturb = dtDisturb
   dt_data_output_av = dtDataOutputAv
   dt_dopr = dtDopr
   end_time = endTime
   dt_dots = dtDots
   ideal_solar_division = fac
   ideal_solar_efolding1 = dep1
   ideal_solar_efolding2 = dep2

  ! nVertLevels = 50
  ! nzLES = 128
   nz = nzLES
  ! dzLES = 1.0_wp
   allocate(T_mpas(nVertLevels),S_mpas(nVertLevels),U_mpas(nVertLevels),V_mpas(nVertLevels))
   allocate(tIncrementLES(nVertLevels),sIncrementLES(nVertLevels))
   allocate(uIncrementLES(nVertLevels),vIncrementLES(nVertLevels))
   allocate(zmid(nVertLevels),zedge(nVertLevels+1),lt_mpas(nVertLevels))

!   lt_mpas(:) = 50.0_wp
   zmid(1) = -0.5_wp*lt_mpas(1)
   zedge(1) = 0

   do i=2,nVertLevels
      zmid(i) = zmid(i-1) - 0.5*(lt_mpas(i-1) + lt_mpas(i))
      zedge(i) = zedge(i-1) - lt_mpas(i-1)
   enddo

   zedge(nvertLevels+1) = zedge(nVertLevels) - lt_mpas(nVertLevels)

!   U_mpas(:) = 0.0_wp
!   V_mpas(:) = 0.0_wp
!   S_mpas(:) = 34.0_wp!
!   do i=1,nVertLevels
!      T_mpas(i) = 293.15 + 0.005*zmid(i)
  ! enddo

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

!   f_mpas = 1e-4
!   uwflux = 0.0
!   vwflux = 0.0
!   wsflux = 1e-4
!   wtflux = -1.78e-5

!TODO add check for right / acceptable range.
    top_momentumflux_u = uwflux
    top_momentumflux_v = vwflux
    top_heatflux = wtflux
    top_salinityflux = wsflux
    f = f_mpas

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

    allocate(zu(nzb:nzt+1),zeLES(nzb-1:nzt+1),Tles(0:nzLES+1),Sles(0:nzLES+1))
    allocate(Ules(0:nzLES+1),Vles(0:nzLES+1))

    nzt = nzLES
    ! construct a stretched stretched grid
    z_cntr = zedge(nVertLevels+1)
    z_frst = -dzLES
    z_fac1 = z_cntr / z_frst
    z_fac2 = 1.0_wp / float(nzt)
    z_fac = 1.10_wp
    tol = 1.0E-10_wp
    test = 10.00_wp
    knt = 0

    do while (test > tol)
      knt = knt + 1
      z_facn = (z_fac1*(z_fac - 1.0_wp) + 1.0_wp)**z_fac2
      test = abs(1.0 - z_facn / z_fac)
      if(knt .gt. 500) THEN
        print *, 'cannot find stretching factor,'
        print *, 'z_fac = ',z_fac, 'z_facn = ',z_facn, 'knt = ',knt
        stop
      ENDIF
      z_fac = z_facn
    enddo

    zeLES(nzt+1) = dzLES
    zeLES(nzt) = 0.0_wp
    zeLES(nzt-1) = -dzLES
    iz = 2
    do i = nzt-2,nzb,-1
      zeLES(i) = zeLES(nzt-1)*(z_fac**(float(iz)) - 1.0_wp) / (z_fac - 1.0_wp)
      iz = iz + 1
    enddo
    zeLES(nzb-1) = max(z_cntr,zeLES(nzb) - (zeLES(nzb+1) - zeLES(nzb)))

    do i = nzt,nzb,-1
      zu(i) = 0.5*(zeLES(i) + zeLES(i-1))
    enddo
    zu(nzt+1) = dzLES*0.5

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
         if(zu(j) < zmid(i+1)) then
           i = i+1
           i = min(i,nVertLevels-1)
         endif

         coeff2 = (T_mpas(i) - T_mpas(i+1)) / (zmid(i) - zmid(i+1))
         coeff1 = T_mpas(i+1) - coeff2*zmid(i+1)
         pt(j,:,:) = coeff2*zu(j) + coeff1

         coeff2 = (S_mpas(i) - S_mpas(i+1)) / (zmid(i) - zmid(i+1))
         coeff1 = S_mpas(i+1) - coeff2*zmid(i+1)
         sa(j,:,:) = coeff2*zu(j) + coeff1

         coeff2 = (U_mpas(i) - U_mpas(i+1)) / (zmid(i) - zmid(i+1))
         coeff1 = U_mpas(i+1) - coeff2*zmid(i+1)
         u(j,:,:) = coeff2*zu(j) + coeff1

         coeff2 = (V_mpas(i) - V_mpas(i+1)) / (zmid(i) - zmid(i+1))
         coeff1 = V_mpas(i+1) - coeff2*zmid(i+1)
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

    ! need tto interpolate back to mpas for fluxes, include sgs terms?
    Tles = hom(:,1,4,statistic_regions) - 273.15
    Sles = hom(:,1,23,statistic_regions)
    Ules = hom(:,1,1,statistic_regions)
    Vles = hom(:,1,2,statistic_regions)

    i = nzt-1
    do while (i > 1)
      do j=1,nVertLevels
         if(zmid(j) < zu(i)) then
           do while (zmid(j) < zu(i))
              i = i-1
           enddo
           i = max(i,nzb)
         endif

         coeff2 = (Tles(i+1) - Tles(i)) / (zu(i+1) - zu(i))
         coeff1 = Tles(i+1) - coeff2*zu(i+1)
         T_mpas2(j) = coeff2*zmid(j) + coeff1
         tIncrementLES(j) = (T_mpas2(j) - T_mpas(j)) / end_time

         coeff2 = (Sles(i+1) - Sles(i)) / (zu(i+1) - zu(i))
         coeff1 = Sles(i+1) - coeff2*zu(i+1)
         S_mpas2(j) = coeff2*zmid(j) + coeff1
         sIncrementLES(j) = (S_mpas2(j) - S_mpas(j)) / end_time

         coeff2 = (Ules(i+1) - Ules(i)) / (zu(i+1) - zu(i))
         coeff1 = Ules(i+1) - coeff2*zu(i+1)
         U_mpas2(j) = coeff2*zmid(j) + coeff1
         uIncrementLES(j) = (U_mpas2(j) - U_mpas(j)) / end_time

         coeff2 = (Vles(i+1) - Vles(i)) / (zu(i+1) - zu(i))
         coeff1 = Vles(i+1) - coeff2*zu(i+1)
         V_mpas2(j) = coeff2*zmid(j) + coeff1
         vIncrementLES(j) = (V_mpas2(j) - V_mpas(j)) / end_time

      enddo
    enddo

#if defined( __parallel )
    CALL MPI_FINALIZE( ierr )
#endif

END subroutine palm
