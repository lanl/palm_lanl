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
             f_mpas,nVertLevels,wtflux,wtflux_solar, wsflux,uwflux, &
             vwflux,fac,dep1,dep2,dzLES,nzLES,        &
             dtDataOutput, dtDisturb, dtDataOutputAv, &
             dtDopr, endTime, dtDots,                &
             tIncrementLES,sIncrementLES,             &
             uIncrementLES,vIncrementLES)


    USE arrays_3d

    USE control_parameters

    USE configure_3D_MODEL

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s, cpu_statistics

    USE indices

    USE netcdf_data_input_mod,                                                 &
        ONLY:  netcdf_data_input_inquire_file, netcdf_data_input_init,         &
               netcdf_data_input_surface_data, netcdf_data_input_topo

    USE kinds

    USE pegrid

    USE random_generator_parallel, ONLY: deallocate_random_generator

    use surface_mod

    use tridia_solver, ONLY: tridia_deallocate

    use turbulence_closure_mod, ONLY: tcm_deallocate_arrays

    USE write_restart_data_mod,                                                &
        ONLY:  wrd_global, wrd_local

    use statistics

    IMPLICIT NONE

!
! -- Variables from MPAS
   integer(iwp) :: nVertLevels, il, jl, kl, knt, nzLES, iz
   Real(wp),dimension(nVertLevels),intent(inout)   :: T_mpas, S_mpas, U_mpas, V_mpas
   Real(wp),dimension(nVertLevels),intent(inout)   :: tIncrementLES, sIncrementLES, &
                                                        uIncrementLES, vIncrementLES
   Real(wp),dimension(nVertLevels),intent(in)      :: lt_mpas
   Real(wp),allocatable,dimension(:)   :: T_mpas2, S_mpas2, U_mpas2, V_mpas2
   Real(wp),allocatable,dimension(:)   :: Tles, Sles, Ules, Vles, zmid, zedge
   real(wp),allocatable,dimension(:)   :: zeLES, wtLES, wsLES, wuLES, wvLES
   Real(wp) :: wtflux, wsflux, uwflux, vwflux, dzLES, z_fac, z_frst, z_cntr
   real(wp) :: z_fac1, z_fac2, z_facn, tol, test, f_mpas, fac, dep1, dep2
   real(wp) :: dtDataOutput, dtDisturb, dtDataOutputAv, dtDopr, endTime, dtDots
   real(wp) :: wtflux_solar
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

   call init_control_parameters

   dt_data_output = dtDataOutput
   dt_disturb = dtDisturb
   dt_data_output_av = dtDataOutputAv
   dt_dopr = dtDopr
   end_time = endTime
   dt_dots = dtDots
   ideal_solar_division = fac
   ideal_solar_efolding1 = dep1
   ideal_solar_efolding2 = dep2
   wb_solar = wtflux_solar
  ! nVertLevels = 50
  ! nzLES = 128
   nz = nzLES
  ! dzLES = 1.0_wp
!   allocate(T_mpas(nVertLevels),S_mpas(nVertLevels),U_mpas(nVertLevels),V_mpas(nVertLevels))
!   allocate(tIncrementLES(nVertLevels),sIncrementLES(nVertLevels))
!   allocate(uIncrementLES(nVertLevels),vIncrementLES(nVertLevels))
   allocate(zmid(nVertLevels),zedge(nVertLevels+1))!,lt_mpas(nVertLevels))
   allocate(T_mpas2(nVertLevels),S_mpas2(nVertLevels),U_mpas2(nVertLevels))
   allocate(V_mpas2(nVertLevels))

!   lt_mpas(:) = 50.0_wp
   zmid(1) = -0.5_wp*lt_mpas(1)
   zedge(1) = 0

   do il=2,nVertLevels
      zmid(il) = zmid(il-1) - 0.5*(lt_mpas(il-1) + lt_mpas(il))
      zedge(il) = zedge(il-1) - lt_mpas(il-1)
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

    print *, top_heatflux
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
    do il = nzt-2,nzb,-1
      zeLES(il) = zeLES(nzt-1)*(z_fac**(float(iz)) - 1.0_wp) / (z_fac - 1.0_wp)
      iz = iz + 1
    enddo
    zeLES(nzb-1) = max(z_cntr,zeLES(nzb) - (zeLES(nzb+1) - zeLES(nzb)))

    do il = nzt,nzb,-1
      zu(il) = 0.5*(zeLES(il) + zeLES(il-1))
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
    il = 1
    do while (il < nVertLevels-1)
      do jl=nzt,nzb,-1
         if(zu(jl) < zmid(il+1)) then
           il = il+1
           il = min(il,nVertLevels-1)
         endif

         coeff2 = (T_mpas(il) - T_mpas(il+1)) / (zmid(il) - zmid(il+1))
         coeff1 = T_mpas(il+1) - coeff2*zmid(il+1)
         pt(jl,:,:) = coeff2*zu(jl) + coeff1 + 273.15_wp

         coeff2 = (S_mpas(il) - S_mpas(il+1)) / (zmid(il) - zmid(il+1))
         coeff1 = S_mpas(il+1) - coeff2*zmid(il+1)
         sa(jl,:,:) = coeff2*zu(jl) + coeff1

         coeff2 = (U_mpas(il) - U_mpas(il+1)) / (zmid(il) - zmid(il+1))
         coeff1 = U_mpas(il+1) - coeff2*zmid(il+1)
         u(jl,:,:) = coeff2*zu(jl) + coeff1

         coeff2 = (V_mpas(il) - V_mpas(il+1)) / (zmid(il) - zmid(il+1))
         coeff1 = V_mpas(il+1) - coeff2*zmid(il+1)
         v(jl,:,:) = coeff2*zu(jl) + coeff1

      enddo
    enddo

    pt(nzt+1,:,:) = pt(nzt,:,:)
    sa(nzt+1,:,:) = sa(nzt,:,:)
    u(nzt+1,:,:) = u(nzt,:,:)
    v(nzt+1,:,:) = v(nzt,:,:)

    !
!-- Output of program header
!    IF ( myid == 0 )  CALL header

    CALL cpu_log( log_point(2), 'initialisation', 'stop' )

!
!-- Set start time in format hh:mm:ss
    simulated_time_chr = time_to_string( time_since_reference_point )

!    IF ( do3d_at_begin )  THEN
!       CALL data_output_3d( 0 )
!    ENDIF

!
!-- Integration of the model equations using timestep-scheme
    CALL time_integration

!-- If required, repeat output of header including the required CPU-time
!    IF ( myid == 0 )  CALL header
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
!    CALL cpu_statistics

    ! need tto interpolate back to mpas for fluxes, include sgs terms?
    Tles = hom(:,1,4,statistic_regions) - 273.15
    Sles = hom(:,1,23,statistic_regions)
    Ules = hom(:,1,1,statistic_regions)
    Vles = hom(:,1,2,statistic_regions)

    il = nzt-1
    do while (il > 1)
      do jl=1,nVertLevels
         if(zmid(jl) < zu(il)) then
           do while (zmid(jl) < zu(il))
              il = il-1
           enddo
           il = max(il,nzb)
         endif

         coeff2 = (Tles(il+1) - Tles(il)) / (zu(il+1) - zu(il))
         coeff1 = Tles(il+1) - coeff2*zu(il+1)

         T_mpas2(jl) = coeff2*zmid(jl) + coeff1 
         tIncrementLES(jl) = (T_mpas2(jl) - T_mpas(jl)) / end_time

         coeff2 = (Sles(il+1) - Sles(il)) / (zu(il+1) - zu(il))
         coeff1 = Sles(il+1) - coeff2*zu(il+1)
         S_mpas2(jl) = coeff2*zmid(jl) + coeff1
         sIncrementLES(jl) = (S_mpas2(jl) - S_mpas(jl)) / end_time

         coeff2 = (Ules(il+1) - Ules(il)) / (zu(il+1) - zu(il))
         coeff1 = Ules(il+1) - coeff2*zu(il+1)
         U_mpas2(jl) = coeff2*zmid(jl) + coeff1
         uIncrementLES(jl) = (U_mpas2(jl) - U_mpas(jl)) / end_time

         coeff2 = (Vles(il+1) - Vles(il)) / (zu(il+1) - zu(il))
         coeff1 = Vles(il+1) - coeff2*zu(il+1)
         V_mpas2(jl) = coeff2*zmid(jl) + coeff1
         vIncrementLES(jl) = (V_mpas2(jl) - V_mpas(jl)) / end_time

      enddo
    enddo

    DEALLOCATE( pt_init, q_init, s_init, ref_state, sa_init, ug,         &
                       u_init, v_init, vg, hom, hom_sum )

   deallocate(hor_index_bounds)

    deallocate(zu,zeLES,Tles,Sles)
    deallocate(hyp, Ules,Vles)
    deallocate(ddzu, ddzw, dd2zu, dzu, dzw, zw, ddzu_pres, nzb_s_inner,  &
               nzb_s_outer, nzb_u_inner, nzb_u_outer, nzb_v_inner,       &
               nzb_v_outer, nzb_w_inner, nzb_w_outer, nzb_diff_s_inner,  &
               nzb_diff_s_outer, wall_flags_0, advc_flags_1, advc_flags_2)

    call deallocate_bc
    call deallocate_3d_variables
    call tcm_deallocate_arrays
    call deallocate_random_generator
    call tridia_deallocate

    close(18)

#if defined( __parallel )
    CALL MPI_FINALIZE( ierr )
#endif

END subroutine palm

subroutine init_control_parameters
    USE arrays_3d

    USE control_parameters

    USE statistics, only: flow_statistics_called
    USE kinds


    openfile = file_status(.FALSE.,.FALSE.)

        poisfft_initialized = .FALSE.
        psolver = 'poisfft'
        momentum_advec = 'ws-scheme'
        loop_optimization = 'vector'
        bc_e_b = 'neumann'
        bc_lr = 'cyclic'
        bc_ns = 'cyclic'
        bc_p_b = 'neumann'
        bc_p_t = 'neumann'
        bc_pt_b = 'neumann'
        bc_pt_t = 'neumann'
        bc_sa_t = 'neumann'
        bc_sa_b = 'neumann'
        bc_uv_b = 'neumann'
        bc_uv_t = 'neumann'
        coupling_mode = 'uncoupled'
        fft_method = 'temperton-algorithm'
        topography = 'flat'
        initializing_actions = 'set_constant_profiles'
        random_generator = 'random-parallel'
        reference_state = 'initial_profile'
        data_output = ' '
        data_output_user = ' '
        doav = ' '
        data_output_masks = ' ' 
        data_output_pr = ' '
        domask = ' '
        do2d = ' '
        do3d = ' '

        do3d_no(0:1) = 0

        abort_mode = 1
        average_count_pr = 0
        average_count_3d = 0
        current_timestep_number = 0
        coupling_topology = 0
        dist_range = 0
        doav_n = 0
        dopr_n = 0
        dopr_time_count = 0
        dopts_time_count = 0
        dots_time_count = 0
        dp_level_ind_b = 0 
        dvrp_filecount = 0
        ensemble_member_nr = 0

        iran = -1234567
        length = 0
        io_group = 0
        io_blocks = 1
        masks = 0
        maximum_parallel_io_streams = -1
        mgcycles = 0
        mg_cycles = 4
        mg_switch_to_pe0_level = -1
        ngsrb = 2
        nr_timesteps_this_run = 0
        nsor = 20
        nsor_ini = 100
        normalizing_region = 0
        num_leg = 0
        num_var_fl_user = 0
        nz_do3 = -9999
        y_shift = 0
        mask_size(max_masks,3) = -1
        mask_size_l(max_masks,3) = -1
        mask_start_l(max_masks,3) = -1
        pt_vertical_gradient_level_ind(10) = -9999
        sa_vertical_gradient_level_ind(10) = -9999
        stokes_drift_method = -9999

        dz(10) = -1.0_wp 
        dzconst = 2.5_wp
        dt_disturb = 20.0_wp 
        dt_do3d = 9999999.9_wp
        dt_3d = 0.01_wp

        simulated_time = 0.0_wp
        flow_statistics_called = .FALSE.
        disturbance_created = .FALSE.
        time_disturb = 0.0_wp
        time_dopr = 0.0_wp
        time_dopr_av = 0.0_wp
        time_dots = 0.0_wp
        time_do2d_xy = 0.0_wp
        time_do2d_xz = 0.0_wp
        time_do2d_yz = 0.0_wp
        time_do3d = 0.0_wp
        time_do_av = 0.0_wp
        time_run_control = 0.0_wp

end subroutine init_control_parameters

subroutine deallocate_memory
        
        use pegrid

        deallocate(hor_index_bounds)

end subroutine deallocate_memory
