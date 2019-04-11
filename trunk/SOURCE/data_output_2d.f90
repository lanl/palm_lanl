!> @file data_output_2d.f90
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
! $Id: data_output_2d.f90 3052 2018-05-31 06:11:20Z raasch $
! Do not open FORTRAN binary files in case of parallel netCDF I/O
! 
! 3049 2018-05-29 13:52:36Z Giersch
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
! 
! 3014 2018-05-09 08:42:38Z maronga
! Added nzb_do and nzt_do for some modules for 2d output
! 
! 3004 2018-04-27 12:33:25Z Giersch
! precipitation_rate removed, case prr*_xy removed, to_be_resorted have to point
! to ql_vp_av and not to ql_vp, allocation checks implemented (averaged data 
! will be assigned to fill values if no allocation happened so far)   
! 
! 2963 2018-04-12 14:47:44Z suehring
! Introduce index for vegetation/wall, pavement/green-wall and water/window 
! surfaces, for clearer access of surface fraction, albedo, emissivity, etc. .
! 
! 2817 2018-02-19 16:32:21Z knoop
! Preliminary gust module interface implemented
! 
! 2805 2018-02-14 17:00:09Z suehring
! Consider also default-type surfaces for surface temperature output. 
! 
! 2797 2018-02-08 13:24:35Z suehring
! Enable output of ground-heat flux also at urban surfaces.
! 
! 2743 2018-01-12 16:03:39Z suehring
! In case of natural- and urban-type surfaces output surfaces fluxes in W/m2.
! 
! 2742 2018-01-12 14:59:47Z suehring
! Enable output of surface temperature
! 
! 2735 2018-01-11 12:01:27Z suehring
! output of r_a moved from land-surface to consider also urban-type surfaces
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Implementation of uv exposure model (FK)
! Implementation of turbulence_closure_mod (TG)
! Set fill values at topography grid points or e.g. non-natural-type surface
! in case of LSM output (MS)
! 
! 2512 2017-10-04 08:26:59Z raasch
! upper bounds of cross section output changed from nx+1,ny+1 to nx,ny
! no output of ghost layer data
! 
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison' 
! includes two more prognostic equations for cloud drop concentration (nc)  
! and cloud water content (qc). 
! 
! 2277 2017-06-12 10:47:51Z kanani
! Removed unused variables do2d_xy_n, do2d_xz_n, do2d_yz_n
! 
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new surface concept
! 
!
! 2190 2017-03-21 12:16:43Z raasch
! bugfix for r2031: string rho replaced by rho_ocean
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean and rho_av to rho_ocean_av
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1980 2016-07-29 15:51:57Z suehring
! Bugfix, in order to steer user-defined output, setting flag found explicitly
! to .F.
! 
! 1976 2016-07-27 13:28:04Z maronga
! Output of radiation quantities is now done directly in the respective module
! 
! 1972 2016-07-26 07:52:02Z maronga
! Output of land surface quantities is now done directly in the respective 
! module
! 
! 1960 2016-07-12 16:34:24Z suehring
! Scalar surface flux added
! Rename INTEGER variable s into s_ind, as s is already assigned to scalar
!
! 1849 2016-04-08 11:33:18Z hoffmann
! precipitation_amount, precipitation_rate, prr moved to arrays_3d
!
! 1822 2016-04-07 07:49:42Z hoffmann
! Output of bulk cloud physics simplified.
!
! 1788 2016-03-10 11:01:04Z maronga
! Added output of z0q
! 
! 1783 2016-03-06 18:36:17Z raasch
! name change of netcdf routines and module + related changes
!
! 1745 2016-02-05 13:06:51Z gronemeier
! Bugfix: test if time axis limit exceeds moved to point after call of check_open
!
! 1703 2015-11-02 12:38:44Z raasch
! bugfix for output of single (*) xy-sections in case of parallel netcdf I/O
!
! 1701 2015-11-02 07:43:04Z maronga
! Bugfix in output of RRTGM data
! 
! 1691 2015-10-26 16:17:44Z maronga
! Added output of Obukhov length (ol) and radiative heating rates  for RRTMG. 
! Formatting corrections.
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1585 2015-04-30 07:05:52Z maronga
! Added support for RRTMG
! 
! 1555 2015-03-04 17:44:27Z maronga
! Added output of r_a and r_s
! 
! 1551 2015-03-03 14:18:16Z maronga
! Added suppport for land surface model and radiation model output. In the course
! of this action, the limits for vertical loops have been changed (from nzb and 
! nzt+1 to nzb_do and nzt_do, respectively in order to allow soil model output).
! Moreover, a new vertical grid zs was introduced.
! 
! 1359 2014-04-11 17:15:14Z hoffmann
! New particle structure integrated. 
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
! 
! 1327 2014-03-21 11:00:16Z raasch
! parts concerning iso2d output removed,
! -netcdf output queries
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1318 2014-03-17 13:35:16Z raasch
! barrier argument removed from cpu_log.
! module interfaces removed
!
! 1311 2014-03-14 12:13:39Z heinze
! bugfix: close #if defined( __netcdf )
!
! 1308 2014-03-13 14:58:42Z fricke
! +local_2d_sections, local_2d_sections_l, ns
! Check, if the limit of the time dimension is exceeded for parallel output
! To increase the performance for parallel output, the following is done:
! - Update of time axis is only done by PE0
! - Cross sections are first stored on a local array and are written
!   collectively to the output file by all PEs.
!
! 1115 2013-03-26 18:16:16Z hoffmann
! ql is calculated by calc_liquid_water_content
!
! 1076 2012-12-05 08:30:18Z hoffmann
! Bugfix in output of ql 
!
! 1065 2012-11-22 17:42:36Z hoffmann
! Bugfix: Output of cross sections of ql
!
! 1053 2012-11-13 17:11:03Z hoffmann
! +qr, nr, qc and cross sections
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1031 2012-10-19 14:35:30Z raasch
! netCDF4 without parallel file support implemented
!
! 1007 2012-09-19 14:30:36Z franke
! Bugfix: missing calculation of ql_vp added
!
! 978 2012-08-09 08:28:32Z fricke
! +z0h
!
! Revision 1.1  1997/08/11 06:24:09  raasch
! Initial revision
!
!
! Description:
! ------------
!> Data output of cross-sections in netCDF format or binary format
!> to be later converted to NetCDF by helper routine combine_plot_fields.
!> Attention: The position of the sectional planes is still not always computed 
!> ---------  correctly. (zu is used always)!
!------------------------------------------------------------------------------!
 SUBROUTINE data_output_2d( mode, av )
 

    USE arrays_3d,                                                             &
        ONLY:  dzw, e, heatflux_output_conversion, nc, nr, p, pt,              &
               precipitation_amount, prr, q, qc, ql, ql_c, ql_v, ql_vp, qr,    &
               rho_ocean, s, sa, tend, u, v, vpt, w, zu, zw,                   &
               salinityflux_output_conversion, waterflux_output_conversion
        
    USE averaging
        
    USE cloud_parameters,                                                      &
        ONLY:  cp, hyrho, l_d_cp, l_v, pt_d_t
               
    USE control_parameters,                                                    &
        ONLY:  cloud_physics, constant_flux_layer, data_output_2d_on_each_pe,  &
               data_output_xy, data_output_xz, data_output_yz, do2d,           &
               do2d_xy_last_time, do2d_xy_time_count,                          &
               do2d_xz_last_time, do2d_xz_time_count,                          &
               do2d_yz_last_time, do2d_yz_time_count,                          &
               ibc_uv_b, io_blocks, io_group, land_surface, message_string,    &
               ntdim_2d_xy, ntdim_2d_xz, ntdim_2d_yz,                          &
               psolver, section, simulated_time, simulated_time_chr,           &
               time_since_reference_point, uv_exposure
        
    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point 
        
    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE gust_mod,                                                              &
        ONLY:  gust_data_output_2d, gust_module_enabled
        
    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg, nz,   &
               nzb, nzt, wall_flags_0
               
    USE kinds
    
    USE land_surface_model_mod,                                                &
        ONLY:  lsm_data_output_2d, zs
    
#if defined( __netcdf )
    USE NETCDF
#endif

    USE netcdf_interface,                                                      &
        ONLY:  fill_value, id_set_xy, id_set_xz, id_set_yz, id_var_do2d,       &
               id_var_time_xy, id_var_time_xz, id_var_time_yz, nc_stat,        &
               netcdf_data_format, netcdf_handle_error

    USE particle_attributes,                                                   &
        ONLY:  grid_particles, number_of_particles, particle_advection_start,  &
               particles, prt_count
    
    USE pegrid

    USE radiation_model_mod,                                                   &
        ONLY:  radiation, radiation_data_output_2d

    USE surface_mod,                                                           &
        ONLY:  ind_pav_green, ind_veg_wall, ind_wat_win, surf_def_h,           &
               surf_lsm_h, surf_usm_h

    USE turbulence_closure_mod,                                                &
        ONLY:  tcm_data_output_2d

    USE uv_exposure_model_mod,                                                 &
        ONLY:  uvem_data_output_2d


    IMPLICIT NONE

    CHARACTER (LEN=2)  ::  do2d_mode    !< 
    CHARACTER (LEN=2)  ::  mode         !< 
    CHARACTER (LEN=4)  ::  grid         !< 
    CHARACTER (LEN=25) ::  section_chr  !< 
    CHARACTER (LEN=50) ::  rtext        !< 
    
    INTEGER(iwp) ::  av        !< 
    INTEGER(iwp) ::  ngp       !< 
    INTEGER(iwp) ::  file_id   !< 
    INTEGER(iwp) ::  flag_nr   !< number of masking flag
    INTEGER(iwp) ::  i         !< 
    INTEGER(iwp) ::  if        !< 
    INTEGER(iwp) ::  is        !< 
    INTEGER(iwp) ::  iis       !< 
    INTEGER(iwp) ::  j         !< 
    INTEGER(iwp) ::  k         !< 
    INTEGER(iwp) ::  l         !< 
    INTEGER(iwp) ::  layer_xy  !< 
    INTEGER(iwp) ::  m         !<
    INTEGER(iwp) ::  n         !< 
    INTEGER(iwp) ::  nis       !<
    INTEGER(iwp) ::  ns        !< 
    INTEGER(iwp) ::  nzb_do    !< lower limit of the data field (usually nzb)
    INTEGER(iwp) ::  nzt_do    !< upper limit of the data field (usually nzt+1)
    INTEGER(iwp) ::  psi       !< 
    INTEGER(iwp) ::  s_ind     !< 
    INTEGER(iwp) ::  sender    !< 
    INTEGER(iwp) ::  ind(4)    !< 
    
    LOGICAL ::  found          !< 
    LOGICAL ::  resorted       !< 
    LOGICAL ::  two_d          !< 
    
    REAL(wp) ::  mean_r        !<
    REAL(wp) ::  s_r2          !<
    REAL(wp) ::  s_r3          !<
    
    REAL(wp), DIMENSION(:), ALLOCATABLE     ::  level_z             !<
    REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  local_2d            !<
    REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  local_2d_l          !<

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  local_pf            !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  local_2d_sections   !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  local_2d_sections_l !<

#if defined( __parallel )
    REAL(wp), DIMENSION(:,:),   ALLOCATABLE ::  total_2d    !<
#endif
    REAL(wp), DIMENSION(:,:,:), POINTER ::  to_be_resorted  !< 

    NAMELIST /LOCAL/  rtext

!
!-- Immediate return, if no output is requested (no respective sections
!-- found in parameter data_output)
    IF ( mode == 'xy'  .AND.  .NOT. data_output_xy(av) )  RETURN
    IF ( mode == 'xz'  .AND.  .NOT. data_output_xz(av) )  RETURN
    IF ( mode == 'yz'  .AND.  .NOT. data_output_yz(av) )  RETURN

    CALL cpu_log (log_point(3),'data_output_2d','start')

    two_d = .FALSE.    ! local variable to distinguish between output of pure 2D
                       ! arrays and cross-sections of 3D arrays.

!
!-- Depending on the orientation of the cross-section, the respective output 
!-- files have to be opened.
    SELECT CASE ( mode )

       CASE ( 'xy' )
          s_ind = 1
          ALLOCATE( level_z(nzb:nzt+1), local_2d(nxl:nxr,nys:nyn) )

          IF ( netcdf_data_format > 4 )  THEN
             ns = 1
             DO WHILE ( section(ns,s_ind) /= -9999  .AND.  ns <= 100 )
                ns = ns + 1
             ENDDO
             ns = ns - 1
             ALLOCATE( local_2d_sections(nxl:nxr,nys:nyn,1:ns) )
             local_2d_sections = 0.0_wp
          ENDIF

!
!--       Parallel netCDF4/HDF5 output is done on all PEs, all other on PE0 only
          IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
             CALL check_open( 101+av*10 )
          ENDIF
          IF ( data_output_2d_on_each_pe  .AND.  netcdf_data_format < 5 )  THEN
             CALL check_open( 21 )
          ELSE
             IF ( myid == 0 )  THEN
#if defined( __parallel )
                ALLOCATE( total_2d(0:nx,0:ny) )
#endif
             ENDIF
          ENDIF

       CASE ( 'xz' )
          s_ind = 2
          ALLOCATE( local_2d(nxl:nxr,nzb:nzt+1) )

          IF ( netcdf_data_format > 4 )  THEN
             ns = 1
             DO WHILE ( section(ns,s_ind) /= -9999  .AND.  ns <= 100 )
                ns = ns + 1
             ENDDO
             ns = ns - 1
             ALLOCATE( local_2d_sections(nxl:nxr,1:ns,nzb:nzt+1) )
             ALLOCATE( local_2d_sections_l(nxl:nxr,1:ns,nzb:nzt+1) )
             local_2d_sections = 0.0_wp; local_2d_sections_l = 0.0_wp
          ENDIF

!
!--       Parallel netCDF4/HDF5 output is done on all PEs, all other on PE0 only
          IF ( myid == 0 .OR. netcdf_data_format > 4 )  THEN
             CALL check_open( 102+av*10 )
          ENDIF

          IF ( data_output_2d_on_each_pe  .AND.  netcdf_data_format < 5 )  THEN
             CALL check_open( 22 )
          ELSE
             IF ( myid == 0 )  THEN
#if defined( __parallel )
                ALLOCATE( total_2d(0:nx,nzb:nzt+1) )
#endif
             ENDIF
          ENDIF

       CASE ( 'yz' )
          s_ind = 3
          ALLOCATE( local_2d(nys:nyn,nzb:nzt+1) )

          IF ( netcdf_data_format > 4 )  THEN
             ns = 1
             DO WHILE ( section(ns,s_ind) /= -9999  .AND.  ns <= 100 )
                ns = ns + 1
             ENDDO
             ns = ns - 1
             ALLOCATE( local_2d_sections(1:ns,nys:nyn,nzb:nzt+1) )
             ALLOCATE( local_2d_sections_l(1:ns,nys:nyn,nzb:nzt+1) )
             local_2d_sections = 0.0_wp; local_2d_sections_l = 0.0_wp
          ENDIF

!
!--       Parallel netCDF4/HDF5 output is done on all PEs, all other on PE0 only
          IF ( myid == 0 .OR. netcdf_data_format > 4 )  THEN
             CALL check_open( 103+av*10 )
          ENDIF

          IF ( data_output_2d_on_each_pe  .AND.  netcdf_data_format < 5 )  THEN
             CALL check_open( 23 )
          ELSE
             IF ( myid == 0 )  THEN
#if defined( __parallel )
                ALLOCATE( total_2d(0:ny,nzb:nzt+1) )
#endif
             ENDIF
          ENDIF

       CASE DEFAULT
          message_string = 'unknown cross-section: ' // TRIM( mode )
          CALL message( 'data_output_2d', 'PA0180', 1, 2, 0, 6, 0 )

    END SELECT

!
!-- For parallel netcdf output the time axis must be limited. Return, if this
!-- limit is exceeded. This could be the case, if the simulated time exceeds 
!-- the given end time by the length of the given output interval.
    IF ( netcdf_data_format > 4 )  THEN
       IF ( mode == 'xy'  .AND.  do2d_xy_time_count(av) + 1 >                  &
            ntdim_2d_xy(av) )  THEN
          WRITE ( message_string, * ) 'Output of xy cross-sections is not ',   &
                          'given at t=', simulated_time, '&because the',       & 
                          ' maximum number of output time levels is exceeded.'
          CALL message( 'data_output_2d', 'PA0384', 0, 1, 0, 6, 0 )
          CALL cpu_log( log_point(3), 'data_output_2d', 'stop' )
          RETURN
       ENDIF
       IF ( mode == 'xz'  .AND.  do2d_xz_time_count(av) + 1 >                  &
            ntdim_2d_xz(av) )  THEN
          WRITE ( message_string, * ) 'Output of xz cross-sections is not ',   &
                          'given at t=', simulated_time, '&because the',       & 
                          ' maximum number of output time levels is exceeded.'
          CALL message( 'data_output_2d', 'PA0385', 0, 1, 0, 6, 0 )
          CALL cpu_log( log_point(3), 'data_output_2d', 'stop' )
          RETURN
       ENDIF
       IF ( mode == 'yz'  .AND.  do2d_yz_time_count(av) + 1 >                  &
            ntdim_2d_yz(av) )  THEN
          WRITE ( message_string, * ) 'Output of yz cross-sections is not ',   &
                          'given at t=', simulated_time, '&because the',       & 
                          ' maximum number of output time levels is exceeded.'
          CALL message( 'data_output_2d', 'PA0386', 0, 1, 0, 6, 0 )
          CALL cpu_log( log_point(3), 'data_output_2d', 'stop' )
          RETURN
       ENDIF
    ENDIF

!
!-- Allocate a temporary array for resorting (kji -> ijk).
    ALLOCATE( local_pf(nxl:nxr,nys:nyn,nzb:nzt+1) )
    local_pf = 0.0

!
!-- Loop of all variables to be written.
!-- Output dimensions chosen
    if = 1
    l = MAX( 2, LEN_TRIM( do2d(av,if) ) )
    do2d_mode = do2d(av,if)(l-1:l)

    DO  WHILE ( do2d(av,if)(1:1) /= ' ' )

       IF ( do2d_mode == mode )  THEN
!
!--       Set flag to steer output of radiation, land-surface, or user-defined
!--       quantities
          found = .FALSE.

          nzb_do = nzb
          nzt_do = nzt+1
!
!--       Before each output, set array local_pf to fill value
          local_pf = fill_value
!
!--       Set masking flag for topography for not resorted arrays
          flag_nr = 0
          
!
!--       Store the array chosen on the temporary array.
          resorted = .FALSE.
          SELECT CASE ( TRIM( do2d(av,if) ) )
             CASE ( 'e_xy', 'e_xz', 'e_yz' )
                IF ( av == 0 )  THEN
                   to_be_resorted => e
                ELSE
                   IF ( .NOT. ALLOCATED( e_av ) ) THEN
                      ALLOCATE( e_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      e_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   to_be_resorted => e_av
                ENDIF
                IF ( mode == 'xy' )  level_z = zu

             CASE ( 'lpt_xy', 'lpt_xz', 'lpt_yz' )
                IF ( av == 0 )  THEN
                   to_be_resorted => pt
                ELSE
                   IF ( .NOT. ALLOCATED( lpt_av ) ) THEN
                      ALLOCATE( lpt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      lpt_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   to_be_resorted => lpt_av
                ENDIF
                IF ( mode == 'xy' )  level_z = zu

             CASE ( 'lwp*_xy' )        ! 2d-array
                IF ( av == 0 )  THEN
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         local_pf(i,j,nzb+1) = SUM( ql(nzb:nzt,j,i) *          &
                                                    dzw(1:nzt+1) )
                      ENDDO
                   ENDDO
                ELSE
                   IF ( .NOT. ALLOCATED( lwp_av ) ) THEN
                      ALLOCATE( lwp_av(nysg:nyng,nxlg:nxrg) )
                      lwp_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         local_pf(i,j,nzb+1) = lwp_av(j,i)
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)

             CASE ( 'melt*_xy' )        ! 2d-array
                IF ( TRIM(constant_flux_layer) == 'top') THEN
                   IF ( av == 0 ) THEN
                      DO  m = 1, surf_def_h(2)%ns
                         i = surf_def_h(2)%i(m)
                         j = surf_def_h(2)%j(m)
                         local_pf(i,j,nzb+1) = surf_def_h(2)%melt(m)
                      ENDDO
                   ELSE
                      IF ( .NOT. ALLOCATED( melt_av ) ) THEN
                         ALLOCATE( melt_av(nysg:nyng,nxlg:nxrg) )
                         melt_av = REAL( fill_value, KIND = wp )
                      ENDIF
                      DO  i = nxl, nxr
                         DO  j = nys, nyn
                            local_pf(i,j,nzb+1) =  melt_av(j,i)
                         ENDDO
                      ENDDO
                   ENDIF
                ENDIF
                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)
             
             CASE ( 'nc_xy', 'nc_xz', 'nc_yz' )
                IF ( av == 0 )  THEN
                   to_be_resorted => nc
                ELSE
                   IF ( .NOT. ALLOCATED( nc_av ) ) THEN
                      ALLOCATE( nc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      nc_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   to_be_resorted => nc_av
                ENDIF
                IF ( mode == 'xy' )  level_z = zu

             CASE ( 'nr_xy', 'nr_xz', 'nr_yz' )
                IF ( av == 0 )  THEN
                   to_be_resorted => nr
                ELSE
                   IF ( .NOT. ALLOCATED( nr_av ) ) THEN
                      ALLOCATE( nr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      nr_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   to_be_resorted => nr_av
                ENDIF
                IF ( mode == 'xy' )  level_z = zu

             CASE ( 'ghf*_xy' )        ! 2d-array
                IF ( av == 0 )  THEN
                   DO  m = 1, surf_lsm_h%ns
                      i                   = surf_lsm_h%i(m)            
                      j                   = surf_lsm_h%j(m)
                      local_pf(i,j,nzb+1) = surf_lsm_h%ghf(m)
                   ENDDO
                   DO  m = 1, surf_usm_h%ns
                      i                   = surf_usm_h%i(m)            
                      j                   = surf_usm_h%j(m)
                      local_pf(i,j,nzb+1) = surf_usm_h%frac(ind_veg_wall,m)  *  &
                                            surf_usm_h%wghf_eb(m)        +      &
                                            surf_usm_h%frac(ind_pav_green,m) *  &
                                            surf_usm_h%wghf_eb_green(m)  +      &
                                            surf_usm_h%frac(ind_wat_win,m)   *  &
                                            surf_usm_h%wghf_eb_window(m)
                   ENDDO
                ELSE
                   IF ( .NOT. ALLOCATED( ghf_av ) ) THEN
                      ALLOCATE( ghf_av(nysg:nyng,nxlg:nxrg) )
                      ghf_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         local_pf(i,j,nzb+1) = ghf_av(j,i)
                      ENDDO
                   ENDDO
                ENDIF

                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)

             CASE ( 'ol*_xy' )        ! 2d-array
                IF ( TRIM(constant_flux_layer) == 'top') THEN
                   IF ( av == 0 ) THEN
                      DO  m = 1, surf_def_h(2)%ns
                         i = surf_def_h(2)%i(m)
                         j = surf_def_h(2)%j(m)
                         local_pf(i,j,nzb+1) = surf_def_h(2)%ol(m)
                      ENDDO
                   ELSE
                      IF ( .NOT. ALLOCATED( ol_av ) ) THEN
                         ALLOCATE( ol_av(nysg:nyng,nxlg:nxrg) )
                         ol_av = REAL( fill_value, KIND = wp )
                      ENDIF
                      DO  i = nxl, nxr
                         DO  j = nys, nyn
                            local_pf(i,j,nzb+1) = ol_av(j,i)
                         ENDDO
                      ENDDO
                   ENDIF
                ELSE
                   IF ( av == 0 ) THEN
                      DO  m = 1, surf_def_h(0)%ns
                         i = surf_def_h(0)%i(m)
                         j = surf_def_h(0)%j(m)
                         local_pf(i,j,nzb+1) = surf_def_h(0)%ol(m)
                      ENDDO
                      DO  m = 1, surf_lsm_h%ns
                         i = surf_lsm_h%i(m)
                         j = surf_lsm_h%j(m)
                         local_pf(i,j,nzb+1) = surf_lsm_h%ol(m)
                      ENDDO
                      DO  m = 1, surf_usm_h%ns
                         i = surf_usm_h%i(m)
                         j = surf_usm_h%j(m)
                         local_pf(i,j,nzb+1) = surf_usm_h%ol(m)
                      ENDDO
                   ELSE
                      IF ( .NOT. ALLOCATED( ol_av ) ) THEN
                         ALLOCATE( ol_av(nysg:nyng,nxlg:nxrg) )
                         ol_av = REAL( fill_value, KIND = wp )
                      ENDIF
                      DO  i = nxl, nxr
                         DO  j = nys, nyn
                            local_pf(i,j,nzb+1) = ol_av(j,i)
                         ENDDO
                      ENDDO
                   ENDIF
                ENDIF
                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)

             CASE ( 'p_xy', 'p_xz', 'p_yz' )
                IF ( av == 0 )  THEN
                   IF ( psolver /= 'sor' )  CALL exchange_horiz( p, nbgp )
                   to_be_resorted => p
                ELSE
                   IF ( .NOT. ALLOCATED( p_av ) ) THEN
                      ALLOCATE( p_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      p_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   IF ( psolver /= 'sor' )  CALL exchange_horiz( p_av, nbgp )
                   to_be_resorted => p_av
                ENDIF
                IF ( mode == 'xy' )  level_z = zu

             CASE ( 'pc_xy', 'pc_xz', 'pc_yz' )  ! particle concentration
                IF ( av == 0 )  THEN
                   IF ( simulated_time >= particle_advection_start )  THEN
                      tend = prt_count
!                      CALL exchange_horiz( tend, nbgp )
                   ELSE
                      tend = 0.0_wp
                   ENDIF
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         DO  k = nzb, nzt+1
                            local_pf(i,j,k) = tend(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                   resorted = .TRUE.
                ELSE
                   IF ( .NOT. ALLOCATED( pc_av ) ) THEN
                      ALLOCATE( pc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      pc_av = REAL( fill_value, KIND = wp )
                   ENDIF
!                   CALL exchange_horiz( pc_av, nbgp )
                   to_be_resorted => pc_av
                ENDIF

             CASE ( 'pr_xy', 'pr_xz', 'pr_yz' )  ! mean particle radius (effective radius)
                IF ( av == 0 )  THEN
                   IF ( simulated_time >= particle_advection_start )  THEN
                      DO  i = nxl, nxr
                         DO  j = nys, nyn
                            DO  k = nzb, nzt+1
                               number_of_particles = prt_count(k,j,i)
                               IF (number_of_particles <= 0)  CYCLE
                               particles => grid_particles(k,j,i)%particles(1:number_of_particles)
                               s_r2 = 0.0_wp
                               s_r3 = 0.0_wp
                               DO  n = 1, number_of_particles
                                  IF ( particles(n)%particle_mask )  THEN
                                     s_r2 = s_r2 + particles(n)%radius**2 * &
                                            particles(n)%weight_factor
                                     s_r3 = s_r3 + particles(n)%radius**3 * &
                                            particles(n)%weight_factor
                                  ENDIF
                               ENDDO
                               IF ( s_r2 > 0.0_wp )  THEN
                                  mean_r = s_r3 / s_r2
                               ELSE
                                  mean_r = 0.0_wp
                               ENDIF
                               tend(k,j,i) = mean_r
                            ENDDO
                         ENDDO
                      ENDDO
!                      CALL exchange_horiz( tend, nbgp )
                   ELSE
                      tend = 0.0_wp
                   ENDIF
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         DO  k = nzb, nzt+1
                            local_pf(i,j,k) = tend(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                   resorted = .TRUE.
                ELSE
                   IF ( .NOT. ALLOCATED( pr_av ) ) THEN
                      ALLOCATE( pr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      pr_av = REAL( fill_value, KIND = wp )
                   ENDIF
!                   CALL exchange_horiz( pr_av, nbgp )
                   to_be_resorted => pr_av
                ENDIF

             CASE ( 'pra*_xy' )        ! 2d-array / integral quantity => no av
!                CALL exchange_horiz_2d( precipitation_amount )
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                      local_pf(i,j,nzb+1) =  precipitation_amount(j,i)
                   ENDDO
                ENDDO
                precipitation_amount = 0.0_wp   ! reset for next integ. interval
                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)

             CASE ( 'prr_xy', 'prr_xz', 'prr_yz' )
                IF ( av == 0 )  THEN
!                   CALL exchange_horiz( prr, nbgp )
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         DO  k = nzb, nzt+1
                            local_pf(i,j,k) = prr(k,j,i) * hyrho(nzb+1)
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
                   IF ( .NOT. ALLOCATED( prr_av ) ) THEN
                      ALLOCATE( prr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      prr_av = REAL( fill_value, KIND = wp )
                   ENDIF
!                   CALL exchange_horiz( prr_av, nbgp )
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         DO  k = nzb, nzt+1
                            local_pf(i,j,k) = prr_av(k,j,i) * hyrho(nzb+1)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
                IF ( mode == 'xy' )  level_z = zu

             CASE ( 'pt1*_xy' )        ! 2d-array
                IF ( TRIM(constant_flux_layer) == 'top') THEN
                   IF ( av == 0 ) THEN
                      DO  m = 1, surf_def_h(2)%ns
                         i = surf_def_h(2)%i(m)
                         j = surf_def_h(2)%j(m)
                         k = surf_def_h(2)%k(m)
                         local_pf(i,j,nzb+1) = surf_def_h(2)%pt1(m)
                      ENDDO
                   ELSE
                      IF ( .NOT. ALLOCATED( pt1_av ) ) THEN
                         ALLOCATE( pt1_av(nysg:nyng,nxlg:nxrg) )
                         pt1_av = REAL( fill_value, KIND = wp )
                      ENDIF
                      DO  i = nxl, nxr
                         DO  j = nys, nyn
                            local_pf(i,j,nzb+1) =  pt1_av(j,i)
                         ENDDO
                      ENDDO
                   ENDIF
                ENDIF
                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)

             CASE ( 'pt_io*_xy' )        ! 2d-array
                IF ( TRIM(constant_flux_layer) == 'top') THEN
                   IF ( av == 0 ) THEN
                      DO  m = 1, surf_def_h(2)%ns
                         i = surf_def_h(2)%i(m)
                         j = surf_def_h(2)%j(m)
                         k = surf_def_h(2)%k(m)
                         local_pf(i,j,nzb+1) = surf_def_h(2)%pt_io(m)
                      ENDDO
                   ELSE
                      IF ( .NOT. ALLOCATED( pt_io_av ) ) THEN
                         ALLOCATE( pt_io_av(nysg:nyng,nxlg:nxrg) )
                         pt_io_av = REAL( fill_value, KIND = wp )
                      ENDIF
                      DO  i = nxl, nxr
                         DO  j = nys, nyn
                            local_pf(i,j,nzb+1) =  pt_io_av(j,i)
                         ENDDO
                      ENDDO
                   ENDIF
                ENDIF
                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)

             CASE ( 'pt_xy', 'pt_xz', 'pt_yz' )
                IF ( av == 0 )  THEN
                   IF ( .NOT. cloud_physics ) THEN
                      to_be_resorted => pt
                   ELSE
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                            DO  k = nzb, nzt+1
                               local_pf(i,j,k) = pt(k,j,i) + l_d_cp *          &
                                                             pt_d_t(k) *       &
                                                             ql(k,j,i)
                            ENDDO
                         ENDDO
                      ENDDO
                      resorted = .TRUE.
                   ENDIF
                ELSE
                   IF ( .NOT. ALLOCATED( pt_av ) ) THEN
                      ALLOCATE( pt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      pt_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   to_be_resorted => pt_av
                ENDIF
                IF ( mode == 'xy' )  level_z = zu

             CASE ( 'q_xy', 'q_xz', 'q_yz' )
                IF ( av == 0 )  THEN
                   to_be_resorted => q
                ELSE
                   IF ( .NOT. ALLOCATED( q_av ) ) THEN
                      ALLOCATE( q_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      q_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   to_be_resorted => q_av
                ENDIF
                IF ( mode == 'xy' )  level_z = zu

             CASE ( 'qc_xy', 'qc_xz', 'qc_yz' )
                IF ( av == 0 )  THEN
                   to_be_resorted => qc
                ELSE
                   IF ( .NOT. ALLOCATED( qc_av ) ) THEN
                      ALLOCATE( qc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      qc_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   to_be_resorted => qc_av
                ENDIF
                IF ( mode == 'xy' )  level_z = zu

             CASE ( 'ql_xy', 'ql_xz', 'ql_yz' )
                IF ( av == 0 )  THEN
                   to_be_resorted => ql
                ELSE
                   IF ( .NOT. ALLOCATED( ql_av ) ) THEN
                      ALLOCATE( ql_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      ql_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   to_be_resorted => ql_av
                ENDIF
                IF ( mode == 'xy' )  level_z = zu

             CASE ( 'ql_c_xy', 'ql_c_xz', 'ql_c_yz' )
                IF ( av == 0 )  THEN
                   to_be_resorted => ql_c
                ELSE
                   IF ( .NOT. ALLOCATED( ql_c_av ) ) THEN
                      ALLOCATE( ql_c_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      ql_c_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   to_be_resorted => ql_c_av
                ENDIF
                IF ( mode == 'xy' )  level_z = zu

             CASE ( 'ql_v_xy', 'ql_v_xz', 'ql_v_yz' )
                IF ( av == 0 )  THEN
                   to_be_resorted => ql_v
                ELSE
                   IF ( .NOT. ALLOCATED( ql_v_av ) ) THEN
                      ALLOCATE( ql_v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      ql_v_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   to_be_resorted => ql_v_av
                ENDIF
                IF ( mode == 'xy' )  level_z = zu

             CASE ( 'ql_vp_xy', 'ql_vp_xz', 'ql_vp_yz' )
                IF ( av == 0 )  THEN
                   IF ( simulated_time >= particle_advection_start )  THEN
                      DO  i = nxl, nxr
                         DO  j = nys, nyn
                            DO  k = nzb, nzt+1
                               number_of_particles = prt_count(k,j,i)
                               IF (number_of_particles <= 0)  CYCLE
                               particles => grid_particles(k,j,i)%particles(1:number_of_particles)
                               DO  n = 1, number_of_particles
                                  IF ( particles(n)%particle_mask )  THEN
                                     tend(k,j,i) =  tend(k,j,i) +                 &
                                                    particles(n)%weight_factor /  &
                                                    prt_count(k,j,i)
                                  ENDIF
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
!                      CALL exchange_horiz( tend, nbgp )
                   ELSE
                      tend = 0.0_wp
                   ENDIF
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         DO  k = nzb, nzt+1
                            local_pf(i,j,k) = tend(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                   resorted = .TRUE.
                ELSE
                   IF ( .NOT. ALLOCATED( ql_vp_av ) ) THEN
                      ALLOCATE( ql_vp_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      ql_vp_av = REAL( fill_value, KIND = wp )
                   ENDIF
!                   CALL exchange_horiz( ql_vp_av, nbgp )
                   to_be_resorted => ql_vp_av
                ENDIF
                IF ( mode == 'xy' )  level_z = zu

             CASE ( 'qr_xy', 'qr_xz', 'qr_yz' )
                IF ( av == 0 )  THEN
                   to_be_resorted => qr
                ELSE
                   IF ( .NOT. ALLOCATED( qr_av ) ) THEN
                      ALLOCATE( qr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      qr_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   to_be_resorted => qr_av
                ENDIF
                IF ( mode == 'xy' )  level_z = zu

             CASE ( 'qsws*_xy' )        ! 2d-array
                IF ( av == 0 ) THEN
!
!--                In case of default surfaces, clean-up flux by density.
!--                In case of land- and urban-surfaces, convert fluxes into
!--                dynamic units
                   DO  m = 1, surf_def_h(0)%ns
                      i = surf_def_h(0)%i(m)
                      j = surf_def_h(0)%j(m)
                      k = surf_def_h(0)%k(m)
                      local_pf(i,j,nzb+1) = surf_def_h(0)%qsws(m) *            &
                                            waterflux_output_conversion(k)
                   ENDDO
                   DO  m = 1, surf_lsm_h%ns
                      i = surf_lsm_h%i(m)
                      j = surf_lsm_h%j(m)
                      k = surf_lsm_h%k(m)
                      local_pf(i,j,nzb+1) = surf_lsm_h%qsws(m) * l_v
                   ENDDO
                   DO  m = 1, surf_usm_h%ns
                      i = surf_usm_h%i(m)
                      j = surf_usm_h%j(m)
                      k = surf_usm_h%k(m)
                      local_pf(i,j,nzb+1) = surf_usm_h%qsws(m) * l_v
                   ENDDO
                ELSE
                   IF ( .NOT. ALLOCATED( qsws_av ) ) THEN
                      ALLOCATE( qsws_av(nysg:nyng,nxlg:nxrg) )
                      qsws_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   DO  i = nxl, nxr
                      DO  j = nys, nyn 
                         local_pf(i,j,nzb+1) =  qsws_av(j,i)
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)

             CASE ( 'qv_xy', 'qv_xz', 'qv_yz' )
                IF ( av == 0 )  THEN
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         DO  k = nzb, nzt+1
                            local_pf(i,j,k) = q(k,j,i) - ql(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                   resorted = .TRUE.
                ELSE
                   IF ( .NOT. ALLOCATED( qv_av ) ) THEN
                      ALLOCATE( qv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      qv_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   to_be_resorted => qv_av
                ENDIF
                IF ( mode == 'xy' )  level_z = zu

             CASE ( 'r_a*_xy' )        ! 2d-array
                IF ( av == 0 )  THEN
                   DO  m = 1, surf_lsm_h%ns
                      i                   = surf_lsm_h%i(m)            
                      j                   = surf_lsm_h%j(m)
                      local_pf(i,j,nzb+1) = surf_lsm_h%r_a(m)
                   ENDDO

                   DO  m = 1, surf_usm_h%ns
                      i   = surf_usm_h%i(m)            
                      j   = surf_usm_h%j(m)
                      local_pf(i,j,nzb+1) =                                          &
                                 ( surf_usm_h%frac(ind_veg_wall,m)  *                &
                                   surf_usm_h%r_a(m)       +                         & 
                                   surf_usm_h%frac(ind_pav_green,m) *                &
                                   surf_usm_h%r_a_green(m) +                         & 
                                   surf_usm_h%frac(ind_wat_win,m)   *                &
                                   surf_usm_h%r_a_window(m) )
                   ENDDO
                ELSE
                   IF ( .NOT. ALLOCATED( r_a_av ) ) THEN
                      ALLOCATE( r_a_av(nysg:nyng,nxlg:nxrg) )
                      r_a_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         local_pf(i,j,nzb+1) = r_a_av(j,i)
                      ENDDO
                   ENDDO
                ENDIF
                resorted       = .TRUE.
                two_d          = .TRUE.
                level_z(nzb+1) = zu(nzb+1)

             CASE ( 'rho_ocean_xy', 'rho_ocean_xz', 'rho_ocean_yz' )
                IF ( av == 0 )  THEN
                   to_be_resorted => rho_ocean
                ELSE
                   IF ( .NOT. ALLOCATED( rho_ocean_av ) ) THEN
                      ALLOCATE( rho_ocean_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      rho_ocean_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   to_be_resorted => rho_ocean_av
                ENDIF

             CASE ( 's_xy', 's_xz', 's_yz' )
                IF ( av == 0 )  THEN
                   to_be_resorted => s
                ELSE
                   IF ( .NOT. ALLOCATED( s_av ) ) THEN
                      ALLOCATE( s_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      s_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   to_be_resorted => s_av
                ENDIF

             CASE ( 'sa1*_xy' )        ! 2d-array
                IF ( TRIM(constant_flux_layer) == 'top') THEN
                   IF ( av == 0 ) THEN
                      DO  m = 1, surf_def_h(2)%ns
                         i = surf_def_h(2)%i(m)
                         j = surf_def_h(2)%j(m)
                         k = surf_def_h(2)%k(m)
                         local_pf(i,j,nzb+1) = surf_def_h(2)%sa1(m)
                      ENDDO
                   ELSE
                      IF ( .NOT. ALLOCATED( sa1_av ) ) THEN
                         ALLOCATE( sa1_av(nysg:nyng,nxlg:nxrg) )
                         sa1_av = REAL( fill_value, KIND = wp )
                      ENDIF
                      DO  i = nxl, nxr
                         DO  j = nys, nyn
                            local_pf(i,j,nzb+1) =  sa1_av(j,i)
                         ENDDO
                      ENDDO
                   ENDIF
                ENDIF
                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)

             CASE ( 'sa_io*_xy' )        ! 2d-array
                IF ( TRIM(constant_flux_layer) == 'top') THEN
                   IF ( av == 0 ) THEN
                      DO  m = 1, surf_def_h(2)%ns
                         i = surf_def_h(2)%i(m)
                         j = surf_def_h(2)%j(m)
                         k = surf_def_h(2)%k(m)
                         local_pf(i,j,nzb+1) = surf_def_h(2)%sa_io(m)
                      ENDDO
                   ELSE
                      IF ( .NOT. ALLOCATED( sa_io_av ) ) THEN
                         ALLOCATE( sa_io_av(nysg:nyng,nxlg:nxrg) )
                         sa_io_av = REAL( fill_value, KIND = wp )
                      ENDIF
                      DO  i = nxl, nxr
                         DO  j = nys, nyn
                            local_pf(i,j,nzb+1) =  sa_io_av(j,i)
                         ENDDO
                      ENDDO
                   ENDIF
                ENDIF
                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)

             CASE ( 'sa_xy', 'sa_xz', 'sa_yz' )
                IF ( av == 0 )  THEN
                   to_be_resorted => sa
                ELSE
                   IF ( .NOT. ALLOCATED( sa_av ) ) THEN
                      ALLOCATE( sa_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      sa_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   to_be_resorted => sa_av
                ENDIF

             CASE ( 'shf_sol*_xy' )        ! 2d-array
                IF ( av == 0 ) THEN
!
!--                In case of default surfaces, clean-up flux by density.
!--                In case of land- and urban-surfaces, convert fluxes into
!--                dynamic units.
                   DO  m = 1, surf_def_h(0)%ns
                      i = surf_def_h(0)%i(m)
                      j = surf_def_h(0)%j(m)
                      k = surf_def_h(0)%k(m)
                      local_pf(i,j,nzb+1) = surf_def_h(0)%shf_sol(m) *             &
                                            heatflux_output_conversion(k)
                   ENDDO
                   DO  m = 1, surf_lsm_h%ns
                      i = surf_lsm_h%i(m)
                      j = surf_lsm_h%j(m)
                      k = surf_lsm_h%k(m)
                      local_pf(i,j,nzb+1) = surf_lsm_h%shf_sol(m) * cp
                   ENDDO
                   DO  m = 1, surf_usm_h%ns
                      i = surf_usm_h%i(m)
                      j = surf_usm_h%j(m)
                      k = surf_usm_h%k(m)
                      local_pf(i,j,nzb+1) = surf_usm_h%shf_sol(m) * cp
                   ENDDO
                ELSE
                   IF ( .NOT. ALLOCATED( shf_sol_av ) ) THEN
                      ALLOCATE( shf_sol_av(nysg:nyng,nxlg:nxrg) )
                      shf_sol_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         local_pf(i,j,nzb+1) =  shf_sol_av(j,i)
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)
 
             CASE ( 'shf*_xy' )        ! 2d-array
                IF ( av == 0 ) THEN
!
!--                In case of default surfaces, clean-up flux by density.
!--                In case of land- and urban-surfaces, convert fluxes into
!--                dynamic units.
                   IF ( TRIM(constant_flux_layer) == 'top') THEN
                      DO  m = 1, surf_def_h(2)%ns
                         i = surf_def_h(2)%i(m)
                         j = surf_def_h(2)%j(m)
                         k = surf_def_h(2)%k(m)
                         local_pf(i,j,nzb+1) = surf_def_h(2)%shf(m) *             &
                                               heatflux_output_conversion(k)
                      ENDDO
                   ELSE
                      DO  m = 1, surf_def_h(0)%ns
                         i = surf_def_h(0)%i(m)
                         j = surf_def_h(0)%j(m)
                         k = surf_def_h(0)%k(m)
                         local_pf(i,j,nzb+1) = surf_def_h(0)%shf(m) *             &
                                               heatflux_output_conversion(k)
                      ENDDO
                      DO  m = 1, surf_lsm_h%ns
                         i = surf_lsm_h%i(m)
                         j = surf_lsm_h%j(m)
                         k = surf_lsm_h%k(m)
                         local_pf(i,j,nzb+1) = surf_lsm_h%shf(m) * cp
                      ENDDO
                      DO  m = 1, surf_usm_h%ns
                         i = surf_usm_h%i(m)
                         j = surf_usm_h%j(m)
                         k = surf_usm_h%k(m)
                         local_pf(i,j,nzb+1) = surf_usm_h%shf(m) * cp
                      ENDDO
                   ENDIF
                ELSE
                   IF ( .NOT. ALLOCATED( shf_av ) ) THEN
                      ALLOCATE( shf_av(nysg:nyng,nxlg:nxrg) )
                      shf_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         local_pf(i,j,nzb+1) =  shf_av(j,i)
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)
                
             CASE ( 'sasws*_xy' )        ! 2d-array
                IF ( av == 0 ) THEN
                   IF ( TRIM(constant_flux_layer) == 'top') THEN
                      DO  m = 1, surf_def_h(2)%ns
                         i = surf_def_h(2)%i(m)
                         j = surf_def_h(2)%j(m)
                         k = surf_def_h(2)%k(m)
                         local_pf(i,j,nzb+1) = surf_def_h(2)%sasws(m) *        &
                                               salinityflux_output_conversion(k)
                      ENDDO
                   ELSE
                      DO  m = 1, surf_def_h(0)%ns
                         i = surf_def_h(0)%i(m)
                         j = surf_def_h(0)%j(m)
                         k = surf_def_h(0)%k(m)
                         local_pf(i,j,nzb+1) = surf_def_h(0)%sasws(m) *        &
                                               salinityflux_output_conversion(k)
                      ENDDO
                   ENDIF
                ELSE
                   IF ( .NOT. ALLOCATED( sasws_av ) ) THEN
                      ALLOCATE( sasws_av(nysg:nyng,nxlg:nxrg) )
                      sasws_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         local_pf(i,j,nzb+1) =  sasws_av(j,i)
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)
             
             CASE ( 'ssws*_xy' )        ! 2d-array
                IF ( av == 0 ) THEN
                   DO  m = 1, surf_def_h(0)%ns
                      i = surf_def_h(0)%i(m)
                      j = surf_def_h(0)%j(m)
                      local_pf(i,j,nzb+1) = surf_def_h(0)%ssws(m)
                   ENDDO
                   DO  m = 1, surf_lsm_h%ns
                      i = surf_lsm_h%i(m)
                      j = surf_lsm_h%j(m)
                      local_pf(i,j,nzb+1) = surf_lsm_h%ssws(m)
                   ENDDO
                   DO  m = 1, surf_usm_h%ns
                      i = surf_usm_h%i(m)
                      j = surf_usm_h%j(m)
                      local_pf(i,j,nzb+1) = surf_usm_h%ssws(m)
                   ENDDO
                ELSE
                   IF ( .NOT. ALLOCATED( ssws_av ) ) THEN
                      ALLOCATE( ssws_av(nysg:nyng,nxlg:nxrg) )
                      ssws_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   DO  i = nxl, nxr
                      DO  j = nys, nyn 
                         local_pf(i,j,nzb+1) =  ssws_av(j,i)
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)                

             CASE ( 't*_xy' )        ! 2d-array
                IF ( av == 0 )  THEN
                   DO  m = 1, surf_def_h(0)%ns
                      i = surf_def_h(0)%i(m)
                      j = surf_def_h(0)%j(m)
                      local_pf(i,j,nzb+1) = surf_def_h(0)%ts(m)
                   ENDDO
                   DO  m = 1, surf_lsm_h%ns
                      i = surf_lsm_h%i(m)
                      j = surf_lsm_h%j(m)
                      local_pf(i,j,nzb+1) = surf_lsm_h%ts(m)
                   ENDDO
                   DO  m = 1, surf_usm_h%ns
                      i = surf_usm_h%i(m)
                      j = surf_usm_h%j(m)
                      local_pf(i,j,nzb+1) = surf_usm_h%ts(m)
                   ENDDO
                ELSE
                   IF ( .NOT. ALLOCATED( ts_av ) ) THEN
                      ALLOCATE( ts_av(nysg:nyng,nxlg:nxrg) )
                      ts_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         local_pf(i,j,nzb+1) = ts_av(j,i)
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)

             CASE ( 'tsurf*_xy' )        ! 2d-array
                IF( TRIM(constant_flux_layer) == 'top') THEN
                   IF ( av == 0 )  THEN
                      DO  m = 1, surf_def_h(2)%ns
                         i                   = surf_def_h(2)%i(m)            
                         j                   = surf_def_h(2)%j(m)
                         local_pf(i,j,nzb+1) = surf_def_h(2)%pt_surface(m)
                      ENDDO
                
                   ELSE
                      IF ( .NOT. ALLOCATED( tsurf_av ) ) THEN
                         ALLOCATE( tsurf_av(nysg:nyng,nxlg:nxrg) )
                         tsurf_av = REAL( fill_value, KIND = wp )
                      ENDIF
                      DO  i = nxl, nxr
                         DO  j = nys, nyn
                            local_pf(i,j,nzb+1) = tsurf_av(j,i)
                         ENDDO
                      ENDDO
                   ENDIF
                ELSE
                   IF ( av == 0 )  THEN
                      DO  m = 1, surf_def_h(0)%ns
                         i                   = surf_def_h(0)%i(m)            
                         j                   = surf_def_h(0)%j(m)
                         local_pf(i,j,nzb+1) = surf_def_h(0)%pt_surface(m)
                      ENDDO

                      DO  m = 1, surf_lsm_h%ns
                         i                   = surf_lsm_h%i(m)            
                         j                   = surf_lsm_h%j(m)
                         local_pf(i,j,nzb+1) = surf_lsm_h%pt_surface(m)
                      ENDDO

                      DO  m = 1, surf_usm_h%ns
                         i   = surf_usm_h%i(m)            
                         j   = surf_usm_h%j(m)
                         local_pf(i,j,nzb+1) = surf_usm_h%pt_surface(m)
                      ENDDO

                   ELSE
                      IF ( .NOT. ALLOCATED( tsurf_av ) ) THEN
                         ALLOCATE( tsurf_av(nysg:nyng,nxlg:nxrg) )
                         tsurf_av = REAL( fill_value, KIND = wp )
                      ENDIF
                      DO  i = nxl, nxr
                         DO  j = nys, nyn
                            local_pf(i,j,nzb+1) = tsurf_av(j,i)
                         ENDDO
                      ENDDO
                   ENDIF
                ENDIF
                resorted       = .TRUE.
                two_d          = .TRUE.
                level_z(nzb+1) = zu(nzb+1)

             CASE ( 'u_xy', 'u_xz', 'u_yz' )
                flag_nr = 1
                IF ( av == 0 )  THEN
                   to_be_resorted => u
                ELSE
                   IF ( .NOT. ALLOCATED( u_av ) ) THEN
                      ALLOCATE( u_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      u_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   to_be_resorted => u_av
                ENDIF
                IF ( mode == 'xy' )  level_z = zu
!
!--             Substitute the values generated by "mirror" boundary condition
!--             at the bottom boundary by the real surface values.
                IF ( do2d(av,if) == 'u_xz'  .OR.  do2d(av,if) == 'u_yz' )  THEN
                   IF ( ibc_uv_b == 0 )  local_pf(:,:,nzb) = 0.0_wp
                ENDIF

             CASE ( 'u*_xy' )        ! 2d-array
                IF ( TRIM(constant_flux_layer) == 'top' ) l = 2
                IF ( TRIM(constant_flux_layer) == 'bottom') l = 0
                IF ( av == 0 )  THEN
                   DO  m = 1, surf_def_h(l)%ns
                      i = surf_def_h(l)%i(m)
                      j = surf_def_h(l)%j(m)
                      local_pf(i,j,nzb+1) = surf_def_h(l)%us(m)
                   ENDDO
                ELSE
                   IF ( .NOT. ALLOCATED( us_av ) ) THEN
                      ALLOCATE( us_av(nysg:nyng,nxlg:nxrg) )
                      us_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         local_pf(i,j,nzb+1) = us_av(j,i)
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)

             CASE ( 'usws*_xy' )        ! 2d-array
                IF ( TRIM(constant_flux_layer) == 'top' ) l = 2
                IF ( TRIM(constant_flux_layer) == 'bottom') l = 0
                IF ( av == 0 )  THEN
                   DO  m = 1, surf_def_h(l)%ns
                      i = surf_def_h(l)%i(m)
                      j = surf_def_h(l)%j(m)
                      local_pf(i,j,nzb+1) = surf_def_h(l)%usws(m)
                   ENDDO
                ELSE
                   IF ( .NOT. ALLOCATED( usws_av ) ) THEN
                      ALLOCATE( usws_av(nysg:nyng,nxlg:nxrg) )
                      usws_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         local_pf(i,j,nzb+1) = usws_av(j,i)
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)

             CASE ( 'v_xy', 'v_xz', 'v_yz' )
                flag_nr = 2
                IF ( av == 0 )  THEN
                   to_be_resorted => v
                ELSE
                   IF ( .NOT. ALLOCATED( v_av ) ) THEN
                      ALLOCATE( v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      v_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   to_be_resorted => v_av
                ENDIF
                IF ( mode == 'xy' )  level_z = zu
!
!--             Substitute the values generated by "mirror" boundary condition
!--             at the bottom boundary by the real surface values.
                IF ( do2d(av,if) == 'v_xz'  .OR.  do2d(av,if) == 'v_yz' )  THEN
                   IF ( ibc_uv_b == 0 )  local_pf(:,:,nzb) = 0.0_wp
                ENDIF

             CASE ( 'vsws*_xy' )        ! 2d-array
                IF ( TRIM(constant_flux_layer) == 'top' ) l = 2
                IF ( TRIM(constant_flux_layer) == 'bottom') l = 0
                IF ( av == 0 )  THEN
                   DO  m = 1, surf_def_h(l)%ns
                      i = surf_def_h(l)%i(m)
                      j = surf_def_h(l)%j(m)
                      local_pf(i,j,nzb+1) = surf_def_h(l)%vsws(m)
                   ENDDO
                ELSE
                   IF ( .NOT. ALLOCATED( vsws_av ) ) THEN
                      ALLOCATE( vsws_av(nysg:nyng,nxlg:nxrg) )
                      vsws_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         local_pf(i,j,nzb+1) = vsws_av(j,i)
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)

             CASE ( 'vpt_xy', 'vpt_xz', 'vpt_yz' )
                IF ( av == 0 )  THEN
                   to_be_resorted => vpt
                ELSE
                   IF ( .NOT. ALLOCATED( vpt_av ) ) THEN
                      ALLOCATE( vpt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      vpt_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   to_be_resorted => vpt_av
                ENDIF
                IF ( mode == 'xy' )  level_z = zu

             CASE ( 'w_xy', 'w_xz', 'w_yz' )
                flag_nr = 3
                IF ( av == 0 )  THEN
                   to_be_resorted => w
                ELSE
                   IF ( .NOT. ALLOCATED( w_av ) ) THEN
                      ALLOCATE( w_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                      w_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   to_be_resorted => w_av
                ENDIF
                IF ( mode == 'xy' )  level_z = zw

             CASE ( 'z0*_xy' )        ! 2d-array
                IF ( TRIM(constant_flux_layer) == 'top' ) THEN
                   IF ( av == 0 ) THEN
                      DO  m = 1, surf_def_h(2)%ns
                         i = surf_def_h(2)%i(m)
                         j = surf_def_h(2)%j(m)
                         local_pf(i,j,nzb+1) = surf_def_h(2)%z0(m)
                      ENDDO
                   ELSE
                      IF ( .NOT. ALLOCATED( z0_av ) ) THEN
                         ALLOCATE( z0_av(nysg:nyng,nxlg:nxrg) )
                         z0_av = REAL( fill_value, KIND = wp )
                      ENDIF
                      DO  i = nxl, nxr
                         DO  j = nys, nyn
                            local_pf(i,j,nzb+1) =  z0_av(j,i)
                         ENDDO
                      ENDDO
                   ENDIF
                ELSE
                   IF ( av == 0 ) THEN
                      DO  m = 1, surf_def_h(0)%ns
                         i = surf_def_h(0)%i(m)
                         j = surf_def_h(0)%j(m)
                         local_pf(i,j,nzb+1) = surf_def_h(0)%z0(m)
                      ENDDO
                      DO  m = 1, surf_lsm_h%ns
                         i = surf_lsm_h%i(m)
                         j = surf_lsm_h%j(m)
                         local_pf(i,j,nzb+1) = surf_lsm_h%z0(m)
                      ENDDO
                      DO  m = 1, surf_usm_h%ns
                         i = surf_usm_h%i(m)
                         j = surf_usm_h%j(m)
                         local_pf(i,j,nzb+1) = surf_usm_h%z0(m)
                      ENDDO
                   ELSE
                      IF ( .NOT. ALLOCATED( z0_av ) ) THEN
                         ALLOCATE( z0_av(nysg:nyng,nxlg:nxrg) )
                         z0_av = REAL( fill_value, KIND = wp )
                      ENDIF
                      DO  i = nxl, nxr
                         DO  j = nys, nyn
                            local_pf(i,j,nzb+1) =  z0_av(j,i)
                         ENDDO
                      ENDDO
                   ENDIF
                ENDIF
                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)

             CASE ( 'z0h*_xy' )        ! 2d-array
                IF ( TRIM(constant_flux_layer) == 'top' ) THEN
                   IF ( av == 0 ) THEN
                      DO  m = 1, surf_def_h(2)%ns
                         i = surf_def_h(2)%i(m)
                         j = surf_def_h(2)%j(m)
                         local_pf(i,j,nzb+1) = surf_def_h(2)%z0h(m)
                      ENDDO
                   ELSE
                      IF ( .NOT. ALLOCATED( z0h_av ) ) THEN
                         ALLOCATE( z0h_av(nysg:nyng,nxlg:nxrg) )
                         z0h_av = REAL( fill_value, KIND = wp )
                      ENDIF
                      DO  i = nxl, nxr
                         DO  j = nys, nyn
                            local_pf(i,j,nzb+1) =  z0h_av(j,i)
                         ENDDO
                      ENDDO
                   ENDIF
                ELSE
                   IF ( av == 0 ) THEN
                      DO  m = 1, surf_def_h(0)%ns
                         i = surf_def_h(0)%i(m)
                         j = surf_def_h(0)%j(m)
                         local_pf(i,j,nzb+1) = surf_def_h(0)%z0h(m)
                      ENDDO
                      DO  m = 1, surf_lsm_h%ns
                         i = surf_lsm_h%i(m)
                         j = surf_lsm_h%j(m)
                         local_pf(i,j,nzb+1) = surf_lsm_h%z0h(m)
                      ENDDO
                      DO  m = 1, surf_usm_h%ns
                         i = surf_usm_h%i(m)
                         j = surf_usm_h%j(m)
                         local_pf(i,j,nzb+1) = surf_usm_h%z0h(m)
                      ENDDO
                   ELSE
                      IF ( .NOT. ALLOCATED( z0h_av ) ) THEN
                         ALLOCATE( z0h_av(nysg:nyng,nxlg:nxrg) )
                         z0h_av = REAL( fill_value, KIND = wp )
                      ENDIF
                      DO  i = nxl, nxr
                         DO  j = nys, nyn
                            local_pf(i,j,nzb+1) =  z0h_av(j,i)
                         ENDDO
                      ENDDO
                   ENDIF
                ENDIF
                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)

             CASE ( 'z0q*_xy' )        ! 2d-array
                IF ( av == 0 ) THEN
                   DO  m = 1, surf_def_h(0)%ns
                      i = surf_def_h(0)%i(m)
                      j = surf_def_h(0)%j(m)
                      local_pf(i,j,nzb+1) = surf_def_h(0)%z0q(m)
                   ENDDO
                   DO  m = 1, surf_lsm_h%ns
                      i = surf_lsm_h%i(m)
                      j = surf_lsm_h%j(m)
                      local_pf(i,j,nzb+1) = surf_lsm_h%z0q(m)
                   ENDDO
                   DO  m = 1, surf_usm_h%ns
                      i = surf_usm_h%i(m)
                      j = surf_usm_h%j(m)
                      local_pf(i,j,nzb+1) = surf_usm_h%z0q(m)
                   ENDDO
                ELSE
                   IF ( .NOT. ALLOCATED( z0q_av ) ) THEN
                      ALLOCATE( z0q_av(nysg:nyng,nxlg:nxrg) )
                      z0q_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         local_pf(i,j,nzb+1) =  z0q_av(j,i)
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
                two_d = .TRUE.
                level_z(nzb+1) = zu(nzb+1)

             CASE DEFAULT

!
!--             Land surface model quantity
                IF ( land_surface )  THEN
                   CALL lsm_data_output_2d( av, do2d(av,if), found, grid, mode,&
                                            local_pf, two_d, nzb_do, nzt_do )
                ENDIF

!
!--             Turbulence closure variables
                IF ( .NOT. found )  THEN
                   CALL tcm_data_output_2d( av, do2d(av,if), found, grid, mode,&
                                             local_pf, two_d, nzb_do, nzt_do )
                ENDIF

!
!--             Radiation quantity
                IF ( .NOT. found  .AND.  radiation )  THEN
                   CALL radiation_data_output_2d( av, do2d(av,if), found, grid,&
                                                  mode, local_pf, two_d,       &
                                                  nzb_do, nzt_do  )
                ENDIF

!
!--             Gust module quantities
                IF ( .NOT. found  .AND.  gust_module_enabled )  THEN
                   CALL gust_data_output_2d( av, do2d(av,if), found, grid,     &
                                             local_pf, two_d, nzb_do, nzt_do )
                ENDIF

!
!--             UV exposure model quantity
                IF ( uv_exposure )  THEN
                   CALL uvem_data_output_2d( av, do2d(av,if), found, grid, mode,&
                                             local_pf, two_d, nzb_do, nzt_do )
                ENDIF

!
!--             User defined quantity
                IF ( .NOT. found )  THEN
                   CALL user_data_output_2d( av, do2d(av,if), found, grid,     &
                                             local_pf, two_d, nzb_do, nzt_do )
                ENDIF

                resorted = .TRUE.

                IF ( grid == 'zu' )  THEN
                   IF ( mode == 'xy' )  level_z = zu
                ELSEIF ( grid == 'zw' )  THEN
                   IF ( mode == 'xy' )  level_z = zw
                ELSEIF ( grid == 'zu1' ) THEN
                   IF ( mode == 'xy' )  level_z(nzb+1) = zu(nzb+1)
                ELSEIF ( grid == 'zs' ) THEN
                   IF ( mode == 'xy' )  level_z = zs
                ENDIF

                IF ( .NOT. found )  THEN
                   message_string = 'no output provided for: ' //              &
                                    TRIM( do2d(av,if) )
                   CALL message( 'data_output_2d', 'PA0181', 0, 0, 0, 6, 0 )
                ENDIF

          END SELECT

!
!--       Resort the array to be output, if not done above. Flag topography 
!--       grid points with fill values, using the corresponding maksing flag.
          IF ( .NOT. resorted )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = MERGE( to_be_resorted(k,j,i),          &
                                               REAL( fill_value, KIND = wp ),  &
                                               BTEST( wall_flags_0(k,j,i),     &
                                                      flag_nr ) ) 
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

!
!--       Output of the individual cross-sections, depending on the cross- 
!--       section mode chosen.
          is = 1
   loop1: DO WHILE ( section(is,s_ind) /= -9999  .OR.  two_d )

             SELECT CASE ( mode )

                CASE ( 'xy' )
!
!--                Determine the cross section index
                   IF ( two_d )  THEN
                      layer_xy = nzb+1
                   ELSE
                      layer_xy = section(is,s_ind)
                   ENDIF

!
!--                Exit the loop for layers beyond the data output domain
!--                (used for soil model)
                   IF ( layer_xy > nzt_do )  THEN 
                      EXIT loop1
                   ENDIF

!
!--                Update the netCDF xy cross section time axis.
!--                In case of parallel output, this is only done by PE0
!--                to increase the performance.
                   IF ( simulated_time /= do2d_xy_last_time(av) )  THEN
                      do2d_xy_time_count(av) = do2d_xy_time_count(av) + 1
                      do2d_xy_last_time(av)  = simulated_time
                      IF ( myid == 0 )  THEN
                         IF ( .NOT. data_output_2d_on_each_pe  &
                              .OR.  netcdf_data_format > 4 )   &
                         THEN
#if defined( __netcdf )
                            nc_stat = NF90_PUT_VAR( id_set_xy(av),             &
                                                    id_var_time_xy(av),        &
                                             (/ time_since_reference_point /), &
                                         start = (/ do2d_xy_time_count(av) /), &
                                                    count = (/ 1 /) )
                            CALL netcdf_handle_error( 'data_output_2d', 53 )
#endif
                         ENDIF
                      ENDIF
                   ENDIF
!
!--                If required, carry out averaging along z
                   IF ( section(is,s_ind) == -1  .AND.  .NOT. two_d )  THEN

                      local_2d = 0.0_wp
!
!--                   Carry out the averaging (all data are on the PE)
                      DO  k = nzb_do, nzt_do
                         DO  j = nys, nyn
                            DO  i = nxl, nxr
                               local_2d(i,j) = local_2d(i,j) + local_pf(i,j,k)
                            ENDDO
                         ENDDO
                      ENDDO

                      local_2d = local_2d / ( nzt_do - nzb_do + 1.0_wp)

                   ELSE
!
!--                   Just store the respective section on the local array
                      local_2d = local_pf(:,:,layer_xy)

                   ENDIF

#if defined( __parallel )
                   IF ( netcdf_data_format > 4 )  THEN
!
!--                   Parallel output in netCDF4/HDF5 format.
                      IF ( two_d ) THEN
                         iis = 1
                      ELSE
                         iis = is
                      ENDIF

#if defined( __netcdf )
!
!--                   For parallel output, all cross sections are first stored
!--                   here on a local array and will be written to the output
!--                   file afterwards to increase the performance.
                      DO  i = nxl, nxr
                         DO  j = nys, nyn
                            local_2d_sections(i,j,iis) = local_2d(i,j)
                         ENDDO
                      ENDDO
#endif
                   ELSE

                      IF ( data_output_2d_on_each_pe )  THEN
!
!--                      Output of partial arrays on each PE
#if defined( __netcdf )
                         IF ( myid == 0 )  THEN
                            WRITE ( 21 )  time_since_reference_point,          &
                                          do2d_xy_time_count(av), av
                         ENDIF
#endif
                         DO  i = 0, io_blocks-1
                            IF ( i == io_group )  THEN
                               WRITE ( 21 )  nxl, nxr, nys, nyn, nys, nyn
                               WRITE ( 21 )  local_2d
                            ENDIF
#if defined( __parallel )
                            CALL MPI_BARRIER( comm2d, ierr )
#endif
                         ENDDO

                      ELSE
!
!--                      PE0 receives partial arrays from all processors and
!--                      then outputs them. Here a barrier has to be set,
!--                      because otherwise "-MPI- FATAL: Remote protocol queue
!--                      full" may occur.
                         CALL MPI_BARRIER( comm2d, ierr )

                         ngp = ( nxr-nxl+1 ) * ( nyn-nys+1 )
                         IF ( myid == 0 )  THEN
!
!--                         Local array can be relocated directly.
                            total_2d(nxl:nxr,nys:nyn) = local_2d
!
!--                         Receive data from all other PEs.
                            DO  n = 1, numprocs-1
!
!--                            Receive index limits first, then array.
!--                            Index limits are received in arbitrary order from
!--                            the PEs.
                               CALL MPI_RECV( ind(1), 4, MPI_INTEGER,          &
                                              MPI_ANY_SOURCE, 0, comm2d,       &
                                              status, ierr )
                               sender = status(MPI_SOURCE)
                               DEALLOCATE( local_2d )
                               ALLOCATE( local_2d(ind(1):ind(2),ind(3):ind(4)) )
                               CALL MPI_RECV( local_2d(ind(1),ind(3)), ngp,    &
                                              MPI_REAL, sender, 1, comm2d,     &
                                              status, ierr )
                               total_2d(ind(1):ind(2),ind(3):ind(4)) = local_2d
                            ENDDO
!
!--                         Relocate the local array for the next loop increment
                            DEALLOCATE( local_2d )
                            ALLOCATE( local_2d(nxl:nxr,nys:nyn) )

#if defined( __netcdf )
                            IF ( two_d ) THEN
                               nc_stat = NF90_PUT_VAR( id_set_xy(av),       &
                                                       id_var_do2d(av,if),  &
                                                       total_2d(0:nx,0:ny), &
                             start = (/ 1, 1, 1, do2d_xy_time_count(av) /), &
                                             count = (/ nx+1, ny+1, 1, 1 /) )
                            ELSE
                               nc_stat = NF90_PUT_VAR( id_set_xy(av),       &
                                                       id_var_do2d(av,if),  &
                                                       total_2d(0:nx,0:ny), &
                            start = (/ 1, 1, is, do2d_xy_time_count(av) /), &
                                             count = (/ nx+1, ny+1, 1, 1 /) )
                            ENDIF
                            CALL netcdf_handle_error( 'data_output_2d', 54 )
#endif

                         ELSE
!
!--                         First send the local index limits to PE0
                            ind(1) = nxl; ind(2) = nxr
                            ind(3) = nys; ind(4) = nyn
                            CALL MPI_SEND( ind(1), 4, MPI_INTEGER, 0, 0,       &
                                           comm2d, ierr )
!
!--                         Send data to PE0
                            CALL MPI_SEND( local_2d(nxl,nys), ngp,             &
                                           MPI_REAL, 0, 1, comm2d, ierr )
                         ENDIF
!
!--                      A barrier has to be set, because otherwise some PEs may
!--                      proceed too fast so that PE0 may receive wrong data on
!--                      tag 0
                         CALL MPI_BARRIER( comm2d, ierr )
                      ENDIF

                   ENDIF
#else
#if defined( __netcdf )
                   IF ( two_d ) THEN
                      nc_stat = NF90_PUT_VAR( id_set_xy(av),                &
                                              id_var_do2d(av,if),           &
                                              local_2d(nxl:nxr,nys:nyn),    &
                             start = (/ 1, 1, 1, do2d_xy_time_count(av) /), &
                                           count = (/ nx+1, ny+1, 1, 1 /) )
                   ELSE
                      nc_stat = NF90_PUT_VAR( id_set_xy(av),                &
                                              id_var_do2d(av,if),           &
                                              local_2d(nxl:nxr,nys:nyn),    &
                            start = (/ 1, 1, is, do2d_xy_time_count(av) /), &
                                           count = (/ nx+1, ny+1, 1, 1 /) )
                   ENDIF
                   CALL netcdf_handle_error( 'data_output_2d', 447 )
#endif
#endif

!
!--                For 2D-arrays (e.g. u*) only one cross-section is available.
!--                Hence exit loop of output levels.
                   IF ( two_d )  THEN
                      IF ( netcdf_data_format < 5 )  two_d = .FALSE.
                      EXIT loop1
                   ENDIF

                CASE ( 'xz' )
!
!--                Update the netCDF xz cross section time axis.
!--                In case of parallel output, this is only done by PE0
!--                to increase the performance.
                   IF ( simulated_time /= do2d_xz_last_time(av) )  THEN
                      do2d_xz_time_count(av) = do2d_xz_time_count(av) + 1
                      do2d_xz_last_time(av)  = simulated_time
                      IF ( myid == 0 )  THEN
                         IF ( .NOT. data_output_2d_on_each_pe  &
                              .OR.  netcdf_data_format > 4 )   &
                         THEN
#if defined( __netcdf )
                            nc_stat = NF90_PUT_VAR( id_set_xz(av),             &
                                                    id_var_time_xz(av),        &
                                             (/ time_since_reference_point /), &
                                         start = (/ do2d_xz_time_count(av) /), &
                                                    count = (/ 1 /) )
                            CALL netcdf_handle_error( 'data_output_2d', 56 )
#endif
                         ENDIF
                      ENDIF
                   ENDIF

!
!--                If required, carry out averaging along y
                   IF ( section(is,s_ind) == -1 )  THEN

                      ALLOCATE( local_2d_l(nxl:nxr,nzb_do:nzt_do) )
                      local_2d_l = 0.0_wp
                      ngp = ( nxr-nxl + 1 ) * ( nzt_do-nzb_do + 1 )
!
!--                   First local averaging on the PE
                      DO  k = nzb_do, nzt_do
                         DO  j = nys, nyn
                            DO  i = nxl, nxr
                               local_2d_l(i,k) = local_2d_l(i,k) +             &
                                                 local_pf(i,j,k)
                            ENDDO
                         ENDDO
                      ENDDO
#if defined( __parallel )
!
!--                   Now do the averaging over all PEs along y
                      IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
                      CALL MPI_ALLREDUCE( local_2d_l(nxl,nzb_do),                &
                                          local_2d(nxl,nzb_do), ngp, MPI_REAL,   &
                                          MPI_SUM, comm1dy, ierr )
#else
                      local_2d = local_2d_l
#endif
                      local_2d = local_2d / ( ny + 1.0_wp )

                      DEALLOCATE( local_2d_l )

                   ELSE
!
!--                   Just store the respective section on the local array
!--                   (but only if it is available on this PE!)
                      IF ( section(is,s_ind) >= nys  .AND.  section(is,s_ind) <= nyn ) &
                      THEN
                         local_2d = local_pf(:,section(is,s_ind),nzb_do:nzt_do)
                      ENDIF

                   ENDIF

#if defined( __parallel )
                   IF ( netcdf_data_format > 4 )  THEN
!
!--                   Output in netCDF4/HDF5 format.
!--                   Output only on those PEs where the respective cross
!--                   sections reside. Cross sections averaged along y are
!--                   output on the respective first PE along y (myidy=0).
                      IF ( ( section(is,s_ind) >= nys  .AND.                   &
                             section(is,s_ind) <= nyn )  .OR.                  &
                           ( section(is,s_ind) == -1  .AND.  myidy == 0 ) )  THEN
#if defined( __netcdf )
!
!--                      For parallel output, all cross sections are first
!--                      stored here on a local array and will be written to the
!--                      output file afterwards to increase the performance.
                         DO  i = nxl, nxr
                            DO  k = nzb_do, nzt_do
                               local_2d_sections_l(i,is,k) = local_2d(i,k)
                            ENDDO
                         ENDDO
#endif
                      ENDIF

                   ELSE

                      IF ( data_output_2d_on_each_pe )  THEN
!
!--                      Output of partial arrays on each PE. If the cross
!--                      section does not reside on the PE, output special
!--                      index values.
#if defined( __netcdf )
                         IF ( myid == 0 )  THEN
                            WRITE ( 22 )  time_since_reference_point,          &
                                          do2d_xz_time_count(av), av
                         ENDIF
#endif
                         DO  i = 0, io_blocks-1
                            IF ( i == io_group )  THEN
                               IF ( ( section(is,s_ind) >= nys  .AND.          &
                                      section(is,s_ind) <= nyn )  .OR.         &
                                    ( section(is,s_ind) == -1  .AND.           &
                                      nys-1 == -1 ) )                          &
                               THEN
                                  WRITE (22)  nxl, nxr, nzb_do, nzt_do, nzb, nzt+1
                                  WRITE (22)  local_2d
                               ELSE
                                  WRITE (22)  -1, -1, -1, -1, -1, -1
                               ENDIF
                            ENDIF
#if defined( __parallel )
                            CALL MPI_BARRIER( comm2d, ierr )
#endif
                         ENDDO

                      ELSE
!
!--                      PE0 receives partial arrays from all processors of the 
!--                      respective cross section and outputs them. Here a 
!--                      barrier has to be set, because otherwise 
!--                      "-MPI- FATAL: Remote protocol queue full" may occur.
                         CALL MPI_BARRIER( comm2d, ierr )

                         ngp = ( nxr-nxl + 1 ) * ( nzt_do-nzb_do + 1 )
                         IF ( myid == 0 )  THEN
!
!--                         Local array can be relocated directly.
                            IF ( ( section(is,s_ind) >= nys  .AND.              &
                                   section(is,s_ind) <= nyn )  .OR.             &
                                 ( section(is,s_ind) == -1  .AND.               &
                                   nys-1 == -1 ) )  THEN
                               total_2d(nxl:nxr,nzb_do:nzt_do) = local_2d
                            ENDIF
!
!--                         Receive data from all other PEs.
                            DO  n = 1, numprocs-1
!
!--                            Receive index limits first, then array.
!--                            Index limits are received in arbitrary order from
!--                            the PEs.
                               CALL MPI_RECV( ind(1), 4, MPI_INTEGER,          &
                                              MPI_ANY_SOURCE, 0, comm2d,       &
                                              status, ierr )
!
!--                            Not all PEs have data for XZ-cross-section.
                               IF ( ind(1) /= -9999 )  THEN
                                  sender = status(MPI_SOURCE)
                                  DEALLOCATE( local_2d )
                                  ALLOCATE( local_2d(ind(1):ind(2),            &
                                                     ind(3):ind(4)) )
                                  CALL MPI_RECV( local_2d(ind(1),ind(3)), ngp, &
                                                 MPI_REAL, sender, 1, comm2d,  &
                                                 status, ierr )
                                  total_2d(ind(1):ind(2),ind(3):ind(4)) =      &
                                                                        local_2d
                               ENDIF
                            ENDDO
!
!--                         Relocate the local array for the next loop increment
                            DEALLOCATE( local_2d )
                            ALLOCATE( local_2d(nxl:nxr,nzb_do:nzt_do) )

#if defined( __netcdf )
                            nc_stat = NF90_PUT_VAR( id_set_xz(av),             &
                                                 id_var_do2d(av,if),           &
                                                 total_2d(0:nx,nzb_do:nzt_do), &
                               start = (/ 1, is, 1, do2d_xz_time_count(av) /), &
                                          count = (/ nx+1, 1, nzt_do-nzb_do+1, 1 /) )
                            CALL netcdf_handle_error( 'data_output_2d', 58 )
#endif

                         ELSE
!
!--                         If the cross section resides on the PE, send the
!--                         local index limits, otherwise send -9999 to PE0.
                            IF ( ( section(is,s_ind) >= nys  .AND.              &
                                   section(is,s_ind) <= nyn )  .OR.             &
                                 ( section(is,s_ind) == -1  .AND.  nys-1 == -1 ) ) &
                            THEN
                               ind(1) = nxl; ind(2) = nxr
                               ind(3) = nzb_do;   ind(4) = nzt_do
                            ELSE
                               ind(1) = -9999; ind(2) = -9999
                               ind(3) = -9999; ind(4) = -9999
                            ENDIF
                            CALL MPI_SEND( ind(1), 4, MPI_INTEGER, 0, 0,       &
                                           comm2d, ierr )
!
!--                         If applicable, send data to PE0.
                            IF ( ind(1) /= -9999 )  THEN
                               CALL MPI_SEND( local_2d(nxl,nzb_do), ngp,         &
                                              MPI_REAL, 0, 1, comm2d, ierr )
                            ENDIF
                         ENDIF
!
!--                      A barrier has to be set, because otherwise some PEs may
!--                      proceed too fast so that PE0 may receive wrong data on
!--                      tag 0
                         CALL MPI_BARRIER( comm2d, ierr )
                      ENDIF

                   ENDIF
#else
#if defined( __netcdf )
                   nc_stat = NF90_PUT_VAR( id_set_xz(av),                   &
                                           id_var_do2d(av,if),              &
                                           local_2d(nxl:nxr,nzb_do:nzt_do), &
                            start = (/ 1, is, 1, do2d_xz_time_count(av) /), &
                                       count = (/ nx+1, 1, nzt_do-nzb_do+1, 1 /) )
                   CALL netcdf_handle_error( 'data_output_2d', 451 )
#endif
#endif

                CASE ( 'yz' )
!
!--                Update the netCDF yz cross section time axis.
!--                In case of parallel output, this is only done by PE0
!--                to increase the performance.
                   IF ( simulated_time /= do2d_yz_last_time(av) )  THEN
                      do2d_yz_time_count(av) = do2d_yz_time_count(av) + 1
                      do2d_yz_last_time(av)  = simulated_time
                      IF ( myid == 0 )  THEN
                         IF ( .NOT. data_output_2d_on_each_pe  &
                              .OR.  netcdf_data_format > 4 )   &
                         THEN
#if defined( __netcdf )
                            nc_stat = NF90_PUT_VAR( id_set_yz(av),             &
                                                    id_var_time_yz(av),        &
                                             (/ time_since_reference_point /), &
                                         start = (/ do2d_yz_time_count(av) /), &
                                                    count = (/ 1 /) )
                            CALL netcdf_handle_error( 'data_output_2d', 59 )
#endif
                         ENDIF
                      ENDIF
                   ENDIF

!
!--                If required, carry out averaging along x
                   IF ( section(is,s_ind) == -1 )  THEN

                      ALLOCATE( local_2d_l(nys:nyn,nzb_do:nzt_do) )
                      local_2d_l = 0.0_wp
                      ngp = ( nyn-nys+1 ) * ( nzt_do-nzb_do+1 )
!
!--                   First local averaging on the PE
                      DO  k = nzb_do, nzt_do
                         DO  j = nys, nyn
                            DO  i = nxl, nxr
                               local_2d_l(j,k) = local_2d_l(j,k) +             &
                                                 local_pf(i,j,k)
                            ENDDO
                         ENDDO
                      ENDDO
#if defined( __parallel )
!
!--                   Now do the averaging over all PEs along x
                      IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
                      CALL MPI_ALLREDUCE( local_2d_l(nys,nzb_do),                &
                                          local_2d(nys,nzb_do), ngp, MPI_REAL,   &
                                          MPI_SUM, comm1dx, ierr )
#else
                      local_2d = local_2d_l
#endif
                      local_2d = local_2d / ( nx + 1.0_wp )

                      DEALLOCATE( local_2d_l )

                   ELSE
!
!--                   Just store the respective section on the local array
!--                   (but only if it is available on this PE!)
                      IF ( section(is,s_ind) >= nxl  .AND.  section(is,s_ind) <= nxr ) &
                      THEN
                         local_2d = local_pf(section(is,s_ind),:,nzb_do:nzt_do)
                      ENDIF

                   ENDIF

#if defined( __parallel )
                   IF ( netcdf_data_format > 4 )  THEN
!
!--                   Output in netCDF4/HDF5 format.
!--                   Output only on those PEs where the respective cross
!--                   sections reside. Cross sections averaged along x are
!--                   output on the respective first PE along x (myidx=0).
                      IF ( ( section(is,s_ind) >= nxl  .AND.                       &
                             section(is,s_ind) <= nxr )  .OR.                      &
                           ( section(is,s_ind) == -1  .AND.  myidx == 0 ) )  THEN
#if defined( __netcdf )
!
!--                      For parallel output, all cross sections are first
!--                      stored here on a local array and will be written to the
!--                      output file afterwards to increase the performance.
                         DO  j = nys, nyn
                            DO  k = nzb_do, nzt_do
                               local_2d_sections_l(is,j,k) = local_2d(j,k)
                            ENDDO
                         ENDDO
#endif
                      ENDIF

                   ELSE

                      IF ( data_output_2d_on_each_pe )  THEN
!
!--                      Output of partial arrays on each PE. If the cross
!--                      section does not reside on the PE, output special
!--                      index values.
#if defined( __netcdf )
                         IF ( myid == 0 )  THEN
                            WRITE ( 23 )  time_since_reference_point,          &
                                          do2d_yz_time_count(av), av
                         ENDIF
#endif
                         DO  i = 0, io_blocks-1
                            IF ( i == io_group )  THEN
                               IF ( ( section(is,s_ind) >= nxl  .AND.          &
                                      section(is,s_ind) <= nxr )  .OR.         &
                                    ( section(is,s_ind) == -1  .AND.           &
                                      nxl-1 == -1 ) )                          &
                               THEN
                                  WRITE (23)  nys, nyn, nzb_do, nzt_do, nzb, nzt+1
                                  WRITE (23)  local_2d
                               ELSE
                                  WRITE (23)  -1, -1, -1, -1, -1, -1
                               ENDIF
                            ENDIF
#if defined( __parallel )
                            CALL MPI_BARRIER( comm2d, ierr )
#endif
                         ENDDO

                      ELSE
!
!--                      PE0 receives partial arrays from all processors of the 
!--                      respective cross section and outputs them. Here a
!--                      barrier has to be set, because otherwise 
!--                      "-MPI- FATAL: Remote protocol queue full" may occur.
                         CALL MPI_BARRIER( comm2d, ierr )

                         ngp = ( nyn-nys+1 ) * ( nzt_do-nzb_do+1 )
                         IF ( myid == 0 )  THEN
!
!--                         Local array can be relocated directly.
                            IF ( ( section(is,s_ind) >= nxl  .AND.             &
                                   section(is,s_ind) <= nxr )   .OR.           &
                                 ( section(is,s_ind) == -1  .AND.  nxl-1 == -1 ) ) &
                            THEN
                               total_2d(nys:nyn,nzb_do:nzt_do) = local_2d
                            ENDIF
!
!--                         Receive data from all other PEs.
                            DO  n = 1, numprocs-1
!
!--                            Receive index limits first, then array.
!--                            Index limits are received in arbitrary order from
!--                            the PEs.
                               CALL MPI_RECV( ind(1), 4, MPI_INTEGER,          &
                                              MPI_ANY_SOURCE, 0, comm2d,       &
                                              status, ierr )
!
!--                            Not all PEs have data for YZ-cross-section.
                               IF ( ind(1) /= -9999 )  THEN
                                  sender = status(MPI_SOURCE)
                                  DEALLOCATE( local_2d )
                                  ALLOCATE( local_2d(ind(1):ind(2),            &
                                                     ind(3):ind(4)) )
                                  CALL MPI_RECV( local_2d(ind(1),ind(3)), ngp, &
                                                 MPI_REAL, sender, 1, comm2d,  &
                                                 status, ierr )
                                  total_2d(ind(1):ind(2),ind(3):ind(4)) =      &
                                                                        local_2d
                               ENDIF
                            ENDDO
!
!--                         Relocate the local array for the next loop increment
                            DEALLOCATE( local_2d )
                            ALLOCATE( local_2d(nys:nyn,nzb_do:nzt_do) )

#if defined( __netcdf )
                            nc_stat = NF90_PUT_VAR( id_set_yz(av),             &
                                                 id_var_do2d(av,if),           &
                                                 total_2d(0:ny,nzb_do:nzt_do), &
                            start = (/ is, 1, 1, do2d_yz_time_count(av) /),    &
                                       count = (/ 1, ny+1, nzt_do-nzb_do+1, 1 /) )
                            CALL netcdf_handle_error( 'data_output_2d', 61 )
#endif

                         ELSE
!
!--                         If the cross section resides on the PE, send the
!--                         local index limits, otherwise send -9999 to PE0.
                            IF ( ( section(is,s_ind) >= nxl  .AND.              &
                                   section(is,s_ind) <= nxr )  .OR.             &
                                 ( section(is,s_ind) == -1  .AND.  nxl-1 == -1 ) ) &
                            THEN
                               ind(1) = nys; ind(2) = nyn
                               ind(3) = nzb_do;   ind(4) = nzt_do
                            ELSE
                               ind(1) = -9999; ind(2) = -9999
                               ind(3) = -9999; ind(4) = -9999
                            ENDIF
                            CALL MPI_SEND( ind(1), 4, MPI_INTEGER, 0, 0,       &
                                           comm2d, ierr )
!
!--                         If applicable, send data to PE0.
                            IF ( ind(1) /= -9999 )  THEN
                               CALL MPI_SEND( local_2d(nys,nzb_do), ngp,         &
                                              MPI_REAL, 0, 1, comm2d, ierr )
                            ENDIF
                         ENDIF
!
!--                      A barrier has to be set, because otherwise some PEs may
!--                      proceed too fast so that PE0 may receive wrong data on
!--                      tag 0
                         CALL MPI_BARRIER( comm2d, ierr )
                      ENDIF

                   ENDIF
#else
#if defined( __netcdf )
                   nc_stat = NF90_PUT_VAR( id_set_yz(av),                   &
                                           id_var_do2d(av,if),              &
                                           local_2d(nys:nyn,nzb_do:nzt_do), &
                            start = (/ is, 1, 1, do2d_xz_time_count(av) /), &
                                           count = (/ 1, ny+1, nzt_do-nzb_do+1, 1 /) )
                   CALL netcdf_handle_error( 'data_output_2d', 452 )
#endif
#endif

             END SELECT

             is = is + 1
          ENDDO loop1

!
!--       For parallel output, all data were collected before on a local array
!--       and are written now to the netcdf file. This must be done to increase
!--       the performance of the parallel output.
#if defined( __netcdf )
          IF ( netcdf_data_format > 4 )  THEN

                SELECT CASE ( mode )

                   CASE ( 'xy' )
                      IF ( two_d ) THEN
                         nis = 1
                         two_d = .FALSE.
                      ELSE
                         nis = ns
                      ENDIF
!
!--                   Do not output redundant ghost point data except for the
!--                   boundaries of the total domain.
!                      IF ( nxr == nx  .AND.  nyn /= ny )  THEN
!                         nc_stat = NF90_PUT_VAR( id_set_xy(av),                &
!                                                 id_var_do2d(av,if),           &
!                                                 local_2d_sections(nxl:nxr+1,  &
!                                                    nys:nyn,1:nis),            &
!                                                 start = (/ nxl+1, nys+1, 1,   &
!                                                    do2d_xy_time_count(av) /), &
!                                                 count = (/ nxr-nxl+2,         &
!                                                            nyn-nys+1, nis, 1  &
!                                                          /) )
!                      ELSEIF ( nxr /= nx  .AND.  nyn == ny )  THEN
!                         nc_stat = NF90_PUT_VAR( id_set_xy(av),                &
!                                                 id_var_do2d(av,if),           &
!                                                 local_2d_sections(nxl:nxr,    &
!                                                    nys:nyn+1,1:nis),          &
!                                                 start = (/ nxl+1, nys+1, 1,   &
!                                                    do2d_xy_time_count(av) /), &
!                                                 count = (/ nxr-nxl+1,         &
!                                                            nyn-nys+2, nis, 1  &
!                                                          /) )
!                      ELSEIF ( nxr == nx  .AND.  nyn == ny )  THEN
!                         nc_stat = NF90_PUT_VAR( id_set_xy(av),                &
!                                                 id_var_do2d(av,if),           &
!                                                 local_2d_sections(nxl:nxr+1,  &
!                                                    nys:nyn+1,1:nis),          &
!                                                 start = (/ nxl+1, nys+1, 1,   &
!                                                    do2d_xy_time_count(av) /), &
!                                                 count = (/ nxr-nxl+2,         &
!                                                            nyn-nys+2, nis, 1  &
!                                                          /) )
!                      ELSE
                         nc_stat = NF90_PUT_VAR( id_set_xy(av),                &
                                                 id_var_do2d(av,if),           &
                                                 local_2d_sections(nxl:nxr,    &
                                                    nys:nyn,1:nis),            &
                                                 start = (/ nxl+1, nys+1, 1,   &
                                                    do2d_xy_time_count(av) /), &
                                                 count = (/ nxr-nxl+1,         &
                                                            nyn-nys+1, nis, 1  &
                                                          /) )
!                      ENDIF    

                      CALL netcdf_handle_error( 'data_output_2d', 55 )

                   CASE ( 'xz' )
!
!--                   First, all PEs get the information of all cross-sections.
!--                   Then the data are written to the output file by all PEs
!--                   while NF90_COLLECTIVE is set in subroutine 
!--                   define_netcdf_header. Although redundant information are
!--                   written to the output file in that case, the performance 
!--                   is significantly better compared to the case where only
!--                   the first row of PEs in x-direction (myidx = 0) is given
!--                   the output while NF90_INDEPENDENT is set.
                      IF ( npey /= 1 )  THEN
                         
#if defined( __parallel )
!
!--                      Distribute data over all PEs along y
                         ngp = ( nxr-nxl+1 ) * ( nzt_do-nzb_do+1 ) * ns
                         IF ( collective_wait ) CALL MPI_BARRIER( comm2d, ierr )
                         CALL MPI_ALLREDUCE( local_2d_sections_l(nxl,1,nzb_do),  &
                                             local_2d_sections(nxl,1,nzb_do),    &
                                             ngp, MPI_REAL, MPI_SUM, comm1dy,  &
                                             ierr )
#else
                         local_2d_sections = local_2d_sections_l
#endif
                      ENDIF
!
!--                   Do not output redundant ghost point data except for the
!--                   boundaries of the total domain.
!                      IF ( nxr == nx )  THEN
!                         nc_stat = NF90_PUT_VAR( id_set_xz(av),                &
!                                             id_var_do2d(av,if),               & 
!                                             local_2d_sections(nxl:nxr+1,1:ns, &
!                                                nzb_do:nzt_do),                &
!                                             start = (/ nxl+1, 1, 1,           &
!                                                do2d_xz_time_count(av) /),     &
!                                             count = (/ nxr-nxl+2, ns, nzt_do-nzb_do+1,  &
!                                                        1 /) )
!                      ELSE
                         nc_stat = NF90_PUT_VAR( id_set_xz(av),                &
                                             id_var_do2d(av,if),               &
                                             local_2d_sections(nxl:nxr,1:ns,   &
                                                nzb_do:nzt_do),                &
                                             start = (/ nxl+1, 1, 1,           &
                                                do2d_xz_time_count(av) /),     &
                                             count = (/ nxr-nxl+1, ns, nzt_do-nzb_do+1,  &
                                                1 /) )
!                      ENDIF

                      CALL netcdf_handle_error( 'data_output_2d', 57 )

                   CASE ( 'yz' )
!
!--                   First, all PEs get the information of all cross-sections.
!--                   Then the data are written to the output file by all PEs
!--                   while NF90_COLLECTIVE is set in subroutine 
!--                   define_netcdf_header. Although redundant information are
!--                   written to the output file in that case, the performance 
!--                   is significantly better compared to the case where only
!--                   the first row of PEs in y-direction (myidy = 0) is given
!--                   the output while NF90_INDEPENDENT is set.
                      IF ( npex /= 1 )  THEN

#if defined( __parallel )
!
!--                      Distribute data over all PEs along x
                         ngp = ( nyn-nys+1 ) * ( nzt-nzb + 2 ) * ns
                         IF ( collective_wait ) CALL MPI_BARRIER( comm2d, ierr )
                         CALL MPI_ALLREDUCE( local_2d_sections_l(1,nys,nzb_do),  &
                                             local_2d_sections(1,nys,nzb_do),    &
                                             ngp, MPI_REAL, MPI_SUM, comm1dx,  &
                                             ierr )
#else
                         local_2d_sections = local_2d_sections_l
#endif
                      ENDIF
!
!--                   Do not output redundant ghost point data except for the
!--                   boundaries of the total domain.
!                      IF ( nyn == ny )  THEN
!                         nc_stat = NF90_PUT_VAR( id_set_yz(av),                &
!                                             id_var_do2d(av,if),               &
!                                             local_2d_sections(1:ns,           &
!                                                nys:nyn+1,nzb_do:nzt_do),      &
!                                             start = (/ 1, nys+1, 1,           &
!                                                do2d_yz_time_count(av) /),     &
!                                             count = (/ ns, nyn-nys+2,         &
!                                                        nzt_do-nzb_do+1, 1 /) )
!                      ELSE
                         nc_stat = NF90_PUT_VAR( id_set_yz(av),                &
                                             id_var_do2d(av,if),               &
                                             local_2d_sections(1:ns,nys:nyn,   &
                                                nzb_do:nzt_do),                &
                                             start = (/ 1, nys+1, 1,           &
                                                do2d_yz_time_count(av) /),     &
                                             count = (/ ns, nyn-nys+1,         &
                                                        nzt_do-nzb_do+1, 1 /) )
!                      ENDIF

                      CALL netcdf_handle_error( 'data_output_2d', 60 )

                   CASE DEFAULT
                      message_string = 'unknown cross-section: ' // TRIM( mode )
                      CALL message( 'data_output_2d', 'PA0180', 1, 2, 0, 6, 0 )

                END SELECT                     

          ENDIF
#endif
       ENDIF

       if = if + 1
       l = MAX( 2, LEN_TRIM( do2d(av,if) ) )
       do2d_mode = do2d(av,if)(l-1:l)

    ENDDO

!
!-- Deallocate temporary arrays.
    IF ( ALLOCATED( level_z ) )  DEALLOCATE( level_z )
    IF ( netcdf_data_format > 4 )  THEN
       DEALLOCATE( local_pf, local_2d, local_2d_sections )
       IF( mode == 'xz' .OR. mode == 'yz' ) DEALLOCATE( local_2d_sections_l )
    ENDIF
#if defined( __parallel )
    IF ( .NOT.  data_output_2d_on_each_pe  .AND.  myid == 0 )  THEN
       DEALLOCATE( total_2d )
    ENDIF
#endif

!
!-- Close plot output file.
    file_id = 20 + s_ind

    IF ( data_output_2d_on_each_pe )  THEN
       DO  i = 0, io_blocks-1
          IF ( i == io_group )  THEN
             CALL close_file( file_id )
          ENDIF
#if defined( __parallel )
          CALL MPI_BARRIER( comm2d, ierr )
#endif
       ENDDO
    ELSE
       IF ( myid == 0 )  CALL close_file( file_id )
    ENDIF

    CALL cpu_log( log_point(3), 'data_output_2d', 'stop' )

 END SUBROUTINE data_output_2d
