!> @file data_output_3d.f90
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
! ------------------
! 
! 
! Former revisions:
! -----------------
! $Id: data_output_3d.f90 3049 2018-05-29 13:52:36Z Giersch $
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error message revised
! 
! 3014 2018-05-09 08:42:38Z maronga
! Added nzb_do and nzt_do for some modules for 3d output
! 
! 3004 2018-04-27 12:33:25Z Giersch
! Allocation checks implemented (averaged data will be assigned to fill values 
! if no allocation happened so far) 
! 
! 2967 2018-04-13 11:22:08Z raasch
! bugfix: missing parallel cpp-directives added
! 
! 2817 2018-02-19 16:32:21Z knoop
! Preliminary gust module interface implemented
! 
! 2766 2018-01-22 17:17:47Z kanani
! Removed preprocessor directive __chem
! 
! 2756 2018-01-16 18:11:14Z suehring
! Fill values for 3D output of chemical species introduced.
! 
! 2746 2018-01-15 12:06:04Z suehring
! Move flag plant canopy to modules
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Implementation of turbulence_closure_mod (TG)
! Implementation of chemistry module (FK)
! Set fill values at topography grid points or e.g. non-natural-type surface
! in case of LSM output (MS)
! 
! 2512 2017-10-04 08:26:59Z raasch
! upper bounds of 3d output changed from nx+1,ny+1 to nx,ny
! no output of ghost layer data
! 
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison' 
! includes two more prognostic equations for cloud drop concentration (nc)  
! and cloud water content (qc). 
! 
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography concept
! 
! 2209 2017-04-19 09:34:46Z kanani
! Added plant canopy model output
! 
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean and rho_av to rho_ocean_av
! 
! 2011 2016-09-19 17:29:57Z kanani
! Flag urban_surface is now defined in module control_parameters,
! changed prefix for urban surface model output to "usm_",
! introduced control parameter varnamelength for LEN of trimvar.
! 
! 2007 2016-08-24 15:47:17Z kanani
! Added support for new urban surface model (temporary modifications of 
! SELECT CASE ( ) necessary, see variable trimvar)
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
! Output of land surface quantities is now done directly in the respective module.
! Unnecessary directive __parallel removed.
! 
! 1960 2016-07-12 16:34:24Z suehring
! Scalar surface flux added
!
! 1849 2016-04-08 11:33:18Z hoffmann
! prr moved to arrays_3d
!
! 1822 2016-04-07 07:49:42Z hoffmann
! prr vertical dimensions set to nzb_do to nzt_do. Unused variables deleted.
!
! 1808 2016-04-05 19:44:00Z raasch
! test output removed
!
! 1783 2016-03-06 18:36:17Z raasch
! name change of netcdf routines and module + related changes
!
! 1745 2016-02-05 13:06:51Z gronemeier
! Bugfix: test if time axis limit exceeds moved to point after call of check_open
! 
! 1691 2015-10-26 16:17:44Z maronga
! Added output of radiative heating rates for RRTMG
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1585 2015-04-30 07:05:52Z maronga
! Added support for RRTMG
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
! parts concerning avs output removed,
! -netcdf output queries
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
! barrier argument removed from cpu_log,
! module interfaces removed
!
! 1308 2014-03-13 14:58:42Z fricke
! Check, if the limit of the time dimension is exceeded for parallel output
! To increase the performance for parallel output, the following is done:
! - Update of time axis is only done by PE0
!
! 1244 2013-10-31 08:16:56Z raasch
! Bugfix for index bounds in case of 3d-parallel output
!
! 1115 2013-03-26 18:16:16Z hoffmann
! ql is calculated by calc_liquid_water_content
!
! 1106 2013-03-04 05:31:38Z raasch
! array_kind renamed precision_kind
!
! 1076 2012-12-05 08:30:18Z hoffmann
! Bugfix in output of ql
!
! 1053 2012-11-13 17:11:03Z hoffmann
! +nr, qr, prr, qc and averaged quantities
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
! Revision 1.1  1997/09/03 06:29:36  raasch
! Initial revision
!
!
! Description:
! ------------
!> Output of the 3D-arrays in netCDF and/or AVS format.
!------------------------------------------------------------------------------!
 SUBROUTINE data_output_3d( av )
 

    USE arrays_3d,                                                             &
        ONLY:  e, nc, nr, p, pt, prr, q, qc, ql, ql_c, ql_v, qr, rho_ocean, s, &
               sa, tend, u, v, vpt, w, alpha_T, beta_S, solar3d
        
    USE averaging
        
    USE control_parameters,                                                    &
        ONLY:  air_chemistry, cloud_physics, do3d, do3d_no, do3d_time_count,   &
               io_blocks, io_group, land_surface, message_string,              &
               ntdim_3d, nz_do3d,  plant_canopy,                               &
               psolver, simulated_time, time_since_reference_point,            &
               urban_surface, varnamelength
        
    USE cpulog,                                                                &
        ONLY:  log_point, cpu_log

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg, nzb,  &
               nzt, wall_flags_0
        
    USE kinds
    
#if defined( __netcdf )
    USE NETCDF
#endif

    USE netcdf_interface,                                                      &
        ONLY:  fill_value, id_set_3d, id_var_do3d, id_var_time_3d, nc_stat,    &
               netcdf_data_format, netcdf_handle_error
        
    USE pegrid
    USE turbulence_closure_mod,                                                &
        ONLY:  tcm_data_output_3d

    IMPLICIT NONE

    INTEGER(iwp) ::  av        !< 
    INTEGER(iwp) ::  flag_nr   !< number of masking flag
    INTEGER(iwp) ::  i         !< 
    INTEGER(iwp) ::  if        !< 
    INTEGER(iwp) ::  j         !< 
    INTEGER(iwp) ::  k         !< 
    INTEGER(iwp) ::  n         !< 
    INTEGER(iwp) ::  nzb_do    !< vertical lower limit for data output
    INTEGER(iwp) ::  nzt_do    !< vertical upper limit for data output

    LOGICAL      ::  found     !< 
    LOGICAL      ::  resorted  !< 

    REAL(wp)     ::  mean_r    !< 
    REAL(wp)     ::  s_r2      !< 
    REAL(wp)     ::  s_r3      !< 

    REAL(sp), DIMENSION(:,:,:), ALLOCATABLE ::  local_pf  !<

    REAL(wp), DIMENSION(:,:,:), POINTER ::  to_be_resorted  !< 

    CHARACTER (LEN=varnamelength) ::  trimvar  !< TRIM of output-variable string

!
!-- Return, if nothing to output
    IF ( do3d_no(av) == 0 )  RETURN

    CALL cpu_log (log_point(14),'data_output_3d','start')

!
!-- Open output file.
!-- For classic or 64bit netCDF output on more than one PE, each PE opens its
!-- own file and writes the data of its subdomain in binary format. After the
!-- run, these files are combined to one NetCDF file by combine_plot_fields.
!-- For netCDF4/HDF5 output, data is written in parallel into one file.
    IF ( netcdf_data_format < 5 )  THEN
#if defined( __parallel )
       CALL check_open( 30 )
#endif
       IF ( myid == 0 )  CALL check_open( 106+av*10 )
    ELSE
       CALL check_open( 106+av*10 )
    ENDIF

!
!-- For parallel netcdf output the time axis must be limited. Return, if this
!-- limit is exceeded. This could be the case, if the simulated time exceeds 
!-- the given end time by the length of the given output interval.
    IF ( netcdf_data_format > 4 )  THEN
       IF ( do3d_time_count(av) + 1 > ntdim_3d(av) )  THEN
          WRITE ( message_string, * ) 'Output of 3d data is not given at t=',  &
                                      simulated_time, '&because the maximum ', & 
                                      'number of output time levels is ',      &
                                      'exceeded.'
          CALL message( 'data_output_3d', 'PA0387', 0, 1, 0, 6, 0 )
          CALL cpu_log( log_point(14), 'data_output_3d', 'stop' )
          RETURN
       ENDIF
    ENDIF

!
!-- Update the netCDF time axis
!-- In case of parallel output, this is only done by PE0 to increase the
!-- performance.
#if defined( __netcdf )
    do3d_time_count(av) = do3d_time_count(av) + 1
    IF ( myid == 0 )  THEN
       nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_time_3d(av),           &
                               (/ time_since_reference_point /),            &
                               start = (/ do3d_time_count(av) /),           &
                               count = (/ 1 /) )
       CALL netcdf_handle_error( 'data_output_3d', 376 )
    ENDIF
#endif

!
!-- Loop over all variables to be written.
    if = 1

    DO  WHILE ( do3d(av,if)(1:1) /= ' ' )

!
!--    Temporary solution to account for data output within the new urban 
!--    surface model (urban_surface_mod.f90), see also SELECT CASE ( trimvar ).
!--    Store the array chosen on the temporary array.
       trimvar = TRIM( do3d(av,if) )
         resorted = .FALSE.
          nzb_do   = nzb
          nzt_do   = nz_do3d
!
!--    Set flag to steer output of radiation, land-surface, or user-defined
!--    quantities
       found = .FALSE.
!
!--    Allocate a temporary array with the desired output dimensions.
       ALLOCATE( local_pf(nxl:nxr,nys:nyn,nzb_do:nzt_do) )
!
!--    Before each output, set array local_pf to fill value
       local_pf = fill_value
!
!--    Set masking flag for topography for not resorted arrays
       flag_nr = 0

       SELECT CASE ( trimvar )

          CASE ( 'e' )
             IF ( av == 0 )  THEN
                to_be_resorted => e
             ELSE
                IF ( .NOT. ALLOCATED( e_av ) ) THEN
                   ALLOCATE( e_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   e_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => e_av
             ENDIF

          CASE ( 'lpt' )
             IF ( av == 0 )  THEN
                to_be_resorted => pt
             ELSE
                IF ( .NOT. ALLOCATED( lpt_av ) ) THEN
                   ALLOCATE( lpt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   lpt_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => lpt_av
             ENDIF

          CASE ( 'nc' )
             IF ( av == 0 )  THEN
                to_be_resorted => nc
             ELSE
                IF ( .NOT. ALLOCATED( nc_av ) ) THEN
                   ALLOCATE( nc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   nc_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => nc_av
             ENDIF

          CASE ( 'nr' )
             IF ( av == 0 )  THEN
                to_be_resorted => nr
             ELSE
                IF ( .NOT. ALLOCATED( nr_av ) ) THEN
                   ALLOCATE( nr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   nr_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => nr_av
             ENDIF

          CASE ( 'p' )
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
          CASE ( 'prr' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = prr(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                IF ( .NOT. ALLOCATED( prr_av ) ) THEN
                   ALLOCATE( prr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   prr_av = REAL( fill_value, KIND = wp )
                ENDIF
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = prr_av(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
             resorted = .TRUE.

          CASE ( 'pt' )
             IF ( av == 0 )  THEN
                to_be_resorted => pt
             ELSE
                IF ( .NOT. ALLOCATED( pt_av ) ) THEN
                   ALLOCATE( pt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   pt_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => pt_av
             ENDIF

          CASE ( 'q' )
             IF ( av == 0 )  THEN
                to_be_resorted => q
             ELSE
                IF ( .NOT. ALLOCATED( q_av ) ) THEN
                   ALLOCATE( q_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   q_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => q_av
             ENDIF

          CASE ( 'qc' )
             IF ( av == 0 )  THEN
                to_be_resorted => qc
             ELSE
                IF ( .NOT. ALLOCATED( qc_av ) ) THEN
                   ALLOCATE( qc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   qc_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => qc_av
             ENDIF

          CASE ( 'ql' )
             IF ( av == 0 )  THEN
                to_be_resorted => ql
             ELSE
                IF ( .NOT. ALLOCATED( ql_av ) ) THEN
                   ALLOCATE( ql_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ql_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => ql_av
             ENDIF

          CASE ( 'ql_c' )
             IF ( av == 0 )  THEN
                to_be_resorted => ql_c
             ELSE
                IF ( .NOT. ALLOCATED( ql_c_av ) ) THEN
                   ALLOCATE( ql_c_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ql_c_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => ql_c_av
             ENDIF

          CASE ( 'ql_v' )
             IF ( av == 0 )  THEN
                to_be_resorted => ql_v
             ELSE
                IF ( .NOT. ALLOCATED( ql_v_av ) ) THEN
                   ALLOCATE( ql_v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ql_v_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => ql_v_av
             ENDIF
          CASE ( 'qr' )
             IF ( av == 0 )  THEN
                to_be_resorted => qr
             ELSE
                IF ( .NOT. ALLOCATED( qr_av ) ) THEN
                   ALLOCATE( qr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   qr_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => qr_av
             ENDIF

          CASE ( 'qv' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
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

           CASE ( 'solar3d' )
             IF ( av == 0 )  THEN
                to_be_resorted => solar3d
             ELSE
                IF ( .NOT. ALLOCATED( solar3d_av ) ) THEN
                   ALLOCATE( solar3d_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   u_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => solar3d_av
             ENDIF


          CASE ( 'rho_ocean' )
             IF ( av == 0 )  THEN
                to_be_resorted => rho_ocean
             ELSE
                IF ( .NOT. ALLOCATED( rho_ocean_av ) ) THEN
                   ALLOCATE( rho_ocean_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   u_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => rho_ocean_av
             ENDIF

          CASE ( 'alpha_T' )
             IF ( av == 0 )  THEN
                to_be_resorted => alpha_T
             ELSE
                IF ( .NOT. ALLOCATED( alpha_T_av ) ) THEN
                   ALLOCATE( alpha_T_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   u_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => alpha_T_av
             ENDIF

          CASE ( 'beta_S' )
             IF ( av == 0 )  THEN
                to_be_resorted => beta_S
             ELSE
                IF ( .NOT. ALLOCATED( beta_S_av ) ) THEN
                   ALLOCATE( beta_S_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   u_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => beta_S_av
             ENDIF


          CASE ( 's' )
             IF ( av == 0 )  THEN
                to_be_resorted => s
             ELSE
                IF ( .NOT. ALLOCATED( s_av ) ) THEN
                   ALLOCATE( s_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   s_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => s_av
             ENDIF

          CASE ( 'sa' )
             IF ( av == 0 )  THEN
                to_be_resorted => sa
             ELSE
                IF ( .NOT. ALLOCATED( sa_av ) ) THEN
                   ALLOCATE( sa_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   sa_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => sa_av
             ENDIF

          CASE ( 'u' )
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

          CASE ( 'v' )
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

          CASE ( 'vpt' )
             IF ( av == 0 )  THEN
                to_be_resorted => vpt
             ELSE
                IF ( .NOT. ALLOCATED( vpt_av ) ) THEN
                   ALLOCATE( vpt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   vpt_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => vpt_av
             ENDIF

          CASE ( 'w' )
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

          CASE DEFAULT
             IF ( .NOT. found )  THEN
                CALL tcm_data_output_3d( av, do3d(av,if), found, local_pf,     &
                                         nzb_do, nzt_do )
                resorted = .TRUE.
             ENDIF
             IF ( .NOT. found )  THEN
                message_string =  'no output available for: ' //               &
                                  TRIM( do3d(av,if) )
                CALL message( 'data_output_3d', 'PA0182', 0, 0, 0, 6, 0 )
             ENDIF

       END SELECT

!
!--    Resort the array to be output, if not done above
       IF ( .NOT. resorted )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb_do, nzt_do
                   local_pf(i,j,k) = MERGE(                                    &
                                      to_be_resorted(k,j,i),                   &
                                      REAL( fill_value, KIND = wp ),           &
                                      BTEST( wall_flags_0(k,j,i), flag_nr ) )
                ENDDO
             ENDDO
          ENDDO
       ENDIF

!
!--    Output of the 3D-array
#if defined( __parallel )
       IF ( netcdf_data_format < 5 )  THEN
!
!--       Non-parallel netCDF output. Data is output in parallel in
!--       FORTRAN binary format here, and later collected into one file by
!--       combine_plot_fields
          IF ( myid == 0 )  THEN
             WRITE ( 30 )  time_since_reference_point,                   &
                           do3d_time_count(av), av
          ENDIF
          DO  i = 0, io_blocks-1
             IF ( i == io_group )  THEN
                WRITE ( 30 )  nxl, nxr, nys, nyn, nzb_do, nzt_do
                WRITE ( 30 )  local_pf(:,:,nzb_do:nzt_do)
             ENDIF

             CALL MPI_BARRIER( comm2d, ierr )

          ENDDO

       ELSE
#if defined( __netcdf )
!
!--       Parallel output in netCDF4/HDF5 format.
!          IF ( nxr == nx  .AND.  nyn /= ny )  THEN
!             nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_do3d(av,if),  &
!                               local_pf(nxl:nxr+1,nys:nyn,nzb_do:nzt_do),    &
!                start = (/ nxl+1, nys+1, nzb_do+1, do3d_time_count(av) /),  &
!                count = (/ nxr-nxl+2, nyn-nys+1, nzt_do-nzb_do+1, 1 /) )
!          ELSEIF ( nxr /= nx  .AND.  nyn == ny )  THEN
!             nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_do3d(av,if),  &
!                               local_pf(nxl:nxr,nys:nyn+1,nzb_do:nzt_do),    &
!                start = (/ nxl+1, nys+1, nzb_do+1, do3d_time_count(av) /),  &
!                count = (/ nxr-nxl+1, nyn-nys+2, nzt_do-nzb_do+1, 1 /) )
!          ELSEIF ( nxr == nx  .AND.  nyn == ny )  THEN
!             nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_do3d(av,if),  &
!                             local_pf(nxl:nxr+1,nys:nyn+1,nzb_do:nzt_do  ),  &
!                start = (/ nxl+1, nys+1, nzb_do+1, do3d_time_count(av) /),  &
!                count = (/ nxr-nxl+2, nyn-nys+2, nzt_do-nzb_do+1, 1 /) )
!          ELSE
             nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_do3d(av,if),  &
                                 local_pf(nxl:nxr,nys:nyn,nzb_do:nzt_do),    &
                start = (/ nxl+1, nys+1, nzb_do+1, do3d_time_count(av) /),  &
                count = (/ nxr-nxl+1, nyn-nys+1, nzt_do-nzb_do+1, 1 /) )
!          ENDIF
          CALL netcdf_handle_error( 'data_output_3d', 386 )
#endif
       ENDIF
#else
#if defined( __netcdf )
       nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_do3d(av,if),        &
                         local_pf(nxl:nxr+1,nys:nyn+1,nzb_do:nzt_do),        &
                         start = (/ 1, 1, 1, do3d_time_count(av) /),     &
                         count = (/ nx+1, ny+1, nzt_do-nzb_do+1, 1 /) )
       CALL netcdf_handle_error( 'data_output_3d', 446 )
#endif
#endif

       if = if + 1

!
!--    Deallocate temporary array
       DEALLOCATE ( local_pf )

    ENDDO

    CALL cpu_log( log_point(14), 'data_output_3d', 'stop' )

!
!-- Formats.
3300 FORMAT ('variable ',I4,'  file=',A,'  filetype=unformatted  skip=',I12/   &
             'label = ',A,A)

 END SUBROUTINE data_output_3d
