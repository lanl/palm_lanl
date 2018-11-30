!> @file netcdf_interface_mod.f90
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
! $Id: netcdf_interface_mod.f90 3049 2018-05-29 13:52:36Z Giersch $
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised, code adjusted to PALMs coding standards, CASE pt_ext
! pt_new disabled, comment revised
! 
! 3004 2018-04-27 12:33:25Z Giersch
! .NOT. found in if-query added to account for variables found in tcm
! 
! 2964 2018-04-12 16:04:03Z Giersch
! Calculation of fixed number of output time levels for parallel netcdf output
! has been moved completely to check_parameters
! 
! 2932 2018-03-26 09:39:22Z maronga
! Renamed inipar to initialization_parameters.
! 
! 2817 2018-02-19 16:32:21Z knoop
! Preliminary gust module interface implemented
! 
! 2769 2018-01-25 09:22:24Z raasch
! bugfix for calculating number of required output time levels in case of output
! at the beginning of a restart run
! 
! 2766 2018-01-22 17:17:47Z kanani
! Removed preprocessor directive __chem
! 
! 2746 2018-01-15 12:06:04Z suehring
! Move flag plant canopy to modules
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Implementation of uv exposure model (FK)
! Implemented checks for turbulence_closure_mod (TG)
! Implementation of chemistry module (FK)
! Bugfix in setting netcdf grids for LSM variables
! Enable setting of _FillValue attribute in output files (MS)
!
! 2512 2017-10-04 08:26:59Z raasch
! upper bounds of cross section and 3d output changed from nx+1,ny+1 to nx,ny
! no output of ghost layer data any more
! 
! 2302 2017-07-03 14:07:20Z suehring
! Reading of 3D topography using NetCDF data type NC_BYTE
! 
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison' 
! includes two more prognostic equations for cloud drop concentration (nc)  
! and cloud water content (qc). 
! 
! 2270 2017-06-09 12:18:47Z maronga
! Removed 2 timeseries (shf_eb + qsws_eb). Removed _eb suffixes
! 
! 2265 2017-06-08 16:58:28Z schwenkel
! Unused variables removed.
! 
! 2239 2017-06-01 12:04:51Z suehring
! Bugfix xy-output of land-surface variables
! 
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography and surface concept
! 
! Topograpyh height arrays (zu_s_inner, zw_w_inner) are defined locally, output
! only if parallel netcdf. 
!
! Build interface for topography input: 
! - open file in read-only mode
! - read global attributes
! - read variables
!
! Bugfix in xy output (land-surface case)
!
! 2209 2017-04-19 09:34:46Z kanani
! Added support for plant canopy model output
! 
! 2189 2017-03-21 09:29:52Z raasch
! bugfix: rho renamed rho_ocean for the cross section output
!
! 2109 2017-01-10 12:18:08Z raasch
! bugfix: length of character string netcdf_var_name extended to avoid problems
!         which appeared in restart runs due to truncation
!
! 2040 2016-10-26 16:58:09Z gronemeier
! Increased number of possible statistic_regions to 99
! 
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
! 
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean
! 
! 2011 2016-09-19 17:29:57Z kanani
! Flag urban_surface is now defined in module control_parameters,
! changed prefix for urban surface model output to "usm_",
! introduced control parameter varnamelength for LEN of trimvar.
! 
! 2007 2016-08-24 15:47:17Z kanani
! Added support for new urban surface model (temporary modifications of 
! SELECT CASE ( ) necessary, see variable trimvar),
! increased DIMENSION of do2d_unit, do3d_unit, id_var_do2d, id_var_do3d,
! increased LEN of char_cross_profiles, var_list, var_list_old
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1990 2016-08-12 09:54:36Z gronemeier
! Bugfix: variable list was not written for time series output
! 
! 1980 2016-07-29 15:51:57Z suehring
! Bugfix, in order to steer user-defined output, setting flag found explicitly
! to .F.
! 
! 1976 2016-07-27 13:28:04Z maronga
! Removed remaining 2D land surface quantities. Definition of radiation 
! quantities is now done directly in the respective module
! 
! 1972 2016-07-26 07:52:02Z maronga
! Bugfix: wrong units for lsm quantities.
! Definition of grids for land surface quantities is now done directly in the
! respective module.
! 
! 1960 2016-07-12 16:34:24Z suehring
! Additional labels and units for timeseries output of passive scalar-related 
! quantities
! 
! 1957 2016-07-07 10:43:48Z suehring
! flight module added
!
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
! 
! 
! 1833 2016-04-07 14:23:03Z raasch
! spectrum renamed spectra_mod
!
! 1786 2016-03-08 05:49:27Z raasch
! Bugfix: id_var_time_sp made public
!
! 1783 2016-03-06 18:36:17Z raasch
! netcdf interface has been modularized, former file netcdf renamed to
! netcdf_interface, creation of netcdf-dimensions and -variables moved to
! specific new subroutines create_netcdf_dim and create_netcdf_var,
! compression (deflation) of variables implemented,
! ibmy special cpp directive removed
!
! 1745 2016-02-05 13:06:51Z gronemeier
! Bugfix: recalculating ntdim_3d, ntdim_2d_xy/xz/yz when checking the
!         extensibility of an existing file (only when using parallel NetCDF). 
! 
! 1691 2015-10-26 16:17:44Z maronga
! Added output of radiative heating rates for RRTMG. Corrected output of 
! radiative fluxes
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1596 2015-05-21 09:34:28Z gronemeier
! Bugfix in masked data output. Read 'zu_3d' when trying to extend masked data
! 
! 1551 2015-03-03 14:18:16Z maronga
! Added support for land surface model and radiation model output. In the course
! of this action a new vertical grid zs (soil) was introduced.
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1322 2014-03-20 16:38:49Z raasch
! Forgotten ONLY-attribute added to USE-statements
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1308 2014-03-13 14:58:42Z fricke
! +ntime_count, oldmode
! Adjust NF90_CREATE and NF90_OPEN statement for parallel output
! To increase the performance for parallel output, the following is done:
! - Limit time dimension
! - Values of axis data are only written by PE0
! - No fill is set for all variables
! Check the number of output time levels for restart jobs
!
! 1206 2013-07-18 12:49:16Z witha
! Bugfix: typo in preprocessor directive in subroutine open_write_netcdf_file
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1053 2012-11-13 17:11:03Z hoffmann
! +qr, nr, prr
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1031 2012-10-19 14:35:30Z raasch
! netCDF4 without parallel file support implemented, new routines
! create_netcdf_file and open_write_netcdf_file at end
!
! 992 2012-09-05 15:08:26Z hoffmann
! Removal of the informative messages PA0352 and PA0353.
!  
! 983 2012-08-21 14:17:57Z hoffmann
! Bugfix in cross_profiles.
!
! 964 2012-07-26 09:14:24Z raasch
! rev 951 and 959 reformatted
!
! 959 2012-07-24 13:13:41Z hoffmann
! Bugfix in cross_profiles. It is not allowed to arrange more than 100
! profiles with cross_profiles. 
!
! 951 2012-07-19 14:22:52Z hoffmann
! cross_profiles, profile_rows, profile_columns are written to netCDF header
!
! Revision 1.1  2005/05/18 15:37:16  raasch
! Initial revision
!
!
! Description:
! ------------
!> In case of extend = .FALSE.:
!> Define all necessary dimensions, axes and variables for the different
!> netCDF datasets. This subroutine is called from check_open after a new
!> dataset is created. It leaves the open netCDF files ready to write.
!>
!> In case of extend = .TRUE.:
!> Find out if dimensions and variables of an existing file match the values
!> of the actual run. If so, get all necessary information (ids, etc.) from
!> this file.
!>
!> Parameter av can assume values 0 (non-averaged data) and 1 (time averaged
!> data)
!>
!> @todo calculation of output time levels for parallel NetCDF still does not
!>       cover every exception (change of dt_do, end_time in restart)
!> @todo timeseries and profile output still needs to be rewritten to allow
!>       modularization
!------------------------------------------------------------------------------!
 MODULE netcdf_interface

    USE control_parameters,                                                    &
        ONLY:  max_masks, fl_max, var_fl_max, varnamelength
    USE kinds
#if defined( __netcdf )
    USE NETCDF
#endif

    PRIVATE

    INTEGER(iwp), PARAMETER ::  dopr_norm_num = 7, dopts_num = 29, dots_max = 100

    CHARACTER (LEN=6), DIMENSION(dopr_norm_num) ::  dopr_norm_names =          &
         (/ 'wpt0  ', 'ws2   ', 'tsw2  ', 'ws3   ', 'ws2tsw', 'wstsw2',        &
            'z_i   ' /)

    CHARACTER (LEN=6), DIMENSION(dopr_norm_num) ::  dopr_norm_longnames =      &
         (/ 'wpt0  ', 'w*2   ', 't*w2  ', 'w*3   ', 'w*2t*w', 'w*t*w2',        &
            'z_i   ' /)

    CHARACTER (LEN=7), DIMENSION(dopts_num) :: dopts_label =                   &
          (/ 'tnpt   ', 'x_     ', 'y_     ', 'z_     ', 'z_abs  ', 'u      ', &
             'v      ', 'w      ', 'u"     ', 'v"     ', 'w"     ', 'npt_up ', &
             'w_up   ', 'w_down ', 'radius ', 'r_min  ', 'r_max  ', 'npt_max', &
             'npt_min', 'x*2    ', 'y*2    ', 'z*2    ', 'u*2    ', 'v*2    ', &
             'w*2    ', 'u"2    ', 'v"2    ', 'w"2    ', 'npt*2  ' /)

    CHARACTER (LEN=7), DIMENSION(dopts_num) :: dopts_unit =                    &
          (/ 'number ', 'm      ', 'm      ', 'm      ', 'm      ', 'm/s    ', &
             'm/s    ', 'm/s    ', 'm/s    ', 'm/s    ', 'm/s    ', 'number ', &
             'm/s    ', 'm/s    ', 'm      ', 'm      ', 'm      ', 'number ', &
             'number ', 'm2     ', 'm2     ', 'm2     ', 'm2/s2  ', 'm2/s2  ', &
             'm2/s2  ', 'm2/s2  ', 'm2/s2  ', 'm2/s2  ', 'number2' /)

    INTEGER(iwp) ::  dots_num  = 31  !< number of timeseries defined by default
    INTEGER(iwp) ::  dots_soil = 26  !< starting index for soil-timeseries
    INTEGER(iwp) ::  dots_rad  = 32  !< starting index for radiation-timeseries

    CHARACTER (LEN=13), DIMENSION(dots_max) :: dots_label =                    &
          (/ 'E            ', 'E*           ', 'dt           ',                &
             'u*           ', 'th*          ', 'umax         ',                &
             'vmax         ', 'wmax         ', 'div_new      ',                &
             'div_old      ', 'z_i_wpt      ', 'z_i_pt       ',                &
             'w*           ', 'w"pt"0       ', 'w"pt"        ',                &
             'wpt          ', 'pt(0)        ', 'pt(z_mo)     ',                &
             'w"u"0        ', 'w"v"0        ', 'w"q"0        ',                &
             'ol           ', 'q*           ', 'w"s"         ',                &
             's*           ', 'ghf          ', 'qsws_liq     ',                &
             'qsws_soil    ', 'shf_sol      ', 'shf          ',                &
             'FW_sfc       ', 'qsws_veg     ', 'r_a          ',                &
             'r_s          ',                                                  &
             'rad_net      ', 'rad_lw_in    ', 'rad_lw_out   ',                &
             'rad_sw_in    ', 'rad_sw_out   ', 'rrtm_aldif   ',                &
             'rrtm_aldir   ', 'rrtm_asdif   ', 'rrtm_asdir   ',                &                                               
             ( 'unknown      ', i9 = 1, dots_max-43 ) /)

    CHARACTER (LEN=13), DIMENSION(dots_max) :: dots_unit =                     &
          (/ 'm2/s2        ', 'm2/s2        ', 's            ',                &
             'm/s          ', 'K            ', 'm/s          ',                &
             'm/s          ', 'm/s          ', 's-1          ',                &
             's-1          ', 'm            ', 'm            ',                &
             'm/s          ', 'K m/s        ', 'K m/s        ',                &
             'K m/s        ', 'K            ', 'K            ',                &
             'm2/s2        ', 'm2/s2        ', 'kg m/s       ',                &
             'm            ', 'kg/kg        ', 'kg m/(kg s)  ',                &
             'kg/kg        ', 'W/m2         ', 'W/m2         ',                &
             'W/m2         ', 'mK/s         ', 'mK / s       ',                &
             'm PSU/s      ', 'W/m2         ', 's/m          ',                &
             's/m          ',                                                  &
             'W/m2         ', 'W/m2         ', 'W/m2         ',                &
             'W/m2         ', 'W/m2         ', '             ',                &
             '             ', '             ', '             ',                &
             ( 'unknown      ', i9 = 1, dots_max-43 ) /)

    CHARACTER (LEN=16) :: heatflux_output_unit     !< unit for heatflux output
    CHARACTER (LEN=16) :: waterflux_output_unit    !< unit for waterflux output
    CHARACTER (LEN=16) :: momentumflux_output_unit !< unit for momentumflux output

    CHARACTER (LEN=9), DIMENSION(300) ::  dopr_unit = 'unknown'

    CHARACTER (LEN=7), DIMENSION(0:1,500) ::  do3d_unit

    CHARACTER (LEN=16), DIMENSION(25) ::  prt_var_names = &
          (/ 'pt_age          ', 'pt_dvrp_size    ', 'pt_origin_x     ', &
             'pt_origin_y     ', 'pt_origin_z     ', 'pt_radius       ', &
             'pt_speed_x      ', 'pt_speed_y      ', 'pt_speed_z      ', &
             'pt_weight_factor', 'pt_x            ', 'pt_y            ', &
             'pt_z            ', 'pt_color        ', 'pt_group        ', &
             'pt_tailpoints   ', 'pt_tail_id      ', 'pt_density_ratio', &
             'pt_exp_arg      ', 'pt_exp_term     ', 'not_used        ', &
             'not_used        ', 'not_used        ', 'not_used        ', &
             'not_used        ' /)

    CHARACTER (LEN=16), DIMENSION(25) ::  prt_var_units = &
          (/ 'seconds         ', 'meters          ', 'meters          ', &
             'meters          ', 'meters          ', 'meters          ', &
             'm/s             ', 'm/s             ', 'm/s             ', &
             'factor          ', 'meters          ', 'meters          ', &
             'meters          ', 'none            ', 'none            ', &
             'none            ', 'none            ', 'ratio           ', &
             'none            ', 'none            ', 'not_used        ', &
             'not_used        ', 'not_used        ', 'not_used        ', &
             'not_used        ' /)

    CHARACTER(LEN=20), DIMENSION(11) ::  netcdf_precision = ' '
    CHARACTER(LEN=40) ::  netcdf_data_format_string

    INTEGER(iwp) ::  id_dim_prtnum, id_dim_time_fl, id_dim_time_pr, id_dim_time_prt, &
                     id_dim_time_pts, id_dim_time_sp, id_dim_time_ts, id_dim_x_sp, &
                     id_dim_y_sp, id_dim_zu_sp, id_dim_zw_sp, id_set_fl, id_set_pr, &
                     id_set_prt, id_set_pts, id_set_sp, id_set_ts, id_var_time_fl, &
                     id_var_prtnum, id_var_rnop_prt, id_var_time_pr, id_var_time_prt, &
                     id_var_time_pts, id_var_time_sp, id_var_time_ts, id_var_x_sp, &
                     id_var_y_sp, id_var_zu_sp, id_var_zw_sp, nc_stat


    INTEGER(iwp), DIMENSION(0:1) ::  id_dim_time_xy, id_dim_time_xz, &
                    id_dim_time_yz, id_dim_time_3d, id_dim_x_xy, id_dim_xu_xy, &
                    id_dim_x_xz, id_dim_xu_xz, id_dim_x_yz, id_dim_xu_yz, &
                    id_dim_x_3d, id_dim_xu_3d, id_dim_y_xy, id_dim_yv_xy, &
                    id_dim_y_xz, id_dim_yv_xz, id_dim_y_yz, id_dim_yv_yz, &
                    id_dim_y_3d, id_dim_yv_3d, id_dim_zs_xy, id_dim_zs_xz, &
                    id_dim_zs_yz, id_dim_zs_3d, id_dim_zu_xy, id_dim_zu1_xy, &
                    id_dim_zu_xz, id_dim_zu_yz, id_dim_zu_3d, id_dim_zw_xy, &
                    id_dim_zw_xz, id_dim_zw_yz, id_dim_zw_3d, id_set_xy, &
                    id_set_xz, id_set_yz, id_set_3d, id_var_ind_x_yz, &
                    id_var_ind_y_xz, id_var_ind_z_xy, id_var_time_xy, &
                    id_var_time_xz, id_var_time_yz, id_var_time_3d, id_var_x_xy, &
                    id_var_xu_xy, id_var_x_xz, id_var_xu_xz, id_var_x_yz, &
                    id_var_xu_yz, id_var_x_3d, id_var_xu_3d, id_var_y_xy, &
                    id_var_yv_xy, id_var_y_xz, id_var_yv_xz, id_var_y_yz, &
                    id_var_yv_yz, id_var_y_3d, id_var_yv_3d, id_var_zs_xy, &
                    id_var_zs_xz, id_var_zs_yz, id_var_zs_3d, id_var_zusi_xy, &
                    id_var_zusi_3d, id_var_zu_xy, id_var_zu1_xy, id_var_zu_xz, &
                    id_var_zu_yz, id_var_zu_3d, id_var_zwwi_xy, id_var_zwwi_3d, &
                    id_var_zw_xy, id_var_zw_xz, id_var_zw_yz, id_var_zw_3d

    INTEGER ::  netcdf_data_format = 2  !< NetCDF3 64bit offset format
    INTEGER ::  netcdf_deflate = 0      !< NetCDF compression, default: no
                                        !< compression

    INTEGER(iwp)                 ::  dofl_time_count
    INTEGER(iwp), DIMENSION(10)  ::  id_var_dospx, id_var_dospy
    INTEGER(iwp), DIMENSION(20)  ::  id_var_prt
    INTEGER(iwp), DIMENSION(11)  ::  nc_precision
    INTEGER(iwp), DIMENSION(dopr_norm_num) ::  id_var_norm_dopr
    
    INTEGER(iwp), DIMENSION(fl_max) ::  id_dim_x_fl, id_dim_y_fl, id_dim_z_fl
    INTEGER(iwp), DIMENSION(fl_max) ::  id_var_x_fl, id_var_y_fl, id_var_z_fl
    
    CHARACTER (LEN=20), DIMENSION(fl_max*var_fl_max) :: dofl_label
    CHARACTER (LEN=20), DIMENSION(fl_max*var_fl_max) :: dofl_unit 
    CHARACTER (LEN=20), DIMENSION(fl_max) :: dofl_dim_label_x
    CHARACTER (LEN=20), DIMENSION(fl_max) :: dofl_dim_label_y
    CHARACTER (LEN=20), DIMENSION(fl_max) :: dofl_dim_label_z

    INTEGER(iwp), DIMENSION(fl_max*var_fl_max) :: id_var_dofl    

    INTEGER(iwp), DIMENSION(dopts_num,0:10) ::  id_var_dopts
    INTEGER(iwp), DIMENSION(0:1,500)        ::  id_var_do3d
    INTEGER(iwp), DIMENSION(100,0:99)       ::  id_dim_z_pr, id_var_dopr, &
                                                id_var_z_pr
    INTEGER(iwp), DIMENSION(dots_max,0:99)  ::  id_var_dots

!
    LOGICAL ::  output_for_t0 = .FALSE.

    INTEGER(iwp), DIMENSION(1:max_masks,0:1) ::  id_dim_time_mask, id_dim_x_mask, &
                   id_dim_xu_mask, id_dim_y_mask, id_dim_yv_mask, id_dim_zs_mask, &
                   id_dim_zu_mask, id_dim_zw_mask, &
                   id_set_mask, &
                   id_var_time_mask, id_var_x_mask, id_var_xu_mask, &
                   id_var_y_mask, id_var_yv_mask, id_var_zs_mask, &
                   id_var_zu_mask, id_var_zw_mask, &
                   id_var_zusi_mask, id_var_zwwi_mask

    REAL(wp) ::  fill_value = -999.0_wp    !< value for the _FillValue attribute


    PUBLIC   &
            dopr_unit, dopts_num,     &
            dots_label, dots_max, dots_num, dots_rad, dots_soil, dots_unit,    &
             do3d_unit, fill_value,                                  & 
            id_set_fl, id_set_mask, id_set_pr,                                 &
            id_set_prt, id_set_pts, id_set_sp, id_set_ts,&
            id_set_3d, id_var_dopr,     &
            id_var_dopts, id_var_dots, &
            id_var_do3d, id_var_norm_dopr,    &
            id_var_time_pr, id_var_time_pts, id_var_time_ts,   &
            id_var_time_3d,    &
            nc_stat,                   &
            netcdf_data_format, netcdf_data_format_string, netcdf_deflate,     &
            netcdf_precision, output_for_t0, heatflux_output_unit,             &
            waterflux_output_unit, momentumflux_output_unit

    SAVE

    INTERFACE netcdf_create_dim
       MODULE PROCEDURE netcdf_create_dim
    END INTERFACE netcdf_create_dim

    INTERFACE netcdf_create_file
       MODULE PROCEDURE netcdf_create_file
    END INTERFACE netcdf_create_file

    INTERFACE netcdf_create_var
       MODULE PROCEDURE netcdf_create_var
    END INTERFACE netcdf_create_var

    INTERFACE netcdf_define_header
       MODULE PROCEDURE netcdf_define_header
    END INTERFACE netcdf_define_header

    INTERFACE netcdf_handle_error
       MODULE PROCEDURE netcdf_handle_error
    END INTERFACE netcdf_handle_error

    INTERFACE netcdf_open_write_file
       MODULE PROCEDURE netcdf_open_write_file
    END INTERFACE netcdf_open_write_file

    PUBLIC netcdf_create_file, netcdf_define_header,                           &
           netcdf_handle_error, netcdf_open_write_file

 CONTAINS

 SUBROUTINE netcdf_define_header( callmode, extend, av )
 
#if defined( __netcdf )

    USE arrays_3d,                                                             &
        ONLY:  zu, zw

    USE constants,                                                             &
        ONLY:  pi

    USE control_parameters,                                                    &
        ONLY:  air_chemistry, averaging_interval, averaging_interval_pr,       &
               data_output_pr, dopr_n,                                 &
               dopr_time_count, dopts_time_count, dots_time_count,             &
               do3d, do3d_at_begin,   &
               dt_data_output_av,   &
               dt_do3d, mask_size,              &
               do3d_time_count, end_time, land_surface,     &
               mask_size_l, mask_i, mask_i_global, mask_j, mask_j_global,      &
               mask_k_global, message_string, mid, ntdim_2d_xy,                &
               ntdim_2d_xz, ntdim_2d_yz, ntdim_3d, nz_do3d, plant_canopy,      &
               prt_time_count, run_description_header, section, simulated_time,&
               simulated_time_at_begin, skip_time_data_output_av,              &
               skip_time_do3d, topography, num_leg, num_var_fl,                &
               urban_surface, uv_exposure

    USE grid_variables,                                                        &
        ONLY:  dx, dy, zu_s_inner, zw_w_inner

    USE indices,                                                               &
        ONLY:  nx, nxl, nxr, ny, nys, nyn, nz ,nzb, nzt

    USE kinds
    USE pegrid

    USE profil_parameter,                                                      &
        ONLY:  crmax, cross_profiles, dopr_index, profile_columns, profile_rows

    USE statistics,                                                            &
        ONLY:  hom, statistic_regions

    USE turbulence_closure_mod,                                                &
        ONLY:  tcm_define_netcdf_grid


    IMPLICIT NONE

    CHARACTER (LEN=3)              ::  suffix                !<
    CHARACTER (LEN=2), INTENT (IN) ::  callmode              !<
    CHARACTER (LEN=4)              ::  grid_x                !<
    CHARACTER (LEN=4)              ::  grid_y                !<
    CHARACTER (LEN=4)              ::  grid_z                !<
    CHARACTER (LEN=6)              ::  mode                  !<
    CHARACTER (LEN=10)             ::  precision             !<
    CHARACTER (LEN=10)             ::  var                   !<
    CHARACTER (LEN=20)             ::  netcdf_var_name       !<
    CHARACTER (LEN=varnamelength)  ::  trimvar               !< TRIM of output-variable string
    CHARACTER (LEN=80)             ::  time_average_text     !<
    CHARACTER (LEN=4000)           ::  char_cross_profiles   !<
    CHARACTER (LEN=4000)           ::  var_list              !<
    CHARACTER (LEN=4000)           ::  var_list_old          !<

    CHARACTER (LEN=100), DIMENSION(1:crmax) ::  cross_profiles_adj   !<
    CHARACTER (LEN=100), DIMENSION(1:crmax) ::  cross_profiles_char  !<

    INTEGER(iwp) ::  av                                      !<
    INTEGER(iwp) ::  cross_profiles_count                    !<
    INTEGER(iwp) ::  cross_profiles_maxi                     !<
    INTEGER(iwp) ::  delim                                   !<
    INTEGER(iwp) ::  delim_old                               !<
    INTEGER(iwp) ::  file_id                                 !<
    INTEGER(iwp) ::  i                                       !<
    INTEGER(iwp) ::  id_last                                 !<
    INTEGER(iwp) ::  id_x                                    !<
    INTEGER(iwp) ::  id_y                                    !<
    INTEGER(iwp) ::  id_z                                    !<
    INTEGER(iwp) ::  j                                       !<
    INTEGER(iwp) ::  k                                       !<
    INTEGER(iwp) ::  kk                                      !<
    INTEGER(iwp) ::  ns                                      !<
    INTEGER(iwp) ::  ns_do                                   !< actual value of ns for soil model data
    INTEGER(iwp) ::  ns_old                                  !<
    INTEGER(iwp) ::  ntime_count                             !< number of time levels found in file
    INTEGER(iwp) ::  nz_old                                  !<
    INTEGER(iwp) ::  l                                       !<

    INTEGER(iwp), SAVE ::  oldmode                           !<

    INTEGER(iwp), DIMENSION(1) ::  id_dim_time_old           !<
    INTEGER(iwp), DIMENSION(1) ::  id_dim_x_yz_old           !<
    INTEGER(iwp), DIMENSION(1) ::  id_dim_y_xz_old           !<
    INTEGER(iwp), DIMENSION(1) ::  id_dim_zu_sp_old          !<
    INTEGER(iwp), DIMENSION(1) ::  id_dim_zu_xy_old          !<
    INTEGER(iwp), DIMENSION(1) ::  id_dim_zu_3d_old          !<
    INTEGER(iwp), DIMENSION(1) ::  id_dim_zu_mask_old        !<


    INTEGER(iwp), DIMENSION(1:crmax) ::  cross_profiles_numb !<

    LOGICAL ::  found                                        !<

    LOGICAL, INTENT (INOUT) ::  extend                       !<

    LOGICAL, SAVE ::  init_netcdf = .FALSE.                  !<

    REAL(wp), DIMENSION(1) ::  last_time_coordinate          !< last time value in file

    REAL(wp), DIMENSION(:), ALLOCATABLE   ::  netcdf_data    !<
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  netcdf_data_2d !<


!
!-- Initializing actions
    IF ( .NOT. init_netcdf )  THEN
!
!--    Check and set accuracy for netCDF output. First set default value
       nc_precision = NF90_REAL4

       i = 1
       DO  WHILE ( netcdf_precision(i) /= ' ' )
          j = INDEX( netcdf_precision(i), '_' )
          IF ( j == 0 )  THEN
             WRITE ( message_string, * ) 'netcdf_precision must contain a ', &
                                         '"_"netcdf_precision(', i, ')="',   &
                                         TRIM( netcdf_precision(i) ),'"'
             CALL message( 'netcdf_define_header', 'PA0241', 2, 2, 0, 6, 0 )
          ENDIF

          var       = netcdf_precision(i)(1:j-1)
          precision = netcdf_precision(i)(j+1:)

          IF ( precision == 'NF90_REAL4' )  THEN
             j = NF90_REAL4
          ELSEIF ( precision == 'NF90_REAL8' )  THEN
             j = NF90_REAL8
          ELSE
             WRITE ( message_string, * ) 'illegal netcdf precision: ',  &
                                         'netcdf_precision(', i, ')="', &
                                         TRIM( netcdf_precision(i) ),'"'
             CALL message( 'netcdf_define_header', 'PA0242', 1, 2, 0, 6, 0 )
          ENDIF

          SELECT CASE ( var )
             CASE ( 'xy' )
                nc_precision(1) = j
             CASE ( 'xz' )
                nc_precision(2) = j
             CASE ( 'yz' )
                nc_precision(3) = j
             CASE ( '2d' )
                nc_precision(1:3) = j
             CASE ( '3d' )
                nc_precision(4) = j
             CASE ( 'pr' )
                nc_precision(5) = j
             CASE ( 'ts' )
                nc_precision(6) = j
             CASE ( 'sp' )
                nc_precision(7) = j
             CASE ( 'prt' )
                nc_precision(8) = j
             CASE ( 'masks' )
                nc_precision(11) = j
             CASE ( 'fl' )
                nc_precision(9) = j
             CASE ( 'all' )
                nc_precision    = j

             CASE DEFAULT
                WRITE ( message_string, * ) 'unknown variable in ' //          &
                                  'initialization_parameters ',                & 
                                  'assignment: netcdf_precision(', i, ')="',   &
                                            TRIM( netcdf_precision(i) ),'"'
                CALL message( 'netcdf_define_header', 'PA0243', 1, 2, 0, 6, 0 )

          END SELECT

          i = i + 1
          IF ( i > 50 )  EXIT
       ENDDO

!
!--    Check for allowed parameter range
       IF ( netcdf_deflate < 0  .OR.  netcdf_deflate > 9 )  THEN
          WRITE ( message_string, '(A,I3,A)' ) 'netcdf_deflate out of ' //     &
                                      'range & given value: ', netcdf_deflate, &
                                      ', allowed range: 0-9'
          CALL message( 'netcdf_define_header', 'PA0355', 2, 2, 0, 6, 0 )
       ENDIF
!
!--    Data compression does not work with parallel NetCDF/HDF5
       IF ( netcdf_deflate > 0  .AND.  netcdf_data_format /= 3 )  THEN
          message_string = 'netcdf_deflate reset to 0'
          CALL message( 'netcdf_define_header', 'PA0356', 0, 1, 0, 6, 0 )

          netcdf_deflate = 0
       ENDIF

       init_netcdf = .TRUE.

    ENDIF

!
!-- Determine the mode to be processed
    IF ( extend )  THEN
       mode = callmode // '_ext'
    ELSE
       mode = callmode // '_new'
    ENDIF

!
!-- Select the mode to be processed. Possibilities are 3d, ma (mask), xy, xz, 
!-- yz, pr (profiles), ps (particle timeseries), fl (flight data), ts 
!-- (timeseries) or sp (spectra)
    SELECT CASE ( mode )

       CASE ( '3d_new' )

!
!--       Define some global attributes of the dataset
          nc_stat = NF90_PUT_ATT( id_set_3d(av), NF90_GLOBAL, 'Conventions',   &
                                  'COARDS' )
          CALL netcdf_handle_error( 'netcdf_define_header', 62 )

          IF ( av == 0 )  THEN
             time_average_text = ' '
          ELSE
             WRITE (time_average_text, '('', '',F7.1,'' s average'')') &
                                                            averaging_interval
          ENDIF
          nc_stat = NF90_PUT_ATT( id_set_3d(av), NF90_GLOBAL, 'title', &
                                  TRIM( run_description_header ) //    &
                                  TRIM( time_average_text ) )
          CALL netcdf_handle_error( 'netcdf_define_header', 63 )
          IF ( av == 1 )  THEN
             WRITE ( time_average_text,'(F7.1,'' s avg'')' )  averaging_interval
             nc_stat = NF90_PUT_ATT( id_set_3d(av), NF90_GLOBAL, 'time_avg', &
                                     TRIM( time_average_text ) )
             CALL netcdf_handle_error( 'netcdf_define_header', 63 )
          ENDIF

!
!--       Define time coordinate for volume data.
!--       For parallel output the time dimensions has to be limited, otherwise
!--       the performance drops significantly.
          IF ( netcdf_data_format < 5 )  THEN
             CALL netcdf_create_dim( id_set_3d(av), 'time', NF90_UNLIMITED,    &
                                     id_dim_time_3d(av), 64 )
          ELSE
             CALL netcdf_create_dim( id_set_3d(av), 'time', ntdim_3d(av),      &
                                     id_dim_time_3d(av), 523 )
          ENDIF

          CALL netcdf_create_var( id_set_3d(av), (/ id_dim_time_3d(av) /),     &
                                  'time', NF90_DOUBLE, id_var_time_3d(av),     &
                                  'seconds', '', 65, 66, 00 )
!
!--       Define spatial dimensions and coordinates:
!--       Define vertical coordinate grid (zu grid)
          CALL netcdf_create_dim( id_set_3d(av), 'zu_3d', nz_do3d-nzb+1,       &
                                  id_dim_zu_3d(av), 67 )
          CALL netcdf_create_var( id_set_3d(av), (/ id_dim_zu_3d(av) /),       &
                                  'zu_3d', NF90_DOUBLE, id_var_zu_3d(av),      &
                                  'meters', '', 68, 69, 00 )
!
!--       Define vertical coordinate grid (zw grid)
          CALL netcdf_create_dim( id_set_3d(av), 'zw_3d', nz_do3d-nzb+1,       &
                                  id_dim_zw_3d(av), 70 )
          CALL netcdf_create_var( id_set_3d(av), (/ id_dim_zw_3d(av) /),       &
                                  'zw_3d', NF90_DOUBLE, id_var_zw_3d(av),      &
                                  'meters', '', 71, 72, 00 )
!
!--       Define x-axis (for scalar position)
          CALL netcdf_create_dim( id_set_3d(av), 'x', nx+1, id_dim_x_3d(av),   &
                                  73 )
          CALL netcdf_create_var( id_set_3d(av), (/ id_dim_x_3d(av) /), 'x',   &
                                  NF90_DOUBLE, id_var_x_3d(av), 'meters', '',  &
                                  74, 75, 00 )
!
!--       Define x-axis (for u position)
          CALL netcdf_create_dim( id_set_3d(av), 'xu', nx+1, id_dim_xu_3d(av), &
                                  358 )
          CALL netcdf_create_var( id_set_3d(av), (/ id_dim_xu_3d(av) /), 'xu', &
                                  NF90_DOUBLE, id_var_xu_3d(av), 'meters', '', &
                                  359, 360, 000 )
!
!--       Define y-axis (for scalar position)
          CALL netcdf_create_dim( id_set_3d(av), 'y', ny+1, id_dim_y_3d(av),   &
                                  76 )
          CALL netcdf_create_var( id_set_3d(av), (/ id_dim_y_3d(av) /), 'y',   &
                                  NF90_DOUBLE, id_var_y_3d(av), 'meters', '',  &
                                  77, 78, 00 )
!
!--       Define y-axis (for v position)
          CALL netcdf_create_dim( id_set_3d(av), 'yv', ny+1, id_dim_yv_3d(av), &
                                  361 )
          CALL netcdf_create_var( id_set_3d(av), (/ id_dim_yv_3d(av) /), 'yv', &
                                  NF90_DOUBLE, id_var_yv_3d(av), 'meters', '', &
                                  362, 363, 000 )
!
!--       In case of non-flat topography define 2d-arrays containing the height
!--       information. Only output 2d topography information in case of parallel
!--       output.
          IF ( TRIM( topography ) /= 'flat'  .AND.                             &
               netcdf_data_format > 4 )  THEN
!
!--          Define zusi = zu(nzb_s_inner)
             CALL netcdf_create_var( id_set_3d(av), (/ id_dim_x_3d(av),        &
                                     id_dim_y_3d(av) /), 'zusi', NF90_DOUBLE,  &
                                     id_var_zusi_3d(av), 'meters',             &
                                     'zu(nzb_s_inner)', 413, 414, 415 )
!             
!--          Define zwwi = zw(nzb_w_inner)
             CALL netcdf_create_var( id_set_3d(av), (/ id_dim_x_3d(av),        &
                                     id_dim_y_3d(av) /), 'zwwi', NF90_DOUBLE,  &
                                     id_var_zwwi_3d(av), 'meters',             &
                                     'zw(nzb_w_inner)', 416, 417, 418 )

          ENDIF             
!
!--       Define the variables
          var_list = ';'
          i = 1

          DO WHILE ( do3d(av,i)(1:1) /= ' ' )
!
!--          Temporary solution to account for data output within the new urban 
!--          surface model (urban_surface_mod.f90), see also SELECT CASE ( trimvar )
             trimvar = TRIM( do3d(av,i) )
             IF ( urban_surface  .AND.  trimvar(1:4) == 'usm_' )  THEN
                trimvar = 'usm_output'
             ENDIF
!
!--          Check for the grid
             found = .FALSE.
             SELECT CASE ( trimvar )
!
!--             Most variables are defined on the scalar grid
                CASE ( 'e', 'lpt', 'nc', 'nr', 'p', 'pc', 'pr', 'prr', 'pt',   &
                       'q', 'qc', 'ql', 'ql_c', 'ql_v', 'ql_vp', 'qr', 'qv',   &
                       'rho_ocean', 'alpha_T', 'beta_S', 'solar3d', 's', 'sa', 'vpt' )

                   grid_x = 'x'
                   grid_y = 'y'
                   grid_z = 'zu'
!
!--             u grid
                CASE ( 'u' )

                   grid_x = 'xu'
                   grid_y = 'y'
                   grid_z = 'zu'
!
!--             v grid
                CASE ( 'v' )

                   grid_x = 'x'
                   grid_y = 'yv'
                   grid_z = 'zu'
!
!--             w grid
                CASE ( 'w' )

                   grid_x = 'x'
                   grid_y = 'y'
                   grid_z = 'zw'

!             

                CASE DEFAULT

                   CALL tcm_define_netcdf_grid( do3d(av,i), found, &
                                                   grid_x, grid_y, grid_z )


                  IF ( .NOT. found )  THEN
                      WRITE ( message_string, * ) 'no grid defined for varia', &
                                                  'ble ', TRIM( do3d(av,i) )
                      CALL message( 'define_netcdf_header', 'PA0244', 0, 1, 0, &
                                    6, 0 )
                   ENDIF

             END SELECT

!
!--          Select the respective dimension ids
             IF ( grid_x == 'x' )  THEN
                id_x = id_dim_x_3d(av)
             ELSEIF ( grid_x == 'xu' )  THEN
                id_x = id_dim_xu_3d(av)
             ENDIF

             IF ( grid_y == 'y' )  THEN
                id_y = id_dim_y_3d(av)
             ELSEIF ( grid_y == 'yv' )  THEN
                id_y = id_dim_yv_3d(av)
             ENDIF

             IF ( grid_z == 'zu' )  THEN
                id_z = id_dim_zu_3d(av)
             ELSEIF ( grid_z == 'zw' )  THEN
                id_z = id_dim_zw_3d(av)
             ELSEIF ( grid_z == 'zs' )  THEN
                id_z = id_dim_zs_3d(av)
             ENDIF

!
!--          Define the grid
             CALL netcdf_create_var( id_set_3d(av),(/ id_x, id_y, id_z,        &
                                     id_dim_time_3d(av) /), do3d(av,i),        &
                                     nc_precision(4), id_var_do3d(av,i),       &
                                     TRIM( do3d_unit(av,i) ), do3d(av,i), 79,  &
                                     80, 357, .TRUE. )
#if defined( __netcdf4_parallel )
             IF ( netcdf_data_format > 4 )  THEN
!
!--             Set no fill for every variable to increase performance.
                nc_stat = NF90_DEF_VAR_FILL( id_set_3d(av),     &
                                             id_var_do3d(av,i), &
                                             1, 0 )
                CALL netcdf_handle_error( 'netcdf_define_header', 532 )
!
!--             Set collective io operations for parallel io
                nc_stat = NF90_VAR_PAR_ACCESS( id_set_3d(av),     &
                                               id_var_do3d(av,i), &
                                               NF90_COLLECTIVE )
                CALL netcdf_handle_error( 'netcdf_define_header', 445 )
             ENDIF
#endif
             var_list = TRIM( var_list ) // TRIM( do3d(av,i) ) // ';'

             i = i + 1

          ENDDO

!
!--       No arrays to output
          IF ( i == 1 )  RETURN

!
!--       Write the list of variables as global attribute (this is used by
!--       restart runs and by combine_plot_fields)
          nc_stat = NF90_PUT_ATT( id_set_3d(av), NF90_GLOBAL, 'VAR_LIST', &
                                  var_list )
          CALL netcdf_handle_error( 'netcdf_define_header', 81 )

!
!--       Set general no fill, otherwise the performance drops significantly for
!--       parallel output.
          nc_stat = NF90_SET_FILL( id_set_3d(av), NF90_NOFILL, oldmode )
          CALL netcdf_handle_error( 'netcdf_define_header', 528 )

!
!--       Leave netCDF define mode
          nc_stat = NF90_ENDDEF( id_set_3d(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 82 )

!
!--       These data are only written by PE0 for parallel output to increase
!--       the performance.
          IF ( myid == 0  .OR.  netcdf_data_format < 5 )  THEN
!
!--          Write data for x (shifted by +dx/2) and xu axis
             ALLOCATE( netcdf_data(0:nx) )

             DO  i = 0, nx
                netcdf_data(i) = ( i + 0.5 ) * dx
             ENDDO

             nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_x_3d(av),  &
                                     netcdf_data, start = (/ 1 /),    &
                                     count = (/ nx+1 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 83 )

             DO  i = 0, nx
                netcdf_data(i) = i * dx
             ENDDO

             nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_xu_3d(av), &
                                     netcdf_data, start = (/ 1 /),    &
                                     count = (/ nx+1 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 385 )

             DEALLOCATE( netcdf_data )

!
!--          Write data for y (shifted by +dy/2) and yv axis
             ALLOCATE( netcdf_data(0:ny) )

             DO  i = 0, ny
                netcdf_data(i) = ( i + 0.5_wp ) * dy
             ENDDO

             nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_y_3d(av),  &
                                     netcdf_data, start = (/ 1 /),    &
                                     count = (/ ny+1 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 84 )

             DO  i = 0, ny
                netcdf_data(i) = i * dy
             ENDDO

             nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_yv_3d(av), &
                                     netcdf_data, start = (/ 1 /),    &
                                     count = (/ ny+1 /))
             CALL netcdf_handle_error( 'netcdf_define_header', 387 )

             DEALLOCATE( netcdf_data )

!
!--          Write zu and zw data (vertical axes)
             nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zu_3d(av),  &
                                     zu(nzb:nz_do3d), start = (/ 1 /), &
                                     count = (/ nz_do3d-nzb+1 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 85 )


             nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zw_3d(av),  &
                                     zw(nzb:nz_do3d), start = (/ 1 /), &
                                     count = (/ nz_do3d-nzb+1 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 86 )
          ENDIF
!
!--       In case of non-flat topography write height information. Only for
!--       parallel netcdf output.
          IF ( TRIM( topography ) /= 'flat'  .AND.                             &
               netcdf_data_format > 4 )  THEN

!             IF ( nxr == nx  .AND.  nyn /= ny )  THEN
!                nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zusi_3d(av),     &
!                                        zu_s_inner(nxl:nxr+1,nys:nyn),         &
!                                        start = (/ nxl+1, nys+1 /),            &
!                                        count = (/ nxr-nxl+2, nyn-nys+1 /) )
!             ELSEIF ( nxr /= nx  .AND.  nyn == ny )  THEN
!                nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zusi_3d(av),     &
!                                        zu_s_inner(nxl:nxr,nys:nyn+1),         &
!                                        start = (/ nxl+1, nys+1 /),            &
!                                        count = (/ nxr-nxl+1, nyn-nys+2 /) )
!             ELSEIF ( nxr == nx  .AND.  nyn == ny )  THEN
!                nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zusi_3d(av),     &
!                                        zu_s_inner(nxl:nxr+1,nys:nyn+1),       &
!                                        start = (/ nxl+1, nys+1 /),            &
!                                        count = (/ nxr-nxl+2, nyn-nys+2 /) )
!             ELSE
                nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zusi_3d(av),     &
                                        zu_s_inner(nxl:nxr,nys:nyn),           &
                                        start = (/ nxl+1, nys+1 /),            &
                                        count = (/ nxr-nxl+1, nyn-nys+1 /) )
!             ENDIF
             CALL netcdf_handle_error( 'netcdf_define_header', 419 )

!             IF ( nxr == nx  .AND.  nyn /= ny )  THEN
!                nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zwwi_3d(av),     &
!                                        zw_w_inner(nxl:nxr+1,nys:nyn),         &
!                                        start = (/ nxl+1, nys+1 /),            &
!                                        count = (/ nxr-nxl+2, nyn-nys+1 /) )
!             ELSEIF ( nxr /= nx  .AND.  nyn == ny )  THEN
!                nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zwwi_3d(av),     &
!                                        zw_w_inner(nxl:nxr,nys:nyn+1),         &
!                                        start = (/ nxl+1, nys+1 /),            &
!                                        count = (/ nxr-nxl+1, nyn-nys+2 /) )
!             ELSEIF ( nxr == nx  .AND.  nyn == ny )  THEN
!                nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zwwi_3d(av),     &
!                                        zw_w_inner(nxl:nxr+1,nys:nyn+1),       &
!                                        start = (/ nxl+1, nys+1 /),            &
!                                        count = (/ nxr-nxl+2, nyn-nys+2 /) )
!             ELSE
                nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zwwi_3d(av),     &
                                        zw_w_inner(nxl:nxr,nys:nyn),           &
                                        start = (/ nxl+1, nys+1 /),            &
                                        count = (/ nxr-nxl+1, nyn-nys+1 /) )
!             ENDIF
             CALL netcdf_handle_error( 'netcdf_define_header', 420 )

          ENDIF

       CASE ( '3d_ext' )

!
!--       Get the list of variables and compare with the actual run.
!--       First var_list_old has to be reset, since GET_ATT does not assign
!--       trailing blanks.
          var_list_old = ' '
          nc_stat = NF90_GET_ATT( id_set_3d(av), NF90_GLOBAL, 'VAR_LIST', &
                                  var_list_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 87 )

          var_list = ';'
          i = 1
          DO WHILE ( do3d(av,i)(1:1) /= ' ' )
             var_list = TRIM(var_list) // TRIM( do3d(av,i) ) // ';'
             i = i + 1
          ENDDO

          IF ( av == 0 )  THEN
             var = '(3d)'
          ELSE
             var = '(3d_av)'
          ENDIF

          IF ( TRIM( var_list ) /= TRIM( var_list_old ) )  THEN
             message_string = 'netCDF file for volume data ' //             &
                              TRIM( var ) // ' from previous run found,' // &
                              '&but this file cannot be extended due to' // &
                              ' variable mismatch.' //                      &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0245', 0, 1, 0, 6, 0 )
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Get and compare the number of vertical gridpoints
          nc_stat = NF90_INQ_VARID( id_set_3d(av), 'zu_3d', id_var_zu_3d(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 88 )

          nc_stat = NF90_INQUIRE_VARIABLE( id_set_3d(av), id_var_zu_3d(av), &
                                           dimids = id_dim_zu_3d_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 89 )
          id_dim_zu_3d(av) = id_dim_zu_3d_old(1)

          nc_stat = NF90_INQUIRE_DIMENSION( id_set_3d(av), id_dim_zu_3d(av), &
                                            len = nz_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 90 )

          IF ( nz_do3d-nzb+1 /= nz_old )  THEN
              message_string = 'netCDF file for volume data ' //             &
                               TRIM( var ) // ' from previous run found,' // &
                               '&but this file cannot be extended due to' // &
                               ' mismatch in number of' //                   &
                               ' vertical grid points (nz_do3d).' //         &
                               '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0246', 0, 1, 0, 6, 0 )
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Get the id of the time coordinate (unlimited coordinate) and its
!--       last index on the file. The next time level is pl3d..count+1.
!--       The current time must be larger than the last output time
!--       on the file.
          nc_stat = NF90_INQ_VARID( id_set_3d(av), 'time', id_var_time_3d(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 91 )

          nc_stat = NF90_INQUIRE_VARIABLE( id_set_3d(av), id_var_time_3d(av), &
                                           dimids = id_dim_time_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 92 )

          id_dim_time_3d(av) = id_dim_time_old(1)

          nc_stat = NF90_INQUIRE_DIMENSION( id_set_3d(av), id_dim_time_3d(av), &
                                            len = ntime_count )
          CALL netcdf_handle_error( 'netcdf_define_header', 93 )

!
!--       For non-parallel output use the last output time level of the netcdf
!--       file because the time dimension is unlimited. In case of parallel
!--       output the variable ntime_count could get the value of 9*10E36 because
!--       the time dimension is limited.
          IF ( netcdf_data_format < 5 ) do3d_time_count(av) = ntime_count

          nc_stat = NF90_GET_VAR( id_set_3d(av), id_var_time_3d(av), &
                                  last_time_coordinate,              &
                                  start = (/ do3d_time_count(av) /), &
                                  count = (/ 1 /) )
          CALL netcdf_handle_error( 'netcdf_define_header', 94 )

          IF ( last_time_coordinate(1) >= simulated_time )  THEN
             message_string = 'netCDF file for volume data ' //             &
                              TRIM( var ) // ' from previous run found,' // &
                              '&but this file cannot be extended becaus' // &
                              'e the current output time' //                &
                              '&is less or equal than the last output t' // &
                              'ime on this file.' //                        &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0247', 0, 1, 0, 6, 0 )
             do3d_time_count(av) = 0
             extend = .FALSE.
             RETURN
          ENDIF

          IF ( netcdf_data_format > 4 )  THEN
!
!--          Check if the needed number of output time levels is increased
!--          compared to the number of time levels in the existing file.
             IF ( ntdim_3d(av) > ntime_count )  THEN
                message_string = 'netCDF file for volume data ' // &
                                 TRIM( var ) // ' from previous run found,' // &
                                 '&but this file cannot be extended becaus' // &
                                 'e the number of output time levels has b' // &
                                 'een increased compared to the previous s' // &
                                 'imulation.' //                               &
                                 '&New file is created instead.'
                CALL message( 'define_netcdf_header', 'PA0388', 0, 1, 0, 6, 0 )
                do3d_time_count(av) = 0
                extend = .FALSE.
!
!--             Recalculate the needed time levels for the new file.
                IF ( av == 0 )  THEN
                   ntdim_3d(0) = CEILING(                               &
                           ( end_time - MAX( skip_time_do3d,            &
                                             simulated_time_at_begin )  &
                           ) / dt_do3d )
                   IF ( do3d_at_begin )  ntdim_3d(0) = ntdim_3d(0) + 1
                ELSE
                   ntdim_3d(1) = CEILING(                               &
                           ( end_time - MAX( skip_time_data_output_av,  &
                                             simulated_time_at_begin )  &
                           ) / dt_data_output_av )
                ENDIF
                RETURN
             ENDIF
          ENDIF

!
!--       Dataset seems to be extendable.
!--       Now get the variable ids.
          i = 1
          DO WHILE ( do3d(av,i)(1:1) /= ' ' )
             nc_stat = NF90_INQ_VARID( id_set_3d(av), TRIM( do3d(av,i) ), &
                                       id_var_do3d(av,i) )
             CALL netcdf_handle_error( 'netcdf_define_header', 95 )
#if defined( __netcdf4_parallel )
!
!--          Set collective io operations for parallel io
             IF ( netcdf_data_format > 4 )  THEN
                nc_stat = NF90_VAR_PAR_ACCESS( id_set_3d(av),     &
                                               id_var_do3d(av,i), &
                                               NF90_COLLECTIVE )
                CALL netcdf_handle_error( 'netcdf_define_header', 453 )
             ENDIF
#endif
             i = i + 1
          ENDDO

!
!--       Update the title attribute on file
!--       In order to avoid 'data mode' errors if updated attributes are larger 
!--       than their original size, NF90_PUT_ATT is called in 'define mode' 
!--       enclosed by NF90_REDEF and NF90_ENDDEF calls. This implies a possible 
!--       performance loss due to data copying; an alternative strategy would be 
!--       to ensure equal attribute size. Maybe revise later.
          IF ( av == 0 )  THEN
             time_average_text = ' '
          ELSE
             WRITE (time_average_text, '('', '',F7.1,'' s average'')') &
                                                            averaging_interval
          ENDIF
          nc_stat = NF90_REDEF( id_set_3d(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 429 )
          nc_stat = NF90_PUT_ATT( id_set_3d(av), NF90_GLOBAL, 'title', &
                                  TRIM( run_description_header ) //    &
                                  TRIM( time_average_text ) )
          CALL netcdf_handle_error( 'netcdf_define_header', 96 )
          nc_stat = NF90_ENDDEF( id_set_3d(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 430 )
          message_string = 'netCDF file for volume data ' //             &
                           TRIM( var ) // ' from previous run found.' // &
                           '&This file will be extended.'
          CALL message( 'define_netcdf_header', 'PA0248', 0, 0, 0, 6, 0 )

       CASE ( 'pr_new' )

!
!--       Define some global attributes of the dataset
          IF ( averaging_interval_pr /= 0.0_wp )  THEN
             WRITE (time_average_text,'('', '',F7.1,'' s average'')')          &
                                                            averaging_interval_pr
             nc_stat = NF90_PUT_ATT( id_set_pr, NF90_GLOBAL, 'title',          &
                                     TRIM( run_description_header ) //         &
                                     TRIM( time_average_text ) )
             CALL netcdf_handle_error( 'netcdf_define_header', 218 )

             WRITE ( time_average_text,'(F7.1,'' s avg'')' ) averaging_interval_pr
             nc_stat = NF90_PUT_ATT( id_set_pr, NF90_GLOBAL, 'time_avg',       &
                                     TRIM( time_average_text ) )
          ELSE
             nc_stat = NF90_PUT_ATT( id_set_pr, NF90_GLOBAL, 'title',          &
                                     TRIM( run_description_header ) )
          ENDIF
          CALL netcdf_handle_error( 'netcdf_define_header', 219 )

!
!--       Write number of columns and rows of coordinate systems to be plotted
!--       on one page to the netcdf header.
!--       This information can be used by palmplot.
          nc_stat = NF90_PUT_ATT( id_set_pr, NF90_GLOBAL,                     &
                                  'no_rows',                                  & 
                                  profile_rows ) 
          CALL netcdf_handle_error( 'netcdf_define_header', 519 )

          nc_stat = NF90_PUT_ATT( id_set_pr, NF90_GLOBAL,                     &
                                  'no_columns',                               & 
                                  profile_columns ) 
          CALL netcdf_handle_error( 'netcdf_define_header', 520 )


          cross_profiles_adj  = ADJUSTL( cross_profiles )
          cross_profiles_numb = 999999
          cross_profiles_char = ''

!
!--       Each profile defined in cross_profiles is written to an array
!--       (cross_profiles_char). The number of the respective coordinate
!--       system is assigned in a second array (cross_profiles_numb).
          k = 1

          DO  i = 1, crmax

             IF ( TRIM( cross_profiles_adj(i) ) == ' ' )  EXIT
             delim_old = 0

             DO   j = 1, crmax
                delim = INDEX( cross_profiles_adj(i)(delim_old+1:), ' ' )
                IF ( delim == 1 )  EXIT
                kk = MIN( crmax, k )
                cross_profiles_char(kk) = cross_profiles_adj(i)(delim_old+1: &
                                                              delim_old+delim-1)
                cross_profiles_numb(kk) = i
                k = k + 1
                cross_profiles_maxi  = i
                delim_old = delim_old + delim
             ENDDO

          ENDDO

          cross_profiles_count = MIN( crmax, k-1 )
!
!--       Check if all profiles defined in cross_profiles are defined in 
!--       data_output_pr. If not, they will be skipped.
          DO  i = 1, cross_profiles_count
             DO  j = 1, dopr_n

                IF ( TRIM(cross_profiles_char(i)) == TRIM(data_output_pr(j)) ) &
                THEN
                   EXIT
                ENDIF

                IF ( j == dopr_n )  THEN
                   cross_profiles_numb(i) = 999999
                ENDIF

             ENDDO
          ENDDO

          DO i = 1, crmax
             IF ( cross_profiles_numb(i) == 999999 ) THEN
                DO j = i + 1, crmax
                   IF ( cross_profiles_numb(j) /= 999999 ) THEN
                      cross_profiles_char(i) = cross_profiles_char(j)
                      cross_profiles_numb(i) = cross_profiles_numb(j)
                      cross_profiles_numb(j) = 999999
                      EXIT
                   ENDIF
                ENDDO
             ENDIF
          ENDDO

          DO i = 1, crmax-1
             IF ( cross_profiles_numb(i + 1) == 999999 ) THEN
                cross_profiles_count = i
                EXIT
             ENDIF
          ENDDO
!
!--       Check if all profiles defined in data_output_pr are defined in 
!--       cross_profiles. If not, they will be added to cross_profiles.
          DO  i = 1, dopr_n
             DO  j = 1, cross_profiles_count

                IF ( TRIM(cross_profiles_char(j)) == TRIM(data_output_pr(i)))  &
                THEN
                   EXIT
                ENDIF

                IF (( j == cross_profiles_count ) .AND.                        &
                    ( cross_profiles_count <= crmax - 1))  THEN
                   cross_profiles_count = cross_profiles_count + 1
                   cross_profiles_maxi  = cross_profiles_maxi  + 1
                   cross_profiles_char(MIN( crmax, cross_profiles_count )) =   &
                                                      TRIM( data_output_pr(i) )
                   cross_profiles_numb(MIN( crmax, cross_profiles_count )) =   &
                                                      cross_profiles_maxi
                ENDIF

             ENDDO
          ENDDO

          IF ( cross_profiles_count >= crmax )  THEN
             message_string = 'It is not allowed to arrange more than '        &
                              // '100 profiles with & cross_profiles. Apart'   &
                              // ' from that, all profiles are saved & to '    &
                              // 'the netCDF file.'
             CALL message( 'define_netcdf_header', 'PA0354', 0, 0, 0, 6, 0 )
          ENDIF

!
!--       Writing cross_profiles to netcdf header. This information can be 
!--       used by palmplot. Each profile is separated by ",", each cross is 
!--       separated by ";".
          char_cross_profiles = ';'
          id_last = 1
          cross_profiles_count = MIN( cross_profiles_count, crmax )

          DO  i = 1, cross_profiles_count

             IF ( cross_profiles_numb(i) /= 999999 )  THEN
                IF ( TRIM( char_cross_profiles ) == ';' )  THEN
                   char_cross_profiles = TRIM( char_cross_profiles ) // &
                                         TRIM( cross_profiles_char(i) )
                ELSEIF ( id_last == cross_profiles_numb(i) )  THEN
                   char_cross_profiles = TRIM( char_cross_profiles ) // &
                                         ',' // TRIM( cross_profiles_char(i) )
                ELSE
                   char_cross_profiles = TRIM( char_cross_profiles ) // &
                                         ';' // TRIM( cross_profiles_char(i) )
                ENDIF
                id_last = cross_profiles_numb(i)
             ENDIF

          ENDDO

          char_cross_profiles = TRIM( char_cross_profiles ) // ';'

          nc_stat = NF90_PUT_ATT( id_set_pr, NF90_GLOBAL, 'cross_profiles',   &
                                  TRIM( char_cross_profiles ) )
          CALL netcdf_handle_error( 'netcdf_define_header', 521 )

!
!--       Define time coordinate for profiles (unlimited dimension)
          CALL netcdf_create_dim( id_set_pr, 'time', NF90_UNLIMITED,           &
                                  id_dim_time_pr, 220 )
          CALL netcdf_create_var( id_set_pr, (/ id_dim_time_pr /), 'time',     &
                                  NF90_DOUBLE, id_var_time_pr, 'seconds', '',  &
                                  221, 222, 000 )
!
!--       Define the variables
          var_list = ';'
          DO  i = 1, dopr_n

             IF ( statistic_regions == 0 )  THEN

!
!--             Define the z-axes (each variable gets its own z-axis)
                CALL netcdf_create_dim( id_set_pr,                             &
                                        'z' // TRIM( data_output_pr(i) ),      &
                                        nzt+2-nzb, id_dim_z_pr(i,0), 223 )
                CALL netcdf_create_var( id_set_pr, (/ id_dim_z_pr(i,0) /),     &
                                        'z' // TRIM( data_output_pr(i) ),      &
                                       NF90_DOUBLE, id_var_z_pr(i,0),          &
                                       'meters', '', 224, 225, 000 )
!
!--             Define the variable
                CALL netcdf_create_var( id_set_pr, (/ id_dim_z_pr(i,0),        &
                                        id_dim_time_pr /), data_output_pr(i),  &
                                        nc_precision(5), id_var_dopr(i,0),     &
                                        TRIM( dopr_unit(i) ),                  &
                                        TRIM( data_output_pr(i) ), 226, 227,   &
                                        228 )

                var_list = TRIM( var_list ) // TRIM( data_output_pr(i) ) //  ';'

             ELSE
!
!--             If statistic regions are defined, add suffix _SR+#SR to the
!--             names
                DO  j = 0, statistic_regions
                   WRITE ( suffix, '(''_'',I2.2)' )  j

!
!--                Define the z-axes (each variable gets it own z-axis)
                   CALL netcdf_create_dim( id_set_pr, 'z' //                   &
                                           TRIM(data_output_pr(i)) // suffix,  &
                                           nzt+2-nzb, id_dim_z_pr(i,j), 229 )
                   CALL netcdf_create_var( id_set_pr, (/ id_dim_z_pr(i,j) /),  &
                                           'z' // TRIM(data_output_pr(i)) //   &
                                           suffix, NF90_DOUBLE,                &
                                           id_var_z_pr(i,j), 'meters', '',     &
                                           230, 231, 000 )
!
!--                Define the variable
                   CALL netcdf_create_var( id_set_pr, (/ id_dim_z_pr(i,j),     &
                                           id_dim_time_pr /),                  &
                                           TRIM(data_output_pr(i)) // suffix,  &
                                           nc_precision(5), id_var_dopr(i,j),  &
                                           TRIM( dopr_unit(i) ),               &
                                           TRIM( data_output_pr(i) ) //        &
                                           ' SR ', 232, 233, 234 )

                   var_list = TRIM( var_list ) // TRIM( data_output_pr(i) ) // &
                              suffix // ';'

                ENDDO

             ENDIF

          ENDDO

!
!--       Write the list of variables as global attribute (this is used by
!--       restart runs)
          nc_stat = NF90_PUT_ATT( id_set_pr, NF90_GLOBAL, 'VAR_LIST', var_list )
          CALL netcdf_handle_error( 'netcdf_define_header', 235 )

!
!--       Define normalization variables (as time series)
          DO  i = 1, dopr_norm_num

             CALL netcdf_create_var( id_set_pr, (/ id_dim_time_pr /),          &
                                     'NORM_' // TRIM( dopr_norm_names(i) ),    &
                                     nc_precision(5), id_var_norm_dopr(i),     &
                                     '', TRIM( dopr_norm_longnames(i) ), 236,  &
                                     000, 237 )

          ENDDO

!
!--       Leave netCDF define mode
          nc_stat = NF90_ENDDEF( id_set_pr )
          CALL netcdf_handle_error( 'netcdf_define_header', 238 )

!
!--       Write z-axes data
          DO  i = 1, dopr_n
             DO  j = 0, statistic_regions

                nc_stat = NF90_PUT_VAR( id_set_pr, id_var_z_pr(i,j),      &
                                        hom(nzb:nzt+1,2,dopr_index(i),0), &
                                        start = (/ 1 /),                  &
                                        count = (/ nzt-nzb+2 /) )
                CALL netcdf_handle_error( 'netcdf_define_header', 239 )

             ENDDO
          ENDDO


       CASE ( 'pr_ext' )

!
!--       Get the list of variables and compare with the actual run.
!--       First var_list_old has to be reset, since GET_ATT does not assign
!--       trailing blanks.
          var_list_old = ' '
          nc_stat = NF90_GET_ATT( id_set_pr, NF90_GLOBAL, 'VAR_LIST', &
                                  var_list_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 240 )

          var_list = ';'
          DO  i = 1, dopr_n

             IF ( statistic_regions == 0 )  THEN
                var_list = TRIM( var_list ) // TRIM( data_output_pr(i) ) // ';'
             ELSE
                DO  j = 0, statistic_regions
                   WRITE ( suffix, '(''_'',I2.2)' )  j
                   var_list = TRIM( var_list ) // TRIM( data_output_pr(i) ) // &
                              suffix // ';'
                ENDDO
             ENDIF

          ENDDO

          IF ( TRIM( var_list ) /= TRIM( var_list_old ) )  THEN
             message_string = 'netCDF file for vertical profiles ' //          &
                              'from previous run found,' //                    &
                              '&but this file cannot be extended due to' //    &
                              ' variable mismatch.' //                         &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0254', 0, 1, 0, 6, 0 )
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Get the id of the time coordinate (unlimited coordinate) and its
!--       last index on the file. The next time level is dopr..count+1.
!--       The current time must be larger than the last output time
!--       on the file.
          nc_stat = NF90_INQ_VARID( id_set_pr, 'time', id_var_time_pr )
          CALL netcdf_handle_error( 'netcdf_define_header', 241 )

          nc_stat = NF90_INQUIRE_VARIABLE( id_set_pr, id_var_time_pr, &
                                           dimids = id_dim_time_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 242 )
          id_dim_time_pr = id_dim_time_old(1)

          nc_stat = NF90_INQUIRE_DIMENSION( id_set_pr, id_dim_time_pr, &
                                            len = dopr_time_count )
          CALL netcdf_handle_error( 'netcdf_define_header', 243 )

          nc_stat = NF90_GET_VAR( id_set_pr, id_var_time_pr,        &
                                  last_time_coordinate,             &
                                  start = (/ dopr_time_count /), &
                                  count = (/ 1 /) )
          CALL netcdf_handle_error( 'netcdf_define_header', 244 )

          IF ( last_time_coordinate(1) >= simulated_time )  THEN
             message_string = 'netCDF file for vertical profiles ' //          &
                              'from previous run found,' //                    &
                              '&but this file cannot be extended becaus' //    &
                              'e the current output time' //                   &
                              '&is less or equal than the last output t' //    &
                              'ime on this file.' //                           &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0255', 0, 1, 0, 6, 0 )
             dopr_time_count = 0
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Dataset seems to be extendable.
!--       Now get the variable ids.
          i = 1
          DO  i = 1, dopr_n
 
             IF ( statistic_regions == 0 )  THEN
                nc_stat = NF90_INQ_VARID( id_set_pr, data_output_pr(i),        &
                                          id_var_dopr(i,0) )
                CALL netcdf_handle_error( 'netcdf_define_header', 245 )
             ELSE
                DO  j = 0, statistic_regions
                   WRITE ( suffix, '(''_'',I2.2)' )  j
                   netcdf_var_name = TRIM( data_output_pr(i) ) // suffix
                   nc_stat = NF90_INQ_VARID( id_set_pr, netcdf_var_name,       &
                                             id_var_dopr(i,j) )
                   CALL netcdf_handle_error( 'netcdf_define_header', 246 )
                ENDDO
             ENDIF

          ENDDO

!
!--       Get ids of the normalization variables
          DO  i = 1, dopr_norm_num
             nc_stat = NF90_INQ_VARID( id_set_pr,                             &
                                       'NORM_' // TRIM( dopr_norm_names(i) ), &
                                       id_var_norm_dopr(i) )
             CALL netcdf_handle_error( 'netcdf_define_header', 247 )
          ENDDO

!
!--       Update the title attribute on file
!--       In order to avoid 'data mode' errors if updated attributes are larger 
!--       than their original size, NF90_PUT_ATT is called in 'define mode' 
!--       enclosed by NF90_REDEF and NF90_ENDDEF calls. This implies a possible 
!--       performance loss due to data copying; an alternative strategy would be 
!--       to ensure equal attribute size in a job chain. Maybe revise later.
          IF ( averaging_interval_pr == 0.0_wp )  THEN
             time_average_text = ' '
          ELSE
             WRITE (time_average_text, '('', '',F7.1,'' s average'')') &
                                                            averaging_interval_pr
          ENDIF
          nc_stat = NF90_REDEF( id_set_pr )
          CALL netcdf_handle_error( 'netcdf_define_header', 437 )
          nc_stat = NF90_PUT_ATT( id_set_pr, NF90_GLOBAL, 'title',             &
                                  TRIM( run_description_header ) //            &
                                  TRIM( time_average_text ) )
          CALL netcdf_handle_error( 'netcdf_define_header', 248 )

          nc_stat = NF90_ENDDEF( id_set_pr )
          CALL netcdf_handle_error( 'netcdf_define_header', 438 )
          message_string = 'netCDF file for vertical profiles ' //             &
                           'from previous run found.' //                       &
                           '&This file will be extended.'
          CALL message( 'define_netcdf_header', 'PA0256', 0, 0, 0, 6, 0 )


       CASE ( 'ts_new' )

!
!--       Define some global attributes of the dataset
          nc_stat = NF90_PUT_ATT( id_set_ts, NF90_GLOBAL, 'title',             &
                                  TRIM( run_description_header ) )
          CALL netcdf_handle_error( 'netcdf_define_header', 249 )

!
!--       Define time coordinate for time series (unlimited dimension)
          CALL netcdf_create_dim( id_set_ts, 'time', NF90_UNLIMITED,           &
                                  id_dim_time_ts, 250 )
          CALL netcdf_create_var( id_set_ts, (/ id_dim_time_ts /), 'time',     &
                                  NF90_DOUBLE, id_var_time_ts, 'seconds', '',  &
                                  251, 252, 000 )
!
!--       Define the variables
          var_list = ';'
          DO  i = 1, dots_num

             IF ( statistic_regions == 0 )  THEN

                CALL netcdf_create_var( id_set_ts, (/ id_dim_time_ts /),       &
                                        dots_label(i), nc_precision(6),        &
                                        id_var_dots(i,0),                      &
                                        TRIM( dots_unit(i) ),                  &
                                        TRIM( dots_label(i) ), 253, 254, 255 )

                var_list = TRIM( var_list ) // TRIM( dots_label(i) ) // ';'

             ELSE
!
!--             If statistic regions are defined, add suffix _SR+#SR to the
!--             names
                DO  j = 0, statistic_regions
                   WRITE ( suffix, '(''_'',I2.2)' )  j

                   CALL netcdf_create_var( id_set_ts, (/ id_dim_time_ts /),    &
                                           TRIM( dots_label(i) ) // suffix,    &
                                           nc_precision(6), id_var_dots(i,j),  &
                                           TRIM( dots_unit(i) ),               &
                                           TRIM( dots_label(i) ) // ' SR ' //  &
                                           suffix(2:2), 256, 257, 347)

                   var_list = TRIM( var_list ) // TRIM( dots_label(i) ) // &
                              suffix // ';'

                ENDDO

             ENDIF

          ENDDO

!
!--       Write the list of variables as global attribute (this is used by
!--       restart runs)
          nc_stat = NF90_PUT_ATT( id_set_ts, NF90_GLOBAL, 'VAR_LIST', var_list )
          CALL netcdf_handle_error( 'netcdf_define_header', 258 )

!
!--       Leave netCDF define mode
          nc_stat = NF90_ENDDEF( id_set_ts )
          CALL netcdf_handle_error( 'netcdf_define_header', 259 )


       CASE ( 'ts_ext' )

!
!--       Get the list of variables and compare with the actual run.
!--       First var_list_old has to be reset, since GET_ATT does not assign
!--       trailing blanks.
          var_list_old = ' '
          nc_stat = NF90_GET_ATT( id_set_ts, NF90_GLOBAL, 'VAR_LIST', &
                                  var_list_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 260 )

          var_list = ';'
          i = 1
          DO  i = 1, dots_num

             IF ( statistic_regions == 0 )  THEN
                var_list = TRIM( var_list ) // TRIM( dots_label(i) ) // ';'
             ELSE
                DO  j = 0, statistic_regions
                   WRITE ( suffix, '(''_'',I2.2)' )  j
                   var_list = TRIM( var_list ) // TRIM( dots_label(i) ) //     &
                              suffix // ';'
                ENDDO
             ENDIF

          ENDDO

          IF ( TRIM( var_list ) /= TRIM( var_list_old ) )  THEN
             message_string = 'netCDF file for time series ' //                &
                              'from previous run found,' //                    &
                              '&but this file cannot be extended due to' //    &
                              ' variable mismatch.' //                         &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0257', 0, 1, 0, 6, 0 )
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Get the id of the time coordinate (unlimited coordinate) and its
!--       last index on the file. The next time level is dots..count+1.
!--       The current time must be larger than the last output time
!--       on the file.
          nc_stat = NF90_INQ_VARID( id_set_ts, 'time', id_var_time_ts )
          CALL netcdf_handle_error( 'netcdf_define_header', 261 )

          nc_stat = NF90_INQUIRE_VARIABLE( id_set_ts, id_var_time_ts,          &
                                           dimids = id_dim_time_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 262 )
          id_dim_time_ts = id_dim_time_old(1)

          nc_stat = NF90_INQUIRE_DIMENSION( id_set_ts, id_dim_time_ts,         &
                                            len = dots_time_count )
          CALL netcdf_handle_error( 'netcdf_define_header', 263 )

          nc_stat = NF90_GET_VAR( id_set_ts, id_var_time_ts,                   &
                                  last_time_coordinate,                        &
                                  start = (/ dots_time_count /),               &
                                  count = (/ 1 /) )
          CALL netcdf_handle_error( 'netcdf_define_header', 264 )

          IF ( last_time_coordinate(1) >= simulated_time )  THEN
             message_string = 'netCDF file for time series ' //                &
                              'from previous run found,' //                    &
                              '&but this file cannot be extended becaus' //    &
                              'e the current output time' //                   &
                              '&is less or equal than the last output t' //    &
                              'ime on this file.' //                           &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0258', 0, 1, 0, 6, 0 )
             dots_time_count = 0
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Dataset seems to be extendable.
!--       Now get the variable ids
          i = 1
          DO  i = 1, dots_num
 
             IF ( statistic_regions == 0 )  THEN
                nc_stat = NF90_INQ_VARID( id_set_ts, dots_label(i), &
                                          id_var_dots(i,0) )
                CALL netcdf_handle_error( 'netcdf_define_header', 265 )
             ELSE
                DO  j = 0, statistic_regions
                   WRITE ( suffix, '(''_'',I2.2)' )  j
                   netcdf_var_name = TRIM( dots_label(i) ) // suffix
                   nc_stat = NF90_INQ_VARID( id_set_ts, netcdf_var_name, &
                                             id_var_dots(i,j) )
                   CALL netcdf_handle_error( 'netcdf_define_header', 266 )
                ENDDO
             ENDIF

          ENDDO

!
!--       Update the title attribute on file
!--       In order to avoid 'data mode' errors if updated attributes are larger 
!--       than their original size, NF90_PUT_ATT is called in 'define mode' 
!--       enclosed by NF90_REDEF and NF90_ENDDEF calls. This implies a possible 
!--       performance loss due to data copying; an alternative strategy would be 
!--       to ensure equal attribute size in a job chain. Maybe revise later.
          nc_stat = NF90_REDEF( id_set_ts )
          CALL netcdf_handle_error( 'netcdf_define_header', 439 )
          nc_stat = NF90_PUT_ATT( id_set_ts, NF90_GLOBAL, 'title',             &
                                  TRIM( run_description_header ) )
          CALL netcdf_handle_error( 'netcdf_define_header', 267 )
          nc_stat = NF90_ENDDEF( id_set_ts )
          CALL netcdf_handle_error( 'netcdf_define_header', 440 )
          message_string = 'netCDF file for time series ' //                   &
                           'from previous run found.' //                       &
                           '&This file will be extended.'
          CALL message( 'define_netcdf_header', 'PA0259', 0, 0, 0, 6, 0 )
          
       CASE DEFAULT

          message_string = 'mode "' // TRIM( mode) // '" not supported'
          CALL message( 'netcdf_define_header', 'PA0270', 0, 0, 0, 6, 0 )

    END SELECT

#endif
 END SUBROUTINE netcdf_define_header


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Creates a netCDF file and give back the id. The parallel flag has to be TRUE
!> for parallel netCDF output support.
!------------------------------------------------------------------------------!
 
 SUBROUTINE netcdf_create_file( filename , id, parallel, errno )
#if defined( __netcdf )

    USE pegrid

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN)           :: errno
    INTEGER, INTENT(OUT)          :: id
    LOGICAL, INTENT(IN)           :: parallel


!
!-- Create a new netCDF output file with requested netCDF format
    IF ( netcdf_data_format == 1 )  THEN
!
!--    Classic netCDF format
       nc_stat = NF90_CREATE( filename, NF90_NOCLOBBER, id )

    ELSEIF ( netcdf_data_format == 2 )  THEN
!
!--    64bit-offset format
       nc_stat = NF90_CREATE( filename,                                        &
                              IOR( NF90_NOCLOBBER, NF90_64BIT_OFFSET ), id )

#if defined( __netcdf4 )
    ELSEIF ( netcdf_data_format == 3  .OR.                                     &
             ( .NOT. parallel  .AND.  netcdf_data_format == 5 ) )  THEN
!
!--    netCDF4/HDF5 format
       nc_stat = NF90_CREATE( filename, IOR( NF90_NOCLOBBER, NF90_NETCDF4 ), id )

    ELSEIF ( netcdf_data_format == 4  .OR.                                     &
             ( .NOT. parallel  .AND.  netcdf_data_format == 6 ) )  THEN
!
!--    netCDF4/HDF5 format with classic model flag
       nc_stat = NF90_CREATE( filename,                                        &
                              IOR( NF90_NOCLOBBER,                             &
                              IOR( NF90_CLASSIC_MODEL, NF90_HDF5 ) ), id )

#if defined( __netcdf4_parallel )
    ELSEIF ( netcdf_data_format == 5  .AND.  parallel )  THEN
!
!--    netCDF4/HDF5 format, parallel
       nc_stat = NF90_CREATE( filename,                                        &
                              IOR( NF90_NOCLOBBER,                             &
                              IOR( NF90_NETCDF4, NF90_MPIIO ) ),               &
                              id, COMM = comm2d, INFO = MPI_INFO_NULL )

    ELSEIF ( netcdf_data_format == 6  .AND.  parallel )  THEN
!
!--    netCDF4/HDF5 format with classic model flag, parallel
       nc_stat = NF90_CREATE( filename,                                        &
                              IOR( NF90_NOCLOBBER,                             &
                              IOR( NF90_MPIIO,                                 &
                              IOR( NF90_CLASSIC_MODEL, NF90_HDF5 ) ) ),        &
                              id, COMM = comm2d, INFO = MPI_INFO_NULL )

#endif
#endif
    ENDIF

    CALL netcdf_handle_error( 'netcdf_create_file', errno )
#endif
 END SUBROUTINE netcdf_create_file

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Opens an existing netCDF file for writing and gives back the id.
!> The parallel flag has to be TRUE for parallel netCDF output support.
!------------------------------------------------------------------------------!
 SUBROUTINE netcdf_open_write_file( filename, id, parallel, errno )
#if defined( __netcdf )

    USE pegrid

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN)           :: errno
    INTEGER, INTENT(OUT)          :: id
    LOGICAL, INTENT(IN)           :: parallel


    IF ( netcdf_data_format < 5  .OR.  .NOT. parallel )  THEN
       nc_stat = NF90_OPEN( filename, NF90_WRITE, id )
#if defined( __netcdf4 )
#if defined( __netcdf4_parallel )
    ELSEIF ( netcdf_data_format > 4  .AND.  parallel )  THEN
       nc_stat = NF90_OPEN( filename, IOR( NF90_WRITE, NF90_MPIIO ), id,  &
                            COMM = comm2d, INFO = MPI_INFO_NULL )
#endif
#endif
    ENDIF

    CALL netcdf_handle_error( 'netcdf_open_write_file', errno )
#endif
 END SUBROUTINE netcdf_open_write_file


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Prints out a text message corresponding to the current status.
!------------------------------------------------------------------------------!
 
 SUBROUTINE netcdf_handle_error( routine_name, errno )
#if defined( __netcdf )


    USE control_parameters,                                                    &
        ONLY:  message_string

    IMPLICIT NONE

    CHARACTER(LEN=6) ::  message_identifier
    CHARACTER(LEN=*) ::  routine_name

    INTEGER(iwp) ::  errno

    IF ( nc_stat /= NF90_NOERR )  THEN

       WRITE( message_identifier, '(''NC'',I4.4)' )  errno
       message_string = TRIM( NF90_STRERROR( nc_stat ) )

       CALL message( routine_name, message_identifier, 2, 2, 0, 6, 1 )

    ENDIF

#endif
 END SUBROUTINE netcdf_handle_error


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Create a dimension in NetCDF file
!------------------------------------------------------------------------------!

 SUBROUTINE netcdf_create_dim(ncid, dim_name, ncdim_type, ncdim_id, error_no)

#if defined( __netcdf )

    USE kinds

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  dim_name

    INTEGER, INTENT(IN)  ::  error_no
    INTEGER, INTENT(IN)  ::  ncid
    INTEGER, INTENT(OUT) ::  ncdim_id
    INTEGER, INTENT(IN)  ::  ncdim_type

!
!-- Define time coordinate for volume data (unlimited dimension)
    nc_stat = NF90_DEF_DIM( ncid, dim_name, ncdim_type, ncdim_id )
    CALL netcdf_handle_error( 'netcdf_create_dim', error_no )

#endif

 END SUBROUTINE netcdf_create_dim


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Create a one dimensional variable in specific units in NetCDF file
!------------------------------------------------------------------------------!

 SUBROUTINE netcdf_create_var( ncid, dim_id, var_name, var_type, var_id,       &
                               unit_name, long_name, error_no1, error_no2,     &
                               error_no3, fill )

#if defined( __netcdf )
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  long_name
    CHARACTER(LEN=*), INTENT(IN) ::  unit_name
    CHARACTER(LEN=*), INTENT(IN) ::  var_name

    LOGICAL, OPTIONAL ::  fill  !< indicates setting of _FillValue attribute 

    INTEGER, INTENT(IN)  ::  error_no1
    INTEGER, INTENT(IN)  ::  error_no2
    INTEGER, INTENT(IN)  ::  error_no3
    INTEGER, INTENT(IN)  ::  ncid
    INTEGER, INTENT(OUT) ::  var_id
    INTEGER, INTENT(IN)  ::  var_type

    INTEGER, DIMENSION(:), INTENT(IN) ::  dim_id

!
!-- Define variable
    nc_stat = NF90_DEF_VAR( ncid, var_name, var_type, dim_id, var_id )
    CALL netcdf_handle_error( 'netcdf_create_var', error_no1 )

#if defined( __netcdf4 )
!
!-- Check if variable should be deflate (including shuffling)
!-- and if it is possible (only NetCDF4 with HDF5 supports compression)
    IF ( netcdf_data_format > 2  .AND.  netcdf_deflate > 0 )  THEN
       nc_stat = NF90_DEF_VAR_DEFLATE( ncid, var_id, 1, 1, netcdf_deflate )
       CALL netcdf_handle_error( 'netcdf_create_var_deflate', error_no1 )
    ENDIF
#endif
!
!-- Set unit name if set
    IF ( unit_name /= '' )  THEN
       nc_stat = NF90_PUT_ATT( ncid, var_id, 'units', unit_name )
       CALL netcdf_handle_error( 'netcdf_create_var', error_no2 )
    ENDIF

!
!-- Set long name if set
    IF ( long_name /= '' )  THEN
       nc_stat = NF90_PUT_ATT( ncid, var_id, 'long_name', long_name )
       CALL netcdf_handle_error( 'netcdf_create_var', error_no3 )
    ENDIF

!
!-- Set _FillValue for all variables, except for dimension variables.
!-- Set the fill values accordingly to the corresponding output precision. 
    IF ( PRESENT( fill ) )  THEN
       IF ( var_type == NF90_REAL4 )  THEN
          nc_stat = NF90_PUT_ATT( ncid, var_id, '_FillValue',                  &
                                  REAL( fill_value, KIND = 4 ) )
          CALL netcdf_handle_error( 'netcdf_create_var', 0 )
       ELSE
          nc_stat = NF90_PUT_ATT( ncid, var_id, '_FillValue',                  &
                                  REAL( fill_value, KIND = 8 ) )
          CALL netcdf_handle_error( 'netcdf_create_var', 0 )
       ENDIF
    ENDIF

#endif
 END SUBROUTINE netcdf_create_var

 END MODULE netcdf_interface
