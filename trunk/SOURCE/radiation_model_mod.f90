!> @file radiation_model_mod.f90
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
! Copyright 2015-2018 Czech Technical University in Prague
! Copyright 2015-2018 Institute of Computer Science of the 
!                     Czech Academy of Sciences, Prague
! Copyright 1997-2018 Leibniz Universitaet Hannover
!------------------------------------------------------------------------------!
!
! Current revisions:
! ------------------
! 
! 
! Former revisions:
! -----------------
! $Id: radiation_model_mod.f90 3066 2018-06-12 08:55:55Z Giersch $
! Error message revised 
! 
! 3065 2018-06-12 07:03:02Z Giersch
! dz was replaced by dz(1), error message concerning vertical stretching was 
! added  
! 
! 3049 2018-05-29 13:52:36Z Giersch
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error message revised
! 
! 3026 2018-05-22 10:30:53Z schwenkel
! Changed the name specific humidity to mixing ratio, since we are computing
! mixing ratios.
! 
! 3016 2018-05-09 10:53:37Z Giersch
! Revised structure of reading svf data according to PALM coding standard: 
! svf_code_field/len and fsvf removed, error messages PA0493 and PA0494 added,
! allocation status of output arrays checked. 
! 
! 3014 2018-05-09 08:42:38Z maronga
! Introduced plant canopy height similar to urban canopy height to limit
! the memory requirement to allocate lad.
! Deactivated automatic setting of minimum raytracing distance.
! 
! 3004 2018-04-27 12:33:25Z Giersch
! Further allocation checks implemented (averaged data will be assigned to fill
! values if no allocation happened so far) 
! 
! 2995 2018-04-19 12:13:16Z Giersch
! IF-statement in radiation_init removed so that the calculation of radiative 
! fluxes at model start is done in any case, bugfix in 
! radiation_presimulate_solar_pos (end_time is the sum of end_time and the 
! spinup_time specified in the p3d_file ), list of variables/fields that have
! to be written out or read in case of restarts has been extended
! 
! 2977 2018-04-17 10:27:57Z kanani
! Implement changes from branch radiation (r2948-2971) with minor modifications,
! plus some formatting.
! (moh.hefny):
! - replaced plant_canopy by npcbl to check tree existence to avoid weird 
!   allocation of related arrays (after domain decomposition some domains 
!   contains no trees although plant_canopy (global parameter) is still TRUE).
! - added a namelist parameter to force RTM settings
! - enabled the option to switch radiation reflections off
! - renamed surf_reflections to surface_reflections 
! - removed average_radiation flag from the namelist (now it is implicitly set 
!   in init_3d_model according to RTM)
! - edited read and write sky view factors and CSF routines to account for
!   the sub-domains which may not contain any of them
! 
! 2967 2018-04-13 11:22:08Z raasch
! bugfix: missing parallel cpp-directives added
! 
! 2964 2018-04-12 16:04:03Z Giersch
! Error message PA0491 has been introduced which could be previously found in 
! check_open. The variable numprocs_previous_run is only known in case of 
! initializing_actions == read_restart_data
! 
! 2963 2018-04-12 14:47:44Z suehring
! - Introduce index for vegetation/wall, pavement/green-wall and water/window 
!   surfaces, for clearer access of surface fraction, albedo, emissivity, etc. .
! - Minor bugfix in initialization of albedo for window surfaces
! 
! 2944 2018-04-03 16:20:18Z suehring
! Fixed bad commit
! 
! 2943 2018-04-03 16:17:10Z suehring
! No read of nsurfl from SVF file since it is calculated in
! radiation_interaction_init,
! allocation of arrays in radiation_read_svf only if not yet allocated,
! update of 2920 revision comment.
! 
! 2932 2018-03-26 09:39:22Z maronga
! renamed radiation_par to radiation_parameters
! 
! 2930 2018-03-23 16:30:46Z suehring
! Remove default surfaces from radiation model, does not make much sense to 
! apply radiation model without energy-balance solvers; Further, add check for 
! this. 
! 
! 2920 2018-03-22 11:22:01Z kanani
! - Bugfix: Initialize pcbl array (=-1)
! RTM version 2.0 (Jaroslav Resler, Pavel Krc, Mohamed Salim):
! - new major version of radiation interactions
! - substantially enhanced performance and scalability
! - processing of direct and diffuse solar radiation separated from reflected 
!   radiation, removed virtual surfaces
! - new type of sky discretization by azimuth and elevation angles
! - diffuse radiation processed cumulatively using sky view factor
! - used precalculated apparent solar positions for direct irradiance
! - added new 2D raytracing process for processing whole vertical column at once
!   to increase memory efficiency and decrease number of MPI RMA operations
! - enabled limiting the number of view factors between surfaces by the distance
!   and value
! - fixing issues induced by transferring radiation interactions from 
!   urban_surface_mod to radiation_mod
! - bugfixes and other minor enhancements
! 
! 2906 2018-03-19 08:56:40Z Giersch
! NAMELIST paramter read/write_svf_on_init have been removed, functions 
! check_open and close_file are used now for opening/closing files related to
! svf data, adjusted unit number and error numbers 
! 
! 2894 2018-03-15 09:17:58Z Giersch
! Calculations of the index range of the subdomain on file which overlaps with
! the current subdomain are already done in read_restart_data_mod
! radiation_read_restart_data was renamed to radiation_rrd_local and 
! radiation_last_actions was renamed to radiation_wrd_local, variable named 
! found has been introduced for checking if restart data was found, reading 
! of restart strings has been moved completely to read_restart_data_mod, 
! radiation_rrd_local is already inside the overlap loop programmed in 
! read_restart_data_mod, the marker *** end rad *** is not necessary anymore,
! strings and their respective lengths are written out and read now in case of 
! restart runs to get rid of prescribed character lengths (Giersch)
!
! 2809 2018-02-15 09:55:58Z suehring
! Bugfix for gfortran: Replace the function C_SIZEOF with STORAGE_SIZE
! 
! 2753 2018-01-16 14:16:49Z suehring
! Tile approach for spectral albedo implemented. 
! 
! 2746 2018-01-15 12:06:04Z suehring
! Move flag plant canopy to modules
! 
! 2724 2018-01-05 12:12:38Z maronga
! Set default of average_radiation to .FALSE.
! 
! 2723 2018-01-05 09:27:03Z maronga
! Bugfix in calculation of rad_lw_out (clear-sky). First grid level was used
! instead of the surface value
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2707 2017-12-18 18:34:46Z suehring
! Changes from last commit documented
! 
! 2706 2017-12-18 18:33:49Z suehring
! Bugfix, in average radiation case calculate exner function before using it.
!
! 2701 2017-12-15 15:40:50Z suehring
! Changes from last commit documented
! 
! 2698 2017-12-14 18:46:24Z suehring
! Bugfix in get_topography_top_index
!
! 2696 2017-12-14 17:12:51Z kanani
! - Change in file header (GPL part)
! - Improved reading/writing of SVF from/to file (BM)
! - Bugfixes concerning RRTMG as well as average_radiation options (M. Salim)
! - Revised initialization of surface albedo and some minor bugfixes (MS)
! - Update net radiation after running radiation interaction routine (MS)
! - Revisions from M Salim included
! - Adjustment to topography and surface structure (MS)
! - Initialization of albedo and surface emissivity via input file (MS)
! - albedo_pars extended (MS)
! 
! 2604 2017-11-06 13:29:00Z schwenkel
! bugfix for calculation of effective radius using morrison microphysics
! 
! 2601 2017-11-02 16:22:46Z scharf
! added emissivity to namelist
! 
! 2575 2017-10-24 09:57:58Z maronga
! Bugfix: calculation of shortwave and longwave albedos for RRTMG swapped
! 
! 2547 2017-10-16 12:41:56Z schwenkel
! extended by cloud_droplets option, minor bugfix and correct calculation of
! cloud droplet number concentration
! 
! 2544 2017-10-13 18:09:32Z maronga
! Moved date and time quantitis to separate module date_and_time_mod
! 
! 2512 2017-10-04 08:26:59Z raasch
! upper bounds of cross section and 3d output changed from nx+1,ny+1 to nx,ny
! no output of ghost layer data
! 
! 2504 2017-09-27 10:36:13Z maronga
! Updates pavement types and albedo parameters
! 
! 2328 2017-08-03 12:34:22Z maronga
! Emissivity can now be set individually for each pixel.
! Albedo type can be inferred from land surface model.
! Added default albedo type for bare soil
! 
! 2318 2017-07-20 17:27:44Z suehring
! Get topography top index via Function call 
! 
! 2317 2017-07-20 17:27:19Z suehring
! Improved syntax layout
! 
! 2298 2017-06-29 09:28:18Z raasch
! type of write_binary changed from CHARACTER to LOGICAL
! 
! 2296 2017-06-28 07:53:56Z maronga
! Added output of rad_sw_out for radiation_scheme = 'constant'
! 
! 2270 2017-06-09 12:18:47Z maronga
! Numbering changed (2 timeseries removed)
! 
! 2249 2017-06-06 13:58:01Z sward
! Allow for RRTMG runs without humidity/cloud physics
! 
! 2248 2017-06-06 13:52:54Z sward
! Error no changed
! 
! 2233 2017-05-30 18:08:54Z suehring
! 
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography concept
! Bugfix in read restart
! 
! 2200 2017-04-11 11:37:51Z suehring
! Bugfix in call of exchange_horiz_2d and read restart data
! 
! 2163 2017-03-01 13:23:15Z schwenkel
! Bugfix in radiation_check_data_output
! 
! 2157 2017-02-22 15:10:35Z suehring
! Bugfix in read_restart data
! 
! 2011 2016-09-19 17:29:57Z kanani
! Removed CALL of auxiliary SUBROUTINE get_usm_info,
! flag urban_surface is now defined in module control_parameters.
! 
! 2007 2016-08-24 15:47:17Z kanani
! Added calculation of solar directional vector for new urban surface 
! model,
! accounted for urban_surface model in radiation_check_parameters,
! correction of comments for zenith angle.
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1976 2016-07-27 13:28:04Z maronga
! Output of 2D/3D/masked data is now directly done within this module. The
! radiation schemes have been simplified for better usability so that
! rad_lw_in, rad_lw_out, rad_sw_in, and rad_sw_out are available independent of
! the radiation code used.
! 
! 1856 2016-04-13 12:56:17Z maronga
! Bugfix: allocation of rad_lw_out for radiation_scheme = 'clear-sky'
! 
! 1853 2016-04-11 09:00:35Z maronga
! Added routine for radiation_scheme = constant.
!  
! 1849 2016-04-08 11:33:18Z hoffmann 
! Adapted for modularization of microphysics
!
! 1826 2016-04-07 12:01:39Z maronga
! Further modularization.
! 
! 1788 2016-03-10 11:01:04Z maronga
! Added new albedo class for pavements / roads.
!
! 1783 2016-03-06 18:36:17Z raasch
! palm-netcdf-module removed in order to avoid a circular module dependency,
! netcdf-variables moved to netcdf-module, new routine netcdf_handle_error_rad
! added
!
! 1757 2016-02-22 15:49:32Z maronga
! Added parameter unscheduled_radiation_calls. Bugfix: interpolation of sounding
! profiles for pressure and temperature above the LES domain.
! 
! 1709 2015-11-04 14:47:01Z maronga
! Bugfix: set initial value for rrtm_lwuflx_dt to zero, small formatting
! corrections
! 
! 1701 2015-11-02 07:43:04Z maronga
! Bugfixes: wrong index for output of timeseries, setting of nz_snd_end
! 
! 1691 2015-10-26 16:17:44Z maronga
! Added option for spin-up runs without radiation (skip_time_do_radiation). Bugfix
! in calculation of pressure profiles. Bugfix in calculation of trace gas profiles.
! Added output of radiative heating rates.
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1606 2015-06-29 10:43:37Z maronga
! Added preprocessor directive __netcdf to allow for compiling without netCDF.
! Note, however, that RRTMG cannot be used without netCDF.
! 
! 1590 2015-05-08 13:56:27Z maronga
! Bugfix: definition of character strings requires same length for all elements
! 
! 1587 2015-05-04 14:19:01Z maronga
! Added albedo class for snow
! 
! 1585 2015-04-30 07:05:52Z maronga
! Added support for RRTMG
! 
! 1571 2015-03-12 16:12:49Z maronga
! Added missing KIND attribute. Removed upper-case variable names
! 
! 1551 2015-03-03 14:18:16Z maronga
! Added support for data output. Various variables have been renamed. Added
! interface for different radiation schemes (currently: clear-sky, constant, and
! RRTM (not yet implemented).
! 
! 1496 2014-12-02 17:25:50Z maronga
! Initial revision
! 
!
! Description:
! ------------
!> Radiation models and interfaces
!> @todo Replace dz(1) appropriatly to account for grid stretching
!> @todo move variable definitions used in radiation_init only to the subroutine
!>       as they are no longer required after initialization.
!> @todo Output of full column vertical profiles used in RRTMG
!> @todo Output of other rrtm arrays (such as volume mixing ratios)
!> @todo Adapt for use with topography
!> @todo Optimize radiation_tendency routines
!>
!> @note Many variables have a leading dummy dimension (0:0) in order to
!>       match the assume-size shape expected by the RRTMG model.
!------------------------------------------------------------------------------!
 MODULE radiation_model_mod
 
    USE arrays_3d,                                                             &
        ONLY:  dzw, hyp, nc, pt, q, ql, zu, zw

    USE calc_mean_profile_mod,                                                 &
        ONLY:  calc_mean_profile

    USE cloud_parameters,                                                      &
        ONLY:  cp, l_d_cp, l_v, r_d, rho_l

    USE constants,                                                             &
        ONLY:  pi

    USE control_parameters,                                                    &
        ONLY:  cloud_droplets, cloud_physics, coupling_char, dz, g,            &
               initializing_actions, io_blocks, io_group,                      &
               latitude, longitude, large_scale_forcing, lsf_surf,             &
               message_string, microphysics_morrison, plant_canopy, pt_surface,&
               rho_surface, surface_pressure, time_since_reference_point,      &
               urban_surface, land_surface, end_time, spinup_time, dt_spinup

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s

    USE grid_variables,                                                        &
         ONLY:  ddx, ddy, dx, dy 

    USE date_and_time_mod,                                                     &
        ONLY:  calc_date_and_time, d_hours_day, d_seconds_hour, day_of_year,   &
               d_seconds_year, day_of_year_init, time_utc_init, time_utc

    USE indices,                                                               &
        ONLY:  nnx, nny, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg,   &
               nzb, nzt

    USE, INTRINSIC :: iso_c_binding

    USE kinds

    USE microphysics_mod,                                                      &
        ONLY:  na_init, nc_const, sigma_gc

#if defined ( __netcdf )
    USE NETCDF
#endif

    USE netcdf_data_input_mod,                                                 &
        ONLY:  albedo_type_f, albedo_pars_f, building_type_f, pavement_type_f, &
               vegetation_type_f, water_type_f

    USE plant_canopy_model_mod,                                                &
        ONLY:  lad_s, pc_heating_rate, pc_transpiration_rate

    USE pegrid

#if defined ( __rrtmg )
    USE parrrsw,                                                               &
        ONLY:  naerec, nbndsw

    USE parrrtm,                                                               &
        ONLY:  nbndlw

    USE rrtmg_lw_init,                                                         &
        ONLY:  rrtmg_lw_ini

    USE rrtmg_sw_init,                                                         &
        ONLY:  rrtmg_sw_ini

    USE rrtmg_lw_rad,                                                          &
        ONLY:  rrtmg_lw

    USE rrtmg_sw_rad,                                                          &
        ONLY:  rrtmg_sw
#endif
    USE statistics,                                                            &
        ONLY:  hom

    USE surface_mod,                                                           &
        ONLY:  get_topography_top_index, get_topography_top_index_ji,          &
               ind_pav_green, ind_veg_wall, ind_wat_win,                       &
               surf_lsm_h, surf_lsm_v, surf_type, surf_usm_h, surf_usm_v

    IMPLICIT NONE

    CHARACTER(10) :: radiation_scheme = 'clear-sky' ! 'constant', 'clear-sky', or 'rrtmg'

!
!-- Predefined Land surface classes (albedo_type) after Briegleb (1992)
    CHARACTER(37), DIMENSION(0:33), PARAMETER :: albedo_type_name = (/      &
                                   'user defined                         ', & !  0 
                                   'ocean                                ', & !  1
                                   'mixed farming, tall grassland        ', & !  2
                                   'tall/medium grassland                ', & !  3 
                                   'evergreen shrubland                  ', & !  4
                                   'short grassland/meadow/shrubland     ', & !  5
                                   'evergreen needleleaf forest          ', & !  6
                                   'mixed deciduous evergreen forest     ', & !  7
                                   'deciduous forest                     ', & !  8
                                   'tropical evergreen broadleaved forest', & !  9
                                   'medium/tall grassland/woodland       ', & ! 10
                                   'desert, sandy                        ', & ! 11 
                                   'desert, rocky                        ', & ! 12 
                                   'tundra                               ', & ! 13
                                   'land ice                             ', & ! 14 
                                   'sea ice                              ', & ! 15 
                                   'snow                                 ', & ! 16
                                   'bare soil                            ', & ! 17
                                   'asphalt/concrete mix                 ', & ! 18
                                   'asphalt (asphalt concrete)           ', & ! 19
                                   'concrete (Portland concrete)         ', & ! 20
                                   'sett                                 ', & ! 21
                                   'paving stones                        ', & ! 22
                                   'cobblestone                          ', & ! 23
                                   'metal                                ', & ! 24
                                   'wood                                 ', & ! 25
                                   'gravel                               ', & ! 26
                                   'fine gravel                          ', & ! 27
                                   'pebblestone                          ', & ! 28
                                   'woodchips                            ', & ! 29
                                   'tartan (sports)                      ', & ! 30
                                   'artifical turf (sports)              ', & ! 31
                                   'clay (sports)                        ', & ! 32
                                   'building (dummy)                     '  & ! 33
                                                         /)

    INTEGER(iwp) :: albedo_type  = 9999999, & !< Albedo surface type
                    dots_rad     = 0          !< starting index for timeseries output 

    LOGICAL ::  unscheduled_radiation_calls = .TRUE., & !< flag parameter indicating whether additional calls of the radiation code are allowed
                constant_albedo = .FALSE.,            & !< flag parameter indicating whether the albedo may change depending on zenith
                force_radiation_call = .FALSE.,       & !< flag parameter for unscheduled radiation calls
                lw_radiation = .TRUE.,                & !< flag parameter indicating whether longwave radiation shall be calculated
                radiation = .FALSE.,                  & !< flag parameter indicating whether the radiation model is used
                sun_up    = .TRUE.,                   & !< flag parameter indicating whether the sun is up or down
                sw_radiation = .TRUE.,                & !< flag parameter indicating whether shortwave radiation shall be calculated
                sun_direction = .FALSE.,              & !< flag parameter indicating whether solar direction shall be calculated
                average_radiation = .FALSE.,          & !< flag to set the calculation of radiation averaging for the domain
                radiation_interactions = .FALSE.,     & !< flag to activiate RTM (TRUE only if vertical urban/land surface and trees exist)
                surface_reflections = .TRUE.,         & !< flag to switch the calculation of radiation interaction between surfaces.
                                                        !< When it switched off, only the effect of buildings and trees shadow will
                                                        !< will be considered. However fewer SVFs are expected.
                radiation_interactions_on = .TRUE.      !< namelist flag to force RTM activiation regardless to vertical urban/land surface and trees


    REAL(wp), PARAMETER :: sigma_sb       = 5.67037321E-8_wp,       & !< Stefan-Boltzmann constant
                           solar_constant = 1368.0_wp                 !< solar constant at top of atmosphere

    REAL(wp) :: albedo = 9999999.9_wp,           & !< NAMELIST alpha
                albedo_lw_dif = 9999999.9_wp,    & !< NAMELIST aldif
                albedo_lw_dir = 9999999.9_wp,    & !< NAMELIST aldir
                albedo_sw_dif = 9999999.9_wp,    & !< NAMELIST asdif
                albedo_sw_dir = 9999999.9_wp,    & !< NAMELIST asdir
                decl_1,                          & !< declination coef. 1
                decl_2,                          & !< declination coef. 2
                decl_3,                          & !< declination coef. 3
                dt_radiation = 0.0_wp,           & !< radiation model timestep
                emissivity = 9999999.9_wp,       & !< NAMELIST surface emissivity
                lon = 0.0_wp,                    & !< longitude in radians
                lat = 0.0_wp,                    & !< latitude in radians
                net_radiation = 0.0_wp,          & !< net radiation at surface
                skip_time_do_radiation = 0.0_wp, & !< Radiation model is not called before this time
                sky_trans,                       & !< sky transmissivity
                time_radiation = 0.0_wp            !< time since last call of radiation code


    REAL(wp), DIMENSION(0:0) ::  zenith,         & !< cosine of solar zenith angle
                                 sun_dir_lat,    & !< solar directional vector in latitudes
                                 sun_dir_lon       !< solar directional vector in longitudes

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rad_net_av   !< average of rad_net
!
!-- Land surface albedos for solar zenith angle of 60° after Briegleb (1992)     
!-- (shortwave, longwave, broadband):   sw,      lw,      bb,
    REAL(wp), DIMENSION(0:2,1:33), PARAMETER :: albedo_pars = RESHAPE( (/& 
                                   0.06_wp, 0.06_wp, 0.06_wp,            & !  1
                                   0.09_wp, 0.28_wp, 0.19_wp,            & !  2
                                   0.11_wp, 0.33_wp, 0.23_wp,            & !  3
                                   0.11_wp, 0.33_wp, 0.23_wp,            & !  4
                                   0.14_wp, 0.34_wp, 0.25_wp,            & !  5
                                   0.06_wp, 0.22_wp, 0.14_wp,            & !  6
                                   0.06_wp, 0.27_wp, 0.17_wp,            & !  7
                                   0.06_wp, 0.31_wp, 0.19_wp,            & !  8
                                   0.06_wp, 0.22_wp, 0.14_wp,            & !  9
                                   0.06_wp, 0.28_wp, 0.18_wp,            & ! 10
                                   0.35_wp, 0.51_wp, 0.43_wp,            & ! 11
                                   0.24_wp, 0.40_wp, 0.32_wp,            & ! 12
                                   0.10_wp, 0.27_wp, 0.19_wp,            & ! 13
                                   0.90_wp, 0.65_wp, 0.77_wp,            & ! 14
                                   0.90_wp, 0.65_wp, 0.77_wp,            & ! 15
                                   0.95_wp, 0.70_wp, 0.82_wp,            & ! 16
                                   0.08_wp, 0.08_wp, 0.08_wp,            & ! 17
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 18
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 19
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 20
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 21
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 22
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 23
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 24
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 25
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 26
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 27
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 28
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 29
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 30
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 31
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 32
                                   0.17_wp, 0.17_wp, 0.17_wp             & ! 33
                                 /), (/ 3, 33 /) )

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: &
                        rad_lw_cs_hr,                  & !< longwave clear sky radiation heating rate (K/s)
                        rad_lw_cs_hr_av,               & !< average of rad_lw_cs_hr
                        rad_lw_hr,                     & !< longwave radiation heating rate (K/s)
                        rad_lw_hr_av,                  & !< average of rad_sw_hr
                        rad_lw_in,                     & !< incoming longwave radiation (W/m2)
                        rad_lw_in_av,                  & !< average of rad_lw_in
                        rad_lw_out,                    & !< outgoing longwave radiation (W/m2)
                        rad_lw_out_av,                 & !< average of rad_lw_out
                        rad_sw_cs_hr,                  & !< shortwave clear sky radiation heating rate (K/s)
                        rad_sw_cs_hr_av,               & !< average of rad_sw_cs_hr
                        rad_sw_hr,                     & !< shortwave radiation heating rate (K/s)
                        rad_sw_hr_av,                  & !< average of rad_sw_hr
                        rad_sw_in,                     & !< incoming shortwave radiation (W/m2)
                        rad_sw_in_av,                  & !< average of rad_sw_in
                        rad_sw_out,                    & !< outgoing shortwave radiation (W/m2)
                        rad_sw_out_av                    !< average of rad_sw_out


!
!-- Variables and parameters used in RRTMG only
#if defined ( __rrtmg )
    CHARACTER(LEN=12) :: rrtm_input_file = "RAD_SND_DATA" !< name of the NetCDF input file (sounding data)


!
!-- Flag parameters for RRTMGS (should not be changed)
    INTEGER(iwp), PARAMETER :: rrtm_idrv     = 1, & !< flag for longwave upward flux calculation option (0,1)
                               rrtm_inflglw  = 2, & !< flag for lw cloud optical properties (0,1,2)
                               rrtm_iceflglw = 0, & !< flag for lw ice particle specifications (0,1,2,3)
                               rrtm_liqflglw = 1, & !< flag for lw liquid droplet specifications
                               rrtm_inflgsw  = 2, & !< flag for sw cloud optical properties (0,1,2)
                               rrtm_iceflgsw = 0, & !< flag for sw ice particle specifications (0,1,2,3)
                               rrtm_liqflgsw = 1    !< flag for sw liquid droplet specifications

!
!-- The following variables should be only changed with care, as this will 
!-- require further setting of some variables, which is currently not
!-- implemented (aerosols, ice phase).
    INTEGER(iwp) :: nzt_rad,           & !< upper vertical limit for radiation calculations
                    rrtm_icld = 0,     & !< cloud flag (0: clear sky column, 1: cloudy column)
                    rrtm_iaer = 0        !< aerosol option flag (0: no aerosol layers, for lw only: 6 (requires setting of rrtm_sw_ecaer), 10: one or more aerosol layers (not implemented)

    INTEGER(iwp) :: nc_stat !< local variable for storin the result of netCDF calls for error message handling

    LOGICAL :: snd_exists = .FALSE.      !< flag parameter to check whether a user-defined input files exists

    REAL(wp), PARAMETER :: mol_mass_air_d_wv = 1.607793_wp !< molecular weight dry air / water vapor

    REAL(wp), DIMENSION(:), ALLOCATABLE :: hyp_snd,     & !< hypostatic pressure from sounding data (hPa)
                                           q_snd,       & !< mixing ratio from sounding data (kg/kg) - dummy at the moment
                                           rrtm_tsfc,   & !< dummy array for storing surface temperature
                                           t_snd          !< actual temperature from sounding data (hPa)

    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: rrtm_ccl4vmr,   & !< CCL4 volume mixing ratio (g/mol)
                                             rrtm_cfc11vmr,  & !< CFC11 volume mixing ratio (g/mol)
                                             rrtm_cfc12vmr,  & !< CFC12 volume mixing ratio (g/mol)
                                             rrtm_cfc22vmr,  & !< CFC22 volume mixing ratio (g/mol)
                                             rrtm_ch4vmr,    & !< CH4 volume mixing ratio
                                             rrtm_cicewp,    & !< in-cloud ice water path (g/m²)
                                             rrtm_cldfr,     & !< cloud fraction (0,1)
                                             rrtm_cliqwp,    & !< in-cloud liquid water path (g/m²)
                                             rrtm_co2vmr,    & !< CO2 volume mixing ratio (g/mol)
                                             rrtm_emis,      & !< surface emissivity (0-1)  
                                             rrtm_h2ovmr,    & !< H2O volume mixing ratio
                                             rrtm_n2ovmr,    & !< N2O volume mixing ratio
                                             rrtm_o2vmr,     & !< O2 volume mixing ratio
                                             rrtm_o3vmr,     & !< O3 volume mixing ratio
                                             rrtm_play,      & !< pressure layers (hPa, zu-grid)
                                             rrtm_plev,      & !< pressure layers (hPa, zw-grid)
                                             rrtm_reice,     & !< cloud ice effective radius (microns)
                                             rrtm_reliq,     & !< cloud water drop effective radius (microns)
                                             rrtm_tlay,      & !< actual temperature (K, zu-grid)
                                             rrtm_tlev,      & !< actual temperature (K, zw-grid)
                                             rrtm_lwdflx,    & !< RRTM output of incoming longwave radiation flux (W/m2)
                                             rrtm_lwdflxc,   & !< RRTM output of outgoing clear sky longwave radiation flux (W/m2) 
                                             rrtm_lwuflx,    & !< RRTM output of outgoing longwave radiation flux (W/m2)
                                             rrtm_lwuflxc,   & !< RRTM output of incoming clear sky longwave radiation flux (W/m2)
                                             rrtm_lwuflx_dt, & !< RRTM output of incoming clear sky longwave radiation flux (W/m2)
                                             rrtm_lwuflxc_dt,& !< RRTM output of outgoing clear sky longwave radiation flux (W/m2)
                                             rrtm_lwhr,      & !< RRTM output of longwave radiation heating rate (K/d)
                                             rrtm_lwhrc,     & !< RRTM output of incoming longwave clear sky radiation heating rate (K/d)
                                             rrtm_swdflx,    & !< RRTM output of incoming shortwave radiation flux (W/m2)
                                             rrtm_swdflxc,   & !< RRTM output of outgoing clear sky shortwave radiation flux (W/m2) 
                                             rrtm_swuflx,    & !< RRTM output of outgoing shortwave radiation flux (W/m2)
                                             rrtm_swuflxc,   & !< RRTM output of incoming clear sky shortwave radiation flux (W/m2)
                                             rrtm_swhr,      & !< RRTM output of shortwave radiation heating rate (K/d)
                                             rrtm_swhrc        !< RRTM output of incoming shortwave clear sky radiation heating rate (K/d) 


    REAL(wp), DIMENSION(1) ::                rrtm_aldif,     & !< surface albedo for longwave diffuse radiation
                                             rrtm_aldir,     & !< surface albedo for longwave direct radiation
                                             rrtm_asdif,     & !< surface albedo for shortwave diffuse radiation
                                             rrtm_asdir        !< surface albedo for shortwave direct radiation

!
!-- Definition of arrays that are currently not used for calling RRTMG (due to setting of flag parameters)
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  rad_lw_cs_in,   & !< incoming clear sky longwave radiation (W/m2) (not used)
                                                rad_lw_cs_out,  & !< outgoing clear sky longwave radiation (W/m2) (not used)
                                                rad_sw_cs_in,   & !< incoming clear sky shortwave radiation (W/m2) (not used)
                                                rad_sw_cs_out,  & !< outgoing clear sky shortwave radiation (W/m2) (not used)
                                                rrtm_lw_tauaer, & !< lw aerosol optical depth
                                                rrtm_lw_taucld, & !< lw in-cloud optical depth
                                                rrtm_sw_taucld, & !< sw in-cloud optical depth
                                                rrtm_sw_ssacld, & !< sw in-cloud single scattering albedo
                                                rrtm_sw_asmcld, & !< sw in-cloud asymmetry parameter
                                                rrtm_sw_fsfcld, & !< sw in-cloud forward scattering fraction
                                                rrtm_sw_tauaer, & !< sw aerosol optical depth
                                                rrtm_sw_ssaaer, & !< sw aerosol single scattering albedo
                                                rrtm_sw_asmaer, & !< sw aerosol asymmetry parameter
                                                rrtm_sw_ecaer     !< sw aerosol optical detph at 0.55 microns (rrtm_iaer = 6 only)

#endif
!
!-- Parameters of urban and land surface models
    INTEGER(iwp)                                   ::  nzu                                !< number of layers of urban surface (will be calculated)
    INTEGER(iwp)                                   ::  nzp                                !< number of layers of plant canopy (will be calculated)
    INTEGER(iwp)                                   ::  nzub,nzut                          !< bottom and top layer of urban surface (will be calculated)
    INTEGER(iwp)                                   ::  nzpt                               !< top layer of plant canopy (will be calculated)
!-- parameters of urban and land surface models
    INTEGER(iwp), PARAMETER                        ::  nzut_free = 3                      !< number of free layers above top of of topography
    INTEGER(iwp), PARAMETER                        ::  ndsvf = 2                          !< number of dimensions of real values in SVF
    INTEGER(iwp), PARAMETER                        ::  idsvf = 2                          !< number of dimensions of integer values in SVF
    INTEGER(iwp), PARAMETER                        ::  ndcsf = 2                          !< number of dimensions of real values in CSF
    INTEGER(iwp), PARAMETER                        ::  idcsf = 2                          !< number of dimensions of integer values in CSF
    INTEGER(iwp), PARAMETER                        ::  kdcsf = 4                          !< number of dimensions of integer values in CSF calculation array
    INTEGER(iwp), PARAMETER                        ::  id = 1                             !< position of d-index in surfl and surf
    INTEGER(iwp), PARAMETER                        ::  iz = 2                             !< position of k-index in surfl and surf
    INTEGER(iwp), PARAMETER                        ::  iy = 3                             !< position of j-index in surfl and surf
    INTEGER(iwp), PARAMETER                        ::  ix = 4                             !< position of i-index in surfl and surf

    INTEGER(iwp), PARAMETER                        ::  nsurf_type = 16                    !< number of surf types incl. phys.(land+urban) & (atm.,sky,boundary) surfaces - 1

    INTEGER(iwp), PARAMETER                        ::  iup_u    = 0                       !< 0 - index of urban upward surface (ground or roof)
    INTEGER(iwp), PARAMETER                        ::  idown_u  = 1                       !< 1 - index of urban downward surface (overhanging)
    INTEGER(iwp), PARAMETER                        ::  inorth_u = 2                       !< 2 - index of urban northward facing wall
    INTEGER(iwp), PARAMETER                        ::  isouth_u = 3                       !< 3 - index of urban southward facing wall
    INTEGER(iwp), PARAMETER                        ::  ieast_u  = 4                       !< 4 - index of urban eastward facing wall
    INTEGER(iwp), PARAMETER                        ::  iwest_u  = 5                       !< 5 - index of urban westward facing wall

    INTEGER(iwp), PARAMETER                        ::  iup_l    = 6                       !< 6 - index of land upward surface (ground or roof)
    INTEGER(iwp), PARAMETER                        ::  inorth_l = 7                       !< 7 - index of land northward facing wall
    INTEGER(iwp), PARAMETER                        ::  isouth_l = 8                       !< 8 - index of land southward facing wall
    INTEGER(iwp), PARAMETER                        ::  ieast_l  = 9                       !< 9 - index of land eastward facing wall
    INTEGER(iwp), PARAMETER                        ::  iwest_l  = 10                      !< 10- index of land westward facing wall

    INTEGER(iwp), PARAMETER                        ::  iup_a    = 11                      !< 11- index of atm. cell ubward virtual surface
    INTEGER(iwp), PARAMETER                        ::  idown_a  = 12                      !< 12- index of atm. cell downward virtual surface
    INTEGER(iwp), PARAMETER                        ::  inorth_a = 13                      !< 13- index of atm. cell northward facing virtual surface
    INTEGER(iwp), PARAMETER                        ::  isouth_a = 14                      !< 14- index of atm. cell southward facing virtual surface
    INTEGER(iwp), PARAMETER                        ::  ieast_a  = 15                      !< 15- index of atm. cell eastward facing virtual surface
    INTEGER(iwp), PARAMETER                        ::  iwest_a  = 16                      !< 16- index of atm. cell westward facing virtual surface

    INTEGER(iwp), DIMENSION(0:nsurf_type), PARAMETER ::  idir = (/0, 0,0, 0,1,-1,0,0, 0,1,-1,0, 0,0, 0,1,-1/)   !< surface normal direction x indices
    INTEGER(iwp), DIMENSION(0:nsurf_type), PARAMETER ::  jdir = (/0, 0,1,-1,0, 0,0,1,-1,0, 0,0, 0,1,-1,0, 0/)   !< surface normal direction y indices
    INTEGER(iwp), DIMENSION(0:nsurf_type), PARAMETER ::  kdir = (/1,-1,0, 0,0, 0,1,0, 0,0, 0,1,-1,0, 0,0, 0/)   !< surface normal direction z indices
                                                                                          !< parameter but set in the code


!-- indices and sizes of urban and land surface models
    INTEGER(iwp)                                   ::  startland        !< start index of block of land and roof surfaces
    INTEGER(iwp)                                   ::  endland          !< end index of block of land and roof surfaces
    INTEGER(iwp)                                   ::  nlands           !< number of land and roof surfaces in local processor
    INTEGER(iwp)                                   ::  startwall        !< start index of block of wall surfaces
    INTEGER(iwp)                                   ::  endwall          !< end index of block of wall surfaces
    INTEGER(iwp)                                   ::  nwalls           !< number of wall surfaces in local processor

!-- indices and sizes of urban and land surface models
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  surfl            !< coordinates of i-th local surface in local grid - surfl[:,k] = [d, z, y, x]
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  surf             !< coordinates of i-th surface in grid - surf[:,k] = [d, z, y, x]
    INTEGER(iwp)                                   ::  nsurfl           !< number of all surfaces in local processor
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE        ::  nsurfs           !< array of number of all surfaces in individual processors
    INTEGER(iwp)                                   ::  nsurf            !< global number of surfaces in index array of surfaces (nsurf = proc nsurfs)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE        ::  surfstart        !< starts of blocks of surfaces for individual processors in array surf
                                                                        !< respective block for particular processor is surfstart[iproc]+1 : surfstart[iproc+1]

!-- block variables needed for calculation of the plant canopy model inside the urban surface model
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  pct              !< top layer of the plant canopy
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  pch              !< heights of the plant canopy
    INTEGER(iwp)                                   ::  npcbl = 0        !< number of the plant canopy gridboxes in local processor
    INTEGER(wp), DIMENSION(:,:), ALLOCATABLE       ::  pcbl             !< k,j,i coordinates of l-th local plant canopy box pcbl[:,l] = [k, j, i]
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  pcbinsw          !< array of absorbed sw radiation for local plant canopy box
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  pcbinswdir       !< array of absorbed direct sw radiation for local plant canopy box
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  pcbinswdif       !< array of absorbed diffusion sw radiation for local plant canopy box
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  pcbinlw          !< array of absorbed lw radiation for local plant canopy box

!-- configuration parameters (they can be setup in PALM config)
    LOGICAL                                        ::  split_diffusion_radiation = .TRUE. !< split direct and diffusion dw radiation
                                                                                          !< (.F. in case the radiation model already does it)    
    LOGICAL                                        ::  rma_lad_raytrace = .FALSE.         !< use MPI RMA to access LAD for raytracing (instead of global array)
    LOGICAL                                        ::  mrt_factors = .FALSE.              !< whether to generate MRT factor files during init
    INTEGER(iwp)                                   ::  nrefsteps = 0                      !< number of reflection steps to perform
    REAL(wp), PARAMETER                            ::  ext_coef = 0.6_wp                  !< extinction coefficient (a.k.a. alpha)
    INTEGER(iwp), PARAMETER                        ::  rad_version_len = 10               !< length of identification string of rad version
    CHARACTER(rad_version_len), PARAMETER          ::  rad_version = 'RAD v. 1.1'         !< identification of version of binary svf and restart files
    INTEGER(iwp)                                   ::  raytrace_discrete_elevs = 40       !< number of discretization steps for elevation (nadir to zenith)
    INTEGER(iwp)                                   ::  raytrace_discrete_azims = 80       !< number of discretization steps for azimuth (out of 360 degrees)
    REAL(wp)                                       ::  max_raytracing_dist = -999.0_wp    !< maximum distance for raytracing (in metres)
    REAL(wp)                                       ::  min_irrf_value = 1e-6_wp           !< minimum potential irradiance factor value for raytracing
    REAL(wp), DIMENSION(1:30)                      ::  svfnorm_report_thresh = 1e21_wp    !< thresholds of SVF normalization values to report
    INTEGER(iwp)                                   ::  svfnorm_report_num                 !< number of SVF normalization thresholds to report

!-- radiation related arrays to be used in radiation_interaction routine
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  rad_sw_in_dir    !< direct sw radiation
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  rad_sw_in_diff   !< diffusion sw radiation
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  rad_lw_in_diff   !< diffusion lw radiation

!-- parameters required for RRTMG lower boundary condition
    REAL(wp)                   :: albedo_urb      !< albedo value retuned to RRTMG boundary cond.
    REAL(wp)                   :: emissivity_urb  !< emissivity value retuned to RRTMG boundary cond.
    REAL(wp)                   :: t_rad_urb       !< temperature value retuned to RRTMG boundary cond.

!-- type for calculation of svf
    TYPE t_svf
        INTEGER(iwp)                               :: isurflt           !< 
        INTEGER(iwp)                               :: isurfs            !< 
        REAL(wp)                                   :: rsvf              !< 
        REAL(wp)                                   :: rtransp           !< 
    END TYPE

!-- type for calculation of csf
    TYPE t_csf
        INTEGER(iwp)                               :: ip                !< 
        INTEGER(iwp)                               :: itx               !< 
        INTEGER(iwp)                               :: ity               !< 
        INTEGER(iwp)                               :: itz               !< 
        INTEGER(iwp)                               :: isurfs            !< 
        REAL(wp)                                   :: rsvf              !< 
        REAL(wp)                                   :: rtransp           !< 
    END TYPE

!-- arrays storing the values of USM
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  svfsurf          !< svfsurf[:,isvf] = index of source and target surface for svf[isvf]
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  svf              !< array of shape view factors+direct irradiation factors for local surfaces
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfins          !< array of sw radiation falling to local surface after i-th reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinl          !< array of lw radiation for local surface after i-th reflection

    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  skyvf            !< array of sky view factor for each local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  skyvft           !< array of sky view factor including transparency for each local surface
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  dsitrans         !< dsidir[isvfl,i] = path transmittance of i-th
                                                                        !< direction of direct solar irradiance per target surface
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  dsitransc        !< dtto per plant canopy box
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  dsidir           !< dsidir[:,i] = unit vector of i-th
                                                                        !< direction of direct solar irradiance
    INTEGER(iwp)                                   ::  ndsidir          !< number of apparent solar directions used
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  dsidir_rev       !< dsidir_rev[ielev,iazim] = i for dsidir or -1 if not present

    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinsw         !< array of sw radiation falling to local surface including radiation from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinlw         !< array of lw radiation falling to local surface including radiation from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinswdir      !< array of direct sw radiation falling to local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinswdif      !< array of diffuse sw radiation from sky and model boundary falling to local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinlwdif      !< array of diffuse lw radiation from sky and model boundary falling to local surface
    
                                                                        !< Outward radiation is only valid for nonvirtual surfaces
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutsl        !< array of reflected sw radiation for local surface in i-th reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutll        !< array of reflected + emitted lw radiation for local surface in i-th reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfouts         !< array of reflected sw radiation for all surfaces in i-th reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutl         !< array of reflected + emitted lw radiation for all surfaces in i-th reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutsw        !< array of total sw radiation outgoing from nonvirtual surfaces surfaces after all reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutlw        !< array of total lw radiation outgoing from nonvirtual surfaces surfaces after all reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfhf           !< array of total radiation flux incoming to minus outgoing from local surface

!-- block variables needed for calculation of the plant canopy model inside the urban surface model
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  csfsurf          !< csfsurf[:,icsf] = index of target surface and csf grid index for csf[icsf]
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  csf              !< array of plant canopy sink fators + direct irradiation factors (transparency)
    REAL(wp), DIMENSION(:,:,:), POINTER            ::  sub_lad          !< subset of lad_s within urban surface, transformed to plain Z coordinate
    REAL(wp), DIMENSION(:), POINTER                ::  sub_lad_g        !< sub_lad globalized (used to avoid MPI RMA calls in raytracing)
    REAL(wp)                                       ::  prototype_lad    !< prototype leaf area density for computing effective optical depth
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE        ::  nzterr, plantt   !< temporary global arrays for raytracing
    INTEGER(iwp)                                   ::  plantt_max

!-- arrays and variables for calculation of svf and csf
    TYPE(t_svf), DIMENSION(:), POINTER             ::  asvf             !< pointer to growing svc array
    TYPE(t_csf), DIMENSION(:), POINTER             ::  acsf             !< pointer to growing csf array
    TYPE(t_svf), DIMENSION(:), ALLOCATABLE, TARGET ::  asvf1, asvf2     !< realizations of svf array
    TYPE(t_csf), DIMENSION(:), ALLOCATABLE, TARGET ::  acsf1, acsf2     !< realizations of csf array
    INTEGER(iwp)                                   ::  nsvfla           !< dimmension of array allocated for storage of svf in local processor
    INTEGER(iwp)                                   ::  ncsfla           !< dimmension of array allocated for storage of csf in local processor
    INTEGER(iwp)                                   ::  msvf, mcsf       !< mod for swapping the growing array
    INTEGER(iwp), PARAMETER                        ::  gasize = 10000   !< initial size of growing arrays
    REAL(wp)                                       ::  dist_max_svf = -9999.0 !< maximum distance to calculate the minimum svf to be considered. It is
                                                                        !< used to avoid very small SVFs resulting from too far surfaces with mutual visibility
    INTEGER(iwp)                                   ::  nsvfl            !< number of svf for local processor
    INTEGER(iwp)                                   ::  ncsfl            !< no. of csf in local processor
                                                                        !< needed only during calc_svf but must be here because it is
                                                                        !< shared between subroutines calc_svf and raytrace
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE    ::  gridpcbl         !< index of local pcb[k,j,i]

!-- temporary arrays for calculation of csf in raytracing
    INTEGER(iwp)                                   ::  maxboxesg        !< max number of boxes ray can cross in the domain
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  boxes            !< coordinates of gridboxes being crossed by ray
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  crlens           !< array of crossing lengths of ray for particular grid boxes
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE        ::  lad_ip           !< array of numbers of process where lad is stored 
#if defined( __parallel )
    INTEGER(kind=MPI_ADDRESS_KIND), &
                  DIMENSION(:), ALLOCATABLE        ::  lad_disp         !< array of displaycements of lad in local array of proc lad_ip
#endif
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  lad_s_ray        !< array of received lad_s for appropriate gridboxes crossed by ray
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  rt2_track
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  rt2_track_lad
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  rt2_track_dist
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  rt2_dist



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-- Energy balance variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-- parameters of the land, roof and wall surfaces
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: albedo_surf        !< albedo of the surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: emiss_surf         !< emissivity of the wall surface


    INTERFACE radiation_check_data_output
       MODULE PROCEDURE radiation_check_data_output
    END INTERFACE radiation_check_data_output

    INTERFACE radiation_check_data_output_pr
       MODULE PROCEDURE radiation_check_data_output_pr
    END INTERFACE radiation_check_data_output_pr
  
    INTERFACE radiation_check_parameters
       MODULE PROCEDURE radiation_check_parameters
    END INTERFACE radiation_check_parameters
  
    INTERFACE radiation_clearsky
       MODULE PROCEDURE radiation_clearsky
    END INTERFACE radiation_clearsky
 
    INTERFACE radiation_constant
       MODULE PROCEDURE radiation_constant
    END INTERFACE radiation_constant
 
    INTERFACE radiation_control
       MODULE PROCEDURE radiation_control
    END INTERFACE radiation_control

    INTERFACE radiation_3d_data_averaging
       MODULE PROCEDURE radiation_3d_data_averaging
    END INTERFACE radiation_3d_data_averaging

    INTERFACE radiation_data_output_2d
       MODULE PROCEDURE radiation_data_output_2d
    END INTERFACE radiation_data_output_2d

    INTERFACE radiation_data_output_3d
       MODULE PROCEDURE radiation_data_output_3d
    END INTERFACE radiation_data_output_3d

    INTERFACE radiation_data_output_mask
       MODULE PROCEDURE radiation_data_output_mask
    END INTERFACE radiation_data_output_mask

    INTERFACE radiation_define_netcdf_grid
       MODULE PROCEDURE radiation_define_netcdf_grid
    END INTERFACE radiation_define_netcdf_grid

    INTERFACE radiation_header
       MODULE PROCEDURE radiation_header
    END INTERFACE radiation_header 
  
    INTERFACE radiation_init
       MODULE PROCEDURE radiation_init
    END INTERFACE radiation_init

    INTERFACE radiation_parin
       MODULE PROCEDURE radiation_parin
    END INTERFACE radiation_parin
    
    INTERFACE radiation_rrtmg
       MODULE PROCEDURE radiation_rrtmg
    END INTERFACE radiation_rrtmg

    INTERFACE radiation_tendency
       MODULE PROCEDURE radiation_tendency
       MODULE PROCEDURE radiation_tendency_ij
    END INTERFACE radiation_tendency

    INTERFACE radiation_rrd_local
       MODULE PROCEDURE radiation_rrd_local
    END INTERFACE radiation_rrd_local

    INTERFACE radiation_wrd_local
       MODULE PROCEDURE radiation_wrd_local
    END INTERFACE radiation_wrd_local

    INTERFACE radiation_interaction
       MODULE PROCEDURE radiation_interaction
    END INTERFACE radiation_interaction

    INTERFACE radiation_interaction_init
       MODULE PROCEDURE radiation_interaction_init
    END INTERFACE radiation_interaction_init
 
    INTERFACE radiation_presimulate_solar_pos
       MODULE PROCEDURE radiation_presimulate_solar_pos
    END INTERFACE radiation_presimulate_solar_pos

    INTERFACE radiation_radflux_gridbox
       MODULE PROCEDURE radiation_radflux_gridbox
    END INTERFACE radiation_radflux_gridbox

    INTERFACE radiation_calc_svf
       MODULE PROCEDURE radiation_calc_svf
    END INTERFACE radiation_calc_svf

    INTERFACE radiation_write_svf
       MODULE PROCEDURE radiation_write_svf
    END INTERFACE radiation_write_svf

    INTERFACE radiation_read_svf
       MODULE PROCEDURE radiation_read_svf
    END INTERFACE radiation_read_svf


    SAVE

    PRIVATE

!
!-- Public functions / NEEDS SORTING
    PUBLIC radiation_check_data_output, radiation_check_data_output_pr,        &
           radiation_check_parameters, radiation_control,                      &
           radiation_header, radiation_init, radiation_parin,                  &
           radiation_3d_data_averaging, radiation_tendency,                    &
           radiation_data_output_2d, radiation_data_output_3d,                 &
           radiation_define_netcdf_grid, radiation_wrd_local,                  &
           radiation_rrd_local, radiation_data_output_mask,                    &
           radiation_radflux_gridbox, radiation_calc_svf, radiation_write_svf, &
           radiation_interaction, radiation_interaction_init,                  &
           radiation_read_svf, radiation_presimulate_solar_pos
            

    
!
!-- Public variables and constants / NEEDS SORTING
    PUBLIC albedo, albedo_type, decl_1, decl_2, decl_3, dots_rad, dt_radiation,&
           emissivity, force_radiation_call,                                   &
           lat, lon, rad_net_av, radiation, radiation_scheme, rad_lw_in,       &
           rad_lw_in_av, rad_lw_out, rad_lw_out_av,                            &
           rad_lw_cs_hr, rad_lw_cs_hr_av, rad_lw_hr, rad_lw_hr_av, rad_sw_in,  &
           rad_sw_in_av, rad_sw_out, rad_sw_out_av, rad_sw_cs_hr,              &
           rad_sw_cs_hr_av, rad_sw_hr, rad_sw_hr_av, sigma_sb, solar_constant, &
           skip_time_do_radiation, time_radiation, unscheduled_radiation_calls,&
           zenith, calc_zenith, sun_direction, sun_dir_lat, sun_dir_lon,       &
           split_diffusion_radiation,                                          &
           nrefsteps, mrt_factors, dist_max_svf, nsvfl, svf,                   &
           svfsurf, surfinsw, surfinlw, surfins, surfinl, surfinswdir,         &
           surfinswdif, surfoutsw, surfoutlw, surfinlwdif, rad_sw_in_dir,      &
           rad_sw_in_diff, rad_lw_in_diff, surfouts, surfoutl, surfoutsl,      &
           surfoutll, idir, jdir, kdir, id, iz, iy, ix, nsurfs, surfstart,     &
           surf, surfl, nsurfl, pcbinswdir, pcbinswdif, pcbinsw, pcbinlw,      &
           iup_u, inorth_u, isouth_u, ieast_u, iwest_u,           &
           iup_l, inorth_l, isouth_l, ieast_l, iwest_l,                        &
           nsurf_type, nzub, nzut, nzu, pch, nsurf,                            &
           iup_a, idown_a, inorth_a, isouth_a, ieast_a, iwest_a,               &
           idsvf, ndsvf, idcsf, ndcsf, kdcsf, pct,                             &
           radiation_interactions, startwall, startland, endland, endwall,     &
           skyvf, skyvft, radiation_interactions_on, average_radiation

#if defined ( __rrtmg )
    PUBLIC rrtm_aldif, rrtm_aldir, rrtm_asdif, rrtm_asdir
#endif

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine controls the calls of the radiation schemes
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_control
 
 
       IMPLICIT NONE


       SELECT CASE ( TRIM( radiation_scheme ) )

          CASE ( 'constant' )
             CALL radiation_constant
          
          CASE ( 'clear-sky' ) 
             CALL radiation_clearsky
       
          CASE ( 'rrtmg' )
             CALL radiation_rrtmg

          CASE DEFAULT

       END SELECT


    END SUBROUTINE radiation_control

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output for radiation model
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_check_data_output( var, unit, i, ilen, k )
 
 
       USE control_parameters,                                                 &
           ONLY: data_output, message_string

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  unit     !< 
       CHARACTER (LEN=*) ::  var      !<

       INTEGER(iwp) :: i
       INTEGER(iwp) :: ilen
       INTEGER(iwp) :: k

       SELECT CASE ( TRIM( var ) )

         CASE ( 'rad_lw_cs_hr', 'rad_lw_hr', 'rad_sw_cs_hr', 'rad_sw_hr' )
             IF (  .NOT.  radiation  .OR.  radiation_scheme /= 'rrtmg' )  THEN
                message_string = '"output of "' // TRIM( var ) // '" requi' // &
                                 'res radiation = .TRUE. and ' //              &
                                 'radiation_scheme = "rrtmg"'
                CALL message( 'check_parameters', 'PA0406', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'K/h'     

          CASE ( 'rad_net*', 'rrtm_aldif*', 'rrtm_aldir*', 'rrtm_asdif*',      &
                 'rrtm_asdir*' )
             IF ( k == 0  .OR.  data_output(i)(ilen-2:ilen) /= '_xy' )  THEN
                message_string = 'illegal value for data_output: "' //         &
                                 TRIM( var ) // '" & only 2d-horizontal ' //   &
                                 'cross sections are allowed for this value'
                CALL message( 'check_parameters', 'PA0111', 1, 2, 0, 6, 0 )
             ENDIF
             IF (  .NOT.  radiation  .OR.  radiation_scheme /= "rrtmg" )  THEN
                IF ( TRIM( var ) == 'rrtm_aldif*'  .OR.                        &
                     TRIM( var ) == 'rrtm_aldir*'  .OR.                        &
                     TRIM( var ) == 'rrtm_asdif*'  .OR.                        &
                     TRIM( var ) == 'rrtm_asdir*'      )                       &
                THEN
                   message_string = 'output of "' // TRIM( var ) // '" require'&
                                    // 's radiation = .TRUE. and radiation_sch'&
                                    // 'eme = "rrtmg"'
                   CALL message( 'check_parameters', 'PA0409', 1, 2, 0, 6, 0 )
                ENDIF
             ENDIF

             IF ( TRIM( var ) == 'rad_net*'      ) unit = 'W/m2'   
             IF ( TRIM( var ) == 'rrtm_aldif*'   ) unit = ''   
             IF ( TRIM( var ) == 'rrtm_aldir*'   ) unit = '' 
             IF ( TRIM( var ) == 'rrtm_asdif*'   ) unit = '' 
             IF ( TRIM( var ) == 'rrtm_asdir*'   ) unit = '' 

          CASE DEFAULT
             unit = 'illegal'

       END SELECT


    END SUBROUTINE radiation_check_data_output

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output of profiles for radiation model
!------------------------------------------------------------------------------!  
    SUBROUTINE radiation_check_data_output_pr( variable, var_count, unit,      &
               dopr_unit )
 
       USE arrays_3d,                                                          &
           ONLY: zu

       USE control_parameters,                                                 &
           ONLY: data_output_pr, message_string

       USE indices

       USE profil_parameter

       USE statistics

       IMPLICIT NONE
   
       CHARACTER (LEN=*) ::  unit      !< 
       CHARACTER (LEN=*) ::  variable  !< 
       CHARACTER (LEN=*) ::  dopr_unit !< local value of dopr_unit
 
       INTEGER(iwp) ::  user_pr_index !< 
       INTEGER(iwp) ::  var_count     !< 

       SELECT CASE ( TRIM( variable ) )
       
         CASE ( 'rad_net' )
             IF ( (  .NOT.  radiation )  .OR.  radiation_scheme == 'constant' )&
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme = "constant"'
                CALL message( 'check_parameters', 'PA0408', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 99
                dopr_unit  = 'W/m2'
                hom(:,2,99,:)  = SPREAD( zw, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF

          CASE ( 'rad_lw_in' )
             IF ( (  .NOT.  radiation)  .OR.  radiation_scheme == 'constant' ) &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme = "constant"'
                CALL message( 'check_parameters', 'PA0408', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 100
                dopr_unit  = 'W/m2'
                hom(:,2,100,:)  = SPREAD( zw, 2, statistic_regions+1 )
                unit = dopr_unit 
             ENDIF

          CASE ( 'rad_lw_out' )
             IF ( (  .NOT. radiation )  .OR.  radiation_scheme == 'constant' ) &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme = "constant"'
                CALL message( 'check_parameters', 'PA0408', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 101
                dopr_unit  = 'W/m2'
                hom(:,2,101,:)  = SPREAD( zw, 2, statistic_regions+1 )
                unit = dopr_unit   
             ENDIF

          CASE ( 'rad_sw_in' )
             IF ( (  .NOT. radiation )  .OR.  radiation_scheme == 'constant' ) &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme = "constant"'
                CALL message( 'check_parameters', 'PA0408', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 102
                dopr_unit  = 'W/m2'
                hom(:,2,102,:)  = SPREAD( zw, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF

          CASE ( 'rad_sw_out')
             IF ( (  .NOT.  radiation )  .OR.  radiation_scheme == 'constant' )&
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme = "constant"'
                CALL message( 'check_parameters', 'PA0408', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 103
                dopr_unit  = 'W/m2'
                hom(:,2,103,:)  = SPREAD( zw, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF

          CASE ( 'rad_lw_cs_hr' )
             IF ( (  .NOT.  radiation )  .OR.  radiation_scheme /= 'rrtmg' )   &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme /= "rrtmg"'
                CALL message( 'check_parameters', 'PA0413', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 104
                dopr_unit  = 'K/h'
                hom(:,2,104,:)  = SPREAD( zu, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF

          CASE ( 'rad_lw_hr' )
             IF ( (  .NOT.  radiation )  .OR.  radiation_scheme /= 'rrtmg' )   &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme /= "rrtmg"'
                CALL message( 'check_parameters', 'PA0413', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 105
                dopr_unit  = 'K/h'
                hom(:,2,105,:)  = SPREAD( zu, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF

          CASE ( 'rad_sw_cs_hr' )
             IF ( (  .NOT.  radiation )  .OR.  radiation_scheme /= 'rrtmg' )   &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme /= "rrtmg"'
                CALL message( 'check_parameters', 'PA0413', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 106
                dopr_unit  = 'K/h'
                hom(:,2,106,:)  = SPREAD( zu, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF

          CASE ( 'rad_sw_hr' )
             IF ( (  .NOT.  radiation )  .OR.  radiation_scheme /= 'rrtmg' )   &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme /= "rrtmg"'
                CALL message( 'check_parameters', 'PA0413', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 107
                dopr_unit  = 'K/h'
                hom(:,2,107,:)  = SPREAD( zu, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF


          CASE DEFAULT
             unit = 'illegal'

       END SELECT


    END SUBROUTINE radiation_check_data_output_pr
 
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for radiation model
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_check_parameters

       USE control_parameters,                                                 &
           ONLY: land_surface, message_string, topography, urban_surface

       USE netcdf_data_input_mod,                                              &
           ONLY:  input_pids_static                 
    
       IMPLICIT NONE
       
!
!--    In case no urban-surface or land-surface model is applied, usage of 
!--    a radiation model make no sense.         
       IF ( .NOT. land_surface  .AND.  .NOT. urban_surface )  THEN
          message_string = 'Usage of radiation module is only allowed if ' //  &
                           'land-surface and/or urban-surface model is applied.'
          CALL message( 'check_parameters', 'PA0486', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( radiation_scheme /= 'constant'   .AND.                             &
            radiation_scheme /= 'clear-sky'  .AND.                             &
            radiation_scheme /= 'rrtmg' )  THEN
          message_string = 'unknown radiation_scheme = '//                     &
                           TRIM( radiation_scheme )
          CALL message( 'check_parameters', 'PA0405', 1, 2, 0, 6, 0 )
       ELSEIF ( radiation_scheme == 'rrtmg' )  THEN
#if ! defined ( __rrtmg )
          message_string = 'radiation_scheme = "rrtmg" requires ' //           & 
                           'compilation of PALM with pre-processor ' //        &
                           'directive -D__rrtmg'
          CALL message( 'check_parameters', 'PA0407', 1, 2, 0, 6, 0 )
#endif
#if defined ( __rrtmg ) && ! defined( __netcdf )
          message_string = 'radiation_scheme = "rrtmg" requires ' //           & 
                           'the use of NetCDF (preprocessor directive ' //     &
                           '-D__netcdf'
          CALL message( 'check_parameters', 'PA0412', 1, 2, 0, 6, 0 )
#endif

       ENDIF
!
!--    Checks performed only if data is given via namelist only. 
       IF ( .NOT. input_pids_static )  THEN
          IF ( albedo_type == 0  .AND.  albedo == 9999999.9_wp  .AND.          &
               radiation_scheme == 'clear-sky')  THEN
             message_string = 'radiation_scheme = "clear-sky" in combination'//& 
                              'with albedo_type = 0 requires setting of'//     &
                              'albedo /= 9999999.9'
             CALL message( 'check_parameters', 'PA0410', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( albedo_type == 0  .AND.  radiation_scheme == 'rrtmg'  .AND.     &
             ( albedo_lw_dif == 9999999.9_wp .OR. albedo_lw_dir == 9999999.9_wp&
          .OR. albedo_sw_dif == 9999999.9_wp .OR. albedo_sw_dir == 9999999.9_wp& 
             ) ) THEN
             message_string = 'radiation_scheme = "rrtmg" in combination' //   & 
                              'with albedo_type = 0 requires setting of ' //   &
                              'albedo_lw_dif /= 9999999.9' //                  &
                              'albedo_lw_dir /= 9999999.9' //                  &
                              'albedo_sw_dif /= 9999999.9 and' //              &
                              'albedo_sw_dir /= 9999999.9'
             CALL message( 'check_parameters', 'PA0411', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

!
!--    Incialize svf normalization reporting histogram
       svfnorm_report_num = 1
       DO WHILE ( svfnorm_report_thresh(svfnorm_report_num) < 1e20_wp          &
                   .AND. svfnorm_report_num <= 30 )
          svfnorm_report_num = svfnorm_report_num + 1
       ENDDO
       svfnorm_report_num = svfnorm_report_num - 1


 
    END SUBROUTINE radiation_check_parameters 
 
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the radiation model
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_init
    
       IMPLICIT NONE

       INTEGER(iwp) ::  i         !< running index x-direction 
       INTEGER(iwp) ::  ind_type  !< running index for subgrid-surface tiles 
       INTEGER(iwp) ::  ioff      !< offset in x between surface element reference grid point in atmosphere and actual surface 
       INTEGER(iwp) ::  j         !< running index y-direction 
       INTEGER(iwp) ::  joff      !< offset in y between surface element reference grid point in atmosphere and actual surface 
       INTEGER(iwp) ::  l         !< running index for orientation of vertical surfaces 
       INTEGER(iwp) ::  m         !< running index for surface elements  

!
!--    Allocate array for storing the surface net radiation
       IF ( .NOT. ALLOCATED ( surf_lsm_h%rad_net )  .AND.                      &
                  surf_lsm_h%ns > 0  )   THEN
          ALLOCATE( surf_lsm_h%rad_net(1:surf_lsm_h%ns) )
          surf_lsm_h%rad_net = 0.0_wp 
       ENDIF
       IF ( .NOT. ALLOCATED ( surf_usm_h%rad_net )  .AND.                      &
                  surf_usm_h%ns > 0  )  THEN
          ALLOCATE( surf_usm_h%rad_net(1:surf_usm_h%ns) )
          surf_usm_h%rad_net = 0.0_wp 
       ENDIF
       DO  l = 0, 3
          IF ( .NOT. ALLOCATED ( surf_lsm_v(l)%rad_net )  .AND.                &
                     surf_lsm_v(l)%ns > 0  )  THEN
             ALLOCATE( surf_lsm_v(l)%rad_net(1:surf_lsm_v(l)%ns) )
             surf_lsm_v(l)%rad_net = 0.0_wp 
          ENDIF
          IF ( .NOT. ALLOCATED ( surf_usm_v(l)%rad_net )  .AND.                &
                     surf_usm_v(l)%ns > 0  )  THEN
             ALLOCATE( surf_usm_v(l)%rad_net(1:surf_usm_v(l)%ns) )
             surf_usm_v(l)%rad_net = 0.0_wp 
          ENDIF
       ENDDO


!
!--    Allocate array for storing the surface longwave (out) radiation change
       IF ( .NOT. ALLOCATED ( surf_lsm_h%rad_lw_out_change_0 )  .AND.          &
                  surf_lsm_h%ns > 0  )   THEN
          ALLOCATE( surf_lsm_h%rad_lw_out_change_0(1:surf_lsm_h%ns) )
          surf_lsm_h%rad_lw_out_change_0 = 0.0_wp 
       ENDIF
       IF ( .NOT. ALLOCATED ( surf_usm_h%rad_lw_out_change_0 )  .AND.          &
                  surf_usm_h%ns > 0  )  THEN
          ALLOCATE( surf_usm_h%rad_lw_out_change_0(1:surf_usm_h%ns) )
          surf_usm_h%rad_lw_out_change_0 = 0.0_wp 
       ENDIF
       DO  l = 0, 3
          IF ( .NOT. ALLOCATED ( surf_lsm_v(l)%rad_lw_out_change_0 )  .AND.    &
                     surf_lsm_v(l)%ns > 0  )  THEN
             ALLOCATE( surf_lsm_v(l)%rad_lw_out_change_0(1:surf_lsm_v(l)%ns) )
             surf_lsm_v(l)%rad_lw_out_change_0 = 0.0_wp 
          ENDIF
          IF ( .NOT. ALLOCATED ( surf_usm_v(l)%rad_lw_out_change_0 )  .AND.    &
                     surf_usm_v(l)%ns > 0  )  THEN
             ALLOCATE( surf_usm_v(l)%rad_lw_out_change_0(1:surf_usm_v(l)%ns) )
             surf_usm_v(l)%rad_lw_out_change_0 = 0.0_wp 
          ENDIF
       ENDDO

!
!--    Allocate surface arrays for incoming/outgoing short/longwave radiation
       IF ( .NOT. ALLOCATED ( surf_lsm_h%rad_sw_in )  .AND.                    &
                  surf_lsm_h%ns > 0  )   THEN
          ALLOCATE( surf_lsm_h%rad_sw_in(1:surf_lsm_h%ns)  )
          ALLOCATE( surf_lsm_h%rad_sw_out(1:surf_lsm_h%ns) )
          ALLOCATE( surf_lsm_h%rad_lw_in(1:surf_lsm_h%ns)  )
          ALLOCATE( surf_lsm_h%rad_lw_out(1:surf_lsm_h%ns) )
          surf_lsm_h%rad_sw_in  = 0.0_wp 
          surf_lsm_h%rad_sw_out = 0.0_wp 
          surf_lsm_h%rad_lw_in  = 0.0_wp 
          surf_lsm_h%rad_lw_out = 0.0_wp 
       ENDIF
       IF ( .NOT. ALLOCATED ( surf_usm_h%rad_sw_in )  .AND.                    &
                  surf_usm_h%ns > 0  )  THEN
          ALLOCATE( surf_usm_h%rad_sw_in(1:surf_usm_h%ns)  )
          ALLOCATE( surf_usm_h%rad_sw_out(1:surf_usm_h%ns) )
          ALLOCATE( surf_usm_h%rad_lw_in(1:surf_usm_h%ns)  )
          ALLOCATE( surf_usm_h%rad_lw_out(1:surf_usm_h%ns) )
          surf_usm_h%rad_sw_in  = 0.0_wp 
          surf_usm_h%rad_sw_out = 0.0_wp 
          surf_usm_h%rad_lw_in  = 0.0_wp 
          surf_usm_h%rad_lw_out = 0.0_wp 
       ENDIF
       DO  l = 0, 3
          IF ( .NOT. ALLOCATED ( surf_lsm_v(l)%rad_sw_in )  .AND.              &
                     surf_lsm_v(l)%ns > 0  )  THEN
             ALLOCATE( surf_lsm_v(l)%rad_sw_in(1:surf_lsm_v(l)%ns)  )
             ALLOCATE( surf_lsm_v(l)%rad_sw_out(1:surf_lsm_v(l)%ns) )
             ALLOCATE( surf_lsm_v(l)%rad_lw_in(1:surf_lsm_v(l)%ns)  )
             ALLOCATE( surf_lsm_v(l)%rad_lw_out(1:surf_lsm_v(l)%ns) )
             surf_lsm_v(l)%rad_sw_in  = 0.0_wp 
             surf_lsm_v(l)%rad_sw_out = 0.0_wp 
             surf_lsm_v(l)%rad_lw_in  = 0.0_wp 
             surf_lsm_v(l)%rad_lw_out = 0.0_wp 
          ENDIF
          IF ( .NOT. ALLOCATED ( surf_usm_v(l)%rad_sw_in )  .AND.              &
                     surf_usm_v(l)%ns > 0  )  THEN
             ALLOCATE( surf_usm_v(l)%rad_sw_in(1:surf_usm_v(l)%ns)  )
             ALLOCATE( surf_usm_v(l)%rad_sw_out(1:surf_usm_v(l)%ns) )
             ALLOCATE( surf_usm_v(l)%rad_lw_in(1:surf_usm_v(l)%ns)  )
             ALLOCATE( surf_usm_v(l)%rad_lw_out(1:surf_usm_v(l)%ns) )
             surf_usm_v(l)%rad_sw_in  = 0.0_wp 
             surf_usm_v(l)%rad_sw_out = 0.0_wp 
             surf_usm_v(l)%rad_lw_in  = 0.0_wp 
             surf_usm_v(l)%rad_lw_out = 0.0_wp 
          ENDIF
       ENDDO
!
!--    Fix net radiation in case of radiation_scheme = 'constant'
       IF ( radiation_scheme == 'constant' )  THEN
          IF ( ALLOCATED( surf_lsm_h%rad_net ) )                               &
             surf_lsm_h%rad_net    = net_radiation
          IF ( ALLOCATED( surf_usm_h%rad_net ) )                               &
             surf_usm_h%rad_net    = net_radiation
!
!--       Todo: weight with inclination angle
          DO  l = 0, 3
             IF ( ALLOCATED( surf_lsm_v(l)%rad_net ) )                         &
                surf_lsm_v(l)%rad_net = net_radiation
             IF ( ALLOCATED( surf_usm_v(l)%rad_net ) )                         &
                surf_usm_v(l)%rad_net = net_radiation
          ENDDO
!          radiation = .FALSE.
!
!--    Calculate orbital constants
       ELSE
          decl_1 = SIN(23.45_wp * pi / 180.0_wp)
          decl_2 = 2.0_wp * pi / 365.0_wp
          decl_3 = decl_2 * 81.0_wp
          lat    = latitude * pi / 180.0_wp
          lon    = longitude * pi / 180.0_wp
       ENDIF

       IF ( radiation_scheme == 'clear-sky'  .OR.                              &
            radiation_scheme == 'constant')  THEN


!
!--       Allocate arrays for incoming/outgoing short/longwave radiation
          IF ( .NOT. ALLOCATED ( rad_sw_in ) )  THEN
             ALLOCATE ( rad_sw_in(0:0,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( .NOT. ALLOCATED ( rad_sw_out ) )  THEN
             ALLOCATE ( rad_sw_out(0:0,nysg:nyng,nxlg:nxrg) )
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_in ) )  THEN
             ALLOCATE ( rad_lw_in(0:0,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( .NOT. ALLOCATED ( rad_lw_out ) )  THEN
             ALLOCATE ( rad_lw_out(0:0,nysg:nyng,nxlg:nxrg) )
          ENDIF

!
!--       Allocate average arrays for incoming/outgoing short/longwave radiation
          IF ( .NOT. ALLOCATED ( rad_sw_in_av ) )  THEN
             ALLOCATE ( rad_sw_in_av(0:0,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( .NOT. ALLOCATED ( rad_sw_out_av ) )  THEN
             ALLOCATE ( rad_sw_out_av(0:0,nysg:nyng,nxlg:nxrg) )
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_in_av ) )  THEN
             ALLOCATE ( rad_lw_in_av(0:0,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( .NOT. ALLOCATED ( rad_lw_out_av ) )  THEN
             ALLOCATE ( rad_lw_out_av(0:0,nysg:nyng,nxlg:nxrg) )
          ENDIF
!
!--       Allocate arrays for broadband albedo, and level 1 initialization
!--       via namelist paramter, unless already allocated.
          IF ( .NOT. ALLOCATED(surf_lsm_h%albedo) )  THEN
             ALLOCATE( surf_lsm_h%albedo(0:2,1:surf_lsm_h%ns)     )
             surf_lsm_h%albedo    = albedo
          ENDIF
          IF ( .NOT. ALLOCATED(surf_usm_h%albedo) )  THEN
             ALLOCATE( surf_usm_h%albedo(0:2,1:surf_usm_h%ns)     )
             surf_usm_h%albedo    = albedo
          ENDIF

          DO  l = 0, 3
             IF ( .NOT. ALLOCATED( surf_lsm_v(l)%albedo ) )  THEN
                ALLOCATE( surf_lsm_v(l)%albedo(0:2,1:surf_lsm_v(l)%ns) )
                surf_lsm_v(l)%albedo = albedo
             ENDIF
             IF ( .NOT. ALLOCATED( surf_usm_v(l)%albedo ) )  THEN
                ALLOCATE( surf_usm_v(l)%albedo(0:2,1:surf_usm_v(l)%ns) )
                surf_usm_v(l)%albedo = albedo
             ENDIF
          ENDDO
!
!--       Level 2 initialization of broadband albedo via given albedo_type. 
!--       Only if albedo_type is non-zero
          DO  m = 1, surf_lsm_h%ns
             IF ( surf_lsm_h%albedo_type(ind_veg_wall,m) /= 0 )                &
                surf_lsm_h%albedo(ind_veg_wall,m) =                            &
                           albedo_pars(2,surf_lsm_h%albedo_type(ind_veg_wall,m))
             IF ( surf_lsm_h%albedo_type(ind_pav_green,m) /= 0 )               &
                surf_lsm_h%albedo(ind_pav_green,m) =                           &
                           albedo_pars(2,surf_lsm_h%albedo_type(ind_pav_green,m))
             IF ( surf_lsm_h%albedo_type(ind_wat_win,m) /= 0 )                 &
                surf_lsm_h%albedo(ind_wat_win,m) =                             &
                           albedo_pars(2,surf_lsm_h%albedo_type(ind_wat_win,m))
          ENDDO
          DO  m = 1, surf_usm_h%ns
             IF ( surf_usm_h%albedo_type(ind_veg_wall,m) /= 0 )                &
                surf_usm_h%albedo(ind_veg_wall,m) =                            &
                           albedo_pars(2,surf_usm_h%albedo_type(ind_veg_wall,m))
             IF ( surf_usm_h%albedo_type(ind_pav_green,m) /= 0 )               &
                surf_usm_h%albedo(ind_pav_green,m) =                           &
                           albedo_pars(2,surf_usm_h%albedo_type(ind_pav_green,m))
             IF ( surf_usm_h%albedo_type(ind_wat_win,m) /= 0 )                 &
                surf_usm_h%albedo(ind_wat_win,m) =                             &
                           albedo_pars(2,surf_usm_h%albedo_type(ind_wat_win,m))
          ENDDO

          DO  l = 0, 3
             DO  m = 1, surf_lsm_v(l)%ns
                IF ( surf_lsm_v(l)%albedo_type(ind_veg_wall,m) /= 0 )          &
                   surf_lsm_v(l)%albedo(ind_veg_wall,m) =                      &
                        albedo_pars(2,surf_lsm_v(l)%albedo_type(ind_veg_wall,m))
                IF ( surf_lsm_v(l)%albedo_type(ind_pav_green,m) /= 0 )         &
                   surf_lsm_v(l)%albedo(ind_pav_green,m) =                     &
                        albedo_pars(2,surf_lsm_v(l)%albedo_type(ind_pav_green,m))
                IF ( surf_lsm_v(l)%albedo_type(ind_wat_win,m) /= 0 )           &
                   surf_lsm_v(l)%albedo(ind_wat_win,m) =                       &
                        albedo_pars(2,surf_lsm_v(l)%albedo_type(ind_wat_win,m))
             ENDDO
             DO  m = 1, surf_usm_v(l)%ns
                IF ( surf_usm_v(l)%albedo_type(ind_veg_wall,m) /= 0 )          &
                   surf_usm_v(l)%albedo(ind_veg_wall,m) =                      &
                        albedo_pars(2,surf_usm_v(l)%albedo_type(ind_veg_wall,m))
                IF ( surf_usm_v(l)%albedo_type(ind_pav_green,m) /= 0 )         &
                   surf_usm_v(l)%albedo(ind_pav_green,m) =                     &
                        albedo_pars(2,surf_usm_v(l)%albedo_type(ind_pav_green,m))
                IF ( surf_usm_v(l)%albedo_type(ind_wat_win,m) /= 0 )           &
                   surf_usm_v(l)%albedo(ind_wat_win,m) =                       &
                        albedo_pars(2,surf_usm_v(l)%albedo_type(ind_wat_win,m))
             ENDDO
          ENDDO

!
!--       Level 3 initialization at grid points where albedo type is zero.
!--       This case, albedo is taken from file. In case of constant radiation
!--       or clear sky, only broadband albedo is given.
          IF ( albedo_pars_f%from_file )  THEN
!
!--          Horizontal surfaces
             DO  m = 1, surf_lsm_h%ns
                i = surf_lsm_h%i(m)
                j = surf_lsm_h%j(m)
                IF ( albedo_pars_f%pars_xy(0,j,i) /= albedo_pars_f%fill )  THEN
                   IF ( surf_lsm_h%albedo_type(ind_veg_wall,m) == 0 )          &
                      surf_lsm_h%albedo(ind_veg_wall,m) = albedo_pars_f%pars_xy(0,j,i)
                   IF ( surf_lsm_h%albedo_type(ind_pav_green,m) == 0 )         &
                      surf_lsm_h%albedo(ind_pav_green,m) = albedo_pars_f%pars_xy(0,j,i)
                   IF ( surf_lsm_h%albedo_type(ind_wat_win,m) == 0 )           &
                      surf_lsm_h%albedo(ind_wat_win,m) = albedo_pars_f%pars_xy(0,j,i)
                ENDIF
             ENDDO
             DO  m = 1, surf_usm_h%ns
                i = surf_usm_h%i(m)
                j = surf_usm_h%j(m)
                IF ( albedo_pars_f%pars_xy(0,j,i) /= albedo_pars_f%fill )  THEN
                   IF ( surf_usm_h%albedo_type(ind_veg_wall,m) == 0 )          &
                      surf_usm_h%albedo(ind_veg_wall,m) = albedo_pars_f%pars_xy(0,j,i)
                   IF ( surf_usm_h%albedo_type(ind_pav_green,m) == 0 )         &
                      surf_usm_h%albedo(ind_pav_green,m) = albedo_pars_f%pars_xy(0,j,i)
                   IF ( surf_usm_h%albedo_type(ind_wat_win,m) == 0 )           &
                      surf_usm_h%albedo(ind_wat_win,m) = albedo_pars_f%pars_xy(0,j,i)
                ENDIF
             ENDDO 
!
!--          Vertical surfaces            
             DO  l = 0, 3

                ioff = surf_lsm_v(l)%ioff
                joff = surf_lsm_v(l)%joff
                DO  m = 1, surf_lsm_v(l)%ns
                   i = surf_lsm_v(l)%i(m) + ioff
                   j = surf_lsm_v(l)%j(m) + joff
                   IF ( albedo_pars_f%pars_xy(0,j,i) /= albedo_pars_f%fill )  THEN
                      IF ( surf_lsm_v(l)%albedo_type(ind_veg_wall,m) == 0 )    &
                         surf_lsm_v(l)%albedo(ind_veg_wall,m) = albedo_pars_f%pars_xy(0,j,i)
                      IF ( surf_lsm_v(l)%albedo_type(ind_pav_green,m) == 0 )   &
                         surf_lsm_v(l)%albedo(ind_pav_green,m) = albedo_pars_f%pars_xy(0,j,i)
                      IF ( surf_lsm_v(l)%albedo_type(ind_wat_win,m) == 0 )     &
                         surf_lsm_v(l)%albedo(ind_wat_win,m) = albedo_pars_f%pars_xy(0,j,i)
                   ENDIF
                ENDDO

                ioff = surf_usm_v(l)%ioff
                joff = surf_usm_v(l)%joff
                DO  m = 1, surf_usm_h%ns
                   i = surf_usm_h%i(m) + joff
                   j = surf_usm_h%j(m) + joff
                   IF ( albedo_pars_f%pars_xy(0,j,i) /= albedo_pars_f%fill )  THEN
                      IF ( surf_usm_v(l)%albedo_type(ind_veg_wall,m) == 0 )    &
                         surf_usm_v(l)%albedo(ind_veg_wall,m) = albedo_pars_f%pars_xy(0,j,i)
                      IF ( surf_usm_v(l)%albedo_type(ind_pav_green,m) == 0 )   &
                         surf_usm_v(l)%albedo(ind_pav_green,m) = albedo_pars_f%pars_xy(0,j,i)
                      IF ( surf_usm_v(l)%albedo_type(ind_wat_win,m) == 0 )     &
                         surf_lsm_v(l)%albedo(ind_wat_win,m) = albedo_pars_f%pars_xy(0,j,i)
                   ENDIF
                ENDDO
             ENDDO

          ENDIF  
!
!--    Initialization actions for RRTMG
       ELSEIF ( radiation_scheme == 'rrtmg' )  THEN
#if defined ( __rrtmg )
!
!--       Allocate albedos for short/longwave radiation, horizontal surfaces
!--       for wall/green/window (USM) or vegetation/pavement/water surfaces
!--       (LSM). 
          ALLOCATE ( surf_lsm_h%aldif(0:2,1:surf_lsm_h%ns)       )
          ALLOCATE ( surf_lsm_h%aldir(0:2,1:surf_lsm_h%ns)       )
          ALLOCATE ( surf_lsm_h%asdif(0:2,1:surf_lsm_h%ns)       )
          ALLOCATE ( surf_lsm_h%asdir(0:2,1:surf_lsm_h%ns)       )
          ALLOCATE ( surf_lsm_h%rrtm_aldif(0:2,1:surf_lsm_h%ns)  )
          ALLOCATE ( surf_lsm_h%rrtm_aldir(0:2,1:surf_lsm_h%ns)  )
          ALLOCATE ( surf_lsm_h%rrtm_asdif(0:2,1:surf_lsm_h%ns)  )
          ALLOCATE ( surf_lsm_h%rrtm_asdir(0:2,1:surf_lsm_h%ns)  )

          ALLOCATE ( surf_usm_h%aldif(0:2,1:surf_usm_h%ns)       )
          ALLOCATE ( surf_usm_h%aldir(0:2,1:surf_usm_h%ns)       )
          ALLOCATE ( surf_usm_h%asdif(0:2,1:surf_usm_h%ns)       )
          ALLOCATE ( surf_usm_h%asdir(0:2,1:surf_usm_h%ns)       )
          ALLOCATE ( surf_usm_h%rrtm_aldif(0:2,1:surf_usm_h%ns)  )
          ALLOCATE ( surf_usm_h%rrtm_aldir(0:2,1:surf_usm_h%ns)  )
          ALLOCATE ( surf_usm_h%rrtm_asdif(0:2,1:surf_usm_h%ns)  )
          ALLOCATE ( surf_usm_h%rrtm_asdir(0:2,1:surf_usm_h%ns)  )

!
!--       Allocate broadband albedo (temporary for the current radiation 
!--       implementations)
          IF ( .NOT. ALLOCATED(surf_lsm_h%albedo) )                            &
             ALLOCATE( surf_lsm_h%albedo(0:2,1:surf_lsm_h%ns)     )
          IF ( .NOT. ALLOCATED(surf_usm_h%albedo) )                            &
             ALLOCATE( surf_usm_h%albedo(0:2,1:surf_usm_h%ns)     )

!
!--       Allocate albedos for short/longwave radiation, vertical surfaces 
          DO  l = 0, 3

             ALLOCATE ( surf_lsm_v(l)%aldif(0:2,1:surf_lsm_v(l)%ns)      )
             ALLOCATE ( surf_lsm_v(l)%aldir(0:2,1:surf_lsm_v(l)%ns)      )
             ALLOCATE ( surf_lsm_v(l)%asdif(0:2,1:surf_lsm_v(l)%ns)      )
             ALLOCATE ( surf_lsm_v(l)%asdir(0:2,1:surf_lsm_v(l)%ns)      )

             ALLOCATE ( surf_lsm_v(l)%rrtm_aldif(0:2,1:surf_lsm_v(l)%ns) )
             ALLOCATE ( surf_lsm_v(l)%rrtm_aldir(0:2,1:surf_lsm_v(l)%ns) )
             ALLOCATE ( surf_lsm_v(l)%rrtm_asdif(0:2,1:surf_lsm_v(l)%ns) )
             ALLOCATE ( surf_lsm_v(l)%rrtm_asdir(0:2,1:surf_lsm_v(l)%ns) )

             ALLOCATE ( surf_usm_v(l)%aldif(0:2,1:surf_usm_v(l)%ns)      )
             ALLOCATE ( surf_usm_v(l)%aldir(0:2,1:surf_usm_v(l)%ns)      )
             ALLOCATE ( surf_usm_v(l)%asdif(0:2,1:surf_usm_v(l)%ns)      )
             ALLOCATE ( surf_usm_v(l)%asdir(0:2,1:surf_usm_v(l)%ns)      )

             ALLOCATE ( surf_usm_v(l)%rrtm_aldif(0:2,1:surf_usm_v(l)%ns) )
             ALLOCATE ( surf_usm_v(l)%rrtm_aldir(0:2,1:surf_usm_v(l)%ns) )
             ALLOCATE ( surf_usm_v(l)%rrtm_asdif(0:2,1:surf_usm_v(l)%ns) )
             ALLOCATE ( surf_usm_v(l)%rrtm_asdir(0:2,1:surf_usm_v(l)%ns) )
!
!--          Allocate broadband albedo (temporary for the current radiation
!--          implementations)
             IF ( .NOT. ALLOCATED( surf_lsm_v(l)%albedo ) )                    &
                ALLOCATE( surf_lsm_v(l)%albedo(0:2,1:surf_lsm_v(l)%ns) )
             IF ( .NOT. ALLOCATED( surf_usm_v(l)%albedo ) )                    &
                ALLOCATE( surf_usm_v(l)%albedo(0:2,1:surf_usm_v(l)%ns) )

          ENDDO
!
!--       Level 1 initialization of spectral albedos via namelist 
!--       paramters. Please note, this case all surface tiles are initialized 
!--       the same. 
          IF ( surf_lsm_h%ns > 0 )  THEN
             surf_lsm_h%aldif  = albedo_lw_dif
             surf_lsm_h%aldir  = albedo_lw_dir
             surf_lsm_h%asdif  = albedo_sw_dif
             surf_lsm_h%asdir  = albedo_sw_dir
             surf_lsm_h%albedo = albedo_sw_dif
          ENDIF
          IF ( surf_usm_h%ns > 0 )  THEN
             surf_usm_h%aldif  = albedo_lw_dif
             surf_usm_h%aldir  = albedo_lw_dir
             surf_usm_h%asdif  = albedo_sw_dif
             surf_usm_h%asdir  = albedo_sw_dir
             surf_usm_h%albedo = albedo_sw_dif
          ENDIF

          DO  l = 0, 3

             IF ( surf_lsm_v(l)%ns > 0 )  THEN
                surf_lsm_v(l)%aldif  = albedo_lw_dif
                surf_lsm_v(l)%aldir  = albedo_lw_dir
                surf_lsm_v(l)%asdif  = albedo_sw_dif
                surf_lsm_v(l)%asdir  = albedo_sw_dir
                surf_lsm_v(l)%albedo = albedo_sw_dif
             ENDIF

             IF ( surf_usm_v(l)%ns > 0 )  THEN
                surf_usm_v(l)%aldif  = albedo_lw_dif
                surf_usm_v(l)%aldir  = albedo_lw_dir
                surf_usm_v(l)%asdif  = albedo_sw_dif
                surf_usm_v(l)%asdir  = albedo_sw_dir
                surf_usm_v(l)%albedo = albedo_sw_dif
             ENDIF
          ENDDO

!
!--       Level 2 initialization of spectral albedos via albedo_type. 
!--       Please note, for natural- and urban-type surfaces, a tile approach 
!--       is applied so that the resulting albedo is calculated via the weighted 
!--       average of respective surface fractions. 
          DO  m = 1, surf_lsm_h%ns
!
!--          Spectral albedos for vegetation/pavement/water surfaces
             DO  ind_type = 0, 2
                IF ( surf_lsm_h%albedo_type(ind_type,m) /= 0 )  THEN
                   surf_lsm_h%aldif(ind_type,m) =                              &
                               albedo_pars(0,surf_lsm_h%albedo_type(ind_type,m))
                   surf_lsm_h%asdif(ind_type,m) =                              &
                               albedo_pars(1,surf_lsm_h%albedo_type(ind_type,m))
                   surf_lsm_h%aldir(ind_type,m) =                              &
                               albedo_pars(0,surf_lsm_h%albedo_type(ind_type,m))
                   surf_lsm_h%asdir(ind_type,m) =                              &
                               albedo_pars(1,surf_lsm_h%albedo_type(ind_type,m))
                   surf_lsm_h%albedo(ind_type,m) =                             &
                               albedo_pars(2,surf_lsm_h%albedo_type(ind_type,m))
                ENDIF
             ENDDO

          ENDDO

          DO  m = 1, surf_usm_h%ns
!
!--          Spectral albedos for wall/green/window surfaces
             DO  ind_type = 0, 2
                IF ( surf_usm_h%albedo_type(ind_type,m) /= 0 )  THEN
                   surf_usm_h%aldif(ind_type,m) =                              &
                               albedo_pars(0,surf_usm_h%albedo_type(ind_type,m))
                   surf_usm_h%asdif(ind_type,m) =                              &
                               albedo_pars(1,surf_usm_h%albedo_type(ind_type,m))
                   surf_usm_h%aldir(ind_type,m) =                              &
                               albedo_pars(0,surf_usm_h%albedo_type(ind_type,m))
                   surf_usm_h%asdir(ind_type,m) =                              &
                               albedo_pars(1,surf_usm_h%albedo_type(ind_type,m))
                   surf_usm_h%albedo(ind_type,m) =                             &
                               albedo_pars(2,surf_usm_h%albedo_type(ind_type,m))
                ENDIF
             ENDDO

          ENDDO

          DO l = 0, 3

             DO  m = 1, surf_lsm_v(l)%ns
!
!--             Spectral albedos for vegetation/pavement/water surfaces
                DO  ind_type = 0, 2
                   IF ( surf_lsm_v(l)%albedo_type(ind_type,m) /= 0 )  THEN
                      surf_lsm_v(l)%aldif(ind_type,m) =                        &
                            albedo_pars(0,surf_lsm_v(l)%albedo_type(ind_type,m))
                      surf_lsm_v(l)%asdif(ind_type,m) =                        &
                            albedo_pars(1,surf_lsm_v(l)%albedo_type(ind_type,m))
                      surf_lsm_v(l)%aldir(ind_type,m) =                        &
                            albedo_pars(0,surf_lsm_v(l)%albedo_type(ind_type,m))
                      surf_lsm_v(l)%asdir(ind_type,m) =                        &
                            albedo_pars(1,surf_lsm_v(l)%albedo_type(ind_type,m))
                      surf_lsm_v(l)%albedo(ind_type,m) =                       &
                            albedo_pars(2,surf_lsm_v(l)%albedo_type(ind_type,m))
                   ENDIF
                ENDDO
             ENDDO

             DO  m = 1, surf_usm_v(l)%ns
!
!--             Spectral albedos for wall/green/window surfaces
                DO  ind_type = 0, 2
                   IF ( surf_usm_v(l)%albedo_type(ind_type,m) /= 0 )  THEN
                      surf_usm_v(l)%aldif(ind_type,m) =                        &
                            albedo_pars(0,surf_usm_v(l)%albedo_type(ind_type,m))
                      surf_usm_v(l)%asdif(ind_type,m) =                        &
                            albedo_pars(1,surf_usm_v(l)%albedo_type(ind_type,m))
                      surf_usm_v(l)%aldir(ind_type,m) =                        &
                            albedo_pars(0,surf_usm_v(l)%albedo_type(ind_type,m))
                      surf_usm_v(l)%asdir(ind_type,m) =                        &
                            albedo_pars(1,surf_usm_v(l)%albedo_type(ind_type,m))
                      surf_usm_v(l)%albedo(ind_type,m) =                       &
                            albedo_pars(2,surf_usm_v(l)%albedo_type(ind_type,m))
                   ENDIF
                ENDDO

             ENDDO
          ENDDO
!
!--       Level 3 initialization at grid points where albedo type is zero.
!--       This case, spectral albedos are taken from file if available
          IF ( albedo_pars_f%from_file )  THEN
!
!--          Horizontal
             DO  m = 1, surf_lsm_h%ns
                i = surf_lsm_h%i(m)
                j = surf_lsm_h%j(m)
!
!--             Spectral albedos for vegetation/pavement/water surfaces
                DO  ind_type = 0, 2
                   IF ( surf_lsm_h%albedo_type(ind_type,m) == 0 )  THEN
                      IF ( albedo_pars_f%pars_xy(1,j,i) /= albedo_pars_f%fill )&
                         surf_lsm_h%albedo(ind_type,m) =                       &
                                                albedo_pars_f%pars_xy(1,j,i)
                      IF ( albedo_pars_f%pars_xy(1,j,i) /= albedo_pars_f%fill )&
                         surf_lsm_h%aldir(ind_type,m) =                        &
                                                albedo_pars_f%pars_xy(1,j,i)
                      IF ( albedo_pars_f%pars_xy(2,j,i) /= albedo_pars_f%fill )&
                         surf_lsm_h%aldif(ind_type,m) =                        &
                                                albedo_pars_f%pars_xy(2,j,i)
                      IF ( albedo_pars_f%pars_xy(3,j,i) /= albedo_pars_f%fill )&
                         surf_lsm_h%asdir(ind_type,m) =                        &
                                                albedo_pars_f%pars_xy(3,j,i)
                      IF ( albedo_pars_f%pars_xy(4,j,i) /= albedo_pars_f%fill )&
                         surf_lsm_h%asdif(ind_type,m) =                        &
                                                albedo_pars_f%pars_xy(4,j,i)
                   ENDIF
                ENDDO
             ENDDO

             DO  m = 1, surf_usm_h%ns
                i = surf_usm_h%i(m)
                j = surf_usm_h%j(m)
!
!--             Spectral albedos for wall/green/window surfaces
                DO  ind_type = 0, 2
                   IF ( surf_usm_h%albedo_type(ind_type,m) == 0 )  THEN
                      IF ( albedo_pars_f%pars_xy(1,j,i) /= albedo_pars_f%fill )&
                         surf_usm_h%albedo(ind_type,m) =                       &
                                                albedo_pars_f%pars_xy(1,j,i)
                      IF ( albedo_pars_f%pars_xy(1,j,i) /= albedo_pars_f%fill )&
                         surf_usm_h%aldir(ind_type,m) =                        &
                                                albedo_pars_f%pars_xy(1,j,i)
                      IF ( albedo_pars_f%pars_xy(2,j,i) /= albedo_pars_f%fill )&
                         surf_usm_h%aldif(ind_type,m) =                        &
                                                albedo_pars_f%pars_xy(2,j,i)
                      IF ( albedo_pars_f%pars_xy(3,j,i) /= albedo_pars_f%fill )&
                         surf_usm_h%asdir(ind_type,m) =                        &
                                                albedo_pars_f%pars_xy(3,j,i)
                      IF ( albedo_pars_f%pars_xy(4,j,i) /= albedo_pars_f%fill )&
                         surf_usm_h%asdif(ind_type,m) =                        &
                                                albedo_pars_f%pars_xy(4,j,i)
                   ENDIF
                ENDDO

             ENDDO
!
!--          Vertical
             DO  l = 0, 3
                ioff = surf_lsm_v(l)%ioff
                joff = surf_lsm_v(l)%joff

                DO  m = 1, surf_lsm_v(l)%ns
                   i = surf_lsm_v(l)%i(m)
                   j = surf_lsm_v(l)%j(m)
!
!--                Spectral albedos for vegetation/pavement/water surfaces
                   DO  ind_type = 0, 2
                      IF ( surf_lsm_v(l)%albedo_type(ind_type,m) == 0 )  THEN
                         IF ( albedo_pars_f%pars_xy(1,j+joff,i+ioff) /=        &
                              albedo_pars_f%fill )                             &
                            surf_lsm_v(l)%albedo(ind_type,m) =                 &
                                          albedo_pars_f%pars_xy(1,j+joff,i+ioff)
                         IF ( albedo_pars_f%pars_xy(1,j+joff,i+ioff) /=        &
                              albedo_pars_f%fill )                             &
                            surf_lsm_v(l)%aldir(ind_type,m) =                  &
                                          albedo_pars_f%pars_xy(1,j+joff,i+ioff)
                         IF ( albedo_pars_f%pars_xy(2,j+joff,i+ioff) /=        &
                              albedo_pars_f%fill )                             &
                            surf_lsm_v(l)%aldif(ind_type,m) =                  &
                                          albedo_pars_f%pars_xy(2,j+joff,i+ioff)
                         IF ( albedo_pars_f%pars_xy(3,j+joff,i+ioff) /=        &
                              albedo_pars_f%fill )                             &
                            surf_lsm_v(l)%asdir(ind_type,m) =                  &
                                          albedo_pars_f%pars_xy(3,j+joff,i+ioff)
                         IF ( albedo_pars_f%pars_xy(4,j+joff,i+ioff) /=        &
                              albedo_pars_f%fill )                             &
                            surf_lsm_v(l)%asdif(ind_type,m) =                  &
                                          albedo_pars_f%pars_xy(4,j+joff,i+ioff)
                      ENDIF
                   ENDDO
                ENDDO

                ioff = surf_usm_v(l)%ioff
                joff = surf_usm_v(l)%joff

                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
!
!--                Spectral albedos for wall/green/window surfaces
                   DO  ind_type = 0, 2
                      IF ( surf_usm_v(l)%albedo_type(ind_type,m) == 0 )  THEN
                         IF ( albedo_pars_f%pars_xy(1,j+joff,i+ioff) /=        &
                              albedo_pars_f%fill )                             &
                            surf_usm_v(l)%albedo(ind_type,m) =                 &
                                          albedo_pars_f%pars_xy(1,j+joff,i+ioff)
                         IF ( albedo_pars_f%pars_xy(1,j+joff,i+ioff) /=        &
                              albedo_pars_f%fill )                             &
                            surf_usm_v(l)%aldir(ind_type,m) =                  &
                                          albedo_pars_f%pars_xy(1,j+joff,i+ioff)
                         IF ( albedo_pars_f%pars_xy(2,j+joff,i+ioff) /=        &
                              albedo_pars_f%fill )                             &
                            surf_usm_v(l)%aldif(ind_type,m) =                  &
                                          albedo_pars_f%pars_xy(2,j+joff,i+ioff)
                         IF ( albedo_pars_f%pars_xy(3,j+joff,i+ioff) /=        &
                              albedo_pars_f%fill )                             &
                            surf_usm_v(l)%asdir(ind_type,m) =                  &
                                          albedo_pars_f%pars_xy(3,j+joff,i+ioff)
                         IF ( albedo_pars_f%pars_xy(4,j+joff,i+ioff) /=        &
                              albedo_pars_f%fill )                             &
                            surf_usm_v(l)%asdif(ind_type,m) =                  &
                                          albedo_pars_f%pars_xy(4,j+joff,i+ioff)
                      ENDIF
                   ENDDO

                ENDDO
             ENDDO

          ENDIF

!
!--       Calculate initial values of current (cosine of) the zenith angle and 
!--       whether the sun is up
          CALL calc_zenith     
!
!--       Calculate initial surface albedo for different surfaces
          IF ( .NOT. constant_albedo )  THEN
!
!--          Horizontally aligned natural and urban surfaces
             CALL calc_albedo( surf_lsm_h    )
             CALL calc_albedo( surf_usm_h    )
!
!--          Vertically aligned natural and urban surfaces
             DO  l = 0, 3
                CALL calc_albedo( surf_lsm_v(l) )
                CALL calc_albedo( surf_usm_v(l) )
             ENDDO
          ELSE
!
!--          Initialize sun-inclination independent spectral albedos
!--          Horizontal surfaces
             IF ( surf_lsm_h%ns > 0 )  THEN
                surf_lsm_h%rrtm_aldir = surf_lsm_h%aldir
                surf_lsm_h%rrtm_asdir = surf_lsm_h%asdir
                surf_lsm_h%rrtm_aldif = surf_lsm_h%aldif
                surf_lsm_h%rrtm_asdif = surf_lsm_h%asdif
             ENDIF
             IF ( surf_usm_h%ns > 0 )  THEN
                surf_usm_h%rrtm_aldir = surf_usm_h%aldir
                surf_usm_h%rrtm_asdir = surf_usm_h%asdir
                surf_usm_h%rrtm_aldif = surf_usm_h%aldif
                surf_usm_h%rrtm_asdif = surf_usm_h%asdif
             ENDIF
!
!--          Vertical surfaces
             DO  l = 0, 3
                IF ( surf_lsm_v(l)%ns > 0 )  THEN
                   surf_lsm_v(l)%rrtm_aldir = surf_lsm_v(l)%aldir
                   surf_lsm_v(l)%rrtm_asdir = surf_lsm_v(l)%asdir
                   surf_lsm_v(l)%rrtm_aldif = surf_lsm_v(l)%aldif
                   surf_lsm_v(l)%rrtm_asdif = surf_lsm_v(l)%asdif
                ENDIF
                IF ( surf_usm_v(l)%ns > 0 )  THEN
                   surf_usm_v(l)%rrtm_aldir = surf_usm_v(l)%aldir
                   surf_usm_v(l)%rrtm_asdir = surf_usm_v(l)%asdir
                   surf_usm_v(l)%rrtm_aldif = surf_usm_v(l)%aldif
                   surf_usm_v(l)%rrtm_asdif = surf_usm_v(l)%asdif
                ENDIF
             ENDDO

          ENDIF

!
!--       Allocate 3d arrays of radiative fluxes and heating rates
          IF ( .NOT. ALLOCATED ( rad_sw_in ) )  THEN
             ALLOCATE ( rad_sw_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_sw_in = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_in_av ) )  THEN
             ALLOCATE ( rad_sw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_out ) )  THEN
             ALLOCATE ( rad_sw_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_sw_out = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_out_av ) )  THEN
             ALLOCATE ( rad_sw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_hr ) )  THEN
             ALLOCATE ( rad_sw_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_sw_hr = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_hr_av ) )  THEN
             ALLOCATE ( rad_sw_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_sw_hr_av = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_cs_hr ) )  THEN
             ALLOCATE ( rad_sw_cs_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_sw_cs_hr = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_cs_hr_av ) )  THEN
             ALLOCATE ( rad_sw_cs_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_sw_cs_hr_av = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_in ) )  THEN
             ALLOCATE ( rad_lw_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_lw_in     = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_in_av ) )  THEN
             ALLOCATE ( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_out ) )  THEN
             ALLOCATE ( rad_lw_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
            rad_lw_out    = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_out_av ) )  THEN
             ALLOCATE ( rad_lw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_hr ) )  THEN
             ALLOCATE ( rad_lw_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_lw_hr = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_hr_av ) )  THEN
             ALLOCATE ( rad_lw_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_lw_hr_av = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_cs_hr ) )  THEN
             ALLOCATE ( rad_lw_cs_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_lw_cs_hr = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_cs_hr_av ) )  THEN
             ALLOCATE ( rad_lw_cs_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_lw_cs_hr_av = 0.0_wp
          ENDIF

          ALLOCATE ( rad_sw_cs_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ALLOCATE ( rad_sw_cs_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          rad_sw_cs_in  = 0.0_wp
          rad_sw_cs_out = 0.0_wp

          ALLOCATE ( rad_lw_cs_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ALLOCATE ( rad_lw_cs_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          rad_lw_cs_in  = 0.0_wp
          rad_lw_cs_out = 0.0_wp

!
!--       Allocate 1-element array for surface temperature 
!--       (RRTMG anticipates an array as passed argument).
          ALLOCATE ( rrtm_tsfc(1) )
!
!--       Allocate surface emissivity. 
!--       Values will be given directly before calling rrtm_lw.
          ALLOCATE ( rrtm_emis(0:0,1:nbndlw+1) )

!
!--       Initialize RRTMG
          IF ( lw_radiation )  CALL rrtmg_lw_ini ( cp )
          IF ( sw_radiation )  CALL rrtmg_sw_ini ( cp )

!
!--       Set input files for RRTMG
          INQUIRE(FILE="RAD_SND_DATA", EXIST=snd_exists) 
          IF ( .NOT. snd_exists )  THEN
             rrtm_input_file = "rrtmg_lw.nc"
          ENDIF

!
!--       Read vertical layers for RRTMG from sounding data
!--       The routine provides nzt_rad, hyp_snd(1:nzt_rad),
!--       t_snd(nzt+2:nzt_rad), rrtm_play(1:nzt_rad), rrtm_plev(1_nzt_rad+1), 
!--       rrtm_tlay(nzt+2:nzt_rad), rrtm_tlev(nzt+2:nzt_rad+1)
          CALL read_sounding_data

!
!--       Read trace gas profiles from file. This routine provides
!--       the rrtm_ arrays (1:nzt_rad+1)
          CALL read_trace_gas_data
#endif
       ENDIF

!
!--    Perform user actions if required
       CALL user_init_radiation

!
!--    Calculate radiative fluxes at model start
       SELECT CASE ( TRIM( radiation_scheme ) )

          CASE ( 'rrtmg' )
             CALL radiation_rrtmg

          CASE ( 'clear-sky' )
             CALL radiation_clearsky

          CASE ( 'constant' )
             CALL radiation_constant

          CASE DEFAULT

       END SELECT

       RETURN

    END SUBROUTINE radiation_init


!------------------------------------------------------------------------------!
! Description:
! ------------
!> A simple clear sky radiation model
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_clearsky


       IMPLICIT NONE

       INTEGER(iwp) ::  l         !< running index for surface orientation

       REAL(wp)     ::  exn       !< Exner functions at surface
       REAL(wp)     ::  exn1      !< Exner functions at first grid level or at urban layer top 
       REAL(wp)     ::  pt1       !< potential temperature at first grid level or mean value at urban layer top 
       REAL(wp)     ::  pt1_l     !< potential temperature at first grid level or mean value at urban layer top at local subdomain 
       REAL(wp)     ::  ql1       !< liquid water mixing ratio at first grid level or mean value at urban layer top 
       REAL(wp)     ::  ql1_l     !< liquid water mixing ratio at first grid level or mean value at urban layer top at local subdomain 

       TYPE(surf_type), POINTER ::  surf !< pointer on respective surface type, used to generalize routine   

!
!--    Calculate current zenith angle
       CALL calc_zenith

!
!--    Calculate sky transmissivity
       sky_trans = 0.6_wp + 0.2_wp * zenith(0)

!
!--    Calculate value of the Exner function at model surface
       exn = (surface_pressure / 1000.0_wp )**0.286_wp
!
!--    In case averaged radiation is used, calculate mean temperature and 
!--    liquid water mixing ratio at the urban-layer top.
       IF ( average_radiation ) THEN   
          pt1   = 0.0_wp
          IF ( cloud_physics )  ql1   = 0.0_wp

          pt1_l = SUM( pt(nzut,nys:nyn,nxl:nxr) )
          IF ( cloud_physics )  ql1_l = SUM( ql(nzut,nys:nyn,nxl:nxr) )

#if defined( __parallel )      
          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_ALLREDUCE( pt1_l, pt1, 1, MPI_REAL, MPI_SUM, comm2d, ierr )
          IF ( cloud_physics )                                                 &
             CALL MPI_ALLREDUCE( ql1_l, ql1, 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
          pt1 = pt1_l 
          IF ( cloud_physics )  ql1 = ql1_l
#endif
 
          exn1 = ( hyp(nzut) / 100000.0_wp )**0.286_wp
          IF ( cloud_physics )  pt1 = pt1 + l_d_cp / exn1 * ql1
!
!--       Finally, divide by number of grid points
          pt1 = pt1 / REAL( ( nx + 1 ) * ( ny + 1 ), KIND=wp )
       ENDIF
!
!--    Call clear-sky calculation for each surface orientation. 
!--    First, horizontal surfaces
       surf => surf_lsm_h
       CALL radiation_clearsky_surf
       surf => surf_usm_h
       CALL radiation_clearsky_surf
!
!--    Vertical surfaces
       DO  l = 0, 3
          surf => surf_lsm_v(l)
          CALL radiation_clearsky_surf
          surf => surf_usm_v(l)
          CALL radiation_clearsky_surf
       ENDDO

       CONTAINS

          SUBROUTINE radiation_clearsky_surf

             IMPLICIT NONE

             INTEGER(iwp) ::  i         !< index x-direction
             INTEGER(iwp) ::  j         !< index y-direction
             INTEGER(iwp) ::  k         !< index z-direction
             INTEGER(iwp) ::  m         !< running index for surface elements

             IF ( surf%ns < 1 )  RETURN

!
!--          Calculate radiation fluxes and net radiation (rad_net) assuming 
!--          homogeneous urban radiation conditions. 
             IF ( average_radiation ) THEN       

                k = nzut

                exn1 = ( hyp(k+1) / 100000.0_wp )**0.286_wp

                surf%rad_sw_in  = solar_constant * sky_trans * zenith(0)
                surf%rad_sw_out = albedo_urb * surf%rad_sw_in
                
                surf%rad_lw_in  = 0.8_wp * sigma_sb * (pt1 * exn1)**4

                surf%rad_lw_out = emissivity_urb * sigma_sb * (t_rad_urb)**4   &
                                    + (1.0_wp - emissivity_urb) * surf%rad_lw_in

                surf%rad_net = surf%rad_sw_in - surf%rad_sw_out                &
                             + surf%rad_lw_in - surf%rad_lw_out

                surf%rad_lw_out_change_0 = 3.0_wp * emissivity_urb * sigma_sb  &
                                           * (t_rad_urb)**3

!
!--          Calculate radiation fluxes and net radiation (rad_net) for each surface 
!--          element.
             ELSE

                DO  m = 1, surf%ns
                   i = surf%i(m)
                   j = surf%j(m)
                   k = surf%k(m)

                   exn1 = (hyp(k) / 100000.0_wp )**0.286_wp

                   surf%rad_sw_in(m) = solar_constant * sky_trans * zenith(0)

!
!--                Weighted average according to surface fraction. 
!--                ATTENTION: when radiation interactions are switched on the 
!--                calculated fluxes below are not actually used as they are 
!--                overwritten in radiation_interaction.
                   surf%rad_sw_out(m) = ( surf%frac(ind_veg_wall,m)  *         &
                                          surf%albedo(ind_veg_wall,m)          &
                                        + surf%frac(ind_pav_green,m) *         &
                                          surf%albedo(ind_pav_green,m)         &
                                        + surf%frac(ind_wat_win,m)   *         &
                                          surf%albedo(ind_wat_win,m) )         &
                                        * surf%rad_sw_in(m)

                   surf%rad_lw_out(m) = ( surf%frac(ind_veg_wall,m)  *         &
                                          surf%emissivity(ind_veg_wall,m)      &
                                        + surf%frac(ind_pav_green,m) *         &
                                          surf%emissivity(ind_pav_green,m)     &
                                        + surf%frac(ind_wat_win,m)   *         &
                                          surf%emissivity(ind_wat_win,m)       &
                                        )                                      &
                                        * sigma_sb                             &
                                        * ( surf%pt_surface(m) * exn )**4

                   surf%rad_lw_out_change_0(m) =                               &
                                      ( surf%frac(ind_veg_wall,m)  *           &
                                        surf%emissivity(ind_veg_wall,m)        &
                                      + surf%frac(ind_pav_green,m) *           &
                                        surf%emissivity(ind_pav_green,m)       &
                                      + surf%frac(ind_wat_win,m)   *           &
                                        surf%emissivity(ind_wat_win,m)         &
                                      ) * 3.0_wp * sigma_sb                    &
                                      * ( surf%pt_surface(m) * exn )** 3


                   IF ( cloud_physics )  THEN
                      pt1 = pt(k,j,i) + l_d_cp / exn1 * ql(k,j,i)
                      surf%rad_lw_in(m)  = 0.8_wp * sigma_sb * (pt1 * exn1)**4
                   ELSE
                      surf%rad_lw_in(m)  = 0.8_wp * sigma_sb * (pt(k,j,i) * exn1)**4
                   ENDIF

                   surf%rad_net(m) = surf%rad_sw_in(m) - surf%rad_sw_out(m)    &
                                   + surf%rad_lw_in(m) - surf%rad_lw_out(m)

                ENDDO

             ENDIF

!
!--          Fill out values in radiation arrays
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                rad_sw_in(0,j,i) = surf%rad_sw_in(m)
                rad_sw_out(0,j,i) = surf%rad_sw_out(m)
                rad_lw_in(0,j,i) = surf%rad_lw_in(m)
                rad_lw_out(0,j,i) = surf%rad_lw_out(m)
             ENDDO
 
          END SUBROUTINE radiation_clearsky_surf

    END SUBROUTINE radiation_clearsky


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This scheme keeps the prescribed net radiation constant during the run
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_constant


       IMPLICIT NONE

       INTEGER(iwp) ::  l         !< running index for surface orientation

       REAL(wp)     ::  exn       !< Exner functions at surface
       REAL(wp)     ::  exn1      !< Exner functions at first grid level
       REAL(wp)     ::  pt1       !< potential temperature at first grid level or mean value at urban layer top 
       REAL(wp)     ::  pt1_l     !< potential temperature at first grid level or mean value at urban layer top at local subdomain 
       REAL(wp)     ::  ql1       !< liquid water mixing ratio at first grid level or mean value at urban layer top 
       REAL(wp)     ::  ql1_l     !< liquid water mixing ratio at first grid level or mean value at urban layer top at local subdomain 

       TYPE(surf_type), POINTER ::  surf !< pointer on respective surface type, used to generalize routine   

!
!--    Calculate value of the Exner function
       exn = (surface_pressure / 1000.0_wp )**0.286_wp
!
!--    In case averaged radiation is used, calculate mean temperature and 
!--    liquid water mixing ratio at the urban-layer top.
       IF ( average_radiation ) THEN   
          pt1   = 0.0_wp
          IF ( cloud_physics )  ql1   = 0.0_wp

          pt1_l = SUM( pt(nzut,nys:nyn,nxl:nxr) )
          IF ( cloud_physics )  ql1_l = SUM( ql(nzut,nys:nyn,nxl:nxr) )

#if defined( __parallel )      
          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_ALLREDUCE( pt1_l, pt1, 1, MPI_REAL, MPI_SUM, comm2d, ierr )
          IF ( cloud_physics )                                                 &
             CALL MPI_ALLREDUCE( ql1_l, ql1, 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
          pt1 = pt1_l
          IF ( cloud_physics )  ql1 = ql1_l
#endif
          IF ( cloud_physics )  pt1 = pt1 + l_d_cp / exn1 * ql1
!
!--       Finally, divide by number of grid points
          pt1 = pt1 / REAL( ( nx + 1 ) * ( ny + 1 ), KIND=wp )
       ENDIF

!
!--    First, horizontal surfaces
       surf => surf_lsm_h
       CALL radiation_constant_surf
       surf => surf_usm_h
       CALL radiation_constant_surf
!
!--    Vertical surfaces
       DO  l = 0, 3
          surf => surf_lsm_v(l)
          CALL radiation_constant_surf
          surf => surf_usm_v(l)
          CALL radiation_constant_surf
       ENDDO

       CONTAINS

          SUBROUTINE radiation_constant_surf

             IMPLICIT NONE

             INTEGER(iwp) ::  i         !< index x-direction
             INTEGER(iwp) ::  ioff      !< offset between surface element and adjacent grid point along x
             INTEGER(iwp) ::  j         !< index y-direction
             INTEGER(iwp) ::  joff      !< offset between surface element and adjacent grid point along y
             INTEGER(iwp) ::  k         !< index z-direction
             INTEGER(iwp) ::  koff      !< offset between surface element and adjacent grid point along z
             INTEGER(iwp) ::  m         !< running index for surface elements

             IF ( surf%ns < 1 )  RETURN

!--          Calculate homogenoeus urban radiation fluxes
             IF ( average_radiation ) THEN

                ! set height above canopy
                k = nzut

                surf%rad_net = net_radiation
! MS: Wyh k + 1 ?
                exn1 = (hyp(k+1) / 100000.0_wp )**0.286_wp

                surf%rad_lw_in  = 0.8_wp * sigma_sb * (pt1 * exn1)**4

                surf%rad_lw_out = emissivity_urb * sigma_sb * (t_rad_urb)**4   &
                                    + ( 10.0_wp - emissivity_urb )             & ! shouldn't be this a bulk value -- emissivity_urb?
                                    * surf%rad_lw_in

                surf%rad_lw_out_change_0 = 3.0_wp * emissivity_urb * sigma_sb  &
                                           * t_rad_urb**3

                surf%rad_sw_in = ( surf%rad_net - surf%rad_lw_in               &
                                     + surf%rad_lw_out )                       &
                                     / ( 1.0_wp - albedo_urb )

                surf%rad_sw_out =  albedo_urb * surf%rad_sw_in

!
!--          Calculate radiation fluxes for each surface element
             ELSE
!
!--             Determine index offset between surface element and adjacent 
!--             atmospheric grid point
                ioff = surf%ioff
                joff = surf%joff
                koff = surf%koff

!
!--             Prescribe net radiation and estimate the remaining radiative fluxes
                DO  m = 1, surf%ns
                   i = surf%i(m)
                   j = surf%j(m)
                   k = surf%k(m)

                   surf%rad_net(m) = net_radiation

                   exn1 = (hyp(k) / 100000.0_wp )**0.286_wp

                   IF ( cloud_physics )  THEN
                      pt1 = pt(k,j,i) + l_d_cp / exn1 * ql(k,j,i)
                      surf%rad_lw_in(m)  = 0.8_wp * sigma_sb * (pt1 * exn1)**4
                   ELSE
                      surf%rad_lw_in(m)  = 0.8_wp * sigma_sb *                 &
                                             ( pt(k,j,i) * exn1 )**4
                   ENDIF

!
!--                Weighted average according to surface fraction. 
                   surf%rad_lw_out(m) = ( surf%frac(ind_veg_wall,m)  *         &
                                          surf%emissivity(ind_veg_wall,m)      &
                                        + surf%frac(ind_pav_green,m) *         &
                                          surf%emissivity(ind_pav_green,m)     &
                                        + surf%frac(ind_wat_win,m)   *         &
                                          surf%emissivity(ind_wat_win,m)       &
                                        )                                      &
                                      * sigma_sb                               &
                                      * ( surf%pt_surface(m) * exn )**4

                   surf%rad_sw_in(m) = ( surf%rad_net(m) - surf%rad_lw_in(m)   &
                                       + surf%rad_lw_out(m) )                  &
                                       / ( 1.0_wp -                            &
                                          ( surf%frac(ind_veg_wall,m)  *       &
                                            surf%albedo(ind_veg_wall,m)        &
                                         +  surf%frac(ind_pav_green,m) *       &
                                            surf%albedo(ind_pav_green,m)       &
                                         +  surf%frac(ind_wat_win,m)   *       &
                                            surf%albedo(ind_wat_win,m) )       &
                                         )

                   surf%rad_sw_out(m) = ( surf%frac(ind_veg_wall,m)  *         &
                                          surf%albedo(ind_veg_wall,m)          &
                                        + surf%frac(ind_pav_green,m) *         &
                                          surf%albedo(ind_pav_green,m)         &
                                        + surf%frac(ind_wat_win,m)   *         &
                                          surf%albedo(ind_wat_win,m) )         &
                                      * surf%rad_sw_in(m)

                ENDDO

             ENDIF

!
!--          Fill out values in radiation arrays
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                rad_sw_in(0,j,i) = surf%rad_sw_in(m)
                rad_sw_out(0,j,i) = surf%rad_sw_out(m)
                rad_lw_in(0,j,i) = surf%rad_lw_in(m)
                rad_lw_out(0,j,i) = surf%rad_lw_out(m)
             ENDDO

          END SUBROUTINE radiation_constant_surf
          

    END SUBROUTINE radiation_constant

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for radiation model
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_header ( io )


       IMPLICIT NONE
 
       INTEGER(iwp), INTENT(IN) ::  io            !< Unit of the output file
    

       
!
!--    Write radiation model header
       WRITE( io, 3 )

       IF ( radiation_scheme == "constant" )  THEN
          WRITE( io, 4 ) net_radiation
       ELSEIF ( radiation_scheme == "clear-sky" )  THEN
          WRITE( io, 5 )
       ELSEIF ( radiation_scheme == "rrtmg" )  THEN
          WRITE( io, 6 )
          IF ( .NOT. lw_radiation )  WRITE( io, 10 )
          IF ( .NOT. sw_radiation )  WRITE( io, 11 )
       ENDIF 

       IF ( albedo_type_f%from_file  .OR.  vegetation_type_f%from_file  .OR.   &
            pavement_type_f%from_file  .OR.  water_type_f%from_file  .OR.      &
            building_type_f%from_file )  THEN
             WRITE( io, 13 )
       ELSE  
          IF ( albedo_type == 0 )  THEN
             WRITE( io, 7 ) albedo
          ELSE
             WRITE( io, 8 ) TRIM( albedo_type_name(albedo_type) )
          ENDIF
       ENDIF
       IF ( constant_albedo )  THEN
          WRITE( io, 9 )
       ENDIF
       
       WRITE( io, 12 ) dt_radiation
 

 3 FORMAT (//' Radiation model information:'/                                  &
              ' ----------------------------'/)
 4 FORMAT ('    --> Using constant net radiation: net_radiation = ', F6.2,     &
           // 'W/m**2')
 5 FORMAT ('    --> Simple radiation scheme for clear sky is used (no clouds,',&
                   ' default)')
 6 FORMAT ('    --> RRTMG scheme is used')
 7 FORMAT (/'    User-specific surface albedo: albedo =', F6.3)
 8 FORMAT (/'    Albedo is set for land surface type: ', A)
 9 FORMAT (/'    --> Albedo is fixed during the run')
10 FORMAT (/'    --> Longwave radiation is disabled')
11 FORMAT (/'    --> Shortwave radiation is disabled.')
12 FORMAT  ('    Timestep: dt_radiation = ', F6.2, '  s')
13 FORMAT (/'    Albedo is set individually for each xy-location, according '  &
                 'to given surface type.')


    END SUBROUTINE radiation_header
   

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &radiation_parameters for radiation model
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_parin


       IMPLICIT NONE

       CHARACTER (LEN=80) ::  line  !< dummy string that contains the current line of the parameter file 
       
       NAMELIST /radiation_par/   albedo, albedo_type, albedo_lw_dir,          &
                                  albedo_lw_dif, albedo_sw_dir, albedo_sw_dif, &
                                  constant_albedo, dt_radiation, emissivity,   &
                                  lw_radiation, net_radiation,                 &
                                  radiation_scheme, skip_time_do_radiation,    &
                                  sw_radiation, unscheduled_radiation_calls,   &
                                  split_diffusion_radiation,                   &
                                  max_raytracing_dist, min_irrf_value,         &
                                  nrefsteps, mrt_factors, rma_lad_raytrace,    &
                                  dist_max_svf,                                &
                                  surface_reflections, svfnorm_report_thresh,  &
                                  radiation_interactions_on
   
       NAMELIST /radiation_parameters/   albedo, albedo_type, albedo_lw_dir,   &
                                  albedo_lw_dif, albedo_sw_dir, albedo_sw_dif, &
                                  constant_albedo, dt_radiation, emissivity,   &
                                  lw_radiation, net_radiation,                 &
                                  radiation_scheme, skip_time_do_radiation,    &
                                  sw_radiation, unscheduled_radiation_calls,   &
                                  split_diffusion_radiation,                   &
                                  max_raytracing_dist, min_irrf_value,         &
                                  nrefsteps, mrt_factors, rma_lad_raytrace,    &
                                  dist_max_svf,                                &
                                  surface_reflections, svfnorm_report_thresh,  &
                                  radiation_interactions_on
   
       line = ' '
       
!
!--    Try to find radiation model namelist
       REWIND ( 11 )
       line = ' '
       DO   WHILE ( INDEX( line, '&radiation_parameters' ) == 0 )
          READ ( 11, '(A)', END=10 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read user-defined namelist
       READ ( 11, radiation_parameters )

!
!--    Set flag that indicates that the radiation model is switched on
       radiation = .TRUE.
       
       GOTO 12
!
!--    Try to find old namelist
 10    REWIND ( 11 )
       line = ' '
       DO   WHILE ( INDEX( line, '&radiation_par' ) == 0 )
          READ ( 11, '(A)', END=12 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read user-defined namelist
       READ ( 11, radiation_par )
       
       message_string = 'namelist radiation_par is deprecated and will be ' // &
                     'removed in near future. Please use namelist ' //         &
                     'radiation_parameters instead'
       CALL message( 'radiation_parin', 'PA0487', 0, 1, 0, 6, 0 )

       
!
!--    Set flag that indicates that the radiation model is switched on
       radiation = .TRUE.

       IF ( .NOT.  radiation_interactions_on  .AND.  surface_reflections )  THEN
          message_string = 'surface_reflections is allowed only when '      // &
               'radiation_interactions_on is set to TRUE'
          CALL message( 'radiation_parin', 'PA0293',1, 2, 0, 6, 0 )
       ENDIF

 12    CONTINUE
       
    END SUBROUTINE radiation_parin


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Implementation of the RRTMG radiation_scheme
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_rrtmg

       USE indices,                                                            &
           ONLY:  nbgp

       USE particle_attributes,                                                &
           ONLY:  grid_particles, number_of_particles, particles,              &
                  particle_advection_start, prt_count

       IMPLICIT NONE

#if defined ( __rrtmg )

       INTEGER(iwp) ::  i, j, k, l, m, n !< loop indices
       INTEGER(iwp) ::  k_topo     !< topography top index

       REAL(wp)     ::  nc_rad, &    !< number concentration of cloud droplets
                        s_r2,   &    !< weighted sum over all droplets with r^2
                        s_r3         !< weighted sum over all droplets with r^3

       REAL(wp), DIMENSION(0:nzt+1) :: pt_av, q_av, ql_av
!
!--    Just dummy arguments
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: rrtm_lw_taucld_dum,          &
                                                  rrtm_lw_tauaer_dum,          &
                                                  rrtm_sw_taucld_dum,          &
                                                  rrtm_sw_ssacld_dum,          &
                                                  rrtm_sw_asmcld_dum,          &
                                                  rrtm_sw_fsfcld_dum,          &
                                                  rrtm_sw_tauaer_dum,          &
                                                  rrtm_sw_ssaaer_dum,          &
                                                  rrtm_sw_asmaer_dum,          &
                                                  rrtm_sw_ecaer_dum

!
!--    Calculate current (cosine of) zenith angle and whether the sun is up
       CALL calc_zenith     
!
!--    Calculate surface albedo. In case average radiation is applied, 
!--    this is not required.
       IF ( .NOT. constant_albedo )  THEN
!
!--       Horizontally aligned default, natural and urban surfaces
          CALL calc_albedo( surf_lsm_h    )
          CALL calc_albedo( surf_usm_h    )
!
!--       Vertically aligned default, natural and urban surfaces
          DO  l = 0, 3
             CALL calc_albedo( surf_lsm_v(l) )
             CALL calc_albedo( surf_usm_v(l) )
          ENDDO
       ENDIF

!
!--    Prepare input data for RRTMG

!
!--    In case of large scale forcing with surface data, calculate new pressure
!--    profile. nzt_rad might be modified by these calls and all required arrays
!--    will then be re-allocated
       IF ( large_scale_forcing  .AND.  lsf_surf )  THEN
          CALL read_sounding_data
          CALL read_trace_gas_data
       ENDIF


       IF ( average_radiation ) THEN

          rrtm_asdir(1)  = albedo_urb
          rrtm_asdif(1)  = albedo_urb
          rrtm_aldir(1)  = albedo_urb
          rrtm_aldif(1)  = albedo_urb

          rrtm_emis = emissivity_urb
!
!--       Calculate mean pt profile. Actually, only one height level is required.
          CALL calc_mean_profile( pt, 4 )
          pt_av = hom(:, 1, 4, 0)

!
!--       Prepare profiles of temperature and H2O volume mixing ratio
          rrtm_tlev(0,nzb+1) = t_rad_urb

          IF ( cloud_physics )  THEN
             CALL calc_mean_profile( q, 41 )
             ! average  q is now in hom(:, 1, 41, 0)
             q_av = hom(:, 1, 41, 0)
             CALL calc_mean_profile( ql, 54 )
             ! average ql is now in hom(:, 1, 54, 0)
             ql_av = hom(:, 1, 54, 0)
             
             DO k = nzb+1, nzt+1
                rrtm_tlay(0,k) = pt_av(k) * ( (hyp(k) ) / 100000._wp       &
                                 )**.286_wp + l_d_cp * ql_av(k)
                rrtm_h2ovmr(0,k) = mol_mass_air_d_wv * (q_av(k) - ql_av(k))
             ENDDO
          ELSE
             DO k = nzb+1, nzt+1
                rrtm_tlay(0,k) = pt_av(k) * ( (hyp(k) ) / 100000._wp       &
                                 )**.286_wp
                rrtm_h2ovmr(0,k) = 0._wp
              ENDDO
          ENDIF

!
!--       Avoid temperature/humidity jumps at the top of the LES domain by 
!--       linear interpolation from nzt+2 to nzt+7
          DO k = nzt+2, nzt+7
             rrtm_tlay(0,k) = rrtm_tlay(0,nzt+1)                            &
                           + ( rrtm_tlay(0,nzt+8) - rrtm_tlay(0,nzt+1) )    &
                           / ( rrtm_play(0,nzt+8) - rrtm_play(0,nzt+1) )    &
                           * ( rrtm_play(0,k) - rrtm_play(0,nzt+1) )

             rrtm_h2ovmr(0,k) = rrtm_h2ovmr(0,nzt+1)                        &
                           + ( rrtm_h2ovmr(0,nzt+8) - rrtm_h2ovmr(0,nzt+1) )&
                           / ( rrtm_play(0,nzt+8)   - rrtm_play(0,nzt+1)   )&
                           * ( rrtm_play(0,k) - rrtm_play(0,nzt+1) )

          ENDDO

!--       Linear interpolate to zw grid
          DO k = nzb+2, nzt+8
             rrtm_tlev(0,k)   = rrtm_tlay(0,k-1) + (rrtm_tlay(0,k) -        &
                                rrtm_tlay(0,k-1))                           &
                                / ( rrtm_play(0,k) - rrtm_play(0,k-1) )     &
                                * ( rrtm_plev(0,k) - rrtm_play(0,k-1) )
          ENDDO


!
!--       Calculate liquid water path and cloud fraction for each column.
!--       Note that LWP is required in g/m² instead of kg/kg m.
          rrtm_cldfr  = 0.0_wp
          rrtm_reliq  = 0.0_wp
          rrtm_cliqwp = 0.0_wp
          rrtm_icld   = 0

          IF ( cloud_physics )  THEN
             DO k = nzb+1, nzt+1
                rrtm_cliqwp(0,k) =  ql_av(k) * 1000._wp *                  &
                                    (rrtm_plev(0,k) - rrtm_plev(0,k+1))     &
                                    * 100._wp / g 

                IF ( rrtm_cliqwp(0,k) > 0._wp )  THEN
                   rrtm_cldfr(0,k) = 1._wp
                   IF ( rrtm_icld == 0 )  rrtm_icld = 1

!
!--                Calculate cloud droplet effective radius
                   IF ( cloud_physics )  THEN
                      rrtm_reliq(0,k) = 1.0E6_wp * ( 3._wp * ql_av(k)      &
                                        * rho_surface                       &
                                        / ( 4._wp * pi * nc_const * rho_l )&
                                        )**.33333333333333_wp              &
                                        * EXP( LOG( sigma_gc )**2 )

                   ENDIF

!
!--                Limit effective radius
                   IF ( rrtm_reliq(0,k) > 0.0_wp )  THEN
                      rrtm_reliq(0,k) = MAX(rrtm_reliq(0,k),2.5_wp)
                      rrtm_reliq(0,k) = MIN(rrtm_reliq(0,k),60.0_wp)
                   ENDIF
                ENDIF
             ENDDO
          ENDIF

!
!--       Set surface temperature
          rrtm_tsfc = t_rad_urb

          IF ( lw_radiation )  THEN
             CALL rrtmg_lw( 1, nzt_rad      , rrtm_icld    , rrtm_idrv      ,&
             rrtm_play       , rrtm_plev    , rrtm_tlay    , rrtm_tlev      ,&
             rrtm_tsfc       , rrtm_h2ovmr  , rrtm_o3vmr   , rrtm_co2vmr    ,&
             rrtm_ch4vmr     , rrtm_n2ovmr  , rrtm_o2vmr   , rrtm_cfc11vmr  ,&
             rrtm_cfc12vmr   , rrtm_cfc22vmr, rrtm_ccl4vmr , rrtm_emis      ,&
             rrtm_inflglw    , rrtm_iceflglw, rrtm_liqflglw, rrtm_cldfr     ,&
             rrtm_lw_taucld  , rrtm_cicewp  , rrtm_cliqwp  , rrtm_reice     ,& 
             rrtm_reliq      , rrtm_lw_tauaer,                               &
             rrtm_lwuflx     , rrtm_lwdflx  , rrtm_lwhr  ,                   &
             rrtm_lwuflxc    , rrtm_lwdflxc , rrtm_lwhrc ,                   &
             rrtm_lwuflx_dt  ,  rrtm_lwuflxc_dt )

!
!--          Save fluxes
             DO k = nzb, nzt+1
                rad_lw_in(k,:,:)  = rrtm_lwdflx(0,k)
                rad_lw_out(k,:,:) = rrtm_lwuflx(0,k)
             ENDDO

!
!--          Save heating rates (convert from K/d to K/h)
             DO k = nzb+1, nzt+1
                rad_lw_hr(k,:,:)     = rrtm_lwhr(0,k)  * d_hours_day
                rad_lw_cs_hr(k,:,:)  = rrtm_lwhrc(0,k) * d_hours_day
             ENDDO

!
!--          Save surface radiative fluxes and change in LW heating rate 
!--          onto respective surface elements
!--          Horizontal surfaces
             IF ( surf_lsm_h%ns > 0 )  THEN
                surf_lsm_h%rad_lw_in           = rrtm_lwdflx(0,nzb)
                surf_lsm_h%rad_lw_out          = rrtm_lwuflx(0,nzb)
                surf_lsm_h%rad_lw_out_change_0 = rrtm_lwuflx_dt(0,nzb)
             ENDIF             
             IF ( surf_usm_h%ns > 0 )  THEN
                surf_usm_h%rad_lw_in           = rrtm_lwdflx(0,nzb)
                surf_usm_h%rad_lw_out          = rrtm_lwuflx(0,nzb)
                surf_usm_h%rad_lw_out_change_0 = rrtm_lwuflx_dt(0,nzb)
             ENDIF
!
!--          Vertical surfaces. 
             DO  l = 0, 3
                IF ( surf_lsm_v(l)%ns > 0 )  THEN
                   surf_lsm_v(l)%rad_lw_in           = rrtm_lwdflx(0,nzb)
                   surf_lsm_v(l)%rad_lw_out          = rrtm_lwuflx(0,nzb)
                   surf_lsm_v(l)%rad_lw_out_change_0 = rrtm_lwuflx_dt(0,nzb)
                ENDIF
                IF ( surf_usm_v(l)%ns > 0 )  THEN
                   surf_usm_v(l)%rad_lw_in           = rrtm_lwdflx(0,nzb)
                   surf_usm_v(l)%rad_lw_out          = rrtm_lwuflx(0,nzb)
                   surf_usm_v(l)%rad_lw_out_change_0 = rrtm_lwuflx_dt(0,nzb)
                ENDIF
             ENDDO

          ENDIF

          IF ( sw_radiation .AND. sun_up )  THEN
             CALL rrtmg_sw( 1, nzt_rad      , rrtm_icld  , rrtm_iaer        ,&
             rrtm_play       , rrtm_plev    , rrtm_tlay  , rrtm_tlev        ,&
             rrtm_tsfc       , rrtm_h2ovmr  , rrtm_o3vmr , rrtm_co2vmr      ,&
             rrtm_ch4vmr     , rrtm_n2ovmr  , rrtm_o2vmr , rrtm_asdir       ,&
             rrtm_asdif      , rrtm_aldir   , rrtm_aldif , zenith,           &
             0.0_wp          , day_of_year  , solar_constant,   rrtm_inflgsw,&
             rrtm_iceflgsw   , rrtm_liqflgsw, rrtm_cldfr , rrtm_sw_taucld   ,&
             rrtm_sw_ssacld  , rrtm_sw_asmcld, rrtm_sw_fsfcld, rrtm_cicewp  ,&
             rrtm_cliqwp     , rrtm_reice   , rrtm_reliq , rrtm_sw_tauaer   ,&
             rrtm_sw_ssaaer  , rrtm_sw_asmaer  , rrtm_sw_ecaer ,             &
             rrtm_swuflx     , rrtm_swdflx  , rrtm_swhr  ,                   &
             rrtm_swuflxc    , rrtm_swdflxc , rrtm_swhrc )
  
!
!--          Save fluxes
             DO k = nzb, nzt+1
                rad_sw_in(k,:,:)  = rrtm_swdflx(0,k)
                rad_sw_out(k,:,:) = rrtm_swuflx(0,k)
             ENDDO

!
!--          Save heating rates (convert from K/d to K/s)
             DO k = nzb+1, nzt+1
                rad_sw_hr(k,:,:)     = rrtm_swhr(0,k)  * d_hours_day
                rad_sw_cs_hr(k,:,:)  = rrtm_swhrc(0,k) * d_hours_day
             ENDDO

!
!--          Save surface radiative fluxes onto respective surface elements
!--          Horizontal surfaces
             IF ( surf_lsm_h%ns > 0 )  THEN
                   surf_lsm_h%rad_sw_in     = rrtm_swdflx(0,nzb)
                   surf_lsm_h%rad_sw_out    = rrtm_swuflx(0,nzb)
             ENDIF
             IF ( surf_usm_h%ns > 0 )  THEN
                   surf_usm_h%rad_sw_in     = rrtm_swdflx(0,nzb)
                   surf_usm_h%rad_sw_out    = rrtm_swuflx(0,nzb)
             ENDIF
!
!--          Vertical surfaces. Fluxes are obtain at respective vertical 
!--          level of the surface element
             DO  l = 0, 3
                IF ( surf_lsm_v(l)%ns > 0 )  THEN
                      surf_lsm_v(l)%rad_sw_in  = rrtm_swdflx(0,nzb)
                      surf_lsm_v(l)%rad_sw_out = rrtm_swuflx(0,nzb)
                ENDIF             
                IF ( surf_usm_v(l)%ns > 0 )  THEN
                      surf_usm_v(l)%rad_sw_in  = rrtm_swdflx(0,nzb)
                      surf_usm_v(l)%rad_sw_out = rrtm_swuflx(0,nzb)
                ENDIF       
             ENDDO

          ENDIF
!
!--    RRTMG is called for each (j,i) grid point separately, starting at the 
!--    highest topography level
       ELSE
!
!--       Loop over all grid points
          DO i = nxl, nxr
             DO j = nys, nyn

!
!--             Prepare profiles of temperature and H2O volume mixing ratio
                rrtm_tlev(0,nzb+1) = pt(nzb,j,i) * ( surface_pressure          &
                                                     / 1000.0_wp )**0.286_wp


                IF ( cloud_physics )  THEN
                   DO k = nzb+1, nzt+1
                      rrtm_tlay(0,k) = pt(k,j,i) * ( (hyp(k) ) / 100000.0_wp   &
                                       )**0.286_wp + l_d_cp * ql(k,j,i)
                      rrtm_h2ovmr(0,k) = mol_mass_air_d_wv * (q(k,j,i) - ql(k,j,i))
                   ENDDO
                ELSE
                   DO k = nzb+1, nzt+1
                      rrtm_tlay(0,k) = pt(k,j,i) * ( (hyp(k) ) / 100000.0_wp   &
                                       )**0.286_wp
                      rrtm_h2ovmr(0,k) = 0.0_wp
                   ENDDO
                ENDIF

!
!--             Avoid temperature/humidity jumps at the top of the LES domain by 
!--             linear interpolation from nzt+2 to nzt+7
                DO k = nzt+2, nzt+7
                   rrtm_tlay(0,k) = rrtm_tlay(0,nzt+1)                         &
                                 + ( rrtm_tlay(0,nzt+8) - rrtm_tlay(0,nzt+1) ) &
                                 / ( rrtm_play(0,nzt+8) - rrtm_play(0,nzt+1) ) &
                                 * ( rrtm_play(0,k)     - rrtm_play(0,nzt+1) )

                   rrtm_h2ovmr(0,k) = rrtm_h2ovmr(0,nzt+1)                     &
                              + ( rrtm_h2ovmr(0,nzt+8) - rrtm_h2ovmr(0,nzt+1) )&
                              / ( rrtm_play(0,nzt+8)   - rrtm_play(0,nzt+1)   )&
                              * ( rrtm_play(0,k)       - rrtm_play(0,nzt+1) )

                ENDDO

!--             Linear interpolate to zw grid
                DO k = nzb+2, nzt+8
                   rrtm_tlev(0,k)   = rrtm_tlay(0,k-1) + (rrtm_tlay(0,k) -     &
                                      rrtm_tlay(0,k-1))                        &
                                      / ( rrtm_play(0,k) - rrtm_play(0,k-1) )  &
                                      * ( rrtm_plev(0,k) - rrtm_play(0,k-1) )
                ENDDO


!
!--             Calculate liquid water path and cloud fraction for each column.
!--             Note that LWP is required in g/m² instead of kg/kg m.
                rrtm_cldfr  = 0.0_wp
                rrtm_reliq  = 0.0_wp
                rrtm_cliqwp = 0.0_wp
                rrtm_icld   = 0

                IF ( cloud_physics  .OR.  cloud_droplets )  THEN
                   DO k = nzb+1, nzt+1
                      rrtm_cliqwp(0,k) =  ql(k,j,i) * 1000.0_wp *              &
                                          (rrtm_plev(0,k) - rrtm_plev(0,k+1))  &
                                          * 100.0_wp / g 

                      IF ( rrtm_cliqwp(0,k) > 0.0_wp )  THEN
                         rrtm_cldfr(0,k) = 1.0_wp
                         IF ( rrtm_icld == 0 )  rrtm_icld = 1

!
!--                      Calculate cloud droplet effective radius
                         IF ( cloud_physics )  THEN
!
!--                         Calculete effective droplet radius. In case of using 
!--                         cloud_scheme = 'morrison' and a non reasonable number 
!--                         of cloud droplets the inital aerosol number  
!--                         concentration is considered.
                            IF ( microphysics_morrison )  THEN
                               IF ( nc(k,j,i) > 1.0E-20_wp )  THEN
                                  nc_rad = nc(k,j,i)
                               ELSE
                                  nc_rad = na_init
                               ENDIF 
                            ELSE
                               nc_rad = nc_const
                            ENDIF  

                            rrtm_reliq(0,k) = 1.0E6_wp * ( 3.0_wp * ql(k,j,i)     &
                                              * rho_surface                       &
                                              / ( 4.0_wp * pi * nc_rad * rho_l )  &
                                              )**0.33333333333333_wp              &
                                              * EXP( LOG( sigma_gc )**2 )

                         ELSEIF ( cloud_droplets )  THEN
                            number_of_particles = prt_count(k,j,i)

                            IF (number_of_particles <= 0)  CYCLE
                            particles => grid_particles(k,j,i)%particles(1:number_of_particles)
                            s_r2 = 0.0_wp
                            s_r3 = 0.0_wp

                            DO  n = 1, number_of_particles
                               IF ( particles(n)%particle_mask )  THEN
                                  s_r2 = s_r2 + particles(n)%radius**2 *       &
                                         particles(n)%weight_factor
                                  s_r3 = s_r3 + particles(n)%radius**3 *       &
                                         particles(n)%weight_factor
                               ENDIF
                            ENDDO

                            IF ( s_r2 > 0.0_wp )  rrtm_reliq(0,k) = s_r3 / s_r2

                         ENDIF

!
!--                      Limit effective radius
                         IF ( rrtm_reliq(0,k) > 0.0_wp )  THEN
                            rrtm_reliq(0,k) = MAX(rrtm_reliq(0,k),2.5_wp)
                            rrtm_reliq(0,k) = MIN(rrtm_reliq(0,k),60.0_wp)
                        ENDIF
                      ENDIF
                   ENDDO
                ENDIF

!
!--             Write surface emissivity and surface temperature at current 
!--             surface element on RRTMG-shaped array.
!--             Please note, as RRTMG is a single column model, surface attributes
!--             are only obtained from horizontally aligned surfaces (for 
!--             simplicity). Taking surface attributes from horizontal and 
!--             vertical walls would lead to multiple solutions.  
!--             Moreover, for natural- and urban-type surfaces, several surface
!--             classes can exist at a surface element next to each other. 
!--             To obtain bulk parameters, apply a weighted average for these 
!--             surfaces. 
                DO  m = surf_lsm_h%start_index(j,i), surf_lsm_h%end_index(j,i)
                   rrtm_emis = surf_lsm_h%frac(ind_veg_wall,m)  *              &
                               surf_lsm_h%emissivity(ind_veg_wall,m)  +        &
                               surf_lsm_h%frac(ind_pav_green,m) *              &
                               surf_lsm_h%emissivity(ind_pav_green,m) +        & 
                               surf_lsm_h%frac(ind_wat_win,m)   *              &
                               surf_lsm_h%emissivity(ind_wat_win,m)
                   rrtm_tsfc = pt(surf_lsm_h%k(m)+surf_lsm_h%koff,j,i) *       &
                                       (surface_pressure / 1000.0_wp )**0.286_wp
                ENDDO             
                DO  m = surf_usm_h%start_index(j,i), surf_usm_h%end_index(j,i)
                   rrtm_emis = surf_usm_h%frac(ind_veg_wall,m)  *              &
                               surf_usm_h%emissivity(ind_veg_wall,m)  +        &
                               surf_usm_h%frac(ind_pav_green,m) *              &
                               surf_usm_h%emissivity(ind_pav_green,m) +        & 
                               surf_usm_h%frac(ind_wat_win,m)   *              &
                               surf_usm_h%emissivity(ind_wat_win,m)
                   rrtm_tsfc = pt(surf_usm_h%k(m)+surf_usm_h%koff,j,i) *       &
                                       (surface_pressure / 1000.0_wp )**0.286_wp
                ENDDO
!
!--             Obtain topography top index (lower bound of RRTMG)
                k_topo = get_topography_top_index_ji( j, i, 's' )

                IF ( lw_radiation )  THEN
!
!--                Due to technical reasons, copy optical depth to dummy arguments
!--                which are allocated on the exact size as the rrtmg_lw is called.
!--                As one dimesion is allocated with zero size, compiler complains
!--                that rank of the array does not match that of the 
!--                assumed-shaped arguments in the RRTMG library. In order to 
!--                avoid this, write to dummy arguments and give pass the entire 
!--                dummy array. Seems to be the only existing work-around.  
                   ALLOCATE( rrtm_lw_taucld_dum(1:nbndlw+1,0:0,k_topo+1:nzt_rad+1) )
                   ALLOCATE( rrtm_lw_tauaer_dum(0:0,k_topo+1:nzt_rad+1,1:nbndlw+1) )

                   rrtm_lw_taucld_dum =                                        &
                               rrtm_lw_taucld(1:nbndlw+1,0:0,k_topo+1:nzt_rad+1)
                   rrtm_lw_tauaer_dum =                                        &
                               rrtm_lw_tauaer(0:0,k_topo+1:nzt_rad+1,1:nbndlw+1)

                   CALL rrtmg_lw( 1,                                           &                                        
                                  nzt_rad-k_topo,                              &
                                  rrtm_icld,                                   &
                                  rrtm_idrv,                                   &
                                  rrtm_play(:,k_topo+1:nzt_rad+1),             &
                                  rrtm_plev(:,k_topo+1:nzt_rad+2),             &
                                  rrtm_tlay(:,k_topo+1:nzt_rad+1),             &
                                  rrtm_tlev(:,k_topo+1:nzt_rad+2),             &
                                  rrtm_tsfc,                                   &
                                  rrtm_h2ovmr(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_o3vmr(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_co2vmr(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_ch4vmr(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_n2ovmr(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_o2vmr(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_cfc11vmr(:,k_topo+1:nzt_rad+1),         &
                                  rrtm_cfc12vmr(:,k_topo+1:nzt_rad+1),         &
                                  rrtm_cfc22vmr(:,k_topo+1:nzt_rad+1),         &
                                  rrtm_ccl4vmr(:,k_topo+1:nzt_rad+1),          &
                                  rrtm_emis,                                   &
                                  rrtm_inflglw,                                &
                                  rrtm_iceflglw,                               &
                                  rrtm_liqflglw,                               &
                                  rrtm_cldfr(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_lw_taucld_dum,                          &
                                  rrtm_cicewp(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_cliqwp(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_reice(:,k_topo+1:nzt_rad+1),            & 
                                  rrtm_reliq(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_lw_tauaer_dum,                          &
                                  rrtm_lwuflx(:,k_topo:nzt_rad+1),             &
                                  rrtm_lwdflx(:,k_topo:nzt_rad+1),             &
                                  rrtm_lwhr(:,k_topo+1:nzt_rad+1),             &
                                  rrtm_lwuflxc(:,k_topo:nzt_rad+1),            &
                                  rrtm_lwdflxc(:,k_topo:nzt_rad+1),            &
                                  rrtm_lwhrc(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_lwuflx_dt(:,k_topo:nzt_rad+1),          &
                                  rrtm_lwuflxc_dt(:,k_topo:nzt_rad+1) )

                   DEALLOCATE ( rrtm_lw_taucld_dum )
                   DEALLOCATE ( rrtm_lw_tauaer_dum )
!
!--                Save fluxes
                   DO k = k_topo, nzt+1
                      rad_lw_in(k,j,i)  = rrtm_lwdflx(0,k)
                      rad_lw_out(k,j,i) = rrtm_lwuflx(0,k)
                   ENDDO

!
!--                Save heating rates (convert from K/d to K/h)
                   DO k = k_topo+1, nzt+1
                      rad_lw_hr(k,j,i)     = rrtm_lwhr(0,k)  * d_hours_day
                      rad_lw_cs_hr(k,j,i)  = rrtm_lwhrc(0,k) * d_hours_day
                   ENDDO

!
!--                Save surface radiative fluxes and change in LW heating rate 
!--                onto respective surface elements
!--                Horizontal surfaces
                   DO  m = surf_lsm_h%start_index(j,i),                        &
                           surf_lsm_h%end_index(j,i)
                      surf_lsm_h%rad_lw_in(m)           = rrtm_lwdflx(0,k_topo)
                      surf_lsm_h%rad_lw_out(m)          = rrtm_lwuflx(0,k_topo)
                      surf_lsm_h%rad_lw_out_change_0(m) = rrtm_lwuflx_dt(0,k_topo)
                   ENDDO             
                   DO  m = surf_usm_h%start_index(j,i),                        &
                           surf_usm_h%end_index(j,i)
                      surf_usm_h%rad_lw_in(m)           = rrtm_lwdflx(0,k_topo)
                      surf_usm_h%rad_lw_out(m)          = rrtm_lwuflx(0,k_topo)
                      surf_usm_h%rad_lw_out_change_0(m) = rrtm_lwuflx_dt(0,k_topo)
                   ENDDO 
!
!--                Vertical surfaces. Fluxes are obtain at vertical level of the 
!--                respective surface element
                   DO  l = 0, 3
                      DO  m = surf_lsm_v(l)%start_index(j,i),                  &
                              surf_lsm_v(l)%end_index(j,i)
                         k                                    = surf_lsm_v(l)%k(m)
                         surf_lsm_v(l)%rad_lw_in(m)           = rrtm_lwdflx(0,k)
                         surf_lsm_v(l)%rad_lw_out(m)          = rrtm_lwuflx(0,k)
                         surf_lsm_v(l)%rad_lw_out_change_0(m) = rrtm_lwuflx_dt(0,k)
                      ENDDO             
                      DO  m = surf_usm_v(l)%start_index(j,i),                  &
                              surf_usm_v(l)%end_index(j,i)
                         k                                    = surf_usm_v(l)%k(m)
                         surf_usm_v(l)%rad_lw_in(m)           = rrtm_lwdflx(0,k)
                         surf_usm_v(l)%rad_lw_out(m)          = rrtm_lwuflx(0,k)
                         surf_usm_v(l)%rad_lw_out_change_0(m) = rrtm_lwuflx_dt(0,k)
                      ENDDO 
                   ENDDO

                ENDIF

                IF ( sw_radiation .AND. sun_up )  THEN
!
!--                Get albedo for direct/diffusive long/shortwave radiation at 
!--                current (y,x)-location from surface variables. 
!--                Only obtain it from horizontal surfaces, as RRTMG is a single 
!--                column model
!--                (Please note, only one loop will entered, controlled by 
!--                start-end index.)
                   DO  m = surf_lsm_h%start_index(j,i),                        &
                           surf_lsm_h%end_index(j,i)
                      rrtm_asdir(1)  = SUM( surf_lsm_h%frac(:,m) *             &
                                            surf_lsm_h%rrtm_asdir(:,m) )
                      rrtm_asdif(1)  = SUM( surf_lsm_h%frac(:,m) *             &
                                            surf_lsm_h%rrtm_asdif(:,m) )
                      rrtm_aldir(1)  = SUM( surf_lsm_h%frac(:,m) *             &
                                            surf_lsm_h%rrtm_aldir(:,m) )
                      rrtm_aldif(1)  = SUM( surf_lsm_h%frac(:,m) *             &
                                            surf_lsm_h%rrtm_aldif(:,m) )
                   ENDDO             
                   DO  m = surf_usm_h%start_index(j,i),                        &
                           surf_usm_h%end_index(j,i)
                      rrtm_asdir(1)  = SUM( surf_usm_h%frac(:,m) *             &
                                            surf_usm_h%rrtm_asdir(:,m) )
                      rrtm_asdif(1)  = SUM( surf_usm_h%frac(:,m) *             &
                                            surf_usm_h%rrtm_asdif(:,m) )
                      rrtm_aldir(1)  = SUM( surf_usm_h%frac(:,m) *             &
                                            surf_usm_h%rrtm_aldir(:,m) )
                      rrtm_aldif(1)  = SUM( surf_usm_h%frac(:,m) *             &
                                            surf_usm_h%rrtm_aldif(:,m) )
                   ENDDO
!
!--                Due to technical reasons, copy optical depths and other 
!--                to dummy arguments which are allocated on the exact size as the 
!--                rrtmg_sw is called.
!--                As one dimesion is allocated with zero size, compiler complains
!--                that rank of the array does not match that of the 
!--                assumed-shaped arguments in the RRTMG library. In order to 
!--                avoid this, write to dummy arguments and give pass the entire 
!--                dummy array. Seems to be the only existing work-around.  
                   ALLOCATE( rrtm_sw_taucld_dum(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1) )
                   ALLOCATE( rrtm_sw_ssacld_dum(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1) )
                   ALLOCATE( rrtm_sw_asmcld_dum(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1) )
                   ALLOCATE( rrtm_sw_fsfcld_dum(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1) )
                   ALLOCATE( rrtm_sw_tauaer_dum(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1) )
                   ALLOCATE( rrtm_sw_ssaaer_dum(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1) )
                   ALLOCATE( rrtm_sw_asmaer_dum(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1) )
                   ALLOCATE( rrtm_sw_ecaer_dum(0:0,k_topo+1:nzt_rad+1,1:naerec+1)  )
     
                   rrtm_sw_taucld_dum = rrtm_sw_taucld(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1)
                   rrtm_sw_ssacld_dum = rrtm_sw_ssacld(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1)
                   rrtm_sw_asmcld_dum = rrtm_sw_asmcld(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1)
                   rrtm_sw_fsfcld_dum = rrtm_sw_fsfcld(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1)
                   rrtm_sw_tauaer_dum = rrtm_sw_tauaer(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1)
                   rrtm_sw_ssaaer_dum = rrtm_sw_ssaaer(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1)
                   rrtm_sw_asmaer_dum = rrtm_sw_asmaer(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1)
                   rrtm_sw_ecaer_dum  = rrtm_sw_ecaer(0:0,k_topo+1:nzt_rad+1,1:naerec+1)

                   CALL rrtmg_sw( 1,                                           &
                                  nzt_rad-k_topo,                              &
                                  rrtm_icld,                                   &
                                  rrtm_iaer,                                   &
                                  rrtm_play(:,k_topo+1:nzt_rad+1),             &
                                  rrtm_plev(:,k_topo+1:nzt_rad+2),             &
                                  rrtm_tlay(:,k_topo+1:nzt_rad+1),             &
                                  rrtm_tlev(:,k_topo+1:nzt_rad+2),             &
                                  rrtm_tsfc,                                   &
                                  rrtm_h2ovmr(:,k_topo+1:nzt_rad+1),           &                               
                                  rrtm_o3vmr(:,k_topo+1:nzt_rad+1),            &        
                                  rrtm_co2vmr(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_ch4vmr(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_n2ovmr(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_o2vmr(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_asdir,                                  & 
                                  rrtm_asdif,                                  &
                                  rrtm_aldir,                                  &
                                  rrtm_aldif,                                  &
                                  zenith,                                      &
                                  0.0_wp,                                      &
                                  day_of_year,                                 &
                                  solar_constant,                              &
                                  rrtm_inflgsw,                                &
                                  rrtm_iceflgsw,                               &
                                  rrtm_liqflgsw,                               &
                                  rrtm_cldfr(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_sw_taucld_dum,                          &
                                  rrtm_sw_ssacld_dum,                          &
                                  rrtm_sw_asmcld_dum,                          &
                                  rrtm_sw_fsfcld_dum,                          &
                                  rrtm_cicewp(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_cliqwp(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_reice(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_reliq(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_sw_tauaer_dum,                          &
                                  rrtm_sw_ssaaer_dum,                          &
                                  rrtm_sw_asmaer_dum,                          &
                                  rrtm_sw_ecaer_dum,                           &
                                  rrtm_swuflx(:,k_topo:nzt_rad+1),             &  
                                  rrtm_swdflx(:,k_topo:nzt_rad+1),             & 
                                  rrtm_swhr(:,k_topo+1:nzt_rad+1),             & 
                                  rrtm_swuflxc(:,k_topo:nzt_rad+1),            & 
                                  rrtm_swdflxc(:,k_topo:nzt_rad+1),            &
                                  rrtm_swhrc(:,k_topo+1:nzt_rad+1) )

                   DEALLOCATE( rrtm_sw_taucld_dum )
                   DEALLOCATE( rrtm_sw_ssacld_dum )
                   DEALLOCATE( rrtm_sw_asmcld_dum )
                   DEALLOCATE( rrtm_sw_fsfcld_dum )
                   DEALLOCATE( rrtm_sw_tauaer_dum )
                   DEALLOCATE( rrtm_sw_ssaaer_dum )
                   DEALLOCATE( rrtm_sw_asmaer_dum )
                   DEALLOCATE( rrtm_sw_ecaer_dum )
!
!--                Save fluxes
                   DO k = nzb, nzt+1
                      rad_sw_in(k,j,i)  = rrtm_swdflx(0,k)
                      rad_sw_out(k,j,i) = rrtm_swuflx(0,k)
                   ENDDO
!
!--                Save heating rates (convert from K/d to K/s)
                   DO k = nzb+1, nzt+1
                      rad_sw_hr(k,j,i)     = rrtm_swhr(0,k)  * d_hours_day
                      rad_sw_cs_hr(k,j,i)  = rrtm_swhrc(0,k) * d_hours_day
                   ENDDO

!
!--                Save surface radiative fluxes onto respective surface elements
!--                Horizontal surfaces
                   DO  m = surf_lsm_h%start_index(j,i),                        &
                           surf_lsm_h%end_index(j,i)
                      surf_lsm_h%rad_sw_in(m)     = rrtm_swdflx(0,k_topo)
                      surf_lsm_h%rad_sw_out(m)    = rrtm_swuflx(0,k_topo)
                   ENDDO             
                   DO  m = surf_usm_h%start_index(j,i),                        &
                           surf_usm_h%end_index(j,i)
                      surf_usm_h%rad_sw_in(m)     = rrtm_swdflx(0,k_topo)
                      surf_usm_h%rad_sw_out(m)    = rrtm_swuflx(0,k_topo)
                   ENDDO 
!
!--                Vertical surfaces. Fluxes are obtain at respective vertical 
!--                level of the surface element
                   DO  l = 0, 3
                      DO  m = surf_lsm_v(l)%start_index(j,i),                  &
                              surf_lsm_v(l)%end_index(j,i)
                         k                           = surf_lsm_v(l)%k(m)
                         surf_lsm_v(l)%rad_sw_in(m)  = rrtm_swdflx(0,k)
                         surf_lsm_v(l)%rad_sw_out(m) = rrtm_swuflx(0,k)
                      ENDDO             
                      DO  m = surf_usm_v(l)%start_index(j,i),                  &
                              surf_usm_v(l)%end_index(j,i)
                         k                           = surf_usm_v(l)%k(m)
                         surf_usm_v(l)%rad_sw_in(m)  = rrtm_swdflx(0,k)
                         surf_usm_v(l)%rad_sw_out(m) = rrtm_swuflx(0,k)
                      ENDDO 
                   ENDDO

                ENDIF

             ENDDO
          ENDDO

       ENDIF
!
!--    Finally, calculate surface net radiation for surface elements.
!--    First, for horizontal surfaces    
       DO  m = 1, surf_lsm_h%ns
          surf_lsm_h%rad_net(m) = surf_lsm_h%rad_sw_in(m)                      &
                                - surf_lsm_h%rad_sw_out(m)                     &
                                + surf_lsm_h%rad_lw_in(m)                      &
                                - surf_lsm_h%rad_lw_out(m)
       ENDDO
       DO  m = 1, surf_usm_h%ns
          surf_usm_h%rad_net(m) = surf_usm_h%rad_sw_in(m)                      &
                                - surf_usm_h%rad_sw_out(m)                     &
                                + surf_usm_h%rad_lw_in(m)                      &
                                - surf_usm_h%rad_lw_out(m)
       ENDDO
!
!--    Vertical surfaces. 
!--    Todo: weight with azimuth and zenith angle according to their orientation!
       DO  l = 0, 3     
          DO  m = 1, surf_lsm_v(l)%ns
             surf_lsm_v(l)%rad_net(m) = surf_lsm_v(l)%rad_sw_in(m)             &
                                      - surf_lsm_v(l)%rad_sw_out(m)            &
                                      + surf_lsm_v(l)%rad_lw_in(m)             &
                                      - surf_lsm_v(l)%rad_lw_out(m)
          ENDDO
          DO  m = 1, surf_usm_v(l)%ns
             surf_usm_v(l)%rad_net(m) = surf_usm_v(l)%rad_sw_in(m)             &
                                      - surf_usm_v(l)%rad_sw_out(m)            &
                                      + surf_usm_v(l)%rad_lw_in(m)             &
                                      - surf_usm_v(l)%rad_lw_out(m)
          ENDDO
       ENDDO


       CALL exchange_horiz( rad_lw_in,  nbgp )
       CALL exchange_horiz( rad_lw_out, nbgp )
       CALL exchange_horiz( rad_lw_hr,    nbgp )
       CALL exchange_horiz( rad_lw_cs_hr, nbgp )

       CALL exchange_horiz( rad_sw_in,  nbgp )
       CALL exchange_horiz( rad_sw_out, nbgp ) 
       CALL exchange_horiz( rad_sw_hr,    nbgp )
       CALL exchange_horiz( rad_sw_cs_hr, nbgp )

#endif

    END SUBROUTINE radiation_rrtmg


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the cosine of the zenith angle (variable is called zenith)
!------------------------------------------------------------------------------!
    SUBROUTINE calc_zenith

       IMPLICIT NONE

       REAL(wp) ::  declination,  & !< solar declination angle
                    hour_angle      !< solar hour angle
!
!--    Calculate current day and time based on the initial values and simulation
!--    time
       CALL calc_date_and_time

!
!--    Calculate solar declination and hour angle   
       declination = ASIN( decl_1 * SIN(decl_2 * REAL(day_of_year, KIND=wp) - decl_3) )
       hour_angle  = 2.0_wp * pi * (time_utc / 86400.0_wp) + lon - pi

!
!--    Calculate cosine of solar zenith angle
       zenith(0) = SIN(lat) * SIN(declination) + COS(lat) * COS(declination)   &
                                            * COS(hour_angle)
       zenith(0) = MAX(0.0_wp,zenith(0))

!
!--    Calculate solar directional vector
       IF ( sun_direction )  THEN

!
!--       Direction in longitudes equals to sin(solar_azimuth) * sin(zenith)
          sun_dir_lon(0) = -SIN(hour_angle) * COS(declination)

!
!--       Direction in latitues equals to cos(solar_azimuth) * sin(zenith)
          sun_dir_lat(0) = SIN(declination) * COS(lat) - COS(hour_angle) &
                              * COS(declination) * SIN(lat)
       ENDIF

!
!--    Check if the sun is up (otheriwse shortwave calculations can be skipped)
       IF ( zenith(0) > 0.0_wp )  THEN
          sun_up = .TRUE.
       ELSE
          sun_up = .FALSE.
       END IF

    END SUBROUTINE calc_zenith

#if defined ( __rrtmg ) && defined ( __netcdf )
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates surface albedo components based on Briegleb (1992) and 
!> Briegleb et al. (1986)
!------------------------------------------------------------------------------!
    SUBROUTINE calc_albedo( surf )

        IMPLICIT NONE

        INTEGER(iwp)    ::  ind_type !< running index surface tiles
        INTEGER(iwp)    ::  m        !< running index surface elements

        TYPE(surf_type) ::  surf !< treated surfaces

        IF ( sun_up  .AND.  .NOT. average_radiation )  THEN

           DO  m = 1, surf%ns
!
!--           Loop over surface elements
              DO  ind_type = 0, SIZE( surf%albedo_type, 1 ) - 1
           
!
!--              Ocean
                 IF ( surf%albedo_type(ind_type,m) == 1 )  THEN
                    surf%rrtm_aldir(ind_type,m) = 0.026_wp /                    &
                                                ( zenith(0)**1.7_wp + 0.065_wp )&
                                     + 0.15_wp * ( zenith(0) - 0.1_wp )         &
                                               * ( zenith(0) - 0.5_wp )         &
                                               * ( zenith(0) - 1.0_wp )
                    surf%rrtm_asdir(ind_type,m) = surf%rrtm_aldir(ind_type,m)
!
!--              Snow
                 ELSEIF ( surf%albedo_type(ind_type,m) == 16 )  THEN
                    IF ( zenith(0) < 0.5_wp )  THEN
                       surf%rrtm_aldir(ind_type,m) =                           &
                                 0.5_wp * ( 1.0_wp - surf%aldif(ind_type,m) )  &
                                        * ( 3.0_wp / ( 1.0_wp + 4.0_wp         &
                                        * zenith(0) ) ) - 1.0_wp
                       surf%rrtm_asdir(ind_type,m) =                           &
                                 0.5_wp * ( 1.0_wp - surf%asdif(ind_type,m) )  &
                                        * ( 3.0_wp / ( 1.0_wp + 4.0_wp         &
                                        * zenith(0) ) ) - 1.0_wp

                       surf%rrtm_aldir(ind_type,m) =                           &
                                       MIN(0.98_wp, surf%rrtm_aldir(ind_type,m))
                       surf%rrtm_asdir(ind_type,m) =                           &
                                       MIN(0.98_wp, surf%rrtm_asdir(ind_type,m))
                    ELSE
                       surf%rrtm_aldir(ind_type,m) = surf%aldif(ind_type,m)
                       surf%rrtm_asdir(ind_type,m) = surf%asdif(ind_type,m)
                    ENDIF
!
!--              Sea ice
                 ELSEIF ( surf%albedo_type(ind_type,m) == 15 )  THEN
                    surf%rrtm_aldir(ind_type,m) = surf%aldif(ind_type,m)
                    surf%rrtm_asdir(ind_type,m) = surf%asdif(ind_type,m)

!
!--              Asphalt
                 ELSEIF ( surf%albedo_type(ind_type,m) == 17 )  THEN
                    surf%rrtm_aldir(ind_type,m) = surf%aldif(ind_type,m)
                    surf%rrtm_asdir(ind_type,m) = surf%asdif(ind_type,m)


!
!--              Bare soil
                 ELSEIF ( surf%albedo_type(ind_type,m) == 18 )  THEN
                    surf%rrtm_aldir(ind_type,m) = surf%aldif(ind_type,m)
                    surf%rrtm_asdir(ind_type,m) = surf%asdif(ind_type,m)

!
!--              Land surfaces
                 ELSE
                    SELECT CASE ( surf%albedo_type(ind_type,m) )

!
!--                    Surface types with strong zenith dependence
                       CASE ( 1, 2, 3, 4, 11, 12, 13 )
                          surf%rrtm_aldir(ind_type,m) =                        &
                                surf%aldif(ind_type,m) * 1.4_wp /              &
                                           ( 1.0_wp + 0.8_wp * zenith(0) )
                          surf%rrtm_asdir(ind_type,m) =                        &
                                surf%asdif(ind_type,m) * 1.4_wp /              &
                                           ( 1.0_wp + 0.8_wp * zenith(0) )
!
!--                    Surface types with weak zenith dependence
                       CASE ( 5, 6, 7, 8, 9, 10, 14 )
                          surf%rrtm_aldir(ind_type,m) =                        &
                                surf%aldif(ind_type,m) * 1.1_wp /              &
                                           ( 1.0_wp + 0.2_wp * zenith(0) )
                          surf%rrtm_asdir(ind_type,m) =                        &
                                surf%asdif(ind_type,m) * 1.1_wp /              &
                                           ( 1.0_wp + 0.2_wp * zenith(0) )

                       CASE DEFAULT

                    END SELECT
                 ENDIF
!
!--              Diffusive albedo is taken from Table 2
                 surf%rrtm_aldif(ind_type,m) = surf%aldif(ind_type,m)
                 surf%rrtm_asdif(ind_type,m) = surf%asdif(ind_type,m)
              ENDDO
           ENDDO
!
!--     Set albedo in case of average radiation
        ELSEIF ( sun_up  .AND.  average_radiation )  THEN
           surf%rrtm_asdir = albedo_urb
           surf%rrtm_asdif = albedo_urb
           surf%rrtm_aldir = albedo_urb
           surf%rrtm_aldif = albedo_urb  
!
!--     Darkness
        ELSE
           surf%rrtm_aldir = 0.0_wp
           surf%rrtm_asdir = 0.0_wp
           surf%rrtm_aldif = 0.0_wp
           surf%rrtm_asdif = 0.0_wp
        ENDIF

    END SUBROUTINE calc_albedo

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Read sounding data (pressure and temperature) from RADIATION_DATA.
!------------------------------------------------------------------------------!
    SUBROUTINE read_sounding_data

       IMPLICIT NONE

       INTEGER(iwp) :: id,           & !< NetCDF id of input file
                       id_dim_zrad,  & !< pressure level id in the NetCDF file
                       id_var,       & !< NetCDF variable id
                       k,            & !< loop index
                       nz_snd,       & !< number of vertical levels in the sounding data
                       nz_snd_start, & !< start vertical index for sounding data to be used
                       nz_snd_end      !< end vertical index for souding data to be used

       REAL(wp) :: t_surface           !< actual surface temperature

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  hyp_snd_tmp, & !< temporary hydrostatic pressure profile (sounding)
                                               t_snd_tmp      !< temporary temperature profile (sounding)

!
!--    In case of updates, deallocate arrays first (sufficient to check one
!--    array as the others are automatically allocated). This is required
!--    because nzt_rad might change during the update
       IF ( ALLOCATED ( hyp_snd ) )  THEN
          DEALLOCATE( hyp_snd )
          DEALLOCATE( t_snd )
          DEALLOCATE( q_snd  )
          DEALLOCATE ( rrtm_play )
          DEALLOCATE ( rrtm_plev )
          DEALLOCATE ( rrtm_tlay )
          DEALLOCATE ( rrtm_tlev )

          DEALLOCATE ( rrtm_h2ovmr )
          DEALLOCATE ( rrtm_cicewp )
          DEALLOCATE ( rrtm_cldfr )
          DEALLOCATE ( rrtm_cliqwp )
          DEALLOCATE ( rrtm_reice )
          DEALLOCATE ( rrtm_reliq )
          DEALLOCATE ( rrtm_lw_taucld )
          DEALLOCATE ( rrtm_lw_tauaer )

          DEALLOCATE ( rrtm_lwdflx  )
          DEALLOCATE ( rrtm_lwdflxc )
          DEALLOCATE ( rrtm_lwuflx  )
          DEALLOCATE ( rrtm_lwuflxc )
          DEALLOCATE ( rrtm_lwuflx_dt )
          DEALLOCATE ( rrtm_lwuflxc_dt )
          DEALLOCATE ( rrtm_lwhr  )
          DEALLOCATE ( rrtm_lwhrc )

          DEALLOCATE ( rrtm_sw_taucld )
          DEALLOCATE ( rrtm_sw_ssacld )
          DEALLOCATE ( rrtm_sw_asmcld )
          DEALLOCATE ( rrtm_sw_fsfcld )
          DEALLOCATE ( rrtm_sw_tauaer )
          DEALLOCATE ( rrtm_sw_ssaaer )
          DEALLOCATE ( rrtm_sw_asmaer ) 
          DEALLOCATE ( rrtm_sw_ecaer )   
 
          DEALLOCATE ( rrtm_swdflx  )
          DEALLOCATE ( rrtm_swdflxc )
          DEALLOCATE ( rrtm_swuflx  )
          DEALLOCATE ( rrtm_swuflxc )
          DEALLOCATE ( rrtm_swhr  )
          DEALLOCATE ( rrtm_swhrc )

       ENDIF

!
!--    Open file for reading
       nc_stat = NF90_OPEN( rrtm_input_file, NF90_NOWRITE, id )
       CALL netcdf_handle_error_rad( 'read_sounding_data', 549 )

!
!--    Inquire dimension of z axis and save in nz_snd
       nc_stat = NF90_INQ_DIMID( id, "Pressure", id_dim_zrad )
       nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim_zrad, len = nz_snd )
       CALL netcdf_handle_error_rad( 'read_sounding_data', 551 )

!
! !--    Allocate temporary array for storing pressure data
       ALLOCATE( hyp_snd_tmp(1:nz_snd) )
       hyp_snd_tmp = 0.0_wp


!--    Read pressure from file
       nc_stat = NF90_INQ_VARID( id, "Pressure", id_var )
       nc_stat = NF90_GET_VAR( id, id_var, hyp_snd_tmp(:), start = (/1/),      &
                               count = (/nz_snd/) )
       CALL netcdf_handle_error_rad( 'read_sounding_data', 552 )

!
!--    Allocate temporary array for storing temperature data
       ALLOCATE( t_snd_tmp(1:nz_snd) )
       t_snd_tmp = 0.0_wp

!
!--    Read temperature from file
       nc_stat = NF90_INQ_VARID( id, "ReferenceTemperature", id_var )
       nc_stat = NF90_GET_VAR( id, id_var, t_snd_tmp(:), start = (/1/),        &
                               count = (/nz_snd/) )
       CALL netcdf_handle_error_rad( 'read_sounding_data', 553 )

!
!--    Calculate start of sounding data
       nz_snd_start = nz_snd + 1
       nz_snd_end   = nz_snd + 1

!
!--    Start filling vertical dimension at 10hPa above the model domain (hyp is
!--    in Pa, hyp_snd in hPa).
       DO  k = 1, nz_snd
          IF ( hyp_snd_tmp(k) < ( hyp(nzt+1) - 1000.0_wp) * 0.01_wp )  THEN
             nz_snd_start = k
             EXIT
          END IF
       END DO

       IF ( nz_snd_start <= nz_snd )  THEN
          nz_snd_end = nz_snd
       END IF


!
!--    Calculate of total grid points for RRTMG calculations
       nzt_rad = nzt + nz_snd_end - nz_snd_start + 1

!
!--    Save data above LES domain in hyp_snd, t_snd and q_snd
!--    Note: q_snd_tmp is not calculated at the moment (dry residual atmosphere)
       ALLOCATE( hyp_snd(nzb+1:nzt_rad) )
       ALLOCATE( t_snd(nzb+1:nzt_rad)   )
       ALLOCATE( q_snd(nzb+1:nzt_rad)   )
       hyp_snd = 0.0_wp
       t_snd = 0.0_wp
       q_snd = 0.0_wp

       hyp_snd(nzt+2:nzt_rad) = hyp_snd_tmp(nz_snd_start+1:nz_snd_end)
       t_snd(nzt+2:nzt_rad)   = t_snd_tmp(nz_snd_start+1:nz_snd_end)

       nc_stat = NF90_CLOSE( id )

!
!--    Calculate pressure levels on zu and zw grid. Sounding data is added at 
!--    top of the LES domain. This routine does not consider horizontal or 
!--    vertical variability of pressure and temperature
       ALLOCATE ( rrtm_play(0:0,nzb+1:nzt_rad+1)   )
       ALLOCATE ( rrtm_plev(0:0,nzb+1:nzt_rad+2)   )

       t_surface = pt_surface * ( surface_pressure / 1000.0_wp )**0.286_wp
       DO k = nzb+1, nzt+1
          rrtm_play(0,k) = hyp(k) * 0.01_wp
          rrtm_plev(0,k) = surface_pressure * ( (t_surface - g/cp * zw(k-1)) / &
                         t_surface )**(1.0_wp/0.286_wp)
       ENDDO

       DO k = nzt+2, nzt_rad
          rrtm_play(0,k) = hyp_snd(k)
          rrtm_plev(0,k) = 0.5_wp * ( rrtm_play(0,k) + rrtm_play(0,k-1) )
       ENDDO
       rrtm_plev(0,nzt_rad+1) = MAX( 0.5 * hyp_snd(nzt_rad),                   &
                                   1.5 * hyp_snd(nzt_rad)                      &
                                 - 0.5 * hyp_snd(nzt_rad-1) )
       rrtm_plev(0,nzt_rad+2)  = MIN( 1.0E-4_wp,                               &
                                      0.25_wp * rrtm_plev(0,nzt_rad+1) )

       rrtm_play(0,nzt_rad+1) = 0.5 * rrtm_plev(0,nzt_rad+1)

!
!--    Calculate temperature/humidity levels at top of the LES domain. 
!--    Currently, the temperature is taken from sounding data (might lead to a 
!--    temperature jump at interface. To do: Humidity is currently not 
!--    calculated above the LES domain.
       ALLOCATE ( rrtm_tlay(0:0,nzb+1:nzt_rad+1)   )
       ALLOCATE ( rrtm_tlev(0:0,nzb+1:nzt_rad+2)   )
       ALLOCATE ( rrtm_h2ovmr(0:0,nzb+1:nzt_rad+1) )

       DO k = nzt+8, nzt_rad
          rrtm_tlay(0,k)   = t_snd(k)
          rrtm_h2ovmr(0,k) = q_snd(k)
       ENDDO
       rrtm_tlay(0,nzt_rad+1) = 2.0_wp * rrtm_tlay(0,nzt_rad)                  &
                                - rrtm_tlay(0,nzt_rad-1)
       DO k = nzt+9, nzt_rad+1
          rrtm_tlev(0,k)   = rrtm_tlay(0,k-1) + (rrtm_tlay(0,k)                &
                             - rrtm_tlay(0,k-1))                               &
                             / ( rrtm_play(0,k) - rrtm_play(0,k-1) )           &
                             * ( rrtm_plev(0,k) - rrtm_play(0,k-1) )
       ENDDO
       rrtm_h2ovmr(0,nzt_rad+1) = rrtm_h2ovmr(0,nzt_rad)

       rrtm_tlev(0,nzt_rad+2)   = 2.0_wp * rrtm_tlay(0,nzt_rad+1)              &
                                  - rrtm_tlev(0,nzt_rad)
!
!--    Allocate remaining RRTMG arrays
       ALLOCATE ( rrtm_cicewp(0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_cldfr(0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_cliqwp(0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_reice(0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_reliq(0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_lw_taucld(1:nbndlw+1,0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_lw_tauaer(0:0,nzb+1:nzt_rad+1,1:nbndlw+1) )
       ALLOCATE ( rrtm_sw_taucld(1:nbndsw+1,0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_sw_ssacld(1:nbndsw+1,0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_sw_asmcld(1:nbndsw+1,0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_sw_fsfcld(1:nbndsw+1,0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_sw_tauaer(0:0,nzb+1:nzt_rad+1,1:nbndsw+1) )
       ALLOCATE ( rrtm_sw_ssaaer(0:0,nzb+1:nzt_rad+1,1:nbndsw+1) )
       ALLOCATE ( rrtm_sw_asmaer(0:0,nzb+1:nzt_rad+1,1:nbndsw+1) ) 
       ALLOCATE ( rrtm_sw_ecaer(0:0,nzb+1:nzt_rad+1,1:naerec+1) )    

!
!--    The ice phase is currently not considered in PALM
       rrtm_cicewp = 0.0_wp
       rrtm_reice  = 0.0_wp

!
!--    Set other parameters (move to NAMELIST parameters in the future)
       rrtm_lw_tauaer = 0.0_wp
       rrtm_lw_taucld = 0.0_wp
       rrtm_sw_taucld = 0.0_wp
       rrtm_sw_ssacld = 0.0_wp
       rrtm_sw_asmcld = 0.0_wp
       rrtm_sw_fsfcld = 0.0_wp
       rrtm_sw_tauaer = 0.0_wp
       rrtm_sw_ssaaer = 0.0_wp
       rrtm_sw_asmaer = 0.0_wp
       rrtm_sw_ecaer  = 0.0_wp


       ALLOCATE ( rrtm_swdflx(0:0,nzb:nzt_rad+1)  )
       ALLOCATE ( rrtm_swuflx(0:0,nzb:nzt_rad+1)  )
       ALLOCATE ( rrtm_swhr(0:0,nzb+1:nzt_rad+1)  )
       ALLOCATE ( rrtm_swuflxc(0:0,nzb:nzt_rad+1) )
       ALLOCATE ( rrtm_swdflxc(0:0,nzb:nzt_rad+1) )
       ALLOCATE ( rrtm_swhrc(0:0,nzb+1:nzt_rad+1) )

       rrtm_swdflx  = 0.0_wp
       rrtm_swuflx  = 0.0_wp
       rrtm_swhr    = 0.0_wp  
       rrtm_swuflxc = 0.0_wp
       rrtm_swdflxc = 0.0_wp
       rrtm_swhrc   = 0.0_wp

       ALLOCATE ( rrtm_lwdflx(0:0,nzb:nzt_rad+1)  )
       ALLOCATE ( rrtm_lwuflx(0:0,nzb:nzt_rad+1)  )
       ALLOCATE ( rrtm_lwhr(0:0,nzb+1:nzt_rad+1)  )
       ALLOCATE ( rrtm_lwuflxc(0:0,nzb:nzt_rad+1) )
       ALLOCATE ( rrtm_lwdflxc(0:0,nzb:nzt_rad+1) )
       ALLOCATE ( rrtm_lwhrc(0:0,nzb+1:nzt_rad+1) )

       rrtm_lwdflx  = 0.0_wp
       rrtm_lwuflx  = 0.0_wp
       rrtm_lwhr    = 0.0_wp  
       rrtm_lwuflxc = 0.0_wp
       rrtm_lwdflxc = 0.0_wp
       rrtm_lwhrc   = 0.0_wp

       ALLOCATE ( rrtm_lwuflx_dt(0:0,nzb:nzt_rad+1) )
       ALLOCATE ( rrtm_lwuflxc_dt(0:0,nzb:nzt_rad+1) )

       rrtm_lwuflx_dt = 0.0_wp
       rrtm_lwuflxc_dt = 0.0_wp

    END SUBROUTINE read_sounding_data


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Read trace gas data from file
!------------------------------------------------------------------------------!
    SUBROUTINE read_trace_gas_data

       USE rrsw_ncpar

       IMPLICIT NONE

       INTEGER(iwp), PARAMETER :: num_trace_gases = 9 !< number of trace gases (absorbers)

       CHARACTER(LEN=5), DIMENSION(num_trace_gases), PARAMETER ::              & !< trace gas names
           trace_names = (/'O3   ', 'CO2  ', 'CH4  ', 'N2O  ', 'O2   ',        &
                           'CFC11', 'CFC12', 'CFC22', 'CCL4 '/)

       INTEGER(iwp) :: id,     & !< NetCDF id
                       k,      & !< loop index
                       m,      & !< loop index
                       n,      & !< loop index
                       nabs,   & !< number of absorbers
                       np,     & !< number of pressure levels
                       id_abs, & !< NetCDF id of the respective absorber
                       id_dim, & !< NetCDF id of asborber's dimension
                       id_var    !< NetCDf id ot the absorber

       REAL(wp) :: p_mls_l, p_mls_u, p_wgt_l, p_wgt_u, p_mls_m


       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  p_mls,          & !< pressure levels for the absorbers
                                                 rrtm_play_tmp,  & !< temporary array for pressure zu-levels
                                                 rrtm_plev_tmp,  & !< temporary array for pressure zw-levels
                                                 trace_path_tmp    !< temporary array for storing trace gas path data

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  trace_mls,      & !< array for storing the absorber amounts
                                                 trace_mls_path, & !< array for storing trace gas path data
                                                 trace_mls_tmp     !< temporary array for storing trace gas data


!
!--    In case of updates, deallocate arrays first (sufficient to check one
!--    array as the others are automatically allocated)
       IF ( ALLOCATED ( rrtm_o3vmr ) )  THEN
          DEALLOCATE ( rrtm_o3vmr  )
          DEALLOCATE ( rrtm_co2vmr )
          DEALLOCATE ( rrtm_ch4vmr )
          DEALLOCATE ( rrtm_n2ovmr )
          DEALLOCATE ( rrtm_o2vmr  )
          DEALLOCATE ( rrtm_cfc11vmr )
          DEALLOCATE ( rrtm_cfc12vmr )
          DEALLOCATE ( rrtm_cfc22vmr )
          DEALLOCATE ( rrtm_ccl4vmr  )
       ENDIF

!
!--    Allocate trace gas profiles
       ALLOCATE ( rrtm_o3vmr(0:0,1:nzt_rad+1)  )
       ALLOCATE ( rrtm_co2vmr(0:0,1:nzt_rad+1) )
       ALLOCATE ( rrtm_ch4vmr(0:0,1:nzt_rad+1) )
       ALLOCATE ( rrtm_n2ovmr(0:0,1:nzt_rad+1) )
       ALLOCATE ( rrtm_o2vmr(0:0,1:nzt_rad+1)  )
       ALLOCATE ( rrtm_cfc11vmr(0:0,1:nzt_rad+1) )
       ALLOCATE ( rrtm_cfc12vmr(0:0,1:nzt_rad+1) )
       ALLOCATE ( rrtm_cfc22vmr(0:0,1:nzt_rad+1) )
       ALLOCATE ( rrtm_ccl4vmr(0:0,1:nzt_rad+1)  )

!
!--    Open file for reading
       nc_stat = NF90_OPEN( rrtm_input_file, NF90_NOWRITE, id )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 549 )
!
!--    Inquire dimension ids and dimensions
       nc_stat = NF90_INQ_DIMID( id, "Pressure", id_dim )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )
       nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim, len = np) 
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )

       nc_stat = NF90_INQ_DIMID( id, "Absorber", id_dim )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )
       nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim, len = nabs ) 
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )
   

!
!--    Allocate pressure, and trace gas arrays     
       ALLOCATE( p_mls(1:np) )
       ALLOCATE( trace_mls(1:num_trace_gases,1:np) ) 
       ALLOCATE( trace_mls_tmp(1:nabs,1:np) ) 


       nc_stat = NF90_INQ_VARID( id, "Pressure", id_var )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )
       nc_stat = NF90_GET_VAR( id, id_var, p_mls )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )

       nc_stat = NF90_INQ_VARID( id, "AbsorberAmountMLS", id_var )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )
       nc_stat = NF90_GET_VAR( id, id_var, trace_mls_tmp )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )


!
!--    Write absorber amounts (mls) to trace_mls
       DO n = 1, num_trace_gases
          CALL getAbsorberIndex( TRIM( trace_names(n) ), id_abs )

          trace_mls(n,1:np) = trace_mls_tmp(id_abs,1:np)

!
!--       Replace missing values by zero
          WHERE ( trace_mls(n,:) > 2.0_wp )  
             trace_mls(n,:) = 0.0_wp
          END WHERE
       END DO

       DEALLOCATE ( trace_mls_tmp )

       nc_stat = NF90_CLOSE( id )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 551 )

!
!--    Add extra pressure level for calculations of the trace gas paths
       ALLOCATE ( rrtm_play_tmp(1:nzt_rad+1) )
       ALLOCATE ( rrtm_plev_tmp(1:nzt_rad+2) )

       rrtm_play_tmp(1:nzt_rad)   = rrtm_play(0,1:nzt_rad) 
       rrtm_plev_tmp(1:nzt_rad+1) = rrtm_plev(0,1:nzt_rad+1)
       rrtm_play_tmp(nzt_rad+1)   = rrtm_plev(0,nzt_rad+1) * 0.5_wp
       rrtm_plev_tmp(nzt_rad+2)   = MIN( 1.0E-4_wp, 0.25_wp                    &
                                         * rrtm_plev(0,nzt_rad+1) )
 
!
!--    Calculate trace gas path (zero at surface) with interpolation to the
!--    sounding levels
       ALLOCATE ( trace_mls_path(1:nzt_rad+2,1:num_trace_gases) )

       trace_mls_path(nzb+1,:) = 0.0_wp
       
       DO k = nzb+2, nzt_rad+2
          DO m = 1, num_trace_gases
             trace_mls_path(k,m) = trace_mls_path(k-1,m)

!
!--          When the pressure level is higher than the trace gas pressure
!--          level, assume that 
             IF ( rrtm_plev_tmp(k-1) > p_mls(1) )  THEN             
                
                trace_mls_path(k,m) = trace_mls_path(k,m) + trace_mls(m,1)     &
                                      * ( rrtm_plev_tmp(k-1)                   &
                                          - MAX( p_mls(1), rrtm_plev_tmp(k) )  &
                                        ) / g
             ENDIF

!
!--          Integrate for each sounding level from the contributing p_mls 
!--          levels
             DO n = 2, np
!
!--             Limit p_mls so that it is within the model level
                p_mls_u = MIN( rrtm_plev_tmp(k-1),                             &
                          MAX( rrtm_plev_tmp(k), p_mls(n) ) )
                p_mls_l = MIN( rrtm_plev_tmp(k-1),                             &
                          MAX( rrtm_plev_tmp(k), p_mls(n-1) ) )

                IF ( p_mls_l > p_mls_u )  THEN

!
!--                Calculate weights for interpolation
                   p_mls_m = 0.5_wp * (p_mls_l + p_mls_u)
                   p_wgt_u = (p_mls(n-1) - p_mls_m) / (p_mls(n-1) - p_mls(n))
                   p_wgt_l = (p_mls_m - p_mls(n))   / (p_mls(n-1) - p_mls(n))

!
!--                Add level to trace gas path
                   trace_mls_path(k,m) = trace_mls_path(k,m)                   &
                                         +  ( p_wgt_u * trace_mls(m,n)         &
                                            + p_wgt_l * trace_mls(m,n-1) )     &
                                         * (p_mls_l - p_mls_u) / g
                ENDIF
             ENDDO

             IF ( rrtm_plev_tmp(k) < p_mls(np) )  THEN 
                trace_mls_path(k,m) = trace_mls_path(k,m) + trace_mls(m,np)    &
                                      * ( MIN( rrtm_plev_tmp(k-1), p_mls(np) ) &
                                          - rrtm_plev_tmp(k)                   &
                                        ) / g 
             ENDIF  
          ENDDO
       ENDDO


!
!--    Prepare trace gas path profiles
       ALLOCATE ( trace_path_tmp(1:nzt_rad+1) )

       DO m = 1, num_trace_gases

          trace_path_tmp(1:nzt_rad+1) = ( trace_mls_path(2:nzt_rad+2,m)        &
                                       - trace_mls_path(1:nzt_rad+1,m) ) * g   &
                                       / ( rrtm_plev_tmp(1:nzt_rad+1)          &
                                       - rrtm_plev_tmp(2:nzt_rad+2) )

!
!--       Save trace gas paths to the respective arrays
          SELECT CASE ( TRIM( trace_names(m) ) )

             CASE ( 'O3' )

                rrtm_o3vmr(0,:) = trace_path_tmp(:)

             CASE ( 'CO2' )

                rrtm_co2vmr(0,:) = trace_path_tmp(:)

             CASE ( 'CH4' )

                rrtm_ch4vmr(0,:) = trace_path_tmp(:)

             CASE ( 'N2O' )

                rrtm_n2ovmr(0,:) = trace_path_tmp(:)

             CASE ( 'O2' )

                rrtm_o2vmr(0,:) = trace_path_tmp(:)

             CASE ( 'CFC11' )

                rrtm_cfc11vmr(0,:) = trace_path_tmp(:)

             CASE ( 'CFC12' )

                rrtm_cfc12vmr(0,:) = trace_path_tmp(:)

             CASE ( 'CFC22' )

                rrtm_cfc22vmr(0,:) = trace_path_tmp(:)

             CASE ( 'CCL4' )

                rrtm_ccl4vmr(0,:) = trace_path_tmp(:)

             CASE DEFAULT

          END SELECT

       ENDDO

       DEALLOCATE ( trace_path_tmp )
       DEALLOCATE ( trace_mls_path )
       DEALLOCATE ( rrtm_play_tmp )
       DEALLOCATE ( rrtm_plev_tmp )
       DEALLOCATE ( trace_mls )
       DEALLOCATE ( p_mls )

    END SUBROUTINE read_trace_gas_data


    SUBROUTINE netcdf_handle_error_rad( routine_name, errno )

       USE control_parameters,                                                 &
           ONLY:  message_string

       USE NETCDF

       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=6) ::  message_identifier
       CHARACTER(LEN=*) ::  routine_name

       INTEGER(iwp) ::  errno

       IF ( nc_stat /= NF90_NOERR )  THEN

          WRITE( message_identifier, '(''NC'',I4.4)' )  errno
          message_string = TRIM( NF90_STRERROR( nc_stat ) )

          CALL message( routine_name, message_identifier, 2, 2, 0, 6, 1 )

       ENDIF

    END SUBROUTINE netcdf_handle_error_rad
#endif


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate temperature tendency due to radiative cooling/heating.
!> Cache-optimized version.
!------------------------------------------------------------------------------!
 SUBROUTINE radiation_tendency_ij ( i, j, tend )

    USE cloud_parameters,                                                      &
        ONLY:  pt_d_t

    IMPLICIT NONE

    INTEGER(iwp) :: i, j, k !< loop indices

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) :: tend !< pt tendency term

    IF ( radiation_scheme == 'rrtmg' )  THEN
#if defined  ( __rrtmg )
!
!--    Calculate tendency based on heating rate
       DO k = nzb+1, nzt+1
          tend(k,j,i) = tend(k,j,i) + (rad_lw_hr(k,j,i) + rad_sw_hr(k,j,i))    &
                                         * pt_d_t(k) * d_seconds_hour
       ENDDO
#endif
    ENDIF

    END SUBROUTINE radiation_tendency_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate temperature tendency due to radiative cooling/heating.
!> Vector-optimized version
!------------------------------------------------------------------------------!
 SUBROUTINE radiation_tendency ( tend )

    USE cloud_parameters,                                                      &
        ONLY:  pt_d_t

    USE indices,                                                               &
        ONLY:  nxl, nxr, nyn, nys

    IMPLICIT NONE

    INTEGER(iwp) :: i, j, k !< loop indices

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) :: tend !< pt tendency term

    IF ( radiation_scheme == 'rrtmg' )  THEN
#if defined  ( __rrtmg )
!
!--    Calculate tendency based on heating rate
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO k = nzb+1, nzt+1
                tend(k,j,i) = tend(k,j,i) + ( rad_lw_hr(k,j,i)                 &
                                          +  rad_sw_hr(k,j,i) ) * pt_d_t(k)    &
                                          * d_seconds_hour
             ENDDO
          ENDDO
       ENDDO
#endif
    ENDIF


 END SUBROUTINE radiation_tendency

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine calculates interaction of the solar radiation
!> with urban and land surfaces and updates all surface heatfluxes.
!> It calculates also the required parameters for RRTMG lower BC.
!>
!> For more info. see Resler et al. 2017
!>
!> The new version 2.0 was radically rewriten, the discretization scheme
!> has been changed. This new version significantly improves effectivity
!> of the paralelization and the scalability of the model.
!------------------------------------------------------------------------------!

 SUBROUTINE radiation_interaction

     IMPLICIT NONE

     INTEGER(iwp)                      :: i, j, k, kk, is, js, d, ku, refstep, m, mm, l, ll
     INTEGER(iwp)                      :: nzubl, nzutl, isurf, isurfsrc, isvf, icsf, ipcgb
     INTEGER(iwp)                      :: isd                !< solar direction number
     REAL(wp), DIMENSION(3,3)          :: mrot               !< grid rotation matrix (zyx)
     REAL(wp), DIMENSION(3,0:nsurf_type):: vnorm             !< face direction normal vectors (zyx)
     REAL(wp), DIMENSION(3)            :: sunorig            !< grid rotated solar direction unit vector (zyx)
     REAL(wp), DIMENSION(3)            :: sunorig_grid       !< grid squashed solar direction unit vector (zyx)
     REAL(wp), DIMENSION(0:nsurf_type) :: costheta           !< direct irradiance factor of solar angle
     REAL(wp), DIMENSION(nzub:nzut)    :: pchf_prep          !< precalculated factor for canopy temperature tendency
     REAL(wp), DIMENSION(nzub:nzut)    :: pctf_prep          !< precalculated factor for canopy transpiration tendency
     REAL(wp), PARAMETER               :: alpha = 0._wp      !< grid rotation (TODO: add to namelist or remove)
     REAL(wp)                          :: pc_box_area, pc_abs_frac, pc_abs_eff
     INTEGER(iwp)                      :: pc_box_dimshift    !< transform for best accuracy
     INTEGER(iwp), DIMENSION(0:3)      :: reorder = (/ 1, 0, 3, 2 /)
     REAL(wp), DIMENSION(0:nsurf_type) :: facearea
     REAL(wp)                          :: pabsswl  = 0.0_wp  !< total absorbed SW radiation energy in local processor (W)
     REAL(wp)                          :: pabssw   = 0.0_wp  !< total absorbed SW radiation energy in all processors (W)
     REAL(wp)                          :: pabslwl  = 0.0_wp  !< total absorbed LW radiation energy in local processor (W)
     REAL(wp)                          :: pabslw   = 0.0_wp  !< total absorbed LW radiation energy in all processors (W)
     REAL(wp)                          :: pemitlwl = 0.0_wp  !< total emitted LW radiation energy in all processors (W)
     REAL(wp)                          :: pemitlw  = 0.0_wp  !< total emitted LW radiation energy in all processors (W)
     REAL(wp)                          :: pinswl   = 0.0_wp  !< total received SW radiation energy in local processor (W)
     REAL(wp)                          :: pinsw    = 0.0_wp  !< total received SW radiation energy in all processor (W)
     REAL(wp)                          :: pinlwl   = 0.0_wp  !< total received LW radiation energy in local processor (W)
     REAL(wp)                          :: pinlw    = 0.0_wp  !< total received LW radiation energy in all processor (W)
     REAL(wp)                          :: emiss_sum_surfl    !< sum of emissisivity of surfaces in local processor
     REAL(wp)                          :: emiss_sum_surf     !< sum of emissisivity of surfaces in all processor
     REAL(wp)                          :: area_surfl         !< total area of surfaces in local processor
     REAL(wp)                          :: area_surf          !< total area of surfaces in all processor



#if ! defined( __nopointer )
     IF ( plant_canopy )  THEN
         pchf_prep(:) = r_d * (hyp(nzub:nzut) / 100000.0_wp)**0.286_wp &
                     / (cp * hyp(nzub:nzut) * dx*dy*dz(1)) !< equals to 1 / (rho * c_p * Vbox * T)
         pctf_prep(:) = r_d * (hyp(nzub:nzut) / 100000.0_wp)**0.286_wp &
                     / (l_v * hyp(nzub:nzut) * dx*dy*dz)
     ENDIF
#endif
     sun_direction = .TRUE.
     CALL calc_zenith  !< required also for diffusion radiation

!--     prepare rotated normal vectors and irradiance factor
     vnorm(1,:) = kdir(:)
     vnorm(2,:) = jdir(:)
     vnorm(3,:) = idir(:)
     mrot(1, :) = (/ 1._wp,  0._wp,      0._wp      /)
     mrot(2, :) = (/ 0._wp,  COS(alpha), SIN(alpha) /)
     mrot(3, :) = (/ 0._wp, -SIN(alpha), COS(alpha) /)
     sunorig = (/ zenith(0), sun_dir_lat, sun_dir_lon /)
     sunorig = MATMUL(mrot, sunorig)
     DO d = 0, nsurf_type
         costheta(d) = DOT_PRODUCT(sunorig, vnorm(:,d))
     ENDDO

     IF ( zenith(0) > 0 )  THEN
!--         now we will "squash" the sunorig vector by grid box size in
!--         each dimension, so that this new direction vector will allow us
!--         to traverse the ray path within grid coordinates directly
         sunorig_grid = (/ sunorig(1)/dz(1), sunorig(2)/dy, sunorig(3)/dx /)
!--         sunorig_grid = sunorig_grid / norm2(sunorig_grid)
         sunorig_grid = sunorig_grid / SQRT(SUM(sunorig_grid**2))

         IF ( npcbl > 0 )  THEN
!--            precompute effective box depth with prototype Leaf Area Density
            pc_box_dimshift = MAXLOC(ABS(sunorig), 1) - 1
            CALL box_absorb(CSHIFT((/dz(1),dy,dx/), pc_box_dimshift),      &
                                60, prototype_lad,                          &
                                CSHIFT(ABS(sunorig), pc_box_dimshift),      &
                                pc_box_area, pc_abs_frac)
            pc_box_area = pc_box_area * ABS(sunorig(pc_box_dimshift+1) / sunorig(1))
            pc_abs_eff = LOG(1._wp - pc_abs_frac) / prototype_lad
         ENDIF
     ENDIF

!--     split diffusion and direct part of the solar downward radiation
!--     comming from radiation model and store it in 2D arrays
!--     rad_sw_in_diff, rad_sw_in_dir and rad_lw_in_diff
     IF ( split_diffusion_radiation )  THEN
         CALL calc_diffusion_radiation
     ELSE
         rad_sw_in_diff = 0.0_wp
         rad_sw_in_dir(:,:)  = rad_sw_in(0,:,:)
         rad_lw_in_diff(:,:) = rad_lw_in(0,:,:)
     ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--     First pass: direct + diffuse irradiance + thermal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     surfinswdir   = 0._wp !nsurfl
     surfins       = 0._wp !nsurfl
     surfinl       = 0._wp !nsurfl
     surfoutsl(:)  = 0.0_wp !start-end
     surfoutll(:)  = 0.0_wp !start-end

!--  Set up thermal radiation from surfaces
!--  emiss_surf is defined only for surfaces for which energy balance is calculated
!--  Workaround: reorder surface data type back on 1D array including all surfaces,
!--  which implies to reorder horizontal and vertical surfaces
!
!--  Horizontal walls
     mm = 1
     DO  i = nxl, nxr
        DO  j = nys, nyn
!--           urban
           DO  m = surf_usm_h%start_index(j,i), surf_usm_h%end_index(j,i)
              surfoutll(mm) = SUM ( surf_usm_h%frac(:,m) *                  &
                                    surf_usm_h%emissivity(:,m) )            &
                                  * sigma_sb                                &
                                  * surf_usm_h%pt_surface(m)**4
              albedo_surf(mm) = SUM ( surf_usm_h%frac(:,m) *                &
                                      surf_usm_h%albedo(:,m) )
              emiss_surf(mm)  = SUM ( surf_usm_h%frac(:,m) *                &
                                      surf_usm_h%emissivity(:,m) )
              mm = mm + 1
           ENDDO
!--           land
           DO  m = surf_lsm_h%start_index(j,i), surf_lsm_h%end_index(j,i)
              surfoutll(mm) = SUM ( surf_lsm_h%frac(:,m) *                  &
                                    surf_lsm_h%emissivity(:,m) )            &
                                  * sigma_sb                                &
                                  * surf_lsm_h%pt_surface(m)**4
              albedo_surf(mm) = SUM ( surf_lsm_h%frac(:,m) *                &
                                      surf_lsm_h%albedo(:,m) )
              emiss_surf(mm)  = SUM ( surf_lsm_h%frac(:,m) *                &
                                      surf_lsm_h%emissivity(:,m) )
              mm = mm + 1
           ENDDO
        ENDDO
     ENDDO
!
!--     Vertical walls
     DO  i = nxl, nxr
        DO  j = nys, nyn
           DO  ll = 0, 3
              l = reorder(ll)
!--              urban
              DO  m = surf_usm_v(l)%start_index(j,i),                       &
                      surf_usm_v(l)%end_index(j,i)
                 surfoutll(mm) = SUM ( surf_usm_v(l)%frac(:,m) *            &
                                       surf_usm_v(l)%emissivity(:,m) )      &
                                  * sigma_sb                                &
                                  * surf_usm_v(l)%pt_surface(m)**4
                 albedo_surf(mm) = SUM ( surf_usm_v(l)%frac(:,m) *          &
                                         surf_usm_v(l)%albedo(:,m) )
                 emiss_surf(mm)  = SUM ( surf_usm_v(l)%frac(:,m) *          &
                                         surf_usm_v(l)%emissivity(:,m) )
                 mm = mm + 1
              ENDDO
!--              land
              DO  m = surf_lsm_v(l)%start_index(j,i),                       &
                      surf_lsm_v(l)%end_index(j,i)
                 surfoutll(mm) = SUM ( surf_lsm_v(l)%frac(:,m) *            &
                                       surf_lsm_v(l)%emissivity(:,m) )      &
                                  * sigma_sb                                &
                                  * surf_lsm_v(l)%pt_surface(m)**4
                 albedo_surf(mm) = SUM ( surf_lsm_v(l)%frac(:,m) *          &
                                         surf_lsm_v(l)%albedo(:,m) )
                 emiss_surf(mm)  = SUM ( surf_lsm_v(l)%frac(:,m) *          &
                                         surf_lsm_v(l)%emissivity(:,m) )
                 mm = mm + 1
              ENDDO
           ENDDO
        ENDDO
     ENDDO

#if defined( __parallel )
!--     might be optimized and gather only values relevant for current processor
     CALL MPI_AllGatherv(surfoutll, nsurfl, MPI_REAL, &
                         surfoutl, nsurfs, surfstart, MPI_REAL, comm2d, ierr) !nsurf global
#else
     surfoutl(:) = surfoutll(:) !nsurf global
#endif

     IF ( surface_reflections)  THEN
        DO  isvf = 1, nsvfl
           isurf = svfsurf(1, isvf)
           k     = surfl(iz, isurf)
           j     = surfl(iy, isurf)
           i     = surfl(ix, isurf)
           isurfsrc = svfsurf(2, isvf)
!
!--        For surface-to-surface factors we calculate thermal radiation in 1st pass
           surfinl(isurf) = surfinl(isurf) + svf(1,isvf) * surfoutl(isurfsrc)
        ENDDO
     ENDIF

     !-- diffuse radiation using sky view factor, TODO: homogeneous rad_*w_in_diff because now it depends on no. of processors
     surfinswdif(:) = rad_sw_in_diff(nyn,nxl) * skyvft(:)
     surfinlwdif(:) = rad_lw_in_diff(nyn,nxl) * skyvf(:)

     !-- direct radiation
     IF ( zenith(0) > 0 )  THEN
        !--Identify solar direction vector (discretized number) 1)
        !--
        j = FLOOR(ACOS(zenith(0)) / pi * raytrace_discrete_elevs)
        i = MODULO(NINT(ATAN2(sun_dir_lon(0), sun_dir_lat(0))               &
                        / (2._wp*pi) * raytrace_discrete_azims-.5_wp, iwp), &
                   raytrace_discrete_azims)
        isd = dsidir_rev(j, i)
        DO isurf = 1, nsurfl
           surfinswdir(isurf) = rad_sw_in_dir(nyn,nxl) * costheta(surfl(id, isurf)) * dsitrans(isurf, isd) / zenith(0)
        ENDDO
     ENDIF

     IF ( npcbl > 0 )  THEN

         pcbinswdir(:) = 0._wp
         pcbinswdif(:) = 0._wp
         pcbinlw(:) = 0._wp  !< will stay always 0 since we don't absorb lw anymore
!
!--         pcsf first pass
         DO icsf = 1, ncsfl
             ipcgb = csfsurf(1, icsf)
             i = pcbl(ix,ipcgb)
             j = pcbl(iy,ipcgb)
             k = pcbl(iz,ipcgb)
             isurfsrc = csfsurf(2, icsf)

             IF ( isurfsrc == -1 )  THEN
!--                 Diffuse rad from sky.
                 pcbinswdif(ipcgb) = csf(1,icsf) * csf(2,icsf) * rad_sw_in_diff(j,i)

                 !--Direct rad
                 IF ( zenith(0) > 0 )  THEN
                    !--Estimate directed box absorption
                    pc_abs_frac = 1._wp - exp(pc_abs_eff * lad_s(k,j,i))

                    !--isd has already been established, see 1)
                    pcbinswdir(ipcgb) = rad_sw_in_dir(j, i) * pc_box_area &
                                        * pc_abs_frac * dsitransc(ipcgb, isd)
                 ENDIF

                 EXIT ! only isurfsrc=-1 is processed here
             ENDIF
         ENDDO

         pcbinsw(:) = pcbinswdir(:) + pcbinswdif(:)
     ENDIF
     surfins = surfinswdir + surfinswdif
     surfinl = surfinl + surfinlwdif
     surfinsw = surfins
     surfinlw = surfinl
     surfoutsw = 0.0_wp
     surfoutlw = surfoutll
!        surfhf = surfinsw + surfinlw - surfoutsw - surfoutlw

     IF ( .NOT.  surface_reflections )  THEN
!
!--     Set nrefsteps to 0 to disable reflections        
        nrefsteps = 0
        surfoutsl = albedo_surf * surfins
        surfoutll = (1._wp - emiss_surf) * surfinl
        surfoutsw = surfoutsw + surfoutsl
        surfoutlw = surfoutlw + surfoutll
     ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--     Next passes - reflections
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     DO refstep = 1, nrefsteps

         surfoutsl = albedo_surf * surfins
!--         for non-transparent surfaces, longwave albedo is 1 - emissivity
         surfoutll = (1._wp - emiss_surf) * surfinl

#if defined( __parallel )
         CALL MPI_AllGatherv(surfoutsl, nsurfl, MPI_REAL, &
             surfouts, nsurfs, surfstart, MPI_REAL, comm2d, ierr)
         CALL MPI_AllGatherv(surfoutll, nsurfl, MPI_REAL, &
             surfoutl, nsurfs, surfstart, MPI_REAL, comm2d, ierr)
#else
         surfouts = surfoutsl
         surfoutl = surfoutll
#endif

!--         reset for next pass input
         surfins = 0._wp
         surfinl = 0._wp

!--         reflected radiation
         DO isvf = 1, nsvfl
             isurf = svfsurf(1, isvf)
             isurfsrc = svfsurf(2, isvf)
             surfins(isurf) = surfins(isurf) + svf(1,isvf) * svf(2,isvf) * surfouts(isurfsrc)
             surfinl(isurf) = surfinl(isurf) + svf(1,isvf) * surfoutl(isurfsrc)
         ENDDO

!--         radiation absorbed by plant canopy
         DO icsf = 1, ncsfl
             ipcgb = csfsurf(1, icsf)
             isurfsrc = csfsurf(2, icsf)
             IF ( isurfsrc == -1 )  CYCLE ! sky->face only in 1st pass, not here

             pcbinsw(ipcgb) = pcbinsw(ipcgb) + csf(1,icsf) * csf(2,icsf) * surfouts(isurfsrc)
         ENDDO

         surfinsw = surfinsw  + surfins
         surfinlw = surfinlw  + surfinl
         surfoutsw = surfoutsw + surfoutsl
         surfoutlw = surfoutlw + surfoutll
!            surfhf = surfinsw + surfinlw - surfoutsw - surfoutlw

     ENDDO

!--  push heat flux absorbed by plant canopy to respective 3D arrays
     IF ( npcbl > 0 )  THEN
         pc_heating_rate(:,:,:) = 0.0_wp
         pc_transpiration_rate(:,:,:) = 0.0_wp
         DO ipcgb = 1, npcbl
                 
             j = pcbl(iy, ipcgb)
             i = pcbl(ix, ipcgb)
             k = pcbl(iz, ipcgb)
!
!--             Following expression equals former kk = k - nzb_s_inner(j,i)
             kk = k - get_topography_top_index_ji( j, i, 's' )  !- lad arrays are defined flat
             pc_heating_rate(kk, j, i) = (pcbinsw(ipcgb) + pcbinlw(ipcgb)) &
                 * pchf_prep(k) * pt(k, j, i) !-- = dT/dt

!             pc_transpiration_rate(kk,j,i) = 0.75_wp* (pcbinsw(ipcgb) + pcbinlw(ipcgb)) &
!                 * pctf_prep(k) * pt(k, j, i) !-- = dq/dt

         ENDDO
     ENDIF
!
!--     Transfer radiation arrays required for energy balance to the respective data types
     DO  i = 1, nsurfl
        m  = surfl(5,i)
!
!--     (1) Urban surfaces
!--     upward-facing
        IF ( surfl(1,i) == iup_u )  THEN
           surf_usm_h%rad_sw_in(m)  = surfinsw(i)
           surf_usm_h%rad_sw_out(m) = surfoutsw(i)
           surf_usm_h%rad_lw_in(m)  = surfinlw(i)
           surf_usm_h%rad_lw_out(m) = surfoutlw(i)
           surf_usm_h%rad_net(m)    = surfinsw(i) - surfoutsw(i) +          &
                                      surfinlw(i) - surfoutlw(i)
!
!--     northward-facding
        ELSEIF ( surfl(1,i) == inorth_u )  THEN
           surf_usm_v(0)%rad_sw_in(m)  = surfinsw(i)
           surf_usm_v(0)%rad_sw_out(m) = surfoutsw(i)
           surf_usm_v(0)%rad_lw_in(m)  = surfinlw(i)
           surf_usm_v(0)%rad_lw_out(m) = surfoutlw(i)
           surf_usm_v(0)%rad_net(m)    = surfinsw(i) - surfoutsw(i) +       &
                                         surfinlw(i) - surfoutlw(i)
!
!--     southward-facding
        ELSEIF ( surfl(1,i) == isouth_u )  THEN
           surf_usm_v(1)%rad_sw_in(m)  = surfinsw(i)
           surf_usm_v(1)%rad_sw_out(m) = surfoutsw(i)
           surf_usm_v(1)%rad_lw_in(m)  = surfinlw(i)
           surf_usm_v(1)%rad_lw_out(m) = surfoutlw(i)
           surf_usm_v(1)%rad_net(m)    = surfinsw(i) - surfoutsw(i) +       &
                                         surfinlw(i) - surfoutlw(i)
!
!--     eastward-facing
        ELSEIF ( surfl(1,i) == ieast_u )  THEN
           surf_usm_v(2)%rad_sw_in(m)  = surfinsw(i)
           surf_usm_v(2)%rad_sw_out(m) = surfoutsw(i)
           surf_usm_v(2)%rad_lw_in(m)  = surfinlw(i)
           surf_usm_v(2)%rad_lw_out(m) = surfoutlw(i)
           surf_usm_v(2)%rad_net(m)    = surfinsw(i) - surfoutsw(i) +       &
                                         surfinlw(i) - surfoutlw(i)
!
!--     westward-facding
        ELSEIF ( surfl(1,i) == iwest_u )  THEN
           surf_usm_v(3)%rad_sw_in(m)  = surfinsw(i)
           surf_usm_v(3)%rad_sw_out(m) = surfoutsw(i)
           surf_usm_v(3)%rad_lw_in(m)  = surfinlw(i)
           surf_usm_v(3)%rad_lw_out(m) = surfoutlw(i)
           surf_usm_v(3)%rad_net(m)    = surfinsw(i) - surfoutsw(i) +       &
                                         surfinlw(i) - surfoutlw(i)
!
!--     (2) land surfaces
!--     upward-facing
        ELSEIF ( surfl(1,i) == iup_l )  THEN
           surf_lsm_h%rad_sw_in(m)  = surfinsw(i)
           surf_lsm_h%rad_sw_out(m) = surfoutsw(i)
           surf_lsm_h%rad_lw_in(m)  = surfinlw(i)
           surf_lsm_h%rad_lw_out(m) = surfoutlw(i)
           surf_lsm_h%rad_net(m)    = surfinsw(i) - surfoutsw(i) +          &
                                      surfinlw(i) - surfoutlw(i)
!
!--     northward-facding
        ELSEIF ( surfl(1,i) == inorth_l )  THEN
           surf_lsm_v(0)%rad_sw_in(m)  = surfinsw(i)
           surf_lsm_v(0)%rad_sw_out(m) = surfoutsw(i)
           surf_lsm_v(0)%rad_lw_in(m)  = surfinlw(i)
           surf_lsm_v(0)%rad_lw_out(m) = surfoutlw(i)
           surf_lsm_v(0)%rad_net(m)    = surfinsw(i) - surfoutsw(i) +       &
                                         surfinlw(i) - surfoutlw(i)
!
!--     southward-facding
        ELSEIF ( surfl(1,i) == isouth_l )  THEN
           surf_lsm_v(1)%rad_sw_in(m)  = surfinsw(i)
           surf_lsm_v(1)%rad_sw_out(m) = surfoutsw(i)
           surf_lsm_v(1)%rad_lw_in(m)  = surfinlw(i)
           surf_lsm_v(1)%rad_lw_out(m) = surfoutlw(i)
           surf_lsm_v(1)%rad_net(m)    = surfinsw(i) - surfoutsw(i) +       &
                                         surfinlw(i) - surfoutlw(i)
!
!--     eastward-facing
        ELSEIF ( surfl(1,i) == ieast_l )  THEN
           surf_lsm_v(2)%rad_sw_in(m)  = surfinsw(i)
           surf_lsm_v(2)%rad_sw_out(m) = surfoutsw(i)
           surf_lsm_v(2)%rad_lw_in(m)  = surfinlw(i)
           surf_lsm_v(2)%rad_lw_out(m) = surfoutlw(i)
           surf_lsm_v(2)%rad_net(m)    = surfinsw(i) - surfoutsw(i) +       &
                                         surfinlw(i) - surfoutlw(i)
!
!--     westward-facing
        ELSEIF ( surfl(1,i) == iwest_l )  THEN
           surf_lsm_v(3)%rad_sw_in(m)  = surfinsw(i)
           surf_lsm_v(3)%rad_sw_out(m) = surfoutsw(i)
           surf_lsm_v(3)%rad_lw_in(m)  = surfinlw(i)
           surf_lsm_v(3)%rad_lw_out(m) = surfoutlw(i)
           surf_lsm_v(3)%rad_net(m)    = surfinsw(i) - surfoutsw(i) +       &
                                         surfinlw(i) - surfoutlw(i)
        ENDIF

     ENDDO

     DO  m = 1, surf_usm_h%ns
        surf_usm_h%surfhf(m) = surf_usm_h%rad_sw_in(m)  +                   &
                               surf_usm_h%rad_lw_in(m)  -                   &
                               surf_usm_h%rad_sw_out(m) -                   &
                               surf_usm_h%rad_lw_out(m)
     ENDDO
     DO  m = 1, surf_lsm_h%ns
        surf_lsm_h%surfhf(m) = surf_lsm_h%rad_sw_in(m)  +                   &
                               surf_lsm_h%rad_lw_in(m)  -                   &
                               surf_lsm_h%rad_sw_out(m) -                   &
                               surf_lsm_h%rad_lw_out(m)
     ENDDO

     DO  l = 0, 3
!--     urban
        DO  m = 1, surf_usm_v(l)%ns
           surf_usm_v(l)%surfhf(m) = surf_usm_v(l)%rad_sw_in(m)  +          &
                                     surf_usm_v(l)%rad_lw_in(m)  -          &
                                     surf_usm_v(l)%rad_sw_out(m) -          &
                                     surf_usm_v(l)%rad_lw_out(m)
        ENDDO
!--     land
        DO  m = 1, surf_lsm_v(l)%ns
           surf_lsm_v(l)%surfhf(m) = surf_lsm_v(l)%rad_sw_in(m)  +          &
                                     surf_lsm_v(l)%rad_lw_in(m)  -          &
                                     surf_lsm_v(l)%rad_sw_out(m) -          &
                                     surf_lsm_v(l)%rad_lw_out(m)

        ENDDO
     ENDDO
!
!--  Calculate the average temperature, albedo, and emissivity for urban/land
!--  domain when using average_radiation in the respective radiation model

!--  Precalculate face areas for all face directions using normal vector
     DO d = 0, nsurf_type
        facearea(d) = 1._wp
        IF ( idir(d) == 0 ) facearea(d) = facearea(d) * dx
        IF ( jdir(d) == 0 ) facearea(d) = facearea(d) * dy
        IF ( kdir(d) == 0 ) facearea(d) = facearea(d) * dz(1)
     ENDDO
!
!--  absorbed/received SW & LW and emitted LW energy of all physical
!--  surfaces (land and urban) in local processor
     pinswl = 0._wp
     pinlwl = 0._wp
     pabsswl = 0._wp
     pabslwl = 0._wp
     pemitlwl = 0._wp
     emiss_sum_surfl = 0._wp
     area_surfl = 0._wp
     DO  i = 1, nsurfl
        d = surfl(id, i)
!--  received SW & LW
        pinswl = pinswl + surfinsw(i) * facearea(d)
        pinlwl = pinlwl + surfinlw(i) * facearea(d)
!--   absorbed SW & LW
        pabsswl = pabsswl + (1._wp - albedo_surf(i)) *                   &
                                                surfinsw(i) * facearea(d)
        pabslwl = pabslwl + emiss_surf(i) * surfinlw(i) * facearea(d)
!--   emitted LW
        pemitlwl = pemitlwl + surfoutlw(i) * facearea(d)
!--   emissivity and area sum
        emiss_sum_surfl = emiss_sum_surfl + emiss_surf(i) * facearea(d)
        area_surfl = area_surfl + facearea(d)
     END DO
!
!--  add the absorbed SW energy by plant canopy
     IF ( npcbl > 0 )  THEN
        pabsswl = pabsswl + SUM(pcbinsw)
        pabslwl = pabslwl + SUM(pcbinlw)
     ENDIF
!
!--  gather all rad flux energy in all processors
#if defined( __parallel )
     CALL MPI_ALLREDUCE( pinswl, pinsw, 1, MPI_REAL, MPI_SUM, comm2d, ierr)
     CALL MPI_ALLREDUCE( pinlwl, pinlw, 1, MPI_REAL, MPI_SUM, comm2d, ierr)
     CALL MPI_ALLREDUCE( pabsswl, pabssw, 1, MPI_REAL, MPI_SUM, comm2d, ierr )
     CALL MPI_ALLREDUCE( pabslwl, pabslw, 1, MPI_REAL, MPI_SUM, comm2d, ierr )
     CALL MPI_ALLREDUCE( pemitlwl, pemitlw, 1, MPI_REAL, MPI_SUM, comm2d, ierr )
     CALL MPI_ALLREDUCE( emiss_sum_surfl, emiss_sum_surf, 1, MPI_REAL, MPI_SUM, comm2d, ierr )
     CALL MPI_ALLREDUCE( area_surfl, area_surf, 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
     pinsw = pinswl
     pinlw = pinlwl
     pabssw = pabsswl
     pabslwl = pabslw
     pemitlwl = pemitlw
     emiss_sum_surf = emiss_sum_surfl
     area_surf = area_surfl
#endif

!--  (1) albedo
     IF ( pinsw /= 0.0_wp )  albedo_urb = 1._wp - pabssw / pinsw

!--  (2) average emmsivity
     IF ( area_surf /= 0.0_wp ) emissivity_urb = emiss_sum_surf / area_surf

!--  (3) temperature
     t_rad_urb = ((pemitlw - pabslw + emissivity_urb*pinlw)/(emissivity_urb*sigma_sb*area_surf))**0.25_wp
     

    CONTAINS

!------------------------------------------------------------------------------!
!> Calculates radiation absorbed by box with given size and LAD.
!>
!> Simulates resol**2 rays (by equally spacing a bounding horizontal square
!> conatining all possible rays that would cross the box) and calculates
!> average transparency per ray. Returns fraction of absorbed radiation flux
!> and area for which this fraction is effective.
!------------------------------------------------------------------------------!
    PURE SUBROUTINE box_absorb(boxsize, resol, dens, uvec, area, absorb)
       IMPLICIT NONE

       REAL(wp), DIMENSION(3), INTENT(in) :: &
            boxsize, &      !< z, y, x size of box in m
            uvec            !< z, y, x unit vector of incoming flux
       INTEGER(iwp), INTENT(in) :: &
            resol           !< No. of rays in x and y dimensions
       REAL(wp), INTENT(in) :: &
            dens            !< box density (e.g. Leaf Area Density)
       REAL(wp), INTENT(out) :: &
            area, &         !< horizontal area for flux absorbtion
            absorb          !< fraction of absorbed flux
       REAL(wp) :: &
            xshift, yshift, &
            xmin, xmax, ymin, ymax, &
            xorig, yorig, &
            dx1, dy1, dz1, dx2, dy2, dz2, &
            crdist, &
            transp
       INTEGER(iwp) :: &
            i, j

       xshift = uvec(3) / uvec(1) * boxsize(1)
       xmin = min(0._wp, -xshift)
       xmax = boxsize(3) + max(0._wp, -xshift)
       yshift = uvec(2) / uvec(1) * boxsize(1)
       ymin = min(0._wp, -yshift)
       ymax = boxsize(2) + max(0._wp, -yshift)

       transp = 0._wp
       DO i = 1, resol
          xorig = xmin + (xmax-xmin) * (i-.5_wp) / resol
          DO j = 1, resol
             yorig = ymin + (ymax-ymin) * (j-.5_wp) / resol

             dz1 = 0._wp
             dz2 = boxsize(1)/uvec(1)

             IF ( uvec(2) > 0._wp )  THEN
                dy1 = -yorig             / uvec(2) !< crossing with y=0
                dy2 = (boxsize(2)-yorig) / uvec(2) !< crossing with y=boxsize(2)
             ELSE !uvec(2)==0
                dy1 = -huge(1._wp)
                dy2 = huge(1._wp)
             ENDIF

             IF ( uvec(3) > 0._wp )  THEN
                dx1 = -xorig             / uvec(3) !< crossing with x=0
                dx2 = (boxsize(3)-xorig) / uvec(3) !< crossing with x=boxsize(3)
             ELSE !uvec(3)==0
                dx1 = -huge(1._wp)
                dx2 = huge(1._wp)
             ENDIF

             crdist = max(0._wp, (min(dz2, dy2, dx2) - max(dz1, dy1, dx1)))
             transp = transp + exp(-ext_coef * dens * crdist)
          ENDDO
       ENDDO
       transp = transp / resol**2
       area = (boxsize(3)+xshift)*(boxsize(2)+yshift)
       absorb = 1._wp - transp

    END SUBROUTINE box_absorb

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine splits direct and diffusion dw radiation
!> It sould not be called in case the radiation model already does it
!> It follows <CITATION>
!------------------------------------------------------------------------------!
    SUBROUTINE calc_diffusion_radiation 
    
        REAL(wp), PARAMETER                          :: lowest_solarUp = 0.1_wp  !< limit the sun elevation to protect stability of the calculation
        INTEGER(iwp)                                 :: i, j
        REAL(wp)                                     ::  year_angle              !< angle
        REAL(wp)                                     ::  etr                     !< extraterestrial radiation
        REAL(wp)                                     ::  corrected_solarUp       !< corrected solar up radiation
        REAL(wp)                                     ::  horizontalETR           !< horizontal extraterestrial radiation
        REAL(wp)                                     ::  clearnessIndex          !< clearness index
        REAL(wp)                                     ::  diff_frac               !< diffusion fraction of the radiation

        
!--     Calculate current day and time based on the initial values and simulation time
        year_angle = ( (day_of_year_init * 86400) + time_utc_init              &
                        + time_since_reference_point )  * d_seconds_year       &
                        * 2.0_wp * pi
        
        etr = solar_constant * (1.00011_wp +                                   &
                          0.034221_wp * cos(year_angle) +                      &
                          0.001280_wp * sin(year_angle) +                      &
                          0.000719_wp * cos(2.0_wp * year_angle) +             &
                          0.000077_wp * sin(2.0_wp * year_angle))
        
!--    
!--     Under a very low angle, we keep extraterestrial radiation at
!--     the last small value, therefore the clearness index will be pushed
!--     towards 0 while keeping full continuity.
!--    
        IF ( zenith(0) <= lowest_solarUp )  THEN
            corrected_solarUp = lowest_solarUp
        ELSE
            corrected_solarUp = zenith(0)
        ENDIF
        
        horizontalETR = etr * corrected_solarUp
        
        DO i = nxl, nxr
            DO j = nys, nyn
                clearnessIndex = rad_sw_in(0,j,i) / horizontalETR
                diff_frac = 1.0_wp / (1.0_wp + exp(-5.0033_wp + 8.6025_wp * clearnessIndex))
                rad_sw_in_diff(j,i) = rad_sw_in(0,j,i) * diff_frac
                rad_sw_in_dir(j,i)  = rad_sw_in(0,j,i) * (1.0_wp - diff_frac)
                rad_lw_in_diff(j,i) = rad_lw_in(0,j,i)
            ENDDO
        ENDDO
        
    END SUBROUTINE calc_diffusion_radiation


 END SUBROUTINE radiation_interaction
    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine initializes structures needed for radiative transfer
!> model. This model calculates transformation processes of the
!> radiation inside urban and land canopy layer. The module includes also
!> the interaction of the radiation with the resolved plant canopy.
!>
!> For more info. see Resler et al. 2017
!>
!> The new version 2.0 was radically rewriten, the discretization scheme
!> has been changed. This new version significantly improves effectivity
!> of the paralelization and the scalability of the model.
!>
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_interaction_init

       USE control_parameters,                                                 &
           ONLY:  dz_stretch_level_start
           
       USE netcdf_data_input_mod,                                              &
           ONLY:  leaf_area_density_f

       USE plant_canopy_model_mod,                                             &
           ONLY:  pch_index, pc_heating_rate, lad_s

       IMPLICIT NONE

       INTEGER(iwp) :: i, j, k, d, l, ir, jr, ids, m
       INTEGER(iwp) :: k_topo     !< vertical index indicating topography top for given (j,i)
       INTEGER(iwp) :: k_topo2    !< vertical index indicating topography top for given (j,i)
       INTEGER(iwp) :: nzptl, nzubl, nzutl, isurf, ipcgb
       INTEGER(iwp) :: procid
       REAL(wp)     :: mrl


       !INTEGER(iwp), DIMENSION(1:4,inorth_b:iwest_b)  ::  ijdb                               !< start and end of the local domain border coordinates (set in code)
       !LOGICAL, DIMENSION(inorth_b:iwest_b)           ::  isborder                           !< is PE on the border of the domain in four corresponding directions

!
!--    Find nzub, nzut, nzu via wall_flag_0 array (nzb_s_inner will be
!--    removed later). The following contruct finds the lowest / largest index
!--    for any upward-facing wall (see bit 12).
       nzubl = MINVAL( get_topography_top_index( 's' ) )
       nzutl = MAXVAL( get_topography_top_index( 's' ) )

       nzubl = MAX( nzubl, nzb )

       IF ( plant_canopy )  THEN
!--        allocate needed arrays
           ALLOCATE( pct(nys:nyn,nxl:nxr) )
           ALLOCATE( pch(nys:nyn,nxl:nxr) )

!--        calculate plant canopy height
           npcbl = 0
           pct   = 0
           pch   = 0
           DO i = nxl, nxr
               DO j = nys, nyn
!
!--                Find topography top index
                   k_topo = get_topography_top_index_ji( j, i, 's' )

                   DO k = nzt+1, 0, -1
                       IF ( lad_s(k,j,i) /= 0.0_wp )  THEN
!--                        we are at the top of the pcs
                           pct(j,i) = k + k_topo
                           pch(j,i) = k
                           npcbl = npcbl + pch(j,i)
                           EXIT
                       ENDIF
                   ENDDO
               ENDDO
           ENDDO

           nzutl = MAX( nzutl, MAXVAL( pct ) )
           nzptl = MAXVAL( pct )
!--        code of plant canopy model uses parameter pch_index
!--        we need to setup it here to right value
!--        (pch_index, lad_s and other arrays in PCM are defined flat)
           pch_index = MERGE( leaf_area_density_f%nz - 1, MAXVAL( pch ),       &
                              leaf_area_density_f%from_file )

           prototype_lad = MAXVAL( lad_s ) * .9_wp  !< better be *1.0 if lad is either 0 or maxval(lad) everywhere
           IF ( prototype_lad <= 0._wp ) prototype_lad = .3_wp
           !WRITE(message_string, '(a,f6.3)') 'Precomputing effective box optical ' &
           !    // 'depth using prototype leaf area density = ', prototype_lad
           !CALL message('usm_init_urban_surface', 'PA0520', 0, 0, -1, 6, 0)
       ENDIF

       nzutl = MIN( nzutl + nzut_free, nzt )

#if defined( __parallel )
       CALL MPI_AllReduce(nzubl, nzub, 1, MPI_INTEGER, MPI_MIN, comm2d, ierr )
       CALL MPI_AllReduce(nzutl, nzut, 1, MPI_INTEGER, MPI_MAX, comm2d, ierr )
       CALL MPI_AllReduce(nzptl, nzpt, 1, MPI_INTEGER, MPI_MAX, comm2d, ierr )
#else
       nzub = nzubl
       nzut = nzutl
       nzpt = nzptl
#endif
!
!--    Stretching (non-uniform grid spacing) is not considered in the radiation 
!--    model. Therefore, vertical stretching has to be applied above the area 
!--    where the parts of the radiation model which assume constant grid spacing
!--    are active. ABS (...) is required because the default value of 
!--    dz_stretch_level_start is -9999999.9_wp (negative). 
       IF ( ABS( dz_stretch_level_start(1) ) <= zw(nzut) ) THEN
          WRITE( message_string, * ) 'The lowest level where vertical ',       &
                                     'stretching is applied have to be ',      &
                                     'greater than ', zw(nzut)
          CALL message( 'radiation_interaction_init', 'PA0496', 1, 2, 0, 6, 0 )
       ENDIF 
!
!--    global number of urban and plant layers
       nzu = nzut - nzub + 1
       nzp = nzpt - nzub + 1
!
!--    check max_raytracing_dist relative to urban surface layer height
!        mrl = 2.0_wp * nzu * dz
!        IF ( max_raytracing_dist <= mrl ) THEN
!           IF ( max_raytracing_dist /= -999.0_wp ) THEN
! !--          max_raytracing_dist too low
!              WRITE(message_string, '(a,f6.1)') 'Max_raytracing_dist too low, ' &
!                    // 'override to value ', mrl
!              CALL message('radiation_interaction_init', 'PA0521', 0, 0, -1, 6, 0)
!           ENDIF
!           max_raytracing_dist = mrl
!        ENDIF
!
!--    allocate urban surfaces grid
!--    calc number of surfaces in local proc
       CALL location_message( '    calculation of indices for surfaces', .TRUE. )
       nsurfl = 0
!
!--    Number of horizontal surfaces including land- and roof surfaces in both USM and LSM. Note that
!--    All horizontal surface elements are already counted in surface_mod.
       startland = 1
       nsurfl    = surf_usm_h%ns + surf_lsm_h%ns
       endland   = nsurfl
       nlands    = endland - startland + 1

!
!--    Number of vertical surfaces in both USM and LSM. Note that all vertical surface elements are
!--    already counted in surface_mod.
       startwall = nsurfl+1
       DO  i = 0,3
          nsurfl = nsurfl + surf_usm_v(i)%ns + surf_lsm_v(i)%ns
       ENDDO
       endwall = nsurfl
       nwalls  = endwall - startwall + 1

!--    fill gridpcbl and pcbl
       IF ( npcbl > 0 )  THEN
           ALLOCATE( pcbl(iz:ix, 1:npcbl) )
           ALLOCATE( gridpcbl(nzub:nzpt,nys:nyn,nxl:nxr) )
           pcbl = -1
           gridpcbl(:,:,:) = 0
           ipcgb = 0
           DO i = nxl, nxr
               DO j = nys, nyn
!
!--                Find topography top index
                   k_topo = get_topography_top_index_ji( j, i, 's' )

                   DO k = k_topo + 1, pct(j,i)
                       ipcgb = ipcgb + 1
                       gridpcbl(k,j,i) = ipcgb
                       pcbl(:,ipcgb) = (/ k, j, i /)
                   ENDDO
               ENDDO
           ENDDO
           ALLOCATE( pcbinsw( 1:npcbl ) )
           ALLOCATE( pcbinswdir( 1:npcbl ) )
           ALLOCATE( pcbinswdif( 1:npcbl ) )
           ALLOCATE( pcbinlw( 1:npcbl ) )
       ENDIF

!--    fill surfl (the ordering of local surfaces given by the following
!--    cycles must not be altered, certain file input routines may depend
!--    on it)
       ALLOCATE(surfl(5,nsurfl))  ! is it mecessary to allocate it with (5,nsurfl)?
       isurf = 0

!--    add horizontal surface elements (land and urban surfaces)
!--    TODO: add urban overhanging surfaces (idown_u)
       DO i = nxl, nxr
           DO j = nys, nyn
              DO  m = surf_usm_h%start_index(j,i), surf_usm_h%end_index(j,i)
                 k = surf_usm_h%k(m)

                 isurf = isurf + 1
                 surfl(:,isurf) = (/iup_u,k,j,i,m/)
              ENDDO

              DO  m = surf_lsm_h%start_index(j,i), surf_lsm_h%end_index(j,i)
                 k = surf_lsm_h%k(m)

                 isurf = isurf + 1
                 surfl(:,isurf) = (/iup_l,k,j,i,m/)
              ENDDO

           ENDDO
       ENDDO

!--    add vertical surface elements (land and urban surfaces)
!--    TODO: remove the hard coding of l = 0 to l = idirection
       DO i = nxl, nxr
           DO j = nys, nyn
              l = 0
              DO  m = surf_usm_v(l)%start_index(j,i), surf_usm_v(l)%end_index(j,i)
                 k = surf_usm_v(l)%k(m)

                 isurf          = isurf + 1
                 surfl(:,isurf) = (/inorth_u,k,j,i,m/)
              ENDDO
              DO  m = surf_lsm_v(l)%start_index(j,i), surf_lsm_v(l)%end_index(j,i)
                 k = surf_lsm_v(l)%k(m)

                 isurf          = isurf + 1
                 surfl(:,isurf) = (/inorth_l,k,j,i,m/)
              ENDDO

              l = 1
              DO  m = surf_usm_v(l)%start_index(j,i), surf_usm_v(l)%end_index(j,i)
                 k = surf_usm_v(l)%k(m)

                 isurf          = isurf + 1
                 surfl(:,isurf) = (/isouth_u,k,j,i,m/)
              ENDDO
              DO  m = surf_lsm_v(l)%start_index(j,i), surf_lsm_v(l)%end_index(j,i)
                 k = surf_lsm_v(l)%k(m)

                 isurf          = isurf + 1
                 surfl(:,isurf) = (/isouth_l,k,j,i,m/)
              ENDDO

              l = 2
              DO  m = surf_usm_v(l)%start_index(j,i), surf_usm_v(l)%end_index(j,i)
                 k = surf_usm_v(l)%k(m)

                 isurf          = isurf + 1
                 surfl(:,isurf) = (/ieast_u,k,j,i,m/)
              ENDDO
              DO  m = surf_lsm_v(l)%start_index(j,i), surf_lsm_v(l)%end_index(j,i)
                 k = surf_lsm_v(l)%k(m)

                 isurf          = isurf + 1
                 surfl(:,isurf) = (/ieast_l,k,j,i,m/)
              ENDDO

              l = 3
              DO  m = surf_usm_v(l)%start_index(j,i), surf_usm_v(l)%end_index(j,i)
                 k = surf_usm_v(l)%k(m)

                 isurf          = isurf + 1
                 surfl(:,isurf) = (/iwest_u,k,j,i,m/)
              ENDDO
              DO  m = surf_lsm_v(l)%start_index(j,i), surf_lsm_v(l)%end_index(j,i)
                 k = surf_lsm_v(l)%k(m)

                 isurf          = isurf + 1
                 surfl(:,isurf) = (/iwest_l,k,j,i,m/)
              ENDDO
           ENDDO
       ENDDO

!
!--    broadband albedo of the land, roof and wall surface
!--    for domain border and sky set artifically to 1.0
!--    what allows us to calculate heat flux leaving over
!--    side and top borders of the domain
       ALLOCATE ( albedo_surf(nsurfl) )
       albedo_surf = 1.0_wp
!
!--    Also allocate further array for emissivity with identical order of
!--    surface elements as radiation arrays.
       ALLOCATE ( emiss_surf(nsurfl)  )


!
!--    global array surf of indices of surfaces and displacement index array surfstart
       ALLOCATE(nsurfs(0:numprocs-1))

#if defined( __parallel )
       CALL MPI_Allgather(nsurfl,1,MPI_INTEGER,nsurfs,1,MPI_INTEGER,comm2d,ierr)
#else
       nsurfs(0) = nsurfl
#endif
       ALLOCATE(surfstart(0:numprocs))
       k = 0
       DO i=0,numprocs-1
           surfstart(i) = k
           k = k+nsurfs(i)
       ENDDO
       surfstart(numprocs) = k
       nsurf = k
       ALLOCATE(surf(5,nsurf))

#if defined( __parallel )
       CALL MPI_AllGatherv(surfl, nsurfl*5, MPI_INTEGER, surf, nsurfs*5, &
           surfstart(0:numprocs-1)*5, MPI_INTEGER, comm2d, ierr)
#else
       surf = surfl
#endif

!--
!--    allocation of the arrays for direct and diffusion radiation
       CALL location_message( '    allocation of radiation arrays', .TRUE. )
!--    rad_sw_in, rad_lw_in are computed in radiation model,
!--    splitting of direct and diffusion part is done
!--    in calc_diffusion_radiation for now

       ALLOCATE( rad_sw_in_dir(nysg:nyng,nxlg:nxrg) )
       ALLOCATE( rad_sw_in_diff(nysg:nyng,nxlg:nxrg) )
       ALLOCATE( rad_lw_in_diff(nysg:nyng,nxlg:nxrg) )
       rad_sw_in_dir  = 0.0_wp
       rad_sw_in_diff = 0.0_wp
       rad_lw_in_diff = 0.0_wp

!--    allocate radiation arrays
       ALLOCATE( surfins(nsurfl) )
       ALLOCATE( surfinl(nsurfl) )
       ALLOCATE( surfinsw(nsurfl) )
       ALLOCATE( surfinlw(nsurfl) )
       ALLOCATE( surfinswdir(nsurfl) )
       ALLOCATE( surfinswdif(nsurfl) )
       ALLOCATE( surfinlwdif(nsurfl) )
       ALLOCATE( surfoutsl(nsurfl) )
       ALLOCATE( surfoutll(nsurfl) )
       ALLOCATE( surfoutsw(nsurfl) )
       ALLOCATE( surfoutlw(nsurfl) )
       ALLOCATE( surfouts(nsurf) )
       ALLOCATE( surfoutl(nsurf) )
       ALLOCATE( skyvf(nsurfl) )
       ALLOCATE( skyvft(nsurfl) )

!
!--    In case of average_radiation, aggregated surface albedo and emissivity,
!--    also set initial value for t_rad_urb.
!--    For now set an arbitrary initial value.
       IF ( average_radiation )  THEN
          albedo_urb = 0.5_wp
          emissivity_urb = 0.5_wp
          t_rad_urb = pt_surface
       ENDIF

    END SUBROUTINE radiation_interaction_init

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates shape view factors (SVF), plant sink canopy factors (PCSF),
!> sky-view factors, discretized path for direct solar radiation, MRT factors
!> and other preprocessed data needed for radiation_interaction.
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_calc_svf
    
        IMPLICIT NONE
        
        INTEGER(iwp)                                  :: i, j, k, l, d, ip, jp
        INTEGER(iwp)                                  :: isvf, ksvf, icsf, kcsf, npcsfl, isvf_surflt, imrtt, imrtf, ipcgb
        INTEGER(iwp)                                  :: sd, td, ioln, iproc
        INTEGER(iwp)                                  :: iaz, izn      !< azimuth, zenith counters
        INTEGER(iwp)                                  :: naz, nzn      !< azimuth, zenith num of steps
        REAL(wp)                                      :: az0, zn0      !< starting azimuth/zenith
        REAL(wp)                                      :: azs, zns      !< azimuth/zenith cycle step
        REAL(wp)                                      :: az1, az2      !< relative azimuth of section borders
        REAL(wp)                                      :: azmid         !< ray (center) azimuth
        REAL(wp)                                      :: horizon       !< computed horizon height (tangent of elevation)
        REAL(wp)                                      :: azen          !< zenith angle
        REAL(wp), DIMENSION(:), ALLOCATABLE           :: zdirs         !< directions in z (tangent of elevation)
        REAL(wp), DIMENSION(:), ALLOCATABLE           :: zbdry         !< zenith angle boundaries
        REAL(wp), DIMENSION(:), ALLOCATABLE           :: vffrac        !< view factor fractions for individual rays
        REAL(wp), DIMENSION(:), ALLOCATABLE           :: ztransp       !< array of transparency in z steps
        REAL(wp),     DIMENSION(0:nsurf_type)         :: facearea
        INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE     :: nzterrl
        REAL(wp),     DIMENSION(:,:), ALLOCATABLE     :: csflt, pcsflt
        INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE     :: kcsflt,kpcsflt
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE       :: icsflt,dcsflt,ipcsflt,dpcsflt
        REAL(wp), DIMENSION(3)                        :: uv
        LOGICAL                                       :: visible
        REAL(wp), DIMENSION(3)                        :: sa, ta          !< real coordinates z,y,x of source and target
        REAL(wp)                                      :: transparency, rirrf, sqdist, svfsum
        INTEGER(iwp)                                  :: isurflt, isurfs, isurflt_prev
        INTEGER(iwp)                                  :: itx, ity, itz
        INTEGER(idp)                                  :: ray_skip_maxdist, ray_skip_minval !< skipped raytracing counts
        INTEGER(iwp)                                  :: max_track_len !< maximum 2d track length
        CHARACTER(len=7)                              :: pid_char = ''
        INTEGER(iwp)                                  :: win_lad, minfo
        REAL(wp), DIMENSION(:,:,:), POINTER           :: lad_s_rma       !< fortran pointer, but lower bounds are 1
        TYPE(c_ptr)                                   :: lad_s_rma_p     !< allocated c pointer
#if defined( __parallel )
        INTEGER(kind=MPI_ADDRESS_KIND)                :: size_lad_rma
#endif
!    
        INTEGER(iwp), DIMENSION(0:svfnorm_report_num) :: svfnorm_counts
        CHARACTER(200)                                :: msg

!--     calculation of the SVF
        CALL location_message( '    calculation of SVF and CSF', .TRUE. )
!         CALL radiation_write_debug_log('Start calculation of SVF and CSF')

!--     precalculate face areas for different face directions using normal vector
        DO d = 0, nsurf_type
            facearea(d) = 1._wp
            IF ( idir(d) == 0 ) facearea(d) = facearea(d) * dx
            IF ( jdir(d) == 0 ) facearea(d) = facearea(d) * dy
            IF ( kdir(d) == 0 ) facearea(d) = facearea(d) * dz(1)
        ENDDO

!--     initialize variables and temporary arrays for calculation of svf and csf
        nsvfl  = 0
        ncsfl  = 0
        nsvfla = gasize
        msvf   = 1
        ALLOCATE( asvf1(nsvfla) )
        asvf => asvf1
        IF ( plant_canopy )  THEN
            ncsfla = gasize
            mcsf   = 1
            ALLOCATE( acsf1(ncsfla) )
            acsf => acsf1
        ENDIF
        ray_skip_maxdist = 0
        ray_skip_minval = 0
        
!--     initialize temporary terrain and plant canopy height arrays (global 2D array!)
        ALLOCATE( nzterr(0:(nx+1)*(ny+1)-1) )
#if defined( __parallel )
        ALLOCATE( nzterrl(nys:nyn,nxl:nxr) )
        nzterrl = get_topography_top_index( 's' )
        CALL MPI_AllGather( nzterrl, nnx*nny, MPI_INTEGER, &
                            nzterr, nnx*nny, MPI_INTEGER, comm2d, ierr )
        DEALLOCATE(nzterrl)
#else
        nzterr = RESHAPE( get_topography_top_index( 's' ), (/(nx+1)*(ny+1)/) )
#endif
        IF ( plant_canopy )  THEN
            ALLOCATE( plantt(0:(nx+1)*(ny+1)-1) )
            maxboxesg = nx + ny + nzp + 1
            max_track_len = nx + ny + 1
!--         temporary arrays storing values for csf calculation during raytracing
            ALLOCATE( boxes(3, maxboxesg) )
            ALLOCATE( crlens(maxboxesg) )

#if defined( __parallel )
            CALL MPI_AllGather( pct, nnx*nny, MPI_INTEGER, &
                                plantt, nnx*nny, MPI_INTEGER, comm2d, ierr )
            
!--         temporary arrays storing values for csf calculation during raytracing
            ALLOCATE( lad_ip(maxboxesg) )
            ALLOCATE( lad_disp(maxboxesg) )

            IF ( rma_lad_raytrace )  THEN
                ALLOCATE( lad_s_ray(maxboxesg) )
                
                ! set conditions for RMA communication
                CALL MPI_Info_create(minfo, ierr)
                CALL MPI_Info_set(minfo, 'accumulate_ordering', '', ierr)
                CALL MPI_Info_set(minfo, 'accumulate_ops', 'same_op', ierr)
                CALL MPI_Info_set(minfo, 'same_size', 'true', ierr)
                CALL MPI_Info_set(minfo, 'same_disp_unit', 'true', ierr)

!--             Allocate and initialize the MPI RMA window
!--             must be in accordance with allocation of lad_s in plant_canopy_model
!--             optimization of memory should be done
!--             Argument X of function STORAGE_SIZE(X) needs arbitrary REAL(wp) value, set to 1.0_wp for now
                size_lad_rma = STORAGE_SIZE(1.0_wp)/8*nnx*nny*nzp
                CALL MPI_Win_allocate(size_lad_rma, STORAGE_SIZE(1.0_wp)/8, minfo, comm2d, &
                                        lad_s_rma_p, win_lad, ierr)
                CALL c_f_pointer(lad_s_rma_p, lad_s_rma, (/ nzp, nny, nnx /))
                sub_lad(nzub:, nys:, nxl:) => lad_s_rma(:,:,:)
            ELSE
                ALLOCATE(sub_lad(nzub:nzpt, nys:nyn, nxl:nxr))
            ENDIF
#else
            plantt = RESHAPE( pct(nys:nyn,nxl:nxr), (/(nx+1)*(ny+1)/) )
            ALLOCATE(sub_lad(nzub:nzpt, nys:nyn, nxl:nxr))
#endif
            plantt_max = MAXVAL(plantt)
            ALLOCATE( rt2_track(2, max_track_len), rt2_track_lad(nzub:plantt_max, max_track_len), &
                      rt2_track_dist(0:max_track_len), rt2_dist(plantt_max-nzub+2) )

            sub_lad(:,:,:) = 0._wp
            DO i = nxl, nxr
                DO j = nys, nyn
                    k = get_topography_top_index_ji( j, i, 's' )

                    sub_lad(k:nzpt, j, i) = lad_s(0:nzpt-k, j, i)
                ENDDO
            ENDDO

#if defined( __parallel )
            IF ( rma_lad_raytrace )  THEN
                CALL MPI_Info_free(minfo, ierr)
                CALL MPI_Win_lock_all(0, win_lad, ierr)
            ELSE
                ALLOCATE( sub_lad_g(0:(nx+1)*(ny+1)*nzp-1) )
                CALL MPI_AllGather( sub_lad, nnx*nny*nzp, MPI_REAL, &
                                    sub_lad_g, nnx*nny*nzp, MPI_REAL, comm2d, ierr )
            ENDIF
#endif
        ENDIF

        IF ( mrt_factors )  THEN
            OPEN(153, file='MRT_TARGETS', access='SEQUENTIAL', &
                    action='READ', status='OLD', form='FORMATTED', err=524)
            OPEN(154, file='MRT_FACTORS'//myid_char, access='DIRECT', recl=(5*4+2*8), &
                    action='WRITE', status='REPLACE', form='UNFORMATTED', err=525)
            imrtf = 1
            DO
                READ(153, *, end=526, err=524) imrtt, i, j, k
                IF ( i < nxl  .OR.  i > nxr &
                     .OR.  j < nys  .OR.  j > nyn ) CYCLE
                ta = (/ REAL(k), REAL(j), REAL(i) /)

                DO isurfs = 1, nsurf
                    IF ( .NOT.  surface_facing(i, j, k, -1, &
                        surf(ix, isurfs), surf(iy, isurfs), &
                        surf(iz, isurfs), surf(id, isurfs)) )  THEN
                        CYCLE
                    ENDIF
                      
                    sd = surf(id, isurfs)
                    sa = (/ REAL(surf(iz, isurfs), wp) - 0.5_wp * kdir(sd), &
                            REAL(surf(iy, isurfs), wp) - 0.5_wp * jdir(sd), &
                            REAL(surf(ix, isurfs), wp) - 0.5_wp * idir(sd) /)

!--                 unit vector source -> target
                    uv = (/ (ta(1)-sa(1))*dz(1), (ta(2)-sa(2))*dy, (ta(3)-sa(3))*dx /)
                    sqdist = SUM(uv(:)**2)
                    uv = uv / SQRT(sqdist)

!--                 irradiance factor - see svf. Here we consider that target face is always normal,
!--                 i.e. the second dot product equals 1
                    rirrf = dot_product((/ kdir(sd), jdir(sd), idir(sd) /), uv) &
                        / (pi * sqdist) * facearea(sd)

!--                 raytrace while not creating any canopy sink factors
                    CALL raytrace(sa, ta, isurfs, rirrf, 1._wp, .FALSE., &
                            visible, transparency, win_lad)
                    IF ( .NOT.  visible ) CYCLE

                    !rsvf = rirrf * transparency
                    WRITE(154, rec=imrtf, err=525) INT(imrtt, kind=4), &
                        INT(surf(id, isurfs), kind=4), &
                        INT(surf(iz, isurfs), kind=4), &
                        INT(surf(iy, isurfs), kind=4), &
                        INT(surf(ix, isurfs), kind=4), &
                        REAL(rirrf, kind=8), REAL(transparency, kind=8)
                    imrtf = imrtf + 1

                ENDDO !< isurfs
            ENDDO !< MRT_TARGETS record

524         message_string = 'error reading file MRT_TARGETS'
            CALL message( 'radiation_calc_svf', 'PA0524', 1, 2, 0, 6, 0 )

525         message_string = 'error writing file MRT_FACTORS'//myid_char
            CALL message( 'radiation_calc_svf', 'PA0525', 1, 2, 0, 6, 0 )

526         CLOSE(153)
            CLOSE(154)
        ENDIF  !< mrt_factors
        
        !--Directions opposite to face normals are not even calculated,
        !--they must be preset to 0
        !--
        dsitrans(:,:) = 0._wp
        
        DO isurflt = 1, nsurfl
!--         determine face centers
            td = surfl(id, isurflt)
            ta = (/ REAL(surfl(iz, isurflt), wp) - 0.5_wp * kdir(td),  &
                      REAL(surfl(iy, isurflt), wp) - 0.5_wp * jdir(td),  &
                      REAL(surfl(ix, isurflt), wp) - 0.5_wp * idir(td)  /)

            !--Calculate sky view factor and raytrace DSI paths
            skyvf(isurflt) = 0._wp
            skyvft(isurflt) = 0._wp

            !--Select a proper half-sphere for 2D raytracing
            SELECT CASE ( td )
               CASE ( iup_u, iup_l )
                  az0 = 0._wp
                  naz = raytrace_discrete_azims
                  azs = 2._wp * pi / REAL(naz, wp)
                  zn0 = 0._wp
                  nzn = raytrace_discrete_elevs / 2
                  zns = pi / 2._wp / REAL(nzn, wp)
               CASE ( isouth_u, isouth_l )
                  az0 = pi / 2._wp
                  naz = raytrace_discrete_azims / 2
                  azs = pi / REAL(naz, wp)
                  zn0 = 0._wp
                  nzn = raytrace_discrete_elevs
                  zns = pi / REAL(nzn, wp)
               CASE ( inorth_u, inorth_l )
                  az0 = - pi / 2._wp
                  naz = raytrace_discrete_azims / 2
                  azs = pi / REAL(naz, wp)
                  zn0 = 0._wp
                  nzn = raytrace_discrete_elevs
                  zns = pi / REAL(nzn, wp)
               CASE ( iwest_u, iwest_l )
                  az0 = pi
                  naz = raytrace_discrete_azims / 2
                  azs = pi / REAL(naz, wp)
                  zn0 = 0._wp
                  nzn = raytrace_discrete_elevs
                  zns = pi / REAL(nzn, wp)
               CASE ( ieast_u, ieast_l )
                  az0 = 0._wp
                  naz = raytrace_discrete_azims / 2
                  azs = pi / REAL(naz, wp)
                  zn0 = 0._wp
                  nzn = raytrace_discrete_elevs
                  zns = pi / REAL(nzn, wp)
               CASE DEFAULT
                  WRITE(message_string, *) 'ERROR: the surface type ', td,     &
                                           ' is not supported for calculating',&
                                           ' SVF'
                  CALL message( 'radiation_calc_svf', 'PA0488', 1, 2, 0, 6, 0 )
            END SELECT

            ALLOCATE ( zdirs(1:nzn), zbdry(0:nzn), vffrac(1:nzn), ztransp(1:nzn) )
            zdirs(:) = (/( TAN(pi/2 - (zn0+(REAL(izn,wp)-.5_wp)*zns)), izn=1, nzn )/)
            zbdry(:) = (/( zn0+REAL(izn,wp)*zns, izn=0, nzn )/)
            IF ( td == iup_u  .OR.  td == iup_l )  THEN
               !-- For horizontal target, vf fractions are constant per azimuth
               vffrac(:) = (COS(2 * zbdry(0:nzn-1)) - COS(2 * zbdry(1:nzn))) / 2._wp / REAL(naz, wp)
               !--sum of vffrac for all iaz equals 1, verified
            ENDIF

            !--Calculate sky-view factor and direct solar visibility using 2D raytracing
            DO iaz = 1, naz
               azmid = az0 + (REAL(iaz, wp) - .5_wp) * azs
               IF ( td /= iup_u  .AND.  td /= iup_l )  THEN
                  az2 = REAL(iaz, wp) * azs - pi/2._wp
                  az1 = az2 - azs
                  !TODO precalculate after 1st line
                  vffrac(:) = (SIN(az2) - SIN(az1))                           &
                              * (zbdry(1:nzn) - zbdry(0:nzn-1)                &
                                 + SIN(zbdry(0:nzn-1))*COS(zbdry(0:nzn-1))    &
                                 - SIN(zbdry(1:nzn))*COS(zbdry(1:nzn)))       &
                              / (2._wp * pi)
                  !--sum of vffrac for all iaz equals 1, verified
               ENDIF
               CALL raytrace_2d(ta, (/ COS(azmid), SIN(azmid) /), zdirs,      &
                                    surfstart(myid) + isurflt, facearea(td),  &
                                    vffrac, .TRUE., .FALSE., win_lad, horizon,&
                                    ztransp) !FIXME unit vect in grid units + zdirs

               azen = pi/2 - ATAN(horizon)
               IF ( td == iup_u  .OR.  td == iup_l )  THEN
                  azen = MIN(azen, pi/2) !only above horizontal direction
                  skyvf(isurflt) = skyvf(isurflt) + (1._wp - COS(2*azen)) /   &
                     (2._wp * raytrace_discrete_azims)
               ELSE
                  skyvf(isurflt) = skyvf(isurflt) + (SIN(az2) - SIN(az1)) *   &
                              (azen - SIN(azen)*COS(azen)) / (2._wp*pi)
               ENDIF
               skyvft(isurflt) = skyvft(isurflt) + SUM(ztransp(:) * vffrac(:))
 
               !--Save direct solar transparency
               j = MODULO(NINT(azmid/                                          &
                               (2._wp*pi)*raytrace_discrete_azims-.5_wp, iwp), &
                          raytrace_discrete_azims)

               DO k = 1, raytrace_discrete_elevs/2
                  i = dsidir_rev(k-1, j)
                  IF ( i /= -1 )  dsitrans(isurflt, i) = ztransp(k)
               ENDDO
            ENDDO

            DEALLOCATE ( zdirs, zbdry, vffrac, ztransp )
!
!--         Following calculations only required for surface_reflections
            IF ( surface_reflections )  THEN

               DO  isurfs = 1, nsurf
                  IF ( .NOT.  surface_facing(surfl(ix, isurflt), surfl(iy, isurflt), &
                     surfl(iz, isurflt), surfl(id, isurflt), &
                     surf(ix, isurfs), surf(iy, isurfs), &
                     surf(iz, isurfs), surf(id, isurfs)) )  THEN
                     CYCLE
                  ENDIF
                  
                  sd = surf(id, isurfs)
                  sa = (/ REAL(surf(iz, isurfs), wp) - 0.5_wp * kdir(sd),  &
                          REAL(surf(iy, isurfs), wp) - 0.5_wp * jdir(sd),  &
                          REAL(surf(ix, isurfs), wp) - 0.5_wp * idir(sd)  /)

!--               unit vector source -> target
                  uv = (/ (ta(1)-sa(1))*dz(1), (ta(2)-sa(2))*dy, (ta(3)-sa(3))*dx /)
                  sqdist = SUM(uv(:)**2)
                  uv = uv / SQRT(sqdist)

!--               reject raytracing above max distance
                  IF ( SQRT(sqdist) > max_raytracing_dist ) THEN
                     ray_skip_maxdist = ray_skip_maxdist + 1
                     CYCLE
                  ENDIF
                 
!--               irradiance factor (our unshaded shape view factor) = view factor per differential target area * source area
                  rirrf = dot_product((/ kdir(sd), jdir(sd), idir(sd) /), uv) & ! cosine of source normal and direction
                      * dot_product((/ kdir(td), jdir(td), idir(td) /), -uv) &  ! cosine of target normal and reverse direction
                      / (pi * sqdist) & ! square of distance between centers
                      * facearea(sd)

!--               reject raytracing for potentially too small view factor values
                  IF ( rirrf < min_irrf_value ) THEN
                      ray_skip_minval = ray_skip_minval + 1
                      CYCLE
                  ENDIF

!--               raytrace + process plant canopy sinks within
                  CALL raytrace(sa, ta, isurfs, rirrf, facearea(td), .TRUE., &
                                visible, transparency, win_lad)

                  IF ( .NOT.  visible ) CYCLE
                 ! rsvf = rirrf * transparency

!--               write to the svf array
                  nsvfl = nsvfl + 1
!--               check dimmension of asvf array and enlarge it if needed
                  IF ( nsvfla < nsvfl )  THEN
                     k = nsvfla * 2
                     IF ( msvf == 0 )  THEN
                        msvf = 1
                        ALLOCATE( asvf1(k) )
                        asvf => asvf1
                        asvf1(1:nsvfla) = asvf2
                        DEALLOCATE( asvf2 )
                     ELSE
                        msvf = 0
                        ALLOCATE( asvf2(k) )
                        asvf => asvf2
                        asvf2(1:nsvfla) = asvf1
                        DEALLOCATE( asvf1 )
                     ENDIF

!                      WRITE(msg,'(A,3I12)') 'Grow asvf:',nsvfl,nsvfla,k
!                      CALL radiation_write_debug_log( msg )
                     
                     nsvfla = k
                  ENDIF
!--               write svf values into the array
                  asvf(nsvfl)%isurflt = isurflt
                  asvf(nsvfl)%isurfs = isurfs
                  asvf(nsvfl)%rsvf = rirrf !we postopne multiplication by transparency
                  asvf(nsvfl)%rtransp = transparency !a.k.a. Direct Irradiance Factor
               ENDDO
            ENDIF
        ENDDO

        !--Raytrace to canopy boxes to fill dsitransc TODO optimize
        !--
        dsitransc(:,:) = -999._wp !FIXME
        az0 = 0._wp
        naz = raytrace_discrete_azims
        azs = 2._wp * pi / REAL(naz, wp)
        zn0 = 0._wp
        nzn = raytrace_discrete_elevs / 2
        zns = pi / 2._wp / REAL(nzn, wp)
        ALLOCATE ( zdirs(1:nzn), vffrac(1:nzn), ztransp(1:nzn) )
        zdirs(:) = (/( TAN(pi/2 - (zn0+(REAL(izn,wp)-.5_wp)*zns)), izn=1, nzn )/)
        vffrac(:) = 0._wp

        DO ipcgb = 1, npcbl
           ta = (/ REAL(pcbl(iz, ipcgb), wp),  &
                   REAL(pcbl(iy, ipcgb), wp),  &
                   REAL(pcbl(ix, ipcgb), wp) /)
           !--Calculate sky-view factor and direct solar visibility using 2D raytracing
           DO iaz = 1, naz
              azmid = az0 + (REAL(iaz, wp) - .5_wp) * azs
              CALL raytrace_2d(ta, (/ COS(azmid), SIN(azmid) /), zdirs,     &
                                   -999, -999._wp, vffrac, .FALSE., .TRUE., &
                                   win_lad, horizon, ztransp) !FIXME unit vect in grid units + zdirs

              !--Save direct solar transparency
              j = MODULO(NINT(azmid/                                         &
                             (2._wp*pi)*raytrace_discrete_azims-.5_wp, iwp), &
                         raytrace_discrete_azims)
              DO k = 1, raytrace_discrete_elevs/2
                 i = dsidir_rev(k-1, j)
                 IF ( i /= -1 )  dsitransc(ipcgb, i) = ztransp(k)
              ENDDO
           ENDDO
        ENDDO
        DEALLOCATE ( zdirs, vffrac, ztransp )

!         CALL radiation_write_debug_log( 'End of calculation SVF' )
!         WRITE(msg, *) 'Raytracing skipped for maximum distance of ', &
!             max_raytracing_dist, ' m on ', ray_skip_maxdist, ' pairs.'
!         CALL radiation_write_debug_log( msg )
!         WRITE(msg, *) 'Raytracing skipped for minimum potential value of ', &
!             min_irrf_value , ' on ', ray_skip_minval, ' pairs.'
!         CALL radiation_write_debug_log( msg )

        CALL location_message( '    waiting for completion of SVF and CSF calculation in all processes', .TRUE. )
!--     deallocate temporary global arrays
        DEALLOCATE(nzterr)
        
        IF ( plant_canopy )  THEN
!--         finalize mpi_rma communication and deallocate temporary arrays
#if defined( __parallel )
            IF ( rma_lad_raytrace )  THEN
                CALL MPI_Win_flush_all(win_lad, ierr)
!--             unlock MPI window
                CALL MPI_Win_unlock_all(win_lad, ierr)
!--             free MPI window
                CALL MPI_Win_free(win_lad, ierr)
                
!--             deallocate temporary arrays storing values for csf calculation during raytracing
                DEALLOCATE( lad_s_ray )
!--             sub_lad is the pointer to lad_s_rma in case of rma_lad_raytrace
!--             and must not be deallocated here
            ELSE
                DEALLOCATE(sub_lad)
                DEALLOCATE(sub_lad_g)
            ENDIF
#else
            DEALLOCATE(sub_lad)
#endif
            DEALLOCATE( boxes )
            DEALLOCATE( crlens )
            DEALLOCATE( plantt )
            DEALLOCATE( rt2_track, rt2_track_lad, rt2_track_dist, rt2_dist )
        ENDIF

        CALL location_message( '    calculation of the complete SVF array', .TRUE. )

!         CALL radiation_write_debug_log( 'Start SVF sort' )
!--     sort svf ( a version of quicksort )
        CALL quicksort_svf(asvf,1,nsvfl)

        !< load svf from the structure array to plain arrays
!         CALL radiation_write_debug_log( 'Load svf from the structure array to plain arrays' )
        ALLOCATE( svf(ndsvf,nsvfl) )
        ALLOCATE( svfsurf(idsvf,nsvfl) )
        svfnorm_counts(:) = 0._wp
        isurflt_prev = -1
        ksvf = 1
        svfsum = 0._wp
        DO isvf = 1, nsvfl
!--         normalize svf per target face
            IF ( asvf(ksvf)%isurflt /= isurflt_prev )  THEN
                IF ( isurflt_prev /= -1  .AND.  svfsum /= 0._wp )  THEN
                    !< update histogram of logged svf normalization values
                    i = searchsorted(svfnorm_report_thresh, svfsum / (1._wp-skyvf(isurflt_prev)))
                    svfnorm_counts(i) = svfnorm_counts(i) + 1

                    svf(1, isvf_surflt:isvf-1) = svf(1, isvf_surflt:isvf-1) / svfsum * (1._wp-skyvf(isurflt_prev))
                ENDIF
                isurflt_prev = asvf(ksvf)%isurflt
                isvf_surflt = isvf
                svfsum = asvf(ksvf)%rsvf !?? / asvf(ksvf)%rtransp
            ELSE
                svfsum = svfsum + asvf(ksvf)%rsvf !?? / asvf(ksvf)%rtransp
            ENDIF

            svf(:, isvf) = (/ asvf(ksvf)%rsvf, asvf(ksvf)%rtransp /)
            svfsurf(:, isvf) = (/ asvf(ksvf)%isurflt, asvf(ksvf)%isurfs /)

!--         next element
            ksvf = ksvf + 1
        ENDDO

        IF ( isurflt_prev /= -1  .AND.  svfsum /= 0._wp )  THEN
            i = searchsorted(svfnorm_report_thresh, svfsum / (1._wp-skyvf(isurflt_prev)))
            svfnorm_counts(i) = svfnorm_counts(i) + 1

            svf(1, isvf_surflt:nsvfl) = svf(1, isvf_surflt:nsvfl) / svfsum * (1._wp-skyvf(isurflt_prev))
        ENDIF
        !TODO we should be able to deallocate skyvf, from now on we only need skyvft

!--     deallocate temporary asvf array
!--     DEALLOCATE(asvf) - ifort has a problem with deallocation of allocatable target
!--     via pointing pointer - we need to test original targets
        IF ( ALLOCATED(asvf1) )  THEN
            DEALLOCATE(asvf1)
        ENDIF
        IF ( ALLOCATED(asvf2) )  THEN
            DEALLOCATE(asvf2)
        ENDIF

        npcsfl = 0
        IF ( plant_canopy )  THEN

            CALL location_message( '    calculation of the complete CSF array', .TRUE. )
!             CALL radiation_write_debug_log( 'Calculation of the complete CSF array' )
!--         sort and merge csf for the last time, keeping the array size to minimum
            CALL merge_and_grow_csf(-1)
            
!--         aggregate csb among processors
!--         allocate necessary arrays
            ALLOCATE( csflt(ndcsf,max(ncsfl,ndcsf)) )
            ALLOCATE( kcsflt(kdcsf,max(ncsfl,kdcsf)) )
            ALLOCATE( icsflt(0:numprocs-1) )
            ALLOCATE( dcsflt(0:numprocs-1) )
            ALLOCATE( ipcsflt(0:numprocs-1) )
            ALLOCATE( dpcsflt(0:numprocs-1) )
            
!--         fill out arrays of csf values and 
!--         arrays of number of elements and displacements
!--         for particular precessors
            icsflt = 0
            dcsflt = 0
            ip = -1
            j = -1
            d = 0
            DO kcsf = 1, ncsfl
                j = j+1
                IF ( acsf(kcsf)%ip /= ip )  THEN
!--                 new block of the processor
!--                 number of elements of previous block
                    IF ( ip>=0) icsflt(ip) = j
                    d = d+j
!--                 blank blocks
                    DO jp = ip+1, acsf(kcsf)%ip-1
!--                     number of elements is zero, displacement is equal to previous
                        icsflt(jp) = 0
                        dcsflt(jp) = d
                    ENDDO
!--                 the actual block
                    ip = acsf(kcsf)%ip
                    dcsflt(ip) = d
                    j = 0
                ENDIF
!--             fill out real values of rsvf, rtransp
                csflt(1,kcsf) = acsf(kcsf)%rsvf
                csflt(2,kcsf) = acsf(kcsf)%rtransp
!--             fill out integer values of itz,ity,itx,isurfs
                kcsflt(1,kcsf) = acsf(kcsf)%itz
                kcsflt(2,kcsf) = acsf(kcsf)%ity
                kcsflt(3,kcsf) = acsf(kcsf)%itx
                kcsflt(4,kcsf) = acsf(kcsf)%isurfs
            ENDDO
!--         last blank blocks at the end of array
            j = j+1
            IF ( ip>=0 ) icsflt(ip) = j
            d = d+j
            DO jp = ip+1, numprocs-1
!--             number of elements is zero, displacement is equal to previous
                icsflt(jp) = 0
                dcsflt(jp) = d
            ENDDO
            
!--         deallocate temporary acsf array
!--         DEALLOCATE(acsf) - ifort has a problem with deallocation of allocatable target
!--         via pointing pointer - we need to test original targets
            IF ( ALLOCATED(acsf1) )  THEN
                DEALLOCATE(acsf1)
            ENDIF
            IF ( ALLOCATED(acsf2) )  THEN
                DEALLOCATE(acsf2)
            ENDIF
                    
#if defined( __parallel )
!--         scatter and gather the number of elements to and from all processor
!--         and calculate displacements
!             CALL radiation_write_debug_log( 'Scatter and gather the number of elements to and from all processor' )
            CALL MPI_AlltoAll(icsflt,1,MPI_INTEGER,ipcsflt,1,MPI_INTEGER,comm2d, ierr)
            
            npcsfl = SUM(ipcsflt)
            d = 0
            DO i = 0, numprocs-1
                dpcsflt(i) = d
                d = d + ipcsflt(i)
            ENDDO

!--         exchange csf fields between processors
!             CALL radiation_write_debug_log( 'Exchange csf fields between processors' )
            ALLOCATE( pcsflt(ndcsf,max(npcsfl,ndcsf)) )
            ALLOCATE( kpcsflt(kdcsf,max(npcsfl,kdcsf)) )
            CALL MPI_AlltoAllv(csflt, ndcsf*icsflt, ndcsf*dcsflt, MPI_REAL, &
                pcsflt, ndcsf*ipcsflt, ndcsf*dpcsflt, MPI_REAL, comm2d, ierr)
            CALL MPI_AlltoAllv(kcsflt, kdcsf*icsflt, kdcsf*dcsflt, MPI_INTEGER, &
                kpcsflt, kdcsf*ipcsflt, kdcsf*dpcsflt, MPI_INTEGER, comm2d, ierr)
            
#else
            npcsfl = ncsfl
            ALLOCATE( pcsflt(ndcsf,max(npcsfl,ndcsf)) )
            ALLOCATE( kpcsflt(kdcsf,max(npcsfl,kdcsf)) )
            pcsflt = csflt
            kpcsflt = kcsflt
#endif

!--         deallocate temporary arrays
            DEALLOCATE( csflt )
            DEALLOCATE( kcsflt )
            DEALLOCATE( icsflt )
            DEALLOCATE( dcsflt )
            DEALLOCATE( ipcsflt )
            DEALLOCATE( dpcsflt )

!--         sort csf ( a version of quicksort )
!             CALL radiation_write_debug_log( 'Sort csf' )
            CALL quicksort_csf2(kpcsflt, pcsflt, 1, npcsfl)

!--         aggregate canopy sink factor records with identical box & source
!--         againg across all values from all processors
!             CALL radiation_write_debug_log( 'Aggregate canopy sink factor records with identical box' )

            IF ( npcsfl > 0 )  THEN
                icsf = 1 !< reading index
                kcsf = 1 !< writing index
                DO while (icsf < npcsfl)
!--                 here kpcsf(kcsf) already has values from kpcsf(icsf)
                    IF ( kpcsflt(3,icsf) == kpcsflt(3,icsf+1)  .AND.  &
                         kpcsflt(2,icsf) == kpcsflt(2,icsf+1)  .AND.  &
                         kpcsflt(1,icsf) == kpcsflt(1,icsf+1)  .AND.  &
                         kpcsflt(4,icsf) == kpcsflt(4,icsf+1) )  THEN
!--                     We could simply take either first or second rtransp, both are valid. As a very simple heuristic about which ray
!--                     probably passes nearer the center of the target box, we choose DIF from the entry with greater CSF, since that
!--                     might mean that the traced beam passes longer through the canopy box.
                        IF ( pcsflt(1,kcsf) < pcsflt(1,icsf+1) )  THEN
                            pcsflt(2,kcsf) = pcsflt(2,icsf+1)
                        ENDIF
                        pcsflt(1,kcsf) = pcsflt(1,kcsf) + pcsflt(1,icsf+1)

!--                     advance reading index, keep writing index
                        icsf = icsf + 1
                    ELSE
!--                     not identical, just advance and copy
                        icsf = icsf + 1
                        kcsf = kcsf + 1
                        kpcsflt(:,kcsf) = kpcsflt(:,icsf)
                        pcsflt(:,kcsf) = pcsflt(:,icsf)
                    ENDIF
                ENDDO
!--             last written item is now also the last item in valid part of array
                npcsfl = kcsf
            ENDIF

            ncsfl = npcsfl
            IF ( ncsfl > 0 )  THEN
                ALLOCATE( csf(ndcsf,ncsfl) )
                ALLOCATE( csfsurf(idcsf,ncsfl) )
                DO icsf = 1, ncsfl
                    csf(:,icsf) = pcsflt(:,icsf)
                    csfsurf(1,icsf) =  gridpcbl(kpcsflt(1,icsf),kpcsflt(2,icsf),kpcsflt(3,icsf))
                    csfsurf(2,icsf) =  kpcsflt(4,icsf)
                ENDDO
            ENDIF
            
!--         deallocation of temporary arrays
            DEALLOCATE( pcsflt )
            DEALLOCATE( kpcsflt )
!             CALL radiation_write_debug_log( 'End of aggregate csf' )
            
        ENDIF

#if defined( __parallel )
        CALL MPI_BARRIER( comm2d, ierr )
#endif
!         CALL radiation_write_debug_log( 'End of radiation_calc_svf (after mpi_barrier)' )

        RETURN
        
301     WRITE( message_string, * )  &
            'I/O error when processing shape view factors / ',  &
            'plant canopy sink factors / direct irradiance factors.'
        CALL message( 'init_urban_surface', 'PA0502', 2, 2, 0, 6, 0 )
       
    END SUBROUTINE radiation_calc_svf

    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Raytracing for detecting obstacles and calculating compound canopy sink
!> factors. (A simple obstacle detection would only need to process faces in
!> 3 dimensions without any ordering.)
!> Assumtions:
!> -----------
!> 1. The ray always originates from a face midpoint (only one coordinate equals
!>    *.5, i.e. wall) and doesn't travel parallel to the surface (that would mean
!>    shape factor=0). Therefore, the ray may never travel exactly along a face
!>    or an edge.
!> 2. From grid bottom to urban surface top the grid has to be *equidistant*
!>    within each of the dimensions, including vertical (but the resolution
!>    doesn't need to be the same in all three dimensions).
!------------------------------------------------------------------------------!
    SUBROUTINE raytrace(src, targ, isrc, rirrf, atarg, create_csf, visible, transparency, win_lad)
        IMPLICIT NONE

        REAL(wp), DIMENSION(3), INTENT(in)     :: src, targ    !< real coordinates z,y,x
        INTEGER(iwp), INTENT(in)               :: isrc         !< index of source face for csf
        REAL(wp), INTENT(in)                   :: rirrf        !< irradiance factor for csf
        REAL(wp), INTENT(in)                   :: atarg        !< target surface area for csf
        LOGICAL, INTENT(in)                    :: create_csf   !< whether to generate new CSFs during raytracing
        LOGICAL, INTENT(out)                   :: visible
        REAL(wp), INTENT(out)                  :: transparency !< along whole path
        INTEGER(iwp), INTENT(in)               :: win_lad
        INTEGER(iwp)                           :: i, j, k, d
        INTEGER(iwp)                           :: seldim       !< dimension to be incremented
        INTEGER(iwp)                           :: ncsb         !< no of written plant canopy sinkboxes
        INTEGER(iwp)                           :: maxboxes     !< max no of gridboxes visited
        REAL(wp)                               :: distance     !< euclidean along path
        REAL(wp)                               :: crlen        !< length of gridbox crossing
        REAL(wp)                               :: lastdist     !< beginning of current crossing
        REAL(wp)                               :: nextdist     !< end of current crossing
        REAL(wp)                               :: realdist     !< distance in meters per unit distance
        REAL(wp)                               :: crmid        !< midpoint of crossing
        REAL(wp)                               :: cursink      !< sink factor for current canopy box
        REAL(wp), DIMENSION(3)                 :: delta        !< path vector
        REAL(wp), DIMENSION(3)                 :: uvect        !< unit vector
        REAL(wp), DIMENSION(3)                 :: dimnextdist  !< distance for each dimension increments
        INTEGER(iwp), DIMENSION(3)             :: box          !< gridbox being crossed
        INTEGER(iwp), DIMENSION(3)             :: dimnext      !< next dimension increments along path
        INTEGER(iwp), DIMENSION(3)             :: dimdelta     !< dimension direction = +- 1
        INTEGER(iwp)                           :: px, py       !< number of processors in x and y dir before 
                                                               !< the processor in the question
        INTEGER(iwp)                           :: ip           !< number of processor where gridbox reside
        INTEGER(iwp)                           :: ig           !< 1D index of gridbox in global 2D array
        REAL(wp)                               :: lad_s_target !< recieved lad_s of particular grid box
        REAL(wp), PARAMETER                    :: grow_factor = 1.5_wp !< factor of expansion of grow arrays

!
!--     Maximum number of gridboxes visited equals to maximum number of boundaries crossed in each dimension plus one. That's also
!--     the maximum number of plant canopy boxes written. We grow the acsf array accordingly using exponential factor.
        maxboxes = SUM(ABS(NINT(targ, iwp) - NINT(src, iwp))) + 1
        IF ( plant_canopy  .AND.  ncsfl + maxboxes > ncsfla )  THEN
!--         use this code for growing by fixed exponential increments (equivalent to case where ncsfl always increases by 1)
!--         k = CEILING(grow_factor ** real(CEILING(log(real(ncsfl + maxboxes, kind=wp)) &
!--                                                / log(grow_factor)), kind=wp))
!--         or use this code to simply always keep some extra space after growing
            k = CEILING(REAL(ncsfl + maxboxes, kind=wp) * grow_factor)

            CALL merge_and_grow_csf(k)
        ENDIF
        
        transparency = 1._wp
        ncsb = 0

        delta(:) = targ(:) - src(:)
        distance = SQRT(SUM(delta(:)**2))
        IF ( distance == 0._wp )  THEN
            visible = .TRUE.
            RETURN
        ENDIF
        uvect(:) = delta(:) / distance
        realdist = SQRT(SUM( (uvect(:)*(/dz(1),dy,dx/))**2 ))

        lastdist = 0._wp

!--     Since all face coordinates have values *.5 and we'd like to use
!--     integers, all these have .5 added
        DO d = 1, 3
            IF ( uvect(d) == 0._wp )  THEN
                dimnext(d) = 999999999
                dimdelta(d) = 999999999
                dimnextdist(d) = 1.0E20_wp
            ELSE IF ( uvect(d) > 0._wp )  THEN
                dimnext(d) = CEILING(src(d) + .5_wp)
                dimdelta(d) = 1
                dimnextdist(d) = (dimnext(d) - .5_wp - src(d)) / uvect(d)
            ELSE
                dimnext(d) = FLOOR(src(d) + .5_wp)
                dimdelta(d) = -1
                dimnextdist(d) = (dimnext(d) - .5_wp - src(d)) / uvect(d)
            ENDIF
        ENDDO

        DO
!--         along what dimension will the next wall crossing be?
            seldim = minloc(dimnextdist, 1)
            nextdist = dimnextdist(seldim)
            IF ( nextdist > distance ) nextdist = distance

            crlen = nextdist - lastdist
            IF ( crlen > .001_wp )  THEN
                crmid = (lastdist + nextdist) * .5_wp
                box = NINT(src(:) + uvect(:) * crmid, iwp)

!--             calculate index of the grid with global indices (box(2),box(3))
!--             in the array nzterr and plantt and id of the coresponding processor
                px = box(3)/nnx
                py = box(2)/nny
                ip = px*pdims(2)+py
                ig = ip*nnx*nny + (box(3)-px*nnx)*nny + box(2)-py*nny
                IF ( box(1) <= nzterr(ig) )  THEN
                    visible = .FALSE.
                    RETURN
                ENDIF

                IF ( plant_canopy )  THEN
                    IF ( box(1) <= plantt(ig) )  THEN
                        ncsb = ncsb + 1
                        boxes(:,ncsb) = box
                        crlens(ncsb) = crlen
#if defined( __parallel )
                        lad_ip(ncsb) = ip
                        lad_disp(ncsb) = (box(3)-px*nnx)*(nny*nzp) + (box(2)-py*nny)*nzp + box(1)-nzub
#endif
                    ENDIF
                ENDIF
            ENDIF

            IF ( nextdist >= distance ) EXIT
            lastdist = nextdist
            dimnext(seldim) = dimnext(seldim) + dimdelta(seldim)
            dimnextdist(seldim) = (dimnext(seldim) - .5_wp - src(seldim)) / uvect(seldim)
        ENDDO
        
        IF ( plant_canopy )  THEN
#if defined( __parallel )
            IF ( rma_lad_raytrace )  THEN
!--             send requests for lad_s to appropriate processor
                CALL cpu_log( log_point_s(77), 'rad_init_rma', 'start' )
                DO i = 1, ncsb
                    CALL MPI_Get(lad_s_ray(i), 1, MPI_REAL, lad_ip(i), lad_disp(i), &
                                 1, MPI_REAL, win_lad, ierr)
                    IF ( ierr /= 0 )  THEN
                        WRITE(message_string, *) 'MPI error ', ierr, ' at MPI_Get'
                        CALL message( 'raytrace', 'PA0519', 1, 2, 0, 6, 0 )
                    ENDIF
                ENDDO
                
!--             wait for all pending local requests complete
                CALL MPI_Win_flush_local_all(win_lad, ierr)
                IF ( ierr /= 0 )  THEN
                    WRITE(message_string, *) 'MPI error ', ierr, ' at MPI_Win_flush_local_all'
                    CALL message( 'raytrace', 'PA0519', 1, 2, 0, 6, 0 )
                ENDIF
                CALL cpu_log( log_point_s(77), 'rad_init_rma', 'stop' )
                
            ENDIF
#endif

!--         calculate csf and transparency
            DO i = 1, ncsb
#if defined( __parallel )
                IF ( rma_lad_raytrace )  THEN
                    lad_s_target = lad_s_ray(i)
                ELSE
                    lad_s_target = sub_lad_g(lad_ip(i)*nnx*nny*nzp + lad_disp(i))
                ENDIF
#else
                lad_s_target = sub_lad(boxes(1,i),boxes(2,i),boxes(3,i))
#endif
                cursink = 1._wp - exp(-ext_coef * lad_s_target * crlens(i)*realdist)

                IF ( create_csf )  THEN
!--                 write svf values into the array
                    ncsfl = ncsfl + 1
                    acsf(ncsfl)%ip = lad_ip(i)
                    acsf(ncsfl)%itx = boxes(3,i)
                    acsf(ncsfl)%ity = boxes(2,i)
                    acsf(ncsfl)%itz = boxes(1,i)
                    acsf(ncsfl)%isurfs = isrc
                    acsf(ncsfl)%rsvf = REAL(cursink*rirrf*atarg, wp) !-- we postpone multiplication by transparency
                    acsf(ncsfl)%rtransp = REAL(transparency, wp)
                ENDIF  !< create_csf

                transparency = transparency * (1._wp - cursink)
                
            ENDDO
        ENDIF
        
        visible = .TRUE.

    END SUBROUTINE raytrace
    
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> A new, more efficient version of ray tracing algorithm that processes a whole
!> arc instead of a single ray.
!>
!> In all comments, horizon means tangent of horizon angle, i.e.
!> vertical_delta / horizontal_distance
!------------------------------------------------------------------------------!
   SUBROUTINE raytrace_2d(origin, yxdir, zdirs, iorig, aorig, vffrac, &
                              create_csf, skip_1st_pcb, win_lad, horizon, &
                              transparency)
      IMPLICIT NONE

      REAL(wp), DIMENSION(3), INTENT(IN)     ::  origin        !< z,y,x coordinates of ray origin
      REAL(wp), DIMENSION(2), INTENT(IN)     ::  yxdir         !< y,x *unit* vector of ray direction (in grid units)
      REAL(wp), DIMENSION(:), INTENT(IN)     ::  zdirs         !< list of z directions to raytrace (z/hdist, in grid)
      INTEGER(iwp), INTENT(in)               ::  iorig         !< index of origin face for csf
      REAL(wp), INTENT(in)                   ::  aorig         !< origin face area for csf
      REAL(wp), DIMENSION(LBOUND(zdirs, 1):UBOUND(zdirs, 1)), INTENT(in) ::  vffrac !< 
                                                               !< view factor fractions of each ray for csf
      LOGICAL, INTENT(in)                    ::  create_csf    !< whether to generate new CSFs during raytracing
      LOGICAL, INTENT(in)                    ::  skip_1st_pcb  !< whether to skip first plant canopy box during raytracing
      INTEGER(iwp), INTENT(in)               ::  win_lad       !< leaf area density MPI window
      REAL(wp), INTENT(OUT)                  ::  horizon       !< highest horizon found after raytracing (z/hdist)
      REAL(wp), DIMENSION(LBOUND(zdirs, 1):UBOUND(zdirs, 1)), INTENT(OUT) ::  transparency !<
                                                                !< transparencies of zdirs paths
      !--INTEGER(iwp), DIMENSION(3, LBOUND(zdirs, 1):UBOUND(zdirs, 1)), INTENT(OUT) :: itarget !<
                                                                !< (z,y,x) coordinates of target faces for zdirs
      INTEGER(iwp)                           ::  i, k, l, d
      INTEGER(iwp)                           ::  seldim       !< dimension to be incremented
      REAL(wp), DIMENSION(2)                 ::  yxorigin     !< horizontal copy of origin (y,x)
      REAL(wp)                               ::  distance     !< euclidean along path
      REAL(wp)                               ::  lastdist     !< beginning of current crossing
      REAL(wp)                               ::  nextdist     !< end of current crossing
      REAL(wp)                               ::  crmid        !< midpoint of crossing
      REAL(wp)                               ::  horz_entry   !< horizon at entry to column
      REAL(wp)                               ::  horz_exit    !< horizon at exit from column
      REAL(wp)                               ::  bdydim       !< boundary for current dimension
      REAL(wp), DIMENSION(2)                 ::  crossdist    !< distances to boundary for dimensions
      REAL(wp), DIMENSION(2)                 ::  dimnextdist  !< distance for each dimension increments
      INTEGER(iwp), DIMENSION(2)             ::  column       !< grid column being crossed
      INTEGER(iwp), DIMENSION(2)             ::  dimnext      !< next dimension increments along path
      INTEGER(iwp), DIMENSION(2)             ::  dimdelta     !< dimension direction = +- 1
      INTEGER(iwp)                           ::  px, py       !< number of processors in x and y dir before 
                                                              !< the processor in the question
      INTEGER(iwp)                           ::  ip           !< number of processor where gridbox reside
      INTEGER(iwp)                           ::  ig           !< 1D index of gridbox in global 2D array
      INTEGER(iwp)                           ::  wcount       !< RMA window item count
      INTEGER(iwp)                           ::  maxboxes     !< max no of CSF created
      INTEGER(iwp)                           ::  nly          !< maximum  plant canopy height 
      INTEGER(iwp)                           ::  ntrack
      REAL(wp)                               ::  zbottom, ztop !< urban surface boundary in real numbers
      REAL(wp)                               ::  zorig        !< z coordinate of ray column entry
      REAL(wp)                               ::  zexit        !< z coordinate of ray column exit
      REAL(wp)                               ::  qdist        !< ratio of real distance to z coord difference
      REAL(wp)                               ::  dxxyy        !< square of real horizontal distance
      REAL(wp)                               ::  curtrans     !< transparency of current PC box crossing
      INTEGER(iwp)                           ::  zb0
      INTEGER(iwp)                           ::  zb1
      INTEGER(iwp)                           ::  nz
      INTEGER(iwp)                           ::  iz
      INTEGER(iwp)                           ::  zsgn
      REAL(wp), PARAMETER                    ::  grow_factor = 1.5_wp !< factor of expansion of grow arrays

#if defined( __parallel )
      INTEGER(MPI_ADDRESS_KIND)              ::  wdisp        !< RMA window displacement
#endif
      
      yxorigin(:) = origin(2:3)
      transparency(:) = 1._wp !-- Pre-set the all rays to transparent before reducing
      horizon = -HUGE(1._wp)

      !--Determine distance to boundary (in 2D xy)
      IF ( yxdir(1) > 0._wp )  THEN
         bdydim = ny + .5_wp !< north global boundary
         crossdist(1) = (bdydim - yxorigin(1)) / yxdir(1)
      ELSEIF ( yxdir(1) == 0._wp )  THEN
         crossdist(1) = HUGE(1._wp)
      ELSE
          bdydim = -.5_wp !< south global boundary
          crossdist(1) = (bdydim - yxorigin(1)) / yxdir(1)
      ENDIF

      IF ( yxdir(2) >= 0._wp )  THEN
          bdydim = nx + .5_wp !< east global boundary
          crossdist(2) = (bdydim - yxorigin(2)) / yxdir(2)
      ELSEIF ( yxdir(2) == 0._wp )  THEN
         crossdist(2) = HUGE(1._wp)
      ELSE
          bdydim = -.5_wp !< west global boundary
          crossdist(2) = (bdydim - yxorigin(2)) / yxdir(2)
      ENDIF
      distance = minval(crossdist, 1)

      IF ( plant_canopy )  THEN
         rt2_track_dist(0) = 0._wp
         rt2_track_lad(:,:) = 0._wp
         nly = plantt_max - nzub + 1
      ENDIF

      lastdist = 0._wp

!--   Since all face coordinates have values *.5 and we'd like to use
!--   integers, all these have .5 added
      DO d = 1, 2
          IF ( yxdir(d) == 0._wp )  THEN
              dimnext(d) = HUGE(1_iwp)
              dimdelta(d) = HUGE(1_iwp)
              dimnextdist(d) = HUGE(1._wp)
          ELSE IF ( yxdir(d) > 0._wp )  THEN
              dimnext(d) = FLOOR(yxorigin(d) + .5_wp) + 1
              dimdelta(d) = 1
              dimnextdist(d) = (dimnext(d) - .5_wp - yxorigin(d)) / yxdir(d)
          ELSE
              dimnext(d) = CEILING(yxorigin(d) + .5_wp) - 1
              dimdelta(d) = -1
              dimnextdist(d) = (dimnext(d) - .5_wp - yxorigin(d)) / yxdir(d)
          ENDIF
      ENDDO

      ntrack = 0
      DO
!--      along what dimension will the next wall crossing be?
         seldim = minloc(dimnextdist, 1)
         nextdist = dimnextdist(seldim)
         IF ( nextdist > distance ) nextdist = distance

         IF ( nextdist > lastdist )  THEN
            ntrack = ntrack + 1
            crmid = (lastdist + nextdist) * .5_wp
            column = NINT(yxorigin(:) + yxdir(:) * crmid, iwp)

!--         calculate index of the grid with global indices (column(1),column(2))
!--         in the array nzterr and plantt and id of the coresponding processor
            px = column(2)/nnx
            py = column(1)/nny
            ip = px*pdims(2)+py
            ig = ip*nnx*nny + (column(2)-px*nnx)*nny + column(1)-py*nny

            IF ( lastdist == 0._wp )  THEN
               horz_entry = -HUGE(1._wp)
            ELSE
               horz_entry = (nzterr(ig) - origin(1)) / lastdist
            ENDIF
            horz_exit = (nzterr(ig) - origin(1)) / nextdist
            horizon = MAX(horizon, horz_entry, horz_exit)

            IF ( plant_canopy )  THEN
               rt2_track(:, ntrack) = column(:)
               rt2_track_dist(ntrack) = nextdist
            ENDIF
         ENDIF

         IF ( nextdist >= distance )  EXIT
         lastdist = nextdist
         dimnext(seldim) = dimnext(seldim) + dimdelta(seldim)
         dimnextdist(seldim) = (dimnext(seldim) - .5_wp - yxorigin(seldim)) / yxdir(seldim)
      ENDDO

      IF ( plant_canopy )  THEN
         !--Request LAD WHERE applicable
         !--
#if defined( __parallel )
         IF ( rma_lad_raytrace )  THEN
!--         send requests for lad_s to appropriate processor
            !CALL cpu_log( log_point_s(77), 'usm_init_rma', 'start' )
            DO i = 1, ntrack
               px = rt2_track(2,i)/nnx
               py = rt2_track(1,i)/nny
               ip = px*pdims(2)+py
               ig = ip*nnx*nny + (rt2_track(2,i)-px*nnx)*nny + rt2_track(1,i)-py*nny
               IF ( plantt(ig) <= nzterr(ig) )  CYCLE
               wdisp = (rt2_track(2,i)-px*nnx)*(nny*nzp) + (rt2_track(1,i)-py*nny)*nzp + nzterr(ig)+1-nzub
               wcount = plantt(ig)-nzterr(ig)
               ! TODO send request ASAP - even during raytracing
               CALL MPI_Get(rt2_track_lad(nzterr(ig)+1:plantt(ig), i), wcount, MPI_REAL, ip,    &
                            wdisp, wcount, MPI_REAL, win_lad, ierr)
               IF ( ierr /= 0 )  THEN
                  WRITE(message_string, *) 'MPI error ', ierr, ' at MPI_Get'
                  CALL message( 'raytrace_2d', 'PA0526', 1, 2, 0, 6, 0 )
               ENDIF
            ENDDO

!--         wait for all pending local requests complete
            ! TODO WAIT selectively for each column later when needed
            CALL MPI_Win_flush_local_all(win_lad, ierr)
            IF ( ierr /= 0 )  THEN
               WRITE(message_string, *) 'MPI error ', ierr, ' at MPI_Win_flush_local_all'
               CALL message( 'raytrace', 'PA0527', 1, 2, 0, 6, 0 )
            ENDIF
            !CALL cpu_log( log_point_s(77), 'usm_init_rma', 'stop' )
         ELSE ! rma_lad_raytrace
            DO i = 1, ntrack
               px = rt2_track(2,i)/nnx
               py = rt2_track(1,i)/nny
               ip = px*pdims(2)+py
               ig = ip*nnx*nny*nzp + (rt2_track(2,i)-px*nnx)*(nny*nzp) + (rt2_track(1,i)-py*nny)*nzp
               rt2_track_lad(nzub:plantt_max, i) = sub_lad_g(ig:ig+nly-1)
            ENDDO
         ENDIF
#else
         DO i = 1, ntrack
            rt2_track_lad(nzub:plantt_max, i) = sub_lad(rt2_track(1,i), rt2_track(2,i), nzub:plantt_max)
         ENDDO
#endif

         !--Skip the PCB around origin if requested
         !--
         IF ( skip_1st_pcb )  THEN
            rt2_track_lad(NINT(origin(1), iwp), 1) = 0._wp
         ENDIF

         !--Assert that we have space allocated for CSFs
         !--
         maxboxes = (ntrack + MAX(origin(1) - nzub, nzpt - origin(1))) * SIZE(zdirs, 1)
         IF ( ncsfl + maxboxes > ncsfla )  THEN
!--         use this code for growing by fixed exponential increments (equivalent to case where ncsfl always increases by 1)
!--         k = CEILING(grow_factor ** real(CEILING(log(real(ncsfl + maxboxes, kind=wp)) &
!--                                                / log(grow_factor)), kind=wp))
!--         or use this code to simply always keep some extra space after growing
            k = CEILING(REAL(ncsfl + maxboxes, kind=wp) * grow_factor)
            CALL merge_and_grow_csf(k)
         ENDIF

         !--Calculate transparencies and store new CSFs
         !--
         zbottom = REAL(nzub, wp) - .5_wp
         ztop = REAL(plantt_max, wp) + .5_wp

         !--Reverse direction of radiation (face->sky), only when create_csf
         !--
         IF ( create_csf )  THEN
            DO i = 1, ntrack ! for each column
               dxxyy = ((dy*yxdir(1))**2 + (dx*yxdir(2))**2) * (rt2_track_dist(i)-rt2_track_dist(i-1))**2
               px = rt2_track(2,i)/nnx
               py = rt2_track(1,i)/nny
               ip = px*pdims(2)+py

               DO k = LBOUND(zdirs, 1), UBOUND(zdirs, 1) ! for each ray
                  IF ( zdirs(k) <= horizon )  THEN
                     CYCLE
                  ENDIF

                  zorig = REAL(origin(1), wp) + zdirs(k) * rt2_track_dist(i-1)
                  IF ( zorig <= zbottom .OR. zorig >= ztop )  CYCLE

                  zsgn = INT(SIGN(1._wp, zdirs(k)), iwp)
                  rt2_dist(1) = 0._wp
                  IF ( zdirs(k) == 0._wp )  THEN ! ray is exactly horizontal
                     nz = 2
                     rt2_dist(nz) = SQRT(dxxyy)
                     iz = NINT(zorig, iwp)
                  ELSE
                     zexit = MIN(MAX(REAL(origin(1), wp) + zdirs(k) * rt2_track_dist(i), zbottom), ztop)

                     zb0 = FLOOR(  zorig * zsgn - .5_wp) + 1  ! because it must be greater than orig
                     zb1 = CEILING(zexit * zsgn - .5_wp) - 1  ! because it must be smaller than exit
                     nz = MAX(zb1 - zb0 + 3, 2)
                     rt2_dist(nz) = SQRT(((zexit-zorig)*dz(1))**2 + dxxyy)
                     qdist = rt2_dist(nz) / (zexit-zorig)
                     rt2_dist(2:nz-1) = (/( ((REAL(l, wp) + .5_wp) * zsgn - zorig) * qdist , l = zb0, zb1 )/)
                     iz = zb0 * zsgn
                  ENDIF

                  DO l = 2, nz
                     IF ( rt2_track_lad(iz, i) > 0._wp )  THEN
                        curtrans = exp(-ext_coef * rt2_track_lad(iz, i) * (rt2_dist(l)-rt2_dist(l-1)))

                        ncsfl = ncsfl + 1
                        acsf(ncsfl)%ip = ip
                        acsf(ncsfl)%itx = rt2_track(2,i)
                        acsf(ncsfl)%ity = rt2_track(1,i)
                        acsf(ncsfl)%itz = iz
                        acsf(ncsfl)%isurfs = iorig
                        acsf(ncsfl)%rsvf = REAL((1._wp - curtrans)*aorig*vffrac(k), wp) ! we postpone multiplication by transparency
                        acsf(ncsfl)%rtransp = REAL(transparency(k), wp)

                        transparency(k) = transparency(k) * curtrans
                     ENDIF
                     iz = iz + zsgn
                  ENDDO ! l = 1, nz - 1
               ENDDO ! k = LBOUND(zdirs, 1), UBOUND(zdirs, 1)
            ENDDO ! i = 1, ntrack

            transparency(:) = 1._wp !-- Reset all rays to transparent
         ENDIF

         !-- Forward direction of radiation (sky->face), always
         !--
         DO i = ntrack, 1, -1 ! for each column backwards
            dxxyy = ((dy*yxdir(1))**2 + (dx*yxdir(2))**2) * (rt2_track_dist(i)-rt2_track_dist(i-1))**2
            px = rt2_track(2,i)/nnx
            py = rt2_track(1,i)/nny
            ip = px*pdims(2)+py

            DO k = LBOUND(zdirs, 1), UBOUND(zdirs, 1) ! for each ray
               IF ( zdirs(k) <= horizon )  THEN
                  transparency(k) = 0._wp
                  CYCLE
               ENDIF

               zexit = REAL(origin(1), wp) + zdirs(k) * rt2_track_dist(i-1)
               IF ( zexit <= zbottom .OR. zexit >= ztop )  CYCLE

               zsgn = -INT(SIGN(1._wp, zdirs(k)), iwp)
               rt2_dist(1) = 0._wp
               IF ( zdirs(k) == 0._wp )  THEN ! ray is exactly horizontal
                  nz = 2
                  rt2_dist(nz) = SQRT(dxxyy)
                  iz = NINT(zexit, iwp)
               ELSE
                  zorig = MIN(MAX(REAL(origin(1), wp) + zdirs(k) * rt2_track_dist(i), zbottom), ztop)

                  zb0 = FLOOR(  zorig * zsgn - .5_wp) + 1  ! because it must be greater than orig
                  zb1 = CEILING(zexit * zsgn - .5_wp) - 1  ! because it must be smaller than exit
                  nz = MAX(zb1 - zb0 + 3, 2)
                  rt2_dist(nz) = SQRT(((zexit-zorig)*dz(1))**2 + dxxyy)
                  qdist = rt2_dist(nz) / (zexit-zorig)
                  rt2_dist(2:nz-1) = (/( ((REAL(l, wp) + .5_wp) * zsgn - zorig) * qdist , l = zb0, zb1 )/)
                  iz = zb0 * zsgn
               ENDIF

               DO l = 2, nz
                  IF ( rt2_track_lad(iz, i) > 0._wp )  THEN
                     curtrans = exp(-ext_coef * rt2_track_lad(iz, i) * (rt2_dist(l)-rt2_dist(l-1)))

                     IF ( create_csf )  THEN
                        ncsfl = ncsfl + 1
                        acsf(ncsfl)%ip = ip
                        acsf(ncsfl)%itx = rt2_track(2,i)
                        acsf(ncsfl)%ity = rt2_track(1,i)
                        acsf(ncsfl)%itz = iz
                        acsf(ncsfl)%isurfs = -1 ! a special ID indicating sky
                        acsf(ncsfl)%rsvf = REAL((1._wp - curtrans)*aorig*vffrac(k), wp) ! we postpone multiplication by transparency
                        acsf(ncsfl)%rtransp = REAL(transparency(k), wp)
                     ENDIF  !< create_csf

                     transparency(k) = transparency(k) * curtrans
                  ENDIF
                  iz = iz + zsgn
               ENDDO ! l = 1, nz - 1
            ENDDO ! k = LBOUND(zdirs, 1), UBOUND(zdirs, 1)
         ENDDO ! i = 1, ntrack

      ELSE ! not plant_canopy
         DO k = UBOUND(zdirs, 1), LBOUND(zdirs, 1), -1 ! TODO make more generic
            IF ( zdirs(k) > horizon )  EXIT
            transparency(k) = 0._wp
         ENDDO
      ENDIF

   END SUBROUTINE raytrace_2d
 

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Calculates apparent solar positions for all timesteps and stores discretized
!> positions.
!------------------------------------------------------------------------------!
   SUBROUTINE radiation_presimulate_solar_pos
      IMPLICIT NONE

      INTEGER(iwp)                              ::  it, i, j
      REAL(wp)                                  ::  tsrp_prev
      REAL(wp), DIMENSION(:,:), ALLOCATABLE     ::  dsidir_tmp       !< dsidir_tmp[:,i] = unit vector of i-th
                                                                     !< appreant solar direction

      ALLOCATE ( dsidir_rev(0:raytrace_discrete_elevs/2-1,                 &
                            0:raytrace_discrete_azims-1) )
      dsidir_rev(:,:) = -1
      ALLOCATE ( dsidir_tmp(3,                                             &
                     raytrace_discrete_elevs/2*raytrace_discrete_azims) )
      ndsidir = 0

!
!--   We will artificialy update time_since_reference_point and return to
!--   true value later
      tsrp_prev = time_since_reference_point
      sun_direction = .TRUE.

!
!--   Process spinup time if configured
      IF ( spinup_time > 0._wp )  THEN
         DO  it = 0, CEILING(spinup_time / dt_spinup)
            time_since_reference_point = -spinup_time + REAL(it, wp) * dt_spinup
            CALL simulate_pos
         ENDDO
      ENDIF
!
!--   Process simulation time
      DO  it = 0, CEILING(( end_time - spinup_time ) / dt_radiation)
         time_since_reference_point = REAL(it, wp) * dt_radiation
         CALL simulate_pos
      ENDDO

      time_since_reference_point = tsrp_prev

!--   Allocate global vars which depend on ndsidir
      ALLOCATE ( dsidir ( 3, ndsidir ) )
      dsidir(:,:) = dsidir_tmp(:, 1:ndsidir)
      DEALLOCATE ( dsidir_tmp )

      ALLOCATE ( dsitrans(nsurfl, ndsidir) )
      ALLOCATE ( dsitransc(npcbl, ndsidir) )

      WRITE ( message_string, * ) 'Precalculated', ndsidir, ' solar positions', &
                                  'from', it, ' timesteps.'
      CALL message( 'radiation_presimulate_solar_pos', 'UI0013', 0, 0, 0, 6, 0 )

      CONTAINS

      !------------------------------------------------------------------------!
      ! Description:
      ! ------------
      !> Simuates a single position
      !------------------------------------------------------------------------!
      SUBROUTINE simulate_pos
         IMPLICIT NONE
!
!--      Update apparent solar position based on modified t_s_r_p
         CALL calc_zenith
         IF ( zenith(0) > 0 )  THEN
!--         
!--         Identify solar direction vector (discretized number) 1)
            i = MODULO(NINT(ATAN2(sun_dir_lon(0), sun_dir_lat(0))               &
                            / (2._wp*pi) * raytrace_discrete_azims-.5_wp, iwp), &
                       raytrace_discrete_azims)
            j = FLOOR(ACOS(zenith(0)) / pi * raytrace_discrete_elevs)
            IF ( dsidir_rev(j, i) == -1 )  THEN
               ndsidir = ndsidir + 1
               dsidir_tmp(:, ndsidir) =                                              &
                     (/ COS((REAL(j,wp)+.5_wp) * pi      / raytrace_discrete_elevs), &
                        SIN((REAL(j,wp)+.5_wp) * pi      / raytrace_discrete_elevs)  &
                      * COS((REAL(i,wp)+.5_wp) * 2_wp*pi / raytrace_discrete_azims), &
                        SIN((REAL(j,wp)+.5_wp) * pi      / raytrace_discrete_elevs)  &
                      * SIN((REAL(i,wp)+.5_wp) * 2_wp*pi / raytrace_discrete_azims) /)
               dsidir_rev(j, i) = ndsidir
            ENDIF
         ENDIF
      END SUBROUTINE simulate_pos

   END SUBROUTINE radiation_presimulate_solar_pos



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Determines whether two faces are oriented towards each other. Since the
!> surfaces follow the gird box surfaces, it checks first whether the two surfaces
!> are directed in the same direction, then it checks if the two surfaces are
!> located in confronted direction but facing away from each other, e.g. <--| |-->
!------------------------------------------------------------------------------!
    PURE LOGICAL FUNCTION surface_facing(x, y, z, d, x2, y2, z2, d2)
        IMPLICIT NONE
        INTEGER(iwp),   INTENT(in)  :: x, y, z, d, x2, y2, z2, d2
      
        surface_facing = .FALSE.

!-- first check: are the two surfaces directed in the same direction
        IF ( (d==iup_u  .OR.  d==iup_l  .OR.  d==iup_a )                             &
             .AND. (d2==iup_u  .OR. d2==iup_l) ) RETURN
        IF ( (d==isouth_u  .OR.  d==isouth_l  .OR.  d==isouth_a ) &
             .AND.  (d2==isouth_u  .OR.  d2==isouth_l) ) RETURN
        IF ( (d==inorth_u  .OR.  d==inorth_l  .OR.  d==inorth_a ) &
             .AND.  (d2==inorth_u  .OR.  d2==inorth_l) ) RETURN
        IF ( (d==iwest_u  .OR.  d==iwest_l  .OR.  d==iwest_a )     &
             .AND.  (d2==iwest_u  .OR.  d2==iwest_l ) ) RETURN
        IF ( (d==ieast_u  .OR.  d==ieast_l  .OR.  d==ieast_a )     &
             .AND.  (d2==ieast_u  .OR.  d2==ieast_l ) ) RETURN

!-- second check: are surfaces facing away from each other
        SELECT CASE (d)
            CASE (iup_u, iup_l, iup_a)              !< upward facing surfaces
                IF ( z2 < z ) RETURN
            CASE (idown_a)                          !< downward facing surfaces
                IF ( z2 > z ) RETURN
            CASE (isouth_u, isouth_l, isouth_a)     !< southward facing surfaces
                IF ( y2 > y ) RETURN
            CASE (inorth_u, inorth_l, inorth_a)     !< northward facing surfaces
                IF ( y2 < y ) RETURN
            CASE (iwest_u, iwest_l, iwest_a)        !< westward facing surfaces
                IF ( x2 > x ) RETURN
            CASE (ieast_u, ieast_l, ieast_a)        !< eastward facing surfaces
                IF ( x2 < x ) RETURN
        END SELECT

        SELECT CASE (d2)
            CASE (iup_u)                            !< ground, roof
                IF ( z < z2 ) RETURN
            CASE (isouth_u, isouth_l)               !< south facing
                IF ( y > y2 ) RETURN
            CASE (inorth_u, inorth_l)               !< north facing
                IF ( y < y2 ) RETURN
            CASE (iwest_u, iwest_l)                 !< west facing
                IF ( x > x2 ) RETURN
            CASE (ieast_u, ieast_l)                 !< east facing
                IF ( x < x2 ) RETURN
            CASE (-1)
                CONTINUE
        END SELECT

        surface_facing = .TRUE.
        
    END FUNCTION surface_facing


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Soubroutine reads svf and svfsurf data from saved file
!> SVF means sky view factors and CSF means canopy sink factors
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_read_svf

       IMPLICIT NONE
       
       CHARACTER(rad_version_len)   :: rad_version_field
       
       INTEGER(iwp)                 :: i
       INTEGER(iwp)                 :: ndsidir_from_file = 0
       INTEGER(iwp)                 :: npcbl_from_file = 0
       INTEGER(iwp)                 :: nsurfl_from_file = 0
       
       DO  i = 0, io_blocks-1
          IF ( i == io_group )  THEN

!
!--          numprocs_previous_run is only known in case of reading restart
!--          data. If a new initial run which reads svf data is started the
!--          following query will be skipped
             IF ( initializing_actions == 'read_restart_data' ) THEN

                IF ( numprocs_previous_run /= numprocs ) THEN
                   WRITE( message_string, * ) 'A different number of ',        &
                                              'processors between the run ',   &
                                              'that has written the svf data ',&
                                              'and the one that will read it ',&
                                              'is not allowed' 
                   CALL message( 'check_open', 'PA0491', 1, 2, 0, 6, 0 )
                ENDIF

             ENDIF
              
!
!--          Open binary file 
             CALL check_open( 88 )

!
!--          read and check version
             READ ( 88 ) rad_version_field
             IF ( TRIM(rad_version_field) /= TRIM(rad_version) )  THEN
                 WRITE( message_string, * ) 'Version of binary SVF file "',    &
                             TRIM(rad_version_field), '" does not match ',     &
                             'the version of model "', TRIM(rad_version), '"'
                 CALL message( 'radiation_read_svf', 'PA0482', 1, 2, 0, 6, 0 )
             ENDIF
              
!
!--          read nsvfl, ncsfl, nsurfl
             READ ( 88 ) nsvfl, ncsfl, nsurfl_from_file, npcbl_from_file,      &
                         ndsidir_from_file
             
             IF ( nsvfl < 0  .OR.  ncsfl < 0 )  THEN
                 WRITE( message_string, * ) 'Wrong number of SVF or CSF'
                 CALL message( 'radiation_read_svf', 'PA0483', 1, 2, 0, 6, 0 )
             ELSE
                 WRITE(message_string,*) '    Number of SVF, CSF, and nsurfl ',&
                                         'to read', nsvfl, ncsfl,              &
                                         nsurfl_from_file
                 CALL location_message( message_string, .TRUE. )
             ENDIF
             
             IF ( nsurfl_from_file /= nsurfl )  THEN
                 WRITE( message_string, * ) 'nsurfl from SVF file does not ',  &
                                            'match calculated nsurfl from ',   &
                                            'radiation_interaction_init'
                 CALL message( 'radiation_read_svf', 'PA0490', 1, 2, 0, 6, 0 )
             ENDIF
             
             IF ( npcbl_from_file /= npcbl )  THEN
                 WRITE( message_string, * ) 'npcbl from SVF file does not ',   &
                                            'match calculated npcbl from ',    &
                                            'radiation_interaction_init'
                 CALL message( 'radiation_read_svf', 'PA0493', 1, 2, 0, 6, 0 )
             ENDIF
             
             IF ( ndsidir_from_file /= ndsidir )  THEN
                 WRITE( message_string, * ) 'ndsidir from SVF file does not ', &
                                            'match calculated ndsidir from ',  &
                                            'radiation_presimulate_solar_pos'
                 CALL message( 'radiation_read_svf', 'PA0494', 1, 2, 0, 6, 0 )
             ENDIF
              
!
!--          Arrays skyvf, skyvft, dsitrans and dsitransc are allready 
!--          allocated in radiation_interaction_init and 
!--          radiation_presimulate_solar_pos
             IF ( nsurfl > 0 )  THEN
                READ(88) skyvf
                READ(88) skyvft
                READ(88) dsitrans  
             ENDIF
              
             IF ( plant_canopy  .AND.  npcbl > 0 ) THEN
                READ ( 88 )  dsitransc
             ENDIF
             
!
!--          The allocation of svf, svfsurf, csf and csfsurf happens in routine 
!--          radiation_calc_svf which is not called if the program enters 
!--          radiation_read_svf. Therefore these arrays has to allocate in the 
!--          following
             IF ( nsvfl > 0 )  THEN
                ALLOCATE( svf(ndsvf,nsvfl) )
                ALLOCATE( svfsurf(idsvf,nsvfl) )
                READ(88) svf
                READ(88) svfsurf
             ENDIF

             IF ( plant_canopy  .AND.  ncsfl > 0 )  THEN
                ALLOCATE( csf(ndcsf,ncsfl) )
                ALLOCATE( csfsurf(idcsf,ncsfl) )
                READ(88) csf
                READ(88) csfsurf
             ENDIF
              
!
!--          Close binary file                 
             CALL close_file( 88 )
                
          ENDIF
#if defined( __parallel )
          CALL MPI_BARRIER( comm2d, ierr )
#endif
       ENDDO

    END SUBROUTINE radiation_read_svf


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine stores svf, svfsurf, csf and csfsurf data to a file.
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_write_svf

       IMPLICIT NONE
       
       INTEGER(iwp)        :: i

       DO  i = 0, io_blocks-1
          IF ( i == io_group )  THEN
!
!--          Open binary file 
             CALL check_open( 89 )

             WRITE ( 89 )  rad_version
             WRITE ( 89 )  nsvfl, ncsfl, nsurfl, npcbl, ndsidir
             IF ( nsurfl > 0 ) THEN
                WRITE ( 89 )  skyvf
                WRITE ( 89 )  skyvft
                WRITE ( 89 )  dsitrans
             ENDIF
             IF ( npcbl > 0 ) THEN
                WRITE ( 89 )  dsitransc
             ENDIF
             IF ( nsvfl > 0 ) THEN
                WRITE ( 89 )  svf
                WRITE ( 89 )  svfsurf
             ENDIF
             IF ( plant_canopy  .AND.  ncsfl > 0 )  THEN
                 WRITE ( 89 )  csf
                 WRITE ( 89 )  csfsurf
             ENDIF

!
!--          Close binary file                 
             CALL close_file( 89 )

          ENDIF
#if defined( __parallel )
          CALL MPI_BARRIER( comm2d, ierr )
#endif
       ENDDO
    END SUBROUTINE radiation_write_svf

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Block of auxiliary subroutines:
!> 1. quicksort and corresponding comparison
!> 2. merge_and_grow_csf for implementation of "dynamical growing"
!>    array for csf
!------------------------------------------------------------------------------!
    PURE FUNCTION svf_lt(svf1,svf2) result (res)
      TYPE (t_svf), INTENT(in) :: svf1,svf2
      LOGICAL                  :: res
      IF ( svf1%isurflt < svf2%isurflt  .OR.    &
          (svf1%isurflt == svf2%isurflt  .AND.  svf1%isurfs < svf2%isurfs) )  THEN
          res = .TRUE.
      ELSE
          res = .FALSE.
      ENDIF
    END FUNCTION svf_lt
    
 
!-- quicksort.f -*-f90-*-
!-- Author: t-nissie, adaptation J.Resler
!-- License: GPLv3
!-- Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    RECURSIVE SUBROUTINE quicksort_svf(svfl, first, last)
        IMPLICIT NONE
        TYPE(t_svf), DIMENSION(:), INTENT(INOUT)  :: svfl
        INTEGER(iwp), INTENT(IN)                  :: first, last
        TYPE(t_svf)                               :: x, t
        INTEGER(iwp)                              :: i, j

        IF ( first>=last ) RETURN
        x = svfl( (first+last) / 2 )
        i = first
        j = last
        DO
            DO while ( svf_lt(svfl(i),x) )
               i=i+1
            ENDDO
            DO while ( svf_lt(x,svfl(j)) )
                j=j-1
            ENDDO
            IF ( i >= j ) EXIT
            t = svfl(i);  svfl(i) = svfl(j);  svfl(j) = t
            i=i+1
            j=j-1
        ENDDO
        IF ( first < i-1 ) CALL quicksort_svf(svfl, first, i-1)
        IF ( j+1 < last )  CALL quicksort_svf(svfl, j+1, last)
    END SUBROUTINE quicksort_svf

    
    PURE FUNCTION csf_lt(csf1,csf2) result (res)
      TYPE (t_csf), INTENT(in) :: csf1,csf2
      LOGICAL                  :: res
      IF ( csf1%ip < csf2%ip  .OR.    &
           (csf1%ip == csf2%ip  .AND.  csf1%itx < csf2%itx)  .OR.  &
           (csf1%ip == csf2%ip  .AND.  csf1%itx == csf2%itx  .AND.  csf1%ity < csf2%ity)  .OR.  &
           (csf1%ip == csf2%ip  .AND.  csf1%itx == csf2%itx  .AND.  csf1%ity == csf2%ity  .AND.   &
            csf1%itz < csf2%itz)  .OR.  &
           (csf1%ip == csf2%ip  .AND.  csf1%itx == csf2%itx  .AND.  csf1%ity == csf2%ity  .AND.   &
            csf1%itz == csf2%itz  .AND.  csf1%isurfs < csf2%isurfs) )  THEN
          res = .TRUE.
      ELSE
          res = .FALSE.
      ENDIF
    END FUNCTION csf_lt


!-- quicksort.f -*-f90-*-
!-- Author: t-nissie, adaptation J.Resler
!-- License: GPLv3
!-- Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    RECURSIVE SUBROUTINE quicksort_csf(csfl, first, last)
        IMPLICIT NONE
        TYPE(t_csf), DIMENSION(:), INTENT(INOUT)  :: csfl
        INTEGER(iwp), INTENT(IN)                  :: first, last
        TYPE(t_csf)                               :: x, t
        INTEGER(iwp)                              :: i, j

        IF ( first>=last ) RETURN
        x = csfl( (first+last)/2 )
        i = first
        j = last
        DO
            DO while ( csf_lt(csfl(i),x) )
                i=i+1
            ENDDO
            DO while ( csf_lt(x,csfl(j)) )
                j=j-1
            ENDDO
            IF ( i >= j ) EXIT
            t = csfl(i);  csfl(i) = csfl(j);  csfl(j) = t
            i=i+1
            j=j-1
        ENDDO
        IF ( first < i-1 ) CALL quicksort_csf(csfl, first, i-1)
        IF ( j+1 < last )  CALL quicksort_csf(csfl, j+1, last)
    END SUBROUTINE quicksort_csf

    
    SUBROUTINE merge_and_grow_csf(newsize)
        INTEGER(iwp), INTENT(in)                :: newsize  !< new array size after grow, must be >= ncsfl
                                                            !< or -1 to shrink to minimum
        INTEGER(iwp)                            :: iread, iwrite
        TYPE(t_csf), DIMENSION(:), POINTER      :: acsfnew
        CHARACTER(100)                          :: msg

        IF ( newsize == -1 )  THEN
!--         merge in-place
            acsfnew => acsf
        ELSE
!--         allocate new array
            IF ( mcsf == 0 )  THEN
                ALLOCATE( acsf1(newsize) )
                acsfnew => acsf1
            ELSE
                ALLOCATE( acsf2(newsize) )
                acsfnew => acsf2
            ENDIF
        ENDIF

        IF ( ncsfl >= 1 )  THEN
!--         sort csf in place (quicksort)
            CALL quicksort_csf(acsf,1,ncsfl)

!--         while moving to a new array, aggregate canopy sink factor records with identical box & source
            acsfnew(1) = acsf(1)
            iwrite = 1
            DO iread = 2, ncsfl
!--             here acsf(kcsf) already has values from acsf(icsf)
                IF ( acsfnew(iwrite)%itx == acsf(iread)%itx &
                         .AND.  acsfnew(iwrite)%ity == acsf(iread)%ity &
                         .AND.  acsfnew(iwrite)%itz == acsf(iread)%itz &
                         .AND.  acsfnew(iwrite)%isurfs == acsf(iread)%isurfs )  THEN
!--                 We could simply take either first or second rtransp, both are valid. As a very simple heuristic about which ray
!--                 probably passes nearer the center of the target box, we choose DIF from the entry with greater CSF, since that
!--                 might mean that the traced beam passes longer through the canopy box.
                    IF ( acsfnew(iwrite)%rsvf < acsf(iread)%rsvf )  THEN
                        acsfnew(iwrite)%rtransp = acsf(iread)%rtransp
                    ENDIF
                    acsfnew(iwrite)%rsvf = acsfnew(iwrite)%rsvf + acsf(iread)%rsvf
!--                 advance reading index, keep writing index
                ELSE
!--                 not identical, just advance and copy
                    iwrite = iwrite + 1
                    acsfnew(iwrite) = acsf(iread)
                ENDIF
            ENDDO
            ncsfl = iwrite
        ENDIF

        IF ( newsize == -1 )  THEN
!--         allocate new array and copy shrinked data
            IF ( mcsf == 0 )  THEN
                ALLOCATE( acsf1(ncsfl) )
                acsf1(1:ncsfl) = acsf2(1:ncsfl)
            ELSE
                ALLOCATE( acsf2(ncsfl) )
                acsf2(1:ncsfl) = acsf1(1:ncsfl)
            ENDIF
        ENDIF

!--     deallocate old array
        IF ( mcsf == 0 )  THEN
            mcsf = 1
            acsf => acsf1
            DEALLOCATE( acsf2 )
        ELSE
            mcsf = 0
            acsf => acsf2
            DEALLOCATE( acsf1 )
        ENDIF
        ncsfla = newsize

!         WRITE(msg,'(A,2I12)') 'Grow acsf2:',ncsfl,ncsfla
!         CALL radiation_write_debug_log( msg )

    END SUBROUTINE merge_and_grow_csf

    
!-- quicksort.f -*-f90-*-
!-- Author: t-nissie, adaptation J.Resler
!-- License: GPLv3
!-- Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    RECURSIVE SUBROUTINE quicksort_csf2(kpcsflt, pcsflt, first, last)
        IMPLICIT NONE
        INTEGER(iwp), DIMENSION(:,:), INTENT(INOUT)  :: kpcsflt
        REAL(wp), DIMENSION(:,:), INTENT(INOUT)      :: pcsflt
        INTEGER(iwp), INTENT(IN)                     :: first, last
        REAL(wp), DIMENSION(ndcsf)                   :: t2
        INTEGER(iwp), DIMENSION(kdcsf)               :: x, t1
        INTEGER(iwp)                                 :: i, j

        IF ( first>=last ) RETURN
        x = kpcsflt(:, (first+last)/2 )
        i = first
        j = last
        DO
            DO while ( csf_lt2(kpcsflt(:,i),x) )
                i=i+1
            ENDDO
            DO while ( csf_lt2(x,kpcsflt(:,j)) )
                j=j-1
            ENDDO
            IF ( i >= j ) EXIT
            t1 = kpcsflt(:,i);  kpcsflt(:,i) = kpcsflt(:,j);  kpcsflt(:,j) = t1
            t2 = pcsflt(:,i);  pcsflt(:,i) = pcsflt(:,j);  pcsflt(:,j) = t2
            i=i+1
            j=j-1
        ENDDO
        IF ( first < i-1 ) CALL quicksort_csf2(kpcsflt, pcsflt, first, i-1)
        IF ( j+1 < last )  CALL quicksort_csf2(kpcsflt, pcsflt, j+1, last)
    END SUBROUTINE quicksort_csf2
    

    PURE FUNCTION csf_lt2(item1, item2) result(res)
        INTEGER(iwp), DIMENSION(kdcsf), INTENT(in)  :: item1, item2
        LOGICAL                                     :: res
        res = ( (item1(3) < item2(3))                                                        &
             .OR.  (item1(3) == item2(3)  .AND.  item1(2) < item2(2))                            &
             .OR.  (item1(3) == item2(3)  .AND.  item1(2) == item2(2)  .AND.  item1(1) < item2(1)) &
             .OR.  (item1(3) == item2(3)  .AND.  item1(2) == item2(2)  .AND.  item1(1) == item2(1) &
                 .AND.  item1(4) < item2(4)) )
    END FUNCTION csf_lt2

    PURE FUNCTION searchsorted(athresh, val) result(ind)
        REAL(wp), DIMENSION(:), INTENT(IN)  :: athresh
        REAL(wp), INTENT(IN)                :: val
        INTEGER(iwp)                        :: ind
        INTEGER(iwp)                        :: i

        DO i = LBOUND(athresh, 1), UBOUND(athresh, 1)
            IF ( val < athresh(i) ) THEN
                ind = i - 1
                RETURN
            ENDIF
        ENDDO
        ind = UBOUND(athresh, 1)
    END FUNCTION searchsorted

!------------------------------------------------------------------------------!
! Description:
! ------------
!
!> radiation_radflux_gridbox subroutine gives the sw and lw radiation fluxes at the 
!> faces of a gridbox defined at i,j,k and located in the urban layer.
!> The total sw and the diffuse sw radiation as well as the lw radiation fluxes at
!> the gridbox 6 faces are stored in sw_gridbox, swd_gridbox, and lw_gridbox arrays, 
!> respectively, in the following order:
!>  up_face, down_face, north_face, south_face, east_face, west_face
!> 
!> The subroutine reports also how successful was the search process via the parameter
!> i_feedback as follow:
!> - i_feedback =  1 : successful
!> - i_feedback = -1 : unsuccessful; the requisted point is outside the urban domain
!> - i_feedback =  0 : uncomplete; some gridbox faces fluxes are missing
!> 
!> 
!> It is called outside from usm_urban_surface_mod whenever the radiation fluxes
!> are needed.
!>
!> TODO:
!>    - Compare performance when using some combination of the Fortran intrinsic 
!>      functions, e.g. MINLOC, MAXLOC, ALL, ANY and COUNT functions, which search
!>      surfl array for elements meeting user-specified criterion, i.e. i,j,k
!>    - Report non-found or incomplete radiation fluxes arrays , if any, at the 
!>      gridbox faces in an error message form
!>
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_radflux_gridbox(i,j,k,sw_gridbox,swd_gridbox,lw_gridbox,i_feedback)
        
        IMPLICIT NONE

        INTEGER(iwp),                 INTENT(in)  :: i,j,k                 !< gridbox indices at which fluxes are required
        INTEGER(iwp)                              :: ii,jj,kk,d            !< surface indices and type
        INTEGER(iwp)                              :: l                     !< surface id
        REAL(wp)    , DIMENSION(1:6), INTENT(out) :: sw_gridbox,lw_gridbox !< total sw and lw radiation fluxes of 6 faces of a gridbox, w/m2 
        REAL(wp)    , DIMENSION(1:6), INTENT(out) :: swd_gridbox           !< diffuse sw radiation from sky and model boundary of 6 faces of a gridbox, w/m2
        INTEGER(iwp),                 INTENT(out) :: i_feedback            !< feedback to report how the search was successful 


!-- initialize variables
        i_feedback  = -999999
        sw_gridbox  = -999999.9_wp
        lw_gridbox  = -999999.9_wp
        swd_gridbox = -999999.9_wp
        
!-- check the requisted grid indices
        IF ( k < nzb   .OR.  k > nzut  .OR.   &
             j < nysg  .OR.  j > nyng  .OR.   &
             i < nxlg  .OR.  i > nxrg         &
             ) THEN
           i_feedback = -1
           RETURN
        ENDIF

!-- search for the required grid and formulate the fluxes at the 6 gridbox faces
        DO l = 1, nsurfl
            ii = surfl(ix,l)
            jj = surfl(iy,l)
            kk = surfl(iz,l)

            IF ( ii == i  .AND.  jj == j  .AND.  kk == k ) THEN
               d = surfl(id,l)

               SELECT CASE ( d )

               CASE (iup_u,iup_l,iup_a)                    !- gridbox up_facing face
                  sw_gridbox(1) = surfinsw(l)
                  lw_gridbox(1) = surfinlw(l)
                  swd_gridbox(1) = surfinswdif(l)

               CASE (idown_a)                         !- gridbox down_facing face
                  sw_gridbox(2) = surfinsw(l)
                  lw_gridbox(2) = surfinlw(l)
                  swd_gridbox(2) = surfinswdif(l)

               CASE (inorth_u,inorth_l,inorth_a)  !- gridbox north_facing face
                  sw_gridbox(3) = surfinsw(l)
                  lw_gridbox(3) = surfinlw(l)
                  swd_gridbox(3) = surfinswdif(l)

               CASE (isouth_u,isouth_l,isouth_a)  !- gridbox south_facing face
                  sw_gridbox(4) = surfinsw(l)
                  lw_gridbox(4) = surfinlw(l)
                  swd_gridbox(4) = surfinswdif(l)

               CASE (ieast_u,ieast_l,ieast_a)      !- gridbox east_facing face
                  sw_gridbox(5) = surfinsw(l)
                  lw_gridbox(5) = surfinlw(l)
                  swd_gridbox(5) = surfinswdif(l)

               CASE (iwest_u,iwest_l,iwest_a)      !- gridbox west_facing face
                  sw_gridbox(6) = surfinsw(l)
                  lw_gridbox(6) = surfinlw(l)
                  swd_gridbox(6) = surfinswdif(l)

               END SELECT

            ENDIF

        IF ( ALL( sw_gridbox(:)  /= -999999.9_wp )  ) EXIT
        ENDDO

!-- check the completeness of the fluxes at all gidbox faces        
!-- TODO: report non-found or incomplete rad fluxes arrays in an error message form
        IF ( ANY( sw_gridbox(:)  <= -999999.9_wp )  .OR.   &
             ANY( swd_gridbox(:) <= -999999.9_wp )  .OR.   &
             ANY( lw_gridbox(:)  <= -999999.9_wp ) ) THEN
           i_feedback = 0
        ELSE
           i_feedback = 1
        ENDIF
        
        RETURN
        
    END SUBROUTINE radiation_radflux_gridbox

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine for averaging 3D data
!------------------------------------------------------------------------------!
SUBROUTINE radiation_3d_data_averaging( mode, variable )
 

    USE control_parameters

    USE indices

    USE kinds

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  mode    !< 
    CHARACTER (LEN=*) :: variable !< 

    INTEGER(iwp) ::  i !< 
    INTEGER(iwp) ::  j !< 
    INTEGER(iwp) ::  k !< 
    INTEGER(iwp) ::  m !< index of current surface element 

    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( variable ) )

             CASE ( 'rad_net*' )
                IF ( .NOT. ALLOCATED( rad_net_av ) )  THEN
                   ALLOCATE( rad_net_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_net_av = 0.0_wp

             CASE ( 'rad_lw_in' )
                IF ( .NOT. ALLOCATED( rad_lw_in_av ) )  THEN
                   ALLOCATE( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_lw_in_av = 0.0_wp

             CASE ( 'rad_lw_out' )
                IF ( .NOT. ALLOCATED( rad_lw_out_av ) )  THEN
                   ALLOCATE( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_lw_out_av = 0.0_wp

             CASE ( 'rad_lw_cs_hr' )
                IF ( .NOT. ALLOCATED( rad_lw_cs_hr_av ) )  THEN
                   ALLOCATE( rad_lw_cs_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_lw_cs_hr_av = 0.0_wp

             CASE ( 'rad_lw_hr' )
                IF ( .NOT. ALLOCATED( rad_lw_hr_av ) )  THEN
                   ALLOCATE( rad_lw_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_lw_hr_av = 0.0_wp

             CASE ( 'rad_sw_in' )
                IF ( .NOT. ALLOCATED( rad_sw_in_av ) )  THEN
                   ALLOCATE( rad_sw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_sw_in_av = 0.0_wp

             CASE ( 'rad_sw_out' )
                IF ( .NOT. ALLOCATED( rad_sw_out_av ) )  THEN
                   ALLOCATE( rad_sw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_sw_out_av = 0.0_wp

             CASE ( 'rad_sw_cs_hr' )
                IF ( .NOT. ALLOCATED( rad_sw_cs_hr_av ) )  THEN
                   ALLOCATE( rad_sw_cs_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_sw_cs_hr_av = 0.0_wp

             CASE ( 'rad_sw_hr' )
                IF ( .NOT. ALLOCATED( rad_sw_hr_av ) )  THEN
                   ALLOCATE( rad_sw_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_sw_hr_av = 0.0_wp

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'rad_net*' )
             IF ( ALLOCATED( rad_net_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO m = surf_lsm_h%start_index(j,i), surf_lsm_h%end_index(j,i)
                         rad_net_av(j,i) = rad_net_av(j,i) + surf_lsm_h%rad_net(m)
                      ENDDO
                      DO m = surf_usm_h%start_index(j,i), surf_usm_h%end_index(j,i)
                         rad_net_av(j,i) = rad_net_av(j,i) + surf_usm_h%rad_net(m)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_in' )
             IF ( ALLOCATED( rad_lw_in_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_in_av(k,j,i) = rad_lw_in_av(k,j,i)             &
                                               + rad_lw_in(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_out' )
             IF ( ALLOCATED( rad_lw_out_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_out_av(k,j,i) = rad_lw_out_av(k,j,i)           &
                                                + rad_lw_out(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_cs_hr' )
             IF ( ALLOCATED( rad_lw_cs_hr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_cs_hr_av(k,j,i) = rad_lw_cs_hr_av(k,j,i)       &
                                                  + rad_lw_cs_hr(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_hr' )
             IF ( ALLOCATED( rad_lw_hr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_hr_av(k,j,i) = rad_lw_hr_av(k,j,i)             &
                                               + rad_lw_hr(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_in' )
             IF ( ALLOCATED( rad_sw_in_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_in_av(k,j,i) = rad_sw_in_av(k,j,i)             &
                                               + rad_sw_in(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_out' )
             IF ( ALLOCATED( rad_sw_out_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_out_av(k,j,i) = rad_sw_out_av(k,j,i)           &
                                                + rad_sw_out(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_cs_hr' )
             IF ( ALLOCATED( rad_sw_cs_hr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_cs_hr_av(k,j,i) = rad_sw_cs_hr_av(k,j,i)       &
                                                  + rad_sw_cs_hr(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_hr' )
             IF ( ALLOCATED( rad_sw_hr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_hr_av(k,j,i) = rad_sw_hr_av(k,j,i)             &
                                               + rad_sw_hr(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( variable ) )

         CASE ( 'rad_net*' )
             IF ( ALLOCATED( rad_net_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      rad_net_av(j,i) = rad_net_av(j,i)                        &
                                        / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_in' )
             IF ( ALLOCATED( rad_lw_in_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_in_av(k,j,i) = rad_lw_in_av(k,j,i)             &
                                               / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_out' )
             IF ( ALLOCATED( rad_lw_out_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_out_av(k,j,i) = rad_lw_out_av(k,j,i)           &
                                                / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_cs_hr' )
             IF ( ALLOCATED( rad_lw_cs_hr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_cs_hr_av(k,j,i) = rad_lw_cs_hr_av(k,j,i)       &
                                                / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_hr' )
             IF ( ALLOCATED( rad_lw_hr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_hr_av(k,j,i) = rad_lw_hr_av(k,j,i)             &
                                               / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_in' )
             IF ( ALLOCATED( rad_sw_in_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_in_av(k,j,i) = rad_sw_in_av(k,j,i)             &
                                               / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_out' )
             IF ( ALLOCATED( rad_sw_out_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_out_av(k,j,i) = rad_sw_out_av(k,j,i)           &
                                                / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_cs_hr' )
             IF ( ALLOCATED( rad_sw_cs_hr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_cs_hr_av(k,j,i) = rad_sw_cs_hr_av(k,j,i)       &
                                                / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_hr' )
             IF ( ALLOCATED( rad_sw_hr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_hr_av(k,j,i) = rad_sw_hr_av(k,j,i)             &
                                               / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

       END SELECT

    ENDIF

END SUBROUTINE radiation_3d_data_averaging


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!------------------------------------------------------------------------------!
SUBROUTINE radiation_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )
    
    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN)  ::  var         !< 
    LOGICAL, INTENT(OUT)           ::  found       !< 
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_x      !< 
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_y      !< 
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_z      !< 

    found  = .TRUE.


!
!-- Check for the grid
    SELECT CASE ( TRIM( var ) )

       CASE ( 'rad_lw_cs_hr', 'rad_lw_hr', 'rad_sw_cs_hr', 'rad_sw_hr',        &
              'rad_lw_cs_hr_xy', 'rad_lw_hr_xy', 'rad_sw_cs_hr_xy',            &
              'rad_sw_hr_xy', 'rad_lw_cs_hr_xz', 'rad_lw_hr_xz',               &
              'rad_sw_cs_hr_xz', 'rad_sw_hr_xz', 'rad_lw_cs_hr_yz',            &
              'rad_lw_hr_yz', 'rad_sw_cs_hr_yz', 'rad_sw_hr_yz' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu'

       CASE ( 'rad_lw_in', 'rad_lw_out', 'rad_sw_in', 'rad_sw_out',            &
              'rad_lw_in_xy', 'rad_lw_out_xy', 'rad_sw_in_xy','rad_sw_out_xy', &
              'rad_lw_in_xz', 'rad_lw_out_xz', 'rad_sw_in_xz','rad_sw_out_xz', &
              'rad_lw_in_yz', 'rad_lw_out_yz', 'rad_sw_in_yz','rad_sw_out_yz' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zw'


       CASE DEFAULT
          found  = .FALSE.
          grid_x = 'none'
          grid_y = 'none'
          grid_z = 'none'

        END SELECT

    END SUBROUTINE radiation_define_netcdf_grid

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 3D output variables
!------------------------------------------------------------------------------!
 SUBROUTINE radiation_data_output_2d( av, variable, found, grid, mode,         &
                                      local_pf, two_d, nzb_do, nzt_do )
 
    USE indices

    USE kinds


    IMPLICIT NONE

    CHARACTER (LEN=*) ::  grid     !< 
    CHARACTER (LEN=*) ::  mode     !< 
    CHARACTER (LEN=*) ::  variable !< 

    INTEGER(iwp) ::  av !< 
    INTEGER(iwp) ::  i  !< 
    INTEGER(iwp) ::  j  !< 
    INTEGER(iwp) ::  k  !< 
    INTEGER(iwp) ::  m  !< index of surface element at grid point (j,i)
    INTEGER(iwp) ::  nzb_do   !< 
    INTEGER(iwp) ::  nzt_do   !< 

    LOGICAL      ::  found !< 
    LOGICAL      ::  two_d !< flag parameter that indicates 2D variables (horizontal cross sections)

    REAL(wp) ::  fill_value = -999.0_wp    !< value for the _FillValue attribute

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf !< 

    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'rad_net*_xy' )        ! 2d-array
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Obtain rad_net from its respective surface type
!--                Natural-type surfaces
                   DO  m = surf_lsm_h%start_index(j,i),                        &
                           surf_lsm_h%end_index(j,i)  
                      local_pf(i,j,nzb+1) = surf_lsm_h%rad_net(m)
                   ENDDO
!
!--                Urban-type surfaces
                   DO  m = surf_usm_h%start_index(j,i),                        &
                           surf_usm_h%end_index(j,i)  
                      local_pf(i,j,nzb+1) = surf_usm_h%rad_net(m)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_net_av ) ) THEN
               ALLOCATE( rad_net_av(nysg:nyng,nxlg:nxrg) )
               rad_net_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn 
                   local_pf(i,j,nzb+1) = rad_net_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          two_d = .TRUE.
          grid = 'zu1'

 
       CASE ( 'rad_lw_in_xy', 'rad_lw_in_xz', 'rad_lw_in_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_lw_in(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_lw_in_av ) ) THEN
               ALLOCATE( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_in_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn 
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_lw_in_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'rad_lw_out_xy', 'rad_lw_out_xz', 'rad_lw_out_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_lw_out(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_lw_out_av ) ) THEN
               ALLOCATE( rad_lw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_out_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn 
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_lw_out_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF   
          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'rad_lw_cs_hr_xy', 'rad_lw_cs_hr_xz', 'rad_lw_cs_hr_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_lw_cs_hr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_lw_cs_hr_av ) ) THEN
               ALLOCATE( rad_lw_cs_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_cs_hr_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn 
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_lw_cs_hr_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'rad_lw_hr_xy', 'rad_lw_hr_xz', 'rad_lw_hr_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_lw_hr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_lw_hr_av ) ) THEN
               ALLOCATE( rad_lw_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_hr_av= REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn 
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_lw_hr_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'rad_sw_in_xy', 'rad_sw_in_xz', 'rad_sw_in_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_sw_in(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_sw_in_av ) ) THEN
               ALLOCATE( rad_sw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_in_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn 
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_sw_in_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'rad_sw_out_xy', 'rad_sw_out_xz', 'rad_sw_out_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_sw_out(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_sw_out_av ) ) THEN
               ALLOCATE( rad_sw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_out_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn 
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = rad_sw_out_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'rad_sw_cs_hr_xy', 'rad_sw_cs_hr_xz', 'rad_sw_cs_hr_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_sw_cs_hr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_sw_cs_hr_av ) ) THEN
               ALLOCATE( rad_sw_cs_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_cs_hr_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn 
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_sw_cs_hr_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'rad_sw_hr_xy', 'rad_sw_hr_xz', 'rad_sw_hr_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_sw_hr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_sw_hr_av ) ) THEN
               ALLOCATE( rad_sw_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_hr_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn 
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_sw_hr_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zw'

       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT
 
 END SUBROUTINE radiation_data_output_2d


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 3D output variables
!------------------------------------------------------------------------------!
 SUBROUTINE radiation_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )
 

    USE indices

    USE kinds


    IMPLICIT NONE

    CHARACTER (LEN=*) ::  variable !< 

    INTEGER(iwp) ::  av    !< 
    INTEGER(iwp) ::  i     !< 
    INTEGER(iwp) ::  j     !< 
    INTEGER(iwp) ::  k     !< 
    INTEGER(iwp) ::  nzb_do   !< 
    INTEGER(iwp) ::  nzt_do   !< 

    LOGICAL      ::  found !< 

    REAL(wp) ::  fill_value = -999.0_wp    !< value for the _FillValue attribute

    REAL(sp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf !< 


    found = .TRUE.


    SELECT CASE ( TRIM( variable ) )

      CASE ( 'rad_sw_in' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_sw_in(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_sw_in_av ) ) THEN
               ALLOCATE( rad_sw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_in_av = REAL( fill_value, KIND = wp )
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_sw_in_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_sw_out' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_sw_out(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_sw_out_av ) ) THEN
               ALLOCATE( rad_sw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_out_av = REAL( fill_value, KIND = wp )
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_sw_out_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_sw_cs_hr' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_sw_cs_hr(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_sw_cs_hr_av ) ) THEN
               ALLOCATE( rad_sw_cs_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_cs_hr_av = REAL( fill_value, KIND = wp )
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_sw_cs_hr_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_sw_hr' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_sw_hr(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_sw_hr_av ) ) THEN
               ALLOCATE( rad_sw_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_hr_av = REAL( fill_value, KIND = wp )
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_sw_hr_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_lw_in' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_lw_in(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_lw_in_av ) ) THEN
               ALLOCATE( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_in_av = REAL( fill_value, KIND = wp )
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_lw_in_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_lw_out' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_lw_out(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_lw_out_av ) ) THEN
               ALLOCATE( rad_lw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_out_av = REAL( fill_value, KIND = wp )
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_lw_out_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_lw_cs_hr' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_lw_cs_hr(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_lw_cs_hr_av ) ) THEN
               ALLOCATE( rad_lw_cs_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_cs_hr_av = REAL( fill_value, KIND = wp )
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_lw_cs_hr_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_lw_hr' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_lw_hr(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_lw_hr_av ) ) THEN
               ALLOCATE( rad_lw_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
              rad_lw_hr_av = REAL( fill_value, KIND = wp )
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_lw_hr_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

       CASE DEFAULT
          found = .FALSE.

    END SELECT


 END SUBROUTINE radiation_data_output_3d

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining masked data output
!------------------------------------------------------------------------------!
 SUBROUTINE radiation_data_output_mask( av, variable, found, local_pf )
 
    USE control_parameters
        
    USE indices
    
    USE kinds
    

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  variable   !< 

    INTEGER(iwp) ::  av   !< 
    INTEGER(iwp) ::  i    !< 
    INTEGER(iwp) ::  j    !< 
    INTEGER(iwp) ::  k    !< 

    LOGICAL ::  found     !< 

    REAL(wp),                                                                  &
       DIMENSION(mask_size_l(mid,1),mask_size_l(mid,2),mask_size_l(mid,3)) ::  &
          local_pf   !< 


    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )


       CASE ( 'rad_lw_in' )
          IF ( av == 0 )  THEN
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_lw_in(mask_k(mid,k),              &
                                            mask_j(mid,j),mask_i(mid,i))
                    ENDDO
                 ENDDO
              ENDDO
          ELSE
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_lw_in_av(mask_k(mid,k),           &
                                               mask_j(mid,j),mask_i(mid,i))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'rad_lw_out' )
          IF ( av == 0 )  THEN
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_lw_out(mask_k(mid,k),             &
                                            mask_j(mid,j),mask_i(mid,i))
                    ENDDO
                 ENDDO
              ENDDO
          ELSE
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_lw_out_av(mask_k(mid,k),          &
                                               mask_j(mid,j),mask_i(mid,i))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'rad_lw_cs_hr' )
          IF ( av == 0 )  THEN
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_lw_cs_hr(mask_k(mid,k),           &
                                            mask_j(mid,j),mask_i(mid,i))
                    ENDDO
                 ENDDO
              ENDDO
          ELSE
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_lw_cs_hr_av(mask_k(mid,k),        &
                                               mask_j(mid,j),mask_i(mid,i))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'rad_lw_hr' )
          IF ( av == 0 )  THEN
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_lw_hr(mask_k(mid,k),              &
                                            mask_j(mid,j),mask_i(mid,i))
                    ENDDO
                 ENDDO
              ENDDO
          ELSE
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_lw_hr_av(mask_k(mid,k),           &
                                               mask_j(mid,j),mask_i(mid,i))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'rad_sw_in' )
          IF ( av == 0 )  THEN
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_sw_in(mask_k(mid,k),              &
                                            mask_j(mid,j),mask_i(mid,i))
                    ENDDO
                 ENDDO
              ENDDO
          ELSE
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_sw_in_av(mask_k(mid,k),           &
                                               mask_j(mid,j),mask_i(mid,i))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'rad_sw_out' )
          IF ( av == 0 )  THEN
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_sw_out(mask_k(mid,k),             &
                                            mask_j(mid,j),mask_i(mid,i))
                    ENDDO
                 ENDDO
              ENDDO
          ELSE
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_sw_out_av(mask_k(mid,k),          &
                                               mask_j(mid,j),mask_i(mid,i))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'rad_sw_cs_hr' )
          IF ( av == 0 )  THEN
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_sw_cs_hr(mask_k(mid,k),           &
                                            mask_j(mid,j),mask_i(mid,i))
                    ENDDO
                 ENDDO
              ENDDO
          ELSE
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_sw_cs_hr_av(mask_k(mid,k),        &
                                               mask_j(mid,j),mask_i(mid,i))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'rad_sw_hr' )
          IF ( av == 0 )  THEN
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_sw_hr(mask_k(mid,k),              &
                                            mask_j(mid,j),mask_i(mid,i))
                    ENDDO
                 ENDDO
              ENDDO
          ELSE
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_sw_hr_av(mask_k(mid,k),           &
                                               mask_j(mid,j),mask_i(mid,i))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE DEFAULT
          found = .FALSE.

    END SELECT


 END SUBROUTINE radiation_data_output_mask


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine writes local (subdomain) restart data
!------------------------------------------------------------------------------!
 SUBROUTINE radiation_wrd_local


    IMPLICIT NONE


    IF ( ALLOCATED( rad_net_av ) )  THEN
       CALL wrd_write_string( 'rad_net_av' )
       WRITE ( 14 )  rad_net_av
    ENDIF

    IF ( ALLOCATED( rad_lw_in ) )  THEN
       CALL wrd_write_string( 'rad_lw_in' )
       WRITE ( 14 )  rad_lw_in
    ENDIF

    IF ( ALLOCATED( rad_lw_in_av ) )  THEN
       CALL wrd_write_string( 'rad_lw_in_av' )
       WRITE ( 14 )  rad_lw_in_av
    ENDIF

    IF ( ALLOCATED( rad_lw_out ) )  THEN
       CALL wrd_write_string( 'rad_lw_out' )
       WRITE ( 14 )  rad_lw_out
    ENDIF

    IF ( ALLOCATED( rad_lw_out_av) )  THEN
       CALL wrd_write_string( 'rad_lw_out_av' )
       WRITE ( 14 )  rad_lw_out_av
    ENDIF

    IF ( ALLOCATED( rad_lw_cs_hr) )  THEN
       CALL wrd_write_string( 'rad_lw_cs_hr' )
       WRITE ( 14 )  rad_lw_cs_hr
    ENDIF

    IF ( ALLOCATED( rad_lw_cs_hr_av) )  THEN
       CALL wrd_write_string( 'rad_lw_cs_hr_av' )
       WRITE ( 14 )  rad_lw_cs_hr_av
    ENDIF

    IF ( ALLOCATED( rad_lw_hr) )  THEN
       CALL wrd_write_string( 'rad_lw_hr' )
       WRITE ( 14 )  rad_lw_hr
    ENDIF

    IF ( ALLOCATED( rad_lw_hr_av) )  THEN
       CALL wrd_write_string( 'rad_lw_hr_av' )
       WRITE ( 14 )  rad_lw_hr_av
    ENDIF

    IF ( ALLOCATED( rad_sw_in) )  THEN
       CALL wrd_write_string( 'rad_sw_in' )
       WRITE ( 14 )  rad_sw_in
    ENDIF

    IF ( ALLOCATED( rad_sw_in_av) )  THEN
       CALL wrd_write_string( 'rad_sw_in_av' )
       WRITE ( 14 )  rad_sw_in_av
    ENDIF

    IF ( ALLOCATED( rad_sw_out) )  THEN
       CALL wrd_write_string( 'rad_sw_out' )
       WRITE ( 14 )  rad_sw_out
    ENDIF

    IF ( ALLOCATED( rad_sw_out_av) )  THEN
       CALL wrd_write_string( 'rad_sw_out_av' )
       WRITE ( 14 )  rad_sw_out_av
    ENDIF

    IF ( ALLOCATED( rad_sw_cs_hr) )  THEN
       CALL wrd_write_string( 'rad_sw_cs_hr' )
       WRITE ( 14 )  rad_sw_cs_hr
    ENDIF

    IF ( ALLOCATED( rad_sw_cs_hr_av) )  THEN
       CALL wrd_write_string( 'rad_sw_cs_hr_av' )
       WRITE ( 14 )  rad_sw_cs_hr_av
    ENDIF

    IF ( ALLOCATED( rad_sw_hr) )  THEN
       CALL wrd_write_string( 'rad_sw_hr' )
       WRITE ( 14 )  rad_sw_hr
    ENDIF

    IF ( ALLOCATED( rad_sw_hr_av) )  THEN
       CALL wrd_write_string( 'rad_sw_hr_av' )
       WRITE ( 14 )  rad_sw_hr_av
    ENDIF


 END SUBROUTINE radiation_wrd_local

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine reads local (subdomain) restart data
!------------------------------------------------------------------------------!
 SUBROUTINE radiation_rrd_local( i, k, nxlf, nxlc, nxl_on_file, nxrf, nxrc,    &
                                nxr_on_file, nynf, nync, nyn_on_file, nysf,    &
                                nysc, nys_on_file, tmp_2d, tmp_3d, found )
 

    USE control_parameters
        
    USE indices
    
    USE kinds
    
    USE pegrid


    IMPLICIT NONE

    INTEGER(iwp) ::  i               !< 
    INTEGER(iwp) ::  k               !< 
    INTEGER(iwp) ::  nxlc            !< 
    INTEGER(iwp) ::  nxlf            !< 
    INTEGER(iwp) ::  nxl_on_file     !< 
    INTEGER(iwp) ::  nxrc            !< 
    INTEGER(iwp) ::  nxrf            !< 
    INTEGER(iwp) ::  nxr_on_file     !< 
    INTEGER(iwp) ::  nync            !< 
    INTEGER(iwp) ::  nynf            !< 
    INTEGER(iwp) ::  nyn_on_file     !< 
    INTEGER(iwp) ::  nysc            !< 
    INTEGER(iwp) ::  nysf            !< 
    INTEGER(iwp) ::  nys_on_file     !< 

    LOGICAL, INTENT(OUT)  :: found

    REAL(wp), DIMENSION(nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_2d   !< 

    REAL(wp), DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_3d   !< 

    REAL(wp), DIMENSION(0:0,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_3d2   !< 


    found = .TRUE.


    SELECT CASE ( restart_string(1:length) )

       CASE ( 'rad_net_av' )
          IF ( .NOT. ALLOCATED( rad_net_av ) )  THEN
             ALLOCATE( rad_net_av(nysg:nyng,nxlg:nxrg) )
          ENDIF  
          IF ( k == 1 )  READ ( 13 )  tmp_2d
          rad_net_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =           &
                        tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
       CASE ( 'rad_lw_in' )
          IF ( .NOT. ALLOCATED( rad_lw_in ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant')  THEN
                ALLOCATE( rad_lw_in(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_lw_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF  
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant')  THEN
                READ ( 13 )  tmp_3d2
                rad_lw_in(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =   &
                   tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ ( 13 )  tmp_3d
                rad_lw_in(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =     &
                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_lw_in_av' )
          IF ( .NOT. ALLOCATED( rad_lw_in_av ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant')  THEN
                ALLOCATE( rad_lw_in_av(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF  
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant')  THEN
                READ ( 13 )  tmp_3d2
                rad_lw_in_av(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =&
                    tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ ( 13 )  tmp_3d
                rad_lw_in_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =  &
                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_lw_out' )
          IF ( .NOT. ALLOCATED( rad_lw_out ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant')  THEN
                ALLOCATE( rad_lw_out(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_lw_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF  
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant')  THEN
                READ ( 13 )  tmp_3d2
                rad_lw_out(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =  &
                    tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ ( 13 )  tmp_3d
                rad_lw_out(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =    &
                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_lw_out_av' )
          IF ( .NOT. ALLOCATED( rad_lw_out_av ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant')  THEN
                ALLOCATE( rad_lw_out_av(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_lw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF  
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant')  THEN
                READ ( 13 )  tmp_3d2
                rad_lw_out_av(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) &
                   = tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ ( 13 )  tmp_3d
                rad_lw_out_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_lw_cs_hr' )
          IF ( .NOT. ALLOCATED( rad_lw_cs_hr ) )  THEN
             ALLOCATE( rad_lw_cs_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          rad_lw_cs_hr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =        &
                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_lw_cs_hr_av' )
          IF ( .NOT. ALLOCATED( rad_lw_cs_hr_av ) )  THEN
             ALLOCATE( rad_lw_cs_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          rad_lw_cs_hr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =     &
                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_lw_hr' )
          IF ( .NOT. ALLOCATED( rad_lw_hr ) )  THEN
             ALLOCATE( rad_lw_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          rad_lw_hr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_lw_hr_av' )
          IF ( .NOT. ALLOCATED( rad_lw_hr_av ) )  THEN
             ALLOCATE( rad_lw_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          rad_lw_hr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =        &
                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_sw_in' )
          IF ( .NOT. ALLOCATED( rad_sw_in ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant')  THEN
                ALLOCATE( rad_sw_in(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_sw_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF  
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant')  THEN
                READ ( 13 )  tmp_3d2
                rad_sw_in(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =   &
                    tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ ( 13 )  tmp_3d
                rad_sw_in(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =     &
                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_sw_in_av' )
          IF ( .NOT. ALLOCATED( rad_sw_in_av ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant')  THEN
                ALLOCATE( rad_sw_in_av(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_sw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF  
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant')  THEN
                READ ( 13 )  tmp_3d2
                rad_sw_in_av(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =&
                    tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ ( 13 )  tmp_3d
                rad_sw_in_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =  &
                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_sw_out' )
          IF ( .NOT. ALLOCATED( rad_sw_out ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant')  THEN
                ALLOCATE( rad_sw_out(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_sw_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF  
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant')  THEN
                READ ( 13 )  tmp_3d2
                rad_sw_out(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =  &
                    tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ ( 13 )  tmp_3d
                rad_sw_out(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =    &
                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_sw_out_av' )
          IF ( .NOT. ALLOCATED( rad_sw_out_av ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant')  THEN
                ALLOCATE( rad_sw_out_av(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_sw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF  
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant')  THEN
                READ ( 13 )  tmp_3d2
                rad_sw_out_av(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) &
                   = tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ ( 13 )  tmp_3d
                rad_sw_out_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_sw_cs_hr' )
          IF ( .NOT. ALLOCATED( rad_sw_cs_hr ) )  THEN
             ALLOCATE( rad_sw_cs_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          rad_sw_cs_hr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =        &
                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_sw_cs_hr_av' )
          IF ( .NOT. ALLOCATED( rad_sw_cs_hr_av ) )  THEN
             ALLOCATE( rad_sw_cs_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          rad_sw_cs_hr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =     &
                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_sw_hr' )
          IF ( .NOT. ALLOCATED( rad_sw_hr ) )  THEN
             ALLOCATE( rad_sw_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          rad_sw_hr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_sw_hr_av' )
          IF ( .NOT. ALLOCATED( rad_sw_hr_av ) )  THEN
             ALLOCATE( rad_sw_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          rad_lw_hr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =        &
                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE DEFAULT

          found = .FALSE.

    END SELECT

 END SUBROUTINE radiation_rrd_local

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine writes debug information
!------------------------------------------------------------------------------!
 SUBROUTINE radiation_write_debug_log ( message )
    !> it writes debug log with time stamp
    CHARACTER(*)  :: message
    CHARACTER(15) :: dtc
    CHARACTER(8)  :: date
    CHARACTER(10) :: time
    CHARACTER(5)  :: zone
    CALL date_and_time(date, time, zone)
    dtc = date(7:8)//','//time(1:2)//':'//time(3:4)//':'//time(5:10)
    WRITE(9,'(2A)') dtc, TRIM(message)
    FLUSH(9)
 END SUBROUTINE radiation_write_debug_log

 END MODULE radiation_model_mod
