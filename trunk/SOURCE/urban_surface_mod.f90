!> @file urban_surface_mod.f90
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
! $Id: urban_surface_mod.f90 3091 2018-06-28 16:20:35Z suehring $
! - Limit aerodynamic resistance at vertical walls.
! - Add check for local roughness length not exceeding surface-layer height and
!   limit roughness length where necessary.
! 
! 3065 2018-06-12 07:03:02Z Giersch
! Unused array dxdir was removed, dz was replaced by dzu to consider vertical
! grid stretching
! 
! 3049 2018-05-29 13:52:36Z Giersch
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error message added
! 
! 3029 2018-05-23 12:19:17Z raasch
! bugfix: close unit 151 instead of 90
! 
! 3014 2018-05-09 08:42:38Z maronga
! Added pc_transpiration_rate
! 
! 2977 2018-04-17 10:27:57Z kanani
! Implement changes from branch radiation (r2948-2971) with minor modifications.
! (moh.hefny):
! Extended exn for all model domain height to avoid the need to get nzut.
! 
! 2963 2018-04-12 14:47:44Z suehring
! Introduce index for vegetation/wall, pavement/green-wall and water/window 
! surfaces, for clearer access of surface fraction, albedo, emissivity, etc. .
! 
! 2943 2018-04-03 16:17:10Z suehring
! Calculate exner function at all height levels and remove some un-used
! variables.
! 
! 2932 2018-03-26 09:39:22Z maronga
! renamed urban_surface_par to urban_surface_parameters
! 
! 2921 2018-03-22 15:05:23Z Giersch
! The activation of spinup has been moved to parin
! 
! 2920 2018-03-22 11:22:01Z kanani
! Remove unused pcbl, npcbl from ONLY list
! moh.hefny:
! Fixed bugs introduced by new structures and by moving radiation interaction
! into radiation_model_mod.f90.
! Bugfix: usm data output 3D didn't respect directions
! 
! 2906 2018-03-19 08:56:40Z Giersch
! Local variable ids has to be initialized with a value of -1 in 
! usm_average_3d_data
! 
! 2894 2018-03-15 09:17:58Z Giersch
! Calculations of the index range of the subdomain on file which overlaps with
! the current subdomain are already done in read_restart_data_mod,
! usm_read/write_restart_data have been renamed to usm_r/wrd_local, variable 
! named found has been introduced for checking if restart data was found, 
! reading of restart strings has been moved completely to 
! read_restart_data_mod, usm_rrd_local is already inside the overlap loop 
! programmed in read_restart_data_mod, SAVE attribute added where necessary, 
! deallocation and allocation of some arrays have been changed to take care of 
! different restart files that can be opened (index i), the marker *** end usm 
! *** is not necessary anymore, strings and their respective lengths are 
! written out and read now in case of restart runs to get rid of prescribed 
! character lengths
! 
! 2805 2018-02-14 17:00:09Z suehring
! Initialization of resistances.
! 
! 2797 2018-02-08 13:24:35Z suehring
! Comment concerning output of ground-heat flux added.
! 
! 2766 2018-01-22 17:17:47Z kanani
! Removed redundant commas, added some blanks
! 
! 2765 2018-01-22 11:34:58Z maronga
! Major bugfix in calculation of f_shf. Adjustment of roughness lengths in
! building_pars
! 
! 2750 2018-01-15 16:26:51Z knoop
! Move flag plant canopy to modules
! 
! 2737 2018-01-11 14:58:11Z kanani
! Removed unused variables t_surf_whole...
! 
! 2735 2018-01-11 12:01:27Z suehring
! resistances are saved in surface attributes
! 
! 2723 2018-01-05 09:27:03Z maronga
! Bugfix for spinups (end_time was increased twice in case of LSM + USM runs)
! 
! 2720 2018-01-02 16:27:15Z kanani
! Correction of comment
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2705 2017-12-18 11:26:23Z maronga
! Changes from last commit documented
! 
! 2703 2017-12-15 20:12:38Z maronga
! Workaround for calculation of r_a
!
! 2696 2017-12-14 17:12:51Z kanani
! - Change in file header (GPL part)
! - Bugfix in calculation of pt_surface and related fluxes. (BM)
! - Do not write surface temperatures onto pt array as this might cause 
!   problems with nesting. (MS)
! - Revised calculation of pt1 (now done in surface_layer_fluxes).
!   Bugfix, f_shf_window and f_shf_green were not set at vertical surface 
!   elements. (MS)
! - merged with branch ebsolver
!   green building surfaces do not evaporate yet
!   properties of green wall layers and window layers are taken from wall layers
!   this input data is missing. (RvT)
! - Merged with branch radiation (developed by Mohamed Salim)
! - Revised initialization. (MS)
! - Rename emiss_surf into emissivity, roughness_wall into z0, albedo_surf into 
!   albedo. (MS)
! - Move first call of usm_radiatin from usm_init to init_3d_model
! - fixed problem with near surface temperature
! - added near surface temperature t_surf_10cm_h(m), t_surf_10cm_v(l)%t(m)
! - does not work with temp profile including stability, ol
!   t_surf_10cm = pt1 now
! - merged with 2357 bugfix, error message for nopointer version
! - added indoor model coupling with wall heat flux
! - added green substrate/ dry vegetation layer for buildings
! - merged with 2232 new surface-type structure
! - added transmissivity of window tiles
! - added MOSAIK tile approach for 3 different surfaces (RvT)
! 
! 2583 2017-10-26 13:58:38Z knoop
! Bugfix: reverted MPI_Win_allocate_cptr introduction in last commit 
!
! 2582 2017-10-26 13:19:46Z hellstea
! Workaround for gnufortran compiler added in usm_calc_svf. CALL MPI_Win_allocate is
! replaced by CALL MPI_Win_allocate_cptr if defined ( __gnufortran ).
!
! 2544 2017-10-13 18:09:32Z maronga
! Date and time quantities are now read from date_and_time_mod. Solar constant is
! read from radiation_model_mod
! 
! 2516 2017-10-04 11:03:04Z suehring
! Remove tabs
! 
! 2514 2017-10-04 09:52:37Z suehring
! upper bounds of 3d output changed from nx+1,ny+1 to nx,ny
! no output of ghost layer data
! 
! 2350 2017-08-15 11:48:26Z kanani
! Bugfix and error message for nopointer version. 
! Additional "! defined(__nopointer)" as workaround to enable compilation of 
! nopointer version.
! 
! 2318 2017-07-20 17:27:44Z suehring
! Get topography top index via Function call 
! 
! 2317 2017-07-20 17:27:19Z suehring
! Bugfix: adjust output of shf. Added support for spinups
! 
! 2287 2017-06-15 16:46:30Z suehring
! Bugfix in determination topography-top index
! 
! 2269 2017-06-09 11:57:32Z suehring
! Enable restart runs with different number of PEs
! Bugfixes nopointer branch
! 
! 2258 2017-06-08 07:55:13Z suehring
! Bugfix, add pre-preprocessor directives to enable non-parrallel mode
! 
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments according to new surface-type structure. Remove usm_wall_heat_flux; 
! insteat, heat fluxes are directly applied in diffusion_s.
! 
! 2213 2017-04-24 15:10:35Z kanani
! Removal of output quantities usm_lad and usm_canopy_hr
! 
! 2209 2017-04-19 09:34:46Z kanani
! cpp switch __mpi3 removed,
! minor formatting,
! small bugfix for division by zero (Krc)
! 
! 2113 2017-01-12 13:40:46Z kanani
! cpp switch __mpi3 added for MPI-3 standard code (Ketelsen)
! 
! 2071 2016-11-17 11:22:14Z maronga
! Small bugfix (Resler)
! 
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean
! 
! 2024 2016-10-12 16:42:37Z kanani
! Bugfixes in deallocation of array plantt and reading of csf/csfsurf,
! optimization of MPI-RMA operations,
! declaration of pcbl as integer,
! renamed usm_radnet -> usm_rad_net, usm_canopy_khf -> usm_canopy_hr,
! splitted arrays svf -> svf & csf, svfsurf -> svfsurf & csfsurf,
! use of new control parameter varnamelength,
! added output variables usm_rad_ressw, usm_rad_reslw,
! minor formatting changes,
! minor optimizations.
! 
! 2011 2016-09-19 17:29:57Z kanani
! Major reformatting according to PALM coding standard (comments, blanks,
! alphabetical ordering, etc.),
! removed debug_prints, 
! removed auxiliary SUBROUTINE get_usm_info, instead, USM flag urban_surface is
! defined in MODULE control_parameters (modules.f90) to avoid circular 
! dependencies,
! renamed canopy_heat_flux to pc_heating_rate, as meaning of quantity changed.
! 
! 2007 2016-08-24 15:47:17Z kanani
! Initial revision
!
!
! Description:
! ------------
! 2016/6/9 - Initial version of the USM (Urban Surface Model)
!            authors: Jaroslav Resler, Pavel Krc
!                     (Czech Technical University in Prague and Institute of
!                      Computer Science of the Czech Academy of Sciences, Prague)
!            with contributions: Michal Belda, Nina Benesova, Ondrej Vlcek
!            partly inspired by PALM LSM (B. Maronga)
!            parameterizations of Ra checked with TUF3D (E. S. Krayenhoff)
!> Module for Urban Surface Model (USM)
!> The module includes:
!>    1. radiation model with direct/diffuse radiation, shading, reflections
!>       and integration with plant canopy
!>    2. wall and wall surface model
!>    3. surface layer energy balance
!>    4. anthropogenic heat (only from transportation so far)
!>    5. necessary auxiliary subroutines (reading inputs, writing outputs,
!>       restart simulations, ...)
!> It also make use of standard radiation and integrates it into
!> urban surface model.
!>
!> Further work:
!> -------------
!> 1. Remove global arrays surfouts, surfoutl and only keep track of radiosity
!>    from surfaces that are visible from local surfaces (i.e. there is a SVF
!>    where target is local). To do that, radiosity will be exchanged after each
!>    reflection step using MPI_Alltoall instead of current MPI_Allgather.
!>
!> 2. Temporarily large values of surface heat flux can be observed, up to
!>    1.2 Km/s, which seem to be not realistic.
!>
!> @todo Output of _av variables in case of restarts
!> @todo Revise flux conversion in energy-balance solver
!> @todo Bugfixing in nopointer branch 
!> @todo Check optimizations for RMA operations
!> @todo Alternatives for MPI_WIN_ALLOCATE? (causes problems with openmpi)
!> @todo Check for load imbalances in CPU measures, e.g. for exchange_horiz_prog
!>       factor 3 between min and max time
!> @todo Move setting of flag indoor_model to indoor_model_mod once available
!> @todo Check divisions in wtend (etc.) calculations for possible division
!>       by zero, e.g. in case fraq(0,m) + fraq(1,m) = 0?!
!> @todo Use unit 90 for OPEN/CLOSE of input files (FK)
!> @todo Move plant canopy stuff into plant canopy code
!------------------------------------------------------------------------------!
 MODULE urban_surface_mod

#if ! defined( __nopointer )
    USE arrays_3d,                                                             &
        ONLY:  dzu, hyp, zu, pt, pt_1, pt_2, p, u, v, w, hyp, tend
#endif

    USE cloud_parameters,                                                      &
        ONLY:  cp, r_d

    USE constants,                                                             &
        ONLY:  pi
    
    USE control_parameters,                                                    &
        ONLY:  coupling_start_time, topography, dt_3d,                         &
               intermediate_timestep_count, initializing_actions,              &
               intermediate_timestep_count_max, simulated_time, end_time,      &
               timestep_scheme, tsc, coupling_char, io_blocks, io_group,       &
               message_string, time_since_reference_point, surface_pressure,   &
               g, pt_surface, large_scale_forcing, lsf_surf, spinup,           &
               spinup_pt_mean, spinup_time, time_do3d, dt_do3d,                &
               average_count_3d, varnamelength, urban_surface, kappa,          &
               plant_canopy

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s

    USE date_and_time_mod,                                                     &
        ONLY:  time_utc_init

    USE grid_variables,                                                        &
        ONLY:  dx, dy, ddx, ddy, ddx2, ddy2

    USE indices,                                                               &
        ONLY:  nx, ny, nnx, nny, nnz, nxl, nxlg, nxr, nxrg, nyn, nyng, nys,    &
               nysg, nzb, nzt, nbgp, wall_flags_0

    USE, INTRINSIC :: iso_c_binding 

    USE kinds
              
    USE pegrid
    
    USE plant_canopy_model_mod,                                                &
        ONLY:  pc_heating_rate, pc_transpiration_rate
    
    USE radiation_model_mod,                                                   &
        ONLY:  albedo_type, radiation_interaction, calc_zenith, zenith,        &
               radiation, rad_sw_in, rad_lw_in, rad_sw_out, rad_lw_out,        &
               sigma_sb, solar_constant, sun_direction, sun_dir_lat,           &
               sun_dir_lon,                                                    &
               force_radiation_call, surfinsw, surfinlw, surfinswdir,          &
               surfinswdif, surfoutsw, surfoutlw, surfins,nsvfl, svf, svfsurf, &
               surfinl, surfinlwdif, rad_sw_in_dir, rad_sw_in_diff,            &
               rad_lw_in_diff, surfouts, surfoutl, surfoutsl, surfoutll, surf, &
               surfl, nsurfl, nsurfs, surfstart, pcbinsw, pcbinlw,             &
               iup_u, inorth_u, isouth_u, ieast_u, iwest_u, iup_l,             &
               inorth_l, isouth_l, ieast_l, iwest_l, id,                       &
               iz, iy, ix, idir, jdir, kdir,  nsurf_type, nsurf, idsvf, ndsvf, &
               iup_a, idown_a, inorth_a, isouth_a, ieast_a, iwest_a,           &
               idcsf, ndcsf, kdcsf, pct,                                       &
               startland, endland, startwall, endwall, skyvf, skyvft

    USE statistics,                                                            &
        ONLY:  hom, statistic_regions

    USE surface_mod,                                                           &
        ONLY:  get_topography_top_index_ji, get_topography_top_index,          &
               ind_pav_green, ind_veg_wall, ind_wat_win, surf_usm_h,           &
               surf_usm_v, surface_restore_elements


    IMPLICIT NONE


!-- configuration parameters (they can be setup in PALM config)
    LOGICAL                                        ::  usm_material_model = .TRUE.        !< flag parameter indicating wheather the  model of heat in materials is used
    LOGICAL                                        ::  usm_anthropogenic_heat = .FALSE.   !< flag parameter indicating wheather the anthropogenic heat sources (e.g.transportation) are used
    LOGICAL                                        ::  force_radiation_call_l = .FALSE.   !< flag parameter for unscheduled radiation model calls
    LOGICAL                                        ::  indoor_model = .FALSE.              !< whether to use the indoor model
    LOGICAL                                        ::  read_wall_temp_3d = .FALSE.


    INTEGER(iwp)                                   ::  building_type = 1                  !< default building type (preleminary setting)
    INTEGER(iwp)                                   ::  land_category = 2                  !< default category for land surface
    INTEGER(iwp)                                   ::  wall_category = 2                  !< default category for wall surface over pedestrian zone
    INTEGER(iwp)                                   ::  pedestrian_category = 2            !< default category for wall surface in pedestrian zone
    INTEGER(iwp)                                   ::  roof_category = 2                  !< default category for root surface
    REAL(wp)                                       ::  roughness_concrete = 0.001_wp      !< roughness length of average concrete surface
!
!-- Indices of input attributes for (above) ground floor level 
    INTEGER(iwp) ::  ind_alb_wall          = 38 !< index in input list for albedo_type of wall fraction
    INTEGER(iwp) ::  ind_alb_green         = 39 !< index in input list for albedo_type of green fraction
    INTEGER(iwp) ::  ind_alb_win           = 40 !< index in input list for albedo_type of window fraction
    INTEGER(iwp) ::  ind_emis_wall_agfl    = 14 !< index in input list for wall emissivity, above ground floor level
    INTEGER(iwp) ::  ind_emis_wall_gfl     = 32 !< index in input list for wall emissivity, ground floor level
    INTEGER(iwp) ::  ind_emis_green_agfl   = 15 !< index in input list for green emissivity, above ground floor level
    INTEGER(iwp) ::  ind_emis_green_gfl    = 33 !< index in input list for green emissivity, ground floor level
    INTEGER(iwp) ::  ind_emis_win_agfl     = 16 !< index in input list for window emissivity, above ground floor level
    INTEGER(iwp) ::  ind_emis_win_gfl      = 34 !< index in input list for window emissivity, ground floor level
    INTEGER(iwp) ::  ind_green_frac_w_agfl = 2  !< index in input list for green fraction on wall, above ground floor level
    INTEGER(iwp) ::  ind_green_frac_w_gfl  = 23 !< index in input list for green fraction on wall, ground floor level
    INTEGER(iwp) ::  ind_green_frac_r_agfl = 3  !< index in input list for green fraction on roof, above ground floor level
    INTEGER(iwp) ::  ind_green_frac_r_gfl  = 24 !< index in input list for green fraction on roof, ground floor level
    INTEGER(iwp) ::  ind_hc1_agfl          =  6 !< index in input list for heat capacity at first wall layer, above ground floor level
    INTEGER(iwp) ::  ind_hc1_gfl           = 26 !< index in input list for heat capacity at first wall layer, ground floor level
    INTEGER(iwp) ::  ind_hc2_agfl          = 7  !< index in input list for heat capacity at second wall layer, above ground floor level
    INTEGER(iwp) ::  ind_hc2_gfl           = 27 !< index in input list for heat capacity at second wall layer, ground floor level
    INTEGER(iwp) ::  ind_hc3_agfl          = 8  !< index in input list for heat capacity at third wall layer, above ground floor level
    INTEGER(iwp) ::  ind_hc3_gfl           = 28 !< index in input list for heat capacity at third wall layer, ground floor level
    INTEGER(iwp) ::  ind_gflh              = 20 !< index in input list for ground floor level height
    INTEGER(iwp) ::  ind_lai_r_agfl        = 4  !< index in input list for LAI on roof, above ground floor level
    INTEGER(iwp) ::  ind_lai_r_gfl         = 4  !< index in input list for LAI on roof, ground floor level
    INTEGER(iwp) ::  ind_lai_w_agfl        = 5  !< index in input list for LAI on wall, above ground floor level
    INTEGER(iwp) ::  ind_lai_w_gfl         = 25 !< index in input list for LAI on wall, ground floor level
    INTEGER(iwp) ::  ind_tc1_agfl          = 9  !< index in input list for thermal conductivity at first wall layer, above ground floor level
    INTEGER(iwp) ::  ind_tc1_gfl           = 29 !< index in input list for thermal conductivity at first wall layer, ground floor level
    INTEGER(iwp) ::  ind_tc2_agfl          = 10 !< index in input list for thermal conductivity at second wall layer, above ground floor level
    INTEGER(iwp) ::  ind_tc2_gfl           = 30 !< index in input list for thermal conductivity at second wall layer, ground floor level
    INTEGER(iwp) ::  ind_tc3_agfl          = 11 !< index in input list for thermal conductivity at third wall layer, above ground floor level
    INTEGER(iwp) ::  ind_tc3_gfl           = 31 !< index in input list for thermal conductivity at third wall layer, ground floor level
    INTEGER(iwp) ::  ind_thick_1           = 41 !< index for wall layer thickness - 1st layer
    INTEGER(iwp) ::  ind_thick_2           = 42 !< index for wall layer thickness - 2nd layer
    INTEGER(iwp) ::  ind_thick_3           = 43 !< index for wall layer thickness - 3rd layer
    INTEGER(iwp) ::  ind_thick_4           = 44 !< index for wall layer thickness - 4th layer
    INTEGER(iwp) ::  ind_trans_agfl        = 17 !< index in input list for window transmissivity, above ground floor level
    INTEGER(iwp) ::  ind_trans_gfl         = 35 !< index in input list for window transmissivity, ground floor level
    INTEGER(iwp) ::  ind_wall_frac_agfl    = 0  !< index in input list for wall fraction, above ground floor level
    INTEGER(iwp) ::  ind_wall_frac_gfl     = 21 !< index in input list for wall fraction, ground floor level
    INTEGER(iwp) ::  ind_win_frac_agfl     = 1  !< index in input list for window fraction, above ground floor level
    INTEGER(iwp) ::  ind_win_frac_gfl      = 22 !< index in input list for window fraction, ground floor level
    INTEGER(iwp) ::  ind_z0_agfl           = 18 !< index in input list for z0, above ground floor level
    INTEGER(iwp) ::  ind_z0_gfl            = 36 !< index in input list for z0, ground floor level
    INTEGER(iwp) ::  ind_z0qh_agfl         = 19 !< index in input list for z0h / z0q, above ground floor level
    INTEGER(iwp) ::  ind_z0qh_gfl          = 37 !< index in input list for z0h / z0q, ground floor level


    REAL(wp)  ::  roof_height_limit = 4._wp          !< height for distinguish between land surfaces and roofs
    REAL(wp)  ::  ground_floor_level = 4.0_wp        !< default ground floor level


    CHARACTER(37), DIMENSION(0:6), PARAMETER :: building_type_name = (/     &
                                   'user-defined                         ', & !  0 
                                   'residential - 1950                   ', & !  1 
                                   'residential 1951 - 2000              ', & !  2
                                   'residential 2001 -                   ', & !  3
                                   'office - 1950                        ', & !  4 
                                   'office 1951 - 2000                   ', & !  5
                                   'office 2001 -                        '  & !  6
                                                                     /)
!
!-- building parameters, 4 different types
!-- 0 - wall fraction, 1- window fraction, 2 - green fraction on wall, 3- green fraction 
!-- at roof, 4 - lai of green fraction at roof,  5 - lai of green fraction at wall, 
!-- 6 - heat capacity of wall layer 1, 7 - heat capacity of wall layer 2, 
!-- 8 - heat capacity of wall layer 3, 9 - thermal conductivity of wall layer 1, 
!-- 10 - thermal conductivity of wall layer 2, 11 - thermal conductivity of wall layer 3,  
!-- 12 - indoor target summer temperature ( K ), 13 - indoor target winter temperature (K),
!-- 14 - emissivity of wall fraction, 15 - emissivity of green fraction, 16 - emissivity of window fraction,
!-- 17 - transmissivity of window fraction, 18 - z0, 19 - z0h/z0q, 20 - ground floor height, 
!-- 21 - ground floor wall fraction, 22 - ground floor window fraction, 23 ground floor green fraction,
!-- 24 - ground floor green fraction on roof, 25 - ground floor lai of green fraction, 
!-- 26 - ground floor heat capacity of wall layer 1, 27 - ground floor heat capacity of wall layer 1,
!-- 28 - ground floor heat capacity of wall layer 3, 29 - ground floor thermal conductivity of wall layer 1,
!-- 30 - ground floor thermal conductivity of wall layer 2, 31 - ground floor thermal conductivity of wall layer 3,
!-- 32 - ground floor emissivity of wall fraction, 33 - ground floor emissivity of green fraction,
!-- 34 - ground floor emissivity of window fraction, 35 - ground floor transmissivity of window fraction, 
!-- 36 - ground floor z0, 37 - ground floor z0h/z0q, 38 - albedo type wall fraction
!-- 39 - albedo type green fraction, 40 - albedo type window fraction
!-- 41 - wall layer thickness - 1st layer, 42 - wall layer thickness - 2nd layer,
!-- 43 - wall layer thickness - 3rd layer, 44 - wall layer thickness - 4th layer,
!-- 45 - heat capacity of the wall surface, 46 - heat conductivity
!-- Please note, only preleminary dummy values so far!
    REAL(wp), DIMENSION(0:46,1:6), PARAMETER :: building_pars = RESHAPE( (/    &
        1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp,                        & !parameter 0-5
        1000000.0_wp, 1000000.0_wp, 1000000.0_wp, 0.3_wp, 0.3_wp, 0.3_wp,      & !parameter 6-11
        296.15_wp, 293.15_wp, 0.9_wp, 0.9_wp, 0.01_wp, 0.99_wp,                & !parameter 12-17
        0.001_wp, 0.0001_wp, 4.0_wp,                                             & !parameter 18-20
        1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp,                                & !parameter 21-25
        1000000.0_wp, 1000000.0_wp, 1000000.0_wp,                              & !parameter 26-28                     
        0.3_wp, 0.3_wp, 0.3_wp,                                                & !parameter 29-31       
        0.4_wp, 0.4_wp, 0.4_wp, 0.4_wp, 0.01_wp, 0.001_wp,                     & !parameter 32-37
        24.0_wp, 24.0_wp, 24.0_wp,                                             & !parameter 38-40
        0.0242_wp, 0.0969_wp, 0.346_wp, 1.0_wp,                                & !parameter 41-44
        20000.0_wp, 10.0_wp,                                                   & !parameter 45-46 - end of type 1
        1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp,                        & !parameter 0-5
        1000000.0_wp, 1000000.0_wp, 1000000.0_wp, 0.3_wp, 0.3_wp, 0.3_wp,      & !parameter 6-11
        296.15_wp, 293.15_wp, 0.9_wp, 0.9_wp, 0.01_wp, 0.99_wp,                & !parameter 12-17
        0.001_wp, 0.0001_wp, 4.0_wp,                                             & !parameter 18-20
        1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp,                                & !parameter 21-25
        1000000.0_wp, 1000000.0_wp, 1000000.0_wp,                              & !parameter 26-28                     
        0.3_wp, 0.3_wp, 0.3_wp,                                                & !parameter 29-31       
        0.4_wp, 0.4_wp, 0.4_wp, 0.4_wp, 0.01_wp, 0.001_wp,                     & !parameter 32-37
        24.0_wp, 24.0_wp, 24.0_wp,                                             & !parameter 38-40
        0.0242_wp, 0.0969_wp, 0.346_wp, 1.0_wp,                                & !parameter 41-44
        20000.0_wp, 10.0_wp,                                                   & !parameter 45-46 - end of type 2
        1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp,                        & !parameter 0-5
        1000000.0_wp, 1000000.0_wp, 1000000.0_wp, 0.3_wp, 0.3_wp, 0.3_wp,      & !parameter 6-11
        296.15_wp, 293.15_wp, 0.9_wp, 0.9_wp, 0.01_wp, 0.99_wp,                & !parameter 12-17
        0.001_wp, 0.0001_wp, 4.0_wp,                                             & !parameter 18-20
        1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp,                                & !parameter 21-25
        1000000.0_wp, 1000000.0_wp, 1000000.0_wp,                              & !parameter 26-28                     
        0.3_wp, 0.3_wp, 0.3_wp,                                                & !parameter 29-31       
        0.4_wp, 0.4_wp, 0.4_wp, 0.4_wp, 0.01_wp, 0.001_wp,                     & !parameter 32-37
        24.0_wp, 24.0_wp, 24.0_wp,                                             & !parameter 38-40
        0.0242_wp, 0.0969_wp, 0.346_wp, 1.0_wp,                                & !parameter 41-44
        20000.0_wp, 10.0_wp,                                                   & !parameter 45-46 - end of type 3
        1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp,                        & !parameter 0-5
        1000000.0_wp, 1000000.0_wp, 1000000.0_wp, 0.3_wp, 0.3_wp, 0.3_wp,      & !parameter 6-11
        296.15_wp, 293.15_wp, 0.9_wp, 0.9_wp, 0.01_wp, 0.99_wp,                & !parameter 12-17
        0.01_wp, 0.001_wp, 4.0_wp,                                             & !parameter 18-20
        1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp,                                & !parameter 21-25
        1000000.0_wp, 1000000.0_wp, 1000000.0_wp,                              & !parameter 26-28                     
        0.3_wp, 0.3_wp, 0.3_wp,                                                & !parameter 29-31       
        0.4_wp, 0.4_wp, 0.4_wp, 0.4_wp, 0.01_wp, 0.001_wp,                     & !parameter 32-37
        24.0_wp, 24.0_wp, 24.0_wp,                                             & !parameter 38-40
        0.0242_wp, 0.0969_wp, 0.346_wp, 1.0_wp,                                & !parameter 41-44
        20000.0_wp, 10.0_wp,                                                   & !parameter 45-46 - end of type 4
        1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp,                        & !parameter 0-5
        1000000.0_wp, 1000000.0_wp, 1000000.0_wp, 0.3_wp, 0.3_wp, 0.3_wp,      & !parameter 6-11
        296.15_wp, 293.15_wp, 0.9_wp, 0.9_wp, 0.01_wp, 0.99_wp,                & !parameter 12-17
        0.001_wp, 0.0001_wp, 4.0_wp,                                             & !parameter 18-20
        1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp,                                & !parameter 21-25
        1000000.0_wp, 1000000.0_wp, 1000000.0_wp,                              & !parameter 26-28                     
        0.3_wp, 0.3_wp, 0.3_wp,                                                & !parameter 29-31       
        0.4_wp, 0.4_wp, 0.4_wp, 0.4_wp, 0.01_wp, 0.001_wp,                     & !parameter 32-37
        24.0_wp, 24.0_wp, 24.0_wp,                                             & !parameter 38-40
        0.0242_wp, 0.0969_wp, 0.346_wp, 1.0_wp,                                & !parameter 41-44
        20000.0_wp, 10.0_wp,                                                   & !parameter 45-46 - end of type 5
        1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp,                        & !parameter 0-5
        1000000.0_wp, 1000000.0_wp, 1000000.0_wp, 0.3_wp, 0.3_wp, 0.3_wp,      & !parameter 6-11
        296.15_wp, 293.15_wp, 0.9_wp, 0.9_wp, 0.01_wp, 0.99_wp,                & !parameter 12-17
        0.001_wp, 0.0001_wp, 4.0_wp,                                             & !parameter 18-20
        1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp,                                & !parameter 21-25
        1000000.0_wp, 1000000.0_wp, 1000000.0_wp,                              & !parameter 26-28                     
        0.3_wp, 0.3_wp, 0.3_wp,                                                & !parameter 29-31       
        0.4_wp, 0.4_wp, 0.4_wp, 0.4_wp, 0.01_wp, 0.001_wp,                     & !parameter 32-37
        24.0_wp, 24.0_wp, 24.0_wp,                                             & !parameter 38-40
        0.0242_wp, 0.0969_wp, 0.346_wp, 1.0_wp,                                & !parameter 41-44
        20000.0_wp, 10.0_wp                                                    & !parameter 45-46 - end of type 6
                                                                          /),  &
                                                               (/47, 6/) )

!
!-- Type for surface temperatures at vertical walls. Is not necessary for horizontal walls. 
    TYPE t_surf_vertical
       REAL(wp), DIMENSION(:), ALLOCATABLE         :: t
    END TYPE t_surf_vertical
!
!-- Type for wall temperatures at vertical walls. Is not necessary for horizontal walls. 
    TYPE t_wall_vertical
       REAL(wp), DIMENSION(:,:), ALLOCATABLE       :: t
    END TYPE t_wall_vertical


!-- arrays for time averages
!-- Attention: the variable rad_net_av is also used in the 3d field variable in radiation_model_mod.f90. It may be better to rename it 
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  rad_net_av       !< average of rad_net_l
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinsw_av      !< average of sw radiation falling to local surface including radiation from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinlw_av      !< average of lw radiation falling to local surface including radiation from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinswdir_av   !< average of direct sw radiation falling to local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinswdif_av   !< average of diffuse sw radiation from sky and model boundary falling to local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinlwdif_av   !< average of diffuse lw radiation from sky and model boundary falling to local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinswref_av   !< average of sw radiation falling to surface from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinlwref_av   !< average of lw radiation falling to surface from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutsw_av     !< average of total sw radiation outgoing from nonvirtual surfaces surfaces after all reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutlw_av     !< average of total lw radiation outgoing from nonvirtual surfaces surfaces after all reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfins_av       !< average of array of residua of sw radiation absorbed in surface after last reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinl_av       !< average of array of residua of lw radiation absorbed in surface after last reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfhf_av        !< average of total radiation flux incoming to minus outgoing from local surface  
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  wghf_eb_av       !< average of wghf_eb
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  wshf_eb_av       !< average of wshf_eb
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  t_wall_av        !< Average of t_wall
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  wghf_eb_green_av !< average of wghf_eb_green
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  t_green_av       !< Average of t_green
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  wghf_eb_window_av !< average of wghf_eb_window
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  t_window_av      !< Average of t_window    
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-- anthropogenic heat sources
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE        ::  aheat             !< daily average of anthropogenic heat (W/m2)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  aheatprof         !< diurnal profiles of anthropogenic heat for particular layers
    INTEGER(wp)                                    ::  naheatlayers = 1  !< number of layers of anthropogenic heat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-- wall surface model
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-- wall surface model constants
    INTEGER(iwp), PARAMETER                        :: nzb_wall = 0       !< inner side of the wall model (to be switched)
    INTEGER(iwp), PARAMETER                        :: nzt_wall = 3       !< outer side of the wall model (to be switched)
    INTEGER(iwp), PARAMETER                        :: nzw = 4            !< number of wall layers (fixed for now)

    REAL(wp), DIMENSION(nzb_wall:nzt_wall)         :: zwn_default = (/0.0242_wp, 0.0969_wp, 0.346_wp, 1.0_wp /)
                                                                         !< normalized soil, wall and roof layer depths (m/m)
!    REAL(wp), DIMENSION(nzb_wall:nzt_wall)         :: zwn_default = (/0.33_wp, 0.66_wp, 1.0_wp /)
    REAL(wp), DIMENSION(nzb_wall:nzt_wall)         :: zwn_default_window = (/0.25_wp, 0.5_wp, 0.75_wp, 1.0_wp /)
!    REAL(wp), DIMENSION(nzb_wall:nzt_wall)         :: zwn_default_window = (/0.33_wp, 0.66_wp, 1.0_wp /)
!    REAL(wp), DIMENSION(nzb_wall:nzt_wall)         :: zwn_default_window = (/0.0242_wp, 0.0969_wp, 0.346_wp, 1.0_wp /)
                                                                         !< normalized window layer depths (m/m)
!    REAL(wp), DIMENSION(nzb_wall:nzt_wall)         :: zwn_default_green = (/0.0242_wp, 0.0969_wp, 0.346_wp, 1.0_wp /)
                                                                         !< normalized green layer depths (m/m)
    REAL(wp), DIMENSION(nzb_wall:nzt_wall)         :: zwn_default_green = (/0.25_wp, 0.5_wp, 0.75_wp, 1.0_wp /)
!    REAL(wp), DIMENSION(nzb_wall:nzt_wall)         :: zwn_default_green = (/0.33_wp, 0.66_wp, 1.0_wp /)


    REAL(wp)                                       :: wall_inner_temperature = 295.0_wp    !< temperature of the inner wall surface (~22 degrees C) (K)
    REAL(wp)                                       :: roof_inner_temperature = 295.0_wp    !< temperature of the inner roof surface (~22 degrees C) (K)
    REAL(wp)                                       :: soil_inner_temperature = 288.0_wp    !< temperature of the deep soil (~15 degrees C) (K)
    REAL(wp)                                       :: window_inner_temperature = 295.0_wp  !< temperature of the inner window surface (~22 degrees C) (K)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-- surface and material model variables for walls, ground, roofs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: zwn                !< normalized wall layer depths (m)
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: zwn_window         !< normalized window layer depths (m)
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: zwn_green          !< normalized green layer depths (m)

#if defined( __nopointer )
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_h           !< wall surface temperature (K) at horizontal walls
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_h_p         !< progn. wall surface temperature (K) at horizontal walls
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_window_h    !< window surface temperature (K) at horizontal walls
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_window_h_p  !< progn. window surface temperature (K) at horizontal walls
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_green_h     !< green surface temperature (K) at horizontal walls
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_green_h_p   !< progn. green surface temperature (K) at horizontal walls
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_10cm_h      !< near surface temperature (10cm) (K) at horizontal walls
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_10cm_h_p    !< progn. near surface temperature (10cm) (K) at horizontal walls
    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  ::  t_surf_v
    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  ::  t_surf_v_p
    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  ::  t_surf_window_v
    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  ::  t_surf_window_v_p
    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  ::  t_surf_green_v
    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  ::  t_surf_green_v_p
    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  ::  t_surf_10cm_v
    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  ::  t_surf_10cm_v_p
#else
    REAL(wp), DIMENSION(:), POINTER                :: t_surf_h
    REAL(wp), DIMENSION(:), POINTER                :: t_surf_h_p 
    REAL(wp), DIMENSION(:), POINTER                :: t_surf_window_h
    REAL(wp), DIMENSION(:), POINTER                :: t_surf_window_h_p 
    REAL(wp), DIMENSION(:), POINTER                :: t_surf_green_h
    REAL(wp), DIMENSION(:), POINTER                :: t_surf_green_h_p 
    REAL(wp), DIMENSION(:), POINTER                :: t_surf_10cm_h
    REAL(wp), DIMENSION(:), POINTER                :: t_surf_10cm_h_p

    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_h_1
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_h_2
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_window_h_1
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_window_h_2
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_green_h_1
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_green_h_2
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_10cm_h_1
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_10cm_h_2

    TYPE(t_surf_vertical), DIMENSION(:), POINTER ::  t_surf_v
    TYPE(t_surf_vertical), DIMENSION(:), POINTER ::  t_surf_v_p
    TYPE(t_surf_vertical), DIMENSION(:), POINTER ::  t_surf_window_v
    TYPE(t_surf_vertical), DIMENSION(:), POINTER ::  t_surf_window_v_p
    TYPE(t_surf_vertical), DIMENSION(:), POINTER ::  t_surf_green_v
    TYPE(t_surf_vertical), DIMENSION(:), POINTER ::  t_surf_green_v_p
    TYPE(t_surf_vertical), DIMENSION(:), POINTER ::  t_surf_10cm_v
    TYPE(t_surf_vertical), DIMENSION(:), POINTER ::  t_surf_10cm_v_p

    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  :: t_surf_v_1
    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  :: t_surf_v_2
    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  :: t_surf_window_v_1
    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  :: t_surf_window_v_2
    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  :: t_surf_green_v_1
    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  :: t_surf_green_v_2
    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  :: t_surf_10cm_v_1
    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  :: t_surf_10cm_v_2
    
#endif
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_av          !< average of wall surface temperature (K)
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_window_av   !< average of window surface temperature (K)
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_green_av    !< average of green wall surface temperature (K)
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_10cm_av    !< average of whole wall surface temperature (K)

!-- Temporal tendencies for time stepping            
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: tt_surface_m       !< surface temperature tendency of wall (K)
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: tt_surface_window_m !< surface temperature tendency of window (K)
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: tt_surface_green_m !< surface temperature tendency of green wall (K)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-- Energy balance variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-- parameters of the land, roof and wall surfaces

#if defined( __nopointer )
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET    :: t_wall_h             !< Wall temperature (K)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET    :: t_wall_h_av          !< Average of t_wall
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET    :: t_wall_h_p           !< Prog. wall temperature (K)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET    :: t_window_h           !< Window temperature (K)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET    :: t_window_h_av        !< Average of t_window
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET    :: t_window_h_p         !< Prog. window temperature (K)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET    :: t_green_h            !< Green temperature (K)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET    :: t_green_h_av         !< Average of t_green
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET    :: t_green_h_p          !< Prog. green temperature (K)

    TYPE(t_wall_vertical), DIMENSION(0:3), TARGET  :: t_wall_v             !< Wall temperature (K)
    TYPE(t_wall_vertical), DIMENSION(0:3), TARGET  :: t_wall_v_av          !< Average of t_wall
    TYPE(t_wall_vertical), DIMENSION(0:3), TARGET  :: t_wall_v_p           !< Prog. wall temperature (K)
    TYPE(t_wall_vertical), DIMENSION(0:3), TARGET  :: t_window_v           !< Window temperature (K)
    TYPE(t_wall_vertical), DIMENSION(0:3), TARGET  :: t_window_v_av        !< Average of t_window
    TYPE(t_wall_vertical), DIMENSION(0:3), TARGET  :: t_window_v_p         !< Prog. window temperature (K)
    TYPE(t_wall_vertical), DIMENSION(0:3), TARGET  :: t_green_v            !< Green temperature (K)
    TYPE(t_wall_vertical), DIMENSION(0:3), TARGET  :: t_green_v_av         !< Average of t_green
    TYPE(t_wall_vertical), DIMENSION(0:3), TARGET  :: t_green_v_p          !< Prog. green temperature (K)
#else
    REAL(wp), DIMENSION(:,:), POINTER                :: t_wall_h, t_wall_h_p
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET    :: t_wall_h_av, t_wall_h_1, t_wall_h_2
    REAL(wp), DIMENSION(:,:), POINTER                :: t_window_h, t_window_h_p
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET    :: t_window_h_av, t_window_h_1, t_window_h_2
    REAL(wp), DIMENSION(:,:), POINTER                :: t_green_h, t_green_h_p
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET    :: t_green_h_av, t_green_h_1, t_green_h_2

    TYPE(t_wall_vertical), DIMENSION(:), POINTER   :: t_wall_v, t_wall_v_p
    TYPE(t_wall_vertical), DIMENSION(0:3), TARGET  :: t_wall_v_av, t_wall_v_1, t_wall_v_2
    TYPE(t_wall_vertical), DIMENSION(:), POINTER   :: t_window_v, t_window_v_p
    TYPE(t_wall_vertical), DIMENSION(0:3), TARGET  :: t_window_v_av, t_window_v_1, t_window_v_2
    TYPE(t_wall_vertical), DIMENSION(:), POINTER   :: t_green_v, t_green_v_p
    TYPE(t_wall_vertical), DIMENSION(0:3), TARGET  :: t_green_v_av, t_green_v_1, t_green_v_2
#endif

!-- Wall temporal tendencies for time stepping
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          :: tt_wall_m          !< t_wall prognostic array 
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          :: tt_window_m        !< t_window prognostic array 
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          :: tt_green_m         !< t_green prognostic array 

!-- Surface and material parameters classes (surface_type)
!-- albedo, emissivity, lambda_surf, roughness, thickness, volumetric heat capacity, thermal conductivity
    INTEGER(iwp)                                   :: n_surface_types      !< number of the wall type categories
    INTEGER(iwp), PARAMETER                        :: n_surface_params = 9 !< number of parameters for each type of the wall
    INTEGER(iwp), PARAMETER                        :: ialbedo  = 1         !< albedo of the surface
    INTEGER(iwp), PARAMETER                        :: iemiss   = 2         !< emissivity of the surface
    INTEGER(iwp), PARAMETER                        :: ilambdas = 3         !< heat conductivity lambda S between surface and material ( W m-2 K-1 )
    INTEGER(iwp), PARAMETER                        :: irough   = 4         !< roughness length z0 for movements
    INTEGER(iwp), PARAMETER                        :: iroughh  = 5         !< roughness length z0h for scalars (heat, humidity,...)
    INTEGER(iwp), PARAMETER                        :: icsurf   = 6         !< Surface skin layer heat capacity (J m-2 K-1 )
    INTEGER(iwp), PARAMETER                        :: ithick   = 7         !< thickness of the surface (wall, roof, land)  ( m )
    INTEGER(iwp), PARAMETER                        :: irhoC    = 8         !< volumetric heat capacity rho*C of the material ( J m-3 K-1 )
    INTEGER(iwp), PARAMETER                        :: ilambdah = 9         !< thermal conductivity lambda H of the wall (W m-1 K-1 )
    CHARACTER(12), DIMENSION(:), ALLOCATABLE       :: surface_type_names   !< names of wall types (used only for reports)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE        :: surface_type_codes   !< codes of wall types
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          :: surface_params       !< parameters of wall types

    
!-- interfaces of subroutines accessed from outside of this module
    INTERFACE usm_boundary_condition
       MODULE PROCEDURE usm_boundary_condition
    END INTERFACE usm_boundary_condition

    INTERFACE usm_check_data_output
       MODULE PROCEDURE usm_check_data_output
    END INTERFACE usm_check_data_output
    
    INTERFACE usm_check_parameters
       MODULE PROCEDURE usm_check_parameters
    END INTERFACE usm_check_parameters
    
    INTERFACE usm_data_output_3d
       MODULE PROCEDURE usm_data_output_3d
    END INTERFACE usm_data_output_3d
    
    INTERFACE usm_define_netcdf_grid
       MODULE PROCEDURE usm_define_netcdf_grid
    END INTERFACE usm_define_netcdf_grid

    INTERFACE usm_init_urban_surface
       MODULE PROCEDURE usm_init_urban_surface
    END INTERFACE usm_init_urban_surface

    INTERFACE usm_material_heat_model
       MODULE PROCEDURE usm_material_heat_model
    END INTERFACE usm_material_heat_model
    
    INTERFACE usm_green_heat_model
       MODULE PROCEDURE usm_green_heat_model
    END INTERFACE usm_green_heat_model
    
    INTERFACE usm_parin
       MODULE PROCEDURE usm_parin
    END INTERFACE usm_parin
    
    INTERFACE usm_temperature_near_surface
       MODULE PROCEDURE usm_temperature_near_surface
    END INTERFACE usm_temperature_near_surface

    INTERFACE usm_rrd_local 
       MODULE PROCEDURE usm_rrd_local
    END INTERFACE usm_rrd_local

    INTERFACE usm_surface_energy_balance
       MODULE PROCEDURE usm_surface_energy_balance
    END INTERFACE usm_surface_energy_balance
    
    INTERFACE usm_swap_timelevel
       MODULE PROCEDURE usm_swap_timelevel
    END INTERFACE usm_swap_timelevel
        
    INTERFACE usm_wrd_local
       MODULE PROCEDURE usm_wrd_local
    END INTERFACE usm_wrd_local

    INTERFACE usm_allocate_surface
       MODULE PROCEDURE usm_allocate_surface
    END INTERFACE usm_allocate_surface

    INTERFACE usm_average_3d_data
       MODULE PROCEDURE usm_average_3d_data
    END INTERFACE usm_average_3d_data

    
    SAVE

    PRIVATE 
    
!-- Public functions
    PUBLIC usm_boundary_condition, usm_check_parameters, usm_init_urban_surface,&
           usm_rrd_local,                                                      & 
           usm_surface_energy_balance, usm_material_heat_model,                &
           usm_swap_timelevel, usm_check_data_output, usm_average_3d_data,     &
           usm_data_output_3d, usm_define_netcdf_grid, usm_parin,              &
           usm_wrd_local, usm_allocate_surface

!-- Public parameters, constants and initial values
    PUBLIC usm_anthropogenic_heat, usm_material_model,                          &
           usm_green_heat_model, usm_temperature_near_surface



 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine creates the necessary indices of the urban surfaces
!> and plant canopy and it allocates the needed arrays for USM
!------------------------------------------------------------------------------!
    SUBROUTINE usm_allocate_surface
    
        IMPLICIT NONE
       
        INTEGER(iwp) ::  l

!
!--     Allocate radiation arrays which are part of the new data type. 
!--     For horizontal surfaces.
        ALLOCATE( surf_usm_h%surfhf(1:surf_usm_h%ns)    )
        ALLOCATE( surf_usm_h%rad_net_l(1:surf_usm_h%ns) )
!
!--     For vertical surfaces
        DO  l = 0, 3
           ALLOCATE( surf_usm_v(l)%surfhf(1:surf_usm_v(l)%ns)    )
           ALLOCATE( surf_usm_v(l)%rad_net_l(1:surf_usm_v(l)%ns) )
        ENDDO

!--     Wall surface model
!--     allocate arrays for wall surface model and define pointers
       
!--     allocate array of wall types and wall parameters
        ALLOCATE ( surf_usm_h%surface_types(1:surf_usm_h%ns) )
        DO  l = 0, 3
           ALLOCATE( surf_usm_v(l)%surface_types(1:surf_usm_v(l)%ns) )
        ENDDO
!
!--     Allocate albedo_type and albedo. Each surface element
!--     has 3 values, 0: wall fraction, 1: green fraction, 2: window fraction.
        ALLOCATE( surf_usm_h%albedo_type(0:2,1:surf_usm_h%ns) )
        ALLOCATE( surf_usm_h%albedo(0:2,1:surf_usm_h%ns)      )
        surf_usm_h%albedo_type = albedo_type
        DO  l = 0, 3
           ALLOCATE( surf_usm_v(l)%albedo_type(0:2,1:surf_usm_v(l)%ns) )
           ALLOCATE( surf_usm_v(l)%albedo(0:2,1:surf_usm_v(l)%ns)      )
           surf_usm_v(l)%albedo_type = albedo_type
        ENDDO       


!
!--     Allocate indoor target temperature for summer and winter
        ALLOCATE( surf_usm_h%target_temp_summer(1:surf_usm_h%ns) )
        ALLOCATE( surf_usm_h%target_temp_winter(1:surf_usm_h%ns) )
        DO  l = 0, 3
           ALLOCATE( surf_usm_v(l)%target_temp_summer(1:surf_usm_v(l)%ns) )
           ALLOCATE( surf_usm_v(l)%target_temp_winter(1:surf_usm_v(l)%ns) )
        ENDDO   
!
!--     Allocate flag indicating ground floor level surface elements
        ALLOCATE ( surf_usm_h%ground_level(1:surf_usm_h%ns) ) 
        DO  l = 0, 3
           ALLOCATE( surf_usm_v(l)%ground_level(1:surf_usm_v(l)%ns) )
        ENDDO   
!
!--      Allocate arrays for relative surface fraction. 
!--      0 - wall fraction, 1 - green fraction, 2 - window fraction
         ALLOCATE( surf_usm_h%frac(0:2,1:surf_usm_h%ns) )
         surf_usm_h%frac = 0.0_wp
         DO  l = 0, 3
            ALLOCATE( surf_usm_v(l)%frac(0:2,1:surf_usm_v(l)%ns) )
            surf_usm_v(l)%frac = 0.0_wp
         ENDDO
       
!--     wall and roof surface parameters. First for horizontal surfaces
        ALLOCATE ( surf_usm_h%isroof_surf(1:surf_usm_h%ns)     )
        ALLOCATE ( surf_usm_h%lambda_surf(1:surf_usm_h%ns)     )
        ALLOCATE ( surf_usm_h%lambda_surf_window(1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%lambda_surf_green(1:surf_usm_h%ns)  )
        ALLOCATE ( surf_usm_h%c_surface(1:surf_usm_h%ns)       )
        ALLOCATE ( surf_usm_h%c_surface_window(1:surf_usm_h%ns)   )
        ALLOCATE ( surf_usm_h%c_surface_green(1:surf_usm_h%ns)    )
        ALLOCATE ( surf_usm_h%transmissivity(1:surf_usm_h%ns)  )
        ALLOCATE ( surf_usm_h%lai(1:surf_usm_h%ns)             )
        ALLOCATE ( surf_usm_h%emissivity(0:2,1:surf_usm_h%ns)  )
        ALLOCATE ( surf_usm_h%r_a(1:surf_usm_h%ns)             )
        ALLOCATE ( surf_usm_h%r_a_green(1:surf_usm_h%ns)       )
        ALLOCATE ( surf_usm_h%r_a_window(1:surf_usm_h%ns)      )

!
!--     For vertical surfaces.
        DO  l = 0, 3
           ALLOCATE ( surf_usm_v(l)%lambda_surf(1:surf_usm_v(l)%ns)     )
           ALLOCATE ( surf_usm_v(l)%c_surface(1:surf_usm_v(l)%ns)       )
           ALLOCATE ( surf_usm_v(l)%lambda_surf_window(1:surf_usm_v(l)%ns) )
           ALLOCATE ( surf_usm_v(l)%c_surface_window(1:surf_usm_v(l)%ns)   )
           ALLOCATE ( surf_usm_v(l)%lambda_surf_green(1:surf_usm_v(l)%ns)  )
           ALLOCATE ( surf_usm_v(l)%c_surface_green(1:surf_usm_v(l)%ns)    )
           ALLOCATE ( surf_usm_v(l)%transmissivity(1:surf_usm_v(l)%ns)  )
           ALLOCATE ( surf_usm_v(l)%lai(1:surf_usm_v(l)%ns)             )
           ALLOCATE ( surf_usm_v(l)%emissivity(0:2,1:surf_usm_v(l)%ns)  )
           ALLOCATE ( surf_usm_v(l)%r_a(1:surf_usm_v(l)%ns)             )
           ALLOCATE ( surf_usm_v(l)%r_a_green(1:surf_usm_v(l)%ns)       )
           ALLOCATE ( surf_usm_v(l)%r_a_window(1:surf_usm_v(l)%ns)      )
        ENDDO

!       
!--     allocate wall and roof material parameters. First for horizontal surfaces
        ALLOCATE ( surf_usm_h%thickness_wall(1:surf_usm_h%ns)               )
        ALLOCATE ( surf_usm_h%thickness_window(1:surf_usm_h%ns)                  )
        ALLOCATE ( surf_usm_h%thickness_green(1:surf_usm_h%ns)                   )
        ALLOCATE ( surf_usm_h%lambda_h(nzb_wall:nzt_wall,1:surf_usm_h%ns)   )
        ALLOCATE ( surf_usm_h%rho_c_wall(nzb_wall:nzt_wall,1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%lambda_h_window(nzb_wall:nzt_wall,1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%rho_c_window(nzb_wall:nzt_wall,1:surf_usm_h%ns)    )
        ALLOCATE ( surf_usm_h%lambda_h_green(nzb_wall:nzt_wall,1:surf_usm_h%ns)  )
        ALLOCATE ( surf_usm_h%rho_c_green(nzb_wall:nzt_wall,1:surf_usm_h%ns)     )

!
!--     For vertical surfaces.
        DO  l = 0, 3
           ALLOCATE ( surf_usm_v(l)%thickness_wall(1:surf_usm_v(l)%ns)               )
           ALLOCATE ( surf_usm_v(l)%thickness_window(1:surf_usm_v(l)%ns)                  )
           ALLOCATE ( surf_usm_v(l)%thickness_green(1:surf_usm_v(l)%ns)                   )
           ALLOCATE ( surf_usm_v(l)%lambda_h(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)   )
           ALLOCATE ( surf_usm_v(l)%rho_c_wall(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns) )
           ALLOCATE ( surf_usm_v(l)%lambda_h_window(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns) )
           ALLOCATE ( surf_usm_v(l)%rho_c_window(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)    )
           ALLOCATE ( surf_usm_v(l)%lambda_h_green(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)  )
           ALLOCATE ( surf_usm_v(l)%rho_c_green(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)     )
        ENDDO

!--     allocate wall and roof layers sizes. For horizontal surfaces.
        ALLOCATE ( zwn(nzb_wall:nzt_wall) )
        ALLOCATE ( surf_usm_h%dz_wall(nzb_wall:nzt_wall+1,1:surf_usm_h%ns)     )
        ALLOCATE ( zwn_window(nzb_wall:nzt_wall) )
        ALLOCATE ( surf_usm_h%dz_window(nzb_wall:nzt_wall+1,1:surf_usm_h%ns)     )
        ALLOCATE ( zwn_green(nzb_wall:nzt_wall) )
        ALLOCATE ( surf_usm_h%dz_green(nzb_wall:nzt_wall+1,1:surf_usm_h%ns)      )
        ALLOCATE ( surf_usm_h%ddz_wall(nzb_wall:nzt_wall+1,1:surf_usm_h%ns)    )
        ALLOCATE ( surf_usm_h%dz_wall_stag(nzb_wall:nzt_wall,1:surf_usm_h%ns)  )
        ALLOCATE ( surf_usm_h%ddz_wall_stag(nzb_wall:nzt_wall,1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%zw(nzb_wall:nzt_wall,1:surf_usm_h%ns)            )
        ALLOCATE ( surf_usm_h%ddz_window(nzb_wall:nzt_wall+1,1:surf_usm_h%ns)    )
        ALLOCATE ( surf_usm_h%dz_window_stag(nzb_wall:nzt_wall,1:surf_usm_h%ns)  )
        ALLOCATE ( surf_usm_h%ddz_window_stag(nzb_wall:nzt_wall,1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%zw_window(nzb_wall:nzt_wall,1:surf_usm_h%ns)       )
        ALLOCATE ( surf_usm_h%ddz_green(nzb_wall:nzt_wall+1,1:surf_usm_h%ns)     )
        ALLOCATE ( surf_usm_h%dz_green_stag(nzb_wall:nzt_wall,1:surf_usm_h%ns)   )
        ALLOCATE ( surf_usm_h%ddz_green_stag(nzb_wall:nzt_wall,1:surf_usm_h%ns)  )
        ALLOCATE ( surf_usm_h%zw_green(nzb_wall:nzt_wall,1:surf_usm_h%ns)        )
!
!--     For vertical surfaces.
        DO  l = 0, 3
           ALLOCATE ( surf_usm_v(l)%dz_wall(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)     )
           ALLOCATE ( surf_usm_v(l)%dz_window(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)     )
           ALLOCATE ( surf_usm_v(l)%dz_green(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)      )
           ALLOCATE ( surf_usm_v(l)%ddz_wall(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)    )
           ALLOCATE ( surf_usm_v(l)%dz_wall_stag(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)  )
           ALLOCATE ( surf_usm_v(l)%ddz_wall_stag(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns) )
           ALLOCATE ( surf_usm_v(l)%zw(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)            )
           ALLOCATE ( surf_usm_v(l)%ddz_window(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)    )
           ALLOCATE ( surf_usm_v(l)%dz_window_stag(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)  )
           ALLOCATE ( surf_usm_v(l)%ddz_window_stag(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns) )
           ALLOCATE ( surf_usm_v(l)%zw_window(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)       )
           ALLOCATE ( surf_usm_v(l)%ddz_green(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)     )
           ALLOCATE ( surf_usm_v(l)%dz_green_stag(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)   )
           ALLOCATE ( surf_usm_v(l)%ddz_green_stag(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)  )
           ALLOCATE ( surf_usm_v(l)%zw_green(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)        )
        ENDDO

!--     allocate wall and roof temperature arrays, for horizontal walls
#if defined( __nopointer )
        IF ( .NOT. ALLOCATED( t_surf_h ) )                                     &
           ALLOCATE ( t_surf_h(1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( t_surf_h_p ) )                                   &
           ALLOCATE ( t_surf_h_p(1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( t_wall_h ) )                                     &           
           ALLOCATE ( t_wall_h(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) ) 
        IF ( .NOT. ALLOCATED( t_wall_h_p ) )                                   &           
           ALLOCATE ( t_wall_h_p(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) )  
        IF ( .NOT. ALLOCATED( t_surf_window_h ) )                              &
           ALLOCATE ( t_surf_window_h(1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( t_surf_window_h_p ) )                            &
           ALLOCATE ( t_surf_window_h_p(1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( t_window_h ) )                                   &           
           ALLOCATE ( t_window_h(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) ) 
        IF ( .NOT. ALLOCATED( t_window_h_p ) )                                 &           
           ALLOCATE ( t_window_h_p(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) )  
        IF ( .NOT. ALLOCATED( t_surf_green_h ) )                               &
           ALLOCATE ( t_surf_green_h(1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( t_surf_green_h_p ) )                             &
           ALLOCATE ( t_surf_green_h_p(1:surf_usm_h%ns) )           
        IF ( .NOT. ALLOCATED( t_green_h ) )                                    &           
           ALLOCATE ( t_green_h(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) ) 
        IF ( .NOT. ALLOCATED( t_green_h_p ) )                                  &           
           ALLOCATE ( t_green_h_p(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) )  
        IF ( .NOT. ALLOCATED( t_surf_10cm_h ) )                                &
           ALLOCATE ( t_surf_10cm_h(1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( t_surf_10cm_h_p ) )                              &
           ALLOCATE ( t_surf_10cm_h_p(1:surf_usm_h%ns) )
#else
!
!--     Allocate if required. Note, in case of restarts, some of these arrays 
!--     might be already allocated.
        IF ( .NOT. ALLOCATED( t_surf_h_1 ) )                                   &
           ALLOCATE ( t_surf_h_1(1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( t_surf_h_2 ) )                                   &
           ALLOCATE ( t_surf_h_2(1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( t_wall_h_1 ) )                                   &           
           ALLOCATE ( t_wall_h_1(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) ) 
        IF ( .NOT. ALLOCATED( t_wall_h_2 ) )                                   &           
           ALLOCATE ( t_wall_h_2(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) )          
        IF ( .NOT. ALLOCATED( t_surf_window_h_1 ) )                            &
           ALLOCATE ( t_surf_window_h_1(1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( t_surf_window_h_2 ) )                            &
           ALLOCATE ( t_surf_window_h_2(1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( t_window_h_1 ) )                                 &           
           ALLOCATE ( t_window_h_1(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) ) 
        IF ( .NOT. ALLOCATED( t_window_h_2 ) )                                 &           
           ALLOCATE ( t_window_h_2(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) )          
        IF ( .NOT. ALLOCATED( t_surf_green_h_1 ) )                             &
           ALLOCATE ( t_surf_green_h_1(1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( t_surf_green_h_2 ) )                             &
           ALLOCATE ( t_surf_green_h_2(1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( t_green_h_1 ) )                                  &           
           ALLOCATE ( t_green_h_1(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) ) 
        IF ( .NOT. ALLOCATED( t_green_h_2 ) )                                  &           
           ALLOCATE ( t_green_h_2(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) )          
        IF ( .NOT. ALLOCATED( t_surf_10cm_h_1 ) )                              &
           ALLOCATE ( t_surf_10cm_h_1(1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( t_surf_10cm_h_2 ) )                              &
           ALLOCATE ( t_surf_10cm_h_2(1:surf_usm_h%ns) )
!           
!--     initial assignment of the pointers
        t_wall_h    => t_wall_h_1;    t_wall_h_p    => t_wall_h_2
        t_window_h    => t_window_h_1;    t_window_h_p    => t_window_h_2
        t_green_h    => t_green_h_1;    t_green_h_p    => t_green_h_2
        t_surf_h => t_surf_h_1; t_surf_h_p => t_surf_h_2           
        t_surf_window_h => t_surf_window_h_1; t_surf_window_h_p => t_surf_window_h_2  
        t_surf_green_h => t_surf_green_h_1; t_surf_green_h_p => t_surf_green_h_2           
        t_surf_10cm_h => t_surf_10cm_h_1; t_surf_10cm_h_p => t_surf_10cm_h_2  
 
#endif

!--     allocate wall and roof temperature arrays, for vertical walls if required
#if defined( __nopointer )
        DO  l = 0, 3
           IF ( .NOT. ALLOCATED( t_surf_v(l)%t ) )                             &
              ALLOCATE ( t_surf_v(l)%t(1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_surf_v_p(l)%t ) )                           &
              ALLOCATE ( t_surf_v_p(l)%t(1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_wall_v(l)%t ) )                             &
              ALLOCATE ( t_wall_v(l)%t(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_wall_v_p(l)%t ) )                           &                 
              ALLOCATE ( t_wall_v_p(l)%t(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_surf_window_v(l)%t ) )                      &
              ALLOCATE ( t_surf_window_v(l)%t(1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_surf_window_v_p(l)%t ) )                    &
              ALLOCATE ( t_surf_window_v_p(l)%t(1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_window_v(l)%t ) )                           &
              ALLOCATE ( t_window_v(l)%t(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_window_v_p(l)%t ) )                         &                 
              ALLOCATE ( t_window_v_p(l)%t(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_green_v(l)%t ) )                            &
              ALLOCATE ( t_green_v(l)%t(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_green_v_p(l)%t ) )                          &                 
              ALLOCATE ( t_green_v_p(l)%t(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_surf_green_v(l)%t ) )                       &
              ALLOCATE ( t_surf_green_v(l)%t(1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_surf_green_v_p(l)%t ) )                     &
              ALLOCATE ( t_surf_green_v_p(l)%t(1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_surf_10cm_v(l)%t ) )                        &
              ALLOCATE ( t_surf_10cm_v(l)%t(1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_surf_10cm_v_p(l)%t ) )                        &
              ALLOCATE ( t_surf_10cm_v_p(l)%t(1:surf_usm_v(l)%ns) )
        ENDDO
#else
!
!--     Allocate if required. Note, in case of restarts, some of these arrays 
!--     might be already allocated.
        DO  l = 0, 3
           IF ( .NOT. ALLOCATED( t_surf_v_1(l)%t ) )                           &
              ALLOCATE ( t_surf_v_1(l)%t(1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_surf_v_2(l)%t ) )                           &
              ALLOCATE ( t_surf_v_2(l)%t(1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_wall_v_1(l)%t ) )                           &           
              ALLOCATE ( t_wall_v_1(l)%t(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) ) 
           IF ( .NOT. ALLOCATED( t_wall_v_2(l)%t ) )                           &           
              ALLOCATE ( t_wall_v_2(l)%t(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )  
           IF ( .NOT. ALLOCATED( t_surf_window_v_1(l)%t ) )                    &
              ALLOCATE ( t_surf_window_v_1(l)%t(1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_surf_window_v_2(l)%t ) )                    &
              ALLOCATE ( t_surf_window_v_2(l)%t(1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_window_v_1(l)%t ) )                         &           
              ALLOCATE ( t_window_v_1(l)%t(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) ) 
           IF ( .NOT. ALLOCATED( t_window_v_2(l)%t ) )                         &           
              ALLOCATE ( t_window_v_2(l)%t(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )  
           IF ( .NOT. ALLOCATED( t_surf_green_v_1(l)%t ) )                     &
              ALLOCATE ( t_surf_green_v_1(l)%t(1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_surf_green_v_2(l)%t ) )                     &
              ALLOCATE ( t_surf_green_v_2(l)%t(1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_green_v_1(l)%t ) )                          &           
              ALLOCATE ( t_green_v_1(l)%t(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) ) 
           IF ( .NOT. ALLOCATED( t_green_v_2(l)%t ) )                          &           
              ALLOCATE ( t_green_v_2(l)%t(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )  
           IF ( .NOT. ALLOCATED( t_surf_10cm_v_1(l)%t ) )                     &
              ALLOCATE ( t_surf_10cm_v_1(l)%t(1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_surf_10cm_v_2(l)%t ) )                     &
              ALLOCATE ( t_surf_10cm_v_2(l)%t(1:surf_usm_v(l)%ns) )
        ENDDO
!
!--     initial assignment of the pointers
        t_wall_v    => t_wall_v_1;    t_wall_v_p    => t_wall_v_2
        t_surf_v => t_surf_v_1; t_surf_v_p => t_surf_v_2
        t_window_v    => t_window_v_1;    t_window_v_p    => t_window_v_2
        t_green_v    => t_green_v_1;    t_green_v_p    => t_green_v_2
        t_surf_window_v => t_surf_window_v_1; t_surf_window_v_p => t_surf_window_v_2
        t_surf_green_v => t_surf_green_v_1; t_surf_green_v_p => t_surf_green_v_2
        t_surf_10cm_v => t_surf_10cm_v_1; t_surf_10cm_v_p => t_surf_10cm_v_2

#endif
!
!--     Allocate intermediate timestep arrays. For horizontal surfaces.
        ALLOCATE ( surf_usm_h%tt_surface_m(1:surf_usm_h%ns)                  )
        ALLOCATE ( surf_usm_h%tt_wall_m(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%tt_surface_window_m(1:surf_usm_h%ns)             )
        ALLOCATE ( surf_usm_h%tt_window_m(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%tt_green_m(nzb_wall:nzt_wall+1,1:surf_usm_h%ns)  )
        ALLOCATE ( surf_usm_h%tt_surface_green_m(1:surf_usm_h%ns)              )

!
!--     Set inital values for prognostic quantities
        IF ( ALLOCATED( surf_usm_h%tt_surface_m ) )  surf_usm_h%tt_surface_m = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%tt_wall_m    ) )  surf_usm_h%tt_wall_m    = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%tt_surface_window_m ) )  surf_usm_h%tt_surface_window_m = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%tt_window_m    )      )  surf_usm_h%tt_window_m         = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%tt_green_m    )       )  surf_usm_h%tt_green_m          = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%tt_surface_green_m )  )  surf_usm_h%tt_surface_green_m  = 0.0_wp
!
!--     Now, for vertical surfaces
        DO  l = 0, 3
           ALLOCATE ( surf_usm_v(l)%tt_surface_m(1:surf_usm_v(l)%ns)                  )
           ALLOCATE ( surf_usm_v(l)%tt_wall_m(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
           IF ( ALLOCATED( surf_usm_v(l)%tt_surface_m ) )  surf_usm_v(l)%tt_surface_m = 0.0_wp
           IF ( ALLOCATED( surf_usm_v(l)%tt_wall_m    ) )  surf_usm_v(l)%tt_wall_m    = 0.0_wp
           ALLOCATE ( surf_usm_v(l)%tt_surface_window_m(1:surf_usm_v(l)%ns)             )
           ALLOCATE ( surf_usm_v(l)%tt_window_m(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
           IF ( ALLOCATED( surf_usm_v(l)%tt_surface_window_m ) )  surf_usm_v(l)%tt_surface_window_m = 0.0_wp
           IF ( ALLOCATED( surf_usm_v(l)%tt_window_m  ) )  surf_usm_v(l)%tt_window_m    = 0.0_wp
           ALLOCATE ( surf_usm_v(l)%tt_surface_green_m(1:surf_usm_v(l)%ns)              )
           IF ( ALLOCATED( surf_usm_v(l)%tt_surface_green_m ) )  surf_usm_v(l)%tt_surface_green_m = 0.0_wp
           ALLOCATE ( surf_usm_v(l)%tt_green_m(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)  )
           IF ( ALLOCATED( surf_usm_v(l)%tt_green_m   ) )  surf_usm_v(l)%tt_green_m    = 0.0_wp
        ENDDO

!--     allocate wall heat flux output array and set initial values. For horizontal surfaces
!         ALLOCATE ( surf_usm_h%wshf(1:surf_usm_h%ns)    )  !can be removed 
        ALLOCATE ( surf_usm_h%wshf_eb(1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%wghf_eb(1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%wghf_eb_window(1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%wghf_eb_green(1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%iwghf_eb(1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%iwghf_eb_window(1:surf_usm_h%ns) )
        IF ( ALLOCATED( surf_usm_h%wshf    ) )  surf_usm_h%wshf    = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%wshf_eb ) )  surf_usm_h%wshf_eb = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%wghf_eb ) )  surf_usm_h%wghf_eb = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%wghf_eb_window ) )  surf_usm_h%wghf_eb_window = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%wghf_eb_green ) )  surf_usm_h%wghf_eb_green = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%iwghf_eb ) )  surf_usm_h%iwghf_eb = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%iwghf_eb_window ) )  surf_usm_h%iwghf_eb_window = 0.0_wp
!
!--     Now, for vertical surfaces
        DO  l = 0, 3
!            ALLOCATE ( surf_usm_v(l)%wshf(1:surf_usm_v(l)%ns)    )    ! can be removed
           ALLOCATE ( surf_usm_v(l)%wshf_eb(1:surf_usm_v(l)%ns) )
           ALLOCATE ( surf_usm_v(l)%wghf_eb(1:surf_usm_v(l)%ns) )
           ALLOCATE ( surf_usm_v(l)%wghf_eb_window(1:surf_usm_v(l)%ns) )
           ALLOCATE ( surf_usm_v(l)%wghf_eb_green(1:surf_usm_v(l)%ns) )
           ALLOCATE ( surf_usm_v(l)%iwghf_eb(1:surf_usm_v(l)%ns) )
           ALLOCATE ( surf_usm_v(l)%iwghf_eb_window(1:surf_usm_v(l)%ns) )
           IF ( ALLOCATED( surf_usm_v(l)%wshf    ) )  surf_usm_v(l)%wshf    = 0.0_wp
           IF ( ALLOCATED( surf_usm_v(l)%wshf_eb ) )  surf_usm_v(l)%wshf_eb = 0.0_wp
           IF ( ALLOCATED( surf_usm_v(l)%wghf_eb ) )  surf_usm_v(l)%wghf_eb = 0.0_wp
           IF ( ALLOCATED( surf_usm_v(l)%wghf_eb_window ) )  surf_usm_v(l)%wghf_eb_window = 0.0_wp
           IF ( ALLOCATED( surf_usm_v(l)%wghf_eb_green ) )  surf_usm_v(l)%wghf_eb_green = 0.0_wp
           IF ( ALLOCATED( surf_usm_v(l)%iwghf_eb ) )  surf_usm_v(l)%iwghf_eb = 0.0_wp
           IF ( ALLOCATED( surf_usm_v(l)%iwghf_eb_window ) )  surf_usm_v(l)%iwghf_eb_window = 0.0_wp
        ENDDO
        
    END SUBROUTINE usm_allocate_surface


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sum up and time-average urban surface output quantities as well as allocate
!> the array necessary for storing the average.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_average_3d_data( mode, variable )

        IMPLICIT NONE

        CHARACTER (len=*), INTENT(IN) ::  mode
        CHARACTER (len=*), INTENT(IN) :: variable
  
        INTEGER(iwp)                                       :: i, j, k, l, m, ids, idsint, iwl, istat
        CHARACTER (len=varnamelength)                      :: var, surfid
        INTEGER(iwp), PARAMETER                            :: nd = 5
        CHARACTER(len=6), DIMENSION(0:nd-1), PARAMETER     :: dirname = (/ '_roof ', '_south', '_north', '_west ', '_east ' /)
        INTEGER(iwp), DIMENSION(0:nd-1), PARAMETER         :: dirint = (/ iup_u, isouth_u, inorth_u, iwest_u, ieast_u /)

!--     find the real name of the variable
        ids = -1
        var = TRIM(variable)
        DO i = 0, nd-1
            k = len(TRIM(var))
            j = len(TRIM(dirname(i)))
            IF ( var(k-j+1:k) == dirname(i) )  THEN
                ids = i
                idsint = dirint(ids)
                var = var(:k-j)
                EXIT
            ENDIF
        ENDDO
        IF ( ids == -1 )  THEN
            var = TRIM(variable)
        ENDIF
        IF ( var(1:11) == 'usm_t_wall_'  .AND.  len(TRIM(var)) >= 12 )  THEN
!--          wall layers
            READ(var(12:12), '(I1)', iostat=istat ) iwl
            IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
                var = var(1:10)
            ELSE
!--             wrong wall layer index
                RETURN
            ENDIF
        ENDIF
        IF ( var(1:13) == 'usm_t_window_'  .AND.  len(TRIM(var)) >= 14 )  THEN
!--          wall layers
            READ(var(14:14), '(I1)', iostat=istat ) iwl
            IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
                var = var(1:12)
            ELSE
!--             wrong window layer index
                RETURN
            ENDIF
        ENDIF
        IF ( var(1:12) == 'usm_t_green_'  .AND.  len(TRIM(var)) >= 13 )  THEN
!--          wall layers
            READ(var(13:13), '(I1)', iostat=istat ) iwl
            IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
                var = var(1:11)
            ELSE
!--             wrong green layer index
                RETURN
            ENDIF
        ENDIF

        IF ( mode == 'allocate' )  THEN
           
           SELECT CASE ( TRIM( var ) )
                
                CASE ( 'usm_rad_net' )
!--                 array of complete radiation balance
                    IF ( .NOT.  ALLOCATED(surf_usm_h%rad_net_av) )  THEN
                        ALLOCATE( surf_usm_h%rad_net_av(1:surf_usm_h%ns) )
                        surf_usm_h%rad_net_av = 0.0_wp
                    ENDIF
                    DO  l = 0, 3
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%rad_net_av) )  THEN
                           ALLOCATE( surf_usm_v(l)%rad_net_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%rad_net_av = 0.0_wp
                       ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_insw' )
!--                 array of sw radiation falling to surface after i-th reflection
                    IF ( .NOT.  ALLOCATED(surfinsw_av) )  THEN
                        ALLOCATE( surfinsw_av(nsurfl) )
                        surfinsw_av = 0.0_wp
                    ENDIF

                CASE ( 'usm_rad_inlw' )
!--                 array of lw radiation falling to surface after i-th reflection
                    IF ( .NOT.  ALLOCATED(surfinlw_av) )  THEN
                        ALLOCATE( surfinlw_av(nsurfl) )
                        surfinlw_av = 0.0_wp
                    ENDIF

                CASE ( 'usm_rad_inswdir' )
!--                 array of direct sw radiation falling to surface from sun
                    IF ( .NOT.  ALLOCATED(surfinswdir_av) )  THEN
                        ALLOCATE( surfinswdir_av(nsurfl) )
                        surfinswdir_av = 0.0_wp
                    ENDIF

                CASE ( 'usm_rad_inswdif' )
!--                 array of difusion sw radiation falling to surface from sky and borders of the domain
                    IF ( .NOT.  ALLOCATED(surfinswdif_av) )  THEN
                        ALLOCATE( surfinswdif_av(nsurfl) )
                        surfinswdif_av = 0.0_wp
                    ENDIF

                CASE ( 'usm_rad_inswref' )
!--                 array of sw radiation falling to surface from reflections
                    IF ( .NOT.  ALLOCATED(surfinswref_av) )  THEN
                        ALLOCATE( surfinswref_av(nsurfl) )
                        surfinswref_av = 0.0_wp
                    ENDIF

                CASE ( 'usm_rad_inlwdif' )
!--                 array of sw radiation falling to surface after i-th reflection
                   IF ( .NOT.  ALLOCATED(surfinlwdif_av) )  THEN
                        ALLOCATE( surfinlwdif_av(nsurfl) )
                        surfinlwdif_av = 0.0_wp
                    ENDIF

                CASE ( 'usm_rad_inlwref' )
!--                 array of lw radiation falling to surface from reflections
                    IF ( .NOT.  ALLOCATED(surfinlwref_av) )  THEN
                        ALLOCATE( surfinlwref_av(nsurfl) )
                        surfinlwref_av = 0.0_wp
                    ENDIF

                CASE ( 'usm_rad_outsw' )
!--                 array of sw radiation emitted from surface after i-th reflection
                    IF ( .NOT.  ALLOCATED(surfoutsw_av) )  THEN
                        ALLOCATE( surfoutsw_av(nsurfl) )
                        surfoutsw_av = 0.0_wp
                    ENDIF

                CASE ( 'usm_rad_outlw' )
!--                 array of lw radiation emitted from surface after i-th reflection
                    IF ( .NOT.  ALLOCATED(surfoutlw_av) )  THEN
                        ALLOCATE( surfoutlw_av(nsurfl) )
                        surfoutlw_av = 0.0_wp
                    ENDIF
                CASE ( 'usm_rad_ressw' )
!--                 array of residua of sw radiation absorbed in surface after last reflection
                    IF ( .NOT.  ALLOCATED(surfins_av) )  THEN
                        ALLOCATE( surfins_av(nsurfl) )
                        surfins_av = 0.0_wp
                    ENDIF
                                   
                CASE ( 'usm_rad_reslw' )
!--                 array of residua of lw radiation absorbed in surface after last reflection
                    IF ( .NOT.  ALLOCATED(surfinl_av) )  THEN
                        ALLOCATE( surfinl_av(nsurfl) )
                        surfinl_av = 0.0_wp
                    ENDIF
                                    
                CASE ( 'usm_rad_hf' )
!--                 array of heat flux from radiation for surfaces after i-th reflection
                    IF ( .NOT.  ALLOCATED(surf_usm_h%surfhf_av) )  THEN
                        ALLOCATE( surf_usm_h%surfhf_av(1:surf_usm_h%ns) )
                        surf_usm_h%surfhf_av = 0.0_wp
                    ENDIF
                    DO  l = 0, 3
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%surfhf_av) )  THEN
                           ALLOCATE( surf_usm_v(l)%surfhf_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%surfhf_av = 0.0_wp
                       ENDIF
                    ENDDO

                CASE ( 'usm_wshf' )
!--                 array of sensible heat flux from surfaces
!--                 land surfaces
                    IF ( .NOT.  ALLOCATED(surf_usm_h%wshf_eb_av) )  THEN
                        ALLOCATE( surf_usm_h%wshf_eb_av(1:surf_usm_h%ns) )
                        surf_usm_h%wshf_eb_av = 0.0_wp
                    ENDIF
                    DO  l = 0, 3
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%wshf_eb_av) )  THEN
                           ALLOCATE( surf_usm_v(l)%wshf_eb_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%wshf_eb_av = 0.0_wp
                       ENDIF
                    ENDDO
!
!--             Please note, the following output quantities belongs to the 
!--             individual tile fractions - ground heat flux at wall-, window-, 
!--             and green fraction. Aggregated ground-heat flux is treated
!--             accordingly in average_3d_data, sum_up_3d_data, etc..
                CASE ( 'usm_wghf' )
!--                 array of heat flux from ground (wall, roof, land)
                    IF ( .NOT.  ALLOCATED(surf_usm_h%wghf_eb_av) )  THEN
                        ALLOCATE( surf_usm_h%wghf_eb_av(1:surf_usm_h%ns) )
                        surf_usm_h%wghf_eb_av = 0.0_wp
                    ENDIF
                    DO  l = 0, 3
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%wghf_eb_av) )  THEN
                           ALLOCATE( surf_usm_v(l)%wghf_eb_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%wghf_eb_av = 0.0_wp
                       ENDIF
                    ENDDO

                CASE ( 'usm_wghf_window' )
!--                 array of heat flux from window ground (wall, roof, land)
                    IF ( .NOT.  ALLOCATED(surf_usm_h%wghf_eb_window_av) )  THEN
                        ALLOCATE( surf_usm_h%wghf_eb_window_av(1:surf_usm_h%ns) )
                        surf_usm_h%wghf_eb_window_av = 0.0_wp
                    ENDIF
                    DO  l = 0, 3
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%wghf_eb_window_av) )  THEN
                           ALLOCATE( surf_usm_v(l)%wghf_eb_window_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%wghf_eb_window_av = 0.0_wp
                       ENDIF
                    ENDDO

                CASE ( 'usm_wghf_green' )
!--                 array of heat flux from green ground (wall, roof, land)
                    IF ( .NOT.  ALLOCATED(surf_usm_h%wghf_eb_green_av) )  THEN
                        ALLOCATE( surf_usm_h%wghf_eb_green_av(1:surf_usm_h%ns) )
                        surf_usm_h%wghf_eb_green_av = 0.0_wp
                    ENDIF
                    DO  l = 0, 3
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%wghf_eb_green_av) )  THEN
                           ALLOCATE( surf_usm_v(l)%wghf_eb_green_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%wghf_eb_green_av = 0.0_wp
                       ENDIF
                    ENDDO

                CASE ( 'usm_iwghf' )
!--                 array of heat flux from indoor ground (wall, roof, land)
                    IF ( .NOT.  ALLOCATED(surf_usm_h%iwghf_eb_av) )  THEN
                        ALLOCATE( surf_usm_h%iwghf_eb_av(1:surf_usm_h%ns) )
                        surf_usm_h%iwghf_eb_av = 0.0_wp
                    ENDIF
                    DO  l = 0, 3
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%iwghf_eb_av) )  THEN
                           ALLOCATE( surf_usm_v(l)%iwghf_eb_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%iwghf_eb_av = 0.0_wp
                       ENDIF
                    ENDDO

                CASE ( 'usm_iwghf_window' )
!--                 array of heat flux from indoor window ground (wall, roof, land)
                    IF ( .NOT.  ALLOCATED(surf_usm_h%iwghf_eb_window_av) )  THEN
                        ALLOCATE( surf_usm_h%iwghf_eb_window_av(1:surf_usm_h%ns) )
                        surf_usm_h%iwghf_eb_window_av = 0.0_wp
                    ENDIF
                    DO  l = 0, 3
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%iwghf_eb_window_av) )  THEN
                           ALLOCATE( surf_usm_v(l)%iwghf_eb_window_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%iwghf_eb_window_av = 0.0_wp
                       ENDIF
                    ENDDO
                    
                CASE ( 'usm_t_surf' )
!--                 surface temperature for surfaces
                    IF ( .NOT.  ALLOCATED(surf_usm_h%t_surf_av) )  THEN
                        ALLOCATE( surf_usm_h%t_surf_av(1:surf_usm_h%ns) )
                        surf_usm_h%t_surf_av = 0.0_wp
                    ENDIF
                    DO  l = 0, 3
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%t_surf_av) )  THEN
                           ALLOCATE( surf_usm_v(l)%t_surf_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%t_surf_av = 0.0_wp
                       ENDIF
                    ENDDO

                CASE ( 'usm_t_surf_window' )
!--                 surface temperature for window surfaces
                    IF ( .NOT.  ALLOCATED(surf_usm_h%t_surf_window_av) )  THEN
                        ALLOCATE( surf_usm_h%t_surf_window_av(1:surf_usm_h%ns) )
                        surf_usm_h%t_surf_window_av = 0.0_wp
                    ENDIF
                    DO  l = 0, 3
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%t_surf_window_av) )  THEN
                           ALLOCATE( surf_usm_v(l)%t_surf_window_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%t_surf_window_av = 0.0_wp
                       ENDIF
                    ENDDO
                    
                CASE ( 'usm_t_surf_green' )
!--                 surface temperature for green surfaces
                    IF ( .NOT.  ALLOCATED(surf_usm_h%t_surf_green_av) )  THEN
                        ALLOCATE( surf_usm_h%t_surf_green_av(1:surf_usm_h%ns) )
                        surf_usm_h%t_surf_green_av = 0.0_wp
                    ENDIF
                    DO  l = 0, 3
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%t_surf_green_av) )  THEN
                           ALLOCATE( surf_usm_v(l)%t_surf_green_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%t_surf_green_av = 0.0_wp
                       ENDIF
                    ENDDO
                
                CASE ( 'usm_t_surf_10cm' )
!--                 near surface temperature for whole surfaces
                    IF ( .NOT.  ALLOCATED(surf_usm_h%t_surf_10cm_av) )  THEN
                        ALLOCATE( surf_usm_h%t_surf_10cm_av(1:surf_usm_h%ns) )
                        surf_usm_h%t_surf_10cm_av = 0.0_wp
                    ENDIF
                    DO  l = 0, 3
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%t_surf_10cm_av) )  THEN
                           ALLOCATE( surf_usm_v(l)%t_surf_10cm_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%t_surf_10cm_av = 0.0_wp
                       ENDIF
                    ENDDO

                CASE ( 'usm_t_wall' )
!--                 wall temperature for iwl layer of walls and land
                    IF ( .NOT.  ALLOCATED(surf_usm_h%t_wall_av) )  THEN
                        ALLOCATE( surf_usm_h%t_wall_av(nzb_wall:nzt_wall,1:surf_usm_h%ns) )
                        surf_usm_h%t_wall_av = 0.0_wp
                    ENDIF
                    DO  l = 0, 3
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%t_wall_av) )  THEN
                           ALLOCATE( surf_usm_v(l)%t_wall_av(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%t_wall_av = 0.0_wp
                       ENDIF
                    ENDDO

                CASE ( 'usm_t_window' )
!--                 window temperature for iwl layer of walls and land
                    IF ( .NOT.  ALLOCATED(surf_usm_h%t_window_av) )  THEN
                        ALLOCATE( surf_usm_h%t_window_av(nzb_wall:nzt_wall,1:surf_usm_h%ns) )
                        surf_usm_h%t_window_av = 0.0_wp
                    ENDIF
                    DO  l = 0, 3
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%t_window_av) )  THEN
                           ALLOCATE( surf_usm_v(l)%t_window_av(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%t_window_av = 0.0_wp
                       ENDIF
                    ENDDO

                CASE ( 'usm_t_green' )
!--                 green temperature for iwl layer of walls and land
                    IF ( .NOT.  ALLOCATED(surf_usm_h%t_green_av) )  THEN
                        ALLOCATE( surf_usm_h%t_green_av(nzb_wall:nzt_wall,1:surf_usm_h%ns) )
                        surf_usm_h%t_green_av = 0.0_wp
                    ENDIF
                    DO  l = 0, 3
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%t_green_av) )  THEN
                           ALLOCATE( surf_usm_v(l)%t_green_av(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%t_green_av = 0.0_wp
                       ENDIF
                    ENDDO

               CASE DEFAULT
                   CONTINUE

           END SELECT

        ELSEIF ( mode == 'sum' )  THEN
           
           SELECT CASE ( TRIM( var ) )
                
                CASE ( 'usm_rad_net' )
!--                 array of complete radiation balance
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%rad_net_av(m) =                              &
                                          surf_usm_h%rad_net_av(m) +           &
                                          surf_usm_h%rad_net_l(m)
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%rad_net_av(m) =                        &
                                          surf_usm_v(l)%rad_net_av(m) +        &
                                          surf_usm_v(l)%rad_net_l(m)
                       ENDDO
                    ENDDO
                    
                CASE ( 'usm_rad_insw' )
!--                 array of sw radiation falling to surface after i-th reflection
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfinsw_av(l) = surfinsw_av(l) + surfinsw(l)
                        ENDIF
                    ENDDO
                             
                CASE ( 'usm_rad_inlw' )
!--                 array of lw radiation falling to surface after i-th reflection
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfinlw_av(l) = surfinlw_av(l) + surfinlw(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_inswdir' )
!--                 array of direct sw radiation falling to surface from sun
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfinswdir_av(l) = surfinswdir_av(l) + surfinswdir(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_inswdif' )
!--                 array of difusion sw radiation falling to surface from sky and borders of the domain
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfinswdif_av(l) = surfinswdif_av(l) + surfinswdif(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_inswref' )
!--                 array of sw radiation falling to surface from reflections
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfinswref_av(l) = surfinswref_av(l) + surfinsw(l) - &
                                                surfinswdir(l) - surfinswdif(l)
                        ENDIF
                    ENDDO

                    
                CASE ( 'usm_rad_inlwdif' )
!--                 array of sw radiation falling to surface after i-th reflection
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfinlwdif_av(l) = surfinlwdif_av(l) + surfinlwdif(l)
                        ENDIF
                    ENDDO
!                     
                CASE ( 'usm_rad_inlwref' )
!--                 array of lw radiation falling to surface from reflections
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfinlwref_av(l) = surfinlwref_av(l) + &
                                                surfinlw(l) - surfinlwdif(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_outsw' )
!--                 array of sw radiation emitted from surface after i-th reflection
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfoutsw_av(l) = surfoutsw_av(l) + surfoutsw(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_outlw' )
!--                 array of lw radiation emitted from surface after i-th reflection
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfoutlw_av(l) = surfoutlw_av(l) + surfoutlw(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_ressw' )
!--                 array of residua of sw radiation absorbed in surface after last reflection
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfins_av(l) = surfins_av(l) + surfins(l)
                        ENDIF
                    ENDDO
                                    
                CASE ( 'usm_rad_reslw' )
!--                 array of residua of lw radiation absorbed in surface after last reflection
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfinl_av(l) = surfinl_av(l) + surfinl(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_hf' )
!--                 array of heat flux from radiation for surfaces after i-th reflection
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%surfhf_av(m) =                               &
                                          surf_usm_h%surfhf_av(m) +            &
                                          surf_usm_h%surfhf(m)
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%surfhf_av(m) =                         &
                                          surf_usm_v(l)%surfhf_av(m) +         &
                                          surf_usm_v(l)%surfhf(m)
                       ENDDO
                    ENDDO
                    
                CASE ( 'usm_wshf' )
!--                 array of sensible heat flux from surfaces (land, roof, wall)
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%wshf_eb_av(m) =                              &
                                          surf_usm_h%wshf_eb_av(m) +           &
                                          surf_usm_h%wshf_eb(m)
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%wshf_eb_av(m) =                        &
                                          surf_usm_v(l)%wshf_eb_av(m) +        &
                                          surf_usm_v(l)%wshf_eb(m)
                       ENDDO
                    ENDDO
                    
                CASE ( 'usm_wghf' )
!--                 array of heat flux from ground (wall, roof, land)
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%wghf_eb_av(m) =                              &
                                          surf_usm_h%wghf_eb_av(m) +           &
                                          surf_usm_h%wghf_eb(m)
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%wghf_eb_av(m) =                        &
                                          surf_usm_v(l)%wghf_eb_av(m) +        &
                                          surf_usm_v(l)%wghf_eb(m)
                       ENDDO
                    ENDDO
                    
                CASE ( 'usm_wghf_window' )
!--                 array of heat flux from window ground (wall, roof, land)
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%wghf_eb_window_av(m) =                              &
                                          surf_usm_h%wghf_eb_window_av(m) +           &
                                          surf_usm_h%wghf_eb_window(m)
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%wghf_eb_window_av(m) =                        &
                                          surf_usm_v(l)%wghf_eb_window_av(m) +        &
                                          surf_usm_v(l)%wghf_eb_window(m)
                       ENDDO
                    ENDDO

                CASE ( 'usm_wghf_green' )
!--                 array of heat flux from green ground (wall, roof, land)
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%wghf_eb_green_av(m) =                              &
                                          surf_usm_h%wghf_eb_green_av(m) +           &
                                          surf_usm_h%wghf_eb_green(m)
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%wghf_eb_green_av(m) =                        &
                                          surf_usm_v(l)%wghf_eb_green_av(m) +        &
                                          surf_usm_v(l)%wghf_eb_green(m)
                       ENDDO
                    ENDDO
                    
                CASE ( 'usm_iwghf' )
!--                 array of heat flux from indoor ground (wall, roof, land)
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%iwghf_eb_av(m) =                              &
                                          surf_usm_h%iwghf_eb_av(m) +           &
                                          surf_usm_h%iwghf_eb(m)
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%iwghf_eb_av(m) =                        &
                                          surf_usm_v(l)%iwghf_eb_av(m) +        &
                                          surf_usm_v(l)%iwghf_eb(m)
                       ENDDO
                    ENDDO
                    
                CASE ( 'usm_iwghf_window' )
!--                 array of heat flux from indoor window ground (wall, roof, land)
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%iwghf_eb_window_av(m) =                              &
                                          surf_usm_h%iwghf_eb_window_av(m) +           &
                                          surf_usm_h%iwghf_eb_window(m)
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%iwghf_eb_window_av(m) =                        &
                                          surf_usm_v(l)%iwghf_eb_window_av(m) +        &
                                          surf_usm_v(l)%iwghf_eb_window(m)
                       ENDDO
                    ENDDO
                    
                CASE ( 'usm_t_surf' )
!--                 surface temperature for surfaces
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%t_surf_av(m) =                               & 
                                          surf_usm_h%t_surf_av(m) +            &
                                          t_surf_h(m)
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_surf_av(m) =                         &
                                          surf_usm_v(l)%t_surf_av(m) +         &
                                          t_surf_v(l)%t(m)
                       ENDDO
                    ENDDO
                    
                CASE ( 'usm_t_surf_window' )
!--                 surface temperature for window surfaces
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%t_surf_window_av(m) =                               & 
                                          surf_usm_h%t_surf_window_av(m) +            &
                                          t_surf_window_h(m)
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_surf_window_av(m) =                         &
                                          surf_usm_v(l)%t_surf_window_av(m) +         &
                                          t_surf_window_v(l)%t(m)
                       ENDDO
                    ENDDO
                    
                CASE ( 'usm_t_surf_green' )
!--                 surface temperature for green surfaces
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%t_surf_green_av(m) =                               & 
                                          surf_usm_h%t_surf_green_av(m) +            &
                                          t_surf_green_h(m)
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_surf_green_av(m) =                         &
                                          surf_usm_v(l)%t_surf_green_av(m) +         &
                                          t_surf_green_v(l)%t(m)
                       ENDDO
                    ENDDO
                
                CASE ( 'usm_t_surf_10cm' )
!--                 near surface temperature for whole surfaces
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%t_surf_10cm_av(m) =                               & 
                                          surf_usm_h%t_surf_10cm_av(m) +            &
                                          t_surf_10cm_h(m)
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_surf_10cm_av(m) =                         &
                                          surf_usm_v(l)%t_surf_10cm_av(m) +         &
                                          t_surf_10cm_v(l)%t(m)
                       ENDDO
                    ENDDO

                    
                CASE ( 'usm_t_wall' )
!--                 wall temperature for  iwl layer of walls and land
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%t_wall_av(iwl,m) =                           &
                                          surf_usm_h%t_wall_av(iwl,m) +        &
                                          t_wall_h(iwl,m)
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_wall_av(iwl,m) =                     &
                                          surf_usm_v(l)%t_wall_av(iwl,m) +     &
                                          t_wall_v(l)%t(iwl,m)
                       ENDDO
                    ENDDO
                    
                CASE ( 'usm_t_window' )
!--                 window temperature for  iwl layer of walls and land
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%t_window_av(iwl,m) =                           &
                                          surf_usm_h%t_window_av(iwl,m) +        &
                                          t_window_h(iwl,m)
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_window_av(iwl,m) =                     &
                                          surf_usm_v(l)%t_window_av(iwl,m) +     &
                                          t_window_v(l)%t(iwl,m)
                       ENDDO
                    ENDDO

                CASE ( 'usm_t_green' )
!--                 green temperature for  iwl layer of walls and land
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%t_green_av(iwl,m) =                           &
                                          surf_usm_h%t_green_av(iwl,m) +        &
                                          t_green_h(iwl,m)
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_green_av(iwl,m) =                     &
                                          surf_usm_v(l)%t_green_av(iwl,m) +     &
                                          t_green_v(l)%t(iwl,m)
                       ENDDO
                    ENDDO

                CASE DEFAULT
                    CONTINUE

           END SELECT

        ELSEIF ( mode == 'average' )  THEN
           
           SELECT CASE ( TRIM( var ) )
                
                CASE ( 'usm_rad_net' )
!--                 array of complete radiation balance
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%rad_net_av(m) =                              &
                                          surf_usm_h%rad_net_av(m) /           &
                                          REAL( average_count_3d, kind=wp )
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%rad_net_av(m) =                        &
                                          surf_usm_v(l)%rad_net_av(m) /        &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDDO
                    
                CASE ( 'usm_rad_insw' )
!--                 array of sw radiation falling to surface after i-th reflection
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfinsw_av(l) = surfinsw_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO
                             
                CASE ( 'usm_rad_inlw' )
!--                 array of lw radiation falling to surface after i-th reflection
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfinlw_av(l) = surfinlw_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_inswdir' )
!--                 array of direct sw radiation falling to surface from sun
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfinswdir_av(l) = surfinswdir_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_inswdif' )
!--                 array of difusion sw radiation falling to surface from sky and borders of the domain
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfinswdif_av(l) = surfinswdif_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_inswref' )
!--                 array of sw radiation falling to surface from reflections
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfinswref_av(l) = surfinswref_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_inlwdif' )
!--                 array of sw radiation falling to surface after i-th reflection
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfinlwdif_av(l) = surfinlwdif_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_inlwref' )
!--                 array of lw radiation falling to surface from reflections
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfinlwref_av(l) = surfinlwref_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_outsw' )
!--                 array of sw radiation emitted from surface after i-th reflection
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfoutsw_av(l) = surfoutsw_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_outlw' )
!--                 array of lw radiation emitted from surface after i-th reflection
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfoutlw_av(l) = surfoutlw_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_ressw' )
!--                 array of residua of sw radiation absorbed in surface after last reflection
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfins_av(l) = surfins_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO
                                    
                CASE ( 'usm_rad_reslw' )
!--                 array of residua of lw radiation absorbed in surface after last reflection
                    DO l = 1, nsurfl
                        IF ( surfl(id,l) == idsint )  THEN
                            surfinl_av(l) = surfinl_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_hf' )
!--                 array of heat flux from radiation for surfaces after i-th reflection
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%surfhf_av(m) =                               &
                                          surf_usm_h%surfhf_av(m) /            &
                                          REAL( average_count_3d, kind=wp )
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%surfhf_av(m) =                         &
                                          surf_usm_v(l)%surfhf_av(m) /         &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDDO
                    
                CASE ( 'usm_wshf' )
!--                 array of sensible heat flux from surfaces (land, roof, wall)
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%wshf_eb_av(m) =                              &
                                          surf_usm_h%wshf_eb_av(m) /           &
                                          REAL( average_count_3d, kind=wp )
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%wshf_eb_av(m) =                        &
                                          surf_usm_v(l)%wshf_eb_av(m) /        &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDDO
                    
                CASE ( 'usm_wghf' )
!--                 array of heat flux from ground (wall, roof, land)
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%wghf_eb_av(m) =                              &
                                          surf_usm_h%wghf_eb_av(m) /           &
                                          REAL( average_count_3d, kind=wp )
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%wghf_eb_av(m) =                        &
                                          surf_usm_v(l)%wghf_eb_av(m) /        &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDDO
                    
                CASE ( 'usm_wghf_window' )
!--                 array of heat flux from window ground (wall, roof, land)
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%wghf_eb_window_av(m) =                              &
                                          surf_usm_h%wghf_eb_window_av(m) /           &
                                          REAL( average_count_3d, kind=wp )
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%wghf_eb_window_av(m) =                        &
                                          surf_usm_v(l)%wghf_eb_window_av(m) /        &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDDO

                CASE ( 'usm_wghf_green' )
!--                 array of heat flux from green ground (wall, roof, land)
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%wghf_eb_green_av(m) =                              &
                                          surf_usm_h%wghf_eb_green_av(m) /           &
                                          REAL( average_count_3d, kind=wp )
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%wghf_eb_green_av(m) =                        &
                                          surf_usm_v(l)%wghf_eb_green_av(m) /        &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDDO

                CASE ( 'usm_iwghf' )
!--                 array of heat flux from indoor ground (wall, roof, land)
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%iwghf_eb_av(m) =                              &
                                          surf_usm_h%iwghf_eb_av(m) /           &
                                          REAL( average_count_3d, kind=wp )
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%iwghf_eb_av(m) =                        &
                                          surf_usm_v(l)%iwghf_eb_av(m) /        &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDDO
                    
                CASE ( 'usm_iwghf_window' )
!--                 array of heat flux from indoor window ground (wall, roof, land)
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%iwghf_eb_window_av(m) =                              &
                                          surf_usm_h%iwghf_eb_window_av(m) /           &
                                          REAL( average_count_3d, kind=wp )
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%iwghf_eb_window_av(m) =                        &
                                          surf_usm_v(l)%iwghf_eb_window_av(m) /        &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDDO
                    
                CASE ( 'usm_t_surf' )
!--                 surface temperature for surfaces
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%t_surf_av(m) =                               & 
                                          surf_usm_h%t_surf_av(m) /            &
                                          REAL( average_count_3d, kind=wp )
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_surf_av(m) =                         &
                                          surf_usm_v(l)%t_surf_av(m) /         &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDDO
                    
                CASE ( 'usm_t_surf_window' )
!--                 surface temperature for window surfaces
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%t_surf_window_av(m) =                               & 
                                          surf_usm_h%t_surf_window_av(m) /            &
                                          REAL( average_count_3d, kind=wp )
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_surf_window_av(m) =                         &
                                          surf_usm_v(l)%t_surf_window_av(m) /         &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDDO
                    
                CASE ( 'usm_t_surf_green' )
!--                 surface temperature for green surfaces
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%t_surf_green_av(m) =                               & 
                                          surf_usm_h%t_surf_green_av(m) /            &
                                          REAL( average_count_3d, kind=wp )
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_surf_green_av(m) =                         &
                                          surf_usm_v(l)%t_surf_green_av(m) /         &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDDO
                    
                CASE ( 'usm_t_surf_10cm' )
!--                 near surface temperature for whole surfaces
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%t_surf_10cm_av(m) =                               & 
                                          surf_usm_h%t_surf_10cm_av(m) /            &
                                          REAL( average_count_3d, kind=wp )
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_surf_10cm_av(m) =                         &
                                          surf_usm_v(l)%t_surf_10cm_av(m) /         &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDDO
                    
                CASE ( 'usm_t_wall' )
!--                 wall temperature for  iwl layer of walls and land
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%t_wall_av(iwl,m) =                           &
                                          surf_usm_h%t_wall_av(iwl,m) /        &
                                          REAL( average_count_3d, kind=wp )
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_wall_av(iwl,m) =                     &
                                          surf_usm_v(l)%t_wall_av(iwl,m) /     &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDDO

                CASE ( 'usm_t_window' )
!--                 window temperature for  iwl layer of walls and land
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%t_window_av(iwl,m) =                           &
                                          surf_usm_h%t_window_av(iwl,m) /        &
                                          REAL( average_count_3d, kind=wp )
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_window_av(iwl,m) =                     &
                                          surf_usm_v(l)%t_window_av(iwl,m) /     &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDDO

                CASE ( 'usm_t_green' )
!--                 green temperature for  iwl layer of walls and land
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%t_green_av(iwl,m) =                           &
                                          surf_usm_h%t_green_av(iwl,m) /        &
                                          REAL( average_count_3d, kind=wp )
                    ENDDO
                    DO  l = 0, 3
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_green_av(iwl,m) =                     &
                                          surf_usm_v(l)%t_green_av(iwl,m) /     &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDDO


           END SELECT

        ENDIF

    END SUBROUTINE usm_average_3d_data



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Set internal Neumann boundary condition at outer soil grid points 
!> for temperature and humidity. 
!------------------------------------------------------------------------------!
 SUBROUTINE usm_boundary_condition
 
    IMPLICIT NONE

    INTEGER(iwp) :: i      !< grid index x-direction
    INTEGER(iwp) :: ioff   !< offset index x-direction indicating location of soil grid point
    INTEGER(iwp) :: j      !< grid index y-direction
    INTEGER(iwp) :: joff   !< offset index x-direction indicating location of soil grid point
    INTEGER(iwp) :: k      !< grid index z-direction
    INTEGER(iwp) :: koff   !< offset index x-direction indicating location of soil grid point
    INTEGER(iwp) :: l      !< running index surface-orientation
    INTEGER(iwp) :: m      !< running index surface elements

    koff = surf_usm_h%koff
    DO  m = 1, surf_usm_h%ns
       i = surf_usm_h%i(m)
       j = surf_usm_h%j(m)
       k = surf_usm_h%k(m)
       pt(k+koff,j,i) = pt(k,j,i)
    ENDDO

    DO  l = 0, 3
       ioff = surf_usm_v(l)%ioff
       joff = surf_usm_v(l)%joff
       DO  m = 1, surf_usm_v(l)%ns
          i = surf_usm_v(l)%i(m)
          j = surf_usm_v(l)%j(m)
          k = surf_usm_v(l)%k(m)
          pt(k,j+joff,i+ioff) = pt(k,j,i)
       ENDDO
    ENDDO

 END SUBROUTINE usm_boundary_condition


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine checks variables and assigns units.
!> It is called out from subroutine check_parameters.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_check_data_output( variable, unit )
        
        IMPLICIT NONE
 
        CHARACTER (len=*),INTENT(IN)    ::  variable !:
        CHARACTER (len=*),INTENT(OUT)   ::  unit     !:
        
        CHARACTER (len=varnamelength)   :: var

        var = TRIM(variable)
        IF ( var(1:12) == 'usm_rad_net_'  .OR.  var(1:13) == 'usm_rad_insw_'  .OR.        &
             var(1:13) == 'usm_rad_inlw_'  .OR.  var(1:16) == 'usm_rad_inswdir_'  .OR.    &
             var(1:16) == 'usm_rad_inswdif_'  .OR.  var(1:16) == 'usm_rad_inswref_'  .OR. &
             var(1:16) == 'usm_rad_inlwdif_'  .OR.  var(1:16) == 'usm_rad_inlwref_'  .OR. &
             var(1:14) == 'usm_rad_outsw_'  .OR.  var(1:14) == 'usm_rad_outlw_'  .OR.     &
             var(1:14) == 'usm_rad_ressw_'  .OR.  var(1:14) == 'usm_rad_reslw_'  .OR.     &
             var(1:11) == 'usm_rad_hf_'  .OR.                                             &
             var(1:9)  == 'usm_wshf_'  .OR.  var(1:9) == 'usm_wghf_' .OR.                 &
             var(1:16) == 'usm_wghf_window_' .OR. var(1:15) == 'usm_wghf_green_' .OR.     &
             var(1:10) == 'usm_iwghf_' .OR. var(1:17) == 'usm_iwghf_window_' )  THEN
            unit = 'W/m2'
        ELSE IF ( var(1:10) == 'usm_t_surf'   .OR.  var(1:10) == 'usm_t_wall' .OR.         &
                  var(1:12) == 'usm_t_window' .OR. var(1:17) == 'usm_t_surf_window' .OR.  &
                  var(1:16) == 'usm_t_surf_green'  .OR.                                   &
                  var(1:11) == 'usm_t_green' .OR.                                         &
                  var(1:15) == 'usm_t_surf_10cm')  THEN
            unit = 'K'
        ELSE IF ( var(1:9) == 'usm_surfz'  .OR.  var(1:7) == 'usm_svf'  .OR.              & 
                  var(1:7) == 'usm_dif'  .OR.  var(1:11) == 'usm_surfcat'  .OR.           &
                  var(1:11) == 'usm_surfalb'  .OR.  var(1:12) == 'usm_surfemis'  .OR.     &
                  var(1:9) == 'usm_skyvf' .OR. var(1:9) == 'usm_skyvft' )  THEN
            unit = '1'
        ELSE
            unit = 'illegal'
        ENDIF

    END SUBROUTINE usm_check_data_output


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for urban surface model
!------------------------------------------------------------------------------!
    SUBROUTINE usm_check_parameters
    
       USE control_parameters,                                                 &
           ONLY:  bc_pt_b, bc_q_b, constant_flux_layer, large_scale_forcing,   &
                  lsf_surf, topography

!
!--    Dirichlet boundary conditions are required as the surface fluxes are
!--    calculated from the temperature/humidity gradients in the urban surface
!--    model
       IF ( bc_pt_b == 'neumann'   .OR.   bc_q_b == 'neumann' )  THEN
          message_string = 'urban surface model requires setting of '//        &
                           'bc_pt_b = "dirichlet" and '//                      &
                           'bc_q_b  = "dirichlet"'
          CALL message( 'usm_check_parameters', 'PA0590', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( .NOT. TRIM(constant_flux_layer) == 'bottom' )  THEN
          message_string = 'urban surface model requires '//                   &
                           'constant_flux_layer = bottom '
          CALL message( 'usm_check_parameters', 'PA0084', 1, 2, 0, 6, 0 )
       ENDIF

       IF (  .NOT.  radiation )  THEN
          message_string = 'urban surface model requires '//                   &
                           'the radiation model to be switched on'
          CALL message( 'usm_check_parameters', 'PA0084', 1, 2, 0, 6, 0 )
       ENDIF
!        
!--    Surface forcing has to be disabled for LSF in case of enabled 
!--    urban surface module
       IF ( large_scale_forcing )  THEN
          lsf_surf = .FALSE.
       ENDIF
!
!--    Topography
       IF ( topography == 'flat' )  THEN
          message_string = 'topography /= "flat" is required '//               &
                           'when using the urban surface model'
          CALL message( 'check_parameters', 'PA0592', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    naheatlayers
       IF ( naheatlayers > nzt )  THEN
          message_string = 'number of anthropogenic heat layers '//            &
                           '"naheatlayers" can not be larger than'//           &
                           ' number of domain layers "nzt"'
          CALL message( 'check_parameters', 'PA0593', 1, 2, 0, 6, 0 )
       ENDIF

    END SUBROUTINE usm_check_parameters


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Output of the 3D-arrays in netCDF and/or AVS format 
!> for variables of urban_surface model.
!> It resorts the urban surface module output quantities from surf style
!> indexing into temporary 3D array with indices (i,j,k).
!> It is called from subroutine data_output_3d.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )
        
        IMPLICIT NONE

        INTEGER(iwp), INTENT(IN)       ::  av        !< 
        CHARACTER (len=*), INTENT(IN)  ::  variable  !< 
        INTEGER(iwp), INTENT(IN)       ::  nzb_do    !< lower limit of the data output (usually 0)
        INTEGER(iwp), INTENT(IN)       ::  nzt_do    !< vertical upper limit of the data output (usually nz_do3d)
        LOGICAL, INTENT(OUT)           ::  found     !< 
        REAL(sp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf   !< sp - it has to correspond to module data_output_3d
        REAL(wp), DIMENSION(nzb:nzt+1,nys:nyn,nxl:nxr)     ::  temp_pf    !< temp array for urban surface output procedure
        
        CHARACTER (len=varnamelength)                          :: var, surfid
        INTEGER(iwp), PARAMETER                                :: nd = 5
        CHARACTER(len=6), DIMENSION(0:nd-1), PARAMETER         :: dirname = (/ '_roof ', '_south', '_north', '_west ', '_east ' /)
        INTEGER(iwp), DIMENSION(0:nd-1), PARAMETER             :: dirint =  (/    iup_u, isouth_u, inorth_u,  iwest_u,  ieast_u /)
        INTEGER(iwp), DIMENSION(0:nd-1), PARAMETER             :: diridx =  (/       -1,        1,        0,        3,        2 /)
                                                                     !< index for surf_*_v: 0:3 = (North, South, East, West)
        INTEGER(iwp), DIMENSION(0:nd-1)                        :: dirstart
        INTEGER(iwp), DIMENSION(0:nd-1)                        :: dirend
        INTEGER(iwp)                                           :: ids,idsint,idsidx,isurf,isvf,isurfs,isurflt
        INTEGER(iwp)                                           :: is,js,ks,i,j,k,iwl,istat, l, m
        INTEGER(iwp)                                           :: k_topo    !< topography top index

        dirstart = (/ startland, startwall, startwall, startwall, startwall /)
        dirend = (/ endland, endwall, endwall, endwall, endwall /)

        found = .TRUE.
        temp_pf = -1._wp
        
        ids = -1
        var = TRIM(variable)
        DO i = 0, nd-1
            k = len(TRIM(var))
            j = len(TRIM(dirname(i)))
            IF ( var(k-j+1:k) == dirname(i) )  THEN
                ids = i
                idsint = dirint(ids)
                idsidx = diridx(ids)
                var = var(:k-j)
                EXIT
            ENDIF
        ENDDO
        IF ( ids == -1 )  THEN
            var = TRIM(variable)
        ENDIF
        IF ( var(1:11) == 'usm_t_wall_'  .AND.  len(TRIM(var)) >= 12 )  THEN
!--         wall layers
            READ(var(12:12), '(I1)', iostat=istat ) iwl
            IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
                var = var(1:10)
            ENDIF
        ENDIF
        IF ( var(1:13) == 'usm_t_window_'  .AND.  len(TRIM(var)) >= 14 )  THEN
!--         window layers
            READ(var(14:14), '(I1)', iostat=istat ) iwl
            IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
                var = var(1:12)
            ENDIF
        ENDIF
        IF ( var(1:12) == 'usm_t_green_'  .AND.  len(TRIM(var)) >= 13 )  THEN
!--         green layers
            READ(var(13:13), '(I1)', iostat=istat ) iwl
            IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
                var = var(1:11)
            ENDIF
        ENDIF
        IF ( (var(1:8) == 'usm_svf_'  .OR.  var(1:8) == 'usm_dif_')  .AND.  len(TRIM(var)) >= 13 )  THEN
!--         svf values to particular surface
            surfid = var(9:)
            i = index(surfid,'_')
            j = index(surfid(i+1:),'_')
            READ(surfid(1:i-1),*, iostat=istat ) is
            IF ( istat == 0 )  THEN
                READ(surfid(i+1:i+j-1),*, iostat=istat ) js
            ENDIF
            IF ( istat == 0 )  THEN
                READ(surfid(i+j+1:),*, iostat=istat ) ks
            ENDIF
            IF ( istat == 0 )  THEN
                var = var(1:7)
            ENDIF
        ENDIF
        
        SELECT CASE ( TRIM(var) )

          CASE ( 'usm_surfz' )
!--           array of lw radiation falling to local surface after i-th reflection
              IF ( idsint == iup_u )  THEN
                 DO  m = 1, surf_usm_h%ns
                    i = surf_usm_h%i(m)
                    j = surf_usm_h%j(m)
                    k = surf_usm_h%k(m)
                    temp_pf(0,j,i) = MAX( temp_pf(0,j,i), REAL( k, kind=wp) )
                 ENDDO
              ELSE
                 l = idsidx
                 DO  m = 1, surf_usm_v(l)%ns
                    i = surf_usm_v(l)%i(m)
                    j = surf_usm_v(l)%j(m)
                    k = surf_usm_v(l)%k(m)
                    temp_pf(0,j,i) = MAX( temp_pf(0,j,i), REAL( k, kind=wp) + 1.0_wp )
                 ENDDO
              ENDIF

          CASE ( 'usm_surfcat' )
!--           surface category
              IF ( idsint == iup_u )  THEN
                 DO  m = 1, surf_usm_h%ns
                    i = surf_usm_h%i(m)
                    j = surf_usm_h%j(m)
                    k = surf_usm_h%k(m)
                    temp_pf(k,j,i) = surf_usm_h%surface_types(m)
                 ENDDO
              ELSE
                 l = idsidx
                 DO  m = 1, surf_usm_v(l)%ns
                    i = surf_usm_v(l)%i(m)
                    j = surf_usm_v(l)%j(m)
                    k = surf_usm_v(l)%k(m)
                    temp_pf(k,j,i) = surf_usm_v(l)%surface_types(m)
                 ENDDO
              ENDIF
              
          CASE ( 'usm_surfalb' )
!--           surface albedo, weighted average
              IF ( idsint == iup_u )  THEN
                 DO  m = 1, surf_usm_h%ns
                    i = surf_usm_h%i(m)
                    j = surf_usm_h%j(m)
                    k = surf_usm_h%k(m)
                    temp_pf(k,j,i) = surf_usm_h%frac(ind_veg_wall,m)     *     &
                                     surf_usm_h%albedo(ind_veg_wall,m)  +      &
                                     surf_usm_h%frac(ind_pav_green,m)    *     &
                                     surf_usm_h%albedo(ind_pav_green,m) +      &
                                     surf_usm_h%frac(ind_wat_win,m)      *     &
                                     surf_usm_h%albedo(ind_wat_win,m)
                 ENDDO
              ELSE
                 l = idsidx
                 DO  m = 1, surf_usm_v(l)%ns
                    i = surf_usm_v(l)%i(m)
                    j = surf_usm_v(l)%j(m)
                    k = surf_usm_v(l)%k(m)
                    temp_pf(k,j,i) = surf_usm_v(l)%frac(ind_veg_wall,m)     *  &
                                     surf_usm_v(l)%albedo(ind_veg_wall,m)  +   &
                                     surf_usm_v(l)%frac(ind_pav_green,m)    *  &
                                     surf_usm_v(l)%albedo(ind_pav_green,m) +   &
                                     surf_usm_v(l)%frac(ind_wat_win,m)      *  &
                                     surf_usm_v(l)%albedo(ind_wat_win,m)
                 ENDDO
              ENDIF
              
          CASE ( 'usm_surfemis' )
!--           surface emissivity, weighted average
              IF ( idsint == iup_u )  THEN
                 DO  m = 1, surf_usm_h%ns
                    i = surf_usm_h%i(m)
                    j = surf_usm_h%j(m)
                    k = surf_usm_h%k(m)
                    temp_pf(k,j,i) =  surf_usm_h%frac(ind_veg_wall,m)      *   &
                                      surf_usm_h%emissivity(ind_veg_wall,m)  + &
                                      surf_usm_h%frac(ind_pav_green,m)     *   &
                                      surf_usm_h%emissivity(ind_pav_green,m) + &
                                      surf_usm_h%frac(ind_wat_win,m)       *   &
                                      surf_usm_h%emissivity(ind_wat_win,m)
                 ENDDO
              ELSE
                 l = idsidx
                 DO  m = 1, surf_usm_v(l)%ns
                    i = surf_usm_v(l)%i(m)
                    j = surf_usm_v(l)%j(m)
                    k = surf_usm_v(l)%k(m)
                    temp_pf(k,j,i) = surf_usm_v(l)%frac(ind_veg_wall,m)       *&
                                     surf_usm_v(l)%emissivity(ind_veg_wall,m) +&
                                     surf_usm_v(l)%frac(ind_pav_green,m)      *&
                                     surf_usm_v(l)%emissivity(ind_pav_green,m)+&
                                     surf_usm_v(l)%frac(ind_wat_win,m)        *&
                                     surf_usm_v(l)%emissivity(ind_wat_win,m)
                 ENDDO
              ENDIF

          CASE ( 'usm_surfwintrans' )
!--           transmissivity window tiles
              IF ( idsint == iup_u )  THEN
                 DO  m = 1, surf_usm_h%ns
                    i = surf_usm_h%i(m)
                    j = surf_usm_h%j(m)
                    k = surf_usm_h%k(m)
                    temp_pf(k,j,i) = surf_usm_h%transmissivity(m)
                 ENDDO
              ELSE
                 l = idsidx
                 DO  m = 1, surf_usm_v(l)%ns
                    i = surf_usm_v(l)%i(m)
                    j = surf_usm_v(l)%j(m)
                    k = surf_usm_v(l)%k(m)
                    temp_pf(k,j,i) = surf_usm_v(l)%transmissivity(m)
                 ENDDO
              ENDIF

          CASE ( 'usm_skyvf' )
!--           sky view factor
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == idsint )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = skyvf(isurf)
                 ENDIF
              ENDDO
              
          CASE ( 'usm_skyvft' )
!--           sky view factor
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == ids )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = skyvft(isurf)
                 ENDIF
              ENDDO

!
!-- Not adjusted so far              
          CASE ( 'usm_svf', 'usm_dif' )
!--           shape view factors or iradiance factors to selected surface
              IF ( TRIM(var)=='usm_svf' )  THEN
                  k = 1
              ELSE
                  k = 2
              ENDIF
              DO isvf = 1, nsvfl
                  isurflt = svfsurf(1, isvf)
                  isurfs = svfsurf(2, isvf)
                             
                  IF ( surf(ix,isurfs) == is  .AND.  surf(iy,isurfs) == js  .AND.       &
                       surf(iz,isurfs) == ks  .AND.  surf(id,isurfs) == idsint )  THEN
  !--                 correct source surface
                      temp_pf(surfl(iz,isurflt),surfl(iy,isurflt),surfl(ix,isurflt)) = svf(k,isvf)
                  ENDIF
              ENDDO

          CASE ( 'usm_rad_net' )
!--           array of complete radiation balance
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%rad_net_l(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%rad_net_l(m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%rad_net_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%rad_net_av(m)
                    ENDDO
                 ENDIF
              ENDIF

          CASE ( 'usm_rad_insw' )
!--           array of sw radiation falling to surface after i-th reflection
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == idsint )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinsw(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinsw_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_inlw' )
!--           array of lw radiation falling to surface after i-th reflection
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == idsint )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinlw(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinlw_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_inswdir' )
!--           array of direct sw radiation falling to surface from sun
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == idsint )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinswdir(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinswdir_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_inswdif' )
!--           array of difusion sw radiation falling to surface from sky and borders of the domain
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == idsint )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinswdif(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinswdif_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_inswref' )
!--           array of sw radiation falling to surface from reflections
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == idsint )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = &
                       surfinsw(isurf) - surfinswdir(isurf) - surfinswdif(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinswref_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_inlwdif' )
!--           array of difusion lw radiation falling to surface from sky and borders of the domain
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == idsint )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinlwdif(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinlwdif_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_inlwref' )
!--           array of lw radiation falling to surface from reflections
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == idsint )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinlw(isurf) - surfinlwdif(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinlwref_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_outsw' )
!--           array of sw radiation emitted from surface after i-th reflection
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == idsint )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfoutsw(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfoutsw_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_outlw' )
!--           array of lw radiation emitted from surface after i-th reflection
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == idsint )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfoutlw(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfoutlw_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_ressw' )
!--           average of array of residua of sw radiation absorbed in surface after last reflection
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == idsint )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfins(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfins_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_reslw' )
!--           average of array of residua of lw radiation absorbed in surface after last reflection
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == idsint )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinl(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinl_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO
 
          CASE ( 'usm_rad_hf' )
!--           array of heat flux from radiation for surfaces after all reflections
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%surfhf(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%surfhf(m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%surfhf_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%surfhf_av(m)
                    ENDDO
                 ENDIF
              ENDIF
 
          CASE ( 'usm_wshf' )
!--           array of sensible heat flux from surfaces
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%wshf_eb(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%wshf_eb(m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%wshf_eb_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%wshf_eb_av(m)
                    ENDDO
                 ENDIF
              ENDIF


          CASE ( 'usm_wghf' )
!--           array of heat flux from ground (land, wall, roof)
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%wghf_eb(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%wghf_eb(m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%wghf_eb_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%wghf_eb_av(m)
                    ENDDO
                 ENDIF
              ENDIF

          CASE ( 'usm_wghf_window' )
!--           array of heat flux from window ground (land, wall, roof)

              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%wghf_eb_window(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%wghf_eb_window(m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%wghf_eb_window_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%wghf_eb_window_av(m)
                    ENDDO
                 ENDIF
              ENDIF

          CASE ( 'usm_wghf_green' )
!--           array of heat flux from green ground (land, wall, roof)

              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%wghf_eb_green(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%wghf_eb_green(m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%wghf_eb_green_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%wghf_eb_green_av(m)
                    ENDDO
                 ENDIF
              ENDIF

          CASE ( 'usm_iwghf' )
!--           array of heat flux from indoor ground (land, wall, roof)
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%iwghf_eb(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%iwghf_eb(m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%iwghf_eb_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%iwghf_eb_av(m)
                    ENDDO
                 ENDIF
              ENDIF

          CASE ( 'usm_iwghf_window' )
!--           array of heat flux from indoor window ground (land, wall, roof)

              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%iwghf_eb_window(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%iwghf_eb_window(m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%iwghf_eb_window_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%iwghf_eb_window_av(m)
                    ENDDO
                 ENDIF
              ENDIF
              
          CASE ( 'usm_t_surf' )
!--           surface temperature for surfaces
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = t_surf_h(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = t_surf_v(l)%t(m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%t_surf_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%t_surf_av(m)
                    ENDDO
                 ENDIF
              ENDIF
              
          CASE ( 'usm_t_surf_window' )
!--           surface temperature for window surfaces

              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = t_surf_window_h(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = t_surf_window_v(l)%t(m)
                    ENDDO
                 ENDIF

              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%t_surf_window_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%t_surf_window_av(m)
                    ENDDO

                 ENDIF

              ENDIF

          CASE ( 'usm_t_surf_green' )
!--           surface temperature for green surfaces

              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = t_surf_green_h(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = t_surf_green_v(l)%t(m)
                    ENDDO
                 ENDIF

              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%t_surf_green_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%t_surf_green_av(m)
                    ENDDO

                 ENDIF

              ENDIF

          CASE ( 'usm_t_surf_10cm' )
!--           near surface temperature for whole surfaces

              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = t_surf_10cm_h(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = t_surf_10cm_v(l)%t(m)
                    ENDDO
                 ENDIF

              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%t_surf_10cm_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%t_surf_10cm_av(m)
                    ENDDO

                 ENDIF

              ENDIF

              
          CASE ( 'usm_t_wall' )
!--           wall temperature for  iwl layer of walls and land
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = t_wall_h(iwl,m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = t_wall_v(l)%t(iwl,m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%t_wall_av(iwl,m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%t_wall_av(iwl,m)
                    ENDDO
                 ENDIF
              ENDIF
             
          CASE ( 'usm_t_window' )
!--           window temperature for iwl layer of walls and land
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = t_window_h(iwl,m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = t_window_v(l)%t(iwl,m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%t_window_av(iwl,m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%t_window_av(iwl,m)
                    ENDDO
                 ENDIF
              ENDIF

          CASE ( 'usm_t_green' )
!--           green temperature for  iwl layer of walls and land
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = t_green_h(iwl,m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = t_green_v(l)%t(iwl,m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%t_green_av(iwl,m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%t_green_av(iwl,m)
                    ENDDO
                 ENDIF
              ENDIF

             
          CASE DEFAULT
              found = .FALSE.
              
        END SELECT

!
!--     Rearrange dimensions for NetCDF output
!--     FIXME: this may generate FPE overflow upon conversion from DP to SP
        DO  j = nys, nyn
            DO  i = nxl, nxr
                DO  k = nzb_do, nzt_do
                    local_pf(i,j,k) = temp_pf(k,j,i)
                ENDDO
            ENDDO
        ENDDO
        
    END SUBROUTINE usm_data_output_3d
    

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Soubroutine defines appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_define_netcdf_grid( variable, found, grid_x, grid_y, grid_z )
    
        IMPLICIT NONE

        CHARACTER (len=*), INTENT(IN)  ::  variable    !< 
        LOGICAL, INTENT(OUT)           ::  found       !< 
        CHARACTER (len=*), INTENT(OUT) ::  grid_x      !< 
        CHARACTER (len=*), INTENT(OUT) ::  grid_y      !< 
        CHARACTER (len=*), INTENT(OUT) ::  grid_z      !< 

        CHARACTER (len=varnamelength)  :: var

        var = TRIM(variable)
        IF ( var(1:12) == 'usm_rad_net_'  .OR.  var(1:13) == 'usm_rad_insw_'  .OR.          &
             var(1:13) == 'usm_rad_inlw_'  .OR.  var(1:16) == 'usm_rad_inswdir_'  .OR.      &
             var(1:16) == 'usm_rad_inswdif_'  .OR.  var(1:16) == 'usm_rad_inswref_'  .OR.   &
             var(1:16) == 'usm_rad_inlwdif_'  .OR.  var(1:16) == 'usm_rad_inlwref_'  .OR.   &
             var(1:14) == 'usm_rad_outsw_'  .OR.  var(1:14) == 'usm_rad_outlw_'  .OR.       &
             var(1:14) == 'usm_rad_ressw_'  .OR.  var(1:14) == 'usm_rad_reslw_'  .OR.       &
             var(1:11) == 'usm_rad_hf_'  .OR.                                               &
             var(1:9) == 'usm_wshf_'  .OR.  var(1:9) == 'usm_wghf_'  .OR.                   &
             var(1:16) == 'usm_wghf_window_'  .OR. var(1:15) == 'usm_wghf_green_' .OR.      &
             var(1:10) == 'usm_iwghf_'  .OR. var(1:17) == 'usm_iwghf_window_' .OR.          &
             var(1:10) == 'usm_t_surf'  .OR.  var(1:10) == 'usm_t_wall'  .OR.               &
             var(1:17) == 'usm_t_surf_window'  .OR.  var(1:12) == 'usm_t_window'  .OR.      &
             var(1:16) == 'usm_t_surf_green'  .OR.                                          &
             var(1:15) == 'usm_t_surf_10cm' .OR.                                            &
             var(1:9) == 'usm_surfz'  .OR.  var(1:7) == 'usm_svf'  .OR.                     & 
             var(1:7) == 'usm_dif'  .OR.  var(1:11) == 'usm_surfcat'  .OR.                  &
             var(1:11) == 'usm_surfalb'  .OR.  var(1:12) == 'usm_surfemis'  .OR.            &
             var(1:16) == 'usm_surfwintrans'  .OR.                                          &
             var(1:9) == 'usm_skyvf' .OR. var(1:9) == 'usm_skyvft' ) THEN

            found = .TRUE.
            grid_x = 'x'
            grid_y = 'y'
            grid_z = 'zu'
        ELSE
            found  = .FALSE.
            grid_x = 'none'
            grid_y = 'none'
            grid_z = 'none'
        ENDIF

    END SUBROUTINE usm_define_netcdf_grid
    

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the wall surface model
!------------------------------------------------------------------------------!
    SUBROUTINE usm_init_material_model

        IMPLICIT NONE

        INTEGER(iwp) ::  k, l, m            !< running indices
        
        CALL location_message( '    initialization of wall surface model', .TRUE. )
       
!--     Calculate wall grid spacings. 
!--     Temperature is defined at the center of the wall layers,
!--     whereas gradients/fluxes are defined at the edges (_stag)      
!--     apply for all particular surface grids. First for horizontal surfaces
        DO  m = 1, surf_usm_h%ns

           surf_usm_h%dz_wall(nzb_wall,m) = surf_usm_h%zw(nzb_wall,m)
           DO k = nzb_wall+1, nzt_wall
               surf_usm_h%dz_wall(k,m) = surf_usm_h%zw(k,m) -                  &
                                         surf_usm_h%zw(k-1,m)
           ENDDO
           surf_usm_h%dz_window(nzb_wall,m) = surf_usm_h%zw_window(nzb_wall,m)
           DO k = nzb_wall+1, nzt_wall
               surf_usm_h%dz_window(k,m) = surf_usm_h%zw_window(k,m) -         &
                                         surf_usm_h%zw_window(k-1,m)
           ENDDO
           surf_usm_h%dz_green(nzb_wall,m) = surf_usm_h%zw_green(nzb_wall,m)
           DO k = nzb_wall+1, nzt_wall
               surf_usm_h%dz_green(k,m) = surf_usm_h%zw_green(k,m) -           &
                                         surf_usm_h%zw_green(k-1,m)
           ENDDO
           
           surf_usm_h%dz_wall(nzt_wall+1,m) = surf_usm_h%dz_wall(nzt_wall,m)

           DO k = nzb_wall, nzt_wall-1
               surf_usm_h%dz_wall_stag(k,m) = 0.5 * (                          &
                           surf_usm_h%dz_wall(k+1,m) + surf_usm_h%dz_wall(k,m) )
           ENDDO
           surf_usm_h%dz_wall_stag(nzt_wall,m) = surf_usm_h%dz_wall(nzt_wall,m)
           
           surf_usm_h%dz_window(nzt_wall+1,m) = surf_usm_h%dz_window(nzt_wall,m)

           DO k = nzb_wall, nzt_wall-1
               surf_usm_h%dz_window_stag(k,m) = 0.5 * (                        &
                           surf_usm_h%dz_window(k+1,m) + surf_usm_h%dz_window(k,m) )
           ENDDO
           surf_usm_h%dz_window_stag(nzt_wall,m) = surf_usm_h%dz_window(nzt_wall,m)

           surf_usm_h%dz_green(nzt_wall+1,m) = surf_usm_h%dz_green(nzt_wall,m)

           DO k = nzb_wall, nzt_wall-1
               surf_usm_h%dz_green_stag(k,m) = 0.5 * (                         &
                           surf_usm_h%dz_green(k+1,m) + surf_usm_h%dz_green(k,m) )
           ENDDO
           surf_usm_h%dz_green_stag(nzt_wall,m) = surf_usm_h%dz_green(nzt_wall,m)
        ENDDO
        surf_usm_h%ddz_wall        = 1.0_wp / surf_usm_h%dz_wall
        surf_usm_h%ddz_wall_stag   = 1.0_wp / surf_usm_h%dz_wall_stag
        surf_usm_h%ddz_window      = 1.0_wp / surf_usm_h%dz_window
        surf_usm_h%ddz_window_stag = 1.0_wp / surf_usm_h%dz_window_stag
        surf_usm_h%ddz_green       = 1.0_wp / surf_usm_h%dz_green
        surf_usm_h%ddz_green_stag  = 1.0_wp / surf_usm_h%dz_green_stag
!        
!--     For vertical surfaces
        DO  l = 0, 3
           DO  m = 1, surf_usm_v(l)%ns
              surf_usm_v(l)%dz_wall(nzb_wall,m) = surf_usm_v(l)%zw(nzb_wall,m)
              DO k = nzb_wall+1, nzt_wall
                  surf_usm_v(l)%dz_wall(k,m) = surf_usm_v(l)%zw(k,m) -         &
                                               surf_usm_v(l)%zw(k-1,m)
              ENDDO
              surf_usm_v(l)%dz_window(nzb_wall,m) = surf_usm_v(l)%zw_window(nzb_wall,m)
              DO k = nzb_wall+1, nzt_wall
                  surf_usm_v(l)%dz_window(k,m) = surf_usm_v(l)%zw_window(k,m) - &
                                               surf_usm_v(l)%zw_window(k-1,m)
              ENDDO
              surf_usm_v(l)%dz_green(nzb_wall,m) = surf_usm_v(l)%zw_green(nzb_wall,m)
              DO k = nzb_wall+1, nzt_wall
                  surf_usm_v(l)%dz_green(k,m) = surf_usm_v(l)%zw_green(k,m) - &
                                               surf_usm_v(l)%zw_green(k-1,m)
              ENDDO
           
              surf_usm_v(l)%dz_wall(nzt_wall+1,m) =                            &
                                              surf_usm_v(l)%dz_wall(nzt_wall,m)

              DO k = nzb_wall, nzt_wall-1
                  surf_usm_v(l)%dz_wall_stag(k,m) = 0.5 * (                    &
                                                surf_usm_v(l)%dz_wall(k+1,m) + &
                                                surf_usm_v(l)%dz_wall(k,m) )
              ENDDO
              surf_usm_v(l)%dz_wall_stag(nzt_wall,m) =                         &
                                              surf_usm_v(l)%dz_wall(nzt_wall,m)
              surf_usm_v(l)%dz_window(nzt_wall+1,m) =                            &
                                              surf_usm_v(l)%dz_window(nzt_wall,m)

              DO k = nzb_wall, nzt_wall-1
                  surf_usm_v(l)%dz_window_stag(k,m) = 0.5 * (                    &
                                                surf_usm_v(l)%dz_window(k+1,m) + &
                                                surf_usm_v(l)%dz_window(k,m) )
              ENDDO
              surf_usm_v(l)%dz_window_stag(nzt_wall,m) =                         &
                                              surf_usm_v(l)%dz_window(nzt_wall,m)
              surf_usm_v(l)%dz_green(nzt_wall+1,m) =                            &
                                              surf_usm_v(l)%dz_green(nzt_wall,m)

              DO k = nzb_wall, nzt_wall-1
                  surf_usm_v(l)%dz_green_stag(k,m) = 0.5 * (                    &
                                                surf_usm_v(l)%dz_green(k+1,m) + &
                                                surf_usm_v(l)%dz_green(k,m) )
              ENDDO
              surf_usm_v(l)%dz_green_stag(nzt_wall,m) =                         &
                                              surf_usm_v(l)%dz_green(nzt_wall,m)
           ENDDO
           surf_usm_v(l)%ddz_wall        = 1.0_wp / surf_usm_v(l)%dz_wall
           surf_usm_v(l)%ddz_wall_stag   = 1.0_wp / surf_usm_v(l)%dz_wall_stag
           surf_usm_v(l)%ddz_window      = 1.0_wp / surf_usm_v(l)%dz_window
           surf_usm_v(l)%ddz_window_stag = 1.0_wp / surf_usm_v(l)%dz_window_stag
           surf_usm_v(l)%ddz_green       = 1.0_wp / surf_usm_v(l)%dz_green
           surf_usm_v(l)%ddz_green_stag  = 1.0_wp / surf_usm_v(l)%dz_green_stag
        ENDDO      

        
        CALL location_message( '    wall structures filed out', .TRUE. )

        CALL location_message( '    initialization of wall surface model finished', .TRUE. )

    END SUBROUTINE usm_init_material_model

 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the urban surface model
!------------------------------------------------------------------------------!
    SUBROUTINE usm_init_urban_surface

        USE arrays_3d,                                                         &
            ONLY:  zw

        USE netcdf_data_input_mod,                                             &
            ONLY:  building_pars_f, building_type_f, terrain_height_f
    
        IMPLICIT NONE

        INTEGER(iwp) ::  i                   !< loop index x-dirction
        INTEGER(iwp) ::  ind_emis_wall       !< index in input list for wall emissivity
        INTEGER(iwp) ::  ind_emis_green      !< index in input list for green emissivity
        INTEGER(iwp) ::  ind_emis_win        !< index in input list for window emissivity
        INTEGER(iwp) ::  ind_green_frac_w    !< index in input list for green fraction on wall
        INTEGER(iwp) ::  ind_green_frac_r    !< index in input list for green fraction on roof
        INTEGER(iwp) ::  ind_hc1             !< index in input list for heat capacity at first wall layer
        INTEGER(iwp) ::  ind_hc2             !< index in input list for heat capacity at second wall layer
        INTEGER(iwp) ::  ind_hc3             !< index in input list for heat capacity at third wall layer
        INTEGER(iwp) ::  ind_lai_r           !< index in input list for LAI on roof
        INTEGER(iwp) ::  ind_lai_w           !< index in input list for LAI on wall
        INTEGER(iwp) ::  ind_tc1             !< index in input list for thermal conductivity at first wall layer
        INTEGER(iwp) ::  ind_tc2             !< index in input list for thermal conductivity at second wall layer
        INTEGER(iwp) ::  ind_tc3             !< index in input list for thermal conductivity at third wall layer
        INTEGER(iwp) ::  ind_trans           !< index in input list for window transmissivity
        INTEGER(iwp) ::  ind_wall_frac       !< index in input list for wall fraction
        INTEGER(iwp) ::  ind_win_frac        !< index in input list for window fraction
        INTEGER(iwp) ::  ind_z0              !< index in input list for z0
        INTEGER(iwp) ::  ind_z0qh            !< index in input list for z0h / z0q
        INTEGER(iwp) ::  j                   !< loop index y-dirction
        INTEGER(iwp) ::  k                   !< loop index z-dirction
        INTEGER(iwp) ::  l                   !< loop index surface orientation
        INTEGER(iwp) ::  m                   !< loop index surface element
        INTEGER(iwp) ::  st                  !< dummy  

        REAL(wp)     ::  c, d, tin, twin
        REAL(wp)     ::  ground_floor_level_l         !< local height of ground floor level
        REAL(wp)     ::  z_agl                        !< height above ground
        REAL(wp), DIMENSION(nzb:nzt)   ::  exn        !< value of the Exner function in layers

!
!-- NOPOINTER version not implemented yet
#if defined( __nopointer )
    message_string = 'The urban surface module only runs with POINTER version'
    CALL message( 'urban_surface_mod', 'PA0452', 1, 2, 0, 6, 0 )
#endif

        CALL cpu_log( log_point_s(78), 'usm_init', 'start' )
!--     surface forcing have to be disabled for LSF 
!--     in case of enabled urban surface module
        IF ( large_scale_forcing )  THEN
            lsf_surf = .FALSE.
        ENDIF

!
!--     Flag surface elements belonging to the ground floor level. Therefore, 
!--     use terrain height array from file, if available. This flag is later used
!--     to control initialization of surface attributes.
        surf_usm_h%ground_level = .FALSE. 
        DO  m = 1, surf_usm_h%ns
           i = surf_usm_h%i(m)
           j = surf_usm_h%j(m)
           k = surf_usm_h%k(m)
!
!--        Get local ground level. If no ground level is given in input file,
!--        use default value.
           ground_floor_level_l = ground_floor_level
           IF ( building_pars_f%from_file )  THEN
              IF ( building_pars_f%pars_xy(ind_gflh,j,i) /=                    &
                   building_pars_f%fill )  &
                 ground_floor_level_l = building_pars_f%pars_xy(ind_gflh,j,i)          
           ENDIF
!
!--        Determine height of surface element above ground level
           IF (  terrain_height_f%from_file )  THEN
              z_agl = zw(k) - terrain_height_f%var(j,i)
           ELSE
              z_agl = zw(k)
           ENDIF
!
!--        Set flag for ground level
           IF ( z_agl <= ground_floor_level_l )                                &
              surf_usm_h%ground_level(m) = .TRUE.
        ENDDO

        DO  l = 0, 3
           surf_usm_v(l)%ground_level = .FALSE.
           DO  m = 1, surf_usm_v(l)%ns
              i = surf_usm_v(l)%i(m) + surf_usm_v(l)%ioff
              j = surf_usm_v(l)%j(m) + surf_usm_v(l)%joff
              k = surf_usm_v(l)%k(m)
!
!--           Get local ground level. If no ground level is given in input file,
!--           use default value.
              ground_floor_level_l = ground_floor_level
              IF ( building_pars_f%from_file )  THEN
                 IF ( building_pars_f%pars_xy(ind_gflh,j,i) /=                 &
                      building_pars_f%fill ) &
                    ground_floor_level_l = building_pars_f%pars_xy(ind_gflh,j,i)
              ENDIF
!
!--           Determine height of surface element above ground level. Please 
!--           note, height of surface element is determined with respect to
!--           its height of the adjoing atmospheric grid point. 
              IF (  terrain_height_f%from_file )  THEN
                 z_agl = zw(k) - terrain_height_f%var(j-surf_usm_v(l)%joff,    &
                                                      i-surf_usm_v(l)%ioff)
              ELSE
                 z_agl = zw(k)
              ENDIF
!
!--           Set flag for ground level
              IF ( z_agl <= ground_floor_level_l )                                &
                 surf_usm_v(l)%ground_level(m) = .TRUE.

           ENDDO
        ENDDO
!
!--     Initialization of resistances. 
        DO  m = 1, surf_usm_h%ns
           surf_usm_h%r_a(m)        = 50.0_wp
           surf_usm_h%r_a_green(m)  = 50.0_wp
           surf_usm_h%r_a_window(m) = 50.0_wp
        ENDDO
        DO  l = 0, 3
           DO  m = 1, surf_usm_v(l)%ns
              surf_usm_v(l)%r_a(m)        = 50.0_wp
              surf_usm_v(l)%r_a_green(m)  = 50.0_wp
              surf_usm_v(l)%r_a_window(m) = 50.0_wp
           ENDDO
        ENDDO
!
!--     Initialize urban-type surface attribute. According to initialization in 
!--     land-surface model, follow a 3-level approach.
!--     Level 1 - initialization via default attributes 
        DO  m = 1, surf_usm_h%ns
!
!--        Now, all horizontal surfaces are roof surfaces (?)
           surf_usm_h%isroof_surf(m)   = .TRUE.
           surf_usm_h%surface_types(m) = roof_category         !< default category for root surface
!
!--        In order to distinguish between ground floor level and 
!--        above-ground-floor level surfaces, set input indices. 
           ind_wall_frac    = MERGE( ind_wall_frac_gfl,    ind_wall_frac_agfl,    &
                                     surf_usm_h%ground_level(m) )
           ind_win_frac     = MERGE( ind_win_frac_gfl,     ind_win_frac_agfl,     &
                                     surf_usm_h%ground_level(m) )
           ind_green_frac_w = MERGE( ind_green_frac_w_gfl, ind_green_frac_w_agfl, &
                                     surf_usm_h%ground_level(m) )
           ind_green_frac_r = MERGE( ind_green_frac_r_gfl, ind_green_frac_r_agfl, &
                                     surf_usm_h%ground_level(m) )
           ind_lai_r        = MERGE( ind_lai_r_gfl,        ind_lai_r_agfl,        &
                                     surf_usm_h%ground_level(m) )
           ind_lai_w        = MERGE( ind_lai_w_gfl,        ind_lai_w_agfl,        &
                                     surf_usm_h%ground_level(m) )
           ind_hc1          = MERGE( ind_hc1_gfl,          ind_hc1_agfl,          &
                                     surf_usm_h%ground_level(m) )
           ind_hc2          = MERGE( ind_hc2_gfl,          ind_hc2_agfl,          &
                                     surf_usm_h%ground_level(m) )
           ind_hc3          = MERGE( ind_hc3_gfl,          ind_hc3_agfl,          &
                                     surf_usm_h%ground_level(m) )
           ind_tc1          = MERGE( ind_tc1_gfl,          ind_tc1_agfl,          &
                                     surf_usm_h%ground_level(m) )
           ind_tc2          = MERGE( ind_tc2_gfl,          ind_tc2_agfl,          &
                                     surf_usm_h%ground_level(m) )
           ind_tc3          = MERGE( ind_tc3_gfl,          ind_tc3_agfl,          &
                                     surf_usm_h%ground_level(m) )
           ind_emis_wall    = MERGE( ind_emis_wall_gfl,    ind_emis_wall_agfl,    &
                                     surf_usm_h%ground_level(m) )
           ind_emis_green   = MERGE( ind_emis_green_gfl,   ind_emis_green_agfl,   &
                                     surf_usm_h%ground_level(m) )
           ind_emis_win     = MERGE( ind_emis_win_gfl,     ind_emis_win_agfl,     &
                                     surf_usm_h%ground_level(m) )
           ind_trans        = MERGE( ind_trans_gfl,        ind_trans_agfl,        &
                                     surf_usm_h%ground_level(m) )
           ind_z0           = MERGE( ind_z0_gfl,           ind_z0_agfl,           &
                                     surf_usm_h%ground_level(m) )
           ind_z0qh         = MERGE( ind_z0qh_gfl,         ind_z0qh_agfl,         &
                                     surf_usm_h%ground_level(m) )
!
!--        Initialize relatvie wall- (0), green- (1) and window (2) fractions
           surf_usm_h%frac(ind_veg_wall,m)  = building_pars(ind_wall_frac,building_type)   
           surf_usm_h%frac(ind_pav_green,m) = building_pars(ind_green_frac_r,building_type)  
           surf_usm_h%frac(ind_wat_win,m)   = building_pars(ind_win_frac,building_type)  
           surf_usm_h%lai(m)                = building_pars(ind_green_frac_r,building_type)  

           surf_usm_h%rho_c_wall(nzb_wall,m)   = building_pars(ind_hc1,building_type)  
           surf_usm_h%rho_c_wall(nzb_wall+1,m) = building_pars(ind_hc1,building_type)
           surf_usm_h%rho_c_wall(nzb_wall+2,m) = building_pars(ind_hc2,building_type)
           surf_usm_h%rho_c_wall(nzb_wall+3,m) = building_pars(ind_hc3,building_type)    
           surf_usm_h%lambda_h(nzb_wall,m)   = building_pars(ind_tc1,building_type)  
           surf_usm_h%lambda_h(nzb_wall+1,m) = building_pars(ind_tc1,building_type) 
           surf_usm_h%lambda_h(nzb_wall+2,m) = building_pars(ind_tc2,building_type)
           surf_usm_h%lambda_h(nzb_wall+3,m) = building_pars(ind_tc3,building_type)    
           surf_usm_h%rho_c_green(nzb_wall,m)   = building_pars(ind_hc1,building_type)  
           surf_usm_h%rho_c_green(nzb_wall+1,m) = building_pars(ind_hc1,building_type)
           surf_usm_h%rho_c_green(nzb_wall+2,m) = building_pars(ind_hc2,building_type)
           surf_usm_h%rho_c_green(nzb_wall+3,m) = building_pars(ind_hc3,building_type)    
           surf_usm_h%lambda_h_green(nzb_wall,m)   = building_pars(ind_tc1,building_type)  
           surf_usm_h%lambda_h_green(nzb_wall+1,m) = building_pars(ind_tc1,building_type) 
           surf_usm_h%lambda_h_green(nzb_wall+2,m) = building_pars(ind_tc2,building_type)
           surf_usm_h%lambda_h_green(nzb_wall+3,m) = building_pars(ind_tc3,building_type)
           surf_usm_h%rho_c_window(nzb_wall,m)   = building_pars(ind_hc1,building_type)  
           surf_usm_h%rho_c_window(nzb_wall+1,m) = building_pars(ind_hc1,building_type)
           surf_usm_h%rho_c_window(nzb_wall+2,m) = building_pars(ind_hc2,building_type)
           surf_usm_h%rho_c_window(nzb_wall+3,m) = building_pars(ind_hc3,building_type)    
           surf_usm_h%lambda_h_window(nzb_wall,m)   = building_pars(ind_tc1,building_type)  
           surf_usm_h%lambda_h_window(nzb_wall+1,m) = building_pars(ind_tc1,building_type) 
           surf_usm_h%lambda_h_window(nzb_wall+2,m) = building_pars(ind_tc2,building_type)
           surf_usm_h%lambda_h_window(nzb_wall+3,m) = building_pars(ind_tc3,building_type)    

           surf_usm_h%target_temp_summer(m)  = building_pars(12,building_type)    
           surf_usm_h%target_temp_winter(m)  = building_pars(13,building_type)    
!
!--        emissivity of wall-, green- and window fraction 
           surf_usm_h%emissivity(ind_veg_wall,m)  = building_pars(ind_emis_wall,building_type)
           surf_usm_h%emissivity(ind_pav_green,m) = building_pars(ind_emis_green,building_type)
           surf_usm_h%emissivity(ind_wat_win,m)   = building_pars(ind_emis_win,building_type)

           surf_usm_h%transmissivity(m)      = building_pars(ind_trans,building_type)

           surf_usm_h%z0(m)                  = building_pars(ind_z0,building_type)
           surf_usm_h%z0h(m)                 = building_pars(ind_z0qh,building_type)
           surf_usm_h%z0q(m)                 = building_pars(ind_z0qh,building_type)
!
!--        albedo type for wall fraction, green fraction, window fraction
           surf_usm_h%albedo_type(ind_veg_wall,m)  = INT( building_pars(ind_alb_wall,building_type)  )
           surf_usm_h%albedo_type(ind_pav_green,m) = INT( building_pars(ind_alb_green,building_type) )
           surf_usm_h%albedo_type(ind_wat_win,m)   = INT( building_pars(ind_alb_win,building_type)   )

           surf_usm_h%zw(nzb_wall,m)         = building_pars(ind_thick_1,building_type)
           surf_usm_h%zw(nzb_wall+1,m)       = building_pars(ind_thick_2,building_type)
           surf_usm_h%zw(nzb_wall+2,m)       = building_pars(ind_thick_3,building_type)
           surf_usm_h%zw(nzb_wall+3,m)       = building_pars(ind_thick_4,building_type)
           
           surf_usm_h%zw_green(nzb_wall,m)         = building_pars(ind_thick_1,building_type)
           surf_usm_h%zw_green(nzb_wall+1,m)       = building_pars(ind_thick_2,building_type)
           surf_usm_h%zw_green(nzb_wall+2,m)       = building_pars(ind_thick_3,building_type)
           surf_usm_h%zw_green(nzb_wall+3,m)       = building_pars(ind_thick_4,building_type)
           
           surf_usm_h%zw_window(nzb_wall,m)         = building_pars(ind_thick_1,building_type)
           surf_usm_h%zw_window(nzb_wall+1,m)       = building_pars(ind_thick_2,building_type)
           surf_usm_h%zw_window(nzb_wall+2,m)       = building_pars(ind_thick_3,building_type)
           surf_usm_h%zw_window(nzb_wall+3,m)       = building_pars(ind_thick_4,building_type)

           surf_usm_h%c_surface(m)           = building_pars(45,building_type)  
           surf_usm_h%lambda_surf(m)         = building_pars(46,building_type)  
           surf_usm_h%c_surface_green(m)     = building_pars(45,building_type)  
           surf_usm_h%lambda_surf_green(m)   = building_pars(46,building_type)  
           surf_usm_h%c_surface_window(m)    = building_pars(45,building_type)  
           surf_usm_h%lambda_surf_window(m)  = building_pars(46,building_type)  

        ENDDO

        DO  l = 0, 3
           DO  m = 1, surf_usm_v(l)%ns

              surf_usm_v(l)%surface_types(m) = wall_category         !< default category for root surface
!
!--           In order to distinguish between ground floor level and 
!--           above-ground-floor level surfaces, set input indices. 
              ind_wall_frac    = MERGE( ind_wall_frac_gfl,    ind_wall_frac_agfl,    &
                                        surf_usm_v(l)%ground_level(m) )
              ind_win_frac     = MERGE( ind_win_frac_gfl,     ind_win_frac_agfl,     &
                                        surf_usm_v(l)%ground_level(m) )
              ind_green_frac_w = MERGE( ind_green_frac_w_gfl, ind_green_frac_w_agfl, &
                                        surf_usm_v(l)%ground_level(m) )
              ind_green_frac_r = MERGE( ind_green_frac_r_gfl, ind_green_frac_r_agfl, &
                                        surf_usm_v(l)%ground_level(m) )
              ind_lai_r        = MERGE( ind_lai_r_gfl,        ind_lai_r_agfl,        &
                                        surf_usm_v(l)%ground_level(m) )
              ind_lai_w        = MERGE( ind_lai_w_gfl,        ind_lai_w_agfl,        &
                                        surf_usm_v(l)%ground_level(m) )
              ind_hc1          = MERGE( ind_hc1_gfl,          ind_hc1_agfl,          &
                                        surf_usm_v(l)%ground_level(m) )
              ind_hc2          = MERGE( ind_hc2_gfl,          ind_hc2_agfl,          &
                                        surf_usm_v(l)%ground_level(m) )
              ind_hc3          = MERGE( ind_hc3_gfl,          ind_hc3_agfl,          &
                                        surf_usm_v(l)%ground_level(m) )
              ind_tc1          = MERGE( ind_tc1_gfl,          ind_tc1_agfl,          &
                                        surf_usm_v(l)%ground_level(m) )
              ind_tc2          = MERGE( ind_tc2_gfl,          ind_tc2_agfl,          &
                                        surf_usm_v(l)%ground_level(m) )
              ind_tc3          = MERGE( ind_tc3_gfl,          ind_tc3_agfl,          &
                                        surf_usm_v(l)%ground_level(m) )
              ind_emis_wall    = MERGE( ind_emis_wall_gfl,    ind_emis_wall_agfl,    &
                                        surf_usm_v(l)%ground_level(m) )
              ind_emis_green   = MERGE( ind_emis_green_gfl,   ind_emis_green_agfl,   &
                                        surf_usm_v(l)%ground_level(m) )
              ind_emis_win     = MERGE( ind_emis_win_gfl,     ind_emis_win_agfl,     &
                                        surf_usm_v(l)%ground_level(m) )
              ind_trans        = MERGE( ind_trans_gfl,       ind_trans_agfl,         &
                                        surf_usm_v(l)%ground_level(m) )
              ind_z0           = MERGE( ind_z0_gfl,           ind_z0_agfl,           &
                                        surf_usm_v(l)%ground_level(m) )
              ind_z0qh         = MERGE( ind_z0qh_gfl,         ind_z0qh_agfl,         &
                                        surf_usm_v(l)%ground_level(m) )

!
!--           Initialize relatvie wall- (0), green- (1) and window (2) fractions
              surf_usm_v(l)%frac(ind_veg_wall,m)   = building_pars(ind_wall_frac,building_type)   
              surf_usm_v(l)%frac(ind_pav_green,m)  = building_pars(ind_green_frac_w,building_type) 
              surf_usm_v(l)%frac(ind_wat_win,m)    = building_pars(ind_win_frac,building_type)  
              surf_usm_v(l)%lai(m)                 = building_pars(ind_lai_w,building_type)  

              surf_usm_v(l)%rho_c_wall(nzb_wall,m)   = building_pars(ind_hc1,building_type)  
              surf_usm_v(l)%rho_c_wall(nzb_wall+1,m) = building_pars(ind_hc1,building_type)
              surf_usm_v(l)%rho_c_wall(nzb_wall+2,m) = building_pars(ind_hc2,building_type)
              surf_usm_v(l)%rho_c_wall(nzb_wall+3,m) = building_pars(ind_hc3,building_type)    
              
              surf_usm_v(l)%rho_c_green(nzb_wall,m)   = building_pars(ind_hc1,building_type)  
              surf_usm_v(l)%rho_c_green(nzb_wall+1,m) = building_pars(ind_hc1,building_type)
              surf_usm_v(l)%rho_c_green(nzb_wall+2,m) = building_pars(ind_hc2,building_type)
              surf_usm_v(l)%rho_c_green(nzb_wall+3,m) = building_pars(ind_hc3,building_type)    
              
              surf_usm_v(l)%rho_c_window(nzb_wall,m)   = building_pars(ind_hc1,building_type)  
              surf_usm_v(l)%rho_c_window(nzb_wall+1,m) = building_pars(ind_hc1,building_type)
              surf_usm_v(l)%rho_c_window(nzb_wall+2,m) = building_pars(ind_hc2,building_type)
              surf_usm_v(l)%rho_c_window(nzb_wall+3,m) = building_pars(ind_hc3,building_type)    

              surf_usm_v(l)%lambda_h(nzb_wall,m)   = building_pars(ind_tc1,building_type)  
              surf_usm_v(l)%lambda_h(nzb_wall+1,m) = building_pars(ind_tc1,building_type) 
              surf_usm_v(l)%lambda_h(nzb_wall+2,m) = building_pars(ind_tc2,building_type)
              surf_usm_v(l)%lambda_h(nzb_wall+3,m) = building_pars(ind_tc3,building_type)    
              
              surf_usm_v(l)%lambda_h_green(nzb_wall,m)   = building_pars(ind_tc1,building_type)  
              surf_usm_v(l)%lambda_h_green(nzb_wall+1,m) = building_pars(ind_tc1,building_type) 
              surf_usm_v(l)%lambda_h_green(nzb_wall+2,m) = building_pars(ind_tc2,building_type)
              surf_usm_v(l)%lambda_h_green(nzb_wall+3,m) = building_pars(ind_tc3,building_type)    

              surf_usm_v(l)%lambda_h_window(nzb_wall,m)   = building_pars(ind_tc1,building_type)  
              surf_usm_v(l)%lambda_h_window(nzb_wall+1,m) = building_pars(ind_tc1,building_type) 
              surf_usm_v(l)%lambda_h_window(nzb_wall+2,m) = building_pars(ind_tc2,building_type)
              surf_usm_v(l)%lambda_h_window(nzb_wall+3,m) = building_pars(ind_tc3,building_type)    

              surf_usm_v(l)%target_temp_summer(m)  = building_pars(12,building_type)    
              surf_usm_v(l)%target_temp_winter(m)  = building_pars(13,building_type)    
!
!--           emissivity of wall-, green- and window fraction 
              surf_usm_v(l)%emissivity(ind_veg_wall,m)  = building_pars(ind_emis_wall,building_type)
              surf_usm_v(l)%emissivity(ind_pav_green,m) = building_pars(ind_emis_green,building_type)
              surf_usm_v(l)%emissivity(ind_wat_win,m)   = building_pars(ind_emis_win,building_type)

              surf_usm_v(l)%transmissivity(m)      = building_pars(ind_trans,building_type)

              surf_usm_v(l)%z0(m)                  = building_pars(ind_z0,building_type)
              surf_usm_v(l)%z0h(m)                 = building_pars(ind_z0qh,building_type)
              surf_usm_v(l)%z0q(m)                 = building_pars(ind_z0qh,building_type)

              surf_usm_v(l)%albedo_type(ind_veg_wall,m)  = INT( building_pars(ind_alb_wall,building_type) )
              surf_usm_v(l)%albedo_type(ind_pav_green,m) = INT( building_pars(ind_alb_green,building_type) )
              surf_usm_v(l)%albedo_type(ind_wat_win,m)   = INT( building_pars(ind_alb_win,building_type) )

              surf_usm_v(l)%zw(nzb_wall,m)         = building_pars(ind_thick_1,building_type)
              surf_usm_v(l)%zw(nzb_wall+1,m)       = building_pars(ind_thick_2,building_type)
              surf_usm_v(l)%zw(nzb_wall+2,m)       = building_pars(ind_thick_3,building_type)
              surf_usm_v(l)%zw(nzb_wall+3,m)       = building_pars(ind_thick_4,building_type)
              
              surf_usm_v(l)%zw_green(nzb_wall,m)         = building_pars(ind_thick_1,building_type)
              surf_usm_v(l)%zw_green(nzb_wall+1,m)       = building_pars(ind_thick_2,building_type)
              surf_usm_v(l)%zw_green(nzb_wall+2,m)       = building_pars(ind_thick_3,building_type)
              surf_usm_v(l)%zw_green(nzb_wall+3,m)       = building_pars(ind_thick_4,building_type)

              surf_usm_v(l)%zw_window(nzb_wall,m)         = building_pars(ind_thick_1,building_type)
              surf_usm_v(l)%zw_window(nzb_wall+1,m)       = building_pars(ind_thick_2,building_type)
              surf_usm_v(l)%zw_window(nzb_wall+2,m)       = building_pars(ind_thick_3,building_type)
              surf_usm_v(l)%zw_window(nzb_wall+3,m)       = building_pars(ind_thick_4,building_type)

              surf_usm_v(l)%c_surface(m)           = building_pars(45,building_type)  
              surf_usm_v(l)%lambda_surf(m)         = building_pars(46,building_type)
              surf_usm_v(l)%c_surface_green(m)     = building_pars(45,building_type)  
              surf_usm_v(l)%lambda_surf_green(m)   = building_pars(46,building_type)
              surf_usm_v(l)%c_surface_window(m)    = building_pars(45,building_type)  
              surf_usm_v(l)%lambda_surf_window(m)  = building_pars(46,building_type)

           ENDDO
        ENDDO
!
!--     Level 2 - initialization via building type read from file
        IF ( building_type_f%from_file )  THEN
           DO  m = 1, surf_usm_h%ns
              i = surf_usm_h%i(m)
              j = surf_usm_h%j(m)
!
!--           For the moment, limit building type to 6 (to overcome errors in input file).
              st = building_type_f%var(j,i)
              IF ( st /= building_type_f%fill )  THEN

!
!--              In order to distinguish between ground floor level and 
!--              above-ground-floor level surfaces, set input indices. 
                 ind_wall_frac    = MERGE( ind_wall_frac_gfl,    ind_wall_frac_agfl,    &
                                           surf_usm_h%ground_level(m) )
                 ind_win_frac     = MERGE( ind_win_frac_gfl,     ind_win_frac_agfl,     &
                                           surf_usm_h%ground_level(m) )
                 ind_green_frac_w = MERGE( ind_green_frac_w_gfl, ind_green_frac_w_agfl, &
                                           surf_usm_h%ground_level(m) )
                 ind_green_frac_r = MERGE( ind_green_frac_r_gfl, ind_green_frac_r_agfl, &
                                           surf_usm_h%ground_level(m) )
                 ind_lai_r        = MERGE( ind_lai_r_gfl,        ind_lai_r_agfl,        &
                                           surf_usm_h%ground_level(m) )
                 ind_lai_w        = MERGE( ind_lai_w_gfl,        ind_lai_w_agfl,        &
                                           surf_usm_h%ground_level(m) )
                 ind_hc1          = MERGE( ind_hc1_gfl,          ind_hc1_agfl,          &
                                           surf_usm_h%ground_level(m) )
                 ind_hc2          = MERGE( ind_hc2_gfl,          ind_hc2_agfl,          &
                                           surf_usm_h%ground_level(m) )
                 ind_hc3          = MERGE( ind_hc3_gfl,          ind_hc3_agfl,          &
                                           surf_usm_h%ground_level(m) )
                 ind_tc1          = MERGE( ind_tc1_gfl,          ind_tc1_agfl,          &
                                           surf_usm_h%ground_level(m) )
                 ind_tc2          = MERGE( ind_tc2_gfl,          ind_tc2_agfl,          &
                                           surf_usm_h%ground_level(m) )
                 ind_tc3          = MERGE( ind_tc3_gfl,          ind_tc3_agfl,          &
                                           surf_usm_h%ground_level(m) )
                 ind_emis_wall    = MERGE( ind_emis_wall_gfl,    ind_emis_wall_agfl,    &
                                           surf_usm_h%ground_level(m) )
                 ind_emis_green   = MERGE( ind_emis_green_gfl,   ind_emis_green_agfl,   &
                                           surf_usm_h%ground_level(m) )
                 ind_emis_win     = MERGE( ind_emis_win_gfl,     ind_emis_win_agfl,     &
                                           surf_usm_h%ground_level(m) )
                 ind_trans        = MERGE( ind_trans_gfl,        ind_trans_agfl,        &
                                           surf_usm_h%ground_level(m) )
                 ind_z0           = MERGE( ind_z0_gfl,           ind_z0_agfl,           &
                                           surf_usm_h%ground_level(m) )
                 ind_z0qh         = MERGE( ind_z0qh_gfl,         ind_z0qh_agfl,         &
                                           surf_usm_h%ground_level(m) )

!
!--              Initialize relatvie wall- (0), green- (1) and window (2) fractions
                 surf_usm_h%frac(ind_veg_wall,m)  = building_pars(ind_wall_frac,st)   
                 surf_usm_h%frac(ind_pav_green,m) = building_pars(ind_green_frac_r,st)  
                 surf_usm_h%frac(ind_wat_win,m)   = building_pars(ind_win_frac,st)  
                 surf_usm_h%lai(m)                = building_pars(ind_green_frac_r,st)  

                 surf_usm_h%rho_c_wall(nzb_wall,m)   = building_pars(ind_hc1,st)  
                 surf_usm_h%rho_c_wall(nzb_wall+1,m) = building_pars(ind_hc1,st)
                 surf_usm_h%rho_c_wall(nzb_wall+2,m) = building_pars(ind_hc2,st)
                 surf_usm_h%rho_c_wall(nzb_wall+3,m) = building_pars(ind_hc3,st)    
                 surf_usm_h%lambda_h(nzb_wall,m)   = building_pars(ind_tc1,st)  
                 surf_usm_h%lambda_h(nzb_wall+1,m) = building_pars(ind_tc1,st) 
                 surf_usm_h%lambda_h(nzb_wall+2,m) = building_pars(ind_tc2,st)
                 surf_usm_h%lambda_h(nzb_wall+3,m) = building_pars(ind_tc3,st)    
                 
                 surf_usm_h%rho_c_green(nzb_wall,m)   = building_pars(ind_hc1,st)  
                 surf_usm_h%rho_c_green(nzb_wall+1,m) = building_pars(ind_hc1,st)
                 surf_usm_h%rho_c_green(nzb_wall+2,m) = building_pars(ind_hc2,st)
                 surf_usm_h%rho_c_green(nzb_wall+3,m) = building_pars(ind_hc3,st)    
                 surf_usm_h%lambda_h_green(nzb_wall,m)   = building_pars(ind_tc1,st)  
                 surf_usm_h%lambda_h_green(nzb_wall+1,m) = building_pars(ind_tc1,st) 
                 surf_usm_h%lambda_h_green(nzb_wall+2,m) = building_pars(ind_tc2,st)
                 surf_usm_h%lambda_h_green(nzb_wall+3,m) = building_pars(ind_tc3,st)    
                
                 surf_usm_h%rho_c_window(nzb_wall,m)   = building_pars(ind_hc1,st)  
                 surf_usm_h%rho_c_window(nzb_wall+1,m) = building_pars(ind_hc1,st)
                 surf_usm_h%rho_c_window(nzb_wall+2,m) = building_pars(ind_hc2,st)
                 surf_usm_h%rho_c_window(nzb_wall+3,m) = building_pars(ind_hc3,st)    
                 surf_usm_h%lambda_h_window(nzb_wall,m)   = building_pars(ind_tc1,st)  
                 surf_usm_h%lambda_h_window(nzb_wall+1,m) = building_pars(ind_tc1,st) 
                 surf_usm_h%lambda_h_window(nzb_wall+2,m) = building_pars(ind_tc2,st)
                 surf_usm_h%lambda_h_window(nzb_wall+3,m) = building_pars(ind_tc3,st)    

                 surf_usm_h%target_temp_summer(m)  = building_pars(12,st)    
                 surf_usm_h%target_temp_winter(m)  = building_pars(13,st)    
!
!--              emissivity of wall-, green- and window fraction 
                 surf_usm_h%emissivity(ind_veg_wall,m)  = building_pars(ind_emis_wall,st)
                 surf_usm_h%emissivity(ind_pav_green,m) = building_pars(ind_emis_green,st)
                 surf_usm_h%emissivity(ind_wat_win,m)   = building_pars(ind_emis_win,st)

                 surf_usm_h%transmissivity(m)      = building_pars(ind_trans,st)

                 surf_usm_h%z0(m)                  = building_pars(ind_z0,st)
                 surf_usm_h%z0h(m)                 = building_pars(ind_z0qh,st)
                 surf_usm_h%z0q(m)                 = building_pars(ind_z0qh,st)
!
!--              albedo type for wall fraction, green fraction, window fraction
                 surf_usm_h%albedo_type(ind_veg_wall,m)  = INT( building_pars(ind_alb_wall,st) )
                 surf_usm_h%albedo_type(ind_pav_green,m) = INT( building_pars(ind_alb_green,st) )
                 surf_usm_h%albedo_type(ind_wat_win,m)   = INT( building_pars(ind_alb_win,st) )

                 surf_usm_h%zw(nzb_wall,m)         = building_pars(ind_thick_1,st)
                 surf_usm_h%zw(nzb_wall+1,m)       = building_pars(ind_thick_2,st)
                 surf_usm_h%zw(nzb_wall+2,m)       = building_pars(ind_thick_3,st)
                 surf_usm_h%zw(nzb_wall+3,m)       = building_pars(ind_thick_4,st)
                 
                 surf_usm_h%zw_green(nzb_wall,m)         = building_pars(ind_thick_1,st)
                 surf_usm_h%zw_green(nzb_wall+1,m)       = building_pars(ind_thick_2,st)
                 surf_usm_h%zw_green(nzb_wall+2,m)       = building_pars(ind_thick_3,st)
                 surf_usm_h%zw_green(nzb_wall+3,m)       = building_pars(ind_thick_4,st)

                 surf_usm_h%zw_window(nzb_wall,m)         = building_pars(ind_thick_1,st)
                 surf_usm_h%zw_window(nzb_wall+1,m)       = building_pars(ind_thick_2,st)
                 surf_usm_h%zw_window(nzb_wall+2,m)       = building_pars(ind_thick_3,st)
                 surf_usm_h%zw_window(nzb_wall+3,m)       = building_pars(ind_thick_4,st)

                 surf_usm_h%c_surface(m)           = building_pars(45,st)  
                 surf_usm_h%lambda_surf(m)         = building_pars(46,st)
                 surf_usm_h%c_surface_green(m)     = building_pars(45,st)  
                 surf_usm_h%lambda_surf_green(m)   = building_pars(46,st)
                 surf_usm_h%c_surface_window(m)    = building_pars(45,st)  
                 surf_usm_h%lambda_surf_window(m)  = building_pars(46,st)

              ENDIF
           ENDDO

           DO  l = 0, 3
              DO  m = 1, surf_usm_v(l)%ns
                 i = surf_usm_v(l)%i(m) + surf_usm_v(l)%ioff
                 j = surf_usm_v(l)%j(m) + surf_usm_v(l)%joff
!
!--              For the moment, limit building type to 6 (to overcome errors in input file).

                 st = building_type_f%var(j,i)
                 IF ( st /= building_type_f%fill )  THEN

!
!--                 In order to distinguish between ground floor level and 
!--                 above-ground-floor level surfaces, set input indices. 
                    ind_wall_frac    = MERGE( ind_wall_frac_gfl,    ind_wall_frac_agfl,    &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_win_frac     = MERGE( ind_win_frac_gfl,     ind_win_frac_agfl,     &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_green_frac_w = MERGE( ind_green_frac_w_gfl, ind_green_frac_w_agfl, &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_green_frac_r = MERGE( ind_green_frac_r_gfl, ind_green_frac_r_agfl, &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_lai_r        = MERGE( ind_lai_r_gfl,        ind_lai_r_agfl,        &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_lai_w        = MERGE( ind_lai_w_gfl,        ind_lai_w_agfl,        &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_hc1          = MERGE( ind_hc1_gfl,          ind_hc1_agfl,          &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_hc2          = MERGE( ind_hc2_gfl,          ind_hc2_agfl,          &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_hc3          = MERGE( ind_hc3_gfl,          ind_hc3_agfl,          &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_tc1          = MERGE( ind_tc1_gfl,          ind_tc1_agfl,          &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_tc2          = MERGE( ind_tc2_gfl,          ind_tc2_agfl,          &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_tc3          = MERGE( ind_tc3_gfl,          ind_tc3_agfl,          &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_emis_wall    = MERGE( ind_emis_wall_gfl,    ind_emis_wall_agfl,    &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_emis_green   = MERGE( ind_emis_green_gfl,   ind_emis_green_agfl,   &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_emis_win     = MERGE( ind_emis_win_gfl,     ind_emis_win_agfl,     &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_trans        = MERGE( ind_trans_gfl,       ind_trans_agfl,         &
                                        surf_usm_v(l)%ground_level(m) )
                    ind_z0           = MERGE( ind_z0_gfl,           ind_z0_agfl,           &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_z0qh         = MERGE( ind_z0qh_gfl,         ind_z0qh_agfl,         &
                                              surf_usm_v(l)%ground_level(m) )

!
!--                 Initialize relatvie wall- (0), green- (1) and window (2) fractions
                    surf_usm_v(l)%frac(ind_veg_wall,m)  = building_pars(ind_wall_frac,st)   
                    surf_usm_v(l)%frac(ind_pav_green,m) = building_pars(ind_green_frac_w,st)  
                    surf_usm_v(l)%frac(ind_wat_win,m)   = building_pars(ind_win_frac,st)   
                    surf_usm_v(l)%lai(m)                = building_pars(ind_lai_w,st)  

                    surf_usm_v(l)%rho_c_wall(nzb_wall,m)   = building_pars(ind_hc1,st)  
                    surf_usm_v(l)%rho_c_wall(nzb_wall+1,m) = building_pars(ind_hc1,st)
                    surf_usm_v(l)%rho_c_wall(nzb_wall+2,m) = building_pars(ind_hc2,st)
                    surf_usm_v(l)%rho_c_wall(nzb_wall+3,m) = building_pars(ind_hc3,st)
                    
                    surf_usm_v(l)%rho_c_green(nzb_wall,m)   = building_pars(ind_hc1,st)  
                    surf_usm_v(l)%rho_c_green(nzb_wall+1,m) = building_pars(ind_hc1,st)
                    surf_usm_v(l)%rho_c_green(nzb_wall+2,m) = building_pars(ind_hc2,st)
                    surf_usm_v(l)%rho_c_green(nzb_wall+3,m) = building_pars(ind_hc3,st)
                    
                    surf_usm_v(l)%rho_c_window(nzb_wall,m)   = building_pars(ind_hc1,st)  
                    surf_usm_v(l)%rho_c_window(nzb_wall+1,m) = building_pars(ind_hc1,st)
                    surf_usm_v(l)%rho_c_window(nzb_wall+2,m) = building_pars(ind_hc2,st)
                    surf_usm_v(l)%rho_c_window(nzb_wall+3,m) = building_pars(ind_hc3,st)

                    surf_usm_v(l)%lambda_h(nzb_wall,m)   = building_pars(ind_tc1,st)  
                    surf_usm_v(l)%lambda_h(nzb_wall+1,m) = building_pars(ind_tc1,st) 
                    surf_usm_v(l)%lambda_h(nzb_wall+2,m) = building_pars(ind_tc2,st)
                    surf_usm_v(l)%lambda_h(nzb_wall+3,m) = building_pars(ind_tc3,st) 
                    
                    surf_usm_v(l)%lambda_h_green(nzb_wall,m)   = building_pars(ind_tc1,st)  
                    surf_usm_v(l)%lambda_h_green(nzb_wall+1,m) = building_pars(ind_tc1,st) 
                    surf_usm_v(l)%lambda_h_green(nzb_wall+2,m) = building_pars(ind_tc2,st)
                    surf_usm_v(l)%lambda_h_green(nzb_wall+3,m) = building_pars(ind_tc3,st) 
                    
                    surf_usm_v(l)%lambda_h_window(nzb_wall,m)   = building_pars(ind_tc1,st)  
                    surf_usm_v(l)%lambda_h_window(nzb_wall+1,m) = building_pars(ind_tc1,st) 
                    surf_usm_v(l)%lambda_h_window(nzb_wall+2,m) = building_pars(ind_tc2,st)
                    surf_usm_v(l)%lambda_h_window(nzb_wall+3,m) = building_pars(ind_tc3,st) 

                    surf_usm_v(l)%target_temp_summer(m)  = building_pars(12,st)    
                    surf_usm_v(l)%target_temp_winter(m)  = building_pars(13,st)    
!
!--                 emissivity of wall-, green- and window fraction 
                    surf_usm_v(l)%emissivity(ind_veg_wall,m)  = building_pars(ind_emis_wall,st)
                    surf_usm_v(l)%emissivity(ind_pav_green,m) = building_pars(ind_emis_green,st)
                    surf_usm_v(l)%emissivity(ind_wat_win,m)   = building_pars(ind_emis_win,st)

                    surf_usm_v(l)%transmissivity(m)      = building_pars(ind_trans,st)

                    surf_usm_v(l)%z0(m)                  = building_pars(ind_z0,st)
                    surf_usm_v(l)%z0h(m)                 = building_pars(ind_z0qh,st)
                    surf_usm_v(l)%z0q(m)                 = building_pars(ind_z0qh,st)

                    surf_usm_v(l)%albedo_type(ind_veg_wall,m)  = INT( building_pars(ind_alb_wall,st) )
                    surf_usm_v(l)%albedo_type(ind_pav_green,m) = INT( building_pars(ind_alb_green,st) )
                    surf_usm_v(l)%albedo_type(ind_wat_win,m)   = INT( building_pars(ind_alb_win,st) )

                    surf_usm_v(l)%zw(nzb_wall,m)         = building_pars(ind_thick_1,st)
                    surf_usm_v(l)%zw(nzb_wall+1,m)       = building_pars(ind_thick_2,st)
                    surf_usm_v(l)%zw(nzb_wall+2,m)       = building_pars(ind_thick_3,st)
                    surf_usm_v(l)%zw(nzb_wall+3,m)       = building_pars(ind_thick_4,st)
                    
                    surf_usm_v(l)%zw_green(nzb_wall,m)         = building_pars(ind_thick_1,st)
                    surf_usm_v(l)%zw_green(nzb_wall+1,m)       = building_pars(ind_thick_2,st)
                    surf_usm_v(l)%zw_green(nzb_wall+2,m)       = building_pars(ind_thick_3,st)
                    surf_usm_v(l)%zw_green(nzb_wall+3,m)       = building_pars(ind_thick_4,st)
                    
                    surf_usm_v(l)%zw_window(nzb_wall,m)         = building_pars(ind_thick_1,st)
                    surf_usm_v(l)%zw_window(nzb_wall+1,m)       = building_pars(ind_thick_2,st)
                    surf_usm_v(l)%zw_window(nzb_wall+2,m)       = building_pars(ind_thick_3,st)
                    surf_usm_v(l)%zw_window(nzb_wall+3,m)       = building_pars(ind_thick_4,st)

                    surf_usm_v(l)%c_surface(m)           = building_pars(45,st)  
                    surf_usm_v(l)%lambda_surf(m)         = building_pars(46,st) 
                    surf_usm_v(l)%c_surface_green(m)     = building_pars(45,st)  
                    surf_usm_v(l)%lambda_surf_green(m)   = building_pars(46,st) 
                    surf_usm_v(l)%c_surface_window(m)    = building_pars(45,st)  
                    surf_usm_v(l)%lambda_surf_window(m)  = building_pars(46,st) 


                 ENDIF
              ENDDO
           ENDDO
        ENDIF 

!
!--     Level 3 - initialization via building_pars read from file
        IF ( building_pars_f%from_file )  THEN
           DO  m = 1, surf_usm_h%ns
              i = surf_usm_h%i(m)
              j = surf_usm_h%j(m)

!
!--           In order to distinguish between ground floor level and 
!--           above-ground-floor level surfaces, set input indices. 
              ind_wall_frac    = MERGE( ind_wall_frac_gfl,    ind_wall_frac_agfl,    &
                                        surf_usm_h%ground_level(m) )
              ind_win_frac     = MERGE( ind_win_frac_gfl,     ind_win_frac_agfl,     &
                                        surf_usm_h%ground_level(m) )
              ind_green_frac_w = MERGE( ind_green_frac_w_gfl, ind_green_frac_w_agfl, &
                                        surf_usm_h%ground_level(m) )
              ind_green_frac_r = MERGE( ind_green_frac_r_gfl, ind_green_frac_r_agfl, &
                                        surf_usm_h%ground_level(m) )
              ind_lai_r        = MERGE( ind_lai_r_gfl,        ind_lai_r_agfl,        &
                                        surf_usm_h%ground_level(m) )
              ind_lai_w        = MERGE( ind_lai_w_gfl,        ind_lai_w_agfl,        &
                                        surf_usm_h%ground_level(m) )
              ind_hc1          = MERGE( ind_hc1_gfl,          ind_hc1_agfl,          &
                                        surf_usm_h%ground_level(m) )
              ind_hc2          = MERGE( ind_hc2_gfl,          ind_hc2_agfl,          &
                                        surf_usm_h%ground_level(m) )
              ind_hc3          = MERGE( ind_hc3_gfl,          ind_hc3_agfl,          &
                                        surf_usm_h%ground_level(m) )
              ind_tc1          = MERGE( ind_tc1_gfl,          ind_tc1_agfl,          &
                                        surf_usm_h%ground_level(m) )
              ind_tc2          = MERGE( ind_tc2_gfl,          ind_tc2_agfl,          &
                                        surf_usm_h%ground_level(m) )
              ind_tc3          = MERGE( ind_tc3_gfl,          ind_tc3_agfl,          &
                                        surf_usm_h%ground_level(m) )
              ind_emis_wall    = MERGE( ind_emis_wall_gfl,    ind_emis_wall_agfl,    &
                                        surf_usm_h%ground_level(m) )
              ind_emis_green   = MERGE( ind_emis_green_gfl,   ind_emis_green_agfl,   &
                                        surf_usm_h%ground_level(m) )
              ind_emis_win     = MERGE( ind_emis_win_gfl,     ind_emis_win_agfl,     &
                                        surf_usm_h%ground_level(m) )
              ind_trans        = MERGE( ind_trans_gfl,        ind_trans_agfl,        &
                                        surf_usm_h%ground_level(m) )
              ind_z0           = MERGE( ind_z0_gfl,           ind_z0_agfl,           &
                                        surf_usm_h%ground_level(m) )
              ind_z0qh         = MERGE( ind_z0qh_gfl,         ind_z0qh_agfl,         &
                                        surf_usm_h%ground_level(m) )

!
!--           Initialize relatvie wall- (0), green- (1) and window (2) fractions
              IF ( building_pars_f%pars_xy(ind_wall_frac,j,i) /= building_pars_f%fill )    &
                 surf_usm_h%frac(ind_veg_wall,m)  = building_pars_f%pars_xy(ind_wall_frac,j,i)   
              IF ( building_pars_f%pars_xy(ind_green_frac_r,j,i) /= building_pars_f%fill ) & 
                 surf_usm_h%frac(ind_pav_green,m) = building_pars_f%pars_xy(ind_green_frac_r,j,i) 
              IF ( building_pars_f%pars_xy(ind_win_frac,j,i) /= building_pars_f%fill )     & 
                 surf_usm_h%frac(ind_wat_win,m)   = building_pars_f%pars_xy(ind_win_frac,j,i)

 
              IF ( building_pars_f%pars_xy(ind_lai_r,j,i) /= building_pars_f%fill )        &
                 surf_usm_h%lai(m)             = building_pars_f%pars_xy(ind_lai_r,j,i)

              IF ( building_pars_f%pars_xy(ind_hc1,j,i) /= building_pars_f%fill )  THEN 
                 surf_usm_h%rho_c_wall(nzb_wall,m)   = building_pars_f%pars_xy(ind_hc1,j,i)  
                 surf_usm_h%rho_c_wall(nzb_wall+1,m) = building_pars_f%pars_xy(ind_hc1,j,i)
              ENDIF
              IF ( building_pars_f%pars_xy(ind_hc2,j,i) /= building_pars_f%fill )    &
                 surf_usm_h%rho_c_wall(nzb_wall+2,m) = building_pars_f%pars_xy(ind_hc2,j,i)
              IF ( building_pars_f%pars_xy(ind_hc3,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%rho_c_wall(nzb_wall+3,m) = building_pars_f%pars_xy(ind_hc3,j,i)
              IF ( building_pars_f%pars_xy(ind_hc1,j,i) /= building_pars_f%fill )  THEN 
                 surf_usm_h%rho_c_green(nzb_wall,m)   = building_pars_f%pars_xy(ind_hc1,j,i)  
                 surf_usm_h%rho_c_green(nzb_wall+1,m) = building_pars_f%pars_xy(ind_hc1,j,i)
              ENDIF
              IF ( building_pars_f%pars_xy(ind_hc2,j,i) /= building_pars_f%fill )    &
                 surf_usm_h%rho_c_green(nzb_wall+2,m) = building_pars_f%pars_xy(ind_hc2,j,i)
              IF ( building_pars_f%pars_xy(ind_hc3,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%rho_c_green(nzb_wall+3,m) = building_pars_f%pars_xy(ind_hc3,j,i)
              IF ( building_pars_f%pars_xy(ind_hc1,j,i) /= building_pars_f%fill )  THEN 
                 surf_usm_h%rho_c_window(nzb_wall,m)   = building_pars_f%pars_xy(ind_hc1,j,i)  
                 surf_usm_h%rho_c_window(nzb_wall+1,m) = building_pars_f%pars_xy(ind_hc1,j,i)
              ENDIF
              IF ( building_pars_f%pars_xy(ind_hc2,j,i) /= building_pars_f%fill )    &
                 surf_usm_h%rho_c_window(nzb_wall+2,m) = building_pars_f%pars_xy(ind_hc2,j,i)
              IF ( building_pars_f%pars_xy(ind_hc3,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%rho_c_window(nzb_wall+3,m) = building_pars_f%pars_xy(ind_hc3,j,i)

              IF ( building_pars_f%pars_xy(ind_tc1,j,i) /= building_pars_f%fill )  THEN
                 surf_usm_h%lambda_h(nzb_wall,m)   = building_pars_f%pars_xy(ind_tc1,j,i)         
                 surf_usm_h%lambda_h(nzb_wall+1,m) = building_pars_f%pars_xy(ind_tc1,j,i)        
              ENDIF
              IF ( building_pars_f%pars_xy(ind_tc2,j,i) /= building_pars_f%fill )    &
                 surf_usm_h%lambda_h(nzb_wall+2,m) = building_pars_f%pars_xy(ind_tc2,j,i)
              IF ( building_pars_f%pars_xy(ind_tc3,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%lambda_h(nzb_wall+3,m) = building_pars_f%pars_xy(ind_tc3,j,i)    
              IF ( building_pars_f%pars_xy(ind_tc1,j,i) /= building_pars_f%fill )  THEN
                 surf_usm_h%lambda_h_green(nzb_wall,m)   = building_pars_f%pars_xy(ind_tc1,j,i)         
                 surf_usm_h%lambda_h_green(nzb_wall+1,m) = building_pars_f%pars_xy(ind_tc1,j,i)        
              ENDIF
              IF ( building_pars_f%pars_xy(ind_tc2,j,i) /= building_pars_f%fill )    &
                 surf_usm_h%lambda_h_green(nzb_wall+2,m) = building_pars_f%pars_xy(ind_tc2,j,i)
              IF ( building_pars_f%pars_xy(ind_tc3,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%lambda_h_green(nzb_wall+3,m) = building_pars_f%pars_xy(ind_tc3,j,i)    
              IF ( building_pars_f%pars_xy(ind_tc1,j,i) /= building_pars_f%fill )  THEN
                 surf_usm_h%lambda_h_window(nzb_wall,m)   = building_pars_f%pars_xy(ind_tc1,j,i)         
                 surf_usm_h%lambda_h_window(nzb_wall+1,m) = building_pars_f%pars_xy(ind_tc1,j,i)        
              ENDIF
              IF ( building_pars_f%pars_xy(ind_tc2,j,i) /= building_pars_f%fill )    &
                 surf_usm_h%lambda_h_window(nzb_wall+2,m) = building_pars_f%pars_xy(ind_tc2,j,i)
              IF ( building_pars_f%pars_xy(ind_tc3,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%lambda_h_window(nzb_wall+3,m) = building_pars_f%pars_xy(ind_tc3,j,i)    

              IF ( building_pars_f%pars_xy(12,j,i) /= building_pars_f%fill )         & 
                 surf_usm_h%target_temp_summer(m)  = building_pars_f%pars_xy(12,j,i)   
              IF ( building_pars_f%pars_xy(13,j,i) /= building_pars_f%fill )         & 
                 surf_usm_h%target_temp_winter(m)  = building_pars_f%pars_xy(13,j,i)   

              IF ( building_pars_f%pars_xy(ind_emis_wall,j,i) /= building_pars_f%fill ) & 
                 surf_usm_h%emissivity(ind_veg_wall,m)  = building_pars_f%pars_xy(ind_emis_wall,j,i)
              IF ( building_pars_f%pars_xy(ind_emis_green,j,i) /= building_pars_f%fill )& 
                 surf_usm_h%emissivity(ind_pav_green,m) = building_pars_f%pars_xy(ind_emis_green,j,i)
              IF ( building_pars_f%pars_xy(ind_emis_win,j,i) /= building_pars_f%fill )  & 
                 surf_usm_h%emissivity(ind_wat_win,m)   = building_pars_f%pars_xy(ind_emis_win,j,i)

              IF ( building_pars_f%pars_xy(ind_trans,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%transmissivity(m)      = building_pars_f%pars_xy(ind_trans,j,i)

              IF ( building_pars_f%pars_xy(ind_z0,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%z0(m)                  = building_pars_f%pars_xy(ind_z0,j,i)
              IF ( building_pars_f%pars_xy(ind_z0qh,j,i) /= building_pars_f%fill )  & 
                 surf_usm_h%z0h(m)                 = building_pars_f%pars_xy(ind_z0qh,j,i)
              IF ( building_pars_f%pars_xy(ind_z0qh,j,i) /= building_pars_f%fill )  & 
                 surf_usm_h%z0q(m)                 = building_pars_f%pars_xy(ind_z0qh,j,i)

              IF ( building_pars_f%pars_xy(ind_alb_wall,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%albedo_type(ind_veg_wall,m)  = building_pars_f%pars_xy(ind_alb_wall,j,i)
              IF ( building_pars_f%pars_xy(ind_alb_green,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%albedo_type(ind_pav_green,m) = building_pars_f%pars_xy(ind_alb_green,j,i)
              IF ( building_pars_f%pars_xy(ind_alb_win,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%albedo_type(ind_wat_win,m)   = building_pars_f%pars_xy(ind_alb_win,j,i)

              IF ( building_pars_f%pars_xy(ind_thick_1,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%zw(nzb_wall,m) = building_pars_f%pars_xy(ind_thick_1,j,i)
              IF ( building_pars_f%pars_xy(ind_thick_2,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%zw(nzb_wall+1,m) = building_pars_f%pars_xy(ind_thick_2,j,i)
              IF ( building_pars_f%pars_xy(ind_thick_3,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%zw(nzb_wall+2,m) = building_pars_f%pars_xy(ind_thick_3,j,i)
              IF ( building_pars_f%pars_xy(ind_thick_4,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%zw(nzb_wall+3,m) = building_pars_f%pars_xy(ind_thick_4,j,i)
              IF ( building_pars_f%pars_xy(ind_thick_1,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%zw_green(nzb_wall,m) = building_pars_f%pars_xy(ind_thick_1,j,i)
              IF ( building_pars_f%pars_xy(ind_thick_2,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%zw_green(nzb_wall+1,m) = building_pars_f%pars_xy(ind_thick_2,j,i)
              IF ( building_pars_f%pars_xy(ind_thick_3,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%zw_green(nzb_wall+2,m) = building_pars_f%pars_xy(ind_thick_3,j,i)
              IF ( building_pars_f%pars_xy(ind_thick_4,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%zw_green(nzb_wall+3,m) = building_pars_f%pars_xy(ind_thick_4,j,i)
              IF ( building_pars_f%pars_xy(ind_thick_1,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%zw_window(nzb_wall,m) = building_pars_f%pars_xy(ind_thick_1,j,i)
              IF ( building_pars_f%pars_xy(ind_thick_2,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%zw_window(nzb_wall+1,m) = building_pars_f%pars_xy(ind_thick_2,j,i)
              IF ( building_pars_f%pars_xy(ind_thick_3,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%zw_window(nzb_wall+2,m) = building_pars_f%pars_xy(ind_thick_3,j,i)
              IF ( building_pars_f%pars_xy(ind_thick_4,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%zw_window(nzb_wall+3,m) = building_pars_f%pars_xy(ind_thick_4,j,i)

              IF ( building_pars_f%pars_xy(45,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%c_surface(m)           = building_pars_f%pars_xy(45,j,i)
              IF ( building_pars_f%pars_xy(46,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%lambda_surf(m)         = building_pars_f%pars_xy(46,j,i)
              IF ( building_pars_f%pars_xy(45,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%c_surface_green(m)           = building_pars_f%pars_xy(45,j,i)
              IF ( building_pars_f%pars_xy(46,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%lambda_surf_green(m)         = building_pars_f%pars_xy(46,j,i)
              IF ( building_pars_f%pars_xy(45,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%c_surface_window(m)           = building_pars_f%pars_xy(45,j,i)
              IF ( building_pars_f%pars_xy(46,j,i) /= building_pars_f%fill )    & 
                 surf_usm_h%lambda_surf_window(m)         = building_pars_f%pars_xy(46,j,i)
           ENDDO



           DO  l = 0, 3
              DO  m = 1, surf_usm_v(l)%ns
                 i = surf_usm_v(l)%i(m) + surf_usm_v(l)%ioff
                 j = surf_usm_v(l)%j(m) + surf_usm_v(l)%joff
                
!
!--              In order to distinguish between ground floor level and 
!--              above-ground-floor level surfaces, set input indices. 
                 ind_wall_frac    = MERGE( ind_wall_frac_gfl,    ind_wall_frac_agfl,     &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_win_frac     = MERGE( ind_win_frac_gfl,     ind_win_frac_agfl,     &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_green_frac_w = MERGE( ind_green_frac_w_gfl, ind_green_frac_w_agfl, &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_green_frac_r = MERGE( ind_green_frac_r_gfl, ind_green_frac_r_agfl, &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_lai_r        = MERGE( ind_lai_r_gfl,        ind_lai_r_agfl,        &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_lai_w        = MERGE( ind_lai_w_gfl,        ind_lai_w_agfl,        &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_hc1          = MERGE( ind_hc1_gfl,          ind_hc1_agfl,          &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_hc2          = MERGE( ind_hc2_gfl,          ind_hc2_agfl,          &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_hc3          = MERGE( ind_hc3_gfl,          ind_hc3_agfl,          &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_tc1          = MERGE( ind_tc1_gfl,          ind_tc1_agfl,          &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_tc2          = MERGE( ind_tc2_gfl,          ind_tc2_agfl,          &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_tc3          = MERGE( ind_tc3_gfl,          ind_tc3_agfl,          &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_emis_wall    = MERGE( ind_emis_wall_gfl,    ind_emis_wall_agfl,    &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_emis_green   = MERGE( ind_emis_green_gfl,   ind_emis_green_agfl,   &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_emis_win     = MERGE( ind_emis_win_gfl,     ind_emis_win_agfl,     &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_trans        = MERGE( ind_trans_gfl,       ind_trans_agfl,         &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_z0           = MERGE( ind_z0_gfl,           ind_z0_agfl,           &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_z0qh         = MERGE( ind_z0qh_gfl,         ind_z0qh_agfl,         &
                                              surf_usm_v(l)%ground_level(m) )

!
!--              Initialize relatvie wall- (0), green- (1) and window (2) fractions
                 IF ( building_pars_f%pars_xy(ind_wall_frac,j,i) /=                     &
                      building_pars_f%fill )                                            &
                    surf_usm_v(l)%frac(ind_veg_wall,m)  =                               &
                                      building_pars_f%pars_xy(ind_wall_frac,j,i)   
                 IF ( building_pars_f%pars_xy(ind_green_frac_w,j,i) /=                  &
                      building_pars_f%fill )                                            & 
                    surf_usm_v(l)%frac(ind_pav_green,m) =                               &
                                      building_pars_f%pars_xy(ind_green_frac_w,j,i) 
                 IF ( building_pars_f%pars_xy(ind_win_frac,j,i) /=                      &
                      building_pars_f%fill )                                            & 
                    surf_usm_v(l)%frac(ind_wat_win,m)   =                               &
                                      building_pars_f%pars_xy(ind_win_frac,j,i)
  
                 IF ( building_pars_f%pars_xy(ind_lai_w,j,i) /= building_pars_f%fill )  & 
                    surf_usm_v(l)%lai(m) = building_pars_f%pars_xy(ind_lai_w,j,i)

                 IF ( building_pars_f%pars_xy(ind_hc1,j,i) /= building_pars_f%fill )    &
                 THEN 
                    surf_usm_v(l)%rho_c_wall(nzb_wall,m)   =                            &
                                                    building_pars_f%pars_xy(ind_hc1,j,i)  
                    surf_usm_v(l)%rho_c_wall(nzb_wall+1,m) =                            &
                                                    building_pars_f%pars_xy(ind_hc1,j,i)
                 ENDIF
                 IF ( building_pars_f%pars_xy(ind_hc2,j,i) /= building_pars_f%fill )    &
                    surf_usm_v(l)%rho_c_wall(nzb_wall+2,m) =                            &                            
                                                    building_pars_f%pars_xy(ind_hc2,j,i)
                 IF ( building_pars_f%pars_xy(ind_hc3,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%rho_c_wall(nzb_wall+3,m) =                            &
                                                    building_pars_f%pars_xy(ind_hc3,j,i)
                 IF ( building_pars_f%pars_xy(ind_hc1,j,i) /= building_pars_f%fill )  THEN 
                    surf_usm_v(l)%rho_c_green(nzb_wall,m)   =                           &
                                                    building_pars_f%pars_xy(ind_hc1,j,i)  
                    surf_usm_v(l)%rho_c_green(nzb_wall+1,m) =                           &
                                                    building_pars_f%pars_xy(ind_hc1,j,i)
                 ENDIF
                 IF ( building_pars_f%pars_xy(ind_hc2,j,i) /= building_pars_f%fill )    &
                    surf_usm_v(l)%rho_c_green(nzb_wall+2,m) = building_pars_f%pars_xy(ind_hc2,j,i)
                 IF ( building_pars_f%pars_xy(ind_hc3,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%rho_c_green(nzb_wall+3,m) = building_pars_f%pars_xy(ind_hc3,j,i)    
                 IF ( building_pars_f%pars_xy(ind_hc1,j,i) /= building_pars_f%fill )  THEN 
                    surf_usm_v(l)%rho_c_window(nzb_wall,m)   = building_pars_f%pars_xy(ind_hc1,j,i)  
                    surf_usm_v(l)%rho_c_window(nzb_wall+1,m) = building_pars_f%pars_xy(ind_hc1,j,i)
                 ENDIF
                 IF ( building_pars_f%pars_xy(ind_hc2,j,i) /= building_pars_f%fill )    &
                    surf_usm_v(l)%rho_c_window(nzb_wall+2,m) = building_pars_f%pars_xy(ind_hc2,j,i)
                 IF ( building_pars_f%pars_xy(ind_hc3,j,i) /= building_pars_f%fill )    &
                    surf_usm_v(l)%rho_c_window(nzb_wall+3,m) = building_pars_f%pars_xy(ind_hc3,j,i)

                 IF ( building_pars_f%pars_xy(ind_tc1,j,i) /= building_pars_f%fill )  THEN
                    surf_usm_v(l)%lambda_h(nzb_wall,m)   = building_pars_f%pars_xy(ind_tc1,j,i)         
                    surf_usm_v(l)%lambda_h(nzb_wall+1,m) = building_pars_f%pars_xy(ind_tc1,j,i)        
                 ENDIF
                 IF ( building_pars_f%pars_xy(ind_tc2,j,i) /= building_pars_f%fill )    &
                    surf_usm_v(l)%lambda_h(nzb_wall+2,m) = building_pars_f%pars_xy(ind_tc2,j,i)
                 IF ( building_pars_f%pars_xy(ind_tc3,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%lambda_h(nzb_wall+3,m) = building_pars_f%pars_xy(ind_tc3,j,i)    
                 IF ( building_pars_f%pars_xy(ind_tc1,j,i) /= building_pars_f%fill )  THEN
                    surf_usm_v(l)%lambda_h_green(nzb_wall,m)   = building_pars_f%pars_xy(ind_tc1,j,i)         
                    surf_usm_v(l)%lambda_h_green(nzb_wall+1,m) = building_pars_f%pars_xy(ind_tc1,j,i)        
                 ENDIF
                 IF ( building_pars_f%pars_xy(ind_tc2,j,i) /= building_pars_f%fill )    &
                    surf_usm_v(l)%lambda_h_green(nzb_wall+2,m) = building_pars_f%pars_xy(ind_tc2,j,i)
                 IF ( building_pars_f%pars_xy(ind_tc3,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%lambda_h_green(nzb_wall+3,m) = building_pars_f%pars_xy(ind_tc3,j,i)    
                 IF ( building_pars_f%pars_xy(ind_tc1,j,i) /= building_pars_f%fill )  THEN
                    surf_usm_v(l)%lambda_h_window(nzb_wall,m)   = building_pars_f%pars_xy(ind_tc1,j,i)         
                    surf_usm_v(l)%lambda_h_window(nzb_wall+1,m) = building_pars_f%pars_xy(ind_tc1,j,i)        
                 ENDIF
                 IF ( building_pars_f%pars_xy(ind_tc2,j,i) /= building_pars_f%fill )    &
                    surf_usm_v(l)%lambda_h_window(nzb_wall+2,m) = building_pars_f%pars_xy(ind_tc2,j,i)
                 IF ( building_pars_f%pars_xy(ind_tc3,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%lambda_h_window(nzb_wall+3,m) = building_pars_f%pars_xy(ind_tc3,j,i)

                 IF ( building_pars_f%pars_xy(12,j,i) /= building_pars_f%fill )         & 
                    surf_usm_v(l)%target_temp_summer(m)  = building_pars_f%pars_xy(12,j,i)   
                 IF ( building_pars_f%pars_xy(13,j,i) /= building_pars_f%fill )         & 
                    surf_usm_v(l)%target_temp_winter(m)  = building_pars_f%pars_xy(13,j,i)   

                 IF ( building_pars_f%pars_xy(ind_emis_wall,j,i) /= building_pars_f%fill ) & 
                    surf_usm_v(l)%emissivity(ind_veg_wall,m)  = building_pars_f%pars_xy(ind_emis_wall,j,i)
                 IF ( building_pars_f%pars_xy(ind_emis_green,j,i) /= building_pars_f%fill )& 
                    surf_usm_v(l)%emissivity(ind_pav_green,m) = building_pars_f%pars_xy(ind_emis_green,j,i)
                 IF ( building_pars_f%pars_xy(ind_emis_win,j,i) /= building_pars_f%fill )  & 
                    surf_usm_v(l)%emissivity(ind_wat_win,m)   = building_pars_f%pars_xy(ind_emis_win,j,i)

                 IF ( building_pars_f%pars_xy(ind_trans,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%transmissivity(m)      = building_pars_f%pars_xy(ind_trans,j,i)

                 IF ( building_pars_f%pars_xy(ind_z0,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%z0(m)                  = building_pars_f%pars_xy(ind_z0,j,i)
                 IF ( building_pars_f%pars_xy(ind_z0qh,j,i) /= building_pars_f%fill )  & 
                    surf_usm_v(l)%z0h(m)                 = building_pars_f%pars_xy(ind_z0qh,j,i)
                 IF ( building_pars_f%pars_xy(ind_z0qh,j,i) /= building_pars_f%fill )  & 
                    surf_usm_v(l)%z0q(m)                 = building_pars_f%pars_xy(ind_z0qh,j,i)

                 IF ( building_pars_f%pars_xy(ind_alb_wall,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%albedo_type(ind_veg_wall,m)  = building_pars_f%pars_xy(ind_alb_wall,j,i)
                 IF ( building_pars_f%pars_xy(ind_alb_green,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%albedo_type(ind_pav_green,m) = building_pars_f%pars_xy(ind_alb_green,j,i)
                 IF ( building_pars_f%pars_xy(ind_alb_win,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%albedo_type(ind_wat_win,m)   = building_pars_f%pars_xy(ind_alb_win,j,i)

                 IF ( building_pars_f%pars_xy(ind_thick_1,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%zw(nzb_wall,m) = building_pars_f%pars_xy(ind_thick_1,j,i)
                 IF ( building_pars_f%pars_xy(ind_thick_2,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%zw(nzb_wall+1,m) = building_pars_f%pars_xy(ind_thick_2,j,i)
                 IF ( building_pars_f%pars_xy(ind_thick_3,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%zw(nzb_wall+2,m) = building_pars_f%pars_xy(ind_thick_3,j,i)
                 IF ( building_pars_f%pars_xy(ind_thick_4,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%zw(nzb_wall+3,m) = building_pars_f%pars_xy(ind_thick_4,j,i)
                 IF ( building_pars_f%pars_xy(ind_thick_1,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%zw_green(nzb_wall,m) = building_pars_f%pars_xy(ind_thick_1,j,i)
                 IF ( building_pars_f%pars_xy(ind_thick_2,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%zw_green(nzb_wall+1,m) = building_pars_f%pars_xy(ind_thick_2,j,i)
                 IF ( building_pars_f%pars_xy(ind_thick_3,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%zw_green(nzb_wall+2,m) = building_pars_f%pars_xy(ind_thick_3,j,i)
                 IF ( building_pars_f%pars_xy(ind_thick_4,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%zw_green(nzb_wall+3,m) = building_pars_f%pars_xy(ind_thick_4,j,i)
                 IF ( building_pars_f%pars_xy(ind_thick_1,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%zw_window(nzb_wall,m) = building_pars_f%pars_xy(ind_thick_1,j,i)
                 IF ( building_pars_f%pars_xy(ind_thick_2,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%zw_window(nzb_wall+1,m) = building_pars_f%pars_xy(ind_thick_2,j,i)
                 IF ( building_pars_f%pars_xy(ind_thick_3,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%zw_window(nzb_wall+2,m) = building_pars_f%pars_xy(ind_thick_3,j,i)
                 IF ( building_pars_f%pars_xy(ind_thick_4,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%zw_window(nzb_wall+3,m) = building_pars_f%pars_xy(ind_thick_4,j,i)

                 IF ( building_pars_f%pars_xy(45,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%c_surface(m)           = building_pars_f%pars_xy(45,j,i)
                 IF ( building_pars_f%pars_xy(46,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%lambda_surf(m)         = building_pars_f%pars_xy(46,j,i)
                 IF ( building_pars_f%pars_xy(45,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%c_surface_green(m)     = building_pars_f%pars_xy(45,j,i)
                 IF ( building_pars_f%pars_xy(46,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%lambda_surf_green(m)   = building_pars_f%pars_xy(46,j,i)
                 IF ( building_pars_f%pars_xy(45,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%c_surface_window(m)    = building_pars_f%pars_xy(45,j,i)
                 IF ( building_pars_f%pars_xy(46,j,i) /= building_pars_f%fill )    & 
                    surf_usm_v(l)%lambda_surf_window(m)  = building_pars_f%pars_xy(46,j,i)

              ENDDO
           ENDDO
        ENDIF 
!       
!--     Read the surface_types array.
!--     Please note, here also initialization of surface attributes is done as
!--     long as _urbsurf and _surfpar files are available. Values from above
!--     will be overwritten. This might be removed later, but is still in the 
!--     code to enable compatibility with older model version.
        CALL usm_read_urban_surface_types()
        
!--     init material heat model
        CALL usm_init_material_model()
        
!--     init anthropogenic sources of heat
        IF ( usm_anthropogenic_heat )  THEN
!--         init anthropogenic sources of heat (from transportation for now)
            CALL usm_read_anthropogenic_heat()
        ENDIF

        IF ( plant_canopy )  THEN
            
            IF ( .NOT.  ALLOCATED( pc_heating_rate) )  THEN
!--             then pc_heating_rate is allocated in init_plant_canopy
!--             in case of cthf /= 0 => we need to allocate it for our use here
                ALLOCATE( pc_heating_rate(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

                pc_heating_rate = 0.0_wp

            ENDIF

            IF ( .NOT.  ALLOCATED( pc_transpiration_rate) )  THEN
!--             then pc_heating_rate is allocated in init_plant_canopy
!--             in case of cthf /= 0 => we need to allocate it for our use here
                ALLOCATE( pc_transpiration_rate(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

                pc_transpiration_rate = 0.0_wp


            ENDIF
        ENDIF
!
!--    Check for consistent initialization.
!--    Check if roughness length exceed surface-layer height and decrease 
!--    local roughness length where necessary. 
       DO  m = 1, surf_usm_h%ns
          IF ( surf_usm_h%z0(m) >= surf_usm_h%z_mo(m) )  THEN
          
             surf_usm_h%z0(m) = 0.9_wp * surf_usm_h%z_mo(m)
             
             WRITE( message_string, * ) 'z0 exceeds surface-layer height ' //  &
                            'at horizontal urban surface and is ' //           &
                            'decreased appropriately at grid point (i,j) = ',  &
                            surf_usm_h%i(m), surf_usm_h%j(m)
             CALL message( 'urban_surface_model_mod', 'PA0503',                &
                            0, 0, 0, 6, 0 )
          ENDIF
       ENDDO
       
       DO  l = 0, 3
          DO  m = 1, surf_usm_v(l)%ns
             IF ( surf_usm_v(l)%z0(m) >= surf_usm_v(l)%z_mo(m) )  THEN
          
                surf_usm_v(l)%z0(m) = 0.9_wp * surf_usm_v(l)%z_mo(m)
             
                WRITE( message_string, * ) 'z0 exceeds surface-layer height '//&
                            'at vertical urban surface and is ' //             &
                            'decreased appropriately at grid point (i,j) = ',  &
                            surf_usm_v(l)%i(m)+surf_usm_v(l)%ioff,             &
                            surf_usm_v(l)%j(m)+surf_usm_v(l)%joff
                CALL message( 'urban_surface_model_mod', 'PA0503',             &
                            0, 0, 0, 6, 0 )
             ENDIF
          ENDDO
       ENDDO    

!--     Intitialization of the surface and wall/ground/roof temperature

!--     Initialization for restart runs
        IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.        &
             TRIM( initializing_actions ) /= 'cyclic_fill' )  THEN
       
!--         Calculate initial surface temperature from pt of adjacent gridbox
#if ! defined( __nopointer )
            exn(nzb:nzt) = (hyp(nzb:nzt) / 100000.0_wp )**0.286_wp          !< Exner function
#endif

!
!--         At horizontal surfaces. Please note, t_surf_h is defined on a 
!--         different data type, but with the same dimension.
#if ! defined( __nopointer )
            DO  m = 1, surf_usm_h%ns
               i = surf_usm_h%i(m)            
               j = surf_usm_h%j(m)
               k = surf_usm_h%k(m)

               t_surf_h(m) = pt(k,j,i) * exn(k)
               t_surf_window_h(m) = pt(k,j,i) * exn(k)
               t_surf_green_h(m) = pt(k,j,i) * exn(k)
               surf_usm_h%pt_surface(m) = pt(k,j,i) * exn(k)
            ENDDO
!
!--         At vertical surfaces.
            DO  l = 0, 3
               DO  m = 1, surf_usm_v(l)%ns
                  i = surf_usm_v(l)%i(m)            
                  j = surf_usm_v(l)%j(m)
                  k = surf_usm_v(l)%k(m)

                  t_surf_v(l)%t(m) = pt(k,j,i) * exn(k)
                  t_surf_window_v(l)%t(m) = pt(k,j,i) * exn(k)
                  t_surf_green_v(l)%t(m) = pt(k,j,i) * exn(k)
                  surf_usm_v(l)%pt_surface(m) = pt(k,j,i) * exn(k)
               ENDDO
            ENDDO
#endif
      
!--         initial values for t_wall
!--         outer value is set to surface temperature
!--         inner value is set to wall_inner_temperature
!--         and profile is logaritmic (linear in nz).
!--         Horizontal surfaces
            DO  m = 1, surf_usm_h%ns
!
!--            Roof
               IF ( surf_usm_h%isroof_surf(m) )  THEN
                   tin = roof_inner_temperature
                   twin = window_inner_temperature
!
!--            Normal land surface
               ELSE 
                   tin = soil_inner_temperature
                   twin = window_inner_temperature
               ENDIF

               DO k = nzb_wall, nzt_wall+1
                   c = REAL( k - nzb_wall, wp ) /                              &
                       REAL( nzt_wall + 1 - nzb_wall , wp )

                   t_wall_h(k,m) = ( 1.0_wp - c ) * t_surf_h(m) + c * tin
                   t_window_h(k,m) = ( 1.0_wp - c ) * t_surf_window_h(m) + c * twin
                   t_green_h(k,m) = t_surf_h(m)
               ENDDO
            ENDDO
!
!--         Vertical surfaces
            DO  l = 0, 3
               DO  m = 1, surf_usm_v(l)%ns
!
!--               Inner wall
                  tin = wall_inner_temperature
                  twin = window_inner_temperature

                  DO k = nzb_wall, nzt_wall+1
                     c = REAL( k - nzb_wall, wp ) /                            &
                         REAL( nzt_wall + 1 - nzb_wall , wp )
                     t_wall_v(l)%t(k,m) = ( 1.0_wp - c ) * t_surf_v(l)%t(m) + c * tin
                     t_window_v(l)%t(k,m) = ( 1.0_wp - c ) * t_surf_window_v(l)%t(m) + c * twin
                     t_green_v(l)%t(k,m) = t_surf_v(l)%t(m)
                  ENDDO
               ENDDO
            ENDDO
        ELSE
!--         If specified, replace constant wall temperatures with fully 3D values from file
            IF ( read_wall_temp_3d )  CALL usm_read_wall_temperature()
!
        ENDIF
        
!--    
!--     Possibly DO user-defined actions (e.g. define heterogeneous wall surface)
        CALL user_init_urban_surface

!--     initialize prognostic values for the first timestep
        t_surf_h_p = t_surf_h
        t_surf_v_p = t_surf_v
        t_surf_window_h_p = t_surf_window_h
        t_surf_window_v_p = t_surf_window_v
        t_surf_green_h_p = t_surf_green_h
        t_surf_green_v_p = t_surf_green_v
        t_surf_10cm_h_p = t_surf_10cm_h
        t_surf_10cm_v_p = t_surf_10cm_v

        t_wall_h_p = t_wall_h
        t_wall_v_p = t_wall_v
        t_window_h_p = t_window_h
        t_window_v_p = t_window_v
        t_green_h_p = t_green_h
        t_green_v_p = t_green_v

!--     Adjust radiative fluxes for urban surface at model start
        !CALL radiation_interaction
!--     TODO: interaction should be called once before first output,
!--     that is not yet possible.
        
        CALL cpu_log( log_point_s(78), 'usm_init', 'stop' )

    END SUBROUTINE usm_init_urban_surface


!------------------------------------------------------------------------------!
! Description:
! ------------
!
!> Wall model as part of the urban surface model. The model predicts wall
!> temperature.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_material_heat_model


        IMPLICIT NONE

        INTEGER(iwp) ::  i,j,k,l,kw, m                      !< running indices

        REAL(wp), DIMENSION(nzb_wall:nzt_wall) :: wtend, wintend  !< tendency

!
!--     For horizontal surfaces                                   
        DO  m = 1, surf_usm_h%ns
!
!--        Obtain indices
           i = surf_usm_h%i(m)            
           j = surf_usm_h%j(m)
           k = surf_usm_h%k(m)
!
!--        prognostic equation for ground/roof temperature t_wall_h
           wtend(:) = 0.0_wp
           wtend(nzb_wall) = (1.0_wp / surf_usm_h%rho_c_wall(nzb_wall,m)) *        &
                                       ( surf_usm_h%lambda_h(nzb_wall,m) *         &
                                         ( t_wall_h(nzb_wall+1,m)                  &
                                         - t_wall_h(nzb_wall,m) ) *                &
                                         surf_usm_h%ddz_wall(nzb_wall+1,m)         &
                                       + surf_usm_h%frac(ind_veg_wall,m)           &
                                         / (surf_usm_h%frac(ind_veg_wall,m)        &
                                           + surf_usm_h%frac(ind_pav_green,m) )    &
                                         * surf_usm_h%wghf_eb(m)                   &
                                       - surf_usm_h%frac(ind_pav_green,m)          &
                                          / (surf_usm_h%frac(ind_veg_wall,m)       &
                                            + surf_usm_h%frac(ind_pav_green,m) )   &
                                         * ( surf_usm_h%lambda_h_green(nzt_wall,m) &
                                           * surf_usm_h%ddz_green(nzt_wall,m)      &
                                           + surf_usm_h%lambda_h(nzb_wall,m)       &
                                           * surf_usm_h%ddz_wall(nzb_wall,m) )     &
                                         / ( surf_usm_h%ddz_green(nzt_wall,m)      &
                                           + surf_usm_h%ddz_wall(nzb_wall,m) )     &
                                         * ( t_wall_h(nzb_wall,m)                  &
                                           - t_green_h(nzt_wall,m) ) ) *           &
                                       surf_usm_h%ddz_wall_stag(nzb_wall,m)

!dummy value for testing
surf_usm_h%iwghf_eb(m) = 0.

           IF ( indoor_model ) then 
              DO  kw = nzb_wall+1, nzt_wall-1
                  wtend(kw) = (1.0_wp / surf_usm_h%rho_c_wall(kw,m))              &
                                 * (   surf_usm_h%lambda_h(kw,m)                  &
                                    * ( t_wall_h(kw+1,m) - t_wall_h(kw,m) )       &
                                    * surf_usm_h%ddz_wall(kw+1,m)                 &
                                 - surf_usm_h%lambda_h(kw-1,m)                    &
                                    * ( t_wall_h(kw,m) - t_wall_h(kw-1,m) )       &
                                    * surf_usm_h%ddz_wall(kw,m)                   &
                                   ) * surf_usm_h%ddz_wall_stag(kw,m)
              ENDDO
              wtend(nzt_wall) = (1.0_wp / surf_usm_h%rho_c_wall(nzt_wall,m)) *    &
                                         ( surf_usm_h%lambda_h(nzt_wall-1,m) *    &
                                           ( t_wall_h(nzt_wall,m)                 &
                                           - t_wall_h(nzt_wall-1,m) ) *           &
                                           surf_usm_h%ddz_wall(nzt_wall,m)        &
                                         + surf_usm_h%iwghf_eb(m) ) *             &
                                           surf_usm_h%ddz_wall_stag(nzt_wall,m)
           ELSE
              DO  kw = nzb_wall+1, nzt_wall
                  wtend(kw) = (1.0_wp / surf_usm_h%rho_c_wall(kw,m))              &
                                 * (   surf_usm_h%lambda_h(kw,m)                  &
                                    * ( t_wall_h(kw+1,m) - t_wall_h(kw,m) )       &
                                    * surf_usm_h%ddz_wall(kw+1,m)                 &
                                 - surf_usm_h%lambda_h(kw-1,m)                    &
                                    * ( t_wall_h(kw,m) - t_wall_h(kw-1,m) )       &
                                    * surf_usm_h%ddz_wall(kw,m)                   &
                                   ) * surf_usm_h%ddz_wall_stag(kw,m)
              ENDDO
           ENDIF

           t_wall_h_p(nzb_wall:nzt_wall,m) = t_wall_h(nzb_wall:nzt_wall,m)     &
                                 + dt_3d * ( tsc(2)                            &
                                 * wtend(nzb_wall:nzt_wall) + tsc(3)           &
                                 * surf_usm_h%tt_wall_m(nzb_wall:nzt_wall,m) )   

!--        prognostic equation for ground/roof window temperature t_window_h
           wintend(:) = 0.0_wp
           wintend(nzb_wall) = (1.0_wp / surf_usm_h%rho_c_window(nzb_wall,m)) *   &
                                      ( surf_usm_h%lambda_h_window(nzb_wall,m) *  &
                                        ( t_window_h(nzb_wall+1,m)                &
                                        - t_window_h(nzb_wall,m) ) *              &
                                        surf_usm_h%ddz_window(nzb_wall+1,m)       &
                                      + surf_usm_h%wghf_eb_window(m)              &
                                      + surf_usm_h%rad_sw_in(m)                   &
                                        * (1.0_wp - exp(-surf_usm_h%transmissivity(m) &
                                        * surf_usm_h%zw_window(nzb_wall,m) ) )    &
                                      ) * surf_usm_h%ddz_window_stag(nzb_wall,m)

           IF ( indoor_model ) then 
              DO  kw = nzb_wall+1, nzt_wall-1
                  wintend(kw) = (1.0_wp / surf_usm_h%rho_c_window(kw,m))          &
                                 * (   surf_usm_h%lambda_h_window(kw,m)           &
                                    * ( t_window_h(kw+1,m) - t_window_h(kw,m) )   &
                                    * surf_usm_h%ddz_window(kw+1,m)               &
                                 - surf_usm_h%lambda_h_window(kw-1,m)             &
                                    * ( t_window_h(kw,m) - t_window_h(kw-1,m) )   &
                                    * surf_usm_h%ddz_window(kw,m)                 &
                                 + surf_usm_h%rad_sw_in(m)                        &
                                    * (exp(-surf_usm_h%transmissivity(m)       &
                                        * surf_usm_h%zw_window(kw-1,m) )          &
                                        - exp(-surf_usm_h%transmissivity(m)    &
                                        * surf_usm_h%zw_window(kw,m) ) )          &
                                   ) * surf_usm_h%ddz_window_stag(kw,m)

              ENDDO
              wintend(nzt_wall) = (1.0_wp / surf_usm_h%rho_c_window(nzt_wall,m)) *     &
                                         ( surf_usm_h%lambda_h_window(nzt_wall-1,m) *  &
                                           ( t_window_h(nzt_wall,m)                    &
                                           - t_window_h(nzt_wall-1,m) ) *              &
                                           surf_usm_h%ddz_window(nzt_wall,m)           &
                                         + surf_usm_h%iwghf_eb_window(m)               &
                                         + surf_usm_h%rad_sw_in(m)                     &
                                           * (1.0_wp - exp(-surf_usm_h%transmissivity(m) &
                                           * surf_usm_h%zw_window(nzt_wall,m) ) )      &
                                         ) * surf_usm_h%ddz_window_stag(nzt_wall,m)
           ELSE
              DO  kw = nzb_wall+1, nzt_wall
                  wintend(kw) = (1.0_wp / surf_usm_h%rho_c_window(kw,m))          &
                                 * (   surf_usm_h%lambda_h_window(kw,m)           &
                                    * ( t_window_h(kw+1,m) - t_window_h(kw,m) )   &
                                    * surf_usm_h%ddz_window(kw+1,m)               &
                                 - surf_usm_h%lambda_h_window(kw-1,m)             &
                                    * ( t_window_h(kw,m) - t_window_h(kw-1,m) )   &
                                    * surf_usm_h%ddz_window(kw,m)                 &
                                 + surf_usm_h%rad_sw_in(m)                        &
                                    * (exp(-surf_usm_h%transmissivity(m)       &
                                        * surf_usm_h%zw_window(kw-1,m) )          &
                                        - exp(-surf_usm_h%transmissivity(m)    &
                                        * surf_usm_h%zw_window(kw,m) ) )          &
                                   ) * surf_usm_h%ddz_window_stag(kw,m)

              ENDDO
           ENDIF

           t_window_h_p(nzb_wall:nzt_wall,m) = t_window_h(nzb_wall:nzt_wall,m)    &
                                 + dt_3d * ( tsc(2)                               &
                                 * wintend(nzb_wall:nzt_wall) + tsc(3)            &
                                 * surf_usm_h%tt_window_m(nzb_wall:nzt_wall,m) )   
           
!
!--        calculate t_wall tendencies for the next Runge-Kutta step
           IF ( timestep_scheme(1:5) == 'runge' )  THEN
               IF ( intermediate_timestep_count == 1 )  THEN
                  DO  kw = nzb_wall, nzt_wall
                     surf_usm_h%tt_wall_m(kw,m) = wtend(kw)
                  ENDDO
               ELSEIF ( intermediate_timestep_count <                          &
                        intermediate_timestep_count_max )  THEN
                   DO  kw = nzb_wall, nzt_wall
                      surf_usm_h%tt_wall_m(kw,m) = -9.5625_wp * wtend(kw) +    &
                                         5.3125_wp * surf_usm_h%tt_wall_m(kw,m)
                   ENDDO
               ENDIF
           ENDIF

!--        calculate t_window tendencies for the next Runge-Kutta step
           IF ( timestep_scheme(1:5) == 'runge' )  THEN
               IF ( intermediate_timestep_count == 1 )  THEN
                  DO  kw = nzb_wall, nzt_wall
                     surf_usm_h%tt_window_m(kw,m) = wintend(kw)
                  ENDDO
               ELSEIF ( intermediate_timestep_count <                          &
                        intermediate_timestep_count_max )  THEN
                   DO  kw = nzb_wall, nzt_wall
                      surf_usm_h%tt_window_m(kw,m) = -9.5625_wp * wintend(kw) +    &
                                         5.3125_wp * surf_usm_h%tt_window_m(kw,m)
                   ENDDO
               ENDIF
           ENDIF
        ENDDO

!
!--     For vertical surfaces     
        DO  l = 0, 3                              
           DO  m = 1, surf_usm_v(l)%ns
!
!--           Obtain indices
              i = surf_usm_v(l)%i(m)            
              j = surf_usm_v(l)%j(m)
              k = surf_usm_v(l)%k(m)
!
!--           prognostic equation for wall temperature t_wall_v
              wtend(:) = 0.0_wp

               wtend(nzb_wall) = (1.0_wp / surf_usm_v(l)%rho_c_wall(nzb_wall,m)) *    &
                                       ( surf_usm_v(l)%lambda_h(nzb_wall,m) *         &
                                         ( t_wall_v(l)%t(nzb_wall+1,m)                &
                                         - t_wall_v(l)%t(nzb_wall,m) ) *              &
                                         surf_usm_v(l)%ddz_wall(nzb_wall+1,m)         &
                                       + surf_usm_v(l)%frac(ind_veg_wall,m)           &
                                         / (surf_usm_v(l)%frac(ind_veg_wall,m)        &
                                           + surf_usm_v(l)%frac(ind_pav_green,m) )    &
                                         * surf_usm_v(l)%wghf_eb(m)                   &
                                       - surf_usm_v(l)%frac(ind_pav_green,m)          &
                                         / (surf_usm_v(l)%frac(ind_veg_wall,m)        &
                                           + surf_usm_v(l)%frac(ind_pav_green,m) )    &
                                         * ( surf_usm_v(l)%lambda_h_green(nzt_wall,m) &
                                           * surf_usm_v(l)%ddz_green(nzt_wall,m)      &
                                           + surf_usm_v(l)%lambda_h(nzb_wall,m)       &
                                           * surf_usm_v(l)%ddz_wall(nzb_wall,m) )     &
                                         / ( surf_usm_v(l)%ddz_green(nzt_wall,m)      &
                                           + surf_usm_v(l)%ddz_wall(nzb_wall,m) )     &
                                         * ( t_wall_v(l)%t(nzb_wall,m)                &
                                           - t_green_v(l)%t(nzt_wall,m) ) ) *         &
                                         surf_usm_v(l)%ddz_wall_stag(nzb_wall,m)

!dummy value for testing
surf_usm_v(l)%iwghf_eb(m) = 0.

              IF ( indoor_model ) then 
                 DO  kw = nzb_wall+1, nzt_wall-1
                     wtend(kw) = (1.0_wp / surf_usm_v(l)%rho_c_wall(kw,m))        &
                              * (   surf_usm_v(l)%lambda_h(kw,m)                  &
                                 * ( t_wall_v(l)%t(kw+1,m) - t_wall_v(l)%t(kw,m) )&
                                 * surf_usm_v(l)%ddz_wall(kw+1,m)                 &
                              - surf_usm_v(l)%lambda_h(kw-1,m)                    &
                                 * ( t_wall_v(l)%t(kw,m) - t_wall_v(l)%t(kw-1,m) )&
                                 * surf_usm_v(l)%ddz_wall(kw,m)                   &
                                 ) * surf_usm_v(l)%ddz_wall_stag(kw,m)
                 ENDDO
                 wtend(nzt_wall) = (1.0_wp / surf_usm_v(l)%rho_c_wall(nzt_wall,m)) * &
                                         ( surf_usm_v(l)%lambda_h(nzt_wall-1,m) *    &
                                           ( t_wall_v(l)%t(nzt_wall,m)               &
                                           - t_wall_v(l)%t(nzt_wall-1,m) ) *         &
                                           surf_usm_v(l)%ddz_wall(nzt_wall,m)        &
                                         + surf_usm_v(l)%iwghf_eb(m) ) *             &
                                           surf_usm_v(l)%ddz_wall_stag(nzt_wall,m)
              ELSE
                 DO  kw = nzb_wall+1, nzt_wall
                     wtend(kw) = (1.0_wp / surf_usm_v(l)%rho_c_wall(kw,m))        &
                              * (   surf_usm_v(l)%lambda_h(kw,m)                  &
                                 * ( t_wall_v(l)%t(kw+1,m) - t_wall_v(l)%t(kw,m) )&
                                 * surf_usm_v(l)%ddz_wall(kw+1,m)                 &
                              - surf_usm_v(l)%lambda_h(kw-1,m)                    &
                                 * ( t_wall_v(l)%t(kw,m) - t_wall_v(l)%t(kw-1,m) )&
                                 * surf_usm_v(l)%ddz_wall(kw,m)                   &
                                 ) * surf_usm_v(l)%ddz_wall_stag(kw,m)
                 ENDDO
              ENDIF

              t_wall_v_p(l)%t(nzb_wall:nzt_wall,m) =                           &
                                   t_wall_v(l)%t(nzb_wall:nzt_wall,m)          &
                                 + dt_3d * ( tsc(2)                            &
                                 * wtend(nzb_wall:nzt_wall) + tsc(3)           &
                                 * surf_usm_v(l)%tt_wall_m(nzb_wall:nzt_wall,m) )   

!--           prognostic equation for window temperature t_window_v
              wintend(:) = 0.0_wp
              wintend(nzb_wall) = (1.0_wp / surf_usm_v(l)%rho_c_window(nzb_wall,m)) * &
                                      ( surf_usm_v(l)%lambda_h_window(nzb_wall,m) *   &
                                        ( t_window_v(l)%t(nzb_wall+1,m)               &
                                        - t_window_v(l)%t(nzb_wall,m) ) *             &
                                        surf_usm_v(l)%ddz_window(nzb_wall+1,m)        &
                                      + surf_usm_v(l)%wghf_eb_window(m)               &
                                      + surf_usm_v(l)%rad_sw_in(m)                    &
                                        * (1.0_wp - exp(-surf_usm_v(l)%transmissivity(m) &
                                        * surf_usm_v(l)%zw_window(nzb_wall,m) ) )     &
                                      ) * surf_usm_v(l)%ddz_window_stag(nzb_wall,m)

              IF ( indoor_model ) then 
                 DO  kw = nzb_wall+1, nzt_wall -1
                     wintend(kw) = (1.0_wp / surf_usm_v(l)%rho_c_window(kw,m))         &
                              * (   surf_usm_v(l)%lambda_h_window(kw,m)                &
                                 * ( t_window_v(l)%t(kw+1,m) - t_window_v(l)%t(kw,m) ) &
                                 * surf_usm_v(l)%ddz_window(kw+1,m)                    &
                              - surf_usm_v(l)%lambda_h_window(kw-1,m)                  &
                                 * ( t_window_v(l)%t(kw,m) - t_window_v(l)%t(kw-1,m) ) &
                                 * surf_usm_v(l)%ddz_window(kw,m)                      &
                              + surf_usm_v(l)%rad_sw_in(m)                             &
                                 * (exp(-surf_usm_v(l)%transmissivity(m)            &
                                    * surf_usm_v(l)%zw_window(kw-1,m)       )          &
                                        - exp(-surf_usm_v(l)%transmissivity(m)      &
                                        * surf_usm_v(l)%zw_window(kw,m) ) )            &
                                 ) * surf_usm_v(l)%ddz_window_stag(kw,m)
                  ENDDO
                  wintend(nzt_wall) = (1.0_wp / surf_usm_v(l)%rho_c_window(nzt_wall,m)) * &
                                          ( surf_usm_v(l)%lambda_h_window(nzt_wall-1,m) * &
                                            ( t_window_v(l)%t(nzt_wall,m)                 &
                                            - t_window_v(l)%t(nzt_wall-1,m) ) *           &
                                            surf_usm_v(l)%ddz_window(nzt_wall,m)          &
                                          + surf_usm_v(l)%iwghf_eb_window(m)              &
                                          + surf_usm_v(l)%rad_sw_in(m)                    &
                                            * (1.0_wp - exp(-surf_usm_v(l)%transmissivity(m) &
                                            * surf_usm_v(l)%zw_window(nzt_wall,m) ) )     &
                                          ) * surf_usm_v(l)%ddz_window_stag(nzt_wall,m)
              ELSE
                 DO  kw = nzb_wall+1, nzt_wall
                     wintend(kw) = (1.0_wp / surf_usm_v(l)%rho_c_window(kw,m))         &
                              * (   surf_usm_v(l)%lambda_h_window(kw,m)                &
                                 * ( t_window_v(l)%t(kw+1,m) - t_window_v(l)%t(kw,m) ) &
                                 * surf_usm_v(l)%ddz_window(kw+1,m)                    &
                              - surf_usm_v(l)%lambda_h_window(kw-1,m)                  &
                                 * ( t_window_v(l)%t(kw,m) - t_window_v(l)%t(kw-1,m) ) &
                                 * surf_usm_v(l)%ddz_window(kw,m)                      &
                              + surf_usm_v(l)%rad_sw_in(m)                             &
                                 * (exp(-surf_usm_v(l)%transmissivity(m)            &
                                    * surf_usm_v(l)%zw_window(kw-1,m)       )          &
                                        - exp(-surf_usm_v(l)%transmissivity(m)      &
                                        * surf_usm_v(l)%zw_window(kw,m) ) )            &
                                 ) * surf_usm_v(l)%ddz_window_stag(kw,m)
                 ENDDO
              ENDIF

              t_window_v_p(l)%t(nzb_wall:nzt_wall,m) =                           &
                                   t_window_v(l)%t(nzb_wall:nzt_wall,m)          &
                                 + dt_3d * ( tsc(2)                              &
                                 * wintend(nzb_wall:nzt_wall) + tsc(3)           &
                                 * surf_usm_v(l)%tt_window_m(nzb_wall:nzt_wall,m) )   

!
!--           calculate t_wall tendencies for the next Runge-Kutta step
              IF ( timestep_scheme(1:5) == 'runge' )  THEN
                  IF ( intermediate_timestep_count == 1 )  THEN
                     DO  kw = nzb_wall, nzt_wall
                        surf_usm_v(l)%tt_wall_m(kw,m) = wtend(kw)
                     ENDDO
                  ELSEIF ( intermediate_timestep_count <                       &
                           intermediate_timestep_count_max )  THEN
                      DO  kw = nzb_wall, nzt_wall
                         surf_usm_v(l)%tt_wall_m(kw,m) =                       &
                                     - 9.5625_wp * wtend(kw) +                 &
                                       5.3125_wp * surf_usm_v(l)%tt_wall_m(kw,m)
                      ENDDO
                  ENDIF
              ENDIF
!--           calculate t_window tendencies for the next Runge-Kutta step
              IF ( timestep_scheme(1:5) == 'runge' )  THEN
                  IF ( intermediate_timestep_count == 1 )  THEN
                     DO  kw = nzb_wall, nzt_wall
                        surf_usm_v(l)%tt_window_m(kw,m) = wintend(kw)
                     ENDDO
                  ELSEIF ( intermediate_timestep_count <                       &
                           intermediate_timestep_count_max )  THEN
                      DO  kw = nzb_wall, nzt_wall
                         surf_usm_v(l)%tt_window_m(kw,m) =                     &
                                     - 9.5625_wp * wintend(kw) +               &
                                       5.3125_wp * surf_usm_v(l)%tt_window_m(kw,m)
                      ENDDO
                  ENDIF
              ENDIF
           ENDDO
        ENDDO

    END SUBROUTINE usm_material_heat_model

!------------------------------------------------------------------------------!
! Description:
! ------------
!
!> Green and substrate model as part of the urban surface model. The model predicts ground
!> temperatures.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_green_heat_model


        IMPLICIT NONE

        INTEGER(iwp) ::  i,j,k,l,kw, m                      !< running indices

        REAL(wp), DIMENSION(nzb_wall:nzt_wall) :: gtend  !< tendency

!
!--     For horizontal surfaces                                   
        DO  m = 1, surf_usm_h%ns
!
!--        Obtain indices
           i = surf_usm_h%i(m)            
           j = surf_usm_h%j(m)
           k = surf_usm_h%k(m)

           t_green_h(nzt_wall+1,m) = t_wall_h(nzb_wall,m)
!
!--        prognostic equation for ground/roof temperature t_green_h
           gtend(:) = 0.0_wp
           gtend(nzb_wall) = (1.0_wp / surf_usm_h%rho_c_green(nzb_wall,m)) *    &
                                      ( surf_usm_h%lambda_h_green(nzb_wall,m) * &
                                        ( t_green_h(nzb_wall+1,m)               &
                                        - t_green_h(nzb_wall,m) ) *             &
                                        surf_usm_h%ddz_green(nzb_wall+1,m)      &
                                      + surf_usm_h%wghf_eb_green(m) ) *         &
                                        surf_usm_h%ddz_green_stag(nzb_wall,m)
           
            DO  kw = nzb_wall+1, nzt_wall
                gtend(kw) = (1.0_wp / surf_usm_h%rho_c_green(kw,m))             &
                               * (   surf_usm_h%lambda_h_green(kw,m)            &
                                  * ( t_green_h(kw+1,m) - t_green_h(kw,m) )     &
                                  * surf_usm_h%ddz_green(kw+1,m)                &
                               - surf_usm_h%lambda_h_green(kw-1,m)              &
                                  * ( t_green_h(kw,m) - t_green_h(kw-1,m) )     &
                                  * surf_usm_h%ddz_green(kw,m)                  &
                                 ) * surf_usm_h%ddz_green_stag(kw,m)
            ENDDO

           t_green_h_p(nzb_wall:nzt_wall,m) = t_green_h(nzb_wall:nzt_wall,m)    &
                                 + dt_3d * ( tsc(2)                             &
                                 * gtend(nzb_wall:nzt_wall) + tsc(3)            &
                                 * surf_usm_h%tt_green_m(nzb_wall:nzt_wall,m) )   

          
!
!--        calculate t_green tendencies for the next Runge-Kutta step
           IF ( timestep_scheme(1:5) == 'runge' )  THEN
               IF ( intermediate_timestep_count == 1 )  THEN
                  DO  kw = nzb_wall, nzt_wall
                     surf_usm_h%tt_green_m(kw,m) = gtend(kw)
                  ENDDO
               ELSEIF ( intermediate_timestep_count <                           &
                        intermediate_timestep_count_max )  THEN
                   DO  kw = nzb_wall, nzt_wall
                      surf_usm_h%tt_green_m(kw,m) = -9.5625_wp * gtend(kw) +    &
                                         5.3125_wp * surf_usm_h%tt_green_m(kw,m)
                   ENDDO
               ENDIF
           ENDIF
        ENDDO

!
!--     For vertical surfaces     
        DO  l = 0, 3                              
           DO  m = 1, surf_usm_v(l)%ns
!
!--           Obtain indices
              i = surf_usm_v(l)%i(m)            
              j = surf_usm_v(l)%j(m)
              k = surf_usm_v(l)%k(m)

              t_green_v(l)%t(nzt_wall+1,m) = t_wall_v(l)%t(nzb_wall,m)
!
!--           prognostic equation for green temperature t_green_v
              gtend(:) = 0.0_wp
              gtend(nzb_wall) = (1.0_wp / surf_usm_v(l)%rho_c_green(nzb_wall,m)) * &
                                      ( surf_usm_v(l)%lambda_h_green(nzb_wall,m) * &
                                        ( t_green_v(l)%t(nzb_wall+1,m)             &
                                        - t_green_v(l)%t(nzb_wall,m) ) *           &
                                        surf_usm_v(l)%ddz_green(nzb_wall+1,m)      &
                                      + surf_usm_v(l)%wghf_eb(m) ) *               &
                                        surf_usm_v(l)%ddz_green_stag(nzb_wall,m)
            
              DO  kw = nzb_wall+1, nzt_wall
                 gtend(kw) = (1.0_wp / surf_usm_v(l)%rho_c_green(kw,m))          &
                           * (   surf_usm_v(l)%lambda_h_green(kw,m)              &
                             * ( t_green_v(l)%t(kw+1,m) - t_green_v(l)%t(kw,m) ) &
                             * surf_usm_v(l)%ddz_green(kw+1,m)                   &
                           - surf_usm_v(l)%lambda_h(kw-1,m)                      &
                             * ( t_green_v(l)%t(kw,m) - t_green_v(l)%t(kw-1,m) ) &
                             * surf_usm_v(l)%ddz_green(kw,m) )                   &
                           * surf_usm_v(l)%ddz_green_stag(kw,m)
              ENDDO

              t_green_v_p(l)%t(nzb_wall:nzt_wall,m) =                              &
                                   t_green_v(l)%t(nzb_wall:nzt_wall,m)             &
                                 + dt_3d * ( tsc(2)                                &
                                 * gtend(nzb_wall:nzt_wall) + tsc(3)               &
                                 * surf_usm_v(l)%tt_green_m(nzb_wall:nzt_wall,m) )   

!
!--           calculate t_green tendencies for the next Runge-Kutta step
              IF ( timestep_scheme(1:5) == 'runge' )  THEN
                  IF ( intermediate_timestep_count == 1 )  THEN
                     DO  kw = nzb_wall, nzt_wall
                        surf_usm_v(l)%tt_green_m(kw,m) = gtend(kw)
                     ENDDO
                  ELSEIF ( intermediate_timestep_count <                           &
                           intermediate_timestep_count_max )  THEN
                      DO  kw = nzb_wall, nzt_wall
                         surf_usm_v(l)%tt_green_m(kw,m) =                          &
                                     - 9.5625_wp * gtend(kw) +                     &
                                       5.3125_wp * surf_usm_v(l)%tt_green_m(kw,m)
                      ENDDO
                  ENDIF
              ENDIF

           ENDDO
        ENDDO

    END SUBROUTINE usm_green_heat_model

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &usm_par for urban surface model
!------------------------------------------------------------------------------!
    SUBROUTINE usm_parin

       IMPLICIT NONE

       CHARACTER (LEN=80) ::  line  !< string containing current line of file PARIN

       NAMELIST /urban_surface_par/                                            &
                           building_type,                                      &
                           land_category,                                      &
                           naheatlayers,                                       &
                           pedestrian_category,                                &
                           roughness_concrete,                                 &
                           read_wall_temp_3d,                                  &
                           roof_category,                                      &
                           urban_surface,                                      &
                           usm_anthropogenic_heat,                             &
                           usm_material_model,                                 &
                           wall_category,                                      &
                           indoor_model,                                       &
                           wall_inner_temperature,                             &
                           roof_inner_temperature,                             &
                           soil_inner_temperature,                             &
                           window_inner_temperature

       NAMELIST /urban_surface_parameters/                                     &
                           building_type,                                      &
                           land_category,                                      &
                           naheatlayers,                                       &
                           pedestrian_category,                                &
                           roughness_concrete,                                 &
                           read_wall_temp_3d,                                  &
                           roof_category,                                      &
                           urban_surface,                                      &
                           usm_anthropogenic_heat,                             &
                           usm_material_model,                                 &
                           wall_category,                                      &
                           indoor_model,                                       &
                           wall_inner_temperature,                             &
                           roof_inner_temperature,                             &
                           soil_inner_temperature,                             &
                           window_inner_temperature
!
!--    Try to find urban surface model package
       REWIND ( 11 )
       line = ' '
       DO   WHILE ( INDEX( line, '&urban_surface_parameters' ) == 0 )
          READ ( 11, '(A)', END=10 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read user-defined namelist
       READ ( 11, urban_surface_parameters )
!
!--    Set flag that indicates that the land surface model is switched on
       urban_surface = .TRUE.

       GOTO 12
!
!--    Try to find old namelist
 10    REWIND ( 11 )
       line = ' '
       DO   WHILE ( INDEX( line, '&urban_surface_par' ) == 0 )
          READ ( 11, '(A)', END=12 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read user-defined namelist
       READ ( 11, urban_surface_par )

       message_string = 'namelist urban_surface_par is deprecated and will be ' // &
                     'removed in near future. Please use namelist ' //             &
                     'urban_surface_parameters instead' 
       CALL message( 'usm_parin', 'PA0487', 0, 1, 0, 6, 0 )
!
!--    Set flag that indicates that the land surface model is switched on
       urban_surface = .TRUE.

 12    CONTINUE


    END SUBROUTINE usm_parin

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates temperature near surface (10 cm) for indoor model
!------------------------------------------------------------------------------!
    SUBROUTINE usm_temperature_near_surface

       IMPLICIT NONE

       INTEGER(iwp)                          :: i, j, k, l, m   !< running indices

!        
!--    First, treat horizontal surface elements
       DO  m = 1, surf_usm_h%ns

!--       Get indices of respective grid point
          i = surf_usm_h%i(m)
          j = surf_usm_h%j(m)
          k = surf_usm_h%k(m)

          t_surf_10cm_h(m) = surf_usm_h%pt_surface(m) + surf_usm_h%ts(m) / kappa        &
                             * ( log( 0.1_wp /  surf_usm_h%z0h(m) )              &
                               - psi_h( 0.1_wp / surf_usm_h%ol(m) )              &
                               + psi_h( surf_usm_h%z0h(m) / surf_usm_h%ol(m) ) )

       ENDDO
!
!--    Now, treat vertical surface elements
       DO  l = 0, 3
          DO  m = 1, surf_usm_v(l)%ns

!--          Get indices of respective grid point
             i = surf_usm_v(l)%i(m)
             j = surf_usm_v(l)%j(m)
             k = surf_usm_v(l)%k(m)

             t_surf_10cm_v(l)%t(m) =surf_usm_v(l)%pt_surface(m) + surf_usm_v(l)%ts(m) / kappa &
                                     * ( log( 0.1_wp / surf_usm_v(l)%z0h(m) )             &
                                       - psi_h( 0.1_wp / surf_usm_v(l)%ol(m) )            &
                                       + psi_h( surf_usm_v(l)%z0h(m) / surf_usm_v(l)%ol(m) ) )

          ENDDO

       ENDDO


    END SUBROUTINE usm_temperature_near_surface

    
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!
!> This subroutine is part of the urban surface model.
!> It reads daily heat produced by anthropogenic sources
!> and the diurnal cycle of the heat.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_read_anthropogenic_heat
    
        INTEGER(iwp)                  :: i,j,k,ii
        REAL(wp)                      :: heat

!--     allocation of array of sources of anthropogenic heat and their diural profile
        ALLOCATE( aheat(naheatlayers,nys:nyn,nxl:nxr) )
        ALLOCATE( aheatprof(naheatlayers,0:24) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--     read daily amount of heat and its daily cycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        aheat = 0.0_wp
        DO  ii = 0, io_blocks-1
            IF ( ii == io_group )  THEN

!--             open anthropogenic heat file 
                OPEN( 151, file='ANTHROPOGENIC_HEAT'//TRIM(coupling_char), action='read', &
                           status='old', form='formatted', err=11 )
                i = 0
                j = 0
                DO
                    READ( 151, *, err=12, end=13 )  i, j, k, heat
                    IF ( i >= nxl  .AND.  i <= nxr  .AND.  j >= nys  .AND.  j <= nyn )  THEN
                        IF ( k <= naheatlayers  .AND.  k > get_topography_top_index_ji( j, i, 's' ) )  THEN
!--                         write heat into the array
                            aheat(k,j,i) = heat
                        ENDIF
                    ENDIF
                    CYCLE
 12                 WRITE(message_string,'(a,2i4)') 'error in file ANTHROPOGENIC_HEAT'//TRIM(coupling_char)//' after line ',i,j
                    CALL message( 'usm_read_anthropogenic_heat', 'PA0515', 0, 1, 0, 6, 0 )
                ENDDO
 13             CLOSE(151)
                CYCLE
 11             message_string = 'file ANTHROPOGENIC_HEAT'//TRIM(coupling_char)//' does not exist'
                CALL message( 'usm_read_anthropogenic_heat', 'PA0516', 1, 2, 0, 6, 0 )
            ENDIF
            
#if defined( __parallel ) && ! defined ( __check )
            CALL MPI_BARRIER( comm2d, ierr )
#endif
        ENDDO
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--     read diurnal profiles of heat sources
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        aheatprof = 0.0_wp
        DO  ii = 0, io_blocks-1
            IF ( ii == io_group )  THEN

!--             open anthropogenic heat profile file 
                OPEN( 151, file='ANTHROPOGENIC_HEAT_PROFILE'//TRIM(coupling_char), action='read', &
                           status='old', form='formatted', err=21 )
                i = 0
                DO
                    READ( 151, *, err=22, end=23 )  i, k, heat
                    IF ( i >= 0  .AND.  i <= 24  .AND.  k <= naheatlayers )  THEN
!--                     write heat into the array
                        aheatprof(k,i) = heat
                    ENDIF
                    CYCLE
 22                 WRITE(message_string,'(a,i4)') 'error in file ANTHROPOGENIC_HEAT_PROFILE'// &
                                                     TRIM(coupling_char)//' after line ',i
                    CALL message( 'usm_read_anthropogenic_heat', 'PA0517', 0, 1, 0, 6, 0 )
                ENDDO
                aheatprof(:,24) = aheatprof(:,0)
 23             CLOSE(151)
                CYCLE
 21             message_string = 'file ANTHROPOGENIC_HEAT_PROFILE'//TRIM(coupling_char)//' does not exist'
                CALL message( 'usm_read_anthropogenic_heat', 'PA0518', 1, 2, 0, 6, 0 )
            ENDIF
            
#if defined( __parallel ) && ! defined ( __check )
            CALL MPI_BARRIER( comm2d, ierr )
#endif
        ENDDO
        
    END SUBROUTINE usm_read_anthropogenic_heat
   

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Soubroutine reads t_surf and t_wall data from restart files
!------------------------------------------------------------------------------!
    SUBROUTINE usm_rrd_local( i, k, nxlf, nxlc, nxl_on_file, nxrf, nxrc,       &
                              nxr_on_file, nynf, nync, nyn_on_file, nysf, nysc,&
                              nys_on_file, found )


       USE control_parameters,                                                 &
           ONLY: length, restart_string
           
       IMPLICIT NONE

       CHARACTER (LEN=1)  ::  dum              !< dummy to create correct string for reading input variable

       INTEGER(iwp)       ::  l                !< index variable for surface type
       INTEGER(iwp)       ::  i                !< running index over input files
       INTEGER(iwp)       ::  k                !< running index over previous input files covering current local domain
       INTEGER(iwp)       ::  ns_h_on_file_usm !< number of horizontal surface elements (urban type) on file
       INTEGER(iwp)       ::  nxlc             !< index of left boundary on current subdomain
       INTEGER(iwp)       ::  nxlf             !< index of left boundary on former subdomain 
       INTEGER(iwp)       ::  nxl_on_file      !< index of left boundary on former local domain 
       INTEGER(iwp)       ::  nxrc             !< index of right boundary on current subdomain
       INTEGER(iwp)       ::  nxrf             !< index of right boundary on former subdomain
       INTEGER(iwp)       ::  nxr_on_file      !< index of right boundary on former local domain 
       INTEGER(iwp)       ::  nync             !< index of north boundary on current subdomain
       INTEGER(iwp)       ::  nynf             !< index of north boundary on former subdomain
       INTEGER(iwp)       ::  nyn_on_file      !< index of north boundary on former local domain 
       INTEGER(iwp)       ::  nysc             !< index of south boundary on current subdomain 
       INTEGER(iwp)       ::  nysf             !< index of south boundary on former subdomain
       INTEGER(iwp)       ::  nys_on_file      !< index of south boundary on former local domain
       
       INTEGER(iwp)       ::  ns_v_on_file_usm(0:3) !< number of vertical surface elements (urban type) on file
       
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  start_index_on_file 
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  end_index_on_file

       LOGICAL, INTENT(OUT)  ::  found 
       
       REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE   ::  tmp_surf_h, tmp_surf_window_h, tmp_surf_green_h
       REAL(wp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  tmp_wall_h, tmp_window_h, tmp_green_h
       
       TYPE( t_surf_vertical ), DIMENSION(0:3), SAVE ::  tmp_surf_v, tmp_surf_window_v, tmp_surf_green_v
       TYPE( t_wall_vertical ), DIMENSION(0:3), SAVE ::  tmp_wall_v, tmp_window_v, tmp_green_v


       found = .TRUE.


          SELECT CASE ( restart_string(1:length) ) 

             CASE ( 'ns_h_on_file_usm') 
                IF ( k == 1 )  THEN
                   READ ( 13 ) ns_h_on_file_usm
                
                   IF ( ALLOCATED( tmp_surf_h ) ) DEALLOCATE( tmp_surf_h )
                   IF ( ALLOCATED( tmp_wall_h ) ) DEALLOCATE( tmp_wall_h ) 
                   IF ( ALLOCATED( tmp_surf_window_h ) )                       &
                      DEALLOCATE( tmp_surf_window_h ) 
                   IF ( ALLOCATED( tmp_window_h) ) DEALLOCATE( tmp_window_h ) 
                   IF ( ALLOCATED( tmp_surf_green_h) )                         &
                      DEALLOCATE( tmp_surf_green_h ) 
                   IF ( ALLOCATED( tmp_green_h) ) DEALLOCATE( tmp_green_h )
  
!
!--                Allocate temporary arrays for reading data on file. Note,
!--                the size of allocated surface elements do not necessarily
!--                need  to match the size of present surface elements on 
!--                current processor, as the number of processors between 
!--                restarts can change. 
                   ALLOCATE( tmp_surf_h(1:ns_h_on_file_usm) )
                   ALLOCATE( tmp_wall_h(nzb_wall:nzt_wall+1,                   &
                                        1:ns_h_on_file_usm) )
                   ALLOCATE( tmp_surf_window_h(1:ns_h_on_file_usm) )
                   ALLOCATE( tmp_window_h(nzb_wall:nzt_wall+1,                 &
                                          1:ns_h_on_file_usm) )
                   ALLOCATE( tmp_surf_green_h(1:ns_h_on_file_usm) )
                   ALLOCATE( tmp_green_h(nzb_wall:nzt_wall+1,                  &
                                         1:ns_h_on_file_usm) )

                ENDIF

             CASE ( 'ns_v_on_file_usm')
                IF ( k == 1 )  THEN
                   READ ( 13 ) ns_v_on_file_usm 

                   DO  l = 0, 3
                      IF ( ALLOCATED( tmp_surf_v(l)%t ) )                      &
                         DEALLOCATE( tmp_surf_v(l)%t )
                      IF ( ALLOCATED( tmp_wall_v(l)%t ) )                      &
                         DEALLOCATE( tmp_wall_v(l)%t )
                      IF ( ALLOCATED( tmp_surf_window_v(l)%t ) )               & 
                         DEALLOCATE( tmp_surf_window_v(l)%t )
                      IF ( ALLOCATED( tmp_window_v(l)%t ) )                    &
                         DEALLOCATE( tmp_window_v(l)%t )
                      IF ( ALLOCATED( tmp_surf_green_v(l)%t ) )                &
                         DEALLOCATE( tmp_surf_green_v(l)%t )
                      IF ( ALLOCATED( tmp_green_v(l)%t ) )                     &
                         DEALLOCATE( tmp_green_v(l)%t )
                   ENDDO 

!
!--                Allocate temporary arrays for reading data on file. Note,
!--                the size of allocated surface elements do not necessarily
!--                need to match the size of present surface elements on 
!--                current processor, as the number of processors between 
!--                restarts can change. 
                   DO  l = 0, 3
                      ALLOCATE( tmp_surf_v(l)%t(1:ns_v_on_file_usm(l)) )
                      ALLOCATE( tmp_wall_v(l)%t(nzb_wall:nzt_wall+1,           &
                                                1:ns_v_on_file_usm(l) ) )
                      ALLOCATE( tmp_surf_window_v(l)%t(1:ns_v_on_file_usm(l)) )
                      ALLOCATE( tmp_window_v(l)%t(nzb_wall:nzt_wall+1,         & 
                                                  1:ns_v_on_file_usm(l) ) )
                      ALLOCATE( tmp_surf_green_v(l)%t(1:ns_v_on_file_usm(l)) )
                      ALLOCATE( tmp_green_v(l)%t(nzb_wall:nzt_wall+1,          &
                                                 1:ns_v_on_file_usm(l) ) )
                   ENDDO

                ENDIF    
          
             CASE ( 'usm_start_index_h', 'usm_start_index_v'  )   
                IF ( k == 1 )  THEN

                   IF ( ALLOCATED( start_index_on_file ) )                     &
                      DEALLOCATE( start_index_on_file )

                   ALLOCATE ( start_index_on_file(nys_on_file:nyn_on_file,     &
                                                  nxl_on_file:nxr_on_file) )

                   READ ( 13 )  start_index_on_file

                ENDIF
                
             CASE ( 'usm_end_index_h', 'usm_end_index_v' )   
                IF ( k == 1 )  THEN

                   IF ( ALLOCATED( end_index_on_file ) )                       &
                      DEALLOCATE( end_index_on_file )

                   ALLOCATE ( end_index_on_file(nys_on_file:nyn_on_file,       &
                                                nxl_on_file:nxr_on_file) )

                   READ ( 13 )  end_index_on_file

                ENDIF
          
             CASE ( 't_surf_h' )
#if defined( __nopointer )                   
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_h ) )                         &
                      ALLOCATE( t_surf_h(1:surf_usm_h%ns) )
                   READ ( 13 )  tmp_surf_h
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_h, tmp_surf_h,                  &
                                        surf_usm_h%start_index,                &  
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else                  
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_h_1 ) )                       &
                      ALLOCATE( t_surf_h_1(1:surf_usm_h%ns) )
                   READ ( 13 )  tmp_surf_h
                ENDIF              
                CALL surface_restore_elements(                                 &
                                        t_surf_h_1, tmp_surf_h,                &
                                        surf_usm_h%start_index,                & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif

             CASE ( 't_surf_v(0)' )
#if defined( __nopointer )            
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_v(0)%t ) )                    &
                      ALLOCATE( t_surf_v(0)%t(1:surf_usm_v(0)%ns) )
                   READ ( 13 )  tmp_surf_v(0)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_v(0)%t, tmp_surf_v(0)%t,        &
                                        surf_usm_v(0)%start_index,             &
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else                      
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_v_1(0)%t ) )                  &
                      ALLOCATE( t_surf_v_1(0)%t(1:surf_usm_v(0)%ns) )
                   READ ( 13 )  tmp_surf_v(0)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_v_1(0)%t, tmp_surf_v(0)%t,      &
                                        surf_usm_v(0)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
                      
             CASE ( 't_surf_v(1)' )
#if defined( __nopointer )        
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_v(1)%t ) )                    &
                      ALLOCATE( t_surf_v(1)%t(1:surf_usm_v(1)%ns) )
                   READ ( 13 )  tmp_surf_v(1)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_v(1)%t, tmp_surf_v(1)%t,        &
                                        surf_usm_v(1)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )                 
#else                      
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_v_1(1)%t ) )                  &
                      ALLOCATE( t_surf_v_1(1)%t(1:surf_usm_v(1)%ns) )
                   READ ( 13 )  tmp_surf_v(1)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_v_1(1)%t, tmp_surf_v(1)%t,      &
                                        surf_usm_v(1)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif

             CASE ( 't_surf_v(2)' )
#if defined( __nopointer )          
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_v(2)%t ) )                    &
                      ALLOCATE( t_surf_v(2)%t(1:surf_usm_v(2)%ns) )
                   READ ( 13 )  tmp_surf_v(2)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_v(2)%t, tmp_surf_v(2)%t,        &
                                        surf_usm_v(2)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else                      
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_v_1(2)%t ) )                  &
                      ALLOCATE( t_surf_v_1(2)%t(1:surf_usm_v(2)%ns) )
                   READ ( 13 )  tmp_surf_v(2)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_v_1(2)%t, tmp_surf_v(2)%t,      &
                                        surf_usm_v(2)%start_index,             &  
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
                      
             CASE ( 't_surf_v(3)' )
#if defined( __nopointer )    
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_v(3)%t ) )                    &
                      ALLOCATE( t_surf_v(3)%t(1:surf_usm_v(3)%ns) )
                   READ ( 13 )  tmp_surf_v(3)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_v(3)%t, tmp_surf_v(3)%t,        &
                                        surf_usm_v(3)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else                      
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_v_1(3)%t ) )                  &
                      ALLOCATE( t_surf_v_1(3)%t(1:surf_usm_v(3)%ns) )
                   READ ( 13 )  tmp_surf_v(3)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_v_1(3)%t, tmp_surf_v(3)%t,      &
                                        surf_usm_v(3)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
             CASE ( 't_surf_green_h' )
#if defined( __nopointer )                   
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_green_h ) )                   &
                      ALLOCATE( t_surf_green_h(1:surf_usm_h%ns) )
                   READ ( 13 )  tmp_surf_green_h
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_green_h, tmp_surf_green_h,      &
                                        surf_usm_h%start_index,                & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else                      
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_green_h_1 ) )                 &
                      ALLOCATE( t_surf_green_h_1(1:surf_usm_h%ns) )
                   READ ( 13 )  tmp_surf_green_h
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_green_h_1, tmp_surf_green_h,    &
                                        surf_usm_h%start_index,                &  
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif

             CASE ( 't_surf_green_v(0)' )
#if defined( __nopointer )            
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_green_v(0)%t ) )              &
                      ALLOCATE( t_surf_green_v(0)%t(1:surf_usm_v(0)%ns) )
                   READ ( 13 )  tmp_surf_green_v(0)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_green_v(0)%t,                   &
                                        tmp_surf_green_v(0)%t,                 &
                                        surf_usm_v(0)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else                      
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_green_v_1(0)%t ) )            &
                      ALLOCATE( t_surf_green_v_1(0)%t(1:surf_usm_v(0)%ns) )
                   READ ( 13 )  tmp_surf_green_v(0)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_green_v_1(0)%t,                 &
                                        tmp_surf_green_v(0)%t,                 &
                                        surf_usm_v(0)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
                   
             CASE ( 't_surf_green_v(1)' )
#if defined( __nopointer )        
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_green_v(1)%t ) )              &
                      ALLOCATE( t_surf_green_v(1)%t(1:surf_usm_v(1)%ns) )
                   READ ( 13 )  tmp_surf_green_v(1)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_green_v(1)%t,                   &
                                        tmp_surf_green_v(1)%t,                 &
                                        surf_usm_v(1)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )                 
#else                      
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_green_v_1(1)%t ) )            &
                      ALLOCATE( t_surf_green_v_1(1)%t(1:surf_usm_v(1)%ns) )
                   READ ( 13 )  tmp_surf_green_v(1)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_green_v_1(1)%t,                 &
                                        tmp_surf_green_v(1)%t,                 &
                                        surf_usm_v(1)%start_index,             &  
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif

             CASE ( 't_surf_green_v(2)' )
#if defined( __nopointer )          
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_green_v(2)%t ) )              &
                      ALLOCATE( t_surf_green_v(2)%t(1:surf_usm_v(2)%ns) )
                   READ ( 13 )  tmp_surf_green_v(2)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_green_v(2)%t,                   & 
                                        tmp_surf_green_v(2)%t,                 &
                                        surf_usm_v(2)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else                      
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_green_v_1(2)%t ) )            &
                      ALLOCATE( t_surf_green_v_1(2)%t(1:surf_usm_v(2)%ns) )
                   READ ( 13 )  tmp_surf_green_v(2)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_green_v_1(2)%t,                 &
                                        tmp_surf_green_v(2)%t,                 &
                                        surf_usm_v(2)%start_index,             &  
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
                   
             CASE ( 't_surf_green_v(3)' )
#if defined( __nopointer )    
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_green_v(3)%t ) )              &
                      ALLOCATE( t_surf_green_v(3)%t(1:surf_usm_v(3)%ns) )
                   READ ( 13 )  tmp_surf_green_v(3)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_green_v(3)%t,                   &
                                        tmp_surf_green_v(3)%t,                 &
                                        surf_usm_v(3)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else                      
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_green_v_1(3)%t ) )            &
                      ALLOCATE( t_surf_green_v_1(3)%t(1:surf_usm_v(3)%ns) )
                   READ ( 13 )  tmp_surf_green_v(3)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_green_v_1(3)%t,                 & 
                                        tmp_surf_green_v(3)%t,                 &
                                        surf_usm_v(3)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
             CASE ( 't_surf_window_h' )
#if defined( __nopointer )                   
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_window_h ) )                  &
                      ALLOCATE( t_surf_window_h(1:surf_usm_h%ns) )
                   READ ( 13 )  tmp_surf_window_h
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_window_h, tmp_surf_window_h,    &
                                        surf_usm_h%start_index,                & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else                      
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_window_h_1 ) )                &
                      ALLOCATE( t_surf_window_h_1(1:surf_usm_h%ns) )
                   READ ( 13 )  tmp_surf_window_h
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_window_h_1,                     &
                                        tmp_surf_window_h,                     &
                                        surf_usm_h%start_index,                & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif

             CASE ( 't_surf_window_v(0)' )
#if defined( __nopointer )            
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_window_v(0)%t ) )             &
                      ALLOCATE( t_surf_window_v(0)%t(1:surf_usm_v(0)%ns) )
                   READ ( 13 )  tmp_surf_window_v(0)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_window_v(0)%t,                  &
                                        tmp_surf_window_v(0)%t,                &
                                        surf_usm_v(0)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else                      
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_window_v_1(0)%t ) )           &
                      ALLOCATE( t_surf_window_v_1(0)%t(1:surf_usm_v(0)%ns) )
                   READ ( 13 )  tmp_surf_window_v(0)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_window_v_1(0)%t,                &
                                        tmp_surf_window_v(0)%t,                &
                                        surf_usm_v(0)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
                   
             CASE ( 't_surf_window_v(1)' )
#if defined( __nopointer )        
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_window_v(1)%t ) )             &
                      ALLOCATE( t_surf_window_v(1)%t(1:surf_usm_v(1)%ns) )
                   READ ( 13 )  tmp_surf_window_v(1)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_window_v(1)%t,                  & 
                                        tmp_surf_window_v(1)%t,                &
                                        surf_usm_v(1)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )                 
#else                      
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_window_v_1(1)%t ) )           &
                      ALLOCATE( t_surf_window_v_1(1)%t(1:surf_usm_v(1)%ns) )
                   READ ( 13 )  tmp_surf_window_v(1)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_window_v_1(1)%t,                &
                                        tmp_surf_window_v(1)%t,                &
                                        surf_usm_v(1)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif

             CASE ( 't_surf_window_v(2)' )
#if defined( __nopointer )          
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_window_v(2)%t ) )             &
                      ALLOCATE( t_surf_window_v(2)%t(1:surf_usm_v(2)%ns) )
                   READ ( 13 )  tmp_surf_window_v(2)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_window_v(2)%t,                  &
                                        tmp_surf_window_v(2)%t,                &
                                        surf_usm_v(2)%start_index,             &   
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else                      
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_window_v_1(2)%t ) )           &
                      ALLOCATE( t_surf_window_v_1(2)%t(1:surf_usm_v(2)%ns) )
                   READ ( 13 )  tmp_surf_window_v(2)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_window_v_1(2)%t,                & 
                                        tmp_surf_window_v(2)%t,                &
                                        surf_usm_v(2)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
                   
             CASE ( 't_surf_window_v(3)' )
#if defined( __nopointer )    
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_window_v(3)%t ) )             &
                      ALLOCATE( t_surf_window_v(3)%t(1:surf_usm_v(3)%ns) )
                   READ ( 13 )  tmp_surf_window_v(3)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_window_v(3)%t,                  &
                                        tmp_surf_window_v(3)%t,                &
                                        surf_usm_v(3)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else                      
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_surf_window_v_1(3)%t ) )           &
                      ALLOCATE( t_surf_window_v_1(3)%t(1:surf_usm_v(3)%ns) )
                   READ ( 13 )  tmp_surf_window_v(3)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_surf_window_v_1(3)%t,                &  
                                        tmp_surf_window_v(3)%t,                &
                                        surf_usm_v(3)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
             CASE ( 't_wall_h' )
#if defined( __nopointer )
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_wall_h ) )                         &
                      ALLOCATE( t_wall_h(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) )
                   READ ( 13 )  tmp_wall_h
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_wall_h, tmp_wall_h,                  &
                                        surf_usm_h%start_index,                & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_wall_h_1 ) )                       &
                      ALLOCATE( t_wall_h_1(nzb_wall:nzt_wall+1,                &
                                           1:surf_usm_h%ns) )
                   READ ( 13 )  tmp_wall_h
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_wall_h_1, tmp_wall_h,                &
                                        surf_usm_h%start_index,                & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
             CASE ( 't_wall_v(0)' )
#if defined( __nopointer )
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_wall_v(0)%t ) )                    &
                      ALLOCATE( t_wall_v(0)%t(nzb_wall:nzt_wall+1,             &
                                              1:surf_usm_v(0)%ns) )
                   READ ( 13 )  tmp_wall_v(0)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_wall_v(0)%t, tmp_wall_v(0)%t,        &
                                        surf_usm_v(0)%start_index,             &    
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_wall_v_1(0)%t ) )                  &
                      ALLOCATE( t_wall_v_1(0)%t(nzb_wall:nzt_wall+1,           &
                                                1:surf_usm_v(0)%ns) )
                   READ ( 13 )  tmp_wall_v(0)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_wall_v_1(0)%t, tmp_wall_v(0)%t,      &
                                        surf_usm_v(0)%start_index,             &  
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
             CASE ( 't_wall_v(1)' )
#if defined( __nopointer )
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_wall_v(1)%t ) )                    &
                      ALLOCATE( t_wall_v(1)%t(nzb_wall:nzt_wall+1,             &
                                              1:surf_usm_v(1)%ns) )
                   READ ( 13 )  tmp_wall_v(1)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_wall_v(1)%t, tmp_wall_v(1)%t,        &
                                        surf_usm_v(1)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file ,                    &
                                        nxlc, nysc,                            &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file, nxr_on_file )
#else
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_wall_v_1(1)%t ) )                  &
                      ALLOCATE( t_wall_v_1(1)%t(nzb_wall:nzt_wall+1,           &
                                                1:surf_usm_v(1)%ns) )
                   READ ( 13 )  tmp_wall_v(1)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_wall_v_1(1)%t, tmp_wall_v(1)%t,      &
                                        surf_usm_v(1)%start_index,             &  
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
             CASE ( 't_wall_v(2)' )
#if defined( __nopointer )
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_wall_v(2)%t ) )                    &
                      ALLOCATE( t_wall_v(2)%t(nzb_wall:nzt_wall+1,             &
                                              1:surf_usm_v(2)%ns) )
                   READ ( 13 )  tmp_wall_v(2)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_wall_v(2)%t, tmp_wall_v(2)%t,        &
                                        surf_usm_v(2)%start_index,             &  
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_wall_v_1(2)%t ) )                  &
                      ALLOCATE( t_wall_v_1(2)%t(nzb_wall:nzt_wall+1,           &
                                                1:surf_usm_v(2)%ns) )
                   READ ( 13 )  tmp_wall_v(2)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_wall_v_1(2)%t, tmp_wall_v(2)%t,      &
                                        surf_usm_v(2)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file ,                    &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
             CASE ( 't_wall_v(3)' )
#if defined( __nopointer )
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_wall_v(3)%t ) )                    &
                      ALLOCATE( t_wall_v(3)%t(nzb_wall:nzt_wall+1,             &
                                              1:surf_usm_v(3)%ns) )
                   READ ( 13 )  tmp_wall_v(3)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_wall_v(3)%t, tmp_wall_v(3)%t,        &
                                        surf_usm_v(3)%start_index,             &   
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_wall_v_1(3)%t ) )                  &
                      ALLOCATE( t_wall_v_1(3)%t(nzb_wall:nzt_wall+1,           &
                                                1:surf_usm_v(3)%ns) )
                   READ ( 13 )  tmp_wall_v(3)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_wall_v_1(3)%t, tmp_wall_v(3)%t,      &
                                        surf_usm_v(3)%start_index,             &   
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
             CASE ( 't_green_h' )
#if defined( __nopointer )
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_green_h ) )                        &
                      ALLOCATE( t_green_h(nzb_wall:nzt_wall+1,                 &
                                          1:surf_usm_h%ns) )
                   READ ( 13 )  tmp_green_h
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_green_h, tmp_green_h,                &
                                        surf_usm_h%start_index,                & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_green_h_1 ) )                      &
                      ALLOCATE( t_green_h_1(nzb_wall:nzt_wall+1,               &
                                            1:surf_usm_h%ns) )
                   READ ( 13 )  tmp_green_h
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_green_h_1, tmp_green_h,              &
                                        surf_usm_h%start_index,                & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
             CASE ( 't_green_v(0)' )
#if defined( __nopointer )
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_green_v(0)%t ) )                   &
                      ALLOCATE( t_green_v(0)%t(nzb_wall:nzt_wall+1,            &
                                               1:surf_usm_v(0)%ns) )
                   READ ( 13 )  tmp_green_v(0)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_green_v(0)%t, tmp_green_v(0)%t,      &
                                        surf_usm_v(0)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_green_v_1(0)%t ) )                 &
                      ALLOCATE( t_green_v_1(0)%t(nzb_wall:nzt_wall+1,          &
                                                 1:surf_usm_v(0)%ns) )
                   READ ( 13 )  tmp_green_v(0)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_green_v_1(0)%t, tmp_green_v(0)%t,    &
                                        surf_usm_v(0)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
             CASE ( 't_green_v(1)' )
#if defined( __nopointer )
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_green_v(1)%t ) )                   &
                      ALLOCATE( t_green_v(1)%t(nzb_wall:nzt_wall+1,            &
                                               1:surf_usm_v(1)%ns) )
                   READ ( 13 )  tmp_green_v(1)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_green_v(1)%t, tmp_green_v(1)%t,      &
                                        surf_usm_v(1)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file ,                    &
                                        nxlc, nysc,                            &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_green_v_1(1)%t ) )                 &
                      ALLOCATE( t_green_v_1(1)%t(nzb_wall:nzt_wall+1,          &
                                                 1:surf_usm_v(1)%ns) )
                   READ ( 13 )  tmp_green_v(1)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_green_v_1(1)%t, tmp_green_v(1)%t,    &
                                        surf_usm_v(1)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
             CASE ( 't_green_v(2)' )
#if defined( __nopointer )
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_green_v(2)%t ) )                   &
                      ALLOCATE( t_green_v(2)%t(nzb_wall:nzt_wall+1,            &
                                               1:surf_usm_v(2)%ns) )
                   READ ( 13 )  tmp_green_v(2)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_green_v(2)%t, tmp_green_v(2)%t,      &
                                        surf_usm_v(2)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_green_v_1(2)%t ) )                 &
                      ALLOCATE( t_green_v_1(2)%t(nzb_wall:nzt_wall+1,          &
                                                 1:surf_usm_v(2)%ns) )
                   READ ( 13 )  tmp_green_v(2)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_green_v_1(2)%t, tmp_green_v(2)%t,    &
                                        surf_usm_v(2)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file ,                    &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
             CASE ( 't_green_v(3)' )
#if defined( __nopointer )
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_green_v(3)%t ) )                   &
                      ALLOCATE( t_green_v(3)%t(nzb_wall:nzt_wall+1,            &
                                               1:surf_usm_v(3)%ns) )
                   READ ( 13 )  tmp_green_v(3)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_green_v(3)%t, tmp_green_v(3)%t,      &
                                        surf_usm_v(3)%start_index,             &  
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_green_v_1(3)%t ) )                 &
                      ALLOCATE( t_green_v_1(3)%t(nzb_wall:nzt_wall+1,          &
                                                 1:surf_usm_v(3)%ns) )
                   READ ( 13 )  tmp_green_v(3)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_green_v_1(3)%t, tmp_green_v(3)%t,    &
                                        surf_usm_v(3)%start_index,             &  
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
             CASE ( 't_window_h' )
#if defined( __nopointer )
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_window_h ) )                       &
                      ALLOCATE( t_window_h(nzb_wall:nzt_wall+1,                &
                                           1:surf_usm_h%ns) )
                   READ ( 13 )  tmp_window_h
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_window_h, tmp_window_h,              &
                                        surf_usm_h%start_index,                & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_window_h_1 ) )                     &
                      ALLOCATE( t_window_h_1(nzb_wall:nzt_wall+1,              &
                                             1:surf_usm_h%ns) )
                   READ ( 13 )  tmp_window_h
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_window_h_1, tmp_window_h,            &
                                        surf_usm_h%start_index,                & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file, nxr_on_file )
#endif
             CASE ( 't_window_v(0)' )
#if defined( __nopointer )
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_window_v(0)%t ) )                  &
                      ALLOCATE( t_window_v(0)%t(nzb_wall:nzt_wall+1,           &
                                                1:surf_usm_v(0)%ns) )
                   READ ( 13 )  tmp_window_v(0)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_window_v(0)%t, tmp_window_v(0)%t,    &
                                        surf_usm_v(0)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file, nxr_on_file )
#else
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_window_v_1(0)%t ) )                &
                      ALLOCATE( t_window_v_1(0)%t(nzb_wall:nzt_wall+1,         &
                                                  1:surf_usm_v(0)%ns) )
                   READ ( 13 )  tmp_window_v(0)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_window_v_1(0)%t,                     & 
                                        tmp_window_v(0)%t,                     &
                                        surf_usm_v(0)%start_index,             &
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
             CASE ( 't_window_v(1)' )
#if defined( __nopointer )
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_window_v(1)%t ) )                  &
                      ALLOCATE( t_window_v(1)%t(nzb_wall:nzt_wall+1,           &
                                                1:surf_usm_v(1)%ns) )
                   READ ( 13 )  tmp_window_v(1)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_window_v(1)%t, tmp_window_v(1)%t,    &
                                        surf_usm_v(1)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file ,                    &
                                        nxlc, nysc,                            &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file, nxr_on_file )
#else
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_window_v_1(1)%t ) )                &
                      ALLOCATE( t_window_v_1(1)%t(nzb_wall:nzt_wall+1,         &
                                                  1:surf_usm_v(1)%ns) )
                   READ ( 13 )  tmp_window_v(1)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_window_v_1(1)%t,                     & 
                                        tmp_window_v(1)%t,                     &
                                        surf_usm_v(1)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
             CASE ( 't_window_v(2)' )
#if defined( __nopointer )
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_window_v(2)%t ) )                  &
                      ALLOCATE( t_window_v(2)%t(nzb_wall:nzt_wall+1,           &
                                                1:surf_usm_v(2)%ns) )
                   READ ( 13 )  tmp_window_v(2)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_window_v(2)%t, tmp_window_v(2)%t,    &
                                        surf_usm_v(2)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_window_v_1(2)%t ) )                &
                      ALLOCATE( t_window_v_1(2)%t(nzb_wall:nzt_wall+1,         &
                                                  1:surf_usm_v(2)%ns) )
                   READ ( 13 )  tmp_window_v(2)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_window_v_1(2)%t,                     & 
                                        tmp_window_v(2)%t,                     &
                                        surf_usm_v(2)%start_index,             &  
                                        start_index_on_file,                   &
                                        end_index_on_file ,                    &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
             CASE ( 't_window_v(3)' )
#if defined( __nopointer )
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_window_v(3)%t ) )                  &
                      ALLOCATE( t_window_v(3)%t(nzb_wall:nzt_wall+1,           &
                                                1:surf_usm_v(3)%ns) )
                   READ ( 13 )  tmp_window_v(3)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_window_v(3)%t, tmp_window_v(3)%t,    &
                                        surf_usm_v(3)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#else
                IF ( k == 1 )  THEN
                   IF ( .NOT.  ALLOCATED( t_window_v_1(3)%t ) )                &
                      ALLOCATE( t_window_v_1(3)%t(nzb_wall:nzt_wall+1,1:surf_usm_v(3)%ns) )
                   READ ( 13 )  tmp_window_v(3)%t
                ENDIF
                CALL surface_restore_elements(                                 &
                                        t_window_v_1(3)%t,                     & 
                                        tmp_window_v(3)%t,                     &
                                        surf_usm_v(3)%start_index,             & 
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )
#endif
             CASE DEFAULT

                   found = .FALSE.

          END SELECT

       
    END SUBROUTINE usm_rrd_local
    

    
!------------------------------------------------------------------------------!
! Description:
! ------------
!
!> This subroutine reads walls, roofs and land categories and it parameters
!> from input files.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_read_urban_surface_types
    
        USE netcdf_data_input_mod,                                             &
            ONLY:  building_pars_f, building_type_f

        IMPLICIT NONE

        CHARACTER(12)                                         :: wtn
        INTEGER(iwp)                                          :: wtc
        REAL(wp), DIMENSION(n_surface_params)                 :: wtp
    
        INTEGER(iwp), DIMENSION(0:17, nysg:nyng, nxlg:nxrg)   :: usm_par
        REAL(wp), DIMENSION(1:14, nysg:nyng, nxlg:nxrg)       :: usm_val
        INTEGER(iwp)                                          :: k, l, d, iw, jw, kw, it, ip, ii, ij, m
        INTEGER(iwp)                                          :: i, j
        INTEGER(iwp)                                          :: nz, roof, dirwe, dirsn
        INTEGER(iwp)                                          :: category
        INTEGER(iwp)                                          :: weheight1, wecat1, snheight1, sncat1
        INTEGER(iwp)                                          :: weheight2, wecat2, snheight2, sncat2
        INTEGER(iwp)                                          :: weheight3, wecat3, snheight3, sncat3
        REAL(wp)                                              :: height, albedo, thick
        REAL(wp)                                              :: wealbedo1, wethick1, snalbedo1, snthick1
        REAL(wp)                                              :: wealbedo2, wethick2, snalbedo2, snthick2
        REAL(wp)                                              :: wealbedo3, wethick3, snalbedo3, snthick3
        
        LOGICAL                                               ::  surfpar
        LOGICAL                                               ::  urbsurf

!
!--     If building_pars or building_type are already read from static input 
!--     file, skip reading ASCII file. 
        IF ( building_type_f%from_file  .OR.  building_pars_f%from_file )      &
           RETURN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--     read categories of walls and their parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO  ii = 0, io_blocks-1
            IF ( ii == io_group )  THEN

!--             open urban surface file 
                OPEN( 151, file='SURFACE_PARAMETERS'//coupling_char, action='read', &
                           status='old', form='formatted', err=15 ) 
!--             first test and get n_surface_types
                k = 0
                l = 0
                DO
                    l = l+1
                    READ( 151, *, err=11, end=12 )  wtc, wtp, wtn
                    k = k+1
                    CYCLE
 11                 CONTINUE
                ENDDO
 12             n_surface_types = k
                ALLOCATE( surface_type_names(n_surface_types) )
                ALLOCATE( surface_type_codes(n_surface_types) )
                ALLOCATE( surface_params(n_surface_params, n_surface_types) )
!--             real reading
                rewind( 151 )
                k = 0
                DO
                    READ( 151, *, err=13, end=14 )  wtc, wtp, wtn
                    k = k+1
                    surface_type_codes(k) = wtc
                    surface_params(:,k) = wtp
                    surface_type_names(k) = wtn
                    CYCLE
13                  WRITE(6,'(i3,a,2i5)') myid, 'readparams2 error k=', k
                    FLUSH(6)
                    CONTINUE
                ENDDO
 14             CLOSE(151)
                CYCLE
 15             message_string = 'file SURFACE_PARAMETERS'//TRIM(coupling_char)//' does not exist'
                CALL message( 'usm_read_urban_surface_types', 'PA0513', 1, 2, 0, 6, 0 )
            ENDIF
        ENDDO
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--     read types of surfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        usm_par = 0
        DO  ii = 0, io_blocks-1
            IF ( ii == io_group )  THEN

                !
!--             open csv urban surface file 
                OPEN( 151, file='URBAN_SURFACE'//TRIM(coupling_char), action='read', &
                      status='old', form='formatted', err=23 )
                
                l = 0
                DO
                    l = l+1
!--                 i, j, height, nz, roof, dirwe, dirsn, category, soilcat,
!--                 weheight1, wecat1, snheight1, sncat1, weheight2, wecat2, snheight2, sncat2,
!--                 weheight3, wecat3, snheight3, sncat3
                    READ( 151, *, err=21, end=25 )  i, j, height, nz, roof, dirwe, dirsn,            &
                                            category, albedo, thick,                                 &
                                            weheight1, wecat1, wealbedo1, wethick1,                  &
                                            weheight2, wecat2, wealbedo2, wethick2,                  &
                                            weheight3, wecat3, wealbedo3, wethick3,                  &
                                            snheight1, sncat1, snalbedo1, snthick1,                  &
                                            snheight2, sncat2, snalbedo2, snthick2,                  &
                                            snheight3, sncat3, snalbedo3, snthick3

                    IF ( i >= nxlg  .AND.  i <= nxrg  .AND.  j >= nysg  .AND.  j <= nyng )  THEN
!--                     write integer variables into array
                        usm_par(:,j,i) = (/1, nz, roof, dirwe, dirsn, category,                      &
                                          weheight1, wecat1, weheight2, wecat2, weheight3, wecat3,   &
                                          snheight1, sncat1, snheight2, sncat2, snheight3, sncat3 /)
!--                     write real values into array
                        usm_val(:,j,i) = (/ albedo, thick,                                           &
                                           wealbedo1, wethick1, wealbedo2, wethick2,                 &
                                           wealbedo3, wethick3, snalbedo1, snthick1,                 &
                                           snalbedo2, snthick2, snalbedo3, snthick3 /)
                    ENDIF
                    CYCLE
 21                 WRITE (message_string, "(A,I5)") 'errors in file URBAN_SURFACE'//TRIM(coupling_char)//' on line ', l
                    CALL message( 'usm_read_urban_surface_types', 'PA0512', 0, 1, 0, 6, 0 )
                ENDDO
          
 23             message_string = 'file URBAN_SURFACE'//TRIM(coupling_char)//' does not exist'
                CALL message( 'usm_read_urban_surface_types', 'PA0514', 1, 2, 0, 6, 0 )

 25             CLOSE( 151 )

            ENDIF
#if defined( __parallel ) && ! defined ( __check )
            CALL MPI_BARRIER( comm2d, ierr )
#endif
        ENDDO
        
        !
!--     check completeness and formal correctness of the data
        DO i = nxlg, nxrg
            DO j = nysg, nyng
                IF ( usm_par(0,j,i) /= 0  .AND.  (        &  !< incomplete data,supply default values later
                     usm_par(1,j,i) < nzb  .OR.           &
                     usm_par(1,j,i) > nzt  .OR.           &  !< incorrect height (nz < nzb  .OR.  nz > nzt)
                     usm_par(2,j,i) < 0  .OR.             &
                     usm_par(2,j,i) > 1  .OR.             &  !< incorrect roof sign
                     usm_par(3,j,i) < nzb-nzt  .OR.       & 
                     usm_par(3,j,i) > nzt-nzb  .OR.       &  !< incorrect west-east wall direction sign
                     usm_par(4,j,i) < nzb-nzt  .OR.       &
                     usm_par(4,j,i) > nzt-nzb  .OR.       &  !< incorrect south-north wall direction sign
                     usm_par(6,j,i) < nzb  .OR.           & 
                     usm_par(6,j,i) > nzt  .OR.           &  !< incorrect pedestrian level height for west-east wall
                     usm_par(8,j,i) > nzt  .OR.           &
                     usm_par(10,j,i) > nzt  .OR.          &  !< incorrect wall or roof level height for west-east wall
                     usm_par(12,j,i) < nzb  .OR.          & 
                     usm_par(12,j,i) > nzt  .OR.          &  !< incorrect pedestrian level height for south-north wall
                     usm_par(14,j,i) > nzt  .OR.          &
                     usm_par(16,j,i) > nzt                &  !< incorrect wall or roof level height for south-north wall
                    ) )  THEN
!--                 incorrect input data
                    WRITE (message_string, "(A,2I5)") 'missing or incorrect data in file URBAN_SURFACE'// &
                                                       TRIM(coupling_char)//' for i,j=', i,j
                    CALL message( 'usm_read_urban_surface', 'PA0504', 1, 2, 0, 6, 0 )
                ENDIF
                
            ENDDO
        ENDDO
!        
!--     Assign the surface types to the respective data type. 
!--     First, for horizontal upward-facing surfaces. 
        DO  m = 1, surf_usm_h%ns
           iw = surf_usm_h%i(m)
           jw = surf_usm_h%j(m)
           kw = surf_usm_h%k(m)

           IF ( usm_par(5,jw,iw) == 0 )  THEN
#if ! defined( __nopointer )
              IF ( zu(kw) >= roof_height_limit )  THEN
                 surf_usm_h%isroof_surf(m)   = .TRUE.
                 surf_usm_h%surface_types(m) = roof_category         !< default category for root surface
              ELSE
                 surf_usm_h%isroof_surf(m)   = .FALSE.
                 surf_usm_h%surface_types(m) = land_category         !< default category for land surface
              ENDIF
#endif
              surf_usm_h%albedo(:,m)    = -1.0_wp
              surf_usm_h%thickness_wall(m) = -1.0_wp
              surf_usm_h%thickness_green(m) = -1.0_wp
              surf_usm_h%thickness_window(m) = -1.0_wp
           ELSE
              IF ( usm_par(2,jw,iw)==0 )  THEN
                 surf_usm_h%isroof_surf(m)    = .FALSE.
                 surf_usm_h%thickness_wall(m) = -1.0_wp
                 surf_usm_h%thickness_window(m) = -1.0_wp
                 surf_usm_h%thickness_green(m)  = -1.0_wp
              ELSE
                 surf_usm_h%isroof_surf(m)    = .TRUE.
                 surf_usm_h%thickness_wall(m) = usm_val(2,jw,iw)
                 surf_usm_h%thickness_window(m) = usm_val(2,jw,iw)
                 surf_usm_h%thickness_green(m)  = usm_val(2,jw,iw)
              ENDIF
              surf_usm_h%surface_types(m) = usm_par(5,jw,iw)
              surf_usm_h%albedo(:,m)   = usm_val(1,jw,iw)
              surf_usm_h%transmissivity(m)    = 0.0_wp
           ENDIF
!
!--        Find the type position
           it = surf_usm_h%surface_types(m)
           ip = -99999
           DO k = 1, n_surface_types
              IF ( surface_type_codes(k) == it )  THEN
                 ip = k
                 EXIT
              ENDIF
           ENDDO
           IF ( ip == -99999 )  THEN
!--           wall category not found
              WRITE (message_string, "(A,I5,A,3I5)") 'wall category ', it,     &
                                     ' not found  for i,j,k=', iw,jw,kw
              CALL message( 'usm_read_urban_surface', 'PA0506', 1, 2, 0, 6, 0 )
           ENDIF
!
!--        Albedo
           IF ( surf_usm_h%albedo(ind_veg_wall,m) < 0.0_wp )  THEN
              surf_usm_h%albedo(:,m) = surface_params(ialbedo,ip)
           ENDIF
!--        Albedo type is 0 (custom), others are replaced later
           surf_usm_h%albedo_type(:,m) = 0
!--        Transmissivity
           IF ( surf_usm_h%transmissivity(m) < 0.0_wp )  THEN
              surf_usm_h%transmissivity(m) = 0.0_wp
           ENDIF
!
!--        emissivity of the wall
           surf_usm_h%emissivity(:,m) = surface_params(iemiss,ip)
!            
!--        heat conductivity λS between air and wall ( W m−2 K−1 )
           surf_usm_h%lambda_surf(m) = surface_params(ilambdas,ip)
           surf_usm_h%lambda_surf_window(m) = surface_params(ilambdas,ip)
           surf_usm_h%lambda_surf_green(m)  = surface_params(ilambdas,ip)
!            
!--        roughness length for momentum, heat and humidity
           surf_usm_h%z0(m) = surface_params(irough,ip)
           surf_usm_h%z0h(m) = surface_params(iroughh,ip)
           surf_usm_h%z0q(m) = surface_params(iroughh,ip)
!
!--        Surface skin layer heat capacity (J m−2 K−1 )
           surf_usm_h%c_surface(m) = surface_params(icsurf,ip)
           surf_usm_h%c_surface_window(m) = surface_params(icsurf,ip)
           surf_usm_h%c_surface_green(m)  = surface_params(icsurf,ip)
!            
!--        wall material parameters:
!--        thickness of the wall (m)
!--        missing values are replaced by default value for category
           IF ( surf_usm_h%thickness_wall(m) <= 0.001_wp )  THEN
                surf_usm_h%thickness_wall(m) = surface_params(ithick,ip)
           ENDIF
           IF ( surf_usm_h%thickness_window(m) <= 0.001_wp )  THEN
                surf_usm_h%thickness_window(m) = surface_params(ithick,ip)
           ENDIF
           IF ( surf_usm_h%thickness_green(m) <= 0.001_wp )  THEN
                surf_usm_h%thickness_green(m) = surface_params(ithick,ip)
           ENDIF
!            
!--        volumetric heat capacity rho*C of the wall ( J m−3 K−1 )
           surf_usm_h%rho_c_wall(:,m) = surface_params(irhoC,ip)
           surf_usm_h%rho_c_window(:,m) = surface_params(irhoC,ip)
           surf_usm_h%rho_c_green(:,m)  = surface_params(irhoC,ip)
!            
!--        thermal conductivity λH of the wall (W m−1 K−1 )
           surf_usm_h%lambda_h(:,m) = surface_params(ilambdah,ip)
           surf_usm_h%lambda_h_window(:,m) = surface_params(ilambdah,ip)
           surf_usm_h%lambda_h_green(:,m)  = surface_params(ilambdah,ip)

        ENDDO
!
!--     For vertical surface elements ( 0 -- northward-facing, 1 -- southward-facing,
!--     2 -- eastward-facing, 3 -- westward-facing )
        DO  l = 0, 3
           DO  m = 1, surf_usm_v(l)%ns
              i  = surf_usm_v(l)%i(m)
              j  = surf_usm_v(l)%j(m)
              kw = surf_usm_v(l)%k(m)
              
              IF ( l == 3 )  THEN ! westward facing
                 iw = i
                 jw = j
                 ii = 6
                 ij = 3
              ELSEIF ( l == 2 )  THEN
                 iw = i-1
                 jw = j
                 ii = 6
                 ij = 3
              ELSEIF ( l == 1 )  THEN
                 iw = i
                 jw = j
                 ii = 12
                 ij = 9
              ELSEIF ( l == 0 )  THEN
                 iw = i
                 jw = j-1
                 ii = 12
                 ij = 9
              ENDIF

              IF ( iw < 0 .OR. jw < 0 ) THEN
!--              wall on west or south border of the domain - assign default category
                 IF ( kw <= roof_height_limit ) THEN
                     surf_usm_v(l)%surface_types(m) = wall_category   !< default category for wall surface in wall zone
                 ELSE
                     surf_usm_v(l)%surface_types(m) = roof_category   !< default category for wall surface in roof zone
                 END IF
                 surf_usm_v(l)%albedo(:,m)    = -1.0_wp
                 surf_usm_v(l)%thickness_wall(m) = -1.0_wp
              ELSE IF ( kw <= usm_par(ii,jw,iw) )  THEN
!--                 pedestrian zone
                 IF ( usm_par(ii+1,jw,iw) == 0 )  THEN
                     surf_usm_v(l)%surface_types(m)  = pedestrian_category   !< default category for wall surface in pedestrian zone
                     surf_usm_v(l)%albedo(:,m)    = -1.0_wp
                     surf_usm_v(l)%thickness_wall(m) = -1.0_wp
                     surf_usm_v(l)%thickness_window(m)   = -1.0_wp
                     surf_usm_v(l)%thickness_green(m)    = -1.0_wp
                     surf_usm_v(l)%transmissivity(m)  = -1.0_wp
                 ELSE
                     surf_usm_v(l)%surface_types(m)  = usm_par(ii+1,jw,iw)
                     surf_usm_v(l)%albedo(:,m)    = usm_val(ij,jw,iw)
                     surf_usm_v(l)%thickness_wall(m) = usm_val(ij+1,jw,iw)
                     surf_usm_v(l)%thickness_window(m)   = usm_val(ij+1,jw,iw)
                     surf_usm_v(l)%thickness_green(m)    = usm_val(ij+1,jw,iw)
                     surf_usm_v(l)%transmissivity(m)  = 0.0_wp
                 ENDIF
              ELSE IF ( kw <= usm_par(ii+2,jw,iw) )  THEN
!--              wall zone
                 IF ( usm_par(ii+3,jw,iw) == 0 )  THEN
                     surf_usm_v(l)%surface_types(m)  = wall_category         !< default category for wall surface
                     surf_usm_v(l)%albedo(:,m)    = -1.0_wp
                     surf_usm_v(l)%thickness_wall(m) = -1.0_wp
                     surf_usm_v(l)%thickness_window(m)   = -1.0_wp
                     surf_usm_v(l)%thickness_green(m)    = -1.0_wp
                     surf_usm_v(l)%transmissivity(m)  = -1.0_wp
                 ELSE
                     surf_usm_v(l)%surface_types(m)  = usm_par(ii+3,jw,iw)
                     surf_usm_v(l)%albedo(:,m)    = usm_val(ij+2,jw,iw)
                     surf_usm_v(l)%thickness_wall(m) = usm_val(ij+3,jw,iw)
                     surf_usm_v(l)%thickness_window(m)   = usm_val(ij+3,jw,iw)
                     surf_usm_v(l)%thickness_green(m)    = usm_val(ij+3,jw,iw)
                     surf_usm_v(l)%transmissivity(m)  = 0.0_wp
                 ENDIF
              ELSE IF ( kw <= usm_par(ii+4,jw,iw) )  THEN
!--              roof zone
                 IF ( usm_par(ii+5,jw,iw) == 0 )  THEN
                     surf_usm_v(l)%surface_types(m)  = roof_category         !< default category for roof surface
                     surf_usm_v(l)%albedo(:,m)    = -1.0_wp
                     surf_usm_v(l)%thickness_wall(m) = -1.0_wp
                     surf_usm_v(l)%thickness_window(m)   = -1.0_wp
                     surf_usm_v(l)%thickness_green(m)    = -1.0_wp
                     surf_usm_v(l)%transmissivity(m)  = -1.0_wp
                 ELSE
                     surf_usm_v(l)%surface_types(m)  = usm_par(ii+5,jw,iw)
                     surf_usm_v(l)%albedo(:,m)    = usm_val(ij+4,jw,iw)
                     surf_usm_v(l)%thickness_wall(m) = usm_val(ij+5,jw,iw)
                     surf_usm_v(l)%thickness_window(m)   = usm_val(ij+5,jw,iw)
                     surf_usm_v(l)%thickness_green(m)    = usm_val(ij+5,jw,iw)
                     surf_usm_v(l)%transmissivity(m)  = 0.0_wp
                 ENDIF
              ELSE
!
!--              supply the default category
                 IF ( kw <= roof_height_limit ) THEN
                     surf_usm_v(l)%surface_types(m) = wall_category   !< default category for wall surface in wall zone
                 ELSE
                     surf_usm_v(l)%surface_types(m) = roof_category   !< default category for wall surface in roof zone
                 END IF
                 surf_usm_v(l)%albedo(:,m)    = -1.0_wp
                 surf_usm_v(l)%thickness_wall(m) = -1.0_wp
              ENDIF
!
!--           Find the type position
              it = surf_usm_v(l)%surface_types(m)
              ip = -99999
              DO k = 1, n_surface_types
                 IF ( surface_type_codes(k) == it )  THEN
                    ip = k
                    EXIT
                 ENDIF
              ENDDO
              IF ( ip == -99999 )  THEN
!--              wall category not found
                 WRITE (message_string, "(A,I7,A,3I5)") 'wall category ', it,  &
                                        ' not found  for i,j,k=', iw,jw,kw
                 WRITE(9,*) message_string
              ENDIF
!
!--           Albedo
              IF ( surf_usm_v(l)%albedo(ind_veg_wall,m) < 0.0_wp )  THEN
                 surf_usm_v(l)%albedo(:,m) = surface_params(ialbedo,ip)
              ENDIF
!--           Albedo type is 0 (custom), others are replaced later
              surf_usm_v(l)%albedo_type(:,m) = 0
!--           Transmissivity of the windows
              IF ( surf_usm_v(l)%transmissivity(m) < 0.0_wp )  THEN
                 surf_usm_v(l)%transmissivity(m) = 0.0_wp
              ENDIF
!
!--           emissivity of the wall
              surf_usm_v(l)%emissivity(:,m) = surface_params(iemiss,ip)
!            
!--           heat conductivity lambda S between air and wall ( W m-2 K-1 )
              surf_usm_v(l)%lambda_surf(m) = surface_params(ilambdas,ip)
              surf_usm_v(l)%lambda_surf_window(m) = surface_params(ilambdas,ip)
              surf_usm_v(l)%lambda_surf_green(m) = surface_params(ilambdas,ip)
!            
!--           roughness length
              surf_usm_v(l)%z0(m) = surface_params(irough,ip)
              surf_usm_v(l)%z0h(m) = surface_params(iroughh,ip)
              surf_usm_v(l)%z0q(m) = surface_params(iroughh,ip)
!            
!--           Surface skin layer heat capacity (J m-2 K-1 )
              surf_usm_v(l)%c_surface(m) = surface_params(icsurf,ip)
              surf_usm_v(l)%c_surface_window(m) = surface_params(icsurf,ip)
              surf_usm_v(l)%c_surface_green(m) = surface_params(icsurf,ip)
!            
!--           wall material parameters:
!--           thickness of the wall (m)
!--           missing values are replaced by default value for category
              IF ( surf_usm_v(l)%thickness_wall(m) <= 0.001_wp )  THEN
                   surf_usm_v(l)%thickness_wall(m) = surface_params(ithick,ip)
              ENDIF
              IF ( surf_usm_v(l)%thickness_window(m) <= 0.001_wp )  THEN
                   surf_usm_v(l)%thickness_window(m) = surface_params(ithick,ip)
              ENDIF
              IF ( surf_usm_v(l)%thickness_green(m) <= 0.001_wp )  THEN
                   surf_usm_v(l)%thickness_green(m) = surface_params(ithick,ip)
              ENDIF
!
!--           volumetric heat capacity rho*C of the wall ( J m-3 K-1 )
              surf_usm_v(l)%rho_c_wall(:,m) = surface_params(irhoC,ip)
              surf_usm_v(l)%rho_c_window(:,m) = surface_params(irhoC,ip)
              surf_usm_v(l)%rho_c_green(:,m) = surface_params(irhoC,ip)
!            
!--           thermal conductivity lambda H of the wall (W m-1 K-1 )
              surf_usm_v(l)%lambda_h(:,m) = surface_params(ilambdah,ip)
              surf_usm_v(l)%lambda_h_window(:,m) = surface_params(ilambdah,ip)
              surf_usm_v(l)%lambda_h_green(:,m) = surface_params(ilambdah,ip)

           ENDDO
        ENDDO 
!
!--     Initialize wall layer thicknesses. Please note, this will be removed
!--     after migration to Palm input data standard.  
        DO k = nzb_wall, nzt_wall
           zwn(k) = zwn_default(k)
           zwn_green(k) = zwn_default_green(k)
           zwn_window(k) = zwn_default_window(k)
        ENDDO
!
!--     apply for all particular surface grids. First for horizontal surfaces
        DO  m = 1, surf_usm_h%ns
           surf_usm_h%zw(:,m) = zwn(:) * surf_usm_h%thickness_wall(m)
           surf_usm_h%zw_green(:,m) = zwn_green(:) * surf_usm_h%thickness_green(m)
           surf_usm_h%zw_window(:,m) = zwn_window(:) * surf_usm_h%thickness_window(m)
        ENDDO
        DO  l = 0, 3
           DO  m = 1, surf_usm_v(l)%ns
              surf_usm_v(l)%zw(:,m) = zwn(:) * surf_usm_v(l)%thickness_wall(m)
              surf_usm_v(l)%zw_green(:,m) = zwn_green(:) * surf_usm_v(l)%thickness_green(m)
              surf_usm_v(l)%zw_window(:,m) = zwn_window(:) * surf_usm_v(l)%thickness_window(m)
           ENDDO
        ENDDO

        CALL location_message( '    types and parameters of urban surfaces read', .TRUE. )
   
    END SUBROUTINE usm_read_urban_surface_types


!------------------------------------------------------------------------------!
! Description:
! ------------
!
!> This function advances through the list of local surfaces to find given
!> x, y, d, z coordinates
!------------------------------------------------------------------------------!
    PURE FUNCTION advance_surface(isurfl_start, isurfl_stop, x, y, z, d) &
            result(isurfl)

        INTEGER(iwp), INTENT(in)                :: isurfl_start, isurfl_stop
        INTEGER(iwp), INTENT(in)                :: x, y, z, d
        INTEGER(iwp)                            :: isx, isy, isz, isd
        INTEGER(iwp)                            :: isurfl

        DO isurfl = isurfl_start, isurfl_stop
            isx = surfl(ix, isurfl)
            isy = surfl(iy, isurfl)
            isz = surfl(iz, isurfl)
            isd = surfl(id, isurfl)
            IF ( isx==x .and. isy==y .and. isz==z .and. isd==d )  RETURN
        ENDDO

!--     coordinate not found
        isurfl = -1

    END FUNCTION


!------------------------------------------------------------------------------!
! Description:
! ------------
!
!> This subroutine reads temperatures of respective material layers in walls,
!> roofs and ground from input files. Data in the input file must be in
!> standard order, i.e. horizontal surfaces first ordered by x, y and then
!> vertical surfaces ordered by x, y, direction, z
!------------------------------------------------------------------------------!
    SUBROUTINE usm_read_wall_temperature

        INTEGER(iwp)                                          :: i, j, k, d, ii, iline
        INTEGER(iwp)                                          :: isurfl
        REAL(wp)                                              :: rtsurf
        REAL(wp), DIMENSION(nzb_wall:nzt_wall+1)              :: rtwall




        DO  ii = 0, io_blocks-1
            IF ( ii == io_group )  THEN

!--             open wall temperature file
                OPEN( 152, file='WALL_TEMPERATURE'//coupling_char, action='read', &
                           status='old', form='formatted', err=15 )

                isurfl = 0
                iline = 1
                DO
                    rtwall = -9999.0_wp  !< for incomplete lines
                    READ( 152, *, err=13, end=14 )  i, j, k, d, rtsurf, rtwall

                    IF ( nxl <= i .and. i <= nxr .and. &
                        nys <= j .and. j <= nyn)  THEN  !< local processor
!--                     identify surface id
                        isurfl = advance_surface(isurfl+1, nsurfl, i, j, k, d)
                        IF ( isurfl == -1 )  THEN
                            WRITE(message_string, '(a,4i5,a,i5,a)') 'Coordinates (xyzd) ', i, j, k, d, &
                                ' on line ', iline, &
                                ' in file WALL_TEMPERATURE are either not present or out of standard order of surfaces.'
                            CALL message( 'usm_read_wall_temperature', 'PA0521', 1, 2, 0, 6, 0 )
                        ENDIF

!--                     assign temperatures
                        IF ( d == 0 ) THEN
                           t_surf_h(isurfl) = rtsurf
                           t_wall_h(:,isurfl) = rtwall(:)
                        ELSE
                           t_surf_v(d-1)%t(isurfl) = rtsurf
                           t_wall_v(d-1)%t(:,isurfl) = rtwall(:)
                        ENDIF
                    ENDIF

                    iline = iline + 1
                    CYCLE
 13                 WRITE(message_string, '(a,i5,a)') 'Error reading line ', iline, &
                        ' in file WALL_TEMPERATURE.'
                    CALL message( 'usm_read_wall_temperature', 'PA0522', 1, 2, 0, 6, 0 )
                ENDDO
 14             CLOSE(152)
                CYCLE
 15             message_string = 'file WALL_TEMPERATURE'//TRIM(coupling_char)//' does not exist'
                CALL message( 'usm_read_wall_temperature', 'PA0523', 1, 2, 0, 6, 0 )
            ENDIF
#if defined( __parallel ) && ! defined ( __check )
            CALL MPI_BARRIER( comm2d, ierr )
#endif
        ENDDO

        CALL location_message( '    wall layer temperatures read', .TRUE. )

    END SUBROUTINE usm_read_wall_temperature



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Solver for the energy balance at the ground/roof/wall surface.
!> It follows basic ideas and structure of lsm_energy_balance
!> with many simplifications and adjustments.
!> TODO better description
!------------------------------------------------------------------------------!
    SUBROUTINE usm_surface_energy_balance

        IMPLICIT NONE

        INTEGER(iwp)                          :: i, j, k, l, d, m   !< running indices
        
        REAL(wp)                              :: stend              !< surface tendency
        REAL(wp)                              :: stend_window       !< surface tendency
        REAL(wp)                              :: stend_green        !< surface tendency
        REAL(wp)                              :: coef_1             !< first coeficient for prognostic equation
        REAL(wp)                              :: coef_window_1      !< first coeficient for prognostic window equation
        REAL(wp)                              :: coef_green_1       !< first coeficient for prognostic green wall equation
        REAL(wp)                              :: coef_2             !< second  coeficient for prognostic equation
        REAL(wp)                              :: coef_window_2      !< second  coeficient for prognostic window equation
        REAL(wp)                              :: coef_green_2       !< second  coeficient for prognostic green wall equation
        REAL(wp)                              :: rho_cp             !< rho_wall_surface * cp
        REAL(wp)                              :: f_shf              !< factor for shf_eb
        REAL(wp)                              :: f_shf_window       !< factor for shf_eb window
        REAL(wp)                              :: f_shf_green        !< factor for shf_eb green wall
        REAL(wp)                              :: lambda_surface     !< current value of lambda_surface (heat conductivity between air and wall)
        REAL(wp)                              :: lambda_surface_window !< current value of lambda_surface (heat conductivity between air and window)
        REAL(wp)                              :: lambda_surface_green  !< current value of lambda_surface (heat conductivity between air and greeb wall)
        REAL(wp), DIMENSION(nzb:nzt)          :: exn                !< value of the Exner function in layers
        
        REAL(wp)                              :: dtime              !< simulated time of day (in UTC)
        INTEGER(iwp)                          :: dhour              !< simulated hour of day (in UTC)
        REAL(wp)                              :: acoef              !< actual coefficient of diurnal profile of anthropogenic heat


#if ! defined( __nopointer )
        exn(nzb:nzt) = (hyp(nzb:nzt) / 100000.0_wp )**0.286_wp          !< Exner function
#endif
!        
!--     First, treat horizontal surface elements
        DO  m = 1, surf_usm_h%ns
!
!--        Get indices of respective grid point
           i = surf_usm_h%i(m)
           j = surf_usm_h%j(m)
           k = surf_usm_h%k(m)
!
!--        TODO - how to calculate lambda_surface for horizontal surfaces
!--        (lambda_surface is set according to stratification in land surface model)
!--        MS: ???
           IF ( surf_usm_h%ol(m) >= 0.0_wp )  THEN
              lambda_surface = surf_usm_h%lambda_surf(m)
              lambda_surface_window = surf_usm_h%lambda_surf_window(m)
              lambda_surface_green = surf_usm_h%lambda_surf_green(m)
           ELSE
              lambda_surface = surf_usm_h%lambda_surf(m)
              lambda_surface_window = surf_usm_h%lambda_surf_window(m)
              lambda_surface_green = surf_usm_h%lambda_surf_green(m)
           ENDIF
#if ! defined( __nopointer )
!
!--        calculate rho * cp coefficient at surface layer
           rho_cp  = cp * hyp(k) / ( r_d * surf_usm_h%pt1(m) * exn(k) )
#endif
!
!--        Calculate aerodyamic resistance. 
!--        Calculation for horizontal surfaces follows LSM formulation
!--        pt, us, ts are not available for the prognostic time step,
!--        data from the last time step is used here.

!--        Workaround: use single r_a as stability is only treated for the
!--        average temperature
           surf_usm_h%r_a(m) = ( surf_usm_h%pt1(m) - surf_usm_h%pt_surface(m) ) /&
                               ( surf_usm_h%ts(m) * surf_usm_h%us(m) + 1.0E-20_wp )   
           surf_usm_h%r_a_window(m) = surf_usm_h%r_a(m)
           surf_usm_h%r_a_green(m)  = surf_usm_h%r_a(m)

!            r_a = ( surf_usm_h%pt1(m) - t_surf_h(m) / exn(k) ) /                              &
!                  ( surf_usm_h%ts(m) * surf_usm_h%us(m) + 1.0E-20_wp )
!            r_a_window = ( surf_usm_h%pt1(m) - t_surf_window_h(m) / exn(k) ) /                &
!                  ( surf_usm_h%ts(m) * surf_usm_h%us(m) + 1.0E-20_wp )
!            r_a_green = ( surf_usm_h%pt1(m) - t_surf_green_h(m) / exn(k) ) /                  &
!                  ( surf_usm_h%ts(m) * surf_usm_h%us(m) + 1.0E-20_wp )
                
!--        make sure that the resistance does not drop to zero
           IF ( surf_usm_h%r_a(m)        < 1.0_wp )                            &
               surf_usm_h%r_a(m)        = 1.0_wp
           IF ( surf_usm_h%r_a_green(m)  < 1.0_wp )                            &
               surf_usm_h%r_a_green(m) = 1.0_wp
           IF ( surf_usm_h%r_a_window(m) < 1.0_wp )                            &
               surf_usm_h%r_a_window(m) = 1.0_wp
                
!--        factor for shf_eb
           f_shf  = rho_cp / surf_usm_h%r_a(m)
           f_shf_window  = rho_cp / surf_usm_h%r_a_window(m)
           f_shf_green  = rho_cp / surf_usm_h%r_a_green(m)
        
!--        add LW up so that it can be removed in prognostic equation
           surf_usm_h%rad_net_l(m) = surf_usm_h%rad_sw_in(m)  -                &
                                     surf_usm_h%rad_sw_out(m) +                &
                                     surf_usm_h%rad_lw_in(m)  -                &
                                     surf_usm_h%rad_lw_out(m)

!--        numerator of the prognostic equation
!--     Todo: Adjust to tile approach. So far, emissivity for wall (element 0)
!--           is used
           coef_1 = surf_usm_h%rad_net_l(m) +                                  & 
                 ( 3.0_wp + 1.0_wp ) * surf_usm_h%emissivity(ind_veg_wall,m) * &
                                       sigma_sb * t_surf_h(m) ** 4 +           &  
                                       f_shf * surf_usm_h%pt1(m) +             &
                                       lambda_surface * t_wall_h(nzb_wall,m)
           coef_window_1 = surf_usm_h%rad_net_l(m) +                           & 
                   ( 3.0_wp + 1.0_wp ) * surf_usm_h%emissivity(ind_wat_win,m)  &
                                       * sigma_sb * t_surf_window_h(m) ** 4 +  &  
                                       f_shf_window * surf_usm_h%pt1(m) +      &
                                       lambda_surface_window * t_window_h(nzb_wall,m)
           coef_green_1 = surf_usm_h%rad_net_l(m) +                            & 
                 ( 3.0_wp + 1.0_wp ) * surf_usm_h%emissivity(ind_pav_green,m) *&
                                       sigma_sb * t_surf_green_h(m) ** 4 +     &  
                                       f_shf_green * surf_usm_h%pt1(m) +       &
                                       lambda_surface_green * t_wall_h(nzb_wall,m)

!--        denominator of the prognostic equation
           coef_2 = 4.0_wp * surf_usm_h%emissivity(ind_veg_wall,m) *           &
                             sigma_sb * t_surf_h(m) ** 3                       &
                           + lambda_surface + f_shf / exn(k)
           coef_window_2 = 4.0_wp * surf_usm_h%emissivity(ind_wat_win,m) *     &
                             sigma_sb * t_surf_window_h(m) ** 3                &
                           + lambda_surface_window + f_shf_window / exn(k)
           coef_green_2 = 4.0_wp * surf_usm_h%emissivity(ind_pav_green,m) *    &
                             sigma_sb * t_surf_green_h(m) ** 3                 &
                           + lambda_surface_green + f_shf_green / exn(k)

!--        implicit solution when the surface layer has no heat capacity,
!--        otherwise use RK3 scheme.
           t_surf_h_p(m) = ( coef_1 * dt_3d * tsc(2) +                        &
                             surf_usm_h%c_surface(m) * t_surf_h(m) ) /        & 
                           ( surf_usm_h%c_surface(m) + coef_2 * dt_3d * tsc(2) ) 
           t_surf_window_h_p(m) = ( coef_window_1 * dt_3d * tsc(2) +                        &
                             surf_usm_h%c_surface_window(m) * t_surf_window_h(m) ) /        & 
                           ( surf_usm_h%c_surface_window(m) + coef_window_2 * dt_3d * tsc(2) ) 
           t_surf_green_h_p(m) = ( coef_green_1 * dt_3d * tsc(2) +                        &
                             surf_usm_h%c_surface_green(m) * t_surf_green_h(m) ) /        & 
                           ( surf_usm_h%c_surface_green(m) + coef_green_2 * dt_3d * tsc(2) ) 

!--        add RK3 term
           t_surf_h_p(m) = t_surf_h_p(m) + dt_3d * tsc(3) *                   &
                           surf_usm_h%tt_surface_m(m)
           t_surf_window_h_p(m) = t_surf_window_h_p(m) + dt_3d * tsc(3) *     &
                           surf_usm_h%tt_surface_window_m(m)
           t_surf_green_h_p(m) = t_surf_green_h_p(m) + dt_3d * tsc(3) *       &
                           surf_usm_h%tt_surface_green_m(m)
!
!--        Store surface temperature        
           surf_usm_h%pt_surface(m) = ( surf_usm_h%frac(ind_veg_wall,m) * t_surf_h_p(m)   &
                               + surf_usm_h%frac(ind_wat_win,m) * t_surf_window_h_p(m)   &
                               + surf_usm_h%frac(ind_pav_green,m) * t_surf_green_h_p(m) )  &
                               / exn(k)
                       
!--        calculate true tendency
           stend = ( t_surf_h_p(m) - t_surf_h(m) - dt_3d * tsc(3) *           &
                     surf_usm_h%tt_surface_m(m)) / ( dt_3d  * tsc(2) )
           stend_window = ( t_surf_window_h_p(m) - t_surf_window_h(m) - dt_3d * tsc(3) *        &
                     surf_usm_h%tt_surface_window_m(m)) / ( dt_3d  * tsc(2) )
           stend_green = ( t_surf_green_h_p(m) - t_surf_green_h(m) - dt_3d * tsc(3) *           &
                     surf_usm_h%tt_surface_green_m(m)) / ( dt_3d  * tsc(2) )

!--        calculate t_surf tendencies for the next Runge-Kutta step
           IF ( timestep_scheme(1:5) == 'runge' )  THEN
              IF ( intermediate_timestep_count == 1 )  THEN
                 surf_usm_h%tt_surface_m(m) = stend
                 surf_usm_h%tt_surface_window_m(m) = stend_window
                 surf_usm_h%tt_surface_green_m(m) = stend_green
              ELSEIF ( intermediate_timestep_count <                          &
                        intermediate_timestep_count_max )  THEN
                 surf_usm_h%tt_surface_m(m) = -9.5625_wp * stend +            &
                                     5.3125_wp * surf_usm_h%tt_surface_m(m)
                 surf_usm_h%tt_surface_window_m(m) = -9.5625_wp * stend_window +   &
                                     5.3125_wp * surf_usm_h%tt_surface_window_m(m)
                 surf_usm_h%tt_surface_green_m(m) = -9.5625_wp * stend_green +     &
                                     5.3125_wp * surf_usm_h%tt_surface_green_m(m)
              ENDIF
           ENDIF

!--        in case of fast changes in the skin temperature, it is required to
!--        update the radiative fluxes in order to keep the solution stable
           IF ( ( ABS( t_surf_h_p(m) - t_surf_h(m) ) > 1.0_wp ) .OR. &
                ( ABS( t_surf_green_h_p(m) - t_surf_green_h(m) ) > 1.0_wp ) .OR. &
                ( ABS( t_surf_window_h_p(m) - t_surf_window_h(m) ) > 1.0_wp ) ) THEN
              force_radiation_call_l = .TRUE.
           ENDIF
!
!--        calculate fluxes
!--        rad_net_l is never used!           
           surf_usm_h%rad_net_l(m) = surf_usm_h%rad_net_l(m) +                           &
                                     surf_usm_h%frac(ind_veg_wall,m) *                   &
                                     sigma_sb * surf_usm_h%emissivity(ind_veg_wall,m) *  &
                                     ( t_surf_h_p(m)**4 - t_surf_h(m)**4 )               &
                                    + surf_usm_h%frac(ind_wat_win,m) *                   &
                                     sigma_sb * surf_usm_h%emissivity(ind_wat_win,m) *   &
                                     ( t_surf_window_h_p(m)**4 - t_surf_window_h(m)**4 ) &
                                    + surf_usm_h%frac(ind_pav_green,m) *                 &
                                     sigma_sb * surf_usm_h%emissivity(ind_pav_green,m) * &
                                     ( t_surf_green_h_p(m)**4 - t_surf_green_h(m)**4 )

           surf_usm_h%wghf_eb(m)   = lambda_surface *                         &
                                      ( t_surf_h_p(m) - t_wall_h(nzb_wall,m) )
           surf_usm_h%wghf_eb_green(m)  = lambda_surface_green *                        &
                                          ( t_surf_green_h_p(m) - t_green_h(nzb_wall,m) )
           surf_usm_h%wghf_eb_window(m) = lambda_surface_window *                       &
                                           ( t_surf_window_h_p(m) - t_window_h(nzb_wall,m) )

!
!--        ground/wall/roof surface heat flux
           surf_usm_h%wshf_eb(m)   = - f_shf  * ( surf_usm_h%pt1(m) - t_surf_h_p(m) / exn(k) ) *               &
                                       surf_usm_h%frac(ind_veg_wall,m)         &
                                     - f_shf_window  * ( surf_usm_h%pt1(m) - t_surf_window_h_p(m) / exn(k) ) * &
                                       surf_usm_h%frac(ind_wat_win,m)          &
                                     - f_shf_green  * ( surf_usm_h%pt1(m) - t_surf_green_h_p(m) / exn(k) ) *   &
                                       surf_usm_h%frac(ind_pav_green,m)
!           
!--        store kinematic surface heat fluxes for utilization in other processes
!--        diffusion_s, surface_layer_fluxes,...
           surf_usm_h%shf(m) = surf_usm_h%wshf_eb(m) / cp

       ENDDO
!
!--    Now, treat vertical surface elements
       DO  l = 0, 3
          DO  m = 1, surf_usm_v(l)%ns
!
!--          Get indices of respective grid point
             i = surf_usm_v(l)%i(m)
             j = surf_usm_v(l)%j(m)
             k = surf_usm_v(l)%k(m)

!
!--          TODO - how to calculate lambda_surface for horizontal (??? do you mean verical ???) surfaces
!--          (lambda_surface is set according to stratification in land surface model).
!--          Please note, for vertical surfaces no ol is defined, since 
!--          stratification is not considered in this case.
             lambda_surface = surf_usm_v(l)%lambda_surf(m)
             lambda_surface_window = surf_usm_v(l)%lambda_surf_window(m)
             lambda_surface_green = surf_usm_v(l)%lambda_surf_green(m)

#if ! defined( __nopointer )          
!
!--          calculate rho * cp coefficient at surface layer
             rho_cp  = cp * hyp(k) / ( r_d * surf_usm_v(l)%pt1(m) * exn(k) )
#endif

!--          Calculation of r_a for vertical surfaces
!--
!--          heat transfer coefficient for forced convection along vertical walls
!--          follows formulation in TUF3d model (Krayenhoff & Voogt, 2006)
!--            
!--          H = httc (Tsfc - Tair)
!--          httc = rw * (11.8 + 4.2 * Ueff) - 4.0
!--            
!--                rw: wall patch roughness relative to 1.0 for concrete
!--                Ueff: effective wind speed
!--                - 4.0 is a reduction of Rowley et al (1930) formulation based on
!--                Cole and Sturrock (1977)
!--           
!--                Ucan: Canyon wind speed
!--                wstar: convective velocity
!--                Qs: surface heat flux
!--                zH: height of the convective layer
!--                wstar = (g/Tcan*Qs*zH)**(1./3.)
                
!--          Effective velocity components must always 
!--          be defined at scalar grid point. The wall normal component is 
!--          obtained by simple linear interpolation. ( An alternative would
!--          be an logarithmic interpolation. )
!--          A roughness lenght of 0.001 is assumed for concrete (the inverse,
!--          1000 is used in the nominator for scaling)
             surf_usm_v(l)%r_a(m) = rho_cp / ( surf_usm_v(l)%z0(m) * 1000.0_wp &
                        * ( 11.8_wp + 4.2_wp *                                 &
                        SQRT( MAX( ( ( u(k,j,i) + u(k,j,i+1) ) * 0.5_wp )**2 + &
                                   ( ( v(k,j,i) + v(k,j+1,i) ) * 0.5_wp )**2 + &
                                   ( ( w(k,j,i) + w(k-1,j,i) ) * 0.5_wp )**2,  &
                              0.01_wp ) )                                      &
                           )  - 4.0_wp  ) 
!
!--          Limit aerodynamic resistance
             IF ( surf_usm_v(l)%r_a(m) < 1.0_wp )                              &
                surf_usm_v(l)%r_a(m)        = 1.0_wp          
                           
             f_shf         = rho_cp / surf_usm_v(l)%r_a(m)
             f_shf_window  = rho_cp / surf_usm_v(l)%r_a(m)
             f_shf_green   = rho_cp / surf_usm_v(l)%r_a(m)



!--          add LW up so that it can be removed in prognostic equation
             surf_usm_v(l)%rad_net_l(m) = surf_usm_v(l)%rad_sw_in(m)  -        &
                                          surf_usm_v(l)%rad_sw_out(m) +        &
                                          surf_usm_v(l)%rad_lw_in(m)  -        &
                                          surf_usm_v(l)%rad_lw_out(m)

!--           numerator of the prognostic equation
              coef_1 = surf_usm_v(l)%rad_net_l(m) +                            & ! coef +1 corresponds to -lwout included in calculation of radnet_l
             ( 3.0_wp + 1.0_wp ) * surf_usm_v(l)%emissivity(ind_veg_wall,m) *  &
                                     sigma_sb *  t_surf_v(l)%t(m) ** 4 +       &  
                                     f_shf * surf_usm_v(l)%pt1(m) +            &
                                     lambda_surface * t_wall_v(l)%t(nzb_wall,m)
              coef_window_1 = surf_usm_v(l)%rad_net_l(m) +                     & ! coef +1 corresponds to -lwout included in calculation of radnet_l
               ( 3.0_wp + 1.0_wp ) * surf_usm_v(l)%emissivity(ind_wat_win,m) * &
                                     sigma_sb * t_surf_window_v(l)%t(m) ** 4 + &  
                                     f_shf * surf_usm_v(l)%pt1(m) +            &
                                     lambda_surface_window * t_window_v(l)%t(nzb_wall,m)

              coef_green_1 = surf_usm_v(l)%rad_net_l(m) +                      & ! coef +1 corresponds to -lwout included in calculation of radnet_l
              ( 3.0_wp + 1.0_wp ) * surf_usm_v(l)%emissivity(ind_pav_green,m) *&
                                     sigma_sb * t_surf_green_v(l)%t(m) ** 4 +  &  
                                     f_shf * surf_usm_v(l)%pt1(m) +            &
                                     lambda_surface_green * t_wall_v(l)%t(nzb_wall,m)

!--           denominator of the prognostic equation
              coef_2 = 4.0_wp * surf_usm_v(l)%emissivity(ind_veg_wall,m) *     &
                                sigma_sb * t_surf_v(l)%t(m) ** 3               &
                              + lambda_surface + f_shf / exn(k)
              coef_window_2 = 4.0_wp * surf_usm_v(l)%emissivity(ind_wat_win,m) *&
                                sigma_sb * t_surf_window_v(l)%t(m) ** 3        &
                              + lambda_surface_window + f_shf / exn(k)
              coef_green_2 = 4.0_wp * surf_usm_v(l)%emissivity(ind_pav_green,m) *&
                                sigma_sb * t_surf_green_v(l)%t(m) ** 3         &
                              + lambda_surface_green + f_shf / exn(k)

!--           implicit solution when the surface layer has no heat capacity,
!--           otherwise use RK3 scheme.
              t_surf_v_p(l)%t(m) = ( coef_1 * dt_3d * tsc(2) +                 &
                             surf_usm_v(l)%c_surface(m) * t_surf_v(l)%t(m) ) / & 
                           ( surf_usm_v(l)%c_surface(m) + coef_2 * dt_3d * tsc(2) ) 
              t_surf_window_v_p(l)%t(m) = ( coef_window_1 * dt_3d * tsc(2) +                 &
                             surf_usm_v(l)%c_surface_window(m) * t_surf_window_v(l)%t(m) ) / & 
                           ( surf_usm_v(l)%c_surface_window(m) + coef_window_2 * dt_3d * tsc(2) ) 

              t_surf_green_v_p(l)%t(m) = ( coef_green_1 * dt_3d * tsc(2) +                 &
                             surf_usm_v(l)%c_surface_green(m) * t_surf_green_v(l)%t(m) ) / & 
                           ( surf_usm_v(l)%c_surface_green(m) + coef_green_2 * dt_3d * tsc(2) ) 



!--           add RK3 term
              t_surf_v_p(l)%t(m) = t_surf_v_p(l)%t(m) + dt_3d * tsc(3) *       &
                                surf_usm_v(l)%tt_surface_m(m)
              t_surf_window_v_p(l)%t(m) = t_surf_window_v_p(l)%t(m) + dt_3d * tsc(3) *     &
                                surf_usm_v(l)%tt_surface_window_m(m)
              t_surf_green_v_p(l)%t(m) = t_surf_green_v_p(l)%t(m) + dt_3d * tsc(3) *       &
                                surf_usm_v(l)%tt_surface_green_m(m)
!
!--           Store surface temperature        
              surf_usm_v(l)%pt_surface(m) =  ( surf_usm_v(l)%frac(ind_veg_wall,m) * t_surf_v_p(l)%t(m)  &
                                      + surf_usm_v(l)%frac(ind_wat_win,m) * t_surf_window_v_p(l)%t(m)  &
                                      + surf_usm_v(l)%frac(ind_pav_green,m) * t_surf_green_v_p(l)%t(m) ) &
                                      / exn(k) 
            
!--           calculate true tendency
              stend = ( t_surf_v_p(l)%t(m) - t_surf_v(l)%t(m) - dt_3d * tsc(3) *&
                        surf_usm_v(l)%tt_surface_m(m) ) / ( dt_3d  * tsc(2) )
              stend_window = ( t_surf_window_v_p(l)%t(m) - t_surf_window_v(l)%t(m) - dt_3d * tsc(3) *&
                        surf_usm_v(l)%tt_surface_window_m(m) ) / ( dt_3d  * tsc(2) )
              stend_green = ( t_surf_green_v_p(l)%t(m) - t_surf_green_v(l)%t(m) - dt_3d * tsc(3) *&
                        surf_usm_v(l)%tt_surface_green_m(m) ) / ( dt_3d  * tsc(2) )

!--           calculate t_surf tendencies for the next Runge-Kutta step
              IF ( timestep_scheme(1:5) == 'runge' )  THEN
                 IF ( intermediate_timestep_count == 1 )  THEN
                    surf_usm_v(l)%tt_surface_m(m) = stend
                    surf_usm_v(l)%tt_surface_window_m(m) = stend_window
                    surf_usm_v(l)%tt_surface_green_m(m) = stend_green
                 ELSEIF ( intermediate_timestep_count <                        &
                          intermediate_timestep_count_max )  THEN
                    surf_usm_v(l)%tt_surface_m(m) = -9.5625_wp * stend +       &
                                     5.3125_wp * surf_usm_v(l)%tt_surface_m(m)
                    surf_usm_v(l)%tt_surface_green_m(m) = -9.5625_wp * stend_green +       &
                                     5.3125_wp * surf_usm_v(l)%tt_surface_green_m(m)
                    surf_usm_v(l)%tt_surface_window_m(m) = -9.5625_wp * stend_window +       &
                                     5.3125_wp * surf_usm_v(l)%tt_surface_window_m(m)
                 ENDIF
              ENDIF

!--           in case of fast changes in the skin temperature, it is required to
!--           update the radiative fluxes in order to keep the solution stable

              IF ( ( ABS( t_surf_v_p(l)%t(m) - t_surf_v(l)%t(m) ) > 1.0_wp ) .OR. &
                   ( ABS( t_surf_green_v_p(l)%t(m) - t_surf_green_v(l)%t(m) ) > 1.0_wp ) .OR.  &
                   ( ABS( t_surf_window_v_p(l)%t(m) - t_surf_window_v(l)%t(m) ) > 1.0_wp ) ) THEN
                 force_radiation_call_l = .TRUE.
              ENDIF

!--           calculate fluxes
!--           prognostic rad_net_l is used just for output!           
              surf_usm_v(l)%rad_net_l(m) = surf_usm_v(l)%frac(ind_veg_wall,m) *                      &
                                           ( surf_usm_v(l)%rad_net_l(m) +                            &
                                           3.0_wp * sigma_sb *                                       &
                                           t_surf_v(l)%t(m)**4 - 4.0_wp * sigma_sb *                 &
                                           t_surf_v(l)%t(m)**3 * t_surf_v_p(l)%t(m) )                &
                                         + surf_usm_v(l)%frac(ind_wat_win,m) *                       &
                                           ( surf_usm_v(l)%rad_net_l(m) +                            &
                                           3.0_wp * sigma_sb *                                       &
                                           t_surf_window_v(l)%t(m)**4 - 4.0_wp * sigma_sb *          &
                                           t_surf_window_v(l)%t(m)**3 * t_surf_window_v_p(l)%t(m) )  &
                                         + surf_usm_v(l)%frac(ind_pav_green,m) *                     &
                                           ( surf_usm_v(l)%rad_net_l(m) +                            &
                                           3.0_wp * sigma_sb *                                       &
                                           t_surf_green_v(l)%t(m)**4 - 4.0_wp * sigma_sb *           &
                                           t_surf_green_v(l)%t(m)**3 * t_surf_green_v_p(l)%t(m) )

              surf_usm_v(l)%wghf_eb_window(m) = lambda_surface_window * &
                                                ( t_surf_window_v_p(l)%t(m) - t_window_v(l)%t(nzb_wall,m) )
              surf_usm_v(l)%wghf_eb(m)   = lambda_surface *                    &
                                     ( t_surf_v_p(l)%t(m) - t_wall_v(l)%t(nzb_wall,m) )
              surf_usm_v(l)%wghf_eb_green(m)  = lambda_surface_green *  &
                                                ( t_surf_green_v_p(l)%t(m) - t_green_v(l)%t(nzb_wall,m) )

!--           ground/wall/roof surface heat flux
              surf_usm_v(l)%wshf_eb(m)   =                                     &
                 - f_shf  * ( surf_usm_v(l)%pt1(m) -                           &
                 t_surf_v_p(l)%t(m) / exn(k) ) * surf_usm_v(l)%frac(ind_veg_wall,m)       &
                 - f_shf_window  * ( surf_usm_v(l)%pt1(m) -                    &
                 t_surf_window_v_p(l)%t(m) / exn(k) ) * surf_usm_v(l)%frac(ind_wat_win,m)&
                 - f_shf_green  * ( surf_usm_v(l)%pt1(m) -                     &
                 t_surf_green_v_p(l)%t(m) / exn(k) ) * surf_usm_v(l)%frac(ind_pav_green,m)

!           
!--           store kinematic surface heat fluxes for utilization in other processes
!--           diffusion_s, surface_layer_fluxes,...
              surf_usm_v(l)%shf(m) = surf_usm_v(l)%wshf_eb(m) / cp

           ENDDO

        ENDDO
!
!--     Add-up anthropogenic heat, for now only at upward-facing surfaces
        IF ( usm_anthropogenic_heat  .AND.  &
             intermediate_timestep_count == intermediate_timestep_count_max )  THEN
!--        application of the additional anthropogenic heat sources
!--        we considere the traffic for now so all heat is absorbed
!--        to the first layer, generalization would be worth. 
            
!--        calculation of actual profile coefficient
!--        ??? check time_since_reference_point ???
           dtime = mod(simulated_time + time_utc_init, 24.0_wp*3600.0_wp)
           dhour = INT(dtime/3600.0_wp)
           DO m = 1, naheatlayers
!--           Get indices of respective grid point
              i = surf_usm_h%i(m)
              j = surf_usm_h%j(m)
              k = surf_usm_h%k(m)
              IF ( k > get_topography_top_index_ji( j, i, 's' )  .AND. &
                   k <= naheatlayers )  THEN
!--              increase of pt in box i,j,k in time dt_3d 
!--              given to anthropogenic heat aheat*acoef (W*m-2)
!--              linear interpolation of coeficient
                 acoef = (REAL(dhour+1,wp)-dtime/3600.0_wp)*aheatprof(k,dhour) + &
                         (dtime/3600.0_wp-REAL(dhour,wp))*aheatprof(k,dhour+1)
                 IF ( aheat(k,j,i) > 0.0_wp )  THEN
                    pt(k,j,i) = pt(k,j,i) + aheat(k,j,i)*acoef*dt_3d/(exn(k)*rho_cp*dzu(k))
                 ENDIF
              ENDIF
           ENDDO

        ENDIF
        
!--     pt and shf are defined on nxlg:nxrg,nysg:nyng
!--     get the borders from neighbours
#if ! defined( __nopointer )
        CALL exchange_horiz( pt, nbgp )
#endif

!--     calculation of force_radiation_call:
!--     Make logical OR for all processes.
!--     Force radiation call if at least one processor forces it.
        IF ( intermediate_timestep_count == intermediate_timestep_count_max-1 )&
        THEN
#if defined( __parallel )
          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_ALLREDUCE( force_radiation_call_l, force_radiation_call,    &
                              1, MPI_LOGICAL, MPI_LOR, comm2d, ierr )
#else
          force_radiation_call = force_radiation_call_l
#endif
          force_radiation_call_l = .FALSE.
       ENDIF

    END SUBROUTINE usm_surface_energy_balance


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Swapping of timelevels for t_surf and t_wall
!> called out from subroutine swap_timelevel
!------------------------------------------------------------------------------!
    SUBROUTINE usm_swap_timelevel ( mod_count )

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN) :: mod_count
       INTEGER(iwp)             :: i
      
#if defined( __nopointer )
       t_surf_h    = t_surf_h_p
       t_wall_h    = t_wall_h_p
       t_surf_v    = t_surf_v_p
       t_wall_v    = t_wall_v_p
       t_surf_window_h    = t_surf_window_h_p
       t_window_h    = t_window_h_p
       t_surf_window_v    = t_surf_window_v_p
       t_window_v    = t_window_v_p
       t_surf_green_h    = t_surf_green_h_p
       t_surf_green_v    = t_surf_green_v_p
       t_green_h    = t_green_h_p
       t_green_v    = t_green_v_p
#else
       SELECT CASE ( mod_count )
          CASE ( 0 )
!
!--          Horizontal surfaces
             t_surf_h  => t_surf_h_1; t_surf_h_p  => t_surf_h_2
             t_wall_h     => t_wall_h_1;    t_wall_h_p     => t_wall_h_2
             t_surf_window_h  => t_surf_window_h_1; t_surf_window_h_p  => t_surf_window_h_2
             t_window_h     => t_window_h_1;    t_window_h_p     => t_window_h_2
             t_surf_green_h  => t_surf_green_h_1; t_surf_green_h_p  => t_surf_green_h_2
             t_green_h     => t_green_h_1;    t_green_h_p     => t_green_h_2
!
!--          Vertical surfaces
             t_surf_v  => t_surf_v_1; t_surf_v_p  => t_surf_v_2
             t_wall_v     => t_wall_v_1;    t_wall_v_p     => t_wall_v_2
             t_surf_window_v  => t_surf_window_v_1; t_surf_window_v_p  => t_surf_window_v_2
             t_window_v     => t_window_v_1;    t_window_v_p     => t_window_v_2
             t_surf_green_v  => t_surf_green_v_1; t_surf_green_v_p  => t_surf_green_v_2
             t_green_v     => t_green_v_1;    t_green_v_p     => t_green_v_2
          CASE ( 1 )
!
!--          Horizontal surfaces
             t_surf_h  => t_surf_h_2; t_surf_h_p  => t_surf_h_1
             t_wall_h     => t_wall_h_2;    t_wall_h_p     => t_wall_h_1
             t_surf_window_h  => t_surf_window_h_2; t_surf_window_h_p  => t_surf_window_h_1
             t_window_h     => t_window_h_2;    t_window_h_p     => t_window_h_1
             t_surf_green_h  => t_surf_green_h_2; t_surf_green_h_p  => t_surf_green_h_1
             t_green_h     => t_green_h_2;    t_green_h_p     => t_green_h_1
!
!--          Vertical surfaces
             t_surf_v  => t_surf_v_2; t_surf_v_p  => t_surf_v_1
             t_wall_v     => t_wall_v_2;    t_wall_v_p     => t_wall_v_1
             t_surf_window_v  => t_surf_window_v_2; t_surf_window_v_p  => t_surf_window_v_1
             t_window_v     => t_window_v_2;    t_window_v_p     => t_window_v_1
             t_surf_green_v  => t_surf_green_v_2; t_surf_green_v_p  => t_surf_green_v_1
             t_green_v     => t_green_v_2;    t_green_v_p     => t_green_v_1
       END SELECT
#endif
        
    END SUBROUTINE usm_swap_timelevel

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine writes t_surf and t_wall data into restart files
!------------------------------------------------------------------------------!
    SUBROUTINE usm_wrd_local

    
       IMPLICIT NONE
       
       CHARACTER(LEN=1) ::  dum     !< dummy string to create output-variable name  
       INTEGER(iwp)     ::  l       !< index surface type orientation

       CALL wrd_write_string( 'ns_h_on_file_usm' )
       WRITE ( 14 )  surf_usm_h%ns

       CALL wrd_write_string( 'ns_v_on_file_usm' )
       WRITE ( 14 )  surf_usm_v(0:3)%ns

       CALL wrd_write_string( 'usm_start_index_h' )
       WRITE ( 14 )  surf_usm_h%start_index

       CALL wrd_write_string( 'usm_end_index_h' )
       WRITE ( 14 )  surf_usm_h%end_index

       CALL wrd_write_string( 't_surf_h' )
       WRITE ( 14 )  t_surf_h

       CALL wrd_write_string( 't_surf_window_h' )
       WRITE ( 14 )  t_surf_window_h

       CALL wrd_write_string( 't_surf_green_h' )
       WRITE ( 14 )  t_surf_green_h

       DO  l = 0, 3

          CALL wrd_write_string( 'usm_start_index_v' )
          WRITE ( 14 )  surf_usm_v(l)%start_index

          CALL wrd_write_string( 'usm_end_index_v' )
          WRITE ( 14 )  surf_usm_v(l)%end_index

          WRITE( dum, '(I1)')  l          

          CALL wrd_write_string( 't_surf_v(' // dum // ')' )
          WRITE ( 14 )  t_surf_v(l)%t

          CALL wrd_write_string( 't_surf_window_v(' // dum // ')' )
          WRITE ( 14 ) t_surf_window_v(l)%t     

          CALL wrd_write_string( 't_surf_green_v(' // dum // ')' )
          WRITE ( 14 ) t_surf_green_v(l)%t    
          
       ENDDO

       CALL wrd_write_string( 'usm_start_index_h' )
       WRITE ( 14 )  surf_usm_h%start_index

       CALL wrd_write_string( 'usm_end_index_h' )
       WRITE ( 14 )  surf_usm_h%end_index

       CALL wrd_write_string( 't_wall_h' )
       WRITE ( 14 )  t_wall_h

       CALL wrd_write_string( 't_window_h' )
       WRITE ( 14 )  t_window_h

       CALL wrd_write_string( 't_green_h' )
       WRITE ( 14 )  t_green_h

       DO  l = 0, 3

          CALL wrd_write_string( 'usm_start_index_v' )
          WRITE ( 14 )  surf_usm_v(l)%start_index

          CALL wrd_write_string( 'usm_end_index_v' )
          WRITE ( 14 )  surf_usm_v(l)%end_index

          WRITE( dum, '(I1)')  l     

          CALL wrd_write_string( 't_wall_v(' // dum // ')' )
          WRITE ( 14 )  t_wall_v(l)%t

          CALL wrd_write_string( 't_window_v(' // dum // ')' )
          WRITE ( 14 )  t_window_v(l)%t

          CALL wrd_write_string( 't_green_v(' // dum // ')' )
          WRITE ( 14 )  t_green_v(l)%t
       
       ENDDO

       
    END SUBROUTINE usm_wrd_local

!
!-- Integrated stability function for heat and moisture
    FUNCTION psi_h( zeta )

           USE kinds

       IMPLICIT NONE

       REAL(wp)            :: psi_h !< Integrated similarity function result
       REAL(wp)            :: zeta  !< Stability parameter z/L
       REAL(wp)            :: x     !< dummy variable

       REAL(wp), PARAMETER :: a = 1.0_wp            !< constant
       REAL(wp), PARAMETER :: b = 0.66666666666_wp  !< constant
       REAL(wp), PARAMETER :: c = 5.0_wp            !< constant
       REAL(wp), PARAMETER :: d = 0.35_wp           !< constant
       REAL(wp), PARAMETER :: c_d_d = c / d         !< constant
       REAL(wp), PARAMETER :: bc_d_d = b * c / d    !< constant


      IF ( zeta < 0.0_wp )  THEN
         x = SQRT( 1.0_wp  - 16.0_wp * zeta )
         psi_h = 2.0_wp * LOG( (1.0_wp + x ) / 2.0_wp )
      ELSE
         psi_h = - b * ( zeta - c_d_d ) * EXP( -d * zeta ) - (1.0_wp          &
                 + 0.66666666666_wp * a * zeta )**1.5_wp - bc_d_d             &
                 + 1.0_wp
!
!--       Old version for stable conditions (only valid for z/L < 0.5)
!--       psi_h = - 5.0_wp * zeta
       ENDIF

   END FUNCTION psi_h
    
 END MODULE urban_surface_mod
