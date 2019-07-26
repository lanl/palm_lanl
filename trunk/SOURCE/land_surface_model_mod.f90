!> @file land_surface_model_mod.f90
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
! $Id: land_surface_model_mod.f90 3091 2018-06-28 16:20:35Z suehring $
! Add check for local roughness length not exceeding surface-layer height and
! limit roughness length where necessary.
! 
! 3051 2018-05-30 17:43:55Z suehring
! Bugfix in surface-element loops for pavement surfaces 
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised and added
! 
! 3026 2018-05-22 10:30:53Z schwenkel
! Changed the name specific humidity to mixing ratio, since we are computing
! mixing ratios.
! 
! 3014 2018-05-09 08:42:38Z maronga
! Bugfix: set some initial values
! Bugfix: domain bounds of local_pf corrected
! 
! 3004 2018-04-27 12:33:25Z Giersch
! Further allocation checks implemented (averaged data will be assigned to fill
! values if no allocation happened so far) 
! 
! 2968 2018-04-13 11:52:24Z suehring
! Bugfix in initialization in case of elevated model surface
! 
! 2963 2018-04-12 14:47:44Z suehring
! - In initialization of surface types, consider the case that surface_fractions
!   is not given in static input file.
! - Introduce index for vegetation/wall, pavement/green-wall and water/window 
!   surfaces, for clearer access of surface fraction, albedo, emissivity, etc. .
! 
! 2938 2018-03-27 15:52:42Z suehring
! Initialization of soil moisture and temperature via Inifor-provided data also
! in nested child domains, even if no dynamic input file is available for 
! child domain. 1D soil profiles are received from parent model.  
! 
! 2936 2018-03-27 14:49:27Z suehring
! renamed lsm_par to land_surface_parameters. Bugfix in message calls
! 
! 2921 2018-03-22 15:05:23Z Giersch
! The activation of spinup has been moved to parin
! 
! 2894 2018-03-15 09:17:58Z Giersch
! Calculations of the index range of the subdomain on file which overlaps with
! the current subdomain are already done in read_restart_data_mod,
! lsm_read/write_restart_data was renamed to lsm_r/wrd_local, USE kinds has 
! been removed in several routines, variable named found has been
! introduced for checking if restart data was found, reading of restart strings
! has been moved completely to read_restart_data_mod, lsm_rrd_local is already
! inside the overlap loop programmed in read_restart_data_mod, the marker ***
! end lsm *** is not necessary anymore, strings and their respective lengths 
! are written out and read now in case of restart runs to get rid of prescribed
! character lengths, SAVE attribute added where necessary, deallocation and 
! allocation of some arrays have been changed to take care of different restart
! files that can be opened (index i)
! 
! 2881 2018-03-13 16:24:40Z suehring
! Bugfix: wrong loop structure for soil moisture calculation
! 
! 2805 2018-02-14 17:00:09Z suehring
! Bugfix in initialization of roughness over water surfaces
! 
! 2798 2018-02-09 17:16:39Z suehring
! Minor bugfix for initialization of pt_surface
! 
! 2797 2018-02-08 13:24:35Z suehring
! Move output of ghf to general 2D output to output ghf also at urban-type 
! surfaces. 
! Move restart data of ghf_av to read/write_3d_binary, as this is not a 
! exclusively LSM variable anymore.   
! 
! 2765 2018-01-22 11:34:58Z maronga
! Major bugfix in calculation of f_shf for vertical surfaces
! 
! 2735 2018-01-11 12:01:27Z suehring
! output of r_a moved from land-surface to consider also urban-type surfaces
! 
! 2729 2018-01-09 11:22:28Z maronga
! Separated deep soil temperature from soil_temperature array
! 
! 2724 2018-01-05 12:12:38Z maronga
! Added security check for insufficient soil_temperature values
! 
! 2723 2018-01-05 09:27:03Z maronga
! Bugfix for spinups (end_time was increased twice in case of LSM + USM runs)
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2707 2017-12-18 18:34:46Z suehring
! Changes from last commit documented
! 
! 2706 2017-12-18 18:33:49Z suehring
! Bugfix, read surface temperature in case of restart runs.
!
! 2705 2017-12-18 11:26:23Z maronga
! Bugfix in binary output (wrong sequence)
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Bugfix: missing USE statement for calc_mean_profile
! do not write surface temperatures onto pt array as this might cause 
! problems with nesting (MS)
! Revised calculation of pt1 and qv1 (now done in surface_layer_fluxes). Bugfix
! in calculation of surface density (cannot be done via an surface non-air 
! temperature) (BM)
! Bugfix: g_d was NaN for non-vegetaed surface types (BM)
! Bugfix initialization of c_veg and lai
! Revise data output to enable _FillValues
! Bugfix in calcultion of r_a and rad_net_l (MS)
! Bugfix: rad_net is not updated in case of radiation_interaction and must thu
! be calculated again from the radiative fluxes
! Temporary fix for cases where no soil model is used on some PEs (BM)
! Revised input and initialization of soil and surface paramters
! pavement_depth is variable for each surface element
! radiation quantities belong to surface type now
! surface fractions initialized
! Rename lsm_last_actions into lsm_wrd_subdomain (MS)
! 
! 2608 2017-11-13 14:04:26Z schwenkel
! Calculation of magnus equation in external module (diagnostic_quantities_mod).
! Adjust calculation of vapor pressure and saturation mixing ratio that it is
! consistent with formulations in other parts of PALM.
! 
! 2575 2017-10-24 09:57:58Z maronga
! Pavement parameterization revised
! 
! 2573 2017-10-20 15:57:49Z scharf
! bugfixes in last_actions
! 
! 2548 2017-10-16 13:18:20Z suehring
! extended by cloud_droplets option
! 
! 2532 2017-10-11 16:00:46Z scharf
! bugfixes in data_output_3d
! 
! 2516 2017-10-04 11:03:04Z suehring
! Remove tabs
! 
! 2514 2017-10-04 09:52:37Z suehring
! upper bounds of cross section and 3d output changed from nx+1,ny+1 to nx,ny
! no output of ghost layer data
! 
! 2504 2017-09-27 10:36:13Z maronga
! Support roots and water under pavement. Added several pavement types.
! 
! 2476 2017-09-18 07:54:32Z maronga
! Bugfix for last commit
! 
! 2475 2017-09-18 07:42:36Z maronga
! Bugfix: setting of vegetation_pars for bare soil corrected.
! 
! 2354 2017-08-17 10:49:36Z schwenkel
! minor bugfixes
! 
! 2340 2017-08-07 17:11:13Z maronga
! Revised root_distribution tabel and implemented a pseudo-generic root fraction
! calculation
! 
! 2333 2017-08-04 09:08:26Z maronga
! minor bugfixes
! 
! 2332 2017-08-03 21:15:22Z maronga
! bugfix in pavement_pars
! 
! 2328 2017-08-03 12:34:22Z maronga
! Revised skin layer concept.
! Bugfix for runs with pavement surface and humidity
! Revised some standard values in vegetation_pars
! Added emissivity and default albedo_type as variable to tables
! Changed default surface type to vegetation
! Revised input of soil layer configuration
! 
! 2307 2017-07-07 11:32:10Z suehring
! Bugfix, variable names corrected
! 
! 2299 2017-06-29 10:14:38Z maronga
! Removed pt_p from USE statement. Adjusted call to lsm_soil_model to allow
! spinups without soil moisture prediction
! 
! 2298 2017-06-29 09:28:18Z raasch
! type of write_binary changed from CHARACTER to LOGICAL
! 
! 2296 2017-06-28 07:53:56Z maronga
! Bugfix in calculation of bare soil heat capacity.
! Bugfix in calculation of shf
! Added support for spinups
! 
! 2282 2017-06-13 11:38:46Z schwenkel
! Bugfix for check of saturation moisture
! 
! 2273 2017-06-09 12:46:06Z sward
! Error number changed
! 
! 2270 2017-06-09 12:18:47Z maronga
! Revised parameterization of heat conductivity between skin layer and soil.
! Temperature and moisture are now defined at the center of the layers.
! Renamed veg_type to vegetation_type and pave_type to pavement_type_name
! Renamed and reduced the number of look-up tables (vegetation_pars, soil_pars)
! Revised land surface model initialization
! Removed output of shf_eb and qsws_eb and removed _eb throughout code
! Removed Clapp & Hornberger parameterization
! 
! 2249 2017-06-06 13:58:01Z sward
! 
! 2248 2017-06-06 13:52:54Z sward $
! Error no changed
!
! 2246 2017-06-06 13:09:34Z sward
! Error no changed
!
! Changed soil configuration to 8 layers. The number of soil layers is now
! freely adjustable via the NAMELIST.
! 
! 2237 2017-05-31 10:34:53Z suehring
! Bugfix in write restart data
! 
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography and surface concept
!   - now, also vertical walls are possible
!   - for vertical walls, parametrization of r_a (aerodynamic resisistance) is 
!     implemented. 
!
! Add check for soil moisture, it must not exceed its saturation value.
! 
! 2149 2017-02-09 16:57:03Z scharf
! Land surface parameters II corrected for vegetation_type 18 and 19
! 
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1978 2016-07-29 12:08:31Z maronga
! Bugfix: initial values of pave_surface and water_surface were not set.
! 
! 1976 2016-07-27 13:28:04Z maronga
! Parts of the code have been reformatted. Use of radiation model output is
! generalized and simplified. Added more output quantities due to modularization
! 
! 1972 2016-07-26 07:52:02Z maronga
! Further modularization: output of cross sections and 3D data is now done in this
! module. Moreover, restart data is written and read directly within this module.
!
! 
! 1966 2016-07-18 11:54:18Z maronga
! Bugfix: calculation of m_total in soil model was not set to zero at model start
! 
! 1949 2016-06-17 07:19:16Z maronga
! Bugfix: calculation of qsws_soil_eb with precipitation = .TRUE. gave
! qsws_soil_eb = 0 due to a typo
! 
! 1856 2016-04-13 12:56:17Z maronga
! Bugfix: for water surfaces, the initial water surface temperature is set equal
! to the intital skin temperature. Moreover, the minimum value of r_a is now
! 1.0 to avoid too large fluxes at the first model time step
! 
! 1849 2016-04-08 11:33:18Z hoffmann
! prr moved to arrays_3d
!
! 1826 2016-04-07 12:01:39Z maronga
! Cleanup after modularization
! 
! 1817 2016-04-06 15:44:20Z maronga
! Added interface for lsm_init_arrays. Added subroutines for check_parameters,
! header, and parin. Renamed some subroutines.
! 
! 1788 2016-03-10 11:01:04Z maronga
! Bugfix: calculate lambda_surface based on temperature gradient between skin
! layer and soil layer instead of Obukhov length
! Changed: moved calculation of surface specific humidity to energy balance solver
! New: water surfaces are available by using a fixed sea surface temperature.
! The roughness lengths are calculated dynamically using the Charnock 
! parameterization. This involves the new roughness length for moisture z0q.
! New: modified solution of the energy balance solver and soil model for
! paved surfaces (i.e. asphalt concrete).
! Syntax layout improved.
! Changed: parameter dewfall removed.
! 
! 1783 2016-03-06 18:36:17Z raasch
! netcdf variables moved to netcdf module
!
! 1757 2016-02-22 15:49:32Z maronga
! Bugfix: set tm_soil_m to zero after allocation. Added parameter 
! unscheduled_radiation_calls to control calls of the radiation model based on
! the skin temperature change during one time step (preliminary version). Set 
! qsws_soil_eb to zero at model start (previously set to qsws_eb). Removed MAX
! function as it cannot be vectorized.
! 
! 1709 2015-11-04 14:47:01Z maronga
! Renamed pt_1 and qv_1 to pt1 and qv1.
! Bugfix: set initial values for t_surface_p in case of restart runs
! Bugfix: zero resistance caused crash when using radiation_scheme = 'clear-sky'
! Bugfix: calculation of rad_net when using radiation_scheme = 'clear-sky'
! Added todo action
!
! 1697 2015-10-28 17:14:10Z raasch
! bugfix: misplaced cpp-directive
!
! 1695 2015-10-27 10:03:11Z maronga
! Bugfix: REAL constants provided with KIND-attribute in call of 
! Replaced rif with ol
! 
! 1691 2015-10-26 16:17:44Z maronga
! Added skip_time_do_lsm to allow for spin-ups without LSM. Various bugfixes:
! Soil temperatures are now defined at the edges of the layers, calculation of
! shb_eb corrected, prognostic equation for skin temperature corrected. Surface
! fluxes are now directly transfered to atmosphere
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1590 2015-05-08 13:56:27Z maronga
! Bugfix: definition of character strings requires same length for all elements
! 
! 1585 2015-04-30 07:05:52Z maronga
! Modifications for RRTMG. Changed tables to PARAMETER type. 
! 
! 1571 2015-03-12 16:12:49Z maronga
! Removed upper-case variable names. Corrected distribution of precipitation to 
! the liquid water reservoir and the bare soil fractions.
! 
! 1555 2015-03-04 17:44:27Z maronga
! Added output of r_a and r_s
! 
! 1553 2015-03-03 17:33:54Z maronga
! Improved better treatment of roughness lengths. Added default soil temperature
! profile
! 
! 1551 2015-03-03 14:18:16Z maronga
! Flux calculation is now done in prandtl_fluxes. Added support for data output.
! Vertical indices have been replaced. Restart runs are now possible. Some
! variables have beem renamed. Bugfix in the prognostic equation for the surface
! temperature. Introduced z0_eb and z0h_eb, which overwrite the setting of
! roughness_length and z0_factor. Added Clapp & Hornberger parametrization for
! the hydraulic conductivity. Bugfix for root fraction and extraction
! calculation
! 
! intrinsic function MAX and MIN
!
! 1500 2014-12-03 17:42:41Z maronga
! Corrected calculation of aerodynamic resistance (r_a).
! Precipitation is now added to liquid water reservoir using LE_liq.
! Added support for dry runs.
! 
! 1496 2014-12-02 17:25:50Z maronga
! Initial revision
! 
!
! Description:
! ------------
!> Land surface model, consisting of a solver for the energy balance at the
!> surface and a multi layer soil scheme. The scheme is similar to the TESSEL
!> scheme implemented in the ECMWF IFS model, with modifications according to
!> H-TESSEL. The implementation is based on the formulation implemented in the 
!> DALES and UCLA-LES models.
!>
!> @todo Extensive verification energy-balance solver for vertical surfaces, 
!>       e.g. parametrization of r_a
!> @todo Revise single land-surface processes for vertical surfaces, e.g. 
!>       treatment of humidity, etc. 
!> @todo Consider partial absorption of the net shortwave radiation by the 
!>       skin layer.
!> @todo Improve surface water parameterization
!> @todo Invert indices (running from -3 to 0. Currently: nzb_soil=0, 
!>       nzt_soil=3)).
!> @todo Implement surface runoff model (required when performing long-term LES 
!>       with considerable precipitation.
!> @todo Revise calculation of f2 when wilting point is non-constant in the
!>       soil
!> @todo Allow for zero soil moisture (currently, it is set to wilting point)
!> @note No time step criterion is required as long as the soil layers do not
!>       become too thin.
!> @todo Attention, pavement_subpars_1/2 are hardcoded to 8 levels, in case 
!>       more levels are used this may cause an potential bug
!> @todo Routine calc_q_surface required?
!------------------------------------------------------------------------------!
 MODULE land_surface_model_mod
 
    USE arrays_3d,                                                             &
        ONLY:  hyp, pt, prr, q, q_p, ql, vpt, u, v, w

    USE calc_mean_profile_mod,                                                 &
        ONLY:  calc_mean_profile

    USE cloud_parameters,                                                      &
        ONLY:  cp, hyrho, l_d_cp, l_d_r, l_v, pt_d_t, rho_l, r_d, r_v

    USE control_parameters,                                                    &
        ONLY:  cloud_droplets, cloud_physics, coupling_start_time, dt_3d,      &
               end_time, humidity, intermediate_timestep_count,                &
               initializing_actions, intermediate_timestep_count_max,          &
               land_surface, max_masks, precipitation, pt_surface,             &
               rho_surface, spinup, spinup_pt_mean, spinup_time,               &
               surface_pressure, timestep_scheme, tsc,                         &
               time_since_reference_point

    USE indices,                                                               &
        ONLY:  nbgp, nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb

    USE netcdf_data_input_mod,                                                 &
        ONLY :  building_type_f, init_3d, input_pids_static,                   &
                netcdf_data_input_interpolate,                                 &
                pavement_pars_f, pavement_subsurface_pars_f, pavement_type_f,  &
                root_area_density_lsm_f, soil_pars_f, soil_type_f,             &
                surface_fraction_f, vegetation_pars_f, vegetation_type_f,      &
                water_pars_f, water_type_f

    USE kinds

    USE pegrid

    USE radiation_model_mod,                                                   &
        ONLY:  albedo, albedo_type, emissivity, force_radiation_call,          &
               radiation, radiation_scheme, unscheduled_radiation_calls
        
    USE statistics,                                                            &
        ONLY:  hom, statistic_regions

    USE surface_mod,                                                           &
        ONLY :  ind_pav_green, ind_veg_wall, ind_wat_win, surf_lsm_h,          &
                surf_lsm_v, surf_type, surface_restore_elements

    IMPLICIT NONE

    TYPE surf_type_lsm
       REAL(wp), DIMENSION(:),   ALLOCATABLE ::  var_1d !< 1D prognostic variable
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  var_2d !< 2D prognostic variable
    END TYPE surf_type_lsm

!
!-- LSM model constants

    REAL(wp), PARAMETER  ::                    &
              b_ch               = 6.04_wp,    & ! Clapp & Hornberger exponent
              lambda_h_dry       = 0.19_wp,    & ! heat conductivity for dry soil (W/m/K)  
              lambda_h_sm        = 3.44_wp,    & ! heat conductivity of the soil matrix (W/m/K)
              lambda_h_water     = 0.57_wp,    & ! heat conductivity of water (W/m/K)
              psi_sat            = -0.388_wp,  & ! soil matrix potential at saturation
              rho_c_soil         = 2.19E6_wp,  & ! volumetric heat capacity of soil (J/m3/K)
              rho_c_water        = 4.20E6_wp,  & ! volumetric heat capacity of water (J/m3/K)
              m_max_depth        = 0.0002_wp     ! Maximum capacity of the water reservoir (m)


    REAL(wp), DIMENSION(0:7), PARAMETER  :: dz_soil_default =                  & ! default soil layer configuration
                                            (/ 0.01_wp, 0.02_wp, 0.04_wp,      &
                                               0.06_wp, 0.14_wp, 0.26_wp,      &
                                               0.54_wp, 1.86_wp/)

    REAL(wp), DIMENSION(0:3), PARAMETER  :: dz_soil_ref =                      & ! reference four layer soil configuration used for estimating the root fractions
                                            (/ 0.07_wp, 0.21_wp, 0.72_wp,      &
                                               1.89_wp /)

    REAL(wp), DIMENSION(0:3), PARAMETER  :: zs_ref =                           & ! reference four layer soil configuration used for estimating the root fractions
                                            (/ 0.07_wp, 0.28_wp, 1.0_wp,       &
                                               2.89_wp /)


!
!-- LSM variables
    CHARACTER(10) :: surface_type = 'netcdf'      !< general classification. Allowed are:
                                                  !< 'vegetation', 'pavement', ('building'), 
                                                  !< 'water', and 'netcdf'



    INTEGER(iwp) :: nzb_soil = 0,             & !< bottom of the soil model (Earth's surface)
                    nzt_soil = 7,             & !< top of the soil model
                    nzt_pavement = 0,         & !< top of the pavement within the soil
                    nzs = 8,                  & !< number of soil layers
                    pavement_depth_level = 0, & !< default NAMELIST nzt_pavement
                    pavement_type = 1,        & !< default NAMELIST pavement_type                 
                    soil_type = 3,            & !< default NAMELIST soil_type
                    vegetation_type = 2,      & !< default NAMELIST vegetation_type
                    water_type = 1              !< default NAMELISt water_type
                    
   
        
    LOGICAL :: conserve_water_content = .TRUE.,  & !< open or closed bottom surface for the soil model
               constant_roughness = .FALSE.,     & !< use fixed/dynamic roughness lengths for water surfaces
               force_radiation_call_l = .FALSE., & !< flag to force calling of radiation routine
               aero_resist_kray = .TRUE.           !< flag to control parametrization of aerodynamic resistance at vertical surface elements

!   value 9999999.9_wp -> generic available or user-defined value must be set
!   otherwise -> no generic variable and user setting is optional
    REAL(wp) :: alpha_vangenuchten = 9999999.9_wp,      & !< NAMELIST alpha_vg
                canopy_resistance_coefficient = 9999999.9_wp, & !< NAMELIST g_d
                c_surface = 9999999.9_wp,               & !< Surface (skin) heat capacity (J/m2/K)
                deep_soil_temperature =  9999999.9_wp,  & !< Deep soil temperature (bottom boundary condition)
                drho_l_lv,                              & !< (rho_l * l_v)**-1
                exn,                                    & !< value of the Exner function
                e_s = 0.0_wp,                           & !< saturation water vapour pressure
                field_capacity = 9999999.9_wp,          & !< NAMELIST m_fc
                f_shortwave_incoming = 9999999.9_wp,    & !< NAMELIST f_sw_in
                hydraulic_conductivity = 9999999.9_wp,  & !< NAMELIST gamma_w_sat
                ke = 0.0_wp,                            & !< Kersten number
                lambda_h_sat = 0.0_wp,                  & !< heat conductivity for saturated soil (W/m/K)
                lambda_surface_stable = 9999999.9_wp,   & !< NAMELIST lambda_surface_s (W/m2/K)
                lambda_surface_unstable = 9999999.9_wp, & !< NAMELIST lambda_surface_u (W/m2/K)
                leaf_area_index = 9999999.9_wp,         & !< NAMELIST lai
                l_vangenuchten = 9999999.9_wp,          & !< NAMELIST l_vg
                min_canopy_resistance = 9999999.9_wp,   & !< NAMELIST r_canopy_min
                min_soil_resistance = 50.0_wp,          & !< NAMELIST r_soil_min
                m_total = 0.0_wp,                       & !< weighted total water content of the soil (m3/m3)
                n_vangenuchten = 9999999.9_wp,          & !< NAMELIST n_vg
                pavement_heat_capacity = 9999999.9_wp,  & !< volumetric heat capacity of pavement (e.g. roads) (J/m3/K)
                pavement_heat_conduct  = 9999999.9_wp,  & !< heat conductivity for pavements (e.g. roads) (W/m/K)
                q_s = 0.0_wp,                           & !< saturation water vapor mixing ratio
                residual_moisture = 9999999.9_wp,       & !< NAMELIST m_res
                rho_cp,                                 & !< rho_surface * cp
                rho_lv,                                 & !< rho_ocean * l_v
                rd_d_rv,                                & !< r_d / r_v
                saturation_moisture = 9999999.9_wp,     & !< NAMELIST m_sat
                skip_time_do_lsm = 0.0_wp,              & !< LSM is not called before this time
                vegetation_coverage = 9999999.9_wp,     & !< NAMELIST c_veg
                water_temperature = 9999999.9_wp,       & !< water temperature
                wilting_point = 9999999.9_wp,           & !< NAMELIST m_wilt
                z0_vegetation  = 9999999.9_wp,          & !< NAMELIST z0 (lsm_par)
                z0h_vegetation = 9999999.9_wp,          & !< NAMELIST z0h (lsm_par)
                z0q_vegetation = 9999999.9_wp,          & !< NAMELIST z0q (lsm_par)
                z0_pavement    = 9999999.9_wp,          & !< NAMELIST z0 (lsm_par)
                z0h_pavement   = 9999999.9_wp,          & !< NAMELIST z0h (lsm_par)
                z0q_pavement   = 9999999.9_wp,          & !< NAMELIST z0q (lsm_par)
                z0_water       = 9999999.9_wp,          & !< NAMELIST z0 (lsm_par)
                z0h_water      = 9999999.9_wp,          & !< NAMELIST z0h (lsm_par)
                z0q_water      = 9999999.9_wp             !< NAMELIST z0q (lsm_par)  
                
                
    REAL(wp), DIMENSION(:), ALLOCATABLE  :: ddz_soil_center, & !< 1/dz_soil_center
                                            ddz_soil,        & !< 1/dz_soil
                                            dz_soil_center,  & !< soil grid spacing (center-center)
                                            zs,              & !< depth of the temperature/moisute levels
                                            root_extr          !< root extraction


                                            
    REAL(wp), DIMENSION(0:20)  ::  root_fraction = 9999999.9_wp,     & !< (NAMELIST) distribution of root surface area to the individual soil layers
                                   soil_moisture = 0.0_wp,           & !< NAMELIST soil moisture content (m3/m3)
                                   soil_temperature = 9999999.9_wp,  & !< NAMELIST soil temperature (K) +1
                                   dz_soil  = 9999999.9_wp,          & !< (NAMELIST) soil layer depths (spacing)
                                   zs_layer = 9999999.9_wp         !< soil layer depths (edge)
                                 
#if defined( __nopointer )
    TYPE(surf_type_lsm), TARGET  ::  t_soil_h,    & !< Soil temperature (K), horizontal surface elements
                                     t_soil_h_p,  & !< Prog. soil temperature (K), horizontal surface elements
                                     m_soil_h,    & !< Soil moisture (m3/m3), horizontal surface elements
                                     m_soil_h_p     !< Prog. soil moisture (m3/m3), horizontal surface elements

    TYPE(surf_type_lsm), DIMENSION(0:3), TARGET  ::  &
                                     t_soil_v,       & !< Soil temperature (K), vertical surface elements
                                     t_soil_v_p,     & !< Prog. soil temperature (K), vertical surface elements
                                     m_soil_v,       & !< Soil moisture (m3/m3), vertical surface elements
                                     m_soil_v_p        !< Prog. soil moisture (m3/m3), vertical surface elements

#else
    TYPE(surf_type_lsm), POINTER ::  t_soil_h,    & !< Soil temperature (K), horizontal surface elements
                                     t_soil_h_p,  & !< Prog. soil temperature (K), horizontal surface elements
                                     m_soil_h,    & !< Soil moisture (m3/m3), horizontal surface elements
                                     m_soil_h_p     !< Prog. soil moisture (m3/m3), horizontal surface elements 

    TYPE(surf_type_lsm), TARGET  ::  t_soil_h_1,  & !<
                                     t_soil_h_2,  & !<
                                     m_soil_h_1,  & !<
                                     m_soil_h_2     !<

    TYPE(surf_type_lsm), DIMENSION(:), POINTER :: &
                                     t_soil_v,    & !< Soil temperature (K), vertical surface elements
                                     t_soil_v_p,  & !< Prog. soil temperature (K), vertical surface elements
                                     m_soil_v,    & !< Soil moisture (m3/m3), vertical surface elements
                                     m_soil_v_p     !< Prog. soil moisture (m3/m3), vertical surface elements   

    TYPE(surf_type_lsm), DIMENSION(0:3), TARGET ::&
                                     t_soil_v_1,  & !<
                                     t_soil_v_2,  & !<
                                     m_soil_v_1,  & !<
                                     m_soil_v_2     !<
#endif    

#if defined( __nopointer )
    TYPE(surf_type_lsm), TARGET   ::  t_surface_h,    & !< surface temperature (K), horizontal surface elements 
                                      t_surface_h_p,  & !< progn. surface temperature (K), horizontal surface elements 
                                      m_liq_h,        & !< liquid water reservoir (m), horizontal surface elements 
                                      m_liq_h_p         !< progn. liquid water reservoir (m), horizontal surface elements 

    TYPE(surf_type_lsm), DIMENSION(0:3), TARGET   ::  &
                                      t_surface_v,    & !< surface temperature (K), vertical surface elements 
                                      t_surface_v_p,  & !< progn. surface temperature (K), vertical surface elements 
                                      m_liq_v,        & !< liquid water reservoir (m), vertical surface elements 
                                      m_liq_v_p         !< progn. liquid water reservoir (m), vertical surface elements 
#else
    TYPE(surf_type_lsm), POINTER  ::  t_surface_h,    & !< surface temperature (K), horizontal surface elements 
                                      t_surface_h_p,  & !< progn. surface temperature (K), horizontal surface elements 
                                      m_liq_h,        & !< liquid water reservoir (m), horizontal surface elements 
                                      m_liq_h_p         !< progn. liquid water reservoir (m), horizontal surface elements 

    TYPE(surf_type_lsm), TARGET   ::  t_surface_h_1,  & !<
                                      t_surface_h_2,  & !<
                                      m_liq_h_1,      & !<
                                      m_liq_h_2         !<

    TYPE(surf_type_lsm), DIMENSION(:), POINTER  ::    &
                                      t_surface_v,    & !< surface temperature (K), vertical surface elements 
                                      t_surface_v_p,  & !< progn. surface temperature (K), vertical surface elements 
                                      m_liq_v,        & !< liquid water reservoir (m), vertical surface elements 
                                      m_liq_v_p         !< progn. liquid water reservoir (m), vertical surface elements 

    TYPE(surf_type_lsm), DIMENSION(0:3), TARGET   ::  &
                                      t_surface_v_1,  & !<
                                      t_surface_v_2,  & !<
                                      m_liq_v_1,      & !<
                                      m_liq_v_2         !<
#endif

#if defined( __nopointer )
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET :: m_liq_av
#else
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET :: m_liq_av
#endif

#if defined( __nopointer )
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  t_soil_av, & !< Average of t_soil
                                                        m_soil_av    !< Average of m_soil
#else
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  t_soil_av, & !< Average of t_soil
                                                        m_soil_av    !< Average of m_soil
#endif

    TYPE(surf_type_lsm), TARGET ::  tm_liq_h_m      !< liquid water reservoir tendency (m), horizontal surface elements 
    TYPE(surf_type_lsm), TARGET ::  tt_surface_h_m  !< surface temperature tendency (K), horizontal surface elements 
    TYPE(surf_type_lsm), TARGET ::  tt_soil_h_m     !< t_soil storage array, horizontal surface elements 
    TYPE(surf_type_lsm), TARGET ::  tm_soil_h_m     !< m_soil storage array, horizontal surface elements 

    TYPE(surf_type_lsm), DIMENSION(0:3), TARGET ::  tm_liq_v_m      !< liquid water reservoir tendency (m), vertical surface elements 
    TYPE(surf_type_lsm), DIMENSION(0:3), TARGET ::  tt_surface_v_m  !< surface temperature tendency (K), vertical surface elements 
    TYPE(surf_type_lsm), DIMENSION(0:3), TARGET ::  tt_soil_v_m     !< t_soil storage array, vertical surface elements 
    TYPE(surf_type_lsm), DIMENSION(0:3), TARGET ::  tm_soil_v_m     !< m_soil storage array, vertical surface elements 

!
!-- Energy balance variables                
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: &
              c_liq_av,         & !< average of c_liq
              c_soil_av,        & !< average of c_soil
              c_veg_av,         & !< average of c_veg
              lai_av,           & !< average of lai
              qsws_liq_av,      & !< average of qsws_liq
              qsws_soil_av,     & !< average of qsws_soil
              qsws_veg_av,      & !< average of qsws_veg
              r_s_av              !< average of r_s
                   

!
!-- Predefined Land surface classes (vegetation_type)
    CHARACTER(26), DIMENSION(0:18), PARAMETER :: vegetation_type_name = (/ &
                                   'user defined              ',           & !  0 
                                   'bare soil                 ',           & !  1                           
                                   'crops, mixed farming      ',           & !  2
                                   'short grass               ',           & !  3
                                   'evergreen needleleaf trees',           & !  4
                                   'deciduous needleleaf trees',           & !  5
                                   'evergreen broadleaf trees ',           & !  6
                                   'deciduous broadleaf trees ',           & !  7
                                   'tall grass                ',           & !  8
                                   'desert                    ',           & !  9
                                   'tundra                    ',           & ! 10
                                   'irrigated crops           ',           & ! 11
                                   'semidesert                ',           & ! 12
                                   'ice caps and glaciers     ',           & ! 13
                                   'bogs and marshes          ',           & ! 14
                                   'evergreen shrubs          ',           & ! 15
                                   'deciduous shrubs          ',           & ! 16
                                   'mixed forest/woodland     ',           & ! 17
                                   'interrupted forest        '            & ! 18
                                                                 /)

!
!-- Soil model classes (soil_type)
    CHARACTER(12), DIMENSION(0:6), PARAMETER :: soil_type_name = (/ &
                                   'user defined',                  & ! 0 
                                   'coarse      ',                  & ! 1
                                   'medium      ',                  & ! 2
                                   'medium-fine ',                  & ! 3
                                   'fine        ',                  & ! 4
                                   'very fine   ',                  & ! 5
                                   'organic     '                   & ! 6
                                                                 /)

!
!-- Pavement classes
    CHARACTER(29), DIMENSION(0:15), PARAMETER :: pavement_type_name = (/ &
                                   'user defined                 ', & ! 0 
                                   'asphalt/concrete mix         ', & ! 1
                                   'asphalt (asphalt concrete)   ', & ! 2
                                   'concrete (Portland concrete) ', & ! 3
                                   'sett                         ', & ! 4
                                   'paving stones                ', & ! 5
                                   'cobblestone                  ', & ! 6
                                   'metal                        ', & ! 7
                                   'wood                         ', & ! 8
                                   'gravel                       ', & ! 9
                                   'fine gravel                  ', & ! 10
                                   'pebblestone                  ', & ! 11
                                   'woodchips                    ', & ! 12
                                   'tartan (sports)              ', & ! 13
                                   'artifical turf (sports)      ', & ! 14
                                   'clay (sports)                '  & ! 15
                                                                 /)                                                              
                                                                 
!
!-- Water classes
    CHARACTER(12), DIMENSION(0:5), PARAMETER :: water_type_name = (/ &
                                   'user defined',                   & ! 0 
                                   'lake        ',                   & ! 1
                                   'river       ',                   & ! 2
                                   'ocean       ',                   & ! 3
                                   'pond        ',                   & ! 4
                                   'fountain    '                    & ! 5
                                                                  /)                                                                                  
                    
!
!-- Land surface parameters according to the respective classes (vegetation_type)
    INTEGER(iwp) ::  ind_v_rc_min = 0    !< index for r_canopy_min in vegetation_pars
    INTEGER(iwp) ::  ind_v_rc_lai = 1    !< index for LAI in vegetation_pars
    INTEGER(iwp) ::  ind_v_c_veg   = 2   !< index for c_veg in vegetation_pars
    INTEGER(iwp) ::  ind_v_gd  = 3       !< index for g_d in vegetation_pars
    INTEGER(iwp) ::  ind_v_z0 = 4        !< index for z0 in vegetation_pars
    INTEGER(iwp) ::  ind_v_z0qh = 5      !< index for z0h / z0q in vegetation_pars
    INTEGER(iwp) ::  ind_v_lambda_s = 6  !< index for lambda_s_s in vegetation_pars
    INTEGER(iwp) ::  ind_v_lambda_u = 7  !< index for lambda_s_u in vegetation_pars
    INTEGER(iwp) ::  ind_v_f_sw_in = 8   !< index for f_sw_in in vegetation_pars
    INTEGER(iwp) ::  ind_v_c_surf = 9    !< index for c_surface in vegetation_pars
    INTEGER(iwp) ::  ind_v_at = 10       !< index for albedo_type in vegetation_pars
    INTEGER(iwp) ::  ind_v_emis = 11     !< index for emissivity in vegetation_pars

    INTEGER(iwp) ::  ind_w_temp     = 0    !< index for temperature in water_pars
    INTEGER(iwp) ::  ind_w_z0       = 1    !< index for z0 in water_pars
    INTEGER(iwp) ::  ind_w_z0h      = 2    !< index for z0h in water_pars
    INTEGER(iwp) ::  ind_w_lambda_s = 3    !< index for lambda_s_s in water_pars
    INTEGER(iwp) ::  ind_w_lambda_u = 4    !< index for lambda_s_u in water_pars
    INTEGER(iwp) ::  ind_w_at       = 5    !< index for albedo type in water_pars
    INTEGER(iwp) ::  ind_w_emis     = 6    !< index for emissivity in water_pars

    INTEGER(iwp) ::  ind_p_z0       = 0    !< index for z0 in pavement_pars
    INTEGER(iwp) ::  ind_p_z0h      = 1    !< index for z0h in pavement_pars
    INTEGER(iwp) ::  ind_p_at       = 2    !< index for albedo type in pavement_pars
    INTEGER(iwp) ::  ind_p_emis     = 3    !< index for emissivity in pavement_pars
    INTEGER(iwp) ::  ind_p_lambda_h = 0    !< index for lambda_h in pavement_subsurface_pars
    INTEGER(iwp) ::  ind_p_rho_c    = 1    !< index for rho_c in pavement_pars
!
!-- Land surface parameters
!-- r_canopy_min,     lai,   c_veg,     g_d         z0,         z0h, lambda_s_s, lambda_s_u, f_sw_in,  c_surface, albedo_type, emissivity
    REAL(wp), DIMENSION(0:11,1:18), PARAMETER :: vegetation_pars = RESHAPE( (/ &
          0.0_wp, 0.00_wp, 0.00_wp, 0.00_wp,  0.005_wp,   0.5E-4_wp,     0.0_wp,    0.0_wp, 0.00_wp, 0.00_wp, 17.0_wp, 0.94_wp, & !  1
        180.0_wp, 3.00_wp, 1.00_wp, 0.00_wp,   0.10_wp,    0.001_wp,    10.0_wp,   10.0_wp, 0.05_wp, 0.00_wp,  2.0_wp, 0.95_wp, & !  2
        110.0_wp, 2.00_wp, 1.00_wp, 0.00_wp,   0.03_wp,   0.3E-4_wp,    10.0_wp,   10.0_wp, 0.05_wp, 0.00_wp,  2.0_wp, 0.95_wp, & !  3
        500.0_wp, 5.00_wp, 1.00_wp, 0.03_wp,   2.00_wp,     2.00_wp,    20.0_wp,   15.0_wp, 0.03_wp, 0.00_wp,  5.0_wp, 0.97_wp, & !  4
        500.0_wp, 5.00_wp, 1.00_wp, 0.03_wp,   2.00_wp,     2.00_wp,    20.0_wp,   15.0_wp, 0.03_wp, 0.00_wp,  6.0_wp, 0.97_wp, & !  5
        175.0_wp, 5.00_wp, 1.00_wp, 0.03_wp,   2.00_wp,     2.00_wp,    20.0_wp,   15.0_wp, 0.03_wp, 0.00_wp,  8.0_wp, 0.97_wp, & !  6
        240.0_wp, 6.00_wp, 0.99_wp, 0.13_wp,   2.00_wp,     2.00_wp,    20.0_wp,   15.0_wp, 0.03_wp, 0.00_wp,  9.0_wp, 0.97_wp, & !  7
        100.0_wp, 2.00_wp, 0.70_wp, 0.00_wp,   0.47_wp,  0.47E-2_wp,    10.0_wp,   10.0_wp, 0.05_wp, 0.00_wp,  8.0_wp, 0.97_wp, & !  8
        250.0_wp, 0.05_wp, 0.00_wp, 0.00_wp,  0.013_wp, 0.013E-2_wp,    15.0_wp,   15.0_wp, 0.00_wp, 0.00_wp,  3.0_wp, 0.94_wp, & !  9
         80.0_wp, 1.00_wp, 0.50_wp, 0.00_wp,  0.034_wp, 0.034E-2_wp,    10.0_wp,   10.0_wp, 0.05_wp, 0.00_wp, 11.0_wp, 0.97_wp, & ! 10
        180.0_wp, 3.00_wp, 1.00_wp, 0.00_wp,    0.5_wp,  0.50E-2_wp,    10.0_wp,   10.0_wp, 0.05_wp, 0.00_wp, 13.0_wp, 0.97_wp, & ! 11
        150.0_wp, 0.50_wp, 0.10_wp, 0.00_wp,   0.17_wp,  0.17E-2_wp,    10.0_wp,   10.0_wp, 0.05_wp, 0.00_wp,  2.0_wp, 0.97_wp, & ! 12
          0.0_wp, 0.00_wp, 0.00_wp, 0.00_wp, 1.3E-3_wp,   1.3E-4_wp,    58.0_wp,   58.0_wp, 0.00_wp, 0.00_wp, 11.0_wp, 0.97_wp, & ! 13
        240.0_wp, 4.00_wp, 0.60_wp, 0.00_wp,   0.83_wp,  0.83E-2_wp,    10.0_wp,   10.0_wp, 0.05_wp, 0.00_wp,  4.0_wp, 0.97_wp, & ! 14
        225.0_wp, 3.00_wp, 0.50_wp, 0.00_wp,   0.10_wp,  0.10E-2_wp,    10.0_wp,   10.0_wp, 0.05_wp, 0.00_wp,  4.0_wp, 0.97_wp, & ! 15
        225.0_wp, 1.50_wp, 0.50_wp, 0.00_wp,   0.25_wp,  0.25E-2_wp,    10.0_wp,   10.0_wp, 0.05_wp, 0.00_wp,  4.0_wp, 0.97_wp, & ! 16
        250.0_wp, 5.00_wp, 1.00_wp, 0.03_wp,   2.00_wp,     2.00_wp,    20.0_wp,   15.0_wp, 0.03_wp, 0.00_wp,  7.0_wp, 0.97_wp, & ! 17
        175.0_wp, 2.50_wp, 1.00_wp, 0.03_wp,   1.10_wp,     1.10_wp,    20.0_wp,   15.0_wp, 0.03_wp, 0.00_wp,  8.0_wp, 0.97_wp  & ! 18
                                                               /), (/ 12, 18 /) )

                                    
!
!-- Root distribution for default soil layer configuration (sum = 1)
!--                                level 1 - level 4 according to zs_ref
    REAL(wp), DIMENSION(0:3,1:18), PARAMETER :: root_distribution = RESHAPE( (/ &
                                 1.00_wp, 0.00_wp, 0.00_wp, 0.00_wp,            & !  1
                                 0.24_wp, 0.41_wp, 0.31_wp, 0.04_wp,            & !  2
                                 0.35_wp, 0.38_wp, 0.23_wp, 0.04_wp,            & !  3
                                 0.26_wp, 0.39_wp, 0.29_wp, 0.06_wp,            & !  4
                                 0.26_wp, 0.38_wp, 0.29_wp, 0.07_wp,            & !  5
                                 0.24_wp, 0.38_wp, 0.31_wp, 0.07_wp,            & !  6
                                 0.25_wp, 0.34_wp, 0.27_wp, 0.14_wp,            & !  7
                                 0.27_wp, 0.27_wp, 0.27_wp, 0.09_wp,            & !  8
                                 1.00_wp, 0.00_wp, 0.00_wp, 0.00_wp,            & !  9
                                 0.47_wp, 0.45_wp, 0.08_wp, 0.00_wp,            & ! 10
                                 0.24_wp, 0.41_wp, 0.31_wp, 0.04_wp,            & ! 11
                                 0.17_wp, 0.31_wp, 0.33_wp, 0.19_wp,            & ! 12
                                 0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp,            & ! 13
                                 0.25_wp, 0.34_wp, 0.27_wp, 0.11_wp,            & ! 14
                                 0.23_wp, 0.36_wp, 0.30_wp, 0.11_wp,            & ! 15 
                                 0.23_wp, 0.36_wp, 0.30_wp, 0.11_wp,            & ! 16 
                                 0.19_wp, 0.35_wp, 0.36_wp, 0.10_wp,            & ! 17
                                 0.19_wp, 0.35_wp, 0.36_wp, 0.10_wp             & ! 18
                                 /), (/ 4, 18 /) )

!
!-- Soil parameters according to the following porosity classes (soil_type)

!
!-- Soil parameters  alpha_vg,      l_vg,    n_vg, gamma_w_sat,    m_sat,     m_fc,   m_wilt,    m_res 
    REAL(wp), DIMENSION(0:7,1:6), PARAMETER :: soil_pars = RESHAPE( (/     &
                      3.83_wp,  1.250_wp, 1.38_wp,  6.94E-6_wp, 0.403_wp, 0.244_wp, 0.059_wp, 0.025_wp,& ! 1
                      3.14_wp, -2.342_wp, 1.28_wp,  1.16E-6_wp, 0.439_wp, 0.347_wp, 0.151_wp, 0.010_wp,& ! 2
                      0.83_wp, -0.588_wp, 1.25_wp,  0.26E-6_wp, 0.430_wp, 0.383_wp, 0.133_wp, 0.010_wp,& ! 3
                      3.67_wp, -1.977_wp, 1.10_wp,  2.87E-6_wp, 0.520_wp, 0.448_wp, 0.279_wp, 0.010_wp,& ! 4
                      2.65_wp,  2.500_wp, 1.10_wp,  1.74E-6_wp, 0.614_wp, 0.541_wp, 0.335_wp, 0.010_wp,& ! 5
                      1.30_wp,  0.400_wp, 1.20_wp,  0.93E-6_wp, 0.766_wp, 0.663_wp, 0.267_wp, 0.010_wp & ! 6
                                                                     /), (/ 8, 6 /) )


!
!-- TO BE FILLED
!-- Pavement parameters      z0,       z0h, albedo_type, emissivity  
    REAL(wp), DIMENSION(0:3,1:15), PARAMETER :: pavement_pars = RESHAPE( (/ &
                      1.0E-4_wp, 1.0E-5_wp,     18.0_wp,    0.97_wp,  & !  1
                      1.0E-4_wp, 1.0E-5_wp,     19.0_wp,    0.94_wp,  & !  2
                      1.0E-4_wp, 1.0E-5_wp,     20.0_wp,    0.98_wp,  & !  3                                  
                      1.0E-4_wp, 1.0E-5_wp,     21.0_wp,    0.93_wp,  & !  4
                      1.0E-4_wp, 1.0E-5_wp,     22.0_wp,    0.97_wp,  & !  5
                      1.0E-4_wp, 1.0E-5_wp,     23.0_wp,    0.97_wp,  & !  6
                      1.0E-4_wp, 1.0E-5_wp,     24.0_wp,    0.97_wp,  & !  7
                      1.0E-4_wp, 1.0E-5_wp,     25.0_wp,    0.94_wp,  & !  8
                      1.0E-4_wp, 1.0E-5_wp,     26.0_wp,    0.98_wp,  & !  9                                  
                      1.0E-4_wp, 1.0E-5_wp,     27.0_wp,    0.93_wp,  & ! 10
                      1.0E-4_wp, 1.0E-5_wp,     28.0_wp,    0.97_wp,  & ! 11
                      1.0E-4_wp, 1.0E-5_wp,     29.0_wp,    0.97_wp,  & ! 12
                      1.0E-4_wp, 1.0E-5_wp,     30.0_wp,    0.97_wp,  & ! 13
                      1.0E-4_wp, 1.0E-5_wp,     31.0_wp,    0.94_wp,  & ! 14
                      1.0E-4_wp, 1.0E-5_wp,     32.0_wp,    0.98_wp   & ! 15
                      /), (/ 4, 15 /) )                             
!
!-- Pavement subsurface parameters part 1: thermal conductivity (W/m/K)
!--   0.0-0.01, 0.01-0.03, 0.03-0.07, 0.07-0.15, 0.15-0.30, 0.30-0.50,    0.50-1.25,    1.25-3.00
    REAL(wp), DIMENSION(0:7,1:15), PARAMETER :: pavement_subsurface_pars_1 = RESHAPE( (/ &
       1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp, 9999999.9_wp, 9999999.9_wp, & !  1
       1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp, 9999999.9_wp, 9999999.9_wp, & !  2
       1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp, 9999999.9_wp, 9999999.9_wp, & !  3
       1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp, 9999999.9_wp, 9999999.9_wp, & !  4
       1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp, 9999999.9_wp, 9999999.9_wp, & !  5
       1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp, 9999999.9_wp, 9999999.9_wp, & !  6
       1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp, 9999999.9_wp, 9999999.9_wp, & !  7
       1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp, 9999999.9_wp, 9999999.9_wp, & !  8
       1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp, 9999999.9_wp, 9999999.9_wp, & !  9
       1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp, 9999999.9_wp, 9999999.9_wp, & ! 10
       1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp, 9999999.9_wp, 9999999.9_wp, & ! 11
       1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp, 9999999.9_wp, 9999999.9_wp, & ! 12
       1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp, 9999999.9_wp, 9999999.9_wp, & ! 13
       1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp, 9999999.9_wp, 9999999.9_wp, & ! 14
       1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp,   1.00_wp, 9999999.9_wp, 9999999.9_wp  & ! 15
       /), (/ 8, 15 /) )

!
!-- Pavement subsurface parameters part 2: volumetric heat capacity (J/m3/K)
!--     0.0-0.01, 0.01-0.03, 0.03-0.07, 0.07-0.15, 0.15-0.30, 0.30-0.50,    0.50-1.25,    1.25-3.00
    REAL(wp), DIMENSION(0:7,1:15), PARAMETER :: pavement_subsurface_pars_2 = RESHAPE( (/ &
       1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 9999999.9_wp, 9999999.9_wp, & !  1
       1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 9999999.9_wp, 9999999.9_wp, & !  2
       1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 9999999.9_wp, 9999999.9_wp, & !  3
       1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 9999999.9_wp, 9999999.9_wp, & !  4
       1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 9999999.9_wp, 9999999.9_wp, & !  5
       1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 9999999.9_wp, 9999999.9_wp, & !  6
       1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 9999999.9_wp, 9999999.9_wp, & !  7
       1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 9999999.9_wp, 9999999.9_wp, & !  8
       1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 9999999.9_wp, 9999999.9_wp, & !  9
       1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 9999999.9_wp, 9999999.9_wp, & ! 10
       1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 9999999.9_wp, 9999999.9_wp, & ! 11
       1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 9999999.9_wp, 9999999.9_wp, & ! 12
       1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 9999999.9_wp, 9999999.9_wp, & ! 13
       1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 9999999.9_wp, 9999999.9_wp, & ! 14
       1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 1.94E6_wp, 9999999.9_wp, 9999999.9_wp  & ! 15
                           /), (/ 8, 15 /) )
 
!
!-- TO BE FILLED
!-- Water parameters                    temperature,     z0,      z0h, albedo_type, emissivity,
    REAL(wp), DIMENSION(0:6,1:5), PARAMETER :: water_pars = RESHAPE( (/ &
       283.0_wp, 0.01_wp, 0.001_wp, 1.0E10_wp, 1.0E10_wp, 1.0_wp, 0.99_wp, & ! 1
       283.0_wp, 0.01_wp, 0.001_wp, 1.0E10_wp, 1.0E10_wp, 1.0_wp, 0.99_wp, & ! 2
       283.0_wp, 0.01_wp, 0.001_wp, 1.0E10_wp, 1.0E10_wp, 1.0_wp, 0.99_wp, & ! 3
       283.0_wp, 0.01_wp, 0.001_wp, 1.0E10_wp, 1.0E10_wp, 1.0_wp, 0.99_wp, & ! 4
       283.0_wp, 0.01_wp, 0.001_wp, 1.0E10_wp, 1.0E10_wp, 1.0_wp, 0.99_wp  & ! 5
                                                                     /), (/ 7, 5 /) )                                                                   
                                                                                                                                      
    SAVE


    PRIVATE

    
!
!-- Public functions
    PUBLIC lsm_boundary_condition, lsm_check_data_output,                      &
           lsm_check_data_output_pr,                                           &
           lsm_check_parameters, lsm_define_netcdf_grid, lsm_3d_data_averaging,& 
           lsm_data_output_2d, lsm_data_output_3d, lsm_energy_balance,         &
           lsm_header, lsm_init, lsm_init_arrays, lsm_parin, lsm_soil_model,   &
           lsm_swap_timelevel, lsm_rrd_local, lsm_wrd_local
! !vegetat
!-- Public parameters, constants and initial values
    PUBLIC aero_resist_kray, skip_time_do_lsm

!
!-- Public grid variables
    PUBLIC nzb_soil, nzs, nzt_soil, zs

!
!-- Public prognostic variables
    PUBLIC m_soil_h, t_soil_h

    INTERFACE lsm_boundary_condition
       MODULE PROCEDURE lsm_boundary_condition
    END INTERFACE lsm_boundary_condition

    INTERFACE lsm_check_data_output
       MODULE PROCEDURE lsm_check_data_output
    END INTERFACE lsm_check_data_output
    
    INTERFACE lsm_check_data_output_pr
       MODULE PROCEDURE lsm_check_data_output_pr
    END INTERFACE lsm_check_data_output_pr
    
    INTERFACE lsm_check_parameters
       MODULE PROCEDURE lsm_check_parameters
    END INTERFACE lsm_check_parameters
    
    INTERFACE lsm_3d_data_averaging
       MODULE PROCEDURE lsm_3d_data_averaging
    END INTERFACE lsm_3d_data_averaging

    INTERFACE lsm_data_output_2d
       MODULE PROCEDURE lsm_data_output_2d
    END INTERFACE lsm_data_output_2d

    INTERFACE lsm_data_output_3d
       MODULE PROCEDURE lsm_data_output_3d
    END INTERFACE lsm_data_output_3d

    INTERFACE lsm_define_netcdf_grid
       MODULE PROCEDURE lsm_define_netcdf_grid
    END INTERFACE lsm_define_netcdf_grid

    INTERFACE lsm_energy_balance
       MODULE PROCEDURE lsm_energy_balance
    END INTERFACE lsm_energy_balance

    INTERFACE lsm_header
       MODULE PROCEDURE lsm_header
    END INTERFACE lsm_header
    
    INTERFACE lsm_init
       MODULE PROCEDURE lsm_init
    END INTERFACE lsm_init

    INTERFACE lsm_init_arrays
       MODULE PROCEDURE lsm_init_arrays
    END INTERFACE lsm_init_arrays
    
    INTERFACE lsm_parin
       MODULE PROCEDURE lsm_parin
    END INTERFACE lsm_parin
    
    INTERFACE lsm_soil_model
       MODULE PROCEDURE lsm_soil_model
    END INTERFACE lsm_soil_model

    INTERFACE lsm_swap_timelevel
       MODULE PROCEDURE lsm_swap_timelevel
    END INTERFACE lsm_swap_timelevel

    INTERFACE lsm_rrd_local
       MODULE PROCEDURE lsm_rrd_local
    END INTERFACE lsm_rrd_local

    INTERFACE lsm_wrd_local
       MODULE PROCEDURE lsm_wrd_local
    END INTERFACE lsm_wrd_local

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Set internal Neumann boundary condition at outer soil grid points 
!> for temperature and humidity. 
!------------------------------------------------------------------------------!
 SUBROUTINE lsm_boundary_condition
 
    IMPLICIT NONE

    INTEGER(iwp) :: i      !< grid index x-direction
    INTEGER(iwp) :: ioff   !< offset index x-direction indicating location of soil grid point
    INTEGER(iwp) :: j      !< grid index y-direction
    INTEGER(iwp) :: joff   !< offset index x-direction indicating location of soil grid point
    INTEGER(iwp) :: k      !< grid index z-direction
    INTEGER(iwp) :: koff   !< offset index x-direction indicating location of soil grid point
    INTEGER(iwp) :: l      !< running index surface-orientation
    INTEGER(iwp) :: m      !< running index surface elements

    koff = surf_lsm_h%koff
    DO  m = 1, surf_lsm_h%ns
       i = surf_lsm_h%i(m)
       j = surf_lsm_h%j(m)
       k = surf_lsm_h%k(m)
       pt(k+koff,j,i) = pt(k,j,i)
    ENDDO

    DO  l = 0, 3
       ioff = surf_lsm_v(l)%ioff
       joff = surf_lsm_v(l)%joff
       DO  m = 1, surf_lsm_v(l)%ns
          i = surf_lsm_v(l)%i(m)
          j = surf_lsm_v(l)%j(m)
          k = surf_lsm_v(l)%k(m)
          pt(k,j+joff,i+ioff) = pt(k,j,i)
       ENDDO
    ENDDO
!
!-- In case of humidity, set boundary conditions also for q and vpt.
    IF ( humidity )  THEN
       koff = surf_lsm_h%koff
       DO  m = 1, surf_lsm_h%ns
          i = surf_lsm_h%i(m)
          j = surf_lsm_h%j(m)
          k = surf_lsm_h%k(m)
          q(k+koff,j,i)   = q(k,j,i)
          vpt(k+koff,j,i) = vpt(k,j,i)
       ENDDO

       DO  l = 0, 3
          ioff = surf_lsm_v(l)%ioff
          joff = surf_lsm_v(l)%joff
          DO  m = 1, surf_lsm_v(l)%ns
             i = surf_lsm_v(l)%i(m)
             j = surf_lsm_v(l)%j(m)
             k = surf_lsm_v(l)%k(m)
             q(k,j+joff,i+ioff)   = q(k,j,i)
             vpt(k,j+joff,i+ioff) = vpt(k,j,i)
          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE lsm_boundary_condition

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output for land surface model
!------------------------------------------------------------------------------!
 SUBROUTINE lsm_check_data_output( var, unit, i, ilen, k )
 
 
    USE control_parameters,                                                    &
        ONLY:  data_output, message_string

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  unit  !< 
    CHARACTER (LEN=*) ::  var   !<

    INTEGER(iwp) :: i
    INTEGER(iwp) :: ilen   
    INTEGER(iwp) :: k

    SELECT CASE ( TRIM( var ) )

       CASE ( 'm_soil' )
          IF (  .NOT.  land_surface )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //     &
                      'res land_surface = .TRUE.'
             CALL message( 'lsm_check_data_output', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          unit = 'm3/m3'
           
       CASE ( 't_soil' )
          IF (  .NOT.  land_surface )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //     &
                      'res land_surface = .TRUE.'
             CALL message( 'lsm_check_data_output', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          unit = 'K'   
             
       CASE ( 'lai*', 'c_liq*', 'c_soil*', 'c_veg*', 'm_liq*',                 &
              'qsws_liq*', 'qsws_soil*', 'qsws_veg*', 'r_s*' )
          IF ( k == 0  .OR.  data_output(i)(ilen-2:ilen) /= '_xy' )  THEN
             message_string = 'illegal value for data_output: "' //            &
                              TRIM( var ) // '" & only 2d-horizontal ' //      &
                              'cross sections are allowed for this value'
             CALL message( 'lsm_check_data_output', 'PA0111', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'lai*'  .AND.  .NOT.  land_surface )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //     &
                              'res land_surface = .TRUE.'
             CALL message( 'lsm_check_data_output', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'c_liq*'  .AND.  .NOT.  land_surface )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //     &
                              'res land_surface = .TRUE.'
             CALL message( 'lsm_check_data_output', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'c_soil*'  .AND.  .NOT.  land_surface )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //     &
                              'res land_surface = .TRUE.'
             CALL message( 'lsm_check_data_output', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'c_veg*'  .AND.  .NOT. land_surface )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //     &
                              'res land_surface = .TRUE.'
             CALL message( 'lsm_check_data_output', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'm_liq*'  .AND.  .NOT.  land_surface )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //     &
                              'res land_surface = .TRUE.'
             CALL message( 'lsm_check_data_output', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'qsws_liq*'  .AND.  .NOT. land_surface )         &
          THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //     &
                              'res land_surface = .TRUE.'
             CALL message( 'lsm_check_data_output', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'qsws_soil*'  .AND.  .NOT.  land_surface )       &
          THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //     &
                              'res land_surface = .TRUE.'
             CALL message( 'lsm_check_data_output', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'qsws_veg*'  .AND.  .NOT. land_surface )         &
          THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //     &
                              'res land_surface = .TRUE.'
             CALL message( 'lsm_check_data_output', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( TRIM( var ) == 'r_s*'  .AND.  .NOT.  land_surface )             &
          THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //     &
                              'res land_surface = .TRUE.'
             CALL message( 'lsm_check_data_output', 'PA0404', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( TRIM( var ) == 'lai*'   )      unit = 'none' 
          IF ( TRIM( var ) == 'c_liq*' )      unit = 'none'
          IF ( TRIM( var ) == 'c_soil*')      unit = 'none'
          IF ( TRIM( var ) == 'c_veg*' )      unit = 'none'
          IF ( TRIM( var ) == 'm_liq*'     )  unit = 'm'
          IF ( TRIM( var ) == 'qsws_liq*'  )  unit = 'W/m2'
          IF ( TRIM( var ) == 'qsws_soil*' )  unit = 'W/m2'
          IF ( TRIM( var ) == 'qsws_veg*'  )  unit = 'W/m2'
          IF ( TRIM( var ) == 'r_s*')         unit = 's/m' 
             
       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE lsm_check_data_output



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output of profiles for land surface model
!------------------------------------------------------------------------------!
 SUBROUTINE lsm_check_data_output_pr( variable, var_count, unit, dopr_unit )
 
    USE control_parameters,                                                    &
        ONLY:  data_output_pr, message_string

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
       
       CASE ( 't_soil', '#t_soil' )
          IF (  .NOT.  land_surface )  THEN
             message_string = 'data_output_pr = ' //                           &
                              TRIM( data_output_pr(var_count) ) // ' is' //    &
                              'not implemented for land_surface = .FALSE.'
             CALL message( 'lsm_check_data_output_pr', 'PA0402', 1, 2, 0, 6, 0 )
          ELSE
             dopr_index(var_count) = 89
             dopr_unit     = 'K'
             hom(0:nzs-1,2,89,:)  = SPREAD( - zs(nzb_soil:nzt_soil), 2, statistic_regions+1 )
             IF ( data_output_pr(var_count)(1:1) == '#' )  THEN
                dopr_initial_index(var_count) = 90
                hom(0:nzs-1,2,90,:)   = SPREAD( - zs(nzb_soil:nzt_soil), 2, statistic_regions+1 )
                data_output_pr(var_count)     = data_output_pr(var_count)(2:)
             ENDIF
             unit = dopr_unit
          ENDIF

       CASE ( 'm_soil', '#m_soil' )
          IF (  .NOT.  land_surface )  THEN
             message_string = 'data_output_pr = ' //                           &
                              TRIM( data_output_pr(var_count) ) // ' is' //    &
                              ' not implemented for land_surface = .FALSE.'
             CALL message( 'lsm_check_data_output_pr', 'PA0402', 1, 2, 0, 6, 0 )
          ELSE
             dopr_index(var_count) = 91
             dopr_unit     = 'm3/m3'
             hom(0:nzs-1,2,91,:)  = SPREAD( - zs(nzb_soil:nzt_soil), 2, statistic_regions+1 )
             IF ( data_output_pr(var_count)(1:1) == '#' )  THEN
                dopr_initial_index(var_count) = 92
                hom(0:nzs-1,2,92,:)   = SPREAD( - zs(nzb_soil:nzt_soil), 2, statistic_regions+1 )
                data_output_pr(var_count)     = data_output_pr(var_count)(2:)
             ENDIF
             unit = dopr_unit
          ENDIF


       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE lsm_check_data_output_pr
 
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for land surface model
!------------------------------------------------------------------------------!
 SUBROUTINE lsm_check_parameters

    USE control_parameters,                                                    &
        ONLY:  bc_pt_b, bc_q_b, constant_flux_layer, message_string,           &
               most_method
                     
    
    IMPLICIT NONE

    INTEGER(iwp) ::  k        !< running index, z-dimension

!
!-- Check for a valid setting of surface_type. The default value is 'netcdf'.
!-- In that case, the surface types are read from NetCDF file
    IF ( TRIM( surface_type ) /= 'vegetation'  .AND.                           &
         TRIM( surface_type ) /= 'pavement'    .AND.                           &
         TRIM( surface_type ) /= 'water'       .AND.                           &
         TRIM( surface_type ) /= 'netcdf' )  THEN  
       message_string = 'unknown surface type: surface_type = "' //            &
                        TRIM( surface_type ) // '"'
       CALL message( 'lsm_check_parameters', 'PA0019', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Dirichlet boundary conditions are required as the surface fluxes are
!-- calculated from the temperature/humidity gradients in the land surface
!-- model
    IF ( bc_pt_b == 'neumann'  .OR.  bc_q_b == 'neumann' )  THEN
       message_string = 'lsm requires setting of'//                            &
                        'bc_pt_b = "dirichlet" and '//                         &
                        'bc_q_b  = "dirichlet"'
       CALL message( 'lsm_check_parameters', 'PA0399', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( .NOT. TRIM(constant_flux_layer) == 'bottom' )  THEN
       message_string = 'lsm requires '//                                      &
                        'constant_flux_layer = bottom'
       CALL message( 'lsm_check_parameters', 'PA0400', 1, 2, 0, 6, 0 )
    ENDIF
    
    IF (  .NOT.  radiation )  THEN
       message_string = 'lsm requires '//                                      &
                        'the radiation model to be switched on'
       CALL message( 'lsm_check_parameters', 'PA0400', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( TRIM( surface_type ) == 'vegetation' )  THEN
    
       IF ( vegetation_type == 0 )  THEN
          IF ( min_canopy_resistance == 9999999.9_wp )  THEN
             message_string = 'vegetation_type = 0 (user defined)'//           &
                              'requires setting of min_canopy_resistance'//    &
                              '/= 9999999.9'
             CALL message( 'lsm_check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( leaf_area_index == 9999999.9_wp )  THEN
             message_string = 'vegetation_type = 0 (user_defined)'//           &
                              'requires setting of leaf_area_index'//          &
                              '/= 9999999.9'
             CALL message( 'lsm_check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( vegetation_coverage == 9999999.9_wp )  THEN
             message_string = 'vegetation_type = 0 (user_defined)'//           &
                              'requires setting of vegetation_coverage'//      &
                              '/= 9999999.9'
                CALL message( 'lsm_check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( canopy_resistance_coefficient == 9999999.9_wp)  THEN
             message_string = 'vegetation_type = 0 (user_defined)'//           &
                              'requires setting of'//                          &
                              'canopy_resistance_coefficient /= 9999999.9'
             CALL message( 'lsm_check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( lambda_surface_stable == 9999999.9_wp )  THEN
             message_string = 'vegetation_type = 0 (user_defined)'//           &
                              'requires setting of lambda_surface_stable'//    &
                              '/= 9999999.9'
             CALL message( 'lsm_check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( lambda_surface_unstable == 9999999.9_wp )  THEN
             message_string = 'vegetation_type = 0 (user_defined)'//           &
                              'requires setting of lambda_surface_unstable'//  &
                              '/= 9999999.9'
             CALL message( 'lsm_check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( f_shortwave_incoming == 9999999.9_wp )  THEN
             message_string = 'vegetation_type = 0 (user_defined)'//           &
                              'requires setting of f_shortwave_incoming'//     &
                              '/= 9999999.9'
             CALL message( 'lsm_check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( z0_vegetation == 9999999.9_wp )  THEN
             message_string = 'vegetation_type = 0 (user_defined)'//           &
                              'requires setting of z0_vegetation'//            &
                              '/= 9999999.9'
             CALL message( 'lsm_check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( z0h_vegetation == 9999999.9_wp )  THEN
             message_string = 'vegetation_type = 0 (user_defined)'//           &
                              'requires setting of z0h_vegetation'//           &
                              '/= 9999999.9'
             CALL message( 'lsm_check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       IF ( vegetation_type == 1 )  THEN
          IF ( vegetation_coverage /= 9999999.9_wp  .AND.  vegetation_coverage &
               /= 0.0_wp )  THEN
             message_string = 'vegetation_type = 1 (bare soil)'//              &
                              ' requires vegetation_coverage = 0'
             CALL message( 'lsm_check_parameters', 'PA0294', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
  
    ENDIF
    
    IF ( TRIM( surface_type ) == 'water' )  THEN

       IF ( TRIM( most_method ) == 'lookup' )  THEN    
          WRITE( message_string, * ) 'surface_type = ', surface_type,          &
                                     ' is not allowed in combination with ',   &
                                     'most_method = ', most_method
          CALL message( 'lsm_check_parameters', 'PA0414', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( water_type == 0 )  THEN  
       
          IF ( z0_water == 9999999.9_wp )  THEN
             message_string = 'water_type = 0 (user_defined)'//                &
                              'requires setting of z0_water'//                 &
                              '/= 9999999.9'
             CALL message( 'lsm_check_parameters', 'PA0415', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( z0h_water == 9999999.9_wp )  THEN
             message_string = 'water_type = 0 (user_defined)'//                &
                              'requires setting of z0h_water'//                &
                              '/= 9999999.9'
             CALL message( 'lsm_check_parameters', 'PA0392', 1, 2, 0, 6, 0 )
          ENDIF
          
          IF ( water_temperature == 9999999.9_wp )  THEN
             message_string = 'water_type = 0 (user_defined)'//                &
                              'requires setting of water_temperature'//        &
                              '/= 9999999.9'
             CALL message( 'lsm_check_parameters', 'PA0379', 1, 2, 0, 6, 0 )
          ENDIF       
          
       ENDIF
       
    ENDIF
    
    IF ( TRIM( surface_type ) == 'pavement' )  THEN

       IF ( ANY( dz_soil /= 9999999.9_wp )  .AND.  pavement_type /= 0 )  THEN
          message_string = 'non-default setting of dz_soil '//                  &
                           'does not allow to use pavement_type /= 0)'
             CALL message( 'lsm_check_parameters', 'PA0341', 1, 2, 0, 6, 0 )
          ENDIF

       IF ( pavement_type == 0 )  THEN  
       
          IF ( z0_pavement == 9999999.9_wp )  THEN
             message_string = 'pavement_type = 0 (user_defined)'//             &
                              'requires setting of z0_pavement'//              &
                              '/= 9999999.9'
             CALL message( 'lsm_check_parameters', 'PA0352', 1, 2, 0, 6, 0 )
          ENDIF
          
          IF ( z0h_pavement == 9999999.9_wp )  THEN
             message_string = 'pavement_type = 0 (user_defined)'//             &
                              'requires setting of z0h_pavement'//             &
                              '/= 9999999.9'
             CALL message( 'lsm_check_parameters', 'PA0353', 1, 2, 0, 6, 0 )
          ENDIF
          
          IF ( pavement_heat_conduct == 9999999.9_wp )  THEN
             message_string = 'pavement_type = 0 (user_defined)'//             &
                              'requires setting of pavement_heat_conduct'//    &
                              '/= 9999999.9'
             CALL message( 'lsm_check_parameters', 'PA0342', 1, 2, 0, 6, 0 )
          ENDIF 
          
           IF ( pavement_heat_capacity == 9999999.9_wp )  THEN
             message_string = 'pavement_type = 0 (user_defined)'//             &
                              'requires setting of pavement_heat_capacity'//   &
                              '/= 9999999.9'
             CALL message( 'lsm_check_parameters', 'PA0139', 1, 2, 0, 6, 0 )
          ENDIF 

          IF ( pavement_depth_level == 0 )  THEN
             message_string = 'pavement_type = 0 (user_defined)'//             &
                              'requires setting of pavement_depth_level'//     &
                              '/= 0'
             CALL message( 'lsm_check_parameters', 'PA0474', 1, 2, 0, 6, 0 )
          ENDIF 

       ENDIF 
    
    ENDIF

    IF ( TRIM( surface_type ) == 'netcdf' )  THEN
       IF ( ANY( water_type_f%var /= water_type_f%fill )  .AND.                &
            TRIM( most_method ) == 'lookup' )  THEN    
          WRITE( message_string, * ) 'water-surfaces are not allowed in ' //   &
                                     'combination with most_method = ',        &
                                     TRIM( most_method )
          CALL message( 'lsm_check_parameters', 'PA0999', 2, 2, 0, 6, 0 )
       ENDIF
!
!--    MS: Some problme here, after calling message everythings stucks at
!--        MPI_FINALIZE call. 
       IF ( ANY( pavement_type_f%var /= pavement_type_f%fill )  .AND.           &
            ANY( dz_soil /= 9999999.9_wp ) )  THEN
          message_string = 'pavement-surfaces are not allowed in ' //           &
                           'combination with a non-default setting of dz_soil'
          CALL message( 'lsm_check_parameters', 'PA0999', 2, 2, 0, 6, 0 )
       ENDIF
    ENDIF
    
!
!-- Temporary message as long as NetCDF input is not available
    IF ( TRIM( surface_type ) == 'netcdf'  .AND.  .NOT.  input_pids_static )   &
    THEN
       message_string = 'surface_type = netcdf requires static input file.'
       CALL message( 'lsm_check_parameters', 'PA0465', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( soil_type == 0 )  THEN

       IF ( alpha_vangenuchten == 9999999.9_wp )  THEN
          message_string = 'soil_type = 0 (user_defined)'//                    &
                           'requires setting of alpha_vangenuchten'//          &
                           '/= 9999999.9'
          CALL message( 'lsm_check_parameters', 'PA0403', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( l_vangenuchten == 9999999.9_wp )  THEN
          message_string = 'soil_type = 0 (user_defined)'//                    &
                           'requires setting of l_vangenuchten'//              &
                           '/= 9999999.9'
          CALL message( 'lsm_check_parameters', 'PA0403', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( n_vangenuchten == 9999999.9_wp )  THEN
          message_string = 'soil_type = 0 (user_defined)'//                    &
                           'requires setting of n_vangenuchten'//              &
                           '/= 9999999.9'
          CALL message( 'lsm_check_parameters', 'PA0403', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( hydraulic_conductivity == 9999999.9_wp )  THEN
          message_string = 'soil_type = 0 (user_defined)'//                    &
                           'requires setting of hydraulic_conductivity'//      &
                           '/= 9999999.9'
          CALL message( 'lsm_check_parameters', 'PA0403', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( saturation_moisture == 9999999.9_wp )  THEN
          message_string = 'soil_type = 0 (user_defined)'//                    &
                           'requires setting of saturation_moisture'//         &
                           '/= 9999999.9'
          CALL message( 'lsm_check_parameters', 'PA0403', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( field_capacity == 9999999.9_wp )  THEN
          message_string = 'soil_type = 0 (user_defined)'//                    &
                           'requires setting of field_capacity'//              &
                           '/= 9999999.9'
          CALL message( 'lsm_check_parameters', 'PA0403', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( wilting_point == 9999999.9_wp )  THEN
          message_string = 'soil_type = 0 (user_defined)'//                    &
                           'requires setting of wilting_point'//               &
                           '/= 9999999.9'
          CALL message( 'lsm_check_parameters', 'PA0403', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( residual_moisture == 9999999.9_wp )  THEN
          message_string = 'soil_type = 0 (user_defined)'//                    &
                           'requires setting of residual_moisture'//           &
                           '/= 9999999.9'
          CALL message( 'lsm_check_parameters', 'PA0403', 1, 2, 0, 6, 0 )
       ENDIF

    ENDIF


!!! these checks are not needed for water surfaces??

!
!-- Determine number of soil layers to be used and check whether an appropriate
!-- root fraction is prescribed
    nzb_soil = 0
    nzt_soil = -1
    IF ( ALL( dz_soil == 9999999.9_wp ) )  THEN
       nzt_soil = 7
       dz_soil(nzb_soil:nzt_soil) = dz_soil_default
    ELSE
       DO k = 0, 19
          IF ( dz_soil(k) /= 9999999.9_wp )  THEN
             nzt_soil = nzt_soil + 1
          ENDIF
       ENDDO    
    ENDIF
    nzs = nzt_soil + 1

!
!-- Check whether valid soil temperatures are prescribed
    IF ( COUNT( soil_temperature /= 9999999.9_wp ) /= nzs )  THEN
       WRITE( message_string, * ) 'number of soil layers (', nzs, ') does not',&
                                  ' match to the number of layers specified',  &
                                  ' in soil_temperature (', COUNT(             &
                                   soil_temperature /= 9999999.9_wp ), ')'
          CALL message( 'lsm_check_parameters', 'PA0471', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( deep_soil_temperature == 9999999.9_wp ) THEN
          message_string = 'deep_soil_temperature is not set but must be'//    &
                           '/= 9999999.9'
          CALL message( 'lsm_check_parameters', 'PA0472', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check whether the sum of all root fractions equals one
    IF ( vegetation_type == 0 )  THEN
       IF ( SUM( root_fraction(nzb_soil:nzt_soil) ) /= 1.0_wp )  THEN
          message_string = 'vegetation_type = 0 (user_defined)'//              &
                           'requires setting of root_fraction'//               &
                           '/= 9999999.9 and SUM(root_fraction) = 1'
          CALL message( 'lsm_check_parameters', 'PA0401', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF    
    
    
!
!-- Check for proper setting of soil moisture, must not be larger than its 
!-- saturation value. 
    DO  k = nzb_soil, nzt_soil
       IF ( soil_moisture(k) > saturation_moisture )  THEN 
          message_string = 'soil_moisture must not exceed its saturation' //    &
                            ' value'
          CALL message( 'lsm_check_parameters', 'PA0458', 1, 2, 0, 6, 0 )
       ENDIF
    ENDDO
 
!
!-- Calculate grid spacings. Temperature and moisture are defined at
!-- the center of the soil layers, whereas gradients/fluxes are 
!-- defined at the edges (_layer)
!
!-- Allocate global 1D arrays
    ALLOCATE ( ddz_soil_center(nzb_soil:nzt_soil) )
    ALLOCATE ( ddz_soil(nzb_soil:nzt_soil+1) )
    ALLOCATE ( dz_soil_center(nzb_soil:nzt_soil) )
    ALLOCATE ( zs(nzb_soil:nzt_soil+1) )


    zs(nzb_soil) = 0.5_wp * dz_soil(nzb_soil)
    zs_layer(nzb_soil) = dz_soil(nzb_soil)

    DO  k = nzb_soil+1, nzt_soil
       zs_layer(k) = zs_layer(k-1) + dz_soil(k)
       zs(k) = (zs_layer(k) +  zs_layer(k-1)) * 0.5_wp
    ENDDO

    dz_soil(nzt_soil+1) = zs_layer(nzt_soil) + dz_soil(nzt_soil)
    zs(nzt_soil+1) = zs_layer(nzt_soil) + 0.5_wp * dz_soil(nzt_soil)
 
    DO  k = nzb_soil, nzt_soil-1
       dz_soil_center(k) = zs(k+1) - zs(k)
       IF ( dz_soil_center(k) <= 0.0_wp )  THEN
          message_string = 'invalid soil layer configuration found ' //        &
                           '(dz_soil_center(k) <= 0.0)'
          CALL message( 'lsm_rrd_local', 'PA0140', 1, 2, 0, 6, 0 )
       ENDIF  
    ENDDO
  
    dz_soil_center(nzt_soil) = zs_layer(k-1) + dz_soil(k) - zs(nzt_soil)
       
    ddz_soil_center = 1.0_wp / dz_soil_center
    ddz_soil(nzb_soil:nzt_soil) = 1.0_wp / dz_soil(nzb_soil:nzt_soil)



 END SUBROUTINE lsm_check_parameters
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Solver for the energy balance at the surface.
!------------------------------------------------------------------------------!
 SUBROUTINE lsm_energy_balance( horizontal, l )

    USE diagnostic_quantities_mod,                                             &
        ONLY:  magnus 

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i         !< running index
    INTEGER(iwp) ::  i_off     !< offset to determine index of surface element, seen from atmospheric grid point, for x
    INTEGER(iwp) ::  j         !< running index
    INTEGER(iwp) ::  j_off     !< offset to determine index of surface element, seen from atmospheric grid point, for y
    INTEGER(iwp) ::  k         !< running index
    INTEGER(iwp) ::  k_off     !< offset to determine index of surface element, seen from atmospheric grid point, for z
    INTEGER(iwp) ::  ks        !< running index
    INTEGER(iwp) ::  l         !< surface-facing index
    INTEGER(iwp) ::  m         !< running index concerning wall elements

    LOGICAL      ::  horizontal !< Flag indicating horizontal or vertical surfaces

    REAL(wp) :: c_surface_tmp,& !< temporary variable for storing the volumetric heat capacity of the surface
                f1,          & !< resistance correction term 1
                f2,          & !< resistance correction term 2
                f3,          & !< resistance correction term 3
                m_min,       & !< minimum soil moisture
                e,           & !< water vapour pressure
                e_s,         & !< water vapour saturation pressure
                e_s_dt,      & !< derivate of e_s with respect to T
                tend,        & !< tendency
                dq_s_dt,     & !< derivate of q_s with respect to T
                coef_1,      & !< coef. for prognostic equation
                coef_2,      & !< coef. for prognostic equation
                f_qsws,      & !< factor for qsws
                f_qsws_veg,  & !< factor for qsws_veg
                f_qsws_soil, & !< factor for qsws_soil
                f_qsws_liq,  & !< factor for qsws_liq
                f_shf,       & !< factor for shf
                lambda_soil, & !< Thermal conductivity of the uppermost soil layer (W/m2/K)
                lambda_surface, & !< Current value of lambda_surface (W/m2/K)
                m_liq_max      !< maxmimum value of the liq. water reservoir

    TYPE(surf_type_lsm), POINTER ::  surf_t_surface
    TYPE(surf_type_lsm), POINTER ::  surf_t_surface_p
    TYPE(surf_type_lsm), POINTER ::  surf_tt_surface_m
    TYPE(surf_type_lsm), POINTER ::  surf_m_liq
    TYPE(surf_type_lsm), POINTER ::  surf_m_liq_p
    TYPE(surf_type_lsm), POINTER ::  surf_tm_liq_m

    TYPE(surf_type_lsm), POINTER ::  surf_m_soil
    TYPE(surf_type_lsm), POINTER ::  surf_t_soil

    TYPE(surf_type), POINTER  ::  surf  !< surface-date type variable

    IF ( horizontal )  THEN
       surf              => surf_lsm_h

       surf_t_surface    => t_surface_h
       surf_t_surface_p  => t_surface_h_p
       surf_tt_surface_m => tt_surface_h_m
       surf_m_liq        => m_liq_h
       surf_m_liq_p      => m_liq_h_p
       surf_tm_liq_m     => tm_liq_h_m
       surf_m_soil       => m_soil_h
       surf_t_soil       => t_soil_h
    ELSE
       surf              => surf_lsm_v(l)

       surf_t_surface    => t_surface_v(l)
       surf_t_surface_p  => t_surface_v_p(l)
       surf_tt_surface_m => tt_surface_v_m(l)
       surf_m_liq        => m_liq_v(l)
       surf_m_liq_p      => m_liq_v_p(l)
       surf_tm_liq_m     => tm_liq_v_m(l)
       surf_m_soil       => m_soil_v(l)
       surf_t_soil       => t_soil_v(l)
    ENDIF

!
!-- Index offset of surface element point with respect to adjoining 
!-- atmospheric grid point 
    k_off = surf%koff
    j_off = surf%joff
    i_off = surf%ioff

!
!-- Calculate the exner function for the current time step
    exn = ( surface_pressure / 1000.0_wp )**0.286_wp

    DO  m = 1, surf%ns

       i   = surf%i(m)            
       j   = surf%j(m)
       k   = surf%k(m)

!
!--    Define heat conductivity between surface and soil depending on surface
!--    type. For vegetation, a skin layer parameterization is used. The new
!--    parameterization uses a combination of two conductivities: a constant
!--    conductivity for the skin layer, and a conductivity according to the
!--    uppermost soil layer. For bare soil and pavements, no skin layer is 
!--    applied. In these cases, the temperature is assumed to be constant 
!--    between the surface and the first soil layer. The heat conductivity is 
!--    then derived from the soil/pavement properties. 
!--    For water surfaces, the conductivity is already set to 1E10.
!--    Moreover, the heat capacity is set. For bare soil the heat capacity is
!--    the capacity of the uppermost soil layer, for pavement it is that of
!--    the material involved.

!
!--    for vegetation type surfaces, the thermal conductivity of the soil is
!--    needed

       IF ( surf%vegetation_surface(m) )  THEN

          lambda_h_sat = lambda_h_sm**(1.0_wp - surf%m_sat(nzb_soil,m)) *      &
                         lambda_h_water ** surf_m_soil%var_2d(nzb_soil,m)
                         
          ke = 1.0_wp + LOG10( MAX( 0.1_wp, surf_m_soil%var_2d(nzb_soil,m) /   &
                                                     surf%m_sat(nzb_soil,m) ) )                   
                         
          lambda_soil = (ke * (lambda_h_sat - lambda_h_dry) + lambda_h_dry )   &
                           * ddz_soil(nzb_soil) * 2.0_wp

!
!--       When bare soil is set without a thermal conductivity (no skin layer),
!--       a heat capacity is that of the soil layer, otherwise it is a
!--       combination of the conductivities from the skin and the soil layer
          IF ( surf%lambda_surface_s(m) == 0.0_wp )  THEN 
            surf%c_surface(m) = (rho_c_soil * (1.0_wp - surf%m_sat(nzb_soil,m))&
                              + rho_c_water * surf_m_soil%var_2d(nzb_soil,m) ) &
                              * dz_soil(nzb_soil) * 0.5_wp    
            lambda_surface = lambda_soil

          ELSE IF ( surf_t_surface%var_1d(m) >= surf_t_soil%var_2d(nzb_soil,m))&
          THEN
             lambda_surface = surf%lambda_surface_s(m) * lambda_soil           &
                              / ( surf%lambda_surface_s(m) + lambda_soil )
          ELSE

             lambda_surface = surf%lambda_surface_u(m) * lambda_soil           &
                              / ( surf%lambda_surface_u(m) + lambda_soil )
          ENDIF
       ELSE
          lambda_surface = surf%lambda_surface_s(m)
       ENDIF

!
!--    Set heat capacity of the skin/surface. It is ususally zero when a skin
!--    layer is used, and non-zero otherwise.
       c_surface_tmp = surf%c_surface(m) 

!
!--    First step: calculate aerodyamic resistance. As pt, us, ts
!--    are not available for the prognostic time step, data from the last
!--    time step is used here. Note that this formulation is the
!--    equivalent to the ECMWF formulation using drag coefficients
!        IF ( cloud_physics )  THEN
!           pt1 = pt(k,j,i) + l_d_cp * pt_d_t(k) * ql(k,j,i)
!           qv1 = q(k,j,i) - ql(k,j,i)
!        ELSEIF ( cloud_droplets ) THEN
!           pt1 = pt(k,j,i) + l_d_cp * pt_d_t(k) * ql(k,j,i)
!           qv1 = q(k,j,i) 
!        ELSE
!           pt1 = pt(k,j,i)
!           IF ( humidity )  THEN
!              qv1 = q(k,j,i)
!           ELSE
!              qv1 = 0.0_wp
!           ENDIF
!        ENDIF
!
!--     Calculation of r_a for vertical surfaces
!--
!--     heat transfer coefficient for forced convection along vertical walls
!--     follows formulation in TUF3d model (Krayenhoff & Voogt, 2006)
!--            
!--       H = httc (Tsfc - Tair)
!--       httc = rw * (11.8 + 4.2 * Ueff) - 4.0
!--            
!--             rw: wall patch roughness relative to 1.0 for concrete
!--             Ueff: effective wind speed
!--             - 4.0 is a reduction of Rowley et al (1930) formulation based on
!--             Cole and Sturrock (1977)
!--           
!--             Ucan: Canyon wind speed
!--             wstar: convective velocity
!--             Qs: surface heat flux
!--             zH: height of the convective layer
!--             wstar = (g/Tcan*Qs*zH)**(1./3.)
                
!--    Effective velocity components must always 
!--    be defined at scalar grid point. The wall normal component is 
!--    obtained by simple linear interpolation. ( An alternative would
!--    be an logarithmic interpolation. )
!--    A roughness lenght of 0.001 is assumed for concrete (the inverse,
!--    1000 is used in the nominator for scaling)
!--    To do: detailed investigation which approach gives more reliable results!
!--    Please note, in case of very small friction velocity, e.g. in little 
!--    holes, the resistance can become negative. For this reason, limit r_a 
!--    to positive values. 
       IF ( horizontal  .OR.  .NOT. aero_resist_kray )  THEN
          surf%r_a(m) = ABS( ( surf%pt1(m) - surf%pt_surface(m) ) /            &
                             ( surf%ts(m) * surf%us(m) + 1.0E-20_wp ) )
       ELSE
          surf%r_a(m) = rho_cp / ( surf%z0(m) * 1000.0_wp                      &
                        * ( 11.8_wp + 4.2_wp *                                 &
                        SQRT( MAX( ( ( u(k,j,i) + u(k,j,i+1) ) * 0.5_wp )**2 + &
                                   ( ( v(k,j,i) + v(k,j+1,i) ) * 0.5_wp )**2 + &
                                   ( ( w(k,j,i) + w(k-1,j,i) ) * 0.5_wp )**2,  &
                              0.01_wp ) )                                      &
                           )  - 4.0_wp  ) 
       ENDIF
!
!--    Make sure that the resistance does not drop to zero for neutral 
!--    stratification. 
       IF ( surf%r_a(m) < 1.0_wp )  surf%r_a(m) = 1.0_wp
!
!--    Second step: calculate canopy resistance r_canopy
!--    f1-f3 here are defined as 1/f1-f3 as in ECMWF documentation
 
!--    f1: correction for incoming shortwave radiation (stomata close at 
!--    night)
       f1 = MIN( 1.0_wp, ( 0.004_wp * surf%rad_sw_in(m) + 0.05_wp ) /          &
                        (0.81_wp * (0.004_wp * surf%rad_sw_in(m)               &
                         + 1.0_wp)) )

!
!--    f2: correction for soil moisture availability to plants (the 
!--    integrated soil moisture must thus be considered here)
!--    f2 = 0 for very dry soils
       m_total = 0.0_wp
       DO  ks = nzb_soil, nzt_soil
           m_total = m_total + surf%root_fr(ks,m)                              &
                     * MAX( surf_m_soil%var_2d(ks,m), surf%m_wilt(ks,m) )
       ENDDO 

!
!--    The calculation of f2 is based on only one wilting point value for all
!--    soil layers. The value at k=nzb_soil is used here as a proxy but might
!--    need refinement in the future.
       IF ( m_total > surf%m_wilt(nzb_soil,m)  .AND.                           &
            m_total < surf%m_fc(nzb_soil,m) )  THEN
          f2 = ( m_total - surf%m_wilt(nzb_soil,m) ) /                         &
               ( surf%m_fc(nzb_soil,m) - surf%m_wilt(nzb_soil,m) )
       ELSEIF ( m_total >= surf%m_fc(nzb_soil,m) )  THEN
          f2 = 1.0_wp
       ELSE
          f2 = 1.0E-20_wp
       ENDIF

!
!--    Calculate water vapour pressure at saturation and convert to hPa
       e_s = 0.01_wp * magnus( surf_t_surface%var_1d(m) )

!
!--    f3: correction for vapour pressure deficit
       IF ( surf%g_d(m) /= 0.0_wp )  THEN
!
!--       Calculate vapour pressure
          e  = surf%qv1(m) * surface_pressure / ( surf%qv1(m) + 0.622_wp )
          f3 = EXP ( - surf%g_d(m) * (e_s - e) )
       ELSE
          f3 = 1.0_wp
       ENDIF
!
!--    Calculate canopy resistance. In case that c_veg is 0 (bare soils),
!--    this calculation is obsolete, as r_canopy is not used below.
!--    To do: check for very dry soil -> r_canopy goes to infinity
       surf%r_canopy(m) = surf%r_canopy_min(m) /                               &
                              ( surf%lai(m) * f1 * f2 * f3 + 1.0E-20_wp )
!
!--    Third step: calculate bare soil resistance r_soil.
       m_min = surf%c_veg(m) * surf%m_wilt(nzb_soil,m) +                       &
                         ( 1.0_wp - surf%c_veg(m) ) * surf%m_res(nzb_soil,m)


       f2 = ( surf_m_soil%var_2d(nzb_soil,m) - m_min ) /                       &
            ( surf%m_fc(nzb_soil,m) - m_min )
       f2 = MAX( f2, 1.0E-20_wp )
       f2 = MIN( f2, 1.0_wp     )

       surf%r_soil(m) = surf%r_soil_min(m) / f2
       
!
!--    Calculate the maximum possible liquid water amount on plants and
!--    bare surface. For vegetated surfaces, a maximum depth of 0.2 mm is
!--    assumed, while paved surfaces might hold up 1 mm of water. The 
!--    liquid water fraction for paved surfaces is calculated after 
!--    Noilhan & Planton (1989), while the ECMWF formulation is used for
!--    vegetated surfaces and bare soils.
       IF ( surf%pavement_surface(m) )  THEN
          m_liq_max = m_max_depth * 5.0_wp
          surf%c_liq(m) = MIN( 1.0_wp, ( surf_m_liq%var_1d(m) / m_liq_max)**0.67 )
       ELSE
          m_liq_max = m_max_depth * ( surf%c_veg(m) * surf%lai(m)              &
                      + ( 1.0_wp - surf%c_veg(m) ) )
          surf%c_liq(m) = MIN( 1.0_wp, surf_m_liq%var_1d(m) / m_liq_max )
       ENDIF
!
!--    Calculate saturation water vapor mixing ratio
       q_s = 0.622_wp * e_s / ( surface_pressure - e_s )
!
!--    In case of dewfall, set evapotranspiration to zero
!--    All super-saturated water is then removed from the air
       IF ( humidity  .AND.  q_s <= surf%qv1(m) )  THEN
          surf%r_canopy(m) = 0.0_wp
          surf%r_soil(m)   = 0.0_wp
       ENDIF

!
!--    Calculate coefficients for the total evapotranspiration 
!--    In case of water surface, set vegetation and soil fluxes to zero.
!--    For pavements, only evaporation of liquid water is possible.
       IF ( surf%water_surface(m) )  THEN
          f_qsws_veg  = 0.0_wp
          f_qsws_soil = 0.0_wp
          f_qsws_liq  = rho_lv / surf%r_a(m)
       ELSEIF ( surf%pavement_surface (m) )  THEN
          f_qsws_veg  = 0.0_wp
          f_qsws_soil = 0.0_wp
          f_qsws_liq  = rho_lv * surf%c_liq(m) / surf%r_a(m)
       ELSE
          f_qsws_veg  = rho_lv * surf%c_veg(m) *                               &
                            ( 1.0_wp        - surf%c_liq(m)    ) /             &
                            ( surf%r_a(m) + surf%r_canopy(m) )
          f_qsws_soil = rho_lv * (1.0_wp    - surf%c_veg(m)    ) /             &
                            ( surf%r_a(m) + surf%r_soil(m)   )
          f_qsws_liq  = rho_lv * surf%c_veg(m) * surf%c_liq(m)   /             &
                              surf%r_a(m)
       ENDIF

       f_shf  = rho_cp / surf%r_a(m)
       f_qsws = f_qsws_veg + f_qsws_soil + f_qsws_liq
!
!--    Calculate derivative of q_s for Taylor series expansion
       e_s_dt = e_s * ( 17.62_wp / ( surf_t_surface%var_1d(m) - 29.65_wp) -   &
                        17.62_wp*( surf_t_surface%var_1d(m) - 273.15_wp)      &
                       / ( surf_t_surface%var_1d(m) - 29.65_wp)**2 )

       dq_s_dt = 0.622_wp * e_s_dt / ( surface_pressure - e_s_dt )
!
!--    Calculate net radiation radiation without longwave outgoing flux because
!--    it has a dependency on surface temperature and thus enters the prognostic 
!--    equations directly
       surf%rad_net_l(m) = surf%rad_sw_in(m) - surf%rad_sw_out(m)              &
                           + surf%rad_lw_in(m)
!
!--    Calculate new skin temperature
       IF ( humidity )  THEN
!
!--       Numerator of the prognostic equation
          coef_1 = surf%rad_net_l(m) + surf%rad_lw_out_change_0(m)             &
                   * surf_t_surface%var_1d(m) - surf%rad_lw_out(m)             &
                   + f_shf * surf%pt1(m) + f_qsws * ( surf%qv1(m) - q_s        &
                   + dq_s_dt * surf_t_surface%var_1d(m) ) + lambda_surface     &
                   * surf_t_soil%var_2d(nzb_soil,m)

!
!--       Denominator of the prognostic equation
          coef_2 = surf%rad_lw_out_change_0(m) + f_qsws * dq_s_dt              &
                   + lambda_surface + f_shf / exn
       ELSE
!
!--       Numerator of the prognostic equation
          coef_1 = surf%rad_net_l(m) + surf%rad_lw_out_change_0(m)             &
                   * surf_t_surface%var_1d(m) - surf%rad_lw_out(m)             &
                   + f_shf * surf%pt1(m)  + lambda_surface                     &
                   * surf_t_soil%var_2d(nzb_soil,m)
!
!--       Denominator of the prognostic equation
          coef_2 = surf%rad_lw_out_change_0(m) + lambda_surface + f_shf / exn

       ENDIF

       tend = 0.0_wp

!
!--    Implicit solution when the surface layer has no heat capacity,
!--    otherwise use RK3 scheme.
       surf_t_surface_p%var_1d(m) = ( coef_1 * dt_3d * tsc(2) + c_surface_tmp *&
                          surf_t_surface%var_1d(m) ) / ( c_surface_tmp + coef_2&
                                             * dt_3d * tsc(2) ) 

!
!--    Add RK3 term
       IF ( c_surface_tmp /= 0.0_wp )  THEN

          surf_t_surface_p%var_1d(m) = surf_t_surface_p%var_1d(m) + dt_3d *    &
                                       tsc(3) * surf_tt_surface_m%var_1d(m)

!
!--       Calculate true tendency
          tend = ( surf_t_surface_p%var_1d(m) - surf_t_surface%var_1d(m) -     &
                   dt_3d * tsc(3) * surf_tt_surface_m%var_1d(m)) / (dt_3d  * tsc(2))
!
!--       Calculate t_surface tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                surf_tt_surface_m%var_1d(m) = tend
             ELSEIF ( intermediate_timestep_count <                            &
                      intermediate_timestep_count_max )  THEN
                surf_tt_surface_m%var_1d(m) = -9.5625_wp * tend +              &
                                               5.3125_wp * surf_tt_surface_m%var_1d(m)
             ENDIF
          ENDIF
       ENDIF

!
!--    In case of fast changes in the skin temperature, it is possible to
!--    update the radiative fluxes independently from the prescribed
!--    radiation call frequency. This effectively prevents oscillations,
!--    especially when setting skip_time_do_radiation /= 0. The threshold
!--    value of 0.2 used here is just a first guess. This method should be
!--    revised in the future as tests have shown that the threshold is 
!--    often reached, when no oscillations would occur (causes immense
!--    computing time for the radiation code).
       IF ( ABS( surf_t_surface_p%var_1d(m) - surf_t_surface%var_1d(m) )       &
            > 0.2_wp  .AND. &
            unscheduled_radiation_calls )  THEN
          force_radiation_call_l = .TRUE.
       ENDIF


!        pt(k+k_off,j+j_off,i+i_off) = surf_t_surface_p%var_1d(m) / exn  !is actually no air temperature
       surf%pt_surface(m)          = surf_t_surface_p%var_1d(m) / exn

!
!--    Calculate fluxes
       surf%rad_net_l(m) = surf%rad_net_l(m) +                                 &
                            surf%rad_lw_out_change_0(m)                        &
                          * surf_t_surface%var_1d(m) - surf%rad_lw_out(m)      &
                          - surf%rad_lw_out_change_0(m) * surf_t_surface_p%var_1d(m)

       surf%rad_net(m) = surf%rad_net_l(m)
       surf%rad_lw_out(m) = surf%rad_lw_out(m) + surf%rad_lw_out_change_0(m) * &
                     ( surf_t_surface_p%var_1d(m) - surf_t_surface%var_1d(m) )

       surf%ghf(m) = lambda_surface * ( surf_t_surface_p%var_1d(m)             &
                                             - surf_t_soil%var_2d(nzb_soil,m) )

       surf%shf(m) = - f_shf * ( surf%pt1(m) - surf%pt_surface(m) ) / cp

       IF ( humidity )  THEN
          surf%qsws(m)  = - f_qsws * ( surf%qv1(m) - q_s + dq_s_dt             &
                          * surf_t_surface%var_1d(m) - dq_s_dt *               &
                            surf_t_surface_p%var_1d(m) )

          surf%qsws_veg(m)  = - f_qsws_veg  * ( surf%qv1(m) - q_s              &
                              + dq_s_dt * surf_t_surface%var_1d(m) - dq_s_dt   &
                              * surf_t_surface_p%var_1d(m) )

          surf%qsws_soil(m) = - f_qsws_soil * ( surf%qv1(m) - q_s              &
                              + dq_s_dt * surf_t_surface%var_1d(m) - dq_s_dt   &
                              * surf_t_surface_p%var_1d(m) )

          surf%qsws_liq(m)  = - f_qsws_liq  * ( surf%qv1(m) - q_s              &
                              + dq_s_dt * surf_t_surface%var_1d(m) - dq_s_dt   &
                              * surf_t_surface_p%var_1d(m) )
       ENDIF

!
!--    Calculate the true surface resistance
       IF ( .NOT.  humidity )  THEN
          surf%r_s(m) = 1.0E10_wp
       ELSE
          surf%r_s(m) = - rho_lv * ( surf%qv1(m) - q_s + dq_s_dt               &
                          * surf_t_surface%var_1d(m) - dq_s_dt *               &
                            surf_t_surface_p%var_1d(m) ) /                     &
                            (surf%qsws(m) + 1.0E-20)  - surf%r_a(m)
       ENDIF

!
!--    Calculate change in liquid water reservoir due to dew fall or 
!--    evaporation of liquid water
       IF ( humidity )  THEN
!
!--       If precipitation is activated, add rain water to qsws_liq
!--       and qsws_soil according the the vegetation coverage.
!--       precipitation_rate is given in mm.
          IF ( precipitation )  THEN

!
!--          Add precipitation to liquid water reservoir, if possible.
!--          Otherwise, add the water to soil. In case of
!--          pavements, the exceeding water amount is implicitely removed 
!--          as runoff as qsws_soil is then not used in the soil model
             IF ( surf_m_liq%var_1d(m) /= m_liq_max )  THEN
                surf%qsws_liq(m) = surf%qsws_liq(m)                            &
                                 + surf%c_veg(m) * prr(k+k_off,j+j_off,i+i_off)&
                                 * hyrho(k+k_off)                              &
                                 * 0.001_wp * rho_l * l_v
             ELSE
                surf%qsws_soil(m) = surf%qsws_soil(m)                          &
                                 + surf%c_veg(m) * prr(k+k_off,j+j_off,i+i_off)&
                                 * hyrho(k+k_off)                              &
                                 * 0.001_wp * rho_l * l_v
             ENDIF

!--          Add precipitation to bare soil according to the bare soil
!--          coverage.
             surf%qsws_soil(m) = surf%qsws_soil(m) + ( 1.0_wp                  &
                               - surf%c_veg(m) ) * prr(k+k_off,j+j_off,i+i_off)&
                               * hyrho(k+k_off)                                &
                               * 0.001_wp * rho_l * l_v
          ENDIF

!
!--       If the air is saturated, check the reservoir water level
          IF ( surf%qsws(m) < 0.0_wp )  THEN
!
!--          Check if reservoir is full (avoid values > m_liq_max)
!--          In that case, qsws_liq goes to qsws_soil. In this 
!--          case qsws_veg is zero anyway (because c_liq = 1),       
!--          so that tend is zero and no further check is needed
             IF ( surf_m_liq%var_1d(m) == m_liq_max )  THEN
                surf%qsws_soil(m) = surf%qsws_soil(m) + surf%qsws_liq(m)

                surf%qsws_liq(m)  = 0.0_wp
             ENDIF

!
!--          In case qsws_veg becomes negative (unphysical behavior), 
!--          let the water enter the liquid water reservoir as dew on the
!--          plant
             IF ( surf%qsws_veg(m) < 0.0_wp )  THEN
                surf%qsws_liq(m) = surf%qsws_liq(m) + surf%qsws_veg(m)
                surf%qsws_veg(m) = 0.0_wp
             ENDIF
          ENDIF                    
 
          surf%qsws(m) = surf%qsws(m) / l_v
 
          tend = - surf%qsws_liq(m) * drho_l_lv
          surf_m_liq_p%var_1d(m) = surf_m_liq%var_1d(m) + dt_3d *              &
                                        ( tsc(2) * tend +                      &
                                          tsc(3) * surf_tm_liq_m%var_1d(m) )
!
!--       Check if reservoir is overfull -> reduce to maximum
!--       (conservation of water is violated here)
          surf_m_liq_p%var_1d(m) = MIN( surf_m_liq_p%var_1d(m),m_liq_max )

!
!--       Check if reservoir is empty (avoid values < 0.0)
!--       (conservation of water is violated here)
          surf_m_liq_p%var_1d(m) = MAX( surf_m_liq_p%var_1d(m), 0.0_wp )
!
!--       Calculate m_liq tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                surf_tm_liq_m%var_1d(m) = tend
             ELSEIF ( intermediate_timestep_count <                            &
                      intermediate_timestep_count_max )  THEN
                surf_tm_liq_m%var_1d(m) = -9.5625_wp * tend +                  &
                                           5.3125_wp * surf_tm_liq_m%var_1d(m)
             ENDIF
          ENDIF

       ENDIF

    ENDDO

!
!-- Make a logical OR for all processes. Force radiation call if at
!-- least one processor reached the threshold change in skin temperature
    IF ( unscheduled_radiation_calls  .AND.  intermediate_timestep_count       &
         == intermediate_timestep_count_max-1 )  THEN
#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( force_radiation_call_l, force_radiation_call,       &
                           1, MPI_LOGICAL, MPI_LOR, comm2d, ierr )
#else
       force_radiation_call = force_radiation_call_l
#endif
       force_radiation_call_l = .FALSE.
    ENDIF

!
!-- Calculate surface water vapor mixing ratio
    IF ( humidity )  THEN
       CALL calc_q_surface
    ENDIF

!
!-- Calculate new roughness lengths (for water surfaces only)
    IF ( horizontal  .AND.  .NOT. constant_roughness )  CALL calc_z0_water_surface

    CONTAINS
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of mixing ratio of the skin layer (surface). It is assumend
!> that the skin is always saturated.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_q_surface

       USE diagnostic_quantities_mod

       IMPLICIT NONE

       REAL(wp) :: resistance    !< aerodynamic and soil resistance term

       DO  m = 1, surf%ns

          i   = surf%i(m)            
          j   = surf%j(m)
          k   = surf%k(m)
!
!--       Calculate water vapour pressure at saturation and convert to hPa
          e_s = 0.01_wp * magnus( surf_t_surface_p%var_1d(m) )                   

!
!--       Calculate mixing ratio at saturation
          q_s = 0.622_wp * e_s / ( surface_pressure - e_s )

          resistance = surf%r_a(m) / ( surf%r_a(m) + surf%r_s(m) + 1E-5_wp )

!
!--       Calculate mixing ratio at surface
          IF ( cloud_physics )  THEN
             q(k+k_off,j+j_off,i+i_off) = resistance * q_s +                   &
                                        ( 1.0_wp - resistance ) *              &
                                        ( q(k,j,i) - ql(k,j,i) )
          ELSE
             q(k+k_off,j+j_off,i+i_off) = resistance * q_s +                   &
                                        ( 1.0_wp - resistance ) *              &
                                          q(k,j,i)
          ENDIF
!
!--       Update virtual potential temperature
          vpt(k+k_off,j+j_off,i+i_off) = pt(k+k_off,j+j_off,i+i_off) *         &
                     ( 1.0_wp + 0.61_wp * q(k+k_off,j+j_off,i+i_off) )

       ENDDO

    END SUBROUTINE calc_q_surface



 END SUBROUTINE lsm_energy_balance


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for land surface model
!------------------------------------------------------------------------------!
    SUBROUTINE lsm_header ( io )


       IMPLICIT NONE

       CHARACTER (LEN=86) ::  t_soil_chr          !< String for soil temperature profile
       CHARACTER (LEN=86) ::  roots_chr           !< String for root profile
       CHARACTER (LEN=86) ::  vertical_index_chr  !< String for the vertical index
       CHARACTER (LEN=86) ::  m_soil_chr          !< String for soil moisture
       CHARACTER (LEN=86) ::  soil_depth_chr      !< String for soil depth
       CHARACTER (LEN=10) ::  coor_chr            !< Temporary string
    
       INTEGER(iwp) ::  i                         !< Loop index over soil layers
 
       INTEGER(iwp), INTENT(IN) ::  io            !< Unit of the output file
 
       t_soil_chr = ''
       m_soil_chr    = ''
       soil_depth_chr  = '' 
       roots_chr        = '' 
       vertical_index_chr   = ''

       i = 1
       DO i = nzb_soil, nzt_soil
          WRITE (coor_chr,'(F10.2,7X)') soil_temperature(i)
          t_soil_chr = TRIM( t_soil_chr ) // ' ' // TRIM( coor_chr )

          WRITE (coor_chr,'(F10.2,7X)') soil_moisture(i)
          m_soil_chr = TRIM( m_soil_chr ) // ' ' // TRIM( coor_chr )

          WRITE (coor_chr,'(F10.2,7X)')  - zs(i)
          soil_depth_chr = TRIM( soil_depth_chr ) // ' '  // TRIM( coor_chr )

          WRITE (coor_chr,'(F10.2,7X)')  root_fraction(i)
          roots_chr = TRIM( roots_chr ) // ' '  // TRIM( coor_chr )

          WRITE (coor_chr,'(I10,7X)')  i
          vertical_index_chr = TRIM( vertical_index_chr ) // ' '  //           &
                               TRIM( coor_chr )
       ENDDO

!
!--    Write land surface model header
       WRITE( io,  1 )
       IF ( conserve_water_content )  THEN
          WRITE( io, 2 )
       ELSE
          WRITE( io, 3 )
       ENDIF

       IF ( vegetation_type_f%from_file )  THEN
          WRITE( io, 5 )
       ELSE
          WRITE( io, 4 ) TRIM( vegetation_type_name(vegetation_type) ),        &
                         TRIM (soil_type_name(soil_type) )
       ENDIF
       WRITE( io, 6 ) TRIM( soil_depth_chr ), TRIM( t_soil_chr ),              &
                        TRIM( m_soil_chr ), TRIM( roots_chr ),                 &
                        TRIM( vertical_index_chr )

1   FORMAT (//' Land surface model information:'/                              &
              ' ------------------------------'/)
2   FORMAT ('    --> Soil bottom is closed (water content is conserved',       &
            ', default)')
3   FORMAT ('    --> Soil bottom is open (water content is not conserved)')         
4   FORMAT ('    --> Land surface type  : ',A,/                                &
            '    --> Soil porosity type : ',A)
5   FORMAT ('    --> Land surface type  : read from file' /                    &
            '    --> Soil porosity type : read from file' )
6   FORMAT (/'    Initial soil temperature and moisture profile:'//            &
            '       Height:        ',A,'  m'/                                  &
            '       Temperature:   ',A,'  K'/                                  &
            '       Moisture:      ',A,'  m**3/m**3'/                          &
            '       Root fraction: ',A,'  '/                                   &
            '       Grid point:    ',A)


    END SUBROUTINE lsm_header


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the land surface model
!------------------------------------------------------------------------------!
    SUBROUTINE lsm_init
    
       USE control_parameters,                                                 &
           ONLY:  message_string

       USE indices,                                                            &
           ONLY:  nx, ny, topo_min_level

       USE pmc_interface,                                                      &
           ONLY:  nested_run
    
       IMPLICIT NONE

       LOGICAL      ::  init_soil_dynamically_in_child !< flag controlling initialization of soil in child domains

       INTEGER(iwp) ::  i                       !< running index
       INTEGER(iwp) ::  i_off                   !< index offset of surface element, seen from atmospheric grid point 
       INTEGER(iwp) ::  j                       !< running index
       INTEGER(iwp) ::  j_off                   !< index offset of surface element, seen from atmospheric grid point 
       INTEGER(iwp) ::  k                       !< running index
       INTEGER(iwp) ::  kn                      !< running index
       INTEGER(iwp) ::  ko                      !< running index
       INTEGER(iwp) ::  kroot                   !< running index
       INTEGER(iwp) ::  kzs                     !< running index
       INTEGER(iwp) ::  l                       !< running index surface facing
       INTEGER(iwp) ::  m                       !< running index
       INTEGER(iwp) ::  st                      !< soil-type index
       INTEGER(iwp) ::  n_soil_layers_total     !< temperature variable, stores the total number of soil layers + 4 
       INTEGER(iwp) ::  n_surf                  !< number of surface types of given surface element 

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  bound, bound_root_fr  !< temporary arrays for storing index bounds
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pr_soil_init !< temporary array used for averaging soil profiles

!
!--    Calculate Exner function
       exn = ( surface_pressure / 1000.0_wp )**0.286_wp
!
!--    If no cloud physics is used, rho_surface has not been calculated before
       IF (  .NOT.  cloud_physics )  THEN
          CALL calc_mean_profile( pt, 4 )
          rho_surface = surface_pressure * 100.0_wp / ( r_d * hom(topo_min_level+1,1,4,0) * exn )
       ENDIF

!
!--    Calculate frequently used parameters
       rho_cp    = cp * rho_surface
       rd_d_rv   = r_d / r_v
       rho_lv    = rho_surface * l_v
       drho_l_lv = 1.0_wp / (rho_l * l_v)

!
!--    Set initial values for prognostic quantities
!--    Horizontal surfaces
       tt_surface_h_m%var_1d = 0.0_wp
       tt_soil_h_m%var_2d    = 0.0_wp
       tm_soil_h_m%var_2d    = 0.0_wp
       tm_liq_h_m%var_1d     = 0.0_wp
       surf_lsm_h%c_liq      = 0.0_wp

       surf_lsm_h%ghf = 0.0_wp

       surf_lsm_h%qsws_liq  = 0.0_wp
       surf_lsm_h%qsws_soil = 0.0_wp
       surf_lsm_h%qsws_veg  = 0.0_wp

       surf_lsm_h%r_a        = 50.0_wp
       surf_lsm_h%r_s        = 50.0_wp
       surf_lsm_h%r_canopy   = 0.0_wp
       surf_lsm_h%r_soil     = 0.0_wp
!
!--    Do the same for vertical surfaces
       DO  l = 0, 3
          tt_surface_v_m(l)%var_1d = 0.0_wp
          tt_soil_v_m(l)%var_2d    = 0.0_wp
          tm_soil_v_m(l)%var_2d    = 0.0_wp
          tm_liq_v_m(l)%var_1d     = 0.0_wp
          surf_lsm_v(l)%c_liq      = 0.0_wp

          surf_lsm_v(l)%ghf = 0.0_wp

          surf_lsm_v(l)%qsws_liq  = 0.0_wp
          surf_lsm_v(l)%qsws_soil = 0.0_wp
          surf_lsm_v(l)%qsws_veg  = 0.0_wp

          surf_lsm_v(l)%r_a        = 50.0_wp
          surf_lsm_v(l)%r_s        = 50.0_wp
          surf_lsm_v(l)%r_canopy   = 0.0_wp
          surf_lsm_v(l)%r_soil     = 0.0_wp
       ENDDO

!
!--    Set initial values for prognostic soil quantities
       IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
          t_soil_h%var_2d = 0.0_wp
          m_soil_h%var_2d = 0.0_wp
          m_liq_h%var_1d  = 0.0_wp

          DO  l = 0, 3
             t_soil_v(l)%var_2d = 0.0_wp
             m_soil_v(l)%var_2d = 0.0_wp
             m_liq_v(l)%var_1d  = 0.0_wp
          ENDDO
       ENDIF
!
!--    Allocate 3D soil model arrays
!--    First, for horizontal surfaces
       ALLOCATE ( surf_lsm_h%alpha_vg(nzb_soil:nzt_soil,1:surf_lsm_h%ns)    )
       ALLOCATE ( surf_lsm_h%gamma_w_sat(nzb_soil:nzt_soil,1:surf_lsm_h%ns) )
       ALLOCATE ( surf_lsm_h%lambda_h(nzb_soil:nzt_soil,1:surf_lsm_h%ns)    )
       ALLOCATE ( surf_lsm_h%lambda_h_def(nzb_soil:nzt_soil,1:surf_lsm_h%ns))
       ALLOCATE ( surf_lsm_h%l_vg(nzb_soil:nzt_soil,1:surf_lsm_h%ns)        )
       ALLOCATE ( surf_lsm_h%m_fc(nzb_soil:nzt_soil,1:surf_lsm_h%ns)        )
       ALLOCATE ( surf_lsm_h%m_res(nzb_soil:nzt_soil,1:surf_lsm_h%ns)       )
       ALLOCATE ( surf_lsm_h%m_sat(nzb_soil:nzt_soil,1:surf_lsm_h%ns)       )
       ALLOCATE ( surf_lsm_h%m_wilt(nzb_soil:nzt_soil,1:surf_lsm_h%ns)      )
       ALLOCATE ( surf_lsm_h%n_vg(nzb_soil:nzt_soil,1:surf_lsm_h%ns)        )
       ALLOCATE ( surf_lsm_h%rho_c_total(nzb_soil:nzt_soil,1:surf_lsm_h%ns) )
       ALLOCATE ( surf_lsm_h%rho_c_total_def(nzb_soil:nzt_soil,1:surf_lsm_h%ns) )
       ALLOCATE ( surf_lsm_h%root_fr(nzb_soil:nzt_soil,1:surf_lsm_h%ns)     )
    
       surf_lsm_h%lambda_h     = 0.0_wp
!
!--    If required, allocate humidity-related variables for the soil model
       IF ( humidity )  THEN
          ALLOCATE ( surf_lsm_h%lambda_w(nzb_soil:nzt_soil,1:surf_lsm_h%ns) )
          ALLOCATE ( surf_lsm_h%gamma_w(nzb_soil:nzt_soil,1:surf_lsm_h%ns)  )  

          surf_lsm_h%lambda_w = 0.0_wp 
       ENDIF
!
!--    For vertical surfaces
       DO  l = 0, 3
          ALLOCATE ( surf_lsm_v(l)%alpha_vg(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns)    )
          ALLOCATE ( surf_lsm_v(l)%gamma_w_sat(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns) )
          ALLOCATE ( surf_lsm_v(l)%lambda_h(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns)    )
          ALLOCATE ( surf_lsm_v(l)%lambda_h_def(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns))
          ALLOCATE ( surf_lsm_v(l)%l_vg(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns)        )
          ALLOCATE ( surf_lsm_v(l)%m_fc(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns)        )
          ALLOCATE ( surf_lsm_v(l)%m_res(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns)       )
          ALLOCATE ( surf_lsm_v(l)%m_sat(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns)       )
          ALLOCATE ( surf_lsm_v(l)%m_wilt(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns)      )
          ALLOCATE ( surf_lsm_v(l)%n_vg(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns)        )
          ALLOCATE ( surf_lsm_v(l)%rho_c_total(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns) )  
          ALLOCATE ( surf_lsm_v(l)%rho_c_total_def(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns) )  
          ALLOCATE ( surf_lsm_v(l)%root_fr(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns)     )

          surf_lsm_v(l)%lambda_h     = 0.0_wp 
          
!
!--       If required, allocate humidity-related variables for the soil model
          IF ( humidity )  THEN
             ALLOCATE ( surf_lsm_v(l)%lambda_w(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns) )
             ALLOCATE ( surf_lsm_v(l)%gamma_w(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns)  )  

             surf_lsm_v(l)%lambda_w = 0.0_wp 
          ENDIF      
       ENDDO
!
!--    Allocate albedo type and emissivity for vegetation, water and pavement 
!--    fraction.
!--    Set default values at each surface element. 
       ALLOCATE ( surf_lsm_h%albedo_type(0:2,1:surf_lsm_h%ns) )
       ALLOCATE ( surf_lsm_h%emissivity(0:2,1:surf_lsm_h%ns) )
       surf_lsm_h%albedo_type = 0      
       surf_lsm_h%emissivity  = emissivity
       DO  l = 0, 3
          ALLOCATE ( surf_lsm_v(l)%albedo_type(0:2,1:surf_lsm_v(l)%ns) )
          ALLOCATE ( surf_lsm_v(l)%emissivity(0:2,1:surf_lsm_v(l)%ns)  )
          surf_lsm_v(l)%albedo_type = 0
          surf_lsm_v(l)%emissivity  = emissivity
       ENDDO
!
!--    Allocate arrays for relative surface fraction. 
!--    0 - vegetation fraction, 2 - water fraction, 1 - pavement fraction
       ALLOCATE( surf_lsm_h%frac(0:2,1:surf_lsm_h%ns) )
       surf_lsm_h%frac = 0.0_wp
       DO  l = 0, 3
          ALLOCATE( surf_lsm_v(l)%frac(0:2,1:surf_lsm_v(l)%ns) )
          surf_lsm_v(l)%frac = 0.0_wp
       ENDDO
!
!--    For vertical walls only - allocate special flag indicating if any building is on
!--    top of any natural surfaces. Used for initialization only. 
       DO  l = 0, 3
          ALLOCATE( surf_lsm_v(l)%building_covered(1:surf_lsm_v(l)%ns) )
       ENDDO
!
!--    Set flag parameter for the prescribed surface type depending on user 
!--    input. Set surface fraction to 1 for the respective type.
       SELECT CASE ( TRIM( surface_type ) )
          
          CASE ( 'vegetation' )
         
             surf_lsm_h%vegetation_surface = .TRUE.
             surf_lsm_h%frac(ind_veg_wall,:) = 1.0_wp
             DO  l = 0, 3
                surf_lsm_v(l)%vegetation_surface = .TRUE.
                surf_lsm_v(l)%frac(ind_veg_wall,:) = 1.0_wp
             ENDDO
   
          CASE ( 'water' )
             
             surf_lsm_h%water_surface = .TRUE.
             surf_lsm_h%frac(ind_wat_win,:) = 1.0_wp
!
!--          Note, vertical water surface does not really make sense.
             DO  l = 0, 3  
                surf_lsm_v(l)%water_surface   = .TRUE.
                surf_lsm_v(l)%frac(ind_wat_win,:) = 1.0_wp
             ENDDO

          CASE ( 'pavement' )
             
             surf_lsm_h%pavement_surface = .TRUE.
                surf_lsm_h%frac(ind_pav_green,:) = 1.0_wp
             DO  l = 0, 3
                surf_lsm_v(l)%pavement_surface   = .TRUE.
                surf_lsm_v(l)%frac(ind_pav_green,:) = 1.0_wp
             ENDDO

          CASE ( 'netcdf' )

             DO  m = 1, surf_lsm_h%ns
                i = surf_lsm_h%i(m)
                j = surf_lsm_h%j(m)
                IF ( vegetation_type_f%var(j,i) /= vegetation_type_f%fill )        &
                   surf_lsm_h%vegetation_surface(m) = .TRUE.
                IF ( pavement_type_f%var(j,i)   /= pavement_type_f%fill )          &
                   surf_lsm_h%pavement_surface(m) = .TRUE.
                IF ( water_type_f%var(j,i)      /= water_type_f%fill )             &
                   surf_lsm_h%water_surface(m) = .TRUE.
             ENDDO
             DO  l = 0, 3
                DO  m = 1, surf_lsm_v(l)%ns
!
!--                Only for vertical surfaces. Check if natural walls at reference
!--                grid point are covered by any building. This case, problems 
!--                with initialization will aris if index offsets are used. 
!--                In order to deal with this, set special flag. 
                   surf_lsm_v(l)%building_covered(m) = .FALSE.
                   IF ( building_type_f%from_file )  THEN
                      i = surf_lsm_v(l)%i(m) + surf_lsm_v(l)%ioff
                      j = surf_lsm_v(l)%j(m) + surf_lsm_v(l)%joff
                      IF ( building_type_f%var(j,i) /= 0 )                     &
                         surf_lsm_v(l)%building_covered(m) = .TRUE.
                   ENDIF
!
!--                Normally proceed with setting surface types.
                   i = surf_lsm_v(l)%i(m) + MERGE( 0, surf_lsm_v(l)%ioff,     &
                                            surf_lsm_v(l)%building_covered(m) )
                   j = surf_lsm_v(l)%j(m) + MERGE( 0, surf_lsm_v(l)%joff,     &
                                            surf_lsm_v(l)%building_covered(m) )
                   IF ( vegetation_type_f%var(j,i) /= vegetation_type_f%fill ) &
                      surf_lsm_v(l)%vegetation_surface(m) = .TRUE.
                   IF ( pavement_type_f%var(j,i)   /= pavement_type_f%fill )   &
                      surf_lsm_v(l)%pavement_surface(m) = .TRUE.
                   IF ( water_type_f%var(j,i)      /= water_type_f%fill )      &
                      surf_lsm_v(l)%water_surface(m) = .TRUE.
                ENDDO
             ENDDO

       END SELECT
!
!--    In case of netcdf input file, further initialize surface fractions. 
!--    At the moment only 1 surface is given at a location, so that the fraction
!--    is either 0 or 1. This will be revised later. If surface fraction 
!--    is not given in static input file, relative fractions will be derived
!--    from given surface type. In this case, only 1 type is given at a certain
!--    location (already checked).  
       IF ( input_pids_static  .AND.  surface_fraction_f%from_file )  THEN
          DO  m = 1, surf_lsm_h%ns
             i = surf_lsm_h%i(m)
             j = surf_lsm_h%j(m)
!
!--          0 - vegetation fraction, 1 - pavement fraction, 2 - water fraction             
             surf_lsm_h%frac(ind_veg_wall,m)  =                                &
                                    surface_fraction_f%frac(ind_veg_wall,j,i)         
             surf_lsm_h%frac(ind_pav_green,m) =                                &
                                    surface_fraction_f%frac(ind_pav_green,j,i)        
             surf_lsm_h%frac(ind_wat_win,m)   =                                &
                                    surface_fraction_f%frac(ind_wat_win,j,i)

          ENDDO
          DO  l = 0, 3
             DO  m = 1, surf_lsm_v(l)%ns
                i = surf_lsm_v(l)%i(m) + MERGE( 0, surf_lsm_v(l)%ioff,         &
                                                surf_lsm_v(l)%building_covered(m) ) 
                j = surf_lsm_v(l)%j(m) + MERGE( 0, surf_lsm_v(l)%joff,         &
                                                surf_lsm_v(l)%building_covered(m) ) 
!
!--             0 - vegetation fraction, 1 - pavement fraction, 2 - water fraction        
                surf_lsm_v(l)%frac(ind_veg_wall,m)  =                          &
                                    surface_fraction_f%frac(ind_veg_wall,j,i)         
                surf_lsm_v(l)%frac(ind_pav_green,m) =                          &
                                    surface_fraction_f%frac(ind_pav_green,j,i)        
                surf_lsm_v(l)%frac(ind_wat_win,m)   =                          &
                                    surface_fraction_f%frac(ind_wat_win,j,i)

             ENDDO
          ENDDO
       ELSEIF ( input_pids_static  .AND.  .NOT. surface_fraction_f%from_file ) &
       THEN
          DO  m = 1, surf_lsm_h%ns
             i = surf_lsm_h%i(m)
             j = surf_lsm_h%j(m)

             IF ( vegetation_type_f%var(j,i) /= vegetation_type_f%fill )       &        
                surf_lsm_h%frac(ind_veg_wall,m)  = 1.0_wp
             IF ( pavement_type_f%var(j,i)   /= pavement_type_f%fill   )       &        
                surf_lsm_h%frac(ind_pav_green,m) = 1.0_wp 
             IF ( water_type_f%var(j,i)      /= water_type_f%fill      )       &        
                surf_lsm_h%frac(ind_wat_win,m)   = 1.0_wp       
          ENDDO
          DO  l = 0, 3
             DO  m = 1, surf_lsm_v(l)%ns
                i = surf_lsm_v(l)%i(m) + MERGE( 0, surf_lsm_v(l)%ioff,         &
                                                surf_lsm_v(l)%building_covered(m) ) 
                j = surf_lsm_v(l)%j(m) + MERGE( 0, surf_lsm_v(l)%joff,         &
                                                surf_lsm_v(l)%building_covered(m) ) 
     
                IF ( vegetation_type_f%var(j,i) /= vegetation_type_f%fill )    &        
                   surf_lsm_v(l)%frac(ind_veg_wall,m)  = 1.0_wp
                IF ( pavement_type_f%var(j,i)   /= pavement_type_f%fill   )    &        
                   surf_lsm_v(l)%frac(ind_pav_green,m) = 1.0_wp 
                IF ( water_type_f%var(j,i)      /= water_type_f%fill      )    &        
                   surf_lsm_v(l)%frac(ind_wat_win,m)   = 1.0_wp     
             ENDDO
          ENDDO
       ENDIF
!
!--    Level 1, initialization of soil parameters.
!--    It is possible to overwrite each parameter by setting the respecticy 
!--    NAMELIST variable to a value /= 9999999.9. 
       IF ( soil_type /= 0 )  THEN  
 
          IF ( alpha_vangenuchten == 9999999.9_wp )  THEN
             alpha_vangenuchten = soil_pars(0,soil_type)
          ENDIF

          IF ( l_vangenuchten == 9999999.9_wp )  THEN
             l_vangenuchten = soil_pars(1,soil_type)
          ENDIF

          IF ( n_vangenuchten == 9999999.9_wp )  THEN
             n_vangenuchten = soil_pars(2,soil_type)            
          ENDIF

          IF ( hydraulic_conductivity == 9999999.9_wp )  THEN
             hydraulic_conductivity = soil_pars(3,soil_type)            
          ENDIF

          IF ( saturation_moisture == 9999999.9_wp )  THEN
             saturation_moisture = soil_pars(4,soil_type)           
          ENDIF

          IF ( field_capacity == 9999999.9_wp )  THEN
             field_capacity = soil_pars(5,soil_type)           
          ENDIF

          IF ( wilting_point == 9999999.9_wp )  THEN
             wilting_point = soil_pars(6,soil_type)            
          ENDIF

          IF ( residual_moisture == 9999999.9_wp )  THEN
             residual_moisture = soil_pars(7,soil_type)       
          ENDIF

       ENDIF
!
!--    Map values to the respective 2D/3D arrays
!--    Horizontal surfaces
       surf_lsm_h%alpha_vg      = alpha_vangenuchten
       surf_lsm_h%l_vg          = l_vangenuchten
       surf_lsm_h%n_vg          = n_vangenuchten 
       surf_lsm_h%gamma_w_sat   = hydraulic_conductivity
       surf_lsm_h%m_sat         = saturation_moisture
       surf_lsm_h%m_fc          = field_capacity
       surf_lsm_h%m_wilt        = wilting_point
       surf_lsm_h%m_res         = residual_moisture
       surf_lsm_h%r_soil_min    = min_soil_resistance
!
!--    Vertical surfaces
       DO  l = 0, 3
          surf_lsm_v(l)%alpha_vg      = alpha_vangenuchten
          surf_lsm_v(l)%l_vg          = l_vangenuchten
          surf_lsm_v(l)%n_vg          = n_vangenuchten 
          surf_lsm_v(l)%gamma_w_sat   = hydraulic_conductivity
          surf_lsm_v(l)%m_sat         = saturation_moisture
          surf_lsm_v(l)%m_fc          = field_capacity
          surf_lsm_v(l)%m_wilt        = wilting_point
          surf_lsm_v(l)%m_res         = residual_moisture
          surf_lsm_v(l)%r_soil_min    = min_soil_resistance
       ENDDO
!
!--    Level 2, initialization of soil parameters via soil_type read from file.
!--    Soil parameters are initialized for each (y,x)-grid point 
!--    individually using default paramter settings according to the given
!--    soil type.
       IF ( soil_type_f%from_file )  THEN
!
!--       Level of detail = 1, i.e. a homogeneous soil distribution along the 
!--       vertical dimension is assumed.
          IF ( soil_type_f%lod == 1 )  THEN
!
!--          Horizontal surfaces
             DO  m = 1, surf_lsm_h%ns
                i = surf_lsm_h%i(m)
                j = surf_lsm_h%j(m)
             
                st = soil_type_f%var_2d(j,i)
                IF ( st /= soil_type_f%fill )  THEN
                   surf_lsm_h%alpha_vg(:,m)    = soil_pars(0,st)
                   surf_lsm_h%l_vg(:,m)        = soil_pars(1,st)
                   surf_lsm_h%n_vg(:,m)        = soil_pars(2,st)
                   surf_lsm_h%gamma_w_sat(:,m) = soil_pars(3,st)
                   surf_lsm_h%m_sat(:,m)       = soil_pars(4,st)
                   surf_lsm_h%m_fc(:,m)        = soil_pars(5,st)
                   surf_lsm_h%m_wilt(:,m)      = soil_pars(6,st)
                   surf_lsm_h%m_res(:,m)       = soil_pars(7,st)
                ENDIF
             ENDDO
!
!--          Vertical surfaces ( assumes the soil type given at respective (x,y) 
             DO  l = 0, 3
                DO  m = 1, surf_lsm_v(l)%ns
                   i = surf_lsm_v(l)%i(m) + MERGE( 0, surf_lsm_v(l)%ioff,      &
                                                   surf_lsm_v(l)%building_covered(m) ) 
                   j = surf_lsm_v(l)%j(m) + MERGE( 0, surf_lsm_v(l)%joff,      &
                                                   surf_lsm_v(l)%building_covered(m) ) 

                   st = soil_type_f%var_2d(j,i)
                   IF ( st /= soil_type_f%fill )  THEN
                      surf_lsm_v(l)%alpha_vg(:,m)    = soil_pars(0,st)
                      surf_lsm_v(l)%l_vg(:,m)        = soil_pars(1,st)
                      surf_lsm_v(l)%n_vg(:,m)        = soil_pars(2,st)
                      surf_lsm_v(l)%gamma_w_sat(:,m) = soil_pars(3,st)
                      surf_lsm_v(l)%m_sat(:,m)       = soil_pars(4,st)
                      surf_lsm_v(l)%m_fc(:,m)        = soil_pars(5,st)
                      surf_lsm_v(l)%m_wilt(:,m)      = soil_pars(6,st)
                      surf_lsm_v(l)%m_res(:,m)       = soil_pars(7,st)
                   ENDIF
                ENDDO
             ENDDO
!
!--       Level of detail = 2, i.e. soil type and thus the soil parameters
!--       can be heterogeneous along the vertical dimension. 
          ELSE
!
!--          Horizontal surfaces
             DO  m = 1, surf_lsm_h%ns
                i = surf_lsm_h%i(m)
                j = surf_lsm_h%j(m)
             
                DO  k = nzb_soil, nzt_soil
                   st = soil_type_f%var_3d(k,j,i)
                   IF ( st /= soil_type_f%fill )  THEN
                      surf_lsm_h%alpha_vg(k,m)    = soil_pars(0,st)
                      surf_lsm_h%l_vg(k,m)        = soil_pars(1,st)
                      surf_lsm_h%n_vg(k,m)        = soil_pars(2,st)
                      surf_lsm_h%gamma_w_sat(k,m) = soil_pars(3,st)
                      surf_lsm_h%m_sat(k,m)       = soil_pars(4,st)
                      surf_lsm_h%m_fc(k,m)        = soil_pars(5,st)
                      surf_lsm_h%m_wilt(k,m)      = soil_pars(6,st)
                      surf_lsm_h%m_res(k,m)       = soil_pars(7,st)
                   ENDIF
                ENDDO
             ENDDO
!
!--          Vertical surfaces ( assumes the soil type given at respective (x,y) 
             DO  l = 0, 3
                DO  m = 1, surf_lsm_v(l)%ns
                   i = surf_lsm_v(l)%i(m) + MERGE( 0, surf_lsm_v(l)%ioff,      &
                                                   surf_lsm_v(l)%building_covered(m) ) 
                   j = surf_lsm_v(l)%j(m) + MERGE( 0, surf_lsm_v(l)%joff,      &
                                                   surf_lsm_v(l)%building_covered(m) ) 

                   DO  k = nzb_soil, nzt_soil
                      st = soil_type_f%var_3d(k,j,i)
                      IF ( st /= soil_type_f%fill )  THEN
                         surf_lsm_v(l)%alpha_vg(k,m)    = soil_pars(0,st)
                         surf_lsm_v(l)%l_vg(k,m)        = soil_pars(1,st)
                         surf_lsm_v(l)%n_vg(k,m)        = soil_pars(2,st)
                         surf_lsm_v(l)%gamma_w_sat(k,m) = soil_pars(3,st)
                         surf_lsm_v(l)%m_sat(k,m)       = soil_pars(4,st)
                         surf_lsm_v(l)%m_fc(k,m)        = soil_pars(5,st)
                         surf_lsm_v(l)%m_wilt(k,m)      = soil_pars(6,st)
                         surf_lsm_v(l)%m_res(k,m)       = soil_pars(7,st)
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF
!
!--    Level 3, initialization of single soil parameters at single z,x,y 
!--    position via soil_pars read from file.
       IF ( soil_pars_f%from_file )  THEN
!
!--       Level of detail = 1, i.e. a homogeneous vertical distribution of soil 
!--       parameters is assumed.
!--       Horizontal surfaces
          IF ( soil_pars_f%lod == 1 )  THEN
!
!--          Horizontal surfaces 
             DO  m = 1, surf_lsm_h%ns
                i = surf_lsm_h%i(m)
                j = surf_lsm_h%j(m)

                IF ( soil_pars_f%pars_xy(0,j,i) /= soil_pars_f%fill )              &
                   surf_lsm_h%alpha_vg(:,m)    = soil_pars_f%pars_xy(0,j,i)
                IF ( soil_pars_f%pars_xy(1,j,i) /= soil_pars_f%fill )              &
                   surf_lsm_h%l_vg(:,m)        = soil_pars_f%pars_xy(1,j,i)
                IF ( soil_pars_f%pars_xy(2,j,i) /= soil_pars_f%fill )              &
                   surf_lsm_h%n_vg(:,m)        = soil_pars_f%pars_xy(2,j,i)
                IF ( soil_pars_f%pars_xy(3,j,i) /= soil_pars_f%fill )              &
                   surf_lsm_h%gamma_w_sat(:,m) = soil_pars_f%pars_xy(3,j,i)
                IF ( soil_pars_f%pars_xy(4,j,i) /= soil_pars_f%fill )              &
                   surf_lsm_h%m_sat(:,m)       = soil_pars_f%pars_xy(4,j,i)
                IF ( soil_pars_f%pars_xy(5,j,i) /= soil_pars_f%fill )              &
                   surf_lsm_h%m_fc(:,m)        = soil_pars_f%pars_xy(5,j,i)
                IF ( soil_pars_f%pars_xy(6,j,i) /= soil_pars_f%fill )              &
                   surf_lsm_h%m_wilt(:,m)      = soil_pars_f%pars_xy(6,j,i)
                IF ( soil_pars_f%pars_xy(7,j,i) /= soil_pars_f%fill )              &
                   surf_lsm_h%m_res(:,m)       = soil_pars_f%pars_xy(7,j,i)

             ENDDO
!
!--          Vertical surfaces 
             DO  l = 0, 3
                DO  m = 1, surf_lsm_v(l)%ns
                   i = surf_lsm_v(l)%i(m) + MERGE( 0, surf_lsm_v(l)%ioff,      &
                                                   surf_lsm_v(l)%building_covered(m) ) 
                   j = surf_lsm_v(l)%j(m) + MERGE( 0, surf_lsm_v(l)%joff,      &
                                                   surf_lsm_v(l)%building_covered(m) ) 

                   IF ( soil_pars_f%pars_xy(0,j,i) /= soil_pars_f%fill )           &
                      surf_lsm_v(l)%alpha_vg(:,m)    = soil_pars_f%pars_xy(0,j,i)
                   IF ( soil_pars_f%pars_xy(1,j,i) /= soil_pars_f%fill )           &
                      surf_lsm_v(l)%l_vg(:,m)        = soil_pars_f%pars_xy(1,j,i)
                   IF ( soil_pars_f%pars_xy(2,j,i) /= soil_pars_f%fill )           &
                      surf_lsm_v(l)%n_vg(:,m)        = soil_pars_f%pars_xy(2,j,i)
                   IF ( soil_pars_f%pars_xy(3,j,i) /= soil_pars_f%fill )           &
                      surf_lsm_v(l)%gamma_w_sat(:,m) = soil_pars_f%pars_xy(3,j,i)
                   IF ( soil_pars_f%pars_xy(4,j,i) /= soil_pars_f%fill )           &
                      surf_lsm_v(l)%m_sat(:,m)       = soil_pars_f%pars_xy(4,j,i)
                   IF ( soil_pars_f%pars_xy(5,j,i) /= soil_pars_f%fill )           &
                      surf_lsm_v(l)%m_fc(:,m)        = soil_pars_f%pars_xy(5,j,i)
                   IF ( soil_pars_f%pars_xy(6,j,i) /= soil_pars_f%fill )           &
                      surf_lsm_v(l)%m_wilt(:,m)      = soil_pars_f%pars_xy(6,j,i)
                   IF ( soil_pars_f%pars_xy(7,j,i) /= soil_pars_f%fill )           &
                      surf_lsm_v(l)%m_res(:,m)       = soil_pars_f%pars_xy(7,j,i)

                ENDDO
             ENDDO
!
!--       Level of detail = 2, i.e. soil parameters can be set at each soil 
!--       layer individually. 
          ELSE
!
!--          Horizontal surfaces
             DO  m = 1, surf_lsm_h%ns
                i = surf_lsm_h%i(m)
                j = surf_lsm_h%j(m)

                DO  k = nzb_soil, nzt_soil
                   IF ( soil_pars_f%pars_xyz(0,k,j,i) /= soil_pars_f%fill )        &
                      surf_lsm_h%alpha_vg(k,m)    = soil_pars_f%pars_xyz(0,k,j,i)
                   IF ( soil_pars_f%pars_xyz(1,k,j,i) /= soil_pars_f%fill )        &
                      surf_lsm_h%l_vg(k,m)        = soil_pars_f%pars_xyz(1,k,j,i)
                   IF ( soil_pars_f%pars_xyz(2,k,j,i) /= soil_pars_f%fill )        &
                      surf_lsm_h%n_vg(k,m)        = soil_pars_f%pars_xyz(2,k,j,i)
                   IF ( soil_pars_f%pars_xyz(3,k,j,i) /= soil_pars_f%fill )        &
                      surf_lsm_h%gamma_w_sat(k,m) = soil_pars_f%pars_xyz(3,k,j,i)
                   IF ( soil_pars_f%pars_xyz(4,k,j,i) /= soil_pars_f%fill )        &
                      surf_lsm_h%m_sat(k,m)       = soil_pars_f%pars_xyz(4,k,j,i)
                   IF ( soil_pars_f%pars_xyz(5,k,j,i) /= soil_pars_f%fill )        &
                      surf_lsm_h%m_fc(k,m)        = soil_pars_f%pars_xyz(5,k,j,i)
                   IF ( soil_pars_f%pars_xyz(6,k,j,i) /= soil_pars_f%fill )        &
                      surf_lsm_h%m_wilt(k,m)      = soil_pars_f%pars_xyz(6,k,j,i)
                   IF ( soil_pars_f%pars_xyz(7,k,j,i) /= soil_pars_f%fill )        &
                      surf_lsm_h%m_res(k,m)       = soil_pars_f%pars_xyz(7,k,j,i)
                ENDDO

             ENDDO
!
!--          Vertical surfaces 
             DO  l = 0, 3
                DO  m = 1, surf_lsm_v(l)%ns
                   i = surf_lsm_v(l)%i(m) + MERGE( 0, surf_lsm_v(l)%ioff,      &
                                                   surf_lsm_v(l)%building_covered(m) ) 
                   j = surf_lsm_v(l)%j(m) + MERGE( 0, surf_lsm_v(l)%joff,      &
                                                   surf_lsm_v(l)%building_covered(m) ) 

                   DO  k = nzb_soil, nzt_soil
                      IF ( soil_pars_f%pars_xyz(0,k,j,i) /= soil_pars_f%fill )        &
                         surf_lsm_v(l)%alpha_vg(k,m)    = soil_pars_f%pars_xyz(0,k,j,i)
                      IF ( soil_pars_f%pars_xyz(1,k,j,i) /= soil_pars_f%fill )        &
                         surf_lsm_v(l)%l_vg(k,m)        = soil_pars_f%pars_xyz(1,k,j,i)
                      IF ( soil_pars_f%pars_xyz(2,k,j,i) /= soil_pars_f%fill )        &
                         surf_lsm_v(l)%n_vg(k,m)        = soil_pars_f%pars_xyz(2,k,j,i)
                      IF ( soil_pars_f%pars_xyz(3,k,j,i) /= soil_pars_f%fill )        &
                         surf_lsm_v(l)%gamma_w_sat(k,m) = soil_pars_f%pars_xyz(3,k,j,i)
                      IF ( soil_pars_f%pars_xyz(4,k,j,i) /= soil_pars_f%fill )        &
                         surf_lsm_v(l)%m_sat(k,m)       = soil_pars_f%pars_xyz(4,k,j,i)
                      IF ( soil_pars_f%pars_xyz(5,k,j,i) /= soil_pars_f%fill )        &
                         surf_lsm_v(l)%m_fc(k,m)        = soil_pars_f%pars_xyz(5,k,j,i)
                      IF ( soil_pars_f%pars_xyz(6,k,j,i) /= soil_pars_f%fill )        &
                         surf_lsm_v(l)%m_wilt(k,m)      = soil_pars_f%pars_xyz(6,k,j,i)
                      IF ( soil_pars_f%pars_xyz(7,k,j,i) /= soil_pars_f%fill )        &
                         surf_lsm_v(l)%m_res(k,m)       = soil_pars_f%pars_xyz(7,k,j,i)
                   ENDDO

                ENDDO
             ENDDO

          ENDIF
       ENDIF

!
!--    Level 1, initialization of vegetation parameters. A horizontally 
!--    homogeneous distribution is assumed here. 
       IF ( vegetation_type /= 0 )  THEN

          IF ( min_canopy_resistance == 9999999.9_wp )  THEN
             min_canopy_resistance = vegetation_pars(ind_v_rc_min,vegetation_type)
          ENDIF

          IF ( leaf_area_index == 9999999.9_wp )  THEN
             leaf_area_index = vegetation_pars(ind_v_rc_lai,vegetation_type)         
          ENDIF

          IF ( vegetation_coverage == 9999999.9_wp )  THEN
             vegetation_coverage = vegetation_pars(ind_v_c_veg,vegetation_type)     
          ENDIF

          IF ( canopy_resistance_coefficient == 9999999.9_wp )  THEN
              canopy_resistance_coefficient= vegetation_pars(ind_v_gd,vegetation_type)     
          ENDIF

          IF ( z0_vegetation == 9999999.9_wp )  THEN
             z0_vegetation  = vegetation_pars(ind_v_z0,vegetation_type) 
          ENDIF

          IF ( z0h_vegetation == 9999999.9_wp )  THEN
             z0h_vegetation = vegetation_pars(ind_v_z0qh,vegetation_type)
          ENDIF   
          
          IF ( lambda_surface_stable == 9999999.9_wp )  THEN
             lambda_surface_stable = vegetation_pars(ind_v_lambda_s,vegetation_type) 
          ENDIF

          IF ( lambda_surface_unstable == 9999999.9_wp )  THEN
             lambda_surface_unstable = vegetation_pars(ind_v_lambda_u,vegetation_type)           
          ENDIF

          IF ( f_shortwave_incoming == 9999999.9_wp )  THEN
             f_shortwave_incoming = vegetation_pars(ind_v_f_sw_in,vegetation_type)        
          ENDIF

          IF ( c_surface == 9999999.9_wp )  THEN
             c_surface = vegetation_pars(ind_v_c_surf,vegetation_type)        
          ENDIF

          IF ( albedo_type == 9999999  .AND.  albedo == 9999999.9_wp )  THEN
             albedo_type = INT(vegetation_pars(ind_v_at,vegetation_type))       
          ENDIF
   
          IF ( emissivity == 9999999.9_wp )  THEN
             emissivity = vegetation_pars(ind_v_emis,vegetation_type)     
          ENDIF

       ENDIF
!
!--    Map values onto horizontal elemements
       DO  m = 1, surf_lsm_h%ns
          IF ( surf_lsm_h%vegetation_surface(m) )  THEN
             surf_lsm_h%r_canopy_min(m)     = min_canopy_resistance
             surf_lsm_h%lai(m)              = leaf_area_index
             surf_lsm_h%c_veg(m)            = vegetation_coverage
             surf_lsm_h%g_d(m)              = canopy_resistance_coefficient
             surf_lsm_h%z0(m)               = z0_vegetation
             surf_lsm_h%z0h(m)              = z0h_vegetation
             surf_lsm_h%z0q(m)              = z0h_vegetation
             surf_lsm_h%lambda_surface_s(m) = lambda_surface_stable
             surf_lsm_h%lambda_surface_u(m) = lambda_surface_unstable
             surf_lsm_h%f_sw_in(m)          = f_shortwave_incoming
             surf_lsm_h%c_surface(m)        = c_surface
             surf_lsm_h%albedo_type(ind_veg_wall,m) = albedo_type
             surf_lsm_h%emissivity(ind_veg_wall,m)  = emissivity
          ELSE
             surf_lsm_h%lai(m)   = 0.0_wp
             surf_lsm_h%c_veg(m) = 0.0_wp
             surf_lsm_h%g_d(m)   = 0.0_wp
          ENDIF
  
       ENDDO
!
!--    Map values onto vertical elements, even though this does not make
!--    much sense.
       DO  l = 0, 3
          DO  m = 1, surf_lsm_v(l)%ns
             IF ( surf_lsm_v(l)%vegetation_surface(m) )  THEN
                surf_lsm_v(l)%r_canopy_min(m)     = min_canopy_resistance
                surf_lsm_v(l)%lai(m)              = leaf_area_index
                surf_lsm_v(l)%c_veg(m)            = vegetation_coverage
                surf_lsm_v(l)%g_d(m)              = canopy_resistance_coefficient
                surf_lsm_v(l)%z0(m)               = z0_vegetation
                surf_lsm_v(l)%z0h(m)              = z0h_vegetation
                surf_lsm_v(l)%z0q(m)              = z0h_vegetation
                surf_lsm_v(l)%lambda_surface_s(m) = lambda_surface_stable
                surf_lsm_v(l)%lambda_surface_u(m) = lambda_surface_unstable
                surf_lsm_v(l)%f_sw_in(m)          = f_shortwave_incoming
                surf_lsm_v(l)%c_surface(m)        = c_surface
                surf_lsm_v(l)%albedo_type(ind_veg_wall,m) = albedo_type
                surf_lsm_v(l)%emissivity(ind_veg_wall,m)  = emissivity
             ELSE
                surf_lsm_v(l)%lai(m)   = 0.0_wp
                surf_lsm_v(l)%c_veg(m) = 0.0_wp
                surf_lsm_v(l)%g_d(m)   = 0.0_wp
             ENDIF
          ENDDO
       ENDDO

!
!--    Level 2, initialization of vegation parameters via vegetation_type read
!--    from file. Vegetation parameters are initialized for each (y,x)-grid point 
!--    individually using default paramter settings according to the given
!--    vegetation type.
       IF ( vegetation_type_f%from_file )  THEN
!
!--       Horizontal surfaces
          DO  m = 1, surf_lsm_h%ns
             i = surf_lsm_h%i(m)
             j = surf_lsm_h%j(m)
             
             st = vegetation_type_f%var(j,i)
             IF ( st /= vegetation_type_f%fill  .AND.  st /= 0 )  THEN
                surf_lsm_h%r_canopy_min(m)     = vegetation_pars(ind_v_rc_min,st)
                surf_lsm_h%lai(m)              = vegetation_pars(ind_v_rc_lai,st)
                surf_lsm_h%c_veg(m)            = vegetation_pars(ind_v_c_veg,st)
                surf_lsm_h%g_d(m)              = vegetation_pars(ind_v_gd,st)
                surf_lsm_h%z0(m)               = vegetation_pars(ind_v_z0,st)
                surf_lsm_h%z0h(m)              = vegetation_pars(ind_v_z0qh,st)
                surf_lsm_h%z0q(m)              = vegetation_pars(ind_v_z0qh,st)
                surf_lsm_h%lambda_surface_s(m) = vegetation_pars(ind_v_lambda_s,st)
                surf_lsm_h%lambda_surface_u(m) = vegetation_pars(ind_v_lambda_u,st)
                surf_lsm_h%f_sw_in(m)          = vegetation_pars(ind_v_f_sw_in,st)
                surf_lsm_h%c_surface(m)        = vegetation_pars(ind_v_c_surf,st)
                surf_lsm_h%albedo_type(ind_veg_wall,m) = INT( vegetation_pars(ind_v_at,st) )
                surf_lsm_h%emissivity(ind_veg_wall,m)  = vegetation_pars(ind_v_emis,st)
             ENDIF
          ENDDO
!
!--       Vertical surfaces
          DO  l = 0, 3
             DO  m = 1, surf_lsm_v(l)%ns
                i = surf_lsm_v(l)%i(m) + MERGE( 0, surf_lsm_v(l)%ioff,      &
                                                surf_lsm_v(l)%building_covered(m) ) 
                j = surf_lsm_v(l)%j(m) + MERGE( 0, surf_lsm_v(l)%joff,      &
                                                   surf_lsm_v(l)%building_covered(m) ) 
             
                st = vegetation_type_f%var(j,i)
                IF ( st /= vegetation_type_f%fill  .AND.  st /= 0 )  THEN
                   surf_lsm_v(l)%r_canopy_min(m)     = vegetation_pars(ind_v_rc_min,st)
                   surf_lsm_v(l)%lai(m)              = vegetation_pars(ind_v_rc_lai,st)
                   surf_lsm_v(l)%c_veg(m)            = vegetation_pars(ind_v_c_veg,st)
                   surf_lsm_v(l)%g_d(m)              = vegetation_pars(ind_v_gd,st)
                   surf_lsm_v(l)%z0(m)               = vegetation_pars(ind_v_z0,st)
                   surf_lsm_v(l)%z0h(m)              = vegetation_pars(ind_v_z0qh,st)
                   surf_lsm_v(l)%z0q(m)              = vegetation_pars(ind_v_z0qh,st)
                   surf_lsm_v(l)%lambda_surface_s(m) = vegetation_pars(ind_v_lambda_s,st)
                   surf_lsm_v(l)%lambda_surface_u(m) = vegetation_pars(ind_v_lambda_u,st)
                   surf_lsm_v(l)%f_sw_in(m)          = vegetation_pars(ind_v_f_sw_in,st)
                   surf_lsm_v(l)%c_surface(m)        = vegetation_pars(ind_v_c_surf,st)
                   surf_lsm_v(l)%albedo_type(ind_veg_wall,m) = INT( vegetation_pars(ind_v_at,st) )
                   surf_lsm_v(l)%emissivity(ind_veg_wall,m)  = vegetation_pars(ind_v_emis,st)
                ENDIF
             ENDDO
          ENDDO
       ENDIF
!
!--    Level 3, initialization of vegation parameters at single (x,y) 
!--    position via vegetation_pars read from file.
       IF ( vegetation_pars_f%from_file )  THEN
!
!--       Horizontal surfaces 
          DO  m = 1, surf_lsm_h%ns

             i = surf_lsm_h%i(m)
             j = surf_lsm_h%j(m)
!
!--          If surface element is not a vegetation surface and any value in 
!--          vegetation_pars is given, neglect this information and give an
!--          informative message that this value will not be used.   
             IF ( .NOT. surf_lsm_h%vegetation_surface(m)  .AND.                &
                   ANY( vegetation_pars_f%pars_xy(:,j,i) /=                    &
                   vegetation_pars_f%fill ) )  THEN
                WRITE( message_string, * )                                     &
                                 'surface element at grid point (j,i) = (',    &
                                 j, i, ') is not a vegation surface, ',        &
                                 'so that information given in ',              &
                                 'vegetation_pars at this point is neglected.' 
                CALL message( 'land_surface_model_mod', 'PA0999', 0, 0, 0, 6, 0 )
             ELSE

                IF ( vegetation_pars_f%pars_xy(ind_v_rc_min,j,i) /=            &
                     vegetation_pars_f%fill )                                  &
                   surf_lsm_h%r_canopy_min(m)  =                               &
                                   vegetation_pars_f%pars_xy(ind_v_rc_min,j,i)
                IF ( vegetation_pars_f%pars_xy(ind_v_rc_lai,j,i) /=            &
                     vegetation_pars_f%fill )                                  &
                   surf_lsm_h%lai(m)           =                               &
                                   vegetation_pars_f%pars_xy(ind_v_rc_lai,j,i)
                IF ( vegetation_pars_f%pars_xy(ind_v_c_veg,j,i) /=             &
                     vegetation_pars_f%fill )                                  &
                   surf_lsm_h%c_veg(m)         =                               &
                                   vegetation_pars_f%pars_xy(ind_v_c_veg,j,i)
                IF ( vegetation_pars_f%pars_xy(ind_v_gd,j,i) /=                &
                     vegetation_pars_f%fill )                                  &
                   surf_lsm_h%g_d(m)           =                               &
                                   vegetation_pars_f%pars_xy(ind_v_gd,j,i)
                IF ( vegetation_pars_f%pars_xy(ind_v_z0,j,i) /=                &
                     vegetation_pars_f%fill )                                  &
                   surf_lsm_h%z0(m)            =                               &
                                   vegetation_pars_f%pars_xy(ind_v_z0,j,i)
                IF ( vegetation_pars_f%pars_xy(ind_v_z0qh,j,i) /=              &
                     vegetation_pars_f%fill )  THEN
                   surf_lsm_h%z0h(m)           =                               &
                                   vegetation_pars_f%pars_xy(ind_v_z0qh,j,i)
                   surf_lsm_h%z0q(m)           =                               &
                                   vegetation_pars_f%pars_xy(ind_v_z0qh,j,i)
                ENDIF
                IF ( vegetation_pars_f%pars_xy(ind_v_lambda_s,j,i) /=          &
                     vegetation_pars_f%fill )                                  &
                   surf_lsm_h%lambda_surface_s(m) =                            &
                                   vegetation_pars_f%pars_xy(ind_v_lambda_s,j,i)
                IF ( vegetation_pars_f%pars_xy(ind_v_lambda_u,j,i) /=          &
                     vegetation_pars_f%fill )                                  &
                   surf_lsm_h%lambda_surface_u(m) =                            &
                                   vegetation_pars_f%pars_xy(ind_v_lambda_u,j,i)
                IF ( vegetation_pars_f%pars_xy(ind_v_f_sw_in,j,i) /=           &
                     vegetation_pars_f%fill )                                  &
                   surf_lsm_h%f_sw_in(m)          =                            &
                                   vegetation_pars_f%pars_xy(ind_v_f_sw_in,j,i)
                IF ( vegetation_pars_f%pars_xy(ind_v_c_surf,j,i) /=            &
                     vegetation_pars_f%fill )                                  &
                   surf_lsm_h%c_surface(m)        =                            &
                                   vegetation_pars_f%pars_xy(ind_v_c_surf,j,i)
                IF ( vegetation_pars_f%pars_xy(ind_v_at,j,i) /=                &
                     vegetation_pars_f%fill )                                  &
                   surf_lsm_h%albedo_type(ind_veg_wall,m) =                    &
                                   INT( vegetation_pars_f%pars_xy(ind_v_at,j,i) )
                IF ( vegetation_pars_f%pars_xy(ind_v_emis,j,i) /=              &
                     vegetation_pars_f%fill )                                  &
                   surf_lsm_h%emissivity(ind_veg_wall,m)  =                    &
                                   vegetation_pars_f%pars_xy(ind_v_emis,j,i)
             ENDIF
          ENDDO
!
!--       Vertical surfaces 
          DO  l = 0, 3
             DO  m = 1, surf_lsm_v(l)%ns
                i = surf_lsm_v(l)%i(m) + MERGE( 0, surf_lsm_v(l)%ioff,         &
                                                surf_lsm_v(l)%building_covered(m) ) 
                j = surf_lsm_v(l)%j(m) + MERGE( 0, surf_lsm_v(l)%joff,         &
                                                surf_lsm_v(l)%building_covered(m) ) 
!
!--             If surface element is not a vegetation surface and any value in 
!--             vegetation_pars is given, neglect this information and give an
!--             informative message that this value will not be used.   
                IF ( .NOT. surf_lsm_v(l)%vegetation_surface(m)  .AND.          &
                      ANY( vegetation_pars_f%pars_xy(:,j,i) /=                 &
                      vegetation_pars_f%fill ) )  THEN
                   WRITE( message_string, * )                                  &
                                 'surface element at grid point (j,i) = (',    &
                                 j, i, ') is not a vegation surface, ',        &
                                 'so that information given in ',              &
                                 'vegetation_pars at this point is neglected.' 
                   CALL message( 'land_surface_model_mod', 'PA0999', 0, 0, 0, 6, 0 )
                ELSE

                   IF ( vegetation_pars_f%pars_xy(ind_v_rc_min,j,i) /=         &
                        vegetation_pars_f%fill )                               &
                      surf_lsm_v(l)%r_canopy_min(m)  =                         &
                                   vegetation_pars_f%pars_xy(ind_v_rc_min,j,i)
                   IF ( vegetation_pars_f%pars_xy(ind_v_rc_lai,j,i) /=         &
                        vegetation_pars_f%fill )                               &
                      surf_lsm_v(l)%lai(m)           =                         &
                                   vegetation_pars_f%pars_xy(ind_v_rc_lai,j,i)
                   IF ( vegetation_pars_f%pars_xy(ind_v_c_veg,j,i) /=          &
                        vegetation_pars_f%fill )                               &
                      surf_lsm_v(l)%c_veg(m)         =                         &
                                   vegetation_pars_f%pars_xy(ind_v_c_veg,j,i)
                   IF ( vegetation_pars_f%pars_xy(ind_v_gd,j,i) /=             &
                        vegetation_pars_f%fill )                               &
                     surf_lsm_v(l)%g_d(m)            =                         &
                                   vegetation_pars_f%pars_xy(ind_v_gd,j,i)
                   IF ( vegetation_pars_f%pars_xy(ind_v_z0,j,i) /=             &
                        vegetation_pars_f%fill )                               &
                      surf_lsm_v(l)%z0(m)            =                         &
                                   vegetation_pars_f%pars_xy(ind_v_z0,j,i)
                   IF ( vegetation_pars_f%pars_xy(ind_v_z0qh,j,i) /=           &
                        vegetation_pars_f%fill )  THEN
                      surf_lsm_v(l)%z0h(m)           =                         &
                                   vegetation_pars_f%pars_xy(ind_v_z0qh,j,i)
                      surf_lsm_v(l)%z0q(m)           =                         &
                                   vegetation_pars_f%pars_xy(ind_v_z0qh,j,i)
                   ENDIF
                   IF ( vegetation_pars_f%pars_xy(ind_v_lambda_s,j,i) /=       &
                        vegetation_pars_f%fill )                               &
                      surf_lsm_v(l)%lambda_surface_s(m)  =                     &
                                   vegetation_pars_f%pars_xy(ind_v_lambda_s,j,i)
                   IF ( vegetation_pars_f%pars_xy(ind_v_lambda_u,j,i) /=       &
                        vegetation_pars_f%fill )                               &
                      surf_lsm_v(l)%lambda_surface_u(m)  =                     &
                                   vegetation_pars_f%pars_xy(ind_v_lambda_u,j,i)
                   IF ( vegetation_pars_f%pars_xy(ind_v_f_sw_in,j,i) /=        &
                        vegetation_pars_f%fill )                               &
                      surf_lsm_v(l)%f_sw_in(m)           =                     &
                                   vegetation_pars_f%pars_xy(ind_v_f_sw_in,j,i)
                   IF ( vegetation_pars_f%pars_xy(ind_v_c_surf,j,i) /=         &
                        vegetation_pars_f%fill )                               &
                      surf_lsm_v(l)%c_surface(m)         =                     &
                                   vegetation_pars_f%pars_xy(ind_v_c_surf,j,i)
                   IF ( vegetation_pars_f%pars_xy(ind_v_at,j,i) /=             &
                        vegetation_pars_f%fill )                               &
                      surf_lsm_v(l)%albedo_type(ind_veg_wall,m) =              &
                                   INT( vegetation_pars_f%pars_xy(ind_v_at,j,i) )
                   IF ( vegetation_pars_f%pars_xy(ind_v_emis,j,i) /=           &
                        vegetation_pars_f%fill )                               &
                      surf_lsm_v(l)%emissivity(ind_veg_wall,m)  =              &
                                   vegetation_pars_f%pars_xy(ind_v_emis,j,i)
                ENDIF

             ENDDO
          ENDDO
       ENDIF 

!
!--    Level 1, initialization of water parameters. A horizontally 
!--    homogeneous distribution is assumed here. 
       IF ( water_type /= 0 )  THEN

          IF ( water_temperature == 9999999.9_wp )  THEN
             water_temperature = water_pars(ind_w_temp,water_type)       
          ENDIF

          IF ( z0_water == 9999999.9_wp )  THEN
             z0_water = water_pars(ind_w_z0,water_type)       
          ENDIF       

          IF ( z0h_water == 9999999.9_wp )  THEN
             z0h_water = water_pars(ind_w_z0h,water_type)       
          ENDIF  

          IF ( albedo_type == 9999999  .AND.  albedo == 9999999.9_wp )  THEN
             albedo_type = INT(water_pars(ind_w_at,water_type))        
          ENDIF
   
          IF ( emissivity == 9999999.9_wp )  THEN
             emissivity = water_pars(ind_w_emis,water_type)        
          ENDIF

       ENDIF 
!
!--    Map values onto horizontal elemements
       DO  m = 1, surf_lsm_h%ns
          IF ( surf_lsm_h%water_surface(m) )  THEN
             IF ( TRIM( initializing_actions ) /= 'read_restart_data' )        &
                t_soil_h%var_2d(:,m)           = water_temperature
             surf_lsm_h%z0(m)               = z0_water
             surf_lsm_h%z0h(m)              = z0h_water
             surf_lsm_h%z0q(m)              = z0h_water
             surf_lsm_h%lambda_surface_s(m) = 1.0E10_wp
             surf_lsm_h%lambda_surface_u(m) = 1.0E10_wp               
             surf_lsm_h%c_surface(m)        = 0.0_wp
             surf_lsm_h%albedo_type(ind_wat_win,m) = albedo_type
             surf_lsm_h%emissivity(ind_wat_win,m)  = emissivity
          ENDIF
       ENDDO
!
!--    Map values onto vertical elements, even though this does not make
!--    much sense.
       DO  l = 0, 3
          DO  m = 1, surf_lsm_v(l)%ns
             IF ( surf_lsm_v(l)%water_surface(m) )  THEN
                IF ( TRIM( initializing_actions ) /= 'read_restart_data' )     &
                   t_soil_v(l)%var_2d(:,m)           = water_temperature
                surf_lsm_v(l)%z0(m)               = z0_water
                surf_lsm_v(l)%z0h(m)              = z0h_water
                surf_lsm_v(l)%z0q(m)              = z0h_water
                surf_lsm_v(l)%lambda_surface_s(m) = 1.0E10_wp
                surf_lsm_v(l)%lambda_surface_u(m) = 1.0E10_wp               
                surf_lsm_v(l)%c_surface(m)        = 0.0_wp
                surf_lsm_v(l)%albedo_type(ind_wat_win,m) = albedo_type
                surf_lsm_v(l)%emissivity(ind_wat_win,m)  = emissivity
             ENDIF 
          ENDDO
       ENDDO
!
!
!--    Level 2, initialization of water parameters via water_type read
!--    from file. Water surfaces are initialized for each (y,x)-grid point 
!--    individually using default paramter settings according to the given
!--    water type.
!--    Note, parameter 3/4 of water_pars are albedo and emissivity,
!--    whereas paramter 3/4 of water_pars_f are heat conductivities!
       IF ( water_type_f%from_file )  THEN
!
!--       Horizontal surfaces
          DO  m = 1, surf_lsm_h%ns
             i = surf_lsm_h%i(m)
             j = surf_lsm_h%j(m)
             
             st = water_type_f%var(j,i)
             IF ( st /= water_type_f%fill  .AND.  st /= 0 )  THEN
                IF ( TRIM( initializing_actions ) /= 'read_restart_data' )     &
                   t_soil_h%var_2d(:,m) = water_pars(ind_w_temp,st)
                surf_lsm_h%z0(m)     = water_pars(ind_w_z0,st)
                surf_lsm_h%z0h(m)    = water_pars(ind_w_z0h,st)
                surf_lsm_h%z0q(m)    = water_pars(ind_w_z0h,st)
                surf_lsm_h%lambda_surface_s(m) = water_pars(ind_w_lambda_s,st)
                surf_lsm_h%lambda_surface_u(m) = water_pars(ind_w_lambda_u,st)              
                surf_lsm_h%c_surface(m)        = 0.0_wp
                surf_lsm_h%albedo_type(ind_wat_win,m) = INT( water_pars(ind_w_at,st) )
                surf_lsm_h%emissivity(ind_wat_win,m)  = water_pars(ind_w_emis,st)
             ENDIF
          ENDDO
!
!--       Vertical surfaces
          DO  l = 0, 3
             DO  m = 1, surf_lsm_v(l)%ns
                i = surf_lsm_v(l)%i(m) + MERGE( 0, surf_lsm_v(l)%ioff,         &
                                                surf_lsm_v(l)%building_covered(m) ) 
                j = surf_lsm_v(l)%j(m) + MERGE( 0, surf_lsm_v(l)%joff,         &
                                                surf_lsm_v(l)%building_covered(m) ) 
             
                st = water_type_f%var(j,i)
                IF ( st /= water_type_f%fill  .AND.  st /= 0 )  THEN
                   IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  &
                      t_soil_v(l)%var_2d(:,m) = water_pars(ind_w_temp,st)
                   surf_lsm_v(l)%z0(m)     = water_pars(ind_w_z0,st)
                   surf_lsm_v(l)%z0h(m)    = water_pars(ind_w_z0h,st)
                   surf_lsm_v(l)%z0q(m)    = water_pars(ind_w_z0h,st)
                   surf_lsm_v(l)%lambda_surface_s(m) =                         &
                                                   water_pars(ind_w_lambda_s,st)
                   surf_lsm_v(l)%lambda_surface_u(m) =                         &
                                                   water_pars(ind_w_lambda_u,st)           
                   surf_lsm_v(l)%c_surface(m)     = 0.0_wp
                   surf_lsm_v(l)%albedo_type(ind_wat_win,m) =                  &
                                                  INT( water_pars(ind_w_at,st) )
                   surf_lsm_v(l)%emissivity(ind_wat_win,m)  =                  &
                                                  water_pars(ind_w_emis,st)
                ENDIF
             ENDDO
          ENDDO
       ENDIF      

!
!--    Level 3, initialization of water parameters at single (x,y) 
!--    position via water_pars read from file.
       IF ( water_pars_f%from_file )  THEN
!
!--       Horizontal surfaces 
          DO  m = 1, surf_lsm_h%ns
             i = surf_lsm_h%i(m)
             j = surf_lsm_h%j(m)
!
!--          If surface element is not a water surface and any value in 
!--          water_pars is given, neglect this information and give an
!--          informative message that this value will not be used.   
             IF ( .NOT. surf_lsm_h%water_surface(m)  .AND.                     &
                   ANY( water_pars_f%pars_xy(:,j,i) /= water_pars_f%fill ) )  THEN
                WRITE( message_string, * )                                     &
                              'surface element at grid point (j,i) = (',       &
                              j, i, ') is not a water surface, ',              &
                              'so that information given in ',                 &
                              'water_pars at this point is neglected.' 
                CALL message( 'land_surface_model_mod', 'PA0999', 0, 0, 0, 6, 0 )
             ELSE
                IF ( water_pars_f%pars_xy(ind_w_temp,j,i) /=                   &
                     water_pars_f%fill  .AND.                                  &
                     TRIM( initializing_actions ) /= 'read_restart_data' )     &
                      t_soil_h%var_2d(:,m) = water_pars_f%pars_xy(ind_w_temp,j,i)

                IF ( water_pars_f%pars_xy(ind_w_z0,j,i) /= water_pars_f%fill ) &
                   surf_lsm_h%z0(m)     = water_pars_f%pars_xy(ind_w_z0,j,i)

                IF ( water_pars_f%pars_xy(ind_w_z0h,j,i) /= water_pars_f%fill )&
                THEN
                   surf_lsm_h%z0h(m)    = water_pars_f%pars_xy(ind_w_z0h,j,i)
                   surf_lsm_h%z0q(m)    = water_pars_f%pars_xy(ind_w_z0h,j,i)
                ENDIF
                IF ( water_pars_f%pars_xy(ind_w_lambda_s,j,i) /=               &
                     water_pars_f%fill )                                       &
                   surf_lsm_h%lambda_surface_s(m) =                            &
                                        water_pars_f%pars_xy(ind_w_lambda_s,j,i)

                IF ( water_pars_f%pars_xy(ind_w_lambda_u,j,i) /=               &
                      water_pars_f%fill )                                      &
                   surf_lsm_h%lambda_surface_u(m) =                            &
                                        water_pars_f%pars_xy(ind_w_lambda_u,j,i)     
       
                IF ( water_pars_f%pars_xy(ind_w_at,j,i) /=                     &
                     water_pars_f%fill )                                       &
                   surf_lsm_h%albedo_type(ind_wat_win,m) =                     &
                                       INT( water_pars_f%pars_xy(ind_w_at,j,i) )

                IF ( water_pars_f%pars_xy(ind_w_emis,j,i) /=                   &
                     water_pars_f%fill )                                       &
                   surf_lsm_h%emissivity(ind_wat_win,m) =                      &
                   water_pars_f%pars_xy(ind_w_emis,j,i)  
             ENDIF
          ENDDO
!
!--       Vertical surfaces 
          DO  l = 0, 3
             DO  m = 1, surf_lsm_v(l)%ns
                i = surf_lsm_v(l)%i(m) + MERGE( 0, surf_lsm_v(l)%ioff,         &
                                                surf_lsm_v(l)%building_covered(m) ) 
                j = surf_lsm_v(l)%j(m) + MERGE( 0, surf_lsm_v(l)%joff,         &
                                                surf_lsm_v(l)%building_covered(m) ) 
!
!--             If surface element is not a water surface and any value in 
!--             water_pars is given, neglect this information and give an
!--             informative message that this value will not be used.   
                IF ( .NOT. surf_lsm_v(l)%water_surface(m)  .AND.               &
                      ANY( water_pars_f%pars_xy(:,j,i) /=                      &
                      water_pars_f%fill ) )  THEN
                   WRITE( message_string, * )                                  &
                              'surface element at grid point (j,i) = (',       &
                              j, i, ') is not a water surface, ',              &
                              'so that information given in ',                 &
                              'water_pars at this point is neglected.' 
                   CALL message( 'land_surface_model_mod', 'PA0999',           &
                                  0, 0, 0, 6, 0 )
                ELSE

                   IF ( water_pars_f%pars_xy(ind_w_temp,j,i) /=                &
                     water_pars_f%fill  .AND.                                  &
                     TRIM( initializing_actions ) /= 'read_restart_data' )     &
                      t_soil_v(l)%var_2d(:,m) = water_pars_f%pars_xy(ind_w_temp,j,i)

                   IF ( water_pars_f%pars_xy(ind_w_z0,j,i) /=                  &
                        water_pars_f%fill )                                    &
                      surf_lsm_v(l)%z0(m)   = water_pars_f%pars_xy(ind_w_z0,j,i)

                   IF ( water_pars_f%pars_xy(ind_w_z0h,j,i) /=                 &
                       water_pars_f%fill )  THEN
                      surf_lsm_v(l)%z0h(m)  = water_pars_f%pars_xy(ind_w_z0h,j,i)
                      surf_lsm_v(l)%z0q(m)  = water_pars_f%pars_xy(ind_w_z0h,j,i)
                   ENDIF

                   IF ( water_pars_f%pars_xy(ind_w_lambda_s,j,i) /=            &
                        water_pars_f%fill )                                    &
                      surf_lsm_v(l)%lambda_surface_s(m) =                      &
                                      water_pars_f%pars_xy(ind_w_lambda_s,j,i)

                   IF ( water_pars_f%pars_xy(ind_w_lambda_u,j,i) /=            &
                        water_pars_f%fill )                                    &
                      surf_lsm_v(l)%lambda_surface_u(m) =                      &
                                      water_pars_f%pars_xy(ind_w_lambda_u,j,i)    
 
                   IF ( water_pars_f%pars_xy(ind_w_at,j,i) /=                  &
                        water_pars_f%fill )                                    &
                      surf_lsm_v(l)%albedo_type(ind_wat_win,m) =               &
                                      INT( water_pars_f%pars_xy(ind_w_at,j,i) )

                   IF ( water_pars_f%pars_xy(ind_w_emis,j,i) /=                &
                        water_pars_f%fill )                                    &
                      surf_lsm_v(l)%emissivity(ind_wat_win,m)  =               &
                                      water_pars_f%pars_xy(ind_w_emis,j,i)  
                ENDIF
             ENDDO
          ENDDO

       ENDIF
!
!--    Initialize pavement-type surfaces, level 1
       IF ( pavement_type /= 0 )  THEN  

!
!--       When a pavement_type is used, overwrite a possible setting of 
!--       the pavement depth as it is already defined by the pavement type
          pavement_depth_level = 0

          IF ( z0_pavement == 9999999.9_wp )  THEN
             z0_pavement  = pavement_pars(ind_p_z0,pavement_type) 
          ENDIF

          IF ( z0h_pavement == 9999999.9_wp )  THEN
             z0h_pavement = pavement_pars(ind_p_z0h,pavement_type)
          ENDIF

          IF ( pavement_heat_conduct == 9999999.9_wp )  THEN
             pavement_heat_conduct = pavement_subsurface_pars_1(0,pavement_type)
          ENDIF

          IF ( pavement_heat_capacity == 9999999.9_wp )  THEN
             pavement_heat_capacity = pavement_subsurface_pars_2(0,pavement_type)
          ENDIF   
    
          IF ( albedo_type == 9999999  .AND.  albedo == 9999999.9_wp )  THEN
             albedo_type = INT(pavement_pars(ind_p_at,pavement_type))        
          ENDIF
   
          IF ( emissivity == 9999999.9_wp )  THEN
             emissivity = pavement_pars(ind_p_emis,pavement_type)        
          ENDIF

!
!--       If the depth level of the pavement is not set, determine it from
!--       lookup table.
          IF ( pavement_depth_level == 0 )  THEN
             DO  k = nzb_soil, nzt_soil  
                IF ( pavement_subsurface_pars_1(k,pavement_type) == 9999999.9_wp &
                .OR. pavement_subsurface_pars_2(k,pavement_type) == 9999999.9_wp)&
                THEN
                   nzt_pavement = k-1
                   EXIT
                ENDIF
             ENDDO
          ELSE
             nzt_pavement = pavement_depth_level
          ENDIF

       ENDIF
!
!--    Level 1 initialization of pavement type surfaces. Horizontally 
!--    homogeneous characteristics are assumed
       surf_lsm_h%nzt_pavement = pavement_depth_level
       DO  m = 1, surf_lsm_h%ns
          IF ( surf_lsm_h%pavement_surface(m) )  THEN
             surf_lsm_h%nzt_pavement(m)        = nzt_pavement
             surf_lsm_h%z0(m)                  = z0_pavement
             surf_lsm_h%z0h(m)                 = z0h_pavement
             surf_lsm_h%z0q(m)                 = z0h_pavement
             surf_lsm_h%lambda_surface_s(m)    = pavement_heat_conduct         &
                                                  * ddz_soil(nzb_soil)         &
                                                  * 2.0_wp    
             surf_lsm_h%lambda_surface_u(m)    = pavement_heat_conduct         &
                                                  * ddz_soil(nzb_soil)         &
                                                  * 2.0_wp           
             surf_lsm_h%c_surface(m)           = pavement_heat_capacity        &
                                                        * dz_soil(nzb_soil)    &
                                                        * 0.25_wp                                   

             surf_lsm_h%albedo_type(ind_pav_green,m) = albedo_type
             surf_lsm_h%emissivity(ind_pav_green,m)  = emissivity      
      
             IF ( pavement_type /= 0 )  THEN
                DO  k = nzb_soil, surf_lsm_h%nzt_pavement(m)
                   surf_lsm_h%lambda_h_def(k,m)    =                           &
                                     pavement_subsurface_pars_1(k,pavement_type)                       
                   surf_lsm_h%rho_c_total_def(k,m) =                           &
                                     pavement_subsurface_pars_2(k,pavement_type) 
                ENDDO
             ELSE
                surf_lsm_v(l)%lambda_h_def(:,m)     = pavement_heat_conduct
                surf_lsm_v(l)%rho_c_total_def(:,m)  = pavement_heat_capacity
             ENDIF        
          ENDIF
       ENDDO                               

       DO  l = 0, 3
          surf_lsm_v(l)%nzt_pavement = pavement_depth_level
          DO  m = 1, surf_lsm_v(l)%ns
             IF ( surf_lsm_v(l)%pavement_surface(m) )  THEN
                surf_lsm_v(l)%nzt_pavement(m)        = nzt_pavement
                surf_lsm_v(l)%z0(m)                  = z0_pavement
                surf_lsm_v(l)%z0h(m)                 = z0h_pavement
                surf_lsm_v(l)%z0q(m)                 = z0h_pavement
                surf_lsm_v(l)%lambda_surface_s(m)    = pavement_heat_conduct   &
                                                  * ddz_soil(nzb_soil)         &
                                                  * 2.0_wp    
                surf_lsm_v(l)%lambda_surface_u(m)    = pavement_heat_conduct   &
                                                  * ddz_soil(nzb_soil)         &
                                                  * 2.0_wp           
                surf_lsm_v(l)%c_surface(m)           = pavement_heat_capacity  &
                                                        * dz_soil(nzb_soil)    &
                                                        * 0.25_wp                                      

                surf_lsm_v(l)%albedo_type(ind_pav_green,m) = albedo_type
                surf_lsm_v(l)%emissivity(ind_pav_green,m)  = emissivity

                IF ( pavement_type /= 0 )  THEN
                   DO  k = nzb_soil, surf_lsm_v(l)%nzt_pavement(m)
                      surf_lsm_v(l)%lambda_h_def(k,m)    =                     &
                                     pavement_subsurface_pars_1(k,pavement_type)                       
                      surf_lsm_v(l)%rho_c_total_def(k,m) =                     &
                                     pavement_subsurface_pars_2(k,pavement_type) 
                   ENDDO
                ELSE
                   surf_lsm_v(l)%lambda_h_def(:,m)     = pavement_heat_conduct
                   surf_lsm_v(l)%rho_c_total_def(:,m)  = pavement_heat_capacity
                ENDIF     
             ENDIF
          ENDDO
       ENDDO
!
!--    Level 2 initialization of pavement type surfaces via pavement_type read
!--    from file. Pavement surfaces are initialized for each (y,x)-grid point 
!--    individually.
       IF ( pavement_type_f%from_file )  THEN
!
!--       Horizontal surfaces
          DO  m = 1, surf_lsm_h%ns
             i = surf_lsm_h%i(m)
             j = surf_lsm_h%j(m)
             
             st = pavement_type_f%var(j,i)
             IF ( st /= pavement_type_f%fill  .AND.  st /= 0 )  THEN
!
!--             Determine deepmost index of pavement layer
                DO  k = nzb_soil, nzt_soil  
                   IF ( pavement_subsurface_pars_1(k,st) == 9999999.9_wp       &
                   .OR. pavement_subsurface_pars_2(k,st) == 9999999.9_wp)      &
                   THEN
                      surf_lsm_h%nzt_pavement(m) = k-1
                      EXIT
                   ENDIF
                ENDDO

                surf_lsm_h%z0(m)                = pavement_pars(ind_p_z0,st)
                surf_lsm_h%z0h(m)               = pavement_pars(ind_p_z0h,st)
                surf_lsm_h%z0q(m)               = pavement_pars(ind_p_z0h,st)

                surf_lsm_h%lambda_surface_s(m)  =                              &
                                              pavement_subsurface_pars_1(0,st) &
                                                  * ddz_soil(nzb_soil)         &
                                                  * 2.0_wp    
                surf_lsm_h%lambda_surface_u(m)  =                              &
                                              pavement_subsurface_pars_1(0,st) &
                                                  * ddz_soil(nzb_soil)         &
                                                  * 2.0_wp        
                surf_lsm_h%c_surface(m)         =                              &
                                               pavement_subsurface_pars_2(0,st)&
                                                        * dz_soil(nzb_soil)    &
                                                        * 0.25_wp                               
                surf_lsm_h%albedo_type(ind_pav_green,m) = INT( pavement_pars(ind_p_at,st) )
                surf_lsm_h%emissivity(ind_pav_green,m)  = pavement_pars(ind_p_emis,st)  

                DO  k = nzb_soil, surf_lsm_h%nzt_pavement(m)
                   surf_lsm_h%lambda_h_def(k,m)    =                           &
                                     pavement_subsurface_pars_1(k,pavement_type)                       
                   surf_lsm_h%rho_c_total_def(k,m) =                           &
                                     pavement_subsurface_pars_2(k,pavement_type) 
                ENDDO   
             ENDIF
          ENDDO
!
!--       Vertical surfaces
          DO  l = 0, 3
             DO  m = 1, surf_lsm_v(l)%ns
                i = surf_lsm_v(l)%i(m) + MERGE( 0, surf_lsm_v(l)%ioff,         &
                                                surf_lsm_v(l)%building_covered(m) ) 
                j = surf_lsm_v(l)%j(m) + MERGE( 0, surf_lsm_v(l)%joff,         &
                                                surf_lsm_v(l)%building_covered(m) ) 
             
                st = pavement_type_f%var(j,i)
                IF ( st /= pavement_type_f%fill  .AND.  st /= 0 )  THEN
!
!--                Determine deepmost index of pavement layer
                   DO  k = nzb_soil, nzt_soil  
                      IF ( pavement_subsurface_pars_1(k,st) == 9999999.9_wp    &
                      .OR. pavement_subsurface_pars_2(k,st) == 9999999.9_wp)   &
                      THEN
                         surf_lsm_v(l)%nzt_pavement(m) = k-1
                         EXIT
                      ENDIF
                   ENDDO

                   surf_lsm_v(l)%z0(m)  = pavement_pars(ind_p_z0,st)
                   surf_lsm_v(l)%z0h(m) = pavement_pars(ind_p_z0h,st)
                   surf_lsm_v(l)%z0q(m) = pavement_pars(ind_p_z0h,st)

                   surf_lsm_v(l)%lambda_surface_s(m)  =                        &
                                              pavement_subsurface_pars_1(0,st) &
                                                  * ddz_soil(nzb_soil)         & 
                                                  * 2.0_wp    
                   surf_lsm_v(l)%lambda_surface_u(m)  =                        &
                                              pavement_subsurface_pars_1(0,st) &
                                                  * ddz_soil(nzb_soil)         &
                                                  * 2.0_wp     

                   surf_lsm_v(l)%c_surface(m)    =                             &
                                           pavement_subsurface_pars_2(0,st)    &
                                                        * dz_soil(nzb_soil)    &
                                                        * 0.25_wp                                    
                   surf_lsm_v(l)%albedo_type(ind_pav_green,m) =                &
                                              INT( pavement_pars(ind_p_at,st) )
                   surf_lsm_v(l)%emissivity(ind_pav_green,m)  =                &
                                              pavement_pars(ind_p_emis,st)    
                                              
                   DO  k = nzb_soil, surf_lsm_v(l)%nzt_pavement(m)
                      surf_lsm_v(l)%lambda_h_def(k,m)    =                     &
                                    pavement_subsurface_pars_1(k,pavement_type)                       
                      surf_lsm_v(l)%rho_c_total_def(k,m) =                     &
                                    pavement_subsurface_pars_2(k,pavement_type) 
                   ENDDO   
                ENDIF
             ENDDO
          ENDDO
       ENDIF 
!
!--    Level 3, initialization of pavement parameters at single (x,y) 
!--    position via pavement_pars read from file.
       IF ( pavement_pars_f%from_file )  THEN
!
!--       Horizontal surfaces 
          DO  m = 1, surf_lsm_h%ns
             i = surf_lsm_h%i(m)
             j = surf_lsm_h%j(m)
!
!--          If surface element is not a pavement surface and any value in 
!--          pavement_pars is given, neglect this information and give an
!--          informative message that this value will not be used.   
             IF ( .NOT. surf_lsm_h%pavement_surface(m)  .AND.                  &
                   ANY( pavement_pars_f%pars_xy(:,j,i) /=                      &
                   pavement_pars_f%fill ) )  THEN
                WRITE( message_string, * )                                     &
                              'surface element at grid point (j,i) = (',       &
                              j, i, ') is not a pavement surface, ',           &
                              'so that information given in ',                 &
                              'pavement_pars at this point is neglected.' 
                CALL message( 'land_surface_model_mod', 'PA0999', 0, 0, 0, 6, 0 )
             ELSE
                IF ( pavement_pars_f%pars_xy(ind_p_z0,j,i) /=                  &
                     pavement_pars_f%fill )                                    &
                   surf_lsm_h%z0(m)  = pavement_pars_f%pars_xy(ind_p_z0,j,i)
                IF ( pavement_pars_f%pars_xy(ind_p_z0h,j,i) /=                 &
                     pavement_pars_f%fill )  THEN
                   surf_lsm_h%z0h(m) = pavement_pars_f%pars_xy(ind_p_z0h,j,i)
                   surf_lsm_h%z0q(m) = pavement_pars_f%pars_xy(ind_p_z0h,j,i)
                ENDIF
                IF ( pavement_subsurface_pars_f%pars_xyz(ind_p_lambda_h,0,j,i) &
                     /= pavement_subsurface_pars_f%fill )  THEN
                   surf_lsm_h%lambda_surface_s(m)  =                           &
                      pavement_subsurface_pars_f%pars_xyz(ind_p_lambda_h,0,j,i)&
                                                  * ddz_soil(nzb_soil)         &
                                                  * 2.0_wp
                   surf_lsm_h%lambda_surface_u(m)  =                           &
                      pavement_subsurface_pars_f%pars_xyz(ind_p_lambda_h,0,j,i)&
                                                  * ddz_soil(nzb_soil)         &
                                                  * 2.0_wp    
                ENDIF
                IF ( pavement_subsurface_pars_f%pars_xyz(ind_p_rho_c,0,j,i) /= &
                     pavement_subsurface_pars_f%fill )  THEN
                   surf_lsm_h%c_surface(m)     =                               &
                      pavement_subsurface_pars_f%pars_xyz(ind_p_rho_c,0,j,i)   &
                                                  * dz_soil(nzb_soil)          &
                                                  * 0.25_wp                                    
                ENDIF
                IF ( pavement_pars_f%pars_xy(ind_p_at,j,i) /=                  &
                     pavement_pars_f%fill )                                    &
                   surf_lsm_h%albedo_type(ind_pav_green,m) =                   &
                                              INT( pavement_pars(ind_p_at,st) )
                IF ( pavement_pars_f%pars_xy(ind_p_emis,j,i) /=                &
                     pavement_pars_f%fill )                                    &
                   surf_lsm_h%emissivity(ind_pav_green,m)  =                   &
                                              pavement_pars(ind_p_emis,st) 
             ENDIF 

          ENDDO
!
!--       Vertical surfaces 
          DO  l = 0, 3
             DO  m = 1, surf_lsm_v(l)%ns
                i = surf_lsm_v(l)%i(m) + MERGE( 0, surf_lsm_v(l)%ioff,         &
                                                surf_lsm_v(l)%building_covered(m) ) 
                j = surf_lsm_v(l)%j(m) + MERGE( 0, surf_lsm_v(l)%joff,         &
                                                surf_lsm_v(l)%building_covered(m) ) 
!
!--             If surface element is not a pavement surface and any value in 
!--             pavement_pars is given, neglect this information and give an
!--             informative message that this value will not be used.   
                IF ( .NOT. surf_lsm_v(l)%pavement_surface(m)  .AND.            &
                      ANY( pavement_pars_f%pars_xy(:,j,i) /=                   &
                      pavement_pars_f%fill ) )  THEN
                   WRITE( message_string, * )                                  &
                                 'surface element at grid point (j,i) = (',    &
                                 j, i, ') is not a pavement surface, ',        &
                                 'so that information given in ',              &
                                 'pavement_pars at this point is neglected.' 
                   CALL message( 'land_surface_model_mod', 'PA0999', 0, 0, 0, 6, 0 )
                ELSE

                   IF ( pavement_pars_f%pars_xy(ind_p_z0,j,i) /=               &
                        pavement_pars_f%fill )                                 &
                      surf_lsm_v(l)%z0(m) = pavement_pars_f%pars_xy(ind_p_z0,j,i)
                   IF ( pavement_pars_f%pars_xy(ind_p_z0h,j,i) /=              &
                        pavement_pars_f%fill )  THEN
                      surf_lsm_v(l)%z0h(m) = pavement_pars_f%pars_xy(ind_p_z0h,j,i)
                      surf_lsm_v(l)%z0q(m) = pavement_pars_f%pars_xy(ind_p_z0h,j,i)
                   ENDIF
                   IF ( pavement_subsurface_pars_f%pars_xyz(ind_p_lambda_h,0,j,i)&
                        /= pavement_subsurface_pars_f%fill )  THEN
                      surf_lsm_v(l)%lambda_surface_s(m) =                      &
                      pavement_subsurface_pars_f%pars_xyz(ind_p_lambda_h,0,j,i)&
                                                  * ddz_soil(nzb_soil)         &
                                                  * 2.0_wp
                      surf_lsm_v(l)%lambda_surface_u(m) =                      &
                      pavement_subsurface_pars_f%pars_xyz(ind_p_lambda_h,0,j,i)&
                                                  * ddz_soil(nzb_soil)         &
                                                  * 2.0_wp    
                   ENDIF
                   IF ( pavement_subsurface_pars_f%pars_xyz(ind_p_rho_c,0,j,i) &
                        /= pavement_subsurface_pars_f%fill )  THEN
                      surf_lsm_v(l)%c_surface(m)    =                          &
                         pavement_subsurface_pars_f%pars_xyz(ind_p_rho_c,0,j,i)&
                                                  * dz_soil(nzb_soil)          &
                                                  * 0.25_wp                                  
                   ENDIF
                   IF ( pavement_pars_f%pars_xy(ind_p_at,j,i) /=               &
                        pavement_pars_f%fill )                                 &
                      surf_lsm_v(l)%albedo_type(ind_pav_green,m) =             &
                                            INT( pavement_pars(ind_p_at,st) )

                   IF ( pavement_pars_f%pars_xy(ind_p_emis,j,i) /=             &
                        pavement_pars_f%fill )                                 &
                      surf_lsm_v(l)%emissivity(ind_pav_green,m)  =             &
                                            pavement_pars(ind_p_emis,st)  
                ENDIF 
             ENDDO
          ENDDO
       ENDIF
!
!--    Moreover, for grid points which are flagged with pavement-type 0 or whre 
!--    pavement_subsurface_pars_f is provided, soil heat conductivity and 
!--    capacity are initialized with parameters given in        
!--    pavement_subsurface_pars read from file. 
       IF ( pavement_subsurface_pars_f%from_file )  THEN
!
!--       Set pavement depth to nzt_soil. Please note, this is just a 
!--       workaround at the moment. 
          DO  m = 1, surf_lsm_h%ns
             IF ( surf_lsm_h%pavement_surface(m) )  THEN

                i = surf_lsm_h%i(m)
                j = surf_lsm_h%j(m)

                surf_lsm_h%nzt_pavement(m) = nzt_soil

                DO  k = nzb_soil, nzt_soil  
                   surf_lsm_h%lambda_h_def(k,m) =                              &
                       pavement_subsurface_pars_f%pars_xyz(ind_p_lambda_h,k,j,i)
                   surf_lsm_h%rho_c_total_def(k,m) =                           &
                       pavement_subsurface_pars_f%pars_xyz(ind_p_rho_c,k,j,i)
                ENDDO

             ENDIF
          ENDDO
          DO  l = 0, 3
             DO  m = 1, surf_lsm_v(l)%ns
                IF ( surf_lsm_v(l)%pavement_surface(m) )  THEN

                   i = surf_lsm_v(l)%i(m) + MERGE( 0, surf_lsm_v(l)%ioff,      &
                                                surf_lsm_v(l)%building_covered(m) ) 
                   j = surf_lsm_v(l)%j(m) + MERGE( 0, surf_lsm_v(l)%joff,      &
                                                surf_lsm_v(l)%building_covered(m) ) 

                   surf_lsm_v(l)%nzt_pavement(m) = nzt_soil

                   DO  k = nzb_soil, nzt_soil  
                      surf_lsm_v(l)%lambda_h_def(k,m) =                        &
                       pavement_subsurface_pars_f%pars_xyz(ind_p_lambda_h,k,j,i)
                      surf_lsm_v(l)%rho_c_total_def(k,m) =                     &
                       pavement_subsurface_pars_f%pars_xyz(ind_p_rho_c,k,j,i)
                   ENDDO

                ENDIF
             ENDDO
          ENDDO
       ENDIF

!
!--    Initial run actions
       IF (  TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
!
!--       In case of nested runs, check if soil is initialized via dynamic
!--       input file in root domain and distribute this information
!--       to all embedded child domains. This case, soil moisture and
!--       temperature will be initialized via root domain.  
          init_soil_dynamically_in_child = .FALSE. 
          IF ( nested_run )  THEN
#if defined( __parallel )
             CALL MPI_ALLREDUCE( init_3d%from_file_tsoil  .OR.                 &
                                 init_3d%from_file_msoil,                      &
                                 init_soil_dynamically_in_child,               &
                                 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr )
#endif
          ENDIF

!
!--       First, initialize soil temperature and moisture. 
!--       According to the initialization for surface and soil parameters, 
!--       initialize soil moisture and temperature via a level approach. This 
!--       is to assure that all surface elements are initialized, even if 
!--       data provided from input file contains fill values at some locations.
!--       Level 1, initialization via profiles given in parameter file
          DO  m = 1, surf_lsm_h%ns
             IF ( surf_lsm_h%vegetation_surface(m)  .OR.                       &
                  surf_lsm_h%pavement_surface(m) )  THEN
                DO  k = nzb_soil, nzt_soil  
                   t_soil_h%var_2d(k,m) = soil_temperature(k)
                   m_soil_h%var_2d(k,m) = soil_moisture(k)
                ENDDO
                t_soil_h%var_2d(nzt_soil+1,m) = deep_soil_temperature
             ENDIF
          ENDDO
          DO  l = 0, 3
             DO  m = 1, surf_lsm_v(l)%ns
                IF ( surf_lsm_v(l)%vegetation_surface(m)  .OR.                 &
                     surf_lsm_v(l)%pavement_surface(m) )  THEN
                   DO  k = nzb_soil, nzt_soil  
                      t_soil_v(l)%var_2d(k,m) = soil_temperature(k)
                      m_soil_v(l)%var_2d(k,m) = soil_moisture(k)
                   ENDDO
                   t_soil_v(l)%var_2d(nzt_soil+1,m) = deep_soil_temperature
                ENDIF
             ENDDO
          ENDDO
!
!--       Initialization of soil moisture and temperature from file. 
!--       In case of nested runs, only root parent reads dynamic input file. 
!--       This case, transfer respective soil information provide
!--       by dynamic input file from root parent domain onto all other domains.
          IF ( init_soil_dynamically_in_child )  THEN
!
!--          Child domains will be only initialized with horizontally 
!--          averaged soil profiles in parent domain (for sake of simplicity). 
!--          If required, average soil data on root parent domain before 
!--          distribute onto child domains.
             IF ( init_3d%from_file_msoil  .AND.  init_3d%lod_msoil == 2 )     &
             THEN
                ALLOCATE( pr_soil_init(0:init_3d%nzs-1) )

                DO  k = 0, init_3d%nzs-1
                   pr_soil_init(k) = SUM( init_3d%msoil(k,nys:nyn,nxl:nxr)  )
                ENDDO
!
!--             Allocate 1D array for soil-moisture profile (will not be 
!--             allocated in lod==2 case). 
                ALLOCATE( init_3d%msoil_init(0:init_3d%nzs-1) )
                init_3d%msoil_init = 0.0_wp
#if defined( __parallel )
                CALL MPI_ALLREDUCE( pr_soil_init(0), init_3d%msoil_init(0),    &
                                    SIZE(pr_soil_init),                        &
                                    MPI_REAL, MPI_SUM, comm2d, ierr )
#endif
                init_3d%msoil_init = init_3d%msoil_init /                      &
                                        REAL( ( nx + 1 ) * ( ny + 1), KIND=wp )
                DEALLOCATE( pr_soil_init )
             ENDIF
             IF ( init_3d%from_file_tsoil  .AND.  init_3d%lod_tsoil == 2 )  THEN
                ALLOCATE( pr_soil_init(0:init_3d%nzs-1) )

                DO  k = 0, init_3d%nzs-1
                   pr_soil_init(k) = SUM( init_3d%tsoil(k,nys:nyn,nxl:nxr)  )
                ENDDO
!
!--             Allocate 1D array for soil-temperature profile (will not be 
!--             allocated in lod==2 case). 
                ALLOCATE( init_3d%tsoil_init(0:init_3d%nzs-1) )
                init_3d%tsoil_init = 0.0_wp
#if defined( __parallel )
                CALL MPI_ALLREDUCE( pr_soil_init(0), init_3d%tsoil_init(0),    &
                                    SIZE(pr_soil_init),                        &
                                    MPI_REAL, MPI_SUM, comm2d, ierr )
#endif
                init_3d%tsoil_init = init_3d%tsoil_init /                      &
                                        REAL( ( nx + 1 ) * ( ny + 1), KIND=wp )
                DEALLOCATE( pr_soil_init )

             ENDIF

#if defined( __parallel )
!
!--          Distribute soil grid information on file from root to all childs.
!--          Only process with rank 0 sends the information. 
             CALL MPI_BCAST( init_3d%nzs,    1,                                &
                             MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

             IF ( .NOT.  ALLOCATED( init_3d%z_soil ) )                         &
                ALLOCATE( init_3d%z_soil(1:init_3d%nzs) )

             CALL MPI_BCAST( init_3d%z_soil, SIZE(init_3d%z_soil),             &
                             MPI_REAL, 0, MPI_COMM_WORLD, ierr )
#endif
!
!--          ALLOCATE arrays on child domains and set control attributes.
             IF ( .NOT. init_3d%from_file_msoil )  THEN
                ALLOCATE( init_3d%msoil_init(0:init_3d%nzs-1) )
                init_3d%lod_msoil = 1
                init_3d%from_file_msoil = .TRUE.
             ENDIF
             IF ( .NOT. init_3d%from_file_tsoil )  THEN
                ALLOCATE( init_3d%tsoil_init(0:init_3d%nzs-1) )
                init_3d%lod_tsoil = 1
                init_3d%from_file_tsoil = .TRUE.
             ENDIF


#if defined( __parallel )
!
!--          Distribute soil profiles from root to all childs
             CALL MPI_BCAST( init_3d%msoil_init, SIZE(init_3d%msoil_init),     &
                             MPI_REAL, 0, MPI_COMM_WORLD, ierr )
             CALL MPI_BCAST( init_3d%tsoil_init, SIZE(init_3d%tsoil_init),     &
                             MPI_REAL, 0, MPI_COMM_WORLD, ierr )
#endif

          ENDIF
!
!--       Proceed with Level 2 initialization. Information from dynamic input
!--       is now available on all processes. 
          IF ( init_3d%from_file_msoil )  THEN

             IF ( init_3d%lod_msoil == 1 )  THEN
                DO  m = 1, surf_lsm_h%ns

                   CALL netcdf_data_input_interpolate(                         &
                                       m_soil_h%var_2d(nzb_soil:nzt_soil,m),   &
                                       init_3d%msoil_init(:),                  &
                                       zs(nzb_soil:nzt_soil), init_3d%z_soil,  &
                                       nzb_soil, nzt_soil,                     &
                                       nzb_soil, init_3d%nzs-1 )
                ENDDO
                DO  l = 0, 3
                   DO  m = 1, surf_lsm_v(l)%ns

                      CALL netcdf_data_input_interpolate(                      &
                                       m_soil_v(l)%var_2d(nzb_soil:nzt_soil,m),&
                                       init_3d%msoil_init(:),                  &
                                       zs(nzb_soil:nzt_soil), init_3d%z_soil,  &
                                       nzb_soil, nzt_soil,                     &
                                       nzb_soil, init_3d%nzs-1 )
                   ENDDO
                ENDDO
             ELSE

                DO  m = 1, surf_lsm_h%ns
                   i = surf_lsm_h%i(m)
                   j = surf_lsm_h%j(m)

                   IF ( init_3d%msoil(0,j,i) /= init_3d%fill_msoil )           &
                      CALL netcdf_data_input_interpolate(                      &
                                       m_soil_h%var_2d(nzb_soil:nzt_soil,m),   &
                                       init_3d%msoil(:,j,i),                   &
                                       zs(nzb_soil:nzt_soil), init_3d%z_soil,  &
                                       nzb_soil, nzt_soil,                     &
                                       nzb_soil, init_3d%nzs-1 )
                ENDDO
                DO  l = 0, 3
                   DO  m = 1, surf_lsm_v(l)%ns
                      i = surf_lsm_v(l)%i(m) + MERGE( 0, surf_lsm_v(l)%ioff,   &
                                             surf_lsm_v(l)%building_covered(m) ) 
                      j = surf_lsm_v(l)%j(m) + MERGE( 0, surf_lsm_v(l)%joff,   &
                                             surf_lsm_v(l)%building_covered(m) ) 

                      IF ( init_3d%msoil(0,j,i) /= init_3d%fill_msoil )        &
                         CALL netcdf_data_input_interpolate(                   &
                                       m_soil_v(l)%var_2d(nzb_soil:nzt_soil,m),&
                                       init_3d%msoil(:,j,i),                   &
                                       zs(nzb_soil:nzt_soil), init_3d%z_soil,  &
                                       nzb_soil, nzt_soil,                     &
                                       nzb_soil, init_3d%nzs-1 )
                   ENDDO
                ENDDO
             ENDIF
          ENDIF
!
!--       Soil temperature
          IF ( init_3d%from_file_tsoil )  THEN

             IF ( init_3d%lod_tsoil == 1 )  THEN ! change to 1 if provided correctly by INIFOR
                DO  m = 1, surf_lsm_h%ns

                   CALL netcdf_data_input_interpolate(                         &
                                       t_soil_h%var_2d(nzb_soil:nzt_soil,m),   &
                                       init_3d%tsoil_init(:),                  &
                                       zs(nzb_soil:nzt_soil), init_3d%z_soil,  &
                                       nzb_soil, nzt_soil,                     &
                                       nzb_soil, init_3d%nzs-1 )
                   t_soil_h%var_2d(nzt_soil+1,m) = t_soil_h%var_2d(nzt_soil,m)
                ENDDO
                DO  l = 0, 3
                   DO  m = 1, surf_lsm_v(l)%ns

                      CALL netcdf_data_input_interpolate(                      &
                                       t_soil_v(l)%var_2d(nzb_soil:nzt_soil,m),&
                                       init_3d%tsoil_init(:),                  &
                                       zs(nzb_soil:nzt_soil), init_3d%z_soil,  &
                                       nzb_soil, nzt_soil,                     &
                                       nzb_soil, init_3d%nzs-1 )
                      t_soil_v(l)%var_2d(nzt_soil+1,m) =                       &
                                                 t_soil_v(l)%var_2d(nzt_soil,m)
                   ENDDO
                ENDDO
             ELSE

                DO  m = 1, surf_lsm_h%ns
                   i = surf_lsm_h%i(m)
                   j = surf_lsm_h%j(m)

                   IF ( init_3d%msoil(0,j,i) /= init_3d%fill_msoil )           &
                      CALL netcdf_data_input_interpolate(                      &
                                       t_soil_h%var_2d(nzb_soil:nzt_soil,m),   &
                                       init_3d%tsoil(:,j,i),                   &
                                       zs(nzb_soil:nzt_soil), init_3d%z_soil,  &
                                       nzb_soil, nzt_soil,                     &
                                       nzb_soil, init_3d%nzs-1 )
                   t_soil_h%var_2d(nzt_soil+1,m) = t_soil_h%var_2d(nzt_soil,m)
                ENDDO
                DO  l = 0, 3
                   DO  m = 1, surf_lsm_v(l)%ns
                      i = surf_lsm_v(l)%i(m) + MERGE( 0, surf_lsm_v(l)%ioff,   &
                                                surf_lsm_v(l)%building_covered(m) ) 
                      j = surf_lsm_v(l)%j(m) + MERGE( 0, surf_lsm_v(l)%joff,   &
                                                surf_lsm_v(l)%building_covered(m) ) 

                      IF ( init_3d%msoil(0,j,i) /= init_3d%fill_msoil )        &
                         CALL netcdf_data_input_interpolate(                   &
                                       t_soil_v(l)%var_2d(nzb_soil:nzt_soil,m),&
                                       init_3d%tsoil(:,j,i),                   &
                                       zs(nzb_soil:nzt_soil), init_3d%z_soil,  &
                                       nzb_soil, nzt_soil,                     &
                                       nzb_soil, init_3d%nzs-1 )
                      t_soil_v(l)%var_2d(nzt_soil+1,m) =                       &
                                                 t_soil_v(l)%var_2d(nzt_soil,m)
                   ENDDO
                ENDDO
             ENDIF
          ENDIF
!
!--       Further initialization
          DO  m = 1, surf_lsm_h%ns

             i   = surf_lsm_h%i(m)            
             j   = surf_lsm_h%j(m)
             k   = surf_lsm_h%k(m)
!
!--          Calculate surface temperature. In case of bare soil, the surface 
!--          temperature must be reset to the soil temperature in the first soil 
!--          layer
             IF ( surf_lsm_h%lambda_surface_s(m) == 0.0_wp )  THEN
                t_surface_h%var_1d(m)    = t_soil_h%var_2d(nzb_soil,m)
                surf_lsm_h%pt_surface(m) = t_soil_h%var_2d(nzb_soil,m) / exn
             ELSE
                t_surface_h%var_1d(m)    = pt(k-1,j,i) * exn
                surf_lsm_h%pt_surface(m) = pt(k-1,j,i)  
             ENDIF
             
             IF ( cloud_physics  .OR. cloud_droplets ) THEN
                surf_lsm_h%pt1(m) = pt(k,j,i) + l_d_cp * pt_d_t(k) * ql(k,j,i)
             ELSE
                surf_lsm_h%pt1(m) = pt(k,j,i)
             ENDIF 


!
!--          Assure that r_a cannot be zero at model start
             IF ( surf_lsm_h%pt1(m) == surf_lsm_h%pt_surface(m) )              &
                surf_lsm_h%pt1(m) = surf_lsm_h%pt1(m) + 1.0E-20_wp

             surf_lsm_h%us(m)   = 0.1_wp
             surf_lsm_h%ts(m)   = ( surf_lsm_h%pt1(m) - surf_lsm_h%pt_surface(m) )&
                                  / surf_lsm_h%r_a(m)
             surf_lsm_h%shf(m)  = - surf_lsm_h%us(m) * surf_lsm_h%ts(m)        &
                                  * rho_surface
         ENDDO
!
!--      Vertical surfaces
         DO  l = 0, 3
             DO  m = 1, surf_lsm_v(l)%ns
                i   = surf_lsm_v(l)%i(m)            
                j   = surf_lsm_v(l)%j(m)
                k   = surf_lsm_v(l)%k(m)         
!
!--             Calculate surface temperature. In case of bare soil, the surface 
!--             temperature must be reset to the soil temperature in the first soil 
!--             layer
                IF ( surf_lsm_v(l)%lambda_surface_s(m) == 0.0_wp )  THEN
                   t_surface_v(l)%var_1d(m)      = t_soil_v(l)%var_2d(nzb_soil,m)
                   surf_lsm_v(l)%pt_surface(m)   = t_soil_v(l)%var_2d(nzb_soil,m) / exn
                ELSE
                   j_off = surf_lsm_v(l)%joff
                   i_off = surf_lsm_v(l)%ioff

                   t_surface_v(l)%var_1d(m)      = pt(k,j+j_off,i+i_off) * exn
                   surf_lsm_v(l)%pt_surface(m)   = pt(k,j+j_off,i+i_off)           
                ENDIF


                IF ( cloud_physics  .OR. cloud_droplets ) THEN
                   surf_lsm_v(l)%pt1(m) = pt(k,j,i) + l_d_cp * pt_d_t(k) * ql(k,j,i)
                ELSE
                   surf_lsm_v(l)%pt1(m) = pt(k,j,i)
                ENDIF 

!
!--             Assure that r_a cannot be zero at model start
                IF ( surf_lsm_v(l)%pt1(m) == surf_lsm_v(l)%pt_surface(m) )     &
                     surf_lsm_v(l)%pt1(m) = surf_lsm_v(l)%pt1(m) + 1.0E-20_wp
!
!--             Set artifical values for ts and us so that r_a has its initial value 
!--             for the first time step. Only for interior core domain, not for ghost points
                surf_lsm_v(l)%us(m)   = 0.1_wp
                surf_lsm_v(l)%ts(m)   = ( surf_lsm_v(l)%pt1(m) - surf_lsm_v(l)%pt_surface(m) ) /&
                                          surf_lsm_v(l)%r_a(m)
                surf_lsm_v(l)%shf(m)  = - surf_lsm_v(l)%us(m) *                &
                                          surf_lsm_v(l)%ts(m) * rho_surface

             ENDDO
          ENDDO
       ENDIF
!
!--    Level 1 initialization of root distribution - provided by the user via
!--    via namelist.
       DO  m = 1, surf_lsm_h%ns
          DO  k = nzb_soil, nzt_soil
             surf_lsm_h%root_fr(k,m) = root_fraction(k)
          ENDDO
       ENDDO

       DO  l = 0, 3
          DO  m = 1, surf_lsm_v(l)%ns
             DO  k = nzb_soil, nzt_soil
                surf_lsm_v(l)%root_fr(k,m) = root_fraction(k)
             ENDDO
          ENDDO
       ENDDO

!
!--    Level 2 initialization of root distribution.
!--    When no root distribution is given by the user, use look-up table to prescribe
!--    the root fraction in the individual soil layers.
       IF ( ALL( root_fraction == 9999999.9_wp ) )  THEN
!
!--       First, calculate the index bounds for integration
          n_soil_layers_total = nzt_soil - nzb_soil + 6
          ALLOCATE ( bound(0:n_soil_layers_total) )
          ALLOCATE ( bound_root_fr(0:n_soil_layers_total) )

          kn = 0
          ko = 0
          bound(0) = 0.0_wp
          DO k = 1, n_soil_layers_total-1
             IF ( zs_layer(kn) <= zs_ref(ko) )  THEN
                bound(k) = zs_layer(kn)
                bound_root_fr(k) = ko
                kn = kn + 1
                IF ( kn > nzt_soil+1 )  THEN
                   kn = nzt_soil
                ENDIF
             ELSE
                bound(k) = zs_ref(ko)
                bound_root_fr(k) = ko
                ko = ko + 1
                IF ( ko > 3 )  THEN
                   ko = 3
                ENDIF
             ENDIF

          ENDDO

!
!--       Integrate over all soil layers based on the four-layer root fraction
          kzs = 1
          root_fraction = 0.0_wp
          DO k = 0, n_soil_layers_total-2
             kroot = bound_root_fr(k+1)
             root_fraction(kzs-1) = root_fraction(kzs-1)                       &
                                + root_distribution(kroot,vegetation_type)     &
                                / dz_soil_ref(kroot) * ( bound(k+1) - bound(k) )

             IF ( bound(k+1) == zs_layer(kzs-1) )  THEN
                kzs = kzs+1
             ENDIF
          ENDDO


!
!--       Normalize so that the sum of all fractions equals one
          root_fraction = root_fraction / SUM(root_fraction)

          DEALLOCATE ( bound )
          DEALLOCATE ( bound_root_fr )

!
!--       Map calculated root fractions
          DO  m = 1, surf_lsm_h%ns
             DO  k = nzb_soil, nzt_soil 
                IF ( surf_lsm_h%pavement_surface(m)  .AND.                     &
                     k <= surf_lsm_h%nzt_pavement(m) )  THEN
                   surf_lsm_h%root_fr(k,m) = 0.0_wp
                ELSE
                   surf_lsm_h%root_fr(k,m) = root_fraction(k)
                ENDIF

             ENDDO
!
!--          Normalize so that the sum = 1. Only relevant when the root         
!--          distribution was set to zero due to pavement at some layers.
             IF ( SUM( surf_lsm_h%root_fr(:,m) ) > 0.0_wp )  THEN
                DO k = nzb_soil, nzt_soil
                   surf_lsm_h%root_fr(k,m) = surf_lsm_h%root_fr(k,m)           &
                   / SUM( surf_lsm_h%root_fr(:,m) )
                ENDDO
             ENDIF
          ENDDO
          DO  l = 0, 3
             DO  m = 1, surf_lsm_v(l)%ns
                DO  k = nzb_soil, nzt_soil
                   IF ( surf_lsm_v(l)%pavement_surface(m)  .AND.               &
                        k <= surf_lsm_v(l)%nzt_pavement(m) )  THEN
                      surf_lsm_v(l)%root_fr(k,m) = 0.0_wp
                   ELSE
                      surf_lsm_v(l)%root_fr(k,m) = root_fraction(k)
                   ENDIF
                ENDDO
!
!--             Normalize so that the sum = 1. Only relevant when the root      
!--             distribution was set to zero due to pavement at some layers.
                IF ( SUM( surf_lsm_v(l)%root_fr(:,m) ) > 0.0_wp )  THEN
                   DO  k = nzb_soil, nzt_soil  
                      surf_lsm_v(l)%root_fr(k,m) = surf_lsm_v(l)%root_fr(k,m)  &
                      / SUM( surf_lsm_v(l)%root_fr(:,m) )
                   ENDDO
                ENDIF
             ENDDO
           ENDDO
       ENDIF
!
!--    Level 3 initialization of root distribution.
!--    Take value from file
       IF ( root_area_density_lsm_f%from_file )  THEN
          DO  m = 1, surf_lsm_h%ns
             IF ( surf_lsm_h%vegetation_surface(m) )  THEN
                i = surf_lsm_h%i(m) + MERGE( 0, surf_lsm_v(l)%ioff,            &
                                             surf_lsm_v(l)%building_covered(m) ) 
                j = surf_lsm_h%j(m) + MERGE( 0, surf_lsm_v(l)%joff,            &
                                             surf_lsm_v(l)%building_covered(m) ) 
                DO  k = nzb_soil, nzt_soil 
                   surf_lsm_h%root_fr(k,m) = root_area_density_lsm_f%var(k,j,i) 
                ENDDO

             ENDIF
          ENDDO

          DO  l = 0, 3
             DO  m = 1, surf_lsm_v(l)%ns
                IF ( surf_lsm_v(l)%vegetation_surface(m) )  THEN
                   i = surf_lsm_v(l)%i(m) + MERGE( 0, surf_lsm_v(l)%ioff,      &
                                                   surf_lsm_v(l)%building_covered(m) ) 
                   j = surf_lsm_v(l)%j(m) + MERGE( 0, surf_lsm_v(l)%joff,      &
                                                   surf_lsm_v(l)%building_covered(m) ) 

                   DO  k = nzb_soil, nzt_soil 
                      surf_lsm_v(l)%root_fr(k,m) = root_area_density_lsm_f%var(k,j,i) 
                   ENDDO

                ENDIF
             ENDDO
          ENDDO

       ENDIF
 
!
!--    Possibly do user-defined actions (e.g. define heterogeneous land surface)
       CALL user_init_land_surface


!
!--    Calculate new roughness lengths (for water surfaces only, i.e. only 
!-     horizontal surfaces)
       IF ( .NOT. constant_roughness )  CALL calc_z0_water_surface

       t_soil_h_p    = t_soil_h
       m_soil_h_p    = m_soil_h
       m_liq_h_p     = m_liq_h
       t_surface_h_p = t_surface_h

       t_soil_v_p    = t_soil_v
       m_soil_v_p    = m_soil_v
       m_liq_v_p     = m_liq_v
       t_surface_v_p = t_surface_v



!--    Store initial profiles of t_soil and m_soil (assuming they are 
!--    horizontally homogeneous on this PE)
!--    DEACTIVATED FOR NOW - leads to error when number of locations with
!--    soil model is zero on a PE.
!        hom(nzb_soil:nzt_soil,1,90,:)  = SPREAD( t_soil_h%var_2d(nzb_soil:nzt_soil,1),  &
!                                                 2, statistic_regions+1 )
!        hom(nzb_soil:nzt_soil,1,92,:)  = SPREAD( m_soil_h%var_2d(nzb_soil:nzt_soil,1),  & 
!                                                 2, statistic_regions+1 )

!
!--    Finally, make some consistency checks. 
!--    Ceck for illegal combination of LAI and vegetation coverage.
       IF ( ANY( .NOT. surf_lsm_h%pavement_surface  .AND.                      &
                 surf_lsm_h%lai == 0.0_wp  .AND.  surf_lsm_h%c_veg == 1.0_wp ) &
          )  THEN
          message_string = 'For non-pavement surfaces the combination ' //     &
                           ' lai = 0.0 and c_veg = 1.0 is not allowed.'
          CALL message( 'lsm_rrd_local', 'PA0999', 2, 2, 0, 6, 0 )
       ENDIF

       DO  l = 0, 3
          IF ( ANY( .NOT. surf_lsm_v(l)%pavement_surface  .AND.                &
                    surf_lsm_v(l)%lai == 0.0_wp  .AND.                         &
                    surf_lsm_v(l)%c_veg == 1.0_wp ) )  THEN
             message_string = 'For non-pavement surfaces the combination ' //  &
                              ' lai = 0.0 and c_veg = 1.0 is not allowed.'
             CALL message( 'lsm_rrd_local', 'PA0999', 2, 2, 0, 6, 0 )
          ENDIF
       ENDDO
!
!--    Check if roughness length exceed surface-layer height and decrease 
!--    local roughness length where necessary. 
       DO  m = 1, surf_lsm_h%ns
          IF ( surf_lsm_h%z0(m) >= surf_lsm_h%z_mo(m) )  THEN
          
             surf_lsm_h%z0(m) = 0.9_wp * surf_lsm_h%z_mo(m)
             
             WRITE( message_string, * ) 'z0 exceeds surface-layer height ' //  &
                            'at horizontal natural surface and is ' //         &
                            'decreased appropriately at grid point (i,j) = ',  &
                            surf_lsm_h%i(m), surf_lsm_h%j(m)
             CALL message( 'land_surface_model_mod', 'PA0503',                 &
                            0, 0, 0, 6, 0 )
          ENDIF
       ENDDO
       
       DO  l = 0, 3
          DO  m = 1, surf_lsm_v(l)%ns
             IF ( surf_lsm_v(l)%z0(m) >= surf_lsm_v(l)%z_mo(m) )  THEN
          
                surf_lsm_v(l)%z0(m) = 0.9_wp * surf_lsm_v(l)%z_mo(m)
             
                WRITE( message_string, * ) 'z0 exceeds surface-layer height '//&
                            'at vertical natural surface and is ' //           &
                            'decreased appropriately at grid point (i,j) = ',  &
                            surf_lsm_v(l)%i(m)+surf_lsm_v(l)%ioff,             &
                            surf_lsm_v(l)%j(m)+surf_lsm_v(l)%joff
                CALL message( 'land_surface_model_mod', 'PA0503',              &
                            0, 0, 0, 6, 0 )
             ENDIF
          ENDDO
       ENDDO      



    END SUBROUTINE lsm_init


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate land surface model arrays and define pointers
!------------------------------------------------------------------------------!
    SUBROUTINE lsm_init_arrays
    

       IMPLICIT NONE

       INTEGER(iwp) ::  l !< index indicating facing of surface array 
   
       ALLOCATE ( root_extr(nzb_soil:nzt_soil) )
       root_extr = 0.0_wp 
       
!
!--    Allocate surface and soil temperature / humidity. Please note, 
!--    these arrays are allocated according to surface-data structure,
!--    even if they do not belong to the data type due to the 
!--    pointer arithmetric (TARGET attribute is not allowed in a data-type). 
#if defined( __nopointer )
!
!--    Horizontal surfaces
       ALLOCATE ( m_liq_h_p%var_1d(1:surf_lsm_h%ns)                      )
       ALLOCATE ( t_surface_h%var_1d(1:surf_lsm_h%ns)                    )
       ALLOCATE ( t_surface_h_p%var_1d(1:surf_lsm_h%ns)                  )
       ALLOCATE ( m_soil_h_p%var_2d(nzb_soil:nzt_soil,1:surf_lsm_h%ns)   )
       ALLOCATE ( t_soil_h_p%var_2d(nzb_soil:nzt_soil+1,1:surf_lsm_h%ns) )

!
!--    Vertical surfaces
       DO  l = 0, 3
          ALLOCATE ( m_liq_v(l)%var_1d(1:surf_lsm_v(l)%ns)                        )
          ALLOCATE ( m_liq_v_p(l)%var_1d(1:surf_lsm_v(l)%ns)                      )
          ALLOCATE ( t_surface_v(l)%var_1d(1:surf_lsm_v(l)%ns)                    )
          ALLOCATE ( t_surface_v_p(l)%var_1d(1:surf_lsm_v(l)%ns)                  )
          ALLOCATE ( m_soil_v(l)%var_2d(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns)     )
          ALLOCATE ( m_soil_v_p(l)%var_2d(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns)   )
          ALLOCATE ( t_soil_v(l)%var_2d(nzb_soil:nzt_soil+1,1:surf_lsm_v(l)%ns)   )
          ALLOCATE ( t_soil_v_p(l)%var_2d(nzb_soil:nzt_soil+1,1:surf_lsm_v(l)%ns) )
       ENDDO
!
!--    Allocate soil temperature and moisture. As these variables might be 
!--    already allocated in case of restarts, check this. 
       IF ( .NOT. ALLOCATED( m_liq_h%var_1d ) )                                &
          ALLOCATE ( m_liq_h%var_1d(1:surf_lsm_h%ns) )
       IF ( .NOT. ALLOCATED( m_soil_h%var_2d ) )                               &
          ALLOCATE ( m_soil_h%var_2d(nzb_soil:nzt_soil,1:surf_lsm_h%ns) )
       IF ( .NOT. ALLOCATED( t_soil_h%var_2d ) )                               &
          ALLOCATE ( t_soil_h%var_2d(nzb_soil:nzt_soil,1:surf_lsm_h%ns) )

       DO  l = 0, 3
          IF ( .NOT. ALLOCATED( m_liq_v(l)%var_1d ) )                          &
             ALLOCATE ( m_liq_v(l)%var_1d(1:surf_lsm_v(l)%ns) )
          IF ( .NOT. ALLOCATED( m_soil_v(l)%var_2d ) )                         &
             ALLOCATE ( m_soil_v(l)%var_2d(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns) )
          IF ( .NOT. ALLOCATED( t_soil_v(l)%var_2d ) )                         &
             ALLOCATE ( t_soil_v(l)%var_2d(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns) )
       ENDDO
#else
!
!--    Horizontal surfaces
       ALLOCATE ( m_liq_h_1%var_1d(1:surf_lsm_h%ns)                      )
       ALLOCATE ( m_liq_h_2%var_1d(1:surf_lsm_h%ns)                      )
       ALLOCATE ( t_surface_h_1%var_1d(1:surf_lsm_h%ns)                  )
       ALLOCATE ( t_surface_h_2%var_1d(1:surf_lsm_h%ns)                  )
       ALLOCATE ( m_soil_h_1%var_2d(nzb_soil:nzt_soil,1:surf_lsm_h%ns)   )
       ALLOCATE ( m_soil_h_2%var_2d(nzb_soil:nzt_soil,1:surf_lsm_h%ns)   )
       ALLOCATE ( t_soil_h_1%var_2d(nzb_soil:nzt_soil+1,1:surf_lsm_h%ns) )
       ALLOCATE ( t_soil_h_2%var_2d(nzb_soil:nzt_soil+1,1:surf_lsm_h%ns) )
!
!--    Vertical surfaces
       DO  l = 0, 3
          ALLOCATE ( m_liq_v_1(l)%var_1d(1:surf_lsm_v(l)%ns)                      )
          ALLOCATE ( m_liq_v_2(l)%var_1d(1:surf_lsm_v(l)%ns)                      )
          ALLOCATE ( t_surface_v_1(l)%var_1d(1:surf_lsm_v(l)%ns)                  )
          ALLOCATE ( t_surface_v_2(l)%var_1d(1:surf_lsm_v(l)%ns)                  )
          ALLOCATE ( m_soil_v_1(l)%var_2d(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns)   )
          ALLOCATE ( m_soil_v_2(l)%var_2d(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns)   )
          ALLOCATE ( t_soil_v_1(l)%var_2d(nzb_soil:nzt_soil+1,1:surf_lsm_v(l)%ns) )
          ALLOCATE ( t_soil_v_2(l)%var_2d(nzb_soil:nzt_soil+1,1:surf_lsm_v(l)%ns) )
       ENDDO
#endif
!
!--    Allocate array for heat flux in W/m2, required for radiation?
!--    Consider to remove this array
       ALLOCATE( surf_lsm_h%surfhf(1:surf_lsm_h%ns) )
       DO  l = 0, 3
          ALLOCATE( surf_lsm_v(l)%surfhf(1:surf_lsm_v(l)%ns) )
       ENDDO


!
!--    Allocate intermediate timestep arrays
!--    Horizontal surfaces
       ALLOCATE ( tm_liq_h_m%var_1d(1:surf_lsm_h%ns)                     )
       ALLOCATE ( tt_surface_h_m%var_1d(1:surf_lsm_h%ns)                 )
       ALLOCATE ( tm_soil_h_m%var_2d(nzb_soil:nzt_soil,1:surf_lsm_h%ns)  )
       ALLOCATE ( tt_soil_h_m%var_2d(nzb_soil:nzt_soil,1:surf_lsm_h%ns)  ) 
!
!--    Horizontal surfaces
       DO  l = 0, 3
          ALLOCATE ( tm_liq_v_m(l)%var_1d(1:surf_lsm_v(l)%ns)                     )
          ALLOCATE ( tt_surface_v_m(l)%var_1d(1:surf_lsm_v(l)%ns)                 )
          ALLOCATE ( tm_soil_v_m(l)%var_2d(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns)  )
          ALLOCATE ( tt_soil_v_m(l)%var_2d(nzb_soil:nzt_soil,1:surf_lsm_v(l)%ns)  )
       ENDDO 

!
!--    Allocate 2D vegetation model arrays
!--    Horizontal surfaces
       ALLOCATE ( surf_lsm_h%building_surface(1:surf_lsm_h%ns)    )
       ALLOCATE ( surf_lsm_h%c_liq(1:surf_lsm_h%ns)               )
       ALLOCATE ( surf_lsm_h%c_surface(1:surf_lsm_h%ns)           )
       ALLOCATE ( surf_lsm_h%c_veg(1:surf_lsm_h%ns)               )
       ALLOCATE ( surf_lsm_h%f_sw_in(1:surf_lsm_h%ns)             )
       ALLOCATE ( surf_lsm_h%ghf(1:surf_lsm_h%ns)                 )
       ALLOCATE ( surf_lsm_h%g_d(1:surf_lsm_h%ns)                 )
       ALLOCATE ( surf_lsm_h%lai(1:surf_lsm_h%ns)                 )
       ALLOCATE ( surf_lsm_h%lambda_surface_u(1:surf_lsm_h%ns)    )
       ALLOCATE ( surf_lsm_h%lambda_surface_s(1:surf_lsm_h%ns)    )
       ALLOCATE ( surf_lsm_h%nzt_pavement(1:surf_lsm_h%ns)        )
       ALLOCATE ( surf_lsm_h%pavement_surface(1:surf_lsm_h%ns)    )
       ALLOCATE ( surf_lsm_h%qsws_soil(1:surf_lsm_h%ns)           )
       ALLOCATE ( surf_lsm_h%qsws_liq(1:surf_lsm_h%ns)            )
       ALLOCATE ( surf_lsm_h%qsws_veg(1:surf_lsm_h%ns)            )
       ALLOCATE ( surf_lsm_h%rad_net_l(1:surf_lsm_h%ns)           ) 
       ALLOCATE ( surf_lsm_h%r_a(1:surf_lsm_h%ns)                 )
       ALLOCATE ( surf_lsm_h%r_canopy(1:surf_lsm_h%ns)            )
       ALLOCATE ( surf_lsm_h%r_soil(1:surf_lsm_h%ns)              )
       ALLOCATE ( surf_lsm_h%r_soil_min(1:surf_lsm_h%ns)          )
       ALLOCATE ( surf_lsm_h%r_s(1:surf_lsm_h%ns)                 )
       ALLOCATE ( surf_lsm_h%r_canopy_min(1:surf_lsm_h%ns)        )
       ALLOCATE ( surf_lsm_h%vegetation_surface(1:surf_lsm_h%ns)  )
       ALLOCATE ( surf_lsm_h%water_surface(1:surf_lsm_h%ns)       )

       surf_lsm_h%water_surface        = .FALSE.
       surf_lsm_h%pavement_surface     = .FALSE.
       surf_lsm_h%vegetation_surface   = .FALSE. 

!
!--    Set default values
       surf_lsm_h%r_canopy_min = 0.0_wp

!
!--    Vertical surfaces
       DO  l = 0, 3
          ALLOCATE ( surf_lsm_v(l)%building_surface(1:surf_lsm_v(l)%ns)    )
          ALLOCATE ( surf_lsm_v(l)%c_liq(1:surf_lsm_v(l)%ns)               )
          ALLOCATE ( surf_lsm_v(l)%c_surface(1:surf_lsm_v(l)%ns)           )
          ALLOCATE ( surf_lsm_v(l)%c_veg(1:surf_lsm_v(l)%ns)               )
          ALLOCATE ( surf_lsm_v(l)%f_sw_in(1:surf_lsm_v(l)%ns)             )
          ALLOCATE ( surf_lsm_v(l)%ghf(1:surf_lsm_v(l)%ns)                 )
          ALLOCATE ( surf_lsm_v(l)%g_d(1:surf_lsm_v(l)%ns)                 )
          ALLOCATE ( surf_lsm_v(l)%lai(1:surf_lsm_v(l)%ns)                 )
          ALLOCATE ( surf_lsm_v(l)%lambda_surface_u(1:surf_lsm_v(l)%ns)    )
          ALLOCATE ( surf_lsm_v(l)%lambda_surface_s(1:surf_lsm_v(l)%ns)    )
          ALLOCATE ( surf_lsm_v(l)%nzt_pavement(1:surf_lsm_v(l)%ns)        )
          ALLOCATE ( surf_lsm_v(l)%pavement_surface(1:surf_lsm_v(l)%ns)    )
          ALLOCATE ( surf_lsm_v(l)%qsws_soil(1:surf_lsm_v(l)%ns)           )
          ALLOCATE ( surf_lsm_v(l)%qsws_liq(1:surf_lsm_v(l)%ns)            )
          ALLOCATE ( surf_lsm_v(l)%qsws_veg(1:surf_lsm_v(l)%ns)            )
          ALLOCATE ( surf_lsm_v(l)%rad_net_l(1:surf_lsm_v(l)%ns)           )
          ALLOCATE ( surf_lsm_v(l)%r_a(1:surf_lsm_v(l)%ns)                 )
          ALLOCATE ( surf_lsm_v(l)%r_canopy(1:surf_lsm_v(l)%ns)            )
          ALLOCATE ( surf_lsm_v(l)%r_soil(1:surf_lsm_v(l)%ns)              )
          ALLOCATE ( surf_lsm_v(l)%r_soil_min(1:surf_lsm_v(l)%ns)          )
          ALLOCATE ( surf_lsm_v(l)%r_s(1:surf_lsm_v(l)%ns)                 )
          ALLOCATE ( surf_lsm_v(l)%r_canopy_min(1:surf_lsm_v(l)%ns)        )
          ALLOCATE ( surf_lsm_v(l)%vegetation_surface(1:surf_lsm_v(l)%ns)  )
          ALLOCATE ( surf_lsm_v(l)%water_surface(1:surf_lsm_v(l)%ns)       )

          surf_lsm_v(l)%water_surface       = .FALSE.
          surf_lsm_v(l)%pavement_surface    = .FALSE.
          surf_lsm_v(l)%vegetation_surface  = .FALSE. 
          

!
!--       Set default values
          surf_lsm_v(l)%r_canopy_min = 0.0_wp
        
       ENDDO

   
#if ! defined( __nopointer )
!
!--    Initial assignment of the pointers
!--    Horizontal surfaces
       t_soil_h    => t_soil_h_1;    t_soil_h_p    => t_soil_h_2
       t_surface_h => t_surface_h_1; t_surface_h_p => t_surface_h_2
       m_soil_h    => m_soil_h_1;    m_soil_h_p    => m_soil_h_2
       m_liq_h     => m_liq_h_1;     m_liq_h_p     => m_liq_h_2
!
!--    Vertical surfaces
       t_soil_v    => t_soil_v_1;    t_soil_v_p    => t_soil_v_2
       t_surface_v => t_surface_v_1; t_surface_v_p => t_surface_v_2
       m_soil_v    => m_soil_v_1;    m_soil_v_p    => m_soil_v_2
       m_liq_v     => m_liq_v_1;     m_liq_v_p     => m_liq_v_2

#endif


    END SUBROUTINE lsm_init_arrays


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &lsmpar for land surface model
!------------------------------------------------------------------------------!
    SUBROUTINE lsm_parin

       USE control_parameters,                                                 &
           ONLY:  message_string

       IMPLICIT NONE

       CHARACTER (LEN=80) ::  line  !< dummy string that contains the current line of the parameter file 
       
       NAMELIST /lsm_par/         alpha_vangenuchten, c_surface,               &
                                  canopy_resistance_coefficient,               &
                                  constant_roughness,                          &
                                  conserve_water_content,                      &
                                  deep_soil_temperature,                       &
                                  dz_soil,                                     &
                                  f_shortwave_incoming, field_capacity,        & 
                                  aero_resist_kray, hydraulic_conductivity,    &
                                  lambda_surface_stable,                       &
                                  lambda_surface_unstable, leaf_area_index,    &
                                  l_vangenuchten, min_canopy_resistance,       &
                                  min_soil_resistance, n_vangenuchten,         &
                                  pavement_depth_level,                        &
                                  pavement_heat_capacity,                      &
                                  pavement_heat_conduct, pavement_type,        &
                                  residual_moisture, root_fraction,            &
                                  saturation_moisture, skip_time_do_lsm,       &
                                  soil_moisture, soil_temperature,             &
                                  soil_type,                                   &
                                  surface_type,                                &
                                  vegetation_coverage, vegetation_type,        &
                                  water_temperature, water_type,               &
                                  wilting_point, z0_vegetation,                &
                                  z0h_vegetation, z0q_vegetation, z0_water,    &
                                  z0h_water, z0q_water, z0_pavement,           &
                                  z0h_pavement, z0q_pavement

       NAMELIST /land_surface_parameters/                                      &
                                  alpha_vangenuchten, c_surface,               &
                                  canopy_resistance_coefficient,               &
                                  constant_roughness,                          &
                                  conserve_water_content,                      &
                                  deep_soil_temperature,                       &
                                  dz_soil,                                     &
                                  f_shortwave_incoming, field_capacity,        & 
                                  aero_resist_kray, hydraulic_conductivity,    &
                                  lambda_surface_stable,                       &
                                  lambda_surface_unstable, leaf_area_index,    &
                                  l_vangenuchten, min_canopy_resistance,       &
                                  min_soil_resistance, n_vangenuchten,         &
                                  pavement_depth_level,                        &
                                  pavement_heat_capacity,                      &
                                  pavement_heat_conduct, pavement_type,        &
                                  residual_moisture, root_fraction,            &
                                  saturation_moisture, skip_time_do_lsm,       &
                                  soil_moisture, soil_temperature,             &
                                  soil_type,                                   &
                                  surface_type,                                &
                                  vegetation_coverage, vegetation_type,        &
                                  water_temperature, water_type,               &
                                  wilting_point, z0_vegetation,                &
                                  z0h_vegetation, z0q_vegetation, z0_water,    &
                                  z0h_water, z0q_water, z0_pavement,           &
                                  z0h_pavement, z0q_pavement
                                  
       line = ' '
  
!
!--    Try to find land surface model package
       REWIND ( 11 )
       line = ' '
       DO   WHILE ( INDEX( line, '&land_surface_parameters' ) == 0 )
          READ ( 11, '(A)', END=10 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read user-defined namelist
       READ ( 11, land_surface_parameters )

!
!--    Set flag that indicates that the land surface model is switched on
       land_surface = .TRUE.
       
       GOTO 12
!
!--    Try to find old namelist
 10    REWIND ( 11 )
       line = ' '
       DO   WHILE ( INDEX( line, '&lsm_par' ) == 0 )
          READ ( 11, '(A)', END=12 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read user-defined namelist
       READ ( 11, lsm_par )

       message_string = 'namelist lsm_par is deprecated and will be ' // &
                     'removed in near future. Please use namelist ' //   &
                     'land_surface_parameters instead'
       CALL message( 'lsm_parin', 'PA0487', 0, 1, 0, 6, 0 )
       
!
!--    Set flag that indicates that the land surface model is switched on
       land_surface = .TRUE.


 12    CONTINUE
       

    END SUBROUTINE lsm_parin


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Soil model as part of the land surface model. The model predicts soil
!> temperature and water content.
!------------------------------------------------------------------------------!
    SUBROUTINE lsm_soil_model( horizontal, l, calc_soil_moisture )


       IMPLICIT NONE

       INTEGER(iwp) ::  k       !< running index
       INTEGER(iwp) ::  l       !< surface-data type index indication facing
       INTEGER(iwp) ::  m       !< running index

       LOGICAL, INTENT(IN) ::  calc_soil_moisture !< flag indicating whether soil moisture shall be calculated or not.

       LOGICAL      ::  horizontal !< flag indication horizontal wall, required to set pointer accordingly

       REAL(wp)     ::  h_vg !< Van Genuchten coef. h

       REAL(wp), DIMENSION(nzb_soil:nzt_soil) :: gamma_temp,  & !< temp. gamma
                                                 lambda_temp, & !< temp. lambda
                                                 tend           !< tendency

       TYPE(surf_type_lsm), POINTER ::  surf_m_soil
       TYPE(surf_type_lsm), POINTER ::  surf_m_soil_p
       TYPE(surf_type_lsm), POINTER ::  surf_t_soil
       TYPE(surf_type_lsm), POINTER ::  surf_t_soil_p
       TYPE(surf_type_lsm), POINTER ::  surf_tm_soil_m
       TYPE(surf_type_lsm), POINTER ::  surf_tt_soil_m

       TYPE(surf_type), POINTER  ::  surf  !< surface-date type variable

       IF ( horizontal )  THEN
          surf           => surf_lsm_h

          surf_m_soil    => m_soil_h
          surf_m_soil_p  => m_soil_h_p
          surf_t_soil    => t_soil_h
          surf_t_soil_p  => t_soil_h_p
          surf_tm_soil_m => tm_soil_h_m
          surf_tt_soil_m => tt_soil_h_m
       ELSE
          surf           => surf_lsm_v(l)

          surf_m_soil    => m_soil_v(l)
          surf_m_soil_p  => m_soil_v_p(l)
          surf_t_soil    => t_soil_v(l)
          surf_t_soil_p  => t_soil_v_p(l)
          surf_tm_soil_m => tm_soil_v_m(l)
          surf_tt_soil_m => tt_soil_v_m(l)
       ENDIF

       DO  m = 1, surf%ns

          IF (  .NOT.  surf%water_surface(m) )  THEN
             DO  k = nzb_soil, nzt_soil

                IF ( surf%pavement_surface(m)  .AND.                           &
                     k <= surf%nzt_pavement(m) )  THEN
                   
                   surf%rho_c_total(k,m) = surf%rho_c_total_def(k,m)
                   lambda_temp(k)        = surf%lambda_h_def(k,m) 

                ELSE            
!
!--                Calculate volumetric heat capacity of the soil, taking 
!--                into account water content 
                   surf%rho_c_total(k,m) = (rho_c_soil *                       &
                                               ( 1.0_wp - surf%m_sat(k,m) )    &
                                               + rho_c_water * surf_m_soil%var_2d(k,m) )

!
!--                Calculate soil heat conductivity at the center of the soil
!--                layers
                   lambda_h_sat = lambda_h_sm**(1.0_wp - surf%m_sat(k,m)) *    &
                                  lambda_h_water ** surf_m_soil%var_2d(k,m)

                   ke = 1.0_wp + LOG10( MAX( 0.1_wp, surf_m_soil%var_2d(k,m) / &
                                                     surf%m_sat(k,m) ) )

                   lambda_temp(k) = ke * (lambda_h_sat - lambda_h_dry) +       &
                                    lambda_h_dry
                ENDIF
             ENDDO

!
!--          Calculate soil heat conductivity (lambda_h) at the _layer level 
!--          using linear interpolation. For pavement surface, the
!--          true pavement depth is considered
             DO  k = nzb_soil, nzt_soil-1
                   surf%lambda_h(k,m) = ( lambda_temp(k+1) + lambda_temp(k) )  &
                                        * 0.5_wp
             ENDDO
             surf%lambda_h(nzt_soil,m) = lambda_temp(nzt_soil)

!
!--          Prognostic equation for soil temperature t_soil
             tend(:) = 0.0_wp

             tend(nzb_soil) = ( 1.0_wp / surf%rho_c_total(nzb_soil,m) ) *            &
                    ( surf%lambda_h(nzb_soil,m) * ( surf_t_soil%var_2d(nzb_soil+1,m) &
                      - surf_t_soil%var_2d(nzb_soil,m) ) * ddz_soil_center(nzb_soil) &
                      + surf%ghf(m) ) * ddz_soil(nzb_soil)

             DO  k = nzb_soil+1, nzt_soil
                tend(k) = ( 1.0_wp / surf%rho_c_total(k,m) )                   &
                          * (   surf%lambda_h(k,m)                             &
                     * ( surf_t_soil%var_2d(k+1,m) - surf_t_soil%var_2d(k,m) ) &
                     * ddz_soil_center(k)                                      &
                     - surf%lambda_h(k-1,m)                                    &
                     * ( surf_t_soil%var_2d(k,m) - surf_t_soil%var_2d(k-1,m) ) &
                     * ddz_soil_center(k-1)                                    &
                            ) * ddz_soil(k)

             ENDDO

             surf_t_soil_p%var_2d(nzb_soil:nzt_soil,m) =                       &
                                       surf_t_soil%var_2d(nzb_soil:nzt_soil,m) &
                                               + dt_3d * ( tsc(2)              &
                                               * tend(nzb_soil:nzt_soil)       & 
                                               + tsc(3)                        &
                                               * surf_tt_soil_m%var_2d(nzb_soil:nzt_soil,m) )

!
!--          Calculate t_soil tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   DO  k = nzb_soil, nzt_soil
                      surf_tt_soil_m%var_2d(k,m) = tend(k)
                   ENDDO
                ELSEIF ( intermediate_timestep_count <                         &
                         intermediate_timestep_count_max )  THEN
                   DO  k = nzb_soil, nzt_soil
                      surf_tt_soil_m%var_2d(k,m) = -9.5625_wp * tend(k) +      &
                                                    5.3125_wp *                &
                                                      surf_tt_soil_m%var_2d(k,m)
                   ENDDO
                ENDIF
             ENDIF


             DO  k = nzb_soil, nzt_soil

!
!--             In order to prevent water tranport through paved surfaces, 
!--             conductivity and diffusivity are set to zero
                IF ( surf%pavement_surface(m)  .AND.                           &
                     k <= surf%nzt_pavement(m) )  THEN
                   lambda_temp(k) = 0.0_wp
                   gamma_temp(k)  = 0.0_wp
   
                ELSE 
    
!
!--                Calculate soil diffusivity at the center of the soil layers
                   lambda_temp(k) = (- b_ch * surf%gamma_w_sat(k,m) * psi_sat  &
                                    / surf%m_sat(k,m) ) * (                    &
                                    MAX( surf_m_soil%var_2d(k,m),              &
                                    surf%m_wilt(k,m) ) / surf%m_sat(k,m) )**(  &
                                    b_ch + 2.0_wp )

!
!--                Parametrization of Van Genuchten
!--                Calculate the hydraulic conductivity after Van Genuchten (1980)
                   h_vg = ( ( ( surf%m_res(k,m) - surf%m_sat(k,m) ) /          &
                              ( surf%m_res(k,m) -                              &
                                MAX( surf_m_soil%var_2d(k,m), surf%m_wilt(k,m) )&
                              )                                                &
                            )**(                                               &
                          surf%n_vg(k,m) / ( surf%n_vg(k,m) - 1.0_wp )         &
                               ) - 1.0_wp                                      &
                          )**( 1.0_wp / surf%n_vg(k,m) ) / surf%alpha_vg(k,m)

                   gamma_temp(k) = surf%gamma_w_sat(k,m) * ( ( ( 1.0_wp +      &
                          ( surf%alpha_vg(k,m) * h_vg )**surf%n_vg(k,m)        &
                                                                  )**(         &
                              1.0_wp - 1.0_wp / surf%n_vg(k,m)) - (            &
                          surf%alpha_vg(k,m) * h_vg )**( surf%n_vg(k,m)        &
                              - 1.0_wp) )**2 )                                 &
                              / ( ( 1.0_wp + ( surf%alpha_vg(k,m) * h_vg       &
                              )**surf%n_vg(k,m) )**( ( 1.0_wp  - 1.0_wp        &
                              / surf%n_vg(k,m) ) *                             &
                              ( surf%l_vg(k,m) + 2.0_wp) ) )

                ENDIF

             ENDDO


             IF (  calc_soil_moisture )  THEN

!
!--             Prognostic equation for soil moisture content. Only performed,
!--             when humidity is enabled in the atmosphere.
                IF ( humidity )  THEN
!
!--                Calculate soil diffusivity (lambda_w) at the _layer level 
!--                using linear interpolation. To do: replace this with
!--                ECMWF-IFS Eq. 8.81
                   DO  k = nzb_soil, nzt_soil-1
                
                      surf%lambda_w(k,m) = ( lambda_temp(k+1) + lambda_temp(k) )  &
                                           * 0.5_wp
                      surf%gamma_w(k,m)  = ( gamma_temp(k+1)  +  gamma_temp(k) )  &
                                           * 0.5_wp
                                            
                   ENDDO
!
!
!--                In case of a closed bottom (= water content is conserved), 
!--                set hydraulic conductivity to zero to that no water will be 
!--                lost in the bottom layer. As gamma_w is always a positive value,
!--                it cannot be set to zero in case of purely dry soil since this
!--                would cause accumulation of (non-existing) water in the lowest
!--                soil layer
                   IF ( conserve_water_content .AND.                           &
                        surf_m_soil%var_2d(nzt_soil,m) /= 0.0_wp )  THEN

                      surf%gamma_w(nzt_soil,m) = 0.0_wp
                   ELSE
                      surf%gamma_w(nzt_soil,m) = gamma_temp(nzt_soil)
                   ENDIF     

!--                The root extraction (= root_extr * qsws_veg / (rho_l     
!--                * l_v)) ensures the mass conservation for water. The         
!--                transpiration of plants equals the cumulative withdrawals by 
!--                the roots in the soil. The scheme takes into account the 
!--                availability of water in the soil layers as well as the root 
!--                fraction in the respective layer. Layer with moisture below 
!--                wilting point will not contribute, which reflects the 
!--                preference of plants to take water from moister layers.
!
!--                Calculate the root extraction (ECMWF 7.69, the sum of 
!--                root_extr = 1). The energy balance solver guarantees a 
!--                positive transpiration, so that there is no need for an 
!--                additional check.
                   m_total = 0.0_wp
                   DO  k = nzb_soil, nzt_soil
                      IF ( surf_m_soil%var_2d(k,m) > surf%m_wilt(k,m) )  THEN
                         m_total = m_total + surf%root_fr(k,m)                 &
                                * surf_m_soil%var_2d(k,m)
                      ENDIF
                   ENDDO  
                   IF ( m_total > 0.0_wp )  THEN
                      DO  k = nzb_soil, nzt_soil
                         IF ( surf_m_soil%var_2d(k,m) > surf%m_wilt(k,m) )  THEN
                            root_extr(k) = surf%root_fr(k,m)                   &
                                           * surf_m_soil%var_2d(k,m) / m_total
                         ELSE
                            root_extr(k) = 0.0_wp
                         ENDIF
                      ENDDO
                   ENDIF
!
!--                Prognostic equation for soil water content m_soil_h.
                   tend(:) = 0.0_wp

                   tend(nzb_soil) = ( surf%lambda_w(nzb_soil,m) *   (          &
                         surf_m_soil%var_2d(nzb_soil+1,m)                      &
                         - surf_m_soil%var_2d(nzb_soil,m) )                    &
                         * ddz_soil_center(nzb_soil) - surf%gamma_w(nzb_soil,m)&
                         - ( root_extr(nzb_soil) * surf%qsws_veg(m)            &
                            + surf%qsws_soil(m) ) * drho_l_lv )                &
                            * ddz_soil(nzb_soil)

                   DO  k = nzb_soil+1, nzt_soil-1
                      tend(k) = ( surf%lambda_w(k,m) * ( surf_m_soil%var_2d(k+1,m)  &
                             - surf_m_soil%var_2d(k,m) ) * ddz_soil_center(k)    &
                             - surf%gamma_w(k,m)                                 &
                             - surf%lambda_w(k-1,m) * ( surf_m_soil%var_2d(k,m)  &
                             - surf_m_soil%var_2d(k-1,m)) * ddz_soil_center(k-1) &
                             + surf%gamma_w(k-1,m) - (root_extr(k)               &
                             * surf%qsws_veg(m) * drho_l_lv)                     &
                             ) * ddz_soil(k)
                   ENDDO
                   tend(nzt_soil) = ( - surf%gamma_w(nzt_soil,m)               &
                                   - surf%lambda_w(nzt_soil-1,m)               &
                                   * ( surf_m_soil%var_2d(nzt_soil,m)          &
                                   - surf_m_soil%var_2d(nzt_soil-1,m))         &
                                   * ddz_soil_center(nzt_soil-1)               &
                                   + surf%gamma_w(nzt_soil-1,m) - (            &
                                   root_extr(nzt_soil)                         &
                                   * surf%qsws_veg(m) * drho_l_lv )            &
                                  ) * ddz_soil(nzt_soil)             

                   surf_m_soil_p%var_2d(nzb_soil:nzt_soil,m) =                 &
                                       surf_m_soil%var_2d(nzb_soil:nzt_soil,m) &
                                         + dt_3d * ( tsc(2) * tend(:)          &
                                         + tsc(3) * surf_tm_soil_m%var_2d(:,m) )   
   
!
!--                Account for dry soils (find a better solution here!)
                   DO  k = nzb_soil, nzt_soil
                      IF ( surf_m_soil_p%var_2d(k,m) < 0.0_wp )  surf_m_soil_p%var_2d(k,m) = 0.0_wp
                   ENDDO
  
!
!--                Calculate m_soil tendencies for the next Runge-Kutta step
                   IF ( timestep_scheme(1:5) == 'runge' )  THEN
                      IF ( intermediate_timestep_count == 1 )  THEN
                         DO  k = nzb_soil, nzt_soil
                            surf_tm_soil_m%var_2d(k,m) = tend(k)
                         ENDDO
                      ELSEIF ( intermediate_timestep_count <                   &
                               intermediate_timestep_count_max )  THEN
                         DO  k = nzb_soil, nzt_soil
                            surf_tm_soil_m%var_2d(k,m) = -9.5625_wp * tend(k)  &
                                                    + 5.3125_wp                &
                                                    * surf_tm_soil_m%var_2d(k,m)
                         ENDDO

                      ENDIF
                      
                   ENDIF
                   
                ENDIF

             ENDIF

          ENDIF

       ENDDO

    END SUBROUTINE lsm_soil_model

 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Swapping of timelevels
!------------------------------------------------------------------------------!
    SUBROUTINE lsm_swap_timelevel ( mod_count )

       IMPLICIT NONE

       INTEGER, INTENT(IN) :: mod_count

#if defined( __nopointer )
!
!--    Horizontal surfaces
       t_surface_h  = t_surface_h_p
       t_soil_h     = t_soil_h_p
       IF ( humidity )  THEN
          m_soil_h    = m_soil_h_p
          m_liq_h  = m_liq_h_p
       ENDIF
!
!--    Vertical surfaces
       t_surface_v  = t_surface_v_p
       t_soil_v     = t_soil_v_p
       IF ( humidity )  THEN
          m_soil_v    = m_soil_v_p
          m_liq_v  = m_liq_v_p
       ENDIF

#else
    
       SELECT CASE ( mod_count )

          CASE ( 0 )
!
!--          Horizontal surfaces
             t_surface_h  => t_surface_h_1; t_surface_h_p  => t_surface_h_2
             t_soil_h     => t_soil_h_1;    t_soil_h_p     => t_soil_h_2
             IF ( humidity )  THEN
                m_soil_h  => m_soil_h_1;    m_soil_h_p     => m_soil_h_2
                m_liq_h   => m_liq_h_1;     m_liq_h_p      => m_liq_h_2
             ENDIF

!
!--          Vertical surfaces
             t_surface_v  => t_surface_v_1; t_surface_v_p  => t_surface_v_2
             t_soil_v     => t_soil_v_1;    t_soil_v_p     => t_soil_v_2
             IF ( humidity )  THEN
                m_soil_v  => m_soil_v_1;    m_soil_v_p     => m_soil_v_2
                m_liq_v   => m_liq_v_1;     m_liq_v_p      => m_liq_v_2

             ENDIF



          CASE ( 1 )
!
!--          Horizontal surfaces
             t_surface_h  => t_surface_h_2; t_surface_h_p  => t_surface_h_1
             t_soil_h     => t_soil_h_2;    t_soil_h_p     => t_soil_h_1
             IF ( humidity )  THEN
                m_soil_h  => m_soil_h_2;    m_soil_h_p     => m_soil_h_1
                m_liq_h   => m_liq_h_2;     m_liq_h_p      => m_liq_h_1

             ENDIF
!
!--          Vertical surfaces
             t_surface_v  => t_surface_v_2; t_surface_v_p  => t_surface_v_1
             t_soil_v     => t_soil_v_2;    t_soil_v_p     => t_soil_v_1
             IF ( humidity )  THEN
                m_soil_v  => m_soil_v_2;    m_soil_v_p     => m_soil_v_1
                m_liq_v   => m_liq_v_2;     m_liq_v_p      => m_liq_v_1
             ENDIF

       END SELECT
#endif

    END SUBROUTINE lsm_swap_timelevel




!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine for averaging 3D data
!------------------------------------------------------------------------------!
SUBROUTINE lsm_3d_data_averaging( mode, variable )
 

    USE control_parameters

    USE indices

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  mode    !< 
    CHARACTER (LEN=*) :: variable !< 

    INTEGER(iwp) ::  i       !< 
    INTEGER(iwp) ::  j       !< 
    INTEGER(iwp) ::  k       !< 
    INTEGER(iwp) ::  m       !< running index

    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( variable ) )

             CASE ( 'c_liq*' )
                IF ( .NOT. ALLOCATED( c_liq_av ) )  THEN
                   ALLOCATE( c_liq_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                c_liq_av = 0.0_wp

             CASE ( 'c_soil*' )
                IF ( .NOT. ALLOCATED( c_soil_av ) )  THEN
                   ALLOCATE( c_soil_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                c_soil_av = 0.0_wp

             CASE ( 'c_veg*' )
                IF ( .NOT. ALLOCATED( c_veg_av ) )  THEN
                   ALLOCATE( c_veg_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                c_veg_av = 0.0_wp

             CASE ( 'lai*' )
                IF ( .NOT. ALLOCATED( lai_av ) )  THEN
                   ALLOCATE( lai_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                lai_av = 0.0_wp

             CASE ( 'm_liq*' )
                IF ( .NOT. ALLOCATED( m_liq_av ) )  THEN
                   ALLOCATE( m_liq_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                m_liq_av = 0.0_wp

             CASE ( 'm_soil' )
                IF ( .NOT. ALLOCATED( m_soil_av ) )  THEN
                   ALLOCATE( m_soil_av(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )
                ENDIF
                m_soil_av = 0.0_wp

             CASE ( 'qsws_liq*' )
                IF ( .NOT. ALLOCATED( qsws_liq_av ) )  THEN
                   ALLOCATE( qsws_liq_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                qsws_liq_av = 0.0_wp

             CASE ( 'qsws_soil*' )
                IF ( .NOT. ALLOCATED( qsws_soil_av ) )  THEN
                   ALLOCATE( qsws_soil_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                qsws_soil_av = 0.0_wp

             CASE ( 'qsws_veg*' )
                IF ( .NOT. ALLOCATED( qsws_veg_av ) )  THEN
                   ALLOCATE( qsws_veg_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                qsws_veg_av = 0.0_wp

             CASE ( 'r_s*' )
                IF ( .NOT. ALLOCATED( r_s_av ) )  THEN
                   ALLOCATE( r_s_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                r_s_av = 0.0_wp

             CASE ( 't_soil' )
                IF ( .NOT. ALLOCATED( t_soil_av ) )  THEN
                   ALLOCATE( t_soil_av(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )
                ENDIF
                t_soil_av = 0.0_wp

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'c_liq*' )
             IF ( ALLOCATED( c_liq_av ) ) THEN 
                DO  m = 1, surf_lsm_h%ns
                   i   = surf_lsm_h%i(m)            
                   j   = surf_lsm_h%j(m)
                   c_liq_av(j,i) = c_liq_av(j,i) + surf_lsm_h%c_liq(m)
                ENDDO
             ENDIF   

          CASE ( 'c_soil*' )
             IF ( ALLOCATED( c_soil_av ) ) THEN 
                DO  m = 1, surf_lsm_h%ns
                   i   = surf_lsm_h%i(m)            
                   j   = surf_lsm_h%j(m)
                   c_soil_av(j,i) = c_soil_av(j,i) + (1.0 - surf_lsm_h%c_veg(m))
                ENDDO
             ENDIF

          CASE ( 'c_veg*' )
             IF ( ALLOCATED( c_veg_av ) ) THEN 
                DO  m = 1, surf_lsm_h%ns
                   i   = surf_lsm_h%i(m)            
                   j   = surf_lsm_h%j(m)
                   c_veg_av(j,i) = c_veg_av(j,i) + surf_lsm_h%c_veg(m)
                ENDDO
             ENDIF

          CASE ( 'lai*' )
             IF ( ALLOCATED( lai_av ) ) THEN 
                DO  m = 1, surf_lsm_h%ns
                   i   = surf_lsm_h%i(m)            
                   j   = surf_lsm_h%j(m)
                   lai_av(j,i) = lai_av(j,i) + surf_lsm_h%lai(m)
                ENDDO
             ENDIF

          CASE ( 'm_liq*' )
             IF ( ALLOCATED( m_liq_av ) ) THEN 
                DO  m = 1, surf_lsm_h%ns
                   i   = surf_lsm_h%i(m)            
                   j   = surf_lsm_h%j(m)
                   m_liq_av(j,i) = m_liq_av(j,i) + m_liq_h%var_1d(m)
                ENDDO
             ENDIF

          CASE ( 'm_soil' )
             IF ( ALLOCATED( m_soil_av ) ) THEN 
                DO  m = 1, surf_lsm_h%ns
                   i   = surf_lsm_h%i(m)            
                   j   = surf_lsm_h%j(m)
                   DO  k = nzb_soil, nzt_soil
                      m_soil_av(k,j,i) = m_soil_av(k,j,i) + m_soil_h%var_2d(k,m)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'qsws_liq*' )
             IF ( ALLOCATED( qsws_liq_av ) ) THEN 
                DO  m = 1, surf_lsm_h%ns
                   i   = surf_lsm_h%i(m)            
                   j   = surf_lsm_h%j(m)
                   qsws_liq_av(j,i) = qsws_liq_av(j,i) +                       &
                                         surf_lsm_h%qsws_liq(m)
                ENDDO
             ENDIF

          CASE ( 'qsws_soil*' )
             IF ( ALLOCATED( qsws_soil_av ) ) THEN 
                DO  m = 1, surf_lsm_h%ns
                   i   = surf_lsm_h%i(m)            
                   j   = surf_lsm_h%j(m)
                   qsws_soil_av(j,i) = qsws_soil_av(j,i) +                     &
                                          surf_lsm_h%qsws_soil(m)
                ENDDO
             ENDIF

          CASE ( 'qsws_veg*' )
             IF ( ALLOCATED(qsws_veg_av ) ) THEN 
                DO  m = 1, surf_lsm_h%ns
                   i   = surf_lsm_h%i(m)            
                   j   = surf_lsm_h%j(m)
                   qsws_veg_av(j,i) = qsws_veg_av(j,i) +                       &
                                         surf_lsm_h%qsws_veg(m)
                ENDDO
             ENDIF

          CASE ( 'r_s*' )
             IF ( ALLOCATED( r_s_av) ) THEN 
                DO  m = 1, surf_lsm_h%ns
                   i   = surf_lsm_h%i(m)            
                   j   = surf_lsm_h%j(m)
                   r_s_av(j,i) = r_s_av(j,i) + surf_lsm_h%r_s(m)
                ENDDO
             ENDIF

          CASE ( 't_soil' )
             IF ( ALLOCATED( t_soil_av ) ) THEN 
                DO  m = 1, surf_lsm_h%ns
                   i   = surf_lsm_h%i(m)            
                   j   = surf_lsm_h%j(m)
                   DO  k = nzb_soil, nzt_soil
                      t_soil_av(k,j,i) = t_soil_av(k,j,i) + t_soil_h%var_2d(k,m)
                   ENDDO
                ENDDO
             ENDIF

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'c_liq*' )
             IF ( ALLOCATED( c_liq_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      c_liq_av(j,i) = c_liq_av(j,i)                            &
                                      / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'c_soil*' )
             IF ( ALLOCATED( c_soil_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      c_soil_av(j,i) = c_soil_av(j,i)                          &
                                       / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'c_veg*' )
             IF ( ALLOCATED( c_veg_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      c_veg_av(j,i) = c_veg_av(j,i)                            &
                                      / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

         CASE ( 'lai*' )
             IF ( ALLOCATED( lai_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      lai_av(j,i) = lai_av(j,i)                                &
                                    / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'm_liq*' )
             IF ( ALLOCATED( m_liq_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      m_liq_av(j,i) = m_liq_av(j,i)                            &
                                      / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'm_soil' )
             IF ( ALLOCATED( m_soil_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_soil, nzt_soil
                         m_soil_av(k,j,i) = m_soil_av(k,j,i)                   &
                                            / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'qsws_liq*' )
             IF ( ALLOCATED( qsws_liq_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      qsws_liq_av(j,i) = qsws_liq_av(j,i)                      &
                                         / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'qsws_soil*' )
             IF ( ALLOCATED( qsws_soil_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      qsws_soil_av(j,i) = qsws_soil_av(j,i)                    &
                                          / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'qsws_veg*' )
             IF ( ALLOCATED( qsws_veg_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      qsws_veg_av(j,i) = qsws_veg_av(j,i)                      &
                                         / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'r_s*' )
             IF ( ALLOCATED( r_s_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      r_s_av(j,i) = r_s_av(j,i)                                & 
                                    / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 't_soil' )
             IF ( ALLOCATED( t_soil_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_soil, nzt_soil
                         t_soil_av(k,j,i) = t_soil_av(k,j,i)                   &
                                            / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
!
!-- 

       END SELECT

    ENDIF

END SUBROUTINE lsm_3d_data_averaging


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!------------------------------------------------------------------------------!
 SUBROUTINE lsm_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )
    
     IMPLICIT NONE

     CHARACTER (LEN=*), INTENT(IN)  ::  var         !< 
     LOGICAL, INTENT(OUT)           ::  found       !< 
     CHARACTER (LEN=*), INTENT(OUT) ::  grid_x      !< 
     CHARACTER (LEN=*), INTENT(OUT) ::  grid_y      !< 
     CHARACTER (LEN=*), INTENT(OUT) ::  grid_z      !< 

     found  = .TRUE.

!
!--  Check for the grid
     SELECT CASE ( TRIM( var ) )

        CASE ( 'm_soil', 't_soil', 'm_soil_xy', 't_soil_xy', 'm_soil_xz',      &
               't_soil_xz', 'm_soil_yz', 't_soil_yz' )
           grid_x = 'x'
           grid_y = 'y'
           grid_z = 'zs'

        CASE DEFAULT
           found  = .FALSE.
           grid_x = 'none'
           grid_y = 'none'
           grid_z = 'none'
     END SELECT

 END SUBROUTINE lsm_define_netcdf_grid

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 3D output variables
!------------------------------------------------------------------------------!
 SUBROUTINE lsm_data_output_2d( av, variable, found, grid, mode, local_pf,     &
                                two_d, nzb_do, nzt_do )
 
    USE indices


    IMPLICIT NONE

    CHARACTER (LEN=*) ::  grid     !< 
    CHARACTER (LEN=*) ::  mode     !< 
    CHARACTER (LEN=*) ::  variable !< 

    INTEGER(iwp) ::  av      !< 
    INTEGER(iwp) ::  i       !< running index 
    INTEGER(iwp) ::  j       !< running index
    INTEGER(iwp) ::  k       !< running index
    INTEGER(iwp) ::  m       !< running index
    INTEGER(iwp) ::  nzb_do  !< 
    INTEGER(iwp) ::  nzt_do  !< 

    LOGICAL      ::  found !< 
    LOGICAL      ::  two_d !< flag parameter that indicates 2D variables (horizontal cross sections)

    REAL(wp) ::  fill_value = -999.0_wp    !< value for the _FillValue attribute

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf !< 


    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )
!
!--    Before data is transfered to local_pf, transfer is it 2D dummy variable and exchange ghost points therein. 
!--    However, at this point this is only required for instantaneous arrays, time-averaged quantities are already exchanged. 
       CASE ( 'c_liq*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  m = 1, surf_lsm_h%ns
                i                   = surf_lsm_h%i(m)            
                j                   = surf_lsm_h%j(m)
                local_pf(i,j,nzb+1) = surf_lsm_h%c_liq(m) * surf_lsm_h%c_veg(m)
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( c_liq_av ) ) THEN
               ALLOCATE( c_liq_av(nysg:nyng,nxlg:nxrg) )
               c_liq_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = c_liq_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'c_soil*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  m = 1, surf_lsm_h%ns
                i                   = surf_lsm_h%i(m)            
                j                   = surf_lsm_h%j(m)
                local_pf(i,j,nzb+1) = 1.0_wp - surf_lsm_h%c_veg(m)
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( c_soil_av ) ) THEN
               ALLOCATE( c_soil_av(nysg:nyng,nxlg:nxrg) )
               c_soil_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = c_soil_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'c_veg*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  m = 1, surf_lsm_h%ns
                i                   = surf_lsm_h%i(m)            
                j                   = surf_lsm_h%j(m)
                local_pf(i,j,nzb+1) = surf_lsm_h%c_veg(m)
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( c_veg_av ) ) THEN
               ALLOCATE( c_veg_av(nysg:nyng,nxlg:nxrg) )
               c_veg_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = c_veg_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'lai*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  m = 1, surf_lsm_h%ns
                i                   = surf_lsm_h%i(m)            
                j                   = surf_lsm_h%j(m)
                local_pf(i,j,nzb+1) = surf_lsm_h%lai(m)
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( lai_av ) ) THEN
               ALLOCATE( lai_av(nysg:nyng,nxlg:nxrg) )
               lai_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = lai_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'm_liq*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  m = 1, surf_lsm_h%ns
                i                   = surf_lsm_h%i(m)            
                j                   = surf_lsm_h%j(m)
                local_pf(i,j,nzb+1) = m_liq_h%var_1d(m)
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( m_liq_av ) ) THEN
               ALLOCATE( m_liq_av(nysg:nyng,nxlg:nxrg) )
               m_liq_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = m_liq_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'm_soil_xy', 'm_soil_xz', 'm_soil_yz' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf_lsm_h%ns
                i   = surf_lsm_h%i(m)            
                j   = surf_lsm_h%j(m)
                DO k = nzb_soil, nzt_soil
                   local_pf(i,j,k) = m_soil_h%var_2d(k,m)
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( m_soil_av ) ) THEN
               ALLOCATE( m_soil_av(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )
               m_soil_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO k = nzb_soil, nzt_soil
                      local_pf(i,j,k) = m_soil_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          nzb_do = nzb_soil
          nzt_do = nzt_soil

          IF ( mode == 'xy' ) grid = 'zs'

       CASE ( 'qsws_liq*_xy' )        ! 2d-array
          IF ( av == 0 ) THEN
             DO  m = 1, surf_lsm_h%ns
                i                   = surf_lsm_h%i(m)            
                j                   = surf_lsm_h%j(m)
                local_pf(i,j,nzb+1) = surf_lsm_h%qsws_liq(m)
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( qsws_liq_av ) ) THEN
               ALLOCATE( qsws_liq_av(nysg:nyng,nxlg:nxrg) )
               qsws_liq_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn 
                   local_pf(i,j,nzb+1) =  qsws_liq_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'qsws_soil*_xy' )        ! 2d-array
          IF ( av == 0 ) THEN
             DO  m = 1, surf_lsm_h%ns
                i                   = surf_lsm_h%i(m)            
                j                   = surf_lsm_h%j(m)
                local_pf(i,j,nzb+1) =  surf_lsm_h%qsws_soil(m)
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( qsws_soil_av ) ) THEN
               ALLOCATE( qsws_soil_av(nysg:nyng,nxlg:nxrg) )
               qsws_soil_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn 
                   local_pf(i,j,nzb+1) =  qsws_soil_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'qsws_veg*_xy' )        ! 2d-array
          IF ( av == 0 ) THEN
             DO  m = 1, surf_lsm_h%ns
                i                   = surf_lsm_h%i(m)            
                j                   = surf_lsm_h%j(m)
                local_pf(i,j,nzb+1) =  surf_lsm_h%qsws_veg(m)
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( qsws_veg_av ) ) THEN
               ALLOCATE( qsws_veg_av(nysg:nyng,nxlg:nxrg) )
               qsws_veg_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn 
                   local_pf(i,j,nzb+1) =  qsws_veg_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'


       CASE ( 'r_s*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  m = 1, surf_lsm_h%ns
                i                   = surf_lsm_h%i(m)            
                j                   = surf_lsm_h%j(m)
                local_pf(i,j,nzb+1) = surf_lsm_h%r_s(m)
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( r_s_av ) ) THEN
               ALLOCATE( r_s_av(nysg:nyng,nxlg:nxrg) )
               r_s_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = r_s_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 't_soil_xy', 't_soil_xz', 't_soil_yz' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf_lsm_h%ns
                i   = surf_lsm_h%i(m)            
                j   = surf_lsm_h%j(m)
                DO k = nzb_soil, nzt_soil
                   local_pf(i,j,k) = t_soil_h%var_2d(k,m)
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( t_soil_av ) ) THEN
               ALLOCATE( t_soil_av(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )
               t_soil_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO k = nzb_soil, nzt_soil
                      local_pf(i,j,k) = t_soil_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          nzb_do = nzb_soil
          nzt_do = nzt_soil

          IF ( mode == 'xy' )  grid = 'zs'

       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT
 
 END SUBROUTINE lsm_data_output_2d


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 3D output variables
!------------------------------------------------------------------------------!
 SUBROUTINE lsm_data_output_3d( av, variable, found, local_pf )
 

    USE indices


    IMPLICIT NONE

    CHARACTER (LEN=*) ::  variable !< 

    INTEGER(iwp) ::  av    !< 
    INTEGER(iwp) ::  i     !< 
    INTEGER(iwp) ::  j     !< 
    INTEGER(iwp) ::  k     !< 
    INTEGER(iwp) ::  m     !< running index

    LOGICAL      ::  found !< 

    REAL(wp) ::  fill_value = -999.0_wp    !< value for the _FillValue attribute

    REAL(sp), DIMENSION(nxl:nxr,nys:nyn,nzb_soil:nzt_soil) ::  local_pf !< 


    found = .TRUE.


    SELECT CASE ( TRIM( variable ) )
!
!--   Requires 3D exchange

      CASE ( 'm_soil' )

         IF ( av == 0 )  THEN
            DO  m = 1, surf_lsm_h%ns
                i   = surf_lsm_h%i(m)            
                j   = surf_lsm_h%j(m)
                DO  k = nzb_soil, nzt_soil
                   local_pf(i,j,k) = m_soil_h%var_2d(k,m)
                ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( m_soil_av ) ) THEN
               ALLOCATE( m_soil_av(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )
               m_soil_av = REAL( fill_value, KIND = wp )
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_soil, nzt_soil
                     local_pf(i,j,k) = m_soil_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 't_soil' )

         IF ( av == 0 )  THEN
            DO  m = 1, surf_lsm_h%ns
               i   = surf_lsm_h%i(m)            
               j   = surf_lsm_h%j(m)
               DO  k = nzb_soil, nzt_soil
                  local_pf(i,j,k) = t_soil_h%var_2d(k,m)
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( t_soil_av ) ) THEN
               ALLOCATE( t_soil_av(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )
               t_soil_av = REAL( fill_value, KIND = wp )
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_soil, nzt_soil
                     local_pf(i,j,k) = t_soil_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF


       CASE DEFAULT
          found = .FALSE.

    END SELECT


 END SUBROUTINE lsm_data_output_3d


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Write restart data for land surface model. It is necessary to write 
!> start_index and end_index several times.
!------------------------------------------------------------------------------!
 SUBROUTINE lsm_wrd_local
        

    IMPLICIT NONE

    CHARACTER (LEN=1) ::  dum    !< dummy to create correct string for creating variable string
    INTEGER(iwp)      ::  l      !< index variable for surface orientation

    CALL wrd_write_string( 'ns_h_on_file_lsm' )
    WRITE ( 14 )  surf_lsm_h%ns

    CALL wrd_write_string( 'ns_v_on_file_lsm' )
    WRITE ( 14 )  surf_lsm_v(0:3)%ns


    IF ( ALLOCATED( c_liq_av ) )  THEN
       CALL wrd_write_string( 'c_liq_av' )
       WRITE ( 14 )  c_liq_av
    ENDIF

    IF ( ALLOCATED( c_soil_av ) )  THEN
       CALL wrd_write_string( 'c_soil_av' )
       WRITE ( 14 )  c_soil_av
    ENDIF

    IF ( ALLOCATED( c_veg_av ) )  THEN
       CALL wrd_write_string( 'c_veg_av' )
       WRITE ( 14 )  c_veg_av
    ENDIF

    IF ( ALLOCATED( lai_av ) )  THEN
       CALL wrd_write_string( 'lai_av' )
       WRITE ( 14 )  lai_av
    ENDIF

    IF ( ALLOCATED( m_liq_av ) )  THEN
       CALL wrd_write_string( 'm_liq_av' )
       WRITE ( 14 )  m_liq_av
    ENDIF

    IF ( ALLOCATED( m_soil_av ) )  THEN
       CALL wrd_write_string( 'm_soil_av' )
       WRITE ( 14 )  m_soil_av
    ENDIF

    IF ( ALLOCATED( qsws_liq_av ) )  THEN
       CALL wrd_write_string( 'qsws_liq_av' )
       WRITE ( 14 )  qsws_liq_av
    ENDIF

    IF ( ALLOCATED( qsws_soil_av ) )  THEN
       CALL wrd_write_string( 'qsws_soil_av' )
       WRITE ( 14 )  qsws_soil_av
    ENDIF

    IF ( ALLOCATED( qsws_veg_av ) )  THEN
       CALL wrd_write_string( 'qsws_veg_av' )
       WRITE ( 14 )  qsws_veg_av
    ENDIF

    IF ( ALLOCATED( t_soil_av ) )  THEN
       CALL wrd_write_string( 't_soil_av' )
       WRITE ( 14 )  t_soil_av
    ENDIF

    CALL wrd_write_string( 'lsm_start_index_h' )
    WRITE ( 14 )  surf_lsm_h%start_index

    CALL wrd_write_string( 'lsm_end_index_h' )
    WRITE ( 14 )  surf_lsm_h%end_index

    CALL wrd_write_string( 't_soil_h' )
    WRITE ( 14 )  t_soil_h%var_2d
       

       
    DO  l = 0, 3

       CALL wrd_write_string( 'lsm_start_index_v' )
       WRITE ( 14 )  surf_lsm_v(l)%start_index

       CALL wrd_write_string( 'lsm_end_index_v' )
       WRITE ( 14 )  surf_lsm_v(l)%end_index

       WRITE( dum, '(I1)')  l   

       CALL wrd_write_string( 't_soil_v(' // dum // ')' )
       WRITE ( 14 )  t_soil_v(l)%var_2d
             
    ENDDO

    CALL wrd_write_string( 'lsm_start_index_h' )
    WRITE ( 14 )  surf_lsm_h%start_index

    CALL wrd_write_string( 'lsm_end_index_h' )
    WRITE ( 14 )  surf_lsm_h%end_index

    CALL wrd_write_string( 'm_soil_h' )
    WRITE ( 14 )  m_soil_h%var_2d

    DO  l = 0, 3

       CALL wrd_write_string( 'lsm_start_index_v' )
       WRITE ( 14 )  surf_lsm_v(l)%start_index

       CALL wrd_write_string( 'lsm_end_index_v' )
       WRITE ( 14 )  surf_lsm_v(l)%end_index

       WRITE( dum, '(I1)')  l   

       CALL wrd_write_string( 'm_soil_v(' // dum // ')' )
       WRITE ( 14 )  m_soil_v(l)%var_2d 
      
    ENDDO

    CALL wrd_write_string( 'lsm_start_index_h' )
    WRITE ( 14 )  surf_lsm_h%start_index

    CALL wrd_write_string( 'lsm_end_index_h' )
    WRITE ( 14 )  surf_lsm_h%end_index

    CALL wrd_write_string( 'm_liq_h' )
    WRITE ( 14 )  m_liq_h%var_1d
       
    DO  l = 0, 3

       CALL wrd_write_string( 'lsm_start_index_v' )
       WRITE ( 14 )  surf_lsm_v(l)%start_index

       CALL wrd_write_string( 'lsm_end_index_v' )
       WRITE ( 14 )  surf_lsm_v(l)%end_index

       WRITE( dum, '(I1)')  l   

       CALL wrd_write_string( 'm_liq_v(' // dum // ')' )
       WRITE ( 14 )  m_liq_v(l)%var_1d     
                
    ENDDO

    CALL wrd_write_string( 'lsm_start_index_h' )
    WRITE ( 14 )  surf_lsm_h%start_index

    CALL wrd_write_string( 'lsm_end_index_h' )
    WRITE ( 14 )  surf_lsm_h%end_index

    CALL wrd_write_string( 't_surface_h' )
    WRITE ( 14 )  t_surface_h%var_1d

    DO  l = 0, 3

       CALL wrd_write_string( 'lsm_start_index_v' )
       WRITE ( 14 )  surf_lsm_v(l)%start_index

       CALL wrd_write_string( 'lsm_end_index_v' )
       WRITE ( 14 )  surf_lsm_v(l)%end_index

       WRITE( dum, '(I1)')  l   

       CALL wrd_write_string( 't_surface_v(' // dum // ')' )
       WRITE ( 14 )  t_surface_v(l)%var_1d     
       
    ENDDO


 END SUBROUTINE lsm_wrd_local


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Soubroutine reads lsm data from restart file(s)
!------------------------------------------------------------------------------!
SUBROUTINE lsm_rrd_local( i, k, nxlf, nxlc, nxl_on_file, nxrf, nxrc,           &
                          nxr_on_file, nynf, nync, nyn_on_file, nysf, nysc,    &
                          nys_on_file, tmp_2d, found )
 

    USE control_parameters
        
    USE indices
    
    USE pegrid


    IMPLICIT NONE

    INTEGER(iwp) ::  i                 !< 
    INTEGER(iwp) ::  k                 !< 
    INTEGER(iwp) ::  l                 !< running index surface orientation
    INTEGER(iwp) ::  ns_h_on_file_lsm  !< number of horizontal surface elements (natural type) on file
    INTEGER(iwp) ::  nxlc              !< 
    INTEGER(iwp) ::  nxlf              !< 
    INTEGER(iwp) ::  nxl_on_file       !< index of left boundary on former local domain
    INTEGER(iwp) ::  nxrc              !< 
    INTEGER(iwp) ::  nxrf              !< 
    INTEGER(iwp) ::  nxr_on_file       !< index of right boundary on former local domain
    INTEGER(iwp) ::  nync              !< 
    INTEGER(iwp) ::  nynf              !< 
    INTEGER(iwp) ::  nyn_on_file       !< index of north boundary on former local domain
    INTEGER(iwp) ::  nysc              !< 
    INTEGER(iwp) ::  nysf              !< 
    INTEGER(iwp) ::  nys_on_file       !< index of south boundary on former local domain

    INTEGER(iwp) ::  ns_v_on_file_lsm(0:3) !< number of vertical surface elements (natural type) on file

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  start_index_on_file 
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  end_index_on_file

    LOGICAL, INTENT(OUT)  :: found

    REAL(wp), DIMENSION(nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_2d   !< 

    REAL(wp), DIMENSION(nzb_soil:nzt_soil+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_3d   !< 

    REAL(wp), DIMENSION(nzb_soil:nzt_soil,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_3d2   !< 

    TYPE(surf_type_lsm), SAVE :: tmp_walltype_h_1d   !< temporary 1D array containing the respective surface variable stored on file, horizontal surfaces
    TYPE(surf_type_lsm), SAVE :: tmp_walltype_h_2d   !< temporary 2D array containing the respective surface variable stored on file, horizontal surfaces
    TYPE(surf_type_lsm), SAVE :: tmp_walltype_h_2d2  !< temporary 2D array containing the respective surface variable stored on file, horizontal surfaces 

    TYPE(surf_type_lsm), DIMENSION(0:3), SAVE :: tmp_walltype_v_1d   !< temporary 1D array containing the respective surface variable stored on file, vertical surfaces 
    TYPE(surf_type_lsm), DIMENSION(0:3), SAVE :: tmp_walltype_v_2d   !< temporary 2D array containing the respective surface variable stored on file, vertical surfaces 
    TYPE(surf_type_lsm), DIMENSION(0:3), SAVE :: tmp_walltype_v_2d2  !< temporary 2D array containing the respective surface variable stored on file, vertical surfaces


    found = .TRUE.


       SELECT CASE ( restart_string(1:length) )

           CASE ( 'ns_h_on_file_lsm' )
              IF ( k == 1 )  THEN  
                 READ ( 13 ) ns_h_on_file_lsm

                 IF ( ALLOCATED( tmp_walltype_h_1d%var_1d ) )                  &
                    DEALLOCATE( tmp_walltype_h_1d%var_1d )
                 IF ( ALLOCATED( tmp_walltype_h_2d%var_2d ) )                  &   
                    DEALLOCATE( tmp_walltype_h_2d%var_2d )
                 IF ( ALLOCATED( tmp_walltype_h_2d2%var_2d ) )                 &
                    DEALLOCATE( tmp_walltype_h_2d2%var_2d ) 

!
!--              Allocate temporary arrays to store surface data
                 ALLOCATE( tmp_walltype_h_1d%var_1d(1:ns_h_on_file_lsm) )
                 ALLOCATE( tmp_walltype_h_2d%var_2d(nzb_soil:nzt_soil+1,       &
                                                    1:ns_h_on_file_lsm) )
                 ALLOCATE( tmp_walltype_h_2d2%var_2d(nzb_soil:nzt_soil,        &
                           1:ns_h_on_file_lsm)  )

              ENDIF

           CASE ( 'ns_v_on_file_lsm' )
              IF ( k == 1 )  THEN
                 READ ( 13 ) ns_v_on_file_lsm

                 DO  l = 0, 3
                    IF ( ALLOCATED( tmp_walltype_v_1d(l)%var_1d ) )            &
                       DEALLOCATE( tmp_walltype_v_1d(l)%var_1d )
                    IF ( ALLOCATED( tmp_walltype_v_2d(l)%var_2d ) )            &
                       DEALLOCATE( tmp_walltype_v_2d(l)%var_2d )
                    IF ( ALLOCATED( tmp_walltype_v_2d2(l)%var_2d ) )           &
                       DEALLOCATE( tmp_walltype_v_2d2(l)%var_2d )
                 ENDDO

!
!--              Allocate temporary arrays to store surface data
                 DO  l = 0, 3
                    ALLOCATE( tmp_walltype_v_1d(l)                             &
                                 %var_1d(1:ns_v_on_file_lsm(l)) )
                    ALLOCATE( tmp_walltype_v_2d(l)                             &
                                 %var_2d(nzb_soil:nzt_soil+1,                  &
                                         1:ns_v_on_file_lsm(l)) )
                    ALLOCATE( tmp_walltype_v_2d2(l)                            &
                                 %var_2d(nzb_soil:nzt_soil,                    &
                                         1:ns_v_on_file_lsm(l))  )
                 ENDDO

              ENDIF


           CASE ( 'c_liq_av' )
              IF ( .NOT. ALLOCATED( c_liq_av ) )  THEN
                 ALLOCATE( c_liq_av(nysg:nyng,nxlg:nxrg) )
              ENDIF
              IF ( k == 1 )  READ ( 13 )  tmp_2d
              c_liq_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =              &
                 tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

           CASE ( 'c_soil_av' )
              IF ( .NOT. ALLOCATED( c_soil_av ) )  THEN
                 ALLOCATE( c_soil_av(nysg:nyng,nxlg:nxrg) )
              ENDIF
              IF ( k == 1 )  READ ( 13 )  tmp_2d
              c_soil_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =             &
                 tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

           CASE ( 'c_veg_av' )
              IF ( .NOT. ALLOCATED( c_veg_av ) )  THEN
                 ALLOCATE( c_veg_av(nysg:nyng,nxlg:nxrg) )
              ENDIF
              IF ( k == 1 )  READ ( 13 )  tmp_2d
              c_veg_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =              &
                 tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

           CASE ( 'lai_av' )
              IF ( .NOT. ALLOCATED( lai_av ) )  THEN
                 ALLOCATE( lai_av(nysg:nyng,nxlg:nxrg) )
              ENDIF
              IF ( k == 1 )  READ ( 13 )  tmp_2d
              lai_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                &
                 tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

           CASE ( 'm_liq_av' )
              IF ( .NOT. ALLOCATED( m_liq_av ) )  THEN
                 ALLOCATE( m_liq_av(nysg:nyng,nxlg:nxrg) )
              ENDIF
              IF ( k == 1 )  READ ( 13 )  tmp_2d
              m_liq_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =              &
                 tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

           CASE ( 'm_soil_av' )
              IF ( .NOT. ALLOCATED( m_soil_av ) )  THEN
                 ALLOCATE( m_soil_av(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )
              ENDIF
              IF ( k == 1 )  READ ( 13 )  tmp_3d2(:,:,:)
              m_soil_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                 tmp_3d2(nzb_soil:nzt_soil,nysf                                &
                         -nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

           CASE ( 'qsws_liq_av' )
              IF ( .NOT. ALLOCATED( qsws_liq_av ) )  THEN
                 ALLOCATE( qsws_liq_av(nysg:nyng,nxlg:nxrg) )
              ENDIF  
              IF ( k == 1 )  READ ( 13 )  tmp_2d
              qsws_liq_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =          &
                 tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
           CASE ( 'qsws_soil_av' )
              IF ( .NOT. ALLOCATED( qsws_soil_av ) )  THEN
                 ALLOCATE( qsws_soil_av(nysg:nyng,nxlg:nxrg) )
              ENDIF  
              IF ( k == 1 )  READ ( 13 )  tmp_2d
              qsws_soil_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =         &
                 tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

           CASE ( 'qsws_veg_av' )
              IF ( .NOT. ALLOCATED( qsws_veg_av ) )  THEN
                 ALLOCATE( qsws_veg_av(nysg:nyng,nxlg:nxrg) )
              ENDIF  
              IF ( k == 1 )  READ ( 13 )  tmp_2d
              qsws_veg_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =          &
                 tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

           CASE ( 't_soil_av' )
              IF ( .NOT. ALLOCATED( t_soil_av ) )  THEN
                 ALLOCATE( t_soil_av(nzb_soil:nzt_soil,nysg:nyng,nxlg:nxrg) )
              ENDIF
              IF ( k == 1 )  READ ( 13 )  tmp_3d2(:,:,:)
              t_soil_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                 tmp_3d2(:,nysf-nbgp:nynf+nbgp,                                &
                         nxlf-nbgp:nxrf+nbgp)

           CASE ( 'lsm_start_index_h', 'lsm_start_index_v'  )   
                IF ( k == 1 )  THEN

                   IF ( ALLOCATED( start_index_on_file ) )                     &
                      DEALLOCATE( start_index_on_file )

                   ALLOCATE ( start_index_on_file(nys_on_file:nyn_on_file,     &
                   nxl_on_file:nxr_on_file) )

                   READ ( 13 )  start_index_on_file

                ENDIF
                 
           CASE ( 'lsm_end_index_h', 'lsm_end_index_v' )   
                IF ( k == 1 )  THEN

                   IF ( ALLOCATED( end_index_on_file ) )                       &
                      DEALLOCATE( end_index_on_file )

                   ALLOCATE ( end_index_on_file(nys_on_file:nyn_on_file,       &
                      nxl_on_file:nxr_on_file) )

                   READ ( 13 )  end_index_on_file

                ENDIF
           
           CASE ( 't_soil_h' )
            
              IF ( k == 1 )  THEN
                 IF ( .NOT.  ALLOCATED( t_soil_h%var_2d ) )                    &
                    ALLOCATE( t_soil_h%var_2d(nzb_soil:nzt_soil+1,             &
                                              1:surf_lsm_h%ns) )
                 READ ( 13 )  tmp_walltype_h_2d%var_2d
              ENDIF
              CALL surface_restore_elements(                                   &
                                         t_soil_h%var_2d,                      &
                                         tmp_walltype_h_2d%var_2d,             &
                                         surf_lsm_h%start_index,               & 
                                         start_index_on_file,                  &
                                         end_index_on_file,                    &
                                         nxlc, nysc,                           &
                                         nxlf, nxrf, nysf, nynf,               &
                                         nys_on_file, nyn_on_file,             &
                                         nxl_on_file,nxr_on_file )

           CASE ( 't_soil_v(0)' )
            
              IF ( k == 1 )  THEN
                 IF ( .NOT.  ALLOCATED( t_soil_v(0)%var_2d ) )                 &
                    ALLOCATE( t_soil_v(0)%var_2d(nzb_soil:nzt_soil+1,          &
                                                 1:surf_lsm_v(0)%ns) )
                 READ ( 13 )  tmp_walltype_v_2d(0)%var_2d
              ENDIF
              CALL surface_restore_elements(                                   &
                                      t_soil_v(0)%var_2d,                      &
                                      tmp_walltype_v_2d(0)%var_2d,             &
                                      surf_lsm_v(0)%start_index,               &  
                                      start_index_on_file,                     &
                                      end_index_on_file,                       &
                                      nxlc, nysc,                              &
                                      nxlf, nxrf, nysf, nynf,                  &
                                      nys_on_file, nyn_on_file,                &
                                      nxl_on_file,nxr_on_file )

           CASE ( 't_soil_v(1)' )
            
              IF ( k == 1 )  THEN
                 IF ( .NOT.  ALLOCATED( t_soil_v(1)%var_2d ) )                 &
                    ALLOCATE( t_soil_v(1)%var_2d(nzb_soil:nzt_soil+1,          &
                                                 1:surf_lsm_v(1)%ns) )
                 READ ( 13 )  tmp_walltype_v_2d(1)%var_2d
              ENDIF
              CALL surface_restore_elements(                                   &
                                      t_soil_v(1)%var_2d,                      &
                                      tmp_walltype_v_2d(1)%var_2d,             &
                                      surf_lsm_v(1)%start_index,               &   
                                      start_index_on_file,                     &
                                      end_index_on_file,                       &
                                      nxlc, nysc,                              &
                                      nxlf, nxrf, nysf, nynf,                  &
                                      nys_on_file, nyn_on_file,                &
                                      nxl_on_file,nxr_on_file )

           CASE ( 't_soil_v(2)' )
            
              IF ( k == 1 )  THEN
                 IF ( .NOT.  ALLOCATED( t_soil_v(2)%var_2d ) )                 &
                    ALLOCATE( t_soil_v(2)%var_2d(nzb_soil:nzt_soil+1,          &
                                                 1:surf_lsm_v(2)%ns) )
                 READ ( 13 )  tmp_walltype_v_2d(2)%var_2d
              ENDIF
              CALL surface_restore_elements(                                   &
                                      t_soil_v(2)%var_2d,                      &
                                      tmp_walltype_v_2d(2)%var_2d,             &
                                      surf_lsm_v(2)%start_index,               & 
                                      start_index_on_file,                     &
                                      end_index_on_file,                       &
                                      nxlc, nysc,                              &
                                      nxlf, nxrf, nysf, nynf,                  &
                                      nys_on_file, nyn_on_file,                &
                                      nxl_on_file,nxr_on_file )

           CASE ( 't_soil_v(3)' )
            
              IF ( k == 1 )  THEN
                 IF ( .NOT.  ALLOCATED( t_soil_v(3)%var_2d ) )                 &
                    ALLOCATE( t_soil_v(1)%var_2d(nzb_soil:nzt_soil+1,          &
                                                 1:surf_lsm_v(3)%ns) )
                 READ ( 13 )  tmp_walltype_v_2d(3)%var_2d
              ENDIF
              CALL surface_restore_elements(                                   &
                                      t_soil_v(3)%var_2d,                      &
                                      tmp_walltype_v_2d(3)%var_2d,             &
                                      surf_lsm_v(3)%start_index,               & 
                                      start_index_on_file,                     &
                                      end_index_on_file,                       &
                                      nxlc, nysc,                              &
                                      nxlf, nxrf, nysf, nynf,                  &
                                      nys_on_file, nyn_on_file,                &
                                      nxl_on_file,nxr_on_file )

           CASE ( 'm_soil_h' )
            
              IF ( k == 1 )  THEN
                 IF ( .NOT.  ALLOCATED( m_soil_h%var_2d ) )                    &
                    ALLOCATE( m_soil_h%var_2d(nzb_soil:nzt_soil+1,             &
                                              1:surf_lsm_h%ns) )
                 READ ( 13 )  tmp_walltype_h_2d2%var_2d
              ENDIF
              CALL surface_restore_elements(                                   &
                                        m_soil_h%var_2d,                       &
                                        tmp_walltype_h_2d2%var_2d,             &
                                        surf_lsm_h%start_index,                &  
                                        start_index_on_file,                   &
                                        end_index_on_file,                     &
                                        nxlc, nysc,                            &
                                        nxlf, nxrf, nysf, nynf,                &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file,nxr_on_file )

           CASE ( 'm_soil_v(0)' )
            
              IF ( k == 1 )  THEN
                 IF ( .NOT.  ALLOCATED( m_soil_v(0)%var_2d ) )                 &
                    ALLOCATE( m_soil_v(0)%var_2d(nzb_soil:nzt_soil+1,          &
                                                 1:surf_lsm_v(0)%ns) )
                 READ ( 13 )  tmp_walltype_v_2d2(0)%var_2d
              ENDIF
              CALL surface_restore_elements(                                   &
                                     m_soil_v(0)%var_2d,                       & 
                                     tmp_walltype_v_2d2(0)%var_2d,             &
                                     surf_lsm_v(0)%start_index,                & 
                                     start_index_on_file,                      &
                                     end_index_on_file,                        &
                                     nxlc, nysc,                               &
                                     nxlf, nxrf, nysf, nynf,                   &
                                     nys_on_file, nyn_on_file,                 &
                                     nxl_on_file,nxr_on_file )

           CASE ( 'm_soil_v(1)' )
            
              IF ( k == 1 )  THEN
                 IF ( .NOT.  ALLOCATED( m_soil_v(1)%var_2d ) )                 &
                    ALLOCATE( m_soil_v(1)%var_2d(nzb_soil:nzt_soil+1,          &
                                                 1:surf_lsm_v(1)%ns) )
                 READ ( 13 )  tmp_walltype_v_2d2(1)%var_2d
              ENDIF
              CALL surface_restore_elements(                                   &
                                     m_soil_v(1)%var_2d,                       &    
                                     tmp_walltype_v_2d2(1)%var_2d,             &
                                     surf_lsm_v(1)%start_index,                &  
                                     start_index_on_file,                      &
                                     end_index_on_file,                        &
                                     nxlc, nysc,                               &
                                     nxlf, nxrf, nysf, nynf,                   &
                                     nys_on_file, nyn_on_file,                 &
                                     nxl_on_file,nxr_on_file )


           CASE ( 'm_soil_v(2)' )
            
              IF ( k == 1 )  THEN
                 IF ( .NOT.  ALLOCATED( m_soil_v(2)%var_2d ) )                 &
                    ALLOCATE( m_soil_v(2)%var_2d(nzb_soil:nzt_soil+1,          &
                                                 1:surf_lsm_v(2)%ns) )
                 READ ( 13 )  tmp_walltype_v_2d2(2)%var_2d
              ENDIF
              CALL surface_restore_elements(                                   &
                                     m_soil_v(2)%var_2d,                       & 
                                     tmp_walltype_v_2d2(2)%var_2d,             &
                                     surf_lsm_v(2)%start_index,                &   
                                     start_index_on_file,                      &
                                     end_index_on_file,                        &
                                     nxlc, nysc,                               &
                                     nxlf, nxrf, nysf, nynf,                   &
                                     nys_on_file, nyn_on_file,                 &
                                     nxl_on_file,nxr_on_file )


           CASE ( 'm_soil_v(3)' )
            
              IF ( k == 1 )  THEN
                 IF ( .NOT.  ALLOCATED( m_soil_v(3)%var_2d ) )                 &
                    ALLOCATE( m_soil_v(1)%var_2d(nzb_soil:nzt_soil+1,          &
                                                 1:surf_lsm_v(3)%ns) )
                 READ ( 13 )  tmp_walltype_v_2d2(3)%var_2d
              ENDIF
              CALL surface_restore_elements(                                   &
                                     m_soil_v(3)%var_2d,                       & 
                                     tmp_walltype_v_2d2(3)%var_2d,             &
                                     surf_lsm_v(3)%start_index,                &  
                                     start_index_on_file,                      &
                                     end_index_on_file,                        &
                                     nxlc, nysc,                               &
                                     nxlf, nxrf, nysf, nynf,                   &
                                     nys_on_file, nyn_on_file,                 &
                                     nxl_on_file,nxr_on_file )


           CASE ( 'm_liq_h' )
            
              IF ( k == 1 )  THEN
                 IF ( .NOT.  ALLOCATED( m_liq_h%var_1d ) )                     &
                    ALLOCATE( m_liq_h%var_1d(1:surf_lsm_h%ns) )
                 READ ( 13 )  tmp_walltype_h_1d%var_1d
              ENDIF
              CALL surface_restore_elements(                                   &
                                         m_liq_h%var_1d,                       &
                                         tmp_walltype_h_1d%var_1d,             &
                                         surf_lsm_h%start_index,               &  
                                         start_index_on_file,                  &
                                         end_index_on_file,                    &
                                         nxlc, nysc,                           &
                                         nxlf, nxrf, nysf, nynf,               &
                                         nys_on_file, nyn_on_file,             &
                                         nxl_on_file,nxr_on_file )


           CASE ( 'm_liq_v(0)' )
            
              IF ( k == 1 )  THEN
                 IF ( .NOT.  ALLOCATED( m_liq_v(0)%var_1d ) )                  &
                    ALLOCATE( m_liq_v(0)%var_1d(1:surf_lsm_v(0)%ns) )
                 READ ( 13 )  tmp_walltype_v_1d(0)%var_1d
              ENDIF
              CALL surface_restore_elements(                                   &
                                      m_liq_v(0)%var_1d,                       &
                                      tmp_walltype_v_1d(0)%var_1d,             &
                                      surf_lsm_v(0)%start_index,               & 
                                      start_index_on_file,                     &
                                      end_index_on_file,                       &
                                      nxlc, nysc,                              &
                                      nxlf, nxrf, nysf, nynf,                  &
                                      nys_on_file, nyn_on_file,                &
                                      nxl_on_file,nxr_on_file )


           CASE ( 'm_liq_v(1)' )
            
              IF ( k == 1 )  THEN
                 IF ( .NOT.  ALLOCATED( m_liq_v(1)%var_1d ) )                  &
                    ALLOCATE( m_liq_v(1)%var_1d(1:surf_lsm_v(1)%ns) )
                 READ ( 13 )  tmp_walltype_v_1d(1)%var_1d
              ENDIF
              CALL surface_restore_elements(                                   &
                                      m_liq_v(1)%var_1d,                       &
                                      tmp_walltype_v_1d(1)%var_1d,             &
                                      surf_lsm_v(1)%start_index,               & 
                                      start_index_on_file,                     &
                                      end_index_on_file,                       &
                                      nxlc, nysc,                              &
                                      nxlf, nxrf, nysf, nynf,                  &
                                      nys_on_file, nyn_on_file,                &
                                      nxl_on_file,nxr_on_file )


           CASE ( 'm_liq_v(2)' )
            
              IF ( k == 1 )  THEN
                 IF ( .NOT.  ALLOCATED( m_liq_v(2)%var_1d ) )                  &
                    ALLOCATE( m_liq_v(2)%var_1d(1:surf_lsm_v(2)%ns) )
                 READ ( 13 )  tmp_walltype_v_1d(2)%var_1d
              ENDIF
              CALL surface_restore_elements(                                   &
                                      m_liq_v(2)%var_1d,                       &
                                      tmp_walltype_v_1d(2)%var_1d,             &
                                      surf_lsm_v(2)%start_index,               & 
                                      start_index_on_file,                     &
                                      end_index_on_file,                       &
                                      nxlc, nysc,                              &
                                      nxlf, nxrf, nysf, nynf,                  &
                                      nys_on_file, nyn_on_file,                &
                                      nxl_on_file,nxr_on_file )

           CASE ( 'm_liq_v(3)' )
            
              IF ( k == 1 )  THEN
                 IF ( .NOT.  ALLOCATED( m_liq_v(3)%var_1d ) )                  &
                    ALLOCATE( m_liq_v(3)%var_1d(1:surf_lsm_v(3)%ns) )
                 READ ( 13 )  tmp_walltype_v_1d(3)%var_1d
              ENDIF
              CALL surface_restore_elements(                                   &
                                      m_liq_v(3)%var_1d,                       &
                                      tmp_walltype_v_1d(3)%var_1d,             &
                                      surf_lsm_v(3)%start_index,               & 
                                      start_index_on_file,                     &
                                      end_index_on_file,                       &
                                      nxlc, nysc,                              &
                                      nxlf, nxrf, nysf, nynf,                  &
                                      nys_on_file, nyn_on_file,                &
                                      nxl_on_file,nxr_on_file )


           CASE ( 't_surface_h' )
            
              IF ( k == 1 )  THEN
                 IF ( .NOT.  ALLOCATED( t_surface_h%var_1d ) )                 &
                    ALLOCATE( t_surface_h%var_1d(1:surf_lsm_h%ns) )
                 READ ( 13 )  tmp_walltype_h_1d%var_1d
              ENDIF
              CALL surface_restore_elements(                                   &
                                         t_surface_h%var_1d,                   &
                                         tmp_walltype_h_1d%var_1d,             &
                                         surf_lsm_h%start_index,               & 
                                         start_index_on_file,                  &
                                         end_index_on_file,                    &
                                         nxlc, nysc,                           &
                                         nxlf, nxrf, nysf, nynf,               &
                                         nys_on_file, nyn_on_file,             &
                                         nxl_on_file,nxr_on_file )

           CASE ( 't_surface_v(0)' )
            
              IF ( k == 1 )  THEN
                 IF ( .NOT.  ALLOCATED( t_surface_v(0)%var_1d ) )              &
                    ALLOCATE( t_surface_v(0)%var_1d(1:surf_lsm_v(0)%ns) )
                 READ ( 13 )  tmp_walltype_v_1d(0)%var_1d
              ENDIF
              CALL surface_restore_elements(                                   &
                                      t_surface_v(0)%var_1d,                   &
                                      tmp_walltype_v_1d(0)%var_1d,             &
                                      surf_lsm_v(0)%start_index,               & 
                                      start_index_on_file,                     &
                                      end_index_on_file,                       &
                                      nxlc, nysc,                              &
                                      nxlf, nxrf, nysf, nynf,                  &
                                      nys_on_file, nyn_on_file,                &
                                      nxl_on_file,nxr_on_file )

           CASE ( 't_surface_v(1)' )
            
              IF ( k == 1 )  THEN
                 IF ( .NOT.  ALLOCATED( t_surface_v(1)%var_1d ) )              &
                    ALLOCATE( t_surface_v(1)%var_1d(1:surf_lsm_v(1)%ns) )
                 READ ( 13 )  tmp_walltype_v_1d(1)%var_1d
              ENDIF
              CALL surface_restore_elements(                                   &
                                      t_surface_v(1)%var_1d,                   &
                                      tmp_walltype_v_1d(1)%var_1d,             &
                                      surf_lsm_v(1)%start_index,               & 
                                      start_index_on_file,                     &
                                      end_index_on_file,                       &
                                      nxlc, nysc,                              &
                                      nxlf, nxrf, nysf, nynf,                  &
                                      nys_on_file, nyn_on_file,                & 
                                      nxl_on_file,nxr_on_file )

           CASE ( 't_surface_v(2)' )
            
              IF ( k == 1 )  THEN
                 IF ( .NOT.  ALLOCATED( t_surface_v(2)%var_1d ) )              &
                    ALLOCATE( t_surface_v(2)%var_1d(1:surf_lsm_v(2)%ns) )
                 READ ( 13 )  tmp_walltype_v_1d(2)%var_1d
              ENDIF
              CALL surface_restore_elements(                                   &
                                      t_surface_v(2)%var_1d,                   &
                                      tmp_walltype_v_1d(2)%var_1d,             &
                                      surf_lsm_v(2)%start_index,               & 
                                      start_index_on_file,                     &
                                      end_index_on_file,                       &
                                      nxlc, nysc,                              &
                                      nxlf, nxrf, nysf, nynf,                  &
                                      nys_on_file, nyn_on_file,                &
                                      nxl_on_file,nxr_on_file )

           CASE ( 't_surface_v(3)' )
            
              IF ( k == 1 )  THEN
                 IF ( .NOT.  ALLOCATED( t_surface_v(3)%var_1d ) )              &
                    ALLOCATE( t_surface_v(3)%var_1d(1:surf_lsm_v(3)%ns) )
                 READ ( 13 )  tmp_walltype_v_1d(3)%var_1d
              ENDIF
              CALL surface_restore_elements(                                   &
                                      t_surface_v(3)%var_1d,                   &
                                      tmp_walltype_v_1d(3)%var_1d,             &
                                      surf_lsm_v(3)%start_index,               &   
                                      start_index_on_file,                     &
                                      end_index_on_file,                       &
                                      nxlc, nysc,                              &
                                      nxlf, nxrf, nysf, nynf,                  &
                                      nys_on_file, nyn_on_file,                &
                                      nxl_on_file,nxr_on_file )

          CASE DEFAULT

                found = .FALSE.

       END SELECT


 END SUBROUTINE lsm_rrd_local

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of roughness length for open water (lakes, ocean). The
!> parameterization follows Charnock (1955). Two different implementations
!> are available: as in ECMWF-IFS (Beljaars 1994) or as in FLake (Subin et al.
!> 2012)
!------------------------------------------------------------------------------!
    SUBROUTINE calc_z0_water_surface

       USE control_parameters,                                                 &
           ONLY: g, kappa, molecular_viscosity

       IMPLICIT NONE

       INTEGER(iwp) ::  i       !< running index
       INTEGER(iwp) ::  j       !< running index
       INTEGER(iwp) ::  m       !< running index

       REAL(wp), PARAMETER :: alpha_ch  = 0.018_wp !< Charnock constant (0.01-0.11). Use 0.01 for FLake and 0.018 for ECMWF
!       REAL(wp), PARAMETER :: pr_number = 0.71_wp !< molecular Prandtl number in the Charnock parameterization (differs from prandtl_number)
!       REAL(wp), PARAMETER :: sc_number = 0.66_wp !< molecular Schmidt number in the Charnock parameterization
!       REAL(wp) :: re_0 !< near-surface roughness Reynolds number

       DO  m = 1, surf_lsm_h%ns

          i   = surf_lsm_h%i(m)            
          j   = surf_lsm_h%j(m)
          
          IF ( surf_lsm_h%water_surface(m) )  THEN

!
!--          Disabled: FLake parameterization. Ideally, the Charnock 
!--          coefficient should depend on the water depth and the fetch
!--          length
!             re_0 = z0(j,i) * us(j,i) / molecular_viscosity
!        
!             z0(j,i) = MAX( 0.1_wp * molecular_viscosity / us(j,i),            &
!                           alpha_ch * us(j,i) / g )
!
!             z0h(j,i) = z0(j,i) * EXP( - kappa / pr_number * ( 4.0_wp * SQRT( re_0 ) - 3.2_wp ) )
!             z0q(j,i) = z0(j,i) * EXP( - kappa / pr_number * ( 4.0_wp * SQRT( re_0 ) - 4.2_wp ) )

!
!--           Set minimum roughness length for u* > 0.2
!             IF ( us(j,i) > 0.2_wp )  THEN
!                z0h(j,i) = MAX( 1.0E-5_wp, z0h(j,i) )
!                z0q(j,i) = MAX( 1.0E-5_wp, z0q(j,i) )
!             ENDIF

!
!--          ECMWF IFS model parameterization after Beljaars (1994). At low
!--          wind speed, the sea surface becomes aerodynamically smooth and
!--          the roughness scales with the viscosity. At high wind speed, the
!--          Charnock relation is used.
             surf_lsm_h%z0(m)  = ( 0.11_wp * molecular_viscosity /             &
                                 surf_lsm_h%us(m) )                            &
                               + ( alpha_ch * surf_lsm_h%us(m)**2 / g )

             surf_lsm_h%z0h(m) = 0.40_wp * molecular_viscosity /               &
                                 surf_lsm_h%us(m)
             surf_lsm_h%z0q(m) = 0.62_wp * molecular_viscosity /               &
                                 surf_lsm_h%us(m)

          ENDIF
       ENDDO

    END SUBROUTINE calc_z0_water_surface



 END MODULE land_surface_model_mod
