!> @file init_3d_model.f90
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
! 2018-11-15 cbegeman
! Modify call to init_pt_anomaly
!
! Former revisions:
! -----------------
! $Id: init_3d_model.f90 3083 2018-06-19 14:03:12Z gronemeier $
! Move initialization call for nudging and 1D/3D offline nesting.
! Revise initialization with inifor data.
!
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
!
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
!
! 3042 2018-05-25 10:44:37Z schwenkel
! Changed the name specific humidity to mixing ratio
!
! 3040 2018-05-25 10:22:08Z schwenkel
! Add option to initialize warm air bubble close to surface
!
! 3014 2018-05-09 08:42:38Z maronga
! Bugfix: initialization of ts_value missing
!
! 3011 2018-05-07 14:38:42Z schwenkel
! removed redundant if statement
!
! 3004 2018-04-27 12:33:25Z Giersch
! precipitation_rate removed
!
! 2995 2018-04-19 12:13:16Z Giersch
! CALL radiation_control is not necessary during initialization because
! calculation of radiative fluxes at model start is done in radiation_init
! in any case
!
! 2977 2018-04-17 10:27:57Z kanani
! Implement changes from branch radiation (r2948-2971) with minor modifications
! (moh.hefny):
! - set radiation_interactions according to the existence of urban/land vertical
!   surfaces and trees to activiate RTM
! - set average_radiation to TRUE if RTM is activiated
!
! 2938 2018-03-27 15:52:42Z suehring
! - Revise Inifor initialization for geostrophic wind components
! - Initialize synthetic turbulence generator in case of Inifor initialization
!
! 2936 2018-03-27 14:49:27Z suehring
! Synchronize parent and child models after initialization.
! Remove obsolete masking of topography grid points for Runge-Kutta weighted
! tendency arrays.
!
! 2920 2018-03-22 11:22:01Z kanani
! Add call for precalculating apparent solar positions (moh.hefny)
!
! 2906 2018-03-19 08:56:40Z Giersch
! The variables read/write_svf_on_init have been removed. Instead ENVIRONMENT
! variables read/write_svf have been introduced. Location_message has been
! added.
!
! 2894 2018-03-15 09:17:58Z Giersch
! Renamed routines with respect to reading restart data, file 13 is closed in
! rrd_read_parts_of_global now
!
! 2867 2018-03-09 09:40:23Z suehring
! Further bugfix concerning call of user_init.
!
! 2864 2018-03-08 11:57:45Z suehring
! Bugfix, move call of user_init in front of initialization of grid-point
! arrays
!
! 2817 2018-02-19 16:32:21Z knoop
! Preliminary gust module interface implemented
!
! 2776 2018-01-31 10:44:42Z Giersch
! Variable use_synthetic_turbulence_generator has been abbreviated
!
! 2766 2018-01-22 17:17:47Z kanani
! Removed preprocessor directive __chem
!
! 2758 2018-01-17 12:55:21Z suehring
! In case of spinup of land- and urban-surface model, do not mask wind velocity
! at first computational grid level
!
! 2746 2018-01-15 12:06:04Z suehring
! Move flag plant canopy to modules
!
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
!
! 2705 2017-12-18 11:26:23Z maronga
! Bugfix for reading initial profiles from ls/nuding file
!
! 2701 2017-12-15 15:40:50Z suehring
! Changes from last commit documented
!
! 2700 2017-12-15 14:12:35Z suehring
! Bugfix, missing initialization of surface attributes in case of
! inifor-initialization branch
!
! 2698 2017-12-14 18:46:24Z suehring
! Bugfix in get_topography_top_index
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Implementation of uv exposure model (FK)
! Moved initialisation of diss, e, kh, km to turbulence_closure_mod (TG)
! Added chemical emissions (FK)
! Initialize masking arrays and number-of-grid-points arrays before initialize
! LSM, USM and radiation module
! Initialization with inifor (MS)
!
! 2618 2017-11-16 15:37:30Z suehring
! Reorder calls of init_surfaces.
!
! 2564 2017-10-19 15:56:56Z Giersch
! Variable wind_turbine was added to control_parameters.
!
! 2550 2017-10-16 17:12:01Z boeske
! Modifications to cyclic fill method and turbulence recycling method in case of
! complex terrain simulations
!
! 2513 2017-10-04 09:24:39Z kanani
! Bugfix in storing initial scalar profile (wrong index)
!
! 2350 2017-08-15 11:48:26Z kanani
! Bugfix in nopointer version
!
! 2339 2017-08-07 13:55:26Z gronemeier
! corrected timestamp in header
!
! 2338 2017-08-07 12:15:38Z gronemeier
! Modularize 1D model
!
! 2329 2017-08-03 14:24:56Z knoop
! Removed temporary bugfix (r2327) as bug is properly resolved by this revision
!
! 2327 2017-08-02 07:40:57Z maronga
! Temporary bugfix
!
! 2320 2017-07-21 12:47:43Z suehring
! Modularize large-scale forcing and nudging
!
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison'
! includes two more prognostic equations for cloud drop concentration (nc)
! and cloud water content (qc).
!
! 2277 2017-06-12 10:47:51Z kanani
! Removed unused variable sums_up_fraction_l
!
! 2270 2017-06-09 12:18:47Z maronga
! dots_num must be increased when LSM and/or radiation is used
!
! 2259 2017-06-08 09:09:11Z gronemeier
! Implemented synthetic turbulence generator
!
! 2252 2017-06-07 09:35:37Z knoop
! rho_air now depending on surface_pressure even in Boussinesq mode
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography and surface concept:
!   - Modify passed parameters for disturb_field
!   - Topography representation via flags
!   - Remove unused arrays.
!   - Move initialization of surface-related quantities to surface_mod
!
! 2172 2017-03-08 15:55:25Z knoop
! Bugfix: moved parallel random generator initialization into its module
!
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC directives removed
!
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean
!
! 2011 2016-09-19 17:29:57Z kanani
! Flag urban_surface is now defined in module control_parameters.
!
! 2007 2016-08-24 15:47:17Z kanani
! Added support for urban surface model,
! adjusted location_message in case of plant_canopy
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1992 2016-08-12 15:14:59Z suehring
! Initializaton of scalarflux at model top
! Bugfixes in initialization of surface and top salinity flux, top scalar and
! humidity fluxes
!
! 1960 2016-07-12 16:34:24Z suehring
! Separate humidity and passive scalar
! Increase dimension for mean_inflow_profiles
! Remove inadvertent write-statement
! Bugfix, large-scale forcing is still not implemented for passive scalars
!
! 1957 2016-07-07 10:43:48Z suehring
! flight module added
!
! 1920 2016-05-30 10:50:15Z suehring
! Initialize us with very small number to avoid segmentation fault during
! calculation of Obukhov length
!
! 1918 2016-05-27 14:35:57Z raasch
! intermediate_timestep_count is set 0 instead 1 for first call of pres,
! bugfix: initialization of local sum arrays are moved to the beginning of the
!         routine because otherwise results from pres are overwritten
!
! 1914 2016-05-26 14:44:07Z witha
! Added initialization of the wind turbine model
!
! 1878 2016-04-19 12:30:36Z hellstea
! The zeroth element of weight_pres removed as unnecessary
!
! 1849 2016-04-08 11:33:18Z hoffmann
! Adapted for modularization of microphysics.
! precipitation_amount, precipitation_rate, prr moved to arrays_3d.
! Initialization of nc_1d, nr_1d, pt_1d, qc_1d, qr_1d, q_1d moved to
! microphysics_init.
!
! 1845 2016-04-08 08:29:13Z raasch
! nzb_2d replaced by nzb_u|v_inner
!
! 1833 2016-04-07 14:23:03Z raasch
! initialization of spectra quantities moved to spectra_mod
!
! 1831 2016-04-07 13:15:51Z hoffmann
! turbulence renamed collision_turbulence
!
! 1826 2016-04-07 12:01:39Z maronga
! Renamed radiation calls.
! Renamed canopy model calls.
!
! 1822 2016-04-07 07:49:42Z hoffmann
! icloud_scheme replaced by microphysics_*
!
! 1817 2016-04-06 15:44:20Z maronga
! Renamed lsm calls.
!
! 1815 2016-04-06 13:49:59Z raasch
! zero-settings for velocities inside topography re-activated (was deactivated
! in r1762)
!
! 1788 2016-03-10 11:01:04Z maronga
! Added z0q.
! Syntax layout improved.
!
! 1783 2016-03-06 18:36:17Z raasch
! netcdf module name changed + related changes
!
! 1764 2016-02-28 12:45:19Z raasch
! bugfix: increase size of volume_flow_area_l and volume_flow_initial_l by 1
!
! 1762 2016-02-25 12:31:13Z hellstea
! Introduction of nested domain feature
!
! 1738 2015-12-18 13:56:05Z raasch
! calculate mean surface level height for each statistic region
!
! 1734 2015-12-02 12:17:12Z raasch
! no initial disturbances in case that the disturbance energy limit has been
! set zero
!
! 1707 2015-11-02 15:24:52Z maronga
! Bugfix: transfer of Richardson number from 1D model to Obukhov length caused
! devision by zero in neutral stratification
!
! 1691 2015-10-26 16:17:44Z maronga
! Call to init_surface_layer added. rif is replaced by ol and zeta.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1615 2015-07-08 18:49:19Z suehring
! Enable turbulent inflow for passive_scalar and humidity
!
! 1585 2015-04-30 07:05:52Z maronga
! Initialization of radiation code is now done after LSM initializtion
!
! 1575 2015-03-27 09:56:27Z raasch
! adjustments for psolver-queries
!
! 1551 2015-03-03 14:18:16Z maronga
! Allocation of land surface arrays is now done in the subroutine lsm_init_arrays,
! which is part of land_surface_model.
!
! 1507 2014-12-10 12:14:18Z suehring
! Bugfix: set horizontal velocity components to zero inside topography
!
! 1496 2014-12-02 17:25:50Z maronga
! Added initialization of the land surface and radiation schemes
!
! 1484 2014-10-21 10:53:05Z kanani
! Changes due to new module structure of the plant canopy model:
! canopy-related initialization (e.g. lad and canopy_heat_flux) moved to new
! subroutine init_plant_canopy within the module plant_canopy_model_mod,
! call of subroutine init_plant_canopy added.
!
! 1431 2014-07-15 14:47:17Z suehring
! var_d added, in order to normalize spectra.
!
! 1429 2014-07-15 12:53:45Z knoop
! Ensemble run capability added to parallel random number generator
!
! 1411 2014-05-16 18:01:51Z suehring
! Initial horizontal velocity profiles were not set to zero at the first vertical
! grid level in case of non-cyclic lateral boundary conditions.
!
! 1406 2014-05-16 13:47:01Z raasch
! bugfix: setting of initial velocities at k=1 to zero not in case of a
! no-slip boundary condition for uv
!
! 1402 2014-05-09 14:25:13Z raasch
! location messages modified
!
! 1400 2014-05-09 14:03:54Z knoop
! Parallel random number generator added
!
! 1384 2014-05-02 14:31:06Z raasch
! location messages added
!
! 1361 2014-04-16 15:17:48Z hoffmann
! tend_* removed
! Bugfix: w_subs is not allocated anymore if it is already allocated
!
! 1359 2014-04-11 17:15:14Z hoffmann
! module lpm_init_mod added to use statements, because lpm_init has become a
! module
!
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1340 2014-03-25 19:45:13Z kanani
! REAL constants defined as wp-kind
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL constants defined as wp-kind
! module interfaces removed
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1316 2014-03-17 07:44:59Z heinze
! Bugfix: allocation of w_subs
!
! 1299 2014-03-06 13:15:21Z heinze
! Allocate w_subs due to extension of large scale subsidence in combination
! with large scale forcing data (LSF_DATA)
!
! 1241 2013-10-30 11:36:58Z heinze
! Overwrite initial profiles in case of nudging
! Inititialize shf and qsws in case of large_scale_forcing
!
! 1221 2013-09-10 08:59:13Z raasch
! +rflags_s_inner in copyin statement, use copyin for most arrays instead of
! copy
!
! 1212 2013-08-15 08:46:27Z raasch
! array tri is allocated and included in data copy statement
!
! 1195 2013-07-01 12:27:57Z heinze
! Bugfix: move allocation of ref_state to parin.f90 and read_var_list.f90
!
! 1179 2013-06-14 05:57:58Z raasch
! allocate and set ref_state to be used in buoyancy terms
!
! 1171 2013-05-30 11:27:45Z raasch
! diss array is allocated with full size if accelerator boards are used
!
! 1159 2013-05-21 11:58:22Z fricke
! -bc_lr_dirneu, bc_lr_neudir, bc_ns_dirneu, bc_ns_neudir
!
! 1153 2013-05-10 14:33:08Z raasch
! diss array is allocated with dummy elements even if it is not needed
! (required by PGI 13.4 / CUDA 5.0)
!
! 1115 2013-03-26 18:16:16Z hoffmann
! unused variables removed
!
! 1113 2013-03-10 02:48:14Z raasch
! openACC directive modified
!
! 1111 2013-03-08 23:54:10Z raasch
! openACC directives added for pres
! array diss allocated only if required
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1065 2012-11-22 17:42:36Z hoffmann
! allocation of diss (dissipation rate) in case of turbulence = .TRUE. added
!
! 1053 2012-11-13 17:11:03Z hoffmann
! allocation and initialisation of necessary data arrays for the two-moment
! cloud physics scheme the two new prognostic equations (nr, qr):
! +dr, lambda_r, mu_r, sed_*, xr, *s, *sws, *swst, *, *_p, t*_m, *_1, *_2, *_3,
! +tend_*, prr
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1032 2012-10-21 13:03:21Z letzel
! save memory by not allocating pt_2 in case of neutral = .T.
!
! 1025 2012-10-07 16:04:41Z letzel
! bugfix: swap indices of mask for ghost boundaries
!
! 1015 2012-09-27 09:23:24Z raasch
! mask is set to zero for ghost boundaries
!
! 1010 2012-09-20 07:59:54Z raasch
! cpp switch __nopointer added for pointer free version
!
! 1003 2012-09-14 14:35:53Z raasch
! nxra,nyna, nzta replaced ny nxr, nyn, nzt
!
! 1001 2012-09-13 14:08:46Z raasch
! all actions concerning leapfrog scheme removed
!
! 996 2012-09-07 10:41:47Z raasch
! little reformatting
!
! 978 2012-08-09 08:28:32Z fricke
! outflow damping layer removed
! roughness length for scalar quantites z0h added
! damping zone for the potential temperatur in case of non-cyclic lateral
! boundaries added
! initialization of ptdf_x, ptdf_y
! initialization of c_u_m, c_u_m_l, c_v_m, c_v_m_l, c_w_m, c_w_m_l
!
! 849 2012-03-15 10:35:09Z raasch
! init_particles renamed lpm_init
!
! 825 2012-02-19 03:03:44Z raasch
! wang_collision_kernel renamed wang_kernel
!
! Revision 1.1  1998/03/09 16:22:22  raasch
! Initial revision
!
!
! Description:
! ------------
!> Allocation of arrays and initialization of the 3D model via
!> a) pre-run the 1D model
!> or
!> b) pre-set constant linear profiles
!> or
!> c) read values of a previous run
!------------------------------------------------------------------------------!
 MODULE configure_3D_MODEL

    USE advec_ws

    USE arrays_3d

    USE constants,                                                             &
        ONLY:  pi

    USE control_parameters

    USE grid_variables,                                                        &
        ONLY:  dx, dy, ddx2_mg, ddy2_mg

    USE indices

    USE kinds

    USE netcdf_interface,                                                      &
        ONLY:  dots_max, dots_num, dots_unit, dots_label

    USE netcdf_data_input_mod,                                                 &
        ONLY:  init_3d, netcdf_data_input_interpolate, netcdf_data_input_init_3d

    USE pegrid

    USE random_function_mod

    USE random_generator_parallel,                                             &
        ONLY:  init_parallel_random_generator

    USE read_restart_data_mod,                                                 &
        ONLY:  rrd_read_parts_of_global, rrd_local

    USE statistics,                                                            &
        ONLY:  hom, hom_sum, mean_surface_level_height, pr_palm, rmask,        &
               statistic_regions, sums, sums_divnew_l, sums_divold_l, sums_l,  &
               sums_l_l, sums_wsts_bc_l, ts_value,                             &
               weight_pres, weight_substep

    USE surface_layer_fluxes_mod,                                              &
        ONLY:  init_surface_layer_fluxes

    USE surface_mod,                                                           &
        ONLY :  init_surface_arrays, init_surfaces, surf_def_h,     &
                get_topography_top_index_ji, vertical_surfaces_exist

    USE transpose_indices

    USE turbulence_closure_mod,                                                &
        ONLY:  tcm_init_arrays, tcm_init

    IMPLICIT NONE

    INTEGER(iwp) ::  i             !<
    INTEGER(iwp) ::  ind_array(1)  !<
    INTEGER(iwp) ::  j             !<
    INTEGER(iwp) ::  k             !<
    INTEGER(iwp) ::  k_surf        !< surface level index
    INTEGER(iwp) ::  m             !< index of surface element in surface data type
    INTEGER(iwp) ::  sr            !< index of statistic region

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE   ::  ngp_2dh_l  !<

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  ngp_2dh_outer_l    !<
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  ngp_2dh_s_inner_l  !<

    REAL(wp)     ::  t_surface !< air temperature at the surface

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  p_hydrostatic !< hydrostatic pressure

    INTEGER(iwp) ::  l       !< loop variable
    INTEGER(iwp) ::  nzt_l   !< index of top PE boundary for multigrid level
    REAL(wp) ::  dx_l !< grid spacing along x on different multigrid level
    REAL(wp) ::  dy_l !< grid spacing along y on different multigrid level

    REAL(wp), DIMENSION(1:3) ::  volume_flow_area_l     !<
    REAL(wp), DIMENSION(1:3) ::  volume_flow_initial_l  !<

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  mean_surface_level_height_l    !<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ngp_3d_inner_l    !<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ngp_3d_inner_tmp  !<

    INTEGER(iwp) ::  nz_u_shift   !<
    INTEGER(iwp) ::  nz_v_shift   !<
    INTEGER(iwp) ::  nz_w_shift   !<
    INTEGER(iwp) ::  nz_s_shift   !<
    INTEGER(iwp) ::  nz_u_shift_l !<
    INTEGER(iwp) ::  nz_v_shift_l !<
    INTEGER(iwp) ::  nz_w_shift_l !<
    INTEGER(iwp) ::  nz_s_shift_l !<

    public init_3d_model, deallocate_3d_variables
    contains

   subroutine deallocate_3d_variables

    DEALLOCATE( mean_surface_level_height,         &
              ngp_2dh, ngp_3d, ngp_3d_inner,          &
              ngp_3d_inner_tmp, sums_divold_l, sums_divnew_l)
    DEALLOCATE( dp_smooth_factor, rdf, rdf_sc )
    DEALLOCATE( ngp_2dh_outer, ngp_2dh_s_inner,                 &
              rmask,sums, sums_wsts_bc_l,   &
              ts_value )
    DEALLOCATE( ptdf_x, ptdf_y, weight_pres, weight_substep )

    DEALLOCATE( d, p, tend, sums_l, sums_l_l )

#if defined( __nopointer )
    DEALLOCATE( pt, pt_p, u, u_p, v, v_p, w, w_p, tpt_m, tu_m, tv_m, tw_m)
#else
    DEALLOCATE( pt_1, pt_2, pt_3, u_1, u_2, u_3, v_1,v_2, v_3, w_1, w_2, w_3)
#endif

!
!-- Array for storing constant coeffficients of the tridiagonal solver
    IF ( psolver == 'poisfft' )  THEN
       DEALLOCATE( tri, tric)
    ENDIF

#if defined( __nopointer )
       DEALLOCATE( prho, rho_ocean, alpha_T, beta_S, solar3d,                  &
                 sa, sa_p, tsa_m )
#else
       DEALLOCATE( prho_1,rho_1,alpha_T_1, beta_S_1, solar3d_1, sa_1, sa_2, sa_3 )
#endif

!
!-- Allocation of anelastic and Boussinesq approximation specific arrays
    DEALLOCATE( p_hydrostatic )
    DEALLOCATE( rho_air )
    DEALLOCATE( rho_air_zw )
    DEALLOCATE( drho_air )
    DEALLOCATE( drho_air_zw )

!
!-- Allocation of flux conversion arrays
    DEALLOCATE( heatflux_input_conversion )
    DEALLOCATE( waterflux_input_conversion )
    DEALLOCATE( momentumflux_input_conversion )
    DEALLOCATE( heatflux_output_conversion )
    DEALLOCATE( waterflux_output_conversion )
    DEALLOCATE( momentumflux_output_conversion )

   end subroutine deallocate_3d_variables

   subroutine init_3d_model
    CALL location_message( 'allocating arrays', .FALSE. )

!
!-- Allocate arrays
    ALLOCATE( mean_surface_level_height(0:statistic_regions),                  &
              mean_surface_level_height_l(0:statistic_regions),                &
              ngp_2dh(0:statistic_regions), ngp_2dh_l(0:statistic_regions),    &
              ngp_3d(0:statistic_regions),                                     &
              ngp_3d_inner(0:statistic_regions),                               &
              ngp_3d_inner_l(0:statistic_regions),                             &
              ngp_3d_inner_tmp(0:statistic_regions),                           &
              sums_divnew_l(0:statistic_regions),                              &
              sums_divold_l(0:statistic_regions) )
    ALLOCATE( dp_smooth_factor(nzb:nzt), rdf(nzb+1:nzt), rdf_sc(nzb+1:nzt) )
    ALLOCATE( ngp_2dh_outer(nzb:nzt+1,0:statistic_regions),                    &
              ngp_2dh_outer_l(nzb:nzt+1,0:statistic_regions),                  &
              ngp_2dh_s_inner(nzb:nzt+1,0:statistic_regions),                  &
              ngp_2dh_s_inner_l(nzb:nzt+1,0:statistic_regions),                &
              rmask(nysg:nyng,nxlg:nxrg,0:statistic_regions),                  &
              sums(nzb:nzt+1,pr_palm+max_pr_user),                             &
              sums_l(nzb:nzt+1,pr_palm+max_pr_user,0:threads_per_task-1),      &
              sums_l_l(nzb:nzt+1,0:statistic_regions,0:threads_per_task-1),    &
              sums_wsts_bc_l(nzb:nzt+1,0:statistic_regions),                   &
              ts_value(dots_max,0:statistic_regions) )
    ALLOCATE( ptdf_x(nxlg:nxrg), ptdf_y(nysg:nyng) )

    ALLOCATE( d(nzb+1:nzt,nys:nyn,nxl:nxr),                                    &
              p(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                &
              tend(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

#if defined( __nopointer )
    ALLOCATE( pt(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                               &
              pt_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                             &
              u(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                &
              u_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              v(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                &
              v_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              w(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                &
              w_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              tpt_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                            &
              tu_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                             &
              tv_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                             &
              tw_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
#else
    ALLOCATE( pt_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                             &
              pt_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                             &
              u_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              u_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              u_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              v_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              v_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              v_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              w_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              w_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              w_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    IF (  .NOT.  neutral )  THEN
       ALLOCATE( pt_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ENDIF
#endif

!
!-- Array for storing constant coeffficients of the tridiagonal solver
    IF ( psolver == 'poisfft' )  THEN
       ALLOCATE( tri(nxl_z:nxr_z,nys_z:nyn_z,0:nz-1,2) )
       ALLOCATE( tric(nxl_z:nxr_z,nys_z:nyn_z,0:nz-1) )
    ENDIF

#if defined( __nopointer )
       ALLOCATE( prho(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                          &
                 rho_ocean(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                     &
                 alpha_T(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                       &
                 beta_S(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                        &
                 solar3d(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                       &
                 sa(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                            &
                 sa_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                          &
                 tsa_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
#else
       ALLOCATE( prho_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                        &
                 rho_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                         &
                 alpha_T_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                       &
                 beta_S_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                        &
                 solar3d_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                       &
                 sa_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                          &
                 sa_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                          &
                 sa_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       prho => prho_1
       rho_ocean  => rho_1  ! routines calc_mean_profile and diffusion_e require
                      ! density to be apointer
       alpha_T => alpha_T_1
       beta_S => beta_S_1
       solar3d => solar3d_1
#endif

!
!-- Allocation of anelastic and Boussinesq approximation specific arrays
    ALLOCATE( p_hydrostatic(nzb:nzt+1) )
    ALLOCATE( rho_air(nzb:nzt+1) )
    ALLOCATE( rho_air_zw(nzb:nzt+1) )
    ALLOCATE( drho_air(nzb:nzt+1) )
    ALLOCATE( drho_air_zw(nzb:nzt+1) )

   rho_air(:) = 1.0_wp
   rho_air_zw(:) = 1.0_wp
!
!-- compute the inverse density array in order to avoid expencive divisions
    drho_air    = 1.0_wp / rho_air
    drho_air_zw = 1.0_wp / rho_air_zw

!
!-- Allocation of flux conversion arrays
    ALLOCATE( heatflux_input_conversion(nzb:nzt+1) )
    ALLOCATE( waterflux_input_conversion(nzb:nzt+1) )
    ALLOCATE( momentumflux_input_conversion(nzb:nzt+1) )
    ALLOCATE( heatflux_output_conversion(nzb:nzt+1) )
    ALLOCATE( waterflux_output_conversion(nzb:nzt+1) )
    ALLOCATE( momentumflux_output_conversion(nzb:nzt+1) )

!
!-- calculate flux conversion factors according to approximation and in-/output mode
    DO  k = nzb, nzt+1

        IF ( TRIM( flux_input_mode ) == 'kinematic' )  THEN
            heatflux_input_conversion(k)      = rho_air_zw(k)
            waterflux_input_conversion(k)     = rho_air_zw(k)
            momentumflux_input_conversion(k)  = rho_air_zw(k)
        ENDIF

        IF ( TRIM( flux_output_mode ) == 'kinematic' )  THEN
            heatflux_output_conversion(k)     = drho_air_zw(k)
            waterflux_output_conversion(k)    = drho_air_zw(k)
            momentumflux_output_conversion(k) = drho_air_zw(k)
        ENDIF

        IF ( .NOT. humidity ) THEN
            waterflux_input_conversion(k)  = 1.0_wp
            waterflux_output_conversion(k) = 1.0_wp
        ENDIF

    ENDDO
!
!-- 1D-array for large scale subsidence velocity
    IF ( .NOT. ALLOCATED( w_subs ) )  THEN
       ALLOCATE ( w_subs(nzb:nzt+1) )
       w_subs = 0.0_wp
    ENDIF

#if ! defined( __nopointer )
!
!-- Initial assignment of the pointers
    IF ( .NOT. neutral )  THEN
       pt => pt_1;  pt_p => pt_2;  tpt_m => pt_3
    ELSE
       pt => pt_1;  pt_p => pt_1;  tpt_m => pt_3
    ENDIF
    u  => u_1;   u_p  => u_2;   tu_m  => u_3
    v  => v_1;   v_p  => v_2;   tv_m  => v_3
    w  => w_1;   w_p  => w_2;   tw_m  => w_3
    sa => sa_1;  sa_p => sa_2;  tsa_m => sa_3
#endif
!
!-- Initialize arrays for turbulence closure
    CALL tcm_init_arrays
!
!-- Initialize surface arrays
    CALL init_surface_arrays
!
!-- Allocate arrays containing the RK coefficient for calculation of
!-- perturbation pressure and turbulent fluxes. At this point values are
!-- set for pressure calculation during initialization (where no timestep
!-- is done). Further below the values needed within the timestep scheme
!-- will be set.
    ALLOCATE( weight_substep(1:intermediate_timestep_count_max),               &
              weight_pres(1:intermediate_timestep_count_max) )
    weight_substep = 1.0_wp
    weight_pres    = 1.0_wp
    intermediate_timestep_count = 0  ! needed when simulated_time = 0.0

    CALL location_message( 'finished', .TRUE. )

!
!-- Initialize time series
    ts_value = 0.0_wp

!
!-- Initialize local summation arrays for routine flow_statistics.
!-- This is necessary because they may not yet have been initialized when they
!-- are called from flow_statistics (or - depending on the chosen model run -
!-- are never initialized)
    sums_divnew_l      = 0.0_wp
    sums_divold_l      = 0.0_wp
    sums_l_l           = 0.0_wp
    sums_wsts_bc_l     = 0.0_wp



!
!-- Initialize model variables
    IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.            &
         TRIM( initializing_actions ) /= 'cyclic_fill' )  THEN

       IF ( INDEX(initializing_actions, 'set_constant_profiles') /= 0 )    &
       THEN

          CALL location_message( 'initializing with constant profiles', .FALSE. )
!
!--       Use constructed initial profiles (velocity constant with height,
!--       temperature profile with constant gradient)
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                pt(:,j,i) = pt_init
                u(:,j,i)  = u_init
                v(:,j,i)  = v_init
             ENDDO
          ENDDO
!
!--       Mask topography
          u = MERGE( u, 0.0_wp, BTEST( wall_flags_0, 1 ) )
          v = MERGE( v, 0.0_wp, BTEST( wall_flags_0, 2 ) )
!
!--       Set initial horizontal velocities at the lowest computational grid
!--       levels to zero in order to avoid too small time steps caused by the
!--       diffusion limit in the initial phase of a run (at k=1, dz/2 occurs
!--       in the limiting formula!).
!--       Please note, in case land- or urban-surface model is used and a
!--       spinup is applied, masking the lowest computational level is not
!--       possible as MOST as well as energy-balance parametrizations will not
!--       work with zero wind velocity.
          IF ( ibc_uv_b /= 1  .AND.  .NOT.  spinup )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt
                      u(k,j,i) = MERGE( u(k,j,i), 0.0_wp,                      &
                                        BTEST( wall_flags_0(k,j,i), 20 ) )
                      v(k,j,i) = MERGE( v(k,j,i), 0.0_wp,                      &
                                        BTEST( wall_flags_0(k,j,i), 21 ) )
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   sa(:,j,i) = sa_init
                ENDDO
             ENDDO
!
!--       Compute initial temperature field and other constants used in case
!--       of a sloping surface
          IF ( sloping_surface )  CALL init_slope
!
!--       Initialize surface variables, e.g. friction velocity, momentum
!--       fluxes, etc.
          CALL init_surfaces
       ENDIF

       CALL location_message( 'initializing statistics, boundary conditions, etc.', &
                              .FALSE. )

!
!--    Bottom boundary
       IF ( ibc_uv_b == 0 .OR. ibc_uv_b == 2  )  THEN
          u(nzb,:,:) = 0.0_wp
          v(nzb,:,:) = 0.0_wp
       ENDIF

!
!--    Apply channel flow boundary condition
       IF ( TRIM( bc_uv_t ) == 'dirichlet_0' )  THEN
          u(nzt+1,:,:) = 0.0_wp
          v(nzt+1,:,:) = 0.0_wp
       ENDIF

!
!--    Store initial profiles for output purposes etc.. Please note, in case of
!--    initialization of u, v, w, pt, and q via output data derived from larger
!--    scale models, data will not be horizontally homogeneous. Actually, a mean
!--    profile should be calculated before.
       hom(:,1,5,:) = SPREAD( u(:,nys,nxl), 2, statistic_regions+1 )
       hom(:,1,6,:) = SPREAD( v(:,nys,nxl), 2, statistic_regions+1 )
       IF ( ibc_uv_b == 0 .OR. ibc_uv_b == 2)  THEN
          hom(nzb,1,5,:) = 0.0_wp
          hom(nzb,1,6,:) = 0.0_wp
       ENDIF
       hom(:,1,7,:)  = SPREAD( pt(:,nys,nxl), 2, statistic_regions+1 )


!
!--    Store initial salinity profile
          hom(:,1,26,:)  = SPREAD( sa(:,nys,nxl), 2, statistic_regions+1 )

!
!--    Initialize the random number generators (from numerical recipes)
       CALL random_function_ini

       IF ( random_generator == 'random-parallel' )  THEN
          CALL init_parallel_random_generator(nx, ny, nys, nyn, nxl, nxr)
       ENDIF
!
!--    Set the reference state to be used in the buoyancy terms (for ocean runs
!--    the reference state will be set (overwritten) in init_ocean)
       IF ( use_single_reference_value )  THEN
             ref_state(:) = pt_reference
       ELSE
             ref_state(:) = pt_init(:)
       ENDIF

!
!--    For the moment, vertical velocity is zero
       w = 0.0_wp

!
!--    Initialize array sums (must be defined in first call of pres)
       sums = 0.0_wp

!--    If required, change the surface temperature at the start of the 3D run
       IF ( pt_surface_initial_change /= 0.0_wp )  THEN
          pt(nzb,:,:) = pt(nzb,:,:) + pt_surface_initial_change
       ENDIF

!
!--    Initialize old and new time levels.
       tpt_m = 0.0_wp; tu_m = 0.0_wp; tv_m = 0.0_wp; tw_m = 0.0_wp
       pt_p = pt; u_p = u; v_p = v; w_p = w
       tsa_m = 0.0_wp
       sa_p  = sa

       CALL location_message( 'finished', .TRUE. )

    ELSEIF ( TRIM( initializing_actions ) == 'read_restart_data'  .OR.         &
             TRIM( initializing_actions ) == 'cyclic_fill' )                   &
    THEN

       CALL location_message( 'initializing in case of restart / cyclic_fill', &
                              .FALSE. )
!
!--    Initialize surface elements and its attributes, e.g. heat- and
!--    momentumfluxes, roughness, scaling parameters. As number of surface
!--    elements might be different between runs, e.g. in case of cyclic fill,
!--    and not all surface elements are read, surface elements need to be
!--    initialized before.
       CALL init_surfaces
!
!--    Read processor specific binary data from restart file
       DO  i = 0, io_blocks-1
          IF ( i == io_group )  THEN
             CALL rrd_local
          ENDIF
#if defined( __parallel )
          CALL MPI_BARRIER( comm2d, ierr )
#endif
       ENDDO
!
!
!--    Calculate initial temperature field and other constants used in case
!--    of a sloping surface
       IF ( sloping_surface )  CALL init_slope

!
!--    Initialize new time levels (only done in order to set boundary values
!--    including ghost points)
       pt_p = pt; u_p = u; v_p = v; w_p = w
       sa_p = sa

!
!--    Allthough tendency arrays are set in prognostic_equations, they have
!--    have to be predefined here because they are used (but multiplied with 0)
!--    there before they are set.
       tpt_m = 0.0_wp; tu_m = 0.0_wp; tv_m = 0.0_wp; tw_m = 0.0_wp
       tsa_m = 0.0_wp
!
       CALL location_message( 'finished', .TRUE. )

    ELSE
!
!--    Actually this part of the programm should not be reached
       message_string = 'unknown initializing problem'
       CALL message( 'init_3d_model', 'PA0193', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Initialize TKE, Kh and Km
    CALL tcm_init
!
!-- Before initializing further modules, compute total sum of active mask
!-- grid points and the mean surface level height for each statistic region.
!-- ngp_2dh: number of grid points of a horizontal cross section through the
!--          total domain
!-- ngp_3d:  number of grid points of the total domain
    ngp_2dh_outer_l   = 0
    ngp_2dh_outer     = 0
    ngp_2dh_s_inner_l = 0
    ngp_2dh_s_inner   = 0
    ngp_2dh_l         = 0
    ngp_2dh           = 0
    ngp_3d_inner_l    = 0.0_wp
    ngp_3d_inner      = 0
    ngp_3d            = 0
    ngp_sums          = ( nz + 2 ) * ( pr_palm + max_pr_user )

    mean_surface_level_height   = 0.0_wp
    mean_surface_level_height_l = 0.0_wp
!
!-- Pre-set masks for regional statistics. Default is the total model domain.
!-- Ghost points are excluded because counting values at the ghost boundaries
!-- would bias the statistics
    rmask = 1.0_wp
    rmask(:,nxlg:nxl-1,:) = 0.0_wp;  rmask(:,nxr+1:nxrg,:) = 0.0_wp
    rmask(nysg:nys-1,:,:) = 0.0_wp;  rmask(nyn+1:nyng,:,:) = 0.0_wp
!
!
!-- To do: New concept for these non-topography grid points!
    DO  sr = 0, statistic_regions
       DO  i = nxl, nxr
          DO  j = nys, nyn
             IF ( rmask(j,i,sr) == 1.0_wp )  THEN
!
!--             All xy-grid points
                ngp_2dh_l(sr) = ngp_2dh_l(sr) + 1
!
!--             Determine mean surface-level height. In case of downward-
!--             facing walls are present, more than one surface level exist.
!--             In this case, use the lowest surface-level height.
                IF ( surf_def_h(0)%start_index(j,i) <=                         &
                     surf_def_h(0)%end_index(j,i) )  THEN
                   m = surf_def_h(0)%start_index(j,i)
                   k = surf_def_h(0)%k(m)
                   mean_surface_level_height_l(sr) =                           &
                                       mean_surface_level_height_l(sr) + zw(k-1)
                ENDIF
                k_surf = k - 1

                DO  k = nzb, nzt+1
!
!--                xy-grid points above topography
                   ngp_2dh_outer_l(k,sr) = ngp_2dh_outer_l(k,sr)     +         &
                                  MERGE( 1, 0, BTEST( wall_flags_0(k,j,i), 24 ) )

                   ngp_2dh_s_inner_l(k,sr) = ngp_2dh_s_inner_l(k,sr) +         &
                                  MERGE( 1, 0, BTEST( wall_flags_0(k,j,i), 22 ) )

                ENDDO
!
!--             All grid points of the total domain above topography
                ngp_3d_inner_l(sr) = ngp_3d_inner_l(sr) + ( nz - k_surf + 2 )



             ENDIF
          ENDDO
       ENDDO
    ENDDO
!
!-- Initialize arrays encompassing number of grid-points in inner and outer
!-- domains, statistic regions, etc. Mainly used for horizontal averaging
!-- of turbulence statistics. Please note, user_init must be called before
!-- doing this.
    sr = statistic_regions + 1
#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( ngp_2dh_l(0), ngp_2dh(0), sr, MPI_INTEGER, MPI_SUM,    &
                        comm2d, ierr )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( ngp_2dh_outer_l(0,0), ngp_2dh_outer(0,0), (nz+2)*sr,   &
                        MPI_INTEGER, MPI_SUM, comm2d, ierr )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( ngp_2dh_s_inner_l(0,0), ngp_2dh_s_inner(0,0),          &
                        (nz+2)*sr, MPI_INTEGER, MPI_SUM, comm2d, ierr )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( ngp_3d_inner_l(0), ngp_3d_inner_tmp(0), sr, MPI_REAL,  &
                        MPI_SUM, comm2d, ierr )
    ngp_3d_inner = INT( ngp_3d_inner_tmp, KIND = SELECTED_INT_KIND( 18 ) )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( mean_surface_level_height_l(0),                        &
                        mean_surface_level_height(0), sr, MPI_REAL,            &
                        MPI_SUM, comm2d, ierr )
    mean_surface_level_height = mean_surface_level_height / REAL( ngp_2dh )
#else
    ngp_2dh         = ngp_2dh_l
    ngp_2dh_outer   = ngp_2dh_outer_l
    ngp_2dh_s_inner = ngp_2dh_s_inner_l
    ngp_3d_inner    = INT( ngp_3d_inner_l, KIND = SELECTED_INT_KIND( 18 ) )
    mean_surface_level_height = mean_surface_level_height_l / REAL( ngp_2dh_l )
#endif

    ngp_3d = INT ( ngp_2dh, KIND = SELECTED_INT_KIND( 18 ) ) * &
             INT ( (nz + 2 ), KIND = SELECTED_INT_KIND( 18 ) )

!
!-- Set a lower limit of 1 in order to avoid zero divisions in flow_statistics,
!-- buoyancy, etc. A zero value will occur for cases where all grid points of
!-- the respective subdomain lie below the surface topography
    ngp_2dh_outer   = MAX( 1, ngp_2dh_outer(:,:)   )
    ngp_3d_inner    = MAX( INT(1, KIND = SELECTED_INT_KIND( 18 )),             &
                           ngp_3d_inner(:) )
    ngp_2dh_s_inner = MAX( 1, ngp_2dh_s_inner(:,:) )

!-- Initialize quantities for special advections schemes
    CALL init_advec

!
!-- Impose random perturbation on the horizontal velocity field and then
!-- remove the divergences from the velocity field at the initial stage
    IF ( create_disturbances  .AND.  disturbance_energy_limit /= 0.0_wp  .AND. &
         TRIM( initializing_actions ) /= 'read_restart_data'  .AND.            &
         TRIM( initializing_actions ) /= 'cyclic_fill' )  THEN

       CALL location_message( 'creating initial disturbances', .FALSE. )
       CALL disturb_field( 'u', tend, u )
       CALL disturb_field( 'v', tend, v )
       CALL location_message( 'finished', .TRUE. )

       CALL location_message( 'calling pressure solver', .FALSE. )
       n_sor = nsor_ini
       !$acc data copyin( u, v, w, rho_air, rho_air_zw, ddzw, wall_flags_0, ngp_2dh_outer )
       CALL pres
       !$acc end data
       n_sor = nsor
       CALL location_message( 'finished', .TRUE. )

    ENDIF

!--    Initialize quantities needed for the ocean model
       CALL init_ocean

!-- Initialize surface layer (done after LSM as roughness length are required
!-- for initialization
    IF ( constant_flux_layer )  THEN
       CALL location_message( 'initializing surface layer', .FALSE. )
       CALL init_surface_layer_fluxes
       CALL location_message( 'finished', .TRUE. )
    ENDIF
!
!-- Initialize the ws-scheme.
    IF ( ws_scheme_sca .OR. ws_scheme_mom )  CALL ws_init

!
!-- Setting weighting factors for calculation of perturbation pressure
!-- and turbulent quantities from the RK substeps
    IF ( TRIM(timestep_scheme) == 'runge-kutta-3' )  THEN      ! for RK3-method

       weight_substep(1) = 1._wp/6._wp
       weight_substep(2) = 3._wp/10._wp
       weight_substep(3) = 8._wp/15._wp

       weight_pres(1)    = 1._wp/3._wp
       weight_pres(2)    = 5._wp/12._wp
       weight_pres(3)    = 1._wp/4._wp

    ELSEIF ( TRIM(timestep_scheme) == 'runge-kutta-2' )  THEN  ! for RK2-method

       weight_substep(1) = 1._wp/2._wp
       weight_substep(2) = 1._wp/2._wp

       weight_pres(1)    = 1._wp/2._wp
       weight_pres(2)    = 1._wp/2._wp

    ELSE                                     ! for Euler-method

       weight_substep(1) = 1.0_wp
       weight_pres(1)    = 1.0_wp

    ENDIF

!
!-- Initialize Rayleigh damping factors
    rdf    = 0.0_wp
    rdf_sc = 0.0_wp
    IF ( rayleigh_damping_factor /= 0.0_wp )  THEN
       IF (  .NOT.  ocean )  THEN
          DO  k = nzb+1, nzt
             IF ( zu(k) >= rayleigh_damping_height )  THEN
                rdf(k) = rayleigh_damping_factor *                             &
                      ( SIN( pi * 0.5_wp * ( zu(k) - rayleigh_damping_height ) &
                             / ( zu(nzt) - rayleigh_damping_height ) )         &
                      )**2
             ENDIF
          ENDDO
       ELSE
          DO  k = nzt, nzb+1, -1
             IF ( zu(k) <= rayleigh_damping_height )  THEN
                rdf(k) = rayleigh_damping_factor *                             &
                      ( SIN( pi * 0.5_wp * ( rayleigh_damping_height - zu(k) ) &
                             / ( rayleigh_damping_height - zu(nzb+1) ) )       &
                      )**2
             ENDIF
          ENDDO
       ENDIF
    ENDIF
    IF ( scalar_rayleigh_damping )  rdf_sc = rdf

!
!-- Initialize the starting level and the vertical smoothing factor used for
!-- the external pressure gradient
    dp_smooth_factor = 1.0_wp
    IF ( dp_external )  THEN
!
!--    Set the starting level dp_level_ind_b only if it has not been set before
!--    (e.g. in init_grid).
       IF ( dp_level_ind_b == 0 )  THEN
          ind_array = MINLOC( ABS( dp_level_b - zu ) )
          dp_level_ind_b = ind_array(1) - 1 + nzb
                                        ! MINLOC uses lower array bound 1
       ENDIF
       IF ( dp_smooth )  THEN
          dp_smooth_factor(:dp_level_ind_b) = 0.0_wp
          DO  k = dp_level_ind_b+1, nzt
             dp_smooth_factor(k) = 0.5_wp * ( 1.0_wp + SIN( pi *               &
                        ( REAL( k - dp_level_ind_b, KIND=wp ) /                &
                          REAL( nzt - dp_level_ind_b, KIND=wp ) - 0.5_wp ) ) )
          ENDDO
       ENDIF
    ENDIF

!
!-- Initialize damping zone for the potential temperature in case of
!-- non-cyclic lateral boundaries. The damping zone has the maximum value
!-- at the inflow boundary and decreases to zero at pt_damping_width.
    ptdf_x = 0.0_wp
    ptdf_y = 0.0_wp
    IF ( bc_lr_dirrad )  THEN
       DO  i = nxl, nxr
          IF ( ( i * dx ) < pt_damping_width )  THEN
             ptdf_x(i) = pt_damping_factor * ( SIN( pi * 0.5_wp *              &
                            REAL( pt_damping_width - i * dx, KIND=wp ) / (     &
                            REAL( pt_damping_width, KIND=wp ) ) ) )**2
          ENDIF
       ENDDO
    ELSEIF ( bc_lr_raddir )  THEN
       DO  i = nxl, nxr
          IF ( ( i * dx ) > ( nx * dx - pt_damping_width ) )  THEN
             ptdf_x(i) = pt_damping_factor *                                   &
                         SIN( pi * 0.5_wp *                                    &
                                 ( ( i - nx ) * dx + pt_damping_width ) /      &
                                 REAL( pt_damping_width, KIND=wp ) )**2
          ENDIF
       ENDDO
    ELSEIF ( bc_ns_dirrad )  THEN
       DO  j = nys, nyn
          IF ( ( j * dy ) > ( ny * dy - pt_damping_width ) )  THEN
             ptdf_y(j) = pt_damping_factor *                                   &
                         SIN( pi * 0.5_wp *                                    &
                                 ( ( j - ny ) * dy + pt_damping_width ) /      &
                                 REAL( pt_damping_width, KIND=wp ) )**2
          ENDIF
       ENDDO
    ELSEIF ( bc_ns_raddir )  THEN
       DO  j = nys, nyn
          IF ( ( j * dy ) < pt_damping_width )  THEN
             ptdf_y(j) = pt_damping_factor *                                   &
                         SIN( pi * 0.5_wp *                                    &
                                ( pt_damping_width - j * dy ) /                &
                                REAL( pt_damping_width, KIND=wp ) )**2
          ENDIF
       ENDDO
    ENDIF
!
!-- Check if maximum number of allowed timeseries is exceeded
    IF ( dots_num > dots_max )  THEN
       WRITE( message_string, * ) 'number of time series quantities exceeds',  &
                                  ' its maximum of dots_max = ', dots_max,     &
                                  '&Please increase dots_max in modules.f90.'
       CALL message( 'init_3d_model', 'PA0194', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Input binary data file is not needed anymore. This line must be placed
!-- after call of user_init!
    CALL close_file( 13 )
!
!-- In case of nesting, put an barrier to assure that all parent and child
!-- domains finished initialization.

    Deallocate(mean_surface_level_height_l, ngp_2dh_l, ngp_3d_inner_l)
    DEALLOCATE(ngp_2dh_outer_l, ngp_2dh_s_inner_l)
    CALL location_message( 'leaving init_3d_model', .TRUE. )

 END SUBROUTINE init_3d_model

 END MODULE configure_3D_MODEL
