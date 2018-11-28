!> @file check_parameters.f90
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
! 2018-11-01 cbegeman
! Add checks for bubble initial conditions
!
! 2018-10-31 cbegeman
! Add checks for profile output
!
! 2018-10-25 cbegeman
! Add checks for dirichlet bottom boundary conditions for salinity
!
! Former revisions:
! -----------------
! $Id: check_parameters.f90 2520 2017-10-05 13:50:26Z gronemeier &
! Add inital profile output for e (TG)
!
! 3065 2018-06-12 07:03:02Z Giersch
! dz was replaced by dz(1), error message revised
!
! 3049 2018-05-29 13:52:36Z Giersch
! add variable description
!
! 3046 2018-05-29 08:02:15Z Giersch
! Error messages revised
!
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
!
! 3035 2018-05-24 09:35:20Z schwenkel
! Add option to initialize warm air bubble close to surface
!
! 3034 2018-05-24 08:41:20Z raasch
! bugfix: check that initializing_actions has been set
!
! 2980 2018-04-17 15:19:27Z suehring
! Further improvement for spinup checks.
!
! 2974 2018-04-16 12:59:52Z gronemeier
! Bugfix: check if dt_data_output_av is zero in case of parallel NetCDF output
!
! 2970 2018-04-13 15:09:23Z suehring
! Bugfix in old large-scale forcing mode
!
! 2964 2018-04-12 16:04:03Z Giersch
! Calculation of fixed number of output time levels for parallel netcdf output
! has been revised (based on calculations in netcdf_interface_mod)
!
! 2938 2018-03-27 15:52:42Z suehring
! - Revise start and end indices for imposing random disturbances in case of
!   nesting or large-scale forcing.
! - Remove check for inifor initialization in case of nested runs.
! - Adapt call to check for synthetic turbulence geneartor settings.
!
! 2936 2018-03-27 14:49:27Z suehring
! Check spinup in case of nested runs, spinup time and timestep must be
! identical to assure synchronuous simulations.
!
! 2932 2018-03-26 09:39:22Z maronga
! renamed particles_par to particle_parameters
!
! 2918 2018-03-21 15:52:14Z gronemeier
! Add check for 1D model
!
! 2883 2018-03-14 08:29:10Z Giersch
! dt_dopr_listing is not set to the default value zero anymore
!
! 2851 2018-03-05 14:39:31Z maronga
! Bugfix: calculation of output time levels in case of restart runs (parallel
! NetCDF)
!
! 2836 2018-02-26 13:40:05Z Giersch
! dt_dopr_listing is set to the default value zero
!
! 2817 2018-02-19 16:32:21Z knoop
! Preliminary gust module interface implemented
!
! 2798 2018-02-09 17:16:39Z suehring
! Consider also default-type surfaces for surface temperature output.
!
! 2797 2018-02-08 13:24:35Z suehring
! Enable output of ground-heat flux also at urban surfaces.
!
! 2776 2018-01-31 10:44:42Z Giersch
! Variable synthetic_turbulence_generator has been abbreviated
!
! 2773 2018-01-30 14:12:54Z suehring
! Check for consistent initialization in nesting mode added.
!
! 2766 2018-01-22 17:17:47Z kanani
! Removed preprocessor directive __chem
!
! 2765 2018-01-22 11:34:58Z maronga
! Renamed simulation_time_since_reference to
! time_to_be_simulated_from_reference_point
!
! 2746 2018-01-15 12:06:04Z suehring
! Move flag plant canopy to modules
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
! + new possible value for dissipation_1d
! Added checks for turbulence_closure_mod (TG)
! Implementation of chemistry module (FK)
!
! 2689 2017-12-12 17:46:55Z Giersch
! Bugfix in if query
!
! 2688 2017-12-12 17:27:04Z Giersch
! Check if humidity is set to TRUE in the _p3d file for coupled runs
!
! 2669 2017-12-06 16:03:27Z raasch
! mrun-string replaced by palmrun
!
! 2628 2017-11-20 12:40:38Z schwenkel
! Enabled particle advection with grid stretching -> Removed parameter check
!
! 2575 2017-10-24 09:57:58Z maronga
! Renamed phi --> latitude
!
! 2564 2017-10-19 15:56:56Z Giersch
! Variable wind_turbine was added to control_parameters.
!
! 2550 2017-10-16 17:12:01Z boeske
! Added checks for complex terrain simulations
!
! 2513 2017-10-04 09:24:39Z kanani
! Bugfix for some dopr(_initial)_index values and units connected to
! passive-scalar output
!
! 2508 2017-10-02 08:57:09Z suehring
! Bugfix, change default value of vertical_gradient level in order to consider
! also ocean runs
!
! 2422 2017-09-08 08:25:41Z raasch
! error message in case of missing "restart" file activation string
!
! 2375 2017-08-29 14:10:28Z schwenkel
! Added aerosol for bulk microphysics
!
! 2365 2017-08-21 14:59:59Z kanani
! Vertical grid nesting implemented: Check coupling mode. Generate file header
! (SadiqHuq)
!
! 2354 2017-08-17 10:49:36Z schwenkel
! Bugfix correlated to lsm_check_data_output_pr.
! If-statement for following checks is essential, otherwise units for lsm output
! are set to 'illegal' and palm will be aborted.
!
! 2348 2017-08-10 10:40:10Z kanani
! New: Check for simultaneous use of geostrophic wind and u_profile/v_profile
!
! 2345 2017-08-09 11:50:30Z Giersch
! Remove error message PA0156 and the conserve_volume_flow_mode option
! inflow_profile
!
! 2339 2017-08-07 13:55:26Z raasch
! corrected timestamp in header
!
! 2338 2017-08-07 12:15:38Z gronemeier
! Modularize 1D model
!
! 2329 2017-08-03 14:24:56Z knoop
! Bugfix: index corrected for rho_air and rho_air_zw output
!
! 2320 2017-07-21 12:47:43Z suehring
! Modularize large-scale forcing and nudging
!
! 2312 2017-07-14 20:26:51Z hoffmann
! PA0349 and PA0420 removed.
!
! 2300 2017-06-29 13:31:14Z raasch
! host-specific settings and checks removed
!
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison'
! includes two more prognostic equations for cloud drop concentration (nc)
! and cloud water content (qc).
!
! 2274 2017-06-09 13:27:48Z Giersch
! Changed error messages
!
! 2271 2017-06-09 12:34:55Z sward
! roughness-length check altered
! Error messages fixed
!
! 2270 2017-06-09 12:18:47Z maronga
! Revised numbering (removed 2 timeseries)
!
! 2259 2017-06-08 09:09:11Z gronemeier
! Implemented synthetic turbulence generator
!
! 2251 2017-06-06 15:10:46Z Giersch
!
! 2250 2017-06-06 15:07:12Z Giersch
! Doxygen comment added
!
! 2248 2017-06-06 13:52:54Z sward
! Error message fixed
!
! 2209 2017-04-19 09:34:46Z kanani
! Check for plant canopy model output
!
! 2200 2017-04-11 11:37:51Z suehring
! monotonic_adjustment removed
!
! 2178 2017-03-17 11:07:39Z hellstea
! Index limits for perturbations are now set also in case of nested
! boundary conditions
!
! 2169 2017-03-06 18:16:35Z suehring
! Bugfix, move setting for topography grid convention to init_grid
!
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC related parts of code removed
!
! 2088 2016-12-19 16:30:25Z suehring
! Bugfix in initial salinity profile
!
! 2084 2016-12-09 15:59:42Z knoop
! Removed anelastic + multigrid error message
!
! 2052 2016-11-08 15:14:59Z gronemeier
! Bugfix: remove setting of default value for recycling_width
!
! 2050 2016-11-08 15:00:55Z gronemeier
! Implement turbulent outflow condition
!
! 2044 2016-11-02 16:44:25Z knoop
! Added error code for anelastic approximation
!
! 2042 2016-11-02 13:47:31Z suehring
! Additional checks for wall_heatflux, wall_humidityflux and wall_scalarflux.
! Bugfix, check for constant_scalarflux.
!
! 2040 2016-10-26 16:58:09Z gronemeier
! Removed check for statistic_regions > 9.
!
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean
!
! 2026 2016-10-18 10:27:02Z suehring
! Bugfix, enable output of s*2
!
! 2024 2016-10-12 16:42:37Z kanani
! Added missing CASE, error message and unit for ssws*,
! increased number of possible output quantities to 500.
!
! 2011 2016-09-19 17:29:57Z kanani
! Flag urban_surface is now defined in module control_parameters,
! changed prefix for urban surface model output to "usm_",
! added flag lsf_exception (inipar-Namelist parameter) to allow explicit
! enabling of large scale forcing together with buildings on flat terrain,
! introduced control parameter varnamelength for LEN of var.
!
! 2007 2016-08-24 15:47:17Z kanani
! Added checks for the urban surface model,
! increased counter in DO WHILE loop over data_output (for urban surface output)
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1994 2016-08-15 09:52:21Z suehring
! Add missing check for cloud_physics and cloud_droplets
!
! 1992 2016-08-12 15:14:59Z suehring
! New checks for top_scalarflux
! Bugfixes concerning data output of profiles in case of passive_scalar
!
! 1984 2016-08-01 14:48:05Z suehring
! Bugfix: checking for bottom and top boundary condition for humidity and scalars
!
! 1972 2016-07-26 07:52:02Z maronga
! Removed check of lai* (done in land surface module)
!
! 1970 2016-07-18 14:27:12Z suehring
! Bugfix, check for ibc_q_b and constant_waterflux in case of land-surface scheme
!
! 1962 2016-07-13 09:16:44Z suehring
! typo removed
!
! 1960 2016-07-12 16:34:24Z suehring
! Separate checks and calculations for humidity and passive scalar. Therefore,
! a subroutine to calculate vertical gradients is introduced, in order  to reduce
! redundance.
! Additional check - large-scale forcing is not implemented for passive scalars
! so far.
!
! 1931 2016-06-10 12:06:59Z suehring
! Rename multigrid into multigrid_noopt and multigrid_fast into multigrid
!
! 1929 2016-06-09 16:25:25Z suehring
! Bugfix in check for use_upstream_for_tke
!
! 1918 2016-05-27 14:35:57Z raasch
! setting of a fixed reference state ('single_value') for ocean runs removed
!
! 1914 2016-05-26 14:44:07Z witha
! Added checks for the wind turbine model
!
! 1841 2016-04-07 19:14:06Z raasch
! redundant location message removed
!
! 1833 2016-04-07 14:23:03Z raasch
! check of spectra quantities moved to spectra_mod
!
! 1829 2016-04-07 12:16:29Z maronga
! Bugfix: output of user defined data required reset of the string 'unit'
!
! 1826 2016-04-07 12:01:39Z maronga
! Moved checks for radiation model to the respective module.
! Moved checks for plant canopy model to the respective module.
! Bugfix for check of too large roughness_length
!
! 1824 2016-04-07 09:55:29Z gronemeier
! Check if roughness_length < dz/2. Moved location_message(finished) to the end.
!
! 1822 2016-04-07 07:49:42Z hoffmann
! PALM collision kernel deleted. Collision algorithms are checked for correct
! spelling.
!
! Tails removed.
! !
! Checks for use_sgs_for_particles adopted for the use of droplets with
! use_sgs_for_particles adopted.
!
! Unused variables removed.
!
! 1817 2016-04-06 15:44:20Z maronga
! Moved checks for land_surface model to the respective module
!
! 1806 2016-04-05 18:55:35Z gronemeier
! Check for recycling_yshift
!
! 1804 2016-04-05 16:30:18Z maronga
! Removed code for parameter file check (__check)
!
! 1795 2016-03-18 15:00:53Z raasch
! Bugfix: setting of initial scalar profile in ocean
!
! 1788 2016-03-10 11:01:04Z maronga
! Added check for use of most_method = 'lookup' in combination with water
! surface presribed in the land surface model. Added output of z0q.
! Syntax layout improved.
!
! 1786 2016-03-08 05:49:27Z raasch
! cpp-direktives for spectra removed, check of spectra level removed
!
! 1783 2016-03-06 18:36:17Z raasch
! netcdf variables and module name changed,
! check of netcdf precision removed (is done in the netcdf module)
!
! 1764 2016-02-28 12:45:19Z raasch
! output of nest id in run description header,
! bugfix: check of validity of lateral boundary conditions moved to parin
!
! 1762 2016-02-25 12:31:13Z hellstea
! Introduction of nested domain feature
!
! 1745 2016-02-05 13:06:51Z gronemeier
! Bugfix: check data output intervals to be /= 0.0 in case of parallel NetCDF4
!
! 1701 2015-11-02 07:43:04Z maronga
! Bugfix: definition of rad_net timeseries was missing
!
! 1691 2015-10-26 16:17:44Z maronga
! Added output of Obukhov length (ol) and radiative heating rates for RRTMG.
! Added checks for use of radiation / lsm with topography.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1606 2015-06-29 10:43:37Z maronga
! Added check for use of RRTMG without netCDF.
!
! 1587 2015-05-04 14:19:01Z maronga
! Added check for set of albedos when using RRTMG
!
! 1585 2015-04-30 07:05:52Z maronga
! Added support for RRTMG
!
! 1575 2015-03-27 09:56:27Z raasch
! multigrid_fast added as allowed pressure solver
!
! 1573 2015-03-23 17:53:37Z suehring
! Bugfix: check for advection schemes in case of non-cyclic boundary conditions
!
! 1557 2015-03-05 16:43:04Z suehring
! Added checks for monotonic limiter
!
! 1555 2015-03-04 17:44:27Z maronga
! Added output of r_a and r_s. Renumbering of LSM PA-messages.
!
! 1553 2015-03-03 17:33:54Z maronga
! Removed check for missing soil temperature profile as default values were added.
!
! 1551 2015-03-03 14:18:16Z maronga
! Added various checks for land surface and radiation model. In the course of this
! action, the length of the variable var has to be increased
!
! 1504 2014-12-04 09:23:49Z maronga
! Bugfix: check for large_scale forcing got lost.
!
! 1500 2014-12-03 17:42:41Z maronga
! Boundary conditions changed to dirichlet for land surface model
!
! 1496 2014-12-02 17:25:50Z maronga
! Added checks for the land surface model
!
! 1484 2014-10-21 10:53:05Z kanani
! Changes due to new module structure of the plant canopy model:
!   module plant_canopy_model_mod added,
!   checks regarding usage of new method for leaf area density profile
!   construction added,
!   lad-profile construction moved to new subroutine init_plant_canopy within
!   the module plant_canopy_model_mod,
!   drag_coefficient renamed to canopy_drag_coeff.
! Missing KIND-attribute for REAL constant added
!
! 1455 2014-08-29 10:47:47Z heinze
! empty time records in volume, cross-section and masked data output prevented
! in case of non-parallel netcdf-output in restart runs
!
! 1429 2014-07-15 12:53:45Z knoop
! run_description_header exended to provide ensemble_member_nr if specified
!
! 1425 2014-07-05 10:57:53Z knoop
! bugfix: perturbation domain modified for parallel random number generator
!
! 1402 2014-05-09 14:25:13Z raasch
! location messages modified
!
! 1400 2014-05-09 14:03:54Z knoop
! Check random generator extended by option random-parallel
!
! 1384 2014-05-02 14:31:06Z raasch
! location messages added
!
! 1365 2014-04-22 15:03:56Z boeske
! Usage of large scale forcing for pt and q enabled:
! output for profiles of large scale advection (td_lsa_lpt, td_lsa_q),
! large scale subsidence (td_sub_lpt, td_sub_q)
! and nudging tendencies (td_nud_lpt, td_nud_q, td_nud_u and td_nud_v) added,
! check if use_subsidence_tendencies is used correctly
!
! 1361 2014-04-16 15:17:48Z hoffmann
! PA0363 removed
! PA0362 changed
!
! 1359 2014-04-11 17:15:14Z hoffmann
! Do not allow the execution of PALM with use_particle_tails, since particle
! tails are currently not supported by our new particle structure.
!
! PA0084 not necessary for new particle structure
!
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1330 2014-03-24 17:29:32Z suehring
! In case of SGS-particle velocity advection of TKE is also allowed with
! dissipative 5th-order scheme.
!
! 1327 2014-03-21 11:00:16Z raasch
! "baroclinicity" renamed "baroclinity", "ocean version" replaced by "ocean mode"
! bugfix: duplicate error message 56 removed,
! check of data_output_format and do3d_compress removed
!
! 1322 2014-03-20 16:38:49Z raasch
! some REAL constants defined as wp-kind
!
! 1320 2014-03-20 08:40:49Z raasch
! Kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1308 2014-03-13 14:58:42Z fricke
! +netcdf_data_format_save
! Calculate fixed number of output time levels for parallel netcdf output.
! For masked data, parallel netcdf output is not tested so far, hence
! netcdf_data_format is switched back to non-paralell output.
!
! 1299 2014-03-06 13:15:21Z heinze
! enable usage of large_scale subsidence in combination with large_scale_forcing
! output for profile of large scale vertical velocity w_subs added
!
! 1276 2014-01-15 13:40:41Z heinze
! Use LSF_DATA also in case of Dirichlet bottom boundary condition for scalars
!
! 1241 2013-10-30 11:36:58Z heinze
! output for profiles of ug and vg added
! set surface_heatflux and surface_waterflux also in dependence on
! large_scale_forcing
! checks for nudging and large scale forcing from external file
!
! 1236 2013-09-27 07:21:13Z raasch
! check number of spectra levels
!
! 1216 2013-08-26 09:31:42Z raasch
! check for transpose_compute_overlap (temporary)
!
! 1214 2013-08-21 12:29:17Z kanani
! additional check for simultaneous use of vertical grid stretching
! and particle advection
!
! 1212 2013-08-15 08:46:27Z raasch
! checks for poisfft_hybrid removed
!
! 1210 2013-08-14 10:58:20Z raasch
! check for fftw
!
! 1179 2013-06-14 05:57:58Z raasch
! checks and settings of buoyancy parameters and switches revised,
! initial profile for rho_ocean added to hom (id=77)
!
! 1174 2013-05-31 10:28:08Z gryschka
! Bugfix in computing initial profiles for ug, vg, lad, q in case of Atmosphere
!
! 1159 2013-05-21 11:58:22Z fricke
! bc_lr/ns_dirneu/neudir removed
!
! 1115 2013-03-26 18:16:16Z hoffmann
! unused variables removed
! drizzle can be used without precipitation
!
! 1111 2013-03-08 23:54:10Z raasch
! ibc_p_b = 2 removed
!
! 1103 2013-02-20 02:15:53Z raasch
! Bugfix: turbulent inflow must not require cyclic fill in restart runs
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1069 2012-11-28 16:18:43Z maronga
! allow usage of topography in combination with cloud physics
!
! 1065 2012-11-22 17:42:36Z hoffmann
! Bugfix: It is not allowed to use cloud_scheme = seifert_beheng without
!         precipitation in order to save computational resources.
!
! 1060 2012-11-21 07:19:51Z raasch
! additional check for parameter turbulent_inflow
!
! 1053 2012-11-13 17:11:03Z hoffmann
! necessary changes for the new two-moment cloud physics scheme added:
! - check cloud physics scheme (Kessler or Seifert and Beheng)
! - plant_canopy is not allowed
! - currently, only cache loop_optimization is allowed
! - initial profiles of nr, qr
! - boundary condition of nr, qr
! - check output quantities (qr, nr, prr)
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1031/1034 2012-10-22 11:32:49Z raasch
! check of netcdf4 parallel file support
!
! 1019 2012-09-28 06:46:45Z raasch
! non-optimized version of prognostic_equations not allowed any more
!
! 1015 2012-09-27 09:23:24Z raasch
! acc allowed for loop optimization,
! checks for adjustment of mixing length to the Prandtl mixing length removed
!
! 1003 2012-09-14 14:35:53Z raasch
! checks for cases with unequal subdomain sizes removed
!
! 1001 2012-09-13 14:08:46Z raasch
! all actions concerning leapfrog- and upstream-spline-scheme removed
!
! 996 2012-09-07 10:41:47Z raasch
! little reformatting

! 978 2012-08-09 08:28:32Z fricke
! setting of bc_lr/ns_dirneu/neudir
! outflow damping layer removed
! check for z0h*
! check for pt_damping_width
!
! 964 2012-07-26 09:14:24Z raasch
! check of old profil-parameters removed
!
! 940 2012-07-09 14:31:00Z raasch
! checks for parameter neutral
!
! 924 2012-06-06 07:44:41Z maronga
! Bugfix: preprocessor directives caused error during compilation
!
! 892 2012-05-02 13:51:44Z maronga
! Bugfix for parameter file check ( excluding __netcdf4 )
!
! 866 2012-03-28 06:44:41Z raasch
! use only 60% of the geostrophic wind as translation speed in case of Galilean
! transformation and use_ug_for_galilei_tr = .T. in order to mimimize the
! timestep
!
! 861 2012-03-26 14:18:34Z suehring
! Check for topography and ws-scheme removed.
! Check for loop_optimization = 'vector' and ws-scheme removed.
!
! 845 2012-03-07 10:23:05Z maronga
! Bugfix: exclude __netcdf4 directive part from namelist file check compilation
!
! 828 2012-02-21 12:00:36Z raasch
! check of collision_kernel extended
!
! 825 2012-02-19 03:03:44Z raasch
! check for collision_kernel and curvature_solution_effects
!
! 809 2012-01-30 13:32:58Z maronga
! Bugfix: replaced .AND. and .NOT. with && and ! in the preprocessor directives
!
! 807 2012-01-25 11:53:51Z maronga
! New cpp directive "__check" implemented which is used by check_namelist_files
!
! Revision 1.1  1997/08/26 06:29:23  raasch
! Initial revision
!
!
! Description:
! ------------
!> Check control parameters and deduce further quantities.
!------------------------------------------------------------------------------!
 SUBROUTINE check_parameters


    USE arrays_3d
    USE chemistry_model_mod,                                                   &
        ONLY:  chem_boundary_conds, chem_check_data_output,                    &
               chem_check_data_output_pr, chem_species
    USE chem_modules
    USE cloud_parameters
    USE constants
    USE control_parameters
    USE dvrp_variables
    USE grid_variables
    USE gust_mod,                                                              &
        ONLY: gust_check_data_output, gust_check_data_output_pr,               &
              gust_check_parameters, gust_module_enabled
    USE indices
    USE land_surface_model_mod,                                                &
        ONLY: lsm_check_data_output, lsm_check_data_output_pr,                 &
              lsm_check_parameters

    USE lsf_nudging_mod,                                                       &
        ONLY:  lsf_nudging_check_parameters, lsf_nudging_check_data_output_pr

    USE kinds

    USE model_1d_mod,                                                          &
        ONLY:  damp_level_1d, damp_level_ind_1d

    USE netcdf_data_input_mod,                                                 &
        ONLY:  init_model, input_pids_static, netcdf_data_input_check_dynamic, &
               netcdf_data_input_check_static

    USE netcdf_interface,                                                      &
        ONLY:  dopr_unit, do2d_unit, do3d_unit, netcdf_data_format,            &
               netcdf_data_format_string, dots_unit, heatflux_output_unit,     &
               waterflux_output_unit, momentumflux_output_unit

    USE particle_attributes
    USE pegrid
    USE plant_canopy_model_mod,                                                &
        ONLY:  pcm_check_data_output, pcm_check_parameters

    USE pmc_interface,                                                         &
        ONLY:  cpl_id, nested_run

    USE profil_parameter
    USE radiation_model_mod,                                                   &
        ONLY:  radiation, radiation_check_data_output,                         &
               radiation_check_data_output_pr, radiation_check_parameters
    USE spectra_mod,                                                           &
        ONLY:  calculate_spectra, spectra_check_parameters
    USE statistics
    USE stokes_drift_mod,                                                      &
        ONLY: stokes_drift_check_parameters
    USE subsidence_mod
    USE statistics
    USE synthetic_turbulence_generator_mod,                                    &
        ONLY:  stg_check_parameters
    USE transpose_indices
    USE turbulence_closure_mod,                                                &
        ONLY:  tcm_check_data_output, tcm_check_parameters
    USE urban_surface_mod,                                                     &
        ONLY:  usm_check_data_output, usm_check_parameters
    USE uv_exposure_model_mod,                                                 &
        ONLY:  uvem_check_data_output
    USE wind_turbine_model_mod,                                                &
        ONLY:  wtm_check_parameters
    USE vertical_nesting_mod,                                                  &
        ONLY:  vnested, vnest_check_parameters


    IMPLICIT NONE

    CHARACTER (LEN=varnamelength)  ::  var           !< variable name
    CHARACTER (LEN=7)   ::  unit                     !< unit of variable
    CHARACTER (LEN=8)   ::  date                     !< current date string
    CHARACTER (LEN=10)  ::  time                     !< current time string
    CHARACTER (LEN=20)  ::  ensemble_string          !< string containing number of ensemble member
    CHARACTER (LEN=15)  ::  nest_string              !< string containing id of nested domain
    CHARACTER (LEN=40)  ::  coupling_string          !< string containing type of coupling
    CHARACTER (LEN=100) ::  action                   !< flag string

    INTEGER(iwp) ::  i                               !< loop index
    INTEGER(iwp) ::  ilen                            !< string length
    INTEGER(iwp) ::  j                               !< loop index
    INTEGER(iwp) ::  k                               !< loop index
    INTEGER(iwp) ::  kk                              !< loop index
    INTEGER(iwp) ::  netcdf_data_format_save         !< initial value of netcdf_data_format
    INTEGER(iwp) ::  position                        !< index position of string
    INTEGER(iwp) ::  lsp                             !< running index for chem spcs.

    LOGICAL     ::  found                            !< flag, true if output variable is already marked for averaging

    REAL(wp)    ::  dt_spinup_max                    !< maximum spinup timestep in nested domains
    REAL(wp)    ::  dum                              !< dummy variable
    REAL(wp)    ::  gradient                         !< local gradient
    REAL(wp)    ::  remote = 0.0_wp                  !< MPI id of remote processor
    REAL(wp)    ::  spinup_time_max                  !< maximum spinup time in nested domains
    REAL(wp)    ::  time_to_be_simulated_from_reference_point  !< time to be simulated from reference point


    CALL location_message( 'checking parameters', .FALSE. )

!
!-- At first, check static and dynamic input for consistency
    CALL netcdf_data_input_check_dynamic
    CALL netcdf_data_input_check_static
!
!-- Check for overlap combinations, which are not realized yet
    IF ( transpose_compute_overlap )  THEN
       IF ( numprocs == 1 )  STOP '+++ transpose-compute-overlap not implemented for single PE runs'
    ENDIF

!
!-- Check the coupling mode
    IF ( coupling_mode /= 'uncoupled'            .AND.                         &
         coupling_mode /= 'precursor_atmos'      .AND.                         &
         coupling_mode /= 'precursor_ocean'      .AND.                         &
         coupling_mode /= 'vnested_crse'         .AND.                         &
         coupling_mode /= 'vnested_fine'         .AND.                         &
         coupling_mode /= 'atmosphere_to_ocean'  .AND.                         &
         coupling_mode /= 'ocean_to_atmosphere' )  THEN
       message_string = 'illegal coupling mode: ' // TRIM( coupling_mode )
       CALL message( 'check_parameters', 'PA0002', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check if humidity is set to TRUE in case of the atmospheric run (for coupled runs)
    IF ( coupling_mode == 'atmosphere_to_ocean' .AND. .NOT. humidity) THEN
       message_string = ' Humidity has to be set to .T. in the _p3d file ' //  &
                        'for coupled runs between ocean and atmosphere.'
       CALL message( 'check_parameters', 'PA0476', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check dt_coupling, restart_time, dt_restart, end_time, dx, dy, nx and ny
    IF ( coupling_mode /= 'uncoupled'       .AND.                              &
         coupling_mode(1:8) /= 'vnested_'   .AND.                              &
         coupling_mode /= 'precursor_atmos' .AND.                              &
         coupling_mode /= 'precursor_ocean' )  THEN

       IF ( dt_coupling == 9999999.9_wp )  THEN
          message_string = 'dt_coupling is not set but required for coup' //   &
                           'ling mode "' //  TRIM( coupling_mode ) // '"'
          CALL message( 'check_parameters', 'PA0003', 1, 2, 0, 6, 0 )
       ENDIF

#if defined( __parallel )


       IF ( myid == 0 ) THEN
          CALL MPI_SEND( dt_coupling, 1, MPI_REAL, target_id, 11, comm_inter,  &
                         ierr )
          CALL MPI_RECV( remote, 1, MPI_REAL, target_id, 11, comm_inter,       &
                         status, ierr )
       ENDIF
       CALL MPI_BCAST( remote, 1, MPI_REAL, 0, comm2d, ierr)

       IF ( dt_coupling /= remote )  THEN
          WRITE( message_string, * ) 'coupling mode "', TRIM( coupling_mode ), &
                 '": dt_coupling = ', dt_coupling, '& is not equal to ',       &
                 'dt_coupling_remote = ', remote
          CALL message( 'check_parameters', 'PA0004', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( dt_coupling <= 0.0_wp )  THEN

          IF ( myid == 0  ) THEN
             CALL MPI_SEND( dt_max, 1, MPI_REAL, target_id, 19, comm_inter, ierr )
             CALL MPI_RECV( remote, 1, MPI_REAL, target_id, 19, comm_inter,    &
                            status, ierr )
          ENDIF
          CALL MPI_BCAST( remote, 1, MPI_REAL, 0, comm2d, ierr)

          dt_coupling = MAX( dt_max, remote )
          WRITE( message_string, * ) 'coupling mode "', TRIM( coupling_mode ), &
                 '": dt_coupling <= 0.0 & is not allowed and is reset to ',    &
                 'MAX(dt_max(A,O)) = ', dt_coupling
          CALL message( 'check_parameters', 'PA0005', 0, 1, 0, 6, 0 )
       ENDIF

       IF ( myid == 0 ) THEN
          CALL MPI_SEND( restart_time, 1, MPI_REAL, target_id, 12, comm_inter, &
                         ierr )
          CALL MPI_RECV( remote, 1, MPI_REAL, target_id, 12, comm_inter,       &
                         status, ierr )
       ENDIF
       CALL MPI_BCAST( remote, 1, MPI_REAL, 0, comm2d, ierr)

       IF ( restart_time /= remote )  THEN
          WRITE( message_string, * ) 'coupling mode "', TRIM( coupling_mode ), &
                 '": restart_time = ', restart_time, '& is not equal to ',     &
                 'restart_time_remote = ', remote
          CALL message( 'check_parameters', 'PA0006', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( myid == 0 ) THEN
          CALL MPI_SEND( dt_restart, 1, MPI_REAL, target_id, 13, comm_inter,   &
                         ierr )
          CALL MPI_RECV( remote, 1, MPI_REAL, target_id, 13, comm_inter,       &
                         status, ierr )
       ENDIF
       CALL MPI_BCAST( remote, 1, MPI_REAL, 0, comm2d, ierr)

       IF ( dt_restart /= remote )  THEN
          WRITE( message_string, * ) 'coupling mode "', TRIM( coupling_mode ), &
                 '": dt_restart = ', dt_restart, '& is not equal to ',         &
                 'dt_restart_remote = ', remote
          CALL message( 'check_parameters', 'PA0007', 1, 2, 0, 6, 0 )
       ENDIF

       time_to_be_simulated_from_reference_point = end_time-coupling_start_time

       IF  ( myid == 0 ) THEN
          CALL MPI_SEND( time_to_be_simulated_from_reference_point, 1,         &
                         MPI_REAL, target_id, 14, comm_inter, ierr )
          CALL MPI_RECV( remote, 1, MPI_REAL, target_id, 14, comm_inter,       &
                         status, ierr )
       ENDIF
       CALL MPI_BCAST( remote, 1, MPI_REAL, 0, comm2d, ierr)

       IF ( time_to_be_simulated_from_reference_point /= remote )  THEN
          WRITE( message_string, * ) 'coupling mode "', TRIM( coupling_mode ), &
                 '": time_to_be_simulated_from_reference_point = ',            &
                 time_to_be_simulated_from_reference_point, '& is not equal ', &
                 'to time_to_be_simulated_from_reference_point_remote = ',     &
                 remote
          CALL message( 'check_parameters', 'PA0008', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( myid == 0 ) THEN
          CALL MPI_SEND( dx, 1, MPI_REAL, target_id, 15, comm_inter, ierr )
          CALL MPI_RECV( remote, 1, MPI_REAL, target_id, 15, comm_inter,       &
                                                             status, ierr )
       ENDIF
       CALL MPI_BCAST( remote, 1, MPI_REAL, 0, comm2d, ierr)


       IF ( coupling_mode == 'atmosphere_to_ocean') THEN

          IF ( dx < remote ) THEN
             WRITE( message_string, * ) 'coupling mode "',                     &
                   TRIM( coupling_mode ),                                      &
           '": dx in Atmosphere is not equal to or not larger than dx in ocean'
             CALL message( 'check_parameters', 'PA0009', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( (nx_a+1)*dx /= (nx_o+1)*remote )  THEN
             WRITE( message_string, * ) 'coupling mode "',                     &
                    TRIM( coupling_mode ),                                     &
             '": Domain size in x-direction is not equal in ocean and atmosphere'
             CALL message( 'check_parameters', 'PA0010', 1, 2, 0, 6, 0 )
          ENDIF

       ENDIF

       IF ( myid == 0) THEN
          CALL MPI_SEND( dy, 1, MPI_REAL, target_id, 16, comm_inter, ierr )
          CALL MPI_RECV( remote, 1, MPI_REAL, target_id, 16, comm_inter,       &
                         status, ierr )
       ENDIF
       CALL MPI_BCAST( remote, 1, MPI_REAL, 0, comm2d, ierr)

       IF ( coupling_mode == 'atmosphere_to_ocean') THEN

          IF ( dy < remote )  THEN
             WRITE( message_string, * ) 'coupling mode "',                     &
                    TRIM( coupling_mode ),                                     &
                 '": dy in Atmosphere is not equal to or not larger than dy in ocean'
             CALL message( 'check_parameters', 'PA0011', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( (ny_a+1)*dy /= (ny_o+1)*remote )  THEN
             WRITE( message_string, * ) 'coupling mode "',                     &
                   TRIM( coupling_mode ),                                      &
             '": Domain size in y-direction is not equal in ocean and atmosphere'
             CALL message( 'check_parameters', 'PA0012', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( MOD(nx_o+1,nx_a+1) /= 0 )  THEN
             WRITE( message_string, * ) 'coupling mode "',                     &
                   TRIM( coupling_mode ),                                      &
             '": nx+1 in ocean is not divisible by nx+1 in',                   &
             ' atmosphere without remainder'
             CALL message( 'check_parameters', 'PA0339', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( MOD(ny_o+1,ny_a+1) /= 0 )  THEN
             WRITE( message_string, * ) 'coupling mode "',                     &
                   TRIM( coupling_mode ),                                      &
             '": ny+1 in ocean is not divisible by ny+1 in', &
             ' atmosphere without remainder'

             CALL message( 'check_parameters', 'PA0340', 1, 2, 0, 6, 0 )
          ENDIF

       ENDIF
#else
       WRITE( message_string, * ) 'coupling requires PALM to be compiled with',&
            ' cpp-option "-D__parallel"'
       CALL message( 'check_parameters', 'PA0141', 1, 2, 0, 6, 0 )
#endif
    ENDIF

#if defined( __parallel )
!
!-- Exchange via intercommunicator
    IF ( coupling_mode == 'atmosphere_to_ocean' .AND. myid == 0 )  THEN
       CALL MPI_SEND( humidity, 1, MPI_LOGICAL, target_id, 19, comm_inter,     &
                      ierr )
    ELSEIF ( coupling_mode == 'ocean_to_atmosphere' .AND. myid == 0)  THEN
       CALL MPI_RECV( humidity_remote, 1, MPI_LOGICAL, target_id, 19,          &
                      comm_inter, status, ierr )
    ENDIF
    CALL MPI_BCAST( humidity_remote, 1, MPI_LOGICAL, 0, comm2d, ierr)

#endif

!
!-- User settings for restart times requires that "restart" has been given as
!-- file activation string. Otherwise, binary output would not be saved by
!-- palmrun.
    IF (  ( restart_time /= 9999999.9_wp  .OR.  dt_restart /= 9999999.9_wp )   &
         .AND.  .NOT. write_binary )  THEN
       WRITE( message_string, * ) 'manual restart settings requires file ',    &
                                  'activation string "restart"'
       CALL message( 'check_parameters', 'PA0001', 1, 2, 0, 6, 0 )
    ENDIF


!
!-- Generate the file header which is used as a header for most of PALM's
!-- output files
    CALL DATE_AND_TIME( date, time )
    run_date = date(7:8)//'-'//date(5:6)//'-'//date(3:4)
    run_time = time(1:2)//':'//time(3:4)//':'//time(5:6)
    IF ( coupling_mode == 'uncoupled' )  THEN
       coupling_string = ''
    ELSEIF ( coupling_mode == 'vnested_crse' )  THEN
       coupling_string = ' nested (coarse)'
    ELSEIF ( coupling_mode == 'vnested_fine' )  THEN
       coupling_string = ' nested (fine)'
    ELSEIF ( coupling_mode == 'atmosphere_to_ocean' )  THEN
       coupling_string = ' coupled (atmosphere)'
    ELSEIF ( coupling_mode == 'ocean_to_atmosphere' )  THEN
       coupling_string = ' coupled (ocean)'
    ENDIF
    IF ( ensemble_member_nr /= 0 )  THEN
       WRITE( ensemble_string, '(2X,A,I2.2)' )  'en-no: ', ensemble_member_nr
    ELSE
       ensemble_string = ''
    ENDIF
    IF ( nested_run )  THEN
       WRITE( nest_string, '(2X,A,I2.2)' )  'nest-id: ', cpl_id
    ELSE
       nest_string = ''
    ENDIF

    WRITE ( run_description_header,                                            &
            '(A,2X,A,2X,A,A,A,I2.2,A,A,A,2X,A,A,2X,A,1X,A)' )                  &
          TRIM( version ), TRIM( revision ), 'run: ',                          &
          TRIM( run_identifier ), '.', runnr, TRIM( coupling_string ),         &
          TRIM( nest_string ), TRIM( ensemble_string), 'host: ', TRIM( host ), &
          run_date, run_time

!
!-- Check the general loop optimization method
    SELECT CASE ( TRIM( loop_optimization ) )

       CASE ( 'cache', 'vector' )
          CONTINUE

       CASE DEFAULT
          message_string = 'illegal value given for loop_optimization: "' //   &
                           TRIM( loop_optimization ) // '"'
          CALL message( 'check_parameters', 'PA0013', 1, 2, 0, 6, 0 )

    END SELECT

!
!-- Check topography setting (check for illegal parameter combinations)
    IF ( topography /= 'flat' )  THEN
       action = ' '
       IF ( scalar_advec /= 'pw-scheme' .AND. scalar_advec /= 'ws-scheme'      &
          )  THEN
          WRITE( action, '(A,A)' )  'scalar_advec = ', scalar_advec
       ENDIF
       IF ( momentum_advec /= 'pw-scheme' .AND. momentum_advec /= 'ws-scheme' )&
       THEN
          WRITE( action, '(A,A)' )  'momentum_advec = ', momentum_advec
       ENDIF
       IF ( psolver == 'sor' )  THEN
          WRITE( action, '(A,A)' )  'psolver = ', psolver
       ENDIF
       IF ( sloping_surface )  THEN
          WRITE( action, '(A)' )  'sloping surface = .TRUE.'
       ENDIF
       IF ( galilei_transformation )  THEN
          WRITE( action, '(A)' )  'galilei_transformation = .TRUE.'
       ENDIF
       IF ( cloud_physics )  THEN
          WRITE( action, '(A)' )  'cloud_physics = .TRUE.'
       ENDIF
       IF ( cloud_droplets )  THEN
          WRITE( action, '(A)' )  'cloud_droplets = .TRUE.'
       ENDIF
       IF ( .NOT. constant_flux_layer )  THEN
          WRITE( action, '(A)' )  'constant_flux_layer = .FALSE.'
       ENDIF
       IF ( action /= ' ' )  THEN
          message_string = 'a non-flat topography does not allow ' //          &
                           TRIM( action )
          CALL message( 'check_parameters', 'PA0014', 1, 2, 0, 6, 0 )
       ENDIF

    ENDIF

!
!-- Check turbulence closure setup
    CALL tcm_check_parameters
!
!-- Check approximation
    IF ( TRIM( approximation ) /= 'boussinesq'   .AND.                         &
         TRIM( approximation ) /= 'anelastic' )  THEN
       message_string = 'unknown approximation: approximation = "' //          &
                        TRIM( approximation ) // '"'
       CALL message( 'check_parameters', 'PA0446', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check approximation requirements
    IF ( TRIM( approximation ) == 'anelastic'   .AND.                          &
         TRIM( momentum_advec ) /= 'ws-scheme' )  THEN
       message_string = 'Anelastic approximation requires ' //                 &
                        'momentum_advec = "ws-scheme"'
       CALL message( 'check_parameters', 'PA0447', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( TRIM( approximation ) == 'anelastic'   .AND.                          &
         TRIM( psolver ) == 'multigrid' )  THEN
       message_string = 'Anelastic approximation currently only supports ' //  &
                        'psolver = "poisfft", ' //                             &
                        'psolver = "sor" and ' //                              &
                        'psolver = "multigrid_noopt"'
       CALL message( 'check_parameters', 'PA0448', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( TRIM( approximation ) == 'anelastic'   .AND.                          &
         conserve_volume_flow )  THEN
       message_string = 'Anelastic approximation is not allowed with ' //      &
                        'conserve_volume_flow = .TRUE.'
       CALL message( 'check_parameters', 'PA0449', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check flux input mode
    IF ( TRIM( flux_input_mode ) /= 'dynamic'    .AND.                         &
         TRIM( flux_input_mode ) /= 'kinematic'  .AND.                         &
         TRIM( flux_input_mode ) /= 'approximation-specific' )  THEN
       message_string = 'unknown flux input mode: flux_input_mode = "' //      &
                        TRIM( flux_input_mode ) // '"'
       CALL message( 'check_parameters', 'PA0450', 1, 2, 0, 6, 0 )
    ENDIF
!-- Set flux input mode according to approximation if applicable
    IF ( TRIM( flux_input_mode ) == 'approximation-specific' )  THEN
       IF ( TRIM( approximation ) == 'anelastic' )  THEN
          flux_input_mode = 'dynamic'
       ELSEIF ( TRIM( approximation ) == 'boussinesq' )  THEN
          flux_input_mode = 'kinematic'
       ENDIF
    ENDIF

!
!-- Check flux output mode
    IF ( TRIM( flux_output_mode ) /= 'dynamic'    .AND.                        &
         TRIM( flux_output_mode ) /= 'kinematic'  .AND.                        &
         TRIM( flux_output_mode ) /= 'approximation-specific' )  THEN
       message_string = 'unknown flux output mode: flux_output_mode = "' //    &
                        TRIM( flux_output_mode ) // '"'
       CALL message( 'check_parameters', 'PA0451', 1, 2, 0, 6, 0 )
    ENDIF
!-- Set flux output mode according to approximation if applicable
    IF ( TRIM( flux_output_mode ) == 'approximation-specific' )  THEN
       IF ( TRIM( approximation ) == 'anelastic' )  THEN
          flux_output_mode = 'dynamic'
       ELSEIF ( TRIM( approximation ) == 'boussinesq' )  THEN
          flux_output_mode = 'kinematic'
       ENDIF
    ENDIF


!
!-- When the land- or urban-surface model is used, the flux output must be
!-- dynamic.
    IF ( land_surface  .OR.  urban_surface )  THEN
       flux_output_mode = 'dynamic'
    ENDIF

!
!-- set the flux output units according to flux_output_mode
    IF ( TRIM( flux_output_mode ) == 'kinematic' ) THEN
        heatflux_output_unit              = 'K m/s'
        waterflux_output_unit             = 'kg/kg m/s'
        momentumflux_output_unit          = 'm2/s2'
    ELSEIF ( TRIM( flux_output_mode ) == 'dynamic' ) THEN
        heatflux_output_unit              = 'W/m2'
        waterflux_output_unit             = 'W/m2'
        momentumflux_output_unit          = 'N/m2'
    ENDIF

!-- set time series output units for fluxes
    dots_unit(14:16) = heatflux_output_unit
    dots_unit(21)    = waterflux_output_unit
    dots_unit(19:20) = momentumflux_output_unit

!
!-- Check ocean setting
    IF ( TRIM( coupling_mode ) == 'uncoupled'  .AND.                           &
         TRIM( coupling_char ) == '_O' .AND.                                   &
         .NOT. ocean )  THEN

!
!--    Check whether an (uncoupled) atmospheric run has been declared as an
!--    ocean run (this setting is done via mrun-option -y)
       message_string = 'ocean = .F. does not allow coupling_char = "' //      &
                        TRIM( coupling_char ) // '" set by palmrun-option "-y"'
       CALL message( 'check_parameters', 'PA0317', 1, 2, 0, 6, 0 )

    ENDIF
!
!-- Check cloud scheme
    IF ( cloud_scheme == 'saturation_adjust' )  THEN
       microphysics_sat_adjust = .TRUE.
       microphysics_seifert    = .FALSE.
       microphysics_kessler    = .FALSE.
       precipitation           = .FALSE.
    ELSEIF ( cloud_scheme == 'seifert_beheng' )  THEN
       microphysics_sat_adjust = .FALSE.
       microphysics_seifert    = .TRUE.
       microphysics_kessler    = .FALSE.
       microphysics_morrison  = .FALSE.
       precipitation           = .TRUE.
    ELSEIF ( cloud_scheme == 'kessler' )  THEN
       microphysics_sat_adjust = .FALSE.
       microphysics_seifert    = .FALSE.
       microphysics_kessler    = .TRUE.
       microphysics_morrison   = .FALSE.
       precipitation           = .TRUE.
    ELSEIF ( cloud_scheme == 'morrison' )  THEN
       microphysics_sat_adjust = .FALSE.
       microphysics_seifert    = .TRUE.
       microphysics_kessler    = .FALSE.
       microphysics_morrison   = .TRUE.
       precipitation           = .TRUE.
    ELSE
       message_string = 'unknown cloud microphysics scheme cloud_scheme ="' // &
                        TRIM( cloud_scheme ) // '"'
       CALL message( 'check_parameters', 'PA0357', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check aerosol
    IF ( aerosol_bulk == 'nacl' )  THEN
       aerosol_nacl   = .TRUE.
       aerosol_c3h4o4 = .FALSE.
       aerosol_nh4no3 = .FALSE.
    ELSEIF ( aerosol_bulk == 'c3h4o4' )  THEN
       aerosol_nacl   = .FALSE.
       aerosol_c3h4o4 = .TRUE.
       aerosol_nh4no3 = .FALSE.
    ELSEIF ( aerosol_bulk == 'nh4no3' )  THEN
       aerosol_nacl   = .FALSE.
       aerosol_c3h4o4 = .FALSE.
       aerosol_nh4no3 = .TRUE.
    ELSE
       message_string = 'unknown aerosol = "' // TRIM( aerosol_bulk ) // '"'
       CALL message( 'check_parameters', 'PA0469', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check whether there are any illegal values
!-- Pressure solver:
    IF ( psolver /= 'poisfft'  .AND.  psolver /= 'sor'  .AND.                  &
         psolver /= 'multigrid'  .AND.  psolver /= 'multigrid_noopt' )  THEN
       message_string = 'unknown solver for perturbation pressure: psolver' // &
                        ' = "' // TRIM( psolver ) // '"'
       CALL message( 'check_parameters', 'PA0016', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( psolver(1:9) == 'multigrid' )  THEN
       IF ( cycle_mg == 'w' )  THEN
          gamma_mg = 2
       ELSEIF ( cycle_mg == 'v' )  THEN
          gamma_mg = 1
       ELSE
          message_string = 'unknown multigrid cycle: cycle_mg = "' //          &
                           TRIM( cycle_mg ) // '"'
          CALL message( 'check_parameters', 'PA0020', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( fft_method /= 'singleton-algorithm'  .AND.                            &
         fft_method /= 'temperton-algorithm'  .AND.                            &
         fft_method /= 'fftw'                 .AND.                            &
         fft_method /= 'system-specific' )  THEN
       message_string = 'unknown fft-algorithm: fft_method = "' //             &
                        TRIM( fft_method ) // '"'
       CALL message( 'check_parameters', 'PA0021', 1, 2, 0, 6, 0 )
    ENDIF

    IF( momentum_advec == 'ws-scheme' .AND.                                    &
        .NOT. call_psolver_at_all_substeps  ) THEN
        message_string = 'psolver must be called at each RK3 substep when "'// &
                      TRIM(momentum_advec) // ' "is used for momentum_advec'
        CALL message( 'check_parameters', 'PA0344', 1, 2, 0, 6, 0 )
    END IF
!
!-- Advection schemes:
    IF ( momentum_advec /= 'pw-scheme'  .AND.                                  &
         momentum_advec /= 'ws-scheme'  .AND.                                  &
         momentum_advec /= 'up-scheme' )                                       &
    THEN
       message_string = 'unknown advection scheme: momentum_advec = "' //      &
                        TRIM( momentum_advec ) // '"'
       CALL message( 'check_parameters', 'PA0022', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( ( momentum_advec == 'ws-scheme' .OR.  scalar_advec == 'ws-scheme' )   &
           .AND. ( timestep_scheme == 'euler' .OR.                             &
                   timestep_scheme == 'runge-kutta-2' ) )                      &
    THEN
       message_string = 'momentum_advec or scalar_advec = "'                   &
         // TRIM( momentum_advec ) // '" is not allowed with ' //              &
         'timestep_scheme = "' // TRIM( timestep_scheme ) // '"'
       CALL message( 'check_parameters', 'PA0023', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( scalar_advec /= 'pw-scheme'  .AND.  scalar_advec /= 'ws-scheme' .AND. &
         scalar_advec /= 'bc-scheme' .AND. scalar_advec /= 'up-scheme' )       &
    THEN
       message_string = 'unknown advection scheme: scalar_advec = "' //        &
                        TRIM( scalar_advec ) // '"'
       CALL message( 'check_parameters', 'PA0024', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( scalar_advec == 'bc-scheme'  .AND.  loop_optimization == 'cache' )    &
    THEN
       message_string = 'advection_scheme scalar_advec = "'                    &
         // TRIM( scalar_advec ) // '" not implemented for ' //                &
         'loop_optimization = "' // TRIM( loop_optimization ) // '"'
       CALL message( 'check_parameters', 'PA0026', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( use_sgs_for_particles  .AND.  .NOT. cloud_droplets  .AND.             &
         .NOT. use_upstream_for_tke  .AND.                                     &
         scalar_advec /= 'ws-scheme'                                           &
       )  THEN
       use_upstream_for_tke = .TRUE.
       message_string = 'use_upstream_for_tke is set to .TRUE. because ' //    &
                        'use_sgs_for_particles = .TRUE. '          //          &
                        'and scalar_advec /= ws-scheme'
       CALL message( 'check_parameters', 'PA0025', 0, 1, 0, 6, 0 )
    ENDIF

!
!-- Set LOGICAL switches to enhance performance
    IF ( momentum_advec == 'ws-scheme' )  ws_scheme_mom = .TRUE.
    IF ( scalar_advec   == 'ws-scheme' )  ws_scheme_sca = .TRUE.


!
!-- Timestep schemes:
    SELECT CASE ( TRIM( timestep_scheme ) )

       CASE ( 'euler' )
          intermediate_timestep_count_max = 1

       CASE ( 'runge-kutta-2' )
          intermediate_timestep_count_max = 2

       CASE ( 'runge-kutta-3' )
          intermediate_timestep_count_max = 3

       CASE DEFAULT
          message_string = 'unknown timestep scheme: timestep_scheme = "' //   &
                           TRIM( timestep_scheme ) // '"'
          CALL message( 'check_parameters', 'PA0027', 1, 2, 0, 6, 0 )

    END SELECT

    IF ( (momentum_advec /= 'pw-scheme' .AND. momentum_advec /= 'ws-scheme')   &
         .AND. timestep_scheme(1:5) == 'runge' ) THEN
       message_string = 'momentum advection scheme "' // &
                        TRIM( momentum_advec ) // '" & does not work with ' // &
                        'timestep_scheme "' // TRIM( timestep_scheme ) // '"'
       CALL message( 'check_parameters', 'PA0029', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check for proper settings for microphysics
    IF ( cloud_physics  .AND.  cloud_droplets )  THEN
       message_string = 'cloud_physics = .TRUE. is not allowed with ' //       &
                        'cloud_droplets = .TRUE.'
       CALL message( 'check_parameters', 'PA0442', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Collision kernels:
    SELECT CASE ( TRIM( collision_kernel ) )

       CASE ( 'hall', 'hall_fast' )
          hall_kernel = .TRUE.

       CASE ( 'wang', 'wang_fast' )
          wang_kernel = .TRUE.

       CASE ( 'none' )


       CASE DEFAULT
          message_string = 'unknown collision kernel: collision_kernel = "' // &
                           TRIM( collision_kernel ) // '"'
          CALL message( 'check_parameters', 'PA0350', 1, 2, 0, 6, 0 )

    END SELECT
    IF ( collision_kernel(6:9) == 'fast' )  use_kernel_tables = .TRUE.

!
!-- Initializing actions must have been set by the user

    WRITE(message_string,*) 'initializing_actions = ',initializing_actions
    CALL location_message(adjustl(trim(message_string)),.TRUE.)

    IF ( TRIM( initializing_actions ) == '' )  THEN
       message_string = 'no value specified for initializing_actions'
       CALL message( 'check_parameters', 'PA0149', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.            &
         TRIM( initializing_actions ) /= 'cyclic_fill' )  THEN
!
!--    No restart run: several initialising actions are possible
       action = initializing_actions
       DO  WHILE ( TRIM( action ) /= '' )
          position = INDEX( action, ' ' )
          SELECT CASE ( action(1:position-1) )

             CASE ( 'set_constant_profiles', 'set_1d-model_profiles',          &
                    'by_user', 'initialize_vortex', 			       &
                    'initialize_2D_bubble', 'initialize_3D_bubble', 'inifor' )
                action = action(position+1:)

             CASE DEFAULT
                message_string = 'initializing_action = "' //                  &
                                 TRIM( action ) // '" unknown or not allowed'
                CALL message( 'check_parameters', 'PA0030', 1, 2, 0, 6, 0 )

          END SELECT
       ENDDO
    ENDIF

    IF ( TRIM( initializing_actions ) == 'initialize_vortex'  .AND.            &
         conserve_volume_flow ) THEN
         message_string = 'initializing_actions = "initialize_vortex"' //      &
                        ' is not allowed with conserve_volume_flow = .T.'
       CALL message( 'check_parameters', 'PA0343', 1, 2, 0, 6, 0 )
    ENDIF


    IF ( INDEX( initializing_actions, 'set_constant_profiles' ) /= 0  .AND.    &
         INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN
       message_string = 'initializing_actions = "set_constant_profiles"' //    &
                        ' and "set_1d-model_profiles" are not allowed ' //     &
                        'simultaneously'
       CALL message( 'check_parameters', 'PA0031', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( INDEX( initializing_actions, 'set_constant_profiles' ) /= 0  .AND.    &
         INDEX( initializing_actions, 'by_user' ) /= 0 )  THEN
       message_string = 'initializing_actions = "set_constant_profiles"' //    &
                        ' and "by_user" are not allowed simultaneously'
       CALL message( 'check_parameters', 'PA0032', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( INDEX( initializing_actions, 'by_user' ) /= 0  .AND.                  &
         INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN
       message_string = 'initializing_actions = "by_user" and ' //             &
                        '"set_1d-model_profiles" are not allowed simultaneously'
       CALL message( 'check_parameters', 'PA0033', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- In case of spinup and nested run, spinup end time must be identical
!-- in order to have synchronously running simulations.
    IF ( nested_run )  THEN
#if defined( __parallel )
       CALL MPI_ALLREDUCE( spinup_time, spinup_time_max, 1, MPI_REAL,          &
                           MPI_MAX, MPI_COMM_WORLD, ierr )
       CALL MPI_ALLREDUCE( dt_spinup,   dt_spinup_max,   1, MPI_REAL,          &
                           MPI_MAX, MPI_COMM_WORLD, ierr )
       IF ( spinup_time /= spinup_time_max  .OR.  dt_spinup /= dt_spinup_max ) &
       THEN
          message_string = 'In case of nesting, spinup_time and ' //           &
                           'dt_spinup must be identical in all parent ' //     &
                           'and child domains.'
          CALL message( 'check_parameters', 'PA0489', 3, 2, 0, 6, 0 )
       ENDIF
#endif
    ENDIF

    IF ( cloud_physics  .AND.  .NOT.  humidity )  THEN
       WRITE( message_string, * ) 'cloud_physics = ', cloud_physics, ' is ',   &
              'not allowed with humidity = ', humidity
       CALL message( 'check_parameters', 'PA0034', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( humidity  .AND.  sloping_surface )  THEN
       message_string = 'humidity = .TRUE. and sloping_surface = .TRUE. ' //   &
                        'are not allowed simultaneously'
       CALL message( 'check_parameters', 'PA0036', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check for synthetic turbulence generator settings
    CALL stg_check_parameters

!
!-- Check for Stokes drift settings if Stokes forcing in ocean mode is used
    IF ( ocean .AND. stokes_force ) CALL stokes_drift_check_parameters

!
!-- When plant canopy model is used, peform addtional checks
    IF ( plant_canopy )  CALL pcm_check_parameters

!
!-- Additional checks for spectra
    IF ( calculate_spectra )  CALL spectra_check_parameters
!
!-- When land surface model is used, perform additional checks
    IF ( land_surface )  CALL lsm_check_parameters

!
!-- When urban surface model is used, perform additional checks
    IF ( urban_surface )  CALL usm_check_parameters

!
!-- If wind turbine model is used, perform additional checks
    IF ( wind_turbine )  CALL wtm_check_parameters
!
!-- When radiation model is used, peform addtional checks
    IF ( radiation )  CALL radiation_check_parameters
!
!-- When gust module is used, perform additional checks
    IF ( gust_module_enabled )  CALL gust_check_parameters
!
!-- When large-scale forcing or nudging is used, peform addtional checks
    IF ( large_scale_forcing  .OR.  nudging )  CALL lsf_nudging_check_parameters

!
!-- In case of no model continuation run, check initialising parameters and
!-- deduce further quantities
    IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN

!
!--    Initial profiles for 1D and 3D model, respectively (u,v further below)
       pt_init = pt_surface
       IF ( humidity       )  q_init  = q_surface
       IF ( ocean          )  sa_init = sa_surface
       IF ( passive_scalar )  s_init  = s_surface
       IF ( air_chemistry )  THEN
         DO lsp = 1, nvar
            chem_species(lsp)%conc_pr_init = cs_surface(lsp)
         ENDDO
       ENDIF
!
!--
!--    If required, compute initial profile of the geostrophic wind
!--    (component ug)
       i = 1
       gradient = 0.0_wp

       IF (  .NOT.  ocean )  THEN

          ug_vertical_gradient_level_ind(1) = 0
          ug(0) = ug_surface
          DO  k = 1, nzt+1
             IF ( i < 11 )  THEN
                IF ( ug_vertical_gradient_level(i) < zu(k)  .AND.              &
                     ug_vertical_gradient_level(i) >= 0.0_wp )  THEN
                   gradient = ug_vertical_gradient(i) / 100.0_wp
                   ug_vertical_gradient_level_ind(i) = k - 1
                   i = i + 1
                ENDIF
             ENDIF
             IF ( gradient /= 0.0_wp )  THEN
                IF ( k /= 1 )  THEN
                   ug(k) = ug(k-1) + dzu(k) * gradient
                ELSE
                   ug(k) = ug_surface + dzu(k) * gradient
                ENDIF
             ELSE
                ug(k) = ug(k-1)
             ENDIF
          ENDDO

       ELSE

          ug_vertical_gradient_level_ind(1) = nzt+1
          ug(nzt+1) = ug_surface
          DO  k = nzt, nzb, -1
             IF ( i < 11 )  THEN
                IF ( ug_vertical_gradient_level(i) > zu(k)  .AND.              &
                     ug_vertical_gradient_level(i) <= 0.0_wp )  THEN
                   gradient = ug_vertical_gradient(i) / 100.0_wp
                   ug_vertical_gradient_level_ind(i) = k + 1
                   i = i + 1
                ENDIF
             ENDIF
             IF ( gradient /= 0.0_wp )  THEN
                IF ( k /= nzt )  THEN
                   ug(k) = ug(k+1) - dzu(k+1) * gradient
                ELSE
                   ug(k)   = ug_surface - 0.5_wp * dzu(k+1) * gradient
                   ug(k+1) = ug_surface + 0.5_wp * dzu(k+1) * gradient
                ENDIF
             ELSE
                ug(k) = ug(k+1)
             ENDIF
          ENDDO

       ENDIF

!
!--    In case of no given gradients for ug, choose a zero gradient
       IF ( ug_vertical_gradient_level(1) == -9999999.9_wp )  THEN
          ug_vertical_gradient_level(1) = 0.0_wp
       ENDIF

!
!--
!--    If required, compute initial profile of the geostrophic wind
!--    (component vg)
       i = 1
       gradient = 0.0_wp

       IF (  .NOT.  ocean )  THEN

          vg_vertical_gradient_level_ind(1) = 0
          vg(0) = vg_surface
          DO  k = 1, nzt+1
             IF ( i < 11 )  THEN
                IF ( vg_vertical_gradient_level(i) < zu(k)  .AND.              &
                     vg_vertical_gradient_level(i) >= 0.0_wp )  THEN
                   gradient = vg_vertical_gradient(i) / 100.0_wp
                   vg_vertical_gradient_level_ind(i) = k - 1
                   i = i + 1
                ENDIF
             ENDIF
             IF ( gradient /= 0.0_wp )  THEN
                IF ( k /= 1 )  THEN
                   vg(k) = vg(k-1) + dzu(k) * gradient
                ELSE
                   vg(k) = vg_surface + dzu(k) * gradient
                ENDIF
             ELSE
                vg(k) = vg(k-1)
             ENDIF
          ENDDO

       ELSE

          vg_vertical_gradient_level_ind(1) = nzt+1
          vg(nzt+1) = vg_surface
          DO  k = nzt, nzb, -1
             IF ( i < 11 )  THEN
                IF ( vg_vertical_gradient_level(i) > zu(k)  .AND.              &
                     vg_vertical_gradient_level(i) <= 0.0_wp )  THEN
                   gradient = vg_vertical_gradient(i) / 100.0_wp
                   vg_vertical_gradient_level_ind(i) = k + 1
                   i = i + 1
                ENDIF
             ENDIF
             IF ( gradient /= 0.0_wp )  THEN
                IF ( k /= nzt )  THEN
                   vg(k) = vg(k+1) - dzu(k+1) * gradient
                ELSE
                   vg(k)   = vg_surface - 0.5_wp * dzu(k+1) * gradient
                   vg(k+1) = vg_surface + 0.5_wp * dzu(k+1) * gradient
                ENDIF
             ELSE
                vg(k) = vg(k+1)
             ENDIF
          ENDDO

       ENDIF

!
!--    In case of no given gradients for vg, choose a zero gradient
       IF ( vg_vertical_gradient_level(1) == -9999999.9_wp )  THEN
          vg_vertical_gradient_level(1) = 0.0_wp
       ENDIF

!
!--    Let the initial wind profiles be the calculated ug/vg profiles or
!--    interpolate them from wind profile data (if given)
       IF ( u_profile(1) == 9999999.9_wp  .AND.  v_profile(1) == 9999999.9_wp )  THEN

          u_init = ug
          v_init = vg

       ELSEIF ( u_profile(1) == 0.0_wp  .AND.  v_profile(1) == 0.0_wp )  THEN

          IF ( uv_heights(1) /= 0.0_wp )  THEN
             message_string = 'uv_heights(1) must be 0.0'
             CALL message( 'check_parameters', 'PA0345', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( omega /= 0.0_wp )  THEN
             message_string = 'Coriolis force must be switched off (by setting omega=0.0)' //  &
                              ' when prescribing the forcing by u_profile and v_profile'
             CALL message( 'check_parameters', 'PA0347', 1, 2, 0, 6, 0 )
          ENDIF

          use_prescribed_profile_data = .TRUE.

          kk = 1
          u_init(0) = 0.0_wp
          v_init(0) = 0.0_wp

          DO  k = 1, nz+1

             IF ( kk < 100 )  THEN
                DO  WHILE ( uv_heights(kk+1) <= zu(k) )
                   kk = kk + 1
                   IF ( kk == 100 )  EXIT
                ENDDO
             ENDIF

             IF ( kk < 100  .AND.  uv_heights(kk+1) /= 9999999.9_wp )  THEN
                u_init(k) = u_profile(kk) + ( zu(k) - uv_heights(kk) ) /       &
                                       ( uv_heights(kk+1) - uv_heights(kk) ) * &
                                       ( u_profile(kk+1) - u_profile(kk) )
                v_init(k) = v_profile(kk) + ( zu(k) - uv_heights(kk) ) /       &
                                       ( uv_heights(kk+1) - uv_heights(kk) ) * &
                                       ( v_profile(kk+1) - v_profile(kk) )
             ELSE
                u_init(k) = u_profile(kk)
                v_init(k) = v_profile(kk)
             ENDIF

          ENDDO

       ELSE

          message_string = 'u_profile(1) and v_profile(1) must be 0.0'
          CALL message( 'check_parameters', 'PA0346', 1, 2, 0, 6, 0 )

       ENDIF

!
!--    Compute initial temperature profile using the given temperature gradients
       IF (  .NOT.  neutral )  THEN
          CALL init_vertical_profiles( pt_vertical_gradient_level_ind,          &
                                       pt_vertical_gradient_level,              &
                                       pt_vertical_gradient, pt_init,           &
                                       pt_surface, bc_pt_t_val )
       ENDIF
!
!--    Compute initial humidity profile using the given humidity gradients
       IF ( humidity )  THEN
          CALL init_vertical_profiles( q_vertical_gradient_level_ind,          &
                                       q_vertical_gradient_level,              &
                                       q_vertical_gradient, q_init,            &
                                       q_surface, bc_q_t_val )
       ENDIF
!
!--    Compute initial scalar profile using the given scalar gradients
       IF ( passive_scalar )  THEN
          CALL init_vertical_profiles( s_vertical_gradient_level_ind,          &
                                       s_vertical_gradient_level,              &
                                       s_vertical_gradient, s_init,            &
                                       s_surface, bc_s_t_val )
       ENDIF
!
!--    Compute initial chemistry profile using the given chemical species gradients
       IF ( air_chemistry )  THEN
         DO lsp = 1, nvar
          CALL init_vertical_profiles( cs_vertical_gradient_level_ind(lsp,:),  &
                                       cs_vertical_gradient_level(lsp,:),      &
                                       cs_vertical_gradient(lsp,:),            &
                                       chem_species(lsp)%conc_pr_init,         &
                                       cs_surface(lsp), bc_cs_t_val(lsp) )
         ENDDO
       ENDIF
!
!
!--    If required, compute initial salinity profile using the given salinity
!--    gradients
       IF ( ocean )  THEN
          CALL init_vertical_profiles( sa_vertical_gradient_level_ind,          &
                                       sa_vertical_gradient_level,              &
                                       sa_vertical_gradient, sa_init,           &
                                       sa_surface, dum )
       ENDIF


    ENDIF

!
!-- Check if the control parameter use_subsidence_tendencies is used correctly
    IF ( use_subsidence_tendencies  .AND.  .NOT.  large_scale_subsidence )  THEN
       message_string = 'The usage of use_subsidence_tendencies ' //           &
                            'requires large_scale_subsidence = .T..'
       CALL message( 'check_parameters', 'PA0396', 1, 2, 0, 6, 0 )
    ELSEIF ( use_subsidence_tendencies  .AND.  .NOT. large_scale_forcing )  THEN
       message_string = 'The usage of use_subsidence_tendencies ' //           &
                            'requires large_scale_forcing = .T..'
       CALL message( 'check_parameters', 'PA0397', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Initialize large scale subsidence if required
    If ( large_scale_subsidence )  THEN
       IF ( subs_vertical_gradient_level(1) /= -9999999.9_wp  .AND.            &
                                     .NOT.  large_scale_forcing )  THEN
          CALL init_w_subsidence
       ENDIF
!
!--    In case large_scale_forcing is used, profiles for subsidence velocity
!--    are read in from file LSF_DATA

       IF ( subs_vertical_gradient_level(1) == -9999999.9_wp  .AND.            &
            .NOT.  large_scale_forcing )  THEN
          message_string = 'There is no default large scale vertical ' //      &
                           'velocity profile set. Specify the subsidence ' //  &
                           'velocity profile via subs_vertical_gradient ' //   &
                           'and subs_vertical_gradient_level.'
          CALL message( 'check_parameters', 'PA0380', 1, 2, 0, 6, 0 )
       ENDIF
    ELSE
        IF ( subs_vertical_gradient_level(1) /= -9999999.9_wp )  THEN
           message_string = 'Enable usage of large scale subsidence by ' //    &
                            'setting large_scale_subsidence = .T..'
          CALL message( 'check_parameters', 'PA0381', 1, 2, 0, 6, 0 )
        ENDIF
    ENDIF

!
!-- Overwrite latitude if necessary and compute Coriolis parameter.
!-- To do - move initialization of f and fs to coriolis_mod.
    IF ( input_pids_static )  THEN
       latitude  = init_model%latitude
       longitude = init_model%longitude
    ENDIF

    f  = 2.0_wp * omega * SIN( latitude / 180.0_wp * pi )
    fs = 2.0_wp * omega * COS( latitude / 180.0_wp * pi )

!
!-- Check and set buoyancy related parameters and switches
    IF ( reference_state == 'horizontal_average' )  THEN
       CONTINUE
    ELSEIF ( reference_state == 'initial_profile' )  THEN
       use_initial_profile_as_reference = .TRUE.
    ELSEIF ( reference_state == 'single_value' )  THEN
       use_single_reference_value = .TRUE.
       IF ( pt_reference == 9999999.9_wp )  pt_reference = pt_surface
       vpt_reference = pt_reference * ( 1.0_wp + 0.61_wp * q_surface )
    ELSE
       message_string = 'illegal value for reference_state: "' //              &
                        TRIM( reference_state ) // '"'
       CALL message( 'check_parameters', 'PA0056', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Sign of buoyancy/stability terms
    IF ( ocean )  atmos_ocean_sign = -1.0_wp

!
!-- Ocean version must use flux boundary conditions at the top
    IF ( ocean .AND. .NOT. use_top_fluxes )  THEN
       message_string = 'use_top_fluxes must be .TRUE. in ocean mode'
       CALL message( 'check_parameters', 'PA0042', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- In case of a given slope, compute the relevant quantities
    IF ( alpha_surface /= 0.0_wp )  THEN
       IF ( ABS( alpha_surface ) > 90.0_wp )  THEN
          WRITE( message_string, * ) 'ABS( alpha_surface = ', alpha_surface,   &
                                     ' ) must be < 90.0'
          CALL message( 'check_parameters', 'PA0043', 1, 2, 0, 6, 0 )
       ENDIF
       sloping_surface = .TRUE.
       cos_alpha_surface = COS( alpha_surface / 180.0_wp * pi )
       sin_alpha_surface = SIN( alpha_surface / 180.0_wp * pi )
    ENDIF

!
!-- Check time step and cfl_factor
    IF ( dt /= -1.0_wp )  THEN
       IF ( dt <= 0.0_wp )  THEN
          WRITE( message_string, * ) 'dt = ', dt , ' <= 0.0'
          CALL message( 'check_parameters', 'PA0044', 1, 2, 0, 6, 0 )
       ENDIF
       dt_3d = dt
       dt_fixed = .TRUE.
    ENDIF

    IF ( cfl_factor <= 0.0_wp  .OR.  cfl_factor > 1.0_wp )  THEN
       IF ( cfl_factor == -1.0_wp )  THEN
          IF ( timestep_scheme == 'runge-kutta-2' )  THEN
             cfl_factor = 0.8_wp
          ELSEIF ( timestep_scheme == 'runge-kutta-3' )  THEN
             cfl_factor = 0.9_wp
          ELSE
             cfl_factor = 0.9_wp
          ENDIF
       ELSE
          WRITE( message_string, * ) 'cfl_factor = ', cfl_factor,              &
                 ' out of range &0.0 < cfl_factor <= 1.0 is required'
          CALL message( 'check_parameters', 'PA0045', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Store simulated time at begin
    simulated_time_at_begin = simulated_time

!
!-- Store reference time for coupled runs and change the coupling flag,
!-- if ...
    IF ( simulated_time == 0.0_wp )  THEN
       IF ( coupling_start_time == 0.0_wp )  THEN
          time_since_reference_point = 0.0_wp
       ELSEIF ( time_since_reference_point < 0.0_wp )  THEN
          run_coupled = .FALSE.
       ENDIF
    ENDIF

!
!-- Set wind speed in the Galilei-transformed system
    IF ( galilei_transformation )  THEN
       IF ( use_ug_for_galilei_tr                    .AND.                     &
            ug_vertical_gradient_level(1) == 0.0_wp  .AND.                     &
            ug_vertical_gradient(1) == 0.0_wp        .AND.                     &
            vg_vertical_gradient_level(1) == 0.0_wp  .AND.                     &
            vg_vertical_gradient(1) == 0.0_wp )  THEN
          u_gtrans = ug_surface * 0.6_wp
          v_gtrans = vg_surface * 0.6_wp
       ELSEIF ( use_ug_for_galilei_tr  .AND.                                   &
                ( ug_vertical_gradient_level(1) /= 0.0_wp  .OR.                &
                ug_vertical_gradient(1) /= 0.0_wp ) )  THEN
          message_string = 'baroclinity (ug) not allowed simultaneously' //    &
                           ' with galilei transformation'
          CALL message( 'check_parameters', 'PA0046', 1, 2, 0, 6, 0 )
       ELSEIF ( use_ug_for_galilei_tr  .AND.                                   &
                ( vg_vertical_gradient_level(1) /= 0.0_wp  .OR.                &
                vg_vertical_gradient(1) /= 0.0_wp ) )  THEN
          message_string = 'baroclinity (vg) not allowed simultaneously' //    &
                           ' with galilei transformation'
          CALL message( 'check_parameters', 'PA0047', 1, 2, 0, 6, 0 )
       ELSE
          message_string = 'variable translation speed used for Galilei-' //   &
             'transformation, which may cause & instabilities in stably ' //   &
             'stratified regions'
          CALL message( 'check_parameters', 'PA0048', 0, 1, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- In case of using a prandtl-layer, calculated (or prescribed) surface
!-- fluxes have to be used in the diffusion-terms
    IF ( constant_flux_layer )  use_surface_fluxes = .TRUE.

!-- Check initial conditions for bubble case
    IF ( INDEX( initializing_actions, 'initialize_3D_bubble' ) /= 0 .OR.         &
         INDEX( initializing_actions, 'initialize_2D_bubble' ) /= 0 ) THEN

       IF ( bubble_radius == 9999999.9_wp ) THEN
          message_string = 'initializing_actions includes bubble, so '// &
                           'bubble_radius must be specified in namelist'
          CALL message( message_string, 'PA0562', 1, 2, 0, 6, 0 )
       ELSEIF ( bubble_radius == 0 ) THEN
          message_string = 'bubble_radius is zero, so bubble will not be '//     &
                           'initialized'
          CALL message( message_string, 'PA0562', 0, 2, 0, 6, 0 )
       ELSEIF ( bubble_radius < 0 ) THEN
          message_string = 'bubble_radius is negative. please specify a non-'//  &
                           'negative bubble_radius'
          CALL message( message_string, 'PA0562', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( bubble_center_x /= 9999999.9_wp .AND.                              &
            ( ( bubble_center_x + bubble_radius ) > nx*dx .OR.                 &
              ( bubble_center_x - bubble_radius ) < 0          ) ) THEN
          message_string = 'bubble is outside the x-axis limits'
          CALL message( message_string, 'PA0562', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( bubble_center_y /= 9999999.9_wp .AND.                              &
            ( ( bubble_center_y + bubble_radius ) > ny*dy .OR.                 &
              ( bubble_center_y - bubble_radius ) < 0          ) ) THEN
          message_string = 'bubble is outside the y-axis limits'
          CALL message( message_string, 'PA0562', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( bubble_center_z /= 9999999.9_wp .AND.                              &
            ( ( bubble_center_z + bubble_radius ) > zu(nzt) .OR.               &
              ( bubble_center_z - bubble_radius ) < zu(nzb) ) ) THEN
          message_string = 'bubble is outside the z-axis limits'
          CALL message( message_string, 'PA0562', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( bubble_pt == 9999999.9_wp .AND. bubble_sa == 9999999.9_wp ) THEN
          message_string = 'initializing_actions includes bubble, so either '//   &
                           'bubble_pt or bubble_sa should be specified.'
          CALL message( message_string, 'PA0562', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( .NOT. ocean .AND. bubble_pt == 0.0_wp ) THEN
          message_string = 'bubble_pt is 0 so bubble will not be initialized.'
          CALL message( message_string, 'PA0562', 0, 2, 0, 6, 0 )
       ENDIF

       IF ( ocean .AND. bubble_pt == 0.0_wp .AND. bubble_sa == 0.0_wp ) THEN
          message_string = 'Both bubble_pt and bubble_sa are 0 so bubble will ' // &
                           'not be initialized.'
          CALL message( message_string, 'PA0562', 0, 2, 0, 6, 0 )
       ENDIF

    ENDIF

!
!-- Check boundary conditions and set internal variables:
!-- Attention: the lateral boundary conditions have been already checked in
!-- parin
!
!-- Non-cyclic lateral boundaries require the multigrid method and Piascek-
!-- Willimas or Wicker - Skamarock advection scheme. Several schemes
!-- and tools do not work with non-cyclic boundary conditions.
    IF ( bc_lr /= 'cyclic'  .OR.  bc_ns /= 'cyclic' )  THEN
       IF ( psolver(1:9) /= 'multigrid' )  THEN
          message_string = 'non-cyclic lateral boundaries do not allow ' //    &
                           'psolver = "' // TRIM( psolver ) // '"'
          CALL message( 'check_parameters', 'PA0051', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( momentum_advec /= 'pw-scheme'  .AND.                               &
            momentum_advec /= 'ws-scheme' )  THEN

          message_string = 'non-cyclic lateral boundaries do not allow ' //    &
                           'momentum_advec = "' // TRIM( momentum_advec ) // '"'
          CALL message( 'check_parameters', 'PA0052', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( scalar_advec /= 'pw-scheme'  .AND.                                 &
            scalar_advec /= 'ws-scheme' )  THEN
          message_string = 'non-cyclic lateral boundaries do not allow ' //    &
                           'scalar_advec = "' // TRIM( scalar_advec ) // '"'
          CALL message( 'check_parameters', 'PA0053', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( galilei_transformation )  THEN
          message_string = 'non-cyclic lateral boundaries do not allow ' //    &
                           'galilei_transformation = .T.'
          CALL message( 'check_parameters', 'PA0054', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Bottom boundary condition for the turbulent Kinetic energy
    IF ( bc_e_b == 'neumann' )  THEN
       ibc_e_b = 1
    ELSEIF ( bc_e_b == '(u*)**2+neumann' )  THEN
       ibc_e_b = 2
       IF ( .NOT. constant_flux_layer )  THEN
          bc_e_b = 'neumann'
          ibc_e_b = 1
          message_string = 'boundary condition bc_e_b changed to "' //         &
                           TRIM( bc_e_b ) // '"'
          CALL message( 'check_parameters', 'PA0057', 0, 1, 0, 6, 0 )
       ENDIF
    ELSE
       message_string = 'unknown boundary condition: bc_e_b = "' //            &
                        TRIM( bc_e_b ) // '"'
       CALL message( 'check_parameters', 'PA0058', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Boundary conditions for perturbation pressure
    IF ( bc_p_b == 'dirichlet' )  THEN
       ibc_p_b = 0
    ELSEIF ( bc_p_b == 'neumann' )  THEN
       ibc_p_b = 1
    ELSE
       message_string = 'unknown boundary condition: bc_p_b = "' //            &
                        TRIM( bc_p_b ) // '"'
       CALL message( 'check_parameters', 'PA0059', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( bc_p_t == 'dirichlet' )  THEN
       ibc_p_t = 0
!-- TO_DO: later set bc_p_t to neumann before, in case of nested domain
    ELSEIF ( bc_p_t == 'neumann' .OR. bc_p_t == 'nested'  .OR.                 &
             bc_p_t == 'forcing'  )  THEN
       ibc_p_t = 1
    ELSE
       message_string = 'unknown boundary condition: bc_p_t = "' //            &
                        TRIM( bc_p_t ) // '"'
       CALL message( 'check_parameters', 'PA0061', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Boundary conditions for potential temperature
    IF ( coupling_mode == 'atmosphere_to_ocean' )  THEN
       ibc_pt_b = 2
    ELSE
       IF ( bc_pt_b == 'dirichlet' )  THEN
          ibc_pt_b = 0
       ELSEIF ( bc_pt_b == 'neumann' )  THEN
          ibc_pt_b = 1
       ELSE
          message_string = 'unknown boundary condition: bc_pt_b = "' //        &
                           TRIM( bc_pt_b ) // '"'
          CALL message( 'check_parameters', 'PA0062', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( bc_pt_t == 'dirichlet' )  THEN
       ibc_pt_t = 0
    ELSEIF ( bc_pt_t == 'neumann' )  THEN
       ibc_pt_t = 1
    ELSEIF ( bc_pt_t == 'initial_gradient' )  THEN
       ibc_pt_t = 2
    ELSEIF ( bc_pt_t == 'nested'  .OR.  bc_pt_t == 'forcing' )  THEN
       ibc_pt_t = 3
    ELSE
       message_string = 'unknown boundary condition: bc_pt_t = "' //           &
                        TRIM( bc_pt_t ) // '"'
       CALL message( 'check_parameters', 'PA0063', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( ANY( wall_heatflux /= 0.0_wp )  .AND.                        &
         surface_heatflux == 9999999.9_wp )  THEN
       message_string = 'wall_heatflux additionally requires ' //     &
                        'setting of surface_heatflux'
       CALL message( 'check_parameters', 'PA0443', 1, 2, 0, 6, 0 )
    ENDIF

!
!   This IF clause needs revision, got too complex!!
    IF ( surface_heatflux == 9999999.9_wp  )  THEN
       constant_heatflux = .FALSE.
       IF ( large_scale_forcing  .OR.  land_surface  .OR.  urban_surface )  THEN
          IF ( ibc_pt_b == 0 )  THEN
             constant_heatflux = .FALSE.
          ELSEIF ( ibc_pt_b == 1 )  THEN
             constant_heatflux = .TRUE.
             surface_heatflux = 0.0_wp
          ENDIF
       ENDIF
    ELSE
       constant_heatflux = .TRUE.
    ENDIF

    IF ( top_heatflux     == 9999999.9_wp )  constant_top_heatflux = .FALSE.

    IF ( neutral )  THEN

       IF ( surface_heatflux /= 0.0_wp  .AND.                                  &
            surface_heatflux /= 9999999.9_wp )  THEN
          message_string = 'heatflux must not be set for pure neutral flow'
          CALL message( 'check_parameters', 'PA0351', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( top_heatflux /= 0.0_wp  .AND.  top_heatflux /= 9999999.9_wp )      &
       THEN
          message_string = 'heatflux must not be set for pure neutral flow'
          CALL message( 'check_parameters', 'PA0351', 1, 2, 0, 6, 0 )
       ENDIF

    ENDIF

    IF ( top_momentumflux_u /= 9999999.9_wp  .AND.                             &
         top_momentumflux_v /= 9999999.9_wp )  THEN
       constant_top_momentumflux = .TRUE.
    ELSEIF (  .NOT. ( top_momentumflux_u == 9999999.9_wp  .AND.                &
           top_momentumflux_v == 9999999.9_wp ) )  THEN
       message_string = 'both, top_momentumflux_u AND top_momentumflux_v ' //  &
                        'must be set'
       CALL message( 'check_parameters', 'PA0064', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- A given surface temperature implies Dirichlet boundary condition for
!-- temperature. In this case specification of a constant heat flux is
!-- forbidden.
    IF ( ibc_pt_b == 0  .AND.  constant_heatflux  .AND.                        &
         surface_heatflux /= 0.0_wp )  THEN
       message_string = 'boundary_condition: bc_pt_b = "' // TRIM( bc_pt_b ) //&
                        '& is not allowed with constant_heatflux = .TRUE.'
       CALL message( 'check_parameters', 'PA0065', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( constant_heatflux  .AND.  pt_surface_initial_change /= 0.0_wp )  THEN
       WRITE ( message_string, * )  'constant_heatflux = .TRUE. is not allo',  &
               'wed with pt_surface_initial_change (/=0) = ',                  &
               pt_surface_initial_change
       CALL message( 'check_parameters', 'PA0066', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- A given temperature at the top implies Dirichlet boundary condition for
!-- temperature. In this case specification of a constant heat flux is
!-- forbidden.
    IF ( ibc_pt_t == 0  .AND.  constant_top_heatflux  .AND.                    &
         top_heatflux /= 0.0_wp )  THEN
       message_string = 'boundary_condition: bc_pt_t = "' // TRIM( bc_pt_t ) //&
                        '" is not allowed with constant_top_heatflux = .TRUE.'
       CALL message( 'check_parameters', 'PA0067', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Boundary conditions for salinity
    IF ( ocean )  THEN
       IF ( bc_sa_t == 'dirichlet' )  THEN
          ibc_sa_t = 0
       ELSEIF ( bc_sa_t == 'neumann' )  THEN
          ibc_sa_t = 1
       ELSE
          message_string = 'unknown boundary condition: bc_sa_t = "' //        &
                           TRIM( bc_sa_t ) // '"'
          CALL message( 'check_parameters', 'PA0068', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( bc_sa_b == 'dirichlet' )  THEN
          ibc_sa_b = 0
          CALL location_message('ib_sa_b assigned for dirichlet conditions',.TRUE.) !CB
       ELSEIF ( bc_sa_b == 'neumann' )  THEN
          ibc_sa_b = 1
          CALL location_message('ib_sa_b assigned for neumann conditions',.TRUE.) !CB
       ELSE
          message_string = 'unknown boundary condition: bc_sa_b = "' //        &
                           TRIM( bc_sa_b ) // '"'
          CALL message( 'check_parameters', 'PA0068', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( ibc_sa_t == 1  .AND.  top_salinityflux == 9999999.9_wp )  THEN
          message_string = 'boundary condition: bc_sa_t = "' //                &
                           TRIM( bc_sa_t ) // '" requires to set ' //          &
                           'top_salinityflux'
          CALL message( 'check_parameters', 'PA0069', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( ibc_sa_b == 1  .AND.  bottom_salinityflux == 9999999.9_wp )  THEN
          message_string = 'boundary condition: bc_sa_b = "' //                &
                           TRIM( bc_sa_b ) // '" requires to set ' //          &
                           'bottom_salinityflux'
          CALL message( 'check_parameters', 'PA0069', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    A fixed salinity at the top implies Dirichlet boundary condition for
!--    salinity. In this case specification of a constant salinity flux is
!--    forbidden.
       IF ( top_salinityflux == 9999999.9_wp )  constant_top_salinityflux = .FALSE.
       IF ( ibc_sa_t == 0  .AND.  constant_top_salinityflux  .AND.             &
            top_salinityflux /= 0.0_wp )  THEN
          message_string = 'boundary condition: bc_sa_t = "' //                &
                           TRIM( bc_sa_t ) // '" is not allowed with ' //      &
                           'top_salinityflux /= 0.0'
          CALL message( 'check_parameters', 'PA0070', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( bottom_salinityflux == 9999999.9_wp )  constant_bottom_salinityflux = .FALSE.
       IF ( ibc_sa_b == 0  .AND.  constant_bottom_salinityflux  .AND.             &
            bottom_salinityflux /= 0.0_wp )  THEN
          message_string = 'boundary condition: bc_sa_b = "' //                &
                           TRIM( bc_sa_t ) // '" is not allowed with ' //      &
                           'bottom_salinityflux /= 0.0'
          CALL message( 'check_parameters', 'PA0070', 1, 2, 0, 6, 0 )
       ENDIF

    ENDIF

!
!-- Set boundary conditions for total water content
    IF ( humidity )  THEN

       IF ( ANY( wall_humidityflux /= 0.0_wp )  .AND.                        &
            surface_waterflux == 9999999.9_wp )  THEN
          message_string = 'wall_humidityflux additionally requires ' //     &
                           'setting of surface_waterflux'
          CALL message( 'check_parameters', 'PA0444', 1, 2, 0, 6, 0 )
       ENDIF

       CALL set_bc_scalars( 'q', bc_q_b, bc_q_t, ibc_q_b, ibc_q_t,           &
                            'PA0071', 'PA0072' )

       IF ( surface_waterflux == 9999999.9_wp  )  THEN
          constant_waterflux = .FALSE.
          IF ( large_scale_forcing .OR. land_surface )  THEN
             IF ( ibc_q_b == 0 )  THEN
                constant_waterflux = .FALSE.
             ELSEIF ( ibc_q_b == 1 )  THEN
                constant_waterflux = .TRUE.
             ENDIF
          ENDIF
       ELSE
          constant_waterflux = .TRUE.
       ENDIF

       CALL check_bc_scalars( 'q', bc_q_b, ibc_q_b, 'PA0073', 'PA0074',        &
                              constant_waterflux, q_surface_initial_change )

    ENDIF

    IF ( passive_scalar )  THEN

       IF ( ANY( wall_scalarflux /= 0.0_wp )  .AND.                            &
            surface_scalarflux == 9999999.9_wp )  THEN
          message_string = 'wall_scalarflux additionally requires ' //         &
                           'setting of surface_scalarflux'
          CALL message( 'check_parameters', 'PA0445', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( surface_scalarflux == 9999999.9_wp )  constant_scalarflux = .FALSE.

       CALL set_bc_scalars( 's', bc_s_b, bc_s_t, ibc_s_b, ibc_s_t,             &
                            'PA0071', 'PA0072' )

       CALL check_bc_scalars( 's', bc_s_b, ibc_s_b, 'PA0073', 'PA0074',        &
                              constant_scalarflux, s_surface_initial_change )

       IF ( top_scalarflux == 9999999.9_wp )  constant_top_scalarflux = .FALSE.
!
!--    A fixed scalar concentration at the top implies Dirichlet boundary
!--    condition for scalar. Hence, in this case specification of a constant
!--    scalar flux is forbidden.
       IF ( ( ibc_s_t == 0 .OR. ibc_s_t == 2 )  .AND.  constant_top_scalarflux &
               .AND.  top_scalarflux /= 0.0_wp )  THEN
          message_string = 'boundary condition: bc_s_t = "' //                 &
                           TRIM( bc_sa_t ) // '" is not allowed with ' //      &
                           'top_scalarflux /= 0.0'
          CALL message( 'check_parameters', 'PA0441', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Boundary conditions for chemical species
    IF ( air_chemistry )  CALL chem_boundary_conds( 'init' )

!
!-- Boundary conditions for horizontal components of wind speed
    IF ( bc_uv_b == 'dirichlet' )  THEN
       ibc_uv_b = 0
    ELSEIF ( bc_uv_b == 'neumann' )  THEN
       ibc_uv_b = 1
       IF ( constant_flux_layer )  THEN
          message_string = 'boundary condition: bc_uv_b = "' //                &
               TRIM( bc_uv_b ) // '" is not allowed with constant_flux_layer'  &
               // ' = .TRUE.'
          CALL message( 'check_parameters', 'PA0075', 1, 2, 0, 6, 0 )
       ENDIF
    ELSE
       message_string = 'unknown boundary condition: bc_uv_b = "' //           &
                        TRIM( bc_uv_b ) // '"'
       CALL message( 'check_parameters', 'PA0076', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- In case of coupled simulations u and v at the ground in atmosphere will be
!-- assigned with the u and v values of the ocean surface
    IF ( coupling_mode == 'atmosphere_to_ocean' )  THEN
       ibc_uv_b = 2
    ENDIF

    IF ( coupling_mode == 'ocean_to_atmosphere' )  THEN
       bc_uv_t = 'neumann'
       ibc_uv_t = 1
    ELSE
       IF ( bc_uv_t == 'dirichlet' .OR. bc_uv_t == 'dirichlet_0' )  THEN
          ibc_uv_t = 0
          IF ( bc_uv_t == 'dirichlet_0' )  THEN
!
!--          Velocities for the initial u,v-profiles are set zero at the top
!--          in case of dirichlet_0 conditions
             u_init(nzt+1)    = 0.0_wp
             v_init(nzt+1)    = 0.0_wp
          ENDIF
       ELSEIF ( bc_uv_t == 'neumann' )  THEN
          ibc_uv_t = 1
       ELSEIF ( bc_uv_t == 'nested'  .OR.  bc_uv_t == 'forcing' )  THEN
          ibc_uv_t = 3
       ELSE
          message_string = 'unknown boundary condition: bc_uv_t = "' //        &
                           TRIM( bc_uv_t ) // '"'
          CALL message( 'check_parameters', 'PA0077', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Compute and check, respectively, the Rayleigh Damping parameter
    IF ( rayleigh_damping_factor == -1.0_wp )  THEN
       rayleigh_damping_factor = 0.0_wp
    ELSE
       IF ( rayleigh_damping_factor < 0.0_wp  .OR.                             &
            rayleigh_damping_factor > 1.0_wp )  THEN
          WRITE( message_string, * )  'rayleigh_damping_factor = ',            &
                              rayleigh_damping_factor, ' out of range [0.0,1.0]'
          CALL message( 'check_parameters', 'PA0078', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( rayleigh_damping_height == -1.0_wp )  THEN
       IF (  .NOT.  ocean )  THEN
          rayleigh_damping_height = 0.66666666666_wp * zu(nzt)
       ELSE
          rayleigh_damping_height = 0.66666666666_wp * zu(nzb)
       ENDIF
    ELSE
       IF (  .NOT.  ocean )  THEN
          IF ( rayleigh_damping_height < 0.0_wp  .OR.                          &
               rayleigh_damping_height > zu(nzt) )  THEN
             WRITE( message_string, * )  'rayleigh_damping_height = ',         &
                   rayleigh_damping_height, ' out of range [0.0,', zu(nzt), ']'
             CALL message( 'check_parameters', 'PA0079', 1, 2, 0, 6, 0 )
          ENDIF
       ELSE
          IF ( rayleigh_damping_height > 0.0_wp  .OR.                          &
               rayleigh_damping_height < zu(nzb) )  THEN
             WRITE( message_string, * )  'rayleigh_damping_height = ',         &
                   rayleigh_damping_height, ' out of range [0.0,', zu(nzb), ']'
             CALL message( 'check_parameters', 'PA0079', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
    ENDIF

!
!-- Check number of chosen statistic regions
    IF ( statistic_regions < 0 )  THEN
       WRITE ( message_string, * ) 'number of statistic_regions = ',           &
                   statistic_regions+1, ' is not allowed'
       CALL message( 'check_parameters', 'PA0082', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( normalizing_region > statistic_regions  .OR.                          &
         normalizing_region < 0)  THEN
       WRITE ( message_string, * ) 'normalizing_region = ',                    &
                normalizing_region, ' must be >= 0 and <= ',statistic_regions, &
                ' (value of statistic_regions)'
       CALL message( 'check_parameters', 'PA0083', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set the default intervals for data output, if necessary
!-- NOTE: dt_dosp has already been set in spectra_parin
    IF ( dt_data_output /= 9999999.9_wp )  THEN
       IF ( dt_dopr           == 9999999.9_wp )  dt_dopr           = dt_data_output
       IF ( dt_dopts          == 9999999.9_wp )  dt_dopts          = dt_data_output
       IF ( dt_do2d_xy        == 9999999.9_wp )  dt_do2d_xy        = dt_data_output
       IF ( dt_do2d_xz        == 9999999.9_wp )  dt_do2d_xz        = dt_data_output
       IF ( dt_do2d_yz        == 9999999.9_wp )  dt_do2d_yz        = dt_data_output
       IF ( dt_do3d           == 9999999.9_wp )  dt_do3d           = dt_data_output
       IF ( dt_data_output_av == 9999999.9_wp )  dt_data_output_av = dt_data_output
       DO  mid = 1, max_masks
          IF ( dt_domask(mid) == 9999999.9_wp )  dt_domask(mid)    = dt_data_output
       ENDDO
    ENDIF

!
!-- Set the default skip time intervals for data output, if necessary
    IF ( skip_time_dopr    == 9999999.9_wp )                                   &
                                       skip_time_dopr    = skip_time_data_output
    IF ( skip_time_do2d_xy == 9999999.9_wp )                                   &
                                       skip_time_do2d_xy = skip_time_data_output
    IF ( skip_time_do2d_xz == 9999999.9_wp )                                   &
                                       skip_time_do2d_xz = skip_time_data_output
    IF ( skip_time_do2d_yz == 9999999.9_wp )                                   &
                                       skip_time_do2d_yz = skip_time_data_output
    IF ( skip_time_do3d    == 9999999.9_wp )                                   &
                                       skip_time_do3d    = skip_time_data_output
    IF ( skip_time_data_output_av == 9999999.9_wp )                            &
                                skip_time_data_output_av = skip_time_data_output
    DO  mid = 1, max_masks
       IF ( skip_time_domask(mid) == 9999999.9_wp )                            &
                                skip_time_domask(mid)    = skip_time_data_output
    ENDDO

!
!-- Check the average intervals (first for 3d-data, then for profiles)
    IF ( averaging_interval > dt_data_output_av )  THEN
       WRITE( message_string, * )  'averaging_interval = ',                    &
             averaging_interval, ' must be <= dt_data_output_av = ',           &
             dt_data_output_av
       CALL message( 'check_parameters', 'PA0085', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( averaging_interval_pr == 9999999.9_wp )  THEN
       averaging_interval_pr = averaging_interval
    ENDIF

    IF ( averaging_interval_pr > dt_dopr )  THEN
       WRITE( message_string, * )  'averaging_interval_pr = ',                 &
             averaging_interval_pr, ' must be <= dt_dopr = ', dt_dopr
       CALL message( 'check_parameters', 'PA0086', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set the default interval for profiles entering the temporal average
    IF ( dt_averaging_input_pr == 9999999.9_wp )  THEN
       dt_averaging_input_pr = dt_averaging_input
    ENDIF

!
!-- Set the default interval for the output of timeseries to a reasonable
!-- value (tries to minimize the number of calls of flow_statistics)
    IF ( dt_dots == 9999999.9_wp )  THEN
       IF ( averaging_interval_pr == 0.0_wp )  THEN
          dt_dots = MIN( dt_run_control, dt_dopr )
       ELSE
          dt_dots = MIN( dt_run_control, dt_averaging_input_pr )
       ENDIF
    ENDIF

!
!-- Check the sample rate for averaging (first for 3d-data, then for profiles)
    IF ( dt_averaging_input > averaging_interval )  THEN
       WRITE( message_string, * )  'dt_averaging_input = ',                    &
                dt_averaging_input, ' must be <= averaging_interval = ',       &
                averaging_interval
       CALL message( 'check_parameters', 'PA0088', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( dt_averaging_input_pr > averaging_interval_pr )  THEN
       WRITE( message_string, * )  'dt_averaging_input_pr = ',                 &
                dt_averaging_input_pr, ' must be <= averaging_interval_pr = ', &
                averaging_interval_pr
       CALL message( 'check_parameters', 'PA0089', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set the default value for the integration interval of precipitation amount
    IF ( microphysics_seifert  .OR.  microphysics_kessler )  THEN
       IF ( precipitation_amount_interval == 9999999.9_wp )  THEN
          precipitation_amount_interval = dt_do2d_xy
       ELSE
          IF ( precipitation_amount_interval > dt_do2d_xy )  THEN
             WRITE( message_string, * )  'precipitation_amount_interval = ',   &
                 precipitation_amount_interval, ' must not be larger than ',   &
                 'dt_do2d_xy = ', dt_do2d_xy
             CALL message( 'check_parameters', 'PA0090', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
    ENDIF

!
!-- Determine the number of output profiles and check whether they are
!-- permissible
    DO  WHILE ( data_output_pr(dopr_n+1) /= '          ' )

       dopr_n = dopr_n + 1
       i = dopr_n

!
!--    Determine internal profile number (for hom, homs)
!--    and store height levels
       SELECT CASE ( TRIM( data_output_pr(i) ) )

          CASE ( 'u', '#u' )
             dopr_index(i) = 1
             dopr_unit(i)  = 'm/s'
             hom(:,2,1,:)  = SPREAD( zu, 2, statistic_regions+1 )
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 5
                hom(:,2,5,:)          = SPREAD( zu, 2, statistic_regions+1 )
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'v', '#v' )
             dopr_index(i) = 2
             dopr_unit(i)  = 'm/s'
             hom(:,2,2,:)  = SPREAD( zu, 2, statistic_regions+1 )
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 6
                hom(:,2,6,:)          = SPREAD( zu, 2, statistic_regions+1 )
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'w' )
             dopr_index(i) = 3
             dopr_unit(i)  = 'm/s'
             hom(:,2,3,:)  = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'pt', '#pt' )
             IF ( .NOT. cloud_physics ) THEN
                dopr_index(i) = 4
                dopr_unit(i)  = 'K'
                hom(:,2,4,:)  = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 7
                   hom(:,2,7,:)          = SPREAD( zu, 2, statistic_regions+1 )
                   hom(nzb,2,7,:)        = 0.0_wp    ! because zu(nzb) is negative
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ELSE
                dopr_index(i) = 43
                dopr_unit(i)  = 'K'
                hom(:,2,43,:)  = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 28
                   hom(:,2,28,:)         = SPREAD( zu, 2, statistic_regions+1 )
                   hom(nzb,2,28,:)       = 0.0_wp    ! because zu(nzb) is negative
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ENDIF

          CASE ( 'e', '#e' )
             dopr_index(i)  = 8
             dopr_unit(i)   = 'm2/s2'
             hom(:,2,8,:)   = SPREAD( zu, 2, statistic_regions+1 )
             hom(nzb,2,8,:) = 0.0_wp
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 8
                hom(:,2,8,:)          = SPREAD( zu, 2, statistic_regions+1 )
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'km', '#km' )
             dopr_index(i)  = 9
             dopr_unit(i)   = 'm2/s'
             hom(:,2,9,:)   = SPREAD( zu, 2, statistic_regions+1 )
             hom(nzb,2,9,:) = 0.0_wp
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 23
                hom(:,2,23,:)         = hom(:,2,9,:)
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'kh', '#kh' )
             dopr_index(i)   = 10
             dopr_unit(i)    = 'm2/s'
             hom(:,2,10,:)   = SPREAD( zu, 2, statistic_regions+1 )
             hom(nzb,2,10,:) = 0.0_wp
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 24
                hom(:,2,24,:)         = hom(:,2,10,:)
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'l', '#l' )
             dopr_index(i)   = 11
             dopr_unit(i)    = 'm'
             hom(:,2,11,:)   = SPREAD( zu, 2, statistic_regions+1 )
             hom(nzb,2,11,:) = 0.0_wp
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 25
                hom(:,2,25,:)         = hom(:,2,11,:)
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'w"u"' )
             dopr_index(i) = 12
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,12,:) = SPREAD( zw, 2, statistic_regions+1 )
             IF ( constant_flux_layer )  hom(nzb,2,12,:) = zu(1)

          CASE ( 'w*u*' )
             dopr_index(i) = 13
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,13,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w"v"' )
             dopr_index(i) = 14
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,14,:) = SPREAD( zw, 2, statistic_regions+1 )
             IF ( constant_flux_layer )  hom(nzb,2,14,:) = zu(1)

          CASE ( 'w*v*' )
             dopr_index(i) = 15
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,15,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w"pt"' )
             dopr_index(i) = 16
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,16,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w*pt*' )
             dopr_index(i) = 17
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,17,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'wpt' )
             dopr_index(i) = 18
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,18,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'wu' )
             dopr_index(i) = 19
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,19,:) = SPREAD( zw, 2, statistic_regions+1 )
             IF ( constant_flux_layer )  hom(nzb,2,19,:) = zu(1)

          CASE ( 'wv' )
             dopr_index(i) = 20
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,20,:) = SPREAD( zw, 2, statistic_regions+1 )
             IF ( constant_flux_layer )  hom(nzb,2,20,:) = zu(1)

          CASE ( 'w*pt*BC' )
             dopr_index(i) = 21
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,21,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'wptBC' )
             dopr_index(i) = 22
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,22,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'sa', '#sa' )
             IF ( .NOT. ocean )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for ocean = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 23
                dopr_unit(i)  = 'psu'
                hom(:,2,23,:) = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 26
                   hom(:,2,26,:)         = SPREAD( zu, 2, statistic_regions+1 )
                   hom(nzb,2,26,:)       = 0.0_wp    ! because zu(nzb) is negative
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ENDIF

          CASE ( 'u*2' )
             dopr_index(i) = 30
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,30,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'v*2' )
             dopr_index(i) = 31
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,31,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'w*2' )
             dopr_index(i) = 32
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,32,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'pt*2' )
             dopr_index(i) = 33
             dopr_unit(i)  = 'K2'
             hom(:,2,33,:) = SPREAD( zu, 2, statistic_regions+1 )


          CASE ( 'sa*2' )
             dopr_index(i) = 153
             dopr_unit(i)  = 'PSU^2'
             hom(:,2,153,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'e*' )
             dopr_index(i) = 34
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,34,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'w*2pt*' )
             dopr_index(i) = 35
             dopr_unit(i)  = 'K m2/s2'
             hom(:,2,35,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w*pt*2' )
             dopr_index(i) = 36
             dopr_unit(i)  = 'K2 m/s'
             hom(:,2,36,:) = SPREAD( zw, 2, statistic_regions+1 )
          CASE ( 'w*2sa*' )
             dopr_index(i) = 154
             dopr_unit(i)  = 'PSU m2/s2'
             hom(:,2,154,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w*sa*2' )
             dopr_index(i) = 155
             dopr_unit(i)  = 'sa2 m/s'
             hom(:,2,155,:) = SPREAD( zw, 2, statistic_regions+1 )


          CASE ( 'w*e*' )
             dopr_index(i) = 37
             dopr_unit(i)  = 'm3/s3'
             hom(:,2,37,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w*3' )
             dopr_index(i) = 38
             dopr_unit(i)  = 'm3/s3'
             hom(:,2,38,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'Sw' )
             dopr_index(i) = 39
             dopr_unit(i)  = 'none'
             hom(:,2,39,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'p' )
             dopr_index(i) = 40
             dopr_unit(i)  = 'Pa'
             hom(:,2,40,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'q', '#q' )
             IF ( .NOT. humidity )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for humidity = .FALSE.'
                CALL message( 'check_parameters', 'PA0092', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 41
                dopr_unit(i)  = 'kg/kg'
                hom(:,2,41,:) = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 26
                   hom(:,2,26,:)         = SPREAD( zu, 2, statistic_regions+1 )
                   hom(nzb,2,26,:)       = 0.0_wp    ! because zu(nzb) is negative
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ENDIF

          CASE ( 's', '#s' )
             IF ( .NOT. passive_scalar )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for passive_scalar = .FALSE.'
                CALL message( 'check_parameters', 'PA0093', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 115
                dopr_unit(i)  = 'kg/m3'
                hom(:,2,115,:) = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 121
                   hom(:,2,121,:)        = SPREAD( zu, 2, statistic_regions+1 )
                   hom(nzb,2,121,:)      = 0.0_wp    ! because zu(nzb) is negative
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ENDIF

          CASE ( 'qv', '#qv' )
             IF ( .NOT. cloud_physics ) THEN
                dopr_index(i) = 41
                dopr_unit(i)  = 'kg/kg'
                hom(:,2,41,:) = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 26
                   hom(:,2,26,:)         = SPREAD( zu, 2, statistic_regions+1 )
                   hom(nzb,2,26,:)       = 0.0_wp    ! because zu(nzb) is negative
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ELSE
                dopr_index(i) = 42
                dopr_unit(i)  = 'kg/kg'
                hom(:,2,42,:) = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 27
                   hom(:,2,27,:)         = SPREAD( zu, 2, statistic_regions+1 )
                   hom(nzb,2,27,:)       = 0.0_wp   ! because zu(nzb) is negative
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ENDIF

          CASE ( 'lpt', '#lpt' )
             IF ( .NOT. cloud_physics ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for cloud_physics = .FALSE.'
                CALL message( 'check_parameters', 'PA0094', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 4
                dopr_unit(i)  = 'K'
                hom(:,2,4,:)  = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 7
                   hom(:,2,7,:)          = SPREAD( zu, 2, statistic_regions+1 )
                   hom(nzb,2,7,:)        = 0.0_wp    ! because zu(nzb) is negative
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ENDIF

          CASE ( 'vpt', '#vpt' )
             dopr_index(i) = 44
             dopr_unit(i)  = 'K'
             hom(:,2,44,:) = SPREAD( zu, 2, statistic_regions+1 )
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 29
                hom(:,2,29,:)         = SPREAD( zu, 2, statistic_regions+1 )
                hom(nzb,2,29,:)       = 0.0_wp    ! because zu(nzb) is negative
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'w"vpt"' )
             dopr_index(i) = 45
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,45,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w*vpt*' )
             dopr_index(i) = 46
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,46,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'wvpt' )
             dopr_index(i) = 47
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,47,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w"q"' )
             IF (  .NOT.  humidity )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for humidity = .FALSE.'
                CALL message( 'check_parameters', 'PA0092', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 48
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,48,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'w*q*' )
             IF (  .NOT.  humidity )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for humidity = .FALSE.'
                CALL message( 'check_parameters', 'PA0092', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 49
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,49,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'wq' )
             IF (  .NOT.  humidity )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for humidity = .FALSE.'
                CALL message( 'check_parameters', 'PA0092', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 50
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,50,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'w"s"' )
             IF (  .NOT.  passive_scalar )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for passive_scalar = .FALSE.'
                CALL message( 'check_parameters', 'PA0093', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 117
                dopr_unit(i)  = 'kg/m3 m/s'
                hom(:,2,117,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'w*s*' )
             IF (  .NOT.  passive_scalar )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for passive_scalar = .FALSE.'
                CALL message( 'check_parameters', 'PA0093', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 114
                dopr_unit(i)  = 'kg/m3 m/s'
                hom(:,2,114,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'ws' )
             IF (  .NOT.  passive_scalar )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for passive_scalar = .FALSE.'
                CALL message( 'check_parameters', 'PA0093', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 118
                dopr_unit(i)  = 'kg/m3 m/s'
                hom(:,2,118,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'w"qv"' )
             IF ( humidity  .AND.  .NOT.  cloud_physics )  THEN
                dopr_index(i) = 48
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,48,:) = SPREAD( zw, 2, statistic_regions+1 )
             ELSEIF ( humidity  .AND.  cloud_physics )  THEN
                dopr_index(i) = 51
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,51,:) = SPREAD( zw, 2, statistic_regions+1 )
             ELSE
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for cloud_physics = .FALSE. an'  // &
                                 'd humidity = .FALSE.'
                CALL message( 'check_parameters', 'PA0095', 1, 2, 0, 6, 0 )
             ENDIF

          CASE ( 'w*qv*' )
             IF ( humidity  .AND.  .NOT. cloud_physics )                       &
             THEN
                dopr_index(i) = 49
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,49,:) = SPREAD( zw, 2, statistic_regions+1 )
             ELSEIF( humidity .AND. cloud_physics ) THEN
                dopr_index(i) = 52
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,52,:) = SPREAD( zw, 2, statistic_regions+1 )
             ELSE
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for cloud_physics = .FALSE. an'  // &
                                 'd humidity = .FALSE.'
                CALL message( 'check_parameters', 'PA0095', 1, 2, 0, 6, 0 )
             ENDIF

          CASE ( 'wqv' )
             IF ( humidity  .AND.  .NOT.  cloud_physics )  THEN
                dopr_index(i) = 50
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,50,:) = SPREAD( zw, 2, statistic_regions+1 )
             ELSEIF ( humidity  .AND.  cloud_physics )  THEN
                dopr_index(i) = 53
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,53,:) = SPREAD( zw, 2, statistic_regions+1 )
             ELSE
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for cloud_physics = .FALSE. an'  // &
                                 'd humidity = .FALSE.'
                CALL message( 'check_parameters', 'PA0095', 1, 2, 0, 6, 0 )
             ENDIF

          CASE ( 'ql' )
             IF (  .NOT.  cloud_physics  .AND.  .NOT.  cloud_droplets )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for cloud_physics = .FALSE. and' // &
                                 'cloud_droplets = .FALSE.'
                CALL message( 'check_parameters', 'PA0096', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 54
                dopr_unit(i)  = 'kg/kg'
                hom(:,2,54,:)  = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'w*u*u*:dz' )
             dopr_index(i) = 55
             dopr_unit(i)  = 'm2/s3'
             hom(:,2,55,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'w*p*:dz' )
             dopr_index(i) = 56
             dopr_unit(i)  = 'm2/s3'
             hom(:,2,56,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w"e:dz' )
             dopr_index(i) = 57
             dopr_unit(i)  = 'm2/s3'
             hom(:,2,57,:) = SPREAD( zu, 2, statistic_regions+1 )


          CASE ( 'u"pt"' )
             dopr_index(i) = 58
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,58,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'u*pt*' )
             dopr_index(i) = 59
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,59,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'upt_t' )
             dopr_index(i) = 60
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,60,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'v"pt"' )
             dopr_index(i) = 61
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,61,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'v*pt*' )
             dopr_index(i) = 62
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,62,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'vpt_t' )
             dopr_index(i) = 63
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,63,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'rho_ocean' )
             IF (  .NOT.  ocean ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for ocean = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 64
                dopr_unit(i)  = 'kg/m3'
                hom(:,2,64,:) = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 77
                   hom(:,2,77,:)         = SPREAD( zu, 2, statistic_regions+1 )
                   hom(nzb,2,77,:)       = 0.0_wp    ! because zu(nzb) is negative
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ENDIF

           CASE ( 'solar3d' )
             IF (  .NOT.  ocean ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for ocean = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 152
                dopr_unit(i)  = '^oC^{-1}/s'
                hom(:,2,152,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

           CASE ( 'alpha_T' )
             IF (  .NOT.  ocean ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for ocean = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 150
                dopr_unit(i)  = '^oC^{-1}'
                hom(:,2,150,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'beta_S' )
             IF (  .NOT.  ocean ) THEN
                message_string = 'data_output_pr = ' //                        &
                    TRIM( data_output_pr(i) ) // ' is not imp' // &
                    'lemented for ocean = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 151
                dopr_unit(i)  = 'PSU^{-1}'
                hom(:,2,151,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'w"sa"' )
             IF (  .NOT.  ocean ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for ocean = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 65
                dopr_unit(i)  = 'psu m/s'
                hom(:,2,65,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'w*sa*' )
             IF (  .NOT. ocean  ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for ocean = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 66
                dopr_unit(i)  = 'psu m/s'
                hom(:,2,66,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'wsa' )
             IF (  .NOT.  ocean ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for ocean = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 67
                dopr_unit(i)  = 'psu m/s'
                hom(:,2,67,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'u_stk' )
             IF (  .NOT.  ocean ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for ocean = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSEIF (  .NOT.  stokes_force ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for stokes_force = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 161
                dopr_unit(i)  = 'm/s'
                hom(:,2,161,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'v_stk' )
             IF (  .NOT.  ocean ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for ocean = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSEIF (  .NOT.  stokes_force ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for stokes_force = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 162
                dopr_unit(i)  = 'm/s'
                hom(:,2,162,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'u_stk_zw' )
             IF (  .NOT.  ocean ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for ocean = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSEIF (  .NOT.  stokes_force ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for stokes_force = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 163
                dopr_unit(i)  = 'm/s'
                hom(:,2,163,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'v_stk_zw' )
             IF (  .NOT.  ocean ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for ocean = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSEIF (  .NOT.  stokes_force ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for stokes_force = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 164
                dopr_unit(i)  = 'm/s'
                hom(:,2,164,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'w*p*' )
             dopr_index(i) = 68
             dopr_unit(i)  = 'm3/s3'
             hom(:,2,68,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'w"e' )
             dopr_index(i) = 69
             dopr_unit(i)  = 'm3/s3'
             hom(:,2,69,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'q*2' )
             IF (  .NOT.  humidity )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for humidity = .FALSE.'
                CALL message( 'check_parameters', 'PA0092', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 70
                dopr_unit(i)  = 'kg2/kg2'
                hom(:,2,70,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'prho' )
             IF (  .NOT.  ocean ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for ocean = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 71
                dopr_unit(i)  = 'kg/m3'
                hom(:,2,71,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'hyp' )
             dopr_index(i) = 72
             dopr_unit(i)  = 'hPa'
             hom(:,2,72,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'rho_air' )
             dopr_index(i)  = 119
             dopr_unit(i)   = 'kg/m3'
             hom(:,2,119,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'rho_air_zw' )
             dopr_index(i)  = 120
             dopr_unit(i)   = 'kg/m3'
             hom(:,2,120,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'nc' )
             IF (  .NOT.  cloud_physics )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for cloud_physics = .FALSE.'
                CALL message( 'check_parameters', 'PA0094', 1, 2, 0, 6, 0 )
             ELSEIF ( .NOT.  microphysics_morrison )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for cloud_scheme /= morrison'
                CALL message( 'check_parameters', 'PA0358', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 89
                dopr_unit(i)  = '1/m3'
                hom(:,2,89,:)  = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'nr' )
             IF (  .NOT.  cloud_physics )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for cloud_physics = .FALSE.'
                CALL message( 'check_parameters', 'PA0094', 1, 2, 0, 6, 0 )
             ELSEIF ( .NOT.  microphysics_seifert )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for cloud_scheme /= seifert_beheng'
                CALL message( 'check_parameters', 'PA0358', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 73
                dopr_unit(i)  = '1/m3'
                hom(:,2,73,:)  = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'qr' )
             IF (  .NOT.  cloud_physics )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for cloud_physics = .FALSE.'
                CALL message( 'check_parameters', 'PA0094', 1, 2, 0, 6, 0 )
             ELSEIF ( .NOT.  microphysics_seifert )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for cloud_scheme /= seifert_beheng'
                CALL message( 'check_parameters', 'PA0358', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 74
                dopr_unit(i)  = 'kg/kg'
                hom(:,2,74,:)  = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'qc' )
             IF (  .NOT.  cloud_physics )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for cloud_physics = .FALSE.'
                CALL message( 'check_parameters', 'PA0094', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 75
                dopr_unit(i)  = 'kg/kg'
                hom(:,2,75,:)  = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'prr' )
             IF (  .NOT.  cloud_physics )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for cloud_physics = .FALSE.'
                CALL message( 'check_parameters', 'PA0094', 1, 2, 0, 6, 0 )
             ELSEIF ( microphysics_sat_adjust )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not ava' // &
                                 'ilable for cloud_scheme = saturation_adjust'
                CALL message( 'check_parameters', 'PA0422', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 76
                dopr_unit(i)  = 'kg/kg m/s'
                hom(:,2,76,:)  = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'ug' )
             dopr_index(i) = 78
             dopr_unit(i)  = 'm/s'
             hom(:,2,78,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'vg' )
             dopr_index(i) = 79
             dopr_unit(i)  = 'm/s'
             hom(:,2,79,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'w_subs' )
             IF (  .NOT.  large_scale_subsidence )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for large_scale_subsidence = .FALSE.'
                CALL message( 'check_parameters', 'PA0382', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 80
                dopr_unit(i)  = 'm/s'
                hom(:,2,80,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 's*2' )
             IF (  .NOT.  passive_scalar )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for passive_scalar = .FALSE.'
                CALL message( 'check_parameters', 'PA0185', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 116
                dopr_unit(i)  = 'kg2/m6'
                hom(:,2,116,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF



          CASE DEFAULT

             CALL lsm_check_data_output_pr( data_output_pr(i), i, unit,        &
                                            dopr_unit(i) )

             IF ( unit == 'illegal' )  THEN
                CALL lsf_nudging_check_data_output_pr( data_output_pr(i), i,   &
                                                    unit, dopr_unit(i) )
             ENDIF

             IF ( unit == 'illegal' )  THEN
                CALL radiation_check_data_output_pr( data_output_pr(i), i,     &
                                                     unit, dopr_unit(i) )
             ENDIF
!
!--          Block of gust module profile outputs
             IF ( unit == 'illegal'  .AND.  gust_module_enabled  )  THEN
                   CALL gust_check_data_output_pr( data_output_pr(i), i, unit, &
                                                   dopr_unit(i) )
             ENDIF

             IF ( unit == 'illegal' )  THEN
                CALL chem_check_data_output_pr( data_output_pr(i), i,          &
                                                unit, dopr_unit(i) )
             ENDIF

             IF ( unit == 'illegal' )  THEN
                unit = ''
                CALL user_check_data_output_pr( data_output_pr(i), i, unit )
             ENDIF

             IF ( unit == 'illegal' )  THEN
                IF ( data_output_pr_user(1) /= ' ' )  THEN
                   message_string = 'illegal value for data_output_pr or ' //  &
                                    'data_output_pr_user = "' //               &
                                    TRIM( data_output_pr(i) ) // '"'
                   CALL message( 'check_parameters', 'PA0097', 1, 2, 0, 6, 0 )
                ELSE
                   message_string = 'illegal value for data_output_pr = "' //  &
                                    TRIM( data_output_pr(i) ) // '"'
                   CALL message( 'check_parameters', 'PA0098', 1, 2, 0, 6, 0 )
                ENDIF
             ENDIF

       END SELECT

    ENDDO


!
!-- Append user-defined data output variables to the standard data output
    IF ( data_output_user(1) /= ' ' )  THEN
       i = 1
       DO  WHILE ( data_output(i) /= ' '  .AND.  i <= 500 )
          i = i + 1
       ENDDO
       j = 1
       DO  WHILE ( data_output_user(j) /= ' '  .AND.  j <= 500 )
          IF ( i > 500 )  THEN
             message_string = 'number of output quantitities given by data' // &
                '_output and data_output_user exceeds the limit of 500'
             CALL message( 'check_parameters', 'PA0102', 1, 2, 0, 6, 0 )
          ENDIF
          data_output(i) = data_output_user(j)
          i = i + 1
          j = j + 1
       ENDDO
    ENDIF

!
!-- Check and set steering parameters for 2d/3d data output and averaging
    i   = 1
    DO  WHILE ( data_output(i) /= ' '  .AND.  i <= 500 )
!
!--    Check for data averaging
       ilen = LEN_TRIM( data_output(i) )
       j = 0                                                 ! no data averaging
       IF ( ilen > 3 )  THEN
          IF ( data_output(i)(ilen-2:ilen) == '_av' )  THEN
             j = 1                                           ! data averaging
             data_output(i) = data_output(i)(1:ilen-3)
          ENDIF
       ENDIF
!
!--    Check for cross section or volume data
       ilen = LEN_TRIM( data_output(i) )
       k = 0                                                   ! 3d data
       var = data_output(i)(1:ilen)
       IF ( ilen > 3 )  THEN
          IF ( data_output(i)(ilen-2:ilen) == '_xy'  .OR.                      &
               data_output(i)(ilen-2:ilen) == '_xz'  .OR.                      &
               data_output(i)(ilen-2:ilen) == '_yz' )  THEN

             k = 1                                             ! 2d data
             var = data_output(i)(1:ilen-3)

!--          Make sure that section_nn is defined
             IF ( (data_output(i)(ilen-2:ilen) == '_xy') .AND. ( ALL(section_xy == -9999) ) ) THEN
                message_string = 'to output _xy variables, the depth of at least one xy_section must be specified via the namelist parameter section_xy'
                CALL message( 'check_parameters', 'PA0561', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( (data_output(i)(ilen-2:ilen) == '_xz') .AND. ( ALL(section_xz == -9999) ) ) THEN
                message_string = 'to output _xz variables, the depth of at least one xz_section must be specified via the namelist parameter section_xz'
                CALL message( 'check_parameters', 'PA0561', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( (data_output(i)(ilen-2:ilen) == '_yz') .AND. ( ALL(section_yz == -9999) ) ) THEN
                message_string = 'to output _yz variables, the depth of at least one yz_section must be specified via the namelist parameter section_yz'
                CALL message( 'check_parameters', 'PA0561', 1, 2, 0, 6, 0 )
             ENDIF

          ENDIF

       ENDIF

!
!--    Check for allowed value and set units
       SELECT CASE ( TRIM( var ) )

          CASE ( 'e' )
             IF ( constant_diffusion )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res constant_diffusion = .FALSE.'
                CALL message( 'check_parameters', 'PA0103', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'm2/s2'

          CASE ( 'lpt' )
             IF (  .NOT.  cloud_physics )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                         'res cloud_physics = .TRUE.'
                CALL message( 'check_parameters', 'PA0108', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'K'

          CASE ( 'nc' )
             IF (  .NOT.  cloud_physics )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                         'res cloud_physics = .TRUE.'
                CALL message( 'check_parameters', 'PA0108', 1, 2, 0, 6, 0 )
             ELSEIF ( .NOT.  microphysics_morrison )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                         'res = morrison '
                CALL message( 'check_parameters', 'PA0359', 1, 2, 0, 6, 0 )
             ENDIF
             unit = '1/m3'

          CASE ( 'nr' )
             IF (  .NOT.  cloud_physics )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                         'res cloud_physics = .TRUE.'
                CALL message( 'check_parameters', 'PA0108', 1, 2, 0, 6, 0 )
             ELSEIF ( .NOT.  microphysics_seifert )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                         'res cloud_scheme = seifert_beheng'
                CALL message( 'check_parameters', 'PA0359', 1, 2, 0, 6, 0 )
             ENDIF
             unit = '1/m3'

          CASE ( 'pc', 'pr' )
             IF (  .NOT.  particle_advection )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requir' // &
                   'es a "particle_parameters"-NAMELIST in the parameter ' //  &
                   'file (PARIN)'
                CALL message( 'check_parameters', 'PA0104', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( TRIM( var ) == 'pc' )  unit = 'number'
             IF ( TRIM( var ) == 'pr' )  unit = 'm'

          CASE ( 'prr' )
             IF (  .NOT.  cloud_physics )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                         'res cloud_physics = .TRUE.'
                CALL message( 'check_parameters', 'PA0108', 1, 2, 0, 6, 0 )
             ELSEIF ( microphysics_sat_adjust )  THEN
                message_string = 'output of "' // TRIM( var ) // '" is ' //    &
                         'not available for cloud_scheme = saturation_adjust'
                CALL message( 'check_parameters', 'PA0423', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg m/s'

          CASE ( 'q', 'vpt' )
             IF (  .NOT.  humidity )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res humidity = .TRUE.'
                CALL message( 'check_parameters', 'PA0105', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( TRIM( var ) == 'q'   )  unit = 'kg/kg'
             IF ( TRIM( var ) == 'vpt' )  unit = 'K'

          CASE ( 'qc' )
             IF (  .NOT.  cloud_physics )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                         'res cloud_physics = .TRUE.'
                CALL message( 'check_parameters', 'PA0108', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg'

          CASE ( 'ql' )
             IF ( .NOT.  ( cloud_physics  .OR.  cloud_droplets ) )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                         'res cloud_physics = .TRUE. or cloud_droplets = .TRUE.'
                CALL message( 'check_parameters', 'PA0106', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg'

          CASE ( 'ql_c', 'ql_v', 'ql_vp' )
             IF (  .NOT.  cloud_droplets )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res cloud_droplets = .TRUE.'
                CALL message( 'check_parameters', 'PA0107', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( TRIM( var ) == 'ql_c'  )  unit = 'kg/kg'
             IF ( TRIM( var ) == 'ql_v'  )  unit = 'm3'
             IF ( TRIM( var ) == 'ql_vp' )  unit = 'none'

          CASE ( 'qr' )
             IF (  .NOT.  cloud_physics )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                         'res cloud_physics = .TRUE.'
                CALL message( 'check_parameters', 'PA0108', 1, 2, 0, 6, 0 )
             ELSEIF ( .NOT.  microphysics_seifert ) THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                         'res cloud_scheme = seifert_beheng'
                CALL message( 'check_parameters', 'PA0359', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg'

          CASE ( 'qv' )
             IF (  .NOT.  cloud_physics )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res cloud_physics = .TRUE.'
                CALL message( 'check_parameters', 'PA0108', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg'

          CASE ( 'solar3d' )
             IF (  .NOT.  ocean )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res ocean = .TRUE.'
                CALL message( 'check_parameters', 'PA0109', 1, 2, 0, 6, 0 )
             ENDIF
             unit = '^oC s^{-1}'


          CASE ( 'rho_ocean' )
             IF (  .NOT.  ocean )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res ocean = .TRUE.'
                CALL message( 'check_parameters', 'PA0109', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/m3'

          CASE ( 'alpha_T' )
             IF (  .NOT.  ocean )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res ocean = .TRUE.'
                CALL message( 'check_parameters', 'PA0109', 1, 2, 0, 6, 0 )
             ENDIF
             unit = '^oC^{-1}'

          CASE ( 'beta_S' )
             IF (  .NOT.  ocean )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res ocean = .TRUE.'
                CALL message( 'check_parameters', 'PA0109', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'PSU^{-1}'

          CASE ( 's' )
             IF (  .NOT.  passive_scalar )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res passive_scalar = .TRUE.'
                CALL message( 'check_parameters', 'PA0110', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/m3'

          CASE ( 'sa' )
             IF (  .NOT.  ocean )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res ocean = .TRUE.'
                CALL message( 'check_parameters', 'PA0109', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'psu'

          CASE ( 'p', 'pt', 'u', 'v', 'w' )
             IF ( TRIM( var ) == 'p'  )  unit = 'Pa'
             IF ( TRIM( var ) == 'pt' )  unit = 'K'
             IF ( TRIM( var ) == 'u'  )  unit = 'm/s'
             IF ( TRIM( var ) == 'v'  )  unit = 'm/s'
             IF ( TRIM( var ) == 'w'  )  unit = 'm/s'
             CONTINUE

          CASE ( 'ghf*', 'lwp*', 'ol*', 'pra*', 'prr*', 'qsws*', 'r_a*',       &
                 'shf_sol*', 'shf*', 'ssws*', 't*', 'tsurf*', 'u*', 'z0*', 'z0h*', 'z0q*' )
             IF ( k == 0  .OR.  data_output(i)(ilen-2:ilen) /= '_xy' )  THEN
                message_string = 'illegal value for data_output: "' //         &
                                 TRIM( var ) // '" & only 2d-horizontal ' //   &
                                 'cross sections are allowed for this value'
                CALL message( 'check_parameters', 'PA0111', 1, 2, 0, 6, 0 )
             ENDIF

             IF ( TRIM( var ) == 'lwp*'  .AND.  .NOT. cloud_physics )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res cloud_physics = .TRUE.'
                CALL message( 'check_parameters', 'PA0108', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( TRIM( var ) == 'pra*'  .AND.                                 &
                  .NOT. ( microphysics_kessler .OR. microphysics_seifert ) )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res cloud_scheme = kessler or seifert_beheng'
                CALL message( 'check_parameters', 'PA0112', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( TRIM( var ) == 'pra*'  .AND.  j == 1 )  THEN
                message_string = 'temporal averaging of precipitation ' //     &
                          'amount "' // TRIM( var ) // '" is not possible'
                CALL message( 'check_parameters', 'PA0113', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( TRIM( var ) == 'prr*'  .AND.                                 &
                  .NOT. ( microphysics_kessler .OR. microphysics_seifert ) )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res cloud_scheme = kessler or seifert_beheng'
                CALL message( 'check_parameters', 'PA0112', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( TRIM( var ) == 'qsws*'  .AND.  .NOT.  humidity )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res humidity = .TRUE.'
                CALL message( 'check_parameters', 'PA0322', 1, 2, 0, 6, 0 )
             ENDIF

             IF ( TRIM( var ) == 'ghf*'  .AND.  .NOT.  land_surface )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res land_surface = .TRUE.'
                CALL message( 'check_parameters', 'PA0404', 1, 2, 0, 6, 0 )
             ENDIF

             IF ( ( TRIM( var ) == 'r_a*' .OR.  TRIM( var ) == 'ghf*' )        &
                 .AND.  .NOT.  land_surface  .AND.  .NOT.  urban_surface )     &
             THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res land_surface = .TRUE. or ' //            &
                                 'urban_surface = .TRUE.'
                CALL message( 'check_parameters', 'PA0404', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( TRIM( var ) == 'ssws*'  .AND.  .NOT.  passive_scalar )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res passive_scalar = .TRUE.'
                CALL message( 'check_parameters', 'PA0361', 1, 2, 0, 6, 0 )
             ENDIF

             IF ( TRIM( var ) == 'ghf*'   )  unit = 'W/m2'
             IF ( TRIM( var ) == 'lwp*'   )  unit = 'kg/m2'
             IF ( TRIM( var ) == 'ol*'    )  unit = 'm'
             IF ( TRIM( var ) == 'pra*'   )  unit = 'mm'
             IF ( TRIM( var ) == 'prr*'   )  unit = 'mm/s'
             IF ( TRIM( var ) == 'qsws*'  )  unit = 'kgm/kgs'
             IF ( TRIM( var ) == 'r_a*'   )  unit = 's/m'
             IF ( TRIM( var ) == 'shf*'   )  unit = 'K*m/s'
             IF ( TRIM( var ) == 'shf_sol*' ) unit = 'K*m/s'
             IF ( TRIM( var ) == 'ssws*'  )  unit = 'kg/m2*s'
             IF ( TRIM( var ) == 't*'     )  unit = 'K'
             IF ( TRIM( var ) == 'tsurf*' )  unit = 'K'
             IF ( TRIM( var ) == 'u*'     )  unit = 'm/s'
             IF ( TRIM( var ) == 'z0*'    )  unit = 'm'
             IF ( TRIM( var ) == 'z0h*'   )  unit = 'm'
!
!--          Output of surface latent and sensible heat flux will be in W/m2
!--          in case of natural- and urban-type surfaces, even if
!--          flux_output_mode is set to kinematic units.
             IF ( land_surface  .OR.  urban_surface )  THEN
                IF ( TRIM( var ) == 'shf*'   )  unit = 'W/m2'
                IF ( TRIM( var ) == 'qsws*'  )  unit = 'W/m2'
             ENDIF

          CASE DEFAULT

             CALL tcm_check_data_output ( var, unit, i, ilen, k )

             IF ( unit == 'illegal' )  THEN
                CALL lsm_check_data_output ( var, unit, i, ilen, k )
             ENDIF

             IF ( unit == 'illegal' )  THEN
                CALL radiation_check_data_output( var, unit, i, ilen, k )
             ENDIF

!
!--          Block of gust module outputs
             IF ( unit == 'illegal'  .AND.  gust_module_enabled  )  THEN
                CALL gust_check_data_output ( var, unit )
             ENDIF

!
!--          Block of chemistry model outputs
             IF ( unit == 'illegal'  .AND.  air_chemistry                      &
                                     .AND.  var(1:3) == 'kc_' )  THEN
                CALL chem_check_data_output( var, unit, i, ilen, k )
             ENDIF

!
!--          Block of urban surface model outputs
             IF ( unit == 'illegal' .AND. urban_surface .AND. var(1:4) == 'usm_' ) THEN
                 CALL usm_check_data_output( var, unit )
             ENDIF

!
!--          Block of plant canopy model outputs
             IF ( unit == 'illegal' .AND. plant_canopy .AND. var(1:4) == 'pcm_' ) THEN
                 CALL pcm_check_data_output( var, unit )
             ENDIF

!
!--          Block of uv exposure model outputs
             IF ( unit == 'illegal' .AND. uv_exposure .AND. var(1:5) == 'uvem_' )  THEN
                CALL uvem_check_data_output( var, unit, i, ilen, k )
             ENDIF

             IF ( unit == 'illegal' )  THEN
                unit = ''
                CALL user_check_data_output( var, unit )
             ENDIF

             IF ( unit == 'illegal' )  THEN
                IF ( data_output_user(1) /= ' ' )  THEN
                   message_string = 'illegal value for data_output or ' //     &
                         'data_output_user = "' // TRIM( data_output(i) ) // '"'
                   CALL message( 'check_parameters', 'PA0114', 1, 2, 0, 6, 0 )
                ELSE
                   message_string = 'illegal value for data_output = "' //     &
                                    TRIM( data_output(i) ) // '"'
                   CALL message( 'check_parameters', 'PA0115', 1, 2, 0, 6, 0 )
                ENDIF
             ENDIF

       END SELECT
!
!--    Set the internal steering parameters appropriately
       IF ( k == 0 )  THEN
          do3d_no(j)              = do3d_no(j) + 1
          do3d(j,do3d_no(j))      = data_output(i)
          do3d_unit(j,do3d_no(j)) = unit
       ELSE
          do2d_no(j)              = do2d_no(j) + 1
          do2d(j,do2d_no(j))      = data_output(i)
          do2d_unit(j,do2d_no(j)) = unit
          IF ( data_output(i)(ilen-2:ilen) == '_xy' )  THEN
             data_output_xy(j) = .TRUE.
          ENDIF
          IF ( data_output(i)(ilen-2:ilen) == '_xz' )  THEN
             data_output_xz(j) = .TRUE.
          ENDIF
          IF ( data_output(i)(ilen-2:ilen) == '_yz' )  THEN
             data_output_yz(j) = .TRUE.
          ENDIF
       ENDIF

       IF ( j == 1 )  THEN
!
!--       Check, if variable is already subject to averaging
          found = .FALSE.
          DO  k = 1, doav_n
             IF ( TRIM( doav(k) ) == TRIM( var ) )  found = .TRUE.
          ENDDO

          IF ( .NOT. found )  THEN
             doav_n = doav_n + 1
             doav(doav_n) = var
          ENDIF
       ENDIF

       i = i + 1
    ENDDO

!
!-- Averaged 2d or 3d output requires that an averaging interval has been set
    IF ( doav_n > 0  .AND.  averaging_interval == 0.0_wp )  THEN
       WRITE( message_string, * )  'output of averaged quantity "',            &
                                   TRIM( doav(1) ), '_av" requires to set a ', &
                                   'non-zero averaging interval'
       CALL message( 'check_parameters', 'PA0323', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check sectional planes and store them in one shared array
    IF ( ANY( section_xy > nz + 1 ) )  THEN
       WRITE( message_string, * )  'section_xy must be <= nz + 1 = ', nz + 1
       CALL message( 'check_parameters', 'PA0319', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( ANY( section_xz > ny + 1 ) )  THEN
       WRITE( message_string, * )  'section_xz must be <= ny + 1 = ', ny + 1
       CALL message( 'check_parameters', 'PA0320', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( ANY( section_yz > nx + 1 ) )  THEN
       WRITE( message_string, * )  'section_yz must be <= nx + 1 = ', nx + 1
       CALL message( 'check_parameters', 'PA0321', 1, 2, 0, 6, 0 )
    ENDIF
    section(:,1) = section_xy
    section(:,2) = section_xz
    section(:,3) = section_yz

!
!-- Upper plot limit for 3D arrays
    IF ( nz_do3d == -9999 )  nz_do3d = nzt + 1

!
!-- Set output format string (used in header)
    SELECT CASE ( netcdf_data_format )
       CASE ( 1 )
          netcdf_data_format_string = 'netCDF classic'
       CASE ( 2 )
          netcdf_data_format_string = 'netCDF 64bit offset'
       CASE ( 3 )
          netcdf_data_format_string = 'netCDF4/HDF5'
       CASE ( 4 )
          netcdf_data_format_string = 'netCDF4/HDF5 classic'
       CASE ( 5 )
          netcdf_data_format_string = 'parallel netCDF4/HDF5'
       CASE ( 6 )
          netcdf_data_format_string = 'parallel netCDF4/HDF5 classic'

    END SELECT

!
!-- Check mask conditions
    DO mid = 1, max_masks
       IF ( data_output_masks(mid,1) /= ' '  .OR.                              &
            data_output_masks_user(mid,1) /= ' ' ) THEN
          masks = masks + 1
       ENDIF
    ENDDO

    IF ( masks < 0  .OR.  masks > max_masks )  THEN
       WRITE( message_string, * )  'illegal value: masks must be >= 0 and ',   &
            '<= ', max_masks, ' (=max_masks)'
       CALL message( 'check_parameters', 'PA0325', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( masks > 0 )  THEN
       mask_scale(1) = mask_scale_x
       mask_scale(2) = mask_scale_y
       mask_scale(3) = mask_scale_z
       IF ( ANY( mask_scale <= 0.0_wp ) )  THEN
          WRITE( message_string, * )                                           &
               'illegal value: mask_scale_x, mask_scale_y and mask_scale_z',   &
               'must be > 0.0'
          CALL message( 'check_parameters', 'PA0326', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Generate masks for masked data output
!--    Parallel netcdf output is not tested so far for masked data, hence
!--    netcdf_data_format is switched back to non-parallel output.
       netcdf_data_format_save = netcdf_data_format
       IF ( netcdf_data_format > 4 )  THEN
          IF ( netcdf_data_format == 5 ) netcdf_data_format = 3
          IF ( netcdf_data_format == 6 ) netcdf_data_format = 4
          message_string = 'netCDF file formats '//                            &
                           '5 (parallel netCDF 4) and ' //                     &
                           '6 (parallel netCDF 4 Classic model) '//            &
                           '& are currently not supported (not yet tested) ' //&
                           'for masked data. &Using respective non-parallel' //&
                           ' output for masked data.'
          CALL message( 'check_parameters', 'PA0383', 0, 0, 0, 6, 0 )
       ENDIF
       CALL init_masks
       netcdf_data_format = netcdf_data_format_save
    ENDIF

!
!-- Check the NetCDF data format
    IF ( netcdf_data_format > 2 )  THEN
#if defined( __netcdf4 )
       CONTINUE
#else
       message_string = 'netCDF: netCDF4 format requested but no ' //          &
                        'cpp-directive __netcdf4 given & switch '  //          &
                        'back to 64-bit offset format'
       CALL message( 'check_parameters', 'PA0171', 0, 1, 0, 6, 0 )
       netcdf_data_format = 2
#endif
    ENDIF
    IF ( netcdf_data_format > 4 )  THEN
#if defined( __netcdf4 ) && defined( __netcdf4_parallel )
       CONTINUE
#else
       message_string = 'netCDF: netCDF4 parallel output requested but no ' // &
                        'cpp-directive __netcdf4_parallel given & switch '   //&
                        'back to netCDF4 non-parallel output'
       CALL message( 'check_parameters', 'PA0099', 0, 1, 0, 6, 0 )
       netcdf_data_format = netcdf_data_format - 2
#endif
    ENDIF

!
!-- Calculate fixed number of output time levels for parallel netcdf output.
!-- The time dimension has to be defined as limited for parallel output,
!-- because otherwise the I/O performance drops significantly.
    IF ( netcdf_data_format > 4 )  THEN

!
!--    Check if any of the follwoing data output interval is 0.0s, which is
!--    not allowed for parallel output.
       CALL check_dt_do( dt_do3d,           'dt_do3d'           )
       CALL check_dt_do( dt_do2d_xy,        'dt_do2d_xy'        )
       CALL check_dt_do( dt_do2d_xz,        'dt_do2d_xz'        )
       CALL check_dt_do( dt_do2d_yz,        'dt_do2d_yz'        )
       CALL check_dt_do( dt_data_output_av, 'dt_data_output_av' )

!--    Set needed time levels (ntdim) to
!--    saved time levels + to be saved time levels.
       ntdim_3d(0) = do3d_time_count(0) + CEILING(                             &
                     ( end_time - MAX( skip_time_do3d,                         &
                                       simulated_time_at_begin )               &
                     ) / dt_do3d )
       IF ( do3d_at_begin ) ntdim_3d(0) = ntdim_3d(0) + 1

       ntdim_3d(1) = do3d_time_count(1) + CEILING(                             &
                     ( end_time - MAX( skip_time_data_output_av,               &
                                       simulated_time_at_begin )               &
                     ) / dt_data_output_av )

       ntdim_2d_xy(0) = do2d_xy_time_count(0) + CEILING(                       &
                        ( end_time - MAX( skip_time_do2d_xy,                   &
                                          simulated_time_at_begin )            &
                        ) / dt_do2d_xy )

       ntdim_2d_xz(0) = do2d_xz_time_count(0) + CEILING(                       &
                        ( end_time - MAX( skip_time_do2d_xz,                   &
                                          simulated_time_at_begin )            &
                        ) / dt_do2d_xz )

       ntdim_2d_yz(0) = do2d_yz_time_count(0) + CEILING(                       &
                        ( end_time - MAX( skip_time_do2d_yz,                   &
                                          simulated_time_at_begin )            &
                        ) / dt_do2d_yz )

       IF ( do2d_at_begin )  THEN
          ntdim_2d_xy(0) = ntdim_2d_xy(0) + 1
          ntdim_2d_xz(0) = ntdim_2d_xz(0) + 1
          ntdim_2d_yz(0) = ntdim_2d_yz(0) + 1
       ENDIF

       ntdim_2d_xy(1) = do2d_xy_time_count(1) + CEILING(                       &
                        ( end_time - MAX( skip_time_data_output_av,            &
                                          simulated_time_at_begin )            &
                        ) / dt_data_output_av )

       ntdim_2d_xz(1) = do2d_xz_time_count(1) + CEILING(                       &
                        ( end_time - MAX( skip_time_data_output_av,            &
                                          simulated_time_at_begin )            &
                                 ) / dt_data_output_av )

       ntdim_2d_yz(1) = do2d_yz_time_count(1) + CEILING(                       &
                        ( end_time - MAX( skip_time_data_output_av,            &
                                          simulated_time_at_begin )            &
                        ) / dt_data_output_av )

    ENDIF

!
!-- Check, whether a constant diffusion coefficient shall be used
    IF ( km_constant /= -1.0_wp )  THEN
       IF ( km_constant < 0.0_wp )  THEN
          WRITE( message_string, * )  'km_constant = ', km_constant, ' < 0.0'
          CALL message( 'check_parameters', 'PA0121', 1, 2, 0, 6, 0 )
       ELSE
          IF ( prandtl_number < 0.0_wp )  THEN
             WRITE( message_string, * )  'prandtl_number = ', prandtl_number,  &
                                         ' < 0.0'
             CALL message( 'check_parameters', 'PA0122', 1, 2, 0, 6, 0 )
          ENDIF
          constant_diffusion = .TRUE.

          IF ( constant_flux_layer )  THEN
             message_string = 'constant_flux_layer is not allowed with fixed ' &
                              // 'value of km'
             CALL message( 'check_parameters', 'PA0123', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
    ENDIF

!
!-- In case of non-cyclic lateral boundaries and a damping layer for the
!-- potential temperature, check the width of the damping layer
    IF ( bc_lr /= 'cyclic' ) THEN
       IF ( pt_damping_width < 0.0_wp  .OR.                                    &
            pt_damping_width > REAL( (nx+1) * dx ) )  THEN
          message_string = 'pt_damping_width out of range'
          CALL message( 'check_parameters', 'PA0124', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( bc_ns /= 'cyclic' )  THEN
       IF ( pt_damping_width < 0.0_wp  .OR.                                    &
            pt_damping_width > REAL( (ny+1) * dy ) )  THEN
          message_string = 'pt_damping_width out of range'
          CALL message( 'check_parameters', 'PA0124', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Check value range for zeta = z/L
    IF ( zeta_min >= zeta_max )  THEN
       WRITE( message_string, * )  'zeta_min = ', zeta_min, ' must be less ',  &
                                   'than zeta_max = ', zeta_max
       CALL message( 'check_parameters', 'PA0125', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check random generator
    IF ( (random_generator /= 'system-specific'      .AND.                     &
          random_generator /= 'random-parallel'   )  .AND.                     &
          random_generator /= 'numerical-recipes' )  THEN
       message_string = 'unknown random generator: random_generator = "' //    &
                        TRIM( random_generator ) // '"'
       CALL message( 'check_parameters', 'PA0135', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Determine upper and lower hight level indices for random perturbations
    IF ( disturbance_level_b == -9999999.9_wp )  THEN
       IF ( ocean )  THEN
          disturbance_level_b     = zu((nzt*2)/3)
          disturbance_level_ind_b = ( nzt * 2 ) / 3
       ELSE
          disturbance_level_b     = zu(nzb+3)
          disturbance_level_ind_b = nzb + 3
       ENDIF
    ELSEIF ( disturbance_level_b < zu(3) )  THEN
       WRITE( message_string, * )  'disturbance_level_b = ',                   &
                           disturbance_level_b, ' must be >= ', zu(3), '(zu(3))'
       CALL message( 'check_parameters', 'PA0126', 1, 2, 0, 6, 0 )
    ELSEIF ( disturbance_level_b > zu(nzt-2) )  THEN
       WRITE( message_string, * )  'disturbance_level_b = ',                   &
                   disturbance_level_b, ' must be <= ', zu(nzt-2), '(zu(nzt-2))'
       CALL message( 'check_parameters', 'PA0127', 1, 2, 0, 6, 0 )
    ELSE
       DO  k = 3, nzt-2
          IF ( disturbance_level_b <= zu(k) )  THEN
             disturbance_level_ind_b = k
             EXIT
          ENDIF
       ENDDO
    ENDIF

    IF ( disturbance_level_t == -9999999.9_wp )  THEN
       IF ( ocean )  THEN
          disturbance_level_t     = zu(nzt-3)
          disturbance_level_ind_t = nzt - 3
       ELSE
          disturbance_level_t     = zu(nzt/3)
          disturbance_level_ind_t = nzt / 3
       ENDIF
    ELSEIF ( disturbance_level_t > zu(nzt-2) )  THEN
       WRITE( message_string, * )  'disturbance_level_t = ',                   &
                   disturbance_level_t, ' must be <= ', zu(nzt-2), '(zu(nzt-2))'
       CALL message( 'check_parameters', 'PA0128', 1, 2, 0, 6, 0 )
    ELSEIF ( disturbance_level_t < disturbance_level_b )  THEN
       WRITE( message_string, * )  'disturbance_level_t = ',                   &
                   disturbance_level_t, ' must be >= disturbance_level_b = ',  &
                   disturbance_level_b
       CALL message( 'check_parameters', 'PA0129', 1, 2, 0, 6, 0 )
    ELSE
       DO  k = 3, nzt-2
          IF ( disturbance_level_t <= zu(k) )  THEN
             disturbance_level_ind_t = k
             EXIT
          ENDIF
       ENDDO
    ENDIF

!
!-- Check again whether the levels determined this way are ok.
!-- Error may occur at automatic determination and too few grid points in
!-- z-direction.
    IF ( disturbance_level_ind_t < disturbance_level_ind_b )  THEN
       WRITE( message_string, * )  'disturbance_level_ind_t = ',               &
                disturbance_level_ind_t, ' must be >= ',                       &
                'disturbance_level_ind_b = ', disturbance_level_ind_b
       CALL message( 'check_parameters', 'PA0130', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Determine the horizontal index range for random perturbations.
!-- In case of non-cyclic horizontal boundaries, no perturbations are imposed
!-- near the inflow and the perturbation area is further limited to ...(1)
!-- after the initial phase of the flow.

    IF ( bc_lr /= 'cyclic' )  THEN
       IF ( inflow_disturbance_begin == -1 )  THEN
          inflow_disturbance_begin = MIN( 10, nx/2 )
       ENDIF
       IF ( inflow_disturbance_begin < 0  .OR.  inflow_disturbance_begin > nx )&
       THEN
          message_string = 'inflow_disturbance_begin out of range'
          CALL message( 'check_parameters', 'PA0131', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( inflow_disturbance_end == -1 )  THEN
          inflow_disturbance_end = MIN( 100, 3*nx/4 )
       ENDIF
       IF ( inflow_disturbance_end < 0  .OR.  inflow_disturbance_end > nx )    &
       THEN
          message_string = 'inflow_disturbance_end out of range'
          CALL message( 'check_parameters', 'PA0132', 1, 2, 0, 6, 0 )
       ENDIF
    ELSEIF ( bc_ns /= 'cyclic' )  THEN
       IF ( inflow_disturbance_begin == -1 )  THEN
          inflow_disturbance_begin = MIN( 10, ny/2 )
       ENDIF
       IF ( inflow_disturbance_begin < 0  .OR.  inflow_disturbance_begin > ny )&
       THEN
          message_string = 'inflow_disturbance_begin out of range'
          CALL message( 'check_parameters', 'PA0131', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( inflow_disturbance_end == -1 )  THEN
          inflow_disturbance_end = MIN( 100, 3*ny/4 )
       ENDIF
       IF ( inflow_disturbance_end < 0  .OR.  inflow_disturbance_end > ny )    &
       THEN
          message_string = 'inflow_disturbance_end out of range'
          CALL message( 'check_parameters', 'PA0132', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( random_generator == 'random-parallel' )  THEN
       dist_nxl = nxl;  dist_nxr = nxr
       dist_nys = nys;  dist_nyn = nyn
       IF ( bc_lr == 'radiation/dirichlet' )  THEN
          dist_nxr    = MIN( nx - inflow_disturbance_begin, nxr )
          dist_nxl(1) = MAX( nx - inflow_disturbance_end, nxl )
       ELSEIF ( bc_lr == 'dirichlet/radiation' )  THEN
          dist_nxl    = MAX( inflow_disturbance_begin, nxl )
          dist_nxr(1) = MIN( inflow_disturbance_end, nxr )
       ELSEIF ( bc_lr == 'nested'  .OR.  bc_lr == 'forcing' )  THEN
          dist_nxl    = MAX( inflow_disturbance_begin, nxl )
          dist_nxr    = MIN( nx - inflow_disturbance_begin, nxr )
       ENDIF
       IF ( bc_ns == 'dirichlet/radiation' )  THEN
          dist_nyn    = MIN( ny - inflow_disturbance_begin, nyn )
          dist_nys(1) = MAX( ny - inflow_disturbance_end, nys )
       ELSEIF ( bc_ns == 'radiation/dirichlet' )  THEN
          dist_nys    = MAX( inflow_disturbance_begin, nys )
          dist_nyn(1) = MIN( inflow_disturbance_end, nyn )
       ELSEIF ( bc_ns == 'nested'  .OR.  bc_ns == 'forcing' )  THEN
          dist_nys    = MAX( inflow_disturbance_begin, nys )
          dist_nyn    = MIN( ny - inflow_disturbance_begin, nyn )
       ENDIF
    ELSE
       dist_nxl = 0;  dist_nxr = nx
       dist_nys = 0;  dist_nyn = ny
       IF ( bc_lr == 'radiation/dirichlet' )  THEN
          dist_nxr    = nx - inflow_disturbance_begin
          dist_nxl(1) = nx - inflow_disturbance_end
       ELSEIF ( bc_lr == 'dirichlet/radiation' )  THEN
          dist_nxl    = inflow_disturbance_begin
          dist_nxr(1) = inflow_disturbance_end
       ELSEIF ( bc_lr == 'nested'  .OR.  bc_lr == 'forcing' )  THEN
          dist_nxr    = nx - inflow_disturbance_begin
          dist_nxl    = inflow_disturbance_begin
       ENDIF
       IF ( bc_ns == 'dirichlet/radiation' )  THEN
          dist_nyn    = ny - inflow_disturbance_begin
          dist_nys(1) = ny - inflow_disturbance_end
       ELSEIF ( bc_ns == 'radiation/dirichlet' )  THEN
          dist_nys    = inflow_disturbance_begin
          dist_nyn(1) = inflow_disturbance_end
       ELSEIF ( bc_ns == 'nested'  .OR.  bc_ns == 'forcing' )  THEN
          dist_nyn    = ny - inflow_disturbance_begin
          dist_nys    = inflow_disturbance_begin
       ENDIF
    ENDIF

!
!-- A turbulent inflow requires Dirichlet conditions at the respective inflow
!-- boundary (so far, a turbulent inflow is realized from the left side only)
    IF ( turbulent_inflow  .AND.  bc_lr /= 'dirichlet/radiation' )  THEN
       message_string = 'turbulent_inflow = .T. requires a Dirichlet ' //      &
                        'condition at the inflow boundary'
       CALL message( 'check_parameters', 'PA0133', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Turbulent inflow requires that 3d arrays have been cyclically filled with
!-- data from prerun in the first main run
    IF ( turbulent_inflow  .AND.  initializing_actions /= 'cyclic_fill'  .AND. &
         initializing_actions /= 'read_restart_data' )  THEN
       message_string = 'turbulent_inflow = .T. requires ' //                  &
                        'initializing_actions = ''cyclic_fill'' or ' //        &
                        'initializing_actions = ''read_restart_data'' '
       CALL message( 'check_parameters', 'PA0055', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- In case of turbulent inflow calculate the index of the recycling plane
    IF ( turbulent_inflow )  THEN
       IF ( recycling_width <= dx  .OR.  recycling_width >= nx * dx )  THEN
          WRITE( message_string, * )  'illegal value for recycling_width: ',   &
                                      recycling_width
          CALL message( 'check_parameters', 'PA0134', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Calculate the index
       recycling_plane = recycling_width / dx
!
!--    Because the y-shift is done with a distance of INT( npey / 2 ) no shift
!--    is possible if there is only one PE in y direction.
       IF ( recycling_yshift .AND. pdims(2) < 2 )  THEN
          WRITE( message_string, * )  'recycling_yshift = .T. requires more',  &
                                      ' than one processor in y direction'
          CALL message( 'check_parameters', 'PA0421', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF


    IF ( turbulent_outflow )  THEN
!
!--    Turbulent outflow requires Dirichlet conditions at the respective inflow
!--    boundary (so far, a turbulent outflow is realized at the right side only)
       IF ( bc_lr /= 'dirichlet/radiation' )  THEN
          message_string = 'turbulent_outflow = .T. requires ' //              &
                           'bc_lr = "dirichlet/radiation"'
          CALL message( 'check_parameters', 'PA0038', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    The ouflow-source plane must lay inside the model domain
       IF ( outflow_source_plane < dx  .OR.  &
            outflow_source_plane > nx * dx )  THEN
          WRITE( message_string, * )  'illegal value for outflow_source'//     &
                                      '_plane: ', outflow_source_plane
          CALL message( 'check_parameters', 'PA0145', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Determine damping level index for 1D model
    IF ( INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN
       IF ( damp_level_1d == -1.0_wp )  THEN
          damp_level_1d     = zu(nzt+1)
          damp_level_ind_1d = nzt + 1
       ELSEIF ( damp_level_1d < 0.0_wp  .OR.  damp_level_1d > zu(nzt+1) )  THEN
          WRITE( message_string, * )  'damp_level_1d = ', damp_level_1d,       &
                 ' must be >= 0.0 and <= ', zu(nzt+1), '(zu(nzt+1))'
          CALL message( 'check_parameters', 'PA0136', 1, 2, 0, 6, 0 )
       ELSE
          DO  k = 1, nzt+1
             IF ( damp_level_1d <= zu(k) )  THEN
                damp_level_ind_1d = k
                EXIT
             ENDIF
          ENDDO
       ENDIF
    ENDIF

!
!-- Check some other 1d-model parameters
    IF ( TRIM( mixing_length_1d ) /= 'as_in_3d_model'  .AND.                   &
         TRIM( mixing_length_1d ) /= 'blackadar' )  THEN
       message_string = 'mixing_length_1d = "' // TRIM( mixing_length_1d ) //  &
                        '" is unknown'
       CALL message( 'check_parameters', 'PA0137', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( TRIM( dissipation_1d ) /= 'as_in_3d_model'  .AND.                     &
         TRIM( dissipation_1d ) /= 'detering'  .AND.                           &
         TRIM( dissipation_1d ) /= 'prognostic' )  THEN
       message_string = 'dissipation_1d = "' // TRIM( dissipation_1d ) //      &
                        '" is unknown'
       CALL message( 'check_parameters', 'PA0138', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( TRIM( mixing_length_1d ) /= 'as_in_3d_model'  .AND.                   &
         TRIM( dissipation_1d ) == 'as_in_3d_model' )  THEN
       message_string = 'dissipation_1d = "' // TRIM( dissipation_1d ) //      &
                        '" requires mixing_length_1d = "as_in_3d_model"'
       CALL message( 'check_parameters', 'PA0485', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set time for the next user defined restart (time_restart is the
!-- internal parameter for steering restart events)
    IF ( restart_time /= 9999999.9_wp )  THEN
       IF ( restart_time > time_since_reference_point )  THEN
          time_restart = restart_time
       ENDIF
    ELSE
!
!--    In case of a restart run, set internal parameter to default (no restart)
!--    if the NAMELIST-parameter restart_time is omitted
       time_restart = 9999999.9_wp
    ENDIF

!
!-- Check pressure gradient conditions
    IF ( dp_external  .AND.  conserve_volume_flow )  THEN
       WRITE( message_string, * )  'Both dp_external and conserve_volume_flo', &
            'w are .TRUE. but one of them must be .FALSE.'
       CALL message( 'check_parameters', 'PA0150', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( dp_external )  THEN
       IF ( dp_level_b < zu(nzb)  .OR.  dp_level_b > zu(nzt) )  THEN
          WRITE( message_string, * )  'dp_level_b = ', dp_level_b, ' is out ', &
               ' of range [zu(nzb), zu(nzt)]'
          CALL message( 'check_parameters', 'PA0151', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( .NOT. ANY( dpdxy /= 0.0_wp ) )  THEN
          WRITE( message_string, * )  'dp_external is .TRUE. but dpdxy is ze', &
               'ro, i.e. the external pressure gradient will not be applied'
          CALL message( 'check_parameters', 'PA0152', 0, 1, 0, 6, 0 )
       ENDIF
    ENDIF
    IF ( ANY( dpdxy /= 0.0_wp )  .AND.  .NOT.  dp_external )  THEN
       WRITE( message_string, * )  'dpdxy is nonzero but dp_external is ',     &
            '.FALSE., i.e. the external pressure gradient & will not be applied'
       CALL message( 'check_parameters', 'PA0153', 0, 1, 0, 6, 0 )
    ENDIF
    IF ( conserve_volume_flow )  THEN
       IF ( TRIM( conserve_volume_flow_mode ) == 'default' )  THEN

          conserve_volume_flow_mode = 'initial_profiles'

       ELSEIF ( TRIM( conserve_volume_flow_mode ) /= 'initial_profiles' .AND.  &
            TRIM( conserve_volume_flow_mode ) /= 'bulk_velocity' )  THEN
          WRITE( message_string, * )  'unknown conserve_volume_flow_mode: ',   &
               conserve_volume_flow_mode
          CALL message( 'check_parameters', 'PA0154', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( (bc_lr /= 'cyclic'  .OR.  bc_ns /= 'cyclic')  .AND.                &
          TRIM( conserve_volume_flow_mode ) == 'bulk_velocity' )  THEN
          WRITE( message_string, * )  'non-cyclic boundary conditions ',       &
               'require  conserve_volume_flow_mode = ''initial_profiles'''
          CALL message( 'check_parameters', 'PA0155', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF
    IF ( ( u_bulk /= 0.0_wp  .OR.  v_bulk /= 0.0_wp )  .AND.                   &
         ( .NOT. conserve_volume_flow  .OR.                                    &
         TRIM( conserve_volume_flow_mode ) /= 'bulk_velocity' ) )  THEN
       WRITE( message_string, * )  'nonzero bulk velocity requires ',          &
            'conserve_volume_flow = .T. and ',                                 &
            'conserve_volume_flow_mode = ''bulk_velocity'''
       CALL message( 'check_parameters', 'PA0157', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check particle attributes
    IF ( particle_color /= 'none' )  THEN
       IF ( particle_color /= 'absuv'  .AND.  particle_color /= 'pt*'  .AND.   &
            particle_color /= 'z' )  THEN
          message_string = 'illegal value for parameter particle_color: ' //   &
                           TRIM( particle_color)
          CALL message( 'check_parameters', 'PA0313', 1, 2, 0, 6, 0 )
       ELSE
          IF ( color_interval(2) <= color_interval(1) )  THEN
             message_string = 'color_interval(2) <= color_interval(1)'
             CALL message( 'check_parameters', 'PA0315', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
    ENDIF

    IF ( particle_dvrpsize /= 'none' )  THEN
       IF ( particle_dvrpsize /= 'absw' )  THEN
          message_string = 'illegal value for parameter particle_dvrpsize:' // &
                           ' ' // TRIM( particle_dvrpsize)
          CALL message( 'check_parameters', 'PA0314', 1, 2, 0, 6, 0 )
       ELSE
          IF ( dvrpsize_interval(2) <= dvrpsize_interval(1) )  THEN
             message_string = 'dvrpsize_interval(2) <= dvrpsize_interval(1)'
             CALL message( 'check_parameters', 'PA0316', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
    ENDIF

!
!-- Prevent empty time records in volume, cross-section and masked data in case
!-- of non-parallel netcdf-output in restart runs
    IF ( netcdf_data_format < 5 )  THEN
       IF ( TRIM( initializing_actions ) == 'read_restart_data' )  THEN
          do3d_time_count    = 0
          do2d_xy_time_count = 0
          do2d_xz_time_count = 0
          do2d_yz_time_count = 0
          domask_time_count  = 0
       ENDIF
    ENDIF

!
!-- Check for valid setting of most_method
    IF ( TRIM( most_method ) /= 'circular'  .AND.                              &
         TRIM( most_method ) /= 'newton'    .AND.                              &
         TRIM( most_method ) /= 'lookup' )  THEN
       message_string = 'most_method = "' // TRIM( most_method ) //            &
                        '" is unknown'
       CALL message( 'check_parameters', 'PA0416', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check roughness length, which has to be smaller than dz/2
    IF ( ( constant_flux_layer .OR.  &
           INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )       &
         .AND. roughness_length >= 0.5 * dz(1) )  THEN
       message_string = 'roughness_length must be smaller than dz/2'
       CALL message( 'check_parameters', 'PA0424', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Vertical nesting: check fine and coarse grid compatibility for data exchange
    IF ( vnested )  CALL vnest_check_parameters

!
!-- Check if topography is read from file in case of complex terrain simulations
    IF ( complex_terrain  .AND.  TRIM( topography ) /= 'read_from_file' )  THEN
       message_string = 'complex_terrain requires topography' //               &
                        ' = ''read_from_file'''
       CALL message( 'check_parameters', 'PA0295', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check if vertical grid stretching is switched off in case of complex
!-- terrain simulations
    IF ( complex_terrain  .AND.                                                &
         ANY( dz_stretch_level_start /= -9999999.9_wp ) )  THEN
       message_string = 'Vertical grid stretching is not allowed for ' //      &
                        'complex_terrain = .T.'
       CALL message( 'check_parameters', 'PA0473', 1, 2, 0, 6, 0 )
    ENDIF

    CALL location_message( 'finished', .TRUE. )
!
!-- Check &userpar parameters
    CALL user_check_parameters

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check the length of data output intervals. In case of parallel NetCDF output
!> the time levels of the output files need to be fixed. Therefore setting the
!> output interval to 0.0s (usually used to output each timestep) is not
!> possible as long as a non-fixed timestep is used.
!------------------------------------------------------------------------------!

    SUBROUTINE check_dt_do( dt_do, dt_do_name )

       IMPLICIT NONE

       CHARACTER (LEN=*), INTENT (IN) :: dt_do_name !< parin variable name

       REAL(wp), INTENT (INOUT)       :: dt_do      !< data output interval

       IF ( dt_do == 0.0_wp )  THEN
          IF ( dt_fixed )  THEN
             WRITE( message_string, '(A,F9.4,A)' )  'Output at every '  //     &
                    'timestep is wanted (' // dt_do_name // ' = 0.0).&'//      &
                    'The output interval is set to the fixed timestep dt '//   &
                    '= ', dt, 's.'
             CALL message( 'check_parameters', 'PA0060', 0, 0, 0, 6, 0 )
             dt_do = dt
          ELSE
             message_string = dt_do_name // ' = 0.0 while using a ' //         &
                              'variable timestep and parallel netCDF4 ' //     &
                              'is not allowed.'
             CALL message( 'check_parameters', 'PA0081', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

    END SUBROUTINE check_dt_do




!------------------------------------------------------------------------------!
! Description:
! ------------
!> Inititalizes the vertical profiles of scalar quantities.
!------------------------------------------------------------------------------!

    SUBROUTINE init_vertical_profiles( vertical_gradient_level_ind,            &
                                       vertical_gradient_level,                &
                                       vertical_gradient,                      &
                                       pr_init, surface_value, bc_t_val )


       IMPLICIT NONE

       INTEGER(iwp) ::  i     !< counter
       INTEGER(iwp), DIMENSION(1:10) ::  vertical_gradient_level_ind !< vertical grid indices for gradient levels

       REAL(wp)     ::  bc_t_val      !< model top gradient
       REAL(wp)     ::  gradient      !< vertica gradient of the respective quantity
       REAL(wp)     ::  surface_value !< surface value of the respecitve quantity


       REAL(wp), DIMENSION(0:nz+1) ::  pr_init                 !< initialisation profile
       REAL(wp), DIMENSION(1:10)   ::  vertical_gradient       !< given vertical gradient
       REAL(wp), DIMENSION(1:10)   ::  vertical_gradient_level !< given vertical gradient level

       i = 1
       gradient = 0.0_wp

       IF (  .NOT.  ocean )  THEN

          vertical_gradient_level_ind(1) = 0
          DO  k = 1, nzt+1
             IF ( i < 11 )  THEN
                IF ( vertical_gradient_level(i) < zu(k)  .AND.            &
                     vertical_gradient_level(i) >= 0.0_wp )  THEN
                   gradient = vertical_gradient(i) / 100.0_wp
                   vertical_gradient_level_ind(i) = k - 1
                   i = i + 1
                ENDIF
             ENDIF
             IF ( gradient /= 0.0_wp )  THEN
                IF ( k /= 1 )  THEN
                   pr_init(k) = pr_init(k-1) + dzu(k) * gradient
                ELSE
                   pr_init(k) = pr_init(k-1) + dzu(k) * gradient
                ENDIF
             ELSE
                pr_init(k) = pr_init(k-1)
             ENDIF
   !
   !--       Avoid negative values
             IF ( pr_init(k) < 0.0_wp )  THEN
                pr_init(k) = 0.0_wp
             ENDIF
          ENDDO

       ELSE

          vertical_gradient_level_ind(1) = nzt+1
          DO  k = nzt, 0, -1
             IF ( i < 11 )  THEN
                IF ( vertical_gradient_level(i) >= zu(k)  .AND.            &
                     vertical_gradient_level(i) <= 0.0_wp )  THEN
                   gradient = vertical_gradient(i) / 100.0_wp
                   vertical_gradient_level_ind(i) = k + 1
                   i = i + 1
                ENDIF
             ENDIF
             IF ( gradient /= 0.0_wp )  THEN
                IF ( k /= nzt )  THEN
                   pr_init(k) = pr_init(k+1) - dzu(k+1) * gradient
                ELSE
                   pr_init(k)   = surface_value - 0.5_wp * dzu(k+1) * gradient
                   pr_init(k+1) = surface_value + 0.5_wp * dzu(k+1) * gradient
                ENDIF
             ELSE
                pr_init(k) = pr_init(k+1)
             ENDIF
   !
   !--       Avoid negative humidities
             IF ( pr_init(k) < 0.0_wp )  THEN
                pr_init(k) = 0.0_wp
             ENDIF
          ENDDO

       ENDIF

!
!--    In case of no given gradients, choose zero gradient conditions
       IF ( vertical_gradient_level(1) == -999999.9_wp )  THEN
          vertical_gradient_level(1) = 0.0_wp
       ENDIF
!
!--    Store gradient at the top boundary for possible Neumann boundary condition
       bc_t_val  = ( pr_init(nzt+1) - pr_init(nzt) ) / dzu(nzt+1)

    END SUBROUTINE init_vertical_profiles



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the bottom and top boundary conditions for humidity and scalars.
!------------------------------------------------------------------------------!

    SUBROUTINE set_bc_scalars( sq, bc_b, bc_t, ibc_b, ibc_t, err_nr_b, err_nr_t )


       IMPLICIT NONE

       CHARACTER (LEN=1)   ::  sq         !< name of scalar quantity
       CHARACTER (LEN=*)   ::  bc_b       !< bottom boundary condition
       CHARACTER (LEN=*)   ::  bc_t       !< top boundary condition
       CHARACTER (LEN=*)   ::  err_nr_b   !< error number if bottom bc is unknown
       CHARACTER (LEN=*)   ::  err_nr_t   !< error number if top bc is unknown

       INTEGER(iwp)        ::  ibc_b      !< index for bottom boundary condition
       INTEGER(iwp)        ::  ibc_t      !< index for top boundary condition

!
!--    Set Integer flags and check for possilbe errorneous settings for bottom
!--    boundary condition
       IF ( bc_b == 'dirichlet' )  THEN
          ibc_b = 0
       ELSEIF ( bc_b == 'neumann' )  THEN
          ibc_b = 1
       ELSE
          message_string = 'unknown boundary condition: bc_' // TRIM( sq ) //  &
                           '_b ="' // TRIM( bc_b ) // '"'
          CALL message( 'check_parameters', err_nr_b, 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Set Integer flags and check for possilbe errorneous settings for top
!--    boundary condition
       IF ( bc_t == 'dirichlet' )  THEN
          ibc_t = 0
       ELSEIF ( bc_t == 'neumann' )  THEN
          ibc_t = 1
       ELSEIF ( bc_t == 'initial_gradient' )  THEN
          ibc_t = 2
       ELSEIF ( bc_t == 'nested'  .OR.  bc_t == 'forcing' )  THEN
          ibc_t = 3
       ELSE
          message_string = 'unknown boundary condition: bc_' // TRIM( sq ) //  &
                           '_t ="' // TRIM( bc_t ) // '"'
          CALL message( 'check_parameters', err_nr_t, 1, 2, 0, 6, 0 )
       ENDIF


    END SUBROUTINE set_bc_scalars



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check for consistent settings of bottom boundary conditions for humidity
!> and scalars.
!------------------------------------------------------------------------------!

    SUBROUTINE check_bc_scalars( sq, bc_b, ibc_b,                      &
                                 err_nr_1, err_nr_2,                   &
                                 constant_flux, surface_initial_change )


       IMPLICIT NONE

       CHARACTER (LEN=1)   ::  sq                       !< name of scalar quantity
       CHARACTER (LEN=*)   ::  bc_b                     !< bottom boundary condition
       CHARACTER (LEN=*)   ::  err_nr_1                 !< error number of first error
       CHARACTER (LEN=*)   ::  err_nr_2                 !< error number of second error

       INTEGER(iwp)        ::  ibc_b                    !< index of bottom boundary condition

       LOGICAL             ::  constant_flux            !< flag for constant-flux layer

       REAL(wp)            ::  surface_initial_change   !< value of initial change at the surface

!
!--    A given surface value implies Dirichlet boundary condition for
!--    the respective quantity. In this case specification of a constant flux is
!--    forbidden. However, an exception is made for large-scale forcing as well
!--    as land-surface model.
       IF ( .NOT. land_surface  .AND.  .NOT. large_scale_forcing )  THEN
          IF ( ibc_b == 0  .AND.  constant_flux )  THEN
             message_string = 'boundary condition: bc_' // TRIM( sq ) //       &
                              '_b ' // '= "' // TRIM( bc_b ) //                &
                              '" is not allowed with prescribed surface flux'
             CALL message( 'check_parameters', err_nr_1, 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
       IF ( constant_flux  .AND.  surface_initial_change /= 0.0_wp )  THEN
          WRITE( message_string, * )  'a prescribed surface flux is not allo', &
                 'wed with ', sq, '_surface_initial_change (/=0) = ',          &
                 surface_initial_change
          CALL message( 'check_parameters', err_nr_2, 1, 2, 0, 6, 0 )
       ENDIF


    END SUBROUTINE check_bc_scalars



 END SUBROUTINE check_parameters
