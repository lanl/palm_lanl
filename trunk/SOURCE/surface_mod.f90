!> @file surface_mod.f90
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
!
!------------------------------------------------------------------------------!
!
! Current revisions:
! ------------------
!
!
! Former revisions:
! -----------------
! $Id: surface_mod.f90 3055 2018-06-01 17:15:08Z suehring $
! Initialize Obukhov length, friction velocity and momentum fluxes also
! in case of restarts and cyclic fill
!
! 3026 2018-05-22 10:30:53Z schwenkel
! Changed the name specific humidity to mixing ratio, since we are computing
! mixing ratios.
!
! 2977 2018-04-17 10:27:57Z kanani
! Implement changes from branch radiation (r2948-2971) with minor modifications,
! plus some formatting.
! (moh.hefny)
! Added flag to check the existence of vertical urban/land surfaces, required
! to activate RTM
!
! 2970 2018-04-13 15:09:23Z suehring
! Remove un-necessary initialization of surface elements in old large-scale
! forcing mode
!
! 2963 2018-04-12 14:47:44Z suehring
! Introduce index for vegetation/wall, pavement/green-wall and water/window
! surfaces, for clearer access of surface fraction, albedo, emissivity, etc. .
!
! 2942 2018-04-03 13:51:09Z suehring
! Bugfix in assigning surface element data after restart
!
! 2940 2018-04-03 11:22:42Z suehring
! Bugfix in reading restart data of vertical surface elements
!
! 2920 2018-03-22 11:22:01Z kanani
! Correct comment for surface directions
!
! 2894 2018-03-15 09:17:58Z Giersch
! Calculations of the index range of the subdomain on file which overlaps with
! the current subdomain are already done in read_restart_data_mod,
! surface_read/write_restart_data was renamed to surface_r/wrd_local, variable
! named found has been introduced for checking if restart data was found,
! reading of restart strings has been moved completely to
! read_restart_data_mod, surface_rrd_local is already inside the overlap loop
! programmed in read_restart_data_mod, SAVE attribute added where necessary,
! deallocation and allocation of some arrays have been changed to take care of
! different restart files that can be opened (index i), the marker *** end
! surf *** is not necessary anymore, strings and their respective lengths are
! written out and read now in case of restart runs to get rid of prescribed
! character lengths (Giersch)
!
! 2813 2018-02-16 16:28:14Z suehring
! Some more bugfixes concerning restart runs
!
! 2812 2018-02-16 13:40:25Z hellstea
! Entries 'u_out', 'v_out' and 'w_out' removed from the functions
! get_topography_top_index and get_topography_top_index_ji
!
! 2805 2018-02-14 17:00:09Z suehring
! Bugfix in re-mapping surface elements in case of restart runs
!
! 2766 2018-01-22 17:17:47Z kanani
! Removed preprocessor directive __chem
!
! 2759 2018-01-17 16:24:59Z suehring
! Bugfix, consider density in vertical fluxes of passive scalar as well as
! chemical species.
!
! 2753 2018-01-16 14:16:49Z suehring
! +r_a_green, r_a_window
!
! 2718 2018-01-02 08:49:38Z maronga
! Changes from last commit documented
!
! 2706 2017-12-18 18:33:49Z suehring
! In case of restarts read and write pt_surface
!
! 2698 2017-12-14 18:46:24Z suehring
!
! 2696 2017-12-14 17:12:51Z kanani
! - Change in file header (GPL part)
! - Missing code for css added to surf_*, handling of surface_csflux updated (FK)
! - Bugfixes in reading/writing restart data in case several surface types are
!   present at the same time (MS)
! - Implementation of chemistry module (FK)
! - Allocation of pt1 and qv1 now done for all surface types (MS)
! - Revised classification of surface types
! - Introduce offset values to simplify index determination of surface elements
! - Dimensions of function get_topo_top_index (MS)
! - added variables for USM
! - adapted to changes in USM (RvT)
!
! 2688 2017-12-12 17:27:04Z Giersch
! Allocation and initialization of the latent heat flux (qsws) at the top of
! the ocean domain in case of coupled runs. In addtion, a double if-query has
! been removed.
!
! 2638 2017-11-23 12:44:23Z raasch
! bugfix for constant top momentumflux
!
! 2575 2017-10-24 09:57:58Z maronga
! Pavement parameterization revised
!
! 2547 2017-10-16 12:41:56Z schwenkel
! extended by cloud_droplets option
!
! 2508 2017-10-02 08:57:09Z suehring
! Minor formatting adjustment
!
! 2478 2017-09-18 13:37:24Z suehring
! Bugfixes in initializing model top
!
! 2378 2017-08-31 13:57:27Z suehring
! Bugfix in write restart data
!
! 2339 2017-08-07 13:55:26Z gronemeier
! corrected timestamp in header
!
! 2338 2017-08-07 12:15:38Z gronemeier
! Modularize 1D model
!
! 2318 2017-07-20 17:27:44Z suehring
! New function to obtain topography top index.
!
! 2317 2017-07-20 17:27:19Z suehring
! Implementation of new microphysic scheme: cloud_scheme = 'morrison'
! includes two more prognostic equations for cloud drop concentration (nc)
! and cloud water content (qc).
!
! 2270 2017-06-09 12:18:47Z maronga
! Parameters removed/added due to changes in the LSM
!
! 2269 2017-06-09 11:57:32Z suehring
! Formatting and description adjustments
!
! 2256 2017-06-07 13:58:08Z suehring
! Enable heating at downward-facing surfaces
!
! 2233 2017-05-30 18:08:54Z suehring
! Initial revision
!
!
! Description:
! ------------
!> Surface module defines derived data structures to treat surface-
!> bounded grid cells. Three different types of surfaces are defined:
!> default surfaces, natural surfaces, and urban surfaces. The module
!> encompasses the allocation and initialization of surface arrays, and handles
!> reading and writing restart data.
!> In addition, a further derived data structure is defined, in order to set
!> boundary conditions at surfaces.
!> @todo For the moment, downward-facing surfaces are only classified as
!>        default type
!> @todo Clean up urban-surface variables (some of them are not used any more)
!> @todo Revise chemistry surface flux part (reduce loops?!)
!------------------------------------------------------------------------------!
 MODULE surface_mod

    USE arrays_3d,                                                             &
        ONLY:  heatflux_input_conversion, momentumflux_input_conversion,       &
               alpha_T, beta_S, rho_air, rho_air_zw, zu, zw,                   &
               waterflux_input_conversion

    USE control_parameters

    USE constants

    USE indices,                                                               &
        ONLY:  nxl, nxlg, nxr, nxrg, nys, nysg, nyn, nyng, nzb, nzt, wall_flags_0

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE kinds


    IMPLICIT NONE

!
!-- Data type used to identify grid-points where horizontal boundary conditions
!-- are applied
    TYPE bc_type

       INTEGER(iwp) :: ns                                  !< number of surface elements on the PE

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  i       !< x-index linking to the PALM 3D-grid
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  j       !< y-index linking to the PALM 3D-grid
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  k       !< z-index linking to the PALM 3D-grid

       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE :: start_index !< start index within surface data type for given (j,i)
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE :: end_index   !< end index within surface data type for given (j,i)

    END TYPE bc_type
!
!-- Data type used to identify and treat surface-bounded grid points
    TYPE surf_type

       INTEGER(iwp) :: ioff                                !< offset value in x-direction, used to determine index of surface element
       INTEGER(iwp) :: joff                                !< offset value in y-direction, used to determine index of surface element
       INTEGER(iwp) :: koff                                !< offset value in z-direction, used to determine index of surface element
       INTEGER(iwp) :: ns                                  !< number of surface elements on the PE

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  i       !< x-index linking to the PALM 3D-grid
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  j       !< y-index linking to the PALM 3D-grid
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  k       !< z-index linking to the PALM 3D-grid

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  facing  !< Bit indicating surface orientation

       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE :: start_index !< Start index within surface data type for given (j,i)
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE :: end_index   !< End index within surface data type for given (j,i)

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z_mo      !< surface-layer height
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  uvw_abs   !< absolute surface-parallel velocity
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  us        !< friction velocity
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ts        !< scaling parameter temerature
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qs        !< scaling parameter humidity
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ss        !< scaling parameter passive scalar
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qcs       !< scaling parameter qc
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ncs       !< scaling parameter nc
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qrs       !< scaling parameter qr
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  nrs       !< scaling parameter nr

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ol        !< Obukhov length
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rib       !< Richardson bulk number

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z0        !< roughness length for momentum
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z0h       !< roughness length for heat
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z0q       !< roughness length for humidity

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt1       !< Potential temperature at first grid level
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qv1       !< mixing ratio at first grid level
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  css     !< scaling parameter chemical species
!
!--    Define arrays for surface fluxes
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  usws      !< vertical momentum flux for u-component at horizontal surfaces
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  vsws      !< vertical momentum flux for v-component at horizontal surfaces

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf       !< surface flux sensible heat
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws      !< surface flux latent heat
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ssws      !< surface flux passive scalar
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qcsws     !< surface flux qc
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ncsws     !< surface flux nc
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qrsws     !< surface flux qr
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  nrsws     !< surface flux nr
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  sasws     !< surface flux salinity
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_sol   !< surface solar heat flux -- based on shf and sasws
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  cssws   !< surface flux chemical species
!
!--    Required for horizontal walls in production_e
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  u_0       !< virtual velocity component (see production_e_init for further explanation)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  v_0       !< virtual velocity component (see production_e_init for further explanation)

       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  mom_flux_uv  !< momentum flux usvs and vsus at vertical surfaces (used in diffusion_u and diffusion_v)
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  mom_flux_w   !< momentum flux wsus and wsvs at vertical surfaces (used in diffusion_w)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  mom_flux_tke !< momentum flux usvs, vsus, wsus, wsvs at vertical surfaces at grid center (used in production_e)
!
!--    Variables required for LSM as well as for USM
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE   ::  nzt_pavement  !< top index for pavement in soil
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  albedo_type   !< albedo type, for each fraction (wall,green,window or vegetation,pavement water)

       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  building_surface    !< flag parameter indicating that the surface element is covered by buildings (no LSM actions, not implemented yet)
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  building_covered    !< flag indicating that buildings are on top of orography, only used for vertical surfaces in LSM
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  pavement_surface    !< flag parameter for pavements
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  water_surface       !< flag parameter for water surfaces
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  vegetation_surface  !< flag parameter for natural land surfaces

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  albedo            !< broadband albedo for each surface fraction (LSM: vegetation, water, pavement; USM: wall, green, window)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  emissivity        !< emissivity of the surface, for each fraction  (LSM: vegetation, water, pavement; USM: wall, green, window)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  frac              !< relative surface fraction (LSM: vegetation, water, pavement; USM: wall, green, window)

       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  aldif           !< albedo for longwave diffusive radiation, solar angle of 60°
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  aldir           !< albedo for longwave direct radiation, solar angle of 60°
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  asdif           !< albedo for shortwave diffusive radiation, solar angle of 60°
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  asdir           !< albedo for shortwave direct radiation, solar angle of 60°
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  rrtm_aldif      !< albedo for longwave diffusive radiation, solar angle of 60°
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  rrtm_aldir      !< albedo for longwave direct radiation, solar angle of 60°
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  rrtm_asdif      !< albedo for shortwave diffusive radiation, solar angle of 60°
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  rrtm_asdir      !< albedo for shortwave direct radiation, solar angle of 60°

       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  pt_surface        !< skin-surface temperature
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  rad_net           !< net radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  rad_net_l         !< net radiation, used in USM
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_h          !< heat conductivity of soil/ wall (W/m/K)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_h_green    !< heat conductivity of green soil (W/m/K)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_h_window   !< heat conductivity of windows (W/m/K)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_h_def      !< default heat conductivity of soil (W/m/K)

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_in           !< incoming longwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_out          !< emitted longwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_in           !< incoming shortwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_out          !< emitted shortwave radiation



       REAL(wp), DIMENSION(:), ALLOCATABLE ::  c_liq               !< liquid water coverage (of vegetated area)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  c_veg               !< vegetation coverage
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  f_sw_in             !< fraction of absorbed shortwave radiation by the surface layer (not implemented yet)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ghf                 !< ground heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  g_d                 !< coefficient for dependence of r_canopy on water vapour pressure deficit
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  lai                 !< leaf area index
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  lambda_surface_u    !< coupling between surface and soil (depends on vegetation type) (W/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  lambda_surface_s    !< coupling between surface and soil (depends on vegetation type) (W/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_liq            !< surface flux of latent heat (liquid water portion)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_soil           !< surface flux of latent heat (soil portion)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_veg            !< surface flux of latent heat (vegetation portion)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_a                 !< aerodynamic resistance
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_a_green           !< aerodynamic resistance at green fraction
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_a_window          !< aerodynamic resistance at window fraction
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_canopy            !< canopy resistance
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_soil              !< soil resistance
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_soil_min          !< minimum soil resistance
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_s                 !< total surface resistance (combination of r_soil and r_canopy)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_canopy_min        !< minimum canopy (stomatal) resistance

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  alpha_vg          !< coef. of Van Genuchten
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_w          !< hydraulic diffusivity of soil (?)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  gamma_w           !< hydraulic conductivity of soil (W/m/K)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  gamma_w_sat       !< hydraulic conductivity at saturation
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  l_vg              !< coef. of Van Genuchten
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  m_fc              !< soil moisture at field capacity (m3/m3)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  m_res             !< residual soil moisture
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  m_sat             !< saturation soil moisture (m3/m3)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  m_wilt            !< soil moisture at permanent wilting point (m3/m3)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  n_vg              !< coef. Van Genuchten
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rho_c_total_def   !< default volumetric heat capacity of the (soil) layer (J/m3/K)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rho_c_total       !< volumetric heat capacity of the actual soil matrix (J/m3/K)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  root_fr           !< root fraction within the soil layers
!
!--    Urban surface variables
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  surface_types   !< array of types of wall parameters

       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  isroof_surf          !< flag indicating roof surfaces
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  ground_level         !< flag indicating ground floor level surfaces

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  target_temp_summer  !< indoor target temperature summer
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  target_temp_winter  !< indoor target temperature summer

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  c_surface           !< heat capacity of the wall surface skin (J/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  c_surface_green     !< heat capacity of the green surface skin (J/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  c_surface_window    !< heat capacity of the window surface skin (J/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  lambda_surf         !< heat conductivity between air and surface (W/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  lambda_surf_green   !< heat conductivity between air and green surface (W/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  lambda_surf_window  !< heat conductivity between air and window surface (W/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  thickness_wall      !< thickness of the wall, roof and soil layers
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  thickness_green     !< thickness of the green wall, roof and soil layers
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  thickness_window    !< thickness of the window wall, roof and soil layers
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  transmissivity      !< transmissivity of windows

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutsl           !< reflected shortwave radiation for local surface in i-th reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutll           !< reflected + emitted longwave radiation for local surface in i-th reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfhf              !< total radiation flux incoming to minus outgoing from local surface

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  tt_surface_m        !< surface temperature tendency (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  tt_surface_window_m !< window surface temperature tendency (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  tt_surface_green_m  !< green surface temperature tendency (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wshf                !< kinematic wall heat flux of sensible heat (actually no longer needed)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wshf_eb             !< wall heat flux of sensible heat in wall normal direction

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wghf_eb             !< wall ground heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wghf_eb_window      !< window ground heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wghf_eb_green       !< green ground heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  iwghf_eb            !< indoor wall ground heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  iwghf_eb_window     !< indoor window ground heat flux

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_out_change_0

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinsw            !< shortwave radiation falling to local surface including radiation from reflections
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutsw           !< total shortwave radiation outgoing from nonvirtual surfaces surfaces after all reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinlw            !< longwave radiation falling to local surface including radiation from reflections
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutlw           !< total longwave radiation outgoing from nonvirtual surfaces surfaces after all reflection


       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rho_c_wall        !< volumetric heat capacity of the material ( J m-3 K-1 ) (= 2.19E6)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_wall           !< wall grid spacing (center-center)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ddz_wall          !< 1/dz_wall
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_wall_stag      !< wall grid spacing (edge-edge)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ddz_wall_stag     !< 1/dz_wall_stag
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tt_wall_m         !< t_wall prognostic array
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  zw                !< wall layer depths (m)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rho_c_window      !< volumetric heat capacity of the window material ( J m-3 K-1 ) (= 2.19E6)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_window         !< window grid spacing (center-center)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ddz_window        !< 1/dz_window
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_window_stag    !< window grid spacing (edge-edge)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ddz_window_stag   !< 1/dz_window_stag
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tt_window_m       !< t_window prognostic array
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  zw_window         !< window layer depths (m)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rho_c_green       !< volumetric heat capacity of the green material ( J m-3 K-1 ) (= 2.19E6)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_green          !< green grid spacing (center-center)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ddz_green         !< 1/dz_green
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_green_stag     !< green grid spacing (edge-edge)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ddz_green_stag    !< 1/dz_green_stag
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tt_green_m        !< t_green prognostic array
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  zw_green          !< green layer depths (m)


!-- arrays for time averages
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_net_av       !< average of rad_net_l
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinsw_av      !< average of sw radiation falling to local surface including radiation from reflections
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinlw_av      !< average of lw radiation falling to local surface including radiation from reflections
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinswdir_av   !< average of direct sw radiation falling to local surface
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinswdif_av   !< average of diffuse sw radiation from sky and model boundary falling to local surface
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinlwdif_av   !< average of diffuse lw radiation from sky and model boundary falling to local surface
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinswref_av   !< average of sw radiation falling to surface from reflections
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinlwref_av   !< average of lw radiation falling to surface from reflections
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutsw_av     !< average of total sw radiation outgoing from nonvirtual surfaces surfaces after all reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutlw_av     !< average of total lw radiation outgoing from nonvirtual surfaces surfaces after all reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfins_av       !< average of array of residua of sw radiation absorbed in surface after last reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinl_av       !< average of array of residua of lw radiation absorbed in surface after last reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfhf_av        !< average of total radiation flux incoming to minus outgoing from local surface
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wghf_eb_av       !< average of wghf_eb
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wghf_eb_window_av  !< average of wghf_eb window
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wghf_eb_green_av   !< average of wghf_eb window
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  iwghf_eb_av        !< indoor average of wghf_eb
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  iwghf_eb_window_av !< indoor average of wghf_eb window
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wshf_eb_av       !< average of wshf_eb
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_surf_av        !< average of wall surface temperature (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_surf_window_av !< average of window surface temperature (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_surf_green_av  !< average of green wall surface temperature (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_surf_whole_av  !< average of whole wall surface temperature (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_surf_10cm_av   !< average of the near surface temperature (K)

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  t_wall_av      !< Average of t_wall
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  t_window_av    !< Average of t_window
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  t_green_av     !< Average of t_green

    END TYPE surf_type

    TYPE (bc_type), DIMENSION(0:1)           ::  bc_h        !< boundary condition data type, horizontal upward- and downward facing surfaces

    TYPE (surf_type), DIMENSION(0:2), TARGET ::  surf_def_h  !< horizontal default surfaces (Up, Down, and Top)
    TYPE (surf_type), DIMENSION(0:3), TARGET ::  surf_def_v  !< vertical default surfaces (North, South, East, West)

    INTEGER(iwp), PARAMETER ::  ind_veg_wall  = 0            !< index for vegetation / wall-surface fraction, used for access of albedo, emissivity, etc., for each surface type
    INTEGER(iwp), PARAMETER ::  ind_pav_green = 1            !< index for pavement / green-wall surface fraction, used for access of albedo, emissivity, etc., for each surface type
    INTEGER(iwp), PARAMETER ::  ind_wat_win   = 2            !< index for water / window-surface fraction, used for access of albedo, emissivity, etc., for each surface type

    INTEGER(iwp) ::  ns_h_on_file(0:2)                       !< total number of horizontal surfaces with the same facing, required for writing restart data
    INTEGER(iwp) ::  ns_v_on_file(0:3)                       !< total number of vertical surfaces with the same facing, required for writing restart data

    LOGICAL ::  vertical_surfaces_exist = .FALSE.   !< flag indicating that there are vertical urban/land surfaces
                                                    !< in the domain (required to activiate RTM)


    SAVE

    PRIVATE

    INTERFACE deallocate_bc
       MODULE PROCEDURE deallocate_bc
    END INTERFACE deallocate_bc

    INTERFACE init_bc
       MODULE PROCEDURE init_bc
    END INTERFACE init_bc

    INTERFACE init_surfaces
       MODULE PROCEDURE init_surfaces
    END INTERFACE init_surfaces

    INTERFACE init_surface_arrays
       MODULE PROCEDURE init_surface_arrays
    END INTERFACE init_surface_arrays

    INTERFACE surface_rrd_local
       MODULE PROCEDURE surface_rrd_local
    END INTERFACE surface_rrd_local

    INTERFACE surface_wrd_local
       MODULE PROCEDURE surface_wrd_local
    END INTERFACE surface_wrd_local

    INTERFACE surface_last_actions
       MODULE PROCEDURE surface_last_actions
    END INTERFACE surface_last_actions

    INTERFACE surface_restore_elements
       MODULE PROCEDURE surface_restore_elements_1d
       MODULE PROCEDURE surface_restore_elements_2d
    END INTERFACE surface_restore_elements

!
!-- Public variables
    PUBLIC bc_h, ind_pav_green, ind_veg_wall, ind_wat_win, ns_h_on_file,       &
           ns_v_on_file, surf_def_h, surf_def_v,      &
           surf_type, vertical_surfaces_exist
!
!-- Public subroutines and functions
    PUBLIC get_topography_top_index, get_topography_top_index_ji, init_bc,     &
           init_surfaces, deallocate_bc,                                                     &
           init_surface_arrays, surface_rrd_local,                     &
           surface_restore_elements, surface_wrd_local,               &
           surface_last_actions


 CONTAINS

    SUBROUTINE deallocate_bc

       DEALLOCATE( bc_h(0)%i )
       DEALLOCATE( bc_h(0)%j )
       DEALLOCATE( bc_h(0)%k )
       DEALLOCATE( bc_h(0)%start_index )
       DEALLOCATE( bc_h(0)%end_index )
!
!--    Downward facing
       DEALLOCATE( bc_h(1)%i )
       DeALLOCATE( bc_h(1)%j)
       DEALLOCATE( bc_h(1)%k )
       DEALLOCATE( bc_h(1)%start_index )
       DEALLOCATE( bc_h(1)%end_index )

    end subroutine deallocate_bc
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize data type for setting boundary conditions at horizontal surfaces.
!------------------------------------------------------------------------------!
    SUBROUTINE init_bc

       IMPLICIT NONE

       INTEGER(iwp) ::  i         !<
       INTEGER(iwp) ::  j         !<
       INTEGER(iwp) ::  k         !<

       INTEGER(iwp), DIMENSION(0:1) ::  num_h         !<
       INTEGER(iwp), DIMENSION(0:1) ::  num_h_kji     !<
       INTEGER(iwp), DIMENSION(0:1) ::  start_index_h !<

!
!--    First of all, count the number of upward- and downward-facing surfaces
       num_h = 0
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb+1, nzt
!
!--             Check if current gridpoint belongs to the atmosphere
                IF ( BTEST( wall_flags_0(k,j,i), 0 ) )  THEN
!
!--                Upward-facing
                   IF ( .NOT. BTEST( wall_flags_0(k-1,j,i), 0 ) )              &
                      num_h(0) = num_h(0) + 1
!
!--                Downward-facing
                   IF ( .NOT. BTEST( wall_flags_0(k+1,j,i), 0 ) )              &
                      num_h(1) = num_h(1) + 1
                ENDIF
             ENDDO
          ENDDO
       ENDDO
!
!--    Save the number of surface elements
       bc_h(0)%ns = num_h(0)
       bc_h(1)%ns = num_h(1)
!
!--    ALLOCATE data type variables
!--    Upward facing
        ALLOCATE( bc_h(0)%i(1:bc_h(0)%ns) )
       ALLOCATE( bc_h(0)%j(1:bc_h(0)%ns) )
       ALLOCATE( bc_h(0)%k(1:bc_h(0)%ns) )
       ALLOCATE( bc_h(0)%start_index(nysg:nyng,nxlg:nxrg) )
       ALLOCATE( bc_h(0)%end_index(nysg:nyng,nxlg:nxrg)   )
       bc_h(0)%start_index = 1
       bc_h(0)%end_index   = 0
!
!--    Downward facing
       ALLOCATE( bc_h(1)%i(1:bc_h(1)%ns) )
       ALLOCATE( bc_h(1)%j(1:bc_h(1)%ns) )
       ALLOCATE( bc_h(1)%k(1:bc_h(1)%ns) )
       ALLOCATE( bc_h(1)%start_index(nysg:nyng,nxlg:nxrg) )
       ALLOCATE( bc_h(1)%end_index(nysg:nyng,nxlg:nxrg)   )
       bc_h(1)%start_index = 1
       bc_h(1)%end_index   = 0
!
!--    Store the respective indices on data type
       num_h(0:1)         = 1
       start_index_h(0:1) = 1
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng

             num_h_kji(0:1) = 0
             DO  k = nzb+1, nzt
!
!--             Check if current gridpoint belongs to the atmosphere
                IF ( BTEST( wall_flags_0(k,j,i), 0 ) )  THEN
!
!--                Upward-facing
                   IF ( .NOT. BTEST( wall_flags_0(k-1,j,i), 0 ) )  THEN
                      bc_h(0)%i(num_h(0)) = i
                      bc_h(0)%j(num_h(0)) = j
                      bc_h(0)%k(num_h(0)) = k
                      num_h_kji(0)        = num_h_kji(0) + 1
                      num_h(0)            = num_h(0) + 1
                   ENDIF
!
!--                Downward-facing
                   IF ( .NOT. BTEST( wall_flags_0(k+1,j,i), 0 ) )  THEN
                      bc_h(1)%i(num_h(1)) = i
                      bc_h(1)%j(num_h(1)) = j
                      bc_h(1)%k(num_h(1)) = k
                      num_h_kji(1)        = num_h_kji(1) + 1
                      num_h(1)            = num_h(1) + 1
                   ENDIF
                ENDIF
             ENDDO
             bc_h(0)%start_index(j,i) = start_index_h(0)
             bc_h(0)%end_index(j,i)   = bc_h(0)%start_index(j,i) + num_h_kji(0) - 1
             start_index_h(0)         = bc_h(0)%end_index(j,i) + 1

             bc_h(1)%start_index(j,i) = start_index_h(1)
             bc_h(1)%end_index(j,i)   = bc_h(1)%start_index(j,i) + num_h_kji(1) - 1
             start_index_h(1)         = bc_h(1)%end_index(j,i) + 1
          ENDDO
       ENDDO


    END SUBROUTINE init_bc


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize horizontal and vertical surfaces. Counts the number of default-,
!> natural and urban surfaces and allocates memory, respectively.
!------------------------------------------------------------------------------!
    SUBROUTINE init_surface_arrays


       USE pegrid


       IMPLICIT NONE

       INTEGER(iwp)                 ::  i         !< running index x-direction
       INTEGER(iwp)                 ::  j         !< running index y-direction
       INTEGER(iwp)                 ::  k         !< running index z-direction
       INTEGER(iwp)                 ::  l         !< index variable for surface facing
       INTEGER(iwp), DIMENSION(0:2) ::  num_def_h !< number of horizontally-aligned default surfaces
       INTEGER(iwp), DIMENSION(0:3) ::  num_def_v !< number of vertically-aligned default surfaces

       INTEGER(iwp)              ::  num_surf_v_l !< number of vertically-aligned local urban/land surfaces
       INTEGER(iwp)              ::  num_surf_v   !< number of vertically-aligned total urban/land surfaces


       num_def_h = 0
       num_def_v = 0
!
!--    Surfaces are classified according to the input data read from static
!--    input file. If no input file is present, all surfaces are classified
!--    either as natural, urban, or default, depending on the setting of
!--    land_surface and urban_surface. To control this, use the control
!--    flag topo_no_distinct
!
!--    Count number of horizontal surfaces on local domain
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Check if current gridpoint belongs to the atmosphere
                IF ( BTEST( wall_flags_0(k,j,i), 0 ) )  THEN
!
!--                Check if grid point adjoins to any upward-facing horizontal
!--                surface, e.g. the Earth surface, plane roofs, or ceilings.

                   IF ( .NOT. BTEST( wall_flags_0(k-1,j,i), 0 ) )  THEN
!
                      num_def_h(0) = num_def_h(0) + 1
!
                   ENDIF
!
!--                Check for top-fluxes
                   IF ( k == nzt  .AND.  use_top_fluxes )  THEN
                      num_def_h(2) = num_def_h(2) + 1
!
!--                Check for any other downward-facing surface. So far only for
!--                default surface type.
                   ELSEIF ( .NOT. BTEST( wall_flags_0(k+1,j,i), 0 ) )  THEN
                      num_def_h(1) = num_def_h(1) + 1
                   ENDIF

                ENDIF
             ENDDO
          ENDDO
       ENDDO
!
!--    Count number of vertical surfaces on local domain
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                IF ( BTEST( wall_flags_0(k,j,i), 0 ) )  THEN
!
!--                Northward-facing
                   IF ( .NOT. BTEST( wall_flags_0(k,j-1,i), 0 ) )  THEN
!
                         num_def_v(0) = num_def_v(0) + 1
!
                   ENDIF
!
!--                Southward-facing
                   IF ( .NOT. BTEST( wall_flags_0(k,j+1,i), 0 ) )  THEN
!
                      num_def_v(1) = num_def_v(1) + 1
!
                   ENDIF
!
!--                Eastward-facing
                   IF ( .NOT. BTEST( wall_flags_0(k,j,i-1), 0 ) )  THEN
!
                      num_def_v(2) = num_def_v(2) + 1
!
                   ENDIF
!
!--                Westward-facing
                   IF ( .NOT. BTEST( wall_flags_0(k,j,i+1), 0 ) )  THEN
!
                      num_def_v(3) = num_def_v(3) + 1
!
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO

!
!--    Store number of surfaces per core.
!--    Horizontal surface, default type, upward facing
       surf_def_h(0)%ns = num_def_h(0)
!
!--    Horizontal surface, default type, downward facing
       surf_def_h(1)%ns = num_def_h(1)
!
!--    Horizontal surface, default type, top downward facing
       surf_def_h(2)%ns = num_def_h(2)
!
!--    Vertical surface, default type, northward facing
       surf_def_v(0)%ns = num_def_v(0)
!
!--    Vertical surface, default type, southward facing
       surf_def_v(1)%ns = num_def_v(1)
!
!--    Vertical surface, default type, eastward facing
       surf_def_v(2)%ns = num_def_v(2)
!
!--    Vertical surface, default type, westward facing
       surf_def_v(3)%ns = num_def_v(3)
!
!--    Allocate required attributes for horizontal surfaces - default type.
!--    Upward-facing (l=0) and downward-facing (l=1).
       DO  l = 0, 1
          CALL allocate_surface_attributes_h ( surf_def_h(l), nys, nyn, nxl, nxr )
       ENDDO
!
!--    Allocate required attributes for model top
       CALL allocate_surface_attributes_h_top ( surf_def_h(2), nys, nyn, nxl, nxr )
!
!--    Allocate required attributes for vertical surfaces.
!--    Northward-facing (l=0), southward-facing (l=1), eastward-facing (l=2)
!--    and westward-facing (l=3).
!--    Default type.
       DO  l = 0, 3
          CALL allocate_surface_attributes_v ( surf_def_v(l),                  &
                                               nys, nyn, nxl, nxr )
       ENDDO
!

#if defined( __parallel )
       CALL MPI_ALLREDUCE( num_surf_v_l, num_surf_v, 1, MPI_INTEGER, MPI_SUM, comm2d, ierr)
#else
       num_surf_v = num_surf_v_l
#endif
       IF ( num_surf_v > 0 ) vertical_surfaces_exist = .TRUE.


    END SUBROUTINE init_surface_arrays


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Deallocating memory for upward and downward-facing horizontal surface types,
!> except for top fluxes.
!------------------------------------------------------------------------------!
    SUBROUTINE deallocate_surface_attributes_h( surfaces )

       IMPLICIT NONE


       TYPE(surf_type) ::  surfaces  !< respective surface type


       DEALLOCATE ( surfaces%start_index )
       DEALLOCATE ( surfaces%end_index )
!
!--    Indices to locate surface element
       DEALLOCATE ( surfaces%i )
       DEALLOCATE ( surfaces%j )
       DEALLOCATE ( surfaces%k )
!
!--    Surface-layer height
       DEALLOCATE ( surfaces%z_mo )
!
!--    Surface orientation
       DEALLOCATE ( surfaces%facing )
!
!--    Surface-parallel wind velocity
       DEALLOCATE ( surfaces%uvw_abs )
!
!--    Roughness
       DEALLOCATE ( surfaces%z0 )
       DEALLOCATE ( surfaces%z0h )
       DEALLOCATE ( surfaces%z0q )
!
!--    Friction velocity
       DEALLOCATE ( surfaces%us )
!
!--    Stability parameter
       DEALLOCATE ( surfaces%ol )
!
!--    Bulk Richardson number
       DEALLOCATE ( surfaces%rib )
!
!--    Vertical momentum fluxes of u and v
       DEALLOCATE ( surfaces%usws )
       DEALLOCATE ( surfaces%vsws )
!
!--    Required in production_e
       DEALLOCATE ( surfaces%u_0 )
       DEALLOCATE ( surfaces%v_0 )
!
!--    Characteristic temperature and surface flux of sensible heat
       DEALLOCATE ( surfaces%ts )
       DEALLOCATE ( surfaces%shf )
!
!--    surface temperature
       DEALLOCATE ( surfaces%pt1 )
       DEALLOCATE ( surfaces%pt_surface )
       DEALLOCATE ( surfaces%sasws )
       DEALLOCATE ( surfaces%shf_sol )

    END SUBROUTINE deallocate_surface_attributes_h


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocating memory for upward and downward-facing horizontal surface types,
!> except for top fluxes.
!------------------------------------------------------------------------------!
    SUBROUTINE allocate_surface_attributes_h( surfaces,                        &
                                              nys_l, nyn_l, nxl_l, nxr_l )

       IMPLICIT NONE

       INTEGER(iwp) ::  nyn_l  !< north bound of local 2d array start/end_index, is equal to nyn, except for restart-array
       INTEGER(iwp) ::  nys_l  !< south bound of local 2d array start/end_index, is equal to nyn, except for restart-array
       INTEGER(iwp) ::  nxl_l  !< west bound of local 2d array start/end_index, is equal to nyn, except for restart-array
       INTEGER(iwp) ::  nxr_l  !< east bound of local 2d array start/end_index, is equal to nyn, except for restart-array

       TYPE(surf_type) ::  surfaces  !< respective surface type


       IF ( ALLOCATED( surfaces%start_index ) )                      &
           CALL deallocate_surface_attributes_h( surfaces )


!
!--    Allocate arrays for start and end index of horizontal surface type
!--    for each (j,i)-grid point. This is required e.g. in diffion_x, which is
!--    called for each (j,i). In order to find the location where the
!--    respective flux is store within the surface-type, start- and end-
!--    index are stored for each (j,i). For example, each (j,i) can have
!--    several entries where fluxes for horizontal surfaces might be stored,
!--    e.g. for overhanging structures where several upward-facing surfaces
!--    might exist for given (j,i).
!--    If no surface of respective type exist at current (j,i), set indicies
!--    such that loop in diffusion routines will not be entered.
       ALLOCATE ( surfaces%start_index(nys_l:nyn_l,nxl_l:nxr_l) )
       ALLOCATE ( surfaces%end_index(nys_l:nyn_l,nxl_l:nxr_l)   )
       surfaces%start_index = 0
       surfaces%end_index   = -1
!
!--    Indices to locate surface element
       ALLOCATE ( surfaces%i(1:surfaces%ns)  )
       ALLOCATE ( surfaces%j(1:surfaces%ns)  )
       ALLOCATE ( surfaces%k(1:surfaces%ns)  )
!
!--    Surface-layer height
       ALLOCATE ( surfaces%z_mo(1:surfaces%ns) )
!
!--    Surface orientation
       ALLOCATE ( surfaces%facing(1:surfaces%ns) )
!
!--    Surface-parallel wind velocity
       ALLOCATE ( surfaces%uvw_abs(1:surfaces%ns) )
!
!--    Roughness
       ALLOCATE ( surfaces%z0(1:surfaces%ns)  )
       ALLOCATE ( surfaces%z0h(1:surfaces%ns) )
       ALLOCATE ( surfaces%z0q(1:surfaces%ns) )
!
!--    Friction velocity
       ALLOCATE ( surfaces%us(1:surfaces%ns) )
!
!--    Stability parameter
       ALLOCATE ( surfaces%ol(1:surfaces%ns) )
!
!--    Bulk Richardson number
       ALLOCATE ( surfaces%rib(1:surfaces%ns) )
!
!--    Vertical momentum fluxes of u and v
       ALLOCATE ( surfaces%usws(1:surfaces%ns) )
       ALLOCATE ( surfaces%vsws(1:surfaces%ns) )
!
!--    Required in production_e
       ALLOCATE ( surfaces%u_0(1:surfaces%ns) )
       ALLOCATE ( surfaces%v_0(1:surfaces%ns) )
!
!--    Characteristic temperature and surface flux of sensible heat
       ALLOCATE ( surfaces%ts(1:surfaces%ns)  )
       ALLOCATE ( surfaces%shf(1:surfaces%ns) )
!
!--    surface temperature
       ALLOCATE ( surfaces%pt_surface(1:surfaces%ns) )

!--    Arrays for storing potential temperature and
!--    mixing ratio at first grid level
       ALLOCATE ( surfaces%pt1(1:surfaces%ns) )
       ALLOCATE ( surfaces%sasws(1:surfaces%ns) )
       ALLOCATE ( surfaces%shf_sol(1:surfaces%ns) )

    END SUBROUTINE allocate_surface_attributes_h


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Deallocating memory for model-top fluxes
!------------------------------------------------------------------------------!
    SUBROUTINE deallocate_surface_attributes_h_top( surfaces )

       IMPLICIT NONE


       TYPE(surf_type) ::  surfaces !< respective surface type

       DEALLOCATE ( surfaces%start_index )
       DEALLOCATE ( surfaces%end_index )
!
!--    Indices to locate surface (model-top) element
       DEALLOCATE ( surfaces%i )
       DEALLOCATE ( surfaces%j )
       DEALLOCATE ( surfaces%k )

       DEALLOCATE ( surfaces%u_0 )
       DEALLOCATE ( surfaces%v_0 )
!
!--    Vertical momentum fluxes of u and v
       DEALLOCATE ( surfaces%usws )
       DEALLOCATE ( surfaces%vsws )
!
!--    Sensible heat flux
       DEALLOCATE ( surfaces%shf )
       DEALLOCATE ( surfaces%sasws )
       DEALLOCATE ( surfaces%shf_sol )

    END SUBROUTINE deallocate_surface_attributes_h_top


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocating memory for model-top fluxes
!------------------------------------------------------------------------------!
    SUBROUTINE allocate_surface_attributes_h_top( surfaces,                    &
                                                  nys_l, nyn_l, nxl_l, nxr_l )

       IMPLICIT NONE

       INTEGER(iwp) ::  nyn_l  !< north bound of local 2d array start/end_index, is equal to nyn, except for restart-array
       INTEGER(iwp) ::  nys_l  !< south bound of local 2d array start/end_index, is equal to nyn, except for restart-array
       INTEGER(iwp) ::  nxl_l  !< west bound of local 2d array start/end_index, is equal to nyn, except for restart-array
       INTEGER(iwp) ::  nxr_l  !< east bound of local 2d array start/end_index, is equal to nyn, except for restart-array

       TYPE(surf_type) ::  surfaces !< respective surface type

       IF ( ALLOCATED( surfaces%start_index ) )                      &
           CALL deallocate_surface_attributes_h_top( surfaces )

       ALLOCATE ( surfaces%start_index(nys_l:nyn_l,nxl_l:nxr_l) )
       ALLOCATE ( surfaces%end_index(nys_l:nyn_l,nxl_l:nxr_l)   )
       surfaces%start_index = 0
       surfaces%end_index   = -1
!
!--    Indices to locate surface (model-top) element
       ALLOCATE ( surfaces%i(1:surfaces%ns)  )
       ALLOCATE ( surfaces%j(1:surfaces%ns)  )
       ALLOCATE ( surfaces%k(1:surfaces%ns)  )

       ALLOCATE ( surfaces%u_0(1:surfaces%ns) )
       ALLOCATE ( surfaces%v_0(1:surfaces%ns) )
!
!--    Vertical momentum fluxes of u and v
       ALLOCATE ( surfaces%usws(1:surfaces%ns) )
       ALLOCATE ( surfaces%vsws(1:surfaces%ns) )
!
!--    Sensible heat flux
       ALLOCATE ( surfaces%shf(1:surfaces%ns) )
       ALLOCATE ( surfaces%sasws(1:surfaces%ns) )
       ALLOCATE ( surfaces%shf_sol(1:surfaces%ns) )

    END SUBROUTINE allocate_surface_attributes_h_top


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Deallocating memory for vertical surface types.
!------------------------------------------------------------------------------!
    SUBROUTINE deallocate_surface_attributes_v( surfaces )

       IMPLICIT NONE


       TYPE(surf_type) ::  surfaces !< respective surface type

!
!--    Allocate arrays for start and end index of vertical surface type
!--    for each (j,i)-grid point. This is required in diffion_x, which is
!--    called for each (j,i). In order to find the location where the
!--    respective flux is store within the surface-type, start- and end-
!--    index are stored for each (j,i). For example, each (j,i) can have
!--    several entries where fluxes for vertical surfaces might be stored.
!--    In the flat case, where no vertical walls exit, set indicies such
!--    that loop in diffusion routines will not be entered.
       DEALLOCATE ( surfaces%start_index )
       DEALLOCATE ( surfaces%end_index )
!
!--    Indices to locate surface element.
       DEALLOCATE ( surfaces%i )
       DEALLOCATE ( surfaces%j )
       DEALLOCATE ( surfaces%k )
!
!--    Surface-layer height
       DEALLOCATE ( surfaces%z_mo )
!
!--    Surface orientation
       DEALLOCATE ( surfaces%facing )
!
!--    Surface parallel wind velocity
       DEALLOCATE ( surfaces%uvw_abs )
!
!--    Roughness
       DEALLOCATE ( surfaces%z0 )
       DEALLOCATE ( surfaces%z0h )
       DEALLOCATE ( surfaces%z0q )

!
!--    Friction velocity
       DEALLOCATE ( surfaces%us )
!
!--    Allocate Obukhov length and bulk Richardson number. Actually, at
!--    vertical surfaces these are only required for natural surfaces.
!--    for natural land surfaces
       DEALLOCATE( surfaces%ol )
       DEALLOCATE( surfaces%rib )
!
!--    Allocate arrays for surface momentum fluxes for u and v. For u at north-
!--    and south-facing surfaces, for v at east- and west-facing surfaces.
       DEALLOCATE ( surfaces%mom_flux_uv )
!
!--    Allocate array for surface momentum flux for w - wsus and wsvs
       DEALLOCATE ( surfaces%mom_flux_w )
!
!--    Allocate array for surface momentum flux for subgrid-scale tke wsus and
!--    wsvs; first index usvs or vsws, second index for wsus or wsvs, depending
!--    on surface.
       DEALLOCATE ( surfaces%mom_flux_tke )
!
!--    Characteristic temperature and surface flux of sensible heat
       DEALLOCATE ( surfaces%ts )
       DEALLOCATE ( surfaces%shf )

!--    Arrays for storing potential temperature and
!--    mixing ratio at first grid level
       DEALLOCATE ( surfaces%pt1 )
       DEALLOCATE ( surfaces%qv1 )
       DEALLOCATE ( surfaces%sasws )
       DEALLOCATE ( surfaces%shf_sol )
       DEALLOCATE ( surfaces%pt_surface )
    END SUBROUTINE deallocate_surface_attributes_v


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocating memory for vertical surface types.
!------------------------------------------------------------------------------!
    SUBROUTINE allocate_surface_attributes_v( surfaces,                        &
                                              nys_l, nyn_l, nxl_l, nxr_l )

       IMPLICIT NONE

       INTEGER(iwp) ::  nyn_l  !< north bound of local 2d array start/end_index, is equal to nyn, except for restart-array
       INTEGER(iwp) ::  nys_l  !< south bound of local 2d array start/end_index, is equal to nyn, except for restart-array
       INTEGER(iwp) ::  nxl_l  !< west bound of local 2d array start/end_index, is equal to nyn, except for restart-array
       INTEGER(iwp) ::  nxr_l  !< east bound of local 2d array start/end_index, is equal to nyn, except for restart-array

       TYPE(surf_type) ::  surfaces !< respective surface type

       IF ( ALLOCATED( surfaces%start_index ) )                   &
                      CALL deallocate_surface_attributes_v( surfaces )

!
!--    Allocate arrays for start and end index of vertical surface type
!--    for each (j,i)-grid point. This is required in diffion_x, which is
!--    called for each (j,i). In order to find the location where the
!--    respective flux is store within the surface-type, start- and end-
!--    index are stored for each (j,i). For example, each (j,i) can have
!--    several entries where fluxes for vertical surfaces might be stored.
!--    In the flat case, where no vertical walls exit, set indicies such
!--    that loop in diffusion routines will not be entered.
       ALLOCATE ( surfaces%start_index(nys_l:nyn_l,nxl_l:nxr_l) )
       ALLOCATE ( surfaces%end_index(nys_l:nyn_l,nxl_l:nxr_l)   )
       surfaces%start_index = 0
       surfaces%end_index   = -1
!
!--    Indices to locate surface element.
       ALLOCATE ( surfaces%i(1:surfaces%ns) )
       ALLOCATE ( surfaces%j(1:surfaces%ns) )
       ALLOCATE ( surfaces%k(1:surfaces%ns) )
!
!--    Surface-layer height
       ALLOCATE ( surfaces%z_mo(1:surfaces%ns) )
!
!--    Surface orientation
       ALLOCATE ( surfaces%facing(1:surfaces%ns) )
!
!--    Surface parallel wind velocity
       ALLOCATE ( surfaces%uvw_abs(1:surfaces%ns) )
!
!--    Roughness
       ALLOCATE ( surfaces%z0(1:surfaces%ns)  )
       ALLOCATE ( surfaces%z0h(1:surfaces%ns) )
       ALLOCATE ( surfaces%z0q(1:surfaces%ns) )

!
!--    Friction velocity
       ALLOCATE ( surfaces%us(1:surfaces%ns) )
!
!--    Allocate Obukhov length and bulk Richardson number. Actually, at
!--    vertical surfaces these are only required for natural surfaces.
!--    for natural land surfaces
       ALLOCATE( surfaces%ol(1:surfaces%ns)  )
       ALLOCATE( surfaces%rib(1:surfaces%ns) )
!
!--    Allocate arrays for surface momentum fluxes for u and v. For u at north-
!--    and south-facing surfaces, for v at east- and west-facing surfaces.
       ALLOCATE ( surfaces%mom_flux_uv(1:surfaces%ns) )
!
!--    Allocate array for surface momentum flux for w - wsus and wsvs
       ALLOCATE ( surfaces%mom_flux_w(1:surfaces%ns) )
!
!--    Allocate array for surface momentum flux for subgrid-scale tke wsus and
!--    wsvs; first index usvs or vsws, second index for wsus or wsvs, depending
!--    on surface.
       ALLOCATE ( surfaces%mom_flux_tke(0:1,1:surfaces%ns) )
!
!--    Characteristic temperature and surface flux of sensible heat
       ALLOCATE ( surfaces%ts(1:surfaces%ns)  )
       ALLOCATE ( surfaces%shf(1:surfaces%ns) )
!
!--    surface temperature
       ALLOCATE ( surfaces%pt_surface(1:surfaces%ns) )

!--    Arrays for storing potential temperature and
!--    mixing ratio at first grid level
       ALLOCATE ( surfaces%pt1(1:surfaces%ns) )
       ALLOCATE ( surfaces%qv1(1:surfaces%ns) )
       ALLOCATE ( surfaces%sasws(1:surfaces%ns) )
       ALLOCATE ( surfaces%shf_sol(1:surfaces%ns) )

    END SUBROUTINE allocate_surface_attributes_v


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize surface elements, i.e. set initial values for surface fluxes,
!> friction velocity, calcuation of start/end indices, etc. .
!> Please note, further initialization concerning
!> special surface characteristics, e.g. soil- and vegatation type,
!> building type, etc., is done in the land-surface and urban-surface module,
!> respectively.
!------------------------------------------------------------------------------!
    SUBROUTINE init_surfaces

       IMPLICIT NONE

       INTEGER(iwp) ::  i         !< running index x-direction
       INTEGER(iwp) ::  j         !< running index y-direction
       INTEGER(iwp) ::  k         !< running index z-direction
       INTEGER(iwp) ::  l         !< index variable used to distinguish surface facing
       INTEGER(iwp) ::  m         !< running index surface elements

       INTEGER(iwp), DIMENSION(0:2) ::  num_def_h     !< current number of horizontal surface element, default type
       INTEGER(iwp), DIMENSION(0:2) ::  num_def_h_kji !< dummy to determing local end index in surface type for given (j,i), for horizonal default surfaces
       INTEGER(iwp), DIMENSION(0:2) ::  start_index_def_h !< dummy to determing local start index in surface type for given (j,i), for horizontal default surfaces

       INTEGER(iwp), DIMENSION(0:3) ::  num_def_v     !< current number of vertical surface element, default type
       INTEGER(iwp), DIMENSION(0:3) ::  num_def_v_kji !< dummy to determing local end index in surface type for given (j,i), for vertical default surfaces

       INTEGER(iwp), DIMENSION(0:3) ::  start_index_def_v !< dummy to determing local start index in surface type for given (j,i), for vertical default surfaces

       LOGICAL ::  building     !< flag indicating building grid point
       LOGICAL ::  terrain      !< flag indicating natural terrain grid point

!
!--    Set offset indices, i.e. index difference between surface element and
!--    surface-bounded grid point.
!--    Upward facing - no horizontal offsets
       surf_def_h(0:2)%ioff = 0
       surf_def_h(0:2)%joff = 0

!--    Upward facing vertical offsets
       surf_def_h(0)%koff   = -1
!
!--    Downward facing vertical offset
       surf_def_h(1:2)%koff = 1
!
!--    Vertical surfaces - no vertical offset
       surf_def_v(0:3)%koff = 0
!
!--    North- and southward facing - no offset in x
       surf_def_v(0:1)%ioff = 0
!
!--    Northward facing offset in y
       surf_def_v(0)%joff = -1
!
!--    Southward facing offset in y
       surf_def_v(1)%joff = 1
!
!--    East- and westward facing - no offset in y
       surf_def_v(2:3)%joff = 0
!
!--    Eastward facing offset in x
       surf_def_v(2)%ioff = -1
!
!--    Westward facing offset in y
       surf_def_v(3)%ioff = 1
!
!--    Initialize surface attributes, store indicies, surfaces orientation, etc.,
       num_def_h(0:2) = 1
       num_def_v(0:3) = 1

       start_index_def_h(0:2) = 1
       start_index_def_v(0:3) = 1

       DO  i = nxl, nxr
          DO  j = nys, nyn

             num_def_h_kji = 0
             num_def_v_kji = 0

             DO  k = nzb+1, nzt
!
!--             Check if current gridpoint belongs to the atmosphere
                IF ( BTEST( wall_flags_0(k,j,i), 0 ) )  THEN
!
!--                Upward-facing surface. Distinguish between differet surface types.
!--                To do, think about method to flag natural and non-natural
!--                surfaces.
                   IF ( .NOT. BTEST( wall_flags_0(k-1,j,i), 0 ) )  THEN
!
                     CALL initialize_horizontal_surfaces( k, j, i,         &
                                                           surf_def_h(0),   &
                                                            num_def_h(0),    &
                                                            num_def_h_kji(0),&
                                                            .TRUE., .FALSE. )
                   ENDIF
!
!--                downward-facing surface, first, model top. Please note,
!--                for the moment, downward-facing surfaces are always of
!--                default type
                   IF ( k == nzt  .AND.  use_top_fluxes )  THEN
                      CALL initialize_top( k, j, i, surf_def_h(2),             &
                                           num_def_h(2), num_def_h_kji(2) )
!
!--                Check for any other downward-facing surface. So far only for
!--                default surface type.
                   ELSEIF ( .NOT. BTEST( wall_flags_0(k+1,j,i), 0 ) )  THEN
                      CALL initialize_horizontal_surfaces( k, j, i,            &
                                                           surf_def_h(1),      &
                                                           num_def_h(1),       &
                                                           num_def_h_kji(1),   &
                                                           .FALSE., .TRUE. )
                   ENDIF
!
!--                Check for vertical walls and, if required, initialize it.
!                  Start with northward-facing surface.
                   IF ( .NOT. BTEST( wall_flags_0(k,j-1,i), 0 ) )  THEN
                       CALL initialize_vertical_surfaces( 0, k, j, i,        &
                                                            surf_def_v(0),     &
                                                            num_def_v(0),      &
                                                            num_def_v_kji(0),  &
                                                            .FALSE., .FALSE.,  &
                                                            .FALSE., .TRUE. )
                   ENDIF
!
!--                southward-facing surface
                   IF ( .NOT. BTEST( wall_flags_0(k,j+1,i), 0 ) )  THEN
!
                     CALL initialize_vertical_surfaces( 1, k, j, i,        &
                                                            surf_def_v(1),     &
                                                            num_def_v(1),      &
                                                            num_def_v_kji(1),  &
                                                            .FALSE., .FALSE.,  &
                                                            .TRUE., .FALSE. )
                   ENDIF
!
!--                eastward-facing surface
                   IF ( .NOT. BTEST( wall_flags_0(k,j,i-1), 0 ) )  THEN
!
                        CALL initialize_vertical_surfaces( 2, k, j, i,        &
                                                            surf_def_v(2),     &
                                                            num_def_v(2),      &
                                                            num_def_v_kji(2),  &
                                                            .TRUE., .FALSE.,   &
                                                            .FALSE., .FALSE. )
                   ENDIF
!
!--                westward-facing surface
                   IF ( .NOT. BTEST( wall_flags_0(k,j,i+1), 0 ) )  THEN
                        CALL initialize_vertical_surfaces( 3, k, j, i,        &
                                                            surf_def_v(3),     &
                                                            num_def_v(3),      &
                                                            num_def_v_kji(3),  &
                                                           .FALSE., .TRUE.,    &
                                                           .FALSE., .FALSE. )
                   ENDIF
                ENDIF


             ENDDO
!
!--          Determine start- and end-index at grid point (j,i). Also, for
!--          horizontal surfaces more than 1 horizontal surface element can
!--          exist at grid point (j,i) if overhanging structures are present.
!--          Upward-facing surfaces
             surf_def_h(0)%start_index(j,i) = start_index_def_h(0)
             surf_def_h(0)%end_index(j,i)   = surf_def_h(0)%start_index(j,i) + &
                                                 num_def_h_kji(0) - 1
             start_index_def_h(0)           = surf_def_h(0)%end_index(j,i) + 1
!
!--          Downward-facing surfaces, except model top
             surf_def_h(1)%start_index(j,i) = start_index_def_h(1)
             surf_def_h(1)%end_index(j,i)   = surf_def_h(1)%start_index(j,i) + &
                                                 num_def_h_kji(1) - 1
             start_index_def_h(1)           = surf_def_h(1)%end_index(j,i) + 1
!
!--          Downward-facing surfaces -- model top fluxes
             surf_def_h(2)%start_index(j,i) = start_index_def_h(2)
             surf_def_h(2)%end_index(j,i)   = surf_def_h(2)%start_index(j,i) + &
                                                 num_def_h_kji(2) - 1
             start_index_def_h(2)           = surf_def_h(2)%end_index(j,i) + 1
!
!--          Vertical surfaces - Default type
             surf_def_v(0)%start_index(j,i) = start_index_def_v(0)
             surf_def_v(1)%start_index(j,i) = start_index_def_v(1)
             surf_def_v(2)%start_index(j,i) = start_index_def_v(2)
             surf_def_v(3)%start_index(j,i) = start_index_def_v(3)
             surf_def_v(0)%end_index(j,i)   = start_index_def_v(0) +           &
                                              num_def_v_kji(0) - 1
             surf_def_v(1)%end_index(j,i)   = start_index_def_v(1) +           &
                                              num_def_v_kji(1) - 1
             surf_def_v(2)%end_index(j,i)   = start_index_def_v(2) +           &
                                              num_def_v_kji(2) - 1
             surf_def_v(3)%end_index(j,i)   = start_index_def_v(3) +           &
                                              num_def_v_kji(3) - 1
             start_index_def_v(0)           = surf_def_v(0)%end_index(j,i) + 1
             start_index_def_v(1)           = surf_def_v(1)%end_index(j,i) + 1
             start_index_def_v(2)           = surf_def_v(2)%end_index(j,i) + 1
             start_index_def_v(3)           = surf_def_v(3)%end_index(j,i) + 1

          ENDDO
       ENDDO

       CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize horizontal surface elements, upward- and downward-facing.
!> Note, horizontal surface type alsw comprises model-top fluxes, which are,
!> initialized in a different routine.
!------------------------------------------------------------------------------!
          SUBROUTINE initialize_horizontal_surfaces( k, j, i, surf, num_h,     &
                                                     num_h_kji, upward_facing, &
                                                     downward_facing )

             IMPLICIT NONE

             INTEGER(iwp)  ::  i                !< running index x-direction
             INTEGER(iwp)  ::  j                !< running index y-direction
             INTEGER(iwp)  ::  k                !< running index z-direction
             INTEGER(iwp)  ::  num_h            !< current number of surface element
             INTEGER(iwp)  ::  num_h_kji        !< dummy increment

             LOGICAL       ::  upward_facing    !< flag indicating upward-facing surface
             LOGICAL       ::  downward_facing  !< flag indicating downward-facing surface

             TYPE( surf_type ) :: surf          !< respective surface type
!             REAL(WP)      ::  wb_sfc, tod,arg1      !< surface buoyancy forcing -- only matters for ocean
!
!--          Store indices of respective surface element
             surf%i(num_h) = i
             surf%j(num_h) = j
             surf%k(num_h) = k
!
!--          Surface orientation, bit 0 is set to 1 for upward-facing surfaces,
!--          bit 1 is for downward-facing surfaces.
             IF ( upward_facing   )  surf%facing(num_h) = IBSET( surf%facing(num_h), 0 )
             IF ( downward_facing )  surf%facing(num_h) = IBSET( surf%facing(num_h), 1 )
!
!--          Initialize surface-layer height
             IF ( upward_facing )  THEN
                surf%z_mo(num_h)  = zu(k) - zw(k-1)
             ELSE
                surf%z_mo(num_h)  = zw(k) - zu(k)
             ENDIF

             surf%z0(num_h)    = roughness_length
             surf%z0h(num_h)   = z0h_factor * roughness_length
             surf%z0q(num_h)   = z0h_factor * roughness_length

                surf%ol(num_h)   = surf%z_mo(num_h) / zeta_min
!
!--             Very small number is required for calculation of Obukhov length
!--             at first timestep
             surf%us(num_h)    = 1E-30_wp
             surf%usws(num_h)  = 0.0_wp
             surf%vsws(num_h)  = 0.0_wp

             surf%rib(num_h)     = 0.0_wp
             surf%uvw_abs(num_h) = 0.0_wp

             surf%u_0(num_h)     = 0.0_wp
             surf%v_0(num_h)     = 0.0_wp

             surf%ts(num_h)   = 0.0_wp

!--          Set initial value for surface temperature
             surf%pt_surface(num_h) = pt_surface
!
!--          Inititalize surface fluxes of sensible and latent heat, as well as
!--          passive scalar
             IF ( use_surface_fluxes )  THEN

                IF ( upward_facing )  THEN
                   IF ( constant_heatflux )  THEN
!
!--                   Initialize surface heatflux. However, skip this for now if
!--                   if random_heatflux is set. This case, shf is initialized later.
                      IF ( .NOT. random_heatflux )  THEN
                         surf%shf(num_h) = surface_heatflux *                  &
                                                 heatflux_input_conversion(k-1)
!
!--                      Check if surface heat flux might be replaced by
!--                      prescribed wall heatflux
                         IF ( k-1 /= 0 )  THEN
                            surf%shf(num_h) = wall_heatflux(0) *               &
                                                 heatflux_input_conversion(k-1)
                         ENDIF
                      ENDIF
                   ELSE
                      surf%shf(num_h) = 0.0_wp
                   ENDIF
!
!--             Set heat-flux at downward-facing surfaces
                ELSE
                   surf%shf(num_h) = wall_heatflux(5) *                        &
                                             heatflux_input_conversion(k)
                ENDIF
                   IF ( upward_facing )  THEN
                      surf%sasws(num_h) = bottom_salinityflux * rho_air_zw(k-1)
                   ELSE
                      surf%sasws(num_h) = 0.0_wp
                   ENDIF
             ENDIF
!
!--          Increment surface indices
             num_h     = num_h + 1
             num_h_kji = num_h_kji + 1


          END SUBROUTINE initialize_horizontal_surfaces


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize model-top fluxes. Currently, only the heatflux and salinity flux
!> can be prescribed, latent flux is zero in this case!
!------------------------------------------------------------------------------!
          SUBROUTINE initialize_top( k, j, i, surf, num_h, num_h_kji )

             IMPLICIT NONE

             INTEGER(iwp)  ::  i                !< running index x-direction
             INTEGER(iwp)  ::  j                !< running index y-direction
             INTEGER(iwp)  ::  k                !< running index z-direction
             INTEGER(iwp)  ::  num_h            !< current number of surface element
             INTEGER(iwp)  ::  num_h_kji        !< dummy increment
             REAL(WP)      ::  wb_sfc, tod,arg1      !< surface buoyancy forcing -- only matters for ocean

             TYPE( surf_type ) :: surf          !< respective surface type
!
!--          Store indices of respective surface element
             surf%i(num_h) = i
             surf%j(num_h) = j
             surf%k(num_h) = k
!
             IF ( ocean ) THEN
               surf%shf(num_h) = 0.0_wp
               surf%sasws(num_h) = 0.0_wp
               surf%shf_sol(num_h) = 0.0_wp
             endif

!--          Initialize top heat flux
             IF ( constant_top_heatflux )                                      &
                surf%shf(num_h) = top_heatflux * heatflux_input_conversion(nzt+1)
!--          Prescribe top scalar flux

!--          Prescribe top salinity flux
             IF ( ocean .AND. constant_top_salinityflux)                       &
                surf%sasws(num_h) = top_salinityflux * rho_air_zw(nzt+1)

!--          Top momentum fluxes
             IF ( constant_top_momentumflux )  THEN
                surf%usws(num_h) = top_momentumflux_u *                        &
                                            momentumflux_input_conversion(nzt+1)
                surf%vsws(num_h) = top_momentumflux_v *                        &
                                            momentumflux_input_conversion(nzt+1)
             ENDIF
!
!--          Increment surface indices
             num_h     = num_h + 1
             num_h_kji = num_h_kji + 1


          END SUBROUTINE initialize_top


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize vertical surface elements.
!------------------------------------------------------------------------------!
          SUBROUTINE initialize_vertical_surfaces( l, k, j, i, surf, num_v,    &
                                                num_v_kji, east_facing,        &
                                                west_facing, south_facing,     &
                                                north_facing )

             IMPLICIT NONE

             INTEGER(iwp)  ::  component       !< index of wall_fluxes_ array for respective orientation
             INTEGER(iwp)  ::  i               !< running index x-direction
             INTEGER(iwp)  ::  j               !< running index x-direction
             INTEGER(iwp)  ::  k               !< running index x-direction
             INTEGER(iwp)  ::  l               !< index variable for the surface type, indicating the facing
             INTEGER(iwp)  ::  num_v           !< current number of surface element
             INTEGER(iwp)  ::  num_v_kji       !< current number of surface element at (j,i)

             LOGICAL       ::  east_facing     !< flag indicating east-facing surfaces
             LOGICAL       ::  north_facing    !< flag indicating north-facing surfaces
             LOGICAL       ::  south_facing    !< flag indicating south-facing surfaces
             LOGICAL       ::  west_facing     !< flag indicating west-facing surfaces

             TYPE( surf_type ) :: surf         !< respective surface type

!
!--          Store indices of respective wall element
             surf%i(num_v)   = i
             surf%j(num_v)   = j
             surf%k(num_v)   = k
!
!--          Initialize surface-layer height, or more precisely, distance to surface
             IF ( north_facing  .OR.  south_facing )  THEN
                surf%z_mo(num_v)  = 0.5_wp * dy
             ELSE
                surf%z_mo(num_v)  = 0.5_wp * dx
             ENDIF

             surf%facing(num_v)  = 0
!
!--          Surface orientation. Moreover, set component id to map wall_heatflux,
!--          etc., on surface type (further below)
             IF ( north_facing )  THEN
                surf%facing(num_v) = 5 !IBSET( surf%facing(num_v), 0 )
                component          = 4
             ENDIF

             IF ( south_facing )  THEN
                surf%facing(num_v) = 6 !IBSET( surf%facing(num_v), 1 )
                component          = 3
             ENDIF

             IF ( east_facing )  THEN
                surf%facing(num_v) = 7 !IBSET( surf%facing(num_v), 2 )
                component          = 2
             ENDIF

             IF ( west_facing )  THEN
                surf%facing(num_v) = 8 !IBSET( surf%facing(num_v), 3 )
                component          = 1
             ENDIF


             surf%z0(num_v)  = roughness_length
             surf%z0h(num_v) = z0h_factor * roughness_length
             surf%z0q(num_v) = z0h_factor * roughness_length

             surf%us(num_v)  = 0.0_wp
!
!--          If required, initialize Obukhov length
             IF ( ALLOCATED( surf%ol ) )                                       &
                surf%ol(num_v) = surf%z_mo(num_v) / zeta_min

             surf%uvw_abs(num_v)   = 0.0_wp

             surf%mom_flux_uv(num_v) = 0.0_wp
             surf%mom_flux_w(num_v)  = 0.0_wp
             surf%mom_flux_tke(0:1,num_v) = 0.0_wp

             surf%ts(num_v)    = 0.0_wp
             surf%shf(num_v)   = wall_heatflux(component)
!
!--          Set initial value for surface temperature
             surf%pt_surface(num_v) = pt_surface
!
!--          So far, salinityflux at vertical surfaces is simply zero
!--          at the moment
             surf%sasws(num_v) = wall_salinityflux(component)
!
!--          Increment wall indices
             num_v                 = num_v + 1
             num_v_kji             = num_v_kji + 1

          END SUBROUTINE initialize_vertical_surfaces

    END SUBROUTINE init_surfaces


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Determines topography-top index at given (j,i)-position.
!------------------------------------------------------------------------------!
    FUNCTION get_topography_top_index_ji( j, i, grid )

       IMPLICIT NONE

       CHARACTER(LEN=*) ::  grid                         !< flag to distinquish between staggered grids
       INTEGER(iwp)     ::  i                            !< grid index in x-dimension
       INTEGER(iwp)     ::  ibit                         !< bit position where topography information is stored on respective grid
       INTEGER(iwp)     ::  j                            !< grid index in y-dimension
       INTEGER(iwp)     ::  get_topography_top_index_ji  !< topography top index

       SELECT CASE ( TRIM( grid ) )

          CASE ( 's'     )
             ibit = 12
          CASE ( 'u'     )
             ibit = 14
          CASE ( 'v'     )
             ibit = 16
          CASE ( 'w'     )
             ibit = 18
          CASE ( 's_out' )
             ibit = 24
          CASE DEFAULT
!
!--          Set default to scalar grid
             ibit = 12

       END SELECT

       get_topography_top_index_ji = MAXLOC(                                   &
                                     MERGE( 1, 0,                              &
                                            BTEST( wall_flags_0(:,j,i), ibit ) &
                                          ), DIM = 1                           &
                                           ) - 1

       RETURN

    END FUNCTION get_topography_top_index_ji

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Determines topography-top index at each (j,i)-position.
!------------------------------------------------------------------------------!
    FUNCTION get_topography_top_index( grid )

       IMPLICIT NONE

       CHARACTER(LEN=*) ::  grid                      !< flag to distinquish between staggered grids
       INTEGER(iwp)     ::  ibit                      !< bit position where topography information is stored on respective grid
       INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  get_topography_top_index  !< topography top index

       SELECT CASE ( TRIM( grid ) )

          CASE ( 's'     )
             ibit = 12
          CASE ( 'u'     )
             ibit = 14
          CASE ( 'v'     )
             ibit = 16
          CASE ( 'w'     )
             ibit = 18
          CASE ( 's_out' )
             ibit = 24
          CASE DEFAULT
!
!--          Set default to scalar grid
             ibit = 12

       END SELECT

       get_topography_top_index(nys:nyn,nxl:nxr) = MAXLOC(                     &
                         MERGE( 1, 0,                                          &
                                 BTEST( wall_flags_0(:,nys:nyn,nxl:nxr), ibit )&
                              ), DIM = 1                                       &
                                                         ) - 1

       RETURN

    END FUNCTION get_topography_top_index

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Gathers all surface elements with the same facing (but possibly different
!> type) onto a surface type, and writes binary data into restart files.
!------------------------------------------------------------------------------!
    SUBROUTINE surface_wrd_local


       IMPLICIT NONE

       CHARACTER(LEN=1)             ::  dum           !< dummy string to create output-variable name

       INTEGER(iwp)                 ::  i             !< running index x-direction
       INTEGER(iwp)                 ::  j             !< running index y-direction
       INTEGER(iwp)                 ::  l             !< index surface type orientation
       INTEGER(iwp)                 ::  m             !< running index for surface elements on individual surface array
       INTEGER(iwp), DIMENSION(0:2) ::  start_index_h !< start index for horizontal surface elements on gathered surface array
       INTEGER(iwp), DIMENSION(0:3) ::  mm            !< running index for surface elements on gathered surface array
       INTEGER(iwp), DIMENSION(0:3) ::  start_index_v !< start index for vertical surface elements on gathered surface array

       TYPE(surf_type), DIMENSION(0:2) ::  surf_h     !< gathered horizontal surfaces, contains all surface types
       TYPE(surf_type), DIMENSION(0:3) ::  surf_v     !< gathered vertical surfaces, contains all surface types

!
!--    Determine total number of horizontal and vertical surface elements before
!--    writing var_list
       CALL surface_last_actions
!
!--    Count number of grid points with same facing and allocate attributes respectively
!--    Horizontal upward facing
       surf_h(0)%ns = ns_h_on_file(0)
       CALL allocate_surface_attributes_h( surf_h(0), nys, nyn, nxl, nxr )
!
!--    Horizontal downward facing
       surf_h(1)%ns = ns_h_on_file(1)
       CALL allocate_surface_attributes_h( surf_h(1), nys, nyn, nxl, nxr )
!
!--    Model top
       surf_h(2)%ns = ns_h_on_file(2)
       CALL allocate_surface_attributes_h_top( surf_h(2), nys, nyn, nxl, nxr )
!
!--    Vertical surfaces
       DO  l = 0, 3
          surf_v(l)%ns = ns_v_on_file(l)
          CALL allocate_surface_attributes_v( surf_v(l),                       &
                                              nys, nyn, nxl, nxr )
       ENDDO
!
!--    In the following, gather data from surfaces elements with the same
!--    facing (but possibly differt type) on 1 data-type array.
       mm(0:2) = 1
       DO  l = 0, 2
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  m = surf_def_h(l)%start_index(j,i),                        &
                        surf_def_h(l)%end_index(j,i)
                   IF ( ALLOCATED( surf_def_h(l)%us ) )                        &
                      surf_h(l)%us(mm(l))      = surf_def_h(l)%us(m)
                   IF ( ALLOCATED( surf_def_h(l)%ts ) )                        &
                      surf_h(l)%ts(mm(l))      = surf_def_h(l)%ts(m)
                   IF ( ALLOCATED( surf_def_h(l)%qs ) )                        &
                      surf_h(l)%qs(mm(l))      = surf_def_h(l)%qs(m)
                   IF ( ALLOCATED( surf_def_h(l)%ss ) )                        &
                      surf_h(l)%ss(mm(l))      = surf_def_h(l)%ss(m)
                   IF ( ALLOCATED( surf_def_h(l)%qcs ) )                       &
                      surf_h(l)%qcs(mm(l))     = surf_def_h(l)%qcs(m)
                   IF ( ALLOCATED( surf_def_h(l)%ncs ) )                       &
                      surf_h(l)%ncs(mm(l))     = surf_def_h(l)%ncs(m)
                   IF ( ALLOCATED( surf_def_h(l)%qrs ) )                       &
                      surf_h(l)%qrs(mm(l))     = surf_def_h(l)%qrs(m)
                   IF ( ALLOCATED( surf_def_h(l)%nrs ) )                       &
                      surf_h(l)%nrs(mm(l))     = surf_def_h(l)%nrs(m)
                   IF ( ALLOCATED( surf_def_h(l)%ol ) )                        &
                      surf_h(l)%ol(mm(l))      = surf_def_h(l)%ol(m)
                   IF ( ALLOCATED( surf_def_h(l)%rib ) )                       &
                      surf_h(l)%rib(mm(l))     = surf_def_h(l)%rib(m)
                   IF ( ALLOCATED( surf_def_h(l)%pt_surface ) )                &
                      surf_h(l)%pt_surface(mm(l)) = surf_def_h(l)%pt_surface(m)
                   IF ( ALLOCATED( surf_def_h(l)%usws ) )                      &
                      surf_h(l)%usws(mm(l))    = surf_def_h(l)%usws(m)
                   IF ( ALLOCATED( surf_def_h(l)%vsws ) )                      &
                      surf_h(l)%vsws(mm(l))    = surf_def_h(l)%vsws(m)
                   IF ( ALLOCATED( surf_def_h(l)%shf ) )                       &
                      surf_h(l)%shf(mm(l))     = surf_def_h(l)%shf(m)
                   IF ( ALLOCATED( surf_def_h(l)%qsws ) )                      &
                      surf_h(l)%qsws(mm(l))    = surf_def_h(l)%qsws(m)
                   IF ( ALLOCATED( surf_def_h(l)%ssws ) )                      &
                      surf_h(l)%ssws(mm(l))    = surf_def_h(l)%ssws(m)
                   IF ( ALLOCATED( surf_def_h(l)%ncsws ) )                     &
                      surf_h(l)%ncsws(mm(l))   = surf_def_h(l)%ncsws(m)
                   IF ( ALLOCATED( surf_def_h(l)%nrsws ) )                     &
                      surf_h(l)%nrsws(mm(l))   = surf_def_h(l)%nrsws(m)
                   IF ( ALLOCATED( surf_def_h(l)%sasws ) )                     &
                      surf_h(l)%sasws(mm(l))   = surf_def_h(l)%sasws(m)

                   mm(l) = mm(l) + 1
                ENDDO

             ENDDO

          ENDDO
!
!--       Gather start- and end indices
          start_index_h(l) = 1
          DO  i = nxl, nxr
             DO  j = nys, nyn

                surf_h(l)%start_index(j,i) = start_index_h(l)
                surf_h(l)%end_index(j,i)   = surf_h(l)%start_index(j,i) -1

                DO  m = surf_def_h(l)%start_index(j,i),                        &
                        surf_def_h(l)%end_index(j,i)
                   surf_h(l)%end_index(j,i) = surf_h(l)%end_index(j,i) + 1
                ENDDO
                start_index_h(l) = surf_h(l)%end_index(j,i) + 1

             ENDDO
          ENDDO
       ENDDO


       mm(0:3) = 1
       DO  l = 0, 3
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  m = surf_def_v(l)%start_index(j,i),                        &
                        surf_def_v(l)%end_index(j,i)
                   IF ( ALLOCATED( surf_def_v(l)%us ) )                        &
                      surf_v(l)%us(mm(l))      = surf_def_v(l)%us(m)
                   IF ( ALLOCATED( surf_def_v(l)%ts ) )                        &
                      surf_v(l)%ts(mm(l))      = surf_def_v(l)%ts(m)
                   IF ( ALLOCATED( surf_def_v(l)%qs ) )                        &
                      surf_v(l)%qs(mm(l))      = surf_def_v(l)%qs(m)
                   IF ( ALLOCATED( surf_def_v(l)%ss ) )                        &
                      surf_v(l)%ss(mm(l))      = surf_def_v(l)%ss(m)
                   IF ( ALLOCATED( surf_def_v(l)%qcs ) )                       &
                      surf_v(l)%qcs(mm(l))     = surf_def_v(l)%qcs(m)
                   IF ( ALLOCATED( surf_def_v(l)%ncs ) )                       &
                      surf_v(l)%ncs(mm(l))     = surf_def_v(l)%ncs(m)
                   IF ( ALLOCATED( surf_def_v(l)%qrs ) )                       &
                      surf_v(l)%qrs(mm(l))     = surf_def_v(l)%qrs(m)
                   IF ( ALLOCATED( surf_def_v(l)%nrs ) )                       &
                      surf_v(l)%nrs(mm(l))     = surf_def_v(l)%nrs(m)
                   IF ( ALLOCATED( surf_def_v(l)%ol ) )                        &
                      surf_v(l)%ol(mm(l))      = surf_def_v(l)%ol(m)
                   IF ( ALLOCATED( surf_def_v(l)%rib ) )                       &
                      surf_v(l)%rib(mm(l))     = surf_def_v(l)%rib(m)
                   IF ( ALLOCATED( surf_def_v(l)%pt_surface ) )                &
                      surf_v(l)%pt_surface(mm(l)) = surf_def_v(l)%pt_surface(m)
                   IF ( ALLOCATED( surf_def_v(l)%shf ) )                       &
                      surf_v(l)%shf(mm(l))     = surf_def_v(l)%shf(m)
                   IF ( ALLOCATED( surf_def_v(l)%qsws ) )                      &
                      surf_v(l)%qsws(mm(l))    = surf_def_v(l)%qsws(m)
                   IF ( ALLOCATED( surf_def_v(l)%ssws ) )                      &
                      surf_v(l)%ssws(mm(l))    = surf_def_v(l)%ssws(m)
                   IF ( ALLOCATED( surf_def_v(l)%ncsws ) )                     &
                      surf_v(l)%ncsws(mm(l))   = surf_def_v(l)%ncsws(m)
                   IF ( ALLOCATED( surf_def_v(l)%nrsws ) )                     &
                      surf_v(l)%nrsws(mm(l))   = surf_def_v(l)%nrsws(m)
                   IF ( ALLOCATED( surf_def_v(l)%sasws ) )                     &
                      surf_v(l)%sasws(mm(l))   = surf_def_v(l)%sasws(m)
                   IF ( ALLOCATED( surf_def_v(l)%mom_flux_uv) )                &
                      surf_v(l)%mom_flux_uv(mm(l))  = surf_def_v(l)%mom_flux_uv(m)
                   IF ( ALLOCATED( surf_def_v(l)%mom_flux_w) )                 &
                      surf_v(l)%mom_flux_w(mm(l))   = surf_def_v(l)%mom_flux_w(m)
                   IF ( ALLOCATED( surf_def_v(l)%mom_flux_tke) )               &
                      surf_v(l)%mom_flux_tke(0:1,mm(l)) = surf_def_v(l)%mom_flux_tke(0:1,m)

                   mm(l) = mm(l) + 1
                ENDDO
             ENDDO
          ENDDO
!
!--       Gather start- and end indices
          start_index_v(l) = 1
          DO  i = nxl, nxr
             DO  j = nys, nyn

                surf_v(l)%start_index(j,i) = start_index_v(l)
                surf_v(l)%end_index(j,i)   = surf_v(l)%start_index(j,i) -1

                DO  m = surf_def_v(l)%start_index(j,i),                        &
                        surf_def_v(l)%end_index(j,i)
                   surf_v(l)%end_index(j,i) = surf_v(l)%end_index(j,i) + 1
                ENDDO

                start_index_v(l) = surf_v(l)%end_index(j,i) + 1
             ENDDO
          ENDDO

       ENDDO


       CALL wrd_write_string( 'ns_h_on_file' )
       WRITE ( 14 )  ns_h_on_file

       CALL wrd_write_string( 'ns_v_on_file' )
       WRITE ( 14 )  ns_v_on_file

!
!--    Write required restart data.
!--    Start with horizontal surfaces (upward-, downward-facing, and model top)
       DO  l = 0, 2
          WRITE( dum, '(I1)')  l

          CALL wrd_write_string( 'surf_h(' // dum // ')%start_index' )
          WRITE ( 14 )  surf_h(l)%start_index

          CALL wrd_write_string( 'surf_h(' // dum // ')%end_index' )
          WRITE ( 14 )  surf_h(l)%end_index

          IF ( ALLOCATED ( surf_h(l)%us ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%us' )
             WRITE ( 14 )  surf_h(l)%us
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%ts ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%ts' )
             WRITE ( 14 )  surf_h(l)%ts
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%qs ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%qs' )
             WRITE ( 14 )  surf_h(l)%qs
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%ss ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%ss' )
             WRITE ( 14 )  surf_h(l)%ss
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%qcs ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%qcs' )
             WRITE ( 14 )  surf_h(l)%qcs
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%ncs ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%ncs' )
             WRITE ( 14 )  surf_h(l)%ncs
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%qrs ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%qrs' )
             WRITE ( 14 )  surf_h(l)%qrs
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%nrs ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%nrs' )
             WRITE ( 14 )  surf_h(l)%nrs
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%ol ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%ol' )
             WRITE ( 14 )  surf_h(l)%ol
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%rib ) )  THEN
            CALL wrd_write_string( 'surf_h(' // dum // ')%rib' )
             WRITE ( 14 )  surf_h(l)%rib
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%pt_surface ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%pt_surface' )
             WRITE ( 14 )  surf_h(l)%pt_surface
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%usws ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%usws' )
             WRITE ( 14 )  surf_h(l)%usws
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%vsws ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%vsws' )
             WRITE ( 14 )  surf_h(l)%vsws
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%shf ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%shf' )
             WRITE ( 14 )  surf_h(l)%shf
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%qsws ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%qsws' )
             WRITE ( 14 )  surf_h(l)%qsws
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%ssws ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%ssws' )
             WRITE ( 14 )  surf_h(l)%ssws
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%css ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%css' )
             WRITE ( 14 )  surf_h(l)%css
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%cssws ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%cssws' )
             WRITE ( 14 )  surf_h(l)%cssws
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%qcsws ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%qcsws' )
             WRITE ( 14 )  surf_h(l)%qcsws
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%ncsws ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%ncsws' )
             WRITE ( 14 )  surf_h(l)%ncsws
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%qrsws ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%qrsws' )
             WRITE ( 14 )  surf_h(l)%qrsws
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%nrsws ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%nrsws' )
             WRITE ( 14 )  surf_h(l)%nrsws
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%sasws ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%sasws' )
             WRITE ( 14 )  surf_h(l)%sasws
          ENDIF

       ENDDO
!
!--    Write vertical surfaces
       DO  l = 0, 3
          WRITE( dum, '(I1)')  l

          CALL wrd_write_string( 'surf_v(' // dum // ')%start_index' )
          WRITE ( 14 )  surf_v(l)%start_index

          CALL wrd_write_string( 'surf_v(' // dum // ')%end_index' )
          WRITE ( 14 )   surf_v(l)%end_index

          IF ( ALLOCATED ( surf_v(l)%us ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%us' )
             WRITE ( 14 )  surf_v(l)%us
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%ts ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%ts' )
             WRITE ( 14 )  surf_v(l)%ts
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%qs ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%qs' )
             WRITE ( 14 )  surf_v(l)%qs
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%ss ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%ss' )
             WRITE ( 14 )  surf_v(l)%ss
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%qcs ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%qcs' )
             WRITE ( 14 )  surf_v(l)%qcs
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%ncs ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%ncs' )
             WRITE ( 14 )  surf_v(l)%ncs
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%qrs ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%qrs' )
             WRITE ( 14 )  surf_v(l)%qrs
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%nrs ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%nrs' )
             WRITE ( 14 )  surf_v(l)%nrs
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%ol ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%ol' )
             WRITE ( 14 )  surf_v(l)%ol
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%rib ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%rib' )
             WRITE ( 14 )  surf_v(l)%rib
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%pt_surface ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%pt_surface' )
             WRITE ( 14 )  surf_v(l)%pt_surface
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%shf ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%shf' )
             WRITE ( 14 )  surf_v(l)%shf
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%qsws ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%qsws' )
             WRITE ( 14 )  surf_v(l)%qsws
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%ssws ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%ssws' )
             WRITE ( 14 )  surf_v(l)%ssws
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%css ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%css' )
             WRITE ( 14 )  surf_v(l)%css
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%cssws ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%cssws' )
             WRITE ( 14 )  surf_v(l)%cssws
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%qcsws ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%qcsws' )
             WRITE ( 14 )  surf_v(l)%qcsws
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%ncsws ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%ncsws' )
             WRITE ( 14 )  surf_v(l)%ncsws
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%qrsws ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%qrsws' )
             WRITE ( 14 )  surf_v(l)%qrsws
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%nrsws ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%nrsws' )
             WRITE ( 14 )  surf_v(l)%nrsws
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%sasws ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%sasws' )
             WRITE ( 14 )  surf_v(l)%sasws
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%mom_flux_uv ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%mom_uv' )
             WRITE ( 14 )  surf_v(l)%mom_flux_uv
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%mom_flux_w ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%mom_w' )
             WRITE ( 14 )  surf_v(l)%mom_flux_w
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%mom_flux_tke ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%mom_tke' )
             WRITE ( 14 )  surf_v(l)%mom_flux_tke
          ENDIF

       ENDDO


    END SUBROUTINE surface_wrd_local


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads surface-related restart data. Please note, restart data for a certain
!> surface orientation (e.g. horizontal upward-facing) is stored in one
!> array, even if surface elements may belong to different surface types
!> natural or urban for example). Surface elements are redistributed into its
!> respective surface types within this routine. This allows e.g. changing the
!> surface type after reading the restart data, which might be required in case
!> of cyclic_fill mode.
!------------------------------------------------------------------------------!
    SUBROUTINE surface_rrd_local( ii, kk, nxlf, nxlc, nxl_on_file, nxrf, nxrc, &
                                  nxr_on_file, nynf, nync, nyn_on_file, nysf,  &
                                  nysc, nys_on_file, found )

       USE pegrid,                                                             &
           ONLY: myid, numprocs_previous_run

       IMPLICIT NONE

       CHARACTER (LEN=1)  ::  dum         !< dummy to create correct string for reading input variable

       INTEGER(iwp)       ::  i           !< running index along x-direction, refers to former domain size
       INTEGER(iwp)       ::  ic          !< running index along x-direction, refers to current domain size
       INTEGER(iwp)       ::  j           !< running index along y-direction, refers to former domain size
       INTEGER(iwp)       ::  jc          !< running index along y-direction, refers to former domain size
       INTEGER(iwp)       ::  k           !< running index along z-direction
       INTEGER(iwp)       ::  m           !< running index for surface elements, refers to gathered array encompassing all surface types
       INTEGER(iwp)       ::  mm          !< running index for surface elements, refers to individual surface types
       INTEGER(iwp)       ::  ii               !< running index over input files
       INTEGER(iwp)       ::  kk               !< running index over previous input files covering current local domain
       INTEGER(iwp)       ::  nxlc             !< index of left boundary on current subdomain
       INTEGER(iwp)       ::  nxlf             !< index of left boundary on former subdomain
       INTEGER(iwp)       ::  nxl_on_file      !< index of left boundary on former local domain
       INTEGER(iwp)       ::  nxrc             !< index of right boundary on current subdomain
       INTEGER(iwp)       ::  nxrf             !< index of right boundary on former subdomain
       INTEGER(iwp)       ::  nxr_on_file      !< index of right boundary on former local domain
       INTEGER(iwp)       ::  nync             !< index of north boundary on current subdomain
       INTEGER(iwp)       ::  nynf             !< index of north boundary on former subdomain
       INTEGER(iwp)       ::  nyn_on_file      !< index of norht boundary on former local domain
       INTEGER(iwp)       ::  nysc             !< index of south boundary on current subdomain
       INTEGER(iwp)       ::  nysf             !< index of south boundary on former subdomain
       INTEGER(iwp)       ::  nys_on_file      !< index of south boundary on former local domain

       INTEGER(iwp), SAVE  ::  l           !< index variable for surface type

       LOGICAL                         ::  surf_match_def     !< flag indicating that surface element is of default type

       LOGICAL, INTENT(OUT) ::  found

       LOGICAL, SAVE ::  horizontal_surface !< flag indicating horizontal surfaces
       LOGICAL, SAVE ::  vertical_surface   !< flag indicating vertical surfaces

       TYPE(surf_type), DIMENSION(0:2), SAVE ::  surf_h             !< horizontal surface type on file
       TYPE(surf_type), DIMENSION(0:3), SAVE ::  surf_v             !< vertical surface type on file


       found              = .TRUE.

       SELECT CASE ( restart_string(1:length) )

          CASE ( 'ns_h_on_file' )
             IF ( kk == 1 )  THEN
                READ ( 13 )  ns_h_on_file

                IF ( ALLOCATED( surf_h(0)%start_index ) )                      &
                   CALL deallocate_surface_attributes_h( surf_h(0) )
                IF ( ALLOCATED( surf_h(1)%start_index ) )                      &
                   CALL deallocate_surface_attributes_h( surf_h(1) )
                IF ( ALLOCATED( surf_h(2)%start_index ) )                      &
                   CALL deallocate_surface_attributes_h_top( surf_h(2) )

!--             Allocate memory for number of surface elements on file.
!--             Please note, these number is not necessarily the same as
!--             the final number of surface elements on local domain,
!--             which is the case if processor topology changes during
!--             restart runs.
!--             Horizontal upward facing
                surf_h(0)%ns = ns_h_on_file(0)
                CALL allocate_surface_attributes_h( surf_h(0),                 &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file, nxr_on_file )

!--             Horizontal downward facing
                surf_h(1)%ns = ns_h_on_file(1)
                CALL allocate_surface_attributes_h( surf_h(1),                 &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file, nxr_on_file )

!--             Model top
                surf_h(2)%ns = ns_h_on_file(2)
                CALL allocate_surface_attributes_h_top( surf_h(2),             &
                                            nys_on_file, nyn_on_file,          &
                                            nxl_on_file, nxr_on_file )

!
!--             Initial setting of flags for horizontal and vertical surfaces,
!--             will be set after start- and end-indices are read.
                horizontal_surface = .FALSE.
                vertical_surface   = .FALSE.

             ENDIF

          CASE ( 'ns_v_on_file' )
             IF ( kk == 1 ) THEN
                READ ( 13 )  ns_v_on_file

                DO  l = 0, 3
                   IF ( ALLOCATED( surf_v(l)%start_index ) )                   &
                      CALL deallocate_surface_attributes_v( surf_v(l) )
                ENDDO

!--                Vertical surfaces
                DO  l = 0, 3
                   surf_v(l)%ns = ns_v_on_file(l)
                   CALL allocate_surface_attributes_v( surf_v(l),              &
                                           nys_on_file, nyn_on_file,           &
                                           nxl_on_file, nxr_on_file )
               ENDDO

             ENDIF

          CASE ( 'surf_h(0)%start_index' )
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_h(0)%start_index
             l = 0
          CASE ( 'surf_h(0)%end_index' )
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_h(0)%end_index
             horizontal_surface = .TRUE.
             vertical_surface   = .FALSE.
          CASE ( 'surf_h(0)%us' )
             IF ( ALLOCATED( surf_h(0)%us )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(0)%us
          CASE ( 'surf_h(0)%ts' )
             IF ( ALLOCATED( surf_h(0)%ts )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(0)%ts
          CASE ( 'surf_h(0)%qs' )
             IF ( ALLOCATED( surf_h(0)%qs )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(0)%qs
          CASE ( 'surf_h(0)%ss' )
             IF ( ALLOCATED( surf_h(0)%ss )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(0)%ss
          CASE ( 'surf_h(0)%qcs' )
             IF ( ALLOCATED( surf_h(0)%qcs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(0)%qcs
          CASE ( 'surf_h(0)%ncs' )
             IF ( ALLOCATED( surf_h(0)%ncs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(0)%ncs
          CASE ( 'surf_h(0)%qrs' )
             IF ( ALLOCATED( surf_h(0)%qrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(0)%qrs
          CASE ( 'surf_h(0)%nrs' )
             IF ( ALLOCATED( surf_h(0)%nrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(0)%nrs
          CASE ( 'surf_h(0)%ol' )
             IF ( ALLOCATED( surf_h(0)%ol )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(0)%ol
          CASE ( 'surf_h(0)%rib' )
             IF ( ALLOCATED( surf_h(0)%rib )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(0)%rib
          CASE ( 'surf_h(0)%pt_surface' )
             IF ( ALLOCATED( surf_h(0)%pt_surface )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_h(0)%pt_surface
          CASE ( 'surf_h(0)%usws' )
             IF ( ALLOCATED( surf_h(0)%usws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(0)%usws
          CASE ( 'surf_h(0)%vsws' )
             IF ( ALLOCATED( surf_h(0)%vsws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(0)%vsws
          CASE ( 'surf_h(0)%shf' )
             IF ( ALLOCATED( surf_h(0)%shf )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(0)%shf
          CASE ( 'surf_h(0)%qsws' )
             IF ( ALLOCATED( surf_h(0)%qsws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(0)%qsws
          CASE ( 'surf_h(0)%ssws' )
             IF ( ALLOCATED( surf_h(0)%ssws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(0)%ssws
          CASE ( 'surf_h(0)%css' )
             IF ( ALLOCATED( surf_h(0)%css )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(0)%css
          CASE ( 'surf_h(0)%cssws' )
             IF ( ALLOCATED( surf_h(0)%cssws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(0)%cssws
          CASE ( 'surf_h(0)%qcsws' )
             IF ( ALLOCATED( surf_h(0)%qcsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(0)%qcsws
          CASE ( 'surf_h(0)%ncsws' )
             IF ( ALLOCATED( surf_h(0)%ncsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(0)%ncsws
          CASE ( 'surf_h(0)%qrsws' )
             IF ( ALLOCATED( surf_h(0)%qrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(0)%qrsws
          CASE ( 'surf_h(0)%nrsws' )
             IF ( ALLOCATED( surf_h(0)%nrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(0)%nrsws
          CASE ( 'surf_h(0)%sasws' )
             IF ( ALLOCATED( surf_h(0)%sasws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(0)%sasws

          CASE ( 'surf_h(1)%start_index' )
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_h(1)%start_index
             l = 1
          CASE ( 'surf_h(1)%end_index' )
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_h(1)%end_index
          CASE ( 'surf_h(1)%us' )
             IF ( ALLOCATED( surf_h(1)%us )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(1)%us
          CASE ( 'surf_h(1)%ts' )
             IF ( ALLOCATED( surf_h(1)%ts )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(1)%ts
          CASE ( 'surf_h(1)%qs' )
             IF ( ALLOCATED( surf_h(1)%qs )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(1)%qs
          CASE ( 'surf_h(1)%ss' )
             IF ( ALLOCATED( surf_h(1)%ss )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(1)%ss
          CASE ( 'surf_h(1)%qcs' )
             IF ( ALLOCATED( surf_h(1)%qcs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(1)%qcs
          CASE ( 'surf_h(1)%ncs' )
             IF ( ALLOCATED( surf_h(1)%ncs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(1)%ncs
          CASE ( 'surf_h(1)%qrs' )
             IF ( ALLOCATED( surf_h(1)%qrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(1)%qrs
          CASE ( 'surf_h(1)%nrs' )
             IF ( ALLOCATED( surf_h(1)%nrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(1)%nrs
          CASE ( 'surf_h(1)%ol' )
             IF ( ALLOCATED( surf_h(1)%ol )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(1)%ol
          CASE ( 'surf_h(1)%rib' )
             IF ( ALLOCATED( surf_h(1)%rib )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(1)%rib
          CASE ( 'surf_h(1)%pt_surface' )
             IF ( ALLOCATED( surf_h(1)%pt_surface )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_h(1)%pt_surface
          CASE ( 'surf_h(1)%usws' )
             IF ( ALLOCATED( surf_h(1)%usws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(1)%usws
          CASE ( 'surf_h(1)%vsws' )
             IF ( ALLOCATED( surf_h(1)%vsws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(1)%vsws
          CASE ( 'surf_h(1)%shf' )
             IF ( ALLOCATED( surf_h(1)%shf )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(1)%shf
          CASE ( 'surf_h(1)%qsws' )
             IF ( ALLOCATED( surf_h(1)%qsws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(1)%qsws
          CASE ( 'surf_h(1)%ssws' )
             IF ( ALLOCATED( surf_h(1)%ssws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(1)%ssws
          CASE ( 'surf_h(1)%css' )
             IF ( ALLOCATED( surf_h(1)%css )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(1)%css
          CASE ( 'surf_h(1)%cssws' )
             IF ( ALLOCATED( surf_h(1)%cssws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(1)%cssws
          CASE ( 'surf_h(1)%qcsws' )
             IF ( ALLOCATED( surf_h(1)%qcsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(1)%qcsws
          CASE ( 'surf_h(1)%ncsws' )
             IF ( ALLOCATED( surf_h(1)%ncsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(1)%ncsws
          CASE ( 'surf_h(1)%qrsws' )
             IF ( ALLOCATED( surf_h(1)%qrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(1)%qrsws
          CASE ( 'surf_h(1)%nrsws' )
             IF ( ALLOCATED( surf_h(1)%nrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(1)%nrsws
          CASE ( 'surf_h(1)%sasws' )
             IF ( ALLOCATED( surf_h(1)%sasws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(1)%sasws

          CASE ( 'surf_h(2)%start_index' )
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_h(2)%start_index
             l = 2
          CASE ( 'surf_h(2)%end_index' )
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_h(2)%end_index
          CASE ( 'surf_h(2)%us' )
             IF ( ALLOCATED( surf_h(2)%us )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(2)%us
          CASE ( 'surf_h(2)%ts' )
             IF ( ALLOCATED( surf_h(2)%ts )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(2)%ts
          CASE ( 'surf_h(2)%qs' )
             IF ( ALLOCATED( surf_h(2)%qs )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(2)%qs
          CASE ( 'surf_h(2)%ss' )
             IF ( ALLOCATED( surf_h(2)%ss )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(2)%ss
          CASE ( 'surf_h(2)%qcs' )
             IF ( ALLOCATED( surf_h(2)%qcs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(2)%qcs
          CASE ( 'surf_h(2)%ncs' )
             IF ( ALLOCATED( surf_h(2)%ncs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(2)%ncs
          CASE ( 'surf_h(2)%qrs' )
             IF ( ALLOCATED( surf_h(2)%qrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(2)%qrs
          CASE ( 'surf_h(2)%nrs' )
             IF ( ALLOCATED( surf_h(2)%nrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(2)%nrs
          CASE ( 'surf_h(2)%ol' )
             IF ( ALLOCATED( surf_h(2)%ol )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(2)%ol
          CASE ( 'surf_h(2)%rib' )
             IF ( ALLOCATED( surf_h(2)%rib )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(2)%rib
          CASE ( 'surf_h(2)%pt_surface' )
             IF ( ALLOCATED( surf_h(2)%pt_surface )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_h(2)%pt_surface
          CASE ( 'surf_h(2)%usws' )
             IF ( ALLOCATED( surf_h(2)%usws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(2)%usws
          CASE ( 'surf_h(2)%vsws' )
             IF ( ALLOCATED( surf_h(2)%vsws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(2)%vsws
          CASE ( 'surf_h(2)%shf' )
             IF ( ALLOCATED( surf_h(2)%shf )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(2)%shf
          CASE ( 'surf_h(2)%qsws' )
             IF ( ALLOCATED( surf_h(2)%qsws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(2)%qsws
          CASE ( 'surf_h(2)%ssws' )
             IF ( ALLOCATED( surf_h(2)%ssws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(2)%ssws
          CASE ( 'surf_h(2)%css' )
             IF ( ALLOCATED( surf_h(2)%css )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(2)%css
          CASE ( 'surf_h(2)%cssws' )
             IF ( ALLOCATED( surf_h(2)%cssws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(2)%cssws
          CASE ( 'surf_h(2)%qcsws' )
             IF ( ALLOCATED( surf_h(2)%qcsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(2)%qcsws
          CASE ( 'surf_h(2)%ncsws' )
             IF ( ALLOCATED( surf_h(2)%ncsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(2)%ncsws
          CASE ( 'surf_h(2)%qrsws' )
             IF ( ALLOCATED( surf_h(2)%qrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(2)%qrsws
          CASE ( 'surf_h(2)%nrsws' )
             IF ( ALLOCATED( surf_h(2)%nrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(2)%nrsws
          CASE ( 'surf_h(2)%sasws' )
             IF ( ALLOCATED( surf_h(2)%sasws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(2)%sasws

          CASE ( 'surf_v(0)%start_index' )
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_v(0)%start_index
             l = 0
             horizontal_surface = .FALSE.
             vertical_surface   = .TRUE.
          CASE ( 'surf_v(0)%end_index' )
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_v(0)%end_index
          CASE ( 'surf_v(0)%us' )
             IF ( ALLOCATED( surf_v(0)%us )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(0)%us
          CASE ( 'surf_v(0)%ts' )
             IF ( ALLOCATED( surf_v(0)%ts )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(0)%ts
          CASE ( 'surf_v(0)%qs' )
             IF ( ALLOCATED( surf_v(0)%qs )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(0)%qs
          CASE ( 'surf_v(0)%ss' )
             IF ( ALLOCATED( surf_v(0)%ss )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(0)%ss
          CASE ( 'surf_v(0)%qcs' )
             IF ( ALLOCATED( surf_v(0)%qcs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(0)%qcs
          CASE ( 'surf_v(0)%ncs' )
             IF ( ALLOCATED( surf_v(0)%ncs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(0)%ncs
          CASE ( 'surf_v(0)%qrs' )
             IF ( ALLOCATED( surf_v(0)%qrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(0)%qrs
          CASE ( 'surf_v(0)%nrs' )
             IF ( ALLOCATED( surf_v(0)%nrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(0)%nrs
          CASE ( 'surf_v(0)%ol' )
             IF ( ALLOCATED( surf_v(0)%ol )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(0)%ol
          CASE ( 'surf_v(0)%rib' )
             IF ( ALLOCATED( surf_v(0)%rib )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(0)%rib
          CASE ( 'surf_v(0)%pt_surface' )
             IF ( ALLOCATED( surf_v(0)%pt_surface )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_v(0)%pt_surface
          CASE ( 'surf_v(0)%shf' )
             IF ( ALLOCATED( surf_v(0)%shf )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(0)%shf
          CASE ( 'surf_v(0)%qsws' )
             IF ( ALLOCATED( surf_v(0)%qsws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_v(0)%qsws
          CASE ( 'surf_v(0)%ssws' )
             IF ( ALLOCATED( surf_v(0)%ssws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_v(0)%ssws
          CASE ( 'surf_v(0)%css' )
             IF ( ALLOCATED( surf_v(0)%css )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(0)%css
          CASE ( 'surf_v(0)%cssws' )
             IF ( ALLOCATED( surf_v(0)%cssws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(0)%cssws
          CASE ( 'surf_v(0)%qcsws' )
             IF ( ALLOCATED( surf_v(0)%qcsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(0)%qcsws
          CASE ( 'surf_v(0)%ncsws' )
             IF ( ALLOCATED( surf_v(0)%ncsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(0)%ncsws
          CASE ( 'surf_v(0)%qrsws' )
             IF ( ALLOCATED( surf_v(0)%qrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(0)%qrsws
          CASE ( 'surf_v(0)%nrsws' )
             IF ( ALLOCATED( surf_v(0)%nrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(0)%nrsws
          CASE ( 'surf_v(0)%sasws' )
             IF ( ALLOCATED( surf_v(0)%sasws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(0)%sasws
          CASE ( 'surf_v(0)%mom_uv' )
             IF ( ALLOCATED( surf_v(0)%mom_flux_uv )  .AND.  kk == 1 )         &
                READ ( 13 )  surf_v(0)%mom_flux_uv
          CASE ( 'surf_v(0)%mom_w' )
             IF ( ALLOCATED( surf_v(0)%mom_flux_w )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_v(0)%mom_flux_w
          CASE ( 'surf_v(0)%mom_tke' )
             IF ( ALLOCATED( surf_v(0)%mom_flux_tke )  .AND.  kk == 1 )        &
                READ ( 13 )  surf_v(0)%mom_flux_tke

          CASE ( 'surf_v(1)%start_index' )
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_v(1)%start_index
             l = 1
          CASE ( 'surf_v(1)%end_index' )
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_v(1)%end_index
          CASE ( 'surf_v(1)%us' )
             IF ( ALLOCATED( surf_v(1)%us )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(1)%us
          CASE ( 'surf_v(1)%ts' )
             IF ( ALLOCATED( surf_v(1)%ts )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(1)%ts
          CASE ( 'surf_v(1)%qs' )
             IF ( ALLOCATED( surf_v(1)%qs )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(1)%qs
          CASE ( 'surf_v(1)%ss' )
             IF ( ALLOCATED( surf_v(1)%ss )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(1)%ss
          CASE ( 'surf_v(1)%qcs' )
             IF ( ALLOCATED( surf_v(1)%qcs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(1)%qcs
          CASE ( 'surf_v(1)%ncs' )
             IF ( ALLOCATED( surf_v(1)%ncs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(1)%ncs
          CASE ( 'surf_v(1)%qrs' )
             IF ( ALLOCATED( surf_v(1)%qrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(1)%qrs
          CASE ( 'surf_v(1)%nrs' )
             IF ( ALLOCATED( surf_v(1)%nrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(1)%nrs
          CASE ( 'surf_v(1)%ol' )
             IF ( ALLOCATED( surf_v(1)%ol )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(1)%ol
          CASE ( 'surf_v(1)%rib' )
             IF ( ALLOCATED( surf_v(1)%rib )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(1)%rib
          CASE ( 'surf_v(1)%pt_surface' )
             IF ( ALLOCATED( surf_v(1)%pt_surface )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_v(1)%pt_surface
          CASE ( 'surf_v(1)%shf' )
             IF ( ALLOCATED( surf_v(1)%shf )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(1)%shf
          CASE ( 'surf_v(1)%qsws' )
             IF ( ALLOCATED( surf_v(1)%qsws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_v(1)%qsws
          CASE ( 'surf_v(1)%ssws' )
             IF ( ALLOCATED( surf_v(1)%ssws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_v(1)%ssws
          CASE ( 'surf_v(1)%css' )
             IF ( ALLOCATED( surf_v(1)%css )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(1)%css
          CASE ( 'surf_v(1)%cssws' )
             IF ( ALLOCATED( surf_v(1)%cssws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(1)%cssws
          CASE ( 'surf_v(1)%qcsws' )
             IF ( ALLOCATED( surf_v(1)%qcsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(1)%qcsws
          CASE ( 'surf_v(1)%ncsws' )
             IF ( ALLOCATED( surf_v(1)%ncsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(1)%ncsws
          CASE ( 'surf_v(1)%qrsws' )
             IF ( ALLOCATED( surf_v(1)%qrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(1)%qrsws
          CASE ( 'surf_v(1)%nrsws' )
             IF ( ALLOCATED( surf_v(1)%nrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(1)%nrsws
          CASE ( 'surf_v(1)%sasws' )
             IF ( ALLOCATED( surf_v(1)%sasws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(1)%sasws
          CASE ( 'surf_v(1)%mom_uv' )
             IF ( ALLOCATED( surf_v(1)%mom_flux_uv )  .AND.  kk == 1 )         &
                READ ( 13 )  surf_v(1)%mom_flux_uv
          CASE ( 'surf_v(1)%mom_w' )
             IF ( ALLOCATED( surf_v(1)%mom_flux_w )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_v(1)%mom_flux_w
          CASE ( 'surf_v(1)%mom_tke' )
             IF ( ALLOCATED( surf_v(1)%mom_flux_tke )  .AND.  kk == 1 )        &
                READ ( 13 )  surf_v(1)%mom_flux_tke

          CASE ( 'surf_v(2)%start_index' )
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_v(2)%start_index
             l = 2
          CASE ( 'surf_v(2)%end_index' )
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_v(2)%end_index
          CASE ( 'surf_v(2)%us' )
             IF ( ALLOCATED( surf_v(2)%us )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(2)%us
          CASE ( 'surf_v(2)%ts' )
             IF ( ALLOCATED( surf_v(2)%ts )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(2)%ts
          CASE ( 'surf_v(2)%qs' )
             IF ( ALLOCATED( surf_v(2)%qs )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(2)%qs
          CASE ( 'surf_v(2)%ss' )
             IF ( ALLOCATED( surf_v(2)%ss )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(2)%ss
          CASE ( 'surf_v(2)%qcs' )
             IF ( ALLOCATED( surf_v(2)%qcs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(2)%qcs
          CASE ( 'surf_v(2)%ncs' )
             IF ( ALLOCATED( surf_v(2)%ncs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(2)%ncs
          CASE ( 'surf_v(2)%qrs' )
             IF ( ALLOCATED( surf_v(2)%qrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(2)%qrs
          CASE ( 'surf_v(2)%nrs' )
             IF ( ALLOCATED( surf_v(2)%nrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(2)%nrs
          CASE ( 'surf_v(2)%ol' )
             IF ( ALLOCATED( surf_v(2)%ol )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(2)%ol
          CASE ( 'surf_v(2)%rib' )
             IF ( ALLOCATED( surf_v(2)%rib )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(2)%rib
          CASE ( 'surf_v(2)%pt_surface' )
             IF ( ALLOCATED( surf_v(2)%pt_surface )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_v(2)%pt_surface
          CASE ( 'surf_v(2)%shf' )
             IF ( ALLOCATED( surf_v(2)%shf )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(2)%shf
          CASE ( 'surf_v(2)%qsws' )
             IF ( ALLOCATED( surf_v(2)%qsws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_v(2)%qsws
          CASE ( 'surf_v(2)%ssws' )
             IF ( ALLOCATED( surf_v(2)%ssws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_v(2)%ssws
          CASE ( 'surf_v(2)%css' )
             IF ( ALLOCATED( surf_v(2)%css )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(2)%css
          CASE ( 'surf_v(2)%cssws' )
             IF ( ALLOCATED( surf_v(2)%cssws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(2)%cssws
          CASE ( 'surf_v(2)%qcsws' )
             IF ( ALLOCATED( surf_v(2)%qcsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(2)%qcsws
          CASE ( 'surf_v(2)%ncsws' )
             IF ( ALLOCATED( surf_v(2)%ncsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(2)%ncsws
          CASE ( 'surf_v(2)%qrsws' )
             IF ( ALLOCATED( surf_v(2)%qrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(2)%qrsws
          CASE ( 'surf_v(2)%nrsws' )
             IF ( ALLOCATED( surf_v(2)%nrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(2)%nrsws
          CASE ( 'surf_v(2)%sasws' )
             IF ( ALLOCATED( surf_v(2)%sasws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(2)%sasws
          CASE ( 'surf_v(2)%mom_uv' )
             IF ( ALLOCATED( surf_v(2)%mom_flux_uv )  .AND.  kk == 1 )         &
                READ ( 13 )  surf_v(2)%mom_flux_uv
          CASE ( 'surf_v(2)%mom_w' )
             IF ( ALLOCATED( surf_v(2)%mom_flux_w )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_v(2)%mom_flux_w
          CASE ( 'surf_v(2)%mom_tke' )
             IF ( ALLOCATED( surf_v(2)%mom_flux_tke )  .AND.  kk == 1 )        &
                READ ( 13 )  surf_v(2)%mom_flux_tke

          CASE ( 'surf_v(3)%start_index' )
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_v(3)%start_index
             l = 3
          CASE ( 'surf_v(3)%end_index' )
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_v(3)%end_index
          CASE ( 'surf_v(3)%us' )
             IF ( ALLOCATED( surf_v(3)%us )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(3)%us
          CASE ( 'surf_v(3)%ts' )
             IF ( ALLOCATED( surf_v(3)%ts )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(3)%ts
          CASE ( 'surf_v(3)%qs' )
             IF ( ALLOCATED( surf_v(3)%qs )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(3)%qs
          CASE ( 'surf_v(3)%ss' )
             IF ( ALLOCATED( surf_v(3)%ss )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(3)%ss
          CASE ( 'surf_v(3)%qcs' )
             IF ( ALLOCATED( surf_v(3)%qcs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(3)%qcs
          CASE ( 'surf_v(3)%ncs' )
             IF ( ALLOCATED( surf_v(3)%ncs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(3)%ncs
          CASE ( 'surf_v(3)%qrs' )
             IF ( ALLOCATED( surf_v(3)%qrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(3)%qrs
          CASE ( 'surf_v(3)%nrs' )
             IF ( ALLOCATED( surf_v(3)%nrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(3)%nrs
          CASE ( 'surf_v(3)%ol' )
             IF ( ALLOCATED( surf_v(3)%ol )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(3)%ol
          CASE ( 'surf_v(3)%rib' )
             IF ( ALLOCATED( surf_v(3)%rib )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(3)%rib
          CASE ( 'surf_v(3)%pt_surface' )
             IF ( ALLOCATED( surf_v(3)%pt_surface )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_v(3)%pt_surface
          CASE ( 'surf_v(3)%shf' )
             IF ( ALLOCATED( surf_v(3)%shf )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(3)%shf
          CASE ( 'surf_v(3)%qsws' )
             IF ( ALLOCATED( surf_v(3)%qsws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_v(3)%qsws
          CASE ( 'surf_v(3)%ssws' )
             IF ( ALLOCATED( surf_v(3)%ssws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_v(3)%ssws
          CASE ( 'surf_v(3)%css' )
             IF ( ALLOCATED( surf_v(3)%css )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(3)%css
          CASE ( 'surf_v(3)%cssws' )
             IF ( ALLOCATED( surf_v(3)%cssws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(3)%cssws
          CASE ( 'surf_v(3)%qcsws' )
             IF ( ALLOCATED( surf_v(3)%qcsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(3)%qcsws
          CASE ( 'surf_v(3)%ncsws' )
             IF ( ALLOCATED( surf_v(3)%ncsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(3)%ncsws
          CASE ( 'surf_v(3)%qrsws' )
             IF ( ALLOCATED( surf_v(3)%qrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(3)%qrsws
          CASE ( 'surf_v(3)%nrsws' )
             IF ( ALLOCATED( surf_v(3)%nrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(3)%nrsws
          CASE ( 'surf_v(3)%sasws' )
             IF ( ALLOCATED( surf_v(3)%sasws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(3)%sasws
          CASE ( 'surf_v(3)%mom_uv' )
             IF ( ALLOCATED( surf_v(3)%mom_flux_uv )  .AND.  kk == 1 )         &
                READ ( 13 )  surf_v(3)%mom_flux_uv
          CASE ( 'surf_v(3)%mom_w' )
             IF ( ALLOCATED( surf_v(3)%mom_flux_w )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_v(3)%mom_flux_w
          CASE ( 'surf_v(3)%mom_tke' )
             IF ( ALLOCATED( surf_v(3)%mom_flux_tke )  .AND.  kk == 1 )        &
                READ ( 13 )  surf_v(3)%mom_flux_tke

          CASE DEFAULT

                found = .FALSE.

       END SELECT
!
!--    Redistribute surface elements on its respective type.
       IF ( horizontal_surface  .AND.                                          &
            .NOT. INDEX( restart_string(1:length), '%start_index' ) /= 0 )     &
       THEN

          ic = nxlc
          DO  i = nxlf, nxrf
             jc = nysc
             DO  j = nysf, nynf

                surf_match_def  = surf_def_h(l)%end_index(jc,ic) >=            &
                                  surf_def_h(l)%start_index(jc,ic)

                IF ( surf_match_def )  THEN
                   mm = surf_def_h(l)%start_index(jc,ic)
                   DO  m = surf_h(l)%start_index(j,i),                         &
                           surf_h(l)%end_index(j,i)
                      IF ( surf_def_h(l)%end_index(jc,ic) >= mm )              &
                         CALL restore_surface_elements( surf_def_h(l),         &
                                                        mm, surf_h(l), m )
                      mm = mm + 1
                   ENDDO
                ENDIF
               jc = jc + 1
             ENDDO
             ic = ic + 1
          ENDDO
       ELSEIF ( vertical_surface  .AND.                                        &
            .NOT. INDEX( restart_string(1:length), '%start_index' ) /= 0 )     &
       THEN
          ic = nxlc
          DO  i = nxlf, nxrf
             jc = nysc
             DO  j = nysf, nynf

                surf_match_def  = surf_def_v(l)%end_index(jc,ic) >=            &
                                  surf_def_v(l)%start_index(jc,ic)

                IF ( surf_match_def )  THEN
                   mm = surf_def_v(l)%start_index(jc,ic)
                   DO  m = surf_v(l)%start_index(j,i),                         &
                           surf_v(l)%end_index(j,i)
                      IF ( surf_def_v(l)%end_index(jc,ic) >= mm )              &
                         CALL restore_surface_elements( surf_def_v(l),         &
                                                        mm, surf_v(l), m )
                      mm = mm + 1
                   ENDDO
                ENDIF
                jc = jc + 1
             ENDDO
             ic = ic + 1
          ENDDO
       ENDIF


    CONTAINS
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Restores surfacle elements back on its respective type.
!------------------------------------------------------------------------------!
          SUBROUTINE restore_surface_elements( surf_target, m_target,          &
                                               surf_file,   m_file )

             IMPLICIT NONE

             INTEGER(iwp)      ::  m_file      !< respective surface-element index of current surface array
             INTEGER(iwp)      ::  m_target    !< respecitve surface-element index of surface array on file

             TYPE( surf_type ) ::  surf_target !< target surface type
             TYPE( surf_type ) ::  surf_file   !< surface type on file


             IF ( INDEX( restart_string(1:length), '%us' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%us )  .AND.                        &
                     ALLOCATED( surf_file%us   ) )                             &
                   surf_target%us(m_target) = surf_file%us(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%ol' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%ol )  .AND.                        &
                     ALLOCATED( surf_file%ol   ) )                             &
                   surf_target%ol(m_target) = surf_file%ol(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%pt_surface' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%pt_surface )  .AND.                &
                     ALLOCATED( surf_file%pt_surface   ) )                     &
                   surf_target%pt_surface(m_target) = surf_file%pt_surface(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%usws' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%usws )  .AND.                      &
                     ALLOCATED( surf_file%usws   ) )                           &
                   surf_target%usws(m_target) = surf_file%usws(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%vsws' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%vsws )  .AND.                      &
                     ALLOCATED( surf_file%vsws   ) )                           &
                   surf_target%vsws(m_target) = surf_file%vsws(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%ts' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%ts )  .AND.                        &
                     ALLOCATED( surf_file%ts   ) )                             &
                   surf_target%ts(m_target) = surf_file%ts(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%shf' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%shf )  .AND.                       &
                     ALLOCATED( surf_file%shf   ) )                            &
                   surf_target%shf(m_target) = surf_file%shf(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%qs' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%qs )  .AND.                        &
                     ALLOCATED( surf_file%qs   ) )                             &
                   surf_target%qs(m_target) = surf_file%qs(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%qsws' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%qsws )  .AND.                      &
                     ALLOCATED( surf_file%qsws   ) )                           &
                   surf_target%qsws(m_target) = surf_file%qsws(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%ss' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%ss )  .AND.                        &
                     ALLOCATED( surf_file%ss   ) )                             &
                   surf_target%ss(m_target) = surf_file%ss(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%ssws' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%ssws )  .AND.                      &
                     ALLOCATED( surf_file%ssws   ) )                           &
                   surf_target%ssws(m_target) = surf_file%ssws(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%qcs' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%qcs )  .AND.                       &
                     ALLOCATED( surf_file%qcs   ) )                            &
                  surf_target%qcs(m_target) = surf_file%qcs(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%qcsws' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%qcsws )  .AND.                     &
                     ALLOCATED( surf_file%qcsws   ) )                          &
                   surf_target%qcsws(m_target) = surf_file%qcsws(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%ncs' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%ncs )  .AND.                       &
                     ALLOCATED( surf_file%ncs   ) )                            &
                   surf_target%ncs(m_target) = surf_file%ncs(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%ncsws' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%ncsws )  .AND.                     &
                     ALLOCATED( surf_file%ncsws   ) )                          &
                   surf_target%ncsws(m_target) = surf_file%ncsws(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%qrs' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%qrs )  .AND.                       &
                     ALLOCATED( surf_file%qrs   ) )                            &
                  surf_target%qrs(m_target) = surf_file%qrs(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%qrsws' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%qrsws )  .AND.                     &
                     ALLOCATED( surf_file%qrsws   ) )                          &
                   surf_target%qrsws(m_target) = surf_file%qrsws(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%nrs' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%nrs )  .AND.                       &
                     ALLOCATED( surf_file%nrs   ) )                            &
                   surf_target%nrs(m_target) = surf_file%nrs(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%nrsws' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%nrsws )  .AND.                     &
                     ALLOCATED( surf_file%nrsws   ) )                          &
                   surf_target%nrsws(m_target) = surf_file%nrsws(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%sasws' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%sasws )  .AND.                     &
                     ALLOCATED( surf_file%sasws   ) )                          &
                   surf_target%sasws(m_target) = surf_file%sasws(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%mom_uv' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%mom_flux_uv )  .AND.               &
                     ALLOCATED( surf_file%mom_flux_uv   ) )                    &
                   surf_target%mom_flux_uv(m_target) =                         &
                                           surf_file%mom_flux_uv(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%mom_w' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%mom_flux_w )  .AND.                &
                     ALLOCATED( surf_file%mom_flux_w   ) )                     &
                   surf_target%mom_flux_w(m_target) =                          &
                                           surf_file%mom_flux_w(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%mom_tke' ) /= 0 )  THEN
                IF ( ALLOCATED( surf_target%mom_flux_tke )  .AND.              &
                     ALLOCATED( surf_file%mom_flux_tke   ) )                   &
                   surf_target%mom_flux_tke(0:1,m_target) =                    &
                                           surf_file%mom_flux_tke(0:1,m_file)
             ENDIF


          END SUBROUTINE restore_surface_elements


    END SUBROUTINE surface_rrd_local


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Counts the number of surface elements with the same facing, required for
!> reading and writing restart data.
!------------------------------------------------------------------------------!
    SUBROUTINE surface_last_actions

       IMPLICIT NONE
!
!--    Horizontal surfaces
       ns_h_on_file(0) = surf_def_h(0)%ns
       ns_h_on_file(1) = surf_def_h(1)%ns
       ns_h_on_file(2) = surf_def_h(2)%ns
!
!--    Vertical surfaces
       ns_v_on_file(0) = surf_def_v(0)%ns
       ns_v_on_file(1) = surf_def_v(1)%ns
       ns_v_on_file(2) = surf_def_v(2)%ns
       ns_v_on_file(3) = surf_def_v(3)%ns

    END SUBROUTINE surface_last_actions

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Routine maps surface data read from file after restart - 1D arrays
!------------------------------------------------------------------------------!
    SUBROUTINE surface_restore_elements_1d( surf_target, surf_file,            &
                                            start_index_c,                     &
                                            start_index_on_file,               &
                                            end_index_on_file,                 &
                                            nxlc, nysc, nxlf, nxrf, nysf, nynf,&
                                            nys_on_file, nyn_on_file,          &
                                            nxl_on_file,nxr_on_file )

       IMPLICIT NONE

       INTEGER(iwp) ::  i         !< running index along x-direction, refers to former domain size
       INTEGER(iwp) ::  ic        !< running index along x-direction, refers to current domain size
       INTEGER(iwp) ::  j         !< running index along y-direction, refers to former domain size
       INTEGER(iwp) ::  jc        !< running index along y-direction, refers to former domain size
       INTEGER(iwp) ::  m         !< surface-element index on file
       INTEGER(iwp) ::  mm        !< surface-element index on current subdomain
       INTEGER(iwp) ::  nxlc      !< index of left boundary on current subdomain
       INTEGER(iwp) ::  nxlf      !< index of left boundary on former subdomain
       INTEGER(iwp) ::  nxrf      !< index of right boundary on former subdomain
       INTEGER(iwp) ::  nysc      !< index of north boundary on current subdomain
       INTEGER(iwp) ::  nynf      !< index of north boundary on former subdomain
       INTEGER(iwp) ::  nysf      !< index of south boundary on former subdomain

       INTEGER(iwp) ::  nxl_on_file !< leftmost index on file
       INTEGER(iwp) ::  nxr_on_file !< rightmost index on file
       INTEGER(iwp) ::  nyn_on_file !< northmost index on file
       INTEGER(iwp) ::  nys_on_file !< southmost index on file

       INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  start_index_c
       INTEGER(iwp), DIMENSION(nys_on_file:nyn_on_file,nxl_on_file:nxr_on_file) :: &
                            start_index_on_file   !< start index of surface elements on file
       INTEGER(iwp), DIMENSION(nys_on_file:nyn_on_file,nxl_on_file:nxr_on_file) :: &
                            end_index_on_file     !< end index of surface elements on file

       REAL(wp), DIMENSION(:) ::  surf_target !< target surface type
       REAL(wp), DIMENSION(:) ::  surf_file   !< surface type on file

       ic = nxlc
       DO  i = nxlf, nxrf
          jc = nysc
          DO  j = nysf, nynf

             mm = start_index_c(jc,ic)
             DO  m = start_index_on_file(j,i), end_index_on_file(j,i)
                surf_target(mm) = surf_file(m)
                mm = mm + 1
             ENDDO

             jc = jc + 1
          ENDDO
          ic = ic + 1
       ENDDO


    END SUBROUTINE surface_restore_elements_1d

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Routine maps surface data read from file after restart - 2D arrays
!------------------------------------------------------------------------------!
    SUBROUTINE surface_restore_elements_2d( surf_target, surf_file,            &
                                            start_index_c,                     &
                                            start_index_on_file,               &
                                            end_index_on_file,                 &
                                            nxlc, nysc, nxlf, nxrf, nysf, nynf,&
                                            nys_on_file, nyn_on_file,          &
                                            nxl_on_file,nxr_on_file )

       IMPLICIT NONE

       INTEGER(iwp) ::  i         !< running index along x-direction, refers to former domain size
       INTEGER(iwp) ::  ic        !< running index along x-direction, refers to current domain size
       INTEGER(iwp) ::  j         !< running index along y-direction, refers to former domain size
       INTEGER(iwp) ::  jc        !< running index along y-direction, refers to former domain size
       INTEGER(iwp) ::  m         !< surface-element index on file
       INTEGER(iwp) ::  mm        !< surface-element index on current subdomain
       INTEGER(iwp) ::  nxlc      !< index of left boundary on current subdomain
       INTEGER(iwp) ::  nxlf      !< index of left boundary on former subdomain
       INTEGER(iwp) ::  nxrf      !< index of right boundary on former subdomain
       INTEGER(iwp) ::  nysc      !< index of north boundary on current subdomain
       INTEGER(iwp) ::  nynf      !< index of north boundary on former subdomain
       INTEGER(iwp) ::  nysf      !< index of south boundary on former subdomain

       INTEGER(iwp) ::  nxl_on_file !< leftmost index on file
       INTEGER(iwp) ::  nxr_on_file !< rightmost index on file
       INTEGER(iwp) ::  nyn_on_file !< northmost index on file
       INTEGER(iwp) ::  nys_on_file !< southmost index on file

       INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  start_index_c !< start index of surface type
       INTEGER(iwp), DIMENSION(nys_on_file:nyn_on_file,nxl_on_file:nxr_on_file) :: &
                            start_index_on_file   !< start index of surface elements on file
       INTEGER(iwp), DIMENSION(nys_on_file:nyn_on_file,nxl_on_file:nxr_on_file) :: &
                            end_index_on_file     !< end index of surface elements on file

       REAL(wp), DIMENSION(:,:) ::  surf_target !< target surface type
       REAL(wp), DIMENSION(:,:) ::  surf_file   !< surface type on file

       ic = nxlc
       DO  i = nxlf, nxrf
          jc = nysc
          DO  j = nysf, nynf
             mm = start_index_c(jc,ic)
             DO  m = start_index_on_file(j,i), end_index_on_file(j,i)
                surf_target(:,mm) = surf_file(:,m)
                mm = mm + 1
             ENDDO

             jc = jc + 1
          ENDDO
          ic = ic + 1
       ENDDO

    END SUBROUTINE surface_restore_elements_2d


 END MODULE surface_mod
