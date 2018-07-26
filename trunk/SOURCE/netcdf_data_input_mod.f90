!> @file netcdf_data_input_mod.f90
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
! $Id: netcdf_data_input_mod.f90 3089 2018-06-27 13:20:38Z suehring $
! Revise call for message routine in case of local data inconsistencies.
! 
! 3054 2018-06-01 16:08:59Z gronemeier
! Bugfix: force an MPI abort if errors occur while reading building heights
! from ASCII file
! 
! 3053 2018-06-01 12:59:07Z suehring
! Revise checks for variable surface_fraction
! 
! 3051 2018-05-30 17:43:55Z suehring
! - Speed-up NetCDF input
! - Revise input routines and remove NetCDF input via IO-blocks since this is 
!   not working in parallel mode in case blocking collective read operations 
!   are done 
! - Temporarily revoke renaming of input variables in dynamic driver (tend_ug, 
!   tend_vg, zsoil) in order to keep dynamic input file working with current 
!   model version
! - More detailed error messages created
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
! 
! 3041 2018-05-25 10:39:54Z gronemeier
! Add data type for global file attributes
! Add read of global attributes of static driver
!
! 3037 2018-05-24 10:39:29Z gronemeier
! renamed 'depth' to 'zsoil'
!
! 3036 2018-05-24 10:18:26Z gronemeier
! Revision of input vars according to UC2 data standard
!  - renamed 'orography_2D' to 'zt'
!  - renamed 'buildings_2D' to 'buildings_2d'
!  - renamed 'buildings_3D' to 'buildings_3d'
!  - renamed 'leaf_are_density' to 'lad'
!  - renamed 'basal_are_density' to 'bad'
!  - renamed 'root_are_density_lad' to 'root_area_dens_r'
!  - renamed 'root_are_density_lsm' to 'root_area_dens_s'
!  - renamed 'ls_forcing_ug' to 'tend_ug'
!  - renamed 'ls_forcing_vg' to 'tend_vg'
!
! 3019 2018-05-13 07:05:43Z maronga
! Improved reading speed of large NetCDF files
!
! 2963 2018-04-12 14:47:44Z suehring
! - Revise checks for static input variables.
! - Introduce index for vegetation/wall, pavement/green-wall and water/window
!   surfaces, for clearer access of surface fraction, albedo, emissivity, etc. .
!
! 2958 2018-04-11 15:38:13Z suehring
! Synchronize longitude and latitude between nested model domains, values are
! taken from the root model.
!
! 2955 2018-04-09 15:14:01Z suehring
! Extend checks for consistent setting of buildings, its ID and type.
! Add log-points to measure CPU time of NetCDF data input.
!
! 2953 2018-04-09 11:26:02Z suehring
! Bugfix in checks for initialization data
!
! 2947 2018-04-04 18:01:41Z suehring
! Checks for dynamic input revised
!
! 2946 2018-04-04 17:01:23Z suehring
! Bugfix for revision 2945, perform checks only if dynamic input file is
! available.
!
! 2945 2018-04-04 16:27:14Z suehring
! - Mimic for topography input slightly revised, in order to enable consistency
!   checks
! - Add checks for dimensions in dynamic input file and move already existing
!   checks
!
! 2938 2018-03-27 15:52:42Z suehring
! Initial read of geostrophic wind components from dynamic driver.
!
! 2773 2018-01-30 14:12:54Z suehring
! Revise checks for surface_fraction.
!
! 2925 2018-03-23 14:54:11Z suehring
! Check for further inconsistent settings of surface_fractions.
! Some messages slightly rephrased and error numbers renamed.
!
! 2898 2018-03-15 13:03:01Z suehring
! Check if each building has a type. Further, check if dimensions in static
! input file match the model dimensions.
!
! 2897 2018-03-15 11:47:16Z suehring
! Relax restrictions for topography input, terrain and building heights can be
! input separately and are not mandatory any more.
!
! 2874 2018-03-13 10:55:42Z knoop
! Bugfix: wrong placement of netcdf cpp-macros fixed
!
! 2794 2018-02-07 14:09:43Z knoop
! Check if 3D building input is consistent to numeric grid.
!
! 2773 2018-01-30 14:12:54Z suehring
! - Enable initialization with 3D topography.
! - Move check for correct initialization in nesting mode to check_parameters.
!
! 2772 2018-01-29 13:10:35Z suehring
! Initialization of simulation independent on land-surface model.
!
! 2746 2018-01-15 12:06:04Z suehring
! Read plant-canopy variables independently on land-surface model usage
!
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
!
! 2711 2017-12-20 17:04:49Z suehring
! Rename subroutine close_file to avoid double-naming.
!
! 2700 2017-12-15 14:12:35Z suehring
!
! 2696 2017-12-14 17:12:51Z kanani
! Initial revision (suehring)
!
!
!
!
! Authors:
! --------
! @author Matthias Suehring
!
! Description:
! ------------
!> Modulue contains routines to input data according to Palm input data
!> standart using dynamic and static input files.
!>
!> @todo - Order input alphabetically
!> @todo - Revise error messages and error numbers
!> @todo - Input of missing quantities (chemical species, emission rates)
!> @todo - Defninition and input of still missing variable attributes
!> @todo - Input of initial geostrophic wind profiles with cyclic conditions.
!------------------------------------------------------------------------------!
 MODULE netcdf_data_input_mod

    USE control_parameters,                                                    &
        ONLY:  coupling_char, io_blocks, io_group

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE kinds

#if defined ( __netcdf )
    USE NETCDF
#endif

    USE pegrid

    USE surface_mod,                                                           &
        ONLY:  ind_pav_green, ind_veg_wall, ind_wat_win
!
!-- Define type for dimensions.
    TYPE dims_xy
       INTEGER(iwp) :: nx                             !< dimension length in x
       INTEGER(iwp) :: ny                             !< dimension length in y
       INTEGER(iwp) :: nz                             !< dimension length in z
       REAL(wp), DIMENSION(:), ALLOCATABLE :: x       !< dimension array in x
       REAL(wp), DIMENSION(:), ALLOCATABLE :: y       !< dimension array in y
       REAL(wp), DIMENSION(:), ALLOCATABLE :: z       !< dimension array in z
    END TYPE dims_xy
!
!-- Define data type for nesting in larger-scale models like COSMO.
!-- Data type comprises u, v, w, pt, and q at lateral and top boundaries.
    TYPE force_type

       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names

       INTEGER(iwp) ::  nt     !< number of time levels in dynamic input file
       INTEGER(iwp) ::  nzu    !< number of vertical levels on scalar grid in dynamic input file
       INTEGER(iwp) ::  nzw    !< number of vertical levels on w grid in dynamic input file
       INTEGER(iwp) ::  tind   !< time index for reference time in large-scale forcing data
       INTEGER(iwp) ::  tind_p !< time index for following time in large-scale forcing data

       LOGICAL      ::  init         = .FALSE.
       LOGICAL      ::  interpolated = .FALSE.
       LOGICAL      ::  from_file    = .FALSE.

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surface_pressure !< time dependent surface pressure
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  time             !< time levels in dynamic input file
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zu_atmos         !< vertical levels at scalar grid in dynamic input file
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zw_atmos         !< vertical levels at w grid in dynamic input file

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ug         !< domain-averaged geostrophic component
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vg         !< domain-averaged geostrophic component

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_left   !< u-component at left boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_left   !< v-component at left boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_left   !< w-component at left boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  q_left   !< mixing ratio at left boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pt_left  !< potentital temperautre at left boundary

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_north  !< u-component at north boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_north  !< v-component at north boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_north  !< w-component at north boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  q_north  !< mixing ratio at north boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pt_north !< potentital temperautre at north boundary

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_right  !< u-component at right boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_right  !< v-component at right boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_right  !< w-component at right boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  q_right  !< mixing ratio at right boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pt_right !< potentital temperautre at right boundary

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_south  !< u-component at south boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_south  !< v-component at south boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_south  !< w-component at south boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  q_south  !< mixing ratio at south boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pt_south !< potentital temperautre at south boundary

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_top    !< u-component at top boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_top    !< v-component at top boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_top    !< w-component at top boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  q_top    !< mixing ratio at top boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pt_top   !< potentital temperautre at top boundary

    END TYPE force_type

    TYPE init_type

       CHARACTER(LEN=23) ::  origin_time !< reference time of input data

       INTEGER(iwp) ::  lod_msoil !< level of detail - soil moisture
       INTEGER(iwp) ::  lod_pt    !< level of detail - pt
       INTEGER(iwp) ::  lod_q     !< level of detail - q
       INTEGER(iwp) ::  lod_tsoil !< level of detail - soil temperature
       INTEGER(iwp) ::  lod_u     !< level of detail - u-component
       INTEGER(iwp) ::  lod_v     !< level of detail - v-component
       INTEGER(iwp) ::  lod_w     !< level of detail - w-component
       INTEGER(iwp) ::  nx        !< number of scalar grid points along x in dynamic input file
       INTEGER(iwp) ::  nxu       !< number of u grid points along x in dynamic input file
       INTEGER(iwp) ::  ny        !< number of scalar grid points along y in dynamic input file
       INTEGER(iwp) ::  nyv       !< number of v grid points along y in dynamic input file
       INTEGER(iwp) ::  nzs       !< number of vertical soil levels in dynamic input file
       INTEGER(iwp) ::  nzu       !< number of vertical levels on scalar grid in dynamic input file
       INTEGER(iwp) ::  nzw       !< number of vertical levels on w grid in dynamic input file

       LOGICAL ::  from_file_msoil  = .FALSE. !< flag indicating whether soil moisture is already initialized from file
       LOGICAL ::  from_file_pt     = .FALSE. !< flag indicating whether pt is already initialized from file
       LOGICAL ::  from_file_q      = .FALSE. !< flag indicating whether q is already initialized from file
       LOGICAL ::  from_file_tsoil  = .FALSE. !< flag indicating whether soil temperature is already initialized from file
       LOGICAL ::  from_file_u      = .FALSE. !< flag indicating whether u is already initialized from file
       LOGICAL ::  from_file_ug     = .FALSE. !< flag indicating whether ug is already initialized from file
       LOGICAL ::  from_file_v      = .FALSE. !< flag indicating whether v is already initialized from file
       LOGICAL ::  from_file_vg     = .FALSE. !< flag indicating whether ug is already initialized from file
       LOGICAL ::  from_file_w      = .FALSE. !< flag indicating whether w is already initialized from file

       REAL(wp) ::  fill_msoil       !< fill value for soil moisture
       REAL(wp) ::  fill_pt          !< fill value for pt
       REAL(wp) ::  fill_q           !< fill value for q
       REAL(wp) ::  fill_tsoil       !< fill value for soil temperature
       REAL(wp) ::  fill_u           !< fill value for u
       REAL(wp) ::  fill_v           !< fill value for v
       REAL(wp) ::  fill_w           !< fill value for w
       REAL(wp) ::  latitude         !< latitude of the southern model boundary
       REAL(wp) ::  longitude        !< longitude of the western model boundary
       REAL(wp) ::  origin_x         !< x position of the western model boundary
       REAL(wp) ::  origin_y         !< y position of the northern model boundary
       REAL(wp) ::  origin_z         !< reference height of input data
       REAL(wp) ::  rotation_angle   !< rotation angle of input data

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  msoil_init   !< initial vertical profile of soil moisture
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_init      !< initial vertical profile of pt
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  q_init       !< initial vertical profile of q
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  tsoil_init   !< initial vertical profile of soil temperature
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  u_init       !< initial vertical profile of u
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ug_init      !< initial vertical profile of ug
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  v_init       !< initial vertical profile of v
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  vg_init      !< initial vertical profile of ug
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  w_init       !< initial vertical profile of w
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z_soil       !< vertical levels in soil in dynamic input file, used for interpolation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zu_atmos     !< vertical levels at scalar grid in dynamic input file, used for interpolation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zw_atmos     !< vertical levels at w grid in dynamic input file, used for interpolation


       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  msoil        !< initial 3d soil moisture provide by Inifor and interpolated onto soil grid
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  tsoil        !< initial 3d soil temperature provide by Inifor and interpolated onto soil grid

    END TYPE init_type

!
!-- Define data structures for different input data types.
!-- 8-bit Integer 2D
    TYPE int_2d_8bit
       INTEGER(KIND=1) ::  fill = -127                      !< fill value
       INTEGER(KIND=1), DIMENSION(:,:), ALLOCATABLE ::  var !< respective variable

       LOGICAL ::  from_file = .FALSE.  !< flag indicating whether an input variable is available and read from file or default values are used
    END TYPE int_2d_8bit
!
!-- 32-bit Integer 2D
    TYPE int_2d_32bit
       INTEGER(iwp) ::  fill = -9999                      !< fill value
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  var  !< respective variable

       LOGICAL ::  from_file = .FALSE. !< flag indicating whether an input variable is available and read from file or default values are used
    END TYPE int_2d_32bit

!
!-- Define data type to read 2D real variables
    TYPE real_2d
       LOGICAL ::  from_file = .FALSE.  !< flag indicating whether an input variable is available and read from file or default values are used

       REAL(wp) ::  fill = -9999.9_wp                !< fill value
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  var !< respective variable
    END TYPE real_2d

!
!-- Define data type to read 2D real variables
    TYPE real_3d
       LOGICAL ::  from_file = .FALSE.  !< flag indicating whether an input variable is available and read from file or default values are used

       INTEGER(iwp) ::  nz   !< number of grid points along vertical dimension

       REAL(wp) ::  fill = -9999.9_wp                  !< fill value
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  var !< respective variable
    END TYPE real_3d
!
!-- Define data structure where the dimension and type of the input depends
!-- on the given level of detail.
!-- For buildings, the input is either 2D float, or 3d byte.
    TYPE build_in
       INTEGER(iwp)    ::  lod = 1                               !< level of detail
       INTEGER(KIND=1) ::  fill2 = -127                          !< fill value for lod = 2
       INTEGER(iwp)    ::  nz                                    !< number of vertical layers in file
       INTEGER(KIND=1), DIMENSION(:,:,:), ALLOCATABLE ::  var_3d !< 3d variable (lod = 2)

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z                 !< vertical coordinate for 3D building, used for consistency check

       LOGICAL ::  from_file = .FALSE.  !< flag indicating whether an input variable is available and read from file or default values are used

       REAL(wp)                              ::  fill1 = -9999.9_wp !< fill values for lod = 1
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  var_2d             !< 2d variable (lod = 1)
    END TYPE build_in

!
!-- For soil_type, the input is either 2D or 3D one-byte integer.
    TYPE soil_in
       INTEGER(iwp)                                   ::  lod = 1      !< level of detail
       INTEGER(KIND=1)                                ::  fill = -127  !< fill value for lod = 2
       INTEGER(KIND=1), DIMENSION(:,:), ALLOCATABLE   ::  var_2d       !< 2d variable (lod = 1)
       INTEGER(KIND=1), DIMENSION(:,:,:), ALLOCATABLE ::  var_3d       !< 3d variable (lod = 2)

       LOGICAL ::  from_file = .FALSE.  !< flag indicating whether an input variable is available and read from file or default values are used
    END TYPE soil_in

!
!-- Define data type for fractions between surface types
    TYPE fracs
       INTEGER(iwp)                            ::  nf             !< total number of fractions
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nfracs         !< dimension array for fraction

       LOGICAL ::  from_file = .FALSE. !< flag indicating whether an input variable is available and read from file or default values are used

       REAL(wp)                                ::  fill = -9999.9_wp !< fill value
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  frac              !< respective fraction between different surface types
    END TYPE fracs
!
!-- Data type for parameter lists, Depending on the given level of detail,
!-- the input is 3D or 4D
    TYPE pars
       INTEGER(iwp)                            ::  lod = 1         !< level of detail
       INTEGER(iwp)                            ::  np              !< total number of parameters
       INTEGER(iwp)                            ::  nz              !< vertical dimension - number of soil layers
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  layers          !< dimension array for soil layers
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  pars            !< dimension array for parameters

       LOGICAL ::  from_file = .FALSE.  !< flag indicating whether an input variable is available and read from file or default values are used

       REAL(wp)                                  ::  fill = -9999.9_wp !< fill value
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  pars_xy           !< respective parameters, level of detail = 1
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  pars_xyz          !< respective parameters, level of detail = 2
    END TYPE pars
!
!-- Define type for global file attributes
!-- Please refer to the PALM data standard for a detailed description of each
!-- attribute.
    TYPE global_atts_type
       CHARACTER(LEN=12 ) ::  acronym                            !< acronym of institution
       CHARACTER(LEN=7)   ::  acronym_char = 'acronym'           !< name of attribute
       CHARACTER(LEN=200) ::  author                             !< first name, last name, email adress
       CHARACTER(LEN=6)   ::  author_char = 'author'             !< name of attribute
       CHARACTER(LEN=12 ) ::  campaign                           !< name of campaign
       CHARACTER(LEN=8)   ::  campaign_char = 'campaign'         !< name of attribute
       CHARACTER(LEN=200) ::  comment                            !< comment to data
       CHARACTER(LEN=7)   ::  comment_char = 'comment'           !< name of attribute
       CHARACTER(LEN=200) ::  contact_person                     !< first name, last name, email adress
       CHARACTER(LEN=14)  ::  contact_person_char = 'contact_person'  !< name of attribute
       CHARACTER(LEN=200) ::  conventions = 'CF-1.7'             !< netCDF convention
       CHARACTER(LEN=11)  ::  conventions_char = 'Conventions'   !< name of attribute
       CHARACTER(LEN=23 ) ::  creation_time                      !< creation time of data set
       CHARACTER(LEN=13)  ::  creation_time_char = 'creation_time'  !< name of attribute
       CHARACTER(LEN=16 ) ::  data_content                       !< content of data set
       CHARACTER(LEN=12)  ::  data_content_char = 'data_content' !< name of attribute
       CHARACTER(LEN=200) ::  dependencies                       !< dependencies of data set
       CHARACTER(LEN=12)  ::  dependencies_char = 'dependencies' !< name of attribute
       CHARACTER(LEN=200) ::  history                            !< information about data processing
       CHARACTER(LEN=7)   ::  history_char = 'history'           !< name of attribute
       CHARACTER(LEN=200) ::  institution                        !< name of responsible institution
       CHARACTER(LEN=11)  ::  institution_char = 'institution'   !< name of attribute
       CHARACTER(LEN=200) ::  keywords                           !< keywords of data set
       CHARACTER(LEN=8)   ::  keywords_char = 'keywords'         !< name of attribute
       CHARACTER(LEN=200) ::  license                            !< license of data set
       CHARACTER(LEN=7)   ::  license_char = 'license'           !< name of attribute
       CHARACTER(LEN=200) ::  location                           !< place which refers to data set
       CHARACTER(LEN=8)   ::  location_char = 'location'         !< name of attribute
       CHARACTER(LEN=10)  ::  origin_lat_char = 'origin_lat'     !< name of attribute
       CHARACTER(LEN=10)  ::  origin_lon_char = 'origin_lon'     !< name of attribute
       CHARACTER(LEN=23 ) ::  origin_time                        !< reference time
       CHARACTER(LEN=11)  ::  origin_time_char = 'origin_time'   !< name of attribute
       CHARACTER(LEN=8)   ::  origin_x_char = 'origin_x'         !< name of attribute
       CHARACTER(LEN=8)   ::  origin_y_char = 'origin_y'         !< name of attribute
       CHARACTER(LEN=8)   ::  origin_z_char = 'origin_z'         !< name of attribute
       CHARACTER(LEN=12)  ::  palm_version_char = 'palm_version' !< name of attribute
       CHARACTER(LEN=200) ::  references                         !< literature referring to data set
       CHARACTER(LEN=10)  ::  references_char = 'references'     !< name of attribute
       CHARACTER(LEN=14)  ::  rotation_angle_char = 'rotation_angle'  !< name of attribute
       CHARACTER(LEN=12 ) ::  site                               !< name of model domain
       CHARACTER(LEN=4)   ::  site_char = 'site'                 !< name of attribute
       CHARACTER(LEN=200) ::  source                             !< source of data set
       CHARACTER(LEN=6)   ::  source_char = 'source'             !< name of attribute
       CHARACTER(LEN=200) ::  title                              !< title of data set
       CHARACTER(LEN=5)   ::  title_char = 'title'               !< name of attribute
       CHARACTER(LEN=7)   ::  version_char = 'version'           !< name of attribute

       INTEGER(iwp) ::  version              !< version of data set

       REAL(wp) ::  origin_lat               !< latitude of lower left corner
       REAL(wp) ::  origin_lon               !< longitude of lower left corner
       REAL(wp) ::  origin_x                 !< easting (UTM coordinate) of lower left corner
       REAL(wp) ::  origin_y                 !< northing (UTM coordinate) of lower left corner
       REAL(wp) ::  origin_z                 !< reference height
       REAL(wp) ::  palm_version             !< PALM version of data set
       REAL(wp) ::  rotation_angle           !< rotation angle of coordinate system of data set
    END TYPE global_atts_type
!
!-- Define variables
    TYPE(dims_xy)    ::  dim_static  !< data structure for x, y-dimension in static input file

    TYPE(force_type) ::  force     !< data structure for data input at lateral and top boundaries (provided by Inifor)

    TYPE(init_type) ::  init_3d    !< data structure for the initialization of the 3D flow and soil fields
    TYPE(init_type) ::  init_model !< data structure for the initialization of the model

!
!-- Define 2D variables of type NC_BYTE
    TYPE(int_2d_8bit)  ::  albedo_type_f     !< input variable for albedo type
    TYPE(int_2d_8bit)  ::  building_type_f   !< input variable for building type
    TYPE(int_2d_8bit)  ::  pavement_type_f   !< input variable for pavenment type
    TYPE(int_2d_8bit)  ::  street_crossing_f !< input variable for water type
    TYPE(int_2d_8bit)  ::  street_type_f     !< input variable for water type
    TYPE(int_2d_8bit)  ::  vegetation_type_f !< input variable for vegetation type
    TYPE(int_2d_8bit)  ::  water_type_f      !< input variable for water type

!
!-- Define 2D variables of type NC_INT
    TYPE(int_2d_32bit) ::  building_id_f     !< input variable for building ID
!
!-- Define 2D variables of type NC_FLOAT
    TYPE(real_2d) ::  terrain_height_f       !< input variable for terrain height
!
!-- Define 3D variables of type NC_FLOAT
    TYPE(real_3d) ::  basal_area_density_f    !< input variable for basal area density - resolved vegetation
    TYPE(real_3d) ::  leaf_area_density_f     !< input variable for leaf area density - resolved vegetation
    TYPE(real_3d) ::  root_area_density_lad_f !< input variable for root area density - resolved vegetation
    TYPE(real_3d) ::  root_area_density_lsm_f !< input variable for root area density - parametrized vegetation

!
!-- Define input variable for buildings
    TYPE(build_in) ::  buildings_f           !< input variable for buildings
!
!-- Define input variables for soil_type
    TYPE(soil_in)  ::  soil_type_f           !< input variable for soil type

    TYPE(fracs) ::  surface_fraction_f       !< input variable for surface fraction

    TYPE(pars)  ::  albedo_pars_f              !< input variable for albedo parameters
    TYPE(pars)  ::  building_pars_f            !< input variable for building parameters
    TYPE(pars)  ::  pavement_pars_f            !< input variable for pavement parameters
    TYPE(pars)  ::  pavement_subsurface_pars_f !< input variable for pavement parameters
    TYPE(pars)  ::  soil_pars_f                !< input variable for soil parameters
    TYPE(pars)  ::  vegetation_pars_f          !< input variable for vegetation parameters
    TYPE(pars)  ::  water_pars_f               !< input variable for water parameters


    CHARACTER(LEN=3)  ::  char_lod  = 'lod'         !< name of level-of-detail attribute in NetCDF file

    CHARACTER(LEN=10) ::  char_fill = '_FillValue'        !< name of fill value attribute in NetCDF file

    CHARACTER(LEN=100) ::  input_file_static  = 'PIDS_STATIC'  !< Name of file which comprises static input data
    CHARACTER(LEN=100) ::  input_file_dynamic = 'PIDS_DYNAMIC' !< Name of file which comprises dynamic input data

    INTEGER(iwp) ::  nc_stat         !< return value of nf90 function call

    LOGICAL ::  input_pids_static  = .FALSE.   !< Flag indicating whether Palm-input-data-standard file containing static information exists
    LOGICAL ::  input_pids_dynamic = .FALSE.   !< Flag indicating whether Palm-input-data-standard file containing dynamic information exists

    LOGICAL ::  collective_read = .FALSE.      !< Enable NetCDF collective read

    TYPE(global_atts_type) ::  input_file_atts !< global attributes of input file

    SAVE

    PRIVATE

    INTERFACE netcdf_data_input_interpolate
       MODULE PROCEDURE netcdf_data_input_interpolate_1d
       MODULE PROCEDURE netcdf_data_input_interpolate_1d_soil
       MODULE PROCEDURE netcdf_data_input_interpolate_2d
       MODULE PROCEDURE netcdf_data_input_interpolate_3d
    END INTERFACE netcdf_data_input_interpolate

    INTERFACE netcdf_data_input_check_dynamic
       MODULE PROCEDURE netcdf_data_input_check_dynamic
    END INTERFACE netcdf_data_input_check_dynamic

    INTERFACE netcdf_data_input_check_static
       MODULE PROCEDURE netcdf_data_input_check_static
    END INTERFACE netcdf_data_input_check_static

    INTERFACE netcdf_data_input_inquire_file
       MODULE PROCEDURE netcdf_data_input_inquire_file
    END INTERFACE netcdf_data_input_inquire_file

    INTERFACE netcdf_data_input_init
       MODULE PROCEDURE netcdf_data_input_init
    END INTERFACE netcdf_data_input_init

    INTERFACE netcdf_data_input_init_3d
       MODULE PROCEDURE netcdf_data_input_init_3d
    END INTERFACE netcdf_data_input_init_3d

    INTERFACE netcdf_data_input_lsf
       MODULE PROCEDURE netcdf_data_input_lsf
    END INTERFACE netcdf_data_input_lsf

    INTERFACE netcdf_data_input_surface_data
       MODULE PROCEDURE netcdf_data_input_surface_data
    END INTERFACE netcdf_data_input_surface_data

    INTERFACE netcdf_data_input_topo
       MODULE PROCEDURE netcdf_data_input_topo
    END INTERFACE netcdf_data_input_topo

    INTERFACE get_variable
       MODULE PROCEDURE get_variable_1d_int
       MODULE PROCEDURE get_variable_1d_real
       MODULE PROCEDURE get_variable_2d_int8
       MODULE PROCEDURE get_variable_2d_int32
       MODULE PROCEDURE get_variable_2d_real
       MODULE PROCEDURE get_variable_3d_int8
       MODULE PROCEDURE get_variable_3d_real
       MODULE PROCEDURE get_variable_3d_real_dynamic
!        MODULE PROCEDURE get_variable_3d_real_v
       MODULE PROCEDURE get_variable_4d_real
    END INTERFACE get_variable

    INTERFACE get_variable_pr
       MODULE PROCEDURE get_variable_pr
    END INTERFACE get_variable_pr

    INTERFACE get_attribute
       MODULE PROCEDURE get_attribute_real
       MODULE PROCEDURE get_attribute_int8
       MODULE PROCEDURE get_attribute_int32
       MODULE PROCEDURE get_attribute_string
    END INTERFACE get_attribute

!
!-- Public variables
    PUBLIC albedo_pars_f, albedo_type_f, basal_area_density_f, buildings_f,    &
           building_id_f, building_pars_f, building_type_f, force, init_3d,    &
           init_model, input_file_static, input_pids_static,                   &
           input_pids_dynamic, leaf_area_density_f,                            &
           pavement_pars_f, pavement_subsurface_pars_f, pavement_type_f,       &
           root_area_density_lad_f, root_area_density_lsm_f, soil_pars_f,      &
           soil_type_f, street_crossing_f, street_type_f, surface_fraction_f,  &
           terrain_height_f, vegetation_pars_f, vegetation_type_f,             &
           water_pars_f, water_type_f

!
!-- Public subroutines
    PUBLIC netcdf_data_input_check_dynamic, netcdf_data_input_check_static,    &
           netcdf_data_input_inquire_file,                                     &
           netcdf_data_input_init, netcdf_data_input_init_3d,                  &
           netcdf_data_input_interpolate, netcdf_data_input_lsf,               &
           netcdf_data_input_surface_data, netcdf_data_input_topo

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Inquires whether NetCDF input files according to Palm-input-data standard
!> exist. Moreover, basic checks are performed.
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_inquire_file

       USE control_parameters,                                                 &
           ONLY:  land_surface, message_string, topo_no_distinct, urban_surface

       IMPLICIT NONE

       LOGICAL ::  check_nest  !< flag indicating whether a check passed or not

#if defined ( __netcdf )
       INQUIRE( FILE = TRIM( input_file_static ) // TRIM( coupling_char ),     &
                EXIST = input_pids_static  )
       INQUIRE( FILE = TRIM( input_file_dynamic ) // TRIM( coupling_char ),    &
                EXIST = input_pids_dynamic )
#endif

!
!--    As long as topography can be input via ASCII format, no distinction
!--    between building and terrain can be made. This case, classify all
!--    surfaces as default type. Same in case land-surface and urban-surface
!--    model are not applied.
       IF ( .NOT. input_pids_static )  THEN
          topo_no_distinct = .TRUE.
       ENDIF

    END SUBROUTINE netcdf_data_input_inquire_file

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads global attributes required for initialization of the model.
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_init

       IMPLICIT NONE

       INTEGER(iwp) ::  id_mod   !< NetCDF id of input file
       INTEGER(iwp) ::  ii       !< running index for IO blocks

       IF ( .NOT. input_pids_static )  RETURN

#if defined ( __netcdf )
!
!--    Open file in read-only mode
       CALL open_read_file( TRIM( input_file_static ) //                       &
                            TRIM( coupling_char ), id_mod )
!
!--    Read global attributes
       CALL get_attribute( id_mod, input_file_atts%origin_lat_char,            &
                           input_file_atts%origin_lat, .TRUE. )

       CALL get_attribute( id_mod, input_file_atts%origin_lon_char,            &
                           input_file_atts%origin_lon, .TRUE. )

       CALL get_attribute( id_mod, input_file_atts%origin_time_char,           &
                           input_file_atts%origin_time, .TRUE. )

       CALL get_attribute( id_mod, input_file_atts%origin_x_char,              &
                           input_file_atts%origin_x, .TRUE. )

       CALL get_attribute( id_mod, input_file_atts%origin_y_char,              &
                           input_file_atts%origin_y, .TRUE. )

       CALL get_attribute( id_mod, input_file_atts%origin_z_char,              &
                           input_file_atts%origin_z, .TRUE. )

       CALL get_attribute( id_mod, input_file_atts%rotation_angle_char,        &
                           input_file_atts%rotation_angle, .TRUE. )

!
!--    Finally, close input file
       CALL close_input_file( id_mod )
#endif
!
!--    Copy latitude, longitude, origin_z, rotation angle on init type
       init_model%latitude        = input_file_atts%origin_lat
       init_model%longitude       = input_file_atts%origin_lon
       init_model%origin_time     = input_file_atts%origin_time  
       init_model%origin_x        = input_file_atts%origin_x
       init_model%origin_y        = input_file_atts%origin_y
       init_model%origin_z        = input_file_atts%origin_z  
       init_model%rotation_angle  = input_file_atts%rotation_angle  
            
!
!--    In case of nested runs, each model domain might have different longitude
!--    and latitude, which would result in different Coriolis parameters and
!--    sun-zenith angles. To avoid this, longitude and latitude in each model
!--    domain will be set to the values of the root model. Please note, this
!--    synchronization is required already here.
#if defined( __parallel )
       CALL MPI_BCAST( init_model%latitude,  1, MPI_REAL, 0,                   &
                       MPI_COMM_WORLD, ierr )
       CALL MPI_BCAST( init_model%longitude, 1, MPI_REAL, 0,                   &
                       MPI_COMM_WORLD, ierr )
#endif


    END SUBROUTINE netcdf_data_input_init

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads surface classification data, such as vegetation and soil type, etc. .
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_surface_data

       USE control_parameters,                                                 &
           ONLY:  bc_lr_cyc, bc_ns_cyc, land_surface, message_string,          &
                  plant_canopy, urban_surface

       USE indices,                                                            &
           ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg


       IMPLICIT NONE

       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names  !< variable names in static input file

       INTEGER(iwp) ::  id_surf   !< NetCDF id of input file
       INTEGER(iwp) ::  k         !< running index along z-direction
       INTEGER(iwp) ::  k2        !< running index
       INTEGER(iwp) ::  num_vars  !< number of variables in input file
       INTEGER(iwp) ::  nz_soil   !< number of soil layers in file

       INTEGER(iwp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  var_exchange_int !< dummy variables used to exchange 32-bit Integer arrays
       INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE  ::  var_dum_int_3d !< dummy variables used to exchange real arrays

       REAL(wp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  var_exchange_real !< dummy variables used to exchange real arrays

       REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  var_dum_real_3d !< dummy variables used to exchange real arrays
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  var_dum_real_4d !< dummy variables used to exchange real arrays

!
!--    If not static input file is available, skip this routine
       IF ( .NOT. input_pids_static )  RETURN
!
!--    Measure CPU time
       CALL cpu_log( log_point_s(82), 'NetCDF input', 'start' )
!
!--    Read plant canopy variables.
       IF ( plant_canopy )  THEN
#if defined ( __netcdf )
!
!--       Open file in read-only mode
          CALL open_read_file( TRIM( input_file_static ) //                    &
                               TRIM( coupling_char ) , id_surf )
!
!--       At first, inquire all variable names.
!--       This will be used to check whether an optional input variable
!--       exist or not.
          CALL inquire_num_variables( id_surf, num_vars )

          ALLOCATE( var_names(1:num_vars) )
          CALL inquire_variable_names( id_surf, var_names )

!
!--       Read leaf area density - resolved vegetation
          IF ( check_existence( var_names, 'lad' ) )  THEN
             leaf_area_density_f%from_file = .TRUE.
             CALL get_attribute( id_surf, char_fill,                           &
                                 leaf_area_density_f%fill,                     &
                                 .FALSE., 'lad' )
!
!--          Inquire number of vertical vegetation layer
             CALL get_dimension_length( id_surf, leaf_area_density_f%nz,       &
                                        'zlad' )
!
!--          Allocate variable for leaf-area density
             ALLOCATE( leaf_area_density_f%var( 0:leaf_area_density_f%nz-1,    &
                                                nys:nyn,nxl:nxr) )

             CALL get_variable( id_surf, 'lad', leaf_area_density_f%var,       &
                                nxl, nxr, nys, nyn,                            &
                                0, leaf_area_density_f%nz-1 )

          ELSE
             leaf_area_density_f%from_file = .FALSE.
          ENDIF

!
!--       Read basal area density - resolved vegetation
          IF ( check_existence( var_names, 'bad' ) )  THEN
             basal_area_density_f%from_file = .TRUE.
             CALL get_attribute( id_surf, char_fill,                           &
                                 basal_area_density_f%fill,                    &
                                 .FALSE., 'bad' )
!
!--          Inquire number of vertical vegetation layer
             CALL get_dimension_length( id_surf, basal_area_density_f%nz,      &
                                        'zlad' )
!
!--          Allocate variable
             ALLOCATE( basal_area_density_f%var(0:basal_area_density_f%nz-1,   &
                                                nys:nyn,nxl:nxr) )

             CALL get_variable( id_surf, 'bad', basal_area_density_f%var,      &
                                nxl, nxr, nys, nyn,                            &
                                0,  basal_area_density_f%nz-1 )
          ELSE
             basal_area_density_f%from_file = .FALSE.
          ENDIF

!
!--       Read root area density - resolved vegetation
          IF ( check_existence( var_names, 'root_area_dens_r' ) )  THEN
             root_area_density_lad_f%from_file = .TRUE.
             CALL get_attribute( id_surf, char_fill,                           &
                                 root_area_density_lad_f%fill,                 &
                                 .FALSE., 'root_area_dens_r' )
!
!--          Inquire number of vertical soil layers
             CALL get_dimension_length( id_surf,                               &
                                        root_area_density_lad_f%nz,            &
                                        'zsoil' )
!
!--          Allocate variable
             ALLOCATE( root_area_density_lad_f%var                             &
                                         (0:root_area_density_lad_f%nz-1,      &
                                          nys:nyn,nxl:nxr) )

             CALL get_variable( id_surf, 'root_area_dens_r',                   &
                                root_area_density_lad_f%var,                   &
                                nxl, nxr, nys, nyn,                            &
                                0,  root_area_density_lad_f%nz-1 )
          ELSE
             root_area_density_lad_f%from_file = .FALSE.
          ENDIF
!
!--       Finally, close input file
          CALL close_input_file( id_surf )
#endif
       ENDIF
!
!--    Deallocate variable list. Will be re-allocated in case further
!--    variables are read from file.
       IF ( ALLOCATED( var_names ) )  DEALLOCATE( var_names )
!
!--    Skip the following if no land-surface or urban-surface module are
!--    applied. This case, no one of the following variables is used anyway.
       IF (  .NOT. land_surface  .OR.  .NOT. urban_surface )  RETURN
!
!--    Initialize dummy arrays used for ghost-point exchange
       var_exchange_int  = 0
       var_exchange_real = 0.0_wp

#if defined ( __netcdf )
!
!--    Open file in read-only mode
       CALL open_read_file( TRIM( input_file_static ) //                       &
                            TRIM( coupling_char ) , id_surf )
!
!--    Inquire all variable names.
!--    This will be used to check whether an optional input variable exist
!--    or not.
       CALL inquire_num_variables( id_surf, num_vars )

       ALLOCATE( var_names(1:num_vars) )
       CALL inquire_variable_names( id_surf, var_names )
!
!--    Read vegetation type and required attributes
       IF ( check_existence( var_names, 'vegetation_type' ) )  THEN
          vegetation_type_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              vegetation_type_f%fill,                          &
                              .FALSE., 'vegetation_type' )

          ALLOCATE ( vegetation_type_f%var(nys:nyn,nxl:nxr)  )

          CALL get_variable( id_surf, 'vegetation_type',                       &
                             vegetation_type_f%var, nxl, nxr, nys, nyn )
       ELSE
          vegetation_type_f%from_file = .FALSE.
       ENDIF

!
!--    Read soil type and required attributes
       IF ( check_existence( var_names, 'soil_type' ) )  THEN
             soil_type_f%from_file = .TRUE.
!
!--       Note, lod is currently not on file; skip for the moment
!           CALL get_attribute( id_surf, char_lod,                       &
!                                      soil_type_f%lod,                  &
!                                      .FALSE., 'soil_type' )
          CALL get_attribute( id_surf, char_fill,                              &
                              soil_type_f%fill,                                &
                              .FALSE., 'soil_type' )

          IF ( soil_type_f%lod == 1 )  THEN

             ALLOCATE ( soil_type_f%var_2d(nys:nyn,nxl:nxr)  )

             CALL get_variable( id_surf, 'soil_type', soil_type_f%var_2d,      &
                                nxl, nxr, nys, nyn )

          ELSEIF ( soil_type_f%lod == 2 )  THEN
!
!--          Obtain number of soil layers from file.
             CALL get_dimension_length( id_surf, nz_soil, 'zsoil' )

             ALLOCATE ( soil_type_f%var_3d(0:nz_soil,nys:nyn,nxl:nxr) )

             CALL get_variable( id_surf, 'soil_type', soil_type_f%var_3d,      &
                                nxl, nxr, nys, nyn, 0, nz_soil )
 
          ENDIF
       ELSE
          soil_type_f%from_file = .FALSE.
       ENDIF

!
!--    Read pavement type and required attributes
       IF ( check_existence( var_names, 'pavement_type' ) )  THEN
          pavement_type_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              pavement_type_f%fill, .FALSE.,                   &
                              'pavement_type' )

          ALLOCATE ( pavement_type_f%var(nys:nyn,nxl:nxr)  )

          CALL get_variable( id_surf, 'pavement_type', pavement_type_f%var,    &
                             nxl, nxr, nys, nyn )
       ELSE
          pavement_type_f%from_file = .FALSE.
       ENDIF

!
!--    Read water type and required attributes
       IF ( check_existence( var_names, 'water_type' ) )  THEN
          water_type_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill, water_type_f%fill,           &
                              .FALSE., 'water_type' )

          ALLOCATE ( water_type_f%var(nys:nyn,nxl:nxr)  )

          CALL get_variable( id_surf, 'water_type', water_type_f%var,          &
                             nxl, nxr, nys, nyn )

       ELSE
          water_type_f%from_file = .FALSE.
       ENDIF
!
!--    Read relative surface fractions of vegetation, pavement and water.
       IF ( check_existence( var_names, 'surface_fraction' ) )  THEN
          surface_fraction_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              surface_fraction_f%fill,                         &
                              .FALSE., 'surface_fraction' )
!
!--       Inquire number of surface fractions
          CALL get_dimension_length( id_surf,                                  &
                                     surface_fraction_f%nf,                    &
                                     'nsurface_fraction' )
!
!--       Allocate dimension array and input array for surface fractions
          ALLOCATE( surface_fraction_f%nfracs(0:surface_fraction_f%nf-1) )
          ALLOCATE( surface_fraction_f%frac(0:surface_fraction_f%nf-1,         &
                                            nys:nyn,nxl:nxr) )
!
!--       Get dimension of surface fractions
          CALL get_variable( id_surf, 'nsurface_fraction',                     &
                             surface_fraction_f%nfracs )
!
!--       Read surface fractions
          CALL get_variable( id_surf, 'surface_fraction',                      &
                             surface_fraction_f%frac, nxl, nxr, nys, nyn,      &
                             0, surface_fraction_f%nf-1 )
       ELSE
          surface_fraction_f%from_file = .FALSE.
       ENDIF
!
!--    Read building parameters and related information
       IF ( check_existence( var_names, 'building_pars' ) )  THEN
          building_pars_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              building_pars_f%fill,                            &
                              .FALSE., 'building_pars' )
!
!--       Inquire number of building parameters
          CALL get_dimension_length( id_surf,                                  &
                                     building_pars_f%np,                       &
                                     'nbuilding_pars' )
!
!--       Allocate dimension array and input array for building parameters
          ALLOCATE( building_pars_f%pars(0:building_pars_f%np-1) )
          ALLOCATE( building_pars_f%pars_xy(0:building_pars_f%np-1,            &
                                            nys:nyn,nxl:nxr) )
!
!--       Get dimension of building parameters
          CALL get_variable( id_surf, 'nbuilding_pars',                        &
                             building_pars_f%pars )
!
!--       Read building_pars
          CALL get_variable( id_surf, 'building_pars',                         &
                             building_pars_f%pars_xy, nxl, nxr, nys, nyn,      &
                             0, building_pars_f%np-1 )
       ELSE
          building_pars_f%from_file = .FALSE.
       ENDIF

!
!--    Read albedo type and required attributes
       IF ( check_existence( var_names, 'albedo_type' ) )  THEN
          albedo_type_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill, albedo_type_f%fill,          &
                              .FALSE.,  'albedo_type' )

          ALLOCATE ( albedo_type_f%var(nys:nyn,nxl:nxr)  )
          
          CALL get_variable( id_surf, 'albedo_type', albedo_type_f%var,        &
                             nxl, nxr, nys, nyn )
       ELSE
          albedo_type_f%from_file = .FALSE.
       ENDIF
!
!--    Read albedo parameters and related information
       IF ( check_existence( var_names, 'albedo_pars' ) )  THEN
          albedo_pars_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill, albedo_pars_f%fill,          &
                              .FALSE., 'albedo_pars' )
!
!--       Inquire number of albedo parameters
          CALL get_dimension_length( id_surf, albedo_pars_f%np,                &
                                     'nalbedo_pars' )
!
!--       Allocate dimension array and input array for albedo parameters
          ALLOCATE( albedo_pars_f%pars(0:albedo_pars_f%np-1) )
          ALLOCATE( albedo_pars_f%pars_xy(0:albedo_pars_f%np-1,                &
                                          nys:nyn,nxl:nxr) )
!
!--       Get dimension of albedo parameters
          CALL get_variable( id_surf, 'nalbedo_pars', albedo_pars_f%pars )

          CALL get_variable( id_surf, 'albedo_pars', albedo_pars_f%pars_xy,    &
                             nxl, nxr, nys, nyn,                               &
                             0, albedo_pars_f%np-1 )
       ELSE
          albedo_pars_f%from_file = .FALSE.
       ENDIF

!
!--    Read pavement parameters and related information
       IF ( check_existence( var_names, 'pavement_pars' ) )  THEN
          pavement_pars_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              pavement_pars_f%fill,                            &
                              .FALSE., 'pavement_pars' )
!
!--       Inquire number of pavement parameters
          CALL get_dimension_length( id_surf, pavement_pars_f%np,              &
                                     'npavement_pars' )
!
!--       Allocate dimension array and input array for pavement parameters
          ALLOCATE( pavement_pars_f%pars(0:pavement_pars_f%np-1) )
          ALLOCATE( pavement_pars_f%pars_xy(0:pavement_pars_f%np-1,            &
                                            nys:nyn,nxl:nxr) )
!
!--       Get dimension of pavement parameters
          CALL get_variable( id_surf, 'npavement_pars', pavement_pars_f%pars )

          CALL get_variable( id_surf, 'pavement_pars', pavement_pars_f%pars_xy,&
                             nxl, nxr, nys, nyn,                               &
                             0, pavement_pars_f%np-1 )
       ELSE
          pavement_pars_f%from_file = .FALSE.
       ENDIF

!
!--    Read pavement subsurface parameters and related information
       IF ( check_existence( var_names, 'pavement_subsurface_pars' ) )         &
       THEN
          pavement_subsurface_pars_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              pavement_subsurface_pars_f%fill,                 &
                              .FALSE., 'pavement_subsurface_pars' )
!
!--       Inquire number of parameters
          CALL get_dimension_length( id_surf,                                  &
                                     pavement_subsurface_pars_f%np,            &
                                     'npavement_subsurface_pars' )
!
!--       Inquire number of soil layers
          CALL get_dimension_length( id_surf,                                  &
                                     pavement_subsurface_pars_f%nz,            &
                                     'zsoil' )
!
!--       Allocate dimension array and input array for pavement parameters
          ALLOCATE( pavement_subsurface_pars_f%pars                            &
                            (0:pavement_subsurface_pars_f%np-1) )
          ALLOCATE( pavement_subsurface_pars_f%pars_xyz                        &
                            (0:pavement_subsurface_pars_f%np-1,                &
                             0:pavement_subsurface_pars_f%nz-1,                &
                             nys:nyn,nxl:nxr) )
!
!--       Get dimension of pavement parameters
          CALL get_variable( id_surf, 'npavement_subsurface_pars',             &
                             pavement_subsurface_pars_f%pars )

          CALL get_variable( id_surf, 'pavement_subsurface_pars',              &
                             pavement_subsurface_pars_f%pars_xyz,              &
                             nxl, nxr, nys, nyn,                               &
                             0, pavement_subsurface_pars_f%nz-1,               &
                             0, pavement_subsurface_pars_f%np-1 )
       ELSE
          pavement_subsurface_pars_f%from_file = .FALSE.
       ENDIF


!
!--    Read vegetation parameters and related information
       IF ( check_existence( var_names, 'vegetation_pars' ) )  THEN
          vegetation_pars_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              vegetation_pars_f%fill,                          &
                              .FALSE.,  'vegetation_pars' )
!
!--       Inquire number of vegetation parameters
          CALL get_dimension_length( id_surf, vegetation_pars_f%np,            &
                                     'nvegetation_pars' )
!
!--       Allocate dimension array and input array for surface fractions
          ALLOCATE( vegetation_pars_f%pars(0:vegetation_pars_f%np-1) )
          ALLOCATE( vegetation_pars_f%pars_xy(0:vegetation_pars_f%np-1,        &
                                              nys:nyn,nxl:nxr) )
!
!--       Get dimension of the parameters
          CALL get_variable( id_surf, 'nvegetation_pars',                      &
                             vegetation_pars_f%pars )

          CALL get_variable( id_surf, 'vegetation_pars',                       &
                             vegetation_pars_f%pars_xy, nxl, nxr, nys, nyn,    &
                             0, vegetation_pars_f%np-1 )
       ELSE
          vegetation_pars_f%from_file = .FALSE.
       ENDIF

!
!--    Read root parameters/distribution and related information
       IF ( check_existence( var_names, 'soil_pars' ) )  THEN
          soil_pars_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              soil_pars_f%fill,                                &
                              .FALSE., 'soil_pars' )

          CALL get_attribute( id_surf, char_lod,                               &
                              soil_pars_f%lod,                                 &
                              .FALSE., 'soil_pars' )

!
!--       Inquire number of soil parameters
          CALL get_dimension_length( id_surf,                                  &
                                     soil_pars_f%np,                           &
                                     'nsoil_pars' )
!
!--       Read parameters array
          ALLOCATE( soil_pars_f%pars(0:soil_pars_f%np-1) )
          CALL get_variable( id_surf, 'nsoil_pars', soil_pars_f%pars )

!
!--       In case of level of detail 2, also inquire number of vertical
!--       soil layers, allocate memory and read the respective dimension
          IF ( soil_pars_f%lod == 2 )  THEN
             CALL get_dimension_length( id_surf, soil_pars_f%nz, 'zsoil' )

             ALLOCATE( soil_pars_f%layers(0:soil_pars_f%nz-1) )
             CALL get_variable( id_surf, 'zsoil', soil_pars_f%layers )

          ENDIF

!
!--       Read soil parameters, depending on level of detail
          IF ( soil_pars_f%lod == 1 )  THEN
             ALLOCATE( soil_pars_f%pars_xy(0:soil_pars_f%np-1,                 &
                                           nys:nyn,nxl:nxr) )
                  
             CALL get_variable( id_surf, 'soil_pars', soil_pars_f%pars_xy,     &
                                nxl, nxr, nys, nyn, 0, soil_pars_f%np-1 )

          ELSEIF ( soil_pars_f%lod == 2 )  THEN
             ALLOCATE( soil_pars_f%pars_xyz(0:soil_pars_f%np-1,                &
                                            0:soil_pars_f%nz-1,                &
                                            nys:nyn,nxl:nxr) )
             CALL get_variable( id_surf, 'soil_pars',                          &
                                soil_pars_f%pars_xyz,                          &
                                nxl, nxr, nys, nyn, 0, soil_pars_f%nz-1,       &
                                0, soil_pars_f%np-1 )

          ENDIF
       ELSE
          soil_pars_f%from_file = .FALSE.
       ENDIF

!
!--    Read water parameters and related information
       IF ( check_existence( var_names, 'water_pars' ) )  THEN
          water_pars_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              water_pars_f%fill,                               &
                              .FALSE., 'water_pars' )
!
!--       Inquire number of water parameters
          CALL get_dimension_length( id_surf,                                  &
                                     water_pars_f%np,                          &
                                     'nwater_pars' )
!
!--       Allocate dimension array and input array for water parameters
          ALLOCATE( water_pars_f%pars(0:water_pars_f%np-1) )
          ALLOCATE( water_pars_f%pars_xy(0:water_pars_f%np-1,                  &
                                         nys:nyn,nxl:nxr) )
!
!--       Get dimension of water parameters
          CALL get_variable( id_surf, 'nwater_pars', water_pars_f%pars )

          CALL get_variable( id_surf, 'water_pars', water_pars_f%pars_xy,      &
                             nxl, nxr, nys, nyn, 0, water_pars_f%np-1 )
       ELSE
          water_pars_f%from_file = .FALSE.
       ENDIF
!
!--    Read root area density - parametrized vegetation
       IF ( check_existence( var_names, 'root_area_dens_s' ) )  THEN
          root_area_density_lsm_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              root_area_density_lsm_f%fill,                    &
                              .FALSE., 'root_area_dens_s' )
!
!--       Obtain number of soil layers from file and allocate variable
          CALL get_dimension_length( id_surf, root_area_density_lsm_f%nz,      &
                                     'zsoil' )
          ALLOCATE( root_area_density_lsm_f%var                                &
                                        (0:root_area_density_lsm_f%nz-1,       &
                                         nys:nyn,nxl:nxr) )

!
!--       Read root-area density
          CALL get_variable( id_surf, 'root_area_dens_s',                      &
                             root_area_density_lsm_f%var,                      &
                             nxl, nxr, nys, nyn,                               &
                             0, root_area_density_lsm_f%nz-1 )

       ELSE
          root_area_density_lsm_f%from_file = .FALSE.
       ENDIF
!
!--    Read street type and street crossing
       IF ( check_existence( var_names, 'street_type' ) )  THEN
          street_type_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              street_type_f%fill, .FALSE.,                     &
                              'street_type' )

          ALLOCATE ( street_type_f%var(nys:nyn,nxl:nxr)  )
          
          CALL get_variable( id_surf, 'street_type', street_type_f%var,        &
                             nxl, nxr, nys, nyn )
       ELSE
          street_type_f%from_file = .FALSE.
       ENDIF

       IF ( check_existence( var_names, 'street_crossing' ) )  THEN
          street_crossing_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              street_crossing_f%fill, .FALSE.,                 &
                              'street_crossing' )

          ALLOCATE ( street_crossing_f%var(nys:nyn,nxl:nxr)  )

          CALL get_variable( id_surf, 'street_crossing',                       &
                             street_crossing_f%var, nxl, nxr, nys, nyn )

       ELSE
          street_crossing_f%from_file = .FALSE.
       ENDIF
!
!--    Still missing: root_resolved and building_surface_pars.
!--    Will be implemented as soon as they are available.

!
!--    Finally, close input file
       CALL close_input_file( id_surf )
#endif
!
!--    End of CPU measurement
       CALL cpu_log( log_point_s(82), 'NetCDF input', 'stop' )
!
!--    Exchange 1 ghost points for surface variables. Please note, ghost point
!--    exchange for 3D parameter lists should be revised by using additional
!--    MPI datatypes or rewriting exchange_horiz.
!--    Moreover, varialbes will be resized in the following, including ghost
!--    points.
!--    Start with 2D Integer variables. Please note, for 8-bit integer
!--    variables must be swapt to 32-bit integer before calling exchange_horiz.
       IF ( albedo_type_f%from_file )  THEN
          var_exchange_int                  = INT( albedo_type_f%fill, KIND = 1 )
          var_exchange_int(nys:nyn,nxl:nxr) =                                  &
                            INT( albedo_type_f%var(nys:nyn,nxl:nxr), KIND = 4 )
          CALL exchange_horiz_2d_int( var_exchange_int, nys, nyn, nxl, nxr, nbgp )
          DEALLOCATE( albedo_type_f%var )
          ALLOCATE( albedo_type_f%var(nysg:nyng,nxlg:nxrg) )
          albedo_type_f%var = INT( var_exchange_int, KIND = 1 )
       ENDIF
       IF ( pavement_type_f%from_file )  THEN
          var_exchange_int                  = INT( pavement_type_f%fill, KIND = 1 )
          var_exchange_int(nys:nyn,nxl:nxr) =                                  &
                          INT( pavement_type_f%var(nys:nyn,nxl:nxr), KIND = 4 )
          CALL exchange_horiz_2d_int( var_exchange_int, nys, nyn, nxl, nxr, nbgp )
          DEALLOCATE( pavement_type_f%var )
          ALLOCATE( pavement_type_f%var(nysg:nyng,nxlg:nxrg) )
          pavement_type_f%var = INT( var_exchange_int, KIND = 1 )
       ENDIF
       IF ( soil_type_f%from_file  .AND.  ALLOCATED( soil_type_f%var_2d ) )  THEN
          var_exchange_int                  = INT( soil_type_f%fill, KIND = 1 )
          var_exchange_int(nys:nyn,nxl:nxr) =                                  &
                            INT( soil_type_f%var_2d(nys:nyn,nxl:nxr), KIND = 4 )
          CALL exchange_horiz_2d_int( var_exchange_int, nys, nyn, nxl, nxr, nbgp )
          DEALLOCATE( soil_type_f%var_2d )
          ALLOCATE( soil_type_f%var_2d(nysg:nyng,nxlg:nxrg) )
          soil_type_f%var_2d = INT( var_exchange_int, KIND = 1 )
       ENDIF
       IF ( vegetation_type_f%from_file )  THEN
          var_exchange_int                  = INT( vegetation_type_f%fill, KIND = 1 )
          var_exchange_int(nys:nyn,nxl:nxr) =                                  &
                        INT( vegetation_type_f%var(nys:nyn,nxl:nxr), KIND = 4 )
          CALL exchange_horiz_2d_int( var_exchange_int, nys, nyn, nxl, nxr, nbgp )
          DEALLOCATE( vegetation_type_f%var )
          ALLOCATE( vegetation_type_f%var(nysg:nyng,nxlg:nxrg) )
          vegetation_type_f%var = INT( var_exchange_int, KIND = 1 )
       ENDIF
       IF ( water_type_f%from_file )  THEN
          var_exchange_int                  = INT( water_type_f%fill, KIND = 1 )
          var_exchange_int(nys:nyn,nxl:nxr) =                                  &
                         INT( water_type_f%var(nys:nyn,nxl:nxr), KIND = 4 )
          CALL exchange_horiz_2d_int( var_exchange_int, nys, nyn, nxl, nxr, nbgp )
          DEALLOCATE( water_type_f%var )
          ALLOCATE( water_type_f%var(nysg:nyng,nxlg:nxrg) )
          water_type_f%var = INT( var_exchange_int, KIND = 1 )
       ENDIF
!
!--    Exchange 1 ghost point for 3/4-D variables. For the sake of simplicity,
!--    loop further dimensions to use 2D exchange routines.
!--    This should be revised later by introducing new MPI datatypes.
       IF ( soil_type_f%from_file  .AND.  ALLOCATED( soil_type_f%var_3d ) )    &
       THEN
          ALLOCATE( var_dum_int_3d(0:nz_soil,nys:nyn,nxl:nxr) )
          var_dum_int_3d = soil_type_f%var_3d
          DEALLOCATE( soil_type_f%var_3d )
          ALLOCATE( soil_type_f%var_3d(0:nz_soil,nysg:nyng,nxlg:nxrg) )
          soil_type_f%var_3d = soil_type_f%fill

          DO  k = 0, nz_soil
             var_exchange_int(nys:nyn,nxl:nxr) = var_dum_int_3d(k,nys:nyn,nxl:nxr)
             CALL exchange_horiz_2d_int( var_exchange_int, nys, nyn, nxl, nxr, nbgp )
             soil_type_f%var_3d(k,:,:) = INT( var_exchange_int(:,:), KIND = 1 )
          ENDDO
          DEALLOCATE( var_dum_int_3d )
       ENDIF

       IF ( surface_fraction_f%from_file )  THEN
          ALLOCATE( var_dum_real_3d(0:surface_fraction_f%nf-1,nys:nyn,nxl:nxr) )
          var_dum_real_3d = surface_fraction_f%frac
          DEALLOCATE( surface_fraction_f%frac )
          ALLOCATE( surface_fraction_f%frac(0:surface_fraction_f%nf-1,         &
                                            nysg:nyng,nxlg:nxrg) )
          surface_fraction_f%frac = surface_fraction_f%fill

          DO  k = 0, surface_fraction_f%nf-1
             var_exchange_real(nys:nyn,nxl:nxr) = var_dum_real_3d(k,nys:nyn,nxl:nxr)
             CALL exchange_horiz_2d( var_exchange_real, nbgp )
             surface_fraction_f%frac(k,:,:) = var_exchange_real(:,:)
          ENDDO
          DEALLOCATE( var_dum_real_3d )
       ENDIF

       IF ( building_pars_f%from_file )  THEN
          ALLOCATE( var_dum_real_3d(0:building_pars_f%np-1,nys:nyn,nxl:nxr) )
          var_dum_real_3d = building_pars_f%pars_xy
          DEALLOCATE( building_pars_f%pars_xy )
          ALLOCATE( building_pars_f%pars_xy(0:building_pars_f%np-1,            &
                                            nysg:nyng,nxlg:nxrg) )
          building_pars_f%pars_xy = building_pars_f%fill
          DO  k = 0, building_pars_f%np-1
             var_exchange_real(nys:nyn,nxl:nxr) =                              &
                                              var_dum_real_3d(k,nys:nyn,nxl:nxr)
             CALL exchange_horiz_2d( var_exchange_real, nbgp )
             building_pars_f%pars_xy(k,:,:) = var_exchange_real(:,:)
          ENDDO
          DEALLOCATE( var_dum_real_3d )
       ENDIF

       IF ( albedo_pars_f%from_file )  THEN
          ALLOCATE( var_dum_real_3d(0:albedo_pars_f%np-1,nys:nyn,nxl:nxr) )
          var_dum_real_3d = albedo_pars_f%pars_xy
          DEALLOCATE( albedo_pars_f%pars_xy )
          ALLOCATE( albedo_pars_f%pars_xy(0:albedo_pars_f%np-1,                &
                                          nysg:nyng,nxlg:nxrg) )
          albedo_pars_f%pars_xy = albedo_pars_f%fill
          DO  k = 0, albedo_pars_f%np-1
             var_exchange_real(nys:nyn,nxl:nxr) =                              &
                                              var_dum_real_3d(k,nys:nyn,nxl:nxr)
             CALL exchange_horiz_2d( var_exchange_real, nbgp )
             albedo_pars_f%pars_xy(k,:,:) = var_exchange_real(:,:)
          ENDDO
          DEALLOCATE( var_dum_real_3d )
       ENDIF

       IF ( pavement_pars_f%from_file )  THEN
          ALLOCATE( var_dum_real_3d(0:pavement_pars_f%np-1,nys:nyn,nxl:nxr) )
          var_dum_real_3d = pavement_pars_f%pars_xy
          DEALLOCATE( pavement_pars_f%pars_xy )
          ALLOCATE( pavement_pars_f%pars_xy(0:pavement_pars_f%np-1,            &
                                            nysg:nyng,nxlg:nxrg) )
          pavement_pars_f%pars_xy = pavement_pars_f%fill
          DO  k = 0, pavement_pars_f%np-1
             var_exchange_real(nys:nyn,nxl:nxr) =                              &
                                              var_dum_real_3d(k,nys:nyn,nxl:nxr)
             CALL exchange_horiz_2d( var_exchange_real, nbgp )
             pavement_pars_f%pars_xy(k,:,:) = var_exchange_real(:,:)
          ENDDO
          DEALLOCATE( var_dum_real_3d )
       ENDIF

       IF ( vegetation_pars_f%from_file )  THEN
          ALLOCATE( var_dum_real_3d(0:vegetation_pars_f%np-1,nys:nyn,nxl:nxr) )
          var_dum_real_3d = vegetation_pars_f%pars_xy
          DEALLOCATE( vegetation_pars_f%pars_xy )
          ALLOCATE( vegetation_pars_f%pars_xy(0:vegetation_pars_f%np-1,        &
                                            nysg:nyng,nxlg:nxrg) )
          vegetation_pars_f%pars_xy = vegetation_pars_f%fill
          DO  k = 0, vegetation_pars_f%np-1
             var_exchange_real(nys:nyn,nxl:nxr) =                              &
                                              var_dum_real_3d(k,nys:nyn,nxl:nxr)
             CALL exchange_horiz_2d( var_exchange_real, nbgp )
             vegetation_pars_f%pars_xy(k,:,:) = var_exchange_real(:,:)
          ENDDO
          DEALLOCATE( var_dum_real_3d )
       ENDIF

       IF ( water_pars_f%from_file )  THEN
          ALLOCATE( var_dum_real_3d(0:water_pars_f%np-1,nys:nyn,nxl:nxr) )
          var_dum_real_3d = water_pars_f%pars_xy
          DEALLOCATE( water_pars_f%pars_xy )
          ALLOCATE( water_pars_f%pars_xy(0:water_pars_f%np-1,                  &
                                         nysg:nyng,nxlg:nxrg) )
          water_pars_f%pars_xy = water_pars_f%fill
          DO  k = 0, water_pars_f%np-1
             var_exchange_real(nys:nyn,nxl:nxr) =                              &
                                              var_dum_real_3d(k,nys:nyn,nxl:nxr)
             CALL exchange_horiz_2d( var_exchange_real, nbgp )
             water_pars_f%pars_xy(k,:,:) = var_exchange_real(:,:)
          ENDDO
          DEALLOCATE( var_dum_real_3d )
       ENDIF

       IF ( root_area_density_lsm_f%from_file )  THEN
          ALLOCATE( var_dum_real_3d(0:root_area_density_lsm_f%nz-1,nys:nyn,nxl:nxr) )
          var_dum_real_3d = root_area_density_lsm_f%var
          DEALLOCATE( root_area_density_lsm_f%var )
          ALLOCATE( root_area_density_lsm_f%var(0:root_area_density_lsm_f%nz-1,&
                                                nysg:nyng,nxlg:nxrg) )
          root_area_density_lsm_f%var = root_area_density_lsm_f%fill

          DO  k = 0, root_area_density_lsm_f%nz-1
             var_exchange_real(nys:nyn,nxl:nxr) =                              &
                                              var_dum_real_3d(k,nys:nyn,nxl:nxr)
             CALL exchange_horiz_2d( var_exchange_real, nbgp )
             root_area_density_lsm_f%var(k,:,:) = var_exchange_real(:,:)
          ENDDO
          DEALLOCATE( var_dum_real_3d )
       ENDIF

       IF ( soil_pars_f%from_file )  THEN
          IF ( soil_pars_f%lod == 1 )  THEN

             ALLOCATE( var_dum_real_3d(0:soil_pars_f%np-1,nys:nyn,nxl:nxr) )
             var_dum_real_3d = soil_pars_f%pars_xy
             DEALLOCATE( soil_pars_f%pars_xy )
             ALLOCATE( soil_pars_f%pars_xy(0:soil_pars_f%np-1,                 &
                                            nysg:nyng,nxlg:nxrg) )
             soil_pars_f%pars_xy = soil_pars_f%fill

             DO  k = 0, soil_pars_f%np-1
                var_exchange_real(nys:nyn,nxl:nxr) =                           &
                                              var_dum_real_3d(k,nys:nyn,nxl:nxr)
                CALL exchange_horiz_2d( var_exchange_real, nbgp )
                soil_pars_f%pars_xy(k,:,:) = var_exchange_real(:,:)
             ENDDO
             DEALLOCATE( var_dum_real_3d )
          ELSEIF ( soil_pars_f%lod == 2 )  THEN
             ALLOCATE( var_dum_real_4d(0:soil_pars_f%np-1,                     &
                                       0:soil_pars_f%nz-1,                     &
                                       nys:nyn,nxl:nxr) )
             var_dum_real_4d = soil_pars_f%pars_xyz
             DEALLOCATE( soil_pars_f%pars_xyz )
             ALLOCATE( soil_pars_f%pars_xyz(0:soil_pars_f%np-1,                &
                                            0:soil_pars_f%nz-1,                &
                                            nysg:nyng,nxlg:nxrg) )
             soil_pars_f%pars_xyz = soil_pars_f%fill

             DO  k2 = 0, soil_pars_f%nz-1
                DO  k = 0, soil_pars_f%np-1
                   var_exchange_real(nys:nyn,nxl:nxr) =                        &
                                           var_dum_real_4d(k,k2,nys:nyn,nxl:nxr)
                   CALL exchange_horiz_2d( var_exchange_real, nbgp )

                   soil_pars_f%pars_xyz(k,k2,:,:) = var_exchange_real(:,:)
                ENDDO
             ENDDO
             DEALLOCATE( var_dum_real_4d )
          ENDIF
       ENDIF

       IF ( pavement_subsurface_pars_f%from_file )  THEN
          ALLOCATE( var_dum_real_4d(0:pavement_subsurface_pars_f%np-1,         &
                                    0:pavement_subsurface_pars_f%nz-1,         &
                                    nys:nyn,nxl:nxr) )
          var_dum_real_4d = pavement_subsurface_pars_f%pars_xyz
          DEALLOCATE( pavement_subsurface_pars_f%pars_xyz )
          ALLOCATE( pavement_subsurface_pars_f%pars_xyz                        &
                                          (0:pavement_subsurface_pars_f%np-1,  &
                                           0:pavement_subsurface_pars_f%nz-1,  &
                                           nysg:nyng,nxlg:nxrg) )
          pavement_subsurface_pars_f%pars_xyz = pavement_subsurface_pars_f%fill

          DO  k2 = 0, pavement_subsurface_pars_f%nz-1
             DO  k = 0, pavement_subsurface_pars_f%np-1
                var_exchange_real(nys:nyn,nxl:nxr) =                           &
                                          var_dum_real_4d(k,k2,nys:nyn,nxl:nxr)
                CALL exchange_horiz_2d( var_exchange_real, nbgp )
                pavement_subsurface_pars_f%pars_xyz(k,k2,:,:) =                &
                                                        var_exchange_real(:,:)
             ENDDO
          ENDDO
          DEALLOCATE( var_dum_real_4d )
       ENDIF

!
!--    In case of non-cyclic boundary conditions, set Neumann conditions at the
!--    lateral boundaries.
       IF ( .NOT. bc_ns_cyc )  THEN
          IF ( nys == 0  )  THEN
             IF ( albedo_type_f%from_file )                                    &
                albedo_type_f%var(-1,:) = albedo_type_f%var(0,:)
             IF ( pavement_type_f%from_file )                                  &
                pavement_type_f%var(-1,:) = pavement_type_f%var(0,:)
             IF ( soil_type_f%from_file )  THEN
                IF ( ALLOCATED( soil_type_f%var_2d ) )  THEN
                   soil_type_f%var_2d(-1,:) = soil_type_f%var_2d(0,:)
                ELSE
                   soil_type_f%var_3d(:,-1,:) = soil_type_f%var_3d(:,0,:)
                ENDIF
             ENDIF
             IF ( vegetation_type_f%from_file )                                &
                vegetation_type_f%var(-1,:) = vegetation_type_f%var(0,:)
             IF ( water_type_f%from_file )                                     &
                water_type_f%var(-1,:) = water_type_f%var(0,:)
             IF ( surface_fraction_f%from_file )                               &
                surface_fraction_f%frac(:,-1,:) = surface_fraction_f%frac(:,0,:)
             IF ( building_pars_f%from_file )                                  &
                building_pars_f%pars_xy(:,-1,:) = building_pars_f%pars_xy(:,0,:)
             IF ( albedo_pars_f%from_file )                                    &
                albedo_pars_f%pars_xy(:,-1,:) = albedo_pars_f%pars_xy(:,0,:)
             IF ( pavement_pars_f%from_file )                                  &
                pavement_pars_f%pars_xy(:,-1,:) = pavement_pars_f%pars_xy(:,0,:)
             IF ( vegetation_pars_f%from_file )                                &
                vegetation_pars_f%pars_xy(:,-1,:) =                            &
                                               vegetation_pars_f%pars_xy(:,0,:)
             IF ( water_pars_f%from_file )                                     &
                water_pars_f%pars_xy(:,-1,:) = water_pars_f%pars_xy(:,0,:)
             IF ( root_area_density_lsm_f%from_file )                          &
                root_area_density_lsm_f%var(:,-1,:) =                          &
                                            root_area_density_lsm_f%var(:,0,:)
             IF ( soil_pars_f%from_file )  THEN
                IF ( soil_pars_f%lod == 1 )  THEN
                   soil_pars_f%pars_xy(:,-1,:) = soil_pars_f%pars_xy(:,0,:)
                ELSE
                   soil_pars_f%pars_xyz(:,:,-1,:) = soil_pars_f%pars_xyz(:,:,0,:)
                ENDIF
             ENDIF
             IF ( pavement_subsurface_pars_f%from_file )                       &
                pavement_subsurface_pars_f%pars_xyz(:,:,-1,:) =                &
                                   pavement_subsurface_pars_f%pars_xyz(:,:,0,:)
          ENDIF

          IF ( nyn == ny )  THEN
             IF ( albedo_type_f%from_file )                                    &
                albedo_type_f%var(ny+1,:) = albedo_type_f%var(ny,:)
             IF ( pavement_type_f%from_file )                                  &
                pavement_type_f%var(ny+1,:) = pavement_type_f%var(ny,:)
             IF ( soil_type_f%from_file )  THEN
                IF ( ALLOCATED( soil_type_f%var_2d ) )  THEN
                   soil_type_f%var_2d(ny+1,:) = soil_type_f%var_2d(ny,:)
                ELSE
                   soil_type_f%var_3d(:,ny+1,:) = soil_type_f%var_3d(:,ny,:)
                ENDIF
             ENDIF
             IF ( vegetation_type_f%from_file )                                &
                vegetation_type_f%var(ny+1,:) = vegetation_type_f%var(ny,:)
             IF ( water_type_f%from_file )                                     &
                water_type_f%var(ny+1,:) = water_type_f%var(ny,:)
             IF ( surface_fraction_f%from_file )                               &
                surface_fraction_f%frac(:,ny+1,:) =                            &
                                             surface_fraction_f%frac(:,ny,:)
             IF ( building_pars_f%from_file )                                  &
                building_pars_f%pars_xy(:,ny+1,:) =                            &
                                             building_pars_f%pars_xy(:,ny,:)
             IF ( albedo_pars_f%from_file )                                    &
                albedo_pars_f%pars_xy(:,ny+1,:) = albedo_pars_f%pars_xy(:,ny,:)
             IF ( pavement_pars_f%from_file )                                  &
                pavement_pars_f%pars_xy(:,ny+1,:) =                            &
                                             pavement_pars_f%pars_xy(:,ny,:)
             IF ( vegetation_pars_f%from_file )                                &
                vegetation_pars_f%pars_xy(:,ny+1,:) =                          &
                                               vegetation_pars_f%pars_xy(:,ny,:)
             IF ( water_pars_f%from_file )                                     &
                water_pars_f%pars_xy(:,ny+1,:) = water_pars_f%pars_xy(:,ny,:)
             IF ( root_area_density_lsm_f%from_file )                          &
                root_area_density_lsm_f%var(:,ny+1,:) =                        &
                                            root_area_density_lsm_f%var(:,ny,:)
             IF ( soil_pars_f%from_file )  THEN
                IF ( soil_pars_f%lod == 1 )  THEN
                   soil_pars_f%pars_xy(:,ny+1,:) = soil_pars_f%pars_xy(:,ny,:)
                ELSE
                   soil_pars_f%pars_xyz(:,:,ny+1,:) =                          &
                                              soil_pars_f%pars_xyz(:,:,ny,:)
                ENDIF
             ENDIF
             IF ( pavement_subsurface_pars_f%from_file )                       &
                pavement_subsurface_pars_f%pars_xyz(:,:,ny+1,:) =              &
                                   pavement_subsurface_pars_f%pars_xyz(:,:,ny,:)
          ENDIF
       ENDIF

       IF ( .NOT. bc_lr_cyc )  THEN
          IF ( nxl == 0 )  THEN
            IF ( albedo_type_f%from_file )                                     &
                albedo_type_f%var(:,-1) = albedo_type_f%var(:,0)
             IF ( pavement_type_f%from_file )                                  &
                pavement_type_f%var(:,-1) = pavement_type_f%var(:,0)
             IF ( soil_type_f%from_file )  THEN
                IF ( ALLOCATED( soil_type_f%var_2d ) )  THEN
                   soil_type_f%var_2d(:,-1) = soil_type_f%var_2d(:,0)
                ELSE
                   soil_type_f%var_3d(:,:,-1) = soil_type_f%var_3d(:,:,0)
                ENDIF
             ENDIF
             IF ( vegetation_type_f%from_file )                                &
                vegetation_type_f%var(:,-1) = vegetation_type_f%var(:,0)
             IF ( water_type_f%from_file )                                     &
                water_type_f%var(:,-1) = water_type_f%var(:,0)
             IF ( surface_fraction_f%from_file )                               &
                surface_fraction_f%frac(:,:,-1) = surface_fraction_f%frac(:,:,0)
             IF ( building_pars_f%from_file )                                  &
                building_pars_f%pars_xy(:,:,-1) = building_pars_f%pars_xy(:,:,0)
             IF ( albedo_pars_f%from_file )                                    &
                albedo_pars_f%pars_xy(:,:,-1) = albedo_pars_f%pars_xy(:,:,0)
             IF ( pavement_pars_f%from_file )                                  &
                pavement_pars_f%pars_xy(:,:,-1) = pavement_pars_f%pars_xy(:,:,0)
             IF ( vegetation_pars_f%from_file )                                &
                vegetation_pars_f%pars_xy(:,:,-1) =                            &
                                               vegetation_pars_f%pars_xy(:,:,0)
             IF ( water_pars_f%from_file )                                     &
                water_pars_f%pars_xy(:,:,-1) = water_pars_f%pars_xy(:,:,0)
             IF ( root_area_density_lsm_f%from_file )                          &
                root_area_density_lsm_f%var(:,:,-1) =                          &
                                            root_area_density_lsm_f%var(:,:,0)
             IF ( soil_pars_f%from_file )  THEN
                IF ( soil_pars_f%lod == 1 )  THEN
                   soil_pars_f%pars_xy(:,:,-1) = soil_pars_f%pars_xy(:,:,0)
                ELSE
                   soil_pars_f%pars_xyz(:,:,:,-1) = soil_pars_f%pars_xyz(:,:,:,0)
                ENDIF
             ENDIF
             IF ( pavement_subsurface_pars_f%from_file )                       &
                pavement_subsurface_pars_f%pars_xyz(:,:,:,-1) =                &
                                    pavement_subsurface_pars_f%pars_xyz(:,:,:,0)
          ENDIF

          IF ( nxr == nx )  THEN
             IF ( albedo_type_f%from_file )                                    &
                albedo_type_f%var(:,nx+1) = albedo_type_f%var(:,nx)
             IF ( pavement_type_f%from_file )                                  &
                pavement_type_f%var(:,nx+1) = pavement_type_f%var(:,nx)
             IF ( soil_type_f%from_file )  THEN
                IF ( ALLOCATED( soil_type_f%var_2d ) )  THEN
                   soil_type_f%var_2d(:,nx+1) = soil_type_f%var_2d(:,nx)
                ELSE
                   soil_type_f%var_3d(:,:,nx+1) = soil_type_f%var_3d(:,:,nx)
                ENDIF
             ENDIF
             IF ( vegetation_type_f%from_file )                                &
                vegetation_type_f%var(:,nx+1) = vegetation_type_f%var(:,nx)
             IF ( water_type_f%from_file )                                     &
                water_type_f%var(:,nx+1) = water_type_f%var(:,nx)
             IF ( surface_fraction_f%from_file )                               &
                surface_fraction_f%frac(:,:,nx+1) =                            &
                                             surface_fraction_f%frac(:,:,nx)
             IF ( building_pars_f%from_file )                                  &
                building_pars_f%pars_xy(:,:,nx+1) =                            &
                                             building_pars_f%pars_xy(:,:,nx)
             IF ( albedo_pars_f%from_file )                                    &
                albedo_pars_f%pars_xy(:,:,nx+1) = albedo_pars_f%pars_xy(:,:,nx)
             IF ( pavement_pars_f%from_file )                                  &
                pavement_pars_f%pars_xy(:,:,nx+1) =                            &
                                             pavement_pars_f%pars_xy(:,:,nx)
             IF ( vegetation_pars_f%from_file )                                &
                vegetation_pars_f%pars_xy(:,:,nx+1) =                          &
                                               vegetation_pars_f%pars_xy(:,:,nx)
             IF ( water_pars_f%from_file )                                     &
                water_pars_f%pars_xy(:,:,nx+1) = water_pars_f%pars_xy(:,:,nx)
             IF ( root_area_density_lsm_f%from_file )                          &
                root_area_density_lsm_f%var(:,:,nx+1) =                        &
                                            root_area_density_lsm_f%var(:,:,nx)
             IF ( soil_pars_f%from_file )  THEN
                IF ( soil_pars_f%lod == 1 )  THEN
                   soil_pars_f%pars_xy(:,:,nx+1) = soil_pars_f%pars_xy(:,:,nx)
                ELSE
                   soil_pars_f%pars_xyz(:,:,:,nx+1) =                          &
                                              soil_pars_f%pars_xyz(:,:,:,nx)
                ENDIF
             ENDIF
             IF ( pavement_subsurface_pars_f%from_file )                       &
                pavement_subsurface_pars_f%pars_xyz(:,:,:,nx+1) =              &
                                   pavement_subsurface_pars_f%pars_xyz(:,:,:,nx)
          ENDIF
       ENDIF

    END SUBROUTINE netcdf_data_input_surface_data

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads orography and building information.
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_topo

       USE control_parameters,                                                 &
           ONLY:  bc_lr_cyc, bc_ns_cyc, message_string, topography

       USE indices,                                                            &
           ONLY:  nbgp, nx, nxl, nxr, ny, nyn, nys, nzb, nzt


       IMPLICIT NONE

       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names  !< variable names in static input file


       INTEGER(iwp) ::  i             !< running index along x-direction
       INTEGER(iwp) ::  ii            !< running index for IO blocks
       INTEGER(iwp) ::  id_topo       !< NetCDF id of topograhy input file
       INTEGER(iwp) ::  j             !< running index along y-direction
       INTEGER(iwp) ::  k             !< running index along z-direction
       INTEGER(iwp) ::  num_vars      !< number of variables in netcdf input file
       INTEGER(iwp) ::  skip_n_rows   !< counting variable to skip rows while reading topography file

       INTEGER(iwp), DIMENSION(nys-nbgp:nyn+nbgp,nxl-nbgp:nxr+nbgp) ::  var_exchange_int !< dummy variables used to exchange 32-bit Integer arrays

       REAL(wp) ::  dum           !< dummy variable to skip columns while reading topography file
!
!--    CPU measurement
       CALL cpu_log( log_point_s(83), 'NetCDF/ASCII input topo', 'start' )

!
!--    Input via palm-input data standard
       IF ( input_pids_static )  THEN
#if defined ( __netcdf )
!
!--       Open file in read-only mode
          CALL open_read_file( TRIM( input_file_static ) //                    &
                               TRIM( coupling_char ), id_topo )
!
!--       At first, inquire all variable names.
!--       This will be used to check whether an  input variable exist
!--       or not.
          CALL inquire_num_variables( id_topo, num_vars )
!
!--       Allocate memory to store variable names and inquire them.
          ALLOCATE( var_names(1:num_vars) )
          CALL inquire_variable_names( id_topo, var_names )
!
!--       Read x, y - dimensions. Only required for consistency checks.
          CALL get_dimension_length( id_topo, dim_static%nx, 'x' )
          CALL get_dimension_length( id_topo, dim_static%ny, 'y' )
          ALLOCATE( dim_static%x(0:dim_static%nx-1) )
          ALLOCATE( dim_static%y(0:dim_static%ny-1) )
          CALL get_variable( id_topo, 'x', dim_static%x )
          CALL get_variable( id_topo, 'y', dim_static%y )
!
!--       Terrain height. First, get variable-related _FillValue attribute
          IF ( check_existence( var_names, 'zt' ) )  THEN
             terrain_height_f%from_file = .TRUE.
             CALL get_attribute( id_topo, char_fill, terrain_height_f%fill,    &
                                 .FALSE., 'zt' )
!
!--          Input 2D terrain height.
             ALLOCATE ( terrain_height_f%var(nys:nyn,nxl:nxr)  )
             
             CALL get_variable( id_topo, 'zt', terrain_height_f%var,           &
                                nxl, nxr, nys, nyn )

          ELSE
             terrain_height_f%from_file = .FALSE.
          ENDIF

!
!--       Read building height. First, read its _FillValue attribute,
!--       as well as lod attribute
          buildings_f%from_file = .FALSE.
          IF ( check_existence( var_names, 'buildings_2d' ) )  THEN
             buildings_f%from_file = .TRUE.
             CALL get_attribute( id_topo, char_lod, buildings_f%lod,           &
                                 .FALSE., 'buildings_2d' )

             CALL get_attribute( id_topo, char_fill, buildings_f%fill1,        &
                                 .FALSE., 'buildings_2d' )

!
!--          Read 2D buildings
             IF ( buildings_f%lod == 1 )  THEN
                ALLOCATE ( buildings_f%var_2d(nys:nyn,nxl:nxr) )

                CALL get_variable( id_topo, 'buildings_2d',                    &
                                   buildings_f%var_2d,                         &
                                   nxl, nxr, nys, nyn )
             ELSE
                message_string = 'NetCDF attribute lod ' //                    &
                                 '(level of detail) is not set ' //            &
                                 'properly for buildings_2d.'
                CALL message( 'netcdf_data_input_mod', 'NDI000',               &
                               1, 2, 0, 6, 0 )
             ENDIF
          ENDIF
!
!--       If available, also read 3D building information. If both are
!--       available, use 3D information.
          IF ( check_existence( var_names, 'buildings_3d' ) )  THEN
             buildings_f%from_file = .TRUE.
             CALL get_attribute( id_topo, char_lod, buildings_f%lod,           &
                                 .FALSE., 'buildings_3d' )      

             CALL get_attribute( id_topo, char_fill, buildings_f%fill2,        &
                                 .FALSE., 'buildings_3d' )

             CALL get_dimension_length( id_topo, buildings_f%nz, 'z' )
!
!--          Read 3D buildings
             IF ( buildings_f%lod == 2 )  THEN
                ALLOCATE( buildings_f%z(nzb:buildings_f%nz-1) )
                CALL get_variable( id_topo, 'z', buildings_f%z )

                ALLOCATE( buildings_f%var_3d(nzb:buildings_f%nz-1,             &
                                             nys:nyn,nxl:nxr) )
                buildings_f%var_3d = 0
                
                CALL get_variable( id_topo, 'buildings_3d',                    &
                                   buildings_f%var_3d,                         &
                                   nxl, nxr, nys, nyn, 0, buildings_f%nz-1 )
             ELSE
                message_string = 'NetCDF attribute lod ' //                    &
                                 '(level of detail) is not set ' //            &
                                 'properly for buildings_3d.'
                CALL message( 'netcdf_data_input_mod', 'NDI001',               &
                               1, 2, 0, 6, 0 )
             ENDIF
          ENDIF
!
!--       Read building IDs and its FillValue attribute. Further required
!--       for mapping buildings on top of orography.
          IF ( check_existence( var_names, 'building_id' ) )  THEN
             building_id_f%from_file = .TRUE.
             CALL get_attribute( id_topo, char_fill,                           &
                                 building_id_f%fill, .FALSE.,                  &
                                 'building_id' )

             ALLOCATE ( building_id_f%var(nys:nyn,nxl:nxr) )
             
             CALL get_variable( id_topo, 'building_id', building_id_f%var,     &
                                nxl, nxr, nys, nyn )
          ELSE
             building_id_f%from_file = .FALSE.
          ENDIF
!
!--       Read building_type and required attributes.
          IF ( check_existence( var_names, 'building_type' ) )  THEN
             building_type_f%from_file = .TRUE.
             CALL get_attribute( id_topo, char_fill,                           &
                                 building_type_f%fill, .FALSE.,                &
                                 'building_type' )

             ALLOCATE ( building_type_f%var(nys:nyn,nxl:nxr) )

             CALL get_variable( id_topo, 'building_type', building_type_f%var, &
                                nxl, nxr, nys, nyn )

          ELSE
             building_type_f%from_file = .FALSE.
          ENDIF
!
!--       Close topography input file
          CALL close_input_file( id_topo )
#else
          CONTINUE
#endif
!
!--    ASCII input
       ELSEIF ( TRIM( topography ) == 'read_from_file' )  THEN
             
          DO  ii = 0, io_blocks-1
             IF ( ii == io_group )  THEN

                OPEN( 90, FILE='TOPOGRAPHY_DATA'//TRIM( coupling_char ),       &
                      STATUS='OLD', FORM='FORMATTED', ERR=10 )
!
!--             Read topography PE-wise. Rows are read from nyn to nys, columns
!--             are read from nxl to nxr. At first, ny-nyn rows need to be skipped.
                skip_n_rows = 0
                DO WHILE ( skip_n_rows < ny - nyn )
                   READ( 90, * )
                   skip_n_rows = skip_n_rows + 1
                ENDDO
!
!--             Read data from nyn to nys and nxl to nxr. Therefore, skip
!--             column until nxl-1 is reached
                ALLOCATE ( buildings_f%var_2d(nys:nyn,nxl:nxr) )
                DO  j = nyn, nys, -1
                   READ( 90, *, ERR=11, END=11 )                               &
                                   ( dum, i = 0, nxl-1 ),                      &
                                   ( buildings_f%var_2d(j,i), i = nxl, nxr )
                ENDDO

                GOTO 12

 10             message_string = 'file TOPOGRAPHY_DATA'//                      &
                                 TRIM( coupling_char )// ' does not exist'
                CALL message( 'netcdf_data_input_mod', 'PA0208', 1, 2, 0, 6, 0 )

 11             message_string = 'errors in file TOPOGRAPHY_DATA'//            &
                                 TRIM( coupling_char )
                CALL message( 'netcdf_data_input_mod', 'PA0209', 2, 2, 0, 6, 0 )

 12             CLOSE( 90 )
                buildings_f%from_file = .TRUE.

             ENDIF
#if defined( __parallel )
             CALL MPI_BARRIER( comm2d, ierr )
#endif
          ENDDO

       ENDIF
!
!--    End of CPU measurement
       CALL cpu_log( log_point_s(83), 'NetCDF/ASCII input topo', 'stop' )
!
!--    Check for minimum requirement to setup building topography. If buildings
!--    are provided, also an ID and a type are required.
!--    Note, doing this check in check_parameters
!--    will be too late (data will be used for grid inititialization before).
       IF ( input_pids_static )  THEN
          IF ( buildings_f%from_file  .AND.                                    &
               .NOT. building_id_f%from_file )  THEN
             message_string = 'If building heigths are prescribed in ' //      &
                              'static input file, also an ID is required.'
             CALL message( 'netcdf_data_input_mod', 'NDI002', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
!
!--    In case no terrain height is provided by static input file, allocate
!--    array nevertheless and set terrain height to 0, which simplifies
!--    topography initialization.
       IF ( .NOT. terrain_height_f%from_file )  THEN
          ALLOCATE ( terrain_height_f%var(nys:nyn,nxl:nxr) )
          terrain_height_f%var = 0.0_wp
       ENDIF
!
!--    Finally, exchange 1 ghost point for building ID and type.
!--    In case of non-cyclic boundary conditions set Neumann conditions at the
!--    lateral boundaries.
       IF ( building_id_f%from_file )  THEN
          var_exchange_int                  = building_id_f%fill
          var_exchange_int(nys:nyn,nxl:nxr) = building_id_f%var(nys:nyn,nxl:nxr)
          CALL exchange_horiz_2d_int( var_exchange_int, nys, nyn, nxl, nxr, nbgp )
          DEALLOCATE( building_id_f%var )
          ALLOCATE( building_id_f%var(nys-nbgp:nyn+nbgp,nxl-nbgp:nxr+nbgp) )
          building_id_f%var = var_exchange_int

          IF ( .NOT. bc_ns_cyc )  THEN
             IF ( nys == 0  )  building_id_f%var(-1,:)   = building_id_f%var(0,:)
             IF ( nyn == ny )  building_id_f%var(ny+1,:) = building_id_f%var(ny,:)
          ENDIF
          IF ( .NOT. bc_lr_cyc )  THEN
             IF ( nxl == 0  )  building_id_f%var(:,-1)   = building_id_f%var(:,0)
             IF ( nxr == nx )  building_id_f%var(:,nx+1) = building_id_f%var(:,nx)
          ENDIF
       ENDIF

       IF ( building_type_f%from_file )  THEN
          var_exchange_int                  = INT( building_type_f%fill, KIND = 4 )
          var_exchange_int(nys:nyn,nxl:nxr) =                                  &
                          INT( building_type_f%var(nys:nyn,nxl:nxr), KIND = 4 )
          CALL exchange_horiz_2d_int( var_exchange_int, nys, nyn, nxl, nxr, nbgp )
          DEALLOCATE( building_type_f%var )
          ALLOCATE( building_type_f%var(nys-nbgp:nyn+nbgp,nxl-nbgp:nxr+nbgp) )
          building_type_f%var = INT( var_exchange_int, KIND = 1 )

          IF ( .NOT. bc_ns_cyc )  THEN
             IF ( nys == 0  )  building_type_f%var(-1,:)   = building_type_f%var(0,:)
             IF ( nyn == ny )  building_type_f%var(ny+1,:) = building_type_f%var(ny,:)
          ENDIF
          IF ( .NOT. bc_lr_cyc )  THEN
             IF ( nxl == 0  )  building_type_f%var(:,-1)   = building_type_f%var(:,0)
             IF ( nxr == nx )  building_type_f%var(:,nx+1) = building_type_f%var(:,nx)
          ENDIF
       ENDIF

    END SUBROUTINE netcdf_data_input_topo

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads initialization data of u, v, w, pt, q, geostrophic wind components,
!> as well as soil moisture and soil temperature, derived from larger-scale
!> model (COSMO) by Inifor.
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_init_3d

       USE arrays_3d,                                                          &
           ONLY:  q, pt, u, v, w

       USE control_parameters,                                                 &
           ONLY:  bc_lr_cyc, bc_ns_cyc, forcing, humidity, land_surface,       &
                  message_string, neutral, surface_pressure

       USE indices,                                                            &
           ONLY:  nx, nxl, nxlu, nxr, ny, nyn, nys, nysv, nzb, nz, nzt

       IMPLICIT NONE

       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names

       LOGICAL      ::  dynamic_3d = .TRUE. !< flag indicating that 3D data is read from dynamic file
       
       INTEGER(iwp) ::  i          !< running index along x-direction
       INTEGER(iwp) ::  id_dynamic !< NetCDF id of dynamic input file
       INTEGER(iwp) ::  j          !< running index along y-direction
       INTEGER(iwp) ::  k          !< running index along z-direction
       INTEGER(iwp) ::  num_vars   !< number of variables in netcdf input file

       LOGICAL      ::  check_passed !< flag indicating if a check passed

!
!--    Skip routine if no input file with dynamic input data is available.
       IF ( .NOT. input_pids_dynamic )  RETURN
!
!--    Please note, Inifor is designed to provide initial data for u and v for
!--    the prognostic grid points in case of lateral Dirichlet conditions.
!--    This means that Inifor provides data from nxlu:nxr (for u) and
!--    from nysv:nyn (for v) at the left and south domain boundary, respectively.
!--    However, as work-around for the moment, PALM will run with cyclic
!--    conditions and will be initialized with data provided by Inifor
!--    boundaries in case of Dirichlet.
!--    Hence, simply set set nxlu/nysv to 1 (will be reset to its original value
!--    at the end of this routine.
       IF ( bc_lr_cyc  .AND.  nxl == 0 )  nxlu = 1
       IF ( bc_ns_cyc  .AND.  nys == 0 )  nysv = 1

!
!--    CPU measurement
       CALL cpu_log( log_point_s(85), 'NetCDF input init', 'start' )

#if defined ( __netcdf )
!
!--    Open file in read-only mode
       CALL open_read_file( TRIM( input_file_dynamic ) //                      &
                            TRIM( coupling_char ), id_dynamic )

!
!--    At first, inquire all variable names.
       CALL inquire_num_variables( id_dynamic, num_vars )
!
!--    Allocate memory to store variable names.
       ALLOCATE( var_names(1:num_vars) )
       CALL inquire_variable_names( id_dynamic, var_names )
!
!--    Read vertical dimension of scalar und w grid. Will be used for
!--    inter- and extrapolation in case of stretched numeric grid.
!--    This will be removed when Inifor is able to handle stretched grids.
       CALL get_dimension_length( id_dynamic, init_3d%nzu, 'z'     )
       CALL get_dimension_length( id_dynamic, init_3d%nzw, 'zw'    )
       CALL get_dimension_length( id_dynamic, init_3d%nzs, 'depth' )
!
!--    Read also the horizontal dimensions. These are used just used fo
!--    checking the compatibility with the PALM grid before reading.
       CALL get_dimension_length( id_dynamic, init_3d%nx,  'x'  )
       CALL get_dimension_length( id_dynamic, init_3d%nxu, 'xu' )
       CALL get_dimension_length( id_dynamic, init_3d%ny,  'y'  )
       CALL get_dimension_length( id_dynamic, init_3d%nyv, 'yv' )

!
!--    Check for correct horizontal and vertical dimension. Please note,
!--    checks are performed directly here and not called from
!--    check_parameters as some varialbes are still not allocated there.
!--    Moreover, please note, u- and v-grid has 1 grid point less on
!--    Inifor grid.
       IF ( init_3d%nx-1 /= nx  .OR.  init_3d%nxu-1 /= nx - 1  .OR.            &
            init_3d%ny-1 /= ny  .OR.  init_3d%nyv-1 /= ny - 1 )  THEN
          message_string = 'Number of inifor horizontal grid points  '//       &
                           'does not match the number of numeric grid '//      &
                           'points.'
          CALL message( 'netcdf_data_input_mod', 'NDI003', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( init_3d%nzu-1 /= nz )  THEN
          message_string = 'Number of inifor vertical grid points ' //         &
                           'does not match the number of numeric grid '//      &
                           'points.'
          CALL message( 'netcdf_data_input_mod', 'NDI003', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Read vertical dimensions. Later, these are required for eventual
!--    inter- and extrapolations of the initialization data.
       IF ( check_existence( var_names, 'z' ) )  THEN
          ALLOCATE( init_3d%zu_atmos(1:init_3d%nzu) )
          CALL get_variable( id_dynamic, 'z', init_3d%zu_atmos )
       ENDIF
       IF ( check_existence( var_names, 'zw' ) )  THEN
          ALLOCATE( init_3d%zw_atmos(1:init_3d%nzw) )
          CALL get_variable( id_dynamic, 'zw', init_3d%zw_atmos )
       ENDIF
       IF ( check_existence( var_names, 'depth' ) )  THEN
          ALLOCATE( init_3d%z_soil(1:init_3d%nzs) )
          CALL get_variable( id_dynamic, 'depth', init_3d%z_soil )
       ENDIF
!
!--    Read initial geostrophic wind components at t = 0 (index 1 in file).
!        IF ( check_existence( var_names, 'tend_ug' ) )  THEN
       IF ( check_existence( var_names, 'ls_forcing_ug' ) )  THEN
          ALLOCATE( init_3d%ug_init(nzb:nzt+1) )
!           CALL get_variable_pr( id_dynamic, 'tend_ug', 1,                      &
!                                 init_3d%ug_init )
          CALL get_variable_pr( id_dynamic, 'ls_forcing_ug', 1,                &
                                init_3d%ug_init )
          init_3d%from_file_ug = .TRUE.
       ELSE
          init_3d%from_file_ug = .FALSE.
       ENDIF
!        IF ( check_existence( var_names, 'tend_vg' ) )  THEN
       IF ( check_existence( var_names, 'ls_forcing_vg' ) )  THEN
          ALLOCATE( init_3d%vg_init(nzb:nzt+1) )
!           CALL get_variable_pr( id_dynamic, 'tend_vg', 1,                      &
!                                 init_3d%vg_init )
          CALL get_variable_pr( id_dynamic, 'ls_forcing_vg', 1,                &
                                init_3d%vg_init )
          init_3d%from_file_vg = .TRUE.
       ELSE
          init_3d%from_file_vg = .FALSE.
       ENDIF
!
!--    Read inital 3D data of u, v, w, pt and q,
!--    derived from COSMO model. Read PE-wise yz-slices.
!--    Please note, the u-, v- and w-component are defined on different
!--    grids with one element less in the x-, y-,
!--    and z-direction, respectively. Hence, reading is subdivided
!--    into separate loops.  
!--    Read u-component
       IF ( check_existence( var_names, 'init_u' ) )  THEN
!
!--       Read attributes for the fill value and level-of-detail
          CALL get_attribute( id_dynamic, char_fill, init_3d%fill_u,           &
                              .FALSE., 'init_u' )
          CALL get_attribute( id_dynamic, char_lod, init_3d%lod_u,             &
                              .FALSE., 'init_u' )
!
!--       level-of-detail 1 - read initialization profile
          IF ( init_3d%lod_u == 1 )  THEN
             ALLOCATE( init_3d%u_init(nzb:nzt+1) )
             init_3d%u_init = 0.0_wp

             CALL get_variable( id_dynamic, 'init_u',                          &
                                init_3d%u_init(nzb+1:nzt+1) )
!
!--       level-of-detail 2 - read 3D initialization data
          ELSEIF ( init_3d%lod_u == 2 )  THEN

             CALL get_variable( id_dynamic, 'init_u',                          &
                                u(nzb+1:nzt+1,nys:nyn,nxlu:nxr),               &
                                nxlu, nys+1, nzb+1,                            &
                                nxr-nxlu+1, nyn-nys+1, init_3d%nzu,            &
                                dynamic_3d )
          ENDIF
          init_3d%from_file_u = .TRUE.
       ENDIF
!
!--    Read v-component
       IF ( check_existence( var_names, 'init_v' ) )  THEN
!
!--       Read attributes for the fill value and level-of-detail
          CALL get_attribute( id_dynamic, char_fill, init_3d%fill_v,           &
                              .FALSE., 'init_v' )
          CALL get_attribute( id_dynamic, char_lod, init_3d%lod_v,             &
                              .FALSE., 'init_v' )
!
!--       level-of-detail 1 - read initialization profile
          IF ( init_3d%lod_v == 1 )  THEN
             ALLOCATE( init_3d%v_init(nzb:nzt+1) )
             init_3d%v_init = 0.0_wp

             CALL get_variable( id_dynamic, 'init_v',                          &
                                init_3d%v_init(nzb+1:nzt+1) )

!
!--       level-of-detail 2 - read 3D initialization data
          ELSEIF ( init_3d%lod_v == 2 )  THEN

             CALL get_variable( id_dynamic, 'init_v',                          &
                                v(nzb+1:nzt+1,nysv:nyn,nxl:nxr),               &
                                nxl+1, nysv, nzb+1,                            &
                                nxr-nxl+1, nyn-nysv+1, init_3d%nzu,            &
                                dynamic_3d )
             
          ENDIF
          init_3d%from_file_v = .TRUE.
       ENDIF
!
!--    Read w-component
       IF ( check_existence( var_names, 'init_w' ) )  THEN
!
!--       Read attributes for the fill value and level-of-detail
          CALL get_attribute( id_dynamic, char_fill, init_3d%fill_w,           &
                              .FALSE., 'init_w' )
          CALL get_attribute( id_dynamic, char_lod, init_3d%lod_w,             &
                              .FALSE., 'init_w' )
!
!--       level-of-detail 1 - read initialization profile
          IF ( init_3d%lod_w == 1 )  THEN
             ALLOCATE( init_3d%w_init(nzb:nzt+1) )
             init_3d%w_init = 0.0_wp

             CALL get_variable( id_dynamic, 'init_w',                          &
                                init_3d%w_init(nzb+1:nzt) )

!
!--       level-of-detail 2 - read 3D initialization data
          ELSEIF ( init_3d%lod_w == 2 )  THEN

             CALL get_variable( id_dynamic, 'init_w',                           &
                                w(nzb+1:nzt,nys:nyn,nxl:nxr),                   &
                                nxl+1, nys+1, nzb+1,                            &
                                nxr-nxl+1, nyn-nys+1, init_3d%nzw,              &
                                dynamic_3d )

          ENDIF
          init_3d%from_file_w = .TRUE.
       ENDIF
!
!--    Read potential temperature
       IF ( .NOT. neutral )  THEN
          IF ( check_existence( var_names, 'init_pt' ) )  THEN
!
!--          Read attributes for the fill value and level-of-detail
             CALL get_attribute( id_dynamic, char_fill, init_3d%fill_pt,       &
                                 .FALSE., 'init_pt' )
             CALL get_attribute( id_dynamic, char_lod, init_3d%lod_pt,         &
                                 .FALSE., 'init_pt' )
!
!--          level-of-detail 1 - read initialization profile
             IF ( init_3d%lod_pt == 1 )  THEN
                ALLOCATE( init_3d%pt_init(nzb:nzt+1) )

                CALL get_variable( id_dynamic, 'init_pt',                      &
                                   init_3d%pt_init(nzb+1:nzt+1) )
!
!--             Set Neumann surface boundary condition for initial profil
                init_3d%pt_init(nzb) = init_3d%pt_init(nzb+1)
!
!--          level-of-detail 2 - read 3D initialization data
             ELSEIF ( init_3d%lod_pt == 2 )  THEN

                CALL get_variable( id_dynamic, 'init_pt',                      &
                                   pt(nzb+1:nzt+1,nys:nyn,nxl:nxr),            &
                                   nxl+1, nys+1, nzb+1,                        &
                                   nxr-nxl+1, nyn-nys+1, init_3d%nzu,          &
                                   dynamic_3d )


             ENDIF
             init_3d%from_file_pt = .TRUE.
          ENDIF
       ENDIF
!
!--    Read mixing ratio
       IF ( humidity )  THEN
          IF ( check_existence( var_names, 'init_qv' ) )  THEN
!
!--          Read attributes for the fill value and level-of-detail
             CALL get_attribute( id_dynamic, char_fill, init_3d%fill_q,        &
                                 .FALSE., 'init_qv' )
             CALL get_attribute( id_dynamic, char_lod, init_3d%lod_q,          &
                                 .FALSE., 'init_qv' )
!
!--          level-of-detail 1 - read initialization profile
             IF ( init_3d%lod_q == 1 )  THEN
                ALLOCATE( init_3d%q_init(nzb:nzt+1) )

                CALL get_variable( id_dynamic, 'init_qv',                      &
                                   init_3d%q_init(nzb+1:nzt+1) )
!
!--             Set Neumann surface boundary condition for initial profil
                init_3d%q_init(nzb) = init_3d%q_init(nzb+1)

!
!--          level-of-detail 2 - read 3D initialization data
             ELSEIF ( init_3d%lod_q == 2 )  THEN
             
                CALL get_variable( id_dynamic, 'init_qv',                      &
                                   q(nzb+1:nzt+1,nys:nyn,nxl:nxr),             &
                                   nxl+1, nys+1, nzb+1,                        &
                                   nxr-nxl+1, nyn-nys+1, init_3d%nzu,          &
                                   dynamic_3d )



             ENDIF
             init_3d%from_file_q = .TRUE.
          ENDIF
       ENDIF
!
!--    Read soil moisture
       IF ( land_surface )  THEN

          IF ( check_existence( var_names, 'init_soil_m' ) )  THEN
!
!--          Read attributes for the fill value and level-of-detail
             CALL get_attribute( id_dynamic, char_fill,                        &
                                 init_3d%fill_msoil,                           &
                                 .FALSE., 'init_soil_m' )
             CALL get_attribute( id_dynamic, char_lod,                         &
                                 init_3d%lod_msoil,                            &
                                 .FALSE., 'init_soil_m' )
!
!--          level-of-detail 1 - read initialization profile
             IF ( init_3d%lod_msoil == 1 )  THEN
                ALLOCATE( init_3d%msoil_init(0:init_3d%nzs-1) )

                CALL get_variable( id_dynamic, 'init_soil_m',                  &
                                   init_3d%msoil_init(0:init_3d%nzs-1) )
!
!--          level-of-detail 2 - read 3D initialization data
             ELSEIF ( init_3d%lod_msoil == 2 )  THEN
                ALLOCATE ( init_3d%msoil(0:init_3d%nzs-1,nys:nyn,nxl:nxr) )

               CALL get_variable( id_dynamic, 'init_soil_m',                   &   
                                  init_3d%msoil(0:init_3d%nzs-1,nys:nyn,nxl:nxr),&
                                  nxl, nxr, nys, nyn, 0, init_3d%nzs-1 )

             ENDIF
             init_3d%from_file_msoil = .TRUE.
          ENDIF
!
!--       Read soil temperature
          IF ( check_existence( var_names, 'init_soil_t' ) )  THEN
!
!--          Read attributes for the fill value and level-of-detail
             CALL get_attribute( id_dynamic, char_fill,                        &
                                 init_3d%fill_tsoil,                           &
                                 .FALSE., 'init_soil_t' )
             CALL get_attribute( id_dynamic, char_lod,                         &
                                 init_3d%lod_tsoil,                            &
                                 .FALSE., 'init_soil_t' )
!
!--          level-of-detail 1 - read initialization profile
             IF ( init_3d%lod_tsoil == 1 )  THEN
                ALLOCATE( init_3d%tsoil_init(0:init_3d%nzs-1) )

                CALL get_variable( id_dynamic, 'init_soil_t',                  &
                                   init_3d%tsoil_init(0:init_3d%nzs-1) )

!
!--          level-of-detail 2 - read 3D initialization data
             ELSEIF ( init_3d%lod_tsoil == 2 )  THEN
                ALLOCATE ( init_3d%tsoil(0:init_3d%nzs-1,nys:nyn,nxl:nxr) )
                
                CALL get_variable( id_dynamic, 'init_soil_t',                  &   
                                  init_3d%tsoil(0:init_3d%nzs-1,nys:nyn,nxl:nxr),&
                                  nxl, nxr, nys, nyn, 0, init_3d%nzs-1 )
             ENDIF
             init_3d%from_file_tsoil = .TRUE.
          ENDIF
       ENDIF
!
!--    Close input file
       CALL close_input_file( id_dynamic )
#endif
!
!--    End of CPU measurement
       CALL cpu_log( log_point_s(85), 'NetCDF input init', 'stop' )
!
!--    Finally, check if the input data has any fill values. Please note,
!--    checks depend on the LOD of the input data.
       IF ( init_3d%from_file_u )  THEN
          check_passed = .TRUE.
          IF ( init_3d%lod_u == 1 )  THEN
             IF ( ANY( init_3d%u_init(nzb+1:nzt+1) == init_3d%fill_u ) )       &
                check_passed = .FALSE.
          ELSEIF ( init_3d%lod_u == 2 )  THEN
             IF ( ANY( u(nzb+1:nzt+1,nys:nyn,nxlu:nxr) == init_3d%fill_u ) )   &
                check_passed = .FALSE.
          ENDIF
          IF ( .NOT. check_passed )  THEN
             message_string = 'NetCDF input for u_init must not contain ' //   &
                              'any _FillValues'
             CALL message( 'netcdf_data_input_mod', 'NDI004', 2, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       IF ( init_3d%from_file_v )  THEN
          check_passed = .TRUE.
          IF ( init_3d%lod_v == 1 )  THEN
             IF ( ANY( init_3d%v_init(nzb+1:nzt+1) == init_3d%fill_v ) )       &
                check_passed = .FALSE.
          ELSEIF ( init_3d%lod_v == 2 )  THEN
             IF ( ANY( v(nzb+1:nzt+1,nysv:nyn,nxl:nxr) == init_3d%fill_v ) )   &
                check_passed = .FALSE.
          ENDIF
          IF ( .NOT. check_passed )  THEN
             message_string = 'NetCDF input for v_init must not contain ' //   &
                              'any _FillValues'
             CALL message( 'netcdf_data_input_mod', 'NDI005', 2, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       IF ( init_3d%from_file_w )  THEN
          check_passed = .TRUE.
          IF ( init_3d%lod_w == 1 )  THEN
             IF ( ANY( init_3d%w_init(nzb+1:nzt) == init_3d%fill_w ) )         &
                check_passed = .FALSE.
          ELSEIF ( init_3d%lod_w == 2 )  THEN
             IF ( ANY( w(nzb+1:nzt,nys:nyn,nxl:nxr) == init_3d%fill_w ) )      &
                check_passed = .FALSE.
          ENDIF
          IF ( .NOT. check_passed )  THEN
             message_string = 'NetCDF input for w_init must not contain ' //   &
                              'any _FillValues'
             CALL message( 'netcdf_data_input_mod', 'NDI006', 2, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       IF ( init_3d%from_file_pt )  THEN
          check_passed = .TRUE.
          IF ( init_3d%lod_pt == 1 )  THEN
             IF ( ANY( init_3d%pt_init(nzb+1:nzt+1) == init_3d%fill_pt ) )     &
                check_passed = .FALSE.
          ELSEIF ( init_3d%lod_pt == 2 )  THEN
             IF ( ANY( pt(nzb+1:nzt+1,nys:nyn,nxl:nxr) == init_3d%fill_pt ) )  &
                check_passed = .FALSE.
          ENDIF
          IF ( .NOT. check_passed )  THEN
             message_string = 'NetCDF input for pt_init must not contain ' //  &
                              'any _FillValues'
             CALL message( 'netcdf_data_input_mod', 'NDI007', 2, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       IF ( init_3d%from_file_q )  THEN
          check_passed = .TRUE.
          IF ( init_3d%lod_q == 1 )  THEN
             IF ( ANY( init_3d%q_init(nzb+1:nzt+1) == init_3d%fill_q ) )       &
                check_passed = .FALSE.
          ELSEIF ( init_3d%lod_q == 2 )  THEN
             IF ( ANY( q(nzb+1:nzt+1,nys:nyn,nxl:nxr) == init_3d%fill_q ) )    &
                check_passed = .FALSE.
          ENDIF
          IF ( .NOT. check_passed )  THEN
             message_string = 'NetCDF input for q_init must not contain ' //   &
                              'any _FillValues'
             CALL message( 'netcdf_data_input_mod', 'NDI008', 2, 2, 0, 6, 0 )
          ENDIF
       ENDIF
!
!--    Workaround for cyclic conditions. Please see above for further explanation.
       IF ( bc_lr_cyc  .AND.  nxl == 0 )  nxlu = nxl
       IF ( bc_ns_cyc  .AND.  nys == 0 )  nysv = nys

    END SUBROUTINE netcdf_data_input_init_3d

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads data at lateral and top boundaries derived from larger-scale model
!> (COSMO) by Inifor.
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_lsf

       USE control_parameters,                                                 &
           ONLY:  bc_lr_cyc, bc_ns_cyc, force_bound_l, force_bound_n,          &
                  force_bound_r, force_bound_s,                                &
                  forcing, humidity, message_string, neutral, simulated_time


       USE indices,                                                            &
           ONLY:  nxl, nxlu, nxr, nyn, nys, nysv, nzb, nzt

       IMPLICIT NONE

       LOGICAL      ::  dynamic_3d = .TRUE. !< flag indicating that 3D data is read from dynamic file
       
       INTEGER(iwp) ::  i          !< running index along x-direction
       INTEGER(iwp) ::  id_dynamic !< NetCDF id of dynamic input file
       INTEGER(iwp) ::  j          !< running index along y-direction
       INTEGER(iwp) ::  k          !< running index along z-direction
       INTEGER(iwp) ::  num_vars   !< number of variables in netcdf input file
       INTEGER(iwp) ::  t          !< running index time dimension

       REAL(wp) ::  dum           !< dummy variable to skip columns while reading topography file

       force%from_file = MERGE( .TRUE., .FALSE., input_pids_dynamic )
!
!--    Skip input if no forcing from larger-scale models is applied.
       IF ( .NOT. forcing )  RETURN

!
!--    CPU measurement
       CALL cpu_log( log_point_s(86), 'NetCDF input forcing', 'start' )

#if defined ( __netcdf )
!
!--    Open file in read-only mode
       CALL open_read_file( TRIM( input_file_dynamic ) //                      &
                            TRIM( coupling_char ), id_dynamic )
!
!--    Initialize INIFOR forcing.
       IF ( .NOT. force%init )  THEN
!
!--       At first, inquire all variable names.
          CALL inquire_num_variables( id_dynamic, num_vars )
!
!--       Allocate memory to store variable names.
          ALLOCATE( force%var_names(1:num_vars) )
          CALL inquire_variable_names( id_dynamic, force%var_names )
!
!--       Read time dimension, allocate memory and finally read time array
          CALL get_dimension_length( id_dynamic, force%nt, 'time' )

          IF ( check_existence( force%var_names, 'time' ) )  THEN
             ALLOCATE( force%time(0:force%nt-1) )
             CALL get_variable( id_dynamic, 'time', force%time )
          ENDIF
!
!--       Read vertical dimension of scalar und w grid
          CALL get_dimension_length( id_dynamic, force%nzu, 'z' )
          CALL get_dimension_length( id_dynamic, force%nzw, 'zw' )

          IF ( check_existence( force%var_names, 'z' ) )  THEN
             ALLOCATE( force%zu_atmos(1:force%nzu) )
             CALL get_variable( id_dynamic, 'z', force%zu_atmos )
          ENDIF
          IF ( check_existence( force%var_names, 'zw' ) )  THEN
             ALLOCATE( force%zw_atmos(1:force%nzw) )
             CALL get_variable( id_dynamic, 'zw', force%zw_atmos )
          ENDIF

!
!--       Read surface pressure
          IF ( check_existence( force%var_names,                               &
                            'surface_forcing_surface_pressure' ) )  THEN
             ALLOCATE( force%surface_pressure(0:force%nt-1) )
             CALL get_variable( id_dynamic,                                    &
                                'surface_forcing_surface_pressure',            &
                                force%surface_pressure )
          ENDIF
!
!--       Set control flag to indicate that initialization is already done
          force%init = .TRUE.

       ENDIF

!
!--    Obtain time index for current input starting at 0.
!--    @todo: At the moment time, in INIFOR and simulated time correspond
!--           to each other. If required, adjust to daytime.
       force%tind = MINLOC( ABS( force%time - simulated_time ), DIM = 1 )      &
                    - 1
       force%tind_p = force%tind + 1       
!
!--    Read geostrophic wind components. In case of forcing, this is only
!--    required if cyclic boundary conditions are applied.
       IF ( bc_lr_cyc  .AND.  bc_ns_cyc )  THEN
          DO  t = force%tind, force%tind_p
!              CALL get_variable_pr( id_dynamic, 'tend_ug', t+1,           &
!                                    force%ug(t-force%tind,:) )
!              CALL get_variable_pr( id_dynamic, 'tend_vg', t+1,           &
!                                    force%ug(t-force%tind,:) )
             CALL get_variable_pr( id_dynamic, 'ls_forcing_ug', t+1,           &
                                   force%ug(t-force%tind,:) )
             CALL get_variable_pr( id_dynamic, 'ls_forcing_vg', t+1,           &
                                   force%ug(t-force%tind,:) )
          ENDDO
       ENDIF
!
!--    Read data at lateral and top boundaries. Please note, at left and
!--    right domain boundary, yz-layers are read for u, v, w, pt and q.
!--    For the v-component, the data starts at nysv, while for the other
!--    quantities the data starts at nys. This is equivalent at the north
!--    and south domain boundary for the u-component.
       IF ( force_bound_l )  THEN
          CALL get_variable( id_dynamic, 'ls_forcing_left_u',                  &
                           force%u_left(0:1,nzb+1:nzt+1,nys:nyn),              &
                           nys+1, nzb+1, force%tind+1,                         &
                           nyn-nys+1, force%nzu, 2, dynamic_3d )
          
          CALL get_variable( id_dynamic, 'ls_forcing_left_v',                  &
                           force%v_left(0:1,nzb+1:nzt+1,nysv:nyn),             &
                           nysv, nzb+1, force%tind+1,                          &
                           nyn-nysv+1, force%nzu, 2, dynamic_3d )

          CALL get_variable( id_dynamic, 'ls_forcing_left_w',                  &
                           force%w_left(0:1,nzb+1:nzt,nys:nyn),                &
                           nys+1, nzb+1, force%tind+1,                         &
                           nyn-nys+1, force%nzw, 2, dynamic_3d )

          IF ( .NOT. neutral )  THEN
             CALL get_variable( id_dynamic, 'ls_forcing_left_pt',              &
                           force%pt_left(0:1,nzb+1:nzt+1,nys:nyn),             &
                           nys+1, nzb+1, force%tind+1,                         &
                           nyn-nys+1, force%nzu, 2, dynamic_3d )
          ENDIF
          IF ( humidity )  THEN
             CALL get_variable( id_dynamic, 'ls_forcing_left_qv',              &
                           force%q_left(0:1,nzb+1:nzt+1,nys:nyn),              &
                           nys+1, nzb+1, force%tind+1,                         &
                           nyn-nys+1, force%nzu, 2, dynamic_3d )
          ENDIF
       ENDIF

       IF ( force_bound_r )  THEN
          CALL get_variable( id_dynamic, 'ls_forcing_right_u',                 &
                           force%u_right(0:1,nzb+1:nzt+1,nys:nyn),             &
                           nys+1, nzb+1, force%tind+1,                         &
                           nyn-nys+1, force%nzu, 2, dynamic_3d )
                           
          CALL get_variable( id_dynamic, 'ls_forcing_right_v',                 &
                           force%v_right(0:1,nzb+1:nzt+1,nysv:nyn),            &
                           nysv, nzb+1, force%tind+1,                          &
                           nyn-nysv+1, force%nzu, 2, dynamic_3d )
                           
          CALL get_variable( id_dynamic, 'ls_forcing_right_w',                 &
                           force%w_right(0:1,nzb+1:nzt,nys:nyn),               &
                           nys+1, nzb+1, force%tind+1,                         &
                           nyn-nys+1, force%nzw, 2, dynamic_3d )
                           
          IF ( .NOT. neutral )  THEN
             CALL get_variable( id_dynamic, 'ls_forcing_right_pt',             &
                           force%pt_right(0:1,nzb+1:nzt+1,nys:nyn),            &
                           nys+1, nzb+1, force%tind+1,                         &
                           nyn-nys+1, force%nzu, 2, dynamic_3d )
          ENDIF
          IF ( humidity )  THEN
             CALL get_variable( id_dynamic, 'ls_forcing_right_qv',             &
                           force%q_right(0:1,nzb+1:nzt+1,nys:nyn),             &
                           nys+1, nzb+1, force%tind+1,                         &
                           nyn-nys+1, force%nzu, 2, dynamic_3d )
          ENDIF
       ENDIF

       IF ( force_bound_n )  THEN
       
          CALL get_variable( id_dynamic, 'ls_forcing_north_u',                 &
                           force%u_north(0:1,nzb+1:nzt+1,nxlu:nxr),            &
                           nxlu, nzb+1, force%tind+1,                          &
                           nxr-nxlu+1, force%nzu, 2, dynamic_3d )

          CALL get_variable( id_dynamic, 'ls_forcing_north_v',                 &
                           force%v_north(0:1,nzb+1:nzt+1,nxl:nxr),             &
                           nxl+1, nzb+1, force%tind+1,                         &
                           nxr-nxl+1, force%nzu, 2, dynamic_3d )
                           
          CALL get_variable( id_dynamic, 'ls_forcing_north_w',                 &
                           force%w_north(0:1,nzb+1:nzt,nxl:nxr),               &
                           nxl+1, nzb+1, force%tind+1,                         &
                           nxr-nxl+1, force%nzw, 2, dynamic_3d )
                           
          IF ( .NOT. neutral )  THEN
             CALL get_variable( id_dynamic, 'ls_forcing_north_pt',             &
                           force%pt_north(0:1,nzb+1:nzt+1,nxl:nxr),            &
                           nxl+1, nzb+1, force%tind+1,                         &
                           nxr-nxl+1, force%nzu, 2, dynamic_3d )
          ENDIF
          IF ( humidity )  THEN
             CALL get_variable( id_dynamic, 'ls_forcing_north_qv',             &
                           force%q_north(0:1,nzb+1:nzt+1,nxl:nxr),             &
                           nxl+1, nzb+1, force%tind+1,                         &
                           nxr-nxl+1, force%nzu, 2, dynamic_3d )
          ENDIF
       ENDIF

       IF ( force_bound_s )  THEN
          CALL get_variable( id_dynamic, 'ls_forcing_south_u',                 &
                           force%u_south(0:1,nzb+1:nzt+1,nxlu:nxr),            &
                           nxlu, nzb+1, force%tind+1,                          &
                           nxr-nxlu+1, force%nzu, 2, dynamic_3d )

          CALL get_variable( id_dynamic, 'ls_forcing_south_v',                 &
                           force%v_south(0:1,nzb+1:nzt+1,nxl:nxr),             &
                           nxl+1, nzb+1, force%tind+1,                         &
                           nxr-nxl+1, force%nzu, 2, dynamic_3d )
                           
          CALL get_variable( id_dynamic, 'ls_forcing_south_w',                 &
                           force%w_south(0:1,nzb+1:nzt,nxl:nxr),               &
                           nxl+1, nzb+1, force%tind+1,                         &
                           nxr-nxl+1, force%nzw, 2, dynamic_3d )
                           
          IF ( .NOT. neutral )  THEN
             CALL get_variable( id_dynamic, 'ls_forcing_south_pt',             &
                           force%pt_south(0:1,nzb+1:nzt+1,nxl:nxr),            &
                           nxl+1, nzb+1, force%tind+1,                         &
                           nxr-nxl+1, force%nzu, 2, dynamic_3d )
          ENDIF
          IF ( humidity )  THEN
             CALL get_variable( id_dynamic, 'ls_forcing_south_qv',             &
                           force%q_south(0:1,nzb+1:nzt+1,nxl:nxr),             &
                           nxl+1, nzb+1, force%tind+1,                         &
                           nxr-nxl+1, force%nzu, 2, dynamic_3d )
          ENDIF
       ENDIF
!
!--    Top boundary
       CALL get_variable( id_dynamic, 'ls_forcing_top_u',                      &
                             force%u_top(0:1,nys:nyn,nxlu:nxr),                &
                             nxlu, nys+1, force%tind+1,                        &
                             nxr-nxlu+1, nyn-nys+1, 2, dynamic_3d )

       CALL get_variable( id_dynamic, 'ls_forcing_top_v',                      &
                             force%v_top(0:1,nysv:nyn,nxl:nxr),                &
                             nxl+1, nysv, force%tind+1,                        &
                             nxr-nxl+1, nyn-nysv+1, 2, dynamic_3d )
                             
       CALL get_variable( id_dynamic, 'ls_forcing_top_w',                      &
                             force%w_top(0:1,nys:nyn,nxl:nxr),                 &
                             nxl+1, nys+1, force%tind+1,                       &
                             nxr-nxl+1, nyn-nys+1, 2, dynamic_3d )
                             
       IF ( .NOT. neutral )  THEN
          CALL get_variable( id_dynamic, 'ls_forcing_top_pt',                  &
                                force%pt_top(0:1,nys:nyn,nxl:nxr),             &
                                nxl+1, nys+1, force%tind+1,                    &
                                nxr-nxl+1, nyn-nys+1, 2, dynamic_3d )
       ENDIF
       IF ( humidity )  THEN
          CALL get_variable( id_dynamic, 'ls_forcing_top_qv',                  &
                                force%q_top(0:1,nys:nyn,nxl:nxr),              &
                                nxl+1, nys+1, force%tind+1,                    &
                                nxr-nxl+1, nyn-nys+1, 2, dynamic_3d )
       ENDIF

!
!--    Close input file
       CALL close_input_file( id_dynamic )
#endif
!
!--    End of CPU measurement
       CALL cpu_log( log_point_s(86), 'NetCDF input forcing', 'stop' )

!
!--    Finally, after data input set control flag indicating that vertical
!--    inter- and/or extrapolation is required.
!--    Please note, inter/extrapolation of INIFOR data is only a workaroud,
!--    as long as INIFOR delivers vertically equidistant data.
       force%interpolated = .FALSE.

    END SUBROUTINE netcdf_data_input_lsf


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks input file for consistency and minimum requirements.
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_check_dynamic

       USE control_parameters,                                                 &
           ONLY:  initializing_actions, forcing, message_string

       IMPLICIT NONE

!
!--    In case of forcing, check whether dynamic input file is present
       IF ( .NOT. input_pids_dynamic  .AND.  forcing )  THEN
          message_string = 'forcing = .TRUE. requires dynamic input file ' //  &
                            TRIM( input_file_dynamic ) // TRIM( coupling_char )
          CALL message( 'netcdf_data_input_mod', 'NDI009', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Dynamic input file must also be present if initialization via inifor is
!--    prescribed.
       IF ( .NOT. input_pids_dynamic  .AND.                                    &
            TRIM( initializing_actions ) == 'inifor' )  THEN
          message_string = 'initializing_actions = inifor requires dynamic ' //&
                           'input file ' // TRIM( input_file_dynamic ) //      &
                           TRIM( coupling_char )
          CALL message( 'netcdf_data_input_mod', 'NDI010', 1, 2, 0, 6, 0 )
       ENDIF

    END SUBROUTINE netcdf_data_input_check_dynamic

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks input file for consistency and minimum requirements.
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_check_static

       USE arrays_3d,                                                          &
           ONLY:  zu

       USE control_parameters,                                                 &
           ONLY:  land_surface, message_string, urban_surface

       USE grid_variables,                                                     &
           ONLY:  dx, dy

       USE indices,                                                            &
           ONLY:  nx, nxl, nxr, ny, nyn, nys

       IMPLICIT NONE

       INTEGER(iwp) ::  i      !< loop index along x-direction
       INTEGER(iwp) ::  j      !< loop index along y-direction
       INTEGER(iwp) ::  n_surf !< number of different surface types at given location

       LOGICAL      ::  check_passed !< flag indicating if a check passed

!
!--    Return if no static input file is available
       IF ( .NOT. input_pids_static )  RETURN
!
!--    Check whether dimension size in input file matches the model dimensions
       IF ( dim_static%nx-1 /= nx  .OR.  dim_static%ny-1 /= ny )  THEN
          message_string = 'Static input file: horizontal dimension in ' //    &
                           'x- and/or y-direction ' //                         &
                           'do not match the respective model dimension'
          CALL message( 'netcdf_data_input_mod', 'NDI011', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Check if grid spacing of provided input data matches the respective
!--    grid spacing in the model.
       IF ( dim_static%x(1) - dim_static%x(0) /= dx  .OR.                      &
            dim_static%y(1) - dim_static%y(0) /= dy )  THEN
          message_string = 'Static input file: horizontal grid spacing ' //    &
                           'in x- and/or y-direction ' //                      &
                           'do not match the respective model grid spacing.'
          CALL message( 'netcdf_data_input_mod', 'NDI012', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Check orography for fill-values. For the moment, give an error message.
!--    More advanced methods, e.g. a nearest neighbor algorithm as used in GIS
!--    systems might be implemented later.
!--    Please note, if no terrain height is provided, it is set to 0.
       IF ( ANY( terrain_height_f%var == terrain_height_f%fill ) )  THEN
          message_string = 'NetCDF variable zt is not ' //                     &
                           'allowed to have missing data'
          CALL message( 'netcdf_data_input_mod', 'NDI013', 2, 2, myid, 6, 0 )
       ENDIF
!
!--    If 3D buildings are read, check if building information is consistent
!--    to numeric grid.
       IF ( buildings_f%from_file )  THEN
          IF ( buildings_f%lod == 2 )  THEN
             IF ( buildings_f%nz > SIZE( zu ) )  THEN
                message_string = 'Reading 3D building data - too much ' //     &
                                 'data points along the vertical coordinate.'
                CALL message( 'netcdf_data_input_mod', 'NDI014', 2, 2, 0, 6, 0 )
             ENDIF

             IF ( ANY( buildings_f%z(0:buildings_f%nz-1) /=                    &
                       zu(0:buildings_f%nz-1) ) )  THEN
                message_string = 'Reading 3D building data - vertical ' //     &
                                 'coordinate do not match numeric grid.'
                CALL message( 'netcdf_data_input_mod', 'NDI015', 2, 2, 0, 6, 0 )
             ENDIF
          ENDIF
       ENDIF

!
!--    Skip further checks concerning buildings and natural surface properties
!--    if no urban surface and land surface model are applied.
       IF (  .NOT. land_surface  .OR.  .NOT. urban_surface )  RETURN
!
!--    Check for minimum requirement of surface-classification data in case
!--    static input file is used.
       IF ( ( .NOT. vegetation_type_f%from_file  .OR.                          &
              .NOT. pavement_type_f%from_file    .OR.                          &
              .NOT. water_type_f%from_file       .OR.                          &
              .NOT. soil_type_f%from_file             ) .OR.                   &
             ( urban_surface  .AND.  .NOT. building_type_f%from_file ) )  THEN
          message_string = 'Minimum requirement for surface classification ' //&
                           'is not fulfilled. At least ' //                    &
                           'vegetation_type, pavement_type, ' //               &
                           'soil_type and water_type are '//                   &
                           'required. If urban-surface model is applied, ' //  &
                           'also building_type ist required'
          CALL message( 'netcdf_data_input_mod', 'NDI016', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Check for general availability of input variables.
!--    If vegetation_type is 0 at any location, vegetation_pars as well as
!--    root_area_dens_s are required.
       IF ( vegetation_type_f%from_file )  THEN
          IF ( ANY( vegetation_type_f%var == 0 ) )  THEN
             IF ( .NOT. vegetation_pars_f%from_file )  THEN
                message_string = 'If vegegation_type = 0 at any location, ' // &
                                 'vegetation_pars is required'
                CALL message( 'netcdf_data_input_mod', 'NDI017', 2, 2, -1, 6, 0 )
             ENDIF
             IF ( .NOT. root_area_density_lsm_f%from_file )  THEN
                message_string = 'If vegegation_type = 0 at any location, ' // &
                                 'root_area_dens_s is required'
                CALL message( 'netcdf_data_input_mod', 'NDI018', 2, 2, myid, 6, 0 )
             ENDIF
          ENDIF
       ENDIF
!
!--    If soil_type is zero at any location, soil_pars is required.
       IF ( soil_type_f%from_file )  THEN
          check_passed = .TRUE.
          IF ( ALLOCATED( soil_type_f%var_2d ) )  THEN
             IF ( ANY( soil_type_f%var_2d == 0 ) )  THEN
                IF ( .NOT. soil_pars_f%from_file )  check_passed = .FALSE.
             ENDIF
          ELSE
             IF ( ANY( soil_type_f%var_3d == 0 ) )  THEN
                IF ( .NOT. soil_pars_f%from_file )  check_passed = .FALSE.
             ENDIF
          ENDIF
          IF ( .NOT. check_passed )  THEN
             message_string = 'If soil_type = 0 at any location, ' //          &
                              'soil_pars is required'
             CALL message( 'netcdf_data_input_mod', 'NDI019', 2, 2, myid, 6, 0 )
          ENDIF
       ENDIF
!
!--    If building_type is zero at any location, building_pars is required.
       IF ( building_type_f%from_file )  THEN
          IF ( ANY( building_type_f%var == 0 ) )  THEN
             IF ( .NOT. building_pars_f%from_file )  THEN
                message_string = 'If building_type = 0 at any location, ' //   &
                                 'building_pars is required'
                CALL message( 'netcdf_data_input_mod', 'NDI020', 2, 2, myid, 6, 0 )
             ENDIF
          ENDIF
       ENDIF
!
!--    If albedo_type is zero at any location, albedo_pars is required.
       IF ( albedo_type_f%from_file )  THEN
          IF ( ANY( albedo_type_f%var == 0 ) )  THEN
             IF ( .NOT. albedo_pars_f%from_file )  THEN
                message_string = 'If albedo_type = 0 at any location, ' //     &
                                 'albedo_pars is required'
                CALL message( 'netcdf_data_input_mod', 'NDI021', 2, 2, myid, 6, 0 )
             ENDIF
          ENDIF
       ENDIF
!
!--    If pavement_type is zero at any location, pavement_pars is required.
       IF ( pavement_type_f%from_file )  THEN
          IF ( ANY( pavement_type_f%var == 0 ) )  THEN
             IF ( .NOT. pavement_pars_f%from_file )  THEN
                message_string = 'If pavement_type = 0 at any location, ' //   &
                                 'pavement_pars is required'
                CALL message( 'netcdf_data_input_mod', 'NDI022', 2, 2, myid, 6, 0 )
             ENDIF
          ENDIF
       ENDIF
!
!--    If pavement_type is zero at any location, also pavement_subsurface_pars
!--    is required.
       IF ( pavement_type_f%from_file )  THEN
          IF ( ANY( pavement_type_f%var == 0 ) )  THEN
             IF ( .NOT. pavement_subsurface_pars_f%from_file )  THEN
                message_string = 'If pavement_type = 0 at any location, ' //   &
                                 'pavement_subsurface_pars is required'
                CALL message( 'netcdf_data_input_mod', 'NDI023', 2, 2, myid, 6, 0 )
             ENDIF
          ENDIF
       ENDIF
!
!--    If water_type is zero at any location, water_pars is required.
       IF ( water_type_f%from_file )  THEN
          IF ( ANY( water_type_f%var == 0 ) )  THEN
             IF ( .NOT. water_pars_f%from_file )  THEN
                message_string = 'If water_type = 0 at any location, ' //      &
                                 'water_pars is required'
                CALL message( 'netcdf_data_input_mod', 'NDI024', 2, 2,myid, 6, 0 )
             ENDIF
          ENDIF
       ENDIF
!
!--    Check for local consistency of the input data.
       DO  i = nxl, nxr
          DO  j = nys, nyn
!
!--          For each (y,x)-location at least one of the parameters
!--          vegetation_type, pavement_type, building_type, or water_type
!--          must be set to a non­missing value.
             IF ( vegetation_type_f%var(j,i) == vegetation_type_f%fill  .AND.  &
                  pavement_type_f%var(j,i)   == pavement_type_f%fill    .AND.  &
                  building_type_f%var(j,i)   == building_type_f%fill    .AND.  &
                  water_type_f%var(j,i)      == water_type_f%fill )  THEN
                WRITE( message_string, * ) 'At least one of the parameters '// &
                                 'vegetation_type, pavement_type, '     //     &
                                 'building_type, or water_type must be set '// &
                                 'to a non-missing value. Grid point: ', j, i
                CALL message( 'netcdf_data_input_mod', 'NDI025', 2, 2, myid, 6, 0 )
             ENDIF
!
!--          Note that a soil_type is required for each location (y,x) where
!--          either vegetation_type or pavement_type is a non­missing value.
             IF ( ( vegetation_type_f%var(j,i) /= vegetation_type_f%fill  .OR. &
                    pavement_type_f%var(j,i)   /= pavement_type_f%fill ) )  THEN
                check_passed = .TRUE.
                IF ( ALLOCATED( soil_type_f%var_2d ) )  THEN
                   IF ( soil_type_f%var_2d(j,i) == soil_type_f%fill )          &
                      check_passed = .FALSE.
                ELSE
                   IF ( ANY( soil_type_f%var_3d(:,j,i) == soil_type_f%fill) )  &
                      check_passed = .FALSE.
                ENDIF

                IF ( .NOT. check_passed )  THEN
                   message_string = 'soil_type is required for each '//        &
                                 'location (y,x) where vegetation_type or ' // &
                                 'pavement_type is a non-missing value.'
                   CALL message( 'netcdf_data_input_mod', 'NDI026',            &
                                  2, 2, myid, 6, 0 )
                ENDIF
             ENDIF
!
!--          Check for consistency of surface fraction. If more than one type
!--          is set, surface fraction need to be given and the sum must not
!--          be larger than 1.
             n_surf = 0
             IF ( vegetation_type_f%var(j,i) /= vegetation_type_f%fill )       &
                n_surf = n_surf + 1
             IF ( water_type_f%var(j,i)      /= water_type_f%fill )            &
                n_surf = n_surf + 1
             IF ( pavement_type_f%var(j,i)   /= pavement_type_f%fill )         &
                n_surf = n_surf + 1

             IF ( n_surf > 1 )  THEN
                IF ( .NOT. surface_fraction_f%from_file )  THEN
                   message_string = 'If more than one surface type is ' //     &
                                 'given at a location, surface_fraction ' //   &
                                 'must be provided.'
                   CALL message( 'netcdf_data_input_mod', 'NDI027',            &
                                  2, 2, myid, 6, 0 )
                ELSEIF ( ANY ( surface_fraction_f%frac(:,j,i) ==               &
                               surface_fraction_f%fill ) )  THEN
                   message_string = 'If more than one surface type is ' //     &
                                 'given at a location, surface_fraction ' //   &
                                 'must be provided.'
                   CALL message( 'netcdf_data_input_mod', 'NDI027',            &
                                  2, 2, myid, 6, 0 )
                ENDIF
             ENDIF
!
!--          Check for further mismatches. e.g. relative fractions exceed 1 or
!--          vegetation_type is set but surface vegetation fraction is zero, 
!--          etc..
             IF ( surface_fraction_f%from_file )  THEN
!
!--             Sum of relative fractions must not exceed 1.
                IF ( SUM ( surface_fraction_f%frac(:,j,i) ) > 1.0_wp )  THEN
                   message_string = 'surface_fraction must not exceed 1'
                   CALL message( 'netcdf_data_input_mod', 'NDI028',            &
                                  2, 2, myid, 6, 0 )
                ENDIF
!
!--             Relative fraction for a type must not be zero at locations where 
!--             this type is set. 
                IF (                                                           &
                  ( vegetation_type_f%var(j,i) /= vegetation_type_f%fill  .AND.&
                 ( surface_fraction_f%frac(ind_veg_wall,j,i) == 0.0_wp .OR.    &
                   surface_fraction_f%frac(ind_veg_wall,j,i) ==                &
                                                     surface_fraction_f%fill ) &
                  )  .OR.                                                      &
                  ( pavement_type_f%var(j,i) /= pavement_type_f%fill     .AND. &
                 ( surface_fraction_f%frac(ind_pav_green,j,i) == 0.0_wp .OR.   &
                   surface_fraction_f%frac(ind_pav_green,j,i) ==               &
                                                     surface_fraction_f%fill ) &
                  )  .OR.                                                      &
                  ( water_type_f%var(j,i) /= water_type_f%fill           .AND. &
                 ( surface_fraction_f%frac(ind_wat_win,j,i) == 0.0_wp .OR.     &
                   surface_fraction_f%frac(ind_wat_win,j,i) ==                 &
                                                     surface_fraction_f%fill ) &
                  ) )  THEN
                   WRITE( message_string, * ) 'Mismatch in setting of '     // &
                             'surface_fraction. Vegetation-, pavement-, or '// &
                             'water surface is given at (i,j) = ( ', i, j,     &
                             ' ), but surface fraction is 0 for the given type.'
                   CALL message( 'netcdf_data_input_mod', 'NDI029',            &
                                  2, 2, myid, 6, 0 )
                ENDIF
!
!--             Relative fraction for a type must not contain non-zero values
!--             if this type is not set. 
                IF (                                                           &
                  ( vegetation_type_f%var(j,i) == vegetation_type_f%fill  .AND.&
                 ( surface_fraction_f%frac(ind_veg_wall,j,i) /= 0.0_wp .AND.   &
                   surface_fraction_f%frac(ind_veg_wall,j,i) /=                &
                                                     surface_fraction_f%fill ) &
                  )  .OR.                                                      &
                  ( pavement_type_f%var(j,i) == pavement_type_f%fill     .AND. &
                 ( surface_fraction_f%frac(ind_pav_green,j,i) /= 0.0_wp .AND.  &
                   surface_fraction_f%frac(ind_pav_green,j,i) /=               &
                                                     surface_fraction_f%fill ) &
                  )  .OR.                                                      &
                  ( water_type_f%var(j,i) == water_type_f%fill           .AND. &
                 ( surface_fraction_f%frac(ind_wat_win,j,i) /= 0.0_wp .AND.    &
                   surface_fraction_f%frac(ind_wat_win,j,i) /=                 &
                                                     surface_fraction_f%fill ) &
                  ) )  THEN
                   WRITE( message_string, * ) 'Mismatch in setting of '     // &
                             'surface_fraction. Vegetation-, pavement-, or '// &
                             'water surface is not given at (i,j) = ( ', i, j, &
                             ' ), but surface fraction is not 0 for the ' //   &
                             'given type.'
                   CALL message( 'netcdf_data_input_mod', 'NDI030',            &
                                  2, 2, myid, 6, 0 )
                ENDIF
             ENDIF
!
!--          Check vegetation_pars. If vegetation_type is 0, all parameters
!--          need to be set, otherwise, single parameters set by
!--          vegetation_type can be overwritten.
             IF ( vegetation_type_f%from_file )  THEN
                IF ( vegetation_type_f%var(j,i) == 0 )  THEN
                   IF ( ANY( vegetation_pars_f%pars_xy(:,j,i) ==               &
                             vegetation_pars_f%fill ) )  THEN
                      message_string = 'If vegetation_type(y,x) = 0, all '  // &
                                       'parameters of vegetation_pars at '//   &
                                       'this location must be set.'
                      CALL message( 'netcdf_data_input_mod', 'NDI031',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ENDIF
             ENDIF
!
!--          Check root distribution. If vegetation_type is 0, all levels must
!--          be set.
             IF ( vegetation_type_f%from_file )  THEN
                IF ( vegetation_type_f%var(j,i) == 0 )  THEN
                   IF ( ANY( root_area_density_lsm_f%var(:,j,i) ==             &
                             root_area_density_lsm_f%fill ) )  THEN
                      message_string = 'If vegetation_type(y,x) = 0, all ' //  &
                                       'levels of root_area_dens_s ' //        &
                                       'must be set at this location.'
                      CALL message( 'netcdf_data_input_mod', 'NDI032',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ENDIF
             ENDIF
!
!--          Check soil parameters. If soil_type is 0, all parameters
!--          must be set.
             IF ( soil_type_f%from_file )  THEN
                check_passed = .TRUE.
                IF ( ALLOCATED( soil_type_f%var_2d ) )  THEN
                   IF ( soil_type_f%var_2d(j,i) == 0 )  THEN
                      IF ( ANY( soil_pars_f%pars_xy(:,j,i) ==                  &
                                soil_pars_f%fill ) )  check_passed = .FALSE.
                   ENDIF
                ELSE
                   IF ( ANY( soil_type_f%var_3d(:,j,i) == 0 ) )  THEN
                      IF ( ANY( soil_pars_f%pars_xy(:,j,i) ==                  &
                                soil_pars_f%fill ) )  check_passed = .FALSE.
                   ENDIF
                ENDIF
                IF ( .NOT. check_passed )  THEN
                   message_string = 'If soil_type(y,x) = 0, all levels of '  //&
                                    'soil_pars at this location must be set.'
                   CALL message( 'netcdf_data_input_mod', 'NDI033',            &
                                  2, 2, myid, 6, 0 )
                ENDIF
             ENDIF

!
!--          Check building parameters. If building_type is 0, all parameters
!--          must be set.
             IF ( building_type_f%from_file )  THEN
                IF ( building_type_f%var(j,i) == 0 )  THEN
                   IF ( ANY( building_pars_f%pars_xy(:,j,i) ==                 &
                             building_pars_f%fill ) )  THEN
                      message_string = 'If building_type(y,x) = 0, all ' //    &
                                       'parameters of building_pars at this '//&
                                       'location must be set.'
                      CALL message( 'netcdf_data_input_mod', 'NDI034',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ENDIF
             ENDIF
!
!--          Check if building_type is set at each building and vice versa.
             IF ( building_type_f%from_file  .AND.  buildings_f%from_file )  THEN
                IF ( buildings_f%lod == 1 )  THEN
                   IF ( buildings_f%var_2d(j,i)  /= buildings_f%fill1  .AND.   &
                        building_type_f%var(j,i) == building_type_f%fill )  THEN

                      WRITE( message_string, * ) 'Each location where a ' //   &
                                         'building is set requires a type ' // &
                                         '( and vice versa ) in case the ' //  &
                                         'urban-surface model is applied. ' // &
                                         'i, j = ', i, j
                      CALL message( 'netcdf_data_input_mod', 'NDI035',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ENDIF
                IF ( buildings_f%lod == 2 )  THEN
                   IF ( ANY( buildings_f%var_3d(:,j,i) == 1 )  .AND.           &
                        building_type_f%var(j,i) == building_type_f%fill )  THEN
                      WRITE( message_string, * ) 'Each location where a ' //   &
                                         'building is set requires a type ' // &
                                         '( and vice versa ) in case the ' //  &
                                         'urban-surface model is applied. ' // &
                                         'i, j = ', i, j
                      CALL message( 'netcdf_data_input_mod', 'NDI035',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ENDIF
             ENDIF
!
!--          Check if at each location where a building is present also an ID
!--          is set and vice versa.
             IF ( buildings_f%from_file )  THEN
                IF ( buildings_f%lod == 1 )  THEN
                   IF ( buildings_f%var_2d(j,i) /= buildings_f%fill1  .AND.    &
                        building_id_f%var(j,i)  == building_id_f%fill )  THEN
                      WRITE( message_string, * ) 'Each location where a ' //   &
                                         'building is set requires an ID ' //  &
                                         '( and vice versa ). i, j = ', i, j
                      CALL message( 'netcdf_data_input_mod', 'NDI036',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ELSEIF ( buildings_f%lod == 2 )  THEN
                   IF ( ANY( buildings_f%var_3d(:,j,i) == 1 )  .AND.           &
                        building_id_f%var(j,i) == building_id_f%fill )  THEN
                      WRITE( message_string, * ) 'Each location where a ' //   &
                                         'building is set requires an ID ' //  &
                                         '( and vice versa ). i, j = ', i, j
                      CALL message( 'netcdf_data_input_mod', 'NDI036',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ENDIF
             ENDIF
!
!--          Check if at each location where a building ID or a -type is set
!--          also a bulding is defined.
             IF ( buildings_f%from_file )  THEN
                IF ( buildings_f%lod == 1 )  THEN
                   IF ( buildings_f%var_2d(j,i)  /= buildings_f%fill1  .AND.   &
                        building_id_f%var(j,i) == building_id_f%fill )  THEN
                      WRITE( message_string, * ) 'Each building grid point '// &
                                                 'requires an ID.', i, j
                      CALL message( 'netcdf_data_input_mod', 'NDI036',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ELSEIF ( buildings_f%lod == 2 )  THEN
                   IF ( ANY( buildings_f%var_3d(:,j,i) == 1 )                  &
                  .AND. building_id_f%var(j,i) == building_id_f%fill )  THEN
                      WRITE( message_string, * ) 'Each building grid point '// &
                                                 'requires an ID.', i, j
                      CALL message( 'netcdf_data_input_mod', 'NDI036',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ENDIF
             ENDIF
!
!--          Check albedo parameters. If albedo_type is 0, all parameters
!--          must be set.
             IF ( albedo_type_f%from_file )  THEN
                IF ( albedo_type_f%var(j,i) == 0 )  THEN
                   IF ( ANY( albedo_pars_f%pars_xy(:,j,i) ==                   &
                             albedo_pars_f%fill ) )  THEN
                      message_string = 'If albedo_type(y,x) = 0, all ' //      &
                                       'parameters of albedo_pars at this ' // &
                                       'location must be set.'
                      CALL message( 'netcdf_data_input_mod', 'NDI037',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ENDIF
             ENDIF

!
!--          Check pavement parameters. If pavement_type is 0, all parameters
!--          of pavement_pars must be set at this location.
             IF ( pavement_type_f%from_file )  THEN
                IF ( pavement_type_f%var(j,i) == 0 )  THEN
                   IF ( ANY( pavement_pars_f%pars_xy(:,j,i) ==                 &
                             pavement_pars_f%fill ) )  THEN
                      message_string = 'If pavement_type(y,x) = 0, all ' //    &
                                       'parameters of pavement_pars at this '//&
                                       'location must be set.'
                      CALL message( 'netcdf_data_input_mod', 'NDI038',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ENDIF
             ENDIF
!
!--          Check pavement-subsurface parameters. If pavement_type is 0,
!--          all parameters of pavement_subsurface_pars must be set  at this
!--          location.
             IF ( pavement_type_f%from_file )  THEN
                IF ( pavement_type_f%var(j,i) == 0 )  THEN
                   IF ( ANY( pavement_subsurface_pars_f%pars_xyz(:,:,j,i) ==   &
                             pavement_subsurface_pars_f%fill ) )  THEN
                      message_string = 'If pavement_type(y,x) = 0, all ' //    &
                                       'parameters of '                  //    &
                                       'pavement_subsurface_pars at this '//   &
                                       'location must be set.'
                      CALL message( 'netcdf_data_input_mod', 'NDI039',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ENDIF
             ENDIF

!
!--          Check water parameters. If water_type is 0, all parameters
!--          must be set  at this location.
             IF ( water_type_f%from_file )  THEN
                IF ( water_type_f%var(j,i) == 0 )  THEN
                   IF ( ANY( water_pars_f%pars_xy(:,j,i) ==                    &
                             water_pars_f%fill ) )  THEN
                      message_string = 'If water_type(y,x) = 0, all ' //       &
                                       'parameters of water_pars at this ' //  &
                                       'location must be set.'
                      CALL message( 'netcdf_data_input_mod', 'NDI040',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ENDIF
             ENDIF

          ENDDO
       ENDDO

    END SUBROUTINE netcdf_data_input_check_static

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Vertical interpolation and extrapolation of 1D variables.
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_interpolate_1d( var, z_grid, z_file)

       IMPLICIT NONE

       LOGICAL      ::  top     !< flag indicating extrapolation at model top

       INTEGER(iwp) ::  k       !< running index z-direction file
       INTEGER(iwp) ::  kk      !< running index z-direction stretched model grid
       INTEGER(iwp) ::  kl      !< lower index bound along z-direction
       INTEGER(iwp) ::  ku      !< upper index bound along z-direction
       INTEGER(iwp) ::  nz_file !< number of vertical levels on file


       REAL(wp), DIMENSION(:) ::  z_grid                  !< grid levels on numeric grid
       REAL(wp), DIMENSION(:) ::  z_file                  !< grid levels on file grid
       REAL(wp), DIMENSION(:), INTENT(INOUT) ::  var      !< treated variable
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  var_tmp  !< temporary variable


       kl = LBOUND(var,1)
       ku = UBOUND(var,1)
       ALLOCATE( var_tmp(kl:ku) )

       DO  k = kl, ku

          kk = MINLOC( ABS( z_file - z_grid(k) ), DIM = 1 )

          IF ( kk < ku )  THEN
             IF ( z_file(kk) - z_grid(k) <= 0.0_wp )  THEN
                var_tmp(k) = var(kk) +                                         &
                                       ( var(kk+1)        - var(kk)    ) /     &
                                       ( z_file(kk+1)     - z_file(kk) ) *     &
                                       ( z_grid(k)        - z_file(kk) )

             ELSEIF ( z_file(kk) - z_grid(k) > 0.0_wp )  THEN
                var_tmp(k) = var(kk-1) +                                       &
                                         ( var(kk)     - var(kk-1)    ) /      &
                                         ( z_file(kk)  - z_file(kk-1) ) *      &
                                         ( z_grid(k)   - z_file(kk-1) )
             ENDIF
!
!--       Extrapolate
          ELSE

             var_tmp(k) = var(ku) +   ( var(ku)    - var(ku-1)      ) /        &
                                      ( z_file(ku) - z_file(ku-1)   ) *        &
                                      ( z_grid(k)  - z_file(ku)     )

          ENDIF

       ENDDO
       var(:) = var_tmp(:)

       DEALLOCATE( var_tmp )


    END SUBROUTINE netcdf_data_input_interpolate_1d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Vertical interpolation and extrapolation of 1D variables from Inifor grid
!> onto Palm grid, where both have same dimension. Please note, the passed
!> paramter list in 1D version is different compared to 2D version.
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_interpolate_1d_soil( var, var_file,           &
                                                      z_grid, z_file,          &
                                                      nzb_var, nzt_var,        &
                                                      nzb_file, nzt_file )

       IMPLICIT NONE

       INTEGER(iwp) ::  i        !< running index x-direction
       INTEGER(iwp) ::  j        !< running index y-direction
       INTEGER(iwp) ::  k        !< running index z-direction file
       INTEGER(iwp) ::  kk       !< running index z-direction stretched model grid
       INTEGER(iwp) ::  ku       !< upper index bound along z-direction for varialbe from file
       INTEGER(iwp) ::  nzb_var  !< lower bound of final array
       INTEGER(iwp) ::  nzt_var  !< upper bound of final array
       INTEGER(iwp) ::  nzb_file !< lower bound of file array
       INTEGER(iwp) ::  nzt_file !< upper bound of file array

!        LOGICAL, OPTIONAL ::  zsoil !< flag indicating reverse z-axis, i.e. zsoil instead of height, e.g. in case of ocean or soil

       REAL(wp), DIMENSION(nzb_var:nzt_var)   ::  z_grid   !< grid levels on numeric grid
       REAL(wp), DIMENSION(nzb_file:nzt_file) ::  z_file   !< grid levels on file grid
       REAL(wp), DIMENSION(nzb_var:nzt_var)   ::  var      !< treated variable
       REAL(wp), DIMENSION(nzb_file:nzt_file) ::  var_file !< temporary variable

       ku = nzt_file

       DO  k = nzb_var, nzt_var
!
!--       Determine index on Inifor grid which is closest to the actual height
          kk = MINLOC( ABS( z_file - z_grid(k) ), DIM = 1 )
!
!--       If closest index on Inifor grid is smaller than top index,
!--       interpolate the data
          IF ( kk < nzt_file )  THEN
             IF ( z_file(kk) - z_grid(k) <= 0.0_wp )  THEN
                var(k) = var_file(kk) + ( var_file(kk+1) - var_file(kk) ) /    &
                                        ( z_file(kk+1)   - z_file(kk)   ) *    &
                                        ( z_grid(k)      - z_file(kk)   )

             ELSEIF ( z_file(kk) - z_grid(k) > 0.0_wp )  THEN
                var(k) = var_file(kk-1) + ( var_file(kk) - var_file(kk-1) ) /  &
                                          ( z_file(kk)   - z_file(kk-1)   ) *  &
                                          ( z_grid(k)    - z_file(kk-1)   )
             ENDIF
!
!--       Extrapolate if actual height is above the highest Inifor level
          ELSE
             var(k) = var_file(ku) + ( var_file(ku) - var_file(ku-1) ) /       &
                                     ( z_file(ku)   - z_file(ku-1)   ) *       &
                                     ( z_grid(k)    - z_file(ku)     )

          ENDIF

       ENDDO

    END SUBROUTINE netcdf_data_input_interpolate_1d_soil

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Vertical interpolation and extrapolation of 2D variables at lateral boundaries.
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_interpolate_2d( var, z_grid, z_file)

       IMPLICIT NONE

       LOGICAL      ::  top     !< flag indicating extrapolation at model top

       INTEGER(iwp) ::  i       !< running index x- or y -direction
       INTEGER(iwp) ::  il      !< lower index bound along x- or y-direction
       INTEGER(iwp) ::  iu      !< upper index bound along x- or y-direction
       INTEGER(iwp) ::  k       !< running index z-direction file
       INTEGER(iwp) ::  kk      !< running index z-direction stretched model grid
       INTEGER(iwp) ::  kl      !< lower index bound along z-direction
       INTEGER(iwp) ::  ku      !< upper index bound along z-direction
       INTEGER(iwp) ::  nz_file !< number of vertical levels on file


       REAL(wp), DIMENSION(:) ::  z_grid                  !< grid levels on numeric grid
       REAL(wp), DIMENSION(:) ::  z_file                  !< grid levels on file grid
       REAL(wp), DIMENSION(:,:), INTENT(INOUT) ::  var    !< treated variable
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  var_tmp  !< temporary variable


       il = LBOUND(var,2)
       iu = UBOUND(var,2)
       kl = LBOUND(var,1)
       ku = UBOUND(var,1)
       ALLOCATE( var_tmp(kl:ku) )

       DO  i = il, iu
          DO  k = kl, ku

             kk = MINLOC( ABS( z_file - z_grid(k) ), DIM = 1 )

             IF ( kk < ku )  THEN
                IF ( z_file(kk) - z_grid(k) <= 0.0_wp )  THEN
                   var_tmp(k) = var(kk,i) +                                    &
                                          ( var(kk+1,i)      - var(kk,i)  ) /  &
                                          ( z_file(kk+1)     - z_file(kk) ) *  &
                                          ( z_grid(k)        - z_file(kk) )

                ELSEIF ( z_file(kk) - z_grid(k) > 0.0_wp )  THEN
                   var_tmp(k) = var(kk-1,i) +                                  &
                                            ( var(kk,i)   - var(kk-1,i)  ) /   &
                                            ( z_file(kk)  - z_file(kk-1) ) *   &
                                            ( z_grid(k)   - z_file(kk-1) )
                ENDIF
!
!--          Extrapolate
             ELSE

                var_tmp(k) = var(ku,i) + ( var(ku,i)  - var(ku-1,i)    ) /     &
                                         ( z_file(ku) - z_file(ku-1)   ) *     &
                                         ( z_grid(k)  - z_file(ku)     )

             ENDIF

          ENDDO
          var(:,i) = var_tmp(:)

       ENDDO

       DEALLOCATE( var_tmp )


    END SUBROUTINE netcdf_data_input_interpolate_2d

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Vertical interpolation and extrapolation of 3D variables.
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_interpolate_3d( var, z_grid, z_file )

       IMPLICIT NONE

       INTEGER(iwp) ::  i       !< running index x-direction
       INTEGER(iwp) ::  il      !< lower index bound along x-direction
       INTEGER(iwp) ::  iu      !< upper index bound along x-direction
       INTEGER(iwp) ::  j       !< running index y-direction
       INTEGER(iwp) ::  jl      !< lower index bound along x-direction
       INTEGER(iwp) ::  ju      !< upper index bound along x-direction
       INTEGER(iwp) ::  k       !< running index z-direction file
       INTEGER(iwp) ::  kk      !< running index z-direction stretched model grid
       INTEGER(iwp) ::  kl      !< lower index bound along z-direction
       INTEGER(iwp) ::  ku      !< upper index bound along z-direction
       INTEGER(iwp) ::  nz_file !< number of vertical levels on file

       REAL(wp), DIMENSION(:) ::  z_grid                      !< grid levels on numeric grid
       REAL(wp), DIMENSION(:) ::  z_file                      !< grid levels on file grid
       REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::  var      !< treated variable
       REAL(wp), DIMENSION(:), ALLOCATABLE       ::  var_tmp  !< temporary variable

       il = LBOUND(var,3)
       iu = UBOUND(var,3)
       jl = LBOUND(var,2)
       ju = UBOUND(var,2)
       kl = LBOUND(var,1)
       ku = UBOUND(var,1)

       ALLOCATE( var_tmp(kl:ku) )

       DO  i = il, iu
          DO  j = jl, ju
             DO  k = kl, ku

                kk = MINLOC( ABS( z_file - z_grid(k) ), DIM = 1 )

                IF ( kk < ku )  THEN
                   IF ( z_file(kk) - z_grid(k) <= 0.0_wp )  THEN
                      var_tmp(k) = var(kk,j,i) +                               &
                                             ( var(kk+1,j,i) - var(kk,j,i) ) / &
                                             ( z_file(kk+1)  - z_file(kk)  ) * &
                                             ( z_grid(k)     - z_file(kk)  )

                   ELSEIF ( z_file(kk) - z_grid(k) > 0.0_wp )  THEN
                      var_tmp(k) = var(kk-1,j,i) +                             &
                                             ( var(kk,j,i) - var(kk-1,j,i) ) / &
                                             ( z_file(kk)  - z_file(kk-1)  ) * &
                                             ( z_grid(k)   - z_file(kk-1)  )
                   ENDIF
!
!--             Extrapolate
                ELSE
                   var_tmp(k) = var(ku,j,i) +                                  &
                                       ( var(ku,j,i)  - var(ku-1,j,i)   ) /    &
                                       ( z_file(ku)   - z_file(ku-1)    ) *    &
                                       ( z_grid(k)    - z_file(ku)      )

                ENDIF
             ENDDO
             var(:,j,i) = var_tmp(:)
          ENDDO
       ENDDO

       DEALLOCATE( var_tmp )


    END SUBROUTINE netcdf_data_input_interpolate_3d

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks if a given variables is on file
!------------------------------------------------------------------------------!
    FUNCTION check_existence( vars_in_file, var_name )

       IMPLICIT NONE

       CHARACTER(LEN=*) ::  var_name                   !< variable to be checked
       CHARACTER(LEN=*), DIMENSION(:) ::  vars_in_file !< list of variables in file

       INTEGER(iwp) ::  i                              !< loop variable

       LOGICAL ::  check_existence                     !< flag indicating whether a variable exist or not - actual return value

       i = 1
       check_existence = .FALSE.
       DO  WHILE ( i <= SIZE( vars_in_file ) )
          check_existence = TRIM( vars_in_file(i) ) == TRIM( var_name )  .OR.  &
                            check_existence
          i = i + 1
       ENDDO

       RETURN

    END FUNCTION check_existence


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Closes an existing netCDF file.
!------------------------------------------------------------------------------!
    SUBROUTINE close_input_file( id )
#if defined( __netcdf )

       USE pegrid

       IMPLICIT NONE

       INTEGER(iwp), INTENT(INOUT)        ::  id        !< file id

       nc_stat = NF90_CLOSE( id )
       CALL handle_error( 'close', 540 )
#endif
    END SUBROUTINE close_input_file

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Opens an existing netCDF file for reading only and returns its id.
!------------------------------------------------------------------------------!
    SUBROUTINE open_read_file( filename, id )
#if defined( __netcdf )

       USE pegrid

       IMPLICIT NONE

       CHARACTER (LEN=*), INTENT(IN) ::  filename  !< filename
       INTEGER(iwp), INTENT(INOUT)   ::  id        !< file id
       LOGICAL                       ::  file_open = .FALSE.

#if defined( __netcdf4_parallel )
!      if __netcdf4_parallel is defined, parrallel NetCDF will be used unconditionally
       nc_stat = NF90_OPEN( filename, IOR( NF90_WRITE, NF90_MPIIO ), id,  &
                            COMM = comm2d, INFO = MPI_INFO_NULL )
       IF(nc_stat /= NF90_NOERR )  THEN                                       !possible NetCDF 3 file
           nc_stat = NF90_OPEN( filename, NF90_NOWRITE, id )
           collective_read = .FALSE.
       ELSE
           collective_read = .TRUE.
       END IF
#else
!      All MPI processes open und read
       nc_stat = NF90_OPEN( filename, NF90_NOWRITE, id )
#endif

       CALL handle_error( 'open_read_file', 539 )

#endif
    END SUBROUTINE open_read_file

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads global or variable-related attributes of type INTEGER (32-bit)
!------------------------------------------------------------------------------!
     SUBROUTINE get_attribute_int32( id, attribute_name, value, global,        &
                                     variable_name )

       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)            ::  attribute_name   !< attribute name
       CHARACTER(LEN=*), OPTIONAL  ::  variable_name    !< variable name

       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_var           !< variable id
       INTEGER(iwp), INTENT(INOUT) ::  value            !< read value

       LOGICAL, INTENT(IN) ::  global                   !< flag indicating global attribute
#if defined( __netcdf )

!
!--    Read global attribute
       IF ( global )  THEN
          nc_stat = NF90_GET_ATT( id, NF90_GLOBAL, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_int32 global', 522, attribute_name )
!
!--    Read attributes referring to a single variable. Therefore, first inquire
!--    variable id
       ELSE
          nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
          CALL handle_error( 'get_attribute_int32', 522, attribute_name )
          nc_stat = NF90_GET_ATT( id, id_var, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_int32', 522, attribute_name )
       ENDIF
#endif
    END SUBROUTINE get_attribute_int32

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads global or variable-related attributes of type INTEGER (8-bit)
!------------------------------------------------------------------------------!
     SUBROUTINE get_attribute_int8( id, attribute_name, value, global,         &
                                    variable_name )

       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)            ::  attribute_name   !< attribute name
       CHARACTER(LEN=*), OPTIONAL  ::  variable_name    !< variable name

       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_var           !< variable id
       INTEGER(KIND=1), INTENT(INOUT) ::  value         !< read value

       LOGICAL, INTENT(IN) ::  global                   !< flag indicating global attribute
#if defined( __netcdf )

!
!--    Read global attribute
       IF ( global )  THEN
          nc_stat = NF90_GET_ATT( id, NF90_GLOBAL, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_int8 global', 523, attribute_name )
!
!--    Read attributes referring to a single variable. Therefore, first inquire
!--    variable id
       ELSE
          nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
          CALL handle_error( 'get_attribute_int8', 523, attribute_name )
          nc_stat = NF90_GET_ATT( id, id_var, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_int8', 523, attribute_name )
       ENDIF
#endif
    END SUBROUTINE get_attribute_int8

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads global or variable-related attributes of type REAL
!------------------------------------------------------------------------------!
     SUBROUTINE get_attribute_real( id, attribute_name, value, global,         &
                                    variable_name )

       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)            ::  attribute_name   !< attribute name
       CHARACTER(LEN=*), OPTIONAL  ::  variable_name    !< variable name

       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_var           !< variable id

       LOGICAL, INTENT(IN) ::  global                   !< flag indicating global attribute

       REAL(wp), INTENT(INOUT)     ::  value            !< read value
#if defined( __netcdf )


!
!-- Read global attribute
       IF ( global )  THEN
          nc_stat = NF90_GET_ATT( id, NF90_GLOBAL, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_real global', 524, attribute_name )
!
!-- Read attributes referring to a single variable. Therefore, first inquire
!-- variable id
       ELSE
          nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
          CALL handle_error( 'get_attribute_real', 524, attribute_name )
          nc_stat = NF90_GET_ATT( id, id_var, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_real', 524, attribute_name )
       ENDIF
#endif
    END SUBROUTINE get_attribute_real

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads global or variable-related attributes of type CHARACTER
!> Remark: reading attributes of type NF_STRING return an error code 56 -
!> Attempt to convert between text & numbers.
!------------------------------------------------------------------------------!
     SUBROUTINE get_attribute_string( id, attribute_name, value, global,       &
                                      variable_name )

       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)                ::  attribute_name   !< attribute name
       CHARACTER(LEN=*), OPTIONAL      ::  variable_name    !< variable name
       CHARACTER(LEN=*), INTENT(INOUT) ::  value            !< read value

       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_var           !< variable id

       LOGICAL, INTENT(IN) ::  global                   !< flag indicating global attribute
#if defined( __netcdf )

!
!--    Read global attribute
       IF ( global )  THEN
          nc_stat = NF90_GET_ATT( id, NF90_GLOBAL, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_string global', 525, attribute_name )
!
!--    Read attributes referring to a single variable. Therefore, first inquire
!--    variable id
       ELSE
          nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
          CALL handle_error( 'get_attribute_string', 525, attribute_name )

          nc_stat = NF90_GET_ATT( id, id_var, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_string',525, attribute_name )

       ENDIF
#endif
    END SUBROUTINE get_attribute_string



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Get dimension array for a given dimension
!------------------------------------------------------------------------------!
     SUBROUTINE get_dimension_length( id, dim_len, variable_name )
#if defined( __netcdf )

       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)            ::  variable_name    !< dimension name
       CHARACTER(LEN=100)          ::  dum              !< dummy variable to receive return character

       INTEGER(iwp)                ::  dim_len          !< dimension size
       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_dim           !< dimension id

!
!--    First, inquire dimension ID
       nc_stat = NF90_INQ_DIMID( id, TRIM( variable_name ), id_dim )
       CALL handle_error( 'get_dimension_length', 526, variable_name )
!
!--    Inquire dimension length
       nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim, dum, LEN = dim_len )
       CALL handle_error( 'get_dimension_length', 526, variable_name )

#endif
    END SUBROUTINE get_dimension_length

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 1D integer variable from file.
!------------------------------------------------------------------------------!
     SUBROUTINE get_variable_1d_int( id, variable_name, var )

       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)            ::  variable_name    !< variable name

       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_var           !< dimension id

       INTEGER(iwp), DIMENSION(:), INTENT(INOUT) ::  var  !< variable to be read
#if defined( __netcdf )

!
!--    First, inquire variable ID
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
       CALL handle_error( 'get_variable_1d_int', 527, variable_name )
!
!--    Inquire dimension length
       nc_stat = NF90_GET_VAR( id, id_var, var )
       CALL handle_error( 'get_variable_1d_int', 527, variable_name )

#endif
    END SUBROUTINE get_variable_1d_int

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 1D float variable from file.
!------------------------------------------------------------------------------!
     SUBROUTINE get_variable_1d_real( id, variable_name, var )

       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)            ::  variable_name    !< variable name

       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_var           !< dimension id

       REAL(wp), DIMENSION(:), INTENT(INOUT) ::  var  !< variable to be read
#if defined( __netcdf )

!
!--    First, inquire variable ID
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
       CALL handle_error( 'get_variable_1d_real', 528, variable_name )
!
!--    Inquire dimension length
       nc_stat = NF90_GET_VAR( id, id_var, var )
       CALL handle_error( 'get_variable_1d_real', 528, variable_name )

#endif
    END SUBROUTINE get_variable_1d_real


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a time-dependent 1D float variable from file.
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_pr( id, variable_name, t, var )
#if defined( __netcdf )

       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)                      ::  variable_name    !< variable name

       INTEGER(iwp), INTENT(IN)              ::  id               !< file id
       INTEGER(iwp), DIMENSION(1:2)          ::  id_dim           !< dimension ids
       INTEGER(iwp)                          ::  id_var           !< dimension id
       INTEGER(iwp)                          ::  n_file           !< number of data-points in file along z dimension
       INTEGER(iwp), INTENT(IN)              ::  t                !< timestep number

       REAL(wp), DIMENSION(:), INTENT(INOUT) ::  var  !< variable to be read

!
!--    First, inquire variable ID
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Inquire dimension size of vertical dimension
       nc_stat = NF90_INQUIRE_VARIABLE( id, id_var, DIMIDS = id_dim )
       nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim(1), LEN = n_file )
!
!--    Read variable.
       nc_stat = NF90_GET_VAR( id, id_var, var,                                &
                               start = (/ 1,      t     /),                    &
                               count = (/ n_file, 1     /) )
       CALL handle_error( 'get_variable_pr', 529, variable_name )

#endif
    END SUBROUTINE get_variable_pr


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 2D REAL variable from a file. Reading is done processor-wise,
!> i.e. each core reads its own domain in slices along x.
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_2d_real( id, variable_name, var, is, ie, js, je )

       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       INTEGER(iwp)                  ::  i               !< running index along x direction
       INTEGER(iwp)                  ::  ie              !< start index for subdomain input along x direction
       INTEGER(iwp)                  ::  is              !< end index for subdomain input along x direction
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(iwp)                  ::  j               !< running index along y direction
       INTEGER(iwp)                  ::  je              !< start index for subdomain input along y direction
       INTEGER(iwp)                  ::  js              !< end index for subdomain input along y direction
       
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  tmp   !< temporary variable to read data from file according
                                                         !< to its reverse memory access
       REAL(wp), DIMENSION(:,:), INTENT(INOUT) ::  var   !< variable to be read
#if defined( __netcdf )
!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Check for collective read-operation and set respective NetCDF flags if
!--    required. 
       IF ( collective_read )  THEN
          nc_stat = NF90_VAR_PAR_ACCESS (id, id_var, NF90_COLLECTIVE)
       ENDIF
!
!--    Allocate temporary variable according to memory access on file. 
       ALLOCATE( tmp(is:ie,js:je) )
!
!--    Get variable
       nc_stat = NF90_GET_VAR( id, id_var, tmp,                                &
                               start = (/ is+1,      js+1 /),                  &
                               count = (/ ie-is + 1, je-js+1 /) )   
                               
       CALL handle_error( 'get_variable_2d_real', 530, variable_name )
!
!--    Resort data. Please note, dimension subscripts of var all start at 1. 
       DO  i = is, ie 
          DO  j = js, je 
             var(j-js+1,i-is+1) = tmp(i,j)
          ENDDO
       ENDDO
       
       DEALLOCATE( tmp )

#endif
    END SUBROUTINE get_variable_2d_real

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 2D 32-bit INTEGER variable from file. Reading is done processor-wise,
!> i.e. each core reads its own domain in slices along x.
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_2d_int32( id, variable_name, var, is, ie, js, je )

       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       INTEGER(iwp)                  ::  i               !< running index along x direction
       INTEGER(iwp)                  ::  ie              !< start index for subdomain input along x direction
       INTEGER(iwp)                  ::  is              !< end index for subdomain input along x direction
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(iwp)                  ::  j               !< running index along y direction
       INTEGER(iwp)                  ::  je              !< start index for subdomain input along y direction
       INTEGER(iwp)                  ::  js              !< end index for subdomain input along y direction
       
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE   ::  tmp  !< temporary variable to read data from file according
                                                            !< to its reverse memory access
       INTEGER(iwp), DIMENSION(:,:), INTENT(INOUT) ::  var  !< variable to be read
#if defined( __netcdf )
!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Check for collective read-operation and set respective NetCDF flags if
!--    required. 
       IF ( collective_read )  THEN
          nc_stat = NF90_VAR_PAR_ACCESS (id, id_var, NF90_COLLECTIVE)
       ENDIF
!
!--    Allocate temporary variable according to memory access on file. 
       ALLOCATE( tmp(is:ie,js:je) )
!
!--    Get variable
       nc_stat = NF90_GET_VAR( id, id_var, tmp,                                &
                               start = (/ is+1,      js+1 /),                  &
                               count = (/ ie-is + 1, je-js+1 /) )    
                               
       CALL handle_error( 'get_variable_2d_int32', 531, variable_name )                             
!
!--    Resort data. Please note, dimension subscripts of var all start at 1. 
       DO  i = is, ie 
          DO  j = js, je 
             var(j-js+1,i-is+1) = tmp(i,j)
          ENDDO
       ENDDO
       
       DEALLOCATE( tmp )

#endif
    END SUBROUTINE get_variable_2d_int32

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 2D 8-bit INTEGER variable from file. Reading is done processor-wise,
!> i.e. each core reads its own domain in slices along x.
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_2d_int8( id, variable_name, var, is, ie, js, je )

       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       INTEGER(iwp)                  ::  i               !< running index along x direction
       INTEGER(iwp)                  ::  ie              !< start index for subdomain input along x direction
       INTEGER(iwp)                  ::  is              !< end index for subdomain input along x direction
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(iwp)                  ::  j               !< running index along y direction
       INTEGER(iwp)                  ::  je              !< start index for subdomain input along y direction
       INTEGER(iwp)                  ::  js              !< end index for subdomain input along y direction
       
       INTEGER(KIND=1), DIMENSION(:,:), ALLOCATABLE   ::  tmp  !< temporary variable to read data from file according
                                                               !< to its reverse memory access
       INTEGER(KIND=1), DIMENSION(:,:), INTENT(INOUT) ::  var  !< variable to be read
#if defined( __netcdf )
!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Check for collective read-operation and set respective NetCDF flags if
!--    required. 
       IF ( collective_read )  THEN
          nc_stat = NF90_VAR_PAR_ACCESS (id, id_var, NF90_COLLECTIVE)
       ENDIF
!
!--    Allocate temporary variable according to memory access on file. 
       ALLOCATE( tmp(is:ie,js:je) )
!
!--    Get variable
       nc_stat = NF90_GET_VAR( id, id_var, tmp,                                &
                               start = (/ is+1,      js+1 /),                  &
                               count = (/ ie-is + 1, je-js+1 /) )   
                               
       CALL handle_error( 'get_variable_2d_int8', 532, variable_name )
!
!--    Resort data. Please note, dimension subscripts of var all start at 1. 
       DO  i = is, ie 
          DO  j = js, je 
             var(j-js+1,i-is+1) = tmp(i,j)
          ENDDO
       ENDDO
       
       DEALLOCATE( tmp )

#endif
    END SUBROUTINE get_variable_2d_int8


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 3D 8-bit INTEGER variable from file.
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_3d_int8( id, variable_name, var, is, ie, js, je,   &
                                     ks, ke )

       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       INTEGER(iwp)                  ::  i               !< index along x direction
       INTEGER(iwp)                  ::  ie              !< start index for subdomain input along x direction
       INTEGER(iwp)                  ::  is              !< end index for subdomain input along x direction
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(iwp)                  ::  j               !< index along y direction
       INTEGER(iwp)                  ::  je              !< start index for subdomain input along y direction
       INTEGER(iwp)                  ::  js              !< end index for subdomain input along y direction
       INTEGER(iwp)                  ::  k               !< index along any 3rd dimension
       INTEGER(iwp)                  ::  ke              !< start index of 3rd dimension
       INTEGER(iwp)                  ::  ks              !< end index of 3rd dimension
       
       INTEGER(KIND=1), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp  !< temporary variable to read data from file according 
                                                                 !< to its reverse memory access

       INTEGER(KIND=1), DIMENSION(:,:,:), INTENT(INOUT) ::  var  !< variable to be read
#if defined( __netcdf )

!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )   
!
!--    Check for collective read-operation and set respective NetCDF flags if
!--    required. 
       IF ( collective_read )  THEN
          nc_stat = NF90_VAR_PAR_ACCESS (id, id_var, NF90_COLLECTIVE)
       ENDIF
!
!--    Allocate temporary variable according to memory access on file. 
       ALLOCATE( tmp(is:ie,js:je,ks:ke) )
!
!--    Get variable
       nc_stat = NF90_GET_VAR( id, id_var, tmp,                                &
                               start = (/ is+1,    js+1,    ks+1 /),           &
                               count = (/ ie-is+1, je-js+1, ke-ks+1 /) )                              

       CALL handle_error( 'get_variable_3d_int8', 533, variable_name )                               
!
!--    Resort data. Please note, dimension subscripts of var all start at 1. 
       DO  i = is, ie 
          DO  j = js, je
             DO  k = ks, ke
                var(k-ks+1,j-js+1,i-is+1) = tmp(i,j,k)
             ENDDO
          ENDDO
       ENDDO
       
       DEALLOCATE( tmp )

#endif
    END SUBROUTINE get_variable_3d_int8


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 3D float variable from file.
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_3d_real( id, variable_name, var, is, ie, js, je,   &
                                     ks, ke )

       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       INTEGER(iwp)                  ::  i               !< index along x direction
       INTEGER(iwp)                  ::  ie              !< start index for subdomain input along x direction
       INTEGER(iwp)                  ::  is              !< end index for subdomain input along x direction
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(iwp)                  ::  j               !< index along y direction
       INTEGER(iwp)                  ::  je              !< start index for subdomain input along y direction
       INTEGER(iwp)                  ::  js              !< end index for subdomain input along y direction
       INTEGER(iwp)                  ::  k               !< index along any 3rd dimension
       INTEGER(iwp)                  ::  ke              !< start index of 3rd dimension
       INTEGER(iwp)                  ::  ks              !< end index of 3rd dimension
       
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp !< temporary variable to read data from file according 
                                                         !< to its reverse memory access

       REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::  var  !< variable to be read
#if defined( __netcdf )

!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )  
!
!--    Check for collective read-operation and set respective NetCDF flags if
!--    required. 
       IF ( collective_read )  THEN
          nc_stat = NF90_VAR_PAR_ACCESS (id, id_var, NF90_COLLECTIVE)
       ENDIF
!
!--    Allocate temporary variable according to memory access on file. 
       ALLOCATE( tmp(is:ie,js:je,ks:ke) )
!
!--    Get variable
       nc_stat = NF90_GET_VAR( id, id_var, tmp,                                &
                               start = (/ is+1,    js+1,    ks+1 /),           &
                               count = (/ ie-is+1, je-js+1, ke-ks+1 /) )   
                               
       CALL handle_error( 'get_variable_3d_real', 534, variable_name )
!
!--    Resort data. Please note, dimension subscripts of var all start at 1. 
       DO  i = is, ie 
          DO  j = js, je
             DO  k = ks, ke
                var(k-ks+1,j-js+1,i-is+1) = tmp(i,j,k)
             ENDDO
          ENDDO
       ENDDO
       
       DEALLOCATE( tmp )

#endif
    END SUBROUTINE get_variable_3d_real

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 3D float array from file.
!------------------------------------------------------------------------------!
!     SUBROUTINE get_variable_3d_real_v( id, variable_name, is, ie, js, je, var )
! 
!        USE indices
!        USE pegrid
! 
!        IMPLICIT NONE
! 
!        CHARACTER(LEN=*)              ::  variable_name   !< variable name
! 
!        INTEGER(iwp), INTENT(IN)      ::  is,ie           !< index range along x direction
!        INTEGER(iwp), INTENT(IN)      ::  id              !< file id
!        INTEGER(iwp)                  ::  id_var          !< variable id
!        INTEGER(iwp), INTENT(IN)      ::  js,je           !< index range along y direction
!        INTEGER(iwp)                  ::  n3              !< number of data-points along 3rd dimension
! 
!        INTEGER(iwp)                  ::  i,j,k
!        INTEGER(iwp), DIMENSION(3)    ::  id_dim
! 
!        REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::  var         !< variable to be read
! #if defined( __netcdf )
! !
! !--    Inside the ...static NetCDF files, the array is stored as float.
! !--    Therefore single precision is sufficiant for the temporary array
! 
!        REAL(sp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_var     !< temporary array to read NetCDF data in i,j,k direction
! 
! !kk    Please check, if it is time consuming to do the inquire every time
! !
! !--    Inquire variable id
!        nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
! !
! !--    Get length of third dimension, required for the count parameter.
! !--    Therefore, first inquired dimension ids
!        nc_stat = NF90_INQUIRE_VARIABLE( id, id_var, DIMIDS = id_dim )
!        nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim(3), LEN = n3 )
! 
! !
! !--    Check for collective read-operation and set respective NetCDF flags if
! !--    required. 
!        IF ( collective_read )  THEN
!           nc_stat = NF90_VAR_PAR_ACCESS (id, id_var, NF90_COLLECTIVE)
!        ENDIF
! 
! !
! !--    Allocate temporary array ro read NetCDF data in i,j,k direction
! 
!        ALLOCATE(tmp_var(is:ie,js:je,n3))
! !
! !--    Get variable
! !--    Read complete local 3-D array in oone call
! 
!        nc_stat = NF90_GET_VAR( id, id_var, tmp_var,                            &
!                                start = (/ is+1, js+1, 1 /),                    &
!                                count = (/ ie-is+1, je-js+1, n3 /) )
! 
!        CALL handle_error( 'get_variable_3d_real', 532 )
! 
! !
! !--    Resort data in k,j,i direction
! 
!        DO i=is,ie
!           DO j=js,je
!              DO K=1,n3
!                 var (k,j-js+1,i-is+1) = tmp_var(i,j,k)
!              END DO
!           END DO
!        END DO
! 
!        DEALLOCATE(tmp_var)
! 
! #endif
!     END SUBROUTINE get_variable_3d_real_v


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 4D float variable from file. 
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_4d_real( id, variable_name, var, is, ie, js, je,   &
                                     k1s, k1e, k2s, k2e )

       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       INTEGER(iwp)                  ::  i               !< index along x direction
       INTEGER(iwp)                  ::  ie              !< start index for subdomain input along x direction
       INTEGER(iwp)                  ::  is              !< end index for subdomain input along x direction
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(iwp)                  ::  j               !< index along y direction
       INTEGER(iwp)                  ::  je              !< start index for subdomain input along y direction
       INTEGER(iwp)                  ::  js              !< end index for subdomain input along y direction
       INTEGER(iwp)                  ::  k1              !< index along 3rd direction
       INTEGER(iwp)                  ::  k1e             !< start index for 3rd dimension
       INTEGER(iwp)                  ::  k1s             !< end index for 3rd dimension
       INTEGER(iwp)                  ::  k2              !< index along 4th direction
       INTEGER(iwp)                  ::  k2e             !< start index for 4th dimension
       INTEGER(iwp)                  ::  k2s             !< end index for 4th dimension

       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE   ::  tmp  !< temporary variable to read data from file according
                                                            !< to its reverse memory access
       REAL(wp), DIMENSION(:,:,:,:), INTENT(INOUT) ::  var  !< variable to be read
#if defined( __netcdf )

!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Check for collective read-operation and set respective NetCDF flags if
!--    required. 
       IF ( collective_read )  THEN
          nc_stat = NF90_VAR_PAR_ACCESS (id, id_var, NF90_COLLECTIVE)
       ENDIF
!
!--    Allocate temporary variable according to memory access on file. 
       ALLOCATE( tmp(is:ie,js:je,k1s:k1e,k2s:k2e) )
!
!--    Get variable
       nc_stat = NF90_GET_VAR( id, id_var, tmp,                                &
                               start = (/ is+1,    js+1,   k1s+1, k2s+1 /),    &
                               count = (/ ie-is+1, je-js+1,                    &
                                          k1e-k1s+1, k2e-k2s+1 /) )

       CALL handle_error( 'get_variable_4d_real', 535, variable_name )
!
!--    Resort data. Please note, dimension subscripts of var all start at 1. 
       DO  i = is, ie 
          DO  j = js, je
             DO  k1 = k1s, k1e
                DO  k2 = k2s, k2e
                   var(k2-k2s+1,k1-k1s+1,j-js+1,i-is+1) = tmp(i,j,k1,k2)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       
       DEALLOCATE( tmp )
#endif
    END SUBROUTINE get_variable_4d_real



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 3D float variables from dynamic driver, such as time-dependent xy-, 
!> xz- or yz-boundary data as well as 3D initialization data. Please note, 
!> the passed arguments are start indices and number of elements in each 
!> dimension, which is in contrast to the other 3d versions where start- and 
!> end indices are passed. The different handling of 3D dynamic variables is 
!> due to its asymmetry for the u- and v component.
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_3d_real_dynamic( id, variable_name, var,           &
                            i1s, i2s, i3s, count_1, count_2, count_3, dynamic)
                                
       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       LOGICAL                       ::  dynamic         !< additional flag just used to select correct overloaded routine from interface block 
       
       INTEGER(iwp)                  ::  count_1         !< number of elements to be read along 1st dimension (with respect to file)
       INTEGER(iwp)                  ::  count_2         !< number of elements to be read along 2nd dimension (with respect to file)
       INTEGER(iwp)                  ::  count_3         !< number of elements to be read along 3rd dimension (with respect to file)
       INTEGER(iwp)                  ::  i1              !< running index along 1st dimension on file
       INTEGER(iwp)                  ::  i1s             !< start index for subdomain input along 1st dimension (with respect to file)
       INTEGER(iwp)                  ::  i2              !< running index along 2nd dimension on file       
       INTEGER(iwp)                  ::  i2s             !< start index for subdomain input along 2nd dimension (with respect to file)
       INTEGER(iwp)                  ::  i3              !< running index along 3rd dimension on file 
       INTEGER(iwp)                  ::  i3s             !< start index of 3rd dimension, in dynamic file this is either time (2D boundary) or z (3D)
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(iwp)                  ::  lb1             !< lower bound of 1st dimension (with respect to file)
       INTEGER(iwp)                  ::  lb2             !< lower bound of 2nd dimension (with respect to file)
       INTEGER(iwp)                  ::  lb3             !< lower bound of 3rd dimension (with respect to file)
       INTEGER(iwp)                  ::  ub1             !< upper bound of 1st dimension (with respect to file)
       INTEGER(iwp)                  ::  ub2             !< upper bound of 2nd dimension (with respect to file)
       INTEGER(iwp)                  ::  ub3             !< upper bound of 3rd dimension (with respect to file)

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp !< temporary variable to read data from file according
                                                         !< to its reverse memory access
       
       REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::  var !< input variable
       
#if defined( __netcdf )
!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Check for collective read-operation and set respective NetCDF flags if
!--    required. 
       IF ( collective_read )  THEN
          nc_stat = NF90_VAR_PAR_ACCESS (id, id_var, NF90_COLLECTIVE)
       ENDIF    
!
!--    Allocate temporary variable according to memory access on file. 
!--    Therefore, determine dimension bounds of input array. 
       lb1 = LBOUND(var,3)
       ub1 = UBOUND(var,3)
       lb2 = LBOUND(var,2)
       ub2 = UBOUND(var,2)
       lb3 = LBOUND(var,1)
       ub3 = UBOUND(var,1)
       ALLOCATE( tmp(lb1:ub1,lb2:ub2,lb3:ub3) )
!
!--    Get variable
       nc_stat = NF90_GET_VAR( id, id_var, tmp,                                &
                               start = (/ i1s,     i2s,     i3s /),            &
                               count = (/ count_1, count_2, count_3 /) )

       CALL handle_error( 'get_variable_3d_real_dynamic', 536, variable_name )
!
!--    Resort data. Please note, dimension subscripts of var all start at 1. 
       DO  i3 = lb3, ub3
          DO i2 = lb2, ub2
             DO  i1 = lb1, ub1
                var(i3,i2,i1) = tmp(i1,i2,i3)
             ENDDO
          ENDDO
       ENDDO
       
       DEALLOCATE( tmp )       
#endif
    END SUBROUTINE get_variable_3d_real_dynamic



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Inquires the number of variables in a file
!------------------------------------------------------------------------------!
    SUBROUTINE inquire_num_variables( id, num_vars )

       USE indices
       USE pegrid

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp), INTENT(INOUT)   ::  num_vars        !< number of variables in a file
#if defined( __netcdf )

       nc_stat = NF90_INQUIRE( id, NVARIABLES = num_vars )
       CALL handle_error( 'inquire_num_variables', 537 )

#endif
    END SUBROUTINE inquire_num_variables


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Inquires the variable names belonging to a file.
!------------------------------------------------------------------------------!
    SUBROUTINE inquire_variable_names( id, var_names )

       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*), DIMENSION(:), INTENT(INOUT) ::  var_names   !< return variable - variable names
       INTEGER(iwp)                                  ::  i           !< loop variable
       INTEGER(iwp), INTENT(IN)                      ::  id          !< file id
       INTEGER(iwp)                                  ::  num_vars    !< number of variables (unused return parameter)
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE       ::  varids      !< dummy array to strore variable ids temporarily
#if defined( __netcdf )

       ALLOCATE( varids(1:SIZE(var_names)) )
       nc_stat = NF90_INQ_VARIDS( id, NVARS = num_vars, VARIDS = varids )
       CALL handle_error( 'inquire_variable_names', 538 )

       DO  i = 1, SIZE(var_names)
          nc_stat = NF90_INQUIRE_VARIABLE( id, varids(i), NAME = var_names(i) )
          CALL handle_error( 'inquire_variable_names', 538 )
       ENDDO

       DEALLOCATE( varids )
#endif
    END SUBROUTINE inquire_variable_names

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Prints out a text message corresponding to the current status.
!------------------------------------------------------------------------------!
    SUBROUTINE handle_error( routine_name, errno, name )

       USE control_parameters,                                                 &
           ONLY:  message_string

       IMPLICIT NONE

       CHARACTER(LEN=6) ::  message_identifier
       CHARACTER(LEN=*) ::  routine_name
       CHARACTER(LEN=*), OPTIONAL ::  name

       INTEGER(iwp) ::  errno
#if defined( __netcdf )
       
       IF ( nc_stat /= NF90_NOERR )  THEN

          WRITE( message_identifier, '(''NC'',I4.4)' )  errno
          
          IF ( PRESENT( name ) )  THEN
             message_string = "Problem reading attribute/variable - " //       &
                              TRIM(name) // ": " //                            &
                              TRIM( NF90_STRERROR( nc_stat ) )
          ELSE
             message_string = TRIM( NF90_STRERROR( nc_stat ) )
          ENDIF

          CALL message( routine_name, message_identifier, 2, 2, 0, 6, 1 )

       ENDIF

#endif
    END SUBROUTINE handle_error


 END MODULE netcdf_data_input_mod
