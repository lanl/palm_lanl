!> @file src/grid.f90
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
! Copyright 2017-2018 Leibniz Universitaet Hannover
! Copyright 2017-2018 Deutscher Wetterdienst Offenbach
!------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: grid.f90 2718 2018-01-02 08:49:38Z maronga $
! Initial revision
!
! 
!
! Authors:
! --------
! @author Eckhard Kadasch
!
!------------------------------------------------------------------------------!
! Description:
! ------------
!> The grid module contains all variables and routines to specify and work with
!> the numerical grids in INIFOR. By convention, all angles are stored in
!> radians.
!------------------------------------------------------------------------------!

 MODULE grid

    USE control
    USE defs,                                                                  &
        ONLY:  DATE, EARTH_RADIUS, TO_RADIANS, TO_DEGREES, PI, dp, hp, sp,     &
               SNAME, LNAME, PATH, FORCING_FREQ, WATER_ID, FILL_ITERATIONS,    &
               BETA, P_SL, T_SL, BETA, RD, G, P_REF, RD_PALM, CP_PALM, RHO_L
    USE io,                                                                    &
        ONLY:  get_netcdf_variable, get_netcdf_attribute,                      &
               parse_command_line_arguments
    USE netcdf,                                                                &
        ONLY:  NF90_MAX_NAME, NF90_MAX_VAR_DIMS
    USE transform,                                                             &
        ONLY:  rotate_to_cosmo, find_horizontal_neighbours,                    &
               compute_horizontal_interp_weights,                              &
               find_vertical_neighbours_and_weights, interpolate_2d,           &
               gamma_from_hemisphere, phic_to_phin, lamc_to_lamn, average_2d,  &
               project, centre_velocities, phi2phirot, rla2rlarot, uv2uvrot
    USE types
    USE util
    
    IMPLICIT NONE
    
    SAVE
    
    REAL(dp) ::  phi_equat = 0.0_dp    !< latitude of rotated equator of COSMO-DE grid [rad]
    REAL(dp) ::  phi_n     = 0.0_dp    !< latitude of rotated pole of COSMO-DE grid [rad]
    REAL(dp) ::  lambda_n  = 0.0_dp    !< longitude of rotaded pole of COSMO-DE grid [rad]
    REAL(dp) ::  phi_c     = 0.0_dp    !< rotated-grid latitude of the center of the PALM domain [rad]
    REAL(dp) ::  lambda_c  = 0.0_dp    !< rotated-grid longitude of the centre of the PALM domain [rad]
    REAL(dp) ::  phi_cn    = 0.0_dp    !< latitude of the rotated pole relative to the COSMO-DE grid [rad]
    REAL(dp) ::  lambda_cn = 0.0_dp    !< longitude of the rotated pole relative to the COSMO-DE grid [rad]
    REAL(dp) ::  gam       = 0.0_dp    !< angle for working around phirot2phi/rlarot2rla bug
    REAL(dp) ::  dx        = 0.0_dp    !< PALM-4U grid spacing in x direction [m]
    REAL(dp) ::  dy        = 0.0_dp    !< PALM-4U grid spacing in y direction [m]
    REAL(dp) ::  dz        = 0.0_dp    !< PALM-4U grid spacing in z direction [m]
    REAL(dp) ::  dxi       = 0.0_dp    !< inverse PALM-4U grid spacing in x direction [m^-1]
    REAL(dp) ::  dyi       = 0.0_dp    !< inverse PALM-4U grid spacing in y direction [m^-1]
    REAL(dp) ::  dzi       = 0.0_dp    !< inverse PALM-4U grid spacing in z direction [m^-1]
    REAL(dp) ::  lx        = 0.0_dp    !< PALM-4U domain size in x direction [m]
    REAL(dp) ::  ly        = 0.0_dp    !< PALM-4U domain size in y direction [m]
    REAL(dp) ::  lz        = 0.0_dp    !< PALM-4U domain size in z direction [m]
    REAL(dp) ::  ug        = 0.0_dp    !< geostrophic wind in x direction [m/s]
    REAL(dp) ::  vg        = 0.0_dp    !< geostrophic wind in y direction [m/s]
    REAL(dp) ::  p0        = 0.0_dp    !< PALM-4U surface pressure, at z0 [Pa]
    REAL(dp) ::  x0        = 0.0_dp    !< x coordinate of PALM-4U Earth tangent [m] 
    REAL(dp) ::  y0        = 0.0_dp    !< y coordinate of PALM-4U Earth tangent [m] 
    REAL(dp) ::  z0        = 0.0_dp    !< Elevation of the PALM-4U domain above sea level [m]
    REAL(dp) ::  lonmin    = 0.0_dp    !< Minimunm longitude of COSMO-DE's rotated-pole grid
    REAL(dp) ::  lonmax    = 0.0_dp    !< Maximum longitude of COSMO-DE's rotated-pole grid
    REAL(dp) ::  latmin    = 0.0_dp    !< Minimunm latitude of COSMO-DE's rotated-pole grid
    REAL(dp) ::  latmax    = 0.0_dp    !< Maximum latitude of COSMO-DE's rotated-pole grid
    REAL(dp) ::  latitude  = 0.0_dp    !< geograpohical latitude of the PALM-4U origin, from inipar namelist [deg]
    REAL(dp) ::  longitude = 0.0_dp    !< geograpohical longitude of the PALM-4U origin, from inipar namelist [deg]
    REAL(dp) ::  origin_lat= 0.0_dp    !< geograpohical latitude of the PALM-4U origin, from static driver netCDF file [deg]
    REAL(dp) ::  origin_lon= 0.0_dp    !< geograpohical longitude of the PALM-4U origin, from static driver netCDF file [deg]
    REAL(dp) ::  end_time  = 0.0_dp    !< PALM-4U simulation time [s]

    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  hhl             !< heights of half layers (cell faces) above sea level in COSMO-DE, read in from external file
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  hfl             !< heights of full layers (cell centres) above sea level in COSMO-DE, computed as arithmetic average of hhl
    REAL(dp), DIMENSION(:), ALLOCATABLE, TARGET     ::  depths          !< COSMO-DE's TERRA-ML soil layer depths
    REAL(dp), DIMENSION(:), ALLOCATABLE, TARGET     ::  d_depth_rho_inv !< COSMO-DE's TERRA-ML soil layer thicknesses
    REAL(dp), DIMENSION(:), ALLOCATABLE, TARGET     ::  rlon            !< longitudes of COSMO-DE's rotated-pole grid
    REAL(dp), DIMENSION(:), ALLOCATABLE, TARGET     ::  rlat            !< latitudes of COSMO-DE's rotated-pole grid
    REAL(dp), DIMENSION(:), ALLOCATABLE, TARGET     ::  time

    INTEGER(hp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  soiltyp      !< COSMO-DE soil type map
    INTEGER ::  i     !< indexing variable
    INTEGER ::  imin, imax, jmin, jmax !< index ranges for profile averaging
    INTEGER ::  k     !< indexing variable
    INTEGER ::  nt    !< number of output time steps
    INTEGER ::  nx    !< number of PALM-4U grid points in x direction
    INTEGER ::  ny    !< number of PALM-4U grid points in y direction
    INTEGER ::  nz    !< number of PALM-4U grid points in z direction
    INTEGER ::  nlon  !< number of longitudal points in target grid (COSMO-DE)
    INTEGER ::  nlat  !< number of latitudal points in target grid (COSMO-DE)
    INTEGER ::  nlev  !< number of levels in target grid (COSMO-DE)
    INTEGER ::  layers !< number of COSMO-DE soil layers
    INTEGER ::  start_hour_flow         !< start of flow forcing in number of hours relative to start_date
    INTEGER ::  start_hour_soil         !< start of soil forcing in number of hours relative to start_date, typically equals start_hour_flow
    INTEGER ::  start_hour_radiation    !< start of radiation forcing in number of hours relative to start_date, 0 to 2 hours before start_hour_flow to reconstruct hourly averages from one- to three hourly averages of the input data
    INTEGER ::  start_hour_soilmoisture !< start of forcing for the soil moisture spin-up in number of hours relative to start_date, typically -672 (-4 weeks)
    INTEGER ::  end_hour  !< simulation time in hours
    INTEGER ::  end_hour_soilmoisture  !< end of soil moisture spin-up in hours relative to start_hour_flow
    INTEGER ::  step_hour !< number of hours between forcing time steps

    LOGICAL ::  init_variables_required
    LOGICAL ::  boundary_variables_required

    INTEGER ::  n_invar = 0 !< number of variables in the input variable table
    INTEGER ::  n_outvar = 0 !< number of variables in the output variable table
    TYPE(nc_var), ALLOCATABLE, TARGET ::  input_var_table(:)  !< table of input variables
    TYPE(nc_var), ALLOCATABLE, TARGET ::  output_var_table(:) !< table of input variables
    TYPE(nc_var)                      ::  cosmo_var           !< COSMO-DE dummy variable, used for reading HHL, rlon, rlat

    TYPE(grid_definition), TARGET ::  palm_grid                       !< PALM-4U grid in the target system (COSMO-DE rotated-pole)
    TYPE(grid_definition), TARGET ::  palm_intermediate               !< PALM-4U grid with coarse vertical grid wiht levels interpolated from COSMO-DE grid
    TYPE(grid_definition), TARGET ::  cosmo_grid                      !< target system (COSMO-DE rotated-pole)
    TYPE(grid_definition), TARGET ::  scalars_east_grid               !< 
    TYPE(grid_definition), TARGET ::  scalars_west_grid               !< 
    TYPE(grid_definition), TARGET ::  scalars_north_grid              !< 
    TYPE(grid_definition), TARGET ::  scalars_south_grid              !< 
    TYPE(grid_definition), TARGET ::  scalars_top_grid                !< 
    TYPE(grid_definition), TARGET ::  scalars_east_intermediate       !< 
    TYPE(grid_definition), TARGET ::  scalars_west_intermediate       !< 
    TYPE(grid_definition), TARGET ::  scalars_north_intermediate      !< 
    TYPE(grid_definition), TARGET ::  scalars_south_intermediate      !< 
    TYPE(grid_definition), TARGET ::  scalars_top_intermediate        !< 
    TYPE(grid_definition), TARGET ::  u_initial_grid                  !< 
    TYPE(grid_definition), TARGET ::  u_east_grid                     !< 
    TYPE(grid_definition), TARGET ::  u_west_grid                     !< 
    TYPE(grid_definition), TARGET ::  u_north_grid                    !< 
    TYPE(grid_definition), TARGET ::  u_south_grid                    !< 
    TYPE(grid_definition), TARGET ::  u_top_grid                      !< 
    TYPE(grid_definition), TARGET ::  u_initial_intermediate          !<
    TYPE(grid_definition), TARGET ::  u_east_intermediate             !< 
    TYPE(grid_definition), TARGET ::  u_west_intermediate             !< 
    TYPE(grid_definition), TARGET ::  u_north_intermediate            !< 
    TYPE(grid_definition), TARGET ::  u_south_intermediate            !< 
    TYPE(grid_definition), TARGET ::  u_top_intermediate              !< 
    TYPE(grid_definition), TARGET ::  v_initial_grid                  !< 
    TYPE(grid_definition), TARGET ::  v_east_grid                     !< 
    TYPE(grid_definition), TARGET ::  v_west_grid                     !< 
    TYPE(grid_definition), TARGET ::  v_north_grid                    !< 
    TYPE(grid_definition), TARGET ::  v_south_grid                    !< 
    TYPE(grid_definition), TARGET ::  v_top_grid                      !< 
    TYPE(grid_definition), TARGET ::  v_initial_intermediate          !<
    TYPE(grid_definition), TARGET ::  v_east_intermediate             !< 
    TYPE(grid_definition), TARGET ::  v_west_intermediate             !< 
    TYPE(grid_definition), TARGET ::  v_north_intermediate            !< 
    TYPE(grid_definition), TARGET ::  v_south_intermediate            !< 
    TYPE(grid_definition), TARGET ::  v_top_intermediate              !< 
    TYPE(grid_definition), TARGET ::  w_initial_grid                  !< 
    TYPE(grid_definition), TARGET ::  w_east_grid                     !< 
    TYPE(grid_definition), TARGET ::  w_west_grid                     !< 
    TYPE(grid_definition), TARGET ::  w_north_grid                    !< 
    TYPE(grid_definition), TARGET ::  w_south_grid                    !< 
    TYPE(grid_definition), TARGET ::  w_top_grid                      !< 
    TYPE(grid_definition), TARGET ::  w_initial_intermediate          !<
    TYPE(grid_definition), TARGET ::  w_east_intermediate             !< 
    TYPE(grid_definition), TARGET ::  w_west_intermediate             !< 
    TYPE(grid_definition), TARGET ::  w_north_intermediate            !< 
    TYPE(grid_definition), TARGET ::  w_south_intermediate            !< 
    TYPE(grid_definition), TARGET ::  w_top_intermediate              !< 
    TYPE(grid_definition), TARGET ::  scalar_profile_grid             !< 
    TYPE(grid_definition), TARGET ::  scalar_profile_intermediate     !< 
    TYPE(grid_definition), TARGET ::  w_profile_grid                  !< 
    TYPE(grid_definition), TARGET ::  w_profile_intermediate          !< 

    TYPE(io_group), ALLOCATABLE, TARGET ::  io_group_list(:)  !< List of I/O groups, which group together output variables that share the same input variable
 
    NAMELIST /inipar/ nx, ny, nz, dx, dy, dz, longitude, latitude
    NAMELIST /d3par/  end_time
    
    CHARACTER(LEN=LNAME) ::  nc_source_text     = ''  !< Text describing the source of the output data, e.g. 'COSMO-DE analysis from ...'
    CHARACTER(LEN=DATE)  ::  start_date         = ''  !< String of the FORMAT YYYYMMDDHH indicating the start of the intended PALM-4U simulation
    CHARACTER(LEN=PATH)  ::  hhl_file           = ''  !< Path to the file containing the COSMO-DE HHL variable (height of half layers, i.e. vertical cell faces)
    CHARACTER(LEN=PATH)  ::  namelist_file      = ''  !< Path to the PALM-4U namelist file
    CHARACTER(LEN=PATH)  ::  static_driver_file = ''  !< Path to the file containing the COSMO-DE SOILTYP variable (map of COSMO-DE soil types)
    CHARACTER(LEN=PATH)  ::  soiltyp_file       = ''  !< Path to the file containing the COSMO-DE SOILTYP variable (map of COSMO-DE soil types)
    CHARACTER(LEN=PATH)  ::  input_path         = ''  !< Path to the input data file directory
    CHARACTER(LEN=PATH), ALLOCATABLE, DIMENSION(:) ::  flow_files
    CHARACTER(LEN=PATH), ALLOCATABLE, DIMENSION(:) ::  soil_moisture_files
    CHARACTER(LEN=PATH), ALLOCATABLE, DIMENSION(:) ::  soil_files
    CHARACTER(LEN=PATH), ALLOCATABLE, DIMENSION(:) ::  radiation_files
    CHARACTER(LEN=SNAME) ::  input_prefix         !< Prefix of input files, e.g. 'laf' for COSMO-DE analyses
    CHARACTER(LEN=SNAME) ::  flow_suffix          !< Suffix of flow input files, e.g. 'flow'
    CHARACTER(LEN=SNAME) ::  soil_suffix          !< Suffix of soil input files, e.g. 'soil'
    CHARACTER(LEN=SNAME) ::  radiation_suffix     !< Suffix of radiation input files, e.g. 'radiation'
    CHARACTER(LEN=SNAME) ::  soilmoisture_suffix  !< Suffix of input files for soil moisture spin-up, e.g. 'soilmoisture'
    CHARACTER(LEN=SNAME) ::  mode                 !< INIFOR's initialization mode, 'profile' or 'volume'
                          
    TYPE(nc_file) ::  output_file

 CONTAINS
    
    SUBROUTINE setup_parameters()

!
!------------------------------------------------------------------------------
! Section 1: Define default parameters
!------------------------------------------------------------------------------
       start_date = '2013072100'
       end_hour = 2
       start_hour_soil = -2
       start_hour_soilmoisture = - (4 * 7 * 24) - 2

       lonmin =  -5.0_dp * TO_RADIANS
       lonmax =   5.5_dp * TO_RADIANS
       latmin =  -5.0_dp * TO_RADIANS
       latmax =   6.5_dp * TO_RADIANS

       ! COSMO-DE default rotated pole
       phi_n     =   40.0_dp * TO_RADIANS
       phi_equat =   50.0_dp * TO_RADIANS
       lambda_n  = -170.0_dp * TO_RADIANS

       ! COMSMO-DE soil layers
       layers = 8    
       ALLOCATE( depths(layers), d_depth_rho_inv(layers) )
       depths = (/0.005_dp, 0.02_dp, 0.06_dp, 0.18_dp, 0.54_dp, 1.62_dp, 4.86_dp, 14.58_dp/)
       d_depth_rho_inv = 1.0_dp / &
          ( (/0.01_dp, 0.02_dp, 0.06_dp, 0.18_dp, 0.54_dp, 1.62_dp, 4.86_dp, 14.58_dp/) * RHO_L )

       ! Defaultmain centre (_c) of the PALM-4U grid in the geographical system (_g)
       origin_lat = 52.325079_dp * TO_RADIANS ! south-west of Berlin, origin used for the Dec 2017 showcase simulation
       origin_lon = 13.082744_dp * TO_RADIANS
       z0 = 35.0_dp

       ! Default atmospheric parameters
       ug = 0.0_dp
       vg = 0.0_dp
       p0 = P_SL

       ! Parameters for file names
       start_hour_flow = 0
       start_hour_soil = 0
       start_hour_radiation = 0
       start_hour_soilmoisture = start_hour_flow - 2
       end_hour_soilmoisture = start_hour_flow
       step_hour = 1
       input_prefix = 'laf'  ! 'laf' for COSMO-DE analyses
       flow_suffix = '-flow'
       soil_suffix = '-soil'
       radiation_suffix = '-rad'
       soilmoisture_suffix = '-soilmoisture'
!
!------------------------------------------------------------------------------
! Section 2: Read command-line arguments, namelist, and grid parameters
!------------------------------------------------------------------------------

       ! Set default paths
       input_path         = './'
       hhl_file           = ''
       soiltyp_file       = ''
       namelist_file      = './namelist'
       output_file % name = './palm-4u-input.nc'

       ! Update default file names and domain centre
       CALL parse_command_line_arguments( start_date, hhl_file, soiltyp_file,  &
               static_driver_file, input_path, output_file % name,             &
               namelist_file, ug, vg, p0, z0, mode )

       init_variables_required     = .TRUE.
       boundary_variables_required = (TRIM(mode) .NE. 'profile')

       CALL normalize_path(input_path)
       IF (TRIM(hhl_file) == '')  hhl_file = TRIM(input_path) // 'hhl.nc'
       IF (TRIM(soiltyp_file) == '')  soiltyp_file = TRIM(input_path) // 'soil.nc'

       CALL report('setup_parameters', "       data path: " // TRIM(input_path))
       CALL report('setup_parameters', "        hhl file: " // TRIM(hhl_file))
       CALL report('setup_parameters', "    soiltyp file: " // TRIM(soiltyp_file))
       CALL report('setup_parameters', "   namelist file: " // TRIM(namelist_file))
       CALL report('setup_parameters', "output data file: " // TRIM(output_file % name))

 CALL run_control('time', 'init')
       ! Read in namelist parameters
       OPEN(10, FILE=namelist_file)
       READ(10, NML=inipar) ! nx, ny, nz, dx, dy, dz
       READ(10, NML=d3par)  ! end_time
       CLOSE(10)
 CALL run_control('time', 'read')

       end_hour = CEILING(end_time / FORCING_FREQ)

       ! Generate input file lists
       CALL input_file_list(start_date, start_hour_flow, end_hour, step_hour,  &
          input_path, input_prefix, flow_suffix, flow_files)
       CALL input_file_list(start_date, start_hour_soil, end_hour, step_hour,  &
          input_path, input_prefix, soil_suffix, soil_files)
       CALL input_file_list(start_date, start_hour_radiation, end_hour, step_hour, &
          input_path, input_prefix, radiation_suffix, radiation_files)
       CALL input_file_list(start_date, start_hour_soilmoisture, end_hour_soilmoisture, step_hour, &
          input_path, input_prefix, soilmoisture_suffix, soil_moisture_files)

!
!------------------------------------------------------------------------------
! Section 3: Check for consistency
!------------------------------------------------------------------------------
       IF (dx*dy*dz .EQ. 0.0_dp)  THEN
          message = "Grid cells have zero volume. Grid spacings are probably"//&
             " set incorrectly in namelist file '" // TRIM(namelist_file) // "'."
          CALL abort('setup_parameters', message) 
       END IF
!
!------------------------------------------------------------------------------
! Section 4: Compute additional parameters
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Section 4.1: COSMO-DE parameters
!------------------------------------------------------------------------------


 CALL run_control('time', 'init')
       ! Read COSMO-DE soil type map
       cosmo_var % name = 'SOILTYP'
       soiltyp = NINT(get_netcdf_variable(soiltyp_file, cosmo_var), hp)

       message = 'Reading PALM-4U origin from'
       IF (TRIM(static_driver_file) .NE. '')  THEN

          origin_lon = get_netcdf_attribute(static_driver_file, 'origin_lon')
          origin_lat = get_netcdf_attribute(static_driver_file, 'origin_lat')

          message = TRIM(message) // " static driver file '"                   &
                                  // TRIM(static_driver_file) // "'"


       ELSE

          origin_lon = longitude
          origin_lat = latitude

          message = TRIM(message) // " namlist file '"                         &
                                  // TRIM(namelist_file) // "'"

       END IF
       origin_lon = origin_lon * TO_RADIANS
       origin_lat = origin_lat * TO_RADIANS

       CALL report('setup_parameters', message)


       ! Read in COSMO-DE heights of half layers (vertical cell faces)
       cosmo_var % name = 'HHL'
       hhl = get_netcdf_variable(hhl_file, cosmo_var)
 CALL run_control('time', 'read')

       CALL reverse(hhl)
       nlon = SIZE(hhl, 1)
       nlat = SIZE(hhl, 2)
       nlev = SIZE(hhl, 3)

 CALL run_control('time', 'comp')

       ! Appoximate COSMO-DE heights of full layers (cell centres)
       ALLOCATE( hfl(nlon, nlat, nlev-1) )
 CALL run_control('time', 'alloc')
       DO k = 1, nlev-1
          hfl(:,:,k) = 0.5_dp * ( hhl(:,:,k) +                                 &
                                  hhl(:,:,k+1) )
       END DO

!------------------------------------------------------------------------------
! Section 4.2: PALM-4U parameters
!------------------------------------------------------------------------------
       ! PALM-4U domain extents
       lx = (nx+1) * dx
       ly = (ny+1) * dy
       lz = (nz+1) * dz
       
       ! PALM-4U point of Earth tangency
       x0 = 0.0_dp
       y0 = 0.0_dp

       ! time vector
       nt = CEILING(end_time / FORCING_FREQ) + 1
       ALLOCATE( time(nt) )
       CALL linspace(0.0_dp, 3600.0_dp * (nt-1), time)
       output_file % time => time
 CALL run_control('time', 'init')

! Convert the PALM-4U origin coordinates to COSMO's rotated-pole grid
       phi_c    = TO_RADIANS *                                                 &
                  phi2phirot( origin_lat * TO_DEGREES, origin_lon * TO_DEGREES,&
                              phi_n * TO_DEGREES, lambda_n * TO_DEGREES )
       lambda_c = TO_RADIANS *                                                 &
                  rla2rlarot( origin_lat * TO_DEGREES, origin_lon * TO_DEGREES,&
                              phi_n * TO_DEGREES, lambda_n * TO_DEGREES,     &
                              0.0_dp )

! Set gamma according to whether PALM domain is in the northern or southern
! hemisphere of the COSMO-DE rotated-pole system. Gamma assumes either the
! value 0 or PI and is needed to work around around a bug in the rotated-pole
! coordinate transformations.
       gam = gamma_from_hemisphere(origin_lat, phi_equat)

! Compute the north pole of the rotated-pole grid centred at the PALM-4U domain
! centre. The resulting (phi_cn, lambda_cn) are coordinates in COSMO-DE's
! rotated-pole grid.
       phi_cn    = phic_to_phin(phi_c) 
       lambda_cn = lamc_to_lamn(phi_c, lambda_c) 

       message =   "PALM-4U origin:" // NEW_LINE('') // &
          "           lon (lambda) = " // &
          TRIM(real_to_str_f(origin_lon * TO_DEGREES)) // " deg"// NEW_LINE(' ') //&
          "           lat (phi   ) = " // &
          TRIM(real_to_str_f(origin_lat * TO_DEGREES)) // " deg (geographical)" // NEW_LINE(' ') //&
          "           lon (lambda) = " // &
          TRIM(real_to_str_f(lambda_c * TO_DEGREES)) // " deg" // NEW_LINE(' ') // &
          "           lat (phi   ) = " // &
          TRIM(real_to_str_f(phi_c * TO_DEGREES)) // " deg (COSMO-DE rotated-pole)"
      CALL report ('setup_parameters', message)

       message = "North pole of the rotated COSMO-DE system:" // NEW_LINE(' ') // &
          "           lon (lambda) = " // &
          TRIM(real_to_str_f(lambda_n * TO_DEGREES)) // " deg" // NEW_LINE(' ') //&
          "           lat (phi   ) = " // &
          TRIM(real_to_str_f(phi_n * TO_DEGREES)) // " deg (geographical)"
       CALL report ('setup_parameters', message)
          
       message = "North pole of the rotated palm system:" // NEW_LINE(' ') // &
          "           lon (lambda) = " // &
          TRIM(real_to_str_f(lambda_cn * TO_DEGREES)) // " deg" // NEW_LINE(' ') // &
          "           lat (phi   ) = " // &
          TRIM(real_to_str_f(phi_cn * TO_DEGREES)) // " deg (COSMO-DE rotated-pole)"
       CALL report ('setup_parameters', message)

 CALL run_control('time', 'comp')

    END SUBROUTINE setup_parameters


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializes the COSMO-DE, PALM-4U, PALM-4U boundary grids.
!------------------------------------------------------------------------------!
    SUBROUTINE setup_grids() ! setup_grids(inifor_settings(with nx, ny, nz,...))
       CHARACTER ::  interp_mode
        
!------------------------------------------------------------------------------
! Section 1: Define model and initialization grids
!------------------------------------------------------------------------------
       CALL init_grid_definition('palm', grid=palm_grid,                       &
               xmin=0.0_dp, xmax=lx,                                           &
               ymin=0.0_dp, ymax=ly,                                           &
               zmin=0.0_dp, zmax=lz, x0=x0, y0=y0, z0=z0,                      &
               nx=nx, ny=ny, nz=nz, mode=mode)

       ! Subtracting 1 because arrays will be allocated with nlon + 1 elements.
       CALL init_grid_definition('cosmo-de', grid=cosmo_grid,                  &
               xmin=lonmin, xmax=lonmax,                                       &
               ymin=latmin, ymax=latmax,                                       &
               zmin=0.0_dp, zmax=51.0_dp, x0=x0, y0=y0, z0=0.0_dp,             &
               nx=nlon-1, ny=nlat-1, nz=nlev-1)

! Define intermediate grid. This is the same as palm_grid except with a much
! coarser vertical grid. The vertical levels are interpolated in each PALM-4U
! column from COSMO-DE's secondary levels. The main levels are then computed as
! the averages of the bounding secondary levels.
       CALL init_grid_definition('palm intermediate', grid=palm_intermediate,  &
               xmin=0.0_dp, xmax=lx,                                           &
               ymin=0.0_dp, ymax=ly,                                           &
               zmin=0.0_dp, zmax=lz, x0=x0, y0=y0, z0=z0,                      &
               nx=nx, ny=ny, nz=nlev-2)

       CALL init_grid_definition('boundary', grid=u_initial_grid,              &
               xmin = dx, xmax = lx - dx,                                      &
               ymin = 0.5_dp * dy, ymax = ly - 0.5_dp * dy,                    &
               zmin = 0.5_dp * dz, zmax = lz - 0.5_dp * dz,                    &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx-1, ny = ny, nz = nz,                                    &
               dx = dx, dy = dy, dz = dz, mode=mode)

       CALL init_grid_definition('boundary', grid=v_initial_grid,              &
               xmin = 0.5_dp * dx, xmax = lx - 0.5_dp * dx,                    &
               ymin = dy, ymax = ly - dy,                                      &
               zmin = 0.5_dp * dz, zmax = lz - 0.5_dp * dz,                    &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = ny-1, nz = nz,                                    &
               dx = dx, dy = dy, dz = dz, mode=mode)

       CALL init_grid_definition('boundary', grid=w_initial_grid,              &
               xmin = 0.5_dp * dx, xmax = lx - 0.5_dp * dx,                    &
               ymin = 0.5_dp * dy, ymax = ly - 0.5_dp * dy,                    &
               zmin = dz, zmax = lz - dz,                                      &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = ny, nz = nz-1,                                    &
               dx = dx, dy = dy, dz = dz, mode=mode)

       CALL init_grid_definition('boundary intermediate', grid=u_initial_intermediate,      &
               xmin = dx, xmax = lx - dx,                                      &
               ymin = 0.5_dp * dy, ymax = ly - 0.5_dp * dy,                    &
               zmin = 0.5_dp * dz, zmax = lz - 0.5_dp * dz,                    &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx-1, ny = ny, nz = nlev - 2,                              &
               dx = dx, dy = dy, dz = dz)

       CALL init_grid_definition('boundary intermediate', grid=v_initial_intermediate,      &
               xmin = 0.5_dp * dx, xmax = lx - 0.5_dp * dx,                    &
               ymin = dy, ymax = ly - dy,                                      &
               zmin = 0.5_dp * dz, zmax = lz - 0.5_dp * dz,                    &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = ny-1, nz = nlev - 2,                              &
               dx = dx, dy = dy, dz = dz)

       CALL init_grid_definition('boundary intermediate', grid=w_initial_intermediate,      &
               xmin = 0.5_dp * dx, xmax = lx - 0.5_dp * dx,                    &
               ymin = 0.5_dp * dy, ymax = ly - 0.5_dp * dy,                    &
               zmin = 0.5_dp * dz, zmax = lz - 0.5_dp * dz,                    &
               x0=x0, y0=y0, z0 = z0,                                          &
               nx = nx, ny = ny, nz = nlev - 2,                                &
               dx = dx, dy = dy, dz = dz)

      IF (boundary_variables_required)  THEN
!
!------------------------------------------------------------------------------
! Section 2: Define PALM-4U boundary grids
!------------------------------------------------------------------------------
          CALL init_grid_definition('boundary', grid=scalars_east_grid,           &
                  xmin = lx + 0.5_dp * dx, xmax = lx + 0.5_dp * dx,               &
                  ymin =  0.5_dp * dy, ymax = ly - 0.5_dp * dy,                   &
                  zmin =  0.5_dp * dz, zmax = lz - 0.5_dp * dz,                   &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = 0, ny = ny, nz = nz,                                       &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=scalars_west_grid,           &
                  xmin = -0.5_dp * dx, xmax = -0.5_dp * dx,                       &
                  ymin =  0.5_dp * dy, ymax = ly - 0.5_dp * dy,                   &
                  zmin =  0.5_dp * dz, zmax = lz - 0.5_dp * dz,                   &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = 0, ny = ny, nz = nz,                                       &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=scalars_north_grid,          &
                  xmin = 0.5_dp * dx, xmax = lx - 0.5_dp * dx,                    &
                  ymin = ly + 0.5_dp * dy, ymax = ly + 0.5_dp * dy,               &
                  zmin = 0.5_dp * dz, zmax = lz - 0.5_dp * dz,                    &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx, ny = 0, nz = nz,                                       &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=scalars_south_grid,          &
                  xmin =  0.5_dp * dx, xmax = lx - 0.5_dp * dx,                   &
                  ymin = -0.5_dp * dy, ymax = -0.5_dp * dy,                       &
                  zmin =  0.5_dp * dz, zmax = lz - 0.5_dp * dz,                   &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx, ny = 0, nz = nz,                                       &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=scalars_top_grid,            &
                  xmin =  0.5_dp * dx, xmax = lx - 0.5_dp * dx,                   &
                  ymin =  0.5_dp * dy, ymax = ly - 0.5_dp * dy,                   &
                  zmin =  lz + 0.5_dp * dz, zmax = lz + 0.5_dp * dz,              &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx, ny = ny, nz = 0,                                       &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=u_east_grid,                 &
                  xmin = lx, xmax = lx,                                           &
                  ymin = 0.5_dp * dy, ymax = ly - 0.5_dp * dy,                    &
                  zmin = 0.5_dp * dz, zmax = lz - 0.5_dp * dz,                    &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = 0, ny = ny, nz = nz,                                       &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=u_west_grid,                 &
                  xmin = 0.0_dp, xmax = 0.0_dp,                                   &
                  ymin = 0.5_dp * dy, ymax = ly - 0.5_dp * dy,                    &
                  zmin = 0.5_dp * dz, zmax = lz - 0.5_dp * dz,                    &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = 0, ny = ny, nz = nz,                                       &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=u_north_grid,                &
                  xmin = dx, xmax = lx - dx,                                      &
                  ymin = ly + 0.5_dp * dy, ymax = ly + 0.5_dp * dy,               &
                  zmin = 0.5_dp * dz, zmax = lz - 0.5_dp * dz,                    &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx-1, ny = 0, nz = nz,                                     &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=u_south_grid,                &
                  xmin = dx, xmax = lx - dx,                                      &
                  ymin = -0.5_dp * dy, ymax = -0.5_dp * dy,                       &
                  zmin = 0.5_dp * dz, zmax = lz - 0.5_dp * dz,                    &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx-1, ny = 0, nz = nz,                                     &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=u_top_grid,                  &
                  xmin = dx, xmax = lx - dx,                                      &
                  ymin = 0.5_dp * dy, ymax = ly - 0.5_dp * dy,                    &
                  zmin = lz + 0.5_dp * dz, zmax = lz + 0.5_dp * dz,               &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx-1, ny = ny, nz = 0,                                     &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=v_east_grid,                 &
                  xmin = lx + 0.5_dp * dx, xmax = lx + 0.5_dp * dx,               &
                  ymin = dy, ymax = ly - dy,                                      &
                  zmin = 0.5_dp * dz, zmax = lz - 0.5_dp * dz,                    &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = 0, ny = ny-1, nz = nz,                                     &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=v_west_grid,                 &
                  xmin = -0.5_dp * dx, xmax = -0.5_dp * dx,                       &
                  ymin = dy, ymax = ly - dy,                                      &
                  zmin = 0.5_dp * dz, zmax = lz - 0.5_dp * dz,                    &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = 0, ny = ny-1, nz = nz,                                     &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=v_north_grid,                &
                  xmin = 0.5_dp * dx, xmax = lx - 0.5_dp * dx,                    &
                  ymin = ly, ymax = ly,                                           &
                  zmin = 0.5_dp * dz, zmax = lz - 0.5_dp * dz,                    &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx, ny = 0, nz = nz,                                       &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=v_south_grid,                &
                  xmin = 0.5_dp * dx, xmax = lx - 0.5_dp * dx,                    &
                  ymin = 0.0_dp, ymax = 0.0_dp,                                   &
                  zmin = 0.5_dp * dz, zmax = lz - 0.5_dp * dz,                    &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx, ny = 0, nz = nz,                                       &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=v_top_grid,                  &
                  xmin = 0.5_dp * dx, xmax = lx - 0.5_dp * dx,                    &
                  ymin = dy, ymax = ly - dy,                                      &
                  zmin = lz + 0.5_dp * dz, zmax = lz + 0.5_dp * dz,               &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx, ny = ny-1, nz = 0,                                     &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=w_east_grid,                 &
                  xmin = lx + 0.5_dp * dx, xmax = lx + 0.5_dp * dx,               &
                  ymin =  0.5_dp * dy, ymax = ly - 0.5_dp * dy,                   &
                  zmin =  dz, zmax = lz - dz,                                     &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = 0, ny = ny, nz = nz - 1,                                   &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=w_west_grid,                 &
                  xmin = -0.5_dp * dx, xmax = -0.5_dp * dx,                       &
                  ymin = 0.5_dp * dy, ymax = ly - 0.5_dp * dy,                    &
                  zmin = dz, zmax = lz - dz,                                      &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = 0, ny = ny, nz = nz - 1,                                   &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=w_north_grid,                &
                  xmin = 0.5_dp * dx, xmax = lx - 0.5_dp * dx,                    &
                  ymin = ly + 0.5_dp * dy, ymax = ly + 0.5_dp * dy,               &
                  zmin = dz, zmax = lz - dz,                                      &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx, ny = 0, nz = nz - 1,                                   &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=w_south_grid,                &
                  xmin =  0.5_dp * dx, xmax = lx - 0.5_dp * dx,                   &
                  ymin = -0.5_dp * dy, ymax = -0.5_dp * dy,                       &
                  zmin = dz, zmax = lz - dz,                                      &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx, ny = 0, nz = nz - 1,                                   &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=w_top_grid,                  &
                  xmin = 0.5_dp * dx, xmax = lx - 0.5_dp * dx,                    &
                  ymin = 0.5_dp * dy, ymax = ly - 0.5_dp * dy,                    &
                  zmin = lz, zmax = lz,                                           &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx, ny = ny, nz = 0,                                       &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary intermediate', grid=scalars_east_intermediate,   &
                  xmin = lx + 0.5_dp * dx, xmax = lx + 0.5_dp * dx,               &
                  ymin =  0.5_dp * dy, ymax = ly - 0.5_dp * dy,                   &
                  zmin =  0.0_dp, zmax = 0.0_dp,                                  &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = 0, ny = ny, nz = nlev - 2,                                 &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary intermediate', grid=scalars_west_intermediate,   &
                  xmin = -0.5_dp * dx, xmax = -0.5_dp * dx,                       &
                  ymin =  0.5_dp * dy, ymax = ly - 0.5_dp * dy,                   &
                  zmin =  0.0_dp, zmax = 0.0_dp,                                  &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = 0, ny = ny, nz = nlev - 2,                                 &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary intermediate', grid=scalars_north_intermediate,  &
                  xmin = 0.5_dp * dx, xmax = lx - 0.5_dp * dx,                    &
                  ymin = ly + 0.5_dp * dy, ymax = ly + 0.5_dp * dy,               &
                  zmin =  0.0_dp, zmax = 0.0_dp,                                  &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx, ny = 0, nz = nlev - 2,                                 &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary intermediate', grid=scalars_south_intermediate,  &
                  xmin =  0.5_dp * dx, xmax = lx - 0.5_dp * dx,                   &
                  ymin = -0.5_dp * dy, ymax = -0.5_dp * dy,                       &
                  zmin =  0.0_dp, zmax = 0.0_dp,                                  &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx, ny = 0, nz = nlev - 2,                                 &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary intermediate', grid=scalars_top_intermediate,    &
                  xmin =  0.5_dp * dx, xmax = lx - 0.5_dp * dx,                   &
                  ymin =  0.5_dp * dy, ymax = ly - 0.5_dp * dy,                   &
                  zmin =  0.0_dp, zmax = 0.0_dp,                                  &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx, ny = ny, nz = nlev - 2,                                &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary intermediate', grid=u_east_intermediate,         &
                  xmin = lx, xmax = lx,                                           &
                  ymin = 0.5_dp * dy, ymax = ly - 0.5_dp * dy,                    &
                  zmin =  0.0_dp, zmax = 0.0_dp,                                  &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = 0, ny = ny, nz = nlev - 2,                                 &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary intermediate', grid=u_west_intermediate,         &
                  xmin = 0.0_dp, xmax = 0.0_dp,                                   &
                  ymin = 0.5_dp * dy, ymax = ly - 0.5_dp * dy,                    &
                  zmin =  0.0_dp, zmax = 0.0_dp,                                  &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = 0, ny = ny, nz = nlev - 2,                                 &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary intermediate', grid=u_north_intermediate,        &
                  xmin = dx, xmax = lx - dx,                                      &
                  ymin = ly + 0.5_dp * dy, ymax = ly + 0.5_dp * dy,               &
                  zmin =  0.0_dp, zmax = 0.0_dp,                                  &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx-1, ny = 0, nz = nlev - 2,                               &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary intermediate', grid=u_south_intermediate,        &
                  xmin = dx, xmax = lx - dx,                                      &
                  ymin = -0.5_dp * dy, ymax = -0.5_dp * dy,                       &
                  zmin =  0.0_dp, zmax = 0.0_dp,                                  &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx-1, ny = 0, nz = nlev - 2,                               &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary intermediate', grid=u_top_intermediate,          &
                  xmin = dx, xmax = lx - dx,                                      &
                  ymin = 0.5_dp * dy, ymax = ly - 0.5_dp * dy,                    &
                  zmin =  0.0_dp, zmax = 0.0_dp,                                  &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx-1, ny = ny, nz = nlev - 2,                              &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary intermediate', grid=v_east_intermediate,         &
                  xmin = lx + 0.5_dp * dx, xmax = lx + 0.5_dp * dx,               &
                  ymin = dy, ymax = ly - dy,                                      &
                  zmin =  0.0_dp, zmax = 0.0_dp,                                  &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = 0, ny = ny-1, nz = nlev - 2,                               &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary intermediate', grid=v_west_intermediate,         &
                  xmin = -0.5_dp * dx, xmax = -0.5_dp * dx,                       &
                  ymin = dy, ymax = ly - dy,                                      &
                  zmin =  0.0_dp, zmax = 0.0_dp,                                  &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = 0, ny = ny-1, nz = nlev - 2,                               &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary intermediate', grid=v_north_intermediate,        &
                  xmin = 0.5_dp * dx, xmax = lx - 0.5_dp * dx,                    &
                  ymin = ly, ymax = ly,                                           &
                  zmin =  0.0_dp, zmax = 0.0_dp,                                  &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx, ny = 0, nz = nlev - 2,                                 &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary intermediate', grid=v_south_intermediate,        &
                  xmin = 0.5_dp * dx, xmax = lx - 0.5_dp * dx,                    &
                  ymin = 0.0_dp, ymax = 0.0_dp,                                   &
                  zmin =  0.0_dp, zmax = 0.0_dp,                                  &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx, ny = 0, nz = nlev - 2,                                 &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary intermediate', grid=v_top_intermediate,          &
                  xmin = 0.5_dp * dx, xmax = lx - 0.5_dp * dx,                    &
                  ymin = dy, ymax = ly - dy,                                      &
                  zmin =  0.0_dp, zmax = 0.0_dp,                                  &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx, ny = ny-1, nz = nlev - 2,                              &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary intermediate', grid=w_east_intermediate,         &
                  xmin = lx + 0.5_dp * dx, xmax = lx + 0.5_dp * dx,               &
                  ymin =  0.5_dp * dy, ymax = ly - 0.5_dp * dy,                   &
                  zmin =  0.0_dp, zmax = 0.0_dp,                                  &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = 0, ny = ny, nz = nlev - 2,                                 &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary intermediate', grid=w_west_intermediate,         &
                  xmin = -0.5_dp * dx, xmax = -0.5_dp * dx,                       &
                  ymin =  0.5_dp * dy, ymax = ly - 0.5_dp * dy,                   &
                  zmin =  0.0_dp, zmax = 0.0_dp,                                  &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = 0, ny = ny, nz = nlev - 2,                                 &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary intermediate', grid=w_north_intermediate,        &
                  xmin = 0.5_dp * dx, xmax = lx - 0.5_dp * dx,                    &
                  ymin = ly + 0.5_dp * dy, ymax = ly + 0.5_dp * dy,               &
                  zmin =  0.0_dp, zmax = 0.0_dp,                                  &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx, ny = 0, nz = nlev - 2,                                 &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary intermediate', grid=w_south_intermediate,        &
                  xmin =  0.5_dp * dx, xmax = lx - 0.5_dp * dx,                   &
                  ymin = -0.5_dp * dy, ymax = -0.5_dp * dy,                       &
                  zmin =  0.0_dp, zmax = 0.0_dp,                                  &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx, ny = 0, nz = nlev - 2,                                 &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary intermediate', grid=w_top_intermediate,          &
                  xmin =  0.5_dp * dx, xmax = lx - 0.5_dp * dx,                   &
                  ymin =  0.5_dp * dy, ymax = ly - 0.5_dp * dy,                   &
                  zmin =  0.0_dp, zmax = 0.0_dp,                                  &
                  x0=x0, y0=y0, z0 = z0,                                          &
                  nx = nx, ny = ny, nz = nlev - 2,                                &
                  dx = dx, dy = dy, dz = dz)
       END IF

!                                                                              
!------------------------------------------------------------------------------
! Section 3: Define profile grids
!------------------------------------------------------------------------------

       IF (TRIM(mode) == 'profile')  THEN
          CALL init_grid_definition('boundary', grid=scalar_profile_grid,      &
                  xmin = 0.5_dp * lx, xmax = 0.5_dp * lx,                      &
                  ymin = 0.5_dp * ly, ymax = 0.5_dp * ly,                      &
                  zmin = 0.5_dp * dz, zmax = lz - 0.5_dp * dz,                 &
                  x0=x0, y0=y0, z0 = z0,                                       &
                  nx = 0, ny = 0, nz = nz,                                     &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=w_profile_grid,           &
                  xmin = 0.5_dp * lx, xmax = 0.5_dp * lx,                      &
                  ymin = 0.5_dp * ly, ymax = 0.5_dp * ly,                      &
                  zmin = dz, zmax = lz - dz,                                   &
                  x0=x0, y0=y0, z0 = z0,                                       &
                  nx = 0, ny = 0, nz = nz - 1,                                 &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=scalar_profile_intermediate,&
                  xmin = 0.5_dp * lx, xmax = 0.5_dp * lx,                      &
                  ymin = 0.5_dp * ly, ymax = 0.5_dp * ly,                      &
                  zmin = 0.0_dp, zmax = 0.0_dp,                                &
                  x0=x0, y0=y0, z0 = z0,                                       &
                  nx = 0, ny = 0, nz = nlev - 2,                               &
                  dx = dx, dy = dy, dz = dz)

          CALL init_grid_definition('boundary', grid=w_profile_intermediate,   &
                  xmin = 0.5_dp * lx, xmax = 0.5_dp * lx,                      &
                  ymin = 0.5_dp * ly, ymax = 0.5_dp * ly,                      &
                  zmin = 0.0_dp, zmax = 0.0_dp,                                &
                  x0=x0, y0=y0, z0 = z0,                                       &
                  nx = 0, ny = 0, nz = nlev - 2,                               &
                  dx = dx, dy = dy, dz = dz)
       END IF

!                                                                              
!------------------------------------------------------------------------------
! Section 4: Precompute neighbours and weights for interpolation              
!------------------------------------------------------------------------------
       interp_mode = 's'
       CALL setup_interpolation(cosmo_grid, palm_grid, palm_intermediate, interp_mode, mode=mode)
       IF (boundary_variables_required)  THEN
          CALL setup_interpolation(cosmo_grid, scalars_east_grid, scalars_east_intermediate, interp_mode)
          CALL setup_interpolation(cosmo_grid, scalars_west_grid, scalars_west_intermediate, interp_mode)
          CALL setup_interpolation(cosmo_grid, scalars_north_grid, scalars_north_intermediate, interp_mode)
          CALL setup_interpolation(cosmo_grid, scalars_south_grid, scalars_south_intermediate, interp_mode)
          CALL setup_interpolation(cosmo_grid, scalars_top_grid, scalars_top_intermediate, interp_mode)
       END IF

       interp_mode = 'u'
       CALL setup_interpolation(cosmo_grid, u_initial_grid, u_initial_intermediate, interp_mode, mode=mode)
       IF (boundary_variables_required)  THEN
          CALL setup_interpolation(cosmo_grid, u_east_grid, u_east_intermediate, interp_mode)
          CALL setup_interpolation(cosmo_grid, u_west_grid, u_west_intermediate, interp_mode)
          CALL setup_interpolation(cosmo_grid, u_north_grid, u_north_intermediate, interp_mode)
          CALL setup_interpolation(cosmo_grid, u_south_grid, u_south_intermediate, interp_mode)
          CALL setup_interpolation(cosmo_grid, u_top_grid, u_top_intermediate, interp_mode)
       END IF

       interp_mode = 'v'
       CALL setup_interpolation(cosmo_grid, v_initial_grid, v_initial_intermediate, interp_mode, mode=mode)
       IF (boundary_variables_required)  THEN
          CALL setup_interpolation(cosmo_grid, v_east_grid, v_east_intermediate, interp_mode)
          CALL setup_interpolation(cosmo_grid, v_west_grid, v_west_intermediate, interp_mode)
          CALL setup_interpolation(cosmo_grid, v_north_grid, v_north_intermediate, interp_mode)
          CALL setup_interpolation(cosmo_grid, v_south_grid, v_south_intermediate, interp_mode)
          CALL setup_interpolation(cosmo_grid, v_top_grid, v_top_intermediate, interp_mode)
       END IF

       interp_mode = 'w'
       CALL setup_interpolation(cosmo_grid, w_initial_grid, w_initial_intermediate, interp_mode, mode=mode)
       IF (boundary_variables_required)  THEN
          CALL setup_interpolation(cosmo_grid, w_east_grid, w_east_intermediate, interp_mode)
          CALL setup_interpolation(cosmo_grid, w_west_grid, w_west_intermediate, interp_mode)
          CALL setup_interpolation(cosmo_grid, w_north_grid, w_north_intermediate, interp_mode)
          CALL setup_interpolation(cosmo_grid, w_south_grid, w_south_intermediate, interp_mode)
          CALL setup_interpolation(cosmo_grid, w_top_grid, w_top_intermediate, interp_mode)
       END IF

       IF (TRIM(mode) == 'profile')  THEN
           CALL setup_averaging(cosmo_grid, palm_intermediate, imin, imax, jmin, jmax) 
       END IF
        

    END SUBROUTINE setup_grids


    SUBROUTINE setup_interpolation(cosmo_grid, palm_grid, palm_intermediate, kind, mode)

       TYPE(grid_definition), INTENT(IN), TARGET    ::  cosmo_grid
       TYPE(grid_definition), INTENT(INOUT), TARGET ::  palm_grid, palm_intermediate
       CHARACTER, INTENT(IN)                        ::  kind
       CHARACTER(LEN=*), INTENT(IN), OPTIONAL       ::  mode

       TYPE(grid_definition), POINTER      ::  grid
       REAL(dp), DIMENSION(:), POINTER     ::  lat, lon
       REAL(dp), DIMENSION(:,:,:), POINTER ::  h

       LOGICAL :: setup_vertical = .TRUE.

       IF (PRESENT(mode))  THEN
          IF (TRIM(mode) == 'profile')  setup_vertical = .FALSE.
       ELSE
          setup_vertical = .TRUE.
       END IF

!------------------------------------------------------------------------------
! Section 1: Horizontal interpolation                                        
!------------------------------------------------------------------------------
       ! Select horizontal coordinates according to kind of points (s/w, u, v)
       SELECT CASE(kind)

       CASE('s') ! scalars

          lat => cosmo_grid % lat
          lon => cosmo_grid % lon
          h   => cosmo_grid % hfl

       CASE('w') ! vertical velocity

          lat => cosmo_grid % lat
          lon => cosmo_grid % lon
          h   => cosmo_grid % hhl

       CASE('u') ! x velocity

          lat => cosmo_grid % lat
          lon => cosmo_grid % lonu
          h   => cosmo_grid % hfl

       CASE('v') ! y velocity

          lat => cosmo_grid % latv
          lon => cosmo_grid % lon
          h   => cosmo_grid % hfl

       CASE DEFAULT

          message = "Interpolation mode '" // mode // "' is not supported."
          CALL abort('setup_interpolation', message)

       END SELECT

       grid => palm_intermediate

       CALL find_horizontal_neighbours(lat, lon,                               &
          cosmo_grid % dxi, cosmo_grid % dyi, grid % clat,                     &
          grid % clon, grid % ii, grid % jj)

       CALL compute_horizontal_interp_weights(lat, lon,                        &
          cosmo_grid % dxi, cosmo_grid % dyi, grid % clat,                     &
          grid % clon, grid % ii, grid % jj, grid % w_horiz)

!------------------------------------------------------------------------------
! Section 2: Vertical interpolation
!------------------------------------------------------------------------------

       IF (setup_vertical)  THEN
          ALLOCATE( grid % h(0:grid % nx, 0:grid % ny, 0:grid % nz) ) 
          grid % h(:,:,:) = - EARTH_RADIUS

          ! For w points, use hhl, for scalars use hfl
          ! compute the full heights for the intermediate grids
          CALL interpolate_2d(cosmo_grid % hfl, grid % h, grid)
          CALL find_vertical_neighbours_and_weights(palm_grid, grid)
       END IF
       
    END SUBROUTINE setup_interpolation


    SUBROUTINE setup_averaging(cosmo_grid, palm_intermediate, imin, imax, jmin, jmax)

       TYPE(grid_definition), INTENT(IN) ::  cosmo_grid, palm_intermediate
       INTEGER, INTENT(INOUT)            ::  imin, imax, jmin, jmax

       TYPE(grid_definition), POINTER    ::  grid
       REAL ::  lonmin_pos,lonmax_pos, latmin_pos, latmax_pos

       ! find horizontal index ranges for profile averaging
       lonmin_pos = (MINVAL(palm_intermediate % clon(:,:)) - cosmo_grid % lon(0)) * cosmo_grid % dxi
       lonmax_pos = (MAXVAL(palm_intermediate % clon(:,:)) - cosmo_grid % lon(0)) * cosmo_grid % dxi
       latmin_pos = (MINVAL(palm_intermediate % clat(:,:)) - cosmo_grid % lat(0)) * cosmo_grid % dyi
       latmax_pos = (MAXVAL(palm_intermediate % clat(:,:)) - cosmo_grid % lat(0)) * cosmo_grid % dyi

       imin = FLOOR(lonmin_pos)
       imax = CEILING(lonmax_pos)
       jmin = FLOOR(latmin_pos)
       jmax = CEILING(latmax_pos)
       
       ! average heights for intermediate scalar and w profile grids
       grid => scalar_profile_intermediate
       ALLOCATE( grid % h(0:grid % nx, 0:grid % ny, 0:grid % nz) ) 
       grid % h(:,:,:) = - EARTH_RADIUS
       CALL average_2d(cosmo_grid % hfl, grid % h(0,0,:), imin, imax, jmin, jmax)
       CALL find_vertical_neighbours_and_weights(scalar_profile_grid, grid)

       grid => w_profile_intermediate
       ALLOCATE( grid % h(0:grid % nx, 0:grid % ny, 0:grid % nz) ) 
       grid % h(:,:,:) = - EARTH_RADIUS
       CALL average_2d(cosmo_grid % hhl, grid % h(0,0,:), imin, imax, jmin, jmax)
       CALL find_vertical_neighbours_and_weights(w_profile_grid, grid)

    END SUBROUTINE setup_averaging


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Helper function that computes horizontal domain extend in x or y direction
!> such that the centres of a boundary grid fall at -dx/2 or lx + dx/2.
!> 
!> Input parameters:
!> -----------------
!> dxy : grid spacing in x or y direction
!> lxy : domain length in dxy direction
!> 
!> Output parameters:
!> ------------------
!> boundary_extent : Domain minimum xymin (maximum xymax) if dxy < 0 (> 0)
!------------------------------------------------------------------------------!
    REAL(dp) FUNCTION boundary_extent(dxy, lxy)
        REAL(dp), INTENT(IN) ::  dxy, lxy

        boundary_extent = 0.5_dp * lxy + SIGN(lxy + ABS(dxy), dxy)

    END FUNCTION boundary_extent


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializes grid_definition-type variables.
!> 
!> Input parameters:
!> -----------------
!> mode : Initialization mode, distinguishes between PALM-4U and COSMO-DE grids
!>    as well as grids covering the boundary surfaces. Valid modes are:
!>       - 'palm'
!>       - 'cosmo-de'
!>       - 'eastwest-scalar'
!> 
!> <xyx>min, <xyz>max : Domain minima and maxima in x, y, and z direction. Note
!>    that these values do not necessarily translate to the outmost coordinates
!>    of the generated grid but rather refer to the extent of the underlying
!>    PALM-4U computational domain (i.e. the outer cell faces). The coordinates
!>    of the generated grid will be inferred from this information taking into
!>    account the initialization mode. For example, the coordinates of a
!>    boundary grid initialized using mode 'eastwest-scalar' will be located in
!>    planes one half grid point outwards of xmin and xmax.
!>
!> z0 : Elevation of the PALM-4U domain above sea level [m]
!>
!> n<xyz> : Number of grod points in x, y, and z direction
!>
!> Output parameters:
!> ------------------
!> grid : Grid variable to be initialized.
!------------------------------------------------------------------------------!
    SUBROUTINE init_grid_definition(kind, xmin, xmax, ymin, ymax, zmin, zmax,  &
                                    x0, y0, z0, nx, ny, nz, dx, dy, dz, grid, mode)
        CHARACTER(LEN=*), INTENT(IN)           ::  kind
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL ::  mode
        INTEGER, INTENT(IN)                    ::  nx, ny, nz
        REAL(dp), INTENT(IN)                   ::  xmin, xmax, ymin, ymax, zmin, zmax
        REAL(dp), INTENT(IN)                   ::  x0, y0, z0
        REAL(dp), OPTIONAL, INTENT(IN)         ::  dx, dy, dz
        TYPE(grid_definition), INTENT(INOUT)   ::  grid

        grid % nx = nx
        grid % ny = ny
        grid % nz = nz

        grid % lx = xmax - xmin
        grid % ly = ymax - ymin
        grid % lz = zmax - zmin

        grid % x0 = x0
        grid % y0 = y0
        grid % z0 = z0

        SELECT CASE( TRIM (kind) )

        CASE('boundary')
           IF (.NOT. PRESENT(dx))  THEN
              message = "dx is not present but needed for 'eastwest-scalar' "//&
                        "initializaton."
              CALL abort('init_grid_definition', message)
           END IF
           IF (.NOT. PRESENT(dy))  THEN
              message = "dy is not present but needed for 'eastwest-scalar' "//&
                        "initializaton."
              CALL abort('init_grid_definition', message)
           END IF
           IF (.NOT. PRESENT(dz))  THEN
              message = "dz is not present but needed for 'eastwest-scalar' "//&
                        "initializaton."
              CALL abort('init_grid_definition', message)
           END IF

           grid % dx  = dx
           grid % dy  = dy
           grid % dz  = dz

           grid % dxi = 1.0_dp / grid % dx
           grid % dyi = 1.0_dp / grid % dy
           grid % dzi = 1.0_dp / grid % dz

           ALLOCATE( grid % x(0:nx) )
           CALL linspace(xmin, xmax, grid % x)

           ALLOCATE( grid % y(0:ny) )
           CALL linspace(ymin, ymax, grid % y)

           ALLOCATE( grid % z(0:nz) )
           CALL linspace(zmin, zmax, grid % z)

           ! Allocate neighbour indices and weights
           IF (TRIM(mode) .NE. 'profile')  THEN
              ALLOCATE( grid % kk(0:nx, 0:ny, 0:nz, 2) )
              grid % kk(:,:,:,:) = -1

              ALLOCATE( grid % w_verti(0:nx, 0:ny, 0:nz, 2) )
              grid % w_verti(:,:,:,:) = 0.0_dp
           END IF
        
        CASE('boundary intermediate')
           IF (.NOT. PRESENT(dx))  THEN
              message = "dx is not present but needed for 'eastwest-scalar' "//&
                        "initializaton."
              CALL abort('init_grid_definition', message)
           END IF
           IF (.NOT. PRESENT(dy))  THEN
              message = "dy is not present but needed for 'eastwest-scalar' "//&
                        "initializaton."
              CALL abort('init_grid_definition', message)
           END IF
           IF (.NOT. PRESENT(dz))  THEN
              message = "dz is not present but needed for 'eastwest-scalar' "//&
                        "initializaton."
              CALL abort('init_grid_definition', message)
           END IF

           grid % dx  = dx
           grid % dy  = dy
           grid % dz  = dz

           grid % dxi = 1.0_dp / grid % dx
           grid % dyi = 1.0_dp / grid % dy
           grid % dzi = 1.0_dp / grid % dz

           ALLOCATE( grid % x(0:nx) )
           CALL linspace(xmin, xmax, grid % x)

           ALLOCATE( grid % y(0:ny) )
           CALL linspace(ymin, ymax, grid % y)

           ALLOCATE( grid % z(0:nz) )
           CALL linspace(zmin, zmax, grid % z)

           ALLOCATE( grid % clon(0:nx, 0:ny), grid % clat(0:nx, 0:ny)  )
           CALL rotate_to_cosmo(                                               &
              phir = project( grid % y, y0, EARTH_RADIUS ) , & ! = plate-carree latitude
              lamr = project( grid % x, x0, EARTH_RADIUS ) , & ! = plate-carree longitude
              phip = phi_cn, lamp = lambda_cn,                                 &
              phi  = grid % clat,                                              &
              lam  = grid % clon,                                              &
              gam  = gam                                                       &
           )

           ! Allocate neighbour indices and weights
           ALLOCATE( grid % ii(0:nx, 0:ny, 4),                                 &
                     grid % jj(0:nx, 0:ny, 4) )
           grid % ii(:,:,:)   = -1
           grid % jj(:,:,:)   = -1

           ALLOCATE( grid % w_horiz(0:nx, 0:ny, 4) )
           grid % w_horiz(:,:,:)   = 0.0_dp
        
        ! This mode initializes a Cartesian PALM-4U grid and adds the
        ! corresponding latitudes and longitudes of the rotated pole grid.
        CASE('palm')
           grid % name(1) = 'x and lon'
           grid % name(2) = 'y and lat'
           grid % name(3) = 'z'

           grid % dx  = grid % lx / (nx + 1)
           grid % dy  = grid % ly / (ny + 1)
           grid % dz  = grid % lz / (nz + 1)

           grid % dxi = 1.0_dp / grid % dx
           grid % dyi = 1.0_dp / grid % dy
           grid % dzi = 1.0_dp / grid % dz

           ALLOCATE( grid % x(0:nx),   grid % y(0:ny),  grid % z(0:nz)  )
           ALLOCATE( grid % xu(1:nx),  grid % yv(1:ny), grid % zw(1:nz) )
           CALL linspace(xmin + 0.5_dp*grid % dx, xmax - 0.5_dp*grid % dx, grid % x)
           CALL linspace(ymin + 0.5_dp*grid % dy, ymax - 0.5_dp*grid % dy, grid % y)
           CALL linspace(zmin + 0.5_dp*grid % dz, zmax - 0.5_dp*grid % dz, grid % z)
           CALL linspace(xmin + grid % dx, xmax - grid % dx, grid % xu)
           CALL linspace(ymin + grid % dy, ymax - grid % dy, grid % yv)
           CALL linspace(zmin + grid % dz, zmax - grid % dz, grid % zw)

           grid % depths => depths

           ! Allocate neighbour indices and weights
           IF (TRIM(mode) .NE. 'profile')  THEN
              ALLOCATE( grid % kk(0:nx, 0:ny, 0:nz, 2) )
              grid % kk(:,:,:,:) = -1

              ALLOCATE( grid % w_verti(0:nx, 0:ny, 0:nz, 2) )
              grid % w_verti(:,:,:,:) = 0.0_dp
           END IF

        CASE('palm intermediate')
           grid % name(1) = 'x and lon'
           grid % name(2) = 'y and lat'
           grid % name(3) = 'z'

           grid % dx  = grid % lx / (nx + 1)
           grid % dy  = grid % ly / (ny + 1)
           grid % dz  = grid % lz / (nz + 1)

           grid % dxi = 1.0_dp / grid % dx
           grid % dyi = 1.0_dp / grid % dy
           grid % dzi = 1.0_dp / grid % dz

           ALLOCATE( grid % x(0:nx),   grid % y(0:ny),  grid % z(0:nz)  )
           ALLOCATE( grid % xu(1:nx),  grid % yv(1:ny), grid % zw(1:nz) )
           CALL linspace(xmin + 0.5_dp*grid % dx, xmax - 0.5_dp*grid % dx, grid % x)
           CALL linspace(ymin + 0.5_dp*grid % dy, ymax - 0.5_dp*grid % dy, grid % y)
           CALL linspace(zmin + 0.5_dp*grid % dz, zmax - 0.5_dp*grid % dz, grid % z)
           CALL linspace(xmin + grid % dx, xmax - grid % dx, grid % xu)
           CALL linspace(ymin + grid % dy, ymax - grid % dy, grid % yv)
           CALL linspace(zmin + grid % dz, zmax - grid % dz, grid % zw)

           grid % depths => depths

           ! Allocate rotated-pole coordinates, clon is for (c)osmo-de (lon)gitude
           ALLOCATE( grid % clon(0:nx, 0:ny),   grid % clat(0:nx, 0:ny)  )
           ALLOCATE( grid % clonu(1:nx, 0:ny),  grid % clatu(1:nx, 0:ny) )
           ALLOCATE( grid % clonv(0:nx, 1:ny),  grid % clatv(0:nx, 1:ny) )

           ! Compute rotated-pole coordinates of...
           ! ... PALM-4U centres
           CALL rotate_to_cosmo(                                               &
              phir = project( grid % y, y0, EARTH_RADIUS ) , & ! = plate-carree latitude
              lamr = project( grid % x, x0, EARTH_RADIUS ) , & ! = plate-carree longitude
              phip = phi_cn, lamp = lambda_cn,                                 &
              phi  = grid % clat,                                              &
              lam  = grid % clon,                                              &
              gam  = gam                                                       &
           )

           ! ... PALM-4U u winds
           CALL rotate_to_cosmo(                                               &
              phir = project( grid % y,  y0, EARTH_RADIUS ), & ! = plate-carree latitude
              lamr = project( grid % xu, x0, EARTH_RADIUS ), & ! = plate-carree longitude
              phip = phi_cn, lamp = lambda_cn,                                 &
              phi  = grid % clatu,                                             &
              lam  = grid % clonu,                                             &
              gam  = gam                                                       &
           )

           ! ... PALM-4U v winds
           CALL rotate_to_cosmo(                                               &
              phir = project( grid % yv, y0, EARTH_RADIUS ), & ! = plate-carree latitude
              lamr = project( grid % x,  x0, EARTH_RADIUS ), & ! = plate-carree longitude
              phip = phi_cn, lamp = lambda_cn,                                 &
              phi  = grid % clatv,                                             &
              lam  = grid % clonv,                                             &
              gam  = gam                                                       &
           )

           ! Allocate neighbour indices and weights
           ALLOCATE( grid % ii(0:nx, 0:ny, 4),                                 &
                     grid % jj(0:nx, 0:ny, 4) )
           grid % ii(:,:,:)   = -1
           grid % jj(:,:,:)   = -1

           ALLOCATE( grid % w_horiz(0:nx, 0:ny, 4) )
           grid % w_horiz(:,:,:)   = 0.0_dp

        CASE('cosmo-de')
           grid % name(1) = 'rlon'         ! of COMSO-DE cell centres (scalars)
           grid % name(2) = 'rlat'         ! of COMSO-DE cell centres (scalars)
           grid % name(3) = 'height'

           grid % dx  = grid % lx / nx     ! = 0.025 deg, stored in radians
           grid % dy  = grid % ly / ny     ! = 0.025 deg, stored in radians
           grid % dz  = 0.0_dp             ! not defined yet

           grid % dxi = 1.0_dp / grid % dx ! [rad^-1]
           grid % dyi = 1.0_dp / grid % dy ! [rad^-1]
           grid % dzi = 0.0_dp             ! not defined yet

           ALLOCATE( grid % lon(0:nx),   grid % lat(0:ny),  grid % z(0:nz)  )
           ALLOCATE( grid % lonu(0:nx),  grid % latv(0:ny), grid % zw(0:nz) )

           CALL linspace(xmin, xmax, grid % lon)
           CALL linspace(ymin, ymax, grid % lat)
           grid % lonu(:) = grid % lon + 0.5_dp * grid % dx
           grid % latv(:) = grid % lat + 0.5_dp * grid % dy

           ! Point to heights of half levels (hhl) and compute heights of full
           ! levels (hfl) as arithmetic averages
           grid % hhl => hhl
           grid % hfl => hfl
           grid % depths => depths

        CASE DEFAULT
            message = "Grid kind '" // TRIM(kind) // "' is not recognized."
            CALL abort('init_grid_definition', message)

        END SELECT

    END SUBROUTINE init_grid_definition


    SUBROUTINE setup_io_groups()

       INTEGER ::  ngroups

       ngroups = 16
       ALLOCATE( io_group_list(ngroups) )
       
       !soil temp
       io_group_list(1) = init_io_group(                                       &
          in_files = soil_files,                                               &
          out_vars = output_var_table(1:1),                                    &
          in_var_list = input_var_table(1:1),                                  &
          kind = 'soil-temperature'                                            &
       )

       !soil water
       io_group_list(2) = init_io_group(                                       &
          in_files = soil_files,                                               &
          out_vars = output_var_table(2:2),                                    &
          in_var_list = input_var_table(2:2),                                  &
          kind = 'soil-water'                                                  &
       )

       !potential temperature, surface pressure
       io_group_list(3) = init_io_group(                                       &
          in_files = flow_files,                                               &
          out_vars = [output_var_table(3:8), output_var_table(42:42)],         &
          in_var_list = (/input_var_table(3), input_var_table(17)/),           &
          kind = 'temperature'                                                 &
       )

       !specific humidity
       io_group_list(4) = init_io_group(                                       &
          in_files = flow_files,                                               &
          out_vars = output_var_table(9:14),                                   &
          in_var_list = input_var_table(4:4),                                  &
          kind = 'scalar'                                                      &
       )

       !u and v velocity and geostrophic winds ug, vg
       io_group_list(5) = init_io_group(                                       &
          in_files = flow_files,                                               &
          out_vars = [output_var_table(15:26), output_var_table(43:44)],       &
          !out_vars = output_var_table(15:20),                                  &
          in_var_list = input_var_table(5:6),                                  &
          !in_var_list = input_var_table(5:5),                                  &
          kind = 'velocities'                                                  &
       )
   
       !!v velocity, deprecated!
       !io_group_list(6) = init_io_group(                                       &
       !   in_files = flow_files,                                               &
       !   out_vars = output_var_table(21:26),                                  &
       !   in_var_list = input_var_table(6:6),                                  &
       !   kind = 'horizontal velocity'                                         &
       !)
       !io_group_list(6) % to_be_processed = .FALSE.
   
       !w velocity
       io_group_list(7) = init_io_group(                                       &
          in_files = flow_files,                                               &
          out_vars = output_var_table(27:32),                                  &
          in_var_list = input_var_table(7:7),                                  &
          kind = 'scalar'                                                      &
       )

       !rain
       io_group_list(8) = init_io_group(                                       &
          in_files = soil_moisture_files,                                      &
          out_vars = output_var_table(33:33),                                  &
          in_var_list = input_var_table(8:8),                                  &
          kind = 'accumulated'                                                 &
       )

       !snow
       io_group_list(9) = init_io_group(                                       &
          in_files = soil_moisture_files,                                      &
          out_vars = output_var_table(34:34),                                  &
          in_var_list = input_var_table(9:9),                                  &
          kind = 'accumulated'                                                 &
       )

       !graupel
       io_group_list(10) = init_io_group(                                      &
          in_files = soil_moisture_files,                                      &
          out_vars = output_var_table(35:35),                                  &
          in_var_list = input_var_table(10:10),                                &
          kind = 'accumulated'                                                 &
       )

       !evapotranspiration
       io_group_list(11) = init_io_group(                                      &
          in_files = soil_moisture_files,                                      &
          out_vars = output_var_table(37:37),                                  &
          in_var_list = input_var_table(11:11),                                &
          kind = 'accumulated'                                                 &
       )

       !2m air temperature
       io_group_list(12) = init_io_group(                                      &
          in_files = soil_moisture_files,                                      &
          out_vars = output_var_table(36:36),                                  &
          in_var_list = input_var_table(12:12),                                &
          kind = 'surface'                                                     &
       )

       !incoming diffusive sw flux
       io_group_list(13) = init_io_group(                                      &
          in_files = radiation_files,                                          &
          out_vars = output_var_table(38:38),                                  &
          in_var_list = input_var_table(13:13),                                &
          kind = 'running average'                                             &
       )
       io_group_list(13) % to_be_processed = .FALSE.

       !incoming direct sw flux
       io_group_list(14) = init_io_group(                                      &
          in_files = radiation_files,                                          &
          out_vars = output_var_table(39:39),                                  &
          in_var_list = input_var_table(14:14),                                &
          kind = 'running average'                                             &
       )
       io_group_list(14) % to_be_processed = .FALSE.

       !sw radiation balance
       io_group_list(15) = init_io_group(                                      &
          in_files = radiation_files,                                          &
          out_vars = output_var_table(40:40),                                  &
          in_var_list = input_var_table(15:15),                                &
          kind = 'running average'                                             &
       )

       !lw radiation balance
       io_group_list(16) = init_io_group(                                      &
          in_files = radiation_files,                                          &
          out_vars = output_var_table(41:41),                                  &
          in_var_list = input_var_table(16:16),                                &
          kind = 'running average'                                             &
       )

    END SUBROUTINE setup_io_groups


    FUNCTION init_io_group(in_files, out_vars, in_var_list, kind) RESULT(group)
       CHARACTER(LEN=PATH), INTENT(IN) ::  in_files(:)
       CHARACTER(LEN=*), INTENT(IN)    ::  kind
       TYPE(nc_var), INTENT(IN)        ::  out_vars(:)
       TYPE(nc_var), INTENT(IN)        ::  in_var_list(:)

       TYPE(io_group)                  ::  group

       group % nt = SIZE(in_files)
       group % nv = SIZE(out_vars)
       group % kind = TRIM(kind)

       ALLOCATE(group % in_var_list( SIZE(in_var_list) ))
       ALLOCATE(group % in_files(group % nt))
       ALLOCATE(group % out_vars(group % nv))

       group % in_var_list = in_var_list
       group % in_files = in_files
       group % out_vars = out_vars
       group % to_be_processed = .TRUE.

    END FUNCTION init_io_group


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Deallocates all allocated variables.
!------------------------------------------------------------------------------!
    SUBROUTINE fini_grids()

       CALL report('fini_grids', 'Deallocating grids')
       
       DEALLOCATE(palm_grid%x,  palm_grid%y,  palm_grid%z,                     &
                  palm_grid%xu, palm_grid%yv, palm_grid%zw,                    &
                  palm_grid%clon,  palm_grid%clat,                             &
                  palm_grid%clonu, palm_grid%clatu)

       DEALLOCATE(palm_intermediate%x,  palm_intermediate%y,  palm_intermediate%z, &
                  palm_intermediate%xu, palm_intermediate%yv, palm_intermediate%zw,&
                  palm_intermediate%clon,  palm_intermediate%clat,             &  
                  palm_intermediate%clonu, palm_intermediate%clatu)

       DEALLOCATE(cosmo_grid%lon,  cosmo_grid%lat,  cosmo_grid%z,              &
                  cosmo_grid%lonu, cosmo_grid%latv, cosmo_grid%zw,             &
                  cosmo_grid%hfl)

    END SUBROUTINE fini_grids


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializes the the variable list.
!------------------------------------------------------------------------------!
    SUBROUTINE setup_variable_tables(mode)
       CHARACTER(LEN=*), INTENT(IN) ::  mode
       TYPE(nc_var), POINTER        ::  var

       IF (TRIM(start_date) == '')  THEN
          message = 'Simulation start date has not been set.'
          CALL abort('setup_variable_tables', message)
       END IF

       nc_source_text = 'COSMO-DE analysis from ' // TRIM(start_date)

       n_invar = 17
       n_outvar = 44
       ALLOCATE( input_var_table(n_invar) )
       ALLOCATE( output_var_table(n_outvar) )

!
!------------------------------------------------------------------------------
!- Section 1: NetCDF input variables
!------------------------------------------------------------------------------
       var => input_var_table(1)
       var % name = 'T_SO'
       var % to_be_processed = .TRUE.
       var % is_upside_down = .FALSE.

       var => input_var_table(2)
       var % name = 'W_SO'
       var % to_be_processed = .TRUE.
       var % is_upside_down = .FALSE.

       var => input_var_table(3)
       var % name = 'T'
       var % to_be_processed = .TRUE.
       var % is_upside_down = .TRUE.

       var => input_var_table(4)
       var % name = 'QV'
       var % to_be_processed = .TRUE.
       var % is_upside_down = .TRUE.

       var => input_var_table(5)
       var % name = 'U'
       var % to_be_processed = .TRUE.
       var % is_upside_down = .TRUE.

       var => input_var_table(6)
       var % name = 'V'
       var % to_be_processed = .TRUE.
       var % is_upside_down = .TRUE.

       var => input_var_table(7)
       var % name = 'W'
       var % to_be_processed = .TRUE.
       var % is_upside_down = .TRUE.

       var => input_var_table(8)
       var % name = 'RAIN_GSP'
       var % to_be_processed = .TRUE.
       var % is_upside_down = .FALSE.

       var => input_var_table(9)
       var % name = 'SNOW_GSP'
       var % to_be_processed = .TRUE.
       var % is_upside_down = .FALSE.

       var => input_var_table(10)
       var % name = 'GRAU_GSP'
       var % to_be_processed = .TRUE.
       var % is_upside_down = .FALSE.

       var => input_var_table(11)
       var % name = 'AEVAP_S'
       var % to_be_processed = .TRUE.
       var % is_upside_down = .FALSE.

       var => input_var_table(12)
       var % name = 'T_2M'
       var % to_be_processed = .TRUE.
       var % is_upside_down = .FALSE.

       var => input_var_table(13)
       var % name = 'ASWDIFD_S'
       var % to_be_processed = .TRUE.
       var % is_upside_down = .FALSE.

       var => input_var_table(14)
       var % name = 'ASWDIR_S'
       var % to_be_processed = .TRUE.
       var % is_upside_down = .FALSE.

       var => input_var_table(15)
       var % name = 'ASOB_S'
       var % to_be_processed = .TRUE.
       var % is_upside_down = .FALSE.

       var => input_var_table(16)
       var % name = 'ATHB_S'
       var % to_be_processed = .TRUE.
       var % is_upside_down = .FALSE.

       var => input_var_table(17)
       var % name = 'P'
       var % to_be_processed = .TRUE.
       var % is_upside_down = .TRUE.

!
!------------------------------------------------------------------------------
!- Section 2: NetCDF output variables
!------------------------------------------------------------------------------
       output_var_table(1) = init_nc_var(                                      &
          name              = 'init_soil_t',                                   &
          std_name          = "",                                              &
          long_name         = "initial soil temperature",                      &
          units             = "K",                                             &
          kind              = "init soil",                                     &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = palm_grid,                                       &
          intermediate_grid = palm_intermediate                                &
       )

       output_var_table(2) = init_nc_var(                                      &
          name              = 'init_soil_m',                                   &
          std_name          = "",                                              &
          long_name         = "initial soil moisture",                         &
          units             = "m^3/m^3",                                       &
          kind              = "init soil",                                     &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = palm_grid,                                       &
          intermediate_grid = palm_intermediate                                &
       )

       output_var_table(3) = init_nc_var(                                      &
          name              = 'init_pt',                                       &
          std_name          = "",                                              &
          long_name         = "initial potential temperature",                 &
          units             = "K",                                             &
          kind              = "init scalar",                                   &
          input_id          = 1,                                               & ! first in (T, p) IO group
          output_file       = output_file,                                     &
          grid              = palm_grid,                                       &
          intermediate_grid = palm_intermediate,                               &
          is_profile = (TRIM(mode) == 'profile')                               &
       )
       IF (TRIM(mode) == 'profile')  THEN
          output_var_table(3) % grid => scalar_profile_grid
          output_var_table(3) % intermediate_grid => scalar_profile_intermediate
       END IF

       output_var_table(4) = init_nc_var(                                      &
          name              = 'ls_forcing_left_pt',                            &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for left model boundary for the potential temperature", &
          units             = "K",                                             &
          kind              = "left scalar",                                   &
          input_id          = 1,                                               &
          grid              = scalars_west_grid,                               &
          intermediate_grid = scalars_west_intermediate,                       &
          output_file = output_file                                            &
       )

       output_var_table(5) = init_nc_var(                                      &
          name              = 'ls_forcing_right_pt',                           &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for right model boundary for the potential temperature", &
          units             = "K",                                             &
          kind              = "right scalar",                                  &
          input_id          = 1,                                               &
          grid              = scalars_east_grid,                               &
          intermediate_grid = scalars_east_intermediate,                       &
          output_file = output_file                                            &
       )

       output_var_table(6) = init_nc_var(                                      &
          name              = 'ls_forcing_north_pt',                           &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for north model boundary for the potential temperature", &
          units             = "K",                                             &
          kind              = "north scalar",                                  &
          input_id          = 1,                                               &
          grid              = scalars_north_grid,                              &
          intermediate_grid = scalars_north_intermediate,                      &
          output_file = output_file                                            &
       )

       output_var_table(7) = init_nc_var(                                      &
          name              = 'ls_forcing_south_pt',                           &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for south model boundary for the potential temperature", &
          units             = "K",                                             &
          kind              = "south scalar",                                  &
          input_id          = 1,                                               &
          grid              = scalars_south_grid,                              &
          intermediate_grid = scalars_south_intermediate,                      &
          output_file = output_file                                            &
       )

       output_var_table(8) = init_nc_var(                                      &
          name              = 'ls_forcing_top_pt',                             &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for top model boundary for the potential temperature", &
          units             = "K",                                             &
          kind              = "top scalar",                                    &
          input_id          = 1,                                               &
          grid              = scalars_top_grid,                                &
          intermediate_grid = scalars_top_intermediate,                        &
          output_file = output_file                                            &
       )

       output_var_table(9) = init_nc_var(                                      &
          name              = 'init_qv',                                       &
          std_name          = "",                                              &
          long_name         = "initial specific humidity",                     &
          units             = "kg/kg",                                         &
          kind              = "init scalar",                                   &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = palm_grid,                                       &
          intermediate_grid = palm_intermediate,                               &
          is_profile = (TRIM(mode) == 'profile')                               &
       )
       IF (TRIM(mode) == 'profile')  THEN
          output_var_table(9) % grid => scalar_profile_grid
          output_var_table(9) % intermediate_grid => scalar_profile_intermediate
       END IF

       output_var_table(10) = init_nc_var(                                     &
          name              = 'ls_forcing_left_qv',                            &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for left model boundary for the specific humidity", &
          units             = "kg/kg",                                         &
          kind              = "left scalar",                                   &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = scalars_west_grid,                               &
          intermediate_grid = scalars_west_intermediate                        &
       )

       output_var_table(11) = init_nc_var(                                     &
          name              = 'ls_forcing_right_qv',                           &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for right model boundary for the specific humidity", &
          units             = "kg/kg",                                         &
          kind              = "right scalar",                                  &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = scalars_east_grid,                               &
          intermediate_grid = scalars_east_intermediate                        &
       )

       output_var_table(12) = init_nc_var(                                     &
          name              = 'ls_forcing_north_qv',                           &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for north model boundary for the specific humidity", &
          units             = "kg/kg",                                         &
          kind              = "north scalar",                                  &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = scalars_north_grid,                              &
          intermediate_grid = scalars_north_intermediate                       &
       )

       output_var_table(13) = init_nc_var(                                     &
          name              = 'ls_forcing_south_qv',                           &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for south model boundary for the specific humidity", &
          units             = "kg/kg",                                         &
          kind              = "south scalar",                                  &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = scalars_south_grid,                              &
          intermediate_grid = scalars_south_intermediate                       &
       )

       output_var_table(14) = init_nc_var(                                     &
          name              = 'ls_forcing_top_qv',                             &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for top model boundary for the specific humidity", &
          units             = "kg/kg",                                         &
          kind              = "top scalar",                                    &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = scalars_top_grid,                                &
          intermediate_grid = scalars_top_intermediate                         &
       )

       output_var_table(15) = init_nc_var(                                     &
          name              = 'init_u',                                        &
          std_name          = "",                                              &
          long_name         = "initial wind component in x direction",         &
          units             = "m/s",                                           &
          kind              = "init u",                                        &
          input_id          = 1,                                               & ! first in (U, V) I/O group
          output_file       = output_file,                                     &
          grid              = u_initial_grid,                                  &
          intermediate_grid = u_initial_intermediate,                          &
          is_profile = (TRIM(mode) == 'profile')                               &
       )
       IF (TRIM(mode) == 'profile')  THEN
          output_var_table(15) % grid => scalar_profile_grid
          output_var_table(15) % intermediate_grid => scalar_profile_intermediate
       END IF 

       output_var_table(16) = init_nc_var(                                     &
          name              = 'ls_forcing_left_u',                             &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for left model boundary for the wind component in x direction", &
          units             = "m/s",                                           &
          kind              = "left u",                                        &
          input_id          = 1,                                               & ! first in (U, V) I/O group
          output_file       = output_file,                                     &
          grid              = u_west_grid,                                     &
          intermediate_grid = u_west_intermediate                              &
       )

       output_var_table(17) = init_nc_var(                                     &
          name              = 'ls_forcing_right_u',                            &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for right model boundary for the wind component in x direction", &
          units             = "m/s",                                           &
          kind              = "right u",                                       &
          input_id          = 1,                                               & ! first in (U, V) I/O group
          output_file       = output_file,                                     &
          grid              = u_east_grid,                                     &
          intermediate_grid = u_east_intermediate                              &
       )

       output_var_table(18) = init_nc_var(                                     &
          name              = 'ls_forcing_north_u',                            &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for north model boundary for the wind component in x direction", &
          units             = "m/s",                                           &
          kind              = "north u",                                       &
          input_id          = 1,                                               & ! first in (U, V) I/O group
          output_file       = output_file,                                     &
          grid              = u_north_grid,                                    &
          intermediate_grid = u_north_intermediate                             &
       )

       output_var_table(19) = init_nc_var(                                     &
          name              = 'ls_forcing_south_u',                            &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for south model boundary for the wind component in x direction", &
          units             = "m/s",                                           &
          kind              = "south u",                                       &
          input_id          = 1,                                               & ! first in (U, V) I/O group
          output_file       = output_file,                                     &
          grid              = u_south_grid,                                    &
          intermediate_grid = u_south_intermediate                             &
       )

       output_var_table(20) = init_nc_var(                                     &
          name              = 'ls_forcing_top_u',                              &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for top model boundary for the wind component in x direction", &
          units             = "m/s",                                           &
          kind              = "top u",                                         &
          input_id          = 1,                                               & ! first in (U, V) I/O group
          output_file       = output_file,                                     &
          grid              = u_top_grid,                                      &
          intermediate_grid = u_top_intermediate                               &
       )

       output_var_table(21) = init_nc_var(                                     &
          name              = 'init_v',                                        &
          std_name          = "",                                              &
          long_name         = "initial wind component in y direction",         &
          units             = "m/s",                                           &
          kind              = "init v",                                        &
          input_id          = 2,                                               & ! second in (U, V) I/O group
          output_file       = output_file,                                     &
          grid              = v_initial_grid,                                  &
          intermediate_grid = v_initial_intermediate,                          &
          is_profile = (TRIM(mode) == 'profile')                               &
       )
       IF (TRIM(mode) == 'profile')  THEN
          output_var_table(21) % grid => scalar_profile_grid
          output_var_table(21) % intermediate_grid => scalar_profile_intermediate
       END IF

       output_var_table(22) = init_nc_var(                                     &
          name              = 'ls_forcing_left_v',                             &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for left model boundary for the wind component in y direction", &
          units             = "m/s",                                           &
          kind              = "right v",                                       &
          input_id          = 2,                                               & ! second in (U, V) I/O group
          output_file       = output_file,                                     &
          grid              = v_west_grid,                                     &
          intermediate_grid = v_west_intermediate                              &
       )

       output_var_table(23) = init_nc_var(                                     &
          name              = 'ls_forcing_right_v',                            &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for right model boundary for the wind component in y direction", &
          units             = "m/s",                                           &
          kind              = "right v",                                       &
          input_id          = 2,                                               & ! second in (U, V) I/O group
          output_file       = output_file,                                     &
          grid              = v_east_grid,                                     &
          intermediate_grid = v_east_intermediate                              &
       )

       output_var_table(24) = init_nc_var(                                     &
          name              = 'ls_forcing_north_v',                            &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for north model boundary for the wind component in y direction", &
          units             = "m/s",                                           &
          kind              = "north v",                                       &
          input_id          = 2,                                               & ! second in (U, V) I/O group
          output_file       = output_file,                                     &
          grid              = v_north_grid,                                    &
          intermediate_grid = v_north_intermediate                             &
       )

       output_var_table(25) = init_nc_var(                                     &
          name              = 'ls_forcing_south_v',                            &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for south model boundary for the wind component in y direction", &
          units             = "m/s",                                           &
          kind              = "south v",                                       &
          input_id          = 2,                                               & ! second in (U, V) I/O group
          output_file       = output_file,                                     &
          grid              = v_south_grid,                                    &
          intermediate_grid = v_south_intermediate                             &
       )

       output_var_table(26) = init_nc_var(                                     &
          name              = 'ls_forcing_top_v',                              &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for top model boundary for the wind component in y direction", &
          units             = "m/s",                                           &
          kind              = "top v",                                         &
          input_id          = 2,                                               & ! second in (U, V) I/O group
          output_file       = output_file,                                     &
          grid              = v_top_grid,                                      &
          intermediate_grid = v_top_intermediate                               &
       )

       output_var_table(27) = init_nc_var(                                     &
          name              = 'init_w',                                        &
          std_name          = "",                                              &
          long_name         = "initial wind component in z direction",         &
          units             = "m/s",                                           &
          kind              = "init w",                                        &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = w_initial_grid,                                  &
          intermediate_grid = w_initial_intermediate,                          &
          is_profile = (TRIM(mode) == 'profile')                               &
       )
       IF (TRIM(mode) == 'profile')  THEN
          output_var_table(27) % grid => w_profile_grid
          output_var_table(27) % intermediate_grid => w_profile_intermediate
       END IF

       output_var_table(28) = init_nc_var(                                     &
          name              = 'ls_forcing_left_w',                             &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for left model boundary for the wind component in z direction", &
          units             = "m/s",                                           &
          kind              = "left w",                                        &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = w_west_grid,                                     &
          intermediate_grid = w_west_intermediate                              &
       )

       output_var_table(29) = init_nc_var(                                     &
          name              = 'ls_forcing_right_w',                            &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for right model boundary for the wind component in z direction", &
          units             = "m/s",                                           &
          kind              = "right w",                                       &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = w_east_grid,                                     &
          intermediate_grid = w_east_intermediate                              &
       )

       output_var_table(30) = init_nc_var(                                     &
          name              = 'ls_forcing_north_w',                            &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for north model boundary for the wind component in z direction", &
          units             = "m/s",                                           &
          kind              = "north w",                                       &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = w_north_grid,                                    &
          intermediate_grid = w_north_intermediate                             &
       )

       output_var_table(31) = init_nc_var(                                     &
          name              = 'ls_forcing_south_w',                            &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for south model boundary for the wind component in z direction", &
          units             = "m/s",                                           &
          kind              = "south w",                                       &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = w_south_grid,                                    &
          intermediate_grid = w_south_intermediate                             &
       )

       output_var_table(32) = init_nc_var(                                     &
          name              = 'ls_forcing_top_w',                              &
          std_name          = "",                                              &
          long_name         = "large-scale forcing for top model boundary for the wind component in z direction", &
          units             = "m/s",                                           &
          kind              = "top w",                                         &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = w_top_grid,                                      &
          intermediate_grid = w_top_intermediate                               &
       )

       output_var_table(33) = init_nc_var(                                     &
          name              = 'ls_forcing_soil_rain',                          &
          std_name          = "",                                              &
          long_name         = "large-scale forcing rain",                      &
          units             = "kg/m2",                                         &
          kind              = "surface forcing",                               &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = palm_grid,                                       &
          intermediate_grid = palm_intermediate                                &
       )

       output_var_table(34) = init_nc_var(                                     &
          name              = 'ls_forcing_soil_snow',                          &
          std_name          = "",                                              &
          long_name         = "large-scale forcing snow",                      &
          units             = "kg/m2",                                         &
          kind              = "surface forcing",                               &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = palm_grid,                                       &
          intermediate_grid = palm_intermediate                                &
       )

       output_var_table(35) = init_nc_var(                                     &
          name              = 'ls_forcing_soil_graupel',                       &
          std_name          = "",                                              &
          long_name         = "large-scale forcing graupel",                   &
          units             = "kg/m2",                                         &
          kind              = "surface forcing",                               &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = palm_grid,                                       &
          intermediate_grid = palm_intermediate                                &
       )

       output_var_table(36) = init_nc_var(                                     &
          name              = 'ls_forcing_soil_t_2m',                          &
          std_name          = "",                                              &
          long_name         = "large-scale forcing 2m air temperature",        &
          units             = "kg/m2",                                         &
          kind              = "surface forcing",                               &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = palm_grid,                                       &
          intermediate_grid = palm_intermediate                                &
       )

       output_var_table(37) = init_nc_var(                                     &
          name              = 'ls_forcing_soil_evap',                          &
          std_name          = "",                                              &
          long_name         = "large-scale forcing evapo-transpiration",       &
          units             = "kg/m2",                                         &
          kind              = "surface forcing",                               &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = palm_grid,                                       &
          intermediate_grid = palm_intermediate                                &
       )
       ! Radiation fluxes and balances
       output_var_table(38) = init_nc_var(                                     &
          name              = 'rad_swd_dif_0',                                 &
          std_name          = "",                                              &
          long_name         = "incoming diffuse shortwave radiative flux at the surface", &
          units             = "W/m2",                                          &
          kind              = "surface forcing",                               &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = palm_grid,                                       &
          intermediate_grid = palm_intermediate                                &
       )

       output_var_table(39) = init_nc_var(                                     &
          name              = 'rad_swd_dir_0',                                 &
          std_name          = "",                                              &
          long_name         = "incoming direct shortwave radiative flux at the surface", &
          units             = "W/m2",                                          &
          kind              = "surface forcing",                               &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = palm_grid,                                       &
          intermediate_grid = palm_intermediate                                &
       )

       output_var_table(40) = init_nc_var(                                     &
          name              = 'rad_sw_bal_0',                                  &
          std_name          = "",                                              &
          long_name         = "shortwave radiation balance at the surface",    &
          units             = "W/m2",                                          &
          kind              = "surface forcing",                               &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = palm_grid,                                       &
          intermediate_grid = palm_intermediate                                &
       )

       output_var_table(41) = init_nc_var(                                     &
          name              = 'rad_lw_bal_0',                                  &
          std_name          = "",                                              &
          long_name         = "longwave radiation balance at the surface",     &
          units             = "W/m2",                                          &
          kind              = "surface forcing",                               &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = palm_grid,                                       &
          intermediate_grid = palm_intermediate                                &
       )

       output_var_table(42) = init_nc_var(                                     &
          name              = 'surface_forcing_surface_pressure',              &
          std_name          = "",                                              &
          long_name         = "surface pressure",                              &
          units             = "Pa",                                            &
          kind              = "time series",                                   &
          input_id          = 2,                                               & ! second in (T, p) I/O group
          output_file       = output_file,                                     &
          grid              = palm_grid,                                       &
          intermediate_grid = palm_intermediate                                &
       )

       output_var_table(43) = init_nc_var(                                     &
          name              = 'ls_forcing_ug',                                 &
          std_name          = "",                                              &
          long_name         = "geostrophic wind (u component)",                &
          units             = "m/s",                                           &
          kind              = "profile",                                       &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = palm_grid,                                       &
          intermediate_grid = palm_intermediate                                &
       )

       output_var_table(44) = init_nc_var(                                     &
          name              = 'ls_forcing_vg',                                 &
          std_name          = "",                                              &
          long_name         = "geostrophic wind (v component)",                &
          units             = "m/s",                                           &
          kind              = "profile",                                       &
          input_id          = 1,                                               &
          output_file       = output_file,                                     &
          grid              = palm_grid,                                       &
          intermediate_grid = palm_intermediate                                &
       )

       ! Attributes shared among all variables
       output_var_table(:) % source = nc_source_text

    END SUBROUTINE setup_variable_tables


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializes and nc_var varible with the given parameters. The 'kind'
!> parameter is used to infer the correct netCDF IDs and the level of detail,
!> 'lod', as defined by the PALM-4U input data standard.
!------------------------------------------------------------------------------!
    FUNCTION init_nc_var(name, std_name, long_name, units, kind, input_id,     &
                         grid, intermediate_grid, output_file, is_profile) RESULT(var)

       CHARACTER(LEN=*), INTENT(IN)      ::  name, std_name, long_name, units, kind
       INTEGER, INTENT(IN)               ::  input_id
       TYPE(grid_definition), INTENT(IN), TARGET ::  grid, intermediate_grid
       TYPE(nc_file), INTENT(IN)         ::  output_file
       LOGICAL, INTENT(IN), OPTIONAL     ::  is_profile

       CHARACTER(LEN=LNAME)              ::  out_var_kind 
       TYPE(nc_var)                      ::  var

       out_var_kind = TRIM(kind)

       IF (PRESENT(is_profile))  THEN
          IF (is_profile)  out_var_kind = TRIM(kind) // ' profile'
       END IF

       var % name              = name
       var % standard_name     = std_name
       var % long_name         = long_name
       var % units             = units
       var % kind              = TRIM(out_var_kind)
       var % input_id          = input_id
       var % nt                = SIZE (output_file % time)
       var % grid              => grid
       var % intermediate_grid => intermediate_grid

       SELECT CASE( TRIM(out_var_kind) )

       !TODO: Using global module variables 'init_variables_required' and
       !TODO: 'boundary_variables_required'. Encapsulate in settings type
       !TODO: and pass into init_nc_var.
       CASE( 'init soil' )
          var % nt              = 1
          var % lod             = 2
          var % ndim            = 3
          var % dimids(1:3)     = output_file % dimids_soil
          var % dimvarids(1:3)  = output_file % dimvarids_soil
          var % to_be_processed = init_variables_required
          var % task            = "interpolate_2d"

       CASE( 'init scalar' )
          var % nt              = 1
          var % lod             = 2
          var % ndim            = 3
          var % dimids(1:3)     = output_file % dimids_scl
          var % dimvarids(1:3)  = output_file % dimvarids_scl
          var % to_be_processed = init_variables_required
          var % task            = "interpolate_3d"

       CASE( 'init u' )
          var % nt              = 1
          var % lod             = 2
          var % ndim            = 3
          var % dimids(1)       = output_file % dimids_vel(1)
          var % dimids(2)       = output_file % dimids_scl(2)
          var % dimids(3)       = output_file % dimids_scl(3)
          var % dimvarids(1)    = output_file % dimvarids_vel(1)
          var % dimvarids(2)    = output_file % dimvarids_scl(2)
          var % dimvarids(3)    = output_file % dimvarids_scl(3)
          var % to_be_processed = init_variables_required
          var % task            = "interpolate_3d"

       CASE( 'init v' )
          var % nt              = 1
          var % lod             = 2
          var % ndim            = 3
          var % dimids(1)       = output_file % dimids_scl(1)
          var % dimids(2)       = output_file % dimids_vel(2)
          var % dimids(3)       = output_file % dimids_scl(3)
          var % dimvarids(1)    = output_file % dimvarids_scl(1)
          var % dimvarids(2)    = output_file % dimvarids_vel(2)
          var % dimvarids(3)    = output_file % dimvarids_scl(3)
          var % to_be_processed = init_variables_required
          var % task            = "interpolate_3d"

       CASE( 'init w' )
          var % nt              = 1
          var % lod             = 2
          var % ndim            = 3
          var % dimids(1)       = output_file % dimids_scl(1)
          var % dimids(2)       = output_file % dimids_scl(2)
          var % dimids(3)       = output_file % dimids_vel(3)
          var % dimvarids(1)    = output_file % dimvarids_scl(1)
          var % dimvarids(2)    = output_file % dimvarids_scl(2)
          var % dimvarids(3)    = output_file % dimvarids_vel(3)
          var % to_be_processed = init_variables_required
          var % task            = "interpolate_3d"

       CASE( 'init scalar profile', 'init u profile', 'init v profile')
          var % nt              = 1
          var % lod             = 1
          var % ndim            = 1
          var % dimids(1)       = output_file % dimids_scl(3)    !z
          var % dimvarids(1)    = output_file % dimvarids_scl(3) !z
          var % to_be_processed = init_variables_required
          var % task            = "average profile"

       CASE( 'init w profile')
          var % nt              = 1
          var % lod             = 1
          var % ndim            = 1
          var % dimids(1)       = output_file % dimids_vel(3)    !z
          var % dimvarids(1)    = output_file % dimvarids_vel(3) !z
          var % to_be_processed = init_variables_required
          var % task            = "average profile"

       CASE ( 'surface forcing' )
          var % lod             = 1
          var % ndim            = 3
          var % dimids(3)       = output_file % dimid_time
          var % dimids(1:2)     = output_file % dimids_soil(1:2)
          var % dimvarids(3)    = output_file % dimvarid_time
          var % dimvarids(1:2)  = output_file % dimvarids_soil(1:2)
          var % to_be_processed = boundary_variables_required
          var % task            = "interpolate_2d"

       CASE ( 'left scalar', 'right scalar') ! same as right
          var % lod             = 2
          var % ndim            = 3
          var % dimids(3)       = output_file % dimid_time
          var % dimids(1)       = output_file % dimids_scl(2)
          var % dimids(2)       = output_file % dimids_scl(3)
          var % dimvarids(3)    = output_file % dimvarid_time
          var % dimvarids(1)    = output_file % dimvarids_scl(2)
          var % dimvarids(2)    = output_file % dimvarids_scl(3)
          var % to_be_processed = boundary_variables_required
          var % task            = "interpolate_3d"

       CASE ( 'north scalar', 'south scalar') ! same as south
          var % lod             = 2
          var % ndim            = 3
          var % dimids(3)       = output_file % dimid_time
          var % dimids(1)       = output_file % dimids_scl(1)
          var % dimids(2)       = output_file % dimids_scl(3)
          var % dimvarids(3)    = output_file % dimvarid_time
          var % dimvarids(1)    = output_file % dimvarids_scl(1)
          var % dimvarids(2)    = output_file % dimvarids_scl(3)
          var % to_be_processed = boundary_variables_required
          var % task            = "interpolate_3d"

       CASE ( 'top scalar', 'top w' )
          var % lod             = 2
          var % ndim            = 3
          var % dimids(3)       = output_file % dimid_time
          var % dimids(1)       = output_file % dimids_scl(1)
          var % dimids(2)       = output_file % dimids_scl(2)
          var % dimvarids(3)    = output_file % dimvarid_time
          var % dimvarids(1)    = output_file % dimvarids_scl(1)
          var % dimvarids(2)    = output_file % dimvarids_scl(2)
          var % to_be_processed = boundary_variables_required
          var % task            = "interpolate_3d"

       CASE ( 'left u', 'right u' )
          var % lod             = 2
          var % ndim            = 3
          var % dimids(3)       = output_file % dimid_time
          var % dimids(1)       = output_file % dimids_scl(2)
          var % dimids(2)       = output_file % dimids_scl(3)
          var % dimvarids(3)    = output_file % dimvarid_time
          var % dimvarids(1)    = output_file % dimvarids_scl(2)
          var % dimvarids(2)    = output_file % dimvarids_scl(3)
          var % to_be_processed = boundary_variables_required
          var % task            = "interpolate_3d"

       CASE ( 'north u', 'south u' )
          var % lod             = 2
          var % ndim            = 3
          var % dimids(3)       = output_file % dimid_time    !t
          var % dimids(1)       = output_file % dimids_vel(1) !x
          var % dimids(2)       = output_file % dimids_scl(3) !z
          var % dimvarids(3)    = output_file % dimvarid_time
          var % dimvarids(1)    = output_file % dimvarids_vel(1)
          var % dimvarids(2)    = output_file % dimvarids_scl(3)
          var % to_be_processed = boundary_variables_required
          var % task            = "interpolate_3d"

       CASE ( 'top u' )
          var % lod             = 2
          var % ndim            = 3
          var % dimids(3)       = output_file % dimid_time    !t
          var % dimids(1)       = output_file % dimids_vel(1) !x
          var % dimids(2)       = output_file % dimids_scl(2) !z
          var % dimvarids(3)    = output_file % dimvarid_time
          var % dimvarids(1)    = output_file % dimvarids_vel(1)
          var % dimvarids(2)    = output_file % dimvarids_scl(2)
          var % to_be_processed = boundary_variables_required
          var % task            = "interpolate_3d"

       CASE ( 'left v', 'right v' )
          var % lod             = 2
          var % ndim            = 3
          var % dimids(3)       = output_file % dimid_time
          var % dimids(1)       = output_file % dimids_vel(2)
          var % dimids(2)       = output_file % dimids_scl(3)
          var % dimvarids(3)    = output_file % dimvarid_time
          var % dimvarids(1)    = output_file % dimvarids_vel(2)
          var % dimvarids(2)    = output_file % dimvarids_scl(3)
          var % to_be_processed = boundary_variables_required
          var % task            = "interpolate_3d"

       CASE ( 'north v', 'south v' )
          var % lod             = 2
          var % ndim            = 3
          var % dimids(3)       = output_file % dimid_time    !t
          var % dimids(1)       = output_file % dimids_scl(1) !x
          var % dimids(2)       = output_file % dimids_scl(3) !z
          var % dimvarids(3)    = output_file % dimvarid_time
          var % dimvarids(1)    = output_file % dimvarids_scl(1)
          var % dimvarids(2)    = output_file % dimvarids_scl(3)
          var % to_be_processed = boundary_variables_required
          var % task            = "interpolate_3d"

       CASE ( 'top v' )
          var % lod             = 2
          var % ndim            = 3
          var % dimids(3)       = output_file % dimid_time    !t
          var % dimids(1)       = output_file % dimids_scl(1) !x
          var % dimids(2)       = output_file % dimids_vel(2) !z
          var % dimvarids(3)    = output_file % dimvarid_time
          var % dimvarids(1)    = output_file % dimvarids_scl(1)
          var % dimvarids(2)    = output_file % dimvarids_vel(2)
          var % to_be_processed = boundary_variables_required
          var % task            = "interpolate_3d"

       CASE ( 'left w', 'right w')
          var % lod             = 2
          var % ndim            = 3
          var % dimids(3)       = output_file % dimid_time
          var % dimids(1)       = output_file % dimids_scl(2)
          var % dimids(2)       = output_file % dimids_vel(3)
          var % dimvarids(3)    = output_file % dimvarid_time
          var % dimvarids(1)    = output_file % dimvarids_scl(2)
          var % dimvarids(2)    = output_file % dimvarids_vel(3)
          var % to_be_processed = boundary_variables_required
          var % task            = "interpolate_3d"

       CASE ( 'north w', 'south w' )
          var % lod             = 2
          var % ndim            = 3
          var % dimids(3)       = output_file % dimid_time    !t
          var % dimids(1)       = output_file % dimids_scl(1) !x
          var % dimids(2)       = output_file % dimids_vel(3) !z
          var % dimvarids(3)    = output_file % dimvarid_time
          var % dimvarids(1)    = output_file % dimvarids_scl(1)
          var % dimvarids(2)    = output_file % dimvarids_vel(3)
          var % to_be_processed = boundary_variables_required
          var % task            = "interpolate_3d"

       CASE ( 'time series' )
          var % lod             = 0
          var % ndim            = 1
          var % dimids(1)       = output_file % dimid_time    !t
          var % dimvarids(1)    = output_file % dimvarid_time
          var % to_be_processed = .TRUE.
          var % task            = "average scalar"

       CASE ( 'profile' )
          var % lod             = 2
          var % ndim            = 2
          var % dimids(2)       = output_file % dimid_time    !t
          var % dimids(1)       = output_file % dimids_scl(3) !z
          var % dimvarids(2)    = output_file % dimvarid_time
          var % dimvarids(1)    = output_file % dimvarids_scl(3)
          var % to_be_processed = .TRUE.
          var % task            = "profile"

       CASE DEFAULT
           message = "Variable kind '" // TRIM(kind) // "' not recognized."
           CALL abort ('init_nc_var', message)

       END SELECT

    END FUNCTION init_nc_var


    SUBROUTINE fini_variables()

       CALL report('fini_variables', 'Deallocating variable table')
       DEALLOCATE( input_var_table )

    END SUBROUTINE fini_variables


    SUBROUTINE fini_io_groups()

       CALL report('fini_io_groups', 'Deallocating IO groups')
       DEALLOCATE( io_group_list )

    END SUBROUTINE fini_io_groups


    SUBROUTINE fini_file_lists()
       
       CALL report('fini_file_lists', 'Deallocating file lists')
       DEALLOCATE( flow_files, soil_files, radiation_files, soil_moisture_files )

    END SUBROUTINE fini_file_lists


    SUBROUTINE input_file_list(start_date_string, start_hour, end_hour,        &
                               step_hour, path, prefix, suffix, file_list)

       CHARACTER (LEN=DATE), INTENT(IN) ::  start_date_string
       CHARACTER (LEN=*),    INTENT(IN) ::  prefix, suffix, path
       INTEGER,              INTENT(IN) ::  start_hour, end_hour, step_hour
       CHARACTER(LEN=*), ALLOCATABLE, INTENT(INOUT) ::  file_list(:)

       INTEGER             ::  number_of_files, hour, i
       CHARACTER(LEN=DATE) ::  date_string

       number_of_files = end_hour - start_hour + 1

       ALLOCATE( file_list(number_of_files) )

       DO i = 1, number_of_files
          hour = start_hour + (i-1) * step_hour
          date_string = add_hours_to(start_date_string, hour)

          file_list(i) = TRIM(path) // TRIM(prefix) // TRIM(date_string) //    &
                         TRIM(suffix) // '.nc'
          message = "Set up input file name '" // TRIM(file_list(i)) //"'"
          CALL report('input_file_list', message)
       END DO

    END SUBROUTINE input_file_list


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Carries out any physical conversion of the quantities in the given input
!> buffer needed to obtain the quantity required by PALM-4U. For instance, 
!> velocities are rotated to the PALM-4U coordinate system and the potential
!> temperature is computed from the absolute temperature and pressure.
!>
!> Note, that the preprocessing does not include any grid change. The result
!> array will match a COSMO-DE scalar array.
!------------------------------------------------------------------------------!
    SUBROUTINE preprocess(group, input_buffer, cosmo_grid, iter)

       TYPE(io_group), INTENT(INOUT), TARGET       ::  group
       TYPE(container), INTENT(INOUT), ALLOCATABLE ::  input_buffer(:)
       TYPE(grid_definition), INTENT(IN)           ::  cosmo_grid
       INTEGER, INTENT(IN)                         ::  iter
       
       REAL(dp), ALLOCATABLE                       ::  basic_state_pressure(:)
       TYPE(container), ALLOCATABLE                ::  compute_buffer(:)
       INTEGER                                     ::  hour, dt
       INTEGER                                     ::  i, j, k
       INTEGER                                     ::  nx, ny, nz
       
       input_buffer(:) % is_preprocessed = .FALSE.
        
       SELECT CASE( group % kind )
       
       CASE( 'velocities' )
          ! Allocate a compute puffer with the same number of arrays as the input
          ALLOCATE( compute_buffer( SIZE(input_buffer) ) )

          ! Allocate u and v arrays with scalar dimensions
          nx = SIZE(input_buffer(1) % array, 1)
          ny = SIZE(input_buffer(1) % array, 2)
          nz = SIZE(input_buffer(1) % array, 3)
          ALLOCATE( compute_buffer(1) % array(nx, ny, nz) ) ! u buffer
          ALLOCATE( compute_buffer(2) % array(nx, ny, nz) ) ! v buffer

 CALL run_control('time', 'alloc')

          ! interpolate U and V to centres
          CALL centre_velocities( u_face = input_buffer(1) % array,            &
                                  v_face = input_buffer(2) % array,            &
                                  u_centre = compute_buffer(1) % array,        &
                                  v_centre = compute_buffer(2) % array )
          
          ! rotate U and V to PALM-4U orientation and overwrite U and V with
          ! rotated velocities
          DO k = 1, nz
          DO j = 2, ny
          DO i = 2, nx
             CALL uv2uvrot( urot = compute_buffer(1) % array(i,j,k),           &
                            vrot = compute_buffer(2) % array(i,j,k),           &
                            rlat = cosmo_grid % lat(j-1),                      &
                            rlon = cosmo_grid % lon(i-1),                      &
                            pollat = phi_cn,                                   &
                            pollon = lambda_cn,                                &
                            u = input_buffer(1) % array(i,j,k),                &
                            v = input_buffer(2) % array(i,j,k) )
          END DO
          END DO
          END DO
          input_buffer(1) % array(1,:,:) = 0.0_dp
          input_buffer(2) % array(1,:,:) = 0.0_dp
          input_buffer(1) % array(:,1,:) = 0.0_dp
          input_buffer(2) % array(:,1,:) = 0.0_dp

          input_buffer(1:2) % is_preprocessed = .TRUE.
 CALL run_control('time', 'comp')

          DEALLOCATE( compute_buffer )
 CALL run_control('time', 'alloc')

          message = "Input buffers for group '" // TRIM(group % kind) // "'"//&
             " preprocessed sucessfully."
          CALL report('preprocess', message)
       
       CASE( 'temperature' ) ! P and T
          nx = SIZE(input_buffer(1) % array, 1)
          ny = SIZE(input_buffer(1) % array, 2)
          nz = SIZE(input_buffer(1) % array, 3)

          ! Compute absolute pressure if presure perturbation has been read in.
          IF ( TRIM(group % in_var_list(2) % name) == 'PP' )  THEN
             message = "Absolute pressure, P, not available, " //              &
                       "computing from pressure preturbation PP."
             CALL report('preprocess', message)

             ALLOCATE( basic_state_pressure(1:nz) )
 CALL run_control('time', 'alloc')

             DO j = 1, ny
             DO i = 1, nx

                CALL get_basic_state(cosmo_grid % hfl(i,j,:), BETA, P_SL, T_SL,&
                                     RD, G, basic_state_pressure)

                input_buffer (2) % array(i,j,:) =                              &
                   input_buffer (2) % array(i,j,:) + basic_state_pressure(:)

             END DO
             END DO
 CALL run_control('time', 'comp')

             DEALLOCATE( basic_state_pressure )
 CALL run_control('time', 'alloc')

             group % in_var_list(2) % name = 'P'

          END IF

          ! Convert absolute to potential temperature
          input_buffer(1) % array(:,:,:) = input_buffer(1) % array(:,:,:) *    &
             EXP( RD_PALM/CP_PALM * LOG( P_REF / input_buffer(2) % array(:,:,:) ))

          input_buffer(1:2) % is_preprocessed = .TRUE.
          message = "Input buffers for group '" // TRIM(group % kind) // "'"//&
             " preprocessed sucessfully."
          CALL report('preprocess', message)
 CALL run_control('time', 'comp')
       
       CASE( 'scalar' ) ! S or W
          input_buffer(:) % is_preprocessed = .TRUE.
 CALL run_control('time', 'comp')

       CASE( 'soil-temperature' ) ! 
          
          CALL fill_water_cells(soiltyp, input_buffer(1) % array, &
                                SIZE(input_buffer(1) % array, 3), &
                                FILL_ITERATIONS)
          input_buffer(:) % is_preprocessed = .TRUE.
 CALL run_control('time', 'comp')

       CASE( 'soil-water' ) ! 

          CALL fill_water_cells(soiltyp, input_buffer(1) % array, &
                                SIZE(input_buffer(1) % array, 3), &
                                FILL_ITERATIONS)

          nx = SIZE(input_buffer(1) % array, 1)
          ny = SIZE(input_buffer(1) % array, 2)
          nz = SIZE(input_buffer(1) % array, 3)

          DO k = 1, nz
          DO j = 1, ny
          DO i = 1, nx
             input_buffer(1) % array(i,j,k) =                                  &
                 input_buffer(1) % array(i,j,k) * d_depth_rho_inv(k)
          END DO
          END DO
          END DO

          message = "Converted soil water from [kg/m^2] to [m^3/m^3]"
          CALL report('preprocess', message)

          input_buffer(:) % is_preprocessed = .TRUE.
 CALL run_control('time', 'comp')

       CASE( 'surface' ) ! 
          input_buffer(:) % is_preprocessed = .TRUE.
 CALL run_control('time', 'comp')

       CASE( 'accumulated' ) ! 
          message = "De-accumulating '" // TRIM(group % in_var_list(1) % name) //&
                    "' in iteration " // TRIM(str(iter)) 
          CALL report('preprocess', message)

          hour = iter - 1
          dt = MODULO(hour, 3) + 1 ! averaging period
          SELECT CASE(dt)

          CASE(1)
             !input has been accumulated over one hour. Leave as is
             !input_buffer(1) % array(:,:,:) carrries one-hour integral

          CASE(2)
             !input has been accumulated over two hours. Subtract previous step
             !input_buffer(1) % array(:,:,:) carrries one-hour integral
             !input_buffer(2) % array(:,:,:) carrries two-hour integral
             CALL deaverage(                                                   &
                      avg_1 = input_buffer(1) % array(:,:,:), t1 = 1.0_dp,     &
                      avg_2 = input_buffer(2) % array(:,:,:), t2 = 1.0_dp,     &
                      avg_3 = input_buffer(1) % array(:,:,:), t3 = 1.0_dp )
             !input_buffer(1) % array(:,:,:) carrries one-hour integral of second hour

          CASE(3)
             !input has been accumulated over three hours. Subtract previous step
             !input_buffer(1) % array(:,:,:) carrries three-hour integral
             !input_buffer(2) % array(:,:,:) still carrries two-hour integral
             CALL deaverage(                                                   &
                     avg_1 = input_buffer(2) % array(:,:,:), t1 = 1.0_dp,      &
                     avg_2 = input_buffer(1) % array(:,:,:), t2 = 1.0_dp,      &
                     avg_3 = input_buffer(1) % array(:,:,:), t3 = 1.0_dp )
             !input_buffer(1) % array(:,:,:) carrries one-hour integral of third hourA

          CASE DEFAULT
             message = "Invalid averaging period '" // TRIM(str(dt)) // " hours"
             CALL abort('preprocess', message)

          END SELECT
          input_buffer(:) % is_preprocessed = .TRUE.
 CALL run_control('time', 'comp')

       CASE( 'running average' ) ! 
          message = "De-averaging '" // TRIM(group % in_var_list(1) % name) //   &
                    "' in iteration " // TRIM(str(iter)) 
          CALL report('preprocess', message)

          hour = iter - 1
          dt = MODULO(hour, 3) + 1 ! averaging period
          SELECT CASE(dt)

          CASE(1)
             !input has been accumulated over one hour. Leave as is
             !input_buffer(1) % array(:,:,:) carrries one-hour integral

          CASE(2)
             !input has been accumulated over two hours. Subtract previous step
             !input_buffer(1) % array(:,:,:) carrries one-hour integral
             !input_buffer(2) % array(:,:,:) carrries two-hour integral
             CALL deaverage( input_buffer(1) % array(:,:,:), 1.0_dp,           &
                             input_buffer(2) % array(:,:,:), 2.0_dp,           &
                             input_buffer(1) % array(:,:,:), 1.0_dp)
             !input_buffer(1) % array(:,:,:) carrries one-hour integral of second hour

          CASE(3)
             !input has been accumulated over three hours. Subtract previous step
             !input_buffer(1) % array(:,:,:) carrries three-hour integral
             !input_buffer(2) % array(:,:,:) still carrries two-hour integral
             CALL deaverage( input_buffer(2) % array(:,:,:), 2.0_dp,           &
                             input_buffer(1) % array(:,:,:), 3.0_dp,           &
                             input_buffer(1) % array(:,:,:), 1.0_dp)
             !input_buffer(1) % array(:,:,:) carrries one-hour integral of third hourA

          CASE DEFAULT
             message = "Invalid averaging period '" // TRIM(str(dt)) // " hours"
             CALL abort('preprocess', message)

          END SELECT
          input_buffer(:) % is_preprocessed = .TRUE.

       CASE DEFAULT
          message = "Buffer kind '" // TRIM(group % kind) // "' is not supported."
          CALL abort('prerpocess', message)

       END SELECT
 CALL run_control('time', 'comp')

    END SUBROUTINE preprocess


! Description:
! ------------
!> Computes average soil values in COSMO-DE water cells from neighbouring
!> non-water cells. This is done as a preprocessing step for the COSMO-DE
!> soil input arrays, which contain unphysical values for water cells.
!>
!> This routine computes the average of up to all nine neighbouring cells
!> or keeps the original value, if not at least one non-water neightbour
!> is available.
!>
!> By repeatedly applying this step, soil data can be extrapolated into
!> 'water' regions occupying multiple cells, one cell per iteration.
!> 
!> Input parameters:
!> -----------------
!> soiltyp : 2d map of COSMO-DE soil types
!> nz : number of layers in the COSMO-DE soil
!> niter : number iterations
!> 
!> Output parameters:
!> ------------------
!> array : the soil array (i.e. water content or temperature)
!------------------------------------------------------------------------------!
    SUBROUTINE fill_water_cells(soiltyp, array, nz, niter)
       INTEGER(hp), DIMENSION(:,:,:), INTENT(IN) :: soiltyp
       REAL(dp), DIMENSION(:,:,:), INTENT(INOUT) :: array
       INTEGER, INTENT(IN)                       :: nz, niter

       REAL(dp), DIMENSION(nz)                   :: column
       INTEGER(hp), DIMENSION(:,:), ALLOCATABLE  :: old_soiltyp, new_soiltyp
       INTEGER                                   :: l, i, j, nx, ny, n_cells, ii, jj, iter
       INTEGER, DIMENSION(8)                     :: di, dj

       nx = SIZE(array, 1)
       ny = SIZE(array, 2)
       di = (/ -1, -1, -1, 0,  0,  1, 1, 1 /)
       dj = (/ -1,  0,  1, -1, 1, -1, 0, 1 /)

       ALLOCATE(old_soiltyp(SIZE(soiltyp,1), &
                            SIZE(soiltyp,2) ))

       ALLOCATE(new_soiltyp(SIZE(soiltyp,1), &
                            SIZE(soiltyp,2) ))

       old_soiltyp(:,:) = soiltyp(:,:,1)
       new_soiltyp(:,:) = soiltyp(:,:,1)

       DO iter = 1, niter

          DO j = 1, ny
          DO i = 1, nx
          
             IF (old_soiltyp(i,j) == WATER_ID)  THEN 

                n_cells = 0
                column(:) = 0.0_dp
                DO l = 1, SIZE(di)

                   ii = MIN(nx, MAX(1, i + di(l)))
                   jj = MIN(ny, MAX(1, j + dj(l)))

                   IF (old_soiltyp(ii,jj) .NE. WATER_ID)  THEN
                      n_cells = n_cells + 1
                      column(:) = column(:) + array(ii,jj,:)
                   END IF

                END DO

                ! Overwrite if at least one non-water neighbour cell is available
                IF (n_cells > 0)  THEN
                   array(i,j,:) = column(:) / n_cells
                   new_soiltyp(i,j) = 0
                END IF

             END IF

          END DO
          END DO

          old_soiltyp(:,:) = new_soiltyp(:,:)

       END DO

       DEALLOCATE(old_soiltyp, new_soiltyp)

    END SUBROUTINE fill_water_cells

 END MODULE grid
