!> @file src/types.f90
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
! $Id: types.f90 2718 2018-01-02 08:49:38Z maronga $
! Initial revision
!
! 
!
! Authors:
! --------
! @author Eckhard Kadasch
!
! Description:
! ------------
!> The types module provides derived data types used in INIFOR.
!------------------------------------------------------------------------------!
 MODULE types
 
 USE defs,                                                                     &
    ONLY:  dp, PATH, SNAME, LNAME
 USE netcdf,                                                                   &
    ONLY:  NF90_MAX_VAR_DIMS, NF90_MAX_NAME

 IMPLICIT NONE

 TYPE grid_definition
    CHARACTER(LEN=SNAME)  ::  name(3)       !< names of the grid dimensions, e.g. (/'x', 'y', 'z'/) or (/'latitude', 'longitude', 'height'/)
    INTEGER               ::  nx            !< number of gridpoints in the first dimension
    INTEGER               ::  ny            !< number of gridpoints in the second dimension
    INTEGER               ::  nz            !< number of gridpoints in the third dimension
    INTEGER, ALLOCATABLE  ::  ii(:,:,:)     !< Given a point (i,j,k) in the PALM-4U grid, ii(i,j,l) gives the x index of the l'th horizontl neighbour on the COSMO-DE grid.
    INTEGER, ALLOCATABLE  ::  jj(:,:,:)     !< Given a point (i,j,k) in the PALM-4U grid, jj(i,j,l) gives the y index of the l'th horizontl neighbour on the COSMO-DE grid.
    INTEGER, ALLOCATABLE  ::  kk(:,:,:,:)   !< Given a point (i,j,k) in the PALM-4U grid, kk(i,j,k,l) gives the z index of the l'th vertical neighbour in the intermediate grid.
    REAL(dp)              ::  dx            !< grid spacing in the first dimension [m]
    REAL(dp)              ::  dy            !< grid spacing in the second dimension [m]
    REAL(dp)              ::  dz            !< grid spacing in the third dimension [m]
    REAL(dp)              ::  dxi           !< inverse grid spacing in the first dimension [m^-1]
    REAL(dp)              ::  dyi           !< inverse grid spacing in the second dimension [m^-1]
    REAL(dp)              ::  dzi           !< inverse grid spacing in the third dimension [m^-1]
    REAL(dp)              ::  lx            !< domain length in the first dimension [m]
    REAL(dp)              ::  ly            !< domain length in the second dimension [m]
    REAL(dp)              ::  lz            !< domain length in the third dimension [m]
    REAL(dp)              ::  x0            !< x coordinate of PALM-4U domain projection centre, i.e. location of zero distortion
    REAL(dp)              ::  y0            !< y coordinate of PALM-4U domain projection centre, i.e. location of zwro distortion
    REAL(dp)              ::  z0            !< displacement of the coordinate origin above sea level [m]
    REAL(dp), ALLOCATABLE ::  x(:)          !< coordinates of cell centers in x direction [m]
    REAL(dp), ALLOCATABLE ::  y(:)          !< coordinates of cell centers in y direction [m]
    REAL(dp), ALLOCATABLE ::  z(:)          !< coordinates of cell centers in z direction [m]
    REAL(dp), ALLOCATABLE ::  h(:,:,:)      !< heights grid point for intermediate grids [m]
    REAL(dp), POINTER     ::  hhl(:,:,:)    !< heights of half layers (cell faces) above sea level in COSMO-DE, read in from 
    REAL(dp), POINTER     ::  hfl(:,:,:)    !< heights of full layers (cell centres) above sea level in COSMO-DE, computed as arithmetic average of hhl
    REAL(dp), POINTER     ::  depths(:)     !< depths of output soil layers, equal the depths of the source model (e.g. COSMO-DE)
    REAL(dp), ALLOCATABLE ::  xu(:)         !< coordinates of cell faces in x direction [m]
    REAL(dp), ALLOCATABLE ::  yv(:)         !< coordinates of cell faces in y direction [m]
    REAL(dp), ALLOCATABLE ::  zw(:)         !< coordinates of cell faces in z direction [m]
    REAL(dp), ALLOCATABLE ::  lat(:)        !< rotated-pole latitudes of scalars (cell centers) of the COSMO-DE grid [rad]
    REAL(dp), ALLOCATABLE ::  lon(:)        !< rotated-pole longitudes of scalars (cell centres) of the COSMO-DE grid [rad]
    REAL(dp), ALLOCATABLE ::  latv(:)       !< rotated-pole latitudes of v winds (face centres in latitudal/y direction) [rad]
    REAL(dp), ALLOCATABLE ::  lonu(:)       !< rotated-pole latitudes of u winds (face centres in longitudal/x direction) [rad]
    REAL(dp), ALLOCATABLE ::  clat(:,:)     !< latitudes of PALM-4U cell centres in COSMO-DE's rotated-pole grid [rad]
    REAL(dp), ALLOCATABLE ::  clon(:,:)     !< longitudes of PALM-4U scalars (cell centres) in COSMO-DE's rotated-pole grid [rad]
    REAL(dp), ALLOCATABLE ::  clatu(:,:)    !< latitudes of PALM-4U u winds (cell faces in u direction) in COSMO-DE's rotated-pole grid [rad]
    REAL(dp), ALLOCATABLE ::  clonu(:,:)    !< longitudes of PALM-4U u winds (cell faces in u direction) in COSMO-DE's rotated-pole grid [rad]
    REAL(dp), ALLOCATABLE ::  clatv(:,:)    !< latitudes of PALM-4U v winds (cell faces in v direction) in COSMO-DE's rotated-pole grid [rad]
    REAL(dp), ALLOCATABLE ::  clonv(:,:)    !< longitudes of PALM-4U v winds (cell faces in v direction) in COSMO-DE's rotated-pole grid [rad]
    REAL(dp), ALLOCATABLE ::  w_horiz(:,:,:)   !< weights for bilinear horizontal interpolation
    REAL(dp), ALLOCATABLE ::  w_verti(:,:,:,:) !< weights for linear vertical interpolation
 END TYPE


 TYPE nc_file
    CHARACTER(LEN=PATH)   ::  name          !< file name
    INTEGER               ::  dimid_time    !< NetCDF IDs of the time dimension
    INTEGER               ::  dimids_scl(3) !< NetCDF IDs of the grid dimensions for scalar points x, y, z 
    INTEGER               ::  dimids_vel(3) !< NetCDF IDs of the grid dimensions for velocity points xu, yu, zu
    INTEGER               ::  dimids_soil(3)!< NetCDF IDs of the grid dimensions for soil points x, y, depth
    INTEGER               ::  dimvarid_time !< NetCDF IDs of the time variable
    INTEGER               ::  dimvarids_scl(3) !< NetCDF IDs of the grid coordinates of scalars x, y, z 
    INTEGER               ::  dimvarids_vel(3) !< NetCDF IDs of the grid coordinates of velocities xu, yu, zu. Note that velocities are located at mix of both coordinates, e.g. u(xu, y, z).
    INTEGER               ::  dimvarids_soil(3)!< NetCDF IDs of the grid coordinates for soil points x, y, depth 
    REAL(dp), POINTER     ::  time(:)       ! vector of output time steps
 END TYPE


 TYPE nc_var
    INTEGER                               ::  varid     !< NetCDF ID of the variable
    INTEGER                               ::  input_id  !< ID of the correpsonding input variables, only valid for output variables
    INTEGER                               ::  ndim      !< number of NetCDF dimensions
    INTEGER                               ::  nt        !< number of output time steps
    INTEGER                               ::  lod       !< NetCDF attribute indicating the PALM-4U level of detail
    INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) ::  dimids    !< NetCDF IDs of the dimensions
    INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) ::  dimvarids !< IDs of NetCDF dimension variables
    INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) ::  dimlen    !< length of NetCDF dimensions
    CHARACTER(LEN=NF90_MAX_NAME), DIMENSION(NF90_MAX_VAR_DIMS) ::  dimname !< names of NetCDF dimensions
    CHARACTER(LEN=SNAME)                  ::  name                      !< NetCDF short name of the variable
    CHARACTER(LEN=LNAME)                  ::  standard_name             !< NetCDF standard name of the variable
    CHARACTER(LEN=LNAME)                  ::  long_name                 !< NetCDF long name of the variable
    CHARACTER(LEN=LNAME)                  ::  source                    !< NetCDF attribute indicating the data source for the output
    CHARACTER(LEN=SNAME)                  ::  units                     !< NetCDF units of the variable
    CHARACTER(LEN=SNAME)                  ::  kind                      !< Kind of grid
    CHARACTER(LEN=SNAME)                  ::  task                      !< Processing task that generates this variable, e.g. 'interpolate_2d' or 'average profile'
    LOGICAL                               ::  to_be_processed = .FALSE. !< Inifor flag indicating whether variable shall be processed
    LOGICAL                               ::  is_read = .FALSE.         !< Inifor flag indicating whether variable has been read
    LOGICAL                               ::  is_upside_down  = .FALSE. !< Inifor flag indicating whether variable shall be processed
    TYPE(grid_definition), POINTER        ::  grid                      !< Pointer to the corresponding output grid
    TYPE(grid_definition), POINTER        ::  intermediate_grid         !< Pointer to the corresponding intermediate grid
 END TYPE nc_var


 TYPE io_group                                          !< Input/Output group, groups together output variabels that share their input variables. For instance, all boundary surfaces and initialization fields of the potential temperature are base on T and p.
    INTEGER                          ::  nt             !< maximum number of output time steps across all output variables
    INTEGER                          ::  nv             !< number of output variables
    CHARACTER(LEN=SNAME)             ::  kind           !< kind of I/O group
    CHARACTER(LEN=PATH), ALLOCATABLE ::  in_files(:)    !< list of nt input files
    TYPE(nc_var), ALLOCATABLE        ::  out_vars(:)    !< list of output variables 
    TYPE(nc_var), ALLOCATABLE        ::  in_var_list(:) !< list of input variables
    LOGICAL                          ::  to_be_processed = .FALSE. !< Inifor flag indicating whether I/O group shall be processed
    LOGICAL                          ::  is_accumulated = .FALSE.  !< Flag indicating whether this I/O group contains accumulated variables
    LOGICAL                          ::  is_preprocessed = .FALSE. !< Inifor flag indicating whether the I/O group has been preprocessed
 END TYPE io_group 


 TYPE container
   REAL(dp), ALLOCATABLE ::  array(:,:,:)
   LOGICAL               ::  is_preprocessed = .FALSE.
 END TYPE container

 END MODULE types

