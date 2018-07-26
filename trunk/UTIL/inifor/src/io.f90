!> @file src/io.f90
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
! $Id: io.f90 2718 2018-01-02 08:49:38Z maronga $
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
!> The io module contains the functions needed to read and write netCDF data in
!> INIFOR.
!------------------------------------------------------------------------------!
 MODULE io

    USE control
    USE defs,                                                                  &
        ONLY:  DATE, SNAME, PATH, PI, dp, TO_RADIANS, TO_DEGREES, VERSION
    USE netcdf
    USE types
    USE util,                                                                  &
        ONLY:  reverse, str

    IMPLICIT NONE

 CONTAINS

    SUBROUTINE netcdf_define_variable(var, ncid)

        TYPE(nc_var), INTENT(INOUT) ::  var
        INTEGER, INTENT(IN)         ::  ncid

        CALL check(nf90_def_var(ncid, var % name, NF90_FLOAT,       var % dimids(1:var % ndim), var % varid))
        CALL check(nf90_put_att(ncid, var % varid, "standard_name", var % standard_name))
        CALL check(nf90_put_att(ncid, var % varid, "long_name",     var % long_name))
        CALL check(nf90_put_att(ncid, var % varid, "units",         var % units))
        CALL check(nf90_put_att(ncid, var % varid, "lod",           var % lod))
        CALL check(nf90_put_att(ncid, var % varid, "source",        var % source))
        CALL check(nf90_put_att(ncid, var % varid, "_FillValue",    NF90_FILL_REAL))

    END SUBROUTINE netcdf_define_variable
    

    SUBROUTINE netcdf_get_dimensions(var, ncid)

        TYPE(nc_var), INTENT(INOUT) ::  var
        INTEGER, INTENT(IN)         ::  ncid
        INTEGER                     ::  i
        CHARACTER(SNAME)            ::  null

        ! Remember dimension lenghts for NetCDF writing routine
        DO i = 1, var % ndim
           CALL check(nf90_inquire_dimension(ncid, var % dimids(i), &
                                             name = null, &
                                             len  = var % dimlen(i)  ) )
        END DO

    END SUBROUTINE netcdf_get_dimensions


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine initializes Inifor. This includes parsing command-line
!> arguments, setting the names of the input and output file names as well as
!> the name of the input namelist and, subsequently, reading in and setting grid
!> parameters for the PALM-4U computational grid.
!------------------------------------------------------------------------------!
    SUBROUTINE parse_command_line_arguments( start_date, hhl_file,             &
       soiltyp_file, static_driver_file, input_path, output_file,              &
       namelist_file, ug, vg, p0, z0, mode )

       CHARACTER(LEN=PATH), INTENT(INOUT)  ::  hhl_file, soiltyp_file,         &
           static_driver_file, input_path, output_file, namelist_file
       CHARACTER(LEN=SNAME), INTENT(INOUT) ::  mode
       REAL(dp), INTENT(INOUT)             ::  ug, vg, p0, z0
       CHARACTER(LEN=DATE), INTENT(INOUT)  ::  start_date

       CHARACTER(LEN=PATH)     ::  option, arg
       INTEGER                 ::  arg_count, i

       arg_count = COMMAND_ARGUMENT_COUNT()
       IF (arg_count .GT. 0)  THEN

          ! Every option should have an argument.
          IF ( MOD(arg_count, 2) .NE. 0 )  THEN
             message = "Syntax error in command line."
             CALL abort('parse_command_line_arguments', message)
          END IF
          
          message = "The -clon and -clat command line options are depricated. " // &
             "Please remove them form your inifor command and specify the " // &
             "location of the PALM-4U origin either" // NEW_LINE(' ') // &
             "   - by setting the namelist parameters 'origin_lon' and 'origin_lat, or'" // NEW_LINE(' ') // &
             "   - by providing a static driver netCDF file via the -static command-line option."

          ! Loop through option/argument pairs.
          DO i = 1, arg_count, 2

             CALL GET_COMMAND_ARGUMENT( i, option )
             CALL GET_COMMAND_ARGUMENT( i+1, arg )

             SELECT CASE( TRIM(option) )

             CASE( '-date' )
                start_date = TRIM(arg)

             ! Elevation of the PALM-4U domain above sea level
             CASE( '-z0' )
                READ(arg, *) z0

             ! surface pressure, at z0
             CASE( '-p0' )
                READ(arg, *) p0

             ! surface pressure, at z0
             CASE( '-ug' )
                READ(arg, *) ug

             ! surface pressure, at z0
             CASE( '-vg' )
                READ(arg, *) vg

             ! Domain centre geographical longitude
             CASE( '-clon' )
                CALL abort('parse_command_line_arguments', message)          
                !READ(arg, *) lambda_cg
                !lambda_cg = lambda_cg * TO_RADIANS

             ! Domain centre geographical latitude
             CASE( '-clat' )
                CALL abort('parse_command_line_arguments', message)          
                !READ(arg, *) phi_cg
                !phi_cg = phi_cg * TO_RADIANS

             CASE( '-path' )
                 input_path = TRIM(arg)

             CASE( '-hhl' )
                 hhl_file = TRIM(arg)

             CASE( '-static' )
                 static_driver_file = TRIM(arg)

             CASE( '-soil' )
                 soiltyp_file = TRIM(arg)

             CASE( '-o' )
                output_file = TRIM(arg)

             CASE( '-n' )
                namelist_file = TRIM(arg)

             ! Initialization mode: 'profile' / 'volume'
             CASE( '-mode' )
                mode = TRIM(arg)

                SELECT CASE( TRIM(mode) )

                CASE( 'profile' )

                CASE DEFAULT
                   message = "Mode '" // TRIM(mode) // "' is not supported. " //&
                             "Currently, '-mode profile' is the only supported option. " //&
                             "Select this one or omit the -mode option entirely."
                   CALL abort( 'parse_command_line_arguments', message ) 
                END SELECT

             CASE DEFAULT
                message = "unknown option '" // TRIM(option(2:)) // "'."
                CALL abort('parse_command_line_arguments', message)

             END SELECT

          END DO

       ELSE
            
          message = "No arguments present, using default input and output files"
          CALL report('parse_command_line_arguments', message)

       END IF

   END SUBROUTINE parse_command_line_arguments


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine initializes the Inifor output file, i.e. the PALM-4U
!> initializing and forcing data as a NetCDF file.
!> 
!> Besides writing metadata, such as global attributes, coordinates, variables,
!> in the NetCDF file, various NetCDF IDs are saved for later, when Inifor
!> writes the actual data.
!------------------------------------------------------------------------------!
    SUBROUTINE setup_netcdf_dimensions(output_file, palm_grid)

       TYPE(nc_file), INTENT(INOUT)      ::  output_file
       TYPE(grid_definition), INTENT(IN) ::  palm_grid

       CHARACTER (LEN=SNAME) ::  date
       INTEGER               ::  ncid, nx, ny, nz, nt, dimids(3), dimvarids(3)
       REAL(dp)              ::  z0

       ! Create the NetCDF file. NF90_CLOBBER selects overwrite mode.
       CALL check(nf90_create(TRIM(output_file % name), OR(NF90_CLOBBER, NF90_HDF5), ncid))

!
!------------------------------------------------------------------------------
!- Section 1: Write global NetCDF attributes
!------------------------------------------------------------------------------
       CALL date_and_time(date)
       CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'title',          'PALM input file for scenario ...'))
       CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'institution',    'Deutscher Wetterdienst, Offenbach'))
       CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'author',         'Eckhard Kadasch, eckhard.kadasch@dwd.de'))
       CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'history',        'Created on '//date))
       CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'references',     '--'))
       CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'comment',        '--'))
       CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'origin_lat',     '--'))
       CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'origin_lon',     '--'))
       CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'inifor_version', TRIM(VERSION)))
       CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'palm_version',   '--'))

!
!------------------------------------------------------------------------------
!- Section 2: Define NetCDF dimensions and coordinates
!------------------------------------------------------------------------------
       nt = SIZE(output_file % time)
       nx = palm_grid % nx
       ny = palm_grid % ny
       nz = palm_grid % nz
       z0 = palm_grid % z0

!
!------------------------------------------------------------------------------
!- Section 2a: Define dimensions for cell centers (scalars in soil and atmosph.)
!------------------------------------------------------------------------------
       dimids = (/0, 0, 0/) ! reset dimids
          CALL check( nf90_def_dim(ncid, "x", nx+1, dimids(1)) )
          CALL check( nf90_def_dim(ncid, "y", ny+1, dimids(2)) )
          CALL check( nf90_def_dim(ncid, "z", nz+1, dimids(3)) )
          output_file % dimids_scl = dimids ! save dimids for later

       dimvarids = (/0, 0, 0/) ! reset dimvarids
          CALL check(nf90_def_var(ncid, "x", NF90_FLOAT, dimids(1), dimvarids(1)))
          CALL check(nf90_put_att(ncid, dimvarids(1), "standard_name", "x coordinate of cell centers"))
          CALL check(nf90_put_att(ncid, dimvarids(1), "units", "m"))

          CALL check(nf90_def_var(ncid, "y", NF90_FLOAT, dimids(2), dimvarids(2)))
          CALL check(nf90_put_att(ncid, dimvarids(2), "standard_name", "y coordinate of cell centers"))
          CALL check(nf90_put_att(ncid, dimvarids(2), "units", "m"))

          CALL check(nf90_def_var(ncid, "z", NF90_FLOAT, dimids(3), dimvarids(3)))
          CALL check(nf90_put_att(ncid, dimvarids(3), "standard_name", "z coordinate of cell centers"))
          CALL check(nf90_put_att(ncid, dimvarids(3), "units", "m"))
       output_file % dimvarids_scl = dimvarids ! save dimvarids for later

       ! overwrite third dimid with the one of depth
       CALL check(nf90_def_dim(ncid, "depth", SIZE(palm_grid % depths), dimids(3)) )
       output_file % dimids_soil = dimids ! save dimids for later

       ! overwrite third dimvarid with the one of depth
       CALL check(nf90_def_var(ncid, "depth", NF90_FLOAT, output_file % dimids_soil(3), dimvarids(3)))
       CALL check(nf90_put_att(ncid, dimvarids(3), "standard_name", "depth_below_land"))
       CALL check(nf90_put_att(ncid, dimvarids(3), "positive", "down"))
       CALL check(nf90_put_att(ncid, dimvarids(3), "units", "m"))
       output_file % dimvarids_soil = dimvarids ! save dimvarids for later
!
!------------------------------------------------------------------------------
!- Section 2b: Define dimensions for cell faces/velocities
!------------------------------------------------------------------------------
       dimids = (/0, 0, 0/) ! reset dimids
          CALL check(nf90_def_dim(ncid, "xu", nx, dimids(1)) )
          CALL check(nf90_def_dim(ncid, "yv", ny, dimids(2)) )
          CALL check(nf90_def_dim(ncid, "zw", nz, dimids(3)) )
       output_file % dimids_vel = dimids ! save dimids for later

       dimvarids = (/0, 0, 0/) ! reset dimvarids
          CALL check(nf90_def_var(ncid, "xu", NF90_FLOAT, dimids(1), dimvarids(1)))
          CALL check(nf90_put_att(ncid, dimvarids(1), "standard_name", "x coordinate of cell faces"))
          CALL check(nf90_put_att(ncid, dimvarids(1), "units", "m"))

          CALL check(nf90_def_var(ncid, "yv", NF90_FLOAT, dimids(2), dimvarids(2)))
          CALL check(nf90_put_att(ncid, dimvarids(2), "standard_name", "y coordinate of cell faces"))
          CALL check(nf90_put_att(ncid, dimvarids(2), "units", "m"))

          CALL check(nf90_def_var(ncid, "zw", NF90_FLOAT, dimids(3), dimvarids(3)))
          CALL check(nf90_put_att(ncid, dimvarids(3), "standard_name", "z coordinate of cell faces"))
          CALL check(nf90_put_att(ncid, dimvarids(3), "units", "m"))
       output_file % dimvarids_vel = dimvarids ! save dimvarids for later

!
!------------------------------------------------------------------------------
!- Section 2c: Define time dimension
!------------------------------------------------------------------------------
       CALL check(nf90_def_dim(ncid, "time", nt, output_file % dimid_time) )
       CALL check(nf90_def_var(ncid, "time", NF90_FLOAT, &
                                             output_file % dimid_time, &
                                             output_file % dimvarid_time))
       CALL check(nf90_put_att(ncid, output_file % dimvarid_time, "standard_name", "time"))
       CALL check(nf90_put_att(ncid, output_file % dimvarid_time, "long_name", "time"))
       CALL check(nf90_put_att(ncid, output_file % dimvarid_time, "units", "seconds since..."))

       CALL check(nf90_enddef(ncid))

!
!------------------------------------------------------------------------------
!- Section 3: Write grid coordinates
!------------------------------------------------------------------------------
       CALL check(nf90_put_var(ncid, output_file % dimvarids_scl(1), palm_grid%x))
       CALL check(nf90_put_var(ncid, output_file % dimvarids_scl(2), palm_grid%y))
       CALL check(nf90_put_var(ncid, output_file % dimvarids_scl(3), palm_grid%z))

       CALL check(nf90_put_var(ncid, output_file % dimvarids_vel(1), palm_grid%xu))
       CALL check(nf90_put_var(ncid, output_file % dimvarids_vel(2), palm_grid%yv))
       CALL check(nf90_put_var(ncid, output_file % dimvarids_vel(3), palm_grid%zw))

       ! TODO Read in soil depths from input file before this.
       CALL check(nf90_put_var(ncid, output_file % dimvarids_soil(3), palm_grid%depths))

       ! Write time vector
       CALL check(nf90_put_var(ncid, output_file % dimvarid_time, output_file % time))

       ! Close the file
       CALL check(nf90_close(ncid))

    END SUBROUTINE setup_netcdf_dimensions


    SUBROUTINE setup_netcdf_variables(filename, output_variable_table)

       CHARACTER (LEN=*), INTENT(IN)        ::  filename
       TYPE(nc_var), INTENT(INOUT), TARGET  ::  output_variable_table(:)
       TYPE(nc_var), POINTER                ::  var
       INTEGER                              ::  i, ncid

       message = "Initializing PALM-4U forcing file '" // TRIM(filename) // "'."
       CALL report('setup_netcdf_variables', message)

       CALL check(nf90_open(TRIM(filename), NF90_WRITE, ncid))
       CALL check(nf90_redef(ncid))

       DO i = 1, SIZE(output_variable_table)

          var => output_variable_table(i)

          IF ( var % to_be_processed )  THEN
             message = "Defining variable #" // TRIM(str(i)) // " '" // TRIM(var%name) // "'."
             CALL report('setup_netcdf_variables', message)

             CALL netcdf_define_variable(var, ncid)
             CALL netcdf_get_dimensions(var, ncid)
          END IF
           
       END DO

       CALL check(nf90_enddef(ncid))
       CALL check(nf90_close(ncid))

       message = "Forcing file '" // TRIM(filename) // "' initialized successfully."
       CALL report('setup_netcdf_variables', message)

    END SUBROUTINE setup_netcdf_variables


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine reads and returns all input variables of the given IO group
!> It accomodates the data by allocating a container variable that represents a
!> list of arrays of the same length as the groups variable list. (This list
!> will typically contain one or two items.) After the container, its members
!> are allocated one by one with the appropriate, possibly different,
!> dimensions.
!>
!> The 'group' is an INTENT(INOUT) variable so that 'get_netcdf_variable()' can
!> record netCDF IDs in the 'in_var_list()' member variable.
!------------------------------------------------------------------------------!
    SUBROUTINE read_input_variables(group, iter, buffer)
       TYPE(io_group), INTENT(INOUT), TARGET       ::  group
       INTEGER, INTENT(IN)                         ::  iter
       TYPE(container), ALLOCATABLE, INTENT(INOUT) ::  buffer(:)
       INTEGER                                     ::  hour, buf_id
       TYPE(nc_var), POINTER                       ::  input_var
       CHARACTER(LEN=PATH), POINTER                ::  input_file
       INTEGER                                     ::  ivar, nbuffers

       message = "Reading data for I/O group '" // TRIM(group % in_var_list(1) % name) // "'."
       CALL report('read_input_variables', message)

       input_file => group % in_files(iter)

!
!------------------------------------------------------------------------------
!- Section 1: Load input buffers for accumulated variables
!------------------------------------------------------------------------------
       IF (group % kind == 'running average' .OR.                              &
           group % kind == 'accumulated')  THEN ! radiation budgets, precipitation

          IF (SIZE(group % in_var_list) .GT. 1 ) THEN
             message = "I/O groups may not contain more than one " // & 
                       "accumulated variable. Group '" // TRIM(group % kind) //&
                       "' contains " //                                        &
                       TRIM( str(SIZE(group % in_var_list)) ) // "."
             CALL abort('read_input_variables | accumulation', message)
          END IF

          ! use two buffer arrays
          nbuffers = 2
          IF ( .NOT. ALLOCATED( buffer ) )  ALLOCATE( buffer(nbuffers) )

          ! chose correct buffer array
          hour = iter - 1! hour of the day
          buf_id = select_buffer(hour)

 CALL run_control('time', 'read')
          IF ( ALLOCATED(buffer(buf_id) % array) )  DEALLOCATE(buffer(buf_id) % array)
 CALL run_control('time', 'alloc')

          input_var => group % in_var_list(1)
          buffer(buf_id) % array = get_netcdf_variable( input_file, input_var )  
          CALL report('read_input_variables', "Read accumulated " // TRIM(group % in_var_list(1) % name)) 

          IF ( input_var % is_upside_down )  CALL reverse(buffer(buf_id) % array)
 CALL run_control('time', 'comp')
          
!------------------------------------------------------------------------------
!- Section 2: Load input buffers for normal I/O groups
!------------------------------------------------------------------------------
       ELSE

          nbuffers = SIZE( group % in_var_list )
          ALLOCATE( buffer(nbuffers) )
 CALL run_control('time', 'alloc')
          
          DO ivar = 1, nbuffers

             input_var => group % in_var_list(ivar)

             ! Check wheather P or PP is present in input file
             IF (input_var % name == 'P')  THEN
                input_var % name = TRIM( get_pressure_var(input_file) )
 CALL run_control('time', 'read')
             END IF

             buffer(ivar) % array = get_netcdf_variable( input_file, input_var )  

             IF ( input_var % is_upside_down )  CALL reverse(buffer(ivar) % array)
 CALL run_control('time', 'comp')

          END DO
       END IF

    END SUBROUTINE read_input_variables


    INTEGER FUNCTION select_buffer(hour)
       INTEGER, INTENT(IN) ::  hour
       INTEGER             ::  step

       select_buffer = 0
       step = MODULO(hour, 3) + 1

       SELECT CASE(step)
       CASE(1, 3)
           select_buffer = 1
       CASE(2)
           select_buffer = 2
       CASE DEFAULT
           message = "Invalid step '" // TRIM(str(step))
           CALL abort('select_buffer', message)
       END SELECT
    END FUNCTION select_buffer


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks if the input_file contains the absolute pressure, 'P', or the pressure
!> perturbation, 'PP', and returns the appropriate string.
!------------------------------------------------------------------------------!
    CHARACTER(LEN=2) FUNCTION get_pressure_var(input_file) RESULT(var)
       CHARACTER(LEN=*) ::  input_file
       INTEGER          ::  ncid, varid

       CALL check(nf90_open( TRIM(input_file), NF90_NOWRITE, ncid ))
       IF ( nf90_inq_varid( ncid, 'P', varid ) .EQ. NF90_NOERR )  THEN

          var = 'P'

       ELSE IF ( nf90_inq_varid( ncid, 'PP', varid ) .EQ. NF90_NOERR )  THEN

          var = 'PP'
          CALL report('get_pressure_var', 'Using PP instead of P')

       ELSE

          message = "Failed to read '" // TRIM(var) // &
                    "' from file '" // TRIM(input_file) // "'."
          CALL abort('get_pressure_var', message)

       END IF

       CALL check(nf90_close(ncid))

    END FUNCTION get_pressure_var


    FUNCTION get_netcdf_attribute(filename, attribute) RESULT(attribute_value)

       CHARACTER(LEN=*), INTENT(IN) ::  filename, attribute
       REAL(dp)                     ::  attribute_value

       INTEGER                      :: ncid

       IF ( nf90_open( TRIM(filename), NF90_NOWRITE, ncid ) == NF90_NOERR )  THEN

          CALL check(nf90_get_att(ncid, NF90_GLOBAL, TRIM(attribute), attribute_value))

       ELSE

          message = "Failed to read '" // TRIM(attribute) // &
                    "' from file '" // TRIM(filename) // "'."
          CALL abort('get_netcdf_attribute', message)

       END IF

    END FUNCTION get_netcdf_attribute



    FUNCTION get_netcdf_variable(in_file, in_var) RESULT(buffer)

       CHARACTER(LEN=PATH), INTENT(IN)      ::  in_file
       TYPE(nc_var), INTENT(INOUT)          ::  in_var
       REAL(dp), ALLOCATABLE                ::  buffer(:,:,:)
       INTEGER                              ::  i, ncid, start(3)


       ! Read in_var NetCDF attributes
       IF ( nf90_open( TRIM(in_file), NF90_NOWRITE, ncid ) .EQ. NF90_NOERR .AND. &
            nf90_inq_varid( ncid, in_var % name, in_var % varid ) .EQ. NF90_NOERR )  THEN

          CALL check(nf90_get_att(ncid, in_var % varid, "long_name", in_var % long_name))
          CALL check(nf90_get_att(ncid, in_var % varid, "units", in_var % units))

          ! Read in_var NetCDF dimensions
          CALL check(nf90_inquire_variable( ncid, in_var % varid,              &
                                            ndims  = in_var % ndim,            &
                                            dimids = in_var % dimids ))
                                        
          DO i = 1, in_var % ndim
             CALL check(nf90_inquire_dimension( ncid, in_var % dimids(i),      &
                                                name = in_var % dimname(i),    &
                                                len  = in_var % dimlen(i) ))
          END DO

          start = (/ 1, 1, 1 /)
          IF ( TRIM(in_var % name) .EQ. 'T_SO' )  THEN
             ! Skip depth = 0.0 for T_SO and reduce number of depths from 9 to 8
             in_var % dimlen(3) = in_var % dimlen(3) - 1

             ! Start reading from second level, e.g. depth = 0.005 instead of 0.0
             start(3) = 2
          END IF

          SELECT CASE(in_var % ndim)

          CASE (2)

             ALLOCATE( buffer( in_var % dimlen(1),                             &
                               in_var % dimlen(2),                             &
                               1 ) )

          CASE (3)

             ALLOCATE( buffer( in_var % dimlen(1),                             &
                               in_var % dimlen(2),                             &
                               in_var % dimlen(3) ) )
          CASE (4)

             ALLOCATE( buffer( in_var % dimlen(1),                             &
                               in_var % dimlen(2),                             &
                               in_var % dimlen(3) ) )
          CASE DEFAULT

             message = "Failed reading NetCDF variable " //                    &
                TRIM(in_var % name) // " with " // TRIM(str(in_var%ndim)) //   &
                " dimensions because only two- and and three-dimensional" //   &
                " variables are supported."
             CALL abort('get_netcdf_variable', message) 

          END SELECT
 CALL run_control('time', 'alloc')
          
          ! TODO: Check for matching dimensions of buffer and var
          CALL check(nf90_get_var( ncid, in_var % varid, buffer,               &
                                   start = start,                              &
                                   count = in_var % dimlen(1:3) ) )

 CALL run_control('time', 'read')
       ELSE

          message = "Failed to read '" // TRIM(in_var % name) // &
             "' from file '" // TRIM(in_file) // "'."
          CALL report('get_netcdf_variable', message)

       END IF

       CALL check(nf90_close(ncid))

 CALL run_control('time', 'read')

    END FUNCTION get_netcdf_variable


    SUBROUTINE update_output(var, array, iter, output_file)
       TYPE(nc_var), INTENT(IN)  ::  var
       REAL(dp), INTENT(IN)      ::  array(:,:,:)
       INTEGER, INTENT(IN)       ::  iter
       TYPE(nc_file), INTENT(IN) ::  output_file

       INTEGER ::  ncid, ndim, start(4), count(4)
       LOGICAL ::  var_is_time_dependent

       var_is_time_dependent = (                                               &
          var % dimids( var % ndim ) == output_file % dimid_time               &
       )

       ! Skip time dimension for output
       IF ( var_is_time_dependent )  THEN
          ndim = var % ndim - 1
       ELSE 
          ndim = var % ndim
       END IF

       start(:)      = (/1,1,1,1/)
       start(ndim+1) = iter
       count(1:ndim) = var%dimlen(1:ndim)

       CALL check(nf90_open(output_file % name, NF90_WRITE, ncid))

       ! Reduce dimension of output array according to variable kind
       SELECT CASE (TRIM(var % kind))
       
       CASE ( 'init scalar profile', 'init u profile', 'init v profile',       &
              'init w profile' )

          CALL check(nf90_put_var( ncid, var%varid, array(1,1,:) ) )

       CASE ( 'init soil', 'init scalar', 'init u', 'init v', 'init w' )

          CALL check(nf90_put_var( ncid, var%varid, array(:,:,:) ) )

       CASE( 'left scalar', 'right scalar', 'left w', 'right w' ) 

          CALL check(nf90_put_var( ncid, var%varid, array(1,:,:),              &
                                   start=start(1:ndim+1),                      &
                                   count=count(1:ndim) ) )
          

          IF (.NOT. SIZE(array, 2) .EQ. var % dimlen(1))  THEN
             PRINT *, "inifor: update_output: Dimension ", 1, " of variable ", &
                 TRIM(var % name), " (", var % dimlen(1),                      &
                 ") does not match the dimension of the output array (",       &
                 SIZE(array, 2), ")."
             STOP
          END IF
          

       CASE( 'north scalar', 'south scalar', 'north w', 'south w' )

          CALL check(nf90_put_var( ncid, var%varid, array(:,1,:),              &
                                   start=start(1:ndim+1),                      &
                                   count=count(1:ndim) ) )
          

       CASE( 'surface forcing', 'top scalar', 'top w' )

          CALL check(nf90_put_var( ncid, var%varid, array(:,:,1),              &
                                   start=start(1:ndim+1),                      &
                                   count=count(1:ndim) ) )
          
       CASE ( 'left u', 'right u', 'left v', 'right v' )

          CALL check(nf90_put_var( ncid, var%varid, array(1,:,:),              &
                                   start=start(1:ndim+1),                      &
                                   count=count(1:ndim) ) )

       CASE ( 'north u', 'south u', 'north v', 'south v' )

          CALL check(nf90_put_var( ncid, var%varid, array(:,1,:),              &
                                   start=start(1:ndim+1),                      &
                                   count=count(1:ndim) ) )

       CASE ( 'top u', 'top v' )

          CALL check(nf90_put_var( ncid, var%varid, array(:,:,1),              &
                                   start=start(1:ndim+1),                      &
                                   count=count(1:ndim) ) )

       CASE ( 'time series' )

          CALL check(nf90_put_var( ncid, var%varid, array(1,1,1),              &
                                   start=start(1:ndim+1) ) )

       CASE ( 'profile' )

          CALL check(nf90_put_var( ncid, var%varid, array(1,1,:),              &
                                   start=start(1:ndim+1),                      &
                                   count=count(1:ndim) ) )

       CASE DEFAULT

           message = "Variable kind '" // TRIM(var % kind) //                  &
                    "' not recognized."
           CALL abort('update_output', message)

       END SELECT

       CALL check(nf90_close(ncid))

    END SUBROUTINE update_output


    SUBROUTINE write_netcdf_variable_2d(var, iter, output_file, buffer)
       TYPE(nc_var), INTENT(IN)          ::  var
       INTEGER, INTENT(IN)               ::  iter
       TYPE(nc_file), INTENT(IN)         ::  output_file
       REAL(dp), INTENT(IN)              ::  buffer(:,:,:)

       INTEGER ::  ncid, ndim_out, start(4), count(4)
       LOGICAL ::  last_dimension_is_time

       ndim_out = var % ndim

       last_dimension_is_time = var % dimids( var % ndim ) == output_file % dimid_time
       IF ( last_dimension_is_time )  THEN
          ndim_out = ndim_out - 1
       END IF

       start(:)      = (/1,1,1,iter/)
       count(1:ndim_out) = var%dimlen(1:ndim_out)

       CALL check(nf90_open(output_file % name, NF90_WRITE, ncid))

       IF (TRIM(var % kind) .EQ. 'left/right scalar')  THEN

          CALL check(nf90_put_var( ncid, var%varid, buffer(1,:,:) ) )

       ELSE IF (TRIM(var % kind) .EQ. 'north/south scalar')  THEN

          CALL check(nf90_put_var( ncid, var%varid, buffer(:,1,:) ) )

       ELSE IF (TRIM(var % kind) .EQ. 'top scalar')  THEN

          CALL check(nf90_put_var( ncid, var%varid, buffer(:,:,1) ) )
       ELSE

          CALL check(nf90_put_var( ncid, var%varid, buffer ) )

       END IF
       CALL check(nf90_close(ncid))

    END SUBROUTINE write_netcdf_variable_2d


    SUBROUTINE check(status)

       INTEGER, INTENT(IN) ::  status

       IF (status /= nf90_noerr)  THEN
          message = "NetCDF API call failed with error: " //                     &
                    TRIM( nf90_strerror(status) )
          CALL abort('io.check', message) 
       END IF

    END SUBROUTINE check

 END MODULE io
