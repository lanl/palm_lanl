!> @file src/inifor.f90
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
! $Id: inifor.f90 2718 2018-01-02 08:49:38Z maronga $
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
!> INIFOR is an interpolation tool for generating meteorological initialization
!> and forcing data for the urban climate model PALM-4U. The required 
!> meteorological fields are interpolated from output data of the mesoscale
!> model COSMO-DE. This is the main program file.
!------------------------------------------------------------------------------!
 PROGRAM inifor

    USE control
    USE defs
    USE grid,                                                                  &
        ONLY:  setup_parameters, setup_grids, setup_variable_tables,           &
               setup_io_groups, fini_grids, fini_variables, fini_io_groups,    &
               fini_file_lists, preprocess,                                    &
               output_file, io_group_list, output_var_table,                   &
               cosmo_grid, palm_grid, nx, ny, nz, ug, vg, p0, mode,            &
               imin, imax, jmin, jmax

    USE io
    USE transform,                                                             &
        ONLY:  average_profile, average_2d, interpolate_2d, interpolate_3d
    USE types
    
    IMPLICIT NONE
    
    INTEGER                                 ::  igroup
    INTEGER                                 ::  ivar
    INTEGER                                 ::  iter
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::  output_arr
    TYPE(nc_var), POINTER                   ::  output_var
    TYPE(io_group), POINTER                 ::  group
    TYPE(container), ALLOCATABLE            ::  input_buffer(:)
    
!> \mainpage About INIFOR
!>  ...
!
!------------------------------------------------------------------------------
!- Section 1: Initialization
!------------------------------------------------------------------------------
 CALL run_control('init', 'void')

    ! Initialize INIFOR's parameters from command-line interface and namelists
    CALL setup_parameters()

    ! Initialize all grids, including interpolation neighbours and weights
    CALL setup_grids()
 CALL run_control('time', 'init')

    ! Initialize the netCDF output file and define dimensions
    CALL setup_netcdf_dimensions(output_file, palm_grid)
 CALL run_control('time', 'write')

    ! Set up the tables containing the input and output variables and set
    ! the corresponding netCDF dimensions for each output variable
    CALL setup_variable_tables(mode)
 CALL run_control('time', 'write')

    ! Add the output variables to the netCDF output file
    CALL setup_netcdf_variables(output_file % name, output_var_table)

    CALL setup_io_groups() 
 CALL run_control('time', 'init')

!------------------------------------------------------------------------------
!- Section 2: Main loop
!------------------------------------------------------------------------------

    DO igroup = 1, SIZE(io_group_list)

       group => io_group_list(igroup)
       IF ( group % to_be_processed )  THEN
          
          DO iter = 1, group % nt 

!------------------------------------------------------------------------------
!- Section 2.1: Read and preprocess input data
!------------------------------------------------------------------------------
             CALL read_input_variables(group, iter, input_buffer)
 CALL run_control('time', 'read')

             CALL preprocess(group, input_buffer, cosmo_grid, iter)
 CALL run_control('time', 'comp')

             IF ( .NOT. ALL(input_buffer(:) % is_preprocessed .AND. .TRUE.) )  THEN
                message = "Input buffers for group '" // TRIM(group % kind) // &
                   "' could not be preprocessed sucessfully."
                CALL abort('main loop', message)
             END IF

!------------------------------------------------------------------------------
!- Section 2.2: Interpolate each output variable of the group
!------------------------------------------------------------------------------
             DO ivar = 1, group % nv

                output_var => group % out_vars( ivar )

                IF ( output_var % to_be_processed .AND.                        &
                     iter .LE. output_var % nt )  THEN

                   message = "Processing '" // TRIM(output_var % name) //      &
                             "' (" // TRIM(output_var % kind) //               &
                             "), iteration " // TRIM(str(iter)) //" of " //    &
                             TRIM(str(output_var % nt))
                   CALL report('main loop', message)

                   SELECT CASE( TRIM(output_var % task) )

                   CASE( 'interpolate_2d' ) 
                   
                      SELECT CASE( TRIM(output_var % kind) )
                       
                      CASE( 'init soil' )

                         ALLOCATE( output_arr( 0:output_var % grid % nx,       &
                                               0:output_var % grid % ny,       &
                                               SIZE(output_var % grid % depths) ) )

                      CASE ( 'surface forcing' )

                         ALLOCATE( output_arr( 0:output_var % grid % nx,       &
                                               0:output_var % grid % ny, 1 ) )

                      CASE DEFAULT

                          CALL abort("main loop", 'Not a soil variable')

                      END SELECT
 CALL run_control('time', 'alloc')

                      CALL interpolate_2d(input_buffer(output_var % input_id) % array(:,:,:), &
                              output_arr(:,:,:), output_var % intermediate_grid, output_var)
 CALL run_control('time', 'comp')


                   CASE ( 'interpolate_3d' )

                      ALLOCATE( output_arr( 0:output_var % grid % nx,          &
                                            0:output_var % grid % ny,          &
                                            0:output_var % grid % nz ) )

 CALL run_control('time', 'alloc')
                      CALL interpolate_3d(                                     &
                         input_buffer(output_var % input_id) % array(:,:,:),   &
                         output_arr(:,:,:),                                    &
                         output_var % intermediate_grid,                       &
                         output_var % grid)
 CALL run_control('time', 'comp')

                   CASE ( 'average profile' )

                      ALLOCATE( output_arr( 0:output_var % grid % nx,          &
                                            0:output_var % grid % ny,          &
                                            0:output_var % grid % nz ) )
 CALL run_control('time', 'alloc')
                      

                      CALL average_profile(                                    &
                         input_buffer(output_var % input_id) % array(:,:,:),   &
                         output_arr(:,:,:), imin, imax, jmin, jmax,            &
                         output_var % intermediate_grid,                       &
                         output_var % grid)
 CALL run_control('time', 'comp')

                   CASE ( 'average scalar' )

                      ALLOCATE( output_arr(1,1,1) )
 CALL run_control('time', 'alloc')
                      output_arr(1,1,1) = p0
 CALL run_control('time', 'comp')

                   CASE ( 'profile' )
                      
                      ALLOCATE( output_arr( 1, 1, 0:nz ) )
 CALL run_control('time', 'alloc')

                      SELECT CASE (TRIM(output_var % name))

                      CASE('ls_forcing_ug')
                          output_arr(1, 1, :) = ug

                      CASE('ls_forcing_vg')
                          output_arr(1, 1, :) = vg

                      CASE DEFAULT
                          message = "'" // TRIM(output_var % name) //          &
                             "' is not a valid '" // TRIM(output_var % kind) //&
                             "' variable kind."
                          CALL abort('main loop', message)
                      END SELECT
 CALL run_control('time', 'comp')

                   CASE DEFAULT
                      message = "Processing task '" // TRIM(output_var % task) //&
                               "' not recognized."
                      CALL abort('', message)

                   END SELECT
 CALL run_control('time', 'comp')

!------------------------------------------------------------------------------
!- Section 2.3: Write current time step of current variable
!------------------------------------------------------------------------------
                   message = "Writing variable '" // TRIM(output_var%name) // "'."
                   CALL report('main loop', message)
                   CALL update_output(output_var, output_arr, iter, output_file)
 CALL run_control('time', 'write')

                   DEALLOCATE(output_arr)
 CALL run_control('time', 'alloc')

                END IF

             END DO ! ouput variables

             IF ( group % kind == 'running average' .OR. &
                  group % kind == 'accumulated' )  THEN
                ! Keep input buffer around for averaged (radiation) and 
                ! accumulated COSMO-DE quantities (precipitation).
             ELSE
                CALL report('main loop', 'Deallocating input buffer')
                DEALLOCATE(input_buffer)
             END IF
 CALL run_control('time', 'alloc')

          END DO ! time steps / input files 

          IF (ALLOCATED(input_buffer))  THEN
             CALL report('main loop', 'Deallocating input buffer')
             DEALLOCATE(input_buffer)
          END IF
 CALL run_control('time', 'alloc')

       ELSE 

          message = "Skipping IO group '" // TRIM(group % kind) // "'"
          IF ( ALLOCATED(group % in_var_list) )  THEN
              message = TRIM(message) // " with input variable '" //           &
              TRIM(group % in_var_list(1) % name) // "'."
          END IF

          CALL report('main loop', message)

       END IF ! IO group % to_be_processed

    END DO ! IO groups

!------------------------------------------------------------------------------
!- Section 3: Clean up.
!------------------------------------------------------------------------------
    CALL fini_file_lists()
    CALL fini_io_groups()
    CALL fini_variables()
    !CALL fini_grids()
 CALL run_control('time', 'alloc')
 CALL run_control('report', 'void')

    message = "Finished writing forcing file '" // TRIM(output_file % name) // &
              "' successfully."
    CALL report('main loop', message)


 END PROGRAM inifor
