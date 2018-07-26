 PROGRAM read_palm_netcdf_data

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
! Copyright 1997-2018  Leibniz Universitaet Hannover
!------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: read_palm_netcdf_data.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 1310 2014-03-14 08:01:56Z raasch
! update of GPL copyright
!
! 1046 2012-11-09 14:38:45Z maronga
! code put under GPL (PALM 3.9)
!
! Description:
! ------------                     
! This is an example program for reading PALM 2d/3d NetCDF datasets
!
! The user has to add his own code for further storing and analyzing of
! these data!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The NetCDF include file and library has to be given with the respective
! compiler options. Please find out the respective paths on your system and
! set them appropriately.
!
! Here are some examples how this routine should be compiled:
!
! decalpha:
!    f95 -fast -r8 -I/usr/local/netcdf-3.5.1/include
!    -L/usr/local/netcdf-3.5.1/lib -lnetcdf
! IBM-Regatta:
!    xlf95 -qrealsize=8 -q64 -qmaxmem=-1 -Q
!    -I /aws/dataformats/netcdf-3.6.0-p1/64-32/include
!    -L/aws/dataformats/netcdf-3.6.0-p1/64-32/lib -lnetcdf -O3
! IBM-Regatta KISTI:
!    xlf95 -qrealsize=8 -q64 -qmaxmem=-1 -Q
!    -I /applic/netcdf64/src/f90
!    -L/applic/lib/NETCDF64 -lnetcdf -O3
! IBM-Regatta Yonsei (gfdl5):
!    xlf95 -qrealsize=8 -q64 -qmaxmem=-1 -Q
!    -I /usr1/users/raasch/pub/netcdf-3.6.0-p1/include
!    -L/usr1/users/raasch/pub/netcdf-3.6.0-p1/lib -lnetcdf -O3
! IMUK:
!    ifort read_palm...f90 -o read_palm...x
!    -I /muksoft/packages/netcdf/linux/include -axW -r8 -nbs
!    -Vaxlib -L /muksoft/packages/netcdf/linux/lib -lnetcdf
! NEC-SX6:
!    sxf90 read_palm...f90 -o read_palm...x 
!    -I /pool/SX-6/netcdf/netcdf-3.6.0-p1/include  -C hopt -Wf '-A idbl4'
!    -L/pool/SX-6/netcdf/netcdf-3.6.0-p1/lib -lnetcdf
!------------------------------------------------------------------------------!

    USE netcdf

    IMPLICIT NONE

!
!-- Local variables
    CHARACTER (LEN=10)   ::  dimname(4), var_name
    CHARACTER (LEN=40)   ::  filename

    CHARACTER (LEN=2000) ::  title, var_list

    INTEGER ::  i, j, k, nc_stat, pos, time_step

    INTEGER ::  current_level, current_var, id_set, id_var_time, num_var

    INTEGER, DIMENSION(4) ::  id_dims, id_dims_loc, levels

    INTEGER, DIMENSION(1000) ::  id_var

    REAL ::  time(1)

    REAL, DIMENSION(:,:,:), ALLOCATABLE ::  data_array


    PRINT*, '*** Please type NetCDF filename to be read:'
    READ*, filename

    nc_stat = NF90_OPEN( filename, NF90_WRITE, id_set )
    IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 1 )

!
!-- Get the run description header and output
    title = ' '
    nc_stat = NF90_GET_ATT( id_set, NF90_GLOBAL, 'title', title )
    IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 2 )
    WRITE (*,'(/A/A)')  '*** file created by:', TRIM( title )

!
!-- Get the list of variables (order of variables corresponds with the
!-- order of data on the binary file)
    var_list = ' '    ! GET_ATT does not assign trailing blanks
    nc_stat = NF90_GET_ATT( id_set, NF90_GLOBAL, 'VAR_LIST', var_list )
    IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 3 )

!
!-- Inquire id of the time coordinate variable
    nc_stat = NF90_INQ_VARID( id_set, 'time', id_var_time )
    IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 4 )

!
!-- Count number of variables; there is one more semicolon in the
!-- string than variable names
    num_var = -1
    DO  i = 1, LEN( var_list )
       IF ( var_list(i:i) == ';' )  num_var = num_var + 1
    ENDDO
    WRITE (*,'(/A,I3,A/)')  '*** file contains ', num_var, ' variable(s)'


    pos = INDEX( var_list, ';' )
!
!-- Loop over all variables
    DO  i = 1, num_var

!
!--    Extract variable name from list
       var_list = var_list(pos+1:)
       pos = INDEX( var_list, ';' )
       var_name = var_list(1:pos-1)

!
!--    Get variable ID from name
       nc_stat = NF90_INQ_VARID( id_set, TRIM( var_name ), id_var(i) )
       IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 5 )

!
!--    Inquire the dimension IDs
       nc_stat = NF90_INQUIRE_VARIABLE( id_set, id_var(i), &
                                        dimids = id_dims_loc )
       IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 6 )
       id_dims = id_dims_loc

!
!--    Get number of x/y/z/time levels(gridpoints) for that variable
       DO  j = 1, 4
          nc_stat = NF90_INQUIRE_DIMENSION( id_set, id_dims(j),&
                                            dimname(j), levels(j) )
          IF ( nc_stat /= NF90_NOERR ) CALL handle_netcdf_error( 7 )
       ENDDO

       WRITE (*,100)  '*** reading variable "', TRIM(var_name),         &
                      '", dimensioned as', TRIM(var_name), levels(1)-1, &
                      levels(2)-1, levels(3)-1
100    FORMAT (A,A,A/4X,A,'(0:',I4,',0:',I4,',0:',I4,')   (x/y/z)'/)

!
!--    Allocate the data array to be read
       ALLOCATE( data_array(0:levels(1)-1,0:levels(2)-1,0:levels(3)-1) )

!
!--    Read the data from file for each timestep
       DO  j = 1, levels(4)

!
!--        Get the time of the current timelevel and output
           nc_stat = NF90_GET_VAR( id_set, id_var_time, time, start = (/ j /), &
                                   count = (/ 1 /) )

           IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 7+i )

           WRITE (*,'(A,I3,A,F8.1,A)')  '    reading timelevel ', i, &
                                        '    time = ', time(1), ' s'     

           nc_stat = NF90_GET_VAR( id_set, id_var(i),                    &
                                   data_array, start = (/ 1, 1, 1, j /), &
                             count = (/ levels(1), levels(2), levels(3), 1 /) )

           IF ( nc_stat /= NF90_NOERR )  &
                                       CALL handle_netcdf_error( 7+levels(4)+i )
!
!--        ADD YOUR OWN CODE FOR FURTHER STORING OF THESE DATA HERE
!--        --------------------------------------------------------


       ENDDO

       WRITE (*,'(/)')

       DEALLOCATE( data_array )

    ENDDO



 CONTAINS


    SUBROUTINE handle_netcdf_error( errno )
!
!--    Prints out a text message corresponding to the current NetCDF status

       IMPLICIT NONE

       INTEGER, INTENT(IN) ::  errno

       IF ( nc_stat /= NF90_NOERR )  THEN
          WRITE (*,'(A,1X,I3/4X,A)')                                           &
                                   '+++ read_palm_netcdf_data  error handle:', &
                                   errno, TRIM( nf90_strerror( nc_stat ) )
          STOP
       ENDIF

    END SUBROUTINE handle_netcdf_error


 END PROGRAM read_palm_netcdf_data



