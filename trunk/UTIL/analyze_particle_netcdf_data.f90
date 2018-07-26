 PROGRAM analyze_particle_netcdf_data

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
! $Id: analyze_particle_netcdf_data.f90 2718 2018-01-02 08:49:38Z maronga $
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
! Initial revision 12/07/05
!
!
! Description:
! ------------
! This is an EXAMPLE routine how to read NetCDF particle data output by PALM.
! As an example, the mean heigth and radius of all particles are calculated
! and output for each output time available on the file
!
! 
! Any additonal analyzation requested has to be added by the user by following
! the steps given in this example!
!
!
!
! This routine must be compiled with:
! decalpha:
!    f95 -fast -r8 -I/usr/local/netcdf-3.5.1/include
!    -L/usr/local/netcdf-3.5.1/lib -lnetcdf
! IBM-Regatta:
!    xlf95 -qsuffix=cpp=f90 -qrealsize=8 -q64 -qmaxmem=-1 -Q
!    -I /aws/dataformats/netcdf-3.5.0/netcdf-64-32-3.5.0/include_F90_64
!    -L/aws/dataformats/netcdf-3.5.0/netcdf-64-32-3.5.0/lib -lnetcdf -O3
! IBM-Regatta KISTI:
!    xlf95 -qsuffix=cpp=f90 -qrealsize=8 -q64 -qmaxmem=-1 -Q
!    -I /applic/netcdf64/src/f90
!    -L/applic/lib/NETCDF64 -lnetcdf -O3
! IMUK:
!    ifort analyze_particle_netcdf_data
!    -I /muksoft/packages/netcdf/linux/include -axW -r8 -nbs -Vaxlib
!    -L /muksoft/packages/netcdf/linux/lib -lnetcdf
! NEC-SX6:
!    sxf90 analyze_particle_netcdf_data.f90 
!    -I /pool/SX-6/netcdf/netcdf-3.6.0-p1/include  -C hopt -Wf '-A idbl4'
!    -D__netcdf -L/pool/SX-6/netcdf/netcdf-3.6.0-p1/lib -lnetcdf
!------------------------------------------------------------------------------!

    USE netcdf

    IMPLICIT NONE

!
!-- Local variables
    CHARACTER (LEN=7)    ::  id_char
    CHARACTER (LEN=2000) ::  string


    INTEGER ::  f, fn, i, id_dim_time, id_var_rnop, id_var_r, id_var_time, &
                id_var_x, id_var_y, id_var_z, n, nc_stat, ntl, start

    INTEGER, DIMENSION(1000) ::  id_set
    INTEGER, DIMENSION(1)    ::  id_dim_time_old

    INTEGER, DIMENSION(:), ALLOCATABLE   ::  total_nop
    INTEGER, DIMENSION(:,:), ALLOCATABLE ::  nop

    LOGICAL ::  found

    REAL ::  mean_r, mean_z
    REAL, DIMENSION(:), ALLOCATABLE ::  prt_r, prt_x, prt_y, prt_z, tl

!
!-- Check, if file from PE0 exists. If it does not exist, PALM did not
!-- create any output for this cross-section.
    fn = 0
    WRITE (id_char,'(''_'',I4.4)')  fn
    INQUIRE ( FILE=id_char, EXIST=found )

!
!-- Find out the number of files (equal to the number of PEs which
!-- have been used in PALM) and open them
    IF ( .NOT. found )  THEN
       PRINT*, '+++ no file _0000 found in current working directory'
       STOP
    ENDIF

    DO  WHILE ( found )

!
!--    Open NetCDF dataset
       nc_stat = NF90_OPEN( id_char, NF90_NOWRITE, id_set(fn) )
       IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( nc_stat, 1 )
       fn = fn + 1
       WRITE (id_char,'(''_'',I4.4)')  fn
       INQUIRE ( FILE=id_char, EXIST=found )

    ENDDO
    fn = fn - 1

    PRINT*, '*** ', fn+1, ' files found'

!
!-- Get the run description header and print it out
    string = ' '
    nc_stat = NF90_GET_ATT( id_set(0), NF90_GLOBAL, 'title', string )
    IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( nc_stat, 2 )
    PRINT*, '*** run: ', TRIM( string )

!
!-- Get the available time levels
    nc_stat = NF90_INQ_VARID( id_set(0), 'time', id_var_time )
    IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( nc_stat, 3 )

    nc_stat = NF90_INQUIRE_VARIABLE( id_set(0), id_var_time, &
                                     dimids = id_dim_time_old )
    IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( nc_stat, 4 )
    id_dim_time = id_dim_time_old(1)

    nc_stat = NF90_INQUIRE_DIMENSION( id_set(0), id_dim_time, len = ntl )
    IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( nc_stat, 5 )
    ALLOCATE( tl(1:ntl) )
    print*, 'ntl=',ntl

    nc_stat = NF90_GET_VAR( id_set(0), id_var_time, tl )
    IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( nc_stat, 6 )

    DO  n = 1, ntl
       print*, '*** time_level(', n, ') =', tl(n)
    ENDDO

!
!-- Get the number of particles used
    nc_stat = NF90_INQ_VARID( id_set(0), 'real_num_of_prt', id_var_rnop )
    IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( nc_stat, 7 )

    ALLOCATE( nop(1:ntl,0:fn), total_nop(1:ntl) )

    DO  f = 0, fn

       nc_stat = NF90_GET_VAR( id_set(f), id_var_rnop, nop(1:ntl,f) )
       IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( nc_stat, 8 )

    ENDDO

    total_nop = 0
    DO  n = 1, ntl

       DO  f = 0, fn
          total_nop(n) = total_nop(n) + nop(n,f)
       ENDDO

       PRINT*, '*** time = ', tl(n), ' total # of particles: ', total_nop(n)

    ENDDO

!
!-- Get the particle x and y coordinates
    nc_stat = NF90_INQ_VARID( id_set(0), 'pt_x', id_var_x )
    IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( nc_stat, 9 )

    nc_stat = NF90_INQ_VARID( id_set(0), 'pt_y', id_var_y )
    IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( nc_stat, 10 )

    nc_stat = NF90_INQ_VARID( id_set(0), 'pt_z', id_var_z )
    IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( nc_stat, 11 )

    nc_stat = NF90_INQ_VARID( id_set(0), 'pt_radius', id_var_r )
    IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( nc_stat, 12 )

    PRINT*, ' '
!
!-- Loop over all timelevels
    DO  n = 1, ntl

       ALLOCATE( prt_x(total_nop(n)), prt_y(total_nop(n)), &
                 prt_z(total_nop(n)), prt_r(total_nop(n)) )

       start = 1

!
!--    Read the data from the files (one file per processor)
       DO  f = 0, fn

          nc_stat = NF90_GET_VAR( id_set(f), id_var_x,           &
                                  prt_x(start:start+nop(n,f)-1), &
                                  start = (/ 1, n /),            &
                                  count = (/ nop(n,f), 1 /) )
          IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( nc_stat, 13 )

          nc_stat = NF90_GET_VAR( id_set(f), id_var_y,           &
                                  prt_y(start:start+nop(n,f)-1), &
                                  start = (/ 1, n /),            &
                                  count = (/ nop(n,f), 1 /) )
          IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( nc_stat, 14 )

          nc_stat = NF90_GET_VAR( id_set(f), id_var_z,           &
                                  prt_z(start:start+nop(n,f)-1), &
                                  start = (/ 1, n /),            &
                                  count = (/ nop(n,f), 1 /) )
          IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( nc_stat, 15 )

          nc_stat = NF90_GET_VAR( id_set(f), id_var_r,           &
                                  prt_r(start:start+nop(n,f)-1), &
                                  start = (/ 1, n /),            &
                                  count = (/ nop(n,f), 1 /) )
          IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( nc_stat, 16 )

          start = start + nop(n,f)

       ENDDO

       mean_z = 0.0
       mean_r = 0.0
       DO  i = 1, total_nop(n)
          mean_z = mean_z + prt_z(i)
          mean_r = mean_r + prt_r(i)
       ENDDO
       mean_z = mean_z / total_nop(n)
       mean_r = mean_r / total_nop(n)

       PRINT*, '*** time = ', tl(n), ' mean height = ', mean_z, &
                                     ' mean radius = ', mean_r

!
!--    prt_x, prt_y, prt_z, and prt_r contain the particle coordinates and
!--    radii, respectively. Please output these data or carry out the
!--    requested analyzing in this program before you deallocate the arrays.

       DEALLOCATE( prt_x, prt_y, prt_z, prt_r )

    ENDDO


 END PROGRAM analyze_particle_netcdf_data



 SUBROUTINE handle_netcdf_error( nc_stat, position )
!
!-- Prints out a text message corresponding to the current NetCDF status

    USE netcdf

    IMPLICIT NONE

    INTEGER, INTENT(IN) ::  nc_stat, position

    IF ( nc_stat /= NF90_NOERR )  THEN
       PRINT*, '+++ analyze_particle_netcdf_data'
       PRINT*, '    netcdf error: ', TRIM( nf90_strerror( nc_stat ) )
       PRINT*, '    position = ', position
       STOP
    ENDIF

 END SUBROUTINE handle_netcdf_error
