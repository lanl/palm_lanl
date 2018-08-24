 PROGRAM combine_plot_fields

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
! $Id: combine_plot_fields.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2669 2017-12-06 16:03:27Z raasch
! data of 3d-nest runs are completely processed now
! 
! 2523 2017-10-05 14:42:47Z kanani
! Increased LEN for CHARACTER variable var_name, equal to the value in PALM
! 
! 2512 2017-10-04 08:26:59Z raasch
! PALM output does not contain ghost layer data any more
! avs- and iso2d-related parts removed, handling of compressed data removed
! 
! 2365 2017-08-21 14:59:59Z kanani
! Vertical grid nesting implemented (SadiqHuq)
!
! 1809 2016-04-05 20:13:28Z raasch
!
! 1808 2016-04-05 19:44:00Z raasch
! cpu measurements are done with standard FORTRAN routine on every machine
!
! 1551 2015-03-03 14:18:16Z maronga
! Adjustments for data output of soil model quantities
! 
! 1468 2014-09-24 14:06:57Z maronga
! Adapted for use on up to 6-digit processor cores (not tested)
!
! 1394 2014-05-06 10:17:31Z keck
! KIND-attribute added to all INTEGER and REAL declaration statements
!
! 1046 2012-11-09 14:38:45Z maronga
! code put under GPL (PALM 3.9)
!
! 691 2011-03-04 08:45:30Z maronga
! Bugfix: check for precursor ocean runs, removed typo
!
! 493 2010-03-01 08:30:24Z raasch
! Exit in case of already complete NetCDF data (due to parallel output in PALM)
! cpu measurements included
!
! 210 2008-11-06 08:54:02Z raasch
! Size of pf3d adjusted to the required output size (1 gridpoint less, along
! all three dimensions), because output of a subset of the data
! (pf3d(nxa:nxe...) in the NF90_PUT_VAR statement caused segmentation fault
! with the INTEL compiler.
! Subdomain data are read into temporary arrays pf_tmp/pf3d_tmp in order to
! avoid INTEL compiler warnings about (automatic) creation of temporary arrays
! Bugfix: three misplaced #endif directives
!
! 114 2007-10-10 00:03:15Z raasch
! Bugfix: model_string needed a default value
!
! Aug 07    Loop for processing of output by coupled runs, id_string does not
!           contain modus any longer
!
! 18/01/06  Output of time-averaged data
!
! 25/05/05  Errors removed
!
! 26/04/05  Output in NetCDF format, iso2d and avs output only if parameter
!           file exists
!
! 31/10/01  All comments and messages translated into English
!
! 23/02/99  Keine Bearbeitung komprimierter 3D-Daten
! Ursprungsversion vom 28/07/97
!
!
! Description:
! ------------
! This routine combines data of the PALM-subdomains into one file. In PALM
! every processor element opens its own file and writes 2D- or 3D-binary-data
! into it (different files are opened for xy-, xz-, yz-cross-sections and
! 3D-data). For plotting or analyzing these PE-data have to be collected and
! to be put into single files, which is done by this routine.
! Output format is NetCDF. Additionally, a data are output in a binary format
! readable by ISO2D-software (cross-sections) and by AVS (3D-data).
!------------------------------------------------------------------------------!

#if defined( __netcdf )
    USE netcdf
#endif

    IMPLICIT NONE

!
!-- Local variables
    CHARACTER (LEN=2)    ::  modus
    CHARACTER (LEN=4)    ::  model_string
    CHARACTER (LEN=6)    ::  id_string
    CHARACTER (LEN=30)   ::  dimname, var_name
    CHARACTER (LEN=40)   ::  filename

    CHARACTER (LEN=2000), DIMENSION(0:1) ::  var_list

    INTEGER, PARAMETER ::  iwp = 4, spk = SELECTED_REAL_KIND( 6 ), wp = 8

    INTEGER(iwp) ::  av, danz, i, id, j, k, model, models, nc_stat,            &
                     nxa, nxag, nxe, nxeg, nya, nyag, nye, nyeg,               &
                     nza, nze, pos, time_step, xa, xe, xxa, xxe, ya, ya_do,    &
                     ya_tot, ye, ye_do, ye_tot, yya, yye, za, ze, zza, zze

    INTEGER(8)                        ::  count, count_rate

    INTEGER(iwp), DIMENSION(0:1)      ::  current_level, current_var,          &
                                          fanz, id_set, id_var_time, num_var
    INTEGER(iwp), DIMENSION(4)        ::  id_dims_loc
    INTEGER(iwp), DIMENSION(0:1,4)    ::  id_dims
    INTEGER(iwp), DIMENSION(0:1,1000) ::  id_var, levels

    LOGICAL  ::  found, nest3d, netcdf_output, netcdf_parallel, netcdf_0,      &
                 netcdf_1, vnest

    REAL(wp) ::  cpu_start_time, cpu_end_time, dx, simulated_time
    REAL(wp),  DIMENSION(:,:), ALLOCATABLE   ::  pf, pf_tmp
    REAL(spk), DIMENSION(:,:,:), ALLOCATABLE ::  pf3d, pf3d_tmp



    PRINT*, ''
    PRINT*, ''
    PRINT*, '*** combine_plot_fields ***'

!
!-- Find out if a coupled or nested run has been carried out
    INQUIRE( FILE='COUPLING_PORT_OPENED', EXIST=found )
    INQUIRE( FILE='VNESTING_PORT_OPENED', EXIST=vnest )
    INQUIRE( FILE='3DNESTING', EXIST=nest3d )
    IF ( found )  THEN
       models = 2
       PRINT*, '    coupled run'
    ELSEIF ( vnest )  THEN
       models = 2
       PRINT*, '    Vertically nested grid coupling'
    ELSEIF ( nest3d )  THEN
       OPEN( 90, FILE='3DNESTING', FORM='FORMATTED' )
       READ ( 90, '(I2)' )  models
       CLOSE ( 90 )
       PRINT*, '    3d-nest run'
       PRINT*, '    number of nest domains = ', models
    ELSE
       models = 1
       PRINT*, '    uncoupled run'
    ENDIF

!
!-- Find out if a precursor ocean run has been carried out
    INQUIRE( FILE='PRECURSOR_OCEAN', EXIST=found )
    IF ( found )  THEN
       model_string = '_O'
       PRINT*, '    precursor ocean run'
    ELSE
       model_string = ''
    ENDIF

!
!-- Do everything for each model
    DO model = 1, models
!
!--    Set the model string used to identify the filenames
       IF ( found  .OR.  vnest )  THEN
          PRINT*, ''
          IF ( model == 2 )  THEN
             IF ( vnest )  THEN
                model_string = '_N'
                PRINT*, '    now combining FINE data'
                PRINT*, '    ========================'
             ELSE
                model_string = '_O'
                PRINT*, '    now combining ocean data'
                PRINT*, '    ========================'
             ENDIF
          ELSE
             IF ( vnest )  THEN
                PRINT*, '    now combining COARSE data'
                PRINT*, '    ============================='
             ELSE
                PRINT*, '    now combining atmosphere data'
                PRINT*, '    ============================='
             ENDIF
          ENDIF
       ELSEIF ( nest3d )  THEN
          PRINT*, ''
          PRINT*, '--> processing nest id = ', model
          IF ( model == 1 )  THEN
             model_string = ''
          ELSE
             WRITE( model_string, '(A2,I2.2)' )  '_N', model
          ENDIF
       ENDIF
!
!--    2D-arrays for ISO2D
!--    Main loop for the three different cross-sections, starting with 
!--    xy-section
       modus = 'XY'
       PRINT*, ''
       DO  WHILE ( modus == 'XY'  .OR.  modus == 'XZ'  .OR.  modus == 'YZ' )
!
!--       Take current time
          CALL SYSTEM_CLOCK( count, count_rate )
          cpu_start_time = REAL( count ) / REAL( count_rate )

          netcdf_parallel = .FALSE.
!
!--       Check, if file from PE0 exists. If it does not exist, PALM did not
!--       create any output for this cross-section.
          danz = 0
          WRITE (id_string,'(I6.6)')  danz
          INQUIRE ( &
               FILE='PLOT2D_'//modus//TRIM( model_string )//'_'//id_string, &
               EXIST=found )
!
!--       Find out the number of files (equal to the number of PEs which
!--       have been used in PALM) and open them
          DO  WHILE ( found )

             OPEN ( danz+110, &
                  FILE='PLOT2D_'//modus//TRIM( model_string )//'_'//id_string, &
                  FORM='UNFORMATTED' )
             danz = danz + 1
             WRITE (id_string,'(I6.6)')  danz
             INQUIRE ( &
                  FILE='PLOT2D_'//modus//TRIM( model_string )//'_'//id_string, &
                  EXIST=found )

          ENDDO

!
!--       Inquire whether a NetCDF file exists
          INQUIRE( FILE='DATA_2D_'//modus//'_NETCDF'//TRIM( model_string ), &
               EXIST=netcdf_0 )

!
!--       Inquire whether a NetCDF file for time-averaged data exists
          INQUIRE( FILE='DATA_2D_'//modus//'_AV_NETCDF'//TRIM( model_string ),&
               EXIST=netcdf_1 )

          IF ( netcdf_0  .OR.  netcdf_1 )  THEN
             netcdf_output = .TRUE.
!
!--          Inquire whether the NetCDF file is already complete (parallel
!--          output)
             INQUIRE( FILE='NO_COMBINE_PLOT_FIELDS_'//modus, &
                      EXIST=netcdf_parallel )
             IF ( netcdf_parallel )  THEN
                netcdf_parallel = .TRUE.
             ELSE
                netcdf_parallel = .FALSE.
             ENDIF
          ELSE
             netcdf_output = .FALSE.
          ENDIF

!
!--       Info-output
          PRINT*, ''
#if defined( __netcdf )
          IF ( netcdf_output )  THEN
             IF ( netcdf_parallel )  THEN
             PRINT*, '    NetCDF ' // modus // '-data are in one file ', &
                          '(NetCDF4-format) - merging not neccessary'
             ELSE
                PRINT*, '    NetCDF output enabled'
             ENDIF 
          ENDIF
#else
          IF ( netcdf_output )  THEN
             PRINT*, '--- Sorry, no NetCDF support on this host'
             netcdf_output = .FALSE.
          ENDIF
#endif
          IF ( .NOT. netcdf_parallel )  THEN
             IF ( danz /= 0 )  THEN
                PRINT*, '    ',modus,'-section:  ', danz, ' file(s) found'
             ELSE
                PRINT*, '    no ', modus, '-section data available'
             ENDIF
          ENDIF


          IF ( netcdf_output  .AND. .NOT. netcdf_parallel  .AND.  danz /= 0 ) &
          THEN
#if defined( __netcdf )
             DO  av = 0, 1

                IF ( av == 0  .AND.  .NOT.  netcdf_0 )  CYCLE
                IF ( av == 1  .AND.  .NOT.  netcdf_1 )  CYCLE

!
!--             Open NetCDF dataset
                IF ( av == 0 )  THEN
                   filename = 'DATA_2D_'//modus//'_NETCDF' &
                        //TRIM( model_string )
                ELSE
                   filename = 'DATA_2D_'//modus//'_AV_NETCDF' &
                        //TRIM( model_string )
                ENDIF
                nc_stat = NF90_OPEN( filename, NF90_WRITE, id_set(av) )
                IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 1 )

!
!--             Get the list of variables (order of variables corresponds with
!--             the order of data on the binary file)
                var_list(av) = ' '    ! GET_ATT does not assign trailing blanks
                nc_stat = NF90_GET_ATT( id_set(av), NF90_GLOBAL, 'VAR_LIST', &
                     var_list(av) )
                IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 2 )

!
!--             Inquire id of the time coordinate variable
                nc_stat = NF90_INQ_VARID( id_set(av), 'time', id_var_time(av) )
                IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 3 )

!
!--             Count number of variables; there is one more semicolon in the
!--             string than variable names
                num_var(av) = -1
                DO  i = 1, LEN( var_list(av) )
                   IF ( var_list(av)(i:i) == ';' )  num_var(av) = num_var(av) +1
                ENDDO

!
!--             Extract the variable names from the list and inquire their
!--             NetCDF IDs
                pos = INDEX( var_list(av), ';' )
!
!--             Loop over all variables
                DO  i = 1, num_var(av)

!
!--                Extract variable name from list
                   var_list(av) = var_list(av)(pos+1:)
                   pos = INDEX( var_list(av), ';' )
                   var_name = var_list(av)(1:pos-1)

!
!--                Get variable ID from name
                   nc_stat = NF90_INQ_VARID( id_set(av), TRIM( var_name ), &
                        id_var(av,i) )
                   IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 4 )

!
!--                Get number of x/y/z levels for that variable
                   nc_stat = NF90_INQUIRE_VARIABLE( id_set(av), id_var(av,i), &
                        dimids = id_dims_loc )
                   IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 5 )
                   id_dims(av,:) = id_dims_loc

!
!--                Inquire dimension ID
                   DO  j = 1, 4
                      nc_stat = NF90_INQUIRE_DIMENSION( id_set(av), &
                           id_dims(av,j), dimname, levels(av,i) )
                      IF ( nc_stat /= NF90_NOERR ) CALL handle_netcdf_error( 6 )

                      IF ( modus == 'XY' .AND. INDEX(dimname, 'z') /= 0 )  EXIT
                      IF ( modus == 'XZ' .AND. INDEX(dimname, 'y') /= 0 )  EXIT
                      IF ( modus == 'YZ' .AND. INDEX(dimname, 'x') /= 0 )  EXIT
                   ENDDO

                ENDDO

             ENDDO   ! av = 0, 1

#endif
          ENDIF

!
!--       Read the arrays, as long as the end of the file is reached
          IF ( .NOT. netcdf_parallel )  THEN

             fanz          =         0
             current_level =         1
             current_var   = 999999999

             DO  WHILE ( danz /= 0 )

!
!--             Loop over all files (reading data of the subdomains)
                DO  id = 0, danz-1
!
!--                File from PE0 contains special information at the beginning,
!--                concerning the lower and upper indices of the total-domain
!--                used in PALM (nxa, nxe, nya, nye).
!--                Allocate necessary arrays, open the output file and write
!--                the coordinate informations needed by ISO2D.
                   IF ( id == 0  .AND.  fanz(0) == 0  .AND.  fanz(1) == 0 ) THEN

                      READ ( id+110 )  nxa, nxe, nya, nye
                      ALLOCATE ( pf(nxa:nxe,nya:nye) )
!
!--                   Set actual domain bounds to total domain
                      ya_do = nya
                      ye_do = nye

                   ENDIF
!
!--                Read output time
                   IF ( netcdf_output  .AND.  id == 0 )  THEN
                      READ ( id+110, END=998 )  simulated_time, time_step, av
                   ENDIF
!
!--                Read subdomain indices
                   READ ( id+110, END=998 )  xa, xe, ya, ye, ya_tot, ye_tot

!
!--                IF the PE made no output (in case that no part of the
!--                cross-section is situated on this PE), indices have the
!--                value -1
                   IF ( .NOT. ( xa == -1  .AND.  xe == -1  .AND. &
                                ya == -1  .AND.  ye == -1 ) )  THEN


!
!--                   Read the subdomain grid-point values
                      ALLOCATE( pf_tmp(xa:xe,ya:ye) )
                      READ ( id+110 )  pf_tmp
!
!--                   Calculate indices on atmospheric grid (e.g. for soil model
!--                   quantities)
                      IF ( ya /= ya_tot .OR. ye /= ye_tot )  THEN
                         ye_do = ye - ya
                         ya_do = ya
                         pf(xa:xe,ya_do:ye_do) = pf_tmp
                      ELSE
                         ye_do = nye
                         ya_do = nya
                         pf(xa:xe,ya:ye) = pf_tmp
                      ENDIF

                      DEALLOCATE( pf_tmp )
                   ENDIF
                   IF ( id == 0 )  fanz(av) = fanz(av) + 1

                ENDDO
!
!--             Write data in NetCDF format
                IF ( netcdf_output )  THEN
#if defined( __netcdf )
!
!--                Check if a new time step has begun; if yes write data to
!--                time axis
                   IF ( current_var(av) > num_var(av) )  THEN
                      current_var(av) = 1
                      nc_stat = NF90_PUT_VAR( id_set(av), id_var_time(av), &
                                              (/ simulated_time /),        &
                                              start = (/ time_step /),     &
                                              count = (/ 1 /) )
                      IF ( nc_stat /= NF90_NOERR ) CALL handle_netcdf_error( 7 )
                   ENDIF

!
!--                Now write the data; this is mode dependent
                   SELECT CASE ( modus )

                      CASE ( 'XY' )
                         nc_stat = NF90_PUT_VAR( id_set(av),                   &
                                           id_var(av,current_var(av)),         &
                                           pf(nxa:nxe,nya:nye),                &
                             start = (/ 1, 1, current_level(av), time_step /), &
                                      count = (/ nxe-nxa+1, nye-nya+1, 1, 1 /) )
                         IF ( nc_stat /= NF90_NOERR )  THEN
                            CALL handle_netcdf_error( 8 )
                         ENDIF
                  
                      CASE ( 'XZ' )
                         nc_stat = NF90_PUT_VAR( id_set(av),                   &
                                           id_var(av,current_var(av)),         &
                                           pf(nxa:nxe,ya_do:ye_do),            &
                             start = (/ 1, current_level(av), 1, time_step /), &
                                      count = (/ nxe-nxa+1, 1, ye_do-ya_do+1, 1 /) )
                         IF ( nc_stat /= NF90_NOERR )  THEN
                            CALL handle_netcdf_error( 9 )
                         ENDIF

                      CASE ( 'YZ' )
                         nc_stat = NF90_PUT_VAR( id_set(av),                   &
                                           id_var(av,current_var(av)),         &
                                           pf(nxa:nxe,ya_do:ye_do),            &
                             start = (/ current_level(av), 1, 1, time_step /), &
                                      count = (/ 1, nxe-nxa+1, ye_do-ya_do+1, 1 /) )
                         IF ( nc_stat /= NF90_NOERR )  THEN
                            CALL handle_netcdf_error( 10 )
                         ENDIF

                   END SELECT

!
!--                Data is written, check if max level is reached
                   current_level(av) = current_level(av) + 1
                   IF ( current_level(av) > levels(av,current_var(av)) )  THEN
                      current_level(av) = 1
                      current_var(av)   = current_var(av) + 1
                   ENDIF

#endif
                ENDIF

             ENDDO

          ENDIF

998       IF ( danz /= 0  .AND.  .NOT. netcdf_parallel )  THEN
!
!--          Print the number of the arrays processed
             WRITE (*,'(16X,I4,A)')  fanz(0)+fanz(1), ' array(s) processed'
             IF ( fanz(1) /= 0 )  THEN
                WRITE (*,'(16X,I4,A)')  fanz(1), ' array(s) are time-averaged'
             ENDIF

!
!--          Close all files and deallocate arrays
             DO  id = 0, danz-1
                CLOSE ( id+110 )
             ENDDO
             CLOSE ( 2 )
             DEALLOCATE ( pf )

!
!--          Close the NetCDF file
             IF ( netcdf_output )  THEN
#if defined( __netcdf )
                IF ( netcdf_0 )  THEN
                   nc_stat = NF90_CLOSE( id_set(0) )
                   IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 11 )
                ENDIF
                IF ( netcdf_1 )  THEN
                   nc_stat = NF90_CLOSE( id_set(1) )
                   IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 12 )
                ENDIF
#endif
             ENDIF
          ENDIF

!
!--       Output required cpu time
          IF ( danz /= 0  .AND.  .NOT. netcdf_parallel )  THEN
             CALL SYSTEM_CLOCK( count, count_rate )
             cpu_end_time = REAL( count ) / REAL( count_rate )
             WRITE (*,'(5X,A,F9.3,A)')  'Required cpu-time: ', &
                                        cpu_end_time-cpu_start_time, ' sec'
          ENDIF

!
!--       Choose the next cross-section
          SELECT CASE ( modus )
             CASE ( 'XY' )
                modus = 'XZ'
             CASE ( 'XZ' )
                modus = 'YZ'
             CASE ( 'YZ' )
                modus = 'no'
          END SELECT

       ENDDO


!
!--    Combine the 3D-arrays
       netcdf_parallel = .FALSE.

!
!--    Info-output
       PRINT*, ' '

!
!--    Take current time
       CALL SYSTEM_CLOCK( count, count_rate )
       cpu_start_time = REAL( count ) / REAL( count_rate )

!
!--    Inquire whether a NetCDF file exists
       INQUIRE( FILE='DATA_3D_NETCDF'//TRIM( model_string ), EXIST=netcdf_0 )

!
!--    Inquire whether a NetCDF file for time-averaged data exists
       INQUIRE( FILE='DATA_3D_AV_NETCDF'//TRIM( model_string ), EXIST=netcdf_1 )

       IF ( netcdf_0  .OR.  netcdf_1 )  THEN
          netcdf_output = .TRUE.
!
!--       Inquire whether the NetCDF file is already complete (parallel output)
          INQUIRE( FILE='NO_COMBINE_PLOT_FIELDS_3D', EXIST=netcdf_parallel )
          IF ( netcdf_parallel )  THEN
             netcdf_parallel = .TRUE.
          ELSE
             netcdf_parallel = .FALSE.
          ENDIF
       ELSE
          netcdf_output = .FALSE.
       ENDIF

!
!--    Check, if file from PE0 exists; not neccessary in case of parallel
!--    PALM output
       IF ( .NOT. netcdf_parallel )  THEN
          danz = 0
          WRITE (id_string,'(I6.6)')  danz
          INQUIRE ( &
               FILE='PLOT3D_DATA'//TRIM( model_string )//'_'//TRIM( id_string ),  &
               EXIST=found )
       ELSE
          found = .FALSE.
       ENDIF

!
!--    Find out the number of files and open them
       DO  WHILE ( found )

          OPEN ( danz+110, &
               FILE='PLOT3D_DATA'//TRIM( model_string )//'_'//TRIM(id_string), &
               FORM='UNFORMATTED')
          danz = danz + 1
          WRITE (id_string,'(I6.6)')  danz
          INQUIRE ( &
               FILE='PLOT3D_DATA'//TRIM( model_string )//'_'//TRIM(id_string), &
               EXIST=found )

       ENDDO

#if defined( __netcdf )
       IF ( netcdf_output )  THEN
          IF ( netcdf_parallel )  THEN
             PRINT*, '    NetCDF data are in one file (NetCDF4-format)', &
                          ' - merging not neccessary'
          ELSE
             PRINT*, '    NetCDF output enabled'
          ENDIF
       ENDIF
#else
       IF ( netcdf_output )  THEN
          PRINT*, '--- Sorry, no NetCDF support on this host'
          netcdf_output = .FALSE.
       ENDIF
#endif
       IF ( .NOT. netcdf_parallel )  THEN
          IF ( danz /= 0 )  THEN
             PRINT*, '    3D-data:     ', danz, ' file(s) found'
          ELSE
             PRINT*, '    no 3D-data file available'
          ENDIF
       ENDIF

       IF ( netcdf_output  .AND. .NOT. netcdf_parallel  .AND.  danz /= 0 )  THEN
#if defined( __netcdf )
          DO  av = 0, 1

             IF ( av == 0  .AND.  .NOT.  netcdf_0 )  CYCLE
             IF ( av == 1  .AND.  .NOT.  netcdf_1 )  CYCLE

!
!--          Open NetCDF dataset
             IF ( av == 0 )  THEN
                filename = 'DATA_3D_NETCDF'//TRIM( model_string )
             ELSE
                filename = 'DATA_3D_AV_NETCDF'//TRIM( model_string )
             ENDIF
             nc_stat = NF90_OPEN( filename, NF90_WRITE, id_set(av) )
             IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 13 )


!
!--          Get the list of variables (order of variables corresponds with the
!--          order of data on the binary file)
             var_list(av) = ' '    ! GET_ATT does not assign trailing blanks
             nc_stat = NF90_GET_ATT( id_set(av), NF90_GLOBAL, 'VAR_LIST', &
                  var_list(av) )
             IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 14 )

!
!--          Inquire id of the time coordinate variable
             nc_stat = NF90_INQ_VARID( id_set(av), 'time', id_var_time(av) )
             IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 15 )

!
!--          Count number of variables; there is one more semicolon in the
!--          string than variable names
             num_var(av) = -1
             DO  i = 1, LEN( var_list(av) )
                IF ( var_list(av)(i:i) == ';' )  num_var(av) = num_var(av) + 1
             ENDDO

!
!--          Extract the variable names from the list and inquire their NetCDF
!--          IDs
             pos = INDEX( var_list(av), ';' )
!
!--          Loop over all variables
             DO  i = 1, num_var(av)

!
!--             Extract variable name from list
                var_list(av) = var_list(av)(pos+1:)
                pos = INDEX( var_list(av), ';' )
                var_name = var_list(av)(1:pos-1)

!
!--             Get variable ID from name
                nc_stat = NF90_INQ_VARID( id_set(av), TRIM( var_name ), &
                     id_var(av,i) )
                IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 16 )

             ENDDO

          ENDDO    ! av=0,1

#endif
       ENDIF

!
!--    Read arrays, until the end of the file is reached
       IF ( .NOT. netcdf_parallel )  THEN

          current_var = 999999999
          fanz = 0
          DO  WHILE ( danz /= 0 )

!
!--          Loop over all files
             DO  id = 0, danz-1
!
!--             File from PE0 contains at the beginning the index bounds 
!--             of PALM's total domain.
!--             Allocate the array for storing the total domain data
                IF ( id == 0  .AND.  fanz(0) == 0  .AND.  fanz(1) == 0 )  THEN
!                   READ ( id+110 )  nxag, nxeg, nyag, nyeg
                   READ ( id+110 )  nxa, nxe, nya, nye, nza, nze
                   ALLOCATE ( pf3d(nxa:nxe,nya:nye,nza:nze) )
                ENDIF

!
!--             Read output time
                IF ( netcdf_output  .AND.  id == 0 )  THEN
                   IF ( netcdf_1 )  THEN
                      READ ( id+110, END=999 )  simulated_time, time_step, av
                   ELSE
!
!--                   For compatibility with earlier PALM versions
                      READ ( id+110, END=999 )  simulated_time, time_step
                      av = 0
                   ENDIF
                ENDIF

!
!--             Read subdomain indices and grid point values
                READ ( id+110, END=999 )  xa, xe, ya, ye, za, ze
                ALLOCATE( pf3d_tmp(xa:xe,ya:ye,za:ze) )
                READ ( id+110 )  pf3d_tmp

                xxa = MAX( nxa, xa )
                xxe = MIN( nxe, xe )
                yya = MAX( nya, ya )
                yye = MIN( nye, ye )
                DO  k = za, ze
                   DO  j = yya, yye
                      DO  i = xxa, xxe
                         pf3d(i,j,k) = pf3d_tmp(i,j,k)
                      ENDDO
                   ENDDO
                ENDDO

                DEALLOCATE( pf3d_tmp )
                IF ( id == 0 )  fanz(av) = fanz(av) + 1

             ENDDO

!
!--          Write data of the total domain in NetCDF format
             IF ( netcdf_output )  THEN
#if defined( __netcdf )
!
!--             Check if a new time step has begun; if yes write data to time
!--             axis
                IF ( current_var(av) > num_var(av) )  THEN
                   current_var(av) = 1
                   nc_stat = NF90_PUT_VAR( id_set(av), id_var_time(av), &
                                      (/ simulated_time /),&
                                      start = (/ time_step /), count = (/ 1 /) )
                   IF ( nc_stat /= NF90_NOERR ) CALL handle_netcdf_error( 17 )
                ENDIF

!
!--             Now write the data
                nc_stat = NF90_PUT_VAR( id_set(av), id_var(av,current_var(av)),&
                                        pf3d(nxa:nxe,nya:nye,za:ze), start = (/ 1, 1, 1, time_step /),&
                              count = (/ nxe-nxa+1, nye-nya+1, ze-za+1, 1 /) )
                IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 18 )

                current_var(av) = current_var(av) + 1

#endif
             ENDIF

          ENDDO

       ENDIF


999    IF ( danz /= 0  .AND.  .NOT. netcdf_parallel )  THEN
!
!--       Print the number of arrays processed
          WRITE (*,'(16X,I4,A)')  fanz(0)+fanz(1), ' array(s) processed'
          IF ( fanz(1) /= 0 )  THEN
             WRITE (*,'(16X,I4,A)')  fanz(1), ' array(s) are time-averaged'
          ENDIF
!
!--       Close all files and deallocate array
          DO  id = 0, danz-1
             CLOSE ( id+110 )
          ENDDO
          CLOSE ( 2 )
          DEALLOCATE ( pf3d )
!
!--       Close the NetCDF file
          IF ( netcdf_output )  THEN
#if defined( __netcdf )
             IF ( netcdf_0 )  THEN
                nc_stat = NF90_CLOSE( id_set(0) )
                IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 19 )
             ENDIF
             IF ( netcdf_1 )  THEN
                nc_stat = NF90_CLOSE( id_set(1) )
                IF ( nc_stat /= NF90_NOERR )  CALL handle_netcdf_error( 20 )
             ENDIF
#endif
          ENDIF

!
!--       Output required cpu time
          CALL SYSTEM_CLOCK( count, count_rate )
          cpu_end_time = REAL( count ) / REAL( count_rate )
          WRITE (*,'(5X,A,F9.3,A)')  'Required cpu-time: ', &
                                     cpu_end_time-cpu_start_time, ' sec'

       ENDIF

    ENDDO  ! models


 CONTAINS


    SUBROUTINE handle_netcdf_error( errno )
!
!--    Prints out a text message corresponding to the current NetCDF status

       IMPLICIT NONE

       INTEGER, INTENT(IN) ::  errno

#if defined( __netcdf )
       IF ( nc_stat /= NF90_NOERR )  THEN
          PRINT*, '+++ combine_plot_fields  netcdf: ', av, errno, &
                  TRIM( nf90_strerror( nc_stat ) )
       ENDIF
#endif

    END SUBROUTINE handle_netcdf_error


 END PROGRAM combine_plot_fields



