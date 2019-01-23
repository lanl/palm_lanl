!> @file check_open.f90
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
! Add file_id for runfile which specifies coupling parameters
! 
! Former revisions:
! -----------------
! $Id: check_open.f90 3045 2018-05-28 07:55:41Z Giersch $
! Output case 108 disabled
! 
! 2964 2018-04-12 16:04:03Z Giersch
! Error message moved to radiation_model_mod
! 
! 2957 2018-04-11 08:48:06Z Giersch
! Error message has been corrected in case of writing sky view factors
! 
! 2906 2018-03-19 08:56:40Z Giersch
! CASE 88 and 89 has been added for the input/output of sky view factors
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2669 2017-12-06 16:03:27Z raasch
! file name extension for masked data files is changed to "_M##" and is now
! appended at the end of the filename,
! file ids not used any more have been removed
! 
! 2516 2017-10-04 11:03:04Z suehring
! Remove tabs
! 
! 2514 2017-10-04 09:52:37Z suehring
! upper bounds of cross section and 3d output changed from nx+1,ny+1 to nx,ny
! no output of ghost layer data
! iso2d-related parts removed
! 
! 2300 2017-06-29 13:31:14Z raasch
! -host
! 
! 2298 2017-06-29 09:28:18Z raasch
! -return_addres, return_username, avs_coor_file_..., avs_data_file_...,
! cross_ts_numbers, cross_ts_number_count
!
! 2101 2017-01-05 16:42:31Z suehring
!
! 2063 2016-11-10 17:14:35Z raasch
! bugfix: opening of PROGRESS file moved out of the NetCDF block
!
! 2040 2016-10-26 16:58:09Z gronemeier
! Removed open of file 'PLOTTS_DATA' ( CASE(50:59) ) as it is no longer needed
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1988 2016-08-10 14:49:06Z gronemeier
! informative message added if files cannot be opened in newly created directory
! 
! 1986 2016-08-10 14:07:17Z gronemeier
! Bugfix: check if output can be opened in newly created directory. If not
! wait one second and try again.
! 
! 1974 2016-07-26 08:43:25Z gronemeier
! Bugfix: MPI barriers after deleting non-extendable files must only be called
! in case of parallel I/O
! 
! 1957 2016-07-07 10:43:48Z suehring
! flight module added
!
! 1804 2016-04-05 16:30:18Z maronga
! Removed code for parameter file check (__check)
!
! 1783 2016-03-06 18:36:17Z raasch
! name change of netcdf routines and module + related changes
!
! 1779 2016-03-03 08:01:28Z raasch
! coupling_char is trimmed at every place it occurs, because it can have
! different length now
!
! 1745 2016-02-05 13:06:51Z gronemeier
! Bugfix: added MPI barrier after deleting existing non-extendable file by PE0
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1551 2015-03-03 14:18:16Z maronga
! Removed redundant output for combine_plot_fields
! 
! 1468 2014-09-24 14:06:57Z maronga
! Adapted for use on up to 6-digit processor cores
! Added file unit 117 (PROGRESS)
! 
! 1359 2014-04-11 17:15:14Z hoffmann
! Format of particle exchange statistics extended to reasonable numbers of      
! particles.
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute,
! declaration for unused variables xkoor, ykoor, zkoor removed
! 
! 1327 2014-03-21 11:00:16Z raasch
! parts concerning iso2d and avs output removed
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! old module precision_kind is removed,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1106 2013-03-04 05:31:38Z raasch
! array_kind renamed precision_kind
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1031 2012-10-19 14:35:30Z raasch
! netCDF4 without parallel file support implemented,
! opening of netCDF files are done by new routines create_netcdf_file and
! open_write_netcdf_file
!
! 964 2012-07-26 09:14:24Z raasch
! old profil-units (40:49) removed,
! append feature removed from unit 14
!
! 849 2012-03-15 10:35:09Z raasch
! comment changed
!
! 809 2012-01-30 13:32:58Z maronga
! Bugfix: replaced .AND. and .NOT. with && and ! in the preprocessor directives
!
! 807 2012-01-25 11:53:51Z maronga
! New cpp directive "__check" implemented which is used by check_namelist_files
!
! Revision 1.1  1997/08/11 06:10:55  raasch
! Initial revision
!
!
! Description:
! ------------
!> Check if file unit is open. If not, open file and, if necessary, write a
!> header or start other initializing actions, respectively.
!------------------------------------------------------------------------------!
SUBROUTINE check_open( file_id )
 

    USE arrays_3d,                                                             &
        ONLY:  zu

    USE control_parameters,                                                    &
        ONLY:  coupling_char,                      &
               max_masks, message_string, mid, nz_do3d, openfile,              &
               run_description_header, runnr

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxr, ny, nyn, nyng, nys, nysg, nz, nzb, nzt 

    USE kinds

#if defined( __netcdf )
    USE NETCDF
#endif

    USE netcdf_interface,                                                      &
        ONLY:  id_set_pr, id_set_prt, id_set_pts,      &
               id_set_ts,          &
               id_set_3d, nc_stat, netcdf_create_file, netcdf_data_format,     &
               netcdf_define_header, netcdf_handle_error, netcdf_open_write_file
    USE pegrid

    USE posix_calls_from_fortran,                                              &
        ONLY:  fortran_sleep


    IMPLICIT NONE

    CHARACTER (LEN=4)   ::  mask_char               !<
    CHARACTER (LEN=2)   ::  suffix                  !<
    CHARACTER (LEN=30)  ::  filename                !<
    CHARACTER (LEN=80)  ::  rtext                   !<
    CHARACTER (LEN=100) ::  line                    !<

    INTEGER(iwp) ::  av          !<
    INTEGER(iwp) ::  file_id     !<
    INTEGER(iwp) ::  i           !<
    INTEGER(iwp) ::  ioerr       !< IOSTAT flag for IO-commands ( 0 = no error )
    INTEGER(iwp) ::  j           !<
    INTEGER(iwp) ::  k           !<
    
    LOGICAL ::  netcdf_extend, ISOPEN    !<

!
!-- Immediate return if file already open
    IF ( openfile(file_id)%opened )  RETURN

!
!-- Only certain files are allowed to be re-opened
!-- NOTE: some of the other files perhaps also could be re-opened, but it
!--       has not been checked so far, if it works!
    IF ( openfile(file_id)%opened_before )  THEN
       SELECT CASE ( file_id )
          CASE ( 13, 14, 21, 22, 23, 80, 85, 117 )
             IF ( file_id == 14 .AND. openfile(file_id)%opened_before )  THEN
                message_string = 're-open of unit ' //                         &
                                 '14 is not verified. Please check results!'
                CALL message( 'check_open', 'PA0165', 0, 1, 0, 6, 0 )       
             ENDIF

          CASE DEFAULT
             WRITE( message_string, * ) 're-opening of file-id ', file_id,     &
                                        ' is not allowed'
             CALL message( 'check_open', 'PA0166', 0, 1, 0, 6, 0 )    
               
             RETURN

       END SELECT
    ENDIF

!
!-- Check if file may be opened on the relevant PE
    SELECT CASE ( file_id )

       CASE ( 15, 16, 17, 18, 19, 50:59, 104:105, 107, 109, 117 )
      
          IF ( myid /= 0 )  THEN
             WRITE( message_string, * ) 'opening file-id ',file_id,            &
                                        ' not allowed for PE ',myid
             CALL message( 'check_open', 'PA0167', 2, 2, -1, 6, 1 )
          ENDIF

       CASE ( 101:103, 106, 111:113, 116, 201:200+2*max_masks )

          IF ( netcdf_data_format < 5 )  THEN
          
             IF ( myid /= 0 )  THEN
                WRITE( message_string, * ) 'opening file-id ',file_id,         &
                                           ' not allowed for PE ',myid
                CALL message( 'check_open', 'PA0167', 2, 2, -1, 6, 1 )
             ENDIF
     
          ENDIF

       CASE ( 90:99 )

!
!--       File-ids that are used temporarily in other routines
          WRITE( message_string, * ) 'opening file-id ',file_id,               &
                                    ' is not allowed since it is used otherwise'
          CALL message( 'check_open', 'PA0168', 0, 1, 0, 6, 0 ) 
          
    END SELECT

!
!-- Open relevant files
    SELECT CASE ( file_id )

       CASE ( 11 )

          OPEN ( 11, FILE='PARIN'//TRIM( coupling_char ), FORM='FORMATTED',    &
                     STATUS='OLD' )

       CASE ( 13 )

          IF ( myid_char == '' )  THEN
             OPEN ( 13, FILE='BININ'//TRIM( coupling_char )//myid_char,        &
                        FORM='UNFORMATTED', STATUS='OLD' )
          ELSE
!
!--          First opening of unit 13 openes file _000000 on all PEs because
!--          only this file contains the global variables
             IF ( .NOT. openfile(file_id)%opened_before )  THEN
                OPEN ( 13, FILE='BININ'//TRIM( coupling_char )//'/_000000',    &
                           FORM='UNFORMATTED', STATUS='OLD' )
             ELSE
                OPEN ( 13, FILE='BININ'//TRIM( coupling_char )//'/'//          &
                           myid_char, FORM='UNFORMATTED', STATUS='OLD' )
             ENDIF
          ENDIF

       CASE ( 14 )

          IF ( myid_char == '' )  THEN
             OPEN ( 14, FILE='BINOUT'//TRIM( coupling_char )//myid_char,       &
                        FORM='UNFORMATTED', POSITION='APPEND' )
          ELSE
             IF ( myid == 0  .AND. .NOT. openfile(file_id)%opened_before )  THEN
                CALL local_system( 'mkdir  BINOUT' // TRIM( coupling_char ) )
             ENDIF
#if defined( __parallel )
!
!--          Set a barrier in order to allow that all other processors in the 
!--          directory created by PE0 can open their file
             CALL MPI_BARRIER( comm2d, ierr )
#endif
             ioerr = 1
             DO WHILE ( ioerr /= 0 )
                OPEN ( 14, FILE='BINOUT'//TRIM(coupling_char)//'/'//myid_char, &
                           FORM='UNFORMATTED', IOSTAT=ioerr )
                IF ( ioerr /= 0 )  THEN
                   WRITE( 9, * )  '*** could not open "BINOUT'//         &
                                  TRIM(coupling_char)//'/'//myid_char//  &
                                  '"! Trying again in 1 sec.'
                   CALL fortran_sleep( 1 )
                ENDIF
             ENDDO

          ENDIF

       CASE ( 15 )

          OPEN ( 15, FILE='RUN_CONTROL'//TRIM( coupling_char ),                &
                     FORM='FORMATTED' )

       CASE ( 16 )

          OPEN ( 16, FILE='LIST_PROFIL'//TRIM( coupling_char ),                &
                     FORM='FORMATTED' )

       CASE ( 17 )

          OPEN ( 17, FILE='LIST_PROFIL_1D'//TRIM( coupling_char ),             &
                     FORM='FORMATTED' )

       CASE ( 18 )

          INQUIRE( UNIT=18, OPENED=ISOPEN )
          if( .not. ISOPEN) THEN
             OPEN ( 18, FILE='CPU_MEASURES'//TRIM( coupling_char ),               &
                     FORM='FORMATTED' )
          endif

       CASE ( 19 )

          OPEN ( 19, FILE='HEADER'//TRIM( coupling_char ), FORM='FORMATTED' )

       CASE ( 20 )

          IF ( myid == 0  .AND. .NOT. openfile(file_id)%opened_before )  THEN
             CALL local_system( 'mkdir  DATA_LOG' // TRIM( coupling_char ) )
          ENDIF
          IF ( myid_char == '' )  THEN
             OPEN ( 20, FILE='DATA_LOG'//TRIM( coupling_char )//'/_000000',    &
                        FORM='UNFORMATTED', POSITION='APPEND' )
          ELSE
#if defined( __parallel )
!
!--          Set a barrier in order to allow that all other processors in the 
!--          directory created by PE0 can open their file
             CALL MPI_BARRIER( comm2d, ierr )
#endif
             ioerr = 1
             DO WHILE ( ioerr /= 0 )
                OPEN ( 20, FILE='DATA_LOG'//TRIM( coupling_char )//'/'//       &
                           myid_char, FORM='UNFORMATTED', POSITION='APPEND',   &
                           IOSTAT=ioerr )
                IF ( ioerr /= 0 )  THEN
                   WRITE( 9, * )  '*** could not open "DATA_LOG'//         &
                                  TRIM( coupling_char )//'/'//myid_char//  &
                                  '"! Trying again in 1 sec.'
                   CALL fortran_sleep( 1 )
                ENDIF
             ENDDO

          ENDIF

       CASE ( 30 )

          OPEN ( 30, FILE='PLOT3D_DATA'//TRIM( coupling_char )//myid_char,     &
                     FORM='UNFORMATTED' )
!
!--       Specifications for combine_plot_fields
          IF ( myid == 0 )  THEN
#if defined( __parallel )
             WRITE ( 30 )  0, nx, 0, ny, 0, nz_do3d
#endif
          ENDIF

       CASE ( 50 )

          OPEN ( 50, FILE='runfile')


!--    File where sky-view factors and further required data is stored will be 
!--    read
       CASE ( 88 )

          IF ( myid_char == '' )  THEN
             OPEN ( 88, FILE='SVFIN'//TRIM( coupling_char )//myid_char,        &
                        FORM='UNFORMATTED', STATUS='OLD', IOSTAT=ioerr )
          ELSE

             OPEN ( 88, FILE='SVFIN'//TRIM( coupling_char )//'/'//myid_char,   &
                        FORM='UNFORMATTED', STATUS='OLD', IOSTAT=ioerr )
          ENDIF

!
!--    File where sky-view factors and further required data is stored will be
!--    created
       CASE ( 89 )

          IF ( myid_char == '' )  THEN
             OPEN ( 89, FILE='SVFOUT'//TRIM( coupling_char )//myid_char,       &
                        FORM='UNFORMATTED', STATUS='NEW' )
          ELSE
             IF ( myid == 0  .AND. .NOT. openfile(file_id)%opened_before )  THEN
                CALL local_system( 'mkdir  SVFOUT' // TRIM( coupling_char ) )
             ENDIF
#if defined( __parallel )
!
!--          Set a barrier in order to allow that all other processors in the 
!--          directory created by PE0 can open their file
             CALL MPI_BARRIER( comm2d, ierr )
#endif
             ioerr = 1
             DO WHILE ( ioerr /= 0 )
                OPEN ( 89, FILE='SVFOUT'//TRIM(coupling_char)//'/'//myid_char, &
                           FORM='UNFORMATTED', STATUS='NEW', IOSTAT=ioerr )
                IF ( ioerr /= 0 )  THEN
                   WRITE( 9, * )  '*** could not open "SVFOUT'//               &
                                  TRIM(coupling_char)//'/'//myid_char//        &
                                  '"! Trying again in 1 sec.'
                   CALL fortran_sleep( 1 )
                ENDIF
             ENDDO

          ENDIF

!
!--    Progress file that is used by the PALM watchdog
       CASE ( 117 )

          OPEN ( 117, FILE='PROGRESS'//TRIM( coupling_char ),                  &
                      STATUS='REPLACE', FORM='FORMATTED' )

#if defined( __netcdf )
       CASE ( 104 )
!
!--       Set filename
          filename = 'DATA_1D_PR_NETCDF' // TRIM( coupling_char )

!
!--       Inquire, if there is a netCDF file from a previuos run. This should
!--       be opened for extension, if its variables match the actual run.
          INQUIRE( FILE=filename, EXIST=netcdf_extend )

          IF ( netcdf_extend )  THEN
!
!--          Open an existing netCDF file for output
             CALL netcdf_open_write_file( filename, id_set_pr, .FALSE., 29 )
!
!--          Read header information and set all ids. If there is a mismatch
!--          between the previuos and the actual run, netcdf_extend is returned
!--          as .FALSE.
             CALL netcdf_define_header( 'pr', netcdf_extend, 0 )

!
!--          Remove the local file, if it can not be extended
             IF ( .NOT. netcdf_extend )  THEN
                nc_stat = NF90_CLOSE( id_set_pr )
                CALL netcdf_handle_error( 'check_open', 30 )
                CALL local_system( 'rm ' // TRIM( filename ) )
             ENDIF

          ENDIF          

          IF ( .NOT. netcdf_extend )  THEN
!
!--          Create a new netCDF output file with requested netCDF format
             CALL netcdf_create_file( filename, id_set_pr, .FALSE., 31 )
!
!--          Define the header
             CALL netcdf_define_header( 'pr', netcdf_extend, 0 )

          ENDIF

       CASE ( 105 )
!
!--       Set filename
          filename = 'DATA_1D_TS_NETCDF' // TRIM( coupling_char )

!
!--       Inquire, if there is a netCDF file from a previuos run. This should
!--       be opened for extension, if its variables match the actual run.
          INQUIRE( FILE=filename, EXIST=netcdf_extend )

          IF ( netcdf_extend )  THEN
!
!--          Open an existing netCDF file for output
             CALL netcdf_open_write_file( filename, id_set_ts, .FALSE., 32 )
!
!--          Read header information and set all ids. If there is a mismatch
!--          between the previuos and the actual run, netcdf_extend is returned
!--          as .FALSE.
             CALL netcdf_define_header( 'ts', netcdf_extend, 0 )

!
!--          Remove the local file, if it can not be extended
             IF ( .NOT. netcdf_extend )  THEN
                nc_stat = NF90_CLOSE( id_set_ts )
                CALL netcdf_handle_error( 'check_open', 33 )
                CALL local_system( 'rm ' // TRIM( filename ) )
             ENDIF

          ENDIF          

          IF ( .NOT. netcdf_extend )  THEN
!
!--          Create a new netCDF output file with requested netCDF format
             CALL netcdf_create_file( filename, id_set_ts, .FALSE., 34 )
!
!--          Define the header
             CALL netcdf_define_header( 'ts', netcdf_extend, 0 )

          ENDIF


       CASE ( 106, 116 )
!
!--       Set filename depending on unit number
          IF ( file_id == 106 )  THEN
             filename = 'DATA_3D_NETCDF' // TRIM( coupling_char )
             av = 0
          ELSE
             filename = 'DATA_3D_AV_NETCDF' // TRIM( coupling_char )
             av = 1
          ENDIF
!
!--       Inquire, if there is a netCDF file from a previous run. This should
!--       be opened for extension, if its dimensions and variables match the
!--       actual run.
          INQUIRE( FILE=filename, EXIST=netcdf_extend )
          IF ( netcdf_extend )  THEN
!
!--          Open an existing netCDF file for output
             CALL netcdf_open_write_file( filename, id_set_3d(av), .TRUE., 35 )
!
!--          Read header information and set all ids. If there is a mismatch
!--          between the previuos and the actual run, netcdf_extend is returned
!--          as .FALSE.
             CALL netcdf_define_header( '3d', netcdf_extend, av )

!
!--          Remove the local file, if it can not be extended
             IF ( .NOT. netcdf_extend )  THEN
                nc_stat = NF90_CLOSE( id_set_3d(av) )
                CALL netcdf_handle_error( 'check_open', 36 )
                IF ( myid == 0 )  CALL local_system( 'rm ' // TRIM( filename ) )
#if defined( __parallel )
!
!--             Set a barrier in order to assure that PE0 deleted the old file
!--             before any other processor tries to open a new file
!--             Barrier is only needed in case of parallel I/O
                IF ( netcdf_data_format > 4 )  CALL MPI_BARRIER( comm2d, ierr )
#endif
             ENDIF

          ENDIF

          IF ( .NOT. netcdf_extend )  THEN
!
!--          Create a new netCDF output file with requested netCDF format
             CALL netcdf_create_file( filename, id_set_3d(av), .TRUE., 37 )

!
!--          Define the header
             CALL netcdf_define_header( '3d', netcdf_extend, av )

!
!--          In case of parallel netCDF output, create flag file which tells
!--          combine_plot_fields that nothing is to do.
             IF ( myid == 0  .AND.  netcdf_data_format > 4 )  THEN
                OPEN( 99, FILE='NO_COMBINE_PLOT_FIELDS_3D' )
                WRITE ( 99, '(A)' )  'no combine_plot_fields.x neccessary'
                CLOSE( 99 )
             ENDIF

          ENDIF


       CASE ( 109 )
!
!--       Set filename
          filename = 'DATA_1D_PTS_NETCDF' // TRIM( coupling_char )

!
!--       Inquire, if there is a netCDF file from a previuos run. This should
!--       be opened for extension, if its variables match the actual run.
          INQUIRE( FILE=filename, EXIST=netcdf_extend )

          IF ( netcdf_extend )  THEN
!
!--          Open an existing netCDF file for output
             CALL netcdf_open_write_file( filename, id_set_pts, .FALSE., 393 )
!
!--          Read header information and set all ids. If there is a mismatch
!--          between the previuos and the actual run, netcdf_extend is returned
!--          as .FALSE.
             CALL netcdf_define_header( 'ps', netcdf_extend, 0 )

!
!--          Remove the local file, if it can not be extended
             IF ( .NOT. netcdf_extend )  THEN
                nc_stat = NF90_CLOSE( id_set_pts )
                CALL netcdf_handle_error( 'check_open', 394 )
                CALL local_system( 'rm ' // TRIM( filename ) )
             ENDIF

          ENDIF          

          IF ( .NOT. netcdf_extend )  THEN
!
!--          Create a new netCDF output file with requested netCDF format
             CALL netcdf_create_file( filename, id_set_pts, .FALSE., 395 )
!
!--          Define the header
             CALL netcdf_define_header( 'ps', netcdf_extend, 0 )

          ENDIF


#else


#endif

       CASE DEFAULT

          WRITE( message_string, * ) 'no OPEN-statement for file-id ',file_id
          CALL message( 'check_open', 'PA0172', 2, 2, -1, 6, 1 )

    END SELECT

!
!-- Set open flag
    openfile(file_id)%opened = .TRUE.

!
!-- Formats
3300 FORMAT ('#'/                                                              &
             'coord 1  file=',A,'  filetype=unformatted'/                      &
             'coord 2  file=',A,'  filetype=unformatted  skip=',I6/            &
             'coord 3  file=',A,'  filetype=unformatted  skip=',I6/            &
             '#')
4000 FORMAT ('# ',A)
5000 FORMAT ('# ',A/                                                           &
             '#1 E'/'#2 E*'/'#3 dt'/'#4 u*'/'#5 th*'/'#6 umax'/'#7 vmax'/      &
             '#8 wmax'/'#9 div_new'/'#10 div_old'/'#11 z_i_wpt'/'#12 z_i_pt'/  &
             '#13 w*'/'#14 w''pt''0'/'#15 w''pt'''/'#16 wpt'/'#17 pt(0)'/      &
             '#18 pt(zp)'/'#19 splptx'/'#20 splpty'/'#21 splptz')
8000 FORMAT (A/                                                                &
             '  step    time    # of parts     lPE sent/recv  rPE sent/recv  ',&
             'sPE sent/recv  nPE sent/recv    max # of parts  '/               &
             109('-'))

 END SUBROUTINE check_open
