!> @file data_output_spectra.f90
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
! ------------------
! 
! 
! Former revisions:
! -----------------
! $Id: data_output_spectra.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2101 2017-01-05 16:42:31Z suehring
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1960 2016-07-12 16:34:24Z suehring
! Additional default spectra for passive scalar
! 
! 1833 2016-04-07 14:23:03Z raasch
! spectrum renamed spectra_mod, spectra related variables moved to spectra_mod,
! routines data_output_spectra_x/y removed
!
! 1786 2016-03-08 05:49:27Z raasch
! cpp-directives for spectra removed, immediate return if no spectra levels are
! given
!
! 1783 2016-03-06 18:36:17Z raasch
! name change of netcdf routines and module + related changes
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
! 
! 1327 2014-03-21 11:00:16Z raasch
! -netcdf output queries
!
! 1324 2014-03-21 09:13:16Z suehring
! Bugfix: module statistics and module spectrum added, missing variables in ONLY
! arguments added
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL functions provided with KIND-attribute
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1318 2014-03-17 13:35:16Z raasch
! module interfaces removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 964 2012-07-26 09:14:24Z raasch
! code for profil-output removed
!
! Revision 1.1  2001/01/05 15:14:20  raasch
! Initial revision
!
!
! Description:
! ------------
!> Writing spectra data on file, using a special format which allows
!> plotting of these data with PROFIL-graphic-software
!------------------------------------------------------------------------------!
 SUBROUTINE data_output_spectra
 
#if defined( __netcdf )
    USE control_parameters,                                                    &
        ONLY:  message_string, run_description_header,                         &
               time_since_reference_point

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point

    USE kinds

    USE NETCDF

    USE netcdf_interface,                                                      &
        ONLY:  id_set_sp, id_var_time_sp, nc_stat, netcdf_handle_error

    USE pegrid

    USE spectra_mod,                                                           &
        ONLY:  average_count_sp, averaging_interval_sp, comp_spectra_level,    &
               data_output_sp, dosp_time_count, spectra_direction, spectrum_x, &
               spectrum_y


    IMPLICIT NONE

    INTEGER(iwp) ::  cranz_x !<
    INTEGER(iwp) ::  cranz_y !<
    INTEGER(iwp) ::  m       !<
    INTEGER(iwp) ::  pr      !<
    
    LOGICAL      ::  frame_x !< 
    LOGICAL      ::  frame_y !<

    CALL cpu_log( log_point(31), 'data_output_spectra', 'start' )

!
!-- Check if user gave any levels for spectra to be calculated
    IF ( comp_spectra_level(1) == 999999 )  RETURN

!
!-- Output is only performed on PE0
    IF ( myid == 0 )  THEN

!
!--    Open file for spectra output in NetCDF format
       CALL check_open( 107 )

!
!--    Increment the counter for number of output times
       dosp_time_count = dosp_time_count + 1

!
!--    Update the spectra time axis
       nc_stat = NF90_PUT_VAR( id_set_sp, id_var_time_sp,        &
                               (/ time_since_reference_point /), &
                               start = (/ dosp_time_count /), count = (/ 1 /) )
       CALL netcdf_handle_error( 'data_output_spectra', 47 )

!
!--    If necessary, calculate time average and reset average counter
       IF ( average_count_sp == 0 )  THEN
           message_string = 'no spectra data available'
           CALL message( 'data_output_spectra', 'PA0186', 0, 0, 0, 6, 0 )
       ENDIF
       IF ( average_count_sp /= 1 )  THEN
          spectrum_x = spectrum_x / REAL( average_count_sp, KIND=wp )
          spectrum_y = spectrum_y / REAL( average_count_sp, KIND=wp )
          average_count_sp = 0
       ENDIF

!
!--    Loop over all spectra defined by the user
       m = 1
       DO WHILE ( data_output_sp(m) /= ' '  .AND.  m <= 10 )

          SELECT CASE ( TRIM( data_output_sp(m) ) )

             CASE ( 'u' )
                pr = 1

             CASE ( 'v' )
                pr = 2

             CASE ( 'w' )
                pr = 3

             CASE ( 'pt' )
                pr = 4

             CASE ( 'q' )
                pr = 5

             CASE ( 's' )
                pr = 6

             CASE DEFAULT
!
!--             The DEFAULT case is reached either if the parameter 
!--             data_output_sp(m) contains a wrong character string or if the 
!--             user has coded a special case in the user interface. There, the 
!--             subroutine user_spectra checks which of these two conditions 
!--             applies.
                CALL user_spectra( 'data_output', m, pr )

          END SELECT

!
!--       Output of spectra in NetCDF format
!--       Output of x-spectra
          IF ( INDEX( spectra_direction(m), 'x' ) /= 0 ) THEN
             CALL output_spectra_netcdf( m, 'x' )
          ENDIF
!
!--       Output of y-spectra
          IF ( INDEX( spectra_direction(m), 'y' ) /= 0 ) THEN
             CALL output_spectra_netcdf( m, 'y' )
          ENDIF

!
!--       Increase counter for next spectrum
          m = m + 1

       ENDDO

!
!--    Reset spectra values
       spectrum_x = 0.0_wp; spectrum_y = 0.0_wp

    ENDIF

    CALL cpu_log( log_point(31), 'data_output_spectra', 'stop' )

#if defined( __parallel )
!    CALL MPI_BARRIER( comm2d, ierr )  ! really necessary
#endif

#endif
 END SUBROUTINE data_output_spectra


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
 SUBROUTINE output_spectra_netcdf( nsp, direction )
#if defined( __netcdf )

    USE constants,                                                             &
        ONLY:  pi

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE indices,                                                               &
        ONLY:  nx, ny

    USE kinds

    USE NETCDF

    USE netcdf_interface,                                                      &
        ONLY:  id_set_sp, id_var_dospx, id_var_dospy, nc_stat,                 &
               netcdf_handle_error

    USE spectra_mod,                                                           &
        ONLY:  dosp_time_count, n_sp_x, n_sp_y, spectrum_x, spectrum_y


    IMPLICIT NONE

    CHARACTER (LEN=1), INTENT(IN) ::  direction     !<

    INTEGER(iwp), INTENT(IN)      ::  nsp           !<

    INTEGER(iwp)                  ::  i             !<
    INTEGER(iwp)                  ::  k             !<

    REAL(wp)                      ::  frequency     !<

    REAL(wp), DIMENSION(nx/2)     ::  netcdf_data_x !<
    REAL(wp), DIMENSION(ny/2)     ::  netcdf_data_y !<


    IF ( direction == 'x' )  THEN

       DO  k = 1, n_sp_x

          DO  i = 1, nx/2
             frequency = 2.0_wp * pi * i / ( dx * ( nx + 1 ) )
             netcdf_data_x(i) = frequency * spectrum_x(i,k,nsp)
          ENDDO

          nc_stat = NF90_PUT_VAR( id_set_sp, id_var_dospx(nsp), netcdf_data_x, &
                                  start = (/ 1, k, dosp_time_count /), &
                                  count = (/ nx/2, 1, 1 /) )
          CALL netcdf_handle_error( 'data_output_spectra', 348 )

       ENDDO

    ENDIF

    IF ( direction == 'y' )  THEN

       DO  k = 1, n_sp_y

          DO  i = 1, ny/2
             frequency = 2.0_wp * pi * i / ( dy * ( ny + 1 ) )
             netcdf_data_y(i) = frequency * spectrum_y(i,k,nsp)
          ENDDO

          nc_stat = NF90_PUT_VAR( id_set_sp, id_var_dospy(nsp), netcdf_data_y, &
                                  start = (/ 1, k, dosp_time_count /), &
                                  count = (/ ny/2, 1, 1 /) )
          CALL netcdf_handle_error( 'data_output_spectra', 349 )

       ENDDO

    ENDIF

#endif
 END SUBROUTINE output_spectra_netcdf
