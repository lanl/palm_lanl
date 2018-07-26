!> @file data_output_tseries.f90
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
! $Id: data_output_tseries.f90 2718 2018-01-02 08:49:38Z maronga $
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
! 1783 2016-03-06 18:36:17Z raasch
! name change of netcdf routines and module + related changes
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1524 2015-01-14 13:18:19Z keck
! Bugfix: increment dots_time_count after the call of subroutine check_open
! 
! 1327 2014-03-21 11:00:16Z raasch
! -netcdf output queries
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
! barrier argument removed from cpu_log.
! module interfaces removed
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! Revision 1.1  1998/03/03 08:00:13  raasch
! Initial revision
!
!
! Description:
! ------------
!> Time series output for PROFIL. Always all time series are stored. A selection
!> can be applied via the PROFIL-parameters in close_file.
!------------------------------------------------------------------------------!
 SUBROUTINE data_output_tseries
 

    USE control_parameters,                                                    &
        ONLY:  dots_time_count, time_since_reference_point

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point  

    USE kinds

#if defined( __netcdf )
    USE NETCDF
#endif
    USE netcdf_interface,                                                      &
        ONLY:  dots_num, id_set_ts, id_var_dots, id_var_time_ts, nc_stat,      &
               netcdf_handle_error

    USE pegrid

    USE profil_parameter
    
    USE statistics,                                                            &
        ONLY:  flow_statistics_called, statistic_regions, ts_value

    IMPLICIT NONE


    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  sr !<


!
!-- If required, compute statistics.
    IF ( .NOT. flow_statistics_called )  CALL flow_statistics

!
!-- Flow_statistics has its own cpu-time measuring.
    CALL cpu_log( log_point(21), 'data_output_tseries', 'start' )

    IF ( myid == 0 )  THEN

!
!--    Open file for time series output in NetCDF format
       CALL check_open( 105 )
       
!--    Increment the counter for number of output times
!      CAUTION: The following line has to be after the call of the subroutine
!               check_open, since check_open resets the counter dots_time_count
!               to 0, if a new file is opened
       dots_time_count = dots_time_count + 1
       
#if defined( __netcdf )
!
!--    Update the time series time axis
       nc_stat = NF90_PUT_VAR( id_set_ts, id_var_time_ts,        &
                               (/ time_since_reference_point /), &
                               start = (/ dots_time_count /),    &
                               count = (/ 1 /) )
       CALL netcdf_handle_error( 'data_output_tseries', 350 )
#endif

!
!--    Time series output for the total domain (and each subregion, if
!--    applicable)
       DO  sr = 0, statistic_regions

#if defined( __netcdf )
          DO  i = 1, dots_num
             nc_stat = NF90_PUT_VAR( id_set_ts, id_var_dots(i,sr),  &
                                     (/ ts_value(i,sr) /),          &
                                     start = (/ dots_time_count /), &
                                     count = (/ 1 /) )
             CALL netcdf_handle_error( 'data_output_tseries', 351 )
          ENDDO
#endif

       ENDDO

    ENDIF


    CALL cpu_log( log_point(21), 'data_output_tseries', 'stop' )

!
!-- formats
500 FORMAT (23(E15.7,1X))

 END SUBROUTINE data_output_tseries
