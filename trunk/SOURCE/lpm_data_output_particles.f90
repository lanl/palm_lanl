!> @file lpm_data_output_particles.f90
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
! $Id: lpm_data_output_particles.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2123 2017-01-18 12:34:59Z hoffmann
!
! 2122 2017-01-18 12:22:54Z hoffmann
! Calculation of particle ID
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1822 2016-04-07 07:49:42Z hoffmann
! Tails removed. Unused variables removed.
!
! 1783 2016-03-06 18:36:17Z raasch
! name change of netcdf routines and module + related changes
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1359 2014-04-11 17:15:14Z hoffmann
! netCDF output currently not available
! output of particle data in binary format adopted to new particle structure
! 
! 1327 2014-03-21 11:00:16Z raasch
! -netcdf output queries
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! revision history before 2012 removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 849 2012-03-15 10:35:09Z raasch
! initial revision (former part of advec_particles)
!
! 22/02/12 - Initial version
!
! Description:
! ------------
!> Write particle data in FORTRAN binary and/or netCDF format 
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_data_output_particles
 

    USE control_parameters,                                                    &
        ONLY:  simulated_time

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE indices,                                                               &
        ONLY:  nxl, nxr, nyn, nys, nzb, nzt

    USE kinds

    USE particle_attributes,                                                   &
        ONLY:  grid_particles, number_of_particles,  particles, prt_count

    IMPLICIT NONE

    INTEGER(iwp) ::  ip !<
    INTEGER(iwp) ::  jp !<
    INTEGER(iwp) ::  kp !<

    CALL cpu_log( log_point_s(40), 'lpm_data_output', 'start' )

!
!-- Attention: change version number for unit 85 (in routine check_open)
!--            whenever the output format for this unit is changed!
    CALL check_open( 85 )

    WRITE ( 85 )  simulated_time
    WRITE ( 85 )  prt_count
          
    DO  ip = nxl, nxr
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt
             number_of_particles = prt_count(kp,jp,ip)
             particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
             IF ( number_of_particles <= 0 )  CYCLE
             WRITE ( 85 )  particles
          ENDDO
       ENDDO
    ENDDO

    CALL close_file( 85 )


#if defined( __netcdf )
! !
! !-- Output in netCDF format
!     CALL check_open( 108 )
! 
! !
! !-- Update the NetCDF time axis
!     prt_time_count = prt_time_count + 1
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_time_prt, &
!                             (/ simulated_time /),        &
!                             start = (/ prt_time_count /), count = (/ 1 /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 1 )
! 
! !
! !-- Output the real number of particles used
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_rnop_prt, &
!                             (/ number_of_particles /),   &
!                             start = (/ prt_time_count /), count = (/ 1 /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 2 )
! 
! !
! !-- Output all particle attributes
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(1), particles%age,      &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 3 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(2), particles%user,     &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 4 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(3), particles%origin_x, &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 5 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(4), particles%origin_y, &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 6 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(5), particles%origin_z, &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 7 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(6), particles%radius,   &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 8 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(7), particles%speed_x,  &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 9 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(8), particles%speed_y,  &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 10 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(9), particles%speed_z,  &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 11 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt,id_var_prt(10),                     &
!                             particles%weight_factor,                       &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 12 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(11), particles%x,       &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 13 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(12), particles%y,       &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 14 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(13), particles%z,       &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 15 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(14), particles%class,   &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 16 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(15), particles%group,   &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 17 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(16),                    &
!                             particles%id2,                                 &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 18 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(17), particles%id1,     &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 19 )
! 
#endif

    CALL cpu_log( log_point_s(40), 'lpm_data_output', 'stop' )

 END SUBROUTINE lpm_data_output_particles
