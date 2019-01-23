!> @file run_control.f90
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
! $Id: run_control.f90 2718 2018-01-02 08:49:38Z maronga $
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
! 1808 2016-04-05 19:44:00Z raasch
! routine local_flush replaced by FORTRAN statement
!
! 1697 2015-10-28 17:14:10Z raasch
! small E- and F-FORMAT changes to avoid informative compiler messages about
! insufficient field width
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
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
! 1318 2014-03-17 13:35:16Z raasch
! module interfaces removed
!
! 1257 2013-11-08 15:18:40Z raasch
! output format for theta* changed to avoid output of *****
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1001 2012-09-13 14:08:46Z raasch
! all actions concerning leapfrog scheme removed
!
! Revision 1.1  1997/08/11 06:25:38  raasch
! Initial revision
!
!
! Description:
! ------------
!> Computation and output of run-control quantities
!------------------------------------------------------------------------------!
 SUBROUTINE run_control
 

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point

    USE control_parameters,                                                    &
        ONLY:  advected_distance_x, advected_distance_y,                       &
               current_timestep_number, disturbance_created, dt_3d, mgcycles,  &
               run_control_header, runnr, simulated_time, simulated_time_chr,  &
               timestep_reason

    USE indices,                                                               &
        ONLY:  nzb

    USE kinds

    USE pegrid

    USE statistics,                                                            &
        ONLY:  flow_statistics_called, hom, pr_palm, u_max, u_max_ijk, v_max,  &
               v_max_ijk, w_max, w_max_ijk

    IMPLICIT NONE

    CHARACTER (LEN=1) ::  disturb_chr

!
!-- If required, do statistics
    IF ( .NOT. flow_statistics_called ) THEN
           CALL flow_statistics
    endif
!
!-- Flow_statistics has its own cpu-time measurement
    CALL cpu_log( log_point(11), 'run_control', 'start' )

!
!-- Output
    IF ( myid == 0 )  THEN

!
!--    Check, whether file unit is already open (may have been opened in header
!--    before)
       CALL check_open( 15 )

!
!--    If required, write header
!       IF ( .NOT. run_control_header )  THEN
!          WRITE ( 15, 100 )
!          run_control_header = .TRUE.
!       ENDIF

!
!--    If required, set disturbance flag
       IF ( disturbance_created )  THEN
          disturb_chr = 'D'
       ELSE
          disturb_chr = ' '
       ENDIF
!       WRITE ( 15, 101 )  runnr, current_timestep_number, simulated_time_chr,  &
!                          INT( ( simulated_time-INT( simulated_time ) ) * 100),&
!                          dt_3d, timestep_reason, u_max, disturb_chr,          &
!                          v_max, disturb_chr, w_max, hom(nzb,1,pr_palm,0),     &
!                          hom(nzb+8,1,pr_palm,0), hom(nzb+3,1,pr_palm,0),      &
!                          hom(nzb+6,1,pr_palm,0), hom(nzb+4,1,pr_palm,0),      &
!                          hom(nzb+5,1,pr_palm,0), hom(nzb+9,1,pr_palm,0),      &
!                          hom(nzb+10,1,pr_palm,0), u_max_ijk(1:3),             &
!                          v_max_ijk(1:3), w_max_ijk(1:3),                      &
!                          advected_distance_x/1000.0_wp,                       &
!                          advected_distance_y/1000.0_wp, mgcycles
!
!--    Write buffer contents to disc immediately
!       FLUSH( 15 )

    ENDIF
!
!-- If required, reset disturbance flag. This has to be done outside the above 
!-- IF-loop, because the flag would otherwise only be reset on PE0
    IF ( disturbance_created )  disturbance_created = .FALSE.

    CALL cpu_log( log_point(11), 'run_control', 'stop' )

!
!-- Formats
100 FORMAT (///'Run-control output:'/ &
              &'------------------'// &
          &'RUN  ITER. HH:MM:SS.SS    DT(E)     UMAX     VMAX     WMAX     U', &
          &'*    W*      THETA*     Z_I     ENERG.   DISTENERG    DIVOLD    ', &
          &' DIVNEW     UMAX(KJI)    VMAX(KJI)    WMAX(KJI)   ADVECX   ADVEC', &
          &'Y   MGCYC'/                                                        &
          &'----------------------------------------------------------------', &
          &'----------------------------------------------------------------', &
          &'----------------------------------------------------------------', &
          &'---------')
101 FORMAT (I3,1X,I6,1X,A8,'.',I2.2,1X,F8.4,A1,1X,F8.4,A1,F8.4,A1,F8.4,1X,     &
            F6.3,1X,F5.2, &
            2X,E10.3,2X,F6.0,1X,4(E10.3,1X),3(3(I4),1X),F8.3,1X,F8.3,5X,I3)

 END SUBROUTINE run_control
