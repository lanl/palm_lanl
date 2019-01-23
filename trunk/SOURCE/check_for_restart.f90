!> @file check_for_restart.f90
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
! $Id: check_for_restart.f90 3049 2018-05-29 13:52:36Z Giersch $
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2298 2017-06-29 09:28:18Z raasch
! type of write_binary changed from CHARACTER to LOGICAL
! 
! 2101 2017-01-05 16:42:31Z suehring
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1797 2016-03-21 16:50:28Z raasch
! check accounts for nesting mode now
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1509 2014-12-16 08:56:46Z heinze
! bugfix: prevent infinite loop in case of automatic restarts
!
! 1468 2014-09-24 14:06:57Z maronga
! Added support for unscheduled job termination using the flag files
! DO_STOP_NOW and DO_RESTART_NOW
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1032 2012-10-21 13:03:21Z letzel
! minor reformatting
!
! Revision 1.1  1998/03/18 20:06:51  raasch
! Initial revision
!
!
! Description:
! ------------
!> Set stop flag, if restart is neccessary because of expiring cpu-time or
!> if it is forced by user
!------------------------------------------------------------------------------!
 SUBROUTINE check_for_restart
 

    USE control_parameters,                                                    &
        ONLY:  coupling_mode, dt_restart, end_time, message_string,            &
               run_description_header, simulated_time, terminate_coupled,      &
               terminate_coupled_remote, terminate_run,                        &
               termination_time_needed, time_restart,                          &
               time_since_reference_point, write_binary

    USE kinds

    USE pegrid

    IMPLICIT NONE

    INTEGER ::  global_communicator       !< global communicator to be used here

    LOGICAL ::  terminate_run_l           !<
    LOGICAL ::  do_stop_now = .FALSE.     !<
    LOGICAL ::  do_restart_now = .FALSE.  !<

    REAL(wp) ::  remaining_time !<


!
!-- Check remaining CPU-time
    CALL local_tremain( remaining_time )

!
!-- If necessary set a flag to stop the model run
    terminate_run_l = .FALSE.
    IF ( remaining_time <= termination_time_needed  .AND.  write_binary )  THEN

       terminate_run_l = .TRUE.
    ENDIF

!
!-- Set the global communicator to be used (depends on the mode in which PALM is
!-- running)
       global_communicator = comm2d

#if defined( __parallel )
!
!-- Make a logical OR for all processes. Stop the model run if at least
!-- one process has reached the time limit.
    IF ( collective_wait )  CALL MPI_BARRIER( global_communicator, ierr )
    CALL MPI_ALLREDUCE( terminate_run_l, terminate_run, 1, MPI_LOGICAL,     &
                        MPI_LOR, global_communicator, ierr )
#else
    terminate_run = terminate_run_l
#endif

!
!-- Output that job will be terminated
    IF ( terminate_run  .AND.  myid == 0 )  THEN
       WRITE( message_string, * ) 'run will be terminated because it is ',     &
                       'running out of job cpu limit & ',                      &
                       'remaining time:         ', remaining_time, ' s &',     &
                       'termination time needed:', termination_time_needed, ' s'
       CALL message( 'check_for_restart', 'PA0163', 0, 1, 0, 6, 0 )
    ENDIF

!
!-- In case of coupled runs inform the remote model of the termination 
!-- and its reason, provided the remote model has not already been 
!-- informed of another termination reason (terminate_coupled > 0) before, 
!-- or vice versa (terminate_coupled_remote > 0).
    IF ( terminate_run .AND. TRIM( coupling_mode ) /= 'uncoupled'  .AND.       &
         terminate_coupled == 0  .AND.  terminate_coupled_remote == 0 )  THEN

       terminate_coupled = 3

#if defined( __parallel )
       IF ( myid == 0 ) THEN
          CALL MPI_SENDRECV( terminate_coupled,        1, MPI_INTEGER,         &
                             target_id, 0,                                     &
                             terminate_coupled_remote, 1, MPI_INTEGER,         &
                             target_id, 0,                                     &
                             comm_inter, status, ierr )
       ENDIF
       CALL MPI_BCAST( terminate_coupled_remote, 1, MPI_INTEGER, 0, comm2d,    &
                       ierr )
#endif
    ENDIF


!
!-- Check if a flag file exists that forces a termination of the model
    IF ( myid == 0 )  THEN
       INQUIRE(FILE="DO_STOP_NOW", EXIST=do_stop_now)
       INQUIRE(FILE="DO_RESTART_NOW", EXIST=do_restart_now)

       IF ( do_stop_now .OR. do_restart_now )  THEN

          terminate_run_l = .TRUE.

          WRITE( message_string, * ) 'run will be terminated because user ',   &
                                  'forced a job finalization using a flag',    &
                                  'file:',                                     &
                                  '&DO_STOP_NOW: ', do_stop_now,               &
                                  '&DO_RESTART_NOW: ', do_restart_now 
          CALL message( 'check_for_restart', 'PA0398', 0, 0, 0, 6, 0 )

       ENDIF
    ENDIF


#if defined( __parallel )
!
!-- Make a logical OR for all processes. Stop the model run if a flag file has
!-- been detected above.
    IF ( collective_wait )  CALL MPI_BARRIER( global_communicator, ierr )
    CALL MPI_ALLREDUCE( terminate_run_l, terminate_run, 1, MPI_LOGICAL,        &
                        MPI_LOR, global_communicator, ierr )
#else
    terminate_run = terminate_run_l
#endif

!
!-- In case of coupled runs inform the remote model of the termination 
!-- and its reason, provided the remote model has not already been 
!-- informed of another termination reason (terminate_coupled > 0) before,
!-- or vice versa (terminate_coupled_remote > 0).
    IF ( terminate_run .AND. coupling_mode /= 'uncoupled' .AND.                &
         terminate_coupled == 0 .AND.  terminate_coupled_remote == 0 )  THEN

       terminate_coupled = 6

#if defined( __parallel )
       IF ( myid == 0 ) THEN
          CALL MPI_SENDRECV( terminate_coupled,        1, MPI_INTEGER,      &
                             target_id,  0,                                 &
                             terminate_coupled_remote, 1, MPI_INTEGER,      &
                             target_id,  0,                                 &
                             comm_inter, status, ierr )   
       ENDIF
       CALL MPI_BCAST( terminate_coupled_remote, 1, MPI_INTEGER, 0,         &
                       comm2d, ierr )  
#endif

    ENDIF

!
!-- Set the stop flag also, if restart is forced by user settings
    IF ( time_restart /= 9999999.9_wp  .AND.                                   &
         time_restart < time_since_reference_point )  THEN

!
!--    Restart is not neccessary, if the end time of the run (given by
!--    the user) has been reached
       IF ( simulated_time < end_time )  THEN
          terminate_run = .TRUE.
!
!--       Increment restart time, if forced by user, otherwise set restart
!--       time to default (no user restart)
          IF ( dt_restart /= 9999999.9_wp )  THEN
             time_restart = time_restart + dt_restart
          ELSE
             time_restart = 9999999.9_wp
          ENDIF

          WRITE( message_string, * ) 'run will be terminated due to user ',    &
                                  'settings of ',                              &
                                  'restart_time / dt_restart, ',               &
                                  'new restart time is: ', time_restart, ' s' 
          CALL message( 'check_for_restart', 'PA0164', 0, 0, 0, 6, 0 )
 
!
!--       In case of coupled runs inform the remote model of the termination 
!--       and its reason, provided the remote model has not already been 
!--       informed of another termination reason (terminate_coupled > 0) before,
!--       or vice versa (terminate_coupled_remote > 0).
          IF ( coupling_mode /= 'uncoupled' .AND. terminate_coupled == 0       &
               .AND.  terminate_coupled_remote == 0 )  THEN

             IF ( dt_restart /= 9999999.9_wp )  THEN
                terminate_coupled = 4
             ELSE
                terminate_coupled = 5
             ENDIF
#if defined( __parallel )
             IF ( myid == 0 ) THEN
                CALL MPI_SENDRECV( terminate_coupled,        1, MPI_INTEGER,   &
                                   target_id,  0,                              &
                                   terminate_coupled_remote, 1, MPI_INTEGER,   &
                                   target_id,  0,                              &
                                   comm_inter, status, ierr )   
             ENDIF
             CALL MPI_BCAST( terminate_coupled_remote, 1, MPI_INTEGER, 0,      &
                             comm2d, ierr )  
#endif
          ENDIF
       ELSE
          time_restart = 9999999.9_wp
       ENDIF
    ENDIF

!
!-- If the run is stopped, set a flag file which is necessary to initiate
!-- the start of a continuation run, except if the user forced to stop the
!-- run without restart
    IF ( terminate_run  .AND.  myid == 0  .AND.              &
         .NOT. do_stop_now)  THEN

       OPEN ( 90, FILE='CONTINUE_RUN', FORM='FORMATTED' )
       WRITE ( 90, '(A)' )  TRIM( run_description_header )
       CLOSE ( 90 )

    ENDIF


 END SUBROUTINE check_for_restart
