!> @file local_stop.f90
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
! $Id: local_stop.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2371 2017-08-24 13:01:17Z kanani
! Removed unnecessary USE of vertical_nesting_mod
! 
! 2365 2017-08-21 14:59:59Z kanani
! Vertical nesting implemented (SadiqHuq)
! 
! 2101 2017-01-05 16:42:31Z suehring
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1804 2016-04-05 16:30:18Z maronga
! Removed code for parameter file check (__check)
!
! 1764 2016-02-28 12:45:19Z raasch
! abort with MPI_COMM_WORLD added, nested runs always abort with MPI_ABORT
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! revision history before 2012 removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 809 2012-01-30 13:32:58Z maronga
! Bugfix: replaced .AND. and .NOT. with && and ! in the preprocessor directives
!
! 807 2012-01-25 11:53:51Z maronga
! New cpp directive "__check" implemented which is used by check_namelist_files
!
! Revision 1.1  2002/12/19 15:46:23  raasch
! Initial revision
!
!
! Description:
! ------------
!> Stop program execution
!------------------------------------------------------------------------------!
 SUBROUTINE local_stop
 

    USE control_parameters,                                                    &
        ONLY:  abort_mode, coupling_mode, coupling_mode_remote, dt_restart,    &
               stop_dt, terminate_coupled, terminate_coupled_remote,           &
               terminate_run, time_restart

    USE pegrid

#if defined( __parallel )
         IF ( abort_mode == 1 )  THEN
             CALL MPI_FINALIZE( ierr )
             STOP
          ELSEIF ( abort_mode == 2 )  THEN
             CALL MPI_ABORT( comm2d, 9999, ierr )
          ELSEIF ( abort_mode == 3 )  THEN
             CALL MPI_ABORT( MPI_COMM_WORLD, 9999, ierr )
          ENDIF
#else

    STOP

#endif

 END SUBROUTINE local_stop    
