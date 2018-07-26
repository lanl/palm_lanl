!> @file src/control.f90
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
! Copyright 2017-2018 Leibniz Universitaet Hannover
! Copyright 2017-2018 Deutscher Wetterdienst Offenbach
!------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: control.f90 2718 2018-01-02 08:49:38Z maronga $
! Initial revision
!
! 
!
! Authors:
! --------
! @author Eckhard Kadasch
!
! Description:
! ------------
!> The control module provides routines for timing INIFOR and writing runtime
!> feedback to the terminal and a log file.
!------------------------------------------------------------------------------!
 MODULE control

    USE defs,                                                                  &
        ONLY:  LNAME, dp

    USE util,                                                                  &
        ONLY:  real_to_str, real_to_str_f

    IMPLICIT NONE

    CHARACTER (LEN=5000) ::  message = ''

 CONTAINS

    SUBROUTINE report(routine, message)

       CHARACTER(LEN=*), INTENT(IN) ::  routine
       CHARACTER(LEN=*), INTENT(IN) ::  message
       INTEGER                      ::  u
       LOGICAL, SAVE                ::  is_first_run = .TRUE.

       PRINT *, "inifor: " // TRIM(message) // "  [ " // TRIM(routine) // " ]"

       IF (is_first_run)  THEN
          OPEN( NEWUNIT=u, FILE='inifor.log', STATUS='replace' )
          is_first_run = .FALSE.
       ELSE
          OPEN( NEWUNIT=u, FILE='inifor.log', POSITION='append', STATUS='old' )
       END IF
          
       WRITE(u, *)  TRIM(message) // "  [ " // TRIM(routine) // " ]"

       CLOSE(u)

    END SUBROUTINE report


    SUBROUTINE warn(routine, message)

       CHARACTER(LEN=*), INTENT(IN) ::  routine
       CHARACTER(LEN=*), INTENT(IN) ::  message

       CALL report(routine, "WARNING: " // TRIM(message))

    END SUBROUTINE warn


    SUBROUTINE abort(routine, message)

       CHARACTER(LEN=*), INTENT(IN) ::  routine
       CHARACTER(LEN=*), INTENT(IN) ::  message

       CALL report(routine, "ERROR: " // TRIM(message) // " Stopping.")
       STOP

    END SUBROUTINE abort


    SUBROUTINE run_control(mode, budget)

       CHARACTER(LEN=*), INTENT(IN) ::  mode, budget
       REAL(dp), SAVE               ::  t0, t1
       REAL(dp), SAVE               ::  t_comp=0.0_dp, &
                                        t_alloc=0.0_dp, &
                                        t_init=0.0_dp, &
                                        t_read=0.0_dp, &
                                        t_total=0.0_dp, &
                                        t_write=0.0_dp
       CHARACTER(LEN=*), PARAMETER  ::  fmt='(F6.2)'


       SELECT CASE(TRIM(mode))

       CASE('init')
          CALL CPU_TIME(t0)

       CASE('time')

          CALL CPU_TIME(t1)

          SELECT CASE(TRIM(budget))

             CASE('alloc')
                t_alloc = t_alloc + t1 - t0

             CASE('init')
                t_init = t_init + t1 - t0

             CASE('read')
                t_read = t_read + t1 - t0

             CASE('write')
                t_write = t_write + t1 - t0

             CASE('comp')
                t_comp = t_comp + t1 - t0

             CASE DEFAULT
                CALL abort('run_control', "Time Budget '" // TRIM(mode) // "' is not supported.")

          END SELECT

          t0 = t1

       CASE('report')
           t_total = t_init + t_read + t_write + t_comp

           CALL report('run_control', " *** CPU time ***")

           CALL report('run_control', "Initialization: " // real_to_str(t_init)  // &
                       " s (" // TRIM(real_to_str(100*t_init/t_total, fmt))      // " %)")

           CALL report('run_control', "(De-)Allocation:" // real_to_str(t_alloc)  // &
                       " s (" // TRIM(real_to_str(100*t_alloc/t_total, fmt))      // " %)")

           CALL report('run_control', "Reading data:   " // real_to_str(t_read)  // &
                       " s (" // TRIM(real_to_str(100*t_read/t_total, fmt))      // " %)")

           CALL report('run_control', "Writing data:   " // real_to_str(t_write) // &
                       " s (" // TRIM(real_to_str(100*t_write/t_total, fmt))     // " %)")

           CALL report('run_control', "Computation:    " // real_to_str(t_comp)  // &
                       " s (" // TRIM(real_to_str(100*t_comp/t_total, fmt))      // " %)")

           CALL report('run_control', "Total:          " // real_to_str(t_total) // &
                       " s (" // TRIM(real_to_str(100*t_total/t_total, fmt))     // " %)")

       CASE DEFAULT
          CALL abort('run_control', "Mode '" // TRIM(mode) // "' is not supported.")

       END SELECT

    END SUBROUTINE run_control

 END MODULE

