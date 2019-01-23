!> @file cpulog_mod.f90
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
! $Id: cpulog_mod.f90 3049 2018-05-29 13:52:36Z Giersch $
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
! 
! 3028 2018-05-23 10:27:23Z knoop
! renamed d3par to runtime_parameters
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2266 2017-06-09 09:27:21Z raasch
! bugfix: cpu-time per grid point is now calculated using ngp_3d in order to
! avoid integer overflow,
! some re-formatting
! 
! 2119 2017-01-17 16:51:50Z raasch
!
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC relevant code removed
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1931 2016-06-10 12:06:59Z suehring
! Adjustment in character length and format statement
! 
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
! 
!
! 1808 2016-04-05 19:44:00Z raasch
! cpu measurements are done with standard FORTRAN routine on every machine
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1402 2014-05-09 14:25:13Z raasch
! location messages added
!
! 1369 2014-04-24 05:57:38Z raasch
! routine description added
!
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
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
! former files/routines cpu_log and cpu_statistics combined to one module,
! which also includes the former data module cpulog from the modules-file
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! Revision 1.1  1997/07/24 11:12:29  raasch
! Initial revision
!
!
! Description:
! ------------
!> CPU-time measurements for any program part whatever. Results of the
!> measurements are output at the end of the run in local file CPU_MEASURES.
!>
!> To measure the CPU-time (better to say the wallclock time) of a specific code
!> segment, two calls of cpu_log have to be used as brackets in front and at the
!> end of the segment:
!>
!>     CALL cpu_log( log_point(n), 'any identifier', 'start' )
!>       ... code segment ...
!>     CALL cpu_log( log_point(n), 'any identifier', 'stop' )
!>
!> Parts of the code segment can be excluded from the measurement by additional
!> call of cpu_log:
!>
!>       ... first segment to be measured
!>     CALL cpu_log( log_point(n), 'any identifier', 'pause' )
!>       ... oart of segment to be excluded from measurement
!>     CALL cpu_log( log_point(n), 'any identifier', 'continue' )
!>       ... second segment to be mesasured
!>
!> n is an INTEGER within the interval [1,100] defining the id of the specific
!> code segment, 'any identifier' is a string describing the code segment to be
!> measured. It can be freely chosen and results will appear under this name in
!> file CPU_MEASURES. ids can only be used once. If you like to do a
!> measurement of a new part of the code, please look for an id which is unused
!> ao far.
!>
!> runtime_parameters-parameter cpu_log_barrierwait can be used to set an MPI 
!> barrier at the beginning of the measurement (modus 'start' or 'continue'), 
!> to avoid that idle times (due to MPI calls in the code segment, which are
!> waiting for other processes to be finished) affect the measurements.
!> If barriers shall not be used at all, a fourth, optional parameter has to be
!> given:
!>
!>     CALL cpu_log( ..., ..., 'start', cpu_log_nowait )
!>
!> Variable log_point should be used for non-overlapping code segments, and they
!> should sum up to the total cpu-time required by the complete run.
!> Variable log_point_s can be used for any other special (s) measurements.
!------------------------------------------------------------------------------!
 MODULE cpulog
 

    USE control_parameters,                                                    &
        ONLY: message_string, nr_timesteps_this_run, run_description_header,   &
              synchronous_exchange
               
    USE indices,                                                               &
        ONLY: ngp_3d, nx, ny, nz
        
    USE kinds
    
    USE pegrid

    IMPLICIT NONE

    PRIVATE
    PUBLIC   cpu_log, cpu_log_barrierwait, cpu_log_nowait, cpu_statistics,     &
             initial_wallclock_time, log_point, log_point_s

    INTERFACE cpu_log
       MODULE PROCEDURE cpu_log
    END INTERFACE cpu_log

    INTERFACE cpu_statistics
       MODULE PROCEDURE cpu_statistics
    END INTERFACE cpu_statistics

    INTEGER(iwp), PARAMETER ::  cpu_log_continue = 0  !< 
    INTEGER(iwp), PARAMETER ::  cpu_log_pause = 1     !< 
    INTEGER(iwp), PARAMETER ::  cpu_log_start = 2     !< 
    INTEGER(iwp), PARAMETER ::  cpu_log_stop = 3      !< 

    LOGICAL            ::  cpu_log_barrierwait = .FALSE.  !< 
    LOGICAL, PARAMETER ::  cpu_log_nowait = .FALSE.       !< 

    REAL(wp) ::  initial_wallclock_time  !< 

    TYPE logpoint
       REAL(wp)           ::  isum       !< 
       REAL(wp)           ::  ivect      !< 
       REAL(wp)           ::  mean       !< 
       REAL(wp)           ::  mtime      !< 
       REAL(wp)           ::  mtimevec   !< 
       REAL(wp)           ::  sum        !< 
       REAL(wp)           ::  vector     !< 
       INTEGER(iwp)       ::  counts     !< 
       CHARACTER (LEN=25) ::  place      !< 
    END TYPE logpoint

    TYPE(logpoint), DIMENSION(100) ::  log_point = logpoint( 0.0_wp, 0.0_wp,   &
                                       0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                                       0, ' ' )
    TYPE(logpoint), DIMENSION(100) ::  log_point_s = logpoint( 0.0_wp, 0.0_wp, &
                                       0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                                       0, ' ' )

    SAVE

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE cpu_log( log_event, place, modus, barrierwait )

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  modus              !< 
       CHARACTER (LEN=*) ::  place              !< 
       
       LOGICAL           ::  wait_allowed       !< 
       LOGICAL, OPTIONAL ::  barrierwait        !< 
       LOGICAL, SAVE     ::  first = .TRUE.     !< 
       
       REAL(wp)          ::  mtime = 0.0_wp     !< 
       REAL(wp)          ::  mtimevec = 0.0_wp  !< 
       TYPE(logpoint)    ::  log_event          !< 

       INTEGER(idp)     ::  count        !< 
       INTEGER(idp)     ::  count_rate   !< 


!
!--    Initialize and check, respectively, point of measurement
       IF ( log_event%place == ' ' )  THEN
          log_event%place = place
       ELSEIF ( log_event%place /= place )  THEN
          WRITE( message_string, * ) 'wrong argument expected: ',              &
                            TRIM(log_event%place), ' given: ',  TRIM( place )
          CALL message( 'cpu_log', 'PA0174', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    Determine, if barriers are allowed to set
       IF ( PRESENT( barrierwait ) )  THEN
          wait_allowed = barrierwait
       ELSE
          wait_allowed = .TRUE.
       ENDIF

!
!--    MPI barrier, if requested, in order to avoid measuring wait times
!--    caused by MPI routines waiting for other MPI routines of other
!--    PEs that have not yet finished
#if defined( __parallel )
       IF ( cpu_log_barrierwait  .AND.  wait_allowed  .AND.                    &
            ( modus == 'start'  .OR.  modus == 'continue' ) )  THEN
          CALL MPI_BARRIER( comm2d, ierr )
       ENDIF
#endif

!
!--    Take current time
       CALL SYSTEM_CLOCK( count, count_rate )
       mtime = REAL( count, KIND=wp ) / REAL( count_rate, KIND=wp )

!
!--    Start, stop or pause measurement
       IF ( modus == 'start'  .OR.  modus == 'continue' )  THEN
          log_event%mtime    = mtime
          log_event%mtimevec = mtimevec
       ELSEIF ( modus == 'pause' )  THEN
          IF ( ( mtime - log_event%mtime ) < 0.0  .AND.  first )  THEN
             WRITE( message_string, * ) 'negative time interval occured',      &
                         '&PE',myid,' L=PAUSE "',TRIM(log_event%place),        &
                         '" new=', mtime,' last=',log_event%mtime
             CALL message( 'cpu_log', 'PA0176', 0, 1, -1, 6, 0 )
             first = .FALSE.
          ENDIF
          log_event%isum     = log_event%isum + mtime - log_event%mtime
          log_event%ivect    = log_event%ivect  + mtimevec - log_event%mtimevec
       ELSEIF ( modus == 'stop' )  THEN
          IF ( ( mtime - log_event%mtime + log_event%isum ) < 0.0  .AND.       &
               first )  THEN
             WRITE( message_string, * ) 'negative time interval occured',      &
                       '&PE',myid,' L=STOP "',TRIM(log_event%place),'" new=',  &
                       mtime,' last=',log_event%mtime,' isum=',log_event%isum
             CALL message( 'cpu_log', 'PA0177', 0, 1, -1, 6, 0 )
             first = .FALSE.
          ENDIF
          log_event%mtime    = mtime    - log_event%mtime    + log_event%isum
          log_event%mtimevec = mtimevec - log_event%mtimevec + log_event%ivect
          log_event%sum      = log_event%sum  + log_event%mtime
          IF ( log_event%sum < 0.0  .AND.  first )  THEN
             WRITE( message_string, * ) 'negative time interval occured',      &
                       '&PE',myid,' L=STOP "',TRIM(log_event%place),'" sum=',  &
                       log_event%sum,' mtime=',log_event%mtime
             CALL message( 'cpu_log', 'PA0178', 0, 1, -1, 6, 0 )
             first = .FALSE.
          ENDIF
          log_event%vector   = log_event%vector + log_event%mtimevec
          log_event%counts   = log_event%counts + 1
          log_event%isum     = 0.0_wp
          log_event%ivect    = 0.0_wp
       ELSE
          message_string = 'unknown modus of time measurement: ' //            &
                           TRIM( modus )
          CALL message( 'cpu_log', 'PA0179', 0, 1, -1, 6, 0 )
       ENDIF

    END SUBROUTINE cpu_log


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Analysis and output of the cpu-times measured. All PE results are collected
!> on PE0 in order to calculate the mean cpu-time over all PEs and other
!> statistics. The output is sorted according to the amount of cpu-time consumed
!> and output on PE0.
!------------------------------------------------------------------------------!
 
    SUBROUTINE cpu_statistics

       IMPLICIT NONE

       INTEGER(iwp)    ::  i               !< 
       INTEGER(iwp)    ::  ii(1)           !< 
       INTEGER(iwp)    ::  iii             !< 
       INTEGER(iwp)    ::  sender          !< 
       REAL(wp)       ::  average_cputime  !< 
       REAL(wp), SAVE ::  norm = 1.0_wp    !< 
       REAL(wp), DIMENSION(:),   ALLOCATABLE ::  pe_max        !< 
       REAL(wp), DIMENSION(:),   ALLOCATABLE ::  pe_min        !< 
       REAL(wp), DIMENSION(:),   ALLOCATABLE ::  pe_rms        !< 
       REAL(wp), DIMENSION(:),   ALLOCATABLE ::  sum           !< 
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pe_log_points !< 


       CALL location_message( 'calculating cpu statistics', .FALSE. )

!
!--    Compute cpu-times in seconds
       log_point%mtime  = log_point%mtime  / norm
       log_point%sum    = log_point%sum    / norm
       log_point%vector = log_point%vector / norm
       WHERE ( log_point%counts /= 0 )
          log_point%mean = log_point%sum / log_point%counts
       END WHERE


!
!--    Collect cpu-times from all PEs and calculate statistics
       IF ( myid == 0 )  THEN
!
!--       Allocate and initialize temporary arrays needed for statistics
          ALLOCATE( pe_max( SIZE( log_point ) ), pe_min( SIZE( log_point ) ),  &
                    pe_rms( SIZE( log_point ) ),                               &
                    pe_log_points( SIZE( log_point ), 0:numprocs-1 ) )
          pe_min = log_point%sum
          pe_max = log_point%sum    ! need to be set in case of 1 PE
          pe_rms = 0.0_wp

#if defined( __parallel )
!
!--       Receive data from all PEs
          DO  i = 1, numprocs-1
             CALL MPI_RECV( pe_max(1), SIZE( log_point ), MPI_REAL,            &
                            i, i, comm2d, status, ierr )
             sender = status(MPI_SOURCE)
             pe_log_points(:,sender) = pe_max
          ENDDO
          pe_log_points(:,0) = log_point%sum   ! Results from PE0
!
!--       Calculate mean of all PEs, store it on log_point%sum
!--       and find minimum and maximum
          DO  iii = 1, SIZE( log_point )
             DO  i = 1, numprocs-1
                log_point(iii)%sum = log_point(iii)%sum + pe_log_points(iii,i)
                pe_min(iii) = MIN( pe_min(iii), pe_log_points(iii,i) )
                pe_max(iii) = MAX( pe_max(iii), pe_log_points(iii,i) )
             ENDDO
             log_point(iii)%sum = log_point(iii)%sum / numprocs
!
!--          Calculate rms
             DO  i = 0, numprocs-1
                pe_rms(iii) = pe_rms(iii) + (                                  &
                                    pe_log_points(iii,i) - log_point(iii)%sum  &
                                            )**2
             ENDDO
             pe_rms(iii) = SQRT( pe_rms(iii) / numprocs )
          ENDDO
       ELSE
!
!--       Send data to PE0 (pe_max is used as temporary storage to send
!--       the data in order to avoid sending the data type log)
          ALLOCATE( pe_max( SIZE( log_point ) ) )
          pe_max = log_point%sum
          CALL MPI_SEND( pe_max(1), SIZE( log_point ), MPI_REAL, 0, myid,      &
                         comm2d, ierr )
#endif

       ENDIF

!
!--    Write cpu-times
       IF ( myid == 0 )  THEN
!
!--       Re-store sums
          ALLOCATE( sum( SIZE( log_point ) ) )
          WHERE ( log_point%counts /= 0 )
             sum = log_point%sum
          ELSEWHERE
             sum = -1.0_wp
          ENDWHERE

!
!--       Get total time in order to calculate CPU-time per gridpoint and
!--       timestep
          IF ( nr_timesteps_this_run /= 0 )  THEN
             average_cputime = log_point(1)%sum / REAL( ngp_3d(0), KIND=wp ) / &
                               REAL( nr_timesteps_this_run, KIND=wp ) * 1E6_wp
                               ! in micro-sec
          ELSE
             average_cputime = -1.0_wp
          ENDIF

!
!--       Write cpu-times sorted by size
          CALL check_open( 18 )
#if defined( __parallel )
          WRITE ( 18, 100 )  TRIM( run_description_header ),                   &
                             numprocs * threads_per_task, pdims(1), pdims(2),  &
                             threads_per_task, nx+1, ny+1, nz,                 &
                             nr_timesteps_this_run, average_cputime

          WRITE ( 18, 110 )
#else
          WRITE ( 18, 100 )  TRIM( run_description_header ),                   &
                             numprocs * threads_per_task, 1, 1,                &
                             threads_per_task, nx+1, ny+1, nz,                 &
                             nr_timesteps_this_run, average_cputime

          WRITE ( 18, 110 )
#endif
          DO
             ii = MAXLOC( sum )
             i = ii(1)
             IF ( sum(i) /= -1.0_wp )  THEN
                WRITE ( 18, 102 ) &
              log_point(i)%place, log_point(i)%sum,                            &
                   log_point(i)%sum / log_point(1)%sum * 100.0_wp,             &
                   log_point(i)%counts, pe_min(i), pe_max(i), pe_rms(i)
                sum(i) = -1.0_wp
             ELSE
                EXIT
             ENDIF
          ENDDO
       ENDIF


!
!--    The same procedure again for the individual measurements.
!
!--    Compute cpu-times in seconds
       log_point_s%mtime  = log_point_s%mtime  / norm
       log_point_s%sum    = log_point_s%sum    / norm
       log_point_s%vector = log_point_s%vector / norm
       WHERE ( log_point_s%counts /= 0 )
          log_point_s%mean = log_point_s%sum / log_point_s%counts
       END WHERE

!
!--    Collect cpu-times from all PEs and calculate statistics
#if defined( __parallel )
!
!--    Set barrier in order to avoid that PE0 receives log_point_s-data
!--    while still busy with receiving log_point-data (see above)
       CALL MPI_BARRIER( comm2d, ierr )
#endif
       IF ( myid == 0 )  THEN
!
!--       Initialize temporary arrays needed for statistics
          pe_min = log_point_s%sum
          pe_max = log_point_s%sum    ! need to be set in case of 1 PE
          pe_rms = 0.0_wp

#if defined( __parallel )
!
!--       Receive data from all PEs
          DO  i = 1, numprocs-1
             CALL MPI_RECV( pe_max(1), SIZE( log_point ), MPI_REAL, &
                            MPI_ANY_SOURCE, MPI_ANY_TAG, comm2d, status, ierr )
             sender = status(MPI_SOURCE)
             pe_log_points(:,sender) = pe_max
          ENDDO
          pe_log_points(:,0) = log_point_s%sum   ! Results from PE0
!
!--       Calculate mean of all PEs, store it on log_point_s%sum
!--       and find minimum and maximum
          DO  iii = 1, SIZE( log_point )
             DO  i = 1, numprocs-1
                log_point_s(iii)%sum = log_point_s(iii)%sum +                  &
                                       pe_log_points(iii,i)
                pe_min(iii) = MIN( pe_min(iii), pe_log_points(iii,i) )
                pe_max(iii) = MAX( pe_max(iii), pe_log_points(iii,i) )
             ENDDO
             log_point_s(iii)%sum = log_point_s(iii)%sum / numprocs
!
!--          Calculate rms
             DO  i = 0, numprocs-1
                pe_rms(iii) = pe_rms(iii) + (                                  &
                                   pe_log_points(iii,i) - log_point_s(iii)%sum &
                                            )**2
             ENDDO
             pe_rms(iii) = SQRT( pe_rms(iii) / numprocs )
          ENDDO
       ELSE
!
!--       Send data to PE0 (pe_max is used as temporary storage to send
!--       the data in order to avoid sending the data type log)
          pe_max = log_point_s%sum
          CALL MPI_SEND( pe_max(1), SIZE( log_point ), MPI_REAL, 0, 0, comm2d, &
                         ierr )
#endif

       ENDIF

!
!--    Write cpu-times
       IF ( myid == 0 )  THEN
!
!--       Re-store sums
          WHERE ( log_point_s%counts /= 0 )
             sum = log_point_s%sum
          ELSEWHERE
             sum = -1.0_wp
          ENDWHERE

!
!--       Write cpu-times sorted by size
          WRITE ( 18, 101 )
          DO
             ii = MAXLOC( sum )
             i = ii(1)
             IF ( sum(i) /= -1.0_wp )  THEN
                WRITE ( 18, 102 )                                              &
                   log_point_s(i)%place, log_point_s(i)%sum,                   &
                   log_point_s(i)%sum / log_point(1)%sum * 100.0_wp,           &
                   log_point_s(i)%counts, pe_min(i), pe_max(i), pe_rms(i)
                sum(i) = -1.0_wp
             ELSE
                EXIT
             ENDIF
          ENDDO

!
!--       Output of handling of MPI operations
          IF ( collective_wait )  THEN
             WRITE ( 18, 103 )
          ELSE
             WRITE ( 18, 104 )
          ENDIF
          IF ( cpu_log_barrierwait )  WRITE ( 18, 111 )
          IF ( synchronous_exchange )  THEN
             WRITE ( 18, 105 )
          ELSE
             WRITE ( 18, 106 )
          ENDIF

!
!--       Empty lines in order to create a gap to the results of the model
!--       continuation runs
          WRITE ( 18, 107 )

!
!--       Unit 18 is not needed anymore
          CALL close_file( 18 )

       ENDIF

       DEALLOCATE( sum, pe_max, pe_min, pe_rms, pe_log_points )

       CALL location_message( 'finished', .TRUE. )

   100 FORMAT (A/11('-')//'CPU measures for ',I5,' PEs (',I5,'(x) * ',I5,'(y', &
               &') tasks *',I5,' threads):'//                                  &
               'gridpoints (x/y/z): ',20X,I5,' * ',I5,' * ',I5/                &
               'nr of timesteps: ',22X,I6/                                     &
               'cpu time per grid point and timestep: ',5X,F8.5,' * 10**-6 s')

   101 FORMAT (/'special measures:'/ &
               &'-----------------------------------------------------------', &
               &'------------------------')

   102 FORMAT (A25,2X,F9.3,2X,F7.2,1X,I7,3(1X,F9.3))
   103 FORMAT (/'Barriers are set in front of collective operations')
   104 FORMAT (/'No barriers are set in front of collective operations')
   105 FORMAT (/'Exchange of ghostpoints via MPI_SENDRCV')
   106 FORMAT (/'Exchange of ghostpoints via MPI_ISEND/MPI_IRECV')
   107 FORMAT (//)
   110 FORMAT ('-------------------------------------------------------------',     &
               &'---------'//&
               &'place:                              mean        counts      min  ',&
               &'     max       rms'/ &
               &'                                sec.      %                sec. ', &
               &'     sec.      sec.'/  &
               &'-----------------------------------------------------------',      &
               &'------------------------')
   111 FORMAT (/'Barriers are set at beginning (start/continue) of measurements')

    END SUBROUTINE cpu_statistics

 END MODULE cpulog
