 PROGRAM compare_palm_logs

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
! $Id: compare_palm_logs.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 1310 2014-03-14 08:01:56Z raasch
! update of GPL copyright
!
! 1046 2012-11-09 14:38:45Z maronga
! code put under GPL (PALM 3.9)
! name of data-directories are read from input
!
! Description:
! ------------
! This routine compares the log files from two different PALM runs.
!
! This routine must be compiled with:
! decalpha:
!    f95 -cpp -fast -r8
! IBM-Regatta:
!    xlf95 -qsuffix=cpp=f90 -qrealsize=8 -q64 -qmaxmem=-1 -Q -O3
! IMUK:
!    ifort compare...f90 -o compare...x
!    -cpp -axW -r8 -nbs -Vaxlib
! NEC-SX6:
!    sxf90 compare...f90 -o compare...x 
!    -C hopt -Wf '-A idbl4'
!------------------------------------------------------------------------------!

    IMPLICIT NONE

!
!-- Local variables
    CHARACTER (LEN=5) ::  id_char
    CHARACTER (LEN=80), DIMENSION(2)  ::  directory, log_message
    CHARACTER (LEN=100), DIMENSION(2) ::  filename

    INTEGER ::  count=0, i, id, i1(2), i2(2), j, j1(2), j2(2), k, k1(2), k2(2), &
                n_err, n_files(2)

    LOGICAL ::  found

    REAL    ::  simtime(2)

    INTEGER, DIMENSION(:,:),   ALLOCATABLE ::  array_2d_i_1, array_2d_i_2

    REAL, DIMENSION(:,:),   ALLOCATABLE ::  array_2d_1, array_2d_2
    REAL, DIMENSION(:,:,:), ALLOCATABLE ::  array_1, array_2


!
!-- Read the two data-directories to be compared
    PRINT*, '*** please enter name of first data-directory:'
    READ ( *, * )  directory(1)
    directory(1) = TRIM( directory(1) ) // '/'

    PRINT*, '*** please enter name of second data-directory:'
    READ ( *, * )  directory(2)
    directory(2) = TRIM( directory(2) ) // '/'

!
!-- Check, if file from PE0 exists on directory 1. Stop, if it does not exist.
    n_files(1) = 0

    WRITE (id_char,'(''_'',I4.4)')  n_files(1)
    INQUIRE ( FILE=TRIM( directory(1) )//id_char, EXIST=found )
!
!-- Find out the number of files (equal to the number of PEs which
!-- have been used in PALM) and open them
    DO  WHILE ( found )

       OPEN ( n_files(1)+100, FILE=TRIM( directory(1) )//id_char, &
              FORM='UNFORMATTED' )
       n_files(1) = n_files(1) + 1
       WRITE (id_char,'(''_'',I4.4)')  n_files(1)
       INQUIRE ( FILE=TRIM( directory(1) )//id_char, EXIST=found )

    ENDDO

    IF ( n_files(1) == 0 )  THEN
       PRINT*, '+++ no file _0000 in directory "', TRIM( directory(1) ), '"'
       STOP
    ELSE
       PRINT*, '*** directory "', TRIM( directory(1) ), '": ', n_files(1), &
               ' files found'
    ENDIF

!
!-- Same for the second directory
    n_files(2) = 0

    WRITE (id_char,'(''_'',I4.4)')  n_files(2)
    INQUIRE ( FILE=TRIM( directory(2) )//id_char, EXIST=found )

    DO  WHILE ( found )

       OPEN ( n_files(2)+200, FILE=TRIM( directory(2) )//id_char, &
              FORM='UNFORMATTED' )
       n_files(2) = n_files(2) + 1
       WRITE (id_char,'(''_'',I4.4)')  n_files(2)
       INQUIRE ( FILE=TRIM( directory(2) )//id_char, EXIST=found )

    ENDDO

!
!-- Number of files must be identical
    IF ( n_files(1) /= n_files(2) )  THEN
       PRINT*, '+++ file number mismatch'
       PRINT*, '    ', TRIM( directory(1) ), ': ', n_files(1), ' files'
       PRINT*, '    ', TRIM( directory(2) ), ': ', n_files(2), ' files'
       STOP
    ENDIF

!
!-- Compare the data file by file
    DO  id = 0, n_files(1)-1

       count = 0

       WRITE (filename(1),'(A,''_'',I4.4)')  TRIM( directory(1) ), id
       WRITE (filename(2),'(A,''_'',I4.4)')  TRIM( directory(2) ), id

       PRINT*, '*** comparing files "', TRIM( filename(1) ),'" "', &
               TRIM( filename(2) ), '"'
       DO
          PRINT*,' '
          READ ( id+100, END=100 )  log_message(1)
          PRINT*,'    --- ', TRIM( log_message(1) )
          READ ( id+200, END=900 )  log_message(2)

          IF ( TRIM( log_message(1) ) /= TRIM( log_message(2) ) )  THEN
             PRINT*,'    +++ log message on file 2 does not match:'
             PRINT*,'        ', TRIM( log_message(2) )
          ENDIF

          count = count + 1
          IF ( log_message(1)(1:2) == '3d' )  THEN
             PRINT*,'    *** reading 3d array'
             READ ( id+100, END=901 )  simtime(1), i1(1), i2(1), j1(1), j2(1), &
                                       k1(1), k2(1)
             PRINT*,'        time=', simtime(1)
             PRINT*,'        array size=(',i1(1),':',i2(1), &
                                       ',',j1(1),':',j2(1),',',k1(1),':',k2(1),')'
             READ ( id+200, END=902 )  simtime(2), i1(2), i2(2), j1(2), j2(2), &
                                       k1(2), k2(2)
             IF ( simtime(1) /= simtime(2) .OR. i1(1) /= i1(2) .OR. &
                  i2(1) /= i2(2) .OR. j1(1) /= j1(2) .OR. j2(1) /= j2(2) .OR. &
                  k1(1) /= k1(2) .OR. k2(1) /= k2(2) )  THEN
                PRINT*,'    +++ time/indices on file 2 does not match:'
                PRINT*,'        time=', simtime(2)
                PRINT*,'        array size=(',i1(2),':', &
                                i2(2), ',',j1(2),':',j2(2),',',k1(2),':',k2(2),')'
                STOP
             ENDIF

             ALLOCATE( array_1(i1(1):i2(1),j1(1):j2(1),k1(1):k2(1)), &
                       array_2(i1(2):i2(2),j1(2):j2(2),k1(2):k2(2)) )

             READ ( id+100, END=903 )  array_1
             READ ( id+200, END=904 )  array_2

             n_err = 0
loop:        DO  k = k1(1), k2(1)
loop1:           DO  j = j1(1), j2(1)
                   DO  i = i1(1), i2(1)
                      IF ( array_1(i,j,k) /= array_2(i,j,k) )  THEN
                         PRINT*,'+++ data mismatch on element (',i,',',j,',',k,')'
                         PRINT*,'    array_1: ', array_1(i,j,k)
                         PRINT*,'    array_2: ', array_2(i,j,k)
                         n_err = n_err + 1
                         IF ( n_err > 5 )  EXIT loop
                      ENDIF
                   ENDDO
                ENDDO loop1
             ENDDO loop

             DEALLOCATE( array_1, array_2 )

          ELSEIF ( log_message(1)(1:2) == '2d' )  THEN
             PRINT*,'    *** reading 2d array'
             READ ( id+100, END=901 )  simtime(1), i1(1), i2(1), j1(1), j2(1)
             PRINT*,'        time=', simtime(1)
             PRINT*,'        array size=(',i1(1),':',i2(1), &
                                       ',',j1(1),':',j2(1),')'
             READ ( id+200, END=902 )  simtime(2), i1(2), i2(2), j1(2), j2(2)
             IF ( simtime(1) /= simtime(2) .OR. i1(1) /= i1(2) .OR. &
                  i2(1) /= i2(2) .OR. j1(1) /= j1(2) .OR. j2(1) /= j2(2) )  THEN
                PRINT*,'    +++ time/indices on file 2 does not match:'
                PRINT*,'        time=', simtime(2)
                PRINT*,'        array size=(',i1(2),':', &
                                i2(2), ',',j1(2),':',j2(2),')'
             ENDIF

             ALLOCATE( array_2d_1(i1(1):i2(1),j1(1):j2(1)), &
                       array_2d_2(i1(2):i2(2),j1(2):j2(2)) )

             READ ( id+100, END=903 )  array_2d_1
             READ ( id+200, END=904 )  array_2d_2

             IF ( i1(1) /= i1(2) )  i1(1) = i1(2)
             IF ( i2(1) /= i2(2) )  i2(1) = i2(2)
             IF ( j1(1) /= j1(2) )  j1(1) = j1(2)
             IF ( j2(1) /= j2(2) )  j2(1) = j2(2)

             n_err = 0
loop2:       DO  j = j1(1), j2(1)
                DO  i = i1(1), i2(1)
                   IF ( array_2d_1(i,j) /= array_2d_2(i,j) )  THEN
                      PRINT*,'+++ data mismatch on element (',i,',',j,')'
                      PRINT*,'    array_1: ', array_2d_1(i,j)
                      PRINT*,'    array_2: ', array_2d_2(i,j)
                      n_err = n_err + 1
                      IF ( n_err > 5 )  EXIT loop2
                   ENDIF
                ENDDO
             ENDDO loop2

             DEALLOCATE( array_2d_1, array_2d_2 )

          ELSE
             PRINT*,'    *** reading 2d int array'
             READ ( id+100, END=901 )  simtime(1), i1(1), i2(1), j1(1), j2(1)
             PRINT*,'        time=', simtime(1)
             PRINT*,'        array size=(',i1(1),':',i2(1), &
                                       ',',j1(1),':',j2(1),')'
             READ ( id+200, END=902 )  simtime(2), i1(2), i2(2), j1(2), j2(2)
             IF ( simtime(1) /= simtime(2) .OR. i1(1) /= i1(2) .OR. &
                  i2(1) /= i2(2) .OR. j1(1) /= j1(2) .OR. j2(1) /= j2(2) )  THEN
                PRINT*,'    +++ time/indices on file 2 does not match:'
                PRINT*,'        time=', simtime(2)
                PRINT*,'        array size=(',i1(2),':', &
                                i2(2), ',',j1(2),':',j2(2),')'
             ENDIF

             ALLOCATE( array_2d_i_1(i1(1):i2(1),j1(1):j2(1)), &
                       array_2d_i_2(i1(2):i2(2),j1(2):j2(2)) )

             READ ( id+100, END=903 )  array_2d_i_1
             READ ( id+200, END=904 )  array_2d_i_2

             IF ( i1(1) /= i1(2) )  i1(1) = i1(2)
             IF ( i2(1) /= i2(2) )  i2(1) = i2(2)
             IF ( j1(1) /= j1(2) )  j1(1) = j1(2)
             IF ( j2(1) /= j2(2) )  j2(1) = j2(2)

             n_err = 0
loop3:       DO  j = j1(1), j2(1)
                DO  i = i1(1), i2(1)
                   IF ( array_2d_i_1(i,j) /= array_2d_i_2(i,j) )  THEN
                      PRINT*,'+++ data mismatch on element (',i,',',j,')'
                      PRINT*,'    array_1: ', array_2d_i_1(i,j)
                      PRINT*,'    array_2: ', array_2d_i_2(i,j)
                      n_err = n_err + 1
                      IF ( n_err > 5 )  EXIT loop3
                   ENDIF
                ENDDO
             ENDDO loop3

             DEALLOCATE( array_2d_i_1, array_2d_i_2 )

          ENDIF

!          IF ( count > 8 )  STOP
       ENDDO

100    PRINT*, '*** end of data on file "', TRIM( filename(1) ), '"'
       PRINT*, '*** files seem to be identical'
       PRINT*, ' '
    ENDDO

    STOP

900 PRINT*,'+++ unexpected end on file "', TRIM( filename(2) ), '"'
    STOP
901 PRINT*,'+++ unexpected end on file "', TRIM( filename(1) ), '"'
    PRINT*,'    while reading indices'
    STOP
902 PRINT*,'+++ unexpected end on file "', TRIM( filename(2) ), '"'
    PRINT*,'    while reading indices'
    STOP
903 PRINT*,'+++ unexpected end on file "', TRIM( filename(1) ), '"'
    PRINT*,'    while reading array data'
    STOP
904 PRINT*,'+++ unexpected end on file "', TRIM( filename(2) ), '"'
    PRINT*,'    while reading array data'
    STOP

 END PROGRAM compare_palm_logs



