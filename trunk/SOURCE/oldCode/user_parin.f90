!> @file user_parin.f90
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
! $Id: user_parin.f90 3049 2018-05-29 13:52:36Z Giersch $
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! renamed userpar to user_parameters
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2512 2017-10-04 08:26:59Z raasch
! current interface revision number number set to r2512
! 
! 2298 2017-06-29 09:28:18Z raasch
! user interface current revision updated
! 
! 2101 2017-01-05 16:42:31Z suehring
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1833 2016-04-07 14:23:03Z raasch
! required interface revision changed
!
! 1783 2016-03-06 18:36:17Z raasch
! required interface revision changed
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
!
! 1668 2015-09-23 13:45:36Z raasch
! current interface revision number number set to r1663
!
! 1666 2015-09-23 07:31:10Z raasch
! interface revision number is set to blank
!
! 1320 2014-03-20 08:40:49Z raasch
! kind-parameters added to all INTEGER and REAL declaration statements, 
! kinds are defined in new module kinds, 
! old module precision_kind is removed,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 931 2012-06-08 15:09:28Z maronga
! Re-enabled check for max_pr_user
!
! 841 2012-02-28 12:29:49Z maronga
! Bugfix: disable max_pr_user check during prior namelist file check
!
! 217 2008-12-09 18:00:48Z letzel
! +topography_grid_convention
! Former file user_interface.f90 split into one file per subroutine
!
! Description:
! ------------
!> Interface to read user-defined namelist-parameters.
!------------------------------------------------------------------------------!
 SUBROUTINE user_parin
 

    USE control_parameters
    
    USE kinds
    
    USE pegrid
    
    USE statistics
    
    USE user

    IMPLICIT NONE

    CHARACTER (LEN=80) ::  line   !< 

    INTEGER(iwp) ::  i                 !< 
    INTEGER(iwp) ::  j                 !< 
    INTEGER(iwp) ::  max_pr_user_tmp   !< 


    NAMELIST /userpar/  data_output_pr_user, data_output_user, region,         &
                        data_output_masks_user
                        
                        
    NAMELIST /user_parameters/  data_output_pr_user, data_output_user, region, &
                        data_output_masks_user

!
!-- Set revision number of this default interface version. It will be checked within
!-- the main program (palm). Please change the revision number in case that the
!-- current revision does not match with previous revisions (e.g. if routines
!-- have been added/deleted or if parameter lists in subroutines have been changed).
    user_interface_current_revision = 'r2512'

!
!-- Position the namelist-file at the beginning (it was already opened in
!-- parin), search for user-defined namelist-group ("userpar", but any other
!-- name can be choosed) and position the file at this line.
    REWIND ( 11 )

    line = ' '
    DO   WHILE ( INDEX( line, '&user_parameters' ) == 0 )
       READ ( 11, '(A)', END=10 )  line
    ENDDO
    BACKSPACE ( 11 )

!
!-- Read user-defined namelist
    READ ( 11, userpar )

    user_defined_namelist_found = .TRUE.

    GOTO 12
    
    
10  REWIND ( 11 )

    line = ' ' 
    DO   WHILE ( INDEX( line, '&userpar' ) == 0 )
       READ ( 11, '(A)', END=12 )  line
    ENDDO
    BACKSPACE ( 11 )

!
!-- Read user-defined namelist
    READ ( 11, userpar )
    
    
    message_string = 'namelist userpar is deprecated and will be ' //          &
                     'removed in near future. &Please use namelist ' //        &
                     'user_parameters instead' 
    CALL message( 'user_parin', 'PA0487', 0, 1, 0, 6, 0 )
       
    user_defined_namelist_found = .TRUE.
    
    
 12 CONTINUE

!
!-- Determine the number of user-defined profiles and append them to the
!-- standard data output (data_output_pr)
    IF ( user_defined_namelist_found )  THEN
       max_pr_user_tmp = 0
       IF ( data_output_pr_user(1) /= ' ' )  THEN
          i = 1
          DO  WHILE ( data_output_pr(i) /= ' '  .AND.  i <= 100 )
             i = i + 1
          ENDDO
          j = 1
          DO  WHILE ( data_output_pr_user(j) /= ' '  .AND.  j <= 100 )
             data_output_pr(i) = data_output_pr_user(j)
             max_pr_user_tmp   = max_pr_user_tmp + 1
             i = i + 1
             j = j + 1
          ENDDO
      ENDIF


!
!--    In case of a restart run, the number of user-defined profiles on the
!--    restart file (already stored in max_pr_user) has to match the one given
!--    for the current run
       IF ( TRIM( initializing_actions ) == 'read_restart_data' )  THEN
          IF ( max_pr_user /= max_pr_user_tmp )  THEN
             WRITE( message_string, * ) 'the number of user-defined profiles ',&
                     'given in data_output_pr (', max_pr_user_tmp, ') doe',    &
                     'snot match the one ',                                    &
                     'found in the restart file (', max_pr_user,               &
                                     ')'
             CALL message( 'user_parin', 'UI0009', 1, 2, 0, 6, 0 )
          ENDIF
       ELSE
          max_pr_user = max_pr_user_tmp
       ENDIF

    ENDIF
 
    RETURN

 END SUBROUTINE user_parin

