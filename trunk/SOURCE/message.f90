!> @file message.f90
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
! $Id: message.f90 3019 2018-05-13 07:05:43Z maronga $
! Temporaraly deactivate MPI_BARRIER as long as nested systems freeze due to
! asynchronous calls of location_message.
! 
! 2961 2018-04-12 10:13:48Z suehring
! Synchronize location message between parent and child. Message will be not 
! flushed before all models finished their respective task.  
! 
! 2932 2018-03-26 09:39:22Z maronga
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
! 1764 2016-02-28 12:45:19Z raasch
! nest id added to header string, add linefeed to stdout to get messages better
! seperated from the location messages,
! in case of nested runs, location messages are given only by the root domain
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1664 2015-09-23 06:18:12Z knoop
! updated information_string_2 to meet the new server structure
! 
! 1402 2014-05-09 14:25:13Z raasch
! formatting of messages modified
! 
! 1384 2014-05-02 14:31:06Z raasch
! routine location_message added
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 213 2008-11-13 10:26:18Z raasch
! Initial revision
!
! Description:
! ------------
!> Handling of the different kinds of messages.
!> Meaning of formal parameters:
!> requested_action: 0 - continue, 1 - abort by stop, 2 - abort by mpi_abort
!>                   3 - abort by mpi_abort using MPI_COMM_WORLD
!> message_level: 0 - informative, 1 - warning, 2 - error
!> output_on_pe: -1 - all, else - output on specified PE
!> file_id: 6 - stdout (*)
!> flush_file: 0 - no action, 1 - flush the respective output buffer
!------------------------------------------------------------------------------!
 SUBROUTINE message( routine_name, message_identifier, requested_action, &
                     message_level, output_on_pe, file_id, flush_file )
 
    USE control_parameters,                                                    &
        ONLY:  abort_mode, message_string

    USE kinds

    USE pegrid

    IMPLICIT NONE

    CHARACTER(LEN=6)   ::  message_identifier            !<
    CHARACTER(LEN=20)  ::  nest_string                   !< nest id information
    CHARACTER(LEN=*)   ::  routine_name                  !<
    CHARACTER(LEN=200) ::  header_string                 !<
    CHARACTER(LEN=200) ::  information_string_1          !<
    CHARACTER(LEN=200) ::  information_string_2          !<

    INTEGER(iwp) ::  file_id                             !<
    INTEGER(iwp) ::  flush_file                          !<
    INTEGER(iwp) ::  i                                   !<
    INTEGER(iwp) ::  message_level                       !<
    INTEGER(iwp) ::  output_on_pe                        !<
    INTEGER(iwp) ::  requested_action                    !<

    LOGICAL ::  do_output                                !<
    LOGICAL ::  pe_out_of_range                          !<


    do_output       = .FALSE.
    pe_out_of_range = .FALSE.

!
       nest_string = ''
!
!-- Create the complete output string, starting with the message level
    IF ( message_level == 0 )  THEN
       header_string = '--- informative message' // TRIM(nest_string) //       &
                       ' ---  ID:'
    ELSEIF ( message_level == 1 )  THEN
       header_string = '+++ warning message' // TRIM(nest_string) // ' ---  ID:'
    ELSEIF ( message_level == 2 )  THEN
       header_string = '+++ error message' // TRIM(nest_string) // ' ---  ID:'
    ELSE
       WRITE( header_string,'(A,I2)' )  '+++ unknown message level' //         &
                                        TRIM(nest_string) // ': ',             &
                                        message_level
    ENDIF

!
!-- Add the message identifier and the generating routine
    header_string = TRIM( header_string ) // ' ' // message_identifier // &
                    '   generated by routine: ' // TRIM( routine_name )
  
    information_string_1 = 'Further information can be found at'
    IF(message_identifier(1:2) == 'NC') THEN
       information_string_2 = 'http://palm.muk.uni-hannover.de/trac/wiki/doc' // &
                              '/app/errmsg#NC'
    ELSE
       information_string_2 = 'http://palm.muk.uni-hannover.de/trac/wiki/doc' // &
                              '/app/errmsg#' // message_identifier
    ENDIF
    

!
!-- Output the output string and the corresponding message string which had
!-- been already assigned in the calling subroutine.
!
!-- First find out if output shall be done on this PE.
    IF ( output_on_pe == -1 )  THEN
       do_output = .TRUE.
    ELSEIF ( myid == output_on_pe )  THEN
       do_output = .TRUE.
    ENDIF
#if defined( __parallel )
!
!-- In case of illegal pe number output on pe0
    IF ( output_on_pe > numprocs-1 )  THEN
       pe_out_of_range = .TRUE.
       IF ( myid == 0 )  do_output = .TRUE.
    ENDIF
#endif

!
!-- Now do the output
    IF ( do_output )  THEN

       IF ( file_id == 6 )  THEN
!
!--       Output on stdout
          WRITE( *, '(//A/)' )  TRIM( header_string )
!
!--       Cut message string into pieces and output one piece per line.
!--       Remove leading blanks.
          message_string = ADJUSTL( message_string )
          i = INDEX( message_string, '&' )
          DO WHILE ( i /= 0 )
             WRITE( *, '(4X,A)' )  ADJUSTL( message_string(1:i-1) )
             message_string = ADJUSTL( message_string(i+1:) )
             i = INDEX( message_string, '&' )
          ENDDO
          WRITE( *, '(4X,A)' )  TRIM( message_string )
          WRITE( *, '(4X,A)' )  ''
          WRITE( *, '(4X,A)' )  TRIM( information_string_1 ) 
          WRITE( *, '(4X,A)' )  TRIM( information_string_2 ) 
          WRITE( *, '(4X,A)' )  ''

       ELSE
!
!--       Output on requested file id (file must have been opened elsewhere!)
          WRITE( file_id, '(A/)' )  TRIM( header_string )
!
!--       Cut message string into pieces and output one piece per line.
!--       Remove leading blanks.
          message_string = ADJUSTL( message_string )
          i = INDEX( message_string, '&' )
          DO WHILE ( i /= 0 )
             WRITE( file_id, '(4X,A)' )  ADJUSTL( message_string(1:i-1) )
             message_string = ADJUSTL( message_string(i+1:) )
             i = INDEX( message_string, '&' )
          ENDDO
          WRITE( file_id, '(4X,A)' )  TRIM( message_string )
          WRITE( file_id, '(4X,A)' )  ''
          WRITE( file_id, '(4X,A)' )  TRIM( information_string_1 ) 
          WRITE( file_id, '(4X,A)' )  TRIM( information_string_2 ) 
          WRITE( file_id, '(4X,A)' )  ''
!
!--       Flush buffer, if requested
          IF ( flush_file == 1 )  FLUSH( file_id )
       ENDIF

       IF ( pe_out_of_range )  THEN
          WRITE ( *, '(A)' )  '+++ WARNING from routine message:'
          WRITE ( *, '(A,I6,A)' )  '    PE ', output_on_pe, &
                                   ' choosed for output is larger '
          WRITE ( *, '(A,I6)' )  '    than the maximum number of used PEs', &
                                 numprocs-1
          WRITE ( *, '(A)' )  '    Output is done on PE0 instead'
       ENDIF

    ENDIF

!
!-- Abort execution, if requested
    IF ( requested_action > 0 )  THEN
       abort_mode = requested_action
       CALL local_stop
    ENDIF

 END SUBROUTINE message


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Prints out the given location on stdout
!------------------------------------------------------------------------------!
 
 SUBROUTINE location_message( location, advance )


    USE, INTRINSIC ::  ISO_FORTRAN_ENV,                                        &
        ONLY:  OUTPUT_UNIT

    USE pegrid

    IMPLICIT NONE

    CHARACTER(LEN=*) ::  location !< text to be output on stdout
    LOGICAL          ::  advance  !< switch for advancing/noadvancing I/O
!
!-- Output for nested runs only on the root domain

    IF ( myid == 0 )  THEN
       IF ( advance )  THEN
          WRITE ( OUTPUT_UNIT, '(6X,''--- '',A)' )  TRIM( location )
       ELSE
          WRITE ( OUTPUT_UNIT, '(6X,''... '',A)', ADVANCE='NO' )               &
                TRIM( location )
       ENDIF
       FLUSH( OUTPUT_UNIT )
    ENDIF

 END SUBROUTINE location_message
