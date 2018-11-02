!> @file init_coupling.f90
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
! 2018-11-01 cbegeman
! Read coupling parameters from runfile rather than the standard input
! 
! Former revisions:
! ------------------
! $Id: init_coupling.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2669 2017-12-06 16:03:27Z raasch
! file extension for vertical nesting changed to "_NV"
! 
! 2365 2017-08-21 14:59:59Z kanani
! Vertical nesting implemented (SadiqHuq)
! 
! 2298 2017-06-29 09:28:18Z raasch
! MPI2 coupling removed
! 
! 2101 2017-01-05 16:42:31Z suehring
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1808 2016-04-05 19:44:00Z raasch
! routine local_getenv replaced by standard FORTRAN routine
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
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
! 222 2009-01-12 16:04:16Z letzel
! Initial revision
!
! Description:
! ------------
!> Initializing coupling via MPI-1 or MPI-2 if the coupled version of PALM is
!> called.
!------------------------------------------------------------------------------!
  SUBROUTINE init_coupling
 

    USE control_parameters,                                                    &
        ONLY:  coupling_char, coupling_mode
        
    USE kinds
    
    USE pegrid

    USE vertical_nesting_mod

    IMPLICIT NONE

!
!-- Local variables
    INTEGER(iwp) ::  i            !<
    INTEGER(iwp) ::  inter_color  !<
    
    INTEGER(iwp), DIMENSION(:) ::  bc_data(0:3) = 0  !<

!
!-- Get information about the coupling mode from runfile and
!-- distribute it to the other PEs. Distribute PEs to 2 new communicators.
!-- ATTENTION: numprocs will be reset according to the new communicators
#if defined ( __parallel )

    IF ( myid == 0 )  THEN
       CALL check_open( 50 )
       READ( 50,*,ERR=10,END=10 ) coupling_mode, bc_data(1), bc_data(2)
10     CONTINUE
       IF ( TRIM( coupling_mode ) == 'coupled_run' )  THEN
          i = 1
       ELSEIF ( TRIM( coupling_mode ) == 'vnested_twi' )  THEN
          i = 9
       ELSE
          i = 0
       ENDIF
       bc_data(0) = i

!
!--    Check if '_O' has to be used as file extension in an uncoupled ocean
!--    run. This is required, if this run shall be continued as a coupled run.
       IF ( TRIM( coupling_mode ) == 'precursor_ocean' )  bc_data(3) = 1

    ENDIF

    CALL MPI_BCAST( bc_data(0), 4, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
    i = bc_data(0)

    IF ( i == 0 )  THEN
       coupling_mode = 'uncoupled'
!
!--    In case of a precursor ocean run, an additional flag file is created.
!--    This is necessary for data_output_2d_on_each_pe = .T.
       IF ( bc_data(3) == 1 )  THEN
          OPEN( 90, FILE='PRECURSOR_OCEAN', FORM='FORMATTED' )
          WRITE ( 90, '(''TRUE'')' )
          CLOSE ( 90 )
       ENDIF
    ELSEIF ( i == 9 )  THEN

!
!--    Set a flag to identify runs with vertical nesting
       vnested = .TRUE.
       
       comm_inter = MPI_COMM_WORLD
       
!
!--    Split the total available PE's into two groups
!--    numprocs for Coarse and Fine Grid are specified via mrun argument -N 
       IF ( myid < bc_data(1) )  THEN
          inter_color     = 0
          numprocs        = bc_data(1)
          coupling_mode   = 'vnested_crse'
       ELSE
          inter_color     = 1
          numprocs        = bc_data(2)
          coupling_mode   = 'vnested_fine'
       ENDIF
       
       CALL MPI_COMM_SPLIT( MPI_COMM_WORLD, inter_color, 0, comm_palm, ierr )
       comm2d = comm_palm
       
       OPEN( 90, FILE='VNESTING_PORT_OPENED', FORM='FORMATTED' )
       WRITE ( 90, '(''TRUE'')' )
       CLOSE ( 90 )
      
    ELSE
       comm_inter = MPI_COMM_WORLD

       IF ( myid < bc_data(1) ) THEN
          inter_color     = 0
          numprocs        = bc_data(1)
          coupling_mode   = 'atmosphere_to_ocean'
       ELSE
          inter_color     = 1
          numprocs        = bc_data(2)
          coupling_mode   = 'ocean_to_atmosphere'
       ENDIF

       CALL MPI_COMM_SPLIT( MPI_COMM_WORLD, inter_color, 0, comm_palm, ierr )
       comm2d = comm_palm

!
!--    Write a flag file for the ocean model and the other atmosphere
!--    processes.
       OPEN( 90, FILE='COUPLING_PORT_OPENED', FORM='FORMATTED' )
       WRITE ( 90, '(''TRUE'')' )
       CLOSE ( 90 )
    ENDIF
#endif

!
!-- In case of a precursor ocean run (followed by a coupled run), or a
!-- coupled atmosphere-ocean run, set the file extension for the ocean files
    IF ( TRIM( coupling_mode ) == 'ocean_to_atmosphere' .OR. bc_data(3) == 1 ) &
    THEN
       coupling_char = '_O'
    ENDIF

    IF (  TRIM( coupling_mode ) == 'vnested_fine' )  THEN
!
!-- Set file extension for vertical nesting
       coupling_char = '_NV'
    ENDIF

 END SUBROUTINE init_coupling
