!> @file lpm_write_restart_file.f90
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
! $Id: lpm_write_restart_file.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2265 2017-06-08 16:58:28Z schwenkel
! Unused variables removed.
! 
! 2101 2017-01-05 16:42:31Z suehring
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1822 2016-04-07 07:49:42Z hoffmann
! Tails removed. Unused variables removed.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1359 2014-04-11 17:15:14Z hoffmann
! New particle structure integrated. 
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 849 2012-03-15 10:35:09Z raasch
! initial revision (former part of advec_particles)
!
!
! Description:
! ------------
!> Write particle data in FORTRAN binary format on restart file
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_write_restart_file


    USE indices,                                                               &
        ONLY:  nxl, nxr, nyn, nys, nzb, nzt

    USE kinds

    USE particle_attributes,                                                   &
        ONLY:  bc_par_b, bc_par_lr, bc_par_ns, bc_par_t, grid_particles,       &
               number_of_particles, number_of_particle_groups,                 &
               particles, particle_groups, prt_count, time_prel,               &
               time_write_particle_data

    USE pegrid

    IMPLICIT NONE

    CHARACTER (LEN=10) ::  particle_binary_version   !<

    INTEGER(iwp) ::  ip                              !<
    INTEGER(iwp) ::  jp                              !<
    INTEGER(iwp) ::  kp                              !<

!
!-- First open the output unit.
    IF ( myid_char == '' )  THEN
       OPEN ( 90, FILE='PARTICLE_RESTART_DATA_OUT'//myid_char, &
                  FORM='UNFORMATTED')
    ELSE
       IF ( myid == 0 )  CALL local_system( 'mkdir PARTICLE_RESTART_DATA_OUT' )
#if defined( __parallel )
!
!--    Set a barrier in order to allow that thereafter all other processors
!--    in the directory created by PE0 can open their file
       CALL MPI_BARRIER( comm2d, ierr )
#endif
       OPEN ( 90, FILE='PARTICLE_RESTART_DATA_OUT/'//myid_char, &
                  FORM='UNFORMATTED' )
    ENDIF

!
!-- Write the version number of the binary format.
!-- Attention: After changes to the following output commands the version
!-- ---------  number of the variable particle_binary_version must be
!--            changed! Also, the version number and the list of arrays
!--            to be read in lpm_read_restart_file must be adjusted
!--            accordingly.
    particle_binary_version = '4.0'
    WRITE ( 90 )  particle_binary_version

!
!-- Write some particle parameters, the size of the particle arrays as
!-- well as other dvrp-plot variables.
    WRITE ( 90 )  bc_par_b, bc_par_lr, bc_par_ns, bc_par_t,              &
                  number_of_particle_groups, particle_groups, time_prel, &
                  time_write_particle_data

    WRITE ( 90 )  prt_count
          
    DO  ip = nxl, nxr
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt
             number_of_particles = prt_count(kp,jp,ip)
             particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
             IF ( number_of_particles <= 0 )  CYCLE
             WRITE ( 90 )  particles
          ENDDO
       ENDDO
    ENDDO

    CLOSE ( 90 )

#if defined( __parallel )
       CALL MPI_BARRIER( comm2d, ierr )
#endif


 END SUBROUTINE lpm_write_restart_file
