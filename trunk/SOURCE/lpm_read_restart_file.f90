!> @file lpm_read_restart_file.f90
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
! $Id: lpm_read_restart_file.f90 3049 2018-05-29 13:52:36Z Giersch $
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised, code adjusted according to PALMs coding standard
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2606 2017-11-10 10:36:31Z schwenkel
! Changed particle box locations: center of particle box now coincides 
! with scalar grid point of same index.
! Renamed module and subroutines: lpm_pack_arrays_mod -> lpm_pack_and_sort_mod
! lpm_pack_all_arrays -> lpm_sort_in_subboxes, lpm_pack_arrays -> lpm_pack
! lpm_sort -> lpm_sort_timeloop_done
!
! 2312 2017-07-14 20:26:51Z hoffmann
! Extended particle data type.
!
! 2305 2017-07-06 11:18:47Z hoffmann
! Improved calculation of particle IDs.
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
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 849 2012-03-15 10:35:09Z raasch
! initial revision (former part of init_particles)
!
!
! Description:
! ------------
!> Read particle data from the restart file.
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_read_restart_file


    USE control_parameters,                                                    &
        ONLY:  message_string

    USE indices,                                                               &
        ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt

    USE kinds

    USE lpm_pack_and_sort_mod,                                                   &
        ONLY:  lpm_sort_in_subboxes

    USE particle_attributes,                                                   &
        ONLY:  alloc_factor, bc_par_b, bc_par_lr, bc_par_ns, bc_par_t,         &
               grid_particles, min_nr_particle, number_of_particles,           &
               number_of_particle_groups, particle_groups, particle_type,      &
               prt_count, time_prel, time_write_particle_data, zero_particle

    USE pegrid

    IMPLICIT NONE

    CHARACTER (LEN=10) ::  particle_binary_version    !<
    CHARACTER (LEN=10) ::  version_on_file            !<

    INTEGER(iwp) :: alloc_size !<
    INTEGER(iwp) :: ip         !<
    INTEGER(iwp) :: jp         !<
    INTEGER(iwp) :: kp         !<

    TYPE(particle_type), DIMENSION(:), ALLOCATABLE :: tmp_particles !<

!
!-- Read particle data from previous model run.
!-- First open the input unit.
    IF ( myid_char == '' )  THEN
       OPEN ( 90, FILE='PARTICLE_RESTART_DATA_IN'//myid_char,                  &
                  FORM='UNFORMATTED' )
    ELSE
       OPEN ( 90, FILE='PARTICLE_RESTART_DATA_IN/'//myid_char,                 &
                  FORM='UNFORMATTED' )
    ENDIF

!
!-- First compare the version numbers
    READ ( 90 )  version_on_file
    particle_binary_version = '4.0'
    IF ( TRIM( version_on_file ) /= TRIM( particle_binary_version ) )  THEN
       message_string = 'version mismatch concerning data from prior ' //      &
                        'run &version on file = "' //                          &
                                      TRIM( version_on_file ) //               &
                        '&version in program = "' //                           &
                                      TRIM( particle_binary_version ) // '"'
       CALL message( 'lpm_read_restart_file', 'PA0214', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- If less particles are stored on the restart file than prescribed by
!-- min_nr_particle, the remainder is initialized by zero_particle to avoid
!-- errors.
    zero_particle = particle_type( 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,     &
                                   0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,     &
                                   0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,     &
                                   0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,     &
                                   0, 0, 0_idp, .FALSE., -1 )
!
!-- Read some particle parameters and the size of the particle arrays,
!-- allocate them and read their contents.
    READ ( 90 )  bc_par_b, bc_par_lr, bc_par_ns, bc_par_t,                     &
                 number_of_particle_groups, particle_groups, time_prel,        &
                 time_write_particle_data

    ALLOCATE( prt_count(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                        &
              grid_particles(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

    READ ( 90 )  prt_count

    DO  ip = nxl, nxr
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt

             number_of_particles = prt_count(kp,jp,ip)
             IF ( number_of_particles > 0 )  THEN
                alloc_size = MAX( INT( number_of_particles *                   &
                             ( 1.0_wp + alloc_factor / 100.0_wp ) ),           &
                             min_nr_particle )
             ELSE
                alloc_size = min_nr_particle
             ENDIF

             ALLOCATE( grid_particles(kp,jp,ip)%particles(1:alloc_size) )

             IF ( number_of_particles > 0 )  THEN
                ALLOCATE( tmp_particles(1:number_of_particles) )
                READ ( 90 )  tmp_particles
                grid_particles(kp,jp,ip)%particles(1:number_of_particles) = tmp_particles
                DEALLOCATE( tmp_particles )
                IF ( number_of_particles < alloc_size )  THEN
                   grid_particles(kp,jp,ip)%particles(number_of_particles+1:alloc_size) &
                      = zero_particle
                ENDIF
             ELSE
                grid_particles(kp,jp,ip)%particles(1:alloc_size) = zero_particle
             ENDIF

          ENDDO
       ENDDO
    ENDDO

    CLOSE ( 90 )
!
!-- Must be called to sort particles into blocks, which is needed for a fast
!-- interpolation of the LES fields on the particle position.
    CALL lpm_sort_in_subboxes


 END SUBROUTINE lpm_read_restart_file
