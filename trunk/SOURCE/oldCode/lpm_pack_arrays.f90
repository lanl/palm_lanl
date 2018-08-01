!> @file lpm_pack_arrays.f90
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
! $Id: lpm_pack_arrays.f90 2801 2018-02-14 16:01:55Z thiele $
! Introduce particle transfer in nested models.
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2628 2017-11-20 12:40:38Z schwenkel
! Enabled particle advection with grid stretching.
! 
! 2609 2017-11-14 14:14:44Z schwenkel
! Integrated subroutine pack_and_sort into lpm_sort_in_subboxes
! 
! 2606 2017-11-10 10:36:31Z schwenkel
! Changed particle box locations: center of particle box now coincides 
! with scalar grid point of same index.
! Renamed module and subroutines: lpm_pack_arrays_mod -> lpm_pack_and_sort_mod
! lpm_pack_all_arrays -> lpm_sort_in_subboxes, lpm_pack_arrays -> lpm_pack
! lpm_sort -> lpm_sort_timeloop_done
! 
! 2417 2017-09-06 15:22:27Z suehring
! New routine which sorts particles into particles that completed and not 
! completed the LES timestep. 
! 
! 2101 2017-01-05 16:42:31Z suehring
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1822 2016-04-07 07:49:42Z hoffmann
! Tails removed. Unused variables removed.
!
! 1685 2015-10-08 07:32:13Z raasch
! bugfix concerning vertical index calculation in case of ocean
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1359 2014-04-11 17:15:14Z hoffmann
! New particle structure integrated. 
! Kind definition added to all floating point numbers.
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
!> Pack particle arrays, which means eliminate those elements marked for
!> deletion and move data with higher index values to these free indices.
!> Determine the new number of particles.
!> Moreover, particles are also sorted into groups finished and not finished
!> its timestep. 
!------------------------------------------------------------------------------!
 MODULE lpm_pack_and_sort_mod
 

    USE particle_attributes,                                                   &
        ONLY:  grid_particles, number_of_particles, offset_ocean_nzt,          &
               particles, particle_type, prt_count

    PRIVATE
    PUBLIC lpm_sort_in_subboxes, lpm_pack, lpm_sort_timeloop_done

    INTERFACE lpm_sort_in_subboxes
       MODULE PROCEDURE lpm_sort_in_subboxes
    END INTERFACE lpm_sort_in_subboxes

    INTERFACE lpm_pack
       MODULE PROCEDURE lpm_pack
    END INTERFACE lpm_pack

    INTERFACE lpm_sort_timeloop_done
       MODULE PROCEDURE lpm_sort_timeloop_done
    END INTERFACE lpm_sort_timeloop_done


 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! -----------
!> Routine for the whole processor
!> Sort all particles into the 8 respective subgrid boxes
!------------------------------------------------------------------------------!
    SUBROUTINE lpm_sort_in_subboxes

       USE cpulog,                                                              &
          ONLY:  cpu_log, log_point_s

       USE indices,                                                             &
          ONLY:  nxl, nxr, nys, nyn, nzb, nzt

       USE kinds

       USE control_parameters,                                                  &
          ONLY: dz
       
       USE grid_variables,                                                      &
          ONLY: dx,dy,ddx, ddy
           
       USE arrays_3d,                                                           &
          ONLY:  zu
       IMPLICIT NONE

       INTEGER(iwp) ::  i  !<
       INTEGER(iwp) ::  ip !<
       INTEGER(iwp) ::  is !<
       INTEGER(iwp) ::  j  !<
       INTEGER(iwp) ::  jp !<
       INTEGER(iwp) ::  k  !<
       INTEGER(iwp) ::  kp !<
       INTEGER(iwp) ::  m  !<
       INTEGER(iwp) ::  n  !<
       INTEGER(iwp) ::  nn !<
       INTEGER(iwp) ::  sort_index  !<

       INTEGER(iwp), DIMENSION(0:7) ::  sort_count  !<

       TYPE(particle_type), DIMENSION(:,:), ALLOCATABLE ::  sort_particles  !<

       CALL cpu_log( log_point_s(51), 'lpm_sort_in_subboxes', 'start' )
       DO  ip = nxl, nxr
          DO  jp = nys, nyn
             DO  kp = nzb+1, nzt
                number_of_particles = prt_count(kp,jp,ip)
                IF ( number_of_particles <= 0 )  CYCLE
                particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
                      
                nn = 0
                sort_count = 0
                ALLOCATE( sort_particles(number_of_particles, 0:7) )
                
                DO  n = 1, number_of_particles
                   sort_index = 0

                   IF ( particles(n)%particle_mask )  THEN
                      nn = nn + 1
!
!--                   Sorting particles with a binary scheme
!--                   sort_index=111_2=7_10 -> particle at the left,south,bottom subgridbox
!--                   sort_index=000_2=0_10 -> particle at the right,north,top subgridbox
!--                   For this the center of the gridbox is calculated 
                      i = (particles(n)%x + 0.5_wp * dx) * ddx
                      j = (particles(n)%y + 0.5_wp * dy) * ddy
                      
                      IF ( i == ip )  sort_index = sort_index + 4 
                      IF ( j == jp )  sort_index = sort_index + 2                      
                      IF ( zu(kp) > particles(n)%z ) sort_index = sort_index + 1
                      
                      sort_count(sort_index) = sort_count(sort_index) + 1
                      m = sort_count(sort_index)
                      sort_particles(m,sort_index) = particles(n)
                      sort_particles(m,sort_index)%block_nr = sort_index
                   ENDIF
                ENDDO

                nn = 0
                DO is = 0,7
                   grid_particles(kp,jp,ip)%start_index(is) = nn + 1
                   DO n = 1,sort_count(is)
                      nn = nn + 1
                      particles(nn) = sort_particles(n,is)
                   ENDDO
                   grid_particles(kp,jp,ip)%end_index(is) = nn
                ENDDO

                number_of_particles = nn    
                prt_count(kp,jp,ip) = number_of_particles
                DEALLOCATE(sort_particles)
             ENDDO
          ENDDO
       ENDDO
       CALL cpu_log( log_point_s(51), 'lpm_sort_in_subboxes', 'stop' )

    END SUBROUTINE lpm_sort_in_subboxes

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Move all particles not marked for deletion to lowest indices (packing)
!------------------------------------------------------------------------------!
    SUBROUTINE lpm_pack

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  n       !<
       INTEGER(iwp) ::  nn      !<
!
!--    Find out elements marked for deletion and move data from highest index
!--    values to these free indices
       nn = number_of_particles

       DO WHILE ( .NOT. particles(nn)%particle_mask )
          nn = nn-1
          IF ( nn == 0 )  EXIT
       ENDDO

       IF ( nn > 0 )  THEN
          DO  n = 1, number_of_particles
             IF ( .NOT. particles(n)%particle_mask )  THEN
                particles(n) = particles(nn)
                nn = nn - 1
                DO WHILE ( .NOT. particles(nn)%particle_mask )
                   nn = nn-1
                   IF ( n == nn )  EXIT
                ENDDO
             ENDIF
             IF ( n == nn )  EXIT
          ENDDO
       ENDIF

!
!--    The number of deleted particles has been determined in routines
!--    lpm_boundary_conds, lpm_droplet_collision, and lpm_exchange_horiz
       number_of_particles = nn

    END SUBROUTINE lpm_pack 
                
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sort particles in each sub-grid box into two groups: particles that already
!> completed the LES timestep, and particles that need further timestepping to
!> complete the LES timestep. 
!------------------------------------------------------------------------------!
    SUBROUTINE lpm_sort_timeloop_done

       USE control_parameters,                                                 &
           ONLY:  dt_3d

       USE indices,                                                            &
           ONLY: nxl, nxr, nys, nyn, nzb, nzt

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) :: end_index     !< particle end index for each sub-box
       INTEGER(iwp) :: i             !< index of particle grid box in x-direction
       INTEGER(iwp) :: j             !< index of particle grid box in y-direction
       INTEGER(iwp) :: k             !< index of particle grid box in z-direction
       INTEGER(iwp) :: n             !< running index for number of particles
       INTEGER(iwp) :: nb            !< index of subgrid boux
       INTEGER(iwp) :: nf            !< indices for particles in each sub-box that already finalized their substeps
       INTEGER(iwp) :: nnf           !< indices for particles in each sub-box that need further treatment
       INTEGER(iwp) :: num_finalized !< number of particles in each sub-box that already finalized their substeps
       INTEGER(iwp) :: start_index   !< particle start index for each sub-box

       TYPE(particle_type), DIMENSION(:), ALLOCATABLE :: sort_particles  !< temporary particle array 

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt

                number_of_particles = prt_count(k,j,i)
                IF ( number_of_particles <= 0 )  CYCLE

                particles => grid_particles(k,j,i)%particles(1:number_of_particles)

                DO  nb = 0, 7
!
!--                Obtain start and end index for each subgrid box
                   start_index = grid_particles(k,j,i)%start_index(nb)
                   end_index   = grid_particles(k,j,i)%end_index(nb)
!
!--                Allocate temporary array used for sorting. 
                   ALLOCATE( sort_particles(start_index:end_index) )
!
!--                Determine number of particles already completed the LES 
!--                timestep, and write them into a temporary array. 
                   nf = start_index
                   num_finalized = 0
                   DO  n = start_index, end_index
                      IF ( dt_3d - particles(n)%dt_sum < 1E-8_wp )  THEN
                         sort_particles(nf) = particles(n)
                         nf                 = nf + 1
                         num_finalized      = num_finalized + 1
                      ENDIF
                   ENDDO
!
!--                Determine number of particles that not completed the LES 
!--                timestep, and write them into a temporary array. 
                   nnf = nf
                   DO  n = start_index, end_index
                      IF ( dt_3d - particles(n)%dt_sum > 1E-8_wp )  THEN
                         sort_particles(nnf) = particles(n)
                         nnf                 = nnf + 1
                      ENDIF
                   ENDDO
!
!--                Write back sorted particles
                   particles(start_index:end_index) =                          &
                                           sort_particles(start_index:end_index)
!
!--                Determine updated start_index, used to masked already 
!--                completed particles. 
                   grid_particles(k,j,i)%start_index(nb) =                     &
                                      grid_particles(k,j,i)%start_index(nb)    &
                                    + num_finalized
!
!--                Deallocate dummy array
                   DEALLOCATE ( sort_particles )
!
!--                Finally, if number of non-completed particles is non zero 
!--                in any of the sub-boxes, set control flag appropriately. 
                   IF ( nnf > nf )                                             &
                      grid_particles(k,j,i)%time_loop_done = .FALSE.

                ENDDO
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE lpm_sort_timeloop_done


 END MODULE lpm_pack_and_sort_mod
