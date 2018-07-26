!> @file lpm_merging.f90
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
! $Id: lpm_merging.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 
! Change in file header (GPL part)
!
! Added comments
! 
! 
! 2263 2017-06-08 14:59:01Z schwenkel
! Initial revision
! 
! 
!
! Description:
! ------------
! This routine is a part of the Lagrangian particle model. Two Super droplets 
! which fulfill certain criterion's (e.g. a big weighting factor and a small
! radius) can be merged into one super droplet with a increased number of 
! represented particles of the super droplet. This mechanism ensures an
! improved a feasible amount of computational costs. The limits of particle 
! creation should be chosen carefully! The idea of this algorithm is based on 
! Unterstrasser and Soelch, 2014. 
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_merging


    USE arrays_3d,                                                             &
        ONLY:  ql

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE indices,                                                               &
        ONLY:  nxl, nxr, nyn, nys, nzb, nzt

    USE kinds

    USE particle_attributes,                                                   &
        ONLY:  deleted_particles, grid_particles, initial_weighting_factor,    &     
               merge_drp, merging, number_of_particles, particles, prt_count,  &
               radius_merge, sum_merge_drp, weight_factor_merge

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i         !<      
    INTEGER(iwp) ::  j         !<       
    INTEGER(iwp) ::  k         !<       
    INTEGER(iwp) ::  n         !<   
    
    REAL(wp) ::  ql_crit = 1.0E-5_wp  !< threshold lwc for cloudy grid cells 
                                      !< (e.g. Siebesma et al 2003, JAS, 60)
    
    CALL cpu_log( log_point_s(81), 'lpm_merging', 'start' )

    merge_drp  = 0
    
    IF ( weight_factor_merge == -1.0_wp )  THEN
       weight_factor_merge = 0.5_wp * initial_weighting_factor 
    ENDIF

    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
   
             number_of_particles = prt_count(k,j,i)
             IF ( number_of_particles <= 0  .OR.                               &
                   ql(k,j,i) >= ql_crit )  CYCLE
             particles => grid_particles(k,j,i)%particles(1:number_of_particles)
!
!--          Start merging operations: This routine delete super droplets with
!--          a small radius (radius <= radius_merge) and a low weighting 
!--          factor (weight_factor  <= weight_factor_merge). The number of 
!--          represented particles will be added to the next particle of the 
!--          particle array. Tests showed that this simplified method can be 
!--          used because it will only take place outside of cloudy grid 
!--          boxes where ql <= 1.0E-5 kg/kg. Therefore, especially former cloned   
!--          and subsequent evaporated super droplets will be merged.
             DO  n = 1, number_of_particles-1
                IF ( particles(n)%particle_mask                    .AND.       &
                     particles(n+1)%particle_mask                  .AND.       &
                     particles(n)%radius        <= radius_merge    .AND.       &
                     particles(n)%weight_factor <= weight_factor_merge )       &   
                THEN
                   particles(n+1)%weight_factor  =                             &
                                       particles(n+1)%weight_factor +          &
                                       ( particles(n)%radius**3     /          &
                                         particles(n+1)%radius**3   *          &
                                         particles(n)%weight_factor            &
                                       )
                   particles(n)%particle_mask = .FALSE.
                   deleted_particles          = deleted_particles + 1 
                   merge_drp                  = merge_drp + 1
                
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
       
    sum_merge_drp = sum_merge_drp + merge_drp

    CALL cpu_log( log_point_s(81), 'lpm_merging', 'stop' )

 END SUBROUTINE lpm_merging
