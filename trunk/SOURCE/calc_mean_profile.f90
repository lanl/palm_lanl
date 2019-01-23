!> @file calc_mean_profile.f90
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
! $Id: calc_mean_profile.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Enable usage also during initialization, required for initialization of 
! radiation and land-surface model (MS)
! 
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography concept
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1873 2016-04-18 14:50:06Z maronga
! Module renamed (removed _mod)
! 
! 
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
! 
!
! 1738 2015-12-18 13:56:05Z raasch
! bugfix: if a layer is completely filled with topography, no mean is calculated
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1365 2014-04-22 15:03:56Z boeske
! Initial revision
!
! Description:
! ------------
!> Calculate the horizontally averaged vertical temperature profile (pr=4 in case
!> of potential temperature, 44 in case of virtual potential temperature, and 64
!> in case of density (ocean runs)).
!------------------------------------------------------------------------------!
 MODULE calc_mean_profile_mod
 

    PRIVATE
    PUBLIC calc_mean_profile

    INTERFACE calc_mean_profile
       MODULE PROCEDURE calc_mean_profile
    END INTERFACE calc_mean_profile

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_mean_profile( var, pr )

       USE control_parameters,                                                 &
           ONLY:  intermediate_timestep_count, message_string

       USE indices,                                                            &
           ONLY:  ngp_2dh_s_inner, nxl, nxr, nyn, nys, nzb, nzb, nzt,          &
                  wall_flags_0

       USE kinds

       USE pegrid

       USE statistics,                                                         &
           ONLY:  flow_statistics_called, hom, sums, sums_l


       IMPLICIT NONE
       
       INTEGER(iwp) ::  i                  !<
       INTEGER(iwp) ::  j                  !<
       INTEGER(iwp) ::  k                  !<
       INTEGER(iwp) ::  pr                 !< 
       INTEGER(iwp) ::  omp_get_thread_num !<
       INTEGER(iwp) ::  tn                 !<
       
#if defined( __nopointer )
       REAL(wp), DIMENSION(:,:,:) ::  var  !<
#else
       REAL(wp), DIMENSION(:,:,:), POINTER ::  var
#endif

!
!--    Computation of the horizontally averaged profile of variable var, unless
!--    already done by the relevant call from flow_statistics. The calculation
!--    is done only for the first respective intermediate timestep in order to
!--    spare communication time and to produce identical model results with jobs
!--    which are calling flow_statistics at different time intervals. At 
!--    initialization, intermediate_timestep_count = 0 is considered as well.

       IF ( .NOT. flow_statistics_called  .AND.                                &
            intermediate_timestep_count <= 1 )  THEN

!
!--       Horizontal average of variable var
          tn           =   0  ! Default thread number in case of one thread
          !$OMP PARALLEL PRIVATE( i, j, k, tn )
!$        tn = omp_get_thread_num()
          sums_l(:,pr,tn) = 0.0_wp
          !$OMP DO
          DO  i = nxl, nxr
             DO  j =  nys, nyn
                DO  k = nzb, nzt+1
                   sums_l(k,pr,tn) = sums_l(k,pr,tn) + var(k,j,i)              &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,i), 22 ) )
                ENDDO
             ENDDO
          ENDDO
          !$OMP END PARALLEL

          DO  i = 1, threads_per_task-1
             sums_l(:,pr,0) = sums_l(:,pr,0) + sums_l(:,pr,i)
          ENDDO

#if defined( __parallel )

          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_ALLREDUCE( sums_l(nzb,pr,0), sums(nzb,pr), nzt+2-nzb,       &
                              MPI_REAL, MPI_SUM, comm2d, ierr )

#else

          sums(:,pr) = sums_l(:,pr,0)

#endif


          DO  k = nzb, nzt+1
             IF ( ngp_2dh_s_inner(k,0) /= 0 )  THEN
                hom(k,1,pr,0) = sums(k,pr) / ngp_2dh_s_inner(k,0)
             ENDIF
          ENDDO

       ENDIF


    END SUBROUTINE calc_mean_profile

 END MODULE calc_mean_profile_mod
