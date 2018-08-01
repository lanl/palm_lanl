!> @file lpm_calc_liquid_water_content.f90
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
! $Id: lpm_calc_liquid_water_content.f90 3039 2018-05-24 13:13:11Z schwenkel $
! bugfix for lcm with grid stretching
! 
! 2718 2018-01-02 08:49:38Z maronga
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
! 1822 2016-04-07 07:49:42Z hoffmann
! Unused variables removed.
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
!> Calculate the liquid water content for each grid box.
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_calc_liquid_water_content
 

    USE arrays_3d,                                                             &
        ONLY:  dzw, ql, ql_v, ql_vp

    USE cloud_parameters,                                                      &
        ONLY:  rho_l

    USE constants,                                                             &
        ONLY:  pi

    USE control_parameters,                                                    &
        ONLY:  dz, message_string, rho_surface

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE indices,                                                               &
        ONLY:  nxl, nxr, nyn, nys, nzb, nzt

    USE kinds

    USE particle_attributes,                                                   &
        ONLY:  grid_particles, number_of_particles, particles, prt_count

    IMPLICIT NONE

    INTEGER(iwp) ::  i   !<
    INTEGER(iwp) ::  j   !<
    INTEGER(iwp) ::  k   !<
    INTEGER(iwp) ::  n   !<

    CALL cpu_log( log_point_s(45), 'lpm_calc_ql', 'start' )

!
!-- Set water content initially to zero
    ql = 0.0_wp;  ql_v = 0.0_wp;  ql_vp = 0.0_wp

!
!-- Calculate for each grid box
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt

             number_of_particles = prt_count(k,j,i)
             IF ( number_of_particles <= 0 )  CYCLE
             particles => grid_particles(k,j,i)%particles(1:number_of_particles)

!
!--          Calculate the total volume in the boxes (ql_v, weighting factor
!--          has to beincluded)
             DO  n = 1, prt_count(k,j,i)
                ql_v(k,j,i)  = ql_v(k,j,i)  + particles(n)%weight_factor *  &
                                              particles(n)%radius**3
             ENDDO

!
!--          Calculate the liquid water content
             IF ( ql_v(k,j,i) /= 0.0_wp )  THEN
                ql(k,j,i) = ql(k,j,i) + rho_l * 1.33333333_wp * pi *           &
                                        ql_v(k,j,i) /                          &
                                        ( rho_surface * dx * dy * dzw(k) )

                IF ( ql(k,j,i) < 0.0_wp )  THEN
                   WRITE( message_string, * )  'LWC out of range: ' , &
                                               ql(k,j,i),i,j,k
                   CALL message( 'lpm_calc_liquid_water_content', '', 2, 2, &
                                 -1, 6, 1 )
                ENDIF

             ELSE

                ql(k,j,i) = 0.0_wp

             ENDIF

          ENDDO
       ENDDO
    ENDDO

    CALL cpu_log( log_point_s(45), 'lpm_calc_ql', 'stop' )


 END SUBROUTINE lpm_calc_liquid_water_content
