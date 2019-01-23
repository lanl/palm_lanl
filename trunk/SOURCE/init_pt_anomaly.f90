!> @file init_pt_anomaly.f90
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
! 2018-11-15 cbegeman
! Change ptanom to 2D and 3D bubble with user-specified properties
! 
! Former revisions:
! -----------------
! $Id: init_pt_anomaly.f90 3035 2018-05-24 09:35:20Z schwenkel $
! Add option to initialize warm air bubble close to surface
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
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL constants defined as wp_kind
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
! 861 2012-03-26 14:18:34Z suehring
! Modification of the amplitude to obtain a visible temperature perturbation.
!
! Revision 1.1  1997/08/29 08:58:56  raasch
! Initial revision
!
!
! Description:
! ------------
!> Impose a temperature perturbation for an advection test.
!------------------------------------------------------------------------------!
 SUBROUTINE init_pt_anomaly
 

    USE arrays_3d,                                                             &
        ONLY:  pt, sa, zu

    USE constants,                                                             &
        ONLY:  pi

    USE control_parameters    
        
    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxr, ny, nyn, nys, nzb, nzt
        
    USE kinds

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< grid index along x
    INTEGER(iwp) ::  j  !< grid index along y
    INTEGER(iwp) ::  k  !< grid index along z
    
    REAL(wp)     ::  bubble_dr                            !< distance from the center of the bubble    
    
!
!-- Set default bubble center to the center of the domain
    IF ( bubble_center_x == 9999999.9_wp ) bubble_center_x = dx * ( nx+1 ) / 2
    IF ( bubble_center_y == 9999999.9_wp ) bubble_center_y = dy * ( ny+1 ) / 2
    IF ( bubble_center_z == 9999999.9_wp ) bubble_center_z = (zu(nzt) - zu(nzb)) / 2
    IF ( bubble_pt == 9999999.9_wp ) bubble_pt = 0.
    IF ( bubble_sa == 9999999.9_wp .AND. ocean ) bubble_sa = 0.

!
!--    Compute the perturbation.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             IF ( INDEX( initializing_actions, 'initialize_2D_bubble' ) /= 0   &
                  .AND. bubble_radius /= 0 )  THEN

                bubble_dr = SQRT( ( dy*j - bubble_center_y )**2 +              &
                                  ( zu(k) - bubble_center_z )**2 )

                IF ( bubble_dr <= bubble_radius )  THEN

                   pt(k,j,i) = pt(k,j,i) + bubble_pt *                            &
                            EXP( -0.5 * ( (dy*j  - bubble_center_y) /          &
                                                   bubble_radius )**2) *       &
                            EXP( -0.5 * ( (zu(k) - bubble_center_z) /          &
                                                   bubble_radius)**2)
                   IF ( ocean ) THEN

                      sa(k,j,i) = sa(k,j,i) + bubble_sa *                         &
                               EXP( -0.5 * ( (dy*j  - bubble_center_y) /       &
                                                      bubble_radius )**2) *    &
                               EXP( -0.5 * ( (zu(k) - bubble_center_z) /       &
                                                      bubble_radius )**2)

                   ENDIF
                ENDIF
             ELSE

                bubble_dr = SQRT( ( dx*i - bubble_center_x )**2 +              &
                                  ( dy*j - bubble_center_y )**2 +              &
                                  ( zu(k) - bubble_center_z )**2 )

                IF ( bubble_dr <= bubble_radius )  THEN

                   IF ( INDEX( initializing_actions, 'initialize_3D_bubble' ) /= 0 &
                            .AND. bubble_radius /= 0 )  THEN

                      pt(k,j,i) = pt(k,j,i) + bubble_pt*cos(pi*bubble_dr/(2*bubble_radius))

                      IF ( ocean ) THEN                                           
                         sa(k,j,i) = sa(k,j,i) + bubble_sa*cos(pi*bubble_dr/(2*bubble_radius))
                      ENDIF

                   ENDIF
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO

!
!-- Exchange of boundary values
    CALL exchange_horiz( pt, nbgp )
    IF ( ocean ) CALL exchange_horiz( sa, nbgp )


 END SUBROUTINE init_pt_anomaly
