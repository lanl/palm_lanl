!> @file compute_vpt.f90
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
! $Id: compute_vpt.f90 2718 2018-01-02 08:49:38Z maronga $
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
! 803 2012-01-16 15:48:46Z franke
! Bugfix: wrong factor in calculation of vpt in case of cloud droplets
!
! Revision 1.1  2000/04/13 14:40:53  schroeter
! Initial revision
!
!
! Description:
! -------------
!> Computation of the virtual potential temperature 
!------------------------------------------------------------------------------!
 SUBROUTINE compute_vpt
 

    USE arrays_3d,                                                             &
        ONLY:  pt, q, ql, vpt
        
    USE indices,                                                               &
        ONLY:  nzb, nzt
        
    USE cloud_parameters,                                                      &
        ONLY:  l_d_cp, pt_d_t
        
    USE control_parameters,                                                    &
        ONLY:  cloud_droplets, cloud_physics
        
    USE kinds

    IMPLICIT NONE

    INTEGER(iwp) :: k   !< 

    IF ( .NOT. cloud_physics  .AND.  .NOT. cloud_droplets )  THEN
       vpt = pt * ( 1.0_wp + 0.61_wp * q )
    ELSE IF (cloud_physics)  THEN
       DO  k = nzb, nzt+1
          vpt(k,:,:) = ( pt(k,:,:) + pt_d_t(k) * l_d_cp * ql(k,:,:) ) *        &
                       ( 1.0_wp + 0.61_wp * q(k,:,:) - 1.61_wp * ql(k,:,:) ) 
       ENDDO
    ELSE
       vpt = pt * ( 1.0_wp + 0.61_wp * q - ql ) 
    ENDIF

 END SUBROUTINE compute_vpt
