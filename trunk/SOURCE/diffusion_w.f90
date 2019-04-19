!> @file diffusion_w.f90
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
! $Id: diffusion_w.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography and surface concept
!
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC version of subroutine removed
!
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1873 2016-04-18 14:50:06Z maronga
! Module renamed (removed _mod)
!
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1374 2014-04-25 12:55:07Z raasch
! vsws + vswst removed from acc-present-list
!
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1340 2014-03-25 19:45:13Z kanani
! REAL constants defined as wp-kind
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL constants defined as wp-kind
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1257 2013-11-08 15:18:40Z raasch
! openacc loop and loop vector clauses removed, declare create moved after
! the FORTRAN declaration statement
!
! 1128 2013-04-12 06:19:32Z raasch
! loop index bounds in accelerator version replaced by i_left, i_right, j_south,
! j_north
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1015 2012-09-27 09:23:24Z raasch
! accelerator version (*_acc) added
!
! 1001 2012-09-13 14:08:46Z raasch
! arrays comunicated by module instead of parameter list
!
! 978 2012-08-09 08:28:32Z fricke
! outflow damping layer removed
! kmxm_x/_z and kmxp_x/_z change to kmxm and kmxp
! kmym_y/_z and kmyp_y/_z change to kmym and kmyp
!
! Revision 1.1  1997/09/12 06:24:11  raasch
! Initial revision
!
!
! Description:
! ------------
!> Diffusion term of the w-component
!------------------------------------------------------------------------------!
 MODULE diffusion_w_mod


    PRIVATE
    PUBLIC diffusion_w

    INTERFACE diffusion_w
       MODULE PROCEDURE diffusion_w
    END INTERFACE diffusion_w

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE diffusion_w

       USE arrays_3d,                                                          &
           ONLY :  ddzu, ddzw, km, tend, u, v, w, drho_air_zw, rho_air

       USE grid_variables,                                                     &
           ONLY :  ddx, ddy

       USE indices,                                                            &
           ONLY :  nxl, nxr, nyn, nys, nzb, nzt

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  l             !< running index of surface type, south- or north-facing wall
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp) ::  kmxm              !<
       REAL(wp) ::  kmxp              !<
       REAL(wp) ::  kmym              !<
       REAL(wp) ::  kmyp              !<



       !$acc data copy( tend ) &
       !$acc copyin( u, v, w ) &
       !$acc present( km ) &
       !$acc present( ddzu, ddzw, rho_air, drho_air_zw )

       !$acc parallel
       !$acc loop gang vector collapse(3)
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt-1
!
!--             Predetermine flag to mask topography and wall-bounded grid points.
!
!--             Interpolate eddy diffusivities on staggered gridpoints
                kmxp = 0.25_wp * ( km(k,j,i)   +   km(k,j,i+1) +               &
                                   km(k+1,j,i) + km(k+1,j,i+1) )
                kmxm = 0.25_wp * ( km(k,j,i)   + km(k,j,i-1)   +               &
                                   km(k+1,j,i) + km(k+1,j,i-1) )
                kmyp = 0.25_wp * ( km(k,j,i)   + km(k+1,j,i)   +               &
                                   km(k,j+1,i) + km(k+1,j+1,i) )
                kmym = 0.25_wp * ( km(k,j,i)   + km(k+1,j,i)   +               &
                                   km(k,j-1,i) + km(k+1,j-1,i) )

                tend(k,j,i) = tend(k,j,i)                                      &
                                    + ( kmxp * (                               &
                                   ( w(k,j,i+1)   - w(k,j,i)   ) * ddx         &
                                 + ( u(k+1,j,i+1) - u(k,j,i+1) ) * ddzu(k+1)   &
                                               )                               &
                                     - kmxm *  (                               &
                                   ( w(k,j,i)     - w(k,j,i-1) ) * ddx         &
                                 + ( u(k+1,j,i)   - u(k,j,i)   ) * ddzu(k+1)   &
                                               )                               &
                                      ) * ddx                                  &
                                    + ( kmyp * (                               &
                                   ( w(k,j+1,i)   - w(k,j,i)   ) * ddy         &
                                 + ( v(k+1,j+1,i) - v(k,j+1,i) ) * ddzu(k+1)   &
                                               )                               &
                                      - kmym * (                               &
                                   ( w(k,j,i)     - w(k,j-1,i) ) * ddy         &
                                 + ( v(k+1,j,i)   - v(k,j,i)   ) * ddzu(k+1)   &
                                               )                               &
                                      ) * ddy                                  &
                       + 2.0_wp * (                                            &
                         km(k+1,j,i) * ( w(k+1,j,i) - w(k,j,i) )   * ddzw(k+1) &
                                     * rho_air(k+1)                            &
                       - km(k,j,i)   * ( w(k,j,i)   - w(k-1,j,i) ) * ddzw(k)   &
                                     * rho_air(k)                              &
                                  ) * ddzu(k+1) * drho_air_zw(k)
             ENDDO
          ENDDO
       ENDDO
       !$acc end parallel
       !$acc end data

    END SUBROUTINE diffusion_w


 END MODULE diffusion_w_mod
