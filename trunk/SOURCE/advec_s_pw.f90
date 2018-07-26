!> @file advec_s_pw.f90
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
! $Id: advec_s_pw.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! topography representation via flags
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
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1374 2014-04-25 12:55:07Z raasch
! missing variables added to ONLY list
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
! 1010 2012-09-20 07:59:54Z raasch
! cpp switch __nopointer added for pointer free version
!
! Revision 1.1  1997/08/29 08:54:20  raasch
! Initial revision
!
!
! Description:
! ------------
!> Advection term for scalar variables using the Piacsek and Williams scheme
!> (form C3). Contrary to PW itself, for reasons of accuracy their scheme is
!> slightly modified as follows: the values of those scalars that are used for
!> the computation of the flux divergence are reduced by the value of the
!> relevant scalar at the location where the difference is computed (sk(k,j,i)).
!> NOTE: at the first grid point above the surface computation still takes place!
!------------------------------------------------------------------------------!
 MODULE advec_s_pw_mod
 

    PRIVATE
    PUBLIC advec_s_pw

    INTERFACE advec_s_pw
       MODULE PROCEDURE advec_s_pw
       MODULE PROCEDURE advec_s_pw_ij
    END INTERFACE
 
 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE advec_s_pw( sk )

       USE arrays_3d,                                                          &
           ONLY:  dd2zu, tend, u, v, w

       USE control_parameters,                                                 &
           ONLY:  u_gtrans, v_gtrans

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb,             &
                  nzt, wall_flags_0

       USE kinds


       IMPLICIT NONE

       INTEGER(iwp) ::  i !<
       INTEGER(iwp) ::  j !<
       INTEGER(iwp) ::  k !<

#if defined( __nopointer )
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  sk !<
#else
       REAL(wp), DIMENSION(:,:,:), POINTER ::  sk
#endif
 

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                tend(k,j,i) = tend(k,j,i) +                                       &
            ( -0.5_wp * ( ( u(k,j,i+1) - u_gtrans ) * ( sk(k,j,i+1) - sk(k,j,i) ) &
                        - ( u(k,j,i)   - u_gtrans ) * ( sk(k,j,i-1) - sk(k,j,i) ) &
                        ) * ddx                                                   &
              -0.5_wp * ( ( v(k,j+1,i) - v_gtrans ) * ( sk(k,j+1,i) - sk(k,j,i) ) &
                        - ( v(k,j,i)   - v_gtrans ) * ( sk(k,j-1,i) - sk(k,j,i) ) &
                        ) * ddy                                                   &
              -         (   w(k,j,i)                * ( sk(k+1,j,i) - sk(k,j,i) ) &
                        -   w(k-1,j,i)              * ( sk(k-1,j,i) - sk(k,j,i) ) &
                        ) * dd2zu(k)                                              &
            ) * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE advec_s_pw


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE advec_s_pw_ij( i, j, sk )

       USE arrays_3d,                                                          &
           ONLY:  dd2zu, tend, u, v, w

       USE control_parameters,                                                 &
           ONLY:  u_gtrans, v_gtrans

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nyng, nysg, nzb, nzt, wall_flags_0

       USE kinds


       IMPLICIT NONE

       INTEGER(iwp) ::  i !<
       INTEGER(iwp) ::  j !<
       INTEGER(iwp) ::  k !<

#if defined( __nopointer )
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  sk !<
#else
       REAL(wp), DIMENSION(:,:,:), POINTER ::  sk
#endif


       DO  k = nzb+1, nzt
          tend(k,j,i) = tend(k,j,i) +                                             &
            ( -0.5_wp * ( ( u(k,j,i+1) - u_gtrans ) * ( sk(k,j,i+1) - sk(k,j,i) ) &
                        - ( u(k,j,i)   - u_gtrans ) * ( sk(k,j,i-1) - sk(k,j,i) ) &
                        ) * ddx                                                   &
              -0.5_wp * ( ( v(k,j+1,i) - v_gtrans ) * ( sk(k,j+1,i) - sk(k,j,i) ) &
                        - ( v(k,j,i)   - v_gtrans ) * ( sk(k,j-1,i) - sk(k,j,i) ) &
                        ) * ddy                                                   &
              -         (   w(k,j,i)                * ( sk(k+1,j,i) - sk(k,j,i) ) &
                        -   w(k-1,j,i)              * ( sk(k-1,j,i) - sk(k,j,i) ) &
                        ) * dd2zu(k)                                              &
            ) * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )
       ENDDO

    END SUBROUTINE advec_s_pw_ij

 END MODULE advec_s_pw_mod
