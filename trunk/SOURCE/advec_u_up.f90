!> @file advec_u_up.f90
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
! $Id: advec_u_up.f90 2718 2018-01-02 08:49:38Z maronga $
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
! Revision 1.1  1997/08/29 08:55:25  raasch
! Initial revision
!
!
! Description:
! ------------
!> Advection term for the u velocity-component using upstream scheme.
!> NOTE: vertical advection at k=1 still has wrong grid spacing for w>0!
!>       The same problem occurs for all topography boundaries!
!------------------------------------------------------------------------------!
 MODULE advec_u_up_mod
 

    PRIVATE
    PUBLIC advec_u_up

    INTERFACE advec_u_up
       MODULE PROCEDURE advec_u_up
       MODULE PROCEDURE advec_u_up_ij
    END INTERFACE advec_u_up

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE advec_u_up

       USE arrays_3d,                                                          &
           ONLY:  ddzu, tend, u, v, w

       USE control_parameters,                                                 &
           ONLY:  u_gtrans, v_gtrans

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxlu, nxr, nyn, nys, nzb, nzt, wall_flags_0

       USE kinds


       IMPLICIT NONE

       INTEGER(iwp) ::  i !<
       INTEGER(iwp) ::  j !<
       INTEGER(iwp) ::  k !<

       REAL(wp) ::  ukomp !<
       REAL(wp) ::  vkomp !<
       REAL(wp) ::  wkomp !<

       
       DO  i = nxlu, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             x-direction
                ukomp = u(k,j,i) - u_gtrans
                IF ( ukomp > 0.0_wp )  THEN
                   tend(k,j,i) = tend(k,j,i) - ukomp *                         &
                                         ( u(k,j,i) - u(k,j,i-1) ) * ddx       &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 1 ) )
                ELSE
                   tend(k,j,i) = tend(k,j,i) - ukomp *                         &
                                          ( u(k,j,i+1) - u(k,j,i) ) * ddx      &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 1 ) )
                ENDIF
!
!--             y-direction
                vkomp = 0.25_wp * ( v(k,j,i)   + v(k,j+1,i) +                  &
                                 v(k,j,i-1) + v(k,j+1,i-1) ) - v_gtrans
                IF ( vkomp > 0.0_wp )  THEN
                   tend(k,j,i) = tend(k,j,i) - vkomp *                         &
                                         ( u(k,j,i) - u(k,j-1,i) ) * ddy       &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 1 ) )
                ELSE
                   tend(k,j,i) = tend(k,j,i) - vkomp *                         &
                                         ( u(k,j+1,i) - u(k,j,i) ) * ddy       &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 1 ) )
                ENDIF
!
!--             z-direction
                wkomp = 0.25_wp * ( w(k,j,i)   + w(k-1,j,i) +                  &
                                 w(k,j,i-1) + w(k-1,j,i-1) )
                IF ( wkomp > 0.0_wp )  THEN
                   tend(k,j,i) = tend(k,j,i) - wkomp *                         &
                                         ( u(k,j,i) - u(k-1,j,i) ) * ddzu(k)   &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 1 ) )
                ELSE
                   tend(k,j,i) = tend(k,j,i) - wkomp *                         &
                                         ( u(k+1,j,i) - u(k,j,i) ) * ddzu(k+1) &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 1 ) )
                ENDIF

             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE advec_u_up


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE advec_u_up_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  ddzu, tend, u, v, w

       USE control_parameters,                                                 &
           ONLY:  u_gtrans, v_gtrans

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nzb, nzt, wall_flags_0

       USE kinds


       IMPLICIT NONE

       INTEGER(iwp) ::  i !<
       INTEGER(iwp) ::  j !<
       INTEGER(iwp) ::  k !<

       REAL(wp) ::  ukomp !<
       REAL(wp) ::  vkomp !<
       REAL(wp) ::  wkomp !<


       DO  k = nzb+1, nzt
!
!--       x-direction
          ukomp = u(k,j,i) - u_gtrans
          IF ( ukomp > 0.0_wp )  THEN
             tend(k,j,i) = tend(k,j,i) - ukomp *                               &
                                         ( u(k,j,i) - u(k,j,i-1) ) * ddx       &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 1 ) )
          ELSE
             tend(k,j,i) = tend(k,j,i) - ukomp *                               &
                                         ( u(k,j,i+1) - u(k,j,i) ) * ddx       &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 1 ) )
          ENDIF
!
!--       y-direction
          vkomp = 0.25_wp * ( v(k,j,i) + v(k,j+1,i) + v(k,j,i-1) + v(k,j+1,i-1) &
                         ) - v_gtrans
          IF ( vkomp > 0.0_wp )  THEN
             tend(k,j,i) = tend(k,j,i) - vkomp *                               &
                                         ( u(k,j,i) - u(k,j-1,i) ) * ddy       &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 1 ) )
          ELSE
             tend(k,j,i) = tend(k,j,i) - vkomp *                               &
                                         ( u(k,j+1,i) - u(k,j,i) ) * ddy       &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 1 ) )
          ENDIF
!
!--       z-direction
          wkomp = 0.25_wp * ( w(k,j,i) + w(k-1,j,i) + w(k,j,i-1) + w(k-1,j,i-1) )
          IF ( wkomp > 0.0_wp )  THEN
             tend(k,j,i) = tend(k,j,i) - wkomp *                               &
                                         ( u(k,j,i) - u(k-1,j,i) ) * ddzu(k)   &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 1 ) )
          ELSE
             tend(k,j,i) = tend(k,j,i) - wkomp *                               &
                                         ( u(k+1,j,i) - u(k,j,i) ) * ddzu(k+1) &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 1 ) )
          ENDIF

       ENDDO

    END SUBROUTINE advec_u_up_ij

 END MODULE advec_u_up_mod
