!> @file stokes_force_mod.f90
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
!
! Qing Li, 20180925
! Initial version
!
!
! Description:
! ------------
!> Computation of all Stokes force terms in the momentum equations
!> (Stokes-vortex force and Stokes-Coriolis force) and tracer euqations
!> (Stokes-advection)
!------------------------------------------------------------------------------!
 MODULE stokes_force_mod


    PRIVATE
    PUBLIC stokes_force_uvw, stokes_force_s

    INTERFACE stokes_force_uvw
       MODULE PROCEDURE stokes_force_uvw
       MODULE PROCEDURE stokes_force_uvw_ij
    END INTERFACE stokes_force_uvw

    INTERFACE stokes_force_s
       MODULE PROCEDURE stokes_force_s
       MODULE PROCEDURE stokes_force_s_ij
    END INTERFACE stokes_force_s

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Stokes-votex force and Stokes-Coriolis force in momentum equations
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE stokes_force_uvw( component )

       USE arrays_3d,                                                          &
           ONLY:  tend, u, v, u_stk, v_stk

       USE control_parameters,                                                 &
           ONLY:  f, fs, message_string

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxl, nxlu, nxr, nyn, nys, nysv, nzb, nzt, wall_flags_0

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  component  !< u-, v-, and w-component
       INTEGER(iwp) ::  i          !< running index x direction
       INTEGER(iwp) ::  j          !< running index y direction
       INTEGER(iwp) ::  k          !< running index z direction

       REAL(wp)     ::  flag       !< flag to mask topography
!
!--    Compute Stokes forces for the three velocity components
       SELECT CASE ( component )

!
!--       u-component
          CASE ( 1 )
             DO  i = nxlu, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
!
!--                   Predetermine flag to mask topography
                      flag = MERGE( 1.0_wp, 0.0_wp,                            &
                                    BTEST( wall_flags_0(k,j,i), 1 ) )
!
!--                   Stokes-Coriolis force
                      tend(k,j,i) = tend(k,j,i) + f * v_stk(k) * flag
!
!--                   Stokes-vortex force
                      tend(k,j,i) = tend(k,j,i) + v_stk(k) * (                 &
                                    0.5_wp * ( v(k,j,i) - v(k,j,i-1) +         &
                                    v(k,j+1,i) - v(k,j+1,i-1) ) * ddx -        &
                                    0.5_wp *( u(k,j+1,i) - u(k,j-1,i) ) * ddy  &
                                                             ) * flag
                   ENDDO
                ENDDO
             ENDDO

!
!--       v-component
          CASE ( 2 )
             DO  i = nxl, nxr
                DO  j = nysv, nyn
                   DO  k = nzb+1, nzt
!
!--                   Predetermine flag to mask topography
                      flag = MERGE( 1.0_wp, 0.0_wp,                            &
                                    BTEST( wall_flags_0(k,j,i), 2 ) )
!
!--                   Stokes-Coriolis force
                      tend(k,j,i) = tend(k,j,i) - f * u_stk(k) * flag
!
!--                   Stokes-vortex force
                      tend(k,j,i) = tend(k,j,i) + u_stk(k) * (                 &
                                    0.5_wp * ( u(k,j,i) - u(k,j-1,i) +         &
                                    u(k,j,i+1) - u(k,j-1,i+1) ) * ddy -        &
                                    0.5_wp *( v(k,j,i+1) - v(k,j,i-1) ) * ddx  &
                                                             ) * flag
                   ENDDO
                ENDDO
             ENDDO

!
!--       w-component
          CASE ( 3 )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
!
!--                   Predetermine flag to mask topography
                      flag = MERGE( 1.0_wp, 0.0_wp,                            &
                                    BTEST( wall_flags_0(k,j,i), 3 ) )
!
!--                   Stokes-Coriolis force
                      tend(k,j,i) = tend(k,j,i) + fs * u_stk(k) * flag
                   ENDDO
                ENDDO
             ENDDO

          CASE DEFAULT

             WRITE( message_string, * ) ' wrong component: ', component
             CALL message( 'Stokes_force_uvw', 'PA0601', 1, 2, 0, 6, 0 )

      END SELECT

    END SUBROUTINE stokes_force_uvw


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Stokes-votex force and Stokes-Coriolis force in momentum equations
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE stokes_force_uvw_ij( i, j, component )

       USE arrays_3d,                                                          &
           ONLY:  tend, u, v, u_stk, v_stk

       USE control_parameters,                                                 &
           ONLY:  f, fs, message_string

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nzb, nzt, wall_flags_0

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  component  !< u-, v-, and w-component
       INTEGER(iwp) ::  i          !< running index x direction
       INTEGER(iwp) ::  j          !< running index y direction
       INTEGER(iwp) ::  k          !< running index z direction

       REAL(wp)     ::  flag       !< flag to mask topography

!
!--    Compute Stokes forces for the three velocity components
       SELECT CASE ( component )

!
!--       u-component
          CASE ( 1 )
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 1 ) )
!
!--             Stokes-Coriolis force
                tend(k,j,i) = tend(k,j,i) + f * v_stk(k) * flag
!
!--             Stokes-vortex force
                tend(k,j,i) = tend(k,j,i) + v_stk(k) * (                       &
                              0.5_wp * ( v(k,j,i) - v(k,j,i-1) +               &
                              v(k,j+1,i) - v(k,j+1,i-1) ) * ddx -              &
                              0.5_wp *( u(k,j+1,i) - u(k,j-1,i) ) * ddy        &
                                                       ) * flag
             ENDDO

!
!--       v-component
          CASE ( 2 )
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 2 ) )
!
!--             Stokes-Coriolis force
                tend(k,j,i) = tend(k,j,i) - f * u_stk(k) * flag
!
!--             Stokes-vortex force
                tend(k,j,i) = tend(k,j,i) + u_stk(k) * (                       &
                              0.5_wp * ( u(k,j,i) - u(k,j-1,i) +               &
                              u(k,j,i+1) - u(k,j-1,i+1) ) * ddy -              &
                              0.5_wp *( v(k,j,i+1) - v(k,j,i-1) ) * ddx        &
                                                       ) * flag
             ENDDO

!
!--       w-component
          CASE ( 3 )
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 3 ) )
!
!--             Stokes-Coriolis force
                tend(k,j,i) = tend(k,j,i) + fs * u_stk(k) * flag
             ENDDO

          CASE DEFAULT

             WRITE( message_string, * ) ' wrong component: ', component
             CALL message( 'stokes_force_uvw', 'PA0601', 1, 2, 0, 6, 0 )

       END SELECT

    END SUBROUTINE stokes_force_uvw_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Stokes-advection term in tracer equations
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE stokes_force_s( sk )

       USE arrays_3d,                                                          &
           ONLY:  tend, u_stk, v_stk

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb,             &
                  nzt, wall_flags_0

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i          !< running index x direction
       INTEGER(iwp) ::  j          !< running index y direction
       INTEGER(iwp) ::  k          !< running index z direction

       REAL(wp)     ::  flag       !< flag to mask topography

#if defined( __nopointer )
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  sk !<
#else
       REAL(wp), DIMENSION(:,:,:), POINTER ::  sk
#endif

!--    Compute Stokes-advection term for the tracer equation
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp,                                  &
                              BTEST( wall_flags_0(k,j,i), 0 ) )
!
!--             Stokes-advection term
                tend(k,j,i) = tend(k,j,i) - 0.5_wp * (                         &
                              u_stk(k) * ( sk(k,j,i+1) - sk(k,j,i-1) ) * ddx + &
                              v_stk(k) * ( sk(k,j+1,i) - sk(k,j-1,i) ) * ddy   &
                                                     ) * flag
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE stokes_force_s


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Stokes-advection term in tracer equations
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE stokes_force_s_ij( i, j, sk )

       USE arrays_3d,                                                          &
           ONLY:  tend, u_stk, v_stk

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nyng, nysg, nzb, nzt, wall_flags_0

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i          !< running index x direction
       INTEGER(iwp) ::  j          !< running index y direction
       INTEGER(iwp) ::  k          !< running index z direction

       REAL(wp)     ::  flag       !< flag to mask topography

#if defined( __nopointer )
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  sk !<
#else
       REAL(wp), DIMENSION(:,:,:), POINTER ::  sk
#endif

!--    Compute Stokes-advection term for the tracer equation
       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp,                                        &
                        BTEST( wall_flags_0(k,j,i), 0 ) )
!
!--       Stokes-advection term
          tend(k,j,i) = tend(k,j,i) - 0.5_wp * (                               &
                        u_stk(k) * ( sk(k,j,i+1) - sk(k,j,i-1) ) * ddx +       &
                        v_stk(k) * ( sk(k,j+1,i) - sk(k,j-1,i) ) * ddy         &
                                               ) * flag
       ENDDO

    END SUBROUTINE stokes_force_s_ij

 END MODULE stokes_force_mod
