!> @file tide_eqs_mod.f90
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
! L Conlon, 20181031
! Initial version
!
!
!------------------------------------------------------------------------------!
 MODULE tide_eqs_mod


    PRIVATE
    PUBLIC tide_eqs_uvw, tide_eqs_s, tide_production_e,              &
           tide_pressure_head

    INTERFACE tide_eqs_uvw
       MODULE PROCEDURE tide_eqs_uvw
       MODULE PROCEDURE tide_eqs_uvw_ij
    END INTERFACE tide_eqs_uvw

    INTERFACE tide_eqs_s
       MODULE PROCEDURE tide_eqs_s
       MODULE PROCEDURE tide_eqs_s_ij
    END INTERFACE tide_eqs_s

    INTERFACE tide_production_e
       MODULE PROCEDURE tide_production_e
       MODULE PROCEDURE tide_production_e_ij
    END INTERFACE tide_production_e

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Tide in momentum equations (coriolis)
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE tide_eqs_uvw( component )

       USE arrays_3d,                                                          &
           ONLY:  tend, u, v, u_tide, v_tide

       USE control_parameters,                                                 &
           ONLY:  f, fs, message_string

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxl, nxlu, nxr, nyn, nys, nysv, nzb, nzt, wall_flags_0

       USE kinds

       USE tide_mod
       IMPLICIT NONE

       INTEGER(iwp) ::  component  !< u-, v-, and w-component
       INTEGER(iwp) ::  i          !< running index x direction
       INTEGER(iwp) ::  j          !< running index y direction
       INTEGER(iwp) ::  k          !< running index z direction

       REAL(wp)     ::  flag       !< flag to mask topography
!

!--    Compute tidal forces for the three velocity components
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
!--                   Tidal Coriolis force
                      tend(k,j,i) = tend(k,j,i) + f * v_tide(k) * flag
!
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
!--                   Tidal Coriolis force
                      tend(k,j,i) = tend(k,j,i) - f * u_tide(k) * flag
!
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
!--                   Tidal Coriolis force
                      tend(k,j,i) = tend(k,j,i) + fs * u_tide(k) * flag
                   ENDDO
                ENDDO
             ENDDO

          CASE DEFAULT

             WRITE( message_string, * ) ' wrong component: ', component
             CALL message( 'tide_eqs_uvw', 'PA0601', 1, 2, 0, 6, 0 )

      END SELECT

    END SUBROUTINE tide_eqs_uvw


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Tidal Coriolis force in momentum equations
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE tide_eqs_uvw_ij( i, j, component )

       USE arrays_3d,                                                          &
           ONLY:  tend, u, v, u_tide, v_tide

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

!       CALL tide_simple


!
!--    Compute tidal forces for the three velocity components
       SELECT CASE ( component )

!
!--       u-component
          CASE ( 1 )
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 1 ) )
!
!--             Tidal Coriolis force
                tend(k,j,i) = tend(k,j,i) + f * v_tide(k) * flag
!
             ENDDO

!
!--       v-component
          CASE ( 2 )
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 2 ) )
!
!--             Tidal Coriolis force
                tend(k,j,i) = tend(k,j,i) - f * u_tide(k) * flag
!
             ENDDO

!
!--       w-component
          CASE ( 3 )
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 3 ) )
!
!--             Tidal Coriolis force
                tend(k,j,i) = tend(k,j,i) + fs * u_tide(k) * flag
             ENDDO

          CASE DEFAULT

             WRITE( message_string, * ) ' wrong component: ', component
             CALL message( 'tide_eqs_uvw', 'PA0601', 1, 2, 0, 6, 0 )

       END SELECT

    END SUBROUTINE tide_eqs_uvw_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Tidal advection term in tracer equations
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE tide_eqs_s( sk )

       USE arrays_3d,                                                          &
           ONLY:  tend, u_tide, v_tide

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

!--    Compute advection term for the tracer equation
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp,                                  &
                              BTEST( wall_flags_0(k,j,i), 0 ) )
!
!--             tidal advection term
                tend(k,j,i) = tend(k,j,i) - 0.5_wp * (                         &
                              u_tide(k) * ( sk(k,j,i+1) - sk(k,j,i-1) ) * ddx + &
                              v_tide(k) * ( sk(k,j+1,i) - sk(k,j-1,i) ) * ddy   &
                                                     ) * flag
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE tide_eqs_s


!------------------------------------------------------------------------------!
! Description:
! ------------
!> tidal-advection term in tracer equations
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE tide_eqs_s_ij( i, j, sk )

       USE arrays_3d,                                                          &
           ONLY:  tend, u_tide, v_tide

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


!--    Compute tidal advection term for the tracer equation
       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp,                                        &
                        BTEST( wall_flags_0(k,j,i), 0 ) )
!
!--       tidal-advection term
          tend(k,j,i) = tend(k,j,i) - 0.5_wp * (                               &
                        u_tide(k) * ( sk(k,j,i+1) - sk(k,j,i-1) ) * ddx +       &
                        v_tide(k) * ( sk(k,j+1,i) - sk(k,j-1,i) ) * ddy         &
                                               ) * flag
       ENDDO

    END SUBROUTINE tide_eqs_s_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> tidal production term in TKE equations
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE tide_production_e

       USE arrays_3d,                                                          &
           ONLY:  tend, u, v, w, u_tide, v_tide, dd2zu, km

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxl, nxr, nyn, nys, nzb, nzt, wall_flags_0

       USE kinds


       IMPLICIT NONE

       INTEGER(iwp) ::  i          !< running index x direction
       INTEGER(iwp) ::  j          !< running index y direction
       INTEGER(iwp) ::  k          !< running index z direction

       REAL(wp)     ::  flag       !< flag to mask topography
       REAL(wp)     ::  dudz, dvdz, dwdx, dwdy


!--    Compute tidal advection term for the tracer equation
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp,                                  &
                              BTEST( wall_flags_0(k,j,i), 29 ) )
!
!--             Stokes-production term
                dudz = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) -                 &
                                   u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)
                dwdx = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) -                 &
                                   w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
                dvdz = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) -                 &
                                   v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)
                dwdy = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) -                 &
                                   w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
                tend(k,j,i) = tend(k,j,i) + km(k,j,i) * (                      &
                              ( u_tide(k+1) - u_tide(k-1) ) * dd2zu(k) *         &
                              ( dudz + dwdx ) +                                &
                              ( v_tide(k+1) - v_tide(k-1) ) * dd2zu(k) *         &
                              ( dvdz + dwdy )           ) * flag
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE tide_production_e


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE tide_production_e_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  tend, u, v, w, u_tide, v_tide, dd2zu, km

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxl, nxr, nyn, nys, nzb, nzt, wall_flags_0

       USE kinds


       IMPLICIT NONE

       INTEGER(iwp) ::  i          !< running index x direction
       INTEGER(iwp) ::  j          !< running index y direction
       INTEGER(iwp) ::  k          !< running index z direction

       REAL(wp)     ::  flag       !< flag to mask topography
       REAL(wp)     ::  dudz, dvdz, dwdx, dwdy

!--    Compute Stokes-advection term for the tracer equation
       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp,                                  &
                        BTEST( wall_flags_0(k,j,i), 29 ) )
!
!--       Stokes-production term
          dudz = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) -                 &
                             u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)
          dwdx = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) -                 &
                             w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
          dvdz = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) -                 &
                             v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)
          dwdy = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) -                 &
                             w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
          tend(k,j,i) = tend(k,j,i) + km(k,j,i) * (                      &
                        ( u_tide(k+1) - u_tide(k-1) ) * dd2zu(k) *            &
                        ( dudz + dwdx ) +                                &
                        ( v_tide(k+1) - v_tide(k-1) ) * dd2zu(k) *            &
                        ( dvdz + dwdy )           ) * flag
       ENDDO

    END SUBROUTINE tide_production_e_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!------------------------------------------------------------------------------!
    SUBROUTINE tide_pressure_head

       USE arrays_3d,                                                          &
           ONLY:  p, u, v, u_tide, v_tide, rho_air

       USE indices,                                                            &
           ONLY:  nxl, nxr, nyn, nys, nzb, nzt, wall_flags_0

       USE kinds


       IMPLICIT NONE

       INTEGER(iwp) ::  i          !< running index x direction
       INTEGER(iwp) ::  j          !< running index y direction
       INTEGER(iwp) ::  k          !< running index z direction

       REAL(wp)     ::  flag       !< flag to mask topography

! TODO: Check if the stokes pressure head is correctly scaled by rho_air,
!       as rho_air is used in pres.f90 to calculate p
!       <20180926, Qing Li> !

!
!--    Update perturbation pressure
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb, nzt+1
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp,                                  &
                           BTEST( wall_flags_0(k,j,i), 0 ) )
                p(k,j,i) = p(k,j,i) - 0.5_wp * rho_air(k) *                    &
                           ( ( u(k,j,i) + u(k,j,i+1) ) * u_tide(k) +            &
                             ( v(k,j,i) + v(k,j+1,i) ) * v_tide(k) +            &
                             u_tide(k)**2 + v_tide(k)**2                         &
                           ) * flag
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE tide_pressure_head

 END MODULE tide_eqs_mod
