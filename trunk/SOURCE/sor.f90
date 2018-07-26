!> @file sor.f90
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
! $Id: sor.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! - Change in file header (GPL part)
! - Large-scale forcing implemented (MS)
! 
! 2101 2017-01-05 16:42:31Z suehring
!
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1762 2016-02-25 12:31:13Z hellstea
! Introduction of nested domain feature
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
! old module precision_kind is removed,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! Revision 1.1  1997/08/11 06:25:56  raasch
! Initial revision
!
!
! Description:
! ------------
!> Solve the Poisson-equation with the SOR-Red/Black-scheme.
!------------------------------------------------------------------------------!
 SUBROUTINE sor( d, ddzu, ddzw, p )

    USE arrays_3d,                                                             &
        ONLY:  rho_air, rho_air_zw

    USE grid_variables,                                                        &
        ONLY:  ddx2, ddy2

    USE indices,                                                               &
        ONLY:  nbgp, nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nz, nzb, nzt

    USE kinds

    USE control_parameters,                                                    &
        ONLY:  bc_lr_cyc, bc_ns_cyc, force_bound_l, force_bound_n,             &
               force_bound_r, force_bound_s, ibc_p_b, ibc_p_t, inflow_l,       &
               inflow_n, inflow_r, inflow_s, nest_bound_l, nest_bound_n,       &
               nest_bound_r, nest_bound_s, n_sor, omega_sor, outflow_l,        &
               outflow_n, outflow_r, outflow_s

    IMPLICIT NONE

    INTEGER(iwp) ::  i              !<
    INTEGER(iwp) ::  j              !<
    INTEGER(iwp) ::  k              !<
    INTEGER(iwp) ::  n              !<
    INTEGER(iwp) ::  nxl1           !<
    INTEGER(iwp) ::  nxl2           !<
    INTEGER(iwp) ::  nys1           !<
    INTEGER(iwp) ::  nys2           !<

    REAL(wp)     ::  ddzu(1:nz+1)   !<
    REAL(wp)     ::  ddzw(1:nzt+1)  !<

    REAL(wp)     ::  d(nzb+1:nzt,nys:nyn,nxl:nxr)      !<
    REAL(wp)     ::  p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  !<

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  f1         !<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  f2         !<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  f3         !<

    ALLOCATE( f1(1:nz), f2(1:nz), f3(1:nz) )

!
!-- Compute pre-factors.
    DO  k = 1, nz
         f2(k) = ddzu(k+1) * ddzw(k) * rho_air_zw(k)
         f3(k) = ddzu(k)   * ddzw(k) * rho_air_zw(k-1)
         f1(k) = 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k) + f2(k) + f3(k)
    ENDDO

!
!-- Limits for RED- and BLACK-part.
    IF ( MOD( nxl , 2 ) == 0 )  THEN
       nxl1 = nxl
       nxl2 = nxl + 1
    ELSE
       nxl1 = nxl + 1
       nxl2 = nxl
    ENDIF
    IF ( MOD( nys , 2 ) == 0 )  THEN
       nys1 = nys
       nys2 = nys + 1
    ELSE
       nys1 = nys + 1
       nys2 = nys
    ENDIF

    DO  n = 1, n_sor

!
!--    RED-part
       DO  i = nxl1, nxr, 2
          DO  j = nys2, nyn, 2
             DO  k = nzb+1, nzt
                p(k,j,i) = p(k,j,i) + omega_sor / f1(k) * (            &
                           rho_air(k) * ddx2 * ( p(k,j,i+1) + p(k,j,i-1) ) +   &
                           rho_air(k) * ddy2 * ( p(k,j+1,i) + p(k,j-1,i) ) +   &
                           f2(k) * p(k+1,j,i)                              +   &
                           f3(k) * p(k-1,j,i)                              -   &
                           d(k,j,i)                                        -   &
                           f1(k) * p(k,j,i)           )
             ENDDO
          ENDDO
       ENDDO

       DO  i = nxl2, nxr, 2
          DO  j = nys1, nyn, 2
             DO  k = nzb+1, nzt
                p(k,j,i) = p(k,j,i) + omega_sor / f1(k) * (                    &
                           rho_air(k) * ddx2 * ( p(k,j,i+1) + p(k,j,i-1) ) +   &
                           rho_air(k) * ddy2 * ( p(k,j+1,i) + p(k,j-1,i) ) +   &
                           f2(k) * p(k+1,j,i)                              +   &
                           f3(k) * p(k-1,j,i)                              -   &
                           d(k,j,i)                                        -   &
                           f1(k) * p(k,j,i)           )
             ENDDO
          ENDDO
       ENDDO

!
!--    Exchange of boundary values for p.
       CALL exchange_horiz( p, nbgp )

!
!--    Horizontal (Neumann) boundary conditions in case of non-cyclic boundaries
       IF ( .NOT. bc_lr_cyc )  THEN
          IF ( inflow_l      .OR.  outflow_l  .OR.                                 &
               nest_bound_l  .OR.  force_bound_l )  p(:,:,nxl-1) = p(:,:,nxl)
          IF ( inflow_r      .OR.  outflow_r  .OR.                                 &
               nest_bound_r  .OR.  force_bound_r )  p(:,:,nxr+1) = p(:,:,nxr)
       ENDIF
       IF ( .NOT. bc_ns_cyc )  THEN
          IF ( inflow_n      .OR.  outflow_n  .OR.                                 &
               nest_bound_n  .OR.  force_bound_n )  p(:,nyn+1,:) = p(:,nyn,:)
          IF ( inflow_s      .OR.  outflow_s  .OR.                                 &
               nest_bound_s  .OR.  force_bound_s )  p(:,nys-1,:) = p(:,nys,:)
       ENDIF

!
!--    BLACK-part
       DO  i = nxl1, nxr, 2
          DO  j = nys1, nyn, 2
             DO  k = nzb+1, nzt
                p(k,j,i) = p(k,j,i) + omega_sor / f1(k) * (            &
                           rho_air(k) * ddx2 * ( p(k,j,i+1) + p(k,j,i-1) ) +   &
                           rho_air(k) * ddy2 * ( p(k,j+1,i) + p(k,j-1,i) ) +   &
                           f2(k) * p(k+1,j,i)                              +   &
                           f3(k) * p(k-1,j,i)                              -   &
                           d(k,j,i)                                        -   &
                           f1(k) * p(k,j,i)           )
             ENDDO
          ENDDO
       ENDDO

       DO  i = nxl2, nxr, 2
          DO  j = nys2, nyn, 2
             DO  k = nzb+1, nzt
                p(k,j,i) = p(k,j,i) + omega_sor / f1(k) * (            &
                           rho_air(k) * ddx2 * ( p(k,j,i+1) + p(k,j,i-1) ) +   &
                           rho_air(k) * ddy2 * ( p(k,j+1,i) + p(k,j-1,i) ) +   &
                           f2(k) * p(k+1,j,i)                              +   &
                           f3(k) * p(k-1,j,i)                              -   &
                           d(k,j,i)                                        -   &
                           f1(k) * p(k,j,i)           )
             ENDDO
          ENDDO
       ENDDO

!
!--    Exchange of boundary values for p.
       CALL exchange_horiz( p, nbgp )

!
!--    Boundary conditions top/bottom.
!--    Bottom boundary
       IF ( ibc_p_b == 1 )  THEN       !       Neumann
          p(nzb,:,:) = p(nzb+1,:,:)
       ELSE                            !       Dirichlet
          p(nzb,:,:) = 0.0_wp
       ENDIF

!
!--    Top boundary
       IF ( ibc_p_t == 1 )  THEN                 !  Neumann
          p(nzt+1,:,:) = p(nzt,:,:)
       ELSE                      !  Dirichlet
          p(nzt+1,:,:) = 0.0_wp
       ENDIF

!
!--    Horizontal (Neumann) boundary conditions in case of non-cyclic boundaries
       IF ( .NOT. bc_lr_cyc )  THEN
          IF ( inflow_l      .OR.  outflow_l  .OR.                             &
               nest_bound_l  .OR.  force_bound_l )  p(:,:,nxl-1) = p(:,:,nxl)
          IF ( inflow_r      .OR.  outflow_r  .OR.                             &
               nest_bound_r  .OR.  force_bound_r )  p(:,:,nxr+1) = p(:,:,nxr)
       ENDIF
       IF ( .NOT. bc_ns_cyc )  THEN
          IF ( inflow_n      .OR.  outflow_n  .OR.                             &
               nest_bound_n  .OR.  force_bound_n )  p(:,nyn+1,:) = p(:,nyn,:)
          IF ( inflow_s      .OR.  outflow_s  .OR.                             &
               nest_bound_s  .OR.  force_bound_s )  p(:,nys-1,:) = p(:,nys,:)
       ENDIF


    ENDDO

    DEALLOCATE( f1, f2, f3 )

 END SUBROUTINE sor
