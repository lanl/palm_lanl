!> @file stokes_drift_mod.f90
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
! Qing Li, 20180926
! Initial version
!
!
! Description:
! ------------
!> Subroutines to compute Stokes drift profile from various sources
!------------------------------------------------------------------------------!
 MODULE stokes_drift_mod


   PRIVATE
   ! PUBLIC stokes_drift_check_parameters, init_stokes_dirft
   PUBLIC init_stokes_drift

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Stokes-votex force and Stokes-Coriolis force in momentum equations
!> Call for all grid points
!------------------------------------------------------------------------------!
   ! SUBROUTINE stokes_drfit_check_parameters
   ! END SUBROUTINE stokes_drift_check_parameters


   SUBROUTINE init_stokes_drift

      USE arrays_3d,                                                           &
         ONLY:  u, v, u_stk, v_stk, zu, zw

      USE control_parameters,                                                  &
         ONLY:  g

      USE indices,                                                             &
         ONLY:  nxlg, nxrg, nysg, nyng, nzb, nzt

      USE kinds

      IMPLICIT NONE

      INTEGER(iwp) ::  i
      INTEGER(iwp) ::  j
      INTEGER(iwp) ::  k

      REAL(wp)     ::  u_stk0, v_stk0, d_stk

      ALLOCATE( u_stk(nzb:nzt+1), v_stk(nzb:nzt+1) )
      u_stk0 = 0.068_wp
      v_stk0 = 0.0_wp
      d_stk  = 1.0_wp/4.8_wp

      DO  k = nzt, nzb, -1
         u_stk(k) = u_stk0 * EXP( zu(k) * d_stk )
         v_stk(k) = v_stk0 * EXP( zu(k) * d_stk )
      ENDDO
      u_stk(nzt+1) = u_stk(nzt)
      v_stk(nzt+1) = v_stk(nzt)

      ! update initial condition for u and v
      DO  i = nxlg, nxrg
         DO  j = nysg, nyng
            DO  k = nzb, nzt
               u(k,j,i) = u(k,j,i) - u_stk(k)
               v(k,j,i) = v(k,j,i) - v_stk(k)
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE init_stokes_drift

 END MODULE stokes_drift_mod
