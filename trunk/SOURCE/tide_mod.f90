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
! L Conlon, 20181024
! Initial version
!
!
!------------------------------------------------------------------------------!
 MODULE tide_mod

   USE arrays_3d,                                                              &
      ONLY:  u, v, u_tide, v_tide, zu, zw, dzu

   USE control_parameters,                                                     &
      ONLY:  g,  message_string, simulated_time

   USE indices,                                                                &
      ONLY:  nxlg, nxrg, nysg, nyng, nzb, nzt

   USE constants,                                                              &
      ONLY:  pi

   USE kinds

   IMPLICIT NONE



   PRIVATE
   PUBLIC init_tide, tide_simple

 CONTAINS




!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize tide according to the control parameters
!------------------------------------------------------------------------------!
   SUBROUTINE init_tide
      USE control_parameters
      IMPLICIT NONE

      INTEGER(iwp) ::  i
      INTEGER(iwp) ::  j
      INTEGER(iwp) ::  k

      ALLOCATE( u_tide(nzb:nzt+1), v_tide(nzb:nzt+1) )

         IF ( tide_amp == -9999999.9_wp .OR. tide_dir == -9999999.9_wp)  THEN
            WRITE( message_string, * )  'either tide_amp or tide_dir ',    &
               'is not set but required if tide == T'
         ENDIF

      u_tide = 0.0_wp
      v_tide = 0.0_wp
!
!--   Compute tide
         CALL tide_simple
!
!--   Update initial condition for u and v
      DO  i = nxlg, nxrg
         DO  j = nysg, nyng
            DO  k = nzb, nzt
               u(k,j,i) = u(k,j,i) - u_tide(k)
               v(k,j,i) = v(k,j,i) -v_tide(k)
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE init_tide


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute M2 tide
!------------------------------------------------------------------------------!
   SUBROUTINE tide_simple

      USE control_parameters

      IMPLICIT NONE

      INTEGER(iwp) ::  k
      REAL(wp)     ::  kdz, xcomp, ycomp, d2r


!--   tide direction
      d2r = pi / 180.0_wp
      xcomp = COS( tide_dir * d2r )
      ycomp = SIN( tide_dir * d2r )

! Calculate u and v components of M2 tide
        DO  k = nzt, nzb+1, -1
         u_tide(k) = xcomp*tide_amp*cos(.0001404*simulated_time)
         v_tide(k) = ycomp*tide_amp*cos(.0001404*simulated_time)
      ENDDO

      u_tide(nzt+1) = u_tide(nzt)
      v_tide(nzt+1) = v_tide(nzt)
      u_tide(nzb) = u_tide(nzb+1)
      v_tide(nzb) = v_tide(nzb+1)

   END SUBROUTINE tide_simple



 END MODULE tide_mod
