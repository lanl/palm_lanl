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

   USE arrays_3d,                                                              &
      ONLY:  u, v, u_stk, v_stk, zu, zw, dzu

   USE control_parameters,                                                     &
      ONLY:  g, stokes_drift_method, message_string

   USE indices,                                                                &
      ONLY:  nxlg, nxrg, nysg, nyng, nzb, nzt

   USE constants,                                                              &
      ONLY:  pi

   USE kinds

   IMPLICIT NONE

   INTEGER(iwp), PARAMETER :: NOTHING = 0       !< No Stokes drift
   INTEGER(iwp), PARAMETER :: FROMUSDELTA = 1   !< Stokes drift from surface value and decay depth
   INTEGER(iwp), PARAMETER :: FROMSPECDHH85 = 2 !< Stokes drift from DHH85 spectrum

   PRIVATE
   ! PUBLIC stokes_drift_check_parameters, init_stokes_dirft
   PUBLIC init_stokes_drift

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!>
!------------------------------------------------------------------------------!
   ! SUBROUTINE stokes_drfit_check_parameters
   ! END SUBROUTINE stokes_drift_check_parameters


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize Stokes drift profile according to the control parameters
!------------------------------------------------------------------------------!
   SUBROUTINE init_stokes_drift

      IMPLICIT NONE

      INTEGER(iwp) ::  i
      INTEGER(iwp) ::  j
      INTEGER(iwp) ::  k

      ALLOCATE( u_stk(nzb:nzt+1), v_stk(nzb:nzt+1) )

      u_stk = 0.0_wp
      v_stk = 0.0_wp

      ! compute Stokes drift
      SELECT CASE ( stokes_drift_method )
      CASE ( NOTHING )
         u_stk = 0.0_wp
         v_stk = 0.0_wp
      CASE ( FROMUSDELTA )
         CALL stokes_drift_usdelta
      CASE ( FROMSPECDHH85 )
         CALL stokes_drift_spec_dhh85
      CASE DEFAULT
         WRITE( message_string, * ) ' unknown method for Stokes drift: ', stokes_drift_method
         CALL message( 'init_stokes_drift', 'PA0602', 1, 2, 0, 6, 0 )
      END SELECT

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


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute grid cell-averaged Stokes drift profile from surface Stokes drift
!> and the exponential decay depth scale
!------------------------------------------------------------------------------!
   SUBROUTINE stokes_drift_usdelta

      USE control_parameters,                                                  &
         ONLY:  u0_stk, v0_stk, d_stk

      IMPLICIT NONE

      INTEGER(iwp) ::  k
      REAL(wp)     ::  kdz, dd_stk, tmp

!     Stokes drift averaged over the grid cell
      dd_stk = 1.0_wp / d_stk
      DO  k = nzt, nzb+1, -1
         kdz = 0.5 * dzu(k) * dd_stk
         IF ( kdz .LT. 10.0_wp ) THEN
             tmp = SINH(kdz) / kdz * EXP( zu(k) * dd_stk )
         ELSE
             tmp = EXP( zu(k) * dd_stk )
         ENDIF
         u_stk(k) = u0_stk * tmp
         v_stk(k) = v0_stk * tmp
      ENDDO
      u_stk(nzt+1) = u_stk(nzt)
      v_stk(nzt+1) = v_stk(nzt)
      u_stk(nzb) = u_stk(nzb+1)
      v_stk(nzb) = v_stk(nzb+1)

   END SUBROUTINE stokes_drift_usdelta


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute grid cell-averaged Stokes drift profile from the empirical wave
!> spectrum of Donelan et al., 1985
!------------------------------------------------------------------------------!
   SUBROUTINE stokes_drift_spec_dhh85

      USE control_parameters,                                                  &
         ONLY:  wind_dir

      IMPLICIT NONE

      REAL(wp), PARAMETER     :: min_omega = 0.1_wp, max_omega = 10.0_wp
      INTEGER(iwp), PARAMETER :: nomega = 1000

      INTEGER(iwp) :: i, k
      REAL(wp)     :: d2r, xcomp, ycomp, domega, sd_omega, tmp

      ! wind direction
      d2r = pi / 180.0_wp
      xcomp = COS( wind_dir * d2r )
      ycomp = SIN( wind_dir * d2r )
      ! integral over frequency
      domega = ( max_omega - min_omega ) / REAL(nomega, KIND=wp)
      DO  k = nzt, nzb+1, -1
         tmp = 0.0_wp
         sd_omega = min_omega + 0.5_wp * domega
         DO  i = 1, nomega
            tmp = tmp +                                                       &
               domega * stokes_drift_kernel_dhh85(sd_omega, zu(k), dzu(k))
            sd_omega = sd_omega + domega
         ENDDO
         u_stk(k) = xcomp * tmp
         v_stk(k) = ycomp * tmp
      ENDDO
      u_stk(nzt+1) = u_stk(nzt)
      v_stk(nzt+1) = v_stk(nzt)
      u_stk(nzb) = u_stk(nzb+1)
      v_stk(nzb) = v_stk(nzb+1)

   END SUBROUTINE stokes_drift_spec_dhh85


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Evaluate kernel of the Stokes integral for the Donelan et al., 1985 spectrum
!------------------------------------------------------------------------------!
   FUNCTION stokes_drift_kernel_dhh85(sd_omega, sd_z, sd_dz)

      USE control_parameters,                                                  &
         ONLY:  wind_speed, wave_age

      IMPLICIT NONE

      REAL(wp)  :: sd_omega, sd_z, sd_dz
      REAL(wp)  :: dhh_omega_p, dhh_alpha, dhh_sigma, dhh_gamma1, dhh_gamma2
      REAL(wp)  :: wave_spec, sd_filter, kdz, iwa
      REAL(wp)  :: stokes_drift_kernel_dhh85

      ! DHH 85 spectrum
      iwa = 1.0_wp / wave_age !< inverse wave age
      dhh_omega_p = g * iwa / wind_speed !< peak frequency
      dhh_alpha  = 0.006_wp * iwa**(0.55_wp)
      dhh_sigma  = 0.08_wp * ( 1.0_wp + 4.0_wp * wave_age**3 )
      IF ( iwa .LE. 1.0_wp) THEN
         dhh_gamma1 = 1.7_wp
      ELSE
         dhh_gamma1 = 1.7_wp + 6.0_wp * LOG10( iwa )
      ENDIF
      dhh_gamma2 = EXP( -0.5_wp * ( sd_omega - dhh_omega_p )**2 /              &
                       dhh_sigma**2 / dhh_omega_p**2 )
      wave_spec  = dhh_alpha * g**2 / (dhh_omega_p * sd_omega**4 ) *           &
                   EXP( -( dhh_omega_p / sd_omega )**4 ) *                     &
                   dhh_gamma1**dhh_gamma2
      ! Stokes drift integral kernel
      kdz = sd_omega**2 * sd_dz / g
      IF ( kdz .LT. 10.0_wp ) THEN
         sd_filter = SINH(kdz) / kdz
      ELSE
         sd_filter = 1.0_wp
      ENDIF
      stokes_drift_kernel_dhh85 = 2.0_wp * ( wave_spec * sd_omega**3 ) *       &
                     sd_filter * EXP( 2.0_wp * sd_omega**2 * sd_z / g ) / g

      RETURN

   END FUNCTION stokes_drift_kernel_dhh85


 END MODULE stokes_drift_mod
