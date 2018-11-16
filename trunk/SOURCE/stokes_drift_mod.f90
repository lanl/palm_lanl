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

   INTEGER(iwp), PARAMETER :: FROMUSDELTA = 1   !< Stokes drift from surface value and decay depth
   INTEGER(iwp), PARAMETER :: FROMSPECDHH85 = 2 !< Stokes drift from DHH85 spectrum
   INTEGER(iwp), PARAMETER :: COASTAL = 3

   PRIVATE

   PUBLIC init_stokes_drift, stokes_drift_check_parameters

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!>
!------------------------------------------------------------------------------!
   SUBROUTINE stokes_drift_check_parameters

      USE control_parameters,                                                  &
         ONLY: stokes_drift_method, u0_stk, v0_stk, d_stk,                     &
               wind_speed, wind_dir, wave_age, wave_length

      IMPLICIT NONE
!
!--   Check Stokes drift method
      SELECT CASE ( stokes_drift_method )
      CASE ( -9999 )
         WRITE( message_string, * )  'stokes_drift_method ',                   &
               'is not set but required if stokes_force == .TRUE.'
         CALL message( 'stokes_drift_check_parameters', 'PA0603', 1, 2, 0, 6, 0)
      CASE ( 1 )
         IF ( u0_stk == -9999999.9_wp .OR. v0_stk == -9999999.9_wp .OR.        &
              d_stk == -9999999.9_wp ) THEN
            WRITE( message_string, * )  'either u0_stk, v0_skt, or d_stk ',    &
               'is not set but required if stokes_drift_method == 1'
            CALL message( 'stokes_drift_check_parameters',                     &
                          'PA0603', 1, 2, 0, 6, 0)
         ELSEIF ( d_stk < 0.0_wp ) THEN
            WRITE( message_string, * )  'd_stk = ',                            &
                   d_stk, ' must be > 0 (m)'
            CALL message( 'stokes_drift_check_parameters',                     &
                          'PA0603', 1, 2, 0, 6, 0)
         ENDIF
      CASE ( 2 )
         IF ( wind_speed == -9999999.9_wp .OR.                                 &
              wind_dir == -9999999.9_wp .OR.                                   &
              wave_age == -9999999.9_wp ) THEN
            WRITE( message_string, * )  'either wind_speed, wind_dir, or ',    &
               'wave_age is not set but required if stokes_drift_method == 2'
            CALL message( 'stokes_drift_check_parameters',                     &
                          'PA0603', 1, 2, 0, 6, 0)
         ELSEIF ( wave_age < 0.2_wp .OR. wave_age > 1.2_wp ) THEN
            WRITE( message_string, * )  'wave_age = ',                         &
                   wave_age, ' must be between 0.2 and 1.2'
            CALL message( 'stokes_drift_check_parameters',                     &
                          'PA0603', 1, 2, 0, 6, 0)
         ENDIF
       CASE ( 3 )
         IF ( wind_dir == -9999999.9_wp .OR. wave_length == -9999999.9_wp ) THEN
            WRITE( message_string, * )  'either wind_dir or wave_length ',    &
               'is not set but required if stokes_drift_method == 3'
            CALL message( 'stokes_drift_check_parameters',                     &
                          'PA0603', 1, 2, 0, 6, 0)
         ELSEIF ( wave_length < 0.0_wp ) THEN
            WRITE( message_string, * )  'wave_length = ',                            &
                   wave_length, ' must be > 0 (m)'
            CALL message( 'stokes_drift_check_parameters',                     &
                          'PA0603', 1, 2, 0, 6, 0)

         ENDIF
      CASE DEFAULT
         WRITE( message_string, * )  'invalid stokes_drift_mehtod = ',         &
                stokes_drift_method, ', must be 1, 2, or 3'
         CALL message( 'stokes_drift_check_parameters', 'PA0603', 1, 2, 0, 6, 0)
      END SELECT

   END SUBROUTINE stokes_drift_check_parameters


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
!
!--   Compute Stokes drift
      SELECT CASE ( stokes_drift_method )
      CASE ( FROMUSDELTA )
         CALL stokes_drift_usdelta
      CASE ( FROMSPECDHH85 )
         CALL stokes_drift_spec_dhh85
      CASE (COASTAL)
         CALL stokes_drift_coastal
      CASE DEFAULT
         WRITE( message_string, * )  'invalid stokes_drift_mehtod = ',         &
                stokes_drift_method, ', must be 1, 2, or 3'
         CALL message( 'init_stokes_drift', 'PA0602', 1, 2, 0, 6, 0 )
      END SELECT
!
!--   Update initial condition for u and v
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
!
!--   Stokes drift averaged over the grid cell
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

!
!--   Frequency range
      REAL(wp), PARAMETER     :: min_omega = 0.1_wp, max_omega = 10.0_wp
      INTEGER(iwp), PARAMETER :: nomega = 1000

      INTEGER(iwp) :: i, k
      REAL(wp)     :: d2r, xcomp, ycomp, domega, sd_omega, tmp
!
!--   Wind direction
      d2r = pi / 180.0_wp
      xcomp = COS( wind_dir * d2r )
      ycomp = SIN( wind_dir * d2r )
!
!--   Integral over frequency
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
!> Compute Stokes drift profile following Sinha 2015 for coastal
!  applications
!------------------------------------------------------------------------------!
   SUBROUTINE stokes_drift_coastal

      USE control_parameters,                                                  &
         ONLY:  wind_dir, wave_length

      IMPLICIT NONE


      INTEGER(iwp) :: k
      REAL(wp)     :: d2r, xcomp, ycomp, domega, sd_omega, tmp, w_num, depth_s
      REAL(wp) :: w_num_d, half_depth_stokes, num_depths, depth, dzdepth
      REAL(wp), dimension (nzt) :: deps
!

!     Depth
      num_depths=SIZE(zw)
      depth=MAXVAL(abs(zw))
      half_depth_stokes=depth/2
      deps=-((abs(zw*2)/depth)-1)


!     Wave number
      w_num=2*pi/ wave_length
      w_num_d=w_num*half_depth_stokes


!--   Wind direction
      d2r = pi / 180.0_wp
      xcomp = COS( wind_dir * d2r )
      ycomp = SIN( wind_dir * d2r )

      DO  k = nzb+1, nzt, 1
           u_stk(k)=(cosh(2*w_num_d*(deps(k)+1))/(2*sinh(2*w_num_d)**2))*xcomp
           v_stk(k)=(cosh(2*w_num_d*(deps(k)+1))/(2*sinh(2*w_num_d)**2))*ycomp
      ENDDO


      u_stk(nzt+1) = u_stk(nzt)
      v_stk(nzt+1) = v_stk(nzt)
      u_stk(nzb) = u_stk(nzb+1)
      v_stk(nzb) = v_stk(nzb+1)


   END SUBROUTINE stokes_drift_coastal


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
!
!--   DHH 85 spectrum
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
!
!--   Stokes drift integral kernel
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
