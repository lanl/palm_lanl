!> @file diffusion_s.f90
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
! ------------------
!
!
! Former revisions:
! -----------------
! $Id: diffusion_s.f90 2759 2018-01-17 16:24:59Z suehring $
! Major bugfix, horizontal diffusion at vertical surfaces corrected.
!
! 2718 2018-01-02 08:49:38Z maronga
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
! 1691 2015-10-26 16:17:44Z maronga
! Formatting corrections.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1374 2014-04-25 12:55:07Z raasch
! missing variables added to ONLY list
!
! 1340 2014-03-25 19:45:13Z kanani
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
! openacc loop and loop vector clauses removed
!
! 1128 2013-04-12 06:19:32Z raasch
! loop index bounds in accelerator version replaced by i_left, i_right, j_south,
! j_north
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1015 2012-09-27 09:23:24Z raasch
! accelerator version (*_acc) added
!
! 1010 2012-09-20 07:59:54Z raasch
! cpp switch __nopointer added for pointer free version
!
! 1001 2012-09-13 14:08:46Z raasch
! some arrays comunicated by module instead of parameter list
!
! Revision 1.1  2000/04/13 14:54:02  schroeter
! Initial revision
!
!
! Description:
! ------------
!> Diffusion term of scalar quantities (temperature and water content)
!------------------------------------------------------------------------------!
 MODULE diffusion_s_mod


    PRIVATE
    PUBLIC diffusion_s

    INTERFACE diffusion_s
       MODULE PROCEDURE diffusion_s
    END INTERFACE diffusion_s

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE diffusion_s( s, s_flux_t, s_flux_solar_t)

       USE arrays_3d,                                                          &
           ONLY:  dzw, ddzu, ddzw, kh, tend, drho_air, rho_air_zw, solar3d

       USE control_parameters,                                                 &
           ONLY: use_top_fluxes, ideal_solar_division,     &
                 ideal_solar_efolding1, ideal_solar_efolding2

       USE grid_variables,                                                     &
           ONLY:  ddx, ddx2, ddy, ddy2

       USE indices,                                                            &
           ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb,             &
                  nzt, wall_flags_0

       USE kinds

       USE surface_mod,                                                        &
           ONLY :  surf_def_h

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp) ::  zval              !< depth_variable for solar penetration
       REAL(wp) ::  flux1             !< solar flux temp variable
       REAL(wp) ::  flux2             !< solar flux temp variable
       REAL(wp) ::  flag
       REAL(wp) ::  mask_bottom       !< flag to mask vertical upward-facing surface
       REAL(wp) ::  mask_top          !< flag to mask vertical downward-facing surface

       REAL(wp), DIMENSION(1:surf_def_h(2)%ns) ::  s_flux_t           !< flux at model top

#if defined( __nopointer )
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  s  !<
#else
       REAL(wp), DIMENSION(:,:,:), POINTER ::  s  !<
#endif

       REAL(wp), DIMENSION(1:surf_def_h(2)%ns),INTENT(IN),OPTIONAL :: s_flux_solar_t  !<solar flux at sfc

!-- Compute horizontal diffusion

       !$acc data copy( tend ) &
       !$acc copyin( s, s_flux_t, s_flux_solar_t ) &
       !$acc copyout( solar3d ) &
       !$acc present( kh ) &
       !$acc present( ddzu, ddzw, rho_air_zw, drho_air, wall_flags_0 )

       !$acc parallel
       !$acc loop gang vector collapse(3)
       DO  i = nxl, nxr
          DO  j = nys,nyn
             DO  k = nzb+1, nzt
                tend(k,j,i) = tend(k,j,i)                                      &
                                          + 0.5_wp * (                         &
                                     ( kh(k,j,i) + kh(k,j,i+1) )               &
                                   * ( s(k,j,i+1) - s(k,j,i)   )               &
                                   - ( kh(k,j,i) + kh(k,j,i-1) )               &
                                   * ( s(k,j,i)   - s(k,j,i-1) )               &
                                                     ) * ddx2                  &
                                          + 0.5_wp * (                         &
                                     ( kh(k,j,i) + kh(k,j+1,i) )               &
                                   * ( s(k,j+1,i) - s(k,j,i)   )               &
                                   - ( kh(k,j,i) + kh(k,j-1,i) )               &
                                   * ( s(k,j,i)   - s(k,j-1,i) )               &
                                                     ) * ddy2
             ENDDO
          ENDDO
       ENDDO

       !$acc loop gang vector collapse(3)
       DO i = nxl, nxr
          DO j = nys, nyn
!
!--          Compute vertical diffusion. In case that surface fluxes have been
!--          prescribed or computed at bottom and/or top, index k starts/ends at
!--          nzb+2 or nzt-1, respectively. Model top is also mask if top flux
!--          is given.
             DO  k = nzb+1, nzt
!
!--             Determine flags to mask topography below and above. Flag 0 is
!--             used to mask topography in general, and flag 8 implies
!--             information about use_surface_fluxes. Flag 9 is used to control
!--             flux at model top.
                mask_bottom = MERGE( 1.0_wp, 0.0_wp,                           &
                                     BTEST( wall_flags_0(k-1,j,i), 8 ) )
                mask_top    = MERGE( 1.0_wp, 0.0_wp,                           &
                                     BTEST( wall_flags_0(k+1,j,i), 8 ) ) *     &
                              MERGE( 1.0_wp, 0.0_wp,                           &
                                     BTEST( wall_flags_0(k+1,j,i), 9 ) )
                flag        = MERGE( 1.0_wp, 0.0_wp,                           &
                                     BTEST( wall_flags_0(k,j,i), 0 ) )

                tend(k,j,i) = tend(k,j,i)                                      &
                                       + 0.5_wp * (                            &
                                      ( kh(k,j,i) + kh(k+1,j,i) ) *            &
                                          ( s(k+1,j,i)-s(k,j,i) ) * ddzu(k+1)  &
                                                            * rho_air_zw(k)    &
                                                            * mask_top         &
                                    - ( kh(k,j,i) + kh(k-1,j,i) ) *            &
                                          ( s(k,j,i)-s(k-1,j,i) ) * ddzu(k)    &
                                                            * rho_air_zw(k-1)  &
                                                            * mask_bottom      &
                                                  ) * ddzw(k) * drho_air(k)    &
                                                              * flag
             ENDDO
          ENDDO
       ENDDO
       !$acc end parallel

       !$acc parallel
       !$acc loop gang vector collapse(2)
       DO i=nxl,nxr
          DO j=nys,nyn
             !LPV adding solar forcing with depth
             IF ( PRESENT(s_flux_solar_t )) THEN
                m = surf_def_h(2)%start_index(j,i)
                zval = 0.0_wp
                !$acc loop vector
                DO k = nzt,nzb+1,-1
                   flux1 = (1.0_wp - ideal_solar_division)*exp(ideal_solar_efolding2*zval) + &
                             ideal_solar_division*exp(ideal_solar_efolding1*zval)
                   zval = zval - dzw(k)
                   flux2 = (1.0_wp - ideal_solar_division)*exp(ideal_solar_efolding2*zval) + &
                   ideal_solar_division*exp(ideal_solar_efolding1*zval)

                   tend(k,j,i) = tend(k,j,i) - s_flux_solar_t(m)*(flux1 - flux2) / dzw(k)

                   solar3d(k,j,i) = -s_flux_solar_t(m)*(flux1 - flux2) / dzw(k)
                ENDDO
             ENDIF

!
!--          Vertical diffusion at the last computational gridpoint along z-direction
             IF ( use_top_fluxes )  THEN
                surf_s = surf_def_h(2)%start_index(j,i)
                surf_e = surf_def_h(2)%end_index(j,i)
                !$acc loop vector
                DO  m = surf_s, surf_e

                   k   = surf_def_h(2)%k(m)
                   tend(k,j,i) = tend(k,j,i)                                   &
                           + ( - s_flux_t(m) ) * ddzw(k) * drho_air(k)
                ENDDO
             ENDIF
          ENDDO
       ENDDO
       !$acc end parallel
       !$acc end data

    END SUBROUTINE diffusion_s

 END MODULE diffusion_s_mod
