!> @file eqn_state_seawater.f90
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
! $Id: eqn_state_seawater.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2369 2017-08-22 15:20:37Z suehring
! Bugfix, do not mask topography here, since density becomes zero, leading to
! division by zero in production_e
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography and surface concept
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean
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
! 97 2007-06-21 08:23:15Z raasch
! Initial revision
!
!
! Description:
! ------------
!> Equation of state for seawater as a function of potential temperature,
!> salinity, and pressure.
!> For coefficients see Jackett et al., 2006: J. Atm. Ocean Tech.
!> eqn_state_seawater calculates the potential density referred at hyp(0).
!> eqn_state_seawater_func calculates density.
!------------------------------------------------------------------------------!
 MODULE eqn_state_seawater_mod


    USE kinds

    IMPLICIT NONE

    PRIVATE
    PUBLIC eqn_state_seawater, eqn_state_seawater_func

    REAL(wp), DIMENSION(12), PARAMETER ::  nom =                               &
                          (/ 9.9984085444849347D2,   7.3471625860981584D0,     &
                            -5.3211231792841769D-2,  3.6492439109814549D-4,    &
                             2.5880571023991390D0,  -6.7168282786692354D-3,    &
                             1.9203202055760151D-3,  1.1798263740430364D-2,    &
                             9.8920219266399117D-8,  4.6996642771754730D-6,    &
                            -2.5862187075154352D-8, -3.2921414007960662D-12 /)
                          !<

    REAL(wp), DIMENSION(13), PARAMETER ::  den =                               &
                          (/ 1.0D0,                  7.2815210113327091D-3,    &
                            -4.4787265461983921D-5,  3.3851002965802430D-7,    &
                             1.3651202389758572D-10, 1.7632126669040377D-3,    &
                            -8.8066583251206474D-6, -1.8832689434804897D-10,   &
                             5.7463776745432097D-6,  1.4716275472242334D-9,    &
                             6.7103246285651894D-6, -2.4461698007024582D-17,   &
                            -9.1534417604289062D-18 /)
                          !<

    INTERFACE eqn_state_seawater
       MODULE PROCEDURE eqn_state_seawater
    END INTERFACE eqn_state_seawater

    INTERFACE eqn_state_seawater_func
       MODULE PROCEDURE eqn_state_seawater_func
    END INTERFACE eqn_state_seawater_func

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE eqn_state_seawater

       USE arrays_3d,                                                          &
           ONLY:  hyp, prho, pt_p, alpha_T, beta_S, rho_ocean, sa_p
       USE indices,                                                            &
           ONLY:  nxl, nxr, nyn, nys, nzb, nzt

       USE surface_mod,                                                        &
          ONLY :  bc_h

       USE control_parameters,                                                 &
          ONLY :  linear_eqnOfState, alpha_const, beta_const, fixed_alpha,     &
              rho_ref, pt_ref, sa_ref

       IMPLICIT NONE

       INTEGER(iwp) ::  i       !< running index x direction
       INTEGER(iwp) ::  j       !< running index y direction
       INTEGER(iwp) ::  k       !< running index z direction
       INTEGER(iwp) ::  l       !< running index of surface type, south- or north-facing wall
       INTEGER(iwp) ::  m       !< running index surface elements
       INTEGER(iwp) ::  surf_e  !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s  !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp) ::  pden   !<
       REAL(wp) ::  pnom   !<
       REAL(wp) ::  dpdendT1 !<
       REAL(wp) ::  dpdendS1 !<
       REAL(wp) ::  dpnomdT1 !<
       REAL(wp) ::  dpnomdS1 !<
       REAL(wp) ::  dpdendT2 !<
       REAL(wp) ::  dpdendS2 !<
       REAL(wp) ::  dpnomdT2 !<
       REAL(wp) ::  dpnomdS2 !<
       REAL(wp) ::  p1     !<
       REAL(wp) ::  p2     !<
       REAL(wp) ::  p3     !<
       REAL(wp) ::  pt1    !<
       REAL(wp) ::  pt2    !<
       REAL(wp) ::  pt3    !<
       REAL(wp) ::  pt4    !<
       REAL(wp) ::  sa1    !<
       REAL(wp) ::  sa15   !<
       REAL(wp) ::  sa2    !<


       !$acc parallel present( pt_p, sa_p ) &
       !$acc present( hyp, rho_ocean, prho, alpha_T, beta_S )

       !$acc loop collapse(2)
       DO  i = nxl, nxr
          DO  j = nys, nyn
             !$acc loop seq
             DO  k = nzb+1, nzt
             !
!--             Pressure is needed in dbar
                p1 = hyp(k) * 1E-4_wp
                p2 = p1 * p1
                p3 = p2 * p1

!
!--             Temperature needed in degree Celsius
                pt1 = pt_p(k,j,i) - 273.15_wp
                pt2 = pt1 * pt1
                pt3 = pt1 * pt2
                pt4 = pt2 * pt2

                sa1  = sa_p(k,j,i)
                sa15 = sa1 * SQRT( sa1 )
                sa2  = sa1 * sa1

                pnom = nom(1)           + nom(2)*pt1     + nom(3)*pt2     +    &
                       nom(4)*pt3       + nom(5)*sa1     + nom(6)*sa1*pt1 +    &
                       nom(7)*sa2
!-- LPV         Compute pieces for alpha_T and beta_S in steps

                dpnomdT1 = nom(2)       + 2.0*nom(3)*pt1 + 3.0*nom(4)*pt2 +    &
                       nom(6)*sa1
                dpnomdS1 = nom(5)       + nom(6)*pt1     + 2.0*nom(7)*sa1

                pden = den(1)           + den(2)*pt1     + den(3)*pt2     +    &
                       den(4)*pt3       + den(5)*pt4     + den(6)*sa1     +    &
                       den(7)*sa1*pt1   + den(8)*sa1*pt3 + den(9)*sa15    +    &
                       den(10)*sa15*pt2
!
                dpdendT1 = den(2) + 2.0*den(3)*pt1 + 3.0*den(4)*pt2 + 4.0*den(5)*pt3 + &
                       den(7)*sa1 + 3.0*den(8)*sa1*pt2 + 2.0*den(10)*sa15*pt1
                dpdendS1 = den(6) + den(7)*pt1 + den(8)*pt3 + 1.5*den(9)*sqrt(sa1) + &
                       1.5*den(10)*pt2*sqrt(sa1)

!--             Potential density (without pressure terms)
                prho(k,j,i) = pnom / pden

                pnom = pnom +             nom(8)*p1      + nom(9)*p1*pt2  +    &
                       nom(10)*p1*sa1   + nom(11)*p2     + nom(12)*p2*pt2

                pden = pden +             den(11)*p1     + den(12)*p2*pt3 +    &
                       den(13)*p3*pt1

                dpnomdT2 = dpnomdT1     + 2.0*nom(9)*p1*pt1  + 2.0*nom(12)*p2*pt1
                dpnomdS2 = dpnomdS1     + nom(10)*p1

                dpdendT2 = dpdendT1 + 3.0*den(12)*p2*pt2 + den(13)*p3
                dpdendS2 = dpdendS1

!
!--             In-situ density

                rho_ocean(k,j,i) = pnom / pden

                alpha_T(k,j,i) = -1.0/rho_ocean(k,j,i)*((dpnomdT2*pden -         &
                                    pnom*dpdendT2) / (pden*pden))
                beta_S(k,j,i) = 1.0/rho_ocean(k,j,i)*(dpnomdS2*pden - pnom*dpdendS2) / (pden*pden)
                rho_ocean(k,j,i) = pnom / pden

                if (linear_eqnOfState) THEN
                  if (fixed_alpha) THEN
                    rho_ocean(k,j,i) = rho_ref*(1.0 - alpha_const*(pt1 - pt_ref) + &
                        beta_const*(sa1 - sa_ref))
                  ELSE
                    rho_ocean(k,j,i) = rho_ref*(1.0 - alpha_T(k,j,i)*(pt1 - pt_ref) + &
                        beta_S(k,j,i)*(sa1 - sa_ref))
                  ENDIF
                ENDIF
             ENDDO
!
!--          Neumann conditions are assumed at top boundary
             prho(nzt+1,j,i)      = prho(nzt,j,i)
             rho_ocean(nzt+1,j,i) = rho_ocean(nzt,j,i)

          ENDDO
       ENDDO
       !$acc end parallel
!
!--    Neumann conditions at up/downward-facing surfaces
       !$OMP PARALLEL DO PRIVATE( i, j, k )
       !$acc parallel present( bc_h, prho, rho_ocean )
       !$acc loop
       DO  m = 1, bc_h(0)%ns
          i = bc_h(0)%i(m)
          j = bc_h(0)%j(m)
          k = bc_h(0)%k(m)
          prho(k-1,j,i)      = prho(k,j,i)
          rho_ocean(k-1,j,i) = rho_ocean(k,j,i)
       ENDDO
!
!--    Downward facing surfaces
       !$OMP PARALLEL DO PRIVATE( i, j, k )
       !$acc loop
       DO  m = 1, bc_h(1)%ns
          i = bc_h(1)%i(m)
          j = bc_h(1)%j(m)
          k = bc_h(1)%k(m)
          prho(k+1,j,i)      = prho(k,j,i)
          rho_ocean(k+1,j,i) = rho_ocean(k,j,i)
       ENDDO
       !$acc end parallel

    END SUBROUTINE eqn_state_seawater


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Equation of state as a function
!------------------------------------------------------------------------------!
    REAL(wp) FUNCTION eqn_state_seawater_func( p, pt, sa )

       IMPLICIT NONE

       REAL(wp) ::  p      !<
       REAL(wp) ::  p1     !<
       REAL(wp) ::  p2     !<
       REAL(wp) ::  p3     !<
       REAL(wp) ::  pt     !<
       REAL(wp) ::  pt1    !<
       REAL(wp) ::  pt2    !<
       REAL(wp) ::  pt3    !<
       REAL(wp) ::  pt4    !<
       REAL(wp) ::  sa     !<
       REAL(wp) ::  sa15   !<
       REAL(wp) ::  sa2    !<

!
!--    Pressure is needed in dbar
       p1 = p  * 1E-4_wp
       p2 = p1 * p1
       p3 = p2 * p1

!
!--    Temperature needed in degree Celsius
       pt1 = pt - 273.15_wp
       pt2 = pt1 * pt1
       pt3 = pt1 * pt2
       pt4 = pt2 * pt2

       sa15 = sa * SQRT( sa )
       sa2  = sa * sa

       eqn_state_seawater_func =                                               &
         ( nom(1)        + nom(2)*pt1       + nom(3)*pt2    + nom(4)*pt3     + &
           nom(5)*sa     + nom(6)*sa*pt1    + nom(7)*sa2    + nom(8)*p1      + &
           nom(9)*p1*pt2 + nom(10)*p1*sa    + nom(11)*p2    + nom(12)*p2*pt2   &
         ) /                                                                   &
         ( den(1)        + den(2)*pt1       + den(3)*pt2    + den(4)*pt3     + &
           den(5)*pt4    + den(6)*sa        + den(7)*sa*pt1 + den(8)*sa*pt3  + &
           den(9)*sa15   + den(10)*sa15*pt2 + den(11)*p1    + den(12)*p2*pt3 + &
           den(13)*p3*pt1                                                      &
         )


    END FUNCTION eqn_state_seawater_func

 END MODULE eqn_state_seawater_mod
