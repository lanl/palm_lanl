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
! 2018-11-21 cbegeman
! Added functions for freezing point and derivative of freezing point with 
! respect to salinity
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
    PUBLIC eqn_state_seawater, eqn_state_seawater_func,                        &
           pt_freezing, pt_freezing_SA, T_freezing, T_freezing_SA

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

    REAL(wp), DIMENSION(11), PARAMETER ::  ptfrnom =                           &
                          (/ 2.5180516744541290e-3, -5.8545863698926184e-2,    &
                             2.2979985780124325e-3, -3.0086338218235500e-4,    &
                            -7.0023530029351803e-4,  8.4149607219833806e-9,    &
                             1.1845857563107403e-11, 1.                   ,    &
                             1.3632481944285909e-6, -3.8493266309172074e-5,    &
                             9.1686537446749641e-10 /)
    
    REAL(wp), DIMENSION(23), PARAMETER ::  tfrnom =                            &
                          (/ 0.002519             , -5.946302841607319    ,    &
                             4.136051661346983    , -1.115150523403847e1  ,    &
                             1.476878746184548e1  , -1.088873263630961e1  ,    &
                             2.961018839640730    , -7.433320943962606    ,    &
                             -1.561578562479883   ,  4.073774363480365e-2 ,    &
                             1.158414435887717e-2 , -4.122639292422863e-1 ,    &
                             -1.123186915628260e-1,  5.715012685553502e-1 ,    &
                             2.021682115652684e-1 ,  4.140574258089767e-2 ,    &
                             -6.034228641903586e-1, -1.205825928146808e-2 ,    &
                             -2.812172968619369e-1,  1.877244474023750e-2 ,    &
                             -1.204395563789007e-1,  2.349147739749606e-1 ,    &
                             2.748444541144219e-3  /)
       
    INTERFACE eqn_state_seawater
       MODULE PROCEDURE eqn_state_seawater
       MODULE PROCEDURE eqn_state_seawater_ij
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
          ONLY :  drho_dp_const, linear_eqnOfState, alpha_const, beta_const, fixed_alpha,     &
                  message_string,surface_pressure, rho_ref, pt_ref, sa_ref

       IMPLICIT NONE

       INTEGER(iwp) ::  i       !< running index x direction
       INTEGER(iwp) ::  j       !< running index y direction
       INTEGER(iwp) ::  k       !< running index z direction
       INTEGER(iwp) ::  l       !< running index of surface type, south- or north-facing wall
       INTEGER(iwp) ::  m       !< running index surface elements
       INTEGER(iwp) ::  surf_e  !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s  !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp) ::  pden   !<
       REAL(wp) ::  pden_surface = 0.0_wp   !<
       REAL(wp) ::  pnom   !<
       REAL(wp) ::  pnom_surface = 0.0_wp   !<
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
       REAL(wp) ::  p1_surface     !<
       REAL(wp) ::  p2_surface     !<
       REAL(wp) ::  p3_surface     !<
       REAL(wp) ::  pt1    !<
       REAL(wp) ::  pt2    !<
       REAL(wp) ::  pt3    !<
       REAL(wp) ::  pt4    !<
       REAL(wp) ::  sa1    !<
       REAL(wp) ::  sa15   !<
       REAL(wp) ::  sa2    !<


       IF ( surface_pressure > 1014.0_wp ) THEN
          p1_surface = surface_pressure * 1E-2_wp
          p2_surface = p1_surface * p1_surface
          p3_surface = p2_surface * p1_surface
       ENDIF

       IF (linear_eqnOfState .AND. fixed_alpha) THEN
          DO  k = nzb+1, nzt
             DO  j = nys, nyn
                DO  i = nxl, nxr
                   prho(k,j,i)      = rho_ref*(1.0 -                           &
                                         alpha_const*(pt_p(k,j,i) - pt_ref) +  &
                                         beta_const*(sa_p(k,j,i) - sa_ref) )
                   rho_ocean(k,j,i) = rho_ref*(1.0 -                           &
                                         alpha_const*(pt_p(k,j,i) - pt_ref) +  &
                                         beta_const*(sa_p(k,j,i) - sa_ref) +   & 
                                         drho_dp_const * hyp(k) * 1E-4_wp     )
                ENDDO
             ENDDO
          ENDDO
       ELSE
          DO  k = nzb+1, nzt
!
!--          Pressure is needed in dbar
             p1 = hyp(k) * 1E-4_wp
             p2 = p1 * p1
             p3 = p2 * p1

             DO  j = nys, nyn
                DO  i = nxl, nxr
!
!--                Temperature needed in degree Celsius
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
!-- LPV            Compute pieces for alpha_T and beta_S in steps

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

                   IF ( surface_pressure > 1014.0_wp ) THEN
                      pnom_surface =                     nom(8)*p1_surface      + &
                                 nom(9)*p1_surface*pt2 + nom(10)*p1_surface*sa1 + &
                                 nom(11)*p2_surface    + nom(12)*p2_surface*pt2

                      pden_surface =                  den(11)*p1_surface      +   &
                             den(12)*p2_surface*pt3 + den(13)*p3_surface*pt1
                   ENDIF

!--                Potential density referenced to surface_pressure
                   prho(k,j,i) = (pnom + pnom_surface) / (pden + pden_surface)

                   pnom = pnom +             nom(8)*p1      + nom(9)*p1*pt2  +    &
                          nom(10)*p1*sa1   + nom(11)*p2     + nom(12)*p2*pt2

                   pden = pden +             den(11)*p1     + den(12)*p2*pt3 +    &
                          den(13)*p3*pt1

                   dpnomdT2 = dpnomdT1     + 2.0*nom(9)*p1*pt1  + 2.0*nom(12)*p2*pt1
                   dpnomdS2 = dpnomdS1     + nom(10)*p1

                   dpdendT2 = dpdendT1 + 3.0*den(12)*p2*pt2 + den(13)*p3
                   dpdendS2 = dpdendS1

!
!--                In-situ density

                   rho_ocean(k,j,i) = pnom / pden

                   alpha_T(k,j,i) = -1.0/rho_ocean(k,j,i)*((dpnomdT2*pden -         &
                                       pnom*dpdendT2) / (pden*pden))
                   beta_S(k,j,i) = 1.0/rho_ocean(k,j,i)*(dpnomdS2*pden - pnom*dpdendS2) / (pden*pden)

                   IF (linear_eqnOfState) THEN
                      prho(k,j,i)      = rho_ref*(1.0 -                        &
                                         alpha_T(k,j,i)*(pt1 - pt_ref) +       &
                                         beta_S(k,j,i)*(sa1 - sa_ref) )
                      rho_ocean(k,j,i) = rho_ref*(1.0 -                        &
                                         alpha_T(k,j,i)*(pt1 - pt_ref) +       &
                                         beta_S(k,j,i)*(sa1 - sa_ref) +        &
                                         drho_dp_const * hyp(k) * 1E-4_wp)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDIF
!
!--    Neumann conditions are assumed at top boundary
       DO  j = nys, nyn
          DO  i = nxl, nxr
             prho(nzt+1,j,i)      = prho(nzt,j,i)
             rho_ocean(nzt+1,j,i) = rho_ocean(nzt,j,i)
          ENDDO
       ENDDO
!
!--    Neumann conditions at up/downward-facing surfaces
       !$OMP PARALLEL DO PRIVATE( i, j, k )
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
       DO  m = 1, bc_h(1)%ns
          i = bc_h(1)%i(m)
          j = bc_h(1)%j(m)
          k = bc_h(1)%k(m)
          prho(k+1,j,i)      = prho(k,j,i)
          rho_ocean(k+1,j,i) = rho_ocean(k,j,i)
       ENDDO

    END SUBROUTINE eqn_state_seawater


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE eqn_state_seawater_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  hyp, prho, pt_p, alpha_T, beta_S, rho_ocean, sa_p

       USE control_parameters,                                                 &
          ONLY :  drho_dp_const, linear_eqnOfState, alpha_const, beta_const, fixed_alpha,     &
                  message_string,surface_pressure, rho_ref, pt_ref, sa_ref
       
       USE indices,                                                            &
           ONLY:  nzb, nzt

       USE surface_mod,                                                        &
          ONLY :  bc_h

       IMPLICIT NONE

       INTEGER(iwp) ::  i       !< running index x direction
       INTEGER(iwp) ::  j       !< running index y direction
       INTEGER(iwp) ::  k       !< running index z direction
       INTEGER(iwp) ::  l       !< running index of surface type, south- or north-facing wall
       INTEGER(iwp) ::  m       !< running index surface elements
       INTEGER(iwp) ::  surf_e  !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s  !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp) ::  pden   !<
       REAL(wp) ::  pden_surface = 0.0_wp   !<
       REAL(wp) ::  pnom   !<
       REAL(wp) ::  pnom_surface = 0.0_wp   !<
       REAL(wp) ::  p1     !<
       REAL(wp) ::  p2     !<
       REAL(wp) ::  p3     !<
       REAL(wp) ::  p1_surface     !<
       REAL(wp) ::  p2_surface     !<
       REAL(wp) ::  p3_surface     !<
       REAL(wp) ::  pt1    !<
       REAL(wp) ::  pt2    !<
       REAL(wp) ::  pt3    !<
       REAL(wp) ::  pt4    !<
       REAL(wp) ::  sa1    !<
       REAL(wp) ::  sa15   !<
       REAL(wp) ::  sa2    !<

       IF ( surface_pressure > 1014.0_wp ) THEN
          p1_surface = surface_pressure * 1E-2_wp
          p2_surface = p1_surface * p1_surface
          p3_surface = p2_surface * p1_surface
       ENDIF

       IF (linear_eqnOfState .AND. fixed_alpha) THEN
          DO  k = nzb+1, nzt
             prho(k,j,i)      = rho_ref*(1.0 -                                 &
                                   alpha_const*(pt_p(k,j,i) - pt_ref) +        &
                                   beta_const*(sa_p(k,j,i) - sa_ref)  )      
             rho_ocean(k,j,i) = rho_ref*(1.0 -                                 &
                                   alpha_const*(pt_p(k,j,i) - pt_ref) +        &
                                   beta_const*(sa_p(k,j,i) - sa_ref) +         & 
                                   drho_dp_const * hyp(k) * 1E-4_wp     )
          ENDDO
       ELSE
          DO  k = nzb+1, nzt
!
!--          Pressure is needed in dbar
             p1 = hyp(k) * 1E-4_wp
             p2 = p1 * p1
             p3 = p2 * p1

!
!--          Temperature needed in degree Celsius
             pt1 = pt_p(k,j,i) - 273.15_wp
             pt2 = pt1 * pt1
             pt3 = pt1 * pt2
             pt4 = pt2 * pt2

             sa1  = sa_p(k,j,i)
             sa15 = sa1 * SQRT( sa1 )
             sa2  = sa1 * sa1

             pnom = nom(1)           + nom(2)*pt1     + nom(3)*pt2     +          &
                    nom(4)*pt3       + nom(5)*sa1     + nom(6)*sa1*pt1 +          &
                    nom(7)*sa2

             pden = den(1)           + den(2)*pt1     + den(3)*pt2     +          &
                    den(4)*pt3       + den(5)*pt4     + den(6)*sa1     +          &
                    den(7)*sa1*pt1   + den(8)*sa1*pt3 + den(9)*sa15    +          &
                    den(10)*sa15*pt2
             
             IF ( surface_pressure > 1014.0_wp ) THEN
                pnom_surface =                     nom(8)*p1_surface      + &
                           nom(9)*p1_surface*pt2 + nom(10)*p1_surface*sa1 + &
                           nom(11)*p2_surface    + nom(12)*p2_surface*pt2

                pden_surface =                  den(11)*p1_surface      +   &
                       den(12)*p2_surface*pt3 + den(13)*p3_surface*pt1
             ENDIF

!--          Potential density referenced to surface_pressure
             prho(k,j,i) = (pnom + pnom_surface) / (pden + pden_surface)

             pnom = pnom +             nom(8)*p1      + nom(9)*p1*pt2  +          &
                    nom(10)*p1*sa1   + nom(11)*p2     + nom(12)*p2*pt2
             pden = pden +             den(11)*p1     + den(12)*p2*pt3 +          &
                    den(13)*p3*pt1

!
!--          In-situ density
             rho_ocean(k,j,i) = pnom / pden

             IF (linear_eqnOfState) THEN
                prho(k,j,i)      = rho_ref*(1.0 - alpha_T(k,j,i)*(pt1 - pt_ref) + &
                                            beta_S(k,j,i)*(sa1 - sa_ref) )
                rho_ocean(k,j,i) = rho_ref*(1.0 - alpha_T(k,j,i)*(pt1 - pt_ref) + &
                                            beta_S(k,j,i)*(sa1 - sa_ref) +        &
                                            drho_dp_const * hyp(k) * 1E-4_wp)
             ENDIF

          ENDDO
       ENDIF
!
!--    Neumann conditions at up/downward-facing walls
       surf_s = bc_h(0)%start_index(j,i)
       surf_e = bc_h(0)%end_index(j,i)
       DO  m = surf_s, surf_e
          k                  = bc_h(0)%k(m)
          prho(k-1,j,i)      = prho(k,j,i)
          rho_ocean(k-1,j,i) = rho_ocean(k,j,i)
       ENDDO
!
!--    Downward facing surfaces
       surf_s = bc_h(1)%start_index(j,i)
       surf_e = bc_h(1)%end_index(j,i)
       DO  m = surf_s, surf_e
          k                  = bc_h(1)%k(m)
          prho(k+1,j,i)      = prho(k,j,i)
          rho_ocean(k+1,j,i) = rho_ocean(k,j,i)
       ENDDO
!
!--    Neumann condition are assumed at top boundary
       prho(nzt+1,j,i)      = prho(nzt,j,i)
       rho_ocean(nzt+1,j,i) = rho_ocean(nzt,j,i)

    END SUBROUTINE eqn_state_seawater_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Equation of state as a function returning in situ density
!> To return potential density, set p=0
!------------------------------------------------------------------------------!
    REAL(wp) FUNCTION eqn_state_seawater_func( p, pt, sa )
       
       USE control_parameters,                                                 &
          ONLY :  alpha_const, beta_const, drho_dp_const, fixed_alpha,         &
                  linear_eqnOfState, rho_ref, pt_ref, sa_ref

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

       IF (linear_eqnOfState .AND. fixed_alpha) THEN
          eqn_state_seawater_func = rho_ref*(1.0 -                             &
                                             alpha_const*(pt - pt_ref) +       &
                                             beta_const*(sa - sa_ref) +        & 
                                             drho_dp_const * p * 1E-4_wp     )
       ELSE
!
!--       Pressure is needed in dbar
          p1 = p  * 1E-4_wp
          p2 = p1 * p1
          p3 = p2 * p1

!
!--       Temperature needed in degree Celsius
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
       ENDIF

    END FUNCTION eqn_state_seawater_func

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate derivative of in situ temperature freezing point with respect
!> to absolute salinity at a given pressure and absolute salinity according to 
!> TEOS10 polynomial function.
!------------------------------------------------------------------------------!
    REAL (wp) FUNCTION T_freezing_SA( p, SA )
    
       IMPLICIT NONE

       REAL(wp) :: SA, p, p_r, SA_r, x

       T_freezing_SA =                                                         &
           (  tfrnom(2) + x * ( 1.5 * tfrnom(3) +                              &
                                x * ( 2.0 * tfrnom(4) +                        &
                                      x * ( 2.5 * tfrnom(5) +                  &
                                            x * ( 3.0 * tfrnom(6) +            &
                                                  3.5 * tfrnom(7) * x )        &
                                          )                                    &
                                    )                                          &
                              )                                                &
            + p_r * (  tfrnom(11) + x * ( 1.5 * tfrnom(12) +                   &
                                          x * ( 2.0 * tfrnom(14) +             &
                                                x * ( 2.5 * tfrnom(17) +       &
                                                      x * ( 3.0 * tfrnom(20) + &
                                                        3.5 * tfrnom(23) * x ) &
                                                    )                          &
                                              )                                &
                                        )                                      &
                      + p_r * ( tfrnom(13) + x * ( 1.5 * tfrnom(15) +          &
                                                   x * ( 2.0 * tfrnom(18) +    &
                                                         2.5 * tfrnom(21) * x )&
                                                 )                             &
                                + p_r * ( tfrnom(16) + x * ( 1.5 * tfrnom(19) +&
                                                             2.0 * tfrnom(22) * x )&
                                        )                                      &
                              )                                                &
                    )                                                          &
           ) * 1e-2_wp

    END FUNCTION T_freezing_SA

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate potential temperature freezing point at a given pressure and 
!> absolute salinity according to Jackett et al. (2006).
!------------------------------------------------------------------------------!
    REAL(wp) FUNCTION pt_freezing( p, SA )

       IMPLICIT NONE

       REAL(wp) ::  p       !< given in dbar
       REAL(wp) ::  SA
       
       pt_freezing  = ( ( ptfrnom(1)        + ptfrnom(2)*SA   +                    &
                          ptfrnom(3)*SA**1.5 +                                     &
                          ptfrnom(4)*SA**2. + ptfrnom(5)*p    +                    & 
                          ptfrnom(6)*p**2   + ptfrnom(7)*SA*p**2. ) /              &
                        ( ptfrnom(8)        + ptfrnom(9)*SA**2.5 +                 &
                          ptfrnom(10)*p + ptfrnom(11)*p**2.           )            &
                      )

       pt_freezing = pt_freezing + 273.15

    END FUNCTION pt_freezing

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate potential temperature freezing point at a given pressure and 
!> absolute salinity according to Jackett et al. (2006).
!------------------------------------------------------------------------------!
    REAL(wp) FUNCTION pt_freezing_SA( p, SA )

       IMPLICIT NONE

       REAL(wp) ::  p       !< given in dbar
       REAL(wp) ::  SA
       REAL(wp) ::  ptnum,ptden,dnum_dSA,dden_dSA

       ptnum = ( ptfrnom(1)        + ptfrnom(2)*SA + ptfrnom(3)*SA**1.5 +        &
                 ptfrnom(4)*SA**2. + ptfrnom(5)*p +                              &
                 ptfrnom(6)*p**2   + ptfrnom(7)*SA*p**2. )
       ptden = ( ptfrnom(8)        + ptfrnom(9)*SA**2.5 +                        &
                 ptfrnom(10)*p     + ptfrnom(11)*p**2.   )    
       dnum_dSA = ptfrnom(2)     + 1.5*ptfrnom(3)*SA**0.5 +                    &
                  ptfrnom(4)*SA  + ptfrnom(7)*p**2.
       dden_dSA = 2.5*ptfrnom(9)*SA**1.5

       pt_freezing_SA  = ( ( ptden*dnum_dSA - ptnum*dden_dSA ) / ptden**2. )

    END FUNCTION pt_freezing_SA

!------------------------------------------------------------------------------!
! Description:
! ------------
!
!  Calculates the in-situ temperature at which seawater freezes from a 
!  computationally efficient polynomial following TEOS10.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!
!  t_freezing = in-situ temperature at which seawater freezes.    [ deg C ]
!               (ITS-90)                
!------------------------------------------------------------------------------!
    REAL(wp) FUNCTION T_freezing( p, SA )

       IMPLICIT NONE

       REAL(wp) :: SA, p, p_r, SA_r, x

       SA_r = SA * 1e-2_wp
       x = SQRT(SA_r)
       p_r = p * 1e-4_wp

       T_freezing = tfrnom(1) +                                                &
                    SA_r * (   tfrnom(2) + x*(tfrnom(3) + x*(tfrnom(4) +       &
                              x*(tfrnom(5) + x*(tfrnom(6) + tfrnom(7)*x))))) + &
                    p_r  * (tfrnom(8) + p_r*(tfrnom(9) + tfrnom(10)*p_r))    + &
                    SA_r * p_r * (tfrnom(11) +                                 &
                                  p_r*(tfrnom(13) +                            &
                                       p_r*(tfrnom(16) + tfrnom(22)*sa_r)) +   &
                                  sa_r * (tfrnom(14) + tfrnom(18)*p_r +        &
                                          tfrnom(20)*sa_r) +                   &
                                  x    * (tfrnom(12) +                         &
                                          p_r*(tfrnom(15) + tfrnom(19)*p_r) +  &
                                          sa_r * (tfrnom(17) +                 &
                                                  tfrnom(21)*p_r +             &
                                                  tfrnom(23)*sa_r)             &
                                         )                                     &
                                 )
       T_freezing = T_freezing + 273.15
   
   END FUNCTION T_freezing

 END MODULE eqn_state_seawater_mod
