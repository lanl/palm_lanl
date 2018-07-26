!> @file calc_radiation.f90
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
! $Id: calc_radiation.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography concept
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
! 1322 2014-03-20 16:38:49Z raasch
! exponent 4.0 changed to integer
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
! Revision 1.1  2000/04/13 14:42:45  schroeter
! Initial revision
!
!
! Description:
! -------------
!> Calculation of the vertical divergences of the long-wave radiation-fluxes
!> based on the parameterization of the cloud effective emissivity
!------------------------------------------------------------------------------!
 MODULE calc_radiation_mod
 
    USE kinds
    
    PRIVATE
    PUBLIC calc_radiation
    
    LOGICAL, SAVE ::  first_call = .TRUE. !<
    REAL(wp), SAVE ::  sigma = 5.67E-08_wp   !<

    REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::  lwp_ground         !<
    REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::  lwp_top            !<
    REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::  blackbody_emission !<

    INTERFACE calc_radiation
       MODULE PROCEDURE calc_radiation
       MODULE PROCEDURE calc_radiation_ij
    END INTERFACE calc_radiation
 
 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE calc_radiation

       USE arrays_3d,                                                          &
           ONLY:  dzw, pt, ql, tend

       USE cloud_parameters,                                                   &
           ONLY:  cp, l_d_cp, pt_d_t, t_d_pt

       USE control_parameters,                                                 &
           ONLY:  rho_surface

       USE indices,                                                            &
           ONLY:  nxl, nxr, nyn, nys, nzb, nzt, wall_flags_0

       USE kinds

       USE pegrid


       IMPLICIT NONE

       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  k      !<
       INTEGER(iwp) ::  k_help !<
 
       REAL(wp) :: df_p                      !<
       REAL(wp) :: df_m                      !<
       REAL(wp) :: effective_emission_up_m   !<
       REAL(wp) :: effective_emission_up_p   !<
       REAL(wp) :: effective_emission_down_m !<
       REAL(wp) :: effective_emission_down_p !<
       REAL(wp) :: f_up_m                    !<
       REAL(wp) :: f_up_p                    !<
       REAL(wp) :: f_down_m                  !<
       REAL(wp) :: f_down_p                  !<
       REAL(wp) :: impinging_flux_at_top     !<
       REAL(wp) :: temperature               !<


!
!--    On first call, allocate temporary arrays
       IF ( first_call )  THEN
          ALLOCATE( blackbody_emission(nzb:nzt+1), lwp_ground(nzb:nzt+1),      &
                    lwp_top(nzb:nzt+1) )
          first_call = .FALSE.
       ENDIF


       DO  i = nxl, nxr
          DO  j = nys, nyn
!
!--          Compute the liquid water path (LWP) and blackbody_emission
!--          at all vertical levels
             lwp_ground(nzb) = 0.0_wp
             lwp_top(nzt+1)  = rho_surface * ql(nzt+1,j,i) * dzw(nzt+1)

             temperature     = pt(nzb,j,i) * t_d_pt(nzb) + l_d_cp * ql(nzb,j,i)
             blackbody_emission(nzb) = sigma * temperature**4

             DO  k = nzb+1, nzt

                k_help = ( nzt+nzb+1 ) - k
                lwp_ground(k)   = lwp_ground(k-1) + rho_surface * ql(k,j,i) *  &
                                  dzw(k)

                lwp_top(k_help) = lwp_top(k_help+1) +                          &
                                  rho_surface * ql(k_help,j,i) * dzw(k_help)

                temperature     = pt(k,j,i) * t_d_pt(k) + l_d_cp * ql(k,j,i)
                blackbody_emission(k) = sigma * temperature**4                 &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 0 ) )

             ENDDO

             lwp_ground(nzt+1) = lwp_ground(nzt) +                             &
                                 rho_surface * ql(nzt+1,j,i) * dzw(nzt+1)
             lwp_top(nzb)      = lwp_top(nzb+1)

             temperature       = pt(nzt+1,j,i) * t_d_pt(nzt+1) + l_d_cp *      &
                                 ql(nzt+1,j,i)
             blackbody_emission(nzt+1) = sigma * temperature**4

!
!--          See Chlond '92, this is just a first guess
             impinging_flux_at_top = blackbody_emission(nzb) - 100.0_wp

             DO  k = nzb+1, nzt
!
!--             Save some computational time, but this may cause load
!--             imbalances if ql is not distributed uniformly
                IF ( ql(k,j,i) /= 0.0_wp ) THEN
!
!--                Compute effective emissivities
                   effective_emission_up_p   = 1.0_wp -                        &
                                               EXP( -130.0_wp * lwp_ground(k+1) )
                   effective_emission_up_m   = 1.0_wp -                        &
                                               EXP( -130.0_wp * lwp_ground(k-1) )
                   effective_emission_down_p = 1.0_wp -                        &
                                               EXP( -158.0_wp * lwp_top(k+1) )
                   effective_emission_down_m = 1.0_wp -                        &
                                               EXP( -158.0_wp * lwp_top(k-1) )  

!
!--                Compute vertical long wave radiation fluxes
                   f_up_p = blackbody_emission(nzb) +                          &
                            effective_emission_up_p *                          &
                           ( blackbody_emission(k) - blackbody_emission(nzb) )

                   f_up_m = blackbody_emission(nzb) +                          &
                            effective_emission_up_m *                          &
                           ( blackbody_emission(k-1) - blackbody_emission(nzb) )

                   f_down_p = impinging_flux_at_top +                          &
                              effective_emission_down_p *                      &
                             ( blackbody_emission(k) - impinging_flux_at_top )

                   f_down_m = impinging_flux_at_top +                          &
                              effective_emission_down_m *                      &
                             ( blackbody_emission(k-1) - impinging_flux_at_top )

!
!--                Divergence of vertical long wave radiation fluxes
                   df_p = f_up_p - f_down_p
                   df_m = f_up_m - f_down_m

!
!--                Compute tendency term         
                   tend(k,j,i) = tend(k,j,i) -                                 &
                                ( pt_d_t(k) / ( rho_surface * cp ) *           &
                                  ( df_p - df_m ) / dzw(k) )                   &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 0 ) )

                ENDIF

             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE calc_radiation


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE calc_radiation_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  dzw, pt, ql, tend

       USE cloud_parameters,                                                   &
           ONLY:  cp, l_d_cp, pt_d_t, t_d_pt

       USE control_parameters,                                                 &
           ONLY:  rho_surface

       USE indices,                                                            &
           ONLY:  nzb, nzt, wall_flags_0

       USE kinds

       USE pegrid

   
       IMPLICIT NONE

       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  k      !<
       INTEGER(iwp) ::  k_help !<

       REAL(wp) :: df_p                      !<
       REAL(wp) :: df_m                      !<
       REAL(wp) :: effective_emission_up_m   !<
       REAL(wp) :: effective_emission_up_p   !<
       REAL(wp) :: effective_emission_down_m !<
       REAL(wp) :: effective_emission_down_p !<
       REAL(wp) :: f_up_m                    !<
       REAL(wp) :: f_up_p                    !<
       REAL(wp) :: f_down_m                  !<
       REAL(wp) :: f_down_p                  !<
       REAL(wp) :: impinging_flux_at_top     !<
       REAL(wp) :: temperature               !<

       
!
!--    On first call, allocate temporary arrays
       IF ( first_call )  THEN
          ALLOCATE( blackbody_emission(nzb:nzt+1), lwp_ground(nzb:nzt+1),      &
                    lwp_top(nzb:nzt+1) )
          first_call = .FALSE.
       ENDIF

!
!--    Compute the liquid water path (LWP) and blackbody_emission
!--    at all vertical levels
       lwp_ground(nzb) = 0.0_wp
       lwp_top(nzt+1)  = rho_surface * ql(nzt+1,j,i) * dzw(nzt+1)

       temperature     = pt(nzb,j,i) * t_d_pt(nzb) + l_d_cp * ql(nzb,j,i)
       blackbody_emission(nzb) = sigma * temperature**4

       DO  k = nzb+1, nzt
          k_help = ( nzt+nzb+1 ) - k
          lwp_ground(k)   = lwp_ground(k-1) + rho_surface * ql(k,j,i) * dzw(k)

          lwp_top(k_help) = lwp_top(k_help+1) +                                &
                            rho_surface * ql(k_help,j,i) * dzw(k_help)

          temperature     = pt(k,j,i) * t_d_pt(k) + l_d_cp * ql(k,j,i)
          blackbody_emission(k) = sigma * temperature**4                       &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 0 ) )

       ENDDO
       lwp_ground(nzt+1) = lwp_ground(nzt) +                                   &
                           rho_surface * ql(nzt+1,j,i) * dzw(nzt+1)
       lwp_top(nzb)      = lwp_top(nzb+1)

       temperature       = pt(nzt+1,j,i) * t_d_pt(nzt+1) + l_d_cp *            &
                           ql(nzt+1,j,i)
       blackbody_emission(nzt+1) = sigma * temperature**4

!
!--    See Chlond '92, this is just a first guess
       impinging_flux_at_top = blackbody_emission(nzb) - 100.0_wp

       DO  k = nzb+1, nzt
!
!--       Store some computational time,
!--       this may cause load imbalances if ql is not distributed uniformly
          IF ( ql(k,j,i) /= 0.0_wp ) THEN
!
!--          Compute effective emissivities
             effective_emission_up_p   = 1.0_wp -                              &
                                         EXP( -130.0_wp * lwp_ground(k+1) )
             effective_emission_up_m   = 1.0_wp -                              &
                                         EXP( -130.0_wp * lwp_ground(k-1) )
             effective_emission_down_p = 1.0_wp -                              &
                                         EXP( -158.0_wp * lwp_top(k+1) )
             effective_emission_down_m = 1.0_wp -                              &
                                         EXP( -158.0_wp * lwp_top(k-1) )  
             
!
!--          Compute vertical long wave radiation fluxes
             f_up_p = blackbody_emission(nzb) + effective_emission_up_p *      &
                     ( blackbody_emission(k) - blackbody_emission(nzb) )

             f_up_m = blackbody_emission(nzb) + effective_emission_up_m *      &
                     ( blackbody_emission(k-1) - blackbody_emission(nzb) )

             f_down_p = impinging_flux_at_top + effective_emission_down_p *    &
                       ( blackbody_emission(k) - impinging_flux_at_top )

             f_down_m = impinging_flux_at_top + effective_emission_down_m *    &
                       ( blackbody_emission(k-1) - impinging_flux_at_top )

!
!-           Divergence of vertical long wave radiation fluxes
             df_p = f_up_p - f_down_p
             df_m = f_up_m - f_down_m

!
!--          Compute tendency term         
             tend(k,j,i) = tend(k,j,i) - ( pt_d_t(k) / ( rho_surface * cp ) *  &
                                         ( df_p - df_m ) / dzw(k) )            &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 0 ) )

          ENDIF

       ENDDO

    END SUBROUTINE calc_radiation_ij

 END MODULE calc_radiation_mod
