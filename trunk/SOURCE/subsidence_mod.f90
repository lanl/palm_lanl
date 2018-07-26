!> @file subsidence_mod.f90
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
! $Id: subsidence_mod.f90 3045 2018-05-28 07:55:41Z Giersch $
! Error message revised
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
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1862 2016-04-14 09:07:42Z hoffmann
! Bugfix: In case of vector machine optimized model runs, sums_ls_l should not 
! be addressed if large_scale_forcing is false because sums_ls_l is not 
! allocated. 
!
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
! 
! 1729 2015-11-20 11:01:48Z gronemeier
! Bugfix: shifting of initial profiles only once each time step
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
!
! 1489 2014-10-30 08:09:12Z raasch
! bugfix: sums_ls_l can only be used if large_scale_forcing is switched on
!
! 1382 2014-04-30 12:15:41Z boeske
! Changed the weighting factor that is used in the summation of subsidence 
! tendencies for profile data output from weight_pres to weight_substep
! added Neumann boundary conditions for profile data output of subsidence terms
! at nzt+1
! 
! 1380 2014-04-28 12:40:45Z heinze
! Shifting only necessary in case of scalar Rayleigh damping
! 
! 1365 2014-04-22 15:03:56Z boeske
! Summation of subsidence tendencies for data output added 
! +ls_index, sums_ls_l, tmp_tend
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
! Revision 3.7 2009-12-11 14:15:58Z heinze
! Initial revision 
!
! Description:
! ------------
!> Impact of large-scale subsidence or ascent as tendency term for use 
!> in the prognostic equation of potential temperature. This enables the 
!> construction of a constant boundary layer height z_i with time.
!-----------------------------------------------------------------------------!
 MODULE subsidence_mod
 


    IMPLICIT NONE

    PRIVATE
    PUBLIC  init_w_subsidence, subsidence

    INTERFACE init_w_subsidence
       MODULE PROCEDURE init_w_subsidence
    END INTERFACE init_w_subsidence

    INTERFACE subsidence
       MODULE PROCEDURE subsidence
       MODULE PROCEDURE subsidence_ij
    END INTERFACE subsidence

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE init_w_subsidence 

       USE arrays_3d,                                                          &
           ONLY:  dzu, w_subs, zu

       USE control_parameters,                                                 &
           ONLY:  message_string, ocean, subs_vertical_gradient,               &
                  subs_vertical_gradient_level, subs_vertical_gradient_level_i

       USE indices,                                                            &
           ONLY:  nzb, nzt

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i !< 
       INTEGER(iwp) ::  k !< 

       REAL(wp)     ::  gradient   !< 
       REAL(wp)     ::  ws_surface !< 

       IF ( .NOT. ALLOCATED( w_subs ) )  THEN
          ALLOCATE( w_subs(nzb:nzt+1) )
          w_subs = 0.0_wp
       ENDIF 

       IF ( ocean )  THEN
          message_string = 'Applying large scale vertical motion is not ' //   &
                           'allowed for ocean runs'
          CALL message( 'init_w_subsidence', 'PA0324', 2, 2, 0, 6, 0 )
       ENDIF

!
!--   Compute the profile of the subsidence/ascent velocity 
!--   using the given gradients
      i = 1
      gradient = 0.0_wp
      ws_surface = 0.0_wp
      

      subs_vertical_gradient_level_i(1) = 0
      DO  k = 1, nzt+1
         IF ( i < 11 )  THEN
            IF ( subs_vertical_gradient_level(i) < zu(k)  .AND. &
                 subs_vertical_gradient_level(i) >= 0.0_wp )  THEN
               gradient = subs_vertical_gradient(i) / 100.0_wp
               subs_vertical_gradient_level_i(i) = k - 1
               i = i + 1
            ENDIF
         ENDIF
         IF ( gradient /= 0.0_wp )  THEN
            IF ( k /= 1 )  THEN
               w_subs(k) = w_subs(k-1) + dzu(k) * gradient
            ELSE
               w_subs(k) = ws_surface   + 0.5_wp * dzu(k) * gradient
            ENDIF
         ELSE
            w_subs(k) = w_subs(k-1)
         ENDIF
      ENDDO

!
!--   In case of no given gradients for the subsidence/ascent velocity,
!--   choose zero gradient
      IF ( subs_vertical_gradient_level(1) == -9999999.9_wp )  THEN
         subs_vertical_gradient_level(1) = 0.0_wp
      ENDIF

    END SUBROUTINE init_w_subsidence


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE subsidence( tendency, var, var_init, ls_index ) 

       USE arrays_3d,                                                          &
           ONLY:  ddzu, w_subs

       USE control_parameters,                                                 &
           ONLY:  dt_3d, intermediate_timestep_count, large_scale_forcing,     &
                  scalar_rayleigh_damping

       USE indices,                                                            &
           ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt,        &
                  wall_flags_0

       USE kinds

       USE statistics,                                                         &
           ONLY:  sums_ls_l, weight_substep

       IMPLICIT NONE
 
       INTEGER(iwp) ::  i !< 
       INTEGER(iwp) ::  j !< 
       INTEGER(iwp) ::  k !< 
       INTEGER(iwp) ::  ls_index !<

       REAL(wp)     ::  tmp_tend !<
       REAL(wp)     ::  tmp_grad !< 
    
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  var      !< 
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  tendency !< 
       REAL(wp), DIMENSION(nzb:nzt+1) ::  var_init                     !< 
       REAL(wp), DIMENSION(nzb:nzt+1) ::  var_mod                      !< 

       var_mod = var_init

!
!--    Influence of w_subsidence on the current tendency term
       DO  i = nxl, nxr
          DO  j = nys, nyn

             DO  k = nzb+1, nzt 
                IF ( w_subs(k) < 0.0_wp )  THEN    ! large-scale subsidence
                   tmp_tend = - w_subs(k) *                                    &
                              ( var(k+1,j,i) - var(k,j,i) ) * ddzu(k+1) *      &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 0 ) )
                ELSE                               ! large-scale ascent
                   tmp_tend = - w_subs(k) *                                    &
                              ( var(k,j,i) - var(k-1,j,i) ) * ddzu(k) *        &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 0 ) )
                ENDIF

                tendency(k,j,i) = tendency(k,j,i) + tmp_tend

                IF ( large_scale_forcing )  THEN
                   sums_ls_l(k,ls_index) = sums_ls_l(k,ls_index) + tmp_tend    &
                                 * weight_substep(intermediate_timestep_count) &
                                 * MERGE( 1.0_wp, 0.0_wp,                      &
                                          BTEST( wall_flags_0(k,j,i), 0 ) )
                ENDIF
             ENDDO

             IF ( large_scale_forcing )  THEN
                sums_ls_l(nzt+1,ls_index) = sums_ls_l(nzt,ls_index)
             ENDIF

          ENDDO
       ENDDO

!
!--    Shifting of the initial profile is especially necessary with Rayleigh 
!--    damping switched on
       IF ( scalar_rayleigh_damping .AND.                                      &
            intermediate_timestep_count == 1 )  THEN
          DO  k = nzb, nzt
             IF ( w_subs(k) < 0.0_wp )  THEN      ! large-scale subsidence
                var_mod(k) = var_init(k) - dt_3d * w_subs(k) *  &
                                  ( var_init(k+1) - var_init(k) ) * ddzu(k+1)
             ENDIF
          ENDDO
!
!--      At the upper boundary, the initial profile is shifted with aid of 
!--      the gradient tmp_grad. (This is ok if the gradients are linear.)
         IF ( w_subs(nzt) < 0.0_wp )  THEN
            tmp_grad = ( var_init(nzt+1) - var_init(nzt) ) * ddzu(nzt+1)
            var_mod(nzt+1) = var_init(nzt+1) -  &
                                 dt_3d * w_subs(nzt+1) * tmp_grad
         ENDIF
       

         DO  k = nzt+1, nzb+1, -1
            IF ( w_subs(k) >= 0.0_wp )  THEN  ! large-scale ascent
               var_mod(k) = var_init(k) - dt_3d * w_subs(k) *  &
                                  ( var_init(k) - var_init(k-1) ) * ddzu(k) 
            ENDIF
         ENDDO
!
!--      At the lower boundary shifting is not necessary because the 
!--      subsidence velocity w_subs(nzb) vanishes.
         IF ( w_subs(nzb+1) >= 0.0_wp )  THEN
            var_mod(nzb) = var_init(nzb)
         ENDIF

         var_init = var_mod
      ENDIF


 END SUBROUTINE subsidence

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
 SUBROUTINE subsidence_ij( i, j, tendency, var, var_init, ls_index ) 

       USE arrays_3d,                                                          &
           ONLY:  ddzu, w_subs

       USE control_parameters,                                                 &
           ONLY:  dt_3d, intermediate_timestep_count, large_scale_forcing,     &
                  scalar_rayleigh_damping

       USE indices,                                                            &
           ONLY:  nxl, nxlg, nxrg, nyng, nys, nysg, nzb, nzt, wall_flags_0

       USE kinds

       USE statistics,                                                         &
           ONLY:  sums_ls_l, weight_substep

       IMPLICIT NONE
 
       INTEGER(iwp) ::  i !< 
       INTEGER(iwp) ::  j !< 
       INTEGER(iwp) ::  k !< 
       INTEGER(iwp) ::  ls_index !<

       REAL(wp)     ::  tmp_tend !<
       REAL(wp)     ::  tmp_grad !< 
    
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  var      !< 
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  tendency !< 
       REAL(wp), DIMENSION(nzb:nzt+1) ::  var_init                     !< 
       REAL(wp), DIMENSION(nzb:nzt+1) ::  var_mod                      !< 

       var_mod = var_init

!
!--    Influence of w_subsidence on the current tendency term
       DO  k = nzb+1, nzt 
          IF ( w_subs(k) < 0.0_wp )  THEN      ! large-scale subsidence
             tmp_tend = - w_subs(k) * ( var(k+1,j,i) - var(k,j,i) )            &
                                    * ddzu(k+1)                                &
                                    * MERGE( 1.0_wp, 0.0_wp,                   &
                                             BTEST( wall_flags_0(k,j,i), 0 ) )
          ELSE                                 ! large-scale ascent
             tmp_tend = - w_subs(k) * ( var(k,j,i) - var(k-1,j,i) ) * ddzu(k)  &
                                    * MERGE( 1.0_wp, 0.0_wp,                   &
                                             BTEST( wall_flags_0(k,j,i), 0 ) )
          ENDIF

          tendency(k,j,i) = tendency(k,j,i) + tmp_tend

          IF ( large_scale_forcing )  THEN
             sums_ls_l(k,ls_index) = sums_ls_l(k,ls_index) + tmp_tend          &
                                  * weight_substep(intermediate_timestep_count)&
                                  * MERGE( 1.0_wp, 0.0_wp,                     &
                                           BTEST( wall_flags_0(k,j,i), 0 ) )
          ENDIF
       ENDDO

       IF ( large_scale_forcing )  THEN
          sums_ls_l(nzt+1,ls_index) = sums_ls_l(nzt,ls_index)
       ENDIF

!
!--    Shifting of the initial profile is especially necessary with Rayleigh 
!--    damping switched on
       IF ( scalar_rayleigh_damping .AND.                                      &
            intermediate_timestep_count == 1 )  THEN
          IF ( i == nxl .AND. j == nys )  THEN ! shifting only once per PE 

             DO  k = nzb, nzt
                IF ( w_subs(k) < 0.0_wp )  THEN      ! large-scale subsidence
                   var_mod(k) = var_init(k) - dt_3d * w_subs(k) *  &
                                     ( var_init(k+1) - var_init(k) ) * ddzu(k+1)
                ENDIF
             ENDDO
!
!--          At the upper boundary, the initial profile is shifted with aid of 
!--          the gradient tmp_grad. (This is ok if the gradients are linear.)
             IF ( w_subs(nzt) < 0.0_wp )  THEN
                tmp_grad = ( var_init(nzt+1) - var_init(nzt) ) * ddzu(nzt+1)
                var_mod(nzt+1) = var_init(nzt+1) -  &
                                     dt_3d * w_subs(nzt+1) * tmp_grad
             ENDIF
       

             DO  k = nzt+1, nzb+1, -1
                IF ( w_subs(k) >= 0.0_wp )  THEN  ! large-scale ascent
                   var_mod(k) = var_init(k) - dt_3d * w_subs(k) *  &
                                      ( var_init(k) - var_init(k-1) ) * ddzu(k)
                ENDIF
             ENDDO
!
!--          At the lower boundary shifting is not necessary because the 
!--          subsidence velocity w_subs(nzb) vanishes.
             IF ( w_subs(nzb+1) >= 0.0_wp )  THEN
                var_mod(nzb) = var_init(nzb)
             ENDIF

             var_init = var_mod 

          ENDIF
       ENDIF

 END SUBROUTINE subsidence_ij


 END MODULE subsidence_mod
