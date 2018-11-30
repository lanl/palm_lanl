!> @file diffusion_v.f90
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
! $Id: diffusion_v.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2638 2017-11-23 12:44:23Z raasch
! bugfix for constant top momentumflux
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
! 1740 2016-01-13 08:19:40Z raasch
! unnecessary calculations of kmzm and kmzp in wall bounded parts removed
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
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
! openacc loop and loop vector clauses removed, declare create moved after
! the FORTRAN declaration statement
!
! 1128 2013-04-12 06:19:32Z raasch
! loop index bounds in accelerator version replaced by i_left, i_right, j_south,
! j_north
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1015 2012-09-27 09:23:24Z raasch
! accelerator version (*_acc) added
!
! 1001 2012-09-13 14:08:46Z raasch
! arrays comunicated by module instead of parameter list
!
! 978 2012-08-09 08:28:32Z fricke
! outflow damping layer removed
! kmxm_x/_y and kmxp_x/_y change to kmxm and kmxp
!
! Revision 1.1  1997/09/12 06:24:01  raasch
! Initial revision
!
!
! Description:
! ------------
!> Diffusion term of the v-component
!------------------------------------------------------------------------------!
 MODULE diffusion_v_mod
 

    PRIVATE
    PUBLIC diffusion_v

    INTERFACE diffusion_v
       MODULE PROCEDURE diffusion_v
       MODULE PROCEDURE diffusion_v_ij
    END INTERFACE diffusion_v

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE diffusion_v

       USE arrays_3d,                                                          &
           ONLY:  ddzu, ddzw, km, tend, u, v, w, drho_air, rho_air_zw
       
       USE control_parameters,                                                 &
           ONLY:  constant_top_momentumflux, use_surface_fluxes,               &
                  use_top_fluxes
       
       USE grid_variables,                                                     &
           ONLY:  ddx, ddy, ddy2
       
       USE indices,                                                            &
           ONLY:  nxl, nxr, nyn, nys, nysv, nzb, nzt, wall_flags_0
       
       USE kinds

       USE surface_mod,                                                        &
           ONLY :  surf_def_h, surf_def_v

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  l             !< running index of surface type, south- or north-facing wall
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp)     ::  flag          !< flag to mask topography grid points
       REAL(wp)     ::  kmxm          !< 
       REAL(wp)     ::  kmxp          !< 
       REAL(wp)     ::  kmzm          !< 
       REAL(wp)     ::  kmzp          !< 
       REAL(wp)     ::  mask_bottom   !< flag to mask vertical upward-facing surface  
       REAL(wp)     ::  mask_east     !< flag to mask vertical surface south of the grid point 
       REAL(wp)     ::  mask_west     !< flag to mask vertical surface north of the grid point 
       REAL(wp)     ::  mask_top      !< flag to mask vertical downward-facing surface      

       DO  i = nxl, nxr
          DO  j = nysv, nyn
!
!--          Compute horizontal diffusion
             DO  k = nzb+1, nzt

!
!--             Predetermine flag to mask topography and wall-bounded grid points. 
!--             It is sufficient to masked only east- and west-facing surfaces, which
!--             need special treatment for the v-component. 
                flag      = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i),   2 ) ) 
                mask_east = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i+1), 2 ) )
                mask_west = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i-1), 2 ) )
!
!--             Interpolate eddy diffusivities on staggered gridpoints
                kmxp = 0.25_wp * ( km(k,j,i)+km(k,j,i+1)+km(k,j-1,i)+km(k,j-1,i+1) )
                kmxm = 0.25_wp * ( km(k,j,i)+km(k,j,i-1)+km(k,j-1,i)+km(k,j-1,i-1) )

                tend(k,j,i) = tend(k,j,i) +    (                             &
                          mask_east * kmxp * (                               &
                                 ( v(k,j,i+1) - v(k,j,i)     ) * ddx         &
                               + ( u(k,j,i+1) - u(k,j-1,i+1) ) * ddy         &
                                             )                               &
                        - mask_west * kmxm * (                               &
                                 ( v(k,j,i) - v(k,j,i-1) ) * ddx             &
                               + ( u(k,j,i) - u(k,j-1,i) ) * ddy             &
                                             )                               &
                                               ) * ddx  * flag               &
                                    + 2.0_wp * (                             &
                                  km(k,j,i)   * ( v(k,j+1,i) - v(k,j,i)   )  &
                                - km(k,j-1,i) * ( v(k,j,i)   - v(k,j-1,i) )  &
                                               ) * ddy2 * flag

             ENDDO

!
!--          Add horizontal momentum flux v'u' at east- (l=2) and west-facing (l=3)
!--          surfaces. Note, in the the flat case, loops won't be entered as 
!--          start_index > end_index. Furtermore, note, no vertical natural surfaces
!--          so far.           
!--          Default-type surfaces
             DO  l = 2, 3
                surf_s = surf_def_v(l)%start_index(j,i)
                surf_e = surf_def_v(l)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k           = surf_def_v(l)%k(m)
                   tend(k,j,i) = tend(k,j,i) +                                 &
                                    surf_def_v(l)%mom_flux_uv(m) * ddx
                ENDDO   
             ENDDO
!
!--          Compute vertical diffusion. In case of simulating a surface layer,
!--          respective grid diffusive fluxes are masked (flag 10) within this 
!--          loop, and added further below, else, simple gradient approach is
!--          applied. Model top is also mask if top-momentum flux is given. 
             DO  k = nzb+1, nzt
!
!--             Determine flags to mask topography below and above. Flag 2 is 
!--             used to mask topography in general, while flag 8 implies also
!--             information about use_surface_fluxes. Flag 9 is used to control 
!--             momentum flux at model top.  
                mask_bottom = MERGE( 1.0_wp, 0.0_wp,                           &
                                     BTEST( wall_flags_0(k-1,j,i), 8 ) ) 
                mask_top    = MERGE( 1.0_wp, 0.0_wp,                           &
                                     BTEST( wall_flags_0(k+1,j,i), 8 ) ) *     &
                              MERGE( 1.0_wp, 0.0_wp,                           &
                                     BTEST( wall_flags_0(k+1,j,i), 9 ) ) 
                flag        = MERGE( 1.0_wp, 0.0_wp,                           &
                                     BTEST( wall_flags_0(k,j,i), 2 ) ) 
!
!--             Interpolate eddy diffusivities on staggered gridpoints
                kmzp = 0.25_wp * &
                       ( km(k,j,i)+km(k+1,j,i)+km(k,j-1,i)+km(k+1,j-1,i) )
                kmzm = 0.25_wp * &
                       ( km(k,j,i)+km(k-1,j,i)+km(k,j-1,i)+km(k-1,j-1,i) )

                tend(k,j,i) = tend(k,j,i)                                      &
                      & + ( kmzp * ( ( v(k+1,j,i) - v(k,j,i) ) * ddzu(k+1)     &
                      &            + ( w(k,j,i) - w(k,j-1,i) ) * ddy           &
                      &            ) * rho_air_zw(k)   * mask_top              &
                      &   - kmzm * ( ( v(k,j,i)   - v(k-1,j,i)   ) * ddzu(k)   &
                      &            + ( w(k-1,j,i) - w(k-1,j-1,i) ) * ddy       &
                      &            ) * rho_air_zw(k-1) * mask_bottom           &
                      &   ) * ddzw(k) * drho_air(k) * flag
             ENDDO

!
!--          Vertical diffusion at the first grid point above the surface,
!--          if the momentum flux at the bottom is given by the Prandtl law
!--          or if it is prescribed by the user.
!--          Difference quotient of the momentum flux is not formed over
!--          half of the grid spacing (2.0*ddzw(k)) any more, since the
!--          comparison with other (LES) models showed that the values of
!--          the momentum flux becomes too large in this case.
             IF ( use_surface_fluxes )  THEN
!
!--             Default-type surfaces, upward-facing
                surf_s = surf_def_h(0)%start_index(j,i)
                surf_e = surf_def_h(0)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k   = surf_def_h(0)%k(m)

                   tend(k,j,i) = tend(k,j,i)                                   &
                        + ( - ( - surf_def_h(0)%vsws(m) )                      &
                          ) * ddzw(k) * drho_air(k)
                ENDDO
!
!--             Default-type surfaces, dowward-facing
                surf_s = surf_def_h(1)%start_index(j,i)
                surf_e = surf_def_h(1)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k   = surf_def_h(1)%k(m)

                   tend(k,j,i) = tend(k,j,i)                                   &
                        + ( - surf_def_h(1)%vsws(m)                            &
                          ) * ddzw(k) * drho_air(k)
                ENDDO
             ENDIF
!
!--          Add momentum flux at model top
             IF ( use_top_fluxes  .AND.  constant_top_momentumflux )  THEN
                surf_s = surf_def_h(2)%start_index(j,i)
                surf_e = surf_def_h(2)%end_index(j,i)
                DO  m = surf_s, surf_e

                   k   = surf_def_h(2)%k(m)

                   tend(k,j,i) = tend(k,j,i)                                   &
                           + ( - surf_def_h(2)%vsws(m) ) * ddzw(k) * drho_air(k)
                ENDDO
             ENDIF

          ENDDO
       ENDDO

    END SUBROUTINE diffusion_v


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE diffusion_v_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  ddzu, ddzw, km, tend, u, v, w, drho_air, rho_air_zw
       
       USE control_parameters,                                                 &
           ONLY:  constant_top_momentumflux, use_surface_fluxes,               &
                  use_top_fluxes
       
       USE grid_variables,                                                     &
           ONLY:  ddx, ddy, ddy2
       
       USE indices,                                                            &
           ONLY:  nzb, nzt, wall_flags_0
       
       USE kinds

       USE surface_mod,                                                        &
           ONLY :  surf_def_h, surf_def_v

       IMPLICIT NONE


       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  l             !< running index of surface type, south- or north-facing wall
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp)     ::  flag          !< flag to mask topography grid points
       REAL(wp)     ::  kmxm          !< 
       REAL(wp)     ::  kmxp          !< 
       REAL(wp)     ::  kmzm          !< 
       REAL(wp)     ::  kmzp          !< 
       REAL(wp)     ::  mask_bottom   !< flag to mask vertical upward-facing surface  
       REAL(wp)     ::  mask_east     !< flag to mask vertical surface south of the grid point 
       REAL(wp)     ::  mask_west     !< flag to mask vertical surface north of the grid point 
       REAL(wp)     ::  mask_top      !< flag to mask vertical downward-facing surface

!
!--    Compute horizontal diffusion
       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography and wall-bounded grid points. 
!--       It is sufficient to masked only east- and west-facing surfaces, which
!--       need special treatment for the v-component. 
          flag      = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i),   2 ) ) 
          mask_east = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i+1), 2 ) )
          mask_west = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i-1), 2 ) )
!
!--       Interpolate eddy diffusivities on staggered gridpoints
          kmxp = 0.25_wp * ( km(k,j,i)+km(k,j,i+1)+km(k,j-1,i)+km(k,j-1,i+1) )
          kmxm = 0.25_wp * ( km(k,j,i)+km(k,j,i-1)+km(k,j-1,i)+km(k,j-1,i-1) )

          tend(k,j,i) = tend(k,j,i) +          (                             &
                          mask_east * kmxp * (                               &
                                 ( v(k,j,i+1) - v(k,j,i)     ) * ddx         &
                               + ( u(k,j,i+1) - u(k,j-1,i+1) ) * ddy         &
                                             )                               &
                        - mask_west * kmxm * (                               &
                                 ( v(k,j,i) - v(k,j,i-1) ) * ddx             &
                               + ( u(k,j,i) - u(k,j-1,i) ) * ddy             &
                                             )                               &
                                               ) * ddx  * flag               &
                                    + 2.0_wp * (                             &
                                  km(k,j,i)   * ( v(k,j+1,i) - v(k,j,i)   )  &
                                - km(k,j-1,i) * ( v(k,j,i)   - v(k,j-1,i) )  &
                                               ) * ddy2 * flag
       ENDDO

!
!--    Add horizontal momentum flux v'u' at east- (l=2) and west-facing (l=3)
!--    surfaces. Note, in the the flat case, loops won't be entered as 
!--    start_index > end_index. Furtermore, note, no vertical natural surfaces
!--    so far.           
!--    Default-type surfaces
       DO  l = 2, 3
          surf_s = surf_def_v(l)%start_index(j,i)
          surf_e = surf_def_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_def_v(l)%k(m)
             tend(k,j,i) = tend(k,j,i) + surf_def_v(l)%mom_flux_uv(m) * ddx
          ENDDO   
       ENDDO
!
!--    Compute vertical diffusion. In case of simulating a surface layer,
!--    respective grid diffusive fluxes are masked (flag 8) within this 
!--    loop, and added further below, else, simple gradient approach is
!--    applied. Model top is also mask if top-momentum flux is given.
       DO  k = nzb+1, nzt
!
!--       Determine flags to mask topography below and above. Flag 2 is 
!--       used to mask topography in general, while flag 10 implies also
!--       information about use_surface_fluxes. Flag 9 is used to control 
!--       momentum flux at model top.  
          mask_bottom = MERGE( 1.0_wp, 0.0_wp,                                 &
                               BTEST( wall_flags_0(k-1,j,i), 8 ) ) 
          mask_top    = MERGE( 1.0_wp, 0.0_wp,                                 &
                               BTEST( wall_flags_0(k+1,j,i), 8 ) ) *           &
                        MERGE( 1.0_wp, 0.0_wp,                                 &
                               BTEST( wall_flags_0(k+1,j,i), 9 ) ) 
          flag        = MERGE( 1.0_wp, 0.0_wp,                                 &
                               BTEST( wall_flags_0(k,j,i), 2 ) )
!
!--       Interpolate eddy diffusivities on staggered gridpoints
          kmzp = 0.25_wp * ( km(k,j,i)+km(k+1,j,i)+km(k,j-1,i)+km(k+1,j-1,i) )
          kmzm = 0.25_wp * ( km(k,j,i)+km(k-1,j,i)+km(k,j-1,i)+km(k-1,j-1,i) )

          tend(k,j,i) = tend(k,j,i)                                            &
                      & + ( kmzp * ( ( v(k+1,j,i) - v(k,j,i) ) * ddzu(k+1)     &
                      &            + ( w(k,j,i) - w(k,j-1,i) ) * ddy           &
                      &            ) * rho_air_zw(k)   * mask_top              &
                      &   - kmzm * ( ( v(k,j,i)   - v(k-1,j,i)   ) * ddzu(k)   &
                      &            + ( w(k-1,j,i) - w(k-1,j-1,i) ) * ddy       &
                      &            ) * rho_air_zw(k-1) * mask_bottom           &
                      &   ) * ddzw(k) * drho_air(k) * flag
       ENDDO

!
!--    Vertical diffusion at the first grid point above the surface, if the
!--    momentum flux at the bottom is given by the Prandtl law or if it is
!--    prescribed by the user.
!--    Difference quotient of the momentum flux is not formed over half of
!--    the grid spacing (2.0*ddzw(k)) any more, since the comparison with 
!--    other (LES) models showed that the values of the momentum flux becomes
!--    too large in this case.
       IF ( use_surface_fluxes )  THEN
!
!--       Default-type surfaces, upward-facing
          surf_s = surf_def_h(0)%start_index(j,i)
          surf_e = surf_def_h(0)%end_index(j,i)
          DO  m = surf_s, surf_e
             k   = surf_def_h(0)%k(m)

             tend(k,j,i) = tend(k,j,i)                                         &
                        + ( - ( - surf_def_h(0)%vsws(m) )                      &
                          ) * ddzw(k) * drho_air(k)
          ENDDO
!
!--       Default-type surfaces, dowward-facing
          surf_s = surf_def_h(1)%start_index(j,i)
          surf_e = surf_def_h(1)%end_index(j,i)
          DO  m = surf_s, surf_e
             k   = surf_def_h(1)%k(m)

             tend(k,j,i) = tend(k,j,i)                                         &
                        + ( - surf_def_h(1)%vsws(m)                            &
                          ) * ddzw(k) * drho_air(k)
          ENDDO
       ENDIF
!
!--    Add momentum flux at model top
       IF ( use_top_fluxes  .AND.  constant_top_momentumflux )  THEN
          surf_s = surf_def_h(2)%start_index(j,i)
          surf_e = surf_def_h(2)%end_index(j,i)
          DO  m = surf_s, surf_e

             k   = surf_def_h(2)%k(m)

             tend(k,j,i) = tend(k,j,i)                                         &
                           + ( - surf_def_h(2)%vsws(m) ) * ddzw(k) * drho_air(k)
          ENDDO
       ENDIF

    END SUBROUTINE diffusion_v_ij

 END MODULE diffusion_v_mod
