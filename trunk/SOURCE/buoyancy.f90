!> @file buoyancy.f90
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
! $Id: buoyancy.f90 2718 2018-01-02 08:49:38Z maronga $
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
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC version of subroutine removed
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
! 1374 2014-04-25 12:55:07Z raasch
! missing variables added to ONLY list
! 
! 1365 2014-04-22 15:03:56Z boeske
! Calculation of reference state in subroutine calc_mean_profile moved to
! subroutine time_integration,
! subroutine calc_mean_profile moved to new file calc_mean_profile.f90
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
! 1257 2013-11-08 15:18:40Z raasch
! vector length (32) removed from openacc clause
!
! 1241 2013-10-30 11:36:58Z heinze
! Generalize calc_mean_profile for wider use: use additional steering 
! character loc 
!
! 1179 2013-06-14 05:57:58Z raasch
! steering of reference state revised (var_reference and pr removed from
! parameter list), use_reference renamed use_single_reference_value
!
! 1171 2013-05-30 11:27:45Z raasch
! openacc statements added to use_reference-case in accelerator version
!
! 1153 2013-05-10 14:33:08Z raasch
! code adjustments of accelerator version required by PGI 12.3 / CUDA 5.0
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
! 1010 2012-09-20 07:59:54Z raasch
! cpp switch __nopointer added for pointer free version
!
! Revision 1.1  1997/08/29 08:56:48  raasch
! Initial revision
!
!
! Description:
! ------------
!> Buoyancy term of the third component of the equation of motion.
!> @attention Humidity is not regarded when using a sloping surface!
!------------------------------------------------------------------------------!
 MODULE buoyancy_mod
 

    PRIVATE
    PUBLIC buoyancy

    INTERFACE buoyancy
       MODULE PROCEDURE buoyancy
       MODULE PROCEDURE buoyancy_ij
    END INTERFACE buoyancy

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE buoyancy( var, wind_component )

       USE arrays_3d,                                                          &
           ONLY:  pt, pt_slope_ref, ref_state, tend

       USE control_parameters,                                                 &
           ONLY:  atmos_ocean_sign, cos_alpha_surface, g, message_string,      &
                  pt_surface, sin_alpha_surface, sloping_surface

       USE indices,                                                            &
           ONLY:  nxl, nxlg, nxlu, nxr, nxrg, nyn, nyng, nys, nysg, nzb,       &
                  nzt, wall_flags_0

       USE kinds

       USE pegrid


       IMPLICIT NONE

       INTEGER(iwp) ::  i              !<
       INTEGER(iwp) ::  j              !<
       INTEGER(iwp) ::  k              !<
       INTEGER(iwp) ::  wind_component !<
       
#if defined( __nopointer )
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  var !<
#else
       REAL(wp), DIMENSION(:,:,:), POINTER ::  var
#endif

       IF ( .NOT. sloping_surface )  THEN
!
!--       Normal case: horizontal surface
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt-1

                   tend(k,j,i) = tend(k,j,i) + atmos_ocean_sign * g * 0.5_wp *  &
                          (                                                     &
                             ( var(k,j,i)   - ref_state(k) )   / ref_state(k) + &
                             ( var(k+1,j,i) - ref_state(k+1) ) / ref_state(k+1) &
                          ) * MERGE( 1.0_wp, 0.0_wp,                            &
                                     BTEST( wall_flags_0(k,j,i), 0 ) )
                ENDDO
             ENDDO
          ENDDO

       ELSE
!
!--       Buoyancy term for a surface with a slope in x-direction. The equations
!--       for both the u and w velocity-component contain proportionate terms.
!--       Temperature field at time t=0 serves as environmental temperature.
!--       Reference temperature (pt_surface) is the one at the lower left corner
!--       of the total domain.
          IF ( wind_component == 1 )  THEN

             DO  i = nxlu, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt-1
                      tend(k,j,i) = tend(k,j,i) + g * sin_alpha_surface *         &
                           0.5_wp * ( ( pt(k,j,i-1)         + pt(k,j,i)         ) &
                                    - ( pt_slope_ref(k,i-1) + pt_slope_ref(k,i) ) &
                                    ) / pt_surface                                &
                                      * MERGE( 1.0_wp, 0.0_wp,                    &
                                               BTEST( wall_flags_0(k,j,i), 0 ) )
                   ENDDO
                ENDDO
             ENDDO

          ELSEIF ( wind_component == 3 )  THEN

             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt-1
                      tend(k,j,i) = tend(k,j,i) + g * cos_alpha_surface *         &
                           0.5_wp * ( ( pt(k,j,i)         + pt(k+1,j,i)         ) &
                                    - ( pt_slope_ref(k,i) + pt_slope_ref(k+1,i) ) &
                                    ) / pt_surface                                &
                                      * MERGE( 1.0_wp, 0.0_wp,                    &
                                               BTEST( wall_flags_0(k,j,i), 0 ) )
                   ENDDO
                ENDDO
            ENDDO

          ELSE
             
             WRITE( message_string, * ) 'no term for component "',             &
                                       wind_component,'"'
             CALL message( 'buoyancy', 'PA0159', 1, 2, 0, 6, 0 )

          ENDIF

       ENDIF

    END SUBROUTINE buoyancy


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!> @attention PGI-compiler creates SIGFPE if opt>1 is used! Therefore, opt=1 is
!>            forced by compiler-directive.
!------------------------------------------------------------------------------!
!pgi$r opt=1
    SUBROUTINE buoyancy_ij( i, j, var, wind_component )

       USE arrays_3d,                                                          &
           ONLY:  pt, pt_slope_ref, ref_state, tend

       USE control_parameters,                                                 &
           ONLY:  atmos_ocean_sign, cos_alpha_surface, g, message_string,      &
                  pt_surface, sin_alpha_surface, sloping_surface

       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nyng, nysg, nzb, nzt, wall_flags_0

       USE kinds

       USE pegrid


       IMPLICIT NONE

       INTEGER(iwp) ::  i              !<
       INTEGER(iwp) ::  j              !<
       INTEGER(iwp) ::  k              !<
       INTEGER(iwp) ::  pr             !<
       INTEGER(iwp) ::  wind_component !<
       
#if defined( __nopointer )
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  var !<
#else
       REAL(wp), DIMENSION(:,:,:), POINTER ::  var
#endif


       IF ( .NOT. sloping_surface )  THEN
!
!--       Normal case: horizontal surface
          DO  k = nzb+1, nzt-1
              tend(k,j,i) = tend(k,j,i) + atmos_ocean_sign * g * 0.5_wp * (    &
                        ( var(k,j,i)   - ref_state(k)   ) / ref_state(k)   +   &
                        ( var(k+1,j,i) - ref_state(k+1) ) / ref_state(k+1)     &
                                                                          )    &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 0 ) )
          ENDDO

       ELSE
!
!--       Buoyancy term for a surface with a slope in x-direction. The equations
!--       for both the u and w velocity-component contain proportionate terms.
!--       Temperature field at time t=0 serves as environmental temperature.
!--       Reference temperature (pt_surface) is the one at the lower left corner
!--       of the total domain.
          IF ( wind_component == 1 )  THEN

             DO  k = nzb+1, nzt-1
                tend(k,j,i) = tend(k,j,i) + g * sin_alpha_surface *               &
                           0.5_wp * ( ( pt(k,j,i-1)         + pt(k,j,i)         ) &
                                    - ( pt_slope_ref(k,i-1) + pt_slope_ref(k,i) ) &
                                    ) / pt_surface                                &
                                      * MERGE( 1.0_wp, 0.0_wp,                    &
                                               BTEST( wall_flags_0(k,j,i), 0 ) )
             ENDDO

          ELSEIF ( wind_component == 3 )  THEN

             DO  k = nzb+1, nzt-1
                tend(k,j,i) = tend(k,j,i) + g * cos_alpha_surface *               &
                           0.5_wp * ( ( pt(k,j,i)         + pt(k+1,j,i)         ) &
                                    - ( pt_slope_ref(k,i) + pt_slope_ref(k+1,i) ) &
                                    ) / pt_surface                                &
                                      * MERGE( 1.0_wp, 0.0_wp,                    &
                                               BTEST( wall_flags_0(k,j,i), 0 ) )
             ENDDO

          ELSE

             WRITE( message_string, * ) 'no term for component "',             &
                                       wind_component,'"'
             CALL message( 'buoyancy', 'PA0159', 1, 2, 0, 6, 0 )

          ENDIF

       ENDIF

    END SUBROUTINE buoyancy_ij

 END MODULE buoyancy_mod
