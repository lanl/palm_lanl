!> @file data_output_dvrp.f90
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
! $Id: data_output_dvrp.f90 3045 2018-05-28 07:55:41Z Giersch $
! Code adjusted according to PALM coding standards
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2274 2017-06-09 13:27:48Z Giersch
! Changed error messages
! 
! 2101 2017-01-05 16:42:31Z suehring
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1960 2016-07-12 16:34:24Z suehring
! Separate humidity and passive scalar
! 
! 1873 2016-04-18 14:50:06Z maronga
! Module renamed (removed _mod)
! 
! 
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
!
! 
! 1822 2016-04-07 07:49:42Z hoffmann
! Particles and tails removed.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1346 2014-03-27 13:18:20Z heinze
! Bugfix: REAL constants provided with KIND-attribute especially in call of 
! intrinsic function like MAX, MIN, SIGN
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1318 2014-03-17 13:35:16Z raasch
! module interfaces removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 828 2012-02-21 12:00:36Z raasch
! particle feature color renamed class
!
! Revision 1.1  2000/04/27 06:27:17  raasch
! Initial revision
!
!
! Description:
! ------------
!> Plot of isosurface and slicers with dvrp-software
!------------------------------------------------------------------------------!
 MODULE dvrp_color
 

    USE dvrp_variables
    
    USE kinds

    IMPLICIT NONE

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE color_dvrp( value, color )

       REAL(wp), INTENT(IN)  ::  value    !< 
       REAL(wp), INTENT(OUT) ::  color(4) !< 

       REAL(wp)              ::  scale    !< 

       scale = ( value - slicer_range_limits_dvrp(1,islice_dvrp) ) /           &
               ( slicer_range_limits_dvrp(2,islice_dvrp) -                     &
                 slicer_range_limits_dvrp(1,islice_dvrp) )

       scale = MODULO( 180.0_wp + 180.0_wp * scale, 360.0_wp )

       color = (/ scale, 0.5_wp, 1.0_wp, 0.0_wp /)

    END SUBROUTINE color_dvrp

 END MODULE dvrp_color


 RECURSIVE SUBROUTINE data_output_dvrp

#if defined( __dvrp_graphics )

    USE arrays_3d,                                                             &
        ONLY:  p, pt, q, ql, s, ts, u, us, v, w, zu
        
    USE cloud_parameters,                                                      &
        ONLY:  l_d_cp, pt_d_t
        
    USE constants,                                                             &
        ONLY:  pi
        
    USE control_parameters,                                                    &
        ONLY:  cloud_droplets, cloud_physics, do2d, do3d, humidity, ibc_uv_b,  &
               message_string, nz_do3d, passive_scalar, simulated_time,        &
               threshold
        
    USE cpulog,                                                                &
        ONLY:  log_point, log_point_s, cpu_log
        
    USE DVRP
    
    USE dvrp_color
        
    USE dvrp_variables
        
    USE grid_variables,                                                        &
        ONLY:  dx, dy
        
    USE indices,                                                               &
        ONLY:  nxl, nxr, nyn, nys, nzb
        
    USE kinds

    USE pegrid

    IMPLICIT NONE

    CHARACTER (LEN=2) ::  section_chr      !< 
    CHARACTER (LEN=6) ::  output_variable  !< 
    
    INTEGER(iwp) ::  c_mode           !<  
    INTEGER(iwp) ::  c_size_x         !< 
    INTEGER(iwp) ::  c_size_y         !< 
    INTEGER(iwp) ::  c_size_z         !< 
    INTEGER(iwp) ::  dvrp_nop         !< 
    INTEGER(iwp) ::  dvrp_not         !< 
    INTEGER(iwp) ::  gradient_normals !< 
    INTEGER(iwp) ::  i                !< 
    INTEGER(iwp) ::  ip               !< 
    INTEGER(iwp) ::  j                !< 
    INTEGER(iwp) ::  jp               !< 
    INTEGER(iwp) ::  k                !< 
    INTEGER(iwp) ::  l                !< 
    INTEGER(iwp) ::  m                !< 
    INTEGER(iwp) ::  n                !< 
    INTEGER(iwp) ::  n_isosurface     !< 
    INTEGER(iwp) ::  n_slicer         !< 
    INTEGER(iwp) ::  nn               !< 
    INTEGER(iwp) ::  section_mode     !< 
    INTEGER(iwp) ::  vn               !< 
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  p_c  !< 
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  p_t  !< 

    LOGICAL, DIMENSION(:), ALLOCATABLE ::  dvrp_mask  !< 

    REAL(sp) ::  slicer_position  !< 
    REAL(sp) ::  tmp_alpha        !< 
    REAL(sp) ::  tmp_alpha_w      !< 
    REAL(sp) ::  tmp_b            !< 
    REAL(sp) ::  tmp_c_alpha      !< 
    REAL(sp) ::  tmp_g            !< 
    REAL(sp) ::  tmp_norm         !< 
    REAL(sp) ::  tmp_pos          !< 
    REAL(sp) ::  tmp_r            !< 
    REAL(sp) ::  tmp_t            !< 
    REAL(sp) ::  tmp_th           !< 
    REAL(sp), DIMENSION(:),     ALLOCATABLE   ::  psize  !< 
    REAL(sp), DIMENSION(:),     ALLOCATABLE   ::  p_x    !< 
    REAL(sp), DIMENSION(:),     ALLOCATABLE   ::  p_y    !< 
    REAL(sp), DIMENSION(:),     ALLOCATABLE   ::  p_z    !< 
    REAL(sp), DIMENSION(:,:,:), ALLOCATABLE   ::  local_pf  !< 
    REAL(sp), DIMENSION(:,:,:,:), ALLOCATABLE ::  local_pfi !< 


    CALL cpu_log( log_point(27), 'data_output_dvrp', 'start' )

!
!-- Loop over all output modes choosed
    m            = 1
    n_isosurface = 0  ! isosurface counter (for threshold values and color)
    n_slicer     = 0  ! slice plane counter (for range of values)
    DO WHILE ( mode_dvrp(m) /= ' ' )
!
!--    Update of the steering variables
       IF ( .NOT. lock_steering_update )  THEN
!
!--       Set lock to avoid recursive calls of DVRP_STEERING_UPDATE
          lock_steering_update = .TRUE.
!          CALL DVRP_STEERING_UPDATE( m-1, data_output_dvrp )
          lock_steering_update = .FALSE.
       ENDIF

!
!--    Determine the variable which shall be plotted (in case of slicers or
!--    isosurfaces)
       IF ( mode_dvrp(m)(1:10) == 'isosurface' )  THEN
          READ ( mode_dvrp(m), '(10X,I2)' )  vn
          output_variable = do3d(0,vn)
          n_isosurface = n_isosurface + 1
       ELSEIF ( mode_dvrp(m)(1:6) == 'slicer' )  THEN
          READ ( mode_dvrp(m), '(6X,I2)' )  vn
          output_variable = do2d(0,vn)
          l = MAX( 2, LEN_TRIM( do2d(0,vn) ) )
          section_chr = do2d(0,vn)(l-1:l)
          SELECT CASE ( section_chr )
             CASE ( 'xy' )
                section_mode = 2
                slicer_position = zu(MIN( slicer_position_dvrp(m), nz_do3d ))
             CASE ( 'xz' )
                section_mode = 1
                slicer_position = slicer_position_dvrp(m) * dy
             CASE ( 'yz' )
                section_mode = 0
                slicer_position = slicer_position_dvrp(m) * dx
          END SELECT
       ENDIF

!
!--    Select the plot mode (in case of isosurface or slicer only if user has
!--    defined a variable which shall be plotted; otherwise do nothing)
       IF ( ( mode_dvrp(m)(1:10) == 'isosurface'  .OR.                         &
                  mode_dvrp(m)(1:6)  == 'slicer'           )                   &
                  .AND.  output_variable /= ' ' )  THEN

!
!--       Create an intermediate array, properly dimensioned for plot-output 
          ALLOCATE( local_pf(nxl_dvrp:nxr_dvrp+1,nys_dvrp:nyn_dvrp+1,          &
                             nzb:nz_do3d) )

!
!--       Move original array to intermediate array
          IF ( dvrp_overlap )  THEN

             SELECT CASE ( output_variable )

                CASE ( 'u', 'u_xy', 'u_xz', 'u_yz' )
                   DO  i = nxl_dvrp, nxr_dvrp+1
                      DO  j = nys_dvrp, nyn_dvrp+1
                         DO  k = nzb, nz_do3d
                            local_pf(i,j,k) = u(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
!
!--                Replace mirrored values at lower surface by real surface
!--                values
                   IF ( output_variable == 'u_xz'  .OR.                        &
                        output_variable == 'u_yz' )  THEN
                      IF ( ibc_uv_b == 0 )  local_pf(:,:,nzb) = 0.0_wp
                   ENDIF


                CASE ( 'v', 'v_xy', 'v_xz', 'v_yz' )
                   DO  i = nxl_dvrp, nxr_dvrp+1
                      DO  j = nys_dvrp, nyn_dvrp+1
                         DO  k = nzb, nz_do3d
                            local_pf(i,j,k) = v(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
!
!--                Replace mirrored values at lower surface by real surface
!--                values
                   IF ( output_variable == 'v_xz'  .OR.                        &
                        output_variable == 'v_yz' )  THEN
                      IF ( ibc_uv_b == 0 )  local_pf(:,:,nzb) = 0.0_wp
                   ENDIF

                CASE ( 'w', 'w_xy', 'w_xz', 'w_yz' )
                   DO  i = nxl_dvrp, nxr_dvrp+1
                      DO  j = nys_dvrp, nyn_dvrp+1
                         DO  k = nzb, nz_do3d
                            local_pf(i,j,k) = w(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
! Averaging for Langmuir circulation
!                   DO  k = nzb, nz_do3d
!                      DO  j = nys_dvrp+1, nyn_dvrp
!                         DO  i = nxl_dvrp, nxr_dvrp+1
!                            local_pf(i,j,k) = 0.25 * local_pf(i,j-1,k) + &
!                                              0.50 * local_pf(i,j,k)   + &
!                                              0.25 * local_pf(i,j+1,k)
!                         ENDDO
!                      ENDDO
!                   ENDDO

                CASE ( 'p', 'p_xy', 'p_xz', 'p_yz' )
                   DO  i = nxl_dvrp, nxr_dvrp+1
                      DO  j = nys_dvrp, nyn_dvrp+1
                         DO  k = nzb, nz_do3d
                            local_pf(i,j,k) = p(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO

                CASE ( 'pt', 'pt_xy', 'pt_xz', 'pt_yz' )
                   IF ( .NOT. cloud_physics ) THEN
                      DO  i = nxl_dvrp, nxr_dvrp+1
                         DO  j = nys_dvrp, nyn_dvrp+1
                            DO  k = nzb, nz_do3d
                               local_pf(i,j,k) = pt(k,j,i)
                            ENDDO
                         ENDDO
                      ENDDO
                   ELSE
                      DO  i = nxl_dvrp, nxr_dvrp+1
                         DO  j = nys_dvrp, nyn_dvrp+1
                            DO  k = nzb, nz_do3d
                               local_pf(i,j,k) = pt(k,j,i) + l_d_cp *          &
                                                 pt_d_t(k) * ql(k,j,i)
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDIF

                CASE ( 'q', 'q_xy', 'q_xz', 'q_yz' )
                   IF ( humidity )  THEN
                      DO  i = nxl_dvrp, nxr_dvrp+1
                         DO  j = nys_dvrp, nyn_dvrp+1
                            DO  k = nzb, nz_do3d
                               local_pf(i,j,k) = q(k,j,i)
                            ENDDO
                         ENDDO
                      ENDDO            
                   ELSE                    
                      message_string = 'if humidity = '    //                  & 
                            '.FALSE. output of ' // TRIM( output_variable ) // &
                            'is not provided' 
                      CALL message( 'data_output_dvrp', 'PA0183',&
                                                                 0, 0, 0, 6, 0 )
                   ENDIF
             
                CASE ( 'ql', 'ql_xy', 'ql_xz', 'ql_yz' )
                   IF ( cloud_physics  .OR.  cloud_droplets )  THEN
                      DO  i = nxl_dvrp, nxr_dvrp+1
                         DO  j = nys_dvrp, nyn_dvrp+1
                            DO  k = nzb, nz_do3d
                               local_pf(i,j,k) = ql(k,j,i)
                            ENDDO
                         ENDDO
                      ENDDO
                   ELSE                     
                      message_string = 'if cloud_physics = .FALSE. and ' //    & 
                                  'cloud_droplets = .FALSE. '
                                  'output of ' // TRIM( output_variable) //    &
                                  'is not provided' 
                      CALL message( 'data_output_dvrp', 'PA0184',&
                                                                 0, 0, 0, 6, 0 )
                   ENDIF

                CASE ( 's', 's_xy', 's_xz', 's_yz' )
                   IF ( passive_scalar )  THEN
                      DO  i = nxl_dvrp, nxr_dvrp+1
                         DO  j = nys_dvrp, nyn_dvrp+1
                            DO  k = nzb, nz_do3d
                               local_pf(i,j,k) = s(k,j,i)
                            ENDDO
                         ENDDO
                      ENDDO            
                   ELSE                    
                      message_string = 'if passive_scalar = '    //            & 
                            '.FALSE. output of ' // TRIM( output_variable ) //   &
                            'is not provided' 
                      CALL message( 'data_output_dvrp', 'PA0183',&
                                                                 0, 0, 0, 6, 0 )
                   ENDIF

                CASE ( 'u*_xy' )
                   DO  i = nxl_dvrp, nxr_dvrp+1
                      DO  j = nys_dvrp, nyn_dvrp+1
                         local_pf(i,j,nzb+1) = us(j,i)
                      ENDDO
                   ENDDO
                   slicer_position = zu(nzb+1)

                CASE ( 't*_xy' )
                   DO  i = nxl_dvrp, nxr_dvrp+1
                      DO  j = nys_dvrp, nyn_dvrp+1
                         local_pf(i,j,nzb+1) = ts(j,i)
                      ENDDO
                   ENDDO
                   slicer_position = zu(nzb+1)


                CASE DEFAULT
!
!--                The DEFAULT case is reached either if output_variable
!--                contains unsupported variable or if the user has coded a
!--                special case in the user interface. There, the subroutine
!--                user_data_output_dvrp checks which of these two conditions
!--                applies.
                   CALL user_data_output_dvrp( output_variable, local_pf )


             END SELECT

          ELSE
!
!--          No overlap of clipping domain with the current subdomain
             DO  i = nxl_dvrp, nxr_dvrp+1
                DO  j = nys_dvrp, nyn_dvrp+1
                   DO  k = nzb, nz_do3d
                      local_pf(i,j,k) = 0.0_wp
                   ENDDO
                ENDDO
             ENDDO

          ENDIF

          IF ( mode_dvrp(m)(1:10) == 'isosurface' )  THEN

!
!--          DVRP-Calls for plotting isosurfaces:
             CALL cpu_log( log_point_s(26), 'dvrp_isosurface', 'start' )

!
!--          Definition of isosurface color
             tmp_r = isosurface_color(1,n_isosurface)
             tmp_g = isosurface_color(2,n_isosurface)
             tmp_b = isosurface_color(3,n_isosurface)
             tmp_t = 0.0_wp
             CALL DVRP_MATERIAL_RGB( m-1, 1, tmp_r, tmp_g, tmp_b, tmp_t )

!
!--          Compute and plot isosurface in dvr-format
             CALL DVRP_DATA( m-1, local_pf, 1, nx_dvrp, ny_dvrp, nz_dvrp,      &
                             cyclic_dvrp, cyclic_dvrp, cyclic_dvrp )

             c_size_x = vc_size_x;  c_size_y = vc_size_y;  c_size_z = vc_size_z
             CALL DVRP_CLUSTER_SIZE( m-1, c_size_x, c_size_y, c_size_z )

             c_mode   = vc_mode 
             CALL DVRP_CLUSTERING_MODE( m-1, c_mode )

             gradient_normals = vc_gradient_normals
             CALL DVRP_GRADIENTNORMALS( m-1, gradient_normals )

!
!--          A seperate procedure for setting vc_alpha will be in the next
!--          version of libDVRP
             tmp_c_alpha = vc_alpha 
             CALL DVRP_THRESHOLD( -(m-1)-1, tmp_c_alpha )

             IF ( dvrp_overlap )  THEN
                tmp_th = threshold(n_isosurface)
             ELSE
                tmp_th = 1.0_wp  ! nothing is plotted because array values are 0
             ENDIF

             CALL DVRP_THRESHOLD( m-1, tmp_th )

             CALL DVRP_VISUALIZE( m-1, 21, dvrp_filecount )

             CALL cpu_log( log_point_s(26), 'dvrp_isosurface', 'stop' )

          ELSEIF ( mode_dvrp(m)(1:6) == 'slicer' )  THEN

!
!--          DVRP-Calls for plotting slicers:
             CALL cpu_log( log_point_s(27), 'dvrp_slicer', 'start' )

!
!--          Material and color definitions
             tmp_r = 0.0_wp;  tmp_g = 0.0_wp;  tmp_b = 0.0_wp;  tmp_t = 0.0_wp
             CALL DVRP_MATERIAL_RGB( m-1, 1, tmp_r, tmp_g, tmp_b, tmp_t )

             n_slicer = n_slicer + 1

!
!--           Using dolorfunction has not been properly tested
!             islice_dvrp = n_slicer
!             CALL DVRP_COLORFUNCTION( m-1, DVRP_CM_HLS, 25,                 &
!                                      slicer_range_limits_dvrp(:,n_slicer), &
!                                      color_dvrp )

!
!--          Set interval of values defining the colortable
             CALL set_slicer_attributes_dvrp( n_slicer )

!
!--          Create user-defined colortable
             CALL user_dvrp_coltab( 'slicer', output_variable )

             CALL DVRP_COLORTABLE_HLS( m-1, 1, interval_values_dvrp,           &
                                       interval_h_dvrp, interval_l_dvrp,       &
                                       interval_s_dvrp, interval_a_dvrp )

!
!--          Compute and plot slicer in dvr-format
             CALL DVRP_DATA( m-1, local_pf, 1, nx_dvrp, ny_dvrp, nz_dvrp,      &
                             cyclic_dvrp, cyclic_dvrp, cyclic_dvrp )
             tmp_pos = slicer_position
             CALL DVRP_SLICER( m-1, section_mode, tmp_pos )

             CALL DVRP_VISUALIZE( m-1, 2, dvrp_filecount )

             CALL cpu_log( log_point_s(27), 'dvrp_slicer', 'stop' )

          ENDIF

          DEALLOCATE( local_pf )

       ELSEIF ( mode_dvrp(m)(1:9) == 'pathlines' ) THEN

          ALLOCATE( local_pfi(4,nxl:nxr+1,nys:nyn+1,nzb:nz_do3d) )
          DO  i = nxl, nxr+1
             DO  j = nys, nyn+1
                DO  k = nzb, nz_do3d
                   local_pfi(1,i,j,k) = u(k,j,i)
                   local_pfi(2,i,j,k) = v(k,j,i)
                   local_pfi(3,i,j,k) = w(k,j,i)
                   tmp_norm           = SQRT( u(k,j,i) * u(k,j,i) +            &
                                              v(k,j,i) * v(k,j,i) +            &
                                              w(k,j,i) * w(k,j,i) )
                   tmp_alpha          = ACOS( 0.0_wp * u(k,j,i) / tmp_norm +   &
                                              0.0_wp * v(k,j,i) / tmp_norm -   &
                                              1.0_wp * w(k,j,i) / tmp_norm )
                   tmp_alpha_w        = tmp_alpha / pi * 180.0_wp
                   local_pfi(4,i,j,k) = tmp_alpha_w
                ENDDO
             ENDDO
          ENDDO

          CALL cpu_log( log_point_s(31), 'dvrp_pathlines', 'start' )

          CALL DVRP_DATA( m-1, local_pfi, 4, nx_dvrp, ny_dvrp, nz_dvrp,        &
                          cyclic_dvrp, cyclic_dvrp, cyclic_dvrp )
          CALL DVRP_VISUALIZE( m-1, 20, dvrp_filecount )

          CALL cpu_log( log_point_s(31), 'dvrp_pathlines', 'stop' )

          DEALLOCATE( local_pfi )

       ENDIF

       m = m + 1

    ENDDO

    dvrp_filecount = dvrp_filecount + 1

    CALL cpu_log( log_point(27), 'data_output_dvrp', 'stop' )

#endif
 END SUBROUTINE data_output_dvrp
