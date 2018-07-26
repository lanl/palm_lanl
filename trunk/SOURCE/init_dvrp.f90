!> @file init_dvrp.f90
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
! $Id: init_dvrp.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2516 2017-10-04 11:03:04Z suehring
! Remove tabs
! 
! 2514 2017-10-04 09:52:37Z suehring
! NEC related cpp directives removed
! 
! 2298 2017-06-29 09:28:18Z raasch
! MPI2 related part removed
! 
! 2101 2017-01-05 16:42:31Z suehring
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1808 2016-04-05 19:44:00Z raasch
! routine local_getenv replaced by standard FORTRAN routine
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute 
!
! 1322 2014-03-20 16:38:49Z raasch
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
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! Revision 1.1  2000/04/27 06:24:39  raasch
! Initial revision
!
!
! Description:
! ------------
!> Initializing actions needed when using dvrp-software
!------------------------------------------------------------------------------!
  SUBROUTINE init_dvrp
 
#if defined( __dvrp_graphics )

    USE arrays_3d,                                                             &
        ONLY:  zu
        
    USE DVRP
    
    USE dvrp_variables
    
    USE grid_variables,                                                        &
        ONLY:  dx, dy
        
    USE indices,                                                               &
        ONLY:  nx, nxl, nxr, ny, nyn, nys, nzb, nzb_s_inner
        
    USE kinds
    
    USE pegrid
    
    USE control_parameters,                                                               &
        ONLY:  message_string, nz_do3d, run_identifier, topography

    IMPLICIT NONE

    CHARACTER (LEN=2)  ::  section_chr      !<
    CHARACTER (LEN=3)  ::  prefix_chr       !<
    CHARACTER (LEN=80) ::  dvrp_file_local  !<
    
    INTEGER(iwp) ::  cluster_mode      !<
    INTEGER(iwp) ::  cluster_size_x    !<
    INTEGER(iwp) ::  cluster_size_y    !<
    INTEGER(iwp) ::  cluster_size_z    !<
    INTEGER(iwp) ::  gradient_normals  !<
    INTEGER(iwp) ::  i                 !<
    INTEGER(iwp) ::  j                 !<
    INTEGER(iwp) ::  k                 !<
    INTEGER(iwp) ::  l                 !<
    INTEGER(iwp) ::  m                 !<
    INTEGER(iwp) ::  nx_dvrp_l         !<
    INTEGER(iwp) ::  nx_dvrp_r         !<
    INTEGER(iwp) ::  ny_dvrp_n         !<
    INTEGER(iwp) ::  ny_dvrp_s         !<
    INTEGER(iwp) ::  pn                !<
    INTEGER(iwp) ::  tv                !<
    INTEGER(iwp) ::  vn                !<
                     
    LOGICAL  ::  allocated  !<
    
    REAL(sp) ::  center(3)      !<
    REAL(sp) ::  cluster_alpha  !<
    REAL(sp) ::  distance       !<
    REAL(sp) ::  tmp_b          !<
    REAL(sp) ::  tmp_g          !<
    REAL(sp) ::  tmp_r          !<
    REAL(sp) ::  tmp_t          !<
    REAL(sp) ::  tmp_th         !<
    REAL(sp) ::  tmp_thr        !<
    REAL(sp) ::  tmp_x1         !<
    REAL(sp) ::  tmp_x2         !<
    REAL(sp) ::  tmp_y1         !<
    REAL(sp) ::  tmp_y2         !<
    REAL(sp) ::  tmp_z1         !<
    REAL(sp) ::  tmp_z2         !<
    REAL(sp) ::  tmp_1          !<
    REAL(sp) ::  tmp_2          !<
    REAL(sp) ::  tmp_3          !<
    REAL(sp) ::  tmp_4          !<
    REAL(sp) ::  tmp_5          !<
    REAL(sp) ::  tmp_6          !<
    REAL(sp) ::  tmp_7          !<

    REAL(sp), DIMENSION(:,:,:), ALLOCATABLE ::  local_pf  !<

    TYPE(CSTRING), SAVE ::  dvrp_directory_c   !<
    TYPE(CSTRING), SAVE ::  dvrp_file_c        !<
    TYPE(CSTRING), SAVE ::  dvrp_file_local_c  !<
    TYPE(CSTRING), SAVE ::  dvrp_host_c        !<
    TYPE(CSTRING), SAVE ::  dvrp_password_c    !<
    TYPE(CSTRING), SAVE ::  dvrp_username_c    !<
    TYPE(CSTRING), SAVE ::  name_c             !<

!
!-- Set clipping to default (total domain), if not set by user
    IF ( clip_dvrp_l == 9999999.9_wp )  clip_dvrp_l = 0.0_wp
    IF ( clip_dvrp_r == 9999999.9_wp )  clip_dvrp_r = ( nx + 1 ) * dx
    IF ( clip_dvrp_s == 9999999.9_wp )  clip_dvrp_s = 0.0_wp
    IF ( clip_dvrp_n == 9999999.9_wp )  clip_dvrp_n = ( ny + 1 ) * dy

!
!-- Calculate the clipping index limits
    nx_dvrp_l = clip_dvrp_l / dx
    nx_dvrp_r = clip_dvrp_r / dx
    ny_dvrp_s = clip_dvrp_s / dy
    ny_dvrp_n = clip_dvrp_n / dy

    IF ( nx_dvrp_l < nxr  .AND.  nx_dvrp_r > nxl  .AND. &
         ny_dvrp_s < nyn  .AND.  ny_dvrp_n > nys )  THEN

       dvrp_overlap = .TRUE.
       nxl_dvrp = MAX( nxl, nx_dvrp_l )
       nxr_dvrp = MIN( nxr, nx_dvrp_r )
       nys_dvrp = MAX( nys, ny_dvrp_s )
       nyn_dvrp = MIN( nyn, ny_dvrp_n )

       IF ( nxl_dvrp == nxl  .AND.  nxr_dvrp == nxr  .AND.  &
            nys_dvrp == nys  .AND.  nyn_dvrp == nyn )  THEN
          dvrp_total_overlap = .TRUE.
       ELSE
          dvrp_total_overlap = .FALSE.
       ENDIF

    ELSE
!
!--    This subdomain does not overlap with the clipping area. Define an
!--    arbitrary (small) domain within in the clipping area.
       dvrp_overlap       = .FALSE.
       dvrp_total_overlap = .FALSE.
!       nxl_dvrp = nx_dvrp_l
!       nxr_dvrp = nxl_dvrp + 4
!       nys_dvrp = ny_dvrp_s
!       nyn_dvrp = nys_dvrp + 4
       nxl_dvrp = nxl
       nxr_dvrp = MIN( nxl+4, nxr )
       nys_dvrp = nys
       nyn_dvrp = MIN( nys+4, nyn )

    ENDIF

!
!-- Set the maximum time the program can be suspended on user request (by
!-- dvrp steering). This variable is defined in module DVRP.
    DVRP_MAX_SUSPEND_TIME = 7200

!
!-- Allocate array holding the names and limits of the steering variables
!-- (must have the same number of elements as array mode_dvrp!)
    ALLOCATE( steering_dvrp(10) )

!
!-- Check, if output parameters are given and/or allowed
!-- and set default-values, where necessary
    IF ( dvrp_username == ' ' )  THEN
        message_string = 'dvrp_username is undefined'
        CALL message( 'init_dvrp', 'PA0195', 1, 2, 0, 6, 0 )         
    ENDIF

    IF ( dvrp_output /= 'ftp'  .AND.  dvrp_output /= 'rtsp'  .AND. &
         dvrp_output /= 'local' )  THEN
       message_string = 'dvrp_output="' // TRIM( dvrp_output ) // &
                        '" not allowed'
       CALL message( 'init_dvrp', 'PA0196', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( dvrp_directory == 'default' )  THEN
       dvrp_directory = TRIM( dvrp_username ) // '/' // TRIM( run_identifier )
    ENDIF

!
!-- A local dvrserver running always outputs on temporary directory DATA_DVR
    IF ( local_dvrserver_running )  THEN
       dvrp_directory = 'DATA_DVR'
    ENDIF

    IF ( dvrp_output /= 'local' )  THEN
       IF ( dvrp_file /= 'default'  .AND.  dvrp_file /= '/dev/null' )  THEN
          message_string = 'dvrp_file="' // TRIM( dvrp_file ) // '" not allowed'
          CALL message( 'init_dvrp', 'PA0197', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Strings are assigned to strings of special type which have a CHAR( 0 )
!-- (C end-of-character symbol) at their end. This is needed when strings are
!-- passed to C routines.
    dvrp_directory_c = dvrp_directory
    dvrp_file_c      = dvrp_file
    dvrp_host_c      = dvrp_host
    dvrp_password_c  = dvrp_password
    dvrp_username_c  = dvrp_username

!
!-- Loop over all output modes choosed
    m = 1
    allocated = .FALSE.
    DO WHILE ( mode_dvrp(m) /= ' ' )
    
!
!--    Check, if mode is allowed
       IF ( mode_dvrp(m)(1:10) /= 'isosurface'  .AND. &
            mode_dvrp(m)(1:6)  /= 'slicer'      .AND. &
            mode_dvrp(m)(1:9)  /= 'particles'   .AND. &
            mode_dvrp(m)(1:9)  /= 'pathlines' )  THEN

          message_string = 'mode_dvrp="' // TRIM( mode_dvrp(m) ) // &
                           '" not allowed'
          CALL message( 'init_dvrp', 'PA0198', 1, 2, 0, 6, 0 )
          CALL local_stop

       ENDIF
!
!--    Determine prefix for dvrp_file
       WRITE ( prefix_chr, '(I2.2,''_'')' )  m
!
!--    Camera position must be computed and written on file when no dvrp-output
!--    has been generated so far (in former runs)
!       IF ( dvrp_filecount == 0 )  THEN
!
!--       Compute center of domain and distance of camera from center
          center(1) = ( clip_dvrp_l + clip_dvrp_r ) * 0.5_wp * superelevation_x
          center(2) = ( clip_dvrp_s + clip_dvrp_n ) * 0.5_wp * superelevation_y
          center(3) = ( zu(nz_do3d) - zu(nzb) ) * 0.5_wp * superelevation
          distance  = 1.5_wp * MAX(                                            &
                           (clip_dvrp_r-clip_dvrp_l) * superelevation_x,       &
                           (clip_dvrp_n-clip_dvrp_s) * superelevation_y,       &
                           ( zu(nz_do3d) - zu(nzb) ) * superelevation          &
                                  )

!
!--       Write camera position on file
          CALL DVRP_INIT( m-1, 0 )

!
!--       Create filename for camera
          IF ( dvrp_output == 'rtsp' )  THEN

             dvrp_file   = prefix_chr // TRIM( mode_dvrp(m) ) // '/camera.dvr'
             dvrp_file_c = dvrp_file
             CALL DVRP_OUTPUT_RTSP( m-1, dvrp_host_c, dvrp_username_c, &
                                    dvrp_password_c, dvrp_directory_c, &
                                    dvrp_file_c )

          ELSEIF ( dvrp_output == 'ftp' )  THEN

             dvrp_file   = prefix_chr // TRIM( mode_dvrp(m) ) // '.camera.dvr'
             dvrp_file_c = dvrp_file
!             CALL DVRP_OUTPUT_FTP( m-1, 0, dvrp_host_c, dvrp_username_c, &
!                                   dvrp_password_c, dvrp_directory_c,    &
!                                   dvrp_file_c )

          ELSE

             IF ( dvrp_file(1:9) /= '/dev/null' )  THEN
                dvrp_file_local   = prefix_chr // TRIM( mode_dvrp(m) )  &
                     // '.camera.dvr'
                dvrp_file_local_c = dvrp_file_local
             ELSE
                dvrp_file_local_c = dvrp_file_c
             ENDIF
             CALL DVRP_OUTPUT_LOCAL( m-1, 0, dvrp_file_local_c )

          ENDIF

          CALL DVRP_CAMERA( m-1, center, distance )

!
!--       Define bounding box material and create a bounding box
          tmp_r = 0.5_wp;  tmp_g = 0.5_wp;  tmp_b = 0.5_wp;  tmp_t = 0.0_wp
          CALL DVRP_MATERIAL_RGB( m-1, 1, tmp_r, tmp_g, tmp_b, tmp_t )

          tmp_1 = 0.01_wp;
          tmp_2 = clip_dvrp_l * superelevation_x
          tmp_3 = clip_dvrp_s * superelevation_y
          tmp_4 = 0.0_wp
          tmp_5 = (clip_dvrp_r+dx) * superelevation_x
          tmp_6 = (clip_dvrp_n+dy) * superelevation_y
          tmp_7 = zu(nz_do3d) * superelevation
          CALL DVRP_BOUNDINGBOX( m-1, 1, tmp_1, tmp_2, tmp_3, tmp_4, tmp_5, &
                                 tmp_6, tmp_7 )

          CALL DVRP_VISUALIZE( m-1, 0, 0 )
          CALL DVRP_EXIT( m-1 )

!
!--       Write topography isosurface on file
          IF ( TRIM( topography ) /= 'flat' )  THEN

             CALL DVRP_INIT( m-1, 0 )

!
!--          Create filename for topography
             IF ( dvrp_output == 'rtsp' )  THEN

                dvrp_file   = prefix_chr // TRIM( mode_dvrp(m) )  &
                              // '/topography.dvr'
                dvrp_file_c = dvrp_file
                CALL DVRP_OUTPUT_RTSP( m-1, dvrp_host_c, dvrp_username_c, &
                                       dvrp_password_c, dvrp_directory_c, &
                                       dvrp_file_c )

             ELSEIF ( dvrp_output == 'ftp' )  THEN

                dvrp_file   = prefix_chr // TRIM( mode_dvrp(m) )  &
                              // '.topography.dvr'
                dvrp_file_c = dvrp_file
!                CALL DVRP_OUTPUT_FTP( m-1, 0, dvrp_host_c, dvrp_username_c, &
!                                      dvrp_password_c, dvrp_directory_c,    &
!                                      dvrp_file_c )

             ELSE

                IF ( dvrp_file(1:9) /= '/dev/null' )  THEN
                   dvrp_file_local   = prefix_chr // TRIM( mode_dvrp(m) )  &
                                       // '.topography.dvr'
                   dvrp_file_local_c = dvrp_file_local
                ELSE
                   dvrp_file_local_c = dvrp_file_c
                ENDIF
                CALL DVRP_OUTPUT_LOCAL( m-1, 0, dvrp_file_local_c )

             ENDIF

!
!--          Determine local gridpoint coordinates
             IF ( .NOT. allocated )  THEN
                ALLOCATE( xcoor_dvrp(nxl_dvrp:nxr_dvrp+1), &
                          ycoor_dvrp(nys_dvrp:nyn_dvrp+1), &
                          zcoor_dvrp(nzb:nz_do3d) )
                allocated = .TRUE.

                DO  i = nxl_dvrp, nxr_dvrp+1
                   xcoor_dvrp(i) = i * dx * superelevation_x
                ENDDO
                DO  j = nys_dvrp, nyn_dvrp+1
                   ycoor_dvrp(j) = j * dy * superelevation_y
                ENDDO
                zcoor_dvrp = zu(nzb:nz_do3d) * superelevation
                nx_dvrp    = nxr_dvrp+1 - nxl_dvrp + 1
                ny_dvrp    = nyn_dvrp+1 - nys_dvrp + 1
                nz_dvrp    = nz_do3d - nzb + 1
             ENDIF

!
!--          Define the grid used by dvrp
             CALL DVRP_NO_GLOBAL_GRID( m-1, 1 )
             CALL DVRP_GRID( m-1, nx_dvrp, ny_dvrp, nz_dvrp, xcoor_dvrp, &
                             ycoor_dvrp, zcoor_dvrp )

             tmp_r = topography_color(1)
             tmp_g = topography_color(2)
             tmp_b = topography_color(3)
             tmp_t = 0.0_wp
             CALL DVRP_MATERIAL_RGB( m-1, 1, tmp_r, tmp_g, tmp_b, tmp_t )

!
!--          Compute and plot isosurface in dvr-format
             ALLOCATE( local_pf(nxl_dvrp:nxr_dvrp+1,nys_dvrp:nyn_dvrp+1, &
                                nzb:nz_do3d) )
             local_pf = 0.0_wp
             IF ( dvrp_overlap )  THEN
                DO  i = nxl_dvrp, nxr_dvrp+1
                   DO  j = nys_dvrp, nyn_dvrp+1
                      IF ( nzb_s_inner(j,i) > 0 )  THEN
                         local_pf(i,j,nzb:nzb_s_inner(j,i)) = 1.0_wp
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF

             CALL DVRP_DATA( m-1, local_pf, 1, nx_dvrp, ny_dvrp, nz_dvrp, &
                             cyclic_dvrp, cyclic_dvrp, cyclic_dvrp )

             tmp_th = 1.0_wp
             CALL DVRP_THRESHOLD( m-1, tmp_th )

!
!--          Reduce the number of polygones, if required
             IF ( cluster_size > 1 )  THEN

                cluster_size_x = cluster_size
                cluster_size_y = cluster_size
                cluster_size_z = cluster_size
                cluster_mode     = 4    ! vertex clustering mode
                gradient_normals = 0    ! use flat-shading

                CALL DVRP_CLUSTER_SIZE( m-1, cluster_size_x, cluster_size_y, &
                                        cluster_size_z )
                CALL DVRP_CLUSTERING_MODE( m-1, cluster_mode )
                CALL DVRP_GRADIENTNORMALS( m-1, gradient_normals )
!
!--             Set parameter for vertex clustering mode 4.
!--             ATTENTION: A seperate procedure for setting cluster_alpha will
!--                        be in the next version of libDVRP (Feb 09)
                cluster_alpha = 38.0_wp
                CALL DVRP_THRESHOLD( -(m-1)-1, cluster_alpha )

                CALL DVRP_VISUALIZE( m-1, 21, 0 )

             ELSE
!
!--             No polygon reduction
                CALL DVRP_VISUALIZE( m-1, 1, 0 )

             ENDIF

             DEALLOCATE( local_pf )

             CALL DVRP_EXIT( m-1 )

          ENDIF

!
!--       Write the ground plate (z=0) isosurface on file
          CALL DVRP_INIT( m-1, 0 )

!
!--       Create filename for surface
          IF ( dvrp_output == 'rtsp' )  THEN

             dvrp_file   = prefix_chr // TRIM( mode_dvrp(m) ) // &
                           '/groundplate.dvr'
             dvrp_file_c = dvrp_file
             CALL DVRP_OUTPUT_RTSP( m-1, dvrp_host_c, dvrp_username_c, &
                                    dvrp_password_c, dvrp_directory_c, &
                                    dvrp_file_c )

          ELSEIF ( dvrp_output == 'ftp' )  THEN

             dvrp_file   = prefix_chr // TRIM( mode_dvrp(m) ) // &
                           '.groundplate.dvr'
             dvrp_file_c = dvrp_file
!             CALL DVRP_OUTPUT_FTP( m-1, 0, dvrp_host_c, dvrp_username_c, &
!                                   dvrp_password_c, dvrp_directory_c,    &
!                                   dvrp_file_c )

          ELSE

             IF ( dvrp_file(1:9) /= '/dev/null' )  THEN
                dvrp_file_local   = prefix_chr // TRIM( mode_dvrp(m) )  &
                     // '.groundplate.dvr'
                dvrp_file_local_c = dvrp_file_local
             ELSE
                dvrp_file_local_c = dvrp_file_c
             ENDIF
             CALL DVRP_OUTPUT_LOCAL( m-1, 0, dvrp_file_local_c )

          ENDIF

!
!--       Determine local gridpoint coordinates
          IF ( .NOT. allocated )  THEN
             ALLOCATE( xcoor_dvrp(nxl_dvrp:nxr_dvrp+1), &
                       ycoor_dvrp(nys_dvrp:nyn_dvrp+1), &
                       zcoor_dvrp(nzb:nz_do3d) )
             allocated = .TRUE.

             DO  i = nxl_dvrp, nxr_dvrp+1
                xcoor_dvrp(i) = i * dx * superelevation_x
             ENDDO
             DO  j = nys_dvrp, nyn_dvrp+1
                ycoor_dvrp(j) = j * dy * superelevation_y
             ENDDO
             zcoor_dvrp = zu(nzb:nz_do3d) * superelevation
             nx_dvrp    = nxr_dvrp+1 - nxl_dvrp + 1
             ny_dvrp    = nyn_dvrp+1 - nys_dvrp + 1
             nz_dvrp    = nz_do3d - nzb + 1
          ENDIF

!
!--       Define the grid used by dvrp
          CALL DVRP_NO_GLOBAL_GRID( m-1, 1 )
          CALL DVRP_GRID( m-1, nx_dvrp, ny_dvrp, nz_dvrp, xcoor_dvrp, &
                          ycoor_dvrp, zcoor_dvrp )

          tmp_r = groundplate_color(1)
          tmp_g = groundplate_color(2)
          tmp_b = groundplate_color(3)
          tmp_t = 0.0_wp
          CALL DVRP_MATERIAL_RGB( m-1, 1, tmp_r, tmp_g, tmp_b, tmp_t )

!
!--       Compute and plot isosurface in dvr-format
          ALLOCATE( local_pf(nxl_dvrp:nxr_dvrp+1,nys_dvrp:nyn_dvrp+1, &
                             nzb:nz_do3d) )
          local_pf = 0.0_wp
          IF (dvrp_overlap )  local_pf(:,:,0) = 1.0_wp

          CALL DVRP_DATA( m-1, local_pf, 1, nx_dvrp, ny_dvrp, nz_dvrp, &
                          cyclic_dvrp, cyclic_dvrp, cyclic_dvrp )
          tmp_th = 1.0_wp
          CALL DVRP_THRESHOLD( m-1, tmp_th )

!
!--       Always reduce the number of polygones as much as possible
          cluster_size_x = 5
          cluster_size_y = 5
          cluster_size_z = 5
          cluster_mode     = 4    ! vertex clustering mode
          gradient_normals = 0    ! use flat-shading

          CALL DVRP_CLUSTER_SIZE( m-1, cluster_size_x, cluster_size_y, &
                                  cluster_size_z )
          CALL DVRP_CLUSTERING_MODE( m-1, cluster_mode )
          CALL DVRP_GRADIENTNORMALS( m-1, gradient_normals )
!
!--       Set parameter for vertex clustering mode 4.
!--       ATTENTION: A seperate procedure for setting cluster_alpha will be in
!--                  the next version of libDVRP (Feb 09)
          cluster_alpha = 38.0_wp
          CALL DVRP_THRESHOLD( -(m-1)-1, cluster_alpha )

          CALL DVRP_VISUALIZE( m-1, 21, 0 )

          DEALLOCATE( local_pf )

          CALL DVRP_EXIT( m-1 )
    
!       ENDIF


!
!--    Initialize dvrp for all dvrp-calls during the run
       CALL DVRP_INIT( m-1, 0 )

!
!--    Preliminary definition of filename for dvrp-output
       IF ( dvrp_output == 'rtsp' )  THEN

!
!--       First initialize parameters for possible interactive steering.
!--       Every parameter has to be passed to the respective stream.
          pn = 1
!
!--       Initialize threshold counter needed for initialization of the
!--       isosurface steering variables
          tv = 0

          DO WHILE ( mode_dvrp(pn) /= ' ' )

             IF ( mode_dvrp(pn)(1:10) == 'isosurface' )  THEN

                READ ( mode_dvrp(pn), '(10X,I2)' )  vn
                steering_dvrp(pn)%name = do3d(0,vn)
                tv = tv + 1

                IF ( do3d(0,vn)(1:1) == 'w' )  THEN
                   steering_dvrp(pn)%min  = -4.0_wp
                   steering_dvrp(pn)%max  =  5.0_wp
                ELSE
                   steering_dvrp(pn)%min  = 288.0_wp
                   steering_dvrp(pn)%max  = 292.0_wp
                ENDIF

                name_c  = TRIM( do3d(0,vn) )
                tmp_thr = threshold(tv)
                CALL DVRP_STEERING_INIT( m-1, name_c, steering_dvrp(pn)%min, &
                                         steering_dvrp(pn)%max, tmp_thr )

             ELSEIF ( mode_dvrp(pn)(1:6) == 'slicer' )  THEN

                READ ( mode_dvrp(pn), '(6X,I2)' )  vn
                steering_dvrp(pn)%name = do2d(0,vn)
                name_c = TRIM( do2d(0,vn) )

                l = MAX( 2, LEN_TRIM( do2d(0,vn) ) )
                section_chr = do2d(0,vn)(l-1:l)
                SELECT CASE ( section_chr )
                   CASE ( 'xy' )
                      steering_dvrp(pn)%imin   = 0
                      steering_dvrp(pn)%imax   = nz_do3d
                      slicer_position_dvrp(pn) = section(1,1)
                      CALL DVRP_STEERING_INIT( m-1, name_c,            &
                                               steering_dvrp(pn)%imin, &
                                               steering_dvrp(pn)%imax, &
                                               slicer_position_dvrp(pn) )
                   CASE ( 'xz' )
                      steering_dvrp(pn)%imin   = 0
                      steering_dvrp(pn)%imax   = ny
                      slicer_position_dvrp(pn) = section(1,2)
                      CALL DVRP_STEERING_INIT( m-1, name_c,            &
                                               steering_dvrp(pn)%imin, &
                                               steering_dvrp(pn)%imax, &
                                               slicer_position_dvrp(pn) )
                   CASE ( 'yz' )
                      steering_dvrp(pn)%imin = 0
                      steering_dvrp(pn)%imax = nx
                      slicer_position_dvrp(pn) = section(1,3)
                      CALL DVRP_STEERING_INIT( m-1, name_c,            &
                                               steering_dvrp(pn)%imin, &
                                               steering_dvrp(pn)%imax, &
                                               slicer_position_dvrp(pn) )
                END SELECT

             ENDIF

             pn = pn + 1

          ENDDO

          dvrp_file = prefix_chr // TRIM( mode_dvrp(m) ) // '/*****.dvr'
          dvrp_file_c = dvrp_file
          CALL DVRP_OUTPUT_RTSP( m-1, dvrp_host_c, dvrp_username_c, &
                                 dvrp_password_c, dvrp_directory_c, &
                                 dvrp_file_c )

       ELSEIF ( dvrp_output == 'ftp' )  THEN

          dvrp_file   = prefix_chr // TRIM( mode_dvrp(m) ) // '.%05d.dvr'
          dvrp_file_c = dvrp_file
!          CALL DVRP_OUTPUT_FTP( m-1, 0, dvrp_host_c, dvrp_username_c, &
!                                dvrp_password_c, dvrp_directory_c, dvrp_file_c )

       ELSE

          IF ( dvrp_file(1:9) /= '/dev/null' )  THEN
             dvrp_file_local   = prefix_chr // TRIM( mode_dvrp(m) )  &
                  // '_%05d.dvr'
             dvrp_file_local_c = dvrp_file_local
          ELSE
             dvrp_file_local_c = dvrp_file_c
          ENDIF
          CALL DVRP_OUTPUT_LOCAL( m-1, 0, dvrp_file_local_c )

       ENDIF

!
!--    Determine local gridpoint coordinates
       IF ( .NOT. allocated )  THEN
          ALLOCATE( xcoor_dvrp(nxl_dvrp:nxr_dvrp+1), &
                    ycoor_dvrp(nys_dvrp:nyn_dvrp+1), &
                    zcoor_dvrp(nzb:nz_do3d) )
          allocated = .TRUE.

          DO  i = nxl_dvrp, nxr_dvrp+1
             xcoor_dvrp(i) = i * dx * superelevation_x
          ENDDO
          DO  j = nys_dvrp, nyn_dvrp+1
             ycoor_dvrp(j) = j * dy * superelevation_y
          ENDDO
          zcoor_dvrp = zu(nzb:nz_do3d) * superelevation
          nx_dvrp    = nxr_dvrp+1 - nxl_dvrp + 1
          ny_dvrp    = nyn_dvrp+1 - nys_dvrp + 1
          nz_dvrp    = nz_do3d - nzb + 1
       ENDIF

!
!--    Define the grid used by dvrp
       IF ( mode_dvrp(m) /= 'pathlines' )  THEN
          CALL DVRP_NO_GLOBAL_GRID( m-1, 1 )
       ENDIF
       CALL DVRP_GRID( m-1, nx_dvrp, ny_dvrp, nz_dvrp, xcoor_dvrp, ycoor_dvrp, &
                       zcoor_dvrp )

       IF ( mode_dvrp(m) == 'pathlines' )  THEN

          tmp_x1 = 0.0_wp;  tmp_y1 = 0.0_wp;  tmp_z1 = 0.0_wp
          tmp_x2 = 1.0_wp;  tmp_y2 = 1.0_wp;  tmp_z2 = 0.3_wp
          CALL DVRP_CUBIC_SEEDING( m-1, tmp_x1, tmp_y1, tmp_z1, tmp_x2, tmp_y2,&
                                   tmp_z2, pathlines_linecount, 2, 0 )
!
!--       Set wavecount and wavetime
          CALL DVRP_PATHLINES_BEHAVIOUR_WAVE( m-1, pathlines_wavecount, &
                                              pathlines_wavetime,       &
                                              pathlines_fadeintime,     &
                                              pathlines_fadeouttime )
!
!--       Set pathline length
          CALL DVRP_PATHLINES_SETMAXHISTORY( m-1, pathlines_maxhistory )
          CALL DVRP_PATHLINES_SETFADING( m-1, 1, 0.0_wp )

          CALL DVRP_INIT_PATHLINES( m-1, 0 )

       ENDIF

       IF ( mode_dvrp(m)(1:9) == 'particles' )  THEN
!
!--       Define a default colourtable for particles
          DO  i = 1, 11
             interval_values_dvrp_prt(1,i) = i - 1.0_wp
             interval_values_dvrp_prt(2,i) = REAL( i, KIND=wp )
             interval_h_dvrp_prt(:,i) = 270.0_wp - ( i - 1.0_wp ) * 9.0_wp
          ENDDO

          DO  i = 12, 22
             interval_values_dvrp_prt(1,i) = i - 1.0_wp
             interval_values_dvrp_prt(2,i) = REAL( i, KIND=wp )
             interval_h_dvrp_prt(:,i) = 70.0_wp - ( i - 12.0_wp ) * 9.5_wp
          ENDDO

          dvrp_colortable_entries_prt = 22

       ENDIF

       m = m + 1

    ENDDO

#endif
 END SUBROUTINE init_dvrp

 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializes logging events for time measurement with dvrp software
!> and splits one PE from the global communicator in case that dvrp output
!> shall be done by one single PE.
!------------------------------------------------------------------------------!
 
 SUBROUTINE init_dvrp_logging

#if defined( __dvrp_graphics )

    USE dvrp_variables,                                                        &
        ONLY:  use_seperate_pe_for_dvrp_output
    
    USE kinds
    
    USE pegrid

    IMPLICIT NONE

    CHARACTER (LEN=4) ::  chr  !<
    
    INTEGER(iwp) ::  idummy    !<

!
!-- Initialize logging of calls by DVRP graphic software
    CALL DVRP_LOG_INIT( 'DVRP_LOG' // CHAR( 0 ), 0 )

!
!-- User-defined logging events: #1 (total time needed by PALM)
    CALL DVRP_LOG_SYMBOL( 1, 'PALM_total' // CHAR( 0 ) )
    CALL DVRP_LOG_SYMBOL( 2, 'PALM_timestep' // CHAR( 0 ) )
    CALL DVRP_LOG_EVENT( 1, 1 )

#if defined( __parallel )
!
!-- Find out, if dvrp output shall be done by a dedicated PE
    CALL GET_ENVIRONMENT_VARIABLE( 'use_seperate_pe_for_dvrp_output', chr,     &
                                   idummy )
    IF ( chr == 'true' )  THEN

       use_seperate_pe_for_dvrp_output = .TRUE.

!
!--    Adjustment for new MPI-1 coupling. This might be unnecessary.
       IF ( coupling_mode /= 'uncoupled' ) THEN
          message_string = 'split of communicator not realized with' // &
                          ' MPI1 coupling atmosphere-ocean'
          CALL message( 'init_dvrp_logging', 'PA0199', 1, 2, 0, 6, 0 )
 
          CALL DVRP_SPLIT( comm_inter, comm_palm )
       ELSE
          CALL DVRP_SPLIT( MPI_COMM_WORLD, comm_palm )
       ENDIF

       CALL MPI_COMM_SIZE( comm_palm, numprocs, ierr )

    ENDIF
#endif

#endif
 END SUBROUTINE init_dvrp_logging


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Exit of dvrp software and finish dvrp logging
!------------------------------------------------------------------------------!
 
 SUBROUTINE close_dvrp

#if defined( __dvrp_graphics )
                                               
    USE DVRP
    
    USE dvrp_variables,                                                        &
        ONLY: use_seperate_pe_for_dvrp_output
    
    USE kinds

    INTEGER(iwp) ::  m  !<

!
!-- If required, close dvrp-software and logging of dvrp-calls
    IF ( dt_dvrp /= 9999999.9_wp )  THEN
       m = 1
       DO WHILE ( mode_dvrp(m) /= ' ' )
          CALL DVRP_EXIT( m-1 )
          m = m + 1
       ENDDO
       CALL DVRP_LOG_EVENT( -1, 1 )   ! Logging of total cpu-time used by PALM
       IF ( use_seperate_pe_for_dvrp_output )  THEN
          CALL DVRP_SPLIT_EXIT( 1 )      ! Argument 0: reduced output
       ELSE
          CALL DVRP_LOG_EXIT( 1 )        ! Argument 0: reduced output
       ENDIF
    ENDIF

#endif
 END SUBROUTINE close_dvrp
