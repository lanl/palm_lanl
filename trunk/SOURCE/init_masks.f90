!> @file init_masks.f90
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
! $Id: init_masks.f90 3065 2018-06-12 07:03:02Z Giersch $
! dz_stretch_level was replaced by dz_stretch_level_start
! 
! 3049 2018-05-29 13:52:36Z Giersch
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2344 2017-08-09 11:00:34Z raasch
! explicit setting of initial values for array domask required due to a bug
! in the Cray compiler (appears only if option -eD is used)
! 
! 2301 2017-06-29 16:41:21Z gronemeier
! Bugfix: set variable name length to global value
! 
! 
! Implementation of new microphysic scheme: cloud_scheme = 'morrison' 
! includes two more prognostic equations for cloud drop concentration (nc)  
! and cloud water content (qc). 
! 
! 2271 2017-06-09 12:34:55Z sward
! Changed error message
! 
! 2101 2017-01-05 16:42:31Z suehring
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1822 2016-04-07 07:49:42Z hoffmann
! icloud_scheme replaced by microphysics_*
!
! 1804 2016-04-05 16:30:18Z maronga
! Removed code for parameter file check (__check)
!
! 1783 2016-03-06 18:36:17Z raasch
! netcdf module name changed + related changes
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1438 2014-07-22 14:14:06Z heinze
! +nr, qc, qr
! 
! 1414 2014-05-31 11:19:48Z gronemeier
! Bugfix: first and last grid points as they appear in 3D volume data can now
! be included to masked data output
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1324 2014-03-21 09:13:16Z suehring
! Bugfix: ONLY statement for module netcdf_control removed 
!
! 1320 2014-03-20 08:40:49Z raasch 
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1186 2013-06-18 06:22:52Z raasch
! bugfix: 0.0 replaced by zu(nzb) as the lowest default height level for masks,
! because a zero value does not work in case of ocean runs
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1032 2012-10-21 13:03:21Z letzel
! mask locations determined based on scalar positions
!
! 1031 2012-10-19 14:35:30Z raasch
! netCDF4 without parallel file support implemented
!
! 996 2012-09-07 10:41:47Z raasch
! little reformatting
!
! 978 2012-08-09 08:28:32Z fricke
! +z0h*
!
! 809 2012-01-30 13:32:58Z maronga
! Bugfix: replaced .AND. and .NOT. with && and ! in the preprocessor directives 
!
! 807 2012-01-25 11:53:51Z maronga
! New cpp directive "__check" implemented which is used by check_namelist_files
!
! 410 2009-12-04 17:05:40Z letzel
! Initial revision
!
!
! Description:
! ------------
!> Initialize masked data output
!------------------------------------------------------------------------------!
 SUBROUTINE init_masks

    USE arrays_3d,                                                             &
        ONLY:  zu, zw

    USE control_parameters,                                                    &
        ONLY:  constant_diffusion, cloud_droplets, cloud_physics,              &
               data_output_masks, data_output_masks_user,                      &
               doav, doav_n, domask, domask_no, dz, dz_stretch_level_start,    &
               humidity, mask, masks, mask_scale, mask_i,                      &
               mask_i_global, mask_j, mask_j_global, mask_k, mask_k_global,    &
               mask_loop, mask_size, mask_size_l, mask_start_l, mask_x,        &
               mask_x_loop, mask_xyz_dimension, mask_y, mask_y_loop, mask_z,   &
               mask_z_loop, max_masks,  message_string, mid,                   &
               microphysics_morrison, microphysics_seifert, passive_scalar,    &
               ocean, varnamelength
               

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE indices,                                                               &
        ONLY:  nx, nxl, nxr, ny, nyn, nys, nz, nzb, nzt

    USE kinds

    USE netcdf_interface,                                                      &
        ONLY:  domask_unit, netcdf_data_format

    USE pegrid

    IMPLICIT NONE

    CHARACTER (LEN=varnamelength) ::  var  !<
    CHARACTER (LEN=7) ::  unit !<
    
    CHARACTER (LEN=varnamelength), DIMENSION(max_masks,100) ::  do_mask      !<
    CHARACTER (LEN=varnamelength), DIMENSION(max_masks,100) ::  do_mask_user !<

    INTEGER(iwp) ::  i            !< 
    INTEGER(iwp) ::  ilen         !<
    INTEGER(iwp) ::  ind(6)       !<
    INTEGER(iwp) ::  ind_array(1) !<
    INTEGER(iwp) ::  j            !<
    INTEGER(iwp) ::  k            !<
    INTEGER(iwp) ::  n            !<
    INTEGER(iwp) ::  sender       !<
    
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  tmp_array !<

    LOGICAL ::  found !<

!
!-- Initial values are explicitly set here due to a bug in the Cray compiler
!-- in case of assignments of initial values in declaration statements for
!-- arrays with more than 9999 elements (appears with -eD only)
    domask = ' '

!
!-- Allocation and initialization
    ALLOCATE( tmp_array( MAX(nx,ny,nz)+2 ) )

    ALLOCATE( mask_i(max_masks,nxr-nxl+2), &
              mask_j(max_masks,nyn-nys+2), &
              mask_k(max_masks,nzt-nzb+2) )
!
!-- internal mask arrays ("mask,dimension,selection")
    ALLOCATE( mask(max_masks,3,mask_xyz_dimension), &
              mask_loop(max_masks,3,3) )
    
!
!-- Parallel mask output not yet supported. In check_parameters data format 
!-- is restricted and is switched back to non-parallel output. Therefore the
!-- following error can not occur at the moment.
    IF ( netcdf_data_format > 4 )  THEN
       message_string = 'netCDF file formats '//                               &
                        '5 and 6 (with parallel I/O support)'//                &
                        ' are currently not supported.'
       CALL message( 'init_masks', 'PA0328', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Store data output parameters for masked data output in few shared arrays
    DO  mid = 1, masks
    
       do_mask     (mid,:) = data_output_masks(mid,:)
       do_mask_user(mid,:) = data_output_masks_user(mid,:)
       mask      (mid,1,:) = mask_x(mid,:) 
       mask      (mid,2,:) = mask_y(mid,:)
       mask      (mid,3,:) = mask_z(mid,:) 
       
       IF ( mask_x_loop(mid,1) == -1.0_wp  .AND.  mask_x_loop(mid,2) == -1.0_wp&
            .AND.  mask_x_loop(mid,3) == -1.0_wp )  THEN
          mask_loop(mid,1,1:2) = -1.0_wp
          mask_loop(mid,1,3)   =  0.0_wp
       ELSE
          mask_loop(mid,1,:) = mask_x_loop(mid,:)
       ENDIF
       IF ( mask_y_loop(mid,1) == -1.0_wp  .AND.  mask_y_loop(mid,2) == -1.0_wp&
            .AND.  mask_y_loop(mid,3) == -1.0_wp )  THEN
          mask_loop(mid,2,1:2) = -1.0_wp
          mask_loop(mid,2,3)   =  0.0_wp
       ELSE
          mask_loop(mid,2,:) = mask_y_loop(mid,:)
       ENDIF
       IF ( mask_z_loop(mid,1) == -1.0_wp  .AND.  mask_z_loop(mid,2) == -1.0_wp&
            .AND.  mask_z_loop(mid,3) == -1.0_wp )  THEN
          mask_loop(mid,3,1:2) = -1.0_wp
          mask_loop(mid,3,3)   =  0.0_wp
       ELSE
          mask_loop(mid,3,:) = mask_z_loop(mid,:)
       ENDIF
       
    ENDDO
    
    mask_i = -1; mask_j = -1; mask_k = -1
    
!
!-- Global arrays are required by define_netcdf_header.
    IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
       ALLOCATE( mask_i_global(max_masks,nx+2), &
                 mask_j_global(max_masks,ny+2), &
                 mask_k_global(max_masks,nz+2) )
       mask_i_global = -1; mask_j_global = -1; mask_k_global = -1
    ENDIF

!
!-- Determine variable names for each mask
    DO  mid = 1, masks
!
!--    Append user-defined data output variables to the standard data output
       IF ( do_mask_user(mid,1) /= ' ' )  THEN
          i = 1
          DO  WHILE ( do_mask(mid,i) /= ' '  .AND.  i <= 100 )
             i = i + 1
          ENDDO
          j = 1
          DO  WHILE ( do_mask_user(mid,j) /= ' '  .AND.  j <= 100 )
             IF ( i > 100 )  THEN
                WRITE ( message_string, * ) 'number of output quantitities ',  &
                     'given by data_output_mask and data_output_mask_user ',   &
                     'exceeds the limit of 100'
                CALL message( 'init_masks', 'PA0329', 1, 2, 0, 6, 0 )
             ENDIF
             do_mask(mid,i) = do_mask_user(mid,j)
             i = i + 1
             j = j + 1
          ENDDO
       ENDIF

!
!--    Check and set steering parameters for mask data output and averaging
       i   = 1
       DO WHILE ( do_mask(mid,i) /= ' '  .AND.  i <= 100 )
!
!--       Check for data averaging
          ilen = LEN_TRIM( do_mask(mid,i) )
          j = 0                                              ! no data averaging
          IF ( ilen > 3 )  THEN
             IF ( do_mask(mid,i)(ilen-2:ilen) == '_av' )  THEN
                j = 1                                           ! data averaging
                do_mask(mid,i) = do_mask(mid,i)(1:ilen-3)
             ENDIF
          ENDIF
          var = TRIM( do_mask(mid,i) )
!
!--       Check for allowed value and set units
          SELECT CASE ( TRIM( var ) )

             CASE ( 'e' )
                IF ( constant_diffusion )  THEN
                   WRITE ( message_string, * ) 'output of "', TRIM( var ),     &
                        '" requires constant_diffusion = .FALSE.'
                   CALL message( 'init_masks', 'PA0103', 1, 2, 0, 6, 0 )
                ENDIF
                unit = 'm2/s2'

             CASE ( 'rho_ocean' )
                IF ( .NOT. ocean )  THEN
                   WRITE ( message_string, * ) 'output of "', TRIM( var ),     &
                        '" requires ocean = .TRUE.'
                   CALL message( 'init_masks', 'PA0109', 1, 2, 0, 6, 0 )
                ENDIF
                unit = 'kg/m3'

             CASE ( 's' )
                IF ( .NOT. passive_scalar )  THEN
                   WRITE ( message_string, * ) 'output of "', TRIM( var ),     &
                        '" requires passive_scalar = .TRUE.'
                   CALL message( 'init_masks', 'PA0110', 1, 2, 0, 6, 0 )
                ENDIF
                unit = 'conc'

             CASE ( 'sa' )
                IF ( .NOT. ocean )  THEN
                   WRITE ( message_string, * ) 'output of "', TRIM( var ),     &
                        '" requires ocean = .TRUE.'
                   CALL message( 'init_masks', 'PA0109', 1, 2, 0, 6, 0 )
                ENDIF
                unit = 'psu'

             CASE ( 'u*', 't*', 'lwp*', 'pra*', 'prr*', 'z0*', 'z0h*' )
                WRITE ( message_string, * ) 'illegal value for data_',         &
                     'output: "', TRIM( var ), '" is only allowed',            &
                     'for horizontal cross section'
                CALL message( 'init_masks', 'PA0111', 1, 2, 0, 6, 0 )

             CASE ( 'p', 'pt', 'u', 'v', 'w' )
                IF ( TRIM( var ) == 'p'  )  unit = 'Pa'
                IF ( TRIM( var ) == 'pt' )  unit = 'K'
                IF ( TRIM( var ) == 'u'  )  unit = 'm/s'
                IF ( TRIM( var ) == 'v'  )  unit = 'm/s'
                IF ( TRIM( var ) == 'w'  )  unit = 'm/s'
                CONTINUE

             CASE DEFAULT
                IF ( unit == 'illegal' )  THEN
                   IF ( do_mask_user(mid,1) /= ' ' )  THEN
                      WRITE ( message_string, * ) 'illegal value for data_',   &
                           'output_masks or data_output_masks_user: "',        &
                           TRIM( do_mask(mid,i) ), '"'
                      CALL message( 'init_masks', 'PA0018', 1, 2, 0, 6, 0 )
                   ELSE
                      WRITE ( message_string, * ) 'illegal value for data_',   &
                           ' output_masks : "', TRIM( do_mask(mid,i) ), '"'
                      CALL message( 'init_masks', 'PA0330', 1, 2, 0, 6, 0 )
                   ENDIF
                ENDIF

          END SELECT
!
!--       Set the internal steering parameters appropriately
          domask_no(mid,j)                    = domask_no(mid,j) + 1
          domask(mid,j,domask_no(mid,j))      = do_mask(mid,i)
          domask_unit(mid,j,domask_no(mid,j)) = unit

          IF ( j == 1 )  THEN
!
!--          Check, if variable is already subject to averaging
             found = .FALSE.
             DO  k = 1, doav_n
                IF ( TRIM( doav(k) ) == TRIM( var ) )  found = .TRUE.
             ENDDO

             IF ( .NOT. found )  THEN
                doav_n = doav_n + 1
                doav(doav_n) = var
             ENDIF
          ENDIF

          i = i + 1

       ENDDO   ! do_mask(mid,i)
    ENDDO   ! mid


!
!-- Determine mask locations for each mask
    DO  mid = 1, masks
!
!--    Set local masks for each subdomain along all three dimensions
       CALL set_mask_locations( 1, dx, 'dx', nx, 'nx', nxl, nxr )
       CALL set_mask_locations( 2, dy, 'dy', ny, 'ny', nys, nyn )
       CALL set_mask_locations( 3, dz(1), 'dz', nz, 'nz', nzb, nzt )
!
!--    Set global masks along all three dimensions (required by
!--    define_netcdf_header).
#if defined( __parallel )
!
!--    PE0 receives partial arrays from all processors of the respective mask
!--    and outputs them. Here a barrier has to be set, because otherwise 
!--    "-MPI- FATAL: Remote protocol queue full" may occur.

       CALL MPI_BARRIER( comm2d, ierr )

       IF ( myid == 0 )  THEN
!
!--       Local arrays can be relocated directly.
          mask_i_global(mid,mask_start_l(mid,1): &
                       mask_start_l(mid,1)+mask_size_l(mid,1)-1) = &
                       mask_i(mid,:mask_size_l(mid,1))
          mask_j_global(mid,mask_start_l(mid,2): &
                       mask_start_l(mid,2)+mask_size_l(mid,2)-1) = &
                       mask_j(mid,:mask_size_l(mid,2))
          mask_k_global(mid,mask_start_l(mid,3): &
                       mask_start_l(mid,3)+mask_size_l(mid,3)-1) = &
                       mask_k(mid,:mask_size_l(mid,3))
!
!--       Receive data from all other PEs.
          DO  n = 1, numprocs-1
!
!--          Receive index limits first, then arrays.
!--          Index limits are received in arbitrary order from the PEs.
             CALL MPI_RECV( ind(1), 6, MPI_INTEGER, MPI_ANY_SOURCE, 0,  &
                  comm2d, status, ierr )
!
!--          Not all PEs have data for the mask.
             IF ( ind(1) /= -9999 )  THEN
                sender = status(MPI_SOURCE)
                CALL MPI_RECV( tmp_array(ind(1)), ind(2)-ind(1)+1,  &
                               MPI_INTEGER, sender, 1, comm2d, status, ierr )
                mask_i_global(mid,ind(1):ind(2)) = tmp_array(ind(1):ind(2))
                CALL MPI_RECV( tmp_array(ind(3)), ind(4)-ind(3)+1,  &
                               MPI_INTEGER, sender, 2, comm2d, status, ierr )
                mask_j_global(mid,ind(3):ind(4)) = tmp_array(ind(3):ind(4))
                CALL MPI_RECV( tmp_array(ind(5)), ind(6)-ind(5)+1,  &
                               MPI_INTEGER, sender, 3, comm2d, status, ierr )
                mask_k_global(mid,ind(5):ind(6)) = tmp_array(ind(5):ind(6))
             ENDIF
          ENDDO

       ELSE
!
!--       If at least part of the mask resides on the PE, send the index limits
!--       for the target array, otherwise send -9999 to PE0.
          IF ( mask_size_l(mid,1) > 0  .AND.  mask_size_l(mid,2) > 0  .AND.  &
               mask_size_l(mid,3) > 0  )  THEN
             ind(1) = mask_start_l(mid,1)
             ind(2) = mask_start_l(mid,1) + mask_size_l(mid,1) - 1
             ind(3) = mask_start_l(mid,2)
             ind(4) = mask_start_l(mid,2) + mask_size_l(mid,2) - 1
             ind(5) = mask_start_l(mid,3)
             ind(6) = mask_start_l(mid,3) + mask_size_l(mid,3) - 1
          ELSE
             ind(1) = -9999; ind(2) = -9999
             ind(3) = -9999; ind(4) = -9999
             ind(5) = -9999; ind(6) = -9999
          ENDIF
          CALL MPI_SEND( ind(1), 6, MPI_INTEGER, 0, 0, comm2d, ierr )
!
!--       If applicable, send data to PE0.
          IF ( ind(1) /= -9999 )  THEN
             tmp_array(:mask_size_l(mid,1)) = mask_i(mid,:mask_size_l(mid,1))
             CALL MPI_SEND( tmp_array(1), mask_size_l(mid,1),  &
                            MPI_INTEGER, 0, 1, comm2d, ierr )
             tmp_array(:mask_size_l(mid,2)) = mask_j(mid,:mask_size_l(mid,2))
             CALL MPI_SEND( tmp_array(1), mask_size_l(mid,2),  &
                            MPI_INTEGER, 0, 2, comm2d, ierr )
             tmp_array(:mask_size_l(mid,3)) = mask_k(mid,:mask_size_l(mid,3))
             CALL MPI_SEND( tmp_array(1), mask_size_l(mid,3),  &
                            MPI_INTEGER, 0, 3, comm2d, ierr )
          ENDIF
       ENDIF
!
!--    A barrier has to be set, because otherwise some PEs may proceed too fast
!--    so that PE0 may receive wrong data on tag 0.
       CALL MPI_BARRIER( comm2d, ierr )
       
       IF ( netcdf_data_format > 4 )  THEN
         
          CALL MPI_BCAST( mask_i_global(mid,:), nx+2, MPI_INTEGER, 0, comm2d, &
                          ierr )
          CALL MPI_BCAST( mask_j_global(mid,:), ny+2, MPI_INTEGER, 0, comm2d, &
                          ierr )
          CALL MPI_BCAST( mask_k_global(mid,:), nz+2, MPI_INTEGER, 0, comm2d, &
                          ierr ) 
     
       ENDIF

#else
!
!--    Local arrays can be relocated directly.
       mask_i_global(mid,:) = mask_i(mid,:)
       mask_j_global(mid,:) = mask_j(mid,:)
       mask_k_global(mid,:) = mask_k(mid,:)
#endif
    ENDDO   ! mid

    DEALLOCATE( tmp_array )
!
!-- Internal mask arrays cannot be deallocated on PE 0 because they are
!-- required for header output on PE 0.
    IF ( myid /= 0 )  DEALLOCATE( mask, mask_loop )

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Set local mask for each subdomain along 'dim' direction.
!------------------------------------------------------------------------------!
    SUBROUTINE set_mask_locations( dim, dxyz, dxyz_string, nxyz, nxyz_string, &
                                   lb, ub )

       IMPLICIT NONE

       CHARACTER (LEN=2) ::  dxyz_string !<
       CHARACTER (LEN=2) ::  nxyz_string !<
       
       INTEGER(iwp)  ::  count       !<
       INTEGER(iwp)  ::  count_l     !<
       INTEGER(iwp)  ::  dim         !<
       INTEGER(iwp)  ::  m           !<
       INTEGER(iwp)  ::  loop_begin  !<
       INTEGER(iwp)  ::  loop_end    !<
       INTEGER(iwp)  ::  loop_stride !<
       INTEGER(iwp)  ::  lb          !<
       INTEGER(iwp)  ::  nxyz        !<
       INTEGER(iwp)  ::  ub          !<
       
       REAL(wp)      ::  dxyz  !<
       REAL(wp)      ::  ddxyz !<
       REAL(wp)      ::  tmp1  !<
       REAL(wp)      ::  tmp2  !<

       count = 0;  count_l = 0  
       ddxyz = 1.0_wp / dxyz  
       tmp1  = 0.0_wp
       tmp2  = 0.0_wp

       IF ( mask(mid,dim,1) >= 0.0_wp )  THEN
!
!--       use predefined mask_* array
          DO  WHILE ( mask(mid,dim,count+1) >= 0.0_wp )
             count = count + 1
             IF ( dim == 1 .OR. dim == 2 )  THEN
                m = NINT( mask(mid,dim,count) * mask_scale(dim) * ddxyz - 0.5_wp )
                IF ( m < 0 )  m = 0  ! avoid negative values
             ELSEIF ( dim == 3 )  THEN
                ind_array =  &
                     MINLOC( ABS( mask(mid,dim,count) * mask_scale(dim) - zu ) )
                m = ind_array(1) - 1 + nzb  ! MINLOC uses lower array bound 1
             ENDIF
             IF ( m > (nxyz+1) )  THEN
                WRITE ( message_string, '(I3,A,I3,A,I1,3A,I3)' )               &
                     m,' in mask ',mid,' along dimension ',dim,                &
                     ' exceeds (',nxyz_string,'+1) = ',nxyz+1
                CALL message( 'init_masks', 'PA0331', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( ( m >= lb .AND. m <= ub ) .OR.     &
                  ( m == (nxyz+1) .AND. ub == nxyz )  )  THEN
                IF ( count_l == 0 )  mask_start_l(mid,dim) = count
                count_l = count_l + 1
                IF ( dim == 1 )  THEN
                   mask_i(mid,count_l) = m
                ELSEIF ( dim == 2 )  THEN
                   mask_j(mid,count_l) = m
                ELSEIF ( dim == 3 )  THEN
                   mask_k(mid,count_l) = m
                ENDIF
             ENDIF
             IF ( count == mask_xyz_dimension )  EXIT
          ENDDO
          mask_size(mid,dim)   = count
          mask_size_l(mid,dim) = count_l

       ELSE
!
!--       use predefined mask_loop_* array, or use the default (all grid points
!--       along this direction)
          IF ( mask_loop(mid,dim,1) < 0.0_wp )  THEN
             tmp1 = mask_loop(mid,dim,1)
             mask_loop(mid,dim,1) = zw(nzb)  !   lowest level  (default)
          ENDIF
          IF ( dim == 1 .OR. dim == 2 )  THEN
             IF ( mask_loop(mid,dim,2) < 0.0_wp )  THEN
                tmp2 = mask_loop(mid,dim,2)
                mask_loop(mid,dim,2) = (nxyz+1)*dxyz / mask_scale(dim)   ! (default)
             ENDIF
             IF ( MAXVAL( mask_loop(mid,dim,1:2) )  &
                  > (nxyz+1) * dxyz / mask_scale(dim) )  THEN
                WRITE ( message_string, '(2(A,I3,A,I1,A,F9.3),5A,I1,A,F9.3)' ) &
                     'mask_loop(',mid,',',dim,',1)=',mask_loop(mid,dim,1),     &
                     ' and/or mask_loop(',mid,',',dim,',2)=', &
                     mask_loop(mid,dim,2),' exceed (', &
                     nxyz_string,'+1)*',dxyz_string,'/mask_scale(',dim,')=',   &
                     (nxyz+1)*dxyz/mask_scale(dim)
                CALL message( 'init_masks', 'PA0332', 1, 2, 0, 6, 0 )
             ENDIF
             loop_begin  = NINT( mask_loop(mid,dim,1) * mask_scale(dim)        &
                  * ddxyz - 0.5_wp )
             loop_end    = NINT( mask_loop(mid,dim,2) * mask_scale(dim)        &
                  * ddxyz - 0.5_wp )
             loop_stride = NINT( mask_loop(mid,dim,3) * mask_scale(dim)        &
                  * ddxyz )
             IF ( loop_begin == -1 )  loop_begin = 0  ! avoid negative values
          ELSEIF ( dim == 3 )  THEN
             IF ( mask_loop(mid,dim,2) < 0.0_wp )  THEN
                tmp2 = mask_loop(mid,dim,2)
                mask_loop(mid,dim,2) = zu(nz+1) / mask_scale(dim)   ! (default)
             ENDIF
             IF ( MAXVAL( mask_loop(mid,dim,1:2) )  &
                  > zu(nz+1) / mask_scale(dim) )  THEN
                WRITE ( message_string, '(2(A,I3,A,I1,A,F9.3),A,I1,A,F9.3)' )  &
                     'mask_loop(',mid,',',dim,',1)=',mask_loop(mid,dim,1),     &
                     ' and/or mask_loop(',mid,',',dim,',2)=', &
                     mask_loop(mid,dim,2),' exceed zu(nz+1)/mask_scale(',dim,  &
                     ')=',zu(nz+1)/mask_scale(dim)
                CALL message( 'init_masks', 'PA0333', 1, 2, 0, 6, 0 )
             ENDIF
             ind_array =  &
                  MINLOC( ABS( mask_loop(mid,dim,1) * mask_scale(dim) - zu ) )
             loop_begin =  &
                  ind_array(1) - 1 + nzb ! MINLOC uses lower array bound 1
             ind_array =  &
                  MINLOC( ABS( mask_loop(mid,dim,2) * mask_scale(dim) - zu ) )
             loop_end = ind_array(1) - 1 + nzb ! MINLOC uses lower array bound 1
!
!--          The following line assumes a constant vertical grid spacing within
!--          the vertical mask range; it fails for vertical grid stretching.
!--          Maybe revise later. Issue warning but continue execution. ABS(...)
!--          within the IF statement is necessary because the default value of
!--          dz_stretch_level_start is -9999999.9_wp.
             loop_stride = NINT( mask_loop(mid,dim,3) * mask_scale(dim) * ddxyz )

             IF ( mask_loop(mid,dim,2) * mask_scale(dim) >                     &
                  ABS( dz_stretch_level_start(1) ) )  THEN
                WRITE ( message_string, '(A,I3,A,I1,A,F9.3,A,F8.2,3A)' )       &
                     'mask_loop(',mid,',',dim,',2)=', mask_loop(mid,dim,2),    &
                     ' exceeds dz_stretch_level=',dz_stretch_level_start(1),   &
                     '.&Vertical mask locations will not ',                    &
                     'match the desired heights within the stretching ',       &
                     'region.'
                CALL message( 'init_masks', 'PA0334', 0, 1, 0, 6, 0 )
             ENDIF

          ENDIF
!
!--       If necessary, reset mask_loop(mid,dim,1) and mask_loop(mid,dim,2).
          IF ( tmp1 < 0.0_wp )  mask_loop(mid,dim,1) = tmp1
          IF ( tmp2 < 0.0_wp )  mask_loop(mid,dim,2) = tmp2
!
!--       The default stride +/-1 (every grid point) applies if 
!--       mask_loop(mid,dim,3) is not specified (its default is zero).
          IF ( loop_stride == 0 )  THEN
             IF ( loop_end >= loop_begin )  THEN
                loop_stride =  1
             ELSE
                loop_stride = -1
             ENDIF
          ENDIF
          DO  m = loop_begin, loop_end, loop_stride
             count = count + 1
             IF ( ( m >= lb  .AND.  m <= ub ) .OR.   &
                  ( m == (nxyz+1) .AND. ub == nxyz )  )  THEN
                IF ( count_l == 0 )  mask_start_l(mid,dim) = count
                count_l = count_l + 1
                IF ( dim == 1 )  THEN
                   mask_i(mid,count_l) = m
                ELSEIF ( dim == 2 )  THEN
                   mask_j(mid,count_l) = m
                ELSEIF ( dim == 3 )  THEN
                   mask_k(mid,count_l) = m
                ENDIF
             ENDIF
          ENDDO
          mask_size(mid,dim)   = count
          mask_size_l(mid,dim) = count_l

       ENDIF

    END SUBROUTINE set_mask_locations

 END SUBROUTINE init_masks
