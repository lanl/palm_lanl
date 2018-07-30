!> @file data_output_mask.f90
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
! $Id: data_output_mask.f90 3083 2018-06-19 14:03:12Z gronemeier $
! 
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
! 
! 3030 2018-05-23 14:37:00Z raasch
! variable if renamed ivar
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! 
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison' 
! includes two more prognostic equations for cloud drop concentration (nc)  
! and cloud water content (qc). 
! 
! 2101 2017-01-05 16:42:31Z suehring
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean and rho_av to rho_ocean_av
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1980 2016-07-29 15:51:57Z suehring
! Bugfix, in order to steer user-defined output, setting flag found explicitly
! to .F.
! 
! 1976 2016-07-27 13:28:04Z maronga
! Output of radiation quantities is now done directly in the respective module
! 
! 1960 2016-07-12 16:34:24Z suehring
! Separate humidity and passive scalar
!
! 2016-03-06 18:36:17Z raasch
! name change of netcdf routines and module + related changes,
! switch back of netcdf data format moved from time integration routine to here
!
! 1691 2015-10-26 16:17:44Z maronga
! Added output of radiative heating rates for RRTMG
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1585 2015-04-30 07:05:52Z maronga
! Added support for RRTMG
! 
! 1438 2014-07-22 14:14:06Z heinze
! +nr, qc, qr
! 
! 1359 2014-04-11 17:15:14Z hoffmann
! New particle structure integrated. 
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
! 
! 1327 2014-03-21 11:00:16Z raasch
! 
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
! barrier argument removed from cpu_log,
! module interfaces removed
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1031 2012-10-19 14:35:30Z raasch
! netCDF4 without parallel file support implemented
!
! 1007 2012-09-19 14:30:36Z franke
! Bugfix: calculation of pr must depend on the particle weighting factor,
! missing calculation of ql_vp added
!
! 410 2009-12-04 17:05:40Z letzel
! Initial version
!
! Description:
! ------------
!> Masked data output in netCDF format for current mask (current value of mid).
!------------------------------------------------------------------------------!
 SUBROUTINE data_output_mask( av )

 

#if defined( __netcdf )
    USE arrays_3d,                                                             &
        ONLY:  e, nc, nr, p, pt, q, qc, ql, ql_c, ql_v, qr, rho_ocean, s, sa,  &
               tend, u, v, vpt, w, alpha_T, beta_S, solar3d
    
    USE averaging,                                                             &
        ONLY:  e_av, lpt_av, nc_av, nr_av, p_av, pc_av, pr_av, pt_av, q_av,    &
               qc_av, ql_av, ql_c_av, ql_v_av, ql_vp_av, qv_av, qr_av,         &
               rho_ocean_av, s_av, sa_av, u_av, v_av, vpt_av, w_av,            &
               alpha_T_av, beta_S_av, solar3d_av
    
    USE cloud_parameters,                                                      &
        ONLY:  l_d_cp, pt_d_t
    
    USE control_parameters,                                                    &
        ONLY:  cloud_physics, domask, domask_no, domask_time_count, mask_i,    &
               mask_j, mask_k, mask_size, mask_size_l, mask_start_l,           &
               max_masks, message_string, mid, nz_do3d, simulated_time
    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point

    USE indices,                                                               &
        ONLY:  nbgp, nxl, nxr, nyn, nys, nzb, nzt
        
    USE kinds
    
    USE NETCDF
    
    USE netcdf_interface,                                                      &
        ONLY:  id_set_mask, id_var_domask, id_var_time_mask, nc_stat,          &
               netcdf_data_format, netcdf_handle_error
    
    USE particle_attributes,                                                   &
        ONLY:  grid_particles, number_of_particles, particles,                 &
               particle_advection_start, prt_count
    
    USE pegrid

    USE radiation_model_mod,                                                   &
        ONLY:  radiation, radiation_data_output_mask

    IMPLICIT NONE

    INTEGER(iwp) ::  av       !< 
    INTEGER(iwp) ::  ngp      !< 
    INTEGER(iwp) ::  i        !< 
    INTEGER(iwp) ::  ivar     !<
    INTEGER(iwp) ::  j        !< 
    INTEGER(iwp) ::  k        !< 
    INTEGER(iwp) ::  n        !< 
    INTEGER(iwp) ::  netcdf_data_format_save !<
    INTEGER(iwp) ::  psi      !< 
    INTEGER(iwp) ::  sender   !< 
    INTEGER(iwp) ::  ind(6)   !< 
    
    LOGICAL ::  found         !< 
    LOGICAL ::  resorted      !< 
    
    REAL(wp) ::  mean_r       !< 
    REAL(wp) ::  s_r2         !< 
    REAL(wp) ::  s_r3         !< 
    
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  local_pf    !<
#if defined( __parallel )
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  total_pf    !<
#endif
    REAL(wp), DIMENSION(:,:,:), POINTER ::  to_be_resorted  !<

!
!-- Return, if nothing to output
    IF ( domask_no(mid,av) == 0 )  RETURN

    CALL cpu_log (log_point(49),'data_output_mask','start')

!
!-- Parallel netcdf output is not tested so far for masked data, hence
!-- netcdf_data_format is switched back to non-paralell output.
    netcdf_data_format_save = netcdf_data_format
    IF ( netcdf_data_format == 5 ) netcdf_data_format = 3
    IF ( netcdf_data_format == 6 ) netcdf_data_format = 4

!
!-- Open output file.
    IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
       CALL check_open( 200+mid+av*max_masks )
    ENDIF 

!
!-- Allocate total and local output arrays.
#if defined( __parallel )
    IF ( myid == 0 )  THEN
       ALLOCATE( total_pf(mask_size(mid,1),mask_size(mid,2),mask_size(mid,3)) )
    ENDIF
#endif
    ALLOCATE( local_pf(mask_size_l(mid,1),mask_size_l(mid,2), &
                       mask_size_l(mid,3)) )

!
!-- Update the netCDF time axis.
    domask_time_count(mid,av) = domask_time_count(mid,av) + 1
    IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
       nc_stat = NF90_PUT_VAR( id_set_mask(mid,av), id_var_time_mask(mid,av), &
                               (/ simulated_time /),                          &
                               start = (/ domask_time_count(mid,av) /),       &
                               count = (/ 1 /) )
       CALL netcdf_handle_error( 'data_output_mask', 460 )
    ENDIF

!
!-- Loop over all variables to be written.
    ivar = 1

    DO  WHILE ( domask(mid,av,ivar)(1:1) /= ' ' )
!
!--    Reallocate local_pf on PE 0 since its shape changes during MPI exchange
       IF ( netcdf_data_format < 5   .AND.  myid == 0  .AND.  ivar > 1 )  THEN
          DEALLOCATE( local_pf )
          ALLOCATE( local_pf(mask_size_l(mid,1),mask_size_l(mid,2), &
                             mask_size_l(mid,3)) )
       ENDIF
!
!--    Set flag to steer output of radiation, land-surface, or user-defined
!--    quantities
       found = .FALSE.
!
!--    Store the variable chosen.
       resorted = .FALSE.
       SELECT CASE ( TRIM( domask(mid,av,ivar) ) )

          CASE ( 'e' )
             IF ( av == 0 )  THEN
                to_be_resorted => e
             ELSE
                to_be_resorted => e_av
             ENDIF

          CASE ( 'lpt' )
             IF ( av == 0 )  THEN
                to_be_resorted => pt
             ELSE
                to_be_resorted => lpt_av
             ENDIF

          CASE ( 'nc' )
             IF ( av == 0 )  THEN
                to_be_resorted => nc
             ELSE
                to_be_resorted => nc_av
             ENDIF

          CASE ( 'nr' )
             IF ( av == 0 )  THEN
                to_be_resorted => nr
             ELSE
                to_be_resorted => nr_av
             ENDIF

          CASE ( 'p' )
             IF ( av == 0 )  THEN
                to_be_resorted => p
             ELSE
                to_be_resorted => p_av
             ENDIF

          CASE ( 'pc' )  ! particle concentration (requires ghostpoint exchange)
             IF ( av == 0 )  THEN
                tend = prt_count
                CALL exchange_horiz( tend, nbgp )
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
                      DO  k = 1, mask_size_l(mid,3)
                         local_pf(i,j,k) =  tend(mask_k(mid,k), &
                                   mask_j(mid,j),mask_i(mid,i))
                      ENDDO
                   ENDDO
                ENDDO
                resorted = .TRUE.
             ELSE
                CALL exchange_horiz( pc_av, nbgp )
                to_be_resorted => pc_av
             ENDIF

          CASE ( 'pr' )  ! mean particle radius (effective radius)
             IF ( av == 0 )  THEN
                IF ( simulated_time >= particle_advection_start )  THEN
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         DO  k = nzb, nz_do3d
                            number_of_particles = prt_count(k,j,i)
                            IF (number_of_particles <= 0)  CYCLE
                            particles => grid_particles(k,j,i)%particles(1:number_of_particles)
                            s_r2 = 0.0_wp
                            s_r3 = 0.0_wp
                            DO  n = 1, number_of_particles
                               IF ( particles(n)%particle_mask )  THEN
                                  s_r2 = s_r2 + grid_particles(k,j,i)%particles(n)%radius**2 * &
                                         grid_particles(k,j,i)%particles(n)%weight_factor
                                  s_r3 = s_r3 + grid_particles(k,j,i)%particles(n)%radius**3 * &
                                         grid_particles(k,j,i)%particles(n)%weight_factor
                               ENDIF
                            ENDDO
                            IF ( s_r2 > 0.0_wp )  THEN
                               mean_r = s_r3 / s_r2
                            ELSE
                               mean_r = 0.0_wp
                            ENDIF
                            tend(k,j,i) = mean_r
                         ENDDO
                      ENDDO
                   ENDDO
                   CALL exchange_horiz( tend, nbgp )
                ELSE
                   tend = 0.0_wp
                ENDIF
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
                      DO  k = 1, mask_size_l(mid,3)
                         local_pf(i,j,k) =  tend(mask_k(mid,k), &
                                   mask_j(mid,j),mask_i(mid,i))
                      ENDDO
                   ENDDO
                ENDDO
                resorted = .TRUE.
             ELSE
                CALL exchange_horiz( pr_av, nbgp )
                to_be_resorted => pr_av
             ENDIF

          CASE ( 'pt' )
             IF ( av == 0 )  THEN
                IF ( .NOT. cloud_physics ) THEN
                   to_be_resorted => pt
                ELSE
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
                         DO  k = 1, mask_size_l(mid,3)
                            local_pf(i,j,k) =  &
                                 pt(mask_k(mid,k),mask_j(mid,j),mask_i(mid,i)) &
                                 + l_d_cp * pt_d_t(mask_k(mid,k)) * &
                                   ql(mask_k(mid,k),mask_j(mid,j),mask_i(mid,i))
                         ENDDO
                      ENDDO
                   ENDDO
                   resorted = .TRUE.
                ENDIF
             ELSE
                to_be_resorted => pt_av
             ENDIF

          CASE ( 'q' )
             IF ( av == 0 )  THEN
                to_be_resorted => q
             ELSE
                to_be_resorted => q_av
             ENDIF

          CASE ( 'qc' )
             IF ( av == 0 )  THEN
                to_be_resorted => qc
             ELSE
                to_be_resorted => qc_av
             ENDIF

          CASE ( 'ql' )
             IF ( av == 0 )  THEN
                to_be_resorted => ql
             ELSE
                to_be_resorted => ql_av
             ENDIF

          CASE ( 'ql_c' )
             IF ( av == 0 )  THEN
                to_be_resorted => ql_c
             ELSE
                to_be_resorted => ql_c_av
             ENDIF

          CASE ( 'ql_v' )
             IF ( av == 0 )  THEN
                to_be_resorted => ql_v
             ELSE
                to_be_resorted => ql_v_av
             ENDIF

          CASE ( 'ql_vp' )
             IF ( av == 0 )  THEN
                IF ( simulated_time >= particle_advection_start )  THEN
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         DO  k = nzb, nz_do3d
                            number_of_particles = prt_count(k,j,i)
                            IF (number_of_particles <= 0)  CYCLE
                            particles => grid_particles(k,j,i)%particles(1:number_of_particles)
                            DO  n = 1, number_of_particles
                               IF ( particles(n)%particle_mask )  THEN
                                  tend(k,j,i) = tend(k,j,i) + &
                                          particles(n)%weight_factor / &
                                          prt_count(k,j,i)
                               ENDIF
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                   CALL exchange_horiz( tend, nbgp )
                ELSE
                   tend = 0.0_wp
                ENDIF
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
                      DO  k = 1, mask_size_l(mid,3)
                         local_pf(i,j,k) =  tend(mask_k(mid,k), &
                                   mask_j(mid,j),mask_i(mid,i))
                      ENDDO
                   ENDDO
                ENDDO
                resorted = .TRUE.
             ELSE
                CALL exchange_horiz( ql_vp_av, nbgp )
                to_be_resorted => ql_vp_av
             ENDIF

          CASE ( 'qv' )
             IF ( av == 0 )  THEN
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
                      DO  k = 1, mask_size_l(mid,3)
                         local_pf(i,j,k) =  &
                              q(mask_k(mid,k),mask_j(mid,j),mask_i(mid,i)) -  &
                              ql(mask_k(mid,k),mask_j(mid,j),mask_i(mid,i))
                      ENDDO
                   ENDDO
                ENDDO
                resorted = .TRUE.
             ELSE
                to_be_resorted => qv_av
             ENDIF

          CASE ( 'qr' )
             IF ( av == 0 )  THEN
                to_be_resorted => qr
             ELSE
                to_be_resorted => qr_av
             ENDIF

          CASE ( 'rho_ocean' )
             IF ( av == 0 )  THEN
                to_be_resorted => rho_ocean
             ELSE
                to_be_resorted => rho_ocean_av
             ENDIF

          CASE ( 'solar3d' )
             IF ( av == 0 )  THEN
                to_be_resorted => solar3d
             ELSE
                to_be_resorted => solar3d_av
             ENDIF

          CASE ( 'alpha_T' )
             IF ( av == 0 )  THEN
                to_be_resorted => alpha_T
             ELSE
                to_be_resorted => alpha_T_av
             ENDIF

          CASE ( 'beta_S' )
             IF ( av == 0 )  THEN
                to_be_resorted => beta_S
             ELSE
                to_be_resorted => beta_S_av
             ENDIF


          CASE ( 's' )
             IF ( av == 0 )  THEN
                to_be_resorted => s
             ELSE
                to_be_resorted => s_av
             ENDIF

          CASE ( 'sa' )
             IF ( av == 0 )  THEN
                to_be_resorted => sa
             ELSE
                to_be_resorted => sa_av
             ENDIF

          CASE ( 'u' )
             IF ( av == 0 )  THEN
                to_be_resorted => u
             ELSE
                to_be_resorted => u_av
             ENDIF

          CASE ( 'v' )
             IF ( av == 0 )  THEN
                to_be_resorted => v
             ELSE
                to_be_resorted => v_av
             ENDIF

          CASE ( 'vpt' )
             IF ( av == 0 )  THEN
                to_be_resorted => vpt
             ELSE
                to_be_resorted => vpt_av
             ENDIF

          CASE ( 'w' )
             IF ( av == 0 )  THEN
                to_be_resorted => w
             ELSE
                to_be_resorted => w_av
             ENDIF

          CASE DEFAULT

!
!--          Radiation quantity
             IF ( radiation )  THEN
                CALL radiation_data_output_mask(av, domask(mid,av,ivar), found,&
                                                local_pf )
             ENDIF

!
!--          User defined quantity
             IF ( .NOT. found )  THEN
                CALL user_data_output_mask(av, domask(mid,av,ivar), found,     &
                                           local_pf )
             ENDIF

             resorted = .TRUE.

             IF ( .NOT. found )  THEN
                WRITE ( message_string, * ) 'no masked output available for: ',&
                                            TRIM( domask(mid,av,ivar) )
                CALL message( 'data_output_mask', 'PA0327', 0, 0, 0, 6, 0 )
             ENDIF

       END SELECT

!
!--    Resort the array to be output, if not done above
       IF ( .NOT. resorted )  THEN
          DO  i = 1, mask_size_l(mid,1)
             DO  j = 1, mask_size_l(mid,2)
                DO  k = 1, mask_size_l(mid,3)
                   local_pf(i,j,k) =  to_be_resorted(mask_k(mid,k), &
                                      mask_j(mid,j),mask_i(mid,i))
                ENDDO
             ENDDO
          ENDDO
       ENDIF

!
!--    I/O block. I/O methods are implemented
!--    (1) for parallel execution
!--     a. with netCDF 4 parallel I/O-enabled library
!--     b. with netCDF 3 library
!--    (2) for serial execution.
!--    The choice of method depends on the correct setting of preprocessor
!--    directives __parallel and __netcdf4_parallel as well as on the parameter
!--    netcdf_data_format.
#if defined( __parallel )
#if defined( __netcdf4_parallel )
       IF ( netcdf_data_format > 4 )  THEN
!
!--       (1) a. Parallel I/O using netCDF 4 (not yet tested)
          nc_stat = NF90_PUT_VAR( id_set_mask(mid,av),                         &
               id_var_domask(mid,av,ivar), local_pf,                           &
               start = (/ mask_start_l(mid,1), mask_start_l(mid,2),            &
                          mask_start_l(mid,3), domask_time_count(mid,av) /),   &
               count = (/ mask_size_l(mid,1), mask_size_l(mid,2),              &
                          mask_size_l(mid,3), 1 /) )
          CALL netcdf_handle_error( 'data_output_mask', 461 )
       ELSE
#endif
!
!--       (1) b. Conventional I/O only through PE0
!--       PE0 receives partial arrays from all processors of the respective mask
!--       and outputs them. Here a barrier has to be set, because otherwise 
!--       "-MPI- FATAL: Remote protocol queue full" may occur.
          CALL MPI_BARRIER( comm2d, ierr )

          ngp = mask_size_l(mid,1) * mask_size_l(mid,2) * mask_size_l(mid,3)
          IF ( myid == 0 )  THEN
!
!--          Local array can be relocated directly.
             total_pf( &
               mask_start_l(mid,1):mask_start_l(mid,1)+mask_size_l(mid,1)-1, &
               mask_start_l(mid,2):mask_start_l(mid,2)+mask_size_l(mid,2)-1, &
               mask_start_l(mid,3):mask_start_l(mid,3)+mask_size_l(mid,3)-1 ) &
               = local_pf
!
!--          Receive data from all other PEs.
             DO  n = 1, numprocs-1
!
!--             Receive index limits first, then array.
!--             Index limits are received in arbitrary order from the PEs.
                CALL MPI_RECV( ind(1), 6, MPI_INTEGER, MPI_ANY_SOURCE, 0,  &
                     comm2d, status, ierr )
!
!--             Not all PEs have data for the mask
                IF ( ind(1) /= -9999 )  THEN
                   ngp = ( ind(2)-ind(1)+1 ) * (ind(4)-ind(3)+1 ) *  &
                         ( ind(6)-ind(5)+1 )
                   sender = status(MPI_SOURCE)
                   DEALLOCATE( local_pf )
                   ALLOCATE(local_pf(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)))
                   CALL MPI_RECV( local_pf(ind(1),ind(3),ind(5)), ngp,  &
                        MPI_REAL, sender, 1, comm2d, status, ierr )
                   total_pf(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) &
                        = local_pf
                ENDIF
             ENDDO

             nc_stat = NF90_PUT_VAR( id_set_mask(mid,av),                      &
                  id_var_domask(mid,av,ivar), total_pf,                        &
                  start = (/ 1, 1, 1, domask_time_count(mid,av) /),            &
                  count = (/ mask_size(mid,1), mask_size(mid,2),               &
                             mask_size(mid,3), 1 /) )
             CALL netcdf_handle_error( 'data_output_mask', 462 )

          ELSE
!
!--          If at least part of the mask resides on the PE, send the index
!--          limits for the target array, otherwise send -9999 to PE0.
             IF ( mask_size_l(mid,1) > 0 .AND.  mask_size_l(mid,2) > 0 .AND. &
                  mask_size_l(mid,3) > 0  ) &
                  THEN
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
!--          If applicable, send data to PE0.
             IF ( ind(1) /= -9999 )  THEN
                CALL MPI_SEND( local_pf(1,1,1), ngp, MPI_REAL, 0, 1, comm2d, &
                     ierr )
             ENDIF
          ENDIF
!
!--       A barrier has to be set, because otherwise some PEs may proceed too
!--       fast so that PE0 may receive wrong data on tag 0.
          CALL MPI_BARRIER( comm2d, ierr )
#if defined( __netcdf4_parallel )
       ENDIF
#endif
#else
!
!--    (2) For serial execution of PALM, the single processor (PE0) holds all
!--    data and writes them directly to file.
       nc_stat = NF90_PUT_VAR( id_set_mask(mid,av),                            &
                               id_var_domask(mid,av,ivar), local_pf,           &
                             start = (/ 1, 1, 1, domask_time_count(mid,av) /), &
                             count = (/ mask_size_l(mid,1), mask_size_l(mid,2),&
                               mask_size_l(mid,3), 1 /) )
       CALL netcdf_handle_error( 'data_output_mask', 463 )
#endif

       ivar = ivar + 1

    ENDDO

!
!-- Deallocate temporary arrays.
    DEALLOCATE( local_pf )
#if defined( __parallel )
    IF ( myid == 0 )  THEN
       DEALLOCATE( total_pf )
    ENDIF
#endif

!
!-- Switch back to original format given by user (see beginning of this routine)
    netcdf_data_format = netcdf_data_format_save

    CALL cpu_log( log_point(49), 'data_output_mask', 'stop' )
#endif

 END SUBROUTINE data_output_mask
