!> @file large_scale_forcing_nudging_mod.f90
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
! $Id: large_scale_forcing_nudging_mod.f90 3049 2018-05-29 13:52:36Z Giersch $
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised 
! 
! 3026 2018-05-22 10:30:53Z schwenkel
! Changed the name specific humidity to mixing ratio, since we are computing
! mixing ratios.
! 
! 2970 2018-04-13 15:09:23Z suehring
! Bugfix in old large-scale forcing mode
! 
! 2938 2018-03-27 15:52:42Z suehring
! Further improvements for nesting in larger-scale model
! 
! 2863 2018-03-08 11:36:25Z suehring
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Forcing with larger-scale models implemented (MS)
! 
! 2342 2017-08-08 11:00:43Z boeske
! fixed check if surface forcing data is available until end of simulation
! 
! 2320 2017-07-21 12:47:43Z suehring
! initial revision
!
! Description:
! ------------
!> Calculates large scale forcings (geostrophic wind and subsidence velocity) as 
!> well as surfaces fluxes dependent on time given in an external file (LSF_DATA).
!> Moreover, module contains nudging routines, where u, v, pt and q are nudged 
!> to given profiles on a relaxation timescale tnudge. 
!> Profiles are read in from NUDGING_DATA. 
!> Code is based on Neggers et al. (2012) and also in parts on DALES and UCLA-LES.
!> @todo: Revise reading of ASCII-files
!> @todo: Remove unused variables and control flags
!> @todo: Revise large-scale facing of surface variables
!> @todo: Revise control flags lsf_exception, lsf_surf, lsf_vert, etc. 
!--------------------------------------------------------------------------------!
 MODULE lsf_nudging_mod

    USE arrays_3d,                                                             &
        ONLY:  dzw, e, heatflux_input_conversion, pt, pt_init, q, q_init, s,   &
               tend, u, u_init, ug, v, v_init, vg, w, w_subs,                  &
               waterflux_input_conversion, zu, zw                  

    USE control_parameters,                                                    &
        ONLY:  bc_lr, bc_ns, bc_pt_b, bc_q_b, constant_diffusion,              &
               constant_heatflux, constant_waterflux,                          &
               data_output_pr, dt_3d, end_time, forcing,                       &
               force_bound_l, force_bound_n, force_bound_r, force_bound_s,     &  
               humidity, initializing_actions, intermediate_timestep_count,    &
               ibc_pt_b, ibc_q_b,                                              &
               large_scale_forcing, large_scale_subsidence, lsf_surf, lsf_vert,&
               lsf_exception, message_string, neutral, nudging, passive_scalar,&
               pt_surface, ocean, q_surface, surface_heatflux,                 &
               surface_pressure, surface_waterflux, topography,                &
               use_subsidence_tendencies

    USE grid_variables

    USE pegrid

    USE indices,                                                               &
        ONLY:  nbgp, ngp_sums_ls, nx, nxl, nxlg, nxlu, nxr, nxrg, ny, nys,     &
               nysv, nysg, nyn, nyng, nzb, nz, nzt, wall_flags_0

    USE kinds

    USE surface_mod,                                                           &
        ONLY:  surf_def_h, surf_lsm_h, surf_usm_h

    USE statistics,                                                            &
        ONLY:  hom, statistic_regions, sums_ls_l, weight_substep

    USE netcdf_data_input_mod,                                                 &
        ONLY:  force, netcdf_data_input_interpolate

    INTEGER(iwp) ::  nlsf = 1000                       !< maximum number of profiles in LSF_DATA (large scale forcing)
    INTEGER(iwp) ::  ntnudge = 1000                    !< maximum number of profiles in NUDGING_DATA (nudging)

    REAL(wp) ::  d_area_t

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ptnudge     !< vertical profile of pot. temperature interpolated to vertical grid (nudging)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  qnudge      !< vertical profile of water vapor mixing ratio interpolated to vertical grid (nudging) 
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tnudge      !< vertical profile of nudging time scale interpolated to vertical grid (nudging)  
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  td_lsa_lpt  !< temperature tendency due to large scale advection (large scale forcing)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  td_lsa_q    !< water vapor mixing ratio tendency due to large scale advection (large scale forcing)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  td_sub_lpt  !< temperature tendency due to subsidence/ascent (large scale forcing)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  td_sub_q    !< water vapor mixing ratio tendency due to subsidence/ascent (large scale forcing)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ug_vert     !< vertical profile of geostrophic wind component in x-direction interpolated to vertical grid (large scale forcing)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  unudge      !< vertical profile of wind component in x-direction interpolated to vertical grid (nudging) 
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vnudge      !< vertical profile of wind component in y-direction interpolated to vertical grid (nudging) 
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vg_vert     !< vertical profile of geostrophic wind component in y-direction interpolated to vertical grid (large scale forcing)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  wnudge      !< vertical profile of subsidence/ascent velocity interpolated to vertical grid (nudging) ???
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  wsubs_vert  !< vertical profile of wind component in z-direction interpolated to vertical grid (nudging) ???

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_surf      !< time-dependent surface sensible heat flux (large scale forcing)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  timenudge     !< times at which vertical profiles are defined in NUDGING_DATA (nudging)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  time_surf     !< times at which surface values/fluxes are defined in LSF_DATA (large scale forcing)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  time_vert     !< times at which vertical profiles are defined in LSF_DATA (large scale forcing)

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  tmp_tnudge    !< current nudging time scale

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  p_surf        !< time-dependent surface pressure (large scale forcing)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_surf       !< time-dependent surface temperature (large scale forcing)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_surf     !< time-dependent surface latent heat flux (large scale forcing)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  q_surf        !< time-dependent surface water vapor mixing ratio (large scale forcing)

    SAVE
    PRIVATE
!
!-- Public subroutines
    PUBLIC ls_forcing_surf, ls_forcing_vert, ls_advec, lsf_init,               &
           lsf_nudging_check_parameters, nudge_init,                           &
           lsf_nudging_check_data_output_pr, lsf_nudging_header,               &
           calc_tnudge, nudge, nudge_ref, forcing_bc_mass_conservation,        &
           forcing_bc
!
!-- Public variables
    PUBLIC qsws_surf, shf_surf, td_lsa_lpt, td_lsa_q, td_sub_lpt,              &
           td_sub_q, time_vert, force


    INTERFACE ls_advec
       MODULE PROCEDURE ls_advec
       MODULE PROCEDURE ls_advec_ij
    END INTERFACE ls_advec

    INTERFACE nudge
       MODULE PROCEDURE nudge
       MODULE PROCEDURE nudge_ij
    END INTERFACE nudge

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE forcing_bc_mass_conservation

       USE control_parameters,                                                 &
           ONLY:  volume_flow

       IMPLICIT NONE

       INTEGER(iwp) ::  i !<
       INTEGER(iwp) ::  j !<
       INTEGER(iwp) ::  k !<

       REAL(wp) ::  w_correct !< 
       REAL(wp), DIMENSION(1:3) ::  volume_flow_l   !<

       volume_flow   = 0.0_wp
       volume_flow_l = 0.0_wp

       d_area_t = 1.0_wp / ( ( nx + 1 ) * dx * ( ny + 1 ) * dy )

       IF ( force_bound_l )  THEN
          i = nxl
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                volume_flow_l(1) = volume_flow_l(1) + u(k,j,i) * dzw(k) * dy   &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,i), 1 ) )
             ENDDO
          ENDDO
       ENDIF
       IF ( force_bound_r )  THEN
          i = nxr+1
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                volume_flow_l(1) = volume_flow_l(1) - u(k,j,i) * dzw(k) * dy   &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,i), 1 ) )
             ENDDO
          ENDDO
       ENDIF
       IF ( force_bound_s )  THEN
          j = nys
          DO  i = nxl, nxr
             DO  k = nzb+1, nzt
                volume_flow_l(2) = volume_flow_l(2) + v(k,j,i) * dzw(k) * dx   &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,i), 2 ) )
             ENDDO
          ENDDO
       ENDIF
       IF ( force_bound_n )  THEN
          j = nyn+1
          DO  i = nxl, nxr
             DO  k = nzb+1, nzt
                volume_flow_l(2) = volume_flow_l(2) - v(k,j,i) * dzw(k) * dx   &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,i), 2 ) )
             ENDDO
          ENDDO
       ENDIF
!
!--    Top boundary
       k = nzt
       DO  i = nxl, nxr
          DO  j = nys, nyn
             volume_flow_l(3) = volume_flow_l(3) - w(k,j,i) * dx * dy
          ENDDO
       ENDDO

#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( volume_flow_l, volume_flow, 3, MPI_REAL, MPI_SUM,      &
                           comm2d, ierr )
#else
       volume_flow = volume_flow_l
#endif

       w_correct = SUM( volume_flow ) * d_area_t

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzt, nzt + 1
                w(k,j,i) = w(k,j,i) + w_correct
             ENDDO
          ENDDO
       ENDDO

write(9,*) "w correction", w_correct
flush(9)

    END SUBROUTINE forcing_bc_mass_conservation


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE forcing_bc

       USE control_parameters,                                                 &
           ONLY:  force_bound_l, force_bound_n, force_bound_r, force_bound_s,  &
                  humidity, neutral, passive_scalar, simulated_time

       USE netcdf_data_input_mod,                                              &
           ONLY:  force                              

       IMPLICIT NONE

       INTEGER(iwp) ::  i !< running index x-direction
       INTEGER(iwp) ::  j !< running index y-direction
       INTEGER(iwp) ::  k !< running index z-direction
       INTEGER(iwp) ::  t !< running index for time levels

       REAL(wp) ::  ddt_lsf !< inverse value of time resolution of forcing data
       REAL(wp) ::  t_ref   !< time past since last reference step

!
!--    If required, interpolate and/or extrapolate data vertically. This is 
!--    required as Inifor outputs only equidistant vertical data.
       IF ( ANY( zu(1:nzt+1) /= force%zu_atmos(1:force%nzu) ) )  THEN
          IF ( .NOT. force%interpolated )  THEN

             DO  t = 0, 1
                IF ( force_bound_l )  THEN
                   CALL netcdf_data_input_interpolate( force%u_left(t,:,:),    &
                                                       zu(nzb+1:nzt+1),        &
                                                       force%zu_atmos )
                   CALL netcdf_data_input_interpolate( force%v_left(t,:,:),    &
                                                       zu(nzb+1:nzt+1),        &
                                                       force%zu_atmos )
                   CALL netcdf_data_input_interpolate( force%w_left(t,:,:),    &
                                                       zw(nzb+1:nzt+1),        &
                                                       force%zw_atmos )
                   IF ( .NOT. neutral )                                        &
                      CALL netcdf_data_input_interpolate( force%pt_left(t,:,:),&
                                                          zu(nzb+1:nzt+1),     &
                                                          force%zu_atmos )
                   IF ( humidity )                                             &
                      CALL netcdf_data_input_interpolate( force%q_left(t,:,:), &
                                                          zu(nzb+1:nzt+1),     &
                                                          force%zu_atmos )
                ENDIF
                IF ( force_bound_r )  THEN
                   CALL netcdf_data_input_interpolate( force%u_right(t,:,:),   &
                                                       zu(nzb+1:nzt+1),        &
                                                       force%zu_atmos )
                   CALL netcdf_data_input_interpolate( force%v_right(t,:,:),   &
                                                       zu(nzb+1:nzt+1),        &
                                                       force%zu_atmos )
                   CALL netcdf_data_input_interpolate( force%w_right(t,:,:),   &
                                                       zw(nzb+1:nzt+1),        &
                                                       force%zw_atmos )
                   IF ( .NOT. neutral )                                        &
                      CALL netcdf_data_input_interpolate( force%pt_right(t,:,:),&
                                                          zu(nzb+1:nzt+1),     &
                                                          force%zu_atmos )
                   IF ( humidity )                                             &
                      CALL netcdf_data_input_interpolate( force%q_right(t,:,:),&
                                                          zu(nzb+1:nzt+1),     &
                                                          force%zu_atmos )
                ENDIF
                IF ( force_bound_n )  THEN
                   CALL netcdf_data_input_interpolate( force%u_north(t,:,:),   &
                                                       zu(nzb+1:nzt+1),        &
                                                       force%zu_atmos )
                   CALL netcdf_data_input_interpolate( force%v_north(t,:,:),   &
                                                       zu(nzb+1:nzt+1),        &
                                                       force%zu_atmos )
                   CALL netcdf_data_input_interpolate( force%w_north(t,:,:),   &
                                                       zw(nzb+1:nzt+1),        &
                                                       force%zw_atmos )
                   IF ( .NOT. neutral )                                        &
                      CALL netcdf_data_input_interpolate( force%pt_north(t,:,:),&
                                                          zu(nzb+1:nzt+1),     &
                                                          force%zu_atmos )
                   IF ( humidity )                                             &
                      CALL netcdf_data_input_interpolate( force%q_north(t,:,:),&
                                                          zu(nzb+1:nzt+1),     &
                                                          force%zu_atmos )
                ENDIF
                IF ( force_bound_s )  THEN
                   CALL netcdf_data_input_interpolate( force%u_south(t,:,:),   &
                                                       zu(nzb+1:nzt+1),        &
                                                       force%zu_atmos )
                   CALL netcdf_data_input_interpolate( force%v_south(t,:,:),   &
                                                       zu(nzb+1:nzt+1),        &
                                                       force%zu_atmos )
                   CALL netcdf_data_input_interpolate( force%w_south(t,:,:),   &
                                                       zw(nzb+1:nzt+1),        &
                                                       force%zw_atmos )
                   IF ( .NOT. neutral )                                        &
                      CALL netcdf_data_input_interpolate( force%pt_south(t,:,:),&
                                                          zu(nzb+1:nzt+1),     &
                                                          force%zu_atmos )
                   IF ( humidity )                                             &
                      CALL netcdf_data_input_interpolate( force%q_south(t,:,:),&
                                                          zu(nzb+1:nzt+1),     &
                                                          force%zu_atmos )
                ENDIF
             ENDDO
!
!--          Note, no interpolation of top boundary. Just use initial value. 
!--          No physical meaningful extrapolation possible if only one layer is
!--          given. 

             force%interpolated = .TRUE. 
          ENDIF
       ENDIF
      
!
!--    Calculate time interval of forcing data       
       ddt_lsf = 1.0_wp / ( force%time(force%tind_p) - force%time(force%tind) )
!
!--    Calculate reziproke time past since last reference step. Please note, 
!--    as simulated time is still not updated, the actual time here is 
!--    simulated time + dt_3d
       t_ref = simulated_time + dt_3d - force%time(force%tind)

       IF ( force_bound_l )  THEN

          DO  j = nys, nyn
             DO  k = nzb+1, nzt+1
                u(k,j,nxlg:nxl) = force%u_left(0,k,j) + ddt_lsf * t_ref *      &
                         ( force%u_left(1,k,j) - force%u_left(0,k,j) ) *       &
                           MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_0(k,j,nxlg:nxl), 1 ) )
             ENDDO
          ENDDO

          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                w(k,j,nxlg:nxl-1) = force%w_left(0,k,j) + ddt_lsf * t_ref *    &
                         ( force%w_left(1,k,j) - force%w_left(0,k,j) ) *       &
                           MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_0(k,j,nxlg:nxl-1), 3 ) )
             ENDDO
          ENDDO

          DO  j = nysv, nyn
             DO  k = nzb+1, nzt+1
                v(k,j,nxlg:nxl-1) = force%v_left(0,k,j) + ddt_lsf * t_ref *    &
                         ( force%v_left(1,k,j) - force%v_left(0,k,j) ) *       &
                           MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_0(k,j,nxlg:nxl-1), 2 ) )
             ENDDO
          ENDDO

          IF ( .NOT. neutral )  THEN
             DO  j = nys, nyn
                DO  k = nzb+1, nzt+1
                   pt(k,j,nxlg:nxl-1) = force%pt_left(0,k,j) + ddt_lsf *       &
                                                                  t_ref   *    &
                       ( force%pt_left(1,k,j) - force%pt_left(0,k,j) )
 
                ENDDO
             ENDDO
          ENDIF

       ENDIF

       IF ( force_bound_r )  THEN

          DO  j = nys, nyn
             DO  k = nzb+1, nzt+1
                u(k,j,nxr+1:nxrg) = force%u_right(0,k,j) + ddt_lsf * t_ref *   &
                        ( force%u_right(1,k,j) - force%u_right(0,k,j) ) *      &
                           MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_0(k,j,nxr+1:nxrg), 1 ) )

             ENDDO
          ENDDO
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                w(k,j,nxr+1:nxrg) = force%w_right(0,k,j) + ddt_lsf * t_ref *   &
                        ( force%w_right(1,k,j) - force%w_right(0,k,j) ) *      &
                           MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_0(k,j,nxr+1:nxrg), 3 ) )
             ENDDO
          ENDDO

          DO  j = nysv, nyn
             DO  k = nzb+1, nzt+1
                v(k,j,nxr+1:nxrg) = force%v_right(0,k,j) + ddt_lsf * t_ref *   &
                        ( force%v_right(1,k,j) - force%v_right(0,k,j) ) *      &
                           MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_0(k,j,nxr+1:nxrg), 2 ) )
             ENDDO
          ENDDO

          IF ( .NOT. neutral )  THEN
             DO  j = nys, nyn
                DO  k = nzb+1, nzt+1
                   pt(k,j,nxr+1:nxrg) = force%pt_right(0,k,j) + ddt_lsf *      &
                                                                  t_ref   *    &
                     ( force%pt_right(1,k,j) - force%pt_right(0,k,j) )
 
                ENDDO
             ENDDO
          ENDIF

       ENDIF

       IF ( force_bound_s )  THEN

          DO  i = nxl, nxr
             DO  k = nzb+1, nzt+1
                v(k,nysg:nys,i)   = force%v_south(0,k,i) + ddt_lsf * t_ref *   &
                        ( force%v_south(1,k,i) - force%v_south(0,k,i) ) *      &
                           MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_0(k,nysg:nys,i), 2 ) )
             ENDDO
          ENDDO

          DO  i = nxl, nxr
             DO  k = nzb+1, nzt
                w(k,nysg:nys-1,i) = force%w_south(0,k,i) + ddt_lsf * t_ref *   &
                        ( force%w_south(1,k,i) - force%w_south(0,k,i) ) *      &
                           MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_0(k,nysg:nys-1,i), 3 ) )
             ENDDO
          ENDDO

          DO  i = nxlu, nxr
             DO  k = nzb+1, nzt+1
                u(k,nysg:nys-1,i) = force%u_south(0,k,i) + ddt_lsf * t_ref *   &
                        ( force%u_south(1,k,i) - force%u_south(0,k,i) ) *      &
                           MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_0(k,nysg:nys-1,i), 1 ) )
             ENDDO
          ENDDO

          IF ( .NOT. neutral )  THEN
             DO  i = nxl, nxr
                DO  k = nzb+1, nzt+1
                   pt(k,nysg:nys-1,i) = force%pt_south(0,k,i) + ddt_lsf *      &
                                                                  t_ref   *    &
                     ( force%pt_south(1,k,i) - force%pt_south(0,k,i) )
 
                ENDDO
             ENDDO
          ENDIF

       ENDIF

       IF ( force_bound_n )  THEN

          DO  i = nxl, nxr
             DO  k = nzb+1, nzt+1
                v(k,nyn+1:nyng,i)   = force%v_north(0,k,i) + ddt_lsf * t_ref * &
                        ( force%v_north(1,k,i) - force%v_north(0,k,i) ) *      &
                           MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_0(k,nyn+1:nyng,i), 2 ) )
             ENDDO
          ENDDO
          DO  i = nxl, nxr
             DO  k = nzb+1, nzt
                w(k,nyn+1:nyng,i) = force%w_north(0,k,i) + ddt_lsf * t_ref *   &
                        ( force%w_north(1,k,i) - force%w_north(0,k,i) ) *      &
                           MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_0(k,nyn+1:nyng,i), 3 ) )
             ENDDO
          ENDDO

          DO  i = nxlu, nxr
             DO  k = nzb+1, nzt+1
                u(k,nyn+1:nyng,i) = force%u_north(0,k,i) + ddt_lsf * t_ref *   &
                        ( force%u_north(1,k,i) - force%u_north(0,k,i) ) *      &
                           MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_0(k,nyn+1:nyng,i), 1 ) )

             ENDDO
          ENDDO

          IF ( .NOT. neutral )  THEN
             DO  i = nxl, nxr
                DO  k = nzb+1, nzt+1
                   pt(k,nyn+1:nyng,i) = force%pt_north(0,k,i) + ddt_lsf *      &
                                                                  t_ref   *    &
                     ( force%pt_north(1,k,i) - force%pt_north(0,k,i) )
 
                ENDDO
             ENDDO
          ENDIF

       ENDIF
!
!--    Top boundary. 
!--    Please note, only map Inifor data on model top in case the numeric is
!--    identical to the Inifor grid. At the top boundary an extrapolation is 
!--    not possible. 
       DO  i = nxlu, nxr
          DO  j = nys, nyn
             u(nzt+1,j,i) = force%u_top(0,j,i) + ddt_lsf * t_ref *             &
                        ( force%u_top(1,j,i) - force%u_top(0,j,i) ) *          &
                           MERGE( 1.0_wp, 0.0_wp,                             &
                                  BTEST( wall_flags_0(nzt+1,j,i), 1 ) )
          ENDDO
       ENDDO

       DO  i = nxl, nxr
          DO  j = nysv, nyn
             v(nzt+1,j,i) = force%v_top(0,j,i) + ddt_lsf * t_ref *             &
                        ( force%v_top(1,j,i) - force%v_top(0,j,i) ) *          &
                           MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_0(nzt+1,j,i), 2 ) )
          ENDDO
       ENDDO

       DO  i = nxl, nxr
          DO  j = nys, nyn
             w(nzt:nzt+1,j,i) = force%w_top(0,j,i) + ddt_lsf * t_ref *         &
                        ( force%w_top(1,j,i) - force%w_top(0,j,i) ) *          &
                           MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_0(nzt:nzt+1,j,i), 3 ) )
          ENDDO
       ENDDO


       IF ( .NOT. neutral )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                pt(nzt+1,j,i) = force%pt_top(0,j,i) + ddt_lsf * t_ref *        &
                        ( force%pt_top(1,j,i) - force%pt_top(0,j,i) )
             ENDDO
          ENDDO
       ENDIF

       IF ( force_bound_l  .AND.  force_bound_s )  THEN
          DO  i = 1, nbgp
             u(:,nys-i,nxlg:nxl)   = u(:,nys,nxlg:nxl)
             w(:,nys-i,nxlg:nxl-1) = w(:,nys,nxlg:nxl-1)
             IF ( .NOT. neutral )  pt(:,nys-i,nxlg:nxl-1) = pt(:,nys,nxlg:nxl-1)
          ENDDO
          DO  i = 1, nbgp+1
             v(:,nysv-i,nxlg:nxl-1) = v(:,nysv,nxlg:nxl-1)
          ENDDO
       ENDIF
       IF ( force_bound_l  .AND.  force_bound_n )  THEN
          DO  i = 1, nbgp
             u(:,nyn+i,nxlg:nxl)   = u(:,nyn,nxlg:nxl)
             v(:,nyn+i,nxlg:nxl-1) = v(:,nyn,nxlg:nxl-1)
             w(:,nyn+i,nxlg:nxl-1) = w(:,nyn,nxlg:nxl-1)
             IF ( .NOT. neutral )  pt(:,nyn+i,nxlg:nxl-1) = pt(:,nyn,nxlg:nxl-1)
          ENDDO
       ENDIF
       IF ( force_bound_r  .AND.  force_bound_s )  THEN
          DO  i = 1, nbgp
             u(:,nys-i,nxr+1:nxrg) = u(:,nys,nxr+1:nxrg)
             w(:,nys-i,nxr+1:nxrg) = w(:,nys,nxr+1:nxrg)
             IF ( .NOT. neutral )  pt(:,nys-i,nxr+1:nxrg) = pt(:,nys,nxr+1:nxrg)
          ENDDO
          DO  i = 1, nbgp+1
             v(:,nysv-i,nxr+1:nxrg) = v(:,nysv,nxr+1:nxrg)
          ENDDO
       ENDIF
       IF ( force_bound_r  .AND.  force_bound_n )  THEN
          DO  i = 1, nbgp
             u(:,nyn+i,nxr+1:nxrg) = u(:,nyn,nxr+1:nxrg)
             v(:,nyn+i,nxr+1:nxrg) = v(:,nyn,nxr+1:nxrg)
             w(:,nyn+i,nxr+1:nxrg) = w(:,nyn,nxr+1:nxrg)
             IF ( .NOT. neutral )  pt(:,nyn+i,nxr+1:nxrg) = pt(:,nyn,nxr+1:nxrg)
          ENDDO
       ENDIF
!
!--    Moreover, set Neumann boundary condition for subgrid-scale TKE and 
!--    passive scalar
       IF ( .NOT. constant_diffusion )  THEN
          IF (  force_bound_l )  e(:,:,nxl-1) = e(:,:,nxl)
          IF (  force_bound_r )  e(:,:,nxr+1) = e(:,:,nxr)
          IF (  force_bound_s )  e(:,nys-1,:) = e(:,nys,:)
          IF (  force_bound_n )  e(:,nyn+1,:) = e(:,nyn,:)
          e(nzt+1,:,:) = e(nzt,:,:)
       ENDIF
       IF ( passive_scalar )  THEN
          IF (  force_bound_l )  s(:,:,nxl-1) = s(:,:,nxl)
          IF (  force_bound_r )  s(:,:,nxr+1) = s(:,:,nxr)
          IF (  force_bound_s )  s(:,nys-1,:) = s(:,nys,:)
          IF (  force_bound_n )  s(:,nyn+1,:) = s(:,nyn,:)
       ENDIF


       CALL exchange_horiz( u, nbgp )
       CALL exchange_horiz( v, nbgp )
       CALL exchange_horiz( w, nbgp )
       IF ( .NOT. neutral )  CALL exchange_horiz( pt, nbgp )

!
!--    Set surface pressure. Please note, time-dependent surface
!--    pressure would require changes in anelastic approximation and 
!--    treatment of fluxes. 
!--    For the moment, comment this out! 
!      surface_pressure = force%surface_pressure(force%tind) +                 &
!                                                      ddt_lsf * t_ref *       &
!                                    ( force%surface_pressure(force%tind_p)    &
!                                    - force%surface_pressure(force%tind) )

    END SUBROUTINE forcing_bc

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE lsf_nudging_check_parameters

       IMPLICIT NONE
!
!--    Check nudging and large scale forcing from external file
       IF ( nudging  .AND.  (  .NOT.  large_scale_forcing ) )  THEN
          message_string = 'Nudging requires large_scale_forcing = .T.. &'//   &
                        'Surface fluxes and geostrophic wind should be &'//    &
                        'prescribed in file LSF_DATA'
          CALL message( 'check_parameters', 'PA0374', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( large_scale_forcing  .AND.  ( bc_lr /= 'cyclic'  .OR.              &
                                          bc_ns /= 'cyclic' ) )  THEN
          message_string = 'Non-cyclic lateral boundaries do not allow for &'//&
                        'the usage of large scale forcing from external file.'
          CALL message( 'check_parameters', 'PA0375', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( large_scale_forcing  .AND.  (  .NOT.  humidity ) )  THEN
          message_string = 'The usage of large scale forcing from external &'//&
                        'file LSF_DATA requires humidity = .T..'
          CALL message( 'check_parameters', 'PA0376', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( large_scale_forcing  .AND.  passive_scalar )  THEN
          message_string = 'The usage of large scale forcing from external &'// &
                        'file LSF_DATA is not implemented for passive scalars'
          CALL message( 'check_parameters', 'PA0440', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( large_scale_forcing  .AND.  topography /= 'flat'                   &
                              .AND.  .NOT.  lsf_exception )  THEN
          message_string = 'The usage of large scale forcing from external &'//&
                        'file LSF_DATA is not implemented for non-flat topography'
          CALL message( 'check_parameters', 'PA0377', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( large_scale_forcing  .AND.  ocean )  THEN
          message_string = 'The usage of large scale forcing from external &'//&
                        'file LSF_DATA is not implemented for ocean runs'
          CALL message( 'check_parameters', 'PA0378', 1, 2, 0, 6, 0 )
       ENDIF

    END SUBROUTINE lsf_nudging_check_parameters

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output of profiles for land surface model
!------------------------------------------------------------------------------!
    SUBROUTINE lsf_nudging_check_data_output_pr( variable, var_count, unit,    &
                                                 dopr_unit )
 
       USE profil_parameter

       IMPLICIT NONE
   
       CHARACTER (LEN=*) ::  unit      !< 
       CHARACTER (LEN=*) ::  variable  !< 
       CHARACTER (LEN=*) ::  dopr_unit !< local value of dopr_unit
 
       INTEGER(iwp) ::  user_pr_index !<    
       INTEGER(iwp) ::  var_count     !< 

       SELECT CASE ( TRIM( variable ) )
       

          CASE ( 'td_lsa_lpt' )
             IF (  .NOT.  large_scale_forcing )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) //          &
                                 ' is not implemented for ' //                 &
                                 'large_scale_forcing = .FALSE.'
                CALL message( 'lsf_nudging_check_data_output_pr', 'PA0393',    &
                               1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 81
                dopr_unit             = 'K/s'
                unit                  = 'K/s'
                hom(:,2,81,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'td_lsa_q' )
             IF (  .NOT.  large_scale_forcing )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) //          &
                                 ' is not implemented for ' //                 &
                                 'large_scale_forcing = .FALSE.'
                CALL message( 'lsf_nudging_check_data_output_pr', 'PA0393',    &
                               1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 82
                dopr_unit             = 'kg/kgs'
                unit                  = 'kg/kgs'
                hom(:,2,82,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF
          CASE ( 'td_sub_lpt' )
             IF (  .NOT.  large_scale_forcing )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) //          &
                                 ' is not implemented for ' //                 &
                                 'large_scale_forcing = .FALSE.'
                CALL message( 'lsf_nudging_check_data_output_pr', 'PA0393',    &
                               1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 83
                dopr_unit             = 'K/s'
                unit                  = 'K/s'
                hom(:,2,83,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'td_sub_q' )
             IF (  .NOT.  large_scale_forcing )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) //          &
                                 ' is not implemented for ' //                 &
                                 'large_scale_forcing = .FALSE.'
                CALL message( 'lsf_nudging_check_data_output_pr', 'PA0393',    &
                               1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 84
                dopr_unit             = 'kg/kgs'
                unit                  = 'kg/kgs'
                hom(:,2,84,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'td_nud_lpt' )
             IF (  .NOT.  nudging )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) //          &
                                 ' is not implemented for ' //                 &
                                 'nudging = .FALSE.'
                CALL message( 'lsf_nudging_check_data_output_pr', 'PA0394',    &
                               1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 85
                dopr_unit             = 'K/s'
                unit                  = 'K/s'
                hom(:,2,85,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'td_nud_q' )
             IF (  .NOT.  nudging )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) //          &
                                 ' is not implemented for ' //                 &
                                 'nudging = .FALSE.'
                CALL message( 'lsf_nudging_check_data_output_pr', 'PA0394',    &
                               1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 86
                dopr_unit             = 'kg/kgs'
                unit                  = 'kg/kgs'
                hom(:,2,86,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'td_nud_u' )
             IF (  .NOT.  nudging )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) //          &
                                 ' is not implemented for ' //                 &
                                 'nudging = .FALSE.'
                CALL message( 'lsf_nudging_check_data_output_pr', 'PA0394',    &
                               1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 87
                dopr_unit             = 'm/s2'
                unit                  = 'm/s2'
                hom(:,2,87,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'td_nud_v' )
             IF (  .NOT.  nudging )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) //          &
                                 ' is not implemented for ' //                 &
                                 'nudging = .FALSE.'
                CALL message( 'lsf_nudging_check_data_output_pr', 'PA0394',    &
                               1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 88
                dopr_unit             = 'm/s2'
                unit                  = 'm/s2'
                hom(:,2,88,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF


          CASE DEFAULT
             unit = 'illegal'
   
       END SELECT

    END SUBROUTINE lsf_nudging_check_data_output_pr

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE lsf_nudging_header ( io )

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN) ::  io !< Unit of the output file

       WRITE ( io, 1 )
       IF ( large_scale_forcing )  THEN
          WRITE ( io, 3 )
          WRITE ( io, 4 )

          IF ( large_scale_subsidence )  THEN
             IF ( .NOT. use_subsidence_tendencies )  THEN
                WRITE ( io, 5 )
             ELSE
                WRITE ( io, 6 )
             ENDIF 
          ENDIF

          IF ( bc_pt_b == 'dirichlet' )  THEN
             WRITE ( io, 12 )
          ELSEIF ( bc_pt_b == 'neumann' )  THEN
             WRITE ( io, 13 )
          ENDIF

          IF ( bc_q_b == 'dirichlet' )  THEN
             WRITE ( io, 14 )
          ELSEIF ( bc_q_b == 'neumann' )  THEN
             WRITE ( io, 15 )
          ENDIF

          WRITE ( io, 7 )
          IF ( nudging )  THEN
             WRITE ( io, 10 )
          ENDIF
       ELSE
          WRITE ( io, 2 )
          WRITE ( io, 11 )
       ENDIF
       IF ( large_scale_subsidence )  THEN
          WRITE ( io, 8 )
          WRITE ( io, 9 )
       ENDIF


1 FORMAT (//' Large scale forcing and nudging:'/ &
              ' -------------------------------'/)
2 FORMAT (' --> No large scale forcing from external is used (default) ')
3 FORMAT (' --> Large scale forcing from external file LSF_DATA is used: ')
4 FORMAT ('     - large scale advection tendencies ')
5 FORMAT ('     - large scale subsidence velocity w_subs ')
6 FORMAT ('     - large scale subsidence tendencies ')
7 FORMAT ('     - and geostrophic wind components ug and vg')
8 FORMAT (' --> Large-scale vertical motion is used in the ', &
                  'prognostic equation(s) for')
9 FORMAT ('     the scalar(s) only')
10 FORMAT (' --> Nudging is used')
11 FORMAT (' --> No nudging is used (default) ')
12 FORMAT ('     - prescribed surface values for temperature')
13 FORMAT ('     - prescribed surface fluxes for temperature')
14 FORMAT ('     - prescribed surface values for humidity')
15 FORMAT ('     - prescribed surface fluxes for humidity')

    END SUBROUTINE lsf_nudging_header 

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE lsf_init

       USE control_parameters,                                                 &
           ONLY:  bc_lr_cyc, bc_ns_cyc

       USE netcdf_data_input_mod,                                              &
           ONLY:  netcdf_data_input_lsf  

       IMPLICIT NONE

       CHARACTER(100) ::  chmess      !<
       CHARACTER(1)   ::  hash        !<

       INTEGER(iwp) ::  ierrn         !<
       INTEGER(iwp) ::  finput = 90   !<
       INTEGER(iwp) ::  k             !<
       INTEGER(iwp) ::  nt            !<
       INTEGER(iwp) ::  t             !< running index for time levels

       REAL(wp) ::  fac               !<
       REAL(wp) ::  highheight        !<
       REAL(wp) ::  highug_vert       !<
       REAL(wp) ::  highvg_vert       !<
       REAL(wp) ::  highwsubs_vert    !<
       REAL(wp) ::  lowheight         !<
       REAL(wp) ::  lowug_vert        !<
       REAL(wp) ::  lowvg_vert        !<
       REAL(wp) ::  lowwsubs_vert     !<
       REAL(wp) ::  high_td_lsa_lpt   !<
       REAL(wp) ::  low_td_lsa_lpt    !<
       REAL(wp) ::  high_td_lsa_q     !<
       REAL(wp) ::  low_td_lsa_q      !<
       REAL(wp) ::  high_td_sub_lpt   !<
       REAL(wp) ::  low_td_sub_lpt    !<
       REAL(wp) ::  high_td_sub_q     !<
       REAL(wp) ::  low_td_sub_q      !<
       REAL(wp) ::  r_dummy           !<

       IF ( forcing )  THEN
!
!--       Allocate arrays for geostrophic wind components. Arrays will 
!--       incorporate 2 time levels in order to interpolate in between. Please
!--       note, forcing using geostrophic wind components is only required in 
!--       case of cyclic boundary conditions.
          IF ( bc_lr_cyc  .AND.  bc_ns_cyc )  THEN
             ALLOCATE( force%ug(0:1,nzb:nzt+1) )
             ALLOCATE( force%vg(0:1,nzb:nzt+1) )
          ENDIF
!
!--       Allocate arrays for reading boundary values. Arrays will incorporate 2 
!--       time levels in order to interpolate in between. 
          IF ( force_bound_l )  THEN
             ALLOCATE( force%u_left(0:1,nzb+1:nzt+1,nys:nyn)  )
             ALLOCATE( force%v_left(0:1,nzb+1:nzt+1,nysv:nyn) )
             ALLOCATE( force%w_left(0:1,nzb+1:nzt,nys:nyn)    )
             IF ( .NOT. neutral ) ALLOCATE( force%pt_left(0:1,nzb+1:nzt+1,nys:nyn) )
          ENDIF
          IF ( force_bound_r )  THEN
             ALLOCATE( force%u_right(0:1,nzb+1:nzt+1,nys:nyn)  )
             ALLOCATE( force%v_right(0:1,nzb+1:nzt+1,nysv:nyn) )
             ALLOCATE( force%w_right(0:1,nzb+1:nzt,nys:nyn)    )
             IF ( .NOT. neutral ) ALLOCATE( force%pt_right(0:1,nzb+1:nzt+1,nys:nyn) )
          ENDIF
          IF ( force_bound_n )  THEN
             ALLOCATE( force%u_north(0:1,nzb+1:nzt+1,nxlu:nxr) )
             ALLOCATE( force%v_north(0:1,nzb+1:nzt+1,nxl:nxr)  )
             ALLOCATE( force%w_north(0:1,nzb+1:nzt,nxl:nxr)    )
             IF ( .NOT. neutral ) ALLOCATE( force%pt_north(0:1,nzb+1:nzt+1,nxl:nxr) )
          ENDIF
          IF ( force_bound_s )  THEN
             ALLOCATE( force%u_south(0:1,nzb+1:nzt+1,nxlu:nxr) )
             ALLOCATE( force%v_south(0:1,nzb+1:nzt+1,nxl:nxr)  )
             ALLOCATE( force%w_south(0:1,nzb+1:nzt,nxl:nxr)    )
             IF ( .NOT. neutral ) ALLOCATE( force%pt_south(0:1,nzb+1:nzt+1,nxl:nxr) )
          ENDIF
          
          ALLOCATE( force%u_top(0:1,nys:nyn,nxlu:nxr) )
          ALLOCATE( force%v_top(0:1,nysv:nyn,nxl:nxr) )
          ALLOCATE( force%w_top(0:1,nys:nyn,nxl:nxr)  )
          IF ( .NOT. neutral ) ALLOCATE( force%pt_top(0:1,nys:nyn,nxl:nxr) )

!
!--       Initial call of input. Time array, initial 3D data of u, v, w, 
!--       potential temperature, as well as mixing ratio, will be read. 
!--       Moreover, data at lateral and top boundary will be read.  
          CALL netcdf_data_input_lsf
!
!--       Please note, at the moment INIFOR assumes only an equidistant vertical
!--       grid. In case of vertical grid stretching, input of inital 3D data 
!--       need to be inter- and/or extrapolated. 
!--       Therefore, check if zw grid on file is identical to numeric zw grid.   
          IF ( ANY( zu(1:nzt+1) /= force%zu_atmos(1:force%nzu) ) )  THEN
!
!--          Also data at the boundaries need to be inter/extrapolated at both
!--          time levels
             DO  t = 0, 1
                IF ( force_bound_l )  THEN
                   CALL netcdf_data_input_interpolate( force%u_left(t,:,:),    &
                                                       zu(1:nzt+1),            &
                                                       force%zu_atmos )
                   CALL netcdf_data_input_interpolate( force%v_left(t,:,:),    &
                                                       zu(1:nzt+1),            &
                                                       force%zu_atmos )
                   CALL netcdf_data_input_interpolate( force%w_left(t,:,:),    &
                                                       zw(1:nzt+1),            &
                                                       force%zw_atmos )
                   IF ( .NOT. neutral )                                        &
                      CALL netcdf_data_input_interpolate( force%pt_left(t,:,:),&
                                                          zu(1:nzt+1),         &
                                                          force%zu_atmos )
                ENDIF
                IF ( force_bound_r )  THEN
                   CALL netcdf_data_input_interpolate( force%u_right(t,:,:),   &
                                                       zu(1:nzt+1),            &
                                                       force%zu_atmos )
                   CALL netcdf_data_input_interpolate( force%v_right(t,:,:),   &
                                                       zu(1:nzt+1),            &
                                                       force%zu_atmos )
                   CALL netcdf_data_input_interpolate( force%w_right(t,:,:),   &
                                                       zw(1:nzt+1),            &
                                                       force%zw_atmos )
                   IF ( .NOT. neutral )                                        &
                      CALL netcdf_data_input_interpolate( force%pt_right(t,:,:),&
                                                          zu(1:nzt+1),          &
                                                          force%zu_atmos )
                ENDIF
                IF ( force_bound_n )  THEN
                   CALL netcdf_data_input_interpolate( force%u_north(t,:,:),   &
                                                       zu(1:nzt+1),            &
                                                       force%zu_atmos )
                   CALL netcdf_data_input_interpolate( force%v_north(t,:,:),   &
                                                       zu(1:nzt+1),            &
                                                       force%zu_atmos )
                   CALL netcdf_data_input_interpolate( force%w_north(t,:,:),   &
                                                       zw(1:nzt+1),            &
                                                       force%zw_atmos )
                   IF ( .NOT. neutral )                                        &
                      CALL netcdf_data_input_interpolate( force%pt_north(t,:,:),&
                                                          zu(1:nzt+1),          &
                                                          force%zu_atmos )
                ENDIF
                IF ( force_bound_s )  THEN
                   CALL netcdf_data_input_interpolate( force%u_south(t,:,:),   &
                                                       zu(1:nzt+1),            &
                                                       force%zu_atmos )
                   CALL netcdf_data_input_interpolate( force%v_south(t,:,:),   &
                                                       zu(1:nzt+1),            &
                                                       force%zu_atmos )
                   CALL netcdf_data_input_interpolate( force%w_south(t,:,:),   &
                                                       zw(1:nzt+1),            &
                                                       force%zw_atmos )
                   IF ( .NOT. neutral )                                        &
                      CALL netcdf_data_input_interpolate( force%pt_south(t,:,:),&
                                                          zu(1:nzt+1),          &
                                                          force%zu_atmos )
                ENDIF
             ENDDO
          ENDIF

!
!--       Exchange ghost points
          CALL exchange_horiz( u, nbgp )
          CALL exchange_horiz( v, nbgp )
          CALL exchange_horiz( w, nbgp )
          IF ( .NOT. neutral )  CALL exchange_horiz( pt, nbgp )
!
!--       At lateral boundaries, set also initial boundary conditions
          IF ( force_bound_l )  THEN
             u(:,:,nxl)   = u(:,:,nxlu)
             v(:,:,nxl-1) = v(:,:,nxl)
             w(:,:,nxl-1) = w(:,:,nxl)
             IF ( .NOT. neutral )  pt(:,:,nxl-1) = pt(:,:,nxl)
          ENDIF
          IF ( force_bound_r )  THEN
             u(:,:,nxr+1) = u(:,:,nxr)
             v(:,:,nxr+1) = v(:,:,nxr)
             w(:,:,nxr+1) = w(:,:,nxr)
             IF ( .NOT. neutral )  pt(:,:,nxr+1) = pt(:,:,nxr)
          ENDIF
          IF ( force_bound_s )  THEN
             u(:,nys-1,:) = u(:,nys,:)
             v(:,nys,:)   = v(:,nysv,:)
             w(:,nys-1,:) = w(:,nys,:)
             IF ( .NOT. neutral )  pt(:,nys-1,:) = pt(:,nys,:)
          ENDIF
          IF ( force_bound_n )  THEN
             u(:,nyn+1,:) = u(:,nyn,:)
             v(:,nyn+1,:) = v(:,nyn,:)
             w(:,nyn+1,:) = w(:,nyn,:)
             IF ( .NOT. neutral )  pt(:,nyn+1,:) = pt(:,nyn,:)
          ENDIF

!
!--       After 3D data is initialized, ensure mass conservation
          CALL forcing_bc_mass_conservation
!
!--       Initialize surface pressure. Please note, time-dependent surface
!--       pressure would require changes in anelastic approximation and 
!--       treatment of fluxes. 
!--       For the moment, comment this out! 
!         surface_pressure = force%surface_pressure(0)

       ELSE

          ALLOCATE( p_surf(0:nlsf), pt_surf(0:nlsf), q_surf(0:nlsf),           &
                    qsws_surf(0:nlsf), shf_surf(0:nlsf),                       &
                    td_lsa_lpt(nzb:nzt+1,0:nlsf), td_lsa_q(nzb:nzt+1,0:nlsf),  &
                    td_sub_lpt(nzb:nzt+1,0:nlsf), td_sub_q(nzb:nzt+1,0:nlsf),  &
                    time_vert(0:nlsf), time_surf(0:nlsf),                      &
                    ug_vert(nzb:nzt+1,0:nlsf), vg_vert(nzb:nzt+1,0:nlsf),      &
                    wsubs_vert(nzb:nzt+1,0:nlsf) )

          p_surf = 0.0_wp; pt_surf = 0.0_wp; q_surf = 0.0_wp; qsws_surf = 0.0_wp
          shf_surf = 0.0_wp; time_vert = 0.0_wp; td_lsa_lpt = 0.0_wp
          td_lsa_q = 0.0_wp; td_sub_lpt = 0.0_wp; td_sub_q = 0.0_wp
          time_surf = 0.0_wp; ug_vert = 0.0_wp; vg_vert = 0.0_wp
          wsubs_vert = 0.0_wp

!
!--       Array for storing large scale forcing and nudging tendencies at each 
!--       timestep for data output
          ALLOCATE( sums_ls_l(nzb:nzt+1,0:7) )
          sums_ls_l = 0.0_wp

          ngp_sums_ls = (nz+2)*6

          OPEN ( finput, FILE='LSF_DATA', STATUS='OLD', &
                 FORM='FORMATTED', IOSTAT=ierrn )

          IF ( ierrn /= 0 )  THEN
             message_string = 'file LSF_DATA does not exist'
             CALL message( 'ls_forcing', 'PA0368', 1, 2, 0, 6, 0 )
          ENDIF

          ierrn = 0
!
!--       First three lines of LSF_DATA contain header
          READ ( finput, FMT='(a100)', IOSTAT=ierrn ) chmess
          READ ( finput, FMT='(a100)', IOSTAT=ierrn ) chmess
          READ ( finput, FMT='(a100)', IOSTAT=ierrn ) chmess

          IF ( ierrn /= 0 )  THEN
             message_string = 'errors in file LSF_DATA'
             CALL message( 'ls_forcing', 'PA0369', 1, 2, 0, 6, 0 )
          ENDIF

!
!--       Surface values are read in
          nt     = 0
          ierrn = 0

          DO WHILE ( time_surf(nt) < end_time )
             nt = nt + 1
             READ ( finput, *, IOSTAT = ierrn ) time_surf(nt), shf_surf(nt),   &
                                                qsws_surf(nt), pt_surf(nt),    &
                                                q_surf(nt), p_surf(nt)

             IF ( ierrn /= 0 )  THEN
               WRITE ( message_string, * ) 'No time dependent surface ' //     &
                                 'variables in & LSF_DATA for end of run found'

                CALL message( 'ls_forcing', 'PA0363', 1, 2, 0, 6, 0 )
             ENDIF
          ENDDO

          IF ( time_surf(1) > end_time )  THEN
             WRITE ( message_string, * ) 'Time dependent surface variables in ' // &
                                         '&LSF_DATA set in after end of ' ,        &
                                         'simulation - lsf_surf is set to FALSE'
             CALL message( 'ls_forcing', 'PA0371', 0, 0, 0, 6, 0 )
             lsf_surf = .FALSE.
          ENDIF

!
!--       Go to the end of the list with surface variables
          DO WHILE ( ierrn == 0 )
             READ ( finput, *, IOSTAT = ierrn ) r_dummy
          ENDDO

!
!--       Profiles of ug, vg and w_subs are read in (large scale forcing)

          nt = 0
          DO WHILE ( time_vert(nt) < end_time )
             nt = nt + 1
             hash = "#"
             ierrn = 1 ! not zero
!
!--          Search for the next line consisting of "# time", 
!--          from there onwards the profiles will be read
             DO WHILE ( .NOT. ( hash == "#" .AND. ierrn == 0 ) ) 
                READ ( finput, *, IOSTAT=ierrn ) hash, time_vert(nt)
                IF ( ierrn < 0 )  THEN 
                   WRITE( message_string, * ) 'No time dependent vertical profiles',&
                                    ' in & LSF_DATA for end of run found'
                   CALL message( 'ls_forcing', 'PA0372', 1, 2, 0, 6, 0 )
                ENDIF
             ENDDO

             IF ( nt == 1 .AND. time_vert(nt) > end_time ) EXIT

             READ ( finput, *, IOSTAT=ierrn ) lowheight, lowug_vert, lowvg_vert,&
                                              lowwsubs_vert, low_td_lsa_lpt,    &
                                              low_td_lsa_q, low_td_sub_lpt,     &
                                              low_td_sub_q
             IF ( ierrn /= 0 )  THEN
                message_string = 'errors in file LSF_DATA'
                CALL message( 'ls_forcing', 'PA0369', 1, 2, 0, 6, 0 )
             ENDIF

             READ ( finput, *, IOSTAT=ierrn ) highheight, highug_vert,         &
                                              highvg_vert, highwsubs_vert,     &
                                              high_td_lsa_lpt, high_td_lsa_q,  &
                                              high_td_sub_lpt, high_td_sub_q
         
             IF ( ierrn /= 0 )  THEN
                message_string = 'errors in file LSF_DATA'
                CALL message( 'ls_forcing', 'PA0369', 1, 2, 0, 6, 0 )
             ENDIF


             DO  k = nzb, nzt+1
                IF ( highheight < zu(k) )  THEN
                   lowheight      = highheight
                   lowug_vert     = highug_vert
                   lowvg_vert     = highvg_vert
                   lowwsubs_vert  = highwsubs_vert
                   low_td_lsa_lpt = high_td_lsa_lpt
                   low_td_lsa_q   = high_td_lsa_q
                   low_td_sub_lpt = high_td_sub_lpt
                   low_td_sub_q   = high_td_sub_q

                   ierrn = 0
                   READ ( finput, *, IOSTAT=ierrn ) highheight, highug_vert,    &
                                                    highvg_vert, highwsubs_vert,&
                                                    high_td_lsa_lpt,            &
                                                    high_td_lsa_q,              &
                                                    high_td_sub_lpt, high_td_sub_q

                   IF ( ierrn /= 0 )  THEN
                      WRITE( message_string, * ) 'zu(',k,') = ', zu(k), 'm ',  &
                           'is higher than the maximum height in LSF_DATA ',   &
                           'which is ', lowheight, 'm. Interpolation on PALM ',&
                           'grid is not possible.'
                      CALL message( 'ls_forcing', 'PA0395', 1, 2, 0, 6, 0 )
                   ENDIF

                ENDIF

!
!--             Interpolation of prescribed profiles in space 
                fac = (highheight-zu(k))/(highheight - lowheight)

                ug_vert(k,nt)    = fac * lowug_vert                            &
                                   + ( 1.0_wp - fac ) * highug_vert
                vg_vert(k,nt)    = fac * lowvg_vert                            &
                                   + ( 1.0_wp - fac ) * highvg_vert
                wsubs_vert(k,nt) = fac * lowwsubs_vert                         &
                                   + ( 1.0_wp - fac ) * highwsubs_vert

                td_lsa_lpt(k,nt) = fac * low_td_lsa_lpt                        &
                                   + ( 1.0_wp - fac ) * high_td_lsa_lpt
                td_lsa_q(k,nt)   = fac * low_td_lsa_q                          &
                                   + ( 1.0_wp - fac ) * high_td_lsa_q
                td_sub_lpt(k,nt) = fac * low_td_sub_lpt                        &
                                   + ( 1.0_wp - fac ) * high_td_sub_lpt
                td_sub_q(k,nt)   = fac * low_td_sub_q                          &
                                   + ( 1.0_wp - fac ) * high_td_sub_q

             ENDDO

          ENDDO 

!
!--       Large scale vertical velocity has to be zero at the surface
          wsubs_vert(nzb,:) = 0.0_wp
    
          IF ( time_vert(1) > end_time )  THEN
             WRITE ( message_string, * ) 'Time dependent large scale profile ',&
                                'forcing from&LSF_DATA sets in after end of ' ,&
                                'simulation - lsf_vert is set to FALSE'
             CALL message( 'ls_forcing', 'PA0373', 0, 0, 0, 6, 0 )
             lsf_vert = .FALSE.
          ENDIF

          CLOSE( finput )

       ENDIF

    END SUBROUTINE lsf_init

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE ls_forcing_surf ( time )

       IMPLICIT NONE

       INTEGER(iwp) ::  nt                     !<

       REAL(wp)             :: dum_surf_flux  !<
       REAL(wp)             :: fac            !<
       REAL(wp), INTENT(in) :: time           !<

!
!--    Interpolation in time of LSF_DATA at the surface
       nt = 1
       DO WHILE ( time > time_surf(nt) )
          nt = nt + 1
       ENDDO
       IF ( time /= time_surf(nt) )  THEN
         nt = nt - 1
       ENDIF

       fac = ( time -time_surf(nt) ) / ( time_surf(nt+1) - time_surf(nt) )

       IF ( ibc_pt_b == 0 )  THEN
!
!--       In case of Dirichlet boundary condition shf must not 
!--       be set - it is calculated via MOST in prandtl_fluxes
          pt_surface = pt_surf(nt) + fac * ( pt_surf(nt+1) - pt_surf(nt) )

       ELSEIF ( ibc_pt_b == 1 )  THEN
!
!--       In case of Neumann boundary condition pt_surface is needed for 
!--       calculation of reference density
          dum_surf_flux = ( shf_surf(nt) + fac *                               &
                            ( shf_surf(nt+1) - shf_surf(nt) )                  &
                          ) * heatflux_input_conversion(nzb)
!
!--       Save surface sensible heat flux on default, natural and urban surface
!--       type, if required 
          IF ( surf_def_h(0)%ns >= 1 )  surf_def_h(0)%shf(:) = dum_surf_flux
          IF ( surf_lsm_h%ns    >= 1 )  surf_lsm_h%shf(:)    = dum_surf_flux
          IF ( surf_usm_h%ns    >= 1 )  surf_usm_h%shf(:)    = dum_surf_flux

          pt_surface    = pt_surf(nt) + fac * ( pt_surf(nt+1) - pt_surf(nt) )

       ENDIF

       IF ( ibc_q_b == 0 )  THEN
!
!--       In case of Dirichlet boundary condition qsws must not 
!--       be set - it is calculated via MOST in prandtl_fluxes
          q_surface = q_surf(nt) + fac * ( q_surf(nt+1) - q_surf(nt) )

       ELSEIF ( ibc_q_b == 1 )  THEN
          dum_surf_flux = ( qsws_surf(nt) + fac *                              &
                             ( qsws_surf(nt+1) - qsws_surf(nt) )               &
                             ) * waterflux_input_conversion(nzb)
!
!--       Save surface sensible heat flux on default, natural and urban surface
!--       type, if required 
          IF ( surf_def_h(0)%ns >= 1 )  surf_def_h(0)%qsws(:) = dum_surf_flux
          IF ( surf_lsm_h%ns    >= 1 )  surf_lsm_h%qsws(:)    = dum_surf_flux
          IF ( surf_usm_h%ns    >= 1 )  surf_usm_h%qsws(:)    = dum_surf_flux

       ENDIF
!
!--    Surface heat- and waterflux will be written later onto surface elements 
       IF ( .NOT.  neutral  .AND.  constant_heatflux  .AND.                    &
            TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
             surface_heatflux = shf_surf(1)
       ENDIF

       surface_pressure = p_surf(nt) + fac * ( p_surf(nt+1) - p_surf(nt) )

    END SUBROUTINE ls_forcing_surf 




!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE ls_forcing_vert ( time )


       IMPLICIT NONE

       INTEGER(iwp) ::  nt                     !<

       REAL(wp)             ::  fac           !<
       REAL(wp), INTENT(in) ::  time          !<

!
!--    Interpolation in time of LSF_DATA for ug, vg and w_subs
       nt = 1
       DO WHILE ( time > time_vert(nt) )
          nt = nt + 1
       ENDDO
       IF ( time /= time_vert(nt) )  THEN
         nt = nt - 1
       ENDIF

       fac = ( time-time_vert(nt) ) / ( time_vert(nt+1)-time_vert(nt) )

       ug     = ug_vert(:,nt) + fac * ( ug_vert(:,nt+1) - ug_vert(:,nt) )
       vg     = vg_vert(:,nt) + fac * ( vg_vert(:,nt+1) - vg_vert(:,nt) )

    END SUBROUTINE ls_forcing_vert


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE ls_advec ( time, prog_var )
      

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  prog_var   !< 

       REAL(wp), INTENT(in)  :: time    !< 
       REAL(wp) :: fac                  !<  

       INTEGER(iwp) ::  i               !< 
       INTEGER(iwp) ::  j               !< 
       INTEGER(iwp) ::  k               !< 
       INTEGER(iwp) ::  nt               !< 

!
!--    Interpolation in time of LSF_DATA 
       nt = 1
       DO WHILE ( time > time_vert(nt) )
          nt = nt + 1
       ENDDO
       IF ( time /= time_vert(nt) )  THEN
         nt = nt - 1
       ENDIF

       fac = ( time-time_vert(nt) ) / ( time_vert(nt+1)-time_vert(nt) )

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                tend(k,j,i) = tend(k,j,i) + td_lsa_lpt(k,nt) + fac *     &
                              ( td_lsa_lpt(k,nt+1) - td_lsa_lpt(k,nt) ) *&
                                  MERGE( 1.0_wp, 0.0_wp,                 &
                                         BTEST( wall_flags_0(k,j,i), 0 ) )
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE ls_advec


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE ls_advec_ij ( i, j, time, prog_var )

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  prog_var   !< 

       REAL(wp), INTENT(in)  :: time    !< 
       REAL(wp) :: fac                  !< 

       INTEGER(iwp) ::  i               !< 
       INTEGER(iwp) ::  j               !< 
       INTEGER(iwp) ::  k               !< 
       INTEGER(iwp) ::  nt               !< 

!
!--    Interpolation in time of LSF_DATA 
       nt = 1
       DO WHILE ( time > time_vert(nt) )
          nt = nt + 1
       ENDDO
       IF ( time /= time_vert(nt) )  THEN
         nt = nt - 1
       ENDIF

       fac = ( time-time_vert(nt) ) / ( time_vert(nt+1)-time_vert(nt) )


       DO  k = nzb+1, nzt
          tend(k,j,i) = tend(k,j,i) + td_lsa_lpt(k,nt)                   &
                       + fac * ( td_lsa_lpt(k,nt+1) - td_lsa_lpt(k,nt) )*&
                                  MERGE( 1.0_wp, 0.0_wp,                 &
                                         BTEST( wall_flags_0(k,j,i), 0 ) )
       ENDDO

    END SUBROUTINE ls_advec_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE nudge_init

       IMPLICIT NONE


       INTEGER(iwp) ::  finput = 90  !<
       INTEGER(iwp) ::  ierrn        !<
       INTEGER(iwp) ::  k            !<
       INTEGER(iwp) ::  nt            !<

       CHARACTER(1) ::  hash     !<

       REAL(wp) ::  highheight   !<
       REAL(wp) ::  highqnudge   !<
       REAL(wp) ::  highptnudge  !<
       REAL(wp) ::  highunudge   !<
       REAL(wp) ::  highvnudge   !<
       REAL(wp) ::  highwnudge   !<
       REAL(wp) ::  hightnudge   !<

       REAL(wp) ::  lowheight    !<
       REAL(wp) ::  lowqnudge    !<
       REAL(wp) ::  lowptnudge   !<
       REAL(wp) ::  lowunudge    !<
       REAL(wp) ::  lowvnudge    !<
       REAL(wp) ::  lowwnudge    !<
       REAL(wp) ::  lowtnudge    !<

       REAL(wp) ::  fac          !<

       ALLOCATE( ptnudge(nzb:nzt+1,1:ntnudge), qnudge(nzb:nzt+1,1:ntnudge), &
                 tnudge(nzb:nzt+1,1:ntnudge), unudge(nzb:nzt+1,1:ntnudge),  &
                 vnudge(nzb:nzt+1,1:ntnudge), wnudge(nzb:nzt+1,1:ntnudge)  )

       ALLOCATE( tmp_tnudge(nzb:nzt) )

       ALLOCATE( timenudge(0:ntnudge) )

       ptnudge = 0.0_wp; qnudge = 0.0_wp; tnudge = 0.0_wp; unudge = 0.0_wp
       vnudge = 0.0_wp; wnudge = 0.0_wp; timenudge = 0.0_wp
!
!--    Initialize array tmp_nudge with a current nudging time scale of 6 hours
       tmp_tnudge = 21600.0_wp

       nt = 0
       OPEN ( finput, FILE='NUDGING_DATA', STATUS='OLD', &
              FORM='FORMATTED', IOSTAT=ierrn )

       IF ( ierrn /= 0 )  THEN
          message_string = 'file NUDGING_DATA does not exist'
          CALL message( 'nudging', 'PA0365', 1, 2, 0, 6, 0 )
       ENDIF

       ierrn = 0

 rloop:DO
          nt = nt + 1
          hash = "#"
          ierrn = 1 ! not zero
!
!--       Search for the next line consisting of "# time", 
!--       from there onwards the profiles will be read
          DO WHILE ( .NOT. ( hash == "#" .AND. ierrn == 0 ) ) 
          
            READ ( finput, *, IOSTAT=ierrn ) hash, timenudge(nt)
            IF ( ierrn < 0 )  EXIT rloop

          ENDDO

          ierrn = 0
          READ ( finput, *, IOSTAT=ierrn ) lowheight, lowtnudge, lowunudge,   &
                                           lowvnudge, lowwnudge , lowptnudge, &
                                           lowqnudge

          IF ( ierrn /= 0 )  THEN
             message_string = 'errors in file NUDGING_DATA'
             CALL message( 'nudging', 'PA0366', 1, 2, 0, 6, 0 )
          ENDIF

          ierrn = 0
          READ ( finput, *, IOSTAT=ierrn ) highheight, hightnudge, highunudge,   &
                                           highvnudge, highwnudge , highptnudge, &
                                           highqnudge

          IF ( ierrn /= 0 )  THEN
             message_string = 'errors in file NUDGING_DATA'
             CALL message( 'nudging', 'PA0366', 1, 2, 0, 6, 0 )
          ENDIF

          DO  k = nzb, nzt+1
             DO WHILE ( highheight < zu(k) )
                lowheight  = highheight
                lowtnudge  = hightnudge
                lowunudge  = highunudge
                lowvnudge  = highvnudge
                lowwnudge  = highwnudge
                lowptnudge = highptnudge
                lowqnudge  = highqnudge
 
                ierrn = 0
                READ ( finput, *, IOSTAT=ierrn )  highheight , hightnudge ,    &
                                                  highunudge , highvnudge ,    &
                                                  highwnudge , highptnudge,    &
                                                  highqnudge
                IF (ierrn /= 0 )  THEN
                   WRITE( message_string, * ) 'zu(',k,') = ', zu(k), 'm is ',  &
                        'higher than the maximum height in NUDING_DATA which ',&
                        'is ', lowheight, 'm. Interpolation on PALM ',         &
                        'grid is not possible.'
                   CALL message( 'nudging', 'PA0364', 1, 2, 0, 6, 0 )
                ENDIF
             ENDDO

!
!--          Interpolation of prescribed profiles in space 

             fac = ( highheight - zu(k) ) / ( highheight - lowheight )

             tnudge(k,nt)  = fac * lowtnudge  + ( 1.0_wp - fac ) * hightnudge
             unudge(k,nt)  = fac * lowunudge  + ( 1.0_wp - fac ) * highunudge
             vnudge(k,nt)  = fac * lowvnudge  + ( 1.0_wp - fac ) * highvnudge
             wnudge(k,nt)  = fac * lowwnudge  + ( 1.0_wp - fac ) * highwnudge
             ptnudge(k,nt) = fac * lowptnudge + ( 1.0_wp - fac ) * highptnudge
             qnudge(k,nt)  = fac * lowqnudge  + ( 1.0_wp - fac ) * highqnudge
          ENDDO

       ENDDO rloop

       CLOSE ( finput )

!
!--    Overwrite initial profiles in case of nudging
       IF ( nudging )  THEN
          pt_init = ptnudge(:,1)
          u_init  = unudge(:,1)
          v_init  = vnudge(:,1)

          WRITE( message_string, * ) 'Initial profiles of u, v, pt and q ',    &
                                     'from NUDGING_DATA are used.'
          CALL message( 'large_scale_forcing_nudging', 'PA0370', 0, 0, 0, 6, 0 )
       ENDIF


    END SUBROUTINE nudge_init

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_tnudge ( time )

       IMPLICIT NONE


       REAL(wp) ::  dtm         !<
       REAL(wp) ::  dtp         !<
       REAL(wp) ::  time        !<

       INTEGER(iwp) ::  k   !<
       INTEGER(iwp) ::  nt  !<

       nt = 1
       DO WHILE ( time > timenudge(nt) )
         nt = nt+1
       ENDDO
       IF ( time /= timenudge(1) ) THEN
         nt = nt-1
       ENDIF

       dtm = ( time - timenudge(nt) ) / ( timenudge(nt+1) - timenudge(nt) )
       dtp = ( timenudge(nt+1) - time ) / ( timenudge(nt+1) - timenudge(nt) )

       DO  k = nzb, nzt
          tmp_tnudge(k) = MAX( dt_3d, tnudge(k,nt) * dtp + tnudge(k,nt+1) * dtm )
       ENDDO

    END SUBROUTINE calc_tnudge

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE nudge ( time, prog_var )

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  prog_var  !<

       REAL(wp) ::  tmp_tend    !<
       REAL(wp) ::  dtm         !<
       REAL(wp) ::  dtp         !<
       REAL(wp) ::  time        !<

       INTEGER(iwp) ::  i  !<
       INTEGER(iwp) ::  j  !<
       INTEGER(iwp) ::  k  !<
       INTEGER(iwp) ::  nt  !<


       nt = 1
       DO WHILE ( time > timenudge(nt) )
         nt = nt+1
       ENDDO
       IF ( time /= timenudge(1) ) THEN
         nt = nt-1
       ENDIF

       dtm = ( time - timenudge(nt) ) / ( timenudge(nt+1) - timenudge(nt) )
       dtp = ( timenudge(nt+1) - time ) / ( timenudge(nt+1) - timenudge(nt) )

       SELECT CASE ( prog_var )

          CASE ( 'u' )

             DO  i = nxl, nxr
                DO  j = nys, nyn

                   DO  k = nzb+1, nzt

                      tmp_tend = - ( hom(k,1,1,0) - ( unudge(k,nt) * dtp +     &
                                     unudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                      tend(k,j,i) = tend(k,j,i) + tmp_tend *                   &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 1 ) )

                      sums_ls_l(k,6) = sums_ls_l(k,6) + tmp_tend *             &
                                     weight_substep(intermediate_timestep_count)
                   ENDDO
                  
                   sums_ls_l(nzt+1,6) = sums_ls_l(nzt,6)
 
                ENDDO
            ENDDO

          CASE ( 'v' )

             DO  i = nxl, nxr
                DO  j = nys, nyn

                   DO  k = nzb+1, nzt

                      tmp_tend = - ( hom(k,1,2,0) - ( vnudge(k,nt) * dtp +     &
                                     vnudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                      tend(k,j,i) = tend(k,j,i) + tmp_tend *                   &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 2 ) )

                      sums_ls_l(k,7) = sums_ls_l(k,7) + tmp_tend *             &
                                     weight_substep(intermediate_timestep_count)
                   ENDDO
                  
                   sums_ls_l(nzt+1,7) = sums_ls_l(nzt,7)

                ENDDO
            ENDDO

          CASE ( 'pt' )

             DO  i = nxl, nxr
                DO  j = nys, nyn

                   DO  k = nzb+1, nzt

                      tmp_tend = - ( hom(k,1,4,0) - ( ptnudge(k,nt) * dtp +    &
                                     ptnudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                      tend(k,j,i) = tend(k,j,i) + tmp_tend *                   &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 0 ) )

                      sums_ls_l(k,4) = sums_ls_l(k,4) + tmp_tend *             &
                                     weight_substep(intermediate_timestep_count)
                   ENDDO

                   sums_ls_l(nzt+1,4) = sums_ls_l(nzt,4)

                ENDDO
            ENDDO

          CASE DEFAULT
             message_string = 'unknown prognostic variable "' // prog_var // '"'
             CALL message( 'nudge', 'PA0367', 1, 2, 0, 6, 0 )

       END SELECT

    END SUBROUTINE nudge


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!

    SUBROUTINE nudge_ij( i, j, time, prog_var )

       IMPLICIT NONE


       CHARACTER (LEN=*) ::  prog_var  !<

       REAL(wp) ::  tmp_tend    !<
       REAL(wp) ::  dtm         !<
       REAL(wp) ::  dtp         !<
       REAL(wp) ::  time        !<

       INTEGER(iwp) ::  i  !<
       INTEGER(iwp) ::  j  !<
       INTEGER(iwp) ::  k  !<
       INTEGER(iwp) ::  nt  !<


       nt = 1
       DO WHILE ( time > timenudge(nt) )
         nt = nt+1
       ENDDO
       IF ( time /= timenudge(1) )  THEN
         nt = nt-1
       ENDIF

       dtm = ( time - timenudge(nt) ) / ( timenudge(nt+1) - timenudge(nt) )
       dtp = ( timenudge(nt+1) - time ) / ( timenudge(nt+1) - timenudge(nt) )

       SELECT CASE ( prog_var )

          CASE ( 'u' )

             DO  k = nzb+1, nzt

                tmp_tend = - ( hom(k,1,1,0) - ( unudge(k,nt) * dtp +           &
                               unudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                tend(k,j,i) = tend(k,j,i) + tmp_tend *                         &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 1 ) )

                sums_ls_l(k,6) = sums_ls_l(k,6) + tmp_tend                     &
                                 * weight_substep(intermediate_timestep_count)
             ENDDO

             sums_ls_l(nzt+1,6) = sums_ls_l(nzt,6)

          CASE ( 'v' )

             DO  k = nzb+1, nzt

                tmp_tend = - ( hom(k,1,2,0) - ( vnudge(k,nt) * dtp +           &
                               vnudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                tend(k,j,i) = tend(k,j,i) + tmp_tend *                         &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 2 ) )

                sums_ls_l(k,7) = sums_ls_l(k,7) + tmp_tend                     &
                                 * weight_substep(intermediate_timestep_count)
             ENDDO

             sums_ls_l(nzt+1,7) = sums_ls_l(nzt,7)

          CASE ( 'pt' )

             DO  k = nzb+1, nzt

                tmp_tend = - ( hom(k,1,4,0) - ( ptnudge(k,nt) * dtp +          &
                               ptnudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                tend(k,j,i) = tend(k,j,i) + tmp_tend *                         &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 0 ) )

                sums_ls_l(k,4) = sums_ls_l(k,4) + tmp_tend                     &
                                 * weight_substep(intermediate_timestep_count)
             ENDDO

             sums_ls_l(nzt+1,4) = sums_ls_l(nzt,4)

          CASE DEFAULT
             message_string = 'unknown prognostic variable "' // prog_var // '"'
             CALL message( 'nudge', 'PA0367', 1, 2, 0, 6, 0 )

       END SELECT


    END SUBROUTINE nudge_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE nudge_ref ( time )

       IMPLICIT NONE

       INTEGER(iwp) ::  nt                    !<

       REAL(wp)             ::  fac           !<
       REAL(wp), INTENT(in) ::  time          !<

!
!--    Interpolation in time of NUDGING_DATA for pt_init and q_init. This is 
!--    needed for correct upper boundary conditions for pt and q and in case that 
!      large scale subsidence as well as scalar Rayleigh-damping are used
       nt = 1
       DO WHILE ( time > time_vert(nt) )
          nt = nt + 1
       ENDDO
       IF ( time /= time_vert(nt) )  THEN
        nt = nt - 1
       ENDIF

       fac = ( time-time_vert(nt) ) / ( time_vert(nt+1)-time_vert(nt) )

       pt_init = ptnudge(:,nt) + fac * ( ptnudge(:,nt+1) - ptnudge(:,nt) )
       u_init  = unudge(:,nt) + fac * ( unudge(:,nt+1) - unudge(:,nt) )
       v_init  = vnudge(:,nt) + fac * ( vnudge(:,nt+1) - vnudge(:,nt) )

    END SUBROUTINE nudge_ref


 END MODULE lsf_nudging_mod
