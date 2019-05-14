!> @file time_integration.f90
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
! $Id: time_integration.f90 3042 2018-05-25 10:44:37Z schwenkel $
! Changed the name specific humidity to mixing ratio
!
! 3040 2018-05-25 10:22:08Z schwenkel
! Fixed bug in IF statement
! Ensure that the time when calling the radiation to be the time step of the
! pre-calculated time when first calculate the positions of the sun
!
! 3004 2018-04-27 12:33:25Z Giersch
! First call of flow_statistics has been removed. It is already called in
! run_control itself
!
! 2984 2018-04-18 11:51:30Z hellstea
! CALL pmci_ensure_nest_mass_conservation is removed (so far only commented out)
! as seemingly unnecessary.
!
! 2941 2018-04-03 11:54:58Z kanani
! Deduct spinup_time from RUN_CONTROL output of main 3d run
! (use time_since_reference_point instead of simulated_time)
!
! 2938 2018-03-27 15:52:42Z suehring
! Nesting of dissipation rate in case of RANS mode and TKE-e closure is applied
!
! 2936 2018-03-27 14:49:27Z suehring
! Little formatting adjustment.
!
! 2817 2018-02-19 16:32:21Z knoop
! Preliminary gust module interface implemented
!
! 2801 2018-02-14 16:01:55Z thiele
! Changed lpm from subroutine to module.
! Introduce particle transfer in nested models.
!
! 2776 2018-01-31 10:44:42Z Giersch
! Variable use_synthetic_turbulence_generator has been abbreviated
!
! 2773 2018-01-30 14:12:54Z suehring
! - Nesting for chemical species
!
! 2766 2018-01-22 17:17:47Z kanani
! Removed preprocessor directive __chem
!
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! - Change in file header (GPL part)
! - Implementation of uv exposure model (FK)
! - Moved vnest_boundary_conds_khkm from tcm_diffusivities to here (TG)
! - renamed diffusivities to tcm_diffusivities (TG)
! - implement prognostic equation for diss (TG)
! - Moved/commented CALL to chem_emissions (FK)
! - Added CALL to chem_emissions (FK)
! - Implementation of chemistry module (FK)
! - Calls for setting boundary conditions in USM and LSM (MS)
! - Large-scale forcing with larger-scale models implemented (MS)
! - Rename usm_radiation into radiation_interactions; merge with branch
!   radiation (MS)
! - added call for usm_green_heat_model for green building surfaces (RvT)
! - added call for usm_temperature_near_surface for use in indoor model (RvT)
!
! 2617 2017-11-16 12:47:24Z suehring
! Bugfix, assure that the reference state does not become zero.
!
! 2563 2017-10-19 15:36:10Z Giersch
! Variable wind_turbine moved to module control_parameters
!
! 2365 2017-08-21 14:59:59Z kanani
! Vertical grid nesting implemented (SadiqHuq)
!
! 2320 2017-07-21 12:47:43Z suehring
! Set bottom boundary conditions after nesting interpolation and anterpolation
!
! 2299 2017-06-29 10:14:38Z maronga
! Call of soil model adjusted
!
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison'
! includes two more prognostic equations for cloud drop concentration (nc)
! and cloud water content (qc).
!
! 2271 2017-06-09 12:34:55Z sward
! Start timestep message changed
!
! 2259 2017-06-08 09:09:11Z gronemeier
! Implemented synthetic turbulence generator
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography and surface concept
! Modify passed parameters for disturb_field
!
! 2178 2017-03-17 11:07:39Z hellstea
! Setting perturbations at all times near inflow boundary is removed
! in case of nested boundaries
!
! 2174 2017-03-13 08:18:57Z maronga
! Added support for nesting with cloud microphysics
!
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC directives and related code removed
!
! 2050 2016-11-08 15:00:55Z gronemeier
! Implement turbulent outflow condition
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean
!
! 2011 2016-09-19 17:29:57Z kanani
! Flag urban_surface is now defined in module control_parameters,
! removed commented CALLs of global_min_max.
!
! 2007 2016-08-24 15:47:17Z kanani
! Added CALLs for new urban surface model
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1976 2016-07-27 13:28:04Z maronga
! Simplified calls to radiation model
!
! 1960 2016-07-12 16:34:24Z suehring
! Separate humidity and passive scalar
!
! 1957 2016-07-07 10:43:48Z suehring
! flight module added
!
! 1919 2016-05-27 14:51:23Z raasch
! Initial version of purely vertical nesting introduced.
!
! 1918 2016-05-27 14:35:57Z raasch
! determination of time step moved to the end of the time step loop,
! the first time step is now always calculated before the time step loop (i.e.
! also in case of restart runs)
!
! 1914 2016-05-26 14:44:07Z witha
! Added call for wind turbine model
!
! 1878 2016-04-19 12:30:36Z hellstea
! Synchronization for nested runs rewritten
!
! 1853 2016-04-11 09:00:35Z maronga
! Adjusted for use with radiation_scheme = constant
!
! 1849 2016-04-08 11:33:18Z hoffmann
! Adapted for modularization of microphysics
!
! 1833 2016-04-07 14:23:03Z raasch
! spectrum renamed spectra_mod, spectra related variables moved to spectra_mod
!
! 1831 2016-04-07 13:15:51Z hoffmann
! turbulence renamed collision_turbulence
!
! 1822 2016-04-07 07:49:42Z hoffmann
! icloud_scheme replaced by microphysics_*
!
! 1808 2016-04-05 19:44:00Z raasch
! output message in case unscheduled radiation calls removed
!
! 1797 2016-03-21 16:50:28Z raasch
! introduction of different datatransfer modes
!
! 1791 2016-03-11 10:41:25Z raasch
! call of pmci_update_new removed
!
! 1786 2016-03-08 05:49:27Z raasch
! +module spectrum
!
! 1783 2016-03-06 18:36:17Z raasch
! switch back of netcdf data format for mask output moved to the mask output
! routine
!
! 1781 2016-03-03 15:12:23Z raasch
! some pmc calls removed at the beginning (before timeloop),
! pmc initialization moved to the main program
!
! 1764 2016-02-28 12:45:19Z raasch
! PMC_ACTIVE flags removed,
! bugfix: nest synchronization after first call of timestep
!
! 1762 2016-02-25 12:31:13Z hellstea
! Introduction of nested domain feature
!
! 1736 2015-12-04 08:56:33Z raasch
! no perturbations added to total domain if energy limit has been set zero
!
! 1691 2015-10-26 16:17:44Z maronga
! Added option for spin-ups without land surface and radiation models. Moved calls
! for radiation and lan surface schemes.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1671 2015-09-25 03:29:37Z raasch
! bugfix: ghostpoint exchange for array diss in case that sgs velocities are used
! for particles
!
! 1585 2015-04-30 07:05:52Z maronga
! Moved call of radiation scheme. Added support for RRTM
!
! 1551 2015-03-03 14:18:16Z maronga
! Added interface for different radiation schemes.
!
! 1496 2014-12-02 17:25:50Z maronga
! Added calls for the land surface model and radiation scheme
!
! 1402 2014-05-09 14:25:13Z raasch
! location messages modified
!
! 1384 2014-05-02 14:31:06Z raasch
! location messages added
!
! 1380 2014-04-28 12:40:45Z heinze
! CALL of nudge_ref added
! bc_pt_t_val and bc_q_t_val are updated in case nudging is used
!
! 1365 2014-04-22 15:03:56Z boeske
! Reset sums_ls_l to zero at each timestep
! +sums_ls_l
! Calculation of reference state (previously in subroutine calc_mean_profile)

! 1342 2014-03-26 17:04:47Z kanani
! REAL constants defined as wp-kind
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! old module precision_kind is removed,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
! 1318 2014-03-17 13:35:16Z raasch
! module interfaces removed
!
! 1308 2014-03-13 14:58:42Z fricke
! +netcdf_data_format_save
! For masked data, parallel netcdf output is not tested so far, hence
! netcdf_data_format is switched back to non-paralell output.
!
! 1276 2014-01-15 13:40:41Z heinze
! Use LSF_DATA also in case of Dirichlet bottom boundary condition for scalars
!
! 1257 2013-11-08 15:18:40Z raasch
! acc-update-host directive for timestep removed
!
! 1241 2013-10-30 11:36:58Z heinze
! Generalize calc_mean_profile for wider use
! Determine shf and qsws in dependence on data from LSF_DATA
! Determine ug and vg in dependence on data from LSF_DATA
! 1221 2013-09-10 08:59:13Z raasch
! host update of arrays before timestep is called
!
! 1179 2013-06-14 05:57:58Z raasch
! mean profiles for reference state are only calculated if required,
! small bugfix for background communication
!
! 1171 2013-05-30 11:27:45Z raasch
! split of prognostic_equations deactivated (comment lines), for the time being
!
! 1128 2013-04-12 06:19:32Z raasch
! asynchronous transfer of ghost point data realized for acc-optimized version:
! prognostic_equations are first called two times for those points required for
! the left-right and north-south exchange, respectively, and then for the
! remaining points,
! those parts requiring global communication moved from prognostic_equations to
! here
!
! 1115 2013-03-26 18:16:16Z hoffmann
! calculation of qr and nr is restricted to precipitation
!
! 1113 2013-03-10 02:48:14Z raasch
! GPU-porting of boundary conditions,
! openACC directives updated
! formal parameter removed from routine boundary_conds
!
! 1111 2013-03-08 23:54:10Z raasch
! +internal timestep counter for cpu statistics added,
! openACC directives updated
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1065 2012-11-22 17:42:36Z hoffmann
! exchange of diss (dissipation rate) in case of turbulence = .TRUE. added
!
! 1053 2012-11-13 17:11:03Z hoffmann
! exchange of ghost points for nr, qr added
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1019 2012-09-28 06:46:45Z raasch
! non-optimized version of prognostic_equations removed
!
! 1015 2012-09-27 09:23:24Z raasch
! +call of prognostic_equations_acc
!
! 1001 2012-09-13 14:08:46Z raasch
! all actions concerning leapfrog- and upstream-spline-scheme removed
!
! 849 2012-03-15 10:35:09Z raasch
! advec_particles renamed lpm, first_call_advec_particles renamed first_call_lpm
!
! 825 2012-02-19 03:03:44Z raasch
! wang_collision_kernel renamed wang_kernel
!
! Revision 1.1  1997/08/11 06:19:04  raasch
! Initial revision
!
!
! Description:
! ------------
!> Integration in time of the model equations, statistical analysis and graphic
!> output
!------------------------------------------------------------------------------!
 SUBROUTINE time_integration


    USE advec_ws,                                                              &
        ONLY:  ws_statistics, ws_finalize

    USE arrays_3d,                                                             &
        ONLY:  p, hyp, dzu, e, e_p, nc, nc_p, nr, nr_p, prho, pt, pt_p, pt_init, sa_init, &
               q_init, q, qc, qc_p, ql, ql_c, ql_v, ql_vp, qr, qr_p, q_p,      &
               ref_state, rho_ocean, s, s_p, sa_p, tend, u, u_p, v,            &
               v_p, w, w_p, alpha_T, beta_S, solar3d, sa, &
               ddzu, ddzw, dzw, dd2zu, drho_air, drho_air_zw,                  &
               rho_air, rho_air_zw, kh, km,                                    &
               te_m, tu_m, tv_m, tw_m, tpt_m, tsa_m,     &
               u_stk, v_stk, ug, vg, u_init, v_init, rdf, rdf_sc,              &
               ptdf_x, ptdf_y

    USE calc_mean_profile_mod,                                                 &
        ONLY:  calc_mean_profile

    USE control_parameters,                                                    &
        ONLY:  advected_distance_x, advected_distance_y,        &
               average_count_3d, averaging_interval, averaging_interval_pr,    &
               bc_lr_cyc, bc_ns_cyc, bc_pt_t_val,                              &
               bc_q_t_val, call_psolver_at_all_substeps,       &
               constant_flux_layer, constant_heatflux,          &
               create_disturbances, dopr_n, coupling_mode, &
               coupling_start_time, current_timestep_number,                   &
               disturbance_created, disturbance_energy_limit, dist_range,      &
               do_sum, dt_3d, dt_averaging_input, dt_averaging_input_pr,       &
               dt_coupling, dt_data_output_av, dt_disturb, dt_do2d_xy,         &
               dt_do2d_xz, dt_do2d_yz, dt_do3d, dt_domask,dt_dopts, dt_dopr,   &
               dt_dopr_listing, dt_dots, dt_run_control, end_time,    &
               forcing, g,      &
               intermediate_timestep_count, intermediate_timestep_count_max,   &
               masks,                   &
               mid,  &
               neutral, nr_timesteps_this_run, nudging,                        &
               passive_scalar, pt_reference,            &
               pt_slope_offset, random_heatflux,                    &
               run_coupled, simulated_time, simulated_time_chr,    &
               skip_time_do2d_xy, skip_time_do2d_xz, skip_time_do2d_yz,        &
               skip_time_do3d, skip_time_domask, skip_time_dopr,               &
               skip_time_data_output_av, sloping_surface,                      &
               stop_dt, terminate_coupled, terminate_run, timestep_scheme,     &
               time_coupling, time_do2d_xy, time_do2d_xz, time_do2d_yz,        &
               time_do3d, time_domask, time_dopr, time_dopr_av,                &
               time_dopr_listing, time_dopts, time_dosp, time_dosp_av,         &
               time_dots, time_do_av, time_do_sla, time_disturb, time_dvrp,    &
               time_run_control, time_since_reference_point, tsc,              &
               turbulent_inflow, turbulent_outflow, urban_surface,             &
               use_initial_profile_as_reference,                               &
               use_single_reference_value, uv_exposure, u_gtrans, v_gtrans,    &
               virtual_flight, wind_turbine, ws_scheme_mom, ws_scheme_sca,     &
               stokes_force, dp_smooth_factor, dpdxy

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt, &
               wall_flags_0, advc_flags_1, advc_flags_2

    USE interfaces

    USE kinds

    USE netcdf_data_input_mod,                                                 &
        ONLY:  force, netcdf_data_input_lsf

    USE pegrid

    USE progress_bar,                                                          &
        ONLY:  finish_progress_bar, output_progress_bar

    USE prognostic_equations_mod,                                              &
        ONLY:  prognostic_equations_vector

    USE statistics,                                                            &
        ONLY:  flow_statistics_called, hom, rmask

    USE surface_layer_fluxes_mod,                                              &
        ONLY:  surface_layer_fluxes

    USE surface_mod,                                                           &
        ONLY:  surf_def_h, bc_h

    USE turbulence_closure_mod,                                                &
        ONLY:  tcm_diffusivities, production_e_init, &
               l_grid, l_wall

    USE stokes_force_mod,                                                      &
        ONLY:  stokes_pressure_head

    IMPLICIT NONE

    CHARACTER (LEN=9) ::  time_to_string          !<
    INTEGER(iwp)      ::  it
    INTEGER(iwp)      ::  lsp
    INTEGER(iwp)      ::  n

    #define MY_DEBUG print *,"DEBUG",__LINE__,__FILE__

    REAL(wp) ::  dt_3d_old  !< temporary storage of timestep to be used for
                            !< steering of run control output interval
    REAL(wp) ::  tsrp_org   !< original value of time_since_reference_point

! !$acc data    copy( tend ) &
! !$acc         copy( tu_m ) &
! !$acc         copy( tv_m ) &
! !$acc         copy( tw_m ) &
! !$acc         copy( tpt_m ) &
! !$acc         copy( tsa_m ) &
! !$acc         copy( u ) &
! !$acc         copy( v ) &
! !$acc         copy( w ) &
! !$acc         copy( pt ) &
! !$acc         copy( sa ) &
! !$acc         copy( u_p ) &
! !$acc         copy( v_p ) &
! !$acc         copy( w_p ) &
! !$acc         copy( pt_p ) &
! !$acc         copy( sa_p ) &
! !$acc       copyin( ddzu ) &
! !$acc       copyin( ddzw ) &
! !$acc       copyin( dd2zu ) &
! !$acc       copyin( u_init ) &
! !$acc       copyin( v_init ) &
! !$acc       copyin( ug ) &
! !$acc       copyin( vg ) &
! !$acc       copyin( pt_init ) &
! !$acc       copyin( sa_init ) &
! !$acc       copyin( tsc ) &
! !$acc       copyin( etdf_x ) &
! !$acc       copyin( wall_flags_0 ) &
! !$acc       copyin( timestep_scheme ) &
! !$acc       copyin( rdf ) &
! !$acc       copyin( rdf_sc )

!$acc data copyin( g ) &
!$acc      copyin( drho_air ) &
!$acc      copyin( drho_air_zw ) &
!$acc      copyin( rho_air ) &
!$acc      copyin( rho_air_zw ) &
!$acc      copyin( dd2zu ) &
!$acc      copyin( ddzu ) &
!$acc      copyin( ddzw ) &
!$acc      copyin( dzw ) &
!$acc      copyin( bc_h ) &
!$acc      copyin( l_grid ) &
!$acc      copyin( l_wall ) &
!$acc      copyin( surf_def_h ) &
!$acc      copyin( rmask ) &
!$acc      copyin( wall_flags_0 ) &
!$acc      copyin( advc_flags_1, advc_flags_2 ) &
!$acc      copyin( dp_smooth_factor, dpdxy ) &
!$acc      copyin( ptdf_x, ptdf_y ) &
!$acc      copyin( hyp ) &
!$acc      copyin( tsc ) &
!$acc       copyin( rdf, rdf_sc ) &
!!$acc      copyin( u_stk, v_stk ) &
!$acc      copyin( u_init, v_init ) &
!$acc      copyin( pt_init, sa_init ) &
!$acc      copyin( u, u_p, tu_m ) &
!$acc      copyin( v, v_p, tv_m ) &
!$acc      copyin( w, w_p, tw_m ) &
!$acc      copyin( e, e_p, te_m ) &
!$acc      copyin( pt, pt_p, tpt_m ) &
!$acc      copyin( sa, sa_p, tsa_m ) &
!$acc      copyin( kh, km ) &
!$acc      copyin( prho, rho_ocean ) &
!$acc      copyin( alpha_T, beta_S, solar3d ) &
!$acc      copyin( ug, vg )

!
!-- At beginning determine the first time step
    CALL timestep
!
!
!-- Determine and print out the run control quantities before the first time
!-- step of this run. For the initial run, some statistics (e.g. divergence)
!-- need to be determined first --> CALL flow_statistics at the beginning of
!-- run_control
    CALL run_control
!
    CALL location_message( 'starting timestep-sequence', .TRUE. )

    !
!-- Start of the time loop
    DO  WHILE ( simulated_time < end_time  .AND.  .NOT. stop_dt  .AND. &
                .NOT. terminate_run )
print *, simulated_time
       CALL cpu_log( log_point_s(10), 'timesteps', 'start' )
!
!--    Start of intermediate step loop
       intermediate_timestep_count = 0
       DO  WHILE ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )

          intermediate_timestep_count = intermediate_timestep_count + 1

!
!--       Set the steering factors for the prognostic equations which depend
!--       on the timestep scheme
          CALL timestep_scheme_steering
          !$acc update device(tsc)

!--          Horizontally averaged profiles to be used as reference state in
!--          buoyancy terms (WARNING: only the respective last call of
!--          calc_mean_profile defines the reference state!)
          CALL calc_mean_profile( rho_ocean, 64 )

                ref_state(:)  = hom(:,1,64,0)
!--          Assure that ref_state does not become zero at any level
!--          ( might be the case if a vertical level is completely occupied
!--            with topography ).
          ref_state = MERGE( MAXVAL(ref_state), ref_state,                  &
                                ref_state == 0.0_wp )

          CALL production_e_init
          IF ( ( ws_scheme_mom .OR. ws_scheme_sca )  .AND.  &
               intermediate_timestep_count == 1 )  CALL ws_statistics
!
          !$acc update device( u, v, w, e, pt, sa, rho_ocean, prho, alpha_T, beta_S)
          CALL prognostic_equations_vector
          !$acc update self( u_p, v_p, w_p, e_p, pt_p, sa_p, tu_m, tv_m, tw_m, tpt_m, tsa_m ) &
          !$acc self( rho_ocean, prho, alpha_T, beta_S, solar3d )
            !
!
!--       Exchange of ghost points (lateral boundary conditions)
          CALL cpu_log( log_point(26), 'exchange-horiz-progn', 'start' )

          CALL exchange_horiz( u_p, nbgp )
          CALL exchange_horiz( v_p, nbgp )
          CALL exchange_horiz( w_p, nbgp )
          CALL exchange_horiz( pt_p, nbgp )
          CALL exchange_horiz( e_p, nbgp )
          CALL exchange_horiz( sa_p, nbgp )
          CALL exchange_horiz( rho_ocean, nbgp )
          CALL exchange_horiz( prho, nbgp )
          CALL exchange_horiz( alpha_T, nbgp )
          CALL exchange_horiz( beta_S, nbgp )
          call exchange_horiz( solar3d, nbgp )
          CALL cpu_log( log_point(26), 'exchange-horiz-progn', 'stop' )

!
!--       Boundary conditions for the prognostic quantities (except of the
!--       velocities at the outflow in case of a non-cyclic lateral wall)
          CALL boundary_conds
!
!--       Swap the time levels in preparation for the next time step.
          CALL swap_timelevel
!
!--       Temperature offset must be imposed at cyclic boundaries in x-direction
!--       when a sloping surface is used
          IF ( sloping_surface )  THEN
             IF ( nxl ==  0 )  pt(:,:,nxlg:nxl-1) = pt(:,:,nxlg:nxl-1) - &
                                                    pt_slope_offset
             IF ( nxr == nx )  pt(:,:,nxr+1:nxrg) = pt(:,:,nxr+1:nxrg) + &
                                                    pt_slope_offset
          ENDIF
!
!--       Impose a random perturbation on the horizontal velocity field
          IF ( create_disturbances  .AND.                                      &
               ( call_psolver_at_all_substeps  .AND.                           &
               intermediate_timestep_count == intermediate_timestep_count_max )&
          .OR. ( .NOT. call_psolver_at_all_substeps  .AND.                     &
               intermediate_timestep_count == 1 ) )                            &
          THEN
             time_disturb = time_disturb + dt_3d
             IF ( time_disturb <= dt_disturb ) then
                  CALL disturb_field( 'u', tend, u )
                  CALL disturb_field( 'v', tend, v )
                  call disturb_field('pt', tend, pt )
          !        call disturb_field('sa', tend, sa )
            ENDIF
          ENDIF
!
!--       Reduce the velocity divergence via the equation for perturbation
!--       pressure.
          CALL pres
!
!--       Compute the diffusion quantities

!--       First the vertical (and horizontal) fluxes in the surface
!--       (constant flux) layer are computed
          IF ( constant_flux_layer )  THEN
             CALL cpu_log( log_point(19), 'surface_layer_fluxes', 'start' )
             CALL surface_layer_fluxes
             CALL cpu_log( log_point(19), 'surface_layer_fluxes', 'stop' )
          ENDIF

!--       Compute the diffusion coefficients
          CALL cpu_log( log_point(17), 'diffusivities', 'start' )
          !$acc update device(e, prho)
          CALL tcm_diffusivities
          !$acc update self(km, kh, e)
          CALL cpu_log( log_point(17), 'diffusivities', 'stop' )
!

       ENDDO   ! Intermediate step loop
!
!--    Update perturbation pressure to account for the Stokes pressure
!--    head, if required
       IF ( stokes_force ) THEN
          CALL stokes_pressure_head
       ENDIF
!
!--    Increase simulation time and output times
       nr_timesteps_this_run      = nr_timesteps_this_run + 1
       current_timestep_number    = current_timestep_number + 1
       simulated_time             = simulated_time   + dt_3d
       time_since_reference_point = simulated_time - coupling_start_time
       simulated_time_chr         = time_to_string( time_since_reference_point )



       IF ( simulated_time >= skip_time_data_output_av )  THEN
          time_do_av         = time_do_av       + dt_3d
       ENDIF
       IF ( simulated_time >= skip_time_do2d_xy )  THEN
          time_do2d_xy       = time_do2d_xy     + dt_3d
       ENDIF
       IF ( simulated_time >= skip_time_do2d_xz )  THEN
          time_do2d_xz       = time_do2d_xz     + dt_3d
       ENDIF
       IF ( simulated_time >= skip_time_do2d_yz )  THEN
          time_do2d_yz       = time_do2d_yz     + dt_3d
       ENDIF
       IF ( simulated_time >= skip_time_do3d    )  THEN
          time_do3d          = time_do3d        + dt_3d
       ENDIF
       DO  mid = 1, masks
          IF ( simulated_time >= skip_time_domask(mid) )  THEN
             time_domask(mid)= time_domask(mid) + dt_3d
          ENDIF
       ENDDO
       time_dots          = time_dots        + dt_3d
       IF ( simulated_time >= skip_time_dopr )  THEN
          time_dopr       = time_dopr        + dt_3d
       ENDIF
       time_dopr_listing          = time_dopr_listing        + dt_3d
       time_run_control   = time_run_control + dt_3d
!
!--    Check, if restart is necessary (because cpu-time is expiring or
!--    because it is forced by user) and set stop flag
!--    This call is skipped if the remote model has already initiated a restart.
!       IF ( .NOT. terminate_run )  CALL check_for_restart
!
!--    Set a flag indicating that so far no statistics have been created
!--    for this time step
       flow_statistics_called = .FALSE.
!
!--    If required, call flow_statistics for averaging in time
       IF ( averaging_interval_pr /= 0.0_wp  .AND.  &
            ( dt_dopr - time_dopr ) <= averaging_interval_pr  .AND.  &
            simulated_time >= skip_time_dopr )  THEN
          time_dopr_av = time_dopr_av + dt_3d
          IF ( time_dopr_av >= dt_averaging_input_pr )  THEN
             do_sum = .TRUE.
             time_dopr_av = MOD( time_dopr_av, &
                                    MAX( dt_averaging_input_pr, dt_3d ) )
          ENDIF
       ENDIF
       IF ( do_sum )  CALL flow_statistics

!
!--    Sum-up 3d-arrays for later output of time-averaged 2d/3d/masked data
       IF ( averaging_interval /= 0.0_wp  .AND.                                &
            ( dt_data_output_av - time_do_av ) <= averaging_interval  .AND. &
            simulated_time >= skip_time_data_output_av )                    &
       THEN
          time_do_sla = time_do_sla + dt_3d
          IF ( time_do_sla >= dt_averaging_input )  THEN
             CALL sum_up_3d_data
             average_count_3d = average_count_3d + 1
             time_do_sla = MOD( time_do_sla, MAX( dt_averaging_input, dt_3d ) )
          ENDIF
       ENDIF
!
!--    Graphic output for PROFIL
       IF ( time_dopr >= dt_dopr )  THEN
          IF ( dopr_n /= 0 )  CALL data_output_profiles
          time_dopr = MOD( time_dopr, MAX( dt_dopr, dt_3d ) )
          time_dopr_av = 0.0_wp    ! due to averaging (see above)
       ENDIF

!
!--    Graphic output for time series
       IF ( time_dots >= dt_dots )  THEN
          CALL data_output_tseries
          time_dots = MOD( time_dots, MAX( dt_dots, dt_3d ) )
       ENDIF
!
!--    3d-data output (volume data)
       IF ( time_do3d >= dt_do3d )  THEN
          CALL data_output_3d( 0 )
          time_do3d = MOD( time_do3d, MAX( dt_do3d, dt_3d ) )
       ENDIF

!
!--    Output of time-averaged 2d/3d/masked data
       IF ( time_do_av >= dt_data_output_av )  THEN
          CALL average_3d_data
         CALL data_output_3d( 1 )
         time_do_av = MOD( time_do_av, MAX( dt_data_output_av, dt_3d ) )
       ENDIF
!--    Determine size of next time step. Save timestep dt_3d because it is
!--    newly calculated in routine timestep, but required further below for
!--    steering the run control output interval
       dt_3d_old = dt_3d
       CALL timestep
!--    Computation and output of run control parameters.
!--    This is also done whenever perturbations have been imposed
       IF ( time_run_control >= dt_run_control  .OR.                     &
            timestep_scheme(1:5) /= 'runge'  .OR.  disturbance_created ) &
       THEN
          CALL run_control
          IF ( time_run_control >= dt_run_control )  THEN
             time_run_control = MOD( time_run_control, &
                                     MAX( dt_run_control, dt_3d_old ) )
          ENDIF
       ENDIF

!
!--    Output elapsed simulated time in form of a progress bar on stdout
!       IF ( myid == 0 )  CALL output_progress_bar

       CALL cpu_log( log_point_s(10), 'timesteps', 'stop' )

    ENDDO   ! time loop

!$acc end data

    call ws_finalize

    IF ( myid == 0 )  CALL finish_progress_bar
    CALL location_message( 'finished time-stepping', .TRUE. )

 END SUBROUTINE time_integration
