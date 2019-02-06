!> @file prognostic_equations.f90
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
! 2019-01-22 cbegeman
! Add timers to prognostic_equations_cache and to tcm_prognostic
!
! 2019-01-22 cbegeman
! Add timers to prognostic_equations_cache and to tcm_prognostic
!
! 2018-10-19 cbegeman
! Change buoyancy call for sloped ocean cases
!
! Former revisions:
! -----------------
! $Id: prognostic_equations.f90 3022 2018-05-18 11:12:35Z suehring $
! Revise recent bugfix for nesting
!
! 3021 2018-05-16 08:14:20Z maronga
! Bugfix in IF clause for nesting
!
! 3014 2018-05-09 08:42:38Z maronga
! Fixed a bug in the IF condition to call pcm_tendency in case of
! potential temperature
!
! 2815 2018-02-19 11:29:57Z kanani
! Rename chem_tendency to chem_prognostic_equations,
! implement vector version for air chemistry
!
! 2766 2018-01-22 17:17:47Z kanani
! Removed preprocessor directive __chem
!
! 2746 2018-01-15 12:06:04Z suehring
! Move flag plant canopy to modules
!
! 2719 2018-01-02 09:02:06Z maronga
! Bugfix for last change.
!
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! - Change in file header (GPL part)
! - Moved TKE equation to tcm_prognostic (TG)
! - Added switch for chemical reactions (RF, FK)
! - Implementation of chemistry module (RF, BK, FK)
!
! 2563 2017-10-19 15:36:10Z Giersch
! Variable wind_turbine moved to module control_parameters
!
! 2320 2017-07-21 12:47:43Z suehring
! Modularize large-scale forcing and nudging
!
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison'
! includes two more prognostic equations for cloud drop concentration (nc)
! and cloud water content (qc).
!
! 2261 2017-06-08 14:25:57Z raasch
! bugfix for r2232: openmp directives removed
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjutst to new surface-type structure. Remove call for usm_wall_heat_flux,
! which is realized directly in diffusion_s now.
!
! 2192 2017-03-22 04:14:10Z raasch
! Bugfix for misplaced and missing openMP directives from r2155
!
! 2155 2017-02-21 09:57:40Z hoffmann
! Bugfix in the calculation of microphysical quantities on ghost points.
!
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC version of subroutine removed
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean
!
! 2011 2016-09-19 17:29:57Z kanani
! Flag urban_surface is now defined in module control_parameters.
!
! 2007 2016-08-24 15:47:17Z kanani
! Added pt tendency calculation based on energy balance at urban surfaces
! (new urban surface model)
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1976 2016-07-27 13:28:04Z maronga
! Simplied calls to radiation model
!
! 1960 2016-07-12 16:34:24Z suehring
! Separate humidity and passive scalar
!
! 1914 2016-05-26 14:44:07Z witha
! Added calls for wind turbine model
!
! 1873 2016-04-18 14:50:06Z maronga
! Module renamed (removed _mod)
!
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
!
! 1826 2016-04-07 12:01:39Z maronga
! Renamed canopy model calls.
!
! 1822 2016-04-07 07:49:42Z hoffmann
! Kessler microphysics scheme moved to microphysics.
!
! 1757 2016-02-22 15:49:32Z maronga
!
! 1691 2015-10-26 16:17:44Z maronga
! Added optional model spin-up without radiation / land surface model calls.
! Formatting corrections.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1585 2015-04-30 07:05:52Z maronga
! Added call for temperature tendency calculation due to radiative flux divergence
!
! 1517 2015-01-07 19:12:25Z hoffmann
! advec_s_bc_mod addded, since advec_s_bc is now a module
!
! 1496 2014-12-02 17:25:50Z maronga
! Renamed "radiation" -> "cloud_top_radiation"
!
! 1484 2014-10-21 10:53:05Z kanani
! Changes due to new module structure of the plant canopy model:
! parameters cthf and plant_canopy moved to module plant_canopy_model_mod.
! Removed double-listing of use_upstream_for_tke in ONLY-list of module
! control_parameters
!
! 1409 2014-05-23 12:11:32Z suehring
! Bugfix: i_omp_start changed for advec_u_ws at left inflow and outflow boundary.
! This ensures that left-hand side fluxes are also calculated for nxl in that
! case, even though the solution at nxl is overwritten in boundary_conds()
!
! 1398 2014-05-07 11:15:00Z heinze
! Rayleigh-damping for horizontal velocity components changed: instead of damping
! against ug and vg, damping against u_init and v_init is used to allow for a
! homogenized treatment in case of nudging
!
! 1380 2014-04-28 12:40:45Z heinze
! Change order of calls for scalar prognostic quantities:
! ls_advec -> nudging -> subsidence since initial profiles
!
! 1374 2014-04-25 12:55:07Z raasch
! missing variables added to ONLY lists
!
! 1365 2014-04-22 15:03:56Z boeske
! Calls of ls_advec for large scale advection added,
! subroutine subsidence is only called if use_subsidence_tendencies = .F.,
! new argument ls_index added to the calls of subsidence
! +ls_index
!
! 1361 2014-04-16 15:17:48Z hoffmann
! Two-moment microphysics moved to the start of prognostic equations. This makes
! the 3d arrays for tend_q, tend_qr, tend_pt and tend_pt redundant.
! Additionally, it is allowed to call the microphysics just once during the time
! step (not at each sub-time step).
!
! Two-moment cloud physics added for vector and accelerator optimization.
!
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1337 2014-03-25 15:11:48Z heinze
! Bugfix: REAL constants provided with KIND-attribute
!
! 1332 2014-03-25 11:59:43Z suehring
! Bugfix: call advec_ws or advec_pw for TKE only if NOT use_upstream_for_tke
!
! 1330 2014-03-24 17:29:32Z suehring
! In case of SGS-particle velocity advection of TKE is also allowed with
! dissipative 5th-order scheme.
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
! 1318 2014-03-17 13:35:16Z raasch
! module interfaces removed
!
! 1257 2013-11-08 15:18:40Z raasch
! openacc loop vector clauses removed, independent clauses added
!
! 1246 2013-11-01 08:59:45Z heinze
! enable nudging also for accelerator version
!
! 1241 2013-10-30 11:36:58Z heinze
! usage of nudging enabled (so far not implemented for accelerator version)
!
! 1179 2013-06-14 05:57:58Z raasch
! two arguments removed from routine buoyancy, ref_state updated on device
!
! 1128 2013-04-12 06:19:32Z raasch
! those parts requiring global communication moved to time_integration,
! loop index bounds in accelerator version replaced by i_left, i_right, j_south,
! j_north
!
! 1115 2013-03-26 18:16:16Z hoffmann
! optimized cloud physics: calculation of microphysical tendencies transfered
! to microphysics.f90; qr and nr are only calculated if precipitation is required
!
! 1111 2013-03-08 23:54:10Z raasch
! update directives for prognostic quantities removed
!
! 1106 2013-03-04 05:31:38Z raasch
! small changes in code formatting
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1053 2012-11-13 17:11:03Z hoffmann
! implementation of two new prognostic equations for rain drop concentration (nr)
! and rain water content (qr)
!
! currently, only available for cache loop optimization
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1019 2012-09-28 06:46:45Z raasch
! non-optimized version of prognostic_equations removed
!
! 1015 2012-09-27 09:23:24Z raasch
! new branch prognostic_equations_acc
! OpenACC statements added + code changes required for GPU optimization
!
! 1001 2012-09-13 14:08:46Z raasch
! all actions concerning leapfrog- and upstream-spline-scheme removed
!
! 978 2012-08-09 08:28:32Z fricke
! km_damp_x and km_damp_y removed in calls of diffusion_u and diffusion_v
! add ptdf_x, ptdf_y for damping the potential temperature at the inflow
! boundary in case of non-cyclic lateral boundaries
! Bugfix: first thread index changes for WS-scheme at the inflow
!
! 940 2012-07-09 14:31:00Z raasch
! temperature equation can be switched off
!
! Revision 1.1  2000/04/13 14:56:27  schroeter
! Initial revision
!
!
! Description:
! ------------
!> Solving the prognostic equations.
!------------------------------------------------------------------------------!
 MODULE prognostic_equations_mod



    USE arrays_3d,                                                             &
        ONLY:  diss_l_e, diss_l_nc, diss_l_nr, diss_l_pt, diss_l_q, diss_l_qc, &
               diss_l_qr, diss_l_s, diss_l_sa, diss_s_e, diss_s_nc, diss_s_nr, &
               diss_s_pt, diss_s_q, diss_s_qc, diss_s_qr, diss_s_s, diss_s_sa, &
               e, e_p, flux_s_e, flux_s_nc, flux_s_nr, flux_s_pt, flux_s_q,    &
               flux_s_qc, flux_s_qr, flux_s_s, flux_s_sa, flux_l_e, flux_l_nc, &
               flux_l_nr, flux_l_pt, flux_l_q, flux_l_qc, flux_l_qr, flux_l_s, &
               flux_l_sa, nc, nc_p, nr, nr_p, pt, ptdf_x, ptdf_y, pt_init,     &
               pt_p, prho, q, q_init, q_p, qc, qc_p, qr, qr_p, rdf, rdf_sc,    &
               ref_state, rho_ocean, s,  s_init, s_p, sa, sa_init, sa_p, tend, &
               te_m, tnc_m,  tnr_m, tpt_m, tq_m, tqc_m, tqr_m, ts_m, tsa_m,    &
               tu_m, tv_m, tw_m, u, ug, u_init, u_p, v, vg, vpt, v_init, v_p,  &
               w, w_p, alpha_T, beta_S, drho_ref_zu

    USE chemistry_model_mod,                                                   &
        ONLY:  chem_integrate, chem_prognostic_equations,                      &
               chem_species, nspec, nvar, spc_names

    USE chem_photolysis_mod,                                                   &
        ONLY:  photolysis_control

    USE chem_modules,                                                          &
        ONLY:  call_chem_at_all_substeps, chem_gasphase_on

    USE kinds

    USE control_parameters,                                                    &
        ONLY:  air_chemistry, call_microphysics_at_all_substeps,               &
               cloud_physics, cloud_top_radiation, constant_diffusion,         &
               dp_external, dp_level_ind_b, dp_smooth_factor, dpdxy, dt_3d,    &
               humidity, idealized_diurnal, g,                                 &
               inflow_l, intermediate_timestep_count,                          &
               intermediate_timestep_count_max, large_scale_forcing,           &
               large_scale_subsidence, message_string, microphysics_morrison,  &
               microphysics_seifert, microphysics_sat_adjust, neutral, nudging,&
               ocean, outflow_l, outflow_s, passive_scalar, plant_canopy,      &
               prho_reference, pt_reference, pt_reference, pt_reference,       &
               scalar_advec, scalar_advec, simulated_time, sloping_surface,    &
               rho_reference, time_prog_terms, timestep_scheme, tsc,           &
               use_single_reference_value, use_subsidence_tendencies,          &
               use_upstream_for_tke, wind_turbine, ws_scheme_mom,              &
               ws_scheme_sca, urban_surface, land_surface, wb_solar,           &
               stokes_force

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s


    USE eqn_state_seawater_mod,                                                &
        ONLY:  eqn_state_seawater

    USE indices,                                                               &
        ONLY:  nbgp, nxl, nxlg, nxlu, nxr, nxrg, nyn, nyng, nys, nysg, nysv,   &
               nzb, nzt, wall_flags_0

    USE advec_ws,                                                              &
        ONLY:  advec_s_ws, advec_u_ws, advec_v_ws, advec_w_ws

    USE advec_s_bc_mod,                                                        &
        ONLY:  advec_s_bc

    USE advec_s_pw_mod,                                                        &
        ONLY:  advec_s_pw

    USE advec_s_up_mod,                                                        &
        ONLY:  advec_s_up

    USE advec_u_pw_mod,                                                        &
        ONLY:  advec_u_pw

    USE advec_u_up_mod,                                                        &
        ONLY:  advec_u_up

    USE advec_v_pw_mod,                                                        &
        ONLY:  advec_v_pw

    USE advec_v_up_mod,                                                        &
        ONLY:  advec_v_up

    USE advec_w_pw_mod,                                                        &
        ONLY:  advec_w_pw

    USE advec_w_up_mod,                                                        &
        ONLY:  advec_w_up

    USE buoyancy_mod,                                                          &
        ONLY:  buoyancy

    USE calc_radiation_mod,                                                    &
        ONLY:  calc_radiation

    USE coriolis_mod,                                                          &
        ONLY:  coriolis

    USE diffusion_s_mod,                                                       &
        ONLY:  diffusion_s

    USE diffusion_u_mod,                                                       &
        ONLY:  diffusion_u

    USE diffusion_v_mod,                                                       &
        ONLY:  diffusion_v

    USE diffusion_w_mod,                                                       &
        ONLY:  diffusion_w

    USE kinds

    USE lsf_nudging_mod,                                                       &
        ONLY:  ls_advec, nudge

    USE microphysics_mod,                                                      &
        ONLY:  microphysics_control

    USE plant_canopy_model_mod,                                                &
        ONLY:  cthf, pcm_tendency

    USE radiation_model_mod,                                                   &
        ONLY:  radiation, radiation_tendency,                                  &
               skip_time_do_radiation

    USE statistics,                                                            &
        ONLY:  hom

    USE subsidence_mod,                                                        &
        ONLY:  subsidence

    USE turbulence_closure_mod,                                                &
        ONLY:  tcm_prognostic

    USE user_actions_mod,                                                      &
        ONLY:  user_actions

    USE surface_mod,                                                           &
        ONLY :  surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v, surf_usm_h,    &
                surf_usm_v

    USE constants

    USE wind_turbine_model_mod,                                                &
        ONLY:  wtm_tendencies

    USE stokes_force_mod,                                                      &
        ONLY: stokes_force_uvw, stokes_force_s

    PRIVATE
    PUBLIC prognostic_equations_cache, prognostic_equations_vector

    INTERFACE prognostic_equations_cache
       MODULE PROCEDURE prognostic_equations_cache
    END INTERFACE prognostic_equations_cache

    INTERFACE prognostic_equations_vector
       MODULE PROCEDURE prognostic_equations_vector
    END INTERFACE prognostic_equations_vector


 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Version with one optimized loop over all equations. It is only allowed to
!> be called for the Wicker and Skamarock or Piascek-Williams advection scheme.
!>
!> Here the calls of most subroutines are embedded in two DO loops over i and j,
!> so communication between CPUs is not allowed (does not make sense) within
!> these loops.
!>
!> (Optimized to avoid cache missings, i.e. for Power4/5-architectures.)
!------------------------------------------------------------------------------!

 SUBROUTINE prognostic_equations_cache


    IMPLICIT NONE

    INTEGER(iwp) ::  i                   !<
    INTEGER(iwp) ::  i_omp_start         !<
    INTEGER(iwp) ::  j                   !<
    INTEGER(iwp) ::  k                   !<
    INTEGER(iwp) ::  omp_get_thread_num  !<
    INTEGER(iwp) ::  tn = 0              !<
    INTEGER(iwp) ::  m

    LOGICAL      ::  loop_start          !<
    INTEGER      ::  n, lsp              !< lsp running index for chem spcs
    REAL(WP)      ::  wb_sfc, tod,arg1      !< surface buoyancy forcing -- only matters for ocean

!
!-- Time measurement can only be performed for the whole set of equations
    CALL cpu_log( log_point(32), 'all progn.equations', 'start' )

!
!-- Calculation of chemical reactions. This is done outside of main loop,
!-- since exchange of ghost points is required after this update of the
!-- concentrations of chemical species
    IF ( air_chemistry )  THEN
!
!--    If required, calculate photolysis frequencies -
!--    UNFINISHED: Why not before the intermediate timestep loop?
       IF ( intermediate_timestep_count ==  1 )  THEN
          CALL photolysis_control
       ENDIF
!
!--    Chemical reactions
       CALL cpu_log( log_point(82), '(chem react + exch_h)', 'start' )

       IF ( chem_gasphase_on ) THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn

                IF ( intermediate_timestep_count == 1 .OR.                        &
                     call_chem_at_all_substeps ) THEN
                   CALL chem_integrate (i,j)
                ENDIF
             ENDDO
          ENDDO
       ENDIF
!
!--    Loop over chemical species
       CALL cpu_log( log_point_s(84), 'chemistry exch-horiz ', 'start' )
       DO  n = 1, nspec
          CALL exchange_horiz( chem_species(n)%conc, nbgp )
       ENDDO
       CALL cpu_log( log_point_s(84), 'chemistry exch-horiz ', 'stop' )

       CALL cpu_log( log_point(82), '(chem react + exch_h)', 'stop' )

    ENDIF

!
!-- If required, calculate cloud microphysics
    IF ( cloud_physics  .AND.  .NOT. microphysics_sat_adjust  .AND.            &
         ( intermediate_timestep_count == 1  .OR.                              &
           call_microphysics_at_all_substeps ) )                               &
    THEN
       CALL cpu_log( log_point(51), 'microphysics', 'start')
       !$OMP PARALLEL PRIVATE (i,j)
       !$OMP DO
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             CALL microphysics_control( i, j )
           ENDDO
       ENDDO
       !$OMP END PARALLEL
       CALL cpu_log( log_point(51), 'microphysics', 'stop' )
    ENDIF

!
!-- Loop over all prognostic equations
    !$OMP PARALLEL PRIVATE (i,i_omp_start,j,k,loop_start,tn)

    !$  tn = omp_get_thread_num()
    loop_start = .TRUE.

    k = nzt
    DO i = nxl, nxr
       DO j = nys,nyn
              IF (idealized_diurnal) THEN
                 m = surf_def_h(2)%start_index(j,i)
                 wb_sfc = g*(surf_def_h(2)%shf(m)*alpha_T(k,j,i) -        &
                                     surf_def_h(2)%sasws(m)*beta_S(k,j,i))
                 tod = simulated_time / 86400.0_wp
                 arg1 = cos(2.0_wp*pi*(tod - 0.5_wp))
                 surf_def_h(2)%shf_sol(m) = wb_solar*max(arg1,0.0_wp)
             ENDIF
       enddo
    enddo
    !$OMP DO


    DO  i = nxl, nxr

!
!--    Store the first loop index. It differs for each thread and is required
!--    later in advec_ws
       IF ( loop_start )  THEN
          loop_start  = .FALSE.
          i_omp_start = i
       ENDIF

       DO  j = nys, nyn
!
!--       Tendency terms for u-velocity component. Please note, in case of
!--       non-cyclic boundary conditions the grid point i=0 is excluded from
!--       the prognostic equations for the u-component.
          CALL cpu_log( log_point(5), 'u-equation', 'start' )
          IF ( i >= nxlu )  THEN

             IF ( time_prog_terms ) CALL cpu_log( log_point(43), 'advec-uvw', 'start' )
             tend(:,j,i) = 0.0_wp

             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_mom )  THEN
                   CALL advec_u_ws( i, j, i_omp_start, tn )
                ELSE
                   CALL advec_u_pw( i, j )
                ENDIF
             ELSE
                CALL advec_u_up( i, j )
             ENDIF
             IF ( time_prog_terms ) CALL cpu_log( log_point(43), 'advec-uvw', 'stop' )
             
             IF ( time_prog_terms ) CALL cpu_log( log_point(44), 'diffusion', 'start' )
             CALL diffusion_u( i, j )
             IF ( time_prog_terms ) CALL cpu_log( log_point(44), 'diffusion', 'stop' )

             IF ( time_prog_terms ) CALL cpu_log( log_point(45), 'coriolis', 'start' )
             CALL coriolis( i, j, 1 )
             IF ( time_prog_terms ) CALL cpu_log( log_point(45), 'coriolis', 'stop' )
             
             IF ( sloping_surface  .AND.  .NOT. neutral )  THEN
                IF ( time_prog_terms ) CALL cpu_log( log_point(46), 'buoyancy', 'start' )
                IF ( ocean ) THEN
                   CALL buoyancy( i, j, rho_ocean, 1 )
                ELSE
                   CALL buoyancy( i, j, pt, 1 )
                ENDIF
                IF ( time_prog_terms ) CALL cpu_log( log_point(46), 'buoyancy', 'stop' )
             ENDIF
!
!--          If required, compute Stokes forces
             IF ( ocean .AND. stokes_force ) THEN
                IF ( time_prog_terms ) CALL cpu_log( log_point(47), 'stokes', 'start' )
                CALL stokes_force_uvw( i, j, 1 )
                IF ( time_prog_terms ) CALL cpu_log( log_point(47), 'stokes', 'stop' )
             ENDIF

!
!--          Drag by plant canopy
             IF ( plant_canopy )  CALL pcm_tendency( i, j, 1 )

!
!--          External pressure gradient
             IF ( dp_external )  THEN
                DO  k = dp_level_ind_b+1, nzt
                   tend(k,j,i) = tend(k,j,i) - dpdxy(1) * dp_smooth_factor(k)   &
                                 * MERGE( 1.0_wp/rho_reference, drho_ref_zu(k), &
                                          use_single_reference_value )
                ENDDO
             ENDIF

!
!--          Nudging
             IF ( nudging )  CALL nudge( i, j, simulated_time, 'u' )

!
!--          Forces by wind turbines
             IF ( wind_turbine )  CALL wtm_tendencies( i, j, 1 )

             CALL user_actions( i, j, 'u-tendency' )

!
!--          Prognostic equation for u-velocity component
             CALL cpu_log( log_point(48), 'prog-u', 'start' )
             DO  k = nzb+1, nzt

                u_p(k,j,i) = u(k,j,i) + ( dt_3d *                               &
                                            ( tsc(2) * tend(k,j,i) +            &
                                              tsc(3) * tu_m(k,j,i) )            &
                                            - tsc(5) * rdf(k)                   &
                                                     * ( u(k,j,i) - u_init(k) ) &
                                        ) * MERGE( 1.0_wp, 0.0_wp,              &
                                                 BTEST( wall_flags_0(k,j,i), 1 )&
                                                 )
             ENDDO

!
!--          Calculate tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   DO  k = nzb+1, nzt
                      tu_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ELSEIF ( intermediate_timestep_count < &
                         intermediate_timestep_count_max )  THEN
                   DO  k = nzb+1, nzt
                      tu_m(k,j,i) =   -9.5625_wp * tend(k,j,i)                &
                                     + 5.3125_wp * tu_m(k,j,i)
                   ENDDO
                ENDIF
             ENDIF

          ENDIF
          CALL cpu_log( log_point(5), 'u-equation', 'stop' )
!
!--       Tendency terms for v-velocity component. Please note, in case of
!--       non-cyclic boundary conditions the grid point j=0 is excluded from
!--       the prognostic equations for the v-component. !--
          CALL cpu_log( log_point(6), 'v-equation', 'start' )
          IF ( j >= nysv )  THEN

             IF ( time_prog_terms ) CALL cpu_log( log_point(43), 'advec-uvw', 'start' )
             tend(:,j,i) = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_mom )  THEN
                    CALL advec_v_ws( i, j, i_omp_start, tn )
                ELSE
                    CALL advec_v_pw( i, j )
                ENDIF
             ELSE
                CALL advec_v_up( i, j )
             ENDIF
             IF ( time_prog_terms ) CALL cpu_log( log_point(43), 'advec-uvw', 'stop' )
             
             IF ( time_prog_terms ) CALL cpu_log( log_point(44), 'diffusion', 'start' )
             CALL diffusion_v( i, j )
             IF ( time_prog_terms ) CALL cpu_log( log_point(44), 'diffusion', 'stop' )
             
             IF ( time_prog_terms ) CALL cpu_log( log_point(45), 'coriolis', 'start' )
             CALL coriolis( i, j, 2 )
             IF ( time_prog_terms ) CALL cpu_log( log_point(45), 'coriolis', 'stop' )

!
!--          If required, compute Stokes forces
             IF ( ocean .AND. stokes_force ) THEN
                IF ( time_prog_terms ) CALL cpu_log( log_point(47), 'stokes', 'start' )
                CALL stokes_force_uvw( i, j, 2 )
                IF ( time_prog_terms ) CALL cpu_log( log_point(47), 'stokes', 'stop' )
             ENDIF

!
!--          Drag by plant canopy
             IF ( plant_canopy )  CALL pcm_tendency( i, j, 2 )

!
!--          External pressure gradient
             IF ( dp_external )  THEN
                DO  k = dp_level_ind_b+1, nzt
                   tend(k,j,i) = tend(k,j,i) - dpdxy(2) * dp_smooth_factor(k)   &
                                 * MERGE( 1.0_wp/rho_reference, drho_ref_zu(k), &
                                          use_single_reference_value )
                ENDDO
             ENDIF

!
!--          Nudging
             IF ( nudging )  CALL nudge( i, j, simulated_time, 'v' )

!
!--          Forces by wind turbines
             IF ( wind_turbine )  CALL wtm_tendencies( i, j, 2 )

             CALL user_actions( i, j, 'v-tendency' )

!
!--          Prognostic equation for v-velocity component
             DO  k = nzb+1, nzt
                v_p(k,j,i) = v(k,j,i) + ( dt_3d *                              &
                                            ( tsc(2) * tend(k,j,i) +           &
                                              tsc(3) * tv_m(k,j,i) )           &
                                            - tsc(5) * rdf(k)                  &
                                                     * ( v(k,j,i) - v_init(k) )&
                                        ) * MERGE( 1.0_wp, 0.0_wp,             &
                                                   BTEST( wall_flags_0(k,j,i), 2 )&
                                                 )
             ENDDO

!
!--          Calculate tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   DO  k = nzb+1, nzt
                      tv_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ELSEIF ( intermediate_timestep_count < &
                         intermediate_timestep_count_max )  THEN
                   DO  k = nzb+1, nzt
                      tv_m(k,j,i) =   -9.5625_wp * tend(k,j,i)                 &
                                     + 5.3125_wp * tv_m(k,j,i)
                   ENDDO
                ENDIF
             ENDIF

          ENDIF
          CALL cpu_log( log_point(6), 'v-equation', 'stop' )

!
!--       Tendency terms for w-velocity component
          CALL cpu_log( log_point(7), 'w-equation', 'start' )
          
          IF ( time_prog_terms ) CALL cpu_log( log_point(43), 'advec-uvw', 'start' )
          tend(:,j,i) = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_mom )  THEN
                CALL advec_w_ws( i, j, i_omp_start, tn )
             ELSE
                CALL advec_w_pw( i, j )
             END IF
          ELSE
             CALL advec_w_up( i, j )
          ENDIF
          IF ( time_prog_terms ) CALL cpu_log( log_point(43), 'advec-uvw', 'stop' )
          
          IF ( time_prog_terms ) CALL cpu_log( log_point(44), 'diffusion', 'start' )
          CALL diffusion_w( i, j )
          IF ( time_prog_terms ) CALL cpu_log( log_point(44), 'diffusion', 'stop' )
          
          IF ( time_prog_terms ) CALL cpu_log( log_point(45), 'coriolis', 'start' )
          CALL coriolis( i, j, 3 )
          IF ( time_prog_terms ) CALL cpu_log( log_point(45), 'coriolis', 'stop' )

          IF ( .NOT. neutral )  THEN
             IF ( time_prog_terms ) CALL cpu_log( log_point(46), 'buoyancy', 'start' )
             IF ( ocean )  THEN
                CALL buoyancy( i, j, rho_ocean, 3 )
             ELSE
                IF ( .NOT. humidity )  THEN
                   CALL buoyancy( i, j, pt, 3 )
                ELSE
                   CALL buoyancy( i, j, vpt, 3 )
                ENDIF
             ENDIF
             IF ( time_prog_terms ) CALL cpu_log( log_point(46), 'buoyancy', 'stop' )
          ENDIF

!
!--       If required, compute Stokes forces
          IF ( ocean .AND. stokes_force ) THEN
             IF ( time_prog_terms ) CALL cpu_log( log_point(47), 'stokes', 'start' )
             CALL stokes_force_uvw( i, j, 3 )
             IF ( time_prog_terms ) CALL cpu_log( log_point(47), 'stokes', 'stop' )
          ENDIF

!
!--       Drag by plant canopy
          IF ( plant_canopy )  CALL pcm_tendency( i, j, 3 )

!
!--       Forces by wind turbines
          IF ( wind_turbine )  CALL wtm_tendencies( i, j, 3 )

          CALL user_actions( i, j, 'w-tendency' )

!
!--       Prognostic equation for w-velocity component
          DO  k = nzb+1, nzt-1
             w_p(k,j,i) = w(k,j,i) + ( dt_3d *                                 &
                                             ( tsc(2) * tend(k,j,i) +          &
                                               tsc(3) * tw_m(k,j,i) )          &
                                             - tsc(5) * rdf(k) * w(k,j,i)      &
                                     ) * MERGE( 1.0_wp, 0.0_wp,                &
                                                BTEST( wall_flags_0(k,j,i), 3 )&
                                              )
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  k = nzb+1, nzt-1
                   tw_m(k,j,i) = tend(k,j,i)
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  k = nzb+1, nzt-1
                   tw_m(k,j,i) =   -9.5625_wp * tend(k,j,i)                    &
                                  + 5.3125_wp * tw_m(k,j,i)
                ENDDO
             ENDIF
          ENDIF
          CALL cpu_log( log_point(7), 'w-equation', 'stop' )

!
!--       If required, compute prognostic equation for potential temperature
          IF ( .NOT. neutral )  THEN
             CALL cpu_log( log_point(13), 'pt-equation', 'start' )
!
!--          Tendency terms for potential temperature
             tend(:,j,i) = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                   IF ( ws_scheme_sca )  THEN
                      CALL advec_s_ws( i, j, pt, 'pt', flux_s_pt, diss_s_pt,   &
                                       flux_l_pt, diss_l_pt, i_omp_start, tn )
                   ELSE
                      CALL advec_s_pw( i, j, pt )
                   ENDIF
             ELSE
                CALL advec_s_up( i, j, pt )
             ENDIF
             CALL diffusion_s( i, j, pt,                                       &
                               surf_def_h(0)%shf, surf_def_h(1)%shf,           &
                               surf_def_h(2)%shf,                              &
                               surf_lsm_h%shf,    surf_usm_h%shf,              &
                               surf_def_v(0)%shf, surf_def_v(1)%shf,           &
                               surf_def_v(2)%shf, surf_def_v(3)%shf,           &
                               surf_lsm_v(0)%shf, surf_lsm_v(1)%shf,           &
                               surf_lsm_v(2)%shf, surf_lsm_v(3)%shf,           &
                               surf_usm_v(0)%shf, surf_usm_v(1)%shf,           &
                               surf_usm_v(2)%shf, surf_usm_v(3)%shf,           &
                               surf_def_h(2)%shf_sol )

!
!--          If required, compute Stokes-advection term
             IF ( ocean .AND. stokes_force ) THEN
                CALL stokes_force_s( i, j, pt )
             ENDIF

!
!--          If required compute heating/cooling due to long wave radiation
!--          processes
             IF ( cloud_top_radiation )  THEN
                CALL calc_radiation( i, j )
             ENDIF

!
!--          Consideration of heat sources within the plant canopy
             IF ( plant_canopy  .AND.                                          &
                (cthf /= 0.0_wp  .OR. urban_surface  .OR.  land_surface) )  THEN
                CALL pcm_tendency( i, j, 4 )
             ENDIF

!
!--          Large scale advection
             IF ( large_scale_forcing )  THEN
                CALL ls_advec( i, j, simulated_time, 'pt' )
             ENDIF

!
!--          Nudging
             IF ( nudging )  CALL nudge( i, j, simulated_time, 'pt' )

!
!--          If required, compute effect of large-scale subsidence/ascent
             IF ( large_scale_subsidence  .AND.                                &
                  .NOT. use_subsidence_tendencies )  THEN
                CALL subsidence( i, j, tend, pt, pt_init, 2 )
             ENDIF

!
!--          If required, add tendency due to radiative heating/cooling
             IF ( radiation  .AND.                                             &
                  simulated_time > skip_time_do_radiation )  THEN
                CALL radiation_tendency ( i, j, tend )
             ENDIF


             CALL user_actions( i, j, 'pt-tendency' )

!
!--          Prognostic equation for potential temperature
             DO  k = nzb+1, nzt
                   pt_p(k,j,i) = pt(k,j,i) + ( dt_3d *                            &
                                                  ( tsc(2) * tend(k,j,i) +     &
                                                    tsc(3) * tpt_m(k,j,i) )    &
                                                  - tsc(5) * ( pt(k,j,i) - pt_init(k) ) &
                                                  * ( rdf_sc(k) + ptdf_x(i)    &
                                                                + ptdf_y(j) )  &
                                          )                                    &
                                       * MERGE( 1.0_wp, 0.0_wp,                &
                                                BTEST( wall_flags_0(k,j,i), 0 )&
                                              )
             ENDDO

!
!--          Calculate tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   DO  k = nzb+1, nzt
                      tpt_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ELSEIF ( intermediate_timestep_count < &
                         intermediate_timestep_count_max )  THEN
                   DO  k = nzb+1, nzt
                      tpt_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +              &
                                        5.3125_wp * tpt_m(k,j,i)
                   ENDDO
                ENDIF
             ENDIF

             CALL cpu_log( log_point(13), 'pt-equation', 'stop' )

          ENDIF

!
!--       If required, compute prognostic equation for salinity
          IF ( ocean )  THEN

             CALL cpu_log( log_point(37), 'sa-equation', 'start' )
!
!--          Tendency-terms for salinity
             tend(:,j,i) = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' ) &
             THEN
                IF ( ws_scheme_sca )  THEN
                    CALL advec_s_ws( i, j, sa, 'sa', flux_s_sa,  &
                                diss_s_sa, flux_l_sa, diss_l_sa, i_omp_start, tn  )
                ELSE
                    CALL advec_s_pw( i, j, sa )
                ENDIF
             ELSE
                CALL advec_s_up( i, j, sa )
             ENDIF

             CALL diffusion_s( i, j, sa,                                       &
                               surf_def_h(0)%sasws, surf_def_h(1)%sasws,       &
                               surf_def_h(2)%sasws,                            &
                               surf_lsm_h%sasws,    surf_usm_h%sasws,          &
                               surf_def_v(0)%sasws, surf_def_v(1)%sasws,       &
                               surf_def_v(2)%sasws, surf_def_v(3)%sasws,       &
                               surf_lsm_v(0)%sasws, surf_lsm_v(1)%sasws,       &
                               surf_lsm_v(2)%sasws, surf_lsm_v(3)%sasws,       &
                               surf_usm_v(0)%sasws, surf_usm_v(1)%sasws,       &
                               surf_usm_v(2)%sasws, surf_usm_v(3)%sasws )

!
!--          If required, compute Stokes-advection term
             IF ( stokes_force ) THEN
                CALL stokes_force_s( i, j, sa )
             ENDIF


             CALL user_actions( i, j, 'sa-tendency' )

!
!--          Prognostic equation for salinity
             DO  k = nzb+1, nzt
                sa_p(k,j,i) = sa(k,j,i) + ( dt_3d *                            &
                                                  ( tsc(2) * tend(k,j,i) +     &
                                                    tsc(3) * tsa_m(k,j,i) )    &
                                                  - tsc(5) * rdf_sc(k)         &
                                                   * ( sa(k,j,i) - sa_init(k) )&
                                          ) * MERGE( 1.0_wp, 0.0_wp,           &
                                                BTEST( wall_flags_0(k,j,i), 0 )&
                                                   )
                IF ( sa_p(k,j,i) < 0.0_wp )  sa_p(k,j,i) = 0.1_wp * sa(k,j,i)
             ENDDO

!
!--          Calculate tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   DO  k = nzb+1, nzt
                      tsa_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ELSEIF ( intermediate_timestep_count < &
                         intermediate_timestep_count_max )  THEN
                   DO  k = nzb+1, nzt
                      tsa_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +              &
                                        5.3125_wp * tsa_m(k,j,i)
                   ENDDO
                ENDIF
             ENDIF

             CALL cpu_log( log_point(37), 'sa-equation', 'stop' )
!
!--          Calculate density by the equation of state for seawater
             CALL cpu_log( log_point(38), 'eqns-seawater', 'start' )
             CALL eqn_state_seawater( i, j )
             CALL cpu_log( log_point(38), 'eqns-seawater', 'stop' )

          ENDIF

!
!--       If required, compute prognostic equation for total water content
          IF ( humidity )  THEN

             CALL cpu_log( log_point(29), 'q-equation', 'start' )
!
!--          Tendency-terms for total water content / scalar
             tend(:,j,i) = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' ) &
             THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws( i, j, q, 'q', flux_s_q, &
                                diss_s_q, flux_l_q, diss_l_q, i_omp_start, tn )
                ELSE
                   CALL advec_s_pw( i, j, q )
                ENDIF
             ELSE
                CALL advec_s_up( i, j, q )
             ENDIF
             CALL diffusion_s( i, j, q,                                        &
                               surf_def_h(0)%qsws, surf_def_h(1)%qsws,         &
                               surf_def_h(2)%qsws,                             &
                               surf_lsm_h%qsws,    surf_usm_h%qsws,            &
                               surf_def_v(0)%qsws, surf_def_v(1)%qsws,         &
                               surf_def_v(2)%qsws, surf_def_v(3)%qsws,         &
                               surf_lsm_v(0)%qsws, surf_lsm_v(1)%qsws,         &
                               surf_lsm_v(2)%qsws, surf_lsm_v(3)%qsws,         &
                               surf_usm_v(0)%qsws, surf_usm_v(1)%qsws,         &
                               surf_usm_v(2)%qsws, surf_usm_v(3)%qsws )

!
!--          Sink or source of humidity due to canopy elements
             IF ( plant_canopy )  CALL pcm_tendency( i, j, 5 )

!
!--          Large scale advection
             IF ( large_scale_forcing )  THEN
                CALL ls_advec( i, j, simulated_time, 'q' )
             ENDIF

!
!--          Nudging
             IF ( nudging )  CALL nudge( i, j, simulated_time, 'q' )

!
!--          If required compute influence of large-scale subsidence/ascent
             IF ( large_scale_subsidence  .AND.                                &
                  .NOT. use_subsidence_tendencies )  THEN
                CALL subsidence( i, j, tend, q, q_init, 3 )
             ENDIF

             CALL user_actions( i, j, 'q-tendency' )

!
!--          Prognostic equation for total water content / scalar
             DO  k = nzb+1, nzt
                q_p(k,j,i) = q(k,j,i) + ( dt_3d *                              &
                                                ( tsc(2) * tend(k,j,i) +       &
                                                  tsc(3) * tq_m(k,j,i) )       &
                                                - tsc(5) * rdf_sc(k) *         &
                                                      ( q(k,j,i) - q_init(k) ) &
                                        )                                      &
                                       * MERGE( 1.0_wp, 0.0_wp,                &
                                                BTEST( wall_flags_0(k,j,i), 0 )&
                                              )
                IF ( q_p(k,j,i) < 0.0_wp )  q_p(k,j,i) = 0.1_wp * q(k,j,i)
             ENDDO

!
!--          Calculate tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   DO  k = nzb+1, nzt
                      tq_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ELSEIF ( intermediate_timestep_count < &
                         intermediate_timestep_count_max )  THEN
                   DO  k = nzb+1, nzt
                      tq_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +               &
                                       5.3125_wp * tq_m(k,j,i)
                   ENDDO
                ENDIF
             ENDIF

             CALL cpu_log( log_point(29), 'q-equation', 'stop' )
!
!--          If required, calculate prognostic equations for cloud water content
!--          and cloud drop concentration
             IF ( cloud_physics  .AND.  microphysics_morrison )  THEN
!
                CALL cpu_log( log_point(67), 'qc-equation', 'start' )

!--             Calculate prognostic equation for cloud water content
                tend(:,j,i) = 0.0_wp
                IF ( timestep_scheme(1:5) == 'runge' ) &
                THEN
                   IF ( ws_scheme_sca )  THEN
                      CALL advec_s_ws( i, j, qc, 'qc', flux_s_qc,       &
                                       diss_s_qc, flux_l_qc, diss_l_qc, &
                                       i_omp_start, tn )
                   ELSE
                      CALL advec_s_pw( i, j, qc )
                   ENDIF
                ELSE
                   CALL advec_s_up( i, j, qc )
                ENDIF
                CALL diffusion_s( i, j, qc,                                   &
                                  surf_def_h(0)%qcsws, surf_def_h(1)%qcsws,   &
                                  surf_def_h(2)%qcsws,                        &
                                  surf_lsm_h%qcsws,    surf_usm_h%qcsws,      &
                                  surf_def_v(0)%qcsws, surf_def_v(1)%qcsws,   &
                                  surf_def_v(2)%qcsws, surf_def_v(3)%qcsws,   &
                                  surf_lsm_v(0)%qcsws, surf_lsm_v(1)%qcsws,   &
                                  surf_lsm_v(2)%qcsws, surf_lsm_v(3)%qcsws,   &
                                  surf_usm_v(0)%qcsws, surf_usm_v(1)%qcsws,   &
                                  surf_usm_v(2)%qcsws, surf_usm_v(3)%qcsws )

!
!--             Prognostic equation for cloud water content
                DO  k = nzb+1, nzt
                   qc_p(k,j,i) = qc(k,j,i) + ( dt_3d *                         &
                                                      ( tsc(2) * tend(k,j,i) + &
                                                        tsc(3) * tqc_m(k,j,i) )&
                                                      - tsc(5) * rdf_sc(k)     &
                                                               * qc(k,j,i)     &
                                             )                                 &
                                       * MERGE( 1.0_wp, 0.0_wp,                &
                                                BTEST( wall_flags_0(k,j,i), 0 )&
                                              )
                   IF ( qc_p(k,j,i) < 0.0_wp )  qc_p(k,j,i) = 0.0_wp
                ENDDO
!
!--             Calculate tendencies for the next Runge-Kutta step
                IF ( timestep_scheme(1:5) == 'runge' )  THEN
                   IF ( intermediate_timestep_count == 1 )  THEN
                      DO  k = nzb+1, nzt
                         tqc_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ELSEIF ( intermediate_timestep_count < &
                            intermediate_timestep_count_max )  THEN
                      DO  k = nzb+1, nzt
                         tqc_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +           &
                                           5.3125_wp * tqc_m(k,j,i)
                      ENDDO
                   ENDIF
                ENDIF

                CALL cpu_log( log_point(67), 'qc-equation', 'stop' )
                CALL cpu_log( log_point(68), 'nc-equation', 'start' )
!
!--             Calculate prognostic equation for cloud drop concentration.
                tend(:,j,i) = 0.0_wp
                IF ( timestep_scheme(1:5) == 'runge' )  THEN
                   IF ( ws_scheme_sca )  THEN
                      CALL advec_s_ws( i, j, nc, 'nc', flux_s_nc,    &
                                    diss_s_nc, flux_l_nc, diss_l_nc, &
                                    i_omp_start, tn )
                   ELSE
                      CALL advec_s_pw( i, j, nc )
                   ENDIF
                ELSE
                   CALL advec_s_up( i, j, nc )
                ENDIF
                CALL diffusion_s( i, j, nc,                                    &
                                  surf_def_h(0)%ncsws, surf_def_h(1)%ncsws,    &
                                  surf_def_h(2)%ncsws,                         &
                                  surf_lsm_h%ncsws,    surf_usm_h%ncsws,       &
                                  surf_def_v(0)%ncsws, surf_def_v(1)%ncsws,    &
                                  surf_def_v(2)%ncsws, surf_def_v(3)%ncsws,    &
                                  surf_lsm_v(0)%ncsws, surf_lsm_v(1)%ncsws,    &
                                  surf_lsm_v(2)%ncsws, surf_lsm_v(3)%ncsws,    &
                                  surf_usm_v(0)%ncsws, surf_usm_v(1)%ncsws,    &
                                  surf_usm_v(2)%ncsws, surf_usm_v(3)%ncsws )

!
!--             Prognostic equation for cloud drop concentration
                DO  k = nzb+1, nzt
                   nc_p(k,j,i) = nc(k,j,i) + ( dt_3d *                         &
                                                      ( tsc(2) * tend(k,j,i) + &
                                                        tsc(3) * tnc_m(k,j,i) )&
                                                      - tsc(5) * rdf_sc(k)     &
                                                               * nc(k,j,i)     &
                                             )                                 &
                                       * MERGE( 1.0_wp, 0.0_wp,                &
                                                BTEST( wall_flags_0(k,j,i), 0 )&
                                              )
                   IF ( nc_p(k,j,i) < 0.0_wp )  nc_p(k,j,i) = 0.0_wp
                ENDDO
!
!--             Calculate tendencies for the next Runge-Kutta step
                IF ( timestep_scheme(1:5) == 'runge' )  THEN
                   IF ( intermediate_timestep_count == 1 )  THEN
                      DO  k = nzb+1, nzt
                         tnc_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ELSEIF ( intermediate_timestep_count < &
                            intermediate_timestep_count_max )  THEN
                      DO  k = nzb+1, nzt
                         tnc_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +           &
                                           5.3125_wp * tnc_m(k,j,i)
                      ENDDO
                   ENDIF
                ENDIF

                CALL cpu_log( log_point(68), 'nc-equation', 'stop' )

             ENDIF
!
!--          If required, calculate prognostic equations for rain water content
!--          and rain drop concentration
             IF ( cloud_physics  .AND.  microphysics_seifert )  THEN

                CALL cpu_log( log_point(52), 'qr-equation', 'start' )
!
!--             Calculate prognostic equation for rain water content
                tend(:,j,i) = 0.0_wp
                IF ( timestep_scheme(1:5) == 'runge' ) &
                THEN
                   IF ( ws_scheme_sca )  THEN
                      CALL advec_s_ws( i, j, qr, 'qr', flux_s_qr,       &
                                       diss_s_qr, flux_l_qr, diss_l_qr, &
                                       i_omp_start, tn )
                   ELSE
                      CALL advec_s_pw( i, j, qr )
                   ENDIF
                ELSE
                   CALL advec_s_up( i, j, qr )
                ENDIF
                CALL diffusion_s( i, j, qr,                                   &
                                  surf_def_h(0)%qrsws, surf_def_h(1)%qrsws,   &
                                  surf_def_h(2)%qrsws,                        &
                                  surf_lsm_h%qrsws,    surf_usm_h%qrsws,      &
                                  surf_def_v(0)%qrsws, surf_def_v(1)%qrsws,   &
                                  surf_def_v(2)%qrsws, surf_def_v(3)%qrsws,   &
                                  surf_lsm_v(0)%qrsws, surf_lsm_v(1)%qrsws,   &
                                  surf_lsm_v(2)%qrsws, surf_lsm_v(3)%qrsws,   &
                                  surf_usm_v(0)%qrsws, surf_usm_v(1)%qrsws,   &
                                  surf_usm_v(2)%qrsws, surf_usm_v(3)%qrsws )

!
!--             Prognostic equation for rain water content
                DO  k = nzb+1, nzt
                   qr_p(k,j,i) = qr(k,j,i) + ( dt_3d *                         &
                                                      ( tsc(2) * tend(k,j,i) + &
                                                        tsc(3) * tqr_m(k,j,i) )&
                                                      - tsc(5) * rdf_sc(k)     &
                                                               * qr(k,j,i)     &
                                             )                                 &
                                       * MERGE( 1.0_wp, 0.0_wp,                &
                                                BTEST( wall_flags_0(k,j,i), 0 )&
                                              )
                   IF ( qr_p(k,j,i) < 0.0_wp )  qr_p(k,j,i) = 0.0_wp
                ENDDO
!
!--             Calculate tendencies for the next Runge-Kutta step
                IF ( timestep_scheme(1:5) == 'runge' )  THEN
                   IF ( intermediate_timestep_count == 1 )  THEN
                      DO  k = nzb+1, nzt
                         tqr_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ELSEIF ( intermediate_timestep_count < &
                            intermediate_timestep_count_max )  THEN
                      DO  k = nzb+1, nzt
                         tqr_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +           &
                                           5.3125_wp * tqr_m(k,j,i)
                      ENDDO
                   ENDIF
                ENDIF

                CALL cpu_log( log_point(52), 'qr-equation', 'stop' )
                CALL cpu_log( log_point(53), 'nr-equation', 'start' )
!
!--             Calculate prognostic equation for rain drop concentration.
                tend(:,j,i) = 0.0_wp
                IF ( timestep_scheme(1:5) == 'runge' )  THEN
                   IF ( ws_scheme_sca )  THEN
                      CALL advec_s_ws( i, j, nr, 'nr', flux_s_nr,    &
                                    diss_s_nr, flux_l_nr, diss_l_nr, &
                                    i_omp_start, tn )
                   ELSE
                      CALL advec_s_pw( i, j, nr )
                   ENDIF
                ELSE
                   CALL advec_s_up( i, j, nr )
                ENDIF
                CALL diffusion_s( i, j, nr,                                    &
                                  surf_def_h(0)%nrsws, surf_def_h(1)%nrsws,    &
                                  surf_def_h(2)%nrsws,                         &
                                  surf_lsm_h%nrsws,    surf_usm_h%nrsws,       &
                                  surf_def_v(0)%nrsws, surf_def_v(1)%nrsws,    &
                                  surf_def_v(2)%nrsws, surf_def_v(3)%nrsws,    &
                                  surf_lsm_v(0)%nrsws, surf_lsm_v(1)%nrsws,    &
                                  surf_lsm_v(2)%nrsws, surf_lsm_v(3)%nrsws,    &
                                  surf_usm_v(0)%nrsws, surf_usm_v(1)%nrsws,    &
                                  surf_usm_v(2)%nrsws, surf_usm_v(3)%nrsws )

!
!--             Prognostic equation for rain drop concentration
                DO  k = nzb+1, nzt
                   nr_p(k,j,i) = nr(k,j,i) + ( dt_3d *                         &
                                                      ( tsc(2) * tend(k,j,i) + &
                                                        tsc(3) * tnr_m(k,j,i) )&
                                                      - tsc(5) * rdf_sc(k)     &
                                                               * nr(k,j,i)     &
                                             )                                 &
                                       * MERGE( 1.0_wp, 0.0_wp,                &
                                                BTEST( wall_flags_0(k,j,i), 0 )&
                                              )
                   IF ( nr_p(k,j,i) < 0.0_wp )  nr_p(k,j,i) = 0.0_wp
                ENDDO
!
!--             Calculate tendencies for the next Runge-Kutta step
                IF ( timestep_scheme(1:5) == 'runge' )  THEN
                   IF ( intermediate_timestep_count == 1 )  THEN
                      DO  k = nzb+1, nzt
                         tnr_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ELSEIF ( intermediate_timestep_count < &
                            intermediate_timestep_count_max )  THEN
                      DO  k = nzb+1, nzt
                         tnr_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +           &
                                           5.3125_wp * tnr_m(k,j,i)
                      ENDDO
                   ENDIF
                ENDIF

                CALL cpu_log( log_point(53), 'nr-equation', 'stop' )

             ENDIF

          ENDIF

!
!--       If required, compute prognostic equation for scalar
          IF ( passive_scalar )  THEN
       
             CALL cpu_log( log_point(66), 's-equation', 'start' )
!
!--          Tendency-terms for total water content / scalar
             tend(:,j,i) = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' ) &
             THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws( i, j, s, 's', flux_s_s, &
                                diss_s_s, flux_l_s, diss_l_s, i_omp_start, tn )
                ELSE
                   CALL advec_s_pw( i, j, s )
                ENDIF
             ELSE
                CALL advec_s_up( i, j, s )
             ENDIF
             CALL diffusion_s( i, j, s,                                        &
                               surf_def_h(0)%ssws, surf_def_h(1)%ssws,         &
                               surf_def_h(2)%ssws,                             &
                               surf_lsm_h%ssws,    surf_usm_h%ssws,            &
                               surf_def_v(0)%ssws, surf_def_v(1)%ssws,         &
                               surf_def_v(2)%ssws, surf_def_v(3)%ssws,         &
                               surf_lsm_v(0)%ssws, surf_lsm_v(1)%ssws,         &
                               surf_lsm_v(2)%ssws, surf_lsm_v(3)%ssws,         &
                               surf_usm_v(0)%ssws, surf_usm_v(1)%ssws,         &
                               surf_usm_v(2)%ssws, surf_usm_v(3)%ssws )

!
!--          If required, compute Stokes-advection term
             IF ( ocean .AND. stokes_force ) THEN
                CALL stokes_force_s( i, j, s )
             ENDIF

!
!--          Sink or source of scalar concentration due to canopy elements
             IF ( plant_canopy )  CALL pcm_tendency( i, j, 7 )

!
!--          Large scale advection, still need to be extended for scalars
!              IF ( large_scale_forcing )  THEN
!                 CALL ls_advec( i, j, simulated_time, 's' )
!              ENDIF

!
!--          Nudging, still need to be extended for scalars
!              IF ( nudging )  CALL nudge( i, j, simulated_time, 's' )

!
!--          If required compute influence of large-scale subsidence/ascent.
!--          Note, the last argument is of no meaning in this case, as it is
!--          only used in conjunction with large_scale_forcing, which is to
!--          date not implemented for scalars.
             IF ( large_scale_subsidence  .AND.                                &
                  .NOT. use_subsidence_tendencies  .AND.                       &
                  .NOT. large_scale_forcing )  THEN
                CALL subsidence( i, j, tend, s, s_init, 3 )
             ENDIF

             CALL user_actions( i, j, 's-tendency' )

!
!--          Prognostic equation for scalar
             DO  k = nzb+1, nzt
                s_p(k,j,i) = s(k,j,i) + (  dt_3d *                             &
                                            ( tsc(2) * tend(k,j,i) +           &
                                              tsc(3) * ts_m(k,j,i) )           &
                                            - tsc(5) * rdf_sc(k)               &
                                                     * ( s(k,j,i) - s_init(k) )&
                                        )                                      &
                                       * MERGE( 1.0_wp, 0.0_wp,                &
                                                BTEST( wall_flags_0(k,j,i), 0 )&
                                              )
                IF ( s_p(k,j,i) < 0.0_wp )  s_p(k,j,i) = 0.1_wp * s(k,j,i)
             ENDDO

!
!--          Calculate tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   DO  k = nzb+1, nzt
                      ts_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ELSEIF ( intermediate_timestep_count < &
                         intermediate_timestep_count_max )  THEN
                   DO  k = nzb+1, nzt
                      ts_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +               &
                                       5.3125_wp * ts_m(k,j,i)
                   ENDDO
                ENDIF
             ENDIF

             CALL cpu_log( log_point(66), 's-equation', 'stop' )

          ENDIF

!
!--       Calculate prognostic equations for turbulence closure
          CALL cpu_log( log_point(41), 'tcm-equation', 'start' )
          CALL tcm_prognostic( i, j, i_omp_start, tn )
          CALL cpu_log( log_point(41), 'tcm-equation', 'stop' )

!
!--       If required, compute prognostic equation for chemical quantites
          IF ( air_chemistry )  THEN
             CALL cpu_log( log_point(83), '(chem advec+diff+prog)', 'start' )
!
!--          Loop over chemical species
             DO  lsp = 1, nvar
                CALL chem_prognostic_equations ( chem_species(lsp)%conc_p,     &
                                     chem_species(lsp)%conc,                   &
                                     chem_species(lsp)%tconc_m,                &
                                     chem_species(lsp)%conc_pr_init,           &
                                     i, j, i_omp_start, tn, lsp,               &
                                     chem_species(lsp)%flux_s_cs,              &
                                     chem_species(lsp)%diss_s_cs,              &
                                     chem_species(lsp)%flux_l_cs,              &
                                     chem_species(lsp)%diss_l_cs )
             ENDDO

             CALL cpu_log( log_point(83), '(chem advec+diff+prog)', 'stop' )
          ENDIF   ! Chemicals equations

       ENDDO
    ENDDO
    !$OMP END PARALLEL


    CALL cpu_log( log_point(32), 'all progn.equations', 'stop' )


 END SUBROUTINE prognostic_equations_cache


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Version for vector machines
!------------------------------------------------------------------------------!

 SUBROUTINE prognostic_equations_vector


    IMPLICIT NONE

    INTEGER(iwp) ::  i     !<
    INTEGER(iwp) ::  j     !<
    INTEGER(iwp) ::  k     !<
    INTEGER(iwp) ::  lsp   !< running index for chemical species
    INTEGER(iwp) ::  m

    REAL(wp)     ::  sbt  !<
    REAL(WP)      ::  wb_sfc, tod,arg1      !< surface buoyancy forcing -- only matters for ocean


!
!-- If required, calculate cloud microphysical impacts
    IF ( cloud_physics  .AND.  .NOT. microphysics_sat_adjust  .AND.            &
         ( intermediate_timestep_count == 1  .OR.                              &
           call_microphysics_at_all_substeps )                                 &
       )  THEN
       CALL cpu_log( log_point(51), 'microphysics', 'start' )
       CALL microphysics_control
       CALL cpu_log( log_point(51), 'microphysics', 'stop' )
    ENDIF

!
!-- u-velocity component
    CALL cpu_log( log_point(5), 'u-equation', 'start' )

    tend = 0.0_wp
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( ws_scheme_mom )  THEN
          CALL advec_u_ws
       ELSE
          CALL advec_u_pw
       ENDIF
    ELSE
       CALL advec_u_up
    ENDIF
    CALL diffusion_u
    CALL coriolis( 1 )
    IF ( sloping_surface  .AND.  .NOT. neutral )  THEN
       IF ( ocean ) THEN
          CALL buoyancy( rho_ocean, 1 )
       ELSE
          CALL buoyancy( pt, 1 )
       ENDIF
    ENDIF

!
!-- If required, compute Stokes forces
    IF ( ocean .AND. stokes_force ) THEN
       CALL stokes_force_uvw( 1 )
    ENDIF

!
!-- Drag by plant canopy
    IF ( plant_canopy )  CALL pcm_tendency( 1 )

!
!-- External pressure gradient
    IF ( dp_external )  THEN
       DO  i = nxlu, nxr
          DO  j = nys, nyn
             DO  k = dp_level_ind_b+1, nzt
                tend(k,j,i) = tend(k,j,i) - dpdxy(1) * dp_smooth_factor(k)   &
                              * MERGE( 1.0_wp/rho_reference, drho_ref_zu(k), &
                                       use_single_reference_value )
             ENDDO
          ENDDO
       ENDDO
    ENDIF

!
!-- Nudging
    IF ( nudging )  CALL nudge( simulated_time, 'u' )

!
!-- Forces by wind turbines
    IF ( wind_turbine )  CALL wtm_tendencies( 1 )

    CALL user_actions( 'u-tendency' )

!
!-- Prognostic equation for u-velocity component
    DO  i = nxlu, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             u_p(k,j,i) = u(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) +          &
                                                 tsc(3) * tu_m(k,j,i) )          &
                                               - tsc(5) * rdf(k) *               &
                                                        ( u(k,j,i) - u_init(k) ) &
                                     ) * MERGE( 1.0_wp, 0.0_wp,                  &
                                                 BTEST( wall_flags_0(k,j,i), 1 ) &
                                              )
          ENDDO
       ENDDO
    ENDDO

!
!-- Calculate tendencies for the next Runge-Kutta step
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( intermediate_timestep_count == 1 )  THEN
          DO  i = nxlu, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   tu_m(k,j,i) = tend(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ELSEIF ( intermediate_timestep_count < &
                intermediate_timestep_count_max )  THEN
          DO  i = nxlu, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   tu_m(k,j,i) =    -9.5625_wp * tend(k,j,i)                   &
                                   + 5.3125_wp * tu_m(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    CALL cpu_log( log_point(5), 'u-equation', 'stop' )

!
!-- v-velocity component
    CALL cpu_log( log_point(6), 'v-equation', 'start' )

    tend = 0.0_wp
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( ws_scheme_mom )  THEN
          CALL advec_v_ws
       ELSE
          CALL advec_v_pw
       END IF
    ELSE
       CALL advec_v_up
    ENDIF
    CALL diffusion_v
    CALL coriolis( 2 )

!
!-- If required, compute Stokes forces
    IF ( ocean .AND. stokes_force ) THEN
       CALL stokes_force_uvw( 2 )
    ENDIF

!
!-- Drag by plant canopy
    IF ( plant_canopy )  CALL pcm_tendency( 2 )

!
!-- External pressure gradient
    IF ( dp_external )  THEN
       DO  i = nxl, nxr
          DO  j = nysv, nyn
             DO  k = dp_level_ind_b+1, nzt
                tend(k,j,i) = tend(k,j,i) - dpdxy(2) * dp_smooth_factor(k)   &
                              * MERGE( 1.0_wp/rho_reference, drho_ref_zu(k), &
                                       use_single_reference_value )
             ENDDO
          ENDDO
       ENDDO
    ENDIF

!
!-- Nudging
    IF ( nudging )  CALL nudge( simulated_time, 'v' )

!
!-- Forces by wind turbines
    IF ( wind_turbine )  CALL wtm_tendencies( 2 )

    CALL user_actions( 'v-tendency' )

!
!-- Prognostic equation for v-velocity component
    DO  i = nxl, nxr
       DO  j = nysv, nyn
          DO  k = nzb+1, nzt
             v_p(k,j,i) = v(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) +        &
                                                 tsc(3) * tv_m(k,j,i) )        &
                                               - tsc(5) * rdf(k) *             &
                                                      ( v(k,j,i) - v_init(k) ) &
                                     ) * MERGE( 1.0_wp, 0.0_wp,                &
                                                BTEST( wall_flags_0(k,j,i), 2 )&
                                              )
          ENDDO
       ENDDO
    ENDDO

!
!-- Calculate tendencies for the next Runge-Kutta step
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( intermediate_timestep_count == 1 )  THEN
          DO  i = nxl, nxr
             DO  j = nysv, nyn
                DO  k = nzb+1, nzt
                   tv_m(k,j,i) = tend(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ELSEIF ( intermediate_timestep_count < &
                intermediate_timestep_count_max )  THEN
          DO  i = nxl, nxr
             DO  j = nysv, nyn
                DO  k = nzb+1, nzt
                   tv_m(k,j,i) =   -9.5625_wp * tend(k,j,i)                    &
                                  + 5.3125_wp * tv_m(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    CALL cpu_log( log_point(6), 'v-equation', 'stop' )

!
!-- w-velocity component
    CALL cpu_log( log_point(7), 'w-equation', 'start' )

    tend = 0.0_wp
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( ws_scheme_mom )  THEN
          CALL advec_w_ws
       ELSE
          CALL advec_w_pw
       ENDIF
    ELSE
       CALL advec_w_up
    ENDIF
    CALL diffusion_w
    CALL coriolis( 3 )

    IF ( .NOT. neutral )  THEN
       IF ( ocean )  THEN
          CALL buoyancy( rho_ocean, 3 )
       ELSE
          IF ( .NOT. humidity )  THEN
             CALL buoyancy( pt, 3 )
          ELSE
             CALL buoyancy( vpt, 3 )
          ENDIF
       ENDIF
    ENDIF

!
!-- If required, compute Stokes forces
    IF ( ocean .AND. stokes_force ) THEN
       CALL stokes_force_uvw( 3 )
    ENDIF

!
!-- Drag by plant canopy
    IF ( plant_canopy )  CALL pcm_tendency( 3 )

!
!-- Forces by wind turbines
    IF ( wind_turbine )  CALL wtm_tendencies( 3 )

    CALL user_actions( 'w-tendency' )

!
!-- Prognostic equation for w-velocity component
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt-1
             w_p(k,j,i) = w(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) +        &
                                                 tsc(3) * tw_m(k,j,i) )        &
                                               - tsc(5) * rdf(k) * w(k,j,i)    &
                                     ) * MERGE( 1.0_wp, 0.0_wp,                &
                                                BTEST( wall_flags_0(k,j,i), 3 )&
                                              )
          ENDDO
       ENDDO
    ENDDO

!
!-- Calculate tendencies for the next Runge-Kutta step
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( intermediate_timestep_count == 1 )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt-1
                   tw_m(k,j,i) = tend(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ELSEIF ( intermediate_timestep_count < &
                intermediate_timestep_count_max )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt-1
                   tw_m(k,j,i) =   -9.5625_wp * tend(k,j,i)                    &
                                  + 5.3125_wp * tw_m(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    CALL cpu_log( log_point(7), 'w-equation', 'stop' )


!
!-- If required, compute prognostic equation for potential temperature
    IF ( .NOT. neutral )  THEN

       CALL cpu_log( log_point(13), 'pt-equation', 'start' )

!
!--    pt-tendency terms with communication
       sbt = tsc(2)
       IF ( scalar_advec == 'bc-scheme' )  THEN

          IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--          Bott-Chlond scheme always uses Euler time step. Thus:
             sbt = 1.0_wp
          ENDIF
          tend = 0.0_wp
          CALL advec_s_bc( pt, 'pt' )

       ENDIF

!
!--    pt-tendency terms with no communication
       IF ( scalar_advec /= 'bc-scheme' )  THEN
          tend = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( pt, 'pt' )
             ELSE
                CALL advec_s_pw( pt )
             ENDIF
          ELSE
             CALL advec_s_up( pt )
          ENDIF
       ENDIF

       k = nzt
       DO i = nxl, nxr
          DO j = nys,nyn
                 IF (idealized_diurnal) THEN
                    m = surf_def_h(2)%start_index(j,i)
                    wb_sfc = g*(surf_def_h(2)%shf(m)*alpha_T(k,j,i) -        &
                                     surf_def_h(2)%sasws(m)*beta_S(k,j,i))
                    tod = simulated_time / 86400.0_wp
                    arg1 = cos(2.0_wp*pi*(tod - 0.75_wp))
                    surf_def_h(2)%shf_sol(m) = wb_solar*max(arg1,0.0_wp)
                ENDIF
          enddo
       enddo

       CALL diffusion_s( pt,                                                   &
                         surf_def_h(0)%shf, surf_def_h(1)%shf,                 &
                         surf_def_h(2)%shf,                                    &
                         surf_lsm_h%shf,    surf_usm_h%shf,                    &
                         surf_def_v(0)%shf, surf_def_v(1)%shf,                 &
                         surf_def_v(2)%shf, surf_def_v(3)%shf,                 &
                         surf_lsm_v(0)%shf, surf_lsm_v(1)%shf,                 &
                         surf_lsm_v(2)%shf, surf_lsm_v(3)%shf,                 &
                         surf_usm_v(0)%shf, surf_usm_v(1)%shf,                 &
                         surf_usm_v(2)%shf, surf_usm_v(3)%shf,                 &
                         surf_def_h(2)%shf_sol )

!
!--    If required, compute Stokes-advection term
       IF ( ocean .AND. stokes_force ) THEN
          CALL stokes_force_s( pt )
       ENDIF

!
!--    If required compute heating/cooling due to long wave radiation processes
       IF ( cloud_top_radiation )  THEN
          CALL calc_radiation
       ENDIF

!
!--    Consideration of heat sources within the plant canopy
       IF ( plant_canopy  .AND.                                          &
            (cthf /= 0.0_wp  .OR. urban_surface  .OR.  land_surface) )  THEN
          CALL pcm_tendency( 4 )
       ENDIF

!
!--    Large scale advection
       IF ( large_scale_forcing )  THEN
          CALL ls_advec( simulated_time, 'pt' )
       ENDIF

!
!--    Nudging
       IF ( nudging )  CALL nudge( simulated_time, 'pt' )

!
!--    If required compute influence of large-scale subsidence/ascent
       IF ( large_scale_subsidence  .AND.                                      &
            .NOT. use_subsidence_tendencies )  THEN
          CALL subsidence( tend, pt, pt_init, 2 )
       ENDIF

!
!--    If required, add tendency due to radiative heating/cooling
       IF ( radiation  .AND.                                                   &
            simulated_time > skip_time_do_radiation )  THEN
            CALL radiation_tendency ( tend )
       ENDIF

       CALL user_actions( 'pt-tendency' )

!
!--    Prognostic equation for potential temperature
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                pt_p(k,j,i) = pt(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +      &
                                                   tsc(3) * tpt_m(k,j,i) )     &
                                                 - tsc(5) *                    &
                                                   ( pt(k,j,i) - pt_init(k) ) *&
                                          ( rdf_sc(k) + ptdf_x(i) + ptdf_y(j) )&
                                          )                                    &
                                   * MERGE( 1.0_wp, 0.0_wp,                    &
                                             BTEST( wall_flags_0(k,j,i), 0 )   &
                                          )
             ENDDO
          ENDDO
       ENDDO
!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      tpt_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      tpt_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +              &
                                        5.3125_wp * tpt_m(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       CALL cpu_log( log_point(13), 'pt-equation', 'stop' )

    ENDIF

!
!-- If required, compute prognostic equation for salinity
    IF ( ocean )  THEN

       CALL cpu_log( log_point(37), 'sa-equation', 'start' )

!
!--    sa-tendency terms with communication
       sbt = tsc(2)
       IF ( scalar_advec == 'bc-scheme' )  THEN

          IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--          Bott-Chlond scheme always uses Euler time step. Thus:
             sbt = 1.0_wp
          ENDIF
          tend = 0.0_wp
          CALL advec_s_bc( sa, 'sa' )

       ENDIF

!
!--    sa-tendency terms with no communication
       IF ( scalar_advec /= 'bc-scheme' )  THEN
          tend = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                 CALL advec_s_ws( sa, 'sa' )
             ELSE
                 CALL advec_s_pw( sa )
             ENDIF
          ELSE
             CALL advec_s_up( sa )
          ENDIF
       ENDIF

       CALL diffusion_s( sa,                                                   &
                         surf_def_h(0)%sasws, surf_def_h(1)%sasws,             &
                         surf_def_h(2)%sasws,                                  &
                         surf_lsm_h%sasws,    surf_usm_h%sasws,                &
                         surf_def_v(0)%sasws, surf_def_v(1)%sasws,             &
                         surf_def_v(2)%sasws, surf_def_v(3)%sasws,             &
                         surf_lsm_v(0)%sasws, surf_lsm_v(1)%sasws,             &
                         surf_lsm_v(2)%sasws, surf_lsm_v(3)%sasws,             &
                         surf_usm_v(0)%sasws, surf_usm_v(1)%sasws,             &
                         surf_usm_v(2)%sasws, surf_usm_v(3)%sasws )

!
!--    If required, compute Stokes-advection term
       IF ( stokes_force ) THEN
          CALL stokes_force_s( sa )
       ENDIF

       CALL user_actions( 'sa-tendency' )

!
!--    Prognostic equation for salinity
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                sa_p(k,j,i) = sa(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +      &
                                                   tsc(3) * tsa_m(k,j,i) )     &
                                                 - tsc(5) * rdf_sc(k) *        &
                                                 ( sa(k,j,i) - sa_init(k) )    &
                                          )                                    &
                                   * MERGE( 1.0_wp, 0.0_wp,                    &
                                             BTEST( wall_flags_0(k,j,i), 0 )   &
                                          )
                IF ( sa_p(k,j,i) < 0.0_wp )  sa_p(k,j,i) = 0.1_wp * sa(k,j,i)
             ENDDO
          ENDDO
       ENDDO

!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      tsa_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      tsa_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +              &
                                        5.3125_wp * tsa_m(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       CALL cpu_log( log_point(37), 'sa-equation', 'stop' )

!
!--    Calculate density by the equation of state for seawater
       CALL cpu_log( log_point(38), 'eqns-seawater', 'start' )
       CALL eqn_state_seawater
       CALL cpu_log( log_point(38), 'eqns-seawater', 'stop' )

    ENDIF

!
!-- If required, compute prognostic equation for total water content
    IF ( humidity )  THEN

       CALL cpu_log( log_point(29), 'q-equation', 'start' )

!
!--    Scalar/q-tendency terms with communication
       sbt = tsc(2)
       IF ( scalar_advec == 'bc-scheme' )  THEN

          IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--          Bott-Chlond scheme always uses Euler time step. Thus:
             sbt = 1.0_wp
          ENDIF
          tend = 0.0_wp
          CALL advec_s_bc( q, 'q' )

       ENDIF

!
!--    Scalar/q-tendency terms with no communication
       IF ( scalar_advec /= 'bc-scheme' )  THEN
          tend = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( q, 'q' )
             ELSE
                CALL advec_s_pw( q )
             ENDIF
          ELSE
             CALL advec_s_up( q )
          ENDIF
       ENDIF

       CALL diffusion_s( q,                                                    &
                         surf_def_h(0)%qsws, surf_def_h(1)%qsws,               &
                         surf_def_h(2)%qsws,                                   &
                         surf_lsm_h%qsws,    surf_usm_h%qsws,                  &
                         surf_def_v(0)%qsws, surf_def_v(1)%qsws,               &
                         surf_def_v(2)%qsws, surf_def_v(3)%qsws,               &
                         surf_lsm_v(0)%qsws, surf_lsm_v(1)%qsws,               &
                         surf_lsm_v(2)%qsws, surf_lsm_v(3)%qsws,               &
                         surf_usm_v(0)%qsws, surf_usm_v(1)%qsws,               &
                         surf_usm_v(2)%qsws, surf_usm_v(3)%qsws )

!
!--    Sink or source of humidity due to canopy elements
       IF ( plant_canopy ) CALL pcm_tendency( 5 )

!
!--    Large scale advection
       IF ( large_scale_forcing )  THEN
          CALL ls_advec( simulated_time, 'q' )
       ENDIF

!
!--    Nudging
       IF ( nudging )  CALL nudge( simulated_time, 'q' )

!
!--    If required compute influence of large-scale subsidence/ascent
       IF ( large_scale_subsidence  .AND.                                      &
            .NOT. use_subsidence_tendencies )  THEN
         CALL subsidence( tend, q, q_init, 3 )
       ENDIF

       CALL user_actions( 'q-tendency' )

!
!--    Prognostic equation for total water content
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                q_p(k,j,i) = q(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +        &
                                                 tsc(3) * tq_m(k,j,i) )        &
                                               - tsc(5) * rdf_sc(k) *          &
                                                      ( q(k,j,i) - q_init(k) ) &
                                        ) * MERGE( 1.0_wp, 0.0_wp,             &
                                                   BTEST( wall_flags_0(k,j,i), 0 ) &
                                                 )
                IF ( q_p(k,j,i) < 0.0_wp )  q_p(k,j,i) = 0.1_wp * q(k,j,i)
             ENDDO
          ENDDO
       ENDDO

!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      tq_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      tq_m(k,j,i) =   -9.5625_wp * tend(k,j,i)                 &
                                     + 5.3125_wp * tq_m(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       CALL cpu_log( log_point(29), 'q-equation', 'stop' )

!
!--    If required, calculate prognostic equations for cloud water content
!--    and cloud drop concentration
       IF ( cloud_physics  .AND.  microphysics_morrison )  THEN

          CALL cpu_log( log_point(67), 'qc-equation', 'start' )

!
!--       Calculate prognostic equation for cloud water content
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( qc, 'qc' )

          ENDIF

!
!--       qc-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws( qc, 'qc' )
                ELSE
                   CALL advec_s_pw( qc )
                ENDIF
             ELSE
                CALL advec_s_up( qc )
             ENDIF
          ENDIF

          CALL diffusion_s( qc,                                                &
                            surf_def_h(0)%qcsws, surf_def_h(1)%qcsws,          &
                            surf_def_h(2)%qcsws,                               &
                            surf_lsm_h%qcsws,    surf_usm_h%qcsws,             &
                            surf_def_v(0)%qcsws, surf_def_v(1)%qcsws,          &
                            surf_def_v(2)%qcsws, surf_def_v(3)%qcsws,          &
                            surf_lsm_v(0)%qcsws, surf_lsm_v(1)%qcsws,          &
                            surf_lsm_v(2)%qcsws, surf_lsm_v(3)%qcsws,          &
                            surf_usm_v(0)%qcsws, surf_usm_v(1)%qcsws,          &
                            surf_usm_v(2)%qcsws, surf_usm_v(3)%qcsws )

!
!--       Prognostic equation for cloud water content
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   qc_p(k,j,i) = qc(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +   &
                                                      tsc(3) * tqc_m(k,j,i) )  &
                                                    - tsc(5) * rdf_sc(k) *     &
                                                               qc(k,j,i)       &
                                             )                                 &
                                    * MERGE( 1.0_wp, 0.0_wp,                   &
                                             BTEST( wall_flags_0(k,j,i), 0 )   &
                                          )
                   IF ( qc_p(k,j,i) < 0.0_wp )  qc_p(k,j,i) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tqc_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tqc_m(k,j,i) =   -9.5625_wp * tend(k,j,i)             &
                                         + 5.3125_wp * tqc_m(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF

          CALL cpu_log( log_point(67), 'qc-equation', 'stop' )
          CALL cpu_log( log_point(68), 'nc-equation', 'start' )

!
!--       Calculate prognostic equation for cloud drop concentration
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( nc, 'nc' )

          ENDIF

!
!--       nc-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws( nc, 'nc' )
                ELSE
                   CALL advec_s_pw( nc )
                ENDIF
             ELSE
                CALL advec_s_up( nc )
             ENDIF
          ENDIF

          CALL diffusion_s( nc,                                                &
                            surf_def_h(0)%ncsws, surf_def_h(1)%ncsws,          &
                            surf_def_h(2)%ncsws,                               &
                            surf_lsm_h%ncsws,    surf_usm_h%ncsws,             &
                            surf_def_v(0)%ncsws, surf_def_v(1)%ncsws,          &
                            surf_def_v(2)%ncsws, surf_def_v(3)%ncsws,          &
                            surf_lsm_v(0)%ncsws, surf_lsm_v(1)%ncsws,          &
                            surf_lsm_v(2)%ncsws, surf_lsm_v(3)%ncsws,          &
                            surf_usm_v(0)%ncsws, surf_usm_v(1)%ncsws,          &
                            surf_usm_v(2)%ncsws, surf_usm_v(3)%ncsws )

!
!--       Prognostic equation for cloud drop concentration
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   nc_p(k,j,i) = nc(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +   &
                                                      tsc(3) * tnc_m(k,j,i) )  &
                                                    - tsc(5) * rdf_sc(k) *     &
                                                               nc(k,j,i)       &
                                             )                                 &
                                   * MERGE( 1.0_wp, 0.0_wp,                    &
                                             BTEST( wall_flags_0(k,j,i), 0 )   &
                                          )
                   IF ( nc_p(k,j,i) < 0.0_wp )  nc_p(k,j,i) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tnc_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tnc_m(k,j,i) =  -9.5625_wp * tend(k,j,i)             &
                                         + 5.3125_wp * tnc_m(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF

          CALL cpu_log( log_point(68), 'nc-equation', 'stop' )

       ENDIF
!
!--    If required, calculate prognostic equations for rain water content
!--    and rain drop concentration
       IF ( cloud_physics  .AND.  microphysics_seifert )  THEN

          CALL cpu_log( log_point(52), 'qr-equation', 'start' )

!
!--       Calculate prognostic equation for rain water content
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( qr, 'qr' )

          ENDIF

!
!--       qr-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws( qr, 'qr' )
                ELSE
                   CALL advec_s_pw( qr )
                ENDIF
             ELSE
                CALL advec_s_up( qr )
             ENDIF
          ENDIF

          CALL diffusion_s( qr,                                                &
                            surf_def_h(0)%qrsws, surf_def_h(1)%qrsws,          &
                            surf_def_h(2)%qrsws,                               &
                            surf_lsm_h%qrsws,    surf_usm_h%qrsws,             &
                            surf_def_v(0)%qrsws, surf_def_v(1)%qrsws,          &
                            surf_def_v(2)%qrsws, surf_def_v(3)%qrsws,          &
                            surf_lsm_v(0)%qrsws, surf_lsm_v(1)%qrsws,          &
                            surf_lsm_v(2)%qrsws, surf_lsm_v(3)%qrsws,          &
                            surf_usm_v(0)%qrsws, surf_usm_v(1)%qrsws,          &
                            surf_usm_v(2)%qrsws, surf_usm_v(3)%qrsws )

!
!--       Prognostic equation for rain water content
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   qr_p(k,j,i) = qr(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +   &
                                                      tsc(3) * tqr_m(k,j,i) )  &
                                                    - tsc(5) * rdf_sc(k) *     &
                                                               qr(k,j,i)       &
                                             )                                 &
                                    * MERGE( 1.0_wp, 0.0_wp,                   &
                                             BTEST( wall_flags_0(k,j,i), 0 )   &
                                          )
                   IF ( qr_p(k,j,i) < 0.0_wp )  qr_p(k,j,i) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tqr_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tqr_m(k,j,i) =   -9.5625_wp * tend(k,j,i)             &
                                         + 5.3125_wp * tqr_m(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF

          CALL cpu_log( log_point(52), 'qr-equation', 'stop' )
          CALL cpu_log( log_point(53), 'nr-equation', 'start' )

!
!--       Calculate prognostic equation for rain drop concentration
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( nr, 'nr' )

          ENDIF

!
!--       nr-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws( nr, 'nr' )
                ELSE
                   CALL advec_s_pw( nr )
                ENDIF
             ELSE
                CALL advec_s_up( nr )
             ENDIF
          ENDIF

          CALL diffusion_s( nr,                                                &
                            surf_def_h(0)%nrsws, surf_def_h(1)%nrsws,          &
                            surf_def_h(2)%nrsws,                               &
                            surf_lsm_h%nrsws,    surf_usm_h%nrsws,             &
                            surf_def_v(0)%nrsws, surf_def_v(1)%nrsws,          &
                            surf_def_v(2)%nrsws, surf_def_v(3)%nrsws,          &
                            surf_lsm_v(0)%nrsws, surf_lsm_v(1)%nrsws,          &
                            surf_lsm_v(2)%nrsws, surf_lsm_v(3)%nrsws,          &
                            surf_usm_v(0)%nrsws, surf_usm_v(1)%nrsws,          &
                            surf_usm_v(2)%nrsws, surf_usm_v(3)%nrsws )

!
!--       Prognostic equation for rain drop concentration
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   nr_p(k,j,i) = nr(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +   &
                                                      tsc(3) * tnr_m(k,j,i) )  &
                                                    - tsc(5) * rdf_sc(k) *     &
                                                               nr(k,j,i)       &
                                             )                                 &
                                   * MERGE( 1.0_wp, 0.0_wp,                    &
                                             BTEST( wall_flags_0(k,j,i), 0 )   &
                                          )
                   IF ( nr_p(k,j,i) < 0.0_wp )  nr_p(k,j,i) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tnr_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tnr_m(k,j,i) =  -9.5625_wp * tend(k,j,i)             &
                                         + 5.3125_wp * tnr_m(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF

          CALL cpu_log( log_point(53), 'nr-equation', 'stop' )

       ENDIF

    ENDIF
!
!-- If required, compute prognostic equation for scalar
    IF ( passive_scalar )  THEN

       CALL cpu_log( log_point(66), 's-equation', 'start' )

!
!--    Scalar/q-tendency terms with communication
       sbt = tsc(2)
       IF ( scalar_advec == 'bc-scheme' )  THEN

          IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--          Bott-Chlond scheme always uses Euler time step. Thus:
             sbt = 1.0_wp
          ENDIF
          tend = 0.0_wp
          CALL advec_s_bc( s, 's' )

       ENDIF

!
!--    Scalar/q-tendency terms with no communication
       IF ( scalar_advec /= 'bc-scheme' )  THEN
          tend = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( s, 's' )
             ELSE
                CALL advec_s_pw( s )
             ENDIF
          ELSE
             CALL advec_s_up( s )
          ENDIF
       ENDIF

       CALL diffusion_s( s,                                                    &
                         surf_def_h(0)%ssws, surf_def_h(1)%ssws,               &
                         surf_def_h(2)%ssws,                                   &
                         surf_lsm_h%ssws,    surf_usm_h%ssws,                  &
                         surf_def_v(0)%ssws, surf_def_v(1)%ssws,               &
                         surf_def_v(2)%ssws, surf_def_v(3)%ssws,               &
                         surf_lsm_v(0)%ssws, surf_lsm_v(1)%ssws,               &
                         surf_lsm_v(2)%ssws, surf_lsm_v(3)%ssws,               &
                         surf_usm_v(0)%ssws, surf_usm_v(1)%ssws,               &
                         surf_usm_v(2)%ssws, surf_usm_v(3)%ssws )

!
!--    If required, compute Stokes-advection term
       IF ( ocean .AND. stokes_force ) THEN
          CALL stokes_force_s( s )
       ENDIF

!
!--    Sink or source of humidity due to canopy elements
       IF ( plant_canopy ) CALL pcm_tendency( 7 )

!
!--    Large scale advection. Not implemented for scalars so far.
!        IF ( large_scale_forcing )  THEN
!           CALL ls_advec( simulated_time, 'q' )
!        ENDIF

!
!--    Nudging. Not implemented for scalars so far.
!        IF ( nudging )  CALL nudge( simulated_time, 'q' )

!
!--    If required compute influence of large-scale subsidence/ascent.
!--    Not implemented for scalars so far.
       IF ( large_scale_subsidence  .AND.                                      &
            .NOT. use_subsidence_tendencies  .AND.                             &
            .NOT. large_scale_forcing )  THEN
         CALL subsidence( tend, s, s_init, 3 )
       ENDIF

       CALL user_actions( 's-tendency' )

!
!--    Prognostic equation for total water content
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                s_p(k,j,i) = s(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +        &
                                                 tsc(3) * ts_m(k,j,i) )        &
                                               - tsc(5) * rdf_sc(k) *          &
                                                    ( s(k,j,i) - s_init(k) )   &
                                        )                                      &
                                   * MERGE( 1.0_wp, 0.0_wp,                    &
                                             BTEST( wall_flags_0(k,j,i), 0 )   &
                                          )
                IF ( s_p(k,j,i) < 0.0_wp )  s_p(k,j,i) = 0.1_wp * s(k,j,i)
             ENDDO
          ENDDO
       ENDDO

!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      ts_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      ts_m(k,j,i) =   -9.5625_wp * tend(k,j,i)                 &
                                     + 5.3125_wp * ts_m(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       CALL cpu_log( log_point(66), 's-equation', 'stop' )

    ENDIF

    CALL cpu_log( log_point(41), 'tcm-equation', 'start' )
    CALL tcm_prognostic()
    CALL cpu_log( log_point(41), 'tcm-equation', 'stop' )

!
!-- If required, compute prognostic equation for chemical quantites
    IF ( air_chemistry )  THEN
       CALL cpu_log( log_point(83), '(chem advec+diff+prog)', 'start' )
!
!--    Loop over chemical species
       DO  lsp = 1, nvar
          CALL chem_prognostic_equations ( chem_species(lsp)%conc_p,           &
                                           chem_species(lsp)%conc,             &
                                           chem_species(lsp)%tconc_m,          &
                                           chem_species(lsp)%conc_pr_init,     &
                                           lsp )
       ENDDO

       CALL cpu_log( log_point(83), '(chem advec+diff+prog)', 'stop' )
    ENDIF   ! Chemicals equations


 END SUBROUTINE prognostic_equations_vector


 END MODULE prognostic_equations_mod
