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
               w, w_p, alpha_T, beta_S
    USE kinds

    USE control_parameters,                                                    &
        ONLY:  air_chemistry, call_microphysics_at_all_substeps,               &
               cloud_physics, cloud_top_radiation,         &
               dp_external, dp_level_ind_b, dp_smooth_factor, dpdxy, dt_3d,    &
               humidity, idealized_diurnal, g,                                 &
               inflow_l, intermediate_timestep_count,                          &
               intermediate_timestep_count_max, large_scale_forcing,           &
               large_scale_subsidence, microphysics_morrison,                  &
               microphysics_seifert, microphysics_sat_adjust, neutral, nudging,&
               ocean, outflow_l, outflow_s, passive_scalar, plant_canopy,      &
               prho_reference, prho_reference,                                 &
               prho_reference, pt_reference, pt_reference, pt_reference,       &
               scalar_advec, scalar_advec, simulated_time, sloping_surface,    &
               timestep_scheme, tsc, use_subsidence_tendencies,                &
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

    USE buoyancy_mod,                                                          &
        ONLY:  buoyancy

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
    USE statistics,                                                            &
        ONLY:  hom

    USE turbulence_closure_mod,                                                &
        ONLY:  tcm_prognostic

    USE surface_mod,                                                           &
        ONLY :  surf_def_h, surf_def_v

    USE constants

    USE stokes_force_mod,                                                      &
        ONLY: stokes_force_uvw, stokes_force_s

    PRIVATE
    PUBLIC prognostic_equations_vector
    INTERFACE prognostic_equations_vector
       MODULE PROCEDURE prognostic_equations_vector
    END INTERFACE prognostic_equations_vector


 CONTAINS


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
!-- u-velocity component
    CALL cpu_log( log_point(5), 'u-equation', 'start' )

    tend = 0.0_wp

    call cpu_log( log_point(43), 'u advec', 'start')
    CALL advec_u_ws
    call cpu_log( log_point(43), 'u advec', 'stop')

    call cpu_log( log_point(44), 'u diffusion', 'start')
    CALL diffusion_u
    call cpu_log( log_point(44), 'u diffusion', 'stop')

    call cpu_log( log_point(45), 'u coriolis', 'start')
    CALL coriolis( 1 )
    call cpu_log( log_point(45), 'u coriolis', 'stop')

    call cpu_log( log_point(46), 'u stokes', 'start')
    !
!-- If required, compute Stokes forces
    IF ( ocean .AND. stokes_force ) THEN
       CALL stokes_force_uvw( 1 )
    ENDIF
    call cpu_log( log_point(46), 'u stokes', 'stop')
!
!-- External pressure gradient
    IF ( dp_external )  THEN
       DO  i = nxlu, nxr
          DO  j = nys, nyn
             DO  k = dp_level_ind_b+1, nzt
                tend(k,j,i) = tend(k,j,i) - dpdxy(1) * dp_smooth_factor(k)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

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
    CALL advec_v_ws
    CALL diffusion_v
    CALL coriolis( 2 )

!
!-- If required, compute Stokes forces
    IF ( ocean .AND. stokes_force ) THEN
       CALL stokes_force_uvw( 2 )
    ENDIF
!
!-- External pressure gradient
    IF ( dp_external )  THEN
       DO  i = nxl, nxr
          DO  j = nysv, nyn
             DO  k = dp_level_ind_b+1, nzt
                tend(k,j,i) = tend(k,j,i) - dpdxy(2) * dp_smooth_factor(k)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
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
    CALL advec_w_ws
    CALL diffusion_w
    CALL coriolis( 3 )

    call cpu_log( log_point(47), 'w buoy', 'start')
    CALL buoyancy( rho_ocean, 3 )
    call cpu_log( log_point(47), 'w buoy', 'stop')
    !
!-- If required, compute Stokes forces
    IF ( ocean .AND. stokes_force ) THEN
       CALL stokes_force_uvw( 3 )
    ENDIF
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

    CALL cpu_log( log_point(13), 'pt-equation', 'start' )

!
!--    pt-tendency terms with communication
       sbt = tsc(2)
!
       tend = 0.0_wp
       call cpu_log( log_point(40), 'pt advection', 'start')
       CALL advec_s_ws( pt, 'pt' )
       call cpu_log( log_point(40), 'pt advection', 'stop')

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

       call cpu_log( log_point(41), 'pt diffusion' , 'start')

       CALL diffusion_s( pt,                                                   &
                         surf_def_h(0)%shf, surf_def_h(1)%shf,                 &
                         surf_def_h(2)%shf,                                    &
                         surf_def_v(0)%shf, surf_def_v(1)%shf,                 &
                         surf_def_v(2)%shf, surf_def_v(3)%shf,                 &
                         surf_def_h(2)%shf_sol )

       call cpu_log( log_point(41), 'pt diffusion' , 'stop')

       call cpu_log( log_point(42), 'pt stokes' ,'start')
!
!--    If required, compute Stokes-advection term
       IF ( ocean .AND. stokes_force ) THEN
          CALL stokes_force_s( pt )
       ENDIF

       call cpu_log( log_point(42), 'pt stokes', 'stop' )

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

       CALL cpu_log( log_point(37), 'sa-equation', 'start' )

!
!--    sa-tendency terms with communication
       sbt = tsc(2)
!
       tend = 0.0_wp
       CALL advec_s_ws( sa, 'sa' )

       CALL diffusion_s( sa,                                                   &
                         surf_def_h(0)%sasws, surf_def_h(1)%sasws,             &
                         surf_def_h(2)%sasws,                                  &
                         surf_def_v(0)%sasws, surf_def_v(1)%sasws,             &
                         surf_def_v(2)%sasws, surf_def_v(3)%sasws)

!
!--    If required, compute Stokes-advection term
       IF ( stokes_force ) THEN
          CALL stokes_force_s( sa )
       ENDIF

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

    CALL tcm_prognostic()

 END SUBROUTINE prognostic_equations_vector


 END MODULE prognostic_equations_mod
