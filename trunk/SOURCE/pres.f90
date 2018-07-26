!> @file pres.f90
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
! $Id: pres.f90 3016 2018-05-09 10:53:37Z Giersch $
! Dollar sign added before Id
! 
! 2696 2017-12-14 17:12:51Z kanani
! To avoid jumps while plotting w-profiles w level nzt+1 is set to w level nzt 
! after velocity modifications through the pressure solver were carried out
! 
! 2696 2017-12-14 17:12:51Z kanani
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! poismg_noopt modularized (MS)
! 
! 2298 2017-06-29 09:28:18Z raasch
! comment changed + some formatting
! 
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography and surface concept
! 
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC directives and related code removed
! 
! 2073 2016-11-30 14:34:05Z raasch
! openmp bugfix for calculation of new divergence
! 
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1932 2016-06-10 12:09:21Z suehring
! Initial version of purely vertical nesting introduced.
!
! 1931 2016-06-10 12:06:59Z suehring
! Rename multigrid into multigrid_noopt and multigrid_fast into multigrid
!
! 1929 2016-06-09 16:25:25Z suehring
! Bugfix: weight_substep for initial call, replace by local variable
! 
! 1918 2016-05-27 14:35:57Z raasch
! sum of divergence is also calculated when pres is called before the initial
! first time step,
! stearing is modified by using intermediate_timestep_count = 0 in order to
! determine if pres is called before the first initial timestep,
! bugfix: in case of Neumann conditions for pressure both at bottom and top,
!         mean vertical velocity is also removed for the first time step
! bugfix for calculating divergences
!
! 1908 2016-05-25 17:22:32Z suehring
! New divergence for output into RUN_CONTROL file is calculated only at last 
! Runge-Kutta step
!
! 1845 2016-04-08 08:29:13Z raasch
! nzb_2d replace by nzb_u|v_inner
! 
! 1799 2016-04-05 08:35:55Z gronemeier
! Bugfix: excluded third dimension from horizontal volume flow calculation
! 
! 1762 2016-02-25 12:31:13Z hellstea
! Introduction of nested domain feature
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
!
! 1575 2015-03-27 09:56:27Z raasch
! poismg_fast + respective module added, adjustments for psolver-queries
!
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
!
! 1318 2014-03-17 13:35:16Z raasch
! module interfaces removed
!
! 1306 2014-03-13 14:30:59Z raasch
! second argument removed from routine poisfft
!
! 1257 2013-11-08 15:18:40Z raasch
! openacc loop and loop vector clauses removed, independent clauses added,
! end parallel replaced by end parallel loop
!
! 1221 2013-09-10 08:59:13Z raasch
! openACC porting of reduction operations, loops for calculating d are
! using the rflags_s_inner multiply flag instead of the nzb_s_inner loop index
!
! 1212 2013-08-15 08:46:27Z raasch
! call of poisfft_hybrid removed
!
! 1117 2013-03-27 11:15:36Z suehring
! Bugfix in OpenMP parallelization.
!
! 1113 2013-03-10 02:48:14Z raasch
! GPU-porting of several loops, some loops rearranged
!
! 1111 2013-03-08 23:54:10Z
! openACC statements added,
! ibc_p_b = 2 removed
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1003 2012-09-14 14:35:53Z raasch
! adjustment of array tend for cases with unequal subdomain sizes removed
!
! Revision 1.1  1997/07/24 11:24:44  raasch
! Initial revision
!
!
! Description:
! ------------
!> Compute the divergence of the provisional velocity field. Solve the Poisson
!> equation for the perturbation pressure. Compute the final velocities using
!> this perturbation pressure. Compute the remaining divergence.
!------------------------------------------------------------------------------!
 SUBROUTINE pres
 

    USE arrays_3d,                                                             &
        ONLY:  d, ddzu, ddzu_pres, ddzw, dzw, p, p_loc, rho_air, rho_air_zw,   &
               tend, u, v, w

    USE control_parameters,                                                    &
        ONLY:  bc_lr_cyc, bc_ns_cyc, conserve_volume_flow, coupling_mode,      &
               dt_3d, gathered_size, ibc_p_b, ibc_p_t,                         &
               intermediate_timestep_count, intermediate_timestep_count_max,   &
               mg_switch_to_pe0_level, nest_domain, outflow_l, outflow_n,      &
               outflow_r, outflow_s, psolver, subdomain_size, topography,      &
               volume_flow, volume_flow_area, volume_flow_initial

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s

    USE grid_variables,                                                        &
        ONLY:  ddx, ddy

    USE indices,                                                               &
        ONLY:  nbgp, ngp_2dh_outer, nx, nxl, nxlg, nxl_mg, nxr, nxrg, nxr_mg,  &
               ny, nys, nysg, nys_mg, nyn, nyng, nyn_mg, nzb, nzt, nzt_mg,     &
               wall_flags_0

    USE kinds

    USE pegrid
    
    USE pmc_interface,                                                         &
        ONLY:  nesting_mode 

    USE poisfft_mod,                                                           &
        ONLY:  poisfft

    USE poismg_mod

    USE poismg_noopt_mod

    USE statistics,                                                            &
        ONLY:  statistic_regions, sums_divnew_l, sums_divold_l, weight_pres,   &
               weight_substep

    USE surface_mod,                                                           &
        ONLY :  bc_h

    IMPLICIT NONE

    INTEGER(iwp) ::  i              !<
    INTEGER(iwp) ::  j              !<
    INTEGER(iwp) ::  k              !<
    INTEGER(iwp) ::  m              !<

    REAL(wp)     ::  ddt_3d         !<
    REAL(wp)     ::  d_weight_pres  !<
    REAL(wp)     ::  localsum       !<
    REAL(wp)     ::  threadsum      !<
    REAL(wp)     ::  weight_pres_l  !<
    REAL(wp)     ::  weight_substep_l !<

    REAL(wp), DIMENSION(1:3)   ::  volume_flow_l       !<
    REAL(wp), DIMENSION(1:3)   ::  volume_flow_offset  !<
    REAL(wp), DIMENSION(1:nzt) ::  w_l                 !<
    REAL(wp), DIMENSION(1:nzt) ::  w_l_l               !<

    LOGICAL :: nest_domain_nvn      !<


    CALL cpu_log( log_point(8), 'pres', 'start' )


!
!-- Calculate quantities to be used locally
    ddt_3d = 1.0_wp / dt_3d
    IF ( intermediate_timestep_count == 0 )  THEN
!
!--    If pres is called before initial time step
       weight_pres_l    = 1.0_wp
       d_weight_pres    = 1.0_wp
       weight_substep_l = 1.0_wp
    ELSE
       weight_pres_l    = weight_pres(intermediate_timestep_count)
       d_weight_pres    = 1.0_wp / weight_pres(intermediate_timestep_count)
       weight_substep_l = weight_substep(intermediate_timestep_count)
    ENDIF

!
!-- Multigrid method expects array d to have one ghost layer.
!-- 
    IF ( psolver(1:9) == 'multigrid' )  THEN
     
       DEALLOCATE( d )
       ALLOCATE( d(nzb:nzt+1,nys-1:nyn+1,nxl-1:nxr+1) ) 

!
!--    Since p is later used to hold the weighted average of the substeps, it
!--    cannot be used in the iterative solver. Therefore, its initial value is
!--    stored on p_loc, which is then iteratively advanced in every substep.
       IF ( intermediate_timestep_count <= 1 )  THEN
          DO  i = nxl-1, nxr+1
             DO  j = nys-1, nyn+1
                DO  k = nzb, nzt+1
                   p_loc(k,j,i) = p(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       
    ELSEIF ( psolver == 'sor'  .AND.  intermediate_timestep_count <= 1 )  THEN

!
!--    Since p is later used to hold the weighted average of the substeps, it
!--    cannot be used in the iterative solver. Therefore, its initial value is
!--    stored on p_loc, which is then iteratively advanced in every substep.
       p_loc = p

    ENDIF

!
!-- Conserve the volume flow at the outflow in case of non-cyclic lateral
!-- boundary conditions
!-- WARNING: so far, this conservation does not work at the left/south
!--          boundary if the topography at the inflow differs from that at the
!--          outflow! For this case, volume_flow_area needs adjustment!
!
!-- Left/right
    IF ( conserve_volume_flow  .AND.  ( outflow_l .OR. outflow_r ) )  THEN

       volume_flow(1)   = 0.0_wp
       volume_flow_l(1) = 0.0_wp

       IF ( outflow_l )  THEN
          i = 0
       ELSEIF ( outflow_r )  THEN
          i = nx+1
       ENDIF

       DO  j = nys, nyn
!
!--       Sum up the volume flow through the south/north boundary
          DO  k = nzb+1, nzt
             volume_flow_l(1) = volume_flow_l(1) + u(k,j,i) * dzw(k)           &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,i), 1 )  &
                                            )
          ENDDO
       ENDDO

#if defined( __parallel )   
       IF ( collective_wait )  CALL MPI_BARRIER( comm1dy, ierr )
       CALL MPI_ALLREDUCE( volume_flow_l(1), volume_flow(1), 1, MPI_REAL, &
                           MPI_SUM, comm1dy, ierr )   
#else
       volume_flow = volume_flow_l  
#endif
       volume_flow_offset(1) = ( volume_flow_initial(1) - volume_flow(1) ) &
                               / volume_flow_area(1)

       DO  j = nysg, nyng
          DO  k = nzb+1, nzt
             u(k,j,i) = u(k,j,i) + volume_flow_offset(1)                       &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,i), 1 )  &
                                            )
          ENDDO
       ENDDO

    ENDIF

!
!-- South/north
    IF ( conserve_volume_flow  .AND.  ( outflow_n .OR. outflow_s ) )  THEN

       volume_flow(2)   = 0.0_wp
       volume_flow_l(2) = 0.0_wp

       IF ( outflow_s )  THEN
          j = 0
       ELSEIF ( outflow_n )  THEN
          j = ny+1
       ENDIF

       DO  i = nxl, nxr
!
!--       Sum up the volume flow through the south/north boundary
          DO  k = nzb+1, nzt
             volume_flow_l(2) = volume_flow_l(2) + v(k,j,i) * dzw(k)           &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,i), 2 )  &
                                            )
          ENDDO
       ENDDO

#if defined( __parallel )   
       IF ( collective_wait )  CALL MPI_BARRIER( comm1dx, ierr )
       CALL MPI_ALLREDUCE( volume_flow_l(2), volume_flow(2), 1, MPI_REAL, &
                           MPI_SUM, comm1dx, ierr )   
#else
       volume_flow = volume_flow_l  
#endif
       volume_flow_offset(2) = ( volume_flow_initial(2) - volume_flow(2) )    &
                               / volume_flow_area(2)

       DO  i = nxlg, nxrg
          DO  k = nzb+1, nzt
             v(k,j,i) = v(k,j,i) + volume_flow_offset(2)                       &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,i), 2 )  &
                                            )
          ENDDO
       ENDDO

    ENDIF

!
!-- Remove mean vertical velocity in case that Neumann conditions are
!-- used both at bottom and top boundary, and if not a nested domain in a 
!-- normal nesting run. In case of vertical nesting, this must be done.
!-- Therefore an auxiliary logical variable nest_domain_nvn is used here, and 
!-- nvn stands for non-vertical nesting. 
!-- This cannot be done before the first initial time step because ngp_2dh_outer
!-- is not yet known then.
    nest_domain_nvn = nest_domain
    IF ( nest_domain .AND. nesting_mode == 'vertical' )  THEN
       nest_domain_nvn = .FALSE.
    ENDIF

    IF ( ibc_p_b == 1  .AND.  ibc_p_t == 1  .AND.                               &
         .NOT. nest_domain_nvn  .AND. intermediate_timestep_count /= 0 )        &
    THEN
       w_l = 0.0_wp;  w_l_l = 0.0_wp
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                w_l_l(k) = w_l_l(k) + w(k,j,i)                                 &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,i), 3 )  &
                                            )
             ENDDO
          ENDDO
       ENDDO
#if defined( __parallel )   
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( w_l_l(1), w_l(1), nzt, MPI_REAL, MPI_SUM, &
                           comm2d, ierr )
#else
       w_l = w_l_l
#endif
       DO  k = 1, nzt
          w_l(k) = w_l(k) / ngp_2dh_outer(k,0)
       ENDDO
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb+1, nzt
                w(k,j,i) = w(k,j,i) - w_l(k)                                   &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,i), 3 )  &
                                            )
             ENDDO
          ENDDO
       ENDDO
    ENDIF

!
!-- Compute the divergence of the provisional velocity field.
    CALL cpu_log( log_point_s(1), 'divergence', 'start' )

    IF ( psolver(1:9) == 'multigrid' )  THEN
       !$OMP PARALLEL DO SCHEDULE( STATIC ) PRIVATE (i,j,k)
       DO  i = nxl-1, nxr+1
          DO  j = nys-1, nyn+1
             DO  k = nzb, nzt+1
                d(k,j,i) = 0.0_wp
             ENDDO
          ENDDO
       ENDDO
    ELSE
       !$OMP PARALLEL DO SCHEDULE( STATIC ) PRIVATE (i,j,k)
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                d(k,j,i) = 0.0_wp
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    localsum  = 0.0_wp
    threadsum = 0.0_wp

#if defined( __ibm )
    !$OMP PARALLEL PRIVATE (i,j,k) FIRSTPRIVATE(threadsum) REDUCTION(+:localsum)
    !$OMP DO SCHEDULE( STATIC )
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             d(k,j,i) = ( ( u(k,j,i+1) - u(k,j,i) ) * rho_air(k) * ddx +       &
                          ( v(k,j+1,i) - v(k,j,i) ) * rho_air(k) * ddy +       &
                          ( w(k,j,i)   * rho_air_zw(k) -                       &
                            w(k-1,j,i) * rho_air_zw(k-1) ) * ddzw(k)           &
                        ) * ddt_3d * d_weight_pres                             &
                                   * MERGE( 1.0_wp, 0.0_wp,                    &
                                            BTEST( wall_flags_0(k,j,i), 0 )    &
                                          )
          ENDDO
!
!--       Compute possible PE-sum of divergences for flow_statistics
          DO  k = nzb+1, nzt
             threadsum = threadsum + ABS( d(k,j,i) )                           &
                                   * MERGE( 1.0_wp, 0.0_wp,                    &
                                            BTEST( wall_flags_0(k,j,i), 0 )    &
                                          )
          ENDDO

       ENDDO
    ENDDO

    IF ( intermediate_timestep_count == intermediate_timestep_count_max  .OR.  &
         intermediate_timestep_count == 0 )  THEN
       localsum = localsum + threadsum * dt_3d * weight_pres_l
    ENDIF
    !$OMP END PARALLEL
#else

    !$OMP PARALLEL PRIVATE (i,j,k)
    !$OMP DO SCHEDULE( STATIC )
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = 1, nzt
             d(k,j,i) = ( ( u(k,j,i+1) - u(k,j,i) ) * rho_air(k) * ddx +       &
                          ( v(k,j+1,i) - v(k,j,i) ) * rho_air(k) * ddy +       &
                          ( w(k,j,i)   * rho_air_zw(k) -                       &
                            w(k-1,j,i) * rho_air_zw(k-1) ) * ddzw(k)           &
                        ) * ddt_3d * d_weight_pres                             &
                                   * MERGE( 1.0_wp, 0.0_wp,                    &
                                            BTEST( wall_flags_0(k,j,i), 0 )    &
                                          )     
          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL

!
!-- Compute possible PE-sum of divergences for flow_statistics. Carry out
!-- computation only at last Runge-Kutta substep.
    IF ( intermediate_timestep_count == intermediate_timestep_count_max  .OR.  &
         intermediate_timestep_count == 0 )  THEN
       !$OMP PARALLEL PRIVATE (i,j,k) FIRSTPRIVATE(threadsum) REDUCTION(+:localsum)
       !$OMP DO SCHEDULE( STATIC )
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                threadsum = threadsum + ABS( d(k,j,i) )
             ENDDO
          ENDDO
       ENDDO
       localsum = localsum + threadsum * dt_3d * weight_pres_l
       !$OMP END PARALLEL
    ENDIF
#endif

!
!-- For completeness, set the divergence sum of all statistic regions to those
!-- of the total domain
    IF ( intermediate_timestep_count == intermediate_timestep_count_max  .OR.  &
         intermediate_timestep_count == 0 )  THEN
       sums_divold_l(0:statistic_regions) = localsum
    ENDIF

    CALL cpu_log( log_point_s(1), 'divergence', 'stop' )

!
!-- Compute the pressure perturbation solving the Poisson equation
    IF ( psolver == 'poisfft' )  THEN

!
!--    Solve Poisson equation via FFT and solution of tridiagonal matrices
       CALL poisfft( d )

!
!--    Store computed perturbation pressure and set boundary condition in
!--    z-direction
       !$OMP PARALLEL DO PRIVATE (i,j,k)
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                tend(k,j,i) = d(k,j,i)
             ENDDO
          ENDDO
       ENDDO

!
!--    Bottom boundary:
!--    This condition is only required for internal output. The pressure
!--    gradient (dp(nzb+1)-dp(nzb))/dz is not used anywhere else.
       IF ( ibc_p_b == 1 )  THEN
!
!--       Neumann (dp/dz = 0). Using surfae data type, first for non-natural 
!--       surfaces, then for natural and urban surfaces
!--       Upward facing
          !$OMP PARALLEL DO PRIVATE( i, j, k )
          DO  m = 1, bc_h(0)%ns
             i = bc_h(0)%i(m)
             j = bc_h(0)%j(m)
             k = bc_h(0)%k(m)
             tend(k-1,j,i) = tend(k,j,i)
          ENDDO
!
!--       Downward facing
          !$OMP PARALLEL DO PRIVATE( i, j, k )
          DO  m = 1, bc_h(1)%ns
             i = bc_h(1)%i(m)
             j = bc_h(1)%j(m)
             k = bc_h(1)%k(m)
             tend(k+1,j,i) = tend(k,j,i)
          ENDDO

       ELSE
!
!--       Dirichlet. Using surface data type, first for non-natural 
!--       surfaces, then for natural and urban surfaces
!--       Upward facing
          !$OMP PARALLEL DO PRIVATE( i, j, k )
          DO  m = 1, bc_h(0)%ns
             i = bc_h(0)%i(m)
             j = bc_h(0)%j(m)
             k = bc_h(0)%k(m)
             tend(k-1,j,i) = 0.0_wp
          ENDDO
!
!--       Downward facing
          !$OMP PARALLEL DO PRIVATE( i, j, k )
          DO  m = 1, bc_h(1)%ns
             i = bc_h(1)%i(m)
             j = bc_h(1)%j(m)
             k = bc_h(1)%k(m)
             tend(k+1,j,i) = 0.0_wp
          ENDDO

       ENDIF

!
!--    Top boundary
       IF ( ibc_p_t == 1 )  THEN
!
!--       Neumann
          !$OMP PARALLEL DO PRIVATE (i,j,k)
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                tend(nzt+1,j,i) = tend(nzt,j,i)
             ENDDO
          ENDDO

       ELSE
!
!--       Dirichlet
          !$OMP PARALLEL DO PRIVATE (i,j,k)
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                tend(nzt+1,j,i) = 0.0_wp
             ENDDO
          ENDDO

       ENDIF

!
!--    Exchange boundaries for p
       CALL exchange_horiz( tend, nbgp )
      
    ELSEIF ( psolver == 'sor' )  THEN

!
!--    Solve Poisson equation for perturbation pressure using SOR-Red/Black
!--    scheme
       CALL sor( d, ddzu_pres, ddzw, p_loc )
       tend = p_loc

    ELSEIF ( psolver(1:9) == 'multigrid' )  THEN

!
!--    Solve Poisson equation for perturbation pressure using Multigrid scheme,
!--    array tend is used to store the residuals.

!--    If the number of grid points of the gathered grid, which is collected
!--    on PE0, is larger than the number of grid points of an PE, than array
!--    tend will be enlarged. 
       IF ( gathered_size > subdomain_size )  THEN
          DEALLOCATE( tend )
          ALLOCATE( tend(nzb:nzt_mg(mg_switch_to_pe0_level)+1,nys_mg(          &
                    mg_switch_to_pe0_level)-1:nyn_mg(mg_switch_to_pe0_level)+1,&
                    nxl_mg(mg_switch_to_pe0_level)-1:nxr_mg(                   &
                    mg_switch_to_pe0_level)+1) )
       ENDIF

       IF ( psolver == 'multigrid' )  THEN
          CALL poismg( tend )
       ELSE
          CALL poismg_noopt( tend )
       ENDIF

       IF ( gathered_size > subdomain_size )  THEN
          DEALLOCATE( tend )
          ALLOCATE( tend(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ENDIF

!
!--    Restore perturbation pressure on tend because this array is used
!--    further below to correct the velocity fields
       DO  i = nxl-1, nxr+1
          DO  j = nys-1, nyn+1
             DO  k = nzb, nzt+1
                tend(k,j,i) = p_loc(k,j,i)
             ENDDO
          ENDDO
       ENDDO

    ENDIF

!
!-- Store perturbation pressure on array p, used for pressure data output.
!-- Ghost layers are added in the output routines (except sor-method: see below)
    IF ( intermediate_timestep_count <= 1 )  THEN
       !$OMP PARALLEL PRIVATE (i,j,k)
       !$OMP DO
       DO  i = nxl-1, nxr+1
          DO  j = nys-1, nyn+1
             DO  k = nzb, nzt+1
                p(k,j,i) = tend(k,j,i) * &
                           weight_substep_l
             ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL

    ELSEIF ( intermediate_timestep_count > 1 )  THEN
       !$OMP PARALLEL PRIVATE (i,j,k)
       !$OMP DO
       DO  i = nxl-1, nxr+1
          DO  j = nys-1, nyn+1
             DO  k = nzb, nzt+1
                p(k,j,i) = p(k,j,i) + tend(k,j,i) * &
                           weight_substep_l
             ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL

    ENDIF
       
!
!-- SOR-method needs ghost layers for the next timestep
    IF ( psolver == 'sor' )  CALL exchange_horiz( p, nbgp )

!
!-- Correction of the provisional velocities with the current perturbation 
!-- pressure just computed
    IF ( conserve_volume_flow  .AND.  ( bc_lr_cyc .OR. bc_ns_cyc ) )  THEN
       volume_flow_l(1) = 0.0_wp
       volume_flow_l(2) = 0.0_wp
    ENDIF

    !$OMP PARALLEL PRIVATE (i,j,k)
    !$OMP DO
    DO  i = nxl, nxr   
       DO  j = nys, nyn

          DO  k = nzb+1, nzt
             w(k,j,i) = w(k,j,i) - dt_3d *                                     &
                           ( tend(k+1,j,i) - tend(k,j,i) ) * ddzu(k+1)         &
                                     * weight_pres_l                           &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,i), 3 )  &
                                            )
          ENDDO

          DO  k = nzb+1, nzt
             u(k,j,i) = u(k,j,i) - dt_3d *                                     &
                           ( tend(k,j,i) - tend(k,j,i-1) ) * ddx               &
                                     * weight_pres_l                           &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,i), 1 )  &
                                            )
          ENDDO

          DO  k = nzb+1, nzt
             v(k,j,i) = v(k,j,i) - dt_3d *                                     &
                           ( tend(k,j,i) - tend(k,j-1,i) ) * ddy               &
                                     * weight_pres_l                           &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,i), 2 )  &
                                            )
          ENDDO                                                         

       ENDDO
    ENDDO
    !$OMP END PARALLEL

!
!-- The vertical velocity is not set to zero at nzt + 1 for nested domains
!-- Instead it is set to the values of nzt (see routine vnest_boundary_conds 
!-- or pmci_interp_tril_t) BEFORE calling the pressure solver. To avoid jumps
!-- while plotting profiles w at the top has to be set to the values in the 
!-- height nzt after above modifications. Hint: w level nzt+1 does not impact
!-- results. 
    IF (nest_domain .OR. coupling_mode == 'vnested_fine') THEN
       w(nzt+1,:,:) = w(nzt,:,:)
    ENDIF

!
!-- Sum up the volume flow through the right and north boundary
    IF ( conserve_volume_flow  .AND.  bc_lr_cyc  .AND.  bc_ns_cyc  .AND.       &
         nxr == nx )  THEN

       !$OMP PARALLEL PRIVATE (j,k)
       !$OMP DO
       DO  j = nys, nyn
          !$OMP CRITICAL
          DO  k = nzb+1, nzt
             volume_flow_l(1) = volume_flow_l(1) + u(k,j,nxr) * dzw(k)         &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,nxr), 1 )&
                                            )
          ENDDO
          !$OMP END CRITICAL
       ENDDO
       !$OMP END PARALLEL

    ENDIF

    IF ( conserve_volume_flow  .AND.  bc_ns_cyc  .AND.  bc_lr_cyc  .AND.       &
         nyn == ny )  THEN

       !$OMP PARALLEL PRIVATE (i,k)
       !$OMP DO
       DO  i = nxl, nxr
          !$OMP CRITICAL
          DO  k = nzb+1, nzt
             volume_flow_l(2) = volume_flow_l(2) + v(k,nyn,i) * dzw(k)         &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,nyn,i), 2 )&
                                            )
           ENDDO
          !$OMP END CRITICAL
       ENDDO
       !$OMP END PARALLEL

    ENDIF
    
!
!-- Conserve the volume flow
    IF ( conserve_volume_flow  .AND.  ( bc_lr_cyc  .AND.  bc_ns_cyc ) )  THEN

#if defined( __parallel )   
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( volume_flow_l(1), volume_flow(1), 2, MPI_REAL, &
                           MPI_SUM, comm2d, ierr )  
#else
       volume_flow = volume_flow_l  
#endif   

       volume_flow_offset(1:2) = ( volume_flow_initial(1:2) - volume_flow(1:2) ) / &
                            volume_flow_area(1:2)

       !$OMP PARALLEL PRIVATE (i,j,k)
       !$OMP DO
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                u(k,j,i) = u(k,j,i) + volume_flow_offset(1)                    &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,i), 1 )  &
                                            )
             ENDDO
             DO  k = nzb+1, nzt
                v(k,j,i) = v(k,j,i) + volume_flow_offset(2)                    &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,i), 2 )  &
                                            )
             ENDDO
          ENDDO
       ENDDO

       !$OMP END PARALLEL

    ENDIF

!
!-- Exchange of boundaries for the velocities
    CALL exchange_horiz( u, nbgp )
    CALL exchange_horiz( v, nbgp )
    CALL exchange_horiz( w, nbgp )

!
!-- Compute the divergence of the corrected velocity field,
!-- A possible PE-sum is computed in flow_statistics. Carry out computation
!-- only at last Runge-Kutta step.
    IF ( intermediate_timestep_count == intermediate_timestep_count_max  .OR.  &
         intermediate_timestep_count == 0 )  THEN
       CALL cpu_log( log_point_s(1), 'divergence', 'start' )
       sums_divnew_l = 0.0_wp

!
!--    d must be reset to zero because it can contain nonzero values below the
!--    topography
       IF ( topography /= 'flat' )  d = 0.0_wp

       localsum  = 0.0_wp
       threadsum = 0.0_wp

       !$OMP PARALLEL PRIVATE (i,j,k) FIRSTPRIVATE(threadsum) REDUCTION(+:localsum)
#if defined( __ibm )
       !$OMP DO SCHEDULE( STATIC )
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
             d(k,j,i) = ( ( u(k,j,i+1) - u(k,j,i) ) * rho_air(k) * ddx +       &
                          ( v(k,j+1,i) - v(k,j,i) ) * rho_air(k) * ddy +       &
                          ( w(k,j,i)   * rho_air_zw(k) -                       &
                            w(k-1,j,i) * rho_air_zw(k-1) ) * ddzw(k)           &
                        ) * MERGE( 1.0_wp, 0.0_wp,                             &
                                   BTEST( wall_flags_0(k,j,i), 0 )             &
                                 )
             ENDDO
             DO  k = nzb+1, nzt
                threadsum = threadsum + ABS( d(k,j,i) )
             ENDDO
          ENDDO
       ENDDO
#else
       !$OMP DO SCHEDULE( STATIC )
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                d(k,j,i) = ( ( u(k,j,i+1) - u(k,j,i) ) * rho_air(k) * ddx +    &
                             ( v(k,j+1,i) - v(k,j,i) ) * rho_air(k) * ddy +    &
                             ( w(k,j,i)   * rho_air_zw(k) -                    &
                               w(k-1,j,i) * rho_air_zw(k-1) ) * ddzw(k)        &
                           ) * MERGE( 1.0_wp, 0.0_wp,                          &
                                   BTEST( wall_flags_0(k,j,i), 0 )             &
                                    )
             ENDDO
          ENDDO
       ENDDO
!
!--    Compute possible PE-sum of divergences for flow_statistics
       !$OMP DO SCHEDULE( STATIC )
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                threadsum = threadsum + ABS( d(k,j,i) )
             ENDDO
          ENDDO
       ENDDO
#endif

       localsum = localsum + threadsum
       !$OMP END PARALLEL

!
!--    For completeness, set the divergence sum of all statistic regions to those
!--    of the total domain
       sums_divnew_l(0:statistic_regions) = localsum

       CALL cpu_log( log_point_s(1), 'divergence', 'stop' )

    ENDIF

    CALL cpu_log( log_point(8), 'pres', 'stop' )


 END SUBROUTINE pres
