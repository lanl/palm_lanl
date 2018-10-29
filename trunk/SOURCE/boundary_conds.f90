!> @file boundary_conds.f90
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
! 2018-10-25 cbegeman
! Treatment of dirichlet bottom boundary conditions for salinity
! 
! Former revisions:
! -----------------
! $Id: boundary_conds.f90 2938 2018-03-27 15:52:42Z suehring $
! Set boundary condition for TKE and TKE dissipation rate in case of nesting
! and if parent model operates in RANS mode but child model in LES mode.
! mode
! 
! 2793 2018-02-07 10:54:33Z suehring
! Removed preprocessor directive __chem
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Adjust boundary conditions for e and diss in case of TKE-e closure (TG)
! Implementation of chemistry module (FK)
! 
! 2569 2017-10-20 11:54:42Z kanani
! Removed redundant code for ibc_s_b=1 and ibc_q_b=1
! 
! 2365 2017-08-21 14:59:59Z kanani
! Vertical grid nesting implemented: exclude setting vertical velocity to zero 
! on fine grid (SadiqHuq)
! 
! 2320 2017-07-21 12:47:43Z suehring
! Remove unused control parameter large_scale_forcing from only-list
! 
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison' 
! includes two more prognostic equations for cloud drop concentration (nc)  
! and cloud water content (qc). 
! 
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Set boundary conditions on topography top using flag method.
! 
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC directives removed
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1992 2016-08-12 15:14:59Z suehring
! Adjustments for top boundary condition for passive scalar
! 
! 1960 2016-07-12 16:34:24Z suehring
! Treat humidity and passive scalar separately
! 
! 1823 2016-04-07 08:57:52Z hoffmann
! Initial version of purely vertical nesting introduced.
!
! 1822 2016-04-07 07:49:42Z hoffmann
! icloud_scheme removed. microphyisics_seifert added.
!
! 1764 2016-02-28 12:45:19Z raasch
! index bug for u_p at left outflow removed
!
! 1762 2016-02-25 12:31:13Z hellstea
! Introduction of nested domain feature
!
! 1742 2016-01-13 09:50:06Z raasch
! bugfix for outflow Neumann boundary conditions at bottom and top
!
! 1717 2015-11-11 15:09:47Z raasch
! Bugfix: index error in outflow conditions for left boundary
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
! 
! 1410 2014-05-23 12:16:18Z suehring
! Bugfix: set dirichlet boundary condition for passive_scalar at model domain 
! top 
!
! 1399 2014-05-07 11:16:25Z heinze
! Bugfix: set inflow boundary conditions also if no humidity or passive_scalar
! is used.
! 
! 1398 2014-05-07 11:15:00Z heinze
! Dirichlet-condition at the top for u and v changed to u_init and v_init also 
! for large_scale_forcing
! 
! 1380 2014-04-28 12:40:45Z heinze
! Adjust Dirichlet-condition at the top for pt in case of nudging
! 
! 1361 2014-04-16 15:17:48Z hoffmann
! Bottom and top boundary conditions of rain water content (qr) and 
! rain drop concentration (nr) changed to Dirichlet
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
! loop independent clauses added
!
! 1241 2013-10-30 11:36:58Z heinze
! Adjust ug and vg at each timestep in case of large_scale_forcing
!
! 1159 2013-05-21 11:58:22Z fricke
! Bugfix: Neumann boundary conditions for the velocity components at the
! outflow are in fact radiation boundary conditions using the maximum phase
! velocity that ensures numerical stability (CFL-condition).
! Hence, logical operator use_cmax is now used instead of bc_lr_dirneu/_neudir.
! Bugfix: In case of use_cmax at the outflow, u, v, w are replaced by
! u_p, v_p, w_p  
!
! 1115 2013-03-26 18:16:16Z hoffmann
! boundary conditions of two-moment cloud scheme are restricted to Neumann-
! boundary-conditions
!
! 1113 2013-03-10 02:48:14Z raasch
! GPU-porting
! dummy argument "range" removed
! Bugfix: wrong index in loops of radiation boundary condition
!
! 1053 2012-11-13 17:11:03Z hoffmann
! boundary conditions for the two new prognostic equations (nr, qr) of the 
! two-moment cloud scheme
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 996 2012-09-07 10:41:47Z raasch
! little reformatting
!
! 978 2012-08-09 08:28:32Z fricke
! Neumann boudnary conditions are added at the inflow boundary for the SGS-TKE.
! Outflow boundary conditions for the velocity components can be set to Neumann
! conditions or to radiation conditions with a horizontal averaged phase
! velocity. 
!
! 875 2012-04-02 15:35:15Z gryschka
! Bugfix in case of dirichlet inflow bc at the right or north boundary
!
! Revision 1.1  1997/09/12 06:21:34  raasch
! Initial revision
!
!
! Description:
! ------------
!> Boundary conditions for the prognostic quantities.
!> One additional bottom boundary condition is applied for the TKE (=(u*)**2)
!> in prandtl_fluxes. The cyclic lateral boundary conditions are implicitly
!> handled in routine exchange_horiz. Pressure boundary conditions are
!> explicitly set in routines pres, poisfft, poismg and sor.
!------------------------------------------------------------------------------!
 SUBROUTINE boundary_conds
 

    USE arrays_3d,                                                             &
        ONLY:  c_u, c_u_m, c_u_m_l, c_v, c_v_m, c_v_m_l, c_w, c_w_m, c_w_m_l,  &
               diss_p, dzu, e_p, nc_p, nr_p, pt, pt_p, q, q_p, qc_p, qr_p, s,  & 
               s_p, sa, sa_p, u, ug, u_init, u_m_l, u_m_n, u_m_r, u_m_s, u_p,  &
               v, vg, v_init, v_m_l, v_m_n, v_m_r, v_m_s, v_p,                 &
               w, w_p, w_m_l, w_m_n, w_m_r, w_m_s, pt_init

    USE chemistry_model_mod,                                                   &
        ONLY:  chem_boundary_conds 
              
    USE control_parameters,                                                    &
        ONLY:  air_chemistry, bc_pt_t_val, bc_q_t_val, bc_s_t_val,             &
               constant_diffusion, cloud_physics, coupling_mode, dt_3d,        &
               force_bound_l, force_bound_s, forcing, humidity,                &
               ibc_pt_b, ibc_pt_t, ibc_q_b, ibc_q_t, ibc_s_b, ibc_s_t,         &
               ibc_sa_b, ibc_sa_t, ibc_uv_b, ibc_uv_t, inflow_l, inflow_n,     &
               inflow_r, inflow_s, intermediate_timestep_count,                &
               microphysics_morrison, microphysics_seifert, nest_domain,       &
               nest_bound_l, nest_bound_n, nest_bound_r, nest_bound_s, nudging,&
               ocean, outflow_l, outflow_n, outflow_r, outflow_s,              &
               passive_scalar, rans_mode, rans_tke_e, tsc, use_cmax

    USE grid_variables,                                                        &
        ONLY:  ddx, ddy, dx, dy

    USE indices,                                                               &
        ONLY:  nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg,             &
               nzb, nzt, wall_flags_0

    USE kinds

    USE pegrid

    USE pmc_interface,                                                         &
        ONLY : nesting_mode, rans_mode_parent

    USE surface_mod,                                                           &
        ONLY :  bc_h

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< grid index x direction
    INTEGER(iwp) ::  j  !< grid index y direction
    INTEGER(iwp) ::  k  !< grid index z direction
    INTEGER(iwp) ::  kb !< variable to set respective boundary value, depends on facing. 
    INTEGER(iwp) ::  l  !< running index boundary type, for up- and downward-facing walls
    INTEGER(iwp) ::  m  !< running index surface elements

    REAL(wp)    ::  c_max !<
    REAL(wp)    ::  denom !<


!
!-- Bottom boundary 
    IF ( ibc_uv_b == 1 )  THEN
       u_p(nzb,:,:) = u_p(nzb+1,:,:)
       v_p(nzb,:,:) = v_p(nzb+1,:,:)
    ENDIF
!
!-- Set zero vertical velocity at topography top (l=0), or bottom (l=1) in case
!-- of downward-facing surfaces. 
    DO  l = 0, 1
!
!--    Set kb, for upward-facing surfaces value at topography top (k-1) is set,
!--    for downward-facing surfaces at topography bottom (k+1). 
       kb = MERGE( -1, 1, l == 0 )
       !$OMP PARALLEL DO PRIVATE( i, j, k )
       DO  m = 1, bc_h(l)%ns
          i = bc_h(l)%i(m)            
          j = bc_h(l)%j(m)
          k = bc_h(l)%k(m)
          w_p(k+kb,j,i) = 0.0_wp
       ENDDO
    ENDDO

!
!-- Top boundary. A nested domain ( ibc_uv_t = 3 ) does not require settings.
    IF ( ibc_uv_t == 0 )  THEN
        u_p(nzt+1,:,:) = u_init(nzt+1)
        v_p(nzt+1,:,:) = v_init(nzt+1)
    ELSEIF ( ibc_uv_t == 1 )  THEN
        u_p(nzt+1,:,:) = u_p(nzt,:,:)
        v_p(nzt+1,:,:) = v_p(nzt,:,:)
    ENDIF

!
!-- Vertical nesting: Vertical velocity not zero at the top of the fine grid
    IF (  .NOT.  nest_domain  .AND.                                            &
                 TRIM(coupling_mode) /= 'vnested_fine' )  THEN
       w_p(nzt:nzt+1,:,:) = 0.0_wp  !< nzt is not a prognostic level (but cf. pres)
    ENDIF

!
!-- Temperature at bottom and top boundary.
!-- In case of coupled runs (ibc_pt_b = 2) the temperature is given by
!-- the sea surface temperature of the coupled ocean model.
!-- Dirichlet
    IF ( ibc_pt_b == 0 )  THEN
       DO  l = 0, 1
!
!--       Set kb, for upward-facing surfaces value at topography top (k-1) is set,
!--       for downward-facing surfaces at topography bottom (k+1). 
          kb = MERGE( -1, 1, l == 0 )
          !$OMP PARALLEL DO PRIVATE( i, j, k )
          DO  m = 1, bc_h(l)%ns
             i = bc_h(l)%i(m)            
             j = bc_h(l)%j(m)
             k = bc_h(l)%k(m)
             pt_p(k+kb,j,i) = pt(k+kb,j,i)
          ENDDO
       ENDDO
!
!-- Neumann, zero-gradient
    ELSEIF ( ibc_pt_b == 1 )  THEN
       DO  l = 0, 1
!
!--       Set kb, for upward-facing surfaces value at topography top (k-1) is set,
!--       for downward-facing surfaces at topography bottom (k+1). 
          kb = MERGE( -1, 1, l == 0 )
          !$OMP PARALLEL DO PRIVATE( i, j, k )
          DO  m = 1, bc_h(l)%ns
             i = bc_h(l)%i(m)            
             j = bc_h(l)%j(m)
             k = bc_h(l)%k(m)
             pt_p(k+kb,j,i) = pt_p(k,j,i)
          ENDDO
       ENDDO
    ENDIF

!
!-- Temperature at top boundary
    IF ( ibc_pt_t == 0 )  THEN
        pt_p(nzt+1,:,:) = pt(nzt+1,:,:)
!
!--     In case of nudging adjust top boundary to pt which is
!--     read in from NUDGING-DATA
        IF ( nudging )  THEN
           pt_p(nzt+1,:,:) = pt_init(nzt+1)
        ENDIF
    ELSEIF ( ibc_pt_t == 1 )  THEN
        pt_p(nzt+1,:,:) = pt_p(nzt,:,:)
    ELSEIF ( ibc_pt_t == 2 )  THEN
        pt_p(nzt+1,:,:) = pt_p(nzt,:,:) + bc_pt_t_val * dzu(nzt+1)
    ENDIF

!
!-- Boundary conditions for TKE.
!-- Generally Neumann conditions with de/dz=0 are assumed.
    IF ( .NOT. constant_diffusion )  THEN

       IF ( .NOT. rans_tke_e )  THEN
          DO  l = 0, 1
!
!--         Set kb, for upward-facing surfaces value at topography top (k-1) is set,
!--         for downward-facing surfaces at topography bottom (k+1). 
             kb = MERGE( -1, 1, l == 0 )
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, bc_h(l)%ns
                i = bc_h(l)%i(m)            
                j = bc_h(l)%j(m)
                k = bc_h(l)%k(m)
                e_p(k+kb,j,i) = e_p(k,j,i)
             ENDDO
          ENDDO
       ENDIF

       IF ( .NOT. nest_domain )  THEN
          e_p(nzt+1,:,:) = e_p(nzt,:,:)
!
!--    Nesting case: if parent operates in RANS mode and child in LES mode,
!--    no TKE is transfered. This case, set Neumann conditions at lateral and 
!--    top child boundaries. 
!--    If not ( both either in RANS or in LES mode ), TKE boundary condition
!--    is treated in the nesting. 
       ELSE

          IF ( rans_mode_parent  .AND.  .NOT. rans_mode )  THEN



             e_p(nzt+1,:,:) = e_p(nzt,:,:)
             IF ( nest_bound_l )  e_p(:,:,nxl-1) = e_p(:,:,nxl)
             IF ( nest_bound_r )  e_p(:,:,nxr+1) = e_p(:,:,nxr)
             IF ( nest_bound_s )  e_p(:,nys-1,:) = e_p(:,nys,:)
             IF ( nest_bound_n )  e_p(:,nyn+1,:) = e_p(:,nyn,:)

          ENDIF
       ENDIF
    ENDIF

!
!-- Boundary conditions for TKE dissipation rate. 
    IF ( rans_tke_e .AND. .NOT. nest_domain )  THEN
       diss_p(nzt+1,:,:) = diss_p(nzt,:,:)
    ENDIF

!
!-- Boundary conditions for salinity
    IF ( ocean )  THEN

!-- Bottom boundary: Dirichlet condition.
       IF ( ibc_sa_b == 0 )  THEN
          DO  l = 0, 1
!
!--          Set kb, for upward-facing surfaces value at topography top (k-1) is set,
!--          for downward-facing surfaces at topography bottom (k+1). 
             kb = MERGE( -1, 1, l == 0 )
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, bc_h(l)%ns
                i = bc_h(l)%i(m)            
                j = bc_h(l)%j(m)
                k = bc_h(l)%k(m)
                sa_p(k+kb,j,i) = sa(k+kb,j,i)
             ENDDO
          ENDDO

!
!--    Bottom boundary: Neumann condition.
       ELSEIF ( ibc_sa_b == 1 )  THEN
          DO  l = 0, 1
!
!--          Set kb, for upward-facing surfaces value at topography top (k-1) is set,
!--          for downward-facing surfaces at topography bottom (k+1). 
             kb = MERGE( -1, 1, l == 0 )
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, bc_h(l)%ns
                i = bc_h(l)%i(m)            
                j = bc_h(l)%j(m)
                k = bc_h(l)%k(m)
                sa_p(k+kb,j,i) = sa_p(k,j,i)
             ENDDO
          ENDDO
       ENDIF
!
!--    Top boundary: Dirichlet or Neumann
       IF ( ibc_sa_t == 0 )  THEN
           sa_p(nzt+1,:,:) = sa(nzt+1,:,:)
       ELSEIF ( ibc_sa_t == 1 )  THEN
           sa_p(nzt+1,:,:) = sa_p(nzt,:,:)
       ENDIF

    ENDIF

!
!-- Boundary conditions for total water content,
!-- bottom and top boundary (see also temperature)
    IF ( humidity )  THEN
!
!--    Surface conditions for constant_humidity_flux
!--    Run loop over all non-natural and natural walls. Note, in wall-datatype
!--    the k coordinate belongs to the atmospheric grid point, therefore, set
!--    q_p at k-1
       IF ( ibc_q_b == 0 ) THEN

          DO  l = 0, 1
!
!--          Set kb, for upward-facing surfaces value at topography top (k-1) is set,
!--          for downward-facing surfaces at topography bottom (k+1). 
             kb = MERGE( -1, 1, l == 0 )
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, bc_h(l)%ns
                i = bc_h(l)%i(m)            
                j = bc_h(l)%j(m)
                k = bc_h(l)%k(m)
                q_p(k+kb,j,i) = q(k+kb,j,i)
             ENDDO
          ENDDO
          
       ELSE
          
          DO  l = 0, 1
!
!--          Set kb, for upward-facing surfaces value at topography top (k-1) is set,
!--          for downward-facing surfaces at topography bottom (k+1). 
             kb = MERGE( -1, 1, l == 0 )
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, bc_h(l)%ns
                i = bc_h(l)%i(m)            
                j = bc_h(l)%j(m)
                k = bc_h(l)%k(m)
                q_p(k+kb,j,i) = q_p(k,j,i)
             ENDDO
          ENDDO
       ENDIF
!
!--    Top boundary
       IF ( ibc_q_t == 0 ) THEN
          q_p(nzt+1,:,:) = q(nzt+1,:,:)
       ELSEIF ( ibc_q_t == 1 ) THEN
          q_p(nzt+1,:,:) = q_p(nzt,:,:) + bc_q_t_val * dzu(nzt+1)
       ENDIF

       IF ( cloud_physics  .AND.  microphysics_morrison )  THEN
!             
!--       Surface conditions cloud water (Dirichlet)
!--       Run loop over all non-natural and natural walls. Note, in wall-datatype
!--       the k coordinate belongs to the atmospheric grid point, therefore, set
!--       qr_p and nr_p at k-1
          !$OMP PARALLEL DO PRIVATE( i, j, k )
          DO  m = 1, bc_h(0)%ns
             i = bc_h(0)%i(m)            
             j = bc_h(0)%j(m)
             k = bc_h(0)%k(m)
             qc_p(k-1,j,i) = 0.0_wp
             nc_p(k-1,j,i) = 0.0_wp
          ENDDO
!
!--       Top boundary condition for cloud water (Dirichlet)
          qc_p(nzt+1,:,:) = 0.0_wp
          nc_p(nzt+1,:,:) = 0.0_wp
           
       ENDIF

       IF ( cloud_physics  .AND.  microphysics_seifert )  THEN
!             
!--       Surface conditions rain water (Dirichlet)
!--       Run loop over all non-natural and natural walls. Note, in wall-datatype
!--       the k coordinate belongs to the atmospheric grid point, therefore, set
!--       qr_p and nr_p at k-1
          !$OMP PARALLEL DO PRIVATE( i, j, k )
          DO  m = 1, bc_h(0)%ns
             i = bc_h(0)%i(m)            
             j = bc_h(0)%j(m)
             k = bc_h(0)%k(m)
             qr_p(k-1,j,i) = 0.0_wp
             nr_p(k-1,j,i) = 0.0_wp
          ENDDO
!
!--       Top boundary condition for rain water (Dirichlet)
          qr_p(nzt+1,:,:) = 0.0_wp
          nr_p(nzt+1,:,:) = 0.0_wp
           
       ENDIF
    ENDIF
!
!-- Boundary conditions for scalar,
!-- bottom and top boundary (see also temperature)
    IF ( passive_scalar )  THEN
!
!--    Surface conditions for constant_humidity_flux
!--    Run loop over all non-natural and natural walls. Note, in wall-datatype
!--    the k coordinate belongs to the atmospheric grid point, therefore, set
!--    s_p at k-1
       IF ( ibc_s_b == 0 ) THEN
          
          DO  l = 0, 1
!
!--          Set kb, for upward-facing surfaces value at topography top (k-1) is set,
!--          for downward-facing surfaces at topography bottom (k+1). 
             kb = MERGE( -1, 1, l == 0 )
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, bc_h(l)%ns
                i = bc_h(l)%i(m)            
                j = bc_h(l)%j(m)
                k = bc_h(l)%k(m)
                s_p(k+kb,j,i) = s(k+kb,j,i)
             ENDDO
          ENDDO
          
       ELSE
          
          DO  l = 0, 1
!
!--          Set kb, for upward-facing surfaces value at topography top (k-1) is set,
!--          for downward-facing surfaces at topography bottom (k+1). 
             kb = MERGE( -1, 1, l == 0 )
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, bc_h(l)%ns
                i = bc_h(l)%i(m)            
                j = bc_h(l)%j(m)
                k = bc_h(l)%k(m)
                s_p(k+kb,j,i) = s_p(k,j,i)
             ENDDO
          ENDDO
       ENDIF
!
!--    Top boundary condition
       IF ( ibc_s_t == 0 )  THEN
          s_p(nzt+1,:,:) = s(nzt+1,:,:)
       ELSEIF ( ibc_s_t == 1 )  THEN
          s_p(nzt+1,:,:) = s_p(nzt,:,:)
       ELSEIF ( ibc_s_t == 2 )  THEN
          s_p(nzt+1,:,:) = s_p(nzt,:,:) + bc_s_t_val * dzu(nzt+1)
       ENDIF

    ENDIF    
!
!-- Top/bottom boundary conditions for chemical species
    IF ( air_chemistry )  CALL chem_boundary_conds( 'set_bc_bottomtop' )
!
!-- In case of inflow or nest boundary at the south boundary the boundary for v
!-- is at nys and in case of inflow or nest boundary at the left boundary the
!-- boundary for u is at nxl. Since in prognostic_equations (cache optimized
!-- version) these levels are handled as a prognostic level, boundary values
!-- have to be restored here.
!-- For the SGS-TKE, Neumann boundary conditions are used at the inflow.
    IF ( inflow_s )  THEN
       v_p(:,nys,:) = v_p(:,nys-1,:)
       IF ( .NOT. constant_diffusion ) e_p(:,nys-1,:) = e_p(:,nys,:)
    ELSEIF ( inflow_n )  THEN
       IF ( .NOT. constant_diffusion ) e_p(:,nyn+1,:) = e_p(:,nyn,:)
    ELSEIF ( inflow_l ) THEN
       u_p(:,:,nxl) = u_p(:,:,nxl-1)
       IF ( .NOT. constant_diffusion ) e_p(:,:,nxl-1) = e_p(:,:,nxl)
    ELSEIF ( inflow_r )  THEN
       IF ( .NOT. constant_diffusion ) e_p(:,:,nxr+1) = e_p(:,:,nxr)
    ENDIF

!
!-- The same restoration for u at i=nxl and v at j=nys as above must be made
!-- in case of nest boundaries. This must not be done in case of vertical nesting 
!-- mode as in that case the lateral boundaries are actually cyclic.
    IF ( nesting_mode /= 'vertical'  .OR.  forcing )  THEN
       IF ( nest_bound_s  .OR.  force_bound_s )  THEN
          v_p(:,nys,:) = v_p(:,nys-1,:)
       ENDIF
       IF ( nest_bound_l  .OR.  force_bound_l )  THEN
          u_p(:,:,nxl) = u_p(:,:,nxl-1)
       ENDIF
    ENDIF

!
!-- Lateral boundary conditions for scalar quantities at the outflow
    IF ( outflow_s )  THEN
       pt_p(:,nys-1,:)     = pt_p(:,nys,:)
       IF ( .NOT. constant_diffusion )  e_p(:,nys-1,:) = e_p(:,nys,:)
       IF ( rans_tke_e )  diss_p(:,nys-1,:) = diss_p(:,nys,:)
       IF ( humidity )  THEN
          q_p(:,nys-1,:) = q_p(:,nys,:)
          IF ( cloud_physics  .AND.  microphysics_morrison )  THEN
             qc_p(:,nys-1,:) = qc_p(:,nys,:)
             nc_p(:,nys-1,:) = nc_p(:,nys,:)
          ENDIF
          IF ( cloud_physics  .AND.  microphysics_seifert )  THEN
             qr_p(:,nys-1,:) = qr_p(:,nys,:)
             nr_p(:,nys-1,:) = nr_p(:,nys,:)
          ENDIF
       ENDIF
       IF ( passive_scalar )  s_p(:,nys-1,:) = s_p(:,nys,:)
    ELSEIF ( outflow_n )  THEN
       pt_p(:,nyn+1,:)     = pt_p(:,nyn,:)
       IF ( .NOT. constant_diffusion )  e_p(:,nyn+1,:) = e_p(:,nyn,:)
       IF ( rans_tke_e )  diss_p(:,nyn+1,:) = diss_p(:,nyn,:)
       IF ( humidity )  THEN
          q_p(:,nyn+1,:) = q_p(:,nyn,:)
          IF ( cloud_physics  .AND.  microphysics_morrison )  THEN
             qc_p(:,nyn+1,:) = qc_p(:,nyn,:)
             nc_p(:,nyn+1,:) = nc_p(:,nyn,:)
          ENDIF
          IF ( cloud_physics  .AND.  microphysics_seifert )  THEN
             qr_p(:,nyn+1,:) = qr_p(:,nyn,:)
             nr_p(:,nyn+1,:) = nr_p(:,nyn,:)
          ENDIF
       ENDIF
       IF ( passive_scalar )  s_p(:,nyn+1,:) = s_p(:,nyn,:)
    ELSEIF ( outflow_l )  THEN
       pt_p(:,:,nxl-1)     = pt_p(:,:,nxl)
       IF ( .NOT. constant_diffusion )  e_p(:,:,nxl-1) = e_p(:,:,nxl)
       IF ( rans_tke_e )  diss_p(:,:,nxl-1) = diss_p(:,:,nxl)
       IF ( humidity )  THEN
          q_p(:,:,nxl-1) = q_p(:,:,nxl)
          IF ( cloud_physics  .AND.  microphysics_morrison )  THEN
             qc_p(:,:,nxl-1) = qc_p(:,:,nxl)
             nc_p(:,:,nxl-1) = nc_p(:,:,nxl)
          ENDIF
          IF ( cloud_physics  .AND.  microphysics_seifert )  THEN
             qr_p(:,:,nxl-1) = qr_p(:,:,nxl)
             nr_p(:,:,nxl-1) = nr_p(:,:,nxl)
          ENDIF
       ENDIF
       IF ( passive_scalar )  s_p(:,:,nxl-1) = s_p(:,:,nxl)
    ELSEIF ( outflow_r )  THEN
       pt_p(:,:,nxr+1)     = pt_p(:,:,nxr)
       IF ( .NOT. constant_diffusion )  e_p(:,:,nxr+1) = e_p(:,:,nxr)
       IF ( rans_tke_e )  diss_p(:,:,nxr+1) = diss_p(:,:,nxr)
       IF ( humidity )  THEN
          q_p(:,:,nxr+1) = q_p(:,:,nxr)
          IF ( cloud_physics  .AND.  microphysics_morrison )  THEN
             qc_p(:,:,nxr+1) = qc_p(:,:,nxr)
             nc_p(:,:,nxr+1) = nc_p(:,:,nxr)
          ENDIF
          IF ( cloud_physics  .AND.  microphysics_seifert )  THEN
             qr_p(:,:,nxr+1) = qr_p(:,:,nxr)
             nr_p(:,:,nxr+1) = nr_p(:,:,nxr)
          ENDIF
       ENDIF
       IF ( passive_scalar )  s_p(:,:,nxr+1) = s_p(:,:,nxr)
    ENDIF

!
!-- Lateral boundary conditions for chemical species
    IF ( air_chemistry )  CALL chem_boundary_conds( 'set_bc_lateral' )    

!
!-- Radiation boundary conditions for the velocities at the respective outflow.
!-- The phase velocity is either assumed to the maximum phase velocity that
!-- ensures numerical stability (CFL-condition) or calculated after
!-- Orlanski(1976) and averaged along the outflow boundary.
    IF ( outflow_s )  THEN

       IF ( use_cmax )  THEN
          u_p(:,-1,:) = u(:,0,:)
          v_p(:,0,:)  = v(:,1,:)
          w_p(:,-1,:) = w(:,0,:)          
       ELSEIF ( .NOT. use_cmax )  THEN

          c_max = dy / dt_3d

          c_u_m_l = 0.0_wp 
          c_v_m_l = 0.0_wp
          c_w_m_l = 0.0_wp

          c_u_m = 0.0_wp 
          c_v_m = 0.0_wp
          c_w_m = 0.0_wp

!
!--       Calculate the phase speeds for u, v, and w, first local and then
!--       average along the outflow boundary.
          DO  k = nzb+1, nzt+1
             DO  i = nxl, nxr

                denom = u_m_s(k,0,i) - u_m_s(k,1,i)

                IF ( denom /= 0.0_wp )  THEN
                   c_u(k,i) = -c_max * ( u(k,0,i) - u_m_s(k,0,i) ) / ( denom * tsc(2) )
                   IF ( c_u(k,i) < 0.0_wp )  THEN
                      c_u(k,i) = 0.0_wp
                   ELSEIF ( c_u(k,i) > c_max )  THEN
                      c_u(k,i) = c_max
                   ENDIF
                ELSE
                   c_u(k,i) = c_max
                ENDIF

                denom = v_m_s(k,1,i) - v_m_s(k,2,i)

                IF ( denom /= 0.0_wp )  THEN
                   c_v(k,i) = -c_max * ( v(k,1,i) - v_m_s(k,1,i) ) / ( denom * tsc(2) )
                   IF ( c_v(k,i) < 0.0_wp )  THEN
                      c_v(k,i) = 0.0_wp
                   ELSEIF ( c_v(k,i) > c_max )  THEN
                      c_v(k,i) = c_max
                   ENDIF
                ELSE
                   c_v(k,i) = c_max
                ENDIF

                denom = w_m_s(k,0,i) - w_m_s(k,1,i)

                IF ( denom /= 0.0_wp )  THEN
                   c_w(k,i) = -c_max * ( w(k,0,i) - w_m_s(k,0,i) ) / ( denom * tsc(2) )
                   IF ( c_w(k,i) < 0.0_wp )  THEN
                      c_w(k,i) = 0.0_wp
                   ELSEIF ( c_w(k,i) > c_max )  THEN
                      c_w(k,i) = c_max
                   ENDIF
                ELSE
                   c_w(k,i) = c_max
                ENDIF

                c_u_m_l(k) = c_u_m_l(k) + c_u(k,i)
                c_v_m_l(k) = c_v_m_l(k) + c_v(k,i)
                c_w_m_l(k) = c_w_m_l(k) + c_w(k,i)

             ENDDO
          ENDDO

#if defined( __parallel )   
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dx, ierr )
          CALL MPI_ALLREDUCE( c_u_m_l(nzb+1), c_u_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dx, ierr )   
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dx, ierr )
          CALL MPI_ALLREDUCE( c_v_m_l(nzb+1), c_v_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dx, ierr )  
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dx, ierr )
          CALL MPI_ALLREDUCE( c_w_m_l(nzb+1), c_w_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dx, ierr )  
#else
          c_u_m = c_u_m_l
          c_v_m = c_v_m_l
          c_w_m = c_w_m_l
#endif

          c_u_m = c_u_m / (nx+1)
          c_v_m = c_v_m / (nx+1)
          c_w_m = c_w_m / (nx+1)

!
!--       Save old timelevels for the next timestep
          IF ( intermediate_timestep_count == 1 )  THEN
             u_m_s(:,:,:) = u(:,0:1,:)
             v_m_s(:,:,:) = v(:,1:2,:)
             w_m_s(:,:,:) = w(:,0:1,:)
          ENDIF

!
!--       Calculate the new velocities
          DO  k = nzb+1, nzt+1
             DO  i = nxlg, nxrg
                u_p(k,-1,i) = u(k,-1,i) - dt_3d * tsc(2) * c_u_m(k) *          &
                                       ( u(k,-1,i) - u(k,0,i) ) * ddy

                v_p(k,0,i)  = v(k,0,i)  - dt_3d * tsc(2) * c_v_m(k) *          &
                                       ( v(k,0,i) - v(k,1,i) ) * ddy

                w_p(k,-1,i) = w(k,-1,i) - dt_3d * tsc(2) * c_w_m(k) *          &
                                       ( w(k,-1,i) - w(k,0,i) ) * ddy
             ENDDO
          ENDDO

!
!--       Bottom boundary at the outflow
          IF ( ibc_uv_b == 0 )  THEN
             u_p(nzb,-1,:) = 0.0_wp 
             v_p(nzb,0,:)  = 0.0_wp  
          ELSE                    
             u_p(nzb,-1,:) =  u_p(nzb+1,-1,:)
             v_p(nzb,0,:)  =  v_p(nzb+1,0,:)
          ENDIF
          w_p(nzb,-1,:) = 0.0_wp

!
!--       Top boundary at the outflow
          IF ( ibc_uv_t == 0 )  THEN
             u_p(nzt+1,-1,:) = u_init(nzt+1)
             v_p(nzt+1,0,:)  = v_init(nzt+1)
          ELSE
             u_p(nzt+1,-1,:) = u_p(nzt,-1,:)
             v_p(nzt+1,0,:)  = v_p(nzt,0,:)
          ENDIF
          w_p(nzt:nzt+1,-1,:) = 0.0_wp

       ENDIF

    ENDIF

    IF ( outflow_n )  THEN

       IF ( use_cmax )  THEN
          u_p(:,ny+1,:) = u(:,ny,:)
          v_p(:,ny+1,:) = v(:,ny,:)
          w_p(:,ny+1,:) = w(:,ny,:)          
       ELSEIF ( .NOT. use_cmax )  THEN

          c_max = dy / dt_3d

          c_u_m_l = 0.0_wp 
          c_v_m_l = 0.0_wp
          c_w_m_l = 0.0_wp

          c_u_m = 0.0_wp 
          c_v_m = 0.0_wp
          c_w_m = 0.0_wp

!
!--       Calculate the phase speeds for u, v, and w, first local and then
!--       average along the outflow boundary.
          DO  k = nzb+1, nzt+1
             DO  i = nxl, nxr

                denom = u_m_n(k,ny,i) - u_m_n(k,ny-1,i)

                IF ( denom /= 0.0_wp )  THEN
                   c_u(k,i) = -c_max * ( u(k,ny,i) - u_m_n(k,ny,i) ) / ( denom * tsc(2) )
                   IF ( c_u(k,i) < 0.0_wp )  THEN
                      c_u(k,i) = 0.0_wp
                   ELSEIF ( c_u(k,i) > c_max )  THEN
                      c_u(k,i) = c_max
                   ENDIF
                ELSE
                   c_u(k,i) = c_max
                ENDIF

                denom = v_m_n(k,ny,i) - v_m_n(k,ny-1,i)

                IF ( denom /= 0.0_wp )  THEN
                   c_v(k,i) = -c_max * ( v(k,ny,i) - v_m_n(k,ny,i) ) / ( denom * tsc(2) )
                   IF ( c_v(k,i) < 0.0_wp )  THEN
                      c_v(k,i) = 0.0_wp
                   ELSEIF ( c_v(k,i) > c_max )  THEN
                      c_v(k,i) = c_max
                   ENDIF
                ELSE
                   c_v(k,i) = c_max
                ENDIF

                denom = w_m_n(k,ny,i) - w_m_n(k,ny-1,i)

                IF ( denom /= 0.0_wp )  THEN
                   c_w(k,i) = -c_max * ( w(k,ny,i) - w_m_n(k,ny,i) ) / ( denom * tsc(2) )
                   IF ( c_w(k,i) < 0.0_wp )  THEN
                      c_w(k,i) = 0.0_wp
                   ELSEIF ( c_w(k,i) > c_max )  THEN
                      c_w(k,i) = c_max
                   ENDIF
                ELSE
                   c_w(k,i) = c_max
                ENDIF

                c_u_m_l(k) = c_u_m_l(k) + c_u(k,i)
                c_v_m_l(k) = c_v_m_l(k) + c_v(k,i)
                c_w_m_l(k) = c_w_m_l(k) + c_w(k,i)

             ENDDO
          ENDDO

#if defined( __parallel )   
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dx, ierr )
          CALL MPI_ALLREDUCE( c_u_m_l(nzb+1), c_u_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dx, ierr )   
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dx, ierr )
          CALL MPI_ALLREDUCE( c_v_m_l(nzb+1), c_v_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dx, ierr )  
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dx, ierr )
          CALL MPI_ALLREDUCE( c_w_m_l(nzb+1), c_w_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dx, ierr )  
#else
          c_u_m = c_u_m_l
          c_v_m = c_v_m_l
          c_w_m = c_w_m_l
#endif

          c_u_m = c_u_m / (nx+1)
          c_v_m = c_v_m / (nx+1)
          c_w_m = c_w_m / (nx+1)

!
!--       Save old timelevels for the next timestep
          IF ( intermediate_timestep_count == 1 )  THEN
                u_m_n(:,:,:) = u(:,ny-1:ny,:)
                v_m_n(:,:,:) = v(:,ny-1:ny,:)
                w_m_n(:,:,:) = w(:,ny-1:ny,:)
          ENDIF

!
!--       Calculate the new velocities
          DO  k = nzb+1, nzt+1
             DO  i = nxlg, nxrg
                u_p(k,ny+1,i) = u(k,ny+1,i) - dt_3d * tsc(2) * c_u_m(k) *      &
                                       ( u(k,ny+1,i) - u(k,ny,i) ) * ddy

                v_p(k,ny+1,i) = v(k,ny+1,i)  - dt_3d * tsc(2) * c_v_m(k) *     &
                                       ( v(k,ny+1,i) - v(k,ny,i) ) * ddy

                w_p(k,ny+1,i) = w(k,ny+1,i) - dt_3d * tsc(2) * c_w_m(k) *      &
                                       ( w(k,ny+1,i) - w(k,ny,i) ) * ddy
             ENDDO
          ENDDO

!
!--       Bottom boundary at the outflow
          IF ( ibc_uv_b == 0 )  THEN
             u_p(nzb,ny+1,:) = 0.0_wp
             v_p(nzb,ny+1,:) = 0.0_wp   
          ELSE                    
             u_p(nzb,ny+1,:) =  u_p(nzb+1,ny+1,:)
             v_p(nzb,ny+1,:) =  v_p(nzb+1,ny+1,:)
          ENDIF
          w_p(nzb,ny+1,:) = 0.0_wp

!
!--       Top boundary at the outflow
          IF ( ibc_uv_t == 0 )  THEN
             u_p(nzt+1,ny+1,:) = u_init(nzt+1)
             v_p(nzt+1,ny+1,:) = v_init(nzt+1)
          ELSE
             u_p(nzt+1,ny+1,:) = u_p(nzt,nyn+1,:)
             v_p(nzt+1,ny+1,:) = v_p(nzt,nyn+1,:)
          ENDIF
          w_p(nzt:nzt+1,ny+1,:) = 0.0_wp

       ENDIF

    ENDIF

    IF ( outflow_l )  THEN

       IF ( use_cmax )  THEN
          u_p(:,:,0)  = u(:,:,1)
          v_p(:,:,-1) = v(:,:,0)
          w_p(:,:,-1) = w(:,:,0)          
       ELSEIF ( .NOT. use_cmax )  THEN

          c_max = dx / dt_3d

          c_u_m_l = 0.0_wp 
          c_v_m_l = 0.0_wp
          c_w_m_l = 0.0_wp

          c_u_m = 0.0_wp 
          c_v_m = 0.0_wp
          c_w_m = 0.0_wp

!
!--       Calculate the phase speeds for u, v, and w, first local and then
!--       average along the outflow boundary.
          DO  k = nzb+1, nzt+1
             DO  j = nys, nyn

                denom = u_m_l(k,j,1) - u_m_l(k,j,2)

                IF ( denom /= 0.0_wp )  THEN
                   c_u(k,j) = -c_max * ( u(k,j,1) - u_m_l(k,j,1) ) / ( denom * tsc(2) )
                   IF ( c_u(k,j) < 0.0_wp )  THEN
                      c_u(k,j) = 0.0_wp
                   ELSEIF ( c_u(k,j) > c_max )  THEN
                      c_u(k,j) = c_max
                   ENDIF
                ELSE
                   c_u(k,j) = c_max
                ENDIF

                denom = v_m_l(k,j,0) - v_m_l(k,j,1)

                IF ( denom /= 0.0_wp )  THEN
                   c_v(k,j) = -c_max * ( v(k,j,0) - v_m_l(k,j,0) ) / ( denom * tsc(2) )
                   IF ( c_v(k,j) < 0.0_wp )  THEN
                      c_v(k,j) = 0.0_wp
                   ELSEIF ( c_v(k,j) > c_max )  THEN
                      c_v(k,j) = c_max
                   ENDIF
                ELSE
                   c_v(k,j) = c_max
                ENDIF

                denom = w_m_l(k,j,0) - w_m_l(k,j,1)

                IF ( denom /= 0.0_wp )  THEN
                   c_w(k,j) = -c_max * ( w(k,j,0) - w_m_l(k,j,0) ) / ( denom * tsc(2) )
                   IF ( c_w(k,j) < 0.0_wp )  THEN
                      c_w(k,j) = 0.0_wp
                   ELSEIF ( c_w(k,j) > c_max )  THEN
                      c_w(k,j) = c_max
                   ENDIF
                ELSE
                   c_w(k,j) = c_max
                ENDIF

                c_u_m_l(k) = c_u_m_l(k) + c_u(k,j)
                c_v_m_l(k) = c_v_m_l(k) + c_v(k,j)
                c_w_m_l(k) = c_w_m_l(k) + c_w(k,j)

             ENDDO
          ENDDO

#if defined( __parallel )   
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dy, ierr )
          CALL MPI_ALLREDUCE( c_u_m_l(nzb+1), c_u_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dy, ierr )   
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dy, ierr )
          CALL MPI_ALLREDUCE( c_v_m_l(nzb+1), c_v_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dy, ierr )  
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dy, ierr )
          CALL MPI_ALLREDUCE( c_w_m_l(nzb+1), c_w_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dy, ierr )  
#else
          c_u_m = c_u_m_l
          c_v_m = c_v_m_l
          c_w_m = c_w_m_l
#endif

          c_u_m = c_u_m / (ny+1)
          c_v_m = c_v_m / (ny+1)
          c_w_m = c_w_m / (ny+1)

!
!--       Save old timelevels for the next timestep
          IF ( intermediate_timestep_count == 1 )  THEN
                u_m_l(:,:,:) = u(:,:,1:2)
                v_m_l(:,:,:) = v(:,:,0:1)
                w_m_l(:,:,:) = w(:,:,0:1)
          ENDIF

!
!--       Calculate the new velocities
          DO  k = nzb+1, nzt+1
             DO  j = nysg, nyng
                u_p(k,j,0) = u(k,j,0) - dt_3d * tsc(2) * c_u_m(k) *            &
                                       ( u(k,j,0) - u(k,j,1) ) * ddx

                v_p(k,j,-1) = v(k,j,-1) - dt_3d * tsc(2) * c_v_m(k) *          &
                                       ( v(k,j,-1) - v(k,j,0) ) * ddx

                w_p(k,j,-1) = w(k,j,-1) - dt_3d * tsc(2) * c_w_m(k) *          &
                                       ( w(k,j,-1) - w(k,j,0) ) * ddx
             ENDDO
          ENDDO

!
!--       Bottom boundary at the outflow
          IF ( ibc_uv_b == 0 )  THEN
             u_p(nzb,:,0)  = 0.0_wp 
             v_p(nzb,:,-1) = 0.0_wp
          ELSE                    
             u_p(nzb,:,0)  =  u_p(nzb+1,:,0)
             v_p(nzb,:,-1) =  v_p(nzb+1,:,-1)
          ENDIF
          w_p(nzb,:,-1) = 0.0_wp

!
!--       Top boundary at the outflow
          IF ( ibc_uv_t == 0 )  THEN
             u_p(nzt+1,:,0)  = u_init(nzt+1)
             v_p(nzt+1,:,-1) = v_init(nzt+1)
          ELSE
             u_p(nzt+1,:,0)  = u_p(nzt,:,0)
             v_p(nzt+1,:,-1) = v_p(nzt,:,-1)
          ENDIF
          w_p(nzt:nzt+1,:,-1) = 0.0_wp

       ENDIF

    ENDIF

    IF ( outflow_r )  THEN

       IF ( use_cmax )  THEN
          u_p(:,:,nx+1) = u(:,:,nx)
          v_p(:,:,nx+1) = v(:,:,nx)
          w_p(:,:,nx+1) = w(:,:,nx)          
       ELSEIF ( .NOT. use_cmax )  THEN

          c_max = dx / dt_3d

          c_u_m_l = 0.0_wp 
          c_v_m_l = 0.0_wp
          c_w_m_l = 0.0_wp

          c_u_m = 0.0_wp 
          c_v_m = 0.0_wp
          c_w_m = 0.0_wp

!
!--       Calculate the phase speeds for u, v, and w, first local and then
!--       average along the outflow boundary.
          DO  k = nzb+1, nzt+1
             DO  j = nys, nyn

                denom = u_m_r(k,j,nx) - u_m_r(k,j,nx-1)

                IF ( denom /= 0.0_wp )  THEN
                   c_u(k,j) = -c_max * ( u(k,j,nx) - u_m_r(k,j,nx) ) / ( denom * tsc(2) )
                   IF ( c_u(k,j) < 0.0_wp )  THEN
                      c_u(k,j) = 0.0_wp
                   ELSEIF ( c_u(k,j) > c_max )  THEN
                      c_u(k,j) = c_max
                   ENDIF
                ELSE
                   c_u(k,j) = c_max
                ENDIF

                denom = v_m_r(k,j,nx) - v_m_r(k,j,nx-1)

                IF ( denom /= 0.0_wp )  THEN
                   c_v(k,j) = -c_max * ( v(k,j,nx) - v_m_r(k,j,nx) ) / ( denom * tsc(2) )
                   IF ( c_v(k,j) < 0.0_wp )  THEN
                      c_v(k,j) = 0.0_wp
                   ELSEIF ( c_v(k,j) > c_max )  THEN
                      c_v(k,j) = c_max
                   ENDIF
                ELSE
                   c_v(k,j) = c_max
                ENDIF

                denom = w_m_r(k,j,nx) - w_m_r(k,j,nx-1)

                IF ( denom /= 0.0_wp )  THEN
                   c_w(k,j) = -c_max * ( w(k,j,nx) - w_m_r(k,j,nx) ) / ( denom * tsc(2) )
                   IF ( c_w(k,j) < 0.0_wp )  THEN
                      c_w(k,j) = 0.0_wp
                   ELSEIF ( c_w(k,j) > c_max )  THEN
                      c_w(k,j) = c_max
                   ENDIF
                ELSE
                   c_w(k,j) = c_max
                ENDIF

                c_u_m_l(k) = c_u_m_l(k) + c_u(k,j)
                c_v_m_l(k) = c_v_m_l(k) + c_v(k,j)
                c_w_m_l(k) = c_w_m_l(k) + c_w(k,j)

             ENDDO
          ENDDO

#if defined( __parallel )   
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dy, ierr )
          CALL MPI_ALLREDUCE( c_u_m_l(nzb+1), c_u_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dy, ierr )   
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dy, ierr )
          CALL MPI_ALLREDUCE( c_v_m_l(nzb+1), c_v_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dy, ierr )  
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dy, ierr )
          CALL MPI_ALLREDUCE( c_w_m_l(nzb+1), c_w_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dy, ierr )  
#else
          c_u_m = c_u_m_l
          c_v_m = c_v_m_l
          c_w_m = c_w_m_l
#endif

          c_u_m = c_u_m / (ny+1)
          c_v_m = c_v_m / (ny+1)
          c_w_m = c_w_m / (ny+1)

!
!--       Save old timelevels for the next timestep
          IF ( intermediate_timestep_count == 1 )  THEN
                u_m_r(:,:,:) = u(:,:,nx-1:nx)
                v_m_r(:,:,:) = v(:,:,nx-1:nx)
                w_m_r(:,:,:) = w(:,:,nx-1:nx)
          ENDIF

!
!--       Calculate the new velocities
          DO  k = nzb+1, nzt+1
             DO  j = nysg, nyng
                u_p(k,j,nx+1) = u(k,j,nx+1) - dt_3d * tsc(2) * c_u_m(k) *      &
                                       ( u(k,j,nx+1) - u(k,j,nx) ) * ddx

                v_p(k,j,nx+1) = v(k,j,nx+1) - dt_3d * tsc(2) * c_v_m(k) *      &
                                       ( v(k,j,nx+1) - v(k,j,nx) ) * ddx

                w_p(k,j,nx+1) = w(k,j,nx+1) - dt_3d * tsc(2) * c_w_m(k) *      &
                                       ( w(k,j,nx+1) - w(k,j,nx) ) * ddx
             ENDDO
          ENDDO

!
!--       Bottom boundary at the outflow
          IF ( ibc_uv_b == 0 )  THEN
             u_p(nzb,:,nx+1) = 0.0_wp
             v_p(nzb,:,nx+1) = 0.0_wp 
          ELSE                    
             u_p(nzb,:,nx+1) =  u_p(nzb+1,:,nx+1)
             v_p(nzb,:,nx+1) =  v_p(nzb+1,:,nx+1)
          ENDIF
          w_p(nzb,:,nx+1) = 0.0_wp

!
!--       Top boundary at the outflow
          IF ( ibc_uv_t == 0 )  THEN
             u_p(nzt+1,:,nx+1) = u_init(nzt+1)
             v_p(nzt+1,:,nx+1) = v_init(nzt+1)
          ELSE
             u_p(nzt+1,:,nx+1) = u_p(nzt,:,nx+1)
             v_p(nzt+1,:,nx+1) = v_p(nzt,:,nx+1)
          ENDIF
          w_p(nzt:nzt+1,:,nx+1) = 0.0_wp

       ENDIF

    ENDIF

 END SUBROUTINE boundary_conds
