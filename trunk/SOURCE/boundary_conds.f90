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

          e_p(nzt+1,:,:) = e_p(nzt,:,:)
!
!--    Nesting case: if parent operates in RANS mode and child in LES mode,
!--    no TKE is transfered. This case, set Neumann conditions at lateral and 
!--    top child boundaries. 
!--    If not ( both either in RANS or in LES mode ), TKE boundary condition
!--    is treated in the nesting. 
   ENDIF

!
!-- Boundary conditions for TKE dissipation rate. 
    IF ( rans_tke_e .AND. .NOT. nest_domain )  THEN
       diss_p(nzt+1,:,:) = diss_p(nzt,:,:)
    ENDIF

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


 END SUBROUTINE boundary_conds
