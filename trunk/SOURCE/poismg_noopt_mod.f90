!> @file poismg_noopt_mod.f90
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
! $Id: poismg_noopt_mod.f90 2939 2018-03-29 18:20:00Z suehring $
! Set lateral boundary conditions for divergence
! 
! 2793 2018-02-07 10:54:33Z suehring
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! - Change in file header (GPL part)
! - Rename file into poismg_noopt_mod
! - Modularize poismg_noopt - initialization of masking flags moved from 
!   init_grid
! - Masking of downward-facing walls walls 
! - Large-scale forcing implemented (MS) 
! 
! 2591 2017-10-30 08:07:25Z maronga
! 
! 2586 2017-10-27 13:03:50Z kanani
!
! 2232 2017-05-30 17:47:52Z suehring
! Bugfixes OpenMP
! 
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
! 
! 2021 2016-10-07 14:08:57Z suehring
! Bugfix: restore nest_bound_(l/r/s/n) in case of mg_switch_to_pe0
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1934 2016-06-13 09:46:57Z hellstea
! Rename subroutines and cpu-measure log points to indicate _noopt version
! 
! 1762 2016-02-25 12:31:13Z hellstea
! Introduction of nested domain feature
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
! old module precision_kind is removed,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1318 2014-03-17 13:35:16Z raasch
! module interfaces removed
!
! 1159 2013-05-21 11:58:22Z fricke
! bc_lr/ns_dirneu/neudir removed
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1056 2012-11-16 15:28:04Z raasch
! Bugfix: all ghost points have to be used for allocating p3
! arrays p2, f2, and f2_l changed from allocatable to automatic
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 996 2012-09-07 10:41:47Z raasch
! little reformatting
!
! 978 2012-08-09 08:28:32Z fricke
! bc_lr/ns_dirneu/neudir added
!
! 880 2012-04-13 06:28:59Z raasch
! Bugfix: preprocessor statements for parallel execution added
!
! 778 2011-11-07 14:18:25Z fricke
! Allocation of p3 changes when multigrid is used and the collected field on PE0
! has more grid points than the subdomain of an PE.
!
! 707 2011-03-29 11:39:40Z raasch
! p_loc is used instead of p in the main routine (poismg).
! On coarse grid levels, gathered data are identically processed on all PEs
! (before, on PE0 only), so that the subsequent scattering of data is not
! neccessary any more.
! bc_lr/ns replaced by bc_lr/ns_cyc/dirrad/raddir
! Bugfix: bottom (nzb) and top (nzt+1) boundary conditions set in routines
! resid and restrict. They were missed before which may have led to
! unpredictable results.
!
! 667 2010-12-23 12:06:00Z suehring/gryschka
! Calls of exchange_horiz are modified.
!
! 622 2010-12-10 08:08:13Z raasch
! optional barriers included in order to speed up collective operations
!
! 257 2009-03-11 15:17:42Z heinze
! Output of messages replaced by message handling routine.
!
! 181 2008-07-30 07:07:47Z raasch
! Bugfix: grid_level+1 has to be used in restrict for flags-array
!
! 114 2007-10-10 00:03:15Z raasch
! Boundary conditions at walls are implicitly set using flag arrays. Only
! Neumann BC is allowed. Upper walls are still not realized.
! Bottom and top BCs for array f_mg in restrict removed because boundary
! values are not needed (right hand side of SOR iteration).
!
! 75 2007-03-22 09:54:05Z raasch
! 2nd+3rd argument removed from exchange horiz
!
! RCS Log replace by Id keyword, revision history cleaned up
!
! Revision 1.6  2005/03/26 20:55:54  raasch
! Implementation of non-cyclic (Neumann) horizontal boundary conditions,
! routine prolong simplified (one call of exchange_horiz spared)
!
! Revision 1.1  2001/07/20 13:10:51  raasch
! Initial revision
!
!
! Description:
! ------------
!> Solves the Poisson equation for the perturbation pressure with a multigrid
!> V- or W-Cycle scheme.
!>
!> This multigrid method was originally developed for PALM by Joerg Uhlenbrock,
!> September 2000 - July 2001.
!> 
!> @attention Loop unrolling and cache optimization in SOR-Red/Black method
!>            still does not give the expected speedup! 
!>
!> @todo Further work required.
!> @todo Formatting adjustments required (indention after modularization) 
!------------------------------------------------------------------------------!
 MODULE poismg_noopt_mod
 
    USE control_parameters,                                                    &
        ONLY:  grid_level, force_bound_l, force_bound_n, force_bound_r,        &
               force_bound_s, forcing, inflow_l, inflow_n, inflow_r, inflow_s, &
               outflow_l, outflow_n, outflow_r, outflow_s

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE kinds

    USE pegrid

    PRIVATE

    INTERFACE poismg_noopt
       MODULE PROCEDURE poismg_noopt
    END INTERFACE poismg_noopt

    INTERFACE poismg_noopt_init
       MODULE PROCEDURE poismg_noopt_init
    END INTERFACE poismg_noopt_init

    PUBLIC poismg_noopt, poismg_noopt_init

 CONTAINS

    SUBROUTINE poismg_noopt( r )
    

       USE arrays_3d,                                                          &
           ONLY:  d, p_loc

       USE control_parameters,                                                 &
           ONLY:  bc_lr_cyc, bc_ns_cyc, gathered_size, grid_level_count,       &
                  ibc_p_t, maximum_grid_level, message_string, mgcycles,       &
                  mg_cycles, mg_switch_to_pe0_level, residual_limit,           &
                  subdomain_size

       USE control_parameters,                                                 &
           ONLY:  bc_lr_cyc, bc_ns_cyc, gathered_size, grid_level,             &
                  grid_level_count, ibc_p_t,                                   &
                  maximum_grid_level, message_string, mgcycles, mg_cycles,     &
                  mg_switch_to_pe0_level, residual_limit, subdomain_size

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE indices,                                                            &
           ONLY:  nxl, nxlg, nxl_mg, nxr, nxrg, nxr_mg, nys, nysg, nys_mg, nyn,&
                  nyng, nyn_mg, nzb, nzt, nzt_mg

       USE kinds

       USE pegrid

       IMPLICIT NONE

       REAL(wp) ::  maxerror          !<
       REAL(wp) ::  maximum_mgcycles  !<
       REAL(wp) ::  residual_norm     !<

       REAL(wp), DIMENSION(nzb:nzt+1,nys-1:nyn+1,nxl-1:nxr+1) ::  r  !<

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  p3  !<


       CALL cpu_log( log_point_s(29), 'poismg_noopt', 'start' )
!
!--    Initialize arrays and variables used in this subroutine

!--    If the number of grid points of the gathered grid, which is collected
!--    on PE0, is larger than the number of grid points of an PE, than array
!--    p3 will be enlarged.
       IF ( gathered_size > subdomain_size )  THEN 
          ALLOCATE( p3(nzb:nzt_mg(mg_switch_to_pe0_level)+1,nys_mg(               &
                       mg_switch_to_pe0_level)-1:nyn_mg(mg_switch_to_pe0_level)+1,&
                       nxl_mg(mg_switch_to_pe0_level)-1:nxr_mg(                   &
                       mg_switch_to_pe0_level)+1) )
       ELSE
          ALLOCATE ( p3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ENDIF

       p3 = 0.0_wp
    
!
!--    Ghost boundaries have to be added to divergence array.
!--    Exchange routine needs to know the grid level!
       grid_level = maximum_grid_level
       CALL exchange_horiz( d, 1)
!
!--    Set bottom and top boundary conditions
       d(nzb,:,:) = d(nzb+1,:,:)
       IF ( ibc_p_t == 1 )  d(nzt+1,:,: ) = d(nzt,:,:)
!
!
!--    Initiation of the multigrid scheme. Does n cycles until the 
!--    residual is smaller than the given limit. The accuracy of the solution 
!--    of the poisson equation will increase with the number of cycles.
!--    If the number of cycles is preset by the user, this number will be
!--    carried out regardless of the accuracy.
       grid_level_count =  0
       mgcycles         =  0
       IF ( mg_cycles == -1 )  THEN
          maximum_mgcycles = 0
          residual_norm    = 1.0_wp
       ELSE
          maximum_mgcycles = mg_cycles
          residual_norm    = 0.0_wp
       ENDIF

       DO WHILE ( residual_norm > residual_limit  .OR. &
                  mgcycles < maximum_mgcycles )
    
          CALL next_mg_level_noopt( d, p_loc, p3, r)

!
!--       Calculate the residual if the user has not preset the number of
!--       cycles to be performed
          IF ( maximum_mgcycles == 0 )  THEN
             CALL resid_noopt( d, p_loc, r )
             maxerror = SUM( r(nzb+1:nzt,nys:nyn,nxl:nxr)**2 )

#if defined( __parallel )
             IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
                CALL MPI_ALLREDUCE( maxerror, residual_norm, 1, MPI_REAL,      &
                                    MPI_SUM, comm2d, ierr)
#else
                residual_norm = maxerror
#endif
             residual_norm = SQRT( residual_norm )
          ENDIF

          mgcycles = mgcycles + 1

!
!--       If the user has not limited the number of cycles, stop the run in case
!--       of insufficient convergence
          IF ( mgcycles > 1000  .AND.  mg_cycles == -1 )  THEN
             message_string = 'no sufficient convergence within 1000 cycles'
             CALL message( 'poismg_noopt', 'PA0283', 1, 2, 0, 6, 0 )
          ENDIF

       ENDDO

       DEALLOCATE( p3 )

!
!--    Unset the grid level. Variable is used to determine the MPI datatypes for
!--    ghost point exchange
       grid_level = 0 

       CALL cpu_log( log_point_s(29), 'poismg_noopt', 'stop' )

    END SUBROUTINE poismg_noopt


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the residual of the perturbation pressure.
!------------------------------------------------------------------------------!
    SUBROUTINE resid_noopt( f_mg, p_mg, r )


       USE arrays_3d,                                                          &
           ONLY:  f1_mg, f2_mg, f3_mg, rho_air_mg

       USE control_parameters,                                                 &
           ONLY:  bc_lr_cyc, bc_ns_cyc, ibc_p_b, ibc_p_t

       USE grid_variables,                                                     &
           ONLY:  ddx2_mg, ddy2_mg

       USE indices,                                                            &
           ONLY:  flags, wall_flags_1, wall_flags_2, wall_flags_3, wall_flags_4,&
                  wall_flags_5, wall_flags_6, wall_flags_7, wall_flags_8,      &
                  wall_flags_9, wall_flags_10, nxl_mg, nxr_mg, nys_mg, nyn_mg, &
                  nzb, nzt_mg

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i
       INTEGER(iwp) ::  j
       INTEGER(iwp) ::  k
       INTEGER(iwp) ::  l

       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                       nys_mg(grid_level)-1:nyn_mg(grid_level)+1,              &
                       nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  f_mg  !<
       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                       nys_mg(grid_level)-1:nyn_mg(grid_level)+1,              &
                       nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  p_mg  !<
       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                       nys_mg(grid_level)-1:nyn_mg(grid_level)+1,              &
                       nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  r     !<

!
!--    Calculate the residual
       l = grid_level

!
!--    Choose flag array of this level
       SELECT CASE ( l )
          CASE ( 1 )
             flags => wall_flags_1
          CASE ( 2 )
             flags => wall_flags_2
          CASE ( 3 )
             flags => wall_flags_3
          CASE ( 4 )
             flags => wall_flags_4
          CASE ( 5 )
             flags => wall_flags_5
          CASE ( 6 )
             flags => wall_flags_6
          CASE ( 7 )
             flags => wall_flags_7
          CASE ( 8 )
             flags => wall_flags_8
          CASE ( 9 )
             flags => wall_flags_9
          CASE ( 10 )
             flags => wall_flags_10
       END SELECT

!$OMP PARALLEL PRIVATE (i,j,k)
!$OMP DO
       DO  i = nxl_mg(l), nxr_mg(l)
          DO  j = nys_mg(l), nyn_mg(l) 
             DO  k = nzb+1, nzt_mg(l)
                r(k,j,i) = f_mg(k,j,i)                                         &
                           - rho_air_mg(k,l) * ddx2_mg(l) *                    &
                               ( p_mg(k,j,i+1) + IBITS( flags(k,j,i), 5, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i+1) ) + &
                                 p_mg(k,j,i-1) + IBITS( flags(k,j,i), 4, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i-1) ) ) &
                           - rho_air_mg(k,l) * ddy2_mg(l) *                    &
                               ( p_mg(k,j+1,i) + IBITS( flags(k,j,i), 3, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j+1,i) ) + &
                                 p_mg(k,j-1,i) + IBITS( flags(k,j,i), 2, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j-1,i) ) ) &
                           - f2_mg(k,l) *                                      &
                               ( p_mg(k+1,j,i) + IBITS( flags(k,j,i), 0, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k+1,j,i) ) ) &
                           - f3_mg(k,l) *                                      &
                               ( p_mg(k-1,j,i) + IBITS( flags(k,j,i), 0, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k-1,j,i) ) ) &
                           + f1_mg(k,l) * p_mg(k,j,i)
!
!--             Residual within topography should be zero
                r(k,j,i) = r(k,j,i) * ( 1.0_wp - IBITS( flags(k,j,i), 6, 1 ) )
             ENDDO
          ENDDO
       ENDDO
!$OMP END PARALLEL

!
!--    Horizontal boundary conditions
       CALL exchange_horiz( r, 1)

!--    Boundary conditions at bottom and top of the domain.
!--    These points are not handled by the above loop. Points may be within
!--    buildings, but that doesn't matter. 
       IF ( ibc_p_b == 1 )  THEN
          r(nzb,:,: ) = r(nzb+1,:,:)
       ELSE
          r(nzb,:,: ) = 0.0_wp
       ENDIF

       IF ( ibc_p_t == 1 )  THEN
          r(nzt_mg(l)+1,:,: ) = r(nzt_mg(l),:,:)
       ELSE
          r(nzt_mg(l)+1,:,: ) = 0.0_wp
       ENDIF


    END SUBROUTINE resid_noopt


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolates the residual on the next coarser grid with "full weighting"
!> scheme.
!------------------------------------------------------------------------------!
    SUBROUTINE restrict_noopt( f_mg, r )


       USE control_parameters,                                                 &
           ONLY:  bc_lr_cyc, bc_ns_cyc, grid_level, ibc_p_b, ibc_p_t

       USE indices,                                                            &
           ONLY:  flags, wall_flags_1, wall_flags_2, wall_flags_3, wall_flags_4,&
                  wall_flags_5, wall_flags_6, wall_flags_7, wall_flags_8,      &
                  wall_flags_9, wall_flags_10, nxl_mg, nxr_mg, nys_mg, nyn_mg, &
                  nzb, nzt_mg

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i    !<
       INTEGER(iwp) ::  ic   !<
       INTEGER(iwp) ::  j    !<
       INTEGER(iwp) ::  jc   !<
       INTEGER(iwp) ::  k    !<
       INTEGER(iwp) ::  kc   !<
       INTEGER(iwp) ::  l    !<

       REAL(wp) ::  rkjim    !<
       REAL(wp) ::  rkjip    !<
       REAL(wp) ::  rkjmi    !<
       REAL(wp) ::  rkjmim   !<
       REAL(wp) ::  rkjmip   !<
       REAL(wp) ::  rkjpi    !<
       REAL(wp) ::  rkjpim   !<
       REAL(wp) ::  rkjpip   !<
       REAL(wp) ::  rkmji    !<
       REAL(wp) ::  rkmjim   !<
       REAL(wp) ::  rkmjip   !<
       REAL(wp) ::  rkmjmi   !<
       REAL(wp) ::  rkmjmim  !<
       REAL(wp) ::  rkmjmip  !<
       REAL(wp) ::  rkmjpi   !<
       REAL(wp) ::  rkmjpim  !<
       REAL(wp) ::  rkmjpip  !<

       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                           nys_mg(grid_level)-1:nyn_mg(grid_level)+1,          &
                           nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  f_mg  !<

       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level+1)+1,                         &
                           nys_mg(grid_level+1)-1:nyn_mg(grid_level+1)+1,      &
                           nxl_mg(grid_level+1)-1:nxr_mg(grid_level+1)+1) ::  r !<

!
!--    Interpolate the residual
       l = grid_level

!
!--    Choose flag array of the upper level
       SELECT CASE ( l+1 )
          CASE ( 1 )
             flags => wall_flags_1
          CASE ( 2 )
             flags => wall_flags_2
          CASE ( 3 )
             flags => wall_flags_3
          CASE ( 4 )
             flags => wall_flags_4
          CASE ( 5 )
             flags => wall_flags_5
          CASE ( 6 )
             flags => wall_flags_6
          CASE ( 7 )
             flags => wall_flags_7
          CASE ( 8 )
             flags => wall_flags_8
          CASE ( 9 )
             flags => wall_flags_9
          CASE ( 10 )
             flags => wall_flags_10
       END SELECT
       
!$OMP PARALLEL PRIVATE (i,j,k,ic,jc,kc, rkjim,rkjip,rkjpi,rkjmi,rkjmim,rkjpim, &
!$OMP rkjmip, rkjpip,rkmji,rkmjim,rkmjip,rkmjpi,rkmjmi,rkmjmim,rkmjpim,rkmjmip,&
!$OMP rkmjpip          )
!$OMP DO
       DO  ic = nxl_mg(l), nxr_mg(l)    
          i = 2*ic 
          DO  jc = nys_mg(l), nyn_mg(l)  
             j = 2*jc
             DO  kc = nzb+1, nzt_mg(l)
                k = 2*kc-1
!
!--             Use implicit Neumann BCs if the respective gridpoint is inside
!--             the building
                rkjim   = r(k,j,i-1)     + IBITS( flags(k,j,i-1), 6, 1 ) *     &
                                           ( r(k,j,i) - r(k,j,i-1) )
                rkjip   = r(k,j,i+1)     + IBITS( flags(k,j,i+1), 6, 1 ) *     &
                                           ( r(k,j,i) - r(k,j,i+1) )
                rkjpi   = r(k,j+1,i)     + IBITS( flags(k,j+1,i), 6, 1 ) *     &
                                           ( r(k,j,i) - r(k,j+1,i) )
                rkjmi   = r(k,j-1,i)     + IBITS( flags(k,j-1,i), 6, 1 ) *     &
                                           ( r(k,j,i) - r(k,j-1,i) )
                rkjmim  = r(k,j-1,i-1)   + IBITS( flags(k,j-1,i-1), 6, 1 ) *   &
                                           ( r(k,j,i) - r(k,j-1,i-1) )
                rkjpim  = r(k,j+1,i-1)   + IBITS( flags(k,j+1,i-1), 6, 1 ) *   &
                                           ( r(k,j,i) - r(k,j+1,i-1) )
                rkjmip  = r(k,j-1,i+1)   + IBITS( flags(k,j-1,i+1), 6, 1 ) *   &
                                           ( r(k,j,i) - r(k,j-1,i+1) )
                rkjpip  = r(k,j+1,i+1)   + IBITS( flags(k,j+1,i+1), 6, 1 ) *   &
                                           ( r(k,j,i) - r(k,j+1,i+1) )
                rkmji   = r(k-1,j,i)     + IBITS( flags(k-1,j,i), 6, 1 ) *     &
                                           ( r(k,j,i) - r(k-1,j,i) )
                rkmjim  = r(k-1,j,i-1)   + IBITS( flags(k-1,j,i-1), 6, 1 ) *   &
                                           ( r(k,j,i) - r(k-1,j,i-1) )
                rkmjip  = r(k-1,j,i+1)   + IBITS( flags(k-1,j,i+1), 6, 1 ) *   &
                                           ( r(k,j,i) - r(k-1,j,i+1) )
                rkmjpi  = r(k-1,j+1,i)   + IBITS( flags(k-1,j+1,i), 6, 1 ) *   &
                                           ( r(k,j,i) - r(k-1,j+1,i) )
                rkmjmi  = r(k-1,j-1,i)   + IBITS( flags(k-1,j-1,i), 6, 1 ) *   &
                                           ( r(k,j,i) - r(k-1,j-1,i) )
                rkmjmim = r(k-1,j-1,i-1) + IBITS( flags(k-1,j-1,i-1), 6, 1 ) * &
                                           ( r(k,j,i) - r(k-1,j-1,i-1) )
                rkmjpim = r(k-1,j+1,i-1) + IBITS( flags(k-1,j+1,i-1), 6, 1 ) * &
                                           ( r(k,j,i) - r(k-1,j+1,i-1) )
                rkmjmip = r(k-1,j-1,i+1) + IBITS( flags(k-1,j-1,i+1), 6, 1 ) * &
                                           ( r(k,j,i) - r(k-1,j-1,i+1) )
                rkmjpip = r(k-1,j+1,i+1) + IBITS( flags(k-1,j+1,i+1), 6, 1 ) * &
                                           ( r(k,j,i) - r(k-1,j+1,i+1) )

                f_mg(kc,jc,ic) = 1.0_wp / 64.0_wp * (                         &
                                 8.0_wp * r(k,j,i)                            &
                               + 4.0_wp * ( rkjim   + rkjip   +               &
                                            rkjpi   + rkjmi   )               &
                               + 2.0_wp * ( rkjmim  + rkjpim  +               &
                                            rkjmip  + rkjpip  )               &
                               + 4.0_wp * rkmji                               &
                               + 2.0_wp * ( rkmjim  + rkmjim  +               &
                                            rkmjpi  + rkmjmi  )               &
                               +          ( rkmjmim + rkmjpim +               &
                                            rkmjmip + rkmjpip )               &
                               + 4.0_wp * r(k+1,j,i)                          &
                               + 2.0_wp * ( r(k+1,j,i-1)   + r(k+1,j,i+1)   + &
                                            r(k+1,j+1,i)   + r(k+1,j-1,i)   ) &
                               +          ( r(k+1,j-1,i-1) + r(k+1,j+1,i-1) + &
                                            r(k+1,j-1,i+1) + r(k+1,j+1,i+1) ) &
                                                    )

!             f_mg(kc,jc,ic) = 1.0_wp / 64.0_wp * (                         &
!                              8.0_wp * r(k,j,i)                            &
!                            + 4.0_wp * ( r(k,j,i-1)     + r(k,j,i+1)     + &
!                                         r(k,j+1,i)     + r(k,j-1,i)     ) &
!                            + 2.0_wp * ( r(k,j-1,i-1)   + r(k,j+1,i-1)   + &
!                                         r(k,j-1,i+1)   + r(k,j+1,i+1)   ) &
!                            + 4.0_wp * r(k-1,j,i)                          &
!                            + 2.0_wp * ( r(k-1,j,i-1)   + r(k-1,j,i+1)   + &
!                                         r(k-1,j+1,i)   + r(k-1,j-1,i)   ) &
!                            +          ( r(k-1,j-1,i-1) + r(k-1,j+1,i-1) + &
!                                         r(k-1,j-1,i+1) + r(k-1,j+1,i+1) ) &
!                            + 4.0_wp * r(k+1,j,i)                          &
!                            + 2.0_wp * ( r(k+1,j,i-1)   + r(k+1,j,i+1)   + &
!                                         r(k+1,j+1,i)   + r(k+1,j-1,i)   ) &
!                            +          ( r(k+1,j-1,i-1) + r(k+1,j+1,i-1) + &
!                                         r(k+1,j-1,i+1) + r(k+1,j+1,i+1) ) &
!                                                )
             ENDDO
          ENDDO
       ENDDO
!$OMP END PARALLEL

!
!--    Horizontal boundary conditions
       CALL exchange_horiz( f_mg, 1)

!--    Boundary conditions at bottom and top of the domain.
!--    These points are not handled by the above loop. Points may be within
!--    buildings, but that doesn't matter. 
       IF ( ibc_p_b == 1 )  THEN
          f_mg(nzb,:,: ) = f_mg(nzb+1,:,:)
       ELSE
          f_mg(nzb,:,: ) = 0.0_wp
       ENDIF

       IF ( ibc_p_t == 1 )  THEN
          f_mg(nzt_mg(l)+1,:,: ) = f_mg(nzt_mg(l),:,:)
       ELSE
          f_mg(nzt_mg(l)+1,:,: ) = 0.0_wp
       ENDIF


   END SUBROUTINE restrict_noopt


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolates the correction of the perturbation pressure 
!> to the next finer grid.
!------------------------------------------------------------------------------!
 SUBROUTINE prolong_noopt( p, temp )


    USE control_parameters,                                                    &
        ONLY:  bc_lr_cyc, bc_ns_cyc, ibc_p_b, ibc_p_t
    USE indices,                                                               &
        ONLY:  nxl_mg, nxr_mg, nys_mg, nyn_mg, nzb, nzt_mg

    USE kinds

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  l  !<

    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level-1)+1,                           &
                        nys_mg(grid_level-1)-1:nyn_mg(grid_level-1)+1,        &
                        nxl_mg(grid_level-1)-1:nxr_mg(grid_level-1)+1 ) ::  p  !<

    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                        nys_mg(grid_level)-1:nyn_mg(grid_level)+1,          &
                        nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  temp  !<


!
!-- First, store elements of the coarser grid on the next finer grid
    l = grid_level

!$OMP PARALLEL PRIVATE (i,j,k)
!$OMP DO
    DO  i = nxl_mg(l-1), nxr_mg(l-1)
       DO  j = nys_mg(l-1), nyn_mg(l-1)
!CDIR NODEP
          DO  k = nzb+1, nzt_mg(l-1)
!
!--          Points of the coarse grid are directly stored on the next finer
!--          grid
             temp(2*k-1,2*j,2*i) = p(k,j,i) 
! 
!--          Points between two coarse-grid points
             temp(2*k-1,2*j,2*i+1) = 0.5_wp * ( p(k,j,i) + p(k,j,i+1) )
             temp(2*k-1,2*j+1,2*i) = 0.5_wp * ( p(k,j,i) + p(k,j+1,i) )
             temp(2*k,2*j,2*i)     = 0.5_wp * ( p(k,j,i) + p(k+1,j,i) )
!
!--          Points in the center of the planes stretched by four points
!--          of the coarse grid cube
             temp(2*k-1,2*j+1,2*i+1) = 0.25_wp * ( p(k,j,i)   + p(k,j,i+1) + &
                                                   p(k,j+1,i) + p(k,j+1,i+1) )
             temp(2*k,2*j,2*i+1)     = 0.25_wp * ( p(k,j,i)   + p(k,j,i+1) + &
                                                   p(k+1,j,i) + p(k+1,j,i+1) )
             temp(2*k,2*j+1,2*i)     = 0.25_wp * ( p(k,j,i)   + p(k,j+1,i) + &
                                                   p(k+1,j,i) + p(k+1,j+1,i) )
!
!--          Points in the middle of coarse grid cube
             temp(2*k,2*j+1,2*i+1) = 0.125_wp * ( p(k,j,i)     + p(k,j,i+1)   + &
                                                  p(k,j+1,i)   + p(k,j+1,i+1) + &
                                                  p(k+1,j,i)   + p(k+1,j,i+1) + &
                                                  p(k+1,j+1,i) + p(k+1,j+1,i+1) )
          ENDDO
       ENDDO
    ENDDO
!$OMP END PARALLEL 
                          
!
!-- Horizontal boundary conditions
    CALL exchange_horiz( temp, 1)

!-- Bottom and top boundary conditions
    IF ( ibc_p_b == 1 )  THEN
       temp(nzb,:,: ) = temp(nzb+1,:,:)
    ELSE
       temp(nzb,:,: ) = 0.0_wp
    ENDIF

    IF ( ibc_p_t == 1 )  THEN
       temp(nzt_mg(l)+1,:,: ) = temp(nzt_mg(l),:,:)
    ELSE
       temp(nzt_mg(l)+1,:,: ) = 0.0_wp
    ENDIF

  
 END SUBROUTINE prolong_noopt


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Relaxation method for the multigrid scheme. A Gauss-Seidel iteration with
!> 3D-Red-Black decomposition (GS-RB) is used.
!------------------------------------------------------------------------------!
    SUBROUTINE redblack_noopt( f_mg, p_mg )


       USE arrays_3d,                                                          &
           ONLY:  f1_mg, f2_mg, f3_mg, rho_air_mg

       USE control_parameters,                                                 &
           ONLY:  bc_lr_cyc, bc_ns_cyc, ibc_p_b, ibc_p_t, ngsrb

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE grid_variables,                                                     &
           ONLY:  ddx2_mg, ddy2_mg

       USE indices,                                                            &
           ONLY:  flags, wall_flags_1, wall_flags_2, wall_flags_3, wall_flags_4,&
                  wall_flags_5, wall_flags_6, wall_flags_7, wall_flags_8,      &
                  wall_flags_9, wall_flags_10, nxl_mg, nxr_mg, nys_mg, nyn_mg, &
                  nzb, nzt_mg

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) :: color    !<
       INTEGER(iwp) :: i        !<
       INTEGER(iwp) :: ic       !<
       INTEGER(iwp) :: j        !<
       INTEGER(iwp) :: jc       !<
       INTEGER(iwp) :: jj       !<
       INTEGER(iwp) :: k        !<
       INTEGER(iwp) :: l        !<
       INTEGER(iwp) :: n        !<

       LOGICAL :: unroll        !<

       REAL(wp) ::  wall_left   !<
       REAL(wp) ::  wall_north  !<
       REAL(wp) ::  wall_right  !<
       REAL(wp) ::  wall_south  !<
       REAL(wp) ::  wall_total  !<
       REAL(wp) ::  wall_top    !<

       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                           nys_mg(grid_level)-1:nyn_mg(grid_level)+1,          &
                           nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  f_mg  !<
       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                           nys_mg(grid_level)-1:nyn_mg(grid_level)+1,          &
                           nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  p_mg  !<

       l = grid_level

!
!--    Choose flag array of this level
       SELECT CASE ( l )
          CASE ( 1 )
             flags => wall_flags_1
          CASE ( 2 )
             flags => wall_flags_2
          CASE ( 3 )
             flags => wall_flags_3
          CASE ( 4 )
             flags => wall_flags_4
          CASE ( 5 )
             flags => wall_flags_5
          CASE ( 6 )
             flags => wall_flags_6
          CASE ( 7 )
             flags => wall_flags_7
          CASE ( 8 )
             flags => wall_flags_8
          CASE ( 9 )
             flags => wall_flags_9
          CASE ( 10 )
             flags => wall_flags_10
       END SELECT

       unroll = ( MOD( nyn_mg(l)-nys_mg(l)+1, 4 ) == 0  .AND. &
                  MOD( nxr_mg(l)-nxl_mg(l)+1, 2 ) == 0 )

       DO  n = 1, ngsrb
          
          DO  color = 1, 2

             IF ( .NOT. unroll )  THEN

                CALL cpu_log( log_point_s(36), 'redblack_no_unroll_noopt', 'start' )

!
!--             Without unrolling of loops, no cache optimization
                DO  i = nxl_mg(l), nxr_mg(l), 2
                   DO  j = nys_mg(l) + 2 - color, nyn_mg(l), 2 
                      DO  k = nzb+1, nzt_mg(l), 2
!                      p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) * (                    &
!                                ddx2_mg(l) * ( p_mg(k,j,i+1) + p_mg(k,j,i-1) ) &
!                              + ddy2_mg(l) * ( p_mg(k,j+1,i) + p_mg(k,j-1,i) ) &
!                              + f2_mg(k,l) * p_mg(k+1,j,i)                     &
!                              + f3_mg(k,l) * p_mg(k-1,j,i) - f_mg(k,j,i)       &
!                                                       )

                         p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) * (                 &
                             rho_air_mg(k,l) * ddx2_mg(l) *                    &
                               ( p_mg(k,j,i+1) + IBITS( flags(k,j,i), 5, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i+1) ) + &
                                 p_mg(k,j,i-1) + IBITS( flags(k,j,i), 4, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i-1) ) ) &
                           + rho_air_mg(k,l) * ddy2_mg(l) *                    &
                               ( p_mg(k,j+1,i) + IBITS( flags(k,j,i), 3, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j+1,i) ) + &
                                 p_mg(k,j-1,i) + IBITS( flags(k,j,i), 2, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j-1,i) ) ) &
                           + f2_mg(k,l) *                                      &
                               ( p_mg(k+1,j,i) + IBITS( flags(k,j,i), 7, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k+1,j,i) ) ) &
                           + f3_mg(k,l) *                                      &
                               ( p_mg(k-1,j,i) + IBITS( flags(k,j,i), 0, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k-1,j,i) ) ) &
                           - f_mg(k,j,i)                    )
                      ENDDO
                   ENDDO
                ENDDO
      
                DO  i = nxl_mg(l)+1, nxr_mg(l), 2
                   DO  j = nys_mg(l) + (color-1), nyn_mg(l), 2
                      DO  k = nzb+1, nzt_mg(l), 2 
                         p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) * (                 &
                               rho_air_mg(k,l) * ddx2_mg(l) *                  &
                               ( p_mg(k,j,i+1) + IBITS( flags(k,j,i), 5, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i+1) ) + &
                                 p_mg(k,j,i-1) + IBITS( flags(k,j,i), 4, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i-1) ) ) &
                           + rho_air_mg(k,l) * ddy2_mg(l) *                    &
                               ( p_mg(k,j+1,i) + IBITS( flags(k,j,i), 3, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j+1,i) ) + &
                                 p_mg(k,j-1,i) + IBITS( flags(k,j,i), 2, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j-1,i) ) ) &
                           + f2_mg(k,l) *                                      &
                               ( p_mg(k+1,j,i) + IBITS( flags(k,j,i), 7, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k+1,j,i) ) ) &
                           + f3_mg(k,l) *                                      &
                               ( p_mg(k-1,j,i) + IBITS( flags(k,j,i), 0, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k-1,j,i) ) ) &
                           - f_mg(k,j,i)                     )
                      ENDDO
                   ENDDO
                ENDDO
     
                DO  i = nxl_mg(l), nxr_mg(l), 2
                   DO  j = nys_mg(l) + (color-1), nyn_mg(l), 2
                      DO  k = nzb+2, nzt_mg(l), 2
                        p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) * (                  &
                             rho_air_mg(k,l) * ddx2_mg(l) *                    &
                               ( p_mg(k,j,i+1) + IBITS( flags(k,j,i), 5, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i+1) ) + &
                                 p_mg(k,j,i-1) + IBITS( flags(k,j,i), 4, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i-1) ) ) &
                           + rho_air_mg(k,l) * ddy2_mg(l) *                    &
                               ( p_mg(k,j+1,i) + IBITS( flags(k,j,i), 3, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j+1,i) ) + &
                                 p_mg(k,j-1,i) + IBITS( flags(k,j,i), 2, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j-1,i) ) ) &
                           + f2_mg(k,l) *                                      &
                               ( p_mg(k+1,j,i) + IBITS( flags(k,j,i), 7, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k+1,j,i) ) ) &
                           + f3_mg(k,l) *                                      &
                               ( p_mg(k-1,j,i) + IBITS( flags(k,j,i), 0, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k-1,j,i) ) ) &
                           - f_mg(k,j,i)                    )
                      ENDDO
                   ENDDO
                ENDDO

                DO  i = nxl_mg(l)+1, nxr_mg(l), 2
                   DO  j = nys_mg(l) + 2 - color, nyn_mg(l), 2
                      DO  k = nzb+2, nzt_mg(l), 2
                         p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) * (                 &
                             rho_air_mg(k,l) * ddx2_mg(l) *                    &
                               ( p_mg(k,j,i+1) + IBITS( flags(k,j,i), 5, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i+1) ) + &
                                 p_mg(k,j,i-1) + IBITS( flags(k,j,i), 4, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i-1) ) ) &
                           + rho_air_mg(k,l) * ddy2_mg(l) *                    &
                               ( p_mg(k,j+1,i) + IBITS( flags(k,j,i), 3, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j+1,i) ) + &
                                 p_mg(k,j-1,i) + IBITS( flags(k,j,i), 2, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j-1,i) ) ) &
                           + f2_mg(k,l) *                                      &
                               ( p_mg(k+1,j,i) + IBITS( flags(k,j,i), 7, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k+1,j,i) ) ) &
                           + f3_mg(k,l) *                                      &
                               ( p_mg(k-1,j,i) + IBITS( flags(k,j,i), 0, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k-1,j,i) ) ) &
                           - f_mg(k,j,i)                    )
                      ENDDO
                   ENDDO
                ENDDO
                CALL cpu_log( log_point_s(36), 'redblack_no_unroll_noopt', 'stop' )

             ELSE

!
!--             Loop unrolling along y, only one i loop for better cache use
                CALL cpu_log( log_point_s(38), 'redblack_unroll_noopt', 'start' )
                DO  ic = nxl_mg(l), nxr_mg(l), 2
                   DO  jc = nys_mg(l), nyn_mg(l), 4
                      i  = ic
                      jj = jc+2-color
                      DO  k = nzb+1, nzt_mg(l), 2
                         j = jj
                         p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) * (                 &
                             rho_air_mg(k,l) * ddx2_mg(l) *                    &
                               ( p_mg(k,j,i+1) + IBITS( flags(k,j,i), 5, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i+1) ) + &
                                 p_mg(k,j,i-1) + IBITS( flags(k,j,i), 4, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i-1) ) ) &
                           + rho_air_mg(k,l) * ddy2_mg(l) *                    &
                               ( p_mg(k,j+1,i) + IBITS( flags(k,j,i), 3, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j+1,i) ) + &
                                 p_mg(k,j-1,i) + IBITS( flags(k,j,i), 2, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j-1,i) ) ) &
                           + f2_mg(k,l) *                                      &
                               ( p_mg(k+1,j,i) + IBITS( flags(k,j,i), 7, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k+1,j,i) ) ) &
                           + f3_mg(k,l) *                                      &
                               ( p_mg(k-1,j,i) + IBITS( flags(k,j,i), 0, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k-1,j,i) ) ) &
                           - f_mg(k,j,i)               )
                         j = jj+2
                         p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) * (                 &
                             rho_air_mg(k,l) * ddx2_mg(l) *                    &
                               ( p_mg(k,j,i+1) + IBITS( flags(k,j,i), 5, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i+1) ) + &
                                 p_mg(k,j,i-1) + IBITS( flags(k,j,i), 4, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i-1) ) ) &
                           + rho_air_mg(k,l) * ddy2_mg(l) *                    &
                               ( p_mg(k,j+1,i) + IBITS( flags(k,j,i), 3, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j+1,i) ) + &
                                 p_mg(k,j-1,i) + IBITS( flags(k,j,i), 2, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j-1,i) ) ) &
                           + f2_mg(k,l) *                                      &
                               ( p_mg(k+1,j,i) + IBITS( flags(k,j,i), 7, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k+1,j,i) ) ) &
                           + f3_mg(k,l) *                                      &
                               ( p_mg(k-1,j,i) + IBITS( flags(k,j,i), 0, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k-1,j,i) ) ) &
                           - f_mg(k,j,i)                     )
                      ENDDO
      
                      i  = ic+1
                      jj = jc+color-1
                      DO  k = nzb+1, nzt_mg(l), 2 
                         j =jj
                         p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) * (                 &
                             rho_air_mg(k,l) * ddx2_mg(l) *                    &
                               ( p_mg(k,j,i+1) + IBITS( flags(k,j,i), 5, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i+1) ) + &
                                 p_mg(k,j,i-1) + IBITS( flags(k,j,i), 4, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i-1) ) ) &
                           + rho_air_mg(k,l) * ddy2_mg(l) *                    &
                               ( p_mg(k,j+1,i) + IBITS( flags(k,j,i), 3, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j+1,i) ) + &
                                 p_mg(k,j-1,i) + IBITS( flags(k,j,i), 2, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j-1,i) ) ) &
                           + f2_mg(k,l) *                                      &
                               ( p_mg(k+1,j,i) + IBITS( flags(k,j,i), 7, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k+1,j,i) ) ) &
                           + f3_mg(k,l) *                                      &
                               ( p_mg(k-1,j,i) + IBITS( flags(k,j,i), 0, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k-1,j,i) ) ) &
                           - f_mg(k,j,i)                     )
                         j = jj+2
                         p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) * (                &
                            rho_air_mg(k,l) * ddx2_mg(l) *                    &
                              ( p_mg(k,j,i+1) + IBITS( flags(k,j,i), 5, 1 ) * &
                                            ( p_mg(k,j,i) - p_mg(k,j,i+1) ) + &
                                p_mg(k,j,i-1) + IBITS( flags(k,j,i), 4, 1 ) * &
                                            ( p_mg(k,j,i) - p_mg(k,j,i-1) ) ) &
                          + rho_air_mg(k,l) * ddy2_mg(l) *                    &
                              ( p_mg(k,j+1,i) + IBITS( flags(k,j,i), 3, 1 ) * &
                                            ( p_mg(k,j,i) - p_mg(k,j+1,i) ) + &
                                p_mg(k,j-1,i) + IBITS( flags(k,j,i), 2, 1 ) * &
                                            ( p_mg(k,j,i) - p_mg(k,j-1,i) ) ) &
                          + f2_mg(k,l) *                                      &
                              ( p_mg(k+1,j,i) + IBITS( flags(k,j,i), 7, 1 ) * &
                                            ( p_mg(k,j,i) - p_mg(k+1,j,i) ) ) &
                          + f3_mg(k,l) *                                      &
                              ( p_mg(k-1,j,i) + IBITS( flags(k,j,i), 0, 1 ) * &
                                            ( p_mg(k,j,i) - p_mg(k-1,j,i) ) ) &
                          - f_mg(k,j,i)                    )
                      ENDDO

                      i  = ic
                      jj = jc+color-1
                      DO  k = nzb+2, nzt_mg(l), 2
                         j =jj
                         p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) * (                 &
                             rho_air_mg(k,l) * ddx2_mg(l) *                    &
                               ( p_mg(k,j,i+1) + IBITS( flags(k,j,i), 5, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i+1) ) + &
                                 p_mg(k,j,i-1) + IBITS( flags(k,j,i), 4, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i-1) ) ) &
                           + rho_air_mg(k,l) * ddy2_mg(l) *                    &
                               ( p_mg(k,j+1,i) + IBITS( flags(k,j,i), 3, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j+1,i) ) + &
                                 p_mg(k,j-1,i) + IBITS( flags(k,j,i), 2, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j-1,i) ) ) &
                           + f2_mg(k,l) *                                      &
                               ( p_mg(k+1,j,i) + IBITS( flags(k,j,i), 7, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k+1,j,i) ) ) &
                           + f3_mg(k,l) *                                      &
                               ( p_mg(k-1,j,i) + IBITS( flags(k,j,i), 0, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k-1,j,i) ) ) &
                           - f_mg(k,j,i)                    )
                         j = jj+2
                         p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) * (                 &
                             rho_air_mg(k,l) * ddx2_mg(l) *                    &
                               ( p_mg(k,j,i+1) + IBITS( flags(k,j,i), 5, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i+1) ) + &
                                 p_mg(k,j,i-1) + IBITS( flags(k,j,i), 4, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i-1) ) ) &
                           + rho_air_mg(k,l) * ddy2_mg(l) *                    &
                               ( p_mg(k,j+1,i) + IBITS( flags(k,j,i), 3, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j+1,i) ) + &
                                 p_mg(k,j-1,i) + IBITS( flags(k,j,i), 2, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j-1,i) ) ) &
                           + f2_mg(k,l) *                                      &
                               ( p_mg(k+1,j,i) + IBITS( flags(k,j,i), 7, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k+1,j,i) ) ) &
                           + f3_mg(k,l) *                                      &
                               ( p_mg(k-1,j,i) + IBITS( flags(k,j,i), 0, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k-1,j,i) ) ) &
                           - f_mg(k,j,i)                     )
                      ENDDO

                      i  = ic+1
                      jj = jc+2-color
                      DO  k = nzb+2, nzt_mg(l), 2
                         j =jj
                         p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) * (                 &
                             rho_air_mg(k,l) * ddx2_mg(l) *                    &
                               ( p_mg(k,j,i+1) + IBITS( flags(k,j,i), 5, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i+1) ) + &
                                 p_mg(k,j,i-1) + IBITS( flags(k,j,i), 4, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i-1) ) ) &
                           + rho_air_mg(k,l) * ddy2_mg(l) *                    &
                               ( p_mg(k,j+1,i) + IBITS( flags(k,j,i), 3, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j+1,i) ) + &
                                 p_mg(k,j-1,i) + IBITS( flags(k,j,i), 2, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j-1,i) ) ) &
                           + f2_mg(k,l) *                                      &
                               ( p_mg(k+1,j,i) + IBITS( flags(k,j,i), 7, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k+1,j,i) ) ) &
                           + f3_mg(k,l) *                                      &
                               ( p_mg(k-1,j,i) + IBITS( flags(k,j,i), 0, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k-1,j,i) ) ) &
                           - f_mg(k,j,i)                     )
                         j = jj+2
                         p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) * (                 &
                             rho_air_mg(k,l) * ddx2_mg(l) *                    &
                               ( p_mg(k,j,i+1) + IBITS( flags(k,j,i), 5, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i+1) ) + &
                                 p_mg(k,j,i-1) + IBITS( flags(k,j,i), 4, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j,i-1) ) ) &
                           + rho_air_mg(k,l) * ddy2_mg(l) *                    &
                               ( p_mg(k,j+1,i) + IBITS( flags(k,j,i), 3, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j+1,i) ) + &
                                 p_mg(k,j-1,i) + IBITS( flags(k,j,i), 2, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k,j-1,i) ) ) &
                           + f2_mg(k,l) *                                      &
                               ( p_mg(k+1,j,i) + IBITS( flags(k,j,i), 7, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k+1,j,i) ) ) &
                           + f3_mg(k,l) *                                      &
                               ( p_mg(k-1,j,i) + IBITS( flags(k,j,i), 0, 1 ) * &
                                             ( p_mg(k,j,i) - p_mg(k-1,j,i) ) ) &
                           - f_mg(k,j,i)                     )
                      ENDDO

                   ENDDO
                ENDDO
                CALL cpu_log( log_point_s(38), 'redblack_unroll_noopt', 'stop' )

             ENDIF

!
!--          Horizontal boundary conditions
             CALL exchange_horiz( p_mg, 1 )

!
!--          Bottom and top boundary conditions
             IF ( ibc_p_b == 1 )  THEN
                p_mg(nzb,:,: ) = p_mg(nzb+1,:,:)
             ELSE
                p_mg(nzb,:,: ) = 0.0_wp
             ENDIF

             IF ( ibc_p_t == 1 )  THEN
                p_mg(nzt_mg(l)+1,:,: ) = p_mg(nzt_mg(l),:,:)
             ELSE
                p_mg(nzt_mg(l)+1,:,: ) = 0.0_wp
             ENDIF

          ENDDO

       ENDDO

!
!--    Set pressure within topography and at the topography surfaces
!$OMP PARALLEL PRIVATE (i,j,k,wall_left,wall_north,wall_right,wall_south,wall_top,wall_total)
!$OMP DO
       DO  i = nxl_mg(l), nxr_mg(l)
          DO  j = nys_mg(l), nyn_mg(l) 
             DO  k = nzb, nzt_mg(l)
!
!--             First, set pressure inside topography to zero
                p_mg(k,j,i) = p_mg(k,j,i) * ( 1.0_wp - IBITS( flags(k,j,i), 6, 1 ) )
!
!--             Second, determine if the gridpoint inside topography is adjacent
!--             to a wall and set its value to a value given by the average of
!--             those values obtained from Neumann boundary condition
                wall_left  = IBITS( flags(k,j,i-1), 5, 1 )
                wall_right = IBITS( flags(k,j,i+1), 4, 1 )
                wall_south = IBITS( flags(k,j-1,i), 3, 1 )
                wall_north = IBITS( flags(k,j+1,i), 2, 1 )
                wall_top   = IBITS( flags(k+1,j,i), 0, 1 )
                wall_total = wall_left + wall_right + wall_south + wall_north + &
                             wall_top

                IF ( wall_total > 0.0_wp )  THEN
                   p_mg(k,j,i) = 1.0_wp / wall_total *                 &
                                        ( wall_left  * p_mg(k,j,i-1) + &
                                          wall_right * p_mg(k,j,i+1) + &
                                          wall_south * p_mg(k,j-1,i) + &
                                          wall_north * p_mg(k,j+1,i) + &
                                          wall_top   * p_mg(k+1,j,i) )
                ENDIF
             ENDDO
          ENDDO
       ENDDO
!$OMP END PARALLEL

!
!--    One more time horizontal boundary conditions
       CALL exchange_horiz( p_mg, 1)


    END SUBROUTINE redblack_noopt



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Gather subdomain data from all PEs.
!------------------------------------------------------------------------------!
    SUBROUTINE mg_gather_noopt( f2, f2_sub )

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE indices,                                                            &
           ONLY:  mg_loc_ind, nxl_mg, nxr_mg, nys_mg, nyn_mg, nzb, nzt_mg

       USE kinds

       USE pegrid

       IMPLICIT NONE

       INTEGER(iwp) ::  i       !<
       INTEGER(iwp) ::  il      !<
       INTEGER(iwp) ::  ir      !<
       INTEGER(iwp) ::  j       !<
       INTEGER(iwp) ::  jn      !<
       INTEGER(iwp) ::  js      !<
       INTEGER(iwp) ::  k       !<
       INTEGER(iwp) ::  nwords  !<

       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                       nys_mg(grid_level)-1:nyn_mg(grid_level)+1,              &
                       nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  f2    !<
       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                       nys_mg(grid_level)-1:nyn_mg(grid_level)+1,              &
                       nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  f2_l  !<

       REAL(wp), DIMENSION(nzb:mg_loc_ind(5,myid)+1,                           &
                           mg_loc_ind(3,myid)-1:mg_loc_ind(4,myid)+1,          &
                           mg_loc_ind(1,myid)-1:mg_loc_ind(2,myid)+1) ::  f2_sub  !<


#if defined( __parallel )
       CALL cpu_log( log_point_s(34), 'mg_gather_noopt', 'start' )

       f2_l = 0.0_wp

!
!--    Store the local subdomain array on the total array
       js = mg_loc_ind(3,myid)
       IF ( south_border_pe )  js = js - 1
       jn = mg_loc_ind(4,myid)
       IF ( north_border_pe )  jn = jn + 1
       il = mg_loc_ind(1,myid)
       IF ( left_border_pe )   il = il - 1
       ir = mg_loc_ind(2,myid)
       IF ( right_border_pe )  ir = ir + 1
       DO  i = il, ir
          DO  j = js, jn
             DO  k = nzb, nzt_mg(grid_level)+1
                f2_l(k,j,i) = f2_sub(k,j,i)
             ENDDO
          ENDDO
       ENDDO

!
!--    Find out the number of array elements of the total array
       nwords = SIZE( f2 )

!
!--    Gather subdomain data from all PEs
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( f2_l(nzb,nys_mg(grid_level)-1,nxl_mg(grid_level)-1),&
                           f2(nzb,nys_mg(grid_level)-1,nxl_mg(grid_level)-1),  &
                           nwords, MPI_REAL, MPI_SUM, comm2d, ierr )

       CALL cpu_log( log_point_s(34), 'mg_gather_noopt', 'stop' )
#endif
       
    END SUBROUTINE mg_gather_noopt



!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo It may be possible to improve the speed of this routine by using
!>       non-blocking communication
!------------------------------------------------------------------------------!
    SUBROUTINE mg_scatter_noopt( p2, p2_sub )

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE indices,                                                            &
           ONLY:  mg_loc_ind, nxl_mg, nxr_mg, nys_mg, nyn_mg, nzb, nzt_mg

       USE kinds

       USE pegrid

       IMPLICIT NONE

       INTEGER(iwp) ::  nwords  !<

       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level-1)+1,                         &
                           nys_mg(grid_level-1)-1:nyn_mg(grid_level-1)+1,      &
                           nxl_mg(grid_level-1)-1:nxr_mg(grid_level-1)+1) ::  p2  !<

       REAL(wp), DIMENSION(nzb:mg_loc_ind(5,myid)+1,                           &
                           mg_loc_ind(3,myid)-1:mg_loc_ind(4,myid)+1,          &
                           mg_loc_ind(1,myid)-1:mg_loc_ind(2,myid)+1) ::  p2_sub  !<

!
!--    Find out the number of array elements of the subdomain array
       nwords = SIZE( p2_sub )

#if defined( __parallel )
       CALL cpu_log( log_point_s(35), 'mg_scatter_noopt', 'start' )

       p2_sub = p2(:,mg_loc_ind(3,myid)-1:mg_loc_ind(4,myid)+1, &
                     mg_loc_ind(1,myid)-1:mg_loc_ind(2,myid)+1)

       CALL cpu_log( log_point_s(35), 'mg_scatter_noopt', 'stop' )
#endif
       
    END SUBROUTINE mg_scatter_noopt


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This is where the multigrid technique takes place. V- and W- Cycle are
!> implemented and steered by the parameter "gamma". Parameter "nue" determines
!> the convergence of the multigrid iterative solution. There are nue times
!> RB-GS iterations. It should be set to "1" or "2", considering the time effort
!> one would like to invest. Last choice shows a very good converging factor,
!> but leads to an increase in computing time.
!------------------------------------------------------------------------------!
    RECURSIVE SUBROUTINE next_mg_level_noopt( f_mg, p_mg, p3, r )

       USE control_parameters,                                                 &
           ONLY:  bc_lr_dirrad, bc_lr_raddir, bc_ns_dirrad, bc_ns_raddir,      &
                  gamma_mg, grid_level_count, ibc_p_b, ibc_p_t,                &
                  maximum_grid_level,                                          &
                  mg_switch_to_pe0_level, mg_switch_to_pe0, ngsrb


       USE indices,                                                            &
           ONLY:  mg_loc_ind, nxl, nxl_mg, nxr, nxr_mg, nys, nys_mg, nyn,      &
                  nyn_mg, nzb, nzt, nzt_mg

       USE kinds

       USE pegrid

       IMPLICIT NONE

       INTEGER(iwp) ::  i            !<
       INTEGER(iwp) ::  j            !<
       INTEGER(iwp) ::  k            !<
       INTEGER(iwp) ::  nxl_mg_save  !<
       INTEGER(iwp) ::  nxr_mg_save  !<
       INTEGER(iwp) ::  nyn_mg_save  !<
       INTEGER(iwp) ::  nys_mg_save  !<
       INTEGER(iwp) ::  nzt_mg_save  !<

       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                           nys_mg(grid_level)-1:nyn_mg(grid_level)+1,          &
                           nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) :: f_mg  !<
       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                           nys_mg(grid_level)-1:nyn_mg(grid_level)+1,          &
                           nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) :: p_mg  !<
       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                           nys_mg(grid_level)-1:nyn_mg(grid_level)+1,          &
                           nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) :: p3    !<
       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                           nys_mg(grid_level)-1:nyn_mg(grid_level)+1,          &
                           nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) :: r     !<

       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level-1)+1,                         &
                           nys_mg(grid_level-1)-1:nyn_mg(grid_level-1)+1,      &
                           nxl_mg(grid_level-1)-1:nxr_mg(grid_level-1)+1) ::  f2  !<
       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level-1)+1,                         &
                           nys_mg(grid_level-1)-1:nyn_mg(grid_level-1)+1,      &
                           nxl_mg(grid_level-1)-1:nxr_mg(grid_level-1)+1) ::  p2  !<

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  f2_sub  !<
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  p2_sub  !<

!
!--    Restriction to the coarsest grid
    10 IF ( grid_level == 1 )  THEN

!
!--       Solution on the coarsest grid. Double the number of Gauss-Seidel
!--       iterations in order to get a more accurate solution.
          ngsrb = 2 * ngsrb

          CALL redblack_noopt( f_mg, p_mg )

          ngsrb = ngsrb / 2


       ELSEIF ( grid_level /= 1 )  THEN

          grid_level_count(grid_level) = grid_level_count(grid_level) + 1

!
!--       Solution on the actual grid level
          CALL redblack_noopt( f_mg, p_mg )

!
!--       Determination of the actual residual
          CALL resid_noopt( f_mg, p_mg, r )

!
!--       Restriction of the residual (finer grid values!) to the next coarser
!--       grid. Therefore, the grid level has to be decremented now. nxl..nzt have
!--       to be set to the coarse grid values, because these variables are needed
!--       for the exchange of ghost points in routine exchange_horiz
          grid_level = grid_level - 1
          nxl = nxl_mg(grid_level)
          nys = nys_mg(grid_level)
          nxr = nxr_mg(grid_level)
          nyn = nyn_mg(grid_level)
          nzt = nzt_mg(grid_level)

          IF ( grid_level == mg_switch_to_pe0_level )  THEN

!
!--          From this level on, calculations are done on PE0 only.
!--          First, carry out restriction on the subdomain.
!--          Therefore, indices of the level have to be changed to subdomain values
!--          in between (otherwise, the restrict routine would expect
!--          the gathered array)

             nxl_mg_save = nxl_mg(grid_level)
             nxr_mg_save = nxr_mg(grid_level)
             nys_mg_save = nys_mg(grid_level)
             nyn_mg_save = nyn_mg(grid_level)
             nzt_mg_save = nzt_mg(grid_level)
             nxl_mg(grid_level) = mg_loc_ind(1,myid)
             nxr_mg(grid_level) = mg_loc_ind(2,myid)
             nys_mg(grid_level) = mg_loc_ind(3,myid)
             nyn_mg(grid_level) = mg_loc_ind(4,myid)
             nzt_mg(grid_level) = mg_loc_ind(5,myid)
             nxl = mg_loc_ind(1,myid)
             nxr = mg_loc_ind(2,myid)
             nys = mg_loc_ind(3,myid)
             nyn = mg_loc_ind(4,myid)
             nzt = mg_loc_ind(5,myid)

             ALLOCATE( f2_sub(nzb:nzt_mg(grid_level)+1,                    &
                              nys_mg(grid_level)-1:nyn_mg(grid_level)+1,   &
                              nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) )

             CALL restrict_noopt( f2_sub, r )

!
!--          Restore the correct indices of this level
             nxl_mg(grid_level) = nxl_mg_save
             nxr_mg(grid_level) = nxr_mg_save
             nys_mg(grid_level) = nys_mg_save
             nyn_mg(grid_level) = nyn_mg_save
             nzt_mg(grid_level) = nzt_mg_save
             nxl = nxl_mg(grid_level)
             nxr = nxr_mg(grid_level)
             nys = nys_mg(grid_level)
             nyn = nyn_mg(grid_level)
             nzt = nzt_mg(grid_level)
!
!--          Gather all arrays from the subdomains on PE0
             CALL mg_gather_noopt( f2, f2_sub )

!
!--          Set switch for routine exchange_horiz, that no ghostpoint exchange
!--          has to be carried out from now on
             mg_switch_to_pe0 = .TRUE.

!
!--          In case of non-cyclic lateral boundary conditions, both in- and
!--          outflow conditions have to be used on all PEs after the switch,
!--          because then they have the total domain.
             IF ( bc_lr_dirrad )  THEN
                inflow_l  = .TRUE.
                inflow_r  = .FALSE.
                outflow_l = .FALSE.
                outflow_r = .TRUE.
             ELSEIF ( bc_lr_raddir )  THEN
                inflow_l  = .FALSE.
                inflow_r  = .TRUE.
                outflow_l = .TRUE.
                outflow_r = .FALSE.
             ELSEIF ( forcing )  THEN
                force_bound_l = .TRUE.
                force_bound_r = .TRUE.
             ENDIF

             IF ( bc_ns_dirrad )  THEN
                inflow_n  = .TRUE.
                inflow_s  = .FALSE.
                outflow_n = .FALSE.
                outflow_s = .TRUE.
             ELSEIF ( bc_ns_raddir )  THEN
                inflow_n  = .FALSE.
                inflow_s  = .TRUE.
                outflow_n = .TRUE.
                outflow_s = .FALSE.
             ELSEIF ( forcing )  THEN
                force_bound_s = .TRUE.
                force_bound_n = .TRUE.
             ENDIF

             DEALLOCATE( f2_sub )

          ELSE

             CALL restrict_noopt( f2, r )

          ENDIF

          p2 = 0.0_wp

!
!--       Repeat the same procedure till the coarsest grid is reached
          CALL next_mg_level_noopt( f2, p2, p3, r )

       ENDIF

!
!--    Now follows the prolongation
       IF ( grid_level >= 2 )  THEN

!
!--       Prolongation of the new residual. The values are transferred 
!--       from the coarse to the next finer grid.
          IF ( grid_level == mg_switch_to_pe0_level+1 )  THEN

#if defined( __parallel )
!
!--          At this level, the new residual first has to be scattered from
!--          PE0 to the other PEs
             ALLOCATE( p2_sub(nzb:mg_loc_ind(5,myid)+1,             &
                       mg_loc_ind(3,myid)-1:mg_loc_ind(4,myid)+1,   &
                       mg_loc_ind(1,myid)-1:mg_loc_ind(2,myid)+1) )

             CALL mg_scatter_noopt( p2, p2_sub )

!
!--          Therefore, indices of the previous level have to be changed to
!--          subdomain values in between (otherwise, the prolong routine would
!--          expect the gathered array)
             nxl_mg_save = nxl_mg(grid_level-1)
             nxr_mg_save = nxr_mg(grid_level-1)
             nys_mg_save = nys_mg(grid_level-1)
             nyn_mg_save = nyn_mg(grid_level-1)
             nzt_mg_save = nzt_mg(grid_level-1)
             nxl_mg(grid_level-1) = mg_loc_ind(1,myid)
             nxr_mg(grid_level-1) = mg_loc_ind(2,myid)
             nys_mg(grid_level-1) = mg_loc_ind(3,myid)
             nyn_mg(grid_level-1) = mg_loc_ind(4,myid)
             nzt_mg(grid_level-1) = mg_loc_ind(5,myid)

!
!--          Set switch for routine exchange_horiz, that ghostpoint exchange
!--          has to be carried again out from now on
             mg_switch_to_pe0 = .FALSE.

!
!--          For non-cyclic lateral boundary conditions and in case of nesting, 
!--          restore the in-/outflow conditions. 
             inflow_l  = .FALSE.;  inflow_r  = .FALSE.
             inflow_n  = .FALSE.;  inflow_s  = .FALSE.
             outflow_l = .FALSE.;  outflow_r = .FALSE.
             outflow_n = .FALSE.;  outflow_s = .FALSE.
!
             IF ( forcing )  THEN
                force_bound_l = .FALSE.
                force_bound_r = .FALSE.
                force_bound_s = .FALSE.
                force_bound_n = .FALSE.      
             ENDIF

             IF ( pleft == MPI_PROC_NULL )  THEN
                IF ( bc_lr_dirrad )  THEN
                   inflow_l  = .TRUE.
                ELSEIF ( bc_lr_raddir )  THEN
                   outflow_l = .TRUE.
                ELSEIF ( forcing )  THEN
                   force_bound_l = .TRUE.
                ENDIF
             ENDIF

             IF ( pright == MPI_PROC_NULL )  THEN
                IF ( bc_lr_dirrad )  THEN
                   outflow_r = .TRUE.
                ELSEIF ( bc_lr_raddir )  THEN
                   inflow_r  = .TRUE.
                ELSEIF ( forcing )  THEN
                   force_bound_r = .TRUE.
                ENDIF
             ENDIF

             IF ( psouth == MPI_PROC_NULL )  THEN
                IF ( bc_ns_dirrad )  THEN
                   outflow_s = .TRUE.
                ELSEIF ( bc_ns_raddir )  THEN
                   inflow_s  = .TRUE.
                ELSEIF ( forcing )  THEN
                   force_bound_s = .TRUE.
                ENDIF
             ENDIF

             IF ( pnorth == MPI_PROC_NULL )  THEN
                IF ( bc_ns_dirrad )  THEN
                   inflow_n  = .TRUE.
                ELSEIF ( bc_ns_raddir )  THEN
                   outflow_n = .TRUE.
                ELSEIF ( forcing )  THEN
                   force_bound_n = .TRUE.
                ENDIF
             ENDIF

             CALL prolong_noopt( p2_sub, p3 )

!
!--          Restore the correct indices of the previous level
             nxl_mg(grid_level-1) = nxl_mg_save
             nxr_mg(grid_level-1) = nxr_mg_save
             nys_mg(grid_level-1) = nys_mg_save
             nyn_mg(grid_level-1) = nyn_mg_save
             nzt_mg(grid_level-1) = nzt_mg_save

             DEALLOCATE( p2_sub )
#endif

          ELSE

             CALL prolong_noopt( p2, p3 )

          ENDIF

!
!--       Computation of the new pressure correction. Therefore, 
!--       values from prior grids are added up automatically stage by stage.
          DO  i = nxl_mg(grid_level)-1, nxr_mg(grid_level)+1
             DO  j = nys_mg(grid_level)-1, nyn_mg(grid_level)+1
                DO  k = nzb, nzt_mg(grid_level)+1 
                   p_mg(k,j,i) = p_mg(k,j,i) + p3(k,j,i)
                ENDDO
             ENDDO
          ENDDO

!
!--       Relaxation of the new solution
          CALL redblack_noopt( f_mg, p_mg )

       ENDIF


!
!--    The following few lines serve the steering of the multigrid scheme 
       IF ( grid_level == maximum_grid_level )  THEN

          GOTO 20

       ELSEIF ( grid_level /= maximum_grid_level  .AND.  grid_level /= 1  .AND. &
                grid_level_count(grid_level) /= gamma_mg )  THEN

          GOTO 10

       ENDIF

!
!--    Reset counter for the next call of poismg_noopt
       grid_level_count(grid_level) = 0

!
!--    Continue with the next finer level. nxl..nzt have to be
!--    set to the finer grid values, because these variables are needed for the
!--    exchange of ghost points in routine exchange_horiz
       grid_level = grid_level + 1
       nxl = nxl_mg(grid_level)
       nxr = nxr_mg(grid_level)
       nys = nys_mg(grid_level)
       nyn = nyn_mg(grid_level)
       nzt = nzt_mg(grid_level)

    20 CONTINUE

    END SUBROUTINE next_mg_level_noopt



    SUBROUTINE poismg_noopt_init

       USE control_parameters,                                                 &
           ONLY:  bc_lr_cyc, bc_ns_cyc, masking_method, maximum_grid_level

       USE indices,                                                            &
           ONLY:  flags, nxl_mg, nxr_mg, nyn_mg, nys_mg, nzb, nzt_mg,          &
                  wall_flags_0, wall_flags_1,                                  &
                  wall_flags_10, wall_flags_2, wall_flags_3,  wall_flags_4,    &
                  wall_flags_5, wall_flags_6, wall_flags_7, wall_flags_8,      &
                  wall_flags_9

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< index variable along x 
       INTEGER(iwp) ::  inc           !< incremental parameter for coarsening grid level
       INTEGER(iwp) ::  j             !< index variable along y
       INTEGER(iwp) ::  k             !< index variable along z
       INTEGER(iwp) ::  l             !< loop variable indication current grid level
       INTEGER(iwp) ::  nxl_l         !< index of left PE boundary for multigrid level
       INTEGER(iwp) ::  nxr_l         !< index of right PE boundary for multigrid level
       INTEGER(iwp) ::  nyn_l         !< index of north PE boundary for multigrid level
       INTEGER(iwp) ::  nys_l         !< index of south PE boundary for multigrid level
       INTEGER(iwp) ::  nzt_l         !< index of top PE boundary for multigrid level

       INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  topo_tmp
!
!--    Gridpoint increment of the current level. 
       inc = 1
       DO  l = maximum_grid_level, 1 , -1
!
!--       Set grid_level as it is required for exchange_horiz_2d_int
          grid_level = l

          nxl_l = nxl_mg(l)
          nxr_l = nxr_mg(l)
          nys_l = nys_mg(l)
          nyn_l = nyn_mg(l)
          nzt_l = nzt_mg(l)
!
!--       Assign the flag level to be calculated
          SELECT CASE ( l )
             CASE ( 1 )
                flags => wall_flags_1
             CASE ( 2 )
                flags => wall_flags_2
             CASE ( 3 )
                flags => wall_flags_3
             CASE ( 4 )
                flags => wall_flags_4
             CASE ( 5 )
                flags => wall_flags_5
             CASE ( 6 )
                flags => wall_flags_6
             CASE ( 7 )
                flags => wall_flags_7
             CASE ( 8 )
                flags => wall_flags_8
             CASE ( 9 )
                flags => wall_flags_9
             CASE ( 10 )
                flags => wall_flags_10
          END SELECT

!
!--       Depending on the grid level, set the respective bits in case of
!--       neighbouring walls
!--       Bit 0:  wall to the bottom
!--       Bit 1:  wall to the top (not realized in remaining PALM code so far)
!--       Bit 2:  wall to the south
!--       Bit 3:  wall to the north
!--       Bit 4:  wall to the left
!--       Bit 5:  wall to the right
!--       Bit 6:  inside building

          flags = 0
!
!--       In case of masking method, flags are not set and multigrid method
!--       works like FFT-solver
          IF ( .NOT. masking_method )  THEN

!
!--          Allocate temporary array for topography heights on coarser grid
!--          level. Please note, 2 ghoist points are required, in order to
!--          calculate flags() on the interior ghost point. 
             ALLOCATE( topo_tmp(nzb:nzt_l+1,nys_l-1:nyn_l+1,nxl_l-1:nxr_l+1) )
             topo_tmp = 0
             
             DO  i = nxl_l, nxr_l
                DO  j = nys_l, nyn_l
                   DO  k = nzb, nzt_l
                      topo_tmp(k,j,i) = wall_flags_0(k*inc,j*inc,i*inc)
                   ENDDO
                ENDDO
             ENDDO
             topo_tmp(nzt_l+1,:,:) = topo_tmp(nzt_l,:,:)
!
!--          Exchange ghost points on respective multigrid level. 2 ghost points
!--          are required, in order to calculate flags on 
!--          nys_l-1 / nyn_l+1 / nxl_l-1 / nxr_l+1. 
             CALL exchange_horiz_int( topo_tmp, nys_l, nyn_l, nxl_l, nxr_l, nzt_l, 1 )
!
             DO  i = nxl_l, nxr_l
                DO  j = nys_l, nyn_l
                   DO  k = nzb, nzt_l     
!
!--                   Inside/outside building (inside building does not need
!--                   further tests for walls)
                      IF ( .NOT. BTEST( topo_tmp(k,j,i), 0 ) )  THEN

                         flags(k,j,i) = IBSET( flags(k,j,i), 6 )

                      ELSE
!
!--                      Bottom wall
                         IF ( .NOT. BTEST( topo_tmp(k-1,j,i), 0 ) )  THEN
                            flags(k,j,i) = IBSET( flags(k,j,i), 0 )
                         ENDIF
!
!--                      South wall
                         IF ( .NOT. BTEST( topo_tmp(k,j-1,i), 0 ) )  THEN
                            flags(k,j,i) = IBSET( flags(k,j,i), 2 )
                         ENDIF
!
!--                      North wall
                         IF ( .NOT. BTEST( topo_tmp(k,j+1,i), 0 ) )  THEN
                            flags(k,j,i) = IBSET( flags(k,j,i), 3 )
                         ENDIF
!
!--                      Left wall
                         IF ( .NOT. BTEST( topo_tmp(k,j,i-1), 0 ) )  THEN
                            flags(k,j,i) = IBSET( flags(k,j,i), 4 )
                         ENDIF
!
!--                      Right wall
                         IF ( .NOT. BTEST( topo_tmp(k,j,i+1), 0 ) )  THEN
                            flags(k,j,i) = IBSET( flags(k,j,i), 5 )
                         ENDIF
!
!--                      Top wall
                         IF ( .NOT. BTEST( topo_tmp(k+1,j,i), 0 ) )  THEN
                            flags(k,j,i) = IBSET( flags(k,j,i), 7 )
                         ENDIF

                      ENDIF
                            
                   ENDDO
                ENDDO
             ENDDO
             flags(nzt_l+1,:,:) = flags(nzt_l,:,:)

             CALL exchange_horiz_int( flags, nys_l, nyn_l, nxl_l, nxr_l, nzt_l, 1 )

             DEALLOCATE( topo_tmp )

          ENDIF

          inc = inc * 2

       ENDDO
!
!--    Reset grid_level to "normal" grid
       grid_level = 0

    END SUBROUTINE poismg_noopt_init

 END MODULE poismg_noopt_mod
