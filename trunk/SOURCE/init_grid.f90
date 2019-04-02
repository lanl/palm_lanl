!> @file init_grid.f90
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
! $Id: init_grid.f90 3068 2018-06-12 14:49:41Z Giersch $
! New warning message concerning grid stretching has been introduced
! 
! 3066 2018-06-12 08:55:55Z Giersch
! Bugfix in IF statement before error message
! 
! 3065 2018-06-12 07:03:02Z Giersch
! New vertical stretching mechanism introduced
! 
! 3051 2018-05-30 17:43:55Z suehring
! Minor bugfix concerning mapping 3D buildings on top of terrain
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
! 
! 2968 2018-04-13 11:52:24Z suehring
! Bugfix in initialization in case of elevated model surface. Introduce
! index for minimum topography-top. 
! 
! 2955 2018-04-09 15:14:01Z suehring
! Improve topography filter routine and add ghost-point exchange for building 
! ID and building type. 
! 
! 2927 2018-03-23 15:13:00Z suehring
! Bugfix, setting boundary conditions for topography index array.
! 
! 2918 2018-03-21 15:52:14Z gronemeier
! Moved init_mixing_length to turbulence_closure_mod.f90
! 
! 2897 2018-03-15 11:47:16Z suehring
! Relax restrictions for topography input, terrain and building heights can be
! input separately and are not mandatory any more. 
! 
! 2893 2018-03-14 16:20:52Z suehring
! Revise informative message concerning filtered topography (1 grid-point 
! holes). 
! 
! 2892 2018-03-14 15:06:29Z suehring
! Bugfix, uninitialized array in case of single_building.
! 
! 2867 2018-03-09 09:40:23Z suehring
! Revise mapping of 3D buildings onto onto orography.
! 
! 2823 2018-02-20 15:31:45Z Giersch
! Set boundary conditions for 3D topography in case of non-cyclic boundary 
! conditions
! 
! 2796 2018-02-08 12:25:39Z suehring
! Bugfix in 3D building initialization
! 
! 2747 2018-01-15 12:44:17Z suehring
! Bugfix, topography height is rounded to the nearest discrete grid level
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2701 2017-12-15 15:40:50Z suehring
! Changes from last commit documented
! 
! 2698 2017-12-14 18:46:24Z suehring
! Bugfix in get_topography_top_index
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Revised topography input
! Set nzb_max not for the entire nest domain, only for boundary PEs
! Re-organize routine, split-up into several subroutines
! Modularize poismg_noopt
! Remove setting bit 26, 27, 28 in wall_flags_0, indicating former '_outer'
! arrays (not required any more).  
! Bugfix in generic tunnel setup (MS)
! 
! 2550 2017-10-16 17:12:01Z boeske
! Set lateral boundary conditions for topography on all three ghost layers
! 
! 2478 2017-09-18 13:37:24Z suehring
! Bugfix, correct flag for use_top
! 
! 2365 2017-08-21 14:59:59Z kanani
! Vertical nesting implemented (SadiqHuq)
! 
! 2319 2017-07-20 17:33:17Z suehring
! Remove print statements
! 
! 2318 2017-07-20 17:27:44Z suehring
! Get topography top index via Function call 
! 
! 2317 2017-07-20 17:27:19Z suehring
! Bugfixes in reading 3D topography from file
! 
! 2274 2017-06-09 13:27:48Z Giersch
! Changed error messages
! 
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! - Adjustments according to new topography representation
! - Bugfix: Move determination of nzb_max behind topography modification in 
!   cell-edge case 
! - Get rid off global arrays required for topography output
! - Enable topography input via netcdf 
! - Generic tunnel set-up added
! 
! 2200 2017-04-11 11:37:51Z suehring
! monotonic_adjustment removed
! 
! 2169 2017-03-06 18:16:35Z suehring
! Bugfix, move setting for topography grid convention to init_grid, else, if no
! value is set, the simulation may abort in case of restarts
! 
! 2128 2017-01-23 15:00:03Z suehring
! Bugfix in setting topography from file in case of ocean simulations
! 
! 2088 2016-12-19 16:30:25Z suehring
! Bugfix in generic topography in case of ocean simulations
! 
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
! 
! 2021 2016-10-07 14:08:57Z suehring
! Bugfix: setting Neumann boundary conditions for topography required for 
! topography flags in multigrid_noopt solver
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1994 2016-08-15 09:52:21Z suehring
! Bugfix in definition of generic topography
! 
! 1982 2016-08-01 11:04:48Z suehring
! Bugfix concering consistency check for topography
! 
! 1968 2016-07-18 12:01:49Z suehring
! Changed: PE-wise reading of topography file in order to avoid global definition
! of arrays nzb_local and nzb_tmp. Thereby, topography definition for single
! buildings and street canyons has changed, as well as flag setting for 
! multigrid scheme.
!
! Bugfix in checking l_grid anisotropy. 
! Simplify initial computation of lwall and vertical_influence, i.e. remove 
! nzb_s_inner as it is still zero at this point.
! 
! 1942 2016-06-14 12:18:18Z suehring
! Topography filter implemented to fill holes resolved by only one grid point.
! Initialization of flags for ws-scheme moved to advec_ws.  
! 
! 1931 2016-06-10 12:06:59Z suehring
! Rename multigrid into multigrid_noopt and multigrid_fast into multigrid
!
! 1910 2016-05-26 06:49:46Z raasch
! Bugfix: if topography is read from file, Neumann conditions are used for the
! nzb_local array (instead of cyclic conditions) in case that non-cyclic
! boundary conditions are switched on for the run
!
! 1902 2016-05-09 11:18:56Z suehring
! Set topography flags for multigrid solver only (not for multigrid_fast)
!
! 1886 2016-04-21 11:20:47Z suehring
! Bugfix: setting advection flags near walls
! reformulated index values for nzb_v_inner
! variable discriptions added in declaration block
!
! 1845 2016-04-08 08:29:13Z raasch
! nzb_2d removed
! 
! 1804 2016-04-05 16:30:18Z maronga
! Removed code for parameter file check (__check)
!
! 1779 2016-03-03 08:01:28Z raasch
! coupling_char is trimmed at every place it occurs, because it can have
! different length now
!
! 1762 2016-02-25 12:31:13Z hellstea
! Introduction of nested domain feature
!
! 1743 2016-01-13 10:23:51Z raasch
! Bugfix for calculation of nzb_s_outer and nzb_u_outer at north boundary of
! total domain
!
! 1691 2015-10-26 16:17:44Z maronga
! Renamed prandtl_layer to constant_flux_layer.
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1677 2015-10-02 13:25:23Z boeske
! Bugfix: Ghost points are included in wall_flags_0 and wall_flags_00
! 
! 1675 2015-10-02 08:28:59Z gronemeier
! Bugfix: Definition of topography grid levels
! 
! 1660 2015-09-21 08:15:16Z gronemeier
! Bugfix: Definition of topography grid levels if vertical grid stretching
!         starts below the maximum topography height.
!
! 1580 2015-04-10 13:43:49Z suehring
! Bugfix: setting flags for 5th order scheme near buildings
!
! 1575 2015-03-27 09:56:27Z raasch
! adjustments for psolver-queries
!
! 1557 2015-03-05 16:43:04Z suehring
! Adjustment for monotoinic limiter
!
! 1418 2014-06-06 13:05:08Z fricke
! Bugfix: Change if-condition for stretched grid in the ocean, with the old
!          condition and a negative value for dz_stretch_level the condition
!          was always true for the whole model domain
!
! 1409 2014-05-23 12:11:32Z suehring
! Bugfix: set wall_flags_0 at inflow and outflow boundary also for i <= nxlu 
! j <= nysv
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
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
! 
! 1221 2013-09-10 08:59:13Z raasch
! wall_flags_00 introduced to hold bits 32-63,
! additional 3D-flag arrays for replacing the 2D-index array nzb_s_inner in
! loops optimized for openACC (pres + flow_statistics)
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1069 2012-11-28 16:18:43Z maronga
! bugfix: added coupling_char to TOPOGRAPHY_DATA to allow topography in the
!         ocean model in case of coupled runs
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1015 2012-09-27 09:23:24Z raasch
! lower index for calculating wall_flags_0 set to nzb_w_inner instead of
! nzb_w_inner+1
!
! 996 2012-09-07 10:41:47Z raasch
! little reformatting
!
! 978 2012-08-09 08:28:32Z fricke
! Bugfix: nzb_max is set to nzt at non-cyclic lateral boundaries
! Bugfix: Set wall_flags_0 for inflow boundary
!
! 927 2012-06-06 19:15:04Z raasch
! Wall flags are not set for multigrid method in case of masking method
!
! 864 2012-03-27 15:10:33Z gryschka 
! In case of ocean and Dirichlet bottom bc for u and v dzu_mg and ddzu_pres
! were not correctly defined for k=1.
!
! 861 2012-03-26 14:18:34Z suehring
! Set wall_flags_0. The array is needed for degradation in ws-scheme near walls,
! inflow and outflow boundaries as well as near the bottom and the top of the
! model domain.!
! Initialization of nzb_s_inner and nzb_w_inner.
! gls has to be at least nbgp to do not exceed the array bounds of nzb_local
! while setting wall_flags_0
!
! 843 2012-02-29 15:16:21Z gryschka
! In case of ocean and dirichlet bc for u and v at the bottom
! the first u-level ist defined at same height as the first w-level
!
! 818 2012-02-08 16:11:23Z maronga
! Bugfix: topo_height is only required if topography is used. It is thus now
! allocated in the topography branch
!
! 809 2012-01-30 13:32:58Z maronga
! Bugfix: replaced .AND. and .NOT. with && and ! in the preprocessor directives
!
! 807 2012-01-25 11:53:51Z maronga
! New cpp directive "__check" implemented which is used by check_namelist_files
!
! Revision 1.1  1997/08/11 06:17:45  raasch
! Initial revision (Testversion)
!
!
! Description:
! -----------------------------------------------------------------------------!
!> Creating grid depending constants
!> @todo: Rearrange topo flag list
!> @todo: reference 3D buildings on top of orography is not tested and may need
!>        further improvement for steep slopes
!> @todo: Use more advanced setting of building type at filled holes 
!------------------------------------------------------------------------------!
 SUBROUTINE init_grid
 
    USE advec_ws,                                                              &
        ONLY:  ws_init_flags

    USE arrays_3d,                                                             &
        ONLY:  dd2zu, ddzu, ddzu_pres, ddzw, dzu, dzw, zu, zw
        
    USE control_parameters,                                                    &
        ONLY:  bc_lr_cyc, bc_ns_cyc, building_height, building_length_x,       &
               building_length_y, building_wall_left, building_wall_south,     &
               canyon_height, canyon_wall_left, canyon_wall_south,             &
               canyon_width_x, canyon_width_y, constant_flux_layer,            &
               dp_level_ind_b, dz, dz_max, dz_stretch_factor,                  &   
               dz_stretch_factor_array, dz_stretch_level, dz_stretch_level_end,&
               dz_stretch_level_end_index, dz_stretch_level_start_index,       &
               dz_stretch_level_start, grid_level,                             &
               force_bound_l, force_bound_r, force_bound_n, force_bound_s,     &
               ibc_uv_b, inflow_l, inflow_n, inflow_r, inflow_s,               &
               masking_method, maximum_grid_level, message_string,             &
               momentum_advec, nest_domain, nest_bound_l,                      &
               nest_bound_n, nest_bound_r, nest_bound_s,                       &
               number_stretch_level_end, number_stretch_level_start, ocean,    &
               outflow_l, outflow_n, outflow_r, outflow_s, psolver,            & 
               scalar_advec, topography, topography_grid_convention,           &
               tunnel_height, tunnel_length, tunnel_width_x, tunnel_width_y,   &
               tunnel_wall_depth, use_surface_fluxes, use_top_fluxes,          &
               wall_adjustment_factor
         
    USE grid_variables,                                                        &
        ONLY:  ddx, ddx2, ddy, ddy2, dx, dx2, dy, dy2, zu_s_inner, zw_w_inner
        
    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg, nz,   &
               nzb, nzb_max, nzb_s_inner, nzb_s_outer, nzb_u_inner,            &
               nzb_u_outer, nzb_v_inner, nzb_v_outer, nzb_w_inner,             &
               nzb_w_outer, nzt, topo_min_level
    
    USE kinds

    USE pegrid

    USE poismg_noopt_mod

    USE surface_mod,                                                           &
        ONLY:  get_topography_top_index, get_topography_top_index_ji, init_bc

    USE vertical_nesting_mod,                                                  &
        ONLY:  vnested, vnest_init_grid

    IMPLICIT NONE

    INTEGER(iwp) ::  i                           !< index variable along x 
    INTEGER(iwp) ::  j                           !< index variable along y
    INTEGER(iwp) ::  k                           !< index variable along z
    INTEGER(iwp) ::  k_top                       !< topography top index on local PE
    INTEGER(iwp) ::  n                           !< loop variable for stretching
    INTEGER(iwp) ::  number_dz                   !< number of user-specified dz values       
    INTEGER(iwp) ::  nzb_local_max               !< vertical grid index of maximum topography height
    INTEGER(iwp) ::  nzb_local_min               !< vertical grid index of minimum topography height
                                      
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  nzb_local  !< index for topography top at cell-center
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  nzb_tmp    !< dummy to calculate topography indices on u- and v-grid

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  topo !< input array for 3D topography and dummy array for setting "outer"-flags

    REAL(wp) ::  dz_level_end  !< distance between calculated height level for u/v-grid and user-specified end level for stretching
    REAL(wp) ::  dz_stretched  !< stretched vertical grid spacing
    
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  min_dz_stretch_level_end !< Array that contains all minimum heights where the stretching can end


!
!-- Calculation of horizontal array bounds including ghost layers
    nxlg = nxl - nbgp
    nxrg = nxr + nbgp
    nysg = nys - nbgp
    nyng = nyn + nbgp

!
!-- Allocate grid arrays
    ALLOCATE( ddzu(1:nzt+1), ddzw(1:nzt+1), dd2zu(1:nzt), dzu(1:nzt+1),        &
              dzw(1:nzt+1), zu(nzb:nzt+1), zw(nzb:nzt+1) )

!
!-- Compute height of u-levels from constant grid length and dz stretch factors
    IF ( dz(1) == -1.0_wp )  THEN
       message_string = 'missing dz'
       CALL message( 'init_grid', 'PA0200', 1, 2, 0, 6, 0 ) 
    ELSEIF ( dz(1) <= 0.0_wp )  THEN
       WRITE( message_string, * ) 'dz=',dz(1),' <= 0.0'
       CALL message( 'init_grid', 'PA0201', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Initialize dz_stretch_level_start with the value of dz_stretch_level
!-- if it was set by the user
    IF ( dz_stretch_level /= -9999999.9_wp ) THEN
       dz_stretch_level_start(1) = dz_stretch_level
    ENDIF
       
!
!-- Determine number of dz values and stretching levels specified by the
!-- user to allow right controlling of the stretching mechanism and to
!-- perform error checks
    number_dz = COUNT( dz /= -1.0_wp )
    number_stretch_level_start = COUNT( dz_stretch_level_start /=              &
                                       -9999999.9_wp )
    number_stretch_level_end = COUNT( dz_stretch_level_end /=                  &
                                      9999999.9_wp )

!
!-- The number of specified end levels +1 has to be the same than the number 
!-- of specified dz values
    IF ( number_dz /= number_stretch_level_end + 1 ) THEN
       WRITE( message_string, * ) 'The number of values for dz = ',         &
                                   number_dz, 'has to be the same than& ',  &
                                   'the number of values for ',             &
                                   'dz_stretch_level_end + 1 = ',           &
                                   number_stretch_level_end+1
          CALL message( 'init_grid', 'PA0156', 1, 2, 0, 6, 0 )
    ENDIF
    
!
!--    The number of specified start levels has to be the same or one less than 
!--    the number of specified dz values
    IF ( number_dz /= number_stretch_level_start + 1 .AND.                  &
         number_dz /= number_stretch_level_start ) THEN
       WRITE( message_string, * ) 'The number of values for dz = ',         &
                                   number_dz, 'has to be the same or one ', &
                                   'more than& the number of values for ',  &
                                   'dz_stretch_level_start = ',             &
                                   number_stretch_level_start
          CALL message( 'init_grid', 'PA0211', 1, 2, 0, 6, 0 )
    ENDIF
    
!--    The number of specified start levels has to be the same or one more than 
!--    the number of specified end levels
    IF ( number_stretch_level_start /= number_stretch_level_end + 1 .AND.   &
         number_stretch_level_start /= number_stretch_level_end ) THEN
       WRITE( message_string, * ) 'The number of values for ',              &
                                  'dz_stretch_level_start = ',              &
                                   dz_stretch_level_start, 'has to be the ',&
                                   'same or one more than& the number of ', &
                                   'values for dz_stretch_level_end = ',    &
                                   number_stretch_level_end
          CALL message( 'init_grid', 'PA0216', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Initialize dz for the free atmosphere with the value of dz_max
    IF ( dz(number_stretch_level_start+1) == -1.0_wp .AND.                     &
         number_stretch_level_start /= 0 ) THEN 
       dz(number_stretch_level_start+1) = dz_max
    ENDIF
       
!
!-- Initialize the stretching factor if (infinitely) stretching in the free 
!-- atmosphere is desired (dz_stretch_level_end was not specified for the 
!-- free atmosphere)
    IF ( number_stretch_level_start == number_stretch_level_end + 1 ) THEN 
       dz_stretch_factor_array(number_stretch_level_start) =                   &
       dz_stretch_factor
    ENDIF
    
!
!-- Allocation of arrays for stretching
    ALLOCATE( min_dz_stretch_level_end(number_stretch_level_start) )

!
!-- Define the vertical grid levels
    IF ( .NOT. ocean )  THEN
    
!
!--    The stretching region has to be large enough to allow for a smooth
!--    transition between two different grid spacings
       DO n = 1, number_stretch_level_start
          min_dz_stretch_level_end(n) = dz_stretch_level_start(n) +            &
                                        4 * MAX( dz(n),dz(n+1) )
       ENDDO

       IF ( ANY( min_dz_stretch_level_end(1:number_stretch_level_start) >      &
                 dz_stretch_level_end(1:number_stretch_level_start) ) ) THEN
             message_string= 'Eeach dz_stretch_level_end has to be larger ' // &
                             'than its corresponding value for &' //           &
                             'dz_stretch_level_start + 4*MAX(dz(n),dz(n+1)) '//&
                             'to allow for smooth grid stretching'
             CALL message( 'init_grid', 'PA0224', 1, 2, 0, 6, 0 )
       ENDIF
       
!
!--    Stretching must not be applied within the prandtl_layer 
!--    (first two grid points). For the default case dz_stretch_level_start 
!--    is negative. Therefore the absolut value is checked here.
       IF ( ANY( ABS( dz_stretch_level_start ) < dz(1) * 1.5_wp ) ) THEN
          WRITE( message_string, * ) 'Eeach dz_stretch_level_start has to be ',&
                                     'larger than ', dz(1) * 1.5
             CALL message( 'init_grid', 'PA0226', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    The stretching has to start and end on a grid level. Therefore 
!--    user-specified values have to ''interpolate'' to the next lowest level
       IF ( number_stretch_level_start /= 0 ) THEN
          dz_stretch_level_start(1) = INT( (dz_stretch_level_start(1) -        &
                                            dz(1)/2.0) / dz(1) )               &
                                      * dz(1) + dz(1)/2.0
       ENDIF
       
       IF ( number_stretch_level_start > 1 ) THEN
          DO n = 2, number_stretch_level_start
             dz_stretch_level_start(n) = INT( dz_stretch_level_start(n) /      &
                                              dz(n) ) * dz(n)
          ENDDO
       ENDIF
       
       IF ( number_stretch_level_end /= 0 ) THEN
          DO n = 1, number_stretch_level_end
             dz_stretch_level_end(n) = INT( dz_stretch_level_end(n) /          &
                                            dz(n+1) ) * dz(n+1)
          ENDDO
       ENDIF
 
!
!--    Determine stretching factor if necessary
       IF ( number_stretch_level_end >= 1 ) THEN 
          CALL calculate_stretching_factor( number_stretch_level_end )
       ENDIF

!
!--    Grid for atmosphere with surface at z=0 (k=0, w-grid).
!--    First compute the u- and v-levels. In case of dirichlet bc for u and v
!--    the first u/v- and w-level (k=0) are defined at same height (z=0). 
!--    The second u-level (k=1) corresponds to the top of the 
!--    Prandtl-layer.
       IF ( ibc_uv_b == 0 .OR. ibc_uv_b == 2 ) THEN
          zu(0) = 0.0_wp
       ELSE
          zu(0) = - dz(1) * 0.5_wp
       ENDIF
          
       zu(1) =   dz(1) * 0.5_wp
       
!
!--    Determine u and v height levels considering the possibility of grid
!--    stretching in several heights.
       n = 1
       dz_stretch_level_start_index = nzt+1
       dz_stretch_level_end_index = nzt+1
       dz_stretched = dz(1)

!--    The default value of dz_stretch_level_start is negative, thus the first
!--    condition is always true. Hence, the second condition is necessary.
       DO  k = 2, nzt+1
          IF ( dz_stretch_level_start(n) <= zu(k-1) .AND.                      &
               dz_stretch_level_start(n) /= -9999999.9_wp ) THEN
             dz_stretched = dz_stretched * dz_stretch_factor_array(n)
             
             IF ( dz(n) > dz(n+1) ) THEN
                dz_stretched = MAX( dz_stretched, dz(n+1) ) !Restrict dz_stretched to the user-specified (higher) dz
             ELSE
                dz_stretched = MIN( dz_stretched, dz(n+1) ) !Restrict dz_stretched to the user-specified (lower) dz
             ENDIF
             
             IF ( dz_stretch_level_start_index(n) == nzt+1 )                         &
             dz_stretch_level_start_index(n) = k-1
             
          ENDIF
          
          zu(k) = zu(k-1) + dz_stretched
          
!
!--       Make sure that the stretching ends exactly at dz_stretch_level_end 
          dz_level_end = ABS( zu(k) - dz_stretch_level_end(n) ) 
          
          IF ( dz_level_end  < dz(n+1)/3.0 ) THEN
             zu(k) = dz_stretch_level_end(n)
             dz_stretched = dz(n+1)
             dz_stretch_level_end_index(n) = k
             n = n + 1             
          ENDIF
       ENDDO

!
!--    Compute the w-levels. They are always staggered half-way between the 
!--    corresponding u-levels. In case of dirichlet bc for u and v at the 
!--    ground the first u- and w-level (k=0) are defined at same height (z=0). 
!--    The top w-level is extrapolated linearly.
       zw(0) = 0.0_wp
       DO  k = 1, nzt
          zw(k) = ( zu(k) + zu(k+1) ) * 0.5_wp
       ENDDO
       zw(nzt+1) = zw(nzt) + 2.0_wp * ( zu(nzt+1) - zw(nzt) )

    ELSE

!
!--    The stretching region has to be large enough to allow for a smooth
!--    transition between two different grid spacings
       DO n = 1, number_stretch_level_start
          min_dz_stretch_level_end(n) = dz_stretch_level_start(n) -            &
                                        4 * MAX( dz(n),dz(n+1) )
       ENDDO
       
       IF ( ANY( min_dz_stretch_level_end (1:number_stretch_level_start) <     &
                 dz_stretch_level_end(1:number_stretch_level_start) ) ) THEN
             message_string= 'Eeach dz_stretch_level_end has to be less ' //   &
                             'than its corresponding value for &' //           &
                             'dz_stretch_level_start - 4*MAX(dz(n),dz(n+1)) '//&
                             'to allow for smooth grid stretching'
             CALL message( 'init_grid', 'PA0224', 1, 2, 0, 6, 0 )
       ENDIF
       
!
!--    Stretching must not be applied close to the surface (last two grid 
!--    points). For the default case dz_stretch_level_start is negative. 
       IF ( ANY( dz_stretch_level_start > - dz(1) * 1.5_wp ) ) THEN
          WRITE( message_string, * ) 'Eeach dz_stretch_level_start has to be ',&
                                     'less than ', dz(1) * 1.5
             CALL message( 'init_grid', 'PA0226', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    The stretching has to start and end on a grid level. Therefore 
!--    user-specified values have to ''interpolate'' to the next highest level
       IF ( number_stretch_level_start /= 0 ) THEN
          dz_stretch_level_start(1) = INT( (dz_stretch_level_start(1) +        &
                                            dz(1)/2.0) / dz(1) )               &
                                      * dz(1) - dz(1)/2.0
       ENDIF
       
       IF ( number_stretch_level_start > 1 ) THEN
          DO n = 2, number_stretch_level_start
             dz_stretch_level_start(n) = INT( dz_stretch_level_start(n) /      &
                                              dz(n) ) * dz(n)
          ENDDO
       ENDIF
       
       IF ( number_stretch_level_end /= 0 ) THEN
          DO n = 1, number_stretch_level_end
             dz_stretch_level_end(n) = INT( dz_stretch_level_end(n) /          &
                                            dz(n+1) ) * dz(n+1)
          ENDDO
       ENDIF
       
!
!--    Determine stretching factor if necessary
       IF ( number_stretch_level_end >= 1 ) THEN 
          CALL calculate_stretching_factor( number_stretch_level_end )
       ENDIF

!
!--    Grid for ocean with free water surface is at k=nzt (w-grid).
!--    In case of neumann bc at the ground the first first u-level (k=0) lies
!--    below the first w-level (k=0). In case of dirichlet bc the first u- and
!--    w-level are defined at same height, but staggered from the second level.
!--    The second u-level (k=1) corresponds to the top of the Prandtl-layer.
!--    z values are negative starting from z=0 (surface)
       zu(nzt+1) =   dz(1) * 0.5_wp
       zu(nzt)   = - dz(1) * 0.5_wp

!
!--    Determine u and v height levels considering the possibility of grid
!--    stretching in several heights.
       n = 1
       dz_stretch_level_start_index = 0
       dz_stretch_level_end_index = 0
       dz_stretched = dz(1)

       DO  k = nzt-1, 0, -1
          
          IF ( dz_stretch_level_start(n) >= zu(k+1) ) THEN
             dz_stretched = dz_stretched * dz_stretch_factor_array(n)

             IF ( dz(n) > dz(n+1) ) THEN
                dz_stretched = MAX( dz_stretched, dz(n+1) ) !Restrict dz_stretched to the user-specified (higher) dz
             ELSE
                dz_stretched = MIN( dz_stretched, dz(n+1) ) !Restrict dz_stretched to the user-specified (lower) dz
             ENDIF
             
             IF ( dz_stretch_level_start_index(n) == 0 )                             &
             dz_stretch_level_start_index(n) = k+1
             
          ENDIF
          
          zu(k) = zu(k+1) - dz_stretched
          
!
!--       Make sure that the stretching ends exactly at dz_stretch_level_end 
          dz_level_end = ABS( zu(k) - dz_stretch_level_end(n) ) 
          
          IF ( dz_level_end  < dz(n+1)/3.0 ) THEN
             zu(k) = dz_stretch_level_end(n)
             dz_stretched = dz(n+1)
             dz_stretch_level_end_index(n) = k
             n = n + 1             
          ENDIF
       ENDDO
       
!
!--    Compute the w-levels. They are always staggered half-way between the 
!--    corresponding u-levels, except in case of dirichlet bc for u and v
!--    at the ground. In this case the first u- and w-level are defined at 
!--    same height. The top w-level (nzt+1) is not used but set for 
!--    consistency, since w and all scalar variables are defined up tp nzt+1.
       zw(nzt+1) = dz(1)
       zw(nzt)   = 0.0_wp
       DO  k = 0, nzt
          zw(k) = ( zu(k) + zu(k+1) ) * 0.5_wp
       ENDDO

!
!--    In case of dirichlet bc for u and v the first u- and w-level are defined
!--    at same height.
       IF ( ibc_uv_b == 0 ) THEN
          zu(0) = zw(0)
       ENDIF

    ENDIF

!
!-- Compute grid lengths.
    DO  k = 1, nzt+1
       dzu(k)  = zu(k) - zu(k-1)
       ddzu(k) = 1.0_wp / dzu(k)
       dzw(k)  = zw(k) - zw(k-1)
       ddzw(k) = 1.0_wp / dzw(k)
    ENDDO

    DO  k = 1, nzt
       dd2zu(k) = 1.0_wp / ( dzu(k) + dzu(k+1) )
    ENDDO
    
!    
!-- The FFT- SOR-pressure solvers assume grid spacings of a staggered grid
!-- everywhere. For the actual grid, the grid spacing at the lowest level
!-- is only dz/2, but should be dz. Therefore, an additional array
!-- containing with appropriate grid information is created for these
!-- solvers.
    IF ( psolver(1:9) /= 'multigrid' )  THEN
       ALLOCATE( ddzu_pres(1:nzt+1) )
       ddzu_pres = ddzu
       ddzu_pres(1) = ddzu_pres(2)  ! change for lowest level
    ENDIF

!
!-- Compute the reciprocal values of the horizontal grid lengths.
    ddx = 1.0_wp / dx
    ddy = 1.0_wp / dy
    dx2 = dx * dx
    dy2 = dy * dy
    ddx2 = 1.0_wp / dx2
    ddy2 = 1.0_wp / dy2

!
!-- Allocate 3D array to set topography
    ALLOCATE( topo(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    topo = 0
!
!-- Initialize topography by generic topography or read from topography from file.  
    CALL init_topo( topo )
!
!-- Set flags to mask topography on the grid. 
    CALL set_topo_flags( topo )    
!
!-- Calculate wall flag arrays for the multigrid method.
!-- Please note, wall flags are only applied in the non-optimized version.
    IF ( psolver == 'multigrid_noopt' )  CALL poismg_noopt_init 

!
!-- Init flags for ws-scheme to degrade order of the numerics near walls, i.e. 
!-- to decrease the numerical stencil appropriately.
    IF ( momentum_advec == 'ws-scheme'  .OR.  scalar_advec == 'ws-scheme' )    &
       CALL ws_init_flags

!
!-- Determine the maximum level of topography. It is used for
!-- steering the degradation of order of the applied advection scheme, 
!-- as well in the lpm. 
!-- In case of non-cyclic lateral boundaries, the order of the advection
!-- scheme has to be reduced up to nzt (required at the lateral boundaries).
    k_top = 0
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb, nzt + 1
             k_top = MAX( k_top, MERGE( k, 0,                                  &
                                        .NOT. BTEST( topo(k,j,i), 0 ) ) )
          ENDDO
       ENDDO
    ENDDO
#if defined( __parallel )
    CALL MPI_ALLREDUCE( k_top + 1, nzb_max, 1, MPI_INTEGER,                    & !is +1 really necessary here?
                        MPI_MAX, comm2d, ierr )
#else
    nzb_max = k_top + 1
#endif
    IF ( inflow_l  .OR.  outflow_l  .OR.  force_bound_l  .OR.  nest_bound_l  .OR.&
         inflow_r  .OR.  outflow_r  .OR.  force_bound_r  .OR.  nest_bound_r  .OR.&
         inflow_n  .OR.  outflow_n  .OR.  force_bound_n  .OR.  nest_bound_n  .OR.&
         inflow_s  .OR.  outflow_s  .OR.  force_bound_s  .OR.  nest_bound_s )    &
         nzb_max = nzt
!   
!-- Finally, if topography extents up to the model top, limit nzb_max to nzt.
    nzb_max = MIN( nzb_max, nzt )
!
!-- Determine minimum index of topography. Usually, this will be nzb. In case
!-- there is elevated topography, however, the lowest topography will be higher. 
!-- This index is e.g. used to calculate mean first-grid point atmosphere 
!-- temperature, surface pressure and density, etc. .
    topo_min_level   = 0
#if defined( __parallel )
    CALL MPI_ALLREDUCE( MINVAL( get_topography_top_index( 's' ) ),             &
                        topo_min_level, 1, MPI_INTEGER, MPI_MIN, comm2d, ierr )
#else
    topo_min_level = MINVAL( get_topography_top_index( 's' ) )
#endif
!
!-- Initialize boundary conditions via surface type
    CALL init_bc
!
!-- Allocate and set topography height arrays required for data output
    IF ( TRIM( topography ) /= 'flat' )  THEN
!
!--    Allocate and set the arrays containing the topography height
       IF ( nxr == nx  .AND.  nyn /= ny )  THEN
          ALLOCATE( zu_s_inner(nxl:nxr+1,nys:nyn),                             &
                    zw_w_inner(nxl:nxr+1,nys:nyn) )
       ELSEIF ( nxr /= nx  .AND.  nyn == ny )  THEN
          ALLOCATE( zu_s_inner(nxl:nxr,nys:nyn+1),                             &
                    zw_w_inner(nxl:nxr,nys:nyn+1) )
       ELSEIF ( nxr == nx  .AND.  nyn == ny )  THEN
          ALLOCATE( zu_s_inner(nxl:nxr+1,nys:nyn+1),                           &
                    zw_w_inner(nxl:nxr+1,nys:nyn+1) )
       ELSE
          ALLOCATE( zu_s_inner(nxl:nxr,nys:nyn),                               &
                    zw_w_inner(nxl:nxr,nys:nyn) )
       ENDIF

       zu_s_inner   = 0.0_wp
       zw_w_inner   = 0.0_wp
!
!--    Determine local topography height on scalar and w-grid. Note, setting
!--    lateral boundary values is not necessary, realized via wall_flags_0 
!--    array. Further, please note that loop bounds are different from 
!--    nxl to nxr and nys to nyn on south and right model boundary, hence, 
!--    use intrinsic lbound and ubound functions to infer array bounds.
       DO  i = LBOUND(zu_s_inner, 1), UBOUND(zu_s_inner, 1)
          DO  j = LBOUND(zu_s_inner, 2), UBOUND(zu_s_inner, 2)
!
!--          Topography height on scalar grid. Therefore, determine index of 
!--          upward-facing surface element on scalar grid. 
             zu_s_inner(i,j) = zu( get_topography_top_index_ji( j, i, 's' ) )
!
!--          Topography height on w grid. Therefore, determine index of 
!--          upward-facing surface element on w grid.
             zw_w_inner(i,j) = zw( get_topography_top_index_ji( j, i, 's' ) )
          ENDDO
       ENDDO
    ENDIF

!
!-- In the following, calculate 2D index arrays. Note, these will be removed
!-- soon. 
!-- Allocate outer and inner index arrays for topography and set
!-- defaults.                    
    ALLOCATE( nzb_s_inner(nysg:nyng,nxlg:nxrg),                                &
              nzb_s_outer(nysg:nyng,nxlg:nxrg),                                &
              nzb_u_inner(nysg:nyng,nxlg:nxrg),                                &
              nzb_u_outer(nysg:nyng,nxlg:nxrg),                                &
              nzb_v_inner(nysg:nyng,nxlg:nxrg),                                &
              nzb_v_outer(nysg:nyng,nxlg:nxrg),                                &
              nzb_w_inner(nysg:nyng,nxlg:nxrg),                                &
              nzb_w_outer(nysg:nyng,nxlg:nxrg),                                &
              nzb_local(nysg:nyng,nxlg:nxrg),                                  &
              nzb_tmp(nysg:nyng,nxlg:nxrg) )
!
!-- Initialize 2D-index arrays. Note, these will be removed soon! 
    nzb_local(nys:nyn,nxl:nxr) = get_topography_top_index( 's' )
    CALL exchange_horiz_2d_int( nzb_local, nys, nyn, nxl, nxr, nbgp )
!
!-- Check topography for consistency with model domain. Therefore, use
!-- maximum and minium topography-top indices. Note, minimum topography top
!-- index is already calculated.  
    IF ( TRIM( topography ) /= 'flat' )  THEN
#if defined( __parallel )
       CALL MPI_ALLREDUCE( MAXVAL( get_topography_top_index( 's' ) ),          &
                           nzb_local_max, 1, MPI_INTEGER, MPI_MAX, comm2d, ierr )             
#else
       nzb_local_max = MAXVAL( get_topography_top_index( 's' ) )
#endif
       nzb_local_min = topo_min_level
!
!--    Consistency checks
       IF ( nzb_local_min < 0  .OR.  nzb_local_max  > nz + 1 )  THEN
          WRITE( message_string, * ) 'nzb_local values are outside the',       &
                                ' model domain',                               &
                                '&MINVAL( nzb_local ) = ', nzb_local_min,      &
                                '&MAXVAL( nzb_local ) = ', nzb_local_max
          CALL message( 'init_grid', 'PA0210', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    nzb_s_inner = nzb;  nzb_s_outer = nzb
    nzb_u_inner = nzb;  nzb_u_outer = nzb
    nzb_v_inner = nzb;  nzb_v_outer = nzb
    nzb_w_inner = nzb;  nzb_w_outer = nzb


!
!-- Set Neumann conditions for topography. Will be removed soon.
    IF ( .NOT. bc_ns_cyc )  THEN
       IF ( nys == 0  )  THEN
          DO  i = 1, nbgp  
             nzb_local(nys-i,:)   = nzb_local(nys,:)
          ENDDO
       ELSEIF ( nyn == ny )  THEN
          DO  i = 1, nbgp  
             nzb_local(ny+i,:) = nzb_local(ny,:)
          ENDDO
       ENDIF
    ENDIF

    IF ( .NOT. bc_lr_cyc )  THEN
       IF ( nxl == 0  )  THEN
          DO  i = 1, nbgp  
             nzb_local(:,nxl-i)   = nzb_local(:,nxl)
          ENDDO
       ELSEIF ( nxr == nx )  THEN
          DO  i = 1, nbgp  
             nzb_local(:,nx+i) = nzb_local(:,nx)
          ENDDO 
       ENDIF          
    ENDIF
!
!-- Initialization of 2D index arrays, will be removed soon!
!-- Initialize nzb_s_inner and nzb_w_inner
    nzb_s_inner = nzb_local
    nzb_w_inner = nzb_local

!
!-- Initialize remaining index arrays:
!-- first pre-initialize them with nzb_s_inner...
    nzb_u_inner = nzb_s_inner
    nzb_u_outer = nzb_s_inner
    nzb_v_inner = nzb_s_inner
    nzb_v_outer = nzb_s_inner
    nzb_w_outer = nzb_s_inner
    nzb_s_outer = nzb_s_inner

!
!-- nzb_s_outer: 
!-- extend nzb_local east-/westwards first, then north-/southwards
    nzb_tmp = nzb_local
    DO  j = nys, nyn
       DO  i = nxl, nxr
          nzb_tmp(j,i) = MAX( nzb_local(j,i-1), nzb_local(j,i),             &
                              nzb_local(j,i+1) )
       ENDDO
    ENDDO
       
    CALL exchange_horiz_2d_int( nzb_tmp, nys, nyn, nxl, nxr, nbgp )
     
    DO  i = nxl, nxr
       DO  j = nys, nyn
          nzb_s_outer(j,i) = MAX( nzb_tmp(j-1,i), nzb_tmp(j,i),             &
                                  nzb_tmp(j+1,i) )
       ENDDO
!
!--    non-cyclic boundary conditions (overwritten by call of
!--    exchange_horiz_2d_int below in case of cyclic boundary conditions)
       IF ( nys == 0 )  THEN
          j = -1
          nzb_s_outer(j,i) = MAX( nzb_tmp(j+1,i), nzb_tmp(j,i) )
       ENDIF
       IF ( nyn == ny )  THEN
          j = ny + 1
          nzb_s_outer(j,i) = MAX( nzb_tmp(j-1,i), nzb_tmp(j,i) )
       ENDIF
    ENDDO
!
!-- nzb_w_outer: 
!-- identical to nzb_s_outer
    nzb_w_outer = nzb_s_outer
!
!-- nzb_u_inner: 
!-- extend nzb_local rightwards only
    nzb_tmp = nzb_local
    DO  j = nys, nyn
       DO  i = nxl, nxr
          nzb_tmp(j,i) = MAX( nzb_local(j,i-1), nzb_local(j,i) )
       ENDDO
    ENDDO
       
    CALL exchange_horiz_2d_int( nzb_tmp, nys, nyn, nxl, nxr, nbgp )
       
    nzb_u_inner = nzb_tmp
!
!-- nzb_u_outer: 
!-- extend current nzb_tmp (nzb_u_inner) north-/southwards
    DO  i = nxl, nxr
       DO  j = nys, nyn
          nzb_u_outer(j,i) = MAX( nzb_tmp(j-1,i), nzb_tmp(j,i),             &
                                  nzb_tmp(j+1,i) )
       ENDDO
!
!--    non-cyclic boundary conditions (overwritten by call of
!--    exchange_horiz_2d_int below in case of cyclic boundary conditions)
       IF ( nys == 0 )  THEN
          j = -1
          nzb_u_outer(j,i) = MAX( nzb_tmp(j+1,i), nzb_tmp(j,i) )
       ENDIF
       IF ( nyn == ny )  THEN
          j = ny + 1
          nzb_u_outer(j,i) = MAX( nzb_tmp(j-1,i), nzb_tmp(j,i) )
       ENDIF
    ENDDO

!
!-- nzb_v_inner:
!-- extend nzb_local northwards only
    nzb_tmp = nzb_local
    DO  i = nxl, nxr
       DO  j = nys, nyn
          nzb_tmp(j,i) = MAX( nzb_local(j-1,i), nzb_local(j,i) )
       ENDDO
    ENDDO
       
    CALL exchange_horiz_2d_int( nzb_tmp, nys, nyn, nxl, nxr, nbgp )      
    nzb_v_inner = nzb_tmp

!
!-- nzb_v_outer: 
!-- extend current nzb_tmp (nzb_v_inner) right-/leftwards
    DO  j = nys, nyn
       DO  i = nxl, nxr
          nzb_v_outer(j,i) = MAX( nzb_tmp(j,i-1), nzb_tmp(j,i),                &
                                  nzb_tmp(j,i+1) )
       ENDDO
!
!--    non-cyclic boundary conditions (overwritten by call of
!--    exchange_horiz_2d_int below in case of cyclic boundary conditions)
       IF ( nxl == 0 )  THEN
          i = -1
          nzb_v_outer(j,i) = MAX( nzb_tmp(j,i+1), nzb_tmp(j,i) )
       ENDIF
       IF ( nxr == nx )  THEN
          i = nx + 1
          nzb_v_outer(j,i) = MAX( nzb_tmp(j,i-1), nzb_tmp(j,i) )
       ENDIF
    ENDDO

!
!-- Exchange of lateral boundary values (parallel computers) and cyclic
!-- boundary conditions, if applicable.
!-- Since nzb_s_inner and nzb_w_inner are derived directly from nzb_local
!-- they do not require exchange and are not included here.
    CALL exchange_horiz_2d_int( nzb_u_inner, nys, nyn, nxl, nxr, nbgp )
    CALL exchange_horiz_2d_int( nzb_u_outer, nys, nyn, nxl, nxr, nbgp )
    CALL exchange_horiz_2d_int( nzb_v_inner, nys, nyn, nxl, nxr, nbgp )
    CALL exchange_horiz_2d_int( nzb_v_outer, nys, nyn, nxl, nxr, nbgp )
    CALL exchange_horiz_2d_int( nzb_w_outer, nys, nyn, nxl, nxr, nbgp )
    CALL exchange_horiz_2d_int( nzb_s_outer, nys, nyn, nxl, nxr, nbgp )

!
!-- Vertical nesting: communicate vertical grid level arrays between fine and
!-- coarse grid
    IF ( vnested )  CALL vnest_init_grid

 END SUBROUTINE init_grid


! Description:
! -----------------------------------------------------------------------------!
!> Calculation of the stretching factor through an iterative method. Ideas were 
!> taken from the paper "Regional stretched grid generation and its application
!> to the NCAR RegCM (1999)". Normally, no analytic solution exists because the
!> system of equations has two variables (r,l) but four requirements 
!> (l=integer, r=[0,88;1,2], Eq(6), Eq(5) starting from index j=1) which
!> results into an overdetermined system. 
!------------------------------------------------------------------------------!
 SUBROUTINE calculate_stretching_factor( number_end )
 
    USE control_parameters,                                                    &
        ONLY:  dz, dz_stretch_factor, dz_stretch_factor_array,                 &   
               dz_stretch_level_end, dz_stretch_level_start, message_string
  
    USE kinds
    
    IMPLICIT NONE
    
    INTEGER(iwp) ::  iterations  !< number of iterations until stretch_factor_lower/upper_limit is reached  
    INTEGER(iwp) ::  l_rounded   !< after l_rounded grid levels dz(n) is strechted to dz(n+1) with stretch_factor_2 
    INTEGER(iwp) ::  n           !< loop variable for stretching
    
    INTEGER(iwp), INTENT(IN) ::  number_end !< number of user-specified end levels for stretching
        
    REAL(wp) ::  delta_l               !< absolute difference between l and l_rounded
    REAL(wp) ::  delta_stretch_factor  !< absolute difference between stretch_factor_1 and stretch_factor_2
    REAL(wp) ::  delta_total_new       !< sum of delta_l and delta_stretch_factor for the next iteration (should be as small as possible) 
    REAL(wp) ::  delta_total_old       !< sum of delta_l and delta_stretch_factor for the last iteration 
    REAL(wp) ::  distance              !< distance between dz_stretch_level_start and dz_stretch_level_end (stretching region)
    REAL(wp) ::  l                     !< value that fulfil Eq. (5) in the paper mentioned above together with stretch_factor_1 exactly
    REAL(wp) ::  numerator             !< numerator of the quotient
    REAL(wp) ::  stretch_factor_1      !< stretching factor that fulfil Eq. (5) togehter with l exactly
    REAL(wp) ::  stretch_factor_2      !< stretching factor that fulfil Eq. (6) togehter with l_rounded exactly
    
    REAL(wp) ::  dz_stretch_factor_array_2(9) = 1.08_wp  !< Array that contains all stretch_factor_2 that belongs to stretch_factor_1 
    
    REAL(wp), PARAMETER ::  stretch_factor_interval = 1.0E-06  !< interval for sampling possible stretching factors
    REAL(wp), PARAMETER ::  stretch_factor_lower_limit = 0.88  !< lowest possible stretching factor
    REAL(wp), PARAMETER ::  stretch_factor_upper_limit = 1.12  !< highest possible stretching factor
 
 
    l = 0
    DO  n = 1, number_end
    
       iterations = 1
       stretch_factor_1 = 1.0 
       stretch_factor_2 = 1.0
       delta_total_old = 1.0
       
       IF ( dz(n) > dz(n+1) ) THEN
          DO WHILE ( stretch_factor_1 >= stretch_factor_lower_limit ) 
             
             stretch_factor_1 = 1.0 - iterations * stretch_factor_interval
             distance = ABS( dz_stretch_level_end(n) -                   &
                        dz_stretch_level_start(n) ) 
             numerator = distance*stretch_factor_1/dz(n) +               &
                         stretch_factor_1 - distance/dz(n)
             
             IF ( numerator > 0.0 ) THEN
                l = LOG( numerator ) / LOG( stretch_factor_1 ) - 1.0
                l_rounded = NINT( l )
                delta_l = ABS( l_rounded - l ) / l
             ENDIF
             
             stretch_factor_2 = EXP( LOG( dz(n+1)/dz(n) ) / (l_rounded) )
             
             delta_stretch_factor = ABS( stretch_factor_1 -              &
                                         stretch_factor_2 ) /            &
                                    stretch_factor_2
             
             delta_total_new = delta_l + delta_stretch_factor

!
!--                stretch_factor_1 is taken to guarantee that the stretching
!--                procedure ends as close as possible to dz_stretch_level_end.
!--                stretch_factor_2 would guarantee that the stretched dz(n) is
!--                equal to dz(n+1) after l_rounded grid levels.
             IF (delta_total_new < delta_total_old) THEN
                dz_stretch_factor_array(n) = stretch_factor_1
                dz_stretch_factor_array_2(n) = stretch_factor_2
                delta_total_old = delta_total_new
             ENDIF
             
             iterations = iterations + 1
            
          ENDDO
             
       ELSEIF ( dz(n) < dz(n+1) ) THEN 
          DO WHILE ( stretch_factor_1 <= stretch_factor_upper_limit )
                     
             stretch_factor_1 = 1.0 + iterations * stretch_factor_interval
             distance = ABS( dz_stretch_level_end(n) -                      &
                        dz_stretch_level_start(n) ) 
             numerator = distance*stretch_factor_1/dz(n) +                  &
                         stretch_factor_1 - distance/dz(n)
             
             l = LOG( numerator ) / LOG( stretch_factor_1 ) - 1.0
             l_rounded = NINT( l )
             delta_l = ABS( l_rounded - l ) / l
             
             stretch_factor_2 = EXP( LOG( dz(n+1)/dz(n) ) / (l_rounded) )

             delta_stretch_factor = ABS( stretch_factor_1 -                 &
                                        stretch_factor_2 ) /                &
                                        stretch_factor_2
             
             delta_total_new = delta_l + delta_stretch_factor
             
!
!--                stretch_factor_1 is taken to guarantee that the stretching
!--                procedure ends as close as possible to dz_stretch_level_end.
!--                stretch_factor_2 would guarantee that the stretched dz(n) is
!--                equal to dz(n+1) after l_rounded grid levels.
             IF (delta_total_new < delta_total_old) THEN
                dz_stretch_factor_array(n) = stretch_factor_1
                dz_stretch_factor_array_2(n) = stretch_factor_2
                delta_total_old = delta_total_new
             ENDIF
             
             iterations = iterations + 1
          ENDDO
          
       ELSE
          message_string= 'Two adjacent values of dz must be different'
          CALL message( 'init_grid', 'PA0228', 1, 2, 0, 6, 0 )
          
       ENDIF

!
!--    Check if also the second stretching factor fits into the allowed
!--    interval. If not, print a warning for the user.
       IF ( dz_stretch_factor_array_2(n) < stretch_factor_lower_limit .OR.     & 
            dz_stretch_factor_array_2(n) > stretch_factor_upper_limit ) THEN
          WRITE( message_string, * ) 'stretch_factor_2 = ',                    &
                                     dz_stretch_factor_array_2(n), ' which is',&
                                     ' responsible for exactly reaching& dz =',&
                                      dz(n+1), 'after a specific amount of',   & 
                                     ' grid levels& exceeds the upper',        &
                                     ' limit =', stretch_factor_upper_limit,   &
                                     ' &or lower limit = ',                    &
                                     stretch_factor_lower_limit
          CALL message( 'init_grid', 'PA0499', 0, 1, 0, 6, 0 )
            
       ENDIF
    ENDDO
        
 END SUBROUTINE calculate_stretching_factor
 
 
! Description:
! -----------------------------------------------------------------------------!
!> Set temporary topography flags and reference buildings on top of underlying
!> orography. 
!------------------------------------------------------------------------------!
 SUBROUTINE process_topography( topo_3d )

    USE arrays_3d,                                                             &
        ONLY:  zu, zw

    USE control_parameters,                                                    &
        ONLY:  bc_lr_cyc, bc_ns_cyc, land_surface, ocean, urban_surface

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg, nzb,  &
               nzt

    USE netcdf_data_input_mod,                                                 &
        ONLY:  buildings_f, building_id_f, input_pids_static,                  &
               terrain_height_f

    USE kinds

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i                !< running index along x-direction
    INTEGER(iwp) ::  j                !< running index along y-direction
    INTEGER(iwp) ::  k                !< running index along z-direction with respect to numeric grid
    INTEGER(iwp) ::  k2               !< running index along z-direction with respect to netcdf grid
    INTEGER(iwp) ::  nr               !< index variable indication maximum terrain height for respective building ID
    INTEGER(iwp) ::  num_build        !< counter for number of buildings
    INTEGER(iwp) ::  topo_top_index   !< orography top index, used to map 3D buildings onto terrain

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  displace_dum        !< displacements of start addresses, used for MPI_ALLGATHERV
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_ids           !< building IDs on entire model domain
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_ids_final     !< building IDs on entire model domain, multiple occurences are sorted out 
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_ids_final_tmp !< temporary array used for resizing
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_ids_l         !< building IDs on local subdomain
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_ids_l_tmp     !< temporary array used to resize array of building IDs

    INTEGER(iwp), DIMENSION(0:numprocs-1) ::  num_buildings     !< number of buildings with different ID on entire model domain
    INTEGER(iwp), DIMENSION(0:numprocs-1) ::  num_buildings_l   !< number of buildings with different ID on local subdomain

    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  topo_3d !< input array for 3D topography and dummy array for setting "outer"-flags

    REAL(wp)                            ::  ocean_offset        !< offset to consider inverse vertical coordinate at topography definition
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  oro_max             !< maximum terrain height occupied by an building with certain id
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  oro_max_l           !< maximum terrain height occupied by an building with certain id, on local subdomain

!
!-- In the following, buildings and orography are further preprocessed 
!-- before they are mapped on the LES grid.
!-- Buildings are mapped on top of the orography by maintaining the roof 
!-- shape of the building. This can be achieved by referencing building on 
!-- top of the maximum terrain height within the area occupied by the 
!-- respective building. As buildings and terrain height are defined PE-wise,
!-- parallelization of this referencing is required (a building can be 
!-- distributed between different PEs).  
!-- In a first step, determine the number of buildings with different 
!-- building id on each PE. In a next step, all building ids are gathered
!-- into one array which is present to all PEs. For each building ID, 
!-- the maximum terrain height occupied by the respective building is 
!-- computed and distributed to each PE.  
!-- Finally, for each building id and its respective reference orography, 
!-- builidings are mapped on top.   
!-- 
!-- First, pre-set topography flags, bit 1 indicates orography, bit 2 
!-- buildings 
!-- classify the respective surfaces.
    topo_3d          = IBSET( topo_3d, 0 )
    topo_3d(nzb,:,:) = IBCLR( topo_3d(nzb,:,:), 0 )
!
!-- In order to map topography on PALM grid also in case of ocean simulations,
!-- pre-calculate an offset value.
    ocean_offset = MERGE( zw(0), 0.0_wp, ocean )
!
!-- Reference buildings on top of orography. This is not necessary
!-- if topography is read from ASCII file as no distinction between buildings
!-- and terrain height can be made. Moreover, this is also not necessary if
!-- urban-surface and land-surface model are used at the same time. 
    IF ( input_pids_static )  THEN

       IF ( buildings_f%from_file )  THEN 
          num_buildings_l = 0
          num_buildings   = 0
!
!--       Allocate at least one element for building ids, 
          ALLOCATE( build_ids_l(1) )
          DO  i = nxl, nxr
             DO  j = nys, nyn
                IF ( building_id_f%var(j,i) /= building_id_f%fill )  THEN
                   IF ( num_buildings_l(myid) > 0 )  THEN
                      IF ( ANY( building_id_f%var(j,i) .EQ.  build_ids_l ) )   &
                      THEN
                         CYCLE
                      ELSE
                         num_buildings_l(myid) = num_buildings_l(myid) + 1
!
!--                   Resize array with different local building ids
                      ALLOCATE( build_ids_l_tmp(1:SIZE(build_ids_l)) )
                      build_ids_l_tmp = build_ids_l
                      DEALLOCATE( build_ids_l )
                      ALLOCATE( build_ids_l(1:num_buildings_l(myid)) )
                      build_ids_l(1:num_buildings_l(myid)-1) =                 &
                                  build_ids_l_tmp(1:num_buildings_l(myid)-1)
                      build_ids_l(num_buildings_l(myid)) = building_id_f%var(j,i)
                      DEALLOCATE( build_ids_l_tmp )
                   ENDIF
!
!--                First occuring building id on PE 
                   ELSE 
                      num_buildings_l(myid) = num_buildings_l(myid) + 1
                      build_ids_l(1) = building_id_f%var(j,i)
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
!
!--       Determine number of different building ids for the entire domain 
#if defined( __parallel ) 
          CALL MPI_ALLREDUCE( num_buildings_l, num_buildings, numprocs,              &
                              MPI_INTEGER, MPI_SUM, comm2d, ierr ) 
#else
          num_buildings = num_buildings_l
#endif
!
!--       Gather all buildings ids on each PEs. 
!--       First, allocate array encompassing all building ids in model domain.  
          ALLOCATE( build_ids(1:SUM(num_buildings)) )
#if defined( __parallel ) 
!
!--       Allocate array for displacements. 
!--       As each PE may has a different number of buildings, so that
!--       the block sizes send by each PE may not be equal. Hence, 
!--       information about the respective displacement is required, indicating 
!--       the respective adress where each MPI-task writes into the receive 
!--       buffer array  
          ALLOCATE( displace_dum(0:numprocs-1) )
          displace_dum(0) = 0
          DO i = 1, numprocs-1
             displace_dum(i) = displace_dum(i-1) + num_buildings(i-1)
          ENDDO

          CALL MPI_ALLGATHERV( build_ids_l(1:num_buildings_l(myid)),                 &
                               num_buildings(myid),                                  &
                               MPI_INTEGER,                                          &
                               build_ids,                                            &
                               num_buildings,                                        &
                               displace_dum,                                         & 
                               MPI_INTEGER,                                          &
                               comm2d, ierr )   

          DEALLOCATE( displace_dum )

#else
          build_ids = build_ids_l
#endif

!
!--       Note, in parallel mode building ids can occure mutliple times, as 
!--       each PE has send its own ids. Therefore, sort out building ids which 
!--       appear more than one time. 
          num_build = 0
          DO  nr = 1, SIZE(build_ids)

             IF ( ALLOCATED(build_ids_final) )  THEN
                IF ( ANY( build_ids(nr) .EQ. build_ids_final ) )  THEN
                   CYCLE
                ELSE
                   num_build = num_build + 1
!
!--                Resize
                   ALLOCATE( build_ids_final_tmp(1:num_build) )
                   build_ids_final_tmp(1:num_build-1) = build_ids_final(1:num_build-1)
                   DEALLOCATE( build_ids_final )
                   ALLOCATE( build_ids_final(1:num_build) )
                   build_ids_final(1:num_build-1) = build_ids_final_tmp(1:num_build-1)
                   build_ids_final(num_build) = build_ids(nr)
                   DEALLOCATE( build_ids_final_tmp )
                ENDIF             
             ELSE
                num_build = num_build + 1
                ALLOCATE( build_ids_final(1:num_build) )
                build_ids_final(num_build) = build_ids(nr)
             ENDIF
          ENDDO

!
!--       Determine maximumum terrain height occupied by the respective 
!--       building and temporalily store on oro_max 
          ALLOCATE( oro_max_l(1:SIZE(build_ids_final)) )
          ALLOCATE( oro_max(1:SIZE(build_ids_final))   )
          oro_max_l = 0.0_wp

          DO  nr = 1, SIZE(build_ids_final)
             oro_max_l(nr) = MAXVAL(                                              &
                              MERGE( terrain_height_f%var, 0.0_wp,                &
                                     building_id_f%var(nys:nyn,nxl:nxr) .EQ.      &
                                     build_ids_final(nr) ) )
          ENDDO
   
#if defined( __parallel )    
          IF ( SIZE(build_ids_final) >= 1 ) THEN
             CALL MPI_ALLREDUCE( oro_max_l, oro_max, SIZE( oro_max ), MPI_REAL,   &
                                 MPI_MAX, comm2d, ierr ) 
          ENDIF
#else
          oro_max = oro_max_l
#endif
!
!--       Finally, determine discrete grid height of maximum orography occupied
!--       by a building. Use all-or-nothing approach, i.e. a grid box is either
          oro_max_l = 0.0
          DO  nr = 1, SIZE(build_ids_final)
             DO  k = nzb, nzt
                IF ( zu(k) - ocean_offset <= oro_max(nr) )                     &
                   oro_max_l = zw(k) - ocean_offset
             ENDDO
             oro_max = oro_max_l
          ENDDO
       ENDIF
!
!--    Map orography as well as buildings onto grid. 
       DO  i = nxl, nxr
          DO  j = nys, nyn
             topo_top_index = 0
             DO  k = nzb, nzt
!
!--             In a first step, if grid point is below or equal the given 
!--             terrain height, grid point is flagged to be of type natural. 
!--             Please note, in case there is also a building which is lower
!--             than the vertical grid spacing, initialization of surface
!--             attributes will not be correct as given surface information
!--             will not be in accordance to the classified grid points. 
!--             Hence, in this case, de-flag the grid point and give it 
!--             urban type instead. 
                IF ( zu(k) - ocean_offset <= terrain_height_f%var(j,i) )  THEN
                    topo_3d(k,j,i) = IBCLR( topo_3d(k,j,i), 0 )
                    topo_3d(k,j,i) = IBSET( topo_3d(k,j,i), 1 )
                    topo_top_index = k ! topo_top_index + 1
                ENDIF
!
!--             Set building grid points. Here, only consider 2D buildings. 
!--             3D buildings require separate treatment. 
                IF ( buildings_f%from_file  .AND.  buildings_f%lod == 1 )  THEN
                   IF ( building_id_f%var(j,i) /= building_id_f%fill )  THEN
!
!--                   Determine index where maximum terrain height occupied by 
!--                   the respective building height is stored.
                      nr = MINLOC( ABS( build_ids_final -                      &
                                        building_id_f%var(j,i) ), DIM = 1 )
        
                      IF ( zu(k) - ocean_offset <=                             &
                           oro_max(nr) + buildings_f%var_2d(j,i) )  THEN
                         topo_3d(k,j,i) = IBCLR( topo_3d(k,j,i), 0 )
                         topo_3d(k,j,i) = IBSET( topo_3d(k,j,i), 2 )
!
!--                      De-flag grid point of type natural. See comment above.
                         topo_3d(k,j,i) = IBCLR( topo_3d(k,j,i), 1 ) 
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
!
!--          Map 3D buildings onto terrain height.  
!--          In case of any slopes, map building on top of maximum terrain 
!--          height covered by the building. In other words, extend 
!--          building down to the respective local terrain-surface height. 
             IF ( buildings_f%from_file  .AND.  buildings_f%lod == 2 )  THEN
                IF ( building_id_f%var(j,i) /= building_id_f%fill )  THEN
!
!--                Determine index for maximum-terrain-height array.
                   nr = MINLOC( ABS( build_ids_final -                         &
                                     building_id_f%var(j,i) ), DIM = 1 )
!
!--                Extend building down to the terrain surface, i.e. fill-up
!--                surface irregularities below a building. Note, oro_max
!--                is already a discrete height according to the all-or-nothing
!--                approach, i.e. grid box is either topography or atmosphere, 
!--                terrain top is defined at upper bound of the grid box.
!--                Hence, check for zw in this case. 
                   DO k = topo_top_index + 1, nzt + 1      
                      IF ( zw(k) - ocean_offset <= oro_max(nr) )  THEN
                         topo_3d(k,j,i) = IBCLR( topo_3d(k,j,i), 0 )
                         topo_3d(k,j,i) = IBSET( topo_3d(k,j,i), 2 )
                      ENDIF
                   ENDDO        
!
!--                After surface irregularities are smoothen, determine lower
!--                start index where building starts. 
                   DO  k = nzb, nzt
                      IF ( zw(k) - ocean_offset <= oro_max(nr) )               &
                         topo_top_index = k
                   ENDDO
!
!--                Finally, map building on top.
                   k2 = 0
                   DO k = topo_top_index, nzt + 1
                      IF ( k2 <= buildings_f%nz-1 )  THEN
                         IF ( buildings_f%var_3d(k2,j,i) == 1 )  THEN
                            topo_3d(k,j,i) = IBCLR( topo_3d(k,j,i), 0 )
                            topo_3d(k,j,i) = IBCLR( topo_3d(k,j,i), 1 )
                            topo_3d(k,j,i) = IBSET( topo_3d(k,j,i), 2 )
                         ENDIF
                      ENDIF
                      k2 = k2 + 1
                   ENDDO
                ENDIF
             ENDIF
          ENDDO
       ENDDO
!
!--    Deallocate temporary arrays required for processing and reading data
       IF ( ALLOCATED( oro_max         ) )  DEALLOCATE( oro_max         )
       IF ( ALLOCATED( oro_max_l       ) )  DEALLOCATE( oro_max_l       )
       IF ( ALLOCATED( build_ids_final ) )  DEALLOCATE( build_ids_final )
!
!-- Topography input via ASCII format. 
    ELSE
       ocean_offset     = MERGE( zw(0), 0.0_wp, ocean )
       topo_3d          = IBSET( topo_3d, 0 )
       topo_3d(nzb,:,:) = IBCLR( topo_3d(nzb,:,:), 0 )
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb, nzt
                IF ( zu(k) - ocean_offset <= buildings_f%var_2d(j,i) )  THEN
                    topo_3d(k,j,i) = IBCLR( topo_3d(k,j,i), 0 )
                    topo_3d(k,j,i) = IBSET( topo_3d(k,j,i), 2 ) !indicates terrain
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    CALL exchange_horiz_int( topo_3d, nys, nyn, nxl, nxr, nzt, nbgp )

    IF ( .NOT. bc_ns_cyc )  THEN
       IF ( nys == 0  )  topo_3d(:,-1,:)   = topo_3d(:,0,:)
       IF ( nyn == ny )  topo_3d(:,ny+1,:) = topo_3d(:,ny,:)
    ENDIF

    IF ( .NOT. bc_lr_cyc )  THEN
       IF ( nxl == 0  )  topo_3d(:,:,-1)   = topo_3d(:,:,0)
       IF ( nxr == nx )  topo_3d(:,:,nx+1) = topo_3d(:,:,nx)          
    ENDIF

 END SUBROUTINE process_topography


! Description:
! -----------------------------------------------------------------------------!
!> Filter topography, i.e. fill holes resolved by only one grid point.  
!> Such holes are suspected to lead to velocity blow-ups as continuity 
!> equation on discrete grid cannot be fulfilled in such case.
!------------------------------------------------------------------------------!
 SUBROUTINE filter_topography( topo_3d )

    USE control_parameters,                                                    &
        ONLY:  bc_lr_cyc, bc_ns_cyc, message_string

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg, nzb, nzt

    USE netcdf_data_input_mod,                                                 &
        ONLY:  building_id_f, building_type_f 

    USE  pegrid

    IMPLICIT NONE

    LOGICAL      ::  filled = .FALSE. !< flag indicating if holes were filled

    INTEGER(iwp) ::  i          !< running index along x-direction
    INTEGER(iwp) ::  j          !< running index along y-direction
    INTEGER(iwp) ::  k          !< running index along z-direction
    INTEGER(iwp) ::  num_hole   !< number of holes (in topography) resolved by only one grid point 
    INTEGER(iwp) ::  num_hole_l !< number of holes (in topography) resolved by only one grid point on local PE     
    INTEGER(iwp) ::  num_wall   !< number of surrounding vertical walls for a single grid point

    INTEGER(iwp), DIMENSION(nysg:nyng,nxlg:nxrg)           ::  var_exchange_int  !< dummy array for exchanging ghost-points 
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE            ::  topo_tmp          !< temporary 3D-topography used to fill holes
    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  topo_3d           !< 3D-topography array merging buildings and orography
!
!-- Before checking for holes, set lateral boundary conditions for 
!-- topography. After hole-filling, boundary conditions must be set again.
!-- Several iterations are performed, in order to fill holes which might 
!-- emerge by the filling-algorithm itself.
    ALLOCATE( topo_tmp(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    topo_tmp = 0

    num_hole = 99999
    DO WHILE ( num_hole > 0 )       

       num_hole = 0    
       CALL exchange_horiz_int( topo_3d, nys, nyn, nxl, nxr, nzt, nbgp )
!
!--    Exchange also building ID and type. Note, building_type is an one-byte 
!--    variable.
       IF ( building_id_f%from_file )                                          &
          CALL exchange_horiz_2d_int( building_id_f%var, nys, nyn, nxl, nxr, nbgp )
       IF ( building_type_f%from_file )  THEN
          var_exchange_int = INT( building_type_f%var, KIND = 4 )
          CALL exchange_horiz_2d_int( var_exchange_int, nys, nyn, nxl, nxr, nbgp )
          building_type_f%var = INT( var_exchange_int, KIND = 1 )
       ENDIF

       topo_tmp = topo_3d
!
!--    In case of non-cyclic lateral boundaries, assume lateral boundary to be 
!--    a solid wall. Thus, intermediate spaces of one grid point between 
!--    boundary and some topographic structure will be filled.           
       IF ( .NOT. bc_ns_cyc )  THEN
          IF ( nys == 0  )  topo_tmp(:,-1,:)   = IBCLR( topo_tmp(:,0,:),  0 )
          IF ( nyn == ny )  topo_tmp(:,ny+1,:) = IBCLR( topo_tmp(:,ny,:), 0 )
       ENDIF

       IF ( .NOT. bc_lr_cyc )  THEN
          IF ( nxl == 0  )  topo_tmp(:,:,-1)   = IBCLR( topo_tmp(:,:,0),  0 )
          IF ( nxr == nx )  topo_tmp(:,:,nx+1) = IBCLR( topo_tmp(:,:,nx), 0 )          
       ENDIF

       num_hole_l = 0
       DO i = nxl, nxr
          DO j = nys, nyn
             DO  k = nzb+1, nzt
                IF ( BTEST( topo_tmp(k,j,i), 0 ) )  THEN
                   num_wall = 0
                   IF ( .NOT. BTEST( topo_tmp(k,j-1,i), 0 ) )                  &
                      num_wall = num_wall + 1
                   IF ( .NOT. BTEST( topo_tmp(k,j+1,i), 0 ) )                  &
                      num_wall = num_wall + 1
                   IF ( .NOT. BTEST( topo_tmp(k,j,i-1), 0 ) )                  &
                      num_wall = num_wall + 1
                   IF ( .NOT. BTEST( topo_tmp(k,j,i+1), 0 ) )                  &
                      num_wall = num_wall + 1
                   IF ( .NOT. BTEST( topo_tmp(k-1,j,i), 0 ) )                  &
                      num_wall = num_wall + 1   
                   IF ( .NOT. BTEST( topo_tmp(k+1,j,i), 0 ) )                  &
                      num_wall = num_wall + 1

                   IF ( num_wall >= 4 )  THEN
                      num_hole_l     = num_hole_l + 1
!
!--                   Clear flag 0 and set special flag ( bit 3) to indicate 
!--                   that new topography point is a result of filtering process.
                      topo_3d(k,j,i) = IBCLR( topo_3d(k,j,i), 0 )
                      topo_3d(k,j,i) = IBSET( topo_3d(k,j,i), 3 )
!
!--                   If filled grid point is occupied by a building, classify
!--                   it as building grid point.
                      IF ( building_type_f%from_file )  THEN
                         IF ( building_type_f%var(j,i)   /=                    &  
                              building_type_f%fill            .OR.             &       
                              building_type_f%var(j+1,i) /=                    &  
                              building_type_f%fill            .OR.             &                
                              building_type_f%var(j-1,i) /=                    &                
                              building_type_f%fill            .OR.             &                
                              building_type_f%var(j,i+1) /=                    &                
                              building_type_f%fill            .OR.             &                
                              building_type_f%var(j,i-1) /=                    &                
                              building_type_f%fill )  THEN
!
!--                         Set flag indicating building surfaces
                            topo_3d(k,j,i) = IBSET( topo_3d(k,j,i), 2 )
!
!--                         Set building_type and ID at this position if not 
!--                         already set. This is required for proper 
!--                         initialization of urban-surface energy balance 
!--                         solver.
                            IF ( building_type_f%var(j,i) ==                   &
                                 building_type_f%fill )  THEN

                               IF ( building_type_f%var(j+1,i) /=              &
                                    building_type_f%fill )  THEN
                                  building_type_f%var(j,i) =                   &
                                                    building_type_f%var(j+1,i)
                                  building_id_f%var(j,i) =                     &
                                                    building_id_f%var(j+1,i)
                               ELSEIF ( building_type_f%var(j-1,i) /=          &
                                        building_type_f%fill )  THEN
                                  building_type_f%var(j,i) =                   &
                                                    building_type_f%var(j-1,i)
                                  building_id_f%var(j,i) =                     &
                                                    building_id_f%var(j-1,i)
                               ELSEIF ( building_type_f%var(j,i+1) /=          &
                                        building_type_f%fill )  THEN
                                  building_type_f%var(j,i) =                   &
                                                    building_type_f%var(j,i+1)
                                  building_id_f%var(j,i) =                     &
                                                    building_id_f%var(j,i+1)
                               ELSEIF ( building_type_f%var(j,i-1) /=          &
                                        building_type_f%fill )  THEN
                                  building_type_f%var(j,i) =                   &
                                                    building_type_f%var(j,i-1)
                                  building_id_f%var(j,i) =                     &
                                                    building_id_f%var(j,i-1)
                               ENDIF
                            ENDIF
                         ENDIF
                      ENDIF
!
!--                   If filled grid point is already classified as building
!--                   everything is fine, else classify this grid point as
!--                   natural type grid point. This case, values for the 
!--                   surface type are already set.
                      IF ( .NOT. BTEST( topo_3d(k,j,i), 2 ) )  THEN
                         topo_3d(k,j,i) = IBSET( topo_3d(k,j,i), 1 )
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
!
!--    Count the total number of holes, required for informative message.
#if defined( __parallel )
       CALL MPI_ALLREDUCE( num_hole_l, num_hole, 1, MPI_INTEGER, MPI_SUM,      &
                           comm2d, ierr )
#else
       num_hole = num_hole_l
#endif    
       IF ( num_hole > 0  .AND.  .NOT. filled )  filled = .TRUE.

    ENDDO
!
!-- Create an informative message if any holes were filled.
    IF ( filled )  THEN
       WRITE( message_string, * ) 'Topography was filtered, i.e. holes ' //    &
                                  'resolved by only one grid point '     //    &
                                  'were filled during initialization.'
       CALL message( 'init_grid', 'PA0430', 0, 0, 0, 6, 0 )
    ENDIF

    DEALLOCATE( topo_tmp )
!
!-- Finally, exchange topo_3d array again and if necessary set Neumann boundary
!-- condition in case of non-cyclic lateral boundaries. 
    CALL exchange_horiz_int( topo_3d, nys, nyn, nxl, nxr, nzt, nbgp )

    IF ( .NOT. bc_ns_cyc )  THEN
       IF ( nys == 0  )  topo_3d(:,-1,:)   = topo_3d(:,0,:)
       IF ( nyn == ny )  topo_3d(:,ny+1,:) = topo_3d(:,ny,:)
    ENDIF

    IF ( .NOT. bc_lr_cyc )  THEN
       IF ( nxl == 0  )  topo_3d(:,:,-1)   = topo_3d(:,:,0)
       IF ( nxr == nx )  topo_3d(:,:,nx+1) = topo_3d(:,:,nx)          
    ENDIF
!
!-- Exchange building ID and type. Note, building_type is an one-byte variable.
    IF ( building_id_f%from_file )                                             &
       CALL exchange_horiz_2d_int( building_id_f%var, nys, nyn, nxl, nxr, nbgp )
    IF ( building_type_f%from_file )  THEN
       var_exchange_int = INT( building_type_f%var, KIND = 4 )
       CALL exchange_horiz_2d_int( var_exchange_int, nys, nyn, nxl, nxr, nbgp )
       building_type_f%var = INT( var_exchange_int, KIND = 1 )
    ENDIF

 END SUBROUTINE filter_topography


! Description:
! -----------------------------------------------------------------------------!
!> Reads topography information from file or sets generic topography. Moreover, 
!> all topography-relevant topography arrays are initialized, and grid flags
!> are set.  
!------------------------------------------------------------------------------!
 SUBROUTINE init_topo( topo )

    USE arrays_3d,                                                             &
        ONLY:  zw
        
    USE control_parameters,                                                    &
        ONLY:  bc_lr_cyc, bc_ns_cyc, building_height, building_length_x,       &
               building_length_y, building_wall_left, building_wall_south,     &
               canyon_height, canyon_wall_left, canyon_wall_south,             &
               canyon_width_x, canyon_width_y, constant_flux_layer,            &
               dp_level_ind_b, dz,                                             &
               message_string, ocean, topography, topography_grid_convention,  &
               tunnel_height, tunnel_length, tunnel_width_x, tunnel_width_y,   &
               tunnel_wall_depth
         
    USE grid_variables,                                                        &
        ONLY:  dx, dy
        
    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg, nz,   &
               nzb, nzt
    
    USE kinds

    USE pegrid

    USE surface_mod,                                                           &
        ONLY:  get_topography_top_index, get_topography_top_index_ji

    IMPLICIT NONE

    INTEGER(iwp) ::  bh            !< temporary vertical index of building height
    INTEGER(iwp) ::  blx           !< grid point number of building size along x
    INTEGER(iwp) ::  bly           !< grid point number of building size along y
    INTEGER(iwp) ::  bxl           !< index for left building wall
    INTEGER(iwp) ::  bxr           !< index for right building wall
    INTEGER(iwp) ::  byn           !< index for north building wall
    INTEGER(iwp) ::  bys           !< index for south building wall
    INTEGER(iwp) ::  ch            !< temporary vertical index for canyon height
    INTEGER(iwp) ::  cwx           !< grid point number of canyon size along x
    INTEGER(iwp) ::  cwy           !< grid point number of canyon size along y
    INTEGER(iwp) ::  cxl           !< index for left canyon wall
    INTEGER(iwp) ::  cxr           !< index for right canyon wall
    INTEGER(iwp) ::  cyn           !< index for north canyon wall
    INTEGER(iwp) ::  cys           !< index for south canyon wall
    INTEGER(iwp) ::  i             !< index variable along x
    INTEGER(iwp) ::  j             !< index variable along y
    INTEGER(iwp) ::  k             !< index variable along z
    INTEGER(iwp) ::  hv_in         !< heavyside function to model inner tunnel surface 
    INTEGER(iwp) ::  hv_out        !< heavyside function to model outer tunnel surface 
    INTEGER(iwp) ::  txe_out       !< end position of outer tunnel wall in x
    INTEGER(iwp) ::  txs_out       !< start position of outer tunnel wall in x
    INTEGER(iwp) ::  tye_out       !< end position of outer tunnel wall in y
    INTEGER(iwp) ::  tys_out       !< start position of outer tunnel wall in y
    INTEGER(iwp) ::  txe_in        !< end position of inner tunnel wall in x
    INTEGER(iwp) ::  txs_in        !< start position of inner tunnel wall in x
    INTEGER(iwp) ::  tye_in        !< end position of inner tunnel wall in y
    INTEGER(iwp) ::  tys_in        !< start position of inner tunnel wall in y
    INTEGER(iwp) ::  td            !< tunnel wall depth
    INTEGER(iwp) ::  th            !< height of outer tunnel wall

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  nzb_local         !< index for topography top at cell-center
    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  topo !< input array for 3D topography and dummy array for setting "outer"-flags


!
!-- Set outer and inner index arrays for non-flat topography.
!-- Here consistency checks concerning domain size and periodicity are
!-- necessary.
!-- Within this SELECT CASE structure only nzb_local is initialized
!-- individually depending on the chosen topography type, all other index 
!-- arrays are initialized further below.
    SELECT CASE ( TRIM( topography ) )

       CASE ( 'flat' )
!   
!--       Initialilize 3D topography array, used later for initializing flags
          IF ( TRIM(constant_flux_layer) == 'top' ) THEN
             topo(nzb+1:nzt,:,:) = IBSET( topo(nzb+1:nzt,:,:), 0 ) 
          ELSE
             topo(nzb+1:nzt+1,:,:) = IBSET( topo(nzb+1:nzt+1,:,:), 0 ) 
          ENDIF

       CASE ( 'single_building' )
!
!--       Single rectangular building, by default centered in the middle of the
!--       total domain
          blx = NINT( building_length_x / dx )
          bly = NINT( building_length_y / dy )
          bh  = MINLOC( ABS( zw - building_height ), 1 ) - 1
          IF ( ABS( zw(bh)   - building_height ) == &
               ABS( zw(bh+1) - building_height )    )  bh = bh + 1
          IF ( building_wall_left == 9999999.9_wp )  THEN
             building_wall_left = ( nx + 1 - blx ) / 2 * dx
          ENDIF
          bxl = NINT( building_wall_left / dx )
          bxr = bxl + blx

          IF ( building_wall_south == 9999999.9_wp )  THEN
              building_wall_south = ( ny + 1 - bly ) / 2 * dy
          ENDIF
          bys = NINT( building_wall_south / dy )
          byn = bys + bly

!
!--       Building size has to meet some requirements
          IF ( ( bxl < 1 ) .OR. ( bxr > nx-1 ) .OR. ( bxr < bxl+3 ) .OR.       &
               ( bys < 1 ) .OR. ( byn > ny-1 ) .OR. ( byn < bys+3 ) )  THEN
             WRITE( message_string, * ) 'inconsistent building parameters:',   &
                                      '&bxl=', bxl, 'bxr=', bxr, 'bys=', bys,  &
                                      'byn=', byn, 'nx=', nx, 'ny=', ny
             CALL message( 'init_grid', 'PA0203', 1, 2, 0, 6, 0 )
          ENDIF

          ALLOCATE( nzb_local(nysg:nyng,nxlg:nxrg) )
          nzb_local = 0
!
!--       Define the building. 
          IF ( bxl <= nxr  .AND.  bxr >= nxl  .AND.                            &
               bys <= nyn  .AND.  byn >= nys )                                 & 
             nzb_local(MAX(nys,bys):MIN(nyn,byn),MAX(nxl,bxl):MIN(nxr,bxr)) = bh
!
!--       Set bit array on basis of nzb_local
          DO  i = nxl, nxr
             DO  j = nys, nyn
                topo(nzb_local(j,i)+1:nzt+1,j,i) =                             &
                                 IBSET( topo(nzb_local(j,i)+1:nzt+1,j,i), 0 ) 
             ENDDO
          ENDDO
        
          DEALLOCATE( nzb_local )

          CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )
!
!--       Set boundary conditions also for flags. Can be interpreted as Neumann
!--       boundary conditions for topography. 
          IF ( .NOT. bc_ns_cyc )  THEN
             IF ( nys == 0  )  THEN
                DO  i = 1, nbgp     
                   topo(:,nys-i,:)   = topo(:,nys,:)
                ENDDO
             ENDIF
             IF ( nyn == ny )  THEN
                DO  i = 1, nbgp  
                   topo(:,nyn+i,:) = topo(:,nyn,:)
                ENDDO
             ENDIF
          ENDIF
          IF ( .NOT. bc_lr_cyc )  THEN
             IF ( nxl == 0  )  THEN
                DO  i = 1, nbgp   
                   topo(:,:,nxl-i)   = topo(:,:,nxl)
                ENDDO
             ENDIF
             IF ( nxr == nx )  THEN 
                DO  i = 1, nbgp   
                   topo(:,:,nxr+i) = topo(:,:,nxr)     
                ENDDO
             ENDIF      
          ENDIF

       CASE ( 'single_street_canyon' )
!
!--       Single quasi-2D street canyon of infinite length in x or y direction.
!--       The canyon is centered in the other direction by default.
          IF ( canyon_width_x /= 9999999.9_wp )  THEN
!
!--          Street canyon in y direction
             cwx = NINT( canyon_width_x / dx )
             IF ( canyon_wall_left == 9999999.9_wp )  THEN
                canyon_wall_left = ( nx + 1 - cwx ) / 2 * dx
             ENDIF
             cxl = NINT( canyon_wall_left / dx )
             cxr = cxl + cwx
          ELSEIF ( canyon_width_y /= 9999999.9_wp )  THEN
!
!--          Street canyon in x direction
             cwy = NINT( canyon_width_y / dy )
             IF ( canyon_wall_south == 9999999.9_wp )  THEN
                canyon_wall_south = ( ny + 1 - cwy ) / 2 * dy
             ENDIF
             cys = NINT( canyon_wall_south / dy )
             cyn = cys + cwy
      
          ELSE
             
             message_string = 'no street canyon width given'
             CALL message( 'init_grid', 'PA0204', 1, 2, 0, 6, 0 )
  
          ENDIF

          ch  = MINLOC( ABS( zw - canyon_height ), 1 ) - 1
          IF ( ABS( zw(ch)   - canyon_height ) == &
               ABS( zw(ch+1) - canyon_height )    )  ch = ch + 1
          dp_level_ind_b = ch
!
!--       Street canyon size has to meet some requirements
          IF ( canyon_width_x /= 9999999.9_wp )  THEN
             IF ( ( cxl < 1 ) .OR. ( cxr > nx-1 ) .OR. ( cwx < 3 ) .OR.        &
                  ( ch < 3 ) )  THEN
                WRITE( message_string, * ) 'inconsistent canyon parameters:',  &
                                           '&cxl=', cxl, ' cxr=', cxr,         &
                                           ' cwx=', cwx,                       &
                                           ' ch=', ch, ' nx=', nx, ' ny=', ny
                CALL message( 'init_grid', 'PA0205', 1, 2, 0, 6, 0 ) 
             ENDIF
          ELSEIF ( canyon_width_y /= 9999999.9_wp )  THEN
             IF ( ( cys < 1 ) .OR. ( cyn > ny-1 ) .OR. ( cwy < 3 ) .OR.        &
                  ( ch < 3 ) )  THEN
                WRITE( message_string, * ) 'inconsistent canyon parameters:',  &
                                           '&cys=', cys, ' cyn=', cyn,         &
                                           ' cwy=', cwy,                       &
                                           ' ch=', ch, ' nx=', nx, ' ny=', ny
                CALL message( 'init_grid', 'PA0206', 1, 2, 0, 6, 0 ) 
             ENDIF
          ENDIF
          IF ( canyon_width_x /= 9999999.9_wp .AND.                            &                 
               canyon_width_y /= 9999999.9_wp )  THEN
             message_string = 'inconsistent canyon parameters:' //             &   
                              '&street canyon can only be oriented' //         &
                              ' either in x- or in y-direction'
             CALL message( 'init_grid', 'PA0207', 1, 2, 0, 6, 0 )
          ENDIF

          ALLOCATE( nzb_local(nysg:nyng,nxlg:nxrg) )
          nzb_local = ch
          IF ( canyon_width_x /= 9999999.9_wp )  THEN
             IF ( cxl <= nxr  .AND.  cxr >= nxl )                              &
                nzb_local(:,MAX(nxl,cxl+1):MIN(nxr,cxr-1)) = 0
          ELSEIF ( canyon_width_y /= 9999999.9_wp )  THEN
             IF ( cys <= nyn  .AND.  cyn >= nys )                              &          
                nzb_local(MAX(nys,cys+1):MIN(nyn,cyn-1),:) = 0
          ENDIF
!
!--       Set bit array on basis of nzb_local
          DO  i = nxl, nxr
             DO  j = nys, nyn
                topo(nzb_local(j,i)+1:nzt+1,j,i) =                             &
                                 IBSET( topo(nzb_local(j,i)+1:nzt+1,j,i), 0 ) 
             ENDDO
          ENDDO
          DEALLOCATE( nzb_local )

          CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )
!
!--       Set boundary conditions also for flags. Can be interpreted as Neumann
!--       boundary conditions for topography. 
          IF ( .NOT. bc_ns_cyc )  THEN
             IF ( nys == 0  )  THEN
                DO  i = 1, nbgp     
                   topo(:,nys-i,:)   = topo(:,nys,:)
                ENDDO
             ENDIF
             IF ( nyn == ny )  THEN
                DO  i = 1, nbgp  
                   topo(:,nyn+i,:) = topo(:,nyn,:)
                ENDDO
             ENDIF
          ENDIF
          IF ( .NOT. bc_lr_cyc )  THEN
             IF ( nxl == 0  )  THEN
                DO  i = 1, nbgp   
                   topo(:,:,nxl-i)   = topo(:,:,nxl)
                ENDDO
             ENDIF
             IF ( nxr == nx )  THEN 
                DO  i = 1, nbgp   
                   topo(:,:,nxr+i) = topo(:,:,nxr)     
                ENDDO
             ENDIF      
          ENDIF

       CASE ( 'tunnel' )

!
!--       Tunnel height
          IF ( tunnel_height == 9999999.9_wp )  THEN
             th = zw( INT( 0.2 * nz) )
          ELSE
             th = tunnel_height
          ENDIF
!
!--       Tunnel-wall depth
          IF ( tunnel_wall_depth == 9999999.9_wp )  THEN  
             td = MAX ( dx, dy, dz(1) )
          ELSE
             td = tunnel_wall_depth
          ENDIF
!
!--       Check for tunnel width
          IF ( tunnel_width_x == 9999999.9_wp  .AND.                           &
               tunnel_width_y == 9999999.9_wp  )  THEN
             message_string = 'No tunnel width is given. '
             CALL message( 'init_grid', 'PA0280', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( tunnel_width_x /= 9999999.9_wp  .AND.                           &
               tunnel_width_y /= 9999999.9_wp  )  THEN
             message_string = 'Inconsistent tunnel parameters:' //             &   
                              'tunnel can only be oriented' //                 &
                              'either in x- or in y-direction.'
             CALL message( 'init_grid', 'PA0281', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       Tunnel axis along y
          IF ( tunnel_width_x /= 9999999.9_wp )  THEN
             IF ( tunnel_width_x > ( nx + 1 ) * dx )  THEN
                message_string = 'Tunnel width too large'
                CALL message( 'init_grid', 'PA0282', 1, 2, 0, 6, 0 )
             ENDIF

             txs_out = INT( ( nx + 1 ) * 0.5_wp * dx - tunnel_width_x * 0.5_wp )
             txe_out = INT( ( nx + 1 ) * 0.5_wp * dx + tunnel_width_x * 0.5_wp )
             txs_in  = INT( ( nx + 1 ) * 0.5_wp * dx -                         &
                                      ( tunnel_width_x * 0.5_wp - td ) )
             txe_in  = INT( ( nx + 1 ) * 0.5_wp * dx +                         &
                                   ( tunnel_width_x * 0.5_wp - td ) )

             tys_out = INT( ( ny + 1 ) * 0.5_wp * dy - tunnel_length * 0.5_wp )
             tye_out = INT( ( ny + 1 ) * 0.5_wp * dy + tunnel_length * 0.5_wp )
             tys_in  = tys_out
             tye_in  = tye_out
          ENDIF
          IF ( tunnel_width_x /= 9999999.9_wp  .AND.                           &   
               tunnel_width_x - 2.0_wp * td <= 2.0_wp * dx )                   &
          THEN
             message_string = 'Tunnel width too small'
             CALL message( 'init_grid', 'PA0175', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( tunnel_width_y /= 9999999.9_wp  .AND.                           &
               tunnel_width_y - 2.0_wp * td <= 2.0_wp * dy )                   &
          THEN
             message_string = 'Tunnel width too small'
             CALL message( 'init_grid', 'PA0455', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       Tunnel axis along x
          IF ( tunnel_width_y /= 9999999.9_wp )  THEN
             IF ( tunnel_width_y > ( ny + 1 ) * dy )  THEN
                message_string = 'Tunnel width too large'
                CALL message( 'init_grid', 'PA0456', 1, 2, 0, 6, 0 )
             ENDIF

             txs_out = INT( ( nx + 1 ) * 0.5_wp * dx - tunnel_length * 0.5_wp )
             txe_out = INT( ( nx + 1 ) * 0.5_wp * dx + tunnel_length * 0.5_wp )
             txs_in  = txs_out
             txe_in  = txe_out

             tys_out = INT( ( ny + 1 ) * 0.5_wp * dy - tunnel_width_y * 0.5_wp )
             tye_out = INT( ( ny + 1 ) * 0.5_wp * dy + tunnel_width_y * 0.5_wp )
             tys_in  = INT( ( ny + 1 ) * 0.5_wp * dy -                         &
                                        ( tunnel_width_y * 0.5_wp - td ) )
             tye_in  = INT( ( ny + 1 ) * 0.5_wp * dy +                         &
                                     ( tunnel_width_y * 0.5_wp - td ) )
          ENDIF

          topo = 0
          DO  i = nxl, nxr
             DO  j = nys, nyn
!
!--             Use heaviside function to model outer tunnel surface
                hv_out = th * 0.5_wp *                                         &
                              ( ( SIGN( 1.0_wp, i * dx - txs_out ) + 1.0_wp )  &
                              - ( SIGN( 1.0_wp, i * dx - txe_out ) + 1.0_wp ) )

                hv_out = hv_out * 0.5_wp *                                     &
                            ( ( SIGN( 1.0_wp, j * dy - tys_out ) + 1.0_wp )    &
                            - ( SIGN( 1.0_wp, j * dy - tye_out ) + 1.0_wp ) )
!    
!--             Use heaviside function to model inner tunnel surface
                hv_in  = ( th - td ) * 0.5_wp *                                &
                                ( ( SIGN( 1.0_wp, i * dx - txs_in ) + 1.0_wp ) &
                                - ( SIGN( 1.0_wp, i * dx - txe_in ) + 1.0_wp ) )

                hv_in = hv_in * 0.5_wp *                                       &
                                ( ( SIGN( 1.0_wp, j * dy - tys_in ) + 1.0_wp ) &
                                - ( SIGN( 1.0_wp, j * dy - tye_in ) + 1.0_wp ) )
!
!--             Set flags at x-y-positions without any tunnel surface
                IF ( hv_out - hv_in == 0.0_wp )  THEN
                   topo(nzb+1:nzt+1,j,i) = IBSET( topo(nzb+1:nzt+1,j,i), 0 )
!
!--             Set flags at x-y-positions with tunnel surfaces
                ELSE
                   DO  k = nzb + 1, nzt + 1
!
!--                   Inner tunnel
                      IF ( hv_out - hv_in == th )  THEN
                         IF ( zw(k) <= hv_out )  THEN
                            topo(k,j,i) = IBCLR( topo(k,j,i), 0 )
                         ELSE
                            topo(k,j,i) = IBSET( topo(k,j,i), 0 )
                         ENDIF
                      ENDIF
!
!--                   Lateral tunnel walls
                      IF ( hv_out - hv_in == td )  THEN
                         IF ( zw(k) <= hv_in )  THEN 
                            topo(k,j,i) = IBSET( topo(k,j,i), 0 )
                         ELSEIF ( zw(k) > hv_in  .AND.  zw(k) <= hv_out )  THEN 
                            topo(k,j,i) = IBCLR( topo(k,j,i), 0 )
                         ELSEIF ( zw(k) > hv_out )  THEN 
                            topo(k,j,i) = IBSET( topo(k,j,i), 0 )
                         ENDIF
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO
          ENDDO

          CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )
!
!--       Set boundary conditions also for flags. Can be interpreted as Neumann
!--       boundary conditions for topography. 
          IF ( .NOT. bc_ns_cyc )  THEN
             IF ( nys == 0  )  THEN
                DO  i = 1, nbgp     
                   topo(:,nys-i,:)   = topo(:,nys,:)
                ENDDO
             ENDIF
             IF ( nyn == ny )  THEN
                DO  i = 1, nbgp  
                   topo(:,nyn+i,:) = topo(:,nyn,:)
                ENDDO
             ENDIF
          ENDIF
          IF ( .NOT. bc_lr_cyc )  THEN
             IF ( nxl == 0  )  THEN
                DO  i = 1, nbgp   
                   topo(:,:,nxl-i)   = topo(:,:,nxl)
                ENDDO
             ENDIF
             IF ( nxr == nx )  THEN 
                DO  i = 1, nbgp   
                   topo(:,:,nxr+i) = topo(:,:,nxr)     
                ENDDO
             ENDIF      
          ENDIF

       CASE ( 'read_from_file' )
!
!--       Note, topography information have been already read.  
!--       If required, further process topography, i.e. reference buildings on
!--       top of orography and set temporary 3D topography array, which is 
!--       used later to set grid flags. Calling of this rouinte is also 
!--       required in case of ASCII input, even though no distinction between 
!--       terrain- and building height is made in this case.  
          CALL process_topography( topo )
!
!--       Filter holes resolved by only one grid-point
          CALL filter_topography( topo )
!
!--       Exchange ghost-points, as well as add cyclic or Neumann boundary 
!--       conditions.
          CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )
!
!--       Set lateral boundary conditions for topography on all ghost layers          
          IF ( .NOT. bc_ns_cyc )  THEN
             IF ( nys == 0  )  THEN
                DO  i = 1, nbgp         
                   topo(:,nys-i,:) = topo(:,nys,:)
                ENDDO
             ENDIF
             IF ( nyn == ny )  THEN
                DO  i = 1, nbgp         
                   topo(:,nyn+i,:) = topo(:,nyn,:)
                ENDDO
             ENDIF
          ENDIF

          IF ( .NOT. bc_lr_cyc )  THEN
             IF ( nxl == 0  )  THEN
                DO  i = 1, nbgp 
                   topo(:,:,nxl-i) = topo(:,:,nxl)
                ENDDO
             ENDIF
             IF ( nxr == nx )  THEN
                DO  i = 1, nbgp 
                   topo(:,:,nxr+i) = topo(:,:,nxr)
                ENDDO
             ENDIF
          ENDIF


       CASE DEFAULT
!   
!--       The DEFAULT case is reached either if the parameter topography
!--       contains a wrong character string or if the user has defined a special
!--       case in the user interface. There, the subroutine user_init_grid 
!--       checks which of these two conditions applies.
          CALL user_init_grid( topo )
          CALL filter_topography( topo )

    END SELECT
!
!-- Consistency checks and index array initialization are only required for
!-- non-flat topography.
    IF ( TRIM( topography ) /= 'flat' )  THEN
!
!--    In case of non-flat topography, check whether the convention how to 
!--    define the topography grid has been set correctly, or whether the default
!--    is applicable. If this is not possible, abort.
       IF ( TRIM( topography_grid_convention ) == ' ' )  THEN
          IF ( TRIM( topography ) /= 'single_building' .AND.                   &
               TRIM( topography ) /= 'single_street_canyon' .AND.              &
               TRIM( topography ) /= 'tunnel'  .AND.                           &
               TRIM( topography ) /= 'read_from_file')  THEN
!--          The default value is not applicable here, because it is only valid
!--          for the four standard cases 'single_building', 
!--          'single_street_canyon', 'tunnel' and 'read_from_file'
!--          defined in init_grid.
             WRITE( message_string, * )                                        &
               'The value for "topography_grid_convention" ',                  &
               'is not set. Its default value is & only valid for ',           &
               '"topography" = ''single_building'', ''tunnel'' ',              &
               '''single_street_canyon'' & or ''read_from_file''.',            &
               '& Choose ''cell_edge'' or ''cell_center''.'
             CALL message( 'init_grid', 'PA0239', 1, 2, 0, 6, 0 )
          ELSE
!--          The default value is applicable here.
!--          Set convention according to topography.
             IF ( TRIM( topography ) == 'single_building' .OR.                 &
                  TRIM( topography ) == 'single_street_canyon' )  THEN
                topography_grid_convention = 'cell_edge'
             ELSEIF ( TRIM( topography ) == 'read_from_file'  .OR.             &
                      TRIM( topography ) == 'tunnel')  THEN
                topography_grid_convention = 'cell_center'
             ENDIF
          ENDIF
       ELSEIF ( TRIM( topography_grid_convention ) /= 'cell_edge' .AND.        &
                TRIM( topography_grid_convention ) /= 'cell_center' )  THEN
          WRITE( message_string, * )                                           &
            'The value for "topography_grid_convention" is ',                  &
            'not recognized.& Choose ''cell_edge'' or ''cell_center''.'
          CALL message( 'init_grid', 'PA0240', 1, 2, 0, 6, 0 )
       ENDIF


       IF ( topography_grid_convention == 'cell_edge' )  THEN
! 
!--       The array nzb_local as defined using the 'cell_edge' convention 
!--       describes the actual total size of topography which is defined at the 
!--       cell edges where u=0 on the topography walls in x-direction and v=0 
!--       on the topography walls in y-direction. However, PALM uses individual
!--       arrays nzb_u|v|w|s_inner|outer that are based on nzb_s_inner.
!--       Therefore, the extent of topography in nzb_local is now reduced by 
!--       1dx at the E topography walls and by 1dy at the N topography walls 
!--       to form the basis for nzb_s_inner. 
!--       Note, the reverse memory access (i-j instead of j-i) is absolutely
!--       required at this point.
          DO  j = nys+1, nyn+1
             DO  i = nxl-1, nxr
                DO  k = nzb, nzt+1
                   IF ( BTEST( topo(k,j,i), 0 )  .OR.                          &
                        BTEST( topo(k,j,i+1), 0 ) )                            &
                       topo(k,j,i) = IBSET( topo(k,j,i), 0 )
                ENDDO
             ENDDO
          ENDDO      
          CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )

          DO  i = nxl, nxr+1
             DO  j = nys-1, nyn
                DO  k = nzb, nzt+1
                   IF ( BTEST( topo(k,j,i), 0 )  .OR.                          &
                        BTEST( topo(k,j+1,i), 0 ) )                            &
                      topo(k,j,i) = IBSET( topo(k,j,i), 0 )
                ENDDO
             ENDDO
          ENDDO  
          CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )
   
       ENDIF
    ENDIF

 END SUBROUTINE init_topo

 SUBROUTINE set_topo_flags(topo)

    USE control_parameters,                                                    &
        ONLY:  bc_lr_cyc, bc_ns_cyc, constant_flux_layer, land_surface,        &
               use_surface_fluxes, use_top_fluxes, urban_surface

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg, nz,   &
               nzb, nzt, wall_flags_0

    USE kinds

    IMPLICIT NONE

    INTEGER(iwp) ::  i             !< index variable along x
    INTEGER(iwp) ::  j             !< index variable along y
    INTEGER(iwp) ::  k             !< index variable along z

    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  topo !< input array for 3D topography and dummy array for setting "outer"-flags

    ALLOCATE( wall_flags_0(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    wall_flags_0 = 0
!
!-- Set-up topography flags. First, set flags only for s, u, v and w-grid.
!-- Further special flags will be set in following loops. 
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb, nzt+1
!
!--          scalar grid
             IF ( BTEST( topo(k,j,i), 0 ) )                                    &
                wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 0 )
!
!--          u grid
             IF ( BTEST( topo(k,j,i),   0 )  .AND.                             &
                  BTEST( topo(k,j,i-1), 0 ) )                                  &
                wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 1 )
!
!--          v grid
             IF ( BTEST( topo(k,j,i),   0 )  .AND.                             &
                  BTEST( topo(k,j-1,i), 0 ) )                                  &
                 wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 2 )

          ENDDO

          DO k = nzb, nzt
!
!--          w grid
             IF ( BTEST( topo(k,j,i),   0 )  .AND.                             &
                  BTEST( topo(k+1,j,i), 0 ) )                                  &
                wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 3 )
          ENDDO
          IF ( TRIM(constant_flux_layer) /= 'top' )                            &
             wall_flags_0(nzt+1,j,i) = IBSET( wall_flags_0(nzt+1,j,i), 3 )

       ENDDO
    ENDDO

    CALL exchange_horiz_int( wall_flags_0, nys, nyn, nxl, nxr, nzt, nbgp )
!
!-- Set outer array for scalars to mask near-surface grid points in 
!-- production_e
    DO i = nxl, nxr
       DO j = nys, nyn
          DO k = nzb, nzt+1
             IF ( BTEST( wall_flags_0(k,j-1,i), 0 )  .AND.                     &
                  BTEST( wall_flags_0(k,j+1,i), 0 )  .AND.                     &
                  BTEST( wall_flags_0(k,j,i-1), 0 )  .AND.                     &
                  BTEST( wall_flags_0(k,j-1,i-1), 0 )  .AND.                   &
                  BTEST( wall_flags_0(k,j+1,i-1), 0 )  .AND.                   &
                  BTEST( wall_flags_0(k,j-1,i+1), 0 )  .AND.                   &
                  BTEST( wall_flags_0(k,j+1,i+1), 0 ) )                        &
                wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 24 )
          ENDDO
       ENDDO
    ENDDO
!
!-- Set further special flags
    DO i = nxl, nxr
       DO j = nys, nyn
          DO k = nzb, nzt+1
!
!--          scalar grid, former nzb_diff_s_inner.
!--          Note, use this flag also to mask topography in diffusion_u and 
!--          diffusion_v along the vertical direction. In case of 
!--          use_surface_fluxes, fluxes are calculated via MOST, else, simple
!--          gradient approach is applied. Please note, in case of u- and v-
!--          diffuison, a small error is made at edges (on the east side for u,
!--          at the north side for v), since topography on scalar grid point
!--          is used instead of topography on u/v-grid. As number of topography grid
!--          points on uv-grid is different than s-grid, different number of 
!--          surface elements would be required. In order to avoid this, 
!--          treat edges (u(k,j,i+1)) simply by a gradient approach, i.e. these
!--          points are not masked within diffusion_u. Tests had shown that the
!--          effect on the flow is negligible. 
             IF ( TRIM(constant_flux_layer) /= 'none' .OR. use_surface_fluxes )  THEN
                IF ( BTEST( wall_flags_0(k,j,i), 0 ) ) THEN
                   wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 8 )
                ENDIF
             ELSE
                wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 8 )
             ENDIF

          ENDDO
!
!--       Special flag to control vertical diffusion at model top - former 
!--       nzt_diff
          wall_flags_0(:,j,i) = IBSET( wall_flags_0(:,j,i), 9 )
          IF ( TRIM(constant_flux_layer) == 'top' .OR. use_top_fluxes ) THEN
             wall_flags_0(nzt+1,j,i) = IBCLR( wall_flags_0(nzt+1,j,i), 9 )
          ENDIF

          DO k = nzb+1, nzt
!
!--          Special flag on u grid, former nzb_u_inner + 1, required   
!--          for disturb_field and initialization. Do not disturb directly at
!--          topography, as well as initialize u with zero one grid point outside
!--          of topography.
             IF ( BTEST( wall_flags_0(k-1,j,i), 1 )  .AND.                     &
                  BTEST( wall_flags_0(k,j,i),   1 )  .AND.                     &
                  BTEST( wall_flags_0(k+1,j,i), 1 ) )                          &
                wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 20 )
!
!--          Special flag on v grid, former nzb_v_inner + 1, required   
!--          for disturb_field and initialization. Do not disturb directly at
!--          topography, as well as initialize v with zero one grid point outside
!--          of topography.
             IF ( BTEST( wall_flags_0(k-1,j,i), 2 )  .AND.                     &
                  BTEST( wall_flags_0(k,j,i),   2 )  .AND.                     &
                  BTEST( wall_flags_0(k+1,j,i), 2 ) )                          &
                wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 21 )
!
!--          Special flag on scalar grid, former nzb_s_inner+1. Used for 
!--          lpm_sgs_tke
             IF ( BTEST( wall_flags_0(k,j,i),   0 )  .AND.                     &
                  BTEST( wall_flags_0(k-1,j,i), 0 )  .AND.                     &
                  BTEST( wall_flags_0(k+1,j,i), 0 ) )                          &
                wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 25 )
!
!--          Special flag on scalar grid, nzb_diff_s_outer - 1, required in 
!--          in production_e
!--          CB: but bit 29 isn't even used if constant_flux_layer is on
             IF ( TRIM(constant_flux_layer) /= 'none' .OR.                     &
                  use_surface_fluxes .OR. use_top_fluxes )  THEN
                IF ( BTEST( wall_flags_0(k,j,i),   24 )  .AND.                 &
                     BTEST( wall_flags_0(k-1,j,i), 24 )  .AND.                 &
                     BTEST( wall_flags_0(k+1,j,i), 0 ) )                       &
                   wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 29 )
             ELSE
                IF ( BTEST( wall_flags_0(k,j,i), 0 ) )                         &
                   wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 29 )
             ENDIF
!
!--          Special flag on scalar grid, nzb_diff_s_outer - 1, required in 
!--          in production_e
             IF ( TRIM(constant_flux_layer) /= 'none' .OR.                     &
                  use_surface_fluxes .OR. use_top_fluxes )  THEN
                IF ( BTEST( wall_flags_0(k,j,i),   0 )  .AND.                  &
                     BTEST( wall_flags_0(k-1,j,i), 0 )  .AND.                  &
                     BTEST( wall_flags_0(k+1,j,i), 0 ) )                       &
                   wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 30 )
             ELSE
                IF ( BTEST( wall_flags_0(k,j,i), 0 ) )                         &
                   wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 30 )
             ENDIF
          ENDDO
!
!--       Flags indicating downward facing walls
          DO k = nzb+1, nzt+1
!
!--          Scalar grid
             IF ( BTEST( wall_flags_0(k-1,j,i), 0 )  .AND.                     &
            .NOT. BTEST( wall_flags_0(k,j,i), 0   ) )                          & 
                 wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 13 ) 
!
!--          Downward facing wall on u grid
             IF ( BTEST( wall_flags_0(k-1,j,i), 1 )  .AND.                     &
            .NOT. BTEST( wall_flags_0(k,j,i), 1   ) )                          & 
                wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 15 )
!
!--          Downward facing wall on v grid
             IF ( BTEST( wall_flags_0(k-1,j,i), 2 )  .AND.                     &
            .NOT. BTEST( wall_flags_0(k,j,i), 2   ) )                          & 
                wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 17 )
!
!--          Downward facing wall on w grid
             IF ( BTEST( wall_flags_0(k-1,j,i), 3 )  .AND.                     &
            .NOT. BTEST( wall_flags_0(k,j,i), 3 ) )                            & 
                wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 19 )
          ENDDO
!
!--       Flags indicating upward facing walls
          DO k = nzb, nzt
!
!--          Upward facing wall on scalar grid
             IF ( .NOT. BTEST( wall_flags_0(k,j,i),   0 )  .AND.               &
                        BTEST( wall_flags_0(k+1,j,i), 0 ) )                    & 
                wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 12 )
!
!--          Upward facing wall on u grid
             IF ( .NOT. BTEST( wall_flags_0(k,j,i),   1 )  .AND.               &
                        BTEST( wall_flags_0(k+1,j,i), 1 ) )                    & 
                wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 14 )

!   
!--          Upward facing wall on v grid
             IF ( .NOT. BTEST( wall_flags_0(k,j,i),   2 )  .AND.               &
                        BTEST( wall_flags_0(k+1,j,i), 2 ) )                    & 
                wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 16 )
   
!
!--          Upward facing wall on w grid
             IF ( .NOT. BTEST( wall_flags_0(k,j,i),   3 )  .AND.               &
                        BTEST( wall_flags_0(k+1,j,i), 3 ) )                    & 
                wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 18 )
!
!--          Special flag on scalar grid, former nzb_s_inner
             IF ( BTEST( wall_flags_0(k,j,i), 0 )  .OR.                        &
                  BTEST( wall_flags_0(k,j,i), 12 ) .OR.                        &
                  BTEST( wall_flags_0(k,j,i), 13 ) )                           &
                   wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 22 )
!
!--          Special flag on scalar grid, nzb_diff_s_inner - 1, required for 
!--          flow_statistics
             IF ( TRIM(constant_flux_layer) == 'bottom' .OR.                   &
                  use_surface_fluxes                          )  THEN
                IF ( BTEST( wall_flags_0(k,j,i),   0 )  .AND.                  &
                     BTEST( wall_flags_0(k+1,j,i), 0 ) )                       &
                  wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 23 )
             ELSEIF ( TRIM(constant_flux_layer) == 'top' .OR. use_top_fluxes ) &
                THEN
                IF ( BTEST( wall_flags_0(k,j,i),   0 )  .AND.                  &
                     BTEST( wall_flags_0(k-1,j,i), 0 ) )                       &
                  wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 23 )
             ELSE
                IF ( BTEST( wall_flags_0(k,j,i), 22 ) )                        &
                   wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 23 )
             ENDIF

          ENDDO
          
          IF ( TRIM(constant_flux_layer) /= 'top' ) THEN
             wall_flags_0(nzt+1,j,i) = IBSET( wall_flags_0(nzt+1,j,i), 22 )
             wall_flags_0(nzt+1,j,i) = IBSET( wall_flags_0(nzt+1,j,i), 23 )
          ENDIF

       ENDDO
    ENDDO
!
!-- Finally, set identification flags indicating natural terrain or buildings.
!-- Natural terrain grid points.
    IF ( land_surface )  THEN
       DO i = nxl, nxr
          DO j = nys, nyn
             DO k = nzb, nzt+1
!
!--             Natural terrain grid point
                IF ( BTEST( topo(k,j,i), 1 ) )                                 &
                   wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 5 )
             ENDDO
          ENDDO
       ENDDO
    ENDIF
!
!-- Building grid points.
    IF ( urban_surface )  THEN
       DO i = nxl, nxr
          DO j = nys, nyn
             DO k = nzb, nzt+1
                IF ( BTEST( topo(k,j,i), 2 ) )                                 &
                   wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 6 )
             ENDDO
          ENDDO
       ENDDO
    ENDIF
!
!-- Exchange ghost points for wall flags
    CALL exchange_horiz_int( wall_flags_0, nys, nyn, nxl, nxr, nzt, nbgp )
!
!-- Set boundary conditions also for flags. Can be interpreted as Neumann
!-- boundary conditions for topography. 
    IF ( .NOT. bc_ns_cyc )  THEN
       IF ( nys == 0  )  THEN
          DO  i = 1, nbgp     
             wall_flags_0(:,nys-i,:)   = wall_flags_0(:,nys,:)
          ENDDO
       ENDIF
       IF ( nyn == ny )  THEN
          DO  i = 1, nbgp  
             wall_flags_0(:,nyn+i,:) = wall_flags_0(:,nyn,:)
          ENDDO
       ENDIF
    ENDIF
    IF ( .NOT. bc_lr_cyc )  THEN
       IF ( nxl == 0  )  THEN
          DO  i = 1, nbgp   
             wall_flags_0(:,:,nxl-i)   = wall_flags_0(:,:,nxl)
          ENDDO
       ENDIF
       IF ( nxr == nx )  THEN 
          DO  i = 1, nbgp   
             wall_flags_0(:,:,nxr+i) = wall_flags_0(:,:,nxr)     
          ENDDO
       ENDIF      
       
    ENDIF


 END SUBROUTINE set_topo_flags



