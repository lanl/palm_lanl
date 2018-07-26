!> @file lpm_advec.f90
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
! $Id: lpm_advec.f90 3065 2018-06-12 07:03:02Z Giersch $
! dz values were replaced by dzw or dz(1) to allow for right vertical stretching
! 
! 2969 2018-04-13 11:55:09Z thiele
! Bugfix in Interpolation indices.
! 
! 2886 2018-03-14 11:51:53Z thiele
! Bugfix in passive particle SGS Model:
! Sometimes the added SGS velocities would lead to a violation of the CFL 
! criterion for single particles. For this a check was added after the 
! calculation of SGS velocities.
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2701 2017-12-15 15:40:50Z suehring
! Changes from last commit documented
! 
! 2698 2017-12-14 18:46:24Z suehring
! Particle interpolations at walls in case of SGS velocities revised and not
! required parts are removed. (responsible Philipp Thiele)
! Bugfix in get_topography_top_index
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2630 2017-11-20 12:58:20Z schwenkel
! Removed indices ilog and jlog which are no longer needed since particle box
! locations are identical to scalar boxes and topography.
! 
! 2628 2017-11-20 12:40:38Z raasch
! bugfix in logarithmic interpolation of v-component (usws was used by mistake)
! 
! 2606 2017-11-10 10:36:31Z schwenkel
! Changed particle box locations: center of particle box now coincides 
! with scalar grid point of same index.
! Renamed module and subroutines: lpm_pack_arrays_mod -> lpm_pack_and_sort_mod
! lpm_pack_all_arrays -> lpm_sort_in_subboxes, lpm_pack_arrays -> lpm_pack
! lpm_sort -> lpm_sort_timeloop_done 
! 
! 2417 2017-09-06 15:22:27Z suehring
! Particle loops adapted for sub-box structure, i.e. for each sub-box the 
! particle loop runs from start_index up to end_index instead from 1 to 
! number_of_particles. This way, it is possible to skip unnecessary 
! computations for particles that already completed the LES timestep. 
! 
! 2318 2017-07-20 17:27:44Z suehring
! Get topography top index via Function call 
! 
! 2317 2017-07-20 17:27:19Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography and surface concept
! 
! 2100 2017-01-05 16:40:16Z suehring
! Prevent extremely large SGS-velocities in regions where TKE is zero, e.g. 
! at the begin of simulations and/or in non-turbulent regions. 
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1936 2016-06-13 13:37:44Z suehring
! Formatting adjustments
! 
! 1929 2016-06-09 16:25:25Z suehring
! Put stochastic equation in an extra subroutine.
! Set flag for stochastic equation to communicate whether a particle is near 
! topography. This case, memory and drift term are disabled in the Weil equation. 
!
! Enable vertical logarithmic interpolation also above topography. This case, 
! set a lower limit for the friction velocity, as it can become very small 
! in narrow street canyons, leading to too large particle speeds.
!
! 1888 2016-04-21 12:20:49Z suehring
! Bugfix concerning logarithmic interpolation of particle speed
!
! 1822 2016-04-07 07:49:42Z hoffmann
! Random velocity fluctuations for particles added. Terminal fall velocity
! for droplets is calculated from a parameterization (which is better than
! the previous, physically correct calculation, which demands a very short
! time step that is not used in the model).
!
! Unused variables deleted.
!
! 1691 2015-10-26 16:17:44Z maronga
! Renamed prandtl_layer to constant_flux_layer.
!
! 1685 2015-10-08 07:32:13Z raasch
! TKE check for negative values (so far, only zero value was checked)
! offset_ocean_nzt_m1 removed
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1583 2015-04-15 12:16:27Z suehring
! Bugfix: particle advection within Prandtl-layer in case of Galilei 
! transformation.
!
! 1369 2014-04-24 05:57:38Z raasch
! usage of module interfaces removed
!
! 1359 2014-04-11 17:15:14Z hoffmann
! New particle structure integrated. 
! Kind definition added to all floating point numbers.
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL constants defined as wp_kind
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1314 2014-03-14 18:25:17Z suehring
! Vertical logarithmic interpolation of horizontal particle speed for particles 
! between roughness height and first vertical grid level.
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 849 2012-03-15 10:35:09Z raasch
! initial revision (former part of advec_particles)
!
!
! Description:
! ------------
!> Calculation of new particle positions due to advection using a simple Euler
!> scheme. Particles may feel inertia effects. SGS transport can be included
!> using the stochastic model of Weil et al. (2004, JAS, 61, 2877-2887).
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_advec (ip,jp,kp)
 

    USE arrays_3d,                                                             &
        ONLY:  de_dx, de_dy, de_dz, diss, dzw, e, km, u, v, w, zu, zw

    USE cpulog

    USE pegrid

    USE control_parameters,                                                    &
        ONLY:  atmos_ocean_sign, cloud_droplets, constant_flux_layer, dt_3d,   &
               dt_3d_reached_l, dz, g, kappa, topography, u_gtrans, v_gtrans

    USE grid_variables,                                                        &
        ONLY:  ddx, dx, ddy, dy
        
    USE indices,                                                               &
        ONLY:  nzb, nzt, wall_flags_0
        
    USE kinds
    
    USE particle_attributes,                                                   &
        ONLY:  block_offset, c_0, dt_min_part, grid_particles,                 &
               iran_part, log_z_z0, number_of_particles, number_of_sublayers,  &
               particles, particle_groups, offset_ocean_nzt, sgs_wf_part,      &
               use_sgs_for_particles, vertical_particle_advection, z0_av_global
        
    USE statistics,                                                            &
        ONLY:  hom

    USE surface_mod,                                                           &
        ONLY:  get_topography_top_index_ji, surf_def_h, surf_lsm_h, surf_usm_h

    IMPLICIT NONE

    LOGICAL ::  subbox_at_wall !< flag to see if the current subgridbox is adjacent to a wall

    INTEGER(iwp) ::  agp                         !< loop variable
    INTEGER(iwp) ::  gp_outside_of_building(1:8) !< number of grid points used for particle interpolation in case of topography 
    INTEGER(iwp) ::  i                           !< index variable along x
    INTEGER(iwp) ::  ip                          !< index variable along x 
    INTEGER(iwp) ::  j                           !< index variable along y 
    INTEGER(iwp) ::  jp                          !< index variable along y 
    INTEGER(iwp) ::  k                           !< index variable along z 
    INTEGER(iwp) ::  k_wall                      !< vertical index of topography top 
    INTEGER(iwp) ::  kp                          !< index variable along z 
    INTEGER(iwp) ::  kw                          !< index variable along z
    INTEGER(iwp) ::  n                           !< loop variable over all particles in a grid box
    INTEGER(iwp) ::  nb                          !< block number particles are sorted in
    INTEGER(iwp) ::  num_gp                      !< number of adjacent grid points inside topography 
    INTEGER(iwp) ::  surf_start                  !< Index on surface data-type for current grid box 

    INTEGER(iwp), DIMENSION(0:7) ::  start_index !< start particle index for current block
    INTEGER(iwp), DIMENSION(0:7) ::  end_index   !< start particle index for current block

    REAL(wp) ::  aa                 !< dummy argument for horizontal particle interpolation 
    REAL(wp) ::  bb                 !< dummy argument for horizontal particle interpolation 
    REAL(wp) ::  cc                 !< dummy argument for horizontal particle interpolation 
    REAL(wp) ::  d_sum              !< dummy argument for horizontal particle interpolation in case of topography 
    REAL(wp) ::  d_z_p_z0           !< inverse of interpolation length for logarithmic interpolation 
    REAL(wp) ::  dd                 !< dummy argument for horizontal particle interpolation 
    REAL(wp) ::  de_dx_int_l        !< x/y-interpolated TKE gradient (x) at particle position at lower vertical level
    REAL(wp) ::  de_dx_int_u        !< x/y-interpolated TKE gradient (x) at particle position at upper vertical level
    REAL(wp) ::  de_dy_int_l        !< x/y-interpolated TKE gradient (y) at particle position at lower vertical level
    REAL(wp) ::  de_dy_int_u        !< x/y-interpolated TKE gradient (y) at particle position at upper vertical level
    REAL(wp) ::  de_dt              !< temporal derivative of TKE experienced by the particle
    REAL(wp) ::  de_dt_min          !< lower level for temporal TKE derivative
    REAL(wp) ::  de_dz_int_l        !< x/y-interpolated TKE gradient (z) at particle position at lower vertical level
    REAL(wp) ::  de_dz_int_u        !< x/y-interpolated TKE gradient (z) at particle position at upper vertical level
    REAL(wp) ::  diameter           !< diamter of droplet
    REAL(wp) ::  diss_int_l         !< x/y-interpolated dissipation at particle position at lower vertical level
    REAL(wp) ::  diss_int_u         !< x/y-interpolated dissipation at particle position at upper vertical level
    REAL(wp) ::  dt_particle_m      !< previous particle time step
    REAL(wp) ::  dz_temp            !< 
    REAL(wp) ::  e_int_l            !< x/y-interpolated TKE at particle position at lower vertical level
    REAL(wp) ::  e_int_u            !< x/y-interpolated TKE at particle position at upper vertical level
    REAL(wp) ::  e_mean_int         !< horizontal mean TKE at particle height
    REAL(wp) ::  exp_arg            !<
    REAL(wp) ::  exp_term           !<
    REAL(wp) ::  gg                 !< dummy argument for horizontal particle interpolation 
    REAL(wp) ::  height_p           !< dummy argument for logarithmic interpolation
    REAL(wp) ::  location(1:30,1:3) !< wall locations
    REAL(wp) ::  log_z_z0_int       !< logarithmus used for surface_layer interpolation
    REAL(wp) ::  random_gauss       !<
    REAL(wp) ::  RL                 !< Lagrangian autocorrelation coefficient
    REAL(wp) ::  rg1                !< Gaussian distributed random number
    REAL(wp) ::  rg2                !< Gaussian distributed random number
    REAL(wp) ::  rg3                !< Gaussian distributed random number
    REAL(wp) ::  sigma              !< velocity standard deviation
    REAL(wp) ::  u_int_l            !< x/y-interpolated u-component at particle position at lower vertical level
    REAL(wp) ::  u_int_u            !< x/y-interpolated u-component at particle position at upper vertical level
    REAL(wp) ::  us_int             !< friction velocity at particle grid box
    REAL(wp) ::  usws_int           !< surface momentum flux (u component) at particle grid box
    REAL(wp) ::  v_int_l            !< x/y-interpolated v-component at particle position at lower vertical level
    REAL(wp) ::  v_int_u            !< x/y-interpolated v-component at particle position at upper vertical level
    REAL(wp) ::  vsws_int           !< surface momentum flux (u component) at particle grid box
    REAL(wp) ::  vv_int             !<
    REAL(wp) ::  w_int_l            !< x/y-interpolated w-component at particle position at lower vertical level
    REAL(wp) ::  w_int_u            !< x/y-interpolated w-component at particle position at upper vertical level
    REAL(wp) ::  w_s                !< terminal velocity of droplets
    REAL(wp) ::  x                  !< dummy argument for horizontal particle interpolation 
    REAL(wp) ::  y                  !< dummy argument for horizontal particle interpolation 
    REAL(wp) ::  z_p                !< surface layer height (0.5 dz)

    REAL(wp), PARAMETER ::  a_rog = 9.65_wp      !< parameter for fall velocity
    REAL(wp), PARAMETER ::  b_rog = 10.43_wp     !< parameter for fall velocity
    REAL(wp), PARAMETER ::  c_rog = 0.6_wp       !< parameter for fall velocity
    REAL(wp), PARAMETER ::  k_cap_rog = 4.0_wp   !< parameter for fall velocity
    REAL(wp), PARAMETER ::  k_low_rog = 12.0_wp  !< parameter for fall velocity
    REAL(wp), PARAMETER ::  d0_rog = 0.745_wp    !< separation diameter

    REAL(wp), DIMENSION(1:30) ::  d_gp_pl !< dummy argument for particle interpolation scheme in case of topography
    REAL(wp), DIMENSION(1:30) ::  de_dxi  !< horizontal TKE gradient along x at adjacent wall 
    REAL(wp), DIMENSION(1:30) ::  de_dyi  !< horizontal TKE gradient along y at adjacent wall 
    REAL(wp), DIMENSION(1:30) ::  de_dzi  !< horizontal TKE gradient along z at adjacent wall 
    REAL(wp), DIMENSION(1:30) ::  dissi   !< dissipation at adjacent wall
    REAL(wp), DIMENSION(1:30) ::  ei      !< TKE at adjacent wall

    REAL(wp), DIMENSION(number_of_particles) ::  term_1_2       !< flag to communicate whether a particle is near topography or not
    REAL(wp), DIMENSION(number_of_particles) ::  dens_ratio     !<
    REAL(wp), DIMENSION(number_of_particles) ::  de_dx_int      !< horizontal TKE gradient along x at particle position
    REAL(wp), DIMENSION(number_of_particles) ::  de_dy_int      !< horizontal TKE gradient along y at particle position
    REAL(wp), DIMENSION(number_of_particles) ::  de_dz_int      !< horizontal TKE gradient along z at particle position
    REAL(wp), DIMENSION(number_of_particles) ::  diss_int       !< dissipation at particle position
    REAL(wp), DIMENSION(number_of_particles) ::  dt_gap         !< remaining time until particle time integration reaches LES time
    REAL(wp), DIMENSION(number_of_particles) ::  dt_particle    !< particle time step
    REAL(wp), DIMENSION(number_of_particles) ::  e_int          !< TKE at particle position
    REAL(wp), DIMENSION(number_of_particles) ::  fs_int         !< weighting factor for subgrid-scale particle speed
    REAL(wp), DIMENSION(number_of_particles) ::  lagr_timescale !< Lagrangian timescale
    REAL(wp), DIMENSION(number_of_particles) ::  rvar1_temp     !< 
    REAL(wp), DIMENSION(number_of_particles) ::  rvar2_temp     !< 
    REAL(wp), DIMENSION(number_of_particles) ::  rvar3_temp     !< 
    REAL(wp), DIMENSION(number_of_particles) ::  u_int          !< u-component of particle speed
    REAL(wp), DIMENSION(number_of_particles) ::  v_int          !< v-component of particle speed 
    REAL(wp), DIMENSION(number_of_particles) ::  w_int          !< w-component of particle speed
    REAL(wp), DIMENSION(number_of_particles) ::  xv             !< x-position
    REAL(wp), DIMENSION(number_of_particles) ::  yv             !< y-position
    REAL(wp), DIMENSION(number_of_particles) ::  zv             !< z-position

    REAL(wp), DIMENSION(number_of_particles, 3) ::  rg !< vector of Gaussian distributed random numbers 

    CALL cpu_log( log_point_s(44), 'lpm_advec', 'continue' )

!
!-- Determine height of Prandtl layer and distance between Prandtl-layer
!-- height and horizontal mean roughness height, which are required for 
!-- vertical logarithmic interpolation of horizontal particle speeds 
!-- (for particles below first vertical grid level).
    z_p      = zu(nzb+1) - zw(nzb)
    d_z_p_z0 = 1.0_wp / ( z_p - z0_av_global )

    start_index = grid_particles(kp,jp,ip)%start_index
    end_index   = grid_particles(kp,jp,ip)%end_index

    xv = particles(1:number_of_particles)%x
    yv = particles(1:number_of_particles)%y
    zv = particles(1:number_of_particles)%z

    DO  nb = 0, 7
!
!--    Interpolate u velocity-component       
       i = ip
       j = jp + block_offset(nb)%j_off
       k = kp + block_offset(nb)%k_off

       DO  n = start_index(nb), end_index(nb)
!
!--       Interpolation of the u velocity component onto particle position.  
!--       Particles are interpolation bi-linearly in the horizontal and a 
!--       linearly in the vertical. An exception is made for particles below
!--       the first vertical grid level in case of a prandtl layer. In this 
!--       case the horizontal particle velocity components are determined using
!--       Monin-Obukhov relations (if branch). 
!--       First, check if particle is located below first vertical grid level 
!--       above topography (Prandtl-layer height)
!--       Determine vertical index of topography top
          k_wall = get_topography_top_index_ji( jp, ip, 's' )

          IF ( constant_flux_layer  .AND.  zv(n) - zw(k_wall) < z_p )  THEN
!
!--          Resolved-scale horizontal particle velocity is zero below z0. 
             IF ( zv(n) - zw(k_wall) < z0_av_global )  THEN
                u_int(n) = 0.0_wp
             ELSE
!
!--             Determine the sublayer. Further used as index.
                height_p = ( zv(n) - zw(k_wall) - z0_av_global ) &
                                     * REAL( number_of_sublayers, KIND=wp )    &
                                     * d_z_p_z0 
!
!--             Calculate LOG(z/z0) for exact particle height. Therefore,   
!--             interpolate linearly between precalculated logarithm. 
                log_z_z0_int = log_z_z0(INT(height_p))                         &
                                 + ( height_p - INT(height_p) )                &
                                 * ( log_z_z0(INT(height_p)+1)                 &
                                      - log_z_z0(INT(height_p))                &
                                   ) 
!
!--             Get friction velocity and momentum flux from new surface data 
!--             types.
                IF ( surf_def_h(0)%start_index(jp,ip) <=                   &
                     surf_def_h(0)%end_index(jp,ip) )  THEN
                   surf_start = surf_def_h(0)%start_index(jp,ip)
!--                Limit friction velocity. In narrow canyons or holes the 
!--                friction velocity can become very small, resulting in a too
!--                large particle speed.
                   us_int    = MAX( surf_def_h(0)%us(surf_start), 0.01_wp )  
                   usws_int  = surf_def_h(0)%usws(surf_start)
                ELSEIF ( surf_lsm_h%start_index(jp,ip) <=                  &
                         surf_lsm_h%end_index(jp,ip) )  THEN
                   surf_start = surf_lsm_h%start_index(jp,ip)
                   us_int    = MAX( surf_lsm_h%us(surf_start), 0.01_wp )  
                   usws_int  = surf_lsm_h%usws(surf_start)
                ELSEIF ( surf_usm_h%start_index(jp,ip) <=                  &
                         surf_usm_h%end_index(jp,ip) )  THEN
                   surf_start = surf_usm_h%start_index(jp,ip)
                   us_int    = MAX( surf_usm_h%us(surf_start), 0.01_wp )  
                   usws_int  = surf_usm_h%usws(surf_start)
                ENDIF

!
!--             Neutral solution is applied for all situations, e.g. also for 
!--             unstable and stable situations. Even though this is not exact 
!--             this saves a lot of CPU time since several calls of intrinsic 
!--             FORTRAN procedures (LOG, ATAN) are avoided, This is justified 
!--             as sensitivity studies revealed no significant effect of 
!--             using the neutral solution also for un/stable situations.
                u_int(n) = -usws_int / ( us_int * kappa + 1E-10_wp )           & 
                            * log_z_z0_int - u_gtrans
                
             ENDIF
!
!--       Particle above the first grid level. Bi-linear interpolation in the 
!--       horizontal and linear interpolation in the vertical direction.
          ELSE

             x  = xv(n) - i * dx
             y  = yv(n) + ( 0.5_wp - j ) * dy
             aa = x**2          + y**2
             bb = ( dx - x )**2 + y**2
             cc = x**2          + ( dy - y )**2
             dd = ( dx - x )**2 + ( dy - y )**2
             gg = aa + bb + cc + dd

             u_int_l = ( ( gg - aa ) * u(k,j,i)   + ( gg - bb ) * u(k,j,i+1)   &
                         + ( gg - cc ) * u(k,j+1,i) + ( gg - dd ) *            &
                         u(k,j+1,i+1) ) / ( 3.0_wp * gg ) - u_gtrans

             IF ( k == nzt )  THEN
                u_int(n) = u_int_l
             ELSE
                u_int_u = ( ( gg-aa ) * u(k+1,j,i) + ( gg-bb ) * u(k+1,j,i+1)  &
                            + ( gg-cc ) * u(k+1,j+1,i) + ( gg-dd ) *           &
                            u(k+1,j+1,i+1) ) / ( 3.0_wp * gg ) - u_gtrans
                u_int(n) = u_int_l + ( zv(n) - zu(k) ) / dzw(k) *               &
                           ( u_int_u - u_int_l )
             ENDIF

          ENDIF

       ENDDO
!
!--    Same procedure for interpolation of the v velocity-component 
       i = ip + block_offset(nb)%i_off
       j = jp
       k = kp + block_offset(nb)%k_off

       DO  n = start_index(nb), end_index(nb)

!
!--       Determine vertical index of topography top
          k_wall = get_topography_top_index_ji( jp,ip, 's' )

          IF ( constant_flux_layer  .AND.  zv(n) - zw(k_wall) < z_p )  THEN
             IF ( zv(n) - zw(k_wall) < z0_av_global )  THEN
!
!--             Resolved-scale horizontal particle velocity is zero below z0. 
                v_int(n) = 0.0_wp
             ELSE       
!
!--             Determine the sublayer. Further used as index. Please note, 
!--             logarithmus can not be reused from above, as in in case of
!--             topography particle on u-grid can be above surface-layer height,
!--             whereas it can be below on v-grid. 
                height_p = ( zv(n) - zw(k_wall) - z0_av_global ) &
                                  * REAL( number_of_sublayers, KIND=wp )       &
                                  * d_z_p_z0 
!
!--             Calculate LOG(z/z0) for exact particle height. Therefore,   
!--             interpolate linearly between precalculated logarithm. 
                log_z_z0_int = log_z_z0(INT(height_p))                         &
                                 + ( height_p - INT(height_p) )                &
                                 * ( log_z_z0(INT(height_p)+1)                 &
                                      - log_z_z0(INT(height_p))                &
                                   ) 
!
!--             Get friction velocity and momentum flux from new surface data 
!--             types.
                IF ( surf_def_h(0)%start_index(jp,ip) <=                   &
                     surf_def_h(0)%end_index(jp,ip) )  THEN
                   surf_start = surf_def_h(0)%start_index(jp,ip)
!--                Limit friction velocity. In narrow canyons or holes the 
!--                friction velocity can become very small, resulting in a too
!--                large particle speed.
                   us_int    = MAX( surf_def_h(0)%us(surf_start), 0.01_wp )  
                   vsws_int  = surf_def_h(0)%vsws(surf_start)
                ELSEIF ( surf_lsm_h%start_index(jp,ip) <=                  &
                         surf_lsm_h%end_index(jp,ip) )  THEN
                   surf_start = surf_lsm_h%start_index(jp,ip)
                   us_int    = MAX( surf_lsm_h%us(surf_start), 0.01_wp )  
                   vsws_int  = surf_lsm_h%vsws(surf_start)
                ELSEIF ( surf_usm_h%start_index(jp,ip) <=                  &
                         surf_usm_h%end_index(jp,ip) )  THEN
                   surf_start = surf_usm_h%start_index(jp,ip)
                   us_int    = MAX( surf_usm_h%us(surf_start), 0.01_wp )  
                   vsws_int  = surf_usm_h%vsws(surf_start)
                ENDIF 
!
!--             Neutral solution is applied for all situations, e.g. also for 
!--             unstable and stable situations. Even though this is not exact 
!--             this saves a lot of CPU time since several calls of intrinsic 
!--             FORTRAN procedures (LOG, ATAN) are avoided, This is justified 
!--             as sensitivity studies revealed no significant effect of 
!--             using the neutral solution also for un/stable situations.
                v_int(n) = -vsws_int / ( us_int * kappa + 1E-10_wp )           &
                         * log_z_z0_int - v_gtrans

             ENDIF

          ELSE
             x  = xv(n) + ( 0.5_wp - i ) * dx
             y  = yv(n) - j * dy
             aa = x**2          + y**2
             bb = ( dx - x )**2 + y**2
             cc = x**2          + ( dy - y )**2
             dd = ( dx - x )**2 + ( dy - y )**2
             gg = aa + bb + cc + dd

             v_int_l = ( ( gg - aa ) * v(k,j,i)   + ( gg - bb ) * v(k,j,i+1)   &
                       + ( gg - cc ) * v(k,j+1,i) + ( gg - dd ) * v(k,j+1,i+1) &
                       ) / ( 3.0_wp * gg ) - v_gtrans

             IF ( k == nzt )  THEN
                v_int(n) = v_int_l
             ELSE
                v_int_u = ( ( gg-aa ) * v(k+1,j,i)   + ( gg-bb ) * v(k+1,j,i+1)   &
                          + ( gg-cc ) * v(k+1,j+1,i) + ( gg-dd ) * v(k+1,j+1,i+1) &
                          ) / ( 3.0_wp * gg ) - v_gtrans
                v_int(n) = v_int_l + ( zv(n) - zu(k) ) / dzw(k) *               &
                                  ( v_int_u - v_int_l )
             ENDIF

          ENDIF

       ENDDO
!
!--    Same procedure for interpolation of the w velocity-component
       i = ip + block_offset(nb)%i_off
       j = jp + block_offset(nb)%j_off
       k = kp - 1

       DO  n = start_index(nb), end_index(nb)

          IF ( vertical_particle_advection(particles(n)%group) )  THEN

             x  = xv(n) + ( 0.5_wp - i ) * dx
             y  = yv(n) + ( 0.5_wp - j ) * dy
             aa = x**2          + y**2
             bb = ( dx - x )**2 + y**2
             cc = x**2          + ( dy - y )**2
             dd = ( dx - x )**2 + ( dy - y )**2
             gg = aa + bb + cc + dd

             w_int_l = ( ( gg - aa ) * w(k,j,i)   + ( gg - bb ) * w(k,j,i+1)   &
                       + ( gg - cc ) * w(k,j+1,i) + ( gg - dd ) * w(k,j+1,i+1) &
                       ) / ( 3.0_wp * gg )

             IF ( k == nzt )  THEN
                w_int(n) = w_int_l
             ELSE
                w_int_u = ( ( gg-aa ) * w(k+1,j,i)   + &
                            ( gg-bb ) * w(k+1,j,i+1) + &
                            ( gg-cc ) * w(k+1,j+1,i) + &
                            ( gg-dd ) * w(k+1,j+1,i+1) &
                          ) / ( 3.0_wp * gg )
                w_int(n) = w_int_l + ( zv(n) - zw(k) ) / dzw(k) *               &
                           ( w_int_u - w_int_l )
             ENDIF

          ELSE

             w_int(n) = 0.0_wp

          ENDIF

       ENDDO

    ENDDO

!-- Interpolate and calculate quantities needed for calculating the SGS
!-- velocities
    IF ( use_sgs_for_particles  .AND.  .NOT. cloud_droplets )  THEN
    
       DO  nb = 0,7
         
          subbox_at_wall = .FALSE.         
!
!--       In case of topography check if subbox is adjacent to a wall
          IF ( .NOT. topography == 'flat' ) THEN
             i = ip + MERGE( -1_iwp , 1_iwp, BTEST( nb, 2 ) )
             j = jp + MERGE( -1_iwp , 1_iwp, BTEST( nb, 1 ) )
             k = kp + MERGE( -1_iwp , 1_iwp, BTEST( nb, 0 ) )
             IF ( .NOT. BTEST(wall_flags_0(k,  jp, ip), 0) .OR.                &
                  .NOT. BTEST(wall_flags_0(kp, j,  ip), 0) .OR.                &
                  .NOT. BTEST(wall_flags_0(kp, jp, i ), 0) )                   &
             THEN
                subbox_at_wall = .TRUE.
             ENDIF
          ENDIF
          IF ( subbox_at_wall ) THEN
             e_int(start_index(nb):end_index(nb))     = e(kp,jp,ip) 
             diss_int(start_index(nb):end_index(nb))  = diss(kp,jp,ip)
             de_dx_int(start_index(nb):end_index(nb)) = de_dx(kp,jp,ip)
             de_dy_int(start_index(nb):end_index(nb)) = de_dy(kp,jp,ip)
             de_dz_int(start_index(nb):end_index(nb)) = de_dz(kp,jp,ip)
!
!--          Set flag for stochastic equation.
             term_1_2(start_index(nb):end_index(nb)) = 0.0_wp             
          ELSE
             i = ip + block_offset(nb)%i_off
             j = jp + block_offset(nb)%j_off
             k = kp + block_offset(nb)%k_off

             DO  n = start_index(nb), end_index(nb)
!
!--             Interpolate TKE
                x  = xv(n) + ( 0.5_wp - i ) * dx
                y  = yv(n) + ( 0.5_wp - j ) * dy
                aa = x**2          + y**2
                bb = ( dx - x )**2 + y**2
                cc = x**2          + ( dy - y )**2
                dd = ( dx - x )**2 + ( dy - y )**2
                gg = aa + bb + cc + dd

                e_int_l = ( ( gg-aa ) * e(k,j,i)   + ( gg-bb ) * e(k,j,i+1)   &
                          + ( gg-cc ) * e(k,j+1,i) + ( gg-dd ) * e(k,j+1,i+1) &
                          ) / ( 3.0_wp * gg )

                IF ( k+1 == nzt+1 )  THEN
                   e_int(n) = e_int_l
                ELSE
                   e_int_u = ( ( gg - aa ) * e(k+1,j,i)   + &
                               ( gg - bb ) * e(k+1,j,i+1) + &
                               ( gg - cc ) * e(k+1,j+1,i) + &
                               ( gg - dd ) * e(k+1,j+1,i+1) &
                            ) / ( 3.0_wp * gg )
                   e_int(n) = e_int_l + ( zv(n) - zu(k) ) / dzw(k) *            &
                                     ( e_int_u - e_int_l )
                ENDIF
!
!--             Needed to avoid NaN particle velocities (this might not be
!--             required any more)
                IF ( e_int(n) <= 0.0_wp )  THEN
                   e_int(n) = 1.0E-20_wp
                ENDIF
!
!--             Interpolate the TKE gradient along x (adopt incides i,j,k and
!--             all position variables from above (TKE))
                de_dx_int_l = ( ( gg - aa ) * de_dx(k,j,i)   + &
                                ( gg - bb ) * de_dx(k,j,i+1) + &
                                ( gg - cc ) * de_dx(k,j+1,i) + &
                                ( gg - dd ) * de_dx(k,j+1,i+1) &
                               ) / ( 3.0_wp * gg )

                IF ( ( k+1 == nzt+1 )  .OR.  ( k == nzb ) )  THEN
                   de_dx_int(n) = de_dx_int_l
                ELSE
                   de_dx_int_u = ( ( gg - aa ) * de_dx(k+1,j,i)   + &
                                   ( gg - bb ) * de_dx(k+1,j,i+1) + &
                                   ( gg - cc ) * de_dx(k+1,j+1,i) + &
                                   ( gg - dd ) * de_dx(k+1,j+1,i+1) &
                                  ) / ( 3.0_wp * gg )
                   de_dx_int(n) = de_dx_int_l + ( zv(n) - zu(k) ) / dzw(k) *    &
                                              ( de_dx_int_u - de_dx_int_l )
                ENDIF
!
!--             Interpolate the TKE gradient along y
                de_dy_int_l = ( ( gg - aa ) * de_dy(k,j,i)   + &
                                ( gg - bb ) * de_dy(k,j,i+1) + &
                                ( gg - cc ) * de_dy(k,j+1,i) + &
                                ( gg - dd ) * de_dy(k,j+1,i+1) &
                               ) / ( 3.0_wp * gg )
                IF ( ( k+1 == nzt+1 )  .OR.  ( k == nzb ) )  THEN
                   de_dy_int(n) = de_dy_int_l
                ELSE
                   de_dy_int_u = ( ( gg - aa ) * de_dy(k+1,j,i)   + &
                                   ( gg - bb ) * de_dy(k+1,j,i+1) + &
                                   ( gg - cc ) * de_dy(k+1,j+1,i) + &
                                   ( gg - dd ) * de_dy(k+1,j+1,i+1) &
                                  ) / ( 3.0_wp * gg )
                      de_dy_int(n) = de_dy_int_l + ( zv(n) - zu(k) ) / dzw(k) * &
                                                 ( de_dy_int_u - de_dy_int_l )
                ENDIF

!
!--             Interpolate the TKE gradient along z
                IF ( zv(n) < 0.5_wp * dz(1) )  THEN
                   de_dz_int(n) = 0.0_wp
                ELSE
                   de_dz_int_l = ( ( gg - aa ) * de_dz(k,j,i)   + &
                                   ( gg - bb ) * de_dz(k,j,i+1) + &
                                   ( gg - cc ) * de_dz(k,j+1,i) + &
                                   ( gg - dd ) * de_dz(k,j+1,i+1) &
                                  ) / ( 3.0_wp * gg )

                   IF ( ( k+1 == nzt+1 )  .OR.  ( k == nzb ) )  THEN
                      de_dz_int(n) = de_dz_int_l
                   ELSE
                      de_dz_int_u = ( ( gg - aa ) * de_dz(k+1,j,i)   + &
                                      ( gg - bb ) * de_dz(k+1,j,i+1) + &
                                      ( gg - cc ) * de_dz(k+1,j+1,i) + &
                                      ( gg - dd ) * de_dz(k+1,j+1,i+1) &
                                     ) / ( 3.0_wp * gg )
                      de_dz_int(n) = de_dz_int_l + ( zv(n) - zu(k) ) / dzw(k) * &
                                                 ( de_dz_int_u - de_dz_int_l )
                   ENDIF
                ENDIF

!
!--             Interpolate the dissipation of TKE
                diss_int_l = ( ( gg - aa ) * diss(k,j,i)   + &
                               ( gg - bb ) * diss(k,j,i+1) + &
                               ( gg - cc ) * diss(k,j+1,i) + &
                               ( gg - dd ) * diss(k,j+1,i+1) &
                               ) / ( 3.0_wp * gg )

                IF ( k == nzt )  THEN
                   diss_int(n) = diss_int_l
                ELSE
                   diss_int_u = ( ( gg - aa ) * diss(k+1,j,i)   + &
                                  ( gg - bb ) * diss(k+1,j,i+1) + &
                                  ( gg - cc ) * diss(k+1,j+1,i) + &
                                  ( gg - dd ) * diss(k+1,j+1,i+1) &
                                 ) / ( 3.0_wp * gg )
                   diss_int(n) = diss_int_l + ( zv(n) - zu(k) ) / dzw(k) *      &
                                            ( diss_int_u - diss_int_l )
                ENDIF

!
!--             Set flag for stochastic equation.
                term_1_2(n) = 1.0_wp
             ENDDO
          ENDIF
       ENDDO

       DO nb = 0,7
          i = ip + block_offset(nb)%i_off
          j = jp + block_offset(nb)%j_off
          k = kp + block_offset(nb)%k_off

          DO  n = start_index(nb), end_index(nb)
!
!--          Vertical interpolation of the horizontally averaged SGS TKE and
!--          resolved-scale velocity variances and use the interpolated values
!--          to calculate the coefficient fs, which is a measure of the ratio
!--          of the subgrid-scale turbulent kinetic energy to the total amount
!--          of turbulent kinetic energy.
             IF ( k == 0 )  THEN
                e_mean_int = hom(0,1,8,0)
             ELSE
                e_mean_int = hom(k,1,8,0) +                                    &
                                           ( hom(k+1,1,8,0) - hom(k,1,8,0) ) / &
                                           ( zu(k+1) - zu(k) ) *               &
                                           ( zv(n) - zu(k) )
             ENDIF

             kw = kp - 1

             IF ( k == 0 )  THEN
                aa  = hom(k+1,1,30,0)  * ( zv(n) / &
                                         ( 0.5_wp * ( zu(k+1) - zu(k) ) ) )
                bb  = hom(k+1,1,31,0)  * ( zv(n) / &
                                         ( 0.5_wp * ( zu(k+1) - zu(k) ) ) )
                cc  = hom(kw+1,1,32,0) * ( zv(n) / &
                                         ( 1.0_wp * ( zw(kw+1) - zw(kw) ) ) )
             ELSE
                aa  = hom(k,1,30,0) + ( hom(k+1,1,30,0) - hom(k,1,30,0) ) *    &
                           ( ( zv(n) - zu(k) ) / ( zu(k+1) - zu(k) ) )
                bb  = hom(k,1,31,0) + ( hom(k+1,1,31,0) - hom(k,1,31,0) ) *    &
                           ( ( zv(n) - zu(k) ) / ( zu(k+1) - zu(k) ) )
                cc  = hom(kw,1,32,0) + ( hom(kw+1,1,32,0)-hom(kw,1,32,0) ) *   &
                           ( ( zv(n) - zw(kw) ) / ( zw(kw+1)-zw(kw) ) )
             ENDIF

             vv_int = ( 1.0_wp / 3.0_wp ) * ( aa + bb + cc )
!
!--          Needed to avoid NaN particle velocities. The value of 1.0 is just
!--          an educated guess for the given case.
             IF ( vv_int + ( 2.0_wp / 3.0_wp ) * e_mean_int == 0.0_wp )  THEN
                fs_int(n) = 1.0_wp
             ELSE
                fs_int(n) = ( 2.0_wp / 3.0_wp ) * e_mean_int /                 &
                            ( vv_int + ( 2.0_wp / 3.0_wp ) * e_mean_int )
             ENDIF

          ENDDO
       ENDDO

       DO  nb = 0, 7
          DO  n = start_index(nb), end_index(nb)
             rg(n,1) = random_gauss( iran_part, 5.0_wp )
             rg(n,2) = random_gauss( iran_part, 5.0_wp )
             rg(n,3) = random_gauss( iran_part, 5.0_wp )
          ENDDO
       ENDDO

       DO  nb = 0, 7
          DO  n = start_index(nb), end_index(nb)

!
!--          Calculate the Lagrangian timescale according to Weil et al. (2004).
             lagr_timescale(n) = ( 4.0_wp * e_int(n) + 1E-20_wp ) / &
                              ( 3.0_wp * fs_int(n) * c_0 * diss_int(n) + 1E-20_wp )

!
!--          Calculate the next particle timestep. dt_gap is the time needed to
!--          complete the current LES timestep.
             dt_gap(n) = dt_3d - particles(n)%dt_sum
             dt_particle(n) = MIN( dt_3d, 0.025_wp * lagr_timescale(n), dt_gap(n) )
             particles(n)%aux1 = lagr_timescale(n)
             particles(n)%aux2 = dt_gap(n)
!
!--          The particle timestep should not be too small in order to prevent
!--          the number of particle timesteps of getting too large
             IF ( dt_particle(n) < dt_min_part  .AND.  dt_min_part < dt_gap(n) )  THEN
                dt_particle(n) = dt_min_part
             ENDIF
             rvar1_temp(n) = particles(n)%rvar1
             rvar2_temp(n) = particles(n)%rvar2
             rvar3_temp(n) = particles(n)%rvar3
!
!--          Calculate the SGS velocity components
             IF ( particles(n)%age == 0.0_wp )  THEN
!
!--             For new particles the SGS components are derived from the SGS
!--             TKE. Limit the Gaussian random number to the interval
!--             [-5.0*sigma, 5.0*sigma] in order to prevent the SGS velocities
!--             from becoming unrealistically large.
                rvar1_temp(n) = SQRT( 2.0_wp * sgs_wf_part * e_int(n)          &
                                          + 1E-20_wp ) * ( rg(n,1) - 1.0_wp )
                rvar2_temp(n) = SQRT( 2.0_wp * sgs_wf_part * e_int(n)          &
                                          + 1E-20_wp ) * ( rg(n,2) - 1.0_wp )
                rvar3_temp(n) = SQRT( 2.0_wp * sgs_wf_part * e_int(n)          &
                                          + 1E-20_wp ) * ( rg(n,3) - 1.0_wp )

             ELSE
!
!--             Restriction of the size of the new timestep: compared to the 
!--             previous timestep the increase must not exceed 200%. First,
!--             check if age > age_m, in order to prevent that particles get zero
!--             timestep.
                dt_particle_m = MERGE( dt_particle(n),                         &
                                       particles(n)%age - particles(n)%age_m,  &
                                       particles(n)%age - particles(n)%age_m < &
                                       1E-8_wp )
                IF ( dt_particle(n) > 2.0_wp * dt_particle_m )  THEN
                   dt_particle(n) = 2.0_wp * dt_particle_m
                ENDIF

!--             For old particles the SGS components are correlated with the
!--             values from the previous timestep. Random numbers have also to
!--             be limited (see above).
!--             As negative values for the subgrid TKE are not allowed, the
!--             change of the subgrid TKE with time cannot be smaller than
!--             -e_int(n)/dt_particle. This value is used as a lower boundary
!--             value for the change of TKE 
                de_dt_min = - e_int(n) / dt_particle(n)

                de_dt = ( e_int(n) - particles(n)%e_m ) / dt_particle_m

                IF ( de_dt < de_dt_min )  THEN
                   de_dt = de_dt_min
                ENDIF

                CALL weil_stochastic_eq(rvar1_temp(n), fs_int(n), e_int(n),& 
                                        de_dx_int(n), de_dt, diss_int(n),       &
                                        dt_particle(n), rg(n,1), term_1_2(n) )

                CALL weil_stochastic_eq(rvar2_temp(n), fs_int(n), e_int(n),& 
                                        de_dy_int(n), de_dt, diss_int(n),       &
                                        dt_particle(n), rg(n,2), term_1_2(n) )

                CALL weil_stochastic_eq(rvar3_temp(n), fs_int(n), e_int(n),& 
                                        de_dz_int(n), de_dt, diss_int(n),       &
                                        dt_particle(n), rg(n,3), term_1_2(n) )

             ENDIF

          ENDDO
       ENDDO
!
!--    Check if the added SGS velocities result in a violation of the CFL-
!--    criterion. If yes choose a smaller timestep based on the new velocities
!--    and calculate SGS velocities again
       dz_temp = zw(kp)-zw(kp-1)
       
       DO  nb = 0, 7
          DO  n = start_index(nb), end_index(nb)
             IF ( .NOT. particles(n)%age == 0.0_wp .AND.                       &
                (ABS( u_int(n) + rvar1_temp(n) ) > (dx/dt_particle(n))  .OR.   &
                 ABS( v_int(n) + rvar2_temp(n) ) > (dy/dt_particle(n))  .OR.   &
                 ABS( w_int(n) + rvar3_temp(n) ) > (dz_temp/dt_particle(n)))) THEN
                
                dt_particle(n) = 0.9_wp * MIN(                                 &
                                 ( dx / ABS( u_int(n) + rvar1_temp(n) ) ),     &
                                 ( dy / ABS( v_int(n) + rvar2_temp(n) ) ),     &
                                 ( dz_temp / ABS( w_int(n) + rvar3_temp(n) ) ) )

!
!--             Reset temporary SGS velocites to "current" ones
                rvar1_temp(n) = particles(n)%rvar1
                rvar2_temp(n) = particles(n)%rvar2
                rvar3_temp(n) = particles(n)%rvar3
                
                de_dt_min = - e_int(n) / dt_particle(n)

                de_dt = ( e_int(n) - particles(n)%e_m ) / dt_particle_m

                IF ( de_dt < de_dt_min )  THEN
                   de_dt = de_dt_min
                ENDIF

                CALL weil_stochastic_eq(rvar1_temp(n), fs_int(n), e_int(n),& 
                                        de_dx_int(n), de_dt, diss_int(n),       &
                                        dt_particle(n), rg(n,1), term_1_2(n) )

                CALL weil_stochastic_eq(rvar2_temp(n), fs_int(n), e_int(n),& 
                                        de_dy_int(n), de_dt, diss_int(n),       &
                                        dt_particle(n), rg(n,2), term_1_2(n) )

                CALL weil_stochastic_eq(rvar3_temp(n), fs_int(n), e_int(n),& 
                                        de_dz_int(n), de_dt, diss_int(n),       &
                                        dt_particle(n), rg(n,3), term_1_2(n) )
             ENDIF                           
             
!
!--          Update particle velocites 
             particles(n)%rvar1 = rvar1_temp(n)
             particles(n)%rvar2 = rvar2_temp(n)
             particles(n)%rvar3 = rvar3_temp(n)
             u_int(n) = u_int(n) + particles(n)%rvar1
             v_int(n) = v_int(n) + particles(n)%rvar2
             w_int(n) = w_int(n) + particles(n)%rvar3
!
!--          Store the SGS TKE of the current timelevel which is needed for
!--          for calculating the SGS particle velocities at the next timestep
             particles(n)%e_m = e_int(n)
          ENDDO
       ENDDO
       
    ELSE
!
!--    If no SGS velocities are used, only the particle timestep has to
!--    be set
       dt_particle = dt_3d

    ENDIF

    dens_ratio = particle_groups(particles(1:number_of_particles)%group)%density_ratio

    IF ( ANY( dens_ratio == 0.0_wp ) )  THEN
       DO  nb = 0, 7
          DO  n = start_index(nb), end_index(nb)

!
!--          Particle advection
             IF ( dens_ratio(n) == 0.0_wp )  THEN
!
!--             Pure passive transport (without particle inertia)
                particles(n)%x = xv(n) + u_int(n) * dt_particle(n)
                particles(n)%y = yv(n) + v_int(n) * dt_particle(n)
                particles(n)%z = zv(n) + w_int(n) * dt_particle(n)

                particles(n)%speed_x = u_int(n)
                particles(n)%speed_y = v_int(n)
                particles(n)%speed_z = w_int(n)

             ELSE
!
!--             Transport of particles with inertia
                particles(n)%x = particles(n)%x + particles(n)%speed_x * &
                                                  dt_particle(n)
                particles(n)%y = particles(n)%y + particles(n)%speed_y * &
                                                  dt_particle(n)
                particles(n)%z = particles(n)%z + particles(n)%speed_z * &
                                                  dt_particle(n)

!
!--             Update of the particle velocity
                IF ( cloud_droplets )  THEN
!
!--                Terminal velocity is computed for vertical direction (Rogers et 
!--                al., 1993, J. Appl. Meteorol.)
                   diameter = particles(n)%radius * 2000.0_wp !diameter in mm
                   IF ( diameter <= d0_rog )  THEN
                      w_s = k_cap_rog * diameter * ( 1.0_wp - EXP( -k_low_rog * diameter ) )
                   ELSE
                      w_s = a_rog - b_rog * EXP( -c_rog * diameter )
                   ENDIF

!
!--                If selected, add random velocities following Soelch and Kaercher
!--                (2010, Q. J. R. Meteorol. Soc.)
                   IF ( use_sgs_for_particles )  THEN
                      lagr_timescale(n) = km(kp,jp,ip) / MAX( e(kp,jp,ip), 1.0E-20_wp )
                      RL             = EXP( -1.0_wp * dt_3d / lagr_timescale(n) )
                      sigma          = SQRT( e(kp,jp,ip) )

                      rg1 = random_gauss( iran_part, 5.0_wp ) - 1.0_wp
                      rg2 = random_gauss( iran_part, 5.0_wp ) - 1.0_wp
                      rg3 = random_gauss( iran_part, 5.0_wp ) - 1.0_wp

                      particles(n)%rvar1 = RL * particles(n)%rvar1 +              &
                                           SQRT( 1.0_wp - RL**2 ) * sigma * rg1
                      particles(n)%rvar2 = RL * particles(n)%rvar2 +              &
                                           SQRT( 1.0_wp - RL**2 ) * sigma * rg2
                      particles(n)%rvar3 = RL * particles(n)%rvar3 +              &
                                           SQRT( 1.0_wp - RL**2 ) * sigma * rg3

                      particles(n)%speed_x = u_int(n) + particles(n)%rvar1
                      particles(n)%speed_y = v_int(n) + particles(n)%rvar2
                      particles(n)%speed_z = w_int(n) + particles(n)%rvar3 - w_s
                   ELSE
                      particles(n)%speed_x = u_int(n)
                      particles(n)%speed_y = v_int(n)
                      particles(n)%speed_z = w_int(n) - w_s
                   ENDIF

                ELSE

                   IF ( use_sgs_for_particles )  THEN
                      exp_arg  = particle_groups(particles(n)%group)%exp_arg
                      exp_term = EXP( -exp_arg * dt_particle(n) )
                   ELSE
                      exp_arg  = particle_groups(particles(n)%group)%exp_arg
                      exp_term = particle_groups(particles(n)%group)%exp_term
                   ENDIF
                   particles(n)%speed_x = particles(n)%speed_x * exp_term +         &
                                          u_int(n) * ( 1.0_wp - exp_term )
                   particles(n)%speed_y = particles(n)%speed_y * exp_term +         &
                                          v_int(n) * ( 1.0_wp - exp_term )
                   particles(n)%speed_z = particles(n)%speed_z * exp_term +         &
                                          ( w_int(n) - ( 1.0_wp - dens_ratio(n) ) * &
                                          g / exp_arg ) * ( 1.0_wp - exp_term )
                ENDIF

             ENDIF
          ENDDO
       ENDDO
    
    ELSE

       DO  nb = 0, 7
          DO  n = start_index(nb), end_index(nb)
!
!--          Transport of particles with inertia
             particles(n)%x = xv(n) + particles(n)%speed_x * dt_particle(n)
             particles(n)%y = yv(n) + particles(n)%speed_y * dt_particle(n)
             particles(n)%z = zv(n) + particles(n)%speed_z * dt_particle(n)
!
!--          Update of the particle velocity
             IF ( cloud_droplets )  THEN
!
!--             Terminal velocity is computed for vertical direction (Rogers et al., 
!--             1993, J. Appl. Meteorol.)
                diameter = particles(n)%radius * 2000.0_wp !diameter in mm
                IF ( diameter <= d0_rog )  THEN
                   w_s = k_cap_rog * diameter * ( 1.0_wp - EXP( -k_low_rog * diameter ) )
                ELSE
                   w_s = a_rog - b_rog * EXP( -c_rog * diameter )
                ENDIF

!
!--             If selected, add random velocities following Soelch and Kaercher
!--             (2010, Q. J. R. Meteorol. Soc.)
                IF ( use_sgs_for_particles )  THEN
                    lagr_timescale(n) = km(kp,jp,ip) / MAX( e(kp,jp,ip), 1.0E-20_wp )
                    RL             = EXP( -1.0_wp * dt_3d / lagr_timescale(n) )
                    sigma          = SQRT( e(kp,jp,ip) )

                    rg1 = random_gauss( iran_part, 5.0_wp ) - 1.0_wp
                    rg2 = random_gauss( iran_part, 5.0_wp ) - 1.0_wp
                    rg3 = random_gauss( iran_part, 5.0_wp ) - 1.0_wp

                    particles(n)%rvar1 = RL * particles(n)%rvar1 +                &
                                         SQRT( 1.0_wp - RL**2 ) * sigma * rg1
                    particles(n)%rvar2 = RL * particles(n)%rvar2 +                &
                                         SQRT( 1.0_wp - RL**2 ) * sigma * rg2
                    particles(n)%rvar3 = RL * particles(n)%rvar3 +                &
                                         SQRT( 1.0_wp - RL**2 ) * sigma * rg3

                    particles(n)%speed_x = u_int(n) + particles(n)%rvar1
                    particles(n)%speed_y = v_int(n) + particles(n)%rvar2
                    particles(n)%speed_z = w_int(n) + particles(n)%rvar3 - w_s
                ELSE
                    particles(n)%speed_x = u_int(n)
                    particles(n)%speed_y = v_int(n)
                    particles(n)%speed_z = w_int(n) - w_s
                ENDIF

             ELSE

                IF ( use_sgs_for_particles )  THEN
                   exp_arg  = particle_groups(particles(n)%group)%exp_arg
                   exp_term = EXP( -exp_arg * dt_particle(n) )
                ELSE
                   exp_arg  = particle_groups(particles(n)%group)%exp_arg
                   exp_term = particle_groups(particles(n)%group)%exp_term
                ENDIF
                particles(n)%speed_x = particles(n)%speed_x * exp_term +             &
                                       u_int(n) * ( 1.0_wp - exp_term )
                particles(n)%speed_y = particles(n)%speed_y * exp_term +             &
                                       v_int(n) * ( 1.0_wp - exp_term )
                particles(n)%speed_z = particles(n)%speed_z * exp_term +             &
                                       ( w_int(n) - ( 1.0_wp - dens_ratio(n) ) * g / &
                                       exp_arg ) * ( 1.0_wp - exp_term )
             ENDIF
          ENDDO
       ENDDO

    ENDIF

!
!-- Store the old age of the particle ( needed to prevent that a
!-- particle crosses several PEs during one timestep, and for the
!-- evaluation of the subgrid particle velocity fluctuations )
    particles(1:number_of_particles)%age_m = particles(1:number_of_particles)%age

    DO  nb = 0, 7
       DO  n = start_index(nb), end_index(nb)
!
!--       Increment the particle age and the total time that the particle
!--       has advanced within the particle timestep procedure
          particles(n)%age    = particles(n)%age    + dt_particle(n)
          particles(n)%dt_sum = particles(n)%dt_sum + dt_particle(n)

!
!--       Check whether there is still a particle that has not yet completed 
!--       the total LES timestep
          IF ( ( dt_3d - particles(n)%dt_sum ) > 1E-8_wp )  THEN
             dt_3d_reached_l = .FALSE.
          ENDIF

       ENDDO
    ENDDO

    CALL cpu_log( log_point_s(44), 'lpm_advec', 'pause' )


 END SUBROUTINE lpm_advec

! Description:
! ------------
!> Calculation of subgrid-scale particle speed using the stochastic model 
!> of Weil et al. (2004, JAS, 61, 2877-2887).
!------------------------------------------------------------------------------! 
 SUBROUTINE weil_stochastic_eq( v_sgs, fs_n, e_n, dedxi_n, dedt_n, diss_n,     &
                                dt_n, rg_n, fac )

    USE kinds

    USE particle_attributes,                                                   &
        ONLY:  c_0, sgs_wf_part

    IMPLICIT NONE

    REAL(wp) ::  a1      !< dummy argument
    REAL(wp) ::  dedt_n  !< time derivative of TKE at particle position 
    REAL(wp) ::  dedxi_n !< horizontal derivative of TKE at particle position
    REAL(wp) ::  diss_n  !< dissipation at particle position 
    REAL(wp) ::  dt_n    !< particle timestep
    REAL(wp) ::  e_n     !< TKE at particle position
    REAL(wp) ::  fac     !< flag to identify adjacent topography
    REAL(wp) ::  fs_n    !< weighting factor to prevent that subgrid-scale particle speed becomes too large
    REAL(wp) ::  sgs_w   !< constant (1/3)
    REAL(wp) ::  rg_n    !< random number
    REAL(wp) ::  term1   !< memory term
    REAL(wp) ::  term2   !< drift correction term
    REAL(wp) ::  term3   !< random term
    REAL(wp) ::  v_sgs   !< subgrid-scale velocity component 

!-- At first, limit TKE to a small non-zero number, in order to prevent 
!-- the occurrence of extremely large SGS-velocities in case TKE is zero, 
!-- (could occur at the simulation begin).
    e_n = MAX( e_n, 1E-20_wp )
!
!-- Please note, terms 1 and 2 (drift and memory term, respectively) are 
!-- multiplied by a flag to switch of both terms near topography. 
!-- This is necessary, as both terms may cause a subgrid-scale velocity build up 
!-- if particles are trapped in regions with very small TKE, e.g. in narrow street 
!-- canyons resolved by only a few grid points. Hence, term 1 and term 2 are 
!-- disabled if one of the adjacent grid points belongs to topography. 
!-- Moreover, in this case, the  previous subgrid-scale component is also set
!-- to zero.

    a1 = fs_n * c_0 * diss_n
!
!-- Memory term
    term1 = - a1 * v_sgs * dt_n / ( 4.0_wp * sgs_wf_part * e_n + 1E-20_wp )    &
                 * fac
!
!-- Drift correction term
    term2 = ( ( dedt_n * v_sgs / e_n ) + dedxi_n ) * 0.5_wp * dt_n              &
                 * fac
!
!-- Random term
    term3 = SQRT( MAX( a1, 1E-20 ) ) * ( rg_n - 1.0_wp ) * SQRT( dt_n )
!
!-- In cese one of the adjacent grid-boxes belongs to topograhy, the previous
!-- subgrid-scale velocity component is set to zero, in order to prevent a 
!-- velocity build-up. 
!-- This case, set also previous subgrid-scale component to zero. 
    v_sgs = v_sgs * fac + term1 + term2 + term3

 END SUBROUTINE weil_stochastic_eq
