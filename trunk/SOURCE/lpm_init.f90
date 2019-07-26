!> @file lpm_init.f90
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
! $Id: lpm_init.f90 3065 2018-06-12 07:03:02Z Giersch $
! dz was replaced by dzw or dz(1) to allow for right vertical stretching
! 
! 3049 2018-05-29 13:52:36Z Giersch
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error message revised
! 
! 3039 2018-05-24 13:13:11Z schwenkel
! bugfix for lcm with grid stretching
! 
! 2967 2018-04-13 11:22:08Z raasch
! nesting routine is only called if nesting is switched on
! 
! 2954 2018-04-09 14:35:46Z schwenkel
! Bugfix for particle initialization in case of ocean
! 
! 2801 2018-02-14 16:01:55Z thiele
! Introduce particle transfer in nested models.
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2701 2017-12-15 15:40:50Z suehring
! Changes from last commit documented
! 
! 2698 2017-12-14 18:46:24Z suehring
! Grid indices passed to lpm_boundary_conds. (responsible Philipp Thiele)
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2628 2017-11-20 12:40:38Z schwenkel
! Enabled particle advection with grid stretching.
! 
! 2608 2017-11-13 14:04:26Z schwenkel
! Calculation of magnus equation in external module (diagnostic_quantities_mod).
! 
! 2606 2017-11-10 10:36:31Z schwenkel
! Changed particle box locations: center of particle box now coincides 
! with scalar grid point of same index.
! Renamed module and subroutines: lpm_pack_arrays_mod -> lpm_pack_and_sort_mod
! lpm_pack_all_arrays -> lpm_sort_in_subboxes, lpm_pack_arrays -> lpm_pack
! lpm_sort -> lpm_sort_timeloop_done
! 
! 2375 2017-08-29 14:10:28Z schwenkel
! Initialization of chemical aerosol composition
! 
! 2346 2017-08-09 16:39:17Z suehring
! Bugfix, correct determination of topography top index
! 
! 2318 2017-07-20 17:27:44Z suehring
! Get topography top index via Function call 
!
! 2317 2017-07-20 17:27:19Z suehring
! Extended particle data type. Aerosol initialization improved.
!
! 2305 2017-07-06 11:18:47Z hoffmann
! Improved calculation of particle IDs.
!
! 2274 2017-06-09 13:27:48Z Giersch
!  Changed error messages
!
! 2265 2017-06-08 16:58:28Z schwenkel
! Unused variables removed.
!
! 2263 2017-06-08 14:59:01Z schwenkel
! Implemented splitting and merging algorithm
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments according to new topography realization
!
!
! 2223 2017-05-15 16:38:09Z suehring
! Add check for particle release at model top
!
! 2182 2017-03-17 14:27:40Z schwenkel
! Added parameters for simplified particle initialization.
!
! 2122 2017-01-18 12:22:54Z hoffmann
! Improved initialization of equilibrium aerosol radii
! Calculation of particle ID
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 2016-06-09 16:25:25Z suehring
! Bugfix in determining initial particle height and grid index in case of
! seed_follows_topography.
! Bugfix concerning random positions, ensure that particles do not move more
! than one grid length.
! Bugfix logarithmic interpolation.
! Initial setting of sgs_wf_part.
!
! 1890 2016-04-22 08:52:11Z hoffmann
! Initialization of aerosol equilibrium radius not possible in supersaturated
! environments. Therefore, a maximum supersaturation of -1 % is assumed during
! initialization.
!
! 1873 2016-04-18 14:50:06Z maronga
! Module renamed (removed _mod
!
! 1871 2016-04-15 11:46:09Z hoffmann
! Initialization of aerosols added.
!
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
!
! 1831 2016-04-07 13:15:51Z hoffmann
! curvature_solution_effects moved to particle_attributes
!
! 1822 2016-04-07 07:49:42Z hoffmann
! Unused variables removed.
!
! 1783 2016-03-06 18:36:17Z raasch
! netcdf module added
!
! 1725 2015-11-17 13:01:51Z hoffmann
! Bugfix: Processor-dependent seed for random function is generated before it is
! used.
!
! 1691 2015-10-26 16:17:44Z maronga
! Renamed prandtl_layer to constant_flux_layer.
!
! 1685 2015-10-08 07:32:13Z raasch
! bugfix concerning vertical index offset in case of ocean
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1575 2015-03-27 09:56:27Z raasch
! initial vertical particle position is allowed to follow the topography
!
! 1359 2014-04-11 17:15:14Z hoffmann
! New particle structure integrated.
! Kind definition added to all floating point numbers.
! lpm_init changed form a subroutine to a module.
!
! 1327 2014-03-21 11:00:16Z raasch
! -netcdf_output
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL functions provided with KIND-attribute
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
! bugfix: #if defined( __parallel ) added
!
! 1314 2014-03-14 18:25:17Z suehring
! Vertical logarithmic interpolation of horizontal particle speed for particles
! between roughness height and first vertical grid level.
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 849 2012-03-15 10:35:09Z raasch
! routine renamed: init_particles -> lpm_init
! de_dx, de_dy, de_dz are allocated here (instead of automatic arrays in
! advec_particles),
! sort_particles renamed lpm_sort_arrays, user_init_particles renamed lpm_init
!
! 828 2012-02-21 12:00:36Z raasch
! call of init_kernels, particle feature color renamed class
!
! 824 2012-02-17 09:09:57Z raasch
! particle attributes speed_x|y|z_sgs renamed rvar1|2|3,
! array particles implemented as pointer
!
! 667 2010-12-23 12:06:00Z suehring/gryschka
! nxl-1, nxr+1, nys-1, nyn+1 replaced by nxlg, nxrg, nysg, nyng for allocation
! of arrays.
!
! Revision 1.1  1999/11/25 16:22:38  raasch
! Initial revision
!
!
! Description:
! ------------
!> This routine initializes a set of particles and their attributes (position,
!> radius, ..) which are used by the Lagrangian particle model (see lpm).
!------------------------------------------------------------------------------!
 MODULE lpm_init_mod

    USE, INTRINSIC ::  ISO_C_BINDING

    USE arrays_3d,                                                             &
        ONLY:  de_dx, de_dy, de_dz, dzw, zu, zw

    USE control_parameters,                                                    &
        ONLY:  cloud_droplets, constant_flux_layer, current_timestep_number,   &
               dt_3d, dz, initializing_actions, message_string, ocean,         &
               simulated_time

    USE grid_variables,                                                        &
        ONLY:  ddx, dx, ddy, dy

    USE indices,                                                               &
        ONLY:  nx, nxl, nxlg, nxrg, nxr, ny, nyn, nys, nyng, nysg, nz, nzb,    &
               nzt, wall_flags_0

    USE kinds

    USE lpm_collision_kernels_mod,                                             &
        ONLY:  init_kernels

    USE netcdf_interface,                                                      &
        ONLY:  netcdf_data_format

    USE particle_attributes,                                                   &
        ONLY:   alloc_factor, bc_par_b, bc_par_lr, bc_par_ns, bc_par_t,        &
                block_offset, block_offset_def, collision_kernel,              &
                curvature_solution_effects, density_ratio, grid_particles,     &
                isf,i_splitting_mode, initial_weighting_factor, ibc_par_b,     &
                ibc_par_lr, ibc_par_ns, ibc_par_t, iran_part, log_z_z0,        &
                max_number_of_particle_groups, min_nr_particle,                &
                number_concentration,                                          &
                number_particles_per_gridbox,  number_of_particles,            &
                number_of_particle_groups, number_of_sublayers,                &
                offset_ocean_nzt, offset_ocean_nzt_m1,                         &
                particles, particle_advection_start, particle_groups,          &
                particle_groups_type, particles_per_point,                     &
                particle_type, pdx, pdy, pdz,  prt_count, psb, psl, psn, psr,  &
                pss, pst, radius, random_start_position,                       &
                read_particles_from_restartfile, seed_follows_topography,      &
                sgs_wf_part, sort_count, splitting_function, splitting_mode,   &
                total_number_of_particles, use_sgs_for_particles,              &
                write_particle_statistics, zero_particle, z0_av_global

    USE pegrid

    USE random_function_mod,                                                   &
        ONLY:  random_function

    USE surface_mod,                                                           &
        ONLY:  get_topography_top_index_ji, surf_def_h, surf_lsm_h, surf_usm_h

    USE pmc_particle_interface,                                                &
        ONLY:  pmcp_g_init

    IMPLICIT NONE

    PRIVATE

    INTEGER(iwp), PARAMETER         :: PHASE_INIT    = 1  !<
    INTEGER(iwp), PARAMETER, PUBLIC :: PHASE_RELEASE = 2  !<

    INTERFACE lpm_init
       MODULE PROCEDURE lpm_init
    END INTERFACE lpm_init

    INTERFACE lpm_create_particle
       MODULE PROCEDURE lpm_create_particle
    END INTERFACE lpm_create_particle

    PUBLIC lpm_init, lpm_create_particle

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_init

    USE lpm_collision_kernels_mod,                                             &
        ONLY:  init_kernels

    USE pmc_interface,                                                         &
        ONLY: nested_run

    IMPLICIT NONE

    INTEGER(iwp) ::  i                           !<
    INTEGER(iwp) ::  j                           !<
    INTEGER(iwp) ::  k                           !<

    REAL(wp) ::  div                             !<
    REAL(wp) ::  height_int                      !<
    REAL(wp) ::  height_p                        !<
    REAL(wp) ::  z_p                             !<
    REAL(wp) ::  z0_av_local                     !<


!
!-- In case of oceans runs, the vertical index calculations need an offset,
!-- because otherwise the k indices will become negative
    IF ( ocean )  THEN
       offset_ocean_nzt    = nzt
       offset_ocean_nzt_m1 = nzt - 1
    ENDIF

!
!-- Define block offsets for dividing a gridcell in 8 sub cells
!-- See documentation for List of subgrid boxes
!-- See pack_and_sort in lpm_pack_arrays.f90 for assignment of the subgrid boxes
    block_offset(0) = block_offset_def ( 0, 0, 0)
    block_offset(1) = block_offset_def ( 0, 0,-1)
    block_offset(2) = block_offset_def ( 0,-1, 0)
    block_offset(3) = block_offset_def ( 0,-1,-1)
    block_offset(4) = block_offset_def (-1, 0, 0)
    block_offset(5) = block_offset_def (-1, 0,-1)
    block_offset(6) = block_offset_def (-1,-1, 0)
    block_offset(7) = block_offset_def (-1,-1,-1)
!
!-- Check the number of particle groups.
    IF ( number_of_particle_groups > max_number_of_particle_groups )  THEN
       WRITE( message_string, * ) 'max_number_of_particle_groups =',           &
                                  max_number_of_particle_groups ,              &
                                  '&number_of_particle_groups reset to ',      &
                                  max_number_of_particle_groups
       CALL message( 'lpm_init', 'PA0213', 0, 1, 0, 6, 0 )
       number_of_particle_groups = max_number_of_particle_groups
    ENDIF
!
!-- Check if downward-facing walls exist. This case, reflection boundary
!-- conditions (as well as subgrid-scale velocities) may do not work
!-- propably (not realized so far).
    IF ( surf_def_h(1)%ns >= 1 )  THEN
       WRITE( message_string, * ) 'Overhanging topography do not work '//      &
                                  'with particles'
       CALL message( 'lpm_init', 'PA0212', 0, 1, 0, 6, 0 )

    ENDIF

!
!-- Set default start positions, if necessary
    IF ( psl(1) == 9999999.9_wp )  psl(1) = 0.0_wp
    IF ( psr(1) == 9999999.9_wp )  psr(1) = ( nx +1 ) * dx
    IF ( pss(1) == 9999999.9_wp )  pss(1) = 0.0_wp
    IF ( psn(1) == 9999999.9_wp )  psn(1) = ( ny +1 ) * dy
    IF ( psb(1) == 9999999.9_wp )  psb(1) = zu(nz/2)
    IF ( pst(1) == 9999999.9_wp )  pst(1) = psb(1)

    IF ( pdx(1) == 9999999.9_wp  .OR.  pdx(1) == 0.0_wp )  pdx(1) = dx
    IF ( pdy(1) == 9999999.9_wp  .OR.  pdy(1) == 0.0_wp )  pdy(1) = dy
    IF ( pdz(1) == 9999999.9_wp  .OR.  pdz(1) == 0.0_wp )  pdz(1) = zu(2) - zu(1)

!
!-- If number_particles_per_gridbox is set, the parametres pdx, pdy and pdz are
!-- calculated diagnostically. Therfore an isotropic distribution is prescribed.
    IF ( number_particles_per_gridbox /= -1 .AND.   &
         number_particles_per_gridbox >= 1 )    THEN
       pdx(1) = (( dx * dy * ( zu(2) - zu(1) ) ) /  &
             REAL(number_particles_per_gridbox))**0.3333333_wp
!
!--    Ensure a smooth value (two significant digits) of distance between
!--    particles (pdx, pdy, pdz).
       div = 1000.0_wp
       DO  WHILE ( pdx(1) < div )
          div = div / 10.0_wp
       ENDDO
       pdx(1) = NINT( pdx(1) * 100.0_wp / div ) * div / 100.0_wp
       pdy(1) = pdx(1)
       pdz(1) = pdx(1)

    ENDIF

    DO  j = 2, number_of_particle_groups
       IF ( psl(j) == 9999999.9_wp )  psl(j) = psl(j-1)
       IF ( psr(j) == 9999999.9_wp )  psr(j) = psr(j-1)
       IF ( pss(j) == 9999999.9_wp )  pss(j) = pss(j-1)
       IF ( psn(j) == 9999999.9_wp )  psn(j) = psn(j-1)
       IF ( psb(j) == 9999999.9_wp )  psb(j) = psb(j-1)
       IF ( pst(j) == 9999999.9_wp )  pst(j) = pst(j-1)
       IF ( pdx(j) == 9999999.9_wp  .OR.  pdx(j) == 0.0_wp )  pdx(j) = pdx(j-1)
       IF ( pdy(j) == 9999999.9_wp  .OR.  pdy(j) == 0.0_wp )  pdy(j) = pdy(j-1)
       IF ( pdz(j) == 9999999.9_wp  .OR.  pdz(j) == 0.0_wp )  pdz(j) = pdz(j-1)
    ENDDO

!
!-- Allocate arrays required for calculating particle SGS velocities.
!-- Initialize prefactor required for stoachastic Weil equation.
    IF ( use_sgs_for_particles  .AND.  .NOT. cloud_droplets )  THEN
       ALLOCATE( de_dx(nzb:nzt+1,nysg:nyng,nxlg:nxrg), &
                 de_dy(nzb:nzt+1,nysg:nyng,nxlg:nxrg), &
                 de_dz(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

       sgs_wf_part = 1.0_wp / 3.0_wp
    ENDIF

!
!-- Allocate array required for logarithmic vertical interpolation of
!-- horizontal particle velocities between the surface and the first vertical
!-- grid level. In order to avoid repeated CPU cost-intensive CALLS of
!-- intrinsic FORTRAN procedure LOG(z/z0), LOG(z/z0) is precalculated for
!-- several heights. Splitting into 20 sublayers turned out to be sufficient.
!-- To obtain exact height levels of particles, linear interpolation is applied
!-- (see lpm_advec.f90).
    IF ( TRIM(constant_flux_layer) == 'bottom' )  THEN

       ALLOCATE ( log_z_z0(0:number_of_sublayers) )
       z_p = zu(nzb+1) - zw(nzb)

!
!--    Calculate horizontal mean value of z0 used for logartihmic
!--    interpolation. Note: this is not exact for heterogeneous z0.
!--    However, sensitivity studies showed that the effect is
!--    negligible.
       z0_av_local  = SUM( surf_def_h(0)%z0 ) + SUM( surf_lsm_h%z0 ) +         &
                      SUM( surf_usm_h%z0 )
       z0_av_global = 0.0_wp

#if defined( __parallel )
       CALL MPI_ALLREDUCE(z0_av_local, z0_av_global, 1, MPI_REAL, MPI_SUM, &
                          comm2d, ierr )
#else
       z0_av_global = z0_av_local
#endif

       z0_av_global = z0_av_global  / ( ( ny + 1 ) * ( nx + 1 ) )
!
!--    Horizontal wind speed is zero below and at z0
       log_z_z0(0) = 0.0_wp
!
!--    Calculate vertical depth of the sublayers
       height_int  = ( z_p - z0_av_global ) / REAL( number_of_sublayers, KIND=wp )
!
!--    Precalculate LOG(z/z0)
       height_p    = z0_av_global
       DO  k = 1, number_of_sublayers

          height_p    = height_p + height_int
          log_z_z0(k) = LOG( height_p / z0_av_global )

       ENDDO

    ENDIF

!
!-- Check boundary condition and set internal variables
    SELECT CASE ( bc_par_b )

       CASE ( 'absorb' )
          ibc_par_b = 1

       CASE ( 'reflect' )
          ibc_par_b = 2

       CASE DEFAULT
          WRITE( message_string, * )  'unknown boundary condition ',           &
                                       'bc_par_b = "', TRIM( bc_par_b ), '"'
          CALL message( 'lpm_init', 'PA0217', 1, 2, 0, 6, 0 )

    END SELECT
    SELECT CASE ( bc_par_t )

       CASE ( 'absorb' )
          ibc_par_t = 1

       CASE ( 'reflect' )
          ibc_par_t = 2
          
       CASE ( 'nested' )
          ibc_par_t = 3

       CASE DEFAULT
          WRITE( message_string, * ) 'unknown boundary condition ',            &
                                     'bc_par_t = "', TRIM( bc_par_t ), '"'
          CALL message( 'lpm_init', 'PA0218', 1, 2, 0, 6, 0 )

    END SELECT
    SELECT CASE ( bc_par_lr )

       CASE ( 'cyclic' )
          ibc_par_lr = 0

       CASE ( 'absorb' )
          ibc_par_lr = 1

       CASE ( 'reflect' )
          ibc_par_lr = 2
          
       CASE ( 'nested' )
          ibc_par_lr = 3

       CASE DEFAULT
          WRITE( message_string, * ) 'unknown boundary condition ',   &
                                     'bc_par_lr = "', TRIM( bc_par_lr ), '"'
          CALL message( 'lpm_init', 'PA0219', 1, 2, 0, 6, 0 )

    END SELECT
    SELECT CASE ( bc_par_ns )

       CASE ( 'cyclic' )
          ibc_par_ns = 0

       CASE ( 'absorb' )
          ibc_par_ns = 1

       CASE ( 'reflect' )
          ibc_par_ns = 2
          
       CASE ( 'nested' )
          ibc_par_ns = 3

       CASE DEFAULT
          WRITE( message_string, * ) 'unknown boundary condition ',   &
                                     'bc_par_ns = "', TRIM( bc_par_ns ), '"'
          CALL message( 'lpm_init', 'PA0220', 1, 2, 0, 6, 0 )

    END SELECT
    SELECT CASE ( splitting_mode )

       CASE ( 'const' )
          i_splitting_mode = 1

       CASE ( 'cl_av' )
          i_splitting_mode = 2

       CASE ( 'gb_av' )
          i_splitting_mode = 3

       CASE DEFAULT
          WRITE( message_string, * )  'unknown splitting_mode = "',            &
                                      TRIM( splitting_mode ), '"'
          CALL message( 'lpm_init', 'PA0146', 1, 2, 0, 6, 0 )

    END SELECT
    SELECT CASE ( splitting_function )

       CASE ( 'gamma' )
          isf = 1

       CASE ( 'log' )
          isf = 2

       CASE ( 'exp' )
          isf = 3

       CASE DEFAULT
          WRITE( message_string, * )  'unknown splitting function = "',        &
                                       TRIM( splitting_function ), '"'
          CALL message( 'lpm_init', 'PA0147', 1, 2, 0, 6, 0 )

    END SELECT


!
!-- Initialize collision kernels
    IF ( collision_kernel /= 'none' )  CALL init_kernels

!
!-- For the first model run of a possible job chain initialize the
!-- particles, otherwise read the particle data from restart file.
    IF ( TRIM( initializing_actions ) == 'read_restart_data'  &
         .AND.  read_particles_from_restartfile )  THEN

       CALL lpm_read_restart_file

    ELSE

!
!--    Allocate particle arrays and set attributes of the initial set of
!--    particles, which can be also periodically released at later times.
       ALLOCATE( prt_count(nzb:nzt+1,nysg:nyng,nxlg:nxrg), &
                 grid_particles(nzb+1:nzt,nys:nyn,nxl:nxr) )

       number_of_particles         = 0

       sort_count = 0
       prt_count  = 0

!
!--    initialize counter for particle IDs
       grid_particles%id_counter = 1

!
!--    Initialize all particles with dummy values (otherwise errors may
!--    occur within restart runs). The reason for this is still not clear
!--    and may be presumably caused by errors in the respective user-interface.
       zero_particle = particle_type( 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,  &
                                      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,  &
                                      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,  &
                                      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,  &
                                      0, 0, 0_idp, .FALSE., -1 )

       particle_groups = particle_groups_type( 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp )

!
!--    Set values for the density ratio and radius for all particle
!--    groups, if necessary
       IF ( density_ratio(1) == 9999999.9_wp )  density_ratio(1) = 0.0_wp
       IF ( radius(1)        == 9999999.9_wp )  radius(1) = 0.0_wp
       DO  i = 2, number_of_particle_groups
          IF ( density_ratio(i) == 9999999.9_wp )  THEN
             density_ratio(i) = density_ratio(i-1)
          ENDIF
          IF ( radius(i) == 9999999.9_wp )  radius(i) = radius(i-1)
       ENDDO

       DO  i = 1, number_of_particle_groups
          IF ( density_ratio(i) /= 0.0_wp  .AND.  radius(i) == 0 )  THEN
             WRITE( message_string, * ) 'particle group #', i, ' has a',       &
                                        'density ratio /= 0 but radius = 0'
             CALL message( 'lpm_init', 'PA0215', 1, 2, 0, 6, 0 )
          ENDIF
          particle_groups(i)%density_ratio = density_ratio(i)
          particle_groups(i)%radius        = radius(i)
       ENDDO

!
!--    Set a seed value for the random number generator to be exclusively
!--    used for the particle code. The generated random numbers should be
!--    different on the different PEs.
       iran_part = iran_part + myid

       CALL lpm_create_particle (PHASE_INIT)
!
!--    User modification of initial particles
       CALL user_lpm_init

!
!--    Open file for statistical informations about particle conditions
       IF ( write_particle_statistics )  THEN
          CALL check_open( 80 )
          WRITE ( 80, 8000 )  current_timestep_number, simulated_time,         &
                              number_of_particles
          CALL close_file( 80 )
       ENDIF

    ENDIF

    IF ( nested_run )  CALL pmcp_g_init

!
!-- To avoid programm abort, assign particles array to the local version of
!-- first grid cell
    number_of_particles = prt_count(nzb+1,nys,nxl)
    particles => grid_particles(nzb+1,nys,nxl)%particles(1:number_of_particles)
!
!-- Formats
8000 FORMAT (I6,1X,F7.2,4X,I10,71X,I10)

 END SUBROUTINE lpm_init

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_create_particle (phase)
    
    USE arrays_3d,                                                             &
       ONLY:  zw
    USE lpm_exchange_horiz_mod,                                                &
        ONLY: lpm_exchange_horiz, lpm_move_particle, realloc_particles_array

    USE lpm_pack_and_sort_mod,                                                 &
        ONLY: lpm_sort_in_subboxes

    USE particle_attributes,                                                   &
        ONLY: deleted_particles

    IMPLICIT  NONE

    INTEGER(iwp)               ::  alloc_size  !< relative increase of allocated memory for particles
    INTEGER(iwp)               ::  i           !< loop variable ( particle groups )
    INTEGER(iwp)               ::  ip          !< index variable along x
    INTEGER(iwp)               ::  j           !< loop variable ( particles per point )
    INTEGER(iwp)               ::  jp          !< index variable along y
    INTEGER(iwp)               ::  k           !< index variable along z
    INTEGER(iwp)               ::  k_surf      !< index of surface grid point
    INTEGER(iwp)               ::  kp          !< index variable along z
    INTEGER(iwp)               ::  loop_stride !< loop variable for initialization
    INTEGER(iwp)               ::  n           !< loop variable ( number of particles )
    INTEGER(iwp)               ::  new_size    !< new size of allocated memory for particles

    INTEGER(iwp), INTENT(IN)   ::  phase       !< mode of inititialization

    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  local_count !< start address of new particle
    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  local_start !< start address of new particle

    LOGICAL                    ::  first_stride !< flag for initialization

    REAL(wp)                   ::  pos_x      !< increment for particle position in x
    REAL(wp)                   ::  pos_y      !< increment for particle position in y
    REAL(wp)                   ::  pos_z      !< increment for particle position in z
    REAL(wp)                   ::  rand_contr !< dummy argument for random position

    TYPE(particle_type),TARGET ::  tmp_particle !< temporary particle used for initialization

!
!-- Calculate particle positions and store particle attributes, if
!-- particle is situated on this PE
    DO  loop_stride = 1, 2
       first_stride = (loop_stride == 1)
       IF ( first_stride )   THEN
          local_count = 0           ! count number of particles
       ELSE
          local_count = prt_count   ! Start address of new particles
       ENDIF

!
!--    Calculate initial_weighting_factor diagnostically
       IF ( number_concentration /= -1.0_wp .AND. number_concentration > 0.0_wp ) THEN
          initial_weighting_factor =  number_concentration * 1.0E6_wp *             &
                                      pdx(1) * pdy(1) * pdz(1)
       END IF

       n = 0
       DO  i = 1, number_of_particle_groups

          pos_z = psb(i)

          DO WHILE ( pos_z <= pst(i) )

             IF ( pos_z >= zw(0) .AND.  pos_z < zw(nzt) )  THEN


                pos_y = pss(i)

                DO WHILE ( pos_y <= psn(i) )

                   IF ( pos_y >= nys * dy  .AND.                  &
                        pos_y <  ( nyn + 1 ) * dy  ) THEN

                      pos_x = psl(i)

               xloop: DO WHILE ( pos_x <= psr(i) )

                         IF ( pos_x >= nxl * dx  .AND.            &
                              pos_x <  ( nxr + 1) * dx ) THEN

                            DO  j = 1, particles_per_point


                               n = n + 1
                               tmp_particle%x             = pos_x
                               tmp_particle%y             = pos_y
                               tmp_particle%z             = pos_z
                               tmp_particle%age           = 0.0_wp
                               tmp_particle%age_m         = 0.0_wp
                               tmp_particle%dt_sum        = 0.0_wp
                               tmp_particle%e_m           = 0.0_wp
                               tmp_particle%rvar1         = 0.0_wp
                               tmp_particle%rvar2         = 0.0_wp
                               tmp_particle%rvar3         = 0.0_wp
                               tmp_particle%speed_x       = 0.0_wp
                               tmp_particle%speed_y       = 0.0_wp
                               tmp_particle%speed_z       = 0.0_wp
                               tmp_particle%origin_x      = pos_x
                               tmp_particle%origin_y      = pos_y
                               tmp_particle%origin_z      = pos_z
                               IF ( curvature_solution_effects )  THEN
                                  tmp_particle%aux1      = 0.0_wp    ! dry aerosol radius
                                  tmp_particle%aux2      = dt_3d     ! last Rosenbrock timestep
                               ELSE
                                  tmp_particle%aux1      = 0.0_wp    ! free to use
                                  tmp_particle%aux2      = 0.0_wp    ! free to use
                               ENDIF
                               tmp_particle%radius        = particle_groups(i)%radius
                               tmp_particle%weight_factor = initial_weighting_factor
                               tmp_particle%class         = 1
                               tmp_particle%group         = i
                               tmp_particle%id            = 0_idp
                               tmp_particle%particle_mask = .TRUE.
                               tmp_particle%block_nr      = -1
!
!--                            Determine the grid indices of the particle position
                               ip = tmp_particle%x * ddx
                               jp = tmp_particle%y * ddy
                               kp = tmp_particle%z / dz(1) + 1 + offset_ocean_nzt                               
                               DO WHILE( zw(kp) < tmp_particle%z ) 
                                  kp = kp + 1
                               ENDDO
                               DO WHILE( zw(kp-1) > tmp_particle%z )
                                  kp = kp - 1
                               ENDDO 
!
!--                            Determine surface level. Therefore, check for
!--                            upward-facing wall on w-grid. 
                               k_surf = get_topography_top_index_ji( jp, ip, 'w' )

                               IF ( seed_follows_topography )  THEN
!
!--                               Particle height is given relative to topography
                                  kp = kp + k_surf
                                  tmp_particle%z = tmp_particle%z + zw(k_surf)
!--                               Skip particle release if particle position is
!--                               above model top, or within topography in case
!--                               of overhanging structures.
                                  IF ( kp > nzt  .OR.                          &
                                 .NOT. BTEST( wall_flags_0(kp,jp,ip), 0 ) )  THEN
                                     pos_x = pos_x + pdx(i)
                                     CYCLE xloop
                                  ENDIF
!
!--                            Skip particle release if particle position is
!--                            below surface, or within topography in case
!--                            of overhanging structures.
                               ELSEIF ( .NOT. seed_follows_topography .AND.    &
                                         tmp_particle%z <= zw(k_surf)  .OR.    &
                                        .NOT. BTEST( wall_flags_0(kp,jp,ip), 0 ) )&
                               THEN
                                  pos_x = pos_x + pdx(i)
                                  CYCLE xloop
                               ENDIF

                               local_count(kp,jp,ip) = local_count(kp,jp,ip) + 1

                               IF ( .NOT. first_stride )  THEN
                                  IF ( ip < nxl  .OR.  jp < nys  .OR.  kp < nzb+1 )  THEN
                                     write(6,*) 'xl ',ip,jp,kp,nxl,nys,nzb+1
                                  ENDIF
                                  IF ( ip > nxr  .OR.  jp > nyn  .OR.  kp > nzt )  THEN
                                     write(6,*) 'xu ',ip,jp,kp,nxr,nyn,nzt
                                  ENDIF
                                  grid_particles(kp,jp,ip)%particles(local_count(kp,jp,ip)) = tmp_particle

                               ENDIF
                            ENDDO

                         ENDIF

                         pos_x = pos_x + pdx(i)

                      ENDDO xloop

                   ENDIF

                   pos_y = pos_y + pdy(i)

                ENDDO

             ENDIF

             pos_z = pos_z + pdz(i)

          ENDDO

       ENDDO

       IF ( first_stride )  THEN
          DO  ip = nxl, nxr
             DO  jp = nys, nyn
                DO  kp = nzb+1, nzt
                   IF ( phase == PHASE_INIT )  THEN
                      IF ( local_count(kp,jp,ip) > 0 )  THEN
                         alloc_size = MAX( INT( local_count(kp,jp,ip) *        &
                            ( 1.0_wp + alloc_factor / 100.0_wp ) ),            &
                            min_nr_particle )
                      ELSE
                         alloc_size = min_nr_particle
                      ENDIF
                      ALLOCATE(grid_particles(kp,jp,ip)%particles(1:alloc_size))
                      DO  n = 1, alloc_size
                         grid_particles(kp,jp,ip)%particles(n) = zero_particle
                      ENDDO
                   ELSEIF ( phase == PHASE_RELEASE )  THEN
                      IF ( local_count(kp,jp,ip) > 0 )  THEN
                         new_size   = local_count(kp,jp,ip) + prt_count(kp,jp,ip)
                         alloc_size = MAX( INT( new_size * ( 1.0_wp +          &
                            alloc_factor / 100.0_wp ) ), min_nr_particle )
                         IF( alloc_size > SIZE( grid_particles(kp,jp,ip)%particles) )  THEN
                            CALL realloc_particles_array(ip,jp,kp,alloc_size)
                         ENDIF
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDIF

    ENDDO



    local_start = prt_count+1
    prt_count   = local_count

!
!-- Calculate particle IDs
    DO  ip = nxl, nxr
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt
             number_of_particles = prt_count(kp,jp,ip)
             IF ( number_of_particles <= 0 )  CYCLE
             particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)

             DO  n = local_start(kp,jp,ip), number_of_particles  !only new particles

                particles(n)%id = 10000_idp**3 * grid_particles(kp,jp,ip)%id_counter + &
                                  10000_idp**2 * kp + 10000_idp * jp + ip
!
!--             Count the number of particles that have been released before
                grid_particles(kp,jp,ip)%id_counter =                          &
                                         grid_particles(kp,jp,ip)%id_counter + 1

             ENDDO

          ENDDO
       ENDDO
    ENDDO

!
!-- Initialize aerosol background spectrum
    IF ( curvature_solution_effects )  THEN
       CALL lpm_init_aerosols(local_start)
    ENDIF

!
!-- Add random fluctuation to particle positions.
    IF ( random_start_position )  THEN
       DO  ip = nxl, nxr
          DO  jp = nys, nyn
             DO  kp = nzb+1, nzt
                number_of_particles = prt_count(kp,jp,ip)
                IF ( number_of_particles <= 0 )  CYCLE
                particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
!
!--             Move only new particles. Moreover, limit random fluctuation
!--             in order to prevent that particles move more than one grid box,
!--             which would lead to problems concerning particle exchange
!--             between processors in case pdx/pdy are larger than dx/dy,
!--             respectively.
                DO  n = local_start(kp,jp,ip), number_of_particles
                   IF ( psl(particles(n)%group) /= psr(particles(n)%group) )  THEN
                      rand_contr = ( random_function( iran_part ) - 0.5_wp ) * &
                                     pdx(particles(n)%group)
                      particles(n)%x = particles(n)%x +                        &
                              MERGE( rand_contr, SIGN( dx, rand_contr ),       &
                                     ABS( rand_contr ) < dx                    &
                                   )
                   ENDIF
                   IF ( pss(particles(n)%group) /= psn(particles(n)%group) )  THEN
                      rand_contr = ( random_function( iran_part ) - 0.5_wp ) * &
                                     pdy(particles(n)%group)
                      particles(n)%y = particles(n)%y +                        &
                              MERGE( rand_contr, SIGN( dy, rand_contr ),       &
                                     ABS( rand_contr ) < dy                    &
                                   )
                   ENDIF
                   IF ( psb(particles(n)%group) /= pst(particles(n)%group) )  THEN
                      rand_contr = ( random_function( iran_part ) - 0.5_wp ) * &
                                     pdz(particles(n)%group)
                      particles(n)%z = particles(n)%z +                        &
                              MERGE( rand_contr, SIGN( dzw(kp), rand_contr ),  &
                                     ABS( rand_contr ) < dzw(kp)               &
                                   )
                   ENDIF
                ENDDO
!
!--             Identify particles located outside the model domain and reflect
!--             or absorb them if necessary.
                CALL lpm_boundary_conds( 'bottom/top', i, j, k )
!
!--             Furthermore, remove particles located in topography. Note, as
!--             the particle speed is still zero at this point, wall
!--             reflection boundary conditions will not work in this case.
                particles =>                                                   &
                       grid_particles(kp,jp,ip)%particles(1:number_of_particles)
                DO  n = local_start(kp,jp,ip), number_of_particles
                   i = particles(n)%x * ddx
                   j = particles(n)%y * ddy
                   k = particles(n)%z / dz(1) + 1 + offset_ocean_nzt
                   DO WHILE( zw(k) < particles(n)%z )
                      k = k + 1
                   ENDDO
                   DO WHILE( zw(k-1) > particles(n)%z )
                      k = k - 1
                   ENDDO
!
!--                Check if particle is within topography
                   IF ( .NOT. BTEST( wall_flags_0(k,j,i), 0 ) )  THEN
                      particles(n)%particle_mask = .FALSE.
                      deleted_particles = deleted_particles + 1
                   ENDIF

                ENDDO
             ENDDO
          ENDDO
       ENDDO
!
!--    Exchange particles between grid cells and processors
       CALL lpm_move_particle
       CALL lpm_exchange_horiz

    ENDIF
!
!-- In case of random_start_position, delete particles identified by
!-- lpm_exchange_horiz and lpm_boundary_conds. Then sort particles into blocks,
!-- which is needed for a fast interpolation of the LES fields on the particle
!-- position.
    CALL lpm_sort_in_subboxes

!
!-- Determine the current number of particles
    DO  ip = nxl, nxr
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt
             number_of_particles         = number_of_particles                 &
                                           + prt_count(kp,jp,ip)
          ENDDO
       ENDDO
    ENDDO
!
!-- Calculate the number of particles of the total domain
#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( number_of_particles, total_number_of_particles, 1, &
    MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
    total_number_of_particles = number_of_particles
#endif

    RETURN

 END SUBROUTINE lpm_create_particle

 SUBROUTINE lpm_init_aerosols(local_start)

    USE arrays_3d,                                                             &
        ONLY: hyp, pt, q

    USE cloud_parameters,                                                      &
        ONLY: l_d_rv, molecular_weight_of_solute,                              &
              molecular_weight_of_water, rho_l, r_v, rho_s, vanthoff

    USE constants,                                                             &
        ONLY: pi

    USE diagnostic_quantities_mod,                                             &
        ONLY:  magnus


    USE kinds

    USE particle_attributes,                                                   &
        ONLY: aero_species, aero_type, aero_weight, log_sigma, na, rm

    IMPLICIT NONE

    REAL(wp)  :: afactor            !< curvature effects
    REAL(wp)  :: bfactor            !< solute effects
    REAL(wp)  :: dlogr              !< logarithmic width of radius bin
    REAL(wp)  :: e_a                !< vapor pressure
    REAL(wp)  :: e_s                !< saturation vapor pressure
    REAL(wp)  :: rmin = 0.005e-6_wp !< minimum aerosol radius
    REAL(wp)  :: rmax = 10.0e-6_wp  !< maximum aerosol radius
    REAL(wp)  :: r_mid              !< mean radius of bin
    REAL(wp)  :: r_l                !< left radius of bin
    REAL(wp)  :: r_r                !< right radius of bin
    REAL(wp)  :: sigma              !< surface tension
    REAL(wp)  :: t_int              !< temperature

    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(IN) ::  local_start !<

    INTEGER(iwp)  :: n              !<
    INTEGER(iwp)  :: ip             !<
    INTEGER(iwp)  :: jp             !<
    INTEGER(iwp)  :: kp             !<

!
!-- Set constants for different aerosol species
    IF ( TRIM(aero_species) .EQ. 'nacl' ) THEN
       molecular_weight_of_solute = 0.05844_wp 
       rho_s                      = 2165.0_wp
       vanthoff                   = 2.0_wp
    ELSEIF ( TRIM(aero_species) .EQ. 'c3h4o4' ) THEN
       molecular_weight_of_solute = 0.10406_wp 
       rho_s                      = 1600.0_wp
       vanthoff                   = 1.37_wp
    ELSEIF ( TRIM(aero_species) .EQ. 'nh4o3' ) THEN
       molecular_weight_of_solute = 0.08004_wp 
       rho_s                      = 1720.0_wp
       vanthoff                   = 2.31_wp
    ELSE
       WRITE( message_string, * ) 'unknown aerosol species ',   &
                                'aero_species = "', TRIM( aero_species ), '"'
       CALL message( 'lpm_init', 'PA0470', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- The following typical aerosol spectra are taken from Jaenicke (1993):
!-- Tropospheric aerosols. Published in Aerosol-Cloud-Climate Interactions.
    IF ( TRIM(aero_type) .EQ. 'polar' )  THEN
       na        = (/ 2.17e1, 1.86e-1, 3.04e-4 /) * 1.0E6
       rm        = (/ 0.0689, 0.375, 4.29 /) * 1.0E-6
       log_sigma = (/ 0.245, 0.300, 0.291 /)
    ELSEIF ( TRIM(aero_type) .EQ. 'background' )  THEN
       na        = (/ 1.29e2, 5.97e1, 6.35e1 /) * 1.0E6
       rm        = (/ 0.0036, 0.127, 0.259 /) * 1.0E-6
       log_sigma = (/ 0.645, 0.253, 0.425 /)
    ELSEIF ( TRIM(aero_type) .EQ. 'maritime' )  THEN
       na        = (/ 1.33e2, 6.66e1, 3.06e0 /) * 1.0E6
       rm        = (/ 0.0039, 0.133, 0.29 /) * 1.0E-6
       log_sigma = (/ 0.657, 0.210, 0.396 /)
    ELSEIF ( TRIM(aero_type) .EQ. 'continental' )  THEN
       na        = (/ 3.20e3, 2.90e3, 3.00e-1 /) * 1.0E6
       rm        = (/ 0.01, 0.058, 0.9 /) * 1.0E-6
       log_sigma = (/ 0.161, 0.217, 0.380 /)
    ELSEIF ( TRIM(aero_type) .EQ. 'desert' )  THEN
       na        = (/ 7.26e2, 1.14e3, 1.78e-1 /) * 1.0E6
       rm        = (/ 0.001, 0.0188, 10.8 /) * 1.0E-6
       log_sigma = (/ 0.247, 0.770, 0.438 /)
    ELSEIF ( TRIM(aero_type) .EQ. 'rural' )  THEN
       na        = (/ 6.65e3, 1.47e2, 1.99e3 /) * 1.0E6
       rm        = (/ 0.00739, 0.0269, 0.0419 /) * 1.0E-6
       log_sigma = (/ 0.225, 0.557, 0.266 /)
    ELSEIF ( TRIM(aero_type) .EQ. 'urban' )  THEN
       na        = (/ 9.93e4, 1.11e3, 3.64e4 /) * 1.0E6
       rm        = (/ 0.00651, 0.00714, 0.0248 /) * 1.0E-6
       log_sigma = (/ 0.245, 0.666, 0.337 /)
    ELSEIF ( TRIM(aero_type) .EQ. 'user' )  THEN
       CONTINUE
    ELSE
       WRITE( message_string, * ) 'unknown aerosol type ',   &
                                'aero_type = "', TRIM( aero_type ), '"'
       CALL message( 'lpm_init', 'PA0459', 1, 2, 0, 6, 0 )
    ENDIF

    DO  ip = nxl, nxr
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt

             number_of_particles = prt_count(kp,jp,ip)
             IF ( number_of_particles <= 0 )  CYCLE
             particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)

             dlogr   = ( LOG10(rmax) - LOG10(rmin) ) / ( number_of_particles - local_start(kp,jp,ip) + 1 )
!
!--          Initialize the aerosols with a predefined spectral distribution
!--          of the dry radius (logarithmically increasing bins) and a varying
!--          weighting factor
             DO  n = local_start(kp,jp,ip), number_of_particles  !only new particles

                r_l   = 10.0**( LOG10( rmin ) + (n-1) * dlogr )
                r_r   = 10.0**( LOG10( rmin ) + n * dlogr )
                r_mid = SQRT( r_l * r_r )

                particles(n)%aux1          = r_mid
                particles(n)%weight_factor =                                           &
                   ( na(1) / ( SQRT( 2.0 * pi ) * log_sigma(1) ) *                     &
                     EXP( - LOG10( r_mid / rm(1) )**2 / ( 2.0 * log_sigma(1)**2 ) ) +  &
                     na(2) / ( SQRT( 2.0 * pi ) * log_sigma(2) ) *                     &
                     EXP( - LOG10( r_mid / rm(2) )**2 / ( 2.0 * log_sigma(2)**2 ) ) +  &
                     na(3) / ( SQRT( 2.0 * pi ) * log_sigma(3) ) *                     &
                     EXP( - LOG10( r_mid / rm(3) )**2 / ( 2.0 * log_sigma(3)**2 ) )    &
                   ) * ( LOG10(r_r) - LOG10(r_l) ) * ( dx * dy * dzw(kp) )

!
!--             Multiply weight_factor with the namelist parameter aero_weight
!--             to increase or decrease the number of simulated aerosols
                particles(n)%weight_factor = particles(n)%weight_factor * aero_weight

                IF ( particles(n)%weight_factor - FLOOR(particles(n)%weight_factor,KIND=wp) &
                     .GT. random_function( iran_part ) )  THEN
                   particles(n)%weight_factor = FLOOR(particles(n)%weight_factor,KIND=wp) + 1.0
                ELSE
                   particles(n)%weight_factor = FLOOR(particles(n)%weight_factor,KIND=wp)
                ENDIF
!
!--             Unnecessary particles will be deleted
                IF ( particles(n)%weight_factor .LE. 0.0 )  particles(n)%particle_mask = .FALSE.

             ENDDO
!
!--          Set particle radius to equilibrium radius based on the environmental
!--          supersaturation (Khvorostyanov and Curry, 2007, JGR). This avoids
!--          the sometimes lengthy growth toward their equilibrium radius within
!--          the simulation.
             t_int  = pt(kp,jp,ip) * ( hyp(kp) / 100000.0_wp )**0.286_wp

             e_s = magnus( t_int )
             e_a = q(kp,jp,ip) * hyp(kp) / ( q(kp,jp,ip) + 0.622_wp )

             sigma   = 0.0761_wp - 0.000155_wp * ( t_int - 273.15_wp )
             afactor = 2.0_wp * sigma / ( rho_l * r_v * t_int )

             bfactor = vanthoff * molecular_weight_of_water *    &
                       rho_s / ( molecular_weight_of_solute * rho_l )
!
!--          The formula is only valid for subsaturated environments. For
!--          supersaturations higher than -5 %, the supersaturation is set to -5%.
             IF ( e_a / e_s >= 0.95_wp )  e_a = 0.95_wp * e_s

             DO  n = local_start(kp,jp,ip), number_of_particles  !only new particles
!
!--             For details on this equation, see Eq. (14) of Khvorostyanov and
!--             Curry (2007, JGR)
                particles(n)%radius = bfactor**0.3333333_wp *                  &
                   particles(n)%aux1 / ( 1.0_wp - e_a / e_s )**0.3333333_wp / &
                   ( 1.0_wp + ( afactor / ( 3.0_wp * bfactor**0.3333333_wp *   &
                     particles(n)%aux1 ) ) /                                  &
                     ( 1.0_wp - e_a / e_s )**0.6666666_wp                      &
                   )

             ENDDO

          ENDDO
       ENDDO
    ENDDO

 END SUBROUTINE lpm_init_aerosols

END MODULE lpm_init_mod
