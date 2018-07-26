!> @file mod_particle_attributes.f90
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
! $Id: mod_particle_attributes.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2375 2017-08-29 14:10:28Z schwenkel
! molecular_weight_of_solute, molecular_weight_of_water, vanthoff removed and 
! added in modules. Parameters are also used in bulk-microphysics.
!
! 2312 2017-07-14 20:26:51Z hoffmann
! Aerosol initialization improved.
!
! 2305 2017-07-06 11:18:47Z hoffmann
! Improved calculation of particle IDs.
!
! 2278 2017-06-12 13:08:18Z schwenkel
! Added comments
!
! 2265 2017-06-08 16:58:28Z schwenkel
! Unused variables removed.
!
! 2263 2017-06-08 14:59:01Z schwenkel
! Implemented splitting and merging algorithm
!
! 2183 2017-03-17 14:29:15Z schwenkel
!
! 2182 2017-03-17 14:27:40Z schwenkel
! Added parameters for simplified particle initialization.
!
! 2122 2017-01-18 12:22:54Z hoffmann
! Calculation of particle ID
! Particle attribute dvrp_psize renamed to user: this attribute can be used by
! by the user to store any variable
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1936 2016-06-13 13:37:44Z suehring
! +deallocate_memory, step_dealloc
!
! 1929 2016-06-09 16:25:25Z suehring
! -sgs_wfu_par, sgs_wfv_par, sgs_wfw_par
! + sgs_wf_par
!
! 1871 2016-04-15 11:46:09Z hoffmann
! Initialization of aerosols added.
!
! 1849 2016-04-08 11:33:18Z hoffmann
! bfactor, mass_of_solute, molecular_weight_of_solute, molecular_weight_of_water,
! vanthoff added from modules
!
! 1831 2016-04-07 13:15:51Z hoffmann
! palm_kernel removed, curvature_solution_effects added
!
! 1822 2016-04-07 07:49:42Z hoffmann
! +collision_algorithm, all_or_nothing, average_impact
! Tails removed.
!
! 1727 2015-11-20 07:22:02Z knoop
! Bugfix: Cause of syntax warning gfortran preprocessor removed
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1575 2015-03-27 09:56:27Z raasch
! +seed_follows_topography
!
! 1359 2014-04-11 17:15:14Z hoffmann
! new module containing all particle related variables
! -dt_sort_particles
!
! Description:
! ------------
!> Definition of variables used to compute particle transport
!------------------------------------------------------------------------------!
MODULE particle_attributes


    USE kinds

    CHARACTER(LEN=15) ::  aero_species = 'nacl'                    !< aerosol species
    CHARACTER(LEN=15) ::  aero_type    = 'maritime'                !< aerosol type
    CHARACTER(LEN=15) ::  bc_par_lr    = 'cyclic'                  !< left/right boundary condition
    CHARACTER(LEN=15) ::  bc_par_ns    = 'cyclic'                  !< north/south boundary condition
    CHARACTER(LEN=15) ::  bc_par_b     = 'reflect'                 !< bottom boundary condition
    CHARACTER(LEN=15) ::  bc_par_t     = 'absorb'                  !< top boundary condition
    CHARACTER(LEN=15) ::  collision_kernel   = 'none'              !< collision kernel
    CHARACTER(LEN=5)  ::  splitting_function = 'gamma'             !< function for calculation critical weighting factor
    CHARACTER(LEN=5)  ::  splitting_mode     = 'const'             !< splitting mode

    INTEGER(iwp) ::  deleted_particles = 0                        !< number of deleted particles per time step
    INTEGER(iwp) ::  dissipation_classes = 10                     !< namelist parameter (see documentation)
    INTEGER(iwp) ::  ibc_par_b                                    !< particle bottom boundary condition dummy
    INTEGER(iwp) ::  ibc_par_lr                                   !< particle left/right boundary condition dummy
    INTEGER(iwp) ::  ibc_par_ns                                   !< particle north/south boundary condition dummy
    INTEGER(iwp) ::  ibc_par_t                                    !< particle top boundary condition dummy
    INTEGER(iwp) ::  iran_part = -1234567                         !< number for random generator
    INTEGER(iwp) ::  isf                                          !< dummy for splitting function
    INTEGER(iwp) ::  i_splitting_mode                             !< dummy for splitting mode
    INTEGER(iwp) ::  max_number_particles_per_gridbox = 100       !< namelist parameter (see documentation)
    INTEGER(iwp) ::  merge_drp = 0                                !< number of merged droplets
    INTEGER(iwp) ::  min_nr_particle = 50                         !< namelist parameter (see documentation)
    INTEGER(iwp) ::  new_particles = 0                            !< number of new particles
    INTEGER(iwp) ::  n_max = 100                                  !< number of radii bin for splitting functions
    INTEGER(iwp) ::  number_of_particles = 0                      !< number of particles for each grid box (3d array is saved on prt_count)
    INTEGER(iwp) ::  number_of_particle_groups = 1                !< namelist parameter (see documentation)
    INTEGER(iwp) ::  number_of_sublayers = 20                     !< number of sublayers for particle velocities betwenn surface and first grid level
    INTEGER(iwp) ::  number_particles_per_gridbox = -1            !< namelist parameter (see documentation)
    INTEGER(iwp) ::  offset_ocean_nzt = 0                         !< in case of oceans runs, the vertical index calculations need an offset
    INTEGER(iwp) ::  offset_ocean_nzt_m1 = 0                      !< in case of oceans runs, the vertical index calculations need an offset
    INTEGER(iwp) ::  particles_per_point = 1                      !< namelist parameter (see documentation)
    INTEGER(iwp) ::  radius_classes = 20                          !< namelist parameter (see documentation)
    INTEGER(iwp) ::  sort_count = 0                               !< counter for sorting particles
    INTEGER(iwp) ::  splitting_factor = 2                         !< namelist parameter (see documentation)
    INTEGER(iwp) ::  splitting_factor_max = 5                     !< namelist parameter (see documentation)
    INTEGER(iwp) ::  step_dealloc = 100                           !< namelist parameter (see documentation)
    INTEGER(iwp) ::  sum_merge_drp = 0                            !< sum of merged super droplets
    INTEGER(iwp) ::  sum_new_particles = 0                        !< sum of created particles (in splitting algorithm)
    INTEGER(iwp) ::  total_number_of_particles                    !< total number of particles in the whole model domain
    INTEGER(iwp) ::  trlp_count_sum                               !< parameter for particle exchange of PEs
    INTEGER(iwp) ::  trlp_count_recv_sum                          !< parameter for particle exchange of PEs
    INTEGER(iwp) ::  trrp_count_sum                               !< parameter for particle exchange of PEs
    INTEGER(iwp) ::  trrp_count_recv_sum                          !< parameter for particle exchange of PEs
    INTEGER(iwp) ::  trsp_count_sum                               !< parameter for particle exchange of PEs
    INTEGER(iwp) ::  trsp_count_recv_sum                          !< parameter for particle exchange of PEs
    INTEGER(iwp) ::  trnp_count_sum                               !< parameter for particle exchange of PEs
    INTEGER(iwp) ::  trnp_count_recv_sum                          !< parameter for particle exchange of PEs

    INTEGER(iwp), PARAMETER ::  max_number_of_particle_groups = 10 !< maximum allowed number of particle groups

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  prt_count  !< 3d array of number of particles of every grid box

    LOGICAL ::  curvature_solution_effects = .FALSE.      !< namelist parameter (see documentation)
    LOGICAL ::  deallocate_memory = .TRUE.                !< namelist parameter (see documentation)
    LOGICAL ::  hall_kernel = .FALSE.                     !< flag for collision kernel
    LOGICAL ::  merging = .FALSE.                         !< namelist parameter (see documentation)
    LOGICAL ::  particle_advection = .FALSE.              !< parameter to steer the advection of particles
    LOGICAL ::  random_start_position = .FALSE.           !< namelist parameter (see documentation)
    LOGICAL ::  read_particles_from_restartfile = .TRUE.  !< namelist parameter (see documentation)
    LOGICAL ::  seed_follows_topography = .FALSE.         !< namelist parameter (see documentation)
    LOGICAL ::  splitting = .FALSE.                       !< namelist parameter (see documentation)
    LOGICAL ::  use_kernel_tables = .FALSE.               !< parameter, which turns on the use of precalculated collision kernels
    LOGICAL ::  use_sgs_for_particles = .FALSE.           !< namelist parameter (see documentation)
    LOGICAL ::  wang_kernel = .FALSE.                     !< flag for collision kernel
    LOGICAL ::  write_particle_statistics = .FALSE.       !< namelist parameter (see documentation)

    LOGICAL, DIMENSION(max_number_of_particle_groups) ::                       &
                vertical_particle_advection = .TRUE.              !< Switch on/off vertical particle transport

    REAL(wp) ::  aero_weight = 1.0_wp                      !< namelist parameter (see documentation)
    REAL(wp) ::  alloc_factor = 20.0_wp                    !< namelist parameter (see documentation)
    REAL(wp) ::  c_0 = 3.0_wp                              !< parameter for lagrangian timescale
    REAL(wp) ::  dt_min_part = 0.0002_wp                   !< minimum particle time step when SGS velocities are used (s)
    REAL(wp) ::  dt_prel = 9999999.9_wp                    !< namelist parameter (see documentation)
    REAL(wp) ::  dt_write_particle_data = 9999999.9_wp     !< namelist parameter (see documentation)
    REAL(wp) ::  end_time_prel = 9999999.9_wp              !< namelist parameter (see documentation)
    REAL(wp) ::  initial_weighting_factor = 1.0_wp         !< namelist parameter (see documentation)
    REAL(wp) ::  log_sigma(3) = 1.0_wp                     !< namelist parameter (see documentation)
    REAL(wp) ::  na(3) = 0.0_wp                            !< namelist parameter (see documentation)
    REAL(wp) ::  number_concentration = -1.0_wp            !< namelist parameter (see documentation)
    REAL(wp) ::  particle_advection_start = 0.0_wp         !< namelist parameter (see documentation)
    REAL(wp) ::  radius_merge = 1.0E-7_wp                  !< namelist parameter (see documentation)
    REAL(wp) ::  radius_split = 40.0E-6_wp                 !< namelist parameter (see documentation)
    REAL(wp) ::  rm(3) = 1.0E-6_wp                         !< namelist parameter (see documentation)
    REAL(wp) ::  sgs_wf_part                               !< parameter for sgs
    REAL(wp) ::  time_prel = 0.0_wp                        !< time for particle release
    REAL(wp) ::  time_write_particle_data = 0.0_wp         !< write particle data at current time on file
    REAL(wp) ::  weight_factor_merge = -1.0_wp             !< namelist parameter (see documentation)
    REAL(wp) ::  weight_factor_split = -1.0_wp             !< namelist parameter (see documentation)
    REAL(wp) ::  z0_av_global                              !< horizontal mean value of z0

    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  density_ratio = 9999999.9_wp  !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  pdx = 9999999.9_wp            !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  pdy = 9999999.9_wp            !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  pdz = 9999999.9_wp            !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  psb = 9999999.9_wp            !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  psl = 9999999.9_wp            !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  psn = 9999999.9_wp            !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  psr = 9999999.9_wp            !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  pss = 9999999.9_wp            !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  pst = 9999999.9_wp            !< namelist parameter (see documentation).
    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  radius = 9999999.9_wp         !< namelist parameter (see documentation)

    REAL(wp), DIMENSION(:), ALLOCATABLE     ::  log_z_z0   !< Precalculate LOG(z/z0)


    TYPE particle_type
        REAL(wp)     ::  aux1          !< auxiliary multi-purpose feature
        REAL(wp)     ::  aux2          !< auxiliary multi-purpose feature
        REAL(wp)     ::  radius        !< radius of particle
        REAL(wp)     ::  age           !< age of particle
        REAL(wp)     ::  age_m         !<
        REAL(wp)     ::  dt_sum        !<
        REAL(wp)     ::  e_m           !< interpolated sgs tke
        REAL(wp)     ::  origin_x      !< origin x-position of particle (changed cyclic bc)
        REAL(wp)     ::  origin_y      !< origin y-position of particle (changed cyclic bc)
        REAL(wp)     ::  origin_z      !< origin z-position of particle (changed cyclic bc)
        REAL(wp)     ::  rvar1         !<
        REAL(wp)     ::  rvar2         !<
        REAL(wp)     ::  rvar3         !<
        REAL(wp)     ::  speed_x       !< speed of particle in x
        REAL(wp)     ::  speed_y       !< speed of particle in y
        REAL(wp)     ::  speed_z       !< speed of particle in z
        REAL(wp)     ::  weight_factor !< weighting factor
        REAL(wp)     ::  x             !< x-position
        REAL(wp)     ::  y             !< y-position
        REAL(wp)     ::  z             !< z-position
        INTEGER(iwp) ::  class         !< radius class needed for collision
        INTEGER(iwp) ::  group         !< number of particle group
        INTEGER(idp) ::  id            !< particle ID (64 bit integer)
        LOGICAL      ::  particle_mask !< if this parameter is set to false the particle will be deleted
        INTEGER(iwp) ::  block_nr      !< number for sorting (removable?)
    END TYPE particle_type

    TYPE(particle_type), DIMENSION(:), POINTER ::  particles       !< Particle array for this grid cell
    TYPE(particle_type)                        ::  zero_particle   !< zero particle to avoid weird thinge

    TYPE particle_groups_type
        SEQUENCE
        REAL(wp) ::  density_ratio  !< density ratio of the fluid and the particles
        REAL(wp) ::  radius         !< radius of particle
        REAL(wp) ::  exp_arg        !< exponential term of particle inertia
        REAL(wp) ::  exp_term       !< exponential term of particle inertia
    END TYPE particle_groups_type

    TYPE(particle_groups_type), DIMENSION(max_number_of_particle_groups) ::    &
       particle_groups

    TYPE  grid_particle_def
        INTEGER(iwp), DIMENSION(0:7)               ::  start_index        !< start particle index for current block
        INTEGER(iwp), DIMENSION(0:7)               ::  end_index          !< end particle index for current block
        INTEGER(iwp)                               ::  id_counter         !< particle id counter
        LOGICAL                                    ::  time_loop_done     !< timestep loop for particle advection
        TYPE(particle_type), POINTER, DIMENSION(:) ::  particles          !< Particle array for this grid cell
    END TYPE grid_particle_def

    TYPE(grid_particle_def), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  grid_particles

    TYPE block_offset_def          !<
        INTEGER(iwp) ::  i_off     !<
        INTEGER(iwp) ::  j_off     !<
        INTEGER(iwp) ::  k_off     !<
    END TYPE block_offset_def

    TYPE(block_offset_def), DIMENSION(0:7)         ::  block_offset

    SAVE


END MODULE particle_attributes
