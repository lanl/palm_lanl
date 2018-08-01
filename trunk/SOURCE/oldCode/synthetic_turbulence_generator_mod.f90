!> @synthetic_turbulence_generator_mod.f90
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
! Copyright 2017 Leibniz Universitaet Hannover
!------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: synthetic_turbulence_generator_mod.f90 3065 2018-06-12 07:03:02Z Giersch $
! Error message related to vertical stretching has been added, dz was replaced
! by dz(1)
! 
! 3051 2018-05-30 17:43:55Z suehring
! Bugfix in calculation of initial Reynolds-stress tensor.
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
! 
! 3044 2018-05-25 10:59:41Z gronemeier
! Add missing variable descriptions
!
! 3038 2018-05-24 10:54:00Z gronemeier
! updated variable descriptions
!
! 2967 2018-04-13 11:22:08Z raasch
! bugfix: missing parallel cpp-directives added
!
! 2946 2018-04-04 17:01:23Z suehring
! Remove unused module load
!
! 2945 2018-04-04 16:27:14Z suehring
! - Bugfix in parallelization of synthetic turbulence generator in case the
!   z-dimension is not an integral divisor of the number of processors along
!   the x- and y-dimension
! - Revision in control mimic in case of RAN-LES nesting
!
! 2938 2018-03-27 15:52:42Z suehring
! Apply turbulence generator at all non-cyclic lateral boundaries in case of
! realistic Inifor large-scale forcing or RANS-LES nesting
!
! 2936 2018-03-27 14:49:27Z suehring
! variable named found has been introduced for checking if restart data was found,
! reading of restart strings has been moved completely to read_restart_data_mod,
! redundant skipping function has been removed, stg_read/write_restart_data
! have been renamed to stg_r/wrd_global, stg_rrd_global is called in
! read_restart_data_mod now, flag syn_turb_gen_prerun and marker *** end stg
! *** have been removed (Giersch), strings and their respective lengths are
! written out and read now in case of restart runs to get rid of prescribed
! character lengths (Giersch), CASE DEFAULT was added if restart data is read
!
! 2841 2018-02-27 15:02:57Z suehring
! Bugfix: wrong placement of include 'mpif.h' corrected
!
! 2836 2018-02-26 13:40:05Z Giersch
! The variables synthetic_turbulence_generator and
! use_synthetic_turbulence_generator have been abbreviated + syn_turb_gen_prerun
! flag is used to define if module related parameters were outputted as restart
! data
!
! 2716 2017-12-29 16:35:59Z kanani
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2669 2017-12-06 16:03:27Z raasch
! unit number for file containing turbulence generator data changed to 90
! bugfix: preprocessor directives added for MPI specific code
!
! 2576 2017-10-24 13:49:46Z Giersch
! Definition of a new function called stg_skip_global to skip module
! parameters during reading restart data
!
! 2563 2017-10-19 15:36:10Z Giersch
! stg_read_restart_data is called in stg_parin in the case of a restart run
!
! 2259 2017-06-08 09:09:11Z gronemeier
! Initial revision
!
!
!
! Authors:
! --------
! @author Tobias Gronemeier, Atsushi Inagaki, Micha Gryschka, Christoph Knigge
!
!
! Description:
! ------------
!> The module generates turbulence at the inflow boundary based on a method by
!> Xie and Castro (2008) utilizing a Lund rotation (Lund, 1998) and a mass-flux
!> correction by Kim et al. (2013).
!> The turbulence is correlated based on length scales in y- and z-direction and
!> a time scale for each velocity component. The profiles of length and time
!> scales, mean u, v, w, e and pt, and all components of the Reynolds stress
!> tensor are read from file STG_PROFILES.
!>
!> @todo test restart
!>       enable cyclic_fill
!>       implement turbulence generation for e and pt
!> @todo Input of height-constant length scales via namelist
!> @note <Enter notes on the module>
!> @bug  Height information from input file is not used. Profiles from input
!>       must match with current PALM grid.
!>       Transformation of length scales to number of gridpoints does not
!>       consider grid stretching.
!>       In case of restart, velocity seeds differ from precursor run if a11,
!>       a22, or a33 are zero.
!------------------------------------------------------------------------------!
 MODULE synthetic_turbulence_generator_mod


    USE arrays_3d,                                                             &
        ONLY:  mean_inflow_profiles, u, v, w

    USE constants,                                                             &
        ONLY:  pi

    USE control_parameters,                                                    &
        ONLY:  initializing_actions, message_string, syn_turb_gen

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s

    USE indices,                                                               &
        ONLY:  nbgp, nzb, nzt, nxl, nxlg, nxr, nxrg, nys, nyn, nyng, nysg

    USE kinds

#if defined( __parallel )  &&  !defined( __mpifh )
    USE MPI
#endif

    USE pegrid,                                                                &
        ONLY:  comm1dx, comm1dy, comm2d, ierr, myidx, myidy, pdims

    USE transpose_indices,                                                     &
        ONLY: nzb_x, nzt_x


    IMPLICIT NONE

#if defined( __parallel )  &&  defined( __mpifh )
    INCLUDE "mpif.h"
#endif

    LOGICAL :: velocity_seed_initialized = .FALSE.  !< true after first call of stg_main
    LOGICAL :: use_syn_turb_gen = .FALSE.           !< switch to use synthetic turbulence generator

    INTEGER(iwp) ::  id_stg_left        !< left lateral boundary core id in case of turbulence generator
    INTEGER(iwp) ::  id_stg_north       !< north lateral boundary core id in case of turbulence generator
    INTEGER(iwp) ::  id_stg_right       !< right lateral boundary core id in case of turbulence generator
    INTEGER(iwp) ::  id_stg_south       !< south lateral boundary core id in case of turbulence generator
    INTEGER(iwp) ::  stg_type_xz        !< MPI type for full z range
    INTEGER(iwp) ::  stg_type_xz_small  !< MPI type for small z range
    INTEGER(iwp) ::  stg_type_yz        !< MPI type for full z range
    INTEGER(iwp) ::  stg_type_yz_small  !< MPI type for small z range
    INTEGER(iwp) ::  merg               !< maximum length scale (in gp)
    INTEGER(iwp) ::  mergp              !< merg + nbgp
    INTEGER(iwp) ::  nzb_x_stg          !< lower bound of z coordinate (required for transposing z on PEs along x)
    INTEGER(iwp) ::  nzt_x_stg          !< upper bound of z coordinate (required for transposing z on PEs along x)
    INTEGER(iwp) ::  nzb_y_stg          !< lower bound of z coordinate (required for transposing z on PEs along y)
    INTEGER(iwp) ::  nzt_y_stg          !< upper bound of z coordinate (required for transposing z on PEs along y)

    REAL(wp) :: mc_factor    !< mass flux correction factor

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: displs_xz      !< displacement for MPI_GATHERV
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: recv_count_xz  !< receive count for MPI_GATHERV
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: displs_yz      !< displacement for MPI_GATHERV
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: recv_count_yz  !< receive count for MPI_GATHERV
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: nux            !< length scale of u in x direction (in gp)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: nuy            !< length scale of u in y direction (in gp)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: nuz            !< length scale of u in z direction (in gp)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: nvx            !< length scale of v in x direction (in gp)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: nvy            !< length scale of v in y direction (in gp)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: nvz            !< length scale of v in z direction (in gp)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: nwx            !< length scale of w in x direction (in gp)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: nwy            !< length scale of w in y direction (in gp)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: nwz            !< length scale of w in z direction (in gp)

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: seed        !< seed of random number for rn-generator

    REAL(wp), DIMENSION(:), ALLOCATABLE :: a11             !< coefficient for Lund rotation
    REAL(wp), DIMENSION(:), ALLOCATABLE :: a21             !< coefficient for Lund rotation
    REAL(wp), DIMENSION(:), ALLOCATABLE :: a22             !< coefficient for Lund rotation
    REAL(wp), DIMENSION(:), ALLOCATABLE :: a31             !< coefficient for Lund rotation
    REAL(wp), DIMENSION(:), ALLOCATABLE :: a32             !< coefficient for Lund rotation
    REAL(wp), DIMENSION(:), ALLOCATABLE :: a33             !< coefficient for Lund rotation
    REAL(wp), DIMENSION(:), ALLOCATABLE :: tu              !< Lagrangian time scale of u
    REAL(wp), DIMENSION(:), ALLOCATABLE :: tv              !< Lagrangian time scale of v
    REAL(wp), DIMENSION(:), ALLOCATABLE :: tw              !< Lagrangian time scale of w

    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: bux           !< filter function for u in x direction
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: buy           !< filter function for u in y direction
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: buz           !< filter function for u in z direction
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: bvx           !< filter function for v in x direction
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: bvy           !< filter function for v in y direction
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: bvz           !< filter function for v in z direction
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: bwx           !< filter function for w in y direction
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: bwy           !< filter function for w in y direction
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: bwz           !< filter function for w in z direction
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: fu_xz         !< velocity seed for u at xz plane
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: fuo_xz        !< velocity seed for u at xz plane with new random number
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: fu_yz         !< velocity seed for u at yz plane
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: fuo_yz        !< velocity seed for u at yz plane with new random number
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: fv_xz         !< velocity seed for v at xz plane
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: fvo_xz        !< velocity seed for v at xz plane with new random number
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: fv_yz         !< velocity seed for v at yz plane
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: fvo_yz        !< velocity seed for v at yz plane with new random number
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: fw_xz         !< velocity seed for w at xz plane
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: fwo_xz        !< velocity seed for w at xz plane with new random number
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: fw_yz         !< velocity seed for w at yz plane
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: fwo_yz        !< velocity seed for w at yz plane with new random number


!
!-- PALM interfaces:
!-- Input parameter checks to be done in check_parameters
    INTERFACE stg_check_parameters
       MODULE PROCEDURE stg_check_parameters
    END INTERFACE stg_check_parameters

!
!-- Calculate filter functions
    INTERFACE stg_filter_func
       MODULE PROCEDURE stg_filter_func
    END INTERFACE stg_filter_func

!
!-- Generate velocity seeds at south and north domain boundary
    INTERFACE stg_generate_seed_xz
       MODULE PROCEDURE stg_generate_seed_xz
    END INTERFACE stg_generate_seed_xz
!
!-- Generate velocity seeds at left and/or right domain boundary
    INTERFACE stg_generate_seed_yz
       MODULE PROCEDURE stg_generate_seed_yz
    END INTERFACE stg_generate_seed_yz

!
!-- Output of information to the header file
    INTERFACE stg_header
       MODULE PROCEDURE stg_header
    END INTERFACE stg_header

!
!-- Initialization actions
    INTERFACE stg_init
       MODULE PROCEDURE stg_init
    END INTERFACE stg_init

!
!-- Main procedure of synth. turb. gen.
    INTERFACE stg_main
       MODULE PROCEDURE stg_main
    END INTERFACE stg_main

!
!-- Reading of NAMELIST parameters
    INTERFACE stg_parin
       MODULE PROCEDURE stg_parin
    END INTERFACE stg_parin

!
!-- Reading of parameters for restart runs
    INTERFACE stg_rrd_global
       MODULE PROCEDURE stg_rrd_global
    END INTERFACE stg_rrd_global

!
!-- Writing of binary output for restart runs
    INTERFACE stg_wrd_global
       MODULE PROCEDURE stg_wrd_global
    END INTERFACE stg_wrd_global

    SAVE

    PRIVATE

!
!-- Public interfaces
    PUBLIC  stg_check_parameters, stg_header, stg_init, stg_main, stg_parin,   &
            stg_wrd_global, stg_rrd_global

!
!-- Public variables
    PUBLIC  id_stg_left, id_stg_north, id_stg_right, id_stg_south,             &
            use_syn_turb_gen


 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for synthetic turbulence generator
!------------------------------------------------------------------------------!
 SUBROUTINE stg_check_parameters


    USE control_parameters,                                                    &
        ONLY:  bc_lr, bc_ns, forcing, nest_domain, number_stretch_level_start, &
               rans_mode, turbulent_inflow

    USE pmc_interface,                                                         &
        ONLY : rans_mode_parent


    IMPLICIT NONE

    IF ( .NOT. use_syn_turb_gen  .AND.  .NOT. rans_mode  .AND.  forcing )  THEN
       message_string = 'Synthetic turbulence generator has to be applied ' // &
                        'when forcing is used and model operates in LES mode.'
       CALL message( 'stg_check_parameters', 'PA0000', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( .NOT. use_syn_turb_gen  .AND.  nest_domain                            &
         .AND. rans_mode_parent  .AND.  .NOT. rans_mode )  THEN
       message_string = 'Synthetic turbulence generator has to be applied ' // &
                        'when nesting is applied and parent operates in '  //  &
                        'RANS-mode but current child in LES mode.'
       CALL message( 'stg_check_parameters', 'PA0000', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( use_syn_turb_gen )  THEN

       IF ( .NOT. forcing  .AND.  .NOT. nest_domain )  THEN

          IF ( INDEX( initializing_actions, 'set_constant_profiles' ) == 0     &
        .AND.  INDEX( initializing_actions, 'read_restart_data' ) == 0 )  THEN
             message_string = 'Using synthetic turbulence generator ' //       &
                              'requires %initializing_actions = '         //   &
                              '"set_constant_profiles" or "read_restart_data"'
             CALL message( 'stg_check_parameters', 'PA0015', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( bc_lr /= 'dirichlet/radiation' )  THEN
             message_string = 'Using synthetic turbulence generator ' //       &
                              'requires &bc_lr = "dirichlet/radiation"'
             CALL message( 'stg_check_parameters', 'PA0035', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( bc_ns /= 'cyclic' )  THEN
             message_string = 'Using synthetic turbulence generator ' //       &
                              'requires &bc_ns = "cyclic"'
             CALL message( 'stg_check_parameters', 'PA0037', 1, 2, 0, 6, 0 )
          ENDIF

       ENDIF

       IF ( turbulent_inflow )  THEN
          message_string = 'Using synthetic turbulence generator ' //          &
                           'in combination &with turbulent_inflow = .T. '//    &
                              'is not allowed'
          CALL message( 'stg_check_parameters', 'PA0039', 1, 2, 0, 6, 0 )
       ENDIF
       
       IF ( number_stretch_level_start > 0 )  THEN
          message_string = 'Using synthetic turbulence generator ' //          &
                           'in combination with stretching is not allowed'
          CALL message( 'stg_check_parameters', 'PA0420', 1, 2, 0, 6, 0 )
       ENDIF

    ENDIF

 END SUBROUTINE stg_check_parameters


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for synthetic turbulence generator
!------------------------------------------------------------------------------!
 SUBROUTINE stg_header ( io )


    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  io   !< Unit of the output file

!
!-- Write synthetic turbulence generator Header
    WRITE( io, 1 )
    IF ( use_syn_turb_gen )  THEN
       WRITE( io, 2 )
    ELSE
       WRITE( io, 3 )
    ENDIF

1   FORMAT (//' Synthetic turbulence generator information:'/                  &
              ' ------------------------------------------'/)
2   FORMAT ('    synthetic turbulence generator switched on')
3   FORMAT ('    synthetic turbulence generator switched off')

 END SUBROUTINE stg_header


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the synthetic turbulence generator
!------------------------------------------------------------------------------!
 SUBROUTINE stg_init


    USE arrays_3d,                                                             &
        ONLY:  ddzw, u_init, v_init, zu

    USE control_parameters,                                                    &
        ONLY:  coupling_char, dz, e_init, forcing, nest_domain, rans_mode

    USE grid_variables,                                                        &
        ONLY:  ddx, ddy

    USE indices,                                                               &
        ONLY:  nz

    USE pmc_interface,                                                         &
        ONLY : rans_mode_parent


    IMPLICIT NONE

    LOGICAL ::  file_stg_exist = .FALSE. !< flag indication whether parameter file for Reynolds stress and length scales exist

#if defined( __parallel )
    INTEGER(KIND=MPI_ADDRESS_KIND) :: extent !< extent of new MPI type
    INTEGER(KIND=MPI_ADDRESS_KIND) :: tob=0  !< dummy variable
#endif

    INTEGER(iwp) :: i                        !> grid index in x-direction
    INTEGER(iwp) :: j                        !> loop index
    INTEGER(iwp) :: k                        !< index
    INTEGER(iwp) :: newtype                  !< dummy MPI type
    INTEGER(iwp) :: realsize                 !< size of REAL variables
    INTEGER(iwp) :: nseed                    !< dimension of random number seed
    INTEGER(iwp) :: startseed = 1234567890   !< start seed for random number generator

!
!-- Dummy variables used for reading profiles
    REAL(wp) :: d1      !< u profile
    REAL(wp) :: d2      !< v profile
    REAL(wp) :: d3      !< w profile
    REAL(wp) :: d5      !< e profile
    REAL(wp) :: d11     !< vertical interpolation for a11
    REAL(wp) :: d21     !< vertical interpolation for a21
    REAL(wp) :: d22     !< vertical interpolation for a22
    REAL(wp) :: dum_exp !< dummy variable used for exponential vertical decrease of turbulent length scales
    REAL(wp) :: luy     !< length scale for u in y direction
    REAL(wp) :: luz     !< length scale for u in z direction
    REAL(wp) :: lvy     !< length scale for v in y direction
    REAL(wp) :: lvz     !< length scale for v in z direction
    REAL(wp) :: lwy     !< length scale for w in y direction
    REAL(wp) :: lwz     !< length scale for w in z direction
    REAL(wp) :: nnz     !< increment used to determine processor decomposition of z-axis along x and y direction
    REAL(wp) :: zz      !< height

    REAL(wp) :: length_scale_surface !< typical length scale
    REAL(wp) :: r_ii_0               !< correction factor
    REAL(wp) :: time_scale           !< typical time scale
    REAL(wp) :: length_scale_z       !< typical length scale

    REAL(wp),DIMENSION(nzb:nzt+1) :: r11  !< Reynolds parameter
    REAL(wp),DIMENSION(nzb:nzt+1) :: r21  !< Reynolds parameter
    REAL(wp),DIMENSION(nzb:nzt+1) :: r22  !< Reynolds parameter
    REAL(wp),DIMENSION(nzb:nzt+1) :: r31  !< Reynolds parameter
    REAL(wp),DIMENSION(nzb:nzt+1) :: r32  !< Reynolds parameter
    REAL(wp),DIMENSION(nzb:nzt+1) :: r33  !< Reynolds parameter


#if defined( __parallel )
    CALL MPI_BARRIER( comm2d, ierr )
#endif

    CALL  cpu_log( log_point(57), 'synthetic_turbulence_gen', 'start' )

    IF ( .NOT. ALLOCATED( mean_inflow_profiles ) )                             &
       ALLOCATE( mean_inflow_profiles(nzb:nzt+1,5) )

    ALLOCATE ( a11(nzb:nzt+1), a21(nzb:nzt+1), a22(nzb:nzt+1),                 &
               a31(nzb:nzt+1), a32(nzb:nzt+1), a33(nzb:nzt+1),                 &
               nux(nzb:nzt+1), nuy(nzb:nzt+1), nuz(nzb:nzt+1),                 &
               nvx(nzb:nzt+1), nvy(nzb:nzt+1), nvz(nzb:nzt+1),                 &
               nwx(nzb:nzt+1), nwy(nzb:nzt+1), nwz(nzb:nzt+1),                 &
               tu(nzb:nzt+1),  tv(nzb:nzt+1),  tw(nzb:nzt+1)   )

#if defined( __parallel )
!
!-- Determine processor decomposition of z-axis along x- and y-direction
    nnz = nz / pdims(1)
    nzb_x_stg = 1 + myidx * INT( nnz )
    nzt_x_stg = ( myidx + 1 ) * INT( nnz )

    IF ( MOD( nz , pdims(1) ) /= 0  .AND.  myidx == id_stg_right )             &
       nzt_x_stg = nzt_x_stg + myidx * ( nnz - INT( nnz ) )
!        nzt_x_stg = myidx * nnz + MOD( nz , pdims(1) )

    IF ( forcing  .OR.  ( nest_domain .AND.  rans_mode_parent  .AND.           &
                   .NOT.  rans_mode ) )  THEN
       nnz = nz / pdims(2)
       nzb_y_stg = 1 + myidy * INT( nnz )
       nzt_y_stg = ( myidy + 1 ) * INT( nnz )

       IF ( MOD( nz , pdims(2) ) /= 0  .AND.  myidy == id_stg_north )          &
          nzt_y_stg = nzt_y_stg + myidy * ( nnz - INT( nnz ) )
!           nzt_y_stg = myidy * nnz + MOD( nz , pdims(2) )
    ENDIF

!
!-- Define MPI type used in stg_generate_seed_yz to gather vertical splitted
!-- velocity seeds
    CALL MPI_TYPE_SIZE( MPI_REAL, realsize, ierr )
    extent = 1 * realsize
!
!-- Set-up MPI datatyp to involve all cores for turbulence generation at yz
!-- layer
!-- stg_type_yz: yz-slice with vertical bounds nzb:nzt+1
    CALL MPI_TYPE_CREATE_SUBARRAY( 2, [nzt-nzb+2,nyng-nysg+1],                 &
            [1,nyng-nysg+1], [0,0], MPI_ORDER_FORTRAN, MPI_REAL, newtype, ierr )
    CALL MPI_TYPE_CREATE_RESIZED( newtype, tob, extent, stg_type_yz, ierr )
    CALL MPI_TYPE_COMMIT( stg_type_yz, ierr )
    CALL MPI_TYPE_FREE( newtype, ierr )

    ! stg_type_yz_small: yz-slice with vertical bounds nzb_x_stg:nzt_x_stg+1
    CALL MPI_TYPE_CREATE_SUBARRAY( 2, [nzt_x_stg-nzb_x_stg+2,nyng-nysg+1],     &
            [1,nyng-nysg+1], [0,0], MPI_ORDER_FORTRAN, MPI_REAL, newtype, ierr )
    CALL MPI_TYPE_CREATE_RESIZED( newtype, tob, extent, stg_type_yz_small, ierr )
    CALL MPI_TYPE_COMMIT( stg_type_yz_small, ierr )
    CALL MPI_TYPE_FREE( newtype, ierr )

    ! receive count and displacement for MPI_GATHERV in stg_generate_seed_yz
    ALLOCATE( recv_count_yz(pdims(1)), displs_yz(pdims(1)) )

    recv_count_yz           = nzt_x_stg-nzb_x_stg + 1
    recv_count_yz(pdims(1)) = recv_count_yz(pdims(1)) + 1

    DO  j = 1, pdims(1)
       displs_yz(j) = 0 + (nzt_x_stg-nzb_x_stg+1) * (j-1)
    ENDDO
!
!-- Set-up MPI datatyp to involve all cores for turbulence generation at xz
!-- layer
!-- stg_type_xz: xz-slice with vertical bounds nzb:nzt+1
    IF ( forcing  .OR.  ( nest_domain .AND.  rans_mode_parent  .AND.           &
                   .NOT.  rans_mode ) )  THEN
       CALL MPI_TYPE_CREATE_SUBARRAY( 2, [nzt-nzb+2,nxrg-nxlg+1],              &
               [1,nxrg-nxlg+1], [0,0], MPI_ORDER_FORTRAN, MPI_REAL, newtype, ierr )
       CALL MPI_TYPE_CREATE_RESIZED( newtype, tob, extent, stg_type_xz, ierr )
       CALL MPI_TYPE_COMMIT( stg_type_xz, ierr )
       CALL MPI_TYPE_FREE( newtype, ierr )

       ! stg_type_yz_small: xz-slice with vertical bounds nzb_x_stg:nzt_x_stg+1
       CALL MPI_TYPE_CREATE_SUBARRAY( 2, [nzt_y_stg-nzb_y_stg+2,nxrg-nxlg+1],  &
               [1,nxrg-nxlg+1], [0,0], MPI_ORDER_FORTRAN, MPI_REAL, newtype, ierr )
       CALL MPI_TYPE_CREATE_RESIZED( newtype, tob, extent, stg_type_xz_small, ierr )
       CALL MPI_TYPE_COMMIT( stg_type_xz_small, ierr )
       CALL MPI_TYPE_FREE( newtype, ierr )

       ! receive count and displacement for MPI_GATHERV in stg_generate_seed_yz
       ALLOCATE( recv_count_xz(pdims(2)), displs_xz(pdims(2)) )

       recv_count_xz           = nzt_y_stg-nzb_y_stg + 1
       recv_count_xz(pdims(2)) = recv_count_xz(pdims(2)) + 1

       DO  j = 1, pdims(2)
          displs_xz(j) = 0 + (nzt_y_stg-nzb_y_stg+1) * (j-1)
       ENDDO

    ENDIF

#endif
!
!-- Define seed of random number
    CALL RANDOM_SEED()
    CALL RANDOM_SEED( size=nseed )
    ALLOCATE( seed(1:nseed) )
    DO  j = 1, nseed
       seed(j) = startseed + j
    ENDDO
    CALL RANDOM_SEED(put=seed)

!-- Read inflow profile
!-- Height levels of profiles in input profile is as follows:
!-- zu: luy, luz, tu, lvy, lvz, tv, r11, r21, r22, d1, d2, d5
!-- zw: lwy, lwz, tw, r31, r32, r33, d3
!-- WARNING: zz is not used at the moment
    INQUIRE( FILE = 'STG_PROFILES' // TRIM( coupling_char ),                   &
             EXIST = file_stg_exist  )

    IF ( file_stg_exist )  THEN

       OPEN( 90, FILE='STG_PROFILES'//TRIM( coupling_char ), STATUS='OLD',     &
                      FORM='FORMATTED')
!
!--    Skip header
       READ( 90, * )

       DO  k = nzb, nzt+1
          READ( 90, * ) zz, luy, luz, tu(k), lvy, lvz, tv(k), lwy, lwz, tw(k), &
                        r11(k), r21(k), r22(k), r31(k), r32(k), r33(k),        &
                        d1, d2, d3, d5

!
!--       Convert length scales from meter to number of grid points. Attention:
!--       Does not work if grid stretching is used
          nuy(k) = INT( luy * ddy )
          nuz(k) = INT( luz / dz(1)  )
          nvy(k) = INT( lvy * ddy )
          nvz(k) = INT( lvz / dz(1)  )
          nwy(k) = INT( lwy * ddy )
          nwz(k) = INT( lwz / dz(1)  )
!
!--       Workaround, assume isotropic turbulence
          nwx(k) = nwy(k)
          nvx(k) = nvy(k)
          nux(k) = nuy(k)
!
!--       Save Mean inflow profiles
          IF ( TRIM( initializing_actions ) /= 'read_restart_data' ) THEN
             mean_inflow_profiles(k,1) = d1
             mean_inflow_profiles(k,2) = d2
            !  mean_inflow_profiles(k,4) = d4
             mean_inflow_profiles(k,5) = d5
          ENDIF
       ENDDO

       CLOSE( 90 )

    ELSE
!
!--    Set-up defaul length scales. Assume exponentially decreasing length
!--    scales and isotropic turbulence.
!--    Typical length (time) scales of 100 m (s) should be a good compromise
!--    between all stratrifications. Near-surface variances are fixed to
!--    0.1 m2/s2, vertical fluxes are one order of magnitude smaller.
!--    Vertical fluxes
       length_scale_surface = 100.0_wp
       r_ii_0               = 0.1_wp
       time_scale           = 100.0_wp
       DO  k = nzb+1, nzt+1
          dum_exp        = MERGE( -zu(k) / ( 0.3* zu(nzt) ),                   &
                                  0.0_wp,                                      &
                                  zu(k) > 0.3 * zu(nzt)                        &
                                )
          length_scale_z = length_scale_surface * EXP( dum_exp )

          nux(k) = MAX( INT( length_scale_z * ddx     ), 1 )
          nuy(k) = MAX( INT( length_scale_z * ddy     ), 1 )
          nuz(k) = MAX( INT( length_scale_z * ddzw(k) ), 1 )
          nvx(k) = MAX( INT( length_scale_z * ddx     ), 1 )
          nvy(k) = MAX( INT( length_scale_z * ddy     ), 1 )
          nvz(k) = MAX( INT( length_scale_z * ddzw(k) ), 1 )
          nwx(k) = MAX( INT( length_scale_z * ddx     ), 1 )
          nwy(k) = MAX( INT( length_scale_z * ddy     ), 1 )
          nwz(k) = MAX( INT( length_scale_z * ddzw(k) ), 1 )

          r11(k) = r_ii_0 * EXP( dum_exp )
          r22(k) = r_ii_0 * EXP( dum_exp )
          r33(k) = r_ii_0 * EXP( dum_exp )

          r21(k) = 0.1_wp * r_ii_0 * EXP( dum_exp )
          r31(k) = 0.1_wp * r_ii_0 * EXP( dum_exp )
          r32(k) = 0.1_wp * r_ii_0 * EXP( dum_exp )

          tu(k)  = time_scale
          tv(k)  = time_scale
          tw(k)  = time_scale

       ENDDO
       nux(nzb) = nux(nzb+1)
       nuy(nzb) = nuy(nzb+1)
       nuz(nzb) = nuz(nzb+1)
       nvx(nzb) = nvx(nzb+1)
       nvy(nzb) = nvy(nzb+1)
       nvz(nzb) = nvz(nzb+1)
       nwx(nzb) = nwx(nzb+1)
       nwy(nzb) = nwy(nzb+1)
       nwz(nzb) = nwz(nzb+1)

       r11(nzb) = r11(nzb+1)
       r22(nzb) = r22(nzb+1)
       r33(nzb) = r33(nzb+1)

       r21(nzb) = r11(nzb+1)
       r31(nzb) = r31(nzb+1)
       r32(nzb) = r32(nzb+1)

       tu(nzb)  = time_scale
       tv(nzb)  = time_scale
       tw(nzb)  = time_scale

    ENDIF

!
!-- Assign initial profiles
    IF ( .NOT. forcing  .AND.  .NOT.  nest_domain )  THEN
       u_init = mean_inflow_profiles(:,1)
       v_init = mean_inflow_profiles(:,2)
      !pt_init = mean_inflow_profiles(:,4)
       e_init = MAXVAL( mean_inflow_profiles(:,5) )
    ENDIF
!
!-- Calculate coefficient matrix from Reynolds stress tensor (Lund rotation)
    DO  k = nzb, nzt+1
       IF ( r11(k) > 0.0_wp )  THEN
          a11(k) = SQRT( r11(k) )
          a21(k) = r21(k) / a11(k)
       ELSE
          a11(k) = 0.0_wp
          a21(k) = 0.0_wp
       ENDIF

       a22(k) = r22(k) - a21(k)**2
       IF ( a22(k) > 0.0_wp )  THEN
          a22(k) = SQRT( a22(k) )
       ELSE
          a22(k) = 0.0_wp
       ENDIF

!
!--    a31, a32, a33 must be calculated with interpolated a11, a21, a22 (d11,
!--    d21, d22) because of different vertical grid
       IF ( k .le. nzt )  THEN
          d11 = 0.5_wp * ( r11(k) + r11(k+1) )
          IF ( d11 > 0.0_wp )  THEN
             d11 = SQRT( d11 )
             d21 = ( 0.5_wp * ( r21(k) + r21(k+1) ) ) / d11
             a31(k) = r31(k) / d11
          ELSE
             d21 = 0.0_wp
             a31(k) = 0.0_wp
          ENDIF

          d22 = 0.5_wp * ( r22(k) + r22(k+1) ) - d21 ** 2
          IF ( d22 > 0.0_wp )  THEN
             a32(k) = ( r32(k) - d21 * a31(k) ) / SQRT( d22 )
          ELSE
             a32(k) = 0.0_wp
          ENDIF

          a33(k) = r33(k) - a31(k) ** 2 - a32(k) ** 2
          IF ( a33(k) > 0.0_wp )  THEN
             a33(k) = SQRT( a33(k) )
          ELSE
             a33(k) = 0.0_wp
          ENDIF
       ELSE
          a31(k) = a31(k-1)
          a32(k) = a32(k-1)
          a33(k) = a33(k-1)
       ENDIF

    ENDDO
!
!-- Define the size of the filter functions and allocate them.
    merg = 0

    ! arrays must be large enough to cover the largest length scale
    DO  k = nzb, nzt+1
       j = MAX( ABS(nux(k)), ABS(nuy(k)), ABS(nuz(k)), &
                ABS(nvx(k)), ABS(nvy(k)), ABS(nvz(k)), &
                ABS(nwx(k)), ABS(nwy(k)), ABS(nwz(k))  )
       IF ( j > merg )  merg = j
    ENDDO

    merg  = 2 * merg
    mergp = merg + nbgp

    ALLOCATE ( bux(-merg:merg,nzb:nzt+1),                                      &
               buy(-merg:merg,nzb:nzt+1),                                      &
               buz(-merg:merg,nzb:nzt+1),                                      &
               bvx(-merg:merg,nzb:nzt+1),                                      &
               bvy(-merg:merg,nzb:nzt+1),                                      &
               bvz(-merg:merg,nzb:nzt+1),                                      &
               bwx(-merg:merg,nzb:nzt+1),                                      &
               bwy(-merg:merg,nzb:nzt+1),                                      &
               bwz(-merg:merg,nzb:nzt+1)  )

!
!-- Allocate velocity seeds for turbulence at xz-layer
    ALLOCATE ( fu_xz( nzb:nzt+1,nxlg:nxrg), fuo_xz(nzb:nzt+1,nxlg:nxrg),       &
               fv_xz( nzb:nzt+1,nxlg:nxrg), fvo_xz(nzb:nzt+1,nxlg:nxrg),       &
               fw_xz( nzb:nzt+1,nxlg:nxrg), fwo_xz(nzb:nzt+1,nxlg:nxrg)  )

!
!-- Allocate velocity seeds for turbulence at yz-layer
    ALLOCATE ( fu_yz( nzb:nzt+1,nysg:nyng), fuo_yz(nzb:nzt+1,nysg:nyng),       &
               fv_yz( nzb:nzt+1,nysg:nyng), fvo_yz(nzb:nzt+1,nysg:nyng),       &
               fw_yz( nzb:nzt+1,nysg:nyng), fwo_yz(nzb:nzt+1,nysg:nyng)  )

    fu_xz  = 0.0_wp
    fuo_xz = 0.0_wp
    fv_xz  = 0.0_wp
    fvo_xz = 0.0_wp
    fw_xz  = 0.0_wp
    fwo_xz = 0.0_wp

    fu_yz  = 0.0_wp
    fuo_yz = 0.0_wp
    fv_yz  = 0.0_wp
    fvo_yz = 0.0_wp
    fw_yz  = 0.0_wp
    fwo_yz = 0.0_wp

!
!-- Create filter functions
    CALL stg_filter_func( nux, bux ) !filter ux
    CALL stg_filter_func( nuy, buy ) !filter uy
    CALL stg_filter_func( nuz, buz ) !filter uz
    CALL stg_filter_func( nvx, bvx ) !filter vx
    CALL stg_filter_func( nvy, bvy ) !filter vy
    CALL stg_filter_func( nvz, bvz ) !filter vz
    CALL stg_filter_func( nwx, bwx ) !filter wx
    CALL stg_filter_func( nwy, bwy ) !filter wy
    CALL stg_filter_func( nwz, bwz ) !filter wz

#if defined( __parallel )
    CALL MPI_BARRIER( comm2d, ierr )
#endif

!
!-- In case of restart, calculate velocity seeds fu, fv, fw from former
!   time step.
!   Bug: fu, fv, fw are different in those heights where a11, a22, a33
!        are 0 compared to the prerun. This is mostly for k=nzt+1.
    IF ( TRIM( initializing_actions ) == 'read_restart_data' )  THEN
       IF ( myidx == id_stg_left  .OR.  myidx == id_stg_right )  THEN

          IF ( myidx == id_stg_left  )  i = -1
          IF ( myidx == id_stg_right )  i = nxr+1

          DO  j = nysg, nyng
             DO  k = nzb, nzt+1

                IF  ( a11(k) .NE. 0._wp ) THEN
                   fu_yz(k,j) = ( u(k,j,i) / mc_factor - u_init(k) ) / a11(k)
                ELSE
                   fu_yz(k,j) = 0._wp
                ENDIF

                IF  ( a22(k) .NE. 0._wp ) THEN
                   fv_yz(k,j) = ( v(k,j,i) / mc_factor - a21(k) * fu_yz(k,j) - &
                               v_init(k) ) / a22(k)
                ELSE
                   fv_yz(k,j) = 0._wp
                ENDIF

                IF  ( a33(k) .NE. 0._wp ) THEN
                   fw_yz(k,j) = ( w(k,j,i) / mc_factor - a31(k) * fu_yz(k,j) - &
                               a32(k) * fv_yz(k,j) ) / a33(k)
                ELSE
                   fw_yz = 0._wp
                ENDIF

             ENDDO
          ENDDO
       ENDIF
    ENDIF

    CALL  cpu_log( log_point(57), 'synthetic_turbulence_gen', 'stop' )

 END SUBROUTINE stg_init


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate filter function bxx from length scale nxx following Eg.9 and 10
!> (Xie and Castro, 2008)
!------------------------------------------------------------------------------!
 SUBROUTINE stg_filter_func( nxx, bxx )


    IMPLICIT NONE

    INTEGER(iwp) :: k         !< loop index
    INTEGER(iwp) :: n_k       !< length scale nXX in height k
    INTEGER(iwp) :: n_k2      !< n_k * 2
    INTEGER(iwp) :: nf        !< index for length scales

    REAL(wp) :: bdenom        !< denominator for filter functions bXX
    REAL(wp) :: qsi = 1.0_wp  !< minimization factor

    INTEGER(iwp), DIMENSION(:) :: nxx(nzb:nzt+1)           !< length scale (in gp)

    REAL(wp), DIMENSION(:,:) :: bxx(-merg:merg,nzb:nzt+1)  !< filter function


    bxx = 0.0_wp

    DO  k = nzb, nzt+1
       bdenom = 0.0_wp
       n_k    = nxx(k)
       IF ( n_k /= 0 )  THEN
          n_k2 = n_k * 2

!
!--       ( Eq.10 )^2
          DO  nf = -n_k2, n_k2
             bdenom = bdenom + EXP( -qsi * pi * ABS(nf) / n_k )**2
          ENDDO

!
!--       ( Eq.9 )
          bdenom = SQRT( bdenom )
          DO  nf = -n_k2, n_k2
             bxx(nf,k) = EXP( -qsi * pi * ABS(nf) / n_k ) / bdenom
          ENDDO
       ENDIF
    ENDDO

END SUBROUTINE stg_filter_func


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &stg_par for synthetic turbulence generator
!------------------------------------------------------------------------------!
 SUBROUTINE stg_parin


    IMPLICIT NONE

    CHARACTER (LEN=80) ::  line   !< dummy string that contains the current line of the parameter file


    NAMELIST /stg_par/   use_syn_turb_gen

    line = ' '

!
!-- Try to find stg package
    REWIND ( 11 )
    line = ' '
    DO WHILE ( INDEX( line, '&stg_par' ) == 0 )
       READ ( 11, '(A)', END=10 )  line
    ENDDO
    BACKSPACE ( 11 )

!
!-- Read namelist
    READ ( 11, stg_par )

!
!-- Set flag that indicates that the synthetic turbulence generator is switched
!-- on
    syn_turb_gen = .TRUE.


 10 CONTINUE

 END SUBROUTINE stg_parin


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine reads the respective restart data.
!------------------------------------------------------------------------------!
 SUBROUTINE stg_rrd_global( found )


    USE control_parameters,                                                    &
        ONLY: length, restart_string


    IMPLICIT NONE

    LOGICAL, INTENT(OUT)  ::  found !< flag indicating if variable was found


    found = .TRUE.


    SELECT CASE ( restart_string(1:length) )

       CASE ( 'mc_factor' )
          READ ( 13 )  mc_factor
       CASE ( 'use_syn_turb_gen' )
          READ ( 13 )  use_syn_turb_gen

       CASE DEFAULT

          found = .FALSE.

    END SELECT


 END SUBROUTINE stg_rrd_global


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the respective restart data.
!------------------------------------------------------------------------------!
 SUBROUTINE stg_wrd_global


    IMPLICIT NONE

    CALL wrd_write_string( 'mc_factor' )
    WRITE ( 14 )  mc_factor

    CALL wrd_write_string( 'use_syn_turb_gen' )
    WRITE ( 14 )  use_syn_turb_gen


 END SUBROUTINE stg_wrd_global


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Create turbulent inflow fields for u, v, w with prescribed length scales and
!> Reynolds stress tensor after a method of Xie and Castro (2008), modified
!> following suggestions of Kim et al. (2013), and using a Lund rotation
!> (Lund, 1998).
!------------------------------------------------------------------------------!
 SUBROUTINE stg_main


    USE arrays_3d,                                                             &
        ONLY:  dzw

    USE control_parameters,                                                    &
        ONLY:  dt_3d, forcing, intermediate_timestep_count,  nest_domain,      &
               rans_mode, simulated_time, volume_flow_initial

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE indices,                                                               &
        ONLY:  wall_flags_0

    USE statistics,                                                            &
        ONLY:  weight_substep

    USE pmc_interface,                                                         &
        ONLY : rans_mode_parent


    IMPLICIT NONE

    INTEGER(iwp) :: i           !< grid index in x-direction
    INTEGER(iwp) :: j           !< loop index in y-direction
    INTEGER(iwp) :: k           !< loop index in z-direction

    REAL(wp) :: dt_stg          !< wheighted subtimestep
    REAL(wp) :: mc_factor_l     !< local mass flux correction factor
    REAL(wp) :: volume_flow     !< mass flux through lateral boundary
    REAL(wp) :: volume_flow_l   !< local mass flux through lateral boundary

    REAL(wp), DIMENSION(nzb:nzt+1,nxlg:nxrg,5) :: dist_xz !< imposed disturbances at north/south boundary
    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,5) :: dist_yz !< imposed disturbances at left/right boundary


    CALL  cpu_log( log_point(57), 'synthetic_turbulence_gen', 'start' )

!
!-- Calculate time step which is needed for filter functions
    dt_stg = dt_3d * weight_substep(intermediate_timestep_count)

!
!-- Initial value of fu, fv, fw
    IF ( simulated_time == 0.0_wp .AND. .NOT. velocity_seed_initialized )  THEN
       CALL stg_generate_seed_yz( nuy, nuz, buy, buz, fu_yz, id_stg_left )
       CALL stg_generate_seed_yz( nvy, nvz, bvy, bvz, fv_yz, id_stg_left )
       CALL stg_generate_seed_yz( nwy, nwz, bwy, bwz, fw_yz, id_stg_left )

       IF ( forcing  .OR.  ( nest_domain .AND.  rans_mode_parent  .AND.        &
                      .NOT.  rans_mode ) )  THEN
!
!--       Generate turbulence at right boundary
          CALL stg_generate_seed_yz( nuy, nuz, buy, buz, fu_yz, id_stg_right )
          CALL stg_generate_seed_yz( nvy, nvz, bvy, bvz, fv_yz, id_stg_right )
          CALL stg_generate_seed_yz( nwy, nwz, bwy, bwz, fw_yz, id_stg_right )
!
!--       Generate turbulence at north boundary
          CALL stg_generate_seed_xz( nux, nuz, bux, buz, fu_xz, id_stg_north )
          CALL stg_generate_seed_xz( nvx, nvz, bvx, bvz, fv_xz, id_stg_north )
          CALL stg_generate_seed_xz( nwx, nwz, bwx, bwz, fw_xz, id_stg_north )
!
!--       Generate turbulence at south boundary
          CALL stg_generate_seed_xz( nux, nuz, bux, buz, fu_xz, id_stg_south )
          CALL stg_generate_seed_xz( nvx, nvz, bvx, bvz, fv_xz, id_stg_south )
          CALL stg_generate_seed_xz( nwx, nwz, bwx, bwz, fw_xz, id_stg_south )
       ENDIF
       velocity_seed_initialized = .TRUE.
    ENDIF
!
!-- New set of fu, fv, fw
    CALL stg_generate_seed_yz( nuy, nuz, buy, buz, fuo_yz, id_stg_left )
    CALL stg_generate_seed_yz( nvy, nvz, bvy, bvz, fvo_yz, id_stg_left )
    CALL stg_generate_seed_yz( nwy, nwz, bwy, bwz, fwo_yz, id_stg_left )

    IF ( forcing  .OR.  ( nest_domain .AND.  rans_mode_parent  .AND.           &
                   .NOT.  rans_mode ) )  THEN
!
!--       Generate turbulence at right boundary
          CALL stg_generate_seed_yz( nuy, nuz, buy, buz, fuo_yz, id_stg_right )
          CALL stg_generate_seed_yz( nvy, nvz, bvy, bvz, fvo_yz, id_stg_right )
          CALL stg_generate_seed_yz( nwy, nwz, bwy, bwz, fwo_yz, id_stg_right )
!
!--       Generate turbulence at north boundary
          CALL stg_generate_seed_xz( nux, nuz, bux, buz, fuo_xz, id_stg_north )
          CALL stg_generate_seed_xz( nvx, nvz, bvx, bvz, fvo_xz, id_stg_north )
          CALL stg_generate_seed_xz( nwx, nwz, bwx, bwz, fwo_xz, id_stg_north )
!
!--       Generate turbulence at south boundary
          CALL stg_generate_seed_xz( nux, nuz, bux, buz, fuo_xz, id_stg_south )
          CALL stg_generate_seed_xz( nvx, nvz, bvx, bvz, fvo_xz, id_stg_south )
          CALL stg_generate_seed_xz( nwx, nwz, bwx, bwz, fwo_xz, id_stg_south )
    ENDIF
!
!-- Turbulence generation at left and or right boundary
    IF ( myidx == id_stg_left  .OR.  myidx == id_stg_right )  THEN

       DO  j = nysg, nyng
          DO  k = nzb, nzt + 1
!
!--          Update fu, fv, fw following Eq. 14 of Xie and Castro (2008)
             IF ( tu(k) == 0.0_wp )  THEN
                fu_yz(k,j) = fuo_yz(k,j)
             ELSE
                fu_yz(k,j) = fu_yz(k,j) * EXP( -pi * dt_stg * 0.5_wp / tu(k) ) +     &
                         fuo_yz(k,j) * SQRT( 1.0_wp - EXP( -pi * dt_stg / tu(k) ) )
             ENDIF

             IF ( tv(k) == 0.0_wp )  THEN
                fv_yz(k,j) = fvo_yz(k,j)
             ELSE
                fv_yz(k,j) = fv_yz(k,j) * EXP( -pi * dt_stg * 0.5_wp / tv(k) ) +     &
                         fvo_yz(k,j) * SQRT( 1.0_wp - EXP( -pi * dt_stg / tv(k) ) )
             ENDIF

             IF ( tw(k) == 0.0_wp )  THEN
                fw_yz(k,j) = fwo_yz(k,j)
             ELSE
                fw_yz(k,j) = fw_yz(k,j) * EXP( -pi * dt_stg * 0.5_wp / tw(k) ) +     &
                         fwo_yz(k,j) * SQRT( 1.0_wp - EXP( -pi * dt_stg / tw(k) ) )
             ENDIF
!
!--          Lund rotation following Eq. 17 in Xie and Castro (2008).
!--          Additional factors are added to improve the variance of v and w
             IF( k == 0 )  THEN
                dist_yz(k,j,1) = 0.0_wp
                dist_yz(k,j,2) = 0.0_wp
                dist_yz(k,j,3) = 0.0_wp
!                 dist_yz(k,j,4) = 0.0_wp
!                 dist_yz(k,j,5) = 0.0_wp
             ELSE
                dist_yz(k,j,1) = a11(k) * fu_yz(k,j)
                !experimental test of 1.2
                dist_yz(k,j,2) = ( SQRT( a22(k) / MAXVAL(a22) )                &
                                         * 1.2_wp )                            &
                                       * (   a21(k) * fu_yz(k,j)               &
                                           + a22(k) * fv_yz(k,j) )
                dist_yz(k,j,3) = ( SQRT(a33(k) / MAXVAL(a33) )                 &
                                         * 1.3_wp )                            &
                                       * (   a31(k) * fu_yz(k,j)               &
                                           + a32(k) * fv_yz(k,j)               &
                                           + a33(k) * fw_yz(k,j) )
                ! Calculation for pt and e not yet implemented
!                 dist_yz(k,j,4) = 0.0_wp
!                 dist_yz(k,j,5) = 0.0_wp
             ENDIF

          ENDDO
       ENDDO

!
!--    Mass flux correction following Kim et al. (2013)
!--    This correction factor insures that the mass flux is preserved at the
!--    inflow boundary
       IF ( .NOT. forcing  .AND.  .NOT. nest_domain )  THEN
          mc_factor_l = 0.0_wp
          mc_factor   = 0.0_wp
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                mc_factor_l = mc_factor_l + dzw(k)  *                          &
                              ( mean_inflow_profiles(k,1) + dist_yz(k,j,1) )
             ENDDO
          ENDDO

#if defined( __parallel )
          CALL MPI_ALLREDUCE( mc_factor_l, mc_factor,  &
                              1, MPI_REAL, MPI_SUM, comm1dy, ierr )
#else
          mc_factor = mc_factor_l
#endif

          mc_factor = volume_flow_initial(1) / mc_factor

!
!--       Add disturbance at the inflow
          DO  j = nysg, nyng
             DO  k = nzb, nzt+1
                 u(k,j,-nbgp+1:0) = ( mean_inflow_profiles(k,1) +              &
                                      dist_yz(k,j,1)             ) * mc_factor
                 v(k,j,-nbgp:-1)  = ( mean_inflow_profiles(k,2) +              &
                                      dist_yz(k,j,2)             ) * mc_factor
                 w(k,j,-nbgp:-1)  =   dist_yz(k,j,3)               * mc_factor
             ENDDO
          ENDDO

       ELSE
!
!--       First, calculate volume flow at yz boundary
          IF ( myidx == id_stg_left  )  i = nxl
          IF ( myidx == id_stg_right )  i = nxr+1

          volume_flow_l = 0.0_wp
          volume_flow   = 0.0_wp
          mc_factor_l   = 0.0_wp
          mc_factor     = 0.0_wp
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                volume_flow_l = volume_flow_l + u(k,j,i) * dzw(k) * dy         &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,i), 1 ) )

                mc_factor_l = mc_factor_l     + ( u(k,j,i) + dist_yz(k,j,1) )  &
                                                         * dzw(k) * dy         &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                              BTEST( wall_flags_0(k,j,i), 1 ) )
             ENDDO
          ENDDO
#if defined( __parallel )
          CALL MPI_ALLREDUCE( volume_flow_l, volume_flow,                      &
                              1, MPI_REAL, MPI_SUM, comm1dy, ierr )
          CALL MPI_ALLREDUCE( mc_factor_l, mc_factor,                          &
                              1, MPI_REAL, MPI_SUM, comm1dy, ierr )
#else
          volume_flow = volume_flow_l
          mc_factor   = mc_factor_l
#endif

          mc_factor = volume_flow / mc_factor

!
!--       Add disturbances
          IF ( myidx == id_stg_left  )  THEN

             DO  j = nysg, nyng
                DO  k = nzb, nzt+1
                   u(k,j,-nbgp+1:0) = ( u(k,j,-nbgp+1:0) + dist_yz(k,j,1) )    &
                                        * mc_factor
                   v(k,j,-nbgp:-1)  = ( v(k,j,-nbgp:-1)  + dist_yz(k,j,2)  )   &
                                        * mc_factor
                   w(k,j,-nbgp:-1)  = ( w(k,j,-nbgp:-1)  + dist_yz(k,j,3)  )   &
                                        * mc_factor
                ENDDO
             ENDDO
          ENDIF
          IF ( myidx == id_stg_right  )  THEN

             DO  j = nysg, nyng
                DO  k = nzb, nzt+1
                   u(k,j,nxr+1:nxr+nbgp) = ( u(k,j,nxr+1:nxr+nbgp) +           &
                                             dist_yz(k,j,1) ) * mc_factor
                   v(k,j,nxr+1:nxr+nbgp) = ( v(k,j,nxr+1:nxr+nbgp) +           &
                                             dist_yz(k,j,2) ) * mc_factor
                   w(k,j,nxr+1:nxr+nbgp) = ( w(k,j,nxr+1:nxr+nbgp) +           &
                                             dist_yz(k,j,3) ) * mc_factor
                ENDDO
             ENDDO
          ENDIF
       ENDIF

    ENDIF
!
!-- Turbulence generation at north and south boundary
    IF ( myidy == id_stg_north  .OR.  myidy == id_stg_south )  THEN

       DO  i = nxlg, nxrg
          DO  k = nzb, nzt + 1
!
!--          Update fu, fv, fw following Eq. 14 of Xie and Castro (2008)
             IF ( tu(k) == 0.0_wp )  THEN
                fu_xz(k,i) = fuo_xz(k,i)
             ELSE
                fu_xz(k,i) = fu_xz(k,i) * EXP( -pi * dt_stg * 0.5_wp / tu(k) ) +     &
                         fuo_xz(k,i) * SQRT( 1.0_wp - EXP( -pi * dt_stg / tu(k) ) )
             ENDIF

             IF ( tv(k) == 0.0_wp )  THEN
                fv_xz(k,i) = fvo_xz(k,i)
             ELSE
                fv_xz(k,i) = fv_xz(k,i) * EXP( -pi * dt_stg * 0.5_wp / tv(k) ) +     &
                         fvo_xz(k,i) * SQRT( 1.0_wp - EXP( -pi * dt_stg / tv(k) ) )
             ENDIF

             IF ( tw(k) == 0.0_wp )  THEN
                fw_xz(k,i) = fwo_xz(k,i)
             ELSE
                fw_xz(k,i) = fw_xz(k,i) * EXP( -pi * dt_stg * 0.5_wp / tw(k) ) +     &
                         fwo_xz(k,i) * SQRT( 1.0_wp - EXP( -pi * dt_stg / tw(k) ) )
             ENDIF
!
!--          Lund rotation following Eq. 17 in Xie and Castro (2008).
!--          Additional factors are added to improve the variance of v and w
             IF( k == 0 )  THEN
                dist_xz(k,i,1) = 0.0_wp
                dist_xz(k,i,2) = 0.0_wp
                dist_xz(k,i,3) = 0.0_wp

             ELSE
                dist_xz(k,i,1) = a11(k) * fu_xz(k,i)
                !experimental test of 1.2
                dist_xz(k,i,2) = ( SQRT( a22(k) / MAXVAL(a22) )                &
                                         * 1.2_wp )                            &
                                       * (   a21(k) * fu_xz(k,i)               &
                                           + a22(k) * fv_xz(k,i) )
                dist_xz(k,i,3) = ( SQRT(a33(k) / MAXVAL(a33) )                 &
                                         * 1.3_wp )                            &
                                       * (   a31(k) * fu_xz(k,i)               &
                                           + a32(k) * fv_xz(k,i)               &
                                           + a33(k) * fw_xz(k,i) )
             ENDIF

          ENDDO
       ENDDO
!
!--    Mass flux correction following Kim et al. (2013)
!--    This correction factor insures that the mass flux is preserved at the
!--    inflow boundary.
!--    First, calculate volume flow at xz boundary
       IF ( myidy == id_stg_south  ) j = nys
       IF ( myidy == id_stg_north )  j = nyn+1

       volume_flow_l = 0.0_wp
       volume_flow   = 0.0_wp
       mc_factor_l   = 0.0_wp
       mc_factor     = 0.0_wp
       DO  i = nxl, nxr
          DO  k = nzb+1, nzt
             volume_flow_l = volume_flow_l + v(k,j,i) * dzw(k) * dx            &
                                  * MERGE( 1.0_wp, 0.0_wp,                     &
                                           BTEST( wall_flags_0(k,j,i), 2 ) )

             mc_factor_l = mc_factor_l     + ( v(k,j,i) + dist_xz(k,i,2) )     &
                                                      * dzw(k) * dx            &
                                  * MERGE( 1.0_wp, 0.0_wp,                     &
                                           BTEST( wall_flags_0(k,j,i), 2 ) )
          ENDDO
       ENDDO
#if defined( __parallel )
       CALL MPI_ALLREDUCE( volume_flow_l, volume_flow,                         &
                           1, MPI_REAL, MPI_SUM, comm1dx, ierr )
       CALL MPI_ALLREDUCE( mc_factor_l, mc_factor,                             &
                           1, MPI_REAL, MPI_SUM, comm1dx, ierr )
#else
       volume_flow = volume_flow_l
       mc_factor   = mc_factor_l
#endif

       mc_factor = volume_flow / mc_factor

!
!--    Add disturbances
       IF ( myidy == id_stg_south  )  THEN

          DO  i = nxlg, nxrg
             DO  k = nzb, nzt+1
                u(k,-nbgp:-1,i) = ( u(k,-nbgp:-1,i) + dist_xz(k,i,1) )         &
                                     * mc_factor
                v(k,-nbgp:0,i)  = ( v(k,-nbgp:0,i)  + dist_xz(k,i,2)  )        &
                                     * mc_factor
                w(k,-nbgp:-1,i) = ( w(k,-nbgp:-1,i) + dist_xz(k,i,3)  )        &
                                     * mc_factor
             ENDDO
          ENDDO
       ENDIF
       IF ( myidy == id_stg_north  )  THEN

          DO  i = nxlg, nxrg
             DO  k = nzb, nzt+1
                u(k,nyn+1:nyn+nbgp,i) = ( u(k,nyn+1:nyn+nbgp,i) +              &
                                          dist_xz(k,i,1) ) * mc_factor
                v(k,nyn+1:nyn+nbgp,i) = ( v(k,nyn+1:nyn+nbgp,i) +              &
                                          dist_xz(k,i,2) ) * mc_factor
                w(k,nyn+1:nyn+nbgp,i) = ( w(k,nyn+1:nyn+nbgp,i) +              &
                                          dist_xz(k,i,3) ) * mc_factor
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    CALL  cpu_log( log_point(57), 'synthetic_turbulence_gen', 'stop' )

 END SUBROUTINE stg_main

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Generate a set of random number rand_it wich is equal on each PE
!> and calculate the velocity seed f_n.
!> f_n is splitted in vertical direction by the number of PEs in x-direction and
!> and each PE calculates a vertical subsection of f_n. At the the end, all
!> parts are collected to form the full array.
!------------------------------------------------------------------------------!
 SUBROUTINE stg_generate_seed_yz( n_y, n_z, b_y, b_z, f_n, id )


    USE indices,                                                               &
        ONLY: ny

    IMPLICIT NONE

    INTEGER(iwp) :: id          !< core ids at respective boundaries
    INTEGER(iwp) :: j           !< loop index in y-direction
    INTEGER(iwp) :: jj          !< loop index in y-direction
    INTEGER(iwp) :: k           !< loop index in z-direction
    INTEGER(iwp) :: kk          !< loop index in z-direction
    INTEGER(iwp) :: send_count  !< send count for MPI_GATHERV

    INTEGER(iwp), DIMENSION(nzb:nzt+1) :: n_y    !< length scale in y-direction
    INTEGER(iwp), DIMENSION(nzb:nzt+1) :: n_z    !< length scale in z-direction
    INTEGER(iwp), DIMENSION(nzb:nzt+1) :: n_y2   !< n_y*2
    INTEGER(iwp), DIMENSION(nzb:nzt+1) :: n_z2   !< n_z*2

    REAL(wp) :: nyz_inv         !< inverse of number of grid points in yz-slice
    REAL(wp) :: rand_av         !< average of random number
    REAL(wp) :: rand_sigma_inv  !< inverse of stdev of random number

    REAL(wp), DIMENSION(-merg:merg,nzb:nzt+1)    :: b_y     !< filter func in y-dir
    REAL(wp), DIMENSION(-merg:merg,nzb:nzt+1)    :: b_z     !< filter func in z-dir
    REAL(wp), DIMENSION(nzb_x_stg:nzt_x_stg+1,nysg:nyng) :: f_n_l   !<  local velocity seed
    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng)     :: f_n     !<  velocity seed
    REAL(wp), DIMENSION(:,:), ALLOCATABLE        :: rand_it !<  random number


!
!-- Generate random numbers using a seed generated in stg_init.
!-- The set of random numbers are modified to have an average of 0 and
!-- unit variance.
    ALLOCATE( rand_it(nzb-mergp:nzt+1+mergp,-mergp:ny+mergp) )

    rand_av        = 0.0_wp
    rand_sigma_inv = 0.0_wp
    nyz_inv        = 1.0_wp / REAL( ( nzt+1 - nzb+1 ) * ( ny+1 ), KIND=wp )

    DO  j = 0, ny
       DO  k = nzb, nzt+1
          CALL RANDOM_NUMBER( rand_it(k,j) )
          rand_av = rand_av + rand_it(k,j)
       ENDDO
    ENDDO

    rand_av = rand_av * nyz_inv

    DO  j = 0, ny
       DO  k = nzb, nzt+1
          rand_it(k,j)   = rand_it(k,j) - rand_av
          rand_sigma_inv = rand_sigma_inv + rand_it(k,j) ** 2
       ENDDO
    ENDDO

    rand_sigma_inv = 1.0_wp / SQRT(rand_sigma_inv * nyz_inv)

    DO  j = 0, ny
       DO  k = nzb, nzt+1
          rand_it(k,j) = rand_it(k,j) * rand_sigma_inv
       ENDDO
    ENDDO

!
!-- Periodic fill of random number in space
    DO  j = 0, ny
       DO  k = 1, mergp
          rand_it(nzb  -k,j) = rand_it(nzt+2-k,j)    ! bottom margin
          rand_it(nzt+1+k,j) = rand_it(nzb+k-1,j)    ! top margin
       ENDDO
    ENDDO
    DO  j = 1, mergp
       DO  k = nzb-mergp, nzt+1+mergp
          rand_it(k,  -j) = rand_it(k,ny-j+1)        ! south margin
          rand_it(k,ny+j) = rand_it(k,   j-1)        ! north margin
       ENDDO
    ENDDO

!
!-- Generate velocity seed following Eq.6 of Xie and Castro (2008)
    n_y2 = n_y * 2
    n_z2 = n_z * 2
    f_n_l  = 0.0_wp

    DO  j = nysg, nyng
       DO  k = nzb_x_stg, nzt_x_stg+1
          DO  jj = -n_y2(k), n_y2(k)
             DO  kk = -n_z2(k), n_z2(k)
                f_n_l(k,j) = f_n_l(k,j)                                        &
                           + b_y(jj,k) * b_z(kk,k) * rand_it(k+kk,j+jj)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE( rand_it )
!
!-- Gather velocity seeds of full subdomain
    send_count = nzt_x_stg - nzb_x_stg + 1
    IF ( nzt_x_stg == nzt )  send_count = send_count + 1

#if defined( __parallel )
    CALL MPI_GATHERV( f_n_l(nzb_x_stg,nysg), send_count, stg_type_yz_small,    &
                      f_n(nzb+1,nysg), recv_count_yz, displs_yz, stg_type_yz,  &
                      id, comm1dx, ierr )
#else
    f_n(nzb+1:nzt+1,nysg:nyng) = f_n_l(nzb_x_stg:nzt_x_stg+1,nysg:nyng)
#endif


 END SUBROUTINE stg_generate_seed_yz




!------------------------------------------------------------------------------!
! Description:
! ------------
!> Generate a set of random number rand_it wich is equal on each PE
!> and calculate the velocity seed f_n.
!> f_n is splitted in vertical direction by the number of PEs in y-direction and
!> and each PE calculates a vertical subsection of f_n. At the the end, all
!> parts are collected to form the full array.
!------------------------------------------------------------------------------!
 SUBROUTINE stg_generate_seed_xz( n_x, n_z, b_x, b_z, f_n, id )


    USE indices,                                                               &
        ONLY: nx


    IMPLICIT NONE

    INTEGER(iwp) :: id          !< core ids at respective boundaries
    INTEGER(iwp) :: i           !< loop index in x-direction
    INTEGER(iwp) :: ii          !< loop index in x-direction
    INTEGER(iwp) :: k           !< loop index in z-direction
    INTEGER(iwp) :: kk          !< loop index in z-direction
    INTEGER(iwp) :: send_count  !< send count for MPI_GATHERV

    INTEGER(iwp), DIMENSION(nzb:nzt+1) :: n_x    !< length scale in x-direction
    INTEGER(iwp), DIMENSION(nzb:nzt+1) :: n_z    !< length scale in z-direction
    INTEGER(iwp), DIMENSION(nzb:nzt+1) :: n_x2   !< n_y*2
    INTEGER(iwp), DIMENSION(nzb:nzt+1) :: n_z2   !< n_z*2

    REAL(wp) :: nxz_inv         !< inverse of number of grid points in xz-slice
    REAL(wp) :: rand_av         !< average of random number
    REAL(wp) :: rand_sigma_inv  !< inverse of stdev of random number

    REAL(wp), DIMENSION(-merg:merg,nzb:nzt+1)    :: b_x     !< filter func in y-dir
    REAL(wp), DIMENSION(-merg:merg,nzb:nzt+1)    :: b_z     !< filter func in z-dir
    REAL(wp), DIMENSION(nzb_y_stg:nzt_y_stg+1,nxlg:nxrg) :: f_n_l   !<  local velocity seed
    REAL(wp), DIMENSION(nzb:nzt+1,nxlg:nxrg)     :: f_n     !<  velocity seed
    REAL(wp), DIMENSION(:,:), ALLOCATABLE        :: rand_it !<  random number


!
!-- Generate random numbers using a seed generated in stg_init.
!-- The set of random numbers are modified to have an average of 0 and
!-- unit variance.
    ALLOCATE( rand_it(nzb-mergp:nzt+1+mergp,-mergp:nx+mergp) )

    rand_av        = 0.0_wp
    rand_sigma_inv = 0.0_wp
    nxz_inv        = 1.0_wp / REAL( ( nzt+1 - nzb+1 ) * ( nx+1 ), KIND=wp )

    DO  i = 0, nx
       DO  k = nzb, nzt+1
          CALL RANDOM_NUMBER( rand_it(k,i) )
          rand_av = rand_av + rand_it(k,i)
       ENDDO
    ENDDO

    rand_av = rand_av * nxz_inv

    DO  i = 0, nx
       DO  k = nzb, nzt+1
          rand_it(k,i)   = rand_it(k,i) - rand_av
          rand_sigma_inv = rand_sigma_inv + rand_it(k,i) ** 2
       ENDDO
    ENDDO

    rand_sigma_inv = 1.0_wp / SQRT(rand_sigma_inv * nxz_inv)

    DO  i = 0, nx
       DO  k = nzb, nzt+1
          rand_it(k,i) = rand_it(k,i) * rand_sigma_inv
       ENDDO
    ENDDO

!
!-- Periodic fill of random number in space
    DO  i = 0, nx
       DO  k = 1, mergp
          rand_it(nzb-k,i)   = rand_it(nzt+2-k,i)    ! bottom margin
          rand_it(nzt+1+k,i) = rand_it(nzb+k-1,i)    ! top margin
       ENDDO
    ENDDO
    DO  i = 1, mergp
       DO  k = nzb-mergp, nzt+1+mergp
          rand_it(k,-i)   = rand_it(k,nx-i+1)        ! left margin
          rand_it(k,nx+i) = rand_it(k,i-1)           ! right margin
       ENDDO
    ENDDO

!
!-- Generate velocity seed following Eq.6 of Xie and Castro (2008)
    n_x2 = n_x * 2
    n_z2 = n_z * 2
    f_n_l  = 0.0_wp

    DO  i = nxlg, nxrg
       DO  k = nzb_y_stg, nzt_y_stg+1
          DO  ii = -n_x2(k), n_x2(k)
             DO  kk = -n_z2(k), n_z2(k)
                f_n_l(k,i) = f_n_l(k,i)                                        &
                           + b_x(ii,k) * b_z(kk,k) * rand_it(k+kk,i+ii)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE( rand_it )

!
!-- Gather velocity seeds of full subdomain
    send_count = nzt_y_stg - nzb_y_stg + 1
    IF ( nzt_y_stg == nzt )  send_count = send_count + 1


#if defined( __parallel )
    CALL MPI_GATHERV( f_n_l(nzb_y_stg,nxlg), send_count, stg_type_xz_small,    &
                      f_n(nzb+1,nxlg), recv_count_xz, displs_xz, stg_type_xz,  &
                      id, comm1dy, ierr )
#else
    f_n(nzb+1:nzt+1,nxlg:nxrg) = f_n_l(nzb_y_stg:nzt_y_stg+1,nxlg:nxrg)
#endif


 END SUBROUTINE stg_generate_seed_xz

 END MODULE synthetic_turbulence_generator_mod
