!> @file turbulence_closure_mod.f90
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
! Copyright 2017-2018 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
!
!
! Former revisions:
! -----------------
! $Id: turbulence_closure_mod.f90 3086 2018-06-25 09:08:04Z gronemeier $
! bugfix: set rans_const_sigma(1) = 1.3
!
! 3083 2018-06-19 14:03:12Z gronemeier
! - set limits of diss at the end of prognostic equations
! - call production_e to calculate production term of diss
! - limit change of diss to -90% to +100%
! - remove factor 0.5 from diffusion_diss_ij
! - rename c_m into c_0, and c_h into c_4
! - add rans_const_c and rans_const_sigma as namelist parameters
! - add calculation of mixing length for profile output in case of rans_tke_e
! - changed format of annotations to comply with doxygen standards
! - calculate and save dissipation rate during rans_tke_l mode
! - set bc at vertical walls for e, diss, km, kh
! - bugfix: set l_wall = 0.0 within buildings
! - set l_wall at bottom and top boundary (rans-mode)
! - bugfix in production term for dissipation rate
! - bugfix in diffusion of dissipation rate
! - disable check for 1D model if rans_tke_e is used
! - bugfixes for initialization (rans-mode):
!    - correction of dissipation-rate formula
!    - calculate km based on l_wall
!    - initialize diss if 1D model is not used
!
! 3045 2018-05-28 07:55:41Z Giersch
! Error message revised
!
! 3014 2018-05-09 08:42:38Z maronga
! Bugfix: nzb_do and nzt_do were not used for 3d data output
!
! 3004 2018-04-27 12:33:25Z Giersch
! Further allocation checks implemented
!
! 2938 2018-03-27 15:52:42Z suehring
! Further todo's
!
! 2936 2018-03-27 14:49:27Z gronemeier
! - defined l_grid only within this module
! - Moved l_wall definition from modules.f90
! - Get level of highest topography, used to limit upward distance calculation
! - Consider cyclic boundary conditions for mixing length calculation
! - Moved copy of wall_flags into subarray to subroutine
! - Implemented l_wall calculation in case of RANS simulation
! - Moved init of l_black to tcm_init_mixing_length
! - Moved init_mixing_length from init_grid.f90 and
!   renamed it to tcm_init_mixing_length
!
! 2764 2018-01-22 09:25:36Z gronemeier
! Bugfix: remove duplicate SAVE statements
!
! 2746 2018-01-15 12:06:04Z suehring
! Move flag plant canopy to modules
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
! Initial revision
!
!
!
!
! Authors:
! --------
! @author Tobias Gronemeier
!
!
! Description:
! ------------
!> This module contains the available turbulence closures for PALM.
!>
!>
!> @todo test initialization for all possibilities
!>       add OpenMP directives whereever possible
!>       remove debug output variables (dummy1, dummy2, dummy3)
!> @todo Check for random disturbances
!> @note <Enter notes on the module>
!> @bug  TKE-e closure still crashes due to too small dt
!------------------------------------------------------------------------------!
 MODULE turbulence_closure_mod


#if defined( __nopointer )
    USE arrays_3d,                                                             &
        ONLY:  diss, diss_p, dzu, e, e_p, kh, km,                              &
               mean_inflow_profiles, prho, pt, tdiss_m, te_m, tend, u, v, vpt, w
#else
    USE arrays_3d,                                                             &
        ONLY:  diss, diss_1, diss_2, diss_3, diss_p, dzu, e, e_1, e_2, e_3,    &
               e_p, kh, km, mean_inflow_profiles, prho, pt, tdiss_m,           &
               te_m, tend, u, v, vpt, w
#endif

    USE control_parameters,                                                    &
        ONLY:  constant_diffusion, dt_3d, e_init, humidity, inflow_l,          &
               initializing_actions, intermediate_timestep_count,              &
               intermediate_timestep_count_max, kappa, km_constant, les_mw,    &
               ocean, plant_canopy, prandtl_number, prho_reference,            &
               pt_reference, rans_mode, rans_tke_e, rans_tke_l, simulated_time,&
               timestep_scheme, turbulence_closure, turbulent_inflow,          &
               use_upstream_for_tke, vpt_reference, ws_scheme_sca,             &
               stokes_force

    USE advec_ws,                                                              &
        ONLY:  advec_s_ws
    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s

    USE indices,                                                               &
        ONLY:  nbgp, nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg,               &
               nzb, nzb_s_inner, nzb_u_inner, nzb_v_inner, nzb_w_inner, nzt,   &
               wall_flags_0

    USE kinds

    USE pegrid

    USE statistics,                                                            &
        ONLY:  hom, hom_sum, statistic_regions

    USE stokes_force_mod,                                                      &
        ONLY:  stokes_force_s, stokes_production_e

    IMPLICIT NONE


    REAL(wp) ::  c_0                !< constant used for diffusion coefficient and dissipation (dependent on mode RANS/LES)
    REAL(wp) ::  c_1                !< model constant for RANS mode
    REAL(wp) ::  c_2                !< model constant for RANS mode
    REAL(wp) ::  c_3                !< model constant for RANS mode
    REAL(wp) ::  c_4                !< model constant for RANS mode
    REAL(wp) ::  l_max              !< maximum length scale for Blackadar mixing length
    REAL(wp) ::  dsig_e = 1.0_wp    !< factor to calculate Ke from Km (1/sigma_e)
    REAL(wp) ::  dsig_diss = 1.0_wp !< factor to calculate K_diss from Km (1/sigma_diss)
    INTEGER(iwp) ::  surf_e         !< end index of surface elements at given i-j position
    INTEGER(iwp) ::  surf_s         !< start index of surface elements at given i-j position

    REAL(wp), DIMENSION(0:4) :: rans_const_c = &       !< model constants for RANS mode (namelist param)
       (/ 0.55_wp, 1.44_wp, 1.92_wp, 0.0_wp, 0.0_wp /) !> default values fit for standard-tke-e closure

    REAL(wp), DIMENSION(2) :: rans_const_sigma = &     !< model constants for RANS mode, sigma values (sigma_e, sigma_diss) (namelist param)
       (/ 1.0_wp, 1.30_wp /)

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  l_black    !< mixing length according to Blackadar
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  l_grid     !< geometric mean of grid sizes dx, dy, dz

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  l_wall !< near-wall mixing length

    !> @todo remove debug variables
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: diss_prod1, diss_adve1, diss_diff1, &
                                               diss_prod2, diss_adve2, diss_diff2, &
                                               diss_prod3, diss_adve3, diss_diff3, &
                                               dummy1, dummy2, dummy3


    PUBLIC c_0, rans_const_c, rans_const_sigma

!
!-- PALM interfaces:
!-- Input parameter checks to be done in check_parameters
    INTERFACE tcm_check_parameters
       MODULE PROCEDURE tcm_check_parameters
    END INTERFACE tcm_check_parameters

!
!-- Data output checks for 2D/3D data to be done in check_parameters
    INTERFACE tcm_check_data_output
       MODULE PROCEDURE tcm_check_data_output
    END INTERFACE tcm_check_data_output

!
!-- Definition of data output quantities
    INTERFACE tcm_define_netcdf_grid
       MODULE PROCEDURE tcm_define_netcdf_grid
    END INTERFACE tcm_define_netcdf_grid

!
!-- Averaging of 3D data for output
    INTERFACE tcm_3d_data_averaging
       MODULE PROCEDURE tcm_3d_data_averaging
    END INTERFACE tcm_3d_data_averaging

!
!-- Data output of 2D quantities
    INTERFACE tcm_data_output_2d
       MODULE PROCEDURE tcm_data_output_2d
    END INTERFACE tcm_data_output_2d

!
!-- Data output of 3D data
    INTERFACE tcm_data_output_3d
       MODULE PROCEDURE tcm_data_output_3d
    END INTERFACE tcm_data_output_3d

!
!-- Initialization actions
    INTERFACE tcm_init
       MODULE PROCEDURE tcm_init
    END INTERFACE tcm_init

!
!-- deallocate arrays
    INTERFACE tcm_deallocate_arrays
      MODULE PROCEDURE tcm_deallocate_arrays
    END INTERFACE tcm_deallocate_arrays

!
!-- Initialization of arrays
    INTERFACE tcm_init_arrays
       MODULE PROCEDURE tcm_init_arrays
    END INTERFACE tcm_init_arrays

!
!-- Initialization of TKE production term
    INTERFACE production_e_init
       MODULE PROCEDURE production_e_init
    END INTERFACE production_e_init

!
!-- Prognostic equations for TKE and TKE dissipation rate
    INTERFACE tcm_prognostic
       MODULE PROCEDURE tcm_prognostic
       MODULE PROCEDURE tcm_prognostic_ij
    END INTERFACE tcm_prognostic

!
!-- Production term for TKE
    INTERFACE production_e
       MODULE PROCEDURE production_e
       MODULE PROCEDURE production_e_ij
    END INTERFACE production_e

!
!-- Diffusion term for TKE
    INTERFACE diffusion_e
       MODULE PROCEDURE diffusion_e
       MODULE PROCEDURE diffusion_e_ij
    END INTERFACE diffusion_e

!
!-- Diffusion term for TKE dissipation rate
    INTERFACE diffusion_diss
       MODULE PROCEDURE diffusion_diss
       MODULE PROCEDURE diffusion_diss_ij
    END INTERFACE diffusion_diss

!
!-- Mixing length for LES case
    INTERFACE mixing_length_les
       MODULE PROCEDURE mixing_length_les
    END INTERFACE mixing_length_les

!
!-- Mixing length for RANS case
    INTERFACE mixing_length_rans
       MODULE PROCEDURE mixing_length_rans
    END INTERFACE mixing_length_rans

!
!-- Calculate diffusivities
    INTERFACE tcm_diffusivities
       MODULE PROCEDURE tcm_diffusivities
    END INTERFACE tcm_diffusivities

!
!-- Swapping of time levels (required for prognostic variables)
    INTERFACE tcm_swap_timelevel
       MODULE PROCEDURE tcm_swap_timelevel
    END INTERFACE tcm_swap_timelevel

    SAVE

    PRIVATE
!
!-- Add INTERFACES that must be available to other modules (alphabetical order)
    PUBLIC production_e_init, tcm_3d_data_averaging, tcm_check_data_output,    &
           tcm_check_parameters, tcm_data_output_2d, tcm_data_output_3d,       &
           tcm_define_netcdf_grid, tcm_diffusivities, tcm_init,                &
           tcm_init_arrays, tcm_prognostic, tcm_swap_timelevel,                &
           tcm_deallocate_arrays


 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for turbulence closure module.
!> @todo remove rans_mode from initialization namelist and rework checks
!>   The way it is implemented at the moment, the user has to set two variables
!>   so that the RANS mode is working. It would be better if only one parameter
!>   has to be set.
!>   2018-06-18, gronemeier
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_check_parameters

    USE control_parameters,                                                    &
        ONLY:  message_string, nest_domain, neutral, turbulent_inflow,         &
               turbulent_outflow

    IMPLICIT NONE

!
!-- Define which turbulence closure is going to be used
    IF ( rans_mode )  THEN

!
!--    Assign values to constants for RANS mode
       dsig_e    = 1.0_wp / rans_const_sigma(1)
       dsig_diss = 1.0_wp / rans_const_sigma(2)

       c_0 = rans_const_c(0)
       c_1 = rans_const_c(1)
       c_2 = rans_const_c(2)
       !c_3 = rans_const_c(3)   !> @todo clarify how to switch between different models
       c_4 = rans_const_c(4)

       SELECT CASE ( TRIM( turbulence_closure ) )

          CASE ( 'TKE-l' )
             rans_tke_l = .TRUE.

          CASE ( 'TKE-e' )
             rans_tke_e = .TRUE.

          CASE DEFAULT
             message_string = 'Unknown turbulence closure: ' //                &
                              TRIM( turbulence_closure )
             CALL message( 'tcm_check_parameters', 'PA0500', 1, 2, 0, 6, 0 )

       END SELECT

       IF ( turbulent_inflow .OR. turbulent_outflow )  THEN
          message_string = 'turbulent inflow/outflow is not yet '//            &
                           'implemented for RANS mode'
          CALL message( 'tcm_check_parameters', 'PA0501', 1, 2, 0, 6, 0 )
       ENDIF

       message_string = 'RANS mode is still in development! ' //               &
                        '&Not all features of PALM are yet compatible '//      &
                        'with RANS mode. &Use at own risk!'
       CALL message( 'tcm_check_parameters', 'PA0502', 0, 1, 0, 6, 0 )

    ELSE

       c_0 = 0.1_wp !according to Lilly (1967) and Deardorff (1980)

       dsig_e = 1.0_wp !assure to use K_m to calculate TKE instead
                       !of K_e which is used in RANS mode

       SELECT CASE ( TRIM( turbulence_closure ) )

          CASE ( 'Moeng_Wyngaard' )
             les_mw = .TRUE.

          CASE DEFAULT
             !> @todo rework this part so that only one call of this error exists
             message_string = 'Unknown turbulence closure: ' //                &
                              TRIM( turbulence_closure )
             CALL message( 'tcm_check_parameters', 'PA0500', 1, 2, 0, 6, 0 )

       END SELECT

    ENDIF

 END SUBROUTINE tcm_check_parameters

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_check_data_output( var, unit, i, ilen, k )

    USE control_parameters,                                                    &
        ONLY:  data_output, message_string

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  unit     !< unit of output variable
    CHARACTER (LEN=*) ::  var      !< name of output variable

    INTEGER(iwp) ::  i      !< index of var in data_output
    INTEGER(iwp) ::  ilen   !< length of var string
    INTEGER(iwp) ::  k      !< flag if var contains one of '_xy', '_xz' or '_yz'

    SELECT CASE ( TRIM( var ) )

       CASE ( 'diss' )
          unit = 'm2/s3'

       CASE ( 'diss1', 'diss2',                         &                      !> @todo remove later
              'diss_prod1', 'diss_adve1', 'diss_diff1', &
              'diss_prod2', 'diss_adve2', 'diss_diff2', &
              'diss_prod3', 'diss_adve3', 'diss_diff3', 'dummy3'  )
          unit = 'debug output'

       CASE ( 'kh', 'km' )
          unit = 'm2/s'

       CASE DEFAULT
          unit = 'illegal'

    END SELECT

 END SUBROUTINE tcm_check_data_output


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Define appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(OUT) ::  grid_x   !< x grid of output variable
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_y   !< y grid of output variable
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_z   !< z grid of output variable
    CHARACTER (LEN=*), INTENT(IN)  ::  var      !< name of output variable

    LOGICAL, INTENT(OUT) ::  found   !< flag if output variable is found

    found  = .TRUE.

!
!-- Check for the grid
    SELECT CASE ( TRIM( var ) )

       CASE ( 'diss', 'diss_xy', 'diss_xz', 'diss_yz' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu'

       CASE ( 'diss1', 'diss2',                         &                       !> @todo remove later
              'diss_prod1', 'diss_adve1', 'diss_diff1', &
              'diss_prod2', 'diss_adve2', 'diss_diff2', &
              'diss_prod3', 'diss_adve3', 'diss_diff3', 'dummy3' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu'

       CASE ( 'kh', 'kh_xy', 'kh_xz', 'kh_yz' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu'

       CASE ( 'km', 'km_xy', 'km_xz', 'km_yz' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu'

       CASE DEFAULT
          found  = .FALSE.
          grid_x = 'none'
          grid_y = 'none'
          grid_z = 'none'

    END SELECT

 END SUBROUTINE tcm_define_netcdf_grid


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Average 3D data.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_3d_data_averaging( mode, variable )


    USE averaging,                                                             &
        ONLY:  diss_av, kh_av, km_av

    USE control_parameters,                                                    &
        ONLY:  average_count_3d

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  mode       !< flag defining mode 'allocate', 'sum' or 'average'
    CHARACTER (LEN=*) ::  variable   !< name of variable

    INTEGER(iwp) ::  i   !< loop index
    INTEGER(iwp) ::  j   !< loop index
    INTEGER(iwp) ::  k   !< loop index

    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'diss' )
             IF ( .NOT. ALLOCATED( diss_av ) )  THEN
                ALLOCATE( diss_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             diss_av = 0.0_wp

          CASE ( 'kh' )
             IF ( .NOT. ALLOCATED( kh_av ) )  THEN
                ALLOCATE( kh_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             kh_av = 0.0_wp

          CASE ( 'km' )
             IF ( .NOT. ALLOCATED( km_av ) )  THEN
                ALLOCATE( km_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             km_av = 0.0_wp

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'diss' )
             IF ( ALLOCATED( diss_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         diss_av(k,j,i) = diss_av(k,j,i) + diss(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'kh' )
             IF ( ALLOCATED( kh_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         kh_av(k,j,i) = kh_av(k,j,i) + kh(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'km' )
             IF ( ALLOCATED( km_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         km_av(k,j,i) = km_av(k,j,i) + km(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'diss' )
             IF ( ALLOCATED( diss_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         diss_av(k,j,i) = diss_av(k,j,i)                       &
                                        / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'kh' )
             IF ( ALLOCATED( kh_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         kh_av(k,j,i) = kh_av(k,j,i)                           &
                                        / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'km' )
             IF ( ALLOCATED( km_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         km_av(k,j,i) = km_av(k,j,i)                           &
                                        / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

       END SELECT

    ENDIF

 END SUBROUTINE tcm_3d_data_averaging


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Define 2D output variables.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_data_output_2d( av, variable, found, grid, mode, local_pf,     &
                                two_d, nzb_do, nzt_do )

    USE averaging,                                                             &
        ONLY:  diss_av, kh_av, km_av

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  grid       !< name of vertical grid
    CHARACTER (LEN=*) ::  mode       !< either 'xy', 'xz' or 'yz'
    CHARACTER (LEN=*) ::  variable   !< name of variable

    INTEGER(iwp) ::  av   !< flag for (non-)average output
    INTEGER(iwp) ::  i    !< loop index
    INTEGER(iwp) ::  j    !< loop index
    INTEGER(iwp) ::  k    !< loop index
    INTEGER(iwp) ::  nzb_do   !< vertical output index (bottom)
    INTEGER(iwp) ::  nzt_do   !< vertical output index (top)

    LOGICAL ::  found   !< flag if output variable is found
    LOGICAL ::  two_d   !< flag parameter that indicates 2D variables (horizontal cross sections)

    REAL(wp) ::  fill_value = -999.0_wp  !< value for the _FillValue attribute

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf !< local
       !< array to which output data is resorted to

    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )


       CASE ( 'diss_xy', 'diss_xz', 'diss_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO k = nzb_do, nzt_do
                      local_pf(i,j,k) = diss(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( diss_av ) ) THEN
                ALLOCATE( diss_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                diss_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO k = nzb_do, nzt_do
                      local_pf(i,j,k) = diss_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'kh_xy', 'kh_xz', 'kh_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO k = nzb_do, nzt_do
                      local_pf(i,j,k) = kh(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( diss_av ) ) THEN
                ALLOCATE( diss_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                diss_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO k = nzb_do, nzt_do
                      local_pf(i,j,k) = kh_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'km_xy', 'km_xz', 'km_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO k = nzb_do, nzt_do
                      local_pf(i,j,k) = km(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( diss_av ) ) THEN
                ALLOCATE( diss_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                diss_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO k = nzb_do, nzt_do
                      local_pf(i,j,k) = km_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          IF ( mode == 'xy' ) grid = 'zu'

       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT

 END SUBROUTINE tcm_data_output_2d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Define 3D output variables.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )


    USE averaging,                                                             &
        ONLY:  diss_av, kh_av, km_av

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  variable   !< name of variable

    INTEGER(iwp) ::  av     !< flag for (non-)average output
    INTEGER(iwp) ::  i      !< loop index
    INTEGER(iwp) ::  j      !< loop index
    INTEGER(iwp) ::  k      !< loop index
    INTEGER(iwp) ::  nzb_do !< lower limit of the data output (usually 0)
    INTEGER(iwp) ::  nzt_do !< vertical upper limit of the data output (usually nz_do3d)

    LOGICAL ::  found   !< flag if output variable is found

    REAL(wp) ::  fill_value = -999.0_wp  !< value for the _FillValue attribute

    REAL(sp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf   !< local
       !< array to which output data is resorted to


    found = .TRUE.


    SELECT CASE ( TRIM( variable ) )


       CASE ( 'diss' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = diss(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( diss_av ) ) THEN
                ALLOCATE( diss_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                diss_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = diss_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'kh' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = kh(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( kh_av ) ) THEN
                ALLOCATE( kh_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                kh_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = kh_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'km' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = km(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( km_av ) ) THEN
                ALLOCATE( km_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                km_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = km_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'dummy3' )                                                        !> @todo remove later
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = dummy3(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'diss1' )                                                         !> @todo remove later
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = dummy1(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'diss2' )                                                         !> @todo remove later
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = dummy2(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'diss_prod1' )                                                    !> @todo remove later
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = diss_prod1(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'diss_adve1' )                                                    !> @todo remove later
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = diss_adve1(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'diss_diff1' )                                                    !> @todo remove later
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = diss_diff1(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'diss_prod2' )                                                    !> @todo remove later
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = diss_prod2(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'diss_adve2' )                                                    !> @todo remove later
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = diss_adve2(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'diss_diff2' )                                                    !> @todo remove later
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = diss_diff2(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'diss_prod3' )                                                    !> @todo remove later
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = diss_prod3(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'diss_adve3' )                                                    !> @todo remove later
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = diss_adve3(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'diss_diff3' )                                                    !> @todo remove later
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = diss_diff3(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE DEFAULT
          found = .FALSE.

    END SELECT

 END SUBROUTINE tcm_data_output_3d


 SUBROUTINE tcm_deallocate_arrays
    IMPLICIT NONE

    deallocate(kh, km, dummy1, dummy2, dummy3, diss_adve1)
    deallocate(diss_adve2, diss_adve3, diss_prod1, diss_prod2)
    deallocate(diss_prod3, diss_diff1, diss_diff2, diss_diff3)
    deallocate(l_grid, l_wall) 
#if defined( __nopointer )
    deallocate(e,e_p,te_m)
#else
    deallocate(e_1,e_2,e_3)
#endif

 END SUBROUTINE tcm_deallocate_arrays

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate arrays and assign pointers.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_init_arrays
    IMPLICIT NONE

    ALLOCATE( kh(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( km(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

    ALLOCATE( dummy1(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )                           !> @todo remove later
    ALLOCATE( dummy2(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( dummy3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( diss_adve1(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( diss_adve2(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( diss_adve3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( diss_prod1(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( diss_prod2(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( diss_prod3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( diss_diff1(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( diss_diff2(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( diss_diff3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    dummy1 = 0.0_wp
    dummy2 = 0.0_wp
    dummy3 = 0.0_wp
    diss_adve1 = 0.0_wp
    diss_adve2 = 0.0_wp
    diss_adve3 = 0.0_wp
    diss_prod1 = 0.0_wp
    diss_prod2 = 0.0_wp
    diss_prod3 = 0.0_wp
    diss_diff1 = 0.0_wp
    diss_diff2 = 0.0_wp
    diss_diff3 = 0.0_wp

#if defined( __nopointer )
    ALLOCATE( e(nzb:nzt+1,nysg:nyng,nxlg:nxrg)    )
    ALLOCATE( e_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  )
    ALLOCATE( te_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

#else
    ALLOCATE( e_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( e_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( e_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
#endif
!
!-- Allocate arrays required for dissipation.
!-- Please note, if it is a nested run, arrays need to be allocated even if
!-- they do not necessarily need to be transferred, which is attributed to
!-- the design of the model coupler which allocates memory for each variable.
#if ! defined( __nopointer )
!
!-- Initial assignment of pointers
    e  => e_1;   e_p  => e_2;   te_m  => e_3

#endif

 END SUBROUTINE tcm_init_arrays


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of turbulence closure module.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_init

    USE control_parameters,                                                    &
        ONLY:  complex_terrain, dissipation_1d, topography

    USE surface_mod,                                                           &
        ONLY:  get_topography_top_index_ji

    IMPLICIT NONE

    INTEGER(iwp) :: i            !< loop index
    INTEGER(iwp) :: j            !< loop index
    INTEGER(iwp) :: k            !< loop index
    INTEGER(iwp) :: nz_s_shift   !< lower shift index for scalars
    INTEGER(iwp) :: nz_s_shift_l !< local lower shift index in case of turbulent inflow

!
!-- Initialize mixing length
    CALL tcm_init_mixing_length
    dummy3 = l_wall                 !> @todo remove later

!
!-- Actions for initial runs
    IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.            &
         TRIM( initializing_actions ) /= 'cyclic_fill' )  THEN

       IF ( INDEX(initializing_actions, 'set_constant_profiles') /= 0 .OR. &
                INDEX( initializing_actions, 'inifor' ) /= 0 )  THEN

          IF ( constant_diffusion )  THEN
             km = km_constant
             kh = km / prandtl_number
             e  = 0.0_wp
          ELSEIF ( e_init > 0.0_wp )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb+1, nzt
                      km(k,j,i) = c_0 * l_wall(k,j,i) * SQRT( e_init )
                   ENDDO
                ENDDO
             ENDDO
             km(nzb,:,:)   = km(nzb+1,:,:)
             km(nzt+1,:,:) = km(nzt,:,:)
             kh = km / prandtl_number
             e  = e_init
          ELSE
             IF ( .NOT. ocean )  THEN
                kh   = 0.01_wp   ! there must exist an initial diffusion, because
                km   = 0.01_wp   ! otherwise no TKE would be produced by the
                                 ! production terms, as long as not yet
                                 ! e = (u*/cm)**2 at k=nzb+1
             ELSE
                kh   = 0.00001_wp
                km   = 0.00001_wp
             ENDIF
             e    = 0.0_wp
          ENDIF

          IF ( rans_tke_e )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb+1, nzt
                      diss(k,j,i) = c_0**4 * e(k,j,i)**2 / km(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
             diss(nzb,:,:) = diss(nzb+1,:,:)
             diss(nzt+1,:,:) = diss(nzt,:,:)
          ENDIF

       ENDIF
!
!--    Store initial profiles for output purposes etc.
       hom(:,1,23,:) = SPREAD( km(:,nys,nxl), 2, statistic_regions+1 )
       hom(:,1,24,:) = SPREAD( kh(:,nys,nxl), 2, statistic_regions+1 )
!
!--    Initialize old and new time levels.
       te_m = 0.0_wp
       e_p = e
       IF ( rans_tke_e )  THEN
          tdiss_m = 0.0_wp
          diss_p = diss
       ENDIF

    ELSEIF ( TRIM( initializing_actions ) == 'read_restart_data'  .OR.         &
             TRIM( initializing_actions ) == 'cyclic_fill' )                   &
    THEN

!
!--    In case of complex terrain and cyclic fill method as initialization,
!--    shift initial data in the vertical direction for each point in the
!--    x-y-plane depending on local surface height
       IF ( complex_terrain  .AND.                                             &
            TRIM( initializing_actions ) == 'cyclic_fill' )  THEN
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                nz_s_shift = get_topography_top_index_ji( j, i, 's' )

                e(nz_s_shift:nzt+1,j,i)  =  e(0:nzt+1-nz_s_shift,j,i)
                km(nz_s_shift:nzt+1,j,i) = km(0:nzt+1-nz_s_shift,j,i)
                kh(nz_s_shift:nzt+1,j,i) = kh(0:nzt+1-nz_s_shift,j,i)
             ENDDO
          ENDDO
          IF ( rans_tke_e )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   nz_s_shift = get_topography_top_index_ji( j, i, 's' )

                   diss(nz_s_shift:nzt+1,j,i) = diss(0:nzt+1-nz_s_shift,j,i)
                ENDDO
             ENDDO
          ENDIF
       ENDIF

!
!--    Initialization of the turbulence recycling method
       IF ( TRIM( initializing_actions ) == 'cyclic_fill'  .AND.               &
            turbulent_inflow )  THEN
          mean_inflow_profiles(:,5) = hom_sum(:,8,0)   ! e
!
!--       In case of complex terrain, determine vertical displacement at inflow
!--       boundary and adjust mean inflow profiles
          IF ( complex_terrain )  THEN
             IF ( nxlg <= 0 .AND. nxrg >= 0 .AND.  &
                  nysg <= 0 .AND. nyng >= 0        )  THEN
                nz_s_shift_l = get_topography_top_index_ji( 0, 0, 's' )
             ELSE
                nz_s_shift_l = 0
             ENDIF
#if defined( __parallel )
             CALL MPI_ALLREDUCE(nz_s_shift_l, nz_s_shift, 1, MPI_INTEGER,      &
                                MPI_MAX, comm2d, ierr)
#else
             nz_s_shift = nz_s_shift_l
#endif
             mean_inflow_profiles(nz_s_shift:nzt+1,5) =  &
                hom_sum(0:nzt+1-nz_s_shift,8,0)  ! e
          ENDIF
!
!--       Use these mean profiles at the inflow (provided that Dirichlet
!--       conditions are used)
          IF ( inflow_l )  THEN
             DO  j = nysg, nyng
                DO  k = nzb, nzt+1
                   e(k,j,nxlg:-1)  = mean_inflow_profiles(k,5)
                ENDDO
             ENDDO
          ENDIF
       ENDIF
!
!--    Inside buildings set TKE back to zero
       IF ( TRIM( initializing_actions ) == 'cyclic_fill' .AND.                &
            topography /= 'flat' )  THEN
!
!--       Inside buildings set TKE back to zero.
!--       Other scalars (km, kh,...) are ignored at present,
!--       maybe revise later.
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                DO  k = nzb, nzt
                   e(k,j,i)     = MERGE( e(k,j,i), 0.0_wp,                     &
                                         BTEST( wall_flags_0(k,j,i), 0 ) )
                ENDDO
             ENDDO
          ENDDO

          IF ( rans_tke_e )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt
                      diss(k,j,i)    = MERGE( diss(k,j,i), 0.0_wp,             &
                                              BTEST( wall_flags_0(k,j,i), 0 ) )
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF
!
!--    Initialize new time levels (only done in order to set boundary values
!--    including ghost points)
       e_p = e
!
!--    Allthough tendency arrays are set in prognostic_equations, they have
!--    to be predefined here because there they are used (but multiplied with 0)
!--    before they are set.
       te_m = 0.0_wp

       IF ( rans_tke_e )  THEN
          diss_p = diss
          tdiss_m = 0.0_wp
       ENDIF

    ENDIF

 END SUBROUTINE tcm_init


! Description:
! -----------------------------------------------------------------------------!
!> Pre-computation of grid-dependent and near-wall mixing length.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_init_mixing_length

    USE arrays_3d,                                                             &
        ONLY:  dzw, ug, vg, zu, zw

    USE control_parameters,                                                    &
        ONLY:  bc_lr_cyc, bc_ns_cyc, f, kappa, message_string,                 &
               wall_adjustment_factor

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg, nzb,  &
               nzt, wall_flags_0

    USE kinds


    IMPLICIT NONE

    INTEGER(iwp) :: dist_dx        !< found distance devided by dx
    INTEGER(iwp) :: i              !< index variable along x
    INTEGER(iwp) :: ii             !< index variable along x
    INTEGER(iwp) :: j              !< index variable along y
    INTEGER(iwp) :: jj             !< index variable along y
    INTEGER(iwp) :: k              !< index variable along z
    INTEGER(iwp) :: k_max_topo = 0 !< index of maximum topography height
    INTEGER(iwp) :: kk             !< index variable along z
    INTEGER(iwp) :: rad_i          !< search radius in grid points along x
    INTEGER(iwp) :: rad_i_l        !< possible search radius to the left
    INTEGER(iwp) :: rad_i_r        !< possible search radius to the right
    INTEGER(iwp) :: rad_j          !< search radius in grid points along y
    INTEGER(iwp) :: rad_j_n        !< possible search radius to north
    INTEGER(iwp) :: rad_j_s        !< possible search radius to south
    INTEGER(iwp) :: rad_k          !< search radius in grid points along z
    INTEGER(iwp) :: rad_k_b        !< search radius in grid points along negative z
    INTEGER(iwp) :: rad_k_t        !< search radius in grid points along positive z

    INTEGER(KIND=1), DIMENSION(:,:), ALLOCATABLE :: vic_yz !< contains a quarter of a single yz-slice of vicinity

    INTEGER(KIND=1), DIMENSION(:,:,:), ALLOCATABLE :: vicinity !< contains topography information of the vicinity of (i/j/k)

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE :: wall_flags_0_global !< wall_flags_0 of whole domain
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE :: wall_flags_dummy    !< dummy array required for MPI_ALLREDUCE command

    REAL(wp) :: radius           !< search radius in meter

    ALLOCATE( l_grid(1:nzt) )
    ALLOCATE( l_wall(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!
!-- Initialize the mixing length in case of an LES-simulation
    IF ( .NOT. rans_mode )  THEN
!
!--    Compute the grid-dependent mixing length.
       DO  k = 1, nzt
          l_grid(k)  = ( dx * dy * dzw(k) )**0.33333333333333_wp
       ENDDO
!
!--    Initialize near-wall mixing length l_wall only in the vertical direction
!--    for the moment, multiplication with wall_adjustment_factor further below
       l_wall(nzb,:,:)   = l_grid(1)
       DO  k = nzb+1, nzt
          l_wall(k,:,:)  = l_grid(k)
       ENDDO
       l_wall(nzt+1,:,:) = l_grid(nzt)

       DO  k = 1, nzt
          IF ( l_grid(k) > 1.5_wp * dx * wall_adjustment_factor .OR.            &
               l_grid(k) > 1.5_wp * dy * wall_adjustment_factor )  THEN
             WRITE( message_string, * ) 'grid anisotropy exceeds ',             &
                                        'threshold given by only local',        &
                                        ' &horizontal reduction of near_wall ', &
                                        'mixing length l_wall',                 &
                                        ' &starting from height level k = ', k, &
                                        '.'
             CALL message( 'init_grid', 'PA0202', 0, 1, 0, 6, 0 )
             EXIT
          ENDIF
       ENDDO
!
!--    In case of topography: limit near-wall mixing length l_wall further:
!--    Go through all points of the subdomain one by one and look for the closest
!--    surface.
!--    Is this correct in the ocean case?
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Check if current gridpoint belongs to the atmosphere
                IF ( BTEST( wall_flags_0(k,j,i), 0 ) )  THEN
!
!--                Check for neighbouring grid-points.
!--                Vertical distance, down
                   IF ( .NOT. BTEST( wall_flags_0(k-1,j,i), 0 ) )              &
                      l_wall(k,j,i) = MIN( l_grid(k), zu(k) - zw(k-1) )
!
!--                Vertical distance, up
                   IF ( .NOT. BTEST( wall_flags_0(k+1,j,i), 0 ) )              &
                      l_wall(k,j,i) = MIN( l_grid(k), zw(k) - zu(k) )
!
!--                y-distance
                   IF ( .NOT. BTEST( wall_flags_0(k,j-1,i), 0 )  .OR.          &
                        .NOT. BTEST( wall_flags_0(k,j+1,i), 0 ) )              &
                      l_wall(k,j,i) = MIN( l_wall(k,j,i), l_grid(k), 0.5_wp * dy )
!
!--                x-distance
                   IF ( .NOT. BTEST( wall_flags_0(k,j,i-1), 0 )  .OR.          &
                        .NOT. BTEST( wall_flags_0(k,j,i+1), 0 ) )              &
                      l_wall(k,j,i) = MIN( l_wall(k,j,i), l_grid(k), 0.5_wp * dx )
!
!--                 yz-distance (vertical edges, down)
                    IF ( .NOT. BTEST( wall_flags_0(k-1,j-1,i), 0 )  .OR.       &
                         .NOT. BTEST( wall_flags_0(k-1,j+1,i), 0 )  )          &
                      l_wall(k,j,i) = MIN( l_wall(k,j,i), l_grid(k),           &
                                           SQRT( 0.25_wp * dy**2 +             &
                                          ( zu(k) - zw(k-1) )**2 ) )
!
!--                  yz-distance (vertical edges, up)
                    IF ( .NOT. BTEST( wall_flags_0(k+1,j-1,i), 0 )  .OR.       &
                         .NOT. BTEST( wall_flags_0(k+1,j+1,i), 0 )  )          &
                      l_wall(k,j,i) = MIN( l_wall(k,j,i), l_grid(k),           &
                                           SQRT( 0.25_wp * dy**2 +             &
                                          ( zw(k) - zu(k) )**2 ) )
!
!--                 xz-distance (vertical edges, down)
                    IF ( .NOT. BTEST( wall_flags_0(k-1,j,i-1), 0 )  .OR.       &
                         .NOT. BTEST( wall_flags_0(k-1,j,i+1), 0 )  )          &
                      l_wall(k,j,i) = MIN( l_wall(k,j,i), l_grid(k),           &
                                           SQRT( 0.25_wp * dx**2 +             &
                                          ( zu(k) - zw(k-1) )**2 ) )
!
!--                 xz-distance (vertical edges, up)
                    IF ( .NOT. BTEST( wall_flags_0(k+1,j,i-1), 0 )  .OR.       &
                         .NOT. BTEST( wall_flags_0(k+1,j,i+1), 0 )  )          &
                     l_wall(k,j,i) = MIN( l_wall(k,j,i), l_grid(k),            &
                                           SQRT( 0.25_wp * dx**2 +             &
                                          ( zw(k) - zu(k) )**2 ) )
!
!--                xy-distance (horizontal edges)
                   IF ( .NOT. BTEST( wall_flags_0(k,j-1,i-1), 0 )  .OR.        &
                        .NOT. BTEST( wall_flags_0(k,j+1,i-1), 0 )  .OR.        &
                        .NOT. BTEST( wall_flags_0(k,j-1,i+1), 0 )  .OR.        &
                        .NOT. BTEST( wall_flags_0(k,j+1,i+1), 0 ) )            &
                      l_wall(k,j,i) = MIN( l_wall(k,j,i), l_grid(k),           &
                                           SQRT( 0.25_wp * ( dx**2 + dy**2 ) ) )
!
!--                xyz distance (vertical and horizontal edges, down)
                   IF ( .NOT. BTEST( wall_flags_0(k-1,j-1,i-1), 0 )  .OR.      &
                        .NOT. BTEST( wall_flags_0(k-1,j+1,i-1), 0 )  .OR.      &
                        .NOT. BTEST( wall_flags_0(k-1,j-1,i+1), 0 )  .OR.      &
                        .NOT. BTEST( wall_flags_0(k-1,j+1,i+1), 0 ) )          &
                      l_wall(k,j,i) = MIN( l_wall(k,j,i), l_grid(k),           &
                                           SQRT( 0.25_wp * ( dx**2 + dy**2 )   &
                                                 +  ( zu(k) - zw(k-1) )**2  ) )
!
!--                xyz distance (vertical and horizontal edges, up)
                   IF ( .NOT. BTEST( wall_flags_0(k+1,j-1,i-1), 0 )  .OR.      &
                        .NOT. BTEST( wall_flags_0(k+1,j+1,i-1), 0 )  .OR.      &
                        .NOT. BTEST( wall_flags_0(k+1,j-1,i+1), 0 )  .OR.      &
                        .NOT. BTEST( wall_flags_0(k+1,j+1,i+1), 0 ) )          &
                      l_wall(k,j,i) = MIN( l_wall(k,j,i), l_grid(k),           &
                                           SQRT( 0.25_wp * ( dx**2 + dy**2 )   &
                                                 +  ( zw(k) - zu(k) )**2  ) )

                ENDIF
             ENDDO
          ENDDO
       ENDDO

    ELSE
!
!-- Initialize the mixing length in case of a RANS simulation
       ALLOCATE( l_black(nzb:nzt+1) )

!
!--    Calculate mixing length according to Blackadar (1962)
       IF ( f /= 0.0_wp )  THEN
          l_max = 2.7E-4_wp * SQRT( ug(nzt+1)**2 + vg(nzt+1)**2 ) /            &
                  ABS( f ) + 1.0E-10_wp
       ELSE
          l_max = 30.0_wp
       ENDIF

       DO  k = nzb, nzt
          l_black(k) = kappa * zu(k) / ( 1.0_wp + kappa * zu(k) / l_max )
       ENDDO

       l_black(nzt+1) = l_black(nzt)

!
!--    Gather topography information of whole domain
       !> @todo reduce amount of data sent by MPI call
       !>   By now, a whole global 3D-array is sent and received with
       !>   MPI_ALLREDUCE although most of the array is 0. This can be
       !>   drastically reduced if only the local subarray is sent and stored
       !>   in a global array. For that, an MPI data type or subarray must be
       !>   defined.
       !>   2018-03-19, gronemeier
       ALLOCATE( wall_flags_0_global(nzb:nzt+1,0:ny,0:nx) )

#if defined ( __parallel )
       ALLOCATE( wall_flags_dummy(nzb:nzt+1,0:ny,0:nx) )
       wall_flags_dummy = 0
       wall_flags_dummy(nzb:nzt+1,nys:nyn,nxl:nxr) =  &
           wall_flags_0(nzb:nzt+1,nys:nyn,nxl:nxr)

       CALL MPI_ALLREDUCE( wall_flags_dummy,                  &
                           wall_flags_0_global,               &
                           (nzt-nzb+2)*(ny+1)*(nx+1),         &
                           MPI_INTEGER, MPI_SUM, comm2d, ierr )
       DEALLOCATE( wall_flags_dummy )
#else
       wall_flags_0_global(nzb:nzt+1,nys:nyn,nxl:nxr) =  &
              wall_flags_0(nzb:nzt+1,nys:nyn,nxl:nxr)
#endif
!
!--    Get height level of highest topography
       DO  i = 0, nx
          DO  j = 0, ny
             DO  k = nzb+1, nzt-1
                IF ( .NOT. BTEST( wall_flags_0_global(k,j,i), 0 ) .AND.  &
                     k > k_max_topo )  &
                   k_max_topo = k
             ENDDO
          ENDDO
       ENDDO

       l_wall(nzb,:,:) = l_black(nzb)
       l_wall(nzt+1,:,:) = l_black(nzt+1)
!
!--    Limit mixing length to either nearest wall or Blackadar mixing length.
!--    For that, analyze each grid point (i/j/k) ("analysed grid point") and
!--    search within its vicinity for the shortest distance to a wall by cal-
!--    culating the distance between the analysed grid point and the "viewed
!--    grid point" if it contains a wall (belongs to topography).
       DO  k = nzb+1, nzt

          radius = l_black(k)  ! radius within walls are searched
!
!--       Set l_wall to its default maximum value (l_back)
          l_wall(k,:,:) = radius

!
!--       Compute search radius as number of grid points in all directions
          rad_i = CEILING( radius / dx )
          rad_j = CEILING( radius / dy )

          DO  kk = 0, nzt-k
             rad_k_t = kk
!
!--          Limit upward search radius to height of maximum topography
             IF ( zu(k+kk)-zu(k) >= radius .OR. k+kk >= k_max_topo )  EXIT
          ENDDO

          DO  kk = 0, k
             rad_k_b = kk
             IF ( zu(k)-zu(k-kk) >= radius )  EXIT
          ENDDO

!
!--       Get maximum vertical radius; necessary for defining arrays
          rad_k = MAX( rad_k_b, rad_k_t )
!
!--       When analysed grid point lies above maximum topography, set search
!--       radius to 0 if the distance between the analysed grid point and max
!--       topography height is larger than the maximum search radius
          IF ( zu(k-rad_k_b) > zu(k_max_topo) )  rad_k_b = 0
!
!--       Search within vicinity only if the vertical search radius is >0
          IF ( rad_k_b /= 0 .OR. rad_k_t /= 0 )  THEN

             !> @note shape of vicinity is larger in z direction
             !>   Shape of vicinity is two grid points larger than actual search
             !>   radius in vertical direction. The first and last grid point is
             !>   always set to 1 to asure correct detection of topography. See
             !>   function "shortest_distance" for details.
             !>   2018-03-16, gronemeier
             ALLOCATE( vicinity(-rad_k-1:rad_k+1,-rad_j:rad_j,-rad_i:rad_i) )
             ALLOCATE( vic_yz(0:rad_k+1,0:rad_j) )

             vicinity = 1

             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Start search only if (i/j/k) belongs to atmosphere
                   IF ( BTEST( wall_flags_0(k,j,i), 0 )  )  THEN
!
!--                   Reset topography within vicinity
                      vicinity(-rad_k:rad_k,:,:) = 0
!
!--                   Copy area surrounding analysed grid point into vicinity.
!--                   First, limit size of data copied to vicinity by the domain
!--                   border
                      rad_i_l = MIN( rad_i, i )
                      rad_i_r = MIN( rad_i, nx-i )

                      rad_j_s = MIN( rad_j, j )
                      rad_j_n = MIN( rad_j, ny-j )

                      CALL copy_into_vicinity( k, j, i,           &
                                               -rad_k_b, rad_k_t, &
                                               -rad_j_s, rad_j_n, &
                                               -rad_i_l, rad_i_r  )
!
!--                   In case of cyclic boundaries, copy parts into vicinity
!--                   where vicinity reaches over the domain borders.
                      IF ( bc_lr_cyc )  THEN
!
!--                      Vicinity reaches over left domain boundary
                         IF ( rad_i > rad_i_l )  THEN
                            CALL copy_into_vicinity( k, j, nx+rad_i_l+1, &
                                                     -rad_k_b, rad_k_t,  &
                                                     -rad_j_s, rad_j_n,  &
                                                     -rad_i, -rad_i_l-1  )
!
!--                         ...and over southern domain boundary
                            IF ( bc_ns_cyc .AND. rad_j > rad_j_s )  &
                               CALL copy_into_vicinity( k, ny+rad_j_s+1,    &
                                                        nx+rad_i_l+1,       &
                                                        -rad_k_b, rad_k_t,  &
                                                        -rad_j, -rad_j_s-1, &
                                                        -rad_i, -rad_i_l-1  )
!
!--                         ...and over northern domain boundary
                            IF ( bc_ns_cyc .AND. rad_j > rad_j_n )  &
                               CALL copy_into_vicinity( k, 0-rad_j_n-1,    &
                                                        nx+rad_i_l+1,      &
                                                        -rad_k_b, rad_k_t, &
                                                         rad_j_n+1, rad_j, &
                                                        -rad_i, -rad_i_l-1 )
                         ENDIF
!
!--                      Vicinity reaches over right domain boundary
                         IF ( rad_i > rad_i_r )  THEN
                            CALL copy_into_vicinity( k, j, 0-rad_i_r-1, &
                                                     -rad_k_b, rad_k_t, &
                                                     -rad_j_s, rad_j_n, &
                                                      rad_i_r+1, rad_i  )
!
!--                         ...and over southern domain boundary
                            IF ( bc_ns_cyc .AND. rad_j > rad_j_s )  &
                               CALL copy_into_vicinity( k, ny+rad_j_s+1,    &
                                                        0-rad_i_r-1,        &
                                                        -rad_k_b, rad_k_t,  &
                                                        -rad_j, -rad_j_s-1, &
                                                         rad_i_r+1, rad_i   )
!
!--                         ...and over northern domain boundary
                            IF ( bc_ns_cyc .AND. rad_j > rad_j_n )  &
                               CALL copy_into_vicinity( k, 0-rad_j_n-1,    &
                                                        0-rad_i_r-1,       &
                                                        -rad_k_b, rad_k_t, &
                                                         rad_j_n+1, rad_j, &
                                                         rad_i_r+1, rad_i  )
                         ENDIF
                      ENDIF

                      IF ( bc_ns_cyc )  THEN
!
!--                      Vicinity reaches over southern domain boundary
                         IF ( rad_j > rad_j_s )  &
                            CALL copy_into_vicinity( k, ny+rad_j_s+1, i, &
                                                     -rad_k_b, rad_k_t,  &
                                                     -rad_j, -rad_j_s-1, &
                                                     -rad_i_l, rad_i_r   )
!
!--                      Vicinity reaches over northern domain boundary
                         IF ( rad_j > rad_j_n )  &
                            CALL copy_into_vicinity( k, 0-rad_j_n-1, i, &
                                                     -rad_k_b, rad_k_t, &
                                                      rad_j_n+1, rad_j, &
                                                      rad_i_l, rad_i_r  )
                      ENDIF
!
!--                   Search for walls only if there is any within vicinity
                      IF ( MAXVAL( vicinity(-rad_k:rad_k,:,:) ) /= 0 )  THEN
!
!--                      Search within first half (positive x)
                         dist_dx = rad_i
                         DO  ii = 0, dist_dx
!
!--                         Search along vertical direction only if below
!--                         maximum topography
                            IF ( rad_k_t > 0 ) THEN
!
!--                            Search for walls within octant (+++)
                               vic_yz = vicinity(0:rad_k+1,0:rad_j,ii)
                               l_wall(k,j,i) = MIN( l_wall(k,j,i),             &
                                       shortest_distance( vic_yz, .TRUE., ii ) )
!
!--                            Search for walls within octant (+-+)
!--                            Switch order of array so that the analysed grid
!--                            point is always located at (0/0) (required by
!--                            shortest_distance").
                               vic_yz = vicinity(0:rad_k+1,0:-rad_j:-1,ii)
                               l_wall(k,j,i) = MIN( l_wall(k,j,i),             &
                                       shortest_distance( vic_yz, .TRUE., ii ) )

                            ENDIF
!
!--                         Search for walls within octant (+--)
                            vic_yz = vicinity(0:-rad_k-1:-1,0:-rad_j:-1,ii)
                            l_wall(k,j,i) = MIN( l_wall(k,j,i),                &
                                      shortest_distance( vic_yz, .FALSE., ii ) )
!
!--                         Search for walls within octant (++-)
                            vic_yz = vicinity(0:-rad_k-1:-1,0:rad_j,ii)
                            l_wall(k,j,i) = MIN( l_wall(k,j,i),                &
                                      shortest_distance( vic_yz, .FALSE., ii ) )
!
!--                         Reduce search along x by already found distance
                            dist_dx = CEILING( l_wall(k,j,i) / dx )

                         ENDDO
!
!-                       Search within second half (negative x)
                         DO  ii = 0, -dist_dx, -1
!
!--                         Search along vertical direction only if below
!--                         maximum topography
                            IF ( rad_k_t > 0 ) THEN
!
!--                            Search for walls within octant (-++)
                               vic_yz = vicinity(0:rad_k+1,0:rad_j,ii)
                               l_wall(k,j,i) = MIN( l_wall(k,j,i),             &
                                      shortest_distance( vic_yz, .TRUE., -ii ) )
!
!--                            Search for walls within octant (--+)
!--                            Switch order of array so that the analysed grid
!--                            point is always located at (0/0) (required by
!--                            shortest_distance").
                               vic_yz = vicinity(0:rad_k+1,0:-rad_j:-1,ii)
                               l_wall(k,j,i) = MIN( l_wall(k,j,i),             &
                                      shortest_distance( vic_yz, .TRUE., -ii ) )

                            ENDIF
!
!--                         Search for walls within octant (---)
                            vic_yz = vicinity(0:-rad_k-1:-1,0:-rad_j:-1,ii)
                            l_wall(k,j,i) = MIN( l_wall(k,j,i),                &
                                     shortest_distance( vic_yz, .FALSE., -ii ) )
!
!--                         Search for walls within octant (-+-)
                            vic_yz = vicinity(0:-rad_k-1:-1,0:rad_j,ii)
                            l_wall(k,j,i) = MIN( l_wall(k,j,i),                &
                                     shortest_distance( vic_yz, .FALSE., -ii ) )
!
!--                         Reduce search along x by already found distance
                            dist_dx = CEILING( l_wall(k,j,i) / dx )

                         ENDDO

                      ENDIF  !Check for any walls within vicinity

                   ELSE  !Check if (i,j,k) belongs to atmosphere

                      l_wall(k,j,i) = l_black(k)

                   ENDIF

                ENDDO  !j loop
             ENDDO  !i loop

             DEALLOCATE( vicinity )
             DEALLOCATE( vic_yz )

          ENDIF  !check vertical size of vicinity

       ENDDO  !k loop

       DEALLOCATE( wall_flags_0_global )

    ENDIF  !LES or RANS mode

!
!-- Set lateral boundary conditions for l_wall
    CALL exchange_horiz( l_wall, nbgp )

    CONTAINS
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the shortest distance between position (i/j/k)=(0/0/0) and
!> (pos_i/jj/kk), where (jj/kk) is the position of the maximum of 'array'
!> closest to the origin (0/0) of 'array'.
!> @todo this part of PALM does not reproduce the same results for optimized
!>   and debug options for the compiler. This should be fixed
!------------------------------------------------------------------------------!
    REAL FUNCTION shortest_distance( array, orientation, pos_i )

       IMPLICIT NONE

       LOGICAL, INTENT(IN) :: orientation    !< flag if array represents an array oriented upwards (true) or downwards (false)

       INTEGER(iwp), INTENT(IN) :: pos_i     !< x position of the yz-plane 'array'

       INTEGER(iwp) :: jj                    !< loop index

       INTEGER(iwp), DIMENSION(0:rad_j) :: loc_k !< location of closest wall along vertical dimension

       INTEGER(KIND=1), DIMENSION(0:rad_k+1,0:rad_j), INTENT(IN) :: array !< array containing a yz-plane at position pos_i

!
!--    Get coordinate of first maximum along vertical dimension
!--    at each y position of array.
!--    Substract 1 because indices count from 1 instead of 0 by MAXLOC
       loc_k = MAXLOC( array, DIM = 1) - 1

!
!--    Set distance to the default maximum value (=search radius)
       shortest_distance = radius
!
!--    Calculate distance between position (0/0/0) and
!--    position (pos_i/jj/loc(jj)) and only save the shortest distance.
       IF ( orientation ) THEN  !if array is oriented upwards
          DO  jj = 0, rad_j
             shortest_distance =                                               &
                MIN( shortest_distance,                                        &
                     SQRT( MAX(REAL(pos_i, KIND=wp)*dx-0.5_wp*dx, 0.0_wp)**2   &
                         + MAX(REAL(jj, KIND=wp)*dy-0.5_wp*dy, 0.0_wp)**2      &
                         + MAX(zw(loc_k(jj)+k-1)-zu(k), 0.0_wp)**2             &
                         )                                                     &
                   )
          ENDDO
       ELSE  !if array is oriented downwards
          !> @note MAX within zw required to circumvent error at domain border
          !>   At the domain border, if non-cyclic boundary is present, the
          !>   index for zw could be -1, which will be errorneous (zw(-1) does
          !>   not exist). The MAX function limits the index to be at least 0.
          DO  jj = 0, rad_j
             shortest_distance =                                               &
                MIN( shortest_distance,                                        &
                     SQRT( MAX(REAL(pos_i, KIND=wp)*dx-0.5_wp*dx, 0.0_wp)**2   &
                         + MAX(REAL(jj, KIND=wp)*dy-0.5_wp*dy, 0.0_wp)**2      &
                         + MAX(zu(k)-zw(MAX(k-loc_k(jj),0_iwp)), 0.0_wp)**2    &
                         )                                                     &
                   )
          ENDDO
       ENDIF

    END FUNCTION

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Copy a subarray of size (kb:kt,js:jn,il:ir) centered around grid point
!> (kp,jp,ip) containing the first bit of wall_flags_0 into the array
!> 'vicinity'. Only copy first bit as this indicates the presence of topography.
!------------------------------------------------------------------------------!
    SUBROUTINE copy_into_vicinity( kp, jp, ip, kb, kt, js, jn, il, ir )

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN) :: il !< left loop boundary
       INTEGER(iwp), INTENT(IN) :: ip !< center position in x-direction
       INTEGER(iwp), INTENT(IN) :: ir !< right loop boundary
       INTEGER(iwp), INTENT(IN) :: jn !< northern loop boundary
       INTEGER(iwp), INTENT(IN) :: jp !< center position in y-direction
       INTEGER(iwp), INTENT(IN) :: js !< southern loop boundary
       INTEGER(iwp), INTENT(IN) :: kb !< bottom loop boundary
       INTEGER(iwp), INTENT(IN) :: kp !< center position in z-direction
       INTEGER(iwp), INTENT(IN) :: kt !< top loop boundary

       INTEGER(iwp) :: i   !< loop index
       INTEGER(iwp) :: j   !< loop index
       INTEGER(iwp) :: k   !< loop index


       DO  i = il, ir
          DO  j = js, jn
             DO  k = kb, kt
                vicinity(k,j,i) = MERGE( 0, 1,               &
                       BTEST( wall_flags_0_global(kp+k,jp+j,ip+i), 0 ) )
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE copy_into_vicinity

 END SUBROUTINE tcm_init_mixing_length


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize virtual velocities used later in production_e.
!------------------------------------------------------------------------------!
 SUBROUTINE production_e_init

    USE arrays_3d,                                                             &
        ONLY:  drho_air_zw, zu

    USE control_parameters,                                                    &
        ONLY:  constant_flux_layer

    USE surface_mod,                                                           &
        ONLY :  surf_def_h, surf_def_v, surf_lsm_h, surf_usm_h

    IMPLICIT NONE

    INTEGER(iwp) ::  i   !< grid index x-direction
    INTEGER(iwp) ::  j   !< grid index y-direction
    INTEGER(iwp) ::  k   !< grid index z-direction
    INTEGER(iwp) ::  m   !< running index surface elements

    IF ( constant_flux_layer )  THEN
!
!--    Calculate a virtual velocity at the surface in a way that the
!--    vertical velocity gradient at k = 1 (u(k+1)-u_0) matches the
!--    Prandtl law (-w'u'/km). This gradient is used in the TKE shear
!--    production term at k=1 (see production_e_ij).
!--    The velocity gradient has to be limited in case of too small km
!--    (otherwise the timestep may be significantly reduced by large
!--    surface winds).
!--    not available in case of non-cyclic boundary conditions.
!--    WARNING: the exact analytical solution would require the determination
!--             of the eddy diffusivity by km = u* * kappa * zp / phi_m.
!--    Default surfaces, upward-facing
       !$OMP PARALLEL DO PRIVATE(i,j,k,m)
       DO  m = 1, surf_def_h(0)%ns

          i = surf_def_h(0)%i(m)
          j = surf_def_h(0)%j(m)
          k = surf_def_h(0)%k(m)
!
!--       Note, calculatione of u_0 and v_0 is not fully accurate, as u/v
!--       and km are not on the same grid. Actually, a further
!--       interpolation of km onto the u/v-grid is necessary. However, the
!--       effect of this error is negligible.
          surf_def_h(0)%u_0(m) = u(k+1,j,i) + surf_def_h(0)%usws(m) *          &
                                     drho_air_zw(k-1) *                        &
                                     ( zu(k+1)    - zu(k-1)    )  /            &
                                     ( km(k,j,i)  + 1.0E-20_wp )
          surf_def_h(0)%v_0(m) = v(k+1,j,i) + surf_def_h(0)%vsws(m) *          &
                                     drho_air_zw(k-1) *                        &
                                     ( zu(k+1)    - zu(k-1)    )  /            &
                                     ( km(k,j,i)  + 1.0E-20_wp )

          IF ( ABS( u(k+1,j,i) - surf_def_h(0)%u_0(m) )  >                     &
               ABS( u(k+1,j,i) - u(k-1,j,i)           )                        &
             )  surf_def_h(0)%u_0(m) = u(k-1,j,i)

          IF ( ABS( v(k+1,j,i) - surf_def_h(0)%v_0(m) )  >                     &
               ABS( v(k+1,j,i) - v(k-1,j,i)           )                        &
             )  surf_def_h(0)%v_0(m) = v(k-1,j,i)

       ENDDO
!
!--    Default surfaces, downward-facing surfaces
       !$OMP PARALLEL DO PRIVATE(i,j,k,m)
       DO  m = 1, surf_def_h(1)%ns

          i = surf_def_h(1)%i(m)
          j = surf_def_h(1)%j(m)
          k = surf_def_h(1)%k(m)

          surf_def_h(1)%u_0(m) = u(k-1,j,i) - surf_def_h(1)%usws(m) *          &
                                     drho_air_zw(k-1) *                        &
                                     ( zu(k+1)    - zu(k-1)    )  /            &
                                     ( km(k,j,i)  + 1.0E-20_wp )
          surf_def_h(1)%v_0(m) = v(k-1,j,i) - surf_def_h(1)%vsws(m) *          &
                                     drho_air_zw(k-1) *                        &
                                     ( zu(k+1)    - zu(k-1)    )  /            &
                                     ( km(k,j,i)  + 1.0E-20_wp )

          IF ( ABS( surf_def_h(1)%u_0(m) - u(k-1,j,i) )  >                     &
               ABS( u(k+1,j,i)           - u(k-1,j,i) )                        &
             )  surf_def_h(1)%u_0(m) = u(k+1,j,i)

          IF ( ABS( surf_def_h(1)%v_0(m) - v(k-1,j,i) )  >                     &
               ABS( v(k+1,j,i)           - v(k-1,j,i) )                        &
             )  surf_def_h(1)%v_0(m) = v(k+1,j,i)

       ENDDO
!
!--    Natural surfaces, upward-facing
       !$OMP PARALLEL DO PRIVATE(i,j,k,m)
       DO  m = 1, surf_lsm_h%ns

          i = surf_lsm_h%i(m)
          j = surf_lsm_h%j(m)
          k = surf_lsm_h%k(m)
!
!--       Note, calculatione of u_0 and v_0 is not fully accurate, as u/v
!--       and km are not on the same grid. Actually, a further
!--       interpolation of km onto the u/v-grid is necessary. However, the
!--       effect of this error is negligible.
          surf_lsm_h%u_0(m) = u(k+1,j,i) + surf_lsm_h%usws(m)      *           &
                                        drho_air_zw(k-1) *                     &
                                        ( zu(k+1)   - zu(k-1)    )  /          &
                                        ( km(k,j,i) + 1.0E-20_wp )
          surf_lsm_h%v_0(m) = v(k+1,j,i) + surf_lsm_h%vsws(m)      *           &
                                        drho_air_zw(k-1) *                     &
                                        ( zu(k+1)   - zu(k-1)    )  /          &
                                        ( km(k,j,i) + 1.0E-20_wp )

          IF ( ABS( u(k+1,j,i) - surf_lsm_h%u_0(m) )  >                        &
               ABS( u(k+1,j,i) - u(k-1,j,i)   )                                &
             )  surf_lsm_h%u_0(m) = u(k-1,j,i)

          IF ( ABS( v(k+1,j,i) - surf_lsm_h%v_0(m) )  >                        &
               ABS( v(k+1,j,i) - v(k-1,j,i)   )                                &
             )  surf_lsm_h%v_0(m) = v(k-1,j,i)

       ENDDO
!
!--    Urban surfaces, upward-facing
       !$OMP PARALLEL DO PRIVATE(i,j,k,m)
       DO  m = 1, surf_usm_h%ns

          i = surf_usm_h%i(m)
          j = surf_usm_h%j(m)
          k = surf_usm_h%k(m)
!
!--       Note, calculatione of u_0 and v_0 is not fully accurate, as u/v
!--       and km are not on the same grid. Actually, a further
!--       interpolation of km onto the u/v-grid is necessary. However, the
!--       effect of this error is negligible.
          surf_usm_h%u_0(m) = u(k+1,j,i) + surf_usm_h%usws(m)      *           &
                                        drho_air_zw(k-1) *                     &
                                        ( zu(k+1)   - zu(k-1)    )  /          &
                                        ( km(k,j,i) + 1.0E-20_wp )
          surf_usm_h%v_0(m) = v(k+1,j,i) + surf_usm_h%vsws(m)      *           &
                                        drho_air_zw(k-1) *                     &
                                        ( zu(k+1)   - zu(k-1)    )  /          &
                                        ( km(k,j,i) + 1.0E-20_wp )

          IF ( ABS( u(k+1,j,i) - surf_usm_h%u_0(m) )  >                        &
               ABS( u(k+1,j,i) - u(k-1,j,i)   )                                &
             )  surf_usm_h%u_0(m) = u(k-1,j,i)

          IF ( ABS( v(k+1,j,i) - surf_usm_h%v_0(m) )  >                        &
               ABS( v(k+1,j,i) - v(k-1,j,i)   )                                &
             )  surf_usm_h%v_0(m) = v(k-1,j,i)

       ENDDO

    ENDIF

 END SUBROUTINE production_e_init


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Prognostic equation for subgrid-scale TKE and TKE dissipation rate.
!> Vector-optimized version
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_prognostic

    USE arrays_3d,                                                             &
        ONLY:  ddzu

    USE control_parameters,                                                    &
        ONLY:  f, scalar_advec, tsc

    USE surface_mod,                                                           &
        ONLY :  surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v, surf_usm_h,    &
                surf_usm_v

    IMPLICIT NONE

    INTEGER(iwp) ::  i       !< loop index
    INTEGER(iwp) ::  j       !< loop index
    INTEGER(iwp) ::  k       !< loop index
    INTEGER(iwp) ::  m       !< loop index
    INTEGER(iwp) ::  surf_e  !< end index of surface elements at given i-j position
    INTEGER(iwp) ::  surf_s  !< start index of surface elements at given i-j position

    REAL(wp)     ::  sbt     !< wheighting factor for sub-time step

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) :: advec  !< advection term of TKE tendency
    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) :: produc !< production term of TKE tendency

!
!-- If required, compute prognostic equation for turbulent kinetic
!-- energy (TKE)
    IF ( .NOT. constant_diffusion )  THEN

       CALL cpu_log( log_point(16), 'tke-equation', 'start' )

       sbt = tsc(2)
       tend = 0.0_wp
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
         CALL advec_s_ws( e, 'e' )
       ENDIF

       ! Compute Stokes-advection if required
       IF ( ocean .AND. stokes_force ) THEN
          CALL stokes_force_s( e )
       ENDIF

       IF ( rans_tke_e )  advec = tend

       CALL production_e

       ! Compute Stokes production if required
       IF ( ocean .AND. stokes_force ) THEN
          CALL stokes_production_e
       ENDIF

!
!--    Save production term for prognostic equation of TKE dissipation rate
       IF ( rans_tke_e )  produc = tend - advec

       CALL diffusion_e( prho, prho_reference )
!
!--    Prognostic equation for TKE.
!--    Eliminate negative TKE values, which can occur due to numerical
!--    reasons in the course of the integration. In such cases the old TKE
!--    value is reduced by 90%.
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                e_p(k,j,i) = e(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +        &
                                                 tsc(3) * te_m(k,j,i) )        &
                                        )                                      &
                                   * MERGE( 1.0_wp, 0.0_wp,                    &
                                             BTEST( wall_flags_0(k,j,i), 0 )   &
                                          )
                IF ( e_p(k,j,i) < 0.0_wp )  e_p(k,j,i) = 0.1_wp * e(k,j,i)
             ENDDO
          ENDDO
       ENDDO

!
!--    Use special boundary condition in case of TKE-e closure
       !> @todo do the same for usm and lsm surfaces
       !>   2018-06-05, gronemeier
       IF ( rans_tke_e )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                surf_s = surf_def_h(0)%start_index(j,i)
                surf_e = surf_def_h(0)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k = surf_def_h(0)%k(m)
                   e_p(k,j,i) = surf_def_h(0)%us(m)**2 / c_0**2
                ENDDO
             ENDDO
          ENDDO
       ENDIF

!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      te_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      te_m(k,j,i) =   -9.5625_wp * tend(k,j,i)                 &
                                     + 5.3125_wp * te_m(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       CALL cpu_log( log_point(16), 'tke-equation', 'stop' )

    ENDIF   ! TKE equation

!
!-- If required, compute prognostic equation for TKE dissipation rate
    IF ( rans_tke_e )  THEN

       CALL cpu_log( log_point(33), 'diss-equation', 'start' )

       sbt = tsc(2)
       tend = 0.0_wp
       CALL advec_s_ws( diss, 'diss' )

!
!--    Production of TKE dissipation rate
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!                tend(k,j,i) = tend(k,j,i) + c_1 * diss(k,j,i) / ( e(k,j,i) + 1.0E-20_wp ) * produc(k)
                tend(k,j,i) = tend(k,j,i) + c_1 * c_0**4 * f / c_4               &  !> @todo needs revision
                      / surf_def_h(0)%us(surf_def_h(0)%start_index(j,i))       &
                      * SQRT(e(k,j,i)) * produc(k,j,i)
             ENDDO
          ENDDO
       ENDDO

       CALL diffusion_diss

!
!
!--    Prognostic equation for TKE dissipation.
!--    Eliminate negative dissipation values, which can occur due to numerical
!--    reasons in the course of the integration. In such cases the old
!--    dissipation value is reduced by 90%.
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                diss_p(k,j,i) = diss(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +  &
                                                 tsc(3) * tdiss_m(k,j,i) )     &
                                        )                                      &
                                   * MERGE( 1.0_wp, 0.0_wp,                    &
                                             BTEST( wall_flags_0(k,j,i), 0 )   &
                                          )
                IF ( diss_p(k,j,i) < 0.0_wp )                                  &
                   diss_p(k,j,i) = 0.1_wp * diss(k,j,i)
             ENDDO
          ENDDO
       ENDDO

!
!--    Use special boundary condition in case of TKE-e closure
       DO  i = nxl, nxr
          DO  j = nys, nyn
             surf_s = surf_def_h(0)%start_index(j,i)
             surf_e = surf_def_h(0)%end_index(j,i)
             DO  m = surf_s, surf_e
                k = surf_def_h(0)%k(m)
                diss_p(k,j,i) = surf_def_h(0)%us(m)**3 / kappa * ddzu(k)
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
                      tdiss_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      tdiss_m(k,j,i) =   -9.5625_wp * tend(k,j,i)              &
                                        + 5.3125_wp * tdiss_m(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       CALL cpu_log( log_point(33), 'diss-equation', 'stop' )

    ENDIF

 END SUBROUTINE tcm_prognostic


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Prognostic equation for subgrid-scale TKE and TKE dissipation rate.
!> Cache-optimized version
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_prognostic_ij( i, j, i_omp, tn )

    USE arrays_3d,                                                             &
        ONLY:  ddzu, diss_l_diss, diss_l_e, diss_s_diss, diss_s_e,             &
               flux_l_diss, flux_l_e, flux_s_diss, flux_s_e,&
               u_p,v_p,w_p

    USE control_parameters,                                                    &
        ONLY:  f, tsc

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE surface_mod,                                                           &
        ONLY :  surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v, surf_usm_h,    &
                surf_usm_v

    use indices, only: nx, ny

    IMPLICIT NONE

    INTEGER(iwp) ::  i       !< loop index x direction
    INTEGER(iwp) ::  i_omp   !< first loop index of i-loop in prognostic_equations
    INTEGER(iwp) ::  j       !< loop index y direction
    INTEGER(iwp) ::  k       !< loop index z direction
    INTEGER(iwp) ::  l       !< loop index
    INTEGER(iwp) ::  m       !< loop index
    INTEGER(iwp) ::  surf_e  !< end index of surface elements at given i-j position
    INTEGER(iwp) ::  surf_s  !< start index of surface elements at given i-j position
    INTEGER(iwp) ::  tn      !< task number of openmp task

    INTEGER(iwp) :: pis = 32 !< debug variable, print from i=pis                !> @todo remove later
    INTEGER(iwp) :: pie = 32 !< debug variable, print until i=pie               !> @todo remove later
    INTEGER(iwp) :: pjs = 26 !< debug variable, print from j=pjs                !> @todo remove later
    INTEGER(iwp) :: pje = 26 !< debug variable, print until j=pje               !> @todo remove later
    INTEGER(iwp) :: pkb = 1  !< debug variable, print from k=pkb                !> @todo remove later
    INTEGER(iwp) :: pkt = 7  !< debug variable, print until k=pkt               !> @todo remove later

    REAL(wp), DIMENSION(nzb:nzt+1) :: dum_adv   !< debug variable               !> @todo remove later
    REAL(wp), DIMENSION(nzb:nzt+1) :: dum_pro   !< debug variable               !> @todo remove later
    REAL(wp), DIMENSION(nzb:nzt+1) :: dum_dif   !< debug variable               !> @todo remove later

5555 FORMAT(A,7(1X,E12.5))   !> @todo remove later

!
!-- If required, compute prognostic equation for turbulent kinetic
!-- energy (TKE)
    IF ( .NOT. constant_diffusion )  THEN

!
!--    Tendency-terms for TKE
       tend(:,j,i) = 0.0_wp
       CALL advec_s_ws( i, j, e, 'e', flux_s_e, diss_s_e, &
                       flux_l_e, diss_l_e , i_omp, tn )

       ! Compute Stokes-advection if required
       IF ( ocean .AND. stokes_force ) THEN
          CALL stokes_force_s( i, j, e )
       ENDIF

       dum_adv = tend(:,j,i)                                                    !> @todo remove later

       CALL production_e( i, j, .FALSE. )

       ! Compute Stokes production if required
       IF ( ocean .AND. stokes_force ) THEN
          CALL stokes_production_e( i, j )
       ENDIF

       dum_pro = tend(:,j,i) - dum_adv                                          !> @todo remove later

       CALL diffusion_e( i, j, prho, prho_reference )

       dum_dif = tend(:,j,i) - dum_adv - dum_pro                                !> @todo remove later

!
!--    Prognostic equation for TKE.
!--    Eliminate negative TKE values, which can occur due to numerical
!--    reasons in the course of the integration. In such cases the old
!--    TKE value is reduced by 90%.
       DO  k = nzb+1, nzt
          e_p(k,j,i) = e(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) +           &
                                              tsc(3) * te_m(k,j,i) )           &
                                  )                                            &
                                 * MERGE( 1.0_wp, 0.0_wp,                      &
                                          BTEST( wall_flags_0(k,j,i), 0 )      &
                                        )
          IF ( e_p(k,j,i) <= 0.0_wp )  e_p(k,j,i) = 0.1_wp * e(k,j,i)
       ENDDO

!
!--    Use special boundary condition in case of TKE-e closure
       IF ( rans_tke_e )  THEN
          surf_s = surf_def_h(0)%start_index(j,i)
          surf_e = surf_def_h(0)%end_index(j,i)
          DO  m = surf_s, surf_e
             k = surf_def_h(0)%k(m)
             e_p(k,j,i) = ( surf_def_h(0)%us(m) / c_0 )**2
          ENDDO

          DO  l = 0, 3
             surf_s = surf_def_v(l)%start_index(j,i)
             surf_e = surf_def_v(l)%end_index(j,i)
             DO  m = surf_s, surf_e
                k = surf_def_v(l)%k(m)
                e_p(k,j,i) = ( surf_def_v(l)%us(m) / c_0 )**2
             ENDDO
          ENDDO
       ENDIF

!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  k = nzb+1, nzt
                te_m(k,j,i) = tend(k,j,i)
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )  THEN
             DO  k = nzb+1, nzt
                te_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +                     &
                                 5.3125_wp * te_m(k,j,i)
             ENDDO
          ENDIF
       ENDIF

!        if ( i >= pis .and. i <= pie .and. j >= pjs .and. j <= pje ) then        !> @todo remove later
!           WRITE(9, *) '------'
!           WRITE(9, '(A,F8.3,1X,F8.3,1X,I2)') 't, dt, int_ts:', simulated_time, dt_3d, intermediate_timestep_count
!           WRITE(9, *) 'i:', i
!           WRITE(9, *) 'j:', j
!           WRITE(9, *) 'k:', pkb, ' - ', pkt
!           WRITE(9, *) '---'
!           WRITE(9, *) 'e:'
!           WRITE(9, 5555) 'adv :', dum_adv(pkb:pkt)
!           WRITE(9, 5555) 'pro :', dum_pro(pkb:pkt)
!           WRITE(9, 5555) 'dif :', dum_dif(pkb:pkt)
!           WRITE(9, 5555) 'tend:', tend(pkb:pkt,j,i)
!           WRITE(9, 5555) 'e_p :', e_p(pkb:pkt,j,i)
!           WRITE(9, 5555) 'e   :', e(pkb:pkt,j,i)
!           FLUSH(9)
!        endif

    ENDIF   ! TKE equation

!
!-- If required, compute prognostic equation for TKE dissipation rate
    IF ( rans_tke_e )  THEN

!
!--    Tendency-terms for dissipation
       tend(:,j,i) = 0.0_wp
       CALL advec_s_ws( i, j, diss, 'diss', flux_s_diss, diss_s_diss,  &
                        flux_l_diss, diss_l_diss, i_omp, tn )

       IF ( intermediate_timestep_count == 1 )  diss_adve1(:,j,i) = tend(:,j,i) !> @todo remove later
       IF ( intermediate_timestep_count == 2 )  diss_adve2(:,j,i) = tend(:,j,i)
       IF ( intermediate_timestep_count == 3 )  diss_adve3(:,j,i) = tend(:,j,i)

!
!--    Production of TKE dissipation rate
       CALL production_e( i, j, .TRUE. )

       IF ( intermediate_timestep_count == 1 )  diss_prod1(:,j,i) = tend(:,j,i) - diss_adve1(:,j,i) !> @todo remove later
       IF ( intermediate_timestep_count == 2 )  diss_prod2(:,j,i) = tend(:,j,i) - diss_adve2(:,j,i)
       IF ( intermediate_timestep_count == 3 )  diss_prod3(:,j,i) = tend(:,j,i) - diss_adve3(:,j,i)

       dum_pro = tend(:,j,i) - dum_adv                                          !> @todo remove later

!
!--    Diffusion term of TKE dissipation rate
       CALL diffusion_diss( i, j )

       IF ( intermediate_timestep_count == 1 )  diss_diff1(:,j,i) = tend(:,j,i) - diss_adve1(:,j,i) - diss_prod1(:,j,i) !> @todo remove later
       IF ( intermediate_timestep_count == 2 )  diss_diff2(:,j,i) = tend(:,j,i) - diss_adve2(:,j,i) - diss_prod2(:,j,i)
       IF ( intermediate_timestep_count == 3 )  diss_diff3(:,j,i) = tend(:,j,i) - diss_adve3(:,j,i) - diss_prod3(:,j,i)
       IF ( intermediate_timestep_count == 3 )  dummy3(:,j,i) = km(:,j,i)

       dum_dif = tend(:,j,i) - dum_adv - dum_pro                                !> @todo remove later

!
!
!--    Prognostic equation for TKE dissipation
!--    Eliminate negative dissipation values, which can occur due to
!--    numerical reasons in the course of the integration. In such cases
!--    the old dissipation value is reduced by 90%.
       DO  k = nzb+1, nzt
          diss_p(k,j,i) = diss(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) +     &
                                                    tsc(3) * tdiss_m(k,j,i) )  &
                                        )                                      &
                                        * MERGE( 1.0_wp, 0.0_wp,               &
                                                BTEST( wall_flags_0(k,j,i), 0 )&
                                               )
       ENDDO

!
!--    Use special boundary condition in case of TKE-e closure
       surf_s = surf_def_h(0)%start_index(j,i)
       surf_e = surf_def_h(0)%end_index(j,i)
       DO  m = surf_s, surf_e
          k = surf_def_h(0)%k(m)
          diss_p(k,j,i) = surf_def_h(0)%us(m)**3 / kappa * ddzu(k)
       ENDDO

       DO  l = 0, 1
          surf_s = surf_def_v(l)%start_index(j,i)
          surf_e = surf_def_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k = surf_def_v(l)%k(m)
             diss_p(k,j,i) = surf_def_v(l)%us(m)**3 / ( kappa * 0.5_wp * dy )
          ENDDO
       ENDDO

       DO  l = 2, 3
          surf_s = surf_def_v(l)%start_index(j,i)
          surf_e = surf_def_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k = surf_def_v(l)%k(m)
             diss_p(k,j,i) = surf_def_v(l)%us(m)**3 / ( kappa * 0.5_wp * dx )
          ENDDO
       ENDDO
!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  k = nzb+1, nzt
                tdiss_m(k,j,i) = tend(k,j,i)
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )  THEN
             DO  k = nzb+1, nzt
                tdiss_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +                  &
                                    5.3125_wp * tdiss_m(k,j,i)
             ENDDO
          ENDIF
       ENDIF
!
!--    Limit change of diss to be between -90% and +100%. Also, set an absolute
!--    minimum value
       DO  k = nzb, nzt+1
          diss_p(k,j,i) = MIN( MAX( diss_p(k,j,i),         &
                                    0.1_wp * diss(k,j,i),  &
                                    0.0001_wp ),           &
                               2.0_wp * diss(k,j,i) )
       ENDDO

       IF ( intermediate_timestep_count == 1 )  dummy1(:,j,i) = diss_p(:,j,i)   !> @todo remove later
       IF ( intermediate_timestep_count == 2 )  dummy2(:,j,i) = diss_p(:,j,i)

!        if ( i >= pis .and. i <= pie .and. j >= pjs .and. j <= pje ) then        !> @todo remove later
!           WRITE(9, *) '---'
!           WRITE(9, *) 'diss:'
!           WRITE(9, 5555) 'adv   :', dum_adv(pkb:pkt)
!           WRITE(9, 5555) 'pro   :', dum_pro(pkb:pkt)
!           WRITE(9, 5555) 'dif   :', dum_dif(pkb:pkt)
!           WRITE(9, 5555) 'tend  :', tend(pkb:pkt,j,i)
!           WRITE(9, 5555) 'diss_p:', diss_p(pkb:pkt,j,i)
!           WRITE(9, 5555) 'diss  :', diss(pkb:pkt,j,i)
!           WRITE(9, *) '---'
!           WRITE(9, 5555) 'km    :', km(pkb:pkt,j,i)
!           flush(9)
!        endif

    ENDIF   ! dissipation equation

 END SUBROUTINE tcm_prognostic_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Production terms (shear + buoyancy) of the TKE.
!> Vector-optimized version
!> @warning The case with constant_flux_layer = F and use_surface_fluxes = T is
!>          not considered well!
!> @todo Adjust production term in case of rans_tke_e simulation
!------------------------------------------------------------------------------!
 SUBROUTINE production_e

    USE arrays_3d,                                                             &
        ONLY:  ddzw, dd2zu, drho_air_zw, q, ql

    USE cloud_parameters,                                                      &
        ONLY:  l_d_cp, l_d_r, pt_d_t, t_d_pt

    USE control_parameters,                                                    &
        ONLY:  cloud_droplets, cloud_physics, constant_flux_layer, g, neutral, &
               rho_reference, use_single_reference_value, use_surface_fluxes,  &
               use_top_fluxes

    USE grid_variables,                                                        &
        ONLY:  ddx, dx, ddy, dy

    USE surface_mod,                                                           &
        ONLY :  surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v, surf_usm_h,    &
                surf_usm_v

    IMPLICIT NONE

    INTEGER(iwp) ::  i       !< running index x-direction
    INTEGER(iwp) ::  j       !< running index y-direction
    INTEGER(iwp) ::  k       !< running index z-direction
    INTEGER(iwp) ::  l       !< running index for different surface type orientation
    INTEGER(iwp) ::  m       !< running index surface elements
    INTEGER(iwp) ::  surf_e  !< end index of surface elements at given i-j position
    INTEGER(iwp) ::  surf_s  !< start index of surface elements at given i-j position

    REAL(wp)     ::  def         !<
    REAL(wp)     ::  flag        !< flag to mask topography
    REAL(wp)     ::  k1          !<
    REAL(wp)     ::  k2          !<
    REAL(wp)     ::  km_neutral  !< diffusion coefficient assuming neutral conditions - used to compute shear production at surfaces
    REAL(wp)     ::  theta       !<
    REAL(wp)     ::  temp        !<
    REAL(wp)     ::  sign_dir    !< sign of wall-tke flux, depending on wall orientation
    REAL(wp)     ::  usvs        !< momentum flux u"v"
    REAL(wp)     ::  vsus        !< momentum flux v"u"
    REAL(wp)     ::  wsus        !< momentum flux w"u"
    REAL(wp)     ::  wsvs        !< momentum flux w"v"

    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  dudx  !< Gradient of u-component in x-direction
    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  dudy  !< Gradient of u-component in y-direction
    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  dudz  !< Gradient of u-component in z-direction
    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  dvdx  !< Gradient of v-component in x-direction
    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  dvdy  !< Gradient of v-component in y-direction
    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  dvdz  !< Gradient of v-component in z-direction
    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  dwdx  !< Gradient of w-component in x-direction
    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  dwdy  !< Gradient of w-component in y-direction
    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  dwdz  !< Gradient of w-component in z-direction

    DO  i = nxl, nxr

       IF ( constant_flux_layer )  THEN

!
!--       Calculate TKE production by shear. Calculate gradients at all grid
!--       points first, gradients at surface-bounded grid points will be
!--       overwritten further below.
          DO  j = nys, nyn
             DO  k = nzb+1, nzt

                dudx(k,j) =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
                dudy(k,j) = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) -            &
                                        u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
                dudz(k,j) = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) -            &
                                        u(k-1,j,i) - u(k-1,j,i+1) ) *          &
                                                         dd2zu(k)

                dvdx(k,j) = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) -            &
                                        v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
                dvdy(k,j) =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
                dvdz(k,j) = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) -            &
                                        v(k-1,j,i) - v(k-1,j+1,i) ) *          &
                                                         dd2zu(k)

                dwdx(k,j) = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) -            &
                                        w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
                dwdy(k,j) = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) -            &
                                        w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
                dwdz(k,j) =           ( w(k,j,i)   - w(k-1,j,i)   ) * ddzw(k)

             ENDDO
          ENDDO

!
!--       Position beneath wall
!--       (2) - Will allways be executed.
!--       'bottom and wall: use u_0,v_0 and wall functions'
          DO  j = nys, nyn
!
!--          Compute gradients at north- and south-facing surfaces.
!--          First, for default surfaces, then for urban surfaces.
!--          Note, so far no natural vertical surfaces implemented
             DO  l = 0, 1
                surf_s = surf_def_v(l)%start_index(j,i)
                surf_e = surf_def_v(l)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k           = surf_def_v(l)%k(m)
                   usvs        = surf_def_v(l)%mom_flux_tke(0,m)
                   wsvs        = surf_def_v(l)%mom_flux_tke(1,m)

                   km_neutral = kappa * ( usvs**2 + wsvs**2 )**0.25_wp         &
                                   * 0.5_wp * dy
!
!--                -1.0 for right-facing wall, 1.0 for left-facing wall
                   sign_dir = MERGE( 1.0_wp, -1.0_wp,                          &
                                     BTEST( wall_flags_0(k,j-1,i), 0 ) )
                   dudy(k,j) = sign_dir * usvs / ( km_neutral + 1E-10_wp )
                   dwdy(k,j) = sign_dir * wsvs / ( km_neutral + 1E-10_wp )
                ENDDO
            ENDDO
!
!--          Compute gradients at east- and west-facing walls
             DO  l = 2, 3
                surf_s = surf_def_v(l)%start_index(j,i)
                surf_e = surf_def_v(l)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k     = surf_def_v(l)%k(m)
                   vsus  = surf_def_v(l)%mom_flux_tke(0,m)
                   wsus  = surf_def_v(l)%mom_flux_tke(1,m)

                   km_neutral = kappa * ( vsus**2 + wsus**2 )**0.25_wp         &
                                      * 0.5_wp * dx
!
!--                -1.0 for right-facing wall, 1.0 for left-facing wall
                   sign_dir = MERGE( 1.0_wp, -1.0_wp,                          &
                                     BTEST( wall_flags_0(k,j,i-1), 0 ) )
                   dvdx(k,j) = sign_dir * vsus / ( km_neutral + 1E-10_wp )
                   dwdx(k,j) = sign_dir * wsus / ( km_neutral + 1E-10_wp )
                ENDDO
            ENDDO
!
!--          Compute gradients at upward-facing surfaces
             surf_s = surf_def_h(0)%start_index(j,i)
             surf_e = surf_def_h(0)%end_index(j,i)
             DO  m = surf_s, surf_e
                k = surf_def_h(0)%k(m)
!
!--             Please note, actually, an interpolation of u_0 and v_0
!--             onto the grid center would be required. However, this
!--             would require several data transfers between 2D-grid and
!--             wall type. The effect of this missing interpolation is
!--             negligible. (See also production_e_init).
                dudz(k,j) = ( u(k+1,j,i) - surf_def_h(0)%u_0(m) ) * dd2zu(k)
                dvdz(k,j) = ( v(k+1,j,i) - surf_def_h(0)%v_0(m) ) * dd2zu(k)

             ENDDO

!--          Compute gradients at downward-facing walls, only for
!--          non-natural default surfaces
             surf_s = surf_def_h(1)%start_index(j,i)
             surf_e = surf_def_h(1)%end_index(j,i)
             DO  m = surf_s, surf_e
                k = surf_def_h(1)%k(m)

                dudz(k,j) = ( surf_def_h(1)%u_0(m) - u(k-1,j,i) ) * dd2zu(k)
                dvdz(k,j) = ( surf_def_h(1)%v_0(m) - v(k-1,j,i) ) * dd2zu(k)

             ENDDO
          ENDDO

          DO  j = nys, nyn
             DO  k = nzb+1, nzt

                def = 2.0_wp * ( dudx(k,j)**2 + dvdy(k,j)**2 + dwdz(k,j)**2 ) + &
                                 dudy(k,j)**2 + dvdx(k,j)**2 + dwdx(k,j)**2 +   &
                                 dwdy(k,j)**2 + dudz(k,j)**2 + dvdz(k,j)**2 +   &
                      2.0_wp * ( dvdx(k,j)*dudy(k,j) + dwdx(k,j)*dudz(k,j)  +   &
                                 dwdy(k,j)*dvdz(k,j) )

                IF ( def < 0.0_wp )  def = 0.0_wp

                flag  = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

                tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def * flag

             ENDDO
          ENDDO

       ELSE

          DO  j = nys, nyn
!
!--          Calculate TKE production by shear. Here, no additional
!--          wall-bounded code is considered.
!--          Why?
             DO  k = nzb+1, nzt

                dudx(k,j)  =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
                dudy(k,j)  = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) -           &
                                         u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
                dudz(k,j)  = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) -           &
                                         u(k-1,j,i) - u(k-1,j,i+1) ) *         &
                                                                     dd2zu(k)

                dvdx(k,j)  = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) -           &
                                         v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
                dvdy(k,j)  =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
                dvdz(k,j)  = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) -           &
                                         v(k-1,j,i) - v(k-1,j+1,i) ) *         &
                                                                     dd2zu(k)

                dwdx(k,j)  = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) -           &
                                         w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
                dwdy(k,j)  = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) -           &
                                         w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
                dwdz(k,j)  =           ( w(k,j,i)   - w(k-1,j,i)   ) *         &
                                                                     ddzw(k)

                def = 2.0_wp * (                                               &
                             dudx(k,j)**2 + dvdy(k,j)**2 + dwdz(k,j)**2        &
                               ) +                                             &
                             dudy(k,j)**2 + dvdx(k,j)**2 + dwdx(k,j)**2 +      &
                             dwdy(k,j)**2 + dudz(k,j)**2 + dvdz(k,j)**2 +      &
                      2.0_wp * (                                               &
                             dvdx(k,j)*dudy(k,j) + dwdx(k,j)*dudz(k,j)  +      &
                             dwdy(k,j)*dvdz(k,j)                               &
                               )

                IF ( def < 0.0_wp )  def = 0.0_wp

                flag  = MERGE( 1.0_wp, 0.0_wp,                                 &
                               BTEST( wall_flags_0(k,j,i), 29 ) )
                tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def * flag

             ENDDO
          ENDDO

       ENDIF

!
!--    If required, calculate TKE production by buoyancy
!--             So far in the ocean no special treatment of density flux
!--             in the bottom and top surface layer
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      tend(k,j,i) = tend(k,j,i) +                              &
                                    kh(k,j,i) * g /                            &
                           MERGE( rho_reference, prho(k,j,i),                  &
                                  use_single_reference_value ) *               &
                                    ( prho(k+1,j,i) - prho(k-1,j,i) ) *        &
                                    dd2zu(k) *                                 &
                                MERGE( 1.0_wp, 0.0_wp,                         &
                                       BTEST( wall_flags_0(k,j,i), 30 )        &
                                     )                            *            &
                                MERGE( 1.0_wp, 0.0_wp,                         &
                                       BTEST( wall_flags_0(k,j,i), 9 )         &
                                     )
                   ENDDO
!
!--                Treatment of near-surface grid points, at up- and down-
!--                ward facing surfaces
                   IF ( use_surface_fluxes )  THEN
                      DO  l = 0, 1
                         surf_s = surf_def_h(l)%start_index(j,i)
                         surf_e = surf_def_h(l)%end_index(j,i)
                         DO  m = surf_s, surf_e
                            k = surf_def_h(l)%k(m)
                            tend(k,j,i) = tend(k,j,i) + g /                    &
                                      MERGE( rho_reference, prho(k,j,i),       &
                                             use_single_reference_value ) *    &
                                      drho_air_zw(k-1) *                       &
                                      surf_def_h(l)%shf(m)
                         ENDDO
                      ENDDO

                   ENDIF

                   IF ( use_top_fluxes )  THEN
                      surf_s = surf_def_h(2)%start_index(j,i)
                      surf_e = surf_def_h(2)%end_index(j,i)
                      DO  m = surf_s, surf_e
                         k = surf_def_h(2)%k(m)
                         tend(k,j,i) = tend(k,j,i) + g /                       &
                                      MERGE( rho_reference, prho(k,j,i),       &
                                             use_single_reference_value ) *    &
                                      drho_air_zw(k) *                         &
                                      surf_def_h(2)%shf(m)
                      ENDDO
                   ENDIF

                ENDDO
    ENDDO

 END SUBROUTINE production_e


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Production terms (shear + buoyancy) of the TKE.
!> Cache-optimized version
!> @warning The case with constant_flux_layer = F and use_surface_fluxes = T is
!>          not considered well!
!> @todo non-neutral case is not yet considered for RANS mode
!------------------------------------------------------------------------------!
 SUBROUTINE production_e_ij( i, j, diss_production )

    USE arrays_3d,                                                             &
        ONLY:  ddzw, dd2zu, drho_air_zw, q, ql

    USE cloud_parameters,                                                      &
        ONLY:  l_d_cp, l_d_r, pt_d_t, t_d_pt

    USE control_parameters,                                                    &
        ONLY:  cloud_droplets, cloud_physics, constant_flux_layer, g, neutral, &
               rho_reference, use_single_reference_value, use_surface_fluxes,  &
               use_top_fluxes

    USE grid_variables,                                                        &
        ONLY:  ddx, dx, ddy, dy

    USE surface_mod,                                                           &
        ONLY :  surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v, surf_usm_h,    &
                surf_usm_v

    IMPLICIT NONE

    LOGICAL :: diss_production

    INTEGER(iwp) ::  i       !< running index x-direction
    INTEGER(iwp) ::  j       !< running index y-direction
    INTEGER(iwp) ::  k       !< running index z-direction
    INTEGER(iwp) ::  l       !< running index for different surface type orientation
    INTEGER(iwp) ::  m       !< running index surface elements
    INTEGER(iwp) ::  surf_e  !< end index of surface elements at given i-j position
    INTEGER(iwp) ::  surf_s  !< start index of surface elements at given i-j position

    REAL(wp)     ::  def         !<
    REAL(wp)     ::  flag        !< flag to mask topography
    REAL(wp)     ::  k1          !<
    REAL(wp)     ::  k2          !<
    REAL(wp)     ::  km_neutral  !< diffusion coefficient assuming neutral conditions - used to compute shear production at surfaces
    REAL(wp)     ::  theta       !<
    REAL(wp)     ::  temp        !<
    REAL(wp)     ::  sign_dir    !< sign of wall-tke flux, depending on wall orientation
    REAL(wp)     ::  usvs        !< momentum flux u"v"
    REAL(wp)     ::  vsus        !< momentum flux v"u"
    REAL(wp)     ::  wsus        !< momentum flux w"u"
    REAL(wp)     ::  wsvs        !< momentum flux w"v"


    REAL(wp), DIMENSION(nzb+1:nzt)  ::  dudx        !< Gradient of u-component in x-direction
    REAL(wp), DIMENSION(nzb+1:nzt)  ::  dudy        !< Gradient of u-component in y-direction
    REAL(wp), DIMENSION(nzb+1:nzt)  ::  dudz        !< Gradient of u-component in z-direction
    REAL(wp), DIMENSION(nzb+1:nzt)  ::  dvdx        !< Gradient of v-component in x-direction
    REAL(wp), DIMENSION(nzb+1:nzt)  ::  dvdy        !< Gradient of v-component in y-direction
    REAL(wp), DIMENSION(nzb+1:nzt)  ::  dvdz        !< Gradient of v-component in z-direction
    REAL(wp), DIMENSION(nzb+1:nzt)  ::  dwdx        !< Gradient of w-component in x-direction
    REAL(wp), DIMENSION(nzb+1:nzt)  ::  dwdy        !< Gradient of w-component in y-direction
    REAL(wp), DIMENSION(nzb+1:nzt)  ::  dwdz        !< Gradient of w-component in z-direction
    REAL(wp), DIMENSION(nzb+1:nzt)  ::  tend_temp   !< temporal tendency

    IF ( constant_flux_layer )  THEN
!
!--    Calculate TKE production by shear. Calculate gradients at all grid
!--    points first, gradients at surface-bounded grid points will be
!--    overwritten further below.
       DO  k = nzb+1, nzt

          dudx(k)  =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
          dudy(k)  = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) -                   &
                                 u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
          dudz(k)  = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) -                   &
                                 u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)

          dvdx(k)  = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) -                   &
                                 v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
          dvdy(k)  =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
          dvdz(k)  = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) -                   &
                                 v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)

          dwdx(k)  = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) -                   &
                                 w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
          dwdy(k)  = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) -                   &
                                 w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
          dwdz(k)  =           ( w(k,j,i)   - w(k-1,j,i)   ) * ddzw(k)

       ENDDO
!
!--    Compute gradients at north- and south-facing surfaces.
!--    Note, no vertical natural surfaces so far.
       DO  l = 0, 1
!
!--       Default surfaces
          surf_s = surf_def_v(l)%start_index(j,i)
          surf_e = surf_def_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_def_v(l)%k(m)
             usvs        = surf_def_v(l)%mom_flux_tke(0,m)
             wsvs        = surf_def_v(l)%mom_flux_tke(1,m)

             km_neutral = kappa * ( usvs**2 + wsvs**2 )**0.25_wp               &
                             * 0.5_wp * dy
!
!--          -1.0 for right-facing wall, 1.0 for left-facing wall
             sign_dir = MERGE( 1.0_wp, -1.0_wp,                                &
                               BTEST( wall_flags_0(k,j-1,i), 0 ) )
             dudy(k) = sign_dir * usvs / ( km_neutral + 1E-10_wp )
             dwdy(k) = sign_dir * wsvs / ( km_neutral + 1E-10_wp )
          ENDDO
      ENDDO
!
!--    Compute gradients at east- and west-facing walls
       DO  l = 2, 3
!
!--       Default surfaces
          surf_s = surf_def_v(l)%start_index(j,i)
          surf_e = surf_def_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_def_v(l)%k(m)
             vsus        = surf_def_v(l)%mom_flux_tke(0,m)
             wsus        = surf_def_v(l)%mom_flux_tke(1,m)

             km_neutral = kappa * ( vsus**2 + wsus**2 )**0.25_wp               &
                                * 0.5_wp * dx
!
!--          -1.0 for right-facing wall, 1.0 for left-facing wall
             sign_dir = MERGE( 1.0_wp, -1.0_wp,                                &
                               BTEST( wall_flags_0(k,j,i-1), 0 ) )
             dvdx(k) = sign_dir * vsus / ( km_neutral + 1E-10_wp )
             dwdx(k) = sign_dir * wsus / ( km_neutral + 1E-10_wp )
          ENDDO
      ENDDO
!
!--    Compute gradients at upward-facing walls, first for
!--    non-natural default surfaces
       surf_s = surf_def_h(0)%start_index(j,i)
       surf_e = surf_def_h(0)%end_index(j,i)
       DO  m = surf_s, surf_e
          k = surf_def_h(0)%k(m)
!
!--       Please note, actually, an interpolation of u_0 and v_0
!--       onto the grid center would be required. However, this
!--       would require several data transfers between 2D-grid and
!--       wall type. The effect of this missing interpolation is
!--       negligible. (See also production_e_init).
          dudz(k)     = ( u(k+1,j,i) - surf_def_h(0)%u_0(m) ) * dd2zu(k)
          dvdz(k)     = ( v(k+1,j,i) - surf_def_h(0)%v_0(m) ) * dd2zu(k)

       ENDDO

!--    Compute gradients at downward-facing walls, only for
!--    non-natural default surfaces
       surf_s = surf_def_h(1)%start_index(j,i)
       surf_e = surf_def_h(1)%end_index(j,i)
       DO  m = surf_s, surf_e
          k = surf_def_h(1)%k(m)

          dudz(k)     = ( surf_def_h(1)%u_0(m) - u(k-1,j,i) ) * dd2zu(k)
          dvdz(k)     = ( surf_def_h(1)%v_0(m) - v(k-1,j,i) ) * dd2zu(k)

       ENDDO

!        IF ( .NOT. rans_tke_e )  THEN

       DO  k = nzb+1, nzt

          def = 2.0_wp * ( dudx(k)**2 + dvdy(k)**2 + dwdz(k)**2 ) +         &
                           dudy(k)**2 + dvdx(k)**2 + dwdx(k)**2   +         &
                           dwdy(k)**2 + dudz(k)**2 + dvdz(k)**2   +         &
                2.0_wp * ( dvdx(k)*dudy(k) +                                &
                           dwdx(k)*dudz(k) +                                &
                           dwdy(k)*dvdz(k) )

          IF ( def < 0.0_wp )  def = 0.0_wp

          flag  = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

          tend_temp(k) = km(k,j,i) * def * flag

       ENDDO

!        ELSE
!
!           DO  k = nzb+1, nzt
! !
! !--          Production term according to Kato and Launder (1993)
!              def = SQRT( ( dudy(k)**2 + dvdz(k)**2 + dwdx(k)**2 +              &
!                            dudz(k)**2 + dvdx(k)**2 + dwdy(k)**2 +              &
!                            2.0_wp * ( dudy(k) * dvdx(k) +                      &
!                                       dvdz(k) * dwdy(k) +                      &
!                                       dwdx(k) * dudz(k) )       )              &
!                        * ( dudy(k)**2 + dvdz(k)**2 + dwdx(k)**2 +              &
!                            dudz(k)**2 + dvdx(k)**2 + dwdy(k)**2 -              &
!                            2.0_wp * ( dudy(k) * dvdx(k) +                      &
!                                       dvdz(k) * dwdy(k) +                      &
!                                       dwdx(k) * dudz(k) )       )              &
!                        )
!
!              IF ( def < 0.0_wp )  def = 0.0_wp
!
!              flag  = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )
!
!              tend_temp(k) = km(k,j,i) * def * flag
!
!           ENDDO
!
!        ENDIF

    ELSE  ! not constant_flux_layer

!        IF ( .NOT. rans_tke_e )  THEN
!
!--       Calculate TKE production by shear. Here, no additional
!--       wall-bounded code is considered.
!--       Why?
          DO  k = nzb+1, nzt

             dudx(k)  =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
             dudy(k)  = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) -                &
                                    u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
             dudz(k)  = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) -                &
                                    u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)

             dvdx(k)  = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) -                &
                                    v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
             dvdy(k)  =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
             dvdz(k)  = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) -                &
                                    v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)

             dwdx(k)  = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) -                &
                                    w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
             dwdy(k)  = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) -                &
                                    w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
             dwdz(k)  =           ( w(k,j,i)   - w(k-1,j,i)   ) * ddzw(k)

             def = 2.0_wp * ( dudx(k)**2 + dvdy(k)**2 + dwdz(k)**2 ) +         &
                              dudy(k)**2 + dvdx(k)**2 + dwdx(k)**2   +         &
                              dwdy(k)**2 + dudz(k)**2 + dvdz(k)**2   +         &
                   2.0_wp * ( dvdx(k)*dudy(k) +                                &
                              dwdx(k)*dudz(k) +                                &
                              dwdy(k)*dvdz(k) )

             IF ( def < 0.0_wp )  def = 0.0_wp

             flag        = MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_0(k,j,i), 29 ) )
             tend_temp(k) = km(k,j,i) * def * flag

          ENDDO

!        ELSE
!
!           DO  k = nzb+1, nzt
!
!              dudx(k)  =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
!              dudy(k)  = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) -                &
!                                     u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
!              dudz(k)  = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) -                &
!                                     u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)
!
!              dvdx(k)  = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) -                &
!                                     v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
!              dvdy(k)  =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
!              dvdz(k)  = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) -                &
!                                     v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)
!
!              dwdx(k)  = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) -                &
!                                     w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
!              dwdy(k)  = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) -                &
!                                     w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
!              dwdz(k)  =           ( w(k,j,i)   - w(k-1,j,i)   ) * ddzw(k)
! !
! !--          Production term according to Kato and Launder (1993)
!              def = SQRT( ( dudy(k)**2 + dvdz(k)**2 + dwdx(k)**2 +              &
!                            dudz(k)**2 + dvdx(k)**2 + dwdy(k)**2 +              &
!                            2.0_wp * ( dudy(k) * dvdx(k) +                      &
!                                       dvdz(k) * dwdy(k) +                      &
!                                       dwdx(k) * dudz(k) )       )              &
!                        * ( dudy(k)**2 + dvdz(k)**2 + dwdx(k)**2 +              &
!                            dudz(k)**2 + dvdx(k)**2 + dwdy(k)**2 -              &
!                            2.0_wp * ( dudy(k) * dvdx(k) +                      &
!                                       dvdz(k) * dwdy(k) +                      &
!                                       dwdx(k) * dudz(k) )       )              &
!                        )
!
!              IF ( def < 0.0_wp )  def = 0.0_wp
!
!              flag        = MERGE( 1.0_wp, 0.0_wp,                              &
!                                   BTEST( wall_flags_0(k,j,i), 29 ) )
!              tend_temp(k) = km(k,j,i) * def * flag
!
!           ENDDO
!
!        ENDIF

    ENDIF

    IF ( .NOT. diss_production )  THEN
!
!--    Production term in case of TKE production
       DO  k = nzb+1, nzt
          tend(k,j,i) = tend(k,j,i) + tend_temp(k)
       ENDDO
    ELSE
!
!--    Production term in case of dissipation-rate production (rans_tke_e)
       DO  k = nzb+1, nzt

          ! Standard TKE-e closure
          tend(k,j,i) = tend(k,j,i) + tend_temp(k) * diss(k,j,i)               &
                                    /( e(k,j,i) + 1.0E-20_wp )                 &
                                    * c_1
!           ! Production according to Koblitz (2013)
!           tend(k,j,i) = tend(k,j,i) + tend_temp(k) * diss(k,j,i)               &
!                                     /( e(k,j,i) + 1.0E-20_wp )                 &
!                                     * ( c_1 + ( c_2 - c_1 )                    &
!                                             * l_wall(k,j,i) / l_max )
!           ! Production according to Detering and Etling (1985)
!           !> @todo us is not correct if there are vertical walls
!           tend(k,j,i) = tend(k,j,i) + tend_temp(k) * SQRT(e(k,j,i))            &
!                                     * c_1 * c_0**3 / c_4 * f                   &
!              / surf_def_h(0)%us(surf_def_h(0)%start_index(j,i))
       ENDDO
    ENDIF

!
!-- If required, calculate TKE production by buoyancy
!--          So far in the ocean no special treatment of density flux in
!--          the bottom and top surface layer
             DO  k = nzb+1, nzt

                tend(k,j,i) = tend(k,j,i) +                                    &
                              kh(k,j,i) * g /                                  &
                              MERGE( rho_reference, prho(k,j,i),               &
                                     use_single_reference_value ) *            &
                              ( prho(k+1,j,i) - prho(k-1,j,i) ) *              &
                              dd2zu(k) *                                       &
                                MERGE( 1.0_wp, 0.0_wp,                         &
                                       BTEST( wall_flags_0(k,j,i), 30 )        &
                                     ) *                                       &
                                MERGE( 1.0_wp, 0.0_wp,                         &
                                       BTEST( wall_flags_0(k,j,i), 9 )         &
                                     )
             ENDDO

             IF ( use_surface_fluxes )  THEN
!
!--             Default surfaces, up- and downward-facing
                DO  l = 0, 1
                   surf_s = surf_def_h(l)%start_index(j,i)
                   surf_e = surf_def_h(l)%end_index(j,i)
                   DO  m = surf_s, surf_e
                      k = surf_def_h(l)%k(m)
                      tend(k,j,i) = tend(k,j,i) + g /                          &
                                MERGE( rho_reference, prho(k,j,i),             &
                                       use_single_reference_value ) *          &
                                drho_air_zw(k-1) *                             &
                                surf_def_h(l)%shf(m)
                   ENDDO
                ENDDO

             ENDIF

             IF ( use_top_fluxes )  THEN
                surf_s = surf_def_h(2)%start_index(j,i)
                surf_e = surf_def_h(2)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k = surf_def_h(2)%k(m)
                   tend(k,j,i) = tend(k,j,i) + g /                             &
                                MERGE( rho_reference, prho(k,j,i),             &
                                       use_single_reference_value ) *          &
                                drho_air_zw(k) *                               &
                                surf_def_h(2)%shf(m)
                ENDDO
             ENDIF

  END SUBROUTINE production_e_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Diffusion and dissipation terms for the TKE.
!> Vector-optimized version
!------------------------------------------------------------------------------!
 SUBROUTINE diffusion_e( var, var_reference )

    USE arrays_3d,                                                             &
        ONLY:  ddzu, ddzw, drho_air, rho_air_zw

    USE grid_variables,                                                        &
        ONLY:  ddx2, ddy2
    USE surface_mod,                                                           &
       ONLY :  bc_h

    IMPLICIT NONE

    INTEGER(iwp) ::  i              !< running index x direction
    INTEGER(iwp) ::  j              !< running index y direction
    INTEGER(iwp) ::  k              !< running index z direction
    INTEGER(iwp) ::  m              !< running index surface elements

    REAL(wp)     ::  flag           !< flag to mask topography
    REAL(wp)     ::  l              !< mixing length
    REAL(wp)     ::  ll             !< adjusted l
    REAL(wp)     ::  var_reference  !< reference temperature

#if defined( __nopointer )
    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  var  !< temperature
#else
    REAL(wp), DIMENSION(:,:,:), POINTER ::  var  !< temperature
#endif
    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  dissipation  !< TKE dissipation


!
!-- Calculate the tendency terms
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
!
!--          Predetermine flag to mask topography
             flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

!
!--          Calculate dissipation
             IF ( les_mw )  THEN

                CALL mixing_length_les( i, j, k, l, ll, var, var_reference )

                dissipation(k,j) = ( 0.19_wp + 0.74_wp * l / ll )              &
                                   * e(k,j,i) * SQRT( e(k,j,i) ) / l

             ELSEIF ( rans_tke_l )  THEN

                CALL mixing_length_rans( i, j, k, l, ll, var, var_reference )

                dissipation(k,j) = c_0**3 * e(k,j,i) * SQRT( e(k,j,i) ) / ll

                diss(k,j,i) = dissipation(k,j) * flag

             ELSEIF ( rans_tke_e )  THEN

                dissipation(k,j) = diss(k,j,i)

             ENDIF

             tend(k,j,i) = tend(k,j,i) + (                                     &
                                           (                                   &
                       ( km(k,j,i)+km(k,j,i+1) ) * ( e(k,j,i+1)-e(k,j,i) )     &
                     - ( km(k,j,i)+km(k,j,i-1) ) * ( e(k,j,i)-e(k,j,i-1) )     &
                                           ) * ddx2  * flag                    &
                                         + (                                   &
                       ( km(k,j,i)+km(k,j+1,i) ) * ( e(k,j+1,i)-e(k,j,i) )     &
                     - ( km(k,j,i)+km(k,j-1,i) ) * ( e(k,j,i)-e(k,j-1,i) )     &
                                           ) * ddy2  * flag                    &
                                         + (                                   &
            ( km(k,j,i)+km(k+1,j,i) ) * ( e(k+1,j,i)-e(k,j,i) ) * ddzu(k+1)    &
                                                          * rho_air_zw(k)      &
          - ( km(k,j,i)+km(k-1,j,i) ) * ( e(k,j,i)-e(k-1,j,i) ) * ddzu(k)      &
                                                          * rho_air_zw(k-1)    &
                                           ) * ddzw(k) * drho_air(k)           &
                                         ) * flag * dsig_e                     &
                          - dissipation(k,j) * flag

          ENDDO
       ENDDO

    ENDDO

 END SUBROUTINE diffusion_e


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Diffusion and dissipation terms for the TKE.
!> Cache-optimized version
!------------------------------------------------------------------------------!
 SUBROUTINE diffusion_e_ij( i, j, var, var_reference )

    USE arrays_3d,                                                             &
        ONLY:  ddzu, ddzw, drho_air, rho_air_zw

    USE grid_variables,                                                        &
        ONLY:  ddx2, ddy2

    USE surface_mod,                                                           &
       ONLY :  bc_h

    IMPLICIT NONE

    INTEGER(iwp) ::  i              !< running index x direction
    INTEGER(iwp) ::  j              !< running index y direction
    INTEGER(iwp) ::  k              !< running index z direction
    INTEGER(iwp) ::  m              !< running index surface elements
    INTEGER(iwp) ::  surf_e         !< End index of surface elements at (j,i)-gridpoint
    INTEGER(iwp) ::  surf_s         !< Start index of surface elements at (j,i)-gridpoint

    REAL(wp)     ::  flag           !< flag to mask topography
    REAL(wp)     ::  l              !< mixing length
    REAL(wp)     ::  ll             !< adjusted l
    REAL(wp)     ::  var_reference  !< reference temperature

#if defined( __nopointer )
    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  var  !< temperature
#else
    REAL(wp), DIMENSION(:,:,:), POINTER ::  var     !< temperature
#endif
    REAL(wp), DIMENSION(nzb+1:nzt) ::  dissipation  !< dissipation of TKE

!
!-- Calculate the mixing length (for dissipation)
    DO  k = nzb+1, nzt
!
!--    Predetermine flag to mask topography
       flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

!
!--    Calculate dissipation...
!--    ...in case of LES
       IF ( les_mw )  THEN

          CALL mixing_length_les( i, j, k, l, ll, var, var_reference )

          dissipation(k) = ( 0.19_wp + 0.74_wp * l / ll )                      &
                           * e(k,j,i) * SQRT( e(k,j,i) ) / l

!
!--    ...in case of RANS
       ELSEIF ( rans_tke_l )  THEN

          CALL mixing_length_rans( i, j, k, l, ll, var, var_reference  )

          dissipation(k) = c_0**3 * e(k,j,i) * SQRT( e(k,j,i) ) / ll

          diss(k,j,i) = dissipation(k) * flag

       ELSEIF ( rans_tke_e )  THEN

          dissipation(k) = diss(k,j,i)

       ENDIF

!
!--    Calculate the tendency term
       tend(k,j,i) = tend(k,j,i) + (                                           &
                                      (                                        &
                      ( km(k,j,i)+km(k,j,i+1) ) * ( e(k,j,i+1)-e(k,j,i) )      &
                    - ( km(k,j,i)+km(k,j,i-1) ) * ( e(k,j,i)-e(k,j,i-1) )      &
                                      ) * ddx2                                 &
                                    + (                                        &
                      ( km(k,j,i)+km(k,j+1,i) ) * ( e(k,j+1,i)-e(k,j,i) )      &
                    - ( km(k,j,i)+km(k,j-1,i) ) * ( e(k,j,i)-e(k,j-1,i) )      &
                                      ) * ddy2                                 &
                                    + (                                        &
           ( km(k,j,i)+km(k+1,j,i) ) * ( e(k+1,j,i)-e(k,j,i) ) * ddzu(k+1)     &
                                                         * rho_air_zw(k)       &
         - ( km(k,j,i)+km(k-1,j,i) ) * ( e(k,j,i)-e(k-1,j,i) ) * ddzu(k)       &
                                                         * rho_air_zw(k-1)     &
                                      ) * ddzw(k) * drho_air(k)                &
                                   ) * flag * dsig_e                           &
                                 - dissipation(k) * flag

    ENDDO

 END SUBROUTINE diffusion_e_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Diffusion term for the TKE dissipation rate
!> Vector-optimized version
!------------------------------------------------------------------------------!
 SUBROUTINE diffusion_diss()
    USE arrays_3d,                                                             &
        ONLY:  ddzu, ddzw, drho_air, rho_air_zw

    USE grid_variables,                                                        &
        ONLY:  ddx2, ddy2

    IMPLICIT NONE

    INTEGER(iwp) ::  i              !< running index x direction
    INTEGER(iwp) ::  j              !< running index y direction
    INTEGER(iwp) ::  k              !< running index z direction

    REAL(wp)     ::  flag           !< flag to mask topography

!
!-- Calculate the tendency terms
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt

!
!--          Predetermine flag to mask topography
             flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

             tend(k,j,i) = tend(k,j,i) +                                       &
                         (       (                                             &
                 ( km(k,j,i)+km(k,j,i+1) ) * ( diss(k,j,i+1)-diss(k,j,i) )     &
               - ( km(k,j,i)+km(k,j,i-1) ) * ( diss(k,j,i)-diss(k,j,i-1) )     &
                                 ) * ddx2                                      &
                               + (                                             &
                 ( km(k,j,i)+km(k,j+1,i) ) * ( diss(k,j+1,i)-diss(k,j,i) )     &
               - ( km(k,j,i)+km(k,j-1,i) ) * ( diss(k,j,i)-diss(k,j-1,i) )     &
                                 ) * ddy2                                      &
                               + (                                             &
      ( km(k,j,i)+km(k+1,j,i) ) * ( diss(k+1,j,i)-diss(k,j,i) ) * ddzu(k+1)    &
                                                    * rho_air_zw(k)            &
    - ( km(k,j,i)+km(k-1,j,i) ) * ( diss(k,j,i)-diss(k-1,j,i) ) * ddzu(k)      &
                                                    * rho_air_zw(k-1)          &
                                 ) * ddzw(k) * drho_air(k)                     &
                         ) * flag * dsig_diss                                  &
                         - c_2 * diss(k,j,i)**2                                &
                               / ( e(k,j,i) + 1.0E-20_wp ) * flag

          ENDDO
       ENDDO
    ENDDO

 END SUBROUTINE diffusion_diss


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Diffusion term for the TKE dissipation rate
!> Cache-optimized version
!------------------------------------------------------------------------------!
 SUBROUTINE diffusion_diss_ij( i, j )

    USE arrays_3d,                                                             &
        ONLY:  ddzu, ddzw, drho_air, rho_air_zw

    USE grid_variables,                                                        &
        ONLY:  ddx2, ddy2

    IMPLICIT NONE

    INTEGER(iwp) ::  i         !< running index x direction
    INTEGER(iwp) ::  j         !< running index y direction
    INTEGER(iwp) ::  k         !< running index z direction

    REAL(wp)     ::  flag      !< flag to mask topography

!
!-- Calculate the mixing length (for dissipation)
    DO  k = nzb+1, nzt

!
!--    Predetermine flag to mask topography
       flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

!
!--    Calculate the tendency term
       tend(k,j,i) =  tend(k,j,i) +                                            &
                   (            (                                              &
                ( km(k,j,i)+km(k,j,i+1) ) * ( diss(k,j,i+1)-diss(k,j,i) )      &
              - ( km(k,j,i)+km(k,j,i-1) ) * ( diss(k,j,i)-diss(k,j,i-1) )      &
                                ) * ddx2                                       &
                              + (                                              &
                ( km(k,j,i)+km(k,j+1,i) ) * ( diss(k,j+1,i)-diss(k,j,i) )      &
              - ( km(k,j,i)+km(k,j-1,i) ) * ( diss(k,j,i)-diss(k,j-1,i) )      &
                                ) * ddy2                                       &
                              + (                                              &
     ( km(k,j,i)+km(k+1,j,i) ) * ( diss(k+1,j,i)-diss(k,j,i) ) * ddzu(k+1)     &
                                                   * rho_air_zw(k)             &
   - ( km(k,j,i)+km(k-1,j,i) ) * ( diss(k,j,i)-diss(k-1,j,i) ) * ddzu(k)       &
                                                   * rho_air_zw(k-1)           &
                                ) * ddzw(k) * drho_air(k)                      &
                   ) * flag * dsig_diss                                        &
                   - c_2 * diss(k,j,i)**2 / ( e(k,j,i) + 1.0E-20_wp ) * flag

    ENDDO

 END SUBROUTINE diffusion_diss_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate mixing length for LES mode.
!------------------------------------------------------------------------------!
 SUBROUTINE mixing_length_les( i, j, k, l, ll, var, var_reference )

    USE arrays_3d,                                                             &
        ONLY:  dd2zu

    USE control_parameters,                                                    &
        ONLY:  atmos_ocean_sign, g, use_single_reference_value,                &
               wall_adjustment, wall_adjustment_factor

    IMPLICIT NONE

    INTEGER(iwp) :: i   !< loop index
    INTEGER(iwp) :: j   !< loop index
    INTEGER(iwp) :: k   !< loop index

    REAL(wp)     :: dvar_dz         !< vertical gradient of var
    REAL(wp)     :: l               !< mixing length
    REAL(wp)     :: l_stable        !< mixing length according to stratification
    REAL(wp)     :: ll              !< adjusted l_grid
    REAL(wp)     :: var_reference   !< var at reference height

#if defined( __nopointer )
    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  var  !< temperature
#else
    REAL(wp), DIMENSION(:,:,:), POINTER ::  var     !< temperature
#endif

    dvar_dz = atmos_ocean_sign * ( var(k+1,j,i) - var(k-1,j,i) ) * dd2zu(k)
    IF ( dvar_dz > 0.0_wp ) THEN
       IF ( use_single_reference_value )  THEN
          l_stable = 0.76_wp * SQRT( e(k,j,i) )                                &
                             / SQRT( g / var_reference * dvar_dz ) + 1E-5_wp
       ELSE
          l_stable = 0.76_wp * SQRT( e(k,j,i) )                                &
                             / SQRT( g / var(k,j,i) * dvar_dz ) + 1E-5_wp
       ENDIF
    ELSE
       l_stable = l_grid(k)
    ENDIF
!
!-- Adjustment of the mixing length
    IF ( wall_adjustment )  THEN
       l  = MIN( wall_adjustment_factor * l_wall(k,j,i), l_grid(k), l_stable )
       ll = MIN( wall_adjustment_factor * l_wall(k,j,i), l_grid(k) )
    ELSE
       l  = MIN( l_grid(k), l_stable )
       ll = l_grid(k)
    ENDIF

 END SUBROUTINE mixing_length_les


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate mixing length for RANS mode.
!------------------------------------------------------------------------------!
 SUBROUTINE mixing_length_rans( i, j, k, l, l_diss, var, var_reference  )

    USE arrays_3d,                                                             &
        ONLY:  dd2zu

    USE control_parameters,                                                    &
        ONLY:  atmos_ocean_sign, g, use_single_reference_value

    IMPLICIT NONE

    INTEGER(iwp) :: i   !< loop index
    INTEGER(iwp) :: j   !< loop index
    INTEGER(iwp) :: k   !< loop index

    REAL(wp)     :: duv2_dz2        !< squared vertical gradient of wind vector
    REAL(wp)     :: dvar_dz         !< vertical gradient of var
    REAL(wp)     :: l               !< mixing length
    REAL(wp)     :: l_diss          !< mixing length for dissipation
    REAL(wp)     :: rif             !< Richardson flux number
    REAL(wp)     :: var_reference   !< var at reference height

#if defined( __nopointer )
    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  var  !< temperature
#else
    REAL(wp), DIMENSION(:,:,:), POINTER ::  var     !< temperature
#endif

    dvar_dz = atmos_ocean_sign * ( var(k+1,j,i) - var(k-1,j,i) ) * dd2zu(k)

    duv2_dz2 =   ( ( u(k+1,j,i) - u(k-1,j,i) ) * dd2zu(k) )**2                 &
               + ( ( v(k+1,j,i) - v(k-1,j,i) ) * dd2zu(k) )**2                 &
               + 1E-30_wp

    IF ( use_single_reference_value )  THEN
       rif = g / var_reference * dvar_dz / duv2_dz2
    ELSE
       rif = g / var(k,j,i) * dvar_dz / duv2_dz2
    ENDIF

    rif = MAX( rif, -5.0_wp )
    rif = MIN( rif,  1.0_wp )

!
!-- Calculate diabatic mixing length using Dyer-profile functions
    IF ( rif >= 0.0_wp )  THEN
       l      = MIN( l_black(k) / ( 1.0_wp + 5.0_wp * rif ), l_wall(k,j,i) )
       l_diss = l
    ELSE
!
!--    In case of unstable stratification, use mixing length of neutral case
!--    for l, but consider profile functions for l_diss
       l      = l_wall(k,j,i)
       l_diss = l * SQRT( 1.0_wp - 16.0_wp * rif )
    ENDIF

 END SUBROUTINE mixing_length_rans


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computation of the turbulent diffusion coefficients for momentum and heat
!> according to Prandtl-Kolmogorov.
!> @todo consider non-default surfaces
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_diffusivities( var, var_reference )


    USE control_parameters,                                                    &
        ONLY:  e_min, outflow_l, outflow_n, outflow_r, outflow_s

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE statistics,                                                            &
        ONLY :  rmask, sums_l_l

    USE surface_mod,                                                           &
        ONLY :  bc_h, surf_def_h, surf_def_v

    IMPLICIT NONE

    INTEGER(iwp) ::  i                   !< loop index
    INTEGER(iwp) ::  j                   !< loop index
    INTEGER(iwp) ::  k                   !< loop index
    INTEGER(iwp) ::  m                   !< loop index
    INTEGER(iwp) ::  n                   !< loop index
    INTEGER(iwp) ::  omp_get_thread_num  !< opemmp function to get thread number
    INTEGER(iwp) ::  sr                  !< statistic region
    INTEGER(iwp) ::  tn                  !< thread number

    REAL(wp)     ::  flag                !< topography flag
    REAL(wp)     ::  l                   !< mixing length
    REAL(wp)     ::  ll                  !< adjusted mixing length
    REAL(wp)     ::  var_reference       !< reference temperature

#if defined( __nopointer )
    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  var  !< temperature
#else
    REAL(wp), DIMENSION(:,:,:), POINTER ::  var  !< temperature
#endif

!
!-- Default thread number in case of one thread
    tn = 0

!
!-- Initialization for calculation of the mixing length profile
    sums_l_l = 0.0_wp

!
!-- Compute the turbulent diffusion coefficient for momentum
    !$OMP PARALLEL PRIVATE (i,j,k,l,ll,sr,tn,flag)
!$  tn = omp_get_thread_num()

!
!-- Introduce an optional minimum tke
    IF ( e_min > 0.0_wp )  THEN
       !$OMP DO
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb+1, nzt
                e(k,j,i) = MAX( e(k,j,i), e_min ) *                            &
                        MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    IF ( les_mw )  THEN
       !$OMP DO
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb+1, nzt

                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

!
!--             Determine the mixing length for LES closure
                CALL mixing_length_les( i, j, k, l, ll, var, var_reference )
!
!--             Compute diffusion coefficients for momentum and heat
                km(k,j,i) = c_0 * l * SQRT( e(k,j,i) ) * flag
                kh(k,j,i) = ( 1.0_wp + 2.0_wp * l / ll ) * km(k,j,i) * flag
!
!--             Summation for averaged profile (cf. flow_statistics)
                DO  sr = 0, statistic_regions
                   sums_l_l(k,sr,tn) = sums_l_l(k,sr,tn) + l * rmask(j,i,sr)   &
                                                             * flag
                ENDDO

             ENDDO
          ENDDO
       ENDDO

    ELSEIF ( rans_tke_l )  THEN

       !$OMP DO
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb+1, nzt

                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )
!
!--             Mixing length for RANS mode with TKE-l closure
                CALL mixing_length_rans( i, j, k, l, ll, var, var_reference )
!
!--             Compute diffusion coefficients for momentum and heat
                km(k,j,i) = c_0 * l * SQRT( e(k,j,i) ) * flag
                kh(k,j,i) = km(k,j,i) / prandtl_number * flag
!
!--             Summation for averaged profile (cf. flow_statistics)
                DO  sr = 0, statistic_regions
                   sums_l_l(k,sr,tn) = sums_l_l(k,sr,tn) + l * rmask(j,i,sr)   &
                                                             * flag
                ENDDO

             ENDDO
          ENDDO
       ENDDO

    ELSEIF ( rans_tke_e )  THEN

       !$OMP DO
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb+1, nzt

                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )
!
!--             Compute diffusion coefficients for momentum and heat
                km(k,j,i) = c_0**4 * e(k,j,i)**2 / ( diss(k,j,i) + 1.0E-30_wp ) * flag
                kh(k,j,i) = km(k,j,i) / prandtl_number * flag
!
!--             Summation for averaged profile of mixing length (cf. flow_statistics)
                DO  sr = 0, statistic_regions
                   sums_l_l(k,sr,tn) = sums_l_l(k,sr,tn) +                     &
                      c_0**3 * e(k,j,i) * SQRT(e(k,j,i)) /                     &
                      ( diss(k,j,i) + 1.0E-30_wp ) * rmask(j,i,sr) * flag
                ENDDO

             ENDDO
          ENDDO
       ENDDO

    ENDIF

    sums_l_l(nzt+1,:,tn) = sums_l_l(nzt,:,tn)   ! quasi boundary-condition for
                                                ! data output
!$OMP END PARALLEL

!
!-- Set vertical boundary values (Neumann conditions both at upward- and
!-- downward facing walls. To set wall-boundary values, the surface data type
!-- is applied.
!-- Horizontal boundary conditions at vertical walls are not set because
!-- so far vertical surfaces require usage of a Prandtl-layer where the boundary
!-- values of the diffusivities are not needed.
    IF ( .NOT. rans_tke_e )  THEN
!
!--    Upward-facing
       !$OMP PARALLEL DO PRIVATE( i, j, k )
       DO  m = 1, bc_h(0)%ns
          i = bc_h(0)%i(m)
          j = bc_h(0)%j(m)
          k = bc_h(0)%k(m)
          km(k-1,j,i) = km(k,j,i)
          kh(k-1,j,i) = kh(k,j,i)
       ENDDO
!
!--    Downward facing surfaces
       !$OMP PARALLEL DO PRIVATE( i, j, k )
       DO  m = 1, bc_h(1)%ns
          i = bc_h(1)%i(m)
          j = bc_h(1)%j(m)
          k = bc_h(1)%k(m)
          km(k+1,j,i) = km(k,j,i)
          kh(k+1,j,i) = kh(k,j,i)
       ENDDO
    ELSE
!
!--    Up- and downward facing surfaces
       DO  n = 0, 1
          DO  m = 1, surf_def_h(n)%ns
             i = surf_def_h(n)%i(m)
             j = surf_def_h(n)%j(m)
             k = surf_def_h(n)%k(m)
             km(k,j,i) = kappa * surf_def_h(n)%us(m) * dzu(k)
             kh(k,j,i) = 1.35_wp * km(k,j,i)
          ENDDO
       ENDDO
!
!--    North- and southward facing surfaces
       DO  n = 0, 1
          DO  m = 1, surf_def_v(n)%ns
             i = surf_def_v(n)%i(m)
             j = surf_def_v(n)%j(m)
             k = surf_def_v(n)%k(m)
             km(k,j,i) = kappa * surf_def_v(n)%us(m) * 0.5_wp * dy
             kh(k,j,i) = 1.35_wp * km(k,j,i)
          ENDDO
       ENDDO
!
!--    West- and eastward facing surfaces
       DO  n = 2, 3
          DO  m = 1, surf_def_v(n)%ns
             i = surf_def_v(n)%i(m)
             j = surf_def_v(n)%j(m)
             k = surf_def_v(n)%k(m)
             km(k,j,i) = kappa * surf_def_v(n)%us(m) * 0.5_wp * dx
             kh(k,j,i) = 1.35_wp * km(k,j,i)
          ENDDO
       ENDDO

       CALL exchange_horiz( km, nbgp )
       CALL exchange_horiz( kh, nbgp )

    ENDIF

!
!-- Model top
    !$OMP PARALLEL DO
    DO  i = nxlg, nxrg
       DO  j = nysg, nyng
          km(nzt+1,j,i) = km(nzt,j,i)
          kh(nzt+1,j,i) = kh(nzt,j,i)
       ENDDO
    ENDDO

!
!-- Set Neumann boundary conditions at the outflow boundaries in case of
!-- non-cyclic lateral boundaries
    IF ( outflow_l )  THEN
       km(:,:,nxl-1) = km(:,:,nxl)
       kh(:,:,nxl-1) = kh(:,:,nxl)
    ENDIF
    IF ( outflow_r )  THEN
       km(:,:,nxr+1) = km(:,:,nxr)
       kh(:,:,nxr+1) = kh(:,:,nxr)
    ENDIF
    IF ( outflow_s )  THEN
       km(:,nys-1,:) = km(:,nys,:)
       kh(:,nys-1,:) = kh(:,nys,:)
    ENDIF
    IF ( outflow_n )  THEN
       km(:,nyn+1,:) = km(:,nyn,:)
       kh(:,nyn+1,:) = kh(:,nyn,:)
    ENDIF

 END SUBROUTINE tcm_diffusivities


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Swapping of timelevels.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_swap_timelevel ( mod_count )

    IMPLICIT NONE

    INTEGER(iwp) ::  i      !< loop index x direction
    INTEGER(iwp) ::  j      !< loop index y direction
    INTEGER(iwp) ::  k      !< loop index z direction
    INTEGER, INTENT(IN) ::  mod_count  !< flag defining where pointers point to

#if defined( __nopointer )

    IF ( .NOT. constant_diffusion )  THEN
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb, nzt+1
                e(k,j,i) = e_p(k,j,i)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    IF ( rans_tke_e )  THEN
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb, nzt+1
                diss(k,j,i) = diss_p(k,j,i)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

#else

    SELECT CASE ( mod_count )

       CASE ( 0 )

          IF ( .NOT. constant_diffusion )  THEN
             e => e_1;    e_p => e_2
          ENDIF

          IF ( rans_tke_e )  THEN
             diss => diss_1;    diss_p => diss_2
          ENDIF

       CASE ( 1 )

          IF ( .NOT. constant_diffusion )  THEN
             e => e_2;    e_p => e_1
          ENDIF

          IF ( rans_tke_e )  THEN
             diss => diss_2;    diss_p => diss_1
          ENDIF

    END SELECT
#endif

 END SUBROUTINE tcm_swap_timelevel


 END MODULE turbulence_closure_mod
