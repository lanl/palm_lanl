!> @file time_integration_spinup.f90
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
! $Id: time_integration_spinup.f90 2983 2018-04-18 10:43:40Z suehring $
! Revise limitation of wall-adjacent velocity.
! 
! 2934 2018-03-26 19:13:22Z suehring
! Synchronize parent and child models after spinup.
! 
! 2881 2018-03-13 16:24:40Z maronga
! Added flag for switching on/off calculation of soil moisture
! 
! 2818 2018-02-19 16:42:36Z maronga
! Velocity components near walls/ground are now set to the profiles stored in
! u_init and v_init. Activated soil moisture calculation during spinup.
! 
! 2782 2018-02-02 11:51:10Z maronga
! Bugfix and re-activation of homogeneous setting of velocity components 
! during spinup
! 
! 2758 2018-01-17 12:55:21Z suehring
! Comment out homogeneous setting of wind velocity as this will lead to zero
! friction velocity and cause problems in MOST relationships.
! 
! 2728 2018-01-09 07:03:53Z maronga
! Set velocity componenets to homogeneous values during spinup
! 
! 2724 2018-01-05 12:12:38Z maronga
! Use dt_spinup for all active components during spinup
! 
! 2723 2018-01-05 09:27:03Z maronga
! Bugfix: array rad_sw_in no longer exists and is thus removed from RUN_CONTROL
! output.
! Added output of XY and 3D data during spinup.
! Bugfix: time step in LSM and USM was set to dt_3d instead of dt_spinup
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Added radiation interactions (moved from USM) (MS)
! 
! 2544 2017-10-13 18:09:32Z maronga
! Date and time quantities are now read from date_and_time_mod 
! 
! 2299 2017-06-29 10:14:38Z maronga
! Call of soil model adjusted to avoid prognostic equation for soil moisture
! during spinup.
! Better representation of diurnal cycle of near-surface temperature.
! Excluded prognostic equation for soil moisture during spinup.
! Added output of run control data for spinup.
! 
! 2297 2017-06-28 14:35:57Z scharf
! bugfixes
! 
! 2296 2017-06-28 07:53:56Z maronga
! Initial revision
!
!
! Description:
! ------------
!> Integration in time of the non-atmospheric model components such as land
!> surface model and urban surface model
!------------------------------------------------------------------------------!
 SUBROUTINE time_integration_spinup
 
    USE arrays_3d,                                                             &
        ONLY:  pt, pt_p, u, u_init, v, v_init

    USE control_parameters,                                                    &
        ONLY:  averaging_interval_pr, calc_soil_moisture_during_spinup,        &
               constant_diffusion, constant_flux_layer,                        &
               coupling_start_time, current_timestep_number,                   &
               data_output_during_spinup, disturbance_created, dopr_n, do_sum, &
               dt_averaging_input_pr, dt_dopr, dt_dots, dt_do2d_xy, dt_do3d,   &
               dt_run_control, dt_spinup, dt_3d, humidity,                     &
               intermediate_timestep_count,                                    &
               intermediate_timestep_count_max, land_surface,                  &
               simulated_time, simulated_time_chr,                             &
               skip_time_dopr, skip_time_do2d_xy, skip_time_do3d, spinup,      &
               spinup_pt_amplitude, spinup_pt_mean, spinup_time,               &
               timestep_count, timestep_scheme, time_dopr, time_dopr_av,       &
               time_dots, time_do2d_xy, time_do3d, time_run_control,           &
               time_since_reference_point, ug_surface, vg_surface, urban_surface

    USE constants,                                                             &
        ONLY:  pi

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s

    USE date_and_time_mod,                                                     &
        ONLY: day_of_year_init, time_utc_init

    USE indices,                                                               &
        ONLY:  nbgp, nzb, nzt, nysg, nyng, nxlg, nxrg

    USE pegrid

    USE kinds

    USE statistics,                                                            &
        ONLY:  flow_statistics_called

    USE surface_layer_fluxes_mod,                                              &
        ONLY:  surface_layer_fluxes

    USE surface_mod,                                                           &
        ONLY :  surf_def_h, surf_def_v


    IMPLICIT NONE

    CHARACTER (LEN=9) ::  time_to_string          !< 
  
    INTEGER(iwp) ::  i !< running index
    INTEGER(iwp) ::  j !< running index
    INTEGER(iwp) ::  k !< running index
    INTEGER(iwp) ::  l !< running index
    INTEGER(iwp) ::  m !< running index

    INTEGER(iwp) :: current_timestep_number_spinup = 0  !< number if timestep during spinup
  
    LOGICAL :: run_control_header_spinup = .FALSE.  !< flag parameter for steering whether the header information must be output

    REAL(wp) ::  pt_spinup   !< temporary storage of temperature
    REAL(wp) ::  dt_save     !< temporary storage for time step
                  
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pt_save  !< temporary storage of temperature
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_save   !< temporary storage of u wind component
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_save   !< temporary storage of v wind component


!
!-- Save 3D arrays because they are to be changed for spinup purpose
    ALLOCATE( pt_save(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( u_save(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( v_save(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

    CALL exchange_horiz( pt, nbgp )   
    CALL exchange_horiz( u,  nbgp )  
    CALL exchange_horiz( v,  nbgp )  
 
    pt_save = pt
    u_save  = u
    v_save  = v

!
!-- Set the same wall-adjacent velocity to all grid points. The sign of the 
!-- original velocity field must be preserved because the surface schemes crash 
!-- otherwise. The precise reason is still unknown. A minimum velocity of 0.1
!-- m/s is used to maintain turbulent transfer at the surface.
    CALL exchange_horiz( u,  nbgp )
    CALL exchange_horiz( v,  nbgp )

    dt_save = dt_3d
    dt_3d   = dt_spinup

    CALL location_message( 'starting spinup-sequence', .TRUE. )
!
!-- Start of the time loop
    DO  WHILE ( simulated_time < spinup_time )

       CALL cpu_log( log_point_s(15), 'timesteps spinup', 'start' )
   
!
!--    Start of intermediate step loop
       intermediate_timestep_count = 0
       DO  WHILE ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )

          intermediate_timestep_count = intermediate_timestep_count + 1

!
!--       Set the steering factors for the prognostic equations which depend
!--       on the timestep scheme
          CALL timestep_scheme_steering


!
!--       Estimate a near-surface air temperature based on the position of the 
!--       sun and user input about mean temperature and amplitude. The time is
!--       shifted by one hour to simulate a lag between air temperature and
!--       incoming radiation
          pt_spinup = spinup_pt_mean 


!
          CALL exchange_horiz( pt,  nbgp )    


!
!--       Swap the time levels in preparation for the next time step.
          timestep_count = timestep_count + 1
     
!
!--       Compute the diffusion quantities
          IF ( .NOT. constant_diffusion )  THEN

!
!--          First the vertical (and horizontal) fluxes in the surface 
!--          (constant flux) layer are computed
             IF ( constant_flux_layer )  THEN
                CALL cpu_log( log_point(19), 'surface_layer_fluxes', 'start' )
                CALL surface_layer_fluxes
                CALL cpu_log( log_point(19), 'surface_layer_fluxes', 'stop' )
             ENDIF

!
          ENDIF
       ENDDO   ! Intermediate step loop

!
!--    Increase simulation time and output times
       current_timestep_number_spinup = current_timestep_number_spinup + 1
       simulated_time             = simulated_time   + dt_3d
       simulated_time_chr         = time_to_string( simulated_time )
       time_since_reference_point = simulated_time - coupling_start_time

       IF ( data_output_during_spinup )  THEN
          IF ( simulated_time >= skip_time_do3d    )  THEN
             time_do3d          = time_do3d        + dt_3d
          ENDIF
          time_dots          = time_dots        + dt_3d
          IF ( simulated_time >= skip_time_dopr )  THEN
             time_dopr       = time_dopr        + dt_3d
          ENDIF
          time_run_control   = time_run_control + dt_3d

!
!--       Carry out statistical analysis and output at the requested output times.
!--       The MOD function is used for calculating the output time counters (like
!--       time_dopr) in order to regard a possible decrease of the output time
!--       interval in case of restart runs

!
!--       Set a flag indicating that so far no statistics have been created
!--       for this time step
          flow_statistics_called = .FALSE.

!
!--       If required, call flow_statistics for averaging in time
          IF ( averaging_interval_pr /= 0.0_wp  .AND.                          &
             ( dt_dopr - time_dopr ) <= averaging_interval_pr  .AND.           &
             simulated_time >= skip_time_dopr )  THEN
             time_dopr_av = time_dopr_av + dt_3d
             IF ( time_dopr_av >= dt_averaging_input_pr )  THEN
                do_sum = .TRUE.
                time_dopr_av = MOD( time_dopr_av,                              &
                               MAX( dt_averaging_input_pr, dt_3d ) )
             ENDIF
          ENDIF
          IF ( do_sum )  CALL flow_statistics

!
!--       Output of profiles
          IF ( time_dopr >= dt_dopr )  THEN
             IF ( dopr_n /= 0 )  CALL data_output_profiles
             time_dopr = MOD( time_dopr, MAX( dt_dopr, dt_3d ) )
             time_dopr_av = 0.0_wp    ! due to averaging (see above)
          ENDIF

!
!--       Output of time series
          IF ( time_dots >= dt_dots )  THEN
             CALL data_output_tseries
             time_dots = MOD( time_dots, MAX( dt_dots, dt_3d ) )
          ENDIF
!
!--       3d-data output (volume data)
          IF ( time_do3d >= dt_do3d )  THEN
             CALL data_output_3d( 0 )
             time_do3d = MOD( time_do3d, MAX( dt_do3d, dt_3d ) )
          ENDIF


       ENDIF
       CALL cpu_log( log_point_s(15), 'timesteps spinup', 'stop' )


!
!--    Run control output
       IF ( myid == 0 )  THEN
!
!--       If necessary, write header
          IF ( .NOT. run_control_header_spinup )  THEN
             CALL check_open( 15 )
             WRITE ( 15, 100 )
             run_control_header_spinup = .TRUE.
          ENDIF
!
!--       Write some general information about the spinup in run control file
          WRITE ( 15, 101 )  current_timestep_number_spinup, simulated_time_chr, dt_3d, pt_spinup
!
!--       Write buffer contents to disc immediately
          FLUSH( 15 )
       ENDIF



    ENDDO   ! time loop

!
!-- Write back saved arrays to the 3D arrays
    pt   = pt_save
    pt_p = pt_save
    u    = u_save
    v    = v_save

!
!-- Reset time step
    dt_3d = dt_save

    DEALLOCATE(pt_save)
    DEALLOCATE(u_save)
    DEALLOCATE(v_save)

    CALL location_message( 'finished spinup-sequence', .TRUE. )


!
!-- Formats
100 FORMAT (///'Spinup control output:'/  &
            '--------------------------------'// &
            'ITER.  HH:MM:SS    DT   PT(z_MO)'/   &
            '--------------------------------')
101 FORMAT (I5,2X,A9,1X,F6.2,3X,F6.2,2X,F6.2)
 END SUBROUTINE time_integration_spinup
