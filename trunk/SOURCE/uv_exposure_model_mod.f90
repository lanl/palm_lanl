!> @file uv_exposure_model_mod.f90
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
!------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: uv_exposure_model_mod.f90 3014 2018-05-09 08:42:38Z maronga $
! Bugfix: domain bounds of local_pf corrected
! 
! 3004 2018-04-27 12:33:25Z Giersch
! Further allocation checks implemented (averaged data will be assigned to fill
! values if no allocation happened so far) 
! 
! 2932 2018-03-26 09:39:22Z maronga
! renamed uvexposure_par to biometeorology_parameters
! 
! 2894 2018-03-15 09:17:58Z Giersch
! Routine for skipping global restart data has been removed, uvem_last_actions
! has been renamed to uvem_wrd_global and uvem_read_restart_data has been 
! renamed to uvem_rrd_global, variable named found has been introduced for 
! checking if restart data was found, reading of restart strings has been moved
! completely to read_restart_data_mod, marker *** end new module *** has been
! removed, strings and their respective lengths are written out and read now 
! in case of restart runs to get rid of prescribed character lengths, CASE 
! DEFAULT was added if restart data is read 
! 
! 2848 2018-03-05 10:11:18Z Giersch
! Initial revision
!
! 
!
! Authors:
! --------
! @author Michael Schrempf
!
!
! Description:
! ------------
!> Calculation of vitamin-D weighted UV exposure
!>
!>
!> @todo uv_vitd3dose-->new output type necessary (cumulative)
!> @todo consider upwelling radiation
!>
!> @note <Enter notes on the module>
!>
!> @bug  <Enter known bugs here>
!------------------------------------------------------------------------------!
 MODULE uv_exposure_model_mod
 

    USE kinds

!
!-- Load required variables from existing modules
!-- USE modulename,                                                            &
!--     ONLY: ...


    IMPLICIT NONE


!
!-- Declare all global variables within the module (alphabetical order)
    INTEGER(iwp) ::  ai                      = 0 !< running index
    INTEGER(iwp) ::  clothing                = 3 !< clothing (0=unclothed, 1=Arms,Hands Face free, 3=Hand, Face free)
    INTEGER(iwp) ::  consider_obstructions   = 1 !< 0 = unobstructed scenario, 
                                                 !< 1 = scenario where obstrcutions are considered
    INTEGER(iwp) ::  orientation_angle       = 0 !< orientation of front/face of the human model
    INTEGER(iwp) ::  saa_in_south            = 1 !< 0 = sun always in south direction, 1 = actual sun position  
    INTEGER(iwp) ::  obstruction_direct_beam = 0 !< Obstruction information for direct beam
    INTEGER(iwp) ::  turn_to_sun             = 1 !< 0 = orientation of human as in orientation_angle,
                                                 !< 1 = human is always orientated towards the sun
    INTEGER(iwp) ::  zi                      = 0 !< running index

    INTEGER(iwp), DIMENSION(0:44)  ::  obstruction_temp1 = 0 !< temporary obstruction information
                                                             !< stored with ibset as logical information 
    INTEGER(iwp), DIMENSION(0:359) ::  obstruction_temp2 = 0 !< temporary obstruction information restored values
                                                             !< from logical information, which where stored by ibset  
    
    INTEGER(iwp), DIMENSION(0:35,0:9) ::  obstruction = 0 !< final obstruction information array for all
                                                          !< hemispherical directions
    
    INTEGER(KIND=1), DIMENSION(:,:,:), ALLOCATABLE ::  obstruction_lookup_table !< lookup table of obstruction
                                                                                !< information of entire domain

    REAL(wp) ::  diffuse_exposure            = 0.0_wp !< calculated exposure by diffuse radiation
    REAL(wp) ::  direct_exposure             = 0.0_wp !< calculated exposure by direct beam   
    REAL(wp) ::  projection_area_direct_beam = 0.0_wp !< projection area for  by direct beam
    REAL(wp) ::  saa                         = 180.0_wp !< solar azimuth angle
    REAL(wp) ::  startpos_human              = 0.0_wp !< start position in azimuth direction for the
                                                      !< interpolation of the projection area              
    REAL(wp) ::  startpos_saa_float          = 0.0_wp !< start position in azimuth direction for the
                                                      !< interpolation of the radiance field             
    REAL(wp) ::  sza                         = 20.0_wp !< solar zenith angle
    REAL(wp) ::  xfactor                     = 0.0_wp !< relative x-position used for interpolation
    REAL(wp) ::  yfactor                     = 0.0_wp !< relative y-position used for interpolation
    
    REAL(wp), DIMENSION(0:2)  ::  irradiance  = 0.0_wp !< iradiance values extracted from irradiance lookup table
    
    REAL(wp), DIMENSION(0:2,0:90) ::  irradiance_lookup_table = 0.0_wp !< irradiance lookup table contains values
                                                                       !< for direct, diffuse and global component
    REAL(wp), DIMENSION(0:35,0:9) ::  integration_array            = 0.0_wp
    REAL(wp), DIMENSION(0:35,0:9) ::  projection_area              = 0.0_wp
    REAL(wp), DIMENSION(0:35,0:9) ::  projection_area_lookup_table = 0.0_wp
    REAL(wp), DIMENSION(0:71,0:9) ::  projection_area_temp         = 0.0_wp
    REAL(wp), DIMENSION(0:35,0:9) ::  radiance_array               = 0.0_wp
    REAL(wp), DIMENSION(0:71,0:9) ::  radiance_array_temp          = 0.0_wp
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vitd3_exposure
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vitd3_exposure_av

    REAL(wp), DIMENSION(0:35,0:9,0:90) ::  radiance_lookup_table = 0.0_wp


    SAVE

    PRIVATE

    
!
!-- Add INTERFACES that must be available to other modules (alphabetical order)
    PUBLIC uvem_3d_data_averaging, uvem_calc_exposure, uvem_check_data_output, &
           uvem_data_output_2d, uvem_define_netcdf_grid, uvem_init,            &
           uvem_init_arrays, uvem_parin

!
!-- Add VARIABLES that must be available to other modules (alphabetical order)
!     PUBLIC 

!
!-- Add PROGNOSTIC QUANTITIES that must be available to other modules (alphabetical order)
!-- PUBLIC ...


!
!-- Default procedures for all new modules (not all are necessarily required, 
!-- alphabetical order is not essential)

!
!-- PALM interfaces:
!-- Data output checks for 2D/3D data to be done in check_parameters
    INTERFACE uvem_check_data_output
       MODULE PROCEDURE uvem_check_data_output
    END INTERFACE uvem_check_data_output
    
! !
! !-- Data output checks for profile data to be done in check_parameters
!     INTERFACE uvem_check_data_output_pr
!        MODULE PROCEDURE uvem_check_data_output_pr
!     END INTERFACE uvem_check_data_output_pr
!     
! !
! !-- Input parameter checks to be done in check_parameters
!     INTERFACE uvem_check_parameters
!        MODULE PROCEDURE uvem_check_parameters
!     END INTERFACE uvem_check_parameters
! 
!
!-- Averaging of 3D data for output
    INTERFACE uvem_3d_data_averaging
       MODULE PROCEDURE uvem_3d_data_averaging
    END INTERFACE uvem_3d_data_averaging
! 
! !
!-- Data output of 2D quantities
    INTERFACE uvem_data_output_2d
       MODULE PROCEDURE uvem_data_output_2d
    END INTERFACE uvem_data_output_2d
! 
! !
! !-- Data output of 3D data
!     INTERFACE uvem_data_output_3d
!        MODULE PROCEDURE uvem_data_output_3d
!     END INTERFACE uvem_data_output_3d
! 
! !
! !-- Definition of data output quantities
    INTERFACE uvem_define_netcdf_grid
       MODULE PROCEDURE uvem_define_netcdf_grid
    END INTERFACE uvem_define_netcdf_grid
! 
! !
! !-- Output of information to the header file
!     INTERFACE uvem_header
!        MODULE PROCEDURE uvem_header
!     END INTERFACE uvem_header
  
!
!-- Initialization actions  
    INTERFACE uvem_init
       MODULE PROCEDURE uvem_init
    END INTERFACE uvem_init
 
!
!-- Initialization of arrays
    INTERFACE uvem_init_arrays
       MODULE PROCEDURE uvem_init_arrays
    END INTERFACE uvem_init_arrays
! 
! !
! !-- Writing of binary output for restart runs  !!! renaming?!
!     INTERFACE uvem_wrd_global
!        MODULE PROCEDURE uvem_wrd_global
!     END INTERFACE uvem_wrd_global
!     
!
!-- Reading of NAMELIST parameters
    INTERFACE uvem_parin
       MODULE PROCEDURE uvem_parin
    END INTERFACE uvem_parin
! 
! !
! !-- Reading of parameters for restart runs
!     INTERFACE uvem_rrd_global
!        MODULE PROCEDURE uvem_rrd_global
!     END INTERFACE uvem_rrd_global
! 
! !
! !-- Swapping of time levels (required for prognostic variables)
!     INTERFACE uvem_swap_timelevel
!        MODULE PROCEDURE uvem_swap_timelevel
!     END INTERFACE uvem_swap_timelevel
! 
! !
! !-- New module-specific procedure(s) (alphabetical order):
!     INTERFACE uvem_newprocedure
!        MODULE PROCEDURE uvem_newprocedure
!     END INTERFACE uvem_newprocedure

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output for new module.
!------------------------------------------------------------------------------!
 SUBROUTINE uvem_check_data_output( var, unit, i, ilen, k )
 
    USE control_parameters,                                                    &
        ONLY:  data_output, message_string, uv_exposure

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  unit     !< string for unit of output quantity
    CHARACTER (LEN=*) ::  var      !< string for output quantity

    INTEGER(iwp) ::  i      !<
    INTEGER(iwp) ::  ilen   !<   
    INTEGER(iwp) ::  k      !<

    SELECT CASE ( TRIM( var ) )


       CASE ( 'uvem_vitd3*' )
          IF (  .NOT.  uv_exposure )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                      'res a namelist &uvexposure_par'
             CALL message( 'uvem_check_data_output', 'UV0001', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( k == 0  .OR.  data_output(i)(ilen-2:ilen) /= '_xy' )  THEN
             message_string = 'illegal value for data_output: "' //            &
                              TRIM( var ) // '" & only 2d-horizontal ' //      &
                              'cross sections are allowed for this value'
             CALL message( 'check_parameters', 'PA0111', 1, 2, 0, 6, 0 )
          ENDIF
          unit = 'IU/s'

       CASE ( 'uvem_vitd3dose*' )
          IF (  .NOT.  uv_exposure )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                      'res  a namelist &uvexposure_par'
             CALL message( 'uvem_check_data_output', 'UV0001', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( k == 0  .OR.  data_output(i)(ilen-2:ilen) /= '_xy' )  THEN
             message_string = 'illegal value for data_output: "' //            &
                              TRIM( var ) // '" & only 2d-horizontal ' //      &
                              'cross sections are allowed for this value'
             CALL message( 'check_parameters', 'PA0111', 1, 2, 0, 6, 0 )
          ENDIF
          unit = 'IU/av-h'
             
       CASE DEFAULT
          unit = 'illegal'

    END SELECT

 END SUBROUTINE uvem_check_data_output


!-----------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 2D output variables
!-----------------------------------------------------------------------------!
 SUBROUTINE uvem_data_output_2d( av, variable, found, grid, mode, local_pf,   &
                                 two_d, nzb_do, nzt_do )
 
    USE indices

    USE kinds


    IMPLICIT NONE

    CHARACTER (LEN=*) ::  grid     !< 
    CHARACTER (LEN=*) ::  mode     !< 
    CHARACTER (LEN=*) ::  variable !< 

    INTEGER(iwp) ::  av      !< 
    INTEGER(iwp) ::  i       !< running index 
    INTEGER(iwp) ::  j       !< running index
    INTEGER(iwp) ::  k       !< running index
    INTEGER(iwp) ::  m       !< running index
    INTEGER(iwp) ::  nzb_do  !< 
    INTEGER(iwp) ::  nzt_do  !< 

    LOGICAL      ::  found !< 
    LOGICAL      ::  two_d !< flag parameter that indicates 2D variables (horizontal cross sections)

    REAL(wp) ::  fill_value = -999.0_wp    !< value for the _FillValue attribute

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf !< 


    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )
!
!--    Before data is transfered to local_pf, transfer is it 2D dummy variable and exchange ghost points therein. 
!--    However, at this point this is only required for instantaneous arrays, time-averaged quantities are already exchanged. 
       CASE ( 'uvem_vitd3*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = vitd3_exposure(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'uvem_vitd3dose*_xy' )        ! 2d-array
          IF ( .NOT. ALLOCATED( vitd3_exposure_av ) ) THEN
             ALLOCATE( vitd3_exposure_av(nysg:nyng,nxlg:nxrg) )
             vitd3_exposure_av = REAL( fill_value, KIND = wp )
          ENDIF
          IF ( av == 1 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = vitd3_exposure_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT
 
 END SUBROUTINE uvem_data_output_2d


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!------------------------------------------------------------------------------!
 SUBROUTINE uvem_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )
    
    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN)  ::  var         !< 
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_x      !< 
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_y      !< 
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_z      !< 

    LOGICAL, INTENT(OUT)           ::  found       !<

    found  = .TRUE.

!
!-- Check for the grid
    SELECT CASE ( TRIM( var ) )

       CASE ( 'uvem_vitd3*_xy', 'uvem_vitd3dose*_xy' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu1'

       CASE DEFAULT
          found  = .FALSE.
          grid_x = 'none'
          grid_y = 'none'
          grid_z = 'none'

    END SELECT

 END SUBROUTINE uvem_define_netcdf_grid


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &uvexposure_par for UV exposure model
!------------------------------------------------------------------------------!
 SUBROUTINE uvem_parin

    USE control_parameters,                                                   &
        ONLY:  uv_exposure


    IMPLICIT NONE

    CHARACTER (LEN=80) ::  line  !< dummy string for current line in parameter file 
    
    NAMELIST /biometeorology_parameters/  clothing
    
    line = ' '
    
!
!-- Try to find uv exposure model namelist
    REWIND ( 11 )
    line = ' '
    DO   WHILE ( INDEX( line, '&biometerology_parameters' ) == 0 )
       READ ( 11, '(A)', END=10 )  line
    ENDDO
    BACKSPACE ( 11 )

!
!-- Read user-defined namelist
    READ ( 11, biometeorology_parameters )

!
!-- Set flag that indicates that the uv exposure model is switched on
    uv_exposure = .TRUE.


 10 CONTINUE
       

 END SUBROUTINE uvem_parin


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine for averaging 3D data
!------------------------------------------------------------------------------!
 SUBROUTINE uvem_3d_data_averaging( mode, variable )
 

    USE control_parameters

    USE indices

    USE kinds

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  mode    !< 
    CHARACTER (LEN=*) :: variable !< 

    INTEGER(iwp) ::  i       !< 
    INTEGER(iwp) ::  j       !< 
    INTEGER(iwp) ::  k       !< 
    INTEGER(iwp) ::  m       !< running index

    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'uvem_vitd3dose*' )
             IF ( .NOT. ALLOCATED( vitd3_exposure_av ) )  THEN
                ALLOCATE( vitd3_exposure_av(nysg:nyng,nxlg:nxrg) )
             ENDIF
             vitd3_exposure_av = 0.0_wp


          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'uvem_vitd3dose*' )
             IF ( ALLOCATED( vitd3_exposure_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      vitd3_exposure_av(j,i) = vitd3_exposure_av(j,i) + vitd3_exposure(j,i)
                   ENDDO
                ENDDO
             ENDIF


          CASE DEFAULT
             CONTINUE

       END SELECT

!
!-- No averaging since we are calculating a dose (only sum is calculated and saved to
!-- av.nc file)
!    ELSEIF ( mode == 'average' )  THEN

    ENDIF

 END SUBROUTINE uvem_3d_data_averaging

    


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the new module
!------------------------------------------------------------------------------!
 SUBROUTINE uvem_init
    

    USE control_parameters,                                                   &
        ONLY:  initializing_actions

    IMPLICIT NONE

!
!-- Actions for initial runs
    IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
       OPEN(90, STATUS='old',FILE=&
       'RADIANCE', FORM='UNFORMATTED')
       READ(90) radiance_lookup_table
       CLOSE(90)
       OPEN(90, STATUS='old',FILE=&
       'IRRADIANCE', FORM='UNFORMATTED')
       READ(90) irradiance_lookup_table
       CLOSE(90)
      
      !________________________LOAD Obstruction information _______________________________
       IF ( consider_obstructions == 1 )  THEN
          OPEN(90, STATUS='old',FILE=&
          'OBSTRUCTION', FORM='UNFORMATTED')
          READ(90) obstruction_lookup_table
          CLOSE(90)  
      ELSEIF( consider_obstructions == 0 ) THEN
          obstruction(:,:) = 1
      ENDIF

      !________________________LOAD integration_array ________________________________________________________________________
      OPEN(90, STATUS='old',FILE=&
      'SOLIDANGLE', FORM='UNFORMATTED')
      READ(90) integration_array
      CLOSE(90)

      IF (clothing==1) THEN
          OPEN(90, STATUS='old',FILE='HUMAN', FORM='UNFORMATTED')
      ENDIF
      READ(90) projection_area_lookup_table
      CLOSE(90)

!
!-- Actions for restart runs
    ELSE
! ....

    ENDIF

 END SUBROUTINE uvem_init


!-----------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate new module arrays and define pointers if required
!-----------------------------------------------------------------------------!
 SUBROUTINE uvem_init_arrays
    
    USE indices,                                                              &
        ONLY:  nxlg, nxrg, nyng, nysg


    IMPLICIT NONE

!
!-- Allocate arrays
    ALLOCATE ( vitd3_exposure(nysg:nyng,nxlg:nxrg) )
    ALLOCATE ( vitd3_exposure_av(nysg:nyng,nxlg:nxrg) )
    ALLOCATE ( obstruction_lookup_table(nxlg:nxrg,nysg:nyng,0:44) )

!
!-- Initialize arrays
    vitd3_exposure = 0.0_wp
    vitd3_exposure_av = 0.0_wp
    obstruction_lookup_table = 0


 END SUBROUTINE uvem_init_arrays

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Module-specific routine for new module
!------------------------------------------------------------------------------!
 SUBROUTINE uvem_solar_position 

    USE constants,                                                             &
       ONLY:  pi
    
    USE date_and_time_mod,                                                     &
       ONLY:  calc_date_and_time, day_of_year, time_utc 
    
    USE control_parameters,                                                    &
       ONLY:  latitude, longitude   

    IMPLICIT NONE
    
   
    REAL(wp) ::  alpha       = 0.0_wp   !< solar azimuth angle in radiant       
    REAL(wp) ::  declination = 0.0_wp   !< declination
    REAL(wp) ::  dtor        = 0.0_wp   !< factor to convert degree to radiant
    REAL(wp) ::  js          = 0.0_wp   !< parameter for solar position calculation
    REAL(wp) ::  lat         = 52.39_wp !< latitude
    REAL(wp) ::  lon         = 9.7_wp   !< longitude        
    REAL(wp) ::  thetar      = 0.0_wp   !< angle for solar zenith angle calculation
    REAL(wp) ::  thetasr     = 0.0_wp   !< angle for solar azimuth angle calculation    
    REAL(wp) ::  zgl         = 0.0_wp   !< calculated exposure by direct beam   
    REAL(wp) ::  woz         = 0.0_wp   !< calculated exposure by diffuse radiation
    REAL(wp) ::  wsp         = 0.0_wp   !< calculated exposure by direct beam   
   

    CALL calc_date_and_time
       
    dtor = pi / 180.0_wp
    lat = latitude
    lon = longitude
!
!-- calculation of js :
    js=  72.0_wp * ( day_of_year + ( time_utc / 86400.0_wp ) ) / 73.0_wp 
!
!-- calculation of equation of time (zgl):
    zgl = 0.0066_wp + 7.3525_wp * cos( ( js + 85.9_wp ) * dtor ) + 9.9359_wp *                     &
    cos( ( 2.0_wp * js + 108.9_wp ) * dtor ) + 0.3387_wp * cos( ( 3 * js + 105.2_wp ) * dtor )
!
!-- calculation of apparent solar time woz:
    woz = ( ( time_utc / 3600.0_wp ) - ( 4.0_wp * ( 15.0_wp - lon ) ) / 60.0_wp ) + ( zgl / 60.0_wp )
!
!-- calculation of hour angle (wsp):
    wsp = ( woz - 12.0_wp ) * 15.0_wp
!
!-- calculation of declination:
    declination = 0.3948_wp - 23.2559_wp * cos( ( js + 9.1_wp ) * dtor ) -                         &
    0.3915_wp * cos( ( 2.0_wp * js + 5.4_wp ) * dtor ) - 0.1764_wp * cos( ( 3.0_wp * js + 26.0_wp ) * dtor )
!
!-- calculation of solar zenith angle
    thetar  = acos( sin( lat * dtor) * sin( declination * dtor ) + cos( wsp * dtor ) *             &
    cos( lat * dtor ) * cos( declination * dtor ) )
    thetasr = asin( sin( lat * dtor) * sin( declination * dtor ) + cos( wsp * dtor ) *             & 
    cos( lat * dtor ) * cos( declination * dtor ) )
    sza = thetar / dtor
!
!-- calculation of solar azimuth angle
    IF (woz .LE. 12.0_wp) alpha = pi - acos( sin(thetasr) * sin(lat * dtor) -                           & 
    sin( declination * dtor ) / ( cos(thetasr) * cos( lat * dtor ) ) )    
    IF (woz .GT. 12.0_wp) alpha = pi + acos( sin(thetasr) * sin(lat * dtor) -                           &
    sin( declination * dtor ) / ( cos(thetasr) * cos( lat * dtor ) ) )    
    saa = alpha / dtor

 END SUBROUTINE uvem_solar_position


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Module-specific routine for new module
!------------------------------------------------------------------------------!
 SUBROUTINE uvem_calc_exposure
    
    USE constants,                                                             &
        ONLY:  pi

    USE indices,                                                               &
        ONLY:  nxlg, nxrg, nyng, nysg
    
    
    IMPLICIT NONE   
    
    INTEGER(iwp) ::  i   !< running index
    INTEGER(iwp) ::  j   !< running index


    CALL uvem_solar_position
     
    IF (sza .GE. 90) THEN
       vitd3_exposure(:,:) = 0.0_wp
    ELSE
       
!       
!--    rotate 3D-Model human to desired direction  ----------------------------- 
       projection_area_temp( 0:35,:) = projection_area_lookup_table
       projection_area_temp(36:71,:) = projection_area_lookup_table               
       IF ( Turn_to_Sun .EQ. 0 )  startpos_human = orientation_angle / 10.0_wp
       IF ( Turn_to_Sun .EQ. 1 )  startpos_human = startpos_saa_float        
       DO  ai = 0, 35
           xfactor = ( startpos_human ) - INT( startpos_human )
          DO  zi = 0, 9
              projection_area(ai,zi) = ( projection_area_temp( ai +     INT( startpos_human ), zi) * &
                                       (1.0_wp - xfactor ) ) &
                                      +( projection_area_temp( ai + 1 + INT( startpos_human ), zi) * &
                                        xfactor)
          ENDDO
       ENDDO
!      
!--    calculate Projectionarea for direct beam -----------------------------'
       projection_area_temp( 0:35,:) = projection_area
       projection_area_temp(36:71,:) = projection_area
       yfactor = ( sza / 10.0_wp ) - INT( sza / 10.0_wp )
       xfactor = ( startpos_saa_float ) - INT( startpos_saa_float )
       projection_area_direct_beam = ( projection_area_temp( INT(startpos_saa_float)    ,INT( sza / 10.0_wp ) ) * &
                                     ( 1.0_wp - xfactor ) * ( 1.0_wp - yfactor ) ) + &
                                     ( projection_area_temp( INT(startpos_saa_float) + 1,INT(sza/10.0_wp))  *&
                                     (    xfactor ) * ( 1.0_wp - yfactor ) ) + &
                                     ( projection_area_temp( INT(startpos_saa_float)    ,INT(sza/10.0_wp)+1)*&
                                     ( 1.0_wp - xfactor ) * (   yfactor ) ) + &
                                     ( projection_area_temp( INT(startpos_saa_float) + 1,INT(sza/10.0_wp)+1)*&
                                     (    xfactor ) * (   yfactor ) )
                 
             
!           
!--    interpolate to accurate Solar Zenith Angle  ------------------          
       DO  ai = 0, 35
          xfactor = (sza)-INT(sza)
          DO  zi = 0, 9
              radiance_array(ai,zi) = ( radiance_lookup_table(ai, zi, INT(sza) ) * ( 1.0_wp - xfactor) ) +&
              ( radiance_lookup_table(ai,zi,INT(sza) + 1) * xfactor )
          ENDDO
       ENDDO
       Do  ai = 0, 2            
           irradiance(ai) = ( irradiance_lookup_table(ai, INT(sza) ) * ( 1.0_wp - xfactor)) + &
           (irradiance_lookup_table(ai, INT(sza) + 1) * xfactor )
       ENDDO   
!         
!--    interpolate to accurate Solar Azimuth Angle ------------------
       IF ( saa_in_south .EQ. 0 )  THEN
          startpos_saa_float = 180.0_wp / 10.0_wp
       ELSE 
          startpos_saa_float = saa / 10.0_wp
       ENDIF
       radiance_array_temp( 0:35,:) = radiance_array
       radiance_array_temp(36:71,:) = radiance_array
       xfactor = (startpos_saa_float) - INT(startpos_saa_float)
       DO  ai = 0, 35
          DO  zi = 0, 9
             radiance_array(ai,zi) = ( radiance_array_temp( ai + INT( startpos_saa_float ), zi ) * &
                                     (1.0_wp - xfactor)) &
                                     + ( radiance_array_temp( ai + 1 + INT( startpos_saa_float ), zi ) &
                                     * xfactor)
          ENDDO
       ENDDO 



                                    
       DO  i = nxlg, nxrg    
          DO  j = nysg, nyng
!            
!--          extract obstrcution from IBSET-Integer_Array ------------------'
             IF (consider_obstructions == 1) THEN
                obstruction_temp1 = obstruction_lookup_table(i,j,:)
                IF (obstruction_temp1(0) .NE. 9) THEN 
                   DO  zi = 0, 44 
                      DO  ai = 0, 7 
                         IF ( btest( obstruction_temp1(zi), ai ) .EQV. .TRUE.) THEN
                            obstruction_temp2( ( zi * 8 ) + ai ) = 1
                         ELSE
                            obstruction_temp2( ( zi * 8 ) + ai ) = 0
                         ENDIF
                      ENDDO
                   ENDDO        
                   DO  zi = 0, 9                                          
                      obstruction(:,zi) = obstruction_temp2( zi * 36 :( zi * 36) + 35 )
                   ENDDO
                ELSE 
                   obstruction(:,:) = 0
                ENDIF
             ENDIF
!              
!--          calculated human exposure ------------------'  
             diffuse_exposure = SUM( radiance_array * projection_area * integration_array * obstruction )      
         
             obstruction_direct_beam = obstruction( nint(startpos_saa_float), nint( sza / 10.0_wp ) ) 
             IF (sza .GE. 89.99_wp) THEN
                sza = 89.99999_wp
             ENDIF
!             
!--          calculate direct normal irradiance (direct beam) ------------------'
             direct_exposure = ( irradiance(1) / cos( pi * sza / 180.0_wp ) ) * &
             projection_area_direct_beam * obstruction_direct_beam  
               
             vitd3_exposure(j,i) = ( diffuse_exposure + direct_exposure ) / 1000.0_wp * 70.97_wp 
             ! unit = international units vitamin D per second              
          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE uvem_calc_exposure




 
! ------------------------------------------------------------------------------!
! Description:
! ------------
! > Check parameters routine
! ------------------------------------------------------------------------------!
!  SUBROUTINE uvem_check_parameters
! 
!     USE control_parameters,                                                    &
!         ONLY:  message_string
! 
!        
!     IMPLICIT NONE
! 
! 
! --    Checks go here (cf. check_parameters.f90). 
!            
!  END SUBROUTINE uvem_check_parameters


 
! !------------------------------------------------------------------------------!
! ! Description:
! ! ------------
! !> Header output
! !------------------------------------------------------------------------------!
!  SUBROUTINE uvem_header ( io )
! 
! 
!     IMPLICIT NONE
! 
!  
!     INTEGER(iwp), INTENT(IN) ::  io   !< Unit of the output file
!  
! !
! !-- Header output goes here
! !-- ...
! 
!  END SUBROUTINE uvem_header



! !------------------------------------------------------------------------------!
! ! Description:
! ! ------------
! !> This routine reads the global restart data.
! !------------------------------------------------------------------------------!
!  SUBROUTINE uvem_rrd_global
!
!
!     USE control_parameters,                                                    &
!         ONLY: length, restart_string
!
!
!     IMPLICIT NONE
!
!     LOGICAL, INTENT(OUT)  ::  found 
!
!
!     found = .TRUE. 
!        
!        
!     SELECT CASE ( restart_string(1:length) )
!
!       CASE ( 'param1' )
!          READ ( 13 )  param1
!
!        CASE DEFAULT
!
!          found = .FALSE.   
!
!     END SELECT
! 
!  END SUBROUTINE uvem_rrd_global  
    

! !------------------------------------------------------------------------------!
! ! Description:
! ! ------------
! !> This routine writes the global restart data.
! !------------------------------------------------------------------------------!
!  SUBROUTINE uvem_wrd_global
! 
!
!     IMPLICIT NONE
!
!
!     CALL wrd_write_string( 'param1' )
!     WRITE ( 14 )  param1         
! 
!        
!        
!  END SUBROUTINE uvem_wrd_global   


 END MODULE uv_exposure_model_mod
