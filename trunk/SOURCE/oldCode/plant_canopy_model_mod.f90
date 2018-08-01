!> @file plant_canopy_model_mod.f90
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
! $Id: plant_canopy_model_mod.f90 3065 2018-06-12 07:03:02Z Giersch $
! dz was replaced by the help of zw to allow for vertical stretching
! 
! 3049 2018-05-29 13:52:36Z Giersch
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error message revised
! 
! 3022 2018-05-18 11:12:35Z suehring
! Bugfix in allocation of transpiration rate
! 
! 3014 2018-05-09 08:42:38Z maronga
! Bugfix: nzb_do and nzt_do were not used for 3d data output
! Added pc_transpiration_rate
! 
! 2977 2018-04-17 10:27:57Z kanani
! Implement changes from branch radiation (r2948-2971) with minor modifications,
! plus some formatting.
! (moh.hefny):
! Add plant canopy type to account for changes in LAD (based on the changes
! done by Resler & Pavel) and correct the error message to PALM Standard.
! 
! 2932 2018-03-26 09:39:22Z maronga
! renamed canopy_par to plant_canopy_parameters
! 
! 2920 2018-03-22 11:22:01Z kanani
! Move usm_lad_rma and prototype_lad to radiation_model (moh.hefny)
! 
! 2892 2018-03-14 15:06:29Z suehring
! Bugfix, read separate ASCII LAD files for parent and child model. 
! 
! 2770 2018-01-25 15:10:09Z kanani
! Correction of parameter check
! 
! 2768 2018-01-24 15:38:29Z kanani
! Added check for output quantity pcm_heatrate, some formatting
! 
! 2766 2018-01-22 17:17:47Z kanani
! Increased LEN of canopy mode to 30
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
! Change in file header (GPL part)
! Bugfix for vertical loop index pch_index in case of Netcdf input
! Introduce 2D index array incorporate canopy top index
! Check if canopy on top of topography do not exceed vertical dimension 
! Add check for canopy_mode in case of Netcdf input.
! Enable _FillValue output for 3d quantities
! Bugfix in reading of PIDS leaf area density (MS)
! 
! 2669 2017-12-06 16:03:27Z raasch
! coupling_char removed
! 
! 2512 2017-10-04 08:26:59Z raasch
! upper bounds of 3d output changed from nx+1,ny+1 to nx,ny
! no output of ghost layer data
! 
! 2318 2017-07-20 17:27:44Z suehring
! Get topography top index via Function call 
! 
! 2317 2017-07-20 17:27:19Z suehring
! Changed error messages
! 
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography concept
! 
! 2213 2017-04-24 15:10:35Z kanani
! Bugfix: exchange of ghost points in array pc_heating_rate needed for output 
! of pcm_heatrate, onetime ghost point exchange of lad_s after initialization.
! Formatting and clean-up of subroutine pcm_read_plant_canopy_3d,
! minor re-organization of canopy-heating initialization.
! 
! 2209 2017-04-19 09:34:46Z kanani
! Added 3d output of leaf area density (pcm_lad) and canopy 
! heat rate (pcm_heatrate)
! 
! 2024 2016-10-12 16:42:37Z kanani
! Added missing lad_s initialization
! 
! 2011 2016-09-19 17:29:57Z kanani
! Renamed canopy_heat_flux to pc_heating_rate, since the original meaning/
! calculation of the quantity has changed, related to the urban surface model
! and similar future applications. 
! 
! 2007 2016-08-24 15:47:17Z kanani
! Added SUBROUTINE pcm_read_plant_canopy_3d for reading 3d plant canopy data
! from file (new case canopy_mode=read_from_file_3d) in the course of
! introduction of urban surface model,
! introduced variable ext_coef,
! resorted SUBROUTINEs to alphabetical order
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1960 2016-07-12 16:34:24Z suehring
! Separate humidity and passive scalar
! 
! 1953 2016-06-21 09:28:42Z suehring
! Bugfix, lad_s and lad must be public
! 
! 1826 2016-04-07 12:01:39Z maronga
! Further modularization
!
! 1721 2015-11-16 12:56:48Z raasch
! bugfixes: shf is reduced in areas covered with canopy only,
!           canopy is set on top of topography
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1484 2014-10-21 10:53:05Z kanani
! Changes due to new module structure of the plant canopy model:
!   module plant_canopy_model_mod now contains a subroutine for the 
!   initialization of the canopy model (pcm_init),
!   limitation of the canopy drag (previously accounted for by calculation of 
!   a limiting canopy timestep for the determination of the maximum LES timestep
!   in subroutine timestep) is now realized by the calculation of pre-tendencies
!   and preliminary velocities in subroutine pcm_tendency,
!   some redundant MPI communication removed in subroutine pcm_init
!   (was previously in init_3d_model),
!   unnecessary 3d-arrays lad_u, lad_v, lad_w removed - lad information on the
!   respective grid is now provided only by lad_s (e.g. in the calculation of 
!   the tendency terms or of cum_lai_hf),
!   drag_coefficient, lai, leaf_surface_concentration, 
!   scalar_exchange_coefficient, sec and sls renamed to canopy_drag_coeff, 
!   cum_lai_hf, leaf_surface_conc, leaf_scalar_exch_coeff, lsec and lsc,
!   respectively,
!   unnecessary 3d-arrays cdc, lsc and lsec now defined as single-value constants,
!   USE-statements and ONLY-lists modified accordingly
! 
! 1340 2014-03-25 19:45:13Z kanani
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
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 138 2007-11-28 10:03:58Z letzel
! Initial revision
!
! Description:
! ------------
!> 1) Initialization of the canopy model, e.g. construction of leaf area density 
!> profile (subroutine pcm_init).
!> 2) Calculation of sinks and sources of momentum, heat and scalar concentration 
!> due to canopy elements (subroutine pcm_tendency).
!------------------------------------------------------------------------------!
 MODULE plant_canopy_model_mod
 
    USE arrays_3d,                                                             &
        ONLY:  dzu, dzw, e, q, s, tend, u, v, w, zu, zw 

    USE indices,                                                               &
        ONLY:  nbgp, nxl, nxlg, nxlu, nxr, nxrg, nyn, nyng, nys, nysg, nysv,   &
               nz, nzb, nzt

    USE kinds

    USE surface_mod,                                                           &
        ONLY:  get_topography_top_index_ji


    IMPLICIT NONE


    CHARACTER (LEN=30)   ::  canopy_mode = 'block' !< canopy coverage

    INTEGER(iwp) ::  pch_index = 0                               !< plant canopy height/top index
    INTEGER(iwp) ::  lad_vertical_gradient_level_ind(10) = -9999 !< lad-profile levels (index)

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  pch_index_ji   !< local plant canopy top

    LOGICAL ::  calc_beta_lad_profile = .FALSE. !< switch for calc. of lad from beta func.

    REAL(wp) ::  alpha_lad = 9999999.9_wp         !< coefficient for lad calculation
    REAL(wp) ::  beta_lad = 9999999.9_wp          !< coefficient for lad calculation
    REAL(wp) ::  canopy_drag_coeff = 0.0_wp       !< canopy drag coefficient (parameter)
    REAL(wp) ::  cdc = 0.0_wp                     !< canopy drag coeff. (abbreviation used in equations)
    REAL(wp) ::  cthf = 0.0_wp                    !< canopy top heat flux
    REAL(wp) ::  dt_plant_canopy = 0.0_wp         !< timestep account. for canopy drag
    REAL(wp) ::  ext_coef = 0.6_wp                !< extinction coefficient
    REAL(wp) ::  lad_surface = 0.0_wp             !< lad surface value
    REAL(wp) ::  lai_beta = 0.0_wp                !< leaf area index (lai) for lad calc.
    REAL(wp) ::  leaf_scalar_exch_coeff = 0.0_wp  !< canopy scalar exchange coeff.
    REAL(wp) ::  leaf_surface_conc = 0.0_wp       !< leaf surface concentration
    REAL(wp) ::  lsc = 0.0_wp                     !< leaf surface concentration
    REAL(wp) ::  lsec = 0.0_wp                    !< leaf scalar exchange coeff.

    REAL(wp) ::  lad_vertical_gradient(10) = 0.0_wp              !< lad gradient
    REAL(wp) ::  lad_vertical_gradient_level(10) = -9999999.9_wp !< lad-prof. levels (in m)

    REAL(wp) ::  lad_type_coef(0:10) = 1.0_wp   !< multiplicative coeficients for particular types
                                                !< of plant canopy (e.g. deciduous tree during winter)

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  lad            !< leaf area density
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pre_lad        !< preliminary lad
    
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  cum_lai_hf       !< cumulative lai for heatflux calc.
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  lad_s            !< lad on scalar-grid
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pc_heating_rate  !< plant canopy heating rate
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pc_transpiration_rate  !< plant canopy transpiration rate

    SAVE


    PRIVATE
  
!
!-- Public functions
    PUBLIC pcm_check_data_output, pcm_check_parameters, pcm_data_output_3d,    &
           pcm_define_netcdf_grid, pcm_header, pcm_init, pcm_parin, pcm_tendency

!
!-- Public variables and constants
    PUBLIC pc_heating_rate, pc_transpiration_rate, canopy_mode, cthf, dt_plant_canopy, lad, lad_s,   &
           pch_index
           

    INTERFACE pcm_check_data_output
       MODULE PROCEDURE pcm_check_data_output
    END INTERFACE pcm_check_data_output
    
    INTERFACE pcm_check_parameters
       MODULE PROCEDURE pcm_check_parameters
    END INTERFACE pcm_check_parameters

    INTERFACE pcm_data_output_3d
       MODULE PROCEDURE pcm_data_output_3d
    END INTERFACE pcm_data_output_3d

    INTERFACE pcm_define_netcdf_grid
       MODULE PROCEDURE pcm_define_netcdf_grid
    END INTERFACE pcm_define_netcdf_grid
    
     INTERFACE pcm_header
       MODULE PROCEDURE pcm_header
    END INTERFACE pcm_header        
    
    INTERFACE pcm_init
       MODULE PROCEDURE pcm_init
    END INTERFACE pcm_init

    INTERFACE pcm_parin
       MODULE PROCEDURE pcm_parin
    END INTERFACE pcm_parin

    INTERFACE pcm_read_plant_canopy_3d
       MODULE PROCEDURE pcm_read_plant_canopy_3d
    END INTERFACE pcm_read_plant_canopy_3d
    
    INTERFACE pcm_tendency
       MODULE PROCEDURE pcm_tendency
       MODULE PROCEDURE pcm_tendency_ij
    END INTERFACE pcm_tendency


 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output for plant canopy model
!------------------------------------------------------------------------------!
 SUBROUTINE pcm_check_data_output( var, unit )
 
 
    USE control_parameters,                                                 &
        ONLY:  data_output, message_string, urban_surface

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  unit  !< 
    CHARACTER (LEN=*) ::  var   !<


    SELECT CASE ( TRIM( var ) )

       CASE ( 'pcm_heatrate' )
          IF ( cthf == 0.0_wp  .AND. .NOT.  urban_surface )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                              'res setting of parameter cthf /= 0.0'
             CALL message( 'pcm_check_data_output', 'PA1000', 1, 2, 0, 6, 0 )
          ENDIF
          unit = 'K s-1'
    
       CASE ( 'pcm_transpirationrate' )
          unit = 'kg kg-1 s-1'

       CASE ( 'pcm_lad' )
          unit = 'm2 m-3'


       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE pcm_check_data_output
 
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for plant canopy model
!------------------------------------------------------------------------------!
    SUBROUTINE pcm_check_parameters

       USE control_parameters,                                                 &
           ONLY: cloud_physics, coupling_char, message_string,                 &
                 microphysics_seifert

       USE netcdf_data_input_mod,                                              &
           ONLY:  input_file_static, input_pids_static
                 
    
       IMPLICIT NONE

    
       IF ( canopy_drag_coeff == 0.0_wp )  THEN
          message_string = 'plant_canopy = .TRUE. requires a non-zero drag '// &
                           'coefficient & given value is canopy_drag_coeff = 0.0'
          CALL message( 'pcm_check_parameters', 'PA0041', 1, 2, 0, 6, 0 )
       ENDIF
    
       IF ( ( alpha_lad /= 9999999.9_wp  .AND.  beta_lad == 9999999.9_wp ) .OR.&
              beta_lad /= 9999999.9_wp   .AND.  alpha_lad == 9999999.9_wp )  THEN
          message_string = 'using the beta function for the construction ' //  &
                           'of the leaf area density profile requires '    //  &
                           'both alpha_lad and beta_lad to be /= 9999999.9'
          CALL message( 'pcm_check_parameters', 'PA0118', 1, 2, 0, 6, 0 )
       ENDIF
    
       IF ( calc_beta_lad_profile  .AND.  lai_beta == 0.0_wp )  THEN
          message_string = 'using the beta function for the construction ' //  &
                           'of the leaf area density profile requires '    //  &
                           'a non-zero lai_beta, but given value is '      //  &
                           'lai_beta = 0.0'
          CALL message( 'pcm_check_parameters', 'PA0119', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( calc_beta_lad_profile  .AND.  lad_surface /= 0.0_wp )  THEN
          message_string = 'simultaneous setting of alpha_lad /= 9999999.9 '// &
                           'combined with beta_lad /= 9999999.9 '           // &
                           'and lad_surface /= 0.0 is not possible, '       // &
                           'use either vertical gradients or the beta '     // &
                           'function for the construction of the leaf area '// &
                           'density profile'
          CALL message( 'pcm_check_parameters', 'PA0120', 1, 2, 0, 6, 0 )
       ENDIF 

       IF ( cloud_physics  .AND.  microphysics_seifert )  THEN
          message_string = 'plant_canopy = .TRUE. requires cloud_scheme /=' // &
                          ' seifert_beheng'
          CALL message( 'pcm_check_parameters', 'PA0360', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    If dynamic input file is used, canopy need to be read from file
       IF ( input_pids_static  .AND.                                           &
            TRIM( canopy_mode ) /= 'read_from_file_3d' )  THEN
          message_string = 'Usage of dynamic input file ' //                   &
                           TRIM( input_file_static ) //                        &
                           TRIM( coupling_char ) // ' requires ' //            &
                           'canopy_mode = read_from_file_3d'
          CALL message( 'pcm_check_parameters', 'PA0999', 1, 2, 0, 6, 0 )
       ENDIF

 
    END SUBROUTINE pcm_check_parameters 
 

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 3D output variables
!------------------------------------------------------------------------------!
 SUBROUTINE pcm_data_output_3d( av, variable, found, local_pf, fill_value,     &
                                nzb_do, nzt_do )
  
    USE indices

    USE kinds


    IMPLICIT NONE

    CHARACTER (LEN=*) ::  variable !< 

    INTEGER(iwp) ::  av     !< 
    INTEGER(iwp) ::  i      !< 
    INTEGER(iwp) ::  j      !< 
    INTEGER(iwp) ::  k      !< 
    INTEGER(iwp) ::  k_topo !< topography top index
    INTEGER(iwp) ::  nzb_do !< lower limit of the data output (usually 0)
    INTEGER(iwp) ::  nzt_do !< vertical upper limit of the data output (usually nz_do3d)

    LOGICAL      ::  found !< 

    REAL(wp)     ::  fill_value
    REAL(sp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf !< 


    found = .TRUE.

    local_pf = REAL( fill_value, KIND = 4 )

    SELECT CASE ( TRIM( variable ) )

      CASE ( 'pcm_heatrate' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  IF ( pch_index_ji(j,i) /= 0 )  THEN
                     k_topo = get_topography_top_index_ji( j, i, 's' )
                     DO  k = k_topo, k_topo + pch_index_ji(j,i)
                        local_pf(i,j,k) = pc_heating_rate(k-k_topo,j,i)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
   
       CASE ( 'pcm_transpirationrate' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  IF ( pch_index_ji(j,i) /= 0 )  THEN
                     k_topo = get_topography_top_index_ji( j, i, 's' )
                     DO  k = k_topo, k_topo + pch_index_ji(j,i)
                        local_pf(i,j,k) = pc_transpiration_rate(k-k_topo,j,i)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
    
      CASE ( 'pcm_lad' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  IF ( pch_index_ji(j,i) /= 0 )  THEN
                     k_topo = get_topography_top_index_ji( j, i, 's' )
                     DO  k = k_topo, k_topo + pch_index_ji(j,i)
                        local_pf(i,j,k) = lad_s(k-k_topo,j,i)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
                  
         
       CASE DEFAULT
          found = .FALSE.

    END SELECT


 END SUBROUTINE pcm_data_output_3d 
         
!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!> It is called from subroutine netcdf.
!------------------------------------------------------------------------------!
 SUBROUTINE pcm_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )
    
     IMPLICIT NONE

     CHARACTER (LEN=*), INTENT(IN)  ::  var         !< 
     LOGICAL, INTENT(OUT)           ::  found       !< 
     CHARACTER (LEN=*), INTENT(OUT) ::  grid_x      !< 
     CHARACTER (LEN=*), INTENT(OUT) ::  grid_y      !< 
     CHARACTER (LEN=*), INTENT(OUT) ::  grid_z      !< 

     found  = .TRUE.

!
!--  Check for the grid
     SELECT CASE ( TRIM( var ) )

        CASE ( 'pcm_heatrate', 'pcm_lad', 'pcm_transpirationrate')
           grid_x = 'x'
           grid_y = 'y'
           grid_z = 'zu'

        CASE DEFAULT
           found  = .FALSE.
           grid_x = 'none'
           grid_y = 'none'
           grid_z = 'none'
     END SELECT

 END SUBROUTINE pcm_define_netcdf_grid
 
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for plant canopy model
!------------------------------------------------------------------------------!
    SUBROUTINE pcm_header ( io )

       USE control_parameters,                                                 &
           ONLY: passive_scalar


       IMPLICIT NONE
 
       CHARACTER (LEN=10) ::  coor_chr            !<

       CHARACTER (LEN=86) ::  coordinates         !<
       CHARACTER (LEN=86) ::  gradients           !<
       CHARACTER (LEN=86) ::  leaf_area_density   !<
       CHARACTER (LEN=86) ::  slices              !<
 
       INTEGER(iwp) :: i                !< 
       INTEGER(iwp),  INTENT(IN) ::  io !< Unit of the output file
       INTEGER(iwp) :: k                !<       
   
       REAL(wp) ::  canopy_height       !< canopy height (in m)
       
       canopy_height = zw(pch_index)

       WRITE ( io, 1 )  canopy_mode, canopy_height, pch_index,                 &
                          canopy_drag_coeff
       IF ( passive_scalar )  THEN
          WRITE ( io, 2 )  leaf_scalar_exch_coeff,                             &
                             leaf_surface_conc
       ENDIF

!
!--    Heat flux at the top of vegetation
       WRITE ( io, 3 )  cthf

!
!--    Leaf area density profile, calculated either from given vertical 
!--    gradients or from beta probability density function.
       IF (  .NOT.  calc_beta_lad_profile )  THEN

!--       Building output strings, starting with surface value
          WRITE ( leaf_area_density, '(F7.4)' )  lad_surface
          gradients = '------'
          slices = '     0'
          coordinates = '   0.0'
          i = 1
          DO  WHILE ( i < 11  .AND.  lad_vertical_gradient_level_ind(i)        &
                      /= -9999 )

             WRITE (coor_chr,'(F7.2)')  lad(lad_vertical_gradient_level_ind(i))
             leaf_area_density = TRIM( leaf_area_density ) // ' ' //           &
                                 TRIM( coor_chr )
 
             WRITE (coor_chr,'(F7.2)')  lad_vertical_gradient(i)
             gradients = TRIM( gradients ) // ' ' // TRIM( coor_chr )

             WRITE (coor_chr,'(I7)')  lad_vertical_gradient_level_ind(i)
             slices = TRIM( slices ) // ' ' // TRIM( coor_chr )

             WRITE (coor_chr,'(F7.1)')  lad_vertical_gradient_level(i)
             coordinates = TRIM( coordinates ) // ' '  // TRIM( coor_chr )

             i = i + 1
          ENDDO

          WRITE ( io, 4 )  TRIM( coordinates ), TRIM( leaf_area_density ),     &
                             TRIM( gradients ), TRIM( slices )

       ELSE
       
          WRITE ( leaf_area_density, '(F7.4)' )  lad_surface
          coordinates = '   0.0'
          
          DO  k = 1, pch_index

             WRITE (coor_chr,'(F7.2)')  lad(k)
             leaf_area_density = TRIM( leaf_area_density ) // ' ' //           &
                                 TRIM( coor_chr )
 
             WRITE (coor_chr,'(F7.1)')  zu(k)
             coordinates = TRIM( coordinates ) // ' '  // TRIM( coor_chr )

          ENDDO       

          WRITE ( io, 5 ) TRIM( coordinates ), TRIM( leaf_area_density ),      &
                          alpha_lad, beta_lad, lai_beta

       ENDIF  

1 FORMAT (//' Vegetation canopy (drag) model:'/                                &
              ' ------------------------------'//                              &
              ' Canopy mode: ', A /                                            &
              ' Canopy height: ',F6.2,'m (',I4,' grid points)' /               &
              ' Leaf drag coefficient: ',F6.2 /)
2 FORMAT (/ ' Scalar exchange coefficient: ',F6.2 /                            &
              ' Scalar concentration at leaf surfaces in kg/m**3: ',F6.2 /)
3 FORMAT (' Predefined constant heatflux at the top of the vegetation: ',F6.2, &
          ' K m/s')
4 FORMAT (/ ' Characteristic levels of the leaf area density:'//               &
              ' Height:              ',A,'  m'/                                &
              ' Leaf area density:   ',A,'  m**2/m**3'/                        &
              ' Gradient:            ',A,'  m**2/m**4'/                        &
              ' Gridpoint:           ',A)
5 FORMAT (//' Characteristic levels of the leaf area density and coefficients:'&
          //  ' Height:              ',A,'  m'/                                &
              ' Leaf area density:   ',A,'  m**2/m**3'/                        &
              ' Coefficient alpha: ',F6.2 /                                    &
              ' Coefficient beta: ',F6.2 /                                     &
              ' Leaf area index: ',F6.2,'  m**2/m**2' /)   
       
    END SUBROUTINE pcm_header
 
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the plant canopy model
!------------------------------------------------------------------------------!
    SUBROUTINE pcm_init
    

       USE control_parameters,                                                 &
           ONLY: humidity, io_blocks, io_group, message_string, ocean,         &
                 passive_scalar, urban_surface

       USE netcdf_data_input_mod,                                              &
           ONLY:  leaf_area_density_f

       USE surface_mod,                                                        &
           ONLY: surf_def_h, surf_lsm_h, surf_usm_h

       IMPLICIT NONE

       CHARACTER(10) :: pct
       
       INTEGER(iwp) ::  i   !< running index
       INTEGER(iwp) ::  ii  !< index       
       INTEGER(iwp) ::  j   !< running index
       INTEGER(iwp) ::  k   !< running index
       INTEGER(iwp) ::  m   !< running index

       REAL(wp) ::  int_bpdf        !< vertical integral for lad-profile construction
       REAL(wp) ::  dzh             !< vertical grid spacing in units of canopy height
       REAL(wp) ::  gradient        !< gradient for lad-profile construction
       REAL(wp) ::  canopy_height   !< canopy height for lad-profile construction
       REAL(wp) ::  pcv(nzb:nzt+1)  !<
       
!
!--    Allocate one-dimensional arrays for the computation of the 
!--    leaf area density (lad) profile
       ALLOCATE( lad(0:nz+1), pre_lad(0:nz+1) )
       lad = 0.0_wp
       pre_lad = 0.0_wp

!
!--    Set flag that indicates that the lad-profile shall be calculated by using
!--    a beta probability density function
       IF ( alpha_lad /= 9999999.9_wp  .AND.  beta_lad /= 9999999.9_wp )  THEN
          calc_beta_lad_profile = .TRUE.
       ENDIF
       
       
!
!--    Compute the profile of leaf area density used in the plant
!--    canopy model. The profile can either be constructed from
!--    prescribed vertical gradients of the leaf area density or by
!--    using a beta probability density function (see e.g. Markkanen et al., 
!--    2003: Boundary-Layer Meteorology, 106, 437-459) 
       IF (  .NOT.  calc_beta_lad_profile )  THEN   

!
!--       Use vertical gradients for lad-profile construction    
          i = 1
          gradient = 0.0_wp

          IF (  .NOT.  ocean )  THEN

             lad(0) = lad_surface
             lad_vertical_gradient_level_ind(1) = 0
 
             DO k = 1, pch_index
                IF ( i < 11 )  THEN
                   IF ( lad_vertical_gradient_level(i) < zu(k)  .AND.          &
                        lad_vertical_gradient_level(i) >= 0.0_wp )  THEN
                      gradient = lad_vertical_gradient(i)
                      lad_vertical_gradient_level_ind(i) = k - 1
                      i = i + 1
                   ENDIF
                ENDIF
                IF ( gradient /= 0.0_wp )  THEN
                   IF ( k /= 1 )  THEN
                      lad(k) = lad(k-1) + dzu(k) * gradient
                   ELSE
                      lad(k) = lad_surface + dzu(k) * gradient
                   ENDIF
                ELSE
                   lad(k) = lad(k-1)
                ENDIF
             ENDDO

          ENDIF

!
!--       In case of no given leaf area density gradients, choose a vanishing
!--       gradient. This information is used for the HEADER and the RUN_CONTROL
!--       file.
          IF ( lad_vertical_gradient_level(1) == -9999999.9_wp )  THEN
             lad_vertical_gradient_level(1) = 0.0_wp
          ENDIF

       ELSE

! 
!--       Use beta function for lad-profile construction
          int_bpdf = 0.0_wp
          canopy_height = zw(pch_index)

          DO k = 0, pch_index
             int_bpdf = int_bpdf +                                             &
                      ( ( ( zw(k) / canopy_height )**( alpha_lad-1.0_wp ) ) *  &
                      ( ( 1.0_wp - ( zw(k) / canopy_height ) )**(              &
                          beta_lad-1.0_wp ) )                                  &
                      * ( ( zw(k+1)-zw(k) ) / canopy_height ) )
          ENDDO

!
!--       Preliminary lad profile (defined on w-grid)
          DO k = 0, pch_index
             pre_lad(k) =  lai_beta *                                          &
                        ( ( ( zw(k) / canopy_height )**( alpha_lad-1.0_wp ) )  &
                        * ( ( 1.0_wp - ( zw(k) / canopy_height ) )**(          &
                              beta_lad-1.0_wp ) ) / int_bpdf                   &
                        ) / canopy_height
          ENDDO

!
!--       Final lad profile (defined on scalar-grid level, since most prognostic
!--       quantities are defined there, hence, less interpolation is required 
!--       when calculating the canopy tendencies)
          lad(0) = pre_lad(0)
          DO k = 1, pch_index
             lad(k) = 0.5 * ( pre_lad(k-1) + pre_lad(k) )
          ENDDO          

       ENDIF

!
!--    Allocate 3D-array for the leaf area density (lad_s).
       ALLOCATE( lad_s(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

!
!--    Initialize canopy parameters cdc (canopy drag coefficient), 
!--    lsec (leaf scalar exchange coefficient), lsc (leaf surface concentration)
!--    with the prescribed values
       cdc = canopy_drag_coeff
       lsec = leaf_scalar_exch_coeff
       lsc = leaf_surface_conc

!
!--    Initialization of the canopy coverage in the model domain:
!--    Setting the parameter canopy_mode = 'block' initializes a canopy, which
!--    fully covers the domain surface
       SELECT CASE ( TRIM( canopy_mode ) )

          CASE( 'block' )

             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   lad_s(:,j,i) = lad(:)
                ENDDO
             ENDDO

          CASE ( 'read_from_file_3d' )
!
!--          Initialize LAD with data from file. If LAD is given in NetCDF file,
!--          use these values, else take LAD profiles from ASCII file. 
!--          Please note, in NetCDF file LAD is only given up to the maximum 
!--          canopy top, indicated by leaf_area_density_f%nz.  
             lad_s = 0.0_wp
             IF ( leaf_area_density_f%from_file )  THEN
!
!--             Set also pch_index, used to be the upper bound of the vertical 
!--             loops. Therefore, use the global top of the canopy layer. 
                pch_index = leaf_area_density_f%nz - 1

                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = 0, leaf_area_density_f%nz - 1
                      IF ( leaf_area_density_f%var(k,j,i) /=                   &
                           leaf_area_density_f%fill )                          &
                         lad_s(k,j,i) = leaf_area_density_f%var(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
                CALL exchange_horiz( lad_s, nbgp )
!
!            ASCII file
!--          Initialize canopy parameters cdc (canopy drag coefficient),
!--          lsec (leaf scalar exchange coefficient), lsc (leaf surface concentration)
!--          from file which contains complete 3D data (separate vertical profiles for
!--          each location).
             ELSE
                CALL pcm_read_plant_canopy_3d
             ENDIF

          CASE DEFAULT
!
!--          The DEFAULT case is reached either if the parameter 
!--          canopy mode contains a wrong character string or if the 
!--          user has coded a special case in the user interface. 
!--          There, the subroutine user_init_plant_canopy checks 
!--          which of these two conditions applies.
             CALL user_init_plant_canopy
 
       END SELECT
!
!--    Initialize 2D index array indicating canopy top index. 
       ALLOCATE( pch_index_ji(nysg:nyng,nxlg:nxrg) )
       pch_index_ji = 0

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = 0, pch_index
                IF ( lad_s(k,j,i) /= 0 )  pch_index_ji(j,i) = k
             ENDDO
!
!--          Check whether topography and local vegetation on top exceed 
!--          height of the model domain.
             k = get_topography_top_index_ji( j, i, 's' )
             IF ( k + pch_index_ji(j,i) >= nzt + 1 )  THEN
                message_string =  'Local vegetation height on top of ' //      &
                                  'topography exceeds height of model domain.'
                CALL message( 'pcm_init', 'PA0999', 2, 2, 0, 6, 0 )
             ENDIF

          ENDDO
       ENDDO

       CALL exchange_horiz_2d_int( pch_index_ji, nys, nyn, nxl, nxr, nbgp )

!
!--    Initialization of the canopy heat source distribution due to heating
!--    of the canopy layers by incoming solar radiation, in case that a non-zero
!--    value is set for the canopy top heat flux (cthf), which equals the
!--    available net radiation at canopy top.
!--    The heat source distribution is calculated by a decaying exponential 
!--    function of the downward cumulative leaf area index (cum_lai_hf), 
!--    assuming that the foliage inside the plant canopy is heated by solar
!--    radiation penetrating the canopy layers according to the distribution
!--    of net radiation as suggested by Brown & Covey (1966; Agric. Meteorol. 3, 
!--    73–96). This approach has been applied e.g. by Shaw & Schumann (1992;
!--    Bound.-Layer Meteorol. 61, 47–64).
!--    When using the urban surface model (USM), canopy heating (pc_heating_rate)
!--    by radiation is calculated in the USM.
       IF ( cthf /= 0.0_wp  .AND. .NOT.  urban_surface )  THEN

          ALLOCATE( cum_lai_hf(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                 &
                    pc_heating_rate(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!
!--       Piecewise calculation of the cumulative leaf area index by vertical
!--       integration of the leaf area density
          cum_lai_hf(:,:,:) = 0.0_wp
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                DO  k = pch_index_ji(j,i)-1, 0, -1
                   IF ( k == pch_index_ji(j,i)-1 )  THEN
                      cum_lai_hf(k,j,i) = cum_lai_hf(k+1,j,i) +                &
                         ( 0.5_wp * lad_s(k+1,j,i) *                           &
                           ( zw(k+1) - zu(k+1) ) )  +                          &
                         ( 0.5_wp * ( 0.5_wp * ( lad_s(k+1,j,i) +              &
                                                 lad_s(k,j,i) ) +              &
                                      lad_s(k+1,j,i) ) *                       &
                           ( zu(k+1) - zw(k) ) ) 
                   ELSE
                      cum_lai_hf(k,j,i) = cum_lai_hf(k+1,j,i) +                &
                         ( 0.5_wp * ( 0.5_wp * ( lad_s(k+2,j,i) +              &
                                                 lad_s(k+1,j,i) ) +            &
                                      lad_s(k+1,j,i) ) *                       &
                           ( zw(k+1) - zu(k+1) ) )  +                          &
                         ( 0.5_wp * ( 0.5_wp * ( lad_s(k+1,j,i) +              &
                                                 lad_s(k,j,i) ) +              &
                                      lad_s(k+1,j,i) ) *                       &
                           ( zu(k+1) - zw(k) ) )
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

!            
!--       In areas with canopy the surface value of the canopy heat 
!--       flux distribution overrides the surface heat flux (shf) 
!--       Start with default surface type
          DO  m = 1, surf_def_h(0)%ns
             k = surf_def_h(0)%k(m)
             IF ( cum_lai_hf(0,j,i) /= 0.0_wp )                                &
                surf_def_h(0)%shf(m) = cthf * exp( -ext_coef * cum_lai_hf(0,j,i) )
          ENDDO
!
!--       Natural surfaces
          DO  m = 1, surf_lsm_h%ns
             k = surf_lsm_h%k(m)
             IF ( cum_lai_hf(0,j,i) /= 0.0_wp )                                &
                surf_lsm_h%shf(m) = cthf * exp( -ext_coef * cum_lai_hf(0,j,i) )
          ENDDO
!
!--       Urban surfaces
          DO  m = 1, surf_usm_h%ns
             k = surf_usm_h%k(m)
             IF ( cum_lai_hf(0,j,i) /= 0.0_wp )                                &
                surf_usm_h%shf(m) = cthf * exp( -ext_coef * cum_lai_hf(0,j,i) )
          ENDDO
!
!
!--       Calculation of the heating rate (K/s) within the different layers of 
!--       the plant canopy. Calculation is only necessary in areas covered with
!--       canopy. 
!--       Within the different canopy layers the plant-canopy heating
!--       rate (pc_heating_rate) is calculated as the vertical 
!--       divergence of the canopy heat fluxes at the top and bottom
!--       of the respective layer
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                DO  k = 1, pch_index_ji(j,i)
                   IF ( cum_lai_hf(0,j,i) /= 0.0_wp )  THEN
                      pc_heating_rate(k,j,i) = cthf *                          &
                                ( exp(-ext_coef*cum_lai_hf(k,j,i)) -           &
                                  exp(-ext_coef*cum_lai_hf(k-1,j,i) ) ) / dzw(k)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

       ENDIF
!
!--    Allocate transpiration rate
       IF ( humidity )                                                         &
          ALLOCATE( pc_transpiration_rate(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )



    END SUBROUTINE pcm_init


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &plant_canopy_parameters for plant canopy model
!------------------------------------------------------------------------------!
    SUBROUTINE pcm_parin

       USE control_parameters,                                                 &
           ONLY:  message_string, plant_canopy

       IMPLICIT NONE

       CHARACTER (LEN=80) ::  line  !< dummy string that contains the current line of the parameter file 
       
       NAMELIST /plant_canopy_parameters/                                      &
                                  alpha_lad, beta_lad, canopy_drag_coeff,      &
                                  canopy_mode, cthf,                           &
                                  lad_surface, lad_type_coef,                  & 
                                  lad_vertical_gradient,                       &
                                  lad_vertical_gradient_level,                 &
                                  lai_beta,                                    &
                                  leaf_scalar_exch_coeff,                      &
                                  leaf_surface_conc, pch_index

       NAMELIST /canopy_par/      alpha_lad, beta_lad, canopy_drag_coeff,      &
                                  canopy_mode, cthf,                           &
                                  lad_surface, lad_type_coef,                  & 
                                  lad_vertical_gradient,                       &
                                  lad_vertical_gradient_level,                 &
                                  lai_beta,                                    &
                                  leaf_scalar_exch_coeff,                      &
                                  leaf_surface_conc, pch_index
                                  
       line = ' '
       
!
!--    Try to find radiation model package
       REWIND ( 11 )
       line = ' '
       DO   WHILE ( INDEX( line, '&plant_canopy_parameters' ) == 0 )
          READ ( 11, '(A)', END=10 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read user-defined namelist
       READ ( 11, plant_canopy_parameters )

!
!--    Set flag that indicates that the radiation model is switched on
       plant_canopy = .TRUE.
       
       GOTO 12
!
!--    Try to find old namelist
 10    REWIND ( 11 )
       line = ' '
       DO   WHILE ( INDEX( line, '&canopy_par' ) == 0 )
          READ ( 11, '(A)', END=12 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read user-defined namelist
       READ ( 11, canopy_par )

       message_string = 'namelist canopy_par is deprecated and will be ' // &
                     'removed in near future. Please use namelist ' //      &
                     'plant_canopy_parameters instead'
       CALL message( 'pcm_parin', 'PA0487', 0, 1, 0, 6, 0 )
       
!
!--    Set flag that indicates that the radiation model is switched on
       plant_canopy = .TRUE.

 12    CONTINUE
       

    END SUBROUTINE pcm_parin



!------------------------------------------------------------------------------!
! Description:
! ------------
!
!> Loads 3D plant canopy data from file. File format is as follows:
!>
!> num_levels
!> dtype,x,y,pctype,value(nzb),value(nzb+1), ... ,value(nzb+num_levels-1)
!> dtype,x,y,pctype,value(nzb),value(nzb+1), ... ,value(nzb+num_levels-1)
!> dtype,x,y,pctype,value(nzb),value(nzb+1), ... ,value(nzb+num_levels-1)
!> ...
!>
!> i.e. first line determines number of levels and further lines represent plant
!> canopy data, one line per column and variable. In each data line,
!> dtype represents variable to be set:
!>
!> dtype=1: leaf area density (lad_s)
!> dtype=2....n: some additional plant canopy input data quantity
!>
!> Zeros are added automatically above num_levels until top of domain.  Any
!> non-specified (x,y) columns have zero values as default.
!------------------------------------------------------------------------------!
    SUBROUTINE pcm_read_plant_canopy_3d
    
       USE control_parameters,                                                 &
           ONLY:  coupling_char, message_string, passive_scalar

       USE indices,                                                            &
           ONLY:  nbgp
           
       IMPLICIT NONE

       INTEGER(iwp)                        ::  dtype     !< type of input data (1=lad)
       INTEGER(iwp)                        ::  pctype    !< type of plant canopy (deciduous,non-deciduous,...)
       INTEGER(iwp)                        ::  i, j      !< running index
       INTEGER(iwp)                        ::  nzp       !< number of vertical layers of plant canopy
       INTEGER(iwp)                        ::  nzpltop   !<
       INTEGER(iwp)                        ::  nzpl      !<
       
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  col   !< vertical column of input data

!
!--    Initialize lad_s array
       lad_s = 0.0_wp
       
!
!--    Open and read plant canopy input data
       OPEN(152, FILE='PLANT_CANOPY_DATA_3D' // TRIM( coupling_char ),         &
                 ACCESS='SEQUENTIAL', ACTION='READ', STATUS='OLD',             &
                 FORM='FORMATTED', ERR=515)
       READ(152, *, ERR=516, END=517)  nzp   !< read first line = number of vertical layers
       
       ALLOCATE( col(0:nzp-1) )

       DO
          READ(152, *, ERR=516, END=517) dtype, i, j, pctype, col(:)
          IF ( i < nxlg  .OR.  i > nxrg  .OR.  j < nysg  .OR.  j > nyng )  CYCLE
          
          SELECT CASE (dtype)
             CASE( 1 )   !< leaf area density
!
!--             This is just the pure canopy layer assumed to be grounded to
!--             a flat domain surface. At locations where plant canopy sits 
!--             on top of any kind of topography, the vertical plant column 
!--             must be "lifted", which is done in SUBROUTINE pcm_tendency.           
                IF ( pctype < 0  .OR.  pctype > 10 )  THEN   !< incorrect plant canopy type
                   WRITE( message_string, * ) 'Incorrect type of plant canopy. '   //  &
                                              'Allowed values 0 <= pctype <= 10, ' //  &
                                              'but pctype is ', pctype
                   CALL message( 'pcm_read_plant_canopy_3d', 'PA0349', 1, 2, 0, 6, 0 )
                ENDIF
                lad_s(0:nzp-1,j,i) = col(0:nzp-1) * lad_type_coef(pctype)
                
             CASE DEFAULT
                WRITE(message_string, '(a,i2,a)')   &
                     'Unknown record type in file PLANT_CANOPY_DATA_3D: "', dtype, '"'
                CALL message( 'pcm_read_plant_canopy_3d', 'PA0530', 1, 2, 0, 6, 0 )
          END SELECT
       ENDDO

515    message_string = 'error opening file PLANT_CANOPY_DATA_3D'
       CALL message( 'pcm_read_plant_canopy_3d', 'PA0531', 1, 2, 0, 6, 0 )

516    message_string = 'error reading file PLANT_CANOPY_DATA_3D'
       CALL message( 'pcm_read_plant_canopy_3d', 'PA0532', 1, 2, 0, 6, 0 )

517    CLOSE(152)
       DEALLOCATE( col )
       
       CALL exchange_horiz( lad_s, nbgp )
        
    END SUBROUTINE pcm_read_plant_canopy_3d
    
    

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of the tendency terms, accounting for the effect of the plant
!> canopy on momentum and scalar quantities.
!>
!> The canopy is located where the leaf area density lad_s(k,j,i) > 0.0 
!> (defined on scalar grid), as initialized in subroutine pcm_init. 
!> The lad on the w-grid is vertically interpolated from the surrounding
!> lad_s. The upper boundary of the canopy is defined on the w-grid at 
!> k = pch_index. Here, the lad is zero.
!>
!> The canopy drag must be limited (previously accounted for by calculation of 
!> a limiting canopy timestep for the determination of the maximum LES timestep
!> in subroutine timestep), since it is physically impossible that the canopy
!> drag alone can locally change the sign of a velocity component. This 
!> limitation is realized by calculating preliminary tendencies and velocities.
!> It is subsequently checked if the preliminary new velocity has a different
!> sign than the current velocity. If so, the tendency is limited in a way that
!> the velocity can at maximum be reduced to zero by the canopy drag. 
!>
!>
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE pcm_tendency( component )


       USE control_parameters,                                                 &
           ONLY:  dt_3d, message_string

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  component !< prognostic variable (u,v,w,pt,q,e)
       INTEGER(iwp) ::  i         !< running index
       INTEGER(iwp) ::  j         !< running index
       INTEGER(iwp) ::  k         !< running index
       INTEGER(iwp) ::  k_wall    !< vertical index of topography top
       INTEGER(iwp) ::  kk        !< running index for flat lad arrays

       REAL(wp) ::  ddt_3d    !< inverse of the LES timestep (dt_3d)
       REAL(wp) ::  lad_local !< local lad value
       REAL(wp) ::  pre_tend  !< preliminary tendency
       REAL(wp) ::  pre_u     !< preliminary u-value 
       REAL(wp) ::  pre_v     !< preliminary v-value 
       REAL(wp) ::  pre_w     !< preliminary w-value 


       ddt_3d = 1.0_wp / dt_3d
 
!
!--    Compute drag for the three velocity components and the SGS-TKE:
       SELECT CASE ( component )

!
!--       u-component
          CASE ( 1 )
             DO  i = nxlu, nxr
                DO  j = nys, nyn
!
!--                Determine topography-top index on u-grid
                   k_wall = get_topography_top_index_ji( j, i, 'u' )
                   DO  k = k_wall+1, k_wall + pch_index_ji(j,i)

                      kk = k - k_wall   !- lad arrays are defined flat
!
!--                   In order to create sharp boundaries of the plant canopy, 
!--                   the lad on the u-grid at index (k,j,i) is equal to 
!--                   lad_s(k,j,i), rather than being interpolated from the 
!--                   surrounding lad_s, because this would yield smaller lad 
!--                   at the canopy boundaries than inside of the canopy.
!--                   For the same reason, the lad at the rightmost(i+1)canopy 
!--                   boundary on the u-grid equals lad_s(k,j,i).
                      lad_local = lad_s(kk,j,i)
                      IF ( lad_local == 0.0_wp .AND. lad_s(kk,j,i-1) > 0.0_wp )&
                      THEN
                         lad_local = lad_s(kk,j,i-1)
                      ENDIF

                      pre_tend = 0.0_wp
                      pre_u = 0.0_wp
!
!--                   Calculate preliminary value (pre_tend) of the tendency
                      pre_tend = - cdc *                                       &
                                   lad_local *                                 &
                                   SQRT( u(k,j,i)**2 +                         &
                                         ( 0.25_wp * ( v(k,j,i-1) +            &
                                                       v(k,j,i)   +            &
                                                       v(k,j+1,i) +            &
                                                       v(k,j+1,i-1) )          &
                                         )**2 +                                &
                                         ( 0.25_wp * ( w(k-1,j,i-1) +          &
                                                       w(k-1,j,i)   +          &
                                                       w(k,j,i-1)   +          &
                                                       w(k,j,i) )              &
                                         )**2                                  &
                                       ) *                                     &
                                   u(k,j,i)

!
!--                   Calculate preliminary new velocity, based on pre_tend
                      pre_u = u(k,j,i) + dt_3d * pre_tend
!
!--                   Compare sign of old velocity and new preliminary velocity,
!--                   and in case the signs are different, limit the tendency
                      IF ( SIGN(pre_u,u(k,j,i)) /= pre_u )  THEN
                         pre_tend = - u(k,j,i) * ddt_3d
                      ELSE
                         pre_tend = pre_tend
                      ENDIF
!
!--                   Calculate final tendency
                      tend(k,j,i) = tend(k,j,i) + pre_tend

                   ENDDO
                ENDDO
             ENDDO

!
!--       v-component
          CASE ( 2 )
             DO  i = nxl, nxr
                DO  j = nysv, nyn
!
!--                Determine topography-top index on v-grid
                   k_wall = get_topography_top_index_ji( j, i, 'v' )

                   DO  k = k_wall+1, k_wall + pch_index_ji(j,i)

                      kk = k - k_wall   !- lad arrays are defined flat
!
!--                   In order to create sharp boundaries of the plant canopy, 
!--                   the lad on the v-grid at index (k,j,i) is equal to 
!--                   lad_s(k,j,i), rather than being interpolated from the 
!--                   surrounding lad_s, because this would yield smaller lad 
!--                   at the canopy boundaries than inside of the canopy.
!--                   For the same reason, the lad at the northmost(j+1) canopy 
!--                   boundary on the v-grid equals lad_s(k,j,i).
                      lad_local = lad_s(kk,j,i)
                      IF ( lad_local == 0.0_wp .AND. lad_s(kk,j-1,i) > 0.0_wp )&
                      THEN
                         lad_local = lad_s(kk,j-1,i)
                      ENDIF

                      pre_tend = 0.0_wp
                      pre_v = 0.0_wp
!
!--                   Calculate preliminary value (pre_tend) of the tendency
                      pre_tend = - cdc *                                       &
                                   lad_local *                                 &
                                   SQRT( ( 0.25_wp * ( u(k,j-1,i)   +          &
                                                       u(k,j-1,i+1) +          &
                                                       u(k,j,i)     +          &
                                                       u(k,j,i+1) )            &
                                         )**2 +                                &
                                         v(k,j,i)**2 +                         &
                                         ( 0.25_wp * ( w(k-1,j-1,i) +          &
                                                       w(k-1,j,i)   +          &
                                                       w(k,j-1,i)   +          &
                                                       w(k,j,i) )              &
                                         )**2                                  &
                                       ) *                                     &
                                   v(k,j,i)

!
!--                   Calculate preliminary new velocity, based on pre_tend
                      pre_v = v(k,j,i) + dt_3d * pre_tend
!
!--                   Compare sign of old velocity and new preliminary velocity,
!--                   and in case the signs are different, limit the tendency
                      IF ( SIGN(pre_v,v(k,j,i)) /= pre_v )  THEN
                         pre_tend = - v(k,j,i) * ddt_3d
                      ELSE
                         pre_tend = pre_tend
                      ENDIF
!
!--                   Calculate final tendency
                      tend(k,j,i) = tend(k,j,i) + pre_tend

                   ENDDO
                ENDDO
             ENDDO

!
!--       w-component
          CASE ( 3 )
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Determine topography-top index on w-grid
                   k_wall = get_topography_top_index_ji( j, i, 'w' )

                   DO  k = k_wall+1, k_wall + pch_index_ji(j,i) - 1

                      kk = k - k_wall   !- lad arrays are defined flat

                      pre_tend = 0.0_wp
                      pre_w = 0.0_wp
!
!--                   Calculate preliminary value (pre_tend) of the tendency
                      pre_tend = - cdc *                                       &
                                   (0.5_wp *                                   &
                                      ( lad_s(kk+1,j,i) + lad_s(kk,j,i) )) *   &
                                   SQRT( ( 0.25_wp * ( u(k,j,i)   +            &
                                                       u(k,j,i+1) +            &
                                                       u(k+1,j,i) +            &
                                                       u(k+1,j,i+1) )          &
                                         )**2 +                                &
                                         ( 0.25_wp * ( v(k,j,i)   +            &
                                                       v(k,j+1,i) +            &
                                                       v(k+1,j,i) +            &
                                                       v(k+1,j+1,i) )          &
                                         )**2 +                                &
                                         w(k,j,i)**2                           &
                                       ) *                                     &
                                   w(k,j,i)
!
!--                   Calculate preliminary new velocity, based on pre_tend
                      pre_w = w(k,j,i) + dt_3d * pre_tend
!
!--                   Compare sign of old velocity and new preliminary velocity,
!--                   and in case the signs are different, limit the tendency
                      IF ( SIGN(pre_w,w(k,j,i)) /= pre_w )  THEN
                         pre_tend = - w(k,j,i) * ddt_3d
                      ELSE
                         pre_tend = pre_tend
                      ENDIF
!
!--                   Calculate final tendency
                      tend(k,j,i) = tend(k,j,i) + pre_tend

                   ENDDO
                ENDDO
             ENDDO

!
!--       potential temperature
          CASE ( 4 )
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Determine topography-top index on scalar-grid
                   k_wall = get_topography_top_index_ji( j, i, 's' )

                   DO  k = k_wall+1, k_wall + pch_index_ji(j,i)

                      kk = k - k_wall   !- lad arrays are defined flat
                      tend(k,j,i) = tend(k,j,i) + pc_heating_rate(kk,j,i)
                   ENDDO
                ENDDO
             ENDDO

!
!--       humidity
          CASE ( 5 )
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Determine topography-top index on scalar-grid
                   k_wall = get_topography_top_index_ji( j, i, 's' )

                   DO  k = k_wall+1, k_wall + pch_index_ji(j,i)

                      kk = k - k_wall   !- lad arrays are defined flat
                      pc_transpiration_rate(kk,j,i) =  - lsec                  &
                                       * lad_s(kk,j,i) *                       &
                                       SQRT( ( 0.5_wp * ( u(k,j,i) +           &
                                                          u(k,j,i+1) )         &
                                             )**2 +                            &
                                             ( 0.5_wp * ( v(k,j,i) +           &
                                                          v(k,j+1,i) )         &
                                             )**2 +                            &
                                             ( 0.5_wp * ( w(k-1,j,i) +         & 
                                                          w(k,j,i) )           &
                                             )**2                              &
                                           ) *                                 &
                                       ( q(k,j,i) - lsc )

                      tend(k,j,i) = tend(k,j,i) + pc_transpiration_rate(kk,j,i)
                   ENDDO
                ENDDO
             ENDDO

!
!--       sgs-tke
          CASE ( 6 )
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Determine topography-top index on scalar-grid
                   k_wall = get_topography_top_index_ji( j, i, 's' )

                   DO  k = k_wall+1, k_wall + pch_index_ji(j,i)

                      kk = k - k_wall   !- lad arrays are defined flat
                      tend(k,j,i) = tend(k,j,i) -                              &
                                       2.0_wp * cdc *                          &
                                       lad_s(kk,j,i) *                         &
                                       SQRT( ( 0.5_wp * ( u(k,j,i) +           &
                                                          u(k,j,i+1) )         &
                                             )**2 +                            &
                                             ( 0.5_wp * ( v(k,j,i) +           &
                                                          v(k,j+1,i) )         &
                                             )**2 +                            &
                                             ( 0.5_wp * ( w(k,j,i) +           &
                                                          w(k+1,j,i) )         &
                                             )**2                              &
                                           ) *                                 &
                                       e(k,j,i)
                   ENDDO
                ENDDO
             ENDDO 
!
!--       scalar concentration
          CASE ( 7 )
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Determine topography-top index on scalar-grid
                   k_wall = get_topography_top_index_ji( j, i, 's' )

                   DO  k = k_wall+1, k_wall + pch_index_ji(j,i)

                      kk = k - k_wall   !- lad arrays are defined flat
                      tend(k,j,i) = tend(k,j,i) -                              &
                                       lsec *                                  &
                                       lad_s(kk,j,i) *                         &
                                       SQRT( ( 0.5_wp * ( u(k,j,i) +           &
                                                          u(k,j,i+1) )         &
                                             )**2 +                            &
                                             ( 0.5_wp * ( v(k,j,i) +           &
                                                          v(k,j+1,i) )         &
                                             )**2 +                            &
                                             ( 0.5_wp * ( w(k-1,j,i) +         & 
                                                          w(k,j,i) )           &
                                             )**2                              &
                                           ) *                                 &
                                       ( s(k,j,i) - lsc )
                   ENDDO
                ENDDO
             ENDDO    



          CASE DEFAULT

             WRITE( message_string, * ) 'wrong component: ', component
             CALL message( 'pcm_tendency', 'PA0279', 1, 2, 0, 6, 0 ) 

       END SELECT

    END SUBROUTINE pcm_tendency


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of the tendency terms, accounting for the effect of the plant
!> canopy on momentum and scalar quantities.
!>
!> The canopy is located where the leaf area density lad_s(k,j,i) > 0.0 
!> (defined on scalar grid), as initialized in subroutine pcm_init. 
!> The lad on the w-grid is vertically interpolated from the surrounding
!> lad_s. The upper boundary of the canopy is defined on the w-grid at 
!> k = pch_index. Here, the lad is zero.
!>
!> The canopy drag must be limited (previously accounted for by calculation of 
!> a limiting canopy timestep for the determination of the maximum LES timestep
!> in subroutine timestep), since it is physically impossible that the canopy
!> drag alone can locally change the sign of a velocity component. This 
!> limitation is realized by calculating preliminary tendencies and velocities.
!> It is subsequently checked if the preliminary new velocity has a different
!> sign than the current velocity. If so, the tendency is limited in a way that
!> the velocity can at maximum be reduced to zero by the canopy drag. 
!>
!>
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE pcm_tendency_ij( i, j, component )


       USE control_parameters,                                                 &
           ONLY:  dt_3d, message_string

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  component !< prognostic variable (u,v,w,pt,q,e)
       INTEGER(iwp) ::  i         !< running index
       INTEGER(iwp) ::  j         !< running index
       INTEGER(iwp) ::  k         !< running index
       INTEGER(iwp) ::  k_wall    !< vertical index of topography top
       INTEGER(iwp) ::  kk        !< running index for flat lad arrays

       REAL(wp) ::  ddt_3d    !< inverse of the LES timestep (dt_3d)
       REAL(wp) ::  lad_local !< local lad value
       REAL(wp) ::  pre_tend  !< preliminary tendency
       REAL(wp) ::  pre_u     !< preliminary u-value 
       REAL(wp) ::  pre_v     !< preliminary v-value 
       REAL(wp) ::  pre_w     !< preliminary w-value 


       ddt_3d = 1.0_wp / dt_3d
!
!--    Compute drag for the three velocity components and the SGS-TKE
       SELECT CASE ( component )

!
!--       u-component
          CASE ( 1 )
!
!--          Determine topography-top index on u-grid
             k_wall = get_topography_top_index_ji( j, i, 'u' )
             DO  k = k_wall + 1, k_wall + pch_index_ji(j,i)

                kk = k - k_wall  !- lad arrays are defined flat

!
!--             In order to create sharp boundaries of the plant canopy, 
!--             the lad on the u-grid at index (k,j,i) is equal to lad_s(k,j,i),
!--             rather than being interpolated from the surrounding lad_s, 
!--             because this would yield smaller lad at the canopy boundaries 
!--             than inside of the canopy.
!--             For the same reason, the lad at the rightmost(i+1)canopy 
!--             boundary on the u-grid equals lad_s(k,j,i).
                lad_local = lad_s(kk,j,i)
                IF ( lad_local == 0.0_wp .AND. lad_s(kk,j,i-1) > 0.0_wp )  THEN
                   lad_local = lad_s(kk,j,i-1)
                ENDIF

                pre_tend = 0.0_wp
                pre_u = 0.0_wp
!
!--             Calculate preliminary value (pre_tend) of the tendency
                pre_tend = - cdc *                                             &
                             lad_local *                                       &   
                             SQRT( u(k,j,i)**2 +                               &
                                   ( 0.25_wp * ( v(k,j,i-1)  +                 &
                                                 v(k,j,i)    +                 &
                                                 v(k,j+1,i)  +                 &
                                                 v(k,j+1,i-1) )                &
                                   )**2 +                                      &
                                   ( 0.25_wp * ( w(k-1,j,i-1) +                &
                                                 w(k-1,j,i)   +                &
                                                 w(k,j,i-1)   +                &
                                                 w(k,j,i) )                    &
                                   )**2                                        &
                                 ) *                                           &
                             u(k,j,i)

!
!--             Calculate preliminary new velocity, based on pre_tend
                pre_u = u(k,j,i) + dt_3d * pre_tend
!
!--             Compare sign of old velocity and new preliminary velocity,
!--             and in case the signs are different, limit the tendency
                IF ( SIGN(pre_u,u(k,j,i)) /= pre_u )  THEN
                   pre_tend = - u(k,j,i) * ddt_3d
                ELSE
                   pre_tend = pre_tend
                ENDIF
!
!--             Calculate final tendency
                tend(k,j,i) = tend(k,j,i) + pre_tend
             ENDDO


!
!--       v-component
          CASE ( 2 )
!
!--          Determine topography-top index on v-grid
             k_wall = get_topography_top_index_ji( j, i, 'v' )

             DO  k = k_wall + 1, k_wall + pch_index_ji(j,i)

                kk = k - k_wall  !- lad arrays are defined flat
!
!--             In order to create sharp boundaries of the plant canopy, 
!--             the lad on the v-grid at index (k,j,i) is equal to lad_s(k,j,i),
!--             rather than being interpolated from the surrounding lad_s, 
!--             because this would yield smaller lad at the canopy boundaries 
!--             than inside of the canopy.
!--             For the same reason, the lad at the northmost(j+1)canopy 
!--             boundary on the v-grid equals lad_s(k,j,i).
                lad_local = lad_s(kk,j,i)
                IF ( lad_local == 0.0_wp .AND. lad_s(kk,j-1,i) > 0.0_wp )  THEN
                   lad_local = lad_s(kk,j-1,i)
                ENDIF

                pre_tend = 0.0_wp
                pre_v = 0.0_wp
!
!--             Calculate preliminary value (pre_tend) of the tendency
                pre_tend = - cdc *                                             &
                             lad_local *                                       &
                             SQRT( ( 0.25_wp * ( u(k,j-1,i)   +                &
                                                 u(k,j-1,i+1) +                &
                                                 u(k,j,i)     +                &
                                                 u(k,j,i+1) )                  &
                                   )**2 +                                      &
                                   v(k,j,i)**2 +                               &
                                   ( 0.25_wp * ( w(k-1,j-1,i) +                &
                                                 w(k-1,j,i)   +                &
                                                 w(k,j-1,i)   +                &
                                                 w(k,j,i) )                    &
                                   )**2                                        &
                                 ) *                                           &
                             v(k,j,i)

!
!--             Calculate preliminary new velocity, based on pre_tend
                pre_v = v(k,j,i) + dt_3d * pre_tend
!
!--             Compare sign of old velocity and new preliminary velocity,
!--             and in case the signs are different, limit the tendency
                IF ( SIGN(pre_v,v(k,j,i)) /= pre_v )  THEN
                   pre_tend = - v(k,j,i) * ddt_3d
                ELSE
                   pre_tend = pre_tend
                ENDIF
!
!--             Calculate final tendency
                tend(k,j,i) = tend(k,j,i) + pre_tend
             ENDDO


!
!--       w-component
          CASE ( 3 )
!
!--          Determine topography-top index on w-grid
             k_wall = get_topography_top_index_ji( j, i, 'w' )

             DO  k = k_wall + 1, k_wall + pch_index_ji(j,i) - 1

                kk = k - k_wall  !- lad arrays are defined flat

                pre_tend = 0.0_wp
                pre_w = 0.0_wp
!
!--             Calculate preliminary value (pre_tend) of the tendency
                pre_tend = - cdc *                                             &
                             (0.5_wp *                                         &
                                ( lad_s(kk+1,j,i) + lad_s(kk,j,i) )) *         &
                             SQRT( ( 0.25_wp * ( u(k,j,i)    +                 &  
                                                 u(k,j,i+1)  +                 &
                                                 u(k+1,j,i)  +                 &
                                                 u(k+1,j,i+1) )                &
                                   )**2 +                                      &
                                   ( 0.25_wp * ( v(k,j,i)    +                 &
                                                 v(k,j+1,i)  +                 &
                                                 v(k+1,j,i)  +                 &
                                                 v(k+1,j+1,i) )                &
                                   )**2 +                                      &
                                   w(k,j,i)**2                                 &
                                 ) *                                           &
                             w(k,j,i)
!
!--             Calculate preliminary new velocity, based on pre_tend
                pre_w = w(k,j,i) + dt_3d * pre_tend
!
!--             Compare sign of old velocity and new preliminary velocity,
!--             and in case the signs are different, limit the tendency
                IF ( SIGN(pre_w,w(k,j,i)) /= pre_w )  THEN
                   pre_tend = - w(k,j,i) * ddt_3d
                ELSE
                   pre_tend = pre_tend
                ENDIF
!
!--             Calculate final tendency
                tend(k,j,i) = tend(k,j,i) + pre_tend
             ENDDO

!
!--       potential temperature
          CASE ( 4 )
!
!--          Determine topography-top index on scalar grid
             k_wall = get_topography_top_index_ji( j, i, 's' )

             DO  k = k_wall + 1, k_wall + pch_index_ji(j,i)
                kk = k - k_wall  !- lad arrays are defined flat
                tend(k,j,i) = tend(k,j,i) + pc_heating_rate(kk,j,i)
             ENDDO


!
!--       humidity
          CASE ( 5 )
!
!--          Determine topography-top index on scalar grid
             k_wall = get_topography_top_index_ji( j, i, 's' )

             DO  k = k_wall + 1, k_wall + pch_index_ji(j,i)
                kk = k - k_wall  !- lad arrays are defined flat

                pc_transpiration_rate(kk,j,i) = - lsec                         &
                                 * lad_s(kk,j,i) *                             &
                                 SQRT( ( 0.5_wp * ( u(k,j,i) +                 &
                                                    u(k,j,i+1) )               &
                                       )**2  +                                 &
                                       ( 0.5_wp * ( v(k,j,i) +                 &
                                                    v(k,j+1,i) )               &
                                       )**2 +                                  &
                                       ( 0.5_wp * ( w(k-1,j,i) +               &
                                                    w(k,j,i) )                 &
                                       )**2                                    &
                                     ) *                                       &
                                 ( q(k,j,i) - lsc )

                tend(k,j,i) = tend(k,j,i) + pc_transpiration_rate(kk,j,i)

             ENDDO   

!
!--       sgs-tke
          CASE ( 6 )
!
!--          Determine topography-top index on scalar grid
             k_wall = get_topography_top_index_ji( j, i, 's' )

             DO  k = k_wall + 1, k_wall + pch_index_ji(j,i)

                kk = k - k_wall
                tend(k,j,i) = tend(k,j,i) -                                    &
                                 2.0_wp * cdc *                                &
                                 lad_s(kk,j,i) *                               &
                                 SQRT( ( 0.5_wp * ( u(k,j,i) +                 &
                                                    u(k,j,i+1) )               &
                                       )**2 +                                  &  
                                       ( 0.5_wp * ( v(k,j,i) +                 &
                                                    v(k,j+1,i) )               &
                                       )**2 +                                  &
                                       ( 0.5_wp * ( w(k,j,i) +                 &
                                                    w(k+1,j,i) )               &
                                       )**2                                    &
                                     ) *                                       &
                                 e(k,j,i)
             ENDDO
             
!
!--       scalar concentration 
          CASE ( 7 )
!
!--          Determine topography-top index on scalar grid
             k_wall = get_topography_top_index_ji( j, i, 's' )

             DO  k = k_wall + 1, k_wall + pch_index_ji(j,i)

                kk = k - k_wall
                tend(k,j,i) = tend(k,j,i) -                                    &
                                 lsec *                                        &
                                 lad_s(kk,j,i) *                               &
                                 SQRT( ( 0.5_wp * ( u(k,j,i) +                 &
                                                    u(k,j,i+1) )               &
                                       )**2  +                                 &
                                       ( 0.5_wp * ( v(k,j,i) +                 &
                                                    v(k,j+1,i) )               &
                                       )**2 +                                  &
                                       ( 0.5_wp * ( w(k-1,j,i) +               &
                                                    w(k,j,i) )                 &
                                       )**2                                    &
                                     ) *                                       &
                                 ( s(k,j,i) - lsc )
             ENDDO                

       CASE DEFAULT

          WRITE( message_string, * ) 'wrong component: ', component
          CALL message( 'pcm_tendency', 'PA0279', 1, 2, 0, 6, 0 ) 

       END SELECT

    END SUBROUTINE pcm_tendency_ij



 END MODULE plant_canopy_model_mod
