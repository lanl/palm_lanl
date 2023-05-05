 MODULE ice_surface_mod

    IMPLICIT NONE

    PUBLIC initialize_ice_surface

CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize horizontal surface elements, upward- and downward-facing. 
!> Note, horizontal surface type also comprises model-top fluxes, which are,
!> initialized in a different routine. 
!------------------------------------------------------------------------------!
   SUBROUTINE initialize_ice_surface( surf, num_h )
   
      USE control_parameters, ONLY: ice_cover
   
      USE netcdf_data_input_mod,                       &
          ONLY :  input_pids_static,                   &
                  top_surface_fraction_f
   
      IMPLICIT NONE 
   
      TYPE( surf_type ) :: surf          !< respective surface type
      INTEGER(iwp)  ::  num_h            !< current number of surface element
   
      IF ( TRIM(ice_cover) == 'full' ) THEN
         surf%ice_fraction(num_h) = 1.0_wp
      ELSEIF ( TRIM(ice_cover) == 'read_from_file' ) THEN
         IF ( input_pids_static .AND. top_surface_fraction_f%from_file ) THEN
            !TODO assign from netcdf
            surf%ice_fraction(num_h) = &
               top_surface_fraction_f%frac(0,surf%j(num_h),surf%i(num_h))
         ELSE
            WRITE( message_string,* ) 'top_surface_fraction not properly initialized'
            CALL message( 'surface_mod', 'PA0999', 1, 2, 0, 6, 0 )
         ENDIF
      ENDIF
   
   ENDSUBROUTINE initialize_ice_surface
 
 END MODULE ice_surface_mod
