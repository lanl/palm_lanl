MODULE ice_surface_mod

    IMPLICIT NONE

    SAVE

    PRIVATE

!
!-- Public subroutines and functions
    PUBLIC ice_init

CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize horizontal surface elements, upward- and downward-facing. 
!> Note, horizontal surface type also comprises model-top fluxes, which are,
!> initialized in a different routine. 
!------------------------------------------------------------------------------!
SUBROUTINE ice_init

   USE control_parameters, ONLY: ice_cover, message_string

   USE kinds

   USE netcdf_data_input_mod,                       &
       ONLY :  input_pids_static, top_surface_fraction_f

   USE surface_mod, ONLY: surf_def_h

   IMPLICIT NONE 

   INTEGER(iwp)  ::  i,j,m !< current number of surface element

   IF ( TRIM(ice_cover) == 'full' ) THEN
      write(message_string,*) 'Full ice cover: Assigning ice fraction to 1 everywhere'
      call location_message(message_string, .true.)
      surf_def_h(2)%ice_fraction(:) = 1.0_wp
   ELSEIF ( TRIM(ice_cover) == 'read_from_file' ) THEN
      IF ( input_pids_static  .AND.  top_surface_fraction_f%from_file )  THEN
         write(message_string,*) 'Partial ice cover: Assigning ice fraction from file'
         call location_message(message_string, .true.)
         DO  m = 1, surf_def_h(2)%ns
            i = surf_def_h(2)%i(m)
            j = surf_def_h(2)%j(m)
            surf_def_h(2)%ice_fraction(m) = top_surface_fraction_f%frac(0,j,i)
         ENDDO
      ELSEIF (.not. input_pids_static) then
         write(message_string,*) 'Partial ice cover: input_pids_static not found'
         call location_message(message_string, .true.)
      ELSE
         write(message_string,*) 'Partial ice cover: top_surface_fraction not from_file'
         call location_message(message_string, .true.)
      ENDIF
   ENDIF

END SUBROUTINE ice_init
 
END MODULE ice_surface_mod
