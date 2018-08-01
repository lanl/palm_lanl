!> @file chem_photolysis_mod.f90
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
! Copyright 2017-2018 Karlsruhe Institute of Technology
!------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: chem_photolysis_mod.f90 2766 2018-01-22 17:17:47Z kanani $
! Removed preprocessor directive __chem
! 
! 2718 2018-01-02 08:49:38Z maronga
! Initial revision
!
! 
!
!
! Authors:
! --------
! @author Renate Forkel
!
!
!------------------------------------------------------------------------------!
! Description:
! ------------
!> photolysis models and interfaces (Adapted from photolysis_model_mod.f90)
!> @todo Alles!
!------------------------------------------------------------------------------!
 MODULE chem_photolysis_mod

    USE arrays_3d,                                                             &
        ONLY:  dzw, hyp, pt, q, ql, zu, zw

    USE constants,                                                             &
        ONLY:  pi

    USE control_parameters,                                                    &
        ONLY:  time_since_reference_point

    USE pegrid,             ONLY: myid, threads_per_task

    USE indices,                                                               &
        ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt

    USE radiation_model_mod,                                                   &
        ONLY:  calc_zenith, zenith
               !day, day_init, lat, lambda, lon,                   &
               !time_utc, time_utc_init,    !FKa: This is now handled by date_and_time_mod

    USE control_parameters,                                                    &
        ONLY:  initializing_actions           
!       ONLY:  cloud_droplets, cloud_physics, g, initializing_actions,         &
!              large_scale_forcing, lsf_surf, phi, pt_surface, rho_surface,    &
!              surface_pressure, time_since_reference_point

    USE chem_gasphase_mod,                                                     &
        ONLY: nphot, phot_names, phot 

    USE chemistry_model_mod,                                                   &
        ONLY: phot_frequen, photolysis_scheme

    USE chem_modules,                                                          &
        ONLY: chem_debug2

!   REAL(wp), DIMENSION(0:0) ::  zenith,         & !< cosine of solar zenith angle
!                                sun_dir_lat,    & !< solar directional vector in latitudes
!                                sun_dir_lon       !< solar directional vector in longitudes
    USE kinds

#if defined ( __netcdf )
    USE NETCDF
#endif


    IMPLICIT NONE


!   LOGICAL ::  unscheduled_photolysis_calls = .TRUE., & !< flag parameter indicating whether additional calls of the photolysis code are allowed
!               constant_albedo = .FALSE.,           & !< flag parameter indicating whether the albedo may change depending on zenith
!               force_photolysis_call = .FALSE.,     & !< flag parameter for unscheduled photolysis calls
!               photolysis = .FALSE.,                & !< flag parameter indicating whether the photolysis model is used
!               sun_up    = .TRUE.,                  & !< flag parameter indicating whether the sun is up or down
!               photolysis = .TRUE.,                 & !< flag parameter indicing whether photolysis shall be calculated
!               sun_direction = .FALSE.                !< flag parameter indicing whether solar direction shall be calculated

!-- Parameters for constant photolysis frequencies
     INTEGER,PARAMETER :: nconst  = 15               !< available predefined photolysis prequencies for constant 

! Names for predefined fixed photolysis frequencies at zenith angle 0
     CHARACTER(LEN=10), PARAMETER, DIMENSION(nconst) :: names_c =  (/                    &
                     'J_O31D    ','J_O33P    ','J_NO2     ','J_HNO3    ','J_RCHO    ', &
                     'J         ','J         ','J         ','J         ','J         ', &
                     'J         ','J         ','J         ','J         ','J         ' /)
! Photolysis frequency at zenith angle 0 in 1/s
     REAL(wp), PARAMETER, DIMENSION(nconst) :: phot0 =  (/                             &
                     2.489E-05_wp,3.556E-04_wp, 8.89E-03_wp,5.334E-07_wp,3.734E-05_wp, &
                     0.0000E00_wp,0.0000E00_wp,0.0000E00_wp,0.0000E00_wp,0.0000E00_wp, &
                     0.0000E00_wp,0.0000E00_wp,0.0000E00_wp,0.0000E00_wp,0.0000E00_wp /)


!-- Parameters for simple photolysis frequencies
     INTEGER,PARAMETER :: nsimple = 15               !< available predefined photolysis prequencies for simpel parameterisation
! Names for simple photolysis frequencies parameterisation (
     CHARACTER(LEN=10), PARAMETER, DIMENSION(nsimple) :: names_s =  (/                 &
                     'J_O31D    ','J_O33P    ','J_H2O2    ','J_NO2     ','J_NO3_A   ', &
                     'J_NO3_B   ','J_HONO    ','J_HNO3    ','J_HCHO_A  ','J_HCHO_B  ', &
                     'J_CH3CHO  ','J         ','J         ','J         ','J_RCHO    ' /)

!-- Parameters for simeple photolysis frequencies from MCM (http://mcm.leeds.ac.uk/MCM)
!-- Saunders et al., 2003, Atmos. Chem. Phys., 3, 161-180
     REAL(wp), PARAMETER, DIMENSION(nconst) :: par_l =  (/                             &
                     6.073E-05_wp,4.775E-04_wp,1.041E-05_wp,1.165E-02_wp,2.485E-02_wp, &
                     1.747E-01_wp,2.644E-03_wp,9.312E-07_wp,4.642E-05_wp,6.853E-05_wp, &
                     7.344E-06_wp,0.0000E00_wp,0.0000E00_wp,0.000E00_wp, 6.853E-05_wp /)

     REAL(wp), PARAMETER, DIMENSION(nconst) :: par_m =  (/                             &
                         1.743_wp,    0.298_wp,    0.723_wp,    0.244_wp,    0.168_wp, &
                         0.155_wp,    0.261_wp,    1.230_wp,    0.762_wp,    0.477_wp, &
                         1.202_wp,    0.000_wp,    0.000_wp,    0.000_wp,    0.477_wp /)

     REAL(wp), PARAMETER, DIMENSION(nconst) :: par_n =  (/                             &
                         0.474_wp,    0.080_wp,    0.279_wp,    0.267_wp,    0.108_wp, &
                         0.125_wp,    0.288_wp,    0.307_wp,    0.353_wp,    0.323_wp, &
                         0.417_wp,    0.000_wp,    0.000_wp,    0.000_wp,    0.323_wp /)


    REAL(wp) :: time_photolysis = 0.0_wp,         & !< time since last call of photolysis code
                dt_photolysis = 0.0_wp,           & !< hotolysis model timestep
                skip_time_do_photolysis = 0.0_wp    !< Radiation model is not called before this time

    REAL(wp)     :: cosz = 0.7_wp                   !< cosine of Zenith angle (45 deg, if not specified otherwise)

!
!-- Variables and parameters used in Fast-J only
! ....

!   INTERFACE photolysis_check_parameters
!      MODULE PROCEDURE photolysis_check_parameters
!   END INTERFACE photolysis_check_parameters
  
    INTERFACE photolysis_constant
       MODULE PROCEDURE photolysis_constant
    END INTERFACE photolysis_constant
 
    INTERFACE photolysis_simple
       MODULE PROCEDURE photolysis_simple
    END INTERFACE photolysis_simple
!
!   INTERFACE photolysis_fastj
!      MODULE PROCEDURE photolysis_fastj
!   END INTERFACE photolysis_fastj
!
    INTERFACE photolysis_control
       MODULE PROCEDURE photolysis_control
    END INTERFACE photolysis_control

!   INTERFACE photolysis_header
!      MODULE PROCEDURE photolysis_header
!   END INTERFACE photolysis_header 
! 
    INTERFACE photolysis_init
       MODULE PROCEDURE photolysis_init
    END INTERFACE photolysis_init

!   INTERFACE photolysis_parin
!      MODULE PROCEDURE photolysis_parin
!   END INTERFACE photolysis_parin
    
!   INTERFACE photolysis_read_restart_data
!      MODULE PROCEDURE photolysis_read_restart_data
!   END INTERFACE photolysis_read_restart_data

!   INTERFACE photolysis_last_actions
!      MODULE PROCEDURE photolysis_last_actions
!   END INTERFACE photolysis_last_actions

    SAVE

    PRIVATE

!
!-- Public functions / NEEDS SORTING
!   PUBLIC photolysis_init
!          photolysis_check_parameters, photolysis_control,                      &
!          photolysis_header, photolysis_init, photolysis_parin !,                 &
!          photolysis_define_netcdf_grid, photolysis_last_actions,               &
!          photolysis_read_restart_data, photolysis_data_output_mask

   PUBLIC  photolysis_control, photolysis_init

   PUBLIC  photolysis_scheme
!

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine controls the calls of the photolysis schemes
!------------------------------------------------------------------------------!
    SUBROUTINE photolysis_control
 
 
       IMPLICIT NONE


       SELECT CASE ( TRIM( photolysis_scheme ) )

          CASE ( 'constant' )
             CALL photolysis_constant
          
          CASE ( 'simple' )
             CALL photolysis_simple
       
!         CASE ( 'fastj' )
!            CALL photolysis_fastj

          CASE DEFAULT

       END SELECT


    END SUBROUTINE photolysis_control

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the photolysis model
!------------------------------------------------------------------------------!
    SUBROUTINE photolysis_init    !### not yes used. Call should be placed near call of chem_init
    
       IMPLICIT NONE

!-- Just in case we need anything

       RETURN

    END SUBROUTINE photolysis_init


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This scheme keeps the prescribed net radiation constant during the run
!> Default zenith angle is 45 deg
!------------------------------------------------------------------------------!
    SUBROUTINE photolysis_constant


       IMPLICIT NONE

!      INTEGER(iwp) :: i, j, k   !< loop indices
       INTEGER(iwp) :: iphot,iav !< loop indix for photolysis reaction
!      REAL(wp)     :: exn,   &  !< Exner functions at surface
!                      pt1       !< potential temperature at first grid level

       DO iphot = 1,nphot
          DO iav = 1,nconst
             IF ( TRIM( names_c(iav) ) == TRIM( phot_names(iphot) ) ) then
!--    Prescribe fixed photolysis frequencies  [1/s]
!             IF(myid == 0 .AND. chem_debug2 )  WRITE(06,*) iphot, iav,phot_names(iphot),names_c(iav)

                      phot_frequen(iphot)%freq(nzb+1:nzt,:,:) =    &
                            phot0(iav) * cosz

!              IF(myid == 0 .AND. chem_debug2 )  WRITE(06,*) phot_frequen(iphot)%freq(1,5,5)
            ENDIF
          ENDDO
       ENDDO

    END SUBROUTINE photolysis_constant


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This scheme applies a simple parameterisation for clear sky photolysis frequencies 
!> from the Master Chemical Mechanism, MCM v3.2 (http://mcm.leeds.ac.uk/MCM).
!> Reference: Saunders et al., Atmos. Chem. Phys., 3, 161, 2003
!------------------------------------------------------------------------------!
    SUBROUTINE photolysis_simple

       IMPLICIT NONE

!      INTEGER(iwp) :: i, j, k   !< loop indices
       INTEGER(iwp) :: iphot,iav !< loop indix for photolysis reaction
!      REAL(wp)     :: exn,   &  !< Exner functions at surface
!                      pt1       !< potential temperature at first grid level
       REAL(wp)     :: coszi     !< 1./cosine of zenith angle

    DO iphot = 1,nphot
       phot_frequen(iphot)%freq = 0.0_wp
    ENDDO

    CALL calc_zenith

    IF ( zenith(0) > 0.0_wp ) THEN
       coszi = 1. / zenith(0)

       DO iphot = 1,nphot
          DO iav = 1,nsimple
             IF ( TRIM( names_s(iav) ) == TRIM( phot_names(iphot) ) ) then
!              if(myid == 0 .AND. chem_debug2 )  WRITE(06,*) 'simple', iphot, iav,phot_names(iphot),names_s(iav)

                      phot_frequen(iphot)%freq(nzb+1:nzt,:,:) =    &
                            par_l(iav) * zenith(0)**par_m(iav) *  EXP( -par_n(iav) * coszi ) 

!              if(myid == 0 .AND. chem_debug2 )  WRITE(06,*) 'simple', phot_frequen(iphot)%freq(1,5,5)
            ENDIF
          ENDDO
       ENDDO
    ENDIF

    END SUBROUTINE photolysis_simple

 END MODULE chem_photolysis_mod
