!> @file chemistry_model_mod.f90
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
! Copyright 2017-2018 Freie Universitaet Berlin
!------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: chemistry_model_mod.f90 3045 2018-05-28 07:55:41Z Giersch $
! error messages revised
! 
! 3014 2018-05-09 08:42:38Z maronga
! Bugfix: nzb_do and nzt_do were not used for 3d data output
! 
! 3004 2018-04-27 12:33:25Z Giersch
! Comment concerning averaged data output added 
! 
! 2932 2018-03-26 09:39:22Z maronga
! renamed chemistry_par to chemistry_parameters
! 
! 2894 2018-03-15 09:17:58Z Giersch
! Calculations of the index range of the subdomain on file which overlaps with
! the current subdomain are already done in read_restart_data_mod,
! chem_last_actions was renamed to chem_wrd_local, chem_read_restart_data was
! renamed to chem_rrd_local, chem_write_var_list was renamed to 
! chem_wrd_global, chem_read_var_list was renamed to chem_rrd_global, 
! chem_skip_var_list has been removed, variable named found has been
! introduced for checking if restart data was found, reading of restart strings
! has been moved completely to read_restart_data_mod, chem_rrd_local is already
! inside the overlap loop programmed in read_restart_data_mod, todo list has 
! bees extended, redundant characters in chem_wrd_local have been removed, 
! the marker *** end chemistry *** is not necessary anymore, strings and their 
! respective lengths are written out and read now in case of restart runs to 
! get rid of prescribed character lengths
!
! 2815 2018-02-19 11:29:57Z suehring
! Bugfix in restart mechanism,
! rename chem_tendency to chem_prognostic_equations,
! implement vector-optimized version of chem_prognostic_equations,
! some clean up (incl. todo list)
! 
! 2773 2018-01-30 14:12:54Z suehring
! Declare variables required for nesting as public
! 
! 2772 2018-01-29 13:10:35Z suehring
! Bugfix in string handling
! 
! 2768 2018-01-24 15:38:29Z kanani
! Shorten lines to maximum length of 132 characters
! 
! 2766 2018-01-22 17:17:47Z kanani
! Removed preprocessor directive __chem
! 
! 2756 2018-01-16 18:11:14Z suehring
! Fill values in 3D output introduced.
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
! @author Farah Kanani-Suehring
! @author Klaus Ketelsen
! @author Basit Khan
!
!
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Chemistry model for PALM-4U
!> @todo Adjust chem_rrd_local to CASE structure of others modules. It is not 
!>       allowed to use the chemistry model in a precursor run and additionally 
!>       not using it in a main run
!> @todo Update/clean-up todo list! (FK)
!> @todo Set proper fill values (/= 0) for chem output arrays! (FK)
!> @todo Add routine chem_check_parameters, add checks for inconsistent or
!>       unallowed parameter settings. 
!>       CALL of chem_check_parameters from check_parameters. (FK)
!> @todo Make routine chem_header available, CALL from header.f90
!>       (see e.g. how it is done in routine lsm_header in 
!>        land_surface_model_mod.f90). chem_header should include all setup
!>        info about chemistry parameter settings. (FK)
!> @todo Implement turbulent inflow of chem spcs in inflow_turbulence.
!> @todo Separate boundary conditions for each chem spcs to be implemented
!> @todo Currently only total concentration are calculated. Resolved, parameterized
!>       and chemistry fluxes although partially and some completely coded but
!>       are not operational/activated in this version. bK.
!> @todo slight differences in passive scalar and chem spcs when chem reactions
!>       turned off. Need to be fixed. bK
!> @todo test nesting for chem spcs, was implemented by suehring (kanani)
!> @todo subroutine set_const_initial_values to be taken out from chemistry_model_mod !bK.
!> @todo chemistry error messages
!> @todo Format this module according to PALM coding standard (see e.g. module
!>       template under http://palm.muk.uni-hannover.de/mosaik/downloads/8 or
!>       D3_coding_standard.pdf under https://palm.muk.uni-hannover.de/trac/downloads/16)
! 
!------------------------------------------------------------------------------!

MODULE chemistry_model_mod

   USE kinds,              ONLY: wp, iwp
   USE indices,            ONLY: nz, nzb,nzt,nysg,nyng,nxlg,nxrg,nys,nyn
   USE pegrid,             ONLY: myid, threads_per_task
   USE control_parameters, ONLY: dt_3d, ws_scheme_sca, initializing_actions,   &
                                 intermediate_timestep_count,                  &
                                 intermediate_timestep_count_max,              &
                                 message_string, omega, tsc,                   &
                                 timestep_scheme, use_prescribed_profile_data  
   USE arrays_3d,          ONLY: hyp, kh, pt, rdf_sc, tend, zu                     
   USE chem_gasphase_mod,  ONLY: nspec, spc_names, nkppctrl, nmaxfixsteps,     &
                                 t_steps, fill_temp, chem_gasphase_integrate,  &
                                 nvar, atol, rtol, nphot, phot_names
   USE cpulog,             ONLY: cpu_log, log_point

   USE chem_modules
 

   IMPLICIT   NONE
   PRIVATE
   SAVE

!- Define chemical variables

   TYPE   species_def
      CHARACTER(LEN=8)                                   :: name
      CHARACTER(LEN=16)                                  :: unit
      REAL(kind=wp),POINTER,DIMENSION(:,:,:)             :: conc
      REAL(kind=wp),POINTER,DIMENSION(:,:,:)             :: conc_av
      REAL(kind=wp),POINTER,DIMENSION(:,:,:)             :: conc_p
      REAL(kind=wp),POINTER,DIMENSION(:,:,:)             :: tconc_m
      REAL(kind=wp),ALLOCATABLE,DIMENSION(:,:)           :: cssws_av           
      REAL(kind=wp),ALLOCATABLE,DIMENSION(:,:)           :: flux_s_cs
      REAL(kind=wp),ALLOCATABLE,DIMENSION(:,:)           :: diss_s_cs
      REAL(kind=wp),ALLOCATABLE,DIMENSION(:,:,:)         :: flux_l_cs
      REAL(kind=wp),ALLOCATABLE,DIMENSION(:,:,:)         :: diss_l_cs
      REAL(kind=wp),ALLOCATABLE,DIMENSION(:)             :: conc_pr_init
   END TYPE species_def


   TYPE   photols_def                                                            
      CHARACTER(LEN=8)                                   :: name
      CHARACTER(LEN=16)                                  :: unit
      REAL(kind=wp),POINTER,DIMENSION(:,:,:)             :: freq
   END TYPE photols_def


   PUBLIC  species_def
   PUBLIC  photols_def


   TYPE(species_def),ALLOCATABLE,DIMENSION(:),TARGET, PUBLIC     :: chem_species
   TYPE(photols_def),ALLOCATABLE,DIMENSION(:),TARGET, PUBLIC     :: phot_frequen  

   REAL(kind=wp),ALLOCATABLE,DIMENSION(:,:,:,:),TARGET   :: spec_conc_1
   REAL(kind=wp),ALLOCATABLE,DIMENSION(:,:,:,:),TARGET   :: spec_conc_2
   REAL(kind=wp),ALLOCATABLE,DIMENSION(:,:,:,:),TARGET   :: spec_conc_3
   REAL(kind=wp),ALLOCATABLE,DIMENSION(:,:,:,:),TARGET   :: spec_conc_av       
   REAL(kind=wp),ALLOCATABLE,DIMENSION(:,:,:,:),TARGET   :: freq_1

   INTEGER,DIMENSION(nkppctrl)                           :: icntrl                            ! Fine tuning kpp
   REAL(kind=wp),DIMENSION(nkppctrl)                     :: rcntrl                            ! Fine tuning kpp

   CHARACTER(10), PUBLIC :: photolysis_scheme
                                         ! 'constant',
                                         ! 'simple' (Simple parameterisation from MCM, Saunders et al., 2003, Atmos. Chem. Phys., 3, 161-180
                                         ! 'fastj'  (Wild et al., 2000, J. Atmos. Chem., 37, 245-282) STILL NOT IMPLEMENTED

   PUBLIC nspec
   PUBLIC nvar       
   PUBLIC spc_names
   PUBLIC spec_conc_2  

!- Interface section
  INTERFACE chem_boundary_conds
      MODULE PROCEDURE chem_boundary_conds
  END INTERFACE chem_boundary_conds

   INTERFACE chem_init
      MODULE PROCEDURE chem_init
   END INTERFACE chem_init

   INTERFACE chem_init_profiles
      MODULE PROCEDURE chem_init_profiles
   END INTERFACE chem_init_profiles

   INTERFACE chem_parin
      MODULE PROCEDURE chem_parin
   END INTERFACE chem_parin

   INTERFACE chem_integrate
      MODULE PROCEDURE chem_integrate_ij
   END INTERFACE chem_integrate

   INTERFACE chem_swap_timelevel
      MODULE PROCEDURE chem_swap_timelevel
   END INTERFACE chem_swap_timelevel

   INTERFACE chem_define_netcdf_grid
      MODULE PROCEDURE chem_define_netcdf_grid
   END INTERFACE chem_define_netcdf_grid

   INTERFACE chem_data_output_3d
      MODULE PROCEDURE chem_data_output_3d
   END INTERFACE chem_data_output_3d

   INTERFACE chem_check_data_output
      MODULE PROCEDURE chem_check_data_output
   END INTERFACE chem_check_data_output

   INTERFACE chem_check_data_output_pr
      MODULE PROCEDURE chem_check_data_output_pr
   END INTERFACE chem_check_data_output_pr

   INTERFACE chem_3d_data_averaging
      MODULE PROCEDURE chem_3d_data_averaging 
   END INTERFACE chem_3d_data_averaging

   INTERFACE chem_wrd_local
      MODULE PROCEDURE chem_wrd_local 
   END INTERFACE chem_wrd_local

   INTERFACE chem_rrd_local
      MODULE PROCEDURE chem_rrd_local
   END INTERFACE chem_rrd_local

   INTERFACE chem_prognostic_equations
      MODULE PROCEDURE chem_prognostic_equations
      MODULE PROCEDURE chem_prognostic_equations_ij
   END INTERFACE chem_prognostic_equations

   INTERFACE chem_header
      MODULE PROCEDURE chem_header
   END INTERFACE chem_header

   INTERFACE chem_emissions
      MODULE PROCEDURE chem_emissions
   END INTERFACE chem_emissions

!    INTERFACE chem_wrd_global
!       MODULE PROCEDURE chem_wrd_global
!    END INTERFACE chem_wrd_global
! 
!    INTERFACE chem_rrd_global
!       MODULE PROCEDURE chem_rrd_global
!    END INTERFACE chem_rrd_global


   PUBLIC chem_3d_data_averaging, chem_boundary_conds, chem_check_data_output, &
          chem_check_data_output_pr, chem_data_output_3d,                      &
          chem_define_netcdf_grid, chem_emissions, chem_header, chem_init,     &
          chem_init_profiles, chem_integrate, chem_wrd_local,                  &
          chem_parin, chem_prognostic_equations,                               &
          chem_rrd_local, chem_swap_timelevel
          

 CONTAINS

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine to initialize and set all boundary conditions for chemical species
!------------------------------------------------------------------------------!

 SUBROUTINE chem_boundary_conds( mode )                                           
                                                                                  
    USE control_parameters,                                                    &  
        ONLY: air_chemistry, outflow_l, outflow_n, outflow_r, outflow_s                 
    USE indices,                                                               &  
        ONLY: nxl, nxr,  nxlg, nxrg, nyng, nysg, nzt                              
                                                                                  
!    USE prognostic_equations_mod,                                             &

    USE arrays_3d,                                                             &     
        ONLY: dzu                                                
    USE surface_mod,                                                           &
        ONLY: bc_h                                                           

    CHARACTER (len=*), INTENT(IN) :: mode
    INTEGER(iwp) ::  i                                                            !< grid index x direction.
    INTEGER(iwp) ::  j                                                            !< grid index y direction.
    INTEGER(iwp) ::  k                                                            !< grid index z direction.
    INTEGER(iwp) ::  kb                                                           !< variable to set respective boundary value, depends on facing. 
    INTEGER(iwp) ::  l                                                            !< running index boundary type, for up- and downward-facing walls.
    INTEGER(iwp) ::  m                                                            !< running index surface elements.
    INTEGER(iwp) ::  lsp                                                          !< running index for chem spcs.
    INTEGER(iwp) ::  lph                                                          !< running index for photolysis frequencies


    SELECT CASE ( TRIM( mode ) )        
       CASE ( 'init' )     
          DO lsp = 1, nspec            
            IF ( surface_csflux(lsp) == 9999999.9_wp )  THEN
                 constant_csflux(lsp) = .FALSE.            
            ENDIF
          ENDDO

          IF ( bc_cs_b == 'dirichlet' )    THEN
             ibc_cs_b = 0 
          ELSEIF ( bc_cs_b == 'neumann' )  THEN
             ibc_cs_b = 1 
          ELSE
!             message_string = 'unknown boundary condition: bc_cs_b ="' // TRIM( bc_cs_b ) // '"'  ! bK commented
             CALL message( 'chem_boundary_conds', 'CM0010', 1, 2, 0, 6, 0 )     !< chemistry_model_mod should have special error numbers --> "CHEM###",
          ENDIF                                                                  
!
!--       Set Integer flags and check for possible erroneous settings for top
!--       boundary condition. bK added *_cs_* here.
          IF ( bc_cs_t == 'dirichlet' )             THEN
             ibc_cs_t = 0 
          ELSEIF ( bc_cs_t == 'neumann' )           THEN
             ibc_cs_t = 1
          ELSEIF ( bc_cs_t == 'initial_gradient' )  THEN
             ibc_cs_t = 2
          ELSEIF ( bc_cs_t == 'nested' )            THEN          
             ibc_cs_t = 3
          ELSE
!            message_string = 'unknown boundary condition: bc_c_t ="' // TRIM( bc_cs_t ) // '"'  
             CALL message( 'check_parameters', 'CM0011', 1, 2, 0, 6, 0 )
          ENDIF

      
       CASE ( 'set_bc_bottomtop' )                   
!--      Bottom boundary condtions for chemical species      
          DO lsp = 1, nspec                                                      
             IF ( ibc_cs_b == 0 ) THEN                   
                DO l = 0, 1 
!--             Set index kb: For upward-facing surfaces (l=0), kb=-1, i.e.
!--             the chem_species(nsp)%conc_p value at the topography top (k-1)
!--             is set; for downward-facing surfaces (l=1), kb=1, i.e. the 
!--             value at the topography bottom (k+1) is set. 

                   kb = MERGE( -1, 1, l == 0 )
                   !$OMP PARALLEL DO PRIVATE( i, j, k )
                   DO m = 1, bc_h(l)%ns
                      i = bc_h(l)%i(m)            
                      j = bc_h(l)%j(m)
                      k = bc_h(l)%k(m)
                      chem_species(lsp)%conc_p(k+kb,j,i) = chem_species(lsp)%conc(k+kb,j,i) 
                   ENDDO                                        
                ENDDO                                       
          
             ELSEIF ( ibc_cs_b == 1 ) THEN
         ! in boundary_conds there is som extra loop over m here for passive tracer
                DO l = 0, 1
                   kb = MERGE( -1, 1, l == 0 )
                   !$OMP PARALLEL DO PRIVATE( i, j, k )                                           
                   DO m = 1, bc_h(l)%ns
                      i = bc_h(l)%i(m)            
                      j = bc_h(l)%j(m)
                      k = bc_h(l)%k(m)
                      chem_species(lsp)%conc_p(k+kb,j,i) = chem_species(lsp)%conc_p(k,j,i)

!< @todo: chem_species(nsp)%conc_p(k+kb,j,i) = chem_species(nsp)%conc(k+kb,j,i),    &
!<  pls loop over nsp=1, NSPEC.
!< @todo: We should also think about the possibility to have &
!< individual boundary conditions for each species? This means, bc_cs_b,       &
!< bc_cs_t, ibc_cs_b, ibc_cs_t would need to be added to TYPE chem_species(nsp)%,   &
!< and the loop over nsp would have to be put outside of this IF-clause.
!< i think its better we keep the same bonundary cond i.e. dirichlet or neumann
!< for all chem spcs. ... !bK

                   ENDDO
                ENDDO
             ENDIF
          ENDDO    ! end lsp loop  

!--    Top boundary conditions for chemical species - Should this not be done for all species?
          IF ( ibc_cs_t == 0 )  THEN                    
             DO lsp = 1, nspec
                chem_species(lsp)%conc_p(nzt+1,:,:) = chem_species(lsp)%conc(nzt+1,:,:)        
             ENDDO
          ELSEIF ( ibc_cs_t == 1 )  THEN
             DO lsp = 1, nspec
                chem_species(lsp)%conc_p(nzt+1,:,:) = chem_species(lsp)%conc_p(nzt,:,:)
             ENDDO
          ELSEIF ( ibc_cs_t == 2 )  THEN
             DO lsp = 1, nspec
                chem_species(lsp)%conc_p(nzt+1,:,:) = chem_species(lsp)%conc_p(nzt,:,:) + bc_cs_t_val(lsp) * dzu(nzt+1)
             ENDDO
                                    !@todo: bc_cs_t_val needs to be calculated,    
                                    !see bc_pt_t_val = ( pt_init(nzt+1) - pt_init(nzt) ) / dzu(nzt+1) 
                                    !(in time_integration). pt_init would be the counterpart to
                                    !chem_species(i)%conc_pr_init (see kchem_driver_FKa1408.f90 of my 
                                    !"Hints: prescribing initial concentration.
          ENDIF
!
       CASE ( 'set_bc_lateral' )                       ! bK commented it
!--          Lateral boundary conditions for chem species at inflow boundary
!--          are automatically set when chem_species concentration is
!--          initialized. The initially set value at the inflow boundary is not
!--          touched during time integration, hence, this boundary value remains
!--          at a constant value, which is the concentration that flows into the
!--          domain.                                                           
!--          Lateral boundary conditions for chem species at outflow boundary 

          IF ( outflow_s )  THEN
             DO lsp = 1, nspec
                chem_species(lsp)%conc_p(:,nys-1,:) = chem_species(lsp)%conc_p(:,nys,:)
             ENDDO
          ELSEIF ( outflow_n )  THEN
             DO lsp = 1, nspec
                chem_species(lsp)%conc_p(:,nyn+1,:) = chem_species(lsp)%conc_p(:,nyn,:)
             ENDDO
          ELSEIF ( outflow_l )  THEN
             DO lsp = 1, nspec
                chem_species(lsp)%conc_p(:,:,nxl-1) = chem_species(lsp)%conc_p(:,:,nxl)
             ENDDO
          ELSEIF ( outflow_r )  THEN
             DO lsp = 1, nspec
                chem_species(lsp)%conc_p(:,:,nxr+1) = chem_species(lsp)%conc_p(:,:,nxr)
             ENDDO
          ENDIF
         
    END SELECT

 END SUBROUTINE chem_boundary_conds
!
!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining initial vertical profiles of chemical species (given by
!> namelist parameters chem_profiles and chem_heights)  --> which should work
!> analogue to parameters u_profile, v_profile and uv_heights)
!------------------------------------------------------------------------------!
 SUBROUTINE chem_init_profiles              !< SUBROUTINE is called from chem_init in case of 
                                            !< TRIM( initializing_actions ) /= 'read_restart_data'
                                            !< We still need to see what has to be done in case of restart run
    USE chem_modules

    IMPLICIT NONE

!-- Local variables
    INTEGER ::  lsp        !< running index for number of species in derived data type species_def
    INTEGER ::  lsp_pr     !< running index for number of species in cs_names, cs_profiles etc
    INTEGER ::  lpr_lev    !< running index for profile level for each chem spcs. 
    INTEGER ::  npr_lev    !< the next available profile lev

!-----------------
!-- To prescribe/initialize a vertically constant  'cs_profile', another parameter
!-- "cs_surface" is introduced. If "cs_profile" and "cs_heights" are prescribed,
!-- their values will override the constant profile given by "cs_surface".

!     IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN 
       lsp_pr = 1
       DO  WHILE ( TRIM( cs_name( lsp_pr ) ) /= 'novalue' )   !'novalue' is the default
          DO  lsp = 1, nspec                                !
!--          Create initial profile (conc_pr_init) for each chemical species
             IF ( TRIM( chem_species(lsp)%name ) == TRIM( cs_name(lsp_pr) ) )  THEN   !
!                IF ( cs_profile(1,1) == 9999999.9_wp ) THEN
!--               Set a vertically constant profile based on the surface conc (cs_surface(lsp_pr)) of each species
                  DO lpr_lev = 0, nzt+1
                     chem_species(lsp)%conc_pr_init(lpr_lev) = cs_surface(lsp_pr)
                  ENDDO

!                ELSE
!                   IF ( cs_heights(lsp,1) /= 0.0_wp )  THEN
!                      message_string = 'cs_heights(1,1) must be 0.0'
!                      CALL message( 'chem_check_parameters', 'CM0012', 1, 2, 0, 6, 0 )
!                   ENDIF
!  
!                   IF ( omega /= 0.0_wp )  THEN
!                      message_string = 'Coriolis force must be switched off (by setting omega=0.0)' //  &
!                                       ' when prescribing the forcing by u_profile and v_profile'
!                      CALL message( 'check_parameters', 'PA0347', 1, 2, 0, 6, 0 )
!                   ENDIF
!  
!                   use_prescribed_profile_data = .TRUE.
! 
!                   npr_lev = 1
! !                  chem_sddpecies(lsp)%conc_pr_init(0) = 0.0_wp
!                   DO  lpr_lev = 1, nz+1
!                      IF ( npr_lev < 100 )  THEN
!                         DO  WHILE ( cs_heights(lsp, npr_lev+1) <= zu(lpr_lev) )
!                            npr_lev = npr_lev + 1
!                            IF ( npr_lev == 100 )  EXIT
!                         ENDDO
!                      ENDIF
!  
!                      IF ( npr_lev < 100  .AND.  cs_heights(lsp, npr_lev + 1) /= 9999999.9_wp )  THEN
!                         chem_species(lsp)%conc_pr_init(lpr_lev) = cs_profile(lsp, npr_lev) +         &
!                                                 ( zu(lpr_lev) - cs_heights(lsp, npr_lev) ) /         &
!                                                 ( cs_heights(lsp, npr_lev + 1) - cs_heights(lsp, npr_lev ) ) * &
!                                                 ( cs_profile(lsp, npr_lev + 1) - cs_profile(lsp, npr_lev ) )
!                      ELSE
!                         chem_species(lsp)%conc_pr_init(lpr_lev) = cs_profile(lsp, npr_lev)
!                      ENDIF
!                   ENDDO
!                ENDIF
!-- If a profile is prescribed explicity using cs_profiles and cs_heights, we then have to fill the 
!-- chem_species(lsp)%conc_pr_init for the specific "lsp" based on the cs_profiles(lsp_pr,:) 
!-- and cs_heights(lsp_pr,:). 
             ENDIF
          ENDDO
          lsp_pr = lsp_pr + 1
       ENDDO
!     ELSE
!       IF (chem_debug1 ) print*,'code to be added for initializing_actions == read_restart_data'   !bK
!     ENDIF

!-- Now, go back to chem_init and use the contents of chem_species(lsp)%conc_pr_init to
!-- initialize the 3D conc arrays, as it is currently taken care of in chem_set_constant_values.
!-- After initializing the 3D arrays, these can be used to set the boundary conditions in the
!-- subroutine kchem_boundary_conds, which should be called from boundary_conds.f90.


 END SUBROUTINE chem_init_profiles
!
!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine initializating chemistry_model_mod
!------------------------------------------------------------------------------!
   SUBROUTINE chem_init


      USE control_parameters,                                                  &
         ONLY: message_string, io_blocks, io_group, turbulent_inflow
       
      USE arrays_3d,                                                           &
          ONLY: mean_inflow_profiles

      USE pegrid

    IMPLICIT   none
!-- local variables
    INTEGER                 :: i,j               !< running index for for horiz numerical grid points
    INTEGER                 :: lsp               !< running index for chem spcs
    INTEGER                 :: lpr_lev           !< running index for chem spcs profile level
!
!-- NOPOINTER version not implemented yet
! #if defined( __nopointer )
!     message_string = 'The chemistry module only runs with POINTER version'
!     CALL message( 'chemistry_model_mod', 'CM0001', 1, 2, 0, 6, 0 )      
! #endif
!
!-- Allocate memory for chemical species
    ALLOCATE( chem_species(nspec) )
    ALLOCATE( spec_conc_1 (nzb:nzt+1,nysg:nyng,nxlg:nxrg,nspec) )
    ALLOCATE( spec_conc_2 (nzb:nzt+1,nysg:nyng,nxlg:nxrg,nspec) )
    ALLOCATE( spec_conc_3 (nzb:nzt+1,nysg:nyng,nxlg:nxrg,nspec) )
    ALLOCATE( spec_conc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nspec) ) 
    ALLOCATE( phot_frequen(nphot) ) 
    ALLOCATE( freq_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nphot) )
    ALLOCATE( bc_cs_t_val(nspec) )
!
!-- Initialize arrays
    spec_conc_1 (:,:,:,:) = 0.0_wp
    spec_conc_2 (:,:,:,:) = 0.0_wp
    spec_conc_3 (:,:,:,:) = 0.0_wp
    spec_conc_av(:,:,:,:) = 0.0_wp


    DO lsp = 1, nspec
       chem_species(lsp)%name    = spc_names(lsp)

       chem_species(lsp)%conc   (nzb:nzt+1,nysg:nyng,nxlg:nxrg)       => spec_conc_1 (:,:,:,lsp)
       chem_species(lsp)%conc_p (nzb:nzt+1,nysg:nyng,nxlg:nxrg)       => spec_conc_2 (:,:,:,lsp)
       chem_species(lsp)%tconc_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg)       => spec_conc_3 (:,:,:,lsp)
       chem_species(lsp)%conc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg)       => spec_conc_av(:,:,:,lsp)     

       ALLOCATE (chem_species(lsp)%cssws_av(nysg:nyng,nxlg:nxrg))                   
       chem_species(lsp)%cssws_av    = 0.0_wp
!
!--    (todo (FK): This needs to be revised. This block must go somewhere else)                                         
!      IF ( ws_scheme_sca )  THEN                                                
          ALLOCATE (chem_species(lsp)%flux_s_cs(nzb+1:nzt,0:threads_per_task-1))   
          ALLOCATE (chem_species(lsp)%diss_s_cs(nzb+1:nzt,0:threads_per_task-1))   
          ALLOCATE (chem_species(lsp)%flux_l_cs(nzb+1:nzt,nys:nyn,0:threads_per_task-1)) 
          ALLOCATE (chem_species(lsp)%diss_l_cs(nzb+1:nzt,nys:nyn,0:threads_per_task-1))   
          chem_species(lsp)%flux_s_cs = 0.0_wp                                     
          chem_species(lsp)%flux_l_cs = 0.0_wp                                     
          chem_species(lsp)%diss_s_cs = 0.0_wp                                      
          chem_species(lsp)%diss_l_cs = 0.0_wp                                     
!      ENDIF         
!
!--   Allocate memory for initial concentration profiles
!--   (concentration values come from namelist)
!--   (todo (FK): Because of this, chem_init is called in palm before
!--               check_parameters, since conc_pr_init is used there.
!--               We have to find another solution since chem_init should
!--               eventually be called from init_3d_model!!)
      ALLOCATE ( chem_species(lsp)%conc_pr_init(0:nz+1) )
      chem_species(lsp)%conc_pr_init(:) = 0.0_wp

    ENDDO

!
!-- Set initial concentration of profiles prescribed by parameters cs_profile
!-- and cs_heights in the namelist &chemistry_parameters
!-- (todo (FK): chem_init_profiles not ready yet, has some bugs)
!     CALL chem_init_profiles
!
!-- Initialize model variables


    IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.            &
         TRIM( initializing_actions ) /= 'cyclic_fill' )  THEN


!--    First model run of a possible job queue.
!--    Initial profiles of the variables must be computes.
       IF ( INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN
!            CALL location_message( 'initializing with 1D model profiles', .FALSE. )
!
!          CALL init_1d_model    ...... decide to call or not later     !bK

!--        Transfer initial profiles to the arrays of the 3D model
          DO lsp = 1, nspec
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO lpr_lev = 1, nz + 1
                      chem_species(lsp)%conc(lpr_lev,j,i) = chem_species(lsp)%conc_pr_init(lpr_lev)
                   ENDDO
                ENDDO
             ENDDO   
          ENDDO
       
       ELSEIF ( INDEX(initializing_actions, 'set_constant_profiles') /= 0 )    &
       THEN
!           CALL location_message( 'initializing with constant profiles', .FALSE. )



!--       Set initial horizontal velocities at the lowest computational grid 
!--       levels to zero in order to avoid too small time steps caused by the 
!--       diffusion limit in the initial phase of a run (at k=1, dz/2 occurs
!--       in the limiting formula!).

          DO lsp = 1, nspec 
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   chem_species(lsp)%conc(:,j,i) = chem_species(lsp)%conc_pr_init    !ITS THERE bK
                ENDDO
             ENDDO
          ENDDO

!        ELSEIF ( INDEX(initializing_actions, 'by_user') /= 0 )                  &
!        THEN
!           CALL location_message( 'initializing by user', .FALSE. )
!
!--       Initialization will completely be done by the user
!--       (FK: This should be called only once, in init_3d_model, i.e. remove it here)
!           CALL user_init_3d_model
!           CALL location_message( 'finished', .TRUE. )

       ENDIF

!-- Store initial chem spcs profile
!         DO lsp = 1, nvar
!          hom_cs(:,1,115,:) = SPREAD(  chem_species(lsp)%conc(:,nys,nxl), 2, statistic_regions+1 )
!         ENDDO
!
!-- If required, change the surface chem spcs at the start of the 3D run
       IF ( cs_surface_initial_change(1) /= 0.0_wp ) THEN            
          DO lsp = 1, nspec 
             chem_species(lsp)%conc(nzb,:,:) = chem_species(lsp)%conc(nzb,:,:) +  &
                                               cs_surface_initial_change(lsp)
          ENDDO
       ENDIF 
!
!-- Initiale old and new time levels. 
       DO lsp = 1, nvar
          chem_species(lsp)%tconc_m = 0.0_wp                      
          chem_species(lsp)%conc_p  = chem_species(lsp)%conc     
       ENDDO

    ENDIF



!--- new code add above this line
    DO lsp = 1, nphot
       phot_frequen(lsp)%name = phot_names(lsp)
!        IF( myid == 0 )  THEN
!           WRITE(6,'(a,i4,3x,a)')  'Photolysis: ',lsp,trim(phot_names(lsp))
!        ENDIF 
         phot_frequen(lsp)%freq(nzb:nzt+1,nysg:nyng,nxlg:nxrg)    => freq_1(:,:,:,lsp)
    ENDDO

!--   Set initial values
!    Not required any more, this can now be done with the namelist  by setting cs_surface 
!    and cs_name without specifying cs_profile (Nevertheless, we still want to keep it for a while) 
!   IF (  TRIM( initializing_actions ) /= 'read_restart_data' )  THEN             
!         CALL set_const_initial_values
!   ENDIF

    RETURN

    CONTAINS
!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine setting constant initial values of chemical species
!------------------------------------------------------------------------------!
     SUBROUTINE set_const_initial_values
!    Not required any more, this can now be done with the namelist  by setting cs_surface 
!    and cs_name without specifying cs_profile (Nevertheless, we still want to keep it for a while) 
         IMPLICIT   none

!--      local variables
         INTEGER                  :: lsp

         IF(myid == 0)  THEN
            write(6,'(/,a,/)')  ' chemics >>>> Set constant Initial Values: '
         ENDIF

!        Default values are taken from smog.def from supplied kpp example
         DO lsp = 1, nspec
            IF(trim(chem_species(lsp)%name) == 'NO')           THEN
!              chem_species(lsp)%conc = 8.725*1.0E+08
!              chem_species(lsp)%conc = 0.05_wp                          !added by bK
               chem_species(lsp)%conc = 0.01_wp                          !added by RFo
            ELSEIF (trim(chem_species(lsp)%name) == 'NO2')     THEN
!              chem_species(lsp)%conc = 2.240*1.0E+08             
!              chem_species(lsp)%conc = 0.01_wp                          !added by bK
               chem_species(lsp)%conc = 0.05_wp                          !added by RFo
            ELSEIF( trim( chem_species(lsp)%name ) == 'O3' )   THEN
               chem_species(lsp)%conc = 0.05_wp                          !added by bK
            ELSEIF(trim(chem_species(lsp)%name) == 'H2O')      THEN
!              chem_species(lsp)%conc = 5.326*1.0E+11
               chem_species(lsp)%conc = 1.30*1.0E+4_wp                   !added by bK
            ELSEIF(trim(chem_species(lsp)%name) == 'O2')       THEN
               chem_species(lsp)%conc = 2.0*1.0E+5_wp                    !added by bK
            ELSEIF(trim(chem_species(lsp)%name) == 'RH')       THEN
               chem_species(lsp)%conc = 0.001_wp                         !added by RFo
            ELSEIF(trim(chem_species(lsp)%name) == 'CO')       THEN
               chem_species(lsp)%conc = 0.5_wp                           !added by RFo
            ELSEIF(trim(chem_species(lsp)%name) == 'RCHO')     THEN
!              chem_species(lsp)%conc = 2.0_wp                           !added by bK
               chem_species(lsp)%conc = 0.01_wp                          !added by RFo
!           ELSEIF(trim(chem_species(lsp)%name) == 'OH')       THEN
!              chem_species(lsp)%conc = 1.0*1.0E-07_wp                   !added by bK
!           ELSEIF(trim(chem_species(lsp)%name) == 'HO2')      THEN
!              chem_species(lsp)%conc = 1*1.0E-7_wp                      !added by bK
!           ELSEIF(trim(chem_species(lsp)%name) == 'RCOO2')    THEN      ! corrected RFo
!              chem_species(lsp)%conc = 1.0*1.0E-7_wp                    !added by bK
!           ELSEIF(trim(chem_species(lsp)%name) == 'RCOO2NO2') THEN
!              chem_species(lsp)%conc = 1.0*1.0E-7_wp                    !added by bK
            ELSE
!              H2O = 2.0e+04;
               chem_species(lsp)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg) = 0.0_wp
            ENDIF

            IF(myid == 0)  THEN
               WRITE(6,'(a,3x,a,3x,a,e12.4)')  '   Species:     ',chem_species(lsp)%name(1:7), & 
                                        'Initial Value = ',chem_species(lsp)%conc(nzb,nysg,nxlg)
            ENDIF
         ENDDO

!   #if defined( __nopointer )
!   !kk      Hier mit message abbrechen
!            if(myid == 0)  then
!               write(6,*)  '   KPP does only run with POINTER Version'
!            end if
!            stop 'error'
!   #endif

         RETURN
      END SUBROUTINE set_const_initial_values


   END SUBROUTINE chem_init
!
!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining parin for &chemistry_parameters for chemistry model
!------------------------------------------------------------------------------!
   SUBROUTINE chem_parin
   
      USE control_parameters,                                                  &
          ONLY: air_chemistry
          
      USE chem_modules
      
      USE kinds

      IMPLICIT   none

      CHARACTER (LEN=80) ::  line   !< dummy string that contains the current line of the parameter file
      
      REAL(wp), DIMENSION(nmaxfixsteps) ::   my_steps   !< List of fixed timesteps   my_step(1) = 0.0 automatic stepping

      NAMELIST /chemistry_parameters/   bc_cs_b,                          &
                                        bc_cs_t,                          &
                                        call_chem_at_all_substeps,        &
                                        chem_debug0,                      &
                                        chem_debug1,                      &
                                        chem_debug2,                      &
                                        chem_gasphase_on,                 &
                                        cs_heights,                       &
                                        cs_name,                          &
                                        cs_profile,                       &
                                        cs_profile_name,                  & 
                                        cs_surface,                       &
                                        emiss_factor_main,                &
                                        emiss_factor_side,                &                      
                                        icntrl,                           &
                                        main_street_id,                   &
                                        max_street_id,                    &
                                        my_steps,                         &
                                        rcntrl,                           &
                                        side_street_id,                   &
                                        photolysis_scheme,                &
                                        wall_csflux,                      &
                                        cs_vertical_gradient,             &
                                        top_csflux,                       & 
                                        surface_csflux,                   &
                                        surface_csflux_name,              &
                                        cs_surface_initial_change,        &
                                        cs_vertical_gradient_level
                             
!-- analogue to chem_names(nspj) we could invent chem_surfaceflux(nspj) and chem_topflux(nspj)
!-- so this way we could prescribe a specific flux value for each species
!>  chemistry_parameters for initial profiles
!>  cs_names = 'O3', 'NO2', 'NO', ...   to set initial profiles)
!>  cs_heights(1,:) = 0.0, 100.0, 500.0, 2000.0, .... (height levels where concs will be prescribed for O3)
!>  cs_heights(2,:) = 0.0, 200.0, 400.0, 1000.0, .... (same for NO2 etc.) 
!>  cs_profiles(1,:) = 10.0, 20.0, 20.0, 30.0, .....  (chem spcs conc at height lvls chem_heights(1,:)) etc.
!>  If the respective concentration profile should be constant with height, then use "cs_surface( number of spcs)"
!>  then write these cs_surface values to chem_species(lsp)%conc_pr_init(:)

!
!--   Read chem namelist   
!--   (todo: initialize these parameters in declaration part, do this for
!--          all chemistry_parameters namelist parameters)
      icntrl    = 0
      rcntrl    = 0.0_wp
      my_steps  = 0.0_wp
      icntrl(2) = 1                                  
      photolysis_scheme = 'simple'
      atol = 1.0_wp
      rtol = 0.01_wp
!
!--   Try to find chemistry package
      REWIND ( 11 )
      line = ' '
      DO   WHILE ( INDEX( line, '&chemistry_parameters' ) == 0 )
         READ ( 11, '(A)', END=10 )  line
      ENDDO
      BACKSPACE ( 11 )
!
!--   Read chemistry namelist
      READ ( 11, chemistry_parameters )                          
!
!--   Enable chemistry model
      air_chemistry = .TRUE.                    

      
 10   CONTINUE

      t_steps = my_steps          !(todo: Why not directly make t_steps a
                                  !       namelist parameter?)

   END SUBROUTINE chem_parin

 !
!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine to integrate chemical species in the given chemical mechanism
!------------------------------------------------------------------------------!

   SUBROUTINE chem_integrate_ij (i, j)

      USE cpulog,                                                              &
           ONLY: cpu_log, log_point
      USE statistics,                                                          &   ! ## RFo
           ONLY:  weight_pres
       USE control_parameters,                                                 &   ! ## RFo 
           ONLY:  dt_3d, intermediate_timestep_count

      IMPLICIT   none
      INTEGER,INTENT(IN)       :: i,j

!--   local variables
      INTEGER(iwp) ::  lsp                                                          !< running index for chem spcs.
      INTEGER(iwp) ::  lph                                                          !< running index for photolysis frequencies
      INTEGER                  :: k,m,istatf                                        
      INTEGER,dimension(20)    :: istatus
      REAL(kind=wp),dimension(nzb+1:nzt,nspec)                :: tmp_conc           
      REAL(kind=wp),dimension(nzb+1:nzt)                      :: tmp_temp
      REAL(kind=wp),dimension(nzb+1:nzt,nphot)                :: tmp_phot
      REAL(kind=wp),dimension(nzb+1:nzt)                      :: tmp_fact   
     REAL(kind=wp),dimension(nzb+1:nzt)                      :: tmp_fact_i    !< conversion factor between molecules cm^{-3} and ppm
     REAL(wp)                         ::  conv                                !< conversion factor
     REAL(wp), PARAMETER              ::  ppm2fr  = 1.0e-6_wp                 !< Conversion factor ppm to fraction
     REAL(wp), PARAMETER              ::  pref_i  = 1._wp / 100000.0_wp       !< inverse reference pressure (1/Pa)
     REAL(wp), PARAMETER              ::  t_std   = 273.15_wp                 !< standard pressure (Pa)
     REAL(wp), PARAMETER              ::  p_std   = 101325.0_wp               !< standard pressure (Pa)
     REAL(wp), PARAMETER              ::  r_cp    = 0.286_wp                  !< R / cp (exponent for potential temperature)
     REAL(wp), PARAMETER              ::  vmolcm  = 22.414e3_wp               !< Mole volume (22.414 l) in cm^{-3}
     REAL(wp), PARAMETER              ::  xna     = 6.022e23_wp               !< Avogadro number (molecules/mol)


      REAL(kind=wp)  :: dt_chem                                             

       CALL cpu_log( log_point(80), '[chem_integrate_ij]', 'start' )
!<     Set chem_gasphase_on to .FALSE. if you want to skip computation of gas phase chemistry
    IF (chem_gasphase_on) THEN

       tmp_temp(:) = pt(:,j,i) * ( hyp(nzb+1:nzt) / 100000.0_wp )**0.286_wp
! ppm to molecules/cm**3
!      tmp_fact = 1.e-6_wp*6.022e23_wp/(22.414_wp*1000._wp) * 273.15_wp * hyp(nzb+1:nzt)/( 101300.0_wp * tmp_temp )
       conv = ppm2fr * xna / vmolcm
       tmp_fact(:) = conv * t_std * hyp(nzb+1:nzt) / (tmp_temp(:) * p_std)
       tmp_fact_i = 1.0_wp/tmp_fact

       CALL fill_temp (istatf, tmp_temp)                                      !< Load constant temperature into kpp context
!      CALL fill_temp (istatf, pt(nzb+1:nzt,j,i))                             !< Load temperature into kpp context

       DO lsp = 1,nspec
          tmp_conc(:,lsp) = chem_species(lsp)%conc(nzb+1:nzt,j,i) * tmp_fact(:) 
       ENDDO

       DO lph = 1,nphot
          tmp_phot(:,lph) = phot_frequen(lph)%freq(nzb+1:nzt,j,i)               
       ENDDO

       IF(myid == 0 .AND. chem_debug0 ) THEN
          IF (i == 10 .and. j == 10) WRITE(0,*) 'begin chemics step ',dt_3d
       ENDIF

!--    Compute length of time step     # RFo
       IF ( call_chem_at_all_substeps )  THEN
          dt_chem = dt_3d * weight_pres(intermediate_timestep_count)
       ELSE
          dt_chem = dt_3d
       ENDIF

       CALL cpu_log( log_point(81), '{chem_gasphase_integrate}', 'start' )


       CALL chem_gasphase_integrate (dt_chem, tmp_conc, tmp_temp, tmp_phot, istatus=istatus)


       CALL cpu_log( log_point(81), '{chem_gasphase_integrate}', 'stop' )

       DO lsp = 1,nspec
          chem_species(lsp)%conc (nzb+1:nzt,j,i) = tmp_conc(:,lsp) * tmp_fact_i(:)  ! RFo
       ENDDO

!      IF (myid == 0 .AND. chem_debug2 ) THEN 
!         IF (i == 10 .and. j == 10)   WRITE(6,'(a,8i7)') ' KPP Status ',istatus(1:8)
!         write(6,'(a,8i7)') ' KPP Status ',istatus(1:8)
!      ENDIF

    ENDIF
       CALL cpu_log( log_point(80), '[chem_integrate_ij]', 'stop' )

       RETURN
   END SUBROUTINE chem_integrate_ij
!
!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine for swapping of timelevels for chemical species
!> called out from subroutine swap_timelevel
!------------------------------------------------------------------------------!

   SUBROUTINE chem_swap_timelevel (level)
      IMPLICIT   none

      INTEGER,INTENT(IN)                  :: level

!--   local variables

      INTEGER               :: lsp

!        print*,' *** entering chem_swap_timelevel ***) '
      if(level == 0)  then
         do lsp=1, nvar                                       
            chem_species(lsp)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)    => spec_conc_1(:,:,:,lsp)
            chem_species(lsp)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  => spec_conc_2(:,:,:,lsp)
         end do
      else
         do lsp=1, nvar                                        
            chem_species(lsp)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)    => spec_conc_2(:,:,:,lsp)
            chem_species(lsp)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  => spec_conc_1(:,:,:,lsp)
         end do
      end if

      RETURN
   END SUBROUTINE chem_swap_timelevel

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!------------------------------------------------------------------------------!
   SUBROUTINE chem_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )

      IMPLICIT NONE

      CHARACTER (LEN=*), INTENT(IN)  ::  var         !<
      LOGICAL, INTENT(OUT)           ::  found       !<
      CHARACTER (LEN=*), INTENT(OUT) ::  grid_x      !<
      CHARACTER (LEN=*), INTENT(OUT) ::  grid_y      !<
      CHARACTER (LEN=*), INTENT(OUT) ::  grid_z      !<

      found  = .TRUE.

      if(var(1:3) == 'kc_')   then                    !< always the same grid for chemistry variables
            grid_x = 'x'
            grid_y = 'y'
            grid_z = 'zu'                             !< kk Use same z axis as u variables. Has to be checked if OK
      else
            found  = .FALSE.
            grid_x = 'none'
            grid_y = 'none'
            grid_z = 'none'
      end if

!     write(6,*) 'chem_define_netcdf_grid ',TRIM(var),' ',trim(grid_x),' ',found

   END SUBROUTINE chem_define_netcdf_grid
!
!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine for checking data output for chemical species
!------------------------------------------------------------------------------!

   SUBROUTINE chem_check_data_output( var, unit, i, ilen, k )


      USE control_parameters,                                                 &
         ONLY: data_output, message_string

      IMPLICIT NONE

      CHARACTER (LEN=*) ::  unit     !<
      CHARACTER (LEN=*) ::  var      !<

      INTEGER(iwp) :: i, lsp
      INTEGER(iwp) :: ilen
      INTEGER(iwp) :: k

      CHARACTER(len=16)    ::  spec_name

      unit = 'illegal'

      spec_name = TRIM(var(4:))             !< var 1:3 is 'kc_'.

       DO lsp=1,nspec
         IF (TRIM(spec_name) == TRIM(chem_species(lsp)%name))   THEN
            unit = 'ppm'
         ENDIF
         ! It is possible  to plant PM10 and PM25 into the gasphase chemistry code 
!        ! as passive species (e.g. 'passive' in GASPHASE_PREPROC/mechanisms):
!        ! set unit to micrograms per m**3 for PM10 and PM25 (PM2.5)
         IF (spec_name(1:2) == 'PM')   THEN  
            unit = 'ug m-3'
         ENDIF
       ENDDO

       DO lsp=1,nphot                                                     
         IF (TRIM(spec_name) == TRIM(phot_frequen(lsp)%name))   THEN
            unit = 'sec-1'
         ENDIF
       ENDDO


       RETURN
   END SUBROUTINE chem_check_data_output
!
!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine for checking data output of profiles for chemistry model
!------------------------------------------------------------------------------!

   SUBROUTINE chem_check_data_output_pr( variable, var_count, unit, dopr_unit )


    USE arrays_3d,                                                             &
        ONLY: zu

    USE control_parameters,                                                    &
        ONLY: data_output_pr, message_string, air_chemistry

    USE indices

    USE profil_parameter

    USE statistics


    IMPLICIT NONE

    CHARACTER (LEN=*) ::  unit     !< 
    CHARACTER (LEN=*) ::  variable !< 
    CHARACTER (LEN=*) ::  dopr_unit
    CHARACTER(len=16) ::  spec_name
 
    INTEGER(iwp) ::  var_count, lsp  !<
    

    spec_name = TRIM(variable(4:))    
!             write(9,*) 'fm #32 .. air_chemistry ', air_chemistry

          IF (  .NOT.  air_chemistry )  THEN
                 message_string = 'data_output_pr = ' //                        &
                 TRIM( data_output_pr(var_count) ) // ' is not imp' // &
                          'lemented for air_chemistry = .FALSE.'
!                CALL message( 'check_parameters', 'PA0XXX', 1, 2, 0, 6, 0 )
          ELSE
             DO lsp = 1, nspec
                IF (TRIM( spec_name ) == TRIM( chem_species(lsp)%name ) ) THEN 
                    dopr_index(var_count) = 900 
                    dopr_unit  = 'ppm'
                    hom(:,2,900,:) = SPREAD( zu, 2, statistic_regions+1 )
                ENDIF
             ENDDO
          ENDIF

   END SUBROUTINE chem_check_data_output_pr
!
!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 3D output variables for chemical species
!------------------------------------------------------------------------------!


   SUBROUTINE chem_data_output_3d( av, variable, found, local_pf, fill_value, nzb_do, nzt_do )


      USE indices

      USE kinds


      IMPLICIT NONE

      CHARACTER (LEN=*) ::  variable !<

      INTEGER(iwp) ::  av    !<
      INTEGER(iwp) ::  nzb_do !< lower limit of the data output (usually 0)
      INTEGER(iwp) ::  nzt_do !< vertical upper limit of the data output (usually nz_do3d)

      LOGICAL      ::  found !<

      REAL(wp) ::  fill_value !<
      REAL(sp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf 


      !-- local variables

      INTEGER              ::  i, j, k, lsp
      CHARACTER(len=16)    ::  spec_name


      found = .FALSE.

      spec_name = TRIM(variable(4:))
!av == 0

       DO lsp=1,nspec
          IF (TRIM(spec_name) == TRIM(chem_species(lsp)%name))   THEN
             IF(myid == 0 .AND. chem_debug0 )  WRITE(6,*) 'Output of species ' // TRIM(variable)  //  &
                                                          TRIM(chem_species(lsp)%name)       
             
             IF (av == 0) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                          local_pf(i,j,k) = MERGE(                             &
                                              chem_species(lsp)%conc(k,j,i),   &
                                              REAL( fill_value, KIND = wp ),   &
                                              BTEST( wall_flags_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
           
             ELSE
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                          local_pf(i,j,k) = MERGE(                             &
                                              chem_species(lsp)%conc_av(k,j,i),&
                                              REAL( fill_value, KIND = wp ),   &
                                              BTEST( wall_flags_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

            found = .TRUE.
          ENDIF
       ENDDO

       RETURN
   END SUBROUTINE chem_data_output_3d
!
!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine for averaging 3D data of chemical species. Due to the fact that 
!> the averaged chem arrays are allocated in chem_init, no if-query concerning
!> the allocation is required (in any mode). Attention: If you just specify an 
!> averaged output quantity in the _p3dr file during restarts the first output
!> includes the time between the beginning of the restart run and the first
!> output time (not necessarily the whole averaging_interval you have 
!> specified in your _p3d/_p3dr file )
!------------------------------------------------------------------------------!

    SUBROUTINE chem_3d_data_averaging ( mode, variable )

    USE control_parameters

    USE indices

    USE kinds

    USE surface_mod,                                                         &
        ONLY: surf_def_h, surf_lsm_h, surf_usm_h   
  
    IMPLICIT NONE
 
    CHARACTER (LEN=*) ::  mode    !< 
    CHARACTER (LEN=*) :: variable !< 
 

    INTEGER(iwp) ::  i                  !< grid index x direction
    INTEGER(iwp) ::  j                  !< grid index y direction
    INTEGER(iwp) ::  k                  !< grid index z direction
    INTEGER(iwp) ::  m                  !< running index surface type
    INTEGER(iwp) :: lsp                 !< running index for chem spcs
    INTEGER(iwp) :: lsp_2               !< it looks like redundent .. will be delted ..bK
  
    IF ( mode == 'allocate' )  THEN
       DO lsp = 1, nspec
          IF (TRIM(variable(4:)) == TRIM(chem_species(lsp)%name)) THEN
!                   lsp_2 = lsp 
             chem_species(lsp)%conc_av = 0.0_wp
             
          ENDIF
       ENDDO

    ELSEIF ( mode == 'sum' )  THEN
   
       DO lsp = 1, nspec
          IF (TRIM(variable(4:)) == TRIM(chem_species(lsp)%name)) THEN
!                   lsp_2 = lsp 
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                       chem_species(lsp)%conc_av(k,j,i) = chem_species(lsp)%conc_av(k,j,i) + &
                                                          chem_species(lsp)%conc(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( TRIM(variable(4:)) == TRIM('cssws*') )        THEN
             DO  m = 1, surf_def_h(0)%ns
                 i = surf_def_h(0)%i(m)
                 j = surf_def_h(0)%j(m)
                 chem_species(lsp)%cssws_av(j,i) = chem_species(lsp)%cssws_av(j,i) + surf_def_h(0)%cssws(lsp,m)
             ENDDO 
             DO  m = 1, surf_lsm_h%ns
                 i = surf_lsm_h%i(m)
                 j = surf_lsm_h%j(m)
                 chem_species(lsp)%cssws_av(j,i) = chem_species(lsp)%cssws_av(j,i) + surf_lsm_h%cssws(lsp,m)
             ENDDO
             DO  m = 1, surf_usm_h%ns
                 i = surf_usm_h%i(m)
                 j = surf_usm_h%j(m)
                 chem_species(lsp)%cssws_av(j,i) = chem_species(lsp)%cssws_av(j,i) + surf_usm_h%cssws(lsp,m)
             ENDDO

          ENDIF
       ENDDO
 
    ELSEIF ( mode == 'average' )  THEN
       DO lsp = 1, nspec
          IF (TRIM(variable(4:)) == TRIM(chem_species(lsp)%name)) THEN
!                   lsp_2 = lsp 
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                       chem_species(lsp)%conc_av(k,j,i) = chem_species(lsp)%conc_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          ELSEIF (TRIM(variable(4:)) == TRIM('cssws*') )            THEN
             DO i = nxlg, nxrg
                DO  j = nysg, nyng
                     chem_species(lsp)%cssws_av(j,i) = chem_species(lsp)%cssws_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO
                CALL exchange_horiz_2d( chem_species(lsp)%cssws_av, nbgp )                 
          ENDIF
       ENDDO
       
    ENDIF     

    END SUBROUTINE chem_3d_data_averaging

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine to write restart data for chemistry model
!------------------------------------------------------------------------------!
 SUBROUTINE chem_wrd_local


    IMPLICIT NONE 

    INTEGER(iwp) ::  lsp !< 

!      REAL(kind=wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg)   :: chems_conc


    DO  lsp = 1, nspec

       CALL wrd_write_string( TRIM( chem_species(lsp)%name ))
       WRITE ( 14 )  chem_species(lsp)%conc

       CALL wrd_write_string( TRIM( chem_species(lsp)%name )//'_av' )
       WRITE ( 14 )  chem_species(lsp)%conc_av

    ENDDO


 END SUBROUTINE chem_wrd_local

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine to read restart data of chemical species
!------------------------------------------------------------------------------!

 SUBROUTINE chem_rrd_local( i, k, nxlf, nxlc, nxl_on_file, nxrf, nxrc,         &
                            nxr_on_file, nynf, nync, nyn_on_file, nysf, nysc,  &
                            nys_on_file, tmp_3d, found )    
                                     
    USE control_parameters
            
    USE indices
    
    USE pegrid

    IMPLICIT NONE 

    CHARACTER (LEN=20) :: spc_name_av !<   
       
    INTEGER(iwp) ::  i, lsp          !< 
    INTEGER(iwp) ::  k               !< 
    INTEGER(iwp) ::  nxlc            !< 
    INTEGER(iwp) ::  nxlf            !< 
    INTEGER(iwp) ::  nxl_on_file     !<   
    INTEGER(iwp) ::  nxrc            !< 
    INTEGER(iwp) ::  nxrf            !< 
    INTEGER(iwp) ::  nxr_on_file     !<   
    INTEGER(iwp) ::  nync            !< 
    INTEGER(iwp) ::  nynf            !< 
    INTEGER(iwp) ::  nyn_on_file     !<   
    INTEGER(iwp) ::  nysc            !< 
    INTEGER(iwp) ::  nysf            !< 
    INTEGER(iwp) ::  nys_on_file     !<   
    
    LOGICAL, INTENT(OUT)  :: found 

    REAL(wp), DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_3d   !< 3D array to temp store data


    found = .FALSE.  


       IF ( ALLOCATED(chem_species) )  THEN

          DO  lsp = 1, nspec

              !< for time-averaged chemical conc.
             spc_name_av  =  TRIM(chem_species(lsp)%name)//'_av'

             IF (restart_string(1:length) == TRIM(chem_species(lsp)%name) )    &
             THEN
                !< read data into tmp_3d
                IF ( k == 1 )  READ ( 13 )  tmp_3d  
                !< fill ..%conc in the restart run    
                chem_species(lsp)%conc(:,nysc-nbgp:nync+nbgp,                  &
                                       nxlc-nbgp:nxrc+nbgp) =                  & 
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
                found = .TRUE.
             ELSEIF (restart_string(1:length) == spc_name_av )  THEN
                IF ( k == 1 )  READ ( 13 )  tmp_3d
                chem_species(lsp)%conc_av(:,nysc-nbgp:nync+nbgp,               &
                                          nxlc-nbgp:nxrc+nbgp) =               &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
                found = .TRUE.
             ENDIF

          ENDDO

       ENDIF


 END SUBROUTINE chem_rrd_local


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine calculating prognostic equations for chemical species 
!> (cache-optimized).
!> Routine is called separately for each chemical species over a loop from
!> prognostic_equations.
!------------------------------------------------------------------------------!
 SUBROUTINE chem_prognostic_equations_ij ( cs_scalar_p, cs_scalar, tcs_scalar_m, pr_init_cs,  &
                            i, j, i_omp_start, tn, ilsp, flux_s_cs, diss_s_cs, &
                            flux_l_cs, diss_l_cs )
    USE pegrid         
    USE advec_ws,        ONLY:  advec_s_ws 
    USE advec_s_pw_mod,  ONLY:  advec_s_pw
    USE advec_s_up_mod,  ONLY:  advec_s_up
    USE diffusion_s_mod, ONLY:  diffusion_s
    USE indices,         ONLY:  wall_flags_0
    USE surface_mod,     ONLY:  surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v, surf_usm_h,    &
                                surf_usm_v


    IMPLICIT NONE

    REAL(wp), DIMENSION(:,:,:), POINTER   :: cs_scalar_p, cs_scalar, tcs_scalar_m

    INTEGER(iwp),INTENT(IN) :: i, j, i_omp_start, tn, ilsp
    REAL(wp), DIMENSION(nzb+1:nzt,0:threads_per_task-1)         :: flux_s_cs   !<
    REAL(wp), DIMENSION(nzb+1:nzt,0:threads_per_task-1)         :: diss_s_cs   !<
    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,0:threads_per_task-1) :: flux_l_cs   !<
    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,0:threads_per_task-1) :: diss_l_cs   !<
    REAL(wp), DIMENSION(0:nz+1)                                 :: pr_init_cs  !<

!-- local variables

    INTEGER :: k
    !
    !--    Tendency-terms for chem spcs.
    tend(:,j,i) = 0.0_wp
!    
!-- Advection terms
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( ws_scheme_sca )  THEN
          CALL advec_s_ws( i, j, cs_scalar, 'kc', flux_s_cs, diss_s_cs,                &
             flux_l_cs, diss_l_cs, i_omp_start, tn )
       ELSE
          CALL advec_s_pw( i, j, cs_scalar )
       ENDIF
    ELSE
         CALL advec_s_up( i, j, cs_scalar )
    ENDIF

!

!-- Diffusion terms (the last three arguments are zero)

      CALL diffusion_s( i, j, cs_scalar, kh,                                             &
                        surf_def_h(0)%cssws(ilsp,:), surf_def_h(1)%cssws(ilsp,:),        &
                        surf_def_h(2)%cssws(ilsp,:),                                     &
                        surf_lsm_h%cssws(ilsp,:), surf_usm_h%cssws(ilsp,:),              &
                        surf_def_v(0)%cssws(ilsp,:), surf_def_v(1)%cssws(ilsp,:),        &
                        surf_def_v(2)%cssws(ilsp,:), surf_def_v(3)%cssws(ilsp,:),        &
                        surf_lsm_v(0)%cssws(ilsp,:), surf_lsm_v(1)%cssws(ilsp,:),        &
                        surf_lsm_v(2)%cssws(ilsp,:), surf_lsm_v(3)%cssws(ilsp,:),        &
                        surf_usm_v(0)%cssws(ilsp,:), surf_usm_v(1)%cssws(ilsp,:),        &
                        surf_usm_v(2)%cssws(ilsp,:), surf_usm_v(3)%cssws(ilsp,:) )
 
!    
!-- Prognostic equation for chem spcs
    DO k = nzb+1, nzt
       cs_scalar_p(k,j,i) = cs_scalar(k,j,i) + ( dt_3d  *                       &
                                               ( tsc(2) * tend(k,j,i) +         &
                                                 tsc(3) * tcs_scalar_m(k,j,i) ) & 
                                               - tsc(5) * rdf_sc(k)             &
                                                        * ( cs_scalar(k,j,i) - pr_init_cs(k) )    &   
                                               )                                                  &
                                                        * MERGE( 1.0_wp, 0.0_wp,                  &      
                                                                BTEST( wall_flags_0(k,j,i), 0 )   &             
                                                                 )       

       IF ( cs_scalar_p(k,j,i) < 0.0_wp )  cs_scalar_p(k,j,i) = 0.1_wp * cs_scalar(k,j,i)    !FKS6
    ENDDO

!
!-- Calculate tendencies for the next Runge-Kutta step
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( intermediate_timestep_count == 1 )  THEN
          DO  k = nzb+1, nzt
             tcs_scalar_m(k,j,i) = tend(k,j,i)
          ENDDO
       ELSEIF ( intermediate_timestep_count < &
          intermediate_timestep_count_max )  THEN
          DO  k = nzb+1, nzt
             tcs_scalar_m(k,j,i) = -9.5625_wp * tend(k,j,i) + &
                5.3125_wp * tcs_scalar_m(k,j,i)
          ENDDO
       ENDIF
    ENDIF

 END SUBROUTINE chem_prognostic_equations_ij


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine calculating prognostic equations for chemical species 
!> (vector-optimized).
!> Routine is called separately for each chemical species over a loop from
!> prognostic_equations.
!------------------------------------------------------------------------------!
 SUBROUTINE chem_prognostic_equations ( cs_scalar_p, cs_scalar, tcs_scalar_m,  &
                                        pr_init_cs, ilsp )

    USE advec_s_pw_mod,                                                        &
        ONLY:  advec_s_pw

    USE advec_s_up_mod,                                                        &
        ONLY:  advec_s_up

    USE advec_ws,                                                              &
        ONLY:  advec_s_ws 

    USE diffusion_s_mod,                                                       &
        ONLY:  diffusion_s

    USE indices,                                                               &
        ONLY:  nxl, nxr, nyn, nys, wall_flags_0

    USE pegrid

    USE surface_mod,                                                           &
        ONLY:  surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v, surf_usm_h,     &
               surf_usm_v

    IMPLICIT NONE

    INTEGER ::  i   !< running index
    INTEGER ::  j   !< running index
    INTEGER ::  k   !< running index

    INTEGER(iwp),INTENT(IN) ::  ilsp          !<

    REAL(wp), DIMENSION(0:nz+1) ::  pr_init_cs   !<

    REAL(wp), DIMENSION(:,:,:), POINTER ::  cs_scalar      !<
    REAL(wp), DIMENSION(:,:,:), POINTER ::  cs_scalar_p    !<
    REAL(wp), DIMENSION(:,:,:), POINTER ::  tcs_scalar_m   !<


!
!-- Tendency terms for chemical species
    tend = 0.0_wp
!    
!-- Advection terms
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( ws_scheme_sca )  THEN
          CALL advec_s_ws( cs_scalar, 'kc' )
       ELSE
          CALL advec_s_pw( cs_scalar )
       ENDIF
    ELSE
         CALL advec_s_up( cs_scalar )
    ENDIF
!
!-- Diffusion terms  (the last three arguments are zero)
    CALL diffusion_s( cs_scalar, kh,                                           &
                      surf_def_h(0)%cssws(ilsp,:),                             &
                      surf_def_h(1)%cssws(ilsp,:),                             &
                      surf_def_h(2)%cssws(ilsp,:),                             &
                      surf_lsm_h%cssws(ilsp,:),                                &
                      surf_usm_h%cssws(ilsp,:),                                &
                      surf_def_v(0)%cssws(ilsp,:),                             &
                      surf_def_v(1)%cssws(ilsp,:),                             &
                      surf_def_v(2)%cssws(ilsp,:),                             &
                      surf_def_v(3)%cssws(ilsp,:),                             &
                      surf_lsm_v(0)%cssws(ilsp,:),                             &
                      surf_lsm_v(1)%cssws(ilsp,:),                             &
                      surf_lsm_v(2)%cssws(ilsp,:),                             &
                      surf_lsm_v(3)%cssws(ilsp,:),                             &
                      surf_usm_v(0)%cssws(ilsp,:),                             &
                      surf_usm_v(1)%cssws(ilsp,:),                             &
                      surf_usm_v(2)%cssws(ilsp,:),                             &
                      surf_usm_v(3)%cssws(ilsp,:) )

!    
!-- Prognostic equation for chemical species
    DO  i = nxl, nxr
       DO  j = nys, nyn    
          DO  k = nzb+1, nzt
             cs_scalar_p(k,j,i) =   cs_scalar(k,j,i)                           &
                                  + ( dt_3d  *                                 &
                                      (   tsc(2) * tend(k,j,i)                 &
                                        + tsc(3) * tcs_scalar_m(k,j,i)         &
                                      )                                        & 
                                    - tsc(5) * rdf_sc(k)                       &
                                             * ( cs_scalar(k,j,i) - pr_init_cs(k) )    &   
                                    )                                          &
                                    * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )       

             IF ( cs_scalar_p(k,j,i) < 0.0_wp )  cs_scalar_p(k,j,i) = 0.1_wp * cs_scalar(k,j,i)
          ENDDO
       ENDDO
    ENDDO
!
!-- Calculate tendencies for the next Runge-Kutta step
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( intermediate_timestep_count == 1 )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn   
                DO  k = nzb+1, nzt
                   tcs_scalar_m(k,j,i) = tend(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ELSEIF ( intermediate_timestep_count < &
          intermediate_timestep_count_max )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   tcs_scalar_m(k,j,i) = - 9.5625_wp * tend(k,j,i)             &
                                         + 5.3125_wp * tcs_scalar_m(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF

 END SUBROUTINE chem_prognostic_equations


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining header output for chemistry model
!------------------------------------------------------------------------------!
 SUBROUTINE chem_header ( io )
       
    IMPLICIT NONE 
 
       INTEGER(iwp), INTENT(IN) ::  io            !< Unit of the output file

!      print*,'the header subroutine is still not operational'     
!!
!!--    Write chemistry model header
!       WRITE( io, 3 )
!
!       IF ( radiation_scheme == "constant" )  THEN 
!          WRITE( io, 4 ) net_radiation
!       ELSEIF ( radiation_scheme == "clear-sky" )  THEN 
!          WRITE( io, 5 )
!       ELSEIF ( radiation_scheme == "rrtmg" )  THEN 
!          WRITE( io, 6 )
!          IF ( .NOT. lw_radiation )  WRITE( io, 10 ) 
!          IF ( .NOT. sw_radiation )  WRITE( io, 11 ) 
!       ENDIF 
!
!       IF ( albedo_type == 0 )  THEN 
!          WRITE( io, 7 ) albedo
!       ELSE
!          WRITE( io, 8 ) TRIM( albedo_type_name(albedo_type) )
!       ENDIF
!       IF ( constant_albedo )  THEN 
!          WRITE( io, 9 )
!       ENDIF
!     
!       IF ( radiation .AND. radiation_scheme /= 'constant' )  THEN 
!          WRITE ( io, 1 )  lambda
!          WRITE ( io, 2 )  day_init, time_utc_init
!       ENDIF
!
!       WRITE( io, 12 ) dt_radiation
! 
! 1 FORMAT ('    Geograph. longitude            :   lambda = ',F4.1,' degr')
! 2 FORMAT ('    Day of the year at model start :   day_init = ',I3   &
!            /'    UTC time at model start        :   time_utc_init = ',F7.1' s')
! 3 FORMAT (//' Radiation model information:'/                                  &
!              ' ----------------------------'/)
! 4 FORMAT ('    --> Using constant net radiation: net_radiation = ', F6.2,        &
!           // 'W/m**2')
! 5 FORMAT ('    --> Simple radiation scheme for clear sky is used (no clouds,',   &
!                   ' default)')
! 6 FORMAT ('    --> RRTMG scheme is used')
! 7 FORMAT (/'    User-specific surface albedo: albedo =', F6.3)
! 8 FORMAT (/'    Albedo is set for land surface type: ', A)
! 9 FORMAT (/'    --> Albedo is fixed during the run')
!10 FORMAT (/'    --> Longwave radiation is disabled')
!11 FORMAT (/'    --> Shortwave radiation is disabled.')
!12 FORMAT  ('    Timestep: dt_radiation = ', F6.2, '  s')
!
!
 END SUBROUTINE chem_header

!------------------------------------------------------------------------------
! Description:
! ------------
!> Subroutine reading restart data for chemistry model input parameters
!  (FK: To make restarts work, I had to comment this routine. We actually
!       don't need it, since the namelist parameters are always read in,
!       also in case of a restart run)
!------------------------------------------------------------------------------
!  SUBROUTINE chem_rrd_global
!
!    USE chem_modules
!
!    USE control_parameters,                                                   &
!        ONLY: length, message_string, restart_string
! 
!     
!     IMPLICIT NONE
!     
!     
!     
!     DO
! 
!        SELECT CASE ( restart_string(1:length) )
!        
!           CASE ( 'bc_cs_b' )
!              READ ( 13 )  bc_cs_b
!
!           CASE DEFAULT
!
!              EXIT
!            
!        END SELECT
!        
!!
!!--     Read next string and its length
!        READ ( 13 )  length
!        READ ( 13 )  restart_string(1:length)
!        
!     ENDDO
!     
!  END SUBROUTINE chem_rrd_global


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine writing restart data for chemistry model input parameters
!  (FK: To make restarts work, I had to comment this routine. We actually
!       don't need it, since the namelist parameters are always read in,
!       also in case of a restart run)
!------------------------------------------------------------------------------!
!  SUBROUTINE chem_wrd_global
! 
!     USE chem_modules
!     
!     USE kinds
! 
! 
!     IMPLICIT NONE
! 
!     INTEGER(iwp) ::  lsp  !< running index for chem spcs
! 
! !
! !-- Writing out input parameters that are not part of chemistry_parameters 
! !-- namelist (namelist parameters are anyway read in again in case of restart)
!     DO lsp = 1, nvar
!        CALL wrd_write_string( 'conc_pr_init_'//chem_species(lsp)%name )
!        WRITE ( 14 )  chem_species(lsp)%conc_pr_init
!     ENDDO
! 
! 
!  END SUBROUTINE chem_wrd_global


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine for emission
!------------------------------------------------------------------------------!
 SUBROUTINE chem_emissions

    USE chem_modules
    
    USE netcdf_data_input_mod,                                                 &
        ONLY:  street_type_f
    
    USE surface_mod,                                                           &
        ONLY:  surf_lsm_h


    IMPLICIT NONE

    INTEGER(iwp) ::  i    !< running index for grid in x-direction
    INTEGER(iwp) ::  j    !< running index for grid in y-direction
    INTEGER(iwp) ::  m    !< running index for horizontal surfaces
    INTEGER(iwp) ::  lsp  !< running index for chem spcs

!
!-- Comment??? (todo)
    IF ( street_type_f%from_file )  THEN
!
!--    Streets are lsm surfaces, hence, no usm surface treatment required
       DO  m = 1, surf_lsm_h%ns
          i = surf_lsm_h%i(m)
          j = surf_lsm_h%j(m)
          
          IF ( street_type_f%var(j,i) >= main_street_id  .AND.                 &
               street_type_f%var(j,i) < max_street_id )  THEN
             DO  lsp = 1, nvar
                surf_lsm_h%cssws(lsp,m) = emiss_factor_main * surface_csflux(lsp)
             ENDDO
          ELSEIF ( street_type_f%var(j,i) >= side_street_id  .AND.             &
                   street_type_f%var(j,i) < main_street_id )  THEN
             DO  lsp = 1, nvar
                surf_lsm_h%cssws(lsp,m) = emiss_factor_side * surface_csflux(lsp)
             ENDDO
          ELSE
             surf_lsm_h%cssws(:,m) = 0.0_wp
          ENDIF
       ENDDO
       
    ENDIF

 END SUBROUTINE chem_emissions


 END MODULE chemistry_model_mod

