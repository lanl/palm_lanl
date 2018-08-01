!> @file chem_modules.f90
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
! $Id: chem_modules.f90 2766 2018-01-22 17:17:47Z kanani $
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
! @author Farah Kanani-Suehring
! @author Basit Khan
!
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of global palm-4u chemistry variables
!> (Module written to define global palm-4u chemistry variables. basit 16Nov2017)
!------------------------------------------------------------------------------!
!
 MODULE chem_modules

    USE kinds 
    USE chem_gasphase_mod,                                                     &   
        ONLY: nspec, nvar, spc_names


    IMPLICIT NONE

    PUBLIC nspec
    PUBLIC nvar
    PUBLIC spc_names

    LOGICAL ::  constant_top_csflux(99)            = .TRUE.                       !< chem spcs at the top  orig .TRUE.
    LOGICAL ::  constant_csflux(99)                = .TRUE.                       !< chem spcs at namelist parameter   orig TRUE


    INTEGER(iwp) :: cs_vertical_gradient_level_ind(99,10) = -9999                 !< grid index values of cs_vertical_gradient_level_ind(s)
    REAL(wp),     DIMENSION(:),   ALLOCATABLE :: bc_cs_t_val
    INTEGER(iwp)                              :: ibc_cs_b                         !< integer flag for bc_cs_b
    INTEGER(iwp)                              :: ibc_cs_t                         !< integer flag for bc_cs_t

    REAL(wp),      DIMENSION(:),  ALLOCATABLE ::  css                             !< scaling parameter for chem spcs

!-- Namelist parameters for creating initial chemistry profiles
    CHARACTER (LEN=20)               :: bc_cs_b    = 'dirichlet'                  !< namelist parameter
    CHARACTER (LEN=20)               :: bc_cs_t    = 'initial_gradient'           !< namelist parameter
    REAL(wp) :: wall_csflux (99,0:5)               = 0.0_wp                       !< namelist parameter
    REAL(wp) :: cs_vertical_gradient (99,10)       = 0.0_wp                       !< namelist parameter
    REAL(wp) :: cs_vertical_gradient_level (99,10) = -999999.9_wp                 !< namelist parameter
    REAL(wp) :: top_csflux ( 99 )                  = 0.0_wp                       !< namelist parameter 
    REAL(wp) :: cs_surface_initial_change(99)      = 0.0_wp                       !< namelist parameter
    REAL(wp) :: surface_csflux(99 )                = 0.0_wp                       !< namelist parameter: fluxes where 'surface_csflux_name' is in the namelist
!   RFo: I do not know whether it makes sense to have 'constant_csflux=.TRUE. for only these species where 
!        no flux is given in the namelist. Let's choos surface_csflux=0.0 (and thus 'constant_csflux'=.TRUE.) as default
!       To obtain  constant_csflux=.FALSE., set surface_csflux = 9999999.9 in the namelist
!   @todo: need to think a bit more about constant_csflux for chemistry. 

    LOGICAL :: call_chem_at_all_substeps           = .FALSE.                      !< namelist parameter 
    LOGICAL :: chem_debug0                         = .FALSE.                      !< namelist parameter flag for minimum print output
    LOGICAL :: chem_debug1                         = .FALSE.                      !< namelist parameter flag for print output
    LOGICAL :: chem_debug2                         = .FALSE.                      !< namelist parameter flag for further print output
    LOGICAL :: chem_gasphase_on                    = .TRUE.                       !< namelist parameter 

    CHARACTER (LEN=11), DIMENSION(99)         :: cs_name = 'novalue'              !< Namelist parameter: chem spcs names 
    CHARACTER (LEN=11), DIMENSION(99)         :: cs_profile_name = 'novalue'      !< Namelist parameter: Names of the
    CHARACTER (LEN=11), DIMENSION(99)         :: surface_csflux_name = 'novalue'  !< Namelist parameter: chem species surface fluxes names
                                                                                  !< active chem spcs, default is 'novalue')  ????
    REAL(wp), DIMENSION(99)                   :: cs_surface = 0.0_wp              !< Namelist parameter: Surface conc of chem spcs'
    REAL(wp), DIMENSION(99,100)               :: cs_heights = 9999999.9_wp        !< Namelist parameter: Height lvls(m) for cs_profiles
    REAL(wp), DIMENSION(99,100)               :: cs_profile = 9999999.9_wp        !< Namelist parameter: Chem conc for each spcs defined

#if defined( __nopointer )
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET   :: cs                       !< chem spcs
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET   :: cs_p                     !< prognostic value of chem spc
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET   :: tcs_m                    !< weighted tendency of s for previous sub-timestep (Runge-Kutta)

#else                                                               
! use pointers cs, cs_p and tcs_m to point arrays cs_1, cs_2, and cs_3
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: cs_1                     !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: cs_2                     !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: cs_3                     !< pointer for swapping of timelevels for respective quantity

    REAL(wp), DIMENSION(:,:,:), POINTER               :: cs                       !< pointer: sgs chem spcs)
    REAL(wp), DIMENSION(:,:,:), POINTER               :: cs_p                     !< pointer: prognostic value of sgs chem spcs 
    REAL(wp), DIMENSION(:,:,:), POINTER               :: tcs_m                    !< pointer:

#endif                                                                            
                                                                                  !< by cs_name for each height lvls defined by cs_heights
!
!-- Namelist parameters for chem_emissions
    INTEGER(iwp) ::  main_street_id = 0
    INTEGER(iwp) ::  max_street_id = 0
    INTEGER(iwp) ::  side_street_id = 0
!
!-- Constant emission factors
    REAL(wp) ::  emiss_factor_main = 0.0_wp
    REAL(wp) ::  emiss_factor_side = 0.0_wp
    
!-- Emission factors with daily cycle
!     REAL(wp), DIMENSION(1:24) ::  emiss_factor_main = 0.0_wp
!     REAL(wp), DIMENSION(1:24) ::  emiss_factor_side = 0.0_wp

    SAVE

 END MODULE chem_modules

