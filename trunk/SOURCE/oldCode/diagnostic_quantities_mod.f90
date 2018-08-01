!> @file diagnostic_quantities_mod.f90
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
! -----------------
! 
! 
! Former revisions:
! ----------------- 
! $Id: diagnostic_quantities_mod.f90 3026 2018-05-22 10:30:53Z schwenkel $
! Changed the name specific humidity to mixing ratio, since we are computing
! mixing ratios.
! 
! 2839 2018-02-27 09:49:06Z schwenkel
! Bugfix for Kessler microphysics
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2608 2017-11-13 14:04:26Z schwenkel
! Initial revision
! 
! 
! Description:
! ------------
!> This module contains subroutines and functions for the calculation of 
!> diagnostic quantities. Especially moisture quantities such as the saturation
!> mixining ratio is calculated
!------------------------------------------------------------------------------!
MODULE diagnostic_quantities_mod
 

   USE kinds
   
   IMPLICIT NONE

   REAL(wp) ::  alpha   !< correction factor 
   REAL(wp) ::  e_s     !< saturation water vapor pressure
   REAL(wp) ::  q_s     !< saturation mixing ratio
   REAL(wp) ::  sat     !< supersaturation
   REAL(wp) ::  t_l     !< actual temperature

   PRIVATE
   PUBLIC  e_s, magnus, q_s, sat, supersaturation, t_l

    INTERFACE supersaturation
       MODULE PROCEDURE supersaturation
    END INTERFACE supersaturation

 CONTAINS
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computation of the diagnostic supersaturation sat, actual temperature t_l
!< and saturation water vapor mixing ratio q_
!------------------------------------------------------------------------------!
    SUBROUTINE supersaturation ( i,j,k )

       USE arrays_3d,                                                          &
           ONLY:  hyp, pt, q, qc, qr

       USE cloud_parameters,                                                   &
           ONLY:  l_d_cp, l_d_r, t_d_pt       

       USE control_parameters,                                                 &
           ONLY:  microphysics_kessler                                       

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<
!
!--    Actual liquid water temperature:
       t_l = t_d_pt(k) * pt(k,j,i)
!
!--    Calculate water vapor saturation pressure
       e_s = magnus( t_l )
!
!--    Computation of saturation mixing ratio:
       q_s   = 0.622_wp * e_s / ( hyp(k) - e_s )
!
!--    Correction factor
       alpha = 0.622_wp * l_d_r * l_d_cp / ( t_l * t_l )
!
!--    Correction of the approximated value
!--    (see: Cuijpers + Duynkerke, 1993, JAS, 23)
       q_s   = q_s * ( 1.0_wp + alpha * q(k,j,i) ) / ( 1.0_wp + alpha * q_s )

!
!--    Supersaturation:
!--    Not in case of microphysics_kessler since qr is unallocated
       IF ( .NOT. microphysics_kessler ) THEN
          sat   = ( q(k,j,i) - qr(k,j,i) - qc(k,j,i) ) / q_s - 1.0_wp
       ENDIF

    END SUBROUTINE supersaturation

 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> This function computes the magnus function (Press et al., 1992).
!> The magnus formula is needed to calculate the saturation vapor pressure
!------------------------------------------------------------------------------!

    FUNCTION magnus( temperature )

       IMPLICIT NONE

       REAL(wp)     ::  magnus            !<
       REAL(wp)     ::  temperature       !<

!
!--    Saturation vapor pressure at t_l:
       magnus =  611.2_wp * EXP( 17.62_wp * ( temperature - 273.15_wp ) /      & 
                                            ( temperature - 29.65_wp  ) )

!        magnus = 610.78_wp * EXP( 17.269_wp * ( temperature - 273.16_wp ) /     &
!                                              ( temperature - 35.86_wp )        &
!                                )

    END FUNCTION magnus


END MODULE diagnostic_quantities_mod
