!> @file set_slicer_attributes_dvrp.f90
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
! $Id: set_slicer_attributes_dvrp.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2101 2017-01-05 16:42:31Z suehring
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
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
! 305 2009-04-27 11:58:42Z raasch
! Initial version
!
! Description:
! ------------
!> This routine sets the dvrp-slicer attributes
!------------------------------------------------------------------------------!
 SUBROUTINE set_slicer_attributes_dvrp( n_slicer )
 

#if defined( __dvrp_graphics )

    USE dvrp_variables,                                                        &
        ONLY:  dvrp_colortable_entries, interval_h_dvrp, interval_values_dvrp, &
               slicer_range_limits_dvrp

    USE kinds

    IMPLICIT NONE

    INTEGER(iwp) ::  j         !<
    INTEGER(iwp) ::  n_slicer  !<

    REAL(wp)     ::  maxv      !<
    REAL(wp)     ::  meav      !<
    REAL(wp)     ::  minv      !<


!
!-- Set interval values to user settings.
!-- The middle of this interval defines the change from blue to yellow
    minv = slicer_range_limits_dvrp(1,n_slicer)
    maxv = slicer_range_limits_dvrp(2,n_slicer)
    meav = ( minv + maxv ) * 0.5_wp

!
!-- Create appropriate colortable with 100 entries.
!-- This table ranges from deep blue (min) to deep red (max)
    DO  j = 1, 50
       interval_values_dvrp(1,j) = minv + (meav-minv) * ( j - 1.0_wp ) / 50.0_wp
       interval_values_dvrp(2,j) = minv + (meav-minv) * ( j )          / 50.0_wp
       interval_h_dvrp(:,j) = 270.0_wp - ( j - 1.0_wp ) * 90.0_wp / 49.0_wp
    ENDDO

    DO  j = 51, 100
       interval_values_dvrp(1,j) = meav + (maxv-meav) * ( j - 51.0_wp ) / 50.0_wp
       interval_values_dvrp(2,j) = meav + (maxv-meav) * ( j - 50.0_wp ) / 50.0_wp
       interval_h_dvrp(:,j) = 70.0_wp - ( j - 51.0_wp ) * 90.0_wp / 49.0_wp
    ENDDO

    dvrp_colortable_entries = 100

#endif

 END SUBROUTINE set_slicer_attributes_dvrp
