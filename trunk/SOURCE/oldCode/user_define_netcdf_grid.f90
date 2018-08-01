!> @file user_define_netcdf_grid.f90
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
! $Id: user_define_netcdf_grid.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2271 2017-06-09 12:34:55Z sward
! found = .TRUE. added to grid definition in examples to avoid warning message
! 
! 2101 2017-01-05 16:42:31Z suehring
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
!
! 1320 2014-03-20 08:40:49Z raasch
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 211 2008-11-11 04:46:24Z raasch
! Former file user_interface.f90 split into one file per subroutine
!
! Description:
! ------------
!> Set the grids on which user-defined output quantities are defined.
!> Allowed values for grid_x are "x" and "xu", for grid_y "y" and "yv", and
!> for grid_z "zu" and "zw".
!------------------------------------------------------------------------------!
 SUBROUTINE user_define_netcdf_grid( variable, found, grid_x, grid_y, grid_z )
 

    USE kinds
    
    USE user

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  grid_x     !< 
    CHARACTER (LEN=*) ::  grid_y     !< 
    CHARACTER (LEN=*) ::  grid_z     !< 
    CHARACTER (LEN=*) ::  variable   !< 

    LOGICAL ::  found   !< 


    SELECT CASE ( TRIM( variable ) )

!
!--    Uncomment and extend the following lines, if necessary
!       CASE ( 'u2', 'u2_xy', 'u2_xz', 'u2_yz' )
!          found  = .TRUE.
!          grid_x = 'xu'
!          grid_y = 'y'
!          grid_z = 'zu'

!       CASE ( 'u*v*', 'u*v*_xy', 'u*v*_xz', 'u*v*_yz' )
!          found  = .TRUE.
!          grid_x = 'x'
!          grid_y = 'y'
!          grid_z = 'zu'

       CASE DEFAULT
          found  = .FALSE.
          grid_x = 'none'
          grid_y = 'none'
          grid_z = 'none'

    END SELECT


 END SUBROUTINE user_define_netcdf_grid

