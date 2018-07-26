!> @file user_module.f90
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
! $Id: user_module.f90 2894 2018-03-15 09:17:58Z Giersch $
! An example for a user defined global variable has been added (Giersch)
! 
! 2718 2018-01-02 08:49:38Z suehring
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
! 1873 2016-04-18 14:50:06Z maronga
! Module renamed (removed _mod)
! 
! 
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
! 
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
!
! 1320 2014-03-20 08:40:49Z raasch
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
! Revision 1.1  1998/03/24 15:29:04  raasch
! Initial revision
!
!
! Description:
! ------------
!> Declaration of user-defined variables. This module may only be used
!> in the user-defined routines (contained in user_interface.f90).
!------------------------------------------------------------------------------!
 MODULE user
 

    USE kinds

    IMPLICIT NONE

    INTEGER(iwp) ::  dots_num_palm   !< 
    INTEGER(iwp) ::  user_idummy     !< 
    
    LOGICAL ::  user_defined_namelist_found = .FALSE.   !< 
    
    REAL(wp) ::  user_rdummy   !< 

!
!-- Sample for user-defined output
!    REAL(wp) :: global_parameter !< user defined global parameter 
!
!    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u2       !< user defined array
!    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u2_av    !< user defined array
!    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  ustvst   !< user defined array

    SAVE

 END MODULE user
