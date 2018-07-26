!> @file user_init.f90
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
! $Id: user_init.f90 2718 2018-01-02 08:49:38Z maronga $
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
! 1799 2016-04-05 08:35:55Z gronemeier
! Bugfix: added dots_num to variable list
!
! 1783 2016-03-06 18:36:17Z raasch
! netcdf module name changed + related changes
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute 
!
! 1320 2014-03-20 08:40:49Z raasch
! revision history before 2012 removed,
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
!> Execution of user-defined initializing actions
!------------------------------------------------------------------------------!
 SUBROUTINE user_init
 

    USE control_parameters
    
    USE indices
    
    USE kinds
    
    USE netcdf_interface,                                                      &
        ONLY: dots_label, dots_unit, dots_num
    
    USE pegrid
    
    USE user

    IMPLICIT NONE

    CHARACTER (LEN=20) :: field_char   !< 
!
!-- Here the user-defined initializing actions follow:
!-- Sample for user-defined output
!    ALLOCATE( u2(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!    ALLOCATE( ustvst(nzb:nzt+1,nysg:nyng,nxlg:nxrg) );  ustvst = 0.0_wp

!-- Sample for user-defined time series
!-- For each time series quantity you have to give a label and a unit,
!-- which will be used for the NetCDF file. They must not contain more than
!-- seven characters. The value of dots_num has to be increased by the
!-- number of new time series quantities. Its old value has to be store in
!-- dots_num_palm. See routine user_statistics on how to output calculate
!-- and output these quantities.
!    dots_label(dots_num+1) = 'abs_umx'
!    dots_unit(dots_num+1)  = 'm/s'
!    dots_label(dots_num+2) = 'abs_vmx'
!    dots_unit(dots_num+2)  = 'm/s'
!
!    dots_num_palm = dots_num
!    dots_num = dots_num + 2

 END SUBROUTINE user_init

