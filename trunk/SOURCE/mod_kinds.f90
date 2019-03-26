!> @file mod_kinds.f90
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
! ------------------
!
!
! Former revisions:
! -----------------
! $Id: mod_kinds.f90 2718 2018-01-02 08:49:38Z maronga $
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
! 1306 2014-03-13 14:30:59Z raasch
! Initial revision
!
! Description:
! ------------
!> Standard kind definitions
!> wp (working precision) and iwp (integer working precision) are the kinds
!> used by default in all variable declarations.
!> By default, PALM is using wp = dp (64bit), and iwp = isp (32bit).
!> If you like to switch to other precision, then please set wp/iwp
!> appropriately by assigning other kinds below.
!------------------------------------------------------------------------------!
 MODULE kinds


    IMPLICIT NONE

!
!-- Floating point kinds
    INTEGER, PARAMETER ::  sp = 4           !< single precision (32 bit)
    INTEGER, PARAMETER ::  dp = 8           !< double precision (64 bit)

!
!-- Integer kinds
    INTEGER, PARAMETER ::  isp = SELECTED_INT_KIND(  9 )   !< single precision (32 bit)
    INTEGER, PARAMETER ::  idp = SELECTED_INT_KIND( 14 )   !< double precision (64 bit)

!
!-- Set kinds to be used as defaults
    INTEGER, PARAMETER ::   wp =  dp          !< default real kind
    INTEGER, PARAMETER ::  iwp = isp          !< default integer kind

    SAVE

 END MODULE kinds
