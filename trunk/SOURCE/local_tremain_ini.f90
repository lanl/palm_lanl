!> @file local_tremain_ini.f90
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
! $Id: local_tremain_ini.f90 2718 2018-01-02 08:49:38Z maronga $
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
! 1808 2016-04-05 19:44:00Z raasch
! cpu measurements are done with standard FORTRAN routine on every machine
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL constants defined as wp_kind
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! Revision 1.1  1998/03/18 20:15:05  raasch
! Initial revision
!
!
! Description:
! ------------
!> Initialization of CPU-time measurements for different operating systems
!------------------------------------------------------------------------------!
 SUBROUTINE local_tremain_ini
 
      
    USE cpulog,                                                                &
        ONLY:  initial_wallclock_time
        
    USE kinds

    IMPLICIT NONE

    INTEGER(idp)     ::  count      !<
    INTEGER(idp)     ::  count_rate !<

!
!-- Get initial wall clock time
    CALL SYSTEM_CLOCK( count, count_rate )
    initial_wallclock_time = REAL( count, KIND=wp ) / REAL( count_rate, KIND=wp )

 END SUBROUTINE local_tremain_ini
