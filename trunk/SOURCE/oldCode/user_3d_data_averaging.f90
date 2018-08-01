!> @file user_3d_data_averaging.f90
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
! $Id: user_3d_data_averaging.f90 3004 2018-04-27 12:33:25Z Giersch $
! Further allocation checks implemented
! 
! 2718 2018-01-02 08:49:38Z maronga
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
! 1322 2014-03-20 16:38:49Z raasch
! REAL functions provided with KIND-attribute
!
! 1320 2014-03-20 08:40:49Z raasch
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
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
!> Sum up and time-average user-defined output quantities as well as allocate
!> the array necessary for storing the average.
!------------------------------------------------------------------------------!
 SUBROUTINE user_3d_data_averaging( mode, variable )
 

    USE control_parameters

    USE indices

    USE kinds

    USE user

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  mode    !< 
    CHARACTER (LEN=*) :: variable !< 

    INTEGER(iwp) ::  i !< 
    INTEGER(iwp) ::  j !< 
    INTEGER(iwp) ::  k !< 

    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( variable ) )

!
!--       Uncomment and extend the following lines, if necessary.
!--       The arrays for storing the user defined quantities (here u2_av) have
!--       to be declared and defined by the user!
!--       Sample for user-defined output:
!          CASE ( 'u2' )
!             IF ( .NOT. ALLOCATED( u2_av ) )  THEN
!                ALLOCATE( u2_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!             ENDIF
!             u2_av = 0.0_wp

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( variable ) )

!
!--       Uncomment and extend the following lines, if necessary.
!--       The arrays for storing the user defined quantities (here u2 and
!--       u2_av) have to be declared and defined by the user!
!--       Sample for user-defined output:
!          CASE ( 'u2' )
!             IF ( ALLOCATED( u2_av ) ) THEN
!                DO  i = nxlg, nxrg
!                   DO  j = nysg, nyng
!                      DO  k = nzb, nzt+1
!                         u2_av(k,j,i) = u2_av(k,j,i) + u2(k,j,i)
!                      ENDDO
!                   ENDDO
!                ENDDO
!             ENDIF

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( variable ) )

!
!--       Uncomment and extend the following lines, if necessary.
!--       The arrays for storing the user defined quantities (here u2_av) have
!--       to be declared and defined by the user!
!--       Sample for user-defined output:
!          CASE ( 'u2' )
!             IF ( ALLOCATED( u2_av ) ) THEN
!                DO  i = nxlg, nxrg
!                   DO  j = nysg, nyng
!                      DO  k = nzb, nzt+1
!                         u2_av(k,j,i) = u2_av(k,j,i) / REAL( average_count_3d, KIND=wp )
!                      ENDDO
!                   ENDDO
!                ENDDO
!             ENDIF

       END SELECT

    ENDIF


 END SUBROUTINE user_3d_data_averaging
