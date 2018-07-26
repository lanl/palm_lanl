!> @file init_advec.f90
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
! $Id: init_advec.f90 2718 2018-01-02 08:49:38Z maronga $
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
! 1346 2014-03-27 13:18:20Z heinze
! Bugfix: REAL constants provided with KIND-attribute especially in call of 
! intrinsic function like MAX, MIN, SIGN
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL constants defined as wp-kind
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements, 
! kinds are defined in new module kinds, 
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
! 
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1003 2012-09-14 14:35:53Z raasch
! obsolete variables removed
!
! 1001 2012-09-13 14:08:46Z raasch
! all actions concerning upstream-spline-method removed
!
! Revision 1.1  1999/02/05 09:07:38  raasch
! Initial revision
!
!
! Description:
! ------------
!> Initialize constant coefficients and parameters for certain advection schemes.
!------------------------------------------------------------------------------!
 SUBROUTINE init_advec
 

    USE advection,                                                             &
        ONLY:  aex, bex, dex, eex
        
    USE kinds
    
    USE control_parameters,                                                    &
        ONLY:  scalar_advec

    IMPLICIT NONE

    INTEGER(iwp) ::  i          !<
    INTEGER(iwp) ::  intervals  !< 
    INTEGER(iwp) ::  j          !<
    
    REAL(wp) :: delt   !<
    REAL(wp) :: dn     !<
    REAL(wp) :: dnneu  !<
    REAL(wp) :: ex1    !<
    REAL(wp) :: ex2    !<
    REAL(wp) :: ex3    !<
    REAL(wp) :: ex4    !<
    REAL(wp) :: ex5    !<
    REAL(wp) :: ex6    !<
    REAL(wp) :: sterm  !<


    IF ( scalar_advec == 'bc-scheme' )  THEN

!
!--    Compute exponential coefficients for the Bott-Chlond scheme
       intervals = 1000
       ALLOCATE( aex(intervals), bex(intervals), dex(intervals), eex(intervals) )

       delt  = 1.0_wp / REAL( intervals, KIND=wp )
       sterm = delt * 0.5_wp

       DO  i = 1, intervals

          IF ( sterm > 0.5_wp )  THEN
             dn = -5.0_wp
          ELSE
             dn = 5.0_wp
          ENDIF

          DO  j = 1, 15
             ex1 = dn * EXP( -dn ) - EXP( 0.5_wp * dn ) + EXP( -0.5_wp * dn )
             ex2 = EXP( dn ) - EXP( -dn )
             ex3 = EXP( -dn ) * ( 1.0_wp - dn ) - 0.5_wp * EXP(  0.5_wp * dn ) &
                                                - 0.5_wp * EXP( -0.5_wp * dn )
             ex4 = EXP( dn ) + EXP( -dn )
             ex5 = dn * sterm + ex1 / ex2
             ex6 = sterm + ( ex3 * ex2 - ex4 * ex1 ) / ( ex2 * ex2 )
             dnneu = dn - ex5 / ex6
             dn  = dnneu
          ENDDO

          IF ( sterm < 0.5_wp )  dn = MAX(  2.95E-2_wp, dn )
          IF ( sterm > 0.5_wp )  dn = MIN( -2.95E-2_wp, dn )
          ex1 = EXP( -dn )
          ex2 = EXP( dn ) - ex1
          aex(i) = -ex1 / ex2
          bex(i) = 1.0_wp / ex2
          dex(i) = dn
          eex(i) = EXP( dex(i) * 0.5_wp )
          sterm = sterm + delt

       ENDDO

    ENDIF

 END SUBROUTINE init_advec
