!> @file random_gauss.f90
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
! $Id: random_gauss.f90 2718 2018-01-02 08:49:38Z maronga $
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
! 1342 2014-03-26 17:04:47Z kanani
! REAL constants defined as wp-kind
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
! RCS Log replace by Id keyword, revision history cleaned up
!
! Revision 1.1  1998/03/25 20:09:47  raasch
! Initial revision
!
!
! Description:
! ------------
!> Generates a gaussian distributed random number (mean value 1, sigma = 1)
!> This routine is taken from the "numerical recipies".
!------------------------------------------------------------------------------!
 FUNCTION random_gauss( idum, upper_limit )
 

    USE kinds

    USE random_function_mod,                                                   &
        ONLY:  random_function

    IMPLICIT NONE

    INTEGER(iwp) ::  idum          !<
    INTEGER(iwp) ::  iset          !<

    REAL(wp)     ::  fac           !<
    REAL(wp)     ::  gset          !<
    REAL(wp)     ::  random_gauss  !<
    REAL(wp)     ::  rsq           !<
    REAL(wp)     ::  upper_limit   !<
    REAL(wp)     ::  v1            !<
    REAL(wp)     ::  v2            !<

    SAVE  iset, gset

    DATA  iset /0/

!
!-- Random numbers are created as long as they do not fall below the given
!-- upper limit
    DO

       IF ( iset == 0 )  THEN
          rsq = 0.0_wp
          DO  WHILE ( rsq >= 1.0_wp  .OR.  rsq == 0.0_wp )
             v1  = 2.0_wp * random_function( idum ) - 1.0_wp
             v2  = 2.0_wp * random_function( idum ) - 1.0_wp
             rsq = v1**2 + v2**2
          ENDDO
          fac          = SQRT( -2.0_wp * LOG( rsq ) / rsq )
          gset         = v1 * fac
          random_gauss = v2 * fac + 1.0_wp
          iset         = 1
       ELSE
          random_gauss = gset + 1.0_wp
          iset         = 0
       ENDIF

       IF ( ABS( random_gauss - 1.0_wp ) <= upper_limit )  EXIT

    ENDDO

 END FUNCTION random_gauss
