!> @file random_function_mod.f90
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
! $Id: random_function_mod.f90 2718 2018-01-02 08:49:38Z maronga $
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
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
! 
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
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! RCS Log replace by Id keyword, revision history cleaned up
!
! Revision 1.1  1998/02/04 16:09:45  raasch
! Initial revision
!
!
! Description:
! ------------
!> Random number generator, produces numbers equally distributed in interval [0,1]
!> This routine is taken from the "numerical recipies"
!------------------------------------------------------------------------------!
 MODULE random_function_mod
 

    USE kinds

    IMPLICIT NONE

    PRIVATE

    PUBLIC random_function, random_function_ini

    INTEGER(iwp), PUBLIC, SAVE ::  random_iv(32)  !<
    INTEGER(iwp), PUBLIC, SAVE ::  random_iy      !<

    INTERFACE random_function_ini
       MODULE PROCEDURE random_function_ini
    END INTERFACE random_function_ini

    INTERFACE random_function
       MODULE PROCEDURE random_function
    END INTERFACE random_function

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE random_function_ini

       IMPLICIT NONE

       random_iv = 0
       random_iy = 0

    END SUBROUTINE random_function_ini

    FUNCTION random_function( idum )


       IMPLICIT NONE

       INTEGER(iwp) ::  ia               !<
       INTEGER(iwp) ::  idum             !<
       INTEGER(iwp) ::  im               !<
       INTEGER(iwp) ::  iq               !<
       INTEGER(iwp) ::  ir               !<
       INTEGER(iwp) ::  ndiv             !<
       INTEGER(iwp) ::  ntab             !<

       INTEGER(iwp) ::  j                !<
       INTEGER(iwp) ::  k                !<

       REAL(wp)     ::  am               !<
       REAL(wp)     ::  eps              !<
       REAL(wp)     ::  random_function  !<
       REAL(wp)     ::  rnmx             !<

       PARAMETER ( ia=16807, im=2147483647, am=1.0_wp/im, iq=127773, ir=2836, &
                   ntab=32, ndiv=1+(im-1)/ntab, eps=1.2e-7_wp, rnmx=1.0_wp-eps )

       IF ( idum .le. 0  .or.  random_iy .eq. 0 )  THEN
          idum = max (-idum,1)
          DO  j = ntab+8,1,-1
             k    = idum / iq
             idum = ia * ( idum - k * iq ) - ir * k
             IF ( idum .lt. 0 )  idum = idum + im
             IF ( j .le. ntab )  random_iv(j) = idum
          ENDDO
          random_iy = random_iv(1)
       ENDIF

       k    = idum / iq
       idum = ia * ( idum - k * iq ) - ir * k
       IF ( idum .lt. 0 )  idum = idum + im
       j            = 1 + random_iy / ndiv
       random_iy    = random_iv(j)
       random_iv(j) = idum
       random_function  = min ( am * random_iy , rnmx )

    END FUNCTION random_function

 END MODULE random_function_mod
