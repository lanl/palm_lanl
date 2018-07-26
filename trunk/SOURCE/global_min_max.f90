!> @file global_min_max.f90
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
! $Id: global_min_max.f90 2718 2018-01-02 08:49:38Z maronga $
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
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
! 
! 1188 2013-06-20 12:00:08Z heinze
! Bugfix in modes 'min' and 'max': x and z component were interchanged
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 866 2012-03-28 06:44:41Z raasch
! new mode "absoff" accounts for an offset in the respective array
!
! Revision 1.1  1997/07/24 11:14:03  raasch
! Initial revision
!
!
! Description:
! ------------
!> Determine the array minimum/maximum and the corresponding indices.
!------------------------------------------------------------------------------!
 SUBROUTINE global_min_max( i1, i2, j1, j2, k1, k2, ar, mode, offset, value, &
                            value_ijk, value1, value1_ijk )
 

    USE indices,                                                               &
        ONLY:  nbgp, ny, nx
        
    USE kinds
    
    USE pegrid

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  mode  !<

    INTEGER(iwp) ::  i              !<
    INTEGER(iwp) ::  i1             !<
    INTEGER(iwp) ::  i2             !<
    INTEGER(iwp) ::  id_fmax        !<
    INTEGER(iwp) ::  id_fmin        !<
    INTEGER(iwp) ::  j              !<
    INTEGER(iwp) ::  j1             !<
    INTEGER(iwp) ::  j2             !<
    INTEGER(iwp) ::  k              !<
    INTEGER(iwp) ::  k1             !<
    INTEGER(iwp) ::  k2             !<
    INTEGER(iwp) ::  fmax_ijk(3)    !<
    INTEGER(iwp) ::  fmax_ijk_l(3)  !<
    INTEGER(iwp) ::  fmin_ijk(3)    !<
    INTEGER(iwp) ::  fmin_ijk_l(3)  !<
    INTEGER(iwp) ::  value_ijk(3)   !<
    
    INTEGER(iwp), OPTIONAL ::  value1_ijk(3)  !<
    
    REAL(wp) ::  offset                 !<
    REAL(wp) ::  value                  !<
    REAL(wp) ::  ar(i1:i2,j1:j2,k1:k2)  !<
    
#if defined( __ibm )
    REAL(sp) ::  fmax(2)    !<
    REAL(sp) ::  fmax_l(2)  !<
    REAL(sp) ::  fmin(2)    !<
    REAL(sp) ::  fmin_l(2)  !<
             ! on 32bit-machines MPI_2REAL must not be replaced 
             ! by MPI_2DOUBLE_PRECISION
#else
    REAL(wp) ::  fmax(2)    !<
    REAL(wp) ::  fmax_l(2)  !<
    REAL(wp) ::  fmin(2)    !<
    REAL(wp) ::  fmin_l(2)  !<
#endif
    REAL(wp), OPTIONAL ::  value1  !<


!
!-- Determine array minimum
    IF ( mode == 'min'  .OR.  mode == 'minmax' )  THEN

!
!--    Determine the local minimum
       fmin_ijk_l = MINLOC( ar )
       fmin_ijk_l(1) = i1 + fmin_ijk_l(1) - 1 ! MINLOC assumes lowerbound = 1
       fmin_ijk_l(2) = j1 + fmin_ijk_l(2) - nbgp
       fmin_ijk_l(3) = k1 + fmin_ijk_l(3) - nbgp
       fmin_l(1)  = ar(fmin_ijk_l(1),fmin_ijk_l(2),fmin_ijk_l(3))

#if defined( __parallel )
       fmin_l(2)  = myid
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( fmin_l, fmin, 1, MPI_2REAL, MPI_MINLOC, comm2d, &
                           ierr )

!
!--    Determine the global minimum. Result stored on PE0.
       id_fmin = fmin(2)
       IF ( id_fmin /= 0 )  THEN
          IF ( myid == 0 )  THEN
             CALL MPI_RECV( fmin_ijk, 3, MPI_INTEGER, id_fmin, 0, comm2d, &
                            status, ierr )
          ELSEIF ( myid == id_fmin )  THEN
             CALL MPI_SEND( fmin_ijk_l, 3, MPI_INTEGER, 0, 0, comm2d, ierr )
          ENDIF
       ELSE
          fmin_ijk = fmin_ijk_l
       ENDIF
!
!--    Send the indices of the just determined array minimum to other PEs
       CALL MPI_BCAST( fmin_ijk, 3, MPI_INTEGER, 0, comm2d, ierr )
#else
       fmin(1)  = fmin_l(1)
       fmin_ijk = fmin_ijk_l
#endif

    ENDIF

!
!-- Determine array maximum
    IF ( mode == 'max'  .OR.  mode == 'minmax' )  THEN

!
!--    Determine the local maximum
       fmax_ijk_l = MAXLOC( ar )
       fmax_ijk_l(1) = i1 + fmax_ijk_l(1) - 1 ! MAXLOC assumes lowerbound = 1
       fmax_ijk_l(2) = j1 + fmax_ijk_l(2) - nbgp
       fmax_ijk_l(3) = k1 + fmax_ijk_l(3) - nbgp
       fmax_l(1) = ar(fmax_ijk_l(1),fmax_ijk_l(2),fmax_ijk_l(3))

#if defined( __parallel )
       fmax_l(2)  = myid
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( fmax_l, fmax, 1, MPI_2REAL, MPI_MAXLOC, comm2d, &
                           ierr )

!
!--    Determine the global maximum. Result stored on PE0.
       id_fmax = fmax(2)
       IF ( id_fmax /= 0 )  THEN
          IF ( myid == 0 )  THEN
             CALL MPI_RECV( fmax_ijk, 3, MPI_INTEGER, id_fmax, 0, comm2d, &
                            status, ierr )
          ELSEIF ( myid == id_fmax )  THEN
             CALL MPI_SEND( fmax_ijk_l, 3, MPI_INTEGER, 0, 0, comm2d, ierr )
          ENDIF
       ELSE
          fmax_ijk = fmax_ijk_l
       ENDIF
!
!--    send the indices of the just determined array maximum to other PEs
       CALL MPI_BCAST( fmax_ijk, 3, MPI_INTEGER, 0, comm2d, ierr )
#else
       fmax(1)  = fmax_l(1)
       fmax_ijk = fmax_ijk_l
#endif

    ENDIF

!
!-- Determine absolute array maximum
    IF ( mode == 'abs' )  THEN

!
!--    Determine the local absolut maximum
       fmax_l(1)     = 0.0_wp
       fmax_ijk_l(1) =  i1
       fmax_ijk_l(2) =  j1
       fmax_ijk_l(3) =  k1
       DO  k = k1, k2
          DO  j = j1, j2
             DO  i = i1, i2
                IF ( ABS( ar(i,j,k) ) > fmax_l(1) )  THEN
                   fmax_l(1) = ABS( ar(i,j,k) )
                   fmax_ijk_l(1) = i
                   fmax_ijk_l(2) = j
                   fmax_ijk_l(3) = k
                ENDIF
             ENDDO
          ENDDO
       ENDDO

!
!--    Set a flag in case that the determined value is negative.
!--    A constant offset has to be subtracted in order to handle the special
!--    case i=0 correctly
       IF ( ar(fmax_ijk_l(1),fmax_ijk_l(2),fmax_ijk_l(3)) < 0.0_wp )  THEN
          fmax_ijk_l(1) = -fmax_ijk_l(1) - 10
       ENDIF

#if defined( __parallel )
       fmax_l(2)  = myid
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( fmax_l, fmax, 1, MPI_2REAL, MPI_MAXLOC, comm2d, &
                           ierr )

!
!--    Determine the global absolut maximum. Result stored on PE0.
       id_fmax = fmax(2)
       IF ( id_fmax /= 0 )  THEN
          IF ( myid == 0 )  THEN
             CALL MPI_RECV( fmax_ijk, 3, MPI_INTEGER, id_fmax, 0, comm2d, &
                            status, ierr )
          ELSEIF ( myid == id_fmax )  THEN
             CALL MPI_SEND( fmax_ijk_l, 3, MPI_INTEGER, 0, 0, comm2d, ierr )
          ENDIF
       ELSE
          fmax_ijk = fmax_ijk_l
       ENDIF
!
!--    Send the indices of the just determined absolut maximum to other PEs
       CALL MPI_BCAST( fmax_ijk, 3, MPI_INTEGER, 0, comm2d, ierr )
#else
       fmax(1)  = fmax_l(1)
       fmax_ijk = fmax_ijk_l
#endif

    ENDIF

!
!-- Determine absolute maximum of ( array - offset )
    IF ( mode == 'absoff' )  THEN

!
!--    Determine the local absolut maximum
       fmax_l(1)     = 0.0_wp
       fmax_ijk_l(1) =  i1
       fmax_ijk_l(2) =  j1
       fmax_ijk_l(3) =  k1
       DO  k = k1, k2
          DO  j = j1, j2
!
!--          Attention: the lowest gridpoint is excluded here, because there
!--          ---------  is no advection at nzb=0 and mode 'absoff' is only
!--                     used for calculating u,v extrema for CFL-criteria
             DO  i = i1+1, i2
                IF ( ABS( ar(i,j,k) - offset ) > fmax_l(1) )  THEN
                   fmax_l(1) = ABS( ar(i,j,k) - offset )
                   fmax_ijk_l(1) = i
                   fmax_ijk_l(2) = j
                   fmax_ijk_l(3) = k
                ENDIF
             ENDDO
          ENDDO
       ENDDO

!
!--    Set a flag in case that the determined value is negative.
!--    A constant offset has to be subtracted in order to handle the special
!--    case i=0 correctly
       IF ( ar(fmax_ijk_l(1),fmax_ijk_l(2),fmax_ijk_l(3)) < 0.0_wp )  THEN
          fmax_ijk_l(1) = -fmax_ijk_l(1) - 10
       ENDIF

#if defined( __parallel )
       fmax_l(2)  = myid
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( fmax_l, fmax, 1, MPI_2REAL, MPI_MAXLOC, comm2d, &
                           ierr )

!
!--    Determine the global absolut maximum. Result stored on PE0.
       id_fmax = fmax(2)
       IF ( id_fmax /= 0 )  THEN
          IF ( myid == 0 )  THEN
             CALL MPI_RECV( fmax_ijk, 3, MPI_INTEGER, id_fmax, 0, comm2d, &
                            status, ierr )
          ELSEIF ( myid == id_fmax )  THEN
             CALL MPI_SEND( fmax_ijk_l, 3, MPI_INTEGER, 0, 0, comm2d, ierr )
          ENDIF
       ELSE
          fmax_ijk = fmax_ijk_l
       ENDIF
!
!--    Send the indices of the just determined absolut maximum to other PEs
       CALL MPI_BCAST( fmax_ijk, 3, MPI_INTEGER, 0, comm2d, ierr )
#else
       fmax(1)  = fmax_l(1)
       fmax_ijk = fmax_ijk_l
#endif

    ENDIF

!
!-- Determine output parameters
    SELECT CASE( mode )

       CASE( 'min' )

          value     = fmin(1)
          value_ijk = fmin_ijk

       CASE( 'max' )

          value     = fmax(1)
          value_ijk = fmax_ijk

       CASE( 'minmax' )

          value      = fmin(1)
          value_ijk  = fmin_ijk
          value1     = fmax(1)
          value1_ijk = fmax_ijk

       CASE( 'abs', 'absoff' )

          value     = fmax(1)
          value_ijk = fmax_ijk
          IF ( fmax_ijk(1) < 0 )  THEN
             value        = -value
             value_ijk(1) = -value_ijk(1) - 10         !???
          ENDIF

    END SELECT

!
!-- Limit index values to the range 0..nx, 0..ny
    IF ( value_ijk(3) < 0  ) value_ijk(3) = nx +1 + value_ijk(3)
    IF ( value_ijk(3) > nx ) value_ijk(3) = value_ijk(3) - (nx+1)
    IF ( value_ijk(2) < 0  ) value_ijk(2) = ny +1 + value_ijk(2)
    IF ( value_ijk(2) > ny ) value_ijk(2) = value_ijk(2) - (ny+1)


 END SUBROUTINE global_min_max
