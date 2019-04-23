!> @file transpose.f90
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
! $Id: transpose.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2119 2017-01-17 16:51:50Z raasch
!
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC directives removed
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doygen readable
!
! 1324 2014-03-21 09:13:16Z suehring
! Bugfix: ONLY statement for module pegrid removed
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
! 1318 2014-03-17 13:35:16Z raasch
! cpu_log_nowait parameter added to cpu measurements of the transpositions
! required for solving the Poisson equation (poisfft),
! module interfaces removed
!
! 1257 2013-11-08 15:18:40Z raasch
! openacc loop and loop vector clauses removed
!
! 1216 2013-08-26 09:31:42Z raasch
! re-sorting of the transposed / to be transposed arrays moved to separate
! routines resort_for_...
!
! 1111 2013-03-08 23:54:10Z raasch
! openACC directives added,
! resorting data from/to work changed, work got 4 dimensions instead of 1
!
! 1106 2013-03-04 05:31:38Z raasch
! preprocessor lines rearranged so that routines can also be used in serial
! (non-parallel) mode
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1003 2012-09-14 14:35:53Z raasch
! indices nxa, nya, etc. replaced by nx, ny, etc.
!
! Revision 1.1  1997/07/24 11:25:18  raasch
! Initial revision
!
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorting data for the transposition from x to y. The transposition itself
!> is carried out in transpose_xy
!------------------------------------------------------------------------------!
 SUBROUTINE resort_for_xy( f_in, f_inv )


     USE indices,                                                              &
         ONLY:  nx

     USE kinds

     USE transpose_indices,                                                    &
         ONLY:  nxl_z, nxr_z, nyn_x, nyn_z, nys_x, nys_z, nzb_x, nzt_x

     IMPLICIT NONE

     REAL(wp) ::  f_in(0:nx,nys_x:nyn_x,nzb_x:nzt_x)  !<
     REAL(wp) ::  f_inv(nys_x:nyn_x,nzb_x:nzt_x,0:nx) !<


     INTEGER(iwp) ::  i !<
     INTEGER(iwp) ::  j !<
     INTEGER(iwp) ::  k !<
!
!-- Rearrange indices of input array in order to make data to be send
!-- by MPI contiguous
    !$OMP  PARALLEL PRIVATE ( i, j, k )
    !$OMP  DO
     DO  i = 0, nx
         DO  k = nzb_x, nzt_x
             DO  j = nys_x, nyn_x
                 f_inv(j,k,i) = f_in(i,j,k)
             ENDDO
         ENDDO
     ENDDO
     !$OMP  END PARALLEL

 END SUBROUTINE resort_for_xy


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Transposition of input array (f_in) from x to y. For the input array, all
!> elements along x reside on the same PE, while after transposition, all
!> elements along y reside on the same PE.
!------------------------------------------------------------------------------!
 SUBROUTINE transpose_xy( f_inv, f_out )


    USE cpulog,                                                                &
        ONLY:  cpu_log, cpu_log_nowait, log_point_s

    USE indices,                                                               &
        ONLY:  nx, ny

    USE kinds

    USE pegrid

    USE transpose_indices,                                                     &
        ONLY:  nxl_y, nxr_y, nyn_x, nys_x, nzb_x, nzb_y, nzt_x, nzt_y

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  l  !<
    INTEGER(iwp) ::  ys !<

    REAL(wp) ::  f_inv(nys_x:nyn_x,nzb_x:nzt_x,0:nx) !<
    REAL(wp) ::  f_out(0:ny,nxl_y:nxr_y,nzb_y:nzt_y) !<

    REAL(wp), DIMENSION(nyn_x-nys_x+1,nzb_y:nzt_y,nxl_y:nxr_y,0:pdims(2)-1) ::  work !<


    IF ( numprocs /= 1 )  THEN

#if defined( __parallel )
!
!--    Transpose array
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start', cpu_log_nowait )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLTOALL( f_inv(nys_x,nzb_x,0),  sendrecvcount_xy, MPI_REAL, &
                          work(1,nzb_y,nxl_y,0), sendrecvcount_xy, MPI_REAL, &
                          comm1dy, ierr )
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )

!
!--    Reorder transposed array
!$OMP  PARALLEL PRIVATE ( i, j, k, l, ys )
!$OMP  DO
       DO  l = 0, pdims(2) - 1
          ys = 0 + l * ( nyn_x - nys_x + 1 )
          DO  i = nxl_y, nxr_y
             DO  k = nzb_y, nzt_y
                DO  j = ys, ys + nyn_x - nys_x
                   f_out(j,i,k) = work(j-ys+1,k,i,l)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
!$OMP  END PARALLEL
#endif

    ELSE

!
!--    Reorder transposed array
!$OMP  PARALLEL PRIVATE ( i, j, k )
!$OMP  DO
       DO  k = nzb_y, nzt_y
          DO  i = nxl_y, nxr_y
             DO  j = 0, ny
                f_out(j,i,k) = f_inv(j,k,i)
             ENDDO
          ENDDO
       ENDDO
!$OMP  END PARALLEL

    ENDIF

 END SUBROUTINE transpose_xy


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorting data after the transposition from x to z. The transposition itself
!> is carried out in transpose_xz
!------------------------------------------------------------------------------!
 SUBROUTINE resort_for_xz( f_inv, f_out )


     USE indices,                                                              &
         ONLY:  nxl, nxr, nyn, nys, nz

     USE kinds

     IMPLICIT NONE

     REAL(wp) ::  f_inv(nys:nyn,nxl:nxr,1:nz) !<
     REAL(wp) ::  f_out(1:nz,nys:nyn,nxl:nxr) !<

     INTEGER(iwp) ::  i !<
     INTEGER(iwp) ::  j !<
     INTEGER(iwp) ::  k !<
!
!-- Rearrange indices of input array in order to make data to be send
!-- by MPI contiguous.
!-- In case of parallel fft/transposition, scattered store is faster in
!-- backward direction!!!
    !$OMP  PARALLEL PRIVATE ( i, j, k )
    !$OMP  DO
     DO  k = 1, nz
         DO  i = nxl, nxr
             DO  j = nys, nyn
                 f_out(k,j,i) = f_inv(j,i,k)
             ENDDO
         ENDDO
     ENDDO
     !$OMP  END PARALLEL

 END SUBROUTINE resort_for_xz


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Transposition of input array (f_in) from x to z. For the input array, all
!> elements along x reside on the same PE, while after transposition, all
!> elements along z reside on the same PE.
!------------------------------------------------------------------------------!
 SUBROUTINE transpose_xz( f_in, f_inv )


    USE cpulog,                                                                &
        ONLY:  cpu_log, cpu_log_nowait, log_point_s

    USE indices,                                                               &
        ONLY:  nnx, nx, nxl, nxr, ny, nyn, nys, nz

    USE kinds

    USE pegrid

    USE transpose_indices,                                                     &
        ONLY:  nyn_x, nys_x, nzb_x, nzt_x

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  l  !<
    INTEGER(iwp) ::  xs !<

    REAL(wp) ::  f_in(0:nx,nys_x:nyn_x,nzb_x:nzt_x) !<
    REAL(wp) ::  f_inv(nys:nyn,nxl:nxr,1:nz) !<

    REAL(wp), DIMENSION(nys_x:nyn_x,nnx,nzb_x:nzt_x,0:pdims(1)-1) ::  work !<


!
!-- If the PE grid is one-dimensional along y, the array has only to be
!-- reordered locally and therefore no transposition has to be done.
    IF ( pdims(1) /= 1 )  THEN

#if defined( __parallel )
!
!--    Reorder input array for transposition
!$OMP  PARALLEL PRIVATE ( i, j, k, l, xs )
!$OMP  DO
       DO  l = 0, pdims(1) - 1
          xs = 0 + l * nnx
          DO  k = nzb_x, nzt_x
             DO  i = xs, xs + nnx - 1
                DO  j = nys_x, nyn_x
                   work(j,i-xs+1,k,l) = f_in(i,j,k)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
!$OMP  END PARALLEL

!
!--    Transpose array
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start', cpu_log_nowait )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLTOALL( work(nys_x,1,nzb_x,0), sendrecvcount_zx, MPI_REAL, &
                          f_inv(nys,nxl,1),      sendrecvcount_zx, MPI_REAL, &
                          comm1dx, ierr )
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )
#endif

    ELSE

!
!--    Reorder the array in a way that the z index is in first position
!$OMP  PARALLEL PRIVATE ( i, j, k )
!$OMP  DO
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = 1, nz
                f_inv(j,i,k) = f_in(i,j,k)
             ENDDO
          ENDDO
       ENDDO
!$OMP  END PARALLEL

    ENDIF

 END SUBROUTINE transpose_xz


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorting data after the transposition from y to x. The transposition itself
!> is carried out in transpose_yx
!------------------------------------------------------------------------------!
 SUBROUTINE resort_for_yx( f_inv, f_out )


     USE indices,                                                              &
         ONLY:  nx

     USE kinds

     USE transpose_indices,                                                    &
         ONLY:  nyn_x, nys_x, nzb_x, nzt_x

     IMPLICIT NONE

     REAL(wp) ::  f_inv(nys_x:nyn_x,nzb_x:nzt_x,0:nx) !<
     REAL(wp) ::  f_out(0:nx,nys_x:nyn_x,nzb_x:nzt_x) !<


     INTEGER(iwp) ::  i !<
     INTEGER(iwp) ::  j !<
     INTEGER(iwp) ::  k !<
!
!-- Rearrange indices of input array in order to make data to be send
!-- by MPI contiguous
    !$OMP  PARALLEL PRIVATE ( i, j, k )
    !$OMP  DO
     DO  i = 0, nx
         DO  k = nzb_x, nzt_x
             DO  j = nys_x, nyn_x
                 f_out(i,j,k) = f_inv(j,k,i)
             ENDDO
         ENDDO
     ENDDO
     !$OMP  END PARALLEL

 END SUBROUTINE resort_for_yx


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Transposition of input array (f_in) from y to x. For the input array, all
!> elements along y reside on the same PE, while after transposition, all
!> elements along x reside on the same PE.
!------------------------------------------------------------------------------!
 SUBROUTINE transpose_yx( f_in, f_inv )


    USE cpulog,                                                                &
        ONLY:  cpu_log, cpu_log_nowait, log_point_s

    USE indices,                                                               &
        ONLY:  nx, ny

    USE kinds

    USE pegrid

    USE transpose_indices,                                                     &
        ONLY:  nxl_y, nxr_y, nyn_x, nys_x, nzb_x, nzb_y, nzt_x, nzt_y

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  l  !<
    INTEGER(iwp) ::  ys !<

    REAL(wp) ::  f_in(0:ny,nxl_y:nxr_y,nzb_y:nzt_y)  !<
    REAL(wp) ::  f_inv(nys_x:nyn_x,nzb_x:nzt_x,0:nx) !<

    REAL(wp), DIMENSION(nyn_x-nys_x+1,nzb_y:nzt_y,nxl_y:nxr_y,0:pdims(2)-1) ::  work !<


    IF ( numprocs /= 1 )  THEN

#if defined( __parallel )
!
!--    Reorder input array for transposition
!$OMP  PARALLEL PRIVATE ( i, j, k, l, ys )
!$OMP  DO
       DO  l = 0, pdims(2) - 1
          ys = 0 + l * ( nyn_x - nys_x + 1 )
          DO  i = nxl_y, nxr_y
             DO  k = nzb_y, nzt_y
                DO  j = ys, ys + nyn_x - nys_x
                   work(j-ys+1,k,i,l) = f_in(j,i,k)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
!$OMP  END PARALLEL

!
!--    Transpose array
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start', cpu_log_nowait )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLTOALL( work(1,nzb_y,nxl_y,0), sendrecvcount_xy, MPI_REAL, &
                          f_inv(nys_x,nzb_x,0),  sendrecvcount_xy, MPI_REAL, &
                          comm1dy, ierr )
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )
#endif

    ELSE

!
!--    Reorder array f_in the same way as ALLTOALL did it
!$OMP  PARALLEL PRIVATE ( i, j, k )
!$OMP  DO
       DO  i = nxl_y, nxr_y
          DO  k = nzb_y, nzt_y
             DO  j = 0, ny
                f_inv(j,k,i) = f_in(j,i,k)
             ENDDO
          ENDDO
       ENDDO
!$OMP  END PARALLEL

    ENDIF

 END SUBROUTINE transpose_yx


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Transposition of input array (f_in) from y to x. For the input array, all
!> elements along y reside on the same PE, while after transposition, all
!> elements along x reside on the same PE.
!> This is a direct transposition for arrays with indices in regular order
!> (k,j,i) (cf. transpose_yx).
!------------------------------------------------------------------------------!
 SUBROUTINE transpose_yxd( f_in, f_out )


    USE cpulog,                                                                &
        ONLY:  cpu_log, cpu_log_nowait, log_point_s

    USE indices,                                                               &
        ONLY:  nnx, nny, nnz, nx, nxl, nxr, nyn, nys, nz

    USE kinds

    USE pegrid

    USE transpose_indices,                                                     &
        ONLY:  nyn_x, nys_x, nzb_x, nzt_x

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  l  !<
    INTEGER(iwp) ::  m  !<
    INTEGER(iwp) ::  xs !<

    REAL(wp) ::  f_in(1:nz,nys:nyn,nxl:nxr)          !<
    REAL(wp) ::  f_inv(nxl:nxr,1:nz,nys:nyn)         !<
    REAL(wp) ::  f_out(0:nx,nys_x:nyn_x,nzb_x:nzt_x) !<
    REAL(wp) ::  work(nnx*nny*nnz)                   !<
#if defined( __parallel )

!
!-- Rearrange indices of input array in order to make data to be send
!-- by MPI contiguous
    DO  k = 1, nz
       DO  j = nys, nyn
          DO  i = nxl, nxr
             f_inv(i,k,j) = f_in(k,j,i)
          ENDDO
       ENDDO
    ENDDO

!
!-- Transpose array
    CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start' )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLTOALL( f_inv(nxl,1,nys), sendrecvcount_xy, MPI_REAL, &
                       work(1),          sendrecvcount_xy, MPI_REAL, &
                       comm1dx, ierr )
    CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )

!
!-- Reorder transposed array
    m = 0
    DO  l = 0, pdims(1) - 1
       xs = 0 + l * nnx
       DO  j = nys_x, nyn_x
          DO  k = 1, nz
             DO  i = xs, xs + nnx - 1
                m = m + 1
                f_out(i,j,k) = work(m)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

#endif

 END SUBROUTINE transpose_yxd


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorting data for the transposition from y to z. The transposition itself
!> is carried out in transpose_yz
!------------------------------------------------------------------------------!
 SUBROUTINE resort_for_yz( f_in, f_inv )


     USE indices,                                                              &
         ONLY:  ny

     USE kinds

     USE transpose_indices,                                                    &
         ONLY:  nxl_y, nxr_y, nzb_y, nzt_y

     IMPLICIT NONE

     REAL(wp) ::  f_in(0:ny,nxl_y:nxr_y,nzb_y:nzt_y)  !<
     REAL(wp) ::  f_inv(nxl_y:nxr_y,nzb_y:nzt_y,0:ny) !<

     INTEGER(iwp) ::  i !<
     INTEGER(iwp) ::  j !<
     INTEGER(iwp) ::  k !<

!
!-- Rearrange indices of input array in order to make data to be send
!-- by MPI contiguous
    !$OMP  PARALLEL PRIVATE ( i, j, k )
    !$OMP  DO
     DO  j = 0, ny
         DO  k = nzb_y, nzt_y
             DO  i = nxl_y, nxr_y
                 f_inv(i,k,j) = f_in(j,i,k)
             ENDDO
         ENDDO
     ENDDO
     !$OMP  END PARALLEL

 END SUBROUTINE resort_for_yz


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Transposition of input array (f_in) from y to z. For the input array, all
!> elements along y reside on the same PE, while after transposition, all
!> elements along z reside on the same PE.
!------------------------------------------------------------------------------!
 SUBROUTINE transpose_yz( f_inv, f_out )


    USE cpulog,                                                                &
        ONLY:  cpu_log, cpu_log_nowait, log_point_s

    USE indices,                                                               &
        ONLY:  ny, nz

    USE kinds

    USE pegrid

    USE transpose_indices,                                                     &
        ONLY:  nxl_y, nxl_z, nxr_y, nxr_z, nyn_z, nys_z, nzb_y, nzt_y

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  l  !<
    INTEGER(iwp) ::  zs !<

    REAL(wp) ::  f_inv(nxl_y:nxr_y,nzb_y:nzt_y,0:ny) !<
    REAL(wp) ::  f_out(nxl_z:nxr_z,nys_z:nyn_z,1:nz) !<

    REAL(wp), DIMENSION(nxl_z:nxr_z,nzt_y-nzb_y+1,nys_z:nyn_z,0:pdims(1)-1) ::  work !<


!
!-- If the PE grid is one-dimensional along y, only local reordering
!-- of the data is necessary and no transposition has to be done.
    IF ( pdims(1) == 1 )  THEN

!$OMP  PARALLEL PRIVATE ( i, j, k )
!$OMP  DO
       DO  j = 0, ny
          DO  k = nzb_y, nzt_y
             DO  i = nxl_y, nxr_y
                f_out(i,j,k) = f_inv(i,k,j)
             ENDDO
          ENDDO
       ENDDO
!$OMP  END PARALLEL

    ELSE

#if defined( __parallel )
!
!--    Transpose array
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start', cpu_log_nowait )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLTOALL( f_inv(nxl_y,nzb_y,0),  sendrecvcount_yz, MPI_REAL, &
                          work(nxl_z,1,nys_z,0), sendrecvcount_yz, MPI_REAL, &
                          comm1dx, ierr )
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )

!
!--    Reorder transposed array
!$OMP  PARALLEL PRIVATE ( i, j, k, l, zs )
!$OMP  DO
       DO  l = 0, pdims(1) - 1
          zs = 1 + l * ( nzt_y - nzb_y + 1 )
          DO  j = nys_z, nyn_z
             DO  k = zs, zs + nzt_y - nzb_y
                DO  i = nxl_z, nxr_z
                   f_out(i,j,k) = work(i,k-zs+1,j,l)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
!$OMP  END PARALLEL
#endif

   ENDIF

 END SUBROUTINE transpose_yz


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorting data for the transposition from z to x. The transposition itself
!> is carried out in transpose_zx
!------------------------------------------------------------------------------!
 SUBROUTINE resort_for_zx( f_in, f_inv )


     USE indices,                                                              &
         ONLY:  nxl, nxr, nyn, nys, nz

     USE kinds

     IMPLICIT NONE

     REAL(wp) ::  f_in(1:nz,nys:nyn,nxl:nxr)  !<
     REAL(wp) ::  f_inv(nys:nyn,nxl:nxr,1:nz) !<

     INTEGER(iwp) ::  i !<
     INTEGER(iwp) ::  j !<
     INTEGER(iwp) ::  k !<

!
!-- Rearrange indices of input array in order to make data to be send
!-- by MPI contiguous
    !$OMP  PARALLEL PRIVATE ( i, j, k )
    !$OMP  DO
     DO  k = 1,nz
         DO  i = nxl, nxr
             DO  j = nys, nyn
                 f_inv(j,i,k) = f_in(k,j,i)
             ENDDO
         ENDDO
     ENDDO
     !$OMP  END PARALLEL

 END SUBROUTINE resort_for_zx


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Transposition of input array (f_in) from z to x. For the input array, all
!> elements along z reside on the same PE, while after transposition, all
!> elements along x reside on the same PE.
!------------------------------------------------------------------------------!
 SUBROUTINE transpose_zx( f_inv, f_out )


    USE cpulog,                                                                &
        ONLY:  cpu_log, cpu_log_nowait, log_point_s

    USE indices,                                                               &
        ONLY:  nnx, nx, nxl, nxr, nyn, nys, nz

    USE kinds

    USE pegrid

    USE transpose_indices,                                                     &
        ONLY:  nyn_x, nys_x, nzb_x, nzt_x

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  l  !<
    INTEGER(iwp) ::  xs !<

    REAL(wp) ::  f_inv(nys:nyn,nxl:nxr,1:nz)         !<
    REAL(wp) ::  f_out(0:nx,nys_x:nyn_x,nzb_x:nzt_x) !<

    REAL(wp), DIMENSION(nys_x:nyn_x,nnx,nzb_x:nzt_x,0:pdims(1)-1) ::  work !<


!
!-- If the PE grid is one-dimensional along y, only local reordering
!-- of the data is necessary and no transposition has to be done.
    IF ( pdims(1) == 1 )  THEN

!$OMP  PARALLEL PRIVATE ( i, j, k )
!$OMP  DO
       DO  k = 1, nz
          DO  i = nxl, nxr
             DO  j = nys, nyn
                f_out(i,j,k) = f_inv(j,i,k)
             ENDDO
          ENDDO
       ENDDO
!$OMP  END PARALLEL

    ELSE

#if defined( __parallel )
!
!--    Transpose array
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start', cpu_log_nowait )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLTOALL( f_inv(nys,nxl,1),      sendrecvcount_zx, MPI_REAL, &
                          work(nys_x,1,nzb_x,0), sendrecvcount_zx, MPI_REAL, &
                          comm1dx, ierr )
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )

!
!--    Reorder transposed array
!$OMP  PARALLEL PRIVATE ( i, j, k, l, xs )
!$OMP  DO
       DO  l = 0, pdims(1) - 1
          xs = 0 + l * nnx
          DO  k = nzb_x, nzt_x
             DO  i = xs, xs + nnx - 1
                DO  j = nys_x, nyn_x
                   f_out(i,j,k) = work(j,i-xs+1,k,l)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
!$OMP  END PARALLEL
#endif

    ENDIF

 END SUBROUTINE transpose_zx


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorting data after the transposition from z to y. The transposition itself
!> is carried out in transpose_zy
!------------------------------------------------------------------------------!
 SUBROUTINE resort_for_zy( f_inv, f_out )


     USE indices,                                                              &
         ONLY:  ny

     USE kinds

     USE transpose_indices,                                                    &
         ONLY:  nxl_y, nxr_y, nzb_y, nzt_y

     IMPLICIT NONE

     REAL(wp) ::  f_inv(nxl_y:nxr_y,nzb_y:nzt_y,0:ny) !<
     REAL(wp) ::  f_out(0:ny,nxl_y:nxr_y,nzb_y:nzt_y) !<


     INTEGER(iwp) ::  i !<
     INTEGER(iwp) ::  j !<
     INTEGER(iwp) ::  k !<

!
!-- Rearrange indices of input array in order to make data to be send
!-- by MPI contiguous
    !$OMP  PARALLEL PRIVATE ( i, j, k )
    !$OMP  DO
     DO  k = nzb_y, nzt_y
         DO  j = 0, ny
             DO  i = nxl_y, nxr_y
                 f_out(j,i,k) = f_inv(i,k,j)
             ENDDO
         ENDDO
     ENDDO
     !$OMP  END PARALLEL

 END SUBROUTINE resort_for_zy


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Transposition of input array (f_in) from z to y. For the input array, all
!> elements along z reside on the same PE, while after transposition, all
!> elements along y reside on the same PE.
!------------------------------------------------------------------------------!
 SUBROUTINE transpose_zy( f_in, f_inv )


    USE cpulog,                                                                &
        ONLY:  cpu_log, cpu_log_nowait, log_point_s

    USE indices,                                                               &
        ONLY:  ny, nz

    USE kinds

    USE pegrid

    USE transpose_indices,                                                     &
        ONLY:  nxl_y, nxl_z, nxr_y, nxr_z, nyn_z, nys_z, nzb_y, nzt_y

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  l  !<
    INTEGER(iwp) ::  zs !<

    REAL(wp) ::  f_in(nxl_z:nxr_z,nys_z:nyn_z,1:nz)  !<
    REAL(wp) ::  f_inv(nxl_y:nxr_y,nzb_y:nzt_y,0:ny) !<

    REAL(wp), DIMENSION(nxl_z:nxr_z,nzt_y-nzb_y+1,nys_z:nyn_z,0:pdims(1)-1) ::  work !<

!
!-- If the PE grid is one-dimensional along y, the array has only to be
!-- reordered locally and therefore no transposition has to be done.
    IF ( pdims(1) /= 1 )  THEN

#if defined( __parallel )
!
!--    Reorder input array for transposition
!$OMP  PARALLEL PRIVATE ( i, j, k, l, zs )
!$OMP  DO
       DO  l = 0, pdims(1) - 1
          zs = 1 + l * ( nzt_y - nzb_y + 1 )
          DO  j = nys_z, nyn_z
             DO  k = zs, zs + nzt_y - nzb_y
                DO  i = nxl_z, nxr_z
                   work(i,k-zs+1,j,l) = f_in(i,j,k)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
!$OMP  END PARALLEL

!
!--    Transpose array
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start', cpu_log_nowait )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLTOALL( work(nxl_z,1,nys_z,0), sendrecvcount_yz, MPI_REAL, &
                          f_inv(nxl_y,nzb_y,0),  sendrecvcount_yz, MPI_REAL, &
                          comm1dx, ierr )
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )
#endif

    ELSE
!
!--    Reorder the array in the same way like ALLTOALL did it
!$OMP  PARALLEL PRIVATE ( i, j, k )
!$OMP  DO
       DO  k = nzb_y, nzt_y
          DO  j = 0, ny
             DO  i = nxl_y, nxr_y
                f_inv(i,k,j) = f_in(i,j,k)
             ENDDO
          ENDDO
       ENDDO
!$OMP  END PARALLEL

    ENDIF

 END SUBROUTINE transpose_zy


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Transposition of input array (f_in) from z to y. For the input array, all
!> elements along z reside on the same PE, while after transposition, all
!> elements along y reside on the same PE.
!> This is a direct transposition for arrays with indices in regular order
!> (k,j,i) (cf. transpose_zy).
!------------------------------------------------------------------------------!
 SUBROUTINE transpose_zyd( f_in, f_out )


    USE cpulog,                                                                &
        ONLY:  cpu_log, cpu_log_nowait, log_point_s

    USE indices,                                                               &
        ONLY:  nnx, nny, nnz, nxl, nxr, nyn, nys, ny, nz

    USE kinds

    USE pegrid

    USE transpose_indices,                                                     &
        ONLY:  nxl_y, nxl_yd, nxr_y, nxr_yd, nzb_y, nzb_yd, nzt_y, nzt_yd

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  l  !<
    INTEGER(iwp) ::  m  !<
    INTEGER(iwp) ::  ys !<

    REAL(wp) ::  f_in(1:nz,nys:nyn,nxl:nxr)              !<
    REAL(wp) ::  f_inv(nys:nyn,nxl:nxr,1:nz)             !<
    REAL(wp) ::  f_out(0:ny,nxl_yd:nxr_yd,nzb_yd:nzt_yd) !<
    REAL(wp) ::  work(nnx*nny*nnz)                       !<

#if defined( __parallel )

!
!-- Rearrange indices of input array in order to make data to be send
!-- by MPI contiguous
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = 1, nz
             f_inv(j,i,k) = f_in(k,j,i)
          ENDDO
       ENDDO
    ENDDO

!
!-- Move data to different array, because memory location of work1 is
!-- needed further below (work1 = work2).
!-- If the PE grid is one-dimensional along x, only local reordering
!-- of the data is necessary and no transposition has to be done.
    IF ( pdims(2) == 1 )  THEN
       DO  k = 1, nz
          DO  i = nxl, nxr
             DO  j = nys, nyn
                f_out(j,i,k) = f_inv(j,i,k)
             ENDDO
          ENDDO
       ENDDO
       RETURN
    ENDIF

!
!-- Transpose array
    CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start' )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLTOALL( f_inv(nys,nxl,1), sendrecvcount_zyd, MPI_REAL, &
                       work(1),          sendrecvcount_zyd, MPI_REAL, &
                       comm1dy, ierr )
    CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )

!
!-- Reorder transposed array
    m = 0
    DO  l = 0, pdims(2) - 1
       ys = 0 + l * nny
       DO  k = nzb_yd, nzt_yd
          DO  i = nxl_yd, nxr_yd
             DO  j = ys, ys + nny - 1
                m = m + 1
                f_out(j,i,k) = work(m)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

#endif

 END SUBROUTINE transpose_zyd
