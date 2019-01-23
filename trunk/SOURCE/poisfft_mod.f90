!> @file poisfft_mod.f90
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
! $Id: poisfft_mod.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2300 2017-06-29 13:31:14Z raasch
! settings depending on host variable removed or replaced by loop_optimization
! 
! 2119 2017-01-17 16:51:50Z raasch
!
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC directives and related code removed
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
! 
! 
! 1804 2016-04-05 16:30:18Z maronga
! Removed code for parameter file check (__check)
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
!
! 1482 2014-10-18 12:34:45Z raasch
! use 2d-decomposition, if accelerator boards are used
!
! 1406 2014-05-16 13:47:01Z raasch
! bugfix for pgi 14.4: declare create moved after array declaration
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
! module interfaces removed
!
! 1306 2014-03-13 14:30:59Z raasch
! openmp sections removed from the overlap branch,
! second argument removed from parameter list
!
! 1216 2013-08-26 09:31:42Z raasch
! resorting of arrays moved to separate routines resort_for_...,
! one argument, used as temporary work array, removed from all transpose
! routines
! overlapping fft / transposition implemented
!
! 1212 2013-08-15 08:46:27Z raasch
! tridia routines moved to seperate module tridia_solver
!
! 1208 2013-08-13 06:41:49Z raasch
! acc-update clauses added for "ar" so that ffts other than cufft can also be
! used (although they are not ported and will give a poor performance)
!
! 1111 2013-03-08 23:54:10Z raasch
! further openACC porting of non-parallel (MPI) branch:
! tridiagonal routines split into extermal subroutines (instead using CONTAINS),
! no distinction between parallel/non-parallel in poisfft and tridia any more,
! tridia routines moved to end of file because of probable bug in PGI compiler 12.5
! (otherwise "invalid device function" is indicated during runtime),
! optimization of tridia routines: constant elements and coefficients of tri are
! stored in seperate arrays ddzuw and tric, last dimension of tri reduced from 5
! to 2,
! poisfft_init is now called internally from poisfft, maketri is called from
! poisfft_init,
! ibc_p_b = 2 removed
!
! 1106 2013-03-04 05:31:38Z raasch
! routines fftx, ffty, fftxp, fftyp removed, calls replaced by fft_x, fft_y,
! in the 1D-decomposition routines fft_x, ffty are replaced by fft_x_1d,
! fft_y_1d
!
! 1103 2013-02-20 02:15:53Z raasch
! tri, ar, and ar1 arguments in tridia-routines (2d) are removed because they
! sometimes cause segmentation faults with intel 12.1 compiler
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 2012-09-21 07:03:55Z raasch
! FLOAT type conversion replaced by REAL
!
! 1003 2012-09-14 14:35:53Z raasch
! indices nxa, nya, etc. replaced by nx, ny, etc.
!
! 940 2012-07-09 14:31:00Z raasch
! special handling of tri-array as an argument in tridia_1dd routines switched
! off because it caused segmentation faults with intel 12.1 compiler
!
! 877 2012-04-03 11:21:44Z suehring
! Bugfix: Avoid divisions by zero in case of using a 'neumann' bc for the
! pressure at the top of the model domain.
!
! 809 2012-01-30 13:32:58Z maronga
! Bugfix: replaced .AND. and .NOT. with && and ! in the preprocessor directives
!
! 807 2012-01-25 11:53:51Z maronga
! New cpp directive "__check" implemented which is used by check_namelist_files
! (most of the code is unneeded by check_namelist_files).
!
! Revision 1.1  1997/07/24 11:24:14  raasch
! Initial revision
!
!
! Description:
! ------------
!> Solves the Poisson equation with a 2D spectral method
!>        d^2 p / dx^2 + d^2 p / dy^2 + d^2 p / dz^2 = s
!>
!> Input:
!> real    ar   contains (nnz,nny,nnx) elements of the velocity divergence,
!>              starting from (1,nys,nxl)
!>
!> Output:
!> real    ar   contains the solution for perturbation pressure p
!------------------------------------------------------------------------------!
 MODULE poisfft_mod
 

    USE fft_xy,                                                                &
        ONLY:  fft_init, fft_y, fft_y_1d, fft_y_m, fft_x, fft_x_1d, fft_x_m

    USE indices,                                                               &
        ONLY:  nnx, nny, nx, nxl, nxr, ny, nys, nyn, nz

    USE transpose_indices,                                                     &
        ONLY:  nxl_y, nxl_z, nxr_y, nxr_z, nys_x, nys_z, nyn_x, nyn_z, nzb_x,  &
               nzb_y, nzt_x, nzt_y

    USE tridia_solver,                                                         &
        ONLY:  tridia_1dd, tridia_init, tridia_substi, tridia_substi_overlap

    IMPLICIT NONE

    PRIVATE

    PUBLIC  poisfft, poisfft_init

    INTERFACE poisfft
       MODULE PROCEDURE poisfft
    END INTERFACE poisfft

    INTERFACE poisfft_init
       MODULE PROCEDURE poisfft_init
    END INTERFACE poisfft_init


 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE poisfft_init

       USE arrays_3d,                                                          &
           ONLY:  ddzu_pres, ddzw

       USE kinds

       USE control_parameters, only: poisfft_initialized

       IMPLICIT NONE

       INTEGER(iwp) ::  k  !<


       CALL fft_init

       CALL tridia_init

       poisfft_initialized = .TRUE.

    END SUBROUTINE poisfft_init



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Two-dimensional Fourier Transformation in x- and y-direction.
!------------------------------------------------------------------------------!
    SUBROUTINE poisfft( ar )

       USE control_parameters,                                                 &
           ONLY:  fft_method, transpose_compute_overlap, poisfft_initialized

       USE cpulog,                                                             &
           ONLY:  cpu_log, cpu_log_nowait, log_point_s

       USE kinds

       USE pegrid

       IMPLICIT NONE

       INTEGER(iwp) ::  ii           !<
       INTEGER(iwp) ::  iind         !<
       INTEGER(iwp) ::  inew         !<
       INTEGER(iwp) ::  jj           !<
       INTEGER(iwp) ::  jind         !<
       INTEGER(iwp) ::  jnew         !<
       INTEGER(iwp) ::  ki           !<
       INTEGER(iwp) ::  kk           !<
       INTEGER(iwp) ::  knew         !<
       INTEGER(iwp) ::  n            !<
       INTEGER(iwp) ::  nblk         !<
       INTEGER(iwp) ::  nnx_y        !<
       INTEGER(iwp) ::  nny_z        !<
       INTEGER(iwp) ::  nnz_t        !<
       INTEGER(iwp) ::  nnz_x        !<
       INTEGER(iwp) ::  nxl_y_bound  !<
       INTEGER(iwp) ::  nxr_y_bound  !<

       INTEGER(iwp), DIMENSION(4) ::  isave  !<

       REAL(wp), DIMENSION(1:nz,nys:nyn,nxl:nxr) ::  ar      !<
       REAL(wp), DIMENSION(nys:nyn,nxl:nxr,1:nz) ::  ar_inv  !<

       REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  ar1      !<
       REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  f_in     !<
       REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  f_inv    !<
       REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  f_out_y  !<
       REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  f_out_z  !<

       CALL cpu_log( log_point_s(3), 'poisfft', 'start' )

       IF ( .NOT. poisfft_initialized )  CALL poisfft_init

!
!--    Two-dimensional Fourier Transformation in x- and y-direction.
       IF ( pdims(2) == 1  .AND.  pdims(1) > 1 )  THEN

!
!--       1d-domain-decomposition along x:
!--       FFT along y and transposition y --> x
          CALL ffty_tr_yx( ar, ar )

!
!--       FFT along x, solving the tridiagonal system and backward FFT
          CALL fftx_tri_fftx( ar )

!
!--       Transposition x --> y and backward FFT along y
          CALL tr_xy_ffty( ar, ar )

       ELSEIF ( pdims(1) == 1 .AND. pdims(2) > 1 )  THEN

!
!--       1d-domain-decomposition along y:
!--       FFT along x and transposition x --> y
          CALL fftx_tr_xy( ar, ar )

!
!--       FFT along y, solving the tridiagonal system and backward FFT
          CALL ffty_tri_ffty( ar )

!
!--       Transposition y --> x and backward FFT along x
          CALL tr_yx_fftx( ar, ar )

       ELSEIF ( .NOT. transpose_compute_overlap )  THEN

!
!--       2d-domain-decomposition or no decomposition (1 PE run)
!--       Transposition z --> x
          CALL cpu_log( log_point_s(5), 'transpo forward', 'start' )
          CALL resort_for_zx( ar, ar_inv )
          CALL transpose_zx( ar_inv, ar )
          CALL cpu_log( log_point_s(5), 'transpo forward', 'pause' )

          CALL cpu_log( log_point_s(4), 'fft_x', 'start' )
          CALL fft_x( ar, 'forward' )
          CALL cpu_log( log_point_s(4), 'fft_x', 'pause' )

!
!--       Transposition x --> y
          CALL cpu_log( log_point_s(5), 'transpo forward', 'continue' )
          CALL resort_for_xy( ar, ar_inv )
          CALL transpose_xy( ar_inv, ar )
          CALL cpu_log( log_point_s(5), 'transpo forward', 'pause' )

          CALL cpu_log( log_point_s(7), 'fft_y', 'start' )
          CALL fft_y( ar, 'forward', ar_tr = ar,                &
                      nxl_y_bound = nxl_y, nxr_y_bound = nxr_y, &
                      nxl_y_l = nxl_y, nxr_y_l = nxr_y )
          CALL cpu_log( log_point_s(7), 'fft_y', 'pause' )

!
!--       Transposition y --> z
          CALL cpu_log( log_point_s(5), 'transpo forward', 'continue' )
          CALL resort_for_yz( ar, ar_inv )
          CALL transpose_yz( ar_inv, ar )
          CALL cpu_log( log_point_s(5), 'transpo forward', 'stop' )

!
!--       Solve the tridiagonal equation system along z
          CALL cpu_log( log_point_s(6), 'tridia', 'start' )
          CALL tridia_substi( ar )
          CALL cpu_log( log_point_s(6), 'tridia', 'stop' )

!
!--       Inverse Fourier Transformation
!--       Transposition z --> y
          CALL cpu_log( log_point_s(8), 'transpo invers', 'start' )
          CALL transpose_zy( ar, ar_inv )
          CALL resort_for_zy( ar_inv, ar )
          CALL cpu_log( log_point_s(8), 'transpo invers', 'pause' )

          CALL cpu_log( log_point_s(7), 'fft_y', 'continue' )
          CALL fft_y( ar, 'backward', ar_tr = ar,               &
                      nxl_y_bound = nxl_y, nxr_y_bound = nxr_y, &
                      nxl_y_l = nxl_y, nxr_y_l = nxr_y )
          CALL cpu_log( log_point_s(7), 'fft_y', 'stop' )

!
!--       Transposition y --> x
          CALL cpu_log( log_point_s(8), 'transpo invers', 'continue' )
          CALL transpose_yx( ar, ar_inv )
          CALL resort_for_yx( ar_inv, ar )
          CALL cpu_log( log_point_s(8), 'transpo invers', 'pause' )

          CALL cpu_log( log_point_s(4), 'fft_x', 'continue' )
          CALL fft_x( ar, 'backward' )
          CALL cpu_log( log_point_s(4), 'fft_x', 'stop' )

!
!--       Transposition x --> z
          CALL cpu_log( log_point_s(8), 'transpo invers', 'continue' )
          CALL transpose_xz( ar, ar_inv )
          CALL resort_for_xz( ar_inv, ar )
          CALL cpu_log( log_point_s(8), 'transpo invers', 'stop' )

       ELSE

!
!--       2d-domain-decomposition or no decomposition (1 PE run) with
!--       overlapping transposition / fft
!--       cputime logging must not use barriers, which would prevent overlapping
          ALLOCATE( f_out_y(0:ny,nxl_y:nxr_y,nzb_y:nzt_y), &
                    f_out_z(0:nx,nys_x:nyn_x,nzb_x:nzt_x) )
!
!--       Transposition z --> x + subsequent fft along x
          ALLOCATE( f_inv(nys:nyn,nxl:nxr,1:nz) )
          CALL resort_for_zx( ar, f_inv )
!
!--       Save original indices and gridpoint counter
          isave(1) = nz
          isave(2) = nzb_x
          isave(3) = nzt_x
          isave(4) = sendrecvcount_zx
!
!--       Set new indices for transformation
          nblk  = nz / pdims(1)
          nz    = pdims(1)
          nnz_x = 1
          nzb_x = 1 + myidx * nnz_x
          nzt_x = ( myidx + 1 ) * nnz_x
          sendrecvcount_zx = nnx * nny * nnz_x

          ALLOCATE( ar1(0:nx,nys_x:nyn_x,nzb_x:nzt_x) )
          ALLOCATE( f_in(nys:nyn,nxl:nxr,1:nz) )

          DO  kk = 1, nblk

             IF ( kk == 1 )  THEN
                CALL cpu_log( log_point_s(5), 'transpo forward', 'start', cpu_log_nowait )
             ELSE
                CALL cpu_log( log_point_s(5), 'transpo forward', 'continue', cpu_log_nowait )
             ENDIF

             DO  knew = 1, nz
                ki = kk + nblk * ( knew - 1 )
                f_in(:,:,knew) = f_inv(:,:,ki)
             ENDDO

             CALL transpose_zx( f_in, ar1(:,:,:))
             CALL cpu_log( log_point_s(5), 'transpo forward', 'pause' )

             IF ( kk == 1 )  THEN
                CALL cpu_log( log_point_s(4), 'fft_x', 'start', cpu_log_nowait )
             ELSE
                CALL cpu_log( log_point_s(4), 'fft_x', 'continue', cpu_log_nowait )
             ENDIF

             n = isave(2) + kk - 1
             CALL fft_x( ar1(:,:,:), 'forward',  ar_2d = f_out_z(:,:,n))
             CALL cpu_log( log_point_s(4), 'fft_x', 'pause' )

          ENDDO
!
!--       Restore original indices/counters
          nz               = isave(1)
          nzb_x            = isave(2)
          nzt_x            = isave(3)
          sendrecvcount_zx = isave(4)

          DEALLOCATE( ar1, f_in, f_inv )

!
!--       Transposition x --> y + subsequent fft along y
          ALLOCATE( f_inv(nys_x:nyn_x,nzb_x:nzt_x,0:nx) )
          CALL resort_for_xy( f_out_z, f_inv )
!
!--       Save original indices and gridpoint counter
          isave(1) = nx
          isave(2) = nxl_y
          isave(3) = nxr_y
          isave(4) = sendrecvcount_xy
!
!--       Set new indices for transformation
          nblk  = ( ( nx+1 ) / pdims(2) ) - 1
          nx    = pdims(2)
          nnx_y = 1
          nxl_y = myidy * nnx_y
          nxr_y = ( myidy + 1 ) * nnx_y - 1
          sendrecvcount_xy = nnx_y * ( nyn_x-nys_x+1 ) * ( nzt_x-nzb_x+1 )

          ALLOCATE( ar1(0:ny,nxl_y:nxr_y,nzb_y:nzt_y) )
          ALLOCATE( f_in(nys_x:nyn_x,nzb_x:nzt_x,0:nx) )

          DO  ii = 0, nblk

             CALL cpu_log( log_point_s(5), 'transpo forward', 'continue', cpu_log_nowait )

             DO  inew = 0, nx-1
                iind = ii + ( nblk + 1 ) * inew
                f_in(:,:,inew) = f_inv(:,:,iind)
             ENDDO

             CALL transpose_xy( f_in, ar1(:,:,:) )

             CALL cpu_log( log_point_s(5), 'transpo forward', 'pause' )

             IF ( ii == 1 )  THEN
                CALL cpu_log( log_point_s(7), 'fft_y', 'start', cpu_log_nowait )
             ELSE
                CALL cpu_log( log_point_s(7), 'fft_y', 'continue', cpu_log_nowait )
             ENDIF

             nxl_y_bound = isave(2)
             nxr_y_bound = isave(3)
             n           = isave(2) + ii
             CALL fft_y( ar1(:,:,:), 'forward', ar_tr = f_out_y,               &
                         nxl_y_bound = nxl_y_bound, nxr_y_bound = nxr_y_bound, &
                         nxl_y_l = n, nxr_y_l = n )

             CALL cpu_log( log_point_s(7), 'fft_y', 'pause' )

          ENDDO
!
!--       Restore original indices/counters
          nx               = isave(1)
          nxl_y            = isave(2)
          nxr_y            = isave(3)
          sendrecvcount_xy = isave(4)

          DEALLOCATE( ar1, f_in, f_inv )

!
!--       Transposition y --> z + subsequent tridia + resort for z --> y
          ALLOCATE( f_inv(nxl_y:nxr_y,nzb_y:nzt_y,0:ny) )
          CALL resort_for_yz( f_out_y, f_inv )
!
!--       Save original indices and gridpoint counter
          isave(1) = ny
          isave(2) = nys_z
          isave(3) = nyn_z
          isave(4) = sendrecvcount_yz
!
!--       Set new indices for transformation
          nblk             = ( ( ny+1 ) / pdims(1) ) - 1
          ny               = pdims(1)
          nny_z            = 1
          nys_z            = myidx * nny_z
          nyn_z            = ( myidx + 1 ) * nny_z - 1
          sendrecvcount_yz = ( nxr_y-nxl_y+1 ) * nny_z * ( nzt_y-nzb_y+1 )

          ALLOCATE( ar1(nxl_z:nxr_z,nys_z:nyn_z,1:nz) )
          ALLOCATE( f_in(nxl_y:nxr_y,nzb_y:nzt_y,0:ny) )

          DO  jj = 0, nblk
!
!--          Forward Fourier Transformation
!--          Transposition y --> z
             CALL cpu_log( log_point_s(5), 'transpo forward', 'continue', cpu_log_nowait )

             DO  jnew = 0, ny-1
                jind = jj + ( nblk + 1 ) * jnew
                f_in(:,:,jnew) = f_inv(:,:,jind)
             ENDDO

             CALL transpose_yz( f_in, ar1(:,:,:) )

             IF ( jj == nblk )  THEN
                CALL cpu_log( log_point_s(5), 'transpo forward', 'stop' )
             ELSE
                CALL cpu_log( log_point_s(5), 'transpo forward', 'pause' )
             ENDIF

!
!--          Solve the tridiagonal equation system along z
             CALL cpu_log( log_point_s(6), 'tridia', 'start', cpu_log_nowait )

             n = isave(2) + jj
             CALL tridia_substi_overlap( ar1(:,:,:), n )

             CALL cpu_log( log_point_s(6), 'tridia', 'stop' )

!
!--          Inverse Fourier Transformation
!--          Transposition z --> y
!--          Only one thread should call MPI routines, therefore forward and
!--          backward tranpose are in the same section
             IF ( jj == 0 )  THEN
                CALL cpu_log( log_point_s(8), 'transpo invers', 'start', cpu_log_nowait )
             ELSE
                CALL cpu_log( log_point_s(8), 'transpo invers', 'continue', cpu_log_nowait )
             ENDIF

             CALL transpose_zy( ar1(:,:,:), f_in )

             DO  jnew = 0, ny-1
                jind = jj + ( nblk + 1 ) * jnew
                f_inv(:,:,jind) = f_in(:,:,jnew)
             ENDDO

             CALL cpu_log( log_point_s(8), 'transpo invers', 'pause' )

          ENDDO
!
!--       Restore original indices/counters
          ny               = isave(1)
          nys_z            = isave(2)
          nyn_z            = isave(3)
          sendrecvcount_yz = isave(4)

          CALL resort_for_zy( f_inv, f_out_y )

          DEALLOCATE( ar1, f_in, f_inv )

!
!--       fft along y backward + subsequent transposition y --> x
          ALLOCATE( f_inv(nys_x:nyn_x,nzb_x:nzt_x,0:nx) )
!
!--       Save original indices and gridpoint counter
          isave(1) = nx
          isave(2) = nxl_y
          isave(3) = nxr_y
          isave(4) = sendrecvcount_xy
!
!--       Set new indices for transformation
          nblk             = (( nx+1 ) / pdims(2) ) - 1
          nx               = pdims(2)
          nnx_y            = 1
          nxl_y            = myidy * nnx_y
          nxr_y            = ( myidy + 1 ) * nnx_y - 1
          sendrecvcount_xy = nnx_y * ( nyn_x-nys_x+1 ) * ( nzt_x-nzb_x+1 )

          ALLOCATE( ar1(0:ny,nxl_y:nxr_y,nzb_y:nzt_y) )
          ALLOCATE( f_in(nys_x:nyn_x,nzb_x:nzt_x,0:nx) )

          DO  ii = 0, nblk

             CALL cpu_log( log_point_s(7), 'fft_y', 'continue', cpu_log_nowait )

             n = isave(2) + ii
             nxl_y_bound = isave(2)
             nxr_y_bound = isave(3)

             CALL fft_y( ar1(:,:,:), 'backward', ar_tr = f_out_y,              &
                         nxl_y_bound = nxl_y_bound, nxr_y_bound = nxr_y_bound, &
                         nxl_y_l = n, nxr_y_l = n )

             IF ( ii == nblk )  THEN
                CALL cpu_log( log_point_s(7), 'fft_y', 'stop' )
             ELSE
                CALL cpu_log( log_point_s(7), 'fft_y', 'pause' )
             ENDIF

             CALL cpu_log( log_point_s(8), 'transpo invers', 'continue', cpu_log_nowait )

             CALL transpose_yx( ar1(:,:,:), f_in )

             DO  inew = 0, nx-1
                iind = ii + (nblk+1) * inew
                f_inv(:,:,iind) = f_in(:,:,inew)
             ENDDO

             CALL cpu_log( log_point_s(8), 'transpo invers', 'pause' )

          ENDDO
!
!--       Restore original indices/counters
          nx               = isave(1)
          nxl_y            = isave(2)
          nxr_y            = isave(3)
          sendrecvcount_xy = isave(4)

          CALL resort_for_yx( f_inv, f_out_z )

          DEALLOCATE( ar1, f_in, f_inv )

!
!--       fft along x backward + subsequent final transposition x --> z
          ALLOCATE( f_inv(nys:nyn,nxl:nxr,1:nz) )
!
!--       Save original indices and gridpoint counter
          isave(1) = nz
          isave(2) = nzb_x
          isave(3) = nzt_x
          isave(4) = sendrecvcount_zx
!
!--       Set new indices for transformation
          nblk             = nz / pdims(1)
          nz               = pdims(1)
          nnz_x            = 1
          nzb_x            = 1 + myidx * nnz_x
          nzt_x            = ( myidx + 1 ) * nnz_x
          sendrecvcount_zx = nnx * nny * nnz_x

          ALLOCATE( ar1(0:nx,nys_x:nyn_x,nzb_x:nzt_x) )
          ALLOCATE( f_in(nys:nyn,nxl:nxr,1:nz) )

          DO  kk = 1, nblk

             CALL cpu_log( log_point_s(4), 'fft_x', 'continue', cpu_log_nowait )

             n = isave(2) + kk - 1
             CALL fft_x( ar1(:,:,:), 'backward', f_out_z(:,:,n))

             IF ( kk == nblk )  THEN
                CALL cpu_log( log_point_s(4), 'fft_x', 'stop' )
             ELSE
                CALL cpu_log( log_point_s(4), 'fft_x', 'pause' )
             ENDIF

             CALL cpu_log( log_point_s(8), 'transpo invers', 'continue', cpu_log_nowait )

             CALL transpose_xz( ar1(:,:,:), f_in )

             DO  knew = 1, nz
                ki = kk + nblk * (knew-1)
                f_inv(:,:,ki) = f_in(:,:,knew)
             ENDDO

             IF ( kk == nblk )  THEN
                CALL cpu_log( log_point_s(8), 'transpo invers', 'stop' )
             ELSE
                CALL cpu_log( log_point_s(8), 'transpo invers', 'pause' )
             ENDIF

          ENDDO
!
!--       Restore original indices/counters
          nz               = isave(1)
          nzb_x            = isave(2)
          nzt_x            = isave(3)
          sendrecvcount_zx = isave(4)

          CALL resort_for_xz( f_inv, ar )

          DEALLOCATE( ar1, f_in, f_inv )

       ENDIF

       CALL cpu_log( log_point_s(3), 'poisfft', 'stop' )

    END SUBROUTINE poisfft


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Fourier-transformation along y with subsequent transposition y --> x for
!> a 1d-decomposition along x.
!>
!> @attention The performance of this routine is much faster on the NEC-SX6,
!>            if the first index of work_ffty_vec is odd. Otherwise
!>            memory bank conflicts may occur (especially if the index is a
!>            multiple of 128). That's why work_ffty_vec is dimensioned as
!>            0:ny+1.
!>            Of course, this will not work if users are using an odd number
!>            of gridpoints along y.
!------------------------------------------------------------------------------!
    SUBROUTINE ffty_tr_yx( f_in, f_out )

       USE control_parameters,                                                 &
           ONLY:  loop_optimization

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE kinds

       USE pegrid

       IMPLICIT NONE

       INTEGER(iwp)            ::  i            !<
       INTEGER(iwp)            ::  iend         !<
       INTEGER(iwp)            ::  iouter       !<
       INTEGER(iwp)            ::  ir           !<
       INTEGER(iwp)            ::  j            !<
       INTEGER(iwp)            ::  k            !<

       INTEGER(iwp), PARAMETER ::  stridex = 4  !<

       REAL(wp), DIMENSION(1:nz,0:ny,nxl:nxr)             ::  f_in   !<
       REAL(wp), DIMENSION(nnx,1:nz,nys_x:nyn_x,pdims(1)) ::  f_out  !<
       REAL(wp), DIMENSION(nxl:nxr,1:nz,0:ny)             ::  work   !<

       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  work_ffty      !<
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  work_ffty_vec  !<

!
!--    Carry out the FFT along y, where all data are present due to the
!--    1d-decomposition along x. Resort the data in a way that x becomes
!--    the first index.
       CALL cpu_log( log_point_s(7), 'fft_y_1d', 'start' )

       IF ( loop_optimization == 'vector' )  THEN

          ALLOCATE( work_ffty_vec(0:ny+1,1:nz,nxl:nxr) )
!
!--       Code optimized for vector processors
          !$OMP PARALLEL PRIVATE ( i, j, k )
          !$OMP DO
          DO  i = nxl, nxr

             DO  j = 0, ny
                DO  k = 1, nz
                   work_ffty_vec(j,k,i) = f_in(k,j,i)
                ENDDO
             ENDDO

             CALL fft_y_m( work_ffty_vec(:,:,i), ny+1, 'forward' )

          ENDDO

          !$OMP DO
          DO  k = 1, nz
             DO  j = 0, ny
                DO  i = nxl, nxr
                   work(i,k,j) = work_ffty_vec(j,k,i)
                ENDDO
             ENDDO
          ENDDO
          !$OMP END PARALLEL

          DEALLOCATE( work_ffty_vec )

       ELSE
!
!--       Cache optimized code.
          ALLOCATE( work_ffty(0:ny,stridex) )
!
!--       The i-(x-)direction is split into a strided outer loop and an inner
!--       loop for better cache performance
          !$OMP PARALLEL PRIVATE (i,iend,iouter,ir,j,k,work_ffty)
          !$OMP DO
          DO  iouter = nxl, nxr, stridex

             iend = MIN( iouter+stridex-1, nxr )  ! Upper bound for inner i loop

             DO  k = 1, nz

                DO  i = iouter, iend

                   ir = i-iouter+1  ! counter within a stride
                   DO  j = 0, ny
                      work_ffty(j,ir) = f_in(k,j,i)
                   ENDDO
!
!--                FFT along y
                   CALL fft_y_1d( work_ffty(:,ir), 'forward' )

                ENDDO

!
!--             Resort
                DO  j = 0, ny
                   DO  i = iouter, iend
                      work(i,k,j) = work_ffty(j,i-iouter+1)
                   ENDDO
                ENDDO

             ENDDO

          ENDDO
          !$OMP END PARALLEL

          DEALLOCATE( work_ffty )

       ENDIF

       CALL cpu_log( log_point_s(7), 'fft_y_1d', 'pause' )

!
!--    Transpose array
#if defined( __parallel )
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start' )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLTOALL( work(nxl,1,0),      sendrecvcount_xy, MPI_REAL, &
                          f_out(1,1,nys_x,1), sendrecvcount_xy, MPI_REAL, &
                          comm1dx, ierr )
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )
#endif

    END SUBROUTINE ffty_tr_yx


!------------------------------------------------------------------------------!
! Description:
! ------------
!>  Transposition x --> y with a subsequent backward Fourier transformation for
!>  a 1d-decomposition along x
!------------------------------------------------------------------------------!
    SUBROUTINE tr_xy_ffty( f_in, f_out )

       USE control_parameters,                                                 &
           ONLY:  loop_optimization

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE kinds

       USE pegrid

       IMPLICIT NONE

       INTEGER(iwp)            ::  i            !<
       INTEGER(iwp)            ::  iend         !<
       INTEGER(iwp)            ::  iouter       !<
       INTEGER(iwp)            ::  ir           !<
       INTEGER(iwp)            ::  j            !<
       INTEGER(iwp)            ::  k            !<

       INTEGER(iwp), PARAMETER ::  stridex = 4  !<

       REAL(wp), DIMENSION(nnx,1:nz,nys_x:nyn_x,pdims(1)) ::  f_in   !<
       REAL(wp), DIMENSION(1:nz,0:ny,nxl:nxr)             ::  f_out  !<
       REAL(wp), DIMENSION(nxl:nxr,1:nz,0:ny)             ::  work   !<

       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  work_ffty         !<
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  work_ffty_vec     !<

!
!--    Transpose array
#if defined( __parallel )
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start' )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLTOALL( f_in(1,1,nys_x,1), sendrecvcount_xy, MPI_REAL, &
                          work(nxl,1,0),     sendrecvcount_xy, MPI_REAL, &
                          comm1dx, ierr )
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )
#endif

!
!--    Resort the data in a way that y becomes the first index and carry out the
!--    backward fft along y.
       CALL cpu_log( log_point_s(7), 'fft_y_1d', 'continue' )

       IF ( loop_optimization == 'vector' )  THEN

          ALLOCATE( work_ffty_vec(0:ny+1,1:nz,nxl:nxr) )
!
!--       Code optimized for vector processors
          !$OMP PARALLEL PRIVATE ( i, j, k )
          !$OMP DO
          DO  k = 1, nz
             DO  j = 0, ny
                DO  i = nxl, nxr
                   work_ffty_vec(j,k,i) = work(i,k,j)
                ENDDO
             ENDDO
          ENDDO

          !$OMP DO
          DO  i = nxl, nxr

             CALL fft_y_m( work_ffty_vec(:,:,i), ny+1, 'backward' )

             DO  j = 0, ny
                DO  k = 1, nz
                   f_out(k,j,i) = work_ffty_vec(j,k,i)
                ENDDO
             ENDDO

          ENDDO
          !$OMP END PARALLEL

          DEALLOCATE( work_ffty_vec )

       ELSE
!
!--       Cache optimized code.
          ALLOCATE( work_ffty(0:ny,stridex) )
!
!--       The i-(x-)direction is split into a strided outer loop and an inner
!--       loop for better cache performance
          !$OMP PARALLEL PRIVATE ( i, iend, iouter, ir, j, k, work_ffty )
          !$OMP DO
          DO  iouter = nxl, nxr, stridex

             iend = MIN( iouter+stridex-1, nxr )  ! Upper bound for inner i loop

             DO  k = 1, nz
!
!--             Resort
                DO  j = 0, ny
                   DO  i = iouter, iend
                      work_ffty(j,i-iouter+1) = work(i,k,j)
                   ENDDO
                ENDDO

                DO  i = iouter, iend

!
!--                FFT along y
                   ir = i-iouter+1  ! counter within a stride
                   CALL fft_y_1d( work_ffty(:,ir), 'backward' )

                   DO  j = 0, ny
                      f_out(k,j,i) = work_ffty(j,ir)
                   ENDDO
                ENDDO

             ENDDO

          ENDDO
          !$OMP END PARALLEL

          DEALLOCATE( work_ffty )

       ENDIF

       CALL cpu_log( log_point_s(7), 'fft_y_1d', 'stop' )

    END SUBROUTINE tr_xy_ffty


!------------------------------------------------------------------------------!
! Description:
! ------------
!> FFT along x, solution of the tridiagonal system and backward FFT for
!> a 1d-decomposition along x
!>
!> @warning this subroutine may still not work for hybrid parallelization
!>          with OpenMP (for possible necessary changes see the original
!>          routine poisfft_hybrid, developed by Klaus Ketelsen, May 2002)
!------------------------------------------------------------------------------!
    SUBROUTINE fftx_tri_fftx( ar )

       USE control_parameters,                                                 &
           ONLY:  loop_optimization

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE grid_variables,                                                     &
           ONLY:  ddx2, ddy2

       USE kinds

       USE pegrid

       IMPLICIT NONE

       INTEGER(iwp) ::  i                   !<
       INTEGER(iwp) ::  j                   !<
       INTEGER(iwp) ::  k                   !<
       INTEGER(iwp) ::  m                   !<
       INTEGER(iwp) ::  n                   !<
       INTEGER(iwp) ::  omp_get_thread_num  !<
       INTEGER(iwp) ::  tn                  !<

       REAL(wp), DIMENSION(0:nx)                          ::  work_fftx  !<
       REAL(wp), DIMENSION(0:nx,1:nz)                     ::  work_trix  !<
       REAL(wp), DIMENSION(nnx,1:nz,nys_x:nyn_x,pdims(1)) ::  ar         !<
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE          ::  tri        !<


       CALL cpu_log( log_point_s(33), 'fft_x_1d + tridia', 'start' )

       ALLOCATE( tri(5,0:nx,0:nz-1,0:threads_per_task-1) )

       tn = 0              ! Default thread number in case of one thread
!$OMP  PARALLEL DO PRIVATE ( i, j, k, m, n, tn, work_fftx, work_trix )
       DO  j = nys_x, nyn_x

!$        tn = omp_get_thread_num()

          IF ( loop_optimization == 'vector' )  THEN
!
!--          Code optimized for vector processors
             DO  k = 1, nz

                m = 0
                DO  n = 1, pdims(1)
                   DO  i = 1, nnx
                      work_trix(m,k) = ar(i,k,j,n)
                      m = m + 1
                   ENDDO
                ENDDO

             ENDDO

             CALL fft_x_m( work_trix, 'forward' )

          ELSE
!
!--          Cache optimized code
             DO  k = 1, nz

                m = 0
                DO  n = 1, pdims(1)
                   DO  i = 1, nnx
                      work_fftx(m) = ar(i,k,j,n)
                      m = m + 1
                   ENDDO
                ENDDO

                CALL fft_x_1d( work_fftx, 'forward' )

                DO  i = 0, nx
                   work_trix(i,k) = work_fftx(i)
                ENDDO

             ENDDO

          ENDIF

!
!--       Solve the linear equation system
          CALL tridia_1dd( ddx2, ddy2, nx, ny, j, work_trix, tri(:,:,:,tn) )

          IF ( loop_optimization == 'vector' )  THEN
!
!--          Code optimized for vector processors
             CALL fft_x_m( work_trix, 'backward' )

             DO  k = 1, nz

                m = 0
                DO  n = 1, pdims(1)
                   DO  i = 1, nnx
                      ar(i,k,j,n) = work_trix(m,k)
                      m = m + 1
                   ENDDO
                ENDDO

             ENDDO

          ELSE
!
!--          Cache optimized code
             DO  k = 1, nz

                DO  i = 0, nx
                   work_fftx(i) = work_trix(i,k)
                ENDDO

                CALL fft_x_1d( work_fftx, 'backward' )

                m = 0
                DO  n = 1, pdims(1)
                   DO  i = 1, nnx
                      ar(i,k,j,n) = work_fftx(m)
                      m = m + 1
                   ENDDO
                ENDDO

             ENDDO

          ENDIF

       ENDDO

       DEALLOCATE( tri )

       CALL cpu_log( log_point_s(33), 'fft_x_1d + tridia', 'stop' )

    END SUBROUTINE fftx_tri_fftx


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Fourier-transformation along x with subsequent transposition x --> y for
!> a 1d-decomposition along y.
!>
!> @attention NEC-branch of this routine may significantly profit from
!>            further optimizations. So far, performance is much worse than
!>            for routine ffty_tr_yx (more than three times slower).
!------------------------------------------------------------------------------!
    SUBROUTINE fftx_tr_xy( f_in, f_out )


       USE control_parameters,                                                 &
           ONLY:  loop_optimization

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE kinds

       USE pegrid

       IMPLICIT NONE

       INTEGER(iwp) ::  i  !<
       INTEGER(iwp) ::  j  !<
       INTEGER(iwp) ::  k  !<

       REAL(wp), DIMENSION(0:nx,1:nz,nys:nyn)             ::  work_fftx  !<
       REAL(wp), DIMENSION(1:nz,nys:nyn,0:nx)             ::  f_in       !<
       REAL(wp), DIMENSION(nny,1:nz,nxl_y:nxr_y,pdims(2)) ::  f_out      !<
       REAL(wp), DIMENSION(nys:nyn,1:nz,0:nx)             ::  work       !<

!
!--    Carry out the FFT along x, where all data are present due to the
!--    1d-decomposition along y. Resort the data in a way that y becomes
!--    the first index.
       CALL cpu_log( log_point_s(4), 'fft_x_1d', 'start' )

       IF ( loop_optimization == 'vector' )  THEN
!
!--       Code for vector processors
!$OMP     PARALLEL PRIVATE ( i, j, k )
!$OMP     DO
          DO  i = 0, nx

             DO  j = nys, nyn
                DO  k = 1, nz
                   work_fftx(i,k,j) = f_in(k,j,i)
                ENDDO
             ENDDO

          ENDDO

!$OMP     DO
          DO  j = nys, nyn

             CALL fft_x_m( work_fftx(:,:,j), 'forward' )

             DO  k = 1, nz
                DO  i = 0, nx
                   work(j,k,i) = work_fftx(i,k,j)
                ENDDO
             ENDDO

          ENDDO
!$OMP     END PARALLEL

       ELSE

!
!--       Cache optimized code (there might be still a potential for better
!--       optimization).
!$OMP     PARALLEL PRIVATE (i,j,k)
!$OMP     DO
          DO  i = 0, nx

             DO  j = nys, nyn
                DO  k = 1, nz
                   work_fftx(i,k,j) = f_in(k,j,i)
                ENDDO
             ENDDO

          ENDDO

!$OMP     DO
          DO  j = nys, nyn
             DO  k = 1, nz

                CALL fft_x_1d( work_fftx(0:nx,k,j), 'forward' )

                DO  i = 0, nx
                   work(j,k,i) = work_fftx(i,k,j)
                ENDDO
             ENDDO

          ENDDO
!$OMP     END PARALLEL

       ENDIF
       CALL cpu_log( log_point_s(4), 'fft_x_1d', 'pause' )

!
!--    Transpose array
#if defined( __parallel )
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start' )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLTOALL( work(nys,1,0),      sendrecvcount_xy, MPI_REAL, &
                          f_out(1,1,nxl_y,1), sendrecvcount_xy, MPI_REAL, &
                          comm1dy, ierr )
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )
#endif

    END SUBROUTINE fftx_tr_xy


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Transposition y --> x with a subsequent backward Fourier transformation for
!> a 1d-decomposition along x.
!------------------------------------------------------------------------------!
    SUBROUTINE tr_yx_fftx( f_in, f_out )


       USE control_parameters,                                                 &
           ONLY:  loop_optimization

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE kinds

       USE pegrid

       IMPLICIT NONE

       INTEGER(iwp) ::  i  !<
       INTEGER(iwp) ::  j  !<
       INTEGER(iwp) ::  k  !<

       REAL(wp), DIMENSION(0:nx,1:nz,nys:nyn)             ::  work_fftx  !<
       REAL(wp), DIMENSION(nny,1:nz,nxl_y:nxr_y,pdims(2)) ::  f_in       !<
       REAL(wp), DIMENSION(1:nz,nys:nyn,0:nx)             ::  f_out      !<
       REAL(wp), DIMENSION(nys:nyn,1:nz,0:nx)             ::  work       !<

!
!--    Transpose array
#if defined( __parallel )
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start' )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLTOALL( f_in(1,1,nxl_y,1), sendrecvcount_xy, MPI_REAL, &
                          work(nys,1,0),     sendrecvcount_xy, MPI_REAL, &
                          comm1dy, ierr )
       CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )
#endif

!
!--    Carry out the FFT along x, where all data are present due to the
!--    1d-decomposition along y. Resort the data in a way that y becomes
!--    the first index.
       CALL cpu_log( log_point_s(4), 'fft_x_1d', 'continue' )

       IF ( loop_optimization == 'vector' )  THEN
!
!--       Code optimized for vector processors
!$OMP     PARALLEL PRIVATE ( i, j, k )
!$OMP     DO
          DO  j = nys, nyn

             DO  k = 1, nz
                DO  i = 0, nx
                   work_fftx(i,k,j) = work(j,k,i)
                ENDDO
             ENDDO

             CALL fft_x_m( work_fftx(:,:,j), 'backward' )

          ENDDO

!$OMP     DO
          DO  i = 0, nx
             DO  j = nys, nyn
                DO  k = 1, nz
                   f_out(k,j,i) = work_fftx(i,k,j)
                ENDDO
             ENDDO
          ENDDO
!$OMP     END PARALLEL

       ELSE

!
!--       Cache optimized code (there might be still a potential for better
!--       optimization).
!$OMP     PARALLEL PRIVATE (i,j,k)
!$OMP     DO
          DO  j = nys, nyn
             DO  k = 1, nz

                DO  i = 0, nx
                   work_fftx(i,k,j) = work(j,k,i)
                ENDDO

                CALL fft_x_1d( work_fftx(0:nx,k,j), 'backward' )

             ENDDO
          ENDDO

!$OMP     DO
          DO  i = 0, nx
             DO  j = nys, nyn
                DO  k = 1, nz
                   f_out(k,j,i) = work_fftx(i,k,j)
                ENDDO
             ENDDO
          ENDDO
!$OMP     END PARALLEL

       ENDIF
       CALL cpu_log( log_point_s(4), 'fft_x_1d', 'stop' )

    END SUBROUTINE tr_yx_fftx


!------------------------------------------------------------------------------!
! Description:
! ------------
!> FFT along y, solution of the tridiagonal system and backward FFT for
!> a 1d-decomposition along y.
!>
!> @warning this subroutine may still not work for hybrid parallelization
!>          with OpenMP (for possible necessary changes see the original
!>          routine poisfft_hybrid, developed by Klaus Ketelsen, May 2002)
!------------------------------------------------------------------------------!
    SUBROUTINE ffty_tri_ffty( ar )


       USE control_parameters,                                                 &
           ONLY:  loop_optimization

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE grid_variables,                                                     &
           ONLY:  ddx2, ddy2

       USE kinds

       USE pegrid

       IMPLICIT NONE

       INTEGER(iwp) ::  i                   !<
       INTEGER(iwp) ::  j                   !<
       INTEGER(iwp) ::  k                   !<
       INTEGER(iwp) ::  m                   !<
       INTEGER(iwp) ::  n                   !<
       INTEGER(iwp) ::  omp_get_thread_num  !<
       INTEGER(iwp) ::  tn                  !<

       REAL(wp), DIMENSION(0:ny)                          ::  work_ffty  !<
       REAL(wp), DIMENSION(0:ny,1:nz)                     ::  work_triy  !<
       REAL(wp), DIMENSION(nny,1:nz,nxl_y:nxr_y,pdims(2)) ::  ar         !<
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE          ::  tri        !<


       CALL cpu_log( log_point_s(39), 'fft_y_1d + tridia', 'start' )

       ALLOCATE( tri(5,0:ny,0:nz-1,0:threads_per_task-1) )

       tn = 0           ! Default thread number in case of one thread
!$OMP  PARALLEL DO PRIVATE ( i, j, k, m, n, tn, work_ffty, work_triy )
       DO  i = nxl_y, nxr_y

!$        tn = omp_get_thread_num()

          IF ( loop_optimization == 'vector' )  THEN
!
!--          Code optimized for vector processors
             DO  k = 1, nz

                m = 0
                DO  n = 1, pdims(2)
                   DO  j = 1, nny
                      work_triy(m,k) = ar(j,k,i,n)
                      m = m + 1
                   ENDDO
                ENDDO

             ENDDO

             CALL fft_y_m( work_triy, ny, 'forward' )

          ELSE
!
!--          Cache optimized code
             DO  k = 1, nz

                m = 0
                DO  n = 1, pdims(2)
                   DO  j = 1, nny
                      work_ffty(m) = ar(j,k,i,n)
                      m = m + 1
                   ENDDO
                ENDDO

                CALL fft_y_1d( work_ffty, 'forward' )

                DO  j = 0, ny
                   work_triy(j,k) = work_ffty(j)
                ENDDO

             ENDDO

          ENDIF

!
!--       Solve the linear equation system
          CALL tridia_1dd( ddy2, ddx2, ny, nx, i, work_triy, tri(:,:,:,tn) )

          IF ( loop_optimization == 'vector' )  THEN
!
!--          Code optimized for vector processors
             CALL fft_y_m( work_triy, ny, 'backward' )

             DO  k = 1, nz

                m = 0
                DO  n = 1, pdims(2)
                   DO  j = 1, nny
                      ar(j,k,i,n) = work_triy(m,k)
                      m = m + 1
                   ENDDO
                ENDDO

             ENDDO

          ELSE
!
!--          Cache optimized code
             DO  k = 1, nz

                DO  j = 0, ny
                   work_ffty(j) = work_triy(j,k)
                ENDDO

                CALL fft_y_1d( work_ffty, 'backward' )

                m = 0
                DO  n = 1, pdims(2)
                   DO  j = 1, nny
                      ar(j,k,i,n) = work_ffty(m)
                      m = m + 1
                   ENDDO
                ENDDO

             ENDDO

          ENDIF

       ENDDO

       DEALLOCATE( tri )

       CALL cpu_log( log_point_s(39), 'fft_y_1d + tridia', 'stop' )

    END SUBROUTINE ffty_tri_ffty

 END MODULE poisfft_mod
