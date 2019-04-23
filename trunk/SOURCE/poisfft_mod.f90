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

#ifdef __GPU
    USE cudafor
#endif

#define MY_DEBUG print *,"DEBUG",__LINE__,__FILE__

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

       INTEGER(iwp) ::  i,j,k

       REAL(wp), DIMENSION(1:nz,nys:nyn,nxl:nxr) ::  ar      !<
#ifdef __GPU
       REAL(wp), DEVICE, ALLOCATABLE ::                                        &
                ar_inv_d(:,:,:), ar_x_d(:,:,:),                                &
                ar_inv_x_d(:,:,:), ar_d(:,:,:),                                &
                ar_y_d(:,:,:), ar_inv_y_d(:,:,:),                              &
                ar_z_d(:,:,:)
#else
       REAL(wp), DIMENSION(nys:nyn,nxl:nxr,1:nz) ::  ar_inv  !<
#endif

       CALL cpu_log( log_point_s(3), 'poisfft', 'start' )

       IF ( .NOT. poisfft_initialized )  CALL poisfft_init

#ifdef __GPU
       ALLOCATE(ar_d(1:nz,nys:nyn,nxl:nxr))
       ALLOCATE(ar_y_d(0:ny,nxl_y:nxr_y,nzb_y:nzt_y))
       ALLOCATE(ar_z_d(nxl_z:nxr_z,nys_z:nyn_z,1:nz))
       ALLOCATE(ar_x_d(0:nx,nys_x:nyn_x,nzb_x:nzt_x))
#endif
!

#ifdef __GPU

       CALL cpu_log(log_point_s(5), 'fft routines', 'start')

       ar_d = ar

!--    2d-domain-decomposition or no decomposition (1 PE run)
!--    Transposition z --> x

       !$acc parallel
       !$acc loop collapse(3)
       DO  k = 1,nz
           DO  i = nxl, nxr
               DO  j = nys, nyn
                   ar_x_d(i,j,k) = ar_d(k,j,i)
               ENDDO
           ENDDO
       ENDDO
       !$acc end parallel

       CALL fft_x( ar_x_d, 'forward' )

!--    Transposition x --> y

       !$acc parallel
       !$acc loop collapse(3)
       DO  i = 0, nx
           DO  k = nzb_x, nzt_x
               DO  j = nys_x, nyn_x
                   ar_y_d(j,i,k) = ar_x_d(i,j,k)
               ENDDO
           ENDDO
       ENDDO
       !$acc end parallel

       CALL fft_y( ar_y_d, 'forward' )

       !$acc parallel
       !$acc loop collapse(3)
       DO  j = 0, ny
           DO  k = nzb_y, nzt_y
               DO  i = nxl_y, nxr_y
                   ar_z_d(i,j,k) = ar_y_d(j,i,k)
               ENDDO
           ENDDO
       ENDDO
       !$acc end parallel

       !ar_z = ar_z_d
!
!--    Solve the tridiagonal equation system along z
       CALL tridia_substi( ar_z_d )

!
!--    Inverse Fourier Transformation
!--    Transposition z --> y

       !$acc parallel
       !$acc loop collapse(3)
       DO  k = nzb_y, nzt_y
          DO  j = 0, ny
             DO  i = nxl_y, nxr_y
                ar_y_d(j,i,k) = ar_z_d(i,j,k)
             ENDDO
          ENDDO
       ENDDO
       !$acc end parallel

       CALL fft_y( ar_y_d, 'backward')
!
!--    Transposition y --> x

       !$acc parallel
       !$acc loop collapse(3)
       DO  i = nxl_y, nxr_y
          DO  k = nzb_y, nzt_y
             DO  j = 0, ny
                ar_x_d(i,j,k) = ar_y_d(j,i,k)
             ENDDO
          ENDDO
       ENDDO
       !$acc end parallel

       CALL fft_x( ar_x_d, 'backward' )

!--    Transposition x --> z

       !$acc parallel
       !$acc loop collapse(3)
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = 1, nz
                ar_d(k,j,i) = ar_x_d(i,j,k)
             ENDDO
          ENDDO
       ENDDO
       !$acc end parallel

       ar = ar_d

       CALL cpu_log(log_point_s(5), 'fft routines', 'stop')

#else
!
!--    2d-domain-decomposition or no decomposition (1 PE run)
!--    Transposition z --> x
       CALL cpu_log( log_point_s(5), 'transpo forward', 'start' )
       CALL resort_for_zx( ar, ar_inv )
       CALL transpose_zx( ar_inv, ar )
       CALL cpu_log( log_point_s(5), 'transpo forward', 'pause' )

       CALL cpu_log( log_point_s(4), 'fft_x', 'start' )
       CALL fft_x( ar, 'forward' )
       CALL cpu_log( log_point_s(4), 'fft_x', 'pause' )

!
!--    Transposition x --> y
       CALL cpu_log( log_point_s(5), 'transpo forward', 'continue' )
       CALL resort_for_xy( ar, ar_inv )
       CALL transpose_xy( ar_inv, ar )
       CALL cpu_log( log_point_s(5), 'transpo forward', 'pause' )

       CALL cpu_log( log_point_s(7), 'fft_y', 'start' )
       CALL fft_y( ar, 'forward')
       CALL cpu_log( log_point_s(7), 'fft_y', 'pause' )

!
!--    Transposition y --> z
       CALL cpu_log( log_point_s(5), 'transpo forward', 'continue' )
       CALL resort_for_yz( ar, ar_inv )
       CALL transpose_yz( ar_inv, ar )
       CALL cpu_log( log_point_s(5), 'transpo forward', 'stop' )

!
!--    Solve the tridiagonal equation system along z
       CALL cpu_log( log_point_s(6), 'tridia', 'start' )
       CALL tridia_substi( ar )
       CALL cpu_log( log_point_s(6), 'tridia', 'stop' )

!
!--    Inverse Fourier Transformation
!--    Transposition z --> y
       CALL cpu_log( log_point_s(8), 'transpo invers', 'start' )
       CALL transpose_zy( ar, ar_inv )
       CALL resort_for_zy( ar_inv, ar )
       CALL cpu_log( log_point_s(8), 'transpo invers', 'pause' )

       CALL cpu_log( log_point_s(7), 'fft_y', 'continue' )
       CALL fft_y( ar, 'backward')
       CALL cpu_log( log_point_s(7), 'fft_y', 'stop' )

!
!--    Transposition y --> x
       CALL cpu_log( log_point_s(8), 'transpo invers', 'continue' )
       CALL transpose_yx( ar, ar_inv )
       CALL resort_for_yx( ar_inv, ar )
       CALL cpu_log( log_point_s(8), 'transpo invers', 'pause' )

       CALL cpu_log( log_point_s(4), 'fft_x', 'continue' )
       CALL fft_x( ar, 'backward' )
       CALL cpu_log( log_point_s(4), 'fft_x', 'stop' )

!
!--    Transposition x --> z
       CALL cpu_log( log_point_s(8), 'transpo invers', 'continue' )
       CALL transpose_xz( ar, ar_inv )
       CALL resort_for_xz( ar_inv, ar )
       CALL cpu_log( log_point_s(8), 'transpo invers', 'stop' )

#endif

       CALL cpu_log( log_point_s(3), 'poisfft', 'stop' )

#ifdef __GPU
       DEALLOCATE(ar_inv_d, ar_x_d, ar_inv_x_d,ar_d,ar_y_d, ar_inv_y_d, ar_z_d)
#endif

    END SUBROUTINE poisfft

 END MODULE poisfft_mod
