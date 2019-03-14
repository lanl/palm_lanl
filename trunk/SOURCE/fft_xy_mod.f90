!> @file fft_xy_mod.f90
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
! $Id: fft_xy_mod.f90 3045 2018-05-28 07:55:41Z Giersch $
! Error messages revised
!
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2300 2017-06-29 13:31:14Z raasch
! NEC related code partly removed, host replaced by loop_optimization
!
! 2274 2017-06-09 13:27:48Z Giersch
! Changed error messages
!
! 2119 2017-01-17 16:51:50Z raasch
!
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC directives and CUDA-fft related code removed
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
!
! 1815 2016-04-06 13:49:59Z raasch
! cpp-directives for ibmy removed
!
! 1749 2016-02-09 12:19:56Z raasch
! small OpenACC bugfix
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1600 2015-06-11 15:50:12Z raasch
! bugfix: openMP threadprivate statement moved after variable declaration
!
! 1482 2014-10-18 12:34:45Z raasch
! cudafft workaround for data declaration of ar_tmp because of PGI 14.1 bug
!
! 1402 2014-05-09 14:25:13Z raasch
! fortran bugfix for r1392
!
! 1398 2014-05-07 11:15:00Z heinze
! bugfix: typo removed for KIND in CMPLX function
!
! 1392 2014-05-06 09:10:05Z raasch
! bugfix: KIND attribute added to CMPLX functions
!
! 1374 2014-04-25 12:55:07Z raasch
! bugfixes: missing variables added to ONLY list, dpk renamed dp
!
! 1372 2014-04-24 06:29:32Z raasch
! openMP-bugfix for fftw: some arrays defined as threadprivate
!
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1342 2014-03-26 17:04:47Z kanani
! REAL constants defined as wp-kind
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL functions provided with KIND-attribute
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
! 1304 2014-03-12 10:29:42Z raasch
! openmp bugfix: work1 used in Temperton algorithm must be private
!
! 1257 2013-11-08 15:18:40Z raasch
! openacc loop and loop vector clauses removed, declare create moved after
! the FORTRAN declaration statement
!
! 1219 2013-08-30 09:33:18Z heinze
! bugfix: use own branch for fftw
!
! 1216 2013-08-26 09:31:42Z raasch
! fft_x and fft_y modified for parallel / ovverlapping execution of fft and
! transpositions,
! fftw implemented for 1d-decomposition (fft_x_1d, fft_y_1d)
!
! 1210 2013-08-14 10:58:20Z raasch
! fftw added
!
! 1166 2013-05-24 13:55:44Z raasch
! C_DOUBLE/COMPLEX reset to dpk
!
! 1153 2013-05-10 14:33:08Z raasch
! code adjustment of data types for CUDA fft required by PGI 12.3 / CUDA 5.0
!
! 1111 2013-03-08 23:54:10Z raasch
! further openACC statements added, CUDA branch completely runs on GPU
! bugfix: CUDA fft plans adjusted for domain decomposition (before they always
! used total domain)
!
! 1106 2013-03-04 05:31:38Z raasch
! CUDA fft added
! array_kind renamed precision_kind, 3D- instead of 1D-loops in fft_x and fft_y
! old fft_x, fft_y become fft_x_1d, fft_y_1d and are used for 1D-decomposition
!
! 1092 2013-02-02 11:24:22Z raasch
! variable sizw declared for NEC case only
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! Revision 1.1  2002/06/11 13:00:49  raasch
! Initial revision
!
!
! Description:
! ------------
!> Fast Fourier transformation along x and y for 1d domain decomposition along x.
!> Original version: Klaus Ketelsen (May 2002)
!------------------------------------------------------------------------------!
 MODULE fft_xy


    USE control_parameters,                                                    &
        ONLY:  fft_method, message_string

    USE indices,                                                               &
        ONLY:  nx, ny, nz

#if defined( __fftw )
    USE, INTRINSIC ::  ISO_C_BINDING
#endif

    USE kinds

    USE transpose_indices,                                                     &
        ONLY:  nxl_y, nxr_y, nyn_x, nys_x, nzb_x, nzb_y, nzt_x, nzt_y

#ifdef __GPU
    USE cudafor

    USE cufft
    #endif

    IMPLICIT NONE

    PRIVATE
    PUBLIC fft_x, fft_x_1d, fft_y, fft_y_1d, fft_init, fft_x_m, fft_y_m,    &
            fft_finalize
#define MY_DEBUG print *,"DEBUG",__LINE__,__FILE__

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE, SAVE ::  ifax_x  !<
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE, SAVE ::  ifax_y  !<

    LOGICAL, SAVE ::  init_fft = .FALSE.  !<

    REAL(wp), SAVE ::  dnx      !<
    REAL(wp), SAVE ::  dny      !<
    REAL(wp), SAVE ::  sqr_dnx  !<
    REAL(wp), SAVE ::  sqr_dny  !<

    REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::  trigs_x  !<
    REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::  trigs_y  !<
#if defined( __fftw )
    INCLUDE  'fftw3.f03'
    INTEGER(KIND=C_INT) ::  nx_c  !<
    INTEGER(KIND=C_INT) ::  ny_c  !<

    COMPLEX(KIND=C_DOUBLE_COMPLEX), DIMENSION(:), ALLOCATABLE, SAVE ::  x_out  !<
    COMPLEX(KIND=C_DOUBLE_COMPLEX), DIMENSION(:), ALLOCATABLE, SAVE ::         &
       y_out  !<

    REAL(KIND=C_DOUBLE), DIMENSION(:), ALLOCATABLE, SAVE ::                    &
       x_in   !<
    REAL(KIND=C_DOUBLE), DIMENSION(:), ALLOCATABLE, SAVE ::                    &
       y_in   !<
    !$OMP THREADPRIVATE( x_out, y_out, x_in, y_in )


    TYPE(C_PTR), SAVE ::  plan_xf, plan_xi, plan_yf, plan_yi
#endif


#ifdef __GPU
    INTEGER :: nx_cC, ny_cC

    INTEGER, PARAMETER, SAVE :: batch = 1

    COMPLEX, DEVICE, DIMENSION(:), ALLOCATABLE, SAVE ::       &
            x_out_dev, y_out_dev
    REAL(wp), DEVICE, DIMENSION(:), ALLOCATABLE, SAVE :: x_in_dev, y_in_dev
    INTEGER, SAVE :: plan_xf_dev, plan_xi_dev, plan_yf_dev, plan_yi_dev
#endif
    !
!-- Public interfaces
    INTERFACE fft_init
       MODULE PROCEDURE fft_init
    END INTERFACE fft_init

    INTERFACE fft_x
       MODULE PROCEDURE fft_x
    END INTERFACE fft_x

    INTERFACE fft_x_1d
       MODULE PROCEDURE fft_x_1d
    END INTERFACE fft_x_1d

    INTERFACE fft_y
       MODULE PROCEDURE fft_y
    END INTERFACE fft_y

    INTERFACE fft_y_1d
       MODULE PROCEDURE fft_y_1d
    END INTERFACE fft_y_1d

    INTERFACE fft_x_m
       MODULE PROCEDURE fft_x_m
    END INTERFACE fft_x_m

    INTERFACE fft_y_m
       MODULE PROCEDURE fft_y_m
    END INTERFACE fft_y_m

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE fft_init

       IMPLICIT NONE
!
       INTEGER(iwp) :: ierr
!--    Return, if already called
       IF ( init_fft )  THEN
          RETURN
       ELSE
          init_fft = .TRUE.
       ENDIF
!
!--       FFTW
#if defined( __fftw )
          nx_c = nx+1
          ny_c = ny+1
          !$OMP PARALLEL
          ALLOCATE( x_in(0:nx+2), y_in(0:ny+2), x_out(0:(nx+1)/2),             &
                    y_out(0:(ny+1)/2) )
          !$OMP END PARALLEL
          plan_xf = FFTW_PLAN_DFT_R2C_1D( nx_c, x_in, x_out, FFTW_ESTIMATE )
          plan_xi = FFTW_PLAN_DFT_C2R_1D( nx_c, x_out, x_in, FFTW_ESTIMATE )
          plan_yf = FFTW_PLAN_DFT_R2C_1D( ny_c, y_in, y_out, FFTW_ESTIMATE )
          plan_yi = FFTW_PLAN_DFT_C2R_1D( ny_c, y_out, y_in, FFTW_ESTIMATE )
#else
          message_string = 'preprocessor switch for fftw is missing'
          CALL message( 'fft_init', 'PA0080', 1, 2, 0, 6, 0 )
#endif

#if defined ( __GPU )
       nx_cC = nx+1
       ny_cC = ny+1
       ALLOCATE( x_in_dev(0:nx), y_in_dev(0:ny), x_out_dev(0:(nx+1)/2),    &
               y_out_dev(0:(ny+1)/2) )

       ierr = cufftPlan1d( plan_xf_dev, nx_cC, CUFFT_R2C, batch)
       ierr = cufftPlan1d( plan_xi_dev, nx_cC, CUFFT_C2R, batch)
       ierr = cufftPlan1d( plan_yf_dev, ny_cC, CUFFT_R2C, batch)
       ierr = cufftPlan1d( plan_yi_dev, ny_cC, CUFFT_C2R, batch)
#endif
    END SUBROUTINE fft_init

    SUBROUTINE fft_finalize
            #ifdef __GPU
            DEALLOCATE(x_in_dev, y_in_dev, y_out_dev, x_out_dev)
            #endif
    END SUBROUTINE fft_finalize


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Fourier-transformation along x-direction.
!> Version for 2D-decomposition.
!> It uses internal algorithms (Singleton or Temperton) or
!> system-specific routines, if they are available
!------------------------------------------------------------------------------!

    SUBROUTINE fft_x( ar, direction )


       IMPLICIT NONE

       CHARACTER (LEN=*) ::  direction  !<

       INTEGER(iwp) :: ierr
       INTEGER(iwp) ::  i          !<
       INTEGER(iwp) ::  j          !<
       INTEGER(iwp) ::  k          !<
       LOGICAL ::  forward_fft !<
       #ifdef __GPU
       REAL(wp), DEVICE, DIMENSION(0:nx,nys_x:nyn_x,nzb_x:nzt_x) :: ar      !<
       #else
       REAL(wp), DIMENSION(0:nx,nys_x:nyn_x,nzb_x:nzt_x) :: ar      !<
       #endif

       IF ( direction == 'forward' )  THEN
          forward_fft = .TRUE.
       ELSE
          forward_fft = .FALSE.
       ENDIF

#if defined( __GPU)

       if ( forward_fft )  THEN

         DO k = nzb_x, nzt_x
            DO j = nys_x, nyn_x

               ierr = cufftExecR2C( plan_xf_dev, ar(0:nx,j,k), x_out_dev)

               !$acc kernels deviceptr(x_out_dev, ar)
               DO  i = 0, (nx+1)/2
                   ar(i,j,k) = REAL( x_out_dev(i), KIND=wp ) / ( nx+1 )
               ENDDO
               DO  i = 1, (nx+1)/2 - 1
                   ar(nx+1-i,j,k) = AIMAG( x_out_dev(i) ) / ( nx+1 )
               ENDDO
               !$acc end kernels
           ENDDO
        ENDDO

     ELSE

        DO  k = nzb_x, nzt_x
           DO  j = nys_x, nyn_x

               !$acc kernels deviceptr(ar, x_out_dev)
               x_out_dev(0) = CMPLX( ar(0,j,k), 0.0_wp, KIND=wp )
               DO  i = 1, (nx+1)/2 - 1
                   x_out_dev(i) = CMPLX( ar(i,j,k), ar(nx+1-i,j,k), KIND=wp )
               ENDDO
                   x_out_dev((nx+1)/2) = CMPLX( ar((nx+1)/2,j,k), 0.0_wp,       &
                                             KIND=wp )
               !$acc end kernels
               ierr = cufftExecC2R( plan_xi_dev, x_out_dev, ar(0:nx,j,k) )

           ENDDO
        ENDDO
    ENDIF

#else
          IF ( forward_fft )  THEN

             !$OMP PARALLEL PRIVATE ( work, i, j, k )
             !$OMP DO
             DO  k = nzb_x, nzt_x
                DO  j = nys_x, nyn_x

                   x_in(0:nx) = ar(0:nx,j,k)
                   CALL FFTW_EXECUTE_DFT_R2C( plan_xf, x_in, x_out )

                   DO  i = 0, (nx+1)/2
                      ar(i,j,k) = REAL( x_out(i), KIND=wp ) / ( nx+1 )
                   ENDDO
                   DO  i = 1, (nx+1)/2 - 1
                      ar(nx+1-i,j,k) = AIMAG( x_out(i) ) / ( nx+1 )
                   ENDDO

                ENDDO
             ENDDO
             !$OMP END PARALLEL

          ELSE
             !$OMP PARALLEL PRIVATE ( work, i, j, k )
             !$OMP DO
             DO  k = nzb_x, nzt_x
                DO  j = nys_x, nyn_x


                   x_out(0) = CMPLX( ar(0,j,k), 0.0_wp, KIND=wp )
                   DO  i = 1, (nx+1)/2 - 1
                      x_out(i) = CMPLX( ar(i,j,k), ar(nx+1-i,j,k), KIND=wp )
                   ENDDO
                   x_out((nx+1)/2) = CMPLX( ar((nx+1)/2,j,k), 0.0_wp,       &
                                               KIND=wp )

                   CALL FFTW_EXECUTE_DFT_C2R( plan_xi, x_out, x_in)
                   ar(0:nx,j,k) = x_in(0:nx)

                ENDDO
             ENDDO
             !$OMP END PARALLEL

          ENDIF
#endif
    END SUBROUTINE fft_x

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Fourier-transformation along x-direction.
!> Version for 1D-decomposition.
!> It uses internal algorithms (Singleton or Temperton) or
!> system-specific routines, if they are available
!------------------------------------------------------------------------------!

    SUBROUTINE fft_x_1d( ar, direction )


       IMPLICIT NONE

       CHARACTER (LEN=*) ::  direction  !<

       INTEGER(iwp) ::  i               !<
       INTEGER(iwp) ::  ishape(1)       !<

       LOGICAL ::  forward_fft          !<

       REAL(wp), DIMENSION(0:nx)   ::  ar     !<
       REAL(wp), DIMENSION(0:nx+2) ::  work   !<
       REAL(wp), DIMENSION(nx+2)   ::  work1  !<

       COMPLEX(wp), DIMENSION(:), ALLOCATABLE ::  cwork  !<
       IF ( direction == 'forward' )  THEN
          forward_fft = .TRUE.
       ELSE
          forward_fft = .FALSE.
       ENDIF

      #if defined( __fftw )
          IF ( forward_fft )  THEN

             x_in(0:nx) = ar(0:nx)
             CALL FFTW_EXECUTE_DFT_R2C( plan_xf, x_in, x_out )

             DO  i = 0, (nx+1)/2
                ar(i) = REAL( x_out(i), KIND=wp ) / ( nx+1 )
             ENDDO
             DO  i = 1, (nx+1)/2 - 1
                ar(nx+1-i) = AIMAG( x_out(i) ) / ( nx+1 )
             ENDDO

         ELSE

             x_out(0) = CMPLX( ar(0), 0.0_wp, KIND=wp )
             DO  i = 1, (nx+1)/2 - 1
                x_out(i) = CMPLX( ar(i), ar(nx+1-i), KIND=wp )
             ENDDO
             x_out((nx+1)/2) = CMPLX( ar((nx+1)/2), 0.0_wp, KIND=wp )

             CALL FFTW_EXECUTE_DFT_C2R( plan_xi, x_out, x_in)
             ar(0:nx) = x_in(0:nx)

         ENDIF
#endif

          END SUBROUTINE fft_x_1d

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Fourier-transformation along y-direction.
!> Version for 2D-decomposition.
!> It uses internal algorithms (Singleton or Temperton) or
!> system-specific routines, if they are available.
!>
!> direction:  'forward' or 'backward'
!> ar, ar_tr:  3D data arrays
!>             forward:   ar: before  ar_tr: after transformation
!>             backward:  ar_tr: before  ar: after transfosition
!>
!> In case of non-overlapping transposition/transformation:
!> nxl_y_bound = nxl_y_l = nxl_y
!> nxr_y_bound = nxr_y_l = nxr_y
!>
!> In case of overlapping transposition/transformation
!> - nxl_y_bound  and  nxr_y_bound have the original values of
!>   nxl_y, nxr_y.  ar_tr is dimensioned using these values.
!> - nxl_y_l = nxr_y_r.  ar is dimensioned with these values, so that
!>   transformation is carried out for a 2D-plane only.
!------------------------------------------------------------------------------!

    SUBROUTINE fft_y( ar, direction )


       IMPLICIT NONE

       CHARACTER (LEN=*) ::  direction  !<

       INTEGER(iwp) ::  ierr
       INTEGER(iwp) ::  i            !<
       INTEGER(iwp) ::  j            !<
       INTEGER(iwp) ::  k            !<

       LOGICAL ::  forward_fft  !<


! #if defined( __GPU)
       ! REAL(wp), DEVICE, DIMENSION(0:ny,nxl_y:nxr_y,nzb_y:nzt_y) :: ar
! #else
       REAL(wp), DIMENSION(0:ny,nxl_y:nxr_y,nzb_y:nzt_y) ::  ar
       REAL(wp), DEVICE, DIMENSION(0:ny,nxl_y:nxr_y,nzb_y:nzt_y) :: ar_dev
! #endif

       IF ( direction == 'forward' )  THEN
          forward_fft = .TRUE.
       ELSE
          forward_fft = .FALSE.
       ENDIF

#if defined( __GPU )
       ar_dev = ar
       IF ( forward_fft )  THEN

          !$OMP PARALLEL PRIVATE ( work, i, j, k )
          !$OMP DO
          DO  k = nzb_y, nzt_y
             DO  i = nxl_y, nxr_y

                y_in_dev(0:ny) = ar(0:ny,i,k)
                ierr = cufftExecR2C( plan_yf_dev, y_in_dev, y_out_dev)
                y_out = y_out_dev

                DO  j = 0, (ny+1)/2
                   ar(j,i,k) = REAL( y_out(j), KIND=wp ) / (ny+1)
                ENDDO
                DO  j = 1, (ny+1)/2 - 1
                   ar(ny+1-j,i,k) = AIMAG( y_out(j) ) / (ny+1)
                ENDDO

             ENDDO
          ENDDO
          !$OMP END PARALLEL

       ELSE

          !$OMP PARALLEL PRIVATE ( work, i, j, k )
          !$OMP DO
          DO  k = nzb_y, nzt_y
             DO  i = nxl_y, nxr_y

                y_out(0) = CMPLX( ar(0,i,k), 0.0_wp, KIND=wp )
                DO  j = 1, (ny+1)/2 - 1
                   y_out(j) = CMPLX( ar(j,i,k), ar(ny+1-j,i,k),       &
                                     KIND=wp )
                ENDDO
                y_out((ny+1)/2) = CMPLX( ar((ny+1)/2,i,k), 0.0_wp,       &
                                         KIND=wp )

                y_out_dev = y_out
                ierr = cufftExecC2R( plan_yi_dev, y_out_dev, y_in_dev )
                ar(0:ny,i,k) = y_in_dev(0:ny)

             ENDDO
          ENDDO
          !$OMP END PARALLEL

       ENDIF

#else

          IF ( forward_fft )  THEN

             !$OMP PARALLEL PRIVATE ( work, i, j, k )
             !$OMP DO
             DO  k = nzb_y, nzt_y
                DO  i = nxl_y, nxr_y

                   y_in(0:ny) = ar(0:ny,i,k)
                   CALL FFTW_EXECUTE_DFT_R2C( plan_yf, y_in, y_out )

                   DO  j = 0, (ny+1)/2
                      ar(j,i,k) = REAL( y_out(j), KIND=wp ) / (ny+1)
                   ENDDO
                   DO  j = 1, (ny+1)/2 - 1
                      ar(ny+1-j,i,k) = AIMAG( y_out(j) ) / (ny+1)
                   ENDDO

                ENDDO
             ENDDO
             !$OMP END PARALLEL

          ELSE

             !$OMP PARALLEL PRIVATE ( work, i, j, k )
             !$OMP DO
             DO  k = nzb_y, nzt_y
                DO  i = nxl_y, nxr_y

                   y_out(0) = CMPLX( ar(0,i,k), 0.0_wp, KIND=wp )
                   DO  j = 1, (ny+1)/2 - 1
                      y_out(j) = CMPLX( ar(j,i,k), ar(ny+1-j,i,k),       &
                                        KIND=wp )
                   ENDDO
                   y_out((ny+1)/2) = CMPLX( ar((ny+1)/2,i,k), 0.0_wp,       &
                                            KIND=wp )

                   CALL FFTW_EXECUTE_DFT_C2R( plan_yi, y_out, y_in )
                   ar(0:ny,i,k) = y_in(0:ny)

                ENDDO
             ENDDO
             !$OMP END PARALLEL

          ENDIF
#endif
    END SUBROUTINE fft_y

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Fourier-transformation along y-direction.
!> Version for 1D-decomposition.
!> It uses internal algorithms (Singleton or Temperton) or
!> system-specific routines, if they are available.
!------------------------------------------------------------------------------!

    SUBROUTINE fft_y_1d( ar, direction )


       IMPLICIT NONE

       CHARACTER (LEN=*) ::  direction

       INTEGER(iwp) ::  j          !<
       INTEGER(iwp) ::  jshape(1)  !<

       LOGICAL ::  forward_fft  !<

       REAL(wp), DIMENSION(0:ny)    ::  ar     !<
       REAL(wp), DIMENSION(0:ny+2)  ::  work   !<
       REAL(wp), DIMENSION(ny+2)    ::  work1  !<

       COMPLEX(wp), DIMENSION(:), ALLOCATABLE ::  cwork  !<
             IF ( direction == 'forward' )  THEN
          forward_fft = .TRUE.
       ELSE
          forward_fft = .FALSE.
       ENDIF

      #if defined( __fftw )
          IF ( forward_fft )  THEN

             y_in(0:ny) = ar(0:ny)
             CALL FFTW_EXECUTE_DFT_R2C( plan_yf, y_in, y_out )

             DO  j = 0, (ny+1)/2
                ar(j) = REAL( y_out(j), KIND=wp ) / (ny+1)
             ENDDO
             DO  j = 1, (ny+1)/2 - 1
                ar(ny+1-j) = AIMAG( y_out(j) ) / (ny+1)
             ENDDO

          ELSE

             y_out(0) = CMPLX( ar(0), 0.0_wp, KIND=wp )
             DO  j = 1, (ny+1)/2 - 1
                y_out(j) = CMPLX( ar(j), ar(ny+1-j), KIND=wp )
             ENDDO
             y_out((ny+1)/2) = CMPLX( ar((ny+1)/2), 0.0_wp, KIND=wp )

             CALL FFTW_EXECUTE_DFT_C2R( plan_yi, y_out, y_in )
             ar(0:ny) = y_in(0:ny)

          ENDIF
#endif

          END SUBROUTINE fft_y_1d

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Fourier-transformation along x-direction.
!> Version for 1d domain decomposition
!> using multiple 1D FFT from Math Keisan on NEC or Temperton-algorithm
!> (no singleton-algorithm on NEC because it does not vectorize)
!------------------------------------------------------------------------------!

    SUBROUTINE fft_x_m( ar, direction )


       IMPLICIT NONE

       CHARACTER (LEN=*) ::  direction  !<

       INTEGER(iwp) ::  i     !<
       INTEGER(iwp) ::  k     !<
       INTEGER(iwp) ::  siza  !<
       INTEGER(iwp) ::  sizw  !< required on NEC only

       REAL(wp), DIMENSION(0:nx,nz)       ::  ar     !<
       REAL(wp), DIMENSION(0:nx+3,nz+1)   ::  ai     !<
       REAL(wp), DIMENSION(6*(nx+4),nz+1) ::  work1  !<

       COMPLEX(wp), DIMENSION(:,:), ALLOCATABLE ::  work  !< required on NEC only
    END SUBROUTINE fft_x_m

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Fourier-transformation along y-direction.
!> Version for 1d domain decomposition
!> using multiple 1D FFT from Math Keisan on NEC or Temperton-algorithm
!> (no singleton-algorithm on NEC because it does not vectorize)
!------------------------------------------------------------------------------!

    SUBROUTINE fft_y_m( ar, ny1, direction )


       IMPLICIT NONE

       CHARACTER (LEN=*) ::  direction  !<

       INTEGER(iwp) ::  j     !<
       INTEGER(iwp) ::  k     !<
       INTEGER(iwp) ::  ny1   !<
       INTEGER(iwp) ::  siza  !<
       INTEGER(iwp) ::  sizw  !< required on NEC only

       REAL(wp), DIMENSION(0:ny1,nz)      ::  ar     !<
       REAL(wp), DIMENSION(0:ny+3,nz+1)   ::  ai     !<
       REAL(wp), DIMENSION(6*(ny+4),nz+1) ::  work1  !<

       COMPLEX(wp), DIMENSION(:,:), ALLOCATABLE ::  work !< required on NEC only

    END SUBROUTINE fft_y_m


 END MODULE fft_xy
