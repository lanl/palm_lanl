!> @file tridia_solver_mod.f90
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
! $Id: tridia_solver_mod.f90 2718 2018-01-02 08:49:38Z maronga $
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
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
!
!
! 1815 2016-04-06 13:49:59Z raasch
! cpp-switch intel11 removed
!
! 1808 2016-04-05 19:44:00Z raasch
! test output removed
!
! 1804 2016-04-05 16:30:18Z maronga
! Removed code for parameter file check (__check)
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1406 2014-05-16 13:47:01Z raasch
! bugfix for pgi 14.4: declare create moved after array declaration
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
! 1257 2013-11-08 15:18:40Z raasch
! openacc loop and loop vector clauses removed, declare create moved after
! the FORTRAN declaration statement
!
! 1221 2013-09-10 08:59:13Z raasch
! dummy argument tri in 1d-routines replaced by tri_for_1d because of name
! conflict with arry tri in module arrays_3d
!
! 1216 2013-08-26 09:31:42Z raasch
! +tridia_substi_overlap for handling overlapping fft / transposition
!
! 1212 2013-08-15 08:46:27Z raasch
! Initial revision.
! Routines have been moved to seperate module from former file poisfft to here.
! The tridiagonal matrix coefficients of array tri are calculated only once at
! the beginning, i.e. routine split is called within tridia_init.
!
!
! Description:
! ------------
!> solves the linear system of equations:
!>
!> -(4 pi^2(i^2/(dx^2*nnx^2)+j^2/(dy^2*nny^2))+
!>   1/(dzu(k)*dzw(k))+1/(dzu(k-1)*dzw(k)))*p(i,j,k)+
!> 1/(dzu(k)*dzw(k))*p(i,j,k+1)+1/(dzu(k-1)*dzw(k))*p(i,j,k-1)=d(i,j,k)
!>
!> by using the Thomas algorithm
!------------------------------------------------------------------------------!
 MODULE tridia_solver
 

    USE indices,                                                               &
        ONLY:  nx, ny, nz

    USE kinds

    USE transpose_indices,                                                     &
        ONLY:  nxl_z, nyn_z, nxr_z, nys_z

    IMPLICIT NONE

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ddzuw !< 

    PRIVATE

    INTERFACE tridia_substi
       MODULE PROCEDURE tridia_substi
    END INTERFACE tridia_substi

    INTERFACE tridia_substi_overlap
       MODULE PROCEDURE tridia_substi_overlap
    END INTERFACE tridia_substi_overlap

    PUBLIC  tridia_substi, tridia_substi_overlap, tridia_init, tridia_1dd, &
            tridia_deallocate

 CONTAINS

    subroutine tridia_deallocate

       deallocate(ddzuw)

    end subroutine tridia_deallocate

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE tridia_init

       USE arrays_3d,                                                          &
           ONLY:  ddzu_pres, ddzw, rho_air_zw

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  k !< 

       ALLOCATE( ddzuw(0:nz-1,3) )

       DO  k = 0, nz-1
          ddzuw(k,1) = ddzu_pres(k+1) * ddzw(k+1) * rho_air_zw(k)
          ddzuw(k,2) = ddzu_pres(k+2) * ddzw(k+1) * rho_air_zw(k+1)
          ddzuw(k,3) = -1.0_wp * &
                       ( ddzu_pres(k+2) * ddzw(k+1) * rho_air_zw(k+1) +        &
                         ddzu_pres(k+1) * ddzw(k+1) * rho_air_zw(k) )
       ENDDO
!
!--    Calculate constant coefficients of the tridiagonal matrix
       CALL maketri
       CALL split

    END SUBROUTINE tridia_init


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the i- and j-dependent component of the matrix
!> Provide the constant coefficients of the tridiagonal matrix for solution
!> of the Poisson equation in Fourier space.
!> The coefficients are computed following the method of
!> Schmidt et al. (DFVLR-Mitteilung 84-15), which departs from Stephan
!> Siano's original version by discretizing the Poisson equation,
!> before it is Fourier-transformed.
!------------------------------------------------------------------------------!
    SUBROUTINE maketri


          USE arrays_3d,                                                       &
              ONLY:  tric, rho_air

          USE constants,                                                       &
              ONLY:  pi

          USE control_parameters,                                              &
              ONLY:  ibc_p_b, ibc_p_t

          USE grid_variables,                                                  &
              ONLY:  dx, dy


          USE kinds

          IMPLICIT NONE

          INTEGER(iwp) ::  i    !< 
          INTEGER(iwp) ::  j    !< 
          INTEGER(iwp) ::  k    !< 
          INTEGER(iwp) ::  nnxh !< 
          INTEGER(iwp) ::  nnyh !< 

          REAL(wp)    ::  ll(nxl_z:nxr_z,nys_z:nyn_z) !< 


          nnxh = ( nx + 1 ) / 2
          nnyh = ( ny + 1 ) / 2

          DO  j = nys_z, nyn_z
             DO  i = nxl_z, nxr_z
                IF ( j >= 0  .AND.  j <= nnyh )  THEN
                   IF ( i >= 0  .AND.  i <= nnxh )  THEN
                      ll(i,j) = 2.0_wp * ( 1.0_wp - COS( ( 2.0_wp * pi * i ) / &
                                            REAL( nx+1, KIND=wp ) ) ) / ( dx * dx ) + &
                                2.0_wp * ( 1.0_wp - COS( ( 2.0_wp * pi * j ) / &
                                            REAL( ny+1, KIND=wp ) ) ) / ( dy * dy )
                   ELSE
                      ll(i,j) = 2.0_wp * ( 1.0_wp - COS( ( 2.0_wp * pi * ( nx+1-i ) ) / &
                                            REAL( nx+1, KIND=wp ) ) ) / ( dx * dx ) + &
                                2.0_wp * ( 1.0_wp - COS( ( 2.0_wp * pi * j ) / &
                                            REAL( ny+1, KIND=wp ) ) ) / ( dy * dy )
                   ENDIF
                ELSE
                   IF ( i >= 0  .AND.  i <= nnxh )  THEN
                      ll(i,j) = 2.0_wp * ( 1.0_wp - COS( ( 2.0_wp * pi * i ) / &
                                            REAL( nx+1, KIND=wp ) ) ) / ( dx * dx ) + &
                                2.0_wp * ( 1.0_wp - COS( ( 2.0_wp * pi * ( ny+1-j ) ) / &
                                            REAL( ny+1, KIND=wp ) ) ) / ( dy * dy )
                   ELSE
                      ll(i,j) = 2.0_wp * ( 1.0_wp - COS( ( 2.0_wp * pi * ( nx+1-i ) ) / &
                                            REAL( nx+1, KIND=wp ) ) ) / ( dx * dx ) + &
                                2.0_wp * ( 1.0_wp - COS( ( 2.0_wp * pi * ( ny+1-j ) ) / &
                                            REAL( ny+1, KIND=wp ) ) ) / ( dy * dy )
                   ENDIF
                ENDIF
             ENDDO
          ENDDO

          DO  k = 0, nz-1
             DO  j = nys_z, nyn_z
                DO  i = nxl_z, nxr_z
                   tric(i,j,k) = ddzuw(k,3) - ll(i,j) * rho_air(k+1)
                ENDDO
             ENDDO
          ENDDO

          IF ( ibc_p_b == 1 )  THEN
             DO  j = nys_z, nyn_z
                DO  i = nxl_z, nxr_z
                   tric(i,j,0) = tric(i,j,0) + ddzuw(0,1)
                ENDDO
             ENDDO
          ENDIF
          IF ( ibc_p_t == 1 )  THEN
             DO  j = nys_z, nyn_z
                DO  i = nxl_z, nxr_z
                   tric(i,j,nz-1) = tric(i,j,nz-1) + ddzuw(nz-1,2)
                ENDDO
             ENDDO
          ENDIF

    END SUBROUTINE maketri


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Substitution (Forward and Backward) (Thomas algorithm)
!------------------------------------------------------------------------------!
    SUBROUTINE tridia_substi( ar )


          USE arrays_3d,                                                       & 
              ONLY:  tri

          USE control_parameters,                                              &
              ONLY:  ibc_p_b, ibc_p_t

          USE kinds

          IMPLICIT NONE

          INTEGER(iwp) ::  i !< 
          INTEGER(iwp) ::  j !< 
          INTEGER(iwp) ::  k !< 

          REAL(wp)     ::  ar(nxl_z:nxr_z,nys_z:nyn_z,1:nz) !< 

          REAL(wp), DIMENSION(nxl_z:nxr_z,nys_z:nyn_z,0:nz-1)   ::  ar1 !< 

!
!--       Forward substitution
          DO  k = 0, nz - 1
             DO  j = nys_z, nyn_z
                DO  i = nxl_z, nxr_z

                   IF ( k == 0 )  THEN
                      ar1(i,j,k) = ar(i,j,k+1)
                   ELSE
                      ar1(i,j,k) = ar(i,j,k+1) - tri(i,j,k,2) * ar1(i,j,k-1)
                   ENDIF

                ENDDO
             ENDDO
          ENDDO

!
!--       Backward substitution
!--       Note, the 1.0E-20 in the denominator is due to avoid divisions
!--       by zero appearing if the pressure bc is set to neumann at the top of
!--       the model domain.
          DO  k = nz-1, 0, -1
             DO  j = nys_z, nyn_z
                DO  i = nxl_z, nxr_z

                   IF ( k == nz-1 )  THEN
                      ar(i,j,k+1) = ar1(i,j,k) / ( tri(i,j,k,1) + 1.0E-20_wp )
                   ELSE
                      ar(i,j,k+1) = ( ar1(i,j,k) - ddzuw(k,2) * ar(i,j,k+2) ) &
                              / tri(i,j,k,1)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

!
!--       Indices i=0, j=0 correspond to horizontally averaged pressure.
!--       The respective values of ar should be zero at all k-levels if
!--       acceleration of horizontally averaged vertical velocity is zero.
          IF ( ibc_p_b == 1  .AND.  ibc_p_t == 1 )  THEN
             IF ( nys_z == 0  .AND.  nxl_z == 0 )  THEN
                DO  k = 1, nz
                   ar(nxl_z,nys_z,k) = 0.0_wp
                ENDDO
             ENDIF
          ENDIF

    END SUBROUTINE tridia_substi


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Substitution (Forward and Backward) (Thomas algorithm)
!------------------------------------------------------------------------------!
    SUBROUTINE tridia_substi_overlap( ar, jj )


          USE arrays_3d,                                                       &
              ONLY:  tri

          USE control_parameters,                                              &
              ONLY:  ibc_p_b, ibc_p_t

          USE kinds

          IMPLICIT NONE

          INTEGER(iwp) ::  i  !< 
          INTEGER(iwp) ::  j  !< 
          INTEGER(iwp) ::  jj !< 
          INTEGER(iwp) ::  k  !< 

          REAL(wp)     ::  ar(nxl_z:nxr_z,nys_z:nyn_z,1:nz) !< 

          REAL(wp), DIMENSION(nxl_z:nxr_z,nys_z:nyn_z,0:nz-1) ::  ar1 !<

!
!--       Forward substitution
          DO  k = 0, nz - 1
             DO  j = nys_z, nyn_z
                DO  i = nxl_z, nxr_z

                   IF ( k == 0 )  THEN
                      ar1(i,j,k) = ar(i,j,k+1)
                   ELSE
                      ar1(i,j,k) = ar(i,j,k+1) - tri(i,jj,k,2) * ar1(i,j,k-1)
                   ENDIF

                ENDDO
             ENDDO
          ENDDO

!
!--       Backward substitution
!--       Note, the 1.0E-20 in the denominator is due to avoid divisions
!--       by zero appearing if the pressure bc is set to neumann at the top of
!--       the model domain.
          DO  k = nz-1, 0, -1
             DO  j = nys_z, nyn_z
                DO  i = nxl_z, nxr_z

                   IF ( k == nz-1 )  THEN
                      ar(i,j,k+1) = ar1(i,j,k) / ( tri(i,jj,k,1) + 1.0E-20_wp )
                   ELSE
                      ar(i,j,k+1) = ( ar1(i,j,k) - ddzuw(k,2) * ar(i,j,k+2) ) &
                              / tri(i,jj,k,1)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

!
!--       Indices i=0, j=0 correspond to horizontally averaged pressure.
!--       The respective values of ar should be zero at all k-levels if
!--       acceleration of horizontally averaged vertical velocity is zero.
          IF ( ibc_p_b == 1  .AND.  ibc_p_t == 1 )  THEN
             IF ( nys_z == 0  .AND.  nxl_z == 0 )  THEN
                DO  k = 1, nz
                   ar(nxl_z,nys_z,k) = 0.0_wp
                ENDDO
             ENDIF
          ENDIF

    END SUBROUTINE tridia_substi_overlap


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Splitting of the tridiagonal matrix (Thomas algorithm)
!------------------------------------------------------------------------------!
    SUBROUTINE split


          USE arrays_3d,                                                       & 
              ONLY:  tri, tric

          USE kinds

          IMPLICIT NONE

          INTEGER(iwp) ::  i !< 
          INTEGER(iwp) ::  j !< 
          INTEGER(iwp) ::  k !< 
!
!--       Splitting
          DO  j = nys_z, nyn_z
             DO  i = nxl_z, nxr_z
                tri(i,j,0,1) = tric(i,j,0)
             ENDDO
          ENDDO

          DO  k = 1, nz-1
             DO  j = nys_z, nyn_z
                DO  i = nxl_z, nxr_z
                   tri(i,j,k,2) = ddzuw(k,1) / tri(i,j,k-1,1)
                   tri(i,j,k,1) = tric(i,j,k) - ddzuw(k-1,2) * tri(i,j,k,2)
                ENDDO
             ENDDO
          ENDDO

    END SUBROUTINE split


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Solves the linear system of equations for a 1d-decomposition along x (see
!> tridia)
!>
!> @attention when using the intel compilers older than 12.0, array tri must
!>            be passed as an argument to the contained subroutines. Otherwise
!>            addres faults will occur. This feature can be activated with
!>            cpp-switch __intel11
!>            On NEC, tri should not be passed (except for routine substi_1dd)
!>            because this causes very bad performance.
!------------------------------------------------------------------------------!
 
    SUBROUTINE tridia_1dd( ddx2, ddy2, nx, ny, j, ar, tri_for_1d )


       USE arrays_3d,                                                          &
           ONLY:  ddzu_pres, ddzw, rho_air, rho_air_zw

       USE control_parameters,                                                 &
           ONLY:  ibc_p_b, ibc_p_t

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i                  !< 
       INTEGER(iwp) ::  j                  !< 
       INTEGER(iwp) ::  k                  !< 
       INTEGER(iwp) ::  nnyh               !< 
       INTEGER(iwp) ::  nx                 !< 
       INTEGER(iwp) ::  ny                 !< 
       INTEGER(iwp) ::  omp_get_thread_num !< 
       INTEGER(iwp) ::  tn                 !< 

       REAL(wp)     ::  ddx2 !< 
       REAL(wp)     ::  ddy2 !< 

       REAL(wp), DIMENSION(0:nx,1:nz)     ::  ar         !< 
       REAL(wp), DIMENSION(5,0:nx,0:nz-1) ::  tri_for_1d !< 


       nnyh = ( ny + 1 ) / 2

!
!--    Define constant elements of the tridiagonal matrix.
!--    The compiler on SX6 does loop exchange. If 0:nx is a high power of 2,
!--    the exchanged loops create bank conflicts. The following directive
!--    prohibits loop exchange and the loops perform much better.
!CDIR NOLOOPCHG
       DO  k = 0, nz-1
          DO  i = 0,nx
             tri_for_1d(2,i,k) = ddzu_pres(k+1) * ddzw(k+1) * rho_air_zw(k)
             tri_for_1d(3,i,k) = ddzu_pres(k+2) * ddzw(k+1) * rho_air_zw(k+1)
          ENDDO
       ENDDO

       IF ( j <= nnyh )  THEN
          CALL maketri_1dd( j )
       ELSE
          CALL maketri_1dd( ny+1-j )
       ENDIF

       CALL split_1dd
       CALL substi_1dd( ar, tri_for_1d )

    CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> computes the i- and j-dependent component of the matrix
!------------------------------------------------------------------------------!
       SUBROUTINE maketri_1dd( j )

          USE constants,                                                       &
              ONLY:  pi

          USE kinds

          IMPLICIT NONE

          INTEGER(iwp) ::  i    !< 
          INTEGER(iwp) ::  j    !< 
          INTEGER(iwp) ::  k    !< 
          INTEGER(iwp) ::  nnxh !< 

          REAL(wp)     ::  a !< 
          REAL(wp)     ::  c !< 

          REAL(wp), DIMENSION(0:nx) ::  l !< 


          nnxh = ( nx + 1 ) / 2
!
!--       Provide the tridiagonal matrix for solution of the Poisson equation in
!--       Fourier space. The coefficients are computed following the method of
!--       Schmidt et al. (DFVLR-Mitteilung 84-15), which departs from Stephan
!--       Siano's original version by discretizing the Poisson equation,
!--       before it is Fourier-transformed
          DO  i = 0, nx
             IF ( i >= 0 .AND. i <= nnxh ) THEN
                l(i) = 2.0_wp * ( 1.0_wp - COS( ( 2.0_wp * pi * i ) / &
                                   REAL( nx+1, KIND=wp ) ) ) * ddx2 + &
                       2.0_wp * ( 1.0_wp - COS( ( 2.0_wp * pi * j ) / &
                                   REAL( ny+1, KIND=wp ) ) ) * ddy2
             ELSE
                l(i) = 2.0_wp * ( 1.0_wp - COS( ( 2.0_wp * pi * ( nx+1-i ) ) / &
                                   REAL( nx+1, KIND=wp ) ) ) * ddx2 + &
                       2.0_wp * ( 1.0_wp - COS( ( 2.0_wp * pi * j ) / &
                                   REAL( ny+1, KIND=wp ) ) ) * ddy2
             ENDIF
          ENDDO

          DO  k = 0, nz-1
             DO  i = 0, nx
                a = -1.0_wp * ddzu_pres(k+2) * ddzw(k+1) * rho_air_zw(k+1)
                c = -1.0_wp * ddzu_pres(k+1) * ddzw(k+1) * rho_air_zw(k)
                tri_for_1d(1,i,k) = a + c - l(i) * rho_air(k+1)
             ENDDO
          ENDDO
          IF ( ibc_p_b == 1 )  THEN
             DO  i = 0, nx
                tri_for_1d(1,i,0) = tri_for_1d(1,i,0) + tri_for_1d(2,i,0)
             ENDDO
          ENDIF
          IF ( ibc_p_t == 1 )  THEN
             DO  i = 0, nx
                tri_for_1d(1,i,nz-1) = tri_for_1d(1,i,nz-1) + tri_for_1d(3,i,nz-1)
             ENDDO
          ENDIF

       END SUBROUTINE maketri_1dd


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Splitting of the tridiagonal matrix (Thomas algorithm)
!------------------------------------------------------------------------------!
       SUBROUTINE split_1dd

          IMPLICIT NONE

          INTEGER(iwp) ::  i !< 
          INTEGER(iwp) ::  k !< 


!
!--       Splitting
          DO  i = 0, nx
             tri_for_1d(4,i,0) = tri_for_1d(1,i,0)
          ENDDO
          DO  k = 1, nz-1
             DO  i = 0, nx
                tri_for_1d(5,i,k) = tri_for_1d(2,i,k) / tri_for_1d(4,i,k-1)
                tri_for_1d(4,i,k) = tri_for_1d(1,i,k) - tri_for_1d(3,i,k-1) * tri_for_1d(5,i,k)
             ENDDO
          ENDDO

       END SUBROUTINE split_1dd


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Substitution (Forward and Backward) (Thomas algorithm)
!------------------------------------------------------------------------------!
       SUBROUTINE substi_1dd( ar, tri_for_1d )


          IMPLICIT NONE

          INTEGER(iwp) ::  i !< 
          INTEGER(iwp) ::  k !< 

          REAL(wp), DIMENSION(0:nx,nz)       ::  ar         !< 
          REAL(wp), DIMENSION(0:nx,0:nz-1)   ::  ar1        !< 
          REAL(wp), DIMENSION(5,0:nx,0:nz-1) ::  tri_for_1d !< 

!
!--       Forward substitution
          DO  i = 0, nx
             ar1(i,0) = ar(i,1)
          ENDDO
          DO  k = 1, nz-1
             DO  i = 0, nx
                ar1(i,k) = ar(i,k+1) - tri_for_1d(5,i,k) * ar1(i,k-1)
             ENDDO
          ENDDO

!
!--       Backward substitution
!--       Note, the add of 1.0E-20 in the denominator is due to avoid divisions
!--       by zero appearing if the pressure bc is set to neumann at the top of
!--       the model domain.
          DO  i = 0, nx
             ar(i,nz) = ar1(i,nz-1) / ( tri_for_1d(4,i,nz-1) + 1.0E-20_wp )
          ENDDO
          DO  k = nz-2, 0, -1
             DO  i = 0, nx
                ar(i,k+1) = ( ar1(i,k) - tri_for_1d(3,i,k) * ar(i,k+2) ) &
                            / tri_for_1d(4,i,k)
             ENDDO
          ENDDO

!
!--       Indices i=0, j=0 correspond to horizontally averaged pressure.
!--       The respective values of ar should be zero at all k-levels if
!--       acceleration of horizontally averaged vertical velocity is zero.
          IF ( ibc_p_b == 1  .AND.  ibc_p_t == 1 )  THEN
             IF ( j == 0 )  THEN
                DO  k = 1, nz
                   ar(0,k) = 0.0_wp
                ENDDO
             ENDIF
          ENDIF

       END SUBROUTINE substi_1dd

    END SUBROUTINE tridia_1dd


 END MODULE tridia_solver
