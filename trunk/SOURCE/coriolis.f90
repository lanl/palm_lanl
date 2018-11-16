!> @file coriolis.f90
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
! $Id: coriolis.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Forcing implemented, preliminary (MS)
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography concept
!
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC version of subroutine removed
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1873 2016-04-18 14:50:06Z maronga
! Module renamed (removed _mod)
!
!
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
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
! 1257 2013-11-08 15:18:40Z raasch
! openacc loop and loop vector clauses removed
!
! 1128 2013-04-12 06:19:32Z raasch
! loop index bounds in accelerator version replaced by i_left, i_right, j_south,
! j_north
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1015 2012-09-27 09:23:24Z raasch
! accelerator version (*_acc) added
!
! Revision 1.1  1997/08/29 08:57:38  raasch
! Initial revision
!
!
! Description:
! ------------
!> Computation of all Coriolis terms in the equations of motion.
!------------------------------------------------------------------------------!
 MODULE coriolis_mod


    PRIVATE
    PUBLIC coriolis

    INTERFACE coriolis
       MODULE PROCEDURE coriolis
       MODULE PROCEDURE coriolis_ij
    END INTERFACE coriolis

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE coriolis( component )

       USE arrays_3d,                                                          &
           ONLY:  tend, u, ug, v, vg, w

       USE control_parameters,                                                 &
           ONLY:  f, forcing, fs, message_string

       USE indices,                                                            &
           ONLY:  nxl, nxlu, nxr, nyn, nys, nysv, nzb, nzt, wall_flags_0

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  component  !<
       INTEGER(iwp) ::  i          !< running index x direction
       INTEGER(iwp) ::  j          !< running index y direction
       INTEGER(iwp) ::  k          !< running index z direction

       REAL(wp)     ::  flag       !< flag to mask topography
       REAL(wp)     ::  flag_force !< flag to mask large-scale pressure gradient in case larger-scale forcing is applied

       flag_force = MERGE( 0.0_wp, 1.0_wp, forcing )
!
!--    Compute Coriolis terms for the three velocity components
       SELECT CASE ( component )

!
!--       u-component
          CASE ( 1 )
             !$acc kernels
             !$acc loop independent
             DO  i = nxlu, nxr
                !$acc loop independent
                DO  j = nys, nyn
                   !$acc loop independent
                   DO  k = nzb+1, nzt
!
!--                   Predetermine flag to mask topography
                      flag = MERGE( 1.0_wp, 0.0_wp,                            &
                                    BTEST( wall_flags_0(k,j,i), 1 ) )

                      tend(k,j,i) = tend(k,j,i) + f  *    ( 0.25_wp *          &
                                   ( v(k,j,i-1) + v(k,j,i) + v(k,j+1,i-1) +    &
                                     v(k,j+1,i) ) - vg(k) * flag_force         &
                                                          ) * flag           &
                                                - fs *    ( 0.25_wp *          &
                                   ( w(k-1,j,i-1) + w(k-1,j,i) + w(k,j,i-1) +  &
                                     w(k,j,i)   )                              &
                                                          ) * flag
                   ENDDO
                ENDDO
             ENDDO
             !$acc end kernels

!
!--       v-component
          CASE ( 2 )
             !$acc kernels
             !$acc loop independent
             DO  i = nxl, nxr
                !$acc loop independent
                DO  j = nysv, nyn
                   !$acc loop independent
                   DO  k = nzb+1, nzt
!
!--                   Predetermine flag to mask topography
                      flag = MERGE( 1.0_wp, 0.0_wp,                            &
                                    BTEST( wall_flags_0(k,j,i), 2 ) )

                      tend(k,j,i) = tend(k,j,i) - f *     ( 0.25_wp *          &
                                   ( u(k,j-1,i) + u(k,j,i) + u(k,j-1,i+1) +    &
                                     u(k,j,i+1) ) - ug(k) * flag_force         &
                                                          ) * flag
                   ENDDO
                ENDDO
             ENDDO
             !$acc end kernels

!
!--       w-component
          CASE ( 3 )
             !$acc kernels
             !$acc loop independent
             DO  i = nxl, nxr
                !$acc loop independent
                DO  j = nys, nyn
                   !$acc loop independent
                   DO  k = nzb+1, nzt
!
!--                   Predetermine flag to mask topography
                      flag = MERGE( 1.0_wp, 0.0_wp,                            &
                                    BTEST( wall_flags_0(k,j,i), 3 ) )

                      tend(k,j,i) = tend(k,j,i) + fs * 0.25_wp *               &
                                   ( u(k,j,i) + u(k+1,j,i) + u(k,j,i+1) +      &
                                     u(k+1,j,i+1) ) * flag
                   ENDDO
                ENDDO
             ENDDO
             !$acc end kernels

          CASE DEFAULT

             WRITE( message_string, * ) ' wrong component: ', component
             CALL message( 'coriolis', 'PA0173', 1, 2, 0, 6, 0 )

       END SELECT

    END SUBROUTINE coriolis


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE coriolis_ij( i, j, component )

       USE arrays_3d,                                                          &
           ONLY:  tend, u, ug, v, vg, w

       USE control_parameters,                                                 &
           ONLY:  f, forcing, fs, message_string

       USE indices,                                                            &
           ONLY:  nzb, nzt, wall_flags_0

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  component  !<
       INTEGER(iwp) ::  i          !< running index x direction
       INTEGER(iwp) ::  j          !< running index y direction
       INTEGER(iwp) ::  k          !< running index z direction

       REAL(wp)     ::  flag       !< flag to mask topography
       REAL(wp)     ::  flag_force !< flag to mask large-scale pressure gradient in case larger-scale forcing is applied

       flag_force = MERGE( 0.0_wp, 1.0_wp, forcing )
!
!--    Compute Coriolis terms for the three velocity components
       SELECT CASE ( component )

!
!--       u-component
          CASE ( 1 )
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 1 ) )

                tend(k,j,i) = tend(k,j,i) + f  *     ( 0.25_wp *               &
                                ( v(k,j,i-1) + v(k,j,i) + v(k,j+1,i-1) +       &
                                  v(k,j+1,i) ) - vg(k) * flag_force            &
                                                     ) * flag                  &
                                          - fs *     ( 0.25_wp *               &
                                ( w(k-1,j,i-1) + w(k-1,j,i) + w(k,j,i-1) +     &
                                  w(k,j,i)   )       ) * flag
             ENDDO

!
!--       v-component
          CASE ( 2 )
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 2 ) )

                tend(k,j,i) = tend(k,j,i) - f *        ( 0.25_wp *             &
                                ( u(k,j-1,i) + u(k,j,i) + u(k,j-1,i+1) +       &
                                  u(k,j,i+1) ) - ug(k)   * flag_force          &
                                                       ) * flag
             ENDDO

!
!--       w-component
          CASE ( 3 )
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 3 ) )

                tend(k,j,i) = tend(k,j,i) + fs * 0.25_wp *                     &
                                ( u(k,j,i) + u(k+1,j,i) + u(k,j,i+1) +         &
                                  u(k+1,j,i+1) ) * flag
             ENDDO

          CASE DEFAULT

             WRITE( message_string, * ) ' wrong component: ', component
             CALL message( 'coriolis', 'PA0173', 1, 2, 0, 6, 0 )

       END SELECT

    END SUBROUTINE coriolis_ij

 END MODULE coriolis_mod
