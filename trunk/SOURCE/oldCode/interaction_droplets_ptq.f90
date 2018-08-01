!> @file interaction_droplets_ptq.f90
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
! $Id: interaction_droplets_ptq.f90 3040 2018-05-25 10:22:08Z schwenkel $
! Changed the name specific humidity to mixing ratio
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography concept
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
! 
! 1845 2016-04-08 08:29:13Z raasch
! nzb_2d replaced by nzb_s_inner
!
! 1822 2016-04-07 07:49:42Z hoffmann
! Unused variables removed.
!
! 1779 2016-03-03 08:01:28Z raasch
! bugfix: module procedure names shortened to avoid Intel compiler warnings
! about too long names
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
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
! 799 2011-12-21 17:48:03Z franke
! Bugfix: pt_d_t(k) was missing in calculation of pt_p
!
! RCS Log replace by Id keyword, revision history cleaned up
!
! Revision 1.1  2005/06/26 19:57:47  raasch
! Initial revision
!
!
! Description:
! ------------
!> Release of latent heat and change of mixing ratio due to condensation /
!> evaporation of droplets.
!------------------------------------------------------------------------------!
 MODULE interaction_droplets_ptq_mod
 

    PRIVATE
    PUBLIC interaction_droplets_ptq

    INTERFACE interaction_droplets_ptq
!
!--    Internal names shortened in order ro avoid Intel compiler messages
!--    about too long names
       MODULE PROCEDURE i_droplets_ptq
       MODULE PROCEDURE i_droplets_ptq_ij
    END INTERFACE interaction_droplets_ptq
 
 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE i_droplets_ptq

       USE arrays_3d,                                                          &
           ONLY:  pt_p, ql_c, q_p
           
       USE cloud_parameters,                                                   &
           ONLY:  l_d_cp, pt_d_t
           
       USE indices,                                                            &
           ONLY:  nxl, nxr, nyn, nys, nzb, nzt, wall_flags_0
           
       USE kinds

       USE pegrid

       IMPLICIT NONE

       INTEGER(iwp) ::  i    !< running index x direction
       INTEGER(iwp) ::  j    !< running index y direction 
       INTEGER(iwp) ::  k    !< running index z direction

       REAL(wp) ::  flag     !< flag to mask topography grid points
 
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

                q_p(k,j,i)  = q_p(k,j,i)  - ql_c(k,j,i) * flag
                pt_p(k,j,i) = pt_p(k,j,i) + l_d_cp * ql_c(k,j,i) * pt_d_t(k)   &
                                                        * flag
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE i_droplets_ptq


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE i_droplets_ptq_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  pt_p, ql_c, q_p

       USE cloud_parameters,                                                   &
           ONLY:  l_d_cp, pt_d_t

       USE indices,                                                            &
           ONLY:  nzb, nzt, wall_flags_0

       USE kinds,                                                              &
           ONLY:  iwp, wp

       USE pegrid

       IMPLICIT NONE

       INTEGER(iwp) ::  i    !< running index x direction
       INTEGER(iwp) ::  j    !< running index y direction 
       INTEGER(iwp) ::  k    !< running index z direction

       REAL(wp) ::  flag     !< flag to mask topography grid points


       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 0 ) )

          q_p(k,j,i)  = q_p(k,j,i)  - ql_c(k,j,i) * flag
          pt_p(k,j,i) = pt_p(k,j,i) + l_d_cp * ql_c(k,j,i) * pt_d_t(k) * flag
       ENDDO

    END SUBROUTINE i_droplets_ptq_ij

 END MODULE interaction_droplets_ptq_mod
