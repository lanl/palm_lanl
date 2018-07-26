!> @file user_statistics.f90
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
! $Id: user_statistics.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Ajdustments for new topography masking
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1804 2016-04-05 16:30:18Z maronga
! Removed code for parameter file check (__check)
!
! 1783 2016-03-06 18:36:17Z raasch
! netcdf module name changed + related changes
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute 
!
! 1320 2014-03-20 08:40:49Z raasch
! kind-parameters added to all INTEGER and REAL declaration statements, 
! kinds are defined in new module kinds, 
! old module precision_kind is removed,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
!
! 1046 2012-11-09 14:38:45Z maronga
! added preprocessor directive for parameter file check 
! 
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 211 2008-11-11 04:46:24Z raasch
! Former file user_interface.f90 split into one file per subroutine
!
! Description:
! ------------
!> Calculation of user-defined statistics, i.e. horizontally averaged profiles
!> and time series.
!> This routine is called for every statistic region sr defined by the user,
!> but at least for the region "total domain" (sr=0).
!> See section 3.5.4 on how to define, calculate, and output user defined
!> quantities.
!------------------------------------------------------------------------------!
 SUBROUTINE user_statistics( mode, sr, tn )
 

    USE arrays_3d
    
    USE indices
    
    USE kinds
    
    USE netcdf_interface,                                                      &
        ONLY:  dots_max

    USE pegrid
    
    USE statistics
    
    USE user

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  mode   !< 

    INTEGER(iwp) ::  i    !< 
    INTEGER(iwp) ::  j    !< 
    INTEGER(iwp) ::  k    !< 
    INTEGER(iwp) ::  sr   !< 
    INTEGER(iwp) ::  tn   !< 

    REAL(wp),                                                                  &
       DIMENSION(dots_num_palm+1:dots_max) ::                                  &
          ts_value_l   !< 


    IF ( mode == 'profiles' )  THEN

!
!--    Sample on how to calculate horizontally averaged profiles of user-
!--    defined quantities. Each quantity is identified by the index
!--    "pr_palm+#" where "#" is an integer starting from 1. These
!--    user-profile-numbers must also be assigned to the respective strings
!--    given by data_output_pr_user in routine user_check_data_output_pr.
!       !$OMP DO
!       DO  i = nxl, nxr
!          DO  j = nys, nyn
!             DO  k = nzb+1, nzt
!!
!!--             Sample on how to calculate the profile of the resolved-scale
!!--             horizontal momentum flux u*v*
!                sums_l(k,pr_palm+1,tn) = sums_l(k,pr_palm+1,tn) +             &
!                      ( 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ) - hom(k,1,1,sr) ) *&
!                      ( 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ) - hom(k,1,2,sr) )  &
!                                     * rmask(j,i,sr)                          &
!                                     * MERGE( 1.0_wp, 0.0_wp,                 &
!                                              BTEST( wall_flags_0(k,j,i), 0 ) )
!!
!!--             Further profiles can be defined and calculated by increasing
!!--             the second index of array sums_l (replace ... appropriately)
!                sums_l(k,pr_palm+2,tn) = sums_l(k,pr_palm+2,tn) + ...           &
!                                         * rmask(j,i,sr)
!             ENDDO
!          ENDDO
!       ENDDO

    ELSEIF ( mode == 'time_series' )  THEN

!
!--    Sample on how to add values for the user-defined time series quantities.
!--    These have to be defined before in routine user_init. This sample
!--    creates two time series for the absolut values of the horizontal
!--    velocities u and v.
!       ts_value_l = 0.0_wp
!       ts_value_l(dots_num_palm+1) = ABS( u_max )
!       ts_value_l(dots_num_palm+2) = ABS( v_max )
!
!--     Collect / send values to PE0, because only PE0 outputs the time series.
!--     CAUTION: Collection is done by taking the sum over all processors.
!--              You may have to normalize this sum, depending on the quantity
!--              that you like to calculate. For serial runs, nothing has to be
!--              done.
!--     HINT: If the time series value that you are calculating has the same
!--           value on all PEs, you can omit the MPI_ALLREDUCE call and
!--           assign ts_value(dots_num_palm+1:,sr) = ts_value_l directly.
!#if defined( __parallel )
!       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
!       CALL MPI_ALLREDUCE( ts_value_l(dots_num_palm+1),                         &
!                           ts_value(dots_num_palm+1,sr),                        &
!                           dots_max-dots_num_palm, MPI_REAL, MPI_SUM, comm2d,   &
!                           ierr )
!#else
!       ts_value(dots_num_palm+1:,sr) = ts_value_l
!#endif

    ENDIF

 END SUBROUTINE user_statistics

