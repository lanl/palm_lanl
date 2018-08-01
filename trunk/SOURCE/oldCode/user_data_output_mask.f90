!> @file user_data_output_mask.f90
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
! $Id: user_data_output_mask.f90 2718 2018-01-02 08:49:38Z maronga $
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
! 1320 2014-03-20 08:40:49Z raasch
! kind-parameters added to all INTEGER and REAL declaration statements, 
! kinds are defined in new module kinds, 
! old module precision_kind is removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! Description:
! ------------
!> Resorts the user-defined output quantity with indices (k,j,i) to a
!> temporary array with indices (i,j,k) for masked data output.
!------------------------------------------------------------------------------!
 SUBROUTINE user_data_output_mask( av, variable, found, local_pf )
 

    USE control_parameters
        
    USE indices
    
    USE kinds
    
    USE user

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  variable   !< 

    INTEGER(iwp) ::  av   !< 
    INTEGER(iwp) ::  i    !< 
    INTEGER(iwp) ::  j    !< 
    INTEGER(iwp) ::  k    !< 

    LOGICAL ::  found     !< 

    REAL(wp),                                                                  &
       DIMENSION(mask_size_l(mid,1),mask_size_l(mid,2),mask_size_l(mid,3)) ::  &
          local_pf   !< 


    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

!--    Uncomment and extend the following lines, if necessary.
!--    The arrays for storing the user defined quantities (here u2 and u2_av)
!--    have to be declared and defined by the user!
!--    Sample for user-defined output:
!       CASE ( 'u2' )
!          IF ( av == 0 )  THEN
!            DO  i = 1, mask_size_l(mid,1)
!               DO  j = 1, mask_size_l(mid,2)
!                  DO  k = 1, mask_size_l(mid,3)
!                      local_pf(i,j,k) = u2(mask_k(mid,k),                       &
!                                           mask_j(mid,j),mask_i(mid,i))
!                   ENDDO
!                ENDDO
!             ENDDO
!          ELSE
!            DO  i = 1, mask_size_l(mid,1)
!               DO  j = 1, mask_size_l(mid,2)
!                  DO  k = 1, mask_size_l(mid,3)
!                      local_pf(i,j,k) = u2_av(mask_k(mid,k),                    &
!                                              mask_j(mid,j),mask_i(mid,i))
!                   ENDDO
!                ENDDO
!             ENDDO
!          ENDIF

       CASE DEFAULT
          found = .FALSE.

    END SELECT


 END SUBROUTINE user_data_output_mask
