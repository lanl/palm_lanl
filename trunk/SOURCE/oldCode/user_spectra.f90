!> @file user_spectra.f90
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
! $Id: user_spectra.f90 3049 2018-05-29 13:52:36Z Giersch $
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
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
! 1960 2016-07-12 16:34:24Z suehring
! Consider additional spectra for passive scalar
! 
! 1833 2016-04-07 14:23:03Z raasch
! spectrum renamed spectra_mod
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable  
!
! 1320 2014-03-20 08:40:49Z raasch
! kind-parameters added to all INTEGER and REAL declaration statements, 
! kinds are defined in new module kinds, 
! old module precision_kind is removed,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 211 2008-11-11 04:46:24Z raasch
! Former file user_interface.f90 split into one file per subroutine
!
! Description:
! ------------
!> Calculation of user-defined spectra.
!> See section 3.5.4 on how to define, calculate, and output user defined
!> quantities.
!------------------------------------------------------------------------------!
 SUBROUTINE user_spectra( mode, m, pr )
 

    USE arrays_3d
    
    USE control_parameters
    
    USE indices
    
    USE kinds
    
    USE spectra_mod
    
    USE statistics
    
    USE user

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  mode

    INTEGER(iwp) ::  i    !< 
    INTEGER(iwp) ::  j    !< 
    INTEGER(iwp) ::  k    !< 
    INTEGER(iwp) ::  m    !< 
    INTEGER(iwp) ::  pr   !< 


!
!-- Sample on how to calculate spectra of user-defined quantities. 
!-- Each quantity is identified by the corresponding user profile index
!-- "pr_palm+#" where "#" is an integer starting from 1. These
!-- user-profile-numbers must also be assigned to the respective strings
!-- given by data_output_pr_user in routine user_check_data_output_pr.
    IF ( mode == 'preprocess' )  THEN

       SELECT CASE ( TRIM( data_output_sp(m) ) )
          
          CASE ( 'u', 'v', 'w', 'pt', 'q', 's' )
!--          Not allowed here since these are the standard quantities used in 
!--          preprocess_spectra.
       
!          CASE ( 'u*v*' )
!             pr = pr_palm+1
!             d(nzb+1:nzt,nys:nyn,nxl:nxr) = ustvst(nzb+1:nzt,nys:nyn,nxl:nxr)
       
          CASE DEFAULT
             message_string = 'Spectra of ' //                                 &
                         TRIM( data_output_sp(m) ) // ' can not be calculated'
             CALL message( 'user_spectra', 'UI0010', 0, 1, 0, 6, 0 )
            
       END SELECT

    ELSEIF ( mode == 'data_output' )  THEN

       SELECT CASE ( TRIM( data_output_sp(m) ) )

          CASE ( 'u', 'v', 'w', 'pt', 'q', 's' )
!--          Not allowed here since these are the standard quantities used in 
!--          data_output_spectra.

!          CASE ( 'u*v*' )
!             pr = 6

          CASE DEFAULT
             message_string = 'Spectra of ' //                                 &
                              TRIM( data_output_sp(m) ) // ' are not defined'
             CALL message( 'user_spectra', 'UI0011', 0, 0, 0, 6, 0 )
             
          END SELECT

    ENDIF

 END SUBROUTINE user_spectra

