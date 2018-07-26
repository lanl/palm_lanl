!> @file user_read_restart_data_mod.f90
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
! $Id: user_read_restart_data_mod.f90 2894 2018-03-15 09:17:58Z Giersch $
! Initial revision
!
!
! Description:
! ------------
!> Reads user specific restart data into binary file(s) for restart runs.
!------------------------------------------------------------------------------!
 MODULE user_read_restart_data_mod


    USE user
     

    IMPLICIT NONE


    INTERFACE user_rrd_global
       MODULE PROCEDURE user_rrd_global
    END INTERFACE user_rrd_global

    INTERFACE user_rrd_local
       MODULE PROCEDURE user_rrd_local
    END INTERFACE user_rrd_local


    PUBLIC user_rrd_global, user_rrd_local


 CONTAINS


!-------------
! Description:
! ------------
!> Reading global restart data that has been defined by the user.
!------------------------------------------------------------------------------!
    SUBROUTINE user_rrd_global( found )


       USE control_parameters,                                                 &
           ONLY: length, restart_string


       IMPLICIT NONE

       LOGICAL, INTENT(OUT)  ::  found 


       found = .TRUE.


       SELECT CASE ( restart_string(1:length) )

          CASE ( 'global_paramter' )
!             READ ( 13 )  global_parameter

          CASE DEFAULT
 
             found = .FALSE.

       END SELECT  


    END SUBROUTINE user_rrd_global


! Description:
! ------------
!> Reading processor specific restart data from file(s) that has been defined 
!> by the user.
!> Subdomain index limits on file are given by nxl_on_file, etc.
!> Indices nxlc, etc. indicate the range of gridpoints to be mapped from the
!> subdomain on file (f) to the subdomain of the current PE (c). They have been
!> calculated in routine rrd_local.
!------------------------------------------------------------------------------!
    SUBROUTINE user_rrd_local( i, k, nxlf, nxlc, nxl_on_file, nxrf, nxrc,      &
                               nxr_on_file, nynf, nync, nyn_on_file, nysf,     & 
                               nysc, nys_on_file, tmp_3d, found )
    

       USE control_parameters
           
       USE indices
       
       USE kinds
       
       USE pegrid
 

       IMPLICIT NONE

       INTEGER(iwp) ::  i               !< 
       INTEGER(iwp) ::  k               !< 
       INTEGER(iwp) ::  nxlc            !< 
       INTEGER(iwp) ::  nxlf            !< 
       INTEGER(iwp) ::  nxl_on_file     !< 
       INTEGER(iwp) ::  nxrc            !< 
       INTEGER(iwp) ::  nxrf            !< 
       INTEGER(iwp) ::  nxr_on_file     !< 
       INTEGER(iwp) ::  nync            !< 
       INTEGER(iwp) ::  nynf            !< 
       INTEGER(iwp) ::  nyn_on_file     !< 
       INTEGER(iwp) ::  nysc            !< 
       INTEGER(iwp) ::  nysf            !< 
       INTEGER(iwp) ::  nys_on_file     !< 

       LOGICAL, INTENT(OUT)  ::  found 

       REAL(wp), DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_3d   !< 

!
!-- Here the reading of user-defined restart data follows:
!-- Sample for user-defined output


       found = .TRUE.


          SELECT CASE ( restart_string(1:length) )

             CASE ( 'u2_av' )
!                IF ( .NOT. ALLOCATED( u2_av ) ) THEN
!                     ALLOCATE( u2_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!                ENDIF
!                IF ( k == 1 )  READ ( 13 )  tmp_3d
!                   u2_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =         &
!                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
!
             CASE DEFAULT

                found = .FALSE.

             END SELECT


    END SUBROUTINE user_rrd_local


 END MODULE user_read_restart_data_mod