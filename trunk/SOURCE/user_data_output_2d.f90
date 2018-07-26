!> @file user_data_output_2d.f90
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
! $Id: user_data_output_2d.f90 3014 2018-05-09 08:42:38Z maronga $
! Bugfix: domain bounds of local_pf corrected
! 
! 3004 2018-04-27 12:33:25Z Giersch
! Further allocation checks implemented (averaged data will be assigned to fill
! values if no allocation happened so far)
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2512 2017-10-04 08:26:59Z raasch
! ghost layer points removed from output array local_pf
! 
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Example code added for accessing output quantities stored on surface-data
! types
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1551 2015-03-03 14:18:16Z maronga
! Replaced nzb and nzt+1 with the new array bounds nzb_do and nzt_do.
! 
! 1320 2014-03-20 08:40:49Z raasch
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
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
!> Resorts the user-defined output quantity with indices (k,j,i) to a
!> temporary array with indices (i,j,k) and sets the grid on which it is defined.
!> Allowed values for grid are "zu" and "zw".
!------------------------------------------------------------------------------!
 SUBROUTINE user_data_output_2d( av, variable, found, grid, local_pf, two_d, nzb_do, nzt_do )
 

    USE indices

    USE kinds

    USE surface_mod

    USE user

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  grid     !< 
    CHARACTER (LEN=*) ::  variable !< 

    INTEGER(iwp) ::  av     !< flag to control data output of instantaneous or time-averaged data
    INTEGER(iwp) ::  i      !< grid index along x-direction
    INTEGER(iwp) ::  j      !< grid index along y-direction
    INTEGER(iwp) ::  k      !< grid index along z-direction
    INTEGER(iwp) ::  m      !< running index surface elements
    INTEGER(iwp) ::  nzb_do !< lower limit of the domain (usually nzb)
    INTEGER(iwp) ::  nzt_do !< upper limit of the domain (usually nzt+1)

    LOGICAL      ::  found !< 
    LOGICAL      ::  two_d !< flag parameter that indicates 2D variables (horizontal cross sections)

    REAL(wp) ::  fill_value = -999.0_wp    !< value for the _FillValue attribute

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf !< 


    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

!
!--    Uncomment and extend the following lines, if necessary.
!--    The arrays for storing the user defined quantities (here u2 and u2_av)
!--    have to be declared and defined by the user!
!--    Sample for user-defined output:
!       CASE ( 'u2_xy', 'u2_xz', 'u2_yz' )
!          IF ( av == 0 )  THEN
!             DO  i = nxl, nxr
!                DO  j = nys, nyn
!                   DO  k = nzb_do, nzt_do
!                      local_pf(i,j,k) = u2(k,j,i)
!                   ENDDO
!                ENDDO
!             ENDDO
!          ELSE
!             IF ( .NOT. ALLOCATED( u2_av ) ) THEN
!                ALLOCATE( u2_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!                u2_av = REAL( fill_value, KIND = wp )
!             ENDIF
!             DO  i = nxl, nxr
!                DO  j = nys, nyn
!                   DO  k = nzb_do, nzt_do
!                      local_pf(i,j,k) = u2_av(k,j,i)
!                   ENDDO
!                ENDDO
!             ENDDO
!          ENDIF
!
!          grid = 'zu'
!
!--    In case two-dimensional surface variables are output, the user
!--    has to access related surface-type. Uncomment and extend following lines
!--    appropriately (example output of vertical surface momentum flux of u-
!--    component). Please note, surface elements can be distributed over 
!--    several data type, depending on their respective surface properties.
!       CASE ( 'usws_xy' )
!          IF ( av == 0 )  THEN
!
!--           Horizontal default-type surfaces
!             DO  m = 1, surf_def_h(0)%ns
!                i = surf_def_h(0)%i(m)
!                j = surf_def_h(0)%j(m)
!                local_pf(i,j,1) = surf_def_h(0)%usws(m)
!             ENDDO
!
!--           Horizontal natural-type surfaces
!             DO  m = 1, surf_lsm_h%ns
!                i = surf_lsm_h%i(m)
!                j = surf_lsm_h%j(m)
!                local_pf(i,j,1) = surf_lsm_h%usws(m)
!             ENDDO
!
!--           Horizontal urban-type surfaces
!             DO  m = 1, surf_usm_h%ns
!                i = surf_usm_h%i(m)
!                j = surf_usm_h%j(m)
!                local_pf(i,j,1) = surf_usm_h%usws(m)
!             ENDDO
!          ENDIF
!
!          grid = 'zu'
!--       


       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT


 END SUBROUTINE user_data_output_2d
