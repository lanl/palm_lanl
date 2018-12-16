!> @file exchange_horiz_2d.f90
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
! $Id: exchange_horiz_2d.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Forcing implemented (MS)
! 
! 2101 2017-01-05 16:42:31Z suehring
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1968 2016-07-18 12:01:49Z suehring
! 2D-INTEGER exchange adopted for different multigrid level
! 
! 1818 2016-04-06 15:53:27Z maronga
! Initial version of purely vertical nesting introduced. 
! 
! 1804 2016-04-05 16:30:18Z maronga
! Removed code for parameter file check (__check)
! 
! 1762 2016-02-25 12:31:13Z hellstea
! Introduction of nested domain feature
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
!
! 1348 2014-03-27 18:01:03Z raasch
! bugfix: bc_lr_cyc and bc_ns_cyc added to ONLY-list
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements, 
! kinds are defined in new module kinds, 
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 841 2012-02-28 12:29:49Z maronga
! Excluded routine from compilation of namelist_file_check
!
! Revision 1.1  1998/01/23 09:58:21  raasch
! Initial revision
!
!
! Description:
! ------------
!> Exchange of lateral (ghost) boundaries (parallel computers) and cyclic
!> boundary conditions, respectively, for 2D-arrays.
!------------------------------------------------------------------------------!
 SUBROUTINE exchange_horiz_2d( ar )
 

    USE control_parameters,                                                    &
        ONLY :  bc_lr_cyc, bc_ns_cyc, force_bound_l, force_bound_n,            &
                force_bound_r, force_bound_s,                                  &
                inflow_l, inflow_n, inflow_r, inflow_s,                        &
                nest_bound_l, nest_bound_n, nest_bound_r, nest_bound_s,        &
                outflow_l, outflow_n, outflow_r, outflow_s
                
    USE cpulog,                                                                &
        ONLY :  cpu_log, log_point_s
        
    USE indices,                                                               &
        ONLY :  nbgp, nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg
        
    USE kinds
    
    USE pegrid

    IMPLICIT NONE


    INTEGER(iwp) :: i  !<
    
    REAL(wp) ::  ar(nysg:nyng,nxlg:nxrg)  !<
    

    CALL cpu_log( log_point_s(13), 'exchange_horiz_2d', 'start' )

#if defined( __parallel )

!
!-- Exchange of lateral boundary values for parallel computers
    IF ( pdims(1) == 1 )  THEN

!
!--    One-dimensional decomposition along y, boundary values can be exchanged
!--    within the PE memory
       ar(:,nxlg:nxl-1) = ar(:,nxr-nbgp+1:nxr)
       ar(:,nxr+1:nxrg) = ar(:,nxl:nxl+nbgp-1)

    ELSE
!
!--    Send left boundary, receive right one

       CALL MPI_SENDRECV( ar(nysg,nxl), 1, type_y, pleft,  0,                 &
                          ar(nysg,nxr+1), 1, type_y, pright, 0,               &
                          comm2d, status, ierr )
!
!--    Send right boundary, receive left one
       CALL MPI_SENDRECV( ar(nysg,nxr+1-nbgp), 1, type_y, pright,  1,         &
                          ar(nysg,nxlg), 1, type_y, pleft,   1,               &
                          comm2d, status, ierr )
                          
     
    ENDIF

    IF ( pdims(2) == 1 )  THEN
!
!--    One-dimensional decomposition along x, boundary values can be exchanged
!--    within the PE memory
       ar(nysg:nys-1,:) = ar(nyn-nbgp+1:nyn,:)
       ar(nyn+1:nyng,:) = ar(nys:nys+nbgp-1,:)

    ELSE
!
!--    Send front boundary, receive rear one

       CALL MPI_SENDRECV( ar(nys,nxlg), 1, type_x, psouth, 0,                 &         
                          ar(nyn+1,nxlg), 1, type_x, pnorth, 0,               &
                          comm2d, status, ierr )
!
!--    Send rear boundary, receive front one
       CALL MPI_SENDRECV( ar(nyn+1-nbgp,nxlg), 1, type_x, pnorth, 1,          &
                          ar(nysg,nxlg), 1, type_x, psouth, 1,                &
                          comm2d, status, ierr )

    ENDIF

#else

!
!-- Lateral boundary conditions in the non-parallel case
    IF ( bc_lr_cyc )  THEN
       ar(:,nxlg:nxl-1) = ar(:,nxr-nbgp+1:nxr)
       ar(:,nxr+1:nxrg) = ar(:,nxl:nxl+nbgp-1)
    ENDIF

    IF ( bc_ns_cyc )  THEN
       ar(nysg:nys-1,:) = ar(nyn-nbgp+1:nyn,:)
       ar(nyn+1:nyng,:) = ar(nys:nys+nbgp-1,:)
    ENDIF


#endif
    CALL cpu_log( log_point_s(13), 'exchange_horiz_2d', 'stop' )

 END SUBROUTINE exchange_horiz_2d



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Exchange of lateral (ghost) boundaries (parallel computers) and cyclic
!> boundary conditions, respectively, for 2D integer arrays.
!------------------------------------------------------------------------------!
 
 SUBROUTINE exchange_horiz_2d_int( ar, nys_l, nyn_l, nxl_l, nxr_l, nbgp_local )


    USE control_parameters,                                                    &
        ONLY:  bc_lr_cyc, bc_ns_cyc, grid_level, force_bound_l, force_bound_n, &
               force_bound_r, force_bound_s, nest_bound_l, nest_bound_n,       &
               nest_bound_r, nest_bound_s
        
    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s
               
    USE kinds
    
    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i           !< dummy index to zero-gradient conditions at in/outflow boundaries
    INTEGER(iwp) ::  nxl_l       !< local index bound at current grid level, left side
    INTEGER(iwp) ::  nxr_l       !< local index bound at current grid level, right side
    INTEGER(iwp) ::  nyn_l       !< local index bound at current grid level, north side
    INTEGER(iwp) ::  nys_l       !< local index bound at current grid level, south side
    INTEGER(iwp) ::  nbgp_local  !< number of ghost layers to be exchanged

    INTEGER(iwp), DIMENSION(nys_l-nbgp_local:nyn_l+nbgp_local,                 &
                            nxl_l-nbgp_local:nxr_l+nbgp_local) ::  ar  !< treated array

    CALL cpu_log( log_point_s(13), 'exchange_horiz_2d', 'start' )

#if defined( __parallel )

!
!-- Exchange of lateral boundary values for parallel computers
    IF ( pdims(1) == 1 )  THEN

!
!--    One-dimensional decomposition along y, boundary values can be exchanged
!--    within the PE memory
       ar(:,nxl_l-nbgp_local:nxl_l-1) = ar(:,nxr_l-nbgp_local+1:nxr_l)
       ar(:,nxr_l+1:nxr_l+nbgp_local) = ar(:,nxl_l:nxl_l+nbgp_local-1)

    ELSE
!
!--    Send left boundary, receive right one
       CALL MPI_SENDRECV( ar(nys_l-nbgp_local,nxl_l),   1,                     &
                          type_y_int(grid_level), pleft,  0,                   &
                          ar(nys_l-nbgp_local,nxr_l+1), 1,                     &
                          type_y_int(grid_level), pright, 0,                   &
                          comm2d, status, ierr )
!
!--    Send right boundary, receive left one
       CALL MPI_SENDRECV( ar(nys_l-nbgp_local,nxr_l+1-nbgp_local), 1,          &
                          type_y_int(grid_level), pright, 1,                   &
                          ar(nys_l-nbgp_local,nxl_l-nbgp_local),   1,          & 
                          type_y_int(grid_level), pleft,  1,                   &
                          comm2d, status, ierr )                         

    ENDIF

    IF ( pdims(2) == 1 )  THEN
!
!--    One-dimensional decomposition along x, boundary values can be exchanged
!--    within the PE memory
       ar(nys_l-nbgp_local:nys_l-1,:) = ar(nyn_l+1-nbgp_local:nyn_l,:)
       ar(nyn_l+1:nyn_l+nbgp_local,:) = ar(nys_l:nys_l-1+nbgp_local,:)


    ELSE
!
!--    Send front boundary, receive rear one
       CALL MPI_SENDRECV( ar(nys_l,nxl_l-nbgp_local),   1,                    &
                          type_x_int(grid_level), psouth, 0,                  &
                          ar(nyn_l+1,nxl_l-nbgp_local), 1,                    &
                          type_x_int(grid_level), pnorth, 0,                  &
                          comm2d, status, ierr )                          

!
!--    Send rear boundary, receive front one
       CALL MPI_SENDRECV( ar(nyn_l+1-nbgp_local,nxl_l-nbgp_local), 1,         &
                          type_x_int(grid_level), pnorth, 1,                  &
                          ar(nys_l-nbgp_local,nxl_l-nbgp_local),   1,         &
                          type_x_int(grid_level), psouth, 1,                  &
                          comm2d, status, ierr )

    ENDIF

#else

!
!-- Lateral boundary conditions in the non-parallel case
    IF ( bc_lr_cyc )  THEN
       ar(:,nxl_l-nbgp_local:nxl_l-1) = ar(:,nxr_l-nbgp_local+1:nxr_l)
       ar(:,nxr_l+1:nxr_l+nbgp_local) = ar(:,nxl_l:nxl_l+nbgp_local-1)
    ENDIF

    IF ( bc_ns_cyc )  THEN
       ar(nys_l-nbgp_local:nys_l-1,:) = ar(nyn_l+1-nbgp_local:nyn_l,:)
       ar(nyn_l+1:nyn_l+nbgp_local,:) = ar(nys_l:nys_l-1+nbgp_local,:)
    ENDIF

#endif
!
!-- Neumann-conditions at inflow/outflow/nested boundaries
    IF ( nest_bound_l  .OR.  force_bound_l )  THEN
       DO  i = nbgp_local, 1, -1
         ar(:,nxl_l-i) = ar(:,nxl_l)
       ENDDO
    ENDIF
    IF ( nest_bound_r  .OR.  force_bound_r )  THEN
       DO  i = 1, nbgp_local
          ar(:,nxr_l+i) = ar(:,nxr_l)
       ENDDO
    ENDIF
    IF ( nest_bound_s  .OR.  force_bound_s )  THEN
       DO  i = nbgp_local, 1, -1
         ar(nys_l-i,:) = ar(nys_l,:)
       ENDDO
    ENDIF
    IF ( nest_bound_n  .OR.  force_bound_n )  THEN
       DO  i = 1, nbgp_local
         ar(nyn_l+i,:) = ar(nyn_l,:)
       ENDDO
    ENDIF

    CALL cpu_log( log_point_s(13), 'exchange_horiz_2d', 'stop' )

 END SUBROUTINE exchange_horiz_2d_int
