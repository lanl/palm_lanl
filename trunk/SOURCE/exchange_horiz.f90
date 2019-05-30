!> @file exchange_horiz.f90
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
! $Id: exchange_horiz.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! 3D-Integer exchange on multigrid level (MS)
!
! 2298 2017-06-29 09:28:18Z raasch
! sendrecv_in_background related parts removed
!
! 2119 2017-01-17 16:51:50Z raasch
!
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC directives and related code removed
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1804 2016-04-05 16:30:18Z maronga
! Removed code for parameter file check (__check)
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1677 2015-10-02 13:25:23Z boeske
! Added new routine for exchange of three-dimensional integer arrays
!
! 1569 2015-03-12 07:54:38Z raasch
! bugfix in background communication part
!
! 1348 2014-03-27 18:01:03Z raasch
! bugfix: on_device added to ONLY-list
!
! 1344 2014-03-26 17:33:09Z kanani
! Added missing parameters to ONLY-attribute
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
! modifications for asynchronous transfer,
! local variables req, wait_stat are global now, and have been moved to module
! pegrid
!
! 1113 2013-03-10 02:48:14Z raasch
! GPU-porting for single-core (1PE) mode
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 841 2012-02-28 12:29:49Z maronga
! Excluded routine from compilation of namelist_file_check
!
! Revision 1.1  1997/07/24 11:13:29  raasch
! Initial revision
!
!
! Description:
! ------------
!> Exchange of lateral boundary values (parallel computers) and cyclic
!> lateral boundary conditions, respectively.
!------------------------------------------------------------------------------!
 SUBROUTINE exchange_horiz( ar, nbgp_local)


    USE control_parameters,                                                    &
        ONLY:  bc_lr, bc_lr_cyc, bc_ns, bc_ns_cyc, grid_level,                 &
               mg_switch_to_pe0, synchronous_exchange

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE indices,                                                               &
        ONLY:  nxl, nxr, nyn, nys, nzb, nzt

    USE kinds

    USE pegrid

    IMPLICIT NONE


    INTEGER(iwp) ::  i           !<
    INTEGER(iwp) ::  j           !<
    INTEGER(iwp) ::  k           !<
    INTEGER(iwp) ::  nbgp_local  !<

    REAL(wp), DIMENSION(nzb:nzt+1,nys-nbgp_local:nyn+nbgp_local,               &
                        nxl-nbgp_local:nxr+nbgp_local) ::  ar  !<


    CALL cpu_log( log_point_s(2), 'exchange_horiz', 'start' )

#if defined( __parallel )

!
!-- Exchange in x-direction of lateral boundaries
    IF ( pdims(1) == 1  .OR.  mg_switch_to_pe0 )  THEN
!
!--    One-dimensional decomposition along y, boundary values can be exchanged
!--    within the PE memory
       IF ( bc_lr_cyc )  THEN
          ar(:,:,nxl-nbgp_local:nxl-1) = ar(:,:,nxr-nbgp_local+1:nxr)
          ar(:,:,nxr+1:nxr+nbgp_local) = ar(:,:,nxl:nxl+nbgp_local-1)
       ENDIF

    ELSE

       IF ( synchronous_exchange )  THEN
!
!--       Send left boundary, receive right one (synchronous)
          CALL MPI_SENDRECV(                                                   &
              ar(nzb,nys-nbgp_local,nxl),   1, type_yz(grid_level), pleft,  0, &
              ar(nzb,nys-nbgp_local,nxr+1), 1, type_yz(grid_level), pright, 0, &
              comm2d, status, ierr )
!
!--       Send right boundary, receive left one (synchronous)
          CALL MPI_SENDRECV( ar(nzb,nys-nbgp_local,nxr+1-nbgp_local), 1,       &
                             type_yz(grid_level), pright, 1,                   &
                             ar(nzb,nys-nbgp_local,nxl-nbgp_local), 1,         &
                             type_yz(grid_level), pleft,  1,                   &
                             comm2d, status, ierr )

       ELSE

!
!--       Asynchroneous exchange
          IF ( send_receive == 'lr'  .OR.  send_receive == 'al' )  THEN

             req(1:4)  = 0
             req_count = 0
!
!--          Send left boundary, receive right one (asynchronous)
             CALL MPI_ISEND( ar(nzb,nys-nbgp_local,nxl),   1, type_yz(grid_level), &
                             pleft, req_count, comm2d, req(req_count+1), ierr )
             CALL MPI_IRECV( ar(nzb,nys-nbgp_local,nxr+1), 1, type_yz(grid_level), &
                             pright, req_count, comm2d, req(req_count+2), ierr )
!
!--          Send right boundary, receive left one (asynchronous)
             CALL MPI_ISEND( ar(nzb,nys-nbgp_local,nxr+1-nbgp_local), 1,       &
                             type_yz(grid_level), pright, req_count+1, comm2d, &
                             req(req_count+3), ierr )
             CALL MPI_IRECV( ar(nzb,nys-nbgp_local,nxl-nbgp_local),   1,       &
                             type_yz(grid_level), pleft,  req_count+1, comm2d, &
                             req(req_count+4), ierr )

             CALL MPI_WAITALL( 4, req, wait_stat, ierr )

          ENDIF

       ENDIF

    ENDIF


    IF ( pdims(2) == 1  .OR.  mg_switch_to_pe0 )  THEN
!
!--    One-dimensional decomposition along x, boundary values can be exchanged
!--    within the PE memory
       IF ( bc_ns_cyc )  THEN
          ar(:,nys-nbgp_local:nys-1,:) = ar(:,nyn-nbgp_local+1:nyn,:)
          ar(:,nyn+1:nyn+nbgp_local,:) = ar(:,nys:nys+nbgp_local-1,:)
       ENDIF

    ELSE

       IF ( synchronous_exchange )  THEN
!
!--       Send front boundary, receive rear one (synchronous)
          CALL MPI_SENDRECV(                                                   &
              ar(nzb,nys,nxl-nbgp_local),   1, type_xz(grid_level), psouth, 0, &
              ar(nzb,nyn+1,nxl-nbgp_local), 1, type_xz(grid_level), pnorth, 0, &
              comm2d, status, ierr )
!
!--       Send rear boundary, receive front one (synchronous)
          CALL MPI_SENDRECV( ar(nzb,nyn-nbgp_local+1,nxl-nbgp_local), 1,       &
                             type_xz(grid_level), pnorth, 1,                   &
                             ar(nzb,nys-nbgp_local,nxl-nbgp_local),   1,       &
                             type_xz(grid_level), psouth, 1,                   &
                             comm2d, status, ierr )

       ELSE

!
!--       Asynchroneous exchange
          IF ( send_receive == 'ns'  .OR.  send_receive == 'al' )  THEN

             req(1:4)  = 0
             req_count = 0

!
!--          Send front boundary, receive rear one (asynchronous)
             CALL MPI_ISEND( ar(nzb,nys,nxl-nbgp_local),   1, type_xz(grid_level), &
                             psouth, req_count, comm2d, req(req_count+1), ierr )
             CALL MPI_IRECV( ar(nzb,nyn+1,nxl-nbgp_local), 1, type_xz(grid_level), &
                             pnorth, req_count, comm2d, req(req_count+2), ierr )
!
!--          Send rear boundary, receive front one (asynchronous)
             CALL MPI_ISEND( ar(nzb,nyn-nbgp_local+1,nxl-nbgp_local), 1,       &
                             type_xz(grid_level), pnorth, req_count+1, comm2d, &
                             req(req_count+3), ierr )
             CALL MPI_IRECV( ar(nzb,nys-nbgp_local,nxl-nbgp_local),   1,       &
                             type_xz(grid_level), psouth, req_count+1, comm2d, &
                             req(req_count+4), ierr )

             CALL MPI_WAITALL( 4, req, wait_stat, ierr )

          ENDIF

       ENDIF

    ENDIF

#else

!
!-- Lateral boundary conditions in the non-parallel case.
!-- Case dependent, because in GPU mode still not all arrays are on device. This
!-- workaround has to be removed later. Also, since PGI compiler 12.5 has problems
!-- with array syntax, explicit loops are used.
    !$acc parallel present( ar )
    IF  ( bc_lr == 'cyclic' )  THEN
       ! ar(:,:,nxl-nbgp_local:nxl-1) = ar(:,:,nxr-nbgp_local+1:nxr)
       ! ar(:,:,nxr+1:nxr+nbgp_local) = ar(:,:,nxl:nxl+nbgp_local-1)
       !$acc loop collapse(2)
       DO  j = nys-nbgp_local, nyn+nbgp_local
          DO  k = nzb, nzt+1
             !$acc loop seq
             DO  i = 0, nbgp_local-1
                ar(k,j,nxl-nbgp_local+i) = ar(k,j,nxr-nbgp_local+1+i)
                ar(k,j,nxr+1+i)          = ar(k,j,nxl+i)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    !$acc end parallel

    !$acc parallel present( ar )
    IF ( bc_ns == 'cyclic' )  THEN
       ! ar(:,nys-nbgp_local:nys-1,:) = ar(:,nyn-nbgp_local+1:nyn,:)
       ! ar(:,nyn+1:nyn+nbgp_local,:) = ar(:,nys:nys+nbgp_local-1,:)
       !$acc loop collapse(2)
       DO  i = nxl-nbgp_local, nxr+nbgp_local
          DO  k = nzb, nzt+1
             !$acc loop seq
             DO  j = 0, nbgp_local-1
                ar(k,nys-nbgp_local+j,i) = ar(k,nyn-nbgp_local+1+j,i)
                ar(k,nyn+1+j,i)          = ar(k,nys+j,i)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    !$acc end parallel

#endif
    CALL cpu_log( log_point_s(2), 'exchange_horiz', 'stop' )

 END SUBROUTINE exchange_horiz


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
 SUBROUTINE exchange_horiz_int( ar, nys_l, nyn_l, nxl_l, nxr_l, nzt_l, nbgp_local)

    USE control_parameters,                                                    &
        ONLY:  bc_lr, bc_lr_cyc, bc_ns, bc_ns_cyc, grid_level

    USE indices,                                                               &
        ONLY:  nxl, nxr, nyn, nys, nzb, nzt

    USE kinds

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  nxl_l       !< local index bound at current grid level, left side
    INTEGER(iwp) ::  nxr_l       !< local index bound at current grid level, right side
    INTEGER(iwp) ::  nyn_l       !< local index bound at current grid level, north side
    INTEGER(iwp) ::  nys_l       !< local index bound at current grid level, south side
    INTEGER(iwp) ::  nzt_l       !< local index bound at current grid level, top
    INTEGER(iwp) ::  nbgp_local  !< number of ghost points

    INTEGER(iwp), DIMENSION(nzb:nzt_l+1,nys_l-nbgp_local:nyn_l+nbgp_local,     &
                            nxl_l-nbgp_local:nxr_l+nbgp_local) ::  ar  !< treated array


#if defined( __parallel )
    IF ( pdims(1) == 1 )  THEN
!
!--    One-dimensional decomposition along y, boundary values can be exchanged
!--    within the PE memory
       IF ( bc_lr_cyc )  THEN
          ar(:,:,nxl_l-nbgp_local:nxl_l-1) = ar(:,:,nxr_l-nbgp_local+1:nxr_l)
          ar(:,:,nxr_l+1:nxr_l+nbgp_local) = ar(:,:,nxl_l:nxl_l+nbgp_local-1)
       ENDIF
    ELSE
!
!--    Send left boundary, receive right one (synchronous)
       CALL MPI_SENDRECV(                                                          &
           ar(nzb,nys_l-nbgp_local,nxl_l),   1, type_yz_int(grid_level), pleft,  0,&
           ar(nzb,nys_l-nbgp_local,nxr_l+1), 1, type_yz_int(grid_level), pright, 0,&
           comm2d, status, ierr )
!
!--    Send right boundary, receive left one (synchronous)
       CALL MPI_SENDRECV(                                                          &
           ar(nzb,nys_l-nbgp_local,nxr_l+1-nbgp_local), 1, type_yz_int(grid_level),&
           pright, 1,                                                              &
           ar(nzb,nys_l-nbgp_local,nxl_l-nbgp_local),   1, type_yz_int(grid_level),&
           pleft,  1,                                                              &
           comm2d, status, ierr )
    ENDIF


    IF ( pdims(2) == 1 )  THEN
!
!--    One-dimensional decomposition along x, boundary values can be exchanged
!--    within the PE memory
       IF ( bc_ns_cyc )  THEN
          ar(:,nys_l-nbgp_local:nys_l-1,:) = ar(:,nyn_l-nbgp_local+1:nyn_l,:)
          ar(:,nyn_l+1:nyn_l+nbgp_local,:) = ar(:,nys_l:nys_l+nbgp_local-1,:)
       ENDIF

    ELSE

!
!--    Send front boundary, receive rear one (synchronous)
       CALL MPI_SENDRECV(                                                          &
           ar(nzb,nys_l,nxl_l-nbgp_local),   1, type_xz_int(grid_level), psouth, 0,&
           ar(nzb,nyn_l+1,nxl_l-nbgp_local), 1, type_xz_int(grid_level), pnorth, 0,&
           comm2d, status, ierr )
!
!--    Send rear boundary, receive front one (synchronous)
       CALL MPI_SENDRECV( ar(nzb,nyn_l-nbgp_local+1,nxl_l-nbgp_local), 1,          &
                          type_xz_int(grid_level), pnorth, 1,                      &
                          ar(nzb,nys_l-nbgp_local,nxl_l-nbgp_local),   1,          &
                          type_xz_int(grid_level), psouth, 1,                      &
                          comm2d, status, ierr )

    ENDIF

#else

    IF ( bc_lr == 'cyclic' )  THEN
       ar(:,:,nxl_l-nbgp_local:nxl_l-1) = ar(:,:,nxr_l-nbgp_local+1:nxr_l)
       ar(:,:,nxr_l+1:nxr_l+nbgp_local) = ar(:,:,nxl_l:nxl_l+nbgp_local-1)
    ENDIF

    IF ( bc_ns == 'cyclic' )  THEN
       ar(:,nys_l-nbgp_local:nys_l-1,:) = ar(:,nyn_l-nbgp_local+1:nyn_l,:)
       ar(:,nyn_l+1:nyn_l+nbgp_local,:) = ar(:,nys_l:nys_l+nbgp_local-1,:)
    ENDIF

#endif


 END SUBROUTINE exchange_horiz_int
