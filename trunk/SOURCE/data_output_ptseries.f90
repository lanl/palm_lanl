!> @file data_output_ptseries.f90
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
! $Id: data_output_ptseries.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2312 2017-07-14 20:26:51Z hoffmann
! SGS velocities also possible for curvature_solution_effects = .TRUE. 
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1831 2016-04-07 13:15:51Z hoffmann
! curvature_solution_effects moved to particle_attributes
!
! 1783 2016-03-06 18:36:17Z raasch
! netcdf module name changed + related changes
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1359 2014-04-11 17:15:14Z hoffmann
! New particle structure integrated.
!
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1327 2014-03-21 11:00:16Z raasch
! -netcdf output queries
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1318 2014-03-17 13:35:16Z raasch
! barrier argument removed from cpu_log,
! module interfaces removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 825 2012-02-19 03:03:44Z raasch
! mean/minimum/maximum particle radius added as output quantity,
! particle attributes speed_x|y|z_sgs renamed rvar1|2|3
!
! Revision 1.1  2006/08/04 14:24:18  raasch
! Initial revision
!
!
! Description:
! ------------
!> Output of particle data timeseries in NetCDF format.
!------------------------------------------------------------------------------!
 SUBROUTINE data_output_ptseries


    USE control_parameters,                                                    &
        ONLY:  dopts_time_count, time_since_reference_point

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point

    USE indices,                                                               &
        ONLY: nxl, nxr, nys, nyn, nzb, nzt

    USE kinds

#if defined( __netcdf )
    USE NETCDF
#endif

    USE netcdf_interface,                                                      &
        ONLY:  dopts_num, id_set_pts, id_var_dopts, id_var_time_pts, nc_stat,  &
               netcdf_handle_error

    USE particle_attributes,                                                   &
        ONLY:  grid_particles, number_of_particles, number_of_particle_groups, &
               particles, prt_count

    USE pegrid

    IMPLICIT NONE


    INTEGER(iwp) ::  i    !<
    INTEGER(iwp) ::  inum !<
    INTEGER(iwp) ::  j    !<
    INTEGER(iwp) ::  jg   !<
    INTEGER(iwp) ::  k    !<
    INTEGER(iwp) ::  n    !<

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pts_value   !<
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pts_value_l !<


    CALL cpu_log( log_point(36), 'data_output_ptseries', 'start' )

    IF ( myid == 0 )  THEN
!
!--    Open file for time series output in NetCDF format
       dopts_time_count = dopts_time_count + 1
       CALL check_open( 109 )
#if defined( __netcdf )
!
!--    Update the particle time series time axis
       nc_stat = NF90_PUT_VAR( id_set_pts, id_var_time_pts,      &
                               (/ time_since_reference_point /), &
                               start = (/ dopts_time_count /), count = (/ 1 /) )
       CALL netcdf_handle_error( 'data_output_ptseries', 391 )
#endif

    ENDIF

    ALLOCATE( pts_value(0:number_of_particle_groups,dopts_num), &
              pts_value_l(0:number_of_particle_groups,dopts_num) )

    pts_value_l = 0.0_wp
    pts_value_l(:,16) = 9999999.9_wp    ! for calculation of minimum radius

!
!-- Calculate or collect the particle time series quantities for all particles
!-- and seperately for each particle group (if there is more than one group)
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb, nzt
             number_of_particles = prt_count(k,j,i)
             IF (number_of_particles <= 0)  CYCLE
             particles => grid_particles(k,j,i)%particles(1:number_of_particles)
             DO  n = 1, number_of_particles

                IF ( particles(n)%particle_mask )  THEN  ! Restrict analysis to active particles

                   pts_value_l(0,1)  = pts_value_l(0,1) + 1.0_wp  ! total # of particles
                   pts_value_l(0,2)  = pts_value_l(0,2) +                      &
                          ( particles(n)%x - particles(n)%origin_x )  ! mean x
                   pts_value_l(0,3)  = pts_value_l(0,3) +                      &
                          ( particles(n)%y - particles(n)%origin_y )  ! mean y
                   pts_value_l(0,4)  = pts_value_l(0,4) +                      &
                          ( particles(n)%z - particles(n)%origin_z )  ! mean z
                   pts_value_l(0,5)  = pts_value_l(0,5) + particles(n)%z        ! mean z (absolute)
                   pts_value_l(0,6)  = pts_value_l(0,6) + particles(n)%speed_x  ! mean u
                   pts_value_l(0,7)  = pts_value_l(0,7) + particles(n)%speed_y  ! mean v
                   pts_value_l(0,8)  = pts_value_l(0,8) + particles(n)%speed_z  ! mean w
                   pts_value_l(0,9)  = pts_value_l(0,9)  + particles(n)%rvar1 ! mean sgsu
                   pts_value_l(0,10) = pts_value_l(0,10) + particles(n)%rvar2 ! mean sgsv
                   pts_value_l(0,11) = pts_value_l(0,11) + particles(n)%rvar3 ! mean sgsw
                   IF ( particles(n)%speed_z > 0.0_wp )  THEN
                      pts_value_l(0,12) = pts_value_l(0,12) + 1.0_wp  ! # of upward moving prts
                      pts_value_l(0,13) = pts_value_l(0,13) +                  &
                                              particles(n)%speed_z ! mean w upw.
                   ELSE
                      pts_value_l(0,14) = pts_value_l(0,14) +                  &
                                              particles(n)%speed_z ! mean w down
                   ENDIF
                   pts_value_l(0,15) = pts_value_l(0,15) + particles(n)%radius ! mean rad
                   pts_value_l(0,16) = MIN( pts_value_l(0,16), particles(n)%radius ) ! minrad
                   pts_value_l(0,17) = MAX( pts_value_l(0,17), particles(n)%radius ) ! maxrad
                   pts_value_l(0,18) = pts_value_l(0,18) + 1.0_wp
                   pts_value_l(0,19) = pts_value_l(0,18) + 1.0_wp
!
!--                Repeat the same for the respective particle group
                   IF ( number_of_particle_groups > 1 )  THEN
                      jg = particles(n)%group

                      pts_value_l(jg,1)  = pts_value_l(jg,1) + 1.0_wp
                      pts_value_l(jg,2)  = pts_value_l(jg,2) +                   &
                           ( particles(n)%x - particles(n)%origin_x )
                      pts_value_l(jg,3)  = pts_value_l(jg,3) +                   &
                           ( particles(n)%y - particles(n)%origin_y )
                      pts_value_l(jg,4)  = pts_value_l(jg,4) +                   &
                           ( particles(n)%z - particles(n)%origin_z )
                      pts_value_l(jg,5)  = pts_value_l(jg,5) + particles(n)%z
                      pts_value_l(jg,6)  = pts_value_l(jg,6) + particles(n)%speed_x
                      pts_value_l(jg,7)  = pts_value_l(jg,7) + particles(n)%speed_y
                      pts_value_l(jg,8)  = pts_value_l(jg,8) + particles(n)%speed_z
                      pts_value_l(jg,9)  = pts_value_l(jg,9)  + particles(n)%rvar1
                      pts_value_l(jg,10) = pts_value_l(jg,10) + particles(n)%rvar2
                      pts_value_l(jg,11) = pts_value_l(jg,11) + particles(n)%rvar3
                      IF ( particles(n)%speed_z > 0.0_wp )  THEN
                         pts_value_l(jg,12) = pts_value_l(jg,12) + 1.0_wp
                         pts_value_l(jg,13) = pts_value_l(jg,13) + particles(n)%speed_z
                      ELSE
                         pts_value_l(jg,14) = pts_value_l(jg,14) + particles(n)%speed_z
                      ENDIF
                      pts_value_l(jg,15) = pts_value_l(jg,15) + particles(n)%radius
                      pts_value_l(jg,16) = MIN( pts_value(jg,16), particles(n)%radius )
                      pts_value_l(jg,17) = MAX( pts_value(jg,17), particles(n)%radius )
                      pts_value_l(jg,18) = pts_value_l(jg,18) + 1.0_wp
                      pts_value_l(jg,19) = pts_value_l(jg,19) + 1.0_wp
                   ENDIF

                ENDIF

             ENDDO

          ENDDO
       ENDDO
    ENDDO


#if defined( __parallel )
!
!-- Sum values of the subdomains
    inum = number_of_particle_groups + 1

    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( pts_value_l(0,1), pts_value(0,1), 15*inum, MPI_REAL, &
                        MPI_SUM, comm2d, ierr )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( pts_value_l(0,16), pts_value(0,16), inum, MPI_REAL, &
                        MPI_MIN, comm2d, ierr )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( pts_value_l(0,17), pts_value(0,17), inum, MPI_REAL, &
                        MPI_MAX, comm2d, ierr )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( pts_value_l(0,18), pts_value(0,18), inum, MPI_REAL, &
                        MPI_MAX, comm2d, ierr )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( pts_value_l(0,19), pts_value(0,19), inum, MPI_REAL, &
                        MPI_MIN, comm2d, ierr )
#else
    pts_value(:,1:19) = pts_value_l(:,1:19)
#endif

!
!-- Normalize the above calculated quantities (except min/max values) with the
!-- total number of particles
    IF ( number_of_particle_groups > 1 )  THEN
       inum = number_of_particle_groups
    ELSE
       inum = 0
    ENDIF

    DO  j = 0, inum

       IF ( pts_value(j,1) > 0.0_wp )  THEN

          pts_value(j,2:15) = pts_value(j,2:15) / pts_value(j,1)
          IF ( pts_value(j,12) > 0.0_wp  .AND.  pts_value(j,12) < 1.0_wp )  THEN
             pts_value(j,13) = pts_value(j,13) / pts_value(j,12)
             pts_value(j,14) = pts_value(j,14) / ( 1.0_wp - pts_value(j,12) )
          ELSEIF ( pts_value(j,12) == 0.0_wp )  THEN
             pts_value(j,13) = -1.0_wp
          ELSE
             pts_value(j,14) = -1.0_wp
          ENDIF

       ENDIF

    ENDDO

!
!-- Calculate higher order moments of particle time series quantities,
!-- seperately for each particle group (if there is more than one group)
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb, nzt
             number_of_particles = prt_count(k,j,i)
             IF (number_of_particles <= 0)  CYCLE
             particles => grid_particles(k,j,i)%particles(1:number_of_particles)
             DO  n = 1, number_of_particles

                pts_value_l(0,20) = pts_value_l(0,20) + ( particles(n)%x - &
                                    particles(n)%origin_x - pts_value(0,2) )**2 ! x*2
                pts_value_l(0,21) = pts_value_l(0,21) + ( particles(n)%y - &
                                    particles(n)%origin_y - pts_value(0,3) )**2 ! y*2
                pts_value_l(0,22) = pts_value_l(0,22) + ( particles(n)%z - &
                                    particles(n)%origin_z - pts_value(0,4) )**2 ! z*2
                pts_value_l(0,23) = pts_value_l(0,23) + ( particles(n)%speed_x - &
                                                         pts_value(0,6) )**2   ! u*2
                pts_value_l(0,24) = pts_value_l(0,24) + ( particles(n)%speed_y - &
                                                          pts_value(0,7) )**2   ! v*2
                pts_value_l(0,25) = pts_value_l(0,25) + ( particles(n)%speed_z - &
                                                          pts_value(0,8) )**2   ! w*2
                pts_value_l(0,26) = pts_value_l(0,26) + ( particles(n)%rvar1 - &
                                                          pts_value(0,9) )**2   ! u"2
                pts_value_l(0,27) = pts_value_l(0,27) + ( particles(n)%rvar2 - &
                                                          pts_value(0,10) )**2  ! v"2
                pts_value_l(0,28) = pts_value_l(0,28) + ( particles(n)%rvar3 - &
                                                          pts_value(0,11) )**2  ! w"2
!
!--             Repeat the same for the respective particle group
                IF ( number_of_particle_groups > 1 )  THEN
                   jg = particles(n)%group

                   pts_value_l(jg,20) = pts_value_l(jg,20) + ( particles(n)%x - &
                                       particles(n)%origin_x - pts_value(jg,2) )**2
                   pts_value_l(jg,21) = pts_value_l(jg,21) + ( particles(n)%y - &
                                       particles(n)%origin_y - pts_value(jg,3) )**2
                   pts_value_l(jg,22) = pts_value_l(jg,22) + ( particles(n)%z - &
                                       particles(n)%origin_z - pts_value(jg,4) )**2
                   pts_value_l(jg,23) = pts_value_l(jg,23) + ( particles(n)%speed_x - &
                                                             pts_value(jg,6) )**2
                   pts_value_l(jg,24) = pts_value_l(jg,24) + ( particles(n)%speed_y - &
                                                             pts_value(jg,7) )**2
                   pts_value_l(jg,25) = pts_value_l(jg,25) + ( particles(n)%speed_z - &
                                                             pts_value(jg,8) )**2
                   pts_value_l(jg,26) = pts_value_l(jg,26) + ( particles(n)%rvar1 - &
                                                             pts_value(jg,9) )**2
                   pts_value_l(jg,27) = pts_value_l(jg,27) + ( particles(n)%rvar2 - &
                                                             pts_value(jg,10) )**2
                   pts_value_l(jg,28) = pts_value_l(jg,28) + ( particles(n)%rvar3 - &
                                                             pts_value(jg,11) )**2
                ENDIF

             ENDDO
          ENDDO
       ENDDO
    ENDDO

    pts_value_l(0,29) = ( number_of_particles - pts_value(0,1) / numprocs )**2
                                                 ! variance of particle numbers
    IF ( number_of_particle_groups > 1 )  THEN
       DO  j = 1, number_of_particle_groups
          pts_value_l(j,29) = ( pts_value_l(j,1) - &
                                pts_value(j,1) / numprocs )**2
       ENDDO
    ENDIF

#if defined( __parallel )
!
!-- Sum values of the subdomains
    inum = number_of_particle_groups + 1

    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( pts_value_l(0,20), pts_value(0,20), inum*10, MPI_REAL, &
                        MPI_SUM, comm2d, ierr )
#else
    pts_value(:,20:29) = pts_value_l(:,20:29)
#endif

!
!-- Normalize the above calculated quantities with the total number of
!-- particles
    IF ( number_of_particle_groups > 1 )  THEN
       inum = number_of_particle_groups
    ELSE
       inum = 0
    ENDIF

    DO  j = 0, inum

       IF ( pts_value(j,1) > 0.0_wp )  THEN
          pts_value(j,20:28) = pts_value(j,20:28) / pts_value(j,1)
       ENDIF
       pts_value(j,29) = pts_value(j,29) / numprocs

    ENDDO

#if defined( __netcdf )
!
!-- Output particle time series quantities in NetCDF format
    IF ( myid == 0 )  THEN
       DO  j = 0, inum
          DO  i = 1, dopts_num
             nc_stat = NF90_PUT_VAR( id_set_pts, id_var_dopts(i,j),  &
                                     (/ pts_value(j,i) /),           &
                                     start = (/ dopts_time_count /), &
                                     count = (/ 1 /) )
             CALL netcdf_handle_error( 'data_output_ptseries', 392 )
          ENDDO
       ENDDO
    ENDIF
#endif

    DEALLOCATE( pts_value, pts_value_l )

    CALL cpu_log( log_point(36), 'data_output_ptseries', 'stop' )

 END SUBROUTINE data_output_ptseries
