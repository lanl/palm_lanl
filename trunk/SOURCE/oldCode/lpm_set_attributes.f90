!> @file lpm_set_attributes.f90
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
! $Id: lpm_set_attributes.f90 3065 2018-06-12 07:03:02Z Giersch $
! dz was replaced by dzw to allow for right vertical stretching
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2123 2017-01-18 12:34:59Z hoffmann
!
! 2122 2017-01-18 12:22:54Z hoffmann
! DVRP routine removed
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1822 2016-04-07 07:49:42Z hoffmann
! Unused variables removed.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1359 2014-04-11 17:15:14Z hoffmann
! New particle structure integrated. 
! Kind definition added to all floating point numbers.
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
! module interfaces removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 849 2012-03-15 10:35:09Z raasch
! routine renamed: set_particle_attributes -> lpm_set_attributes
!
! 828 2012-02-21 12:00:36Z raasch
! particle feature color renamed class
!
! 271 2009-03-26 00:47:14Z raasch
! Initial version
!
! Description:
! ------------
!> This routine sets certain particle attributes depending on the values that
!> other PALM variables have at the current particle position.
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_set_attributes
 

    USE arrays_3d,                                                             &
        ONLY:  dzw, pt, u, v, w, zu, zw

    USE control_parameters,                                                    &
        ONLY:  u_gtrans, v_gtrans

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE dvrp_variables,                                                        &
        ONLY:  color_interval, dvrp_colortable_entries_prt, dvrpsize_interval, &
               particle_color, particle_dvrpsize

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE indices,                                                               &
        ONLY:  ngp_2dh, nxl, nxr, nyn, nys, nzb, nzt

    USE kinds

    USE particle_attributes,                                                   &
        ONLY:  block_offset, grid_particles, number_of_particles, particles,   &
               prt_count

    USE pegrid

    USE statistics,                                                            &
        ONLY:  sums, sums_l

    IMPLICIT NONE

    INTEGER(iwp) ::  i        !<
    INTEGER(iwp) ::  ip       !<
    INTEGER(iwp) ::  j        !<
    INTEGER(iwp) ::  jp       !<
    INTEGER(iwp) ::  k        !<
    INTEGER(iwp) ::  kp       !<
    INTEGER(iwp) ::  n        !<
    INTEGER(iwp) ::  nb       !<

    INTEGER(iwp), DIMENSION(0:7) ::  start_index !<
    INTEGER(iwp), DIMENSION(0:7) ::  end_index   !<

    REAL(wp)    ::  aa        !<
    REAL(wp)    ::  absuv     !<
    REAL(wp)    ::  bb        !<
    REAL(wp)    ::  cc        !<
    REAL(wp)    ::  dd        !<
    REAL(wp)    ::  gg        !<
    REAL(wp)    ::  height    !<
    REAL(wp)    ::  pt_int    !<
    REAL(wp)    ::  pt_int_l  !<
    REAL(wp)    ::  pt_int_u  !<
    REAL(wp)    ::  u_int_l   !<
    REAL(wp)    ::  u_int_u   !<
    REAL(wp)    ::  v_int_l   !<
    REAL(wp)    ::  v_int_u   !<
    REAL(wp)    ::  w_int     !<
    REAL(wp)    ::  w_int_l   !<
    REAL(wp)    ::  w_int_u   !<
    REAL(wp)    ::  x         !<
    REAL(wp)    ::  y         !<

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  u_int                  !<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  v_int                  !<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  xv                     !<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  yv                     !<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  zv                     !<

    CALL cpu_log( log_point_s(49), 'lpm_set_attributes', 'start' )

!
!-- Set particle color
    IF ( particle_color == 'absuv' )  THEN

!
!--    Set particle color depending on the absolute value of the horizontal
!--    velocity
       DO  ip = nxl, nxr
          DO  jp = nys, nyn
             DO  kp = nzb+1, nzt

                number_of_particles = prt_count(kp,jp,ip)
                particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
                IF ( number_of_particles <= 0 )  CYCLE
                start_index = grid_particles(kp,jp,ip)%start_index
                end_index   = grid_particles(kp,jp,ip)%end_index

                ALLOCATE( u_int(1:number_of_particles), &
                          v_int(1:number_of_particles), &
                          xv(1:number_of_particles),    &
                          yv(1:number_of_particles),    &
                          zv(1:number_of_particles) )

                xv = particles(1:number_of_particles)%x
                yv = particles(1:number_of_particles)%y
                zv = particles(1:number_of_particles)%z

                DO  nb = 0,7

                   i = ip
                   j = jp + block_offset(nb)%j_off
                   k = kp + block_offset(nb)%k_off

                   DO  n = start_index(nb), end_index(nb)
!
!--                   Interpolation of the velocity components in the xy-plane
                      x  = xv(n) + ( 0.5_wp - i ) * dx
                      y  = yv(n) - j * dy
                      aa = x**2          + y**2
                      bb = ( dx - x )**2 + y**2
                      cc = x**2          + ( dy - y )**2
                      dd = ( dx - x )**2 + ( dy - y )**2
                      gg = aa + bb + cc + dd

                      u_int_l = ( ( gg - aa ) * u(k,j,i)   + ( gg - bb ) *     &
                                  u(k,j,i+1) + ( gg - cc ) * u(k,j+1,i) +      &
                                  ( gg - dd ) * u(k,j+1,i+1)                   &
                                ) / ( 3.0_wp * gg ) - u_gtrans

                      IF ( k+1 == nzt+1 )  THEN
                         u_int(n) = u_int_l
                      ELSE
                         u_int_u = ( ( gg - aa ) * u(k+1,j,i)   + ( gg - bb ) *  &
                                     u(k+1,j,i+1) + ( gg - cc ) * u(k+1,j+1,i) + &
                                     ( gg - dd ) * u(k+1,j+1,i+1)                &
                                   ) / ( 3.0_wp * gg ) - u_gtrans
                         u_int(n) = u_int_l + ( zv(n) - zu(k) ) / dzw(k) *      &
                                           ( u_int_u - u_int_l )
                      ENDIF

                   ENDDO

                   i = ip + block_offset(nb)%i_off
                   j = jp
                   k = kp + block_offset(nb)%k_off

                   DO  n = start_index(nb), end_index(nb)
!
!--                   Same procedure for interpolation of the v velocity-component
                      x  = xv(n) - i * dx
                      y  = yv(n) + ( 0.5_wp - j ) * dy
                      aa = x**2          + y**2
                      bb = ( dx - x )**2 + y**2
                      cc = x**2          + ( dy - y )**2
                      dd = ( dx - x )**2 + ( dy - y )**2
                      gg = aa + bb + cc + dd

                      v_int_l = ( ( gg - aa ) * v(k,j,i)   + ( gg - bb ) *     &
                                  v(k,j,i+1) + ( gg - cc ) * v(k,j+1,i) +      &
                                  ( gg - dd ) * v(k,j+1,i+1)                   &
                                ) / ( 3.0_wp * gg ) - v_gtrans

                      IF ( k+1 == nzt+1 )  THEN
                         v_int(n) = v_int_l
                      ELSE
                         v_int_u  = ( ( gg - aa ) * v(k+1,j,i) + ( gg - bb ) *    &
                                      v(k+1,j,i+1) + ( gg - cc ) * v(k+1,j+1,i) + &
                                      ( gg - dd ) * v(k+1,j+1,i+1)                &
                                    ) / ( 3.0_wp * gg ) - v_gtrans
                         v_int(n) = v_int_l + ( zv(n) - zu(k) ) / dzw(k) *      &
                                           ( v_int_u - v_int_l )
                      ENDIF

                   ENDDO

                ENDDO

                DO  n = 1, number_of_particles

                   absuv = SQRT( u_int(n)**2 + v_int(n)**2 )

!
!--                Limit values by the given interval and normalize to
!--                interval [0,1]
                   absuv = MIN( absuv, color_interval(2) )
                   absuv = MAX( absuv, color_interval(1) )

                   absuv = ( absuv - color_interval(1) ) / &
                           ( color_interval(2) - color_interval(1) )

!
!--                Number of available colors is defined in init_dvrp
                   particles(n)%class = 1 + absuv *                            &
                                            ( dvrp_colortable_entries_prt - 1 )

                ENDDO

                DEALLOCATE( u_int, v_int, xv, yv, zv )

             ENDDO
          ENDDO
       ENDDO

    ELSEIF ( particle_color == 'pt*' )  THEN
!
!--    Set particle color depending on the resolved scale temperature
!--    fluctuation.
!--    First, calculate the horizontal average of the potential temperature
!--    (This is also done in flow_statistics, but flow_statistics is called
!--    after this routine.)
       sums_l(:,4,0) = 0.0_wp
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             DO  k = nzb, nzt+1
                sums_l(k,4,0) = sums_l(k,4,0) + pt(k,j,i)
             ENDDO
          ENDDO
       ENDDO

#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( sums_l(nzb,4,0), sums(nzb,4), nzt+2-nzb, &
                           MPI_REAL, MPI_SUM, comm2d, ierr )
#else
       sums(:,4) = sums_l(:,4,0)
#endif
       sums(:,4) = sums(:,4) / ngp_2dh(0)

       DO  ip = nxl, nxr
          DO  jp = nys, nyn
             DO  kp = nzb+1, nzt

                number_of_particles = prt_count(kp,jp,ip)
                particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
                IF ( number_of_particles <= 0 )  CYCLE
                start_index = grid_particles(kp,jp,ip)%start_index
                end_index   = grid_particles(kp,jp,ip)%end_index

                ALLOCATE( xv(1:number_of_particles), &
                          yv(1:number_of_particles), &
                          zv(1:number_of_particles) )

                xv = particles(1:number_of_particles)%x
                yv = particles(1:number_of_particles)%y
                zv = particles(1:number_of_particles)%z

                DO  nb = 0,7

                   i = ip + block_offset(nb)%i_off
                   j = jp + block_offset(nb)%j_off
                   k = kp + block_offset(nb)%k_off

                   DO  n = start_index(nb), end_index(nb)
!
!--                   Interpolate temperature to the current particle position
                      x  = xv(n) - i * dx
                      y  = yv(n) - j * dy
                      aa = x**2          + y**2
                      bb = ( dx - x )**2 + y**2
                      cc = x**2          + ( dy - y )**2
                      dd = ( dx - x )**2 + ( dy - y )**2
                      gg = aa + bb + cc + dd

                      pt_int_l = ( ( gg - aa ) * pt(k,j,i)   + ( gg - bb ) *   &
                                   pt(k,j,i+1) + ( gg - cc ) * pt(k,j+1,i) +   &
                                   ( gg - dd ) * pt(k,j+1,i+1)                 &
                                 ) / ( 3.0_wp * gg ) - sums(k,4)

                      pt_int_u = ( ( gg - aa ) * pt(k+1,j,i) + ( gg - bb ) *     &
                                   pt(k+1,j,i+1) + ( gg - cc ) * pt(k+1,j+1,i) + &
                                   ( gg - dd ) * pt(k+1,j+1,i+1)                 &
                                 ) / ( 3.0_wp * gg ) - sums(k,4)

                      pt_int = pt_int_l + ( zv(n) - zu(k) ) / dzw(k) *          &
                                          ( pt_int_u - pt_int_l )

!
!--                   Limit values by the given interval and normalize to
!--                   interval [0,1]
                      pt_int = MIN( pt_int, color_interval(2) )
                      pt_int = MAX( pt_int, color_interval(1) )

                      pt_int = ( pt_int - color_interval(1) ) /                &
                               ( color_interval(2) - color_interval(1) )

!
!--                   Number of available colors is defined in init_dvrp
                      particles(n)%class = 1 + pt_int *                        &
                                           ( dvrp_colortable_entries_prt - 1 )

                   ENDDO
                ENDDO

                DEALLOCATE( xv, yv, zv )

             ENDDO
          ENDDO
       ENDDO

    ELSEIF ( particle_color == 'z' )  THEN
!
!--    Set particle color depending on the height above the bottom
!--    boundary (z=0)
       DO  ip = nxl, nxr
          DO  jp = nys, nyn
             DO  kp = nzb+1, nzt

                number_of_particles = prt_count(kp,jp,ip)
                particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
                IF ( number_of_particles <= 0 )  CYCLE
                DO  n = 1, number_of_particles

                   height = particles(n)%z
!
!--                Limit values by the given interval and normalize to
!--                interval [0,1]
                   height = MIN( height, color_interval(2) )
                   height = MAX( height, color_interval(1) )

                   height = ( height - color_interval(1) ) / &
                            ( color_interval(2) - color_interval(1) )

!
!--                Number of available colors is defined in init_dvrp
                   particles(n)%class = 1 + height *                           &
                                            ( dvrp_colortable_entries_prt - 1 )

                ENDDO

             ENDDO
          ENDDO
       ENDDO

    ENDIF

    CALL cpu_log( log_point_s(49), 'lpm_set_attributes', 'stop' )


 END SUBROUTINE lpm_set_attributes
