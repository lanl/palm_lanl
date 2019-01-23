!> @file timestep.f90
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
! $Id: timestep.f90 3084 2018-06-19 15:30:55Z gronemeier $
! limit increase of dt_3d only in case of RANS mode
! 
! 3083 2018-06-19 14:03:12Z gronemeier
! limit dt_3d to be at maximum 2*old_dt; define old_dt at beginning of routine
! Add km/kh_max
! 
! 3049 2018-05-29 13:52:36Z Giersch
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error message revised
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2365 2017-08-21 14:59:59Z kanani
! Vertical grid nesting: Sync fine and coarse grid timestep (SadiqHuq)
! 
! 2258 2017-06-08 07:55:13Z suehring
! Bugfix, add pre-preprocessor directives to enable non-parrallel mode
!
! 2168 2017-03-06 13:08:38Z suehring
!
! 2130 2017-01-24 16:25:39Z raasch
! bugfix: in case of nested runs the stop condition in case of too small
! timesteps is communicated to all parent/child processes,
! some formatting done
!
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC directives and related part of code removed
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1849 2016-04-08 11:33:18Z hoffmann
! Adapted for modularization of microphysics
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1484 2014-10-21 10:53:05Z kanani
! Changes due to new module structure of the plant canopy model:
!   calculations and parameters related to the plant canopy model removed
!   (the limitation of the canopy drag, i.e. that the canopy drag itself should
!   not change the sign of the velocity components, is now assured for in the
!   calculation of the canopy tendency terms in subroutine plant_canopy_model) 
! 
! 1342 2014-03-26 17:04:47Z kanani
! REAL constants defined as wp-kind
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL functions provided with KIND-attribute
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! old module precision_kind is removed,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
! 
! 1257 2013-11-08 15:18:40Z raasch
! openacc porting
! bugfix for calculation of advective timestep in case of vertically stretched
! grids
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1053 2012-11-13 17:11:03Z hoffmann
! timestep is reduced in two-moment cloud scheme according to the maximum 
! terminal velocity of rain drops
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1001 2012-09-13 14:08:46Z raasch
! all actions concerning leapfrog scheme removed
!
! 978 2012-08-09 08:28:32Z fricke
! restriction of the outflow damping layer in the diffusion criterion removed
!
! 866 2012-03-28 06:44:41Z raasch
! bugfix for timestep calculation in case of Galilei transformation,
! special treatment in case of mirror velocity boundary condition removed
!
! Revision 1.1  1997/08/11 06:26:19  raasch
! Initial revision
!
!
! Description:
! ------------
!> Compute the time step under consideration of the FCL and diffusion criterion.
!------------------------------------------------------------------------------!
 SUBROUTINE timestep
 

    USE arrays_3d,                                                             &
        ONLY:  dzu, dzw, kh, km, u, v, w

    USE control_parameters,                                                    &
        ONLY:  cfl_factor, coupling_mode, dt_3d, dt_fixed, dt_max,             &
               old_dt, message_string, stop_dt, timestep_reason

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point

    USE grid_variables,                                                        &
        ONLY:  dx, dx2, dy, dy2

    USE indices,                                                               &
        ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt

    USE interfaces

    USE kinds

    USE pegrid

    USE statistics,                                                            &
        ONLY:  flow_statistics_called, hom, u_max, u_max_ijk, v_max, v_max_ijk,&
               w_max, w_max_ijk

   IMPLICIT NONE

    INTEGER(iwp) ::  i !< 
    INTEGER(iwp) ::  j !< 
    INTEGER(iwp) ::  k !< 
    INTEGER(iwp) ::  km_max_ijk(3) = -1  !< index values (i,j,k) of location where km_max occurs
    INTEGER(iwp) ::  kh_max_ijk(3) = -1  !< index values (i,j,k) of location where kh_max occurs

    LOGICAL ::  stop_dt_local !< local switch for controlling the time stepping

    REAL(wp) ::  div               !< 
    REAL(wp) ::  dt_diff           !< 
    REAL(wp) ::  dt_diff_l         !< 
    REAL(wp) ::  dt_u              !< 
    REAL(wp) ::  dt_u_l            !< 
    REAL(wp) ::  dt_v              !< 
    REAL(wp) ::  dt_v_l            !< 
    REAL(wp) ::  dt_w              !< 
    REAL(wp) ::  dt_w_l            !< 
    REAL(wp) ::  km_max            !< maximum of Km in entire domain
    REAL(wp) ::  kh_max            !< maximum of Kh in entire domain
    REAL(wp) ::  u_max_l           !< 
    REAL(wp) ::  u_min_l           !< 
    REAL(wp) ::  value             !< 
    REAL(wp) ::  v_max_l           !< 
    REAL(wp) ::  v_min_l           !< 
    REAL(wp) ::  w_max_l           !< 
    REAL(wp) ::  w_min_l           !< 
 
    REAL(wp), DIMENSION(3)         ::  reduce      !< 
    REAL(wp), DIMENSION(3)         ::  reduce_l    !< 
    REAL(wp), DIMENSION(nzb+1:nzt) ::  dxyz2_min   !<  


    CALL cpu_log( log_point(12), 'calculate_timestep', 'start' )
!
!--    Save former time step as reference
       old_dt = dt_3d
!
!-- Determine the maxima of the velocity components, including their
!-- grid index positions.
    CALL global_min_max( nzb, nzt+1, nysg, nyng, nxlg, nxrg, u, 'abs', 0.0_wp, &
                         u_max, u_max_ijk )
    CALL global_min_max( nzb, nzt+1, nysg, nyng, nxlg, nxrg, v, 'abs', 0.0_wp, &
                         v_max, v_max_ijk )
    CALL global_min_max( nzb, nzt+1, nysg, nyng, nxlg, nxrg, w, 'abs', 0.0_wp, &
                         w_max, w_max_ijk )

    IF ( .NOT. dt_fixed )  THEN
!
!--    Variable time step:
!--    Calculate the maximum time step according to the CFL-criterion,
!--    individually for each velocity component
       dt_u_l = 999999.9_wp
       dt_v_l = 999999.9_wp
       dt_w_l = 999999.9_wp
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                dt_u_l = MIN( dt_u_l, ( dx     /                               &
                                 ( ABS( u(k,j,i) ) + 1.0E-10_wp ) ) )
                dt_v_l = MIN( dt_v_l, ( dy     /                               &
                                 ( ABS( v(k,j,i) ) + 1.0E-10_wp ) ) )
                dt_w_l = MIN( dt_w_l, ( dzu(k) /                               &
                                 ( ABS( w(k,j,i) )            + 1.0E-10_wp ) ) )
             ENDDO
          ENDDO
       ENDDO

#if defined( __parallel )
       reduce_l(1) = dt_u_l
       reduce_l(2) = dt_v_l
       reduce_l(3) = dt_w_l
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( reduce_l, reduce, 3, MPI_REAL, MPI_MIN, comm2d, ierr )
       dt_u = reduce(1)
       dt_v = reduce(2)
       dt_w = reduce(3)
#else
       dt_u = dt_u_l
       dt_v = dt_v_l
       dt_w = dt_w_l
#endif

!
!--    Compute time step according to the diffusion criterion.
!--    First calculate minimum grid spacing which only depends on index k
       dt_diff_l = 999999.0_wp

       DO  k = nzb+1, nzt
           dxyz2_min(k) = MIN( dx2, dy2, dzw(k)*dzw(k) ) * 0.125_wp
       ENDDO

       !$OMP PARALLEL private(i,j,k,value) reduction(MIN: dt_diff_l)
       !$OMP DO
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                dt_diff_l = MIN( dt_diff_l, dxyz2_min(k) /                     &
                                    ( MAX( kh(k,j,i), km(k,j,i) ) + 1E-20_wp ) )
             ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL
#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( dt_diff_l, dt_diff, 1, MPI_REAL, MPI_MIN, comm2d,   &
                           ierr )
#else
       dt_diff = dt_diff_l
#endif

!
!--    The time step is the minimum of the 3-4 components and the diffusion time
!--    step minus a reduction (cfl_factor) to be on the safe side.
!--    The time step must not exceed the maximum allowed value.
       dt_3d = cfl_factor * MIN( dt_diff, dt_u, dt_v, dt_w )
       dt_3d = MIN( dt_3d, dt_max )
!
!--    Remember the restricting time step criterion for later output.
       IF ( MIN( dt_u, dt_v, dt_w ) < dt_diff )  THEN
          timestep_reason = 'A'
       ELSE
          timestep_reason = 'D'
       ENDIF

!
!--    Set flag if the time step becomes too small.
       IF ( dt_3d < ( 0.00001_wp * dt_max ) )  THEN
          stop_dt = .TRUE.

!
!--       Determine the maxima of the diffusion coefficients, including their
!--       grid index positions.
          CALL global_min_max( nzb, nzt+1, nysg, nyng, nxlg, nxrg, km, 'abs',  &
                               0.0_wp, km_max, km_max_ijk )
          CALL global_min_max( nzb, nzt+1, nysg, nyng, nxlg, nxrg, kh, 'abs',  &
                               0.0_wp, kh_max, kh_max_ijk )

          WRITE( message_string, * ) 'Time step has reached minimum limit.',   &
               '&dt              = ', dt_3d, ' s  Simulation is terminated.',  &
               '&old_dt          = ', old_dt, ' s',                            &
               '&dt_u            = ', dt_u, ' s',                              &
               '&dt_v            = ', dt_v, ' s',                              &
               '&dt_w            = ', dt_w, ' s',                              &
               '&dt_diff         = ', dt_diff, ' s',                           &
               '&u_max           = ', u_max, ' m/s    k=', u_max_ijk(1),       &
               '  j=', u_max_ijk(2), '  i=', u_max_ijk(3),                     &
               '&v_max           = ', v_max, ' m/s    k=', v_max_ijk(1),       &
               '  j=', v_max_ijk(2), '  i=', v_max_ijk(3),                     &
               '&w_max           = ', w_max, ' m/s    k=', w_max_ijk(1),       &
               '  j=', w_max_ijk(2), '  i=', w_max_ijk(3),                     &
               '&km_max          = ', km_max, ' m2/s2  k=', km_max_ijk(1),     &
               '  j=', km_max_ijk(2), '  i=', km_max_ijk(3),                   &
               '&kh_max          = ', kh_max, ' m2/s2  k=', kh_max_ijk(1),     &
                '  j=', kh_max_ijk(2), '  i=', kh_max_ijk(3)
          CALL message( 'timestep', 'PA0312', 0, 1, 0, 6, 0 )
      ENDIF

!
!
!--    Ensure a smooth value (two significant digits) of the timestep.
       div = 1000.0_wp
       DO  WHILE ( dt_3d < div )
          div = div / 10.0_wp
       ENDDO
       dt_3d = NINT( dt_3d * 100.0_wp / div ) * div / 100.0_wp

    ENDIF

    CALL cpu_log( log_point(12), 'calculate_timestep', 'stop' )

 END SUBROUTINE timestep
