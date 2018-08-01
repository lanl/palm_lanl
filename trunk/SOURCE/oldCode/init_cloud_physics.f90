!> @file init_cloud_physics.f90
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
! ------------------
! $Id: init_cloud_physics.f90 2718 2018-01-02 08:49:38Z maronga $
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
! 1849 2016-04-08 11:33:18Z hoffmann
! bfactor removed and microphysics constants moved microphysics_mod
!
! 1822 2016-04-07 07:49:42Z hoffmann
! icloud_scheme replaced by microphysics_*
!
! 1691 2015-10-26 16:17:44Z maronga
! Removed typo
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1361 2014-04-16 15:17:48Z hoffmann
! sed_qc_const is now calculated here (2-moment microphysics)
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
! 
! 1334 2014-03-25 12:21:40Z heinze
! Bugfix: REAL constants provided with KIND-attribute
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL constants defined as wp-kind
! 
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements, 
! kinds are defined in new module mod_kinds, 
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
!
! 1065 2012-11-22 17:42:36Z hoffmann
! The Courant number of sedimentation can be controlled with c_sedimentation. 
!
! 1053 2012-11-13 17:11:03Z hoffmann
! calculation of the maximum timestep according to the terminal velocity of rain 
! drops in the two moment cloud scheme
!
! calculation of frequently used constants (pirho_l, dpirho_l, schmidt_p_1d3, 
! hyrho) 
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 824 2012-02-17 09:09:57Z raasch
! calculation of b_cond replaced by calculation of bfactor
!
! Revision 1.1  2000/04/13 14:37:22  schroeter
! Initial revision
!
!
! Description:
! ------------
!> Initialization of parameters for handling cloud-physics
!------------------------------------------------------------------------------!
 SUBROUTINE init_cloud_physics
 

    USE arrays_3d,                                                             &
        ONLY:  dzu, hyp, pt_init, zu
        
    USE cloud_parameters,                                                      &
        ONLY:  cp, hyrho, l_d_cp, l_d_r, l_d_rv, l_v, pt_d_t, rho_l, r_d, r_v, &
               t_d_pt
                
    USE control_parameters,                                                    &
        ONLY:  g, message_string, pt_surface, rho_surface, surface_pressure
    
    USE indices,                                                               &
        ONLY:  nzb, nzt
    
    USE kinds

    IMPLICIT NONE

    INTEGER(iwp) ::  k      !<
    
    REAL(wp) ::  t_surface  !<

    ALLOCATE( hyp(nzb:nzt+1), pt_d_t(nzb:nzt+1), t_d_pt(nzb:nzt+1),            &
              hyrho(nzb:nzt+1) )

!
!-- Calculate frequently used parameters
    l_d_cp = l_v / cp
    l_d_r  = l_v / r_d
    l_d_rv = l_v / r_v

!
!-- Calculate:
!-- pt / t : ratio of potential and actual temperature (pt_d_t)
!-- t / pt : ratio of actual and potential temperature (t_d_pt)
!-- p_0(z) : vertical profile of the hydrostatic pressure (hyp)
    t_surface = pt_surface * ( surface_pressure / 1000.0_wp )**0.286_wp
    DO  k = nzb, nzt+1
!
!--    Check temperature in case of too large domain height
       IF ( ( t_surface - g/cp * zu(k) ) < 0.0_wp )  THEN
          WRITE( message_string, * )  'absolute temperature < 0.0 at zu(', k, &
                                      ') = ', zu(k)
          CALL message( 'init_cloud_physics', 'PA0142', 1, 2, 0, 6, 0 )
       ENDIF
       hyp(k)    = surface_pressure * 100.0_wp * &
                   ( (t_surface - g/cp * zu(k)) / t_surface )**(1.0_wp/0.286_wp)
       pt_d_t(k) = ( 100000.0_wp / hyp(k) )**0.286_wp
       t_d_pt(k) = 1.0_wp / pt_d_t(k)
       hyrho(k)  = hyp(k) / ( r_d * t_d_pt(k) * pt_init(k) )       
    ENDDO

!
!-- Compute reference density
    rho_surface = surface_pressure * 100.0_wp / ( r_d * t_surface )


 END SUBROUTINE init_cloud_physics
