!> @file swap_timelevel.f90
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
! $Id: swap_timelevel.f90 2817 2018-02-19 16:32:21Z knoop $
! Preliminary gust module interface implemented
!
! 2766 2018-01-22 17:17:47Z kanani
! Removed preprocessor directive __chem
!
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Moved TKE to turbulence_closure_mod (TG)
! Implementation of chemistry module (FK)
!
! 2350 2017-08-15 11:48:26Z kanani
! Bugfix in nopointer version
!
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison'
! includes two more prognostic equations for cloud drop concentration (nc)
! and cloud water content (qc).
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new surface concept
!
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC directives removed
!
! 2011 2016-09-19 17:29:57Z kanani
! Flag urban_surface is now defined in module control_parameters.
!
! 2007 2016-08-24 15:47:17Z kanani
! Added swapping of urban surface model quantities,
! removed redundance for land surface model
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1960 2016-07-12 16:34:24Z suehring
! Separate humidity and passive scalar
!
! 1822 2016-04-07 07:49:42Z hoffmann
! icloud_scheme replaced by microphysics_*
!
! swap_timelevel.f90 1766 2016-02-29 08:37:15Z raasch
! setting the swap level for pmc data transfer
!
! 1747 2016-02-08 12:25:53Z raasch
! explicit loops in nopointer case to omit craypointer option of pgi compiler
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1496 2014-12-02 17:25:50Z maronga
! Added swapping of land surface model quantities
!
! 1374 2014-04-25 12:55:07Z raasch
! bugfix: use-statement for nopointer-case added
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! revision history before 2012 removed,
! 1318 2014-03-17 13:35:16Z raasch
! module interfaces removed
!
! 1115 2013-03-26 18:16:16Z hoffmann
! calculation of qr and nr is restricted to precipitation
!
! 1111 2013-03-08 23:54:10Z raasch
! openACC directives added
!
! 1053 2012-11-13 17:11:03Z hoffmann
! swap of timelevels for nr, qr added
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1032 2012-10-21 13:03:21Z letzel
! save memory by not allocating pt_2 in case of neutral = .T.
!
! 1010 2012-09-20 07:59:54Z raasch
! cpp switch __nopointer added for pointer free version
!
! 1001 2012-09-13 14:08:46Z raasch
! all actions concerning leapfrog scheme removed
!
! Revision 1.1  2000/01/10  10:08:58  10:08:58  raasch (Siegfried Raasch)
! Initial revision
!
!
! Description:
! ------------
!> Swap of timelevels of variables after each timestep
!------------------------------------------------------------------------------!
 SUBROUTINE swap_timelevel


#if defined( __nopointer )
    USE arrays_3d,                                                             &
        ONLY:  e, e_p, nc, nc_p, nr, nr_p, pt, pt_p, q, q_p, qc, qc_p, qr,     &
               qr_p, s, s_p, sa, sa_p, u, u_p, v, v_p, w, w_p
#else
    USE arrays_3d,                                                             &
        ONLY:  e, e_1, e_2, e_p, nc, nc_1, nc_2, nc_p, nr, nr_1, nr_2, nr_p,   &
               pt, pt_1, pt_2, pt_p, q, q_1, q_2, q_p, qc, qc_1, qc_2, qc_p,   &
               qr, qr_1, qr_2, qr_p, s, s_1, s_2, s_p, sa, sa_1, sa_2, sa_p,   &
               u, u_1, u_2, u_p, v, v_1, v_2, v_p, w, w_1, w_2, w_p

#endif

    USE cpulog,                                                                &
        ONLY: cpu_log, log_point

    USE control_parameters,                                                    &
        ONLY:  air_chemistry, cloud_physics,               &
               humidity, land_surface,                                         &
               microphysics_morrison, microphysics_seifert, neutral,    &
               passive_scalar, timestep_count, urban_surface

    USE indices,                                                               &
        ONLY:  nxlg, nxrg, nyng, nysg, nzb, nzt

    USE turbulence_closure_mod,                                                &
        ONLY:  tcm_swap_timelevel



    IMPLICIT NONE

    INTEGER ::  i, j, k     !> loop indices
    INTEGER ::  swap_level  !> swap_level for steering the pmc data transfer

!
!-- Incrementing timestep counter
    timestep_count = timestep_count + 1

!
!-- Swap of variables
#if defined( __nopointer )
    CALL cpu_log( log_point(28), 'swap_timelevel (nop)', 'start' )

    !$acc parallel present( u, v, w, pt, sa ) &
    !$acc present( u_p, v_p, w_p, pt_p, sa_p )
    !$acc loop collapse(3)
    DO  i = nxlg, nxrg
       DO  j = nysg, nyng
          DO  k = nzb, nzt+1
             u(k,j,i)  = u_p(k,j,i)
             v(k,j,i)  = v_p(k,j,i)
             w(k,j,i)  = w_p(k,j,i)
             pt(k,j,i) = pt_p(k,j,i)
             sa(k,j,i) = sa_p(k,j,i)
          ENDDO
       ENDDO
    ENDDO
    !$acc end parallel

    CALL tcm_swap_timelevel ( 0 )

    CALL cpu_log( log_point(28), 'swap_timelevel (nop)', 'stop' )
#else
    CALL cpu_log( log_point(28), 'swap_timelevel', 'start' )

    SELECT CASE ( MOD( timestep_count, 2 ) )

       CASE ( 0 )

          u  => u_1;   u_p  => u_2
          v  => v_1;   v_p  => v_2
          w  => w_1;   w_p  => w_2
          pt => pt_1;  pt_p => pt_2
          sa => sa_1;  sa_p => sa_2
         IF ( passive_scalar )  THEN
             s => s_1;    s_p => s_2
          ENDIF

          swap_level = 1

       CASE ( 1 )

          u  => u_2;   u_p  => u_1
          v  => v_2;   v_p  => v_1
          w  => w_2;   w_p  => w_1
          IF ( .NOT. neutral )  THEN
             pt => pt_2;  pt_p => pt_1
          ENDIF
          sa => sa_2;  sa_p => sa_1
          IF ( humidity )  THEN
             q => q_2;    q_p => q_1
             IF ( cloud_physics  .AND.  microphysics_morrison )  THEN
                qc => qc_2;    qc_p => qc_1
                nc => nc_2;    nc_p => nc_1
             ENDIF
             IF ( cloud_physics  .AND.  microphysics_seifert )  THEN
                qr => qr_2;    qr_p => qr_1
                nr => nr_2;    nr_p => nr_1
             ENDIF
          ENDIF
          IF ( passive_scalar )  THEN
             s => s_2;    s_p => s_1
          ENDIF

          swap_level = 2

    END SELECT

    CALL tcm_swap_timelevel ( MOD( timestep_count, 2) )

    CALL cpu_log( log_point(28), 'swap_timelevel', 'stop' )
#endif

 END SUBROUTINE swap_timelevel


