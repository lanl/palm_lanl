!> @file sum_up_3d_data.f90
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
! $Id: sum_up_3d_data.f90 3004 2018-04-27 12:33:25Z Giersch $
! prr field added to ONLY-list, prr* case/pr* case/precipitation_rate_av 
! removed, further allocation checks implemented
! 
! 2963 2018-04-12 14:47:44Z suehring
! Introduce index for vegetation/wall, pavement/green-wall and water/window 
! surfaces, for clearer access of surface fraction, albedo, emissivity, etc. .
! 
! 2894 2018-03-15 09:17:58Z Giersch
! Changed comment
! 
! 2817 2018-02-19 16:32:21Z suehring
! Preliminary gust module interface implemented
! 
! 2798 2018-02-09 17:16:39Z suehring
! Consider also default-type surfaces for surface temperature output. 
! 
! 2797 2018-02-08 13:24:35Z suehring
! Enable output of ground-heat flux also at urban surfaces.
! 
! 2790 2018-02-06 11:57:19Z suehring
! Bugfix in summation of surface sensible and latent heat flux
! 
! 2766 2018-01-22 17:17:47Z kanani
! Removed preprocessor directive __chem
! 
! 2743 2018-01-12 16:03:39Z suehring
! In case of natural- and urban-type surfaces output surfaces fluxes in W/m2.
! 
! 2742 2018-01-12 14:59:47Z suehring
! Enable output of surface temperature
! 
! 2735 2018-01-11 12:01:27Z suehring
! output of r_a moved from land-surface to consider also urban-type surfaces
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! - Change in file header (GPL part)
! - Implementation of uv exposure model (FK)
! - output of diss_av, kh_av, km_av (turbulence_closure_mod) (TG)
! - Implementation of chemistry module (FK)
! - Workaround for sum-up usm arrays in case of restart runs, to avoid program 
!   crash (MS)
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
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean and rho_av to rho_ocean_av
! 
! 2024 2016-10-12 16:42:37Z kanani
! Added missing CASE for ssws*
! 
! 2011 2016-09-19 17:29:57Z kanani
! Flag urban_surface is now defined in module control_parameters,
! changed prefix for urban surface model output to "usm_",
! introduced control parameter varnamelength for LEN of trimvar.
! 
! 2007 2016-08-24 15:47:17Z kanani
! Added support for new urban surface model (temporary modifications of 
! SELECT CASE ( ) necessary, see variable trimvar),
! added comments in variable declaration section
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1992 2016-08-12 15:14:59Z suehring
! Bugfix in summation of passive scalar
! 
! 1976 2016-07-27 13:28:04Z maronga
! Radiation actions are now done directly in the respective module
! 
! 1972 2016-07-26 07:52:02Z maronga
! Land surface actions are now done directly in the respective module
! 
! 1960 2016-07-12 16:34:24Z suehring
! Scalar surface flux added
! 
! 1949 2016-06-17 07:19:16Z maronga
! Bugfix: calculation of lai_av, c_veg_av and c_liq_av.
! 
! 1849 2016-04-08 11:33:18Z hoffmann
! precipitation_rate moved to arrays_3d
!
! 1788 2016-03-10 11:01:04Z maronga
! Added z0q and z0q_av
! 
! 1693 2015-10-27 08:35:45Z maronga
! Last revision text corrected
! 
! 1691 2015-10-26 16:17:44Z maronga
! Added output of Obukhov length and radiative heating rates for RRTMG. 
! Corrected output of liquid water path.
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1585 2015-04-30 07:05:52Z maronga
! Adapted for RRTMG
! 
! 1555 2015-03-04 17:44:27Z maronga
! Added output of r_a and r_s
! 
! 1551 2015-03-03 14:18:16Z maronga
! Added support for land surface model and radiation model data.
! 
! 1359 2014-04-11 17:15:14Z hoffmann
! New particle structure integrated. 
!
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
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
! 1318 2014-03-17 13:35:16Z raasch
! barrier argument removed from cpu_log,
! module interfaces removed
!
! 1115 2013-03-26 18:16:16Z hoffmann
! ql is calculated by calc_liquid_water_content
!
! 1053 2012-11-13 17:11:03Z hoffmann
! +nr, prr, qr
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1007 2012-09-19 14:30:36Z franke
! Bugfix in calculation of ql_vp
!
! 978 2012-08-09 08:28:32Z fricke
! +z0h*
!
! Revision 1.1  2006/02/23 12:55:23  raasch
! Initial revision
!
!
! Description:
! ------------
!> Sum-up the values of 3d-arrays. The real averaging is later done in routine
!> average_3d_data. 
!------------------------------------------------------------------------------!
 SUBROUTINE sum_up_3d_data
 

    USE arrays_3d,                                                             &
        ONLY:  dzw, e, heatflux_output_conversion, nc, nr, p, prr, pt,         &
               q, qc, ql, ql_c, ql_v, qr, rho_ocean, s, sa, u, v, vpt, w,      &
               waterflux_output_conversion, alpha_T, beta_S, solar3d

    USE averaging,                                                             &
        ONLY:  diss_av, e_av, ghf_av, kh_av, km_av, lpt_av, lwp_av, nc_av,     &
               nr_av,                                                          &
               ol_av, p_av, pc_av, pr_av, prr_av, pt_av, q_av, qc_av, ql_av,   &
               ql_c_av, ql_v_av, ql_vp_av, qr_av, qsws_av, qv_av, r_a_av,      &
               rho_ocean_av, s_av, sa_av, shf_av, ssws_av, ts_av, tsurf_av,    &
               u_av, us_av, v_av, vpt_av, w_av, z0_av, z0h_av, z0q_av,         &
               alpha_T_av, beta_S_av, shf_sol_av, solar3d_av

    USE control_parameters,                                                    &
        ONLY:  air_chemistry, average_count_3d, cloud_physics, doav, doav_n,   &
               land_surface, rho_surface, urban_surface, uv_exposure,          &
               varnamelength

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point

        USE indices,                                                               &
        ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt 

    USE kinds

    USE surface_mod,                                                           &
        ONLY:  surf_def_h

    USE turbulence_closure_mod,                                                &
        ONLY:  tcm_3d_data_averaging

    IMPLICIT NONE

    INTEGER(iwp) ::  i   !< grid index x direction
    INTEGER(iwp) ::  ii  !< running index
    INTEGER(iwp) ::  j   !< grid index y direction
    INTEGER(iwp) ::  k   !< grid index x direction
    INTEGER(iwp) ::  m   !< running index surface type
    INTEGER(iwp) ::  n   !< 

    REAL(wp)     ::  mean_r !< 
    REAL(wp)     ::  s_r2   !< 
    REAL(wp)     ::  s_r3   !< 

    CHARACTER (LEN=varnamelength) ::  trimvar  !< TRIM of output-variable string


    CALL cpu_log (log_point(34),'sum_up_3d_data','start')

!
!-- Allocate and initialize the summation arrays if called for the very first
!-- time or the first time after average_3d_data has been called 
!-- (some or all of the arrays may have been already allocated
!-- in rrd_local)
    IF ( average_count_3d == 0 )  THEN

       DO  ii = 1, doav_n
!
!--       Temporary solution to account for data output within the new urban 
!--       surface model (urban_surface_mod.f90), see also SELECT CASE ( trimvar )
          trimvar = TRIM( doav(ii) )
       
          SELECT CASE ( trimvar )

             CASE ( 'ghf*' )
                IF ( .NOT. ALLOCATED( ghf_av ) )  THEN
                   ALLOCATE( ghf_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                ghf_av = 0.0_wp

             CASE ( 'e' )
                IF ( .NOT. ALLOCATED( e_av ) )  THEN
                   ALLOCATE( e_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                e_av = 0.0_wp

             CASE ( 'lpt' )
                IF ( .NOT. ALLOCATED( lpt_av ) )  THEN
                   ALLOCATE( lpt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                lpt_av = 0.0_wp

             CASE ( 'lwp*' )
                IF ( .NOT. ALLOCATED( lwp_av ) )  THEN
                   ALLOCATE( lwp_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                lwp_av = 0.0_wp

             CASE ( 'nc' )
                IF ( .NOT. ALLOCATED( nc_av ) )  THEN
                   ALLOCATE( nc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                nc_av = 0.0_wp

             CASE ( 'nr' )
                IF ( .NOT. ALLOCATED( nr_av ) )  THEN
                   ALLOCATE( nr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                nr_av = 0.0_wp

             CASE ( 'ol*' )
                IF ( .NOT. ALLOCATED( ol_av ) )  THEN
                   ALLOCATE( ol_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                ol_av = 0.0_wp

             CASE ( 'p' )
                IF ( .NOT. ALLOCATED( p_av ) )  THEN
                   ALLOCATE( p_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                p_av = 0.0_wp

             CASE ( 'pc' )
                IF ( .NOT. ALLOCATED( pc_av ) )  THEN
                   ALLOCATE( pc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                pc_av = 0.0_wp

             CASE ( 'pr' )
                IF ( .NOT. ALLOCATED( pr_av ) )  THEN
                   ALLOCATE( pr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                pr_av = 0.0_wp

             CASE ( 'prr' )
                IF ( .NOT. ALLOCATED( prr_av ) )  THEN
                   ALLOCATE( prr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                prr_av = 0.0_wp

             CASE ( 'pt' )
                IF ( .NOT. ALLOCATED( pt_av ) )  THEN
                   ALLOCATE( pt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                pt_av = 0.0_wp

             CASE ( 'q' )
                IF ( .NOT. ALLOCATED( q_av ) )  THEN
                   ALLOCATE( q_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                q_av = 0.0_wp

             CASE ( 'qc' )
                IF ( .NOT. ALLOCATED( qc_av ) )  THEN
                   ALLOCATE( qc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                qc_av = 0.0_wp

             CASE ( 'ql' )
                IF ( .NOT. ALLOCATED( ql_av ) )  THEN
                   ALLOCATE( ql_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                ql_av = 0.0_wp

             CASE ( 'ql_c' )
                IF ( .NOT. ALLOCATED( ql_c_av ) )  THEN
                   ALLOCATE( ql_c_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                ql_c_av = 0.0_wp

             CASE ( 'ql_v' )
                IF ( .NOT. ALLOCATED( ql_v_av ) )  THEN
                   ALLOCATE( ql_v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                ql_v_av = 0.0_wp

             CASE ( 'ql_vp' )
                IF ( .NOT. ALLOCATED( ql_vp_av ) )  THEN
                   ALLOCATE( ql_vp_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                ql_vp_av = 0.0_wp

             CASE ( 'qr' )
                IF ( .NOT. ALLOCATED( qr_av ) )  THEN
                   ALLOCATE( qr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                qr_av = 0.0_wp

             CASE ( 'qsws*' )
                IF ( .NOT. ALLOCATED( qsws_av ) )  THEN
                   ALLOCATE( qsws_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                qsws_av = 0.0_wp

             CASE ( 'qv' )
                IF ( .NOT. ALLOCATED( qv_av ) )  THEN
                   ALLOCATE( qv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                qv_av = 0.0_wp

             CASE ( 'r_a*' )
                IF ( .NOT. ALLOCATED( r_a_av ) )  THEN
                   ALLOCATE( r_a_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                r_a_av = 0.0_wp

            CASE ( 'solar3d' )
                IF ( .NOT. ALLOCATED( solar3d_av ) )  THEN
                   ALLOCATE( solar3d_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                solar3d_av = 0.0_wp

             CASE ( 'rho_ocean' )
                IF ( .NOT. ALLOCATED( rho_ocean_av ) )  THEN
                   ALLOCATE( rho_ocean_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rho_ocean_av = 0.0_wp

             CASE ( 'alpha_T' )
                IF ( .NOT. ALLOCATED( alpha_T_av ) )  THEN
                   ALLOCATE( alpha_T_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                alpha_T_av = 0.0_wp

             CASE ( 'beta_S' )
                IF ( .NOT. ALLOCATED( beta_S_av ) )  THEN
                   ALLOCATE( beta_S_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                beta_S_av = 0.0_wp


             CASE ( 's' )
                IF ( .NOT. ALLOCATED( s_av ) )  THEN
                   ALLOCATE( s_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                s_av = 0.0_wp

             CASE ( 'sa' )
                IF ( .NOT. ALLOCATED( sa_av ) )  THEN
                   ALLOCATE( sa_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                sa_av = 0.0_wp

             CASE ( 'shf*' )
                IF ( .NOT. ALLOCATED( shf_av ) )  THEN
                   ALLOCATE( shf_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                shf_av = 0.0_wp
 
             CASE ( 'shf_sol*' )
                IF ( .NOT. ALLOCATED( shf_sol_av ) )  THEN
                   ALLOCATE( shf_sol_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                shf_sol_av = 0.0_wp
                
             CASE ( 'ssws*' )
                IF ( .NOT. ALLOCATED( ssws_av ) )  THEN
                   ALLOCATE( ssws_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                ssws_av = 0.0_wp                

             CASE ( 't*' )
                IF ( .NOT. ALLOCATED( ts_av ) )  THEN
                   ALLOCATE( ts_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                ts_av = 0.0_wp

             CASE ( 'tsurf*' )
                IF ( .NOT. ALLOCATED( tsurf_av ) )  THEN
                   ALLOCATE( tsurf_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                tsurf_av = 0.0_wp

             CASE ( 'u' )
                IF ( .NOT. ALLOCATED( u_av ) )  THEN
                   ALLOCATE( u_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                u_av = 0.0_wp

             CASE ( 'u*' )
                IF ( .NOT. ALLOCATED( us_av ) )  THEN
                   ALLOCATE( us_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                us_av = 0.0_wp

             CASE ( 'v' )
                IF ( .NOT. ALLOCATED( v_av ) )  THEN
                   ALLOCATE( v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                v_av = 0.0_wp

             CASE ( 'vpt' )
                IF ( .NOT. ALLOCATED( vpt_av ) )  THEN
                   ALLOCATE( vpt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                vpt_av = 0.0_wp

             CASE ( 'w' )
                IF ( .NOT. ALLOCATED( w_av ) )  THEN
                   ALLOCATE( w_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                w_av = 0.0_wp

             CASE ( 'z0*' )
                IF ( .NOT. ALLOCATED( z0_av ) )  THEN
                   ALLOCATE( z0_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                z0_av = 0.0_wp

             CASE ( 'z0h*' )
                IF ( .NOT. ALLOCATED( z0h_av ) )  THEN
                   ALLOCATE( z0h_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                z0h_av = 0.0_wp

             CASE ( 'z0q*' )
                IF ( .NOT. ALLOCATED( z0q_av ) )  THEN
                   ALLOCATE( z0q_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                z0q_av = 0.0_wp
!             
!--          Block of urban surface model outputs 

             CASE DEFAULT

!
!--             Turbulence closure module
                CALL tcm_3d_data_averaging( 'allocate', doav(ii) )

          END SELECT

       ENDDO

    ENDIF

!
!-- Loop of all variables to be averaged.
    DO  ii = 1, doav_n
!
!--       Temporary solution to account for data output within the new urban 
!--       surface model (urban_surface_mod.f90), see also SELECT CASE ( trimvar )
          trimvar = TRIM( doav(ii) )
!
!--    Store the array chosen on the temporary array.
       SELECT CASE ( trimvar )
          CASE ( 'e' )
             IF ( ALLOCATED( e_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         e_av(k,j,i) = e_av(k,j,i) + e(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'lpt' )
             IF ( ALLOCATED( lpt_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         lpt_av(k,j,i) = lpt_av(k,j,i) + pt(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'lwp*' )
             IF ( ALLOCATED( lwp_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      lwp_av(j,i) = lwp_av(j,i) + SUM( ql(nzb:nzt,j,i)            &
                                                  * dzw(1:nzt+1) ) * rho_surface
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'nc' )
             IF ( ALLOCATED( nc_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         nc_av(k,j,i) = nc_av(k,j,i) + nc(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'nr' )
             IF ( ALLOCATED( nr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         nr_av(k,j,i) = nr_av(k,j,i) + nr(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ol*' )
             IF ( ALLOCATED( ol_av ) ) THEN
                DO  m = 1, surf_def_h(0)%ns
                   i = surf_def_h(0)%i(m)
                   j = surf_def_h(0)%j(m)
                   ol_av(j,i) = ol_av(j,i) + surf_def_h(0)%ol(m)
                ENDDO
             ENDIF

          CASE ( 'p' )
             IF ( ALLOCATED( p_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         p_av(k,j,i) = p_av(k,j,i) + p(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          CASE ( 'prr' )
             IF ( ALLOCATED( prr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         prr_av(k,j,i) = prr_av(k,j,i) + prr(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'pt' )
             IF ( ALLOCATED( pt_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                            pt_av(k,j,i) = pt_av(k,j,i) + pt(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
             ENDIF

          CASE ( 'q' )
             IF ( ALLOCATED( q_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         q_av(k,j,i) = q_av(k,j,i) + q(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'qc' )
             IF ( ALLOCATED( qc_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         qc_av(k,j,i) = qc_av(k,j,i) + qc(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ql' )
             IF ( ALLOCATED( ql_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         ql_av(k,j,i) = ql_av(k,j,i) + ql(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ql_c' )
             IF ( ALLOCATED( ql_c_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         ql_c_av(k,j,i) = ql_c_av(k,j,i) + ql_c(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ql_v' )
             IF ( ALLOCATED( ql_v_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         ql_v_av(k,j,i) = ql_v_av(k,j,i) + ql_v(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'qr' )
             IF ( ALLOCATED( qr_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         qr_av(k,j,i) = qr_av(k,j,i) + qr(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'qsws*' )
!
!--          In case of default surfaces, clean-up flux by density.
!--          In case of land- and urban-surfaces, convert fluxes into
!--          dynamic units.
             IF ( ALLOCATED( qsws_av ) ) THEN 
                DO  m = 1, surf_def_h(0)%ns
                   i = surf_def_h(0)%i(m)
                   j = surf_def_h(0)%j(m)
                   k = surf_def_h(0)%k(m)
                   qsws_av(j,i) = qsws_av(j,i) + surf_def_h(0)%qsws(m) *          &
                                                 waterflux_output_conversion(k)
                ENDDO
             ENDIF

          CASE ( 'qv' )
             IF ( ALLOCATED( qv_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         qv_av(k,j,i) = qv_av(k,j,i) + q(k,j,i) - ql(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          CASE ( 'solar3d' )
             IF ( ALLOCATED( solar3d_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         solar3d_av(k,j,i) = solar3d_av(k,j,i) + solar3d(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF  


          CASE ( 'rho_ocean' )
             IF ( ALLOCATED( rho_ocean_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rho_ocean_av(k,j,i) = rho_ocean_av(k,j,i) + rho_ocean(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF  

          CASE ( 'alpha_T' )
             IF ( ALLOCATED( alpha_T_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         alpha_T_av(k,j,i) = alpha_T_av(k,j,i) + alpha_T(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF  

          CASE ( 'beta_S' )
             IF ( ALLOCATED( beta_S_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         beta_S_av(k,j,i) = beta_S_av(k,j,i) + beta_S(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF  


          CASE ( 's' )
             IF ( ALLOCATED( s_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         s_av(k,j,i) = s_av(k,j,i) + s(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'sa' )
             IF ( ALLOCATED( sa_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         sa_av(k,j,i) = sa_av(k,j,i) + sa(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'shf*' )
!
!--          In case of default surfaces, clean-up flux by density.
!--          In case of land- and urban-surfaces, convert fluxes into
!--          dynamic units.
             IF ( ALLOCATED( shf_av ) ) THEN 
                DO  m = 1, surf_def_h(0)%ns
                   i = surf_def_h(0)%i(m)
                   j = surf_def_h(0)%j(m)
                   k = surf_def_h(0)%k(m)
                   shf_av(j,i) = shf_av(j,i) + surf_def_h(0)%shf(m)  *            &
                                               heatflux_output_conversion(k)
                ENDDO
             ENDIF

          CASE ( 'shf_sol*' )
             IF ( ALLOCATED( shf_sol_av ) ) THEN 
                DO  m = 1, surf_def_h(0)%ns
                   i = surf_def_h(0)%i(m)
                   j = surf_def_h(0)%j(m)
                   k = surf_def_h(0)%k(m)
                   shf_sol_av(j,i) = shf_sol_av(j,i) + surf_def_h(0)%shf_sol(m)  *            &
                                               heatflux_output_conversion(k)
                ENDDO
             ENDIF

          CASE ( 'ssws*' )
             IF ( ALLOCATED( ssws_av ) ) THEN 
                DO  m = 1, surf_def_h(0)%ns
                   i = surf_def_h(0)%i(m)
                   j = surf_def_h(0)%j(m)
                   ssws_av(j,i) = ssws_av(j,i) + surf_def_h(0)%ssws(m)
                ENDDO
             ENDIF

          CASE ( 't*' )
             IF ( ALLOCATED( ts_av ) ) THEN 
                DO  m = 1, surf_def_h(0)%ns
                   i = surf_def_h(0)%i(m)
                   j = surf_def_h(0)%j(m)
                   ts_av(j,i) = ts_av(j,i) + surf_def_h(0)%ts(m)
                ENDDO
             ENDIF

          CASE ( 'tsurf*' )
             IF ( ALLOCATED( tsurf_av ) ) THEN             
                DO  m = 1, surf_def_h(0)%ns
                   i   = surf_def_h(0)%i(m)            
                   j   = surf_def_h(0)%j(m)
                   tsurf_av(j,i) = tsurf_av(j,i) + surf_def_h(0)%pt_surface(m)
                ENDDO

             ENDIF

          CASE ( 'u' )
             IF ( ALLOCATED( u_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         u_av(k,j,i) = u_av(k,j,i) + u(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'u*' )
             IF ( ALLOCATED( us_av ) ) THEN   
                DO  m = 1, surf_def_h(0)%ns
                   i = surf_def_h(0)%i(m)
                   j = surf_def_h(0)%j(m)
                   us_av(j,i) = us_av(j,i) + surf_def_h(0)%us(m)
                ENDDO
             ENDIF

          CASE ( 'v' )
             IF ( ALLOCATED( v_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         v_av(k,j,i) = v_av(k,j,i) + v(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vpt' )
             IF ( ALLOCATED( vpt_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         vpt_av(k,j,i) = vpt_av(k,j,i) + vpt(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'w' )
             IF ( ALLOCATED( w_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         w_av(k,j,i) = w_av(k,j,i) + w(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'z0*' )
             IF ( ALLOCATED( z0_av ) ) THEN 
                DO  m = 1, surf_def_h(0)%ns
                   i = surf_def_h(0)%i(m)
                   j = surf_def_h(0)%j(m)
                   z0_av(j,i) = z0_av(j,i) + surf_def_h(0)%z0(m)
                ENDDO
             ENDIF

          CASE ( 'z0h*' )
             IF ( ALLOCATED( z0h_av ) ) THEN 
                DO  m = 1, surf_def_h(0)%ns
                   i = surf_def_h(0)%i(m)
                   j = surf_def_h(0)%j(m)
                   z0h_av(j,i) = z0h_av(j,i) + surf_def_h(0)%z0h(m)
                ENDDO
             ENDIF
   
          CASE ( 'z0q*' )
             IF ( ALLOCATED( z0q_av ) ) THEN 
                DO  m = 1, surf_def_h(0)%ns
                   i = surf_def_h(0)%i(m)
                   j = surf_def_h(0)%j(m)
                   z0q_av(j,i) = z0q_av(j,i) + surf_def_h(0)%z0q(m)
                ENDDO
             ENDIF
          CASE DEFAULT
!
!--          Turbulence closure module
             CALL tcm_3d_data_averaging( 'sum', doav(ii) )

       END SELECT

    ENDDO

    CALL cpu_log( log_point(34), 'sum_up_3d_data', 'stop' )


 END SUBROUTINE sum_up_3d_data
