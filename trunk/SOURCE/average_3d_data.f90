!> @file average_3d_data.f90
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
! $Id: average_3d_data.f90 3004 2018-04-27 12:33:25Z Giersch $
! Further allocation checks implemented, case z0q* added
! 
! 2885 2018-03-14 11:02:46Z Giersch
! Preliminary gust module interface implemented
! 
! 2797 2018-02-08 13:24:35Z suehring
! Enable output of ground-heat flux also at urban surfaces.
! 
! 2766 2018-01-22 17:17:47Z kanani
! Removed preprocessor directive __chem
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
! Change in file header (GPL part)
! Implement call for turbulence_closure_mod (TG)
! Implementation of chemistry module (FK)
! 
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison' 
! includes two more prognostic equations for cloud drop concentration (nc)  
! and cloud water content (qc). 
! 
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new surface concept - additional ghost point exchange 
! of surface variable required
! 
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean and rho_av to rho_ocean_av
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
! 1972 2016-07-26 07:52:02Z maronga
! Output of land surface quantities is now done directly in the respective module
! 
! 1960 2016-07-12 16:34:24Z suehring
! Treat humidity and passive scalar separately
! 
! 1691 2015-10-26 16:17:44Z maronga
! Added output of Obukhov length and radiative heating rates for RRTMG. 
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
! Added support for land surface and radiation model parameters.
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL functions provided with KIND-attribute
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
! 1115 2013-03-26 18:16:16Z hoffmann
! +qc
!
! 1053 2012-11-13 17:11:03Z hoffmann
! averaging of nr, qr added
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 978 2012-08-09 08:28:32Z fricke
! +z0h_av
!
! Revision 1.1  2006/02/23 09:48:58  raasch
! Initial revision
!
!
! Description:
! ------------
!> Time-averaging of 3d-data-arrays.
!------------------------------------------------------------------------------!
 SUBROUTINE average_3d_data
 

    USE averaging

    USE control_parameters,                                                    &
        ONLY:  air_chemistry, average_count_3d, doav, doav_n, land_surface,    &
               urban_surface, varnamelength

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point

    USE indices,                                                               &
        ONLY:  nbgp, nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt

    USE kinds

    USE turbulence_closure_mod,                                                &
        ONLY:  tcm_3d_data_averaging



    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< running index
    INTEGER(iwp) ::  ii !< running index
    INTEGER(iwp) ::  j  !< running index
    INTEGER(iwp) ::  k  !< running index

    CHARACTER (LEN=varnamelength) ::  trimvar  !< TRIM of output-variable string


    CALL cpu_log (log_point(35),'average_3d_data','start')

!
!-- Check, if averaging is necessary
    IF ( average_count_3d <= 1 )  RETURN

!
!-- Loop of all variables to be averaged.
    DO  ii = 1, doav_n

!
!--    Temporary solution to account for data output within the new urban 
!--    surface model (urban_surface_mod.f90), see also SELECT CASE ( trimvar )
       trimvar = TRIM( doav(ii) )
       IF ( urban_surface  .AND.  trimvar(1:4) == 'usm_' )  THEN
          trimvar = 'usm_output'
       ENDIF

!
!--    Store the array chosen on the temporary array.
       SELECT CASE ( trimvar )

          CASE ( 'e' )
             IF ( ALLOCATED( e_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         e_av(k,j,i) = e_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ghf*' )
             IF ( ALLOCATED( ghf_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      ghf_av(j,i) = ghf_av(j,i)                                   &
                                    / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( ghf_av, nbgp )
             ENDIF

          CASE ( 'qsws*' )
             IF ( ALLOCATED( qsws_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      qsws_av(j,i) = qsws_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( qsws_av, nbgp )
             ENDIF

          CASE ( 'lpt' )
             IF ( ALLOCATED( lpt_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         lpt_av(k,j,i) = lpt_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'lwp*' )
             IF ( ALLOCATED( lwp_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      lwp_av(j,i) = lwp_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'nc' )
             IF ( ALLOCATED( nc_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         nc_av(k,j,i) = nc_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'nr' )
             IF ( ALLOCATED( nr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         nr_av(k,j,i) = nr_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

         CASE ( 'ol*' )
             IF ( ALLOCATED( ol_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      ol_av(j,i) = ol_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( ol_av, nbgp )
             ENDIF

          CASE ( 'p' )
             IF ( ALLOCATED( p_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         p_av(k,j,i) = p_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'pc' )
             IF ( ALLOCATED( pc_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         pc_av(k,j,i) = pc_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'pr' )
             IF ( ALLOCATED( pr_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         pr_av(k,j,i) = pr_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'prr' )
             IF ( ALLOCATED( prr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         prr_av(k,j,i) = prr_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'pt' )
             IF ( ALLOCATED( pt_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         pt_av(k,j,i) = pt_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'q' )
             IF ( ALLOCATED( q_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         q_av(k,j,i) = q_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'qc' )
             IF ( ALLOCATED( qc_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         qc_av(k,j,i) = qc_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ql' )
             IF ( ALLOCATED( ql_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         ql_av(k,j,i) = ql_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ql_c' )
             IF ( ALLOCATED( ql_c_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         ql_c_av(k,j,i) = ql_c_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ql_v' )
             IF ( ALLOCATED( ql_v_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         ql_v_av(k,j,i) = ql_v_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ql_vp' )
             IF ( ALLOCATED( ql_vp_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         ql_vp_av(k,j,i) = ql_vp_av(k,j,i) /                      &
                                           REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'qr' )
             IF ( ALLOCATED( qr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         qr_av(k,j,i) = qr_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'qv' )
             IF ( ALLOCATED( qv_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         qv_av(k,j,i) = qv_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

         CASE ( 'r_a*' )
             IF ( ALLOCATED( r_a_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      r_a_av(j,i) = r_a_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( r_a_av, nbgp )
             ENDIF

          CASE ( 'rho_ocean' )
             IF ( ALLOCATED( rho_ocean_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rho_ocean_av(k,j,i) = rho_ocean_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'beta_S' )
             IF ( ALLOCATED( beta_S_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         beta_S_av(k,j,i) = beta_S_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'alpha_T' )
             IF ( ALLOCATED( alpha_T_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         alpha_T_av(k,j,i) = alpha_T_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'solar3d' )
             IF ( ALLOCATED( solar3d_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         solar3d_av(k,j,i) = solar3d_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
          ENDIF


          CASE ( 's' )
             IF ( ALLOCATED( s_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         s_av(k,j,i) = s_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'sa' )
             IF ( ALLOCATED( sa_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         sa_av(k,j,i) = sa_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

         CASE ( 'shf*' )
             IF ( ALLOCATED( shf_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      shf_av(j,i) = shf_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( shf_av, nbgp )
             ENDIF

          CASE ( 'ssws*' )
             IF ( ALLOCATED( ssws_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      ssws_av(j,i) = ssws_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( ssws_av, nbgp )
             ENDIF

          CASE ( 't*' )
             IF ( ALLOCATED( ts_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      ts_av(j,i) = ts_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( ts_av, nbgp )
             ENDIF

         CASE ( 'tsurf*' )
             IF ( ALLOCATED( tsurf_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      tsurf_av(j,i) = tsurf_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( tsurf_av, nbgp )
             ENDIF

          CASE ( 'u' )
             IF ( ALLOCATED( u_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         u_av(k,j,i) = u_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'u*' )
             IF ( ALLOCATED( us_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      us_av(j,i) = us_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( us_av, nbgp )
             ENDIF

          CASE ( 'v' )
             IF ( ALLOCATED( v_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         v_av(k,j,i) = v_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vpt' )
             IF ( ALLOCATED( vpt_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         vpt_av(k,j,i) = vpt_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'w' )
             IF ( ALLOCATED( w_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         w_av(k,j,i) = w_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'z0*' )
             IF ( ALLOCATED( z0_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      z0_av(j,i) = z0_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( z0_av, nbgp )
             ENDIF

          CASE ( 'z0h*' )
             IF ( ALLOCATED( z0h_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      z0h_av(j,i) = z0h_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( z0h_av, nbgp )
             ENDIF

          CASE ( 'z0q*' )
             IF ( ALLOCATED( z0q_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      z0q_av(j,i) = z0q_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( z0q_av, nbgp )
             ENDIF
          CASE DEFAULT
!
!--          Turbulence closure module
             CALL tcm_3d_data_averaging( 'average', doav(ii) )

       END SELECT

    ENDDO

!
!-- Reset the counter
    average_count_3d = 0.0

    CALL cpu_log( log_point(35), 'average_3d_data', 'stop' )


 END SUBROUTINE average_3d_data
