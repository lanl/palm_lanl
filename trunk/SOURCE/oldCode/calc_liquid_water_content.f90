!> @file calc_liquid_water_content.f90
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
! $Id: calc_liquid_water_content.f90 2719 2018-01-02 09:02:06Z maronga $
! Bugfix for last change.
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2608 2017-11-13 14:04:26Z schwenkel
! Calculation of supersaturation in external module (diagnostic_quantities_mod).
! Change: correct calculation of saturation specific humidity to saturation
! mixing ratio (the factor of 0.378 vanishes).
! 
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison' 
! includes two more prognostic equations for cloud drop concentration (nc)  
! and cloud water content (qc). 
! 
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography concept
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1822 2016-04-07 07:49:42Z hoffmann
! icloud_scheme removed. microphysics_seifert added.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1253 2013-11-07 10:48:12Z fricke
! Bugfix: q is set to qr in case that q is smaller than qr
!
! 1115 2013-03-26 18:16:16Z hoffmann
! drizzle can be used independently from precipitation
!
! 1053 2012-11-13 17:11:03Z hoffmann
! description expanded to the two-moment cloud scheme
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! Revision 1.1  2000/04/13 14:50:45  schroeter
! Initial revision
!
!
!
! Description:
! ------------
!> Calculation of the liquid water content (0%-or-100%-scheme). This scheme is 
!> used by the one and the two moment cloud physics scheme. Using the two moment 
!> scheme, this calculation results in the cloud water content.
!------------------------------------------------------------------------------!
 SUBROUTINE calc_liquid_water_content
 


    USE arrays_3d,                                                             &
        ONLY:  q, qc, ql, qr

    USE control_parameters,                                                    &
        ONLY:  microphysics_morrison, microphysics_seifert

    USE diagnostic_quantities_mod,                                             &
        ONLY:  e_s, magnus, q_s, sat, supersaturation, t_l

    USE indices,                                                               &
        ONLY:  nxlg, nxrg, nyng, nysg, nzb, nzt, wall_flags_0

    USE kinds

    USE pegrid


    IMPLICIT NONE

    INTEGER(iwp) ::  i !<
    INTEGER(iwp) ::  j !<
    INTEGER(iwp) ::  k !<


    DO  i = nxlg, nxrg
       DO  j = nysg, nyng
          DO  k = nzb+1, nzt

!
!--          Call calculation of supersaturation located 
!--          in diagnostic_quantities_mod
             CALL supersaturation ( i, j, k )

!
!--          Compute the liquid water content
             IF ( microphysics_seifert  .AND.  .NOT. microphysics_morrison )   &
             THEN
                IF ( ( q(k,j,i) - q_s - qr(k,j,i) ) > 0.0_wp )  THEN
                   qc(k,j,i) = ( q(k,j,i) - q_s - qr(k,j,i) )                  &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 0 ) )
                   ql(k,j,i) = ( qc(k,j,i) + qr(k,j,i) )                       &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 0 ) )
                ELSE
                   IF ( q(k,j,i) < qr(k,j,i) )  q(k,j,i) = qr(k,j,i)
                   qc(k,j,i) = 0.0_wp 
                   ql(k,j,i) = qr(k,j,i)                                       &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 0 ) )
                ENDIF
             ELSEIF ( microphysics_morrison )  THEN
                ql(k,j,i) = qc(k,j,i) + qr(k,j,i)                              &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 0 ) )
             ELSE
                IF ( ( q(k,j,i) - q_s ) > 0.0_wp )  THEN
                   qc(k,j,i) = ( q(k,j,i) - q_s )                              &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 0 ) )
                   ql(k,j,i) = qc(k,j,i)                                       &
                                      * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_0(k,j,i), 0 ) )
                ELSE
                   qc(k,j,i) = 0.0_wp
                   ql(k,j,i) = 0.0_wp
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    
 END SUBROUTINE calc_liquid_water_content
