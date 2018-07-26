!> @file inflow_turbulence.f90
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
! $Id: inflow_turbulence.f90 2718 2018-01-02 08:49:38Z maronga $
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
! 1960 2016-07-12 16:34:24Z suehring
! Separate humidity and passive scalar
! 
! 1806 2016-04-05 18:55:35Z gronemeier
! Added comments to variables and code segments. Removed code redundancies.
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
! 
! 1615 2015-07-08 18:49:19Z suehring
! Enable turbulent inflow for passive_scalar and humidity
!
! 1560 2015-03-06 10:48:54Z keck
! Option recycling_yshift added. If this option is switched on, the turbulence
! data, which is mapped from the recycling plane to the inflow, is shifted in
! y direction (by ny * dy / 2 )
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute 
!
! 1346 2014-03-27 13:18:20Z heinze
! Bugfix: REAL constants provided with KIND-attribute especially in call of 
! intrinsic function like MAX, MIN, SIGN
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
! Initial version (2008/03/07)
!
! Description:
! ------------
!> Imposing turbulence at the respective inflow using the turbulence
!> recycling method of Kataoka and Mizuno (2002).
!------------------------------------------------------------------------------!
 SUBROUTINE inflow_turbulence
 

    USE arrays_3d,                                                             &
        ONLY:  e, inflow_damping_factor, mean_inflow_profiles, pt, q, s, u, v, w
        
    USE control_parameters,                                                    &
        ONLY:  humidity, passive_scalar, recycling_plane, recycling_yshift
        
    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point
        
    USE indices,                                                               &
        ONLY:  nbgp, nxl, ny, nyn, nys, nyng, nysg, nzb, nzt
        
    USE kinds
    
    USE pegrid


    IMPLICIT NONE

    INTEGER(iwp) ::  i        !< loop index
    INTEGER(iwp) ::  j        !< loop index
    INTEGER(iwp) ::  k        !< loop index
    INTEGER(iwp) ::  l        !< loop index
    INTEGER(iwp) ::  next     !< ID of receiving PE for y-shift
    INTEGER(iwp) ::  ngp_ifd  !< number of grid points stored in avpr
    INTEGER(iwp) ::  ngp_pr   !< number of grid points stored in inflow_dist
    INTEGER(iwp) ::  prev     !< ID of sending PE for y-shift

    REAL(wp), DIMENSION(nzb:nzt+1,7,nbgp)           ::                         &
       avpr               !< stores averaged profiles at recycling plane
    REAL(wp), DIMENSION(nzb:nzt+1,7,nbgp)           ::                         &
       avpr_l             !< auxiliary variable to calculate avpr
    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,7,nbgp) ::                         &
       inflow_dist        !< turbulence signal of vars, added at inflow boundary
    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,7,nbgp) ::                         &
       local_inflow_dist  !< auxiliary variable for inflow_dist, used for yshift

    CALL cpu_log( log_point(40), 'inflow_turbulence', 'start' )

!
!-- Carry out spanwise averaging in the recycling plane
    avpr_l = 0.0_wp
    ngp_pr = ( nzt - nzb + 2 ) * 7 * nbgp
    ngp_ifd = ngp_pr * ( nyn - nys + 1 + 2 * nbgp )

!
!-- First, local averaging within the recycling domain
    i = recycling_plane

#if defined( __parallel )
    IF ( myidx == id_recycling )  THEN
       
       DO  l = 1, nbgp
          DO  j = nys, nyn
             DO  k = nzb, nzt + 1

                avpr_l(k,1,l) = avpr_l(k,1,l) + u(k,j,i)
                avpr_l(k,2,l) = avpr_l(k,2,l) + v(k,j,i)
                avpr_l(k,3,l) = avpr_l(k,3,l) + w(k,j,i)
                avpr_l(k,4,l) = avpr_l(k,4,l) + pt(k,j,i)
                avpr_l(k,5,l) = avpr_l(k,5,l) + e(k,j,i)
                IF ( humidity )                                                &
                   avpr_l(k,6,l) = avpr_l(k,6,l) + q(k,j,i)
                IF ( passive_scalar )                                          &
                   avpr_l(k,7,l) = avpr_l(k,7,l) + s(k,j,i)

             ENDDO
          ENDDO
          i = i + 1
       ENDDO

    ENDIF
!
!-- Now, averaging over all PEs
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( avpr_l(nzb,1,1), avpr(nzb,1,1), ngp_pr, MPI_REAL,      &
                        MPI_SUM, comm2d, ierr )

#else
    DO  l = 1, nbgp
       DO  j = nys, nyn
          DO  k = nzb, nzt + 1

             avpr_l(k,1,l) = avpr_l(k,1,l) + u(k,j,i)
             avpr_l(k,2,l) = avpr_l(k,2,l) + v(k,j,i)
             avpr_l(k,3,l) = avpr_l(k,3,l) + w(k,j,i)
             avpr_l(k,4,l) = avpr_l(k,4,l) + pt(k,j,i)
             avpr_l(k,5,l) = avpr_l(k,5,l) + e(k,j,i)
             IF ( humidity )                                                   &
                avpr_l(k,6,l) = avpr_l(k,6,l) + q(k,j,i)
             IF ( passive_scalar )                                             &
                avpr_l(k,7,l) = avpr_l(k,7,l) + s(k,j,i)

          ENDDO
       ENDDO
       i = i + 1 
    ENDDO
    
    avpr = avpr_l
#endif

    avpr = avpr / ( ny + 1 )
!
!-- Calculate the disturbances at the recycling plane
    i = recycling_plane

#if defined( __parallel )
    IF ( myidx == id_recycling )  THEN
       DO  l = 1, nbgp
          DO  j = nysg, nyng
             DO  k = nzb, nzt + 1

                inflow_dist(k,j,1,l) = u(k,j,i+1) - avpr(k,1,l)
                inflow_dist(k,j,2,l) = v(k,j,i)   - avpr(k,2,l)
                inflow_dist(k,j,3,l) = w(k,j,i)   - avpr(k,3,l)
                inflow_dist(k,j,4,l) = pt(k,j,i)  - avpr(k,4,l)
                inflow_dist(k,j,5,l) = e(k,j,i)   - avpr(k,5,l)
                IF ( humidity )                                                &
                   inflow_dist(k,j,6,l) = q(k,j,i) - avpr(k,6,l)
                IF ( passive_scalar )                                          &
                   inflow_dist(k,j,7,l) = s(k,j,i) - avpr(k,7,l)
            ENDDO
          ENDDO
          i = i + 1
       ENDDO

    ENDIF
#else
    DO  l = 1, nbgp
       DO  j = nysg, nyng
          DO  k = nzb, nzt+1

             inflow_dist(k,j,1,l) = u(k,j,i+1) - avpr(k,1,l)
             inflow_dist(k,j,2,l) = v(k,j,i)   - avpr(k,2,l)
             inflow_dist(k,j,3,l) = w(k,j,i)   - avpr(k,3,l)
             inflow_dist(k,j,4,l) = pt(k,j,i)  - avpr(k,4,l)
             inflow_dist(k,j,5,l) = e(k,j,i)   - avpr(k,5,l)
             IF ( humidity )                                                   &
                inflow_dist(k,j,6,l) = q(k,j,i) - avpr(k,6,l)
             IF ( passive_scalar )                                             &
                inflow_dist(k,j,7,l) = s(k,j,i) - avpr(k,7,l)
              
          ENDDO
       ENDDO
       i = i + 1
    ENDDO
#endif

!
!-- For parallel runs, send the disturbances to the respective inflow PE
#if defined( __parallel )
    IF ( myidx == id_recycling  .AND.  myidx /= id_inflow )  THEN

       CALL MPI_SEND( inflow_dist(nzb,nysg,1,1), ngp_ifd, MPI_REAL,            &
                      id_inflow, 1, comm1dx, ierr )

    ELSEIF ( myidx /= id_recycling  .AND.  myidx == id_inflow )  THEN

       inflow_dist = 0.0_wp
       CALL MPI_RECV( inflow_dist(nzb,nysg,1,1), ngp_ifd, MPI_REAL,            &
                      id_recycling, 1, comm1dx, status, ierr )

    ENDIF

!
!-- y-shift for inflow_dist
!-- Shift inflow_dist in positive y direction by a distance of INT( npey / 2 )
    IF ( recycling_yshift .AND. myidx == id_inflow ) THEN
!
!--    Calculate the ID of the PE which sends data to this PE (prev) and of the
!--    PE which receives data from this PE (next).
       IF ( myidy >= INT( pdims(2) / 2 ) ) THEN
          prev = myidy - INT( pdims(2) / 2 )
       ELSE
          prev = pdims(2) - ( INT( pdims(2) / 2 ) - myidy )
       ENDIF
     
       IF ( myidy < pdims(2) - INT( pdims(2) / 2 ) ) THEN
          next = myidy + INT( pdims(2) / 2 )
       ELSE
          next = INT( pdims(2) / 2 ) - ( pdims(2) - myidy )
       ENDIF

       local_inflow_dist = 0.0_wp

       CALL MPI_SENDRECV( inflow_dist(nzb,nysg,1,1), ngp_ifd, MPI_REAL,        &
                          next, 1, local_inflow_dist(nzb,nysg,1,1), ngp_ifd,   &
                          MPI_REAL, prev, 1, comm1dy, status, ierr )

       inflow_dist = local_inflow_dist

    ENDIF

#endif

!
!-- Add the disturbance at the inflow
    IF ( nxl == 0 )  THEN

       DO  j = nysg, nyng
          DO  k = nzb, nzt + 1

             u(k,j,-nbgp+1:0) = mean_inflow_profiles(k,1) +                 &
                        inflow_dist(k,j,1,1:nbgp) * inflow_damping_factor(k)
             v(k,j,-nbgp:-1)  = mean_inflow_profiles(k,2) +                 &
                        inflow_dist(k,j,2,1:nbgp) * inflow_damping_factor(k)
             w(k,j,-nbgp:-1)  =                                             &
                        inflow_dist(k,j,3,1:nbgp) * inflow_damping_factor(k)
             pt(k,j,-nbgp:-1) = mean_inflow_profiles(k,4) +                 &
                        inflow_dist(k,j,4,1:nbgp) * inflow_damping_factor(k)
             e(k,j,-nbgp:-1)  = mean_inflow_profiles(k,5) +                 &
                        inflow_dist(k,j,5,1:nbgp) * inflow_damping_factor(k)
             e(k,j,-nbgp:-1)  = MAX( e(k,j,-nbgp:-1), 0.0_wp )

             IF ( humidity )                                                &
                q(k,j,-nbgp:-1)  = mean_inflow_profiles(k,6) +              &
                        inflow_dist(k,j,6,1:nbgp) * inflow_damping_factor(k)
             IF ( passive_scalar )                                          &
                s(k,j,-nbgp:-1)  = mean_inflow_profiles(k,7) +              &
                        inflow_dist(k,j,7,1:nbgp) * inflow_damping_factor(k)

          ENDDO
       ENDDO

    ENDIF


    CALL cpu_log( log_point(40), 'inflow_turbulence', 'stop' )


 END SUBROUTINE inflow_turbulence
