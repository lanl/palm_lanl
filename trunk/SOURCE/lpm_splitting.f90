!> @file lpm_splitting.f90
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
! $Id: lpm_splitting.f90 2932 2018-03-26 09:39:22Z maronga $
! renamed particles_par to particle_parameters
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 
! Change in file header (GPL part)
!
! Added comments
! 
! 
! 2263 2017-06-08 14:59:01Z schwenkel
! Initial revision
! 
! 
!
! Description:
! ------------
! This routine is a part of the Lagrangian particle model. Super droplets which 
! fulfill certain criterion's (e.g. a big weighting factor and a large radius) 
! can be split into several super droplets with a reduced number of 
! represented particles of every super droplet. This mechanism ensures an
! improved representation of the right tail of the drop size distribution with 
! a feasible amount of computational costs. The limits of particle creation 
! should be chosen carefully! The idea of this algorithm is based on 
! Unterstrasser and Soelch, 2014. 
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_splitting


    USE arrays_3d,                                                             &
        ONLY:  ql

    USE cloud_parameters,                                                      &
        ONLY:  rho_l

    USE constants,                                                             &
        ONLY:  pi

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE indices,                                                               &
        ONLY:  nxl, nxr, nyn, nys, nzb, nzt

    USE kinds

    USE lpm_exchange_horiz_mod,                                                &
        ONLY:  realloc_particles_array

    USE particle_attributes,                                                   &
        ONLY:  grid_particles, iran_part, initial_weighting_factor, isf,       &
               i_splitting_mode, max_number_particles_per_gridbox,             &  
               new_particles, n_max, number_concentration,                     &
               number_of_particles, number_particles_per_gridbox, particles,   &
               particle_type, prt_count, radius_split, splitting,              &
               splitting_factor, splitting_factor_max, splitting_mode,         &
               sum_new_particles, weight_factor_split                        

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i                !< 
    INTEGER(iwp) ::  j                !<
    INTEGER(iwp) ::  jpp              !<
    INTEGER(iwp) ::  k                !<
    INTEGER(iwp) ::  n                !<
    INTEGER(iwp) ::  new_particles_gb !< counter of created particles within one grid box
    INTEGER(iwp) ::  new_size         !< new particle array size
    INTEGER(iwp) ::  np               !< 
    INTEGER(iwp) ::  old_size         !< old particle array size
    
    LOGICAL ::  first_loop_stride = .TRUE. !< flag to calculate constants only once

    REAL(wp) ::  diameter                 !< diameter of droplet
    REAL(wp) ::  dlog                     !< factor for DSD calculation
    REAL(wp) ::  factor_volume_to_mass    !< pre calculate factor volume to mass
    REAL(wp) ::  lambda                   !< slope parameter of gamma-distribution
    REAL(wp) ::  lwc                      !< liquid water content of grid box
    REAL(wp) ::  lwc_total                !< average liquid water content of cloud
    REAL(wp) ::  m1                       !< first moment of DSD
    REAL(wp) ::  m1_total                 !< average over all PEs of first moment of DSD
    REAL(wp) ::  m2                       !< second moment of DSD
    REAL(wp) ::  m2_total                 !< average average over all PEs second moment of DSD
    REAL(wp) ::  m3                       !< third moment of DSD
    REAL(wp) ::  m3_total                 !< average average over all PEs third moment of DSD
    REAL(wp) ::  mu                       !< spectral shape parameter of gamma distribution
    REAL(wp) ::  nrclgb                   !< number of cloudy grid boxes (ql >= 1.0E-5 kg/kg) 
    REAL(wp) ::  nrclgb_total             !< average over all PEs of number of cloudy grid boxes
    REAL(wp) ::  nr                       !< number concentration of cloud droplets
    REAL(wp) ::  nr_total                 !< average over all PEs of number of cloudy grid boxes
    REAL(wp) ::  nr0                      !< intercept parameter of gamma distribution
    REAL(wp) ::  pirho_l                  !< pi * rho_l / 6.0
    REAL(wp) ::  ql_crit = 1.0E-5_wp      !< threshold lwc for cloudy grid cells 
                                          !< (Siebesma et al 2003, JAS, 60)
    REAL(wp) ::  rm                       !< volume averaged mean radius
    REAL(wp) ::  rm_total                 !< average over all PEs of volume averaged mean radius
    REAL(wp) ::  r_min = 1.0E-6_wp        !< minimum radius of approximated spectra 
    REAL(wp) ::  r_max = 1.0E-3_wp        !< maximum radius of approximated spectra
    REAL(wp) ::  sigma_log = 1.5_wp       !< standard deviation of the LOG-distribution
    REAL(wp) ::  zeta                     !< Parameter for DSD calculation of Seifert

    REAL(wp), DIMENSION(0:n_max-1) ::  an_spl     !< size dependent critical weight factor
    REAL(wp), DIMENSION(0:n_max-1) ::  r_bin_mid  !< mass weighted mean radius of a bin
    REAL(wp), DIMENSION(0:n_max)   ::  r_bin      !< boundaries of a radius bin
    
    TYPE(particle_type) ::  tmp_particle   !< temporary particle TYPE

    CALL cpu_log( log_point_s(80), 'lpm_splitting', 'start' )

    IF ( first_loop_stride )  THEN
       IF ( i_splitting_mode == 2  .OR.  i_splitting_mode == 3 )  THEN
          dlog   = ( LOG10(r_max) - LOG10(r_min) ) / ( n_max - 1 )
          DO  i = 0, n_max-1
             r_bin(i) = 10.0_wp**( LOG10(r_min) + i * dlog - 0.5_wp * dlog )
             r_bin_mid(i) = 10.0_wp**( LOG10(r_min) + i * dlog )
          ENDDO
          r_bin(n_max) = 10.0_wp**( LOG10(r_min) + n_max * dlog - 0.5_wp * dlog )
       ENDIF   
       factor_volume_to_mass =  4.0_wp / 3.0_wp * pi * rho_l
       pirho_l  = pi * rho_l / 6.0_wp
       IF ( weight_factor_split == -1.0_wp )  THEN
          weight_factor_split = 0.1_wp * initial_weighting_factor 
       ENDIF
    ENDIF

    new_particles  = 0

    IF ( i_splitting_mode == 1 )  THEN

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt

                new_particles_gb = 0
                number_of_particles = prt_count(k,j,i)
                IF ( number_of_particles <= 0  .OR.                            &  
                     ql(k,j,i) < ql_crit )  CYCLE
                particles => grid_particles(k,j,i)%particles(1:number_of_particles)
!                
!--             Start splitting operations. Each particle is checked if it
!--             fulfilled the splitting criterion's. In splitting mode 'const'   
!--             a critical radius  (radius_split) a critical weighting factor
!--             (weight_factor_split) and a splitting factor (splitting_factor)
!--             must  be prescribed (see particle_parameters). Super droplets 
!--             which have a larger radius and larger weighting factor are split 
!--             into 'splitting_factor' super droplets. Therefore, the weighting 
!--             factor of  the super droplet and all created clones is reduced 
!--             by the factor of 'splitting_factor'.
                DO  n = 1, number_of_particles
                   IF ( particles(n)%particle_mask  .AND.                      &
                        particles(n)%radius >= radius_split  .AND.             & 
                        particles(n)%weight_factor >= weight_factor_split )    &
                   THEN          
!
!--                   Calculate the new number of particles.
                      new_size = prt_count(k,j,i) + splitting_factor - 1                                    
!
!--                   Cycle if maximum number of particles per grid box
!--                   is greater than the allowed maximum number.
                      IF ( new_size >= max_number_particles_per_gridbox )  CYCLE                      
!
!--                   Reallocate particle array if necessary. 
                      IF ( new_size > SIZE(particles) )  THEN 
                         CALL realloc_particles_array(i,j,k,new_size)
                      ENDIF
                      old_size = prt_count(k,j,i)
!
!--                   Calculate new weighting factor.
                      particles(n)%weight_factor =  & 
                         particles(n)%weight_factor / splitting_factor
                      tmp_particle = particles(n)
!
!--                   Create splitting_factor-1 new particles.
                      DO  jpp = 1, splitting_factor-1
                         grid_particles(k,j,i)%particles(jpp+old_size) =       & 
                            tmp_particle         
                      ENDDO  
                      new_particles_gb = new_particles_gb + splitting_factor - 1
!   
!--                   Save the new number of super droplets for every grid box.
                      prt_count(k,j,i) = prt_count(k,j,i) +                    &
                                         splitting_factor - 1         
                   ENDIF
                ENDDO
                
                new_particles       = new_particles     + new_particles_gb
                sum_new_particles   = sum_new_particles + new_particles_gb 
             ENDDO
          ENDDO
       ENDDO

    ELSEIF ( i_splitting_mode == 2 )  THEN 
!
!--    Initialize summing variables.
       lwc          = 0.0_wp
       lwc_total    = 0.0_wp 
       m1           = 0.0_wp
       m1_total     = 0.0_wp
       m2           = 0.0_wp
       m2_total     = 0.0_wp
       m3           = 0.0_wp
       m3_total     = 0.0_wp
       nr           = 0.0_wp   
       nrclgb       = 0.0_wp
       nrclgb_total = 0.0_wp 
       nr_total     = 0.0_wp
       rm           = 0.0_wp
       rm_total     = 0.0_wp
       
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                number_of_particles = prt_count(k,j,i)
                IF ( number_of_particles <= 0  .OR.                            & 
                     ql(k,j,i) < ql_crit )  CYCLE
                particles => grid_particles(k,j,i)%particles(1:number_of_particles)
                nrclgb = nrclgb + 1.0_wp
!                
!--             Calculate moments of DSD.                
                DO  n = 1, number_of_particles
                   IF ( particles(n)%particle_mask  .AND.                      &
                        particles(n)%radius >= r_min )                         &
                   THEN
                      nr  = nr  + particles(n)%weight_factor
                      rm  = rm  + factor_volume_to_mass  *                     &
                                 particles(n)%radius**3  *                     &
                                 particles(n)%weight_factor
                      IF ( isf == 1 )  THEN           
                         diameter   = particles(n)%radius * 2.0_wp            
                         lwc = lwc + factor_volume_to_mass *                   &
                                     particles(n)%radius**3 *                  & 
                                     particles(n)%weight_factor 
                         m1  = m1  + particles(n)%weight_factor * diameter                                                
                         m2  = m2  + particles(n)%weight_factor * diameter**2           
                         m3  = m3  + particles(n)%weight_factor * diameter**3
                      ENDIF   
                   ENDIF                        
                ENDDO 
             ENDDO
          ENDDO
       ENDDO

#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( nr, nr_total, 1 , &
       MPI_REAL, MPI_SUM, comm2d, ierr )
       CALL MPI_ALLREDUCE( rm, rm_total, 1 , &
       MPI_REAL, MPI_SUM, comm2d, ierr )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( nrclgb, nrclgb_total, 1 , &
       MPI_REAL, MPI_SUM, comm2d, ierr )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( lwc, lwc_total, 1 , &
       MPI_REAL, MPI_SUM, comm2d, ierr )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( m1, m1_total, 1 , &
       MPI_REAL, MPI_SUM, comm2d, ierr )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( m2, m2_total, 1 , &
       MPI_REAL, MPI_SUM, comm2d, ierr )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( m3, m3_total, 1 , &
       MPI_REAL, MPI_SUM, comm2d, ierr )
#endif 

!
!--    Calculate number concentration and mean volume averaged radius.
       nr_total = MERGE( nr_total / nrclgb_total,                              &
                         0.0_wp, nrclgb_total > 0.0_wp                         &
                       )
       rm_total = MERGE( ( rm_total /                                          &
                            ( nr_total * factor_volume_to_mass )               &
                          )**0.3333333_wp, 0.0_wp, nrclgb_total > 0.0_wp       &
                       )                         
!
!--    Check which function should be used to approximate the DSD.
       IF ( isf == 1 )  THEN
          lwc_total = MERGE( lwc_total / nrclgb_total,                         &
                             0.0_wp, nrclgb_total > 0.0_wp                     &
                           )
          m1_total  = MERGE( m1_total / nrclgb_total,                          &
                             0.0_wp, nrclgb_total > 0.0_wp                     &
                           )
          m2_total  = MERGE( m2_total / nrclgb_total,                          &
                             0.0_wp, nrclgb_total > 0.0_wp                     &
                           )
          m3_total  = MERGE( m3_total / nrclgb_total,                          &
                             0.0_wp, nrclgb_total > 0.0_wp                     &
                           )
          zeta = m1_total * m3_total / m2_total**2                             
          mu   = MAX( ( ( 1.0_wp - zeta ) * 2.0_wp + 1.0_wp ) /                &
                        ( zeta - 1.0_wp ), 0.0_wp                              &
                    )

          lambda = ( pirho_l * nr_total / lwc_total *                          &
                     ( mu + 3.0_wp ) * ( mu + 2.0_wp ) * ( mu + 1.0_wp )       &                                          
                   )**0.3333333_wp
          nr0 = nr_total / gamma( mu + 1.0_wp ) * lambda**( mu + 1.0_wp ) 
          
          DO  n = 0, n_max-1
             diameter  = r_bin_mid(n) * 2.0_wp           
             an_spl(n) = nr0 * diameter**mu * EXP( -lambda * diameter ) *      & 
                         ( r_bin(n+1) - r_bin(n) ) * 2.0_wp 
          ENDDO
       ELSEIF ( isf == 2 )  THEN
          DO  n = 0, n_max-1
             an_spl(n) = nr_total / ( SQRT( 2.0_wp * pi ) *                    &
                                     LOG(sigma_log) * r_bin_mid(n)             &
                                     ) *                                       &
                         EXP( -( LOG( r_bin_mid(n) / rm_total )**2 ) /         &
                               ( 2.0_wp * LOG(sigma_log)**2 )                  & 
                             ) *                                               & 
                         ( r_bin(n+1) - r_bin(n) )
          ENDDO
       ELSEIF( isf == 3 )  THEN
          DO  n = 0, n_max-1 
             an_spl(n) = 3.0_wp * nr_total * r_bin_mid(n)**2 / rm_total**3  *  &
                         EXP( - ( r_bin_mid(n)**3 / rm_total**3 ) )         *  &
                         ( r_bin(n+1) - r_bin(n) )
          ENDDO
       ENDIF
!
!--    Criterion to avoid super droplets with a weighting factor < 1.0.
       an_spl = MAX(an_spl, 1.0_wp)
  
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                number_of_particles = prt_count(k,j,i)
                IF ( number_of_particles <= 0  .OR.                            &  
                     ql(k,j,i) < ql_crit )  CYCLE
                particles => grid_particles(k,j,i)%particles(1:number_of_particles)
                new_particles_gb = 0         
!                
!--             Start splitting operations. Each particle is checked if it
!--             fulfilled the splitting criterion's. In splitting mode 'cl_av'   
!--             a critical radius (radius_split) and a splitting function must 
!--             be prescribed (see particles_par). The critical weighting factor
!--             is calculated while approximating a 'gamma', 'log' or 'exp'-
!--             drop size distribution. In this mode the DSD is calculated as 
!--             an average over all cloudy grid boxes. Super droplets which 
!--             have a larger radius and larger weighting factor are split into
!--             'splitting_factor' super droplets. In this case the splitting  
!--             factor is calculated of weighting factor of the super droplet  
!--             and the approximated number concentration for droplet of such 
!--             a size. Due to the splitting, the weighting factor of the  
!--             super droplet and all created clones is reduced by the factor  
!--             of 'splitting_facor'.
                DO  n = 1, number_of_particles
                   DO  np = 0, n_max-1
                      IF ( r_bin(np) >= radius_split  .AND.                    &
                           particles(n)%particle_mask  .AND.                   &
                           particles(n)%radius >= r_bin(np)  .AND.             &
                           particles(n)%radius < r_bin(np+1)  .AND.            &
                           particles(n)%weight_factor >= an_spl(np)  )         &
                      THEN
!
!--                      Calculate splitting factor
                         splitting_factor =                                    & 
                             MIN( INT( particles(n)%weight_factor /            &
                                        an_spl(np)                             &
                                     ), splitting_factor_max                   &
                                )
                         IF ( splitting_factor < 2 )  CYCLE
!
!--                      Calculate the new number of particles.                                                           
                         new_size = prt_count(k,j,i) + splitting_factor - 1
!
!--                      Cycle if maximum number of particles per grid box
!--                      is greater than the allowed maximum number.                         
                         IF ( new_size >= max_number_particles_per_gridbox )   & 
                         CYCLE
!
!--                      Reallocate particle array if necessary. 
                         IF ( new_size > SIZE(particles) )  THEN 
                            CALL realloc_particles_array(i,j,k,new_size)
                         ENDIF
                         old_size  = prt_count(k,j,i)                             
                         new_particles_gb = new_particles_gb +                 &
                                            splitting_factor - 1
!
!--                      Calculate new weighting factor.
                         particles(n)%weight_factor =                          & 
                            particles(n)%weight_factor / splitting_factor
                         tmp_particle = particles(n)
!
!--                      Create splitting_factor-1 new particles.                                                  
                         DO  jpp = 1, splitting_factor-1
                            grid_particles(k,j,i)%particles(jpp+old_size) =    &
                                                                    tmp_particle         
                         ENDDO
!   
!--                      Save the new number of super droplets. 
                         prt_count(k,j,i) = prt_count(k,j,i) +                 &
                                            splitting_factor - 1
                      ENDIF
                   ENDDO
                ENDDO 
                
                new_particles       = new_particles     + new_particles_gb
                sum_new_particles   = sum_new_particles + new_particles_gb                    
             ENDDO
          ENDDO
       ENDDO

    ELSEIF ( i_splitting_mode == 3 )  THEN 

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
             
!
!--             Initialize summing variables.             
                lwc = 0.0_wp
                m1  = 0.0_wp
                m2  = 0.0_wp
                m3  = 0.0_wp
                nr  = 0.0_wp
                rm  = 0.0_wp  
                
                new_particles_gb = 0
                number_of_particles = prt_count(k,j,i)
                IF ( number_of_particles <= 0  .OR.                            & 
                     ql(k,j,i) < ql_crit )  CYCLE
                particles => grid_particles(k,j,i)%particles
!                
!--             Calculate moments of DSD.                
                DO  n = 1, number_of_particles
                   IF ( particles(n)%particle_mask  .AND.                      &
                        particles(n)%radius >= r_min )                         &
                   THEN
                      nr  = nr + particles(n)%weight_factor
                      rm  = rm + factor_volume_to_mass  *                      &
                                 particles(n)%radius**3  *                     &
                                 particles(n)%weight_factor
                      IF ( isf == 1 )  THEN            
                         diameter   = particles(n)%radius * 2.0_wp           
                         lwc = lwc + factor_volume_to_mass *                   &
                                     particles(n)%radius**3 *                  &  
                                     particles(n)%weight_factor 
                         m1  = m1 + particles(n)%weight_factor * diameter
                         m2  = m2 + particles(n)%weight_factor * diameter**2
                         m3  = m3 + particles(n)%weight_factor * diameter**3
                      ENDIF      
                   ENDIF                                            
                ENDDO  

                IF ( nr <= 0.0  .OR.  rm <= 0.0_wp )  CYCLE 
!
!--             Calculate mean volume averaged radius.                
                rm = ( rm / ( nr * factor_volume_to_mass ) )**0.3333333_wp
!
!--             Check which function should be used to approximate the DSD.              
                IF ( isf == 1 )  THEN
!
!--                Gamma size distribution to calculate  
!--                critical weight_factor (e.g. Marshall + Palmer, 1948).
                   zeta = m1 * m3 / m2**2
                   mu   = MAX( ( ( 1.0_wp - zeta ) * 2.0_wp + 1.0_wp ) /       &
                                ( zeta - 1.0_wp ), 0.0_wp                      &
                             )   
                   lambda = ( pirho_l * nr / lwc *                             &
                              ( mu + 3.0_wp ) * ( mu + 2.0_wp ) *              &
                              ( mu + 1.0_wp )                                  &                                  
                            )**0.3333333_wp
                   nr0 =  ( nr / (gamma( mu + 1.0_wp ) ) ) *                   &
                          lambda**( mu + 1.0_wp ) 

                   DO  n = 0, n_max-1
                      diameter         = r_bin_mid(n) * 2.0_wp           
                      an_spl(n) = nr0 * diameter**mu *                         &
                                  EXP( -lambda * diameter ) *                  & 
                                  ( r_bin(n+1) - r_bin(n) ) * 2.0_wp 
                   ENDDO
                ELSEIF ( isf == 2 )  THEN
!
!--                Lognormal size distribution to calculate critical 
!--                weight_factor (e.g. Levin, 1971, Bradley + Stow, 1974).   
                   DO  n = 0, n_max-1
                      an_spl(n) = nr / ( SQRT( 2.0_wp * pi ) *                 &
                                              LOG(sigma_log) * r_bin_mid(n)    &
                                        ) *                                    &
                                  EXP( -( LOG( r_bin_mid(n) / rm )**2 ) /      &
                                        ( 2.0_wp * LOG(sigma_log)**2 )         & 
                                      ) *                                      & 
                                  ( r_bin(n+1) - r_bin(n) )
                   ENDDO
                ELSEIF ( isf == 3 )  THEN
!
!--                Exponential size distribution to calculate critical 
!--                weight_factor (e.g. Berry + Reinhardt, 1974).  
                   DO  n = 0, n_max-1
                      an_spl(n) = 3.0_wp * nr * r_bin_mid(n)**2 / rm**3 *     &
                                  EXP( - ( r_bin_mid(n)**3 / rm**3 ) ) *      &
                                  ( r_bin(n+1) - r_bin(n) )
                   ENDDO
                ENDIF
                
!
!--             Criterion to avoid super droplets with a weighting factor < 1.0.                                   
                an_spl = MAX(an_spl, 1.0_wp)
!                
!--             Start splitting operations. Each particle is checked if it
!--             fulfilled the splitting criterion's. In splitting mode 'gb_av'   
!--             a critical radius (radius_split) and a splitting function must 
!--             be prescribed (see particles_par). The critical weighting factor
!--             is calculated while appoximating a 'gamma', 'log' or 'exp'-
!--             drop size distribution. In this mode a DSD is calculated for 
!--             every cloudy grid box. Super droplets which have a larger 
!--             radius and larger weighting factor are split into
!--             'splitting_factor' super droplets. In this case the splitting  
!--             factor is calculated of weighting factor of the super droplet  
!--             and theapproximated number concentration for droplet of such 
!--             a size. Due to the splitting, the weighting factor of the  
!--             super droplet and all created clones is reduced by the factor  
!--             of 'splitting_facor'.
                DO  n = 1, number_of_particles
                   DO  np = 0, n_max-1
                      IF ( r_bin(np) >= radius_split  .AND.                    &
                           particles(n)%particle_mask  .AND.                   &
                           particles(n)%radius >= r_bin(np)    .AND.           &
                           particles(n)%radius < r_bin(np+1)   .AND.           &
                           particles(n)%weight_factor >= an_spl(np) )          &
                      THEN
!
!--                      Calculate splitting factor.
                         splitting_factor =                                    & 
                             MIN( INT( particles(n)%weight_factor /            &
                                        an_spl(np)                             &
                                     ), splitting_factor_max                   &
                                )
                         IF ( splitting_factor < 2 )  CYCLE

!
!--                      Calculate the new number of particles.                                                                                         
                         new_size = prt_count(k,j,i) + splitting_factor - 1
!
!--                      Cycle if maximum number of particles per grid box
!--                      is greater than the allowed maximum number.                                                  
                         IF ( new_size >= max_number_particles_per_gridbox )   &
                         CYCLE
!
!--                      Reallocate particle array if necessary.                          
                         IF ( new_size > SIZE(particles) )  THEN 
                            CALL realloc_particles_array(i,j,k,new_size)
                         ENDIF
!
!--                      Calculate new weighting factor.
                         particles(n)%weight_factor = & 
                            particles(n)%weight_factor / splitting_factor
                         tmp_particle               = particles(n)
                         old_size                   = prt_count(k,j,i)
!
!--                      Create splitting_factor-1 new particles.
                         DO jpp = 1, splitting_factor-1
                            grid_particles(k,j,i)%particles(jpp+old_size) =    &
                               tmp_particle                  
                         ENDDO
!
!--                      Save the new number of droplets for every grid box.
                         prt_count(k,j,i)    = prt_count(k,j,i) +              &
                                               splitting_factor - 1
                         new_particles_gb    = new_particles_gb +              &
                                               splitting_factor - 1                                       
                      ENDIF
                   ENDDO  
                ENDDO

                new_particles       = new_particles + new_particles_gb
                sum_new_particles   = sum_new_particles + new_particles_gb                                                 
             ENDDO
          ENDDO
       ENDDO
    ENDIF
        
    CALL cpu_log( log_point_s(80), 'lpm_splitting', 'stop' )

 END SUBROUTINE lpm_splitting
 
