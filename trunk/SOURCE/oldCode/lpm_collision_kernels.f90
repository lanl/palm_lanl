!> @file lpm_collision_kernels.f90
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
! $Id: lpm_collision_kernels.f90 2718 2018-01-02 08:49:38Z maronga $
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
! 1880 2016-04-20 09:36:50Z hoffmann
! Bugfix: The index of the larger particle has to be chosen for interpolation.
!
! 1873 2016-04-18 14:50:06Z maronga
! Module renamed (removed _mod)
! 
! 1858 2016-04-13 13:12:11Z hoffmann
! Interpolation of collision kernels adjusted to more reasonable values.
! Reformatting of the code. 
!
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
! 
! 1822 2016-04-07 07:49:42Z hoffmann
! PALM kernel has been deleted.
! Bugfix in the calculation of the turbulent enhancement factor of the
! collection efficiency.
!
! Unused variables removed.
!
! 1776 2016-03-02 17:54:58Z hoffmann
! Bugfix: Collection efficiencies must be calculated for the larger droplet.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
! 
! 1519 2015-01-08 10:20:42Z hoffmann
! Bugfix: Using the new particle structure, particles are not sorted by size.
! Hence, computation of collision efficiencies must ensure that the ratio of 
! two colliding droplets is < 1.
!
! 1359 2014-04-11 17:15:14Z hoffmann
! New particle structure integrated.
! Kind definition added to all floating point numbers.
! 
! 1346 2014-03-27 13:18:20Z heinze
! Bugfix: REAL constants provided with KIND-attribute especially in call of 
! intrinsic function like MAX, MIN, SIGN
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL constants defined as wp_kind
!
! 1320 2014-03-20 08:40:49Z 
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
! 1071 2012-11-29 16:54:55Z franke
! Bugfix: collision efficiencies for Hall kernel should not be < 1.0E-20
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1007 2012-09-19 14:30:36Z franke
! converted all units to SI units and replaced some parameters by corresponding
! PALM parameters
! Bugfix: factor in calculation of enhancement factor for collision efficencies
! changed from 10. to 1.0
!
! 849 2012-03-15 10:35:09Z raasch
! routine collision_efficiency_rogers added (moved from former advec_particles
! to here)
!
! 835 2012-02-22 11:21:19Z raasch $
! Bugfix: array diss can be used only in case of Wang kernel
!
! 828 2012-02-21 12:00:36Z raasch
! code has been completely reformatted, routine colker renamed
! recalculate_kernel,
! routine init_kernels added, radius is now communicated to the collision
! routines by array radclass
!
! Bugfix: transformation factor for dissipation changed from 1E5 to 1E4
!
! 825 2012-02-19 03:03:44Z raasch
! routine renamed from wang_kernel to lpm_collision_kernels,
! turbulence_effects on collision replaced by wang_kernel
!
! 790 2011-11-29 03:11:20Z raasch
! initial revision
!
! Description:
! ------------
!> This module calculates collision efficiencies either due to pure gravitational
!> effects (Hall kernel, see Hall, 1980: J. Atmos. Sci., 2486-2507) or
!> including the effects of turbulence (Wang kernel, see Wang and
!> Grabowski, 2009: Atmos. Sci. Lett., 10, 1-8, and Ayala et al., 2008: 
!> New J. Phys., 10, 075016). The original code has been
!> provided by L.-P. Wang but is substantially reformatted and speed optimized
!> here.
!------------------------------------------------------------------------------!
 MODULE lpm_collision_kernels_mod
 

    USE constants,                                                             &
        ONLY:  pi
        
    USE kinds

    USE particle_attributes,                                                   &
        ONLY:  collision_kernel, dissipation_classes, particles,               &
               radius_classes

    USE pegrid


    IMPLICIT NONE

    PRIVATE

    PUBLIC  ckernel, init_kernels, rclass_lbound, rclass_ubound,               &
            recalculate_kernel

    REAL(wp) ::  epsilon       !<
    REAL(wp) ::  rclass_lbound !<
    REAL(wp) ::  rclass_ubound !<
    REAL(wp) ::  urms          !<

    REAL(wp), DIMENSION(:),   ALLOCATABLE ::  epsclass  !< dissipation rate class
    REAL(wp), DIMENSION(:),   ALLOCATABLE ::  radclass  !< radius class
    REAL(wp), DIMENSION(:),   ALLOCATABLE ::  winf      !<
    
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ec        !<
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ecf       !<
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  gck       !<
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  hkernel   !<
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  hwratio   !<
    
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  ckernel !<

    SAVE

!
!-- Public interfaces
    INTERFACE init_kernels
       MODULE PROCEDURE init_kernels
    END INTERFACE init_kernels

    INTERFACE recalculate_kernel
       MODULE PROCEDURE recalculate_kernel
    END INTERFACE recalculate_kernel


    CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the collision efficiency matrix with fixed radius and
!> dissipation classes, calculated at simulation start only.
!------------------------------------------------------------------------------!
 
    SUBROUTINE init_kernels

       IMPLICIT NONE

       INTEGER(iwp) ::  i !<
       INTEGER(iwp) ::  j !<
       INTEGER(iwp) ::  k !<


!
!--    Calculate collision efficiencies for fixed radius- and dissipation
!--    classes
       IF ( collision_kernel(6:9) == 'fast' )  THEN

          ALLOCATE( ckernel(1:radius_classes,1:radius_classes,                 &
                    0:dissipation_classes), epsclass(1:dissipation_classes),   &
                    radclass(1:radius_classes) )

!
!--       Calculate the radius class bounds with logarithmic distances
!--       in the interval [1.0E-6, 1000.0E-6] m
          rclass_lbound = LOG( 1.0E-6_wp )
          rclass_ubound = LOG( 1000.0E-6_wp )
          radclass(1)   = EXP( rclass_lbound )
          DO  i = 2, radius_classes
             radclass(i) = EXP( rclass_lbound +                                &
                                ( rclass_ubound - rclass_lbound ) *            &
                                ( i - 1.0_wp ) / ( radius_classes - 1.0_wp ) )
          ENDDO

!
!--       Set the class bounds for dissipation in interval [0.0, 600.0] cm**2/s**3
          DO  i = 1, dissipation_classes
             epsclass(i) = 0.06_wp * REAL( i, KIND=wp ) / dissipation_classes
          ENDDO
!
!--       Calculate collision efficiencies of the Wang/ayala kernel
          ALLOCATE( ec(1:radius_classes,1:radius_classes),  &
                    ecf(1:radius_classes,1:radius_classes), &
                    gck(1:radius_classes,1:radius_classes), &
                    winf(1:radius_classes) )

          DO  k = 1, dissipation_classes

             epsilon = epsclass(k)
             urms    = 2.02_wp * ( epsilon / 0.04_wp )**( 1.0_wp / 3.0_wp )

             CALL turbsd
             CALL turb_enhance_eff
             CALL effic

             DO  j = 1, radius_classes
                DO  i = 1, radius_classes
                   ckernel(i,j,k) = ec(i,j) * gck(i,j) * ecf(i,j)
                ENDDO
             ENDDO

          ENDDO

!
!--       Calculate collision efficiencies of the Hall kernel
          ALLOCATE( hkernel(1:radius_classes,1:radius_classes), &
                    hwratio(1:radius_classes,1:radius_classes) )

          CALL fallg
          CALL effic

          DO  j = 1, radius_classes
             DO  i =  1, radius_classes
                hkernel(i,j) = pi * ( radclass(j) + radclass(i) )**2 &
                                  * ec(i,j) * ABS( winf(j) - winf(i) )
                ckernel(i,j,0) = hkernel(i,j)  ! hall kernel stored on index 0
              ENDDO
          ENDDO

!
!--       Test output of efficiencies
          IF ( j == -1 )  THEN

             PRINT*, '*** Hall kernel'
             WRITE ( *,'(5X,20(F4.0,1X))' ) ( radclass(i)*1.0E6_wp, &
                                              i = 1,radius_classes )
             DO  j = 1, radius_classes
                WRITE ( *,'(F4.0,1X,20(F8.4,1X))' ) radclass(j),  &
                                          ( hkernel(i,j), i = 1,radius_classes )
             ENDDO

             DO  k = 1, dissipation_classes
                DO  i = 1, radius_classes
                   DO  j = 1, radius_classes
                      IF ( hkernel(i,j) == 0.0_wp )  THEN
                         hwratio(i,j) = 9999999.9_wp
                      ELSE
                         hwratio(i,j) = ckernel(i,j,k) / hkernel(i,j)
                      ENDIF
                   ENDDO
                ENDDO

                PRINT*, '*** epsilon = ', epsclass(k)
                WRITE ( *,'(5X,20(F4.0,1X))' ) ( radclass(i) * 1.0E6_wp, &
                                                 i = 1,radius_classes )
                DO  j = 1, radius_classes
                   WRITE ( *,'(F4.0,1X,20(F8.4,1X))' ) radclass(j) * 1.0E6_wp, &
                                          ( hwratio(i,j), i = 1,radius_classes )
                ENDDO
             ENDDO

          ENDIF

          DEALLOCATE( ec, ecf, epsclass, gck, hkernel, winf )

       ENDIF

    END SUBROUTINE init_kernels


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of collision kernels during each timestep and for each grid box
!------------------------------------------------------------------------------!
    SUBROUTINE recalculate_kernel( i1, j1, k1 )

       USE arrays_3d,                                                          &
           ONLY:  diss

       USE particle_attributes,                                                &
           ONLY:  number_of_particles, prt_count, radius_classes, wang_kernel

       IMPLICIT NONE

       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  i1     !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  j1     !<
       INTEGER(iwp) ::  k1     !<


       number_of_particles = prt_count(k1,j1,i1)
       radius_classes      = number_of_particles   ! necessary to use the same
                                                   ! subroutines as for 
                                                   ! precalculated kernels

       ALLOCATE( ec(1:number_of_particles,1:number_of_particles), &
                 radclass(1:number_of_particles), winf(1:number_of_particles) )

!
!--    Store particle radii on the radclass array
       radclass(1:number_of_particles) = particles(1:number_of_particles)%radius

       IF ( wang_kernel )  THEN
          epsilon = diss(k1,j1,i1)   ! dissipation rate in m**2/s**3
       ELSE
          epsilon = 0.0_wp
       ENDIF
       urms    = 2.02_wp * ( epsilon / 0.04_wp )**( 0.33333333333_wp )

       IF ( wang_kernel  .AND.  epsilon > 1.0E-7_wp )  THEN
!
!--       Call routines to calculate efficiencies for the Wang kernel
          ALLOCATE( gck(1:number_of_particles,1:number_of_particles), &
                    ecf(1:number_of_particles,1:number_of_particles) )

          CALL turbsd
          CALL turb_enhance_eff
          CALL effic

          DO  j = 1, number_of_particles
             DO  i =  1, number_of_particles
                ckernel(1+i-1,1+j-1,1) = ec(i,j) * gck(i,j) * ecf(i,j)
             ENDDO
          ENDDO

          DEALLOCATE( gck, ecf )

       ELSE
!
!--       Call routines to calculate efficiencies for the Hall kernel
          CALL fallg
          CALL effic

          DO  j = 1, number_of_particles
             DO  i =  1, number_of_particles
                ckernel(i,j,1) = pi * ( radclass(j) + radclass(i) )**2         &
                                    * ec(i,j) * ABS( winf(j) - winf(i) )
             ENDDO
          ENDDO

       ENDIF

       DEALLOCATE( ec, radclass, winf )

    END SUBROUTINE recalculate_kernel


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of effects of turbulence on the geometric collision kernel 
!> (by including the droplets' average radial relative velocities and their 
!> radial distribution function) following the analytic model by Aayala et al. 
!> (2008, New J. Phys.). For details check the second part 2 of the publication,
!> page 37ff.
!>
!> Input parameters, which need to be replaced by PALM parameters:
!>    water density, air density
!------------------------------------------------------------------------------!
    SUBROUTINE turbsd

       USE control_parameters,                                                 &
           ONLY:  g, molecular_viscosity
    
       USE particle_attributes,                                                &
           ONLY:  radius_classes

       IMPLICIT NONE

       INTEGER(iwp) ::  i     !<
       INTEGER(iwp) ::  j     !<

       REAL(wp) ::  ao        !<
       REAL(wp) ::  ao_gr     !<
       REAL(wp) ::  bbb       !<
       REAL(wp) ::  be        !<
       REAL(wp) ::  b1        !<
       REAL(wp) ::  b2        !<
       REAL(wp) ::  ccc       !<
       REAL(wp) ::  c1        !<
       REAL(wp) ::  c1_gr     !<
       REAL(wp) ::  c2        !<
       REAL(wp) ::  d1        !<
       REAL(wp) ::  d2        !<
       REAL(wp) ::  eta       !<
       REAL(wp) ::  e1        !<
       REAL(wp) ::  e2        !<
       REAL(wp) ::  fao_gr    !<
       REAL(wp) ::  fr        !<
       REAL(wp) ::  grfin     !<
       REAL(wp) ::  lambda    !<
       REAL(wp) ::  lambda_re !<
       REAL(wp) ::  lf        !<
       REAL(wp) ::  rc        !<
       REAL(wp) ::  rrp       !<
       REAL(wp) ::  sst       !<
       REAL(wp) ::  tauk      !<
       REAL(wp) ::  tl        !<
       REAL(wp) ::  t2        !<
       REAL(wp) ::  tt        !<
       REAL(wp) ::  t1        !<
       REAL(wp) ::  vk        !<
       REAL(wp) ::  vrms1xy   !<
       REAL(wp) ::  vrms2xy   !<
       REAL(wp) ::  v1        !<
       REAL(wp) ::  v1v2xy    !<
       REAL(wp) ::  v1xysq    !<
       REAL(wp) ::  v2        !<
       REAL(wp) ::  v2xysq    !<
       REAL(wp) ::  wrfin     !<
       REAL(wp) ::  wrgrav2   !<
       REAL(wp) ::  wrtur2xy  !<
       REAL(wp) ::  xx        !<
       REAL(wp) ::  yy        !<
       REAL(wp) ::  z         !<

       REAL(wp), DIMENSION(1:radius_classes) ::  st  !< Stokes number
       REAL(wp), DIMENSION(1:radius_classes) ::  tau !< inertial time scale

       lambda    = urms * SQRT( 15.0_wp * molecular_viscosity / epsilon )
       lambda_re = urms**2 * SQRT( 15.0_wp / epsilon / molecular_viscosity )
       tl        = urms**2 / epsilon
       lf        = 0.5_wp * urms**3 / epsilon
       tauk      = SQRT( molecular_viscosity / epsilon )
       eta       = ( molecular_viscosity**3 / epsilon )**0.25_wp
       vk        = eta / tauk

       ao = ( 11.0_wp + 7.0_wp * lambda_re ) / ( 205.0_wp + lambda_re )
       tt = SQRT( 2.0_wp * lambda_re / ( SQRT( 15.0_wp ) * ao ) ) * tauk

!
!--    Get terminal velocity of droplets
       CALL fallg

       DO  i = 1, radius_classes
          tau(i) = winf(i) / g    ! inertial time scale
          st(i)  = tau(i) / tauk  ! Stokes number
       ENDDO

!
!--    Calculate average radial relative velocity at contact (wrfin)
       z   = tt / tl
       be  = SQRT( 2.0_wp ) * lambda / lf
       bbb = SQRT( 1.0_wp - 2.0_wp * be**2 )
       d1  = ( 1.0_wp + bbb ) / ( 2.0_wp * bbb )
       e1  = lf * ( 1.0_wp + bbb ) * 0.5_wp
       d2  = ( 1.0_wp - bbb ) * 0.5_wp / bbb
       e2  = lf * ( 1.0_wp - bbb ) * 0.5_wp
       ccc = SQRT( 1.0_wp - 2.0_wp * z**2 )
       b1  = ( 1.0_wp + ccc ) * 0.5_wp / ccc
       c1  = tl * ( 1.0_wp + ccc ) * 0.5_wp
       b2  = ( 1.0_wp - ccc ) * 0.5_wp / ccc
       c2  = tl * ( 1.0_wp - ccc ) * 0.5_wp

       DO  i = 1, radius_classes

          v1 = winf(i)
          t1 = tau(i)

          DO  j = 1, i
             rrp = radclass(i) + radclass(j)
             v2  = winf(j)
             t2  = tau(j)

             v1xysq  = b1 * d1 * phi_w(c1,e1,v1,t1) - b1 * d2 * phi_w(c1,e2,v1,t1) &
                     - b2 * d1 * phi_w(c2,e1,v1,t1) + b2 * d2 * phi_w(c2,e2,v1,t1)
             v1xysq  = v1xysq * urms**2 / t1
             vrms1xy = SQRT( v1xysq )

             v2xysq  = b1 * d1 * phi_w(c1,e1,v2,t2) - b1 * d2 * phi_w(c1,e2,v2,t2) &
                     - b2 * d1 * phi_w(c2,e1,v2,t2) + b2 * d2 * phi_w(c2,e2,v2,t2)
             v2xysq  = v2xysq * urms**2 / t2
             vrms2xy = SQRT( v2xysq )

             IF ( winf(i) >= winf(j) )  THEN
                v1 = winf(i)
                t1 = tau(i)
                v2 = winf(j)
                t2 = tau(j)
             ELSE
                v1 = winf(j)
                t1 = tau(j)
                v2 = winf(i)
                t2 = tau(i)
             ENDIF

             v1v2xy   =  b1 * d1 * zhi(c1,e1,v1,t1,v2,t2) - &
                         b1 * d2 * zhi(c1,e2,v1,t1,v2,t2) - &
                         b2 * d1 * zhi(c2,e1,v1,t1,v2,t2) + &
                         b2 * d2* zhi(c2,e2,v1,t1,v2,t2)
             fr       = d1 * EXP( -rrp / e1 ) - d2 * EXP( -rrp / e2 )
             v1v2xy   = v1v2xy * fr * urms**2 / tau(i) / tau(j)
             wrtur2xy = vrms1xy**2 + vrms2xy**2 - 2.0_wp * v1v2xy
             IF ( wrtur2xy < 0.0_wp )  wrtur2xy = 0.0_wp
             wrgrav2  = pi / 8.0_wp * ( winf(j) - winf(i) )**2
             wrfin    = SQRT( ( 2.0_wp / pi ) * ( wrtur2xy + wrgrav2) )

!
!--          Calculate radial distribution function (grfin)
             IF ( st(j) > st(i) )  THEN
                sst = st(j)
             ELSE
                sst = st(i)
             ENDIF

             xx = -0.1988_wp * sst**4 + 1.5275_wp * sst**3 - 4.2942_wp *       &
                   sst**2 + 5.3406_wp * sst
             IF ( xx < 0.0_wp )  xx = 0.0_wp
             yy = 0.1886_wp * EXP( 20.306_wp / lambda_re )

             c1_gr  =  xx / ( g / vk * tauk )**yy

             ao_gr  = ao + ( pi / 8.0_wp) * ( g / vk * tauk )**2
             fao_gr = 20.115_wp * SQRT( ao_gr / lambda_re )
             rc     = SQRT( fao_gr * ABS( st(j) - st(i) ) ) * eta

             grfin  = ( ( eta**2 + rc**2 ) / ( rrp**2 + rc**2) )**( c1_gr*0.5_wp )
             IF ( grfin < 1.0_wp )  grfin = 1.0_wp

!
!--          Calculate general collection kernel (without the consideration of
!--          collection efficiencies)
             gck(i,j) = 2.0_wp * pi * rrp**2 * wrfin * grfin
             gck(j,i) = gck(i,j)

          ENDDO
       ENDDO

    END SUBROUTINE turbsd

    REAL(wp) FUNCTION phi_w( a, b, vsett, tau0 )
!
!--    Function used in the Ayala et al. (2008) analytical model for turbulent
!--    effects on the collision kernel
       IMPLICIT NONE

       REAL(wp) ::  a     !<
       REAL(wp) ::  aa1   !<
       REAL(wp) ::  b     !<
       REAL(wp) ::  tau0  !<
       REAL(wp) ::  vsett !<

       aa1 = 1.0_wp / tau0 + 1.0_wp / a + vsett / b
       phi_w = 1.0_wp / aa1  - 0.5_wp * vsett / b / aa1**2

    END FUNCTION phi_w

    REAL(wp) FUNCTION zhi( a, b, vsett1, tau1, vsett2, tau2 )
!
!--    Function used in the Ayala et al. (2008) analytical model for turbulent
!--    effects on the collision kernel
       IMPLICIT NONE

       REAL(wp) ::  a      !<
       REAL(wp) ::  aa1    !<
       REAL(wp) ::  aa2    !<
       REAL(wp) ::  aa3    !<
       REAL(wp) ::  aa4    !<
       REAL(wp) ::  aa5    !<
       REAL(wp) ::  aa6    !<
       REAL(wp) ::  b      !<
       REAL(wp) ::  tau1   !<
       REAL(wp) ::  tau2   !<
       REAL(wp) ::  vsett1 !<
       REAL(wp) ::  vsett2 !<

       aa1 = vsett2 / b - 1.0_wp / tau2 - 1.0_wp / a
       aa2 = vsett1 / b + 1.0_wp / tau1 + 1.0_wp / a
       aa3 = ( vsett1 - vsett2 ) / b + 1.0_wp / tau1 + 1.0_wp / tau2
       aa4 = ( vsett2 / b )**2 - ( 1.0_wp / tau2 + 1.0_wp / a )**2
       aa5 = vsett2 / b + 1.0_wp / tau2 + 1.0_wp / a
       aa6 = 1.0_wp / tau1 - 1.0_wp / a + ( 1.0_wp / tau2 + 1.0_wp / a) *      &
             vsett1 / vsett2
       zhi = (1.0_wp / aa1 - 1.0_wp / aa2 ) * ( vsett1 - vsett2 ) * 0.5_wp /   &
             b / aa3**2 + ( 4.0_wp / aa4 - 1.0_wp / aa5**2 - 1.0_wp / aa1**2 ) &
             * vsett2 * 0.5_wp / b /aa6 + ( 2.0_wp * ( b / aa2 - b / aa1 ) -   &
             vsett1 / aa2**2 + vsett2 / aa1**2 ) * 0.5_wp / b / aa3

    END FUNCTION zhi


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Parameterization of terminal velocity following Rogers et al. (1993, J. Appl.
!> Meteorol.)
!------------------------------------------------------------------------------!
    SUBROUTINE fallg

       USE particle_attributes,                                                &
           ONLY:  radius_classes

       IMPLICIT NONE

       INTEGER(iwp) ::  j                            !<

       REAL(wp), PARAMETER ::  k_cap_rog = 4.0_wp    !< parameter
       REAL(wp), PARAMETER ::  k_low_rog = 12.0_wp   !< parameter
       REAL(wp), PARAMETER ::  a_rog     = 9.65_wp   !< parameter
       REAL(wp), PARAMETER ::  b_rog     = 10.43_wp  !< parameter
       REAL(wp), PARAMETER ::  c_rog     = 0.6_wp    !< parameter
       REAL(wp), PARAMETER ::  d0_rog    = 0.745_wp  !< seperation diameter

       REAL(wp)            ::  diameter              !< droplet diameter in mm


       DO  j = 1, radius_classes

          diameter = radclass(j) * 2000.0_wp

          IF ( diameter <= d0_rog )  THEN
             winf(j) = k_cap_rog * diameter * ( 1.0_wp -                       &
                                                EXP( -k_low_rog * diameter ) )
          ELSE
             winf(j) = a_rog - b_rog * EXP( -c_rog * diameter )
          ENDIF

       ENDDO

    END SUBROUTINE fallg


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolation of collision efficiencies (Hall, 1980, J. Atmos. Sci.)
!------------------------------------------------------------------------------!
    SUBROUTINE effic
 
       USE particle_attributes,                                                &
           ONLY:  radius_classes

       IMPLICIT NONE

       INTEGER(iwp) ::  i  !<
       INTEGER(iwp) ::  iq !<
       INTEGER(iwp) ::  ir !<
       INTEGER(iwp) ::  j  !<
       INTEGER(iwp) ::  k  !<

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  ira !<

       LOGICAL, SAVE ::  first = .TRUE. !<

       REAL(wp) ::  ek              !<
       REAL(wp) ::  particle_radius !<
       REAL(wp) ::  pp              !<
       REAL(wp) ::  qq              !<
       REAL(wp) ::  rq              !<

       REAL(wp), DIMENSION(1:21), SAVE ::  rat        !<
       
       REAL(wp), DIMENSION(1:15), SAVE ::  r0         !<
       
       REAL(wp), DIMENSION(1:15,1:21), SAVE ::  ecoll !<

!
!--    Initial assignment of constants
       IF ( first )  THEN

         first = .FALSE.
         r0  = (/   6.0_wp,   8.0_wp,  10.0_wp, 15.0_wp,  20.0_wp,  25.0_wp,   &
                   30.0_wp,  40.0_wp,  50.0_wp, 60.0_wp,  70.0_wp, 100.0_wp,   &
                  150.0_wp, 200.0_wp, 300.0_wp /)

         rat = (/ 0.00_wp, 0.05_wp, 0.10_wp, 0.15_wp, 0.20_wp, 0.25_wp,        &
                  0.30_wp, 0.35_wp, 0.40_wp, 0.45_wp, 0.50_wp, 0.55_wp,        &
                  0.60_wp, 0.65_wp, 0.70_wp, 0.75_wp, 0.80_wp, 0.85_wp,        &
                  0.90_wp, 0.95_wp, 1.00_wp /)

         ecoll(:,1)  = (/ 0.001_wp, 0.001_wp, 0.001_wp, 0.001_wp, 0.001_wp,    &
                          0.001_wp, 0.001_wp, 0.001_wp, 0.001_wp, 0.001_wp,    &
                          0.001_wp, 0.001_wp, 0.001_wp, 0.001_wp, 0.001_wp /)
         ecoll(:,2)  = (/ 0.003_wp, 0.003_wp, 0.003_wp, 0.004_wp, 0.005_wp,    &
                          0.005_wp, 0.005_wp, 0.010_wp, 0.100_wp, 0.050_wp,    &
                          0.200_wp, 0.500_wp, 0.770_wp, 0.870_wp, 0.970_wp /)
         ecoll(:,3)  = (/ 0.007_wp, 0.007_wp, 0.007_wp, 0.008_wp, 0.009_wp,    &
                          0.010_wp, 0.010_wp, 0.070_wp, 0.400_wp, 0.430_wp,    &
                          0.580_wp, 0.790_wp, 0.930_wp, 0.960_wp, 1.000_wp /)
         ecoll(:,4)  = (/ 0.009_wp, 0.009_wp, 0.009_wp, 0.012_wp, 0.015_wp,    &
                          0.010_wp, 0.020_wp, 0.280_wp, 0.600_wp, 0.640_wp,    &
                          0.750_wp, 0.910_wp, 0.970_wp, 0.980_wp, 1.000_wp /)
         ecoll(:,5)  = (/ 0.014_wp, 0.014_wp, 0.014_wp, 0.015_wp, 0.016_wp,    &
                          0.030_wp, 0.060_wp, 0.500_wp, 0.700_wp, 0.770_wp,    &
                          0.840_wp, 0.950_wp, 0.970_wp, 1.000_wp, 1.000_wp /)
         ecoll(:,6)  = (/ 0.017_wp, 0.017_wp, 0.017_wp, 0.020_wp, 0.022_wp,    &
                          0.060_wp, 0.100_wp, 0.620_wp, 0.780_wp, 0.840_wp,    &
                          0.880_wp, 0.950_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
         ecoll(:,7)  = (/ 0.030_wp, 0.030_wp, 0.024_wp, 0.022_wp, 0.032_wp,    &
                          0.062_wp, 0.200_wp, 0.680_wp, 0.830_wp, 0.870_wp,    &
                          0.900_wp, 0.950_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
         ecoll(:,8)  = (/ 0.025_wp, 0.025_wp, 0.025_wp, 0.036_wp, 0.043_wp,    &
                          0.130_wp, 0.270_wp, 0.740_wp, 0.860_wp, 0.890_wp,    &
                          0.920_wp, 1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
         ecoll(:,9)  = (/ 0.027_wp, 0.027_wp, 0.027_wp, 0.040_wp, 0.052_wp,    &
                          0.200_wp, 0.400_wp, 0.780_wp, 0.880_wp, 0.900_wp,    &
                          0.940_wp, 1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
         ecoll(:,10) = (/ 0.030_wp, 0.030_wp, 0.030_wp, 0.047_wp, 0.064_wp,    &
                          0.250_wp, 0.500_wp, 0.800_wp, 0.900_wp, 0.910_wp,    &
                          0.950_wp, 1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
         ecoll(:,11) = (/ 0.040_wp, 0.040_wp, 0.033_wp, 0.037_wp, 0.068_wp,    &
                          0.240_wp, 0.550_wp, 0.800_wp, 0.900_wp, 0.910_wp,    &
                          0.950_wp, 1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
         ecoll(:,12) = (/ 0.035_wp, 0.035_wp, 0.035_wp, 0.055_wp, 0.079_wp,    &
                          0.290_wp, 0.580_wp, 0.800_wp, 0.900_wp, 0.910_wp,    &
                          0.950_wp, 1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
         ecoll(:,13) = (/ 0.037_wp, 0.037_wp, 0.037_wp, 0.062_wp, 0.082_wp,    &
                          0.290_wp, 0.590_wp, 0.780_wp, 0.900_wp, 0.910_wp,    &
                          0.950_wp, 1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
         ecoll(:,14) = (/ 0.037_wp, 0.037_wp, 0.037_wp, 0.060_wp, 0.080_wp,    &
                          0.290_wp, 0.580_wp, 0.770_wp, 0.890_wp, 0.910_wp,    &
                          0.950_wp, 1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
         ecoll(:,15) = (/ 0.037_wp, 0.037_wp, 0.037_wp, 0.041_wp, 0.075_wp,    &
                          0.250_wp, 0.540_wp, 0.760_wp, 0.880_wp, 0.920_wp,    &
                          0.950_wp, 1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
         ecoll(:,16) = (/ 0.037_wp, 0.037_wp, 0.037_wp, 0.052_wp, 0.067_wp,    &
                          0.250_wp, 0.510_wp, 0.770_wp, 0.880_wp, 0.930_wp,    &
                          0.970_wp, 1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
         ecoll(:,17) = (/ 0.037_wp, 0.037_wp, 0.037_wp, 0.047_wp, 0.057_wp,    &
                          0.250_wp, 0.490_wp, 0.770_wp, 0.890_wp, 0.950_wp,    &
                          1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
         ecoll(:,18) = (/ 0.036_wp, 0.036_wp, 0.036_wp, 0.042_wp, 0.048_wp,    &
                          0.230_wp, 0.470_wp, 0.780_wp, 0.920_wp, 1.000_wp,    &
                          1.020_wp, 1.020_wp, 1.020_wp, 1.020_wp, 1.020_wp /)
         ecoll(:,19) = (/ 0.040_wp, 0.040_wp, 0.035_wp, 0.033_wp, 0.040_wp,    &
                          0.112_wp, 0.450_wp, 0.790_wp, 1.010_wp, 1.030_wp,    &
                          1.040_wp, 1.040_wp, 1.040_wp, 1.040_wp, 1.040_wp /)
         ecoll(:,20) = (/ 0.033_wp, 0.033_wp, 0.033_wp, 0.033_wp, 0.033_wp,    &
                          0.119_wp, 0.470_wp, 0.950_wp, 1.300_wp, 1.700_wp,    &
                          2.300_wp, 2.300_wp, 2.300_wp, 2.300_wp, 2.300_wp /)
         ecoll(:,21) = (/ 0.027_wp, 0.027_wp, 0.027_wp, 0.027_wp, 0.027_wp,    &
                          0.125_wp, 0.520_wp, 1.400_wp, 2.300_wp, 3.000_wp,    &
                          4.000_wp, 4.000_wp, 4.000_wp, 4.000_wp, 4.000_wp /)
       ENDIF

!
!--    Calculate the radius class index of particles with respect to array r
!--    Radius has to be in microns
       ALLOCATE( ira(1:radius_classes) )
       DO  j = 1, radius_classes
          particle_radius = radclass(j) * 1.0E6_wp
          DO  k = 1, 15
             IF ( particle_radius < r0(k) )  THEN
                ira(j) = k
                EXIT
             ENDIF
          ENDDO
          IF ( particle_radius >= r0(15) )  ira(j) = 16
       ENDDO

!
!--    Two-dimensional linear interpolation of the collision efficiency.
!--    Radius has to be in microns
       DO  j = 1, radius_classes
          DO  i = 1, j

             ir = MAX( ira(i), ira(j) )
             rq = MIN( radclass(i) / radclass(j), radclass(j) / radclass(i) )
             iq = INT( rq * 20 ) + 1
             iq = MAX( iq , 2)

             IF ( ir < 16 )  THEN
                IF ( ir >= 2 )  THEN
                   pp = ( ( MAX( radclass(j), radclass(i) ) * 1.0E6_wp ) -     &
                          r0(ir-1) ) / ( r0(ir) - r0(ir-1) )
                   qq = ( rq - rat(iq-1) ) / ( rat(iq) - rat(iq-1) )
                   ec(j,i) = ( 1.0_wp - pp ) * ( 1.0_wp - qq )                 &
                             * ecoll(ir-1,iq-1)                                &
                             + pp * ( 1.0_wp - qq ) * ecoll(ir,iq-1)           &
                             + qq * ( 1.0_wp - pp ) * ecoll(ir-1,iq)           &
                             + pp * qq * ecoll(ir,iq)
                ELSE
                   qq = ( rq - rat(iq-1) ) / ( rat(iq) - rat(iq-1) )
                   ec(j,i) = ( 1.0_wp - qq ) * ecoll(1,iq-1) + qq * ecoll(1,iq)
                ENDIF
             ELSE
                qq = ( rq - rat(iq-1) ) / ( rat(iq) - rat(iq-1) )
                ek = ( 1.0_wp - qq ) * ecoll(15,iq-1) + qq * ecoll(15,iq)
                ec(j,i) = MIN( ek, 1.0_wp )
             ENDIF

             IF ( ec(j,i) < 1.0E-20_wp )  ec(j,i) = 0.0_wp

             ec(i,j) = ec(j,i)

          ENDDO
       ENDDO

       DEALLOCATE( ira )

    END SUBROUTINE effic


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolation of turbulent enhancement factor for collision efficencies
!> following Wang and Grabowski (2009, Atmos. Sci. Let.)
!------------------------------------------------------------------------------!
    SUBROUTINE turb_enhance_eff

       USE particle_attributes,                                                &
           ONLY:  radius_classes

       IMPLICIT NONE

       INTEGER(iwp) :: i  !<
       INTEGER(iwp) :: iq !<
       INTEGER(iwp) :: ir !<
       INTEGER(iwp) :: j  !<
       INTEGER(iwp) :: k  !<
       INTEGER(iwp) :: kk !<

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  ira !<
       
       LOGICAL, SAVE ::  first = .TRUE. !<

       REAL(wp) ::  particle_radius !<
       REAL(wp) ::  pp              !<
       REAL(wp) ::  qq              !<
       REAL(wp) ::  rq              !<
       REAL(wp) ::  y1              !<
       REAL(wp) ::  y2              !<
       REAL(wp) ::  y3              !<

       REAL(wp), DIMENSION(1:11), SAVE ::  rat           !<
       REAL(wp), DIMENSION(1:7), SAVE  ::  r0            !<
       
       REAL(wp), DIMENSION(1:7,1:11), SAVE ::  ecoll_100 !<
       REAL(wp), DIMENSION(1:7,1:11), SAVE ::  ecoll_400 !<

!
!--    Initial assignment of constants
       IF ( first )  THEN

          first = .FALSE.

          r0  = (/  10.0_wp, 20.0_wp, 30.0_wp, 40.0_wp, 50.0_wp, 60.0_wp,  &
                   100.0_wp /)

          rat = (/ 0.0_wp, 0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp, 0.5_wp, 0.6_wp, &
                   0.7_wp, 0.8_wp, 0.9_wp, 1.0_wp /)
!
!--       Tabulated turbulent enhancement factor at 100 cm**2/s**3
          ecoll_100(:,1)  = (/  1.74_wp,   1.74_wp,   1.773_wp, 1.49_wp,  &
                                1.207_wp,  1.207_wp,  1.0_wp /)
          ecoll_100(:,2)  = (/  1.46_wp,   1.46_wp,   1.421_wp, 1.245_wp, &
                                1.069_wp,  1.069_wp,  1.0_wp /)
          ecoll_100(:,3)  = (/  1.32_wp,   1.32_wp,   1.245_wp, 1.123_wp, &
                                1.000_wp,  1.000_wp,  1.0_wp /)
          ecoll_100(:,4)  = (/  1.250_wp,  1.250_wp,  1.148_wp, 1.087_wp, &
                                1.025_wp,  1.025_wp,  1.0_wp /)
          ecoll_100(:,5)  = (/  1.186_wp,  1.186_wp,  1.066_wp, 1.060_wp, &
                                1.056_wp,  1.056_wp,  1.0_wp /)
          ecoll_100(:,6)  = (/  1.045_wp,  1.045_wp,  1.000_wp, 1.014_wp, &
                                1.028_wp,  1.028_wp,  1.0_wp /)
          ecoll_100(:,7)  = (/  1.070_wp,  1.070_wp,  1.030_wp, 1.038_wp, &
                                1.046_wp,  1.046_wp,  1.0_wp /)
          ecoll_100(:,8)  = (/  1.000_wp,  1.000_wp,  1.054_wp, 1.042_wp, &
                                1.029_wp,  1.029_wp,  1.0_wp /)
          ecoll_100(:,9)  = (/  1.223_wp,  1.223_wp,  1.117_wp, 1.069_wp, &
                                1.021_wp,  1.021_wp,  1.0_wp /)
          ecoll_100(:,10) = (/  1.570_wp,  1.570_wp,  1.244_wp, 1.166_wp, &
                                1.088_wp,  1.088_wp,  1.0_wp /)
          ecoll_100(:,11) = (/ 20.3_wp,   20.3_wp,   14.6_wp,   8.61_wp,  &
                                2.60_wp,   2.60_wp,   1.0_wp /)
!
!--       Tabulated turbulent enhancement factor at 400 cm**2/s**3
          ecoll_400(:,1)  = (/  4.976_wp,  4.976_wp,  3.593_wp,  2.519_wp, &
                                1.445_wp,  1.445_wp,  1.0_wp /)
          ecoll_400(:,2)  = (/  2.984_wp,  2.984_wp,  2.181_wp,  1.691_wp, &
                                1.201_wp,  1.201_wp,  1.0_wp /)
          ecoll_400(:,3)  = (/  1.988_wp,  1.988_wp,  1.475_wp,  1.313_wp, &
                                1.150_wp,  1.150_wp,  1.0_wp /)
          ecoll_400(:,4)  = (/  1.490_wp,  1.490_wp,  1.187_wp,  1.156_wp, &
                                1.126_wp,  1.126_wp,  1.0_wp /)
          ecoll_400(:,5)  = (/  1.249_wp,  1.249_wp,  1.088_wp,  1.090_wp, &
                                1.092_wp,  1.092_wp,  1.0_wp /)
          ecoll_400(:,6)  = (/  1.139_wp,  1.139_wp,  1.130_wp,  1.091_wp, &
                                1.051_wp,  1.051_wp,  1.0_wp /)
          ecoll_400(:,7)  = (/  1.220_wp,  1.220_wp,  1.190_wp,  1.138_wp, &
                                1.086_wp,  1.086_wp,  1.0_wp /)
          ecoll_400(:,8)  = (/  1.325_wp,  1.325_wp,  1.267_wp,  1.165_wp, &
                                1.063_wp,  1.063_wp,  1.0_wp /)
          ecoll_400(:,9)  = (/  1.716_wp,  1.716_wp,  1.345_wp,  1.223_wp, &
                                1.100_wp,  1.100_wp,  1.0_wp /)
          ecoll_400(:,10) = (/  3.788_wp,  3.788_wp,  1.501_wp,  1.311_wp, &
                                1.120_wp,  1.120_wp,  1.0_wp /)
          ecoll_400(:,11) = (/ 36.52_wp,  36.52_wp,  19.16_wp,  22.80_wp,  &
                               26.0_wp,   26.0_wp,    1.0_wp /)

       ENDIF

!
!--    Calculate the radius class index of particles with respect to array r0
!--    The droplet radius has to be given in microns.
       ALLOCATE( ira(1:radius_classes) )

       DO  j = 1, radius_classes
          particle_radius = radclass(j) * 1.0E6_wp
          DO  k = 1, 7
             IF ( particle_radius < r0(k) )  THEN
                ira(j) = k
                EXIT
             ENDIF
          ENDDO
          IF ( particle_radius >= r0(7) )  ira(j) = 8
       ENDDO

!
!--    Two-dimensional linear interpolation of the turbulent enhancement factor.
!--    The droplet radius has to be given in microns.
       DO  j =  1, radius_classes
          DO  i = 1, j

             ir = MAX( ira(i), ira(j) )
             rq = MIN( radclass(i) / radclass(j), radclass(j) / radclass(i) )

             DO  kk = 2, 11
                IF ( rq <= rat(kk) )  THEN
                   iq = kk
                   EXIT
                ENDIF
             ENDDO

             y1 = 1.0_wp  ! turbulent enhancement factor at 0 m**2/s**3

             IF ( ir < 8 )  THEN
                IF ( ir >= 2 )  THEN
                   pp = ( MAX( radclass(j), radclass(i) ) * 1.0E6_wp -  &
                          r0(ir-1) ) / ( r0(ir) - r0(ir-1) )
                   qq = ( rq - rat(iq-1) ) / ( rat(iq) - rat(iq-1) )
                   y2 = ( 1.0_wp - pp ) * ( 1.0_wp - qq ) * ecoll_100(ir-1,iq-1) + &
                                pp * ( 1.0_wp - qq ) * ecoll_100(ir,iq-1)        + &
                                qq * ( 1.0_wp - pp ) * ecoll_100(ir-1,iq)        + &
                                pp * qq              * ecoll_100(ir,iq)
                   y3 = ( 1.0-pp ) * ( 1.0_wp - qq ) * ecoll_400(ir-1,iq-1)      + &
                                pp * ( 1.0_wp - qq ) * ecoll_400(ir,iq-1)        + &
                                qq * ( 1.0_wp - pp ) * ecoll_400(ir-1,iq)        + &
                                pp * qq              * ecoll_400(ir,iq)
                ELSE
                   qq = ( rq - rat(iq-1) ) / ( rat(iq) - rat(iq-1) )
                   y2 = ( 1.0_wp - qq ) * ecoll_100(1,iq-1) + qq * ecoll_100(1,iq)
                   y3 = ( 1.0_wp - qq ) * ecoll_400(1,iq-1) + qq * ecoll_400(1,iq)
                ENDIF
             ELSE
                qq = ( rq - rat(iq-1) ) / ( rat(iq) - rat(iq-1) )
                y2 = ( 1.0_wp - qq ) * ecoll_100(7,iq-1) + qq * ecoll_100(7,iq)
                y3 = ( 1.0_wp - qq ) * ecoll_400(7,iq-1) + qq * ecoll_400(7,iq)
             ENDIF
!
!--          Linear interpolation of turbulent enhancement factor
             IF ( epsilon <= 0.01_wp )  THEN
                ecf(j,i) = ( epsilon - 0.01_wp ) / ( 0.0_wp  - 0.01_wp ) * y1 &
                         + ( epsilon - 0.0_wp  ) / ( 0.01_wp - 0.0_wp  ) * y2
             ELSEIF ( epsilon <= 0.06_wp )  THEN
                ecf(j,i) = ( epsilon - 0.04_wp ) / ( 0.01_wp - 0.04_wp ) * y2 &
                         + ( epsilon - 0.01_wp ) / ( 0.04_wp - 0.01_wp ) * y3
             ELSE
                ecf(j,i) = ( 0.06_wp - 0.04_wp ) / ( 0.01_wp - 0.04_wp ) * y2 &
                         + ( 0.06_wp - 0.01_wp ) / ( 0.04_wp - 0.01_wp ) * y3
             ENDIF

             IF ( ecf(j,i) < 1.0_wp )  ecf(j,i) = 1.0_wp

             ecf(i,j) = ecf(j,i)

          ENDDO
       ENDDO

    END SUBROUTINE turb_enhance_eff

 END MODULE lpm_collision_kernels_mod
