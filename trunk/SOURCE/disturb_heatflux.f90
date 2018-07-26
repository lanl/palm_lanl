!> @file disturb_heatflux.f90
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
! $Id: disturb_heatflux.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustment according new surface data type.
! Implemented parallel random number generator to obtain always the same 
! random number distribution regardless of the processor distribution.
! 
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute 
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL constants defined as wp-kind
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements, 
! kinds are defined in new module kinds, 
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! Revision 1.1  1998/03/25 20:03:47  raasch
! Initial revision
!
!
! Description:
! ------------
!> Generate random, normally distributed heatflux values and store them as the
!> near-surface heatflux.
!------------------------------------------------------------------------------!
 SUBROUTINE disturb_heatflux( surf )
 

    USE arrays_3d,                                                             &
        ONLY:  heatflux_input_conversion
        
    USE control_parameters,                                                    &
        ONLY:  iran, surface_heatflux, random_generator, wall_heatflux
        
    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point
        
    USE kinds
    
    USE indices,                                                               &
        ONLY:  nzb

    USE random_generator_parallel,                                             &
        ONLY:  random_number_parallel, random_seed_parallel, random_dummy,     &
               id_random_array, seq_random_array

    USE surface_mod,                                                           &
        ONLY:  surf_type

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< grid index, x direction
    INTEGER(iwp) ::  j  !< grid index, y direction
    INTEGER(iwp) ::  k  !< grid index, z direction
    INTEGER(iwp) ::  m  !< loop variables over surface elements
    
    REAL(wp) ::  random_gauss  !<
    REAL(wp) ::  randomnumber  !<

    TYPE(surf_type) ::  surf   !< surface-type variable


    CALL cpu_log( log_point(23), 'disturb_heatflux', 'start' )

!
!-- Generate random disturbances and store them. Note, if 
!-- random_generator /= 'random-parallel' it is not guaranteed to obtain  
!-- the same random distribution if the number of processors is changed. 
    IF ( random_generator /= 'random-parallel' )  THEN

       DO  m = 1, surf%ns

          k = surf%k(m)

          randomnumber = random_gauss( iran, 5.0_wp )      
!
!--       k-1 is topography top index. If this is 0, set surface heatflux. Over
!--       topography surface_heatflux is replaced by wall_heatflux(0).
          IF ( k-1 == 0 )  THEN
             surf%shf(m) = randomnumber * surface_heatflux                     &
                              * heatflux_input_conversion(nzb)
          ELSE
             surf%shf(m) = randomnumber * wall_heatflux(0)                     &
                              * heatflux_input_conversion(k-1)
          ENDIF
       ENDDO
    ELSE

       DO  m = 1, surf%ns

          i = surf%i(m)
          j = surf%j(m)
          k = surf%k(m)

          CALL random_seed_parallel( put=seq_random_array(:, j, i) )
          CALL random_number_parallel( random_dummy )    
!
!--       k-1 is topography top index. If this is 0, set surface heatflux. Over
!--       topography surface_heatflux is replaced by wall_heatflux(0).
          IF ( k-1 == 0 )  THEN
             surf%shf(m) = ( random_dummy - 0.5_wp ) * surface_heatflux        &
                              * heatflux_input_conversion(nzb)
          ELSE
             surf%shf(m) = ( random_dummy - 0.5_wp ) * wall_heatflux(0)        &
                              * heatflux_input_conversion(k-1)
          ENDIF

          CALL random_seed_parallel( get=seq_random_array(:, j, i) )

       ENDDO

    ENDIF

    CALL cpu_log( log_point(23), 'disturb_heatflux', 'stop' )


 END SUBROUTINE disturb_heatflux
