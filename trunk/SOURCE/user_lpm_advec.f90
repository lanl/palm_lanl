!> @file user_lpm_advec.f90
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
! $Id: user_lpm_advec.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2417 2017-09-06 15:22:27Z suehring
! Particle loops adapted for sub-box structure, i.e. for each sub-box the 
! particle loop runs from start_index up to end_index instead from 1 to 
! number_of_particles.
! 
! 2101 2017-01-05 16:42:31Z suehring
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1416 2014-06-04 16:04:03Z suehring
! Provide template for new particle structure.
! Moreover, user_lpm_advec is called for each gridpoint.
! 
! 1359 2014-04-11 17:15:14Z hoffmann
! New particle structure integrated. 
!
! 1320 2014-03-20 08:40:49Z raasch
! kind-parameters added to all INTEGER and REAL declaration statements, 
! kinds are defined in new module kinds, 
! old module precision_kind is removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 849 2012-03-15 10:35:09Z raasch
! routine renamed user_lpm_advec
!
! 211 2008-11-11 04:46:24Z raasch
! Former file user_interface.f90 split into one file per subroutine
!
! Description:
! ------------
!> Modification of initial particles by the user.
!------------------------------------------------------------------------------!
 SUBROUTINE user_lpm_advec( ip, jp, kp )
 

    USE kinds
    
    USE particle_attributes
    
    USE user

    IMPLICIT NONE

    INTEGER(iwp) ::  ip   !< index of particle grid box, x-direction
    INTEGER(iwp) ::  jp   !< index of particle grid box, y-direction
    INTEGER(iwp) ::  kp   !< index of particle grid box, z-direction
    INTEGER(iwp) ::  n    !< particle index
    INTEGER(iwp) ::  nb   !< index of sub-box particles are sorted in

    INTEGER(iwp), DIMENSION(0:7)  ::  start_index !< start particle index for current sub-box
    INTEGER(iwp), DIMENSION(0:7)  ::  end_index   !< start particle index for current sub-box

!
!-- Here the user-defined actions follow
!   number_of_particles = prt_count(kp,jp,ip)
!   particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
!
!   start_index = grid_particles(kp,jp,ip)%start_index
!   end_index   = grid_particles(kp,jp,ip)%end_index
!
!   IF ( number_of_particles <= 0 )  CYCLE
!   DO  n = 1, number_of_particles
!   DO  nb = 0, 7
!      DO  n = start_index(nb), end_index(nb)
!         particles(n)%xxx = 
!      ENDDO
!   ENDDO

 END SUBROUTINE user_lpm_advec

