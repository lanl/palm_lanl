!> @file tests/test-grid.f90
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
! Copyright 2017-2018 Leibniz Universitaet Hannover
! Copyright 2017-2018 Deutscher Wetterdienst Offenbach
!------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: test-grid.f90 2718 2018-01-02 08:49:38Z maronga $
! Initial revision
!
! 
!
! Authors:
! --------
! @author Eckhard Kadasch
!
! Description:
! ------------
!> This program tests the PALM grid mode of INIFOR's init_grid_definition()
!> routine.
!------------------------------------------------------------------------------!
 PROGRAM test_grid

    USE grid, ONLY :  grid_definition, init_grid_definition
    USE test_utils
    
    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER ::  title = "grid initialization"
    LOGICAL                     ::  res

    TYPE(grid_definition) ::  mygrid
    INTEGER               ::  i
    INTEGER, PARAMETER    ::  nx = 9,   ny = 19,   nz = 29
    REAL, PARAMETER       ::  lx = 100., ly = 200., lz = 300.
    REAL, DIMENSION(0:nx) ::  x, xu
    REAL, DIMENSION(0:ny) ::  y, yv
    REAL, DIMENSION(0:nz) ::  z, zw

    CALL begin_test(title, res)

    ! Arange
    CALL init_grid_definition('palm', grid = mygrid,                           &
                              xmin = 0., xmax = lx,                            &
                              ymin = 0., ymax = ly,                            &
                              zmin = 0., zmax = lz,                            &
                              x0 = 0.0, y0 = 0.0, z0 = 0.0,                    &
                              nx = nx, ny = ny, nz = nz)

    ! Act
    DO i = 0, nx
       xu(i) = real(i) / (nx+1) * lx
       x(i)  = 0.5*mygrid%dx + xu(i)
    END DO
    DO i = 0, ny
       yv(i) = real(i) / (ny+1) * ly
       y(i)  = 0.5*mygrid%dy + yv(i)
    END DO
    DO i = 0, nz
       zw(i) = real(i) / (nz+1) * lz
       z(i)  = 0.5*mygrid%dz + zw(i)
    END DO

    ! Assert coordinates match
    res = res .AND. assert_equal(x,      mygrid%x,  "x" )
    res = res .AND. assert_equal(xu(1:), mygrid%xu, "xu")
    res = res .AND. assert_equal(y,      mygrid%y,  "y" )
    res = res .AND. assert_equal(yv(1:), mygrid%yv, "yu")
    res = res .AND. assert_equal(z,      mygrid%z,  "z" )
    res = res .AND. assert_equal(zw(1:), mygrid%zw, "zu")

    CALL end_test(title, res)

 END PROGRAM test_grid
