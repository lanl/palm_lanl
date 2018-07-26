!> @file tests/test-boundaries.f90
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
! $Id: test-boundaries.f90 2718 2018-01-02 08:49:38Z maronga $
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
!> This program tests the boundary grid mode of INIFOR's init_grid_definition()
!> routine.
!------------------------------------------------------------------------------!
 PROGRAM test_boundaries

    USE grid, ONLY  :  init_grid_definition
    USE types, ONLY :  grid_definition
    USE test_utils
    
    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER     ::  title = 'boundary initialization'
    CHARACTER(LEN=20), DIMENSION(2) ::  kind_list = (/ 'boundary', 'boundary' /)
    LOGICAL                         ::  res

    INTEGER                         ::  i, nx, ny, nz
    TYPE(grid_definition)           ::  boundary_grid

    REAL ::  dx, dy, dz, lx, ly, lz, x(2), y(10), z(10)

    CALL begin_test(title, res)

    ! Arange
    dx = 1e-3
    dy = 1.0
    dz = 10.
    nx = 9
    ny = 9
    nz = 9
    lx = 1.0
    ly = 1e1
    lz = 1e2
    x =   (/ -0.5*dx, lx + 0.5*dx /)
    y = ( (/0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0/) + 0.5 )
    z = ( (/0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0/) + 0.5 ) * 10

    DO i = 1, SIZE(kind_list)
   
       ! Act
       CALL init_grid_definition(                                              &
          kind = kind_list(i), grid = boundary_grid,                           &
          xmin = x(i), xmax = x(i),                                            &
          ymin =  0.5 * dy, ymax = ly - 0.5 * dy,                              &
          zmin =  0.5 * dz, zmax = lz - 0.5 * dz,                              &
          x0 = 0.0, y0 = 0.0, z0 = 0.0,                                        &
          nx = 0, ny = ny, nz = nz,                                            &
          dx = dx, dy = dy, dz = dz )
   
       ! Assert
       ! asserting that grid % x has exactly two entries and that they match
       ! expected coordinates
       res = res .AND. assert_equal(boundary_grid % x, (/ x(i) /), 'x coordinates')

       ! asserting that grid % y and % z have expected ranges and coordinates
       res = res .AND. assert_equal( boundary_grid % y, y, 'y coordinates')
       res = res .AND. assert_equal( boundary_grid % z, z, 'z coordinates')
   
       CALL fini_grid_definition(boundary_grid)
    END DO

    CALL end_test(title, res)

 CONTAINS

 SUBROUTINE fini_grid_definition(grid)
    TYPE(grid_definition), INTENT(INOUT) ::  grid

    DEALLOCATE( grid % x, grid % y, grid % z )
    DEALLOCATE( grid % kk )
    DEALLOCATE( grid % w_verti )

 END SUBROUTINE fini_grid_definition

 END PROGRAM test_boundaries
