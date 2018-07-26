!> @file tests/test-interpolation.f90
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
! $Id: test-interpolation.f90 2718 2018-01-02 08:49:38Z maronga $
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
!> This program tests INIFOR's horizontal interpolation.
!------------------------------------------------------------------------------!
 PROGRAM test_interpolation

    USE grid, &
        ONLY:  grid_definition, init_grid_definition, TO_RADIANS, TO_DEGREES, &
               linspace, hhl
    USE transform, &
        ONLY:  find_horizontal_neighbours, compute_horizontal_interp_weights
    USE test_utils
    
    IMPLICIT NONE

!
!------------------------------------------------------------------------------
!- Test 1: Find neighbours
!------------------------------------------------------------------------------
    CHARACTER(LEN=30) ::  title = "find neighbours"
    LOGICAL           ::  res

    TYPE(grid_definition) ::  palm_grid, cosmo_grid
    INTEGER               ::  i, j, ii_ref(0:1, 0:1, 4), jj_ref(0:1, 0:1, 4)
    INTEGER, PARAMETER    ::  nlon=3, nlat=3, nlev=2
    REAL                  ::  w_ref(4), lat(0:2), lon(0:2)

    title = "find neighbours"
    CALL begin_test(title, res)

    ! Arange.
    ! Make a COSMO-DE grid with just two horizotal cells/three h. points
    PRINT *, "INIT GRID"

    ! Allocate grid.hhl for use in init_grid_definition. In INIFOR, this is done
    ! in get_netcdf_variable_2d. In this test, grid.hhl is not used and only
    ! defined manually because it is used in init_grid_definition. 
    ALLOCATE (hhl (nlon, nlat, nlev) )
    hhl(:,:,:) = 0.0
    CALL init_grid_definition('cosmo-de', grid = cosmo_grid,                   &
                              xmin = -5.0 * TO_RADIANS, xmax = 5.5 * TO_RADIANS, &
                              ymin = -5.0 * TO_RADIANS, ymax = 6.5 * TO_RADIANS, &
                              zmin =  0.0, zmax = 10.0,                        &
                              x0 = 0.0, y0 = 0.0, z0 = 0.0,                    &
                              nx = nlon-1, ny = nlat-1, nz = nlev-1)

    PRINT *, "GRID DONE"
    PRINT *, "COSMO lats: ", cosmo_grid % lat * TO_DEGREES
    PRINT *, "COSMO lons: ", cosmo_grid % lon * TO_DEGREES
    
    res = assert_equal( (/cosmo_grid%lat(0), cosmo_grid % lon(0),              &
                          cosmo_grid%lat(2), cosmo_grid % lon(2),              &
                          cosmo_grid%dx*TO_DEGREES, cosmo_grid%dy*TO_DEGREES/),&
                        (/-5.0 * TO_RADIANS, -5.0 * TO_RADIANS,                &
                           6.5 * TO_RADIANS,  5.5 * TO_RADIANS,                &
                           5.25, 5.75 /),    &
                        "COSMO grid coordinates" )
    ! Define a PALM-4U grid with only one cell, i.e. four points in the
    ! horizontal plane. The points are located at the centres of
    ! the COSMO-DE cells.
    CALL init_grid_definition('palm intermediate', grid = palm_grid,           &
                              xmin = 0.0, xmax = 1.0,                          &
                              ymin = 0.0, ymax = 1.0,                          &
                              zmin = 0.0, zmax = 1.0,                          &
                              x0 = 0.0, y0 = 0.0, z0 = 0.0,                    &
                              nx = 1, ny = 1, nz = 1)

    palm_grid % clon(0,0) = 0.5 * cosmo_grid % lon(0)
    palm_grid % clat(0,0) = 0.5 * cosmo_grid % lat(0)
    
    palm_grid % clon(0,1) = 0.5 * cosmo_grid % lon(0)
    palm_grid % clat(0,1) = 0.5 * cosmo_grid % lat(2)

    palm_grid % clon(1,1) = 0.5 * cosmo_grid % lon(2)
    palm_grid % clat(1,1) = 0.5 * cosmo_grid % lat(2)

    palm_grid % clon(1,0) = 0.5 * cosmo_grid % lon(2)
    palm_grid % clat(1,0) = 0.5 * cosmo_grid % lat(0)

    ii_ref(0,0,:) = (/0, 0, 1, 1/)
    jj_ref(0,0,:) = (/0, 1, 1, 0/)

    ii_ref(0,1,:) = (/0, 0, 1, 1/)
    jj_ref(0,1,:) = (/1, 2, 2, 1/)
    
    ii_ref(1,1,:) = (/1, 1, 2, 2/)
    jj_ref(1,1,:) = (/1, 2, 2, 1/)

    ii_ref(1,0,:) = (/1, 1, 2, 2/)
    jj_ref(1,0,:) = (/0, 1, 1, 0/)

    ! Act
    CALL find_horizontal_neighbours(cosmo_grid % lat, cosmo_grid % lon,        &
       cosmo_grid % dxi, cosmo_grid % dyi, palm_grid % clat, palm_grid % clon, &
       palm_grid % ii, palm_grid % jj)

    ! Assert
    DO j = 0, 1
    DO i = 0, 1
       res = res .AND. ALL(palm_grid%ii(i,j,:) == ii_ref(i,j,:))
       PRINT *, "ii     : ", palm_grid%ii(i,j,:)
       PRINT *, "ii_ref : ", ii_ref(i,j,:), " indices match? ", res
       res = res .AND. ALL(palm_grid%jj(i,j,:) == jj_ref(i,j,:))
       PRINT *, "jj     : ", palm_grid%jj(i,j,:)
       PRINT *, "jj_ref : ", jj_ref(i,j,:), " indices match? ", res
    END DO
    END DO

    CALL end_test(title, res)
    

!
!------------------------------------------------------------------------------
!- Test 2: Compute weights for linear interpolation
!------------------------------------------------------------------------------
    title = "interpolation weights"
    CALL begin_test(title, res)

    ! Arange
    ! defining some shorthands
    lon(:) = cosmo_grid % lon(:)
    lat(:) = cosmo_grid % lat(:)

    ! set up PALM-4U points at 1/4 and 1/3 of the COSMO grid widths
    palm_grid % clon(0,0) = -0.25  * (lon(1) - lon(0)) + lon(1)
    palm_grid % clat(0,0) = -2./3. * (lat(1) - lat(0)) + lat(1)
    
    palm_grid % clon(0,1) = -2./3. * (lon(1) - lon(0)) + lon(1)
    palm_grid % clat(0,1) = +0.25  * (lat(2) - lat(1)) + lat(1)

    palm_grid % clon(1,1) = +0.25  * (lon(2) - lon(1)) + lon(1)
    palm_grid % clat(1,1) = +2./3. * (lat(2) - lat(1)) + lat(1)

    palm_grid % clon(1,0) = +2./3. * (lon(2) - lon(1)) + lon(1)
    palm_grid % clat(1,0) = -0.25  * (lat(1) - lat(0)) + lat(1)

    DO j = 0, 1
    DO i = 0, 1
       PRINT *, "PALM lon, lat: ", palm_grid % clon(i,j) * TO_DEGREES, palm_grid % clat(i,j)*TO_DEGREES
    END DO
    END DO

    ! Act
    CALL find_horizontal_neighbours(cosmo_grid % lat, cosmo_grid % lon,        &
       cosmo_grid % dxi, cosmo_grid % dyi, palm_grid % clat, palm_grid % clon, &
       palm_grid % ii, palm_grid % jj)

    CALL compute_horizontal_interp_weights(cosmo_grid % lat, cosmo_grid % lon, &
       cosmo_grid % dxi, cosmo_grid % dyi, palm_grid % clat,                   &
       palm_grid % clon, palm_grid % ii, palm_grid % jj, palm_grid % w_horiz)

    ! Assert
    ! asserting that neighbours are still correct
    DO j = 0, 1
    DO i = 0, 1
       res = res .AND. ALL(palm_grid%ii(i,j,:) == ii_ref(i,j,:))
       PRINT *, "ii     : ", palm_grid%ii(i,j,:)
       PRINT *, "ii_ref : ", ii_ref(i,j,:), " indices match? ", res
       res = res .AND. ALL(palm_grid%jj(i,j,:) == jj_ref(i,j,:))
       PRINT *, "jj     : ", palm_grid%jj(i,j,:)
       PRINT *, "jj_ref : ", jj_ref(i,j,:), " indices match? ", res
    END DO
    END DO

    ! asserting that all four weights equal, 0.5, 0.25, 1./6., and 1./12., resp.
    w_ref = (/1./6., 1./12., 0.25, 0.5/)
    res = res .AND. assert_equal(palm_grid % w_horiz(0, 0, :), w_ref(:), "weights at (0,0)")
    !res = res .AND. palm_grid % w_horiz(0, 0, 1) == w_ref(1)
    !res = res .AND. palm_grid % w_horiz(0, 0, 2) == w_ref(2)
    !res = res .AND. palm_grid % w_horiz(0, 0, 3) == w_ref(3)
    !res = res .AND. palm_grid % w_horiz(0, 0, 4) == w_ref(4)

    w_ref = (/0.5, 1./6., 1./12., 0.25/)
    res = res .AND. assert_equal(palm_grid % w_horiz(0, 1, :), w_ref(:), "weights at (0,1)")
    !res = res .AND. palm_grid % w_horiz(0, 1, 1) == w_ref(4)
    !res = res .AND. palm_grid % w_horiz(0, 1, 2) == w_ref(1)
    !res = res .AND. palm_grid % w_horiz(0, 1, 3) == w_ref(2)
    !res = res .AND. palm_grid % w_horiz(0, 1, 4) == w_ref(3)

    w_ref = (/0.25, 0.5, 1./6., 1./12./)
    res = res .AND. assert_equal(palm_grid % w_horiz(1, 1, :), w_ref(:), "weights at (1,1)")
    !res = res .AND. palm_grid % w_horiz(1, 1, 1) == w_ref(3)
    !res = res .AND. palm_grid % w_horiz(1, 1, 2) == w_ref(4)
    !res = res .AND. palm_grid % w_horiz(1, 1, 3) == w_ref(1)
    !res = res .AND. palm_grid % w_horiz(1, 1, 4) == w_ref(2)

    w_ref = (/1./12., 0.25, 0.5, 1./6./)
    res = res .AND. assert_equal(palm_grid % w_horiz(1, 0, :), w_ref(:), "weights at (1,0)")
    !res = res .AND. palm_grid % w_horiz(1, 0, 1) == w_ref(2)
    !res = res .AND. palm_grid % w_horiz(1, 0, 2) == w_ref(3)
    !res = res .AND. palm_grid % w_horiz(1, 0, 3) == w_ref(4)
    !res = res .AND. palm_grid % w_horiz(1, 0, 4) == w_ref(1)

    CALL end_test(title, res)

 END PROGRAM test_interpolation
