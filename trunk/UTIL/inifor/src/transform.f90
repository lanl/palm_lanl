!> @file src/transform.f90
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
! $Id: transform.f90 2718 2018-01-02 08:49:38Z maronga $
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
!> The transform module provides INIFOR's low-level transformation and
!> interpolation routines. The rotated-pole transformation routines phirot2phi,
!> phi2phirot, rlarot2rla, rla2rlarot, uv2uvrot, and uvrot2uv are adapted from
!> int2lm's utility routines.
!------------------------------------------------------------------------------!
 MODULE transform

    USE control
    USE defs,                                                                  &
        ONLY: TO_DEGREES, TO_RADIANS, PI, dp
    USE types
    USE util,                                                                  &       
        ONLY: real_to_str, str

    IMPLICIT NONE

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolates linearly in the vertical direction in very column (i,j) of the
!> output array outvar(i,j,:) using values of the source array invar. In cells
!> that are outside the COSMO-DE domain, indicated by negative interpolation
!> weights, extrapolate constantly from the cell above.
!> 
!> Input parameters:
!> -----------------
!> invar : Array of source data
!> 
!> outgrid % kk : Array of vertical neighbour indices. kk(i,j,k,:) contain the
!>     indices of the two vertical neighbors of PALM-4U point (i,j,k) on the
!>     input grid corresponding to the source data invar.
!> 
!> outgrid % w_verti : Array of weights for vertical linear interpolation
!>     corresponding to neighbour points indexed by kk.
!>
!> Output papameters:
!> ------------------
!> outvar : Array of interpolated data
!------------------------------------------------------------------------------!
    SUBROUTINE interpolate_1d(in_arr, out_arr, outgrid)
       TYPE(grid_definition), INTENT(IN) ::  outgrid
       REAL(dp), INTENT(IN)              ::  in_arr(0:,0:,0:)
       REAL(dp), INTENT(OUT)             ::  out_arr(0:,0:,0:)

       INTEGER :: i, j, k, l, nx, ny, nz

       nx = UBOUND(out_arr, 1)
       ny = UBOUND(out_arr, 2)
       nz = UBOUND(out_arr, 3)

       DO j = 0, ny
       DO i = 0, nx
       DO k = nz, 0, -1

          ! TODO: Remove IF clause and extrapolate based on a critical vertical
          ! TODO: index marking the lower bound of COSMO-DE data coverage.
          ! Check for negative interpolation weights indicating grid points 
          ! below COSMO-DE domain and extrapolate from the top in such cells.
          IF (outgrid % w_verti(i,j,k,1) < -1.0_dp .AND. k < nz)  THEN
             out_arr(i,j,k) = out_arr(i,j,k+1)
          ELSE
             out_arr(i,j,k) = 0.0_dp
             DO l = 1, 2
                out_arr(i,j,k) = out_arr(i,j,k) +                                 &
                    outgrid % w_verti(i,j,k,l) *                                  &
                    in_arr(i,j,outgrid % kk(i,j,k, l) )
             END DO
          END IF
       END DO
       END DO
       END DO
    END SUBROUTINE interpolate_1d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolates bi-linearly in horizontal planes on every k level of the output
!> array outvar(:,:,k) using values of the source array invar(:,:,:). The source
!> (invar) and interpolation array (outvar) need to have matching dimensions.
!>
!> Input parameters:
!> -----------------
!> invar : Array of source data
!> 
!> outgrid % ii, % jj : Array of neighbour indices in x and y direction.
!>     ii(i,j,k,:), and jj(i,j,k,:) contain the four horizontal neighbour points
!>     of PALM-4U point (i,j,k) on the input grid corresponding to the source
!>     data invar. (The outgrid carries the relationship with the ingrid in the
!      form of the interpoaltion weights.)
!> 
!> outgrid % w_horiz: Array of weights for horizontal bi-linear interpolation
!>     corresponding to neighbour points indexed by ii and jj.
!>
!> Output papameters:
!> ------------------
!> outvar : Array of interpolated data
!------------------------------------------------------------------------------!
    SUBROUTINE interpolate_2d(invar, outvar, outgrid, ncvar)
    ! I index 0-based for the indices of the outvar to be consistent with the
    ! outgrid indices and interpolation weights.
       TYPE(grid_definition), INTENT(IN) ::  outgrid
       REAL(dp), INTENT(IN)              ::  invar(0:,0:,0:)
       REAL(dp), INTENT(OUT)             ::  outvar(0:,0:,0:)
       TYPE(nc_var), INTENT(IN), OPTIONAL ::  ncvar

       INTEGER ::  i, j, k, l

       ! TODO: check if input dimensions are consistent, i.e. ranges are correct
       IF (UBOUND(outvar, 3) .GT. UBOUND(invar, 3))  THEN
           message = "Output array for '" // TRIM(ncvar % name) // "' has ' more levels (" // &
              TRIM(str(UBOUND(outvar, 3))) // ") than input variable ("//&
              TRIM(str(UBOUND(invar, 3))) // ")."
           CALL abort('interpolate_2d', message) 
       END IF

       DO k = 0, UBOUND(outvar, 3)
       DO j = 0, UBOUND(outvar, 2)
       DO i = 0, UBOUND(outvar, 1)
          outvar(i,j,k) = 0.0_dp
          DO l = 1, 4
             
             outvar(i,j,k) = outvar(i,j,k) +                                   &
                outgrid % w_horiz(i,j,l) * invar( outgrid % ii(i,j,l),         & 
                                                  outgrid % jj(i,j,l),         &
                                                  k )
          END DO
       END DO
       END DO
       END DO
        
    END SUBROUTINE interpolate_2d


    SUBROUTINE average_2d(in_arr, out_arr, imin, imax, jmin, jmax)
       REAL(dp), INTENT(IN)  ::  in_arr(0:,0:,0:)
       REAL(dp), INTENT(OUT) ::  out_arr(0:)
       INTEGER, INTENT(IN)   ::  imin, imax, jmin, jmax

       INTEGER  ::  i, j, k
       REAL(dp) ::  ni
       
       IF (imin < 0)  CALL abort('average_2d', "imin < 0.")
       IF (jmin < 0)  CALL abort('average_2d', "jmin < 0.")
       IF (imax > UBOUND(in_arr, 1))  CALL abort('average_2d', "imax out of i bound.")
       IF (jmax > UBOUND(in_arr, 2))  CALL abort('average_2d', "jmax out of j bound.")

       DO k = 0, UBOUND(out_arr, 1)

          out_arr(k) = 0.0_dp
          DO j = jmin, jmax
          DO i = imin, imax
             out_arr(k) = out_arr(k) + in_arr(i, j, k)
          END DO
          END DO

       END DO
   
       ! devide by number of grid points
       ni = 1.0_dp / ( (imax - imin + 1) * (jmax - jmin + 1) )
       out_arr(:) = out_arr(:) * ni

    END SUBROUTINE average_2d


    SUBROUTINE interpolate_3d(source_array, palm_array, palm_intermediate, palm_grid)
       TYPE(grid_definition), INTENT(IN) ::  palm_intermediate, palm_grid
       REAL(dp), DIMENSION(:,:,:), INTENT(IN)  ::  source_array
       REAL(dp), DIMENSION(:,:,:), INTENT(OUT) ::  palm_array
       REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::  intermediate_array
       INTEGER ::  nx, ny, nz

       nx = palm_intermediate % nx
       ny = palm_intermediate % ny
       nz = palm_intermediate % nz ! nlev

       ! Interpolate from COSMO-DE to intermediate grid. Allocating with one
       ! less point in the vertical, since scalars like T have 50 instead of 51
       ! points in COSMO-DE.
       ALLOCATE(intermediate_array(0:nx, 0:ny, 0:nz-1)) !

       CALL interpolate_2d(source_array, intermediate_array, palm_intermediate)

       ! Interpolate from intermediate grid to palm_grid grid, includes
       ! extrapolation for cells below COSMO-DE domain.
       CALL interpolate_1d(intermediate_array, palm_array, palm_grid)

       DEALLOCATE(intermediate_array)

    END SUBROUTINE interpolate_3d


    SUBROUTINE average_profile(source_array, profile_array, imin, imax, jmin, jmax,&
                               palm_intermediate, palm_grid)
       TYPE(grid_definition), INTENT(IN)          ::  palm_intermediate, palm_grid
       REAL(dp), DIMENSION(:,:,:), INTENT(IN)     ::  source_array
       INTEGER, INTENT(IN)                        ::  imin, imax, jmin, jmax
       REAL(dp), DIMENSION(0:,0:,0:), INTENT(OUT) ::  profile_array
       REAL(dp), DIMENSION(:,:,:), ALLOCATABLE    ::  intermediate_array
       INTEGER                                    ::  nx, ny, nz

       nx = palm_intermediate % nx
       ny = palm_intermediate % ny
       nz = palm_intermediate % nz
       ALLOCATE(intermediate_array(0:nx, 0:ny, 0:nz-1))
       intermediate_array(:,:,:) = 0.0_dp

       ! average input array to intermediate profile
       CALL average_2d(source_array, intermediate_array(0,0,:), imin, imax, jmin, jmax)

       ! vertically interpolate to ouput array
       CALL interpolate_1d(intermediate_array, profile_array, palm_grid)

       DEALLOCATE(intermediate_array)

    END SUBROUTINE average_profile



!-----------------------------------------------------------------------------!
! Description:
! -----------
!> This routine computes the inverse Plate Carree projection, i.e. in projects
!> Cartesian coordinates (x,y) onto a sphere. It returns the latitude and
!> lngitude of a geographical system centered at x0 and y0.
!-----------------------------------------------------------------------------!
    SUBROUTINE inv_plate_carree(x, y, x0, y0, r, lat, lon)
       REAL(dp), INTENT(IN)  ::  x(:), y(:), x0, y0, r
       REAL(dp), INTENT(OUT) ::  lat(:), lon(:)
       
       REAL(dp) :: ri

       ! TODO check dimensions of lat/lon and y/x match

       ri = 1.0_dp / r
       
       lat(:) = (y(:) - y0) * ri
       lon(:) = (x(:) - x0) * ri
    END SUBROUTINE 


!-----------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the reverse Plate-Carree projection of a x or y position on a
!> Cartesian grid.
!>
!> Input parameters:
!> -----------------
!> xy : x or y coordinate of the Cartasian grid point [m].
!>
!> xy0 : x or y coordinate that coincides with the origin of the underlying
!>     sperical system (crossing point of the equator and prime meridian) [m].
!>
!> r : Radius of the of the underlying sphere, e.g. EARTH_RADIUS [m].
!> 
!> Returns:
!> --------
!> project : Longitude (in case xy = x) or latitude (xy = y) of the given input
!>     coordinate xy.
!------------------------------------------------------------------------------!
    ELEMENTAL REAL(dp) FUNCTION project(xy, xy0, r)
       REAL(dp), INTENT(IN)  ::  xy, xy0, r
       REAL(dp) :: ri

       ! If this elemental function is called with a large array as xy, it is
       ! computationally more efficient to precompute the inverse radius and
       ! then muliply.
       ri = 1.0_dp / r

       project = (xy - xy0) * ri

    END FUNCTION project


    REAL(dp) FUNCTION phic_to_phin(phi_c)
        REAL(dp), INTENT(IN) ::  phi_c

        phic_to_phin = 0.5_dp * PI - ABS(phi_c)

    END FUNCTION phic_to_phin

    
    REAL(dp) FUNCTION lamc_to_lamn(phi_c, lam_c)
       REAL(dp), INTENT(IN) ::  phi_c, lam_c
        
       lamc_to_lamn = lam_c
       IF (phi_c > 0.0_dp)  THEN
          lamc_to_lamn = lam_c - SIGN(PI, lam_c)
       END IF

    END FUNCTION lamc_to_lamn


    REAL(dp) FUNCTION gamma_from_hemisphere(phi_cg, phi_ref)
       REAL(dp), INTENT(IN) ::  phi_cg, phi_ref
       LOGICAL ::  palm_centre_is_south_of_cosmo_origin
       
       palm_centre_is_south_of_cosmo_origin = (phi_cg < phi_ref)

       IF (palm_centre_is_south_of_cosmo_origin)  THEN
           gamma_from_hemisphere = PI
       ELSE
           gamma_from_hemisphere = 0.0_dp
       END IF
    END FUNCTION gamma_from_hemisphere


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the geographical coordinates corresponding to the given rotated-pole
!> coordinates.
!>
!> In INIFOR, this routine is used to convert coordinates between two
!> rotated-pole systems: COSMO-DE's rotated-pole system, and one centred at the
!> PALM-4U domain centre. In this case, the PALM-4U system is thought of as the
!> rotated-pole system and the routine is used to rotate back to COSMO-DE's
!> system which is thought of as the geographical one.
!> 
!> Input parameters:
!> -----------------
!> phir(:), lamr(: ): latitudes and longitudes of the rotated-pole grid
!> 
!> phip, lamp: latitude and longitude of the rotated north pole
!>
!> gam: "angle between the north poles. If [gam] is not present, the other
!>       system is the real geographical system." (original phiro2rot
!>       description)
!> 
!> Output parameters:
!> ------------------
!> phi(:,:), lam(:,:): geographical latitudes and logitudes
!------------------------------------------------------------------------------!
    SUBROUTINE rotate_to_cosmo(phir, lamr, phip, lamp, phi, lam, gam)
       REAL(dp), INTENT(IN)  ::  phir(0:), lamr(0:), phip, lamp, gam
       REAL(dp), INTENT(OUT) ::  phi(0:,0:), lam(0:,0:)

       INTEGER ::  i, j
       
       IF ( SIZE(phi, 1) .NE. SIZE(lam, 1) .OR. &
            SIZE(phi, 2) .NE. SIZE(lam, 2) )  THEN
          PRINT *, "inifor: rotate_to_cosmo: Dimensions of phi and lambda do not match. Dimensions are:"
          PRINT *, "inifor: rotate_to_cosmo: phi: ", SIZE(phi, 1), SIZE(phi, 2)
          PRINT *, "inifor: rotate_to_cosmo: lam: ", SIZE(lam, 1), SIZE(lam, 2)
          STOP
       END IF

       IF ( SIZE(phir) .NE. SIZE(phi, 2) .OR. &
            SIZE(lamr) .NE. SIZE(phi, 1) )  THEN
          PRINT *, "inifor: rotate_to_cosmo: Dimensions of phir and lamr do not match. Dimensions are:"
          PRINT *, "inifor: rotate_to_cosmo: phir: ", SIZE(phir), SIZE(phi, 2)
          PRINT *, "inifor: rotate_to_cosmo: lamr: ", SIZE(lamr), SIZE(phi, 1)
          STOP
       END IF
       
       DO j = 0, UBOUND(phir, 1)
          DO i = 0, UBOUND(lamr, 1)

             phi(i,j) = phirot2phi(phir(j) * TO_DEGREES,                       &
                                   lamr(i) * TO_DEGREES,                       &
                                   phip * TO_DEGREES,                          &
                                   lamp * TO_DEGREES,                          &
                                   gam  * TO_DEGREES) * TO_RADIANS

             lam(i,j) = rlarot2rla(phir(j) * TO_DEGREES,                       &
                                   lamr(i) * TO_DEGREES,                       &
                                   phip * TO_DEGREES,                          &
                                   lamp * TO_DEGREES,                          &
                                   gam  * TO_DEGREES) * TO_RADIANS

          END DO
       END DO

    END SUBROUTINE rotate_to_cosmo


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute indices of PALM-4U grid point neighbours in the target 
!> system (COSMO-DE) by rounding up and down. (i,j) are the indices of
!> the PALM-4U grid and (ii(i,j,1-4), jj(i,j,1-4)) contain the indices
!> of the its four neigbouring points in the COSMO-DE grid.
!>
!>
!>                     COSMO-DE grid
!>                     -------------
!>           jj, lat
!>              ^      j
!>              |       \          i
!>  jj(i,j,2/3) + ... 2 ---\--------/------ 3
!>              |     | ^   \      /        |
!>              |     | |wp  \    /         |
!>              |     | v     \  /          |
!>       latpos + ............ o/ (i,j)     |
!>              |     |        :            |
!>              |     |        :<----wl---->|
!>  jj(i,j,1/4) + ... 1 -------:----------- 4
!>              |     :        :            :
!>              |     :        :            :
!>              |     :      lonpos         :
!>              L-----+--------+------------+------> ii, lon 
!>               ii(i,j,1/2)        ii(i,j,3/4)
!>
!> 
!> Input parameters:
!> -----------------
!> source_lat, source_lon : (rotated-pole) coordinates of the source grid (e.g.
!>    COSMO-DE)
!>
!> source_dxi, source_dyi : inverse grid spacings of the source grid.
!>
!> target_lat, target_lon : (rotated-pole) coordinates of the target grid (e.g.
!>    COSMO-DE)
!> 
!> Output parameters:
!> ------------------
!> palm_ii, palm_jj : x and y index arrays of horizontal neighbour columns
!> 
!------------------------------------------------------------------------------!
    SUBROUTINE find_horizontal_neighbours(cosmo_lat, cosmo_lon, cosmo_dxi,     &
       cosmo_dyi, palm_clat, palm_clon, palm_ii, palm_jj)

       REAL(dp), DIMENSION(0:), INTENT(IN)        ::  cosmo_lat, cosmo_lon
       REAL(dp), DIMENSION(0:,0:), INTENT(IN)     ::  palm_clat, palm_clon
       REAL(dp), INTENT(IN)                       ::  cosmo_dxi, cosmo_dyi
       INTEGER, DIMENSION(0:,0:,1:), INTENT(OUT)  ::  palm_ii, palm_jj

       REAL(dp) ::  lonpos, latpos, lon0, lat0
       INTEGER  ::  i, j

       lon0 = cosmo_lon(0)
       lat0 = cosmo_lat(0)

       DO j = 0, UBOUND(palm_clon, 2)!palm_grid % ny
       DO i = 0, UBOUND(palm_clon, 1)!palm_grid % nx
          ! Compute the floating point index corrseponding to PALM-4U grid point
          ! location along target grid (COSMO-DE) axes.
          lonpos = (palm_clon(i,j) - lon0) * cosmo_dxi
          latpos = (palm_clat(i,j) - lat0) * cosmo_dyi

          IF (lonpos < 0.0 .OR. latpos < 0.0)  THEN
             PRINT *, " Error while finding neighbours: lonpos or latpos out of bounds!"
             PRINT *, "     (i,j) = (", i, ",",j,")"
             PRINT *, "      lonpos ", lonpos*TO_DEGREES, ", latpos ", latpos*TO_DEGREES
             PRINT *, "        lon0 ", lon0  *TO_DEGREES,   ", lat0   ", lat0*TO_DEGREES
             PRINT *, "    PALM lon ", palm_clon(i,j)*TO_DEGREES,   ", PALM lat ",palm_clat(i,j)*TO_DEGREES
             STOP
          END IF

          palm_ii(i,j,1) = FLOOR(lonpos)
          palm_ii(i,j,2) = FLOOR(lonpos)
          palm_ii(i,j,3) = CEILING(lonpos)
          palm_ii(i,j,4) = CEILING(lonpos)

          palm_jj(i,j,1) = FLOOR(latpos)
          palm_jj(i,j,2) = CEILING(latpos)
          palm_jj(i,j,3) = CEILING(latpos)
          palm_jj(i,j,4) = FLOOR(latpos)
       END DO
       END DO

    END SUBROUTINE find_horizontal_neighbours

    
    SUBROUTINE find_vertical_neighbours_and_weights(palm_grid, palm_intermediate)
       TYPE(grid_definition), INTENT(INOUT) ::  palm_grid
       TYPE(grid_definition), INTENT(IN)    ::  palm_intermediate

       INTEGER  ::  i, j, k, nx, ny, nz, nlev, kcur
       LOGICAL  ::  point_is_below_grid, point_is_above_grid,                  &
                    point_is_in_current_cell
       REAL(dp) ::  current_height, column_base, column_top, h_top, h_bottom,  &
                    weight

       nx   = palm_grid % nx
       ny   = palm_grid % ny
       nz   = palm_grid % nz
       nlev = palm_intermediate % nz

       ! in each column of the fine grid, find vertical neighbours of every cell
       DO i = 0, nx
       DO j = 0, ny

          kcur = 0

          column_base = palm_intermediate % h(i,j,0)
          column_top  = palm_intermediate % h(i,j,nlev)

          ! scan through palm_grid column until and set neighbour indices in
          ! case current_height is either below column_base, in the current
          ! cell, or above column_top. Keep increasing current cell index until
          ! the current cell overlaps with the current_height.
          DO k = 0, nz

             ! Memorize the top and bottom boundaries of the coarse cell and the
             ! current height within it
             current_height = palm_grid % z(k) + palm_grid % z0
             h_top    = palm_intermediate % h(i,j,kcur+1)
             h_bottom = palm_intermediate % h(i,j,kcur)

             point_is_above_grid = (current_height > column_top) !22000m, very unlikely
             point_is_below_grid = (current_height < column_base)

             point_is_in_current_cell = (                                      &
                current_height >= h_bottom .AND.                               &
                current_height <  h_top                                        &
             )

             ! set default weights
             palm_grid % w_verti(i,j,k,1:2) = 0.0_dp

             IF (point_is_above_grid)  THEN

                palm_grid % kk(i,j,k,1:2) = nlev
                palm_grid % w_verti(i,j,k,1:2) = - 2.0_dp

             ELSE IF (point_is_below_grid)  THEN

                palm_grid % kk(i,j,k,1:2) = 0
                palm_grid % w_verti(i,j,k,1:2) = - 2.0_dp

             ELSE
                ! cycle through intermediate levels until current
                ! intermediate-grid cell overlaps with current_height
                DO WHILE (.NOT. point_is_in_current_cell .AND. kcur <= nlev-1)
                   kcur = kcur + 1

                   h_top    = palm_intermediate % h(i,j,kcur+1)
                   h_bottom = palm_intermediate % h(i,j,kcur)
                   point_is_in_current_cell = (                                &
                      current_height >= h_bottom .AND.                         &
                      current_height <  h_top                                  &
                   )
                END DO

                ! kcur = 48 indicates the last section (indices 48 and 49), i.e.
                ! kcur = 49 is not the beginning of a valid cell.
                IF (kcur > nlev-1)  THEN
                   message = "Index " // TRIM(str(kcur)) // " is above intermediate grid range."
                   CALL abort('find_vertical_neighbours', message)
                END IF
   
                palm_grid % kk(i,j,k,1) = kcur
                palm_grid % kk(i,j,k,2) = kcur + 1

                ! copmute vertical weights
                weight = (h_top - current_height) / (h_top - h_bottom)
                palm_grid % w_verti(i,j,k,1) = weight
                palm_grid % w_verti(i,j,k,2) = 1.0_dp - weight
             END IF

          END DO

       END DO
       END DO

    END SUBROUTINE find_vertical_neighbours_and_weights

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the four weights for horizontal bilinear interpolation given the
!> coordinates clon(i,j) clat(i,j) of the PALM-4U grid in the COSMO-DE
!> rotated-pole grid and the neightbour indices ii(i,j,1-4) and jj(i,j,1-4).
!>
!> Input parameters:
!> -----------------
!> palm_grid % clon : longitudes of PALM-4U scalars (cell centres) in COSMO-DE's rotated-pole grid [rad]
!>
!> palm_grid % clat : latitudes of PALM-4U cell centres in COSMO-DE's rotated-pole grid [rad]
!>
!> cosmo_grid % lon : rotated-pole longitudes of scalars (cell centres) of the COSMO-DE grid [rad] 
!>
!> cosmo_grid % lat : rotated-pole latitudes of scalars (cell centers) of the COSMO-DE grid [rad]
!>
!> cosmo_grid % dxi : inverse grid spacing in the first dimension [m^-1]
!>
!> cosmo_grid % dyi : inverse grid spacing in the second dimension [m^-1]
!>
!> Output parameters:
!> ------------------
!> palm_grid % w_horiz(:,:,1-4) : weights for bilinear horizontal interpolation
!
!                               COSMO-DE grid
!                               -------------
!                     jj, lat
!                        ^        j
!                        |         \          i
!            jj(i,j,2/3) + ... 2 ---\--------/------ 3
!                        |     | ^   \      /        |
!                        |     | |wp  \    /         |
!                        |     | v     \  /          |
!                 latpos + ............ o/ (i,j)     |
!                        |     |        :            |
!                        |     |        :<----wl---->|
!            jj(i,j,1/4) + ... 1 -------:----------- 4
!                        |     :        :            :
!                        |     :        :            :
!                        |     :      lonpos         :
!                        L-----+--------+------------+------> ii, lon 
!                         ii(i,j,1/2)        ii(i,j,3/4)
!          
    SUBROUTINE compute_horizontal_interp_weights(cosmo_lat, cosmo_lon,         &
       cosmo_dxi, cosmo_dyi, palm_clat, palm_clon, palm_ii, palm_jj, palm_w_horiz)
       
       REAL(dp), DIMENSION(0:), INTENT(IN)        ::  cosmo_lat, cosmo_lon
       REAL(dp), INTENT(IN)                       ::  cosmo_dxi, cosmo_dyi
       REAL(dp), DIMENSION(0:,0:), INTENT(IN)     ::  palm_clat, palm_clon
       INTEGER, DIMENSION(0:,0:,1:), INTENT(IN)   ::  palm_ii, palm_jj

       REAL(dp), DIMENSION(0:,0:,1:), INTENT(OUT) ::  palm_w_horiz

       REAL(dp) ::  wl, wp
       INTEGER  ::  i, j

       DO j = 0, UBOUND(palm_clon, 2)
       DO i = 0, UBOUND(palm_clon, 1)
      
          ! weight in lambda direction
          wl = ( cosmo_lon(palm_ii(i,j,4)) - palm_clon(i,j) ) * cosmo_dxi

          ! weight in phi direction
          wp = ( cosmo_lat(palm_jj(i,j,2)) - palm_clat(i,j) ) * cosmo_dyi

          IF (wl > 1.0_dp .OR. wl < 0.0_dp)  THEN
              message = "Horizontal weight wl = " // TRIM(real_to_str(wl)) //   &
                        " is out bounds."
              CALL abort('compute_horizontal_interp_weights', message)
          END IF
          IF (wp > 1.0_dp .OR. wp < 0.0_dp)  THEN
              message = "Horizontal weight wp = " // TRIM(real_to_str(wp)) //   &
                        " is out bounds."
              CALL abort('compute_horizontal_interp_weights', message)
          END IF

          palm_w_horiz(i,j,1) = wl * wp
          palm_w_horiz(i,j,2) = wl * (1.0_dp - wp)
          palm_w_horiz(i,j,3) = (1.0_dp - wl) * (1.0_dp - wp)
          palm_w_horiz(i,j,4) = 1.0_dp - SUM( palm_w_horiz(i,j,1:3) )

       END DO
       END DO
       
    END SUBROUTINE compute_horizontal_interp_weights


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolates u and v components of velocities located at cell faces to the
!> cell centres by averaging neighbouring values.
!>
!> This routine is designed to be used with COSMO-DE arrays where there are the
!> same number of grid points for scalars (centres) and velocities (faces). In
!> COSMO-DE the velocity points are staggared one half grid spaceing up-grid
!> which means the first centre point has to be omitted and is set to zero.
    SUBROUTINE centre_velocities(u_face, v_face, u_centre, v_centre)
       REAL(dp), DIMENSION(0:,0:,0:), INTENT(IN)  ::  u_face, v_face
       REAL(dp), DIMENSION(0:,0:,0:), INTENT(OUT) ::  u_centre, v_centre
       INTEGER ::  nx, ny

       nx = UBOUND(u_face, 1)
       ny = UBOUND(u_face, 2)

       u_centre(0,:,:)  = 0.0_dp
       u_centre(1:,:,:) = 0.5_dp * ( u_face(0:nx-1,:,:) + u_face(1:,:,:) )

       v_centre(:,0,:)  = 0.0_dp
       v_centre(:,1:,:) = 0.5_dp * ( v_face(:,0:ny-1,:) + v_face(:,1:,:) )
    END SUBROUTINE centre_velocities


    FUNCTION phirot2phi (phirot, rlarot, polphi, pollam, polgam)
    
       REAL(dp), INTENT (IN) ::  polphi      !< latitude of the rotated north pole
       REAL(dp), INTENT (IN) ::  pollam      !< longitude of the rotated north pole
       REAL(dp), INTENT (IN) ::  phirot      !< latitude in the rotated system
       REAL(dp), INTENT (IN) ::  rlarot      !< longitude in the rotated system
       REAL(dp), INTENT (IN) ::  polgam      !< angle between the north poles of the systems

       REAL(dp)              ::  phirot2phi  !< latitude in the geographical system
       
       REAL(dp)              ::  zsinpol, zcospol, zphis, zrlas, zarg, zgam
    
       zsinpol = SIN(polphi * TO_RADIANS)
       zcospol = COS(polphi * TO_RADIANS)
       zphis   = phirot * TO_RADIANS

       IF (rlarot > 180.0_dp)  THEN
          zrlas = rlarot - 360.0_dp
       ELSE
          zrlas = rlarot
       END IF
       zrlas = zrlas * TO_RADIANS
     
       IF (polgam /= 0.0_dp)  THEN
          zgam = polgam * TO_RADIANS
          zarg = zsinpol * SIN (zphis) +                                       &
                 zcospol * COS(zphis) * ( COS(zrlas) * COS(zgam) -             &
                                          SIN(zgam)  * SIN(zrlas) )
       ELSE
          zarg = zcospol * COS (zphis) * COS (zrlas) + zsinpol * SIN (zphis)
       END IF
      
       phirot2phi = ASIN (zarg) * TO_DEGREES
    
    END FUNCTION phirot2phi


    FUNCTION phi2phirot (phi, rla, polphi, pollam)
    
       REAL(dp), INTENT (IN) ::  polphi !< latitude of the rotated north pole
       REAL(dp), INTENT (IN) ::  pollam !< longitude of the rotated north pole
       REAL(dp), INTENT (IN) ::  phi    !< latitude in the geographical system
       REAL(dp), INTENT (IN) ::  rla    !< longitude in the geographical system
       
       REAL(dp) ::  phi2phirot          !< longitude in the rotated system
       
       REAL(dp) ::  zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1
       
       zsinpol = SIN(polphi * TO_RADIANS)
       zcospol = COS(polphi * TO_RADIANS)
       zlampol = pollam * TO_RADIANS
       zphi    = phi * TO_RADIANS

       IF (rla > 180.0_dp)  THEN
          zrla1 = rla - 360.0_dp
       ELSE
          zrla1 = rla
       END IF
       zrla = zrla1 * TO_RADIANS
       
       zarg1 = SIN(zphi) * zsinpol
       zarg2 = COS(zphi) * zcospol * COS(zrla - zlampol)
       
       phi2phirot = ASIN(zarg1 + zarg2) * TO_DEGREES
    
    END FUNCTION phi2phirot


    FUNCTION rlarot2rla(phirot, rlarot, polphi, pollam, polgam)
    
       REAL(dp), INTENT (IN) ::  polphi !< latitude of the rotated north pole
       REAL(dp), INTENT (IN) ::  pollam !< longitude of the rotated north pole
       REAL(dp), INTENT (IN) ::  phirot !< latitude in the rotated system
       REAL(dp), INTENT (IN) ::  rlarot !< longitude in the rotated system
       REAL(dp), INTENT (IN) ::  polgam !< angle between the north poles of the systems
       
       REAL(dp) ::  rlarot2rla          !< latitude in the geographical system
       
       REAL(dp) ::  zsinpol, zcospol, zlampol, zphis, zrlas, zarg1, zarg2, zgam
       
       zsinpol = SIN(TO_RADIANS * polphi)
       zcospol = COS(TO_RADIANS * polphi)
       zlampol = TO_RADIANS * pollam
       zphis   = TO_RADIANS * phirot

       IF (rlarot > 180.0_dp)  THEN
          zrlas = rlarot - 360.0_dp
       ELSE
          zrlas = rlarot
       END IF
       zrlas   = TO_RADIANS * zrlas
      
       IF (polgam /= 0.0_dp)  THEN
          zgam  = TO_RADIANS * polgam
          zarg1 = SIN(zlampol) * (zcospol * SIN(zphis) - zsinpol*COS(zphis) *  &
                  (COS(zrlas) * COS(zgam) - SIN(zrlas) * SIN(zgam)) ) -        &
                  COS(zlampol) * COS(zphis) * ( SIN(zrlas) * COS(zgam) +       &
                                                COS(zrlas) * SIN(zgam) )
       
          zarg2 = COS (zlampol) * (zcospol * SIN(zphis) - zsinpol*COS(zphis) * &
                  (COS(zrlas) * COS(zgam) - SIN(zrlas) * SIN(zgam)) ) +        &
                  SIN(zlampol) * COS(zphis) * ( SIN(zrlas) * COS(zgam) +       &
                                                COS(zrlas) * SIN(zgam) )
       ELSE
          zarg1   = SIN (zlampol) * (-zsinpol * COS(zrlas) * COS(zphis)  +     &
                                      zcospol *              SIN(zphis)) -     &
                    COS (zlampol) *             SIN(zrlas) * COS(zphis)
          zarg2   = COS (zlampol) * (-zsinpol * COS(zrlas) * COS(zphis)  +     &
                                      zcospol *              SIN(zphis)) +     &
                    SIN (zlampol) *             SIN(zrlas) * COS(zphis)
       END IF
      
       IF (zarg2 == 0.0_dp)  zarg2 = 1.0E-20_dp
      
       rlarot2rla = ATAN2(zarg1,zarg2) * TO_DEGREES
        
    END FUNCTION rlarot2rla


    FUNCTION rla2rlarot ( phi, rla, polphi, pollam, polgam )

       REAL(dp), INTENT (IN) ::  polphi !< latitude of the rotated north pole
       REAL(dp), INTENT (IN) ::  pollam !< longitude of the rotated north pole
       REAL(dp), INTENT (IN) ::  phi    !< latitude in geographical system
       REAL(dp), INTENT (IN) ::  rla    !< longitude in geographical system
       REAL(dp), INTENT (IN) ::  polgam !< angle between the north poles of the systems
       
       REAL (KIND=dp) ::  rla2rlarot    !< latitude in the the rotated system
       
       REAL (KIND=dp) ::  zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1
       
       zsinpol = SIN(polphi * TO_RADIANS)
       zcospol = COS(polphi * TO_RADIANS)
       zlampol = pollam * TO_RADIANS
       zphi    = phi * TO_RADIANS

       IF (rla > 180.0_dp)  THEN
          zrla1 = rla - 360.0_dp
       ELSE
          zrla1 = rla
       END IF
       zrla = zrla1 * TO_RADIANS
       
       zarg1 = - SIN (zrla-zlampol) * COS(zphi)
       zarg2 = - zsinpol * COS(zphi) * COS(zrla-zlampol) + zcospol * SIN(zphi)
       
       IF (zarg2 == 0.0_dp)  zarg2 = 1.0E-20_dp
       
       rla2rlarot = ATAN2 (zarg1,zarg2) * TO_DEGREES
       
       IF (polgam /= 0.0_dp )  THEN
          rla2rlarot = polgam + rla2rlarot
          IF (rla2rlarot > 180._dp)  rla2rlarot = rla2rlarot - 360.0_dp
       END IF
       
    END FUNCTION rla2rlarot


    SUBROUTINE uv2uvrot(u, v, rlat, rlon, pollat, pollon, urot, vrot)
    
       REAL(dp), INTENT (IN)  ::  u, v           !< wind components in the true geographical system
       REAL(dp), INTENT (IN)  ::  rlat, rlon     !< coordinates in the true geographical system
       REAL(dp), INTENT (IN)  ::  pollat, pollon !< latitude and longitude of the north pole of the rotated grid
       
       REAL(dp), INTENT (OUT) ::  urot, vrot     !< wind components in the rotated grid             
       
       REAL (dp) ::  zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm
       
       zsinpol = SIN(pollat * TO_RADIANS)
       zcospol = COS(pollat * TO_RADIANS)
       zlonp   = (pollon-rlon) * TO_RADIANS
       zlat    = rlat * TO_RADIANS
       
       zarg1 = zcospol * SIN(zlonp)
       zarg2 = zsinpol * COS(zlat) - zcospol * SIN(zlat) * COS(zlonp)
       znorm = 1.0_dp / SQRT(zarg1*zarg1 + zarg2*zarg2)
       
       urot = u * zarg2 * znorm - v * zarg1 * znorm
       vrot = u * zarg1 * znorm + v * zarg2 * znorm
    
    END SUBROUTINE uv2uvrot


    SUBROUTINE uvrot2uv (urot, vrot, rlat, rlon, pollat, pollon, u, v)
    
       REAL(dp), INTENT(IN) ::  urot, vrot     !< wind components in the rotated grid
       REAL(dp), INTENT(IN) ::  rlat, rlon     !< latitude and longitude in the true geographical system
       REAL(dp), INTENT(IN) ::  pollat, pollon !< latitude and longitude of the north pole of the rotated grid
       
       REAL(dp), INTENT(OUT) ::  u, v          !< wind components in the true geographical system
       
       REAL(dp) ::  zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm
     
       zsinpol = SIN(pollat * TO_RADIANS)
       zcospol = COS(pollat * TO_RADIANS)
       zlonp   = (pollon-rlon) * TO_RADIANS
       zlat    = rlat * TO_RADIANS
     
       zarg1 = zcospol * SIN(zlonp)
       zarg2 = zsinpol * COS(zlat) - zcospol * SIN(zlat) * COS(zlonp)
       znorm = 1.0_dp / SQRT(zarg1*zarg1 + zarg2*zarg2)
     
       u =   urot * zarg2 * znorm + vrot * zarg1 * znorm
       v = - urot * zarg1 * znorm + vrot * zarg2 * znorm
    
    END SUBROUTINE uvrot2uv

 END MODULE

