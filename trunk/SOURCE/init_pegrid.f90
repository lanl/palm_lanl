!> @file init_pegrid.f90
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
! $Id: init_pegrid.f90 3058 2018-06-05 09:21:14Z raasch $
! bugfix: wrong error number in r3057 revised
! 
! 3057 2018-06-05 09:03:41Z raasch
! bugfix: check that nz is even in case that optimized multigrid is used
! 
! 3049 2018-05-29 13:52:36Z Giersch
! Error messages revised
! 
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
! 
! 2938 2018-03-27 15:52:42Z suehring
! - No checks for domain decomposition in case of turbulence generator 
!  (is done in stg module)
! - Introduce ids to indicate lateral processors for turbulence generator
! 
! 2936 2018-03-27 14:49:27Z suehring
! Variable use_synthetic_turbulence_generator has been abbreviated
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! 3D-Integer exchange on multigrid level (MS)
! Forcing implemented (MS)
! 
! 2600 2017-11-01 14:11:20Z raasch
! calculation of block-I/O quantitites removed (is now done in parin)
! 
! 2516 2017-10-04 11:03:04Z suehring
! Remove tabs
! 
! 2514 2017-10-04 09:52:37Z suehring
! Redundant preprocessor directives removed 
! 
! 2372 2017-08-25 12:37:32Z sward
! Shifted cyclic boundary conditions implemented
! 
! 2365 2017-08-21 14:59:59Z kanani
! Vertical nesting implemented (SadiqHuq)
! 
! 2300 2017-06-29 13:31:14Z raasch
! host-specific settings removed
! 
! 2298 2017-06-29 09:28:18Z raasch
! MPI2 related parts removed
! 
! 2271 2017-06-09 12:34:55Z sward
! Error message changed
! 
! 2259 2017-06-08 09:09:11Z gronemeier
! Implemented synthetic turbulence generator
!
! 2238 2017-05-31 16:49:16Z suehring
! Remove unnecessary module load of pmc_interface
! 
! 2231 2017-05-30 16:44:33Z suehring
!
! 2200 2017-04-11 11:37:51Z suehring
! monotonic_adjustment removed
!
! 2197 2017-03-24 02:25:00Z raasch
! bugfix: do not allow odd values for nz at the coarsest grid level in case of
! optimized multigrid solver
!
! 2180 2017-03-17 13:33:05Z hellstea
! Checks to ensure (2178) that pdims match the grid dimensions in the
! automatic determination of pdims are canceled as unnecessary
! 
! 2178 2017-03-17 11:07:39Z hellstea
! Checks to ensure that pdims match the grid dimensions are added in the
! automatic determination of pdims
!
! 2050 2016-11-08 15:00:55Z gronemeier
! Implement turbulent outflow condition
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1968 2016-07-18 12:01:49Z suehring
! Extent MPI-datatypes for exchange of 2D-INTEGER arrays on coarser multigrid
! level  
! 
! 1964 2016-07-14 15:35:18Z hellstea
! Bugfix: erroneous setting of nest_bound_l/r/s/n = .TRUE. for vertical nesting mode removed.
!
! 1923 2016-05-31 16:37:07Z boeske
! Initial version of purely vertical nesting introduced. 
! 
! 1922 2016-05-31 16:36:08Z boeske
! Bugfix: array transposition checks restricted to cases if a fourier 
! transform is used , removed unused variable nnx_z
! 
! 1833 2016-04-07 14:23:03Z raasch
! spectra related variables moved to spectra_mod
!
! 1815 2016-04-06 13:49:59Z raasch
! cpp-directives for intel openmp bug removed
!
! 1804 2016-04-05 16:30:18Z maronga
! Removed code for parameter file check (__check)
!
! 1779 2016-03-03 08:01:28Z raasch
! changes regarding nested domain removed: virtual PE grid will be automatically
! calculated for nested runs too
!
! 1764 2016-02-28 12:45:19Z raasch
! cpp-statements for nesting removed
!
! 1762 2016-02-25 12:31:13Z hellstea
! Introduction of nested domain feature
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1677 2015-10-02 13:25:23Z boeske
! New MPI-data types for exchange of 3D integer arrays.
!
! 1575 2015-03-27 09:56:27Z raasch
! adjustments for psolver-queries, calculation of ngp_xz added
!
! 1565 2015-03-09 20:59:31Z suehring
! Refine if-clause for setting nbgp. 
!
! 1557 2015-03-05 16:43:04Z suehring
! Adjustment for monotonic limiter
!
! 1468 2014-09-24 14:06:57Z maronga
! Adapted for use on up to 6-digit processor cores
! 
! 1435 2014-07-21 10:37:02Z keck
! bugfix: added missing parameter coupling_mode_remote to ONLY-attribute
! 
! 1402 2014-05-09 14:25:13Z raasch
! location messages modified
! 
! 1384 2014-05-02 14:31:06Z raasch
! location messages added
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL functions provided with KIND-attribute
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1304 2014-03-12 10:29:42Z raasch
! bugfix: single core MPI runs missed some settings of transpose indices
!
! 1212 2013-08-15 08:46:27Z raasch
! error message for poisfft_hybrid removed
!
! 1159 2013-05-21 11:58:22Z fricke
! dirichlet/neumann and neumann/dirichlet removed
!
! 1139 2013-04-18 07:25:03Z raasch
! bugfix for calculating the id of the PE carrying the recycling plane
!
! 1111 2013-03-08 23:54:10Z raasch
! initialization of poisfft moved to module poisfft
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1056 2012-11-16 15:28:04Z raasch
! Indices for arrays n.._mg start from zero due to definition of arrays f2 and
! p2 as automatic arrays in recursive subroutine next_mg_level
!
! 1041 2012-11-06 02:36:29Z raasch
! a 2d virtual processor topology is used by default for all machines
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1003 2012-09-14 14:35:53Z raasch
! subdomains must have identical size (grid matching = "match" removed)
!
! 1001 2012-09-13 14:08:46Z raasch
! all actions concerning upstream-spline-method removed
!
! 978 2012-08-09 08:28:32Z fricke
! dirichlet/neumann and neumann/dirichlet added
! nxlu and nysv are also calculated for inflow boundary
!
! 809 2012-01-30 13:32:58Z maronga
! Bugfix: replaced .AND. and .NOT. with && and ! in the preprocessor directives
!
! 807 2012-01-25 11:53:51Z maronga
! New cpp directive "__check" implemented which is used by check_namelist_files
!
! Revision 1.1  1997/07/24 11:15:09  raasch
! Initial revision
!
!
! Description:
! ------------
!> Determination of the virtual processor topology (if not prescribed by the
!> user)and computation of the grid point number and array bounds of the local
!> domains.
!> @todo: remove MPI-data types for 2D exchange on coarse multigrid level (not
!>        used any more) 
!------------------------------------------------------------------------------!
 SUBROUTINE init_pegrid
 

    USE control_parameters,                                                    &
        ONLY:  bc_lr, bc_ns, coupling_mode, coupling_mode_remote,              &
               coupling_topology, force_bound_l, force_bound_n, force_bound_r, &
               force_bound_s, gathered_size, grid_level,                       &
               grid_level_count, inflow_l, inflow_n, inflow_r, inflow_s,       &
               maximum_grid_level, message_string,                             &
               mg_switch_to_pe0_level, momentum_advec, nest_bound_l,           &
               nest_bound_n, nest_bound_r, nest_bound_s, nest_domain, neutral, &
               psolver, outflow_l, outflow_n, outflow_r, outflow_s,            &
               outflow_source_plane, recycling_width, scalar_advec,            &
               subdomain_size, turbulent_outflow, y_shift

    USE grid_variables,                                                        &
        ONLY:  dx
        
    USE indices,                                                               &
        ONLY:  mg_loc_ind, nbgp, nnx, nny, nnz, nx, nx_a, nx_o, nxl, nxl_mg,   &
               nxlu, nxr, nxr_mg, ny, ny_a, ny_o, nyn, nyn_mg, nys, nys_mg,    &
               nysv, nz, nzb, nzt, nzt_mg, wall_flags_1, wall_flags_2,         &
               wall_flags_3, wall_flags_4, wall_flags_5, wall_flags_6,         &
               wall_flags_7, wall_flags_8, wall_flags_9, wall_flags_10

    USE kinds
      
    USE pegrid
      
    USE transpose_indices,                                                     &
        ONLY:  nxl_y, nxl_yd, nxl_z, nxr_y, nxr_yd, nxr_z, nyn_x, nyn_z, nys_x,&
               nys_z, nzb_x, nzb_y, nzb_yd, nzt_x, nzt_yd, nzt_y

    IMPLICIT NONE

    INTEGER(iwp) ::  i                        !<
    INTEGER(iwp) ::  id_inflow_l              !<
    INTEGER(iwp) ::  id_outflow_l             !< local value of id_outflow
    INTEGER(iwp) ::  id_outflow_source_l      !< local value of id_outflow_source
    INTEGER(iwp) ::  id_recycling_l           !<
    INTEGER(iwp) ::  id_stg_left_l            !< left lateral boundary local core id in case of turbulence generator  
    INTEGER(iwp) ::  id_stg_north_l           !< north lateral boundary local core id in case of turbulence generator  
    INTEGER(iwp) ::  id_stg_right_l           !< right lateral boundary local core id in case of turbulence generator  
    INTEGER(iwp) ::  id_stg_south_l           !< south lateral boundary local core id in case of turbulence generator  
    INTEGER(iwp) ::  ind(5)                   !<
    INTEGER(iwp) ::  j                        !<
    INTEGER(iwp) ::  k                        !<
    INTEGER(iwp) ::  maximum_grid_level_l     !<
    INTEGER(iwp) ::  mg_levels_x              !<
    INTEGER(iwp) ::  mg_levels_y              !<
    INTEGER(iwp) ::  mg_levels_z              !<
    INTEGER(iwp) ::  mg_switch_to_pe0_level_l !<
    INTEGER(iwp) ::  nnx_y                    !<
    INTEGER(iwp) ::  nnx_z                    !<
    INTEGER(iwp) ::  nny_x                    !<
    INTEGER(iwp) ::  nny_z                    !<
    INTEGER(iwp) ::  nnz_x                    !<
    INTEGER(iwp) ::  nnz_y                    !<
    INTEGER(iwp) ::  numproc_sqr              !<
    INTEGER(iwp) ::  nxl_l                    !<
    INTEGER(iwp) ::  nxr_l                    !<
    INTEGER(iwp) ::  nyn_l                    !<
    INTEGER(iwp) ::  nys_l                    !<
    INTEGER(iwp) ::  nzb_l                    !<
    INTEGER(iwp) ::  nzt_l                    !<
    INTEGER(iwp) ::  omp_get_num_threads      !<

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  ind_all !<
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nxlf    !<
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nxrf    !<
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nynf    !<
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nysf    !<

    INTEGER(iwp), DIMENSION(2) ::  pdims_remote         !<
    INTEGER(iwp)               ::  lcoord(2)            !< PE coordinates of left neighbor along x and y
    INTEGER(iwp)               ::  rcoord(2)            !< PE coordinates of right neighbor along x and y

!
!-- Get the number of OpenMP threads
    !$OMP PARALLEL
!$  threads_per_task = omp_get_num_threads()
    !$OMP END PARALLEL


#if defined( __parallel )

    CALL location_message( 'creating virtual PE grids + MPI derived data types', &
                           .FALSE. )

!
!-- Determine the processor topology or check it, if prescribed by the user
    IF ( npex == -1  .AND.  npey == -1 )  THEN

!
!--    Automatic determination of the topology
       numproc_sqr = SQRT( REAL( numprocs, KIND=wp ) )
       pdims(1)    = MAX( numproc_sqr , 1 )
       DO  WHILE ( MOD( numprocs , pdims(1) ) /= 0 )
          pdims(1) = pdims(1) - 1
       ENDDO
       pdims(2) = numprocs / pdims(1)

    ELSEIF ( npex /= -1  .AND.  npey /= -1 )  THEN

!
!--    Prescribed by user. Number of processors on the prescribed topology
!--    must be equal to the number of PEs available to the job
       IF ( ( npex * npey ) /= numprocs )  THEN
          WRITE( message_string, * ) 'number of PEs of the prescribed ',       &
              'topology (', npex*npey,') does not match & the number of ',     &
              'PEs available to the job (', numprocs, ')'
          CALL message( 'init_pegrid', 'PA0221', 1, 2, 0, 6, 0 )
       ENDIF
       pdims(1) = npex
       pdims(2) = npey

    ELSE
!
!--    If the processor topology is prescribed by the user, the number of
!--    PEs must be given in both directions
       message_string = 'if the processor topology is prescribed by th' //     &
                'e user & both values of "npex" and "npey" must be given' //   &
                ' in the &NAMELIST-parameter file'
       CALL message( 'init_pegrid', 'PA0222', 1, 2, 0, 6, 0 )

    ENDIF

!
!-- If necessary, set horizontal boundary conditions to non-cyclic
    IF ( bc_lr /= 'cyclic' )  cyclic(1) = .FALSE.
    IF ( bc_ns /= 'cyclic' )  cyclic(2) = .FALSE.


!
!-- Create the virtual processor grid
    CALL MPI_CART_CREATE( comm_palm, ndim, pdims, cyclic, reorder, &
                          comm2d, ierr )
    CALL MPI_COMM_RANK( comm2d, myid, ierr )
    WRITE (myid_char,'(''_'',I6.6)')  myid

    CALL MPI_CART_COORDS( comm2d, myid, ndim, pcoord, ierr )
    CALL MPI_CART_SHIFT( comm2d, 0, 1, pleft, pright, ierr )
    CALL MPI_CART_SHIFT( comm2d, 1, 1, psouth, pnorth, ierr )
!
!-- In case of cyclic boundary conditions, a y-shift at the boundaries in
!-- x-direction can be introduced via parameter y_shift. The shift is done
!-- by modifying the processor grid in such a way that processors located
!-- at the x-boundary communicate across it to processors with y-coordinate
!-- shifted by y_shift relative to their own. This feature can not be used
!-- in combination with an fft pressure solver. It has been implemented to
!-- counter the effect of streak structures in case of cyclic boundary
!-- conditions. For a description of these see Munters
!-- (2016; dx.doi.org/10.1063/1.4941912)
!--
!-- Get coordinates of left and right neighbor on PE grid
    IF ( y_shift /= 0 ) THEN

       IF ( bc_lr /= 'cyclic'  .OR.  bc_ns /= 'cyclic' )  THEN
          message_string = 'y_shift /= 0 is only allowed for cyclic ' //       &
                           'boundary conditions in both directions '
          CALL message( 'check_parameters', 'PA0467', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( TRIM( psolver ) /= 'multigrid' .AND.                               &
            TRIM( psolver ) /= 'multigrid_noopt')                              &
       THEN
          message_string = 'y_shift /= 0 requires a multigrid pressure solver '
          CALL message( 'check_parameters', 'PA0468', 1, 2, 0, 6, 0 )
       ENDIF

       CALL MPI_CART_COORDS( comm2d, pright, ndim, rcoord, ierr )
       CALL MPI_CART_COORDS( comm2d, pleft, ndim, lcoord, ierr )

!
!--    If the x(y)-coordinate of the right (left) neighbor is smaller (greater)
!--    than that of the calling process, then the calling process is located on
!--    the right (left) boundary of the processor grid. In that case,
!--    the y-coordinate of that neighbor is increased (decreased) by y_shift.
!--    The rank of the process with that coordinate is then inquired and the
!--    neighbor rank for MPI_SENDRECV, pright (pleft) is set to it.
!--    In this way, the calling process receives a new right (left) neighbor
!--    for all future MPI_SENDRECV calls. That neighbor has a y-coordinate
!--    of y+(-)y_shift, where y is the original right (left) neighbor's
!--    y-coordinate. The modulo-operation ensures that if the neighbor's
!--    y-coordinate exceeds the grid-boundary, it will be relocated to
!--    the opposite part of the grid cyclicly.
       IF ( rcoord(1) < pcoord(1) ) THEN
          rcoord(2) = MODULO( rcoord(2) + y_shift, pdims(2) )
          CALL MPI_CART_RANK( comm2d, rcoord, pright, ierr )
       ENDIF

       IF ( lcoord(1) > pcoord(1) ) THEN
          lcoord(2) = MODULO( lcoord(2) - y_shift, pdims(2) )
          CALL MPI_CART_RANK( comm2d, lcoord, pleft, ierr )
       ENDIF
    ENDIF
!
!-- Determine sub-topologies for transpositions
!-- Transposition from z to x:
    remain_dims(1) = .TRUE.
    remain_dims(2) = .FALSE.
    CALL MPI_CART_SUB( comm2d, remain_dims, comm1dx, ierr )
    CALL MPI_COMM_RANK( comm1dx, myidx, ierr )
!
!-- Transposition from x to y
    remain_dims(1) = .FALSE.
    remain_dims(2) = .TRUE.
    CALL MPI_CART_SUB( comm2d, remain_dims, comm1dy, ierr )
    CALL MPI_COMM_RANK( comm1dy, myidy, ierr )


!
!-- Calculate array bounds along x-direction for every PE.
    ALLOCATE( nxlf(0:pdims(1)-1), nxrf(0:pdims(1)-1), nynf(0:pdims(2)-1),      &
              nysf(0:pdims(2)-1) )

    IF ( MOD( nx+1 , pdims(1) ) /= 0 )  THEN
       WRITE( message_string, * ) 'x-direction: gridpoint number (',nx+1,') ', &
                               'is not an& integral divisor of the number ',    &
                               'of processors (', pdims(1),')'
       CALL message( 'init_pegrid', 'PA0225', 1, 2, 0, 6, 0 )
    ELSE
       nnx  = ( nx + 1 ) / pdims(1)
    ENDIF    

!
!-- Left and right array bounds, number of gridpoints
    DO  i = 0, pdims(1)-1
       nxlf(i)   = i * nnx
       nxrf(i)   = ( i + 1 ) * nnx - 1
    ENDDO

!
!-- Calculate array bounds in y-direction for every PE.
    IF ( MOD( ny+1 , pdims(2) ) /= 0 )  THEN
       WRITE( message_string, * ) 'y-direction: gridpoint number (',ny+1,') ', &
                           'is not an& integral divisor of the number of',      &
                           'processors (', pdims(2),')'
       CALL message( 'init_pegrid', 'PA0227', 1, 2, 0, 6, 0 )
    ELSE
       nny  = ( ny + 1 ) / pdims(2)
    ENDIF    

!
!-- South and north array bounds
    DO  j = 0, pdims(2)-1
       nysf(j)   = j * nny
       nynf(j)   = ( j + 1 ) * nny - 1
    ENDDO
!
!-- Local array bounds of the respective PEs
    nxl = nxlf(pcoord(1))
    nxr = nxrf(pcoord(1))
    nys = nysf(pcoord(2))
    nyn = nynf(pcoord(2))
    nzb = 0
    nzt = nz
    nnz = nz

!
!-- Set switches to define if the PE is situated at the border of the virtual
!-- processor grid
    IF ( nxl == 0 )   left_border_pe  = .TRUE.
    IF ( nxr == nx )  right_border_pe = .TRUE.
    IF ( nys == 0 )   south_border_pe = .TRUE.
    IF ( nyn == ny )  north_border_pe = .TRUE.

!
!-- Calculate array bounds and gridpoint numbers for the transposed arrays
!-- (needed in the pressure solver)
!-- For the transposed arrays, cyclic boundaries as well as top and bottom
!-- boundaries are omitted, because they are obstructive to the transposition

!
!-- 1. transposition  z --> x
!-- This transposition is not neccessary in case of a 1d-decomposition along x
    IF ( psolver == 'poisfft' )  THEN

       IF ( pdims(2) /= 1 )  THEN
          IF ( MOD( nz , pdims(1) ) /= 0 )  THEN
             WRITE( message_string, * ) 'transposition z --> x:',              &
                       '& nz=',nz,' is not an integral divisior of pdims(1)=', &
                                                                   pdims(1)
             CALL message( 'init_pegrid', 'PA0230', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       nys_x = nys
       nyn_x = nyn
       nny_x = nny
       nnz_x = nz / pdims(1)
       nzb_x = 1 + myidx * nnz_x
       nzt_x = ( myidx + 1 ) * nnz_x
       sendrecvcount_zx = nnx * nny * nnz_x

    ENDIF


    IF ( psolver == 'poisfft' )  THEN 
!
!--    2. transposition  x --> y
       IF ( MOD( nx+1 , pdims(2) ) /= 0 )  THEN
          WRITE( message_string, * ) 'transposition x --> y:',                 &
                            '& nx+1=',nx+1,' is not an integral divisor of ',  &
                            'pdims(2)=',pdims(2)
          CALL message( 'init_pegrid', 'PA0231', 1, 2, 0, 6, 0 )
       ENDIF

       nnz_y = nnz_x
       nzb_y = nzb_x
       nzt_y = nzt_x
       nnx_y = (nx+1) / pdims(2)
       nxl_y = myidy * nnx_y
       nxr_y = ( myidy + 1 ) * nnx_y - 1
       sendrecvcount_xy = nnx_y * nny_x * nnz_y
!
!--    3. transposition  y --> z  
!--    (ELSE:  x --> y  in case of 1D-decomposition along x)
       nxl_z = nxl_y
       nxr_z = nxr_y
       nny_z = (ny+1) / pdims(1)
       nys_z = myidx * nny_z
       nyn_z = ( myidx + 1 ) * nny_z - 1
       sendrecvcount_yz = nnx_y * nny_z * nnz_y

       IF ( pdims(2) /= 1 )  THEN
!
!--       y --> z
!--       This transposition is not neccessary in case of a 1d-decomposition
!--       along x, except that the uptream-spline method is switched on
          IF ( MOD( ny+1 , pdims(1) ) /= 0 )  THEN
             WRITE( message_string, * ) 'transposition y --> z:',              &
                               '& ny+1=',ny+1,' is not an integral divisor of',&
                               ' pdims(1)=',pdims(1)
             CALL message( 'init_pegrid', 'PA0232', 1, 2, 0, 6, 0 )
          ENDIF

       ELSE
!
!--       x --> y
!--       This condition must be fulfilled for a 1D-decomposition along x
          IF ( MOD( ny+1 , pdims(1) ) /= 0 )  THEN
             WRITE( message_string, * ) 'transposition x --> y:',              &
                               '& ny+1=',ny+1,' is not an integral divisor of',&
                               ' pdims(1)=',pdims(1)
             CALL message( 'init_pegrid', 'PA0233', 1, 2, 0, 6, 0 )
          ENDIF

       ENDIF

    ENDIF
!
!-- Indices for direct transpositions z --> y (used for calculating spectra)
    IF ( psolver == 'poisfft' )  THEN
!
!--    Indices for direct transpositions y --> x 
!--    (they are only possible in case of a 1d-decomposition along x)
       IF ( pdims(2) == 1 )  THEN
          nny_x = nny / pdims(1)
          nys_x = myid * nny_x
          nyn_x = ( myid + 1 ) * nny_x - 1
          nzb_x = 1
          nzt_x = nz
          sendrecvcount_xy = nnx * nny_x * nz
       ENDIF

    ENDIF

    IF ( psolver == 'poisfft' )  THEN
!
!--    Indices for direct transpositions x --> y 
!--    (they are only possible in case of a 1d-decomposition along y)
       IF ( pdims(1) == 1 )  THEN
          nnx_y = nnx / pdims(2)
          nxl_y = myid * nnx_y
          nxr_y = ( myid + 1 ) * nnx_y - 1
          nzb_y = 1
          nzt_y = nz
          sendrecvcount_xy = nnx_y * nny * nz
       ENDIF

    ENDIF

!
!-- Arrays for storing the array bounds are needed any more
    DEALLOCATE( nxlf , nxrf , nynf , nysf )


!
!-- Collect index bounds from other PEs (to be written to restart file later)
    ALLOCATE( hor_index_bounds(4,0:numprocs-1) )

    IF ( myid == 0 )  THEN

       hor_index_bounds(1,0) = nxl
       hor_index_bounds(2,0) = nxr
       hor_index_bounds(3,0) = nys
       hor_index_bounds(4,0) = nyn

!
!--    Receive data from all other PEs
       DO  i = 1, numprocs-1
          CALL MPI_RECV( ibuf, 4, MPI_INTEGER, i, MPI_ANY_TAG, comm2d, status, &
                         ierr )
          hor_index_bounds(:,i) = ibuf(1:4)
       ENDDO

    ELSE
!
!--    Send index bounds to PE0
       ibuf(1) = nxl
       ibuf(2) = nxr
       ibuf(3) = nys
       ibuf(4) = nyn
       CALL MPI_SEND( ibuf, 4, MPI_INTEGER, 0, myid, comm2d, ierr )

    ENDIF


#if defined( __print )
!
!-- Control output
    IF ( myid == 0 )  THEN
       PRINT*, '*** processor topology ***'
       PRINT*, ' '
       PRINT*, 'myid   pcoord    left right  south north  idx idy   nxl: nxr',&
               &'   nys: nyn'
       PRINT*, '------------------------------------------------------------',&
               &'-----------'
       WRITE (*,1000)  0, pcoord(1), pcoord(2), pleft, pright, psouth, pnorth, &
                       myidx, myidy, nxl, nxr, nys, nyn
1000   FORMAT (I4,2X,'(',I3,',',I3,')',3X,I4,2X,I4,3X,I4,2X,I4,2X,I3,1X,I3, &
               2(2X,I4,':',I4))

!
!--    Receive data from the other PEs
       DO  i = 1,numprocs-1
          CALL MPI_RECV( ibuf, 12, MPI_INTEGER, i, MPI_ANY_TAG, comm2d, status, &
                         ierr )
          WRITE (*,1000)  i, ( ibuf(j) , j = 1,12 )
       ENDDO
    ELSE

!
!--    Send data to PE0
       ibuf(1) = pcoord(1); ibuf(2) = pcoord(2); ibuf(3) = pleft
       ibuf(4) = pright; ibuf(5) = psouth; ibuf(6) = pnorth; ibuf(7) = myidx
       ibuf(8) = myidy; ibuf(9) = nxl; ibuf(10) = nxr; ibuf(11) = nys
       ibuf(12) = nyn
       CALL MPI_SEND( ibuf, 12, MPI_INTEGER, 0, myid, comm2d, ierr )       
    ENDIF
#endif

! 
!-- Determine the number of ghost point layers
    IF ( ( scalar_advec == 'ws-scheme' .AND. .NOT. neutral ) .OR.             &
         momentum_advec == 'ws-scheme' )  THEN
       nbgp = 3
    ELSE
       nbgp = 1
    ENDIF 
!
!-- Create a new MPI derived datatype for the exchange of surface (xy) data,
!-- which is needed for coupled atmosphere-ocean runs.
!-- First, calculate number of grid points of an xy-plane.
    ngp_xy  = ( nxr - nxl + 1 + 2 * nbgp ) * ( nyn - nys + 1 + 2 * nbgp )
    CALL MPI_TYPE_VECTOR( ngp_xy, 1, nzt-nzb+2, MPI_REAL, type_xy, ierr )
    CALL MPI_TYPE_COMMIT( type_xy, ierr )

#else

!
!-- Array bounds when running on a single PE (respectively a non-parallel
!-- machine)
    nxl = 0
    nxr = nx
    nnx = nxr - nxl + 1
    nys = 0
    nyn = ny
    nny = nyn - nys + 1
    nzb = 0
    nzt = nz
    nnz = nz

    ALLOCATE( hor_index_bounds(4,0:0) )
    hor_index_bounds(1,0) = nxl
    hor_index_bounds(2,0) = nxr
    hor_index_bounds(3,0) = nys
    hor_index_bounds(4,0) = nyn

!
!-- Array bounds for the pressure solver (in the parallel code, these bounds
!-- are the ones for the transposed arrays)
    nys_x = nys
    nyn_x = nyn
    nzb_x = nzb + 1
    nzt_x = nzt

    nxl_y = nxl
    nxr_y = nxr
    nzb_y = nzb + 1
    nzt_y = nzt

    nxl_z = nxl
    nxr_z = nxr
    nys_z = nys
    nyn_z = nyn

#endif

!
!-- Calculate number of grid levels necessary for the multigrid poisson solver
!-- as well as the gridpoint indices on each level
    IF ( psolver(1:9) == 'multigrid' )  THEN

!
!--    First calculate number of possible grid levels for the subdomains
       mg_levels_x = 1
       mg_levels_y = 1
       mg_levels_z = 1

       i = nnx
       DO WHILE ( MOD( i, 2 ) == 0  .AND.  i /= 2 )
          i = i / 2
          mg_levels_x = mg_levels_x + 1
       ENDDO

       j = nny
       DO WHILE ( MOD( j, 2 ) == 0  .AND.  j /= 2 )
          j = j / 2
          mg_levels_y = mg_levels_y + 1
       ENDDO

       k = nz    ! do not use nnz because it might be > nz due to transposition
                 ! requirements
       DO WHILE ( MOD( k, 2 ) == 0  .AND.  k /= 2 )
          k = k / 2
          mg_levels_z = mg_levels_z + 1
       ENDDO
!
!--    The optimized MG-solver does not allow odd values for nz at the coarsest
!--    grid level
       IF ( TRIM( psolver ) /= 'multigrid_noopt' )  THEN
          IF ( MOD( k, 2 ) /= 0 )  mg_levels_z = mg_levels_z - 1
!
!--       An odd value of nz does not work. The finest level must have an even
!--       value.
          IF (  mg_levels_z == 0 )  THEN
             message_string = 'optimized multigrid method requires nz to be even'
             CALL message( 'init_pegrid', 'PA0495', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       maximum_grid_level = MIN( mg_levels_x, mg_levels_y, mg_levels_z )

!
!--    Find out, if the total domain allows more levels. These additional
!--    levels are identically processed on all PEs.
       IF ( numprocs > 1  .AND.  mg_switch_to_pe0_level /= -1 )  THEN

          IF ( mg_levels_z > MIN( mg_levels_x, mg_levels_y ) )  THEN

             mg_switch_to_pe0_level_l = maximum_grid_level

             mg_levels_x = 1
             mg_levels_y = 1

             i = nx+1
             DO WHILE ( MOD( i, 2 ) == 0  .AND.  i /= 2 )
                i = i / 2
                mg_levels_x = mg_levels_x + 1
             ENDDO

             j = ny+1
             DO WHILE ( MOD( j, 2 ) == 0  .AND.  j /= 2 )
                j = j / 2
                mg_levels_y = mg_levels_y + 1
             ENDDO

             maximum_grid_level_l = MIN( mg_levels_x, mg_levels_y, mg_levels_z )

             IF ( maximum_grid_level_l > mg_switch_to_pe0_level_l )  THEN
                mg_switch_to_pe0_level_l = maximum_grid_level_l - &
                                           mg_switch_to_pe0_level_l + 1
             ELSE
                mg_switch_to_pe0_level_l = 0
             ENDIF

          ELSE

             mg_switch_to_pe0_level_l = 0
             maximum_grid_level_l = maximum_grid_level

          ENDIF

!
!--       Use switch level calculated above only if it is not pre-defined
!--       by user
          IF ( mg_switch_to_pe0_level == 0 )  THEN
             IF ( mg_switch_to_pe0_level_l /= 0 )  THEN
                mg_switch_to_pe0_level = mg_switch_to_pe0_level_l
                maximum_grid_level     = maximum_grid_level_l
             ENDIF

          ELSE
!
!--          Check pre-defined value and reset to default, if neccessary
             IF ( mg_switch_to_pe0_level < mg_switch_to_pe0_level_l  .OR.      &
                  mg_switch_to_pe0_level >= maximum_grid_level_l )  THEN
                message_string = 'mg_switch_to_pe0_level ' //                  &
                                 'out of range and reset to 0'
                CALL message( 'init_pegrid', 'PA0235', 0, 1, 0, 6, 0 )
                mg_switch_to_pe0_level = 0
             ELSE
!
!--             Use the largest number of possible levels anyway and recalculate
!--             the switch level to this largest number of possible values
                maximum_grid_level = maximum_grid_level_l

             ENDIF

          ENDIF

       ENDIF

       ALLOCATE( grid_level_count(maximum_grid_level),                       &
                 nxl_mg(0:maximum_grid_level), nxr_mg(0:maximum_grid_level), &
                 nyn_mg(0:maximum_grid_level), nys_mg(0:maximum_grid_level), &
                 nzt_mg(0:maximum_grid_level) )

       grid_level_count = 0
!
!--    Index zero required as dummy due to definition of arrays f2 and p2 in
!--    recursive subroutine next_mg_level
       nxl_mg(0) = 0; nxr_mg(0) = 0; nyn_mg(0) = 0; nys_mg(0) = 0; nzt_mg(0) = 0

       nxl_l = nxl; nxr_l = nxr; nys_l = nys; nyn_l = nyn; nzt_l = nzt

       DO  i = maximum_grid_level, 1 , -1

          IF ( i == mg_switch_to_pe0_level )  THEN
#if defined( __parallel )
!
!--          Save the grid size of the subdomain at the switch level, because
!--          it is needed in poismg.
             ind(1) = nxl_l; ind(2) = nxr_l
             ind(3) = nys_l; ind(4) = nyn_l
             ind(5) = nzt_l
             ALLOCATE( ind_all(5*numprocs), mg_loc_ind(5,0:numprocs-1) )
             CALL MPI_ALLGATHER( ind, 5, MPI_INTEGER, ind_all, 5, &
                                 MPI_INTEGER, comm2d, ierr )
             DO  j = 0, numprocs-1
                DO  k = 1, 5
                   mg_loc_ind(k,j) = ind_all(k+j*5)
                ENDDO
             ENDDO
             DEALLOCATE( ind_all )
!
!--          Calculate the grid size of the total domain
             nxr_l = ( nxr_l-nxl_l+1 ) * pdims(1) - 1
             nxl_l = 0
             nyn_l = ( nyn_l-nys_l+1 ) * pdims(2) - 1
             nys_l = 0
!
!--          The size of this gathered array must not be larger than the
!--          array tend, which is used in the multigrid scheme as a temporary
!--          array. Therefore the subdomain size of an PE is calculated and 
!--          the size of the gathered grid. These values are used in  
!--          routines pres and poismg
             subdomain_size = ( nxr - nxl + 2 * nbgp + 1 ) * &
                              ( nyn - nys + 2 * nbgp + 1 ) * ( nzt - nzb + 2 )
             gathered_size  = ( nxr_l - nxl_l + 3 ) * ( nyn_l - nys_l + 3 ) *  &
                              ( nzt_l - nzb + 2 )

#else
             message_string = 'multigrid gather/scatter impossible ' //        &
                          'in non parallel mode'
             CALL message( 'init_pegrid', 'PA0237', 1, 2, 0, 6, 0 )
#endif
          ENDIF

          nxl_mg(i) = nxl_l
          nxr_mg(i) = nxr_l
          nys_mg(i) = nys_l
          nyn_mg(i) = nyn_l
          nzt_mg(i) = nzt_l

          nxl_l = nxl_l / 2 
          nxr_l = nxr_l / 2
          nys_l = nys_l / 2 
          nyn_l = nyn_l / 2 
          nzt_l = nzt_l / 2 

       ENDDO

!
!--    Temporary problem: Currently calculation of maxerror in routine poismg crashes
!--    if grid data are collected on PE0 already on the finest grid level.
!--    To be solved later.
       IF ( maximum_grid_level == mg_switch_to_pe0_level )  THEN
          message_string = 'grid coarsening on subdomain level cannot be performed'
          CALL message( 'poismg', 'PA0236', 1, 2, 0, 6, 0 )
       ENDIF

    ELSE

       maximum_grid_level = 0

    ENDIF

!
!-- Default level 0 tells exchange_horiz that all ghost planes have to be
!-- exchanged. grid_level is adjusted in poismg, where only one ghost plane
!-- is required.
    grid_level = 0

#if defined( __parallel )
!
!-- Gridpoint number for the exchange of ghost points (y-line for 2D-arrays)
    ngp_y  = nyn - nys + 1 + 2 * nbgp

!
!-- Define new MPI derived datatypes for the exchange of ghost points in
!-- x- and y-direction for 2D-arrays (line)
    CALL MPI_TYPE_VECTOR( nxr-nxl+1+2*nbgp, nbgp, ngp_y, MPI_REAL, type_x,     &
                          ierr )
    CALL MPI_TYPE_COMMIT( type_x, ierr )

    CALL MPI_TYPE_VECTOR( nbgp, ngp_y, ngp_y, MPI_REAL, type_y, ierr )
    CALL MPI_TYPE_COMMIT( type_y, ierr )
!
!-- Define new MPI derived datatypes for the exchange of ghost points in
!-- x- and y-direction for 2D-INTEGER arrays (line) - on normal grid
    ALLOCATE( type_x_int(0:maximum_grid_level),                                &
              type_y_int(0:maximum_grid_level) )

    CALL MPI_TYPE_VECTOR( nxr-nxl+1+2*nbgp, nbgp, ngp_y, MPI_INTEGER,          &
                          type_x_int(0), ierr )
    CALL MPI_TYPE_COMMIT( type_x_int(0), ierr )

    CALL MPI_TYPE_VECTOR( nbgp, ngp_y, ngp_y, MPI_INTEGER, type_y_int(0), ierr )
    CALL MPI_TYPE_COMMIT( type_y_int(0), ierr )
!
!-- Calculate gridpoint numbers for the exchange of ghost points along x
!-- (yz-plane for 3D-arrays) and define MPI derived data type(s) for the
!-- exchange of ghost points in y-direction (xz-plane).
!-- Do these calculations for the model grid and (if necessary) also
!-- for the coarser grid levels used in the multigrid method
    ALLOCATE ( ngp_xz(0:maximum_grid_level),                                   &
               ngp_xz_int(0:maximum_grid_level),                               &
               ngp_yz(0:maximum_grid_level),                                   &
               ngp_yz_int(0:maximum_grid_level),                               &
               type_xz(0:maximum_grid_level),                                  &
               type_xz_int(0:maximum_grid_level),                              &
               type_yz(0:maximum_grid_level),                                  &
               type_yz_int(0:maximum_grid_level) )

    nxl_l = nxl; nxr_l = nxr; nys_l = nys; nyn_l = nyn; nzb_l = nzb; nzt_l = nzt

!
!-- Discern between the model grid, which needs nbgp ghost points and
!-- grid levels for the multigrid scheme. In the latter case only one
!-- ghost point is necessary.
!-- First definition of MPI-datatypes for exchange of ghost layers on normal 
!-- grid. The following loop is needed for data exchange in poismg.f90.
!
!-- Determine number of grid points of yz-layer for exchange
    ngp_yz(0) = (nzt - nzb + 2) * (nyn - nys + 1 + 2 * nbgp)

!
!-- Define an MPI-datatype for the exchange of left/right boundaries.
!-- Although data are contiguous in physical memory (which does not
!-- necessarily require an MPI-derived datatype), the data exchange between
!-- left and right PE's using the MPI-derived type is 10% faster than without.
    CALL MPI_TYPE_VECTOR( nxr-nxl+1+2*nbgp, nbgp*(nzt-nzb+2), ngp_yz(0), &
                          MPI_REAL, type_xz(0), ierr )
    CALL MPI_TYPE_COMMIT( type_xz(0), ierr )

    CALL MPI_TYPE_VECTOR( nbgp, ngp_yz(0), ngp_yz(0), MPI_REAL, type_yz(0), &
                          ierr ) 
    CALL MPI_TYPE_COMMIT( type_yz(0), ierr )

!
!-- Define data types for exchange of 3D Integer arrays.
    ngp_yz_int(0) = (nzt - nzb + 2) * (nyn - nys + 1 + 2 * nbgp)

    CALL MPI_TYPE_VECTOR( nxr-nxl+1+2*nbgp, nbgp*(nzt-nzb+2), ngp_yz_int(0),   &
                          MPI_INTEGER, type_xz_int(0), ierr )
    CALL MPI_TYPE_COMMIT( type_xz_int(0), ierr )

    CALL MPI_TYPE_VECTOR( nbgp, ngp_yz_int(0), ngp_yz_int(0), MPI_INTEGER,     &
                          type_yz_int(0), ierr )
    CALL MPI_TYPE_COMMIT( type_yz_int(0), ierr )

!
!-- Definition of MPI-datatypes for multigrid method (coarser level grids)
    IF ( psolver(1:9) == 'multigrid' )  THEN
!    
!--    Definition of MPI-datatyoe as above, but only 1 ghost level is used
       DO  i = maximum_grid_level, 1 , -1
!
!--       For 3D-exchange on different multigrid level, one ghost point for 
!--       REAL arrays, two ghost points for INTEGER arrays
          ngp_xz(i) = (nzt_l - nzb_l + 2) * (nxr_l - nxl_l + 3)
          ngp_yz(i) = (nzt_l - nzb_l + 2) * (nyn_l - nys_l + 3)

          ngp_xz_int(i) = (nzt_l - nzb_l + 2) * (nxr_l - nxl_l + 3)
          ngp_yz_int(i) = (nzt_l - nzb_l + 2) * (nyn_l - nys_l + 3)
!
!--       MPI data type for REAL arrays, for xz-layers
          CALL MPI_TYPE_VECTOR( nxr_l-nxl_l+3, nzt_l-nzb_l+2, ngp_yz(i),       &
                                MPI_REAL, type_xz(i), ierr )
          CALL MPI_TYPE_COMMIT( type_xz(i), ierr )

!
!--       MPI data type for INTEGER arrays, for xz-layers
          CALL MPI_TYPE_VECTOR( nxr_l-nxl_l+3, nzt_l-nzb_l+2, ngp_yz_int(i),   &
                                MPI_INTEGER, type_xz_int(i), ierr )
          CALL MPI_TYPE_COMMIT( type_xz_int(i), ierr )

!
!--       MPI data type for REAL arrays, for yz-layers
          CALL MPI_TYPE_VECTOR( 1, ngp_yz(i), ngp_yz(i), MPI_REAL, type_yz(i), &
                                ierr )
          CALL MPI_TYPE_COMMIT( type_yz(i), ierr )
!
!--       MPI data type for INTEGER arrays, for yz-layers
          CALL MPI_TYPE_VECTOR( 1, ngp_yz_int(i), ngp_yz_int(i), MPI_INTEGER,  &
                                type_yz_int(i), ierr )
          CALL MPI_TYPE_COMMIT( type_yz_int(i), ierr )


!--       For 2D-exchange of INTEGER arrays on coarser grid level, where 2 ghost
!--       points need to be exchanged.
          CALL MPI_TYPE_VECTOR( nxr_l-nxl_l+5, 2, nyn_l-nys_l+5, MPI_INTEGER,          &
                                type_x_int(i), ierr )
          CALL MPI_TYPE_COMMIT( type_x_int(i), ierr )


          CALL MPI_TYPE_VECTOR( 2, nyn_l-nys_l+5, nyn_l-nys_l+5, MPI_INTEGER,          &
                                type_y_int(i), ierr )
          CALL MPI_TYPE_COMMIT( type_y_int(i), ierr )

          nxl_l = nxl_l / 2
          nxr_l = nxr_l / 2
          nys_l = nys_l / 2
          nyn_l = nyn_l / 2
          nzt_l = nzt_l / 2

       ENDDO

    ENDIF

#endif

#if defined( __parallel )
!
!-- Setting of flags for inflow/outflow/nesting conditions.
    IF ( pleft == MPI_PROC_NULL )  THEN
       IF ( bc_lr == 'dirichlet/radiation' )  THEN
          inflow_l  = .TRUE.
       ELSEIF ( bc_lr == 'radiation/dirichlet' )  THEN
          outflow_l = .TRUE.
       ELSEIF ( bc_lr == 'nested' )  THEN
          nest_bound_l = .TRUE.
       ELSEIF ( bc_lr == 'forcing' )  THEN
          force_bound_l = .TRUE.
       ENDIF
    ENDIF
 
    IF ( pright == MPI_PROC_NULL )  THEN
       IF ( bc_lr == 'dirichlet/radiation' )  THEN
          outflow_r = .TRUE.
       ELSEIF ( bc_lr == 'radiation/dirichlet' )  THEN
          inflow_r  = .TRUE.
       ELSEIF ( bc_lr == 'nested' )  THEN
          nest_bound_r = .TRUE.
       ELSEIF ( bc_lr == 'forcing' )  THEN
          force_bound_r = .TRUE.
       ENDIF
    ENDIF

    IF ( psouth == MPI_PROC_NULL )  THEN
       IF ( bc_ns == 'dirichlet/radiation' )  THEN
          outflow_s = .TRUE.
       ELSEIF ( bc_ns == 'radiation/dirichlet' )  THEN
          inflow_s  = .TRUE.
       ELSEIF ( bc_ns == 'nested' )  THEN
          nest_bound_s = .TRUE.
       ELSEIF ( bc_ns == 'forcing' )  THEN
          force_bound_s = .TRUE.
       ENDIF
    ENDIF

    IF ( pnorth == MPI_PROC_NULL )  THEN
       IF ( bc_ns == 'dirichlet/radiation' )  THEN
          inflow_n  = .TRUE.
       ELSEIF ( bc_ns == 'radiation/dirichlet' )  THEN
          outflow_n = .TRUE.
       ELSEIF ( bc_ns == 'nested' )  THEN
          nest_bound_n = .TRUE.
       ELSEIF ( bc_ns == 'forcing' )  THEN
          force_bound_n = .TRUE.
       ENDIF
    ENDIF
!
!
!-- Broadcast the id of the inflow PE
    IF ( inflow_l )  THEN
       id_inflow_l = myidx
    ELSE
       id_inflow_l = 0
    ENDIF
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( id_inflow_l, id_inflow, 1, MPI_INTEGER, MPI_SUM, &
                        comm1dx, ierr )

!
!-- Broadcast the id of the recycling plane
!-- WARNING: needs to be adjusted in case of inflows other than from left side!
    IF ( NINT( recycling_width / dx ) >= nxl  .AND. &
         NINT( recycling_width / dx ) <= nxr )  THEN
       id_recycling_l = myidx
    ELSE
       id_recycling_l = 0
    ENDIF
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( id_recycling_l, id_recycling, 1, MPI_INTEGER, MPI_SUM, &
                        comm1dx, ierr )
    CALL location_message( 'finished', .TRUE. )

#else
    IF ( bc_lr == 'dirichlet/radiation' )  THEN
       inflow_l  = .TRUE.
       outflow_r = .TRUE.
    ELSEIF ( bc_lr == 'radiation/dirichlet' )  THEN
       outflow_l = .TRUE.
       inflow_r  = .TRUE.
    ENDIF

    IF ( bc_ns == 'dirichlet/radiation' )  THEN
       inflow_n  = .TRUE.
       outflow_s = .TRUE.
    ELSEIF ( bc_ns == 'radiation/dirichlet' )  THEN
       outflow_n = .TRUE.
       inflow_s  = .TRUE.
    ENDIF
#endif

!
!-- At the inflow or outflow, u or v, respectively, have to be calculated for
!-- one more grid point.
    IF ( inflow_l .OR. outflow_l .OR. nest_bound_l .OR. force_bound_l )  THEN
       nxlu = nxl + 1
    ELSE
       nxlu = nxl
    ENDIF
    IF ( inflow_s .OR. outflow_s .OR. nest_bound_s .OR. force_bound_s )  THEN
       nysv = nys + 1
    ELSE
       nysv = nys
    ENDIF

 END SUBROUTINE init_pegrid
