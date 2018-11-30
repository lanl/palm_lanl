!> @file close_file.f90
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
! $Id: close_file.f90 3045 2018-05-28 07:55:41Z Giersch $
! z_max_do2d removed and output case 108 disabled
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2300 2017-06-29 13:31:14Z raasch
! -host
! 
! 2277 2017-06-12 10:47:51Z kanani
! Removed unused variables do2d_xy_n, do2d_xz_n, do2d_yz_n, do3d_avs_n
! 
! 2101 2017-01-05 16:42:31Z suehring
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1992 2016-08-12 15:14:59Z suehring
! -Close file containing flight data
! -Some tabs removed.
!
! 1783 2016-03-06 18:36:17Z raasch
! name change of netcdf routines and module + related changes
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1327 2014-03-21 11:00:16Z raasch
! parts concerning iso2d and avs output removed
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
! 1031 2012-10-19 14:35:30Z raasch
! netCDF4 without parallel file support implemented
!
! 964 2012-07-26 09:14:24Z raasch
! old profil-units (40:49) and respective code removed
!
! Revision 1.1 (close_files) 1997/08/11 06:11:18  raasch
! Initial revision
!
!
! Description:
! ------------
!> Close specified file or all open files, if "0" has been given as the
!> calling argument. In that case, execute last actions for certain unit 
!> numbers, if required.
!------------------------------------------------------------------------------!
 SUBROUTINE close_file( file_id )
 

    USE control_parameters,                                                    &
        ONLY:  max_masks, mid, nz_do3d, openfile, run_description_header
                
    USE grid_variables,                                                        &
        ONLY:  dy
        
    USE indices,                                                               &
        ONLY:  nx, ny, nz
        
    USE kinds
    
#if defined( __netcdf )
    USE NETCDF
#endif

    USE netcdf_interface,                                                      &
        ONLY:  id_set_mask, id_set_pr, id_set_prt, id_set_pts, id_set_sp,      &
               id_set_ts, id_set_3d,          &
               nc_stat, netcdf_data_format, netcdf_handle_error
                
    USE pegrid                                           

    IMPLICIT NONE

    CHARACTER (LEN=10)  ::  datform = 'lit_endian' !< 
    CHARACTER (LEN=80)  ::  title                  !< 

    INTEGER(iwp) ::  av           !< 
    INTEGER(iwp) ::  dimx         !< 
    INTEGER(iwp) ::  dimy         !< 
    INTEGER(iwp) ::  fid          !< 
    INTEGER(iwp) ::  file_id      !< 
    INTEGER(iwp) ::  planz        !< 

    LOGICAL ::  checkuf = .TRUE.  !< 
    LOGICAL ::  datleg = .TRUE.   !< 
    LOGICAL ::  dbp = .FALSE.     !< 

    REAL(wp) ::  sizex            !< 
    REAL(wp) ::  sizey            !< 
    REAL(wp) ::  yright           !< 

    NAMELIST /GLOBAL/  checkuf, datform, dimx, dimy, dbp, planz,               &
                       title
    NAMELIST /RAHMEN/  datleg

!
!-- Close specified unit number (if opened) and set a flag that it has
!-- been opened one time at least
    IF ( file_id /= 0 )  THEN
       IF ( openfile(file_id)%opened )  THEN
          CLOSE ( file_id )
          openfile(file_id)%opened        = .FALSE.
          openfile(file_id)%opened_before = .TRUE.
       ENDIF
       RETURN
    ENDIF

!
!-- Close all open unit numbers
    DO  fid = 1, 200+2*max_masks

       IF ( openfile(fid)%opened .OR. openfile(fid)%opened_before )  THEN
!
!--       Last actions for certain unit numbers
          SELECT CASE ( fid )

#if defined( __netcdf )
             CASE ( 104 )

                IF ( myid == 0 )  THEN
                   nc_stat = NF90_CLOSE( id_set_pr )
                   CALL netcdf_handle_error( 'close_file', 47 )
                ENDIF

             CASE ( 105 )

                IF ( myid == 0 )  THEN
                   nc_stat = NF90_CLOSE( id_set_ts )
                   CALL netcdf_handle_error( 'close_file', 48 )
                ENDIF

             CASE ( 106 )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_CLOSE( id_set_3d(0) )
                   CALL netcdf_handle_error( 'close_file', 49 )
                ENDIF

            CASE ( 109 ) 

                nc_stat = NF90_CLOSE( id_set_pts )
                CALL netcdf_handle_error( 'close_file', 412 )

            CASE ( 116 )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_CLOSE( id_set_3d(1) )
                   CALL netcdf_handle_error( 'close_file', 353 )
                ENDIF

#endif

          END SELECT
!
!--       Close file
          IF ( openfile(fid)%opened )  CLOSE ( fid )

       ENDIF

    ENDDO

!
!-- Formats
3200 FORMAT ('# AVS',A,'field file'/                                           &
             '#'/                                                              &
             '# ',A/                                                           &
             'ndim=3'/                                                         &
             'dim1=',I5/                                                       &
             'dim2=',I5/                                                       &
             'dim3=',I5/                                                       &
             'nspace=3'/                                                       &
             'veclen=',I5/                                                     &
             'data=xdr_float'/                                                 &
             'field=rectilinear')
4000 FORMAT ('time averaged over',F7.1,' s')


 END SUBROUTINE close_file
