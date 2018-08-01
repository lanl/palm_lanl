 MODULE pmc_general

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
! $Id: pmc_general_mod.f90 2841 2018-02-27 15:02:57Z knoop $
! Bugfix: wrong placement of include 'mpif.h' corrected
! 
! 2801 2018-02-14 16:01:55Z thiele
! Introduce particle transfer in nested models.
! 
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2599 2017-11-01 13:18:45Z hellstea
! Some cleanup and commenting improvements only.
! 
! 2101 2017-01-05 16:42:31Z suehring
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1901 2016-05-04 15:39:38Z raasch
! Code clean up. The words server/client changed to parent/child.
!
! 1900 2016-05-04 15:27:53Z raasch
! re-formatted to match PALM style, file renamed again
!
! 1808 2016-04-05 19:44:00Z raasch
! MPI module used by default on all machines
!
! 1786 2016-03-08 05:49:27Z raasch
! change in child-parent data transfer: parent now gets data from child
! instead of that child puts it to the parent
!
! 1779 2016-03-03 08:01:28Z raasch
! PMC_MPI_REAL removed, dim_order removed from type arraydef,
! array management changed from linked list to sequential loop
!
! 1766 2016-02-29 08:37:15Z raasch
! +po_data in type arraydef
!
! 1764 2016-02-28 12:45:19Z raasch
! cpp-statement added (nesting can only be used in parallel mode),
! all kinds given in PALM style
!
! 1762 2016-02-25 12:31:13Z hellstea
! Initial revision by K. Ketelsen
!
! Description:
! ------------
!
! Structure definition and utilities of Palm Model Coupler
!------------------------------------------------------------------------------!

#if defined( __parallel )
    USE, INTRINSIC ::  ISO_C_BINDING

    USE kinds

#if !defined( __mpifh )
    USE MPI
#endif

    IMPLICIT NONE

#if defined( __mpifh )
    INCLUDE "mpif.h"
#endif

    PRIVATE
    SAVE

    INTEGER(iwp), PARAMETER, PUBLIC :: da_desclen       =  8  !<
    INTEGER(iwp), PARAMETER, PUBLIC :: da_namelen       = 16  !<
    INTEGER(iwp), PARAMETER, PUBLIC :: pmc_da_name_err  = 10  !<
    INTEGER(iwp), PARAMETER, PUBLIC :: pmc_max_array    = 32  !< max # of arrays which can be coupled
    INTEGER(iwp), PARAMETER, PUBLIC :: pmc_max_models   = 64  !<
    INTEGER(iwp), PARAMETER, PUBLIC :: pmc_status_ok    =  0  !<
    INTEGER(iwp), PARAMETER, PUBLIC :: pmc_status_error = -1  !<


    TYPE, PUBLIC :: xy_ind  !< pair of indices in horizontal plane
       INTEGER(iwp) ::  i
       INTEGER(iwp) ::  j
    END TYPE

    TYPE, PUBLIC ::  arraydef
       INTEGER(iwp)                   :: coupleindex  !<
       INTEGER(iwp)                   :: nrdims       !< number of dimensions
       INTEGER(iwp)                   :: dimkey       !< key for NR dimensions and array type
       INTEGER(iwp), DIMENSION(4)     :: a_dim        !< size of dimensions
       TYPE(C_PTR)               :: data         !< pointer of data in parent space
       TYPE(C_PTR), DIMENSION(2) :: po_data      !< base pointers,
                                                 !< pmc_s_set_active_data_array
                                                 !< sets active pointer
       INTEGER(idp)              :: SendIndex    !< index in send buffer
       INTEGER(idp)              :: RecvIndex    !< index in receive buffer
       INTEGER(iwp)              :: SendSize     !< size in send buffer
       INTEGER(iwp)              :: RecvSize     !< size in receive buffer
       TYPE(C_PTR)               :: SendBuf      !< data pointer in send buffer
       TYPE(C_PTR)               :: RecvBuf      !< data pointer in receive buffer
       CHARACTER(LEN=8)          :: Name         !< name of array
       TYPE(arraydef), POINTER   :: next
    END TYPE arraydef

    TYPE(arraydef), PUBLIC, POINTER  :: next

    TYPE, PUBLIC ::  pedef
       INTEGER(iwp) :: nr_arrays = 0  !< number of arrays which will be transfered
       INTEGER(iwp) :: nrele          !< number of elements, same for all arrays
       TYPE(xy_ind), POINTER, DIMENSION(:)   ::  locInd      !< xy index local array for remote PE
       TYPE(arraydef), POINTER, DIMENSION(:) ::  array_list  !< list of data arrays to be transfered
    END TYPE pedef

    TYPE, PUBLIC ::  childdef
       INTEGER(idp) ::  totalbuffersize    !<
       INTEGER(iwp) ::  model_comm         !< communicator of this model
       INTEGER(iwp) ::  inter_comm         !< inter communicator model and child
       INTEGER(iwp) ::  intra_comm         !< intra communicator model and child
       INTEGER(iwp) ::  model_rank         !< rank of this model
       INTEGER(iwp) ::  model_npes         !< number of PEs this model
       INTEGER(iwp) ::  inter_npes         !< number of PEs child model
       INTEGER(iwp) ::  intra_rank         !< rank within intra_comm
       INTEGER(iwp) ::  win_parent_child   !< MPI RMA for preparing data on parent AND child side
       TYPE(pedef), DIMENSION(:), POINTER ::  pes  !< list of all child PEs
    END TYPE childdef

    TYPE, PUBLIC ::  da_namedef  !< data array name definition
       INTEGER(iwp)              ::  couple_index  !< unique number of array
       CHARACTER(LEN=da_desclen) ::  parentdesc    !< parent array description
       CHARACTER(LEN=da_namelen) ::  nameonparent  !< name of array within parent
       CHARACTER(LEN=da_desclen) ::  childdesc     !< child array description
       CHARACTER(LEN=da_namelen) ::  nameonchild   !< name of array within child
    END TYPE da_namedef

    INTERFACE pmc_g_setname
       MODULE PROCEDURE pmc_g_setname
    END INTERFACE pmc_g_setname

    INTERFACE pmc_sort
       MODULE PROCEDURE sort_2d_i
    END INTERFACE pmc_sort

    PUBLIC pmc_g_setname, pmc_sort

 CONTAINS


   
 SUBROUTINE pmc_g_setname( mychild, couple_index, aname )

    IMPLICIT NONE

    CHARACTER(LEN=*)              ::  aname         !<
    INTEGER(iwp), INTENT(IN)      ::  couple_index  !<
    TYPE(childdef), INTENT(INOUT) ::  mychild       !<

    INTEGER(iwp) ::  i  !<

    TYPE(arraydef), POINTER ::  ar   !<
    TYPE(pedef), POINTER    ::  ape  !<

!
!-- Assign array to next free index in array list.
!-- Set name of array in arraydef structure
    DO  i = 1, mychild%inter_npes
       ape => mychild%pes(i)
       ape%nr_arrays = ape%nr_arrays + 1
       ape%array_list(ape%nr_arrays)%name        = aname
       ape%array_list(ape%nr_arrays)%coupleindex = couple_index
    ENDDO

 END SUBROUTINE pmc_g_setname



 SUBROUTINE sort_2d_i( array, sort_ind )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN)                    ::  sort_ind
    INTEGER(iwp), DIMENSION(:,:), INTENT(INOUT) ::  array

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  n  !<

    INTEGER(iwp), DIMENSION(SIZE(array,1)) ::  tmp  !<

    n = SIZE(array,2)
    DO  j = 1, n-1
       DO  i = j+1, n
          IF ( array(sort_ind,i) < array(sort_ind,j) )  THEN
             tmp = array(:,i)
             array(:,i) = array(:,j)
             array(:,j) = tmp
          ENDIF
       ENDDO
    ENDDO

 END  SUBROUTINE sort_2d_i

#endif
 END MODULE pmc_general
