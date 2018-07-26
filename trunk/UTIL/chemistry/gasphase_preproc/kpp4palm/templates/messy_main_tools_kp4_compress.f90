module messy_main_tools_kp4_compress
!
! Module to compress data for ros_Integrator using variable time step 
! in vector mode. This module is not dependent on chemistry setup and
! can be use for multiple chemistry setups in one job.
!

!     COPYRIGHT Klaus Ketelsen and MPI-CH              May 2007


  implicit none
  private
  save

! KPP DP - Double precision kind
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)


! variables
  integer                                   :: kpoints
  integer                                   :: kpoints_save
  integer,dimension(:),allocatable,public   :: index_org
  integer,dimension(:),allocatable          :: index_step
  logical,dimension(:),allocatable,public   :: cell_done
  logical,dimension(:),allocatable          :: done_check

! Data arrays
  real(kind=dp),dimension(:,:),allocatable  :: f_done
  integer,      dimension(:),  allocatable  :: Kacc_done, Krej_done
  integer,      dimension(:),  allocatable  :: ierr_done

! interface block

  interface kco_initialize
    module procedure kco_initialize
  end interface kco_initialize

  interface kco_compress
    module procedure kco_compress
  end interface kco_compress

  interface kco_finalize
    module procedure kco_finalize
  end interface kco_finalize

  public kco_initialize,kco_compress
  public kco_finalize

 contains
! Public subroutines

  subroutine kco_initialize (npoints_initial, f, IERR, Kacc, Krej)
    implicit    none
    integer,intent(IN)                                  :: npoints_initial
    real(kind=dp),intent(IN),dimension(:,:)             :: f     !Only to get dimension for f_done
    integer,intent(IN),dimension(:)                     :: IERR
    integer,intent(IN),dimension(:)                     :: Kacc, Krej

!-- local variables
    integer                          :: jl

    kpoints      = npoints_initial
    kpoints_save = kpoints

    allocate(index_org(kpoints))
    allocate(done_check(kpoints))
    allocate(index_step(kpoints))
    allocate(cell_done(kpoints))
    
    allocate(f_done(size(f,1),size(f,2)))
    allocate(Kacc_done(size(Kacc)))
    allocate(Krej_done(size(Krej)))
    allocate(ierr_done(size(ierr)))

    cell_done     = .false.
    done_check    = .false.

    do jl=1,kpoints
      index_org(jl) = jl
    end do

    return
  end subroutine kco_initialize

  subroutine kco_compress (n_points, T, H, Hnew, IERR, f, RCONST, FIX,     &
                                   RejectLastH, RejectMoreH, Kacc, Krej)
    implicit    none
    integer,intent(INOUT)                       :: n_points
    real(kind=dp),intent(INOUT),dimension(:)    :: T, H, Hnew
    INTEGER,intent(INOUT),dimension(:)          :: IERR
    real(kind=dp),intent(INOUT),dimension(:,:)  :: f
    real(kind=dp),intent(INOUT),dimension(:,:)  :: RCONST, FIX
    LOGICAL,intent(INOUT),dimension(:)          :: RejectLastH, RejectMoreH
    integer,intent(INOUT),dimension(:)          :: Kacc, Krej

!-- local variables
    integer                          :: i,jl
    integer                          :: ic

!   Check, if converged variables exist

    if(.not. any(cell_done(1:kpoints)))   return

!   Move converged cells to original position

    do i=1,size(f,2)
!CDIR NODEP
      do jl=1,kpoints
        if(cell_done(jl))   then
          f_done(index_org(jl),i) = f(jl,i) 
        end if
      end do
    end do
    
    do jl=1,kpoints
      if(cell_done(jl))   then
        ierr_done(index_org(jl)) = ierr(jl) 
        Kacc_done(index_org(jl)) = Kacc(jl) 
        Krej_done(index_org(jl)) = Krej(jl) 
      end if
    end do

!CDIR NODEP
    do jl=1,kpoints
      if(cell_done(jl))   then
        done_check(index_org(jl)) = .true.
      end if
    end do
!kk    write(6,*) 'done check ',kpoints,count(done_check)
        
!   Compute new index vector

    ic = 0
    do jl=1,kpoints
      if(.not.cell_done(jl))   then
        ic = ic+1
        index_step(ic) = jl
      end if
    end do

!kk    write(6,*) 'new number of points: ',ic

!   Set new number of Points

    kpoints  = ic
    n_points = ic

!   Compress arrays

    index_org(1:ic) = index_org(index_step(1:ic))

    do i=1,size(f,2)
      f(1:ic,i) = f(index_step(1:ic),i)
    end do

!   compress argument list arrays

    T(1:ic)           = T(index_step(1:ic))
    IERR(1:ic)        = IERR(index_step(1:ic))
    H(1:ic)           = H(index_step(1:ic))
    Hnew(1:ic)        = Hnew(index_step(1:ic))

    do i=1,size(RCONST,2)
       RCONST(1:ic,i) = RCONST(index_step(1:ic),i)
    end do

    do i=1,size(FIX,2)
       FIX(1:ic,i)    = FIX(index_step(1:ic),i)
    end do

    RejectLastH(1:ic) = RejectLastH(index_step(1:ic))
    RejectMoreH(1:ic) = RejectMoreH(index_step(1:ic))
    Kacc(1:ic)        = Kacc(index_step(1:ic))
    Krej(1:ic)        = Krej(index_step(1:ic))

    return
  end subroutine kco_compress

  subroutine kco_finalize ( f, IERR, Kacc, Krej)
    implicit    none
    real(kind=dp),intent(INOUT),dimension(:,:)   :: f 
    integer,      intent(INOUT),dimension(:)     :: IERR, Kacc, Krej

    kpoints = kpoints_save

    f(1:kpoints,:)  = f_done(1:kpoints,:)   ! mz_pj_20070903 (1:kpoints) added
    ierr(1:kpoints) = ierr_done(1:kpoints)  ! mz_pj_20070903 ...
    Kacc(1:kpoints) = Kacc_done(1:kpoints)  ! mz_pj_20070903 ...
    Krej(1:kpoints) = Krej_done(1:kpoints)  ! mz_pj_20070903 ...

    deallocate(index_org)
    deallocate(index_step)
    deallocate(cell_done)
    deallocate(done_check)

    deallocate(f_done)
    deallocate(ierr_done)
    deallocate(Kacc_done)
    deallocate(Krej_done)

    return
  end subroutine kco_finalize

end module messy_main_tools_kp4_compress
