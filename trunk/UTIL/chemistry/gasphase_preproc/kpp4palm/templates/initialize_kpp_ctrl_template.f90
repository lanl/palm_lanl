MODULE initialize_kpp_ctrl_template

  ! HEADER MODULE initialize_kpp_ctrl_template

  ! NOTES:
  ! - L_VECTOR is automatically defined by kp4
  ! - VL_DIM is automatically defined by kp4
  ! - I_LU_DI is automatically defined by kp4
  ! - WANTED is automatically defined by xmecca
  ! - ICNTRL RCNTRL are automatically defined by kpp
  ! - "USE messy_main_tools" is in Module_header of messy_mecca_kpp.f90
  ! - SAVE will be automatically added by kp4

  IMPLICIT NONE
  !SAVE

  ! FOR FIXED TIME STEP CONTROL
  ! ... max. number of fixed time steps (sum must be 1)
  INTEGER, PARAMETER         :: NMAXFIXSTEPS = 50
  ! ... switch for fixed time stepping
  LOGICAL, PUBLIC            :: l_fixed_step = .FALSE.
  INTEGER, PUBLIC            :: nfsteps = 1
  ! ... number of kpp control parameters
  INTEGER, PARAMETER, PUBLIC :: NKPPCTRL = 20
  !
  INTEGER,  DIMENSION(NKPPCTRL), PUBLIC     :: icntrl = 0
  REAL(DP), DIMENSION(NKPPCTRL), PUBLIC     :: rcntrl = 0.0_dp
  REAL(DP), DIMENSION(NMAXFIXSTEPS), PUBLIC :: t_steps = 0.0_dp

  ! END HEADER MODULE initialize_kpp_ctrl_template

CONTAINS

SUBROUTINE initialize_kpp_ctrl(status, iou, modstr)

  IMPLICIT NONE

  ! I/O
  INTEGER,          INTENT(OUT) :: status
  INTEGER,          INTENT(IN)  :: iou     ! logical I/O unit
  CHARACTER(LEN=*), INTENT(IN)  :: modstr  ! read <modstr>.nml

  ! LOCAL
  REAL(DP) :: tsum
  INTEGER  :: i

  ! check fixed time steps
  tsum = 0.0_dp
  DO i=1, NMAXFIXSTEPS
     IF (t_steps(i) < TINY(0.0_DP)) EXIT
     tsum = tsum + t_steps(i)
  END DO

  nfsteps = i-1

  l_fixed_step = (nfsteps > 0) .AND. ( (tsum -1.0) < TINY(0.0_DP) )

  IF (L_VECTOR) THEN
     WRITE(*,*) ' MODE             : VECTOR (LENGTH=',VL_DIM,')'
  ELSE
     WRITE(*,*) ' MODE             : SCALAR'
  END IF
  !
  WRITE(*,*) ' DE-INDEXING MODE :',I_LU_DI
  !
  WRITE(*,*) ' ICNTRL           : ',icntrl
  WRITE(*,*) ' RCNTRL           : ',rcntrl
  !
  ! NOTE: THIS IS ONLY MEANINGFUL FOR VECTORIZED (kp4) ROSENBROCK-METHODS
  IF (L_VECTOR) THEN
     IF (l_fixed_step) THEN
        WRITE(*,*) ' TIME STEPS       : FIXED (',t_steps(1:nfsteps),')'
     ELSE
        WRITE(*,*) ' TIME STEPS       : AUTOMATIC'
     END IF
  ELSE
     WRITE(*,*) ' TIME STEPS       : AUTOMATIC '//&
          &'(t_steps (CTRL_KPP) ignored in SCALAR MODE)'
  END IF
  ! mz_pj_20070531-

  status = 0


END SUBROUTINE initialize_kpp_ctrl

SUBROUTINE error_output(C,ierr,PE)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ierr
  INTEGER, INTENT(IN) :: PE
  REAL(dp), DIMENSION(:),INTENT(IN) :: C

  write(6,*) 'ERROR in chem_gasphase_mod ',ierr,C(1)


END SUBROUTINE error_output

END MODULE initialize_kpp_ctrl_template
