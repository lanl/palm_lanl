! NOTE: This routine is a vectorised version of the routine in
!       rosenbrock_mz.f90 and can only be used with KP4.
!
SUBROUTINE Rosenbrock(Y,Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!    Solves the system y'=F(t,y) using a Rosenbrock method defined by:
!
!     G = 1/(H*gamma(1)) - Jac(t0,Y0)
!     T_i = t0 + Alpha(i)*H
!     Y_i = Y0 + \sum_{j=1}^{i-1} A(i,j)*K_j
!     G * K_i = Fun( T_i, Y_i ) + \sum_{j=1}^S C(i,j)/H * K_j +
!         gamma(i)*dF/dT(t0, Y0)
!     Y1 = Y0 + \sum_{j=1}^S M(j)*K_j
!
!    For details on Rosenbrock methods and their implementation consult:
!      E. Hairer and G. Wanner
!      "Solving ODEs II. Stiff and differential-algebraic problems".
!      Springer series in computational mathematics, Springer-Verlag, 1996.
!    The codes contained in the book inspired this implementation.
!
!    (C)  Adrian Sandu, August 2004
!    Virginia Polytechnic Institute and State University
!    Contact: sandu@cs.vt.edu
!    This implementation is part of KPP - the Kinetic PreProcessor
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>   INPUT ARGUMENTS:
!
!-     Y(NVAR)    = vector of initial conditions (at T=Tstart)
!-    [Tstart,Tend]  = time range of integration
!     (if Tstart>Tend the integration is performed backwards in time)
!-    RelTol, AbsTol = user precribed accuracy
!- SUBROUTINE  Fun( T, Y, Ydot ) = ODE function,
!                       returns Ydot = Y' = F(T,Y)
!- SUBROUTINE  Jac( T, Y, Jcb ) = Jacobian of the ODE function,
!                       returns Jcb = dFun/dY
!-    ICNTRL(1:20)    = integer inputs parameters
!-    RCNTRL(1:20)    = real inputs parameters
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     OUTPUT ARGUMENTS:
!
!-    Y(NVAR)    -> vector of final states (at T->Tend)
!-    ISTATUS(1:20)   -> integer output parameters
!-    RSTATUS(1:20)   -> real output parameters
!-    IERR            -> job status upon return
!                        success (positive value) or
!                        failure (negative value)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     INPUT PARAMETERS:
!
!    Note: For input parameters equal to zero the default values of the
!       corresponding variables are used.
!
!    ICNTRL(1) = 1: F = F(y)   Independent of T (AUTONOMOUS)
!              = 0: F = F(t,y) Depends on T (NON-AUTONOMOUS)
!
!    ICNTRL(2) = 0: AbsTol, RelTol are NVAR-dimensional vectors
!              = 1: AbsTol, RelTol are scalars
!
!    ICNTRL(3)  -> selection of a particular Rosenbrock method
!        = 0 :  default method is Rodas3
!        = 1 :  method is  Ros2
!        = 2 :  method is  Ros3
!        = 3 :  method is  Ros4
!        = 4 :  method is  Rodas3
!        = 5:   method is  Rodas4
!
!    ICNTRL(4)  -> maximum number of integration steps
!        For ICNTRL(4)=0) the default value of 100000 is used
!
!    RCNTRL(1)  -> Hmin, lower bound for the integration step size
!          It is strongly recommended to keep Hmin = ZERO
!    RCNTRL(2)  -> Hmax, upper bound for the integration step size
!    RCNTRL(3)  -> Hstart, starting value for the integration step size
!
!    RCNTRL(4)  -> FacMin, lower bound on step decrease factor (default=0.2)
!    RCNTRL(5)  -> FacMax, upper bound on step increase factor (default=6)
!    RCNTRL(6)  -> FacRej, step decrease factor after multiple rejections
!            (default=0.1)
!    RCNTRL(7)  -> FacSafe, by which the new step is slightly smaller
!         than the predicted value  (default=0.9)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     OUTPUT PARAMETERS:
!
!    Note: each call to Rosenbrock adds the current no. of fcn calls
!      to previous value of ISTATUS(1), and similar for the other params.
!      Set ISTATUS(1:20) = 0 before call to avoid this accumulation.
!
!    ISTATUS(1) = No. of function calls
!    ISTATUS(2) = No. of jacobian calls
!    ISTATUS(3) = No. of steps
!    ISTATUS(4) = No. of accepted steps
!    ISTATUS(5) = No. of rejected steps (except at the beginning)
!    ISTATUS(6) = No. of LU decompositions
!    ISTATUS(7) = No. of forward/backward substitutions
!    ISTATUS(8) = No. of singular matrix decompositions
!
!    RSTATUS(1)  -> Texit, the time corresponding to the
!                   computed Y upon return
!    RSTATUS(2)  -> Hexit, last accepted step before exit
!    For multiple restarts, use Hexit as Hstart in the following run
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~>  Arguments
   REAL(kind=dp), INTENT(INOUT) :: Y(:,:)
   REAL(kind=dp), INTENT(IN)   :: Tstart,Tend
   REAL(kind=dp), INTENT(IN)   :: AbsTol(NVAR),RelTol(NVAR)
   INTEGER, INTENT(IN)    :: ICNTRL(20)
   REAL(kind=dp), INTENT(IN)   :: RCNTRL(20)
   INTEGER, INTENT(INOUT) :: ISTATUS(20)
   REAL(kind=dp), INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR
!~~~>  The method parameters
   INTEGER, PARAMETER :: Smax = 6
   INTEGER  :: Method, ros_S
   REAL(kind=dp), DIMENSION(Smax) :: ros_M, ros_E, ros_Alpha, ros_Gamma
   REAL(kind=dp), DIMENSION(Smax*(Smax-1)/2) :: ros_A, ros_C
   REAL(kind=dp) :: ros_ELO
   LOGICAL, DIMENSION(Smax) :: ros_NewF
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Local variables
   REAL(kind=dp) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   REAL(kind=dp) :: Hmin, Hmax, Hstart, Hexit
   REAL(kind=dp) :: Texit
   INTEGER :: i, UplimTol, Max_no_steps
   LOGICAL :: Autonomous, VectorTol
!~~~>   Parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp

!~~~>  Initialize statistics
   Nfun = ISTATUS(ifun)
   Njac = ISTATUS(ijac)
   Nstp = ISTATUS(istp)
   Nacc = ISTATUS(iacc)
   Nrej = ISTATUS(irej)
   Ndec = ISTATUS(idec)
   Nsol = ISTATUS(isol)
   Nsng = ISTATUS(isng)

!~~~>  Autonomous or time dependent ODE. Default is time dependent.
   Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:NVAR) and RelTol(1:NVAR)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
         UplimTol  = NVAR
   ELSE
      VectorTol = .FALSE.
         UplimTol  = 1
   END IF

!~~~>  The particular Rosenbrock method chosen
   IF (ICNTRL(3) == 0) THEN
      Method = 4
   ELSEIF ( (ICNTRL(3) >= 1).AND.(ICNTRL(3) <= 5) ) THEN
      Method = ICNTRL(3)
   ELSE
      PRINT * , 'User-selected Rosenbrock method: ICNTRL(3)=', Method
      !CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      IERRV(:) = -2
      RETURN
   END IF

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000 ! 200
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      !CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      IERRV(:) = -1
      RETURN
   END IF

!~~~>  Unit roundoff (1+Roundoff>1)
   Roundoff = WLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      !CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      IERRV(:) = -3
      RETURN
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      !CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      IERRV(:) = -3
      RETURN
   END IF
!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      !CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      IERRV(:) = -3
      RETURN
   END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hexit < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      !CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      IERRV(:) = -4
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      !CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      IERRV(:) = -4
      RETURN
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      !CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      IERRV(:) = -4
      RETURN
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      !CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      IERRV(:) = -4
      RETURN
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        !CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        IERRV(:) = -5
        RETURN
      END IF
    END DO


!~~~>   Initialize the particular Rosenbrock method
   SELECT CASE (Method)
     CASE (1)
       CALL Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=', Method
       !CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       IERRV(:) = -2
       RETURN
   END SELECT

!~~~>  CALL Rosenbrock method
   CALL ros_Integrator(Y,Tstart,Tend,Texit,      &
        AbsTol, RelTol,                          &
!  Rosenbrock method coefficients
        ros_S, ros_M, ros_E, ros_A, ros_C,       &
        ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, &
!  Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart, Hexit,     &
        FacMin, FacMax, FacRej, FacSafe,         &
!  Error indicator
        IERR)


!~~~>  Collect run statistics
   ISTATUS(ifun) = Nfun
   ISTATUS(ijac) = Njac
   ISTATUS(istp) = Nstp
   ISTATUS(iacc) = Nacc
   ISTATUS(irej) = Nrej
   ISTATUS(idec) = Ndec
   ISTATUS(isol) = Nsol
   ISTATUS(isng) = Nsng
!~~~> Last T and H
   RSTATUS(itexit) = Texit
   RSTATUS(ihexit) = Hexit

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS !  SUBROUTINES internal to Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   REAL(kind=dp), INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(INOUT) :: IERR

   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:'
   IF ((Code>=-8).AND.(Code<=-1)) THEN
     PRINT *, IERR_NAMES(Code)
   ELSE
     PRINT *, 'Unknown Error code: ', Code
   ENDIF

   PRINT *, "T=", T, "and H=", H

 END SUBROUTINE ros_ErrorMsg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_Integrator (Y, Tstart, Tend, Tout,     &
        AbsTol, RelTol,                          &
!~~~> Rosenbrock method coefficients
        ros_S, ros_M, ros_E, ros_A, ros_C,       &
        ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart, Hexit,     &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method
!      defined by ros_S (no of stages)
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE messy_main_tools_kp4_compress, ONLY: kco_initialize, kco_compress, &
                                           kco_finalize, cell_done
  IMPLICIT NONE

!~~~> Input: the initial condition at Tstart; Output: the solution at T
   REAL(kind=dp), INTENT(INOUT) :: Y(:,:)
!~~~> Input: integration interval
   REAL(kind=dp), INTENT(IN) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   REAL(kind=dp), INTENT(OUT) ::  Tout
!~~~> Input: tolerances
   REAL(kind=dp), INTENT(IN) ::  AbsTol(NVAR), RelTol(NVAR)
!~~~> Input: The Rosenbrock method parameters
   INTEGER, INTENT(IN) ::  ros_S
   REAL(kind=dp), INTENT(IN) :: ros_M(ros_S), ros_E(ros_S),  &
       ros_Alpha(ros_S), ros_A(ros_S*(ros_S-1)/2), &
       ros_Gamma(ros_S), ros_C(ros_S*(ros_S-1)/2), ros_ELO
   LOGICAL, INTENT(IN) :: ros_NewF(ros_S)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   REAL(kind=dp), INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   REAL(kind=dp), INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
!~~~> Output: last accepted step
   REAL(kind=dp), INTENT(OUT) :: Hexit
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables
   REAL(kind=dp) :: Ynew(VL_DIM,NVAR), Fcn0(VL_DIM,NVAR), Fcn(VL_DIM,NVAR)
   REAL(kind=dp) :: K(VL_DIM,NVAR*ros_S)
   REAL(kind=dp) :: Jac0(VL_DIM,LU_NONZERO), Ghimj(VL_DIM,LU_NONZERO)
   REAL(kind=dp) :: Yerr(VL_DIM,NVAR)

!~~~>  Local parameters
   INTEGER :: ioffset, i, j, istage
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp
!~~~>  Locally called functions
!    REAL(kind=dp) WLAMCH
!    EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   real(kind=dp)                  :: deltat
!kk   logical,parameter              :: fixed_step=.true.
   
   integer,save             :: icount=0

!  vector compression algorithm

   REAL(kind=dp), dimension(VL_dim)      :: T
   REAL(kind=dp), dimension(VL_dim)      :: H, Hnew, HC, HC2, HG, Fac, Tau, Err
   LOGICAL, dimension(VL_dim)            :: RejectLastH, RejectMoreH
   LOGICAL, dimension(VL_dim)            :: Accept_step
   integer                               :: vl_save

   REAL(kind=dp), dimension(size(FIX,1),size(FIX,2))      :: FIX_lo

!~~~>  Initial preparations
   T(1:VL)  = Tstart
   Hexit    = 0.0_dp
   H(1:VL)  = MIN(Hstart,Hmax)
   Tout     = Tend
   icount   = icount+1
   where (ABS(H(1:VL)) <= 10.0_dp*Roundoff) H(1:VL) = DeltaMin

   deltat=Tend-Tstart

   RejectLastH(1:VL)=.FALSE.
   RejectMoreH(1:VL)=.FALSE.

!~~~> Time loop begins below

   Nstp       = 0
   IERRV(1:VL) = 1  ! successful integration
   vl_save          = VL
   FIX_lo(1:VL,:)   = FIX(1:VL,:)
   Kacc(1:VL)       = 0
   Krej(1:VL)       = 0

!  intialize compression algorithm

   IF (.not. l_fixed_step) then
      call kco_initialize (VL, y, IERRV, Kacc, Krej)
   end if

TimeLoop: DO 

   IF ( Nstp > Max_no_steps ) THEN  ! Too many steps
      !CALL ros_ErrorMsg(-6,T(1),H(1),IERR)
      IERRV(1:VL) = -6
!      call abort()                                   !Vorlaeufig
      RETURN
   END IF
   
   Where ( (T(1:VL)+0.1_dp*H(1:VL)) == T(1:VL))      ! Step size too small
      IERRV(1:VL) = -7
   endwhere

!!$   if( any(error_flag))  then
!!$      CALL ros_ErrorMsg(-7,T(1),H(1),IERR)
!!$      call abort()                                   !Vorlaeufig
!!$      RETURN
!!$   end if

!~~~>  Limit H if necessary to avoid going beyond Tend
!kk   Hexit = H
   H(1:VL) = MIN(H(1:VL),ABS(Tend-T(1:VL)))

!~~~>   Compute the function at current time
   CALL Fun( Y, FIX_lo, RCONST, Fcn0)


!~~~>   Compute the Jacobian at current time
   CALL Jac_SP( Y, FIX_lo, RCONST, Jac0)


!~~~>  Repeat step calculation until current step accepted

   CALL ros_PrepareMatrix(H,ros_Gamma(1), Jac0,Ghimj)

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S

      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+NVAR)
       ioffset = NVAR*(istage-1)

      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         Fcn(1:VL,1:NVAR) = Fcn0(1:VL,1:NVAR) 
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         select case (istage-1)
          case (1)
           Ynew(1:VL,1:NVAR) = Y(1:VL,1:NVAR)
           Ynew(1:VL,1:NVAR) = Ynew(1:VL,1:NVAR) +  K(1:VL,NVAR*(1- 1)+ 1:NVAR*1) &
                        *ros_A((istage-1)*(istage-2)/2+ 1)
          case (2)
           Ynew(1:VL,1:NVAR) = Y(1:VL,1:NVAR)
           Ynew(1:VL,1:NVAR) = Ynew(1:VL,1:NVAR) +  K(1:VL,NVAR*(2- 1)+ 1:NVAR*2) &
                        *ros_A((istage-1)*(istage-2)/2+ 2)
           Ynew(1:VL,1:NVAR) = Ynew(1:VL,1:NVAR) +  K(1:VL,NVAR*(2- 1)+ 1:NVAR*2) &
                        *ros_A((istage-1)*(istage-2)/2+ 2)
          case DEFAULT
           Ynew(1:VL,1:NVAR) = Y(1:VL,1:NVAR)
           DO j = 1, istage-1
             Ynew(1:VL,1:NVAR) = Ynew(1:VL,1:NVAR) +  K(1:VL,NVAR*(j- 1)+ 1:NVAR*j) &
                          *ros_A((istage-1)*(istage-2)/2+ j)
           END DO
         end select

         CALL Fun( Ynew, FIX_lo, RCONST, Fcn)
       END IF ! if istage == 1 elseif ros_NewF(istage)

       select case (istage-1)
        case (1)
         K(1:VL,ioffset+ 1:ioffset+ NVAR) = Fcn(1:VL,1:NVAR)
         j = 1
         HC(1:VL) = ros_C((istage-1)*(istage-2)/2+ j)/(H(1:VL))
         do i=1,NVAR
           K(1:VL,ioffset+ i) = K(1:VL,ioffset+ i) +  K(1:VL,NVAR*(j- 1)+ i) *HC(1:VL)
         end do
        case (2)
         K(1:VL,ioffset+ 1:ioffset+ NVAR) = Fcn(1:VL,1:NVAR)
         HC(1:VL)  = ros_C((istage-1)*(istage-2)/2+ 1)/(H(1:VL))
         HC2(1:VL) = ros_C((istage-1)*(istage-2)/2+ 2)/(H(1:VL))
         do i=1,NVAR
           K(1:VL,ioffset+ i) = K(1:VL,ioffset+ i) +  K(1:VL,NVAR*(1- 1)+ i) *HC(1:VL)
           K(1:VL,ioffset+ i) = K(1:VL,ioffset+ i) +  K(1:VL,NVAR*(2- 1)+ i) *HC2(1:VL)
         end do
        case DEFAULT
         K(1:VL,ioffset+ 1:ioffset+ NVAR) = Fcn(1:VL,1:NVAR)
         DO j = 1, istage-1
           HC(1:VL) = ros_C((istage-1)*(istage-2)/2+ j)/(H(1:VL))
           do i=1,NVAR
             K(1:VL,ioffset+ i) = K(1:VL,ioffset+ i) +  K(1:VL,NVAR*(j- 1)+ i) *HC(1:VL)
           end do
         END DO
       end select

       CALL ros_Solve(Ghimj, K(:,ioffset+1:ioffset+NVAR))

   END DO Stage


   if(ros_S == 3)   then
     YNEW(1:VL,1:NVAR) = Y(1:VL,1:NVAR)
     Yerr(1:VL,1:NVAR) = ZERO

!~~~>  Compute the new solution
     Ynew(1:VL,1:NVAR) = Ynew(1:VL,1:NVAR) + K(1:VL,NVAR*(1-1)+1:NVAR*1) * ros_m(1)
     Ynew(1:VL,1:NVAR) = Ynew(1:VL,1:NVAR) + K(1:VL,NVAR*(2-1)+1:NVAR*2) * ros_m(2)
     Ynew(1:VL,1:NVAR) = Ynew(1:VL,1:NVAR) + K(1:VL,NVAR*(3-1)+1:NVAR*3) * ros_m(3)

!~~~>  Compute the error estimation

     Yerr(1:VL,1:NVAR) = Yerr(1:VL,1:NVAR) + K(1:VL,NVAR*(1-1)+1:NVAR*1) * ros_E(1)
     Yerr(1:VL,1:NVAR) = Yerr(1:VL,1:NVAR) + K(1:VL,NVAR*(2-1)+1:NVAR*2) * ros_E(2)
     Yerr(1:VL,1:NVAR) = Yerr(1:VL,1:NVAR) + K(1:VL,NVAR*(3-1)+1:NVAR*3) * ros_E(3)
   else
     YNEW(1:VL,1:NVAR) = Y(1:VL,1:NVAR)
     Yerr(1:VL,1:NVAR) = ZERO
     DO j=1,ros_S
!~~~>  Compute the new solution
         Ynew(1:VL,1:NVAR) = Ynew(1:VL,1:NVAR) + &
                        K(1:VL,NVAR*(j-1)+1:NVAR*j) * ros_m(j)

!~~~>  Compute the error estimation

       Yerr(1:VL,1:NVAR) = Yerr(1:VL,1:NVAR) + &
                    K(1:VL,NVAR*(j-1)+1:NVAR*j) * ros_E(j)
     END DO
   end if


!~~~>  Check the error magnitude and adjust step size
   Nstp = Nstp+1

   if(l_fixed_step)   then
     Err  = 0.0
     H    = t_steps(Nstp)*deltat
     Hnew = t_steps(Nstp)*deltat
   else
     Err(1:VL) = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )
!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
     where(Err(1:VL) >= 10.0_dp*Roundoff)
       Fac(1:VL)  = MIN(FacMax,MAX(FacMin,FacSafe/Err(1:VL)**(ONE/ros_ELO)))
     elsewhere
       Fac(1:VL)  = 1000.                        !Force large time step
     endwhere
     Hnew(1:VL) = H(1:VL)*Fac(1:VL)
   end if

   Accept_step(1:VL) = (Err(1:VL) <= ONE).OR.(H(1:VL) <= Hmin)
   where (Accept_step(1:VL))
      Kacc(1:VL) = Kacc(1:VL)+1
      T(1:VL) = T(1:VL) + H(1:VL)
      Hnew(1:VL) = MAX(Hmin,MIN(Hnew(1:VL),Hmax))
   endwhere
   do i=1,NVAR
     where (Accept_step(1:VL))
       Y(1:VL,i) = MAX(Ynew(1:VL,i),ZERO)
     endwhere
   end do
   where (Accept_step(1:VL) .and. RejectLastH(1:VL))
      Hnew(1:VL) = MIN(Hnew(1:VL),H(1:VL))
   endwhere
   where (Accept_step(1:VL))
      RejectLastH(1:VL) = .FALSE.
      RejectMoreH(1:VL) = .FALSE.
      H(1:VL) = Hnew(1:VL)
   endwhere
!kk   ELSE           !~~~> Reject step
   where (.NOT. Accept_step(1:VL) .and. RejectMoreH(1:VL))
      Hnew(1:VL) = H(1:VL)*FacRej
   endwhere
   where (.NOT. Accept_step(1:VL))
      RejectMoreH(1:VL) = RejectLastH(1:VL)
      RejectLastH(1:VL) = .TRUE.
      H(1:VL) = Hnew(1:VL)
   endwhere
   where ( .NOT. Accept_step(1:VL) .and. Kacc(1:VL) >= 1)
      Krej(1:VL) = Krej(1:VL)+1
   endwhere


   if(l_fixed_step)   then
      if(Nstp >= Nfsteps) EXIT
   else
     cell_done(1:VL) = ( .NOT. (T(1:VL)-Tend)+Roundoff <= ZERO  )  .OR. &
          (IERRV(1:VL) < 0) ! mz_pj_20070516
     call kco_compress (VL, T, H, Hnew, IERRV, y, RCONST, FIX_lo,         &
                                     RejectLastH, RejectMoreH, Kacc, Krej)
     if(VL == 0)   EXIT
   end if

!kk   write(0,*) '  aaa ',icount,Nstp,count(cell_done),count(Accept_step),VL,minval(T)
   END DO TimeLoop

!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

   IF (.not. l_fixed_step) then
     call kco_finalize (y, IERRV, Kacc, Krej)
   end if

   VL = VL_save

!kk   write(0,*) 'bbb ',icount,Nstp,maxval(Kacc),maxval(Krej)
!kk   if(icount >= 10)  call abort()

  END SUBROUTINE ros_Integrator


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, &
               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

! Input arguments
   REAL(kind=dp), INTENT(IN) :: Y(:,:), Ynew(:,:), &
          Yerr(:,:), AbsTol(NVAR), RelTol(NVAR)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables
   REAL(kind=dp), dimension(VL)       :: Err, Scale, Ymax
   INTEGER  :: i
   REAL(kind=dp), PARAMETER           :: ZERO = 0.0_dp
   REAL(kind=dp), dimension(VL)       :: ros_ErrorNorm

   Err = ZERO
   DO i=1,NVAR
     Ymax = MAX( ABS(Y(1:VL,i)),ABS(Ynew(1:VL,i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+ RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+ RelTol(1)*Ymax
     END IF
     Err = Err+ (Yerr(1:VL,i)/Scale)**2
   END DO
   Err  = SQRT(Err/NVAR)

   ros_ErrorNorm = Err

  END FUNCTION ros_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE ros_PrepareMatrix ( H, gam, Jac0, Ghimj )

! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   IMPLICIT NONE

!~~~> Input arguments
   REAL(kind=dp), INTENT(IN) ::  Jac0(:,:)
   REAL(kind=dp), INTENT(IN) ::  gam
!~~~> Output arguments
   REAL(kind=dp), INTENT(OUT) :: Ghimj(:,:)
!~~~> Inout arguments
   REAL(kind=dp), INTENT(INOUT) :: H(:)   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: i, ising, Nconsecutive
   REAL(kind=dp), dimension(size(H)) :: ghinv
   REAL(kind=dp), PARAMETER :: ONE  = 1.0_dp, HALF = 0.5_dp

   Nconsecutive = 0


!~~~>    Construct Ghimj = 1/(H*gam) - Jac0
   Ghimj(1:VL,1:LU_NONZERO) = -ONE * JAC0(1:VL,1:LU_NONZERO)
   ghinv(1:VL) = ONE/(H(1:VL)*gam)
   DO i=1,NVAR
     Ghimj(1:VL,LU_DIAG(i)) = Ghimj(1:VL,LU_DIAG(i))+ghinv(1:VL)
   END DO

!~~~>    Compute LU decomposition
!CDIR noiexpand
   CALL KppDecomp( Ghimj, ising )

   Ndec = Ndec + VL

  END SUBROUTINE ros_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Solve( A, b )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution (using pre-computed LU decomposition)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Input variables
   REAL(kind=dp), INTENT(IN)         :: A(:,:)
!~~~> InOut variables
   REAL(kind=dp), INTENT(INOUT)      :: b(:,:)

   CALL KppSolve( A, b )

   Nsol = Nsol+VL

  END SUBROUTINE ros_Solve



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
            ros_Gamma,ros_NewF,ros_ELO,ros_Name)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

   INTEGER, PARAMETER :: S=2
   INTEGER, INTENT(OUT) ::  ros_S
   REAL(kind=dp), DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   REAL(kind=dp), DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   REAL(kind=dp), INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name
   DOUBLE PRECISION g

    g = 1.0_dp + 1.0_dp/SQRT(2.0_dp)

!~~~> Name of the method
    ros_Name = 'ROS-2'
!~~~> Number of stages
    ros_S = S

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = (1.0_dp)/g
    ros_C(1) = (-2.0_dp)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.0_dp)/(2.0_dp*g)
    ros_M(2)= (1.0_dp)/(2.0_dp*g)
! E_i = Coefficients for error estimator
    ros_E(1) = 1.0_dp/(2.0_dp*g)
    ros_E(2) = 1.0_dp/(2.0_dp*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0_dp
    ros_Alpha(2) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = g
    ros_Gamma(2) =-g

 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
           ros_Gamma,ros_NewF,ros_ELO,ros_Name)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

   INTEGER, PARAMETER :: S=3
   INTEGER, INTENT(OUT) ::  ros_S
   REAL(kind=dp), DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   REAL(kind=dp), DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   REAL(kind=dp), INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name

!~~~> Name of the method
   ros_Name = 'ROS-3'
!~~~> Number of stages
   ros_S = S

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1)= 1.0_dp
   ros_A(2)= 1.0_dp
   ros_A(3)= 0.0_dp

   ros_C(1) = -0.10156171083877702091975600115545E+01_dp
   ros_C(2) =  0.40759956452537699824805835358067E+01_dp
   ros_C(3) =  0.92076794298330791242156818474003E+01_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1E+01_dp
   ros_M(2) =  0.61697947043828245592553615689730E+01_dp
   ros_M(3) = -0.42772256543218573326238373806514E+00_dp
! E_i = Coefficients for error estimator
   ros_E(1) =  0.5E+00_dp
   ros_E(2) = -0.29079558716805469821718236208017E+01_dp
   ros_E(3) =  0.22354069897811569627360909276199E+00_dp
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0E+00_dp
   ros_Alpha(2)= 0.43586652150845899941601945119356E+00_dp
   ros_Alpha(3)= 0.43586652150845899941601945119356E+00_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1)= 0.43586652150845899941601945119356E+00_dp
   ros_Gamma(2)= 0.24291996454816804366592249683314E+00_dp
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01_dp

  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
           ros_Gamma,ros_NewF,ros_ELO,ros_Name)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

   INTEGER, PARAMETER :: S=4
   INTEGER, INTENT(OUT) ::  ros_S
   REAL(kind=dp), DIMENSION(4), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   REAL(kind=dp), DIMENSION(6), INTENT(OUT) :: ros_A, ros_C
   REAL(kind=dp), INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(4), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name
   DOUBLE PRECISION g

!~~~> Name of the method
   ros_Name = 'ROS-4'
!~~~> Number of stages
   ros_S = S

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.2000000000000000E+01_dp
   ros_A(2) = 0.1867943637803922E+01_dp
   ros_A(3) = 0.2344449711399156E+00_dp
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0_dp

   ros_C(1) =-0.7137615036412310E+01_dp
   ros_C(2) = 0.2580708087951457E+01_dp
   ros_C(3) = 0.6515950076447975E+00_dp
   ros_C(4) =-0.2137148994382534E+01_dp
   ros_C(5) =-0.3214669691237626E+00_dp
   ros_C(6) =-0.6949742501781779E+00_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735E+01_dp
   ros_M(2) = 0.2870493262186792E+00_dp
   ros_M(3) = 0.4353179431840180E+00_dp
   ros_M(4) = 0.1093502252409163E+01_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) =-0.2815431932141155E+00_dp
   ros_E(2) =-0.7276199124938920E-01_dp
   ros_E(3) =-0.1082196201495311E+00_dp
   ros_E(4) =-0.1093502252409163E+01_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.1145640000000000E+01_dp
   ros_Alpha(3) = 0.6552168638155900E+00_dp
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5728200000000000E+00_dp
   ros_Gamma(2) =-0.1769193891319233E+01_dp
   ros_Gamma(3) = 0.7592633437920482E+00_dp
   ros_Gamma(4) =-0.1049021087100450E+00_dp

  END SUBROUTINE Ros4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
            ros_Gamma,ros_NewF,ros_ELO,ros_Name)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

   INTEGER, PARAMETER :: S=4
   INTEGER, INTENT(OUT) ::  ros_S
   REAL(kind=dp), DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   REAL(kind=dp), DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   REAL(kind=dp), INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name
   DOUBLE PRECISION g

!~~~> Name of the method
   ros_Name = 'RODAS-3'
!~~~> Number of stages
   ros_S = S

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.0E+00_dp
   ros_A(2) = 2.0E+00_dp
   ros_A(3) = 0.0E+00_dp
   ros_A(4) = 2.0E+00_dp
   ros_A(5) = 0.0E+00_dp
   ros_A(6) = 1.0E+00_dp

   ros_C(1) = 4.0E+00_dp
   ros_C(2) = 1.0E+00_dp
   ros_C(3) =-1.0E+00_dp
   ros_C(4) = 1.0E+00_dp
   ros_C(5) =-1.0E+00_dp
   ros_C(6) =-(8.0E+00_dp/3.0E+00_dp)

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0E+00_dp
   ros_M(2) = 0.0E+00_dp
   ros_M(3) = 1.0E+00_dp
   ros_M(4) = 1.0E+00_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) = 0.0E+00_dp
   ros_E(2) = 0.0E+00_dp
   ros_E(3) = 0.0E+00_dp
   ros_E(4) = 1.0E+00_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0E+00_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0E+00_dp
   ros_Alpha(2) = 0.0E+00_dp
   ros_Alpha(3) = 1.0E+00_dp
   ros_Alpha(4) = 1.0E+00_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5E+00_dp
   ros_Gamma(2) = 1.5E+00_dp
   ros_Gamma(3) = 0.0E+00_dp
   ros_Gamma(4) = 0.0E+00_dp

  END SUBROUTINE Rodas3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
             ros_Gamma,ros_NewF,ros_ELO,ros_Name)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

   INTEGER, PARAMETER :: S=6
   INTEGER, INTENT(OUT) ::  ros_S
   REAL(kind=dp), DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   REAL(kind=dp), DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   REAL(kind=dp), INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name
   DOUBLE PRECISION g

!~~~> Name of the method
    ros_Name = 'RODAS-4'
!~~~> Number of stages
    ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000_dp
    ros_Alpha(2) = 0.386_dp
    ros_Alpha(3) = 0.210_dp
    ros_Alpha(4) = 0.630_dp
    ros_Alpha(5) = 1.000_dp
    ros_Alpha(6) = 1.000_dp

!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = 0.2500000000000000E+00_dp
    ros_Gamma(2) =-0.1043000000000000E+00_dp
    ros_Gamma(3) = 0.1035000000000000E+00_dp
    ros_Gamma(4) =-0.3620000000000023E-01_dp
    ros_Gamma(5) = 0.0_dp
    ros_Gamma(6) = 0.0_dp

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = 0.1544000000000000E+01_dp
    ros_A(2) = 0.9466785280815826E+00_dp
    ros_A(3) = 0.2557011698983284E+00_dp
    ros_A(4) = 0.3314825187068521E+01_dp
    ros_A(5) = 0.2896124015972201E+01_dp
    ros_A(6) = 0.9986419139977817E+00_dp
    ros_A(7) = 0.1221224509226641E+01_dp
    ros_A(8) = 0.6019134481288629E+01_dp
    ros_A(9) = 0.1253708332932087E+02_dp
    ros_A(10) =-0.6878860361058950E+00_dp
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0E+00_dp

    ros_C(1) =-0.5668800000000000E+01_dp
    ros_C(2) =-0.2430093356833875E+01_dp
    ros_C(3) =-0.2063599157091915E+00_dp
    ros_C(4) =-0.1073529058151375E+00_dp
    ros_C(5) =-0.9594562251023355E+01_dp
    ros_C(6) =-0.2047028614809616E+02_dp
    ros_C(7) = 0.7496443313967647E+01_dp
    ros_C(8) =-0.1024680431464352E+02_dp
    ros_C(9) =-0.3399990352819905E+02_dp
    ros_C(10) = 0.1170890893206160E+02_dp
    ros_C(11) = 0.8083246795921522E+01_dp
    ros_C(12) =-0.7981132988064893E+01_dp
    ros_C(13) =-0.3152159432874371E+02_dp
    ros_C(14) = 0.1631930543123136E+02_dp
    ros_C(15) =-0.6058818238834054E+01_dp

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0E+00_dp
    ros_M(6) = 1.0E+00_dp

!~~~> E_i  = Coefficients for error estimator
    ros_E(1) = 0.0E+00_dp
    ros_E(2) = 0.0E+00_dp
    ros_E(3) = 0.0E+00_dp
    ros_E(4) = 0.0E+00_dp
    ros_E(5) = 0.0E+00_dp
    ros_E(6) = 1.0E+00_dp

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 4.0_dp

  END SUBROUTINE Rodas4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   End of the set of internal Rosenbrock subroutines
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
END SUBROUTINE Rosenbrock
