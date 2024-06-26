! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! INTEGRATE - Integrator routine
!   Arguments :
!      TIN       - Start Time for Integration
!      TOUT      - End Time for Integration
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Rosenbrock - Implementation of several Rosenbrock methods:             !
!               * Ros2                                                    !
!               * Ros3                                                    !
!               * Ros4                                                    !
!               * Rodas3                                                  !
!               * Rodas4                                                  !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

MODULE ROS_f90_Integrator

  IMPLICIT NONE
  
!~~~>  Statistics on the work performed by the Rosenbrock method
  INTEGER, PARAMETER :: Nfun=1, Njac=2, Nstp=3, Nacc=4, &
                        Nrej=5, Ndec=6, Nsol=7, Nsng=8, &
                        Ntexit=1, Nhexit=2, Nhnew = 3,  &
                        Npassjac=1

CONTAINS

SUBROUTINE INTEGRATE_ROS( TIN, TOUT, N, NNZERO, VAR, RTOL, ATOL, FUN, JAC,&
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: N,NNZERO
   DOUBLE PRECISION, INTENT(INOUT) ::VAR(N)
   DOUBLE PRECISION, INTENT(IN) :: RTOL(N)
   DOUBLE PRECISION, INTENT(IN) :: ATOL(N)
   DOUBLE PRECISION, INTENT(IN) :: TIN  ! Start Time
   DOUBLE PRECISION, INTENT(IN) :: TOUT ! End Time
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U

   DOUBLE PRECISION :: RCNTRL(20), RSTATUS(20)
   INTEGER       :: ICNTRL(20), ISTATUS(20), IERR

   EXTERNAL :: FUN,JAC
   ! INTEGER, SAVE :: Ntotal = 0

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0

    !~~~> fine-tune the integrator:
   ICNTRL(1) = 0	! 0 - non-autonomous, 1 - autonomous
   ICNTRL(2) = 0	! 0 - vector tolerances, 1 - scalars

   ! If optional parameters are given, and if they are >0, 
   ! then they overwrite default settings. 
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF


   CALL Rosenbrock(N,NNZERO,VAR,TIN,TOUT,   &
         ATOL,RTOL, FUN, JAC,          &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)

   !~~~> Debug option: show no of steps
!    Ntotal = Ntotal + ISTATUS(Nstp)
!    PRINT*,'NSTEPS=',ISTATUS(Nstp),' (',Ntotal,')','  O3=', VAR(1)

   !STEPMIN = RSTATUS(Nhexit) seems no use
   ! if optional parameters are given for output they 
   ! are updated with the return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE INTEGRATE_ROS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE Rosenbrock(N,NNZERO,Y,Tstart,Tend, &
           AbsTol,RelTol, FUN, JAC,        &
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
!    Revised by Philipp Miehe and Adrian Sandu, May 2006                  
!    Revided by Hong Zhang and Adrian Sandu, Feb 2011
!    This implementation is part of FATODE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>   INPUT ARGUMENTS:
!
!-     Y(N)    = vector of initial conditions (at T=Tstart)
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
!-    Y(N)    -> vector of final states (at T->Tend)
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
!    ICNTRL(2) = 0: AbsTol, RelTol are N-dimensional vectors
!              = 1: AbsTol, RelTol are scalars
!
!    ICNTRL(3)  -> selection of a particular Rosenbrock method
!        = 0 :    Rodas3 (default)
!        = 1 :    Ros2
!        = 2 :    Ros3
!        = 3 :    Ros4
!        = 4 :    Rodas3
!        = 5 :    Rodas4
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
!                          (default=0.1)
!    RCNTRL(7)  -> FacSafe, by which the new step is slightly smaller
!         than the predicted value  (default=0.9)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!    OUTPUT ARGUMENTS:
!    -----------------
!
!    T           -> T value for which the solution has been computed
!                     (after successful return T=Tend).
!
!    Y(N)        -> Numerical solution at T
!
!    IDID        -> Reports on successfulness upon return:
!                    = 1 for success
!                    < 0 for error (value equals error code)
!
!    ISTATUS(1)  -> No. of function calls
!    ISTATUS(2)  -> No. of jacobian calls
!    ISTATUS(3)  -> No. of steps
!    ISTATUS(4)  -> No. of accepted steps
!    ISTATUS(5)  -> No. of rejected steps (except at very beginning)
!    ISTATUS(6)  -> No. of LU decompositions
!    ISTATUS(7)  -> No. of forward/backward substitutions
!    ISTATUS(8)  -> No. of singular matrix decompositions
!
!    RSTATUS(1)  -> Texit, the time corresponding to the
!                     computed Y upon return
!    RSTATUS(2)  -> Hexit, last accepted step before exit
!    RSTATUS(3)  -> Hnew, last predicted step (not yet taken)
!                   For multiple restarts, use Hnew as Hstart 
!                     in the subsequent run
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   USE LS_Solver_ROS
   IMPLICIT NONE

!~~~>  Arguments
   INTEGER,       INTENT(IN)    :: N,NNZERO
   DOUBLE PRECISION, INTENT(INOUT) :: Y(N)
   DOUBLE PRECISION, INTENT(IN)    :: Tstart,Tend
   DOUBLE PRECISION, INTENT(IN)    :: AbsTol(N),RelTol(N)
   INTEGER,       INTENT(IN)    :: ICNTRL(20)
   DOUBLE PRECISION, INTENT(IN)    :: RCNTRL(20)
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   DOUBLE PRECISION, INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR
   DOUBLE PRECISION, EXTERNAL :: DLAMCH
   EXTERNAL FUN,JAC
   EXTERNAL :: DCOPY, DAXPY, DGETRF, DGETRS, DSCAL
!~~~>  LS Solver data
   TYPE(LSdata) :: lssdata
!~~~>  Parameters of the Rosenbrock method, up to 6 stages
   INTEGER ::  ros_S, rosMethod
   INTEGER, PARAMETER :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5
   DOUBLE PRECISION :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                    ros_Alpha(6), ros_Gamma(6), ros_ELO
   LOGICAL :: ros_NewF(6)
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Local variables
   DOUBLE PRECISION :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   DOUBLE PRECISION :: Hmin, Hmax, Hstart
   DOUBLE PRECISION :: Texit
   INTEGER       :: i, UplimTol, Max_no_steps
   LOGICAL       :: Autonomous, VectorTol
!~~~>   Parameters
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.0, ONE  = 1.0, TOLLINI = 1D-20
   DOUBLE PRECISION, PARAMETER :: DeltaMin = 1.0E-5

   ros_A = ZERO
   ros_C = ZERO
   ros_M = ZERO
   ros_E = ZERO
   ros_Alpha = ZERO
   ros_Gamma = TOLLINI
!~~~>  Initialize statistics
   ISTATUS(1:8) = 0
   RSTATUS(1:3) = ZERO

!~~~>  Autonomous or time dependent ODE. Default is time dependent.
   Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = N
   ELSE
      VectorTol = .FALSE.
      UplimTol  = 1
   END IF

!~~~>   Initialize the particular Rosenbrock method selected
   SELECT CASE (ICNTRL(3))
     CASE (1)
       CALL Ros2
     CASE (2)
       CALL Ros3
     CASE (0,3)
       CALL Ros4
     CASE (4)
       CALL Rodas3
     CASE (5)
       CALL Rodas4
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=',ICNTRL(3) 
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 200000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Unit roundoff (1+Roundoff>1)
   Roundoff = DLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0*Roundoff) &
         .OR. (RelTol(i) >= 1.0) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO

   CALL LSS_Init(lssdata,N,NNZERO)
!~~~>  CALL Rosenbrock method
   CALL ros_Integrator(Y, Tstart, Tend, Texit,   &
        AbsTol, RelTol,                          &
!  Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,lssdata, &
!  Error indicator
        IERR)
   CALL LSS_Free(lssdata)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS !  SUBROUTINES internal to Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
   DOUBLE PRECISION, INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR
   
   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)    
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)    
      PRINT * , '--> Selected Rosenbrock method not implemented'
    CASE (-3)    
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)    
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)    
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT
   
   PRINT *, "T=", T, "and H=", H
     
 END SUBROUTINE ros_ErrorMsg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_Integrator (Y, Tstart, Tend, T,  &
        AbsTol, RelTol,                          &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe, lssdata,&
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method
!      defined by ros_S (no of stages)
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  USE LS_Solver_ROS, only: LSdata
  IMPLICIT NONE

!~~~> Input: the initial condition at Tstart; Output: the solution at T
   DOUBLE PRECISION, INTENT(INOUT) :: Y(N)
!~~~> Input: integration interval
   DOUBLE PRECISION, INTENT(IN) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   DOUBLE PRECISION, INTENT(OUT) ::  T
!~~~> Input: tolerances
   DOUBLE PRECISION, INTENT(IN) ::  AbsTol(N), RelTol(N)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   DOUBLE PRECISION, INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   DOUBLE PRECISION, INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   TYPE(LSdata) :: lssdata
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables
   DOUBLE PRECISION :: Ynew(N), Fcn0(N), Fcn(N)
   DOUBLE PRECISION :: K(N*ros_S), dFdT(N)
!   DOUBLE PRECISION :: Jac0(N,N), Ghimj(N,N)
   DOUBLE PRECISION :: H, Hnew, HC, HG, Fac, Tau
   DOUBLE PRECISION :: Err, Yerr(N)
   INTEGER :: Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular
!~~~>  Local parameters
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.0, ONE  = 1.0
   DOUBLE PRECISION, PARAMETER :: DeltaMin = 1.0E-5
!~~~>  Local data to compite Jacobian
   DOUBLE PRECISION :: DER(N), DER0(N)
   DOUBLE PRECISION, PARAMETER :: UROUND = 1d-19
   DOUBLE PRECISION :: YSAFE, DELY
   INTEGER :: I
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~>  Initial preparations
   T = Tstart
   RSTATUS(Nhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

!~~~> Time loop begins below

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff < ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff < ZERO) )

   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H = MIN(H,ABS(Tend-T))


!~~~>   Compute the function at current time
   CALL FUN(N,T,Y,Fcn0)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1


!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative (T, Roundoff, Y, &
                Fcn0, dFdT )
   END IF
 

!~~~>   Compute the Jacobian at current time
!   CALL JAC(T,Y,Jac0)
   if (Npassjac==0) then
     CALL LSS_Jac(N,lssdata,T,Y,JAC)
   else
      DO I = 1, N
         YSAFE=Y(I)
         CALL FUN(N,T,Y,DER0)
         DELY=DSQRT(UROUND*MAX(1.D-5,ABS(YSAFE)))
         Y(I)=YSAFE+DELY
         CALL FUN(N,T,Y,DER)
         lssdata%fjac(:,I)=(DER-DER0)/DELY
         Y(I)=YSAFE
      ENDDO
   endif
   ISTATUS(Njac) = ISTATUS(Njac) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

   CALL ros_PrepareMatrix( H,Direction,ros_Gamma(1),Singular)

   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


!~~~>   Compute the stages
Stage: DO istage = 1, ros_S

      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)

      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         CALL DCOPY(N,Fcn0,1,Fcn,1)
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL DCOPY(N,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL DAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
            K(N*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FUN(N,Tau,Ynew,Fcn)
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       END IF ! if istage == 1 elseif ros_NewF(istage)

       CALL DCOPY(N,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL DAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL DAXPY(N,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL LSS_Solve(N,lssdata,K(ioffset+1))      
!       CALL ros_Solve(Ghimj, Pivot, K(ioffset+1))

   END DO Stage

!~~~>  Compute the new solution
   CALL DCOPY(N,Y,1,Ynew,1)

   DO j=1,ros_S
         CALL DAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
   END DO

!~~~>  Compute the error estimation
   CALL DSCAL(N,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL DAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
   END DO
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac
!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      CALL DCOPY(N,Ynew,1,Y,1)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      END IF
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1)  ISTATUS(Nrej) = ISTATUS(Nrej) + 1
   END IF ! Err <= 1

   END DO UntilAccepted

   END DO TimeLoop

!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

  END SUBROUTINE ros_Integrator


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  DOUBLE PRECISION FUNCTION ros_ErrorNorm (Y, Ynew, Yerr, &
                               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

! Input arguments
   DOUBLE PRECISION, INTENT(IN) :: Y(N), Ynew(N), &
          Yerr(N), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables
   DOUBLE PRECISION :: Err, Scale, Ymax
   inTEGER  :: i
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.0

   Err = ZERO
   DO i=1,N
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/N)

   ros_ErrorNorm = MAX(Err,1.0d-10)

  END FUNCTION ros_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

!~~~> Input arguments
   DOUBLE PRECISION, INTENT(IN) :: T, Roundoff, Y(N), Fcn0(N)
!~~~> Output arguments
   DOUBLE PRECISION, INTENT(OUT) :: dFdT(N)
!~~~> Local variables
   DOUBLE PRECISION :: Delta
   DOUBLE PRECISION, PARAMETER :: ONE = 1.0, DeltaMin = 1.0E-6

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL FUN(N,T+Delta,Y,dFdT)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1
   CALL DAXPY(N,(-ONE),Fcn0,1,dFdT,1)
   CALL DSCAL(N,(ONE/Delta),dFdT,1)

  END SUBROUTINE ros_FunTimeDerivative


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_PrepareMatrix (  H, Direction, gam, Singular )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   USE LS_Solver_ROS
   IMPLICIT NONE

!~~~> Input arguments
!   DOUBLE PRECISION, INTENT(IN) ::  Jac0(N,N)
   DOUBLE PRECISION, INTENT(IN) :: gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments
!   DOUBLE PRECISION, INTENT(OUT) :: Ghimj(N,N)
   LOGICAL, INTENT(OUT) ::  Singular
!   INTEGER, INTENT(OUT) ::  Pivot(N)
!~~~> Inout arguments
   DOUBLE PRECISION, INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: ising, Nconsecutive
   DOUBLE PRECISION :: ghinv
   DOUBLE PRECISION, PARAMETER :: ONE  = 1.0, HALF = 0.5
   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)

!~~~>    Construct Ghimj = 1/(H*gam) - Jac0
!     CALL DCOPY(N*N,Jac0,1,Ghimj,1)
!     CALL DSCAL(N*N,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     call LSS_Decomp(N,lssdata,ghinv,ising)
!     DO i=1,N
!       Ghimj(i,i) = Ghimj(i,i)+ghinv
!     END DO
!~~~>    Compute LU decomposition
!     CALL ros_Decomp( Ghimj, Pivot, ising )
     IF (ising == 0) THEN
!~~~>    If successful done
        Singular = .FALSE.
     ELSE ! ising .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ising = ',ising
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF    ! ising

   END DO ! WHILE Singular

  END SUBROUTINE ros_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Decomp( A, Pivot, ising )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the LU decomposition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Inout variables CAUTION!!
   DOUBLE PRECISION, INTENT(INOUT) :: A(N,N)
!~~~> Output variables
   INTEGER, INTENT(OUT) :: Pivot(N), ising

   CALL  DGETRF( N, N, A, N, Pivot, ising )
   ISTATUS(Ndec) = ISTATUS(Ndec) + 1

  END SUBROUTINE ros_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Solve( A, Pivot, b )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution (using pre-computed LU decomposition)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Input variables CAUTION!!!
   DOUBLE PRECISION, INTENT(IN) :: A(N,N)
   INTEGER, INTENT(IN) :: Pivot(N)
   INTEGER ::ISING
!~~~> InOut variables
   DOUBLE PRECISION, INTENT(INOUT) :: b(N)

   CALL  DGETRS( 'N', N , 1, A, N, Pivot, b, N, ISING)
   ISTATUS(Nsol) = ISTATUS(Nsol) + 1

  END SUBROUTINE ros_Solve



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION g

    g = 1.0 + 1.0/SQRT(2.0)
    rosMethod = RS2
!~~~> Name of the method
    ros_Name = 'ROS-2'
!~~~> Number of stages
    ros_S = 2

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = (1.0)/g
    ros_C(1) = (-2.0)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.0)/(2.0*g)
    ros_M(2)= (1.0)/(2.0*g)
! E_i = Coefficients for error estimator
    ros_E(1) = 1.0/(2.0*g)
    ros_E(2) = 1.0/(2.0*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0
    ros_Alpha(2) = 1.0
!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = g
    ros_Gamma(2) =-g

 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   rosMethod = RS3
!~~~> Name of the method
   ros_Name = 'ROS-3'
!~~~> Number of stages
   ros_S = 3

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1)= 1.0
   ros_A(2)= 1.0
   ros_A(3)= 0.0

   ros_C(1) = -0.10156171083877702091975600115545E+01
   ros_C(2) =  0.40759956452537699824805835358067E+01
   ros_C(3) =  0.92076794298330791242156818474003E+01
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1E+01
   ros_M(2) =  0.61697947043828245592553615689730E+01
   ros_M(3) = -0.42772256543218573326238373806514
! E_i = Coefficients for error estimator
   ros_E(1) =  0.5
   ros_E(2) = -0.29079558716805469821718236208017E+01
   ros_E(3) =  0.22354069897811569627360909276199
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0
   ros_Alpha(2)= 0.43586652150845899941601945119356
   ros_Alpha(3)= 0.43586652150845899941601945119356
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1)= 0.43586652150845899941601945119356
   ros_Gamma(2)= 0.24291996454816804366592249683314
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01

  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros4
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

   rosMethod = RS4
!~~~> Name of the method
   ros_Name = 'ROS-4'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.2000000000000000E+01
   ros_A(2) = 0.1867943637803922E+01
   ros_A(3) = 0.2344449711399156
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0

   ros_C(1) =-0.7137615036412310E+01
   ros_C(2) = 0.2580708087951457E+01
   ros_C(3) = 0.6515950076447975
   ros_C(4) =-0.2137148994382534E+01
   ros_C(5) =-0.3214669691237626
   ros_C(6) =-0.6949742501781779
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735E+01
   ros_M(2) = 0.2870493262186792
   ros_M(3) = 0.4353179431840180
   ros_M(4) = 0.1093502252409163E+01
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) =-0.2815431932141155
   ros_E(2) =-0.7276199124938920E-01
   ros_E(3) =-0.1082196201495311
   ros_E(4) =-0.1093502252409163E+01
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0
   ros_Alpha(2) = 0.1145640000000000E+01
   ros_Alpha(3) = 0.6552168638155900
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5728200000000000
   ros_Gamma(2) =-0.1769193891319233E+01
   ros_Gamma(3) = 0.7592633437920482
   ros_Gamma(4) =-0.1049021087100450

  END SUBROUTINE Ros4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RD3
!~~~> Name of the method
   ros_Name = 'RODAS-3'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.0
   ros_A(2) = 2.0
   ros_A(3) = 0.0
   ros_A(4) = 2.0
   ros_A(5) = 0.0
   ros_A(6) = 1.0

   ros_C(1) = 4.0
   ros_C(2) = 1.0
   ros_C(3) =-1.0
   ros_C(4) = 1.0
   ros_C(5) =-1.0
   ros_C(6) =-(8.0/3.0)

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0
   ros_M(2) = 0.0
   ros_M(3) = 1.0
   ros_M(4) = 1.0
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) = 0.0
   ros_E(2) = 0.0
   ros_E(3) = 0.0
   ros_E(4) = 1.0
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0
   ros_Alpha(2) = 0.0
   ros_Alpha(3) = 1.0
   ros_Alpha(4) = 1.0
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5
   ros_Gamma(2) = 1.5
   ros_Gamma(3) = 0.0
   ros_Gamma(4) = 0.0

  END SUBROUTINE Rodas3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RD4
!~~~> Name of the method
    ros_Name = 'RODAS-4'
!~~~> Number of stages
    ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000
    ros_Alpha(2) = 0.386
    ros_Alpha(3) = 0.210
    ros_Alpha(4) = 0.630
    ros_Alpha(5) = 1.000
    ros_Alpha(6) = 1.000

!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = 0.2500000000000000
    ros_Gamma(2) =-0.1043000000000000
    ros_Gamma(3) = 0.1035000000000000
    ros_Gamma(4) =-0.3620000000000023E-01
    ros_Gamma(5) = 0.0
    ros_Gamma(6) = 0.0

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = 0.1544000000000000E+01
    ros_A(2) = 0.9466785280815826
    ros_A(3) = 0.2557011698983284
    ros_A(4) = 0.3314825187068521E+01
    ros_A(5) = 0.2896124015972201E+01
    ros_A(6) = 0.9986419139977817
    ros_A(7) = 0.1221224509226641E+01
    ros_A(8) = 0.6019134481288629E+01
    ros_A(9) = 0.1253708332932087E+02
    ros_A(10) =-0.6878860361058950
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0

    ros_C(1) =-0.5668800000000000E+01
    ros_C(2) =-0.2430093356833875E+01
    ros_C(3) =-0.2063599157091915
    ros_C(4) =-0.1073529058151375
    ros_C(5) =-0.9594562251023355E+01
    ros_C(6) =-0.2047028614809616E+02
    ros_C(7) = 0.7496443313967647E+01
    ros_C(8) =-0.1024680431464352E+02
    ros_C(9) =-0.3399990352819905E+02
    ros_C(10) = 0.1170890893206160E+02
    ros_C(11) = 0.8083246795921522E+01
    ros_C(12) =-0.7981132988064893E+01
    ros_C(13) =-0.3152159432874371E+02
    ros_C(14) = 0.1631930543123136E+02
    ros_C(15) =-0.6058818238834054E+01

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0
    ros_M(6) = 1.0

!~~~> E_i  = Coefficients for error estimator
    ros_E(1) = 0.0
    ros_E(2) = 0.0
    ros_E(3) = 0.0
    ros_E(4) = 0.0
    ros_E(5) = 0.0
    ros_E(6) = 1.0

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
    ros_ELO = 4.0

  END SUBROUTINE Rodas4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   End of the set of internal Rosenbrock subroutines
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
END SUBROUTINE Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END MODULE ROS_f90_Integrator
