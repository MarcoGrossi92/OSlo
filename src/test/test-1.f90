! DEMONSTRATION PROGRAM.

! The following is a simple test problem from
! chemical kinetics. It consists of the following three rate
! equations:
!     dy1/dt = -.04d0*y1 + 1.d4*y2*y3
!     dy2/dt = .04d0*y1 - 1.d4*y2*y3 - 3.d7*y2**2
!     dy3/dt = 3.d7*y2**2
! on the interval from t = 0.0d0 to t = 4.d10, with initial
! conditions y1 = 1.0d0, y2 = y3 = 0.0d0. The problem is stiff.
! It uses ITOL = 2 and ATOL much smaller for y2 than y1 or y3 
! because y2 has much smaller values.

MODULE functions

CONTAINS

  SUBROUTINE FEX(NEQ,T,Y,YDOT)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: NEQ
    DOUBLE PRECISION, INTENT (IN) :: T
    DOUBLE PRECISION, INTENT (IN) :: Y(NEQ)
    DOUBLE PRECISION, INTENT (OUT) :: YDOT(NEQ)
    YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
    YDOT(3) = 3.E7*Y(2)*Y(2)
    YDOT(2) = -YDOT(1) - YDOT(3)
    RETURN
  END SUBROUTINE FEX

END MODULE functions


module ode
  use functions
  implicit none

contains

subroutine radau(n,t1,t2,Z)
  implicit none
  integer, intent(in) :: n
  real(8), intent(inout) :: t1, t2
  real(8), intent(inout) :: Z(n)
  real(8) :: h, ATOL(3), RTOL(1)
  integer :: ITOL, IJAC, MLJAC, MUJAC, IOUT, IDID, MUMAS, MLMAS, IMAS
  integer :: LWORK, LIWORK, LRCONT
  real(8) :: WORK(4*(n+1)*(n+1)+12*(n+1)+20)
  integer :: IWORK(3*(n+1)+20)
  real(8) :: RPAR(1)=0d0
  integer :: IPAR(1)=0
  external :: RADAU5

  LWORK=4*n*n+12*n+20
  LIWORK=3*n+20
  LRCONT=5*n+2
  IWORK = 0
  
  h = 0.D0
  RTOL = 1.d-4
  ATOL(1) = 1.D-8
  ATOL(2) = 1.D-14
  ATOL(3) = 1.D-6
  ITOL = 0
  IJAC = 0
  IMAS = 0
  MLJAC = n
  MLMAS = n
  MUJAC = 0
  MUMAS = 0
  IOUT = 0
  WORK = 0d0

  WORK(1)=1.1d-19

  !% Implicit solver with automatic numerical Jacobi matrix computations
  call RADAU5(n,fex,t1,Z,t2,H,                               &
                RTOL,ATOL,ITOL,                                 &
                dummy,IJAC,MLJAC,MUJAC,                          &
                dummy,IMAS,MLMAS,MUMAS,                          &
                dummy,IOUT,                                    &
                WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)

end subroutine radau

subroutine dvode(n,t1,t2,Z)
  use DVODE_F90_M
  implicit none
  integer, intent(in) :: n
  real(8), intent(inout) :: t1, t2
  real(8), intent(inout) :: Z(n)
  real(8) :: ATOL(n), RTOL
  integer :: ITASK, ISTATE
  TYPE (VODE_OPTS) :: OPTIONS

  RTOL = 1.D-4
  ATOL(1) = 1.D-8
  ATOL(2) = 1.D-14
  ATOL(3) = 1.D-6
  ITASK = 1
  ISTATE = 1
  OPTIONS = SET_NORMAL_OPTS(DENSE_J=.TRUE.,ABSERR_VECTOR=ATOL,      &
        RELERR=RTOL)
  CALL DVODE_F90(FEX,n,Z,t1,t2,ITASK,ISTATE,OPTIONS)

end subroutine dvode

subroutine dummy()
endsubroutine

end module ode


PROGRAM RUNEXAMPLE1
  use ode
  implicit none
  integer, parameter :: neq = 3
  real(8) :: Y(neq), t, TOUT

  Y(1) = 1.0D0
  Y(2) = 0.0D0
  Y(3) = 0.0D0
  T = 0.0D0
  TOUT = 40.D0

  call dvode(neq,T,TOUT,Y)
  write(*,'(A,4F15.8)') 'dvode', Y(:), sum(y)

  Y(1) = 1.0D0
  Y(2) = 0.0D0
  Y(3) = 0.0D0
  T = 0.0D0
  TOUT = 40.D0

  call radau(neq,T,TOUT,Y)
  write(*,'(A,4F15.8)') 'radau', Y(:), sum(y)

END PROGRAM RUNEXAMPLE1
