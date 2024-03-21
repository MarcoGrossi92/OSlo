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

  SUBROUTINE Fhairer(NEQ,T,Y,YDOT)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: NEQ
    real(8), INTENT (IN) :: T
    real(8), INTENT (IN) :: Y(NEQ)
    real(8), INTENT (OUT) :: YDOT(NEQ)
    YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
    YDOT(3) = 3.E7*Y(2)*Y(2)
    YDOT(2) = -YDOT(1) - YDOT(3)
  END SUBROUTINE Fhairer

  SUBROUTINE Fdvode(NEQ,T,Y,YDOT,RPAR,IPAR)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: NEQ
    real(8), INTENT (IN) :: T
    real(8), INTENT (IN) :: Y(NEQ)
    real(8), INTENT (OUT) :: YDOT(NEQ)
    real(8), intent(in) :: RPAR(*)
    integer, intent(in) :: IPAR(*)
    YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
    YDOT(3) = 3.E7*Y(2)*Y(2)
    YDOT(2) = -YDOT(1) - YDOT(3)
  END SUBROUTINE Fdvode

  subroutine Fodepack(self, neq, t, y, ydot, ierr)
    use odepack_mod
    class(lsoda_class), intent(inout) :: self
    integer, intent(in) :: neq
    real(8), intent(in) :: t
    real(8), intent(in) :: y(neq)
    real(8), intent(out) :: ydot(neq)
    integer, intent(out) :: ierr
    YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
    YDOT(3) = 3.E7*Y(2)*Y(2)
    YDOT(2) = -YDOT(1) - YDOT(3)
    ierr = 0
  end subroutine Fodepack

END MODULE functions


module ode
  use functions
  implicit none

  real(8) :: RTOL = 1.d-4
  real(8) :: ATOL(3) = [1.D-8, 1.D-14, 1.D-6]

contains

subroutine wrap_sdirk(n,t1,t2,Z)
  implicit none
  integer, intent(in) :: n
  real(8), intent(inout) :: t1, t2
  real(8), intent(inout) :: Z(n)
  real(8) :: h
  integer :: ITOL, IJAC, MLJAC, MUJAC, IOUT, IDID, MUMAS, MLMAS, IMAS
  integer :: LWORK, LIWORK, LRCONT
  real(8) :: WORK(2*N*N+12*N+7)
  integer :: IWORK(2*N+4)
  external :: SDIRK4

  LWORK=2*N*N+12*N+7
  LIWORK=2*N+4
  LRCONT=5*n+2
  IWORK = 0
  
  h = 0.D0
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
  call SDIRK4(n,Fhairer,t1,Z,t2,H,                     &
                RTOL,ATOL,ITOL,                        &
                dummy,IJAC,MLJAC,MUJAC,                &
                dummy,IMAS,MLMAS,MUMAS,                &
                dummy,IOUT,                            &
                WORK,LWORK,IWORK,LIWORK,LRCONT,IDID)

end subroutine wrap_sdirk


subroutine wrap_radau5(n,t1,t2,Z)
  implicit none
  integer, intent(in) :: n
  real(8), intent(inout) :: t1, t2
  real(8), intent(inout) :: Z(n)
  real(8) :: h
  integer :: ITOL, IJAC, MLJAC, MUJAC, IOUT, IDID, MUMAS, MLMAS, IMAS
  integer :: LWORK, LIWORK
  real(8) :: WORK(4*(n)*(n)+12*(n)+20)
  integer :: IWORK(3*(n)+20)
  real(8) :: RPAR(1)=0d0
  integer :: IPAR(1)=0
  external :: RADAU5

  LWORK=4*n*n+12*n+20
  LIWORK=3*n+20
  IWORK = 0
  
  h = 0.D0
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
  call RADAU5(n,Fhairer,t1,Z,t2,H,                              &
                RTOL,ATOL,ITOL,                                 &
                dummy,IJAC,MLJAC,MUJAC,                         &
                dummy,IMAS,MLMAS,MUMAS,                         &
                dummy,IOUT,                                     &
                WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)

end subroutine wrap_radau5


subroutine wrap_radau(n,t1,t2,Z)
  implicit none
  integer, intent(in) :: n
  real(8), intent(inout) :: t1, t2
  real(8), intent(inout) :: Z(n)
  real(8) :: h
  integer :: ITOL, IJAC, MLJAC, MUJAC, IOUT, IDID, MUMAS, MLMAS, IMAS
  integer :: LWORK, LIWORK
  integer, parameter :: NSMAX=7
  real(8) :: WORK((NSMAX+1)*N*N+(3*NSMAX+3)*N+20)
  integer :: IWORK((2+(NSMAX-1)/2)*N+20)
  real(8) :: RPAR(1)=0d0
  integer :: IPAR(1)=0
  external :: RADAU

  LWORK = (NSMAX+1)*N*N+(3*NSMAX+3)*N+20
  LIWORK= (2+(NSMAX-1)/2)*N+20
  IWORK = 0
  
  h = 0.D0
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

  call RADAU(n,Fhairer,t1,Z,t2,H,                              &
                RTOL,ATOL,ITOL,                                 &
                dummy,IJAC,MLJAC,MUJAC,                         &
                dummy,IMAS,MLMAS,MUMAS,                         &
                dummy,IOUT,                                     &
                WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)

end subroutine wrap_radau


subroutine wrap_rodas(n,t1,t2,Z)
  implicit none
  integer, intent(in) :: n
  real(8), intent(inout) :: t1, t2
  real(8), intent(inout) :: Z(n)
  real(8) :: h
  integer :: ITOL, IJAC, MLJAC, MUJAC, IOUT, IDID, MUMAS, MLMAS, IMAS
  integer :: LWORK, LIWORK
  real(8) :: WORK(2*N*N+14*N+20)
  integer :: IWORK(N+20)
  real(8) :: RPAR(1)=0d0
  integer :: IPAR(1)=0, IFCN, IDFX
  external :: RODAS

  LWORK = 2*N*N+14*N+20
  LIWORK= N+20
  IWORK = 0
  
  h = 0.D0
  IFCN = 0
  IDFX = 0
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

  call RODAS(n,Fhairer,IFCN,t1,Z,t2,H,                          &
                RTOL,ATOL,ITOL,                                 &
                dummy,IJAC,MLJAC,MUJAC,dummy,IDFX,              &
                dummy,IMAS,MLMAS,MUMAS,                         &
                dummy,IOUT,                                     &
                WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)

end subroutine wrap_rodas


subroutine wrap_dvodef90OMP(n,t1,t2,Z)
  use DVODE_F90_M
  implicit none
  integer, intent(in) :: n
  real(8), intent(inout) :: t1, t2
  real(8), intent(inout) :: Z(n)
  integer :: ITASK, ISTATE
  TYPE (VODE_OPTS) :: OPTIONS

  ITASK = 1
  ISTATE = 1
  OPTIONS = SET_OPTS(DENSE_J=.TRUE.,ABSERR_VECTOR=ATOL,RELERR=RTOL)
  CALL DVODE_F90(Fhairer,n,z,t1,t2,ITASK,ISTATE,OPTIONS)

end subroutine wrap_dvodef90OMP


subroutine wrap_dvode(n,t1,t2,Z)
  implicit none
  integer, intent(in) :: n
  real(8), intent(inout) :: t1, t2
  real(8), intent(inout) :: Z(n)
  integer :: ITASK, ISTATE, ITOL, MF, LRW, IOPT
  integer :: IWORK(30 + N), LIW, IPAR(1)=0
  real(8) :: RWORK(22 +  9*N + 2*N**2), RPAR(1)=0d0
  external :: dvode

  LRW = 22 +  9*N + 2*N**2
  LIW = 30 + N

  ITASK = 1
  ISTATE = 1
  IOPT = 0
  ITOL = 2
  MF = 22
  call DVODE (Fdvode, n, z, t1, t2, ITOL, RTOL, ATOL, ITASK,     &
                ISTATE, IOPT, RWORK, LRW, IWORK, LIW, dummy, MF, &
                RPAR, IPAR)

end subroutine wrap_dvode


subroutine wrap_odepack(n,t1,t2,Z)
  use odepack_mod
  implicit none
  integer, intent(in) :: n
  real(8), intent(inout) :: t1, t2
  real(8), intent(inout) :: Z(n)
  type(lsoda_class) :: eq
  integer :: itask, istate

  itask = 1
  istate = 1
  call eq%initialize(Fodepack, n, istate=istate)
  call eq%integrate(z, t1, t2, rtol, atol, itask, istate)

end subroutine wrap_odepack

subroutine dummy()
endsubroutine

end module ode


PROGRAM RUNEXAMPLE1
  use ode
  implicit none
  integer, parameter :: neq = 3
  real(8) :: Y(neq), t, TOUT, tlimit
  real(8) :: time1, time2
  character(len=30) :: Format

  tlimit = 4.d+10
  Format = '(A8,4E20.8)'

  call initialize
  call cpu_time(time1)
  call wrap_dvode(neq,T,TOUT,Y)
  call cpu_time(time2)
  write(*,Format) 'dvode', time2-time1, Y(:)

  call initialize
  call cpu_time(time1)
  call wrap_dvodef90OMP(neq,T,TOUT,Y)
  call cpu_time(time2)
  write(*,Format) 'dvodef90', time2-time1, Y(:)

  call initialize
  call cpu_time(time1)
  call wrap_odepack(neq,T,TOUT,Y)
  call cpu_time(time2)
  write(*,Format) 'odepack', time2-time1, Y(:)

  call initialize
  call cpu_time(time1)
  call wrap_radau(neq,T,TOUT,Y)
  call cpu_time(time2)
  write(*,Format) 'radau5', time2-time1, Y(:)

  call initialize
  call cpu_time(time1)
  call wrap_radau(neq,T,TOUT,Y)
  call cpu_time(time2)
  write(*,Format) 'radau', time2-time1, Y(:)

  call initialize
  call cpu_time(time1)
  call wrap_rodas(neq,T,TOUT,Y)
  call cpu_time(time2)
  write(*,Format) 'rodas', time2-time1, Y(:)

  call initialize
  call cpu_time(time1)
  call wrap_sdirk(neq,T,TOUT,Y)
  call cpu_time(time2)
  write(*,Format) 'sdirk4', time2-time1, Y(:)

contains

  subroutine initialize()
    implicit none
    Y(1) = 1.0D0
    Y(2) = 0.0D0
    Y(3) = 0.0D0
    T = 0.0D0
    TOUT = tlimit
  end subroutine

END PROGRAM RUNEXAMPLE1
