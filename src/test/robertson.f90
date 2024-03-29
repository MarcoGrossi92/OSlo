! DEMONSTRATION PROGRAM.
!
! A classic example of a stiff system of ODEs is the kinetic analysis of Robertson's autocatalytic chemical reaction:
! H. H. Robertson, The solution of a set of reaction rate equations, in J. Walsh (Ed.), 
! Numerical Analysis: An Introduction, pp. 178–182, Academic Press, London (1966).
! It consists of the following three rate equations:
!    xdot = -0.04 * x + 1.e4 * y * z
!    ydot = 0.04 * x - 1.e4 * y * z - 3.e7 * y**2
!    zdot = 3.e7 * y**2
! with initial conditions y1 = 1.0d0, y2 = y3 = 0.0d0. The problem is stiff.
! Absolute tolleranvce is set much smaller for y than x or z because y has much smaller values.

MODULE functions

CONTAINS

  SUBROUTINE Fgeneral(NEQ,T,Y,YDOT)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: NEQ
    real(8), INTENT (IN) :: T
    real(8), INTENT (IN) :: Y(NEQ)
    real(8), INTENT (OUT) :: YDOT(NEQ)
    YDOT(1) = -0.04d0*Y(1) + 1.d4*Y(2)*Y(3)
    YDOT(2) = 0.04d0*Y(1) - 1.d4 *Y(2)*Y(3) - 3.d7 * Y(2)**2
    YDOT(3) = 3.E7*Y(2)*Y(2)
  END SUBROUTINE Fgeneral

  SUBROUTINE Fdvode(NEQ,T,Y,YDOT,RPAR,IPAR)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: NEQ
    real(8), INTENT (IN) :: T
    real(8), INTENT (IN) :: Y(NEQ)
    real(8), INTENT (OUT) :: YDOT(NEQ)
    real(8), intent(in) :: RPAR(*)
    integer, intent(in) :: IPAR(*)
    YDOT(1) = -0.04d0*Y(1) + 1.d4*Y(2)*Y(3)
    YDOT(2) = 0.04d0*Y(1) - 1.d4 *Y(2)*Y(3) - 3.d7 * Y(2)**2
    YDOT(3) = 3.E7*Y(2)*Y(2)
  END SUBROUTINE Fdvode

  subroutine Fodepack(self, neq, t, y, ydot, ierr)
    use odepack_mod
    class(lsoda_class), intent(inout) :: self
    integer, intent(in) :: neq
    real(8), intent(in) :: t
    real(8), intent(in) :: y(neq)
    real(8), intent(out) :: ydot(neq)
    integer, intent(out) :: ierr
    YDOT(1) = -0.04d0*Y(1) + 1.d4*Y(2)*Y(3)
    YDOT(2) = 0.04d0*Y(1) - 1.d4 *Y(2)*Y(3) - 3.d7 * Y(2)**2
    YDOT(3) = 3.E7*Y(2)*Y(2)
    ierr = 0
  end subroutine Fodepack

  SUBROUTINE JAC(NVAR,T, V, JF)
    IMPLICIT NONE
    integer :: NVAR
    real(8) :: V(NVAR), T
    real(8) :: JF(NVAR,NVAR)
    real(8) :: UROUND, YSAFE, DELY
    real(8) :: DER(NVAR), DER0(NVAR)
    integer :: I,J

    JF(:,:) = 0.0d0
    UROUND = 1D-19
    
    DO I = 1, NVAR
      YSAFE=V(I)
      CALL Fgeneral(NVAR,T,V,DER0)
      DELY=DSQRT(UROUND*MAX(1.D-5,ABS(YSAFE)))
      V(I)=YSAFE+DELY
      CALL Fgeneral(NVAR,T,V,DER)
      DO J = 1, NVAR
        JF(J,I)=(DER(J)-DER0(J))/DELY
      ENDDO
      V(I)=YSAFE
    ENDDO

  END SUBROUTINE JAC

END MODULE functions


PROGRAM RUNEXAMPLE1
  use oslo
  use functions
  implicit none
  integer, parameter :: neq = 3
  real(8) :: Y(neq), t, TOUT, tlimit
  real(8) :: time1, time2, RT(neq), AT(neq)
  character(len=30) :: Format
  integer :: err

  RT = 1.d-4
  AT = [1.D-8, 1.D-14, 1.D-6]
  tlimit = 4.d+10
  Format = '(A8,4E20.8)'

  ! call setup_odesolver(N=neq,solver='dvode',RT=RT,AT=AT)
  ! call initialize
  ! call cpu_time(time1)
  ! call run_odesolver(neq,T,TOUT,Y,Fdvode)
  ! call cpu_time(time2)
  ! write(*,Format) 'dvode', time2-time1, Y(:)

  ! call setup_odesolver(N=neq,solver='dvodef90',RT=RT,AT=AT)
  ! call initialize
  ! call cpu_time(time1)
  ! call run_odesolver(neq,T,TOUT,Y,Fgeneral)
  ! call cpu_time(time2)
  ! write(*,Format) 'dvodef90', time2-time1, Y(:)

# if defined(__GFORTRAN__)
  call setup_odesolver(N=neq,solver='lsoda',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  call run_odesolver(neq,T,TOUT,Y,Fodepack,JAC,err)
  call cpu_time(time2)
  write(*,Format) 'lsoda', time2-time1, Y(:)
# endif

  call setup_odesolver(N=neq,solver='Hradau5',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  call run_odesolver(neq,T,TOUT,Y,Fgeneral,JAC,err)
  call cpu_time(time2)
  write(*,Format) 'Hradau5', time2-time1, Y(:)

  ! call setup_odesolver(N=neq,solver='radau',RT=RT,AT=AT)
  ! call initialize
  ! call cpu_time(time1)
  ! call run_odesolver(neq,T,TOUT,Y,Fgeneral)
  ! call cpu_time(time2)
  ! write(*,Format) 'radau', time2-time1, Y(:)

  call setup_odesolver(N=neq,solver='rk',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  call run_odesolver(neq,T,TOUT,Y,Fgeneral,JAC,err)
  call cpu_time(time2)
  write(*,Format) 'rk', time2-time1, Y(:)

  call setup_odesolver(N=neq,solver='Hrodas',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  call run_odesolver(neq,T,TOUT,Y,Fgeneral,JAC,err)
  call cpu_time(time2)
  write(*,Format) 'Hrodas', time2-time1, Y(:)

  call setup_odesolver(N=neq,solver='ros',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  call run_odesolver(neq,T,TOUT,Y,Fgeneral,JAC,err)
  call cpu_time(time2)
  write(*,Format) 'ros', time2-time1, Y(:)

  call setup_odesolver(N=neq,solver='Hsdirk4',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  call run_odesolver(neq,T,TOUT,Y,Fgeneral,JAC,err)
  call cpu_time(time2)
  write(*,Format) 'Hsdirk4', time2-time1, Y(:)

  call setup_odesolver(N=neq,solver='sdirk4',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  call run_odesolver(neq,T,TOUT,Y,Fgeneral,JAC,err)
  call cpu_time(time2)
  write(*,Format) 'sdirk', time2-time1, Y(:)

# if defined(INTEL)
  call setup_odesolver(N=neq,solver='dodesol',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  call run_odesolver(neq,T,TOUT,Y,Fgeneral,JAC,err)
  call cpu_time(time2)
  write(*,Format) 'intel', time2-time1, Y(:)
# endif

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
