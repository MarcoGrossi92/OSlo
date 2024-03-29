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

# if defined (__GFORTRAN__)
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
# endif

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
#if defined (_OPENMP)
  use omp_lib
#endif
  implicit none
  integer, parameter :: neq = 3, nc = 10000
  real(8) :: Y(neq,nc), t, TOUT, tlimit
  real(8) :: time1, time2, RT(neq), AT(neq)
  character(len=30) :: Format
  integer :: i, nthreads=1
  character(10) :: try
  integer :: err

#if defined (_OPENMP)
  !$omp parallel
  nthreads = OMP_GET_NUM_THREADS()
  !$omp end parallel
  write(*,*)'  Number of threads:  ',nthreads
#endif

  RT = 1.d-4
  AT = [1.D-8, 1.D-14, 1.D-6]
  tlimit = 4.d+10
  Format = '(A8,E20.8,A15)'

# if defined(__GFORTRAN__)
  call setup_odesolver(N=neq,solver='lsoda',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  !$omp parallel default(none), &
  !$omp shared(y,tlimit,run_odesolver), &
  !$omp private(i,T,TOUT,err)
  !$omp do schedule(dynamic)
  do i = 1, nc
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fodepack,jac,err)
  enddo
  !$omp enddo
  !$omp end parallel
  call cpu_time(time2)
  try = 'fail'
  if (sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))<=1d-20) try = 'success' 
  write(*,Format) 'lsoda', (time2 - time1)/nthreads, try
#endif

  call setup_odesolver(N=neq,solver='Hradau5',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  !$omp parallel default(none), &
  !$omp shared(y,tlimit,run_odesolver), &
  !$omp private(i,T,TOUT,err)
  !$omp do schedule(dynamic)
  do i = 1, nc
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,jac,err)
  enddo
  !$omp enddo
  !$omp end parallel
  call cpu_time(time2)
  try = 'fail'
  if (sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))<=1d-20) try = 'success' 
  write(*,Format) 'Hradau5', (time2 - time1)/nthreads, try

  
  call setup_odesolver(N=neq,solver='rk',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  !$omp parallel default(none), &
  !$omp shared(y,tlimit,run_odesolver), &
  !$omp private(i,T,TOUT,err)
  !$omp do schedule(dynamic)
  do i = 1, nc
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,jac,err)
  enddo
  !$omp enddo
  !$omp end parallel
  call cpu_time(time2)
  try = 'fail'
  if (sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))<=1d-20) try = 'success' 
  write(*,Format) 'rk', (time2 - time1)/nthreads, try
  if (try=='fail')write(*,*)'error = ',sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))

  call setup_odesolver(N=neq,solver='Hrodas',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  !$omp parallel default(none), &
  !$omp shared(y,tlimit,run_odesolver), &
  !$omp private(i,T,TOUT,err)
  !$omp do schedule(dynamic)
  do i = 1, nc
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,jac,err)
  enddo
  !$omp enddo
  !$omp end parallel
  call cpu_time(time2)
  try = 'fail'
  if (sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))<=1d-20) try = 'success' 
  write(*,Format) 'Hrodas', (time2 - time1)/nthreads, try

  call setup_odesolver(N=neq,solver='ros',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  !$omp parallel default(none), &
  !$omp shared(y,tlimit,run_odesolver), &
  !$omp private(i,T,TOUT,err)
  !$omp do schedule(dynamic)
  do i = 1, nc
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,jac,err)
  enddo
  !$omp enddo
  !$omp end parallel
  call cpu_time(time2)
  try = 'fail'
  if (sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))<=1d-20) try = 'success' 
  write(*,Format) 'ros', (time2 - time1)/nthreads, try

  call setup_odesolver(N=neq,solver='Hsdirk4',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  !$omp parallel default(none), &
  !$omp shared(y,tlimit,run_odesolver), &
  !$omp private(i,T,TOUT,err)
  !$omp do schedule(dynamic)
  do i = 1, nc
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,jac,err)
  enddo
  !$omp enddo
  !$omp end parallel
  call cpu_time(time2)
  try = 'fail'
  if (sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))<=1d-20) try = 'success' 
  write(*,Format) 'Hsdirk4', (time2 - time1)/nthreads, try

  call setup_odesolver(N=neq,solver='sdirk4',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  !$omp parallel default(none), &
  !$omp shared(y,tlimit,run_odesolver), &
  !$omp private(i,T,TOUT,err)
  !$omp do schedule(dynamic)
  do i = 1, nc
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,jac,err)
  enddo
  !$omp enddo
  !$omp end parallel
  call cpu_time(time2)
  try = 'fail'
  if (sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))<=1d-20) try = 'success' 
  write(*,Format) 'sdirk', (time2 - time1)/nthreads, try

# if defined(INTEL)
  call setup_odesolver(N=neq,solver='dodesol',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  !$omp parallel default(none), &
  !$omp shared(y,tlimit,run_odesolver), &
  !$omp private(i,T,TOUT,err)
  !$omp do schedule(dynamic)
  do i = 1, nc
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,jac,err)
  enddo
  !$omp enddo
  !$omp end parallel
  call cpu_time(time2)
  try = 'fail'
  if (sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))<=1d-20) try = 'success' 
  write(*,Format) 'intel', (time2 - time1)/nthreads, try
# endif

contains

  subroutine initialize()
    implicit none
    Y(1,:) = 1.0D0
    Y(2,:) = 0.0D0
    Y(3,:) = 0.0D0
  end subroutine

END PROGRAM RUNEXAMPLE1
