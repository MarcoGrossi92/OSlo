! OSLO - ODE SoLver toolbOx
!
! ODE solvers are set to work with the following specifications:
!   - full and unknown Jacobian evaluated by the solver
!   - the mass-matrix is the identity.
! All the parameters in the "setup" subroutine and within each solver procedure are set accordingly.

module oslo
  implicit none

  real(8), allocatable :: RTOL(:)
  real(8), allocatable :: ATOL(:)

  real(8) :: h

  integer :: ITOL
  integer :: IJAC, MLJAC, MUJAC, IOUT, IDID, MUMAS, MLMAS, IMAS
  integer :: LWORK, LIWORK
  integer :: LRCONT     ! see sdirk4.f
  integer :: NSMAX      ! see radau.f

contains

  subroutine setup(N,solver,RT,AT)
    implicit none
    integer, intent(in)           :: N
    character(len=*), intent(in)  :: solver
    real(8), optional, intent(in) :: RT(N), AT(N)

    ! Tollerances are defined as arrays -> ITOL = 1 (Hairer); ITOL = 2 (DVODE)
    ITOL = 1

    if (allocated(RTOL)) deallocate(RTOL)
    if (allocated(ATOL)) deallocate(ATOL)
    allocate(RTOL(1:N)); allocate(ATOL(1:N))
    ! Default values
    RTOL = 1.D-4; ATOL = 1.D-10
    if (present(RT)) RTOL = RT
    if (present(AT)) ATOL = AT

    select case(solver)
    case('dvode')
      ITOL = ITOL+1
      LWORK = 22 +  9*N + 2*N**2
      LIWORK= 30 + N
    case('dvodef90')
    case('lsoda')
    case('radau5')
      h = 0.D0
      LWORK=4*(n)*(n)+12*(n)+20
      LIWORK=3*n+20
    case('radau')
      h = 0.D0
      NSMAX = 7
      LWORK = (NSMAX+1)*N*N+(3*NSMAX+3)*N+20
      LIWORK= (2+(NSMAX-1)/2)*N+20
    case('rodas')
      h = 0.D0
      LWORK = 2*N*N+14*N+20
      LIWORK= N+20
    case('sdirk4')
      h = 0.D0
      LWORK=2*N*N+12*N+7
      LIWORK=2*N+4
      LRCONT=5*n+2
    case('dodesol')
      h = 1.d-15
      LWORK = (7+2*n)*n
      LIWORK = 128
    end select

  end subroutine setup


# if defined(INTEL)
subroutine wrap_dodesol(n,t1,t2,Z,fcn)
  implicit none
  integer, intent(in) :: n
  real(8), intent(inout) :: t1, t2
  real(8), intent(inout) :: Z(n)
  external :: fcn
  ! specific
  external :: dodesol_mk52lfn
  integer :: kd(n), ierr
  real(8) :: hm, ep, tr
  real(8) :: dpar(LWORK)
  integer :: ipar(LIWORK)
  
  hm = 1.d-15
  ep = minval(RTOL)
  tr = minval(ATOL)
  ipar = 0

  call dodesol_mk52lfn(ipar,n,t1,t2,Z,fcn,h,hm,ep,tr,dpar,kd,ierr)

end subroutine wrap_dodesol
# endif


subroutine wrap_sdirk4(n,t1,t2,Z,fcn)
  implicit none
  integer, intent(in) :: n
  real(8), intent(inout) :: t1, t2
  real(8), intent(inout) :: Z(n)
  external :: fcn
  ! specific
  external :: SDIRK4
  real(8) :: WORK(LWORK)
  integer :: IWORK(LIWORK)

  IJAC = 0
  IMAS = 0
  MLJAC = n
  MLMAS = n
  MUJAC = 0
  MUMAS = 0
  IOUT = 0

  IWORK = 0
  WORK = 0d0
  WORK(1)=1.1d-19

  call SDIRK4(n,fcn,t1,Z,t2,H,                         &
                RTOL,ATOL,ITOL,                        &
                dummy,IJAC,MLJAC,MUJAC,                &
                dummy,IMAS,MLMAS,MUMAS,                &
                dummy,IOUT,                            &
                WORK,LWORK,IWORK,LIWORK,LRCONT,IDID)

end subroutine wrap_sdirk4


subroutine wrap_radau5(n,t1,t2,Z,fcn)
  implicit none
  integer, intent(in) :: n
  real(8), intent(inout) :: t1, t2
  real(8), intent(inout) :: Z(n)
  external :: fcn
  ! specific
  external :: RADAU5
  real(8) :: WORK(LWORK)
  integer :: IWORK(LIWORK)
  real(8) :: RPAR(1)=0d0
  integer :: IPAR(1)=0

  IJAC = 0
  IMAS = 0
  MLJAC = n
  MLMAS = n
  MUJAC = 0
  MUMAS = 0
  IOUT = 0

  IWORK = 0
  WORK = 0d0
  WORK(1)=1.1d-19

  call RADAU5(n,fcn,t1,Z,t2,H,                                  &
                RTOL,ATOL,ITOL,                                 &
                dummy,IJAC,MLJAC,MUJAC,                         &
                dummy,IMAS,MLMAS,MUMAS,                         &
                dummy,IOUT,                                     &
                WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)

end subroutine wrap_radau5


subroutine wrap_radau(n,t1,t2,Z,fcn)
  implicit none
  integer, intent(in) :: n
  real(8), intent(inout) :: t1, t2
  real(8), intent(inout) :: Z(n)
  external :: fcn
  ! specific
  external :: RADAU
  real(8) :: WORK(LWORK)
  integer :: IWORK(LIWORK)
  real(8) :: RPAR(1)=0d0
  integer :: IPAR(1)=0

  IJAC = 0
  IMAS = 0
  MLJAC = n
  MLMAS = n
  MUJAC = 0
  MUMAS = 0
  IOUT = 0

  IWORK = 0
  WORK = 0d0
  WORK(1)=1.1d-19

  call RADAU(n,fcn,t1,Z,t2,H,                                   &
                RTOL,ATOL,ITOL,                                 &
                dummy,IJAC,MLJAC,MUJAC,                         &
                dummy,IMAS,MLMAS,MUMAS,                         &
                dummy,IOUT,                                     &
                WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)

end subroutine wrap_radau


subroutine wrap_rodas(n,t1,t2,Z,fcn)
  implicit none
  integer, intent(in) :: n
  real(8), intent(inout) :: t1, t2
  real(8), intent(inout) :: Z(n)
  external :: fcn
  ! specific
  external :: RODAS
  real(8) :: WORK(LWORK)
  integer :: IWORK(LIWORK)
  real(8) :: RPAR(1)=0d0
  integer :: IPAR(1)=0
  integer :: IFCN, IDFX

  IFCN = 0
  IDFX = 0
  IJAC = 0
  IMAS = 0
  MLJAC = n
  MLMAS = n
  MUJAC = 0
  MUMAS = 0
  IOUT = 0

  IWORK = 0
  WORK = 0d0
  WORK(1)=1.1d-19

  call RODAS(n,fcn,IFCN,t1,Z,t2,H,                              &
                RTOL,ATOL,ITOL,                                 &
                dummy,IJAC,MLJAC,MUJAC,dummy,IDFX,              &
                dummy,IMAS,MLMAS,MUMAS,                         &
                dummy,IOUT,                                     &
                WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)

end subroutine wrap_rodas


subroutine wrap_dvodef90OMP(n,t1,t2,Z,fcn)
  use DVODE_F90_M
  implicit none
  integer, intent(in) :: n
  real(8), intent(inout) :: t1, t2
  real(8), intent(inout) :: Z(n)
  external :: fcn
  ! specific
  integer :: ITASK, ISTATE
  TYPE (VODE_OPTS) :: OPTIONS

  ITASK = 1
  ISTATE = 1
  OPTIONS = SET_OPTS(DENSE_J=.TRUE.,ABSERR_VECTOR=ATOL,RELERR_VECTOR=RTOL)
  CALL DVODE_F90(fcn,n,z,t1,t2,ITASK,ISTATE,OPTIONS)

end subroutine wrap_dvodef90OMP


subroutine wrap_dvode(n,t1,t2,Z,fcn)
  implicit none
  integer, intent(in) :: n
  real(8), intent(inout) :: t1, t2
  real(8), intent(inout) :: Z(n)
  external :: fcn
  ! specific
  external :: dvode
  integer :: ITASK, ISTATE, MF, IOPT
  real(8) :: WORK(LWORK)
  integer :: IWORK(LIWORK)
  integer :: IPAR(1)=0
  real(8) :: RPAR(1)=0d0

  ITASK = 1
  ISTATE = 1
  IOPT = 0
  MF = 22

  call DVODE (fcn, n, z, t1, t2, ITOL, RTOL, ATOL, ITASK,          &
              ISTATE, IOPT, WORK, LWORK, IWORK, LIWORK, dummy, MF, &
              RPAR, IPAR)

end subroutine wrap_dvode


subroutine wrap_odepack(n,t1,t2,Z,fcn)
  use odepack_mod
  implicit none
  integer, intent(in) :: n
  real(8), intent(inout) :: t1, t2
  real(8), intent(inout) :: Z(n)
  external :: fcn
  ! specific
  type(lsoda_class) :: eq
  integer :: itask, istate

  itask = 1
  istate = 1
  call eq%initialize(fcn, n, istate=istate)
  call eq%integrate(z, t1, t2, minval(RTOL), ATOL, itask, istate)

end subroutine wrap_odepack

subroutine dummy()
endsubroutine

end module oslo