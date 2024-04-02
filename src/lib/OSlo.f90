! OSLO - ODE SoLver toolbOx
!
! ODE solvers are set to work with the following specifications:
!   - full and unknown Jacobian evaluated by the solver
!   - the mass-matrix is the identity.
! All the parameters in the "setup" subroutine and within each solver procedure are set accordingly.
!
!***************************************************************************************************
!
! INTEGER INPUT
!
! Hairer solvers
!
! radau/radau5 -> 5
! 1) /=0 -> JACOBIAN MATRIX TO HESSENBERG FORM
! 2) max steps number
! 3) max number of newton iterations
! 4) choice of step initial value
! 8) gustafsson vs classic

! rodas -> 3
! 1) max steps number
! 2) choice of method
! 3) gustafsson vs classic

! sdirk4 -> 4
! 1) /=0 -> JACOBIAN MATRIX TO HESSENBERG FORM
! 2) max steps number
! 3) max number of newton iterations
! 4) coefficients

! FATODE solvers

! rk -> 6
! 3) method
! 4) max steps number
! 5) max number of newton iterations
! 6) choice of newton initial value
! 10) error estimation strategy
! 11) gustafsson vs classic

! ros -> 2
! 3) method
! 4) max steps number

! sdirk -> 4
! 3) method
! 4) max steps number
! 5) max number of newton iterations
! 6) choice of newton initial value

module oslo
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  implicit none
  private
  public :: setup_odesolver
  public :: run_odesolver

  !> Concrete procedure pointing to one of the subroutine realizations
  procedure(wrapper_if), pointer :: run_odesolver

  !> Abstract interface relative to the generic procedure
  abstract interface
  subroutine wrapper_if(n,t1,t2,var,fcn,jac,err)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    integer, intent(in)     :: n
    real(R8), intent(inout) :: t1, t2
    real(R8), intent(inout) :: var(n)
    integer, intent(out)    :: err
    external :: fcn, jac
  end subroutine wrapper_if
  end interface

  real(R8), allocatable :: RTOL(:)
  real(R8), allocatable :: ATOL(:)

  integer, parameter :: NIOPT=3
  integer, allocatable :: IWORK_global(:)
  integer :: ITOL
  integer :: IJAC, MLJAC, MUJAC, MUMAS, MLMAS, IMAS
  integer :: LWORK, LIWORK
  integer :: LRCONT     ! see sdirk4.f
  integer :: NSMAX      ! see radau.f
  integer, parameter :: NNZERO=19 ! see FATODE

contains

  subroutine setup_odesolver(N,solver,RT,AT,IOPT)
    implicit none
    integer, intent(in)            :: N
    character(len=*), intent(in)   :: solver
    real(R8), intent(in), optional :: RT(N), AT(N)
    integer, intent(in), optional  :: IOPT(NIOPT)
    integer :: Nsteps, Nnewt, icnt

    ! Tollerances are defined as arrays -> ITOL = 1 (Hairer); ITOL = 2 (DVODE)
    ITOL = 1

    Nsteps = 0; Nnewt = 0; icnt = 0
    if (present(IOPT)) then
      Nsteps = IOPT(1)
      Nnewt = IOPT(2)
      icnt = IOPT(3)
    endif

    if (allocated(IWORK_global)) deallocate(IWORK_global)

    if (allocated(RTOL)) deallocate(RTOL)
    if (allocated(ATOL)) deallocate(ATOL)
    nullify(run_odesolver)
    allocate(RTOL(1:N)); allocate(ATOL(1:N))
    ! Default values
    RTOL = 1.D-4; ATOL = 1.D-10
    if (present(RT)) RTOL = RT
    if (present(AT)) ATOL = AT

    select case(solver)
    ! case('dvode')
    !   ITOL = ITOL+1
    !   LWORK = 22 +  9*N + 2*N**2
    !   LIWORK= 30 + N
    !   run_odesolver => wrap_dvode
    case('dvodef90')
      run_odesolver => wrap_dvodef90OMP
#   if defined(__GFORTRAN__)
    case('lsoda')
      run_odesolver => wrap_odepack
#   endif
    case('H-radau5')
      LWORK=4*(n)*(n)+12*(n)+20
      LIWORK=3*n+20
      allocate(IWORK_global(LIWORK))
      IWORK_global = 0
      IWORK_global(2) = Nsteps
      IWORK_global(3) = Nnewt
      IWORK_global(8) = icnt
      run_odesolver => wrap_radau5
    case('H-radau')
      NSMAX = 7
      LWORK = (NSMAX+1)*N*N+(3*NSMAX+3)*N+20
      LIWORK= (2+(NSMAX-1)/2)*N+20
      allocate(IWORK_global(LIWORK))
      IWORK_global = 0
      IWORK_global(2) = Nsteps
      IWORK_global(3) = Nnewt
      IWORK_global(8) = icnt
      run_odesolver => wrap_radau
    case('radau2a','lobatto3c','gauss','radau1a')
      allocate(IWORK_global(NNZERO+1))
      IWORK_global = 0
      IWORK_global(4) = Nsteps
      IWORK_global(5) = Nnewt
      IWORK_global(11) = icnt
      IWORK_global(10) = 1
      if(solver=='radau2a')IWORK_global(3)=1
      if(solver=='lobatto3c')IWORK_global(3)=2
      if(solver=='gauss')IWORK_global(3)=3
      if(solver=='radau1a')IWORK_global(3)=4
      run_odesolver => wrap_rk_FATODE
    case('H-rodas')
      LWORK = 2*N*N+14*N+20
      LIWORK= N+20
      allocate(IWORK_global(LIWORK))
      IWORK_global = 0
      IWORK_global(1) = Nsteps
      IWORK_global(2) = Nnewt
      IWORK_global(3) = icnt
      run_odesolver => wrap_rodas
    case('ros2','ros3','ros4','rodas3','rodas4')
      allocate(IWORK_global(NNZERO+1))
      IWORK_global = 0
      IWORK_global(4) = Nsteps
      if(solver=='ros2')IWORK_global(3)=1
      if(solver=='ros3')IWORK_global(3)=2
      if(solver=='ros4')IWORK_global(3)=3
      if(solver=='rodas3')IWORK_global(3)=4
      if(solver=='rodas4')IWORK_global(3)=5
      run_odesolver => wrap_ros_FATODE
    case('H-sdirk4')
      LWORK=2*N*N+12*N+7
      LIWORK=2*N+4
      LRCONT=5*n+2
      allocate(IWORK_global(LIWORK))
      IWORK_global = 0
      IWORK_global(2) = Nsteps
      IWORK_global(3) = Nnewt
      run_odesolver => wrap_sdirk4
    case('sdirk2a','sdirk2b','sdirk3a','sdirk4a','sdirk4b')
      allocate(IWORK_global(NNZERO+1))
      IWORK_global = 0
      IWORK_global(4) = Nsteps
      IWORK_global(5) = Nnewt
      if(solver=='sdirk2a')IWORK_global(3)=1
      if(solver=='sdirk2b')IWORK_global(3)=2
      if(solver=='sdirk3a')IWORK_global(3)=3
      if(solver=='sdirk4a')IWORK_global(3)=4
      if(solver=='sdirk4b')IWORK_global(3)=5
      run_odesolver => wrap_sdirk_FATODE
#   if defined(INTEL)
    case('dodesol')
      LWORK = (7+2*n)*n
      LIWORK = 128
      run_odesolver => wrap_dodesol
#   endif
    case default
      write(*,*) trim(solver), ' is not a valid integrator'
      stop
    end select

  end subroutine setup_odesolver


# if defined(INTEL)
subroutine wrap_dodesol(n,t1,t2,var,fcn,jac,err)
  integer, intent(in)     :: n
  real(R8), intent(inout) :: t1, t2
  real(R8), intent(inout) :: var(n)
  integer, intent(out)    :: err
  external :: fcn, jac
  ! specific
  external :: dodesol_mk52lfn
  real(R8) :: h
  integer :: kd(n), ierr
  real(R8) :: hm, ep, tr
  real(R8) :: dpar(LWORK)
  integer :: ipar(LIWORK)
  
  h = 1.d-15
  hm = 1.d-15
  ep = minval(RTOL)
  tr = minval(ATOL)
  ipar = IWORK_global
  err = 0

  call dodesol_mk52lfn(ipar,n,t1,t2,var,fcn,h,hm,ep,tr,dpar,kd,ierr)

end subroutine wrap_dodesol
# endif


subroutine wrap_sdirk4(n,t1,t2,var,fcn,jac,IDID)
  integer, intent(in)     :: n
  real(R8), intent(inout) :: t1, t2
  real(R8), intent(inout) :: var(n)
  integer, intent(out)    :: IDID
  external :: fcn, jac
  ! specific
  external :: SDIRK4
  real(R8) :: h
  real(R8) :: WORK(LWORK)
  integer :: IWORK(LIWORK)
  integer :: NN,NN2,NN3,NN4,IOUT
  real(R8) :: XOLD,HSOL,RCONT(LRCONT-2)
  integer :: NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL

  h = 0.D0
  IJAC = 0
  IMAS = 0
  MLJAC = n
  MLMAS = n
  MUJAC = 0
  MUMAS = 0

  IWORK = IWORK_global
  WORK = 0d0
  WORK(1)=1.1d-19

  call SDIRK4(n,fcn,t1,var,t2,H,                         &
                RTOL,ATOL,ITOL,                        &
                dummy,IJAC,MLJAC,MUJAC,                &
                dummy,IMAS,MLMAS,MUMAS,                &
                dummy,IOUT,                            &
                WORK,LWORK,IWORK,LIWORK,LRCONT,IDID,   &
                NN,NN2,NN3,NN4,XOLD,HSOL,RCONT,        &
                NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL)

end subroutine wrap_sdirk4


subroutine wrap_radau5(n,t1,t2,var,fcn,jac,IDID)
  integer, intent(in)     :: n
  real(R8), intent(inout) :: t1, t2
  real(R8), intent(inout) :: var(n)
  integer, intent(out)    :: IDID
  external :: fcn, jac
  ! specific
  external :: RADAU5
  real(R8) :: h
  real(R8) :: WORK(LWORK)
  integer :: IWORK(LIWORK)
  real(R8) :: RPAR(1)
  integer :: IPAR(1)
  real(R8) :: RTOL_(n), ATOL_(n)
  integer :: NN,NN2,NN3,NN4,IOUT
  real(R8) :: XSOL,HSOL,C2M1,C1M1
  integer :: MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG

  h = 0.D0
  IJAC = 0
  IMAS = 0
  MLJAC = n
  MLMAS = n
  MUJAC = 0
  MUMAS = 0

  IWORK = IWORK_global
  WORK = 0d0
  WORK(1)=1.1d-19

  RPAR=0d0
  IPAR=0
  
  ! Reassignment is necessary since in radau5 tollerances are modified during the calculation
  RTOL_ = RTOL
  ATOL_ = ATOL

  call RADAU5(n,fcn,t1,var,t2,H,                          &
                RTOL_,ATOL_,ITOL,                       &
                dummy,IJAC,MLJAC,MUJAC,                 &
                dummy,IMAS,MLMAS,MUMAS,                 &
                dummy,IOUT,                             &
                WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID, &
                NN,NN2,NN3,NN4,XSOL,HSOL,C2M1,C1M1,     &
                MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG)

end subroutine wrap_radau5


subroutine wrap_radau(n,t1,t2,var,fcn,jac,IDID)
  integer, intent(in)     :: n
  real(R8), intent(inout) :: t1, t2
  real(R8), intent(inout) :: var(n)
  integer, intent(out)    :: IDID
  external :: fcn, jac
  ! specific
  external :: RADAU
  real(R8) :: h
  real(R8) :: WORK(LWORK)
  integer :: IWORK(LIWORK)
  real(R8) :: RPAR(1)=0d0
  integer :: IPAR(1)=0,IOUT

  h = 0.D0
  IJAC = 0
  IMAS = 0
  MLJAC = n
  MLMAS = n
  MUJAC = 0
  MUMAS = 0

  IWORK = 0
  WORK = 0d0
  WORK(1)=1.1d-19

  call RADAU(n,fcn,t1,var,t2,H,                                 &
                RTOL,ATOL,ITOL,                                 &
                dummy,IJAC,MLJAC,MUJAC,                         &
                dummy,IMAS,MLMAS,MUMAS,                         &
                dummy,IOUT,                                     &
                WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)

end subroutine wrap_radau


subroutine wrap_rodas(n,t1,t2,var,fcn,jac,IDID)
  integer, intent(in)     :: n
  real(R8), intent(inout) :: t1, t2
  real(R8), intent(inout) :: var(n)
  integer, intent(out)    :: IDID
  external :: fcn, jac
  ! specific
  external :: RODAS
  real(R8) :: h
  real(R8) :: WORK(LWORK)
  integer :: IWORK(LIWORK)
  real(R8) :: RPAR(1)=0d0
  integer :: IPAR(1)=0
  integer :: IFCN, IDFX, IOUT

  h = 0.D0
  IFCN = 0
  IDFX = 0
  IJAC = 0
  IMAS = 0
  MLJAC = n
  MLMAS = n
  MUJAC = 0
  MUMAS = 0
  
  IWORK = IWORK_global
  WORK = 0d0
  WORK(1)=1.1d-19

  call RODAS(n,fcn,IFCN,t1,var,t2,H,                              &
                RTOL,ATOL,ITOL,                                 &
                dummy,IJAC,MLJAC,MUJAC,dummy,IDFX,              &
                dummy,IMAS,MLMAS,MUMAS,                         &
                dummy,IOUT,                                     &
                WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)

end subroutine wrap_rodas


subroutine wrap_dvodef90OMP(n,t1,t2,var,fcn,jac,err)
  use DVODE_F90_M
  implicit none
  integer, intent(in) :: n
  real(R8), intent(inout) :: t1, t2
  real(R8), intent(inout) :: var(n)
  integer, intent(out) :: err
  external :: fcn, jac
  ! specific
  integer :: ITASK, ISTATE
  TYPE (VODE_OPTS) :: OPTIONS

  ITASK = 1
  ISTATE = 1
  OPTIONS = SET_OPTS(DENSE_J=.TRUE.,ABSERR_VECTOR=ATOL,RELERR_VECTOR=RTOL)
  CALL DVODE_F90(fcn,n,var,t1,t2,ITASK,ISTATE,OPTIONS)
  err = 0

end subroutine wrap_dvodef90OMP


! subroutine wrap_dvode(n,t1,t2,Z,fcn)
!   implicit none
!   integer, intent(in) :: n
!   real(R8), intent(inout) :: t1, t2
!   real(R8), intent(inout) :: Z(n)
!   external :: fcn
!   ! specific
!   external :: dvode
!   integer :: ITASK, ISTATE, MF, IOPT
!   real(R8) :: WORK(LWORK)
!   integer :: IWORK(LIWORK)
!   integer :: IPAR(1)=0
!   real(R8) :: RPAR(1)=0d0

!   ITASK = 1
!   ISTATE = 1
!   IOPT = 0
!   MF = 22

!   call DVODE (fcn, n, z, t1, t2, ITOL, RTOL, ATOL, ITASK,          &
!               ISTATE, IOPT, WORK, LWORK, IWORK, LIWORK, dummy, MF, &
!               RPAR, IPAR)

! end subroutine wrap_dvode

#if defined(__GFORTRAN__)
subroutine wrap_odepack(n,t1,t2,var,fcn,jac,err)
  use odepack_mod
  implicit none
  integer, intent(in) :: n
  real(R8), intent(inout) :: t1, t2
  real(R8), intent(inout) :: var(n)
  integer, intent(out)    :: err
  external :: fcn, jac
  ! specific
  type(lsoda_class) :: eq
  integer :: itask, istate

  itask = 1
  istate = 1
  call eq%initialize(fcn, n, istate=istate)
  call eq%integrate(var, t1, t2, minval(RTOL), ATOL, itask, istate)
  err = 0

end subroutine wrap_odepack
#endif


subroutine wrap_sdirk_FATODE(n,t1,t2,var,fcn,jac,err)
  use SDIRK_f90_Integrator
  implicit none
  integer, intent(in)     :: n
  real(R8), intent(inout) :: t1, t2
  real(R8), intent(inout) :: var(n)
  integer, intent(out)    :: err
  external :: fcn, jac

  real(R8) :: RCNTRL(NNZERO+1), RSTATUS(NNZERO+1)
  integer  :: ICNTRL(NNZERO+1), ISTATUS(NNZERO+1)
   
  ICNTRL  = IWORK_global
  RCNTRL  = 0.0
  
  call SDIRK(N,NNZERO,T1,T2,VAR,RTOL,ATOL,fcn,JAC,  &
             RCNTRL,ICNTRL,RSTATUS,ISTATUS,err)

end subroutine wrap_sdirk_FATODE


subroutine wrap_ros_FATODE(n,t1,t2,var,fcn,jac,err)
  use ROS_f90_Integrator
  implicit none
  integer, intent(in)     :: n
  real(R8), intent(inout) :: t1, t2
  real(R8), intent(inout) :: var(n)
  integer, intent(out)    :: err
  external :: fcn, jac

  real(R8) :: RCNTRL(NNZERO+1), RSTATUS(NNZERO+1)
  integer  :: ICNTRL(NNZERO+1), ISTATUS(NNZERO+1)
   
  ICNTRL  = IWORK_global
  RCNTRL  = 0.0

  call Rosenbrock(N,NNZERO,VAR,T1,T2,   &
          ATOL,RTOL, fcn, JAC,          &
          RCNTRL,ICNTRL,RSTATUS,ISTATUS,err)

end subroutine wrap_ros_FATODE


subroutine wrap_rk_FATODE(n,t1,t2,var,fcn,jac,err)
  use RK_f90_Integrator
  implicit none
  integer, intent(in)     :: n
  real(R8), intent(inout) :: t1, t2
  real(R8), intent(inout) :: var(n)
  integer, intent(out)    :: err
  external :: fcn, jac

  real(R8) :: RCNTRL(NNZERO+1), RSTATUS(NNZERO+1)
  integer  :: ICNTRL(NNZERO+1), ISTATUS(NNZERO+1)
   
  ICNTRL  = IWORK_global
  RCNTRL  = 0.0

  call RungeKutta(  N, NNZERO, T1, T2, VAR, RTOL, ATOL, fcn, JAC, &
                    RCNTRL,ICNTRL,RSTATUS,ISTATUS,err )

  end subroutine wrap_rk_FATODE

subroutine dummy()
endsubroutine

end module oslo