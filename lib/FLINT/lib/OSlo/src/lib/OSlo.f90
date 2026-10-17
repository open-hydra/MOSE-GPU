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

module interface_definitions
  implicit none

  interface
    subroutine solout_if(NR, XOLD, X, Y, N, IRTRN)
      use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
      implicit none
      integer  :: NR, N
      real(R8) :: X, XOLD, Y(N)
      integer  :: IRTRN
    end subroutine solout_if
  end interface

end module interface_definitions

module oslo
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  implicit none
  private
  public :: setup_odesolver
  public :: run_odesolver
  public :: wrap_radau5      !> exposed for a direct (devirtualized) RADAU5 call in radau-only builds
  public :: wrap_radau5_jac  !> exposed so callers with an analytical Jacobian can bypass run_odesolver

  !> Concrete procedure pointing to one of the subroutine realizations
  procedure(solver_if), pointer :: run_odesolver

  abstract interface
    subroutine solver_if(n, t1, t2, var, fcn, err, hmax, solout)
      use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
      use interface_definitions, only : solout_if
      implicit none
      integer, intent(in)            :: n
      real(R8), intent(inout)        :: t1, t2
      real(R8), intent(inout)        :: var(n)
      integer, intent(out)           :: err
      real(R8), intent(in), optional :: hmax
      procedure(solout_if), optional :: solout
      external :: fcn
    end subroutine solver_if
  end interface

  !> Module data

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

  !> DIAGNOSTIC: global Radau5 step statistics, accumulated across every
  !> per-cell integration call (thread-safe via atomic). Reset/print via
  !> radau_stats_reset / radau_stats_print.
  integer(8) :: rad_nfcn=0, rad_njac=0, rad_nstep=0, rad_naccpt=0, rad_nrejct=0
  public :: radau_stats_reset, radau_stats_print

contains

  subroutine radau_stats_reset()
    implicit none
    rad_nfcn=0; rad_njac=0; rad_nstep=0; rad_naccpt=0; rad_nrejct=0
  end subroutine radau_stats_reset

  subroutine radau_stats_print(tag)
    implicit none
    character(len=*), intent(in) :: tag
    write(*,'(A,A)') ' [RADAU-STATS] ', trim(tag)
    write(*,'(A,I12)') '   function evals (NFCN)  = ', rad_nfcn
    write(*,'(A,I12)') '   jacobian evals (NJAC)  = ', rad_njac
    write(*,'(A,I12)') '   computed steps (NSTEP) = ', rad_nstep
    write(*,'(A,I12)') '   accepted steps(NACCPT) = ', rad_naccpt
    write(*,'(A,I12)') '   rejected steps(NREJCT) = ', rad_nrejct
  end subroutine radau_stats_print

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

    case('dvodef90')
      run_odesolver => wrap_dvodef90OMP

    case('H-radau5')
      LWORK=4*(n)*(n)+12*(n)+20
      LIWORK=3*n+20
      allocate(IWORK_global(LIWORK))
      IWORK_global = 0
      IWORK_global(2) = Nsteps
      IWORK_global(3) = Nnewt
      IWORK_global(8) = icnt
      run_odesolver => wrap_radau5

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
      ! RODAS needs LWORK>=4N^2+12N+20, LIWORK>=3N+20 (RADAU5-tier); the old
      ! 2N^2+14N+20 / N+20 under-sized it for n>=8 (ONERA-7) -> INSUFF. STORAGE.
      LWORK = 4*N*N+12*N+40
      LIWORK= 3*N+40
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
      allocate(IWORK_global(LIWORK))
      IWORK_global = 0
      run_odesolver => wrap_dodesol
#   endif

#   if defined(SUNDIALS)
    case('cvode')
      allocate(IWORK_global(1))
      IWORK_global = 10000
      run_odesolver => wrap_cvode
#   endif

    case default

      write(*,*) trim(solver), ' is not a valid integrator'
      stop

    end select

  end subroutine setup_odesolver


# if defined(INTEL)
subroutine wrap_dodesol(n,t1,t2,var,fcn,ierr,hmax,solout)
  use interface_definitions, only : solout_if
  implicit none
  integer, intent(in)            :: n
  real(R8), intent(inout)        :: t1, t2
  real(R8), intent(inout)        :: var(n)
  integer, intent(out)           :: ierr
  real(R8), intent(in), optional :: hmax
  external :: fcn
  procedure(solout_if), optional :: solout
  ! specific
  external :: dodesol_mk52lfn
  real(R8) :: h
  integer :: kd(n)
  real(R8) :: hm, ep, tr
  real(R8) :: dpar(LWORK)
  integer :: ipar(LIWORK)
  
  h = 1.d-15
  hm = 1.d-15
  ep = minval(RTOL)
  tr = minval(ATOL)
  ipar = IWORK_global

  call dodesol_mk52lfn(ipar,n,t1,t2,var,fcn,h,hm,ep,tr,dpar,kd,ierr)

end subroutine wrap_dodesol
# endif


subroutine wrap_sdirk4(n,t1,t2,var,fcn,IDID,hmax,solout)
  use interface_definitions, only : solout_if
  implicit none
  integer, intent(in)            :: n
  real(R8), intent(inout)        :: t1, t2
  real(R8), intent(inout)        :: var(n)
  integer, intent(out)           :: IDID
  real(R8), intent(in), optional :: hmax
  procedure(solout_if), optional :: solout
  external :: fcn
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
  ! IOUT=0: solout is never called
  ! IOUT=1: solout is available for output
  IOUT = 1

  IWORK = IWORK_global
  WORK = 0d0
  WORK(1)=1.1d-19
  if (present(hmax)) WORK(7) = hmax

  if (present(solout)) then
  call SDIRK4(n,fcn,t1,var,t2,H,                       &
                RTOL,ATOL,ITOL,                        &
                dummy,IJAC,MLJAC,MUJAC,                &
                dummy,IMAS,MLMAS,MUMAS,                &
                solout,IOUT,                           &
                WORK,LWORK,IWORK,LIWORK,LRCONT,IDID,   &
                NN,NN2,NN3,NN4,XOLD,HSOL,RCONT,        &
                NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL)
  else
  call SDIRK4(n,fcn,t1,var,t2,H,                       &
                RTOL,ATOL,ITOL,                        &
                dummy,IJAC,MLJAC,MUJAC,                &
                dummy,IMAS,MLMAS,MUMAS,                &
                dummy,IOUT,                           &
                WORK,LWORK,IWORK,LIWORK,LRCONT,IDID,   &
                NN,NN2,NN3,NN4,XOLD,HSOL,RCONT,        &
                NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL)
  endif

end subroutine wrap_sdirk4


subroutine wrap_radau5(n,t1,t2,var,fcn,IDID,hmax,solout)
  use interface_definitions, only : solout_if
  implicit none
  integer, intent(in)            :: n
  real(R8), intent(inout)        :: t1, t2
  real(R8), intent(inout)        :: var(n)
  integer, intent(out)           :: IDID
  real(R8), intent(in), optional :: hmax
  procedure(solout_if), optional :: solout
  external :: fcn
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
  IOUT = 0

  IWORK = IWORK_global
  WORK = 0d0
  WORK(1)=1.1d-19
  if (present(hmax)) WORK(7) = hmax

  RPAR=0d0
  IPAR=0

  ! Reassignment is necessary since in radau5 tollerances are modified during the calculation
  RTOL_ = RTOL
  ATOL_ = ATOL

  call RADAU5(n,fcn,t1,var,t2,H,                        &
                RTOL_,ATOL_,ITOL,                       &
                dummy,IJAC,MLJAC,MUJAC,                 &
                dummy,IMAS,MLMAS,MUMAS,                 &
                dummy,IOUT,                            &
                WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID, &
                NN,NN2,NN3,NN4,XSOL,HSOL,C2M1,C1M1,     &
                MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG)

  !$omp atomic
  rad_nfcn   = rad_nfcn   + IWORK(14)
  !$omp atomic
  rad_njac   = rad_njac   + IWORK(15)
  !$omp atomic
  rad_nstep  = rad_nstep  + IWORK(16)
  !$omp atomic
  rad_naccpt = rad_naccpt + IWORK(17)
  !$omp atomic
  rad_nrejct = rad_nrejct + IWORK(18)

end subroutine wrap_radau5


!------------------------------------------------------------------------------
! Variant of wrap_radau5 that calls RADAU5 with IJAC=1 and a user-supplied
! analytical Jacobian. Same setup as wrap_radau5; only difference is that
! `jac` is required and IJAC=1 is forced. Kept as a separate subroutine so
! the procedure-pointer signature of run_odesolver (solver_if) is untouched.
!------------------------------------------------------------------------------
subroutine wrap_radau5_jac(n,t1,t2,var,fcn,jac,IDID,hmax)
  implicit none
  integer, intent(in)            :: n
  real(R8), intent(inout)        :: t1, t2
  real(R8), intent(inout)        :: var(n)
  integer, intent(out)           :: IDID
  real(R8), intent(in), optional :: hmax
  external :: fcn
  external :: jac
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
  IJAC = 1
  IMAS = 0
  MLJAC = n
  MLMAS = n
  MUJAC = 0
  MUMAS = 0
  IOUT = 0

  IWORK = IWORK_global
  WORK = 0d0
  WORK(1)=1.1d-19
  if (present(hmax)) WORK(7) = hmax

  RPAR = 0d0
  IPAR = 0

  RTOL_ = RTOL
  ATOL_ = ATOL

  call RADAU5(n,fcn,t1,var,t2,H,                        &
                RTOL_,ATOL_,ITOL,                       &
                jac  ,IJAC,MLJAC,MUJAC,                 &
                dummy,IMAS,MLMAS,MUMAS,                 &
                dummy,IOUT,                            &
                WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID, &
                NN,NN2,NN3,NN4,XSOL,HSOL,C2M1,C1M1,     &
                MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG)

  !$omp atomic
  rad_nfcn   = rad_nfcn   + IWORK(14)
  !$omp atomic
  rad_njac   = rad_njac   + IWORK(15)
  !$omp atomic
  rad_nstep  = rad_nstep  + IWORK(16)
  !$omp atomic
  rad_naccpt = rad_naccpt + IWORK(17)
  !$omp atomic
  rad_nrejct = rad_nrejct + IWORK(18)

end subroutine wrap_radau5_jac


subroutine wrap_rodas(n,t1,t2,var,fcn,IDID,hmax,solout)
  use interface_definitions, only : solout_if
  implicit none
  integer, intent(in)            :: n
  real(R8), intent(inout)        :: t1, t2
  real(R8), intent(inout)        :: var(n)
  integer, intent(out)           :: IDID
  real(R8), intent(in), optional :: hmax
  procedure(solout_if), optional :: solout
  external :: fcn
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
  IOUT = 0
  
  IWORK = IWORK_global
  WORK = 0d0
  WORK(1)=1.1d-19
  if (present(hmax)) WORK(2) = hmax

  call RODAS(n,fcn,IFCN,t1,var,t2,H,                   &
                RTOL,ATOL,ITOL,                        &
                dummy,IJAC,MLJAC,MUJAC,dummy,IDFX,     &
                dummy,IMAS,MLMAS,MUMAS,                &
                dummy,IOUT,                           &
                WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)

end subroutine wrap_rodas


subroutine wrap_dvodef90OMP(n,t1,t2,var,fcn,err,hmax,solout)
  use DVODE_F90_M
  use interface_definitions, only : solout_if
  implicit none
  integer, intent(in)            :: n
  real(R8), intent(inout)        :: t1, t2
  real(R8), intent(inout)        :: var(n)
  integer, intent(out)           :: err
  real(R8), intent(in), optional :: hmax
  procedure(solout_if), optional :: solout
  external :: fcn
  ! specific
  integer :: ITASK, ISTATE
  TYPE (VODE_OPTS) :: OPTIONS

  ITASK = 1
  ISTATE = 1
  OPTIONS = SET_OPTS(DENSE_J=.TRUE.,ABSERR_VECTOR=ATOL,RELERR_VECTOR=RTOL)
  CALL DVODE_F90(fcn,n,var,t1,t2,ITASK,ISTATE,OPTIONS)
  err = 0

end subroutine wrap_dvodef90OMP


subroutine wrap_sdirk_FATODE(n,t1,t2,var,fcn,err,hmax,solout)
  use SDIRK_f90_Integrator
  use interface_definitions, only : solout_if
  implicit none
  integer, intent(in)            :: n
  real(R8), intent(inout)        :: t1, t2
  real(R8), intent(inout)        :: var(n)
  integer, intent(out)           :: err
  real(R8), intent(in), optional :: hmax
  procedure(solout_if), optional :: solout
  external :: fcn
  ! specific
  real(R8) :: RCNTRL(NNZERO+1), RSTATUS(NNZERO+1)
  integer  :: ICNTRL(NNZERO+1), ISTATUS(NNZERO+1)
   
  ICNTRL  = IWORK_global
  RCNTRL  = 0.0

  if (present(hmax)) RCNTRL(2) = hmax
  
  call SDIRK(N,NNZERO,T1,T2,VAR,RTOL,ATOL,fcn,dummy,  &
             RCNTRL,ICNTRL,RSTATUS,ISTATUS,err)

end subroutine wrap_sdirk_FATODE


subroutine wrap_ros_FATODE(n,t1,t2,var,fcn,err,hmax,solout)
  use ROS_f90_Integrator, only: Rosenbrock
  use interface_definitions, only : solout_if
  implicit none
  integer, intent(in)            :: n
  real(R8), intent(inout)        :: t1, t2
  real(R8), intent(inout)        :: var(n)
  integer, intent(out)           :: err
  real(R8), intent(in), optional :: hmax
  procedure(solout_if), optional :: solout
  external :: fcn
  ! specific
  real(R8) :: RCNTRL(NNZERO+1), RSTATUS(NNZERO+1)
  integer  :: ICNTRL(NNZERO+1), ISTATUS(NNZERO+1)
   
  ICNTRL  = IWORK_global
  RCNTRL  = 0.0
  RSTATUS = 0.0
  ISTATUS = 0

  if (present(hmax)) RCNTRL(2) = hmax

  call Rosenbrock(N,NNZERO,VAR,T1,T2,   &
          ATOL,RTOL, fcn, dummy,        &
          RCNTRL,ICNTRL,RSTATUS,ISTATUS,err)

end subroutine wrap_ros_FATODE


subroutine wrap_rk_FATODE(n,t1,t2,var,fcn,err,hmax,solout)
  use RK_f90_Integrator
  use interface_definitions, only : solout_if
  implicit none
  integer, intent(in)            :: n
  real(R8), intent(inout)        :: t1, t2
  real(R8), intent(inout)        :: var(n)
  integer, intent(out)           :: err
  real(R8), intent(in), optional :: hmax
  procedure(solout_if), optional :: solout
  external :: fcn
  ! specific
  real(R8) :: RCNTRL(NNZERO+1), RSTATUS(NNZERO+1)
  integer  :: ICNTRL(NNZERO+1), ISTATUS(NNZERO+1)
   
  ICNTRL  = IWORK_global
  RCNTRL  = 0.0
  if (present(hmax)) RCNTRL(2) = hmax

  call RungeKutta(  N, NNZERO, T1, T2, VAR, RTOL, ATOL, fcn, dummy, &
                    RCNTRL,ICNTRL,RSTATUS,ISTATUS,err )

end subroutine wrap_rk_FATODE

# if defined(SUNDIALS)
subroutine wrap_cvode(n,t1,t2,var,fcn,err,hmax,solout)
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod
  use fcvodes_mod                   ! Fortran interface to CVODES
  use fnvector_serial_mod           ! Fortran interface to serial N_Vector
  use fsunmatrix_dense_mod          ! Fortran interface to dense SUNMatrix
  use fsunlinsol_dense_mod          ! Fortran interface to dense SUNLinearSolver
  use fsunnonlinsol_newton_mod
  use interface_definitions, only : solout_if
  implicit none
  integer, intent(in)            :: n
  real(R8), intent(inout)        :: t1, t2
  real(R8), intent(inout)        :: var(n)
  integer, intent(out)           :: err
  real(R8), intent(in), optional :: hmax
  procedure(solout_if), optional :: solout
  external :: fcn
  ! specific
  type(N_Vector), pointer :: sunvec_y      ! sundials solution vector
  type(N_Vector), pointer :: sunvec_f      ! sundials solution vector
  type(N_Vector), pointer :: sunvec_av     ! sundials tolerance vector
  type(SUNMatrix), pointer :: sunmat_A      ! sundials matrix
  type(SUNLinearSolver), pointer :: sunlinsol_LS  ! sundials linear solver
  type(SUNNonLinearSolver), pointer :: sunnonlin_NLS ! sundials nonlinear solver
  type(c_ptr)                       :: cvode_mem     ! CVode memory
  type(c_ptr)                       :: sunctx        ! SUNDIALS simulation context
  real(R8)                          :: fval(n), tret(1)
  integer(c_int)                    :: retval
  integer(c_int64_t)                :: neq

  neq = int(n, c_int64_t)

  retval = FSUNContext_Create(SUN_COMM_NULL, sunctx)

  sunvec_y => FN_VMake_Serial(neq, var, sunctx)
  sunvec_f => FN_VMake_Serial(neq, fval, sunctx)
  sunvec_av => FN_VMake_Serial(neq, ATOL, sunctx)

  ! Call FCVodeCreate and FCVodeInit to create and initialize CVode memory
  cvode_mem = FCVodeCreate(CV_BDF, sunctx)
  retval = FCVodeInit(cvode_mem, c_funloc(cvodefcn), t1, sunvec_y)

  ! Set tolerances and maxstep
  retval = FCVodeSVtolerances(cvode_mem, RTOL(1), sunvec_av)
  retval = FCVodeSetMaxNumSteps(cvode_mem, int(IWORK_global(1), c_long))

  ! Create dense SUNMatrix for use in linear solves
  sunmat_A => FSUNDenseMatrix(neq, neq, sunctx)
  ! Create dense SUNLinearSolver object
  sunlinsol_LS => FSUNLinSol_Dense(sunvec_y, sunmat_A, sunctx)
  ! Attach the matrix and linear solver
  retval = FCVodeSetLinearSolver(cvode_mem, sunlinsol_LS, sunmat_A)
  ! Use Newton non-linear solver
  sunnonlin_NLS => FSUNNonlinSol_Newton(sunvec_y, sunctx)
  retval = FCVodeSetNonlinearSolver(cvode_mem, sunnonlin_NLS)
  ! Assign Jacobian
  !retval = FCVodeSetJacFn(cvode_mem, c_funloc(jac))

  if (present(hmax)) &
    retval =  FCVodeSetMaxStep(cvode_mem, hmax)

  retval = FCVode(cvode_mem, t2, sunvec_y, tret(1), CV_NORMAL)
  if (retval /= CV_SUCCESS) then
    write(*,*) 'Error in FCVode, retval = ', retval, '; halting'
    stop 1
  end if

  !call PrintFinalStats(cvode_mem)

  ! free memory
  call FCVodeFree(cvode_mem)
  retval = FSUNLinSolFree(sunlinsol_LS)
  call FSUNMatDestroy(sunmat_A)
  call FN_VDestroy(sunvec_y)
  call FN_VDestroy(sunvec_av)
  retval = FSUNContext_Free(sunctx)

  err = retval


contains

  ! ----------------------------------------------------------------
  ! cvodefcn: The CVODE RHS operator function
  ! ----------------------------------------------------------------
  integer(c_int) function cvodefcn(t, sunvec_y, sunvec_f, user_data) result(ierr) bind(C)
    use, intrinsic :: iso_c_binding
    use fsundials_core_mod
    use fnvector_serial_mod           ! Fortran interface to serial N_Vector
    use fsunmatrix_dense_mod          ! Fortran interface to dense SUNMatrix
    implicit none
    real(c_double), value :: t         ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_f  ! function N_Vector
    type(c_ptr), value :: user_data ! user-defined data
    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(neq) :: yval(:)
    real(c_double), pointer, dimension(neq) :: fval(:)

    ! get data arrays from SUNDIALS vectors
    yval => FN_VGetArrayPointer(sunvec_y)
    fval => FN_VGetArrayPointer(sunvec_f)

    call fcn(NEQ,t,yval,fval)
    ierr = 0
    if (all(fval==-1.0d0)) ierr = 1

  end function cvodefcn
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  ! PrintFinalStats
  !
  ! Print CVode statstics to standard out
  ! ----------------------------------------------------------------
  subroutine PrintFinalStats(cvode_mem)

    !======= Inclusions ===========
    use iso_c_binding
    use fcvodes_mod

    !======= Declarations =========
    implicit none

    type(c_ptr), intent(in) :: cvode_mem ! solver memory structure

    integer(c_int)  :: retval          ! error flag

    integer(c_long) :: nsteps(1)     ! num steps
    integer(c_long) :: nfe(1)        ! num function evals
    integer(c_long) :: netfails(1)   ! num error test fails
    integer(c_long) :: nniters(1)    ! nonlinear solver iterations
    integer(c_long) :: nncfails(1)   ! nonlinear solver fails
    integer(c_long) :: njacevals(1)  ! number of Jacobian evaluations
    integer(c_long) :: nluevals(1)   ! number of LU evals
    integer(c_long) :: ngevals(1)    ! number of root evals

    !======= Internals ============

    retval = FCVodeGetNumSteps(cvode_mem, nsteps)
    if (retval /= 0) then
      print *, 'Error in FCVodeGetNumSteps, retval = ', retval, '; halting'
      stop 1
    end if

    retval = FCVodeGetNumRhsEvals(cvode_mem, nfe)
    if (retval /= 0) then
      print *, 'Error in FCVodeGetNumRhsEvals, retval = ', retval, '; halting'
      stop 1
    end if

    retval = FCVodeGetNumLinSolvSetups(cvode_mem, nluevals)
    if (retval /= 0) then
      print *, 'Error in FCVodeGetNumLinSolvSetups, retval = ', retval, '; halting'
      stop 1
    end if

    retval = FCVodeGetNumErrTestFails(cvode_mem, netfails)
    if (retval /= 0) then
      print *, 'Error in FCVodeGetNumErrTestFails, retval = ', retval, '; halting'
      stop 1
    end if

    retval = FCVodeGetNumNonlinSolvIters(cvode_mem, nniters)
    if (retval /= 0) then
      print *, 'Error in FCVodeGetNumNonlinSolvIters, retval = ', retval, '; halting'
      stop 1
    end if

    retval = FCVodeGetNumNonlinSolvConvFails(cvode_mem, nncfails)
    if (retval /= 0) then
      print *, 'Error in FCVodeGetNumNonlinSolvConvFails, retval = ', retval, '; halting'
      stop 1
    end if

    retval = FCVodeGetNumJacEvals(cvode_mem, njacevals)
    if (retval /= 0) then
      print *, 'Error in FCVodeGetNumJacEvals, retval = ', retval, '; halting'
      stop 1
    end if

    retval = FCVodeGetNumGEvals(cvode_mem, ngevals)
    if (retval /= 0) then
      print *, 'Error in FCVodeGetNumGEvals, retval = ', retval, '; halting'
      stop 1
    end if

    print *, ' '
    print *, ' General Solver Stats:'
    print '(4x,A,i9)', 'Total internal steps taken    =', nsteps
    print '(4x,A,i9)', 'Total rhs function calls      =', nfe
    print '(4x,A,i9)', 'Total Jacobian function calls =', njacevals
    print '(4x,A,i9)', 'Total LU function calls       =', nluevals
    print '(4x,A,i9)', 'Num error test failures       =', netfails
    print '(4x,A,i9)', 'Num nonlinear solver iters    =', nniters
    print '(4x,A,i9)', 'Num nonlinear solver fails    =', nncfails
    print *, ' '

    return

  end subroutine PrintFinalStats

end subroutine wrap_cvode
# endif


subroutine dummy()
endsubroutine


end module oslo
