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

END MODULE functions


PROGRAM RUNEXAMPLE1
  use oslo
  use functions
# if defined(MPI)
  use mpi
# endif
  implicit none
  integer, parameter :: neq = 3, nc = 1000
  real(8) :: Y(neq,nc), t, TOUT, tlimit
  real(8) :: time1, time2, RT(neq), AT(neq)
  character(len=30) :: Format
  integer :: i
  character(10) :: try
  integer :: err
  integer :: ierr, rank=0, size, start_idx, end_idx, loc_nc
  real(8) :: diff, ref
  real(8) :: local_sum_even, local_sum_odd, global_sum_even, global_sum_odd

  RT = 1.d-4
  AT = [1.D-8, 1.D-14, 1.D-6]
  tlimit = 4.d+10
  Format = '(A8,E20.8,A15)'

# if defined(MPI)
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
# endif

  loc_nc = nc / size
  start_idx = rank * loc_nc + 1
  end_idx = start_idx + loc_nc - 1

  ! dvodef90
  call setup_odesolver(N=neq,solver='dvodef90',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  do i = start_idx, end_idx
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,err)
  enddo
  call cpu_time(time2)

  local_sum_even = sum(Y(2, start_idx:end_idx:2))
  local_sum_odd  = sum(Y(2, start_idx+1:end_idx:2))
# if defined(MPI)
  call MPI_REDUCE(local_sum_even, global_sum_even, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(local_sum_odd,  global_sum_odd,  1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
# endif
  if (rank == 0) then
    diff = abs(global_sum_even - global_sum_odd)
    ref  = max(abs(global_sum_even), abs(global_sum_odd), 1d-30)
    if (diff / ref <= 1d-12) then
      try = 'success'
    else
      try = 'failure'
    end if
    write(*,Format) 'dvodef90', (time2-time1), trim(try)
  end if
  ! -------------------------

  ! Hairer radau5
  call setup_odesolver(N=neq,solver='H-radau5',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  do i = start_idx, end_idx
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,err)
  enddo
  call cpu_time(time2)

  local_sum_even = sum(Y(2, start_idx:end_idx:2))
  local_sum_odd  = sum(Y(2, start_idx+1:end_idx:2))
# if defined(MPI)
  call MPI_REDUCE(local_sum_even, global_sum_even, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(local_sum_odd,  global_sum_odd,  1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
# endif
  if (rank == 0) then
    diff = abs(global_sum_even - global_sum_odd)
    ref  = max(abs(global_sum_even), abs(global_sum_odd), 1d-30)
    if (diff / ref <= 1d-12) then
      try = 'success'
    else
      try = 'failure'
    end if
    write(*,Format) 'H-radau5', (time2-time1), trim(try)
  end if
  ! -------------------------

  ! FATODE radau2a
  call setup_odesolver(N=neq,solver='radau2a',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  do i = start_idx, end_idx
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,err)
  enddo
  call cpu_time(time2)

  local_sum_even = sum(Y(2, start_idx:end_idx:2))
  local_sum_odd  = sum(Y(2, start_idx+1:end_idx:2))
# if defined(MPI)
  call MPI_REDUCE(local_sum_even, global_sum_even, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(local_sum_odd,  global_sum_odd,  1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
# endif
  if (rank == 0) then
    diff = abs(global_sum_even - global_sum_odd)
    ref  = max(abs(global_sum_even), abs(global_sum_odd), 1d-30)
    if (diff / ref <= 1d-12) then
      try = 'success'
    else
      try = 'failure'
    end if
    write(*,Format) 'radau2a', (time2-time1), trim(try)
  end if
  ! -------------------------

  ! FATODE rodas3
  call setup_odesolver(N=neq,solver='rodas3',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  do i = start_idx, end_idx
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,err)
  enddo
  call cpu_time(time2)

  local_sum_even = sum(Y(2, start_idx:end_idx:2))
  local_sum_odd  = sum(Y(2, start_idx+1:end_idx:2))
# if defined(MPI)
  call MPI_REDUCE(local_sum_even, global_sum_even, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(local_sum_odd,  global_sum_odd,  1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
# endif
  if (rank == 0) then
    diff = abs(global_sum_even - global_sum_odd)
    ref  = max(abs(global_sum_even), abs(global_sum_odd), 1d-30)
    if (diff / ref <= 1d-12) then
      try = 'success'
    else
      try = 'failure'
    end if
    write(*,Format) 'rodas3', (time2-time1), trim(try)
  end if
  ! -------------------------

  ! Hairer sdirk4
  call setup_odesolver(N=neq,solver='H-sdirk4',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  do i = start_idx, end_idx
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,err)
  enddo
  call cpu_time(time2)

  local_sum_even = sum(Y(2, start_idx:end_idx:2))
  local_sum_odd  = sum(Y(2, start_idx+1:end_idx:2))
# if defined(MPI)
  call MPI_REDUCE(local_sum_even, global_sum_even, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(local_sum_odd,  global_sum_odd,  1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
# endif
  if (rank == 0) then
    diff = abs(global_sum_even - global_sum_odd)
    ref  = max(abs(global_sum_even), abs(global_sum_odd), 1d-30)
    if (diff / ref <= 1d-12) then
      try = 'success'
    else
      try = 'failure'
    end if
    write(*,Format) 'H-sdirk4', (time2-time1), trim(try)
  end if
  ! -------------------------

  ! FATODE sdirk4b
  call setup_odesolver(N=neq,solver='sdirk4b',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  do i = start_idx, end_idx
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,err)
  enddo
  call cpu_time(time2)

  local_sum_even = sum(Y(2, start_idx:end_idx:2))
  local_sum_odd  = sum(Y(2, start_idx+1:end_idx:2))
# if defined(MPI)
  call MPI_REDUCE(local_sum_even, global_sum_even, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(local_sum_odd,  global_sum_odd,  1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
# endif
  if (rank == 0) then
    diff = abs(global_sum_even - global_sum_odd)
    ref  = max(abs(global_sum_even), abs(global_sum_odd), 1d-30)
    if (diff / ref <= 1d-12) then
      try = 'success'
    else
      try = 'failure'
    end if
    write(*,Format) 'sdirk4b', (time2-time1), trim(try)
  end if
  ! -------------------------

# if defined(INTEL)
  ! Intel dodesol
  call setup_odesolver(N=neq,solver='dodesol',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  do i = start_idx, end_idx
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,err)
  enddo
  call cpu_time(time2)

  local_sum_even = sum(Y(2, start_idx:end_idx:2))
  local_sum_odd  = sum(Y(2, start_idx+1:end_idx:2))
# if defined(MPI)
  call MPI_REDUCE(local_sum_even, global_sum_even, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(local_sum_odd,  global_sum_odd,  1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
# endif
  if (rank == 0) then
    diff = abs(global_sum_even - global_sum_odd)
    ref  = max(abs(global_sum_even), abs(global_sum_odd), 1d-30)
    if (diff / ref <= 1d-12) then
      try = 'success'
    else
      try = 'failure'
    end if
    write(*,Format) 'dodesol', (time2-time1), trim(try)
  end if
  ! -------------------------
# endif

# if defined(MPI)
  call MPI_FINALIZE(ierr)
# endif


contains

  subroutine initialize()
    implicit none
    Y(1,:) = 1.0D0
    Y(2,:) = 0.0D0
    Y(3,:) = 0.0D0
  end subroutine

END PROGRAM RUNEXAMPLE1









! program hello_mpi
!     use mpi
!     implicit none

!     integer :: ierr, rank, size, name_len
!     character(len=MPI_MAX_PROCESSOR_NAME) :: processor_name

!     call MPI_Init(ierr)

!     call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
!     call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
!     call MPI_Get_processor_name(processor_name, name_len, ierr)

!     print *, 'Hello from rank', rank, 'of', size, 'on', trim(processor_name)

!     call MPI_Finalize(ierr)

! end program hello_mpi
