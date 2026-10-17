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

END MODULE functions


PROGRAM RUNEXAMPLE1
  use oslo
  use functions
#if defined (_OPENMP)
  use omp_lib
#endif
  implicit none
  integer, parameter :: neq = 3, nc = 1000
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

  call setup_odesolver(N=neq,solver='dvodef90',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  !$omp parallel private(i,T,TOUT,err)
  !$omp do schedule(dynamic)
  do i = 1, nc
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,err)
  enddo
  !$omp enddo
  !$omp end parallel
  call cpu_time(time2)
  try = 'fail'
  if (sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))<=1d-20) try = 'success' 
  write(*,Format) 'dvodef90', (time2 - time1)/nthreads, try
  if(try=='fail')write(*,*)'error = ',sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))

  call setup_odesolver(N=neq,solver='H-radau5',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  !$omp parallel private(i,T,TOUT,err)
  !$omp do schedule(dynamic)
  do i = 1, nc
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,err)
  enddo
  !$omp enddo
  !$omp end parallel
  call cpu_time(time2)
  try = 'fail'
  if (sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))<=1d-20) try = 'success' 
  write(*,Format) 'H-radau5', (time2 - time1)/nthreads, try
  if(try=='fail')write(*,*)'error = ',sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))

  call setup_odesolver(N=neq,solver='radau2a',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  !$omp parallel private(i,T,TOUT,err)
  !$omp do schedule(dynamic)
  do i = 1, nc
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,err)
  enddo
  !$omp enddo
  !$omp end parallel
  call cpu_time(time2)
  try = 'fail'
  if (sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))<=1d-20) try = 'success' 
  write(*,Format) 'radau2a', (time2 - time1)/nthreads, try
  if(try=='fail')write(*,*)'error = ',sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))

  call setup_odesolver(N=neq,solver='H-rodas',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  !$omp parallel private(i,T,TOUT,err)
  !$omp do schedule(dynamic)
  do i = 1, nc
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,err)
  enddo
  !$omp enddo
  !$omp end parallel
  call cpu_time(time2)
  try = 'fail'
  if (sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))<=1d-20) try = 'success' 
  write(*,Format) 'H-rodas', (time2 - time1)/nthreads, try
  if(try=='fail')write(*,*)'error = ',sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))

  call setup_odesolver(N=neq,solver='rodas3',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  !$omp parallel private(i,T,TOUT,err)
  !$omp do schedule(dynamic)
  do i = 1, nc
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,err)
  enddo
  !$omp enddo
  !$omp end parallel
  call cpu_time(time2)
  try = 'fail'
  if (sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))<=1d-20) try = 'success' 
  write(*,Format) 'rodas3', (time2 - time1)/nthreads, try
    if(try=='fail')write(*,*)'error = ',sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))

  call setup_odesolver(N=neq,solver='H-sdirk4',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  !$omp parallel private(i,T,TOUT,err)
  !$omp do schedule(dynamic)
  do i = 1, nc
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,err)
  enddo
  !$omp enddo
  !$omp end parallel
  call cpu_time(time2)
  try = 'fail'
  if (sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))<=1d-20) try = 'success' 
  write(*,Format) 'H-sdirk4', (time2 - time1)/nthreads, try
  if(try=='fail')write(*,*)'error = ',sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))

  call setup_odesolver(N=neq,solver='sdirk4b',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  !$omp parallel private(i,T,TOUT,err)
  !$omp do schedule(dynamic)
  do i = 1, nc
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,err)
  enddo
  !$omp enddo
  !$omp end parallel
  call cpu_time(time2)
  try = 'fail'
  if (sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))<=1d-20) try = 'success' 
  write(*,Format) 'sdirk4b', (time2 - time1)/nthreads, try
  if(try=='fail')write(*,*)'error = ',sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))

# if defined(INTEL)
  call setup_odesolver(N=neq,solver='dodesol',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  !$omp parallel private(i,T,TOUT,err)
  !$omp do schedule(dynamic)
  do i = 1, nc
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,err)
  enddo
  !$omp enddo
  !$omp end parallel
  call cpu_time(time2)
  try = 'fail'
  if (sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))<=1d-20) try = 'success' 
  write(*,Format) 'dodesol', (time2 - time1)/nthreads, try
  if(try=='fail')write(*,*)'error = ',sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))
# endif

# if defined(SUNDIALS)
  call setup_odesolver(N=neq,solver='cvode',RT=RT,AT=AT)
  call initialize
  call cpu_time(time1)
  !$omp parallel private(i,T,TOUT,err)
  !$omp do schedule(dynamic)
  do i = 1, nc
    T = 0.0D0; TOUT = tlimit
    call run_odesolver(neq,T,TOUT,Y(:,i),Fgeneral,err)
  enddo
  !$omp enddo
  !$omp end parallel
  call cpu_time(time2)
  try = 'fail'
  if (sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))<=1d-20) try = 'success' 
  write(*,Format) 'cvode', (time2 - time1)/nthreads, try
  if(try=='fail')write(*,*)'error = ',sum(Y(2,2:nc:2))-sum(Y(2,1:nc-1:2))
# endif

contains

  subroutine initialize()
    implicit none
    Y(1,:) = 1.0D0
    Y(2,:) = 0.0D0
    Y(3,:) = 0.0D0
  end subroutine

END PROGRAM RUNEXAMPLE1
