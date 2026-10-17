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
  implicit none
  integer, parameter :: neq = 3
  real(8) :: Y(neq), t, TOUT, tlimit
  real(8) :: time1, time2, RT(neq), AT(neq)
  character(len=30) :: Format
  integer :: err, i, ntimes=10000

  RT = 1.d-4
  AT = [1.D-8, 1.D-14, 1.D-6]
  tlimit = 4.d+10
  Format = '(A8,4E20.8)'

  call setup_odesolver(N=neq,solver='dvodef90',RT=RT,AT=AT)
  call cpu_time(time1)
  do i = 1, ntimes
  call initialize
  call run_odesolver(neq,T,TOUT,Y,Fgeneral,err)
  end do
  call cpu_time(time2)
  write(*,Format) 'dvodef90', time2-time1, Y(:)

  call setup_odesolver(N=neq,solver='H-radau5',RT=RT,AT=AT)
  call cpu_time(time1)
  do i = 1, ntimes
  call initialize
  call run_odesolver(neq,T,TOUT,Y,Fgeneral,err)
  end do
  call cpu_time(time2)
  write(*,Format) 'H-radau5', time2-time1, Y(:)

  call setup_odesolver(N=neq,solver='radau2a',RT=RT,AT=AT)
  call cpu_time(time1)
  do i = 1, ntimes
  call initialize
  call run_odesolver(neq,T,TOUT,Y,Fgeneral,err)
  end do
  call cpu_time(time2)
  write(*,Format) 'radau2a', time2-time1, Y(:)

  call setup_odesolver(N=neq,solver='H-rodas',RT=RT,AT=AT)
  call cpu_time(time1)
  do i = 1, ntimes
  call initialize
  call run_odesolver(neq,T,TOUT,Y,Fgeneral,err)
  end do
  call cpu_time(time2)
  write(*,Format) 'H-rodas', time2-time1, Y(:)

  call setup_odesolver(N=neq,solver='rodas3',RT=RT,AT=AT)
  call cpu_time(time1)
  do i = 1, ntimes
  call initialize
  call run_odesolver(neq,T,TOUT,Y,Fgeneral,err)
  end do
  call cpu_time(time2)
  write(*,Format) 'rodas3', time2-time1, Y(:)

  call setup_odesolver(N=neq,solver='H-sdirk4',RT=RT,AT=AT)
  call cpu_time(time1)
  do i = 1, ntimes
  call initialize
  call run_odesolver(neq,T,TOUT,Y,Fgeneral,err)
  end do
  call cpu_time(time2)
  write(*,Format) 'H-sdirk4', time2-time1, Y(:)

  call setup_odesolver(N=neq,solver='sdirk4b',RT=RT,AT=AT)
  call cpu_time(time1)
  do i = 1, ntimes
  call initialize
  call run_odesolver(neq,T,TOUT,Y,Fgeneral,err)
  end do
  call cpu_time(time2)
  write(*,Format) 'sdirk4b', time2-time1, Y(:)

# if defined(INTEL)
  call setup_odesolver(N=neq,solver='dodesol',RT=RT,AT=AT)
  call cpu_time(time1)
  do i = 1, ntimes
  call initialize
  call run_odesolver(neq,T,TOUT,Y,Fgeneral,err)
  end do
  call cpu_time(time2)
  write(*,Format) 'dodesol', time2-time1, Y(:)
# endif

# if defined(SUNDIALS)
  call setup_odesolver(N=neq,solver='cvode',RT=RT,AT=AT)
  call cpu_time(time1)
  do i = 1, ntimes
  call initialize
  call run_odesolver(neq,T,TOUT,Y,Fgeneral,err)
  end do
  call cpu_time(time2)
  write(*,Format) 'cvode', time2-time1, Y(:)
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
