!> \brief Test program to compare native and Cantera thermodynamic property calculations.
!>
!> This program benchmarks and compares the computation of specific heat at constant pressure (cp)
!> using both a native implementation and the Cantera library for a given gas mixture.
!>
!> - Reads ideal gas properties and loads a phase from 'WD.yaml' using Cantera.
!> - Initializes species mass fractions for CH4, O2, CO2, and CO.
!> - Computes mixture density and gas constant.
!> - Benchmarks the native cp calculation over a range of temperatures.
!> - Benchmarks the Cantera cp calculation over the same temperature range.
!> - Outputs the elapsed time for each method and the relative difference in cp values.
!>
!> \author Marco Grossi
!> \date 2025
program test
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
# if defined (CANTERA)
  use cantera
  use FLINT_cantera_load
# endif
  implicit none
# if defined (CANTERA)
  type(phase_t) :: gas
# endif
  real(8) :: rhoi(9), T, rho, R
  real(8) :: cp_native, cp_cantera
  real(8) :: lam_cantera, lam_native, mi_native, nu_cantera
  real(8) :: time1, time2
  integer :: i, j, err

  err = read_idealgas_thermo()
# if defined (CANTERA)
  call load_phase(gas, 'WD/INPUT/WD.yaml')
# endif

  rhoi = 0.d0
  rhoi(1) = 0.25d0 ! CH4
  rhoi(2) = 0.15d0 ! O2
  rhoi(3) = 0.30d0 ! CO2
  rhoi(5) = 0.40d0 ! CO
  R = f_Rtot(rhoi)
  rho = sum(rhoi)

  call cpu_time(time1)
  do j = 1, 1000
    do i = 1,2000
      T = dble(i)
      cp_native = f_cp(rhoi,T,rho)
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) time2-time1

# if defined (CANTERA)

  call setState_TRY(gas, T, rho, rhoi)

  call cpu_time(time1)
  do j = 1, 1000
    do i = 1,2000
      T = dble(i)
      call setTemperature(gas, T)
      cp_cantera = cp_mass(gas)
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) time2-time1

  write(*,*) abs(cp_cantera-cp_native)/cp_cantera*100d0

# endif

end program test