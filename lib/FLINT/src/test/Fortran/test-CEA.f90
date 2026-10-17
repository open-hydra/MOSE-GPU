program test_eq
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
  use FLINT_CEA_setup
  use FLINT_CEA_solver
  implicit none
  integer :: err
  real(kind=8), dimension(:), allocatable :: rhoi, y_eq
  real(kind=8) :: T_, P_, teq, blessed_Teq, blessed_y
  character(len=16) :: verdict

  p_ = 10d0
  T_ = 1000d0

  !-------------------------------------------------------------------------------------------------
  ! WD
  !-------------------------------------------------------------------------------------------------

  err = read_idealgas_thermo('WD/INPUT')

  blessed_Teq = 4919.064001247616 
  blessed_y = 3.26511387e-01

  allocate(rhoi(ns))
  allocate(y_eq, mold=rhoi)

  rhoi = 1d-20
  rhoi(5) = 0.8d0
  rhoi(1) = 0.2d0

  call CEA_initialize_global()
  call CEA_solve(p_, T_, rhoi, teq, y_eq)

  write(*,*)'FLINT equilibrium temperature   = ', teq
  write(*,*)'Cantera equilibrium temperature = ', blessed_Teq
  write(*,*)'FLINT yCO   = ', y_eq(2)
  write(*,*)'Cantera yCO = ', blessed_y

  verdict = 'success'
  if (isnan(teq)) verdict = 'fail'
  if (abs(teq-blessed_Teq)/blessed_Teq*100d0>1d0) verdict = 'fail'
  if (abs(y_eq(2)-blessed_y)/blessed_Teq*100d0>1d0) verdict = 'fail'

  write(*,'(2A20)') 'WD -> ', verdict

  deallocate(rhoi); deallocate(y_eq)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(species_names); deallocate(elements_names)
  deallocate(species_composition); deallocate(el_weight)

  !-------------------------------------------------------------------------------------------------
  ! ZK
  !-------------------------------------------------------------------------------------------------

  err = read_idealgas_thermo('ZK/INPUT')

  blessed_Teq = 3611.151626300349
  blessed_y = 1.05464743e-01

  allocate(rhoi(ns))
  allocate(y_eq, mold=rhoi)

  rhoi = 1d-20
  rhoi(6) = 0.8d0
  rhoi(17) = 0.2d0

  call CEA_initialize_global()
  call CEA_solve(p_, T_, rhoi, teq, y_eq)

  write(*,*)'FLINT equilibrium temperature   = ', teq
  write(*,*)'Cantera equilibrium temperature = ', blessed_Teq
  write(*,*)'FLINT yCO   = ', y_eq(7)
  write(*,*)'Cantera yCO = ', blessed_y

  verdict = 'success'
  if (isnan(teq)) verdict = 'fail'
  if (abs(teq-blessed_Teq)/blessed_Teq*100d0>1d0) verdict = 'fail'
  if (abs(y_eq(7)-blessed_y)/blessed_Teq*100d0>1d0) verdict = 'fail'

  write(*,'(2A20)') 'ZK -> ', verdict

  deallocate(rhoi); deallocate(y_eq)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(species_names); deallocate(elements_names)
  deallocate(species_composition); deallocate(el_weight)

  !-------------------------------------------------------------------------------------------------
  ! TSR-GP-24
  !-------------------------------------------------------------------------------------------------

  err = read_idealgas_thermo('TSR-GP-24/INPUT')

  blessed_Teq = 3616.638618717671
  blessed_y = 1.00013716e-01

  allocate(rhoi(ns))
  allocate(y_eq, mold=rhoi)

  rhoi = 1d-20
  rhoi(19) = 0.8d0
  rhoi(22) = 0.2d0

  call CEA_initialize_global()
  call CEA_solve(p_, T_, rhoi, teq, y_eq)

  write(*,*)'FLINT equilibrium temperature   = ', teq
  write(*,*)'Cantera equilibrium temperature = ', blessed_Teq
  write(*,*)'FLINT yCO   = ', y_eq(4)
  write(*,*)'Cantera yCO = ', blessed_y

  verdict = 'success'
  if (isnan(teq)) verdict = 'fail'
  if (abs(teq-blessed_Teq)/blessed_Teq*100d0>1d0) verdict = 'fail'
  if (abs(y_eq(4)-blessed_y)/blessed_Teq*100d0>1d0) verdict = 'fail'

  write(*,'(2A20)') 'TSR-GP-24 -> ', verdict

  deallocate(rhoi); deallocate(y_eq)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(species_names); deallocate(elements_names)
  deallocate(species_composition); deallocate(el_weight)

  !-------------------------------------------------------------------------------------------------
  ! Ecker
  !-------------------------------------------------------------------------------------------------

  err = read_idealgas_thermo('Ecker/INPUT')

  blessed_Teq = 1571.8416518627969
  blessed_y = 1.93015081e-01

  allocate(rhoi(ns))
  allocate(y_eq, mold=rhoi)

  rhoi = 1d-20
  rhoi(12) = 0.18798856d0
  rhoi(1) = 0.00534534d0
  rhoi(14) = 0.8066661d0

  call CEA_initialize_global()
  call CEA_solve(p_, T_, rhoi, teq, y_eq)

  write(*,*)'FLINT equilibrium temperature   = ', teq
  write(*,*)'Cantera equilibrium temperature = ', blessed_Teq
  write(*,*)'FLINT yOH   = ', y_eq(11)
  write(*,*)'Cantera yOH = ', blessed_y

  verdict = 'success'
  if (isnan(teq)) verdict = 'fail'
  if (abs(teq-blessed_Teq)/blessed_Teq*100d0>1d0) verdict = 'fail'
  if (abs(y_eq(11)-blessed_y)/blessed_Teq*100d0>1d0) verdict = 'fail'

  write(*,'(2A20)') 'Ecker -> ', verdict

  deallocate(rhoi); deallocate(y_eq)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(species_names); deallocate(elements_names)
  deallocate(species_composition); deallocate(el_weight)

end program test_eq
