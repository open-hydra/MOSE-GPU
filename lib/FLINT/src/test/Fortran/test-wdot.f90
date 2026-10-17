! Test for wdot with thrid-body reactions
program test
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
  use FLINT_Load_chemistry
  use FLINT_Lib_Chemistry_data
  use FLINT_Lib_Chemistry_wdot
# if defined (CANTERA)
  use FLINT_Lib_Chemistry_rhs, only: gas
  use cantera
  use FLINT_cantera_load
# endif
  implicit none
  integer, parameter :: Tend=2000, Tstart=100
  real(8) :: T, rho, R
  real(8), allocatable :: droic(:), rhoi(:)
  real(8), allocatable :: wdot_explicit(:,:), wdot_cantera(:,:)
  real(8) :: time1, time2
  integer :: i, j, err
  character(32) :: mech_name

!   !-------------------------------------------------------------------------------------------------
!   ! WD
!   !-------------------------------------------------------------------------------------------------

!   call execute_command_line('mkdir -p WD/OUTPUT')
!   err = read_idealgas_thermo('WD/INPUT')
!   err = read_chemistry( folder='WD/INPUT', mech_name=mech_name )
! # if defined(CANTERA)
!   call load_phase(gas, 'WD/INPUT/WD.yaml')
! # endif
!   call Assign_Mechanism(mech_name)

!   allocate(rhoi(1:ns))
!   allocate(droic(1:ns))
!   allocate(wdot_explicit(ns,Tstart:Tend))
!   allocate(wdot_cantera(ns,Tstart:Tend))

!   rhoi = 1d-20
!   rhoi(1) = 0.2
!   rhoi(2) = 0.8
!   R = f_Rtot(rhoi)
!   rho = sum(rhoi)

!   call cpu_time(time1)
!   do j = 1, 1
!     do i = Tstart, Tend
!       T = dble(i)
!       call Chemistry_Source ( rhoi, T, droic )
!       wdot_explicit(:,i) = droic
!     enddo
!   enddo
!   call cpu_time(time2)

!   write(*,*) 'WD explicit time =', time2-time1

! # if defined(CANTERA)

!   call setState_TRY(gas, T, rho, rhoi)

!   call cpu_time(time1)
!   do j = 1, 1
!     do i = Tstart, Tend
!       T = dble(i)
!       call setTemperature(gas, T)
!       call getNetProductionRates(gas, droic)
!       droic = droic*wm_tab
!       wdot_cantera(:,i) = droic
!     enddo
!   enddo
!   call cpu_time(time2)

!   write(*,*) 'WD Cantera time =', time2-time1

! # endif

!   open(100, file='WD/OUTPUT/wdot-explicit.dat', status='replace', form='formatted')
!   do i = Tstart, Tend
!     write(100,*) dble(i), (wdot_explicit(j,i),j=1,ns)
!   enddo
!   close(100)
! # if defined(CANTERA)
!   open(200, file='WD/OUTPUT/wdot-cantera.dat', status='replace', form='formatted')
!   do i = Tstart, Tend
!     write(200,*) dble(i), (wdot_cantera(j,i),j=1,ns)
!   enddo
!   close(200)
! # endif

!   deallocate(wdot_cantera); deallocate(wdot_explicit)
!   deallocate(droic)
!   deallocate(rhoi)
!   deallocate(wm_tab); deallocate(Ri_tab)
!   deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
!   call free_chemistry_data()

!   !-------------------------------------------------------------------------------------------------
!   ! TROYES
!   !-------------------------------------------------------------------------------------------------

!   call execute_command_line('mkdir -p Troyes/OUTPUT')
!   err = read_idealgas_thermo('Troyes/INPUT')
!   err = read_chemistry( folder='Troyes/INPUT', mech_name=mech_name )
! # if defined (CANTERA)
!   call load_phase(gas, 'Troyes/INPUT/troyes.yaml')
! # endif
!   call Assign_Mechanism(mech_name)

!   allocate(rhoi(1:ns))
!   allocate(droic(1:ns))
!   allocate(wdot_explicit(ns,Tstart:Tend))
!   allocate(wdot_cantera(ns,Tstart:Tend))

!   rhoi = 1d0/ns
!   R = f_Rtot(rhoi)
!   rho = sum(rhoi)

!   call cpu_time(time1)
!   do j = 1, 1
!     do i = Tstart, Tend
!       T = dble(i)
!       call Chemistry_Source ( rhoi, T, droic )
!       wdot_explicit(:,i) = droic
!     enddo
!   enddo
!   call cpu_time(time2)

!   write(*,*) 'Troyes explicit time =', time2-time1

! # if defined (CANTERA)

!   call setState_TRY(gas, T, rho, rhoi)

!   call cpu_time(time1)
!   do j = 1, 1
!     do i = Tstart, Tend
!       T = dble(i)
!       call setTemperature(gas, T)
!       call getNetProductionRates(gas, droic)
!       droic = droic*wm_tab
!       wdot_cantera(:,i) = droic
!     enddo
!   enddo
!   call cpu_time(time2)

!   write(*,*) 'Troyes Cantera time =', time2-time1

! # endif

!   open(100, file='Troyes/OUTPUT/wdot-explicit.dat', status='replace', form='formatted')
!   do i = Tstart, Tend
!     write(100,*) dble(i), (wdot_explicit(j,i),j=1,ns)
!   enddo
!   close(100)
! # if defined (CANTERA)
!   open(200, file='Troyes/OUTPUT/wdot-cantera.dat', status='replace', form='formatted')
!   do i = Tstart, Tend
!     write(200,*) dble(i), (wdot_cantera(j,i),j=1,ns)
!   enddo
!   close(200)
! # endif

!   deallocate(wdot_cantera); deallocate(wdot_explicit)
!   deallocate(droic)
!   deallocate(rhoi)
!   deallocate(wm_tab); deallocate(Ri_tab)
!   deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
!   call free_chemistry_data()

!   !-------------------------------------------------------------------------------------------------
!   ! ECKER
!   !-------------------------------------------------------------------------------------------------

!   call execute_command_line('mkdir -p Ecker/OUTPUT')
!   err = read_idealgas_thermo('Ecker/INPUT')
!   err = read_chemistry( folder='Ecker/INPUT', mech_name=mech_name )
! # if defined (CANTERA)
!   call load_phase(gas, 'Ecker/INPUT/ecker.yaml')
! # endif
!   call Assign_Mechanism(mech_name)

!   allocate(rhoi(1:ns))
!   allocate(droic(1:ns))
!   allocate(wdot_explicit(ns,Tstart:Tend))
!   allocate(wdot_cantera(ns,Tstart:Tend))

!   rhoi = 1d-20
!   rhoi(2) = 0.00606
!   rhoi(7) = 0.00365
!   rhoi(9) = 0.00861
!   rhoi(11) = 0.00025
!   rhoi(14) = 0.98143
!   R = f_Rtot(rhoi)
!   rho = sum(rhoi)

!   call cpu_time(time1)
!   do j = 1, 1
!     do i = Tstart, Tend
!       T = dble(i)
!       call Chemistry_Source ( rhoi, T, droic )
!       wdot_explicit(:,i) = droic
!     enddo
!   enddo
!   call cpu_time(time2)

!   write(*,*) 'Ecker explicit time =', time2-time1

! # if defined (CANTERA)

!   call setState_TRY(gas, T, rho, rhoi)

!   call cpu_time(time1)
!   do j = 1, 1
!     do i = Tstart, Tend
!       T = dble(i)
!       call setTemperature(gas, T)
!       call getNetProductionRates(gas, droic)
!       droic = droic*wm_tab
!       wdot_cantera(:,i) = droic
!     enddo
!   enddo
!   call cpu_time(time2)

!   write(*,*) 'Ecker Cantera time =', time2-time1

! # endif

!   open(100, file='Ecker/OUTPUT/wdot-explicit.dat', status='replace', form='formatted')
!   do i = Tstart, Tend
!     write(100,*) dble(i), (wdot_explicit(j,i),j=1,ns)
!   enddo
!   close(100)
! # if defined (CANTERA)
!   open(200, file='Ecker/OUTPUT/wdot-cantera.dat', status='replace', form='formatted')
!   do i = Tstart, Tend
!     write(200,*) dble(i), (wdot_cantera(j,i),j=1,ns)
!   enddo
!   close(200)
! # endif

!   deallocate(wdot_cantera); deallocate(wdot_explicit)
!   deallocate(droic)
!   deallocate(rhoi)
!   deallocate(wm_tab); deallocate(Ri_tab)
!   deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
!   call free_chemistry_data()

!   !-------------------------------------------------------------------------------------------------
!   ! CROSS
!   !-------------------------------------------------------------------------------------------------

!   call execute_command_line('mkdir -p Cross/OUTPUT')
!   err = read_idealgas_thermo('Cross/INPUT')
!   err = read_chemistry( folder='Cross/INPUT', mech_name=mech_name )
! # if defined (CANTERA)
!   call load_phase(gas, 'Cross/INPUT/cross.yaml')
! # endif
!   call Assign_Mechanism(mech_name)

!   allocate(rhoi(1:ns))
!   allocate(droic(1:ns))
!   allocate(wdot_explicit(ns,Tstart:Tend))
!   allocate(wdot_cantera(ns,Tstart:Tend))

!   rhoi = 1d0/ns
!   R = f_Rtot(rhoi)
!   rho = sum(rhoi)

!   call cpu_time(time1)
!   do j = 1, 1
!     do i = Tstart, Tend
!       T = dble(i)
!       call Chemistry_Source ( rhoi, T, droic )
!       wdot_explicit(:,i) = droic
!     enddo
!   enddo
!   call cpu_time(time2)

!   write(*,*) 'Cross explicit time =', time2-time1

! # if defined (CANTERA)

!   call setState_TRY(gas, T, rho, rhoi)

!   call cpu_time(time1)
!   do j = 1, 1
!     do i = Tstart, Tend
!       T = dble(i)
!       call setTemperature(gas, T)
!       call getNetProductionRates(gas, droic)
!       droic = droic*wm_tab
!       wdot_cantera(:,i) = droic
!     enddo
!   enddo
!   call cpu_time(time2)

!   write(*,*) 'Cross Cantera time =', time2-time1

! # endif

!   open(100, file='Cross/OUTPUT/wdot-explicit.dat', status='replace', form='formatted')
!   do i = Tstart, Tend
!     write(100,*) dble(i), (wdot_explicit(j,i),j=1,ns)
!   enddo
!   close(100)
! # if defined (CANTERA)
!   open(200, file='Cross/OUTPUT/wdot-cantera.dat', status='replace', form='formatted')
!   do i = Tstart, Tend
!     write(200,*) dble(i), (wdot_cantera(j,i),j=1,ns)
!   enddo
!   close(200)
! # endif

!   deallocate(wdot_cantera); deallocate(wdot_explicit)
!   deallocate(droic)
!   deallocate(rhoi)
!   deallocate(wm_tab); deallocate(Ri_tab)
!   deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
!   call free_chemistry_data()

!   !-------------------------------------------------------------------------------------------------
!   ! SMOOKE
!   !-------------------------------------------------------------------------------------------------

!   call execute_command_line('mkdir -p Smooke/OUTPUT')
!   err = read_idealgas_thermo('Smooke/INPUT')
!   err = read_chemistry( folder='Smooke/INPUT', mech_name=mech_name )
! # if defined (CANTERA)
!   call load_phase(gas, 'Smooke/INPUT/smooke.yaml')
! # endif
!   call Assign_Mechanism(mech_name)

!   allocate(rhoi(1:ns))
!   allocate(droic(1:ns))
!   allocate(wdot_explicit(ns,Tstart:Tend))
!   allocate(wdot_cantera(ns,Tstart:Tend))

!   rhoi = 1d0/ns
!   R = f_Rtot(rhoi)
!   rho = sum(rhoi)

!   call cpu_time(time1)
!   do j = 1, 1
!     do i = Tstart, Tend
!       T = dble(i)
!       call Chemistry_Source ( rhoi, T, droic )
!       wdot_explicit(:,i) = droic
!     enddo
!   enddo
!   call cpu_time(time2)

!   write(*,*) 'Smooke explicit time =', time2-time1

! # if defined (CANTERA)

!   call setState_TRY(gas, T, rho, rhoi)

!   call cpu_time(time1)
!   do j = 1, 1
!     do i = Tstart, Tend
!       T = dble(i)
!       call setTemperature(gas, T)
!       call getNetProductionRates(gas, droic)
!       droic = droic*wm_tab
!       wdot_cantera(:,i) = droic
!     enddo
!   enddo
!   call cpu_time(time2)

!   write(*,*) 'Smooke Cantera time =', time2-time1

! # endif

!   open(100, file='Smooke/OUTPUT/wdot-explicit.dat', status='replace', form='formatted')
!   do i = Tstart, Tend
!     write(100,*) dble(i), (wdot_explicit(j,i),j=1,ns)
!   enddo
!   close(100)
! # if defined (CANTERA)
!   open(200, file='Smooke/OUTPUT/wdot-cantera.dat', status='replace', form='formatted')
!   do i = Tstart, Tend
!     write(200,*) dble(i), (wdot_cantera(j,i),j=1,ns)
!   enddo
!   close(200)
! # endif

!   deallocate(wdot_cantera); deallocate(wdot_explicit)
!   deallocate(droic)
!   deallocate(rhoi)
!   deallocate(wm_tab); deallocate(Ri_tab)
!   deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
!   call free_chemistry_data()

!   !-------------------------------------------------------------------------------------------------
!   ! CORIA-CNRS
!   !-------------------------------------------------------------------------------------------------

!   call execute_command_line('mkdir -p CORIA/OUTPUT')
!   err = read_idealgas_thermo('CORIA/INPUT')
!   err = read_chemistry( folder='CORIA/INPUT', mech_name=mech_name )
! # if defined (CANTERA)
!   call load_phase(gas, 'CORIA/INPUT/coria.yaml')
! # endif
!   call Assign_Mechanism(mech_name)

!   allocate(rhoi(1:ns))
!   allocate(droic(1:ns))
!   allocate(wdot_explicit(ns,Tstart:Tend))
!   allocate(wdot_cantera(ns,Tstart:Tend))

!   rhoi = 1d0/ns
!   R = f_Rtot(rhoi)
!   rho = sum(rhoi)

!   call cpu_time(time1)
!   do j = 1, 1
!     do i = Tstart, Tend
!       T = dble(i)
!       call Chemistry_Source ( rhoi, T, droic )
!       wdot_explicit(:,i) = droic
!     enddo
!   enddo
!   call cpu_time(time2)

!   write(*,*) 'CORIA explicit time =', time2-time1

! # if defined (CANTERA)

!   call setState_TRY(gas, T, rho, rhoi)

!   call cpu_time(time1)
!   do j = 1, 1
!     do i = Tstart, Tend
!       T = dble(i)
!       call setTemperature(gas, T)
!       call getNetProductionRates(gas, droic)
!       droic = droic*wm_tab
!       wdot_cantera(:,i) = droic
!     enddo
!   enddo
!   call cpu_time(time2)

!   write(*,*) 'CORIA Cantera time =', time2-time1

! # endif

!   open(100, file='CORIA/OUTPUT/wdot-explicit.dat', status='replace', form='formatted')
!   do i = Tstart, Tend
!     write(100,*) dble(i), (wdot_explicit(j,i),j=1,ns)
!   enddo
!   close(100)
! # if defined (CANTERA)
!   open(200, file='CORIA/OUTPUT/wdot-cantera.dat', status='replace', form='formatted')
!   do i = Tstart, Tend
!     write(200,*) dble(i), (wdot_cantera(j,i),j=1,ns)
!   enddo
!   close(200)
! # endif

!   deallocate(wdot_cantera); deallocate(wdot_explicit)
!   deallocate(droic)
!   deallocate(rhoi)
!   deallocate(wm_tab); deallocate(Ri_tab)
!   deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
!   call free_chemistry_data()

!   !-------------------------------------------------------------------------------------------------
!   ! TSR-CDF-13
!   !-------------------------------------------------------------------------------------------------

!   call execute_command_line('mkdir -p TSR-CDF-13/OUTPUT')
!   err = read_idealgas_thermo('TSR-CDF-13/INPUT')
!   err = read_chemistry( folder='TSR-CDF-13/INPUT', mech_name=mech_name )
! # if defined (CANTERA)
!   call load_phase(gas, 'TSR-CDF-13/INPUT/TSR-CDF-13.yaml')
! # endif
!   call Assign_Mechanism(mech_name)

!   allocate(rhoi(1:ns))
!   allocate(droic(1:ns))
!   allocate(wdot_explicit(ns,Tstart:Tend))
!   allocate(wdot_cantera(ns,Tstart:Tend))

!   rhoi = 1d0/ns
!   R = f_Rtot(rhoi)
!   rho = sum(rhoi)

!   call cpu_time(time1)
!   do j = 1, 1
!     do i = Tstart, Tend
!       T = dble(i)
!       call Chemistry_Source ( rhoi, T, droic )
!       wdot_explicit(:,i) = droic
!     enddo
!   enddo
!   call cpu_time(time2)

!   write(*,*) 'TSR-CDF-13 explicit time =', time2-time1

! # if defined (CANTERA)

!   call setState_TRY(gas, T, rho, rhoi)

!   call cpu_time(time1)
!   do j = 1, 1
!     do i = Tstart, Tend
!       T = dble(i)
!       call setTemperature(gas, T)
!       call getNetProductionRates(gas, droic)
!       droic = droic*wm_tab
!       wdot_cantera(:,i) = droic
!     enddo
!   enddo
!   call cpu_time(time2)

!   write(*,*) 'TSR-CDF-13 Cantera time =', time2-time1

! # endif

!   open(100, file='TSR-CDF-13/OUTPUT/wdot-explicit.dat', status='replace', form='formatted')
!   do i = Tstart, Tend
!     write(100,*) dble(i), (wdot_explicit(j,i),j=1,ns)
!   enddo
!   close(100)
! # if defined (CANTERA)
!   open(200, file='TSR-CDF-13/OUTPUT/wdot-cantera.dat', status='replace', form='formatted')
!   do i = Tstart, Tend
!     write(200,*) dble(i), (wdot_cantera(j,i),j=1,ns)
!   enddo
!   close(200)
! # endif

!   deallocate(wdot_cantera); deallocate(wdot_explicit)
!   deallocate(droic)
!   deallocate(rhoi)
!   deallocate(wm_tab); deallocate(Ri_tab)
!   deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
!   call free_chemistry_data()

!   !-------------------------------------------------------------------------------------------------
!   ! TSR-GP-24
!   !-------------------------------------------------------------------------------------------------

!   call execute_command_line('mkdir -p TSR-GP-24/OUTPUT')
!   err = read_idealgas_thermo('TSR-GP-24/INPUT')
!   err = read_chemistry( folder='TSR-GP-24/INPUT', mech_name=mech_name )
! # if defined (CANTERA)
!   call load_phase(gas, 'TSR-GP-24/INPUT/TSR-GP-24.yaml')
! # endif
!   call Assign_Mechanism(mech_name)

!   allocate(rhoi(1:ns))
!   allocate(droic(1:ns))
!   allocate(wdot_explicit(ns,Tstart:Tend))
!   allocate(wdot_cantera(ns,Tstart:Tend))

!   rhoi = 1d0/ns
!   R = f_Rtot(rhoi)
!   rho = sum(rhoi)

!   call cpu_time(time1)
!   do j = 1, 1
!     do i = Tstart, Tend
!       T = dble(i)
!       call Chemistry_Source ( rhoi, T, droic )
!       wdot_explicit(:,i) = droic
!     enddo
!   enddo
!   call cpu_time(time2)

!   write(*,*) 'TSR-GP-24 explicit time =', time2-time1

! # if defined (CANTERA)

!   call setState_TRY(gas, T, rho, rhoi)

!   call cpu_time(time1)
!   do j = 1, 1
!     do i = Tstart, Tend
!       T = dble(i)
!       call setTemperature(gas, T)
!       call getNetProductionRates(gas, droic)
!       droic = droic*wm_tab
!       wdot_cantera(:,i) = droic
!     enddo
!   enddo
!   call cpu_time(time2)

!   write(*,*) 'TSR-GP-24 Cantera time =', time2-time1

! # endif

!   open(100, file='TSR-GP-24/OUTPUT/wdot-explicit.dat', status='replace', form='formatted')
!   do i = Tstart, Tend
!     write(100,*) dble(i), (wdot_explicit(j,i),j=1,ns)
!   enddo
!   close(100)
! # if defined (CANTERA)
!   open(200, file='TSR-GP-24/OUTPUT/wdot-cantera.dat', status='replace', form='formatted')
!   do i = Tstart, Tend
!     write(200,*) dble(i), (wdot_cantera(j,i),j=1,ns)
!   enddo
!   close(200)
! # endif

!   deallocate(wdot_cantera); deallocate(wdot_explicit)
!   deallocate(droic)
!   deallocate(rhoi)
!   deallocate(wm_tab); deallocate(Ri_tab)
!   deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
!   call free_chemistry_data()

!   !-------------------------------------------------------------------------------------------------
!   ! TSR-Rich-31
!   !-------------------------------------------------------------------------------------------------

!   call execute_command_line('mkdir -p TSR-Rich-31/OUTPUT')
!   err = read_idealgas_thermo('TSR-Rich-31/INPUT')
!   err = read_chemistry( folder='TSR-Rich-31/INPUT', mech_name=mech_name )
! # if defined (CANTERA)
!   call load_phase(gas, 'TSR-Rich-31/INPUT/TSR-Rich-31.yaml')
! # endif
!   call Assign_Mechanism(mech_name)

!   allocate(rhoi(1:ns))
!   allocate(droic(1:ns))
!   allocate(wdot_explicit(ns,Tstart:Tend))
!   allocate(wdot_cantera(ns,Tstart:Tend))

!   rhoi = 1d0/ns
!   R = f_Rtot(rhoi)
!   rho = sum(rhoi)

!   call cpu_time(time1)
!   do j = 1, 1
!     do i = Tstart, Tend
!       T = dble(i)
!       call Chemistry_Source ( rhoi, T, droic )
!       wdot_explicit(:,i) = droic
!     enddo
!   enddo
!   call cpu_time(time2)

!   write(*,*) 'TSR-Rich-31 explicit time =', time2-time1

! # if defined (CANTERA)

!   call setState_TRY(gas, T, rho, rhoi)

!   call cpu_time(time1)
!   do j = 1, 1
!     do i = Tstart, Tend
!       T = dble(i)
!       call setTemperature(gas, T)
!       call getNetProductionRates(gas, droic)
!       droic = droic*wm_tab
!       wdot_cantera(:,i) = droic
!     enddo
!   enddo
!   call cpu_time(time2)

!   write(*,*) 'TSR-Rich-31 Cantera time =', time2-time1

! # endif

!   open(100, file='TSR-Rich-31/OUTPUT/wdot-explicit.dat', status='replace', form='formatted')
!   do i = Tstart, Tend
!     write(100,*) dble(i), (wdot_explicit(j,i),j=1,ns)
!   enddo
!   close(100)
! # if defined (CANTERA)
!   open(200, file='TSR-Rich-31/OUTPUT/wdot-cantera.dat', status='replace', form='formatted')
!   do i = Tstart, Tend
!     write(200,*) dble(i), (wdot_cantera(j,i),j=1,ns)
!   enddo
!   close(200)
! # endif

!   deallocate(wdot_cantera); deallocate(wdot_explicit)
!   deallocate(droic)
!   deallocate(rhoi)
!   deallocate(wm_tab); deallocate(Ri_tab)
!   deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
!   call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! Pelucchi
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p Pelucchi/OUTPUT')
  err = read_idealgas_thermo('Pelucchi/INPUT')
  err = read_chemistry( folder='Pelucchi/INPUT', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, 'Pelucchi/INPUT/pelucchi.yaml')
# endif
  call Assign_Mechanism(mech_name)

  allocate(rhoi(1:ns))
  allocate(droic(1:ns))
  allocate(wdot_explicit(ns,Tstart:Tend))
  allocate(wdot_cantera(ns,Tstart:Tend))

  rhoi = 1d0/ns
  R = f_Rtot(rhoi)
  rho = sum(rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      call Chemistry_Source ( rhoi, T, droic )
      wdot_explicit(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'Pelucchi explicit time =', time2-time1

# if defined (CANTERA)

  call setState_TRY(gas, T, rho, rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      call setTemperature(gas, T)
      call getNetProductionRates(gas, droic)
      droic = droic*wm_tab
      wdot_cantera(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'Pelucchi Cantera time =', time2-time1

# endif

  open(100, file='Pelucchi/OUTPUT/wdot-explicit.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(100,*) dble(i), (wdot_explicit(j,i),j=1,ns)
  enddo
  close(100)
# if defined (CANTERA)
  open(200, file='Pelucchi/OUTPUT/wdot-cantera.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(200,*) dble(i), (wdot_cantera(j,i),j=1,ns)
  enddo
  close(200)
# endif

  deallocate(wdot_cantera); deallocate(wdot_explicit)
  deallocate(droic)
  deallocate(rhoi)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  call free_chemistry_data()

end program test