!> gpu-port-mcdiff unit test: host co_DS_expr against device co_DS_expr_dev on the same states.
!> usage: test_co_DS_dev <dir with phase.txt, thermo.dat, diffusion.dat> <states.bin> <out.bin>
!> states.bin (stream): int32 n; real64 rhoi(7,n), T(n), p(n).
!> Both routines get the same rho, (Tint, Tdiff) and p for each state, so the test isolates the two implementations.
!> The tables reach the device through flint_acc_upload_thermo, the production upload.
!> out.bin (stream): int32 n; real64 Dh(7,n), Dd(7,n). stdout: per species, max relative difference and non-identical count.
program test_co_DS_dev
#ifdef _OPENACC
  use openacc
#endif
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
  use FLINT_Lib_ThermoTransport_dev, only: co_DS_expr_dev
  use FLINT_Lib_Radau5_dev, only: flint_acc_upload_thermo
  implicit none
  integer, parameter :: NSP = 7
  character(512) :: dir, fstates, fout
  integer :: n, i, s, ios, u, Ti(2), ndiff(NSP), nnan_h, nnan_d
  real(8), allocatable :: rhoi(:,:), T(:), p(:), rho(:), Td(:), Dh(:,:), Dd(:,:)
  integer, allocatable :: Tl(:)
  real(8) :: rl(NSP), dl(NSP), relmax(NSP), rel

  call get_command_argument(1, dir); call get_command_argument(2, fstates); call get_command_argument(3, fout)
#ifdef _OPENACC
  write(*,'(A,I0)') ' acc_get_num_devices(acc_device_nvidia) = ', acc_get_num_devices(acc_device_nvidia)
#else
  write(*,'(A)') ' HOST-ONLY BUILD: co_DS_expr_dev runs on the CPU'
#endif
  ios = read_idealgas_thermo(trim(dir))
  if (ios /= 0) error stop 'read_idealgas_thermo failed'
  if (ns /= NSP) error stop 'ns /= 7'
  ios = read_idealgas_diffusion(trim(dir))
  write(*,'(A,I0,A,I0,A,I0,A,I0,A,ES14.6)') ' read_idealgas_diffusion ios ', ios, '  ndij ', ndij, '  Tmin ', Tmin, &
    '  Tmax ', Tmax, '  Pref ', dij_pref
  if (ios /= 0) error stop 'read_idealgas_diffusion failed'
  call flint_acc_upload_thermo()

  open(newunit=u, file=trim(fstates), access='stream', form='unformatted', status='old')
  read(u) n
  allocate(rhoi(NSP,n), T(n), p(n), rho(n), Td(n), Tl(n), Dh(NSP,n), Dd(NSP,n))
  read(u) rhoi, T, p
  close(u)
  write(*,'(A,I0)') ' states ', n
  do i = 1, n
    rho(i) = 0.d0
    do s = 1, NSP
      rho(i) = rho(i) + rhoi(s,i)
    end do
    Tl(i) = idint(T(i)); Td(i) = T(i) - Tl(i)
  end do

  do i = 1, n
    Ti(1) = Tl(i); Ti(2) = Tl(i) + 1
    call co_DS_expr(rhoi(:,i), rho(i), Ti, Td(i), p(i), Dh(:,i))
  end do

  !$acc parallel loop gang vector private(Ti, rl, dl) copyin(rhoi, rho, Tl, Td, p) copyout(Dd)
  do i = 1, n
    Ti(1) = Tl(i); Ti(2) = Tl(i) + 1
    do s = 1, NSP
      rl(s) = rhoi(s,i)
    end do
    call co_DS_expr_dev(rl, rho(i), Ti, Td(i), p(i), dl)
    do s = 1, NSP
      Dd(s,i) = dl(s)
    end do
  end do

  relmax = 0.d0; ndiff = 0; nnan_h = 0; nnan_d = 0
  do i = 1, n
    do s = 1, NSP
      if (Dh(s,i) /= Dh(s,i)) nnan_h = nnan_h + 1
      if (Dd(s,i) /= Dd(s,i)) nnan_d = nnan_d + 1
      if (Dd(s,i) /= Dh(s,i)) then
        ndiff(s) = ndiff(s) + 1
        rel = abs(Dd(s,i) - Dh(s,i))/max(abs(Dh(s,i)), tiny(1.d0))
        relmax(s) = max(relmax(s), rel)
      end if
    end do
  end do
  write(*,'(A,2I10)') ' NaN host / device: ', nnan_h, nnan_d
  do s = 1, NSP
    write(*,'(A,A8,A,I10,A,ES10.3)') ' species ', trim(species_names(s)), '  non-identical ', ndiff(s), '  max rel diff ', relmax(s)
  end do
  write(*,'(A,ES10.3)') ' MAX REL DIFF HOST-DEVICE ', maxval(relmax)
  open(newunit=u, file=trim(fout), access='stream', form='unformatted', status='replace')
  write(u) n
  write(u) Dh, Dd
  close(u)
  write(*,'(A)') ' TEST_CO_DS_DEV_DONE'
end program test_co_DS_dev
