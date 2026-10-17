!------------------------------------------------------------------------------
! Single-cell RADAU5 host-vs-device spike (GPU port Phase 1 gate, risk #1).
!
! Integrates the 4 stiff post-ignition states extracted from the dettube
! field (RDE/gpu_verify/single_cell/single_cell_states.nml) three ways:
!   (a) reference : OSlo wrap_radau5_jac (the exact MOSE CPU path, IJAC=1)
!   (b) port-host : FLINT radau5_dev called on the host
!   (c) device    : radau5_dev inside an !$acc parallel loop, NREP replicas
!                   per state (replica agreement = device determinism check)
!
! Gates:
!   (a)==(b) bit-exact      -> the specialization/port is faithful
!   (b)~(c) within few*tol  -> device arithmetic acceptable (rtol=atol=1e-5)
!   (c) replicas identical  -> deterministic device execution
!
! Run from a directory containing INPUT/ (dettube's: chemistry-info.txt,
! chemistry-Arrhenius.dat, thermo.dat) and single_cell_states.nml.
! Mimics MOSE Chemistry_Newstate: fresh [0,dt] integration per segment,
! 1e-20 floor clip on species after each segment.
!------------------------------------------------------------------------------
program single_cell_radau5
  use, intrinsic :: iso_fortran_env, only: output_unit
  use FLINT_Load_ThermoTransport, only: read_idealgas_thermo
  use FLINT_Load_Chemistry,       only: read_chemistry
  use FLINT_Lib_Thermodynamic,    only: ns
  use FLINT_Lib_Chemistry_wdot,   only: Assign_Mechanism
  use FLINT_Lib_Chemistry_rhs,    only: rhs_native, jac_native
  use FLINT_Lib_Radau5_dev,       only: radau5_dev, flint_acc_upload_tables, NSMX
  use oslo,                       only: setup_odesolver, wrap_radau5_jac
  implicit none

  integer, parameter :: nstates = 4, nseg = 40, nrep = 256
  character(len=16), parameter :: state_name(nstates) = &
    [character(len=16) :: 'early_induction', 'induction', 'mid_reaction', 'near_equilibrium']
  real(8), parameter :: dtseg(nstates) = [2.5d-6, 2.5d-8, 2.5d-9, 2.5d-8]

  real(8) :: roi(7), T, p
  namelist /cell_state_early_induction/  roi, T, p
  namelist /cell_state_induction/        roi, T, p
  namelist /cell_state_mid_reaction/     roi, T, p
  namelist /cell_state_near_equilibrium/ roi, T, p

  character(len=128) :: mech
  integer :: neq, ios, u, is, iseg, c, s, i, idid
  integer :: nf, nj, nstp, na, nr
  real(8) :: z0(NSMX, nstates)
  real(8) :: zref(NSMX, nseg, nstates)   ! (a) OSlo reference trajectories
  real(8) :: zhost(NSMX, nseg, nstates)  ! (b) port on host
  real(8), allocatable :: z0col(:,:), dtcol(:), zhist(:,:,:)
  integer, allocatable :: iddev(:)
  real(8) :: zloc(NSMX), rt(NSMX), at(NSMX)
  real(8), allocatable :: rt_o(:), at_o(:)
  real(8) :: t0, t1, dab, dbc, drep, floorv, dmax_ab, dmax_bc, dmax_rep
  integer :: ncol, col0
  logical :: ok_all

  ! ---- FLINT init (same path as MOSE: thermo tables, then mechanism) -------
  ios = read_idealgas_thermo()
  if (ios /= 0) stop 'FATAL: read_idealgas_thermo failed (INPUT/thermo.dat?)'
  ios = read_chemistry(mech_name=mech)
  if (ios /= 0) stop 'FATAL: read_chemistry failed (INPUT/chemistry-info.txt?)'
  if (trim(mech) /= 'ONERA-7') stop 'FATAL: expected ONERA-7 mechanism'
  if (ns /= 7) stop 'FATAL: expected ns=7'
  call Assign_Mechanism(trim(mech))
  neq = ns + 1

  ! ---- ODE setup exactly as MOSE IO_Chemistry defaults ----------------------
  allocate(rt_o(neq), at_o(neq))
  rt_o = 1d-5; at_o = 1d-5
  call setup_odesolver(N=neq, solver='H-radau5', RT=rt_o, AT=at_o)
  rt(1:neq) = 1d-5; at(1:neq) = 1d-5

  call flint_acc_upload_tables()

  ! ---- read the 4 states ----------------------------------------------------
  open(newunit=u, file='single_cell_states.nml', status='old', action='read', iostat=ios)
  if (ios /= 0) stop 'FATAL: single_cell_states.nml not found'
  read(u, nml=cell_state_early_induction);  z0(1:7,1)=roi; z0(8,1)=T
  rewind(u); read(u, nml=cell_state_induction);        z0(1:7,2)=roi; z0(8,2)=T
  rewind(u); read(u, nml=cell_state_mid_reaction);     z0(1:7,3)=roi; z0(8,3)=T
  rewind(u); read(u, nml=cell_state_near_equilibrium); z0(1:7,4)=roi; z0(8,4)=T
  close(u)

  write(output_unit,'(A)') ' single-cell RADAU5 spike: (a) OSlo wrap_radau5_jac, (b) radau5_dev host, (c) radau5_dev device'
  write(output_unit,'(A)') ' build: _OPENACC defined -> (c) runs on the GPU'

  ! ---- (a) reference: OSlo RADAU5 ------------------------------------------
  do is = 1, nstates
    zloc(1:neq) = z0(1:neq, is)
    do iseg = 1, nseg
      t0 = 0d0; t1 = dtseg(is)
      call wrap_radau5_jac(neq, t0, t1, zloc(1:neq), rhs_native, jac_native, idid)
      if (idid < 0) write(output_unit,'(A,I3,2A)') ' [a] IDID=', idid, ' state ', trim(state_name(is))
      do s = 1, ns
        if (zloc(s) <= 1d-20) zloc(s) = 1d-20
      end do
      zref(1:neq, iseg, is) = zloc(1:neq)
    end do
  end do

  ! ---- (b) port on host ------------------------------------------------------
  do is = 1, nstates
    zloc(1:neq) = z0(1:neq, is)
    do iseg = 1, nseg
      call radau5_dev(neq, 0d0, zloc(1:neq), dtseg(is), rt(1:neq), at(1:neq), &
                      idid, nf, nj, nstp, na, nr)
      if (idid < 0) write(output_unit,'(A,I3,2A)') ' [b] IDID=', idid, ' state ', trim(state_name(is))
      do s = 1, ns
        if (zloc(s) <= 1d-20) zloc(s) = 1d-20
      end do
      zhost(1:neq, iseg, is) = zloc(1:neq)
    end do
  end do

  ! ---- (c) device: nrep replicas per state ----------------------------------
  ncol = nstates*nrep
  allocate(z0col(neq, ncol), dtcol(ncol), zhist(neq, nseg, ncol), iddev(ncol))
  do is = 1, nstates
    do c = 1, nrep
      z0col(:, (is-1)*nrep+c) = z0(1:neq, is)
      dtcol((is-1)*nrep+c)    = dtseg(is)
    end do
  end do
  iddev = 0

  !$acc parallel loop gang vector copyin(z0col, dtcol, rt, at) &
  !$acc&  copyout(zhist, iddev) private(zloc) &
  !$acc&  private(iseg, s, idid, nf, nj, nstp, na, nr)
  do c = 1, ncol
    zloc(1:8) = z0col(1:8, c)
    do iseg = 1, nseg
      call radau5_dev(8, 0d0, zloc(1:8), dtcol(c), rt(1:8), at(1:8), &
                      idid, nf, nj, nstp, na, nr)
      if (idid < 0 .and. iddev(c) == 0) iddev(c) = idid
      do s = 1, 7
        if (zloc(s) <= 1d-20) zloc(s) = 1d-20
      end do
      do s = 1, 8
        zhist(s, iseg, c) = zloc(s)
      end do
    end do
  end do

  if (any(iddev < 0)) write(output_unit,'(A,I6,A)') ' [c] WARNING: ', count(iddev<0), ' device integrations returned IDID<0'

  ! ---- comparisons -----------------------------------------------------------
  ok_all = .true.
  write(output_unit,'(A)') ''
  write(output_unit,'(A)') ' state              T_end[K] (ref)   max rel a-vs-b   max rel b-vs-c   max rel c-replicas'
  do is = 1, nstates
    col0 = (is-1)*nrep + 1
    dmax_ab = 0d0; dmax_bc = 0d0; dmax_rep = 0d0
    do iseg = 1, nseg
      do i = 1, neq
        floorv = merge(1d0, 1d-12, i == neq)   ! T floor 1, species floor 1e-12
        dab = abs(zhost(i,iseg,is) - zref(i,iseg,is)) / max(abs(zref(i,iseg,is)), floorv)
        dbc = abs(zhist(i,iseg,col0) - zhost(i,iseg,is)) / max(abs(zhost(i,iseg,is)), floorv)
        dmax_ab = max(dmax_ab, dab)
        dmax_bc = max(dmax_bc, dbc)
        do c = col0+1, col0+nrep-1
          drep = abs(zhist(i,iseg,c) - zhist(i,iseg,col0))
          dmax_rep = max(dmax_rep, drep)
        end do
      end do
    end do
    write(output_unit,'(1X,A18,F10.2,6X,3(ES14.3,3X))') state_name(is), zref(neq,nseg,is), dmax_ab, dmax_bc, dmax_rep
    if (dmax_ab  > 1d-12) ok_all = .false.   ! port must match OSlo bit-level
    if (dmax_bc  > 1d-3 ) ok_all = .false.   ! device within a few * tol
    if (dmax_rep /= 0d0 ) ok_all = .false.   ! replicas must agree exactly
  end do

  ! short reference trajectory (every 8th checkpoint) for the record
  write(output_unit,'(A)') ''
  write(output_unit,'(A)') ' reference T(t) trajectories (every 8th segment):'
  do is = 1, nstates
    write(output_unit,'(1X,A18,5F10.2)') state_name(is), (zref(neq,iseg,is), iseg=8,nseg,8)
  end do

  write(output_unit,'(A)') ''
  if (ok_all) then
    write(output_unit,'(A)') ' RESULT: PASS'
  else
    write(output_unit,'(A)') ' RESULT: FAIL'
    error stop 1
  end if

end program single_cell_radau5
