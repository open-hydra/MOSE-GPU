!------------------------------------------------------------------------------
! Single-cell RODAS4 host-vs-device V0 gate (rodas_dev port validation).
!
! Integrates 4 stiff states extracted from the dettube field (real MOSE states:
! cold reactant -> induction -> mid-reaction -> near-equilibrium) as fresh
! [0,dt] Strang segments, four ways:
!   (a) reference  : OSlo wrap_rodas (host RODAS4, IJAC=0 numerical Jacobian)
!   (b) port-host  : FLINT rodas_dev on the host (analytic Jacobian)
!   (c) device     : rodas_dev inside an !$acc parallel loop, NREP replicas
!   (d) cross-meth : radau5_dev on device (independent method, same ODE)
!
! Gates:
!   (a)~(b) within a few*tol  -> the ported RODAS4 stepping is faithful to
!                                Hairer's RODAS (analytic vs numerical Jac only)
!   (b)~(c) within a few*tol  -> device arithmetic acceptable
!   (c) replicas identical    -> deterministic device execution
!   (d)~(c) within a few*tol  -> RODAS and RADAU5 land on the same state
!                                (method-independent physical correctness)
!
! Run from a dir with INPUT/ (dettube's ONERA-7) + single_cell_states.nml.
!------------------------------------------------------------------------------
program single_cell_rodas
  use, intrinsic :: iso_fortran_env, only: output_unit
  use FLINT_Load_ThermoTransport, only: read_idealgas_thermo
  use FLINT_Load_Chemistry,       only: read_chemistry
  use FLINT_Lib_Thermodynamic,    only: ns
  use FLINT_Lib_Chemistry_wdot,   only: Assign_Mechanism
  use FLINT_Lib_Chemistry_rhs,    only: rhs_native, jac_native
  use FLINT_Lib_Radau5_dev,       only: radau5_dev, rodas_dev, flint_acc_upload_tables, NSMX, &
                                        set_radau5_dev_tols, rtolt_dev, atolt_dev, fnewt_dev
  use oslo,                       only: setup_odesolver, run_odesolver
  implicit none

  integer, parameter :: nstates = 4, nseg = 40, nrep = 256
  character(len=16), parameter :: state_name(nstates) = &
    [character(len=16) :: 'early_induction', 'induction', 'mid_reaction', 'near_equilibrium']
  real(8), parameter :: dtseg(nstates) = [1.0d-8, 1.0d-9, 1.0d-9, 1.0d-9]

  real(8) :: roi(7), T, p
  namelist /cell_state_early_induction/  roi, T, p
  namelist /cell_state_induction/        roi, T, p
  namelist /cell_state_mid_reaction/     roi, T, p
  namelist /cell_state_near_equilibrium/ roi, T, p

  character(len=128) :: mech
  integer :: neq, ios, u, is, iseg, c, s, i, idid
  integer :: nf, nj, nstp, na, nr
  real(8) :: z0(NSMX, nstates)
  real(8) :: zref(NSMX, nseg, nstates)   ! (a) OSlo H-rodas reference
  real(8) :: zhost(NSMX, nseg, nstates)  ! (b) rodas_dev host
  real(8), allocatable :: z0col(:,:), dtcol(:), zhist(:,:,:), zrad(:,:,:)
  integer, allocatable :: iddev(:), idrad(:)
  real(8) :: zloc(NSMX), rt(NSMX), at(NSMX)
  real(8), allocatable :: rt_o(:), at_o(:)
  real(8) :: t0, t1, dab, dbc, drep, dcd, floorv
  real(8) :: dmax_ab, dmax_bc, dmax_rep, dmax_cd
  real(8) :: fnw
  integer :: ncol, col0
  logical :: ok_all

  ! ---- FLINT init (thermo tables, then mechanism) --------------------------
  ios = read_idealgas_thermo();  if (ios /= 0) stop 'FATAL: read_idealgas_thermo'
  ios = read_chemistry(mech_name=mech); if (ios /= 0) stop 'FATAL: read_chemistry'
  if (trim(mech) /= 'ONERA-7') stop 'FATAL: expected ONERA-7'
  if (ns /= 7) stop 'FATAL: expected ns=7'
  call Assign_Mechanism(trim(mech))
  neq = ns + 1

  ! ---- ODE setup: OSlo H-rodas (reference solver) --------------------------
  allocate(rt_o(neq), at_o(neq)); rt_o = 1d-5; at_o = 1d-5
  call setup_odesolver(N=neq, solver='H-rodas', RT=rt_o, AT=at_o)
  rt(1:neq) = 1d-5; at(1:neq) = 1d-5
  call flint_acc_upload_tables()
  ! CT build: radau5_dev consumes the PRE-transformed device tolerances; rodas_dev
  ! uses raw rt/at (no transform) and ignores fnw. set_radau5_dev_tols populates
  ! rtolt_dev/atolt_dev (on device) + fnewt_dev.
  call set_radau5_dev_tols(neq, rt(1:neq), at(1:neq))
  fnw = fnewt_dev

  ! ---- read the 4 states ----------------------------------------------------
  open(newunit=u, file='single_cell_states.nml', status='old', action='read', iostat=ios)
  if (ios /= 0) stop 'FATAL: single_cell_states.nml not found'
  read(u, nml=cell_state_early_induction);  z0(1:7,1)=roi; z0(8,1)=T
  rewind(u); read(u, nml=cell_state_induction);        z0(1:7,2)=roi; z0(8,2)=T
  rewind(u); read(u, nml=cell_state_mid_reaction);     z0(1:7,3)=roi; z0(8,3)=T
  rewind(u); read(u, nml=cell_state_near_equilibrium); z0(1:7,4)=roi; z0(8,4)=T
  close(u)

  write(output_unit,'(A)') ' single-cell RODAS4 V0: (a) OSlo wrap_rodas, (b) rodas_dev host, (c) rodas_dev device, (d) radau5_dev device'
  write(output_unit,'(A)') ' build: _OPENACC -> (c),(d) run on the GPU'

  ! ---- (a) reference: OSlo H-rodas -----------------------------------------
  do is = 1, nstates
    zloc(1:neq) = z0(1:neq, is)
    do iseg = 1, nseg
      t0 = 0d0; t1 = dtseg(is)
      call run_odesolver(neq, t0, t1, zloc(1:neq), rhs_native, idid)
      if (idid < 0) write(output_unit,'(A,I3,2A)') ' [a] IDID=', idid, ' state ', trim(state_name(is))
      do s = 1, ns
        if (zloc(s) <= 1d-20) zloc(s) = 1d-20
      end do
      zref(1:neq, iseg, is) = zloc(1:neq)
    end do
  end do

  ! ---- (b) rodas_dev on host ------------------------------------------------
  do is = 1, nstates
    zloc(1:neq) = z0(1:neq, is)
    do iseg = 1, nseg
      call rodas_dev(neq, 0d0, zloc(1:neq), dtseg(is), rt(1:neq), at(1:neq), &
                     fnw, idid, nf, nj, nstp, na, nr)
      if (idid < 0) write(output_unit,'(A,I3,2A)') ' [b] IDID=', idid, ' state ', trim(state_name(is))
      do s = 1, ns
        if (zloc(s) <= 1d-20) zloc(s) = 1d-20
      end do
      zhost(1:neq, iseg, is) = zloc(1:neq)
    end do
  end do

  ! ---- (c) rodas_dev device + (d) radau5_dev device (nrep replicas) --------
  ncol = nstates*nrep
  allocate(z0col(neq, ncol), dtcol(ncol), zhist(neq, nseg, ncol), zrad(neq, nseg, ncol))
  allocate(iddev(ncol), idrad(ncol))
  do is = 1, nstates
    do c = 1, nrep
      z0col(:, (is-1)*nrep+c) = z0(1:neq, is)
      dtcol((is-1)*nrep+c)    = dtseg(is)
    end do
  end do
  iddev = 0; idrad = 0

  ! (c) rodas_dev  (raw rt/at; fnw ignored by RODAS)
  !$acc parallel loop gang vector copyin(z0col, dtcol, rt, at, fnw) &
  !$acc&  copyout(zhist, iddev) private(zloc) &
  !$acc&  private(iseg, s, idid, nf, nj, nstp, na, nr)
  do c = 1, ncol
    zloc(1:8) = z0col(1:8, c)
    do iseg = 1, nseg
      call rodas_dev(8, 0d0, zloc(1:8), dtcol(c), rt(1:8), at(1:8), &
                     fnw, idid, nf, nj, nstp, na, nr)
      if (idid < 0 .and. iddev(c) == 0) iddev(c) = idid
      do s = 1, 7
        if (zloc(s) <= 1d-20) zloc(s) = 1d-20
      end do
      do s = 1, 8
        zhist(s, iseg, c) = zloc(s)
      end do
    end do
  end do

  ! (d) radau5_dev (cross-method; CT -> pre-transformed device tolerances rtolt_dev/atolt_dev)
  !$acc parallel loop gang vector copyin(z0col, dtcol, fnw) &
  !$acc&  copyout(zrad, idrad) private(zloc) &
  !$acc&  present(rtolt_dev, atolt_dev) &
  !$acc&  private(iseg, s, idid, nf, nj, nstp, na, nr)
  do c = 1, ncol
    zloc(1:8) = z0col(1:8, c)
    do iseg = 1, nseg
      call radau5_dev(8, 0d0, zloc(1:8), dtcol(c), rtolt_dev(1:8), atolt_dev(1:8), &
                      fnw, idid, nf, nj, nstp, na, nr)
      if (idid < 0 .and. idrad(c) == 0) idrad(c) = idid
      do s = 1, 7
        if (zloc(s) <= 1d-20) zloc(s) = 1d-20
      end do
      do s = 1, 8
        zrad(s, iseg, c) = zloc(s)
      end do
    end do
  end do

  if (any(iddev < 0)) write(output_unit,'(A,I6,A)') ' [c] WARN: ', count(iddev<0), ' rodas_dev device IDID<0'
  if (any(idrad < 0)) write(output_unit,'(A,I6,A)') ' [d] WARN: ', count(idrad<0), ' radau5_dev device IDID<0'

  ! ---- comparisons ----------------------------------------------------------
  ok_all = .true.
  write(output_unit,'(A)') ''
  write(output_unit,'(A)') ' state              T_end[K](ref)  a-vs-b (port)   b-vs-c (device)  c-replicas   c-vs-d (RODAS/RADAU)'
  do is = 1, nstates
    col0 = (is-1)*nrep + 1
    dmax_ab = 0d0; dmax_bc = 0d0; dmax_rep = 0d0; dmax_cd = 0d0
    do iseg = 1, nseg
      do i = 1, neq
        floorv = merge(1d0, 1d-12, i == neq)   ! T floor 1, species floor 1e-12
        dab = abs(zhost(i,iseg,is) - zref(i,iseg,is)) / max(abs(zref(i,iseg,is)), floorv)
        dbc = abs(zhist(i,iseg,col0) - zhost(i,iseg,is)) / max(abs(zhost(i,iseg,is)), floorv)
        dcd = abs(zhist(i,iseg,col0) - zrad(i,iseg,col0)) / max(abs(zrad(i,iseg,col0)), floorv)
        dmax_ab = max(dmax_ab, dab)
        dmax_bc = max(dmax_bc, dbc)
        dmax_cd = max(dmax_cd, dcd)
        do c = col0+1, col0+nrep-1
          drep = abs(zhist(i,iseg,c) - zhist(i,iseg,col0))
          dmax_rep = max(dmax_rep, drep)
        end do
      end do
    end do
    write(output_unit,'(1X,A18,F9.1,5X,4(ES13.3,3X))') state_name(is), zref(neq,nseg,is), &
          dmax_ab, dmax_bc, dmax_rep, dmax_cd
    if (dmax_ab  > 1d-2 ) ok_all = .false.   ! ported RODAS matches Hairer to tol
    if (dmax_bc  > 1d-3 ) ok_all = .false.   ! device within a few * tol
    if (dmax_rep /= 0d0 ) ok_all = .false.   ! replicas must agree exactly
    if (dmax_cd  > 5d-2 ) ok_all = .false.   ! RODAS vs RADAU5 same equilibrium
  end do

  write(output_unit,'(A)') ''
  if (ok_all) then
    write(output_unit,'(A)') ' RESULT: PASS'
  else
    write(output_unit,'(A)') ' RESULT: FAIL'
    error stop 1
  end if

end program single_cell_rodas
