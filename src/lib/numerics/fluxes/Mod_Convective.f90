!> GPU CUDA-Fortran prototype — CONVECTIVE flux pass-1 as a CUF attributes(global)
!> kernel, to break the OpenACC occupancy/launch_bounds ceiling. Pass-2 stays OpenACC.
!>
!> WHY: the OpenACC convective kernel (Mod_Convective_W.f90) is register-bound; OpenACC
!> has no per-kernel launch_bounds and no user-managed shared memory. This module is a
!> SELF-CONTAINED CUF copy: the physics chain (Convective_Flux -> Reconstruction/MINMOD +
!> co_rotot_Rtot/f_ss + riemann_HLLC/roe_averages) is ported VERBATIM to attributes(device)
!> routines. Module SCALARS that the chain reads (indices) become value args. The FLINT
!> thermo TABLES are held as attributes(device) MODULE COPIES (Ri/wm/cp/h _d), populated
!> ONCE from the host tables in cuf_init_tables, and read DIRECTLY by the device routines
!> (gattini's pattern). This avoids NVFORTRAN-S-0520 (device routines cannot read
!> declare-create module data) AND the multi-level host_data table-pointer loss that the
!> first prototype hit. Resident w/fhat/fm/dl/beta are handed over via host_data use_device.
!> The ALGEBRA is byte-identical to the OpenACC path -> gate 0.0. M1 = no shared memory.
!>
!> GATED behind _OPENACC (the var-last GPU flux path). This is now the
!> UNCONDITIONAL production convective kernel; the old OpenACC scatter fallback is gone.
module MOSE_Mod_Convective
  use cudafor
  use iso_fortran_env, only: R8 => real64
  use MOSE_Global_m,   only: NPRIM_DEV
  use FLINT_Lib_Chemistry_data, only: NSMAX

  use MOSE_Global_m,   only: gc, NPRIM_DEV
  use MOSE_Field_m,    only: P_flat, field_block_index_params
  implicit none
  private
  public :: Convective_Flux

  real(R8), parameter :: Runiv = 8314.46261815324d0   ! mirror FLINT_Lib_Thermodynamic Runiv (exact SI value)

  ! Device copies of the FLINT thermo tables, read DIRECTLY by the device routines.
  ! Populated once from the host tables (cuf_init_tables). Bounds mirror the host arrays
  ! so the (Tmin:Tmax) temperature index and (ns) species index are correct.
  integer, parameter :: NTMAX = 10000   ! >= table Tmax (ONERA-7 Tmax=10000)
  include 'mose_ct_dims.inc'
  ! One binary per species count: the tables follow the compile-time nsc, FLINT's NSMAX (8) is only the floor
  ! (MOSE_open's host tables are allocatable (Tmin:Tmax, 1:ns) with no cap). A universal build keeps NSMAX
  ! and stops below when the mechanism is bigger.
  integer, parameter :: NS_TAB = max(NSMAX, nsc)
  real(R8), device, save :: Ri_tab_d(NS_TAB), wm_tab_d(NS_TAB)
  real(R8), device, save :: cp_tab_d(0:NTMAX, NS_TAB), h_tab_d(0:NTMAX, NS_TAB)
  logical,  save :: cuf_tables_ready = .false.
  ! The MOSE_open solvers 4..14 (riemann_open_d): FLINT's table bounds for the clamp and riemann-options-Mco.
  ! (device ARRAYS as the tables: module device SCALARS read by the CUF kernels clash in the NVVM module, pgf90_dev_common)
  integer,  device, save :: ttab_d(2)    ! FLINT Tmin, Tmax
  real(R8), device, save :: mco_d(1)     ! riemann-options-Mco

  public :: Convective_assert_scheme
  public :: convective_alloc_fhat
  public :: fhat, fhat_j, fhat_k   ! exposed for the CUF launcher + the fused pass-2 (Mod_Fluxes)
  real(R8), allocatable, save :: fhat(:,:), fhat_j(:,:), fhat_k(:,:)
  integer,  save :: nface_max = 0
  logical,  save :: fhat_resident = .false.

contains

  !> One-time host->device copy of the 4 thermo tables the convective chain needs.
  subroutine cuf_init_tables ()
    use FLINT_Lib_Thermodynamic, only: ns, Tmin, Tmax, Ri_tab, wm_tab, cp_tab, h_tab
    use MOSE_Lib_Riemann, only: mco_dev
    if (cuf_tables_ready) return
    if (Tmax > NTMAX) error stop 'MOSE_CONV_W_CUF: table Tmax exceeds NTMAX'
    if (ns > NS_TAB) then
      write(*,*) '  [ERROR] ', ns, ' species exceed the device thermo table size NS_TAB = ', NS_TAB, &
                 ' (FLINT NSMAX): the binary is built per species count (mose does it)'
      error stop 'MOSE_CONV_W_CUF: species count exceeds NS_TAB'
    end if
    !$acc update host(Ri_tab, wm_tab, cp_tab, h_tab)
    Ri_tab_d(1:ns) = Ri_tab(1:ns)
    wm_tab_d(1:ns) = wm_tab(1:ns)
    cp_tab_d(Tmin:Tmax, 1:ns) = cp_tab(Tmin:Tmax, 1:ns)
    h_tab_d(Tmin:Tmax, 1:ns)  = h_tab(Tmin:Tmax, 1:ns)
    ttab_d(1:2) = [ Tmin, Tmax ]; mco_d(1:1) = [ mco_dev ]
    cuf_tables_ready = .true.
  end subroutine cuf_init_tables

  !> MOSE_open AUSM+ (4), AUSM+M (5), HLLE++ (6), HLLC+Tramel (7) for the CUF convective kernels: the shared body of
  !> riemann_open_body.inc (the OpenACC twin is Lib_Riemann riemann_open), on the device tables above.
  attributes(device) subroutine riemann_open_d ( rs, bet, mco, runiv_, ns_, &
                                                 dl, ul, vl, wl, pl, al, dltot, dr, ur, vr, wr, pr, ar, drtot, nx, ny, nz, &
                                                 urL_, urR_, F_r, F_u, F_v, F_w, F_E )
#define RO_HTAB(i,s) h_tab_d(i,s)
#define RO_CPTAB(i,s) cp_tab_d(i,s)
#define RO_RITAB(s) Ri_tab_d(s)
#define RO_WMTAB(s) wm_tab_d(s)
#define RO_TLO ttab_d(1)
#define RO_THI ttab_d(2)
#define RO_MCO mco_d(1)
#define RO_INV value
#include "riemann_open_body.inc"
#undef RO_HTAB
#undef RO_CPTAB
#undef RO_RITAB
#undef RO_WMTAB
#undef RO_TLO
#undef RO_THI
#undef RO_MCO
#undef RO_INV
  end subroutine riemann_open_d

  ! ============================================================================
  ! PORTED DEVICE PHYSICS (verbatim bodies; tables read from the device module copies).
  ! ============================================================================

  attributes(device) pure function f_Rtot_d(rhoi, ns) result(result)
    !$acc routine seq
    integer, value :: ns
    real(R8), intent(in) :: rhoi(ns)
    real(R8) :: result, rho
    integer :: s
    rho = sum(rhoi)
    result = 0.d0
    do s = 1, ns
      result = result + Ri_tab_d(s)*rhoi(s)/rho
    enddo
  end function f_Rtot_d

  attributes(device) pure function f_cp_d(rhoi, Tint, Tdiff, rho, ns) result(res)
    !$acc routine seq
    integer, value :: ns
    real(R8), intent(in) :: rhoi(ns), Tdiff, rho
    integer, intent(in)  :: Tint(2)
    real(R8) :: res, cpi
    integer :: s
    res = 0.d0
    do s = 1, ns
      cpi = cp_tab_d(Tint(1),s) + (cp_tab_d(Tint(2),s) - cp_tab_d(Tint(1),s))*Tdiff
      res = res + rhoi(s)/rho*cpi
    enddo
  end function f_cp_d

  attributes(device) pure function f_ss_d(rhoi, p, rho, Rtot, ns) result(res)
    !$acc routine seq
    integer, value :: ns
    real(R8), intent(in) :: rhoi(ns), p, rho, Rtot
    real(R8) :: res, T, cp, gam, Tdiff
    integer :: T_i, Tint(2)
    T = p/(Rtot*rho)
    T_i = idint(T)
    if (T_i < 0) T_i = 0
    if (T_i > NTMAX-1) T_i = NTMAX-1
    Tdiff = T - T_i
    Tint(1) = T_i; Tint(2) = T_i + 1
    cp  = f_cp_d(rhoi, Tint, Tdiff, rho, ns)
    gam = cp/(cp - Rtot)
    res = dsqrt(gam*Rtot*T)
  end function f_ss_d

  attributes(device) pure subroutine co_rotot_Rtot_d(rhoi, rho, R, ns)
    !$acc routine seq
    integer, value :: ns
    real(R8), intent(in)  :: rhoi(ns)
    real(R8), intent(out) :: rho, R
    rho = sum(rhoi(1:ns))
    R   = f_Rtot_d(rhoi, ns)
  end subroutine co_rotot_Rtot_d

  !> The device limiter dispatched by the kernel's limv argument (= lim_dev): 1 minmod (verbatim rlimiter_MINMOD_d),
  !> 2 vanleer, 3 vanalbada, 4 mc, 5 superbee (Lib_Limiters rlimiter).
  attributes(device) pure function rlimiter_d(lim, x, y) result(result)
    !$acc routine seq
    integer,  intent(in) :: lim
    real(R8), intent(in) :: x, y
    real(R8) :: ru, phi, result
    select case (lim)
    case (1)   ! Min-Mod: the statements of rlimiter_MINMOD (bit-identical)
      result = 0.d0
      if (x*y<0.d0) return
      result = x
      if (abs(x)>abs(y)) result = y
    case (2)   ! Van Leer (rlimiter_VANLEER)
      ru = y/(x + 1d-30)
      phi = (ru+abs(ru))/(1d0+abs(ru))
      result = phi*x
    case (3)   ! Van Albada (rlimiter_VANALBADA)
      ru = y/(x+1d-30)
      phi = (ru*ru+ru)/(1d0+ru*ru)
      result = phi*x
    case (4)   ! MC (rlimiter_MC: maxval/minval of the same values)
      ru = y/(x+1d-30)
      phi = max(0d0, min(2d0*ru, 0.5d0*(1d0+ru), 2d0))
      result = phi*x
    case (5)   ! Super-Bee (rlimiter_SB)
      ru = y/(x+1d-30)
      phi = max(0d0, min(2d0*ru,1d0), min(ru,2d0))
      result = phi*x
    case default
      result = 0d0
    end select
  end function rlimiter_d

  attributes(device) pure function rlimiter_MINMOD_d(x, y) result(result)
    !$acc routine seq
    real(R8), intent(in) :: x, y
    real(R8) :: result
    result = 0.d0
    if (x*y < 0.d0) return
    result = x
    if (abs(x) > abs(y)) result = y
  end function rlimiter_MINMOD_d

  attributes(device) pure subroutine roe_averages_d(nsc, dl, dr, dltot, drtot, ul, vl, wl, ur, vr, wr, &
       pl, pr, h0l, h0r, rho_roe, u_roe, v_roe, w_roe, a_roe, h0_roe)
    !$acc routine seq
    integer, value :: nsc
    real(R8), intent(in) :: dl(nsc), dr(nsc), dltot, drtot
    real(R8), intent(in) :: ul, vl, wl, ur, vr, wr, pl, pr
    real(R8), intent(out) :: h0l, h0r, rho_roe, u_roe, v_roe, w_roe, a_roe, h0_roe
    integer :: s, Til, Tir
    real(R8) :: Rl, Rr, Tl, Tr, dTl, dTr, inv_dltot, inv_drtot, srL, srR, inv_sr
    real(R8) :: cv_roe, sum_ei, R_roe, T_roe, gam_roe, vel2, hl, hr, el, er, invW, d_roe
    Rl = f_Rtot_d(dl, nsc)
    Rr = f_Rtot_d(dr, nsc)
    inv_dltot = 1.d0 / dltot
    inv_drtot = 1.d0 / drtot
    Tl = pl * inv_dltot / Rl
    Tr = pr * inv_drtot / Rr
    Til = int(Tl); Tir = int(Tr)
    if (Til < 0) Til = 0
    if (Til > NTMAX-1) Til = NTMAX-1
    if (Tir < 0) Tir = 0
    if (Tir > NTMAX-1) Tir = NTMAX-1
    dTl = Tl - Til; dTr = Tr - Tir
    srL = sqrt(dltot); srR = sqrt(drtot)
    inv_sr = 1.d0 / (srL + srR)
    rho_roe = srL * srR
    u_roe = (srR*ur + srL*ul) * inv_sr
    v_roe = (srR*vr + srL*vl) * inv_sr
    w_roe = (srR*wr + srL*wl) * inv_sr
    R_roe = (srR*Rr + srL*Rl) * inv_sr
    T_roe = (srR*Tr + srL*Tl) * inv_sr
    h0l = 0.d0; h0r = 0.d0; cv_roe = 0.d0; sum_ei = 0.d0
    do s = 1, nsc
      invW = Runiv / wm_tab_d(s)
      hl = h_tab_d(Til, s) + (h_tab_d(Til+1, s) - h_tab_d(Til, s)) * dTl
      hr = h_tab_d(Tir, s) + (h_tab_d(Tir+1, s) - h_tab_d(Tir, s)) * dTr
      el = hl - invW * Tl
      er = hr - invW * Tr
      h0l = h0l + hl * dl(s) * inv_dltot
      h0r = h0r + hr * dr(s) * inv_drtot
      d_roe = (srR*dr(s)*inv_drtot + srL*dl(s)*inv_dltot) * inv_sr
      cv_roe = cv_roe + d_roe * ( &
          0.5d0 * ( &
            cp_tab_d(Til, s) + (cp_tab_d(Til+1, s)-cp_tab_d(Til, s))*dTl + &
            cp_tab_d(Tir, s) + (cp_tab_d(Tir+1, s)-cp_tab_d(Tir, s))*dTr ) &
          - invW )
      sum_ei = sum_ei + d_roe * (srR*er + srL*el) * inv_sr
    end do
    h0l = h0l + 0.5d0*(ul*ul + vl*vl + wl*wl)
    h0r = h0r + 0.5d0*(ur*ur + vr*vr + wr*wr)
    h0_roe = (srR*h0r + srL*h0l) * inv_sr
    gam_roe = 1.d0 + R_roe / cv_roe
    vel2 = u_roe*u_roe + v_roe*v_roe + w_roe*w_roe
    a_roe = sqrt( (gam_roe - 1.d0) * ( h0_roe - 0.5d0*vel2 + cv_roe*T_roe - sum_ei ) )
  end subroutine roe_averages_d

  attributes(device) pure subroutine Reconstruction_d(prev, local, next, next2, dl0, dl1, dl2, dll, dlr, &
       beta, priml, primr, nprim, nsc)
    !$acc routine seq
    integer, value :: nprim, nsc
    real(R8), intent(in)    :: prev(nprim), local(nprim), next(nprim), next2(nprim)
    real(R8), intent(inout) :: priml(nprim), primr(nprim)
    real(R8), intent(in)    :: dl0, dl1, dl2, dll, dlr, beta
    real(R8) :: slope0, slope1, slope2, slopel, sloper, idl0, idl1, idl2
    integer :: i
    idl0 = 1d0/dl0; idl1 = 1d0/dl1; idl2 = 1d0/dl2
    do i = 1, nprim
      slope0 = ( local(i) - prev(i) ) * idl0
      slope1 = ( next(i)  - local(i) ) * idl1
      slope2 = ( next2(i) - next(i) ) * idl2
      slopel = rlimiter_MINMOD_d ( slope1, slope0 )
      sloper = rlimiter_MINMOD_d ( slope2, slope1 )
      priml(i) = local(i) + beta * slopel * dll
      primr(i) =  next(i) - beta * sloper * dlr
    end do
    do i = 1, nsc
      priml(i) = max(priml(i), 1d-20)
      primr(i) = max(primr(i), 1d-20)
    enddo
  end subroutine Reconstruction_d


  ! ============================================================================
  ! CUF PASS-1 KERNELS (one per direction). 2D (i,j) block + in-thread do k.
  ! ============================================================================

  attributes(global) launch_bounds(128) subroutine conv_cuf_i ( &
       w, fhat, fm1, dl, beta, ncell, nfacemax, boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3, &
       SD, limv, rsv, urf, nur)
    integer, value :: ncell, nfacemax, boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3
    logical, value :: SD
    integer, value :: limv   ! the limiter id (lim_dev)
    integer, value :: rsv    ! the Riemann solver id (rs_dev): 1 HLLC, 2 HLLE, 3 SLAU, 4..14 riemann_open_d
    integer, value :: nur    ! the size of urf: ncell with prec, 1 otherwise
    real(R8) :: urf(nur)     ! Ur of the cells (Lib_Preconditioning ur_f), read by HLLC-PC (14) only
    real(R8) :: w(ncell, nprim), fhat(nfacemax, nprim)
    real(R8) :: fm1(4, 0:n1, n2, n3), dl(3, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc), beta(n1, n2, n3)
    integer :: i, j, k, s, f, c0, cm1, cp1, cp2
    real(R8) :: dlc(-1:2), nrm(3), pm1(NPRIM_DEV), p0(NPRIM_DEV), pp1(NPRIM_DEV), pp2(NPRIM_DEV), res0(NPRIM_DEV)
    real(R8) :: area, bet
    integer :: t_, ni_, nj_, rem_
    include 'conv_decls.inc'
    ! flat thread mapping: 1D grid over the (o..h) face range with an
    ! in-kernel decode, i fastest -> coalescing preserved. Kills the per-row partial-warp
    ! waste of the (32,2,2) tile on narrow blocks (THOR: 18-icell injector blocks paid
    ! 32 lanes/row = 43.8% masked; cell-weighted 19.8% across the mesh). Thread mapping
    ! does not touch FP; each face is still computed exactly once -> bit-exact.
    t_  = (blockIdx%x-1)*blockDim%x + threadIdx%x
    ni_ = h1 - o1 + 1
    nj_ = h2 - o2 + 1
    if (t_ > ni_*nj_*(h3-o3+1)) return
    i    = mod(t_-1, ni_) + o1
    rem_ = (t_-1) / ni_
    j    = mod(rem_, nj_) + o2
    k    = rem_ / nj_ + o3
      c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      cm1 = c0 - 1; cp1 = c0 + 1; cp2 = c0 + 2
      dlc(-1)=dl(1,i-1,j,k); dlc(0)=dl(1,i,j,k); dlc(1)=dl(1,i+1,j,k); dlc(2)=dl(1,i+2,j,k)
      nrm(1)=fm1(1,i,j,k); nrm(2)=fm1(2,i,j,k); nrm(3)=fm1(3,i,j,k)
      do s = 1, nprim
        pm1(s)=w(cm1,s); p0(s)=w(c0,s); pp1(s)=w(cp1,s); pp2(s)=w(cp2,s)
      end do
      area = fm1(4,i,j,k); bet = beta(i,j,k)
      include 'conv_physics.inc'
      f = (i-1) + (j-1)*(n1-1) + (k-1)*(n1-1)*n2 + 1
      do s = 1, nprim
        fhat(f, s) = res0(s)
      end do
  end subroutine conv_cuf_i

  attributes(global) launch_bounds(128) subroutine conv_cuf_j ( &
       w, fhat, fm2, dl, beta, ncell, nfacemax, boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3, &
       SD, limv, rsv, urf, nur)
    integer, value :: ncell, nfacemax, boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3
    logical, value :: SD
    integer, value :: limv   ! the limiter id (lim_dev)
    integer, value :: rsv    ! the Riemann solver id (rs_dev): 1 HLLC, 2 HLLE, 3 SLAU, 4..14 riemann_open_d
    integer, value :: nur    ! the size of urf: ncell with prec, 1 otherwise
    real(R8) :: urf(nur)     ! Ur of the cells (Lib_Preconditioning ur_f), read by HLLC-PC (14) only
    real(R8) :: w(ncell, nprim), fhat(nfacemax, nprim)
    real(R8) :: fm2(4, n1, 0:n2, n3), dl(3, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc), beta(n1, n2, n3)
    integer :: i, j, k, s, f, c0, cm1, cp1, cp2
    real(R8) :: dlc(-1:2), nrm(3), pm1(NPRIM_DEV), p0(NPRIM_DEV), pp1(NPRIM_DEV), pp2(NPRIM_DEV), res0(NPRIM_DEV)
    real(R8) :: area, bet
    integer :: t_, ni_, nj_, rem_
    include 'conv_decls.inc'
    ! flat 1D mapping: in-kernel decode, i fastest -> coalescing preserved;
    ! kills the (32,2,2) tile's partial-warp waste on narrow blocks. Bit-exact.
    t_  = (blockIdx%x-1)*blockDim%x + threadIdx%x
    ni_ = h1 - o1 + 1
    nj_ = h2 - o2 + 1
    if (t_ > ni_*nj_*(h3-o3+1)) return
    i    = mod(t_-1, ni_) + o1
    rem_ = (t_-1) / ni_
    j    = mod(rem_, nj_) + o2
    k    = rem_ / nj_ + o3
      c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      cm1 = c0 - pi; cp1 = c0 + pi; cp2 = c0 + 2*pi
      dlc(-1)=dl(2,i,j-1,k); dlc(0)=dl(2,i,j,k); dlc(1)=dl(2,i,j+1,k); dlc(2)=dl(2,i,j+2,k)
      nrm(1)=fm2(1,i,j,k); nrm(2)=fm2(2,i,j,k); nrm(3)=fm2(3,i,j,k)
      do s = 1, nprim
        pm1(s)=w(cm1,s); p0(s)=w(c0,s); pp1(s)=w(cp1,s); pp2(s)=w(cp2,s)
      end do
      area = fm2(4,i,j,k); bet = beta(i,j,k)
      include 'conv_physics.inc'
      f = (i-1) + (j-1)*n1 + (k-1)*n1*(n2-1) + 1
      do s = 1, nprim
        fhat(f, s) = res0(s)
      end do
  end subroutine conv_cuf_j

  attributes(global) launch_bounds(128) subroutine conv_cuf_k ( &
       w, fhat, fm3, dl, beta, ncell, nfacemax, boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3, &
       SD, limv, rsv, urf, nur)
    integer, value :: ncell, nfacemax, boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3
    logical, value :: SD
    integer, value :: limv   ! the limiter id (lim_dev)
    integer, value :: rsv    ! the Riemann solver id (rs_dev): 1 HLLC, 2 HLLE, 3 SLAU, 4..14 riemann_open_d
    integer, value :: nur    ! the size of urf: ncell with prec, 1 otherwise
    real(R8) :: urf(nur)     ! Ur of the cells (Lib_Preconditioning ur_f), read by HLLC-PC (14) only
    real(R8) :: w(ncell, nprim), fhat(nfacemax, nprim)
    real(R8) :: fm3(4, n1, n2, 0:n3), dl(3, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc), beta(n1, n2, n3)
    integer :: i, j, k, s, f, c0, cm1, cp1, cp2
    real(R8) :: dlc(-1:2), nrm(3), pm1(NPRIM_DEV), p0(NPRIM_DEV), pp1(NPRIM_DEV), pp2(NPRIM_DEV), res0(NPRIM_DEV)
    real(R8) :: area, bet
    integer :: t_, ni_, nj_, rem_
    include 'conv_decls.inc'
    ! flat 1D mapping: in-kernel decode, i fastest -> coalescing preserved;
    ! kills the (32,2,2) tile's partial-warp waste on narrow blocks. Bit-exact.
    t_  = (blockIdx%x-1)*blockDim%x + threadIdx%x
    ni_ = h1 - o1 + 1
    nj_ = h2 - o2 + 1
    if (t_ > ni_*nj_*(h3-o3+1)) return
    i    = mod(t_-1, ni_) + o1
    rem_ = (t_-1) / ni_
    j    = mod(rem_, nj_) + o2
    k    = rem_ / nj_ + o3
      c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      cm1 = c0 - pi*pj; cp1 = c0 + pi*pj; cp2 = c0 + 2*pi*pj
      dlc(-1)=dl(3,i,j,k-1); dlc(0)=dl(3,i,j,k); dlc(1)=dl(3,i,j,k+1); dlc(2)=dl(3,i,j,k+2)
      nrm(1)=fm3(1,i,j,k); nrm(2)=fm3(2,i,j,k); nrm(3)=fm3(3,i,j,k)
      do s = 1, nprim
        pm1(s)=w(cm1,s); p0(s)=w(c0,s); pp1(s)=w(cp1,s); pp2(s)=w(cp2,s)
      end do
      area = fm3(4,i,j,k); bet = beta(i,j,k)
      include 'conv_physics.inc'
      f = (i-1) + (j-1)*n1 + (k-1)*n1*n2 + 1
      do s = 1, nprim
        fhat(f, s) = res0(s)
      end do
  end subroutine conv_cuf_k

  ! ============================================================================
  ! HOST LAUNCHER — CUF pass-1 (host_data use_device) + OpenACC pass-2 (verbatim).
  ! Replaces Convective_Flux_w_blk when MOSE_CONV_W_CUF is ON.
  ! ============================================================================
  subroutine Convective_Flux ( b, Res, fm1, fm2, fm3, dl, beta, n)
    use MOSE_Global_m,   only: gc, nprim, nsc, nu, nv, nw, np
    use MOSE_Field_m,    only: P_flat, R_flat, field_block_index_params, ncell_total
    use MOSE_Lib_Limiters, only: SD_space_order, lim_dev
    use MOSE_Lib_Riemann,  only: rs_dev
    use MOSE_Lib_Preconditioning, only: ur_f, Prec_Alloc
    use openacc, only: acc_get_cuda_stream, acc_async_sync, acc_handle_kind
    integer, intent(in) :: b
    integer, dimension(3), intent(in) :: n
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: Res
    real(R8), dimension(4, 0:n(1), n(2), n(3)), intent(in) :: fm1
    real(R8), dimension(4, n(1), 0:n(2), n(3)), intent(in) :: fm2
    real(R8), dimension(4, n(1), n(2), 0:n(3)), intent(in) :: fm3
    real(R8), dimension(3, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: dl
    real(R8), dimension(1:n(1), 1:n(2), 1:n(3)), intent(inout) :: beta
    integer :: i, j, k, s, f, boff, pi, pj, nfacemax, n1, n2, n3, istat, c0
    integer(kind=cuda_stream_kind) :: strm
    type(dim3) :: grid, tblock
    logical :: SD
    integer :: limv, rsv, nur
    integer :: o1, o2, o3, h1, h2, h3, nfc_
    call cuf_init_tables()
    call Prec_Alloc()   ! ur_f: ncell_total with prec, 1 element otherwise
    nur = size(ur_f)
    n1=n(1); n2=n(2); n3=n(3)
    SD = SD_space_order; limv = lim_dev; rsv = rs_dev
    ! fhat is allocated fhat(nface_max_GLOBAL, nprim) (max over ALL blocks). The CUF
    ! kernel dummy fhat(nfacemax, nprim) MUST use that SAME leading dim as the column
    ! stride, else for a block smaller than the largest the per-var (s>1) offset is
    ! wrong -> energy/species written/read at mismatched locations (the 399x col-5 bug).
    nfacemax = size(fhat, 1)
    call field_block_index_params(b, boff, pi, pj)
    tblock = dim3(128, 1, 1)   ! flat 1D launches
    ! Stream-ordered CUF: launch pass-1 on the OpenACC sync-queue's CUDA stream, so
    ! kernel -> pass-2 scatter ordering is guaranteed by the stream and the per-direction
    ! host barriers (`!$acc wait` + cudaDeviceSynchronize, 6/block/substage here) are
    ! compiled out. Same kernels, same execution order -> bit-exact. ASSUMES one host
    ! thread drives the flux phase (omp single + production OMP_NUM_THREADS=1): OpenACC
    ! queues are per host thread, so a multi-thread driver would break the ordering.
    strm = acc_get_cuda_stream(int(acc_async_sync, acc_handle_kind))

    ! ---- i-direction (flat) ----
    o1 = 1; h1 = n1-1; o2 = 1; h2 = n2; o3 = 1; h3 = n3
    nfc_ = (h1-o1+1)*(h2-o2+1)*(h3-o3+1)
    if (nfc_ > 0) then
    grid = dim3((nfc_ + 127)/128, 1, 1)
    !$acc host_data use_device(P_flat, fhat, fm1, dl, beta, ur_f)
    call conv_cuf_i<<<grid,tblock,0,strm>>>(P_flat, fhat, fm1, dl, beta, ncell_total, nfacemax, &
         boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3, SD, limv, rsv, ur_f, nur)
    !$acc end host_data
    ! NO per-direction pass-2 -- the three directions land in fhat/fhat_j/fhat_k
    ! and the ONE fused pass-2 (Flux_Pass2_Fused, Mod_Fluxes) consumes all six buffers.

    ! ---- j-direction (flat) ----
    end if
    o1 = 1; h1 = n1; o2 = 1; h2 = n2-1; o3 = 1; h3 = n3
    nfc_ = (h1-o1+1)*(h2-o2+1)*(h3-o3+1)
    if (nfc_ > 0) then
    grid = dim3((nfc_ + 127)/128, 1, 1)
    !$acc host_data use_device(P_flat, fhat_j, fm2, dl, beta, ur_f)
    call conv_cuf_j<<<grid,tblock,0,strm>>>(P_flat, fhat_j, fm2, dl, beta, ncell_total, nfacemax, &
         boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3, SD, limv, rsv, ur_f, nur)
    !$acc end host_data

    ! ---- k-direction ----
    end if
    o1 = 1; h1 = n1; o2 = 1; h2 = n2; o3 = 1; h3 = n3-1
    nfc_ = (h1-o1+1)*(h2-o2+1)*(h3-o3+1)
    if (nfc_ > 0) then
    grid = dim3((nfc_ + 127)/128, 1, 1)
    !$acc host_data use_device(P_flat, fhat_k, fm3, dl, beta, ur_f)
    call conv_cuf_k<<<grid,tblock,0,strm>>>(P_flat, fhat_k, fm3, dl, beta, ncell_total, nfacemax, &
         boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3, SD, limv, rsv, ur_f, nur)
    !$acc end host_data

    end if
  end subroutine Convective_Flux


  ! ---- from Mod_Convective_W.f90 (MOSE_Mod_Convective) ----

  !> One-time guard: the w-convective kernel only supports the device-default scheme
  !> pair (Min-Mod limiter + HLLC Batten Riemann), exactly the pair Convective_Flux's
  !> _OPENACC branch hardwires. Map the existing scheme procedure pointers -> integer
  !> ids set at setup (Assign_Limiter/Assign_RS) and assert them here, replacing the
  !> device-hostile pointer dispatch with a checked integer flag.
  subroutine Convective_assert_scheme ()
    use MOSE_Lib_Limiters, only: lim_dev
    use MOSE_Lib_Riemann,  only: rs_dev
    if (lim_dev < 1 .or. lim_dev > 6) error stop &
      'w-convective: limiter not on the device (minmod, vanleer, vanalbada, mc, superbee or first order)'
    if (rs_dev < 1 .or. rs_dev > 14) error stop &
      'w-convective: Riemann solver not on the device (HLLC, HLLE, SLAU and the eleven of riemann_open_body.inc)'
  end subroutine Convective_assert_scheme

  !> Size + device-create the fhat scratch from the per-block dims. One buffer reused
  !> across the three directions, so size it to the max single-direction face count:
  !> i-dir interior faces of a block = (ni-1)*nj*nk (and the j/k analogues); take the
  !> max over all of those over all local blocks.
  subroutine convective_alloc_fhat ( domain )
    use MOSE_Types_m, only: MOSE_domain_type
    use MOSE_Mod_MPI, only: is_local_block
    type(MOSE_domain_type), intent(in) :: domain
    integer :: b, ni, nj, nk, nf
    if (fhat_resident) return
    nface_max = 0
    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle
      ni = domain%blk(b)%dim(1); nj = domain%blk(b)%dim(2); nk = domain%blk(b)%dim(3)
      nf = max((ni-1)*nj*nk, ni*(nj-1)*nk, ni*nj*(nk-1))
      nface_max = max(nface_max, nf)
    end do
    if (nface_max < 1) nface_max = 1
    allocate(fhat(nface_max, nprim), fhat_j(nface_max, nprim), fhat_k(nface_max, nprim))
    fhat = 0.0_R8; fhat_j = 0.0_R8; fhat_k = 0.0_R8
    !$acc enter data create(fhat, fhat_j, fhat_k)
    fhat_resident = .true.
  end subroutine convective_alloc_fhat



end module MOSE_Mod_Convective
