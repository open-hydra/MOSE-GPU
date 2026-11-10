!> GPU CUDA-Fortran prototype — DIFFUSIVE flux pass-1 as a CUF attributes(global)
!> kernel (one face per thread, 3D grid), to break the OpenACC launch_bounds ceiling
!> on the register/spill-bound diffusive kernel (#1 GPU cost; OpenACC pass-1 = 163 reg
!> + 3848 B stack spill). Same playbook as the convective CUF: the WHOLE physics chain
!> (Diffusive_Flux_cached + Tangential_Gradient + Compute_Diffusive_Flux_cached + the
!> SST eddy-viscosity / stress / RANS-flux) is FULLY INLINED into diff_cuf_i/j/k via
!> diff_decls.inc + diff_physics.inc — NO device-routine calls, NO array args (that is
!> what corrupts edge threads under nvfortran). h_tab held as an attributes(device)
!> module copy. Pass-2 (Res accumulation, diffusive sign) stays OpenACC. The face thermo
!> comes from the tface fill (Diffusive_Face_Thermo_Fill): rho/T/mil/kl/cp recomputed
!> from the face-averaged state, matching the CPU-original Compute_Diffusive_Flux.
!>
!> The production diffusive kernel; the recompute path without the cache was deleted.
module MOSE_Mod_Diffusive
  use cudafor
  use iso_fortran_env, only: R8 => real64
  use MOSE_Global_m,   only: NPRIM_DEV
  use FLINT_Lib_Chemistry_data, only: NSMAX
  use MOSE_Mod_SST,    only: sigma_k1, sigma_k2, sigma_w1, sigma_w2, beta_star, alpha_1
  use MOSE_Mod_Spalart, only: cv1, sa_sigma => sigma   ! SA eddy viscosity + nu-tilde diffusive flux (eddy_dev = 2)
  use MOSE_Mod_Wilcox2006, only: wx_sigma_k => sigma_k, wx_sigma_w => sigma_w, wx_beta_star => beta_star, &
                                 wx_Clim => Clim     ! Wilcox2006 eddy viscosity + k/omega diffusive flux (eddy_dev = 3)
  use MOSE_Mod_SSGLRR, only: C_mu, sigmaw_e, sigmad_e
  use MOSE_Mod_Fluid,  only: negativetwothird

  use MOSE_Global_m,   only: gc, NPRIM_DEV, VISC
  use MOSE_Field_m,    only: P_flat, field_block_index_params
  implicit none
  private
  public :: Diffusive_Flux, Diffusive_Flux_List, Diffusive_Flux_Pack
  public :: Compute_Diffusive_Flux
  public :: Diffusive_Face_State, Compute_Diffusive_Flux_St   ! mr-conn-fast (Lib_BC BC_Connection_Visc_fast)

  integer, parameter :: NTMAX = 10000
  ! THOR mechanism compile-time dims: per-primitive local arrays scalarize to
  ! registers instead of spilling to local memory. One binary per species count.
  include 'mose_ct_dims.inc'
  integer, parameter :: NS_TAB = max(NSMAX, nsc)   ! the table follows the compile-time nsc (see Mod_Convective_W_cuf)
  real(R8), device, save :: h_tab_d(0:NTMAX, NS_TAB)
  logical,  save :: dcuf_tables_ready = .false.

  public :: diffusive_alloc_fhat
  public :: dfhat, dfhat_j, dfhat_k   ! exposed for the CUF launcher + the fused pass-2 (Mod_Fluxes)
  public :: tface, Diffusive_Face_Thermo_Fill
  real(R8), allocatable, save :: dfhat(:,:), dfhat_j(:,:), dfhat_k(:,:)
  real(R8), allocatable, save :: tface(:,:)
  logical,  save :: tface_resident = .false.
  integer,  save :: nface_max = 0
  logical,  save :: dfhat_resident = .false.

contains

  subroutine dcuf_init_tables ()
    use FLINT_Lib_Thermodynamic, only: ns, Tmin, Tmax, h_tab
    if (dcuf_tables_ready) return
    if (Tmax > NTMAX) error stop 'MOSE_DIFF_W_CUF: table Tmax exceeds NTMAX'
    if (ns > NS_TAB) error stop 'MOSE_DIFF_W_CUF: species count exceeds the device table size NS_TAB: the binary is built per species count'
    !$acc update host(h_tab)
    h_tab_d(Tmin:Tmax, 1:ns) = h_tab(Tmin:Tmax, 1:ns)
    dcuf_tables_ready = .true.
  end subroutine dcuf_init_tables

  ! ============================================================================
  ! CUF PASS-1 KERNELS (one face per thread, 3D grid). 10-cell stencil.
  ! ============================================================================
  ! The three pass-1 kernels live in diff_cuf_kernels.inc and are instantiated twice: the production
  ! kernels diff_cuf_* (MOSE_MCDIFF_KERNEL undefined: the previous source token for token) and the twins diff_cuf_*_mc
  ! (Dm = the face D_k of tface(:,6:5+nsc) + mu_t/(rho Sct), [MOSE-Physics] schmidt <= 0). The launchers pick on mcdiff_dev.
#define DCUF_KI diff_cuf_i
#define DCUF_KJ diff_cuf_j
#define DCUF_KK diff_cuf_k
#include "diff_cuf_kernels.inc"
#undef DCUF_KI
#undef DCUF_KJ
#undef DCUF_KK
#define MOSE_MCDIFF_KERNEL
#define DCUF_KI diff_cuf_i_mc
#define DCUF_KJ diff_cuf_j_mc
#define DCUF_KK diff_cuf_k_mc
#include "diff_cuf_kernels.inc"
#undef MOSE_MCDIFF_KERNEL
#undef DCUF_KI
#undef DCUF_KJ
#undef DCUF_KK
  ! The multirate list twins (diff_cuf_list_kernel.inc): one kernel for the three directions of a level's face lists.
! the block-pack twins (MOSE_Pack_m): diff_cuf_*_pack and diff_cuf_*_pack_mc, same two instantiations
#define DCUF_PKI diff_cuf_i_pack
#define DCUF_PKJ diff_cuf_j_pack
#define DCUF_PKK diff_cuf_k_pack
#include "diff_cuf_pack_kernels.inc"
#undef DCUF_PKI
#undef DCUF_PKJ
#undef DCUF_PKK
#define MOSE_MCDIFF_KERNEL
#define DCUF_PKI diff_cuf_i_pack_mc
#define DCUF_PKJ diff_cuf_j_pack_mc
#define DCUF_PKK diff_cuf_k_pack_mc
#include "diff_cuf_pack_kernels.inc"
#undef MOSE_MCDIFF_KERNEL
#undef DCUF_PKI
#undef DCUF_PKJ
#undef DCUF_PKK

#define DCUF_KL diff_cuf_list
#include "diff_cuf_list_kernel.inc"
#undef DCUF_KL
#define MOSE_MCDIFF_KERNEL
#define DCUF_KL diff_cuf_list_mc
#include "diff_cuf_list_kernel.inc"
#undef MOSE_MCDIFF_KERNEL
#undef DCUF_KL

  ! ============================================================================
  ! HOST LAUNCHER — CUF pass-1 (host_data) + OpenACC pass-2 (diffusive sign).
  ! Replaces Diffusive_Flux_w_blk when MOSE_DIFF_W_CUF is ON.
  ! ============================================================================
  subroutine Diffusive_Flux ( b, Res, fm1, fm2, fm3, M, yn, rho_c, T_c, &
                                     n, lo, hi )
    use MOSE_Global_m,   only: gc, nprim, nsc, nu, nv, nw, np, nt, nrans, RANS, VISC
    use MOSE_Field_m,    only: P_flat, R_flat, field_block_index_params, ncell_total, rho_c_f, T_c_f
    use MOSE_Lib_RANS,   only: sc_dev, sct_dev, prt_dev, stress_dev, eddy_dev, mcdiff_dev
    use MOSE_Lib_Multirate_Base, only: mr_q_flux
    use openacc, only: acc_get_cuda_stream, acc_async_sync, acc_handle_kind
    integer, intent(in) :: b
    integer, dimension(3), intent(in) :: n
    integer, dimension(3), intent(in) :: lo, hi   ! the cells whose faces are computed: the block, or a multirate box
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: Res
    real(R8), dimension(4, 0:n(1), n(2), n(3)), intent(in) :: fm1
    real(R8), dimension(4, n(1), 0:n(2), n(3)), intent(in) :: fm2
    real(R8), dimension(4, n(1), n(2), 0:n(3)), intent(in) :: fm3
    real(R8), dimension(3, 3, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: M
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: yn, rho_c, T_c
    integer :: i, j, k, s, f, boff, pi, pj, nfacemax, n1, n2, n3, istat, c0
    integer(kind=cuda_stream_kind) :: strm
    type(dim3) :: grid, tblock
    integer :: o1, o2, o3, h1, h2, h3, nfc_
    if (.not. VISC) return
    call dcuf_init_tables()
    n1=n(1); n2=n(2); n3=n(3)
    nfacemax = size(dfhat, 1)
    call field_block_index_params(b, boff, pi, pj)
    tblock = dim3(128, 1, 1)   ! flat 1D launches
    ! Stream-ordered CUF: see the twin comment in Convective_Flux -- pass-1 launches
    ! on the OpenACC sync-queue's stream, the 6 host barriers here are compiled out. In a multirate stage the stage
    ! queue's stream (mr-async-flux, Lib_Multirate_Base mr_q_flux), with the face thermo fills on that queue.
    strm = acc_get_cuda_stream(int(mr_q_flux(), acc_handle_kind))

    ! ---- i-direction ----
    o1 = max(lo(1)-1, 1); h1 = min(hi(1), n1-1); o2 = lo(2); h2 = hi(2); o3 = lo(3); h3 = hi(3)
    nfc_ = (h1-o1+1)*(h2-o2+1)*(h3-o3+1)
    if (nfc_ > 0) then
    grid = dim3((nfc_ + 127)/128, 1, 1)
    call Diffusive_Face_Thermo_Fill ( b, 1, n, o1, h1, o2, h2, o3, h3 )   ! same queue as pass 1: stream-ordered before it reads tface
    !$acc host_data use_device(P_flat, dfhat, fm1, M, yn, rho_c_f, T_c_f, tface)
    if (mcdiff_dev) then
    call diff_cuf_i_mc<<<grid,tblock,0,strm>>>(P_flat, dfhat, fm1, M, yn, rho_c_f, T_c_f, tface, &
         ncell_total, nfacemax, boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3, &
         RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
    else
    call diff_cuf_i<<<grid,tblock,0,strm>>>(P_flat, dfhat, fm1, M, yn, rho_c_f, T_c_f, tface, &
         ncell_total, nfacemax, boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3, &
         RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
    end if
    !$acc end host_data
    ! NO per-direction pass-2 -- the three directions land in dfhat/dfhat_j/dfhat_k
    ! and the ONE fused pass-2 (Flux_Pass2_Fused, Mod_Fluxes) consumes all six buffers.
    ! (tface stays ONE buffer: consumed by the pass-1 launch of ITS direction, above.)

    ! ---- j-direction ----
    end if
    o1 = lo(1); h1 = hi(1); o2 = max(lo(2)-1, 1); h2 = min(hi(2), n2-1); o3 = lo(3); h3 = hi(3)
    nfc_ = (h1-o1+1)*(h2-o2+1)*(h3-o3+1)
    if (nfc_ > 0) then
    grid = dim3((nfc_ + 127)/128, 1, 1)
    call Diffusive_Face_Thermo_Fill ( b, 2, n, o1, h1, o2, h2, o3, h3 )
    !$acc host_data use_device(P_flat, dfhat_j, fm2, M, yn, rho_c_f, T_c_f, tface)
    if (mcdiff_dev) then
    call diff_cuf_j_mc<<<grid,tblock,0,strm>>>(P_flat, dfhat_j, fm2, M, yn, rho_c_f, T_c_f, tface, &
         ncell_total, nfacemax, boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3, &
         RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
    else
    call diff_cuf_j<<<grid,tblock,0,strm>>>(P_flat, dfhat_j, fm2, M, yn, rho_c_f, T_c_f, tface, &
         ncell_total, nfacemax, boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3, &
         RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
    end if
    !$acc end host_data

    ! ---- k-direction ----
    end if
    o1 = lo(1); h1 = hi(1); o2 = lo(2); h2 = hi(2); o3 = max(lo(3)-1, 1); h3 = min(hi(3), n3-1)
    nfc_ = (h1-o1+1)*(h2-o2+1)*(h3-o3+1)
    if (nfc_ > 0) then
    grid = dim3((nfc_ + 127)/128, 1, 1)
    call Diffusive_Face_Thermo_Fill ( b, 3, n, o1, h1, o2, h2, o3, h3 )
    !$acc host_data use_device(P_flat, dfhat_k, fm3, M, yn, rho_c_f, T_c_f, tface)
    if (mcdiff_dev) then
    call diff_cuf_k_mc<<<grid,tblock,0,strm>>>(P_flat, dfhat_k, fm3, M, yn, rho_c_f, T_c_f, tface, &
         ncell_total, nfacemax, boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3, &
         RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
    else
    call diff_cuf_k<<<grid,tblock,0,strm>>>(P_flat, dfhat_k, fm3, M, yn, rho_c_f, T_c_f, tface, &
         ncell_total, nfacemax, boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3, &
         RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
    end if
    !$acc end host_data

    end if
  end subroutine Diffusive_Flux






  ! ---- from Mod_Diffusive_W.f90 (MOSE_Mod_Diffusive) ----

  !> Size + device-create the dfhat scratch. Identical sizing to convective_alloc_fhat:
  !> max single-direction interior-face count over all local blocks.
  subroutine diffusive_alloc_fhat ( domain )
    use MOSE_Types_m, only: MOSE_domain_type
    use MOSE_Mod_MPI, only: is_local_block
    type(MOSE_domain_type), intent(in) :: domain
    integer :: b, ni, nj, nk, nf
    if (dfhat_resident) return
    nface_max = 0
    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle
      ni = domain%blk(b)%dim(1); nj = domain%blk(b)%dim(2); nk = domain%blk(b)%dim(3)
      nf = max((ni-1)*nj*nk, ni*(nj-1)*nk, ni*nj*(nk-1))
      nface_max = max(nface_max, nf)
    end do
    if (nface_max < 1) nface_max = 1
    allocate(dfhat(nface_max, nprim), dfhat_j(nface_max, nprim), dfhat_k(nface_max, nprim))
    dfhat = 0.0_R8; dfhat_j = 0.0_R8; dfhat_k = 0.0_R8
    !$acc enter data create(dfhat, dfhat_j, dfhat_k)
        dfhat_resident = .true.
  end subroutine diffusive_alloc_fhat


  !> Fill tface for ONE direction of ONE block: the face thermo (rho, T, mil, kl, cp)
  !> recomputed from the face-averaged primitives -- statements verbatim from
  !> Compute_Diffusive_Flux's Option-B derivation (Mod_Diffusive.f90:311-324), so
  !> feeding these into the cached consumption restores the legacy closure exactly.
  !> Only the species partial densities (1:nsc) and pressure (np) enter the
  !> derivation, so only those are averaged. Faces are strictly interior (no ghost
  !> reads). Lazy alloc: nface_max is set by diffusive_alloc_fhat, which Fluxes
  !> runs before any per-block launcher reaches here.
  subroutine Diffusive_Face_Thermo_Fill ( b, dir, n, o1, h1, o2, h2, o3, h3 )
    use MOSE_Global_m, only: nsc, np
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, f_cp_dev, co_k_mi_lam_Wilke_dev, co_DS_expr_dev
    use MOSE_Lib_RANS, only: mcdiff_dev, prl_dev
    use MOSE_Lib_Multirate_Base, only: mr_q_flux
    integer, intent(in) :: b, dir
    integer, dimension(3), intent(in) :: n
    integer, intent(in) :: o1, h1, o2, h2, o3, h3   ! the faces of direction dir (the launcher's range)
    integer  :: i, j, k, s, f, c1, c2, boff, pi, pj, T_i, Tint(2), q_
    real(R8) :: Prim(NPRIM_DEV), rho, Rgas, T, Tdiff, cp, mil, kl
    real(R8) :: Dk(NPRIM_DEV)   ! DIFF_MC: the face D_k (co_DS_expr_dev), columns 6:5+nsc
    logical  :: mcd_
    ! Dims as scalars -- a dummy array referenced inside the acc region
    ! is implicitly copy()d per launch (upload+download+sync); scalars ride the
    ! launch parameters for free.
    integer :: nb1_, nb2_, nb3_
    nb1_ = n(1); nb2_ = n(2); nb3_ = n(3)
    if (.not. tface_resident) then
      if (nface_max < 1) error stop 'Diffusive_Face_Thermo_Fill: alloc_fhat must run first'
      if (mcdiff_dev) then
        allocate(tface(nface_max, 5 + nsc))
        write(*,'(A,I0,A,I0,A,F9.2,A)') ' [DIFF_MC] tface ', nface_max, ' faces x ', 5 + nsc, ' columns = ', &
              real(nface_max,R8)*(5 + nsc)*8/1048576.0_R8, ' MB (this rank; 5 columns without the key)'
      else
      allocate(tface(nface_max, 5))
      end if
      tface = 0.0_R8
      !$acc enter data create(tface)
      tface_resident = .true.
    end if
    mcd_ = mcdiff_dev
    q_ = mr_q_flux()   ! the queue of the pass-1 launch that reads tface
    call field_block_index_params(b, boff, pi, pj)
    select case (dir)
    case (1)
      !$acc parallel loop collapse(3) gang vector present(P_flat, tface) async(q_) &
      !$acc&  private(Prim, Tint, s, f, c1, c2, T_i, rho, Rgas, T, Tdiff, cp, mil, kl, Dk)
      do k = o3, h3
      do j = o2, h2
      do i = o1, h1
        c1 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
        c2 = c1 + 1
        do s = 1, nsc
          Prim(s) = 0.5_R8 * ( P_flat(c1,s) + P_flat(c2,s) )
        end do
        Prim(np) = 0.5_R8 * ( P_flat(c1,np) + P_flat(c2,np) )
        call co_rotot_Rtot_dev ( Prim(1:nsc), rho, Rgas )
        T = Prim(np) / ( rho * Rgas )
        T_i = int(T)
        Tdiff = T - T_i
        Tint(1) = T_i
        Tint(2) = T_i + 1
        cp = f_cp_dev ( Prim(1:nsc), Tint, Tdiff, rho )
        call co_k_mi_lam_Wilke_dev ( Prim(1:nsc), rho, Tint, Tdiff, mil, kl )
        f = (i-1) + (j-1)*(nb1_-1) + (k-1)*(nb1_-1)*nb2_ + 1
        if (prl_dev > 0d0) kl = mil*cp/prl_dev   ! [MOSE-Physics] prandtl
        tface(f,1) = rho; tface(f,2) = T; tface(f,3) = mil; tface(f,4) = kl; tface(f,5) = cp
        if (mcd_) then
          call co_DS_expr_dev ( Prim(1:nsc), rho, Tint, Tdiff, Prim(np), Dk(1:nsc) )
          do s = 1, nsc
            tface(f,5+s) = Dk(s)
          end do
        end if
      end do ; end do ; end do
    case (2)
      !$acc parallel loop collapse(3) gang vector present(P_flat, tface) async(q_) &
      !$acc&  private(Prim, Tint, s, f, c1, c2, T_i, rho, Rgas, T, Tdiff, cp, mil, kl, Dk)
      do k = o3, h3
      do j = o2, h2
      do i = o1, h1
        c1 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
        c2 = c1 + pi
        do s = 1, nsc
          Prim(s) = 0.5_R8 * ( P_flat(c1,s) + P_flat(c2,s) )
        end do
        Prim(np) = 0.5_R8 * ( P_flat(c1,np) + P_flat(c2,np) )
        call co_rotot_Rtot_dev ( Prim(1:nsc), rho, Rgas )
        T = Prim(np) / ( rho * Rgas )
        T_i = int(T)
        Tdiff = T - T_i
        Tint(1) = T_i
        Tint(2) = T_i + 1
        cp = f_cp_dev ( Prim(1:nsc), Tint, Tdiff, rho )
        call co_k_mi_lam_Wilke_dev ( Prim(1:nsc), rho, Tint, Tdiff, mil, kl )
        f = (i-1) + (j-1)*nb1_ + (k-1)*nb1_*(nb2_-1) + 1
        if (prl_dev > 0d0) kl = mil*cp/prl_dev   ! [MOSE-Physics] prandtl
        tface(f,1) = rho; tface(f,2) = T; tface(f,3) = mil; tface(f,4) = kl; tface(f,5) = cp
        if (mcd_) then
          call co_DS_expr_dev ( Prim(1:nsc), rho, Tint, Tdiff, Prim(np), Dk(1:nsc) )
          do s = 1, nsc
            tface(f,5+s) = Dk(s)
          end do
        end if
      end do ; end do ; end do
    case (3)
      !$acc parallel loop collapse(3) gang vector present(P_flat, tface) async(q_) &
      !$acc&  private(Prim, Tint, s, f, c1, c2, T_i, rho, Rgas, T, Tdiff, cp, mil, kl, Dk)
      do k = o3, h3
      do j = o2, h2
      do i = o1, h1
        c1 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
        c2 = c1 + pi*pj
        do s = 1, nsc
          Prim(s) = 0.5_R8 * ( P_flat(c1,s) + P_flat(c2,s) )
        end do
        Prim(np) = 0.5_R8 * ( P_flat(c1,np) + P_flat(c2,np) )
        call co_rotot_Rtot_dev ( Prim(1:nsc), rho, Rgas )
        T = Prim(np) / ( rho * Rgas )
        T_i = int(T)
        Tdiff = T - T_i
        Tint(1) = T_i
        Tint(2) = T_i + 1
        cp = f_cp_dev ( Prim(1:nsc), Tint, Tdiff, rho )
        call co_k_mi_lam_Wilke_dev ( Prim(1:nsc), rho, Tint, Tdiff, mil, kl )
        f = (i-1) + (j-1)*nb1_ + (k-1)*nb1_*nb2_ + 1
        if (prl_dev > 0d0) kl = mil*cp/prl_dev   ! [MOSE-Physics] prandtl
        tface(f,1) = rho; tface(f,2) = T; tface(f,3) = mil; tface(f,4) = kl; tface(f,5) = cp
        if (mcd_) then
          call co_DS_expr_dev ( Prim(1:nsc), rho, Tint, Tdiff, Prim(np), Dk(1:nsc) )
          do s = 1, nsc
            tface(f,5+s) = Dk(s)
          end do
        end if
      end do ; end do ; end do
    end select
  end subroutine Diffusive_Face_Thermo_Fill


  !> The face thermo of direction dir of batch ib of the block pack (MOSE_Pack_m): the statements of
  !> Diffusive_Face_Thermo_Fill on the faces of every block of the batch, one gang of 128 per CUDA-block entry of the
  !> batch table (the face is the one the pack flux kernels read at the same row of tfp). Sync construct on the default
  !> queue, in order with the CUF launches of Diffusive_Flux_Pack on the sync queue's stream.
  subroutine Diffusive_Face_Thermo_Fill_Pack ( ib, dir )
    use MOSE_Global_m, only: nsc, np, gc
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, f_cp_dev, co_k_mi_lam_Wilke_dev, co_DS_expr_dev
    use MOSE_Lib_RANS, only: mcdiff_dev, prl_dev
    use MOSE_Pack_m
    integer, intent(in) :: ib, dir
    integer  :: e, tt, t_, s, f, c1, c2, T_i, Tint(2), p, n1, n2, n3, boff, pi, pj, i, j, k, ni_, nj_, nk_, rem_, cb0, ncb, sa, of
    real(R8) :: Prim(NPRIM_DEV), rho, Rgas, T, Tdiff, cp, mil, kl
    real(R8) :: Dk(NPRIM_DEV)
    logical  :: mcd_
    mcd_ = mcdiff_dev
    cb0 = cb_n(0, dir, ib); ncb = cb_n(1, dir, ib)
    if ( ncb == 0 ) return
    !$acc parallel loop gang vector_length(128) present(P_flat, tfp, cb_p, cb_t0, bt_n1, bt_n2, bt_n3, bt_boff, bt_pi, bt_pj, bt_ofh) &
    !$acc&  private(p, n1, n2, n3, boff, pi, pj, ni_, nj_, nk_, sa, of)
    do e = 1, ncb
      p = cb_p(cb0+e)
      n1 = bt_n1(p); n2 = bt_n2(p); n3 = bt_n3(p); boff = bt_boff(p); pi = bt_pi(p); pj = bt_pj(p)
      of = bt_ofh(dir, p)
      if ( dir == 1 ) then
        ni_ = n1 - 1; nj_ = n2; nk_ = n3; sa = 1
      else if ( dir == 2 ) then
        ni_ = n1; nj_ = n2 - 1; nk_ = n3; sa = pi
      else
        ni_ = n1; nj_ = n2; nk_ = n3 - 1; sa = pi*pj
      end if
      !$acc loop vector private(Prim, Tint, Dk, t_, s, f, c1, c2, T_i, rho, Rgas, T, Tdiff, cp, mil, kl, i, j, k, rem_)
      do tt = 1, PACK_NT
        t_ = cb_t0(cb0+e) + tt
        if ( t_ <= ni_*nj_*nk_ ) then
          i    = mod(t_-1, ni_) + 1
          rem_ = (t_-1) / ni_
          j    = mod(rem_, nj_) + 1
          k    = rem_ / nj_ + 1
          c1 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
          c2 = c1 + sa
          do s = 1, nsc
            Prim(s) = 0.5_R8 * ( P_flat(c1,s) + P_flat(c2,s) )
          end do
          Prim(np) = 0.5_R8 * ( P_flat(c1,np) + P_flat(c2,np) )
          call co_rotot_Rtot_dev ( Prim(1:nsc), rho, Rgas )
          T = Prim(np) / ( rho * Rgas )
          T_i = int(T)
          Tdiff = T - T_i
          Tint(1) = T_i
          Tint(2) = T_i + 1
          cp = f_cp_dev ( Prim(1:nsc), Tint, Tdiff, rho )
          call co_k_mi_lam_Wilke_dev ( Prim(1:nsc), rho, Tint, Tdiff, mil, kl )
          if ( dir == 1 ) then
            f = of + (i-1) + (j-1)*(n1-1) + (k-1)*(n1-1)*n2 + 1
          else if ( dir == 2 ) then
            f = of + (i-1) + (j-1)*n1 + (k-1)*n1*(n2-1) + 1
          else
            f = of + (i-1) + (j-1)*n1 + (k-1)*n1*n2 + 1
          end if
          if (prl_dev > 0d0) kl = mil*cp/prl_dev   ! [MOSE-Physics] prandtl
          tfp(f,1) = rho; tfp(f,2) = T; tfp(f,3) = mil; tfp(f,4) = kl; tfp(f,5) = cp
          if (mcd_) then
            call co_DS_expr_dev ( Prim(1:nsc), rho, Tint, Tdiff, Prim(np), Dk(1:nsc) )
            do s = 1, nsc
              tfp(f,5+s) = Dk(s)
            end do
          end if
        end if
      end do
    end do
  end subroutine Diffusive_Face_Thermo_Fill_Pack

  !> The diffusive pass 1 of batch ib of the block pack (MOSE_Pack_m): per direction the face thermo of the batch
  !> (Diffusive_Face_Thermo_Fill_Pack), then the pack flux kernel (diff_cuf_*_pack, or the _mc twin), on the sync queue's
  !> stream as Diffusive_Flux. Six launches for the batch.
  subroutine Diffusive_Flux_Pack ( ib )
    use MOSE_Global_m,   only: gc, nprim, RANS, VISC
    use MOSE_Field_m,    only: P_flat, ncell_total, rho_c_f, T_c_f
    use MOSE_Lib_RANS,   only: sc_dev, sct_dev, prt_dev, stress_dev, eddy_dev, mcdiff_dev
    use MOSE_Pack_m
    use openacc, only: acc_get_cuda_stream, acc_async_sync, acc_handle_kind
    integer, intent(in) :: ib
    integer(kind=cuda_stream_kind) :: strm
    type(dim3) :: grid, tblock
    integer :: nl
    if (.not. VISC) return
    call dcuf_init_tables()
    nl = size(bt_n1)
    tblock = dim3(PACK_NT, 1, 1)
    strm = acc_get_cuda_stream(int(acc_async_sync, acc_handle_kind))
    if ( pbat(ib)%ncb(1) > 0 ) then
      grid = dim3(pbat(ib)%ncb(1), 1, 1)
      call Diffusive_Face_Thermo_Fill_Pack ( ib, 1 )
      !$acc host_data use_device(P_flat, dfhp1, fm1p, Mp, ynp, rho_c_f, T_c_f, tfp, cb_p, cb_t0, bt_n1, bt_n2, bt_n3, bt_boff, bt_pi, bt_pj, bt_ofm1, bt_ofh)
      if (mcdiff_dev) then
      call diff_cuf_i_pack_mc<<<grid,tblock,0,strm>>>(P_flat, dfhp1, fm1p, Mp, ynp, rho_c_f, T_c_f, tfp, ncell_total, nfp1, size(fm1p,2), ntfp, &
           cb_p, cb_t0, cb_n(0,1,ib), bt_n1, bt_n2, bt_n3, bt_boff, bt_pi, bt_pj, bt_ofm1, bt_ofh, nl, gc, &
           RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
      else
      call diff_cuf_i_pack<<<grid,tblock,0,strm>>>(P_flat, dfhp1, fm1p, Mp, ynp, rho_c_f, T_c_f, tfp, ncell_total, nfp1, size(fm1p,2), ntfp, &
           cb_p, cb_t0, cb_n(0,1,ib), bt_n1, bt_n2, bt_n3, bt_boff, bt_pi, bt_pj, bt_ofm1, bt_ofh, nl, gc, &
           RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
      end if
      !$acc end host_data
    end if
    if ( pbat(ib)%ncb(2) > 0 ) then
      grid = dim3(pbat(ib)%ncb(2), 1, 1)
      call Diffusive_Face_Thermo_Fill_Pack ( ib, 2 )
      !$acc host_data use_device(P_flat, dfhp2, fm2p, Mp, ynp, rho_c_f, T_c_f, tfp, cb_p, cb_t0, bt_n1, bt_n2, bt_n3, bt_boff, bt_pi, bt_pj, bt_ofm2, bt_ofh)
      if (mcdiff_dev) then
      call diff_cuf_j_pack_mc<<<grid,tblock,0,strm>>>(P_flat, dfhp2, fm2p, Mp, ynp, rho_c_f, T_c_f, tfp, ncell_total, nfp2, size(fm2p,2), ntfp, &
           cb_p, cb_t0, cb_n(0,2,ib), bt_n1, bt_n2, bt_n3, bt_boff, bt_pi, bt_pj, bt_ofm2, bt_ofh, nl, gc, &
           RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
      else
      call diff_cuf_j_pack<<<grid,tblock,0,strm>>>(P_flat, dfhp2, fm2p, Mp, ynp, rho_c_f, T_c_f, tfp, ncell_total, nfp2, size(fm2p,2), ntfp, &
           cb_p, cb_t0, cb_n(0,2,ib), bt_n1, bt_n2, bt_n3, bt_boff, bt_pi, bt_pj, bt_ofm2, bt_ofh, nl, gc, &
           RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
      end if
      !$acc end host_data
    end if
    if ( pbat(ib)%ncb(3) > 0 ) then
      grid = dim3(pbat(ib)%ncb(3), 1, 1)
      call Diffusive_Face_Thermo_Fill_Pack ( ib, 3 )
      !$acc host_data use_device(P_flat, dfhp3, fm3p, Mp, ynp, rho_c_f, T_c_f, tfp, cb_p, cb_t0, bt_n1, bt_n2, bt_n3, bt_boff, bt_pi, bt_pj, bt_ofm3, bt_ofh)
      if (mcdiff_dev) then
      call diff_cuf_k_pack_mc<<<grid,tblock,0,strm>>>(P_flat, dfhp3, fm3p, Mp, ynp, rho_c_f, T_c_f, tfp, ncell_total, nfp3, size(fm3p,2), ntfp, &
           cb_p, cb_t0, cb_n(0,3,ib), bt_n1, bt_n2, bt_n3, bt_boff, bt_pi, bt_pj, bt_ofm3, bt_ofh, nl, gc, &
           RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
      else
      call diff_cuf_k_pack<<<grid,tblock,0,strm>>>(P_flat, dfhp3, fm3p, Mp, ynp, rho_c_f, T_c_f, tfp, ncell_total, nfp3, size(fm3p,2), ntfp, &
           cb_p, cb_t0, cb_n(0,3,ib), bt_n1, bt_n2, bt_n3, bt_boff, bt_pi, bt_pj, bt_ofm3, bt_ofh, nl, gc, &
           RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
      end if
      !$acc end host_data
    end if
  end subroutine Diffusive_Flux_Pack

  !> The diffusive pass 1 of a multirate level on its face lists (Mod_Fluxes Fluxes_List): per direction the face thermo
  !> of the list (Diffusive_Face_Thermo_Fill_List), then the flux kernel, one launch each over every local block.
  subroutine Diffusive_Flux_List ( fl, dsel )
    use MOSE_Lib_Multirate_Base, only: mr_fl_type
    type(mr_fl_type), intent(inout) :: fl(3)
    integer, optional, intent(in) :: dsel   ! one direction only (mr-flux-fanout), else all three
    integer :: d, d0, d1
    if (.not. VISC) return
    call dcuf_init_tables()
    d0 = 1; d1 = 3
    if ( present(dsel) ) then
      d0 = dsel; d1 = dsel
    end if
    do d = d0, d1
      if (fl(d)%n == 0) cycle
      call Diffusive_Face_Thermo_Fill_List ( fl(d)%tf, fl(d)%c, fl(d)%sa, fl(d)%n, size(fl(d)%tf, 2) )
      call diff_list_launch ( fl(d)%dfh, fl(d)%n, d, fl(d)%c, fl(d)%sa, fl(d)%sb, fl(d)%sc, fl(d)%fm, fl(d)%mm, fl(d)%yn, &
                              fl(d)%tf, size(fl(d)%tf, 2) )
    end do
  end subroutine Diffusive_Flux_List

  subroutine diff_list_launch ( dfh, ne, d, c, sa, sb, sc, fmg, mmg, yng, tf, ntf )
    use MOSE_Global_m,   only: nprim, RANS
    use MOSE_Lib_Multirate_Base, only: mr_q
    use MOSE_Field_m,    only: ncell_total, rho_c_f, T_c_f
    use MOSE_Lib_RANS,   only: sc_dev, sct_dev, prt_dev, stress_dev, eddy_dev, mcdiff_dev
    use openacc, only: acc_get_cuda_stream, acc_async_sync, acc_handle_kind
    integer,  intent(in)    :: ne, d, ntf, c(ne), sa(ne), sb(ne), sc(ne)
    real(R8), intent(in)    :: fmg(4, ne), mmg(3, 3, 2, ne), yng(2, ne), tf(ne, ntf)
    real(R8), intent(inout) :: dfh(ne, nprim)
    integer(kind=cuda_stream_kind) :: strm
    type(dim3) :: grid, tblock
    tblock = dim3(128, 1, 1); grid = dim3((ne + 127)/128, 1, 1)
    strm = acc_get_cuda_stream(int(mr_q, acc_handle_kind))   ! the stage queue's stream (mr-async; the sync queue's outside a list stage)
    !$acc host_data use_device(P_flat, dfh, c, sa, sb, sc, fmg, mmg, yng, rho_c_f, T_c_f, tf)
    if (mcdiff_dev) then
    call diff_cuf_list_mc<<<grid,tblock,0,strm>>>(P_flat, dfh, ne, d, c, sa, sb, sc, fmg, mmg, yng, rho_c_f, T_c_f, &
         tf, ncell_total, RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
    else
    call diff_cuf_list<<<grid,tblock,0,strm>>>(P_flat, dfh, ne, d, c, sa, sb, sc, fmg, mmg, yng, rho_c_f, T_c_f, &
         tf, ncell_total, RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
    end if
    !$acc end host_data
  end subroutine diff_list_launch

  !> The face thermo of a face list (Diffusive_Face_Thermo_Fill on the faces of the list: the face of position e between
  !> the cells c(e) and c(e) + sa(e)), the statements of the box fill verbatim.
  subroutine Diffusive_Face_Thermo_Fill_List ( tf, c, sa, ne, ntf )
    use MOSE_Global_m, only: nsc, np
    use MOSE_Lib_Multirate_Base, only: mr_q
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, f_cp_dev, co_k_mi_lam_Wilke_dev, co_DS_expr_dev
    use MOSE_Lib_RANS, only: mcdiff_dev, prl_dev
    integer,  intent(in)    :: ne, ntf, c(ne), sa(ne)
    real(R8), intent(inout) :: tf(ne, ntf)
    integer  :: e, s, c1, c2, T_i, Tint(2)
    real(R8) :: Prim(NPRIM_DEV), rho, Rgas, T, Tdiff, cp, mil, kl
    real(R8) :: Dk(NPRIM_DEV)
    logical  :: mcd_
    mcd_ = mcdiff_dev
    !$acc parallel loop gang vector async(mr_q) present(P_flat, tf, c, sa) &
    !$acc&  private(Prim, Tint, s, c1, c2, T_i, rho, Rgas, T, Tdiff, cp, mil, kl, Dk)
    do e = 1, ne
      c1 = c(e)
      c2 = c1 + sa(e)
      do s = 1, nsc
        Prim(s) = 0.5_R8 * ( P_flat(c1,s) + P_flat(c2,s) )
      end do
      Prim(np) = 0.5_R8 * ( P_flat(c1,np) + P_flat(c2,np) )
      call co_rotot_Rtot_dev ( Prim(1:nsc), rho, Rgas )
      T = Prim(np) / ( rho * Rgas )
      T_i = int(T)
      Tdiff = T - T_i
      Tint(1) = T_i
      Tint(2) = T_i + 1
      cp = f_cp_dev ( Prim(1:nsc), Tint, Tdiff, rho )
      call co_k_mi_lam_Wilke_dev ( Prim(1:nsc), rho, Tint, Tdiff, mil, kl )
      if (prl_dev > 0d0) kl = mil*cp/prl_dev   ! [MOSE-Physics] prandtl
      tf(e,1) = rho; tf(e,2) = T; tf(e,3) = mil; tf(e,4) = kl; tf(e,5) = cp
      if (mcd_) then
        call co_DS_expr_dev ( Prim(1:nsc), rho, Tint, Tdiff, Prim(np), Dk(1:nsc) )
        do s = 1, nsc
          tf(e,5+s) = Dk(s)
        end do
      end if
    end do
  end subroutine Diffusive_Face_Thermo_Fill_List






  ! ---- from Mod_Diffusive.f90 (the host Compute_Diffusive_Flux, an acc routine of the record BC kernels) ----





    subroutine Compute_Diffusive_Flux ( Prim, Gradient, area, normal, waldis, Flux )
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic
      !$acc routine seq
      use MOSE_Lib_RANS
      use MOSE_Mod_Fluid
      use MOSE_IO_RANS, only: obj_rans
      ! Device build: SOOT is unsupported (guarded out below); the diffusion
      ! non-dim numbers, turbulence ids and concrete SST routines come in via
      ! device-resident scalars / direct (devirtualized) calls.
      use MOSE_Mod_SST, only: SST_Eddy_Viscosity, SST_RANS_Diffusive_Flux
      use MOSE_Mod_Spalart, only: Spalart_Eddy_Viscosity, Spalart_RANS_Diffusive_Flux
      use MOSE_Mod_Wilcox2006, only: Wilcox2006_Eddy_Viscosity, Wilcox2006_RANS_Diffusive_Flux
      use MOSE_Lib_RANS, only: eddy_dev, prl_dev
      use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, f_cp_dev, &
                                               co_k_mi_lam_Wilke_dev, co_DS_expr_dev

      implicit none
      real(R8), intent(in) :: Prim(NPRIM_DEV), Gradient(NPRIM_DEV,3), area, normal(3), waldis
      real(R8), intent(out) :: Flux(NPRIM_DEV)
      ! Local
      integer :: s, T_i, Tint(2)
      real(R8) :: rho, Rgas, T, cp, mil, kl
      real(R8) :: Tdiff, mie, mu, kappa
      real(R8) :: Dm(NPRIM_DEV), stress(3), DiffHFlux, DmGradYi, hh, hh1, hh2
      ! Contiguous gathers of the velocity/RANS gradient blocks. Passing the
      ! strided sections Gradient(nu:nw,:) / Gradient(nt:nprim,:) (rows of the
      ! NPRIM_DEV-padded array -> stride-NPRIM_DEV columns) straight into the
      ! `acc routine seq` SST/stress dummies makes nvfortran misread them on the
      ! device (host gathers silently; device codegen does not). Gkw leading dim
      ! is pinned to 2 = nRANS for the device-only SST model (ransdiff_dev=1).
      real(R8) :: Gvel(3,3), Gkw(2,3)
      integer  :: dd

      ! Option B: thermo properties at the face are computed from the
      ! face-averaged primitive Prim — bit-identical to the pre-thermo-cache
      ! (3594321) path. The cell-centered cache values from Mod_Thermo_Cache
      ! are not consulted here; this avoids the O(Δx²) drift introduced by
      ! face-averaging cell-computed thermo properties.
      ! Device: general-ns thermo can't run general-ns on the GPU -> fixed-size
      ! _dev copies. Turbulence model is devirtualized (device supports SST only).
      call co_rotot_Rtot_dev ( Prim(1:nsc), rho, Rgas )
      T = Prim(np) / ( rho * Rgas )

      T_i     = int(T)
      Tdiff   = T - T_i
      Tint(1) = T_i
      Tint(2) = T_i + 1

      cp = f_cp_dev ( Prim(1:nsc), Tint, Tdiff, rho )
      call co_k_mi_lam_Wilke_dev ( Prim(1:nsc), rho, Tint, Tdiff, mil, kl )
      if (prl_dev > 0d0) kl = mil*cp/prl_dev   ! [MOSE-Physics] prandtl (MOSE_open Lib_Diffusive)

      ! Eddy viscosity
      mie = 0d0
      ! Gather the velocity-gradient (3x3) into a contiguous local once.
      do dd = 1, 3
        Gvel(1,dd) = Gradient(nu,dd); Gvel(2,dd) = Gradient(nv,dd); Gvel(3,dd) = Gradient(nw,dd)
      end do
      if (RANS) then
        if (eddy_dev == 2) then
          call Spalart_Eddy_Viscosity ( mie, Prim(nt:nprim), mil, rho, Gvel, waldis )
        else if (eddy_dev == 3) then
          call Wilcox2006_Eddy_Viscosity ( mie, Prim(nt:nprim), mil, rho, Gvel, waldis )
        else
          call SST_Eddy_Viscosity ( mie, Prim(nt:nprim), mil, rho, Gvel, waldis )
        end if
      end if

      if (mcdiff_dev) then
        ! DIFF_MC: mixture-averaged D_k (table, rescaled to the local p) + turbulent constant-Sct part, as MOSE_open
        call co_DS_expr_dev ( Prim(1:nsc), rho, Tint, Tdiff, Prim(np), Dm(1:nsc) )
        if (mie > 0d0) Dm(1:nsc) = Dm(1:nsc) + mie/(rho*sct_dev)
      else
      Dm (1:nsc) = ( mil/sc_dev + mie/sct_dev ) / rho
      end if
      kappa = kl + mie*cp/prt_dev
      ! Devirtualized stress tensor (1 = Std/laminar, 2 = 2-eq k-coupled SST)
      if (stress_dev == 2) then
        Stress = Stress_Vector_2eq ( Gvel, normal, mil, mie, prim(nt:nprim) )
      else
        Stress = Stress_Vector_Std ( Gvel, normal, mil, mie, prim(nt:nprim) )
      end if

      ! Fluxes computation
      DiffHFlux = 0.0d0
      Flux = 0.0d0

      ! Diffusive mass fluxes — h_tab lookup inlined so the species loop
      ! vectorizes (avoids function-call dispatch per species per face).
      do s = 1, nsc
        DmGradYi  = Dm(s) * dot_product ( Gradient(s,:), normal )     ! Dm * dyi/dn
        Flux(s) = rho * area * DmGradYi                               ! rho * Dm *dyi/dn * Area
        hh1 = h_tab ( Tint(1), s )
        hh2 = h_tab ( Tint(2), s )
        hh  = hh1 + Tdiff * ( hh2 - hh1 )
        DiffHFlux = DiffHFlux + DmGradYi * hh
      enddo

      ! Correcting the mass flux for possible errors (sum must be zero)
      Flux(1:nsc) = Flux(1:nsc) - sum(Flux(1:nsc))*Prim(1:nsc)/rho

      ! Momentum flux
      Flux(nu:nw) = area * Stress

      ! Diffusive enthalpy flux
      Flux(np) = area * ( dot_product ( Stress, Prim(nu:nw) ) + &
                 kappa * dot_product ( Gradient(np,:), normal ) + DiffHFlux * rho )

      ! Turbulence variables diffusive flux
      if (RANS) then
        if (eddy_dev == 2) then
          do dd = 1, 3
            Gkw(1,dd) = Gradient(nt,dd)
          end do
          call Spalart_RANS_Diffusive_Flux ( Flux(nt:nprim), Prim(nt:nprim), mil, Gvel, Gkw(1:nrans,1:3), rho, area, normal, waldis )
        else
          do dd = 1, 3
            Gkw(1,dd) = Gradient(nt,dd); Gkw(2,dd) = Gradient(min(nt+1,nprim),dd)   ! min(): no constant out-of-range subscript in an nrans=1 binary
          end do
          if (eddy_dev == 3) then
            call Wilcox2006_RANS_Diffusive_Flux ( Flux(nt:nprim), Prim(nt:nprim), mil, Gvel, Gkw, rho, area, normal, waldis )
          else
            call SST_RANS_Diffusive_Flux ( Flux(nt:nprim), Prim(nt:nprim), mil, &
                                           Gvel, Gkw, &
                                           rho, area, normal, waldis )
          end if
        end if
      end if


    end subroutine Compute_Diffusive_Flux


    !> mr-conn-fast (08/11): Compute_Diffusive_Flux in two parts, for BC_Connection_Visc_fast, which calls it twice with the
    !> same face state Prim. This part depends on Prim only (thermo, transport, the laminar mixture diffusion): the lines of
    !> Compute_Diffusive_Flux verbatim, so the two calls share one evaluation and get the same values.
    subroutine Diffusive_Face_State ( Prim, rho, Tint, Tdiff, cp, mil, kl, Dml )
      !$acc routine seq
      use MOSE_Global_m
      use MOSE_Lib_RANS, only: prl_dev, mcdiff_dev
      use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, f_cp_dev, co_k_mi_lam_Wilke_dev, co_DS_expr_dev
      implicit none
      real(R8), intent(in)  :: Prim(NPRIM_DEV)
      real(R8), intent(out) :: rho, Tdiff, cp, mil, kl, Dml(NPRIM_DEV)
      integer,  intent(out) :: Tint(2)
      integer  :: T_i
      real(R8) :: Rgas, T

      call co_rotot_Rtot_dev ( Prim(1:nsc), rho, Rgas )
      T = Prim(np) / ( rho * Rgas )

      T_i     = int(T)
      Tdiff   = T - T_i
      Tint(1) = T_i
      Tint(2) = T_i + 1

      cp = f_cp_dev ( Prim(1:nsc), Tint, Tdiff, rho )
      call co_k_mi_lam_Wilke_dev ( Prim(1:nsc), rho, Tint, Tdiff, mil, kl )
      if (prl_dev > 0d0) kl = mil*cp/prl_dev   ! [MOSE-Physics] prandtl (MOSE_open Lib_Diffusive)
      if (mcdiff_dev) then
        call co_DS_expr_dev ( Prim(1:nsc), rho, Tint, Tdiff, Prim(np), Dml(1:nsc) )
      else
        Dml(1:nsc) = 0d0
      end if
    end subroutine Diffusive_Face_State


    !> mr-conn-fast: the rest of Compute_Diffusive_Flux (from the eddy viscosity on), verbatim, on the face state of
    !> Diffusive_Face_State.
    subroutine Compute_Diffusive_Flux_St ( Prim, Gradient, area, normal, waldis, rho, Tint, Tdiff, cp, mil, kl, Dml, Flux )
      !$acc routine seq
      use MOSE_Global_m
      use FLINT_Lib_Thermodynamic
      use MOSE_Lib_RANS
      use MOSE_Mod_Fluid
      use MOSE_Mod_SST, only: SST_Eddy_Viscosity, SST_RANS_Diffusive_Flux
      use MOSE_Mod_Spalart, only: Spalart_Eddy_Viscosity, Spalart_RANS_Diffusive_Flux
      use MOSE_Mod_Wilcox2006, only: Wilcox2006_Eddy_Viscosity, Wilcox2006_RANS_Diffusive_Flux
      use MOSE_Lib_RANS, only: eddy_dev
      implicit none
      real(R8), intent(in)  :: Prim(NPRIM_DEV), Gradient(NPRIM_DEV,3), area, normal(3), waldis
      real(R8), intent(in)  :: rho, Tdiff, cp, mil, kl, Dml(NPRIM_DEV)
      integer,  intent(in)  :: Tint(2)
      real(R8), intent(out) :: Flux(NPRIM_DEV)
      integer  :: s, dd
      real(R8) :: mie, kappa
      real(R8) :: Dm(NPRIM_DEV), stress(3), DiffHFlux, DmGradYi, hh, hh1, hh2
      real(R8) :: Gvel(3,3), Gkw(2,3)

      ! Eddy viscosity
      mie = 0d0
      do dd = 1, 3
        Gvel(1,dd) = Gradient(nu,dd); Gvel(2,dd) = Gradient(nv,dd); Gvel(3,dd) = Gradient(nw,dd)
      end do
      if (RANS) then
        if (eddy_dev == 2) then
          call Spalart_Eddy_Viscosity ( mie, Prim(nt:nprim), mil, rho, Gvel, waldis )
        else if (eddy_dev == 3) then
          call Wilcox2006_Eddy_Viscosity ( mie, Prim(nt:nprim), mil, rho, Gvel, waldis )
        else
          call SST_Eddy_Viscosity ( mie, Prim(nt:nprim), mil, rho, Gvel, waldis )
        end if
      end if

      if (mcdiff_dev) then
        Dm(1:nsc) = Dml(1:nsc)
        if (mie > 0d0) Dm(1:nsc) = Dm(1:nsc) + mie/(rho*sct_dev)
      else
      Dm (1:nsc) = ( mil/sc_dev + mie/sct_dev ) / rho
      end if
      kappa = kl + mie*cp/prt_dev
      if (stress_dev == 2) then
        Stress = Stress_Vector_2eq ( Gvel, normal, mil, mie, prim(nt:nprim) )
      else
        Stress = Stress_Vector_Std ( Gvel, normal, mil, mie, prim(nt:nprim) )
      end if

      DiffHFlux = 0.0d0
      Flux = 0.0d0

      do s = 1, nsc
        DmGradYi  = Dm(s) * dot_product ( Gradient(s,:), normal )     ! Dm * dyi/dn
        Flux(s) = rho * area * DmGradYi                               ! rho * Dm *dyi/dn * Area
        hh1 = h_tab ( Tint(1), s )
        hh2 = h_tab ( Tint(2), s )
        hh  = hh1 + Tdiff * ( hh2 - hh1 )
        DiffHFlux = DiffHFlux + DmGradYi * hh
      enddo

      Flux(1:nsc) = Flux(1:nsc) - sum(Flux(1:nsc))*Prim(1:nsc)/rho

      Flux(nu:nw) = area * Stress

      Flux(np) = area * ( dot_product ( Stress, Prim(nu:nw) ) + &
                 kappa * dot_product ( Gradient(np,:), normal ) + DiffHFlux * rho )

      if (RANS) then
        if (eddy_dev == 2) then
          do dd = 1, 3
            Gkw(1,dd) = Gradient(nt,dd)
          end do
          call Spalart_RANS_Diffusive_Flux ( Flux(nt:nprim), Prim(nt:nprim), mil, Gvel, Gkw(1:nrans,1:3), rho, area, normal, waldis )
        else
          do dd = 1, 3
            Gkw(1,dd) = Gradient(nt,dd); Gkw(2,dd) = Gradient(min(nt+1,nprim),dd)
          end do
          if (eddy_dev == 3) then
            call Wilcox2006_RANS_Diffusive_Flux ( Flux(nt:nprim), Prim(nt:nprim), mil, Gvel, Gkw, rho, area, normal, waldis )
          else
            call SST_RANS_Diffusive_Flux ( Flux(nt:nprim), Prim(nt:nprim), mil, &
                                           Gvel, Gkw, &
                                           rho, area, normal, waldis )
          end if
        end if
      end if
    end subroutine Compute_Diffusive_Flux_St


end module MOSE_Mod_Diffusive
