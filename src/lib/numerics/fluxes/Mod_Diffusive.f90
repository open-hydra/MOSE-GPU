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
  public :: Diffusive_Flux
  public :: Compute_Diffusive_Flux

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

  ! ============================================================================
  ! HOST LAUNCHER — CUF pass-1 (host_data) + OpenACC pass-2 (diffusive sign).
  ! Replaces Diffusive_Flux_w_blk when MOSE_DIFF_W_CUF is ON.
  ! ============================================================================
  subroutine Diffusive_Flux ( b, Res, fm1, fm2, fm3, M, yn, rho_c, T_c, &
                                     n)
    use MOSE_Global_m,   only: gc, nprim, nsc, nu, nv, nw, np, nt, nrans, RANS, VISC
    use MOSE_Field_m,    only: P_flat, R_flat, field_block_index_params, ncell_total, rho_c_f, T_c_f
    use MOSE_Lib_RANS,   only: sc_dev, sct_dev, prt_dev, stress_dev, eddy_dev, mcdiff_dev
    use openacc, only: acc_get_cuda_stream, acc_async_sync, acc_handle_kind
    integer, intent(in) :: b
    integer, dimension(3), intent(in) :: n
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
    ! on the OpenACC sync-queue's stream, the 6 host barriers here are compiled out.
    strm = acc_get_cuda_stream(int(acc_async_sync, acc_handle_kind))

    ! ---- i-direction ----
    o1 = 1; h1 = n1-1; o2 = 1; h2 = n2; o3 = 1; h3 = n3
    nfc_ = (h1-o1+1)*(h2-o2+1)*(h3-o3+1)
    if (nfc_ > 0) then
    grid = dim3((nfc_ + 127)/128, 1, 1)
    call Diffusive_Face_Thermo_Fill ( b, 1, n)   ! sync-queue kernel: complete (or stream-ordered) before pass-1 reads tface
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
    o1 = 1; h1 = n1; o2 = 1; h2 = n2-1; o3 = 1; h3 = n3
    nfc_ = (h1-o1+1)*(h2-o2+1)*(h3-o3+1)
    if (nfc_ > 0) then
    grid = dim3((nfc_ + 127)/128, 1, 1)
    call Diffusive_Face_Thermo_Fill ( b, 2, n)
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
    o1 = 1; h1 = n1; o2 = 1; h2 = n2; o3 = 1; h3 = n3-1
    nfc_ = (h1-o1+1)*(h2-o2+1)*(h3-o3+1)
    if (nfc_ > 0) then
    grid = dim3((nfc_ + 127)/128, 1, 1)
    call Diffusive_Face_Thermo_Fill ( b, 3, n)
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
  subroutine Diffusive_Face_Thermo_Fill ( b, dir, n)
    use MOSE_Global_m, only: nsc, np
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, f_cp_dev, co_k_mi_lam_Wilke_dev, co_DS_expr_dev
    use MOSE_Lib_RANS, only: mcdiff_dev, prl_dev
    integer, intent(in) :: b, dir
    integer, dimension(3), intent(in) :: n
    integer  :: i, j, k, s, f, c1, c2, boff, pi, pj, T_i, Tint(2)
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
    call field_block_index_params(b, boff, pi, pj)
    select case (dir)
    case (1)
      !$acc parallel loop collapse(3) gang vector present(P_flat, tface) &
      !$acc&  private(Prim, Tint, s, f, c1, c2, T_i, rho, Rgas, T, Tdiff, cp, mil, kl, Dk)
      do k = 1, n(3)
      do j = 1, n(2)
      do i = 1, nb1_-1
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
      !$acc parallel loop collapse(3) gang vector present(P_flat, tface) &
      !$acc&  private(Prim, Tint, s, f, c1, c2, T_i, rho, Rgas, T, Tdiff, cp, mil, kl, Dk)
      do k = 1, n(3)
      do j = 1, nb2_-1
      do i = 1, n(1)
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
      !$acc parallel loop collapse(3) gang vector present(P_flat, tface) &
      !$acc&  private(Prim, Tint, s, f, c1, c2, T_i, rho, Rgas, T, Tdiff, cp, mil, kl, Dk)
      do k = 1, nb3_-1
      do j = 1, n(2)
      do i = 1, n(1)
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
  
  
end module MOSE_Mod_Diffusive
