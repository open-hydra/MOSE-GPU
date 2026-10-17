module MOSE_Mod_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public  :: Fluxes

contains

  subroutine Fluxes ( domain )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Field_m, only: field_is_allocated, field_alloc, pack_P_to_flat, field_update_P_flat_device, R_flat
    use MOSE_Mod_Convective, only: convective_alloc_fhat, Convective_assert_scheme
    use MOSE_Mod_Diffusive,  only: diffusive_alloc_fhat

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b
    logical, save :: w_setup_done = .false.
    logical, save :: wr_init_done = .false.

    ! Run the CONVECTIVE flux on the var-last P_flat via the
    ! fhat two-pass kernel (Mod_Convective_W). Keep w host-managed for the slice:
    ! pack the legacy var-first P -> w on the host, push H2D once, then each block's
    ! Fluxes_blk(_acc) does zero+SD+DIFFUSIVE only and hands off to Convective_Flux
    ! for the convective contribution (added into the same var-first Res). The rest of
    ! the pipeline (BC/Newstate/diffusive) is unchanged on the P-path. Flag OFF = no-op.
    if (.not. w_setup_done) then
      call Convective_assert_scheme()   ! integer-flag dispatch guard (Min-Mod + HLLC)
      w_setup_done = .true.
    end if
    if (.not. field_is_allocated()) call field_alloc(domain)
    call convective_alloc_fhat(domain)
    call diffusive_alloc_fhat(domain)
    ! Res var-last migration (complete): ALL device
    ! writers accumulate coalesced into the resident var-LAST R_flat; Newstate reads R_flat only.
    ! The var-first Res is a bridge for HOST sources (RotatingFrame/soot) — reset + read
    ! only when one is active (add_r), and rebuilt from R_flat at the IO/diagnostic boundary
    ! (add_wr_to_R_host in Write_Diagnostic).
    ! The per-substage whole-field `R_flat = 0d0` is DEAD — the fused
    ! flux pass-2 (Flux_Pass2_Fused) WRITES (`=`) every interior cell of every local block
    ! each substage, and BC/SST accumulate after it. One-time zero
    ! below keeps the never-written ghost rows deterministic (they are never read).
    if (.not. wr_init_done) then
      !$acc kernels present(R_flat)
      R_flat = 0d0
      !$acc end kernels
      wr_init_done = .true.
    end if
    ! w is kept device-resident + current per substage in Explicit_Step
    ! (BEFORE Refresh_Thermo_Cache and this Fluxes call), so there is NO per-kernel pack
    ! here -- the flux kernels read the resident P_flat directly (the old per-call
    ! update host(P)+pack+H2D is removed). HOW w is made current depends on rank:
    !   n=1  : w interior from Compute_dt's step-start pack + Newstate writes; w ghosts
    !          from Fill_Ghost_Cell_local (the halo runs on w). NO field_sync.
    !   n>1  : field_sync_P_to_flat packs P -> P_flat (the multi-rank halo still fills P ghosts;
    !          halo on P_flat for the MPI path is future work).

    do b = 1, domain % nb ! Loop over blocks
      if (.not. is_local_block(b)) cycle

      call Fluxes_blk ( domain % blk(b) % P,     &
                        domain % blk(b) % r,     &
                        domain % blk(b) % fm1,   &
                        domain % blk(b) % fm2,   &
                        domain % blk(b) % fm3,   &
                        domain % blk(b) % dl_f,  &
                        domain % blk(b) % beta,  &
                        domain % blk(b) % yn,    &
                        domain % blk(b) % m_f,   &
                        domain % blk(b) % rho_c, &
                        domain % blk(b) % T_c,   &
                        domain % blk(b) % dim,   &
                        b )
    end do

  end subroutine Fluxes



  !----------------------------------------------------------------------------
  ! GPU twin of Fluxes_blk. Same algebra / same flux callees (Convective_Flux,
  ! Diffusive_Flux are `acc routine seq`), but each direction is an
  ! `acc parallel loop collapse(2)` over the two cross-directions with the
  ! flux-direction loop kept `seq`. That preserves the host accumulation order
  ! exactly (Res is updated +/- on two neighbours), so no atomics are needed and
  ! the result matches the CPU to FP-reorder tolerance. Whole-block H2D/D2H per
  ! call (correctness-first; residency is a later step).
  !----------------------------------------------------------------------------
  subroutine Fluxes_blk ( Prim, Res, fm1, fm2, fm3, dl, beta, yn, M, rho_c, T_c, &
                              n, b)
    use MOSE_Mod_Convective, only: fhat, fhat_j, fhat_k
    use MOSE_Mod_Diffusive,  only: dfhat, dfhat_j, dfhat_k
    use MOSE_Global_m, only: gc, nprim, VISC, NPRIM_DEV, ROTATING_FRAME, soot
    use MOSE_Lib_Shock_Detector
    use MOSE_Field_m, only: P_flat, field_block_index_params
    use MOSE_Global_m, only: np
    ! Convective flux: always the CUF kernel (production, unconditional).
    use MOSE_Mod_Convective, only: Convective_Flux
    ! Diffusive flux: the CUF kernel reading the resident thermo cache (production,
    ! unconditional).
    use MOSE_Mod_Diffusive, only: Diffusive_Flux

    implicit none
    integer, intent(in) :: b   ! 1-based global block index (w slice, GPU device path)
    integer, dimension(3), intent(in) :: n
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(out) :: Res
    real(R8), dimension(4, 0:n(1), n(2), n(3)), intent(in) :: fm1
    real(R8), dimension(4, n(1), 0:n(2), n(3)), intent(in) :: fm2
    real(R8), dimension(4, n(1), n(2), 0:n(3)), intent(in) :: fm3
    real(R8), dimension(3, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: dl
    real(R8), dimension(1:n(1), 1:n(2), 1:n(3)), intent(inout) :: beta
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: yn
    real(R8), dimension(3, 3, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: M
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: rho_c, T_c
    ! Local
    integer :: i, j, k, s
    integer :: boff_sd, pi_sd, pj_sd, c0_sd, sdid_   ! the MOSE_open shock sensor on w (sd_dev > 0)
    ! Zero the var-first Res only when a host bridge source
    ! (RotatingFrame/soot) will accumulate into it; see the reset kernel below.
    logical :: add_r
    real(R8) :: rho_n10(10), T_n10(10)
    ! Contiguous gather of the convective 4-cell Prim stencil + 2-cell Res, so the
    ! strided j/k face slices Prim(:,i,j-1:j+2,k) are NOT copy-in/out'd per call
    ! (Convective_Flux takes NPRIM_DEV-leading dummies). Without this the j/k flux
    ! kernels cost ~4x the i kernel.
    real(R8) :: ploc(NPRIM_DEV,-1:2), rloc(NPRIM_DEV,0:1)
    ! Contiguous gather of the strided dl(d,·) 4-cell stencil. Passing the
    ! strided section dl(1,i-1:i+2,..) straight into the `acc routine seq`
    ! Convective_Flux (explicit-shape dl(-1:2)) makes nvfortran misread it on the
    ! device (the host compiler silently gathers; device codegen does not) -> dl
    ! garbage -> reconstruction blows up. Gather into a contiguous local first.
    real(R8) :: dlc(-1:2)
    ! Dims as scalars -- a dummy array referenced inside the acc region
    ! is implicitly copy()d per launch (upload+download+sync); scalars ride the
    ! launch parameters for free.
    integer :: nb1_, nb2_, nb3_
    nb1_ = n(1); nb2_ = n(2); nb3_ = n(3)

    ! copy(Res): the host Res ghost cells are 0 and the reset kernel zeros only
    ! the INTERIOR, so copyout would leave the device Res ghosts uninitialized.
    ! copyout(beta): beta is NOT internal — BC_Fluxes (Lib_BC) reads blk%beta for
    ! the boundary-face reconstruction, so the device-set value (1, or the SD
    ! sensor) must be pushed back to the host or BC sees the stale 0 and diverges.
    !$acc data copyin(Prim, fm1, fm2, fm3, dl, yn, M, rho_c, T_c) &
    !$acc&     copy(Res) copyout(beta)

    ! Reset residuals + shock sensor
    ! On the W path all device writers accumulate into R_flat (zeroed
    ! per substage in Fluxes), so the var-first Res is a pure bridge for host sources
    ! (RotatingFrame/soot). Zero it only when such a source is active; otherwise the
    ! per-substage full-field var-first write (the "reset" share of this kernel)
    ! disappears and only the beta init remains.
    add_r = ROTATING_FRAME .or. soot
    if (add_r) then
      !$acc parallel loop collapse(3) present(Res, beta)
      do k = 1, nb3_
      do j = 1, nb2_
      do i = 1, nb1_
        Res(:,i,j,k) = 0d0
        beta(i,j,k)  = 1d0
      end do ; end do ; end do
    else
      !$acc parallel loop collapse(3) present(beta)
      do k = 1, nb3_
      do j = 1, nb2_
      do i = 1, nb1_
        beta(i,j,k)  = 1d0
      end do ; end do ; end do
    end if

    ! Shock-detector
    if (sd_dev > 0) then
      ! MOSE_open (Mod_Fluxes Fluxes_blk): Tramel or Chen from the pressure of the cell and of its six face neighbours,
      ! on the resident P_flat (ghosts filled before the fluxes); beta of the interior cells, the face takes its left cell's.
      call field_block_index_params(b, boff_sd, pi_sd, pj_sd)
      sdid_ = sd_dev
      !$acc parallel loop collapse(3) present(P_flat, beta) private(c0_sd)
      do k = 1, nb3_
      do j = 1, nb2_
      do i = 1, nb1_
        c0_sd = boff_sd + (i+gc-1) + (j+gc-1)*pi_sd + (k+gc-1)*pi_sd*pj_sd + 1
        beta(i,j,k) = SD_beta ( sdid_, P_flat(c0_sd,np), P_flat(c0_sd-1,np), P_flat(c0_sd+1,np), P_flat(c0_sd-pi_sd,np), P_flat(c0_sd+pi_sd,np), &
                                    P_flat(c0_sd-pi_sd*pj_sd,np), P_flat(c0_sd+pi_sd*pj_sd,np) )
      end do ; end do ; end do
    else if (SD) then
      !$acc parallel loop collapse(3) present(Prim, beta)
      do k = 1, nb3_
      do j = 1, nb2_
      do i = 1, nb1_
        beta(i,j,k) = SD_Tramel ( Prim(np,i-1:i+1,j-1:j+1,k-1:k+1) )
      end do ; end do ; end do
    end if


    ! The CONVECTIVE contribution comes from the var-last fhat
    ! two-pass kernel reading the resident P_flat (packed from P by the caller). It ADDS
    ! into the var-first Res just zeroed above; the diffusive contribution (also
    ! added above) and this convective contribution share the one Res. fm/dl/beta/Res
    ! are present in this data region; w/fhat are enter-data resident.
    call Convective_Flux ( b, Res, fm1, fm2, fm3, dl, beta, n )

    ! The DIFFUSIVE contribution comes from the var-last dfhat
    ! two-pass kernel reading the same resident P_flat + the cell caches/metric (present in
    ! this data region). It ADDS into the same var-first Res with the DIFFUSIVE sign
    ! (Res(left) -= Flux, Res(right) += Flux), so convective + diffusive share one Res.
    ! VISC-gated inside the kernel. M (== m_f), yn, rho_c, T_c are present (copyin above).
    call Diffusive_Flux ( b, Res, fm1, fm2, fm3, M, yn, rho_c, T_c, n )

    ! ONE fused pass-2 for the six directional scatters. Same sync queue as the
    ! six pass-1 launches above -> stream-ordered after all of them.
    call Flux_Pass2_Fused ( b, n )

    !$acc end data

  end subroutine Fluxes_blk


  !----------------------------------------------------------------------------
  ! The six
  ! per-direction pass-2 scatters (conv i/j/k + diff i/j/k) fused into ONE kernel.
  ! Per cell, the six +/- face contributions accumulate in a REGISTER in the exact
  ! HISTORICAL order of the six deleted sequential kernels (conv i +,-; j +,-; k +,-;
  ! then, when VISC, diff i -,+; j -,+; k -,+), and R_flat is written ONCE (`=`, from the
  ! same 0 the old whole-field zero provided). Identical per-component FP chain ->
  ! bit-exact; R_flat traffic drops from 1 write + 6 read-modify-writes to 1 write per
  ! cell per substage. Face indexing/edge guards copied verbatim from the deleted
  ! loops. fhat*/dfhat*/R_flat are enter-data resident (present).
  !----------------------------------------------------------------------------

  !> The second pass of the flux phase: the six directional face fluxes of the block (fhat, fhat_j, fhat_k and the
  !> diffusive twins dfhat, written by the pass-1 kernels) scattered into the residual R_flat of every cell, in the
  !> order of pass2_body.inc (the order of the twelve additions fixes the result to the last bit).
  subroutine Flux_Pass2_Fused ( b, n )
    use MOSE_Global_m, only: gc, nprim, VISC, NPRIM_DEV
    use MOSE_Field_m,  only: R_flat, field_block_index_params
    use MOSE_Mod_Convective, only: fhat, fhat_j, fhat_k
    use MOSE_Mod_Diffusive,  only: dfhat, dfhat_j, dfhat_k
    implicit none
    integer, intent(in) :: b, n(3)
    ! Local
    integer  :: i, j, k, s, f, c0, boff, pi, pj, n1, n2, n3
    logical  :: visc_l
    real(R8) :: tmp(NPRIM_DEV)

    call field_block_index_params(b, boff, pi, pj)
    n1 = n(1); n2 = n(2); n3 = n(3)
    visc_l = VISC

    !$acc parallel loop collapse(3) gang vector &
    !$acc&  present(R_flat, fhat, fhat_j, fhat_k, dfhat, dfhat_j, dfhat_k) &
    !$acc&  private(s, f, c0, tmp)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
#include "pass2_body.inc"
    end do; end do; end do

  end subroutine Flux_Pass2_Fused


end module MOSE_Mod_Fluxes
