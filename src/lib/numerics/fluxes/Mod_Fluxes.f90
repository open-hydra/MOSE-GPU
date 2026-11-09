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
    use MOSE_Lib_Multirate_Base, only: mr_block_active, mr_list_on, mr_stage_on
    use MOSE_Pack_m, only: pack_setup, pack_on

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

    ! multirate stage of a list level (mr-list-level): the flux phase on the level's lists, one launch per kernel over
    ! every local block, the same face and cell values as the box sweep below
    if ( mr_list_on() ) then
      call Fluxes_List ()
      return
    end if
    ! the explicit step on the block pack (MOSE_Pack_m): one launch per kernel type over each batch of blocks
    ! (MOSE_FLUX_PACK=0, a shock detector or a host source keep the per-block sweep below)
    if ( .not. mr_stage_on() ) then
      call pack_setup ( domain )
      if ( pack_on() ) then
        call Fluxes_Pack ()
        return
      end if
    end if

    do b = 1, domain % nb ! Loop over blocks
      if (.not. is_local_block(b)) cycle
      if (.not. mr_block_active(b)) cycle   ! multirate stage: no cell of the stage's level in this block

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


  !> The flux phase of a multirate stage on the lists of its level (Lib_Multirate_Base mr_fl_type / mr_cl_type, built by
  !> Mod_Multirate_Partition mr_build_lists): the convective and diffusive pass 1 per direction, the telescoping on the
  !> interface faces, pass 2 on the level's cells. Thirteen launches for the whole level instead of ~14 per (block, box).
  subroutine Fluxes_List ()
    use MOSE_Global_m, only: VISC
    use MOSE_Mod_Convective, only: Convective_Flux_List
    use MOSE_Mod_Diffusive,  only: Diffusive_Flux_List
    use MOSE_Lib_Multirate_Base, only: obj_mr, mrfl, mrcl, MR_Iface_Update_List
    implicit none
    integer :: l
    l = obj_mr%lcur
    call Convective_Flux_List ( mrfl(:,l) )
    call Diffusive_Flux_List ( mrfl(:,l) )
    call MR_Iface_Update_List ( mrfl(:,l), VISC )
    if ( mrcl(l)%n > 0 ) &
      call Flux_Pass2_List ( mrcl(l)%c, mrcl(l)%pf, mrcl(l)%n, mrfl(1,l)%fh, mrfl(2,l)%fh, mrfl(3,l)%fh, &
                             mrfl(1,l)%dfh, mrfl(2,l)%dfh, mrfl(3,l)%dfh, &
                             size(mrfl(1,l)%fh,1), size(mrfl(2,l)%fh,1), size(mrfl(3,l)%fh,1) )
  end subroutine Fluxes_List

  !> Pass 2 on the cell list of a level: per cell the six faces of the lists (positions pf: i+, i-, j+, j-, k+, k-; 0 = a
  !> boundary face) summed into R_flat in the order of pass2_body.inc -- convective i +,-, j +,-, k +,-, then, when VISC,
  !> diffusive i -,+, j -,+, k -,+ -- the per-component chain of the box sweep, so the same R_flat.
  subroutine Flux_Pass2_List ( cc, pf, ne, fh1, fh2, fh3, dfh1, dfh2, dfh3, n1f, n2f, n3f )
    use MOSE_Global_m, only: nprim, VISC, NPRIM_DEV
    use MOSE_Lib_Multirate_Base, only: mr_q
    use MOSE_Field_m,  only: R_flat
    implicit none
    integer,  intent(in) :: ne, n1f, n2f, n3f, cc(ne), pf(6, ne)
    real(R8), intent(in) :: fh1(n1f, nprim), fh2(n2f, nprim), fh3(n3f, nprim)
    real(R8), intent(in) :: dfh1(n1f, nprim), dfh2(n2f, nprim), dfh3(n3f, nprim)
    integer  :: e, s, q, c0
    logical  :: visc_l
    real(R8) :: tmp(NPRIM_DEV)
    visc_l = VISC
    !$acc parallel loop gang vector present(R_flat, cc, pf, fh1, fh2, fh3, dfh1, dfh2, dfh3) private(s, q, c0, tmp) async(mr_q)
    do e = 1, ne
      c0 = cc(e)
      do s = 1, nprim; tmp(s) = 0d0; end do
      q = pf(1,e); if (q > 0) then; do s = 1, nprim; tmp(s) = tmp(s) + fh1(q, s); end do; end if
      q = pf(2,e); if (q > 0) then; do s = 1, nprim; tmp(s) = tmp(s) - fh1(q, s); end do; end if
      q = pf(3,e); if (q > 0) then; do s = 1, nprim; tmp(s) = tmp(s) + fh2(q, s); end do; end if
      q = pf(4,e); if (q > 0) then; do s = 1, nprim; tmp(s) = tmp(s) - fh2(q, s); end do; end if
      q = pf(5,e); if (q > 0) then; do s = 1, nprim; tmp(s) = tmp(s) + fh3(q, s); end do; end if
      q = pf(6,e); if (q > 0) then; do s = 1, nprim; tmp(s) = tmp(s) - fh3(q, s); end do; end if
      if (visc_l) then
        q = pf(1,e); if (q > 0) then; do s = 1, nprim; tmp(s) = tmp(s) - dfh1(q, s); end do; end if
        q = pf(2,e); if (q > 0) then; do s = 1, nprim; tmp(s) = tmp(s) + dfh1(q, s); end do; end if
        q = pf(3,e); if (q > 0) then; do s = 1, nprim; tmp(s) = tmp(s) - dfh2(q, s); end do; end if
        q = pf(4,e); if (q > 0) then; do s = 1, nprim; tmp(s) = tmp(s) + dfh2(q, s); end do; end if
        q = pf(5,e); if (q > 0) then; do s = 1, nprim; tmp(s) = tmp(s) - dfh3(q, s); end do; end if
        q = pf(6,e); if (q > 0) then; do s = 1, nprim; tmp(s) = tmp(s) + dfh3(q, s); end do; end if
      end if
      do s = 1, nprim; R_flat(c0,s) = tmp(s); end do
    end do
  end subroutine Flux_Pass2_List



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
    use MOSE_Lib_Multirate_Base, only: mr_nbox, mr_get_box, MR_Iface_Update

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
    ! The cells this call works on: the whole block, or in a multirate stage the boxes of the stage's level (Lib_Multirate_Base)
    integer :: ib, nbx, lo(3), hi(3), bl1_, bh1_, bl2_, bh2_, bl3_, bh3_
    nb1_ = n(1); nb2_ = n(2); nb3_ = n(3)
    nbx = mr_nbox(b)

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
    do ib = 1, nbx
    call mr_get_box ( b, ib, n, lo, hi )
    lo = max( lo - 1, 1 )   ! the face kernels read beta of the left cell of a face: one cell before the box
    bl1_ = lo(1); bh1_ = hi(1); bl2_ = lo(2); bh2_ = hi(2); bl3_ = lo(3); bh3_ = hi(3)
    if (add_r) then
      !$acc parallel loop collapse(3) present(Res, beta)
      do k = bl3_, bh3_
      do j = bl2_, bh2_
      do i = bl1_, bh1_
        Res(:,i,j,k) = 0d0
        beta(i,j,k)  = 1d0
      end do ; end do ; end do
    else
      !$acc parallel loop collapse(3) present(beta)
      do k = bl3_, bh3_
      do j = bl2_, bh2_
      do i = bl1_, bh1_
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
      do k = bl3_, bh3_
      do j = bl2_, bh2_
      do i = bl1_, bh1_
        c0_sd = boff_sd + (i+gc-1) + (j+gc-1)*pi_sd + (k+gc-1)*pi_sd*pj_sd + 1
        beta(i,j,k) = SD_beta ( sdid_, P_flat(c0_sd,np), P_flat(c0_sd-1,np), P_flat(c0_sd+1,np), P_flat(c0_sd-pi_sd,np), P_flat(c0_sd+pi_sd,np), &
                                    P_flat(c0_sd-pi_sd*pj_sd,np), P_flat(c0_sd+pi_sd*pj_sd,np) )
      end do ; end do ; end do
    else if (SD) then
      !$acc parallel loop collapse(3) present(Prim, beta)
      do k = bl3_, bh3_
      do j = bl2_, bh2_
      do i = bl1_, bh1_
        beta(i,j,k) = SD_Tramel ( Prim(np,i-1:i+1,j-1:j+1,k-1:k+1) )
      end do ; end do ; end do
    end if
    end do


    ! The CONVECTIVE contribution comes from the var-last fhat
    ! two-pass kernel reading the resident P_flat (packed from P by the caller). It ADDS
    ! into the var-first Res just zeroed above; the diffusive contribution (also
    ! added above) and this convective contribution share the one Res. fm/dl/beta/Res
    ! are present in this data region; w/fhat are enter-data resident.
    do ib = 1, nbx
      call mr_get_box ( b, ib, n, lo, hi )
      call Convective_Flux ( b, Res, fm1, fm2, fm3, dl, beta, n, lo, hi )
    end do

    ! The DIFFUSIVE contribution comes from the var-last dfhat
    ! two-pass kernel reading the same resident P_flat + the cell caches/metric (present in
    ! this data region). It ADDS into the same var-first Res with the DIFFUSIVE sign
    ! (Res(left) -= Flux, Res(right) += Flux), so convective + diffusive share one Res.
    ! VISC-gated inside the kernel. M (== m_f), yn, rho_c, T_c are present (copyin above).
    do ib = 1, nbx
      call mr_get_box ( b, ib, n, lo, hi )
      call Diffusive_Flux ( b, Res, fm1, fm2, fm3, M, yn, rho_c, T_c, n, lo, hi )
    end do

    ! Multirate: the flux telescoping at the level interfaces of the block, on the face buffers, before pass 2 reads them
    ! (no-op outside a multirate stage)
    if (VISC) then
      call MR_Iface_Update ( b, fhat, fhat_j, fhat_k, dfhat, dfhat_j, dfhat_k, size(fhat,1), .true. )
    else
      call MR_Iface_Update ( b, fhat, fhat_j, fhat_k, fhat, fhat_j, fhat_k, size(fhat,1), .false. )
    end if

    ! ONE fused pass-2 for the six directional scatters. Same sync queue as the
    ! six pass-1 launches above -> stream-ordered after all of them.
    do ib = 1, nbx
      call mr_get_box ( b, ib, n, lo, hi )
      call Flux_Pass2_Fused ( b, n, lo, hi )
    end do

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
  subroutine Flux_Pass2_Fused ( b, n, lo, hi )
    use MOSE_Global_m, only: gc, nprim, VISC, NPRIM_DEV
    use MOSE_Field_m,  only: R_flat, field_block_index_params
    use MOSE_Mod_Convective, only: fhat, fhat_j, fhat_k
    use MOSE_Mod_Diffusive,  only: dfhat, dfhat_j, dfhat_k
    implicit none
    integer, intent(in) :: b, n(3)
    integer, intent(in) :: lo(3), hi(3)   ! the cells: the whole block, or a box of a multirate stage
    ! Local
    integer  :: i, j, k, s, f, c0, boff, pi, pj, n1, n2, n3, lo1, hi1, lo2, hi2, lo3, hi3
    logical  :: visc_l
    real(R8) :: tmp(NPRIM_DEV)

    call field_block_index_params(b, boff, pi, pj)
    n1 = n(1); n2 = n(2); n3 = n(3)
    lo1 = lo(1); hi1 = hi(1); lo2 = lo(2); hi2 = hi(2); lo3 = lo(3); hi3 = hi(3)
    visc_l = VISC

    !$acc parallel loop collapse(3) gang vector &
    !$acc&  present(R_flat, fhat, fhat_j, fhat_k, dfhat, dfhat_j, dfhat_k) &
    !$acc&  private(s, f, c0, tmp)
    do k = lo3, hi3
    do j = lo2, hi2
    do i = lo1, hi1
#include "pass2_body.inc"
    end do; end do; end do

  end subroutine Flux_Pass2_Fused


  !> Pass 2 of batch ib of the block pack (MOSE_Pack_m): the six face buffers of the batch summed into R_flat per cell,
  !> in the order of pass2_body.inc (convective i +,-, j +,-, k +,-, then, when VISC, diffusive i -,+, j -,+, k -,+: the
  !> per-component chain of the per-block sweep, so the same R_flat), one gang of 128 per CUDA-block entry of the batch's
  !> cell table; the face rows carry the block's offsets bt_ofh(d,p).
  subroutine Flux_Pass2_Pack ( ib )
    use MOSE_Global_m, only: gc, nprim, VISC, NPRIM_DEV
    use MOSE_Field_m,  only: R_flat
    use MOSE_Pack_m
    implicit none
    integer, intent(in) :: ib
    integer  :: e, t, t_, p, n1, n2, n3, boff, pi, pj, i, j, k, rem_, s, f, c0, cb0, ncb, o1, o2, o3
    logical  :: visc_l
    real(R8) :: tmp(NPRIM_DEV)
    visc_l = VISC
    cb0 = cb_n(0, 4, ib); ncb = cb_n(1, 4, ib)
    if ( ncb == 0 ) return
    !$acc parallel loop gang vector_length(128) &
    !$acc&  present(R_flat, fhp1, fhp2, fhp3, dfhp1, dfhp2, dfhp3, cb_p, cb_t0, bt_n1, bt_n2, bt_n3, bt_boff, bt_pi, bt_pj, bt_ofh) &
    !$acc&  private(p, n1, n2, n3, boff, pi, pj, o1, o2, o3)
    do e = 1, ncb
      p = cb_p(cb0+e)
      n1 = bt_n1(p); n2 = bt_n2(p); n3 = bt_n3(p); boff = bt_boff(p); pi = bt_pi(p); pj = bt_pj(p)
      o1 = bt_ofh(1,p); o2 = bt_ofh(2,p); o3 = bt_ofh(3,p)
      !$acc loop vector private(tmp, t_, s, f, c0, i, j, k, rem_)
      do t = 1, PACK_NT
        t_ = cb_t0(cb0+e) + t
        if ( t_ <= n1*n2*n3 ) then
          i    = mod(t_-1, n1) + 1
          rem_ = (t_-1) / n1
          j    = mod(rem_, n2) + 1
          k    = rem_ / n2 + 1
          c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
          do s = 1, nprim; tmp(s) = 0d0; end do
          ! convective i (+ right face, - left face)
          if (i < n1) then
            f = o1 + (i-1) + (j-1)*(n1-1) + (k-1)*(n1-1)*n2 + 1
            do s = 1, nprim; tmp(s) = tmp(s) + fhp1(f, s); end do
          end if
          if (i > 1) then
            f = o1 + (i-2) + (j-1)*(n1-1) + (k-1)*(n1-1)*n2 + 1
            do s = 1, nprim; tmp(s) = tmp(s) - fhp1(f, s); end do
          end if
          ! convective j
          if (j < n2) then
            f = o2 + (i-1) + (j-1)*n1 + (k-1)*n1*(n2-1) + 1
            do s = 1, nprim; tmp(s) = tmp(s) + fhp2(f, s); end do
          end if
          if (j > 1) then
            f = o2 + (i-1) + (j-2)*n1 + (k-1)*n1*(n2-1) + 1
            do s = 1, nprim; tmp(s) = tmp(s) - fhp2(f, s); end do
          end if
          ! convective k
          if (k < n3) then
            f = o3 + (i-1) + (j-1)*n1 + (k-1)*n1*n2 + 1
            do s = 1, nprim; tmp(s) = tmp(s) + fhp3(f, s); end do
          end if
          if (k > 1) then
            f = o3 + (i-1) + (j-1)*n1 + (k-2)*n1*n2 + 1
            do s = 1, nprim; tmp(s) = tmp(s) - fhp3(f, s); end do
          end if
          if (visc_l) then
            ! diffusive i (- right face, + left face)
            if (i < n1) then
              f = o1 + (i-1) + (j-1)*(n1-1) + (k-1)*(n1-1)*n2 + 1
              do s = 1, nprim; tmp(s) = tmp(s) - dfhp1(f, s); end do
            end if
            if (i > 1) then
              f = o1 + (i-2) + (j-1)*(n1-1) + (k-1)*(n1-1)*n2 + 1
              do s = 1, nprim; tmp(s) = tmp(s) + dfhp1(f, s); end do
            end if
            ! diffusive j
            if (j < n2) then
              f = o2 + (i-1) + (j-1)*n1 + (k-1)*n1*(n2-1) + 1
              do s = 1, nprim; tmp(s) = tmp(s) - dfhp2(f, s); end do
            end if
            if (j > 1) then
              f = o2 + (i-1) + (j-2)*n1 + (k-1)*n1*(n2-1) + 1
              do s = 1, nprim; tmp(s) = tmp(s) + dfhp2(f, s); end do
            end if
            ! diffusive k
            if (k < n3) then
              f = o3 + (i-1) + (j-1)*n1 + (k-1)*n1*n2 + 1
              do s = 1, nprim; tmp(s) = tmp(s) - dfhp3(f, s); end do
            end if
            if (k > 1) then
              f = o3 + (i-1) + (j-1)*n1 + (k-2)*n1*n2 + 1
              do s = 1, nprim; tmp(s) = tmp(s) + dfhp3(f, s); end do
            end if
          end if
          do s = 1, nprim; R_flat(c0,s) = tmp(s); end do
        end if
      end do
    end do
  end subroutine Flux_Pass2_Pack

  !> The flux phase of the explicit step on the block pack (MOSE_Pack_m): per batch the convective pass 1 (3 launches),
  !> the diffusive pass 1 (3 face-thermo fills + 3 launches, with viscosity) and pass 2 (1 launch); no per-block beta
  !> init (beta = 1, set once by pack_setup), no shock detector, no host source, no multirate interface (the pack runs
  !> only outside the stages).
  subroutine Fluxes_Pack ()
    use MOSE_Mod_Convective, only: Convective_Flux_Pack
    use MOSE_Mod_Diffusive,  only: Diffusive_Flux_Pack
    use MOSE_Pack_m, only: pack_nbat
    implicit none
    integer :: ib
    do ib = 1, pack_nbat()
      call Convective_Flux_Pack ( ib )
      call Diffusive_Flux_Pack ( ib )
      call Flux_Pass2_Pack ( ib )
    end do
  end subroutine Fluxes_Pack


end module MOSE_Mod_Fluxes
