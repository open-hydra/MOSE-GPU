module MOSE_Lib_Extrapolation
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Types_m
  use MOSE_Parameters_m
  use MOSE_Global_m

  implicit none

! Halo-on-w. The device ghost-fill twins write the var-LAST flat SoA
! w ghosts (read w interior) instead of var-first P. Flat row == cidx; boff/pi/pj are
! host scalars (field_block_index_params), passed as ARGS. WGCELL is the single-block
! index (ghost kernels); WMCELL/WSCELL are the receiver/donor indices for the cross-block
! connection kernel (two different block offsets in the one P_flat). Lvalue-safe (the
! ghost kernels WRITE w(row,s)). Gated behind _OPENACC; the CPU build never
! references these (the W kernels are entirely #if'd out).
#  define WGCELL(s,i,j,k) w( boff   + ((i)+gc-1) + ((j)+gc-1)*pi   + ((k)+gc-1)*pi*pj     + 1, s )
#  define WMCELL(s,i,j,k) w( boff_m + ((i)+gc-1) + ((j)+gc-1)*pi_m + ((k)+gc-1)*pi_m*pj_m + 1, s )
#  define WSCELL(s,i,j,k) w( boff_s + ((i)+gc-1) + ((j)+gc-1)*pi_s + ((k)+gc-1)*pi_s*pj_s + 1, s )

  !> Wall ghost fill mode for the W halo (0 cubic, 1 mirror, 2 zero-gradient): the key wall-ghost = cubic | auto | mirror | zg of
  !> [MOSE-scheme], read by IO_Ini into wall_ghost_key and resolved once at the first ghost fill. cubic (default) is what MOSE_open
  !> does (Lib_Ghost case default, Ghost_Extrapolate); auto = mirror when the case has a viscous wall record (MOSE_open 301-309 =
  !> internal types 5 6 8 9 10 12 13), cubic otherwise -- the GPU default before 2026-10-16 (doc/THOR_BC_WALL_IS_SYMMETRY.md).
  character(len=16), save :: wall_ghost_key = 'cubic'
  integer, save :: wall_ghost_mode = 0
  logical, save :: wgm_checked = .false.
  !> Ghost cells of the extrapolation faces (MOSE_open 400 = internal type 11): the key extrapolation-ghost = cubic | zero-gradient
  !> of [MOSE-scheme], read by IO_Ini. cubic (default) is what MOSE_open does (Lib_Ghost: 400 is not in its `case(0,401:407,420)`
  !> zero-gradient list and falls in `case default`, Ghost_Extrapolate); zero-gradient is the legacy GPU fill (default before
  !> 2026-10-16).
  logical, save :: extrap_ghost_cubic = .true.
  !> Wall ghost of the turbulence variables for the RANS source terms (Extrapolate_Wall): the key rans-wall-ghost =
  !> copy | wall-value of [MOSE-scheme], read by IO_Ini. copy (default) is what MOSE_open does (Lib_Ghost
  !> Ghost_Wall_Extrapolation copies P(np:nprim) of the boundary cell, turbulence included, right before RANS_Source_Terms);
  !> wall-value is the former device fill (SA nu-tilde ghost = -interior, SST k ghost = -interior and omega ghost =
  !> 2 omega_wall - interior; default before 2026-10-16).
  logical, save :: rans_wall_copy = .true.

contains

  !> The ghost cells outside the time loop (setup): P is the state. They are filled on the device, on P_flat, by the
  !> kernels the time loop uses -- the flat-w fill on one rank, the multi-rank driver (inter-rank exchange + local fill)
  !> on more -- and written back into P, ghost cells included.
  !> The ghost cells of every local block on the resident P_flat (MOSE_open's entry): the local fill on one rank, the
  !> multi-rank driver (inter-rank exchange + local fill) on more. Called by every thread of the time loop's parallel
  !> region, and by the setup (Wrap_Setup packs P into w before and unpacks after).
  subroutine Fill_Ghost_Cell ( domain )
    use MOSE_Mod_MPI, only: mpi_size_
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    if ( mpi_size_ == 1 ) then
      !$omp single
      call Fill_Ghost_Cell_local ( domain )
      !$omp end single
    else
      call Fill_Ghost_Cell_mr ( domain )
    end if
  end subroutine Fill_Ghost_Cell













  !> Device twin of Extrapolate_Wall: per-block single-layer wall-ghost reflection on the
  !> resident P (velocity inverted, species/pressure/RANS copied) for the turbulence RHS.
  !> Local (own interior -> own ghost), so it ports without MPI; disjoint ghost targets per
  !> entry -> race-free; sorted bc order is immaterial (same reason as Fill_Ghost_Cell_acc).
  !>
  !> Reflects the velocity (no-slip) into the wall ghost on resident P, and for SST and Wilcox2006 also sets
  !> the k/omega ghost (SST_ / Wilcox2006_Extrapolate_Wall inlined, the same lines; needs the boundary-cell laminar viscosity,
  !> so it computes rho/T/mil there). Replaces the host Extrapolate_Wall + the omega ghost that
  !> host BC_Wall_Heat writes (device BC kernels are P-read-only -> the omega ghost lives here).
  !> Reached only on the RANS device path (Mod_Explicit gates it by `if (RANS)`); validated by
  !> the Flat_Plate_SST CPU-vs-device A/B.
  subroutine Extrapolate_Wall ( domain )
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Field_m, only: field_block_index_params
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: b
    integer :: boff, pim, pjm
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      if (domain % bc_blk_hi(b) < domain % bc_blk_lo(b)) cycle
      ! var-last build: the consumers (SST_Blk + diffusive RANS flux) read the resident P_flat,
      ! NOT P. Every other BC ghost is written to w (Fill_Ghost_Cell_local); the wall ghost must
      ! be too, else the near-wall SST source reads a stale w ghost (writing P is invisible to it).
      call field_block_index_params(b, boff, pim, pjm)
      call Extrapolate_Wall_blk ( boff, pim, pjm, domain % blk(b) % yn, domain % blk(b) % mil_c, &
                                        domain % blk(b) % dim, domain % bc_bi, domain % bc_bj, domain % bc_bk, &
                                        domain % bc_bf, domain % bc_typ, domain % bc_tw, domain % nbound, &
                                        domain % bc_blk_lo(b), domain % bc_blk_hi(b) )
      if ( rans_wall_copy ) &   ! MOSE_open: the turbulence ghost = the boundary cell (a separate launch: the kernel above is untouched)
      call Extrapolate_Wall_RANS_Copy_blk ( boff, pim, pjm, domain % bc_bi, domain % bc_bj, domain % bc_bk, &
                                              domain % bc_bf, domain % bc_typ, domain % nbound, &
                                              domain % bc_blk_lo(b), domain % bc_blk_hi(b) )
    end do
  end subroutine Extrapolate_Wall

  !> Writes the wall ghost into the resident P_flat
  !> (consumers read w). velocity inverted (no-slip), species/p copied, SST k/omega ghost
  !> (k_wall=0, omega ghost = 2*omega_wall - interior; omega_wall = omega_wall_c(1)*mil_c/d^2 from the
  !> resident cache = the boundary-cell mu_l, matching the CPU per-cell recompute). Reads/writes
  !> via the cidx flat row (boff/pim/pjm = field_block_index_params), mirroring WGCELL.
  subroutine Extrapolate_Wall_blk ( boff, pim, pjm, yn, mil_c, n, c_bi, c_bj, c_bk, c_bf, c_typ, c_tw, nbnd, lo, hi )
    use MOSE_Parameters_m, only: guide
    use MOSE_Lib_RANS,     only: eddy_dev
    use MOSE_Field_m,      only: P_flat, mil_c_f
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, co_k_mi_lam_Wilke_dev
    implicit none
    integer,  intent(in) :: boff, pim, pjm, n(3), nbnd, lo, hi
    real(R8), intent(in) :: yn   (1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc)
    real(R8), intent(in) :: mil_c(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc)
    integer,  intent(in) :: c_bi(nbnd), c_bj(nbnd), c_bk(nbnd), c_bf(nbnd), c_typ(nbnd)
    real(R8), intent(in) :: c_tw(nbnd)
    integer  :: idx, Im, Jm, Km, Fm, Ig, Jg, Kg, ci, cg, s, Tint(2)
    real(R8) :: dist_w, om_wall, rho, Rgas, rho_wall, Tdiff, mil_w, kl_w
    real(R8) :: spw(NPRIM_DEV)
    !$acc parallel loop gang vector present(P_flat, mil_c_f) copyin(c_bi, c_bj, c_bk, c_bf, c_typ, c_tw, yn, mil_c) &
    !$acc&   private(Im, Jm, Km, Fm, Ig, Jg, Kg, ci, cg, s, Tint, dist_w, om_wall, &
    !$acc&           rho, Rgas, rho_wall, Tdiff, mil_w, kl_w, spw)
    do idx = lo, hi
      ! Is_Wall(type): walls = {5,6,8,9,10,12,13,1000}
      select case (c_typ(idx))
        case (5, 6, 8, 9, 10, 12, 13, 1000)
          Im = c_bi(idx); Jm = c_bj(idx); Km = c_bk(idx); Fm = c_bf(idx)
          Ig = Im - guide(Fm,1); Jg = Jm - guide(Fm,2); Kg = Km - guide(Fm,3)
          ci = boff + (Im+gc-1) + (Jm+gc-1)*pim + (Km+gc-1)*pim*pjm + 1   ! interior flat row
          cg = boff + (Ig+gc-1) + (Jg+gc-1)*pim + (Kg+gc-1)*pim*pjm + 1   ! ghost flat row
          do s = nu, nw
            P_flat(cg, s) = - P_flat(ci, s)        ! invert velocity (no-slip)
          end do
          do s = 1, nsc
            P_flat(cg, s) =   P_flat(ci, s)        ! copy species
          end do
          P_flat(cg, np) = P_flat(ci, np)          ! copy pressure (zero-gradient)
          if (RANS .and. eddy_dev == 2) then
            P_flat(cg, nt) = - P_flat(ci, nt)                               ! SA nu-tilde ghost (Spalart_Extrapolate_Wall: nu-tilde_wall = 0)
          else if (RANS) then
            dist_w  = yn(Im,Jm,Km)
            if (c_typ(idx) == 6) then
              ! Isothermal wall (host BC_Wall_Temperature + SST_Extrapolate_Wall): the
              ! wall state is (p(cell), T_wall, same mass fractions) -> rho_wall and
              ! mil(T_wall) replace the adiabatic cell values, and the omega ghost keeps
              ! the rho/rho_wall factor: ghost = rho*2*omega_wall/rho_wall - interior.
              do s = 1, nsc
                spw(s) = P_flat(ci, s)
              end do
              call co_rotot_Rtot_dev ( spw(1:nsc), rho, Rgas )
              rho_wall = P_flat(ci, np) / (Rgas * c_tw(idx))
              do s = 1, nsc
                spw(s) = spw(s) * rho_wall / rho
              end do
              Tint(1) = int(c_tw(idx)); Tint(2) = Tint(1) + 1; Tdiff = c_tw(idx) - Tint(1)
              call co_k_mi_lam_Wilke_dev ( spw(1:nsc), rho_wall, Tint, Tdiff, mil_w, kl_w )
              om_wall = omega_wall_c(1) * mil_w / (dist_w*dist_w)             ! SST omega_wall (rho*w), mil at T_wall
              P_flat(cg, nt)   = - P_flat(ci, nt)                           ! k ghost (k_wall = 0)
              P_flat(cg, nt+1) = rho * 2d0*om_wall/rho_wall - P_flat(ci, nt+1)  ! omega ghost
            else
              om_wall = omega_wall_c(1) * mil_c_f(ci) / (dist_w*dist_w)   ! SST omega_wall (the flat twin of mil_c(Im,Jm,Km))
              P_flat(cg, nt)   = - P_flat(ci, nt)                           ! k ghost (k_wall = 0)
              P_flat(cg, nt+1) = 2d0*om_wall - P_flat(ci, nt+1)             ! omega ghost (rho_wall = rho)
            end if
          end if
      end select
    end do
  end subroutine Extrapolate_Wall_blk

  !> rans-wall-ghost = copy (MOSE_open Ghost_Wall_Extrapolation): after Extrapolate_Wall_blk, the turbulence
  !> variables of every wall ghost take the value of the boundary cell (P(np:nprim) copied in MOSE_open; velocity, species
  !> and pressure are already what MOSE_open writes).
  subroutine Extrapolate_Wall_RANS_Copy_blk ( boff, pim, pjm, c_bi, c_bj, c_bk, c_bf, c_typ, nbnd, lo, hi )
    use MOSE_Parameters_m, only: guide
    use MOSE_Field_m,      only: P_flat
    implicit none
    integer,  intent(in) :: boff, pim, pjm, nbnd, lo, hi
    integer,  intent(in) :: c_bi(nbnd), c_bj(nbnd), c_bk(nbnd), c_bf(nbnd), c_typ(nbnd)
    integer  :: idx, Im, Jm, Km, Fm, Ig, Jg, Kg, ci, cg, s
    if (hi < lo .or. .not. RANS) return
    !$acc parallel loop gang vector present(P_flat) copyin(c_bi, c_bj, c_bk, c_bf, c_typ) &
    !$acc&   private(Im, Jm, Km, Fm, Ig, Jg, Kg, ci, cg, s)
    do idx = lo, hi
      select case (c_typ(idx))
        case (5, 6, 8, 9, 10, 12, 13, 1000)   ! Is_Wall
          Im = c_bi(idx); Jm = c_bj(idx); Km = c_bk(idx); Fm = c_bf(idx)
          Ig = Im - guide(Fm,1); Jg = Jm - guide(Fm,2); Kg = Km - guide(Fm,3)
          ci = boff + (Im+gc-1) + (Jm+gc-1)*pim + (Km+gc-1)*pim*pjm + 1
          cg = boff + (Ig+gc-1) + (Jg+gc-1)*pim + (Kg+gc-1)*pim*pjm + 1
          do s = nt, nprim
            P_flat(cg, s) = P_flat(ci, s)
          end do
      end select
    end do
  end subroutine Extrapolate_Wall_RANS_Copy_blk







  !============================================================================
  ! GPU Phase 2 (F4c.3) — device ghost-cell fill (local boundary entries).
  !
  ! Device twins of the local ghost fills + a per-block kernel. Used by the
  ! SERIAL device path (Fill_Ghost_Cell dispatches here when mpi_size_==1; the
  ! multi-rank MPI exchange stays on the host until F5 brings GPU-aware halos).
  ! Routines mirror Ghost_Symmetry / Ghost_ZG_Extrapolate exactly, reading the
  ! flat fm geometry instead of the AoS dir(d)%f. Ghost writes from different
  ! faces target disjoint ghost cells (corner ghosts are never filled by a face
  ! loop), so a single parallel loop over a block's faces is race-free — matching
  ! the host's single `!$omp do`.
  !============================================================================





  subroutine Fill_Ghost_Cell_mr ( domain )
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Compute/comm overlap -- the same-thread invariant. NVHPC maps OpenACC async
    ! queues PER HOST THREAD: the thread that issues `!$acc wait(Q_OVERLAP)` must be the
    ! SAME thread that LAUNCHED the async(Q_OVERLAP) kernels. Both live inside _start,
    ! and the whole driver runs in ONE `!$omp single`, so the invariant holds.
    !$omp single
    call Fill_Ghost_Cell_mr_start ( domain )
    call Fill_Ghost_Cell_mr_finish ( domain )
    !$omp end single
  end subroutine Fill_Ghost_Cell_mr


  !> E1/rk_async split, part 1: post receives, device-pack + post P sends, local
  !> physical + intra-rank fills (async Q_OVERLAP, fenced), transverse Pg post/pack/send.
  !> On return every MPI transfer is IN FLIGHT and the device is free -- Mod_Explicit
  !> could run ghost-independent work in that window.
  !> NO `!$omp single` inside: the caller owns it (start+finish must share one single).
  subroutine Fill_Ghost_Cell_mr_start ( domain )
    use MOSE_Mod_GhostExchange, only: exchange_ghost_P_post_recv, exchange_ghost_P_post_send, &
                                       exchange_ghost_P_pack, &
                                       exchange_ghost_Pg_post_recv, exchange_ghost_Pg_pack, &
                                       exchange_ghost_Pg_post_send
    use openacc, only: acc_async_sync
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Compute/comm overlap: the ACC queue for the INDEPENDENT local/intra-rank ghost
    ! fill so it overlaps the in-flight inter-rank halo MPI. Any positive id works; it is
    ! joined below before the Pg pack reads the fill-written ghosts.
    integer, parameter :: Q_OVERLAP = 1
    integer :: qs
    ! qs: the synchronous queue; the P pack follows the halo interpolation on it, the Pg pack follows the local fill
    ! on Q_OVERLAP, each with one wait.
    qs = acc_async_sync

    ! (1) Post persistent P receives (device- or host-bound per gpu_aware_mpi_runtime).
    call exchange_ghost_P_post_recv(domain)

    ! (2) Pack the P send buffer ON DEVICE from resident P, then post sends.
    ! dev_packed=.true.: the DEVICE buffer was just packed by exchange_ghost_P_pack and is
    ! authoritative. This SUPPRESSES the H2D send-buffer bridge inside post_send, which would
    ! otherwise overwrite the device-packed bytes with the stale (never-packed) HOST buffer on
    ! the device-bound MPI path -> garbage ghosts from iter 1 ( halo bug, 7.288e3 fixed).
    ! The pack runs on the DEFAULT (synchronous) queue: post_send must not start the MPI send
    ! (CUDA-aware -> DMA from the device P_send_buf) until the pack kernel has completed, and
    ! the default-queue host block guarantees exactly that. (Do NOT make pack async here.)
    call exchange_ghost_P_pack(domain, async_q=qs)
    call exchange_ghost_P_post_send(domain, dev_packed=.true.)

    ! (3) Fill LOCAL physical {3,4,11,22} + INTRA-rank type-1 ghosts on device (resident
    ! P + device bc_pg). Overlap: launch these on ASYNC queue Q_OVERLAP so they run on
    ! the GPU CONCURRENTLY with the in-flight inter-rank P MPI exchange (these local fills are
    ! INDEPENDENT of the arriving halo). The `!$acc wait` below -- issued by THIS SAME thread
    ! -- joins Q before the Pg pack reads the fill-written ghosts and before the unpack/flux
    ! kernels consume them. Same kernels as single-rank, just on a non-default queue.
    ! Fill LOCAL physical + intra-rank ghosts on w (inter-rank w ghosts come from
    ! the device w-unpack in _finish). w interior is current (Compute_dt pack + Newstate +
    ! chemistry); P interior stays current (Newstate/chem unpack) so the P-side pack (step 2)
    ! still reads correct donor interior. With SST_Blk now on w, no P-ghost reader remains.
    call Fill_Ghost_Cell_local ( domain, async_q=Q_OVERLAP, skip_chimera=.true. )   ! chimera runs in _finish, after the unpack

    ! Join (correctness): join Q_OVERLAP HERE -- BEFORE the step-(3b) Pg pack. The Pg
    ! transverse gather reads the donor P at Is+-d11, Js+-d12, ... which, for connection faces
    ! on a donor-block EDGE/CORNER, land in the donor GHOST region -- precisely the cells
    ! step (3) writes. An unfenced pack reads those ghosts mid-flight -> a seam race
    ! (rank-variance ~2.8e-4). This SAME join also covers the fill->unpack corner-cell edge
    ! (unpack is in _finish), so ONE wait suffices, reproducing the  order
    ! fill -> pack -> unpack => BIT-IDENTICAL. Same thread as the async launches above.

    ! (3b)  INTER-rank transverse Pg(:,3:6) -- the SEND side. Post the Pg receives,
    ! gather the donor transverse stencil from the resident donor P into the device Pg send
    ! buffer, then post the Pg sends. The Pg buffers are SEPARATE from the P buffers (no
    ! aliasing). recv MUST be posted before send (matched persistent requests, distinct tag).
    call exchange_ghost_Pg_post_recv()
    call exchange_ghost_Pg_pack(domain, async_q=Q_OVERLAP)   ! joins Q_OVERLAP (local fill, then the pack)
    call exchange_ghost_Pg_post_send()

  end subroutine Fill_Ghost_Cell_mr_start


  !> E1/rk_async split, part 2: wait + unpack P and Pg, receiver-normal bc_pg, complete
  !> sends. NO `!$omp single` inside: the caller owns it.
  subroutine Fill_Ghost_Cell_mr_finish ( domain )
    use MOSE_Mod_GhostExchange, only: exchange_ghost_P_unpack, &
                                       exchange_ghost_P_wait_recv_only, &
                                       exchange_ghost_P_wait_send, &
                                       exchange_ghost_Pg_wait_recv, &
                                       exchange_ghost_Pg_unpack, exchange_ghost_Pg_wait_send, &
                                       fill_bc_pg_normal
    use MOSE_Mod_MPI, only: halo_wait_time
    use openacc, only: acc_async_sync   ! R2
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer(8) :: hw0, hw1, hwrate
    integer :: qs
    ! qs: the synchronous queue of the tail (unpack, chimera, receiver-normal bc_pg, Pg unpack), in this order.
    qs = acc_async_sync

    ! (4) Wait for the P receives, then device-unpack into the resident receiver ghost P
    ! (consumed by the flux kernels). With CUDA-aware MPI the recv buffer is device-resident
    ! and filled in place; the unpack reads it directly -- no host staging of P at all now.
    ! Q_OVERLAP was joined in _start (before step 3b), so the step-(3) fills are complete and
    ! the unpack corner-cell write lands after them (the  fill-then-unpack order).
    ! The system_clock brackets accumulate the true MPI-wait share into halo_wait_time
    ! (per-rank [TIMING] line at run end -- multigpu_plan M0 loss decomposition).
    call system_clock(hw0, hwrate)
    call exchange_ghost_P_wait_recv_only(domain) ! MPI wait: P_recv_buf now valid
    call system_clock(hw1); halo_wait_time = halo_wait_time + real(hw1-hw0,8)/real(hwrate,8)
    call exchange_ghost_P_unpack(domain, async_q=qs)     ! device unpack -> resident P ghosts (kernels)
    ! Chimera ghosts LAST (host / MOSE_open order): donors are local interior cells,
    ! the corner copies land after the unpack, and the Pg pack above saw the previous
    ! substage's chimera ghost exactly as the single-rank path does.
    call Fill_Ghost_Chimera ( domain, qs )
    ! Receiver-normal bc_pg(:,1:gc) ON DEVICE from the just-unpacked receiver ghost P
    ! (replaces the host Ghost_Interrank). Inter-rank recv entries only; intra-rank slots
    ! were already filled by Fill_Ghost_Cell_acc in step (3).
    call fill_bc_pg_normal(domain, async_q=qs)

    ! (5) Wait for the Pg receives, then device-unpack the transverse stencil straight into
    ! the resident bc_pg(:,3:6). No host bc%Pg, no per-slot update device(bc_pg).
    call system_clock(hw0)
    call exchange_ghost_Pg_wait_recv()
    call system_clock(hw1); halo_wait_time = halo_wait_time + real(hw1-hw0,8)/real(hwrate,8)
    call exchange_ghost_Pg_unpack(domain, async_q=qs)

    ! (6) Complete the P and Pg sends (buffers free to be repacked next substage).
    call system_clock(hw0)
    call exchange_ghost_P_wait_send(domain)
    call exchange_ghost_Pg_wait_send()
    call system_clock(hw1); halo_wait_time = halo_wait_time + real(hw1-hw0,8)/real(hwrate,8)

  end subroutine Fill_Ghost_Cell_mr_finish





  !============================================================================
  ! Halo-on-w. Var-LAST flat SoA `w` twins of the device ghost
  ! fill (single-rank path). Same algebra as Ghost_Symmetry_dev /
  ! Ghost_ZG_Extrapolate_dev / Ghost_Connection_pair_acc; only the data access
  ! changes (P(s,i,j,k) -> w(flat_row,s) via WGCELL/WMCELL/WSCELL). These let the
  ! flux/BC w-kernels read correct w ghosts WITHOUT the per-substage field_sync
  ! P -> P_flat pack. STUDY the *_dev / *_acc originals above first.
  !============================================================================

  !> Symmetry ghost on w (acc routine seq): reflect the boundary-cell velocity about
  !> the face normal; copy scalars. Mirrors Ghost_Symmetry_dev.
  !> Rotational periodic face (MOSE_open 200, Lib_Ghost Ghost_Axysimmetry): the ghost is the boundary cell (2D wedge)
  !> or the cell of the opposite k-face (3D) with v and w rotated by +/- delthe, the other variables copied.
  subroutine Ghost_Axysimmetry ( w, ncw, boff, pi, pj, n1, n2, n3, Im, Jm, Km, dir, gstep, ndir_l, delthe )
    !$acc routine seq
    implicit none
    integer,  intent(in)    :: ncw, boff, pi, pj, n1, n2, n3, Im, Jm, Km, dir, gstep, ndir_l
    real(R8), intent(in)    :: delthe
    real(R8), intent(inout) :: w(ncw, nprim)
    integer  :: Ig, Jg, Kg, Ks, s, g3

    Ig = Im; Jg = Jm; Kg = Km; g3 = 0
    select case (dir)
      case (1); Ig = Im - gstep
      case (2); Jg = Jm - gstep
      case (3); Kg = Km - gstep; g3 = gstep   ! guide(Fm,3)
    end select
    if ( ndir_l < 3 ) then
      WGCELL(nv,Ig,Jg,Kg) = WGCELL(nv,Im,Jm,Km) * cos(delthe) + g3 * WGCELL(nw,Im,Jm,Km) * sin(delthe)
      WGCELL(nw,Ig,Jg,Kg) = WGCELL(nv,Im,Jm,Km) * sin(delthe) * (-g3) + WGCELL(nw,Im,Jm,Km) * cos(delthe)
      do s = 1, nu
        WGCELL(s,Ig,Jg,Kg) = WGCELL(s,Im,Jm,Km)
      end do
      do s = np, nprim
        WGCELL(s,Ig,Jg,Kg) = WGCELL(s,Im,Jm,Km)
      end do
    else
      Ks = n3 - Km + 1
      WGCELL(nv,Ig,Jg,Kg) = WGCELL(nv,Im,Jm,Ks) * cos(delthe) + g3 * WGCELL(nw,Im,Jm,Ks) * sin(delthe)
      WGCELL(nw,Ig,Jg,Kg) = WGCELL(nv,Im,Jm,Ks) * sin(delthe) - g3 * WGCELL(nw,Im,Jm,Ks) * cos(delthe)
      do s = 1, nu
        WGCELL(s,Ig,Jg,Kg) = WGCELL(s,Im,Jm,Ks)
      end do
      do s = np, nprim
        WGCELL(s,Ig,Jg,Kg) = WGCELL(s,Im,Jm,Ks)
      end do
    end if
  end subroutine Ghost_Axysimmetry


  subroutine Ghost_Symmetry ( w, ncw, boff, pi, pj, normal, n1, n2, n3, Im, Jm, Km, dir, gstep )
    !$acc routine seq
    implicit none
    integer,  intent(in)    :: ncw, boff, pi, pj, n1, n2, n3, Im, Jm, Km, dir, gstep
    real(R8), intent(in)    :: normal(3)
    real(R8), intent(inout) :: w(ncw, nprim)
    integer  :: Ig, Jg, Kg, s
    real(R8) :: vel(3), udotn

    Ig = Im; Jg = Jm; Kg = Km
    select case (dir)
      case (1); Ig = Im - gstep
      case (2); Jg = Jm - gstep
      case (3); Kg = Km - gstep
    end select

    vel(1) = WGCELL(nu,Im,Jm,Km); vel(2) = WGCELL(nv,Im,Jm,Km); vel(3) = WGCELL(nw,Im,Jm,Km)
    udotn  = vel(1)*normal(1) + vel(2)*normal(2) + vel(3)*normal(3)
    vel(1) = vel(1) - 2d0*udotn*normal(1)
    vel(2) = vel(2) - 2d0*udotn*normal(2)
    vel(3) = vel(3) - 2d0*udotn*normal(3)

    do s = 1, nsc
      WGCELL(s,Ig,Jg,Kg) = WGCELL(s,Im,Jm,Km)
    end do
    WGCELL(nu,Ig,Jg,Kg) = vel(1); WGCELL(nv,Ig,Jg,Kg) = vel(2); WGCELL(nw,Ig,Jg,Kg) = vel(3)
    do s = np, nprim
      WGCELL(s,Ig,Jg,Kg) = WGCELL(s,Im,Jm,Km)
    end do
  end subroutine Ghost_Symmetry

  !> Zero-gradient extrapolation ghost on w (acc routine seq). Mirrors Ghost_ZG_Extrapolate_dev.
  subroutine Ghost_ZG_Extrapolate ( w, ncw, boff, pi, pj, n1, n2, n3, Im, Jm, Km, dir, gstep )
    !$acc routine seq
    implicit none
    integer,  intent(in)    :: ncw, boff, pi, pj, n1, n2, n3, Im, Jm, Km, dir, gstep
    real(R8), intent(inout) :: w(ncw, nprim)
    integer :: g, ig, jg, kg, s
    do g = 1, gc
      ig = Im; jg = Jm; kg = Km
      select case (dir)
        case (1); ig = Im - gstep*g
        case (2); jg = Jm - gstep*g
        case (3); kg = Km - gstep*g
      end select
      do s = 1, nprim
        WGCELL(s,ig,jg,kg) = WGCELL(s,Im,Jm,Km)
      end do
    end do
  end subroutine Ghost_ZG_Extrapolate

  !> 2nd-order extrapolation ghost on w (acc routine seq). Device twin of the host
  !> Ghost_Extrapolate (Fill_Ghost_Cell `case default`): ghost = 3*x3 - 3*x2 + x1.
  !> Used for WALL ghosts (types 5/6/8/9/10/12/13/1000) on the laminar path, where the
  !> host fills the wall ghost with this default extrapolation (Extrapolate_Wall, the
  !> no-slip reflection, is RANS-only -> see Mod_Explicit). The g=1..gc loop is
  !> sequential within this seq routine, so layer g reads the freshly written layer g-1,
  !> byte-matching the host do-loop.
  subroutine Ghost_Extrapolate ( w, ncw, boff, pi, pj, n1, n2, n3, Im, Jm, Km, dir, gstep )
    !$acc routine seq
    implicit none
    integer,  intent(in)    :: ncw, boff, pi, pj, n1, n2, n3, Im, Jm, Km, dir, gstep
    real(R8), intent(inout) :: w(ncw, nprim)
    integer :: g, s, step, ig0, jg0, kg0, i1, j1, k1, i2, j2, k2, i3, j3, k3
    step = -gstep                       ! ghost at +g*step; stencil at (g-1),(g-2),(g-3)*step
    do g = 1, gc
      ig0 = Im; jg0 = Jm; kg0 = Km
      i1  = Im; j1  = Jm; k1  = Km
      i2  = Im; j2  = Jm; k2  = Km
      i3  = Im; j3  = Jm; k3  = Km
      select case (dir)
        case (1)
          ig0 = Im + g*step;  i1 = Im + (g-1)*step;  i2 = Im + (g-2)*step;  i3 = Im + (g-3)*step
        case (2)
          jg0 = Jm + g*step;  j1 = Jm + (g-1)*step;  j2 = Jm + (g-2)*step;  j3 = Jm + (g-3)*step
        case (3)
          kg0 = Km + g*step;  k1 = Km + (g-1)*step;  k2 = Km + (g-2)*step;  k3 = Km + (g-3)*step
      end select
      do s = 1, nprim
        WGCELL(s,ig0,jg0,kg0) = 3d0*WGCELL(s,i1,j1,k1) - 3d0*WGCELL(s,i2,j2,k2) + WGCELL(s,i3,j3,k3)
      end do
    end do
  end subroutine Ghost_Extrapolate

  !> Per-block physical ghost fill on w (symmetry/extrap), twin of Fill_Ghost_Cell_blk_acc.
  subroutine Fill_Ghost_Cell_blk ( w, ncw, boff, pi, pj, fm1, fm2, fm3, n, &
                                     bc_bi, bc_bj, bc_bk, bc_bf, bc_typ, nbnd, lo, hi, async_q, wgm, ndir_l, delthe )
    use openacc, only: acc_async_sync
    implicit none
    integer,  intent(in)    :: ncw, boff, pi, pj, n(3), nbnd, lo, hi
    real(R8), intent(inout) :: w(ncw, nprim)
    real(R8), intent(in)    :: fm1(4, 0:n(1), n(2), n(3))
    real(R8), intent(in)    :: fm2(4, n(1), 0:n(2), n(3))
    real(R8), intent(in)    :: fm3(4, n(1), n(2), 0:n(3))
    integer,  intent(in)    :: bc_bi(nbnd), bc_bj(nbnd), bc_bk(nbnd), bc_bf(nbnd), bc_typ(nbnd)
    integer,  intent(in), optional :: async_q
    integer,  intent(in) :: wgm   ! wall ghost mode: 0 cubic, 1 mirror, 2 zero-gradient
    integer,  intent(in) :: ndir_l   ! number of computational directions (2: wedge faces are k-faces)
    real(R8), intent(in) :: delthe   ! the wedge angle of the rotational periodic faces
    integer  :: idx, Im, Jm, Km, Fm, dir, modfm, gstep, Fi, Fj, Fk, q
    real(R8) :: normal(3)
    ! Dims as scalars -- a dummy array referenced inside the acc region
    ! is implicitly copy()d per launch (upload+download+sync); scalars ride the
    ! launch parameters for free.
    integer :: nb1_, nb2_, nb3_
    nb1_ = n(1); nb2_ = n(2); nb3_ = n(3)

    q = acc_async_sync
    if (present(async_q)) q = async_q

    !$acc parallel loop gang vector present(w, fm1, fm2, fm3, bc_bi, bc_bj, bc_bk, bc_bf, bc_typ) &
    !$acc&   private(Im, Jm, Km, Fm, dir, modfm, gstep, Fi, Fj, Fk, normal)
    do idx = lo, hi
      Im = bc_bi(idx); Jm = bc_bj(idx); Km = bc_bk(idx); Fm = bc_bf(idx)
      dir   = (Fm + 1)/2
      modfm = mod(Fm, 2)
      gstep = 2*modfm - 1
      select case (bc_typ(idx))
        case (3)
          Fi = Im; Fj = Jm; Fk = Km
          select case (dir)
            case (1)
              Fi = Im - modfm
              normal(1) = fm1(1,Fi,Fj,Fk); normal(2) = fm1(2,Fi,Fj,Fk); normal(3) = fm1(3,Fi,Fj,Fk)
            case (2)
              Fj = Jm - modfm
              normal(1) = fm2(1,Fi,Fj,Fk); normal(2) = fm2(2,Fi,Fj,Fk); normal(3) = fm2(3,Fi,Fj,Fk)
            case (3)
              Fk = Km - modfm
              normal(1) = fm3(1,Fi,Fj,Fk); normal(2) = fm3(2,Fi,Fj,Fk); normal(3) = fm3(3,Fi,Fj,Fk)
          end select
          call Ghost_Symmetry ( w, ncw, boff, pi, pj, normal, nb1_, nb2_, nb3_, Im, Jm, Km, dir, gstep )
        case (2)   ! rotational periodic (MOSE_open 200): v and w rotated by +/- delthe, as Lib_Ghost Ghost_Axysimmetry
          call Ghost_Axysimmetry ( w, ncw, boff, pi, pj, nb1_, nb2_, nb3_, Im, Jm, Km, dir, gstep, ndir_l, delthe )
        case (0, 4, 11, 22, 24)   ! 24 = MOSE_open 406 outlet: zero-gradient ghost, as MOSE_open case(0,401:407,420)
          ! type 0 = degenerate 2D k-face (host Fill_Ghost_Cell case(0,...)): zero-gradient.
          call Ghost_ZG_Extrapolate ( w, ncw, boff, pi, pj, nb1_, nb2_, nb3_, Im, Jm, Km, dir, gstep )
        case (5, 6, 8, 9, 10, 12, 13, 1000, 1001)
          ! Wall + manifold(1001): host Fill_Ghost_Cell `case default` -> 2nd-order
          ! extrapolation (the no-slip Extrapolate_Wall reflection is RANS-only). Match it
          ! on the w-halo so the near-boundary reconstruction reads the same ghost as the host.
          ! MOSE_WALL_GHOST: 0 cubic (default, bit-identical), 1 mirror (symmetry reflection),
          ! 2 zero-gradient -- the cubic fill has gain 7 on wall-normal odd-even modes
          ! (doc/THOR_BC_WALL_IS_SYMMETRY.md, the viscous-faceplate blow-up).
          if (wgm == 1) then
          Fi = Im; Fj = Jm; Fk = Km
          select case (dir)
            case (1)
              Fi = Im - modfm
              normal(1) = fm1(1,Fi,Fj,Fk); normal(2) = fm1(2,Fi,Fj,Fk); normal(3) = fm1(3,Fi,Fj,Fk)
            case (2)
              Fj = Jm - modfm
              normal(1) = fm2(1,Fi,Fj,Fk); normal(2) = fm2(2,Fi,Fj,Fk); normal(3) = fm2(3,Fi,Fj,Fk)
            case (3)
              Fk = Km - modfm
              normal(1) = fm3(1,Fi,Fj,Fk); normal(2) = fm3(2,Fi,Fj,Fk); normal(3) = fm3(3,Fi,Fj,Fk)
          end select
          call Ghost_Symmetry ( w, ncw, boff, pi, pj, normal, nb1_, nb2_, nb3_, Im, Jm, Km, dir, gstep )
          else if (wgm == 2) then
            call Ghost_ZG_Extrapolate ( w, ncw, boff, pi, pj, nb1_, nb2_, nb3_, Im, Jm, Km, dir, gstep )
          else
            call Ghost_Extrapolate ( w, ncw, boff, pi, pj, nb1_, nb2_, nb3_, Im, Jm, Km, dir, gstep )
          end if
      end select
    end do
  end subroutine Fill_Ghost_Cell_blk

  !> extrapolation-ghost = cubic: the extrapolation faces (400 -> 11) take the cubic extrapolation ghost of MOSE_open Lib_Ghost
  !> (case default) after the standard fill left the legacy zero-gradient one. A separate launch over the block's
  !> records, so the standard kernel's code (and the production gate on it) is untouched. Reads interior cells only.
  subroutine Fill_Ghost_Cubic_Extrapolation_blk ( w, ncw, boff, pi, pj, n, bc_bi, bc_bj, bc_bk, bc_bf, bc_typ, nbnd, lo, hi, async_q )
    use openacc, only: acc_async_sync
    implicit none
    integer,  intent(in)    :: ncw, boff, pi, pj, n(3), nbnd, lo, hi
    real(R8), intent(inout) :: w(ncw, nprim)
    integer,  intent(in)    :: bc_bi(nbnd), bc_bj(nbnd), bc_bk(nbnd), bc_bf(nbnd), bc_typ(nbnd)
    integer,  intent(in), optional :: async_q
    integer :: idx, Im, Jm, Km, Fm, dir, modfm, gstep, q, nb1_, nb2_, nb3_
    q = acc_async_sync
    if (present(async_q)) q = async_q
    nb1_ = n(1); nb2_ = n(2); nb3_ = n(3)
    if (hi < lo) return
    !$acc parallel loop gang vector present(w, bc_bi, bc_bj, bc_bk, bc_bf, bc_typ) &
    !$acc&   private(Im, Jm, Km, Fm, dir, modfm, gstep)
    do idx = lo, hi
      if (bc_typ(idx) /= 11) cycle
      Im = bc_bi(idx); Jm = bc_bj(idx); Km = bc_bk(idx); Fm = bc_bf(idx)
      dir = (Fm + 1)/2; modfm = mod(Fm, 2); gstep = 2*modfm - 1
      call Ghost_Extrapolate ( w, ncw, boff, pi, pj, nb1_, nb2_, nb3_, Im, Jm, Km, dir, gstep )
    end do
  end subroutine Fill_Ghost_Cubic_Extrapolation_blk


  !> Cross-block connection ghost fill on w for ONE ordered (Bm,Bs) pair. Reads donor w
  !> (WSCELL, block-Bs offset), writes receiver w ghost (WMCELL, block-Bm offset) + bc_pg.
  !> Twin of Ghost_Connection_pair_acc; the two block slices live in the one P_flat via
  !> distinct offsets (boff_m/pi_m/pj_m vs boff_s/pi_s/pj_s).
  subroutine Ghost_Connection_pair ( w, ncw, boff_m, pi_m, pj_m, boff_s, pi_s, pj_s, &
                                       bc_pg, nbnd, nconn, c_im, c_jm, c_km, c_fm,      &
                                       c_is, c_js, c_ks, c_fs,                          &
                                       c_d11, c_d12, c_d21, c_d22, c_pg, lo, hi, async_q )
    use MOSE_Parameters_m, only: guide
    use MOSE_Lib_RANS,     only: eddy_dev
    use openacc, only: acc_async_sync
    implicit none
    integer,  intent(in)    :: ncw, boff_m, pi_m, pj_m, boff_s, pi_s, pj_s, nbnd, nconn, lo, hi
    real(R8), intent(inout) :: w(ncw, nprim)
    real(R8), intent(inout) :: bc_pg(nprim, 6, nbnd)
    integer,  intent(in)    :: c_im(nconn), c_jm(nconn), c_km(nconn), c_fm(nconn)
    integer,  intent(in)    :: c_is(nconn), c_js(nconn), c_ks(nconn), c_fs(nconn)
    integer,  intent(in)    :: c_d11(nconn), c_d12(nconn), c_d21(nconn), c_d22(nconn), c_pg(nconn)
    integer,  intent(in), optional :: async_q
    integer :: c, g, Im, Jm, Km, Fm, Is, Js, Ks, Fs, idx, d11, d12, d21, d22, q, s
    integer :: Ig, Jg, Kg, Isr, Jsr, Ksr
    integer :: i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4
    real(R8) :: v1_, v2_, t3_, t4_, t5_, t6_   ! 8c-A

    q = acc_async_sync
    if (present(async_q)) q = async_q

    if ( boff_m /= boff_s ) then   ! 8c-A: one thread per (record, variable), loads before stores
    !$acc parallel loop gang vector collapse(2) present(w, bc_pg, c_im, c_jm, c_km, c_fm, &
    !$acc&   c_is, c_js, c_ks, c_fs, c_d11, c_d12, c_d21, c_d22, c_pg) &
    !$acc&   private(g, Im, Jm, Km, Fm, Is, Js, Ks, Fs, idx, d11, d12, d21, d22, &
    !$acc&           Ig, Jg, Kg, Isr, Jsr, Ksr, i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4, &
    !$acc&           v1_, v2_, t3_, t4_, t5_, t6_)
    do c = lo, hi
    do s = 1, nprim
      Im = c_im(c); Jm = c_jm(c); Km = c_km(c); Fm = c_fm(c)
      Is = c_is(c); Js = c_js(c); Ks = c_ks(c); Fs = c_fs(c)
      idx = c_pg(c)

      d11 = c_d11(c); d12 = c_d12(c); d21 = c_d21(c); d22 = c_d22(c)
      select case (Fs)
        case (1:2)
          i1=Is; j1=Js-d11; k1=Ks-d12;  i2=Is; j2=Js+d11; k2=Ks+d12
          i3=Is; j3=Js-d21; k3=Ks-d22;  i4=Is; j4=Js+d21; k4=Ks+d22
        case (3:4)
          i1=Is-d11; j1=Js; k1=Ks-d12;  i2=Is+d11; j2=Js; k2=Ks+d12
          i3=Is-d21; j3=Js; k3=Ks-d22;  i4=Is+d21; j4=Js; k4=Ks+d22
        case (5:6)
          i1=Is-d11; j1=Js-d12; k1=Ks;  i2=Is+d11; j2=Js+d12; k2=Ks
          i3=Is-d21; j3=Js-d22; k3=Ks;  i4=Is+d21; j4=Js+d22; k4=Ks
      end select
      ! loads first (independent, all in flight together) ...
      v1_ = WSCELL(s,Is,Js,Ks)
      v2_ = 0.0_R8
      if ( gc >= 2 ) v2_ = WSCELL(s,Is+guide(Fs,1),Js+guide(Fs,2),Ks+guide(Fs,3))
      t3_ = WSCELL(s,i1,j1,k1); t4_ = WSCELL(s,i2,j2,k2)
      t5_ = WSCELL(s,i3,j3,k3); t6_ = WSCELL(s,i4,j4,k4)
      ! ... then the stores
      WMCELL(s,Im-guide(Fm,1),Jm-guide(Fm,2),Km-guide(Fm,3)) = v1_
      bc_pg(s,1,idx) = v1_
      if ( gc >= 2 ) then
        WMCELL(s,Im-2*guide(Fm,1),Jm-2*guide(Fm,2),Km-2*guide(Fm,3)) = v2_
        bc_pg(s,2,idx) = v2_
      end if
      do g = 3, gc   ! deeper ghost layers (gc > 2): the original copy, in order
        Ig  = Im - guide(Fm,1)*g;     Jg  = Jm - guide(Fm,2)*g;     Kg  = Km - guide(Fm,3)*g
        Isr = Is + guide(Fs,1)*(g-1); Jsr = Js + guide(Fs,2)*(g-1); Ksr = Ks + guide(Fs,3)*(g-1)
        WMCELL(s,Ig,Jg,Kg) = WSCELL(s,Isr,Jsr,Ksr)
        bc_pg(s,g,idx)     = WMCELL(s,Ig,Jg,Kg)
      end do
      bc_pg(s,3,idx) = t3_; bc_pg(s,4,idx) = t4_
      bc_pg(s,5,idx) = t5_; bc_pg(s,6,idx) = t6_
    end do
    end do
    return
    end if

    ! self pair (Bm == Bs): the original kernel, in-thread write -> read order preserved
    !$acc parallel loop gang vector present(w, bc_pg, c_im, c_jm, c_km, c_fm, &
    !$acc&   c_is, c_js, c_ks, c_fs, c_d11, c_d12, c_d21, c_d22, c_pg) &
    !$acc&   private(g, s, Im, Jm, Km, Fm, Is, Js, Ks, Fs, idx, d11, d12, d21, d22, &
    !$acc&           Ig, Jg, Kg, Isr, Jsr, Ksr, i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4)
    do c = lo, hi
      Im = c_im(c); Jm = c_jm(c); Km = c_km(c); Fm = c_fm(c)
      Is = c_is(c); Js = c_js(c); Ks = c_ks(c); Fs = c_fs(c)
      idx = c_pg(c)
      do g = 1, gc
        Ig  = Im - guide(Fm,1)*g;     Jg  = Jm - guide(Fm,2)*g;     Kg  = Km - guide(Fm,3)*g
        Isr = Is + guide(Fs,1)*(g-1); Jsr = Js + guide(Fs,2)*(g-1); Ksr = Ks + guide(Fs,3)*(g-1)
        do s = 1, nprim
          WMCELL(s,Ig,Jg,Kg) = WSCELL(s,Isr,Jsr,Ksr)
          bc_pg(s,g,idx)     = WMCELL(s,Ig,Jg,Kg)
        end do
      end do
      d11 = c_d11(c); d12 = c_d12(c); d21 = c_d21(c); d22 = c_d22(c)
      select case (Fs)
        case (1:2)
          i1=Is; j1=Js-d11; k1=Ks-d12;  i2=Is; j2=Js+d11; k2=Ks+d12
          i3=Is; j3=Js-d21; k3=Ks-d22;  i4=Is; j4=Js+d21; k4=Ks+d22
        case (3:4)
          i1=Is-d11; j1=Js; k1=Ks-d12;  i2=Is+d11; j2=Js; k2=Ks+d12
          i3=Is-d21; j3=Js; k3=Ks-d22;  i4=Is+d21; j4=Js; k4=Ks+d22
        case (5:6)
          i1=Is-d11; j1=Js-d12; k1=Ks;  i2=Is+d11; j2=Js+d12; k2=Ks
          i3=Is-d21; j3=Js-d22; k3=Ks;  i4=Is+d21; j4=Js+d22; k4=Ks
      end select
      do s = 1, nprim
        bc_pg(s,3,idx) = WSCELL(s,i1,j1,k1); bc_pg(s,4,idx) = WSCELL(s,i2,j2,k2)
        bc_pg(s,5,idx) = WSCELL(s,i3,j3,k3); bc_pg(s,6,idx) = WSCELL(s,i4,j4,k4)
      end do
    end do
  end subroutine Ghost_Connection_pair


  !> R2 phase B: ONE launch per level over the same union: the four transversal neighbours of the donor face cell into
  !> bc_pg(:,3:6) -- MOSE_open's pass 2 -- read from ghosts that phase A, the physical fill and the chimera fill have ALL
  !> written by now. Reads w only, writes bc_pg only: no hazard, no order.
  subroutine Ghost_Connection_pairs_B ( w, ncw, bc_pg, nbnd, nconn, c_is, c_js, c_ks, c_fs, c_d11, c_d12, c_d21, c_d22, c_pg, &
                                          cu, cbs, nu, async_q )
    use MOSE_Field_m, only: dboff, dpi, dpj
    use openacc, only: acc_async_sync
    implicit none
    integer,  intent(in)    :: ncw, nbnd, nconn, nu
    integer,  intent(in)    :: cu(nu), cbs(nu)
    real(R8), intent(in)    :: w(ncw, nprim)
    real(R8), intent(inout) :: bc_pg(nprim, 6, nbnd)
    integer,  intent(in)    :: c_is(nconn), c_js(nconn), c_ks(nconn), c_fs(nconn)
    integer,  intent(in)    :: c_d11(nconn), c_d12(nconn), c_d21(nconn), c_d22(nconn), c_pg(nconn)
    integer,  intent(in), optional :: async_q
    integer :: e, c, s, Is, Js, Ks, Fs, idx, d11, d12, d21, d22, q
    integer :: i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4, boff_s, pi_s, pj_s
    q = acc_async_sync
    if (present(async_q)) q = async_q
    !$acc parallel loop gang vector collapse(2) present(w, bc_pg, c_is, c_js, c_ks, c_fs, c_d11, c_d12, c_d21, c_d22, c_pg, &
    !$acc&   cu, cbs, dboff, dpi, dpj) &
    !$acc&   private(c, Is, Js, Ks, Fs, idx, d11, d12, d21, d22, i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4, boff_s, pi_s, pj_s)
    do e = 1, nu
    do s = 1, nprim
      c = cu(e); boff_s = dboff(cbs(e)); pi_s = dpi(cbs(e)); pj_s = dpj(cbs(e))
      Is = c_is(c); Js = c_js(c); Ks = c_ks(c); Fs = c_fs(c)
      idx = c_pg(c)
      d11 = c_d11(c); d12 = c_d12(c); d21 = c_d21(c); d22 = c_d22(c)
      select case (Fs)
        case (1:2)
          i1=Is; j1=Js-d11; k1=Ks-d12;  i2=Is; j2=Js+d11; k2=Ks+d12
          i3=Is; j3=Js-d21; k3=Ks-d22;  i4=Is; j4=Js+d21; k4=Ks+d22
        case (3:4)
          i1=Is-d11; j1=Js; k1=Ks-d12;  i2=Is+d11; j2=Js; k2=Ks+d12
          i3=Is-d21; j3=Js; k3=Ks-d22;  i4=Is+d21; j4=Js; k4=Ks+d22
        case (5:6)
          i1=Is-d11; j1=Js-d12; k1=Ks;  i2=Is+d11; j2=Js+d12; k2=Ks
          i3=Is-d21; j3=Js-d22; k3=Ks;  i4=Is+d21; j4=Js+d22; k4=Ks
      end select
      bc_pg(s,3,idx) = WSCELL(s,i1,j1,k1); bc_pg(s,4,idx) = WSCELL(s,i2,j2,k2)
      bc_pg(s,5,idx) = WSCELL(s,i3,j3,k3); bc_pg(s,6,idx) = WSCELL(s,i4,j4,k4)
    end do
    end do
  end subroutine Ghost_Connection_pairs_B


  !> Domain dispatch for halo on P_flat (single-rank), twin of Fill_Ghost_Cell_acc. Fetches each
  !> block's flat-index scalars on the host and threads them in; w is field-resident.
  !> Chimera (666) ghost fill on w: per entry, the two ghost rows are volume-fraction
  !> weighted sums of donor cells in CONSERVED variables, converted back with the
  !> receptor cell's temperature as the Newton guess -- the exact device twin of the
  !> host Ghost_Chimera. All cell references are GLOBAL P_flat indices precomputed by
  !> Setup_BC_Device_Flat (cidx), so donors may live in any local block with no
  !> per-pair dispatch. Writes the whole ghost state like the host: flow (1:np)
  !> through the conservative blend, RANS/passive slots (np+1:nprim) through the
  !> density-weighted blend (MOSE_open Ghost_Chimera). Then, as the host
  !> Fill_BC_Ghost_Chimera, the transverse viscous stencil bc_pg(:,3:6) from the
  !> boundary cell's four transverse neighbours (tn) and the block-edge corner
  !> ghost-of-ghost cells (cg) as copies of ghost row 1 -- without these the viscous
  !> connection twin reads zero-density cells (0/0 -> NaN) on every chimera face.
  subroutine Ghost_Chimera ( w, ncw, bc_pg, nbnd, nchim, g1, g2, rc, dlo, dmid, dhi, &
                               slot, ndon, df, vf, tn, cg, async_q, elo, ehi )
    use MOSE_Global_m, only: NPRIM_DEV
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, prim2cons_dev, cons2prim_dev
    use openacc, only: acc_async_sync
    implicit none
    integer,  intent(in)    :: ncw, nbnd, nchim, ndon
    real(R8), intent(inout) :: w(ncw, nprim)
    real(R8), intent(inout) :: bc_pg(nprim, 6, nbnd)
    integer,  intent(in)    :: g1(nchim), g2(nchim), rc(nchim), dlo(nchim), dmid(nchim), dhi(nchim), slot(nchim)
    integer,  intent(in)    :: df(ndon)
    real(R8), intent(in)    :: vf(ndon)
    integer,  intent(in)    :: tn(4,nchim), cg(4,nchim)
    integer,  intent(in), optional :: async_q
    integer,  intent(in), optional :: elo, ehi          ! 5g: entry range to sweep
    integer  :: e, c, s, q, idx, row, clo, chi, gcell, e0, e1
    real(R8) :: rho, Rgas, Tg, rho_d, rho_g
    real(R8) :: pr(NPRIM_DEV), consd(NPRIM_DEV), consg(NPRIM_DEV), pg(NPRIM_DEV)

    q = acc_async_sync
    if (present(async_q)) q = async_q
    e0 = 1; e1 = nchim
    if (present(elo)) e0 = elo
    if (present(ehi)) e1 = ehi

    !$acc parallel loop gang vector present(w, bc_pg, g1, g2, rc, dlo, dmid, dhi, slot, df, vf, tn, cg) &
    !$acc&   private(e, c, s, idx, row, clo, chi, gcell, rho, Rgas, Tg, rho_d, rho_g, pr, consd, consg, pg)
    do e = e0, e1
      ! receptor-cell temperature guess (host: co_rotot_Rtot + EOS on blk%P(...,im,jm,km))
      do s = 1, nsc
        pr(s) = w(rc(e), s)
      end do
      call co_rotot_Rtot_dev ( pr(1:nsc), rho, Rgas )
      Tg  = w(rc(e), np) / (rho*Rgas)
      idx = slot(e)
      do row = 1, 2
        if (row == 1) then
          clo = dlo(e); chi = dmid(e); gcell = g1(e)
        else
          clo = dmid(e) + 1; chi = dhi(e); gcell = g2(e)
        end if
        do s = 1, nprim
          consg(s) = 0d0
        end do
        do c = clo, chi
          do s = 1, nprim
            pr(s) = w(df(c), s)
          end do
          call prim2cons_dev ( pr(1:np), consd(1:np) )
          rho_d = 0d0
          do s = 1, nsc
            rho_d = rho_d + pr(s)
          end do
          do s = 1, np
            consg(s) = consg(s) + consd(s)*vf(c)
          end do
          do s = np+1, nprim                              ! RANS/passive: density-weighted blend
            consg(s) = consg(s) + pr(s)*rho_d*vf(c)
          end do
        end do
        call cons2prim_dev ( consg(1:np), Tg, pg(1:np) )
        rho_g = 0d0
        do s = 1, nsc
          rho_g = rho_g + pg(s)
        end do
        do s = 1, np
          w(gcell, s)      = pg(s)
          bc_pg(s,row,idx) = pg(s)
        end do
        do s = np+1, nprim
          w(gcell, s)      = consg(s) / rho_g
          bc_pg(s,row,idx) = w(gcell, s)
        end do
      end do
      ! transverse viscous stencil: the boundary cell's four transverse neighbours
      do s = 1, nprim
        bc_pg(s,3,idx) = w(tn(1,e), s); bc_pg(s,4,idx) = w(tn(2,e), s)
        bc_pg(s,5,idx) = w(tn(3,e), s); bc_pg(s,6,idx) = w(tn(4,e), s)
      end do
      ! block-edge corner ghosts of ghost row 1 (each written by exactly one entry)
      do c = 1, 4
        if (cg(c,e) > 0) then
          do s = 1, nprim
            w(cg(c,e), s) = w(g1(e), s)
          end do
        end if
      end do
    end do
  end subroutine Ghost_Chimera



  !> [MOSE-scheme] rans-wall-ghost = copy (default) | wall-value: see rans_wall_copy. Called by IO_Ini at start-up and at
  !> every re-read of input.ini.
  subroutine Assign_RANS_Wall_Ghost ( key_in )
    implicit none
    character(len=*), intent(in) :: key_in
    select case ( trim(adjustl(key_in)) )
      case ('wall-value', 'WALL-VALUE'); rans_wall_copy = .false.
      case ('copy', 'COPY');             rans_wall_copy = .true.
      case default
        write(*,'(3A)') ' [ERROR] [MOSE-scheme] rans-wall-ghost = ', trim(key_in), ': copy (default) or wall-value'
        error stop 'Assign_RANS_Wall_Ghost: unknown rans-wall-ghost'
    end select
  end subroutine Assign_RANS_Wall_Ghost


  !> [MOSE-scheme] extrapolation-ghost = cubic (default, MOSE_open) | zero-gradient: the ghost cells of the extrapolation faces
  !> (400 -> 11), see extrap_ghost_cubic. Called by IO_Ini at start-up and at every re-read of input.ini.
  subroutine Assign_Extrap_Ghost ( key_in )
    use MOSE_Mod_MPI, only: mpi_is_root
    implicit none
    character(len=*), intent(in) :: key_in
    character(len=len(key_in)) :: key
    integer :: i, c
    logical, save :: told = .false.
    key = adjustl(key_in)
    do i = 1, len_trim(key)
      c = iachar(key(i:i)); if (c >= iachar('A') .and. c <= iachar('Z')) key(i:i) = achar(c + 32)
    end do
    select case ( trim(key) )
      case ('zero-gradient', 'zg'); extrap_ghost_cubic = .false.
      case ('cubic');               extrap_ghost_cubic = .true.
      case default
        write(*,'(3A)') ' [ERROR] [MOSE-scheme] extrapolation-ghost = ', trim(key_in), ': zero-gradient (default) or cubic'
        error stop 'Assign_Extrap_Ghost: unknown extrapolation-ghost'
    end select
    if ( extrap_ghost_cubic .and. mpi_is_root .and. .not. told ) &
      write(*,'(A)') ' [BC] extrapolation faces (400): cubic extrapolation ghost cells, as MOSE_open Lib_Ghost'
    told = told .or. extrap_ghost_cubic
  end subroutine Assign_Extrap_Ghost


  subroutine wall_ghost_mode_init ( domain )
    use MOSE_Mod_MPI, only: mpi_is_root
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    character(len=16) :: wenv, key
    integer :: wst, i, c, nvisc, nty(7)
    integer, parameter :: vty(7) = [5, 6, 8, 9, 10, 12, 13]
    character(len=160) :: census, prev
    call get_environment_variable('MOSE_WALL_GHOST', wenv, status=wst)
    if (wst == 0 .and. mpi_is_root) write(*,'(3A)') ' [BC] warning: MOSE_WALL_GHOST=', trim(wenv), &
      ' is ignored -- use the key wall-ghost = auto | mirror | zg | cubic in [MOSE-scheme]'
    key = adjustl(wall_ghost_key)
    do i = 1, len_trim(key)
      c = iachar(key(i:i)); if (c >= iachar('A') .and. c <= iachar('Z')) key(i:i) = achar(c + 32)
    end do
    nty = 0       ! domain%bc holds every record of bc.txt on every rank: the same decision everywhere, no collective
    do i = 1, domain % nbound
      do c = 1, size(vty)
        if ( domain % bc(i) % type == vty(c) ) nty(c) = nty(c) + 1
      end do
    end do
    nvisc = sum(nty)
    census = ''
    do c = 1, size(vty)
      if ( nty(c) > 0 ) then
        prev = census     ! the internal file must not appear in its own output list
        write(census, '(A,1X,I0,A,I0)') trim(prev), nty(c), ' x type ', vty(c)
      end if
    end do
    select case ( trim(key) )
      case ('auto');   wall_ghost_mode = merge(1, 0, nvisc > 0)
      case ('mirror'); wall_ghost_mode = 1
      case ('zg');     wall_ghost_mode = 2
      case ('cubic');  wall_ghost_mode = 0
      case default
        write(*,'(3A)') ' [BC] ERROR: wall-ghost = ', trim(key), ' in [MOSE-scheme] -- the values are auto, mirror, zg or cubic'
        error stop 'invalid wall-ghost value (see the [BC] ERROR line)'
    end select
    if (mpi_is_root) write(*,'(3A,I0,A,I0,2A)') ' [BC] wall-ghost ', trim(key), ' -> mode ', wall_ghost_mode, &
      ' (0 cubic, 1 mirror, 2 zero-gradient); viscous wall records ', nvisc, ':', trim(census)
    wgm_checked = .true.
  end subroutine wall_ghost_mode_init


  subroutine Fill_Ghost_Cell_local ( domain, async_q, skip_chimera )
    use openacc, only: acc_async_sync
    use MOSE_Lib_Metrics, only: delthe
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Field_m,  only: P_flat, ncell_total, field_block_index_params
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in), optional :: async_q
    logical, intent(in), optional :: skip_chimera
    logical :: skipc
    integer :: b, p, Bm, Bs, q, ncw, boff, pi, pj, boff_m, pi_m, pj_m, boff_s, pi_s, pj_s

    q = acc_async_sync
    if (.not. wgm_checked) call wall_ghost_mode_init ( domain )
    if (present(async_q)) q = async_q
    skipc = .false.
    if (present(skip_chimera)) skipc = skip_chimera
    ncw = ncell_total

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      if (domain % bc_blk_hi(b) < domain % bc_blk_lo(b)) cycle
      call field_block_index_params(b, boff, pi, pj)
      call Fill_Ghost_Cell_blk ( P_flat, ncw, boff, pi, pj, &
                                   domain % blk(b) % fm1, domain % blk(b) % fm2, domain % blk(b) % fm3, &
                                   domain % blk(b) % dim, &
                                   domain % bc_bi, domain % bc_bj, domain % bc_bk, domain % bc_bf, domain % bc_typ, &
                                   domain % nbound, domain % bc_blk_lo(b), domain % bc_blk_hi(b), q , wall_ghost_mode, ndir, delthe )
      if ( extrap_ghost_cubic ) &   ! 400 faces: cubic extrapolation ghost of MOSE_open (a separate launch: the kernel above is untouched)
      call Fill_Ghost_Cubic_Extrapolation_blk ( P_flat, ncw, boff, pi, pj, domain % blk(b) % dim, &
                                   domain % bc_bi, domain % bc_bj, domain % bc_bk, domain % bc_bf, domain % bc_typ, &
                                   domain % nbound, domain % bc_blk_lo(b), domain % bc_blk_hi(b), q )
    end do

    do p = 1, domain % n_conn_pair
      Bm = domain % conn_pair_bm(p); Bs = domain % conn_pair_bs(p)
      call field_block_index_params(Bm, boff_m, pi_m, pj_m)
      call field_block_index_params(Bs, boff_s, pi_s, pj_s)
      call Ghost_Connection_pair ( P_flat, ncw, boff_m, pi_m, pj_m, boff_s, pi_s, pj_s, &
                                     domain % bc_pg, domain % nbound, domain % n_conn, &
                                     domain % conn_im, domain % conn_jm, domain % conn_km, domain % conn_fm, &
                                     domain % conn_is, domain % conn_js, domain % conn_ks, domain % conn_fs, &
                                     domain % conn_d11, domain % conn_d12, domain % conn_d21, domain % conn_d22, &
                                     domain % conn_pg, domain % conn_pair_lo(p), domain % conn_pair_hi(p), q )
    end do
    ! MOSE_open Fill_Ghost_Cell pass 2: each pair kernel gathers the four transverse neighbours of its donor face cell into
    ! bc_pg(:,3:6) in the launch that writes its own ghosts, so a neighbour in the ghost strip of a LATER pair held the
    ! previous stage's value. Gathered again here, ONE launch over every entry, after all the ghost writes above (physical
    ! fill + every pair) and before the chimera fill -- the order of MOSE_open (and of the host Fill_Ghost_Cell).
    if ( domain % n_conn > 0 ) &
    call Ghost_Connection_pairs_B ( P_flat, ncw, domain % bc_pg, domain % nbound, domain % n_conn, &
                                      domain % conn_is, domain % conn_js, domain % conn_ks, domain % conn_fs, &
                                      domain % conn_d11, domain % conn_d12, domain % conn_d21, domain % conn_d22, &
                                      domain % conn_pg, domain % conn_all, domain % conn_all_bs, domain % n_conn, q )

    if (.not. skipc) call Fill_Ghost_Chimera ( domain, q )
  end subroutine Fill_Ghost_Cell_local


  !> Chimera ghost fill on w (device), split out so the multi-rank halo can run it
  !> AFTER the inter-rank unpack -- the order of the host path and of MOSE_open
  !> (Fill_BC_Ghost_Connection's transverse stencil is gathered BEFORE Ghost_Chimera,
  !> so a donor-edge stencil sees the previous substage's chimera ghost). Running it
  !> before the Pg pack instead made the 2-rank THOR run differ from the 1-rank one
  !> at the ULP level (2026-08-25). Single-rank order unchanged (it is called last).
  subroutine Fill_Ghost_Chimera ( domain, async_q )
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Field_m,  only: P_flat, ncell_total, field_block_index_params
    use openacc, only: acc_async_sync
        implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in), optional :: async_q
    integer :: q, ncw
    q = acc_async_sync
    if (present(async_q)) q = async_q
    ncw = ncell_total
    if (domain % n_chim > 0) then
      call Ghost_Chimera ( P_flat, ncw, domain % bc_pg, domain % nbound, domain % n_chim, &
                             domain % chim_g1, domain % chim_g2, domain % chim_rc, &
                             domain % chim_dlo, domain % chim_dmid, domain % chim_dhi, &
                             domain % chim_slot, size(domain % chim_df), &
                             domain % chim_df, domain % chim_vf, domain % chim_tn, domain % chim_cg, q, 1, domain % n_chim )
    end if
  end subroutine Fill_Ghost_Chimera


end module MOSE_Lib_Extrapolation
