!> The flat field: the var-last device twins of the block arrays.
!>
!> MOSE_open stores the state var-first per block, blk(b)%P(nprim, 1-gc:ni+gc, ...). On the GPU that layout is
!> uncoalesced (consecutive threads sweeping i for one variable read with stride nprim). This module owns the
!> var-last twins
!>
!>     P_flat(cell, var), PO_flat(cell, var), R_flat(cell, var)
!>
!> where `cell` runs over every padded cell (1-gc:ni+gc x ...) of every local block, concatenated through a per-block
!> offset table (blk_off): P_flat(cidx(b,i,j,k), var) over consecutive i is contiguous, hence coalesced (the STREAmS
!> layout). pack_P_to_flat / unpack_flat_to_P convert; P is the truth at setup and at the I/O boundary of a step,
!> P_flat inside the step. The per-cell thermodynamic caches have flat twins too (rho_c_f, T_c_f, Rg_c_f, mil_c_f,
!> hr_f).
module MOSE_Field_m
  use iso_fortran_env, only: R8 => real64
  use MOSE_Global_m,   only: gc, nprim
  use MOSE_Types_m,    only: MOSE_domain_type
  use MOSE_Mod_MPI,    only: is_local_block

  implicit none
  private

  public :: field_alloc, field_is_allocated
  public :: pack_P_to_flat, unpack_flat_to_P, field_roundtrip_check
  public :: field_update_P_flat_device, field_block_index_params, field_sync_P_to_flat
  public :: field_copy_P_flat_to_PO_flat, field_update_PO_flat_host, unpack_flat_to_PO
  public :: field_state_hash
  public :: rho_c_f, T_c_f, n_cache_f, field_cache_flat_alloc
  public :: Rg_c_f, mil_c_f                                    ! the flat twins of the Rg and viscosity caches
  public :: hr_f                                               ! the flat twin of hr_c (|dT|/T of the last chemistry step)
  public :: vol_f, dt_f
  public :: dboff, dpi, dpj
  public :: P_flat, cidx, ncell_total
  !> The var-last twins of PO and R, sharing cidx/blk_off/ncell_total.
  public :: PO_flat, R_flat
  public :: field_update_P_flat_host
  !> R_flat is the sole residual accumulator of the device path; the host sources (RotatingFrame / soot) fold
  !> their var-first r into it through field_add_R_to_R_flat, and the var-first blk%r is rebuilt from it only
  !> at the I/O / diagnostic boundary (unpack_flat_to_R_host).
  public :: field_update_R_flat_host, unpack_flat_to_R_host, field_add_R_to_R_flat

  ! Global flat var-last state SoA. (ncell_total, nprim).
  real(R8), allocatable, save :: P_flat(:,:)
  !> The old state (PO) and the residual (R) in the same var-last flat layout, read coalesced by the Newstate
  !> kernel; device-resident like P_flat.
  real(R8), allocatable, save :: PO_flat(:,:)
  real(R8), allocatable, save :: R_flat(:,:)
  ! 0-based cell offset of block b in P_flat (local blocks only; 0 width for non-local).
  integer,  allocatable, save :: blk_off(:)
  ! Interior dims (ni,nj,nk) per block, cached so cidx/lin are self-contained.
  integer,  allocatable, save :: bdim(:,:)
  !> blk_off/bdim are HOST-only. The device w-kernel does not read them: it fetches
  !> per-block scalars via field_block_index_params on the host and computes the flat
  !> w index inline (module allocatables are unreliable inside acc routine seq -> a
  !> device illegal address; cidx/lin therefore stay host-only, used by pack/unpack).
  integer,  save :: ncell_total = 0
  logical,  save :: allocated_flag = .false.

  !> The flat twins of the per-cell thermo caches, on the same padded arena and the same flat
  !> index as P_flat. The fused diffusive kernel reads rho_c and T_c at STENCIL neighbours, which are
  !> not in the level's list, so those two cannot be gathered by list position like the geometry --
  !> they need a row-addressed array. Allocated with the flat field, written
  !> by the thermo cache kernels alongside the per-block copies.
  !> The per-block flat-index scalars as device tables. A kernel that spans blocks cannot
  !> take boff/pi/pj as launch parameters; it looks them up from its element's own block.
  integer, allocatable, save :: dboff(:), dpi(:), dpj(:)
  real(R8), allocatable, save :: rho_c_f(:), T_c_f(:)
  real(R8), allocatable, save :: Rg_c_f(:), mil_c_f(:)   !< Rg and laminar viscosity twins, written wherever rho_c_f is
  real(R8), allocatable, save :: hr_f(:)                 !< |dT|/T of the last chemistry step, twin of blk%hr_c
  real(R8), allocatable, save :: vol_f(:)   !< 5s-a: cell volume, constant, filled once
  real(R8), allocatable, save :: dt_f(:)    !< 5s-a: the local time step, mirrored by both its writers
  integer,  save :: n_cache_f = 1        !< 1 while the flat caches do not exist: the kernels skip the store

contains

  !> Linear 0-based index of a padded cell within its block, i fastest (column-major
  !> over the padded extents). i in [1-gc, ni+gc], etc. Shift to 0-based: i -> i+gc-1.
  pure integer function lin(i, j, k, ni, nj) result(c)
    integer, intent(in) :: i, j, k, ni, nj
    integer :: pi, pj
    pi = ni + 2*gc                       ! padded i-extent
    pj = nj + 2*gc                       ! padded j-extent
    c = (i + gc - 1) + (j + gc - 1)*pi + (k + gc - 1)*pi*pj
  end function lin

  !> 1-based global cell index of (b,i,j,k) in w. Single point of truth for the map.
  pure integer function cidx(b, i, j, k) result(c)
    integer, intent(in) :: b, i, j, k
    c = blk_off(b) + lin(i, j, k, bdim(1, b), bdim(2, b)) + 1
  end function cidx

  logical function field_is_allocated()
    field_is_allocated = allocated_flag
  end function field_is_allocated

  !> Return block b's flat-index parameters as plain SCALARS so a device kernel can
  !> compute the w row WITHOUT reading the module allocatables blk_off/bdim on the
  !> device (host-resident allocatables are unreliable inside `acc routine seq` ->
  !> device illegal address). The caller fetches these on the HOST before the kernel
  !> and passes them in as firstprivate scalars. boff = 0-based cell offset of block b;
  !> pi/pj = padded i/j extents (ni+2gc, nj+2gc). Row of (i,j,k) is then
  !>   boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1   (1-based, == cidx).
  subroutine field_block_index_params(b, boff, pi, pj)
    integer, intent(in)  :: b
    integer, intent(out) :: boff, pi, pj
    boff = blk_off(b)
    pi   = bdim(1, b) + 2*gc
    pj   = bdim(2, b) + 2*gc
  end subroutine field_block_index_params

  !> Allocate w + the offset/shape tables from the (already-decomposed) domain.
  !> ncell_total = sum over LOCAL blocks of (ni+2gc)(nj+2gc)(nk+2gc).
  !> Allocate the flat cache twins once the arena exists.
  subroutine field_cache_flat_alloc ()
    implicit none
    if (allocated(rho_c_f)) return
    if (.not. allocated_flag) return
    allocate(rho_c_f(ncell_total), T_c_f(ncell_total), vol_f(ncell_total), dt_f(ncell_total))
    rho_c_f = 0.0_R8; T_c_f = 0.0_R8; vol_f = 0.0_R8; dt_f = 0.0_R8
    !$acc enter data copyin(rho_c_f, T_c_f, vol_f, dt_f)
    allocate(Rg_c_f(ncell_total), mil_c_f(ncell_total)); Rg_c_f = 0.0_R8; mil_c_f = 0.0_R8
    !$acc enter data copyin(Rg_c_f, mil_c_f)
    allocate(hr_f(ncell_total)); hr_f = 0.0_R8   ! (blk%hr_c starts at 0 too)
    !$acc enter data copyin(hr_f)
    n_cache_f = ncell_total
  end subroutine field_cache_flat_alloc

  subroutine field_alloc(domain)
    type(MOSE_domain_type), intent(in) :: domain
    integer :: b, ni, nj, nk, off
    if (allocated_flag) return
    allocate(blk_off(domain%nb), bdim(3, domain%nb))
    blk_off = 0; bdim = 0
    off = 0
    do b = 1, domain%nb
      bdim(:, b) = domain%blk(b)%dim
      if (.not. is_local_block(b)) cycle
      blk_off(b) = off
      ni = bdim(1, b); nj = bdim(2, b); nk = bdim(3, b)
      off = off + (ni + 2*gc)*(nj + 2*gc)*(nk + 2*gc)
    end do
    ncell_total = off
    allocate(P_flat(ncell_total, nprim))
    allocate(PO_flat(ncell_total, nprim))
    allocate(R_flat(ncell_total, nprim))
    P_flat = 0.0_R8; PO_flat = 0.0_R8; R_flat = 0.0_R8
    ! R_flat is the device-resident var-last residual accumulator: the flux pass 2 writes it coalesced
    ! (cell-major) and Newstate reads it. ~1.3 GB at THOR scale (THOR ~19.3 GB on an A30).
    allocate(dboff(domain%nb), dpi(domain%nb), dpj(domain%nb))
    do b = 1, domain%nb
      dboff(b) = blk_off(b); dpi(b) = bdim(1,b) + 2*gc; dpj(b) = bdim(2,b) + 2*gc
    end do
    !$acc enter data create(P_flat, PO_flat, R_flat)
    !$acc enter data copyin(dboff, dpi, dpj)
    ! The arena is complete here -- P_flat/PO_flat/R_flat allocated and resident, ncell_total final -- so the
    ! flag goes up here and not at the end of the routine: field_cache_flat_alloc tests it.
    allocated_flag = .true.
    ! The flat cache twin is allocated here. Every writer of
    ! rho_c/T_c mirrors into it, so it is complete only if it exists before the first cache
    ! kernel runs; allocated at the first partition instead, it missed the initial full sweep
    ! and every cell outside the running level's list stayed at zero for the rest of the run.
    call field_cache_flat_alloc ()
    ! NB blk_off/bdim are NOT made device-resident: the device w-kernel computes the
    ! flat index from plain scalars (field_block_index_params), fetched on the host
    ! before each kernel launch, so cidx/lin run host-only (module allocatables are
    ! unreliable inside acc routine seq -> illegal address).
  end subroutine field_alloc

  !> Push the host-packed w to the device (the slice keeps w host-managed: the
  !> caller does pack_P_to_flat on the host then this H2D before the w-convective kernel).
  subroutine field_update_P_flat_device ()
    if (.not. allocated_flag) return
    !$acc update device(P_flat)
  end subroutine field_update_P_flat_device

  !> P -> P_P_flatith the device-filled ghosts: D2H the resident P, pack it on the host, H2D P_flat.
  subroutine field_sync_P_to_flat (domain)
    type(MOSE_domain_type), intent(in) :: domain
    integer :: b
    if (.not. allocated_flag) return
    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle
      !$acc update host(domain%blk(b)%P)
    end do
    call pack_P_to_flat(domain)
    !$acc update device(P_flat)
  end subroutine field_sync_P_to_flat

  !> D2H the resident R_flat for the I/O / diagnostic boundary (Write_Diagnostic); the solve loop never
  !> needs it on the host.
  subroutine field_update_R_flat_host ()
    if (.not. allocated_flag) return
    !$acc update host(R_flat)
  end subroutine field_update_R_flat_host

  !> Rebuild the var-first diagnostic residual on the host:
  !> blk%r := R_flat over the interior of every LOCAL block. R_flat holds the COMPLETE last-
  !> substage residual in both modes (device writers accumulate into it directly; host
  !> bridge sources were folded in by field_add_R_to_R_flat), so a plain assignment is
  !> correct and trivially idempotent. Interior only: the diagnostic output and the
  !> root gather read r(:,1:ni,1:nj,1:nk), and root-side REMOTE blocks are allocated
  !> interior-only.
  subroutine unpack_flat_to_R_host(domain)
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: b, i, j, k, v, ni, nj, nk
    if (.not. allocated_flag) return
    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle
      ni = bdim(1, b); nj = bdim(2, b); nk = bdim(3, b)
      do k = 1, nk
        do j = 1, nj
          do i = 1, ni
            do v = 1, nprim
              domain%blk(b)%r(v, i, j, k) = R_flat(cidx(b, i, j, k), v)
            end do
          end do
        end do
      end do
    end do
  end subroutine unpack_flat_to_R_host

  !> Fold the var-first r (holding ONLY host-source
  !> contributions: RotatingFrame/soot, synced to the device by their update-device
  !> brackets) into the var-LAST R_flat, once per RK substage AFTER all sources and BEFORE
  !> Newstate. R_flat(c,s) + r(s,i,j,k) is the same single addition Newstate's old bridge
  !> read performed -> bit-exact. Per-block kernel with ONE dummy level (r passed
  !> straight in) so the resident present-match holds (the Newstate/Copy_State pattern).
  subroutine field_add_R_to_R_flat(domain)
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: b
    if (.not. allocated_flag) return
    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle
      call R_flat_add_R_blk ( domain%blk(b)%r, b, domain%blk(b)%dim )
    end do
  end subroutine field_add_R_to_R_flat

  subroutine R_flat_add_R_blk ( r, b, n )
    integer,  intent(in) :: b, n(3)
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: r
    integer :: i, j, k, s, c0, boff, pi, pj
    call field_block_index_params(b, boff, pi, pj)
    !$acc parallel loop collapse(3) gang vector present(R_flat, r) private(c0, s)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)
      c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      do s = 1, nprim
        R_flat(c0,s) = R_flat(c0,s) + r(s,i,j,k)
      end do
    end do; end do; end do
  end subroutine R_flat_add_R_blk

  !> State hash: the sums (int64, no overflow below 2**31 elements) of the low and of the high 32 bits of the IEEE
  !> pattern of every interior w(c, v) of the local blocks -- order-independent, so the device reduction and the MPI
  !> sum give the same numbers on every rank count and every launch order; identical states <=> identical numbers.
  subroutine field_state_hash ( nb, h1, h2, ncnt )
    integer, intent(in) :: nb
    integer(8), intent(out) :: h1, h2, ncnt
    integer :: b, boff, pi, pj, n1, n2, n3, i, j, k, v, c
    integer(8) :: x, s1, s2
    h1 = 0_8; h2 = 0_8; ncnt = 0_8
    if (.not. allocated_flag) return
    do b = 1, nb
      if (.not. is_local_block(b)) cycle
      call field_block_index_params(b, boff, pi, pj)
      n1 = bdim(1, b); n2 = bdim(2, b); n3 = bdim(3, b)
      s1 = 0_8; s2 = 0_8
      !$acc parallel loop collapse(4) gang vector present(P_flat) reduction(+:s1, s2) private(c, x)
      do v = 1, nprim
        do k = 1, n3
          do j = 1, n2
            do i = 1, n1
              c = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
              x = transfer( P_flat(c, v), 0_8 )
              s1 = s1 + iand( x, 4294967295_8 )
              s2 = s2 + iand( ishft( x, -32 ), 4294967295_8 )
            end do
          end do
        end do
      end do
      h1 = h1 + s1; h2 = h2 + s2; ncnt = ncnt + int(n1, 8) * int(n2, 8) * int(n3, 8)
    end do
  end subroutine field_state_hash

  !> D2H: pull the kernel-updated P_flat back to host so unpack_flat_to_P reads it.
  subroutine field_update_P_flat_host ()
    if (.not. allocated_flag) return
    !$acc update host(P_flat)
  end subroutine field_update_P_flat_host

  !> End-of-step host syncs: P_flat/PO_flat are the device-resident truth during the step;
  !> P/PO are no longer maintained per-substage. Once per step (for Compute_Residual + IO)
  !> pull P_flat/PO_flat to host and scatter into the legacy var-first P/PO. field_update_PO_flat_host D2Hs
  !> PO_flat; unpack_flat_to_PO scatters host PO_flat -> host PO (twin of field_update_P_flat_host/unpack_flat_to_P).
  subroutine field_update_PO_flat_host ()
    if (.not. allocated_flag) return
    !$acc update host(PO_flat)
  end subroutine field_update_PO_flat_host

  subroutine unpack_flat_to_PO(domain)
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: b, i, j, k, v, ni, nj, nk
    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle
      ni = bdim(1, b); nj = bdim(2, b); nk = bdim(3, b)
      do k = 1-gc, nk+gc
        do j = 1-gc, nj+gc
          do i = 1-gc, ni+gc
            do v = 1, nprim
              domain%blk(b)%PO(v, i, j, k) = PO_flat(cidx(b, i, j, k), v)
            end do
          end do
        end do
      end do
    end do
  end subroutine unpack_flat_to_PO



  subroutine field_copy_P_flat_to_PO_flat ()
    integer :: c, v
    if (.not. allocated_flag) return
    !$acc parallel loop collapse(2) gang vector present(P_flat, PO_flat)
    do v = 1, nprim
      do c = 1, ncell_total
        PO_flat(c, v) = P_flat(c, v)
      end do
    end do
  end subroutine field_copy_P_flat_to_PO_flat

  !> Transpose the legacy var-first blk%P into the var-last P_flat (all padded cells,
  !> all local blocks). The inverse of unpack_flat_to_P.
  subroutine pack_P_to_flat(domain)
    type(MOSE_domain_type), intent(in) :: domain
    integer :: b, i, j, k, v, ni, nj, nk
    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle
      ni = bdim(1, b); nj = bdim(2, b); nk = bdim(3, b)
      do k = 1-gc, nk+gc
        do j = 1-gc, nj+gc
          do i = 1-gc, ni+gc
            do v = 1, nprim
              P_flat(cidx(b, i, j, k), v) = domain%blk(b)%P(v, i, j, k)
            end do
          end do
        end do
      end do
    end do
  end subroutine pack_P_to_flat

  !> Scatter the var-last P_flat back into the legacy var-first blk%P.
  subroutine unpack_flat_to_P(domain)
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: b, i, j, k, v, ni, nj, nk
    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle
      ni = bdim(1, b); nj = bdim(2, b); nk = bdim(3, b)
      do k = 1-gc, nk+gc
        do j = 1-gc, nj+gc
          do i = 1-gc, ni+gc
            do v = 1, nprim
              domain%blk(b)%P(v, i, j, k) = P_flat(cidx(b, i, j, k), v)
            end do
          end do
        end do
      end do
    end do
  end subroutine unpack_flat_to_P

  !> Self-test (risk register: cidx/blk_off over padded multi-block cells is the new
  !> single point of failure). Verifies unpack(pack(P)) == P bit-for-bit and that cidx
  !> is a bijection (no two cells alias, no out-of-range). Returns max|abs| round-trip
  !> deviation (must be exactly 0) and a logical ok. Call once at setup under a flag.
  subroutine field_roundtrip_check(domain, max_dev, ok)
    type(MOSE_domain_type), intent(inout) :: domain
    real(R8), intent(out) :: max_dev
    logical,  intent(out) :: ok
    integer :: b, i, j, k, v, ni, nj, nk, c
    integer, allocatable :: hits(:)
    real(R8), allocatable :: Psave(:,:,:,:)
    max_dev = 0.0_R8; ok = .true.
    ! bijection check: every w cell written exactly once across all local blocks
    allocate(hits(ncell_total)); hits = 0
    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle
      ni = bdim(1, b); nj = bdim(2, b); nk = bdim(3, b)
      do k = 1-gc, nk+gc
        do j = 1-gc, nj+gc
          do i = 1-gc, ni+gc
            c = cidx(b, i, j, k)
            if (c < 1 .or. c > ncell_total) then; ok = .false.
            else; hits(c) = hits(c) + 1; end if
          end do
        end do
      end do
    end do
    if (any(hits /= 1)) ok = .false.
    deallocate(hits)
    ! round-trip check on the first local block's P (cheap, representative)
    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle
      ni = bdim(1, b); nj = bdim(2, b); nk = bdim(3, b)
      allocate(Psave, mold=domain%blk(b)%P)
      Psave = domain%blk(b)%P
      call pack_P_to_flat(domain)
      call unpack_flat_to_P(domain)
      do k = 1-gc, nk+gc
        do j = 1-gc, nj+gc
          do i = 1-gc, ni+gc
            do v = 1, nprim
              max_dev = max(max_dev, abs(domain%blk(b)%P(v,i,j,k) - Psave(v,i,j,k)))
            end do
          end do
        end do
      end do
      deallocate(Psave)
      exit   ! one block is enough for the identity check
    end do
    if (max_dev /= 0.0_R8) ok = .false.
  end subroutine field_roundtrip_check

end module MOSE_Field_m
