!>@brief Ghost cell MPI communication for inter-block boundaries.
!> Builds a communication schedule from the BC connectivity array,
!> then provides non-blocking exchange routines for P and R arrays.
!> Messages are aggregated per remote rank to minimize latency.
!> In serial mode (USE_MPI not defined), all routines are no-ops.
module MOSE_Mod_GhostExchange
  use iso_fortran_env, only: R8 => real64, I4 => int32
  use MOSE_Mod_MPI

  implicit none
  private

! The multi-rank halo on the flat field: the device pack/unpack/Pg/bc_pg-normal kernels read and write
! P_flat (dummy `w`). WGCELL = flat row of (i,j,k) in block b (boff/pi/pj are host scalars from
! field_block_index_params, passed as args). The MPI buffers and the persistent requests are
! layout-independent; only the cell <-> buffer access is flat.
#  define WGCELL(s,i,j,k) w( boff + ((i)+gc-1) + ((j)+gc-1)*pi + ((k)+gc-1)*pi*pj + 1, s )

  !> Single ghost cell exchange entry
  type :: ghost_entry_type
    integer :: bc_idx          !< Index into domain%bc(:)
    integer :: local_block     !< Block on this rank
    integer :: remote_block    !< Block on the remote rank
    integer :: remote_rank     !< MPI rank of the remote block
    integer :: src_i, src_j, src_k, src_f  !< Source cell coords and face
    integer :: dst_i, dst_j, dst_k, dst_f  !< Destination cell coords and face
  end type ghost_entry_type

  !> Per-rank aggregation info
  type :: rank_group_type
    integer :: rank = -1       !< Remote MPI rank
    integer :: offset = 0      !< Start index in entry list (1-based)
    integer :: count  = 0      !< Number of entries for this rank
  end type rank_group_type

  !> Face group: contiguous run of entries sharing (block, face).
  !> Precomputes guide offsets; for rectangular faces enables direct coordinate sweep.
  type :: face_group_type
    integer :: block            !< Local block index
    integer :: face             !< Face number (1-6)
    integer :: entry_start      !< First entry index in send/recv list
    integer :: entry_count      !< Number of entries
    integer :: di, dj, dk       !< Precomputed guide(face,:)
  end type face_group_type

  !> Communication schedule: entries grouped by remote rank for aggregated messaging
  type :: ghost_schedule_type
    integer :: n_send = 0, n_recv = 0
    type(ghost_entry_type), allocatable :: send_list(:)
    type(ghost_entry_type), allocatable :: recv_list(:)
    logical :: built = .false.

    ! Per-rank aggregation for sends
    integer :: n_send_ranks = 0
    type(rank_group_type), allocatable :: send_groups(:)

    ! Per-rank aggregation for receives
    integer :: n_recv_ranks = 0
    type(rank_group_type), allocatable :: recv_groups(:)

    ! Face groups for aggregated pack/unpack (reduce per-cell overhead)
    integer :: n_send_faces = 0
    type(face_group_type), allocatable :: send_faces(:)
    integer :: n_recv_faces = 0
    type(face_group_type), allocatable :: recv_faces(:)

    ! Pre-allocated MPI buffers (one contiguous buffer per rank group).
    ! P_send_buf / P_recv_buf are made DEVICE-RESIDENT (acc enter data
    ! create) so the P halo is packed/unpacked by device kernels and (with CUDA-aware
    ! MPI) handed straight to MPI via host_data use_device. The persistent P requests
    ! bind to these buffers ONCE at init, so the device pointer must be stable for the
    ! whole run -> resident-for-life, freed only at cleanup. Pg_* stay host (the
    ! transverse-stencil exchange is a deferred follow-on increment).
    real(R8), allocatable :: P_send_buf(:)   !< (nprim*gc * n_send) -- device-resident
    real(R8), allocatable :: P_recv_buf(:)   !< (nprim*gc * n_recv) -- device-resident
    ! The transverse Pg(:,3:6) stencil is now ALSO produced
    ! and exchanged ON DEVICE (donor-side gather from resident donor P -> device pack ->
    ! GPU-aware MPI -> device unpack into the resident bc_pg). Pg_send/recv_buf are made
    ! device-resident (acc enter data create) exactly like P_send/recv_buf, and the
    ! persistent Pg requests bind to the DEVICE buffer (host_data use_device) when CUDA-aware
    ! MPI + the device halo path are live. SEPARATE buffers from P (the P exchange can still
    ! be in flight while the Pg gather runs -> no aliasing). Entry i maps to slot
    ! (i-1)*nprim*PG_NCOLS, mirroring the host exchange_Pg_field convention byte-for-byte.
    real(R8), allocatable :: Pg_send_buf(:)  !< (nprim*4  * n_send) -- device-resident
    real(R8), allocatable :: Pg_recv_buf(:)  !< (nprim*4  * n_recv) -- device-resident
    integer, allocatable  :: send_req(:)     !< (n_send_ranks)
    integer, allocatable  :: recv_req(:)     !< (n_recv_ranks)
    integer, allocatable  :: mpi_stat(:,:)   !< (MPI_STATUS_SIZE, max(n_send_ranks, n_recv_ranks))
    ! Persistent MPI requests for P exchange (initialized once, started each step)
    integer, allocatable  :: P_send_req_pers(:)  !< (n_send_ranks)
    integer, allocatable  :: P_recv_req_pers(:)  !< (n_recv_ranks)
    integer, allocatable  :: R_send_req_pers(:)  !< (n_send_ranks)
    integer, allocatable  :: R_recv_req_pers(:)  !< (n_recv_ranks)
    ! Persistent MPI requests for the device-resident Pg exchange.
    integer, allocatable  :: Pg_send_req_pers(:) !< (n_send_ranks)
    integer, allocatable  :: Pg_recv_req_pers(:) !< (n_recv_ranks)
    logical :: persistent_P_init = .false.
    logical :: persistent_R_init = .false.
    logical :: persistent_Pg_init = .false.
    logical :: pg_bufs_resident = .false.

    ! Flat per-ENTRY index tables (plain integers) mirroring send_list /
    ! recv_list in their EXACT rank-sorted order, so device pack/unpack kernels can read
    ! them (the ghost_entry_type AoS arrays do not map to the device; same flatten as the
    ! conn_* tables, ). Buffer slot for entry i is (i-1)*nprim*gc, identical to
    ! the host pack/unpack convention, so the persistent-request rank-group offsets are
    ! unchanged. Grouped by LOCAL block via [blk_lo,blk_hi] so one kernel per block passes
    ! a single present-matched blk(b)%P ( cross-block trick).
    integer, allocatable  :: send_b(:)                 !< local block of send entry i
    integer, allocatable  :: send_is(:), send_js(:), send_ks(:)  !< donor interior cell
    integer, allocatable  :: send_di(:), send_dj(:), send_dk(:)  !< guide(src_f,:)
    integer, allocatable  :: recv_b(:)                 !< local block of recv entry i
    integer, allocatable  :: recv_id(:), recv_jd(:), recv_kd(:)  !< receiver boundary cell
    integer, allocatable  :: recv_di(:), recv_dj(:), recv_dk(:)  !< guide(dst_f,:)
    ! Per-entry tables for the device Pg transverse exchange.
    ! SEND side (donor local): the donor face + the 4 transverse offset pairs, used by the
    ! device pack kernel to gather the 4 stencil cells from the resident donor blk(Bs)%P
    ! (mirrors host Fill_BC_Ghost_Connection). The donor block + donor cell reuse send_b /
    ! send_is/js/ks (already built). RECV side: pg_recv_slot = the flat bc_pg slot
    ! (= bc_perm inverse) the received transverse stencil unpacks into, AND where the
    ! receiver-normal kernel writes bc_pg(:,1:gc) from the device-unpacked receiver ghost P.
    ! pg_recv_slot needs bc_perm, which Setup_BC_Device builds AFTER build_ghost_schedule ->
    ! it (and the resident-make) is filled in finalize_pg_tables, called from there.
    integer, allocatable  :: pg_send_fs(:)             !< donor face Fs of send entry i
    integer, allocatable  :: pg_send_d11(:), pg_send_d12(:), pg_send_d21(:), pg_send_d22(:)
    integer, allocatable  :: pg_recv_slot(:)           !< bc_pg flat slot for recv entry i
    logical :: pg_tables_built = .false.
    logical :: dev_tables_built = .false.
    logical :: bufs_resident = .false.
    ! True iff the device-resident inter-rank halo driver
    ! (Fill_Ghost_Cell_mr) is the active ghost-fill path for this domain, i.e.
    ! every LOCAL boundary type is device-fillable {1,3,4,11,22} (mirrors the
    ! bc_ghost_dev_mr_ok gate, computed here because the persistent-request binding
    ! happens before Setup_BC_Device runs). When .false. the host Fill_Ghost_Cell
    ! path runs and packs/unpacks the HOST P_send/recv_buf, so the persistent MPI
    ! requests MUST stay host-bound -- binding them to the device buffer would make
    ! MPI DMA device memory the host path never touches (-> CUDA illegal address /
    ! coherence break on the inter-rank seam; the  2-rank crash).
    ! Whether the P and Pg persistent requests bind to the DEVICE buffers (GPU-aware MPI; UCX moves them GPU to GPU
    ! with CUDA-IPC) or to the HOST buffers (the halo staged through host memory): follows gpu_aware_mpi_runtime.
    logical :: dev_bound = .true.
  end type ghost_schedule_type

  ! The schedule of all the inter-rank connection ghosts
  type(ghost_schedule_type), public :: ghost_sched

  !> CUDA-aware (GPU-aware) MPI on the P halo: the persistent requests bind to the device buffers.
  !> MOSE_GPU_AWARE_MPI=0 in the environment forces the host-staging fallback (update host -> MPI ->
  !> update device on the slab buffers) on a machine whose MPI is not CUDA-aware, with no rebuild.
  logical, public :: gpu_aware_mpi_runtime = .true.

  public :: build_ghost_schedule
  public :: build_local_bc_index
  public :: cleanup_ghost_schedule
  public :: exchange_ghost_P_post_recv
  public :: exchange_ghost_P_post_send
  public :: exchange_ghost_P_wait_send
  ! GPU Phase 1 device-halo entry points (no-ops / not referenced in the CPU build)
  public :: exchange_ghost_P_pack, exchange_ghost_P_unpack
  public :: exchange_ghost_P_wait_recv_only

  ! Device transverse-Pg pack/exchange/unpack entry points.
  public :: finalize_pg_tables
  public :: exchange_ghost_Pg_pack, exchange_ghost_Pg_post_recv, exchange_ghost_Pg_post_send
  public :: exchange_ghost_Pg_wait_recv, exchange_ghost_Pg_unpack, exchange_ghost_Pg_wait_send
  public :: fill_bc_pg_normal
  public :: gather_P_to_root
  public :: gather_diagnostic_to_root
  public :: mpi_io_barrier


contains


  !> Build the communication schedule by scanning all BC entries of type 1 (connection).
  !> Entries are sorted by remote rank for aggregated MPI messaging.
  !> Must be called after partition_blocks and domain setup.
  subroutine build_ghost_schedule(domain)
    use MOSE_Types_m
    use MOSE_Global_m, only: gc, nprim

    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    ! Local
    integer :: i, ns, nr, bm, bs

    if (mpi_size_ <= 1) then
      ghost_sched%built = .true.
      return
    end if

    ! Runtime override for the CUDA-aware-MPI device-pointer binding (host-staging fallback
    ! for machines without CUDA-aware MPI, and for A/B isolation of the device-bound MPI seam
    ! vs the device pack/unpack logic). Set MOSE_GPU_AWARE_MPI=0 to force host-staging.
    block
      character(len=8) :: envval
      integer :: envlen, envstat
      call get_environment_variable('MOSE_GPU_AWARE_MPI', envval, envlen, envstat)
      if (envstat == 0 .and. envlen > 0) then
        if (envval(1:1) == '0' .or. envval(1:1) == 'n' .or. envval(1:1) == 'N' .or. &
            envval(1:1) == 'f' .or. envval(1:1) == 'F') then
          gpu_aware_mpi_runtime = .false.
          if (mpi_is_root) write(*,'(A)') ' GhostExchange: MOSE_GPU_AWARE_MPI=0 -> host-staging P halo'
        end if
      end if
    end block
    ghost_sched%dev_bound = gpu_aware_mpi_runtime

#ifdef USE_MPI
    ! Count send and recv entries
    ns = 0; nr = 0
    do i = 1, domain%nbound
      if (domain%bc(i)%type /= 1) cycle
      bm = domain%bc(i)%b   ! destination block
      bs = domain%bc(i)%bs  ! source block
      if (.not. allocated(block_owner)) cycle
      ! If I own the destination but not the source => I need to RECV
      if (is_local_block(bm) .and. (.not. is_local_block(bs))) then
        nr = nr + 1
      end if
      ! If I own the source but not the destination => I need to SEND
      if (is_local_block(bs) .and. (.not. is_local_block(bm))) then
        ns = ns + 1
      end if
    end do

    ghost_sched%n_send = ns
    ghost_sched%n_recv = nr
    if (allocated(ghost_sched%send_list)) deallocate(ghost_sched%send_list)
    if (allocated(ghost_sched%recv_list)) deallocate(ghost_sched%recv_list)
    allocate(ghost_sched%send_list(ns))
    allocate(ghost_sched%recv_list(nr))

    ! Fill send and recv lists
    ns = 0; nr = 0
    do i = 1, domain%nbound
      if (domain%bc(i)%type /= 1) cycle
      bm = domain%bc(i)%b
      bs = domain%bc(i)%bs

      if (is_local_block(bm) .and. (.not. is_local_block(bs))) then
        nr = nr + 1
        ghost_sched%recv_list(nr)%bc_idx = i
        ghost_sched%recv_list(nr)%local_block  = bm
        ghost_sched%recv_list(nr)%remote_block = bs
        ghost_sched%recv_list(nr)%remote_rank  = block_owner(bs)
        ghost_sched%recv_list(nr)%src_i = domain%bc(i)%is
        ghost_sched%recv_list(nr)%src_j = domain%bc(i)%js
        ghost_sched%recv_list(nr)%src_k = domain%bc(i)%ks
        ghost_sched%recv_list(nr)%src_f = domain%bc(i)%fs
        ghost_sched%recv_list(nr)%dst_i = domain%bc(i)%i
        ghost_sched%recv_list(nr)%dst_j = domain%bc(i)%j
        ghost_sched%recv_list(nr)%dst_k = domain%bc(i)%k
        ghost_sched%recv_list(nr)%dst_f = domain%bc(i)%f
      end if

      if (is_local_block(bs) .and. (.not. is_local_block(bm))) then
        ns = ns + 1
        ghost_sched%send_list(ns)%bc_idx = i
        ghost_sched%send_list(ns)%local_block  = bs
        ghost_sched%send_list(ns)%remote_block = bm
        ghost_sched%send_list(ns)%remote_rank  = block_owner(bm)
        ghost_sched%send_list(ns)%src_i = domain%bc(i)%is
        ghost_sched%send_list(ns)%src_j = domain%bc(i)%js
        ghost_sched%send_list(ns)%src_k = domain%bc(i)%ks
        ghost_sched%send_list(ns)%src_f = domain%bc(i)%fs
        ghost_sched%send_list(ns)%dst_i = domain%bc(i)%i
        ghost_sched%send_list(ns)%dst_j = domain%bc(i)%j
        ghost_sched%send_list(ns)%dst_k = domain%bc(i)%k
        ghost_sched%send_list(ns)%dst_f = domain%bc(i)%f
      end if
    end do

    ! Sort entries by remote rank and build rank groups
    call sort_entries_by_rank(ghost_sched%send_list, ns)
    call sort_entries_by_rank(ghost_sched%recv_list, nr)
    call build_rank_groups(ghost_sched%send_list, ns, ghost_sched%send_groups, ghost_sched%n_send_ranks)
    call build_rank_groups(ghost_sched%recv_list, nr, ghost_sched%recv_groups, ghost_sched%n_recv_ranks)

    ! Build face groups for aggregated pack/unpack
    call build_face_groups_send(ghost_sched%send_list, ns, &
                                ghost_sched%send_faces, ghost_sched%n_send_faces)
    call build_face_groups_recv(ghost_sched%recv_list, nr, &
                                ghost_sched%recv_faces, ghost_sched%n_recv_faces)

    ! Pre-allocate MPI buffers (reused every RK step)
    ns = ghost_sched%n_send
    nr = ghost_sched%n_recv
    if (ns > 0 .or. nr > 0) then
      block
        use mpi, only: MPI_STATUS_SIZE
        integer :: max_ranks
        max_ranks = max(ghost_sched%n_send_ranks, ghost_sched%n_recv_ranks, 1)
        allocate(ghost_sched%P_send_buf(nprim*gc * max(ns,1)))
        allocate(ghost_sched%P_recv_buf(nprim*gc * max(nr,1)))
        allocate(ghost_sched%Pg_send_buf(nprim*4 * max(ns,1)))
        allocate(ghost_sched%Pg_recv_buf(nprim*4 * max(nr,1)))
        allocate(ghost_sched%send_req(max(ghost_sched%n_send_ranks,1)))
        allocate(ghost_sched%recv_req(max(ghost_sched%n_recv_ranks,1)))
        allocate(ghost_sched%mpi_stat(MPI_STATUS_SIZE, max_ranks))
      end block
    end if

    ! Decide whether the device-resident halo driver (Fill_Ghost_Cell_mr)
    ! will be the active ghost-fill path. It runs ONLY when every LOCAL boundary type is
    ! device-fillable {1,3,4,11,22} -- exactly the bc_ghost_dev_mr_ok criterion, recomputed
    ! here because Setup_BC_Device (which sets that flag) runs AFTER this routine, yet the
    ! persistent-request device binding below must already know the answer. When a local
    ! type falls outside that set (e.g. type-2 rotational periodicity in SpicchiOThor) the
    ! host Fill_Ghost_Cell path runs instead and packs/unpacks the HOST P buffers, so the
    ! persistent MPI requests must stay HOST-bound (binding them to the resident device
    ! buffer made MPI DMA device memory the host path never touches -> the  inter-rank
    ! CUDA_ERROR_ILLEGAL_ADDRESS). Mirror the exhaustive case list of Setup_BC_Device.
    ! NB this list and Setup_BC_Device's bc_ghost_dev_mr_ok MUST agree: if this says
    ! "host" while that says "device", Fill_Ghost_Cell_mr packs on the device and
    ! post_send(dev_packed) skips the host bridge, so host-bound persistent requests ship
    ! a never-packed host buffer -> garbage ghosts on the receivers (chimera, 2026-08-25).

    ! Build the flat per-entry index tables and make the P halo buffers
    ! device-resident. Both MUST precede init_persistent_P_requests: that routine binds
    ! the persistent requests to the resident device buffer pointer when CUDA-aware MPI
    ! is active AND the device halo path is live (see init_persistent_P_requests), and the
    ! pack/unpack kernels read the flat tables.
    call build_dev_index_tables(domain)
    call make_P_buffers_resident()
    call make_Pg_buffers_resident()

    ! Initialize persistent MPI requests for P and R exchanges
    call init_persistent_P_requests()
    call init_persistent_R_requests()

    if (mpi_is_root) then
      write(*,'(A,I0,A,I0,A)') ' Ghost schedule: ', ghost_sched%n_send_ranks, &
        ' send ranks, ', ghost_sched%n_recv_ranks, ' recv ranks'
    end if
#endif

    ghost_sched%built = .true.
  end subroutine build_ghost_schedule


  !> Build pre-filtered list of BC indices where Bm is local to this rank.
  !> Avoids O(nbound) loop scanning in hot paths (Fill_Ghost_Cell, etc.).
  subroutine build_local_bc_index(domain)
    use MOSE_Types_m

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: i, n, ns

    ! Count local BC entries (Bm is local) and local source entries (Bs is local, type-1)
    n = 0; ns = 0
    do i = 1, domain%nbound
      if (is_local_block(domain%bc(i)%b)) n = n + 1
      if (domain%bc(i)%type == 1) then
        if (is_local_block(domain%bc(i)%bs)) ns = ns + 1
      end if
    end do

    ! Build local_bc_idx (Bm is local)
    domain%n_local_bc = n
    if (allocated(domain%local_bc_idx)) deallocate(domain%local_bc_idx)
    allocate(domain%local_bc_idx(n))

    ! Build local_bs_idx (Bs is local, type-1 only — for Pg computation)
    domain%n_local_bs = ns
    if (allocated(domain%local_bs_idx)) deallocate(domain%local_bs_idx)
    allocate(domain%local_bs_idx(ns))

    ! Fill both index arrays in a single pass
    n = 0; ns = 0
    do i = 1, domain%nbound
      if (is_local_block(domain%bc(i)%b)) then
        n = n + 1
        domain%local_bc_idx(n) = i
      end if
      if (domain%bc(i)%type == 1) then
        if (is_local_block(domain%bc(i)%bs)) then
          ns = ns + 1
          domain%local_bs_idx(ns) = i
        end if
      end if
    end do
  end subroutine build_local_bc_index


  !> Free persistent MPI requests and deallocate schedule arrays.
  subroutine cleanup_ghost_schedule()
    implicit none
#ifdef USE_MPI
    call cleanup_persistent_requests()
    ! Release the device-resident P halo buffers. Persistent requests are
    ! freed first (above) so MPI no longer references the device pointer.
    if (ghost_sched%bufs_resident) then
      if (allocated(ghost_sched%P_send_buf)) then
        !$acc exit data delete(ghost_sched%P_send_buf)
      end if
      if (allocated(ghost_sched%P_recv_buf)) then
        !$acc exit data delete(ghost_sched%P_recv_buf)
      end if
      ghost_sched%bufs_resident = .false.
    end if
    ! : release the device-resident transverse-Pg buffers (after the Pg requests freed).
    if (ghost_sched%pg_bufs_resident) then
      if (allocated(ghost_sched%Pg_send_buf)) then
        !$acc exit data delete(ghost_sched%Pg_send_buf)
      end if
      if (allocated(ghost_sched%Pg_recv_buf)) then
        !$acc exit data delete(ghost_sched%Pg_recv_buf)
      end if
      ghost_sched%pg_bufs_resident = .false.
    end if
#endif
  end subroutine cleanup_ghost_schedule


#ifdef USE_MPI
  !> Build flat per-entry index tables that mirror send_list/recv_list in
  !> their exact (rank-sorted) order so the device pack/unpack kernels can read plain
  !> integer arrays instead of the unmappable ghost_entry_type AoS. Entry i still maps to
  !> buffer slot (i-1)*nprim*gc, so the persistent-request rank-group offsets are unchanged.
  subroutine build_dev_index_tables(domain)
    use MOSE_Types_m
    use MOSE_Global_m, only: guide
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer :: i, ns, nr, bidx

    ns = ghost_sched%n_send
    nr = ghost_sched%n_recv

    if (allocated(ghost_sched%send_b)) deallocate(ghost_sched%send_b, &
        ghost_sched%send_is, ghost_sched%send_js, ghost_sched%send_ks, &
        ghost_sched%send_di, ghost_sched%send_dj, ghost_sched%send_dk)
    if (allocated(ghost_sched%recv_b)) deallocate(ghost_sched%recv_b, &
        ghost_sched%recv_id, ghost_sched%recv_jd, ghost_sched%recv_kd, &
        ghost_sched%recv_di, ghost_sched%recv_dj, ghost_sched%recv_dk)
    if (allocated(ghost_sched%pg_send_fs)) deallocate(ghost_sched%pg_send_fs, &
        ghost_sched%pg_send_d11, ghost_sched%pg_send_d12, ghost_sched%pg_send_d21, ghost_sched%pg_send_d22)

    allocate(ghost_sched%send_b(max(ns,1)))
    allocate(ghost_sched%send_is(max(ns,1)), ghost_sched%send_js(max(ns,1)), ghost_sched%send_ks(max(ns,1)))
    allocate(ghost_sched%send_di(max(ns,1)), ghost_sched%send_dj(max(ns,1)), ghost_sched%send_dk(max(ns,1)))
    allocate(ghost_sched%recv_b(max(nr,1)))
    allocate(ghost_sched%recv_id(max(nr,1)), ghost_sched%recv_jd(max(nr,1)), ghost_sched%recv_kd(max(nr,1)))
    allocate(ghost_sched%recv_di(max(nr,1)), ghost_sched%recv_dj(max(nr,1)), ghost_sched%recv_dk(max(nr,1)))
    allocate(ghost_sched%pg_send_fs(max(ns,1)))
    allocate(ghost_sched%pg_send_d11(max(ns,1)), ghost_sched%pg_send_d12(max(ns,1)), &
             ghost_sched%pg_send_d21(max(ns,1)), ghost_sched%pg_send_d22(max(ns,1)))

    do i = 1, ns
      ghost_sched%send_b(i)  = ghost_sched%send_list(i)%local_block
      ghost_sched%send_is(i) = ghost_sched%send_list(i)%src_i
      ghost_sched%send_js(i) = ghost_sched%send_list(i)%src_j
      ghost_sched%send_ks(i) = ghost_sched%send_list(i)%src_k
      ghost_sched%send_di(i) = guide(ghost_sched%send_list(i)%src_f, 1)
      ghost_sched%send_dj(i) = guide(ghost_sched%send_list(i)%src_f, 2)
      ghost_sched%send_dk(i) = guide(ghost_sched%send_list(i)%src_f, 3)
      !  device Pg gather: donor face + transverse offset pairs from the BC record.
      bidx = ghost_sched%send_list(i)%bc_idx
      ghost_sched%pg_send_fs(i)  = domain%bc(bidx)%fs
      ghost_sched%pg_send_d11(i) = domain%bc(bidx)%d11
      ghost_sched%pg_send_d12(i) = domain%bc(bidx)%d12
      ghost_sched%pg_send_d21(i) = domain%bc(bidx)%d21
      ghost_sched%pg_send_d22(i) = domain%bc(bidx)%d22
    end do
    do i = 1, nr
      ghost_sched%recv_b(i)  = ghost_sched%recv_list(i)%local_block
      ghost_sched%recv_id(i) = ghost_sched%recv_list(i)%dst_i
      ghost_sched%recv_jd(i) = ghost_sched%recv_list(i)%dst_j
      ghost_sched%recv_kd(i) = ghost_sched%recv_list(i)%dst_k
      ghost_sched%recv_di(i) = guide(ghost_sched%recv_list(i)%dst_f, 1)
      ghost_sched%recv_dj(i) = guide(ghost_sched%recv_list(i)%dst_f, 2)
      ghost_sched%recv_dk(i) = guide(ghost_sched%recv_list(i)%dst_f, 3)
    end do

    ! Resident-for-life: read by the pack/unpack kernels every substage.
    !$acc enter data copyin(ghost_sched%send_b, ghost_sched%send_is, ghost_sched%send_js, &
    !$acc&  ghost_sched%send_ks, ghost_sched%send_di, ghost_sched%send_dj, ghost_sched%send_dk, &
    !$acc&  ghost_sched%recv_b, ghost_sched%recv_id, ghost_sched%recv_jd, ghost_sched%recv_kd, &
    !$acc&  ghost_sched%recv_di, ghost_sched%recv_dj, ghost_sched%recv_dk, &
    !$acc&  ghost_sched%pg_send_fs, ghost_sched%pg_send_d11, ghost_sched%pg_send_d12, &
    !$acc&  ghost_sched%pg_send_d21, ghost_sched%pg_send_d22)
    ghost_sched%dev_tables_built = .true.
  end subroutine build_dev_index_tables


  !> Make the P send/recv halo buffers device-resident for the whole run.
  !> Done once, before the persistent requests bind to them. No-op in the CPU build.
  subroutine make_P_buffers_resident()
    implicit none
    if (.not. allocated(ghost_sched%P_send_buf)) return
    !$acc enter data create(ghost_sched%P_send_buf, ghost_sched%P_recv_buf)
    ghost_sched%bufs_resident = .true.
  end subroutine make_P_buffers_resident


  !> Make the transverse-Pg send/recv buffers device-resident
  !> for the whole run, exactly like make_P_buffers_resident. The persistent Pg requests bind
  !> to these device pointers, so they must be stable for life. No-op in the CPU build.
  subroutine make_Pg_buffers_resident()
    implicit none
    if (.not. allocated(ghost_sched%Pg_send_buf)) return
    !$acc enter data create(ghost_sched%Pg_send_buf, ghost_sched%Pg_recv_buf)
    ghost_sched%pg_bufs_resident = .true.
  end subroutine make_Pg_buffers_resident
#endif


  !> Build the per-recv-entry bc_pg slot table and initialize
  !> the persistent Pg MPI requests. Must run AFTER Setup_BC_Device has built domain%bc_perm
  !> (the file-order -> sorted-slot permutation that bc_pg is indexed by), which happens after
  !> build_ghost_schedule. Called from Setup_BC_Device. The send-side Pg tables (pg_send_*)
  !> were already built + made resident in build_dev_index_tables (they need no bc_perm).
  !> Public + depth-0 (body USE_MPI-guarded) so the name resolves in the non-MPI build.
  subroutine finalize_pg_tables(domain)
    use MOSE_Types_m
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use MOSE_Global_m, only: nprim
      use mpi
      integer :: i, nr, bidx, slot
      integer :: r, ierr, tag, buf_pos, entry_size
      integer, parameter :: PG_NCOLS = 4
      integer, allocatable :: slot_of(:)

      if (.not. ghost_sched%built) return

      nr = ghost_sched%n_recv

      ! Inverse permutation: bc file index -> sorted bc_pg slot (same as Setup_BC_Device).
      allocate(slot_of(domain%nbound))
      do i = 1, domain%nbound
        slot_of(domain%bc_perm(i)) = i
      end do

      if (allocated(ghost_sched%pg_recv_slot)) deallocate(ghost_sched%pg_recv_slot)
      allocate(ghost_sched%pg_recv_slot(max(nr,1)))
      ghost_sched%pg_recv_slot = 0
      do i = 1, nr
        bidx = ghost_sched%recv_list(i)%bc_idx
        slot = slot_of(bidx)
        ghost_sched%pg_recv_slot(i) = slot
      end do
      deallocate(slot_of)

      !$acc enter data copyin(ghost_sched%pg_recv_slot)
      ghost_sched%pg_tables_built = .true.

      ! ---- persistent Pg requests (mirror init_persistent_P_requests, distinct tags) ----
      if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return
      entry_size = nprim * PG_NCOLS
      if (allocated(ghost_sched%Pg_recv_req_pers)) return  ! already initialised
      allocate(ghost_sched%Pg_recv_req_pers(max(ghost_sched%n_recv_ranks, 1)))
      allocate(ghost_sched%Pg_send_req_pers(max(ghost_sched%n_send_ranks, 1)))

      ! Device-pointer bound when CUDA-aware MPI + the device halo path are live (same gate as
      ! the P requests): MPI_STARTALL then DMAs straight out of / into the resident Pg buffers.
      if (ghost_sched%dev_bound .and. ghost_sched%pg_bufs_resident) then
        do r = 1, ghost_sched%n_recv_ranks
          buf_pos = (ghost_sched%recv_groups(r)%offset - 1) * entry_size + 1
          tag = ghost_sched%recv_groups(r)%rank + 20000
          !$acc host_data use_device(ghost_sched%Pg_recv_buf)
          call MPI_RECV_INIT(ghost_sched%Pg_recv_buf(buf_pos), &
                             ghost_sched%recv_groups(r)%count * entry_size, &
                             MPI_DOUBLE_PRECISION, ghost_sched%recv_groups(r)%rank, tag, &
                             MPI_COMM_WORLD, ghost_sched%Pg_recv_req_pers(r), ierr)
          !$acc end host_data
          call check_mpi_error(ierr)
        end do
        do r = 1, ghost_sched%n_send_ranks
          buf_pos = (ghost_sched%send_groups(r)%offset - 1) * entry_size + 1
          tag = mpi_rank_ + 20000
          !$acc host_data use_device(ghost_sched%Pg_send_buf)
          call MPI_SEND_INIT(ghost_sched%Pg_send_buf(buf_pos), &
                             ghost_sched%send_groups(r)%count * entry_size, &
                             MPI_DOUBLE_PRECISION, ghost_sched%send_groups(r)%rank, tag, &
                             MPI_COMM_WORLD, ghost_sched%Pg_send_req_pers(r), ierr)
          !$acc end host_data
          call check_mpi_error(ierr)
        end do
        ghost_sched%persistent_Pg_init = .true.
        return
      end if

      ! Host-bound persistent Pg requests (CPU build / GPU host-staging fallback).
      do r = 1, ghost_sched%n_recv_ranks
        buf_pos = (ghost_sched%recv_groups(r)%offset - 1) * entry_size + 1
        tag = ghost_sched%recv_groups(r)%rank + 20000
        call MPI_RECV_INIT(ghost_sched%Pg_recv_buf(buf_pos), &
                           ghost_sched%recv_groups(r)%count * entry_size, &
                           MPI_DOUBLE_PRECISION, ghost_sched%recv_groups(r)%rank, tag, &
                           MPI_COMM_WORLD, ghost_sched%Pg_recv_req_pers(r), ierr)
        call check_mpi_error(ierr)
      end do
      do r = 1, ghost_sched%n_send_ranks
        buf_pos = (ghost_sched%send_groups(r)%offset - 1) * entry_size + 1
        tag = mpi_rank_ + 20000
        call MPI_SEND_INIT(ghost_sched%Pg_send_buf(buf_pos), &
                           ghost_sched%send_groups(r)%count * entry_size, &
                           MPI_DOUBLE_PRECISION, ghost_sched%send_groups(r)%rank, tag, &
                           MPI_COMM_WORLD, ghost_sched%Pg_send_req_pers(r), ierr)
        call check_mpi_error(ierr)
      end do
      ghost_sched%persistent_Pg_init = .true.
    end block
#endif
  end subroutine finalize_pg_tables


#ifdef USE_MPI
  !> Initialize persistent MPI requests for P field exchange.
  !> Called once from build_ghost_schedule. Uses MPI_SEND_INIT/RECV_INIT.
  subroutine init_persistent_P_requests()
    use MOSE_Global_m, only: gc, nprim
    use mpi

    implicit none
    integer :: r, ierr, tag, buf_pos, entry_size

    entry_size = nprim * gc

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return

    allocate(ghost_sched%P_recv_req_pers(max(ghost_sched%n_recv_ranks, 1)))
    allocate(ghost_sched%P_send_req_pers(max(ghost_sched%n_send_ranks, 1)))

    ! Bind the persistent P requests to the DEVICE copy of the buffers when
    ! CUDA-aware MPI is in use, so MPI_STARTALL later DMA's straight out of / into device
    ! memory (no host bounce). host_data use_device must cover the *_INIT call because a
    ! persistent request captures the buffer address once and reuses it on every START.
    ! The buffers were made resident just before this call (build_ghost_schedule), so the
    ! device address exists. Outside the GPU-aware path (CPU build, or fallback toggle)
    ! the requests bind to the HOST buffer exactly as before; the exchange then host-stages.
    ! dev_halo_active gates this too: only the device halo driver (Fill_Ghost_Cell_mr)
    ! packs/unpacks the DEVICE buffers, so device-bound requests are valid only then. If the
    ! host Fill_Ghost_Cell path is live (a non-device-fillable local type, e.g. type-2),
    ! it touches the HOST buffers -> the requests MUST stay host-bound or MPI faults DMAing
    ! device memory ( inter-rank CUDA illegal address).
    if (ghost_sched%dev_bound .and. ghost_sched%bufs_resident) then
      ! Persistent receives (device-pointer bound)
      do r = 1, ghost_sched%n_recv_ranks
        buf_pos = (ghost_sched%recv_groups(r)%offset - 1) * entry_size + 1
        tag = ghost_sched%recv_groups(r)%rank
        !$acc host_data use_device(ghost_sched%P_recv_buf)
        call MPI_RECV_INIT(ghost_sched%P_recv_buf(buf_pos), &
                           ghost_sched%recv_groups(r)%count * entry_size, &
                           MPI_DOUBLE_PRECISION, &
                           ghost_sched%recv_groups(r)%rank, tag, &
                           MPI_COMM_WORLD, ghost_sched%P_recv_req_pers(r), ierr)
        !$acc end host_data
        call check_mpi_error(ierr)
      end do
      ! Persistent sends (device-pointer bound)
      do r = 1, ghost_sched%n_send_ranks
        buf_pos = (ghost_sched%send_groups(r)%offset - 1) * entry_size + 1
        tag = mpi_rank_
        !$acc host_data use_device(ghost_sched%P_send_buf)
        call MPI_SEND_INIT(ghost_sched%P_send_buf(buf_pos), &
                           ghost_sched%send_groups(r)%count * entry_size, &
                           MPI_DOUBLE_PRECISION, &
                           ghost_sched%send_groups(r)%rank, tag, &
                           MPI_COMM_WORLD, ghost_sched%P_send_req_pers(r), ierr)
        !$acc end host_data
        call check_mpi_error(ierr)
      end do
      ghost_sched%persistent_P_init = .true.
      return
    end if

    ! Host-bound persistent requests (CPU build, or GPU host-staging fallback).
    ! Initialize persistent receives
    do r = 1, ghost_sched%n_recv_ranks
      buf_pos = (ghost_sched%recv_groups(r)%offset - 1) * entry_size + 1
      tag = ghost_sched%recv_groups(r)%rank
      call MPI_RECV_INIT(ghost_sched%P_recv_buf(buf_pos), &
                         ghost_sched%recv_groups(r)%count * entry_size, &
                         MPI_DOUBLE_PRECISION, &
                         ghost_sched%recv_groups(r)%rank, tag, &
                         MPI_COMM_WORLD, ghost_sched%P_recv_req_pers(r), ierr)
      call check_mpi_error(ierr)
    end do

    ! Initialize persistent sends
    do r = 1, ghost_sched%n_send_ranks
      buf_pos = (ghost_sched%send_groups(r)%offset - 1) * entry_size + 1
      tag = mpi_rank_
      call MPI_SEND_INIT(ghost_sched%P_send_buf(buf_pos), &
                         ghost_sched%send_groups(r)%count * entry_size, &
                         MPI_DOUBLE_PRECISION, &
                         ghost_sched%send_groups(r)%rank, tag, &
                         MPI_COMM_WORLD, ghost_sched%P_send_req_pers(r), ierr)
      call check_mpi_error(ierr)
    end do

    ghost_sched%persistent_P_init = .true.
  end subroutine init_persistent_P_requests


  !> Initialize persistent MPI requests for R field exchange.
  !> Reuses P_send_buf/P_recv_buf (same size). Uses distinct tags (offset 30000).
  !>
  !> ⚠️ R/P BUFFER ALIASING (GPU Phase 1 watch-item): these R requests are ALWAYS
  !> HOST-bound, while on the device-bound P path (gpu_aware_mpi_runtime + dev_halo_active)
  !> the P requests are DEVICE-bound to the SAME P_send_buf/P_recv_buf. This is SAFE in
  !> every current device-runnable case because:
  !>   (1) IRS (the only caller of exchange_ghost_R, via Residual_Smoothing) is gated by
  !>       obj_irs%IRS, which DEFAULTS .false. and is OFF for all device cases (SpicchiOThor,
  !>       dettube, RDE) -> exchange_ghost_R never executes there, AND
  !>   (2) even with IRS on, P (RK-substage Fill_Ghost_Cell) and R (post-Newstate residual
  !>       smoothing) exchange in temporally disjoint phases -> no overlap on the shared buffer.
  !> If IRS is ever enabled together with the device-bound P path, the host-bound R exchange
  !> would read/write the HOST buffer while the DEVICE copy is authoritative for P, with no
  !> H2D/D2H sync bridging the two uses -> add an explicit stage (or device-bind R) then.
  subroutine init_persistent_R_requests()
    use MOSE_Global_m, only: gc, nprim
    use mpi

    implicit none
    integer :: r, ierr, tag, buf_pos, entry_size
    integer, parameter :: TAG_OFFSET = 30000

    entry_size = nprim * gc

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return

    allocate(ghost_sched%R_recv_req_pers(max(ghost_sched%n_recv_ranks, 1)))
    allocate(ghost_sched%R_send_req_pers(max(ghost_sched%n_send_ranks, 1)))

    ! Initialize persistent receives for R
    do r = 1, ghost_sched%n_recv_ranks
      buf_pos = (ghost_sched%recv_groups(r)%offset - 1) * entry_size + 1
      tag = ghost_sched%recv_groups(r)%rank + TAG_OFFSET
      call MPI_RECV_INIT(ghost_sched%P_recv_buf(buf_pos), &
                         ghost_sched%recv_groups(r)%count * entry_size, &
                         MPI_DOUBLE_PRECISION, &
                         ghost_sched%recv_groups(r)%rank, tag, &
                         MPI_COMM_WORLD, ghost_sched%R_recv_req_pers(r), ierr)
      call check_mpi_error(ierr)
    end do

    ! Initialize persistent sends for R
    do r = 1, ghost_sched%n_send_ranks
      buf_pos = (ghost_sched%send_groups(r)%offset - 1) * entry_size + 1
      tag = mpi_rank_ + TAG_OFFSET
      call MPI_SEND_INIT(ghost_sched%P_send_buf(buf_pos), &
                         ghost_sched%send_groups(r)%count * entry_size, &
                         MPI_DOUBLE_PRECISION, &
                         ghost_sched%send_groups(r)%rank, tag, &
                         MPI_COMM_WORLD, ghost_sched%R_send_req_pers(r), ierr)
      call check_mpi_error(ierr)
    end do

    ghost_sched%persistent_R_init = .true.
  end subroutine init_persistent_R_requests


  !> Free all persistent MPI requests.
  subroutine cleanup_persistent_requests()
    use mpi
    implicit none
    integer :: ierr, r

    if (ghost_sched%persistent_P_init) then
      do r = 1, ghost_sched%n_send_ranks
        call MPI_REQUEST_FREE(ghost_sched%P_send_req_pers(r), ierr)
      end do
      do r = 1, ghost_sched%n_recv_ranks
        call MPI_REQUEST_FREE(ghost_sched%P_recv_req_pers(r), ierr)
      end do
      ghost_sched%persistent_P_init = .false.
    end if

    if (ghost_sched%persistent_R_init) then
      do r = 1, ghost_sched%n_send_ranks
        call MPI_REQUEST_FREE(ghost_sched%R_send_req_pers(r), ierr)
      end do
      do r = 1, ghost_sched%n_recv_ranks
        call MPI_REQUEST_FREE(ghost_sched%R_recv_req_pers(r), ierr)
      end do
      ghost_sched%persistent_R_init = .false.
    end if

    !  device transverse-Pg persistent requests.
    if (ghost_sched%persistent_Pg_init) then
      do r = 1, ghost_sched%n_send_ranks
        call MPI_REQUEST_FREE(ghost_sched%Pg_send_req_pers(r), ierr)
      end do
      do r = 1, ghost_sched%n_recv_ranks
        call MPI_REQUEST_FREE(ghost_sched%Pg_recv_req_pers(r), ierr)
      end do
      ghost_sched%persistent_Pg_init = .false.
    end if
  end subroutine cleanup_persistent_requests


  !> Sort ghost entries by remote_rank using insertion sort (lists are small).
  subroutine sort_entries_by_rank(list, n)
    implicit none
    type(ghost_entry_type), intent(inout) :: list(:)
    integer, intent(in) :: n
    ! Local
    integer :: i, j
    type(ghost_entry_type) :: tmp

    do i = 2, n
      tmp = list(i)
      j = i - 1
      do while (j >= 1)
        if (list(j)%remote_rank <= tmp%remote_rank) exit
        list(j+1) = list(j)
        j = j - 1
      end do
      list(j+1) = tmp
    end do
  end subroutine sort_entries_by_rank


  !> Build rank groups from a sorted entry list.
  subroutine build_rank_groups(list, n, groups, n_groups)
    implicit none
    type(ghost_entry_type), intent(in) :: list(:)
    integer, intent(in) :: n
    type(rank_group_type), allocatable, intent(out) :: groups(:)
    integer, intent(out) :: n_groups
    ! Local
    integer :: i, ng, cur_rank

    if (n == 0) then
      n_groups = 0
      allocate(groups(0))
      return
    end if

    ! Count unique ranks
    ng = 1
    do i = 2, n
      if (list(i)%remote_rank /= list(i-1)%remote_rank) ng = ng + 1
    end do
    n_groups = ng

    allocate(groups(ng))
    ng = 1
    groups(1)%rank   = list(1)%remote_rank
    groups(1)%offset = 1
    groups(1)%count  = 1
    do i = 2, n
      if (list(i)%remote_rank /= list(i-1)%remote_rank) then
        ng = ng + 1
        groups(ng)%rank   = list(i)%remote_rank
        groups(ng)%offset = i
        groups(ng)%count  = 1
      else
        groups(ng)%count = groups(ng)%count + 1
      end if
    end do
  end subroutine build_rank_groups


  !> Build face groups from send entry list.
  !> Groups consecutive entries that share (local_block, src_f).
  subroutine build_face_groups_send(list, n, groups, n_groups)
    use MOSE_Global_m, only: guide
    implicit none
    type(ghost_entry_type), intent(in) :: list(:)
    integer, intent(in) :: n
    type(face_group_type), allocatable, intent(out) :: groups(:)
    integer, intent(out) :: n_groups
    integer :: i, ng

    if (n == 0) then
      n_groups = 0
      allocate(groups(0))
      return
    end if

    ! Count face groups
    ng = 1
    do i = 2, n
      if (list(i)%local_block /= list(i-1)%local_block .or. &
          list(i)%src_f /= list(i-1)%src_f) ng = ng + 1
    end do
    n_groups = ng

    allocate(groups(ng))
    ng = 1
    groups(1)%block = list(1)%local_block
    groups(1)%face  = list(1)%src_f
    groups(1)%entry_start = 1
    groups(1)%entry_count = 1
    groups(1)%di = guide(list(1)%src_f, 1)
    groups(1)%dj = guide(list(1)%src_f, 2)
    groups(1)%dk = guide(list(1)%src_f, 3)
    do i = 2, n
      if (list(i)%local_block /= list(i-1)%local_block .or. &
          list(i)%src_f /= list(i-1)%src_f) then
        ng = ng + 1
        groups(ng)%block = list(i)%local_block
        groups(ng)%face  = list(i)%src_f
        groups(ng)%entry_start = i
        groups(ng)%entry_count = 1
        groups(ng)%di = guide(list(i)%src_f, 1)
        groups(ng)%dj = guide(list(i)%src_f, 2)
        groups(ng)%dk = guide(list(i)%src_f, 3)
      else
        groups(ng)%entry_count = groups(ng)%entry_count + 1
      end if
    end do
  end subroutine build_face_groups_send


  !> Build face groups from recv entry list.
  !> Groups consecutive entries that share (local_block, dst_f).
  subroutine build_face_groups_recv(list, n, groups, n_groups)
    use MOSE_Global_m, only: guide
    implicit none
    type(ghost_entry_type), intent(in) :: list(:)
    integer, intent(in) :: n
    type(face_group_type), allocatable, intent(out) :: groups(:)
    integer, intent(out) :: n_groups
    integer :: i, ng

    if (n == 0) then
      n_groups = 0
      allocate(groups(0))
      return
    end if

    ! Count face groups
    ng = 1
    do i = 2, n
      if (list(i)%local_block /= list(i-1)%local_block .or. &
          list(i)%dst_f /= list(i-1)%dst_f) ng = ng + 1
    end do
    n_groups = ng

    allocate(groups(ng))
    ng = 1
    groups(1)%block = list(1)%local_block
    groups(1)%face  = list(1)%dst_f
    groups(1)%entry_start = 1
    groups(1)%entry_count = 1
    groups(1)%di = guide(list(1)%dst_f, 1)
    groups(1)%dj = guide(list(1)%dst_f, 2)
    groups(1)%dk = guide(list(1)%dst_f, 3)
    do i = 2, n
      if (list(i)%local_block /= list(i-1)%local_block .or. &
          list(i)%dst_f /= list(i-1)%dst_f) then
        ng = ng + 1
        groups(ng)%block = list(i)%local_block
        groups(ng)%face  = list(i)%dst_f
        groups(ng)%entry_start = i
        groups(ng)%entry_count = 1
        groups(ng)%di = guide(list(i)%dst_f, 1)
        groups(ng)%dj = guide(list(i)%dst_f, 2)
        groups(ng)%dk = guide(list(i)%dst_f, 3)
      else
        groups(ng)%entry_count = groups(ng)%entry_count + 1
      end if
    end do
  end subroutine build_face_groups_recv
#endif


  !> Post persistent P receives. Call from !$omp single.
  subroutine exchange_ghost_P_post_recv(domain)
    use MOSE_Types_m
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call exchange_P_post_recv()
#endif
  end subroutine exchange_ghost_P_post_recv



  !> Post persistent P sends. Call from !$omp single after pack is complete.
  !> dev_packed: .true. when the caller is the DEVICE halo driver
  !> (Fill_Ghost_Cell_mr), which has already packed the DEVICE P_send_buf via
  !> exchange_ghost_P_pack. In that case the device buffer is authoritative and the
  !> H2D send-buffer bridge MUST NOT run (it would clobber the freshly device-packed data
  !> with stale HOST bytes the device driver never wrote -> garbage ghosts (a halo-correctness bug).
  !> The bridge is for the OTHER caller only: the HOST
  !> Fill_Ghost_Cell path (setup/Compute_Residual), which packs the HOST buffer and, when
  !> the requests are device-bound, needs it staged down. Optional => host path unchanged.
  subroutine exchange_ghost_P_post_send(domain, dev_packed)
    use MOSE_Types_m
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    logical, intent(in), optional :: dev_packed
    logical :: host_packed
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    ! Default (absent) => HOST caller packed the host buffer => bridge allowed.
    host_packed = .true.
    if (present(dev_packed)) host_packed = .not. dev_packed
    !  crash #2: the HOST Fill_Ghost_Cell path packs into the HOST P_send_buf. When the
    ! persistent requests are DEVICE-bound (gpu_aware_mpi_runtime + dev_halo_active), MPI DMAs
    ! the DEVICE buffer, so push the host-packed data down first. No-op when host-bound
    ! (bridge guards on bufs_resident; flag false unless device-bound). Skipped entirely when
    ! dev_packed (device driver) -- the device buffer is already authoritative.
    if (host_packed .and. ghost_sched%dev_bound) &
        call exchange_P_send_buf_to_device()
    call exchange_P_post_send()
#endif
  end subroutine exchange_ghost_P_post_send



  !> Wait for P receives ONLY (no unpack). Lets the device-halo driver
  !> sequence wait -> (stage) -> device-unpack / host-unpack explicitly. !$omp single.
  subroutine exchange_ghost_P_wait_recv_only(domain)
    use MOSE_Types_m
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call exchange_P_wait_recv()
#endif
  end subroutine exchange_ghost_P_wait_recv_only



  !> Wait for P sends to complete. Call from !$omp single.
  subroutine exchange_ghost_P_wait_send(domain)
    use MOSE_Types_m
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call exchange_P_wait_send()
#endif
  end subroutine exchange_ghost_P_wait_send


  !> Pack the P send buffer ON DEVICE from the resident P, then (when CUDA-
  !> aware MPI is OFF) stage it to host so the host-bound MPI send can read it. Call from
  !> !$omp single. CPU build: not referenced. With CUDA-aware MPI no staging happens --
  !> the device-bound persistent send reads device memory directly.
  subroutine exchange_ghost_P_pack(domain, async_q)
    use MOSE_Types_m
    use openacc, only: acc_async_sync
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer, intent(in), optional :: async_q   !< mr-fuse-exchange: one launch on this queue, then ONE wait on it
    if (mpi_size_ <= 1) return
#if defined(USE_MPI)
    call exchange_P_pack(domain, async_q)
    ! T1c: a host-bound schedule copies the packed buffer down ON the pack's queue, before the one wait below: the wait then
    ! covers the pack and the copy, and the host-bound MPI send reads the host buffer after it -- no host synchronisation added
    if (.not. ghost_sched%dev_bound) call exchange_P_send_buf_to_host(async_q)
    ! the wait is unconditional (also when this rank sends nothing at the level): it is the join of
    ! the stage queue that the local fill on Q_OVERLAP relies on
    if (present(async_q)) then
      if (async_q /= acc_async_sync) then
        !$acc wait(async_q)
      end if
    end if
#endif
  end subroutine exchange_ghost_P_pack


  !> (when CUDA-aware MPI is OFF) stage the received host buffer to device,
  !> then unpack it ON DEVICE into the resident P ghost cells. Call from !$omp single
  !> AFTER the P receives have completed. CPU build: not referenced.
  subroutine exchange_ghost_P_unpack(domain, async_q)
    use MOSE_Types_m
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in), optional :: async_q   !< mr-fuse-exchange: one launch on this queue, no wait
    if (mpi_size_ <= 1) return
#if defined(USE_MPI)
    ! T1c: the received host buffer goes up on the unpack's own queue, right before the unpack kernel (no host wait); the host
    ! buffer is not touched again before the stage's final wait (the next MPI_STARTALL of this schedule comes after it)
    if (.not. ghost_sched%dev_bound) call exchange_P_recv_buf_to_device(async_q)
    call exchange_P_unpack(domain, async_q)
#endif
  end subroutine exchange_ghost_P_unpack




  !> Gather all blocks' P interior data to rank 0 for I/O.
  !> Uses non-blocking MPI_ISEND/IRECV + WAITALL for overlap.
  !> Temporarily allocates P on root for remote blocks if not already allocated.
  subroutine gather_P_to_root(domain)
    use MOSE_Types_m
    use MOSE_Global_m, only: nprim, gc

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, n1, n2, n3, ncells

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr, nreq, r, total_buf_size, offset
      integer, allocatable :: reqs(:), blk_offset(:), blk_ncells(:)
      integer, allocatable :: stats(:,:)
      real(R8), allocatable :: buf(:)

      ! Compute per-block sizes and total buffer size
      allocate(blk_offset(domain%nb), blk_ncells(domain%nb))
      total_buf_size = 0
      do b = 1, domain%nb
        n1 = domain%blk(b)%dim(1)
        n2 = domain%blk(b)%dim(2)
        n3 = domain%blk(b)%dim(3)
        blk_ncells(b) = nprim * n1 * n2 * n3
        blk_offset(b) = total_buf_size
        if (mpi_is_root) then
          if (.not. is_local_block(b)) total_buf_size = total_buf_size + blk_ncells(b)
        else
          if (is_local_block(b)) total_buf_size = total_buf_size + blk_ncells(b)
        end if
      end do
      allocate(buf(max(total_buf_size, 1)))
      allocate(reqs(domain%nb))
      allocate(stats(MPI_STATUS_SIZE, domain%nb))

      ! Recompute offsets contiguously
      offset = 0
      do b = 1, domain%nb
        blk_offset(b) = offset
        if (mpi_is_root) then
          if (.not. is_local_block(b)) offset = offset + blk_ncells(b)
        else
          if (is_local_block(b)) offset = offset + blk_ncells(b)
        end if
      end do

      nreq = 0
      if (mpi_is_root) then
        ! Root: post all receives for remote blocks
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            nreq = nreq + 1
            call MPI_IRECV(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           block_owner(b), b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      else
        ! Non-root: pack local blocks and post all sends
        do b = 1, domain%nb
          if (is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            buf(blk_offset(b)+1 : blk_offset(b)+blk_ncells(b)) = &
              reshape(domain%blk(b)%P(:, 1:n1, 1:n2, 1:n3), [blk_ncells(b)])
            nreq = nreq + 1
            call MPI_ISEND(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           0, b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      end if

      ! Wait for all
      if (nreq > 0) then
        call MPI_WAITALL(nreq, reqs(1:nreq), stats(:,1:nreq), ierr)
        call check_mpi_error(ierr)
      end if

      ! Root: unpack received data (allocate P if it was previously freed)
      if (mpi_is_root) then
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            if (.not. allocated(domain%blk(b)%P)) &
              allocate(domain%blk(b)%P(nprim, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc))
            domain%blk(b)%P(:, 1:n1, 1:n2, 1:n3) = &
              reshape(buf(blk_offset(b)+1 : blk_offset(b)+blk_ncells(b)), [nprim, n1, n2, n3])
          end if
        end do
      end if

      deallocate(buf, reqs, stats, blk_offset, blk_ncells)
    end block
#endif
  end subroutine gather_P_to_root


  !> Gather diagnostic arrays (R, dtlocal, beta) from all ranks to root.
  !> On root, remote blocks may have had these arrays deallocated by
  !> deallocate_remote_computation_data, so they are re-allocated here
  !> to receive the data.  All ranks must call this (collective).
  subroutine gather_diagnostic_to_root(domain)
    use MOSE_Types_m
    use MOSE_Global_m, only: nprim

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: b, n1, n2, n3

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr, nreq, offset, total_buf_size
      integer :: ncells_R, ncells_dt
      integer, allocatable :: reqs(:)
      integer, allocatable :: stats(:,:)
      real(R8), allocatable :: buf(:)
      ! Per-block: nprim*n1*n2*n3 for R, n1*n2*n3 for dtlocal, n1*n2*n3 for beta
      integer, allocatable :: blk_offset(:), blk_size(:)

      allocate(blk_offset(domain%nb), blk_size(domain%nb))

      ! Compute buffer sizes
      total_buf_size = 0
      do b = 1, domain%nb
        n1 = domain%blk(b)%dim(1)
        n2 = domain%blk(b)%dim(2)
        n3 = domain%blk(b)%dim(3)
        blk_size(b) = nprim * n1 * n2 * n3 + 2 * n1 * n2 * n3  ! R + dtlocal + beta
      end do

      offset = 0
      do b = 1, domain%nb
        blk_offset(b) = offset
        if (mpi_is_root) then
          if (.not. is_local_block(b)) offset = offset + blk_size(b)
        else
          if (is_local_block(b)) offset = offset + blk_size(b)
        end if
      end do
      total_buf_size = offset
      allocate(buf(max(total_buf_size, 1)))
      allocate(reqs(domain%nb))
      allocate(stats(MPI_STATUS_SIZE, domain%nb))

      ! On root: ensure R, dtlocal, beta are allocated for remote blocks
      if (mpi_is_root) then
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            if (.not. allocated(domain%blk(b)%R)) &
              allocate(domain%blk(b)%R(nprim, n1, n2, n3))
            if (.not. allocated(domain%blk(b)%dtlocal)) &
              allocate(domain%blk(b)%dtlocal(n1, n2, n3))
            if (.not. allocated(domain%blk(b)%beta)) &
              allocate(domain%blk(b)%beta(n1, n2, n3))
          end if
        end do
      end if

      nreq = 0
      if (mpi_is_root) then
        ! Root: post receives for remote blocks
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            nreq = nreq + 1
            call MPI_IRECV(buf(blk_offset(b)+1), blk_size(b), MPI_DOUBLE_PRECISION, &
                           block_owner(b), 1000+b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      else
        ! Non-root: pack local blocks and send
        do b = 1, domain%nb
          if (is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            ncells_R  = nprim * n1 * n2 * n3
            ncells_dt = n1 * n2 * n3
            offset = blk_offset(b)
            buf(offset+1 : offset+ncells_R) = &
              reshape(domain%blk(b)%R(:, 1:n1, 1:n2, 1:n3), [ncells_R])
            offset = offset + ncells_R
            buf(offset+1 : offset+ncells_dt) = &
              reshape(domain%blk(b)%dtlocal(1:n1, 1:n2, 1:n3), [ncells_dt])
            offset = offset + ncells_dt
            buf(offset+1 : offset+ncells_dt) = &
              reshape(domain%blk(b)%beta(1:n1, 1:n2, 1:n3), [ncells_dt])
            nreq = nreq + 1
            call MPI_ISEND(buf(blk_offset(b)+1), blk_size(b), MPI_DOUBLE_PRECISION, &
                           0, 1000+b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      end if

      ! Wait for all
      if (nreq > 0) then
        call MPI_WAITALL(nreq, reqs(1:nreq), stats(:,1:nreq), ierr)
        call check_mpi_error(ierr)
      end if

      ! Root: unpack received data
      if (mpi_is_root) then
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            ncells_R  = nprim * n1 * n2 * n3
            ncells_dt = n1 * n2 * n3
            offset = blk_offset(b)
            domain%blk(b)%R(:, 1:n1, 1:n2, 1:n3) = &
              reshape(buf(offset+1 : offset+ncells_R), [nprim, n1, n2, n3])
            offset = offset + ncells_R
            domain%blk(b)%dtlocal(1:n1, 1:n2, 1:n3) = &
              reshape(buf(offset+1 : offset+ncells_dt), [n1, n2, n3])
            offset = offset + ncells_dt
            domain%blk(b)%beta(1:n1, 1:n2, 1:n3) = &
              reshape(buf(offset+1 : offset+ncells_dt), [n1, n2, n3])
          end if
        end do
      end if

      deallocate(buf, reqs, stats, blk_offset, blk_size)
    end block
#endif
  end subroutine gather_diagnostic_to_root




#ifdef USE_MPI

  !> Post persistent receives for P field. Call from !$omp single.
  subroutine exchange_P_post_recv()
    use mpi
    implicit none
    integer :: ierr

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return
    if (ghost_sched%n_recv_ranks > 0) then
      call MPI_STARTALL(ghost_sched%n_recv_ranks, ghost_sched%P_recv_req_pers, ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_P_post_recv



  !> Post persistent sends for P field. Call from !$omp single after pack.
  subroutine exchange_P_post_send()
    use mpi
    implicit none
    integer :: ierr

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return
    if (ghost_sched%n_send_ranks > 0) then
      call MPI_STARTALL(ghost_sched%n_send_ranks, ghost_sched%P_send_req_pers, ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_P_post_send


  !> Wait for P receives and unpack into ghost cells. Safe to call unpack from !$omp do.
  subroutine exchange_P_wait_recv()
    use mpi
    implicit none
    integer :: ierr

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return
    if (ghost_sched%n_recv_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_recv_ranks, ghost_sched%P_recv_req_pers, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_recv_ranks), ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_P_wait_recv



  !> Wait for P sends to complete. Call from !$omp single.
  subroutine exchange_P_wait_send()
    use mpi
    implicit none
    integer :: ierr

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return
    if (ghost_sched%n_send_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_send_ranks, ghost_sched%P_send_req_pers, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_send_ranks), ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_P_wait_send


  !----------------------------------------------------------------------------
  ! GPU Phase 1 device halo: pack/unpack the P send/recv buffers directly out of /
  ! into the device-RESIDENT blk(b)%P, plus an optional host-staging bridge for the
  ! CUDA-aware-MPI fallback. Mirrors STREAmS bcswap_step_1_cuf (device pack) but in the
  ! current OpenACC style ( Phase 1). One kernel per LOCAL block so a single
  ! present-matched P maps in. The flat send_*/recv_* index tables and
  ! P_send_buf/P_recv_buf are resident; buffer slot (i-1)*nprim*gc matches the host path.
  !----------------------------------------------------------------------------

  !> Pack the resident P halo into the resident P_send_buf, one kernel per local block.
  !> Call from !$omp single (mirrors Fill_Ghost_Cell_acc). The flat index tables and the
  !> send buffer are passed as explicit ARRAY DUMMIES (not through ghost_sched on the
  !> device) -- the proven present-matching pattern that avoids nvfortran's
  !> derived-type-with-allocatable-components device-mapping problem.
  subroutine exchange_P_pack(domain, async_q)
    use MOSE_Types_m
    use MOSE_Global_m, only: gc, nprim
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer, intent(in), optional :: async_q
    integer :: b

    if (ghost_sched%n_send == 0) return
    if (present(async_q)) then
      block
        use MOSE_Field_m, only: P_flat, ncell_total, dboff, dpi, dpj
        call exchange_P_pack_all(P_flat, ncell_total, dboff, dpi, dpj, size(dboff), ghost_sched%n_send, &
                                   size(ghost_sched%P_send_buf), ghost_sched%P_send_buf, ghost_sched%send_b, &
                                   ghost_sched%send_is, ghost_sched%send_js, ghost_sched%send_ks, &
                                   ghost_sched%send_di, ghost_sched%send_dj, ghost_sched%send_dk, async_q)
      end block
      return
    end if
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      block
        use MOSE_Field_m, only: P_flat, ncell_total, field_block_index_params
        integer :: boff, pi, pj
        call field_block_index_params(b, boff, pi, pj)
        call exchange_P_pack_blk(P_flat, ncell_total, boff, pi, pj, b, ghost_sched%n_send, &
                                       size(ghost_sched%P_send_buf), ghost_sched%P_send_buf, ghost_sched%send_b, &
                                       ghost_sched%send_is, ghost_sched%send_js, ghost_sched%send_ks, &
                                       ghost_sched%send_di, ghost_sched%send_dj, ghost_sched%send_dk)
      end block
    end do
  end subroutine exchange_P_pack



  !> Packs the donor interior rows of one block from w into the send buffer.
  subroutine exchange_P_pack_blk(w, ncw, boff, pi, pj, b, ns, nbuf, sbuf, sb, sis, sjs, sks, sdi, sdj, sdk)
    use MOSE_Global_m, only: gc, nprim
    implicit none
    integer,  intent(in)    :: ncw, boff, pi, pj, b, ns, nbuf
    real(R8), intent(in)    :: w(ncw, nprim)
    real(R8), intent(inout) :: sbuf(nbuf)
    integer,  intent(in)    :: sb(ns), sis(ns), sjs(ns), sks(ns), sdi(ns), sdj(ns), sdk(ns)
    integer :: i, g, s, Is, Js, Ks, di, dj, dk, base, ii, jj, kk

    !$acc parallel loop gang vector present(w, sbuf, sb, sis, sjs, sks, sdi, sdj, sdk) &
    !$acc&   private(g, s, Is, Js, Ks, di, dj, dk, base, ii, jj, kk)
    do i = 1, ns
      if (sb(i) /= b) cycle
      Is = sis(i); Js = sjs(i); Ks = sks(i)
      di = sdi(i); dj = sdj(i); dk = sdk(i)
      base = (i - 1) * nprim * gc
      do g = 1, gc
        ii = Is + di*(g-1); jj = Js + dj*(g-1); kk = Ks + dk*(g-1)
        do s = 1, nprim
          sbuf(base + (g-1)*nprim + s) = WGCELL(s, ii, jj, kk)
        end do
      end do
    end do
  end subroutine exchange_P_pack_blk


  !> Unpack the resident P_recv_buf into the resident P ghost cells, one kernel per local
  !> block. Call from !$omp single. Mirrors the host exchange_P_unpack_faces gather.
  subroutine exchange_P_unpack(domain, async_q)
    use MOSE_Types_m
    use MOSE_Global_m, only: gc, nprim
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in), optional :: async_q
    integer :: b

    if (ghost_sched%n_recv == 0) return
    if (present(async_q)) then
      block
        use MOSE_Field_m, only: P_flat, ncell_total, dboff, dpi, dpj
        call exchange_P_unpack_all(P_flat, ncell_total, dboff, dpi, dpj, size(dboff), ghost_sched%n_recv, &
                                     size(ghost_sched%P_recv_buf), ghost_sched%P_recv_buf, ghost_sched%recv_b, &
                                     ghost_sched%recv_id, ghost_sched%recv_jd, ghost_sched%recv_kd, &
                                     ghost_sched%recv_di, ghost_sched%recv_dj, ghost_sched%recv_dk, async_q)
      end block
      return
    end if
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      block
        use MOSE_Field_m, only: P_flat, ncell_total, field_block_index_params
        integer :: boff, pi, pj
        call field_block_index_params(b, boff, pi, pj)
        call exchange_P_unpack_blk(P_flat, ncell_total, boff, pi, pj, b, ghost_sched%n_recv, &
                                         size(ghost_sched%P_recv_buf), ghost_sched%P_recv_buf, ghost_sched%recv_b, &
                                         ghost_sched%recv_id, ghost_sched%recv_jd, ghost_sched%recv_kd, &
                                         ghost_sched%recv_di, ghost_sched%recv_dj, ghost_sched%recv_dk)
      end block
    end do
  end subroutine exchange_P_unpack



  !> Writes the ghost rows of one block from the receive buffer into w.
  subroutine exchange_P_unpack_blk(w, ncw, boff, pi, pj, b, nr, nbuf, rbuf, rb, rid, rjd, rkd, rdi, rdj, rdk)
    use MOSE_Global_m, only: gc, nprim
    implicit none
    integer,  intent(in)    :: ncw, boff, pi, pj, b, nr, nbuf
    real(R8), intent(inout) :: w(ncw, nprim)
    real(R8), intent(in)    :: rbuf(nbuf)
    integer,  intent(in)    :: rb(nr), rid(nr), rjd(nr), rkd(nr), rdi(nr), rdj(nr), rdk(nr)
    integer :: i, g, s, Id, Jd, Kd, di, dj, dk, base, ii, jj, kk

    !$acc parallel loop gang vector present(w, rbuf, rb, rid, rjd, rkd, rdi, rdj, rdk) &
    !$acc&   private(g, s, Id, Jd, Kd, di, dj, dk, base, ii, jj, kk)
    do i = 1, nr
      if (rb(i) /= b) cycle
      Id = rid(i); Jd = rjd(i); Kd = rkd(i)
      di = rdi(i); dj = rdj(i); dk = rdk(i)
      base = (i - 1) * nprim * gc
      do g = 1, gc
        ii = Id - di*g; jj = Jd - dj*g; kk = Kd - dk*g
        do s = 1, nprim
          WGCELL(s, ii, jj, kk) = rbuf(base + (g-1)*nprim + s)
        end do
      end do
    end do
  end subroutine exchange_P_unpack_blk


  !> Host-staging bridge for the CUDA-aware-MPI FALLBACK: copy the packed P_send_buf from
  !> device to host before the (host-bound) MPI send. Used only when gpu_aware_mpi_runtime
  !> is .false. (machine without CUDA-aware MPI). With CUDA-aware MPI this is never called
  !> -- the persistent requests are device-bound and MPI reads device memory directly.
  subroutine exchange_P_send_buf_to_host(async_q)
    use openacc, only: acc_async_sync
    implicit none
    integer, intent(in), optional :: async_q   !< T1c: the stage queue the copy rides (absent = synchronous, as before)
    integer :: q
    if (ghost_sched%n_send == 0) return
    if (.not. ghost_sched%bufs_resident) return
    q = acc_async_sync
    if (present(async_q)) q = async_q
    !$acc update host(ghost_sched%P_send_buf) async(q)
  end subroutine exchange_P_send_buf_to_host


  !> Host-staging bridge (fallback): copy the received P_recv_buf from host to device
  !> after the (host-bound) MPI recv, so the device unpack kernel reads current data.
  subroutine exchange_P_recv_buf_to_device(async_q)
    use openacc, only: acc_async_sync
    implicit none
    integer, intent(in), optional :: async_q   !< T1c: the stage queue the copy rides (absent = synchronous, as before)
    integer :: q
    if (ghost_sched%n_recv == 0) return
    if (.not. ghost_sched%bufs_resident) return
    q = acc_async_sync
    if (present(async_q)) q = async_q
    !$acc update device(ghost_sched%P_recv_buf) async(q)
  end subroutine exchange_P_recv_buf_to_device


  !>  crash #2 bridge: when the persistent P requests are DEVICE-bound
  !> (gpu_aware_mpi_runtime + dev_halo_active) yet the HOST Fill_Ghost_Cell path runs (the
  !> setup-time / Compute_Residual calls, where P is not resident so the device halo driver
  !> is skipped), the host packs into the HOST P_send_buf but MPI DMAs the DEVICE buffer.
  !> Push the host-packed send buffer to the device copy BEFORE the device-bound MPI send.
  !> Mirror image of exchange_P_send_buf_to_host (which serves the opposite, device-driver
  !> + host-bound-MPI fallback). No-op unless the buffers are resident (device-bound case).
  subroutine exchange_P_send_buf_to_device()
    implicit none
    if (ghost_sched%n_send == 0) return
    if (.not. ghost_sched%bufs_resident) return
    !$acc update device(ghost_sched%P_send_buf)
  end subroutine exchange_P_send_buf_to_device



  !----------------------------------------------------------------------------
  ! Device transverse-Pg pack / exchange / unpack.
  ! Removes the host transverse gather (Fill_BC_Ghost_Connection) and therefore the
  ! step-0 whole-P D2H that fed it. The donor-side transverse stencil is gathered on
  ! DEVICE out of the resident donor blk(Bs)%P (mirror of host Fill_BC_Ghost_Connection),
  ! shipped by GPU-aware MPI, and unpacked on DEVICE straight into the resident bc_pg(:,3:6).
  ! The receiver-normal bc_pg(:,1:gc) is filled on DEVICE from the (already device-unpacked)
  ! receiver ghost P, replacing the host Ghost_Interrank + per-slot update device(bc_pg).
  ! One kernel per LOCAL donor block so a single present-matched donor P maps in.
  !----------------------------------------------------------------------------

  !> Pack the resident Pg send buffer from the resident donor P, one kernel per local block.
  subroutine exchange_Pg_pack(domain, async_q)
    use MOSE_Types_m
    use MOSE_Global_m, only: nprim
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer, intent(in), optional :: async_q
    integer :: b

    if (ghost_sched%n_send == 0) return
    if (present(async_q)) then
      block
        use MOSE_Field_m, only: P_flat, ncell_total, dboff, dpi, dpj
        call exchange_Pg_pack_all(P_flat, ncell_total, dboff, dpi, dpj, size(dboff), ghost_sched%n_send, &
                                    size(ghost_sched%Pg_send_buf), ghost_sched%Pg_send_buf, ghost_sched%send_b, &
                                    ghost_sched%send_is, ghost_sched%send_js, ghost_sched%send_ks, &
                                    ghost_sched%pg_send_fs, ghost_sched%pg_send_d11, ghost_sched%pg_send_d12, &
                                    ghost_sched%pg_send_d21, ghost_sched%pg_send_d22, async_q)
      end block
      return
    end if
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      block
        use MOSE_Field_m, only: P_flat, ncell_total, field_block_index_params
        integer :: boff, pi, pj
        call field_block_index_params(b, boff, pi, pj)
        call exchange_Pg_pack_blk(P_flat, ncell_total, boff, pi, pj, b, ghost_sched%n_send, &
                                        size(ghost_sched%Pg_send_buf), ghost_sched%Pg_send_buf, ghost_sched%send_b, &
                                        ghost_sched%send_is, ghost_sched%send_js, ghost_sched%send_ks, &
                                        ghost_sched%pg_send_fs, ghost_sched%pg_send_d11, ghost_sched%pg_send_d12, &
                                        ghost_sched%pg_send_d21, ghost_sched%pg_send_d22)
      end block
    end do
  end subroutine exchange_Pg_pack



  !> Packs the transverse stencil of one block from w. The transverse stencil can reach
  !> donor EDGE GHOSTS (filled in w by Fill_Ghost_Cell_local, fenced before this pack), hence
  !> it MUST read w (the n>1 w-halo no longer fills donor P ghosts).
  subroutine exchange_Pg_pack_blk(w, ncw, boff, pi, pj, b, ns, nbuf, sbuf, sb, sis, sjs, sks, sfs, sd11, sd12, sd21, sd22)
    use MOSE_Global_m, only: gc, nprim
    implicit none
    integer,  intent(in)    :: ncw, boff, pi, pj, b, ns, nbuf
    real(R8), intent(in)    :: w(ncw, nprim)
    real(R8), intent(inout) :: sbuf(nbuf)
    integer,  intent(in)    :: sb(ns), sis(ns), sjs(ns), sks(ns), sfs(ns)
    integer,  intent(in)    :: sd11(ns), sd12(ns), sd21(ns), sd22(ns)
    integer :: i, s, Is, Js, Ks, Fs, d11, d12, d21, d22, base
    integer :: i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4

    !$acc parallel loop gang vector present(w, sbuf, sb, sis, sjs, sks, sfs, sd11, sd12, sd21, sd22) &
    !$acc&   private(s, Is, Js, Ks, Fs, d11, d12, d21, d22, base, &
    !$acc&           i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4)
    do i = 1, ns
      if (sb(i) /= b) cycle
      Is = sis(i); Js = sjs(i); Ks = sks(i); Fs = sfs(i)
      d11 = sd11(i); d12 = sd12(i); d21 = sd21(i); d22 = sd22(i)
      select case (Fs)
        case (1:2)
          i1=Is; j1=Js-d11; k1=Ks-d12;  i2=Is; j2=Js+d11; k2=Ks+d12
          i3=Is; j3=Js-d21; k3=Ks-d22;  i4=Is; j4=Js+d21; k4=Ks+d22
        case (3:4)
          i1=Is-d11; j1=Js; k1=Ks-d12;  i2=Is+d11; j2=Js; k2=Ks+d12
          i3=Is-d21; j3=Js; k3=Ks-d22;  i4=Is+d21; j4=Js; k4=Ks+d22
        case default  ! (5:6)
          i1=Is-d11; j1=Js-d12; k1=Ks;  i2=Is+d11; j2=Js+d12; k2=Ks
          i3=Is-d21; j3=Js-d22; k3=Ks;  i4=Is+d21; j4=Js+d22; k4=Ks
      end select
      base = (i - 1) * nprim * 4
      do s = 1, nprim
        sbuf(base          + s) = WGCELL(s, i1, j1, k1)   ! -> Pg(:,3)
        sbuf(base +   nprim + s) = WGCELL(s, i2, j2, k2)  ! -> Pg(:,4)
        sbuf(base + 2*nprim + s) = WGCELL(s, i3, j3, k3)  ! -> Pg(:,5)
        sbuf(base + 3*nprim + s) = WGCELL(s, i4, j4, k4)  ! -> Pg(:,6)
      end do
    end do
  end subroutine exchange_Pg_pack_blk


  !> Unpack the resident Pg_recv_buf into the resident bc_pg(:,3:6,slot), one kernel call.
  !> The destination slot is the flat bc_pg slot (bc_perm inverse) of each recv entry.
  subroutine exchange_Pg_unpack(domain, async_q)
    use MOSE_Types_m
    use MOSE_Global_m, only: nprim
    use openacc, only: acc_async_sync
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in), optional :: async_q
    integer :: q
    q = acc_async_sync
    if (present(async_q)) q = async_q
    if (ghost_sched%n_recv == 0) return
    call exchange_Pg_unpack_buf(domain % bc_pg, domain % nbound, ghost_sched%n_recv, &
                                    size(ghost_sched%Pg_recv_buf), ghost_sched%Pg_recv_buf, &
                                    ghost_sched%pg_recv_slot, q)
  end subroutine exchange_Pg_unpack


  !> Device unpack kernel: writes bc_pg(:,3:6,slot) from the received Pg buffer. Each recv
  !> entry writes a distinct slot -> race-free. bc_pg is resident (Mod_Explicit F5 block).
  subroutine exchange_Pg_unpack_buf(bc_pg, nbnd, nr, nbuf, rbuf, rslot, q)
    use MOSE_Global_m, only: nprim
    implicit none
    integer,  intent(in)    :: nbnd, nr, nbuf, q
    real(R8), intent(inout) :: bc_pg(nprim, 6, nbnd)
    real(R8), intent(in)    :: rbuf(nbuf)
    integer,  intent(in)    :: rslot(nr)
    integer :: i, s, slot, base

    !$acc parallel loop gang vector async(q) present(bc_pg, rbuf, rslot) private(s, slot, base)
    do i = 1, nr
      slot = rslot(i)
      base = (i - 1) * nprim * 4
      do s = 1, nprim
        bc_pg(s, 3, slot) = rbuf(base          + s)
        bc_pg(s, 4, slot) = rbuf(base +   nprim + s)
        bc_pg(s, 5, slot) = rbuf(base + 2*nprim + s)
        bc_pg(s, 6, slot) = rbuf(base + 3*nprim + s)
      end do
    end do
  end subroutine exchange_Pg_unpack_buf


  !> Fill the receiver-normal stencil bc_pg(:,1:gc,slot) ON DEVICE from the receiver ghost P
  !> (which exchange_P_unpack has already written into the resident blk(Bm)%P), for every
  !> INTER-rank recv entry. Replaces host Ghost_Interrank + the per-slot update device(bc_pg).
  !> One kernel per local receiver block (single present-matched P). bc_pg(:,1:gc) = the
  !> normal ghost layers, exactly Ghost_Interrank's Pg(:,g) = blkm%P(:,Ig,Jg,Kg).
  subroutine fill_bc_pg_normal_recv(domain, async_q)
    use MOSE_Types_m
    use MOSE_Global_m, only: nprim
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in), optional :: async_q
    integer :: b

    if (ghost_sched%n_recv == 0) return
    if (present(async_q)) then
      block
        use MOSE_Field_m, only: P_flat, ncell_total, dboff, dpi, dpj
        call fill_bc_pg_normal_all(P_flat, ncell_total, dboff, dpi, dpj, size(dboff), &
                                     domain % bc_pg, domain % nbound, ghost_sched%n_recv, &
                                     ghost_sched%recv_b, ghost_sched%recv_id, ghost_sched%recv_jd, &
                                     ghost_sched%recv_kd, ghost_sched%recv_di, ghost_sched%recv_dj, &
                                     ghost_sched%recv_dk, ghost_sched%pg_recv_slot, async_q)
      end block
      return
    end if
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      block
        use MOSE_Field_m, only: P_flat, ncell_total, field_block_index_params
        integer :: boff, pi, pj
        call field_block_index_params(b, boff, pi, pj)
        call fill_bc_pg_normal_blk(P_flat, ncell_total, boff, pi, pj, b, &
                                         domain % bc_pg, domain % nbound, ghost_sched%n_recv, &
                                         ghost_sched%recv_b, ghost_sched%recv_id, ghost_sched%recv_jd, &
                                         ghost_sched%recv_kd, ghost_sched%recv_di, ghost_sched%recv_dj, &
                                         ghost_sched%recv_dk, ghost_sched%pg_recv_slot)
      end block
    end do
  end subroutine fill_bc_pg_normal_recv



  !> The receiver-normal bc_pg slots of one block from the just-unpacked w ghosts.
  subroutine fill_bc_pg_normal_blk(w, ncw, boff, pi, pj, b, bc_pg, nbnd, nr, rb, rid, rjd, rkd, rdi, rdj, rdk, rslot)
    use MOSE_Global_m, only: gc, nprim
    implicit none
    integer,  intent(in)    :: ncw, boff, pi, pj, b, nbnd, nr
    real(R8), intent(in)    :: w(ncw, nprim)
    real(R8), intent(inout) :: bc_pg(nprim, 6, nbnd)
    integer,  intent(in)    :: rb(nr), rid(nr), rjd(nr), rkd(nr), rdi(nr), rdj(nr), rdk(nr), rslot(nr)
    integer :: i, g, s, Id, Jd, Kd, di, dj, dk, slot, ii, jj, kk

    !$acc parallel loop gang vector present(w, bc_pg, rb, rid, rjd, rkd, rdi, rdj, rdk, rslot) &
    !$acc&   private(g, s, Id, Jd, Kd, di, dj, dk, slot, ii, jj, kk)
    do i = 1, nr
      if (rb(i) /= b) cycle
      Id = rid(i); Jd = rjd(i); Kd = rkd(i)
      di = rdi(i); dj = rdj(i); dk = rdk(i)
      slot = rslot(i)
      do g = 1, gc
        ii = Id - di*g; jj = Jd - dj*g; kk = Kd - dk*g
        do s = 1, nprim
          bc_pg(s, g, slot) = WGCELL(s, ii, jj, kk)
        end do
      end do
    end do
  end subroutine fill_bc_pg_normal_blk


  !> mr-fuse-exchange: the four per-block exchange kernels as ONE launch each over every entry of the
  !> schedule. Entry i touches only its own buffer / bc_pg slot and cells of its own block (send_b(i) /
  !> recv_b(i), always local), which the per-block launch owning it touches in the same way; the
  !> kernels copy values and do no arithmetic -> bit-identical, whatever the launch partition.
  !> The block's flat-index scalars come from the resident dboff/dpi/dpj tables.
  subroutine exchange_P_pack_all(w, ncw, dboff, dpi, dpj, nb, ns, nbuf, sbuf, sb, sis, sjs, sks, sdi, sdj, sdk, q)
    use MOSE_Global_m, only: gc, nprim
    implicit none
    integer,  intent(in)    :: ncw, nb, ns, nbuf, q
    real(R8), intent(in)    :: w(ncw, nprim)
    real(R8), intent(inout) :: sbuf(nbuf)
    integer,  intent(in)    :: dboff(nb), dpi(nb), dpj(nb)
    integer,  intent(in)    :: sb(ns), sis(ns), sjs(ns), sks(ns), sdi(ns), sdj(ns), sdk(ns)
    integer :: i, g, s, Is, Js, Ks, di, dj, dk, base, ii, jj, kk, bq, boff, pi, pj

    !$acc parallel loop gang vector async(q) present(w, sbuf, sb, sis, sjs, sks, sdi, sdj, sdk, dboff, dpi, dpj) &
    !$acc&   private(g, s, Is, Js, Ks, di, dj, dk, base, ii, jj, kk, bq, boff, pi, pj)
    do i = 1, ns
      bq = sb(i); boff = dboff(bq); pi = dpi(bq); pj = dpj(bq)
      Is = sis(i); Js = sjs(i); Ks = sks(i)
      di = sdi(i); dj = sdj(i); dk = sdk(i)
      base = (i - 1) * nprim * gc
      do g = 1, gc
        ii = Is + di*(g-1); jj = Js + dj*(g-1); kk = Ks + dk*(g-1)
        do s = 1, nprim
          sbuf(base + (g-1)*nprim + s) = WGCELL(s, ii, jj, kk)
        end do
      end do
    end do
  end subroutine exchange_P_pack_all

  subroutine exchange_P_unpack_all(w, ncw, dboff, dpi, dpj, nb, nr, nbuf, rbuf, rb, rid, rjd, rkd, rdi, rdj, rdk, q)
    use MOSE_Global_m, only: gc, nprim
    implicit none
    integer,  intent(in)    :: ncw, nb, nr, nbuf, q
    real(R8), intent(inout) :: w(ncw, nprim)
    real(R8), intent(in)    :: rbuf(nbuf)
    integer,  intent(in)    :: dboff(nb), dpi(nb), dpj(nb)
    integer,  intent(in)    :: rb(nr), rid(nr), rjd(nr), rkd(nr), rdi(nr), rdj(nr), rdk(nr)
    integer :: i, g, s, Id, Jd, Kd, di, dj, dk, base, ii, jj, kk, bq, boff, pi, pj

    !$acc parallel loop gang vector async(q) present(w, rbuf, rb, rid, rjd, rkd, rdi, rdj, rdk, dboff, dpi, dpj) &
    !$acc&   private(g, s, Id, Jd, Kd, di, dj, dk, base, ii, jj, kk, bq, boff, pi, pj)
    do i = 1, nr
      bq = rb(i); boff = dboff(bq); pi = dpi(bq); pj = dpj(bq)
      Id = rid(i); Jd = rjd(i); Kd = rkd(i)
      di = rdi(i); dj = rdj(i); dk = rdk(i)
      base = (i - 1) * nprim * gc
      do g = 1, gc
        ii = Id - di*g; jj = Jd - dj*g; kk = Kd - dk*g
        do s = 1, nprim
          WGCELL(s, ii, jj, kk) = rbuf(base + (g-1)*nprim + s)
        end do
      end do
    end do
  end subroutine exchange_P_unpack_all

  subroutine exchange_Pg_pack_all(w, ncw, dboff, dpi, dpj, nb, ns, nbuf, sbuf, sb, sis, sjs, sks, sfs, sd11, sd12, sd21, sd22, q)
    use MOSE_Global_m, only: gc, nprim
    implicit none
    integer,  intent(in)    :: ncw, nb, ns, nbuf, q
    real(R8), intent(in)    :: w(ncw, nprim)
    real(R8), intent(inout) :: sbuf(nbuf)
    integer,  intent(in)    :: dboff(nb), dpi(nb), dpj(nb)
    integer,  intent(in)    :: sb(ns), sis(ns), sjs(ns), sks(ns), sfs(ns)
    integer,  intent(in)    :: sd11(ns), sd12(ns), sd21(ns), sd22(ns)
    integer :: i, s, Is, Js, Ks, Fs, d11, d12, d21, d22, base, bq, boff, pi, pj
    integer :: i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4

    !$acc parallel loop gang vector async(q) present(w, sbuf, sb, sis, sjs, sks, sfs, sd11, sd12, sd21, sd22, dboff, dpi, dpj) &
    !$acc&   private(s, Is, Js, Ks, Fs, d11, d12, d21, d22, base, bq, boff, pi, pj, &
    !$acc&           i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4)
    do i = 1, ns
      bq = sb(i); boff = dboff(bq); pi = dpi(bq); pj = dpj(bq)
      Is = sis(i); Js = sjs(i); Ks = sks(i); Fs = sfs(i)
      d11 = sd11(i); d12 = sd12(i); d21 = sd21(i); d22 = sd22(i)
      select case (Fs)
        case (1:2)
          i1=Is; j1=Js-d11; k1=Ks-d12;  i2=Is; j2=Js+d11; k2=Ks+d12
          i3=Is; j3=Js-d21; k3=Ks-d22;  i4=Is; j4=Js+d21; k4=Ks+d22
        case (3:4)
          i1=Is-d11; j1=Js; k1=Ks-d12;  i2=Is+d11; j2=Js; k2=Ks+d12
          i3=Is-d21; j3=Js; k3=Ks-d22;  i4=Is+d21; j4=Js; k4=Ks+d22
        case default  ! (5:6)
          i1=Is-d11; j1=Js-d12; k1=Ks;  i2=Is+d11; j2=Js+d12; k2=Ks
          i3=Is-d21; j3=Js-d22; k3=Ks;  i4=Is+d21; j4=Js+d22; k4=Ks
      end select
      base = (i - 1) * nprim * 4
      do s = 1, nprim
        sbuf(base          + s) = WGCELL(s, i1, j1, k1)   ! -> Pg(:,3)
        sbuf(base +   nprim + s) = WGCELL(s, i2, j2, k2)  ! -> Pg(:,4)
        sbuf(base + 2*nprim + s) = WGCELL(s, i3, j3, k3)  ! -> Pg(:,5)
        sbuf(base + 3*nprim + s) = WGCELL(s, i4, j4, k4)  ! -> Pg(:,6)
      end do
    end do
  end subroutine exchange_Pg_pack_all

  subroutine fill_bc_pg_normal_all(w, ncw, dboff, dpi, dpj, nb, bc_pg, nbnd, nr, rb, rid, rjd, rkd, rdi, rdj, rdk, rslot, q)
    use MOSE_Global_m, only: gc, nprim
    implicit none
    integer,  intent(in)    :: ncw, nb, nbnd, nr, q
    real(R8), intent(in)    :: w(ncw, nprim)
    real(R8), intent(inout) :: bc_pg(nprim, 6, nbnd)
    integer,  intent(in)    :: dboff(nb), dpi(nb), dpj(nb)
    integer,  intent(in)    :: rb(nr), rid(nr), rjd(nr), rkd(nr), rdi(nr), rdj(nr), rdk(nr), rslot(nr)
    integer :: i, g, s, Id, Jd, Kd, di, dj, dk, slot, ii, jj, kk, bq, boff, pi, pj

    !$acc parallel loop gang vector async(q) present(w, bc_pg, rb, rid, rjd, rkd, rdi, rdj, rdk, rslot, dboff, dpi, dpj) &
    !$acc&   private(g, s, Id, Jd, Kd, di, dj, dk, slot, ii, jj, kk, bq, boff, pi, pj)
    do i = 1, nr
      bq = rb(i); boff = dboff(bq); pi = dpi(bq); pj = dpj(bq)
      Id = rid(i); Jd = rjd(i); Kd = rkd(i)
      di = rdi(i); dj = rdj(i); dk = rdk(i)
      slot = rslot(i)
      do g = 1, gc
        ii = Id - di*g; jj = Jd - dj*g; kk = Kd - dk*g
        do s = 1, nprim
          bc_pg(s, g, slot) = WGCELL(s, ii, jj, kk)
        end do
      end do
    end do
  end subroutine fill_bc_pg_normal_all


  !> Host-staging bridges for the Pg buffers (CUDA-aware-MPI FALLBACK / device-bound host path).
  subroutine exchange_Pg_send_buf_to_host(async_q)
    use openacc, only: acc_async_sync
    implicit none
    integer, intent(in), optional :: async_q   !< T1c: the stage queue the copy rides (absent = synchronous, as before)
    integer :: q
    if (ghost_sched%n_send == 0) return
    if (.not. ghost_sched%pg_bufs_resident) return
    q = acc_async_sync
    if (present(async_q)) q = async_q
    !$acc update host(ghost_sched%Pg_send_buf) async(q)
  end subroutine exchange_Pg_send_buf_to_host

  subroutine exchange_Pg_recv_buf_to_device(async_q)
    use openacc, only: acc_async_sync
    implicit none
    integer, intent(in), optional :: async_q   !< T1c: the stage queue the copy rides (absent = synchronous, as before)
    integer :: q
    if (ghost_sched%n_recv == 0) return
    if (.not. ghost_sched%pg_bufs_resident) return
    q = acc_async_sync
    if (present(async_q)) q = async_q
    !$acc update device(ghost_sched%Pg_recv_buf) async(q)
  end subroutine exchange_Pg_recv_buf_to_device



#endif


  !> MPI barrier for synchronizing ranks after I/O operations.
  !> Prevents non-root ranks from advancing to the next MPI collective
  !> while root is still writing files.
  subroutine mpi_io_barrier()
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call check_mpi_error(ierr)
    end block
#endif
  end subroutine mpi_io_barrier


  !>  device-Pg exchange orchestration (mirror of the persistent P exchange wrappers).
  !> Public + unconditional (depth-0, like exchange_ghost_P_pack / mpi_io_barrier; bodies
  !> USE_MPI/_OPENACC-guarded) so the names resolve in the non-MPI / non-OPENACC build. Call
  !> all from !$omp single. With CUDA-aware MPI the persistent requests are device-bound, so
  !> MPI DMAs the resident Pg buffers directly; otherwise the host-staging bridges run.

  subroutine exchange_ghost_Pg_post_recv()
    implicit none
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr
      if (.not. ghost_sched%persistent_Pg_init) return
      if (ghost_sched%n_recv_ranks > 0) then
        call MPI_STARTALL(ghost_sched%n_recv_ranks, ghost_sched%Pg_recv_req_pers, ierr)
        call check_mpi_error(ierr)
      end if
    end block
#endif
  end subroutine exchange_ghost_Pg_post_recv


  !> Pack the device Pg send buffer (out of the resident donor P) ready for the persistent
  !> sends. On the fallback (host-bound requests) the packed device buffer is staged to host.
  subroutine exchange_ghost_Pg_pack(domain, async_q)
    use MOSE_Types_m
    use openacc, only: acc_async_sync
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer, intent(in), optional :: async_q   !< mr-fuse-exchange: one launch on this queue, then ONE wait on it
    if (mpi_size_ <= 1) return
#if defined(USE_MPI)
    call exchange_Pg_pack(domain, async_q)
    ! Host-bound requests: the device buffer must reach the host before MPI reads it. T1c: on the pack's queue, before the one
    ! wait below, which then covers the pack and the copy.
    if (.not. (ghost_sched%dev_bound)) &
        call exchange_Pg_send_buf_to_host(async_q)
    ! unconditional join (also with nothing to send): the caller's local fill ran on this queue
    if (present(async_q)) then
      if (async_q /= acc_async_sync) then
        !$acc wait(async_q)
      end if
    end if
#endif
  end subroutine exchange_ghost_Pg_pack


  subroutine exchange_ghost_Pg_post_send()
    implicit none
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr
      if (.not. ghost_sched%persistent_Pg_init) return
      if (ghost_sched%n_send_ranks > 0) then
        call MPI_STARTALL(ghost_sched%n_send_ranks, ghost_sched%Pg_send_req_pers, ierr)
        call check_mpi_error(ierr)
      end if
    end block
#endif
  end subroutine exchange_ghost_Pg_post_send


  !> Wait for the Pg receives, then (fallback) stage the received buffer to device. After this
  !> the device unpack kernel can read Pg_recv_buf.
  subroutine exchange_ghost_Pg_wait_recv()
    implicit none
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr
      if (.not. ghost_sched%persistent_Pg_init) return
      if (ghost_sched%n_recv_ranks > 0) then
        call MPI_WAITALL(ghost_sched%n_recv_ranks, ghost_sched%Pg_recv_req_pers, &
                         ghost_sched%mpi_stat(:,1:ghost_sched%n_recv_ranks), ierr)
        call check_mpi_error(ierr)
      end if
    end block
#endif
    ! T1c: the host-to-device copy of a host-bound Pg receive moved to exchange_ghost_Pg_unpack, on the unpack's queue
    ! (the only caller of this routine, Fill_Ghost_Cell_mr_finish, calls the unpack right after it)
  end subroutine exchange_ghost_Pg_wait_recv


  !> Public unconditional wrapper: device-unpack the received transverse stencil into the
  !> resident bc_pg(:,3:6). Body guarded so the name resolves in the non-OPENACC / non-MPI build.
  subroutine exchange_ghost_Pg_unpack(domain, async_q)
    use MOSE_Types_m
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in), optional :: async_q
    if (mpi_size_ <= 1) return
#if defined(USE_MPI)
    if (.not. (ghost_sched%dev_bound)) &
        call exchange_Pg_recv_buf_to_device(async_q)
    call exchange_Pg_unpack(domain, async_q)
#endif
  end subroutine exchange_ghost_Pg_unpack


  !> Public unconditional wrapper: device-fill the receiver-normal bc_pg(:,1:gc) from the
  !> (already device-unpacked) receiver ghost P. Body guarded for the non-OPENACC / non-MPI build.
  subroutine fill_bc_pg_normal(domain, async_q)
    use MOSE_Types_m
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in), optional :: async_q
    if (mpi_size_ <= 1) return
#if defined(USE_MPI)
    call fill_bc_pg_normal_recv(domain, async_q)
#endif
  end subroutine fill_bc_pg_normal


  subroutine exchange_ghost_Pg_wait_send()
    implicit none
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr
      if (.not. ghost_sched%persistent_Pg_init) return
      if (ghost_sched%n_send_ranks > 0) then
        call MPI_WAITALL(ghost_sched%n_send_ranks, ghost_sched%Pg_send_req_pers, &
                         ghost_sched%mpi_stat(:,1:ghost_sched%n_send_ranks), ierr)
        call check_mpi_error(ierr)
      end if
    end block
#endif
  end subroutine exchange_ghost_Pg_wait_send



end module MOSE_Mod_GhostExchange
