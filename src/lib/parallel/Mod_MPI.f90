!>@brief MPI infrastructure module for MOSE.
!> Provides MPI environment management, block partitioning, and
!> global reduction wrappers. When compiled without USE_MPI, all
!> routines are no-ops or identity operations (serial fallback).
module MOSE_Mod_MPI
#ifdef USE_MPI
  use mpi
#endif
  use iso_fortran_env, only: R8 => real64, I4 => int32

  implicit none
  private

  ! --- Public state ---
  integer, public :: mpi_rank_ = 0        !< This process rank (0 in serial)
  integer, public :: mpi_size_ = 1        !< Total number of MPI processes
  logical, public :: mpi_is_root = .true.  !< True on rank 0
  real(R8), public :: halo_wait_time = 0.0_R8 !< accumulated halo MPI-wait seconds (this rank)
  !> The halo exchange split by phase, seconds accumulated on this rank (Lib_Extrapolation Fill_Ghost_Cell_mr_start/finish):
  !> 1 post the P receives, 2 pack P on device + post the sends, 3 launch the local fills (async queue), 4 Pg: post receives,
  !> pack (joins the local fills), post sends, 5 wait the P receives, 6 unpack P + chimera + receiver-normal bc_pg,
  !> 7 wait the Pg receives, 8 unpack Pg + wait the sends. 5 + 7 + part of 8 = halo_wait_time.
  real(R8), public :: halo_t(8) = 0.0_R8
  !> The same phases per multirate level (halo_level = the level of the stage being filled, -1 outside a stage):
  !> separates the fixed cost of the restricted exchange of the fine levels from their load imbalance.
  real(R8), public :: halo_tl(8, 0:8) = 0.0_R8
  integer,  public :: halo_level = -1

  ! --- Block-to-rank mapping ---
  integer, allocatable, public :: block_owner(:)     !< block_owner(b) = rank owning block b
  integer, allocatable, public :: local_block_ids(:) !< Global block IDs owned by this rank
  integer, public :: n_local_blocks = 0
  logical, public :: block_owner_from_file = .false. !< INPUT/block_owner.txt supplied the map: the connectivity pass is skipped

  ! --- Public procedures ---
  public :: mpi_init_env, mpi_finalize_env
  public :: halo_acc
  public :: mpi_bind_device
  public :: is_local_block
  public :: partition_blocks, partition_blocks_connectivity
  public :: mpi_allreduce_min_r8

  public :: mpi_reduce_sum_r8_array
  public :: mpi_allreduce_sum_r8_array

  public :: mpi_bcast_logical, mpi_bcast_integer
  public :: mpi_abort_all
  public :: check_mpi_error

contains

  !> Add a halo phase time to the run total and to the current multirate level.
  subroutine halo_acc ( k, v )
    integer,  intent(in) :: k
    real(R8), intent(in) :: v
    halo_t(k) = halo_t(k) + v
    if ( halo_level >= 0 .and. halo_level <= 8 ) halo_tl(k, halo_level) = halo_tl(k, halo_level) + v
  end subroutine halo_acc


  !> Initialize MPI environment with thread support. Call early in main program.
  !> Uses MPI_THREAD_FUNNELED: only master thread makes MPI calls.
  subroutine mpi_init_env()
#ifdef USE_MPI
    integer :: ierr, provided
    call MPI_INIT_THREAD(MPI_THREAD_FUNNELED, provided, ierr)
    call check_mpi_error(ierr)
    if (provided < MPI_THREAD_FUNNELED) then
      write(*,'(A)') ' ERROR: MPI does not support the required threading level (MPI_THREAD_FUNNELED)'
      call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank_, ierr)
    call check_mpi_error(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size_, ierr)
    call check_mpi_error(ierr)
    mpi_is_root = (mpi_rank_ == 0)
#else
    mpi_rank_ = 0
    mpi_size_ = 1
    mpi_is_root = .true.
#endif
    call mpi_bind_device()
  end subroutine mpi_init_env


  !> Bind this rank to a node-local GPU, round-robin over the visible devices.
  !> Without this, every rank of an unwrapped multi-GPU launch lands on device 0;
  !> it composes unchanged with a CUDA_VISIBLE_DEVICES wrapper (one visible
  !> device -> mod(lrank,1) = 0 selects it). The current device is a per-thread
  !> OpenACC ICV and device work is launched from inside omp regions, so the
  !> whole thread team binds. MPI calls stay on the initial thread (FUNNELED).
  subroutine mpi_bind_device()
    use openacc
    integer :: ndev, lrank, dev
#ifdef USE_MPI
    integer :: ierr, node_comm
    call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, mpi_rank_, MPI_INFO_NULL, node_comm, ierr)
    call check_mpi_error(ierr)
    call MPI_COMM_RANK(node_comm, lrank, ierr)
    call check_mpi_error(ierr)
    call MPI_COMM_FREE(node_comm, ierr)
#else
    lrank = 0
#endif
    ndev = acc_get_num_devices(acc_device_nvidia)
    if (ndev > 0) then
      dev = mod(lrank, ndev)
      !$omp parallel
      call acc_set_device_num(dev, acc_device_nvidia)
      !$omp end parallel
      write(*,'(A,I4,A,I2,A,I2)') ' [MPI] rank ', mpi_rank_, ' -> GPU ', dev, ' of ', ndev
    end if
  end subroutine mpi_bind_device


  !> Finalize MPI environment. Call at end of main program.
  subroutine mpi_finalize_env()
#ifdef USE_MPI
    integer :: ierr
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_FINALIZE(ierr)
#endif
  end subroutine mpi_finalize_env


  !> Check if block b is owned by this rank.
  logical function is_local_block(b)
    integer, intent(in) :: b

    if (.not. allocated(block_owner)) then
      is_local_block = .true.  ! Serial mode or not yet partitioned
      return
    end if
    is_local_block = (block_owner(b) == mpi_rank_)
  end function is_local_block


  !> Partition nb blocks across MPI ranks.
  !> Strategy: LPT (Longest Processing Time first) -- sort blocks by DESCENDING
  !> cell count and assign each to the rank with the least current load.
  !> (Ascending order was the original strategy; it is the worst greedy variant:
  !> the big blocks land last on evenly-filled ranks, creating the n=3/n=6
  !> imbalance measured 2026-06-26. LPT bounds the makespan at 4/3 OPT.)
  !> Ownership only -- each block's arithmetic is rank-independent, so the field
  !> is bit-identical under any map; only halo pairing/residual-reduce order move.
  !> blk_ncells(b) = total number of cells in block b.
  subroutine partition_blocks(nb, blk_ncells)
    integer, intent(in) :: nb
    integer, intent(in) :: blk_ncells(nb)
    ! Local
    integer :: b, i, r
    integer, allocatable :: rank_load(:)
    integer, allocatable :: order(:)
    integer :: tmp

    if (allocated(block_owner)) deallocate(block_owner)
    if (allocated(rank_load)) deallocate(rank_load)
    allocate(block_owner(nb))
    allocate(rank_load(0:mpi_size_-1))
    rank_load = 0

    ! Optional explicit map: INPUT/block_owner.txt, one 0-based rank per block, in
    ! block order. Used when the physics needs a specific co-location that LPT
    ! cannot know about -- e.g. chimera donors must be local to their receptor on
    ! the device path, so an annular RDE keeps each wedge's blocks on one rank.
    if (read_block_owner_file(nb, rank_load, blk_ncells)) then
      if (mpi_is_root) write(*,'(A)') ' Block partition read from INPUT/block_owner.txt'
      block_owner_from_file = .true.
    else

    ! Build an order index sorted by ascending cell count (insertion sort: nb is small).
    allocate(order(nb))
    do b = 1, nb
      order(b) = b
    end do
    do i = 2, nb
      tmp = order(i)
      b = i - 1
      do while (b >= 1)
        if (blk_ncells(order(b)) >= blk_ncells(tmp)) exit   ! descending (LPT)
        order(b+1) = order(b)
        b = b - 1
      end do
      order(b+1) = tmp
    end do

    ! Walk in size-DESCENDING order (LPT), assigning each block to the least-loaded rank.
    do i = 1, nb
      b = order(i)
      r = minloc(rank_load, dim=1) - 1   ! rank with min load (0-indexed)
      block_owner(b) = r
      rank_load(r) = rank_load(r) + blk_ncells(b)
    end do

    end if

    call rebuild_local_block_list(nb, blk_ncells, ' Block partition (LPT by block size):')

    deallocate(rank_load)
    if (allocated(order)) deallocate(order)   ! not allocated when the map came from a file
  end subroutine partition_blocks


  !> Rebuild local_block_ids / n_local_blocks from block_owner; root prints the load balance.
  subroutine rebuild_local_block_list(nb, blk_ncells, label)
    integer, intent(in)          :: nb
    integer, intent(in)          :: blk_ncells(nb)
    character(len=*), intent(in) :: label
    integer :: b, nloc
    integer(kind=8) :: lmin, lmax, lsum
    integer(kind=8), allocatable :: rank_load(:)

    nloc = count(block_owner == mpi_rank_)
    n_local_blocks = nloc
    if (allocated(local_block_ids)) deallocate(local_block_ids)
    allocate(local_block_ids(nloc))
    nloc = 0
    do b = 1, nb
      if (block_owner(b) == mpi_rank_) then
        nloc = nloc + 1
        local_block_ids(nloc) = b
      end if
    end do

    if (mpi_is_root) then
      allocate(rank_load(0:mpi_size_-1))
      rank_load = 0
      do b = 1, nb
        rank_load(block_owner(b)) = rank_load(block_owner(b)) + int(blk_ncells(b), kind=8)
      end do
      lmin = minval(rank_load)
      lmax = maxval(rank_load)
      lsum = sum(rank_load)
      write(*,'(A)') label
      write(*,'(A,I8,A,I8)')   '  - cells per rank: min = ', lmin, '  max = ', lmax
      write(*,'(A,F8.3)')      '  - imbalance (max/avg) = ', &
        real(lmax, kind=8) / max(real(lsum, kind=8) / real(mpi_size_, kind=8), 1.0_8)
      deallocate(rank_load)
    end if
  end subroutine rebuild_local_block_list


  !> Second pass of the block partition, once the connection records are known (called from
  !> Setup_BC): assign the blocks to the ranks by contiguity along the connection faces, so
  !> that the halo exchange moves the cells of a few cut faces instead of nearly every
  !> connection (LPT by size scatters neighbouring blocks across the ranks: on the 60-block
  !> nozzle ring 92 % of the connection cells crossed ranks and no block was interior).
  !> Search: for each seed block (every stride-th block), the other seeds are farthest-first
  !> in connection hops; each rank grows from its seed taking the unassigned block best
  !> connected to it (least-loaded rank first); overloaded ranks give blocks away; then
  !> blocks move across ranks while the inter-rank connection cells decrease and every rank
  !> stays within tol of the average load (tol = 1.02, MOSE_PART_TOL). The best map within
  !> the tolerance replaces the LPT map when it cuts fewer connection cells. Deterministic
  !> (no random seeds; ties go to the lowest block index), so the map, the halo counts and
  !> the timings are reproducible. Ownership only: the field is bit-identical under any map.
  !> Skipped on 1 rank, with an INPUT/block_owner.txt map, or with MOSE_PART=lpt.
  !> wd(b, bs) = ghost cells of block b fed by block bs (type-1 connection records).
  subroutine partition_blocks_connectivity(nb, blk_ncells, wd)
    integer, intent(in) :: nb
    integer, intent(in) :: blk_ncells(nb)
    integer, intent(in) :: wd(nb, nb)
    ! Local
    integer :: nr, b, c, r, r0, t, i, pass, ntrial, nleft, rbest, gbest, g, nblk, nint, nrep
    integer :: cut_lpt, cut_best, cut_t, wtot
    integer(kind=8) :: lcg
    integer, allocatable :: w(:,:), conn(:,:), own(:), best(:), seed(:), dist(:), order(:)
    integer(kind=8), allocatable :: load(:), recv(:), cint(:)
    integer(kind=8) :: tot, cap
    real(kind=8) :: tol
    logical :: improved
    character(len=16) :: envval
    integer :: envlen, envstat

    nr = mpi_size_
    if (nr <= 1 .or. nb < nr) return
    if (block_owner_from_file) return
    call get_environment_variable('MOSE_PART', envval, envlen, envstat)
    if (envstat == 0 .and. envlen >= 3) then
      if (envval(1:3) == 'lpt' .or. envval(1:3) == 'LPT') then
        if (mpi_is_root) write(*,'(A)') ' Block partition: MOSE_PART=lpt -> the connectivity pass is skipped'
        return
      end if
    end if
    tol = 1.02d0
    call get_environment_variable('MOSE_PART_TOL', envval, envlen, envstat)
    if (envstat == 0 .and. envlen > 0) then
      read(envval(1:envlen), *, iostat=envstat) tol
      if (envstat /= 0 .or. tol < 1.0d0) tol = 1.02d0
    end if

    allocate(w(nb, nb), conn(nb, 0:nr-1), own(nb), best(nb), seed(nr), dist(nb), order(nr))
    allocate(load(0:nr-1), recv(0:nr-1), cint(0:nr-1))
    w = wd + transpose(wd)
    do b = 1, nb
      w(b, b) = 0
    end do
    wtot = sum(wd)
    tot = sum(int(blk_ncells, kind=8))
    cap = int(tol * real(tot, kind=8) / real(nr, kind=8), kind=8)
    best = block_owner
    cut_lpt = cut_cells(block_owner)
    cut_best = cut_lpt

    ! Multi-start: the first nb trials seed rank 0 at block t with the other seeds farthest-first in
    ! connection hops; the remaining trials draw their seeds from a fixed-seed linear congruential
    ! generator (reproducible: the same mesh and rank count give the same map on every run). The
    ! trial count is bounded by the work per trial (~nb^2).
    ntrial = min(50000, max(64, int(2.0d9 / (real(nb,8)**2))))
    lcg = 20261105_8
    do t = 1, ntrial
      if (t <= nb) then
        seed(1) = t
        do r = 2, nr
          call hops_from(seed, r-1, dist)
          seed(r) = 0
          do b = 1, nb
            if (any(seed(1:r-1) == b)) cycle
            if (seed(r) == 0) then
              seed(r) = b
            else if (dist(b) > dist(seed(r))) then
              seed(r) = b
            end if
          end do
        end do
      else
        do r = 1, nr
          do
            lcg = mod(1103515245_8 * lcg + 12345_8, 2147483648_8)
            b = int(mod(lcg / 65536_8, int(nb, 8))) + 1
            if (.not. any(seed(1:r-1) == b)) exit
          end do
          seed(r) = b
        end do
      end if
      ! Growth: the least-loaded rank takes the unassigned block best connected to it
      own = -1; load = 0; conn = 0
      do r = 1, nr
        call assign(seed(r), r-1)
      end do
      nleft = nb - nr
      do while (nleft > 0)
        call ranks_by_load(order)
        c = 0
        do i = 1, nr
          r = order(i)
          gbest = 0
          do b = 1, nb
            if (own(b) >= 0) cycle
            if (conn(b, r) > gbest) then
              gbest = conn(b, r); c = b
            end if
          end do
          if (c > 0) exit
        end do
        if (c == 0) then   ! no unassigned block touches an assigned one: lowest index to the least-loaded rank
          r = order(1)
          do b = 1, nb
            if (own(b) < 0) then
              c = b; exit
            end if
          end do
        end if
        call assign(c, r)
        nleft = nleft - 1
      end do
      ! Repair: an overloaded rank gives its least-attached block to the least-loaded rank with room
      do nrep = 1, nb
        r0 = maxloc(load, dim=1) - 1
        if (load(r0) <= cap) exit
        call ranks_by_load(order)
        r = order(1)
        if (r == r0) exit
        c = 0; gbest = -huge(0)
        do b = 1, nb
          if (own(b) /= r0) cycle
          if (load(r) + blk_ncells(b) > cap) cycle
          if (count(own == r0) <= 1) exit
          g = conn(b, r) - conn(b, r0)
          if (g > gbest) then
            gbest = g; c = b
          end if
        end do
        if (c == 0) exit
        call move(c, r0, r)
      end do
      ! Refinement: move a block to the rank it is best connected to, within the load cap
      do pass = 1, 50
        improved = .false.
        do b = 1, nb
          r0 = own(b)
          if (count(own == r0) <= 1) cycle
          rbest = -1; gbest = 0
          do r = 0, nr-1
            if (r == r0) cycle
            if (load(r) + blk_ncells(b) > cap) cycle
            g = conn(b, r) - conn(b, r0)
            if (g > gbest) then
              gbest = g; rbest = r
            end if
          end do
          if (rbest >= 0) then
            call move(b, r0, rbest)
            improved = .true.
          end if
        end do
        if (.not. improved) exit
      end do
      if (maxval(load) <= cap) then
        cut_t = cut_cells(own)
        if (cut_t < cut_best) then
          cut_best = cut_t; best = own
        end if
      end if
    end do

    if (cut_best < cut_lpt) then
      block_owner = best
      if (mpi_is_root) write(*,'(A,I0,A,I0,A,F5.1,A,I0,A,F5.1,A)') &
        ' Block partition by connectivity: inter-rank connection cells ', cut_best, ' of ', wtot, &
        ' (', 100.0d0*real(cut_best,8)/real(max(wtot,1),8), ' %); LPT by size had ', cut_lpt, &
        ' (', 100.0d0*real(cut_lpt,8)/real(max(wtot,1),8), ' %)'
      call rebuild_local_block_list(nb, blk_ncells, ' Block partition (contiguous along the connection faces):')
    else
      if (mpi_is_root) write(*,'(A,I0,A,I0,A)') ' Block partition by connectivity: no contiguous map within the load tolerance', &
        ' beats LPT by size (inter-rank connection cells ', cut_lpt, ' of ', wtot, '): LPT kept'
    end if

    ! The exchange load of each rank under the final map: cells received across ranks, blocks
    ! without any inter-rank face (the interior batch an exchange/compute overlap can run first).
    if (mpi_is_root) then
      recv = 0; cint = 0
      do r = 0, nr-1
        nblk = 0; nint = 0
        do b = 1, nb
          if (block_owner(b) /= r) cycle
          nblk = nblk + 1
          improved = .true.   ! reused as "interior" flag
          do c = 1, nb
            if (block_owner(c) == r) cycle
            recv(r) = recv(r) + wd(b, c)
            if (w(b, c) > 0) improved = .false.
          end do
          if (improved) then
            nint = nint + 1; cint(r) = cint(r) + blk_ncells(b)
          end if
        end do
        load(r) = 0
        do b = 1, nb
          if (block_owner(b) == r) load(r) = load(r) + blk_ncells(b)
        end do
        write(*,'(A,I4,A,I0,A,I0,A,I0,A,F5.1,A)') '  - rank ', r, ': blocks ', nblk, &
          ', connection cells received from other ranks ', recv(r), &
          ', blocks without inter-rank faces ', nint, ' (', 100.0d0*real(cint(r),8)/real(max(load(r),1_8),8), ' % of its cells)'
      end do
    end if

    deallocate(w, conn, own, best, seed, dist, order, load, recv, cint)

  contains

    !> Directed inter-rank connection cells of a map (= the sum of the receive counts).
    integer function cut_cells(map)
      integer, intent(in) :: map(nb)
      integer :: bb, cc
      cut_cells = 0
      do cc = 1, nb
        do bb = 1, nb
          if (map(bb) /= map(cc)) cut_cells = cut_cells + wd(bb, cc)
        end do
      end do
    end function cut_cells

    !> Multi-source breadth-first hop distance over the connection graph (nb+1 = unreachable).
    subroutine hops_from(src, ns, d)
      integer, intent(in)  :: src(:), ns
      integer, intent(out) :: d(nb)
      integer :: qu(nb), head, tail, bb, cc
      d = nb + 1
      head = 1; tail = 0
      do bb = 1, ns
        d(src(bb)) = 0
        tail = tail + 1; qu(tail) = src(bb)
      end do
      do while (head <= tail)
        bb = qu(head); head = head + 1
        do cc = 1, nb
          if (w(bb, cc) > 0 .and. d(cc) > d(bb) + 1) then
            d(cc) = d(bb) + 1
            tail = tail + 1; qu(tail) = cc
          end if
        end do
      end do
    end subroutine hops_from

    !> Ranks by ascending load (ties: lowest rank).
    subroutine ranks_by_load(ord)
      integer, intent(out) :: ord(nr)
      integer :: ii, jj, tmp
      do ii = 1, nr
        ord(ii) = ii - 1
      end do
      do ii = 2, nr
        tmp = ord(ii); jj = ii - 1
        do while (jj >= 1)
          if (load(ord(jj)) <= load(tmp)) exit
          ord(jj+1) = ord(jj); jj = jj - 1
        end do
        ord(jj+1) = tmp
      end do
    end subroutine ranks_by_load

    subroutine assign(bb, rr)
      integer, intent(in) :: bb, rr
      own(bb) = rr
      load(rr) = load(rr) + blk_ncells(bb)
      conn(:, rr) = conn(:, rr) + w(:, bb)
    end subroutine assign

    subroutine move(bb, rfrom, rto)
      integer, intent(in) :: bb, rfrom, rto
      own(bb) = rto
      load(rfrom) = load(rfrom) - blk_ncells(bb)
      load(rto)   = load(rto)   + blk_ncells(bb)
      conn(:, rfrom) = conn(:, rfrom) - w(:, bb)
      conn(:, rto)   = conn(:, rto)   + w(:, bb)
    end subroutine move

  end subroutine partition_blocks_connectivity


  !> Read INPUT/block_owner.txt (nb lines, 0-based rank per block) into block_owner and
  !> accumulate rank_load. Returns .false. (and touches nothing) if the file is absent;
  !> stops on a malformed file, an out-of-range rank, or a rank left without blocks.
  logical function read_block_owner_file(nb, rank_load, blk_ncells) result(found)
    integer, intent(in)    :: nb, blk_ncells(nb)
    integer, intent(inout) :: rank_load(0:)
    integer :: u, b, r, ios
    logical :: exists

    found = .false.
    inquire(file='INPUT/block_owner.txt', exist=exists)
    if (.not. exists) return
    open(newunit=u, file='INPUT/block_owner.txt', status='old', action='read', iostat=ios)
    if (ios /= 0) return
    do b = 1, nb
      read(u, *, iostat=ios) r
      if (ios /= 0) then
        write(*,*) '[ERROR] INPUT/block_owner.txt: expected ', nb, ' ranks, failed at block ', b
        error stop 'block_owner.txt malformed'
      end if
      if (r < 0 .or. r > mpi_size_-1) then
        write(*,*) '[ERROR] INPUT/block_owner.txt: block ', b, ' -> rank ', r, ' outside 0..', mpi_size_-1
        error stop 'block_owner.txt rank out of range'
      end if
      block_owner(b) = r
      rank_load(r) = rank_load(r) + blk_ncells(b)
    end do
    close(u)
    do r = 0, mpi_size_-1
      if (rank_load(r) == 0) then
        write(*,*) '[ERROR] INPUT/block_owner.txt leaves rank ', r, ' without blocks'
        error stop 'block_owner.txt empty rank'
      end if
    end do
    found = .true.
  end function read_block_owner_file



  !> MPI_ALLREDUCE with MPI_MIN for a scalar real(R8).
  subroutine mpi_allreduce_min_r8(local_val, global_val)
    real(R8), intent(in)  :: local_val
    real(R8), intent(out) :: global_val
#ifdef USE_MPI
    integer :: ierr
    call MPI_ALLREDUCE(local_val, global_val, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#else
    global_val = local_val
#endif
  end subroutine mpi_allreduce_min_r8




  !> MPI_REDUCE with MPI_SUM for an array of real(R8) to root rank 0.
  !> Result is only valid on root. More efficient than ALLREDUCE when
  !> only root needs the result (e.g. for residual output).
  subroutine mpi_reduce_sum_r8_array(arr, n)
    integer, intent(in)     :: n
    real(R8), intent(inout) :: arr(n)
#ifdef USE_MPI
    integer :: ierr
    if (mpi_rank_ == 0) then
      call MPI_REDUCE(MPI_IN_PLACE, arr, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    else
      call MPI_REDUCE(arr, arr, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    end if
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_reduce_sum_r8_array


  !> MPI_ALLREDUCE with MPI_SUM for an array of real(R8), in place: the sum on every rank.
  subroutine mpi_allreduce_sum_r8_array(arr, n)
    integer, intent(in)     :: n
    real(R8), intent(inout) :: arr(n)
#ifdef USE_MPI
    integer :: ierr
    call MPI_ALLREDUCE(MPI_IN_PLACE, arr, n, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_allreduce_sum_r8_array


  !> Broadcast a logical value from root (rank 0) to all ranks.
  subroutine mpi_bcast_logical(val)
    logical, intent(inout) :: val
#ifdef USE_MPI
    integer :: ierr
    call MPI_BCAST(val, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_bcast_logical


  !> Broadcast an integer value from root (rank 0) to all ranks.
  subroutine mpi_bcast_integer(val)
    integer, intent(inout) :: val
#ifdef USE_MPI
    integer :: ierr
    call MPI_BCAST(val, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_bcast_integer


  !> Check MPI return code. Abort all ranks on error.
  subroutine check_mpi_error(ierr)
    integer, intent(in) :: ierr
#ifdef USE_MPI
    integer :: abort_ierr
    if (ierr /= 0) then
      write(*,'(A,I4,A,I6)') ' [RANK ', mpi_rank_, '] MPI error code: ', ierr
      call MPI_ABORT(MPI_COMM_WORLD, 1, abort_ierr)
    end if
#endif
  end subroutine check_mpi_error


  !> MPI-aware abort: prints message, then aborts all ranks.
  subroutine mpi_abort_all(message)
    character(len=*), intent(in) :: message
#ifdef USE_MPI
    integer :: ierr
    write(*,'(A,I4,A,A)') ' [RANK ', mpi_rank_, '] ABORT: ', trim(message)
    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#else
    write(*,'(A,A)') ' ABORT: ', trim(message)
    error stop
#endif
  end subroutine mpi_abort_all


end module MOSE_Mod_MPI
