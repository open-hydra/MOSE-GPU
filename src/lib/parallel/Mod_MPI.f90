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

  ! --- Block-to-rank mapping ---
  integer, allocatable, public :: block_owner(:)     !< block_owner(b) = rank owning block b
  integer, allocatable, public :: local_block_ids(:) !< Global block IDs owned by this rank
  integer, public :: n_local_blocks = 0

  ! --- Public procedures ---
  public :: mpi_init_env, mpi_finalize_env
  public :: mpi_bind_device
  public :: is_local_block
  public :: partition_blocks
  public :: mpi_allreduce_min_r8

  public :: mpi_reduce_sum_r8_array

  public :: mpi_bcast_logical, mpi_bcast_integer
  public :: mpi_abort_all
  public :: check_mpi_error

contains


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
    integer :: b, i, r, nloc
    integer, allocatable :: rank_load(:)
    integer, allocatable :: order(:)
    integer :: tmp
    integer(kind=8) :: lmin, lmax, lsum

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

    ! Build local block list
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

    ! Final balance check (root prints summary)
    if (mpi_is_root) then
      lmin = minval(rank_load)
      lmax = maxval(rank_load)
      lsum = sum(int(rank_load, kind=8))
      write(*,'(A)') ' Block partition (ascending size):'
      write(*,'(A,I8,A,I8)')   '  - cells per rank: min = ', lmin, '  max = ', lmax
      write(*,'(A,F8.3)')      '  - imbalance (max/avg) = ', &
        real(lmax, kind=8) / max(real(lsum, kind=8) / real(mpi_size_, kind=8), 1.0_8)
    end if

    deallocate(rank_load)
    if (allocated(order)) deallocate(order)   ! not allocated when the map came from a file
  end subroutine partition_blocks


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
