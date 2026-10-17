!>@brief Automatic block splitting for MPI load balancing.
!> When the number of blocks is less than the number of MPI ranks,
!> this module splits the largest blocks to ensure every rank has work.
module MOSE_Mod_BlockSplit
  use iso_fortran_env, only: R8 => real64
  use MOSE_Mod_MPI

  implicit none
  private
  public :: split_domain_for_mpi
  public :: rebuild_IOfield_after_split

contains


  !> Main entry point: iteratively split the largest block until nb >= mpi_size_.
  subroutine split_domain_for_mpi(domain)
    use MOSE_Types_m
    use MOSE_Global_m, only: nprim

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, b_max, d_max, ncells_max, ncells, nb_orig
    integer :: dims(3)

    ! No split needed if we already have enough blocks
    if (domain%nb >= mpi_size_) return

    nb_orig = domain%nb
    if (mpi_is_root) write(*,'(A,I0,A,I0,A)') &
      ' Block splitting: ', domain%nb, ' blocks for ', mpi_size_, ' MPI ranks'

    do while (domain%nb < mpi_size_)

      ! Find the largest block (by cell count)
      b_max = 0
      ncells_max = 0
      do b = 1, domain%nb
        ncells = domain%blk(b)%dim(1) * domain%blk(b)%dim(2) * domain%blk(b)%dim(3)
        if (ncells > ncells_max) then
          ncells_max = ncells
          b_max = b
        end if
      end do

      ! Find the largest dimension of this block
      dims = domain%blk(b_max)%dim
      d_max = maxloc(dims, dim=1)

      ! Safety: cannot split a dimension smaller than 4 (need >=2 cells per sub-block)
      if (dims(d_max) < 4) then
        if (mpi_is_root) write(*,'(A)') &
          '  WARNING: cannot split further, block dimensions too small'
        exit
      end if

      ! Perform the split
      call split_block(domain, b_max, d_max)

    end do

    if (mpi_is_root) then
      write(*,'(A,I0,A,I0,A)') '  Split complete: ', nb_orig, ' -> ', domain%nb, ' blocks'
      do b = 1, domain%nb
        write(*,'(A9,I3,A3,3I5)') ' - Block ', b, ' = ', domain%blk(b)%dim
      end do
    end if
   
    ! Rebuild MPI block partitioning with the new block count
    call repartition_blocks(domain)

  end subroutine split_domain_for_mpi


  !> Split block b_idx along direction split_dir at the midpoint.
  !> The first half stays at index b_idx, the second half is appended at nb+1.
  !> All BC entries are remapped accordingly.
  subroutine split_block(domain, b_idx, split_dir)
    use MOSE_Types_m
    use MOSE_Global_m

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: b_idx, split_dir
    ! Local
    integer :: ni, nj, nk, sp, b_new
    integer :: ni1, nj1, nk1, ni2, nj2, nk2
    integer :: i, j, k, d, nbc_new, nbc_old, n_conn
    integer :: face_lo, face_hi
    type(MOSE_block_type), allocatable :: blk_tmp(:)
    type(MOSE_bc_type), allocatable :: bc_tmp(:)
    integer, allocatable :: n_bf_tmp(:,:)

    ni = domain%blk(b_idx)%dim(1)
    nj = domain%blk(b_idx)%dim(2)
    nk = domain%blk(b_idx)%dim(3)

    ! Split position: midpoint
    select case(split_dir)
      case(1); sp = ni / 2; face_lo = 1; face_hi = 2
      case(2); sp = nj / 2; face_lo = 3; face_hi = 4
      case(3); sp = nk / 2; face_lo = 5; face_hi = 6
    end select

    ! Dimensions of the two halves
    ni1 = ni; nj1 = nj; nk1 = nk
    ni2 = ni; nj2 = nj; nk2 = nk
    select case(split_dir)
      case(1); ni1 = sp;      ni2 = ni - sp
      case(2); nj1 = sp;      nj2 = nj - sp
      case(3); nk1 = sp;      nk2 = nk - sp
    end select
    
    b_new = domain%nb + 1

    ! --- 1. Grow the block array by 1 ---
    allocate(blk_tmp(b_new))
    do i = 1, domain%nb
      blk_tmp(i) = domain%blk(i)
    end do
    call move_alloc(blk_tmp, domain%blk)
    domain%nb = b_new
 
    ! --- 2. Allocate the new block (second half) ---
    domain%blk(b_new)%dim = [ni2, nj2, nk2]
    call allocate_split_block(domain%blk(b_new), [ni2, nj2, nk2])

    ! --- 3. Copy data from original block to second half ---
    call copy_block_data(domain%blk(b_idx), domain%blk(b_new), split_dir, sp, &
                         ni, nj, nk, ni2, nj2, nk2)

    ! --- 4. Shrink the original block (first half) ---
    call shrink_block_inplace(domain%blk(b_idx), split_dir, sp, ni, nj, nk, ni1, nj1, nk1)

    ! --- 5. Remap existing BCs ---
    call remap_bcs(domain, b_idx, b_new, split_dir, sp, ni, nj, nk)

    ! --- 6. Generate new type-1 connections at the split interface ---
    ! Count cells on the split face
    select case(split_dir)
      case(1); n_conn = nj * nk
      case(2); n_conn = ni1 * nk  ! ni1 = ni for J/K splits
      case(3); n_conn = ni1 * nj1
    end select

    ! Grow BC array by 2 * n_conn (one pair per face cell)
    nbc_old = domain%nbound
    nbc_new = nbc_old + 2 * n_conn
    allocate(bc_tmp(nbc_new))
    do i = 1, nbc_old
      bc_tmp(i) = domain%bc(i)
    end do
    call move_alloc(bc_tmp, domain%bc)
    domain%nbound = nbc_new

    ! Fill new connection entries
    call generate_split_connections(domain, b_idx, b_new, split_dir, sp, &
                                    ni1, nj1, nk1, ni2, nj2, nk2, &
                                    nbc_old, n_conn)

    ! --- 7. Rebuild n_bf and sort bc array ---
    call rebuild_n_bf(domain)
    call sort_bc_array(domain)

  end subroutine split_block


  !> Allocate arrays for a new block (same pattern as Allocate_Block but simpler).
  subroutine allocate_split_block(blk, nijk)
    use MOSE_Types_m
    use MOSE_Global_m
    use MOSE_IO_RANS, only: obj_rans
    use MOSE_Lib_IRS, only: obj_irs

    implicit none
    type(MOSE_block_type), intent(inout) :: blk
    integer, intent(in) :: nijk(3)
    ! Local
    integer :: n1, n2, n3

    n1 = nijk(1); n2 = nijk(2); n3 = nijk(3)
    blk%dim = nijk

    allocate(blk%node(0:n1, 0:n2, 0:n3))
    allocate(blk%M(1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc))
    allocate(blk%m_f(3, 3, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc))
    allocate(blk%dl(1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc))
    allocate(blk%dl_f(3, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc))
    allocate(blk%vol(1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc))
    allocate(blk%dir(1)%f(0:n1, n2, n3))
    allocate(blk%dir(2)%f(n1, 0:n2, n3))
    allocate(blk%dir(3)%f(n1, n2, 0:n3))
    allocate(blk%fm1(4, 0:n1, n2, n3))
    allocate(blk%fm2(4, n1, 0:n2, n3))
    allocate(blk%fm3(4, n1, n2, 0:n3))
    allocate(blk%yn(1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)); blk%yn = 0.0d0   ! always, zero ghosts included, as MOSE_open (see Mod_Allocate_Data)

    allocate(blk%P(nprim, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc))
    allocate(blk%PO, blk%R, mold=blk%P)
    blk%P = 0.0_R8; blk%PO = 0.0_R8; blk%R = 0.0_R8

    allocate(blk%dtlocal(1:n1, 1:n2, 1:n3))
    allocate(blk%beta(1:n1, 1:n2, 1:n3))
    blk%dtlocal = 0.0_R8; blk%beta = 0.0_R8

    if (obj_irs%IRS) then
      allocate(blk%RS1, blk%RS2, mold=blk%R)
    end if

    blk%no_chem = .false.

  end subroutine allocate_split_block


  !> Copy data from the original block's second half into the new block.
  subroutine copy_block_data(blk_orig, blk_new, sd, sp, ni, nj, nk, ni2, nj2, nk2)
    use MOSE_Types_m
    use MOSE_Global_m, only: nprim

    implicit none
    type(MOSE_block_type), intent(in)    :: blk_orig
    type(MOSE_block_type), intent(inout) :: blk_new
    integer, intent(in) :: sd, sp, ni, nj, nk, ni2, nj2, nk2
    ! Local
    integer :: d

    select case(sd)
      case(1) ! Split along I: second half is i = sp+1..ni
        ! Copy primitive variables (interior only)
        blk_new%P(:, 1:ni2, 1:nj, 1:nk) = blk_orig%P(:, sp+1:ni, 1:nj, 1:nk)
        ! Copy nodes: original node(sp:ni, 0:nj, 0:nk) → new node(0:ni2, 0:nj, 0:nk)
        do d = 1, 3
          blk_new%node(0:ni2, 0:nj, 0:nk)%c(d) = blk_orig%node(sp:ni, 0:nj, 0:nk)%c(d)
        end do

      case(2) ! Split along J: second half is j = sp+1..nj
        blk_new%P(:, 1:ni, 1:nj2, 1:nk) = blk_orig%P(:, 1:ni, sp+1:nj, 1:nk)
        do d = 1, 3
          blk_new%node(0:ni, 0:nj2, 0:nk)%c(d) = blk_orig%node(0:ni, sp:nj, 0:nk)%c(d)
        end do

      case(3) ! Split along K: second half is k = sp+1..nk
        blk_new%P(:, 1:ni, 1:nj, 1:nk2) = blk_orig%P(:, 1:ni, 1:nj, sp+1:nk)
        do d = 1, 3
          blk_new%node(0:ni, 0:nj, 0:nk2)%c(d) = blk_orig%node(0:ni, 0:nj, sp:nk)%c(d)
        end do
    end select

    blk_new%no_chem = blk_orig%no_chem

  end subroutine copy_block_data


  !> Shrink the original block in-place to keep only the first half.
  !> We reallocate the arrays with the smaller dimensions.
  subroutine shrink_block_inplace(blk, sd, sp, ni, nj, nk, ni1, nj1, nk1)
    use MOSE_Types_m
    use MOSE_Global_m
    use MOSE_IO_RANS, only: obj_rans
    use MOSE_Lib_IRS, only: obj_irs

    implicit none
    type(MOSE_block_type), intent(inout) :: blk
    integer, intent(in) :: sd, sp, ni, nj, nk, ni1, nj1, nk1
    ! Local
    type(MOSE_block_type) :: tmp
    integer :: d

    ! Allocate a temporary block with new dimensions
    call allocate_split_block(tmp, [ni1, nj1, nk1])
     
    ! Copy data from first half
    select case(sd)
      case(1)
        tmp%P(:, 1:ni1, 1:nj, 1:nk) = blk%P(:, 1:sp, 1:nj, 1:nk)
        do d = 1, 3
          tmp%node(0:ni1, 0:nj, 0:nk)%c(d) = blk%node(0:sp, 0:nj, 0:nk)%c(d)
        end do
      case(2)
        tmp%P(:, 1:ni, 1:nj1, 1:nk) = blk%P(:, 1:ni, 1:sp, 1:nk)
        do d = 1, 3
          tmp%node(0:ni, 0:nj1, 0:nk)%c(d) = blk%node(0:ni, 0:sp, 0:nk)%c(d)
        end do
      case(3)
        tmp%P(:, 1:ni, 1:nj, 1:nk1) = blk%P(:, 1:ni, 1:nj, 1:sp)
        do d = 1, 3
          tmp%node(0:ni, 0:nj, 0:nk1)%c(d) = blk%node(0:ni, 0:nj, 0:sp)%c(d)
        end do
    end select

    tmp%no_chem = blk%no_chem

    ! Replace block with the shrunk version
    ! Deallocate old arrays and move new ones
    !blk = tmp
   
    call move_block(tmp,blk)

  end subroutine shrink_block_inplace


  !> Deallocate all allocatable arrays in a block.
  subroutine deallocate_block(blk)
    use MOSE_Types_m

    implicit none
    type(MOSE_block_type), intent(inout) :: blk

    if (allocated(blk%node)) deallocate(blk%node)
    if (allocated(blk%M))    deallocate(blk%M)
    if (allocated(blk%dl))   deallocate(blk%dl)
    if (allocated(blk%vol))  deallocate(blk%vol)
    if (allocated(blk%dir(1)%f)) deallocate(blk%dir(1)%f)
    if (allocated(blk%dir(2)%f)) deallocate(blk%dir(2)%f)
    if (allocated(blk%dir(3)%f)) deallocate(blk%dir(3)%f)
    if (allocated(blk%yn))   deallocate(blk%yn)
    if (allocated(blk%P))    deallocate(blk%P)
    if (allocated(blk%PO))   deallocate(blk%PO)
    if (allocated(blk%R))    deallocate(blk%R)
    if (allocated(blk%TE))   deallocate(blk%TE)
    if (allocated(blk%RC))   deallocate(blk%RC)
    if (allocated(blk%RS1))  deallocate(blk%RS1)
    if (allocated(blk%RS2))  deallocate(blk%RS2)
    if (allocated(blk%dtlocal)) deallocate(blk%dtlocal)
    if (allocated(blk%beta))    deallocate(blk%beta)
    if (allocated(blk%vel_gradient)) deallocate(blk%vel_gradient)
    if (allocated(blk%rc_term1))     deallocate(blk%rc_term1)
    if (allocated(blk%rc_term2))     deallocate(blk%rc_term2)

  end subroutine deallocate_block


  subroutine move_block(src, dest)
    use MOSE_Types_m
    implicit none

    type(MOSE_block_type), intent(inout) :: src
    type(MOSE_block_type), intent(inout) :: dest

    call deallocate_block(dest)
    
     dest%dim = src%dim
  
      if (allocated(src%node))     call move_alloc(src%node,    dest%node)
      if (allocated(src%M))        call move_alloc(src%M,    dest%M)
      if (allocated(src%dl))       call move_alloc(src%dl,    dest%dl)
      if (allocated(src%vol))      call move_alloc(src%vol,    dest%vol)
      if (allocated(src%dir(1)%f)) call move_alloc(src%dir(1)%f,    dest%dir(1)%f)
      if (allocated(src%dir(2)%f)) call move_alloc(src%dir(2)%f,    dest%dir(2)%f)
      if (allocated(src%dir(3)%f)) call move_alloc(src%dir(3)%f,    dest%dir(3)%f)
      if (allocated(src%yn))       call move_alloc(src%yn,    dest%yn)
      if (allocated(src%P))        call move_alloc(src%P,    dest%P)
      if (allocated(src%PO))       call move_alloc(src%PO,    dest%PO)
      if (allocated(src%R))        call move_alloc(src%R,    dest%R)
      if (allocated(src%TE))       call move_alloc(src%TE,    dest%TE)
      if (allocated(src%RC))       call move_alloc(src%RC,    dest%RC)
      if (allocated(src%RS1))      call move_alloc(src%RS1,    dest%RS1)
      if (allocated(src%RS2))      call move_alloc(src%RS2,    dest%RS2)
      if (allocated(src%dtlocal))  call move_alloc(src%dtlocal,    dest%dtlocal)
      if (allocated(src%beta))     call move_alloc(src%beta,    dest%beta)
      if (allocated(src%vel_gradient)) call move_alloc(src%vel_gradient,    dest%vel_gradient)
      if (allocated(src%rc_term1))     call move_alloc(src%rc_term1,    dest%rc_term1)
      if (allocated(src%rc_term2))     call move_alloc(src%rc_term2,    dest%rc_term2)

  end subroutine move_block

  !> Remap existing BCs after splitting block b_idx.
  !> Block b_idx now has the first half; b_new has the second half.
  !> nbc_orig is the number of BC entries before this split (before new connections).
  subroutine remap_bcs(domain, b_idx, b_new, sd, sp, ni, nj, nk)
    use MOSE_Types_m
    use MOSE_Global_m, only: nprim

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: b_idx, b_new, sd, sp, ni, nj, nk
    ! Local
    integer :: i

    ! At this point domain%nbound still has the OLD count (before grow).
    ! We iterate over all existing entries.
    do i = 1, domain%nbound

      ! --- Remap destination block ---
      if (domain%bc(i)%b == b_idx) then
        call remap_bc_entry_dest(domain%bc(i), b_idx, b_new, sd, sp)
      end if

      ! --- Remap source block (for type-1 and similar connections) ---
      select case(domain%bc(i)%type)
        case(1, 1000, 1001, 999)
          if (domain%bc(i)%bs == b_idx) then
            call remap_bc_entry_source(domain%bc(i), b_idx, b_new, sd, sp)
          end if
      end select
   
    end do

  end subroutine remap_bcs


  !> Remap the destination (b, i, j, k, f) of a BC entry.
  subroutine remap_bc_entry_dest(bc, b_idx, b_new, sd, sp)
    use MOSE_Types_m

    implicit none
    type(MOSE_bc_type), intent(inout) :: bc
    integer, intent(in) :: b_idx, b_new, sd, sp
    ! Local
    integer :: coord

    ! Get the coordinate along the split direction
    select case(sd)
      case(1); coord = bc%i
      case(2); coord = bc%j
      case(3); coord = bc%k
    end select

    ! Determine which half this BC belongs to
    ! Face on the "high" side of original block → belongs to second half
    ! Interior cells with coord > sp → belongs to second half
    if (is_in_second_half(bc%f, coord, sd, sp)) then
      bc%b = b_new
      ! Shift the coordinate
      select case(sd)
        case(1); bc%i = bc%i - sp
        case(2); bc%j = bc%j - sp
        case(3); bc%k = bc%k - sp
      end select
    end if

  end subroutine remap_bc_entry_dest


  !> Remap the source (bs, is, js, ks, fs) of a type-1 BC entry.
  subroutine remap_bc_entry_source(bc, b_idx, b_new, sd, sp)
    use MOSE_Types_m

    implicit none
    type(MOSE_bc_type), intent(inout) :: bc
    integer, intent(in) :: b_idx, b_new, sd, sp
    ! Local
    integer :: coord

    select case(sd)
      case(1); coord = bc%is
      case(2); coord = bc%js
      case(3); coord = bc%ks
    end select

    if (is_in_second_half(bc%fs, coord, sd, sp)) then
      bc%bs = b_new
      select case(sd)
        case(1); bc%is = bc%is - sp
        case(2); bc%js = bc%js - sp
        case(3); bc%ks = bc%ks - sp
      end select
    end if

  end subroutine remap_bc_entry_source


  !> Determine if a BC cell at (face, coord_along_split) belongs to the second half.
  pure logical function is_in_second_half(face, coord, sd, sp)
    implicit none
    integer, intent(in) :: face, coord, sd, sp
    ! Local
    integer :: face_hi

    ! The "high" face for each split direction
    select case(sd)
      case(1); face_hi = 2  ! I-plus
      case(2); face_hi = 4  ! J-plus
      case(3); face_hi = 6  ! K-plus
      case default; face_hi = 0
    end select

    ! A cell is in the second half if:
    ! 1. It's on the "high" face of the split direction (e.g., face 2 for I-split), OR
    ! 2. Its coordinate along the split direction is > sp
    if (face == face_hi) then
      is_in_second_half = .true.
    else
      is_in_second_half = (coord > sp)
    end if

  end function is_in_second_half


  !> Generate new type-1 connection BCs at the split interface.
  subroutine generate_split_connections(domain, b1, b2, sd, sp, &
                                         ni1, nj1, nk1, ni2, nj2, nk2, &
                                         nbc_start, n_conn)
    use MOSE_Types_m
    use MOSE_Global_m, only: nprim

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: b1, b2, sd, sp
    integer, intent(in) :: ni1, nj1, nk1, ni2, nj2, nk2
    integer, intent(in) :: nbc_start, n_conn
    ! Local
    integer :: j, k, i, idx, face1, face2
    integer :: d11, d12, d21, d22
    integer :: loop1_max, loop2_max
    integer :: i_b1, j_b1, k_b1, i_b2, j_b2, k_b2

    ! Determine face pair and transverse stencil offsets
    select case(sd)
      case(1)
        face1 = 2; face2 = 1  ! I-plus of b1 ↔ I-minus of b2
        d11 = 1; d12 = 0; d21 = 0; d22 = 1
        loop1_max = nj1; loop2_max = nk1
      case(2)
        face1 = 4; face2 = 3  ! J-plus of b1 ↔ J-minus of b2
        d11 = 1; d12 = 0; d21 = 0; d22 = 1
        loop1_max = ni1; loop2_max = nk1
      case(3)
        face1 = 6; face2 = 5  ! K-plus of b1 ↔ K-minus of b2
        d11 = 1; d12 = 0; d21 = 0; d22 = 1
        loop1_max = ni1; loop2_max = nj1
    end select

    idx = nbc_start

    do k = 1, loop2_max
      do j = 1, loop1_max

        ! --- Entry 1: b1 face_hi → b2 face_lo ---
        idx = idx + 1
        select case(sd)
          case(1)
            i_b1 = ni1; j_b1 = j; k_b1 = k
            i_b2 = 1;   j_b2 = j; k_b2 = k
          case(2)
            i_b1 = j; j_b1 = nj1; k_b1 = k
            i_b2 = j; j_b2 = 1;   k_b2 = k
          case(3)
            i_b1 = j; j_b1 = k; k_b1 = nk1
            i_b2 = j; j_b2 = k; k_b2 = 1
        end select

        domain%bc(idx)%b  = b1
        domain%bc(idx)%i  = i_b1
        domain%bc(idx)%j  = j_b1
        domain%bc(idx)%k  = k_b1
        domain%bc(idx)%f  = face1
        domain%bc(idx)%type = 1
        domain%bc(idx)%bs = b2
        domain%bc(idx)%is = i_b2
        domain%bc(idx)%js = j_b2
        domain%bc(idx)%ks = k_b2
        domain%bc(idx)%fs = face2
        domain%bc(idx)%d11 = d11
        domain%bc(idx)%d12 = d12
        domain%bc(idx)%d21 = d21
        domain%bc(idx)%d22 = d22
        allocate(domain%bc(idx)%Pg(nprim, 6))

        ! --- Entry 2: b2 face_lo → b1 face_hi (symmetric) ---
        idx = idx + 1
        domain%bc(idx)%b  = b2
        domain%bc(idx)%i  = i_b2
        domain%bc(idx)%j  = j_b2
        domain%bc(idx)%k  = k_b2
        domain%bc(idx)%f  = face2
        domain%bc(idx)%type = 1
        domain%bc(idx)%bs = b1
        domain%bc(idx)%is = i_b1
        domain%bc(idx)%js = j_b1
        domain%bc(idx)%ks = k_b1
        domain%bc(idx)%fs = face1
        domain%bc(idx)%d11 = d11
        domain%bc(idx)%d12 = d12
        domain%bc(idx)%d21 = d21
        domain%bc(idx)%d22 = d22
        allocate(domain%bc(idx)%Pg(nprim, 6))

      end do
    end do

  end subroutine generate_split_connections


  !> Rebuild n_bf(block, face) from scratch by scanning all BC entries.
  subroutine rebuild_n_bf(domain)
    use MOSE_Types_m

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: i
    integer, allocatable :: nconnect(:), nwall(:), ngsi(:), nSRM(:)
    integer, allocatable :: nio(:), nsym(:), nper(:), nchimera(:)
    integer, allocatable :: nnozzle(:), ncoupled(:), next(:)
    integer, allocatable :: nmanifold(:), nKaffs(:)
   
    if (allocated(domain%n_bf)) deallocate(domain%n_bf)
    allocate(domain%n_bf(domain%nb, 6))
    allocate(nconnect(domain%nb), nwall(domain%nb), ngsi(domain%nb))
    allocate(nSRM(domain%nb), nio(domain%nb), nsym(domain%nb))
    allocate(nper(domain%nb), nchimera(domain%nb), nnozzle(domain%nb))
    allocate(ncoupled(domain%nb), next(domain%nb))
    allocate(nmanifold(domain%nb), nKaffs(domain%nb))
    domain%n_bf = 0
    
    nconnect = 0
    nwall = 0
    ngsi = 0
    nSRM = 0
    nio = 0
    nsym = 0
    nper = 0
    nchimera = 0
    nnozzle = 0
    ncoupled = 0
    nmanifold = 0
    next = 0
    nKaffs = 0

    do i = 1, domain%nbound
      domain%n_bf(domain%bc(i)%b, domain%bc(i)%f) = &
        domain%n_bf(domain%bc(i)%b, domain%bc(i)%f) + 1

        select case( domain%bc(i)%type )
        case(1)
          nconnect(domain%bc(i)%b) = nconnect(domain%bc(i)%b) + 1
        case(2)
          nper(domain%bc(i)%b) = nper(domain%bc(i)%b) + 1
        case(3)
          nsym(domain%bc(i)%b) = nsym(domain%bc(i)%b) + 1
        case(4)
          nio(domain%bc(i)%b) = nio(domain%bc(i)%b) + 1
        case(5)
          nwall(domain%bc(i)%b) = nwall(domain%bc(i)%b) + 1
        case(6)
           nwall(domain%bc(i)%b) = nwall(domain%bc(i)%b) + 1
        case(8)
           nwall(domain%bc(i)%b) = nwall(domain%bc(i)%b) + 1
        case(9,10,12,13)
           nwall(domain%bc(i)%b) = nwall(domain%bc(i)%b) + 1
        case(11)
          next(domain%bc(i)%b) = next(domain%bc(i)%b) + 1
        case(14)
           nwall(domain%bc(i)%b) = nwall(domain%bc(i)%b) + 1
        case(22)
          nnozzle(domain%bc(i)%b) = nnozzle(domain%bc(i)%b) + 1
        case(666)
          nchimera(domain%bc(i)%b) = nchimera(domain%bc(i)%b) + 1
        case(1000)
          ncoupled(domain%bc(i)%b) = ncoupled(domain%bc(i)%b) + 1
        case(1001)
          nmanifold(domain%bc(i)%b) = nmanifold(domain%bc(i)%b) + 1
        case(999)
          nKaffs(domain%bc(i)%b)  = nKaffs(domain%bc(i)%b)  + 1

      end select
    end do
    if (mpi_is_root) then
    write(*,*)
    do i = 1,domain%nb
      write(*,'(A,I8)') ' NEW Boundary Conditions: Block', i
      write(*,'(A,I8)') '  - Connection      --> ', nconnect(i)
      write(*,'(A,I8)') '  - Viscous wall    --> ', nwall(i)
      if (ngsi(i) > 0) then
        write(*,'(A,I8)') '  - GSI (inc. wall) --> ', ngsi(i)
      endif
      write(*,'(A,I8)') '  - Inflow/outflow  --> ', nio(i)
      write(*,'(A,I8)') '  - Symmetry        --> ', nsym(i)
      write(*,'(A,I8)') '  - Periodicity     --> ', nper(i)
      write(*,'(A,I8)') '  - Extrapolation   --> ', next(i)
      if (nchimera(i) > 0) then
        write(*,'(A,I8)') '  - Chimera         --> ', nchimera(i)
      endif
      if (ncoupled(i) > 0) then
        write(*,'(A,I8)') '  - Coupled wall    --> ', ncoupled(i)
      endif
      if (nSRM(i) > 0) then
          write(*,'(A,I8)') '  - SRM wall        --> ', nSRM(i)
      endif
      if (nkaffs(i) > 0) then
        write(*,'(A,I8)') '  - Coupled 1D-MultiD    --> ', nkaffs(i)
      endif
    enddo
   endif
    deallocate(nconnect,nwall,ngsi,nSRM,nio,nsym,nper,nchimera)
    deallocate(nnozzle,ncoupled,next,nmanifold,nKaffs)
  end subroutine rebuild_n_bf


  !> Sort bc array by (block, face) to restore the ordering invariant
  !> assumed by BC_Fluxes, IO_Wall, and other n_bf-based scan loops.
  !> Uses a counting sort based on the already-computed n_bf array.
  subroutine sort_bc_array(domain)
    use MOSE_Types_m

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: i, b, f, pos
    integer, allocatable :: offset(:,:)
    type(MOSE_bc_type), allocatable :: bc_sorted(:)

    ! Compute cumulative offsets from n_bf
    allocate(offset(domain%nb, 6))
    pos = 0
    do b = 1, domain%nb
      do f = 1, 6
        offset(b, f) = pos
        pos = pos + domain%n_bf(b, f)
      end do
    end do

    ! Place each entry in its sorted position
    allocate(bc_sorted(domain%nbound))
    do i = 1, domain%nbound
      b = domain%bc(i)%b
      f = domain%bc(i)%f
      offset(b, f) = offset(b, f) + 1
      bc_sorted(offset(b, f)) = domain%bc(i)
    end do

    ! Replace original array efficiently via move_alloc
    deallocate(domain%bc)
    call move_alloc(bc_sorted, domain%bc)

    deallocate(offset)

  end subroutine sort_bc_array


  !> Rerun MPI block partitioning after splitting.
  subroutine repartition_blocks(domain)
    use MOSE_Types_m

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b
    integer, allocatable :: blk_ncells(:)

    allocate(blk_ncells(domain%nb))
    do b = 1, domain%nb
      blk_ncells(b) = domain%blk(b)%dim(1) * domain%blk(b)%dim(2) * domain%blk(b)%dim(3)
    end do
    call partition_blocks(domain%nb, blk_ncells)
    deallocate(blk_ncells)
  
  end subroutine repartition_blocks




  !> Rebuild IOfield block array to match domain after split.
  !> Must be called after split_domain_for_mpi so that IOfield reflects
  !> the new block count and dimensions.
  subroutine rebuild_IOfield_after_split(IOfield, domain)
    use Lib_ORION_data
    use MOSE_Types_m
    use MOSE_Global_m, only: nprim

    implicit none
    type(ORION_data), intent(inout)       :: IOfield
    type(MOSE_domain_type), intent(in)    :: domain
    ! Local
    integer :: b, d, ni, nj, nk

    ! Nothing to do if block count already matches
    if (size(IOfield%block) == domain%nb) return

    ! Deallocate old blocks
    do b = 1, size(IOfield%block)
      if (allocated(IOfield%block(b)%mesh)) deallocate(IOfield%block(b)%mesh)
      if (allocated(IOfield%block(b)%vars)) deallocate(IOfield%block(b)%vars)
    end do
    deallocate(IOfield%block)
    allocate(IOfield%block(domain%nb))

    do b = 1, domain%nb
      ni = domain%blk(b)%dim(1)
      nj = domain%blk(b)%dim(2)
      nk = domain%blk(b)%dim(3)

      IOfield%block(b)%Ni = ni
      IOfield%block(b)%Nj = nj
      IOfield%block(b)%Nk = nk
      write(IOfield%block(b)%name, '(A,I0)') 'Block', b

      ! Allocate and copy mesh from domain nodes
      allocate(IOfield%block(b)%mesh(1:3, 0:ni, 0:nj, 0:nk))
      do d = 1, 3
        IOfield%block(b)%mesh(d, 0:ni, 0:nj, 0:nk) = domain%blk(b)%node(0:ni, 0:nj, 0:nk)%c(d)
      end do

      ! Allocate vars (will be resized by Output_Solution_Setup)
      allocate(IOfield%block(b)%vars(1:nprim, 1:ni, 1:nj, 1:nk))
      IOfield%block(b)%vars = 0.0_R8
    end do

  end subroutine rebuild_IOfield_after_split


end module MOSE_Mod_BlockSplit
