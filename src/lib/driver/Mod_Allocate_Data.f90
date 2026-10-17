module MOSE_Mod_Allocate_Data
  use, intrinsic :: iso_fortran_env, only : iostat_end

  implicit none

contains


  subroutine Setup_Data_Structure ( domain, IOfield, MGL )
    use MOSE_types_m
    use MOSE_Global_m
    use Lib_ORION_data
    use MOSE_IO_Solution, only: obj_io
    use MOSE_IO_Chemistry, only: obj_chem
    use MOSE_Mod_MPI

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    type(ORION_data), intent(inout)       :: IOfield
    integer, intent(in)                   :: MGL
    ! Local
    integer :: b, d, nblocks
    integer, allocatable :: blk_ncells(:)



    if (nsoot > 0) then
      nc = np + 1
      soot = .true.
    else
      soot = .false.
    end if

    if (npass>0) then
      nps = np + nsoot + 1
    endif
    
    if (nrans > 0) then
      RANS = .true.
    else
      RANS = .false.
    endif
    nres = nprim - nsc + 1

    ! Push the prim-vector indices to the device: the `acc routine seq` flux
    ! callees reference these module scalars directly. RANS is also (re)synced.
    !$acc update device(nsc, nu, nv, nw, np, nt, nrans, nprim, gc, RANS)

    ! Domain is the MOSE-alias of IOfield
    nblocks = size ( IOfield%block )
    allocate( domain%blk( 1:nblocks ) )
    domain%nb = nblocks

    ! Define the ijk dimensions of each block
    do b = 1, nblocks
      domain%blk(b)%dim(1) = IOfield%block(b)%Ni
      domain%blk(b)%dim(2) = IOfield%block(b)%Nj
      domain%blk(b)%dim(3) = IOfield%block(b)%Nk
    enddo

    ! MPI block partitioning (greedy load balancing by cell count)
    allocate(blk_ncells(nblocks))
    do b = 1, nblocks
      blk_ncells(b) = domain%blk(b)%dim(1) * domain%blk(b)%dim(2) * domain%blk(b)%dim(3)
    end do
    call partition_blocks(nblocks, blk_ncells)
    deallocate(blk_ncells)

    if (mpi_is_root) then
      write(*,*)
      write(*,'(A)') ' Mesh Properties: '
    end if
    do b = 1, nblocks
      if (mpi_is_root) &
        write(*,'(A9,I2,A3,3I5,A8,I3)') ' - Block ', b, ' = ', domain%blk(b)%dim, '  rank =', block_owner(b)

      ! Allocate domain block (all blocks on all ranks for ghost cell connectivity)
      call Allocate_Block ( domain%blk(b), domain%blk(b)%dim, MGL )
      
      ! Import domain-block nodes from orion-field
      do d = 1, 3
        domain%blk(b)%node(0:IOfield%block(b)%Ni,0:IOfield%block(b)%Nj,0:IOfield%block(b)%Nk)%c(d) &
        = IOfield%block(b)%mesh(d,0:,0:,0:)
      enddo
      ! Import domain-block primitives from orion-field variables
      do d = 1, nprim
        domain%blk(b)%P(d,1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk) &
        = IOfield%block(b)%vars(d,1:,1:,1:)
      enddo
      domain%time = obj_io%IOtime
      ! Default: chemistry on all blocks
      if (allocated(obj_chem%no_chem_list)) then
        if (any(obj_chem%no_chem_list == b)) then
          domain%blk(b)%no_chem = .true.
        end if
      end if

    enddo

  end subroutine Setup_Data_Structure


  !> Deallocate heavy arrays from non-local blocks to save memory.
  !> On all ranks: deallocates computation arrays (PO, R, TE, RC, etc.)
  !> On non-root ranks: also deallocates metrics (node, M, dl, vol, dir, yn)
  !>   since ghost cell metrics on local blocks were already populated during
  !>   Setup_Metrics and persist (static mesh). Root keeps metrics for wall I/O.
  !> P is deallocated on remote blocks that are NOT needed as chimera donors (666)
  !> or manifold sources (1001). Blocks referenced by these BCs keep P+dir for now.
  subroutine deallocate_remote_computation_data(domain)
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block, mpi_is_root

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: b, d, i, c
    logical, allocatable :: needs_remote_P(:)

    ! Build mask of remote blocks whose P (and dir) must be kept:
    !  - chimera (666): donorID(:,1) can reference remote blocks
    !  - manifold (1001): bc%bs can be remote
    allocate(needs_remote_P(domain%nb))
    needs_remote_P = .false.
    do i = 1, domain%nbound
      select case (domain%bc(i)%type)
        case (666)
          if (allocated(domain%bc(i)%donorID)) then
            do c = 1, size(domain%bc(i)%donorID, 1)
              b = domain%bc(i)%donorID(c, 1)
              if (.not. is_local_block(b)) needs_remote_P(b) = .true.
            end do
          end if
        case (1001)
          b = domain%bc(i)%bs
          if (b > 0) then
            if (.not. is_local_block(b)) needs_remote_P(b) = .true.
          end if
      end select
    end do

    do b = 1, domain%nb
      if (is_local_block(b)) cycle

      ! Computation arrays — free on all ranks
      if (allocated(domain%blk(b)%PO))           deallocate(domain%blk(b)%PO)
      if (allocated(domain%blk(b)%R))            deallocate(domain%blk(b)%R)
      if (allocated(domain%blk(b)%TE))           deallocate(domain%blk(b)%TE)
      if (allocated(domain%blk(b)%RC))           deallocate(domain%blk(b)%RC)
      if (allocated(domain%blk(b)%RS1))          deallocate(domain%blk(b)%RS1)
      if (allocated(domain%blk(b)%RS2))          deallocate(domain%blk(b)%RS2)
      if (allocated(domain%blk(b)%dtlocal))      deallocate(domain%blk(b)%dtlocal)
      if (allocated(domain%blk(b)%beta))         deallocate(domain%blk(b)%beta)
      if (allocated(domain%blk(b)%rho_c))        deallocate(domain%blk(b)%rho_c)
      if (allocated(domain%blk(b)%T_c))          deallocate(domain%blk(b)%T_c)
      if (allocated(domain%blk(b)%hr_c))         deallocate(domain%blk(b)%hr_c)
      if (allocated(domain%blk(b)%Rg_c))         deallocate(domain%blk(b)%Rg_c)
      if (allocated(domain%blk(b)%mil_c))        deallocate(domain%blk(b)%mil_c)
      if (allocated(domain%blk(b)%kl_c))         deallocate(domain%blk(b)%kl_c)
      if (allocated(domain%blk(b)%cp_c))         deallocate(domain%blk(b)%cp_c)

      ! P — free unless this remote block is a chimera donor or manifold source
      if (.not. needs_remote_P(b)) then
        if (allocated(domain%blk(b)%P)) deallocate(domain%blk(b)%P)
      end if

      ! Metrics — free on non-root ranks only (root needs them for wall I/O)
      ! Keep dir on blocks needed for manifold (BC_Manifold reads blk(Bs)%dir%f%A)
      if (.not. mpi_is_root) then
        if (allocated(domain%blk(b)%node)) deallocate(domain%blk(b)%node)
        if (allocated(domain%blk(b)%M))    deallocate(domain%blk(b)%M)
        if (allocated(domain%blk(b)%dl))   deallocate(domain%blk(b)%dl)
        if (allocated(domain%blk(b)%vol))  deallocate(domain%blk(b)%vol)
        if (allocated(domain%blk(b)%yn))   deallocate(domain%blk(b)%yn)
        if (.not. needs_remote_P(b)) then
          do d = 1, 3
            if (allocated(domain%blk(b)%dir(d)%f)) deallocate(domain%blk(b)%dir(d)%f)
          end do
        end if
      end if
    end do

    deallocate(needs_remote_P)

  end subroutine deallocate_remote_computation_data


  subroutine Allocate_Block( blk, nijk, MGL )
    use MOSE_Global_m
    use MOSE_Types_m
    use MOSE_IO_RANS, only: obj_rans
    use MOSE_Lib_IRS, only: obj_irs
    
    implicit none
    integer, intent(in)                  :: nijk(3), MGL
    type(MOSE_block_type), intent(inout) :: blk
    ! Local
    integer :: ni, nj, nk


    ni = nijk(1) ; nj = nijk(2) ; nk = nijk(3)

    ! Metrics
    allocate( blk % node ( 0:ni, 0:nj, 0:nk ) )
    allocate( blk % M    ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    allocate( blk % m_f  ( 3, 3, 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    allocate( blk % dl   ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    allocate( blk % dl_f ( 3, 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    allocate( blk % vol  ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    allocate( blk % dir(1) % f (0:ni, nj, nk) )
    allocate( blk % dir(2) % f (ni, 0:nj, nk) )
    allocate( blk % dir(3) % f (ni, nj, 0:nk) )
    allocate( blk % fm1 (4, 0:ni, nj, nk) )
    allocate( blk % fm2 (4, ni, 0:nj, nk) )
    allocate( blk % fm3 (4, ni, nj, 0:nk) )
    ! yn (wall distance) is allocated for every run, as MOSE_open does (Mod_Allocate_Data): the device BC
    ! record body gathers it for symmetry / connection records whatever the model, and an unallocated yn
    ! is a null device array (illegal address in the 2D Euler tests, 2026-10-15). Zeroed whole, ghost cells
    ! included, as MOSE_open since its allocation refactor (9201c44): a viscous run fills it in Setup_Metrics,
    ! which leaves the ghost cells outside the connections at this zero (before, at whatever memory held).
    allocate ( blk % yn ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    blk % yn = 0.0_R8

    ! Prim and Residuals
    allocate( blk % P (nprim, 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    allocate( blk % PO, blk % R, mold = blk % P )

    ! Thermo cache (cell-center, ghost-inclusive — same shape as a single P slab).
    ! mil_c stays ghost-inclusive because SST_Blk/SSGLRR_Blk/SA_Source_Blk declare
    ! it with explicit ghost-inclusive bounds; consumers still only read interior.
    allocate( blk % rho_c ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    allocate( blk % T_c,   blk % Rg_c, blk % mil_c, mold = blk % rho_c )
    allocate( blk % hr_c, mold = blk % rho_c ); blk % hr_c = 0.0_R8
    ! kl_c, cp_c are consumed only by Compute_Diffusive_Flux at interior cell
    ! centres (averaged at the face). Allocated ghost-inclusive to match the
    ! explicit-shape dummy declaration in Fluxes_blk (Mod_Fluxes.f90); the ghost
    ! layer is never read by the diffusive loop (i = 1..n(1)-1 etc.) and stays
    ! uninitialised — present only to keep dummy/actual strides aligned.
    if (VISC) then
      allocate( blk % kl_c ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
      allocate( blk % cp_c, mold = blk % kl_c )
    end if

    ! Dt cell center with no ghost cells
    allocate( blk % dtlocal ( 1:ni, 1:nj, 1:nk) )

    ! Shock-detector flag with no ghost cells
    allocate( blk % beta ( 1:ni, 1:nj, 1:nk) )

    ! Multigrid forcing function
    if ( MGL > 1 ) then
      allocate( blk % TE, mold = blk % P )
      allocate( blk % RC, mold = blk % P )
    end if

    ! Temp storage for residuals in IRS
    if ( obj_irs%IRS ) then
      allocate( blk % RS1, blk % RS2, mold = blk % R )
    end if

  end subroutine Allocate_Block


end module MOSE_Mod_Allocate_Data
