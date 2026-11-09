!> The block pack of the explicit flux phase (one launch per kernel type over a batch of blocks, the practice of AMReX
!> `AMReX_MFParallelForG.H` / AthenaK `MeshBlockPack` / Parthenon `MeshData`): the per-block geometry copied once into flat
!> device pools indexed like the flat field (dl, metric, wall distance, beta by the flat cell index of P_flat; the face
!> metrics fm1/fm2/fm3 by a per-block chunk offset), a table of the blocks (dims, flat offsets, chunk offsets) and, per
!> batch and per direction, a table of the CUDA blocks of 128 threads: each CUDA block knows its mesh block and the first
!> face (or cell) it covers, so a kernel sweeps the whole batch with the per-block thread -> face mapping of the per-block
!> kernels (same face computed by the same arithmetic: bit-exact). A batch holds as many consecutive local blocks as fit
!> the face-buffer budget (MOSE_FLUX_PACK_MB, default 2048 MB): fhat/dfhat per direction and tface for the batch.
!> MOSE_FLUX_PACK=0 keeps the per-block launches. Off by itself with a shock detector (beta varies per stage) or a host
!> source (RotatingFrame/soot: the var-first r bridge). Only outside a multirate stage (the stages keep boxes and lists).
module MOSE_Pack_m
  use iso_fortran_env, only: R8 => real64
  implicit none
  private
  public :: pack_setup, pack_on, pack_nbat, pack_batch_blocks
  public :: fm1p, fm2p, fm3p, dlp, Mp, ynp, betap
  public :: bt_n1, bt_n2, bt_n3, bt_boff, bt_pi, bt_pj, bt_ofm1, bt_ofm2, bt_ofm3, bt_ofh, bt_ocel, bt_blk
  public :: cb_p, cb_t0, cb_n
  public :: fhp1, fhp2, fhp3, dfhp1, dfhp2, dfhp3, tfp, nfp1, nfp2, nfp3, ntfp
  public :: pack_batch_t, pbat

  integer, parameter, public :: PACK_NT = 128   !< threads per CUDA block of every pack kernel (the per-block kernels use 128 too)

  ! the geometry pools (device resident): fm1p(4, sum of the (n1+1)*n2*n3 chunks), fm2p, fm3p likewise; dlp(3, ncell_total),
  ! Mp(3,3, ncell_total), ynp(ncell_total), betap(ncell_total) by the flat cell index (ghosts included, as P_flat)
  real(R8), allocatable, save :: fm1p(:,:), fm2p(:,:), fm3p(:,:), dlp(:,:), Mp(:,:,:), ynp(:), betap(:)
  ! the block table, p = 1..nloc (the local blocks in block order): dims, flat offsets (boff, pi, pj as
  ! field_block_index_params), chunk offsets into fm1p/fm2p/fm3p (faces), per batch the offsets of the block's faces in the
  ! batch buffers (bt_ofh(d,p)) and of its interior cells in the batch cell range (bt_ocel(p)); bt_blk(p) = global block
  integer, allocatable, save :: bt_n1(:), bt_n2(:), bt_n3(:), bt_boff(:), bt_pi(:), bt_pj(:)
  integer, allocatable, save :: bt_ofm1(:), bt_ofm2(:), bt_ofm3(:), bt_ofh(:,:), bt_ocel(:), bt_blk(:)
  ! the CUDA-block tables, one flat list for every (batch, kind) with kind 1..3 = the faces of direction d, 4 = the cells:
  ! cb_p(e) = the local block p of CUDA block e, cb_t0(e) = the first face/cell of the block it covers (0-based); the range
  ! of a (batch, kind) is cb_n(0:1, kind, batch) = first index - 1 and count
  integer, allocatable, save :: cb_p(:), cb_t0(:), cb_n(:,:,:)
  ! the batch buffers: the face fluxes of the batch per direction (convective fhp*, diffusive dfhp*) and the face thermo
  ! tfp (nf, 5 or 5+nsc), sized to the largest batch
  real(R8), allocatable, save :: fhp1(:,:), fhp2(:,:), fhp3(:,:), dfhp1(:,:), dfhp2(:,:), dfhp3(:,:), tfp(:,:)
  integer, save :: nfp1 = 0, nfp2 = 0, nfp3 = 0, ntfp = 0

  type :: pack_batch_t
    integer :: pA = 0, pB = -1        !< the local blocks p of the batch (consecutive)
    integer :: nf(3) = 0              !< faces of the batch per direction
    integer :: ncel = 0               !< interior cells of the batch
    integer :: ncb(4) = 0             !< CUDA blocks per kind (faces i, j, k, cells)
  end type pack_batch_t
  type(pack_batch_t), allocatable, save :: pbat(:)

  integer, save :: pack_mode = -1     !< -1 not decided, 0 off, 1 on
  integer, save :: nbat = 0, nloc = 0

contains

  logical function pack_on ()
    pack_on = ( pack_mode == 1 )
  end function pack_on

  integer function pack_nbat ()
    pack_nbat = nbat
  end function pack_nbat

  !> The local blocks of batch ib: p = pA..pB (bt_blk(p) the global block)
  subroutine pack_batch_blocks ( ib, pA, pB )
    integer, intent(in)  :: ib
    integer, intent(out) :: pA, pB
    pA = pbat(ib)%pA; pB = pbat(ib)%pB
  end subroutine pack_batch_blocks

  !> Decide once and build everything: the pools, the block table, the batches and their CUDA-block tables, the buffers.
  subroutine pack_setup ( domain )
    use MOSE_Types_m,  only: MOSE_domain_type
    use MOSE_Global_m, only: gc, nprim, nsc, VISC, ROTATING_FRAME, soot
    use MOSE_Mod_MPI,  only: is_local_block, mpi_is_root
    use MOSE_Field_m,  only: ncell_total, field_block_index_params
    use MOSE_Lib_Limiters, only: SD_space_order
    use MOSE_Lib_Shock_Detector, only: sd_dev
    use MOSE_Lib_RANS, only: mcdiff_dev
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    character(32) :: env
    character(64) :: why
    integer :: st, b, p, n1, n2, n3, i, j, k, c, boff, pi, pj, d, ib, e, ncb_tot, t, ntf
    integer(8) :: bytes_face, budget, acc
    integer :: nf(3), chunk
    real(R8) :: mb
    if ( pack_mode >= 0 ) return
    pack_mode = 1; why = ''
    call get_environment_variable ( "MOSE_FLUX_PACK", env, status = st )
    if ( st == 0 ) then
      if ( trim(adjustl(env)) == '0' ) then
        pack_mode = 0; why = 'MOSE_FLUX_PACK=0'
      end if
    end if
    if ( pack_mode == 1 .and. ( SD_space_order .or. sd_dev > 0 ) ) then
      pack_mode = 0; why = 'shock detector on'
    end if
    if ( pack_mode == 1 .and. ( ROTATING_FRAME .or. soot ) ) then
      pack_mode = 0; why = 'host source (rotating frame / soot)'
    end if
    if ( pack_mode == 0 ) then
      if ( mpi_is_root ) write(*,'(A)') ' [PACK] flux phase on per-block launches ('//trim(why)//')'
      return
    end if

    ! ---- the block table and the pool sizes
    nloc = 0
    do b = 1, domain % nb
      if ( is_local_block(b) ) nloc = nloc + 1
    end do
    allocate ( bt_n1(nloc), bt_n2(nloc), bt_n3(nloc), bt_boff(nloc), bt_pi(nloc), bt_pj(nloc), &
               bt_ofm1(nloc), bt_ofm2(nloc), bt_ofm3(nloc), bt_ofh(3, nloc), bt_ocel(nloc), bt_blk(nloc) )
    bt_ofh = 0; bt_ocel = 0
    p = 0; bt_ofm1 = 0; bt_ofm2 = 0; bt_ofm3 = 0
    do b = 1, domain % nb
      if ( .not. is_local_block(b) ) cycle
      p = p + 1
      n1 = domain%blk(b)%dim(1); n2 = domain%blk(b)%dim(2); n3 = domain%blk(b)%dim(3)
      call field_block_index_params ( b, boff, pi, pj )
      bt_blk(p) = b; bt_n1(p) = n1; bt_n2(p) = n2; bt_n3(p) = n3; bt_boff(p) = boff; bt_pi(p) = pi; bt_pj(p) = pj
      if ( p < nloc ) then
        bt_ofm1(p+1) = bt_ofm1(p) + (n1+1)*n2*n3
        bt_ofm2(p+1) = bt_ofm2(p) + n1*(n2+1)*n3
        bt_ofm3(p+1) = bt_ofm3(p) + n1*n2*(n3+1)
      end if
    end do
    p = nloc
    allocate ( fm1p(4, bt_ofm1(p) + (bt_n1(p)+1)*bt_n2(p)*bt_n3(p)) )
    allocate ( fm2p(4, bt_ofm2(p) + bt_n1(p)*(bt_n2(p)+1)*bt_n3(p)) )
    allocate ( fm3p(4, bt_ofm3(p) + bt_n1(p)*bt_n2(p)*(bt_n3(p)+1)) )
    allocate ( dlp(3, ncell_total), Mp(3, 3, ncell_total), ynp(ncell_total), betap(ncell_total) )
    dlp = 0.0_R8; Mp = 0.0_R8; ynp = 0.0_R8; betap = 1.0_R8

    ! ---- fill the pools from the host block arrays (geometry: written once at setup, read-only thereafter)
    do p = 1, nloc
      b = bt_blk(p); n1 = bt_n1(p); n2 = bt_n2(p); n3 = bt_n3(p); boff = bt_boff(p); pi = bt_pi(p); pj = bt_pj(p)
      associate ( bk => domain%blk(b) )
      do k = 1, n3; do j = 1, n2; do i = 0, n1
        fm1p(:, bt_ofm1(p) + i + (n1+1)*((j-1) + n2*(k-1)) + 1) = bk%fm1(:, i, j, k)
      end do; end do; end do
      do k = 1, n3; do j = 0, n2; do i = 1, n1
        fm2p(:, bt_ofm2(p) + (i-1) + n1*(j + (n2+1)*(k-1)) + 1) = bk%fm2(:, i, j, k)
      end do; end do; end do
      do k = 0, n3; do j = 1, n2; do i = 1, n1
        fm3p(:, bt_ofm3(p) + (i-1) + n1*((j-1) + n2*k) + 1) = bk%fm3(:, i, j, k)
      end do; end do; end do
      do k = 1-gc, n3+gc; do j = 1-gc, n2+gc; do i = 1-gc, n1+gc
        c = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
        dlp(:, c) = bk%dl_f(:, i, j, k)
        Mp(:, :, c) = bk%m_f(:, :, i, j, k)
        ynp(c) = bk%yn(i, j, k)
      end do; end do; end do
      end associate
      ! beta = 1 everywhere without a shock detector; the per-block beta of the block path is not refreshed on this path,
      ! so its host and device copies are set to 1 once here (the boundary gather and the per-block BC kernels read it)
      domain%blk(b)%beta = 1.0_R8
      !$acc update device(domain%blk(b)%beta)
    end do
    !$acc enter data copyin(fm1p, fm2p, fm3p, dlp, Mp, ynp, betap)
    !$acc enter data copyin(bt_n1, bt_n2, bt_n3, bt_boff, bt_pi, bt_pj, bt_ofm1, bt_ofm2, bt_ofm3, bt_blk)

    ! ---- the batches: consecutive local blocks within the face-buffer budget
    bytes_face = int(nprim, 8) * 8_8 * 3_8                  ! fhat per direction
    if ( VISC ) bytes_face = bytes_face * 2_8               ! + dfhat
    ntf = 5; if ( mcdiff_dev ) ntf = 5 + nsc
    bytes_face = bytes_face + int(ntf, 8) * 8_8             ! tface (one direction at a time, counted once)
    mb = 2048.0_R8
    call get_environment_variable ( "MOSE_FLUX_PACK_MB", env, status = st )
    if ( st == 0 ) then
      read(env, *, iostat=st) mb
      if ( st /= 0 .or. mb <= 0.0_R8 ) mb = 2048.0_R8
    end if
    budget = int(mb * 1048576.0_R8, 8)
    allocate ( pbat(nloc) ); nbat = 0; acc = 0
    do p = 1, nloc
      n1 = bt_n1(p); n2 = bt_n2(p); n3 = bt_n3(p)
      nf = (/ (n1-1)*n2*n3, n1*(n2-1)*n3, n1*n2*(n3-1) /)
      if ( nbat == 0 .or. acc + int(maxval(nf), 8) * bytes_face > budget ) then
        nbat = nbat + 1; pbat(nbat)%pA = p; pbat(nbat)%pB = p; pbat(nbat)%nf = 0; pbat(nbat)%ncel = 0; acc = 0
      end if
      pbat(nbat)%pB = p
      do d = 1, 3
        bt_ofh(d, p) = pbat(nbat)%nf(d); pbat(nbat)%nf(d) = pbat(nbat)%nf(d) + nf(d)
      end do
      bt_ocel(p) = pbat(nbat)%ncel; pbat(nbat)%ncel = pbat(nbat)%ncel + n1*n2*n3
      acc = acc + int(maxval(nf), 8) * bytes_face
    end do
    !$acc enter data copyin(bt_ofh, bt_ocel)
    nfp1 = 1; nfp2 = 1; nfp3 = 1; ntfp = 1
    do ib = 1, nbat
      nfp1 = max( nfp1, pbat(ib)%nf(1) ); nfp2 = max( nfp2, pbat(ib)%nf(2) ); nfp3 = max( nfp3, pbat(ib)%nf(3) )
    end do
    ntfp = max( nfp1, nfp2, nfp3 )
    allocate ( fhp1(nfp1, nprim), fhp2(nfp2, nprim), fhp3(nfp3, nprim) )
    fhp1 = 0.0_R8; fhp2 = 0.0_R8; fhp3 = 0.0_R8
    !$acc enter data copyin(fhp1, fhp2, fhp3)
    if ( VISC ) then
      allocate ( dfhp1(nfp1, nprim), dfhp2(nfp2, nprim), dfhp3(nfp3, nprim), tfp(ntfp, ntf) )
      dfhp1 = 0.0_R8; dfhp2 = 0.0_R8; dfhp3 = 0.0_R8; tfp = 0.0_R8
      !$acc enter data copyin(dfhp1, dfhp2, dfhp3, tfp)
    else
      allocate ( dfhp1(1, nprim), dfhp2(1, nprim), dfhp3(1, nprim), tfp(1, ntf) )
      dfhp1 = 0.0_R8; dfhp2 = 0.0_R8; dfhp3 = 0.0_R8; tfp = 0.0_R8
      !$acc enter data copyin(dfhp1, dfhp2, dfhp3, tfp)
    end if

    ! ---- the CUDA-block tables: per (batch, kind) the CUDA blocks of 128 threads of every block of the batch in block
    ! order, each with its block p and its first face/cell (0-based, within the block's own face/cell range)
    ncb_tot = 0
    allocate ( cb_n(0:1, 4, nbat) ); cb_n = 0
    do ib = 1, nbat
      do p = pbat(ib)%pA, pbat(ib)%pB
        n1 = bt_n1(p); n2 = bt_n2(p); n3 = bt_n3(p)
        nf = (/ (n1-1)*n2*n3, n1*(n2-1)*n3, n1*n2*(n3-1) /)
        do d = 1, 3
          pbat(ib)%ncb(d) = pbat(ib)%ncb(d) + (nf(d) + PACK_NT - 1) / PACK_NT
        end do
        pbat(ib)%ncb(4) = pbat(ib)%ncb(4) + (n1*n2*n3 + PACK_NT - 1) / PACK_NT
      end do
      ncb_tot = ncb_tot + sum( pbat(ib)%ncb )
    end do
    allocate ( cb_p(max(ncb_tot,1)), cb_t0(max(ncb_tot,1)) ); cb_p = 1; cb_t0 = 0
    e = 0
    do ib = 1, nbat
      do d = 1, 4
        cb_n(0, d, ib) = e; cb_n(1, d, ib) = pbat(ib)%ncb(d)
        do p = pbat(ib)%pA, pbat(ib)%pB
          n1 = bt_n1(p); n2 = bt_n2(p); n3 = bt_n3(p)
          select case ( d )
          case (1); chunk = (n1-1)*n2*n3
          case (2); chunk = n1*(n2-1)*n3
          case (3); chunk = n1*n2*(n3-1)
          case (4); chunk = n1*n2*n3
          end select
          do t = 0, chunk - 1, PACK_NT
            e = e + 1; cb_p(e) = p; cb_t0(e) = t
          end do
        end do
      end do
    end do
    !$acc enter data copyin(cb_p, cb_t0, cb_n)

    if ( mpi_is_root ) then
      write(*,'(A,I0,A,I0,A,I0,A,F8.1,A)') ' [PACK] flux phase on block packs: ', nloc, ' blocks in ', nbat, &
            ' batches (budget ', nint(mb), ' MB), face buffers ', &
            real(nfp1 + nfp2 + nfp3, R8) * real(nprim, R8) * 8.0_R8 * merge(2.0_R8, 1.0_R8, VISC) / 1048576.0_R8 + &
            real(ntfp, R8) * real(ntf, R8) * 8.0_R8 / 1048576.0_R8, ' MB'
      write(*,'(A,F8.1,A,I0,A)') ' [PACK] geometry pools ', &
            ( real(size(fm1p) + size(fm2p) + size(fm3p) + size(dlp) + size(Mp) + size(ynp) + size(betap), R8) * 8.0_R8 ) / 1048576.0_R8, &
            ' MB, CUDA blocks of all batches and kinds ', ncb_tot, ' (128 threads each)'
    end if
  end subroutine pack_setup

end module MOSE_Pack_m
