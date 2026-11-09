module MOSE_Mod_BC
  use MOSE_Lib_Multirate_Base, only: mr_bl_type

  implicit none

  ! The record lists of the boundary fluxes of the explicit step (bc_list_setup): one per block face, every record of
  ! every local block, gathered as the lists of a multirate level; bc_list_mode -1 = not decided yet, 0 = the per-block
  ! launches (BC_Fluxes_blk), 1 = the lists (BC_Fluxes_list). bla_nb / bla_nbound: the domain the lists were built for.
  type(mr_bl_type), allocatable, save :: bla(:)
  integer, save :: bc_list_mode = -1, bla_nb = -1, bla_nbound = -1

! PWARG threads either (w, ncw, boff, pi, pj) or (P) through the
! BC dispatch kernel's SIGNATURE and its per-entry BC_*_dev CALL sites identically
! (the BC twins gather from the var-LAST flat SoA w on the GPU device path (_OPENACC)). Must
! mirror the identical macro in Lib_BC.f90 (cpp macros are file-scoped). Flag OFF ->
! PWARG == P, the committed P-path byte-unchanged.
#  define PWARG w, wr, ncw, nrw, boff, pi, pj, wroff, ur, nur
#  define PWPASS P_flat, R_flat, ncw, nrw, boff, pi, pj, wroff, ur_f, nur

! GPU lever A (BC concurrency): launch each block's BC face kernels on its OWN async
! queue (qb = b). Within a block all face-groups share the queue -> in-order -> the
! corner-cell residual accumulation order is preserved (bit-exact + deterministic).
! Across blocks the queues run CONCURRENTLY, filling the GPU -- each per-face launch is
! only ~nbnd/nb/6 cells (~10K at THOR scale), severely under-filling a serial launch.
! OFF -> BCASYNC expands empty AND the per-block `acc data` region is kept => the sync
! path is byte-unchanged. (Async is unconditional since 2026-07-11; the
! clause only ever lands on the resident-w parallel loop.)
#define BCASYNC async(qb)

contains

  !----------------------------------------------------------------------------
  ! Flatten the boundary-face list into plain int arrays the
  ! device kernel can map (the bc(:) derived type has allocatable components —
  ! the AoS device blocker). bc(:) is laid out in (block, face) order, exactly
  ! as the host BC_Fluxes driver walks it via n_bf, so bc_blk_lo/hi mark each
  ! block's contiguous index range. Runs on both builds (cheap, setup-only); the
  ! arrays are read only by the device path. Idempotent across MG levels.
  !----------------------------------------------------------------------------
  !> The boundary types the device handles, checked right after the bc file is read and before the halo schedules
  !> bind their buffers: the ghost fill and the BC flux kernel know none (0), connection (1), rotational periodic (2),
  !> symmetry (3), inflow / outflow (4), viscous walls (5, 6), extrapolation (11), injector (22), outlet (24), chimera
  !> (666, donors on the same rank) and the manifold (1001, one rank); anything else stops the run here.
  subroutine Check_Device_Boundary_Types ( domain )
    use MOSE_Types_m
    use MOSE_Global_m, only: ROTATING_FRAME
    use MOSE_Mod_MPI, only: is_local_block, mpi_size_
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer :: i, nbnd
    nbnd = domain % nbound
    do i = 1, nbnd
      if (.not. is_local_block(domain % bc(i) % b)) cycle
      select case (domain % bc(i) % type)
        case (0, 1, 2, 3, 4, 5, 6, 11, 22, 24)
          ! type 0 = "none"/degenerate 2D k-faces: device no-op flux (no dispatch case) +
          ! zero-gradient ghost (matches host Fill_Ghost_Cell case(0)); k-gradient is 0 on a
          ! 1-cell-k mesh so there is no spanwise flux to compute.
          if (domain % bc(i) % type == 4 .and. ROTATING_FRAME) &
            error stop 'F4c device BC: rotating-frame inflow (type 4) not ported'
        case (666)
          ! chimera: device gather ghost (Ghost_Chimera) + connection flux twin.
          ! Every donor must be local (the host fallback would read stale P on the
          ! resident path): single rank, or a partition that co-locates receptor and
          ! donors (INPUT/block_owner.txt).
          block
            integer :: cc
            do cc = 1, sum(domain % bc(i) % ni)
              if (.not. is_local_block(domain % bc(i) % donorID(cc,1))) then
                write(*,*) 'F4c device BC: chimera receptor block ', domain % bc(i) % b, &
                           ' has donor block ', domain % bc(i) % donorID(cc,1), ' on another rank'
                error stop 'F4c device BC: chimera (666) donors must be local -- co-locate them via INPUT/block_owner.txt'
              end if
            end do
          end block
        case (667)
          ! MOSE_open 410 (mapped Q2D state): device ghost from the file at the stage time (Ghost_Q2D_Layers/_Edges,
          ! Mod_BC_Q2D) + the connection twins; local to the rank of its block
          if (ROTATING_FRAME) &
            error stop 'F4c device BC: 667 (MOSE_open 410) in a rotating frame not ported'
        case (1001)
          ! manifold mean connection: device source-face reduction + inflow twin.
          ! Single-rank only (the reduction reads the local w arena).
          if (mpi_size_ > 1) &
            error stop 'F4c device BC: manifold (1001) is single-rank only on the device path'
          if (.not. is_local_block(domain % bc(i) % bs)) &
            error stop 'F4c device BC: manifold (1001) source block not local'
        case default
          write(*,*) 'BC_Fluxes device kernel: unsupported BC type ', domain % bc(i) % type
          if (domain % bc(i) % type == 24 .or. domain % bc(i) % type == 25) &
            write(*,*) '  (MOSE_open BC ', domain % bc(i) % open_type, ' -- host closure only, no device twin yet)'
          error stop 'F4c device BC: only none(0), conn(1), periodic(2), symm(3), in/out(4), wall(5), wall-T(6), extrap(11), chimera(666), manifold(1001) ported'
      end select
    end do

  end subroutine Check_Device_Boundary_Types


  subroutine Setup_BC_Device ( domain )
    use MOSE_Types_m
    use MOSE_Global_m, only: nsc, nrans, nprim, ROTATING_FRAME
    use MOSE_Mod_MPI, only: is_local_block, mpi_size_, mpi_is_root
    use MOSE_Mod_GhostExchange, only: finalize_pg_tables
    use MOSE_Lib_BC, only: conn_fast_dev

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer  :: b, f, i, js, upper, nb, nbnd, Bm, Im, Jm, Km, Fm, Fi, Fj, Fk
    real(R8) :: tv(3), tmod

    nb   = domain % nb
    nbnd = domain % nbound
    if (nbnd <= 0) return
    if (allocated(domain % bc_bi)) return        ! already built for this level

    ! mr-conn-fast (Lib_BC conn_fast_dev): MOSE_CONN_FAST = 0 (as before) / 1 (face state once) / 2 (lean, default)
    block
      character(8) :: env
      integer :: ln, st
      call get_environment_variable ( 'MOSE_CONN_FAST', env, ln, st )
      conn_fast_dev = 2
      if ( st == 0 .and. ln > 0 ) then
        if ( env(1:1) == '0' .or. env(1:1) == 'n' .or. env(1:1) == 'N' .or. env(1:1) == 'f' .or. env(1:1) == 'F' ) &
          conn_fast_dev = 0
        if ( env(1:1) == '1' ) conn_fast_dev = 1
      end if
      !$acc update device(conn_fast_dev)
      if ( mpi_is_root ) write(*,'(A,I0)') ' [BC] mr-conn-fast (connection records: 0 as before, 1 face state once, 2 lean): ', &
                                           conn_fast_dev
    end block

    allocate ( domain % bc_bi(nbnd), domain % bc_bj(nbnd), domain % bc_bk(nbnd), &
               domain % bc_bf(nbnd), domain % bc_typ(nbnd) )
    allocate ( domain % bc_blk_lo(nb), domain % bc_blk_hi(nb) )

    ! The host bc(:) list is in bc.txt FILE order, which is NOT (block,face)-
    ! contiguous (e.g. dettube interleaves faces 3/4 one cell at a time; multi-
    ! block ATLAS meshes interleave block1-f2 / block2-f1). The device dispatcher
    ! (BC_Fluxes_blk) slices [lo,hi] windows from bc_blk_lo/hi + n_bf and
    ! ASSUMES each window is one block / one face — so an unordered list makes a
    ! single acc launch mix blocks (wrong P/r/geometry) and faces (corner r+=
    ! race). Build a stable block-major,face-major permutation here (counting
    ! sort using the n_bf bucket sizes) so the flat device arrays satisfy the
    ! documented invariant. The host bc(:) list is left untouched (reordering it
    ! would change the CPU corner-residual summation order -> breaks O1).
    allocate ( domain % bc_perm(nbnd) )
    block
      integer :: bf_off(nb,6), bb, ff
      upper = 0
      do b = 1, nb
        do f = 1, 6
          bf_off(b,f) = upper          ! 0-based start of bucket (b,f)
          upper = upper + domain % n_bf(b,f)
        end do
      end do
      do i = 1, nbnd                   ! stable placement in file order
        bb = domain % bc(i) % b ; ff = domain % bc(i) % f
        bf_off(bb,ff) = bf_off(bb,ff) + 1
        domain % bc_perm(bf_off(bb,ff)) = i
      end do
    end block

    do i = 1, nbnd
      js = domain % bc_perm(i)
      domain % bc_bi(i)  = domain % bc(js) % i
      domain % bc_bj(i)  = domain % bc(js) % j
      domain % bc_bk(i)  = domain % bc(js) % k
      domain % bc_bf(i)  = domain % bc(js) % f
      domain % bc_typ(i) = domain % bc(js) % type
    end do

    upper = 0
    do b = 1, nb
      domain % bc_blk_lo(b) = upper + 1
      do f = 1, 6
        upper = upper + domain % n_bf(b,f)
      end do
      domain % bc_blk_hi(b) = upper
    end do

    ! Flat per-face inflow/outflow (type 4) scalars + tangential face vector.
    allocate ( domain % bc_mach(nbnd), domain % bc_t0(nbnd), domain % bc_p0(nbnd), &
               domain % bc_relfac(nbnd), domain % bc_alpha(nbnd), domain % bc_beta(nbnd), &
               domain % bc_pamb(nbnd), domain % bc_ci(nsc,nbnd), domain % bc_rans(max(nrans,1),nbnd), &
               domain % bc_tvec(3,nbnd), domain % bc_pg(nprim,6,nbnd), domain % bc_tw(nbnd) )
    domain % bc_mach = 0d0; domain % bc_t0 = 0d0; domain % bc_p0 = 0d0; domain % bc_relfac = 0d0
    domain % bc_alpha = 0d0; domain % bc_beta = 0d0; domain % bc_pamb = 0d0
    domain % bc_ci = 0d0; domain % bc_rans = 0d0; domain % bc_tvec = 0d0; domain % bc_pg = 0d0
    domain % bc_tw = 0d0

    do i = 1, nbnd
      js = domain % bc_perm(i)                                                   ! sorted slot i <- file-order bc(js)
      if (domain % bc(js) % type == 6) then                                      ! isothermal wall: only Tw
        domain % bc_tw(i) = domain % bc(js) % Tw
        cycle
      end if
      if (domain % bc(js) % type /= 4 .and. domain % bc(js) % type /= 22 .and. &
          domain % bc(js) % type /= 1001 .and. domain % bc(js) % type /= 24) cycle   ! 22 + manifold(1001) routed as 4; 24 = outlet (pamb, rel_fac)
      Bm = domain % bc(js) % b
      if (.not. is_local_block(Bm)) cycle
      if (domain % bc(js) % type == 1001) then
        ! Manifold mean connection: static inflow params mirror the host BC_Fluxes
        ! case(1001) call (mach=-10 mass-flux mode, rel_fac=1, ci=1, rans=1e-6; the
        ! dynamic T0/g/pstat land in bc_t0/bc_p0/bc_pamb each substage from the
        ! device source-face reduction in BC_Manifold_Update).
        domain % bc_mach(i)   = -10.0d0
        domain % bc_relfac(i) = 1.0d0
        domain % bc_alpha(i)  = 0.0d0
        domain % bc_beta(i)   = 0.0d0
        domain % bc_ci(1:nsc,i) = 1.0d0
        if (nrans > 0) domain % bc_rans(1:nrans,i) = 1.0d-6
      else
      domain % bc_mach(i)   = domain % bc(js) % mach
      domain % bc_t0(i)     = domain % bc(js) % t0
      domain % bc_p0(i)     = domain % bc(js) % p0
      domain % bc_relfac(i) = domain % bc(js) % rel_fac
      domain % bc_alpha(i)  = domain % bc(js) % alpha
      domain % bc_beta(i)   = domain % bc(js) % beta
      domain % bc_pamb(i)   = domain % bc(js) % pamb
      ! bc%ci is laid out [ RANS(1:nrans) , species(1:nsc) ] (see host BC_Fluxes); the outlet (24) has no composition
      if (allocated(domain % bc(js) % ci)) then
        domain % bc_ci(1:nsc,i) = domain % bc(js) % ci(nrans+1:nrans+nsc)
        if (nrans > 0) domain % bc_rans(1:nrans,i) = domain % bc(js) % ci(1:nrans)
      end if
      end if
      ! Unit tangential face vector from the node geometry (mirror BC_Inflow_Outflow)
      Im = domain % bc(js) % i; Jm = domain % bc(js) % j; Km = domain % bc(js) % k
      Fm = domain % bc(js) % f
      select case (Fm)
        case (1:2); Fi = Im - mod(Fm,2); Fj = Jm;            Fk = Km
        case (3:4); Fi = Im;             Fj = Jm - mod(Fm,2); Fk = Km
        case (5:6); Fi = Im;             Fj = Jm;             Fk = Km - mod(Fm,2)
      end select
      select case (Fm)
        case (1:2); tv = domain % blk(Bm) % node(Fi,Fj,Fk) % c - domain % blk(Bm) % node(Fi,Fj-1,Fk-1) % c
        case (3:4); tv = domain % blk(Bm) % node(Fi,Fj,Fk) % c - domain % blk(Bm) % node(Fi-1,Fj,Fk-1) % c
        case (5:6); tv = domain % blk(Bm) % node(Fi,Fj,Fk) % c - domain % blk(Bm) % node(Fi-1,Fj-1,Fk) % c
      end select
      tmod = norm2(tv)
      domain % bc_tvec(:,i) = tv / tmod
    end do

    ! ---- F5: device connection (type-1) ghost-fill tables (single-rank) ----
    ! Flatten type-1 connection entries whose donor block is ALSO local into per-entry
    ! donor/receiver index arrays, grouped by ordered (Bm,Bs) pair so Fill_Ghost_Cell_acc
    ! fills the cross-block ghosts + bc_pg on device with ONE kernel per pair (both block
    ! P arrays present-matched). Inter-rank donors are left to the host path.
    block
      integer :: nc, c, p, np_, js2
      integer, allocatable :: slot_of(:)
      integer, allocatable :: tbm(:), tbs(:), tim(:), tjm(:), tkm(:), tfm(:)
      integer, allocatable :: tis(:), tjs(:), tks(:), tfs(:)
      integer, allocatable :: td11(:), td12(:), td21(:), td22(:), tpg(:), tpair(:)
      integer, allocatable :: pbm(:), pbs(:), cnt(:), pos(:)
      allocate ( slot_of(nbnd) )                 ! inverse perm: file index -> sorted slot
      do i = 1, nbnd
        slot_of(domain % bc_perm(i)) = i
      end do
      nc = 0
      do js2 = 1, nbnd
        if (domain % bc(js2) % type /= 1) cycle
        if (.not. is_local_block(domain % bc(js2) % b))  cycle
        if (.not. is_local_block(domain % bc(js2) % bs)) cycle
        nc = nc + 1
      end do
      domain % n_conn = nc
      if (nc > 0) then
        allocate ( tbm(nc), tbs(nc), tim(nc), tjm(nc), tkm(nc), tfm(nc),     &
                   tis(nc), tjs(nc), tks(nc), tfs(nc), td11(nc), td12(nc),   &
                   td21(nc), td22(nc), tpg(nc), tpair(nc), pbm(nc), pbs(nc), &
                   cnt(nc), pos(nc) )
        c = 0
        do js2 = 1, nbnd
          if (domain % bc(js2) % type /= 1) cycle
          if (.not. is_local_block(domain % bc(js2) % b))  cycle
          if (.not. is_local_block(domain % bc(js2) % bs)) cycle
          c = c + 1
          tbm(c)=domain%bc(js2)%b;   tbs(c)=domain%bc(js2)%bs
          tim(c)=domain%bc(js2)%i;   tjm(c)=domain%bc(js2)%j;   tkm(c)=domain%bc(js2)%k;   tfm(c)=domain%bc(js2)%f
          tis(c)=domain%bc(js2)%is;  tjs(c)=domain%bc(js2)%js;  tks(c)=domain%bc(js2)%ks;  tfs(c)=domain%bc(js2)%fs
          td11(c)=domain%bc(js2)%d11; td12(c)=domain%bc(js2)%d12
          td21(c)=domain%bc(js2)%d21; td22(c)=domain%bc(js2)%d22
          tpg(c)=slot_of(js2)
        end do
        np_ = 0; cnt(:) = 0                       ! distinct ordered pairs + counts
        do c = 1, nc
          p = 0
          do i = 1, np_
            if (pbm(i)==tbm(c) .and. pbs(i)==tbs(c)) then; p = i; exit; end if
          end do
          if (p == 0) then
            np_ = np_ + 1; pbm(np_)=tbm(c); pbs(np_)=tbs(c); p = np_
          end if
          tpair(c) = p; cnt(p) = cnt(p) + 1
        end do
        domain % n_conn_pair = np_
        allocate ( domain%conn_im(nc), domain%conn_jm(nc), domain%conn_km(nc), domain%conn_fm(nc), &
                   domain%conn_is(nc), domain%conn_js(nc), domain%conn_ks(nc), domain%conn_fs(nc), &
                   domain%conn_d11(nc), domain%conn_d12(nc), domain%conn_d21(nc), domain%conn_d22(nc), &
                   domain%conn_pg(nc) )
        allocate ( domain%conn_pair_bm(np_), domain%conn_pair_bs(np_), &
                   domain%conn_pair_lo(np_), domain%conn_pair_hi(np_) )
        allocate ( domain%conn_all(nc), domain%conn_all_bs(nc) )
        upper = 0                                  ! prefix offsets -> per-pair [lo,hi]
        do p = 1, np_
          domain%conn_pair_bm(p) = pbm(p); domain%conn_pair_bs(p) = pbs(p)
          domain%conn_pair_lo(p) = upper + 1
          upper = upper + cnt(p)
          domain%conn_pair_hi(p) = upper
          pos(p) = domain%conn_pair_lo(p)
        end do
        do c = 1, nc                               ! stable placement by pair
          p = tpair(c); i = pos(p); pos(p) = pos(p) + 1
          domain%conn_im(i)=tim(c); domain%conn_jm(i)=tjm(c); domain%conn_km(i)=tkm(c); domain%conn_fm(i)=tfm(c)
          domain%conn_is(i)=tis(c); domain%conn_js(i)=tjs(c); domain%conn_ks(i)=tks(c); domain%conn_fs(i)=tfs(c)
          domain%conn_d11(i)=td11(c); domain%conn_d12(i)=td12(c); domain%conn_d21(i)=td21(c); domain%conn_d22(i)=td22(c)
          domain%conn_pg(i)=tpg(c)
          domain%conn_all(i)=i; domain%conn_all_bs(i)=tbs(c)
        end do
      end if
      deallocate ( slot_of )
    end block

    !  GPU Phase-1 follow-on: now that bc_perm exists, build the per-recv-entry bc_pg
    ! slot table and the persistent transverse-Pg requests for the device Pg exchange. No-op
    ! single-rank / non-MPI. The send-side Pg tables + resident buffers were already built in
    ! build_ghost_schedule (they need no bc_perm).
    if (mpi_size_ > 1) call finalize_pg_tables(domain)

  end subroutine Setup_BC_Device


  !----------------------------------------------------------------------------
  ! Chimera (666) + manifold (1001) flat device tables. Separate from
  ! Setup_BC_Device because every cell reference here is a GLOBAL P_flat index
  ! (MOSE_Field_m cidx), so the field arena must exist. Called once from the
  ! Mod_Explicit geometry-residency block; allocates + enter-data's its tables.
  !----------------------------------------------------------------------------
  subroutine Setup_BC_Device_Flat ( domain )
    use MOSE_Types_m
    use MOSE_Global_m, only: nsc, nprim, np
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Parameters_m, only: guide
    use MOSE_Field_m, only: cidx, field_alloc, field_is_allocated

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: i, js, e, c, nd, ndon, Bm, Fm, Bs, Fs, dirm, dirs, modf
    integer :: im, jm, km, is_, js_, ks_, Fi, Fj, Fk, ig, jg, kg, n1, n2, n3
    integer :: i0, i1, j0, j1, k0, k1, maxc, ncf, key(4), cur(4)
    logical :: open_grp

    if (domain % nbound <= 0) return
    if (allocated(domain % chim_g1) .or. allocated(domain % mf_cf)) return
    if (.not. field_is_allocated()) call field_alloc(domain)

    ! ---- chimera flat tables --------------------------------------------------
    domain % n_chim = 0; ndon = 0
    do i = 1, domain % nbound
      js = domain % bc_perm(i)
      if (domain % bc(js) % type /= 666) cycle
      if (.not. is_local_block(domain % bc(js) % b)) cycle
      domain % n_chim = domain % n_chim + 1
      ndon = ndon + sum(domain % bc(js) % ni)
    end do
    if (domain % n_chim > 0) then
      allocate ( domain % chim_g1(domain % n_chim), domain % chim_g2(domain % n_chim), &
                 domain % chim_rc(domain % n_chim), domain % chim_dlo(domain % n_chim), &
                 domain % chim_dmid(domain % n_chim), domain % chim_dhi(domain % n_chim), &
                 domain % chim_slot(domain % n_chim), domain % chim_df(ndon), domain % chim_vf(ndon) )
      allocate ( domain % chim_tn(4, domain % n_chim), domain % chim_cg(4, domain % n_chim) )
      e = 0; nd = 0
      do i = 1, domain % nbound
        js = domain % bc_perm(i)
        if (domain % bc(js) % type /= 666) cycle
        if (.not. is_local_block(domain % bc(js) % b)) cycle
        e = e + 1
        Bm = domain % bc(js) % b; Fm = domain % bc(js) % f
        im = domain % bc(js) % i; jm = domain % bc(js) % j; km = domain % bc(js) % k
        domain % chim_rc(e)   = cidx(Bm, im, jm, km)
        domain % chim_g1(e)   = cidx(Bm, im -   guide(Fm,1), jm -   guide(Fm,2), km -   guide(Fm,3))
        domain % chim_g2(e)   = cidx(Bm, im - 2*guide(Fm,1), jm - 2*guide(Fm,2), km - 2*guide(Fm,3))
        domain % chim_slot(e) = i
        domain % chim_dlo(e)  = nd + 1
        domain % chim_dmid(e) = nd + domain % bc(js) % ni(1)
        do c = 1, sum(domain % bc(js) % ni)
          if (.not. is_local_block(domain % bc(js) % donorID(c,1))) &
            error stop 'Setup_BC_Device_Flat: chimera donor block not local'
          nd = nd + 1
          domain % chim_df(nd) = cidx(domain % bc(js) % donorID(c,1), domain % bc(js) % donorID(c,2), &
                                      domain % bc(js) % donorID(c,3), domain % bc(js) % donorID(c,4))
          domain % chim_vf(nd) = domain % bc(js) % volume_fraction(c)
        end do
        domain % chim_dhi(e) = nd
        ! Transverse viscous stencil + block-edge corner ghosts, exactly as the host
        ! Fill_BC_Ghost_Chimera: bc_pg(:,3:6) <- the boundary cell's four transverse
        ! neighbours; at a block edge the first-layer ghost's transverse ghost-of-ghost
        ! cell is a copy of the ghost itself. Flat indices precomputed here.
        n1 = domain % blk(Bm) % dim(1); n2 = domain % blk(Bm) % dim(2); n3 = domain % blk(Bm) % dim(3)
        ig = im - guide(Fm,1); jg = jm - guide(Fm,2); kg = km - guide(Fm,3)
        domain % chim_cg(:,e) = 0
        select case (Fm)
          case (1:2)
            domain % chim_tn(:,e) = (/ cidx(Bm,im,jm-1,km), cidx(Bm,im,jm+1,km), cidx(Bm,im,jm,km-1), cidx(Bm,im,jm,km+1) /)
            if (jm == 1)  domain % chim_cg(1,e) = cidx(Bm,ig,jm-1,km)
            if (jm == n2) domain % chim_cg(2,e) = cidx(Bm,ig,jm+1,km)
            if (km == 1)  domain % chim_cg(3,e) = cidx(Bm,ig,jm,km-1)
            if (km == n3) domain % chim_cg(4,e) = cidx(Bm,ig,jm,km+1)
          case (3:4)
            domain % chim_tn(:,e) = (/ cidx(Bm,im-1,jm,km), cidx(Bm,im+1,jm,km), cidx(Bm,im,jm,km-1), cidx(Bm,im,jm,km+1) /)
            if (im == 1)  domain % chim_cg(1,e) = cidx(Bm,im-1,jg,km)
            if (im == n1) domain % chim_cg(2,e) = cidx(Bm,im+1,jg,km)
            if (km == 1)  domain % chim_cg(3,e) = cidx(Bm,im,jg,km-1)
            if (km == n3) domain % chim_cg(4,e) = cidx(Bm,im,jg,km+1)
          case (5:6)
            domain % chim_tn(:,e) = (/ cidx(Bm,im-1,jm,km), cidx(Bm,im+1,jm,km), cidx(Bm,im,jm-1,km), cidx(Bm,im,jm+1,km) /)
            if (im == 1)  domain % chim_cg(1,e) = cidx(Bm,im-1,jm,kg)
            if (im == n1) domain % chim_cg(2,e) = cidx(Bm,im+1,jm,kg)
            if (jm == 1)  domain % chim_cg(3,e) = cidx(Bm,im,jm-1,kg)
            if (jm == n2) domain % chim_cg(4,e) = cidx(Bm,im,jm+1,kg)
        end select
      end do
      !$acc enter data copyin(domain % chim_g1, domain % chim_g2, domain % chim_rc, &
      !$acc&  domain % chim_dlo, domain % chim_dmid, domain % chim_dhi, domain % chim_slot, &
      !$acc&  domain % chim_df, domain % chim_vf, domain % chim_tn, domain % chim_cg)
    end if

    ! ---- manifold flat tables -------------------------------------------------
    ! Group maximal contiguous sorted-slot runs of type 1001 sharing (b,f,bs,fs).
    block
      integer, allocatable :: glo(:), ghi(:)
      integer :: ng
      allocate ( glo(domain % nbound), ghi(domain % nbound) )
      ng = 0; open_grp = .false.; cur = 0
      do i = 1, domain % nbound
        js = domain % bc_perm(i)
        if (domain % bc(js) % type == 1001 .and. is_local_block(domain % bc(js) % b)) then
          key = (/ domain % bc(js) % b, domain % bc(js) % f, domain % bc(js) % bs, domain % bc(js) % fs /)
          if (open_grp .and. all(key == cur) .and. i == ghi(ng) + 1) then
            ghi(ng) = i
          else
            ng = ng + 1; glo(ng) = i; ghi(ng) = i; cur = key; open_grp = .true.
          end if
        else
          open_grp = .false.
        end if
      end do
      domain % n_mfold = ng
      if (ng > 0) then
        allocate ( domain % mf_clo(ng), domain % mf_chi(ng), domain % mf_lo(ng), &
                   domain % mf_hi(ng), domain % mf_aM(ng) )
        ! size the source-cell tables
        nd = 0; maxc = 0
        do e = 1, ng
          js = domain % bc_perm(glo(e))
          Bs = domain % bc(js) % bs
          if (.not. is_local_block(Bs)) error stop 'Setup_BC_Device_Flat: manifold source block not local'
          call mf_face_bounds ( domain, Bs, domain % bc(js) % fs, i0, i1, j0, j1, k0, k1, dirs )
          ncf = (i1-i0+1)*(j1-j0+1)*(k1-k0+1)
          nd = nd + ncf; maxc = max(maxc, ncf)
        end do
        allocate ( domain % mf_cf(nd), domain % mf_af(nd), domain % mf_tmp(maxc, np+2) )
        nd = 0
        do e = 1, ng
          js = domain % bc_perm(glo(e))
          Bm = domain % bc(js) % b;  Fm = domain % bc(js) % f
          Bs = domain % bc(js) % bs; Fs = domain % bc(js) % fs
          domain % mf_lo(e) = glo(e); domain % mf_hi(e) = ghi(e)
          ! receptor-face total area: host BC_Manifold sums the FULL block face
          call mf_face_bounds ( domain, Bm, Fm, i0, i1, j0, j1, k0, k1, dirm )
          modf = mod(Fm, 2)
          domain % mf_aM(e) = 0d0
          do km = k0, k1; do jm = j0, j1; do im = i0, i1
            Fi = im; Fj = jm; Fk = km
            select case (dirm)
              case (1); Fi = im - modf
              case (2); Fj = jm - modf
              case (3); Fk = km - modf
            end select
            domain % mf_aM(e) = domain % mf_aM(e) + domain % blk(Bm) % dir(dirm) % f(Fi,Fj,Fk) % A
          end do; end do; end do
          ! source-face cells + areas, host loop order (ks,js,is; is fastest)
          call mf_face_bounds ( domain, Bs, Fs, i0, i1, j0, j1, k0, k1, dirs )
          modf = mod(Fs, 2)
          domain % mf_clo(e) = nd + 1
          do ks_ = k0, k1; do js_ = j0, j1; do is_ = i0, i1
            nd = nd + 1
            domain % mf_cf(nd) = cidx(Bs, is_, js_, ks_)
            Fi = is_; Fj = js_; Fk = ks_
            select case (dirs)
              case (1); Fi = is_ - modf
              case (2); Fj = js_ - modf
              case (3); Fk = ks_ - modf
            end select
            domain % mf_af(nd) = domain % blk(Bs) % dir(dirs) % f(Fi,Fj,Fk) % A
          end do; end do; end do
          domain % mf_chi(e) = nd
        end do
        !$acc enter data copyin(domain % mf_cf, domain % mf_af) create(domain % mf_tmp)
      end if
      deallocate ( glo, ghi )
    end block


  end subroutine Setup_BC_Device_Flat





  !> Full-face cell bounds + direction for block b face f (host BC_Manifold mapping).
  subroutine mf_face_bounds ( domain, b, f, i0, i1, j0, j1, k0, k1, dirf )
    use MOSE_Types_m
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer, intent(in)  :: b, f
    integer, intent(out) :: i0, i1, j0, j1, k0, k1, dirf
    i0 = 1; j0 = 1; k0 = 1
    i1 = domain % blk(b) % dim(1); j1 = domain % blk(b) % dim(2); k1 = domain % blk(b) % dim(3)
    select case (f)
      case (1); i1 = 1;                        dirf = 1
      case (2); i0 = domain % blk(b) % dim(1); dirf = 1
      case (3); j1 = 1;                        dirf = 2
      case (4); j0 = domain % blk(b) % dim(2); dirf = 2
      case (5); k1 = 1;                        dirf = 3
      case (6); k0 = domain % blk(b) % dim(3); dirf = 3
    end select
  end subroutine mf_face_bounds


  !----------------------------------------------------------------------------
  ! Manifold (1001) per-substage refresh: area-weighted conserved mean over the
  ! source face -> (T0, g, pstat) -> bc_t0/bc_p0/bc_pamb slots, all on device.
  ! Two passes: a parallel per-cell terms kernel (no races) and a SERIAL device
  ! reduction that sums in the host loop order -- deterministic and bit-ordered
  ! like host BC_Manifold's accumulation.
  !----------------------------------------------------------------------------
  subroutine BC_Manifold_Update ( domain )
    use MOSE_Types_m
    use MOSE_Global_m, only: np
    use MOSE_Field_m, only: P_flat, ncell_total

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: m, ncf

    do m = 1, domain % n_mfold
      ncf = domain % mf_chi(m) - domain % mf_clo(m) + 1
      call Manifold_Face_Terms ( P_flat, ncell_total, domain % mf_cf, domain % mf_af, &
                                   size(domain % mf_cf), domain % mf_clo(m), ncf,  &
                                   domain % mf_tmp, size(domain % mf_tmp,1) )
      call Manifold_Reduce ( domain % mf_tmp, size(domain % mf_tmp,1), ncf, &
                               domain % mf_aM(m), domain % mf_lo(m), domain % mf_hi(m), &
                               domain % bc_t0, domain % bc_p0, domain % bc_pamb, domain % nbound )
    end do

  end subroutine BC_Manifold_Update


  subroutine Manifold_Face_Terms ( w, ncw, cf, af, ncf_tot, clo, ncf, tmp, ntmp )
    use iso_fortran_env, only: R8 => real64
    use MOSE_Global_m, only: nsc, nprim, np, NPRIM_DEV
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, prim2cons_dev

    implicit none
    integer,  intent(in)    :: ncw, ncf_tot, clo, ncf, ntmp
    real(R8), intent(in)    :: w(ncw, nprim)
    integer,  intent(in)    :: cf(ncf_tot)
    real(R8), intent(in)    :: af(ncf_tot)
    real(R8), intent(inout) :: tmp(ntmp, np+2)
    integer  :: c, s
    real(R8) :: pr(NPRIM_DEV), cons(NPRIM_DEV), rho, Rgas, T, a

    !$acc parallel loop gang vector present(w, cf, af, tmp) private(c, s, pr, cons, rho, Rgas, T, a)
    do c = 1, ncf
      do s = 1, np
        pr(s) = w(cf(clo+c-1), s)
      end do
      call prim2cons_dev ( pr(1:np), cons(1:np) )
      a = af(clo+c-1)
      do s = 1, np
        tmp(c,s) = cons(s)*a
      end do
      call co_rotot_Rtot_dev ( pr(1:nsc), rho, Rgas )
      T = pr(np) / (rho*Rgas)
      tmp(c,np+1) = T*a
      tmp(c,np+2) = a
    end do

  end subroutine Manifold_Face_Terms


  subroutine Manifold_Reduce ( tmp, ntmp, ncf, areaM, lo, hi, bc_t0, bc_p0, bc_pamb, nbnd )
    use iso_fortran_env, only: R8 => real64
    use MOSE_Global_m, only: nsc, np, nu, nv, nw, NPRIM_DEV
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, cons2prim_dev
    use FLINT_Lib_Thermodynamic, only: f_tabT_expr, cp_tab

    implicit none
    integer,  intent(in)    :: ntmp, ncf, lo, hi, nbnd
    real(R8), intent(in)    :: tmp(ntmp, np+2)
    real(R8), intent(in)    :: areaM
    real(R8), intent(inout) :: bc_t0(nbnd), bc_p0(nbnd), bc_pamb(nbnd)
    integer  :: c, s, i, T_i, Tint(2)
    real(R8) :: mcons(NPRIM_DEV), mprim(NPRIM_DEV), mtroz, aS, g, v, rho, Rgas, T, Tdiff, cp, T0

    !$acc serial present(tmp, bc_t0, bc_p0, bc_pamb) &
    !$acc&   private(c, s, i, T_i, Tint, mcons, mprim, mtroz, aS, g, v, rho, Rgas, T, Tdiff, cp, T0)
    do s = 1, np
      mcons(s) = 0d0
    end do
    mtroz = 0d0; aS = 0d0
    do c = 1, ncf
      do s = 1, np
        mcons(s) = mcons(s) + tmp(c,s)
      end do
      mtroz = mtroz + tmp(c,np+1)
      aS    = aS    + tmp(c,np+2)
    end do
    mtroz = mtroz / aS
    do s = 1, np
      mcons(s) = mcons(s) / aS
    end do
    call cons2prim_dev ( mcons(1:np), mtroz, mprim(1:np) )
    g = sqrt( mcons(nu)**2 + mcons(nv)**2 + mcons(nw)**2 ) * aS / areaM
    v = sqrt( mprim(nu)**2 + mprim(nv)**2 + mprim(nw)**2 )
    call co_rotot_Rtot_dev ( mprim(1:nsc), rho, Rgas )
    T = mprim(np) / (rho*Rgas)
    T_i = int(T); Tdiff = T - T_i
    Tint(1) = T_i; Tint(2) = T_i + 1
    cp = f_tabT_expr ( 1, cp_tab, Tint, Tdiff )
    T0 = T + 0.5d0 * v**2 / cp
    do i = lo, hi
      bc_t0(i)   = T0
      bc_p0(i)   = g
      bc_pamb(i) = mprim(np)
    end do
    !$acc end serial

  end subroutine Manifold_Reduce


  !----------------------------------------------------------------------------
  ! Device driver: one per-block kernel launch over that block's BC faces.
  !----------------------------------------------------------------------------
  subroutine BC_Fluxes ( domain )
    use MOSE_Field_m,            only: dboff, dpi, dpj
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block, mpi_size_
    use MOSE_Lib_Metrics, only: delthe
    use MOSE_Lib_RANS, only: mcdiff_dev
    use MOSE_Field_m, only: P_flat, R_flat, ncell_total, field_is_allocated, field_alloc, &
                            pack_P_to_flat, field_update_P_flat_device, field_block_index_params
    use MOSE_Lib_Preconditioning, only: ur_f, Prec_Alloc
    use MOSE_Lib_Multirate_Base, only: mr_block_active, mr_stage_on, mrb, obj_mr, mr_rec_on, mrbl, mr_q_async
    use MOSE_Global_m, only: gc, nprim
    use MOSE_Mod_BC_Q2D, only: q2d_time, q2d_time_set

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer  :: b, i, js, qarg, nur
    real(R8) :: delthe_l, t_bc
    logical, save :: bc_tv_checked = .false., bc_has_tv_inlet = .false.
    integer  :: ncw, nrw, boff, pi, pj, wroff
    ! the record list of a multirate stage; outside one a resident placeholder (the kernels loop over the ranges)
    integer, allocatable, save :: lst0(:)
    integer, parameter :: nz6(6) = 0
    ! the list path: a resident stand-in for the per-block var-first r the record twins take (they never touch it), as
    ! large as the largest block so the explicit-shape dummies conform
    real(R8), allocatable, save :: r0(:,:)
    integer :: f, l, nmax
    ! Stream count for the per-block async fan-out (unconditional since 2026-07-11).
    ! Default 2 = the measured A30 plateau; runtime env MOSE_BC_NQ overrides.
    integer, save :: bc_nq = -1
    integer       :: nqst
    character(32) :: nqenv

    delthe_l = delthe   ! host scalar -> firstprivate into the kernel (no residency)
    if (bc_nq < 0) then
      bc_nq = 2   ! measured plateau NQ 2-4 on A30; NQ=nb over-subscribes (+0.8%). Env MOSE_BC_NQ overrides.
      call get_environment_variable("MOSE_BC_NQ", nqenv, status=nqst)
      if (nqst == 0) then
        read(nqenv, *, iostat=nqst) bc_nq
        if (nqst /= 0 .or. bc_nq <= 0) bc_nq = 2
      end if
    end if

    ! Time-varying-inlet check ONCE (BC types are static): the per-substage host scan +
    ! synchronous update device(bc_p0) is pure overhead when no inlet is time-varying
    ! (RDE/THOR). bc_p0 is then constant + already device-resident (F5 copyin) -> skip.
    if (.not. bc_tv_checked) then
      bc_has_tv_inlet = .false.
      do i = 1, domain % nbound
        js = domain % bc_perm(i)
        if (domain % bc(js) % type == 4 .and. domain % bc(js) % p0time % exists) then
          bc_has_tv_inlet = .true.; exit
        end if
      end do
      bc_tv_checked = .true.
    end if
    if (bc_has_tv_inlet) then
      t_bc = domain % time
      if ( q2d_time_set ) t_bc = q2d_time
      do i = 1, domain % nbound
        js = domain % bc_perm(i)
        ! p0(t) at the time of the stage (as the boundary data of type 667, Lib_Extrapolation Fill_Ghost_Cell): t^n + theta dt
        ! of the RK stage, or the stage time of the multirate ladder; domain%time is t^n + dt for every stage
        if (domain % bc(js) % type == 4 .and. domain % bc(js) % p0time % exists) &
          domain % bc_p0(i) = domain % bc(js) % p0time % update(t_bc)
      end do
      !$acc update device(domain % bc_p0)
    end if

    ! The BC kernel reads the var-LAST P_flat. w is kept resident +
    ! current per substage in Explicit_Step (n=1: Compute_dt step-start pack + Newstate +
    ! Fill_Ghost_Cell_local; n>1: field_sync_P_to_flat), and nothing writes the relevant w
    ! between there and here (flux only writes r), so NO pack here. boff/pi/pj are host
    ! scalars fetched per block; ncw = ncell_total.
    if (.not. field_is_allocated()) call field_alloc(domain)
    if (.not. allocated(lst0)) then
      allocate ( lst0(1) ); lst0 = 0
      !$acc enter data copyin(lst0)
    end if
    ncw = ncell_total; nrw = ncw
    call Prec_Alloc ()   ! ur_f (HLLC-PC at the boundary faces): ncell_total with prec, 1 element otherwise
    nur = size(ur_f)

    ! Manifold (1001): refresh the per-substage (T0, g, pstat) inflow scalars from the
    ! source-face mean on device (reads the current w; writes bc_t0/bc_p0/bc_pamb slots).
    ! Runs on the default queue -> covered by the blanket wait below.
    if (domain % n_mfold > 0) call BC_Manifold_Update ( domain )

    ! Ghost-fill / bc_pg device writes (default queue) must complete before the async BC
    ! kernels READ bc_pg + the resident ghosts. One blanket wait covers every prior queue.
    ! mr-async: on a list stage the fills and the BC list kernels share the stage queue (in order): no wait, unless a
    ! manifold update ran on the default queue just above.
    if ( domain % n_mfold > 0 .or. .not. mr_q_async() ) then
      !$acc wait
    end if
    ! a list level of a multirate stage (mr-list-level): the level's records of every local block, one launch per block
    ! face in face order (the per-cell order of the per-block launches), the values of the record body gathered by list
    if ( mr_rec_on() ) then
      l = obj_mr%lcur
      call bc_lists ( mrbl(:,l) )
      return
    end if
    ! the explicit step (no multirate stage): every record of every local block on the lists of bc_list_setup (the same
    ! launches as a list level; MOSE_BC_LIST=0, a shock detector or a 3D periodic seam keep the per-block launches below)
    if ( .not. mr_stage_on() ) then
      if ( bc_list_mode < 0 ) call bc_list_setup ( domain )
      if ( bc_list_mode == 1 .and. domain % nb == bla_nb .and. domain % nbound == bla_nbound ) then
        call bc_lists ( bla )
        return
      end if
    end if
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      if (domain % bc_blk_hi(b) < domain % bc_blk_lo(b)) cycle
      if (.not. mr_block_active(b)) cycle   ! multirate stage: no cell of the stage's level in this block
      call field_block_index_params(b, boff, pi, pj); wroff = boff
      qarg = mod(b-1, bc_nq) + 1     ! blocks sharing a stream serialize (still race-free; just less concurrency)
      if ( mr_stage_on() ) then      ! multirate stage: the records of the stage's level only (Mod_Multirate_Partition)
      call BC_Fluxes_blk ( PWPASS,   domain % blk(b) % r,   &
                               domain % blk(b) % fm1, domain % blk(b) % fm2, &
                               domain % blk(b) % fm3, domain % blk(b) % dl_f, &
                               domain % blk(b) % m_f, domain % blk(b) % yn,  &
                               domain % blk(b) % beta, domain % blk(b) % dim, &
                               domain % bc_bi, domain % bc_bj, domain % bc_bk, &
                               domain % bc_bf, domain % bc_typ, domain % nbound, &
                               domain % n_bf(b,1:6), domain % bc_blk_lo(b), &
                               domain % bc_mach, domain % bc_t0, domain % bc_p0, &
                               domain % bc_relfac, domain % bc_alpha, domain % bc_beta, &
                               domain % bc_pamb, domain % bc_ci, domain % bc_rans, domain % bc_tvec, &
                               domain % bc_tw, domain % bc_pg, delthe_l, qarg, &
                               mrb(b)%bcl, size(mrb(b)%bcl), mrb(b)%bcn(:,obj_mr%lcur), mrb(b)%bco(:,obj_mr%lcur), .true. )
      else
      call BC_Fluxes_blk ( PWPASS,   domain % blk(b) % r,   &
                               domain % blk(b) % fm1, domain % blk(b) % fm2, &
                               domain % blk(b) % fm3, domain % blk(b) % dl_f, &
                               domain % blk(b) % m_f, domain % blk(b) % yn,  &
                               domain % blk(b) % beta, domain % blk(b) % dim, &
                               domain % bc_bi, domain % bc_bj, domain % bc_bk, &
                               domain % bc_bf, domain % bc_typ, domain % nbound, &
                               domain % n_bf(b,1:6), domain % bc_blk_lo(b), &
                               domain % bc_mach, domain % bc_t0, domain % bc_p0, &
                               domain % bc_relfac, domain % bc_alpha, domain % bc_beta, &
                               domain % bc_pamb, domain % bc_ci, domain % bc_rans, domain % bc_tvec, &
                               domain % bc_tw, domain % bc_pg, delthe_l, qarg, &
                               lst0, 1, nz6, nz6, .false. )
      end if
    end do
    ! Drain every per-block BC queue before downstream stages (Newstate/chemistry) read r/R_flat.
    !$acc wait

  contains

    !> The launches on the record lists bl(1:6) (a multirate level, or every record for the explicit step): one per block
    !> face in face order, then the drain. r0 stands in for the per-block var-first r the record twins take.
    subroutine bc_lists ( bl )
      type(mr_bl_type), intent(in) :: bl(:)
      if ( .not. allocated(r0) ) then
        nmax = 1
        do b = 1, domain % nb
          if (is_local_block(b)) nmax = max( nmax, product(domain%blk(b)%dim + 2*gc) )
        end do
        allocate ( r0(nprim, nmax) ); r0 = 0.0_R8
        !$acc enter data copyin(r0)
      end if
      do f = 1, 6
        if ( bl(f)%n == 0 ) cycle
        call BC_Fluxes_list ( P_flat, R_flat, ncw, nrw, ur_f, nur, r0, size(r0,2), &
                              bl(f)%rec, bl(f)%rb, bl(f)%gn, bl(f)%gdl, bl(f)%gm, bl(f)%gy, &
                              bl(f)%n, domain % bc_bi, domain % bc_bj, domain % bc_bk, domain % bc_bf, domain % bc_typ, &
                              domain % nbound, domain % bc_mach, domain % bc_t0, domain % bc_p0, &
                              domain % bc_relfac, domain % bc_alpha, domain % bc_beta, &
                              domain % bc_pamb, domain % bc_ci, domain % bc_rans, domain % bc_tvec, &
                              domain % bc_tw, domain % bc_pg, delthe_l )
      end do
      if ( .not. mr_q_async() ) then   ! mr-async: the downstream kernels follow on the same queue
        !$acc wait
      end if
    end subroutine bc_lists

  end subroutine BC_Fluxes


  !> The record lists of the explicit step's boundary fluxes (bla), decided once: MOSE_BC_LIST=0 keeps the per-block
  !> launches (BC_Fluxes_blk, 3 per block, every record of the block in each); otherwise every record of every local block
  !> goes on the lists of mr_gather_bc_list -- one launch per block face over all the blocks (BC_Fluxes_list, measured
  !> 3x faster per record on the multirate stages) -- when the gather allows it: no 3D rotational periodic seam, and no
  !> shock detector (the gather takes beta of the connection records once, from the host: with the detector beta
  !> changes every stage).
  subroutine bc_list_setup ( domain )
    use MOSE_Types_m, only: MOSE_domain_type
    use MOSE_Mod_MPI, only: mpi_is_root
    use MOSE_Lib_Limiters, only: SD_space_order
    use MOSE_Lib_Shock_Detector, only: sd_dev
    use MOSE_Lib_Multirate_Base, only: mr_setup_dn, mr_gather_bc_list
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    character(32) :: env
    character(48) :: why
    integer :: st, f
    logical :: ok
    bc_list_mode = 1; why = ''
    call get_environment_variable ( "MOSE_BC_LIST", env, status = st )
    if ( st == 0 ) then
      if ( trim(adjustl(env)) == '0' ) then
        bc_list_mode = 0; why = 'MOSE_BC_LIST=0'
      end if
    end if
    if ( bc_list_mode == 1 .and. ( SD_space_order .or. sd_dev > 0 ) ) then
      bc_list_mode = 0; why = 'shock detector on'
    end if
    if ( bc_list_mode == 1 ) then
      call mr_setup_dn ( domain )
      allocate ( bla(6) )
      call mr_gather_bc_list ( domain, -1, bla, ok )
      if ( ok ) then
        bla_nb = domain % nb; bla_nbound = domain % nbound
      else
        bc_list_mode = 0; why = '3D rotational periodic records'; deallocate ( bla )
      end if
    end if
    if ( mpi_is_root ) then
      if ( bc_list_mode == 1 ) then
        write(*,'(A,6(1X,I0))') ' [BC] explicit step: boundary fluxes on record lists, records per block face', &
                                 (bla(f)%n, f = 1, 6)
      else
        write(*,'(A)') ' [BC] explicit step: boundary fluxes on per-block launches ('//trim(why)//')'
      end if
    end if
  end subroutine bc_list_setup


  !----------------------------------------------------------------------------
  ! Per-block boundary-flux kernel. One `acc data` region copies the block's P +
  ! flat geometry in and r back; one `acc parallel loop` per FACE (f=1..6) so the
  ! launch boundary acts as the host's per-face barrier — avoids the corner-cell
  ! residual race (two faces += into the same r cell) without atomics. Each face
  ! group's cells are distinct, so the loop over a single face is race-free.
  ! Geometry (normal/area/dl) is gathered here and passed to the `acc routine
  ! seq` BC twins (which take plain arrays + scalars).
  !----------------------------------------------------------------------------
  subroutine BC_Fluxes_blk ( PWARG, r, fm1, fm2, fm3, dl_f, m_f, yn, beta, n, &
                                 bc_bi, bc_bj, bc_bk, bc_bf, bc_typ, nbnd, nbf6, off, &
                                 bc_mach, bc_t0, bc_p0, bc_relfac, bc_alpha, bc_beta, &
                                 bc_pamb, bc_ci, bc_rans, bc_tvec, bc_tw, bc_pg, delthe, qb, lst, nlst, lcnt, loff, uselist )
    use iso_fortran_env, only: R8 => real64
    use MOSE_Global_m, only: gc, VISC, nprim, nsc, nrans, nu, nv, nw
    use MOSE_Lib_BC,   only: BC_Symmetry_Eul, BC_Symmetry_Visc, BC_Extrapolation, &
                             BC_Inflow_Outflow, BC_Outlet, BC_Inflow_Nozzle, BC_Connection_Eul, &
                             BC_Connection_Visc, BC_Connection_Visc_fast, BC_Connection_Eul_lean, BC_Connection_Visc_lean, &
                             conn_fast_dev, BC_Rotational_Periodic_Eul, BC_Wall_Heat_dev, &
                             BC_Wall_Temperature_dev

    implicit none
    integer,  intent(in)    :: n(3), nbnd, nbf6(6), off
    integer,  intent(in)    :: qb         ! async queue id (= block index; per-block stream, ordered within block)
    integer,  intent(in)    :: ncw, nrw, boff, pi, pj, wroff
    real(R8), intent(in)    :: w(ncw, nprim)
    real(R8), intent(inout) :: wr(nrw, nprim)   ! Res var-last: BC residual write target (threaded to BC_*_dev)
    integer,  intent(in)    :: nur
    real(R8), intent(in)    :: ur(nur)          ! Lib_Preconditioning ur_f (1 element without prec): HLLC-PC
    real(R8), intent(inout) :: r(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc)
    real(R8), intent(in)    :: fm1(4, 0:n(1), n(2), n(3))
    real(R8), intent(in)    :: fm2(4, n(1), 0:n(2), n(3))
    real(R8), intent(in)    :: fm3(4, n(1), n(2), 0:n(3))
    real(R8), intent(in)    :: dl_f(3, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc)
    real(R8), intent(in)    :: m_f(3, 3, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc)
    real(R8), intent(in)    :: yn(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc)
    real(R8), intent(in)    :: beta(1:n(1), 1:n(2), 1:n(3))
    integer,  intent(in)    :: bc_bi(nbnd), bc_bj(nbnd), bc_bk(nbnd), bc_bf(nbnd), bc_typ(nbnd)
    real(R8), intent(in)    :: bc_mach(nbnd), bc_t0(nbnd), bc_p0(nbnd), bc_relfac(nbnd)
    real(R8), intent(in)    :: bc_alpha(nbnd), bc_beta(nbnd), bc_pamb(nbnd)
    real(R8), intent(in)    :: bc_ci(nsc,nbnd), bc_rans(max(nrans,1),nbnd), bc_tvec(3,nbnd)
    real(R8), intent(in)    :: bc_tw(nbnd)
    real(R8), intent(in)    :: bc_pg(nprim, 6, nbnd)
    real(R8), intent(in)    :: delthe   ! wedge angle (firstprivate scalar, type-2)
    integer,  intent(in)    :: nlst, lst(nlst)        ! multirate stage: the records of the stage's level (Lib_Multirate_Base bcl)
    integer,  intent(in)    :: lcnt(6), loff(6)       ! ... per launch g: lst(loff(g)+1 : loff(g)+lcnt(g))
    logical,  intent(in)    :: uselist                ! t: the launches sweep the list; f: every record of their range
    ! Local
    integer  :: g, ng, glo(6), ghi(6), idx, lo, hi, e, nl_, ol_, ierr
    integer  :: Im, Jm, Km, Fm, dir, modfm, gstep, Fi, Fj, Fk, Ii, Ij, Ik
    integer  :: Ig_, Jg_, Kg_, o1_, o2_, o3_, o4_     ! 5r-b: derived here, not inside the routines
    real(R8) :: normal(3), area, Bound_Dl, Int_Dl, dlc4_(4)
    real(R8) :: Mrec_(3,3), Mgho_(3,3), ynrec_, yngho_, beta_   ! 5r-c: handed to the shared apply
    ! Dims as scalars -- a dummy array referenced inside the acc region
    ! is implicitly copy()d per launch (upload+download+sync); scalars ride the
    ! launch parameters for free.
    integer :: nb1_, nb2_, nb3_
    nb1_ = n(1); nb2_ = n(2); nb3_ = n(3)
    ! BC FUSE opposite face pairs {1,2}{3,4}{5,6}: when the direction is non-degenerate
    ! (n_d>=2), the two faces write DISJOINT planes (face 2d-1 -> plane 1, face 2d -> plane
    ! n_d), so they run concurrently in ONE launch (6->3 launches, ~2x SM fill) with NO
    ! overlapping r write. The launches are default-SYNC -> hard barrier between pairs -> a
    ! corner on faces {1,3,5} accumulates in launch order 1->3->5, byte-identical to the
    ! 6-launch version (bit-exact + deterministic; golden 0.0). DEGENERATE direction (n_d=1)
    ! collapses opposite faces onto one plane -> would alias -> falls back to a per-face launch.
    ! the launches: opposite faces fused (bc_launch_ranges); a multirate stage sweeps only its level's records of each range,
    ! in the same launches and order (a range with none is not launched)
    call bc_launch_ranges ( n, nbf6, off, ng, glo, ghi )
    do g = 1, ng
      lo = glo(g); hi = ghi(g)
      if ( uselist ) then
      nl_ = lcnt(g); ol_ = loff(g)
      if ( nl_ == 0 ) cycle
      !$acc parallel loop gang vector BCASYNC present(w, wr, ur, r, fm1, fm2, fm3, dl_f, m_f, yn, beta, lst, &
      !$acc&   bc_bi, bc_bj, bc_bk, bc_bf, bc_typ, &
      !$acc&   bc_mach, bc_t0, bc_p0, bc_relfac, bc_alpha, bc_beta, bc_pamb, bc_ci, bc_rans, bc_tvec, bc_tw, bc_pg) &
      !$acc&   private(idx, Im, Jm, Km, Fm, dir, modfm, gstep, Fi, Fj, Fk, Ii, Ij, Ik, ierr, &
      !$acc&           normal, area, Bound_Dl, Int_Dl, Ig_, Jg_, Kg_, o1_, o2_, o3_, o4_, dlc4_, &
      !$acc&           Mrec_, Mgho_, ynrec_, yngho_, beta_)
      do e = 1, nl_
        idx = lst(ol_ + e)
#include "bc_record_body.inc"
      end do
      else
      !$acc parallel loop gang vector BCASYNC present(w, wr, ur, r, fm1, fm2, fm3, dl_f, m_f, yn, beta, &
      !$acc&   bc_bi, bc_bj, bc_bk, bc_bf, bc_typ, &
      !$acc&   bc_mach, bc_t0, bc_p0, bc_relfac, bc_alpha, bc_beta, bc_pamb, bc_ci, bc_rans, bc_tvec, bc_tw, bc_pg) &
      !$acc&   private(Im, Jm, Km, Fm, dir, modfm, gstep, Fi, Fj, Fk, Ii, Ij, Ik, ierr, &
      !$acc&           normal, area, Bound_Dl, Int_Dl, Ig_, Jg_, Kg_, o1_, o2_, o3_, o4_, dlc4_, &
      !$acc&           Mrec_, Mgho_, ynrec_, yngho_, beta_)
      do idx = lo, hi
#include "bc_record_body.inc"
      end do
      end if
    end do


  end subroutine BC_Fluxes_blk


  !> The boundary fluxes of a multirate list level over the records of one block face of every local block (Lib_Multirate_Base
  !> mr_bl_type): per record the block's flat offset, row lengths and dims from the device tables, the values of
  !> bc_record_body.inc from the gather, then bc_record_apply.inc -- the per-block kernel's work on the same record.
  subroutine BC_Fluxes_list ( w, wr, ncw, nrw, ur, nur, r, nr, rec, rb, gn, gdl, gm, gy, ne, &
                              bc_bi, bc_bj, bc_bk, bc_bf, bc_typ, nbnd, &
                              bc_mach, bc_t0, bc_p0, bc_relfac, bc_alpha, bc_beta, &
                              bc_pamb, bc_ci, bc_rans, bc_tvec, bc_tw, bc_pg, delthe )
    use MOSE_Lib_Multirate_Base, only: mr_q
    use iso_fortran_env, only: R8 => real64
    use MOSE_Global_m, only: gc, VISC, nprim, nsc, nrans, nu, nv, nw
    use MOSE_Lib_BC,   only: BC_Symmetry_Eul, BC_Symmetry_Visc, BC_Extrapolation, &
                             BC_Inflow_Outflow, BC_Outlet, BC_Inflow_Nozzle, BC_Connection_Eul, &
                             BC_Connection_Visc, BC_Connection_Visc_fast, BC_Connection_Eul_lean, BC_Connection_Visc_lean, &
                             conn_fast_dev, BC_Rotational_Periodic_Eul, BC_Wall_Heat_dev, &
                             BC_Wall_Temperature_dev
    use MOSE_Field_m,  only: dboff, dpi, dpj
    use MOSE_Lib_Multirate_Base, only: mr_dn
    implicit none
    integer,  intent(in)    :: ncw, nrw, nur, nr, ne, nbnd
    real(R8), intent(in)    :: w(ncw, nprim)
    real(R8), intent(inout) :: wr(nrw, nprim)
    real(R8), intent(in)    :: ur(nur)
    real(R8), intent(inout) :: r(nprim, nr)
    integer,  intent(in)    :: rec(ne), rb(ne)
    real(R8), intent(in)    :: gn(4,ne), gdl(6,ne), gm(3,3,2,ne), gy(3,ne)
    integer,  intent(in)    :: bc_bi(nbnd), bc_bj(nbnd), bc_bk(nbnd), bc_bf(nbnd), bc_typ(nbnd)
    real(R8), intent(in)    :: bc_mach(nbnd), bc_t0(nbnd), bc_p0(nbnd), bc_relfac(nbnd)
    real(R8), intent(in)    :: bc_alpha(nbnd), bc_beta(nbnd), bc_pamb(nbnd)
    real(R8), intent(in)    :: bc_ci(nsc,nbnd), bc_rans(max(nrans,1),nbnd), bc_tvec(3,nbnd)
    real(R8), intent(in)    :: bc_tw(nbnd)
    real(R8), intent(in)    :: bc_pg(nprim, 6, nbnd)
    real(R8), intent(in)    :: delthe
    integer  :: e, idx, b, boff, pi, pj, wroff, nb1_, nb2_, nb3_, Im, Jm, Km, Fm, ierr
    real(R8) :: normal(3), area, Bound_Dl, Int_Dl, dlc4_(4)
    real(R8) :: Mrec_(3,3), Mgho_(3,3), ynrec_, yngho_, beta_
    !$acc parallel loop gang vector async(mr_q) present(w, wr, ur, r, rec, rb, gn, gdl, gm, gy, dboff, dpi, dpj, mr_dn, &
    !$acc&   bc_bi, bc_bj, bc_bk, bc_bf, bc_typ, &
    !$acc&   bc_mach, bc_t0, bc_p0, bc_relfac, bc_alpha, bc_beta, bc_pamb, bc_ci, bc_rans, bc_tvec, bc_tw, bc_pg) &
    !$acc&   private(idx, b, boff, pi, pj, wroff, nb1_, nb2_, nb3_, Im, Jm, Km, Fm, ierr, &
    !$acc&           normal, area, Bound_Dl, Int_Dl, dlc4_, Mrec_, Mgho_, ynrec_, yngho_, beta_)
    do e = 1, ne
      idx = rec(e); b = rb(e)
      boff = dboff(b); pi = dpi(b); pj = dpj(b); wroff = boff
      nb1_ = mr_dn(1,b); nb2_ = mr_dn(2,b); nb3_ = mr_dn(3,b)
      Im = bc_bi(idx); Jm = bc_bj(idx); Km = bc_bk(idx); Fm = bc_bf(idx)
      normal(1) = gn(1,e); normal(2) = gn(2,e); normal(3) = gn(3,e); area = gn(4,e)
      Bound_Dl = gdl(1,e); Int_Dl = gdl(2,e)
      dlc4_(1) = gdl(3,e); dlc4_(2) = gdl(4,e); dlc4_(3) = gdl(5,e); dlc4_(4) = gdl(6,e)
      Mrec_ = gm(:,:,1,e); Mgho_ = gm(:,:,2,e)
      ynrec_ = gy(1,e); yngho_ = gy(2,e); beta_ = gy(3,e)
#define dl_f r
#include "bc_record_apply.inc"
#undef dl_f
    end do
  end subroutine BC_Fluxes_list


  !> The launches of BC_Fluxes_blk for a block: the record ranges glo(g)..ghi(g), g = 1..ng, in launch order. The records
  !> of the block are grouped by face (nbf6, from off); the two faces of a direction share one launch when the direction is
  !> not degenerate (n_d >= 2: they write disjoint planes), else each face is a launch. Also the multirate record lists
  !> (Mod_Multirate_Partition) follow these ranges.
  subroutine bc_launch_ranges ( n, nbf6, off, ng, glo, ghi )
    integer, intent(in)  :: n(3), nbf6(6), off
    integer, intent(out) :: ng, glo(6), ghi(6)
    integer :: f, lo, hi, base, flo
    ng = 0; glo = 0; ghi = -1; base = off - 1; flo = 0
    do f = 1, 6
      lo   = base + 1
      hi   = base + nbf6(f)
      base = hi
      if (mod(f,2) == 1) then            ! odd face: start of a candidate opposite-face pair
        if (n((f+1)/2) >= 2) then          ! non-degenerate -> defer; fuse at the even face
          flo = lo
          cycle
        end if
      else                               ! even face: close a deferred fused pair
        if (n(f/2) >= 2) lo = flo           ! extend the range back over the odd face
      end if
      if (hi < lo) cycle
      ng = ng + 1; glo(ng) = lo; ghi(ng) = hi
    end do
  end subroutine bc_launch_ranges


end module MOSE_Mod_BC
