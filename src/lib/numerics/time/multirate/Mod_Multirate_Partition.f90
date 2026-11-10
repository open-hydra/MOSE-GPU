!> Multirate: the partition of the cells into levels, run by Compute_dt after the reduction of the step and before
!> Set_Global_dt (doc/archive/MULTIRATE_PLAN.md of branch multirate, sec. 1.2 and 2):
!>   dtcap <- dtlocal (each cell's own explicit limit), s = mr-safety (mr-wall-safety on the wall rows),
!>   lev   =  clamp( ceil(log2( s H / dtcap )), 0, L ),         H = mr-ratio * dt_min,
!>   R1    :  the levels of two face neighbours differ by at most one (inside a block),
!>   R2    :  across a block connection the two rows on each side take the largest of the four levels, so a connection
!>            face never separates two levels (its flux is computed by both sides from ghosts; it is not telescoped),
!>   cap   :  H <= min over the cells of dtcap 2**lev / s     (a cell the levels do not resolve shortens the macro step).
!> Between two full partitions (at least every mr-repart macro steps) the levels are frozen and only the cap is refreshed
!> from the new dtcap (a light step); the full partition runs again when the cap falls below mr-repart-hcap of the step
!> the levels were built for.
module MOSE_Mod_Multirate_Partition
  use iso_fortran_env, only: R8 => real64
  use MOSE_Lib_Multirate_Base
  implicit none
  private
  public :: MR_Partition, MR_Set_Local_dt, MR_Write_State, MR_Read_State, MR_Rebuild_For_H

  ! mr-restart (open-hydra MOSE feature/multirate-integrator F2, by persisting the state): the partition state read at a
  ! restart (MR_Read_State), applied by the first MR_Partition (mr_restore_state) before its usual light/full decision
  logical :: rs_pending = .false.
  integer :: rs_L = 0, rs_nsince = 0, rs_nlight = 0, rs_nfull = 0
  real(R8) :: rs_H = 0.0_R8, rs_Hreq = 0.0_R8, rs_Hcap = 0.0_R8, rs_Hbuilt = 0.0_R8, rs_Hreqb = 1.0_R8, rs_ratio_eff = 16.0_R8
  type rs_lev_type
    integer, allocatable :: lev(:,:,:)
  end type rs_lev_type
  type(rs_lev_type), allocatable :: rs_blk(:)

  integer, parameter :: WALL_TYPES(7) = [5, 6, 8, 9, 10, 12, 13]   !< the viscous walls (MOSE_open 301-309)
  ! mr-list-dev scratch on the device (Mod_Multirate_Partition mr_build_lists_dev): flags, positions, chunk sums, face tables
  integer, parameter :: LD_CH = 1024
  integer, allocatable :: ld_flg(:), ld_pos(:), ld_p1(:), ld_p2(:), ld_p3(:), ld_csum(:), ld_coff(:)
  integer :: ld_q = -1   !< the async queue of the list kernels (above the block and stage queues); -1 = synchronous
  ! mr-build-fast (mr_build_fast_on): the scan in sub-chunks of LD_SC (LD_NS per chunk) with their sums ld_ssum (lt_ss for
  ! the thermo lists); the counts of every block of a level (ld_cnt(7, nb)) and the changed cells of the nesting rule (ld_nc)
  ! summed on the device and read back once
  integer, parameter :: LD_SC = 32, LD_NS = LD_CH / LD_SC
  integer, allocatable :: ld_ssum(:), ld_cnt(:,:), ld_nc(:)
  ! mr-thermo-exact scratch (mr_build_thermo_list_dev): the read mark of every flat cell, the interior flag, flags and positions
  integer, allocatable :: lt_mark(:), lt_int(:), lt_fa(:), lt_pa(:), lt_fb(:), lt_pb(:), lt_cs(:), lt_co(:), lt_ss(:)
  ! mr-ghost-exact (mr_mark_reads_dev): 1 on the flat cells whose primitives a stage of the level reads; host copy for the
  ! ghost-record choice of mr_build_ghost_list
  integer, allocatable :: rk_mark(:)
  logical :: rk_ready = .false.
  ! the census of every rank (the root prints it with the census line, multi-rank runs): cells per level, column = rank
  real(R8), allocatable :: mr_histr(:,:)

contains


  subroutine MR_Partition ( domain, config )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block, mpi_is_root, mpi_allreduce_min_r8, mpi_allreduce_sum_r8_array
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    type(config_type),      intent(in)    :: config
    logical, save :: setup_done = .false.
    integer  :: b, L, ipass, nchg, rnd, l_, r_
    integer(8) :: c0_, c1_, cr_, cp_
    real(R8) :: Hreq, Htry, hcap, tmp, sump, hist(0:MR_LMAX_HARD), cnt(2), smax, dtcmax
    integer  :: Lh

    call system_clock ( c0_, cr_ )
    if ( .not. setup_done ) then
      call MR_Setup ( domain, config )
      setup_done = .true.
    end if
    if ( rs_pending ) call mr_restore_state ( domain )   ! mr-restart: the partition of the step the restart file holds
    obj_mr%dtmin = domain % dtglobal
    Hreq = obj_mr%ratio_eff * obj_mr%dtmin   ! mr-ratio, or the ratio of the last full partition with mr-depth = auto

    ! light step: frozen levels, the cap from the new limits. They hold while their cap keeps at least repart_hcap of the
    ! step they gave at the full partition, scaled with the request since (a CFL ramp raises the request every step)
    if ( obj_mr%have_part .and. obj_mr%nsince < obj_mr%repart ) then
      call mr_store_dtcap ( domain )
      hcap = mr_cap ( domain )
      Htry = min( Hreq, obj_mr%Hbuilt * Hreq / obj_mr%Hreqb )
      if ( hcap >= obj_mr%repart_hcap * Htry ) then
        obj_mr%nsince = obj_mr%nsince + 1; obj_mr%nlight = obj_mr%nlight + 1
        obj_mr%Hreq = Hreq; obj_mr%Hcap = hcap
        obj_mr%H = min( Hreq, hcap )
        if ( obj_mr%drive ) domain % dtglobal = obj_mr%H
        if ( mpi_is_root .and. mod(domain%iter, max(obj_mr%nprint,1)) == 0 ) &
          write(*,'(A,I0,A,ES10.3,A,F7.3,A,F6.2,A,I0)') ' [MR] iter ', domain%iter, '  light  dt_min ', obj_mr%dtmin, &
            '  H_cap/H_req ', hcap/Hreq, '  R_eff ', obj_mr%H/obj_mr%dtmin, '  since full ', obj_mr%nsince
        call system_clock ( c1_ )
        obj_mr%tpart = obj_mr%tpart + real(c1_ - c0_, R8) / real(cr_, R8)
        return
      end if
    end if

    ! full partition
    obj_mr%nsince = 0; obj_mr%nfull = obj_mr%nfull + 1
    cp_ = c0_
    smax = max( obj_mr%sigma, obj_mr%sigma_wall )
    L = 0
    if ( smax * obj_mr%ratio > 1.0_R8 + 1.0e-12_R8 ) L = ceiling( log(smax*obj_mr%ratio) / log(2.0_R8) - 1.0e-12_R8 )
    if ( obj_mr%lmax >= 0 ) L = min( L, obj_mr%lmax )
    L = max( 0, min(L, MR_LMAX_HARD) )
    obj_mr%L = L
    call mr_store_dtcap ( domain )
    ! mr-depth = auto (open-hydra MOSE feature/multirate-integrator, Depth_Cap / Htry_Of, Mod_Multirate.f90 3898-3945 there): the
    ! depth from the octaves of the admissible steps, L^ = min( L_max, floor( log2( dtcap_max / dt_min ) ) ), and the step
    ! H = 2**L^ dt_min / mr-safety: the slowest cells are on level 0 (H <= dtcap_max / s) and the finest sub-step stays
    ! dt_min / s. L_max = mr-lmax, at most 8 (MR_LMAX_HARD; the branch default is 10). Kept until the next full partition.
    if ( obj_mr%depth_auto ) then
      dtcmax = mr_dtcap_max ( domain )
      Lh = 0
      if ( dtcmax > obj_mr%dtmin ) Lh = int( floor( log( dtcmax / obj_mr%dtmin ) / log(2.0_R8) + 1.0e-12_R8 ) )
      Lh = max( 0, min( Lh, MR_LMAX_HARD ) )
      if ( obj_mr%lmax >= 0 ) Lh = min( Lh, obj_mr%lmax )
      obj_mr%ratio_eff = 2.0_R8 ** Lh / obj_mr%sigma
      Hreq = obj_mr%ratio_eff * obj_mr%dtmin
      L = 0
      if ( smax * obj_mr%ratio_eff > 1.0_R8 + 1.0e-12_R8 ) &
        L = ceiling( log(smax*obj_mr%ratio_eff) / log(2.0_R8) - 1.0e-12_R8 )
      if ( obj_mr%lmax >= 0 ) L = min( L, obj_mr%lmax )
      L = max( 0, min(L, MR_LMAX_HARD) )
      obj_mr%L = L
    end if
    Htry = Hreq
    passes: do ipass = 1, 3
      ! raw levels for the step Htry
      sump = 0.0_R8; hist = 0.0_R8
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        call mr_levels_blk ( mrb(b)%dtcap, mrb(b)%cls, mrb(b)%lev, domain%blk(b)%dim, Htry, obj_mr%sigma, &
                             obj_mr%sigma_wall, L, sump, obj_mr%margin, obj_mr%kappa )
      end do
      cnt(1) = sump; cnt(2) = real(mr_ncell(domain), R8)
      call mpi_allreduce_sum_r8_array ( cnt, 2 )
      obj_mr%Cp_raw = cnt(1) / max(cnt(2), 1.0_R8)
      ! the rules: R1 (L passes), then rounds of R2 + R1 until nothing moves
      obj_mr%nrounds = 0; obj_mr%nseam = 0
      if ( L > 0 ) then
        nchg = mr_nest ( domain, L )
        do rnd = 1, 2*L + 2
          nchg = mr_seam ( domain )
          if ( rnd == 1 ) obj_mr%nseam = nchg
          nchg = nchg + mr_nest ( domain, L )
          cnt(1) = real(nchg, R8); cnt(2) = 0.0_R8
          call mpi_allreduce_sum_r8_array ( cnt, 2 )
          obj_mr%nrounds = rnd
          if ( cnt(1) < 0.5_R8 ) exit
        end do
        if ( cnt(1) > 0.5_R8 ) error stop 'multirate: the level rules (R1 nesting, R2 connections) did not converge'
      end if
      hcap = mr_cap ( domain )
      obj_mr%npass = ipass
      if ( hcap >= Htry * ( 1.0_R8 - 1.0e-12_R8 ) ) exit passes
      Htry = hcap      ! a cell the levels do not resolve: rebuild them for the step the cap allows
    end do passes

    ! mr-octaves (diagnostic, prints only): the admissible steps of the cells in octaves above the smallest, and the ideal
    ! work of other macro steps, at the first full partition and every 10th after it
    if ( obj_mr%octaves .and. ( obj_mr%nfull == 1 .or. mod(obj_mr%nfull, 10) == 0 ) ) call mr_octave_census ( domain )
    call mr_build_from_levels ( domain, L, cp_, reuse = .true. )   ! census, finest populated level, boxes, lists, halo
    if ( mpi_is_root .and. obj_mr%profile ) write(*,'(A,I0,A,I0,A,I0,A,L1)') ' [MR] iter ', domain%iter, &
      '  full partition ', obj_mr%nfull, '  cells whose level changed since the last build ', obj_mr%nchg, &
      '  structures kept ', ( obj_mr%nchg == 0 .and. mr_reuse_on() )
    tmp = sum(obj_mr%hist)

    obj_mr%Hreq = Hreq; obj_mr%Hcap = hcap
    obj_mr%H = min( Htry, hcap ); obj_mr%Hbuilt = obj_mr%H; obj_mr%Hreqb = Hreq
    obj_mr%have_part = .true.
    if ( obj_mr%drive ) domain % dtglobal = obj_mr%H

    if ( mpi_is_root .and. ( mod(domain%iter, max(obj_mr%nprint,1)) == 0 .or. obj_mr%nfull == 1 ) ) then
      write(*,'(A,I0,A,ES10.3,A,F6.1,A,I0,A,F7.3,A,F6.2,A,F7.4,A,F7.4,A,I0,A,I0)') ' [MR] iter ', domain%iter, &
        '  dt_min ', obj_mr%dtmin, '  R_req ', obj_mr%ratio_eff, '  L ', L, '  H_cap/H_req ', hcap/Hreq, &
        '  R_eff ', obj_mr%H/obj_mr%dtmin, '  C_p ', obj_mr%Cp, '  C_p_raw ', obj_mr%Cp_raw, '  passes ', obj_mr%npass, &
        '  rule rounds ', obj_mr%nrounds
      write(*,'(A)',advance='no') ' [MR]   cells per level %:'
      do l_ = 0, L
        write(*,'(A,I0,A,F8.4)',advance='no') '  L', l_, ' ', 100.0_R8*obj_mr%hist(l_)/max(tmp,1.0_R8)
      end do
      write(*,'(A)') ''
      if ( allocated(mr_histr) ) then   ! multi-rank: each rank's cells per level and its work (sum of cells x 2^level)
        do r_ = 0, size(mr_histr, 2) - 1
          write(*,'(A,I0,A)',advance='no') ' [MR]   rank ', r_, ' cells per level:'
          do l_ = 0, L
            write(*,'(A,I0,A,I0)',advance='no') '  L', l_, ' ', nint(mr_histr(l_, r_))
          end do
          write(*,'(A,ES11.4)') '   work ', sum( (/ ( mr_histr(l_, r_) * 2.0_R8**l_, l_ = 0, L ) /) )
        end do
      end if
      write(*,'(A,I0,A,I0,A,I0)') ' [MR]   connection cells raised by R2 (round 1, this rank) ', obj_mr%nseam, &
        '  interior faces |dl| > 1 after the rules ', obj_mr%nviol, '  interface faces (this rank) ', obj_mr%nifc
      write(*,'(A,I0,A,I0,A,I0)') ' [MR]   mr-reuse: cells whose level changed since the last build (all ranks) ', &
        obj_mr%nchg, '  full partitions that kept the structures ', obj_mr%nreuse, ' of ', obj_mr%nfull
      if ( obj_mr%nviol > 0 ) write(*,'(A)') ' [MR]   WARNING: faces jump more than one level -- the halo and the ' // &
        'telescoping assume |dl| <= 1'
    end if
    call system_clock ( c1_ )
    obj_mr%tpart = obj_mr%tpart + real(c1_ - c0_, R8) / real(cr_, R8)
    obj_mr%tfull = obj_mr%tfull + real(c1_ - c0_, R8) / real(cr_, R8)
  end subroutine MR_Partition

  !> The structures of the levels mrb%lev (device copy current): the census (hist, C_p, face jumps), the finest populated
  !> level (L trimmed), the host copy of the levels, the boxes, the inter-rank masks, the interface faces, the boundary
  !> lists, the face and cell lists, the halo lists. Run by the full partition and by a restart (mr_restore_state).
  !> reuse (the full partition; mr-reuse, MOSE_MR_REUSE=0 turns it off): when no cell of any rank changed level since the
  !> last build and the ladder is the same, the structures are kept -- they depend on the levels and the geometry only, so
  !> the stages are the same to the bit (open-hydra MOSE feature/multirate-integrator, Build_Map, Mod_Multirate.f90
  !> 4039-4057 there: Build_Lists only if a cell changed level in any process).
  subroutine mr_build_from_levels ( domain, L, cp_, reuse )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block, mpi_allreduce_sum_r8_array, mpi_size_, mpi_rank_
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer,    intent(inout) :: L
    integer(8), intent(inout) :: cp_
    logical, optional, intent(in) :: reuse
    integer  :: b
    real(R8) :: tmp, sump, hist(0:MR_LMAX_HARD), hb(0:MR_LMAX_HARD), cnt(2)
    ! census of the final levels
    sump = 0.0_R8; hist = 0.0_R8
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      hb = 0.0_R8
      call mr_census_blk ( mrb(b)%lev, domain%blk(b)%dim, hb, sump )
      hist = hist + hb   ! whole numbers: the same sums as accumulating in place
      if ( mr_block_census_on() ) call mr_block_census_write ( domain, b, hb )
    end do
    if ( mpi_size_ > 1 ) then   ! the census of every rank, for the root's print (load of each level per rank)
      if ( .not. allocated(mr_histr) ) allocate ( mr_histr(0:MR_LMAX_HARD, 0:mpi_size_-1) )
      mr_histr = 0.0_R8; mr_histr(:, mpi_rank_) = hist
      call mpi_allreduce_sum_r8_array ( mr_histr, (MR_LMAX_HARD+1)*mpi_size_ )
    end if
    call mpi_allreduce_sum_r8_array ( hist, MR_LMAX_HARD+1 )
    tmp = sum(hist)
    cnt(1) = sump; cnt(2) = 0.0_R8
    call mpi_allreduce_sum_r8_array ( cnt, 2 )
    obj_mr%hist = hist; obj_mr%Cp = cnt(1) / max(tmp, 1.0_R8)
    obj_mr%nviol = mr_face_jumps ( domain )
    ! the ladder stops at the finest populated level (the rules leave no empty level between populated ones): an empty
    ! finer level would only cost empty stages
    do while ( L > 0 )
      if ( hist(L) > 0.5_R8 ) exit
      L = L - 1
    end do
    obj_mr%L = L

    ! what the driver needs: the host levels, the boxes, the interface faces, the halo lists
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$acc update host(mrb(b)%lev)
    end do
    ! mr-reuse: the cells whose level differs from the build the structures hold (cnt(2) > 0: a rank without a build, or
    ! another ladder)
    obj_mr%nchg = -1
    if ( present(reuse) ) then
      cnt = 0.0_R8
      if ( obj_mr%Lbuilt /= L ) cnt(2) = 1.0_R8
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        if ( .not. allocated(mrb(b)%levb) ) then
          cnt(2) = 1.0_R8
        else
          cnt(1) = cnt(1) + real( count( mrb(b)%lev /= mrb(b)%levb ), R8 )
        end if
      end do
      call mpi_allreduce_sum_r8_array ( cnt, 2 )
      if ( cnt(2) < 0.5_R8 ) obj_mr%nchg = nint(cnt(1))
      if ( reuse .and. mr_reuse_on() .and. obj_mr%nchg == 0 ) then
        obj_mr%nreuse = obj_mr%nreuse + 1
        if ( obj_mr%profile ) call mr_ptoc ( obj_mr%tpb(1), cp_ )
        return
      end if
    end if
    if ( obj_mr%profile ) call mr_ptoc ( obj_mr%tpb(1), cp_ )   ! mr-profile: the parts of the full partitions
    call mr_build_boxes ( domain, L )
    if ( obj_mr%profile ) call mr_ptoc ( obj_mr%tpb(2), cp_ )
    if ( obj_mr%drive ) call mr_build_xch_masks ( domain, L )   ! mr-xch: the inter-rank entries each level reads
    if ( obj_mr%drive ) call mr_build_iface ( domain )
    if ( obj_mr%drive ) call mr_build_bclist ( domain, L )
    if ( obj_mr%profile ) call mr_ptoc ( obj_mr%tpb(3), cp_ )
    if ( obj_mr%drive .and. obj_mr%listlev >= 0 ) call mr_build_lists ( domain, L, cp_ )
    if ( obj_mr%drive ) call mr_build_halo ( domain, L )
    if ( obj_mr%profile ) call mr_ptoc ( obj_mr%tpb(8), cp_ )
    if ( obj_mr%drive .and. obj_mr%listlev >= 0 ) call mr_build_halo_near ( domain, L )
    if ( obj_mr%profile ) call mr_ptoc ( obj_mr%tpb(9), cp_ )
    ! mr-reuse: the levels these structures are built from
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      if ( .not. allocated(mrb(b)%levb) ) allocate ( mrb(b)%levb, mold = mrb(b)%lev )
      mrb(b)%levb = mrb(b)%lev
    end do
    obj_mr%Lbuilt = L
  end subroutine mr_build_from_levels

  !> MOSE_MR_BLOCK_CENSUS=1 (diagnostic, default off): every build of the level structures appends the cells per level of
  !> each local block to OUTPUT/mr_block_census_r<rank>.txt (full partition number, iteration, block, rank, levels 0..8):
  !> the per-block load of the levels, for a block-to-rank map that balances each level (INPUT/block_owner.txt).
  logical function mr_block_census_on ()
    logical, save :: checked = .false., on = .false.
    character(len=8) :: ev
    integer :: st
    if ( .not. checked ) then
      checked = .true.
      call get_environment_variable ( 'MOSE_MR_BLOCK_CENSUS', ev, status = st )
      if ( st == 0 ) on = ( ev(1:1) == '1' .or. ev(1:1) == 'y' .or. ev(1:1) == 'Y' .or. ev(1:1) == 't' .or. ev(1:1) == 'T' )
    end if
    mr_block_census_on = on
  end function mr_block_census_on

  subroutine mr_block_census_write ( domain, b, hb )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: mpi_rank_
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer,  intent(in) :: b
    real(R8), intent(in) :: hb(0:MR_LMAX_HARD)
    integer, save :: u = 0
    logical, save :: opened = .false.
    character(len=64) :: fn
    if ( .not. opened ) then   ! newunit= gives a negative unit number
      opened = .true.
      write(fn,'(A,I0,A)') 'OUTPUT/mr_block_census_r', mpi_rank_, '.txt'
      open ( newunit = u, file = trim(fn), status = 'replace' )
      write(u,'(A)') '# full_partition iter block rank cells_L0 .. cells_L8'
    end if
    write(u,'(I6,I9,I5,I4,9I9)') obj_mr%nfull, domain%iter, b, mpi_rank_, nint(hb)
    flush(u)
  end subroutine mr_block_census_write

  !> mr-reuse (MOSE_MR_REUSE, default on): a full partition with the map of the last build keeps its structures.
  logical function mr_reuse_on ()
    logical, save :: checked = .false., on = .true.
    character(len=8) :: env
    integer :: st, ln
    if ( .not. checked ) then
      checked = .true.
      call get_environment_variable ( 'MOSE_MR_REUSE', env, ln, st )
      if ( st == 0 .and. ln > 0 ) then
        if ( env(1:1) == '0' .or. env(1:1) == 'n' .or. env(1:1) == 'N' .or. env(1:1) == 'f' .or. env(1:1) == 'F' ) on = .false.
      end if
    end if
    mr_reuse_on = on
  end function mr_reuse_on


  !> mr-nonlocal = redo (open-hydra MOSE feature/multirate-integrator, Mod_Multirate.f90 737-757 there): after a violation of
  !> condition (4), dtcap := min(dtcap, dmin) (the smallest admissible step the certificate saw) and the map rebuilt for the
  !> step H of the attempt; the levels only rise (max with the current map) and the depth can grow. Returns the cap of the new
  !> map; the caller halves H while it is above it. Counted as a full partition (the light steps follow from it).
  subroutine MR_Rebuild_For_H ( domain, H, hcap, deeper )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block, mpi_allreduce_min_r8, mpi_allreduce_sum_r8_array
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    real(R8), intent(in)  :: H
    real(R8), intent(out) :: hcap
    logical,  intent(out) :: deeper
    integer    :: b, L, Lold, nchg, rnd, i, j, k, n1, n2, n3
    integer(8) :: cp_
    real(R8)   :: smax, sump, dloc, dg, cnt(2)
    Lold = obj_mr%L
    dloc = huge(1.0_R8)
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      n1 = domain%blk(b)%dim(1); n2 = domain%blk(b)%dim(2); n3 = domain%blk(b)%dim(3)
      associate ( dc => mrb(b)%dtcap, dm => mrb(b)%dmin, lv => mrb(b)%lev, l2 => mrb(b)%lev2 )
      !$acc parallel loop collapse(3) gang vector present(dc, dm, lv, l2) reduction(min:dloc)
      do k = 1, n3
      do j = 1, n2
      do i = 1, n1
        dc(i,j,k) = min( dc(i,j,k), dm(i,j,k) )
        l2(i,j,k) = lv(i,j,k)
        dloc = min( dloc, dc(i,j,k) )
      end do
      end do
      end do
      end associate
    end do
    call mpi_allreduce_min_r8 ( dloc, dg )
    smax = max( obj_mr%sigma, obj_mr%sigma_wall )
    L = 0
    if ( smax * H / dg > 1.0_R8 + 1.0e-12_R8 ) L = ceiling( log(smax*H/dg) / log(2.0_R8) - 1.0e-12_R8 )
    L = max( L, Lold )
    if ( obj_mr%lmax >= 0 ) L = min( L, max(obj_mr%lmax, Lold) )
    L = max( 0, min(L, MR_LMAX_HARD) )
    sump = 0.0_R8
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      call mr_levels_blk ( mrb(b)%dtcap, mrb(b)%cls, mrb(b)%lev, domain%blk(b)%dim, H, obj_mr%sigma, &
                           obj_mr%sigma_wall, L, sump, obj_mr%margin, obj_mr%kappa )
      n1 = domain%blk(b)%dim(1); n2 = domain%blk(b)%dim(2); n3 = domain%blk(b)%dim(3)
      associate ( lv => mrb(b)%lev, l2 => mrb(b)%lev2 )
      !$acc parallel loop collapse(3) gang vector present(lv, l2)
      do k = 1, n3
      do j = 1, n2
      do i = 1, n1
        lv(i,j,k) = max( lv(i,j,k), l2(i,j,k) )   ! the levels only rise
      end do
      end do
      end do
      end associate
    end do
    cnt(1) = 0.0_R8
    if ( L > 0 ) then
      nchg = mr_nest ( domain, L )
      do rnd = 1, 2*L + 2
        nchg = mr_seam ( domain )
        nchg = nchg + mr_nest ( domain, L )
        cnt(1) = real(nchg, R8); cnt(2) = 0.0_R8
        call mpi_allreduce_sum_r8_array ( cnt, 2 )
        if ( cnt(1) < 0.5_R8 ) exit
      end do
      if ( cnt(1) > 0.5_R8 ) error stop 'multirate redo: the level rules (R1 nesting, R2 connections) did not converge'
    end if
    hcap = mr_cap ( domain )
    obj_mr%L = L
    call system_clock ( cp_ )
    call mr_build_from_levels ( domain, L, cp_ )
    deeper = ( obj_mr%L > Lold )
    obj_mr%nsince = 0; obj_mr%nfull = obj_mr%nfull + 1
    obj_mr%Hcap = hcap; obj_mr%Hbuilt = min( H, hcap ); obj_mr%have_part = .true.
  end subroutine MR_Rebuild_For_H


  !> mr-restart: write the partition state of the current step (the levels of the local blocks and the scalars of the
  !> light-step test) to the open unformatted unit u, after the counters of the run (Wrap_Postprocess). The host levels
  !> are current: a full partition updates them and the light steps do not change them.
  subroutine MR_Write_State ( domain, u )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer, intent(in) :: u
    integer :: b, nloc
    write(u) obj_mr%drive .and. obj_mr%have_part
    if ( .not. ( obj_mr%drive .and. obj_mr%have_part ) ) return
    write(u) obj_mr%L, obj_mr%nsince, obj_mr%nlight, obj_mr%nfull
    write(u) obj_mr%H, obj_mr%Hreq, obj_mr%Hcap, obj_mr%Hbuilt, obj_mr%Hreqb, obj_mr%ratio_eff
    nloc = 0
    do b = 1, domain % nb
      if ( is_local_block(b) ) nloc = nloc + 1
    end do
    write(u) nloc
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      write(u) b, domain%blk(b)%dim
      write(u) mrb(b)%lev
    end do
  end subroutine MR_Write_State

  !> mr-restart: read what MR_Write_State wrote; the first MR_Partition applies it. Returns ok = f (and applies nothing)
  !> when the file has no partition or its blocks are not the local blocks of this run.
  subroutine MR_Read_State ( domain, u, ok )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer, intent(in)  :: u
    logical, intent(out) :: ok
    logical :: has
    integer :: b, bb, nloc, e, n(3)
    ok = .false.
    read(u) has
    if ( .not. has ) return
    read(u) rs_L, rs_nsince, rs_nlight, rs_nfull
    read(u) rs_H, rs_Hreq, rs_Hcap, rs_Hbuilt, rs_Hreqb, rs_ratio_eff
    read(u) nloc
    if ( allocated(rs_blk) ) deallocate ( rs_blk )
    allocate ( rs_blk(domain%nb) )
    do e = 1, nloc
      read(u) bb, n
      if ( bb < 1 .or. bb > domain%nb ) return
      if ( .not. is_local_block(bb) .or. any( n /= domain%blk(bb)%dim ) ) return
      allocate ( rs_blk(bb)%lev(n(1),n(2),n(3)) )
      read(u) rs_blk(bb)%lev
    end do
    do b = 1, domain % nb
      if ( is_local_block(b) .and. .not. allocated(rs_blk(b)%lev) ) return
    end do
    rs_pending = .true.; ok = .true.
  end subroutine MR_Read_State

  !> mr-restart: the levels and the scalars read by MR_Read_State, and the structures built from them, as the uninterrupted
  !> run holds them after the step written; the light/full decision of the next step then runs as in that run.
  subroutine mr_restore_state ( domain )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block, mpi_is_root
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer    :: b, L
    integer(8) :: cp_
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      mrb(b)%lev = rs_blk(b)%lev
      !$acc update device(mrb(b)%lev)
    end do
    deallocate ( rs_blk )
    obj_mr%nfull = rs_nfull
    call system_clock ( cp_ )
    L = rs_L
    call mr_build_from_levels ( domain, L, cp_ )
    obj_mr%L = rs_L
    obj_mr%nsince = rs_nsince; obj_mr%nlight = rs_nlight
    obj_mr%H = rs_H; obj_mr%Hreq = rs_Hreq; obj_mr%Hcap = rs_Hcap; obj_mr%Hbuilt = rs_Hbuilt; obj_mr%Hreqb = rs_Hreqb
    obj_mr%ratio_eff = rs_ratio_eff
    obj_mr%have_part = .true.
    rs_pending = .false.
    if ( mpi_is_root ) write(*,'(A,I0,A,I0,A,ES12.5)') ' [MR] restart: partition restored (L ', rs_L, ', since full ', &
      rs_nsince, ', H built ', rs_Hbuilt
  end subroutine mr_restore_state


  !> mr-xch: the inter-rank halo entries each level's stages read -- the entry's two ghost cells inside a box of the level
  !> grown by gc (the rule of mr_build_ghost_list: the cells the thermo cache, the face stencils, the boundary records and
  !> the SST source of the level read) -- as a bit mask per receive entry, to Mod_GhostExchange, which ships it to the
  !> donors and builds the per-level sub-lists on both sides. The levels are frozen between full partitions.
  subroutine mr_build_xch_masks ( domain, L )
    use MOSE_Types_m
    use MOSE_Global_m, only: gc, guide
    use MOSE_Mod_MPI, only: mpi_size_, mpi_is_root
    use MOSE_Mod_GhostExchange, only: ghost_sched, ghost_sched_set_level_masks, xch_level_count, xch_lazy_setup
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer, intent(in) :: L
    integer, allocatable :: rmask(:), cmask(:)
    integer :: nr, e, b, l_, ib, m, c1(3), d(3), g1(3), g2(3), lo(3), hi(3), ns_, nr_, nsg_, nrg_, p, c
    if ( mpi_size_ <= 1 ) return
    if ( .not. ghost_sched%built ) return
    call xch_lazy_setup ( domain )   ! the ghost-owner table of round 2 (and the pg-edge lists) if the flat field exists
    nr = ghost_sched%n_recv
    allocate ( rmask(max(nr,1)) ); rmask = 0
    do e = 1, nr
      b  = ghost_sched%recv_b(e)
      c1 = (/ ghost_sched%recv_id(e), ghost_sched%recv_jd(e), ghost_sched%recv_kd(e) /)
      d  = (/ ghost_sched%recv_di(e), ghost_sched%recv_dj(e), ghost_sched%recv_dk(e) /)
      g1 = c1 - d; g2 = c1 - 2*d
      m = 0
      do l_ = 0, L
        do ib = 1, mrb(b)%nbox(l_)
          lo = (/ mrb(b)%box(1,ib,l_), mrb(b)%box(3,ib,l_), mrb(b)%box(5,ib,l_) /) - gc
          hi = (/ mrb(b)%box(2,ib,l_), mrb(b)%box(4,ib,l_), mrb(b)%box(6,ib,l_) /) + gc
          if ( all( g1 >= lo .and. g1 <= hi ) .or. all( g2 >= lo .and. g2 <= hi ) ) then
            m = ibset( m, l_ ); exit
          end if
        end do
      end do
      ! never unread: at least the level of its own rows
      if ( m == 0 ) m = ibset( 0, max( 0, min( L, mrb(b)%lev(c1(1),c1(2),c1(3)) ) ) )
      rmask(e) = m
    end do
    ! the same-rank connection records: the levels whose pass-2 gather (Ghost_Connection_pairs_B, mr-ghost-near rule: a
    ! ghost of the record in the near set) reads them -- their donor-side edge reads join the sub-lists (round 2)
    allocate ( cmask(max(domain%n_conn,1)) ); cmask = 0
    do p = 1, domain % n_conn_pair
      b = domain % conn_pair_bm(p)
      do c = domain % conn_pair_lo(p), domain % conn_pair_hi(p)
        c1 = (/ domain%conn_im(c), domain%conn_jm(c), domain%conn_km(c) /)
        d  = (/ guide(domain%conn_fm(c),1), guide(domain%conn_fm(c),2), guide(domain%conn_fm(c),3) /)
        g1 = c1 - d; g2 = c1 - 2*d
        m = 0
        do l_ = 0, L
          do ib = 1, mrb(b)%nbox(l_)
            lo = (/ mrb(b)%box(1,ib,l_), mrb(b)%box(3,ib,l_), mrb(b)%box(5,ib,l_) /) - gc
            hi = (/ mrb(b)%box(2,ib,l_), mrb(b)%box(4,ib,l_), mrb(b)%box(6,ib,l_) /) + gc
            if ( all( g1 >= lo .and. g1 <= hi ) .or. all( g2 >= lo .and. g2 <= hi ) ) then
              m = ibset( m, l_ ); exit
            end if
          end do
        end do
        cmask(c) = m
      end do
    end do
    call ghost_sched_set_level_masks ( domain, rmask, cmask, L + 1 )
    deallocate ( cmask )
    if ( mpi_is_root .and. obj_mr%nfull == 1 ) then
      write(*,'(A)',advance='no') ' [MR]   inter-rank halo entries read per level (rank 0, send/recv):'
      do l_ = 0, L
        call xch_level_count ( l_, ns_, nr_, nsg_, nrg_ )
        write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)',advance='no') '  L', l_, ' ', ns_, '/', nr_, ' (Pg ', nsg_, '/', nrg_, ')'
      end do
      write(*,'(A)') ''
    end if
    deallocate ( rmask )
  end subroutine mr_build_xch_masks


  !> Set_Global_dt under multirate: each cell's own sub-step dtlocal = DT / 2**lev (read by Newstate; it is also the level
  !> predicate of the stage kernels, Lib_Multirate_Base).
  subroutine MR_Set_Local_dt ( domain )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: b
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      call mr_set_dt_blk ( domain%blk(b)%dtlocal, mrb(b)%lev, domain%blk(b)%dim, domain%dtglobal, b )
    end do
  end subroutine MR_Set_Local_dt


  !------------------------------------------------------------------------------------------------------------------------
  !> Once: the scheme the driver supports, the per-block arrays (device-resident) and the wall rows.
  subroutine MR_Setup ( domain, config )
    use MOSE_Types_m
    use MOSE_Global_m, only: ROTATING_FRAME, soot, RANS
    use MOSE_Lib_RANS, only: eddy_dev
    use MOSE_Mod_MPI, only: is_local_block, mpi_is_root, mpi_allreduce_sum_r8_array
    use MOSE_Lib_RK, only: obj_rk
    use MOSE_Lib_Strang, only: N_strang
    use MOSE_Lib_IRS, only: obj_irs
    use MOSE_Lib_Preconditioning, only: obj_prec
    use MOSE_Lib_Omega_Relax, only: obj_omega_relax
    use MOSE_Lib_Limiters, only: lim_dev
    use MOSE_Global_m, only: gc, VISC_ => VISC
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    type(config_type),      intent(in)    :: config
    integer :: b, n(3), e, c, r, i, j, k, f, nw
    real(R8) :: cnt(2), theta_l
    if ( obj_mr%drive ) then
      if ( obj_rk%n_RK /= 2 ) error stop 'multirate: time-scheme = RK2 only (the ladder is built on the Heun stages)'
      if ( trim(config%intvar_name) /= 'conservative' ) error stop 'multirate: integration-variables = cons only'
      if ( .not. config%time_accurate .or. config%HYDRA_time_accurate ) error stop 'multirate: time-accurate runs only'
      if ( N_strang /= 1 ) error stop 'multirate: frozen chemistry only (no Strang splitting in this version)'
      if ( obj_irs%IRS ) error stop 'multirate: no implicit residual smoothing'
      if ( obj_prec%enabled ) error stop 'multirate: no preconditioning'
      if ( obj_omega_relax%active ) error stop 'multirate: no omega-source-relax'
      if ( ROTATING_FRAME .or. soot ) error stop 'multirate: no rotating frame and no soot (host sources)'
      if ( RANS .and. eddy_dev /= 1 ) error stop 'multirate: turbulence SST or laminar only (other sources sweep the block)'
      if ( domain%n_chim > 0 ) error stop 'multirate: no chimera faces (666) in this version'
      if ( config%MGL /= 1 ) error stop 'multirate: one multigrid level only'
    end if
    ! mr-margin: kappa = max( cfl (1 + theta), 1.1 ), theta of the limiter (open-hydra MOSE feature/multirate-integrator,
    ! Mod_Multirate.f90 1307-1322 there): 1/2 minmod, 1 van Leer, MC and superbee, sqrt(2)/2 van Albada, 0 first order
    select case ( lim_dev )
      case (1);       theta_l = 0.5_R8
      case (2, 4, 5); theta_l = 1.0_R8
      case (3);       theta_l = sqrt(2.0_R8) / 2.0_R8
      case default;   theta_l = 0.0_R8
    end select
    obj_mr%kappa = max( config%cfl * ( 1.0_R8 + theta_l ), 1.1_R8 )
    allocate ( mrb(domain%nb) )
    nw = 0
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      n = domain%blk(b)%dim
      allocate ( mrb(b)%lev(n(1),n(2),n(3)), mrb(b)%lev2(n(1),n(2),n(3)), mrb(b)%cls(n(1),n(2),n(3)), &
                 mrb(b)%dtcap(n(1),n(2),n(3)) )
      mrb(b)%lev = 0; mrb(b)%lev2 = 0; mrb(b)%cls = 0; mrb(b)%dtcap = 0.0_R8
      ! the wall rows: mr-wall-rows cells from every viscous wall record, along the inward normal of its face
      if ( obj_mr%sigma_wall > 0.0_R8 .and. domain%bc_blk_hi(b) >= domain%bc_blk_lo(b) ) then
        do e = domain%bc_blk_lo(b), domain%bc_blk_hi(b)
          if ( .not. any( WALL_TYPES == domain%bc_typ(e) ) ) cycle
          i = domain%bc_bi(e); j = domain%bc_bj(e); k = domain%bc_bk(e); f = domain%bc_bf(e)
          do r = 0, obj_mr%wall_rows - 1
            select case (f)
              case (1); c = i + r; if (c > n(1)) exit; mrb(b)%cls(c,j,k) = 1
              case (2); c = i - r; if (c < 1)    exit; mrb(b)%cls(c,j,k) = 1
              case (3); c = j + r; if (c > n(2)) exit; mrb(b)%cls(i,c,k) = 1
              case (4); c = j - r; if (c < 1)    exit; mrb(b)%cls(i,c,k) = 1
              case (5); c = k + r; if (c > n(3)) exit; mrb(b)%cls(i,j,c) = 1
              case (6); c = k - r; if (c < 1)    exit; mrb(b)%cls(i,j,c) = 1
            end select
          end do
        end do
        nw = nw + count( mrb(b)%cls == 1 )
      end if
      !$acc enter data copyin(mrb(b)%lev, mrb(b)%lev2, mrb(b)%cls, mrb(b)%dtcap)
      if ( obj_mr%nonlocal > 0 ) then   ! mr-nonlocal: the certificate's arrays
        allocate ( mrb(b)%dtc(n(1),n(2),n(3)), mrb(b)%dmin(n(1),n(2),n(3)), &
                   mrb(b)%crho(1-gc:n(1)+gc,1-gc:n(2)+gc,1-gc:n(3)+gc), mrb(b)%crg(1-gc:n(1)+gc,1-gc:n(2)+gc,1-gc:n(3)+gc), &
                   mrb(b)%ct(1-gc:n(1)+gc,1-gc:n(2)+gc,1-gc:n(3)+gc), mrb(b)%cmil(1-gc:n(1)+gc,1-gc:n(2)+gc,1-gc:n(3)+gc) )
        mrb(b)%dtc = 0.0_R8; mrb(b)%dmin = 0.0_R8
        mrb(b)%crho = 1.0_R8; mrb(b)%crg = 1.0_R8; mrb(b)%ct = 300.0_R8; mrb(b)%cmil = 0.0_R8
        !$acc enter data copyin(mrb(b)%dtc, mrb(b)%dmin, mrb(b)%crho, mrb(b)%crg, mrb(b)%ct, mrb(b)%cmil)
      end if
    end do
    if ( obj_mr%nonlocal == 2 ) then   ! redo: the linear stability of the accepted sub-steps (the branch, BT 1151-1208)
      if ( config%cfl > 1.0_R8 ) error stop 'mr-nonlocal = redo: cfl <= 1'
      if ( VISC_ .and. config%vnn > 0.5_R8 ) error stop 'mr-nonlocal = redo: vnn <= 0.5'
    end if
    cnt(1) = real(nw, R8); cnt(2) = real(mr_ncell(domain), R8)
    call mpi_allreduce_sum_r8_array ( cnt, 2 )
    obj_mr%nwallc = nint(cnt(1))
    if ( mpi_is_root ) then
      write(*,'(A,L1,A,F6.1,A,F5.2,A,I0)') ' [MR] multirate: integrates ', obj_mr%drive, '  mr-ratio ', obj_mr%ratio, &
        '  mr-safety ', obj_mr%sigma, '  mr-lmax ', obj_mr%lmax
      if ( obj_mr%sigma_wall > 0.0_R8 ) write(*,'(A,F5.2,A,I0,A,I0,A,F7.3,A)') ' [MR]   wall margin ', obj_mr%sigma_wall, &
        ' on ', obj_mr%wall_rows, ' rows from the viscous walls: ', obj_mr%nwallc, ' cells (', &
        100.0_R8*cnt(1)/max(cnt(2),1.0_R8), ' %)'
      write(*,'(A,I0,A,F5.2)') ' [MR]   levels rebuilt at least every ', obj_mr%repart, ' macro steps, or when H_cap/H < ', &
        obj_mr%repart_hcap
      write(*,'(A,L1,A,F6.3,A,L1,A,L1,A,L1,A,I0)') ' [MR]   mr-margin ', obj_mr%margin, ' (kappa ', obj_mr%kappa, &
        ')  mr-depth auto ', obj_mr%depth_auto, '  mr-cons-halo ', obj_mr%cons_halo, '  mr-land ', obj_mr%land, &
        '  mr-nonlocal ', obj_mr%nonlocal
    end if
  end subroutine MR_Setup


  integer function mr_ncell ( domain ) result ( nc )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer :: b
    nc = 0
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      nc = nc + domain%blk(b)%dim(1) * domain%blk(b)%dim(2) * domain%blk(b)%dim(3)
    end do
  end function mr_ncell


  !> dtcap <- dtlocal on every local block (the per-cell limits Compute_dt just wrote).
  subroutine mr_store_dtcap ( domain )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: b
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      call mr_copy_dt_blk ( domain%blk(b)%dtlocal, mrb(b)%dtcap, domain%blk(b)%dim )
    end do
  end subroutine mr_store_dtcap

  !> mr-octaves: the admissible steps dtcap of the cells (the values the levels are built from) in octaves
  !> k = floor(log2(dtcap/dt_min)), the share of the wall rows (margin class 1) in each, and for a macro step
  !> H = 2^m dt_min the ideal work per unit time, sum over the cells of 2^lev / H with lev = max(0, ceil(log2(s H/dtcap)))
  !> and no level rule (the rules only raise levels), relative to the work of the requested step H = mr-ratio dt_min.
  !> Host copy of dtcap: a diagnostic of the first full partition and of every 10th.
  subroutine mr_octave_census ( domain )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block, mpi_is_root, mpi_allreduce_sum_r8_array
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, parameter :: NK = 24, NM = 14
    real(R8) :: oct(0:NK), octw(0:NK), work(0:NM), x, s_, sw_, dtmin, hm, ratio_req
    real(R8) :: buf(2*(NK+1) + (NM+1) + 1)
    integer  :: b, i, j, k, n1, n2, n3, kk, m, lev
    dtmin = obj_mr%dtmin
    if ( dtmin <= 0.0_R8 ) return
    s_ = obj_mr%sigma; sw_ = max( obj_mr%sigma_wall, obj_mr%sigma )
    oct = 0.0_R8; octw = 0.0_R8; work = 0.0_R8
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$acc update host(mrb(b)%dtcap, mrb(b)%cls)
      n1 = domain%blk(b)%dim(1); n2 = domain%blk(b)%dim(2); n3 = domain%blk(b)%dim(3)
      do k = 1, n3; do j = 1, n2; do i = 1, n1
        x = mrb(b)%dtcap(i,j,k) / dtmin
        if ( .not. ( x > 0.0_R8 ) ) cycle
        kk = min( max( int( floor( log(x) / log(2.0_R8) + 1.0e-12_R8 ) ), 0 ), NK )
        oct(kk) = oct(kk) + 1.0_R8
        if ( mrb(b)%cls(i,j,k) == 1 ) octw(kk) = octw(kk) + 1.0_R8
        do m = 0, NM   ! H = 2^m dt_min
          hm = 2.0_R8**m
          if ( mrb(b)%cls(i,j,k) == 1 ) then
            lev = max( 0, ceiling( log( sw_ * hm / x ) / log(2.0_R8) - 1.0e-12_R8 ) )
          else
            lev = max( 0, ceiling( log( s_ * hm / x ) / log(2.0_R8) - 1.0e-12_R8 ) )
          end if
          work(m) = work(m) + 2.0_R8**lev / hm
        end do
      end do; end do; end do
    end do
    buf(1:NK+1) = oct; buf(NK+2:2*NK+2) = octw; buf(2*NK+3:2*NK+3+NM) = work; buf(2*NK+4+NM) = 0.0_R8
    call mpi_allreduce_sum_r8_array ( buf, size(buf) )
    oct = buf(1:NK+1); octw = buf(NK+2:2*NK+2); work = buf(2*NK+3:2*NK+3+NM)
    if ( .not. mpi_is_root ) return
    write(*,'(A,I0,A,ES11.4,A,F6.2,A,F6.2)') ' [MR-OCT] full partition ', obj_mr%nfull, '  dt_min ', dtmin, &
      '  safety ', s_, '  wall safety ', sw_
    write(*,'(A)') ' [MR-OCT] cells per octave of dtcap/dt_min (k = floor log2), of which wall rows:'
    do kk = 0, NK
      if ( oct(kk) > 0.0_R8 ) write(*,'(A,I3,A,F14.0,A,F8.4,A,F14.0)') ' [MR-OCT]   k ', kk, '  cells ', oct(kk), &
        '  share ', oct(kk) / sum(oct), '  wall rows ', octw(kk)
    end do
    ratio_req = obj_mr%ratio_eff
    write(*,'(A,F6.1,A)') ' [MR-OCT] ideal work (no level rules) of H = 2^m dt_min over the work of H = ', ratio_req, ' dt_min:'
    m = nint( log(ratio_req) / log(2.0_R8) )
    do kk = 0, NM
      write(*,'(A,I3,A,F8.1,A,F9.4)') ' [MR-OCT]   m ', kk, '  H/dt_min ', 2.0_R8**kk, '  work ', work(kk) / work(max(min(m,NM),0))
    end do
  end subroutine mr_octave_census

  subroutine mr_copy_dt_blk ( dtl, dtcap, n )
    implicit none
    integer,  intent(in)    :: n(3)
    real(R8), intent(in)    :: dtl(n(1),n(2),n(3))
    real(R8), intent(inout) :: dtcap(n(1),n(2),n(3))
    integer :: i, j, k, n1, n2, n3
    n1 = n(1); n2 = n(2); n3 = n(3)
    !$acc parallel loop collapse(3) gang vector copy(dtl) present(dtcap)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      dtcap(i,j,k) = dtl(i,j,k)
    end do
    end do
    end do
  end subroutine mr_copy_dt_blk


  !> The cap of the current levels over all ranks: min over the cells of dtcap 2**lev / s.
  real(R8) function mr_cap ( domain ) result ( hcap )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block, mpi_allreduce_min_r8
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer  :: b
    real(R8) :: h
    h = huge(1.0_R8)
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      call mr_cap_blk ( mrb(b)%dtcap, mrb(b)%cls, mrb(b)%lev, domain%blk(b)%dim, obj_mr%sigma, obj_mr%sigma_wall, h, &
                        obj_mr%margin, obj_mr%kappa )
    end do
    call mpi_allreduce_min_r8 ( h, hcap )
  end function mr_cap

  !> mr-depth = auto: the largest admissible step of the cells, over the local blocks and the ranks (the min of -x).
  real(R8) function mr_dtcap_max ( domain ) result ( dmax )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block, mpi_allreduce_min_r8
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer  :: b, i, j, k, n1, n2, n3
    real(R8) :: h, hl, mh
    h = -huge(1.0_R8)
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      n1 = domain%blk(b)%dim(1); n2 = domain%blk(b)%dim(2); n3 = domain%blk(b)%dim(3)
      hl = -huge(1.0_R8)
      associate ( dc => mrb(b)%dtcap )
      !$acc parallel loop collapse(3) gang vector present(dc) reduction(max:hl)
      do k = 1, n3
      do j = 1, n2
      do i = 1, n1
        hl = max( hl, dc(i,j,k) )
      end do
      end do
      end do
      end associate
      h = max( h, hl )
    end do
    call mpi_allreduce_min_r8 ( -h, mh )
    dmax = -mh
  end function mr_dtcap_max

  subroutine mr_cap_blk ( dtcap, cls, lev, n, sig, sigw, hcap, marg, kap )
    implicit none
    integer,  intent(in)    :: n(3), cls(n(1),n(2),n(3)), lev(n(1),n(2),n(3))
    real(R8), intent(in)    :: dtcap(n(1),n(2),n(3)), sig, sigw, kap
    logical,  intent(in)    :: marg
    real(R8), intent(inout) :: hcap
    integer  :: i, j, k, n1, n2, n3
    real(R8) :: h_, s_, sig_, sigw_, kap_, dnb
    logical  :: marg_
    n1 = n(1); n2 = n(2); n3 = n(3); h_ = huge(1.0_R8); sig_ = sig; sigw_ = max(sigw, sig)
    marg_ = marg; kap_ = kap
    !$acc parallel loop collapse(3) gang vector present(dtcap, cls, lev) reduction(min:h_) private(s_, dnb)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      s_ = sig_
      if ( cls(i,j,k) == 1 ) s_ = sigw_
      h_ = min( h_, dtcap(i,j,k) * ( 2.0_R8 ** lev(i,j,k) ) / s_ )
      if ( marg_ ) then   ! mr-margin: min over the face neighbours of dtcap 2**lev_c / kappa
        dnb = huge(1.0_R8)
        if ( i > 1  ) dnb = min( dnb, dtcap(i-1,j,k) )
        if ( i < n1 ) dnb = min( dnb, dtcap(i+1,j,k) )
        if ( j > 1  ) dnb = min( dnb, dtcap(i,j-1,k) )
        if ( j < n2 ) dnb = min( dnb, dtcap(i,j+1,k) )
        if ( k > 1  ) dnb = min( dnb, dtcap(i,j,k-1) )
        if ( k < n3 ) dnb = min( dnb, dtcap(i,j,k+1) )
        h_ = min( h_, dnb * ( 2.0_R8 ** lev(i,j,k) ) / kap_ )
      end if
    end do
    end do
    end do
    hcap = min( hcap, h_ )
  end subroutine mr_cap_blk


  !> The raw levels of one block for the macro step H; adds sum(2**lev) to sump.
  subroutine mr_levels_blk ( dtcap, cls, lev, n, H, sig, sigw, L, sump, marg, kap )
    implicit none
    integer,  intent(in)    :: n(3), L, cls(n(1),n(2),n(3))
    real(R8), intent(in)    :: dtcap(n(1),n(2),n(3)), H, sig, sigw, kap
    logical,  intent(in)    :: marg
    integer,  intent(inout) :: lev(n(1),n(2),n(3))
    real(R8), intent(inout) :: sump
    integer  :: i, j, k, n1, n2, n3, l_, lcap
    real(R8) :: x, s_, H_, sig_, sigw_, sp, kap_, dnb
    logical  :: marg_
    n1 = n(1); n2 = n(2); n3 = n(3); lcap = L; H_ = H; sig_ = sig; sigw_ = max(sigw, sig); sp = 0.0_R8
    marg_ = marg; kap_ = kap
    !$acc parallel loop collapse(3) gang vector present(dtcap, cls, lev) reduction(+:sp) private(x, s_, l_, dnb)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      s_ = sig_
      if ( cls(i,j,k) == 1 ) s_ = sigw_
      x  = s_ * H_ / dtcap(i,j,k)
      if ( marg_ ) then   ! mr-margin: x = max( s H / dtcap_c, kappa H / min over the face neighbours of dtcap )
        dnb = huge(1.0_R8)
        if ( i > 1  ) dnb = min( dnb, dtcap(i-1,j,k) )
        if ( i < n1 ) dnb = min( dnb, dtcap(i+1,j,k) )
        if ( j > 1  ) dnb = min( dnb, dtcap(i,j-1,k) )
        if ( j < n2 ) dnb = min( dnb, dtcap(i,j+1,k) )
        if ( k > 1  ) dnb = min( dnb, dtcap(i,j,k-1) )
        if ( k < n3 ) dnb = min( dnb, dtcap(i,j,k+1) )
        x = max( x, kap_ * H_ / dnb )
      end if
      l_ = 0
      if ( x > 1.0_R8 ) l_ = ceiling( log(x) / log(2.0_R8) - 1.0e-12_R8 )
      l_ = min( max(l_, 0), lcap )
      lev(i,j,k) = l_
      sp = sp + 2.0_R8 ** l_
    end do
    end do
    end do
    sump = sump + sp
  end subroutine mr_levels_blk


  subroutine mr_census_blk ( lev, n, hist, sump )
    implicit none
    integer,  intent(in)    :: n(3), lev(n(1),n(2),n(3))
    real(R8), intent(inout) :: hist(0:MR_LMAX_HARD), sump
    integer  :: i, j, k, n1, n2, n3, l_
    real(R8) :: s_, h0, h1, h2, h3, h4, h5, h6, h7, h8
    n1 = n(1); n2 = n(2); n3 = n(3); s_ = 0.0_R8
    h0 = 0; h1 = 0; h2 = 0; h3 = 0; h4 = 0; h5 = 0; h6 = 0; h7 = 0; h8 = 0
    !$acc parallel loop collapse(3) gang vector present(lev) &
    !$acc&  reduction(+:s_, h0, h1, h2, h3, h4, h5, h6, h7, h8) private(l_)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      l_ = lev(i,j,k)
      s_ = s_ + 2.0_R8 ** l_
      select case (l_)
        case (0); h0 = h0 + 1
        case (1); h1 = h1 + 1
        case (2); h2 = h2 + 1
        case (3); h3 = h3 + 1
        case (4); h4 = h4 + 1
        case (5); h5 = h5 + 1
        case (6); h6 = h6 + 1
        case (7); h7 = h7 + 1
        case default; h8 = h8 + 1
      end select
    end do
    end do
    end do
    hist = hist + (/ h0, h1, h2, h3, h4, h5, h6, h7, h8 /)
    sump = sump + s_
  end subroutine mr_census_blk


  !> R1, L passes on every local block: lev <- max(lev, max over the six face neighbours of lev - 1). Returns the cells raised.
  integer function mr_nest ( domain, L ) result ( nchg )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: L
    integer :: b, p, nc
    nchg = 0
    if ( mr_build_fast_on() ) then
      ! mr-build-fast: the passes of every block on the list queue (in order: block b's pass p reads its pass p-1), the
      ! changed cells summed on the device, one read back -- not one synchronous reduction per block and pass
      ld_q = domain % nb + 65
      if ( .not. allocated(ld_nc) ) then
        allocate ( ld_nc(1) )
        !$acc enter data create(ld_nc)
      end if
      call ld_set_i ( ld_nc, 1, 0 )
      do p = 1, L
        do b = 1, domain % nb
          if (.not. is_local_block(b)) cycle
          call mr_nest_blk_q ( mrb(b)%lev, mrb(b)%lev2, domain%blk(b)%dim, ld_nc )
        end do
      end do
      !$acc update host(ld_nc) async(ld_q)
      !$acc wait(ld_q)
      nchg = ld_nc(1)
      return
    end if
    do p = 1, L
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        call mr_nest_blk ( mrb(b)%lev, mrb(b)%lev2, domain%blk(b)%dim, nc )
        nchg = nchg + nc
      end do
    end do
  end function mr_nest

  subroutine mr_nest_blk ( lev, lev2, n, nc )
    implicit none
    integer, intent(in)    :: n(3)
    integer, intent(inout) :: lev(n(1),n(2),n(3)), lev2(n(1),n(2),n(3))
    integer, intent(out)   :: nc
    integer :: i, j, k, n1, n2, n3, m, c
    n1 = n(1); n2 = n(2); n3 = n(3); c = 0
    !$acc parallel loop collapse(3) gang vector present(lev, lev2) reduction(+:c) private(m)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      m = lev(i,j,k)
      if (i > 1)  m = max(m, lev(i-1,j,k) - 1)
      if (i < n1) m = max(m, lev(i+1,j,k) - 1)
      if (j > 1)  m = max(m, lev(i,j-1,k) - 1)
      if (j < n2) m = max(m, lev(i,j+1,k) - 1)
      if (k > 1)  m = max(m, lev(i,j,k-1) - 1)
      if (k < n3) m = max(m, lev(i,j,k+1) - 1)
      if (m /= lev(i,j,k)) c = c + 1
      lev2(i,j,k) = m
    end do
    end do
    end do
    !$acc parallel loop collapse(3) gang vector present(lev, lev2)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      lev(i,j,k) = lev2(i,j,k)
    end do
    end do
    end do
    nc = c
  end subroutine mr_nest_blk

  !> mr-build-fast: mr_nest_blk on the list queue; the changed cells of each chunk of LD_CH cells added to nc(1) (one atomic
  !> per chunk; integers, so the order of the additions does not matter). The same levels as mr_nest_blk.
  subroutine mr_nest_blk_q ( lev, lev2, n, nc )
    implicit none
    integer, intent(in)    :: n(3)
    integer, intent(inout) :: lev(n(1),n(2),n(3)), lev2(n(1),n(2),n(3)), nc(1)
    integer :: i, j, k, n1, n2, n3, m, c, ncell, nch, ch, e, r
    n1 = n(1); n2 = n(2); n3 = n(3); ncell = n1*n2*n3; nch = (ncell + LD_CH - 1) / LD_CH
    !$acc parallel loop gang async(ld_q) present(lev, lev2, nc) private(c)
    do ch = 1, nch
      c = 0
      !$acc loop vector reduction(+:c) private(i, j, k, r, m)
      do e = (ch-1)*LD_CH + 1, min(ch*LD_CH, ncell)
        r = e - 1; i = mod(r, n1) + 1; r = r / n1; j = mod(r, n2) + 1; k = r / n2 + 1
        m = lev(i,j,k)
        if (i > 1)  m = max(m, lev(i-1,j,k) - 1)
        if (i < n1) m = max(m, lev(i+1,j,k) - 1)
        if (j > 1)  m = max(m, lev(i,j-1,k) - 1)
        if (j < n2) m = max(m, lev(i,j+1,k) - 1)
        if (k > 1)  m = max(m, lev(i,j,k-1) - 1)
        if (k < n3) m = max(m, lev(i,j,k+1) - 1)
        if (m /= lev(i,j,k)) c = c + 1
        lev2(i,j,k) = m
      end do
      if ( c > 0 ) then
        !$acc atomic update
        nc(1) = nc(1) + c
      end if
    end do
    !$acc parallel loop collapse(3) gang vector async(ld_q) present(lev, lev2)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      lev(i,j,k) = lev2(i,j,k)
    end do
    end do
    end do
  end subroutine mr_nest_blk_q


  !> Interior faces of the local blocks whose two cells differ by more than one level (0 after the rules).
  integer function mr_face_jumps ( domain ) result ( nv )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block, mpi_allreduce_sum_r8_array
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer  :: b, nb_
    real(R8) :: cnt(2)
    nv = 0
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      call mr_jumps_blk ( mrb(b)%lev, domain%blk(b)%dim, nb_ )
      nv = nv + nb_
    end do
    cnt(1) = real(nv, R8); cnt(2) = 0.0_R8
    call mpi_allreduce_sum_r8_array ( cnt, 2 )
    nv = nint(cnt(1))
  end function mr_face_jumps

  subroutine mr_jumps_blk ( lev, n, nv )
    implicit none
    integer, intent(in)  :: n(3), lev(n(1),n(2),n(3))
    integer, intent(out) :: nv
    integer :: i, j, k, n1, n2, n3, c
    n1 = n(1); n2 = n(2); n3 = n(3); c = 0
    !$acc parallel loop collapse(3) gang vector present(lev) reduction(+:c)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      if (i < n1) then; if (abs(lev(i,j,k) - lev(i+1,j,k)) > 1) c = c + 1; end if
      if (j < n2) then; if (abs(lev(i,j,k) - lev(i,j+1,k)) > 1) c = c + 1; end if
      if (k < n3) then; if (abs(lev(i,j,k) - lev(i,j,k+1)) > 1) c = c + 1; end if
    end do
    end do
    end do
    nv = c
  end subroutine mr_jumps_blk


  !> R2 on the host: the first two rows on both sides of every connection face (local pairs, then the rank seams through
  !> the halo schedule) take the largest of their four levels. Returns the cells raised on this rank.
  integer function mr_seam ( domain ) result ( nchg )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: b, p, c, bm, bs, m
    integer :: a1(3), a2(3), b1(3), b2(3)
    nchg = 0
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$acc update host(mrb(b)%lev)
    end do
    do p = 1, domain % n_conn_pair
      bm = domain % conn_pair_bm(p); bs = domain % conn_pair_bs(p)
      do c = domain % conn_pair_lo(p), domain % conn_pair_hi(p)
        a1 = (/ domain%conn_im(c), domain%conn_jm(c), domain%conn_km(c) /)
        b1 = (/ domain%conn_is(c), domain%conn_js(c), domain%conn_ks(c) /)
        a2 = mr_row2 ( a1, domain%conn_fm(c), domain%blk(bm)%dim )
        b2 = mr_row2 ( b1, domain%conn_fs(c), domain%blk(bs)%dim )
        m = max( mrb(bm)%lev(a1(1),a1(2),a1(3)), mrb(bm)%lev(a2(1),a2(2),a2(3)), &
                 mrb(bs)%lev(b1(1),b1(2),b1(3)), mrb(bs)%lev(b2(1),b2(2),b2(3)) )
        call mr_raise ( mrb(bm)%lev(a1(1),a1(2),a1(3)), m, nchg )
        call mr_raise ( mrb(bm)%lev(a2(1),a2(2),a2(3)), m, nchg )
        call mr_raise ( mrb(bs)%lev(b1(1),b1(2),b1(3)), m, nchg )
        call mr_raise ( mrb(bs)%lev(b2(1),b2(2),b2(3)), m, nchg )
      end do
    end do
    call mr_seam_ranks ( domain, nchg )
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$acc update device(mrb(b)%lev)
    end do
  end function mr_seam

  subroutine mr_raise ( lv, m, nchg )
    implicit none
    integer, intent(inout) :: lv, nchg
    integer, intent(in)    :: m
    if ( lv < m ) then
      lv = m; nchg = nchg + 1
    end if
  end subroutine mr_raise

  !> The second row of a boundary cell: one cell inward along the normal of its face (the cell itself in a one-cell direction).
  function mr_row2 ( c1, f, n ) result ( c2 )
    implicit none
    integer, intent(in) :: c1(3), f, n(3)
    integer :: c2(3), d
    c2 = c1
    d = (f + 1) / 2
    if ( mod(f, 2) == 1 ) then
      c2(d) = min( c1(d) + 1, n(d) )
    else
      c2(d) = max( c1(d) - 1, 1 )
    end if
  end function mr_row2

  !> R2 across the ranks: every send entry of the halo schedule carries the larger level of its donor's two rows, every
  !> receive entry raises its own two rows to it (the schedule pairs the entries of the two ranks in the same order).
  subroutine mr_seam_ranks ( domain, nchg )
    use MOSE_Types_m
    use MOSE_Mod_GhostExchange, only: ghost_sched
#ifdef USE_MPI
    use mpi
    use MOSE_Mod_MPI, only: mpi_rank_
#endif
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(inout) :: nchg
#ifdef USE_MPI
    integer, allocatable :: sbuf(:), rbuf(:), sreq(:), rreq(:)
    integer :: e, r, ierr, tag, b, ns, nr, nsr, nrr, m
    integer :: c1(3), c2(3)
    if ( .not. ghost_sched%built ) return
    ns = ghost_sched%n_send; nr = ghost_sched%n_recv
    nsr = ghost_sched%n_send_ranks; nrr = ghost_sched%n_recv_ranks
    if ( ns == 0 .and. nr == 0 ) return
    allocate ( sbuf(max(ns,1)), rbuf(max(nr,1)), sreq(max(nsr,1)), rreq(max(nrr,1)) )
    do r = 1, nrr
      tag = ghost_sched%recv_groups(r)%rank + 41000
      call MPI_IRECV ( rbuf(ghost_sched%recv_groups(r)%offset), ghost_sched%recv_groups(r)%count, MPI_INTEGER, &
                       ghost_sched%recv_groups(r)%rank, tag, MPI_COMM_WORLD, rreq(r), ierr )
    end do
    do e = 1, ns
      b  = ghost_sched%send_list(e)%local_block
      c1 = (/ ghost_sched%send_list(e)%src_i, ghost_sched%send_list(e)%src_j, ghost_sched%send_list(e)%src_k /)
      c2 = mr_row2 ( c1, ghost_sched%send_list(e)%src_f, domain%blk(b)%dim )
      sbuf(e) = max( mrb(b)%lev(c1(1),c1(2),c1(3)), mrb(b)%lev(c2(1),c2(2),c2(3)) )
    end do
    do r = 1, nsr
      tag = mpi_rank_ + 41000
      call MPI_ISEND ( sbuf(ghost_sched%send_groups(r)%offset), ghost_sched%send_groups(r)%count, MPI_INTEGER, &
                       ghost_sched%send_groups(r)%rank, tag, MPI_COMM_WORLD, sreq(r), ierr )
    end do
    if ( nrr > 0 ) call MPI_WAITALL ( nrr, rreq, MPI_STATUSES_IGNORE, ierr )
    do e = 1, nr
      b  = ghost_sched%recv_list(e)%local_block
      c1 = (/ ghost_sched%recv_list(e)%dst_i, ghost_sched%recv_list(e)%dst_j, ghost_sched%recv_list(e)%dst_k /)
      c2 = mr_row2 ( c1, ghost_sched%recv_list(e)%dst_f, domain%blk(b)%dim )
      m  = max( rbuf(e), mrb(b)%lev(c1(1),c1(2),c1(3)), mrb(b)%lev(c2(1),c2(2),c2(3)) )
      call mr_raise ( mrb(b)%lev(c1(1),c1(2),c1(3)), m, nchg )
      call mr_raise ( mrb(b)%lev(c2(1),c2(2),c2(3)), m, nchg )
    end do
    if ( nsr > 0 ) call MPI_WAITALL ( nsr, sreq, MPI_STATUSES_IGNORE, ierr )
    deallocate ( sbuf, rbuf, sreq, rreq )
#endif
  end subroutine mr_seam_ranks


  subroutine mr_set_dt_blk ( dtl, lev, n, DT, b )
    use MOSE_Global_m, only: gc
    use MOSE_Field_m, only: dt_f, n_cache_f, field_block_index_params
    implicit none
    integer,  intent(in)    :: n(3), b, lev(n(1),n(2),n(3))
    real(R8), intent(in)    :: DT
    real(R8), intent(inout) :: dtl(n(1),n(2),n(3))
    integer  :: i, j, k, n1, n2, n3, boff_, pi_, pj_
    real(R8) :: dt_
    n1 = n(1); n2 = n(2); n3 = n(3); dt_ = DT
    call field_block_index_params ( b, boff_, pi_, pj_ )
    !$acc parallel loop collapse(3) gang vector copy(dtl) present(lev, dt_f)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      dtl(i,j,k) = dt_ / ( 2.0_R8 ** lev(i,j,k) )
      if ( n_cache_f > 1 ) dt_f( boff_ + (i+gc-1) + (j+gc-1)*pi_ + (k+gc-1)*pi_*pj_ + 1 ) = dtl(i,j,k)
    end do
    end do
    end do
  end subroutine mr_set_dt_blk


  !------------------------------------------------------------------------------------------------------------------------
  !> The boxes of every (local block, level) from the host levels: the bounding box of the level's cells, split along the
  !> emptiest slab (signature bisection) while that removes enough empty volume, at most mr-boxes boxes. The stage kernels of
  !> a level sweep its boxes; the cells of other levels inside a box are skipped by the predicate on dtlocal. At the first
  !> full partition the census line also gives the sweep the caps 8, 16, 32 and 64 would give (host only, not kept).
  subroutine mr_build_boxes ( domain, L )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block, mpi_is_root, mpi_allreduce_sum_r8_array
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: L
    integer, parameter :: NCAP = 4, CAPS(NCAP) = (/ 8, 16, 32, 64 /)
    integer  :: b, l_, ic, nb_, bx_(6, MR_NBOX_MAX)
    real(R8) :: vol(0:MR_LMAX_HARD), cells(0:MR_LMAX_HARD), volc(0:MR_LMAX_HARD, NCAP), nlch(0:MR_LMAX_HARD, NCAP)
    real(R8) :: buf(2*(MR_LMAX_HARD+1)), bufc(2*(MR_LMAX_HARD+1)*NCAP)
    logical  :: census
    census = ( obj_mr%nfull == 1 )
    vol = 0.0_R8; cells = 0.0_R8; volc = 0.0_R8; nlch = 0.0_R8
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      mrb(b)%nbox = 0
      do l_ = 0, L
        call mr_boxes_level ( mrb(b)%lev, domain%blk(b)%dim, l_, obj_mr%nboxmax, mrb(b)%nbox(l_), mrb(b)%box(:,:,l_) )
        cells(l_) = cells(l_) + real(count(mrb(b)%lev == l_), R8)
        vol(l_) = vol(l_) + boxes_vol ( mrb(b)%box(:,:,l_), mrb(b)%nbox(l_) )
        if ( census ) then
          do ic = 1, NCAP
            call mr_boxes_level ( mrb(b)%lev, domain%blk(b)%dim, l_, CAPS(ic), nb_, bx_ )
            volc(l_,ic) = volc(l_,ic) + boxes_vol ( bx_, nb_ ); nlch(l_,ic) = nlch(l_,ic) + real(nb_, R8)
          end do
        end if
      end do
    end do
    buf(1:MR_LMAX_HARD+1) = vol; buf(MR_LMAX_HARD+2:) = cells
    call mpi_allreduce_sum_r8_array ( buf, 2*(MR_LMAX_HARD+1) )
    if ( census ) then
      bufc(1:(MR_LMAX_HARD+1)*NCAP) = reshape( volc, (/ (MR_LMAX_HARD+1)*NCAP /) )
      bufc((MR_LMAX_HARD+1)*NCAP+1:) = reshape( nlch, (/ (MR_LMAX_HARD+1)*NCAP /) )
      call mpi_allreduce_sum_r8_array ( bufc, 2*(MR_LMAX_HARD+1)*NCAP )
      volc = reshape( bufc(1:(MR_LMAX_HARD+1)*NCAP), (/ MR_LMAX_HARD+1, NCAP /) )
      nlch = reshape( bufc((MR_LMAX_HARD+1)*NCAP+1:), (/ MR_LMAX_HARD+1, NCAP /) )
    end if
    if ( mpi_is_root .and. census ) then
      write(*,'(A,I0,A)',advance='no') ' [MR]   box sweep / level cells (mr-boxes ', obj_mr%nboxmax, '):'
      do l_ = 0, L
        write(*,'(A,I0,A,F7.3)',advance='no') '  L', l_, ' ', buf(l_+1) / max(buf(MR_LMAX_HARD+2+l_), 1.0_R8)
      end do
      write(*,'(A)') ''
      do ic = 1, NCAP
        write(*,'(A,I2,A)',advance='no') ' [MR]     with mr-boxes ', CAPS(ic), ': sweep / cells (boxes)'
        do l_ = 0, L
          write(*,'(A,I0,A,F7.3,A,I0,A)',advance='no') '  L', l_, ' ', volc(l_,ic) / max(buf(MR_LMAX_HARD+2+l_), 1.0_R8), &
            ' (', nint(nlch(l_,ic)), ')'
        end do
        write(*,'(A)') ''
      end do
    end if

  contains

    real(R8) function boxes_vol ( bx, nbx )
      integer, intent(in) :: bx(6, MR_NBOX_MAX), nbx
      integer :: ib
      boxes_vol = 0.0_R8
      do ib = 1, nbx
        boxes_vol = boxes_vol + real(bx(2,ib)-bx(1,ib)+1, R8) * real(bx(4,ib)-bx(3,ib)+1, R8) * real(bx(6,ib)-bx(5,ib)+1, R8)
      end do
    end function boxes_vol

  end subroutine mr_build_boxes

  subroutine mr_boxes_level ( lev, n, l, nmax, nbox, box )
    implicit none
    integer, intent(in)  :: n(3), lev(n(1),n(2),n(3)), l, nmax
    integer, intent(out) :: nbox, box(6, MR_NBOX_MAX)
    integer  :: cand(6, MR_NBOX_MAX), ncnt(MR_NBOX_MAX), done(MR_NBOX_MAX), ncand, c, cbest, nc_, bx(6), b1(6), b2(6), n1_, n2_
    real(R8) :: waste, wbest
    logical  :: ok
    nbox = 0
    bx = (/ 1, n(1), 1, n(2), 1, n(3) /)
    call tight ( bx, nc_ )
    if ( nc_ == 0 ) return
    ncand = 1; cand(:,1) = bx; ncnt(1) = nc_; done(1) = 0
    do
      if ( ncand >= nmax ) exit
      cbest = 0; wbest = 0.0_R8          ! the open box with the most empty volume
      do c = 1, ncand
        if ( done(c) == 1 ) cycle
        waste = boxvol(cand(:,c)) - real(ncnt(c), R8)
        if ( waste > wbest ) then
          wbest = waste; cbest = c
        end if
      end do
      if ( cbest == 0 ) exit
      call split ( cand(:,cbest), ncnt(cbest), b1, n1_, b2, n2_, ok )
      if ( .not. ok ) then
        done(cbest) = 1; cycle
      end if
      cand(:,cbest) = b1; ncnt(cbest) = n1_; done(cbest) = 0
      ncand = ncand + 1; cand(:,ncand) = b2; ncnt(ncand) = n2_; done(ncand) = 0
    end do
    nbox = ncand; box(:,1:ncand) = cand(:,1:ncand)

  contains

    real(R8) function boxvol ( q )
      integer, intent(in) :: q(6)
      boxvol = real(q(2)-q(1)+1, R8) * real(q(4)-q(3)+1, R8) * real(q(6)-q(5)+1, R8)
    end function boxvol

    !> Shrink q to the bounding box of its level-l cells; nc = their number (q unchanged when 0).
    subroutine tight ( q, nc )
      integer, intent(inout) :: q(6)
      integer, intent(out)   :: nc
      integer :: i, j, k, t(6)
      t = (/ q(2)+1, q(1)-1, q(4)+1, q(3)-1, q(6)+1, q(5)-1 /); nc = 0
      do k = q(5), q(6); do j = q(3), q(4); do i = q(1), q(2)
        if ( lev(i,j,k) /= l ) cycle
        nc = nc + 1
        t(1) = min(t(1), i); t(2) = max(t(2), i); t(3) = min(t(3), j); t(4) = max(t(4), j); t(5) = min(t(5), k); t(6) = max(t(6), k)
      end do; end do; end do
      if ( nc > 0 ) q = t
    end subroutine tight

    !> Cut the tight box q (nq level-l cells) in two tight boxes along the axis and slab that leave the least empty volume:
    !> per slab of each axis the count and the extent of the level's cells in the other two axes, then prefix and suffix
    !> unions, so every cut costs O(1). ok when the cut removes at least a fifth of the empty volume of q.
    subroutine split ( q, nq, q1, nq1, q2, nq2, ok )
      integer, intent(in)  :: q(6), nq
      integer, intent(out) :: q1(6), nq1, q2(6), nq2
      logical, intent(out) :: ok
      integer  :: d, da, db, p, i, j, k, lo, hi, s, ca, cb, first, last
      integer, allocatable :: cnt(:), amn(:), amx(:), bmn(:), bmx(:)
      integer, allocatable :: pc(:), pamn(:), pamx(:), pbmn(:), pbmx(:), plast(:)
      integer, allocatable :: sc(:), samn(:), samx(:), sbmn(:), sbmx(:), sfirst(:)
      real(R8) :: w0, w, wb
      w0 = boxvol(q) - real(nq, R8)
      wb = w0; ok = .false.; nq1 = 0; nq2 = 0; q1 = q; q2 = q
      do d = 1, 3
        lo = q(2*d-1); hi = q(2*d)
        if ( hi - lo < 1 ) cycle
        da = mod(d, 3) + 1; db = mod(d+1, 3) + 1
        allocate ( cnt(lo:hi), amn(lo:hi), amx(lo:hi), bmn(lo:hi), bmx(lo:hi) )
        cnt = 0; amn = huge(1); amx = -huge(1); bmn = huge(1); bmx = -huge(1)
        do k = q(5), q(6); do j = q(3), q(4); do i = q(1), q(2)
          if ( lev(i,j,k) /= l ) cycle
          select case (d)
            case (1); s = i; ca = j; cb = k
            case (2); s = j; ca = k; cb = i
            case default; s = k; ca = i; cb = j
          end select
          cnt(s) = cnt(s) + 1
          amn(s) = min(amn(s), ca); amx(s) = max(amx(s), ca); bmn(s) = min(bmn(s), cb); bmx(s) = max(bmx(s), cb)
        end do; end do; end do
        allocate ( pc(lo:hi), pamn(lo:hi), pamx(lo:hi), pbmn(lo:hi), pbmx(lo:hi), plast(lo:hi) )
        allocate ( sc(lo:hi), samn(lo:hi), samx(lo:hi), sbmn(lo:hi), sbmx(lo:hi), sfirst(lo:hi) )
        do p = lo, hi
          if ( p == lo ) then
            pc(p) = cnt(p); pamn(p) = amn(p); pamx(p) = amx(p); pbmn(p) = bmn(p); pbmx(p) = bmx(p); plast(p) = lo - 1
          else
            pc(p) = pc(p-1) + cnt(p); pamn(p) = min(pamn(p-1), amn(p)); pamx(p) = max(pamx(p-1), amx(p))
            pbmn(p) = min(pbmn(p-1), bmn(p)); pbmx(p) = max(pbmx(p-1), bmx(p)); plast(p) = plast(p-1)
          end if
          if ( cnt(p) > 0 ) plast(p) = p
        end do
        do p = hi, lo, -1
          if ( p == hi ) then
            sc(p) = cnt(p); samn(p) = amn(p); samx(p) = amx(p); sbmn(p) = bmn(p); sbmx(p) = bmx(p); sfirst(p) = hi + 1
          else
            sc(p) = sc(p+1) + cnt(p); samn(p) = min(samn(p+1), amn(p)); samx(p) = max(samx(p+1), amx(p))
            sbmn(p) = min(sbmn(p+1), bmn(p)); sbmx(p) = max(sbmx(p+1), bmx(p)); sfirst(p) = sfirst(p+1)
          end if
          if ( cnt(p) > 0 ) sfirst(p) = p
        end do
        first = sfirst(lo); last = plast(hi)
        do p = lo, hi - 1           ! cut between slabs p and p+1
          if ( pc(p) == 0 .or. sc(p+1) == 0 ) cycle
          w = real(plast(p) - first + 1, R8) * real(pamx(p) - pamn(p) + 1, R8) * real(pbmx(p) - pbmn(p) + 1, R8) &
            - real(pc(p), R8) &
            + real(last - sfirst(p+1) + 1, R8) * real(samx(p+1) - samn(p+1) + 1, R8) * real(sbmx(p+1) - sbmn(p+1) + 1, R8) &
            - real(sc(p+1), R8)
          if ( w < wb ) then
            wb = w; ok = .true.; nq1 = pc(p); nq2 = sc(p+1)
            q1(2*d-1) = first;        q1(2*d) = plast(p)
            q1(2*da-1) = pamn(p);     q1(2*da) = pamx(p)
            q1(2*db-1) = pbmn(p);     q1(2*db) = pbmx(p)
            q2(2*d-1) = sfirst(p+1);  q2(2*d) = last
            q2(2*da-1) = samn(p+1);   q2(2*da) = samx(p+1)
            q2(2*db-1) = sbmn(p+1);   q2(2*db) = sbmx(p+1)
          end if
        end do
        deallocate ( cnt, amn, amx, bmn, bmx, pc, pamn, pamx, pbmn, pbmx, plast, sc, samn, samx, sbmn, sbmx, sfirst )
      end do
      if ( ok ) ok = ( wb <= 0.8_R8 * w0 )
    end subroutine split

  end subroutine mr_boxes_level


  !> The interface faces of every local block (faces between levels l and l+1; the rules allow no larger jump) and the
  !> telescoping accumulators, zeroed: a full partition runs between two macro steps, when no accumulation is pending.
  !> The boundary records of each level by launch of BC_Fluxes_blk (Mod_BC bc_launch_ranges): a stage of level l runs the
  !> launches of its block over bcl(bco(g,l)+1 : bco(g,l)+bcn(g,l)), the records whose interior cell is of level l, in the
  !> order of their range (the set and the order of the predicate on dtlocal, without sweeping the other records).
  subroutine mr_build_bclist ( domain, L )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Mod_BC,  only: bc_launch_ranges
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: L
    integer :: b, g, ng, glo(6), ghi(6), idx, l_, m
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      if ( allocated(mrb(b)%bcl) ) then
        !$acc exit data delete(mrb(b)%bcl)
        deallocate ( mrb(b)%bcl )
      end if
      mrb(b)%bcn = 0; mrb(b)%bco = 0
      ng = 0
      if ( domain % bc_blk_hi(b) >= domain % bc_blk_lo(b) ) &
        call bc_launch_ranges ( domain%blk(b)%dim, domain % n_bf(b,1:6), domain % bc_blk_lo(b), ng, glo, ghi )
      allocate ( mrb(b)%bcl(max(sum(ghi(1:ng) - glo(1:ng) + 1), 1)) )
      mrb(b)%bcl = 0; m = 0
      do g = 1, ng
        do l_ = 0, L
          mrb(b)%bco(g,l_) = m
          do idx = glo(g), ghi(g)
            if ( mrb(b)%lev(domain%bc_bi(idx), domain%bc_bj(idx), domain%bc_bk(idx)) /= l_ ) cycle
            m = m + 1; mrb(b)%bcl(m) = idx
          end do
          mrb(b)%bcn(g,l_) = m - mrb(b)%bco(g,l_)
        end do
      end do
      !$acc enter data copyin(mrb(b)%bcl)
    end do
  end subroutine mr_build_bclist


  !> The lists of the list path (mr-list-level, Lib_Multirate_Base mr_fl_type / mr_cl_type) for every level l >= listlev:
  !> per direction the interior faces with a cell of level l on either side -- the faces whose fluxes the level's cells
  !> read --, each once, in block order and in each block in the face order of the kernels (i fastest); per face the flat
  !> index of its left cell, the flat strides of the diffusive stencil and the block geometry the kernels read (the
  !> convective widths dl(d, -1..2), the face metric, the cell metrics and wall distances of the two cells); the cells of
  !> level l with the list positions of their six faces; the interface records of the level (MR_Iface_Update) by list
  !> position. Rebuilt with the levels; beta is set to 1 once (the list kernels take beta = 1: no shock detector).
  subroutine mr_build_lists ( domain, L, cp_ )
    use MOSE_Types_m
    use MOSE_Global_m,   only: gc, nprim, nsc
    use MOSE_Mod_MPI,    only: is_local_block, mpi_is_root
    use MOSE_Field_m,    only: field_block_index_params, n_cache_f
    use MOSE_Lib_RANS,   only: mcdiff_dev
    use MOSE_Lib_Limiters, only: SD_space_order
    use MOSE_Lib_Shock_Detector, only: sd_dev
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: L
    integer(8), intent(inout) :: cp_   ! mr-profile clock: faces and cells (tpb 4), ghost (5), boundary (6), thermo (7) lists
    logical, save :: beta_done = .false.
    integer :: l_, d, b, i, j, k, n1, n2, n3, boff, pi, pj, nf(3), ncl, e, ec, ef(3), q, r, ntf, ni(3), ie(3)
    integer :: la, lb, di, dj, dk, f
    integer, allocatable :: pos1(:), pos2(:), pos3(:)
    if ( SD_space_order .or. sd_dev > 0 ) error stop 'mr-list-level: the list kernels take beta = 1 (no shock detector)'
    if ( n_cache_f <= 1 ) error stop 'mr-list-level: the list Newstate seeds T from the flat thermo cache (viscous runs)'
    if ( .not. allocated(mrfl) ) then
      allocate ( mrfl(3, 0:MR_LMAX_HARD), mrcl(0:MR_LMAX_HARD), mrgl(0:MR_LMAX_HARD), mrbl(6, 0:MR_LMAX_HARD), &
                 mrtl(0:MR_LMAX_HARD) )
      call mr_setup_dn ( domain )
    end if
    if ( .not. beta_done ) then
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        call mr_beta_one ( domain%blk(b)%beta, domain%blk(b)%dim )
      end do
      beta_done = .true.
    end if
    ntf = 5
    if ( mcdiff_dev ) ntf = 5 + nsc
    do l_ = 0, MR_LMAX_HARD
      call mr_free_lists ( l_ )
    end do
    do l_ = max(mr_rec_level(), 0), L
      if ( l_ < obj_mr%listlev ) then   ! mr-list-records-level: the record lists only (faces and cells on the boxes)
        call mr_build_ghost_list ( domain, l_ )
        call mr_build_bc_list ( domain, l_ )
        if ( mpi_is_root .and. obj_mr%nfull == 1 ) write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0)') ' [MR]   record lists of level ', &
          l_, ': ghost records ', mrgl(l_)%np, ' + connection ', mrgl(l_)%nc, ' (gather ', mrgl(l_)%npb, ')  410 slots ', &
          mrgl(l_)%nq
        cycle
      end if
      ! count
      if ( mr_list_dev_on() ) then   ! mr-list-dev: the face and cell lists built on the device (same order)
        call mr_build_lists_dev ( domain, l_, ntf )
      else
      nf = 0; ncl = 0; ni = 0
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        n1 = domain%blk(b)%dim(1); n2 = domain%blk(b)%dim(2); n3 = domain%blk(b)%dim(3)
        do k = 1, n3; do j = 1, n2; do i = 1, n1
          la = mrb(b)%lev(i,j,k)
          if ( la == l_ ) ncl = ncl + 1
          if ( i < n1 ) then; if ( la == l_ .or. mrb(b)%lev(i+1,j,k) == l_ ) nf(1) = nf(1) + 1; end if
          if ( j < n2 ) then; if ( la == l_ .or. mrb(b)%lev(i,j+1,k) == l_ ) nf(2) = nf(2) + 1; end if
          if ( k < n3 ) then; if ( la == l_ .or. mrb(b)%lev(i,j,k+1) == l_ ) nf(3) = nf(3) + 1; end if
        end do; end do; end do
        do d = 1, 3
          do q = 1, mrb(b)%ifl(d)%n
            if ( mrb(b)%ifl(d)%lc(q) == l_ .or. mrb(b)%ifl(d)%lc(q) + 1 == l_ ) ni(d) = ni(d) + 1
          end do
        end do
      end do
      do d = 1, 3
        associate ( fl => mrfl(d, l_) )
          fl%n = nf(d); fl%ni = ni(d)
          r = max(nf(d), 1)
          allocate ( fl%c(r), fl%sa(r), fl%sb(r), fl%sc(r), fl%dl(4,r), fl%fm(4,r), fl%mm(3,3,2,r), fl%yn(2,r), &
                     fl%fh(r,nprim), fl%dfh(r,nprim), fl%tf(r,ntf), &
                     fl%if_(max(ni(d),1)), fl%ilc(max(ni(d),1)), fl%isl(max(ni(d),1)) )
          fl%c = 0; fl%sa = 0; fl%sb = 0; fl%sc = 0; fl%if_ = 1; fl%ilc = 0; fl%isl = 1
        end associate
      end do
      mrcl(l_)%n = ncl
      r = max(ncl,1)
      allocate ( mrcl(l_)%c(r), mrcl(l_)%pf(6,r), mrcl(l_)%ijk(3,r), mrcl(l_)%bo(r), mrcl(l_)%pp(2,r), &
                 mrcl(l_)%vol(r), mrcl(l_)%met(3,3,r), mrcl(l_)%wd(r) )
      mrcl(l_)%c = 0; mrcl(l_)%pf = 0; mrcl(l_)%ijk = 1; mrcl(l_)%bo = 0; mrcl(l_)%pp = 1
      mrcl(l_)%vol = 1.0_R8; mrcl(l_)%met = 0.0_R8; mrcl(l_)%wd = 0.0_R8
      ! fill
      ef = 0; ec = 0; ie = 0
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        n1 = domain%blk(b)%dim(1); n2 = domain%blk(b)%dim(2); n3 = domain%blk(b)%dim(3)
        call field_block_index_params ( b, boff, pi, pj )
        allocate ( pos1(max((n1-1)*n2*n3,1)), pos2(max(n1*(n2-1)*n3,1)), pos3(max(n1*n2*(n3-1),1)) )
        pos1 = 0; pos2 = 0; pos3 = 0
        do d = 1, 3
          di = 0; dj = 0; dk = 0
          if ( d == 1 ) di = 1
          if ( d == 2 ) dj = 1
          if ( d == 3 ) dk = 1
          do k = 1, n3 - dk; do j = 1, n2 - dj; do i = 1, n1 - di
            if ( mrb(b)%lev(i,j,k) /= l_ .and. mrb(b)%lev(i+di,j+dj,k+dk) /= l_ ) cycle
            ef(d) = ef(d) + 1; e = ef(d)
            associate ( fl => mrfl(d, l_), bk => domain%blk(b) )
              fl%c(e) = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
              if ( d == 1 ) then
                f = (i-1) + (j-1)*(n1-1) + (k-1)*(n1-1)*n2 + 1; pos1(f) = e
                fl%sa(e) = 1; fl%sb(e) = pi; fl%sc(e) = pi*pj
                fl%dl(1:4,e) = (/ bk%dl_f(1,i-1,j,k), bk%dl_f(1,i,j,k), bk%dl_f(1,i+1,j,k), bk%dl_f(1,i+2,j,k) /)
                fl%fm(1:4,e) = bk%fm1(1:4,i,j,k)
              else if ( d == 2 ) then
                f = (i-1) + (j-1)*n1 + (k-1)*n1*(n2-1) + 1; pos2(f) = e
                fl%sa(e) = pi; fl%sb(e) = 1; fl%sc(e) = pi*pj
                fl%dl(1:4,e) = (/ bk%dl_f(2,i,j-1,k), bk%dl_f(2,i,j,k), bk%dl_f(2,i,j+1,k), bk%dl_f(2,i,j+2,k) /)
                fl%fm(1:4,e) = bk%fm2(1:4,i,j,k)
              else
                f = (i-1) + (j-1)*n1 + (k-1)*n1*n2 + 1; pos3(f) = e
                fl%sa(e) = pi*pj; fl%sb(e) = 1; fl%sc(e) = pi
                fl%dl(1:4,e) = (/ bk%dl_f(3,i,j,k-1), bk%dl_f(3,i,j,k), bk%dl_f(3,i,j,k+1), bk%dl_f(3,i,j,k+2) /)
                fl%fm(1:4,e) = bk%fm3(1:4,i,j,k)
              end if
              fl%mm(:,:,1,e) = bk%m_f(:,:,i,j,k); fl%mm(:,:,2,e) = bk%m_f(:,:,i+di,j+dj,k+dk)
              fl%yn(1,e) = bk%yn(i,j,k); fl%yn(2,e) = bk%yn(i+di,j+dj,k+dk)
            end associate
          end do; end do; end do
        end do
        do k = 1, n3; do j = 1, n2; do i = 1, n1
          if ( mrb(b)%lev(i,j,k) /= l_ ) cycle
          ec = ec + 1
          mrcl(l_)%c(ec) = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
          mrcl(l_)%ijk(:,ec) = (/ i, j, k /); mrcl(l_)%bo(ec) = boff; mrcl(l_)%pp(:,ec) = (/ pi, pj /)
          mrcl(l_)%vol(ec) = domain%blk(b)%vol(i,j,k); mrcl(l_)%met(:,:,ec) = domain%blk(b)%m_f(:,:,i,j,k)
          mrcl(l_)%wd(ec) = domain%blk(b)%yn(i,j,k)
          if ( i < n1 ) mrcl(l_)%pf(1,ec) = pos1( (i-1) + (j-1)*(n1-1) + (k-1)*(n1-1)*n2 + 1 )
          if ( i > 1  ) mrcl(l_)%pf(2,ec) = pos1( (i-2) + (j-1)*(n1-1) + (k-1)*(n1-1)*n2 + 1 )
          if ( j < n2 ) mrcl(l_)%pf(3,ec) = pos2( (i-1) + (j-1)*n1 + (k-1)*n1*(n2-1) + 1 )
          if ( j > 1  ) mrcl(l_)%pf(4,ec) = pos2( (i-1) + (j-2)*n1 + (k-1)*n1*(n2-1) + 1 )
          if ( k < n3 ) mrcl(l_)%pf(5,ec) = pos3( (i-1) + (j-1)*n1 + (k-1)*n1*n2 + 1 )
          if ( k > 1  ) mrcl(l_)%pf(6,ec) = pos3( (i-1) + (j-1)*n1 + (k-2)*n1*n2 + 1 )
        end do; end do; end do
        do d = 1, 3
          do q = 1, mrb(b)%ifl(d)%n
            if ( mrb(b)%ifl(d)%lc(q) /= l_ .and. mrb(b)%ifl(d)%lc(q) + 1 /= l_ ) cycle
            ie(d) = ie(d) + 1
            f = mrb(b)%ifl(d)%f(q)
            if ( d == 1 ) mrfl(d,l_)%if_(ie(d)) = pos1(f)
            if ( d == 2 ) mrfl(d,l_)%if_(ie(d)) = pos2(f)
            if ( d == 3 ) mrfl(d,l_)%if_(ie(d)) = pos3(f)
            if ( mrfl(d,l_)%if_(ie(d)) <= 0 ) error stop 'mr_build_lists: an interface face of the level is not in its list'
            mrfl(d,l_)%ilc(ie(d)) = mrb(b)%ifl(d)%lc(q)
            mrfl(d,l_)%isl(ie(d)) = mrb(b)%ifl(d)%s(q)
          end do
        end do
        deallocate ( pos1, pos2, pos3 )
      end do
      if ( any(ef /= nf) .or. ec /= ncl .or. any(ie /= ni) ) error stop 'mr_build_lists: count and fill differ'
      if ( any(mrcl(l_)%pf(:,1:ncl) < 0) ) error stop 'mr_build_lists: a face position is negative'
      do d = 1, 3
        !$acc enter data copyin(mrfl(d,l_)%c, mrfl(d,l_)%sa, mrfl(d,l_)%sb, mrfl(d,l_)%sc, mrfl(d,l_)%dl, mrfl(d,l_)%fm) &
        !$acc&           copyin(mrfl(d,l_)%mm, mrfl(d,l_)%yn, mrfl(d,l_)%if_, mrfl(d,l_)%ilc, mrfl(d,l_)%isl) &
        !$acc&           create(mrfl(d,l_)%fh, mrfl(d,l_)%dfh, mrfl(d,l_)%tf)
      end do
      !$acc enter data copyin(mrcl(l_)%c, mrcl(l_)%pf, mrcl(l_)%ijk, mrcl(l_)%bo, mrcl(l_)%pp, mrcl(l_)%vol, mrcl(l_)%met, &
      !$acc&                  mrcl(l_)%wd)
      end if   ! mr-list-dev
      if ( obj_mr%profile ) call mr_ptoc ( obj_mr%tpb(4), cp_ )
      call mr_build_bc_list ( domain, l_ )   ! (independent of the ghost list: built first, its records feed the read mark)
      if ( obj_mr%profile ) call mr_ptoc ( obj_mr%tpb(6), cp_ )
      rk_ready = .false.
      if ( mr_ghost_exact_on() ) call mr_mark_reads_dev ( domain, l_ )   ! mr-ghost-exact: the primitives the stage reads
      call mr_build_ghost_list ( domain, l_ )
      rk_ready = .false.
      if ( obj_mr%profile ) call mr_ptoc ( obj_mr%tpb(5), cp_ )
      if ( mr_thermo_exact_on() ) then   ! mr-thermo-exact: the cells the stage reads, built on the device
        call mr_build_thermo_list_dev ( domain, l_ )
      else
        call mr_build_thermo_list ( domain, l_ )
      end if
      if ( obj_mr%profile ) call mr_ptoc ( obj_mr%tpb(7), cp_ )
      if ( mpi_is_root .and. obj_mr%nfull == 1 ) write(*,'(A,I0,A,I0,A,3(1X,I0),A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
        ' [MR]   list level ', l_, ': cells ', mrcl(l_)%n, '  faces i/j/k', mrfl(1,l_)%n, mrfl(2,l_)%n, mrfl(3,l_)%n, &
        '  interface faces ', mrfl(1,l_)%ni + mrfl(2,l_)%ni + mrfl(3,l_)%ni, &
        '  ghost records ', mrgl(l_)%np, ' + connection ', mrgl(l_)%nc, ' (gather ', mrgl(l_)%npb, ') + self pairs ', &
        mrgl(l_)%nself, '  410 slots ', mrgl(l_)%nq
    end do
  end subroutine mr_build_lists

  !> mr-list-dev: MOSE_MR_LIST_DEV (default on) builds the face and cell lists of a level on the device.
  logical function mr_list_dev_on ()
    logical, save :: checked = .false., on = .true.
    character(len=8) :: env
    integer :: st, ln
    if ( .not. checked ) then
      checked = .true.
      call get_environment_variable ( 'MOSE_MR_LIST_DEV', env, ln, st )
      if ( st == 0 .and. ln > 0 ) then
        if ( env(1:1) == '0' .or. env(1:1) == 'n' .or. env(1:1) == 'N' .or. env(1:1) == 'f' .or. env(1:1) == 'F' ) on = .false.
      end if
    end if
    mr_list_dev_on = on
  end function mr_list_dev_on

  !> mr-build-fast: MOSE_MR_BUILD_FAST (default on; 0/n/N/f/F = off) -- the in-chunk prefix of ld_scan by sub-chunks, the
  !> counts of the list build and the changed cells of the nesting rule summed on the device with one read back. Integers
  !> only: the same levels, lists and positions either way (the switch is for A/B timing).
  logical function mr_build_fast_on ()
    logical, save :: checked = .false., on = .true.
    character(len=8) :: env
    integer :: st, ln
    if ( .not. checked ) then
      checked = .true.
      call get_environment_variable ( 'MOSE_MR_BUILD_FAST', env, ln, st )
      if ( st == 0 .and. ln > 0 ) then
        if ( env(1:1) == '0' .or. env(1:1) == 'n' .or. env(1:1) == 'N' .or. env(1:1) == 'f' .or. env(1:1) == 'F' ) on = .false.
      end if
    end if
    mr_build_fast_on = on
  end function mr_build_fast_on

  !> The face and cell lists of level l_ built on the device, in the order of the host build (blocks ascending, faces and
  !> cells in lexicographic i-j-k order within a block): per block, a count kernel (the sizes, read back once), then per
  !> direction a flag kernel, an exclusive prefix sum (chunks of 1024: sums, their scan, the in-chunk prefix) and a gather
  !> kernel that writes the list entries at their positions and the face -> position table of the block (pos1/2/3, on the
  !> device) for the six-face links of the cells and the interface faces. The list arrays are allocated on the host at
  !> their sizes and created on the device without a copy (the host copies hold nothing; no code reads them). The host
  !> build (MOSE_MR_LIST_DEV=0) stays as the reference: the same lists to the element.
  subroutine mr_build_lists_dev ( domain, l_, ntf )
    use MOSE_Types_m
    use MOSE_Global_m,   only: gc, nprim
    use MOSE_Mod_MPI,    only: is_local_block
    use MOSE_Field_m,    only: field_block_index_params
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: l_, ntf
    integer :: b, d, n1, n2, n3, boff, pi, pj, nmax, nf(3), ncl, ni(3), nface, r, nifl, i
    integer :: cnt(7), ofs_f(3), ofs_c, ofs_i(3)
    integer, allocatable :: cntb(:,:)   ! (7, nb): nf1 nf2 nf3 ncl ni1 ni2 ni3 of each local block

    ! scratch on the device: flags, positions, chunk sums, the three face->position tables of the current block
    nmax = 1
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      nmax = max( nmax, product(domain%blk(b)%dim) )
    end do
    if ( allocated(ld_flg) ) then
      if ( size(ld_flg) < nmax ) then
        !$acc exit data delete(ld_flg, ld_pos, ld_p1, ld_p2, ld_p3, ld_csum, ld_coff, ld_ssum)
        deallocate ( ld_flg, ld_pos, ld_p1, ld_p2, ld_p3, ld_csum, ld_coff, ld_ssum )
      end if
    end if
    if ( .not. allocated(ld_flg) ) then
      allocate ( ld_flg(nmax), ld_pos(nmax), ld_p1(nmax), ld_p2(nmax), ld_p3(nmax), ld_csum(nmax/LD_CH+2), ld_coff(nmax/LD_CH+2), &
                 ld_ssum(nmax/LD_SC+2) )
      !$acc enter data create(ld_flg, ld_pos, ld_p1, ld_p2, ld_p3, ld_csum, ld_coff, ld_ssum)
    end if
    if ( .not. allocated(ld_cnt) ) then
      allocate ( ld_cnt(7, domain % nb) )
      !$acc enter data create(ld_cnt)
    end if

    ld_q = domain % nb + 65   ! one queue for every list kernel of this level: in order, one wait at the end
    ! 1. the sizes: per block (7 reductions), then the totals and the offsets of each block in the level's lists
    allocate ( cntb(7, domain % nb) ); cntb = 0
    if ( mr_build_fast_on() ) then   ! mr-build-fast: every block's counts on the queue into ld_cnt, one read back
      call ld_set_i ( ld_cnt, 7 * domain % nb, 0 )
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        call ld_count_blk_q ( mrb(b)%lev, domain%blk(b)%dim, l_, ld_cnt, domain % nb, b )
        do d = 1, 3
          if ( mrb(b)%ifl(d)%n > 0 ) call ld_count_ifl_q ( mrb(b)%ifl(d)%lc, mrb(b)%ifl(d)%n, l_, ld_cnt, domain % nb, b, 4+d )
        end do
      end do
      !$acc update host(ld_cnt) async(ld_q)
      !$acc wait(ld_q)
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        cntb(:, b) = ld_cnt(:, b)
      end do
    else
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      nifl = max( mrb(b)%ifl(1)%n, mrb(b)%ifl(2)%n, mrb(b)%ifl(3)%n )
      call ld_count_blk ( mrb(b)%lev, domain%blk(b)%dim, l_, cnt )
      do d = 1, 3
        if ( mrb(b)%ifl(d)%n > 0 ) then
          call ld_count_ifl ( mrb(b)%ifl(d)%lc, mrb(b)%ifl(d)%n, l_, cnt(4+d) )
        else
          cnt(4+d) = 0
        end if
      end do
      cntb(:, b) = cnt
    end do
    end if   ! mr-build-fast
    nf = 0; ncl = 0; ni = 0
    do b = 1, domain % nb
      nf = nf + cntb(1:3, b); ncl = ncl + cntb(4, b); ni = ni + cntb(5:7, b)
    end do

    ! 2. the lists: host allocation at size (unfilled), device creation
    do d = 1, 3
      associate ( fl => mrfl(d, l_) )
        fl%n = nf(d); fl%ni = ni(d)
        r = max(nf(d), 1)
        allocate ( fl%c(r), fl%sa(r), fl%sb(r), fl%sc(r), fl%dl(4,r), fl%fm(4,r), fl%mm(3,3,2,r), fl%yn(2,r), &
                   fl%fh(r,nprim), fl%dfh(r,nprim), fl%tf(r,ntf), &
                   fl%if_(max(ni(d),1)), fl%ilc(max(ni(d),1)), fl%isl(max(ni(d),1)) )
      end associate
      !$acc enter data create(mrfl(d,l_)%c, mrfl(d,l_)%sa, mrfl(d,l_)%sb, mrfl(d,l_)%sc, mrfl(d,l_)%dl, mrfl(d,l_)%fm) &
      !$acc&           create(mrfl(d,l_)%mm, mrfl(d,l_)%yn, mrfl(d,l_)%if_, mrfl(d,l_)%ilc, mrfl(d,l_)%isl) &
      !$acc&           create(mrfl(d,l_)%fh, mrfl(d,l_)%dfh, mrfl(d,l_)%tf)
      associate ( fl => mrfl(d, l_) )
        if ( nf(d) == 0 ) then   ! the padding entry of an empty list, as the host build leaves it
          call ld_set_i ( fl%c, 1, 0 ); call ld_set_i ( fl%sa, 1, 0 ); call ld_set_i ( fl%sb, 1, 0 ); call ld_set_i ( fl%sc, 1, 0 )
        end if
        if ( ni(d) == 0 ) then
          call ld_set_i ( fl%if_, 1, 1 ); call ld_set_i ( fl%ilc, 1, 0 ); call ld_set_i ( fl%isl, 1, 1 )
        end if
      end associate
    end do
    mrcl(l_)%n = ncl
    r = max(ncl, 1)
    allocate ( mrcl(l_)%c(r), mrcl(l_)%pf(6,r), mrcl(l_)%ijk(3,r), mrcl(l_)%bo(r), mrcl(l_)%pp(2,r), &
               mrcl(l_)%vol(r), mrcl(l_)%met(3,3,r), mrcl(l_)%wd(r) )
    !$acc enter data create(mrcl(l_)%c, mrcl(l_)%pf, mrcl(l_)%ijk, mrcl(l_)%bo, mrcl(l_)%pp, mrcl(l_)%vol, mrcl(l_)%met, &
    !$acc&                  mrcl(l_)%wd)
    if ( ncl == 0 ) then
      call ld_set_i ( mrcl(l_)%c, 1, 0 ); call ld_set_i ( mrcl(l_)%pf, 6, 0 ); call ld_set_i ( mrcl(l_)%ijk, 3, 1 )
      call ld_set_i ( mrcl(l_)%bo, 1, 0 ); call ld_set_i ( mrcl(l_)%pp, 2, 1 )
      call ld_set_r ( mrcl(l_)%vol, 1, 1.0_R8 ); call ld_set_r ( mrcl(l_)%met, 9, 0.0_R8 ); call ld_set_r ( mrcl(l_)%wd, 1, 0.0_R8 )
    end if

    ! 3. the entries: per block, flags -> prefix -> gather, faces then cells (their links read the face tables)
    ofs_f = 0; ofs_c = 0; ofs_i = 0
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      n1 = domain%blk(b)%dim(1); n2 = domain%blk(b)%dim(2); n3 = domain%blk(b)%dim(3)
      call field_block_index_params ( b, boff, pi, pj )
      do d = 1, 3
        if ( d == 1 ) nface = (n1-1)*n2*n3
        if ( d == 2 ) nface = n1*(n2-1)*n3
        if ( d == 3 ) nface = n1*n2*(n3-1)
        if ( nface <= 0 ) cycle
        call ld_flag_faces ( mrb(b)%lev, n1, n2, n3, d, l_, ld_flg, nface )
        call ld_scan ( ld_flg, ld_pos, nface, ld_csum, ld_coff, size(ld_csum), ld_ssum, size(ld_ssum) )
        if ( d == 1 ) then
          call ld_gather_faces ( mrb(b)%lev, domain%blk(b)%dl_f, domain%blk(b)%fm1, 0, 1, 1, domain%blk(b)%m_f, domain%blk(b)%yn, &
                                 n1, n2, n3, d, boff, pi, pj, gc, ofs_f(d), ld_flg, ld_pos, ld_p1, nface, &
                                 mrfl(d,l_)%c, mrfl(d,l_)%sa, mrfl(d,l_)%sb, mrfl(d,l_)%sc, mrfl(d,l_)%dl, mrfl(d,l_)%fm, &
                                 mrfl(d,l_)%mm, mrfl(d,l_)%yn, max(nf(d),1) )
        else if ( d == 2 ) then
          call ld_gather_faces ( mrb(b)%lev, domain%blk(b)%dl_f, domain%blk(b)%fm2, 1, 0, 1, domain%blk(b)%m_f, domain%blk(b)%yn, &
                                 n1, n2, n3, d, boff, pi, pj, gc, ofs_f(d), ld_flg, ld_pos, ld_p2, nface, &
                                 mrfl(d,l_)%c, mrfl(d,l_)%sa, mrfl(d,l_)%sb, mrfl(d,l_)%sc, mrfl(d,l_)%dl, mrfl(d,l_)%fm, &
                                 mrfl(d,l_)%mm, mrfl(d,l_)%yn, max(nf(d),1) )
        else
          call ld_gather_faces ( mrb(b)%lev, domain%blk(b)%dl_f, domain%blk(b)%fm3, 1, 1, 0, domain%blk(b)%m_f, domain%blk(b)%yn, &
                                 n1, n2, n3, d, boff, pi, pj, gc, ofs_f(d), ld_flg, ld_pos, ld_p3, nface, &
                                 mrfl(d,l_)%c, mrfl(d,l_)%sa, mrfl(d,l_)%sb, mrfl(d,l_)%sc, mrfl(d,l_)%dl, mrfl(d,l_)%fm, &
                                 mrfl(d,l_)%mm, mrfl(d,l_)%yn, max(nf(d),1) )
        end if
        ! the interface faces of this direction: positions in the level's list
        if ( mrb(b)%ifl(d)%n > 0 .and. cntb(4+d, b) > 0 ) then
          call ld_flag_ifl ( mrb(b)%ifl(d)%lc, mrb(b)%ifl(d)%n, l_, ld_flg )
          call ld_scan ( ld_flg, ld_pos, mrb(b)%ifl(d)%n, ld_csum, ld_coff, size(ld_csum), ld_ssum, size(ld_ssum) )
          if ( d == 1 ) then
            call ld_gather_ifl ( mrb(b)%ifl(d)%f, mrb(b)%ifl(d)%lc, mrb(b)%ifl(d)%s, mrb(b)%ifl(d)%n, ld_flg, ld_pos, &
                                 ld_p1, nface, ofs_i(d), mrfl(d,l_)%if_, mrfl(d,l_)%ilc, mrfl(d,l_)%isl, max(ni(d),1) )
          else if ( d == 2 ) then
            call ld_gather_ifl ( mrb(b)%ifl(d)%f, mrb(b)%ifl(d)%lc, mrb(b)%ifl(d)%s, mrb(b)%ifl(d)%n, ld_flg, ld_pos, &
                                 ld_p2, nface, ofs_i(d), mrfl(d,l_)%if_, mrfl(d,l_)%ilc, mrfl(d,l_)%isl, max(ni(d),1) )
          else
            call ld_gather_ifl ( mrb(b)%ifl(d)%f, mrb(b)%ifl(d)%lc, mrb(b)%ifl(d)%s, mrb(b)%ifl(d)%n, ld_flg, ld_pos, &
                                 ld_p3, nface, ofs_i(d), mrfl(d,l_)%if_, mrfl(d,l_)%ilc, mrfl(d,l_)%isl, max(ni(d),1) )
          end if
        end if
        ofs_f(d) = ofs_f(d) + cntb(d, b); ofs_i(d) = ofs_i(d) + cntb(4+d, b)
      end do
      ! the cells of the level in this block, with their six face positions
      if ( cntb(4, b) > 0 ) then
        call ld_flag_cells ( mrb(b)%lev, n1, n2, n3, l_, ld_flg, n1*n2*n3 )
        call ld_scan ( ld_flg, ld_pos, n1*n2*n3, ld_csum, ld_coff, size(ld_csum), ld_ssum, size(ld_ssum) )
        call ld_gather_cells ( mrb(b)%lev, domain%blk(b)%vol, domain%blk(b)%m_f, domain%blk(b)%yn, n1, n2, n3, boff, pi, pj, &
                               gc, ofs_c, ld_flg, ld_pos, ld_p1, ld_p2, ld_p3, nmax, &
                               mrcl(l_)%c, mrcl(l_)%pf, mrcl(l_)%ijk, mrcl(l_)%bo, mrcl(l_)%pp, mrcl(l_)%vol, mrcl(l_)%met, &
                               mrcl(l_)%wd, max(ncl,1) )
      end if
      ofs_c = ofs_c + cntb(4, b)
    end do
    !$acc wait(ld_q)
    deallocate ( cntb )
  end subroutine mr_build_lists_dev

  subroutine ld_set_i ( a, n, v )
    implicit none
    integer, intent(in)  :: n, v
    integer, intent(out) :: a(n)
    integer :: e
    !$acc parallel loop gang vector async(ld_q) present(a)
    do e = 1, n
      a(e) = v
    end do
  end subroutine ld_set_i

  subroutine ld_set_r ( a, n, v )
    implicit none
    integer,  intent(in)  :: n
    real(R8), intent(in)  :: v
    real(R8), intent(out) :: a(n)
    integer :: e
    !$acc parallel loop gang vector async(ld_q) present(a)
    do e = 1, n
      a(e) = v
    end do
  end subroutine ld_set_r

  !> The counts of one block for level l_: faces i/j/k (a cell of the level on either side) and cells.
  subroutine ld_count_blk ( lev, n, l_, cnt )
    implicit none
    integer, intent(in)  :: n(3), l_, lev(n(1),n(2),n(3))
    integer, intent(out) :: cnt(7)
    integer :: i, j, k, n1, n2, n3, c1, c2, c3, cc, la
    n1 = n(1); n2 = n(2); n3 = n(3); c1 = 0; c2 = 0; c3 = 0; cc = 0
    !$acc parallel loop collapse(3) gang vector present(lev) reduction(+:c1, c2, c3, cc) private(la)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      la = lev(i,j,k)
      if ( la == l_ ) cc = cc + 1
      if ( i < n1 ) then; if ( la == l_ .or. lev(i+1,j,k) == l_ ) c1 = c1 + 1; end if
      if ( j < n2 ) then; if ( la == l_ .or. lev(i,j+1,k) == l_ ) c2 = c2 + 1; end if
      if ( k < n3 ) then; if ( la == l_ .or. lev(i,j,k+1) == l_ ) c3 = c3 + 1; end if
    end do
    end do
    end do
    cnt(1) = c1; cnt(2) = c2; cnt(3) = c3; cnt(4) = cc; cnt(5:7) = 0
  end subroutine ld_count_blk

  subroutine ld_count_ifl ( lc, n, l_, cnt )
    implicit none
    integer, intent(in)  :: n, l_, lc(n)
    integer, intent(out) :: cnt
    integer :: q, c
    c = 0
    !$acc parallel loop gang vector present(lc) reduction(+:c)
    do q = 1, n
      if ( lc(q) == l_ .or. lc(q) + 1 == l_ ) c = c + 1
    end do
    cnt = c
  end subroutine ld_count_ifl

  !> mr-build-fast: ld_count_blk on the list queue, the four sums of each chunk of LD_CH cells added to cnt(1:4, b) (one
  !> atomic per chunk and sum; integers, so the order of the additions does not matter)
  subroutine ld_count_blk_q ( lev, n, l_, cnt, nb, b )
    implicit none
    integer, intent(in)    :: n(3), l_, nb, b, lev(n(1),n(2),n(3))
    integer, intent(inout) :: cnt(7, nb)
    integer :: i, j, k, n1, n2, n3, c1, c2, c3, cc, la, ncell, nch, ch, e, r
    n1 = n(1); n2 = n(2); n3 = n(3); ncell = n1*n2*n3; nch = (ncell + LD_CH - 1) / LD_CH
    !$acc parallel loop gang async(ld_q) present(lev, cnt) private(c1, c2, c3, cc)
    do ch = 1, nch
      c1 = 0; c2 = 0; c3 = 0; cc = 0
      !$acc loop vector reduction(+:c1, c2, c3, cc) private(i, j, k, r, la)
      do e = (ch-1)*LD_CH + 1, min(ch*LD_CH, ncell)
        r = e - 1; i = mod(r, n1) + 1; r = r / n1; j = mod(r, n2) + 1; k = r / n2 + 1
        la = lev(i,j,k)
        if ( la == l_ ) cc = cc + 1
        if ( i < n1 ) then; if ( la == l_ .or. lev(i+1,j,k) == l_ ) c1 = c1 + 1; end if
        if ( j < n2 ) then; if ( la == l_ .or. lev(i,j+1,k) == l_ ) c2 = c2 + 1; end if
        if ( k < n3 ) then; if ( la == l_ .or. lev(i,j,k+1) == l_ ) c3 = c3 + 1; end if
      end do
      if ( c1 > 0 ) then
        !$acc atomic update
        cnt(1,b) = cnt(1,b) + c1
      end if
      if ( c2 > 0 ) then
        !$acc atomic update
        cnt(2,b) = cnt(2,b) + c2
      end if
      if ( c3 > 0 ) then
        !$acc atomic update
        cnt(3,b) = cnt(3,b) + c3
      end if
      if ( cc > 0 ) then
        !$acc atomic update
        cnt(4,b) = cnt(4,b) + cc
      end if
    end do
  end subroutine ld_count_blk_q

  !> mr-build-fast: ld_count_ifl on the list queue, added to cnt(m, b) per chunk of LD_CH records
  subroutine ld_count_ifl_q ( lc, n, l_, cnt, nb, b, m )
    implicit none
    integer, intent(in)    :: n, l_, nb, b, m, lc(n)
    integer, intent(inout) :: cnt(7, nb)
    integer :: q, c, nch, ch
    nch = (n + LD_CH - 1) / LD_CH
    !$acc parallel loop gang async(ld_q) present(lc, cnt) private(c)
    do ch = 1, nch
      c = 0
      !$acc loop vector reduction(+:c)
      do q = (ch-1)*LD_CH + 1, min(ch*LD_CH, n)
        if ( lc(q) == l_ .or. lc(q) + 1 == l_ ) c = c + 1
      end do
      if ( c > 0 ) then
        !$acc atomic update
        cnt(m,b) = cnt(m,b) + c
      end if
    end do
  end subroutine ld_count_ifl_q

  !> flg(f) = 1 on the faces of direction d (lexicographic index f) with a cell of the level on either side.
  subroutine ld_flag_faces ( lev, n1, n2, n3, d, l_, flg, nface )
    implicit none
    integer, intent(in)  :: n1, n2, n3, d, l_, nface, lev(n1,n2,n3)
    integer, intent(out) :: flg(nface)
    integer :: i, j, k, f, m1, m2, m3, di, dj, dk
    di = 0; dj = 0; dk = 0; m1 = n1; m2 = n2; m3 = n3
    if ( d == 1 ) then; di = 1; m1 = n1 - 1; end if
    if ( d == 2 ) then; dj = 1; m2 = n2 - 1; end if
    if ( d == 3 ) then; dk = 1; m3 = n3 - 1; end if
    !$acc parallel loop collapse(3) gang vector async(ld_q) present(lev, flg) private(f)
    do k = 1, m3
    do j = 1, m2
    do i = 1, m1
      f = (i-1) + (j-1)*m1 + (k-1)*m1*m2 + 1
      if ( lev(i,j,k) == l_ .or. lev(i+di,j+dj,k+dk) == l_ ) then
        flg(f) = 1
      else
        flg(f) = 0
      end if
    end do
    end do
    end do
  end subroutine ld_flag_faces

  subroutine ld_flag_cells ( lev, n1, n2, n3, l_, flg, ncell )
    implicit none
    integer, intent(in)  :: n1, n2, n3, l_, ncell, lev(n1,n2,n3)
    integer, intent(out) :: flg(ncell)
    integer :: i, j, k, f
    !$acc parallel loop collapse(3) gang vector async(ld_q) present(lev, flg) private(f)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      f = (i-1) + (j-1)*n1 + (k-1)*n1*n2 + 1
      if ( lev(i,j,k) == l_ ) then
        flg(f) = 1
      else
        flg(f) = 0
      end if
    end do
    end do
    end do
  end subroutine ld_flag_cells

  subroutine ld_flag_ifl ( lc, n, l_, flg )
    implicit none
    integer, intent(in)  :: n, l_, lc(n)
    integer, intent(out) :: flg(n)
    integer :: q
    !$acc parallel loop gang vector async(ld_q) present(lc, flg)
    do q = 1, n
      if ( lc(q) == l_ .or. lc(q) + 1 == l_ ) then
        flg(q) = 1
      else
        flg(q) = 0
      end if
    end do
  end subroutine ld_flag_ifl

  !> Exclusive prefix sum of flg(1:n) (0/1) into pos(1:n): chunk sums, their sequential scan, the in-chunk prefix.
  !> Integer arithmetic in a fixed order: the same positions on every run.
  !> mr-build-fast: the in-chunk prefix by one thread per chunk of LD_CH waited on memory (126 us a call on the cone, 2.2 s
  !> of device time in 38 macro steps at one GPU, nsys of job 315305); the chunk is cut in LD_NS sub-chunks of LD_SC: their
  !> sums (ssum), the chunk sums from them, the sequential scan of the chunks, then one thread per sub-chunk adds the earlier
  !> sub-chunks of its chunk and writes its LD_SC positions. Integers: the same positions as the chunk prefix.
  subroutine ld_scan ( flg, pos, n, csum, coff, ncs, ssum, nss )
    implicit none
    integer, intent(in)  :: n, ncs, nss
    integer, intent(in)  :: flg(n)
    integer, intent(out) :: pos(n), csum(ncs), coff(ncs), ssum(nss)
    integer :: nch, ch, e, lo, hi, s, run, nsub, q, q0, p
    logical :: fast
    nch = (n + LD_CH - 1) / LD_CH
    fast = mr_build_fast_on()
    nsub = (n + LD_SC - 1) / LD_SC
    if ( fast ) then
      if ( nsub > nss ) error stop 'ld_scan: the sub-chunk sums do not fit'
      !$acc parallel loop gang vector async(ld_q) present(flg, ssum) private(lo, hi, s)
      do q = 1, nsub
        lo = (q-1)*LD_SC + 1; hi = min(q*LD_SC, n); s = 0
        !$acc loop seq
        do e = lo, hi
          s = s + flg(e)
        end do
        ssum(q) = s
      end do
      !$acc parallel loop gang vector async(ld_q) present(ssum, csum) private(lo, hi, s)
      do ch = 1, nch
        lo = (ch-1)*LD_NS + 1; hi = min(ch*LD_NS, nsub); s = 0
        !$acc loop seq
        do p = lo, hi
          s = s + ssum(p)
        end do
        csum(ch) = s
      end do
    else
    !$acc parallel loop gang async(ld_q) present(flg, csum) private(lo, hi, s)
    do ch = 1, nch
      lo = (ch-1)*LD_CH + 1; hi = min(ch*LD_CH, n); s = 0
      !$acc loop vector reduction(+:s)
      do e = lo, hi
        s = s + flg(e)
      end do
      csum(ch) = s
    end do
    end if
    !$acc parallel async(ld_q) present(csum, coff) num_gangs(1) vector_length(1)
    run = 0
    !$acc loop seq
    do ch = 1, nch
      coff(ch) = run
      run = run + csum(ch)
    end do
    !$acc end parallel
    if ( fast ) then
      !$acc parallel loop gang vector async(ld_q) present(flg, pos, ssum, coff) private(ch, q0, lo, hi, run)
      do q = 1, nsub
        ch = (q-1) / LD_NS + 1; q0 = (ch-1)*LD_NS + 1; run = coff(ch)
        !$acc loop seq
        do p = q0, q - 1
          run = run + ssum(p)
        end do
        lo = (q-1)*LD_SC + 1; hi = min(q*LD_SC, n)
        !$acc loop seq
        do e = lo, hi
          pos(e) = run
          run = run + flg(e)
        end do
      end do
      return
    end if
    !$acc parallel loop gang async(ld_q) present(flg, pos, coff) private(lo, hi, run)
    do ch = 1, nch
      lo = (ch-1)*LD_CH + 1; hi = min(ch*LD_CH, n); run = coff(ch)
      !$acc loop seq
      do e = lo, hi
        pos(e) = run
        run = run + flg(e)
      end do
    end do
  end subroutine ld_scan

  !> The entries of the faces of direction d flagged in this block, at ofs + pos + 1 of the level's list, and the block's
  !> face -> position table (0 where the face is not in the list).
  subroutine ld_gather_faces ( lev, dl_f, fmd, f1, f2, f3, m_f, yn, n1, n2, n3, d, boff, pi, pj, gc, ofs, flg, pos, ptab, nface, &
                               c, sa, sb, sc, dl, fm, mm, ynl, nl )
    implicit none
    integer,  intent(in)  :: n1, n2, n3, d, boff, pi, pj, gc, ofs, nface, nl, f1, f2, f3
    integer,  intent(in)  :: lev(n1,n2,n3), flg(nface), pos(nface)
    real(R8), intent(in)  :: dl_f(3, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    real(R8), intent(in)  :: fmd(4, f1:n1, f2:n2, f3:n3)   ! fm1 (f1 = 0) / fm2 (f2 = 0) / fm3 (f3 = 0)
    real(R8), intent(in)  :: m_f(3,3,1-gc:n1+gc,1-gc:n2+gc,1-gc:n3+gc), yn(1-gc:n1+gc,1-gc:n2+gc,1-gc:n3+gc)
    integer,  intent(out) :: ptab(nface)
    integer,  intent(out) :: c(nl), sa(nl), sb(nl), sc(nl)
    real(R8), intent(out) :: dl(4,nl), fm(4,nl), mm(3,3,2,nl), ynl(2,nl)
    integer :: i, j, k, f, e, m1, m2, m3, di, dj, dk, a1, a2
    di = 0; dj = 0; dk = 0; m1 = n1; m2 = n2; m3 = n3
    if ( d == 1 ) then; di = 1; m1 = n1 - 1; end if
    if ( d == 2 ) then; dj = 1; m2 = n2 - 1; end if
    if ( d == 3 ) then; dk = 1; m3 = n3 - 1; end if
    !$acc parallel loop collapse(3) gang vector async(ld_q) present(lev, dl_f, fmd, m_f, yn, flg, pos, ptab, c, sa, sb, sc, dl, fm, mm, ynl) &
    !$acc&  private(f, e, a1, a2)
    do k = 1, m3
    do j = 1, m2
    do i = 1, m1
      f = (i-1) + (j-1)*m1 + (k-1)*m1*m2 + 1
      if ( flg(f) == 1 ) then
        e = ofs + pos(f) + 1
        ptab(f) = e
        c(e) = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
        if ( d == 1 ) then
          sa(e) = 1; sb(e) = pi; sc(e) = pi*pj
          dl(1,e) = dl_f(1,i-1,j,k); dl(2,e) = dl_f(1,i,j,k); dl(3,e) = dl_f(1,i+1,j,k); dl(4,e) = dl_f(1,i+2,j,k)
        else if ( d == 2 ) then
          sa(e) = pi; sb(e) = 1; sc(e) = pi*pj
          dl(1,e) = dl_f(2,i,j-1,k); dl(2,e) = dl_f(2,i,j,k); dl(3,e) = dl_f(2,i,j+1,k); dl(4,e) = dl_f(2,i,j+2,k)
        else
          sa(e) = pi*pj; sb(e) = 1; sc(e) = pi
          dl(1,e) = dl_f(3,i,j,k-1); dl(2,e) = dl_f(3,i,j,k); dl(3,e) = dl_f(3,i,j,k+1); dl(4,e) = dl_f(3,i,j,k+2)
        end if
        fm(1,e) = fmd(1,i,j,k); fm(2,e) = fmd(2,i,j,k); fm(3,e) = fmd(3,i,j,k); fm(4,e) = fmd(4,i,j,k)
        do a2 = 1, 3
          do a1 = 1, 3
            mm(a1,a2,1,e) = m_f(a1,a2,i,j,k)
            mm(a1,a2,2,e) = m_f(a1,a2,i+di,j+dj,k+dk)
          end do
        end do
        ynl(1,e) = yn(i,j,k); ynl(2,e) = yn(i+di,j+dj,k+dk)
      else
        ptab(f) = 0
      end if
    end do
    end do
    end do
  end subroutine ld_gather_faces

  subroutine ld_gather_ifl ( f, lc, sl, n, flg, pos, ptab, nface, ofs, if_, ilc, isl, nl )
    implicit none
    integer, intent(in)  :: n, nface, ofs, nl, f(n), lc(n), sl(n), flg(n), pos(n), ptab(nface)
    integer, intent(out) :: if_(nl), ilc(nl), isl(nl)
    integer :: q, e
    !$acc parallel loop gang vector async(ld_q) present(f, lc, sl, flg, pos, ptab, if_, ilc, isl) private(e)
    do q = 1, n
      if ( flg(q) == 1 ) then
        e = ofs + pos(q) + 1
        if_(e) = ptab(f(q)); ilc(e) = lc(q); isl(e) = sl(q)
      end if
    end do
  end subroutine ld_gather_ifl

  !> The cells of the level in this block at ofs + pos + 1, with the positions of their six faces from the face tables.
  subroutine ld_gather_cells ( lev, vol, m_f, yn, n1, n2, n3, boff, pi, pj, gc, ofs, flg, pos, p1, p2, p3, nmax, &
                               c, pf, ijk, bo, pp, vl, met, wd, nl )
    implicit none
    integer,  intent(in)  :: n1, n2, n3, boff, pi, pj, gc, ofs, nmax, nl
    integer,  intent(in)  :: lev(n1,n2,n3), flg(n1*n2*n3), pos(n1*n2*n3), p1(nmax), p2(nmax), p3(nmax)
    real(R8), intent(in)  :: vol(1-gc:n1+gc,1-gc:n2+gc,1-gc:n3+gc), m_f(3,3,1-gc:n1+gc,1-gc:n2+gc,1-gc:n3+gc)
    real(R8), intent(in)  :: yn(1-gc:n1+gc,1-gc:n2+gc,1-gc:n3+gc)
    integer,  intent(out) :: c(nl), pf(6,nl), ijk(3,nl), bo(nl), pp(2,nl)
    real(R8), intent(out) :: vl(nl), met(3,3,nl), wd(nl)
    integer :: i, j, k, f, e, a1, a2
    !$acc parallel loop collapse(3) gang vector async(ld_q) present(lev, vol, m_f, yn, flg, pos, p1, p2, p3, c, pf, ijk, bo, pp, vl, met, wd) &
    !$acc&  private(f, e, a1, a2)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      f = (i-1) + (j-1)*n1 + (k-1)*n1*n2 + 1
      if ( flg(f) == 1 ) then
        e = ofs + pos(f) + 1
        c(e) = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
        ijk(1,e) = i; ijk(2,e) = j; ijk(3,e) = k; bo(e) = boff; pp(1,e) = pi; pp(2,e) = pj
        vl(e) = vol(i,j,k); wd(e) = yn(i,j,k)
        do a2 = 1, 3
          do a1 = 1, 3
            met(a1,a2,e) = m_f(a1,a2,i,j,k)
          end do
        end do
        pf(1,e) = 0; pf(2,e) = 0; pf(3,e) = 0; pf(4,e) = 0; pf(5,e) = 0; pf(6,e) = 0
        if ( i < n1 ) pf(1,e) = p1( (i-1) + (j-1)*(n1-1) + (k-1)*(n1-1)*n2 + 1 )
        if ( i > 1  ) pf(2,e) = p1( (i-2) + (j-1)*(n1-1) + (k-1)*(n1-1)*n2 + 1 )
        if ( j < n2 ) pf(3,e) = p2( (i-1) + (j-1)*n1 + (k-1)*n1*(n2-1) + 1 )
        if ( j > 1  ) pf(4,e) = p2( (i-1) + (j-2)*n1 + (k-1)*n1*(n2-1) + 1 )
        if ( k < n3 ) pf(5,e) = p3( (i-1) + (j-1)*n1 + (k-1)*n1*n2 + 1 )
        if ( k > 1  ) pf(6,e) = p3( (i-1) + (j-1)*n1 + (k-2)*n1*n2 + 1 )
      end if
    end do
    end do
    end do
  end subroutine ld_gather_cells

  !> mr-thermo-exact: MOSE_MR_THERMO_EXACT (default on) builds the thermo-cache lists of a list level from the cells read.
  logical function mr_thermo_exact_on ()
    use MOSE_Mod_MPI, only: mpi_is_root
    logical, save :: checked = .false., on = .true.
    character(len=8) :: env
    integer :: st, ln
    if ( .not. checked ) then
      checked = .true.
      call get_environment_variable ( 'MOSE_MR_THERMO_EXACT', env, ln, st )
      if ( st == 0 .and. ln > 0 ) then
        if ( env(1:1) == '0' .or. env(1:1) == 'n' .or. env(1:1) == 'N' .or. env(1:1) == 'f' .or. env(1:1) == 'F' ) on = .false.
      end if
      if ( mpi_is_root ) then
        if ( on ) then
          write(*,'(A)') ' [MR]   thermo cache of the list levels on the cells the stage reads (MOSE_MR_THERMO_EXACT=0: the boxes)'
        else
          write(*,'(A)') ' [MR]   thermo cache of the list levels on the boxes grown by gc'
        end if
      end if
    end if
    mr_thermo_exact_on = on
  end function mr_thermo_exact_on

  !> The thermo-cache lists of list level l_ (Mod_Thermo_Cache Refresh_Thermo_Cache: ca = rho/Rg/T on ghosts and interior,
  !> cb = the interior part, rho/Rg/T + viscosity with the guard) built on the device from the cells the stage's kernels
  !> read in the cache: the ten cells of every face of the level's face lists (the two cells and their four transverse
  !> neighbours each: the diffusive list kernel reads rho and T there), the level's cells (turbulence sources and the T seed
  !> of the update), and the interior cell of every ghost record of the level (the wall ghost of the turbulence reads the
  !> viscosity there). Every cached value is a function of the cell's primitives alone and Compute_dt refreshes every
  !> interior cell at each macro step, so a cell the stage reads holds the same doubles as with the boxes (a superset):
  !> the same result. Order: ascending flat index (the order of the box build). Needs the level's face, cell and ghost
  !> record lists on the device (both list builds put them there).
  subroutine mr_build_thermo_list_dev ( domain, l_ )
    use MOSE_Types_m
    use MOSE_Global_m, only: gc
    use MOSE_Mod_MPI,  only: is_local_block, mpi_is_root
    use MOSE_Field_m,  only: ncell_total, field_block_index_params
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: l_
    integer :: d, b, nct, na, nb, boff, pi, pj
    nct = ncell_total
    if ( allocated(lt_mark) ) then
      if ( size(lt_mark) /= nct ) then
        !$acc exit data delete(lt_mark, lt_int, lt_fa, lt_pa, lt_fb, lt_pb, lt_cs, lt_co, lt_ss)
        deallocate ( lt_mark, lt_int, lt_fa, lt_pa, lt_fb, lt_pb, lt_cs, lt_co, lt_ss )
      end if
    end if
    if ( .not. allocated(lt_mark) ) then
      allocate ( lt_mark(nct), lt_int(nct), lt_fa(nct), lt_pa(nct), lt_fb(nct), lt_pb(nct), lt_cs(nct/LD_CH+2), lt_co(nct/LD_CH+2), &
                 lt_ss(nct/LD_SC+2) )
      !$acc enter data create(lt_mark, lt_int, lt_fa, lt_pa, lt_fb, lt_pb, lt_cs, lt_co, lt_ss)
      call lt_fill ( lt_int, nct, 0 )
      do b = 1, domain % nb   ! the interior flag of the flat field, once
        if (.not. is_local_block(b)) cycle
        call field_block_index_params ( b, boff, pi, pj )
        call lt_interior_blk ( lt_int, nct, boff, pi, pj, domain%blk(b)%dim, gc )
      end do
    end if
    call lt_fill ( lt_mark, nct, 0 )
    do d = 1, 3
      if ( mrfl(d,l_)%n > 0 ) &
        call lt_mark_faces ( mrfl(d,l_)%c, mrfl(d,l_)%sa, mrfl(d,l_)%sb, mrfl(d,l_)%sc, mrfl(d,l_)%n, lt_mark, nct )
    end do
    if ( mrcl(l_)%n > 0 ) call lt_mark_cells ( mrcl(l_)%c, mrcl(l_)%n, lt_mark, nct )
    if ( mrgl(l_)%np > 0 ) &
      call lt_mark_records ( mrgl(l_)%pr, mrgl(l_)%pb, mrgl(l_)%np, domain%bc_bi, domain%bc_bj, domain%bc_bk, &
                             domain%nbound, gc, lt_mark, nct )
    call lt_flags ( lt_mark, lt_int, lt_fa, lt_fb, nct, na, nb )
    mrtl(l_)%na = na; mrtl(l_)%nb = nb
    allocate ( mrtl(l_)%ca(max(na,1)), mrtl(l_)%cb(max(nb,1)) )
    !$acc enter data create(mrtl(l_)%ca, mrtl(l_)%cb)
    if ( na == 0 ) call ld_set_i ( mrtl(l_)%ca, 1, 1 )
    if ( nb == 0 ) call ld_set_i ( mrtl(l_)%cb, 1, 1 )
    if ( na > 0 ) then
      call ld_scan ( lt_fa, lt_pa, nct, lt_cs, lt_co, size(lt_cs), lt_ss, size(lt_ss) )
      call lt_gather ( lt_fa, lt_pa, nct, mrtl(l_)%ca, na )
    end if
    if ( nb > 0 ) then
      call ld_scan ( lt_fb, lt_pb, nct, lt_cs, lt_co, size(lt_cs), lt_ss, size(lt_ss) )
      call lt_gather ( lt_fb, lt_pb, nct, mrtl(l_)%cb, nb )
    end if
    !$acc wait(ld_q)
    if ( mpi_is_root .and. obj_mr%nfull == 1 ) write(*,'(A,I0,A,I0,A,I0,A)') ' [MR]   thermo list level ', l_, &
      ': ', na, ' cells (', nb, ' interior), the cells the stage reads'
  end subroutine mr_build_thermo_list_dev

  subroutine lt_fill ( a, n, v )
    implicit none
    integer, intent(in)  :: n, v
    integer, intent(out) :: a(n)
    integer :: e
    !$acc parallel loop gang vector present(a)
    do e = 1, n
      a(e) = v
    end do
  end subroutine lt_fill

  subroutine lt_interior_blk ( fint, nct, boff, pi, pj, n, gc )
    implicit none
    integer, intent(in)    :: nct, boff, pi, pj, n(3), gc
    integer, intent(inout) :: fint(nct)
    integer :: i, j, k, n1, n2, n3
    n1 = n(1); n2 = n(2); n3 = n(3)
    !$acc parallel loop collapse(3) gang vector present(fint)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      fint( boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1 ) = 1
    end do
    end do
    end do
  end subroutine lt_interior_blk

  !> The ten cells of every face of the list (the cells the diffusive list kernel reads in the cache).
  subroutine lt_mark_faces ( c, sa, sb, sc, n, mark, nct )
    implicit none
    integer, intent(in)    :: n, nct, c(n), sa(n), sb(n), sc(n)
    integer, intent(inout) :: mark(nct)
    integer :: e, c1, c2
    !$acc parallel loop gang vector present(c, sa, sb, sc, mark) private(c1, c2)
    do e = 1, n
      c1 = c(e); c2 = c1 + sa(e)
      mark(c1) = 1; mark(c2) = 1
      mark(c1 - sb(e)) = 1; mark(c1 + sb(e)) = 1; mark(c2 - sb(e)) = 1; mark(c2 + sb(e)) = 1
      mark(c1 - sc(e)) = 1; mark(c1 + sc(e)) = 1; mark(c2 - sc(e)) = 1; mark(c2 + sc(e)) = 1
    end do
  end subroutine lt_mark_faces

  subroutine lt_mark_cells ( c, n, mark, nct )
    implicit none
    integer, intent(in)    :: n, nct, c(n)
    integer, intent(inout) :: mark(nct)
    integer :: e
    !$acc parallel loop gang vector present(c, mark)
    do e = 1, n
      mark(c(e)) = 1
    end do
  end subroutine lt_mark_cells

  !> The interior cell of every ghost record of the level (the wall ghost of the turbulence reads the viscosity there).
  subroutine lt_mark_records ( pr, pb, n, bc_bi, bc_bj, bc_bk, nbnd, gc, mark, nct )
    use MOSE_Field_m, only: dboff, dpi, dpj
    implicit none
    integer, intent(in)    :: n, nbnd, gc, nct, pr(n), pb(n), bc_bi(nbnd), bc_bj(nbnd), bc_bk(nbnd)
    integer, intent(inout) :: mark(nct)
    integer :: e, idx, b
    !$acc parallel loop gang vector present(pr, pb, bc_bi, bc_bj, bc_bk, dboff, dpi, dpj, mark) private(idx, b)
    do e = 1, n
      idx = pr(e); b = pb(e)
      mark( dboff(b) + (bc_bi(idx)+gc-1) + (bc_bj(idx)+gc-1)*dpi(b) + (bc_bk(idx)+gc-1)*dpi(b)*dpj(b) + 1 ) = 1
    end do
  end subroutine lt_mark_records

  !> fa = the marked cells, fb = the marked interior cells; their counts.
  subroutine lt_flags ( mark, fint, fa, fb, nct, na, nb )
    implicit none
    integer, intent(in)  :: nct, mark(nct), fint(nct)
    integer, intent(out) :: fa(nct), fb(nct), na, nb
    integer :: e, sa_, sb_
    sa_ = 0; sb_ = 0
    !$acc parallel loop gang vector present(mark, fint, fa, fb) reduction(+:sa_, sb_)
    do e = 1, nct
      fa(e) = mark(e); fb(e) = mark(e) * fint(e)
      sa_ = sa_ + fa(e); sb_ = sb_ + fb(e)
    end do
    na = sa_; nb = sb_
  end subroutine lt_flags

  subroutine lt_gather ( flg, pos, nct, lst, nl )
    implicit none
    integer, intent(in)  :: nct, nl, flg(nct), pos(nct)
    integer, intent(out) :: lst(nl)
    integer :: e
    !$acc parallel loop gang vector async(ld_q) present(flg, pos, lst)
    do e = 1, nct
      if ( flg(e) == 1 ) lst( pos(e) + 1 ) = e
    end do
  end subroutine lt_gather

  !> mr-ghost-exact: MOSE_MR_GHOST_EXACT (default on) chooses the ghost records of a list level from the cells read.
  logical function mr_ghost_exact_on ()
    use MOSE_Mod_MPI, only: mpi_is_root, mpi_size_
    logical, save :: checked = .false., on = .true.
    character(len=8) :: env
    integer :: st, ln
    if ( .not. checked ) then
      checked = .true.
      call get_environment_variable ( 'MOSE_MR_GHOST_EXACT', env, ln, st )
      if ( st == 0 .and. ln > 0 ) then
        if ( env(1:1) == '0' .or. env(1:1) == 'n' .or. env(1:1) == 'N' .or. env(1:1) == 'f' .or. env(1:1) == 'F' ) on = .false.
      end if
      ! one rank only for now: with more ranks the transverse pack of the exchange reads donor ghosts chosen by the
      ! receiver's stencils, which the read mark of the donor does not hold (the boxes cover them, gated)
      if ( mpi_size_ > 1 ) on = .false.
      if ( mpi_is_root ) then
        if ( on ) then
          write(*,'(A)') ' [MR]   ghost records of the list levels from the cells the stage reads (MOSE_MR_GHOST_EXACT=0: the boxes)'
        else
          write(*,'(A)') ' [MR]   ghost records of the list levels from the boxes grown by gc'
        end if
      end if
    end if
    mr_ghost_exact_on = on
  end function mr_ghost_exact_on

  !> The cells whose primitives a stage of list level l_ reads after its ghost fill, marked on the device and copied to
  !> the host for mr_build_ghost_list: the convective stencil of every face (two cells each side along the normal), the ten
  !> cells of the diffusive stencil, the level's cells and their six neighbours (turbulence sources, update), and around
  !> every boundary record of the level the line of gc cells each side along the normal with the transverse neighbours
  !> (a box that holds what the boundary flux functions read). The thermo list (the cells read in the cache) is a subset.
  !> A ghost value is a function of the current interior state and boundary data alone, so filling exactly the ghosts the
  !> stage reads gives the same doubles as the boxes. The connection transverse cells stay marked by mr_build_ghost_list.
  subroutine mr_mark_reads_dev ( domain, l_ )
    use MOSE_Types_m
    use MOSE_Global_m, only: gc
    use MOSE_Field_m,  only: ncell_total
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: l_
    integer :: d, nct
    nct = ncell_total
    if ( allocated(rk_mark) ) then
      if ( size(rk_mark) /= nct ) then
        !$acc exit data delete(rk_mark)
        deallocate ( rk_mark )
      end if
    end if
    if ( .not. allocated(rk_mark) ) then
      allocate ( rk_mark(nct) )
      !$acc enter data create(rk_mark)
    end if
    call lt_fill ( rk_mark, nct, 0 )
    do d = 1, 3
      if ( mrfl(d,l_)%n > 0 ) then
        call rk_mark_conv ( mrfl(d,l_)%c, mrfl(d,l_)%sa, mrfl(d,l_)%n, rk_mark, nct )
        call lt_mark_faces ( mrfl(d,l_)%c, mrfl(d,l_)%sa, mrfl(d,l_)%sb, mrfl(d,l_)%sc, mrfl(d,l_)%n, rk_mark, nct )
      end if
    end do
    if ( mrcl(l_)%n > 0 ) call rk_mark_cells6 ( mrcl(l_)%c, mrcl(l_)%pp, mrcl(l_)%n, rk_mark, nct )
    do d = 1, 6
      if ( mrbl(d,l_)%n > 0 ) &
        call rk_mark_records ( mrbl(d,l_)%rec, mrbl(d,l_)%rb, mrbl(d,l_)%n, domain%bc_bi, domain%bc_bj, domain%bc_bk, &
                               domain%bc_bf, domain%nbound, gc, rk_mark, nct )
    end do
    !$acc update host(rk_mark)
    rk_ready = .true.
  end subroutine mr_mark_reads_dev

  !> The convective stencil of every face of the list: c - sa, c, c + sa, c + 2 sa.
  subroutine rk_mark_conv ( c, sa, n, mark, nct )
    implicit none
    integer, intent(in)    :: n, nct, c(n), sa(n)
    integer, intent(inout) :: mark(nct)
    integer :: e
    !$acc parallel loop gang vector present(c, sa, mark)
    do e = 1, n
      mark(c(e) - sa(e)) = 1; mark(c(e)) = 1; mark(c(e) + sa(e)) = 1; mark(c(e) + 2*sa(e)) = 1
    end do
  end subroutine rk_mark_conv

  !> The level's cells and their six neighbours (pp = the row lengths of the cell's block).
  subroutine rk_mark_cells6 ( c, pp, n, mark, nct )
    implicit none
    integer, intent(in)    :: n, nct, c(n), pp(2, n)
    integer, intent(inout) :: mark(nct)
    integer :: e, c0, si, sj
    !$acc parallel loop gang vector present(c, pp, mark) private(c0, si, sj)
    do e = 1, n
      c0 = c(e); si = pp(1,e); sj = pp(1,e)*pp(2,e)
      mark(c0) = 1
      mark(c0 - 1) = 1; mark(c0 + 1) = 1; mark(c0 - si) = 1; mark(c0 + si) = 1; mark(c0 - sj) = 1; mark(c0 + sj) = 1
    end do
  end subroutine rk_mark_cells6

  !> Around every boundary record: the cells at -gc..gc along the face normal from the record's cell, each with its
  !> transverse neighbours (the 3 x 3 cross-section), inside the padded block.
  subroutine rk_mark_records ( rec, rb, n, bc_bi, bc_bj, bc_bk, bc_bf, nbnd, gc, mark, nct )
    use MOSE_Field_m, only: dboff, dpi, dpj
    implicit none
    integer, intent(in)    :: n, nbnd, gc, nct, rec(n), rb(n), bc_bi(nbnd), bc_bj(nbnd), bc_bk(nbnd), bc_bf(nbnd)
    integer, intent(inout) :: mark(nct)
    integer :: e, idx, b, c0, sn, s1, s2, t, a1, a2, dir, pi, pj
    !$acc parallel loop gang vector present(rec, rb, bc_bi, bc_bj, bc_bk, bc_bf, dboff, dpi, dpj, mark) &
    !$acc&  private(idx, b, c0, sn, s1, s2, t, a1, a2, dir, pi, pj)
    do e = 1, n
      idx = rec(e); b = rb(e); pi = dpi(b); pj = dpj(b)
      c0 = dboff(b) + (bc_bi(idx)+gc-1) + (bc_bj(idx)+gc-1)*pi + (bc_bk(idx)+gc-1)*pi*pj + 1
      dir = (bc_bf(idx) + 1) / 2
      if ( dir == 1 ) then
        sn = 1; s1 = pi; s2 = pi*pj
      else if ( dir == 2 ) then
        sn = pi; s1 = 1; s2 = pi*pj
      else
        sn = pi*pj; s1 = 1; s2 = pi
      end if
      !$acc loop seq
      do t = -gc, gc
        !$acc loop seq
        do a2 = -1, 1
          !$acc loop seq
          do a1 = -1, 1
            mark( c0 + t*sn + a1*s1 + a2*s2 ) = 1
          end do
        end do
      end do
    end do
  end subroutine rk_mark_records

  subroutine mr_free_lists ( l_ )
    implicit none
    integer, intent(in) :: l_
    integer :: d
    if ( .not. allocated(mrfl) ) return
    do d = 1, 3
      if ( allocated(mrfl(d,l_)%c) ) then
        !$acc exit data delete(mrfl(d,l_)%c, mrfl(d,l_)%sa, mrfl(d,l_)%sb, mrfl(d,l_)%sc, mrfl(d,l_)%dl, mrfl(d,l_)%fm) &
        !$acc&          delete(mrfl(d,l_)%mm, mrfl(d,l_)%yn, mrfl(d,l_)%if_, mrfl(d,l_)%ilc, mrfl(d,l_)%isl) &
        !$acc&          delete(mrfl(d,l_)%fh, mrfl(d,l_)%dfh, mrfl(d,l_)%tf)
        deallocate ( mrfl(d,l_)%c, mrfl(d,l_)%sa, mrfl(d,l_)%sb, mrfl(d,l_)%sc, mrfl(d,l_)%dl, mrfl(d,l_)%fm, mrfl(d,l_)%mm, &
                     mrfl(d,l_)%yn, mrfl(d,l_)%fh, mrfl(d,l_)%dfh, mrfl(d,l_)%tf, mrfl(d,l_)%if_, mrfl(d,l_)%ilc, &
                     mrfl(d,l_)%isl )
      end if
      mrfl(d,l_)%n = 0; mrfl(d,l_)%ni = 0
    end do
    do d = 1, 6
      if ( allocated(mrbl(d,l_)%rec) ) then
        !$acc exit data delete(mrbl(d,l_)%rec, mrbl(d,l_)%rb, mrbl(d,l_)%gn, mrbl(d,l_)%gdl, mrbl(d,l_)%gm, mrbl(d,l_)%gy)
        deallocate ( mrbl(d,l_)%rec, mrbl(d,l_)%rb, mrbl(d,l_)%gn, mrbl(d,l_)%gdl, mrbl(d,l_)%gm, mrbl(d,l_)%gy )
      end if
      mrbl(d,l_)%n = 0
    end do
    if ( allocated(mrtl(l_)%ca) ) then
      !$acc exit data delete(mrtl(l_)%ca, mrtl(l_)%cb)
      deallocate ( mrtl(l_)%ca, mrtl(l_)%cb )
    end if
    mrtl(l_)%na = 0; mrtl(l_)%nb = 0
    if ( allocated(mrgl(l_)%pr) ) then
      !$acc exit data delete(mrgl(l_)%pr, mrgl(l_)%pb, mrgl(l_)%pn, mrgl(l_)%pyn, mrgl(l_)%cr, mrgl(l_)%cm, mrgl(l_)%cs, &
      !$acc&                 mrgl(l_)%pbr, mrgl(l_)%pbs)
      deallocate ( mrgl(l_)%pr, mrgl(l_)%pb, mrgl(l_)%pn, mrgl(l_)%pyn, mrgl(l_)%cr, mrgl(l_)%cm, mrgl(l_)%cs, &
                   mrgl(l_)%sp, mrgl(l_)%pbr, mrgl(l_)%pbs )
    end if
    mrgl(l_)%np = 0; mrgl(l_)%nc = 0; mrgl(l_)%nself = 0
    if ( allocated(mrgl(l_)%qk) ) then
      !$acc exit data delete(mrgl(l_)%qk)
      deallocate ( mrgl(l_)%qk )
    end if
    mrgl(l_)%nq = -1
    if ( allocated(mrcl(l_)%c) ) then
      !$acc exit data delete(mrcl(l_)%c, mrcl(l_)%pf, mrcl(l_)%ijk, mrcl(l_)%bo, mrcl(l_)%pp, mrcl(l_)%vol, mrcl(l_)%met, &
      !$acc&                mrcl(l_)%wd)
      deallocate ( mrcl(l_)%c, mrcl(l_)%pf, mrcl(l_)%ijk, mrcl(l_)%bo, mrcl(l_)%pp, mrcl(l_)%vol, mrcl(l_)%met, mrcl(l_)%wd )
    end if
    mrcl(l_)%n = 0
  end subroutine mr_free_lists

  !> The ghost-fill lists of level l_ (Lib_Multirate_Base mr_gl_type): every physical record of the blocks with cells of
  !> the level -- the records the per-block fill of a stage of the level writes --, with the face normal the per-block
  !> kernel reads from fm<dir>; the connection records of the pairs into those blocks from another block, in pair order;
  !> the pairs of a block with itself apart.
  subroutine mr_build_ghost_list ( domain, l_ )
    use MOSE_Types_m
    use MOSE_Global_m,     only: gc
    use MOSE_Parameters_m, only: guide
    use MOSE_Mod_MPI,      only: is_local_block
    use MOSE_Field_m,      only: ncell_total, field_block_index_params
    use MOSE_Mod_BC_Q2D,   only: q2d_n, q2d_idx, q2d_blk
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: l_
    integer :: b, p, idx, c, np_, nc_, ns_, nb_, Im, Jm, Km, Fm, dir, modfm, Fi, Fj, Fk, pass, Bs, g, ib, i, j, k, q
    integer :: n(3), lo(3), hi(3), boff, pi, pj, nq(3,4)
    logical :: near, keep
    logical, allocatable :: nk(:)
    ! the self pairs (their own launch over the whole pair) and, with them, the gather over every record: no selection
    ns_ = 0
    do p = 1, domain % n_conn_pair
      b = domain % conn_pair_bm(p)
      if ( mrb(b)%nbox(l_) > 0 .and. domain % conn_pair_bs(p) == b ) ns_ = ns_ + 1
    end do
    near = obj_mr%ghost_near .and. ns_ == 0
    ! mr-ghost-near: the ghosts the stage reads, on the flat rows -- the boxes of the level grown by gc (the cells the thermo
    ! cache, the face stencils, the boundary records and the SST source of the level read), and the donor ghosts the
    ! transverse gather (Ghost_Connection_pairs_B) of the connection records with a cell of the level reads
    if ( near ) then
      allocate ( nk(ncell_total) ); nk = .false.
      if ( rk_ready ) then   ! mr-ghost-exact: the cells whose primitives the stage reads (mr_mark_reads_dev)
        nk = ( rk_mark(1:ncell_total) /= 0 )
      else
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        if ( mrb(b)%nbox(l_) <= 0 ) cycle
        n = domain%blk(b)%dim
        call field_block_index_params ( b, boff, pi, pj )
        do ib = 1, mrb(b)%nbox(l_)
          lo = (/ mrb(b)%box(1,ib,l_), mrb(b)%box(3,ib,l_), mrb(b)%box(5,ib,l_) /)
          hi = (/ mrb(b)%box(2,ib,l_), mrb(b)%box(4,ib,l_), mrb(b)%box(6,ib,l_) /)
          lo = max( lo - gc, 1 - gc ); hi = min( hi + gc, n + gc )
          do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
            nk( boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1 ) = .true.
          end do; end do; end do
        end do
      end do
      end if   ! mr-ghost-exact
      do p = 1, domain % n_conn_pair
        b = domain % conn_pair_bm(p); Bs = domain % conn_pair_bs(p)
        if ( mrb(b)%nbox(l_) <= 0 .or. Bs == b ) cycle
        if (.not. is_local_block(Bs)) cycle
        n = domain%blk(Bs)%dim
        call field_block_index_params ( Bs, boff, pi, pj )
        do c = domain % conn_pair_lo(p), domain % conn_pair_hi(p)
          if ( mrb(b)%lev(domain%conn_im(c), domain%conn_jm(c), domain%conn_km(c)) /= l_ ) cycle
          call mr_conn_tnb ( domain, c, nq )
          do q = 1, 4
            if ( all(nq(:,q) >= 1 .and. nq(:,q) <= n) ) cycle   ! interior: always current
            nk( boff + (nq(1,q)+gc-1) + (nq(2,q)+gc-1)*pi + (nq(3,q)+gc-1)*pi*pj + 1 ) = .true.
          end do
        end do
      end do
    end if
    do pass = 1, 2
      np_ = 0; nc_ = 0; nb_ = 0
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        if ( mrb(b)%nbox(l_) <= 0 ) cycle
        call field_block_index_params ( b, boff, pi, pj )
        do idx = domain % bc_blk_lo(b), domain % bc_blk_hi(b)
          Im = domain%bc_bi(idx); Jm = domain%bc_bj(idx); Km = domain%bc_bk(idx); Fm = domain%bc_bf(idx)
          keep = .not. near
          do g = 1, gc
            if ( keep ) exit
            keep = nk( boff + (Im-g*guide(Fm,1)+gc-1) + (Jm-g*guide(Fm,2)+gc-1)*pi + (Km-g*guide(Fm,3)+gc-1)*pi*pj + 1 )
          end do
          if ( .not. keep ) cycle
          np_ = np_ + 1
          if ( pass == 1 ) cycle
          mrgl(l_)%pr(np_) = idx; mrgl(l_)%pb(np_) = b
          mrgl(l_)%pyn(np_) = domain%blk(b)%yn(Im,Jm,Km)
          dir = (Fm + 1)/2; modfm = mod(Fm, 2)
          Fi = Im; Fj = Jm; Fk = Km
          if ( dir == 1 ) then
            Fi = Im - modfm; mrgl(l_)%pn(1:3,np_) = domain%blk(b)%fm1(1:3,Fi,Fj,Fk)
          else if ( dir == 2 ) then
            Fj = Jm - modfm; mrgl(l_)%pn(1:3,np_) = domain%blk(b)%fm2(1:3,Fi,Fj,Fk)
          else
            Fk = Km - modfm; mrgl(l_)%pn(1:3,np_) = domain%blk(b)%fm3(1:3,Fi,Fj,Fk)
          end if
        end do
      end do
      ns_ = 0
      do p = 1, domain % n_conn_pair
        b = domain % conn_pair_bm(p)
        if ( mrb(b)%nbox(l_) <= 0 ) cycle
        if ( domain % conn_pair_bs(p) == b ) then
          ns_ = ns_ + 1
          if ( pass == 2 ) mrgl(l_)%sp(ns_) = p
          cycle
        end if
        call field_block_index_params ( b, boff, pi, pj )
        do c = domain % conn_pair_lo(p), domain % conn_pair_hi(p)
          Im = domain%conn_im(c); Jm = domain%conn_jm(c); Km = domain%conn_km(c); Fm = domain%conn_fm(c)
          keep = .not. near
          do g = 1, gc
            if ( keep ) exit
            keep = nk( boff + (Im-g*guide(Fm,1)+gc-1) + (Jm-g*guide(Fm,2)+gc-1)*pi + (Km-g*guide(Fm,3)+gc-1)*pi*pj + 1 )
          end do
          if ( keep ) then
            nc_ = nc_ + 1
            if ( pass == 2 ) then
              mrgl(l_)%cr(nc_) = c; mrgl(l_)%cm(nc_) = b; mrgl(l_)%cs(nc_) = domain % conn_pair_bs(p)
            end if
          end if
          ! the transverse gather: every record of the pair, or (mr-ghost-near) the records with a cell of the level
          if ( near ) then
            if ( mrb(b)%lev(Im,Jm,Km) /= l_ ) cycle
          else if ( .not. keep ) then
            cycle
          end if
          nb_ = nb_ + 1
          if ( pass == 2 ) then
            mrgl(l_)%pbr(nb_) = c; mrgl(l_)%pbs(nb_) = domain % conn_pair_bs(p)
          end if
        end do
      end do
      if ( pass == 1 ) then
        mrgl(l_)%np = np_; mrgl(l_)%nc = nc_; mrgl(l_)%nself = ns_; mrgl(l_)%npb = nb_
        allocate ( mrgl(l_)%pr(max(np_,1)), mrgl(l_)%pb(max(np_,1)), mrgl(l_)%pn(3,max(np_,1)), mrgl(l_)%pyn(max(np_,1)), &
                   mrgl(l_)%cr(max(nc_,1)), mrgl(l_)%cm(max(nc_,1)), mrgl(l_)%cs(max(nc_,1)), mrgl(l_)%sp(max(ns_,1)), &
                   mrgl(l_)%pbr(max(nb_,1)), mrgl(l_)%pbs(max(nb_,1)) )
        mrgl(l_)%pr = 1; mrgl(l_)%pb = 1; mrgl(l_)%pn = 0.0_R8; mrgl(l_)%pyn = 0.0_R8; mrgl(l_)%cr = 1; mrgl(l_)%cm = 1
        mrgl(l_)%cs = 1; mrgl(l_)%sp = 0; mrgl(l_)%pbr = 1; mrgl(l_)%pbs = 1
      end if
    end do
    ! mr-ghost-near: the 410 slots (Mod_BC_Q2D) of the level's blocks with a ghost the stage reads -- one of its layers, or
    ! at a face edge the ghost Ghost_Q2D_Edges writes beside the first layer (the transverse gather reads it)
    mrgl(l_)%nq = -1
    if ( near .and. q2d_n >= 0 ) then
      do pass = 1, 2
        nb_ = 0
        do q = 1, q2d_n
          idx = q2d_idx(q); b = q2d_blk(q)
          if ( mrb(b)%nbox(l_) <= 0 ) cycle
          n = domain%blk(b)%dim
          call field_block_index_params ( b, boff, pi, pj )
          Im = domain%bc_bi(idx); Jm = domain%bc_bj(idx); Km = domain%bc_bk(idx); Fm = domain%bc_bf(idx)
          keep = .false.
          do g = 1, gc
            if ( keep ) exit
            keep = nk( boff + (Im-g*guide(Fm,1)+gc-1) + (Jm-g*guide(Fm,2)+gc-1)*pi + (Km-g*guide(Fm,3)+gc-1)*pi*pj + 1 )
          end do
          if ( .not. keep ) then
            lo = (/ Im, Jm, Km /) - guide(Fm,1:3)   ! the first layer
            do dir = 1, 3
              if ( keep ) exit
              if ( guide(Fm,dir) /= 0 ) cycle
              do g = -1, 1, 2
                hi = lo; hi(dir) = lo(dir) + g
                if ( hi(dir) /= 0 .and. hi(dir) /= n(dir) + 1 ) cycle   ! an edge ghost only at the face's edge
                if ( nk( boff + (hi(1)+gc-1) + (hi(2)+gc-1)*pi + (hi(3)+gc-1)*pi*pj + 1 ) ) keep = .true.
              end do
            end do
          end if
          if ( .not. keep ) cycle
          nb_ = nb_ + 1
          if ( pass == 2 ) mrgl(l_)%qk(nb_) = q
        end do
        if ( pass == 1 ) then
          mrgl(l_)%nq = nb_
          allocate ( mrgl(l_)%qk(max(nb_,1)) ); mrgl(l_)%qk = 1
        end if
      end do
      !$acc enter data copyin(mrgl(l_)%qk)
    end if
    if ( allocated(nk) ) deallocate ( nk )
    !$acc enter data copyin(mrgl(l_)%pr, mrgl(l_)%pb, mrgl(l_)%pn, mrgl(l_)%pyn, mrgl(l_)%cr, mrgl(l_)%cm, mrgl(l_)%cs, &
    !$acc&                  mrgl(l_)%pbr, mrgl(l_)%pbs)
  end subroutine mr_build_ghost_list

  !> The four transverse neighbours of the donor cell of connection record c (Ghost_Connection_pairs_B), in the donor block.
  subroutine mr_conn_tnb ( domain, c, nq )
    use MOSE_Types_m
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer, intent(in)  :: c
    integer, intent(out) :: nq(3,4)
    integer :: Is, Js, Ks, d11, d12, d21, d22
    Is = domain%conn_is(c); Js = domain%conn_js(c); Ks = domain%conn_ks(c)
    d11 = domain%conn_d11(c); d12 = domain%conn_d12(c); d21 = domain%conn_d21(c); d22 = domain%conn_d22(c)
    select case ( domain%conn_fs(c) )
      case (1:2)
        nq(:,1) = (/ Is, Js-d11, Ks-d12 /); nq(:,2) = (/ Is, Js+d11, Ks+d12 /)
        nq(:,3) = (/ Is, Js-d21, Ks-d22 /); nq(:,4) = (/ Is, Js+d21, Ks+d22 /)
      case (3:4)
        nq(:,1) = (/ Is-d11, Js, Ks-d12 /); nq(:,2) = (/ Is+d11, Js, Ks+d12 /)
        nq(:,3) = (/ Is-d21, Js, Ks-d22 /); nq(:,4) = (/ Is+d21, Js, Ks+d22 /)
      case default
        nq(:,1) = (/ Is-d11, Js-d12, Ks /); nq(:,2) = (/ Is+d11, Js+d12, Ks /)
        nq(:,3) = (/ Is-d21, Js-d22, Ks /); nq(:,4) = (/ Is+d21, Js+d22, Ks /)
    end select
  end subroutine mr_conn_tnb

  !> The boundary-flux lists of level l_ (mrbl(:,l_)): the gather of Lib_Multirate_Base mr_gather_bc_list, shared with the
  !> lists of the explicit step (Mod_BC bc_list_setup).
  subroutine mr_build_bc_list ( domain, l_ )
    use MOSE_Types_m
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: l_
    logical :: ok
    call mr_gather_bc_list ( domain, l_, mrbl(:,l_), ok )
  end subroutine mr_build_bc_list

  !> The thermo-cache lists of level l_ (Lib_Multirate_Base mr_tl_type): per active block the union of its boxes of the
  !> level grown by gc and clipped to the block with its ghosts (ca), and the interior part of that union (cb) -- the cells
  !> Refresh_Thermo_Cache sweeps box by box -- each once, in block and i-fastest order.
  subroutine mr_build_thermo_list ( domain, l_ )
    use MOSE_Types_m
    use MOSE_Global_m, only: gc
    use MOSE_Mod_MPI,  only: is_local_block, mpi_is_root
    use MOSE_Field_m,  only: field_block_index_params
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: l_
    integer :: b, ib, i, j, k, n(3), lo(3), hi(3), boff, pi, pj, na, nb, pass
    logical, allocatable :: mk(:,:,:)
    do pass = 1, 2
      na = 0; nb = 0
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        if ( mrb(b)%nbox(l_) <= 0 ) cycle
        n = domain%blk(b)%dim
        call field_block_index_params ( b, boff, pi, pj )
        allocate ( mk(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc) ); mk = .false.
        do ib = 1, mrb(b)%nbox(l_)
          lo = (/ mrb(b)%box(1,ib,l_), mrb(b)%box(3,ib,l_), mrb(b)%box(5,ib,l_) /)
          hi = (/ mrb(b)%box(2,ib,l_), mrb(b)%box(4,ib,l_), mrb(b)%box(6,ib,l_) /)
          lo = max( lo - gc, 1 - gc ); hi = min( hi + gc, n + gc )
          mk(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) = .true.
        end do
        do k = 1-gc, n(3)+gc; do j = 1-gc, n(2)+gc; do i = 1-gc, n(1)+gc
          if ( .not. mk(i,j,k) ) cycle
          na = na + 1
          if ( pass == 2 ) mrtl(l_)%ca(na) = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
          if ( i < 1 .or. j < 1 .or. k < 1 .or. i > n(1) .or. j > n(2) .or. k > n(3) ) cycle
          nb = nb + 1
          if ( pass == 2 ) mrtl(l_)%cb(nb) = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
        end do; end do; end do
        deallocate ( mk )
      end do
      if ( pass == 1 ) then
        mrtl(l_)%na = na; mrtl(l_)%nb = nb
        allocate ( mrtl(l_)%ca(max(na,1)), mrtl(l_)%cb(max(nb,1)) )
        mrtl(l_)%ca = 1; mrtl(l_)%cb = 1
      end if
    end do
    !$acc enter data copyin(mrtl(l_)%ca, mrtl(l_)%cb)
    if ( mpi_is_root .and. obj_mr%nfull == 1 ) write(*,'(A,I0,A,I0,A,I0,A)') ' [MR]   thermo list level ', l_, &
      ': ', na, ' cells (', nb, ' interior), the boxes grown by gc'
  end subroutine mr_build_thermo_list

  subroutine mr_beta_one ( beta, n )
    implicit none
    integer,  intent(in)    :: n(3)
    real(R8), intent(inout) :: beta(n(1), n(2), n(3))
    integer :: i, j, k
    !$acc parallel loop collapse(3) present(beta)
    do k = 1, n(3); do j = 1, n(2); do i = 1, n(1)
      beta(i,j,k) = 1d0
    end do; end do; end do
    !$acc update host(beta)
  end subroutine mr_beta_one


  subroutine mr_build_iface ( domain )
    use MOSE_Types_m
    use MOSE_Global_m, only: nprim
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: b, d, i, j, k, n1, n2, n3, la, lb, nfd(3), cur(3), slot
    if ( allocated(gacc) ) then
      !$acc exit data delete(gacc, gsave)
      deallocate ( gacc, gsave )
    end if
    slot = 0
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      do d = 1, 3
        if ( allocated(mrb(b)%ifl(d)%f) ) then
          !$acc exit data delete(mrb(b)%ifl(d)%f, mrb(b)%ifl(d)%lc, mrb(b)%ifl(d)%s)
          deallocate ( mrb(b)%ifl(d)%f, mrb(b)%ifl(d)%lc, mrb(b)%ifl(d)%s )
        end if
        mrb(b)%ifl(d)%n = 0
      end do
      n1 = domain%blk(b)%dim(1); n2 = domain%blk(b)%dim(2); n3 = domain%blk(b)%dim(3)
      nfd = 0
      do k = 1, n3; do j = 1, n2; do i = 1, n1
        la = mrb(b)%lev(i,j,k)
        if ( i < n1 ) then; if ( abs(la - mrb(b)%lev(i+1,j,k)) == 1 ) nfd(1) = nfd(1) + 1; end if
        if ( j < n2 ) then; if ( abs(la - mrb(b)%lev(i,j+1,k)) == 1 ) nfd(2) = nfd(2) + 1; end if
        if ( k < n3 ) then; if ( abs(la - mrb(b)%lev(i,j,k+1)) == 1 ) nfd(3) = nfd(3) + 1; end if
      end do; end do; end do
      do d = 1, 3
        mrb(b)%ifl(d)%n = nfd(d)
        if ( nfd(d) > 0 ) allocate ( mrb(b)%ifl(d)%f(nfd(d)), mrb(b)%ifl(d)%lc(nfd(d)), mrb(b)%ifl(d)%s(nfd(d)) )
      end do
      cur = 0
      do k = 1, n3; do j = 1, n2; do i = 1, n1
        la = mrb(b)%lev(i,j,k)
        do d = 1, 3
          if ( d == 1 ) then; if ( i == n1 ) cycle; lb = mrb(b)%lev(i+1,j,k); end if
          if ( d == 2 ) then; if ( j == n2 ) cycle; lb = mrb(b)%lev(i,j+1,k); end if
          if ( d == 3 ) then; if ( k == n3 ) cycle; lb = mrb(b)%lev(i,j,k+1); end if
          if ( abs(la - lb) /= 1 ) cycle
          cur(d) = cur(d) + 1; slot = slot + 1
          ! the face index of the pass-1 buffers (Mod_Convective conv_cuf_i/j/k)
          if ( d == 1 ) mrb(b)%ifl(d)%f(cur(d)) = (i-1) + (j-1)*(n1-1) + (k-1)*(n1-1)*n2 + 1
          if ( d == 2 ) mrb(b)%ifl(d)%f(cur(d)) = (i-1) + (j-1)*n1 + (k-1)*n1*(n2-1) + 1
          if ( d == 3 ) mrb(b)%ifl(d)%f(cur(d)) = (i-1) + (j-1)*n1 + (k-1)*n1*n2 + 1
          mrb(b)%ifl(d)%lc(cur(d)) = min(la, lb)
          mrb(b)%ifl(d)%s(cur(d))  = slot
        end do
      end do; end do; end do
      do d = 1, 3
        if ( nfd(d) > 0 ) then
          !$acc enter data copyin(mrb(b)%ifl(d)%f, mrb(b)%ifl(d)%lc, mrb(b)%ifl(d)%s)
        end if
      end do
    end do
    obj_mr%nifc = slot
    allocate ( gacc(max(slot,1), nprim), gsave(max(slot,1), nprim) )
    gacc = 0.0_R8; gsave = 0.0_R8
    !$acc enter data copyin(gacc, gsave)
  end subroutine mr_build_iface


  !> The halo lists: for every level m < L, the level-m cells within two cells (any direction, the 5x5x5 neighbourhood) of
  !> a finer cell of the same block, as flat rows of P_flat, concatenated over the local blocks. Across a connection the
  !> two rows of each side have one level (R2), so a finer cell reads coarser cells only inside its own block.
  subroutine mr_build_halo ( domain, L )
    use MOSE_Types_m
    use MOSE_Global_m, only: gc, nprim
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Field_m, only: field_block_index_params
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: L
    integer :: m, b, i, j, k, nh, boff, pi, pj, n(3), e
    type :: ivec_type
      integer, allocatable :: v(:)
    end type ivec_type
    type(ivec_type), allocatable :: rows(:)
    integer, allocatable :: nrow(:)
    if ( allocated(mrh) ) then
      do m = lbound(mrh,1), ubound(mrh,1)
        if ( allocated(mrh(m)%c) ) then
          !$acc exit data delete(mrh(m)%c, mrh(m)%ws)
          deallocate ( mrh(m)%c, mrh(m)%ws )
        end if
      end do
      deallocate ( mrh )
    end if
    allocate ( mrh(0:L) )
    allocate ( rows(domain%nb), nrow(domain%nb) )
    do m = 0, L - 1
      nh = 0; nrow = 0
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        call mr_halo_mark_blk ( mrb(b)%lev, mrb(b)%lev2, domain%blk(b)%dim, m )
        !$acc update host(mrb(b)%lev2)
        nrow(b) = count( mrb(b)%lev2 == 1 )
        if ( nrow(b) == 0 ) cycle
        allocate ( rows(b)%v(nrow(b)) )
        call field_block_index_params ( b, boff, pi, pj )
        n = domain%blk(b)%dim; e = 0
        do k = 1, n(3); do j = 1, n(2); do i = 1, n(1)
          if ( mrb(b)%lev2(i,j,k) /= 1 ) cycle
          e = e + 1
          rows(b)%v(e) = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
        end do; end do; end do
        nh = nh + nrow(b)
      end do
      mrh(m)%n = nh
      allocate ( mrh(m)%c(max(nh,1)), mrh(m)%ws(max(nh,1), nprim) )
      mrh(m)%c = 1; mrh(m)%ws = 0.0_R8
      nh = 0
      do b = 1, domain % nb
        if ( nrow(b) == 0 ) cycle
        mrh(m)%c(nh+1:nh+nrow(b)) = rows(b)%v(1:nrow(b))
        nh = nh + nrow(b)
        deallocate ( rows(b)%v )
      end do
      !$acc enter data copyin(mrh(m)%c, mrh(m)%ws)
    end do
    deallocate ( rows, nrow )
    mrh(L)%n = 0
  end subroutine mr_build_halo

  !> The halo cells each list level reads (Lib_Multirate_Base mr_hn_type), with mr-ghost-near: the cells of the halo of
  !> every coarser level m in the rows the stage of level l reads -- the boxes of l grown by gc (thermo cache, face stencils,
  !> boundary records, SST source) and the interior cells the ghost fill of its records reads (Ghost_Extrapolate: the
  !> record's cell and the next two along the normal; Ghost_Connection_list: the donor's gc cells along its normal and its
  !> four transverse neighbours, which Ghost_Connection_pairs_B reads too). The other halo cells are interpolated by the
  !> stages that read them; the saves and the restore of the predictor (MR_Halo_Save, MR_Halo_Interp at 1) keep the whole
  !> halo. Self pairs (their own launch over the whole pair): every halo cell.
  subroutine mr_build_halo_near ( domain, L )
    use MOSE_Types_m
    use MOSE_Global_m,     only: gc
    use MOSE_Parameters_m, only: guide
    use MOSE_Mod_MPI,      only: is_local_block, mpi_is_root
    use MOSE_Field_m,      only: ncell_total, field_block_index_params
    use MOSE_Mod_BC_Q2D,   only: q2d_n, q2d_idx, q2d_blk
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: L
    integer :: l_, m, b, ib, i, j, k, e, d, q, c, idx, nh, boff, pi, pj, Fm
    integer :: n(3), lo(3), hi(3), x(3), nq(3,4)
    logical, allocatable :: rd(:)
    if ( allocated(mrhn) ) then
      do l_ = lbound(mrhn,2), ubound(mrhn,2)
        do m = lbound(mrhn,1), ubound(mrhn,1)
          if ( .not. allocated(mrhn(m,l_)%c) ) cycle
          !$acc exit data delete(mrhn(m,l_)%c, mrhn(m,l_)%pos)
          deallocate ( mrhn(m,l_)%c, mrhn(m,l_)%pos )
        end do
      end do
      deallocate ( mrhn )
    end if
    if ( .not. obj_mr%ghost_near ) return
    allocate ( mrhn(0:L, 0:L) )
    allocate ( rd(ncell_total) )
    do l_ = max(obj_mr%listlev, 1), L
      rd = ( mrgl(l_)%nself > 0 )
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        if ( mrb(b)%nbox(l_) <= 0 ) cycle
        n = domain%blk(b)%dim
        call field_block_index_params ( b, boff, pi, pj )
        do ib = 1, mrb(b)%nbox(l_)
          lo = (/ mrb(b)%box(1,ib,l_), mrb(b)%box(3,ib,l_), mrb(b)%box(5,ib,l_) /)
          hi = (/ mrb(b)%box(2,ib,l_), mrb(b)%box(4,ib,l_), mrb(b)%box(6,ib,l_) /)
          lo = max( lo - gc, 1 - gc ); hi = min( hi + gc, n + gc )
          do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
            rd( boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1 ) = .true.
          end do; end do; end do
        end do
      end do
      do e = 1, mrgl(l_)%np                      ! physical records: the cell and two more along the normal
        idx = mrgl(l_)%pr(e); b = mrgl(l_)%pb(e); Fm = domain%bc_bf(idx)
        n = domain%blk(b)%dim
        call field_block_index_params ( b, boff, pi, pj )
        do d = 0, 2
          x = (/ domain%bc_bi(idx), domain%bc_bj(idx), domain%bc_bk(idx) /) + d*guide(Fm,1:3)
          if ( any(x < 1) .or. any(x > n) ) exit
          rd( boff + (x(1)+gc-1) + (x(2)+gc-1)*pi + (x(3)+gc-1)*pi*pj + 1 ) = .true.
        end do
      end do
      do e = 1, mrgl(l_)%nc + mrgl(l_)%npb       ! connection records: the donor cells, then the transverse neighbours
        if ( e <= mrgl(l_)%nc ) then
          c = mrgl(l_)%cr(e); b = mrgl(l_)%cs(e)
        else
          c = mrgl(l_)%pbr(e-mrgl(l_)%nc); b = mrgl(l_)%pbs(e-mrgl(l_)%nc)
        end if
        if (.not. is_local_block(b)) cycle
        n = domain%blk(b)%dim
        call field_block_index_params ( b, boff, pi, pj )
        if ( e <= mrgl(l_)%nc ) then
          do d = 0, gc - 1
            x = (/ domain%conn_is(c), domain%conn_js(c), domain%conn_ks(c) /) + d*guide(domain%conn_fs(c),1:3)
            if ( any(x < 1) .or. any(x > n) ) exit
            rd( boff + (x(1)+gc-1) + (x(2)+gc-1)*pi + (x(3)+gc-1)*pi*pj + 1 ) = .true.
          end do
        end if
        call mr_conn_tnb ( domain, c, nq )
        do q = 1, 4
          if ( any(nq(:,q) < 1-gc) .or. any(nq(:,q) > n+gc) ) cycle
          rd( boff + (nq(1,q)+gc-1) + (nq(2,q)+gc-1)*pi + (nq(3,q)+gc-1)*pi*pj + 1 ) = .true.
        end do
      end do
      do e = 1, max(mrgl(l_)%nq, 0)              ! 410 slots: the boundary cell (Ghost_Q2D_Layers)
        idx = q2d_idx(mrgl(l_)%qk(e)); b = q2d_blk(mrgl(l_)%qk(e))
        call field_block_index_params ( b, boff, pi, pj )
        rd( boff + (domain%bc_bi(idx)+gc-1) + (domain%bc_bj(idx)+gc-1)*pi + (domain%bc_bk(idx)+gc-1)*pi*pj + 1 ) = .true.
      end do
      if ( mrgl(l_)%nq < 0 ) rd = .true.          ! every 410 slot filled: no selection (their cells are not listed here)
      do m = 0, l_ - 1
        nh = 0
        do e = 1, mrh(m)%n
          if ( rd(mrh(m)%c(e)) ) nh = nh + 1
        end do
        mrhn(m,l_)%n = nh
        allocate ( mrhn(m,l_)%c(max(nh,1)), mrhn(m,l_)%pos(max(nh,1)) )
        mrhn(m,l_)%c = 1; mrhn(m,l_)%pos = 1
        nh = 0
        do e = 1, mrh(m)%n
          if ( .not. rd(mrh(m)%c(e)) ) cycle
          nh = nh + 1
          mrhn(m,l_)%c(nh) = mrh(m)%c(e); mrhn(m,l_)%pos(nh) = e
        end do
        !$acc enter data copyin(mrhn(m,l_)%c, mrhn(m,l_)%pos)
      end do
      if ( mpi_is_root .and. obj_mr%nfull == 1 ) write(*,'(A,I0,A,9(1X,I0))') ' [MR]   halo cells read by level ', l_, &
        ' (of the halos of the coarser levels):', ( mrhn(m,l_)%n, m = 0, l_ - 1 )
    end do
    deallocate ( rd )
  end subroutine mr_build_halo_near

  !> lev2 = 1 on the level-m cells that have a finer cell in their 5x5x5 neighbourhood inside the block, else 0.
  subroutine mr_halo_mark_blk ( lev, lev2, n, m )
    implicit none
    integer, intent(in)    :: n(3), m, lev(n(1),n(2),n(3))
    integer, intent(inout) :: lev2(n(1),n(2),n(3))
    integer :: i, j, k, n1, n2, n3, ii, jj, kk, f
    n1 = n(1); n2 = n(2); n3 = n(3)
    !$acc parallel loop collapse(3) gang vector present(lev, lev2) private(f, ii, jj, kk)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      f = 0
      if ( lev(i,j,k) == m ) then
        do kk = max(k-2,1), min(k+2,n3)
        do jj = max(j-2,1), min(j+2,n2)
        do ii = max(i-2,1), min(i+2,n1)
          if ( lev(ii,jj,kk) > m ) f = 1
        end do
        end do
        end do
      end if
      lev2(i,j,k) = f
    end do
    end do
    end do
  end subroutine mr_halo_mark_blk

end module MOSE_Mod_Multirate_Partition
