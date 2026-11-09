!> Multirate: the driver of the levels, in place of the RK loop of Explicit_Step when multirate = t (doc/archive/
!> MULTIRATE_PLAN.md of branch multirate, sec. 1.3 and 1.4). Over the macro step DT a cell of level l takes 2**l Heun steps
!> of h_l = DT / 2**l:
!>   MR_Advance(l, k0):  wo := w on the level-l cells
!>                       Heun stage 1 of level l at k0
!>                       if l < L: stash the level-l halo, MR_Advance(l+1, k0), MR_Advance(l+1, k0 + h_l/2), put it back
!>                       Heun stage 2 of level l at k0 + h_l
!> At stage 1 of level l every finer level is at the start of the step, at stage 2 at its end: only the coarser levels need a
!> state at the stage time, w = wo + theta (w1 - wo) on their halo (wo the start, w1 the stage-1 predictor; theta = 0 and 1
!> give wo and w1 exactly). Times are counted in ticks of the finest sub-step DT / 2**L, so theta and the stage times of the
!> boundary data (667 = MOSE_open 410) are exact; with one level the stages are those of the explicit RK2 step, bit for bit.
!> The stage kernels work on the boxes of the stage's level (Lib_Multirate_Base) and update only its cells.
module MOSE_Mod_Multirate
  use iso_fortran_env, only: R8 => real64, I8 => int64
  use MOSE_Lib_Multirate_Base
  implicit none
  private
  public :: MR_Macro_Step, MR_Macro_Step_Redo
  real(R8), allocatable, save :: P0_flat(:,:)   !< mr-nonlocal = redo: the state at the start of the macro step

  integer(I8), save :: tau_tick(0:MR_LMAX_HARD) = 0   !< the tick each level's current step started at
  !> per level: the host time to issue a stage (start of MR_Stage to its final wait) and the wait for the device after it
  real(R8), save :: t_issue(0:MR_LMAX_HARD) = 0.0_R8, t_wait(0:MR_LMAX_HARD) = 0.0_R8

contains


  !> One macro step of the ladder over [time_old, time_old + DT], DT = domain%dtglobal (Compute_dt + Set_Global_dt set the
  !> levels and each cell's sub-step). Called outside any OpenMP parallel region.
  subroutine MR_Macro_Step ( domain, config, External_Function, time_old )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: mpi_is_root
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    type(config_type),      intent(inout) :: config
    external :: External_Function
    real(R8), intent(in) :: time_old
    integer    :: l, c
    integer(8) :: c0_, c1_, cr_
    character(len=6), parameter :: cname(MR_NCALL) = (/ '  halo', ' ghost', 'thermo', '  flux', '    bc', '   ext', &
                                                        '  rans', '   new' /)
    call system_clock ( c0_, cr_ )
    do l = 0, obj_mr%L
      obj_mr%h_lev(l) = domain % dtglobal / ( 2.0_R8 ** l )   ! = the dtlocal of the level's cells (MR_Set_Local_dt)
    end do
    if ( obj_mr%nonlocal > 0 ) call mr_cert_open ( domain )   ! mr-nonlocal: dmin := dtcap (the steps at t_n), count 0
    call MR_Advance ( domain, config, External_Function, 0, 0_I8, time_old )
    if ( obj_mr%nonlocal > 0 ) call mr_cert_close ()          ! the violations of the macro step, over the ranks
    call system_clock ( c1_ )
    obj_mr%tmacro = obj_mr%tmacro + real(c1_ - c0_, R8) / real(cr_, R8)
    if ( mpi_is_root .and. obj_mr%nprint > 0 ) then
      if ( mod(domain%iter, obj_mr%nprint) == 0 ) then
        write(*,'(A)',advance='no') ' [MR]   stage wall s / stages per level (rank 0):'
        do l = 0, obj_mr%L
          write(*,'(A,I0,A,F9.2,A,I0)',advance='no') '  L', l, ' ', obj_mr%tlev(l), ' / ', obj_mr%nlev(l)
        end do
        write(*,'(A)') ''
        write(*,'(A)',advance='no') ' [MR]   host issue s / device wait s per level (rank 0):'
        do l = 0, obj_mr%L
          write(*,'(A,I0,A,F8.2,A,F8.2)',advance='no') '  L', l, ' ', t_issue(l), ' / ', t_wait(l)
        end do
        write(*,'(A)') ''
        write(*,'(A,F10.2,A,F10.2,A,F9.2,A,F9.2)') ' [MR]   wall s (rank 0): macro steps ', obj_mr%tmacro, '  stages ', &
          sum(obj_mr%tlev(0:obj_mr%L)), '  partition ', obj_mr%tpart, '  of which full ', obj_mr%tfull
        if ( obj_mr%nonlocal > 0 ) write(*,'(A,I0,A,I0)') ' [MR]   certificate (mr-nonlocal): violations ', &
          obj_mr%ncert_viol, '  in macro steps ', obj_mr%ncert_steps
        if ( obj_mr%profile ) then
          write(*,'(A,8(A,A))') ' [MR]   profile wall ms per stage (rank 0):', ('  ', cname(c), c = 1, MR_NCALL)
          do l = 0, obj_mr%L
            write(*,'(A,I0,A,8(F8.2))') ' [MR]     L', l, '                               ', &
              ( 1.0e3_R8 * obj_mr%tcall(c,l) / real(max(obj_mr%nlev(l),1), R8), c = 1, MR_NCALL )
          end do
          write(*,'(A)') ' [MR]   profile ghost ms per stage (rank 0):     phys  q2dval  q2dlay    conn  q2dedg   passB   chim' &
                         // '  410 stages'
          do l = 0, obj_mr%L
            write(*,'(A,I0,A,7(F8.2),I8)') ' [MR]     L', l, '                               ', &
              ( 1.0e3_R8 * obj_mr%tgh(c,l) / real(max(obj_mr%nlev(l),1), R8), c = 1, 7 ), nint(obj_mr%tgh(8,l))
          end do
          write(*,'(A,F9.2,A,F9.2,A,F9.2,A,F9.2)') ' [MR]   profile wall s (rank 0): Compute_dt ', obj_mr%tit(1), &
            ' (partition ', obj_mr%tpart, ')  rest of the step ', obj_mr%tit(2), '  between steps ', obj_mr%tit(3)
          write(*,'(A,I0,A,9(F8.2))') ' [MR]   profile full partitions (', obj_mr%nfull, ') wall s: levels boxes ' // &
            'iface+bc faces+cells ghost bclist thermo halo halo_near', obj_mr%tpb(1:9)
        end if
      end if
    end if
  end subroutine MR_Macro_Step


  recursive subroutine MR_Advance ( domain, config, External_Function, l, k0, time_old )
    use MOSE_Types_m
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    type(config_type),      intent(inout) :: config
    external :: External_Function
    integer,     intent(in) :: l
    integer(I8), intent(in) :: k0
    real(R8),    intent(in) :: time_old
    integer(I8) :: hk
    hk = 2_I8 ** ( obj_mr%L - l )          ! the level's step in ticks
    tau_tick(l) = k0
    call MR_Copy_State ( domain, l )
    call MR_Stage ( domain, config, External_Function, l, 1, k0, time_old )
    if ( l < obj_mr%L ) then
      call MR_Halo_Save ( l )
      call MR_Advance ( domain, config, External_Function, l + 1, k0, time_old )
      call MR_Advance ( domain, config, External_Function, l + 1, k0 + hk/2, time_old )
      call MR_Halo_Interp ( l, 1.0_R8 )     ! the stage-1 predictor back on the halo for stage 2
    end if
    call MR_Stage ( domain, config, External_Function, l, 2, k0 + hk, time_old )
  end subroutine MR_Advance


  !> One Heun stage of level l at tick k: the stage sequence of Explicit_Step with the level's boxes armed.
  subroutine MR_Stage ( domain, config, External_Function, l, irk, k, time_old )
    use MOSE_Types_m
    use MOSE_Global_m,          only: RANS, VISC
    use MOSE_Lib_Extrapolation, only: Fill_Ghost_Cell, Extrapolate_Wall
    use MOSE_Mod_Thermo_Cache,  only: Refresh_Thermo_Cache, thermo_guard_report
    use MOSE_Mod_Fluxes,        only: Fluxes
    use MOSE_Mod_BC,            only: BC_Fluxes
    use MOSE_Lib_RANS,          only: RANS_Source_Terms
    use MOSE_Lib_Newstate,      only: RK_Newstate
    use MOSE_Mod_BC_Q2D,        only: q2d_time, q2d_time_set
    use MOSE_Mod_GhostExchange, only: xch_set_level
    use MOSE_Mod_MPI,           only: halo_level
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    type(config_type),      intent(inout) :: config
    external :: External_Function
    integer,     intent(in) :: l, irk
    integer(I8), intent(in) :: k
    real(R8),    intent(in) :: time_old
    integer    :: m
    integer(8) :: c0_, c1_, cr_, ct_, cw_, cx_
    real(R8)   :: theta
    call system_clock ( c0_, cr_ )
    ct_ = c0_
    ! mr-async: the queue of this stage's kernels (the async queue on a list level); a stage off the list path (the block
    ! kernels, synchronous queue) first joins the chain, so nothing a list kernel still writes is read by a block kernel
    call mr_async_setup ( domain % nb )
    call mr_set_queue ( l )
    if ( mr_async .and. mr_q /= mr_qid ) then
      !$acc wait(mr_qid)
    end if
    ! the coarser levels' halo at this stage's time
    do m = 0, l - 1
      theta = real( k - tau_tick(m), R8 ) / real( 2_I8 ** (obj_mr%L - m), R8 )
      if ( mr_halo_near(l) ) then
        call MR_Halo_Interp_Near ( m, l, theta )   ! mr-ghost-near: only the halo cells the stage reads
      else
        call MR_Halo_Interp ( m, theta )
      end if
    end do
    call mr_toc ( 1, l, ct_ )
    obj_mr%lcur = l; obj_mr%irkcur = irk; obj_mr%in_stage = .true.
    ! boundary data of type 667 (MOSE_open 410): the state of the file at the stage time
    q2d_time = time_old + real(k, R8) * domain % dtglobal / ( 2.0_R8 ** obj_mr%L )
    q2d_time_set = .true.
    call xch_set_level ( l )              ! mr-xch: the inter-rank exchange of this stage ships the level's sub-list
    halo_level = l                        ! the halo phase clocks of this stage go to the level's row (Mod_MPI halo_tl)
    call Fill_Ghost_Cell ( domain );      call mr_toc ( 2, l, ct_ )
    halo_level = -1
    call xch_set_level ( -1 )
    call Refresh_Thermo_Cache ( domain ); call mr_toc ( 3, l, ct_ )
    if ( mr_q_async() .and. .not. mr_list_on() ) then
      ! mr-async-blk: the block-path fluxes run on the synchronous stream (CUDA Fortran kernels, the interface correction,
      ! the pass 2): join the queue, run them synchronously, back to the queue for the boundary fluxes and the rest
      !$acc wait(mr_q)
      call mr_set_queue_sync ()
      call Fluxes ( domain )
      call mr_set_queue ( l )
    else
      call Fluxes ( domain )
    end if
    call mr_toc ( 4, l, ct_ )
    call BC_Fluxes ( domain );            call mr_toc ( 5, l, ct_ )
    call External_Function ( domain );    call mr_toc ( 6, l, ct_ )
    if ( RANS ) then
      call Extrapolate_Wall ( domain )
      call RANS_Source_Terms ( domain )
    end if
    call mr_toc ( 7, l, ct_ )
    ! mr-nonlocal: the certificate of condition (4) at the start of the sub-step, on the start state (stage 1)
    if ( obj_mr%nonlocal > 0 .and. irk == 1 .and. obj_mr%L > 0 ) call mr_certificate ( domain, config, l, k )
    call RK_Newstate ( domain, irk )
    obj_mr%in_stage = .false.
    call system_clock ( cw_ )                 ! the host has issued the whole stage: the rest is the device catching up
    !$acc wait
    call system_clock ( cx_ )
    t_issue(l) = t_issue(l) + real(cw_ - c0_, R8) / real(cr_, R8)
    t_wait(l)  = t_wait(l)  + real(cx_ - cw_, R8) / real(cr_, R8)
    call mr_toc ( 8, l, ct_ )
    ! mr-async: the thermodynamic guard (a device counter, read on the host) after the stage's wait instead of mid-stage
    if ( VISC .and. mr_q_async() ) call thermo_guard_report ( -1, domain%blk(1)%dim, domain%iter )
    call system_clock ( c1_ )
    obj_mr%tlev(l) = obj_mr%tlev(l) + real(c1_ - c0_, R8) / real(cr_, R8)
    obj_mr%nlev(l) = obj_mr%nlev(l) + 1
  end subroutine MR_Stage

  !> mr-profile: the wall since ct_ (after a device wait) goes to call c of level l; ct_ moves to now.
  subroutine mr_toc ( c, l, ct_ )
    implicit none
    integer,    intent(in)    :: c, l
    integer(8), intent(inout) :: ct_
    integer(8) :: c1_, cr_
    if ( .not. obj_mr%profile ) return
    !$acc wait
    call system_clock ( c1_, cr_ )
    obj_mr%tcall(c,l) = obj_mr%tcall(c,l) + real(c1_ - ct_, R8) / real(cr_, R8)
    ct_ = c1_
  end subroutine mr_toc


  !> wo := w on the cells of level l (the start of their step). One level: the whole-field copy of the explicit step.
  subroutine MR_Copy_State ( domain, l )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Field_m, only: field_copy_P_flat_to_PO_flat
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: l
    integer :: b, ib
    if ( obj_mr%L == 0 ) then
      call field_copy_P_flat_to_PO_flat ()
      return
    end if
    ! a list level (mr-list-level): the level's cells of every block in one launch, the cells the box copy takes
    if ( obj_mr%listlev >= 0 .and. l >= obj_mr%listlev ) then
      if ( mrcl(l)%n > 0 ) call mr_copy_level_list ( mrcl(l)%c, mrcl(l)%n )
      return
    end if
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      do ib = 1, mrb(b)%nbox(l)
        call mr_copy_level_blk ( b, domain%blk(b)%dtlocal, domain%blk(b)%dim, mrb(b)%box(:,ib,l), obj_mr%h_lev(l) )
      end do
    end do
  end subroutine MR_Copy_State

  subroutine mr_copy_level_blk ( b, dtl, n, bx, hsel )
    use MOSE_Global_m, only: gc, nprim
    use MOSE_Field_m,  only: P_flat, PO_flat, field_block_index_params
    implicit none
    integer,  intent(in) :: b, n(3), bx(6)
    real(R8), intent(in) :: dtl(n(1),n(2),n(3)), hsel
    integer  :: i, j, k, v, c0, boff, pi, pj, i0, i1, j0, j1, k0, k1
    real(R8) :: h_
    call field_block_index_params ( b, boff, pi, pj )
    i0 = bx(1); i1 = bx(2); j0 = bx(3); j1 = bx(4); k0 = bx(5); k1 = bx(6); h_ = hsel
    !$acc parallel loop collapse(3) gang vector present(P_flat, PO_flat, dtl) private(c0, v)
    do k = k0, k1
    do j = j0, j1
    do i = i0, i1
      if ( dtl(i,j,k) /= h_ ) cycle
      c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      do v = 1, nprim
        PO_flat(c0, v) = P_flat(c0, v)
      end do
    end do
    end do
    end do
  end subroutine mr_copy_level_blk


  subroutine mr_copy_level_list ( cl, n )
    use MOSE_Global_m, only: nprim
    use MOSE_Field_m,  only: P_flat, PO_flat
    implicit none
    integer, intent(in) :: n, cl(n)
    integer :: e, v
    !$acc parallel loop collapse(2) gang vector present(P_flat, PO_flat, cl) async(mr_q)
    do v = 1, nprim
      do e = 1, n
        PO_flat(cl(e), v) = P_flat(cl(e), v)
      end do
    end do
  end subroutine mr_copy_level_list

  !> The stage-1 predictor of the level-m halo cells, stashed before the finer levels run.
  subroutine MR_Halo_Save ( m )
    implicit none
    integer, intent(in) :: m
    if ( .not. allocated(mrh) ) return
    if ( mrh(m)%n == 0 ) return
    call mr_halo_save_list ( mrh(m)%c, mrh(m)%ws, mrh(m)%n )
  end subroutine MR_Halo_Save

  subroutine mr_halo_save_list ( c, ws, n )
    use MOSE_Global_m, only: nprim
    use MOSE_Field_m,  only: P_flat
    implicit none
    integer,  intent(in)    :: n, c(n)
    real(R8), intent(inout) :: ws(n, nprim)
    integer :: e, v
    !$acc parallel loop collapse(2) gang vector present(P_flat, c, ws) async(mr_q)
    do v = 1, nprim
      do e = 1, n
        ws(e, v) = P_flat(c(e), v)
      end do
    end do
  end subroutine mr_halo_save_list

  !> w = wo + theta (w1 - wo) on the level-m halo cells: wo exactly at theta = 0, the stashed predictor w1 at theta = 1.
  subroutine MR_Halo_Interp ( m, theta )
    implicit none
    integer,  intent(in) :: m
    real(R8), intent(in) :: theta
    if ( .not. allocated(mrh) ) return
    if ( mrh(m)%n == 0 ) return
    call mr_halo_interp_list ( mrh(m)%c, mrh(m)%ws, mrh(m)%n, theta )
  end subroutine MR_Halo_Interp

  !> A list level with mr-ghost-near reads its own subset of each coarser halo (Mod_Multirate_Partition mr_build_halo_near).
  logical function mr_halo_near ( l )
    integer, intent(in) :: l
    mr_halo_near = .false.
    if ( .not. allocated(mrhn) ) return
    if ( obj_mr%listlev < 0 .or. l < max(obj_mr%listlev, 1) .or. l > ubound(mrhn,2) ) return
    mr_halo_near = .true.
  end function mr_halo_near

  !> MR_Halo_Interp on the cells of the halo of m that the stage of level l reads: the same expression, the predictor of
  !> each cell from its row of the whole halo's ws.
  subroutine MR_Halo_Interp_Near ( m, l, theta )
    implicit none
    integer,  intent(in) :: m, l
    real(R8), intent(in) :: theta
    if ( mrhn(m,l)%n == 0 ) return
    call mr_halo_interp_near_list ( mrhn(m,l)%c, mrhn(m,l)%pos, mrh(m)%ws, mrhn(m,l)%n, mrh(m)%n, theta )
  end subroutine MR_Halo_Interp_Near

  subroutine mr_halo_interp_near_list ( c, pos, ws, n, nw, theta )
    use MOSE_Global_m, only: nprim, NPRIM_DEV
    use FLINT_Lib_Thermodynamic, only: Tmin, Tmax
    use MOSE_Field_m,  only: P_flat, PO_flat
    implicit none
    integer,  intent(in) :: n, nw, c(n), pos(n)
    real(R8), intent(in) :: ws(nw, nprim), theta
    integer  :: e, v, tmin_, tmax_
    real(R8) :: th, po_(NPRIM_DEV), pw_(NPRIM_DEV), pout_(NPRIM_DEV)
    th = theta
    if ( th <= 0.0_R8 ) then
      !$acc parallel loop collapse(2) gang vector present(P_flat, PO_flat, c) async(mr_q)
      do v = 1, nprim
        do e = 1, n
          P_flat(c(e), v) = PO_flat(c(e), v)
        end do
      end do
    else if ( th >= 1.0_R8 ) then
      !$acc parallel loop collapse(2) gang vector present(P_flat, c, pos, ws) async(mr_q)
      do v = 1, nprim
        do e = 1, n
          P_flat(c(e), v) = ws(pos(e), v)
        end do
      end do
    else if ( obj_mr%cons_halo ) then   ! mr-cons-halo: in the conservative variables (mr_cons_interp)
      tmin_ = Tmin; tmax_ = Tmax
      !$acc parallel loop gang vector present(P_flat, PO_flat, c, pos, ws) async(mr_q) private(po_, pw_, pout_, v)
      do e = 1, n
        do v = 1, nprim
          po_(v) = PO_flat(c(e), v); pw_(v) = ws(pos(e), v)
        end do
        call mr_cons_interp ( po_, pw_, th, pout_, tmin_, tmax_ )
        do v = 1, nprim
          P_flat(c(e), v) = pout_(v)
        end do
      end do
    else
      !$acc parallel loop collapse(2) gang vector present(P_flat, PO_flat, c, pos, ws) async(mr_q)
      do v = 1, nprim
        do e = 1, n
          P_flat(c(e), v) = PO_flat(c(e), v) + th * ( ws(pos(e), v) - PO_flat(c(e), v) )
        end do
      end do
    end if
  end subroutine mr_halo_interp_near_list

  subroutine mr_halo_interp_list ( c, ws, n, theta )
    use MOSE_Global_m, only: nprim, NPRIM_DEV
    use FLINT_Lib_Thermodynamic, only: Tmin, Tmax
    use MOSE_Field_m,  only: P_flat, PO_flat
    implicit none
    integer,  intent(in) :: n, c(n)
    real(R8), intent(in) :: ws(n, nprim), theta
    integer  :: e, v, tmin_, tmax_
    real(R8) :: th, po_(NPRIM_DEV), pw_(NPRIM_DEV), pout_(NPRIM_DEV)
    th = theta
    if ( th <= 0.0_R8 ) then
      !$acc parallel loop collapse(2) gang vector present(P_flat, PO_flat, c) async(mr_q)
      do v = 1, nprim
        do e = 1, n
          P_flat(c(e), v) = PO_flat(c(e), v)
        end do
      end do
    else if ( th >= 1.0_R8 ) then
      !$acc parallel loop collapse(2) gang vector present(P_flat, c, ws) async(mr_q)
      do v = 1, nprim
        do e = 1, n
          P_flat(c(e), v) = ws(e, v)
        end do
      end do
    else if ( obj_mr%cons_halo ) then   ! mr-cons-halo: in the conservative variables (mr_cons_interp)
      tmin_ = Tmin; tmax_ = Tmax
      !$acc parallel loop gang vector present(P_flat, PO_flat, c, ws) async(mr_q) private(po_, pw_, pout_, v)
      do e = 1, n
        do v = 1, nprim
          po_(v) = PO_flat(c(e), v); pw_(v) = ws(e, v)
        end do
        call mr_cons_interp ( po_, pw_, th, pout_, tmin_, tmax_ )
        do v = 1, nprim
          P_flat(c(e), v) = pout_(v)
        end do
      end do
    else
      !$acc parallel loop collapse(2) gang vector present(P_flat, PO_flat, c, ws) async(mr_q)
      do v = 1, nprim
        do e = 1, n
          P_flat(c(e), v) = PO_flat(c(e), v) + th * ( ws(e, v) - PO_flat(c(e), v) )
        end do
      end do
    end if
  end subroutine mr_halo_interp_list


  !> mr-cons-halo (open-hydra MOSE feature/multirate-integrator, Mod_Multirate.f90 1929-1958 there; Hoang et al. 2020, eq.
  !> 3.2): the coarse halo state at 0 < theta < 1 interpolated in the conservative variables, U = Ub + theta (Us - Ub), with
  !> Ub and Us the conservative states of the start po and of the stashed predictor pw (prim2cons of newstate_cons_body.inc),
  !> then cons2prim with the Newton iteration of newstate_cons_body.inc seeded with T of the predictor (as the branch). The
  !> species densities and the turbulence variables are the same quantities in both sets: their value is the expression of
  !> the primitive interpolation, to the bit; u, v, w and p change. No clamp: U is a convex combination of two states.
  subroutine mr_cons_interp ( po, pw, th, pout, tmin_, tmax_ )
    !$acc routine seq
    use FLINT_Lib_Thermodynamic, only: h_tab, cp_tab, f_Rtot
    use MOSE_Global_m, only: nprim, NPRIM_DEV, nsc, nu, nw, np
    implicit none
    real(R8), intent(in)  :: po(NPRIM_DEV), pw(NPRIM_DEV), th
    real(R8), intent(out) :: pout(NPRIM_DEV)
    integer,  intent(in)  :: tmin_, tmax_
    real(R8) :: ub(NPRIM_DEV), us(NPRIM_DEV), u(NPRIM_DEV)
    real(R8) :: rho, Rg, vmag, Tloc, Td, hmix, h_s, cp_s, E0, energy, diff_en, T_it, p_it, cv, enth
    integer  :: s, T_i, it

    ! prim2cons(po) -> ub
    rho  = sum( po(1:nsc) )
    Rg   = f_Rtot( po(1:nsc) )
    vmag = 0d0; do s = nu, nw; vmag = vmag + po(s)*po(s); end do
    Tloc = po(np) / (Rg*rho)
    T_i = int(Tloc); if (T_i < tmin_) T_i = tmin_; if (T_i > tmax_-1) T_i = tmax_-1; Td = Tloc - T_i
    hmix = 0d0
    do s = 1, nsc
      h_s  = h_tab(T_i,s) + Td*( h_tab(T_i+1,s) - h_tab(T_i,s) )
      hmix = hmix + h_s * po(s)/rho
    end do
    E0 = hmix - po(np)/rho + 0.5d0*vmag
    ub(1:nsc) = po(1:nsc)
    ub(nu:nw) = rho * po(nu:nw)
    ub(np)    = rho * E0
    if (nprim>np) ub(np+1:nprim) = po(np+1:nprim)

    ! prim2cons(pw) -> us, and the Newton seed T of the predictor
    rho  = sum( pw(1:nsc) )
    Rg   = f_Rtot( pw(1:nsc) )
    vmag = 0d0; do s = nu, nw; vmag = vmag + pw(s)*pw(s); end do
    Tloc = pw(np) / (Rg*rho)
    T_it = Tloc
    T_i = int(Tloc); if (T_i < tmin_) T_i = tmin_; if (T_i > tmax_-1) T_i = tmax_-1; Td = Tloc - T_i
    hmix = 0d0
    do s = 1, nsc
      h_s  = h_tab(T_i,s) + Td*( h_tab(T_i+1,s) - h_tab(T_i,s) )
      hmix = hmix + h_s * pw(s)/rho
    end do
    E0 = hmix - pw(np)/rho + 0.5d0*vmag
    us(1:nsc) = pw(1:nsc)
    us(nu:nw) = rho * pw(nu:nw)
    us(np)    = rho * E0
    if (nprim>np) us(np+1:nprim) = pw(np+1:nprim)

    do s = 1, nprim
      u(s) = ub(s) + th * ( us(s) - ub(s) )
    end do

    ! cons2prim (the Newton-Raphson of newstate_cons_body.inc)
    rho  = sum( u(1:nsc) )
    Rg   = f_Rtot( u(1:nsc) )
    pout(1:nsc) = u(1:nsc)
    pout(nu:nw) = u(nu:nw)/rho
    vmag   = 0d0; do s = nu, nw; vmag = vmag + (u(s)/rho)*(u(s)/rho); end do
    energy = u(np)/rho - 0.5d0*vmag
    diff_en = 1d3; it = 1
    do while ( (abs(diff_en/energy) >= 1d-6) .and. (it <= 20) )
      T_i = int(T_it); if (T_i < tmin_) T_i = tmin_; if (T_i > tmax_-1) T_i = tmax_-1; Td = T_it - T_i
      cv = -Rg; enth = 0d0
      do s = 1, nsc
        cp_s = cp_tab(T_i,s) + Td*( cp_tab(T_i+1,s) - cp_tab(T_i,s) )
        h_s  = h_tab(T_i,s)  + Td*( h_tab(T_i+1,s)  - h_tab(T_i,s)  )
        cv   = cv   + u(s)/rho * cp_s
        enth = enth + u(s)/rho * h_s
      end do
      p_it    = rho*Rg*T_it
      diff_en = abs( (enth - p_it/rho) - energy )
      T_it    = T_it + (energy - (enth - p_it/rho))/cv
      it      = it + 1
    end do
    pout(np) = rho*Rg*T_it
    if (nprim>np) pout(np+1:nprim) = u(np+1:nprim)
  end subroutine mr_cons_interp

  !> mr-nonlocal = redo (open-hydra MOSE feature/multirate-integrator, the attempt loop, Mod_Multirate.f90 663-771 there): the
  !> macro step is run from a copy of the start state; if the certificate saw a violation of condition (4), the state is put
  !> back, dtcap is lowered to the smallest admissible step seen, the map is rebuilt for the step of the attempt (levels only
  !> rise) and H is halved while it is above the cap of the new map; a state that is not physical at the end of the attempt
  !> halves H with the same map. Stops below 2**-10 of the first H or after 21 rejected attempts. A halved H loses the landing
  !> on an output instant (F1). MOSE_MR_FORCE_REDO=1 redoes the first accepted attempt once with the same H and map: the result
  !> must not change (the test of the copy and the restore).
  subroutine MR_Macro_Step_Redo ( domain, config, External_Function, time_old )
    use MOSE_Types_m
    use MOSE_Field_m, only: P_flat, PO_flat
    use MOSE_Mod_MPI, only: mpi_is_root, mpi_allreduce_sum_r8_array
    use MOSE_Mod_Multirate_Partition, only: MR_Rebuild_For_H, MR_Set_Local_dt
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    type(config_type),      intent(inout) :: config
    external :: External_Function
    real(R8), intent(in) :: time_old
    integer, parameter :: MAX_REJ = 21
    logical, save :: force_checked = .false., force_on = .false.
    character(len=8) :: env
    integer  :: nrej, ln, st
    real(R8) :: H, H_first, hcap, cnt(2)
    logical  :: forced_done, nonphys, deeper, landed
    if ( .not. force_checked ) then
      force_checked = .true.
      call get_environment_variable ( 'MOSE_MR_FORCE_REDO', env, ln, st )
      force_on = ( st == 0 .and. ln > 0 .and. env(1:1) == '1' )
      if ( force_on .and. mpi_is_root ) write(*,'(A)') ' [MR] MOSE_MR_FORCE_REDO=1: every macro step is redone once'
    end if
    if ( .not. allocated(P0_flat) ) then
      allocate ( P0_flat(size(P_flat,1), size(P_flat,2)) )
      !$acc enter data create(P0_flat)
    end if
    !$acc wait
    call mr_field_copy ( P_flat, P0_flat )
    H_first = domain % dtglobal; H = H_first; nrej = 0; forced_done = .false.; landed = mr_landed
    attempts: do
      obj_mr%nonphys = .false.
      call MR_Macro_Step ( domain, config, External_Function, time_old )
      !$acc wait
      cnt(1) = merge(1.0_R8, 0.0_R8, obj_mr%nonphys); cnt(2) = 0.0_R8
      call mpi_allreduce_sum_r8_array ( cnt, 2 )
      nonphys = ( cnt(1) > 0.5_R8 )
      if ( obj_mr%viol_step == 0 .and. .not. nonphys ) then
        if ( .not. ( force_on .and. .not. forced_done ) ) exit attempts
        forced_done = .true.                       ! the forced redo: same H, same map
      else
        nrej = nrej + 1
        if ( nonphys ) then
          H = 0.5_R8 * H                           ! same map, half the step
        else
          call MR_Rebuild_For_H ( domain, H, hcap, deeper )
          if ( deeper ) obj_mr%nredo_deeper = obj_mr%nredo_deeper + 1
          do while ( H > hcap * ( 1.0_R8 + 1.0e-12_R8 ) )
            H = 0.5_R8 * H
          end do
        end if
      end if
      if ( H < H_first * 2.0_R8**(-10) ) error stop 'multirate redo: the macro step fell below 2**-10 of its first value'
      if ( nrej >= MAX_REJ ) error stop 'multirate redo: 21 rejected attempts in one macro step'
      call mr_field_copy ( P0_flat, P_flat )       ! back to the start state
      call mr_field_copy ( P0_flat, PO_flat )
      domain % dtglobal = H
      call MR_Set_Local_dt ( domain )
      if ( landed .and. time_old + H < mr_t_land - 1.0e-12_R8 * H ) then   ! a halved H no longer lands (F1)
        landed = .false.; mr_landed = .false.; mr_output_landed = .false.
      end if
      domain % time = time_old + H
      if ( landed ) domain % time = mr_t_land
    end do attempts
    if ( nrej > 0 ) then
      obj_mr%nredo_steps = obj_mr%nredo_steps + 1
      obj_mr%nredo_att = obj_mr%nredo_att + nrej
    end if
    if ( mpi_is_root .and. obj_mr%nprint > 0 ) then
      if ( mod(domain%iter, obj_mr%nprint) == 0 ) write(*,'(A,I0,A,I0,A,I0)') ' [MR]   redo (mr-nonlocal): macro steps ' // &
        'redone ', obj_mr%nredo_steps, '  rejected attempts ', obj_mr%nredo_att, '  deeper maps ', obj_mr%nredo_deeper
    end if
  end subroutine MR_Macro_Step_Redo

  subroutine mr_field_copy ( a, b )
    implicit none
    real(R8), intent(in)    :: a(:,:)
    real(R8), intent(inout) :: b(:,:)
    integer :: i, v
    !$acc parallel loop collapse(2) gang vector present(a, b)
    do v = 1, size(a,2)
      do i = 1, size(a,1)
        b(i,v) = a(i,v)
      end do
    end do
  end subroutine mr_field_copy


  !> mr-nonlocal (open-hydra MOSE feature/multirate-integrator, Cert_Stage, Mod_Multirate.f90 1047-1119 there): condition (4)
  !> of Gnedin, Semenov and Kravtsov (2018): at the start tau of a sub-step of a cell of level l, its step h^ = h_l (2 h_l if a
  !> face neighbour of level l-1 in its block starts a sub-step at tau) must not exceed its admissible step at that state,
  !> recomputed by Compute_dt_cert (the rule of Compute_dt); dmin keeps the smallest admissible step seen in the macro step.
  !> Counted on the device; reduced over the ranks once per macro step (mr_cert_close). Reads the field, changes nothing.
  subroutine mr_certificate ( domain, config, l, k )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Mod_dt,  only: Compute_dt_cert
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    type(config_type),      intent(inout) :: config
    integer,     intent(in) :: l
    integer(I8), intent(in) :: k
    integer  :: b, nv
    logical  :: cstart
    cstart = .false.
    if ( l > 0 ) cstart = ( tau_tick(l-1) == k )
    !$acc wait
    nv = 0
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      if (.not. mr_block_active(b)) cycle
      call Compute_dt_cert ( domain, config, b, mrb(b)%crho, mrb(b)%crg, mrb(b)%ct, mrb(b)%cmil, mrb(b)%dtc )
      call mr_cert_blk ( mrb(b)%lev, mrb(b)%dtc, mrb(b)%dmin, domain%blk(b)%dim, l, obj_mr%h_lev(l), cstart, nv )
    end do
    obj_mr%viol_step = obj_mr%viol_step + nv
  end subroutine mr_certificate

  subroutine mr_cert_blk ( lev, dtc, dmin, n, l, hl, cstart, nv )
    implicit none
    integer,  intent(in)    :: n(3), l, lev(n(1),n(2),n(3))
    real(R8), intent(in)    :: dtc(n(1),n(2),n(3)), hl
    real(R8), intent(inout) :: dmin(n(1),n(2),n(3))
    logical,  intent(in)    :: cstart
    integer,  intent(inout) :: nv
    integer  :: i, j, k, n1, n2, n3, nv_, l_, lm
    real(R8) :: hh, hl_
    logical  :: cs_, nb
    n1 = n(1); n2 = n(2); n3 = n(3); nv_ = 0; l_ = l; lm = l - 1; hl_ = hl; cs_ = cstart
    !$acc parallel loop collapse(3) gang vector present(lev, dtc, dmin) reduction(+:nv_) private(hh, nb)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      if ( lev(i,j,k) /= l_ ) cycle
      hh = hl_
      if ( cs_ ) then
        nb = .false.
        if ( i > 1  ) nb = nb .or. ( lev(i-1,j,k) == lm )
        if ( i < n1 ) nb = nb .or. ( lev(i+1,j,k) == lm )
        if ( j > 1  ) nb = nb .or. ( lev(i,j-1,k) == lm )
        if ( j < n2 ) nb = nb .or. ( lev(i,j+1,k) == lm )
        if ( k > 1  ) nb = nb .or. ( lev(i,j,k-1) == lm )
        if ( k < n3 ) nb = nb .or. ( lev(i,j,k+1) == lm )
        if ( nb ) hh = 2.0_R8 * hl_
      end if
      dmin(i,j,k) = min( dmin(i,j,k), dtc(i,j,k) )
      if ( hh > dtc(i,j,k) * ( 1.0_R8 + 1.0e-12_R8 ) ) nv_ = nv_ + 1
    end do
    end do
    end do
    nv = nv + nv_
  end subroutine mr_cert_blk

  !> mr-nonlocal: the start of the macro step: dmin := dtcap (the admissible steps at t_n), no violation yet.
  subroutine mr_cert_open ( domain )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: b, i, j, k, n1, n2, n3
    obj_mr%viol_step = 0
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      n1 = domain%blk(b)%dim(1); n2 = domain%blk(b)%dim(2); n3 = domain%blk(b)%dim(3)
      associate ( dm => mrb(b)%dmin, dc => mrb(b)%dtcap )
      !$acc parallel loop collapse(3) gang vector present(dm, dc)
      do k = 1, n3
      do j = 1, n2
      do i = 1, n1
        dm(i,j,k) = dc(i,j,k)
      end do
      end do
      end do
      end associate
    end do
  end subroutine mr_cert_open

  !> mr-nonlocal: the violations of the macro step summed over the ranks, into the run's counters.
  subroutine mr_cert_close ()
    use MOSE_Mod_MPI, only: mpi_allreduce_sum_r8_array
    implicit none
    real(R8) :: cnt(2)
    cnt(1) = real(obj_mr%viol_step, R8); cnt(2) = 0.0_R8
    call mpi_allreduce_sum_r8_array ( cnt, 2 )
    obj_mr%viol_step = nint(cnt(1))
    if ( obj_mr%viol_step > 0 ) then
      obj_mr%ncert_viol = obj_mr%ncert_viol + obj_mr%viol_step
      obj_mr%ncert_steps = obj_mr%ncert_steps + 1
    end if
  end subroutine mr_cert_close

end module MOSE_Mod_Multirate
