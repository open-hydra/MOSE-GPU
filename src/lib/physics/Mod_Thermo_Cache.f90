module MOSE_Mod_Thermo_Cache
  ! Per-cell thermo scratchpad, refreshed once per RK stage and once before
  ! Compute_dt. Replaces redundant co_rotot_Rtot / T = p/(rho R) /
  ! f_laminarViscosity calls scattered across Mod_Diffusive, Mod_dt, Mod_SST,
  ! Lib_Newstate. Cache is filled including ghost cells so interface kernels
  ! can read the four neighbours without bounds care.
  use iso_fortran_env, only: I4 => int32, R8 => real64

  use MOSE_Global_m,   only: gc, nprim, nsc, np, nu
  use MOSE_Field_m,    only: P_flat, field_block_index_params
  implicit none
  private
  public :: Refresh_Thermo_Cache, Refresh_Thermo_Cache_Block

  public :: Thermo_Cache_rrt_blk, Thermo_Cache_visc_blk
  public :: Thermo_Cache_fused_blk, Thermo_Cache_rrt_shell_blk, thermo_guard_report, thermo_guard_on
  integer  :: tg_nbad = 0, tg_c0 = 0, tg_b = 0
  real(R8) :: tg_T = 0.0_R8, tg_rho = 0.0_R8, tg_p = 0.0_R8
  !$acc declare create(tg_nbad, tg_c0, tg_b, tg_T, tg_rho, tg_p)
  logical, save :: thermo_guard_on = .false., tg_checked = .false.

contains

  subroutine Refresh_Thermo_Cache ( domain, shell_only )
    ! Domain-wide refresh. Caller must already be inside an !$omp parallel.
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Global_m, only: gc, VISC

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Shell_only skips the INTERIOR pass. The caller
    ! asserts the interior cache is already current and bit-identical (substage (1,1):
    ! Compute_dt just ran the split rrt+visc interior kernels on the SAME w, and the
    ! fused interior kernel is gate-proven bit-equal to that pair); only the
    ! halo-refreshed GHOST shell needs the rrt pass. Ignored on the host path (full
    ! refresh recomputes identical values — correct, just no saving).
    logical, intent(in), optional :: shell_only
    integer :: b
    logical :: shell_l, guard_once

    shell_l = .false.; guard_once = .false.
    if (present(shell_only)) shell_l = shell_only

    ! Per-stage thermo cache on the var-LAST flat SoA w. w is kept
    ! resident + current per substage in Explicit_Step BEFORE this cache call (n=1:
    ! Compute_dt step-start pack + Newstate + Fill_Ghost_Cell_local; n>1: field_sync_P_to_flat)
    ! -- so NO pack here. rrt covers interior+ghosts (matches the P-path); visc is
    ! interior-only (ghost prims can be vacuum -> OOB Wilke int(T)).
    ! Tier-1.3 sweep fusion (VISC): ONE interior pass does rrt+Wilke+cp from a single
    ! species gather, then a cheap ghost-shell rrt pass covers the padding -- replaces
    ! the rrt(padded)+visc(interior) pair, dropping a full redundant read of w per
    ! substage. Bit-exact: same cells, same algebra, cache stores are IEEE round-trips.
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      if (VISC) then
        if (.not. shell_l) &
        call Thermo_Cache_fused_blk ( b, domain%blk(b)%dim, domain%blk(b)%rho_c,      &
               domain%blk(b)%Rg_c, domain%blk(b)%T_c, domain%blk(b)%mil_c )
        call thermo_guard_report ( b, domain%blk(b)%dim, domain%iter )
        call Thermo_Cache_rrt_shell_blk ( b, domain%blk(b)%dim, domain%blk(b)%rho_c,  &
               domain%blk(b)%Rg_c, domain%blk(b)%T_c )
      else if (shell_l) then
        call Thermo_Cache_rrt_shell_blk ( b, domain%blk(b)%dim, domain%blk(b)%rho_c,  &
               domain%blk(b)%Rg_c, domain%blk(b)%T_c )
      else
        call Thermo_Cache_rrt_blk ( b, domain%blk(b)%dim,                             &
               1-gc, domain%blk(b)%dim(1)+gc, 1-gc, domain%blk(b)%dim(2)+gc,            &
               1-gc, domain%blk(b)%dim(3)+gc,                                           &
               domain%blk(b)%rho_c, domain%blk(b)%Rg_c, domain%blk(b)%T_c )
      end if
    enddo
    if ( guard_once ) call thermo_guard_report ( -1, domain%blk(1)%dim, domain%iter )   ! 5i: one device read per list stage
    return
  end subroutine Refresh_Thermo_Cache


  subroutine Refresh_Thermo_Cache_Block ( blk, interior_only )
    ! Single-block refresh. Caller must already be inside an !$omp parallel.
    ! interior_only=.true. skips ghost cells — used by Compute_dt, which runs
    ! before any Fill_Ghost_Cell and therefore cannot trust ghost prims.
    ! f_laminarViscosity dispatches on int(T) into a tabulated array, so feeding
    ! it garbage temperatures from uninitialised ghosts causes OOB segfaults.
    use MOSE_Types_m
    use MOSE_Global_m, only: nsc, np, gc, VISC
    use FLINT_Lib_Thermodynamic, only: co_rotot_Rtot, co_k_mi_lam_Wilke_expr, f_cp_expr
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, co_k_mi_lam_Wilke_dev, f_cp_dev

    implicit none
    type(MOSE_block_type), intent(inout) :: blk
    logical, intent(in), optional        :: interior_only
    integer  :: i, j, k, n1, n2, n3, lo1, hi1, lo2, hi2, lo3, hi3
    integer  :: T_i, Tint(2)
    real(R8) :: rho, Rg, T, Tdiff, mil, kl
    logical  :: int_only

    int_only = .false.
    if (present(interior_only)) int_only = interior_only

    n1 = blk % dim(1) ; n2 = blk % dim(2) ; n3 = blk % dim(3)
    if (int_only) then
      lo1 = 1 ; hi1 = n1 ; lo2 = 1 ; hi2 = n2 ; lo3 = 1 ; hi3 = n3
    else
      lo1 = 1-gc ; hi1 = n1+gc ; lo2 = 1-gc ; hi2 = n2+gc ; lo3 = 1-gc ; hi3 = n3+gc
    end if

    !> Device path: Per-cell thermo on the GPU. The acc loops live
    !> in helper subroutines that take PLAIN ARRAY dummies (blk%P etc. passed as
    !> actual args) — accessing derived-type allocatable members directly in an
    !> acc region needs manual deep-copy/attach and faults (ILLEGAL_ADDRESS);
    !> the chemistry kernel works for the same reason (plain-array dummy).
    !> P copyin / caches copy per call (correctness; residency removes the
    !> transfer later). omp single: one thread launches inside an OpenMP team.
    !$omp single
    call Thermo_Cache_rrt ( blk%dim(1), blk%dim(2), blk%dim(3), &
                                lo1, hi1, lo2, hi2, lo3, hi3, &
                                blk%P, blk%rho_c, blk%Rg_c, blk%T_c )
    if (VISC) &
      call Thermo_Cache_visc ( blk%dim(1), blk%dim(2), blk%dim(3), &
                                   blk%P, blk%rho_c, blk%T_c, &
                                   blk%mil_c, blk%kl_c, blk%cp_c )
    !$omp end single
    return
  end subroutine Refresh_Thermo_Cache_Block


  !> Device kernel: rho_c/Rg_c/T_c over [lo:hi]. Plain-array dummies (explicit
  !> block bounds) so acc copyin(P)/copy(caches) work without derived-type
  !> deep-copy. caches use copy() not copyout so cells outside [lo:hi] (ghosts
  !> in the interior_only call) keep their host values, matching the CPU path.
  subroutine Thermo_Cache_rrt ( n1, n2, n3, lo1, hi1, lo2, hi2, lo3, hi3, &
                                    P, rho_c, Rg_c, T_c )
    use MOSE_Global_m, only: nprim, nsc, np, gc
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev
    implicit none
    integer,  intent(in)    :: n1, n2, n3, lo1, hi1, lo2, hi2, lo3, hi3
    real(R8), intent(in)    :: P    (nprim, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    real(R8), intent(inout) :: rho_c(       1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    real(R8), intent(inout) :: Rg_c (       1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    real(R8), intent(inout) :: T_c  (       1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    integer  :: i, j, k
    real(R8) :: rho, Rg

    !$acc data copyin(P) copy(rho_c, Rg_c, T_c)
    !$acc parallel loop collapse(3) private(rho, Rg)
    do k = lo3, hi3
    do j = lo2, hi2
    do i = lo1, hi1
      call co_rotot_Rtot_dev ( P(1:nsc,i,j,k), rho, Rg )
      rho_c(i,j,k) = rho
      Rg_c (i,j,k) = Rg
      T_c  (i,j,k) = P(np,i,j,k) / ( rho * Rg )
    enddo; enddo; enddo
    !$acc end data
  end subroutine Thermo_Cache_rrt


  !> Device kernel: transport (mil_c/kl_c) + mixture cp (cp_c) at interior cells.
  subroutine Thermo_Cache_visc ( n1, n2, n3, P, rho_c, T_c, mil_c, kl_c, cp_c )
    use MOSE_Global_m, only: nprim, nsc, gc
    use FLINT_Lib_ThermoTransport_dev, only: co_k_mi_lam_Wilke_dev, f_cp_dev
    implicit none
    integer,  intent(in)    :: n1, n2, n3
    real(R8), intent(in)    :: P    (nprim, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    real(R8), intent(in)    :: rho_c(       1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    real(R8), intent(in)    :: T_c  (       1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    real(R8), intent(inout) :: mil_c(       1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    real(R8), intent(inout) :: kl_c (       1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    real(R8), intent(inout) :: cp_c (       1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    integer  :: i, j, k, T_i, Tint(2)
    real(R8) :: T, Tdiff, mil, kl

    !$acc data copyin(P, rho_c, T_c) copy(mil_c, kl_c, cp_c)
    !$acc parallel loop collapse(3) private(T, T_i, Tdiff, Tint, mil, kl)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      T = T_c(i,j,k)
      T_i = idint(T)
      Tdiff   = T - T_i
      Tint(1) = T_i
      Tint(2) = T_i + 1
      call co_k_mi_lam_Wilke_dev ( P(1:nsc,i,j,k), rho_c(i,j,k), Tint, Tdiff, mil, kl )
      mil_c(i,j,k) = mil
      kl_c (i,j,k) = kl
      cp_c (i,j,k) = f_cp_dev ( P(1:nsc,i,j,k), Tint, Tdiff, rho_c(i,j,k) )
    enddo; enddo; enddo
    !$acc end data
  end subroutine Thermo_Cache_visc


  ! ---- from Mod_Thermo_Cache_W.f90 (MOSE_Mod_Thermo_Cache) ----

  !> FUSED interior sweep: rho/Rg/T (rrt) AND transport+cp (visc) from ONE gather of
  !> the species row. Bit-exact vs the split kernels: T_c/rho_c are IEEE doubles, so
  !> using the just-computed rho/T locally is bitwise the store+reload the split visc
  !> kernel performs; Wilke/cp receive identical inputs in identical order.
  subroutine Thermo_Cache_fused_blk ( b, n, rho_c, Rg_c, T_c, mil_c)
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, co_mi_lam_Wilke_dev
    use MOSE_Field_m, only: rho_c_f, T_c_f, Rg_c_f, mil_c_f, n_cache_f
    integer,  intent(in) :: b, n(3)
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: rho_c, Rg_c, T_c
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: mil_c
    integer  :: i, j, k, s, c0, boff, pi, pj, T_i, Tint(2)
    real(R8) :: rho, Rg, T, Tdiff, mil, rl(nsc)
    ! Dims as scalars -- a dummy array referenced inside the acc region
    ! is implicitly copy()d per launch (upload+download+sync); scalars ride the
    ! launch parameters for free.
    integer :: ib
    logical :: guard_
    integer :: nb1_, nb2_, nb3_
    nb1_ = n(1); nb2_ = n(2); nb3_ = n(3)
    call field_block_index_params(b, boff, pi, pj)

    !$acc data copy(rho_c, Rg_c, T_c, mil_c)
    guard_ = thermo_guard_on
    !$acc parallel loop collapse(3) gang vector present(P_flat, rho_c, Rg_c, T_c, mil_c, rho_c_f, T_c_f, Rg_c_f, mil_c_f) &
    !$acc&  private(rl, rho, Rg, T, T_i, Tdiff, Tint, mil, c0, s, ib)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)
      c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      do s = 1, nsc
        rl(s) = P_flat(c0,s)
      end do
      call co_rotot_Rtot_dev ( rl(1:nsc), rho, Rg )
      T = P_flat(c0,np) / ( rho * Rg )
      if ( guard_ ) then
        if ( .not. ( T > 1.0_R8 .and. T < 9990.0_R8 ) ) then   ! NaN fails both tests
          !$acc atomic capture
          tg_nbad = tg_nbad + 1
          ib = tg_nbad
          !$acc end atomic
          if ( ib == 1 ) then
            tg_c0 = c0; tg_b = b; tg_T = T; tg_rho = rho; tg_p = P_flat(c0,np)
          end if
          if ( T /= T ) T = 300.0_R8
          T = min( max(T, 1.0_R8), 9990.0_R8 )
        end if
      end if
      rho_c(i,j,k) = rho
      if ( n_cache_f > 1 ) then; rho_c_f(c0) = rho; T_c_f(c0) = T; Rg_c_f(c0) = Rg; end if   ! row-addressed twin (+ Rg)
      Rg_c (i,j,k) = Rg
      T_c  (i,j,k) = T
      T_i = idint(T); Tdiff = T - T_i
      Tint(1) = T_i; Tint(2) = T_i + 1
      ! face-exact production: the cell conductivity/cp caches (kl_c/cp_c) are dead
      ! (the diffusive flux recomputes transport at faces from face-averaged Prim), so
      ! only mil_c is produced, via the mil-only Wilke twin (milam bit-identical).
      call co_mi_lam_Wilke_dev ( rl(1:nsc), rho, Tint, Tdiff, mil )
      mil_c(i,j,k) = mil
      if ( n_cache_f > 1 ) mil_c_f(c0) = mil   ! the flat twin
    end do; end do; end do
    !$acc end data
  end subroutine Thermo_Cache_fused_blk


  !> Ghost-SHELL rrt: rho/Rg/T over the padded range MINUS the interior (the fused
  !> kernel covered that). Interior threads fall through -- rrt is cheap and the shell
  !> is ~15% of cells, so the idle lanes cost less than a second full species gather.
  subroutine Thermo_Cache_rrt_shell_blk ( b, n, rho_c, Rg_c, T_c )
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev
    use MOSE_Field_m, only: rho_c_f, T_c_f, Rg_c_f, mil_c_f, n_cache_f
    integer,  intent(in) :: b, n(3)
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: rho_c, Rg_c, T_c
    integer  :: i, j, k, s, c0, boff, pi, pj
    real(R8) :: rho, Rg, rl(nsc)
    ! Dims as scalars -- a dummy array referenced inside the acc region
    ! is implicitly copy()d per launch (upload+download+sync); scalars ride the
    ! launch parameters for free.
    integer :: nb1_, nb2_, nb3_
    nb1_ = n(1); nb2_ = n(2); nb3_ = n(3)
    call field_block_index_params(b, boff, pi, pj)

    !$acc data copy(rho_c, Rg_c, T_c)
    !$acc parallel loop collapse(3) gang vector present(P_flat, rho_c, Rg_c, T_c, rho_c_f, T_c_f, Rg_c_f) &
    !$acc&  private(rl, rho, Rg, c0, s)
    do k = 1-gc, nb3_+gc
    do j = 1-gc, nb2_+gc
    do i = 1-gc, nb1_+gc
      if ( i >= 1 .and. i <= nb1_ .and. j >= 1 .and. j <= nb2_ .and. &
           k >= 1 .and. k <= nb3_ ) cycle
      c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      do s = 1, nsc
        rl(s) = P_flat(c0,s)
      end do
      call co_rotot_Rtot_dev ( rl(1:nsc), rho, Rg )
      rho_c(i,j,k) = rho
      Rg_c (i,j,k) = Rg
      T_c  (i,j,k) = P_flat(c0,np) / ( rho * Rg )
      if ( n_cache_f > 1 ) then   ! row-addressed twin, read by the fused diffusive flux
        rho_c_f(c0) = rho; T_c_f(c0) = P_flat(c0,np) / ( rho * Rg ); Rg_c_f(c0) = Rg   ! (+ Rg)
      end if
    end do; end do; end do
    !$acc end data
  end subroutine Thermo_Cache_rrt_shell_blk




  !> rho_c/Rg_c/T_c over [lo:hi] for block b, reading the var-last `w` primitives.
  subroutine Thermo_Cache_rrt_blk ( b, n, lo1, hi1, lo2, hi2, lo3, hi3, rho_c, Rg_c, T_c )
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev
    use MOSE_Field_m, only: rho_c_f, T_c_f, Rg_c_f, mil_c_f, n_cache_f
    integer,  intent(in) :: b, n(3), lo1, hi1, lo2, hi2, lo3, hi3
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: rho_c, Rg_c, T_c
    integer  :: i, j, k, s, c0, boff, pi, pj
    real(R8) :: rho, Rg, rl(nsc)
    call field_block_index_params(b, boff, pi, pj)

    !$acc data copy(rho_c, Rg_c, T_c)
    !$acc parallel loop collapse(3) gang vector present(P_flat, rho_c, Rg_c, T_c, rho_c_f, T_c_f, Rg_c_f) &
    !$acc&  private(rl, rho, Rg, c0, s)
    do k = lo3, hi3
    do j = lo2, hi2
    do i = lo1, hi1
      c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      do s = 1, nsc
        rl(s) = P_flat(c0,s)
      end do
      call co_rotot_Rtot_dev ( rl(1:nsc), rho, Rg )
      rho_c(i,j,k) = rho
      Rg_c (i,j,k) = Rg
      T_c  (i,j,k) = P_flat(c0,np) / ( rho * Rg )
      if ( n_cache_f > 1 ) then   ! row-addressed twin, read by the fused diffusive flux
        rho_c_f(c0) = rho; T_c_f(c0) = P_flat(c0,np) / ( rho * Rg ); Rg_c_f(c0) = Rg   ! (+ Rg)
      end if
    end do; end do; end do
    !$acc end data
  end subroutine Thermo_Cache_rrt_blk

  !> Transport (mil_c) at INTERIOR cells for block b, reading the var-last `w` species
  !> + the rho_c/T_c just computed by the rrt kernel. face-exact: only mil_c (dt reads
  !> it); the conductivity/cp caches are dead (recomputed at faces).
  subroutine Thermo_Cache_visc_blk ( b, n, rho_c, T_c, mil_c )
    use FLINT_Lib_ThermoTransport_dev, only: co_mi_lam_Wilke_dev
    use MOSE_Field_m, only: mil_c_f, n_cache_f   ! every writer of mil_c mirrors into the twin
    integer,  intent(in) :: b, n(3)
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in)    :: rho_c, T_c
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: mil_c
    integer  :: i, j, k, s, c0, boff, pi, pj, T_i, Tint(2)
    real(R8) :: T, Tdiff, mil, rl(nsc)
    call field_block_index_params(b, boff, pi, pj)

    !$acc data copyin(rho_c, T_c) copy(mil_c)
    !$acc parallel loop collapse(3) gang vector present(P_flat, rho_c, T_c, mil_c, mil_c_f) &
    !$acc&  private(rl, T, T_i, Tdiff, Tint, mil, c0, s)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)
      c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      do s = 1, nsc
        rl(s) = P_flat(c0,s)
      end do
      T = T_c(i,j,k)
      T_i = idint(T); Tdiff = T - T_i
      Tint(1) = T_i; Tint(2) = T_i + 1
      call co_mi_lam_Wilke_dev ( rl(1:nsc), rho_c(i,j,k), Tint, Tdiff, mil )
      mil_c(i,j,k) = mil
      if ( n_cache_f > 1 ) mil_c_f(c0) = mil   ! the flat twin
    end do; end do; end do
    !$acc end data
  end subroutine Thermo_Cache_visc_blk
  !> Host side of the guard: enable once from MOSE_THERMO_GUARD, and after a block's refresh
  !> report the first off-table cell with its full primitive state, then stop.
  subroutine thermo_guard_report ( b, n, iter )
    use MOSE_Mod_MPI, only: mpi_rank_
    implicit none
    integer, intent(in) :: b, n(3), iter
    integer :: boff, pi, pj, c0, i, j, k, rem, st, s_, q, cq, di, dj, dk
    character(len=8) :: env
    character(len=32) :: wenv
    real(R8) :: wl(nprim), rhoq
    integer, save :: wb = 0, wi = 0, wj = 0, wk = 0, wfrom = 0
    logical, save :: wchecked = .false., won = .false.
    integer, parameter :: dq(3,7) = reshape([0,0,0, -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1], [3,7])
    if ( .not. wchecked ) then
      call get_environment_variable('MOSE_WATCH', wenv, status=st)
      if ( st == 0 ) then
        read(wenv, *, iostat=st) wb, wi, wj, wk
        won = ( st == 0 .and. wb > 0 )
      end if
      call get_environment_variable('MOSE_WATCH_FROM', wenv, status=st)
      if ( st == 0 ) read(wenv, *, iostat=st) wfrom
      wchecked = .true.
    end if
    if ( won .and. b == wb .and. iter >= wfrom .and. mpi_rank_ >= 0 ) then
      call field_block_index_params(b, boff, pi, pj)
      do q = 1, 7
        di = wi + dq(1,q); dj = wj + dq(2,q); dk = wk + dq(3,q)
        cq = boff + (di+gc-1) + (dj+gc-1)*pi + (dk+gc-1)*pi*pj + 1
        !$acc update host(P_flat(cq:cq, 1:nprim))
        wl = P_flat(cq, 1:nprim); rhoq = sum(wl(1:nsc))
        write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES11.4,A,3ES11.3,A,ES12.4,A,ES10.3,A,ES10.3)') &
          ' [WATCH] it ', iter, ' b', b, ' (', di, ',', dj, ',', dk, ')  rho ', rhoq, '  u,v,w ', wl(nu), wl(nu+1), wl(nu+2), &
          '  p ', wl(np), '  rk ', wl(np+1), '  rw ', wl(np+2)
      end do
    end if
    if ( .not. tg_checked ) then
      call get_environment_variable('MOSE_THERMO_GUARD', env, status=st)
      thermo_guard_on = .true.     ! on by default (it only checks and stops); MOSE_THERMO_GUARD=0 or off disables it
      if ( st == 0 ) then
        if ( trim(env) == '0' .or. trim(env) == 'off' ) thermo_guard_on = .false.
      end if
      tg_checked = .true.
    end if
    if ( .not. thermo_guard_on ) return
    !$acc update host(tg_nbad)
    if ( tg_nbad == 0 ) return
    call thermo_guard_dump ()
  end subroutine thermo_guard_report

  !> The report itself.
  subroutine thermo_guard_dump ()
    use MOSE_Mod_MPI, only: mpi_rank_
    implicit none
    integer :: boff, pi, pj, c0, i, j, k, rem, s_
    real(R8) :: wl(nprim)
    !$acc update host(tg_c0, tg_b, tg_T, tg_rho, tg_p)
    call field_block_index_params(tg_b, boff, pi, pj)
    c0 = tg_c0 - boff - 1; k = c0 / (pi*pj); rem = c0 - k*pi*pj; j = rem / pi; i = rem - j*pi
    i = i - gc + 1; j = j - gc + 1; k = k - gc + 1
    !$acc update host(P_flat(tg_c0:tg_c0, 1:nprim))
    wl = P_flat(tg_c0, 1:nprim)
    write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,ES12.4,A,ES12.4,A,ES12.4)') &
      ' [THERMO-GUARD] rank ', mpi_rank_, ' block ', tg_b, ' cell (', i, ',', j, ',', k, ') offenders ', tg_nbad, &
      '  T ', tg_T, '  rho ', tg_rho, '  p ', tg_p
    write(*,'(A)') ' [THERMO-GUARD] primitive state P_flat(1:nprim) (species rho_s, u, v, w, p, rho*k, rho*omega, ...):'
    do s_ = 1, nprim
      write(*,'(A,I0,A,ES14.6)') '   P_flat(', s_, ') = ', wl(s_)
    end do
    error stop 'THERMO-GUARD: off-table temperature -- state printed above'
  end subroutine thermo_guard_dump

end module MOSE_Mod_Thermo_Cache
