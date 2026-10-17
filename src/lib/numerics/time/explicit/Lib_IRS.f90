!>@brief: Implicit Residual Smoothing related subroutines.
module MOSE_Lib_IRS

  use iso_fortran_env, only: R8 => real64
  use MOSE_Types_m
  use MOSE_Global_m, only: gc, nprim, np, ndir
  use MOSE_Field_m,  only: R_flat, ncell_total, field_block_index_params
  implicit none

  type :: irs_type
    integer :: njacobi = 2
    real(8) :: beta
    logical :: IRS
    contains
      procedure, pass(self) :: ini_setup => read_ini
  end type

  type(irs_type), public :: obj_irs
  private :: read_ini

  public :: IRS_apply
  real(R8), allocatable, save :: rs1(:,:), rs2(:,:)   ! Res_Star / Res_Star_New: var-last flat, like R_flat
  logical, save :: irs_ready = .false.

contains


  subroutine read_ini(self, fini)
    use Finer, only: file_ini
    use MOSE_Parameters_m, only: codename
    
    implicit none
    class(irs_type) :: self
    type(file_ini) :: fini
    ! Local
    integer :: error

    call fini%get(section_name=trim(codename)//'-scheme', option_name='IRS', val=self%IRS, error=error)
    if (error/=0) self%IRS = .false.
    call fini%get(section_name=trim(codename)//'-scheme', option_name='beta', val=self%beta, error=error)
    if (error/=0) self%beta = 0.5d0

  end subroutine read_ini


  !% Compute the Implicit Residual Smooting operator: R* = R + beta*LAPLACE(R*) with Jacobi iterations.



  ! ---- from Mod_IRS_W.f90 (MOSE_Lib_IRS) ----

  !> MOSE_open Newstate_Conservative, IRS branch, phases 1 and 2.
  subroutine IRS_apply ( domain, beta, njacobi, strang, scale )
    use MOSE_Mod_MPI, only: is_local_block, mpi_size_
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    real(R8), intent(in) :: beta, strang
    integer,  intent(in) :: njacobi
    logical,  intent(in), optional :: scale   ! .false.: smooth the raw residual (MOSE_open Newstate_Primitive)
    integer :: b, d, boff, pi, pj, off
    logical :: scale_

    if ( mpi_size_ > 1 ) error stop 'IRS on the GPU: one rank only (the inter-rank ghost residual exchange is not ported)'
    if ( .not. irs_ready ) then
      allocate( rs1(ncell_total, np), rs2(ncell_total, np) )
      rs1 = 0.0_R8; rs2 = 0.0_R8
      !$acc enter data copyin(rs1, rs2)
      irs_ready = .true.
    end if
    !$acc wait
    ! Phase 1: Scale_Residual_Cons on every interior cell (conservative path only). Two statements: Fortran does not
    ! short-circuit .or., so `.not. present(scale) .or. scale` reads the absent argument
    scale_ = .true.
    if ( present(scale) ) scale_ = scale
    if ( scale_ ) then
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        call irs_scale_blk ( b, domain % blk(b) % dtlocal, domain % blk(b) % vol, domain % blk(b) % dim, strang )
      end do
    end if
    ! Phase 2: Residual_Smoothing, direction by direction
    do d = 1, ndir
      call irs_ghost_residual ( domain )
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        call field_block_index_params(b, boff, pi, pj)
        off = 1
        if (d == 2) off = pi
        if (d == 3) off = pi*pj
        call irs_smooth_blk ( boff, pi, pj, domain%blk(b)%dim(1), domain%blk(b)%dim(2), domain%blk(b)%dim(3), &
                              off, beta, njacobi )
      end do
    end do
    !$acc wait
  end subroutine IRS_apply

  !> residual = -residual / volume * dt * strangcoeff (Scale_Residual_Cons), interior cells, all the variables.
  subroutine irs_scale_blk ( b, dtlocal, vol, n, strang )
    implicit none
    integer,  intent(in) :: b, n(3)
    real(R8), intent(in) :: dtlocal(1:n(1), 1:n(2), 1:n(3))
    real(R8), intent(in) :: vol(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc)
    real(R8), intent(in) :: strang
    integer  :: i, j, k, s, c0, boff, pi, pj, n1_, n2_, n3_
    real(R8) :: dtl_, vol_
    call field_block_index_params(b, boff, pi, pj)
    n1_ = n(1); n2_ = n(2); n3_ = n(3)
    !$acc data copyin(dtlocal, vol)
    !$acc parallel loop collapse(3) gang vector present(R_flat, dtlocal, vol) private(c0, s, dtl_, vol_)
    do k = 1, n3_
    do j = 1, n2_
    do i = 1, n1_
      c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      dtl_ = dtlocal(i,j,k); vol_ = vol(i,j,k)
      do s = 1, nprim
        R_flat(c0,s) = -R_flat(c0,s)/vol_*dtl_*strang
      end do
    end do; end do; end do
    !$acc end data
  end subroutine irs_scale_blk

  !> Ghost_Residual_Extrapolation: the first ghost layer of every block face takes the adjacent cell's
  !> residual (zero gradient: wall, inlet, outlet, symmetry, periodicity, no-flux faces), then every
  !> connection record overwrites its ghost with the donor cell's residual. Flow equations 1:np.
  subroutine irs_ghost_residual ( domain )
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer :: b, p, boff, pi, pj, n1, n2, n3, i, j, k, s, cg, cm
    integer :: boff_m, pi_m, pj_m, boff_s, pi_s, pj_s, Bm, Bs, lo, hi
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      call field_block_index_params(b, boff, pi, pj)
      n1 = domain%blk(b)%dim(1); n2 = domain%blk(b)%dim(2); n3 = domain%blk(b)%dim(3)
      !$acc parallel loop collapse(2) gang vector present(R_flat) private(cg, cm, s)
      do k = 1, n3
      do j = 1, n2
        cm = boff + (1+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1;  cg = cm - 1      ! face 1
        do s = 1, np; R_flat(cg,s) = R_flat(cm,s); end do
        cm = boff + (n1+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1; cg = cm + 1      ! face 2
        do s = 1, np; R_flat(cg,s) = R_flat(cm,s); end do
      end do; end do
      !$acc parallel loop collapse(2) gang vector present(R_flat) private(cg, cm, s)
      do k = 1, n3
      do i = 1, n1
        cm = boff + (i+gc-1) + (1+gc-1)*pi + (k+gc-1)*pi*pj + 1;  cg = cm - pi     ! face 3
        do s = 1, np; R_flat(cg,s) = R_flat(cm,s); end do
        cm = boff + (i+gc-1) + (n2+gc-1)*pi + (k+gc-1)*pi*pj + 1; cg = cm + pi     ! face 4
        do s = 1, np; R_flat(cg,s) = R_flat(cm,s); end do
      end do; end do
      !$acc parallel loop collapse(2) gang vector present(R_flat) private(cg, cm, s)
      do j = 1, n2
      do i = 1, n1
        cm = boff + (i+gc-1) + (j+gc-1)*pi + (1+gc-1)*pi*pj + 1;  cg = cm - pi*pj  ! face 5
        do s = 1, np; R_flat(cg,s) = R_flat(cm,s); end do
        cm = boff + (i+gc-1) + (j+gc-1)*pi + (n3+gc-1)*pi*pj + 1; cg = cm + pi*pj  ! face 6
        do s = 1, np; R_flat(cg,s) = R_flat(cm,s); end do
      end do; end do
    end do
    do p = 1, domain % n_conn_pair
      Bm = domain % conn_pair_bm(p); Bs = domain % conn_pair_bs(p)
      if (.not. (is_local_block(Bm) .and. is_local_block(Bs))) cycle
      call field_block_index_params(Bm, boff_m, pi_m, pj_m)
      call field_block_index_params(Bs, boff_s, pi_s, pj_s)
      lo = domain % conn_pair_lo(p); hi = domain % conn_pair_hi(p)
      if (hi < lo) cycle
      call irs_ghost_conn ( lo, hi, domain % conn_im, domain % conn_jm, domain % conn_km, domain % conn_fm, &
                            domain % conn_is, domain % conn_js, domain % conn_ks, domain % n_conn, &
                            boff_m, pi_m, pj_m, boff_s, pi_s, pj_s )
    end do
  end subroutine irs_ghost_residual

  !> r(bm)(ig,jg,kg) = r(bs)(is,js,ks): the ghost behind the receiver face cell (one layer) from the donor cell.
  subroutine irs_ghost_conn ( lo, hi, c_im, c_jm, c_km, c_fm, c_is, c_js, c_ks, nc, boff_m, pi_m, pj_m, boff_s, pi_s, pj_s )
    implicit none
    integer, intent(in) :: lo, hi, nc, boff_m, pi_m, pj_m, boff_s, pi_s, pj_s
    integer, intent(in) :: c_im(nc), c_jm(nc), c_km(nc), c_fm(nc), c_is(nc), c_js(nc), c_ks(nc)
    integer :: e, s, ig, jg, kg, fm, gstep, cg, cd
    !$acc parallel loop gang vector present(R_flat, c_im, c_jm, c_km, c_fm, c_is, c_js, c_ks) private(s, ig, jg, kg, fm, gstep, cg, cd)
    do e = lo, hi
      fm = c_fm(e); gstep = 2*mod(fm,2) - 1               ! +1 on the odd (low) faces, -1 on the even ones
      ig = c_im(e); jg = c_jm(e); kg = c_km(e)
      if ( fm <= 2 ) then; ig = ig - gstep; else if ( fm <= 4 ) then; jg = jg - gstep; else; kg = kg - gstep; end if
      cg = boff_m + (ig+gc-1) + (jg+gc-1)*pi_m + (kg+gc-1)*pi_m*pj_m + 1
      cd = boff_s + (c_is(e)+gc-1) + (c_js(e)+gc-1)*pi_s + (c_ks(e)+gc-1)*pi_s*pj_s + 1
      do s = 1, np
        R_flat(cg,s) = R_flat(cd,s)
      end do
    end do
  end subroutine irs_ghost_conn

  !> Residual_Smoothing_Blk in one direction: Res_Star = Res (ghosts included), njacobi Jacobi sweeps on the
  !> interior, Res(1:np) = Res_Star.
  subroutine irs_smooth_blk ( boff, pi, pj, n1, n2, n3, off, beta, njacobi )
    implicit none
    integer,  intent(in) :: boff, pi, pj, n1, n2, n3, off, njacobi
    real(R8), intent(in) :: beta
    integer :: i, j, k, s, c, c0, jac, nflat
    nflat = pi*pj*(n3 + 2*gc)
    !$acc parallel loop collapse(2) gang vector present(R_flat, rs1)
    do c = 1, nflat
    do s = 1, np
      rs1(boff+c, s) = R_flat(boff+c, s)
    end do; end do
    do jac = 1, njacobi
      !$acc parallel loop collapse(3) gang vector present(R_flat, rs1, rs2) private(c0, s)
      do k = 1, n3
      do j = 1, n2
      do i = 1, n1
        c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
        do s = 1, np
          rs2(c0,s) = ( R_flat(c0,s) + beta * ( rs1(c0-off,s) + rs1(c0+off,s) ) ) / ( 1d0 + 2d0*beta )
        end do
      end do; end do; end do
      !$acc parallel loop collapse(3) gang vector present(rs1, rs2) private(c0, s)
      do k = 1, n3
      do j = 1, n2
      do i = 1, n1
        c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
        do s = 1, np
          rs1(c0,s) = rs2(c0,s)
        end do
      end do; end do; end do
    end do
    !$acc parallel loop collapse(3) gang vector present(R_flat, rs1) private(c0, s)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      do s = 1, np
        R_flat(c0,s) = rs1(c0,s)
      end do
    end do; end do; end do
  end subroutine irs_smooth_blk

end module MOSE_Lib_IRS
