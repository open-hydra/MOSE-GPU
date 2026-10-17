!> @brief The rotation/curvature correction of SA-RC (Spalart & Shur 1997; Shur et al. 2000) and SST-RC (Smirnov & Menter
!> 2009) on the device, as MOSE_open b29130b Lib_SpalartShur: r* = S/Omega and the numerator of r~, 2 W_ik S_jk DS_ij/Dt,
!> shared by the two models (fr1, with its constants and the denominator of r~, is computed in each source).
!> DS_ij/Dt is the steady convective part u.grad(S_ij); the frame-rotation term (eps_imn S_jn + eps_jmn S_in) Omega_m of
!> the published r~ is not included (as MOSE_open), so in a rotating frame the correction sees the frame only through the
!> absolute vorticity in W_ik and Omega.
!>
!> Device layout, on the flat rows of P_flat (WR below): gv_f(row, a + 3*(c-1)) = dU_a/dx_c, the cell-centred velocity
!> gradient of Compute_Velocity_Gradient on the interior cells and on the first ghost layer of the six faces, which takes
!> the value of the adjacent interior cell as in MOSE_open (no exchange between blocks: the strain gradient of a boundary
!> cell is one-sided); rc1_f = r* and rc2_f = 2 W_ik S_jk DS_ij/Dt on the interior rows. Allocated at the first call:
!> only the -RC models use them.
module MOSE_Lib_SpalartShur
  use iso_fortran_env, only: R8 => real64
  implicit none
  private
  public :: RC_Prepare, gv_f, rc1_f, rc2_f

  real(R8), allocatable, save :: gv_f(:,:), rc1_f(:), rc2_f(:)

#  define WR(i,j,k) ( boff + ((i)+gc-1) + ((j)+gc-1)*pi + ((k)+gc-1)*pi*pj + 1 )

contains


  !> MOSE_open Compute_Velocity_Gradient and Compute_RC_Terms on every local block, before the -RC sources.
  subroutine RC_Prepare ( domain )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Global_m, only: ROTATING_FRAME
    use MOSE_Lib_RotatingFrame, only: obj_rot
    use MOSE_Field_m, only: ncell_total

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b

    if ( .not. allocated(gv_f) ) then
      allocate ( gv_f(ncell_total, 9), rc1_f(ncell_total), rc2_f(ncell_total) )
      gv_f = 0d0; rc1_f = 0d0; rc2_f = 0d0
      !$acc enter data copyin(gv_f, rc1_f, rc2_f)
    end if
    do b = 1, domain % nb
      if ( .not. is_local_block(b) ) cycle
      call Velocity_Gradient_Blk ( b, domain % blk(b) % m_f, domain % blk(b) % dim )
    end do
    do b = 1, domain % nb
      if ( .not. is_local_block(b) ) cycle
      call RC_Terms_Blk ( b, domain % blk(b) % m_f, ROTATING_FRAME, obj_rot % Omega_vec, domain % blk(b) % dim )
    end do

  end subroutine RC_Prepare


  !> Compute_Velocity_Gradient on one block: dU/dx = matmul(central differences + 1d-40, M) on the interior cells, then the
  !> ghost layer of each face from the adjacent interior cell (the other two indices on the interior range, as MOSE_open).
  subroutine Velocity_Gradient_Blk ( b, Metric, n )
    use MOSE_Global_m, only: gc, nu
    use MOSE_Field_m, only: P_flat, field_block_index_params

    implicit none
    integer, intent(in) :: b, n(3)
    real(R8), dimension(3, 3, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Metric
    ! Local
    integer  :: boff, pi, pj, i, j, k, a, c, s, n1, n2, n3
    real(R8) :: vd(3,3)

    call field_block_index_params(b, boff, pi, pj)
    n1 = n(1); n2 = n(2); n3 = n(3)
    !$acc data copyin(Metric)
    !$acc parallel loop collapse(3) gang vector present(P_flat, gv_f, Metric) private(vd)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      do a = 1, 3
        vd(a,1) = ( P_flat(WR(i+1,j,k), nu-1+a) - P_flat(WR(i-1,j,k), nu-1+a) )/2d0
        vd(a,2) = ( P_flat(WR(i,j+1,k), nu-1+a) - P_flat(WR(i,j-1,k), nu-1+a) )/2d0
        vd(a,3) = ( P_flat(WR(i,j,k+1), nu-1+a) - P_flat(WR(i,j,k-1), nu-1+a) )/2d0
      end do
      do c = 1, 3
      do a = 1, 3
        gv_f(WR(i,j,k), a+3*(c-1)) = ( vd(a,1) + 1d-40 )*Metric(1,c,i,j,k) + ( vd(a,2) + 1d-40 )*Metric(2,c,i,j,k) &
                                   + ( vd(a,3) + 1d-40 )*Metric(3,c,i,j,k)
      end do
      end do
    end do
    end do
    end do
    ! i faces
    !$acc parallel loop collapse(2) gang vector present(gv_f)
    do k = 1, n3
    do j = 1, n2
      do s = 1, 9
        gv_f(WR(0,j,k), s)    = gv_f(WR(1,j,k), s)
        gv_f(WR(n1+1,j,k), s) = gv_f(WR(n1,j,k), s)
      end do
    end do
    end do
    ! j faces
    !$acc parallel loop collapse(2) gang vector present(gv_f)
    do k = 1, n3
    do i = 1, n1
      do s = 1, 9
        gv_f(WR(i,0,k), s)    = gv_f(WR(i,1,k), s)
        gv_f(WR(i,n2+1,k), s) = gv_f(WR(i,n2,k), s)
      end do
    end do
    end do
    ! k faces
    !$acc parallel loop collapse(2) gang vector present(gv_f)
    do j = 1, n2
    do i = 1, n1
      do s = 1, 9
        gv_f(WR(i,j,0), s)    = gv_f(WR(i,j,1), s)
        gv_f(WR(i,j,n3+1), s) = gv_f(WR(i,j,n3), s)
      end do
    end do
    end do
    !$acc end data

  end subroutine Velocity_Gradient_Blk


  !> Compute_RC_Terms on one block: the strain gradient by central differences of the strain tensors of the neighbours,
  !> to the physical frame component by component; r* = S/Omega (Omega = |absolute vorticity|, at least 1d-30) and
  !> 2 W_ik S_jk (u . grad S_ij), the numerator of r~ (its denominator differs between SA and SST).
  subroutine RC_Terms_Blk ( b, Metric, rotating_frame, omega_vec, n )
    use MOSE_Global_m, only: gc, nu, nv, nw
    use MOSE_Field_m, only: P_flat, field_block_index_params
    use MOSE_Mod_Fluid, only: Strain_Tensor, Vorticity_Vector, Vorticity_Tensor

    implicit none
    integer, intent(in) :: b, n(3)
    real(R8), dimension(3, 3, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Metric
    logical, intent(in) :: rotating_frame
    real(R8), intent(in) :: omega_vec(3)
    ! Local
    integer  :: boff, pi, pj, i, j, k, a, c, d, ii, jj, kk, rp, rm
    real(R8) :: Gp(3,3), Gm(3,3), Sp(3,3), Sm(3,3), sg(3,3,3), sgm(3,3,3), Sij(3,3), Wij(3,3), st(3,3), omega(3)
    real(R8) :: S, O, rtilde, ov1_, ov2_, ov3_

    call field_block_index_params(b, boff, pi, pj)
    ov1_ = omega_vec(1); ov2_ = omega_vec(2); ov3_ = omega_vec(3)
    !$acc data copyin(Metric)
    !$acc parallel loop collapse(3) gang vector present(P_flat, gv_f, rc1_f, rc2_f, Metric) &
    !$acc&  private(Gp, Gm, Sp, Sm, sg, sgm, Sij, Wij, st, omega)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)

      ! Finite difference strain rate tensor
      do d = 1, 3
        if ( d == 1 ) then
          rp = WR(i+1,j,k); rm = WR(i-1,j,k)
        else if ( d == 2 ) then
          rp = WR(i,j+1,k); rm = WR(i,j-1,k)
        else
          rp = WR(i,j,k+1); rm = WR(i,j,k-1)
        end if
        do c = 1, 3
        do a = 1, 3
          Gp(a,c) = gv_f(rp, a+3*(c-1))
          Gm(a,c) = gv_f(rm, a+3*(c-1))
        end do
        end do
        Sp = Strain_Tensor ( Gp )
        Sm = Strain_Tensor ( Gm )
        do jj = 1, 3
        do ii = 1, 3
          sg(ii,jj,d) = ( Sp(ii,jj) - Sm(ii,jj) )*0.5d0
        end do
        end do
      end do

      ! Transformation component by component: sgm(ii,jj,:) = matmul(sg(ii,jj,:), M)
      do jj = 1, 3
      do ii = 1, 3
        do c = 1, 3
          sgm(ii,jj,c) = sg(ii,jj,1)*Metric(1,c,i,j,k) + sg(ii,jj,2)*Metric(2,c,i,j,k) + sg(ii,jj,3)*Metric(3,c,i,j,k)
        end do
      end do
      end do

      do c = 1, 3
      do a = 1, 3
        Gp(a,c) = gv_f(WR(i,j,k), a+3*(c-1))
      end do
      end do
      Sij = Strain_Tensor ( Gp )
      omega = Vorticity_Vector ( Gp )
      ! the absolute vorticity in a rotating frame: omega_abs = omega_rel + 2*Omega_frame
      if ( rotating_frame ) then
        omega(1) = omega(1) + 2.0d0 * ov1_
        omega(2) = omega(2) + 2.0d0 * ov2_
        omega(3) = omega(3) + 2.0d0 * ov3_
      end if
      Wij = Vorticity_Tensor ( omega )
      S = sqrt ( 2d0*sum (Sij**2) )
      O = max ( sqrt ( sum(omega**2) ), 1d-30 )
      rc1_f(WR(i,j,k)) = S/O

      ! Material derivative of the strain rate tensor, the time term neglected
      do jj = 1, 3
      do ii = 1, 3
        st(ii,jj) = P_flat(WR(i,j,k), nu)*sgm(ii,jj,1) + P_flat(WR(i,j,k), nv)*sgm(ii,jj,2) &
                  + P_flat(WR(i,j,k), nw)*sgm(ii,jj,3)
      end do
      end do
      rtilde = 0d0
      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        rtilde = rtilde + Wij(ii,kk)*Sij(jj,kk)*st(ii,jj)
      end do
      end do
      end do
      rc2_f(WR(i,j,k)) = 2d0*rtilde   ! the denominator of r~ is in each source

    end do
    end do
    end do
    !$acc end data

  end subroutine RC_Terms_Blk

end module MOSE_Lib_SpalartShur
