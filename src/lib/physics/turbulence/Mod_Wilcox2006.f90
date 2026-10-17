!> @brief Wilcox (2006) k-omega model, as specified by the NASA Turbulence Modeling Resource (TMR):
!> https://tmbwg.github.io/turbmodels/wilcox.html. MOSE_open b29130b Lib_Wilcox2006 on the device.
!> Exact production P = tau_ij du_i/dx_j with the k-equation limiter 20 beta* rho omega k (TMR Wilcox2006-klim), the
!> stress limiter Clim on the eddy viscosity, the cross diffusion sigma_d rho/omega grad k . grad omega where positive,
!> beta = beta_0 f_beta of the vortex-stretching invariant. Without k-energy-coupling, -2/3 rho k delta_ij is left out of
!> the stress of the momentum and energy equations (TMR Wilcox2006-klim-e); with it, k is in the total energy
!> (stress_dev = 2 and the k source out of the energy equation, as SST).
!>
!> Device twin, 2026-10-17: the source runs on the resident var-last w (wilcox_body.inc through Wilcox_Blk, one cell per
!> thread as SST_Blk); the eddy viscosity and the k/omega diffusive flux are `acc routine seq`, dispatched by eddy_dev = 3
!> in the dt kernels, the symmetry and wall BC kernels and the record BC kernels, and inlined in the face kernels
!> (diff_physics.inc). The wall treatment is that of SST in both codes: k_wall = 0, rho omega_wall = omega_wall_c(1) mu_l /
!> d^2 (omega-wall-bc), the ghost 2 omega_wall - omega (Extrapolate_Wall_blk); the floors 1d-10 of the state update.
!>
!> One difference from MOSE_open, on purpose: S^_ij = S_ij - div(u)/2 delta_ij on all nine entries in the invariant
!> X_w of f_beta. MOSE_open b29130b assigns only the diagonal of Sij_hat and leaves the off-diagonal entries undefined.
!> In 2-D they multiply zero entries of W (the vorticity has the z component only), so the two codes agree; in 3-D
!> MOSE_open reads undefined values there.
!>
!> References
!>   Wilcox, D. C., Formulation of the k-omega Turbulence Model Revisited, AIAA J. 46(11), 2008,
!>     pp. 2823-2838, doi:10.2514/1.36541.
!>   Wilcox, D. C., Turbulence Modeling for CFD, 3rd ed., DCW Industries, 2006.
module MOSE_Mod_Wilcox2006
  implicit none
  real(kind=8), parameter, public :: &
  sigma_k = 0.6d0, &
  sigma_w = 0.5d0, &
  beta_star = 0.09d0, &
  beta_0 = 0.0708d0, &
  gamma = 13d0/25d0, &
  Clim = 7d0/8d0

! The device source (wilcox_body.inc) reads the state from the flat field P_flat, as the SST and SA sources: PS(s,i,j,k)
! scalar, PV(i,j,k) the nu:nw velocity slice, PSUM(i,j,k) the species sum, WR the flat row (boff/pi/pj from
! field_block_index_params). The residual is written into R_flat at the same row.
#  define WR(i,j,k) ( boff + ((i)+gc-1) + ((j)+gc-1)*pi + ((k)+gc-1)*pi*pj + 1 )
#  define PS(s,i,j,k) P_flat( WR(i,j,k), s )
#  define PV(i,j,k)   P_flat( WR(i,j,k), nu:nw )
#  define PSUM(i,j,k) sum( P_flat( WR(i,j,k), 1:nsc ) )

contains


  subroutine Wilcox2006_Source_Terms ( domain )
    use MOSE_Types_m
    use MOSE_IO_RANS, only: obj_rans
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Global_m, only: ROTATING_FRAME
    use MOSE_Lib_RotatingFrame, only: obj_rot

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b

    ! One launch (`!$omp single`): the caller runs inside the Mod_Explicit `!$omp parallel` region.
    !$omp single
    do b = 1, domain % nb
      if ( .not. ( is_local_block(b) ) ) cycle
        call Wilcox_Blk ( b,                             &
                          domain % blk(b) % m_f,         &
                          domain % blk(b) % vol,         &
                          domain % blk(b) % dtlocal, obj_rans % point_implicit, &
                          obj_rans % k_energy_coupling,  &
                          ROTATING_FRAME,                &
                          obj_rot % Omega_vec,           &
                          domain % blk(b) % dim )
    end do
    !$omp end single

  end subroutine Wilcox2006_Source_Terms


  !> Device Wilcox (2006) source over the block box: one cell per thread, wilcox_body.inc.
  subroutine Wilcox_Blk ( b, Metric, Volume, dtl, pimp, k_energy_coupling, rotating_frame, omega_vec, n )
    use MOSE_Global_m
    use MOSE_Mod_Fluid, only: Strain_Tensor, Vorticity_Vector, Vorticity_Tensor
    use MOSE_Field_m, only: P_flat, R_flat, field_block_index_params

    implicit none
    integer, intent(in) :: n(3)
    integer, intent(in) :: b
    integer :: boff, pi, pj
    real(kind=8), dimension(3, 3, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Metric
    real(kind=8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Volume
    real(kind=8), dimension(n(1), n(2), n(3)), intent(in) :: dtl   !< local time step (point-implicit destruction)
    logical, intent(in) :: pimp, k_energy_coupling, rotating_frame
    real(kind=8), intent(in) :: omega_vec(3)
    ! Local (per cell, fixed-size)
    integer      :: i, j, k, a, c, kk
    real(kind=8) :: rho, kap, ome, Gtmp(3,3), Gradvel(3,3), Divel, Sij(3,3), diag, Sb2, ome_hat, mi_t, Tij_, Prod(2)
    real(kind=8) :: Grw(2,3), Grad(2,3), dkDotdw, sigma_d, Diff, Sij_hat(3,3), vort(3), Wij(3,3), Num, X_w, f_beta, beta
    real(kind=8) :: Diss(2), Source(2), fk, fw
    real(kind=8) :: sum_im, sum_ip, sum_jm, sum_jp, sum_km, sum_kp
    real(kind=8) :: ov1_, ov2_, ov3_
    ov1_ = omega_vec(1); ov2_ = omega_vec(2); ov3_ = omega_vec(3)

    call field_block_index_params(b, boff, pi, pj)
    !$acc data copyin(Metric, Volume, dtl)
    !$acc parallel loop collapse(3) gang vector present(P_flat, R_flat, Metric, Volume, dtl) &
    !$acc&  private(Gtmp, Gradvel, Sij, Prod, Grw, Grad, Sij_hat, vort, Wij, Diss, Source)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)

#include "wilcox_body.inc"

    end do
    end do
    end do
    !$acc end data

  end subroutine Wilcox_Blk


  !> Eddy viscosity mu_t = rho k / max(omega, Clim sqrt(2 Sbar_ij Sbar_ij / beta*)), Sbar = S - div(u)/3 I. `acc routine
  !> seq`: called on the device from the dt kernels, the symmetry BC kernel and the record BC kernels (eddy_dev = 3);
  !> diff_physics.inc inlines the same lines in the face kernels.
  subroutine Wilcox2006_Eddy_Viscosity ( mi_t, rkw, mi_l, rho, Gradvel, dist )
    !$acc routine seq
    use MOSE_Global_m
    use MOSE_Mod_Fluid

    implicit none
    real(kind=8), intent(in), dimension(nRANS) :: rkw        ! : rkw = [rho*k rho*w]
    real(kind=8), intent(in)                   :: mi_l       ! : Molecular viscosity
    real(kind=8), intent(in)                   :: rho        ! : Density
    real(kind=8), intent(in)                   :: dist       ! : Wall distance (not used)
    real(kind=8), intent(in), dimension(3,3)   :: Gradvel    ! : Velocity gradient
    real(kind=8), intent(out)                  :: mi_t       ! : Eddy viscosity
    ! Local
    real(kind=8) :: w, Divel, Sij(3,3), diag, w_hat

    w = rkw(2) / rho
    Divel = Gradvel(1,1) + Gradvel(2,2) + Gradvel(3,3)
    Sij = Strain_Tensor ( Gradvel )
    diag = 1d0/3d0 * Divel
    Sij(1,1) = Sij(1,1) - diag
    Sij(2,2) = Sij(2,2) - diag
    Sij(3,3) = Sij(3,3) - diag
    w_hat = Max ( w, Clim * Sqrt( 2d0 * Sum( Sij**2 ) / beta_star ) )
    mi_t = rkw(1) / w_hat

  end subroutine Wilcox2006_Eddy_Viscosity


  !> Diffusive flux of rho k and rho omega: area (mu_l + sigma rho k / omega) grad . n, with the unlimited omega (not the
  !> eddy viscosity of the stress limiter). `acc routine seq`.
  subroutine Wilcox2006_RANS_Diffusive_Flux ( Flux, rkw, mi_l, vGrad, Gradkw, rho, Area, Normal, Dist )
    !$acc routine seq
    use MOSE_Global_m

    implicit none
    real(kind=8), intent(in), dimension(nRANS)   :: rkw                    ! : [ rho*k rho*w ]
    real(kind=8), intent(in)                     :: mi_l                   ! : laminar viscosity
    real(kind=8), intent(in), dimension(3,3)     :: vGrad                  ! : Velocity gradient (not used)
    real(kind=8), intent(in), dimension(nRANS,3) :: Gradkw                 ! : Gradient of [ k w ]
    real(kind=8), intent(in)                     :: rho                    ! : Density
    real(kind=8), intent(in)                     :: Area, Normal(3), Dist  ! : Metrics
    real(kind=8), intent(out), dimension(nRANS)  :: Flux                   ! : Diffusive flux
    ! Local
    real(kind=8) :: rhok, w

    rhok = rkw(1)
    w = rkw(2)/rho
    ! explicit dot products: a strided row of the (nRANS,3) dummy inside an acc routine is misread by nvfortran
    Flux(1) = ( mi_l + sigma_k * rhok / w ) * ( Gradkw(1,1)*Normal(1) + Gradkw(1,2)*Normal(2) + Gradkw(1,3)*Normal(3) )
    Flux(2) = ( mi_l + sigma_w * rhok / w ) * ( Gradkw(2,1)*Normal(1) + Gradkw(2,2)*Normal(2) + Gradkw(2,3)*Normal(3) )
    Flux(1) = Flux(1) * Area
    Flux(2) = Flux(2) * Area

  end subroutine Wilcox2006_RANS_Diffusive_Flux


  subroutine Wilcox2006_Set_Wall_Values ( mi_l, rkw, dist )
    use MOSE_Global_m

    implicit none
    real(kind=8), intent(out), dimension(nRANS)  :: rkw              ! : [ rho*k rho*w ] at wall
    real(kind=8), intent(in)                     :: mi_l             ! : laminar viscosity at wall
    real(kind=8), intent(in)                     :: dist             ! : cell center distance from wall

    rkw(1) = 0d0 ! Solid surface condition on k
    rkw(2) = omega_wall_c(1) * mi_l / ( dist**2 ) ! smooth surface: Menter practical (800) or asymptotic (80), omega-wall-bc

  end subroutine Wilcox2006_Set_Wall_Values


  subroutine Wilcox2006_Extrapolate_Wall ( prim, wall_prim, rho, wall_rho, ghost_prim )
    use MOSE_Global_m

    implicit none
    real(kind=8), intent(in),  dimension(nRANS) :: prim, wall_prim  ! : boundary cell and wall variables
    real(kind=8), intent(in) :: rho, wall_rho ! : boundary cell density and wall density
    real(kind=8), intent(out), dimension(nRANS) :: ghost_prim ! : ghost cell RANS variables

    ghost_prim(1) = - prim(1) ! enforcing k_wall == 0
    ghost_prim(2) = rho * 2d0 * wall_prim(2)/wall_rho - prim(2) ! extrapolating from omega_wall

  end subroutine Wilcox2006_Extrapolate_Wall


  subroutine Wilcox2006_Enforce_Realizability ( rhokrhow )

    implicit none
    real(kind=8), intent(inout),  dimension(:) :: rhokrhow
    if ( rhokrhow(1) < 1d-10 ) rhokrhow(1) = 1d-10
    if ( rhokrhow(2) < 1d-10 ) rhokrhow(2) = 1d-10

  end subroutine Wilcox2006_Enforce_Realizability


end module MOSE_Mod_Wilcox2006
