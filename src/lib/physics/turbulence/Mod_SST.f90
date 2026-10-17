!> @brief Module for Kappa-Omega SST (2003) model. see https://turbmodels.larc.nasa.gov/sst.html.
module MOSE_Mod_SST
  implicit none
  real(kind=8), parameter, public :: &
  sigma_k1 = 0.85d0, &
  sigma_k2 = 1.0d0, &
  sigma_w1 = 0.5d0, &
  sigma_w2 = 0.856d0, &
  beta_1 = 0.075d0, &
  beta_2 = 0.0828d0, &
  beta_star = 0.09d0, &
  karman = 0.41d0, &
  alpha_1 = 0.31d0, &
  gamma_1 = 5d0/9d0, &
  gamma_2 = 0.44d0, &
  cr1 = 1.d0, &
  cr2 = 2.d0, &
  cr3 = 1.d0

! The device SST source (SST_Blk) reads state from the
! var-LAST flat SoA w on the GPU device path (_OPENACC), else var-first Prim(=P). It reads NEIGHBOUR
! cells (i+-1,j+-1,k+-1) for the velocity + k/omega gradients -> at the inter-rank seam
! those neighbours are GHOSTS, so this reader MUST be on w once the n>1 halo fills w ghosts
! (not P). Accessors: PS(s,i,j,k) scalar, PV(i,j,k) the nu:nw velocity slice, PSUM(i,j,k)
! the species-sum; WR = flat row (boff/pi/pj from field_block_index_params). Res stays
! var-first (NOT switched). Flag OFF -> PS/PV/PSUM == Prim(...), byte-identical.
#  define SSTARG  b
#  define SSTPASS b
#  define WR(i,j,k) ( boff + ((i)+gc-1) + ((j)+gc-1)*pi + ((k)+gc-1)*pi*pj + 1 )
#  define PS(s,i,j,k) P_flat( WR(i,j,k), s )
#  define PV(i,j,k)   P_flat( WR(i,j,k), nu:nw )
#  define PSUM(i,j,k) sum( P_flat( WR(i,j,k), 1:nsc ) )

contains


  subroutine SST_Source_Terms ( domain )
    use MOSE_Types_m
    use MOSE_IO_RANS, only: obj_rans
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Global_m, only: ROTATING_FRAME
    use MOSE_Lib_RotatingFrame, only: obj_rot
    use MOSE_Lib_SpalartShur, only: RC_Prepare

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b

    ! One launch (`!$omp single`): SST_Source_Terms runs inside the Mod_Explicit `!$omp parallel` region.
    !$omp single
    if ( obj_rans % SpalartShur ) then   ! SST-RC: the velocity gradient and the correction terms first (MOSE_open order)
      call RC_Prepare ( domain )
      do b = 1, domain % nb
        if ( .not. ( is_local_block(b) ) ) cycle
        call SST_RC_Blk ( SSTPASS, domain % blk(b) % r, domain % blk(b) % m_f, domain % blk(b) % vol, &
                          domain % blk(b) % yn, domain % blk(b) % rho_c, domain % blk(b) % mil_c, &
                          domain % blk(b) % dtlocal, obj_rans % point_implicit, obj_rans % k_energy_coupling, &
                          ROTATING_FRAME, obj_rot % Omega_vec, domain % blk(b) % dim )
      end do
    else
    do b = 1, domain % nb
      if ( .not. ( is_local_block(b) ) ) cycle
        call SST_Blk ( SSTPASS,                  &
                           domain % blk(b) % r,      &
                           domain % blk(b) % m_f,    &
                           domain % blk(b) % vol,    &
                           domain % blk(b) % yn,     &
                           domain % blk(b) % rho_c,  &
                           domain % blk(b) % mil_c,  &
                           domain % blk(b) % dtlocal, obj_rans % point_implicit, &
                           obj_rans % k_energy_coupling, &
                           ROTATING_FRAME,           &
                           obj_rot % Omega_vec,      &
                           domain % blk(b) % dim )
    end do
    end if
    !$omp end single

  end subroutine SST_Source_Terms



  !> Device twin of SST_Blk (the RANS/SST source term); SST-RC runs its twin SST_RC_Blk below (the same body with
  !> MOSE_RC_KERNEL: the stored gradient of Lib_SpalartShur and fr1).
  !> Mirrors the host loop math; only the per-cell arithmetic is reordered (sums
  !> written explicitly to avoid device array temporaries), so it is validated to
  !> the device-FP oracle (~1e-8) vs CPU, not bit-identical.
  !>  - the derived-type metric M(i,j,k)%c is replaced by the flat, already-
  !>    device-resident Metric(:,:,i,j,k) (== blk%m_f).
  !>  - Grad rows are gathered into contiguous gradk_c/gradw_c before compute_F
  !>    (the nvfortran strided-section gotcha, same as SST_RANS_Diffusive_Flux).
  subroutine SST_Blk ( SSTARG, Res, Metric, Volume, WDist, rho_c, mil_c, dtl, pimp, &
                           k_energy_coupling, rotating_frame, omega_vec, n)
    use MOSE_Global_m
    use MOSE_Mod_Fluid, only: Strain_Tensor, Vorticity_Vector, Vorticity_Tensor
    use MOSE_Field_m, only: P_flat, R_flat, ncell_total, field_block_index_params, rho_c_f, mil_c_f

    implicit none
    integer, intent(in) :: n(3)
    integer, intent(in) :: b
    integer :: boff, pi, pj
    real(kind=8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: Res
    real(kind=8), dimension(3, 3, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Metric
    real(kind=8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Volume, WDist
    real(kind=8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: rho_c, mil_c
    real(kind=8), dimension(n(1), n(2), n(3)), intent(in) :: dtl   !< local time step (point-implicit destruction)
    logical, intent(in) :: pimp
    real(kind=8) :: fk_, fw_
    logical, intent(in) :: k_energy_coupling, rotating_frame
    real(kind=8), intent(in) :: omega_vec(3)
    ! Local (per cell, fixed-size)
    integer      :: i, j, k, a, c
    real(kind=8) :: rho, mil, kap, ome, dist, Gtmp(3,3), Gradvel(3,3), Divel
    real(kind=8) :: Sij(3,3), S, S2, Wij(3,3), O, O2, mi_t, F(2)
    real(kind=8) :: vort(3), Grw(2,3), Grad(2,3), gradk_c(3), gradw_c(3)
    real(kind=8) :: Prod(2), dkDotdw, Diff, beta, gamma, Diss(2), Source(2)
    real(kind=8) :: sum_im, sum_ip, sum_jm, sum_jp, sum_km, sum_kp
    ! Dims as scalars -- a dummy array referenced inside the acc region
    ! is implicitly copy()d per launch (upload+download+sync); scalars ride the
    ! launch parameters for free.
    integer :: nb1_, nb2_, nb3_
    real(kind=8) :: ov1_, ov2_, ov3_
    nb1_ = n(1); nb2_ = n(2); nb3_ = n(3)
    ov1_ = omega_vec(1); ov2_ = omega_vec(2); ov3_ = omega_vec(3)

    call field_block_index_params(b, boff, pi, pj)
    !$acc data copy(Res) copyin(Metric, Volume, WDist, rho_c, mil_c, dtl)
    !$acc parallel loop collapse(3) gang vector &
    !$acc&  present(P_flat, R_flat, Res, Metric, Volume, WDist, rho_c, mil_c, dtl, rho_c_f, mil_c_f) &
    !$acc&  private(Gtmp, Gradvel, Sij, Wij, vort, Grw, Grad, gradk_c, gradw_c, &
    !$acc&          Prod, Diss, Source, F)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)

#include "sst_body.inc"

    end do
    end do
    end do
    !$acc end data

  end subroutine SST_Blk


  !> The -RC twin (MOSE_RC_KERNEL): the stored gradient and the correction factor fr1 on the production.
  subroutine SST_RC_Blk ( SSTARG, Res, Metric, Volume, WDist, rho_c, mil_c, dtl, pimp, &
                           k_energy_coupling, rotating_frame, omega_vec, n)
    use MOSE_Global_m
    use MOSE_Mod_Fluid, only: Strain_Tensor, Vorticity_Vector, Vorticity_Tensor
    use MOSE_Lib_SpalartShur, only: gv_f, rc1_f, rc2_f
    use MOSE_Field_m, only: P_flat, R_flat, ncell_total, field_block_index_params, rho_c_f, mil_c_f

    implicit none
    integer, intent(in) :: n(3)
    integer, intent(in) :: b
    integer :: boff, pi, pj
    real(kind=8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: Res
    real(kind=8), dimension(3, 3, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Metric
    real(kind=8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Volume, WDist
    real(kind=8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: rho_c, mil_c
    real(kind=8), dimension(n(1), n(2), n(3)), intent(in) :: dtl   !< local time step (point-implicit destruction)
    logical, intent(in) :: pimp
    real(kind=8) :: fk_, fw_
    logical, intent(in) :: k_energy_coupling, rotating_frame
    real(kind=8), intent(in) :: omega_vec(3)
    ! Local (per cell, fixed-size)
    real(kind=8) :: rstar, Drc, rtilde, frot, fr1   ! SST-RC (Smirnov & Menter 2009)
    integer      :: i, j, k, a, c
    real(kind=8) :: rho, mil, kap, ome, dist, Gtmp(3,3), Gradvel(3,3), Divel
    real(kind=8) :: Sij(3,3), S, S2, Wij(3,3), O, O2, mi_t, F(2)
    real(kind=8) :: vort(3), Grw(2,3), Grad(2,3), gradk_c(3), gradw_c(3)
    real(kind=8) :: Prod(2), dkDotdw, Diff, beta, gamma, Diss(2), Source(2)
    real(kind=8) :: sum_im, sum_ip, sum_jm, sum_jp, sum_km, sum_kp
    ! Dims as scalars -- a dummy array referenced inside the acc region
    ! is implicitly copy()d per launch (upload+download+sync); scalars ride the
    ! launch parameters for free.
    integer :: nb1_, nb2_, nb3_
    real(kind=8) :: ov1_, ov2_, ov3_
    nb1_ = n(1); nb2_ = n(2); nb3_ = n(3)
    ov1_ = omega_vec(1); ov2_ = omega_vec(2); ov3_ = omega_vec(3)

    call field_block_index_params(b, boff, pi, pj)
    !$acc data copy(Res) copyin(Metric, Volume, WDist, rho_c, mil_c, dtl)
    !$acc parallel loop collapse(3) gang vector &
    !$acc&  present(P_flat, R_flat, gv_f, rc1_f, rc2_f, Res, Metric, Volume, WDist, rho_c, mil_c, dtl, rho_c_f, mil_c_f) &
    !$acc&  private(Gtmp, Gradvel, Sij, Wij, vort, Grw, Grad, gradk_c, gradw_c, &
    !$acc&          Prod, Diss, Source, F)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)

#define MOSE_RC_KERNEL
#include "sst_body.inc"
#undef MOSE_RC_KERNEL

    end do
    end do
    end do
    !$acc end data

  end subroutine SST_RC_Blk




  subroutine SST_Eddy_Viscosity ( mi_t, rkw, mi_l, rho, Gradvel, dist )
    !$acc routine seq
    use MOSE_Global_m
    use MOSE_Mod_Fluid

    implicit none
    real(kind=8), intent(in), dimension(nRANS) :: rkw        ! : rkw = [rho*k rho*w]
    real(kind=8), intent(in)                   :: mi_l       ! : Molecular viscosity
    real(kind=8), intent(in)                   :: rho        ! : Density
    real(kind=8), intent(in)                   :: dist       ! : Wall distance
    real(kind=8), intent(in), dimension(3,3)   :: Gradvel    ! : Velocity gradient
    real(kind=8), intent(out)                  :: mi_t       ! : Eddy viscosity
    ! Local
    real(kind=8) :: kap, ome, Sij(3,3), S, F2, vort(3), Wij(3,3), O

    kap = rkw(1) / rho
    ome = rkw(2) / rho
    F2 = compute_F2 ( rho=rho, w=ome, k=kap, mil=mi_l, d=dist )
    Sij = Strain_Tensor ( Gradvel )
    S = Sqrt ( 2d0 * sum ( Sij**2 ) )
    vort = vorticity_vector ( gradvel )
    Wij = Vorticity_Tensor ( vort )
    O = Sqrt ( 2d0 * sum ( Wij**2 ) )
    mi_t = rho*kap*alpha_1 / Max ( alpha_1*ome, S*F2 )

  end subroutine SST_Eddy_Viscosity


  subroutine SST_RANS_Diffusive_Flux ( Flux, rkw, mi_l, Gradvel, Gradkw, rho, Area, Normal, Dist )
    !$acc routine seq
    use MOSE_Global_m

    implicit none
    real(kind=8), intent(in), dimension(nRANS)   :: rkw                    ! : [ rho*k rho*w ]
    real(kind=8), intent(in)                     :: mi_l                   ! : laminar viscosity
    real(kind=8), intent(in), dimension(3,3)     :: Gradvel                ! : Velocity gradient
    real(kind=8), intent(in), dimension(nRANS,3) :: Gradkw                 ! : Gradient of [ k w ]
    real(kind=8), intent(in)                     :: rho                    ! : Density
    real(kind=8), intent(in)                     :: Area, Normal(3), Dist  ! : Metrics
    real(kind=8), intent(out), dimension(nRANS)  :: Flux                   ! : Diffusive flux
    ! Local
    real(kind=8) :: kap, ome, F1, sigma_k, sigma_w, mi_t
    ! Contiguous gathers of the strided Gradkw rows — passing Gradkw(1,:) (a row
    ! of the (nRANS,3) dummy, stride nRANS) straight into compute_F1's Gradk(3)
    ! is misread by nvfortran on the device. (CPU path unchanged.)
    real(kind=8) :: gradk_c(3), gradw_c(3)

    kap = rkw(1)/rho
    ome = rkw(2)/rho
    gradk_c(1)=Gradkw(1,1); gradk_c(2)=Gradkw(1,2); gradk_c(3)=Gradkw(1,3)
    gradw_c(1)=Gradkw(2,1); gradw_c(2)=Gradkw(2,2); gradw_c(3)=Gradkw(2,3)
    F1 = compute_F1 (rho=rho, w=ome, k=kap, mil=mi_l, &
                     gradk=gradk_c, gradw=gradw_c, d=dist)
    sigma_k = sigma_k1*F1 + sigma_k2*(1d0-F1)
    sigma_w = sigma_w1*F1 + sigma_w2*(1d0-F1)
    call SST_Eddy_Viscosity( mi_t=mi_t, rkw=rkw, mi_l=mi_l, rho=rho, &
                             Gradvel=Gradvel, dist=dist )
    Flux(1) = ( mi_l + sigma_k*mi_t ) * Dot_Product ( Gradkw(1,:), Normal )
    Flux(2) = ( mi_l + sigma_w*mi_t ) * Dot_Product ( Gradkw(2,:), Normal )
    Flux = Flux * Area

  end subroutine SST_RANS_Diffusive_Flux


  subroutine SST_Set_Wall_Values ( mi_l, rkw, dist )
    use MOSE_Global_m

    implicit none
    real(kind=8), intent(out), dimension(nRANS)  :: rkw              ! : [ rho*k rho*w ] at wall
    real(kind=8), intent(in)                     :: mi_l             ! : laminar viscosity at wall
    real(kind=8), intent(in)                     :: dist             ! : cell center distance from wall
    
    rkw(1) = 0d0 ! Solid surface condition on k
    rkw(2) = omega_wall_c(1) * mi_l / ( dist**2 ) ! smooth surface: Menter practical (800) or asymptotic (80), omega-wall-bc

  end subroutine SST_Set_Wall_Values


  subroutine SST_Blowing_Correction ( rho, mil, tau, rkw, mdot, dist )
    use MOSE_Global_m

    implicit none
    real(kind=8), intent(out), dimension(nRANS)  :: rkw
    real(kind=8), intent(in)                     :: rho, mil, mdot, dist
    real(kind=8), dimension(3), intent(in)       :: tau
    ! Local
    real(kind=8) :: tauw, uT, v_w, v_plus, SB
    
    rkw(1) = 0d0 ! Solid surface condition on k
    tauw = sqrt( sum(tau**2) )
    uT = sqrt ( tauw / rho )
    v_w = abs(mdot)/rho
    v_plus = v_w / uT
    SB = 25d0 / ( v_plus*( 1d0 + 5d0 * v_plus ) )
    rkw(2) = (rho*uT)**2 * SB / mil ! Wilcox BC for mass injection

  end subroutine SST_Blowing_Correction


  subroutine SST_Blowing_noCorrection ( rho, mil, tau, rkw, mdot, dist )
    use MOSE_Global_m

    implicit none
    real(kind=8), intent(out), dimension(nRANS)  :: rkw
    real(kind=8), intent(in)                     :: rho, mil, mdot, dist
    real(kind=8), dimension(3), intent(in)       :: tau
    
    rkw(1) = 0d0
    rkw(2) = omega_wall_c(1) * mil / ( dist**2 )

  end subroutine SST_Blowing_noCorrection


  subroutine SST_Extrapolate_Wall ( prim, wall_prim, rho, wall_rho, ghost_prim )
    use MOSE_Global_m
    implicit none
    real(kind=8), intent(in),  dimension(nRANS) :: prim, wall_prim  ! : boundary cell and wall variables
    real(kind=8), intent(in) :: rho, wall_rho ! : boundary cell density and wall density
    real(kind=8), intent(out), dimension(nRANS) :: ghost_prim ! : ghost cell RANS variables

    ghost_prim(1) = - prim(1) ! enforcing k_wall == 0
    ghost_prim(2) = rho * 2d0 * wall_prim(2)/wall_rho - prim(2) ! extrapolating from omega_wall

  end subroutine SST_Extrapolate_Wall


  subroutine SST_Enforce_Realizability ( rhokrhow )

    implicit none
    real(kind=8), intent(inout),  dimension(:) :: rhokrhow
    if ( rhokrhow(1) < 1d-10 ) rhokrhow(1) = 1d-10
    if ( rhokrhow(2) < 1d-10 ) rhokrhow(2) = 1d-10

  end subroutine SST_Enforce_Realizability


  function compute_F ( rho, w, k, mil, Gradw, Gradk, d ) result ( F )
    !$acc routine seq
    implicit none
    real(kind=8), intent(in) :: rho, w, k, mil, Gradw(3), Gradk(3), d
    real(kind=8) :: CD, arg1, arg2, F(2)

    CD = Max ( 2d0*rho*sigma_w2/w*Dot_Product( Gradw, Gradk ), 1d-10 )
    arg1 = Min ( Max ( Sqrt(k)/(beta_star*w*d), &
                       5d2*mil/(rho*w*d**2) ), &
                 4d0*sigma_w2*rho*k/(CD*d**2) )
    arg2 = Max ( 2d0*Sqrt(k)/(beta_star*w*d), 5d2*mil/(rho*w*d**2) )
    F(1) = Tanh ( arg1**4 )
    F(2) = Tanh ( arg2**2 )
  end function compute_F

  function compute_F1 ( rho, w, k, mil, Gradw, Gradk, d ) result ( F1 )
    !$acc routine seq
    implicit none
    real(kind=8), intent(in) :: rho, w, k, mil, Gradw(3), Gradk(3), d
    real(kind=8) :: CD, arg1, F1

    CD = Max ( 2d0*rho*sigma_w2/w*Dot_Product( Gradw, Gradk ), 1d-10 )
    arg1 = Min ( Max ( Sqrt(k)/(beta_star*w*d), &
                       5d2*mil/(rho*w*d**2) ), &
                 4d0*sigma_w2*rho*k/(CD*d**2) )
    F1 = Tanh ( arg1**4 )
  end function compute_F1

  function compute_F2 ( rho, w, k, mil, d ) result ( F2 )
    !$acc routine seq
    implicit none
    real(kind=8), intent(in) :: rho, w, k, mil, d
    real(kind=8) :: arg2, F2

    arg2 = Max ( 2d0*Sqrt(k)/(beta_star*w*d), 5d2*mil/(rho*w*d**2) )
    F2 = Tanh ( arg2**2 )

  end function compute_F2
  
  
end module MOSE_Mod_SST
