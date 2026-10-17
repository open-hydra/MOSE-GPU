!> @brief Module for Spalart-Allmaras constants and subroutines, see https://turbmodels.larc.nasa.gov/spalart.html. Note: no ft2 term.
!>
!> Device (OpenACC) twin, 2026-10-15: the source runs on the resident var-last w like the SST one (sa_body.inc through
!> SA_Source_Blk, the block sweep), the eddy viscosity and
!> the nu-tilde diffusive flux are `acc routine seq` and are dispatched by eddy_dev = 2 in the diffusive kernels
!> (diff_physics.inc), the dt kernels and the wall BC kernels. SA, SA-R, SA-RC (SA_RC_Source_Blk, the twin with the
!> stored gradient and fr1 of Lib_SpalartShur) and SAcomp (SA_Comp_Source_Blk, the twin with the Paciorri-Sabetta
!> correction), 2026-10-17.
module MOSE_Mod_Spalart
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  real(R8), parameter, public :: &
    cb1 = 0.1355d0, &
    cb2 = 0.622d0, &
    cw2 = 0.3d0, &
    cw3 = 2.0d0, &
    cv1 = 7.1d0, &
    cK = 0.41d0, &
    ct1 = 1.d0, &
    ct2 = 2.d0, &
    ct3 = 1.2d0, &
    ct4 = 0.5d0, &
    sigma = 2d0/3d0, &
    c2 = 0.7d0, &
    c3 = 0.9d0, &
    cr1 = 1.d0, &
    cr2 = 12.d0, &
    cr3 = 1.d0, &
    Crot = 2.d0, &
    cw1 = cb1/cK**2 + (1d0+cb2)/sigma, &
    cK_sq = cK**2, &   ! for sa_body.inc: the cell-list kernels have a dummy ck(:) that hides cK (case-insensitive names)
    C_cr1 = 0.3d0, &
    tolerance = 1.d-2

! The device SA source (sa_body.inc) reads the state from the flat field P_flat on the GPU device path (_OPENACC), else
! the var-first Prim (= P). Same accessors as the SST source: PS(s,i,j,k) scalar, PV(i,j,k) the nu:nw velocity slice,
! PSUM(i,j,k) the species sum, WR the flat row (boff/pi/pj from field_block_index_params). Res stays var-first.
#  define SAARG  b
#  define SAPASS b
#  define WR(i,j,k) ( boff + ((i)+gc-1) + ((j)+gc-1)*pi + ((k)+gc-1)*pi*pj + 1 )
#  define PS(s,i,j,k) P_flat( WR(i,j,k), s )
#  define PV(i,j,k)   P_flat( WR(i,j,k), nu:nw )
#  define PSUM(i,j,k) sum( P_flat( WR(i,j,k), 1:nsc ) )

contains

  subroutine Spalart_Source_Terms ( domain )
    use MOSE_Types_m
    use MOSE_Mod_Fluid, only : vorticity_vector, vorticity_tensor
    use MOSE_IO_RANS, only: obj_rans
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Global_m, only: ROTATING_FRAME
    use MOSE_Lib_RotatingFrame, only: obj_rot
    use MOSE_Lib_SpalartShur, only: RC_Prepare

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b

    ! One launch (`!$omp single`): the caller runs inside the Mod_Explicit `!$omp parallel` region.
    !$omp single
    if ( obj_rans % SpalartShur ) then   ! SA-RC: the velocity gradient and the correction terms first (MOSE_open order)
      call RC_Prepare ( domain )
      do b = 1, domain % nb
        if ( .not. ( is_local_block(b) ) ) cycle
        call SA_RC_Source_Blk ( SAPASS, domain % blk(b) % r, domain % blk(b) % m_f, domain % blk(b) % vol, &
                                domain % blk(b) % yn, domain % blk(b) % rho_c, domain % blk(b) % mil_c, &
                                domain % blk(b) % dtlocal, obj_rans % point_implicit, &
                                ROTATING_FRAME, obj_rot % Omega_vec, domain % blk(b) % dim )
      end do
    else if ( obj_rans % SAcomp ) then   ! SAcomp: the Paciorri-Sabetta correction of the production
      do b = 1, domain % nb
        if ( .not. ( is_local_block(b) ) ) cycle
        call SA_Comp_Source_Blk ( SAPASS, domain % blk(b) % r, domain % blk(b) % m_f, domain % blk(b) % vol, &
                                  domain % blk(b) % yn, domain % blk(b) % rho_c, domain % blk(b) % mil_c, &
                                  domain % blk(b) % dtlocal, obj_rans % point_implicit, obj_rans % SAR, &
                                  ROTATING_FRAME, obj_rot % Omega_vec, domain % blk(b) % dim )
      end do
    else
    do b = 1, domain % nb
      if ( .not. ( is_local_block(b) ) ) cycle
        call SA_Source_Blk ( SAPASS,                        &
                          domain % blk(b) % r,           &
                          domain % blk(b) % m_f,         &
                          domain % blk(b) % vol,         &
                          domain % blk(b) % yn,          &
                          domain % blk(b) % rho_c,       &
                          domain % blk(b) % mil_c,       &
                          domain % blk(b) % dtlocal, obj_rans % point_implicit, &
                          obj_rans % SAR,                &
                          ROTATING_FRAME,                &
                          obj_rot % Omega_vec,           &
                          domain % blk(b) % dim )
    end do
    end if
    !$omp end single

    end subroutine Spalart_Source_Terms



  !> Device SA source over the block box (the mirror of SST_Blk): one cell per thread, sa_body.inc.
  subroutine SA_Source_Blk ( SAARG, Res, Metric, Volume, WDist, rho_c, mil_c, dtl, pimp, &
                          sar, rotating_frame, omega_vec, n)
    use MOSE_Global_m
    use MOSE_Mod_Fluid, only: Strain_Tensor, Vorticity_Vector
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
    real(kind=8) :: fnu_
    logical, intent(in) :: sar, rotating_frame
    real(kind=8), intent(in) :: omega_vec(3)
    ! Local (per cell, fixed-size)
    integer      :: i, j, k, a, c
    real(kind=8) :: rho, mil, nit, chi, dist, Gtmp(3,3), Gradvel(3,3), vort(3), Om
    real(kind=8) :: fv1_, fv2_, Sbar, Stilde, Sij(3,3), S2, abs_Eij, Production
    real(kind=8) :: Gn(3), Gradnit(3), Gradnit2, Diffusion, r_Fun, g_Fun, fw, Destruction, Source
    real(kind=8) :: sum_im, sum_ip, sum_jm, sum_jp, sum_km, sum_kp
    integer :: nb1_, nb2_, nb3_
    real(kind=8) :: ov1_, ov2_, ov3_
    nb1_ = n(1); nb2_ = n(2); nb3_ = n(3)
    ov1_ = omega_vec(1); ov2_ = omega_vec(2); ov3_ = omega_vec(3)

    call field_block_index_params(b, boff, pi, pj)
    !$acc data copy(Res) copyin(Metric, Volume, WDist, rho_c, mil_c, dtl)
    !$acc parallel loop collapse(3) gang vector &
    !$acc&  present(P_flat, R_flat, Res, Metric, Volume, WDist, rho_c, mil_c, dtl, rho_c_f, mil_c_f) &
    !$acc&  private(Gtmp, Gradvel, vort, Sij, Gn, Gradnit)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)

#include "sa_body.inc"

    end do
    end do
    end do
    !$acc end data

  end subroutine SA_Source_Blk


  !> The SAcomp twin (MOSE_SACOMP_KERNEL): the production times the Paciorri-Sabetta factor.
  subroutine SA_Comp_Source_Blk ( SAARG, Res, Metric, Volume, WDist, rho_c, mil_c, dtl, pimp, &
                          sar, rotating_frame, omega_vec, n)
    use MOSE_Global_m
    use MOSE_Mod_Fluid, only: Strain_Tensor, Vorticity_Vector
    use FLINT_Lib_ThermoTransport_dev, only: f_ss_dev
    use MOSE_Field_m, only: P_flat, R_flat, ncell_total, field_block_index_params, rho_c_f, mil_c_f, Rg_c_f

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
    real(kind=8) :: fnu_
    logical, intent(in) :: sar, rotating_frame
    real(kind=8), intent(in) :: omega_vec(3)
    ! Local (per cell, fixed-size)
    integer      :: s
    real(kind=8) :: rhoi_(nsc), sound_   ! SAcomp: the species densities and the sound speed of the cell
    integer      :: i, j, k, a, c
    real(kind=8) :: rho, mil, nit, chi, dist, Gtmp(3,3), Gradvel(3,3), vort(3), Om
    real(kind=8) :: fv1_, fv2_, Sbar, Stilde, Sij(3,3), S2, abs_Eij, Production
    real(kind=8) :: Gn(3), Gradnit(3), Gradnit2, Diffusion, r_Fun, g_Fun, fw, Destruction, Source
    real(kind=8) :: sum_im, sum_ip, sum_jm, sum_jp, sum_km, sum_kp
    integer :: nb1_, nb2_, nb3_
    real(kind=8) :: ov1_, ov2_, ov3_
    nb1_ = n(1); nb2_ = n(2); nb3_ = n(3)
    ov1_ = omega_vec(1); ov2_ = omega_vec(2); ov3_ = omega_vec(3)

    call field_block_index_params(b, boff, pi, pj)
    !$acc data copy(Res) copyin(Metric, Volume, WDist, rho_c, mil_c, dtl)
    !$acc parallel loop collapse(3) gang vector &
    !$acc&  present(P_flat, R_flat, Res, Metric, Volume, WDist, rho_c, mil_c, dtl, rho_c_f, mil_c_f, Rg_c_f) &
    !$acc&  private(Gtmp, Gradvel, vort, Sij, Gn, Gradnit, rhoi_)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)

#define MOSE_SACOMP_KERNEL
#include "sa_body.inc"
#undef MOSE_SACOMP_KERNEL

    end do
    end do
    end do
    !$acc end data

  end subroutine SA_Comp_Source_Blk


  !> The -RC twin (MOSE_RC_KERNEL): the stored gradient and the correction factor fr1 on the production.
  subroutine SA_RC_Source_Blk ( SAARG, Res, Metric, Volume, WDist, rho_c, mil_c, dtl, pimp, &
                          rotating_frame, omega_vec, n)
    use MOSE_Global_m
    use MOSE_Mod_Fluid, only: Strain_Tensor, Vorticity_Vector
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
    real(kind=8) :: fnu_
    logical, intent(in) :: rotating_frame
    real(kind=8), intent(in) :: omega_vec(3)
    ! Local (per cell, fixed-size)
    real(kind=8) :: rstar, D2, rtilde, fr1   ! SA-RC (Spalart & Shur 1997)
    integer      :: i, j, k, a, c
    real(kind=8) :: rho, mil, nit, chi, dist, Gtmp(3,3), Gradvel(3,3), vort(3), Om
    real(kind=8) :: fv1_, fv2_, Sbar, Stilde, Sij(3,3), S2, abs_Eij, Production
    real(kind=8) :: Gn(3), Gradnit(3), Gradnit2, Diffusion, r_Fun, g_Fun, fw, Destruction, Source
    real(kind=8) :: sum_im, sum_ip, sum_jm, sum_jp, sum_km, sum_kp
    integer :: nb1_, nb2_, nb3_
    real(kind=8) :: ov1_, ov2_, ov3_
    nb1_ = n(1); nb2_ = n(2); nb3_ = n(3)
    ov1_ = omega_vec(1); ov2_ = omega_vec(2); ov3_ = omega_vec(3)

    call field_block_index_params(b, boff, pi, pj)
    !$acc data copy(Res) copyin(Metric, Volume, WDist, rho_c, mil_c, dtl)
    !$acc parallel loop collapse(3) gang vector &
    !$acc&  present(P_flat, R_flat, gv_f, rc1_f, rc2_f, Res, Metric, Volume, WDist, rho_c, mil_c, dtl, rho_c_f, mil_c_f) &
    !$acc&  private(Gtmp, Gradvel, vort, Sij, Gn, Gradnit)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)

#define MOSE_RC_KERNEL
#include "sa_body.inc"
#undef MOSE_RC_KERNEL

    end do
    end do
    end do
    !$acc end data

  end subroutine SA_RC_Source_Blk






  pure function fv1 (chi) result (result)
    !$acc routine seq
    implicit none
    real(R8), intent(in) :: chi
    real(R8) :: result
    result = ( chi**3 ) / ( chi**3 + cv1**3 )
  end function fv1

  function fv2 (chi) result (result)
    !$acc routine seq
    implicit none
    real(R8) :: result
    real(R8), intent(in) :: chi
    result = 1d0 - chi / ( 1d0 + chi * fv1(chi) )
  end function fv2


  !> The residual of the Mach equation of the correction, y^2 (0.44/(1 + 14 y^5) + 0.56) - rr (MOSE_open func).
  function func (yy, rr) result (result)
    !$acc routine seq
    implicit none
    real(R8), intent(in) :: yy, rr
    real(R8) :: result, ff
    ff = 0.44d0 / ( 1d0 + 14d0 * yy**5d0 ) + 0.56d0
    result = yy * yy * ff - rr
  end function func


  !> Compressibility correction Paciorri-Sabetta https://doi.org/10.2514/2.3967 (MOSE_open Lib_Spalart
  !> Compressibility_Correction, the sound speed computed by the caller): Production times f_cor1 f_cor2 of the
  !> turbulent Mach number. The iteration is MOSE_open's (secant derivative, |step| < tolerance = 1d-2); on the device it
  !> stops after 1000 steps as well, where MOSE_open has no cap (a cell that does not converge would hang the kernel).
  subroutine SA_Compressibility_Correction ( nit, chi, sound, om, prod )
    !$acc routine seq
    implicit none
    real(R8), intent(in) :: nit, chi, sound, om
    real(R8), intent(inout) :: prod
    ! Local
    integer  :: it
    real(R8) :: c_mach, staux, err, deltamach
    real(R8) :: aux_num, aux_den, c_plus, c_minus, res, f_cor1, f_cor2, f_cor

    c_mach = 1.d0
    staux = 25.d0*(om*nit*fv1(chi))/sound
    err = 1.d0
    it = 0
    do while ( err > tolerance .and. it < 1000 )
      deltamach = .1d0*c_mach
      aux_num = func(c_mach,staux)*2.*deltamach
      c_plus = c_mach+deltamach
      c_minus = c_mach-deltamach
      aux_den = func(c_plus,staux)-func(c_minus,staux)
      res = aux_num/aux_den
      err = dabs(res)
      c_mach = c_mach-res
      it = it + 1
    end do
    ! correction function
    f_cor1 = 0.6d0*(1.d0/(1.d0+9.d0*(c_mach)**6.))+0.4d0
    f_cor2 = 0.44d0*(1.d0/(1.d0+14.d0*(c_mach)**5.))+0.56d0
    f_cor = f_cor1*f_cor2
    prod = prod*f_cor

  end subroutine SA_Compressibility_Correction

  subroutine Spalart_Set_Wall_Values ( mi_l, mitilde, dist )
    use MOSE_Global_m

    implicit none
    real(R8), intent(out), dimension(nRANS)  :: mitilde          ! : rho*nit at wall
    real(R8), intent(in)                     :: mi_l             ! : laminar viscosity at wall
    real(R8), intent(in)                     :: dist             ! : cell center distance wall

    mitilde = 0d0 ! Smooth solid surface condition

  end subroutine Spalart_Set_Wall_Values


  subroutine Spalart_Set_Blowing_Wall ( rho, mil, tau, mitilde, mdot, dist )
    use MOSE_Global_m

    implicit none
    real(R8), intent(out), dimension(nRANS)  :: mitilde
    real(R8), intent(in)                     :: rho, mil, mdot, dist
    real(R8), dimension(3), intent(in)       :: tau

    mitilde = 0d0 ! Smooth solid surface condition

  end subroutine Spalart_Set_Blowing_Wall


  subroutine Spalart_Extrapolate_Wall ( prim, wall_prim, rho, wall_rho, ghost_prim )
    use MOSE_Global_m

    implicit none
    real(R8), intent(in),  dimension(nRANS) :: prim, wall_prim  ! : boundary cell and wall variables
    real(R8), intent(in) :: rho, wall_rho ! : boundary cell density and wall density
    real(R8), intent(out), dimension(nRANS) :: ghost_prim ! : ghost cell RANS variables

    ghost_prim(1) = - prim(1) ! enforcing mit_wall == 0

  end subroutine Spalart_Extrapolate_Wall


  subroutine Spalart_Enforce_Realizability ( mitilde )

    implicit none
    real(R8), intent(inout),  dimension(:) :: mitilde
    if ( mitilde(1) < 1d-10 ) mitilde(1) = 1d-10

  end subroutine Spalart_Enforce_Realizability


  !> Eddy viscosity mu_t = fv1(chi) * rho*nu-tilde. `acc routine seq`: called from the dt kernels and the wall BC
  !> kernels on the device (dispatched by eddy_dev = 2 next to SST_Eddy_Viscosity).
  subroutine Spalart_Eddy_Viscosity ( mut, mitilde, milam, rho, Gradvel, dist )
    !$acc routine seq
    use MOSE_Global_m

    implicit none
    real(R8), intent(in), dimension(nRANS) :: mitilde        ! : vector of RANS eqs variables in form (rho*nitilde in this case)
    real(R8), intent(in)                   :: milam          ! : Laminar viscosity
    real(R8), intent(in)                   :: rho            ! : Density
    real(R8), intent(in)                   :: dist           ! : wall distance (not used)
    real(R8), intent(in), dimension(3,3)   :: Gradvel        ! : Velocity gradient
    real(R8), intent(out)                  :: mut            ! : Eddy viscosity
    real(R8) :: chi

    chi = mitilde(1)/milam
    mut = fv1( chi )*mitilde(1)

  end subroutine Spalart_Eddy_Viscosity


  !> Diffusive flux of rho*nu-tilde: area*(mu_l + rho*nu-tilde)*(grad nu-tilde . n)/sigma. `acc routine seq`.
  subroutine Spalart_RANS_Diffusive_Flux (flux, mitilde, milam, gradvel, gradnit, rho, area, normal, dist)
    !$acc routine seq
    use MOSE_Global_m

    implicit none
    real(R8), intent(in), dimension(nrans)   :: mitilde               ! : vector of RANS eqs, in this case (rho*nitilde)
    real(R8), intent(in)                     :: milam                 ! : laminar viscosity
    real(R8), intent(in), dimension(3,3)     :: gradvel               ! : Velocity gradient (not used)
    real(R8), intent(in), dimension(nrans,3) :: gradnit               ! : Gradient in form dnitilde/dxj
    real(R8), intent(in)                     :: rho                   ! : Density
    real(R8), intent(in)                     :: area, normal(3), dist ! : Metrics
    real(R8), intent(out), dimension(nrans)  :: flux                  ! : Diffusive flux

    ! explicit dot product: a strided row of the (nrans,3) dummy inside an acc routine is misread by nvfortran
    flux(1) = area*(milam + mitilde(1))*( gradnit(1,1)*normal(1) + gradnit(1,2)*normal(2) + gradnit(1,3)*normal(3) )/sigma

  end subroutine Spalart_RANS_Diffusive_Flux


end module MOSE_Mod_Spalart
