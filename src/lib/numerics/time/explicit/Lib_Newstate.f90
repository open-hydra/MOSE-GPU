module MOSE_Lib_Newstate
  ! ieee_is_nan instead of the non-standard isnan() extension (nvfortran has no isnan)
  use, intrinsic :: ieee_arithmetic
  use iso_fortran_env, only: R8 => real64
  use MOSE_Global_m,   only: gc, nprim, NPRIM_DEV, nsc, nu, nv, nw, np, nt, nrans, RANS
  use MOSE_Field_m,    only: P_flat, PO_flat, R_flat, field_block_index_params
  implicit none

  !> Concrete procedure pointing to one of the function realizations
  procedure(newstate_if), pointer, public :: RK_Newstate

  !> [MOSE-scheme] state-check = t | f (default t, as MOSE_open): after every state update stop on a negative pressure or a
  !> NaN primitive, as MOSE_open Check_And_Fix_State does. Device path only.
  logical, public :: state_check = .true.

  !> Abstract interface relative to the "Newstate" procedure
  abstract interface
    subroutine Newstate_if ( domain, irk )
      use MOSE_Types_m
      use MOSE_Global_m
      use FLINT_Lib_Thermodynamic
      use MOSE_Lib_RK
      use MOSE_Lib_Strang
        use MOSE_Lib_RANS
      
      implicit none
      type(MOSE_domain_type), intent(inout) :: domain 
      integer, intent(in)                   :: irk
      ! Local
      integer :: i, j, k, s, b

    end subroutine Newstate_if
  end interface

  public :: Newstate_Conservative_blk
  public :: Newstate_Primitive_blk
  public :: State_Check_blk, State_Check_report

contains

  
  subroutine Assign_Integration_Variables ( intvar )

    implicit none
    character(len=*), intent(inout)  :: intvar

    select case (intvar)
    case ('prim')
      RK_Newstate => Newstate_Primitive
      intvar = 'primitive'
    case ('cons')
      RK_Newstate => Newstate_Conservative
      intvar = 'conservative'
    case ('prec')   ! MOSE_open Assign_Integration_Variables: the primitive update with the preconditioned Theta
      RK_Newstate => Newstate_Primitive
      intvar = 'preconditioned'
    end select

  end subroutine Assign_Integration_Variables


  !> [MOSE-scheme] state-check = t | f, read by IO_Ini at start-up and at every re-read of input.ini.
  subroutine Assign_State_Check ( key_in )
    implicit none
    character(len=*), intent(in) :: key_in
    select case ( trim(adjustl(key_in)) )
      case ('t', 'T', 'true', '.true.', 'TRUE', 'yes'); state_check = .true.
      case ('f', 'F', 'false', '.false.', 'FALSE', 'no'); state_check = .false.
      case default
        write(*,'(3A)') ' [ERROR] [MOSE-scheme] state-check = ', trim(key_in), ': t or f'
        error stop 'Assign_State_Check: unknown state-check'
    end select
  end subroutine Assign_State_Check


  subroutine Newstate_Conservative ( domain, irk )
    use MOSE_Types_m
    use MOSE_Global_m
    use MOSE_Lib_RK
    use MOSE_Lib_Strang
    use MOSE_Lib_IRS
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Lib_IRS, only: IRS_apply

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in)                   :: irk
    ! Local
    integer :: b
    integer :: nbad_, ib_, jb_, kb_

    !> The state update on the resident P_flat, one kernel per local block (Mod_Newstate_W); IRS on the device too.
    ! IRS on the device (MOSE_open Newstate_Conservative, IRS branch): the resident P_flatr is scaled by
    ! -dt/vol*strang and smoothed (Mod_IRS_W) first, and the state update below takes it pre-scaled.
    if ( obj_irs%IRS ) then
      !$omp single
      call IRS_apply ( domain, obj_irs%beta, obj_irs%njacobi, strangcoeff )
      !$omp end single
    end if
    !$omp single
    do b = 1, domain % nb
      if ( .not. ( is_local_block(b) ) ) cycle
        call Newstate_Conservative_blk ( b, domain%blk(b)%dtlocal,      &
               domain%blk(b)%vol, domain%blk(b)%T_c,                        &
               domain%blk(b)%dim, strangcoeff, rk_coeff(irk), prescaled = obj_irs%IRS )
    end do
    if ( state_check ) then   ! the state check of MOSE_open (state-check key): a separate launch
      !$acc wait
      do b = 1, domain % nb
        if ( .not. is_local_block(b) ) cycle
        call State_Check_blk ( b, domain%blk(b)%dim, nbad_, ib_, jb_, kb_ )
        if ( nbad_ > 0 ) call State_Check_report ( b, ib_, jb_, kb_, nbad_ )
      end do
    end if
    !$omp end single

  end subroutine Newstate_Conservative


  subroutine Newstate_Primitive ( domain, irk )
    use MOSE_Types_m
    use MOSE_Lib_Preconditioning, only: obj_prec, ur_f   ! allocated by Compute_dt (Prec_Alloc), first in every step
    use MOSE_Global_m
    use MOSE_Lib_RK
    use MOSE_Lib_Strang
    use MOSE_Lib_IRS
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Lib_IRS, only: IRS_apply

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in)                   :: irk
    ! Local
    integer :: b
    integer :: nbad_, ib_, jb_, kb_

    !> MOSE_open Newstate_Primitive on the resident P_flat, one kernel per local block. With residual smoothing MOSE_open
    !> smooths the raw conservative residual (Residual_Smoothing, no scaling) and then converts it (Compute_Prim_Residual).
    if ( obj_irs%IRS ) then
      !$omp single
      call IRS_apply ( domain, obj_irs%beta, obj_irs%njacobi, strangcoeff, scale = .false. )
      !$omp end single
    end if
    !$omp single
    do b = 1, domain % nb
      if ( .not. is_local_block(b) ) cycle
      call Newstate_Primitive_blk ( b, domain%blk(b)%dtlocal, domain%blk(b)%vol, domain%blk(b)%dim, &
                                    strangcoeff, rk_coeff(irk), ur_f, size(ur_f), obj_prec%enabled )
    end do
    if ( state_check ) then   ! the state check of MOSE_open (Check_And_Fix_State: NaN or p < 0 stops), a separate launch
      !$acc wait
      do b = 1, domain % nb
        if ( .not. is_local_block(b) ) cycle
        call State_Check_blk ( b, domain%blk(b)%dim, nbad_, ib_, jb_, kb_ )
        if ( nbad_ > 0 ) call State_Check_report ( b, ib_, jb_, kb_, nbad_ )
      end do
    end if
    !$omp end single

  end subroutine Newstate_Primitive






  ! ---- from Mod_Newstate_W.f90 (MOSE_Lib_Newstate) ----

  !> [MOSE-scheme] state-check = t (set by the MOSE_open input dialect): MOSE_open Check_And_Fix_State after every state
  !> update. Counts the interior cells of block b with a negative pressure or a NaN (or infinite) primitive and returns
  !> the first one in memory order; a separate launch, the update kernels are untouched. NaN is caught by ordered
  !> comparisons (false for NaN), not by x /= x, which fast math may fold away. The species floor 1e-20 of the MOSE_open
  !> routine is already in the update (newstate_cons_body.inc).
  subroutine State_Check_blk ( b, n, nbad, ibad, jbad, kbad )
    integer, intent(in)  :: b, n(3)
    integer, intent(out) :: nbad, ibad, jbad, kbad
    integer :: i, j, k, s, c0, cmin, boff, pi, pj, n1, n2, n3
    logical :: bad
    call field_block_index_params(b, boff, pi, pj)
    n1 = n(1); n2 = n(2); n3 = n(3)
    nbad = 0; cmin = huge(1)
    !$acc parallel loop collapse(3) gang vector present(P_flat) reduction(+:nbad) reduction(min:cmin) private(c0, s, bad)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      bad = .not. ( P_flat(c0,np) >= 0.0_R8 )
      do s = 1, nprim
        if ( .not. ( P_flat(c0,s) < huge(1.0_R8) .and. P_flat(c0,s) > -huge(1.0_R8) ) ) bad = .true.
      end do
      if ( bad ) then
        nbad = nbad + 1
        cmin = min(cmin, c0)
      end if
    end do
    end do
    end do
    ibad = 0; jbad = 0; kbad = 0
    if ( nbad > 0 ) then
      c0 = cmin - boff - 1; kbad = c0 / (pi*pj); c0 = c0 - kbad*pi*pj; jbad = c0 / pi; ibad = c0 - jbad*pi
      ibad = ibad - gc + 1; jbad = jbad - gc + 1; kbad = kbad - gc + 1
    end if
  end subroutine State_Check_blk

  !> The stop of MOSE_open Check_And_Fix_State, with its message, for the cell found by State_Check_blk.
  subroutine State_Check_report ( b, i, j, k, nbad )
    integer, intent(in) :: b, i, j, k, nbad
    integer :: c0, boff, pi, pj
    real(R8) :: wl(nprim)
    call field_block_index_params(b, boff, pi, pj)
    c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
    !$acc update host(P_flat(c0:c0, 1:nprim))
    wl = P_flat(c0, 1:nprim)
    write(*,'(A,4I4)') "Integration failed at b, i, j, k:", b, i, j, k
    write(*,*) wl(1:nprim)
    write(*,'(A,I0,A)') " (", nbad, " cells of this block with p < 0 or a NaN primitive; the first one is printed)"
    error stop "NaN or p<0 detected"
  end subroutine State_Check_report

  !> The conservative Newstate for ONE block on var-last P_flat/PO_flat/R_flat. Reads PO_flat (old state)
  !> and R_flat (conservation-form residual), updates w (new primitive state) in place.
  !>   b        : 1-based global block index (selects the P_flat/PO_flat/R_flat slice via cidx).
  !>   dtlocal  : local timestep per interior cell.
  !>   vol      : cell volume (padded).
  !>   T_c      : start-of-stage temperature cache (padded) -> Newton-Raphson seed.
  !>   n        : block interior dims (ni,nj,nk).
  !>   strang   : Strang operator-splitting coefficient.
  !>   rkc      : RK-stage coefficient rk_coeff(irk).
  subroutine Newstate_Conservative_blk ( b, dtlocal, vol, T_c, n, strang, rkc, prescaled )
    use FLINT_Lib_Thermodynamic, only: h_tab, cp_tab, Ri_tab, f_Rtot, Tmin, Tmax
    integer,  intent(in) :: b
    integer,  intent(in) :: n(3)
    real(R8), intent(in) :: strang, rkc
    real(R8), dimension(1:n(1), 1:n(2), 1:n(3)), intent(in) :: dtlocal
    logical,  intent(in), optional :: prescaled   ! IRS: R_flat already holds -dt/vol*strang*R smoothed
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: vol, T_c
    ! Local
    integer  :: i, j, k, s, T_i, it, c0
    integer  :: boff, pi, pj
    real(R8) :: rho, Rg, vmag, Tloc, Td, hmix, E0, dt
    real(R8) :: energy, p_it, cv, enth, diff_en, T_it, h_s, cp_s
    real(R8) :: dt_l, vol_l, tc_l   ! 5s-a: the per-cell inputs the shared body reads
    logical  :: ps_              ! IRS: the residual arrives scaled and smoothed (Mod_IRS_W = MOSE_open Update_State_Cons_IRS)
    integer  :: tmin_, tmax_     ! table bounds: the T index is clamped, an off-table state stops in the thermo guard
    ! Per-cell gathered state (contiguous): nprim is a compile-time parameter, exact-size locals.
    real(R8) :: p_l(nprim), po_l(nprim)
    real(R8) :: consold(nprim), cons(nprim), res_l(nprim)
    ! Dims as scalars -- a dummy array referenced inside the
    ! acc region is implicitly copied per launch; scalars ride the launch parameters.
    integer :: nb1_, nb2_, nb3_
    nb1_ = n(1); nb2_ = n(2); nb3_ = n(3)
    ! The resident var-LAST R_flat is the SOLE residual input. Device
    ! writers (conv/diff pass-2, BC, SST) accumulate into it directly; host bridge
    ! sources (RotatingFrame/soot) are folded in by field_add_R_to_R_flat BEFORE this
    ! kernel (Mod_Explicit), so the var-first r never enters the state update.
    call field_block_index_params(b, boff, pi, pj)

    ps_ = .false.; if (present(prescaled)) ps_ = prescaled
    tmin_ = Tmin; tmax_ = Tmax
    !$acc data copyin(dtlocal, vol, T_c)
    !$acc parallel loop collapse(3) gang vector present(P_flat, PO_flat, R_flat, dtlocal, vol, T_c) &
    !$acc&  private(p_l, po_l, consold, cons, res_l, c0, T_i, it, s, dt_l, vol_l, tc_l)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)
      c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      dt_l = dtlocal(i,j,k); vol_l = vol(i,j,k); tc_l = T_c(i,j,k)
#include "newstate_cons_body.inc"
    end do; end do; end do
    !$acc end data

  end subroutine Newstate_Conservative_blk





  !> MOSE_open Newstate_Primitive for ONE block on the resident P_flat / PO_flat / R_flat, cell by cell:
  !> Compute_Prim_Residual (the conservative residual R_flat turned into the residual of the primitive variables, Weiss-Smith
  !> form, Theta = (drho/dp)_T with prim; with prec, Theta = 1/(beta2 a^2) + 1/(cp T), beta2 = max(Ur^2/a^2, 1d-6)), the
  !> scaling dt/vol*strang, RK_stage, and the floors of Check_And_Fix_State (species 1e-20, the RANS realizability 1e-10).
  !> The thermodynamics as MOSE_open (co_rotot_Rtot, f_cp_expr, H0, f_tabT_expr) on the device tables, T unclamped.
  subroutine Newstate_Primitive_blk ( b, dtlocal, vol, n, strang, rkc, urs, nur, prec )
    use FLINT_Lib_Thermodynamic,       only: h_tab, Ri_tab
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, f_cp_dev, H0_dev
    integer,  intent(in) :: b
    integer,  intent(in) :: n(3)
    real(R8), intent(in) :: strang, rkc
    integer,  intent(in) :: nur
    real(R8), intent(in) :: urs(nur)   ! Lib_Preconditioning ur_f: Ur of the cells (prec), 1 element otherwise
    logical,  intent(in) :: prec
    real(R8), dimension(1:n(1), 1:n(2), 1:n(3)), intent(in) :: dtlocal
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: vol
    integer, parameter :: ns1 = max(nsc-1, 1)
    ! Local
    integer  :: i, j, k, s, r_, c_, T_i, Tint(2), c0, boff, pi, pj
    real(R8) :: rho, Rgas, temperature, Tdiff, inv_rho, modvel, q2, cp, gam, sound, a2, rho_T, rho_p, enthalpy
    real(R8) :: Theta, Denom, drho, drhoN, dyN, dt_, vol_, sYd, sdY, sR1, sR2, sv, sres, sdrdy, beta2
    real(R8) :: p_l(nprim), po_l(nprim), res_l(nprim), vel(3), invG(5,5), cres(5), pres(5)
    real(R8) :: Yi(ns1), drdyi(ns1), dhdyi(ns1), dyi(ns1), Rs(np), dPa(np)
    call field_block_index_params(b, boff, pi, pj)

    !$acc data copyin(dtlocal, vol)
    !$acc parallel loop collapse(3) gang vector present(P_flat, PO_flat, R_flat, dtlocal, vol, urs) &
    !$acc&  private(p_l, po_l, res_l, vel, invG, cres, pres, Yi, drdyi, dhdyi, dyi, Rs, dPa, Tint, c0, T_i, s, r_, c_)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)
      c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      do s = 1, nprim
        p_l(s)   = P_flat(c0,s)
        po_l(s)  = PO_flat(c0,s)
        res_l(s) = R_flat(c0,s)
      end do
      dt_ = dtlocal(i,j,k); vol_ = vol(i,j,k)

      ! ---- Compute_Prim_Residual
      call co_rotot_Rtot_dev ( p_l(1:nsc), rho, Rgas )
      temperature = p_l(np)/(rho*Rgas)
      T_i     = int(temperature)
      Tdiff   = temperature - T_i
      Tint(1) = T_i
      Tint(2) = T_i + 1
      inv_rho  = 1.d0/rho
      vel(1) = p_l(nu); vel(2) = p_l(nv); vel(3) = p_l(nw)
      modvel   = norm2(vel)
      q2       = modvel*modvel
      cp       = f_cp_dev ( p_l(1:nsc), Tint, Tdiff, rho )
      gam      = cp/(cp - Rgas)
      sound    = sqrt( gam*Rgas*temperature )
      a2       = sound*sound
      rho_T    = -rho/temperature
      rho_p    = 1.d0/(Rgas*temperature)
      enthalpy = H0_dev ( p_l(np), p_l(1:nsc), modvel )
      if ( prec ) then   ! Weiss-Smith preconditioning: Theta from the pseudo-acoustic speed Ur of the cell
        beta2 = max( urs(c0)*urs(c0)/a2, 1.d-6 )
        Theta = 1.d0/(beta2*a2) + 1.d0/(cp*temperature)
      else
        Theta = rho_p
      endif
      Denom = 1.d0/( rho*cp*Theta + rho_T )

      if ( nsc == 1 ) then
        ! invGamma : (drho, drhoV, drhoE) -> (dp, du, dv, dw, dT)
        invG(1,1) = ( rho_T*(enthalpy - q2) + rho*cp )*Denom
        invG(1,2) = rho_T*vel(1)*Denom
        invG(1,3) = rho_T*vel(2)*Denom
        invG(1,4) = rho_T*vel(3)*Denom
        invG(1,5) = -rho_T*Denom
        invG(2,1) = -vel(1)*inv_rho
        invG(2,2) =  inv_rho
        invG(2,3) =  0.d0
        invG(2,4) =  0.d0
        invG(2,5) =  0.d0
        invG(3,1) = -vel(2)*inv_rho
        invG(3,2) =  0.d0
        invG(3,3) =  inv_rho
        invG(3,4) =  0.d0
        invG(3,5) =  0.d0
        invG(4,1) = -vel(3)*inv_rho
        invG(4,2) =  0.d0
        invG(4,3) =  0.d0
        invG(4,4) =  inv_rho
        invG(4,5) =  0.d0
        invG(5,1) = -( Theta*(enthalpy - q2) - 1.d0 )*Denom
        invG(5,2) = -vel(1)*Theta*Denom
        invG(5,3) = -vel(2)*Theta*Denom
        invG(5,4) = -vel(3)*Theta*Denom
        invG(5,5) =  Theta*Denom
        ! Cons_Res = -[ sum(residual(1:nsc)), residual(nu:nw), residual(np) ];  Prec_Res = matmul(invGamma, Cons_Res)
        cres(1) = -res_l(1); cres(2) = -res_l(nu); cres(3) = -res_l(nv); cres(4) = -res_l(nw); cres(5) = -res_l(np)
        do r_ = 1, 5
          pres(r_) = 0.d0
          do c_ = 1, 5
            pres(r_) = pres(r_) + invG(r_,c_)*cres(c_)
          end do
        end do
        res_l(nu) = pres(2); res_l(nv) = pres(3); res_l(nw) = pres(4)
        res_l(np) = pres(1)
        ! density from the dp and dT (physical EOS derivatives)
        res_l(1)  = rho_p*pres(1) + rho_T*pres(5)
      else
        ! Multi-species: T and Rgas depend on the composition, the mass-fraction sensitivities enter the p and T rows
        do s = 1, nsc-1
          Yi(s) = p_l(s)*inv_rho
          if ( Yi(s) <= 1.d-9 ) Yi(s) = 0.d0
          drdyi(s) = -rho*( Ri_tab(s) - Ri_tab(nsc) )/Rgas
          dhdyi(s) = ( h_tab(Tint(1),s) + (h_tab(Tint(2),s) - h_tab(Tint(1),s))*Tdiff ) &
                   - ( h_tab(Tint(1),nsc) + (h_tab(Tint(2),nsc) - h_tab(Tint(1),nsc))*Tdiff )
        end do
        ! Res = -[ sum(residual(1:nsc)), residual(nu:np), residual(1:nsc-1) ]
        sres = 0.d0
        do s = 1, nsc
          sres = sres + res_l(s)
        end do
        Rs(1) = -sres
        Rs(2) = -res_l(nu); Rs(3) = -res_l(nv); Rs(4) = -res_l(nw); Rs(5) = -res_l(np)
        do s = 1, nsc-1
          Rs(5+s) = -res_l(s)
        end do
        sYd = 0.d0; sdY = 0.d0; sR1 = 0.d0; sR2 = 0.d0
        do s = 1, nsc-1
          sYd = sYd + Yi(s)*drdyi(s)                                  ! sum(Yi*drdyi)
          sdY = sdY + dhdyi(s)*Yi(s)                                  ! sum(dhdyi*Yi)
          sR1 = sR1 + Rs(5+s)*(rho_T*dhdyi(s) - cp*drdyi(s))          ! sum(Res(6:np)*(rho_T*dhdyi - cp*drdyi))
          sR2 = sR2 + Rs(5+s)*( Theta*dhdyi(s) + drdyi(s)*inv_rho )   ! sum(Res(6:np)*(Theta*dhdyi + drdyi*inv_rho))
        end do
        sv = vel(1)*Rs(2) + vel(2)*Rs(3) + vel(3)*Rs(4)                ! sum(velocity*Res(2:4))
        dPa(1) = ( Rs(1)*( cp*(rho + sYd) - rho_T*( sdY + q2 - enthalpy ) ) &
                 + sR1 - rho_T*Rs(5) + rho_T*sv )*Denom
        dPa(2) = -Rs(1)*vel(1)*inv_rho + Rs(2)*inv_rho
        dPa(3) = -Rs(1)*vel(2)*inv_rho + Rs(3)*inv_rho
        dPa(4) = -Rs(1)*vel(3)*inv_rho + Rs(4)*inv_rho
        dPa(5) = ( Rs(1)*( Theta*( sdY - enthalpy + q2 ) + sYd*inv_rho + 1.d0 ) &
                 - sR2 + Theta*( Rs(5) - sv ) )*Denom
        do s = 1, nsc-1
          dPa(5+s) = ( Rs(s+5) - Rs(1)*Yi(s) )*inv_rho
        end do
        sdrdy = 0.d0; dyN = 0.d0
        do s = 1, nsc-1
          dyi(s) = dPa(5+s)
          sdrdy = sdrdy + drdyi(s)*dyi(s)                             ! sum(drdyi*dyi)
          dyN   = dyN + dyi(s)                                        ! sum(dyi)
        end do
        drho  = rho_p*dPa(1) + rho_T*dPa(5) + sdrdy
        dyN   = -dyN
        drhoN = rho*dyN + p_l(nsc)*inv_rho*drho
        do s = 1, nsc-1
          res_l(s) = rho*dyi(s) + Yi(s)*drho                          ! drhoi
        end do
        res_l(nsc) = drhoN
        res_l(nu) = dPa(2); res_l(nv) = dPa(3); res_l(nw) = dPa(4)
        res_l(np) = dPa(1)
      end if
      ! Turbulence variables (not preconditioned)
      if (nprim > np) then
        do s = np+1, nprim
          res_l(s) = -res_l(s)
        end do
      end if
      ! Time scaling: residual = residual * dt / volume * strangcoeff
      do s = 1, nprim
        res_l(s) = res_l(s) * dt_ / vol_ * strang
      end do

      ! ---- RK_stage: newstate = oldstate + (state - oldstate + residual) * rkcoeff
      do s = 1, nprim
        p_l(s) = po_l(s) + ( p_l(s) - po_l(s) + res_l(s) ) * rkc
      end do
      ! ---- Check_And_Fix_State: the floors (the NaN / p < 0 stop is State_Check_blk, after the launch)
      do s = 1, nsc
        if ( p_l(s) < 1d-20 ) p_l(s) = 1d-20
      end do
      if (RANS) then   ! SST_ / Wilcox2006_ / Spalart_Enforce_Realizability on the nrans turbulence variables
        do s = nt, nt + nrans - 1
          if ( p_l(s) < 1d-10 ) p_l(s) = 1d-10
        end do
      end if

      do s = 1, nprim
        P_flat(c0,s) = p_l(s)
      end do
    end do; end do; end do
    !$acc end data

  end subroutine Newstate_Primitive_blk
end module MOSE_Lib_Newstate
