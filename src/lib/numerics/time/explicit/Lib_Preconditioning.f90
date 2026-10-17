!> @brief Low-Mach preconditioning (integration-variables = prec) on the device, as MOSE_open b29130b Lib_Preconditioning.
!> The reference velocity Ur of every cell (Update_Derived_Variables, comp_Ur) enters the preconditioned primitive update
!> (Theta of Compute_Prim_Residual, Lib_Newstate), the time step (c* = sqrt(alpha^2 u^2 + Ur^2), Mod_dt) and HLLC-PC.
!>
!> Ur modes, in the order of MOSE_open:
!>   1. preconditioning-Uref > 0: Ur = Uref, the sound speed where the velocity exceeds it;
!>   2. preconditioning-Mach > 0: the Ur of the target preconditioned Mach number (the quadratic of MOSE_open), the sound
!>      speed where the target cannot be reached, at least eps_min*a, at most a, at least the viscous Ur;
!>   3. otherwise Weiss-Smith: Ur = |u| Ur-factor, eps_min*a below eps_min*a, a above a, at least the viscous Ur.
!> Then preconditioning-Ur-min as a floor, and preconditioning-Ur-smooth passes of max over the six face neighbours.
!> The viscous Ur (Navier-Stokes) is max(mu_eff/(rho dx_min), k_eff/(rho cp dx_min)), mu_eff = mu_l + mu_t and
!> k_eff = k_l + mu_t cp/Prt with RANS; mu_l and k_l from the Wilke mixture rule (not [MOSE-Physics] prandtl), mu_t from
!> the eddy viscosity of the model with the stored velocity gradient of the -RC models and a zero gradient otherwise:
!> MOSE_open passes blk%vel_gradient, which only Compute_Velocity_Gradient (the -RC sources) writes.
!>
!> Device layout: ur_f(row) on the flat rows of P_flat, written on the interior cells only. The ghost rows stay 0 as
!> MOSE_open blk%Ur, which no boundary condition and no exchange fills: HLLC-PC at a connection face reads 0 on the ghost
!> side in both codes (a MOSE_open defect, reported to the owner).
module MOSE_Lib_Preconditioning
  use iso_fortran_env, only: R8 => real64
  implicit none
  private
  public :: obj_prec, ur_f, Prec_Alloc, Update_Derived_Variables, comp_Ur

  type :: prec_type
    logical  :: enabled     = .false.   !< integration-variables = prec
    real(R8) :: eps_min     = 0.05_R8   !< [MOSE-scheme] preconditioning-eps-min (MOSE_open default 0.05)
    real(R8) :: Ur_ref      = 0.0_R8    !< preconditioning-Uref
    real(R8) :: Ur_min      = 0.0_R8    !< preconditioning-Ur-min
    real(R8) :: Mach_target = 0.0_R8    !< preconditioning-Mach
    real(R8) :: Ur_factor   = 1.0_R8    !< preconditioning-Ur-factor
    integer  :: n_smooth_Ur = 0         !< preconditioning-Ur-smooth
    contains
      procedure, pass(self) :: ini_setup => read_ini
  end type
  type(prec_type), save :: obj_prec

  real(R8), allocatable, save :: ur_f(:), urb_f(:)
  real(R8), allocatable, save :: g0_f(:,:)   !< the zero gradient (1 row) for the models without the -RC gradient
  logical, save :: ur_ready = .false.   !< .true. once Ur has been computed (MOSE_open Wrap_Setup does it before step 1)

#  define WR(i,j,k) ( boff + ((i)+gc-1) + ((j)+gc-1)*pi + ((k)+gc-1)*pi*pj + 1 )

contains


  subroutine read_ini ( self, fini )
    use Finer, only: file_ini
    use MOSE_Parameters_m, only: codename
    implicit none
    class(prec_type) :: self
    type(file_ini) :: fini
    ! Local
    integer :: error
    character(len=32) :: iv
    logical :: ta

    iv = 'cons'
    call fini%get(section_name=trim(codename)//'-scheme', option_name='integration-variables', val=iv, error=error)
    self%enabled = ( trim(adjustl(iv)) == 'prec' )
    call fini%get(section_name=trim(codename)//'-scheme', option_name='preconditioning-eps-min', val=self%eps_min, error=error)
    if (error/=0) self%eps_min = 0.05_R8
    call fini%get(section_name=trim(codename)//'-scheme', option_name='preconditioning-Uref', val=self%Ur_ref, error=error)
    if (error/=0) self%Ur_ref = 0.0_R8
    call fini%get(section_name=trim(codename)//'-scheme', option_name='preconditioning-Ur-min', val=self%Ur_min, error=error)
    if (error/=0) self%Ur_min = 0.0_R8
    call fini%get(section_name=trim(codename)//'-scheme', option_name='preconditioning-Mach', val=self%Mach_target, error=error)
    if (error/=0) self%Mach_target = 0.0_R8
    call fini%get(section_name=trim(codename)//'-scheme', option_name='preconditioning-Ur-factor', val=self%Ur_factor, error=error)
    if (error/=0) self%Ur_factor = 1.0_R8
    call fini%get(section_name=trim(codename)//'-scheme', option_name='preconditioning-Ur-smooth', val=self%n_smooth_Ur, error=error)
    if (error/=0) self%n_smooth_Ur = 0
    if ( self%enabled .and. .not. ( self%eps_min > 0d0 .and. self%eps_min < 1d0 ) ) &
      error stop 'preconditioning-eps-min must be in (0, 1) (MOSE_open Register_Numerics)'
    ta = .false.
    call fini%get(section_name=trim(codename)//'-parameters', option_name='time-accurate', val=ta, error=error)
    if ( self%enabled .and. ta ) &   ! MOSE_open Assign_Setup
      error stop '[ERROR] integration-variables=prec is currently supported only for steady/pseudo-time runs.'

  end subroutine read_ini


  !> Ur of one cell (MOSE_open comp_Ur), the parameters passed as scalars (device routine).
  pure function comp_Ur ( vel, sound, Ur_visc, Ur_ref, Mach_target, eps_min, Ur_factor, Ur_min ) result ( Ur )
    !$acc routine seq
    implicit none
    real(R8), intent(in) :: vel, sound, Ur_visc, Ur_ref, Mach_target, eps_min, Ur_factor, Ur_min
    real(R8) :: Ur
    ! Local
    real(R8) :: Mt, m, A_q, B_q, S_q, disc, w

    if (Ur_ref > 0.0_R8) then
      !-- Mode 1: Uniform reference velocity
      if (vel > Ur_ref) then
        Ur = sound
      else
        Ur = Ur_ref
      endif

    elseif (Mach_target > 0.0_R8) then
      !-- Mode 2: Target preconditioned Mach number (exact): A w^2 + B w + A = 0, w = (Ur/c)^2, physical root w < 1
      Mt = Mach_target
      m  = vel / sound
      A_q = 0.25_R8 * m*m * (Mt*Mt - 1.0_R8)
      B_q = Mt*Mt - 0.5_R8 * m*m * (Mt*Mt + 1.0_R8)

      if (abs(A_q) > 1.0e-30_R8) then
        S_q  = -B_q / A_q
        disc = S_q*S_q - 4.0_R8
        if (disc > 0.0_R8 .and. S_q > 2.0_R8) then
          w = 2.0_R8 / (S_q + sqrt(disc))   ! numerically stable form: w = 2/(S + sqrt(S^2-4))
          Ur = sound * sqrt(w)
        else
          Ur = sound                        ! target not achievable (M >= Mt): no preconditioning
        endif
      else
        Ur = sound * eps_min                ! m ~ 0 (stagnation): the quadratic is degenerate
      endif

      Ur = max(Ur, eps_min * sound)    ! stagnation floor
      Ur = min(Ur, sound)              ! supersonic cap
      Ur = max(Ur, Ur_visc)            ! viscous floor

    else
      !-- Mode 3: Standard Weiss-Smith
      if (vel < eps_min * sound) then
        Ur = sound * eps_min
      elseif (vel >= sound) then
        Ur = sound
      else
        Ur = vel * Ur_factor
      end if
      Ur = max(Ur, Ur_visc)            ! viscous floor

    endif

    !-- Global user-specified floor
    if (Ur_min > 0.0_R8) then
      Ur = max(Ur, Ur_min)
    endif

  end function comp_Ur


  !> The flat Ur, all rows 0 (the ghost rows are never written, as MOSE_open blk%Ur).
  subroutine Prec_Alloc ()
    use MOSE_Field_m, only: ncell_total
    implicit none
    if ( allocated(ur_f) ) return
    allocate ( ur_f(ncell_total), urb_f(ncell_total), g0_f(1,9) )
    ur_f = 0d0; urb_f = 0d0; g0_f = 0d0
    !$acc enter data copyin(ur_f, urb_f, g0_f)
  end subroutine Prec_Alloc


  !> MOSE_open Update_Derived_Variables: Ur of the interior cells of every local block, then the smoothing passes.
  !> `first`: the call before the first time step (MOSE_open Wrap_Setup), done once.
  subroutine Update_Derived_Variables ( domain, first )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_IO_RANS, only: obj_rans
    use MOSE_Lib_SpalartShur, only: gv_f
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    logical, intent(in) :: first
    ! Local
    integer :: b

    if ( first .and. ur_ready ) return
    call Prec_Alloc ()
    do b = 1, domain % nb
      if ( .not. is_local_block(b) ) cycle
      ! the velocity gradient of the eddy viscosity: the one of the last -RC source, zero before it and for the other
      ! models (MOSE_open blk%vel_gradient: zero-initialised, written by Compute_Velocity_Gradient only)
      if ( obj_rans % SpalartShur .and. allocated(gv_f) ) then
        call Ur_Blk ( b, domain % blk(b) % dl_f, domain % blk(b) % yn, domain % blk(b) % dim, gv_f, size(gv_f,1), .true. )
      else
        call Ur_Blk ( b, domain % blk(b) % dl_f, domain % blk(b) % yn, domain % blk(b) % dim, g0_f, 1, .false. )
      end if
      if ( obj_prec % n_smooth_Ur > 0 ) call Smooth_Ur_Blk ( b, domain % blk(b) % dim, obj_prec % n_smooth_Ur )
    end do
    ur_ready = .true.

  end subroutine Update_Derived_Variables


  subroutine Ur_Blk ( b, dl_f, yn, n, gsrc, ng, rcgrad )
    use MOSE_Global_m, only: gc, nprim, nsc, nu, nv, nw, np, nt, nrans, ndir, VISC, RANS
    use MOSE_Field_m, only: P_flat, field_block_index_params
    use MOSE_Lib_RANS, only: eddy_dev, prt_dev
    use MOSE_Mod_SST, only: SST_Eddy_Viscosity
    use MOSE_Mod_Spalart, only: Spalart_Eddy_Viscosity
    use MOSE_Mod_Wilcox2006, only: Wilcox2006_Eddy_Viscosity
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, f_cp_dev, f_ss_dev, co_k_mi_lam_Wilke_dev
    implicit none
    integer, intent(in) :: b, n(3), ng
    real(R8), intent(in) :: gsrc(ng, 9)   !< gv_f of Lib_SpalartShur (rcgrad) or the one zero row g0_f
    logical, intent(in) :: rcgrad
    real(R8), dimension(3, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: dl_f
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: yn
    ! Local
    integer  :: boff, pi, pj, i, j, k, s, a, c, T_i, Tint(2)
    real(R8) :: pl(nprim), rkw(2), Gvel(3,3), rho, R, T, Tdiff, vel_mag, Ur_visc, dx_min, cp, mu, lambda, mut, mu_eff, k_eff
    real(R8) :: sound, ur_ref_, mach_, eps_, fac_, urmin_
    call field_block_index_params(b, boff, pi, pj)
    ur_ref_ = obj_prec%Ur_ref; mach_ = obj_prec%Mach_target; eps_ = obj_prec%eps_min
    fac_ = obj_prec%Ur_factor; urmin_ = obj_prec%Ur_min
    !$acc data copyin(dl_f, yn)
    !$acc parallel loop collapse(3) gang vector present(P_flat, ur_f, dl_f, yn, gsrc) &
    !$acc&  private(pl, rkw, Gvel, Tint)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)
      do s = 1, nprim
        pl(s) = P_flat(WR(i,j,k), s)
      end do
      call co_rotot_Rtot_dev ( pl(1:nsc), rho, R )
      T = pl(np)/(rho*R)
      vel_mag = norm2( pl(nu:nw) )

      Ur_visc = 0.0_R8
      if ( VISC ) then
        dx_min = dl_f(1,i,j,k)
        do s = 2, ndir
          dx_min = min( dx_min, dl_f(s,i,j,k) )
        end do
        T_i = int(T); Tdiff = T - T_i; Tint(1) = T_i; Tint(2) = T_i + 1
        cp = f_cp_dev ( pl(1:nsc), Tint, Tdiff, rho )
        call co_k_mi_lam_Wilke_dev ( pl(1:nsc), rho, Tint, Tdiff, mu, lambda )
        mu_eff = mu
        k_eff  = lambda
        if ( RANS ) then
          if ( rcgrad ) then
            do c = 1, 3
            do a = 1, 3
              Gvel(a,c) = gsrc(WR(i,j,k), a+3*(c-1))
            end do
            end do
          else
            Gvel = 0d0
          end if
          if ( eddy_dev == 2 ) then
            call Spalart_Eddy_Viscosity ( mut, pl(nt:nt+nrans-1), mu, rho, Gvel, yn(i,j,k) )
          else
            rkw(1) = pl(nt); rkw(2) = pl(min(nt+1,nprim))   ! min(): no constant out-of-range subscript in an nrans=1 binary
            if ( eddy_dev == 3 ) then
              call Wilcox2006_Eddy_Viscosity ( mut, rkw, mu, rho, Gvel, yn(i,j,k) )
            else
              call SST_Eddy_Viscosity ( mut, rkw, mu, rho, Gvel, yn(i,j,k) )
            end if
          end if
          mu_eff = mu_eff + mut
          k_eff  = k_eff  + mut * cp / prt_dev
        end if
        Ur_visc = max( mu_eff / (rho * dx_min), k_eff / (rho * cp * dx_min) )
      end if

      sound = f_ss_dev ( pl(1:nsc), pl(np), rho, R )
      ur_f(WR(i,j,k)) = comp_Ur ( vel_mag, sound, Ur_visc, ur_ref_, mach_, eps_, fac_, urmin_ )
    end do
    end do
    end do
    !$acc end data
  end subroutine Ur_Blk


  !> MOSE_open Smooth_Ur_Block: n passes of max(self, the six face neighbours) on the interior, double-buffered; the
  !> ghost neighbours are the zeros of ur_f.
  subroutine Smooth_Ur_Blk ( b, n, n_smooth )
    use MOSE_Global_m, only: gc
    use MOSE_Field_m, only: field_block_index_params
    implicit none
    integer, intent(in) :: b, n(3), n_smooth
    integer :: boff, pi, pj, i, j, k, is
    call field_block_index_params(b, boff, pi, pj)
    do is = 1, n_smooth
      !$acc parallel loop collapse(3) gang vector present(ur_f, urb_f)
      do k = 1, n(3)
      do j = 1, n(2)
      do i = 1, n(1)
        urb_f(WR(i,j,k)) = max( ur_f(WR(i,j,k)), ur_f(WR(i-1,j,k)), ur_f(WR(i+1,j,k)), ur_f(WR(i,j-1,k)), &
                                ur_f(WR(i,j+1,k)), ur_f(WR(i,j,k-1)), ur_f(WR(i,j,k+1)) )
      end do
      end do
      end do
      !$acc parallel loop collapse(3) gang vector present(ur_f, urb_f)
      do k = 1, n(3)
      do j = 1, n(2)
      do i = 1, n(1)
        ur_f(WR(i,j,k)) = urb_f(WR(i,j,k))
      end do
      end do
      end do
    end do
  end subroutine Smooth_Ur_Blk

end module MOSE_Lib_Preconditioning
