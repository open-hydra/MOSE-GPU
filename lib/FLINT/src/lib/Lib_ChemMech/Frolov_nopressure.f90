!>@brief Single-step H2/O2 combustion without pressure correction.
!>
!> Mechanism (from Frolov_nopressure.yaml):
!>     2 H2 + O2  <=>  2 H2O
!> Forward Arrhenius:  A = 8e11,  b = 0,  Ea = 8.314e7 J/kmol
!> Reverse rate computed from equilibrium constant by the chemistry loader
!> and supplied through kb_tab (same convention as ONERA-7, general, etc.).
!>
!> Species ordering (matches the YAML phase definition `species: [O2,H2O,H2,N2]`):
!>     1: O2
!>     2: H2O
!>     3: H2
!>     4: N2     (inert; row is all zero)
!>
!> Hand-rolled to bypass the data-driven `general` dispatcher in
!> FLINT_Lib_Chemistry_wdot. On the THOR Scalability case the `general`
!> dispatcher was the #2 hotspot at ~9000 s of CPU per rank — dominated by
!> dispatch overhead (outer reaction loop, branch on `ni1 /= 0`, integer
!> exponent via `**nint(...)`). This straight-line version eliminates all of
!> that for the single reaction.
module Frolov_nopressure_mod
  implicit none
contains


  !--------------------------------------------------------------------------
  ! RHS: species mass production rate (kg/m^3/s)
  !--------------------------------------------------------------------------
  subroutine Frolov_nopressure(roi, temp, omegadot)
    !$acc routine seq
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    real(8), intent(inout) :: roi(ns)
    real(8), intent(in)    :: temp
    real(8), intent(out)   :: omegadot(ns)
    ! Local (fixed-size: device builds avoid per-call heap malloc; 1:ns used)
    real(8) :: coi(NSCHEM), Tdiff
    integer :: is, T_i, Tint(2)
    real(8) :: kf, kb, net_rate

    do is = 1, ns
      coi(is) = roi(is) / Wm_tab(is)   ! kmol/m^3
    end do

    T_i     = int(temp)
    Tdiff   = temp - T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    kf = f_kf(1, Tint, Tdiff)
    kb = f_kb(1, Tint, Tdiff)

    ! Reaction 1:  2 H2 + O2 <=> 2 H2O
    ! prodf = kf * c(3)^2 * c(1) ;  prodb = kb * c(2)^2
    net_rate = kf * coi(3) * coi(3) * coi(1) - kb * coi(2) * coi(2)

    ! Stoichiometric coefficients (signed, product − reactant):
    !   O2 (1): -1     H2O (2): +2     H2 (3): -2     N2 (4): 0
    omegadot(1) = -      Wm_tab(1) * net_rate
    omegadot(2) =  2.d0 * Wm_tab(2) * net_rate
    omegadot(3) = -2.d0 * Wm_tab(3) * net_rate
    omegadot(4) = 0.d0
  end subroutine Frolov_nopressure


  !--------------------------------------------------------------------------
  ! Analytical Jacobian.
  !
  !   dwdr(i,j) = d omegadot(i) / d roi(j)     [1/s]
  !   dwdT(i)   = d omegadot(i) / d T          [kg/(m^3 s K)]
  !
  ! Same table-derivative trick as ONERA_7_jac: kf and kb are linearly
  ! interpolated in T from kf_tab / kb_tab, so dk/dT is a free table lookup.
  !--------------------------------------------------------------------------
  subroutine Frolov_nopressure_jac(roi, temp, dwdr, dwdT)
    !$acc routine seq
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    real(8), intent(in)  :: roi(ns), temp
    ! dwdr: fixed leading dimension NSCHEM to match jac_native/chemjac_if storage
    ! (sequence association). Only (1:ns,1:ns) written.
    real(8), intent(out) :: dwdr(NSCHEM, ns)
    real(8), intent(out) :: dwdT(ns)
    ! Local (fixed-size NSCHEM: no device-heap malloc; sliced 1:ns throughout)
    real(8) :: coi(NSCHEM), Tdiff
    integer :: is, T_i, j
    real(8) :: kf, kb, dkf_dT, dkb_dT
    real(8) :: dnet_dc(NSCHEM), dnet_dT_r
    real(8) :: dwdr_c(NSCHEM, NSCHEM)

    do is = 1, ns
      coi(is) = roi(is) / Wm_tab(is)
    end do

    T_i   = int(temp)
    Tdiff = temp - T_i

    kf     = kf_tab(T_i  , 1) + Tdiff * (kf_tab(T_i+1, 1) - kf_tab(T_i, 1))
    kb     = kb_tab(T_i  , 1) + Tdiff * (kb_tab(T_i+1, 1) - kb_tab(T_i, 1))
    dkf_dT = kf_tab(T_i+1, 1) - kf_tab(T_i, 1)
    dkb_dT = kb_tab(T_i+1, 1) - kb_tab(T_i, 1)

    ! ∂(net_rate)/∂coi(j) for reaction 1:
    !   net_rate = kf * c(3)^2 * c(1) - kb * c(2)^2
    !   ∂/∂c(1) =        kf * c(3)^2
    !   ∂/∂c(2) = - 2 * kb * c(2)
    !   ∂/∂c(3) =   2 * kf * c(3) * c(1)
    !   ∂/∂c(4) = 0  (N2 inert)
    dnet_dc(1:ns) = 0.d0
    dnet_dc(1) =          kf * coi(3) * coi(3)
    dnet_dc(2) = -2.d0 *  kb * coi(2)
    dnet_dc(3) =  2.d0 *  kf * coi(3) * coi(1)
    dnet_dT_r  = dkf_dT * coi(3) * coi(3) * coi(1) &
               - dkb_dT * coi(2) * coi(2)

    ! Apply stoichiometry to get ∂omegadot/∂coi(j) and ∂omegadot/∂T.
    dwdr_c(1:ns,1:ns) = 0.d0
    dwdr_c(1, 1:ns) = -      Wm_tab(1) * dnet_dc(1:ns)
    dwdr_c(2, 1:ns) =  2.d0 * Wm_tab(2) * dnet_dc(1:ns)
    dwdr_c(3, 1:ns) = -2.d0 * Wm_tab(3) * dnet_dc(1:ns)
    ! row 4 (N2) stays zero by construction.

    dwdT(:) = 0.d0
    dwdT(1) = -      Wm_tab(1) * dnet_dT_r
    dwdT(2) =  2.d0 * Wm_tab(2) * dnet_dT_r
    dwdT(3) = -2.d0 * Wm_tab(3) * dnet_dT_r

    ! Convert ∂/∂coi(j) → ∂/∂roi(j): divide column j by Wm_tab(j).
    do j = 1, ns
      dwdr(1:ns, j) = dwdr_c(1:ns, j) / Wm_tab(j)
    end do
  end subroutine Frolov_nopressure_jac


end module Frolov_nopressure_mod
