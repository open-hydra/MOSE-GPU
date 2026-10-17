
module FLINT_Lib_Chemistry_wdot
  implicit none

  !> Concrete procedure pointing to one of the subroutine realizations
  procedure(chemsource_if), pointer, public :: chemistry_source

  !> Concrete procedure pointing to the analytical Jacobian of `chemistry_source`
  !> (null when the active mechanism has no analytical Jacobian yet; callers
  !> must fall back to finite-difference Jacobian in that case).
  procedure(chemjac_if), pointer, public :: chemistry_jacobian => null()

  !> Device-build mechanism id: procedure pointers cannot be called on the
  !> GPU, so rhs_native/jac_native dispatch on this id with direct calls.
  !> 1 = ONERA-7, 2 = Frolov_nopressure, 0 = unsupported on device.
  integer, public :: mech_dev = 0
  !$acc declare create(mech_dev)

  !> Abstract interface relative to the finite-rate reactions source procedure
  abstract interface
  subroutine chemsource_if(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    integer :: is, T_i, Tint(2)
    real(8), intent(inout) :: roi(ns)
    real(8), intent(in) :: temp
    real(8), intent(out) :: omegadot(ns)
    real(8) :: coi(ns+1), Tdiff
  end subroutine chemsource_if
  end interface

  !> Abstract interface for the analytical chemistry Jacobian.
  !>   dwdr(i,j) = d omegadot(i) / d roi(j)
  !>   dwdT(i)   = d omegadot(i) / d T
  abstract interface
  subroutine chemjac_if(roi,temp,dwdr,dwdT)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data, only: NSCHEM
    implicit none
    real(8), intent(in)  :: roi(ns), temp
    ! dwdr carries a fixed leading dimension NSCHEM (not ns) so the column
    ! stride agrees with the caller's fixed-size array (jac_native): with
    ! sequence association the leading extent must match on both sides, else
    ! columns misalign. Only the (1:ns,1:ns) block is written/read.
    real(8), intent(out) :: dwdr(NSCHEM,ns)
    real(8), intent(out) :: dwdT(ns)
  end subroutine chemjac_if
  end interface

contains

  subroutine Assign_Mechanism(mad_world)
    use WD_mod
    use globH2_mod
    use JLRs_mod
    use smooke_mod
    use coria_mod
    use TSRCDF13_mod
    use TSRGP24_mod
    use TSRRich31_mod
    use ZK_mod
    use coronetti_mod
    use singh_mod
    use troyes_mod
    use ecker_mod
    use cross_mod
    use pelucchi_mod
    use ONERA7_mod
    use sandiego_mod
    use FFCMy_12_mod
    use Frolov_nopressure_mod
    implicit none
    character(*), intent(in) :: mad_world

    ! Device-build dispatch id for rhs_native/jac_native (direct calls).
    select case(mad_world)
    case('ONERA-7')
      mech_dev = 1
    case('Frolov_nopressure')
      mech_dev = 2
    case default
      mech_dev = 0
    end select
# if defined (_OPENACC)
    ! Device build: rhs_native/jac_native are devirtualized to direct calls
    ! (procedure pointers cannot run on the GPU). Any unsupported mechanism
    ! would silently integrate the wrong kinetics -> hard stop.
    if (mech_dev == 0) then
      write(*,*) '[OPENACC] device build supports only ONERA-7 / Frolov_nopressure, got: '//trim(mad_world)
      stop 1
    end if
    !$acc update device(mech_dev)
# endif

    ! Default: no analytical Jacobian available. Each mechanism that has one
    ! overrides this below.
    chemistry_jacobian => null()

    select case(mad_world)
    case('WD')
      chemistry_source => WD
    case('JLR-Nasuti')
      chemistry_source => JLR
    case('Frassoldati')
      chemistry_source => Frassoldati
    case('Smooke')
      chemistry_source => smooke
    case('CORIA-CNRS')
      chemistry_source => coria
    case('TSR-CDF-13')
      chemistry_source => TSRCDF13
    case('TSR-GP-24')
      chemistry_source => TSRGP24
    case('TSR-Rich-31')
      chemistry_source => TSRRich31
    case('ZK')
      chemistry_source => ZK
    case('CoronettiC4H6')
      chemistry_source => Coronetti
    case('CKJLR-10sp')
      chemistry_source => CKJLR10sp
    case('Singh')
      chemistry_source => Singh
    case('Singh-WC32')
      chemistry_source => Singh_WC32
    case('Frolov')
      chemistry_source => Frolov
    case('Nassini')
      chemistry_source => Nassini_4
    case('Troyes')
      chemistry_source => troyes
    case('Ecker')
      chemistry_source => ecker
    case('Cross')
      chemistry_source => cross
    case('Pelucchi')
      chemistry_source => pelucchi
    case('WD-Andersen')
      chemistry_source => Andersen
    case('OSK')
      chemistry_source => OSK
    case('ONERA-7')
      chemistry_source   => ONERA_7
      chemistry_jacobian => ONERA_7_jac
    case('Frolov_nopressure')
      chemistry_source   => Frolov_nopressure
      chemistry_jacobian => Frolov_nopressure_jac
    case('SanDiego')
      chemistry_source => sandiego20161214
    case('FFCMy-12')
      chemistry_source => FFCMy_12

    case default
      write(*,*) "[WARNING] Explicit procedure for "//trim(mad_world)//" not found, defaulting to the general procedure"
      chemistry_source => general
    end select

  end subroutine Assign_Mechanism


  ! General mechanism
  subroutine general(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    use FLINT_Lib_Chemistry_falloff
    implicit none
    real(8), intent(inout) :: roi(ns)
    real(8), intent(in) :: temp 
    real(8), intent(out) :: omegadot(ns)
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(ns), Tdiff
    real(8) :: prod_fwd, prod_rev, deltani, k(2)
    real(8) :: rate_fwd, rate_rev, net_rate, coM
    integer :: ir

    do is = 1, ns
      coi(is) = roi(is)/Wm_tab(is)  ! kmol/m^3
    enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! Initialize mass source terms
    omegadot(:) = 0.d0

    ! Loop over Arrhenius reactions
    do ir = 1, nrc_arrh

      ! Compute third-body effective concentration
      if (ni1_arrh_tab(ns+1, ir)>0d0) then
        coM = 0.d0
        do is = 1, ns
          coM = coM + coi(is) * epsch_arrh_tab(is, ir)
        enddo
      else 
        coM = 1.d0
      endif

      ! Compute forward and reverse rate-of-progress
      prod_fwd = 1.0d0
      prod_rev = 1.0d0
      do is = 1, ns
        if (ni1_arrh_tab(is, ir) /= 0) prod_fwd = prod_fwd * coi(is)**nint(ni1_arrh_tab(is, ir))
        if (ni2_arrh_tab(is, ir) /= 0) prod_rev = prod_rev * coi(is)**nint(ni2_arrh_tab(is, ir))
      enddo
      rate_fwd = f_kf(ir,Tint,Tdiff) * prod_fwd * coM
      rate_rev = f_kb(ir,Tint,Tdiff) * prod_rev * coM
      net_rate = rate_fwd - rate_rev

      ! Sum up net production rate for each species
      do is = 1, ns
        deltani = ni2_arrh_tab(is, ir) - ni1_arrh_tab(is, ir)
        omegadot(is) = omegadot(is) + Wm_tab(is) * deltani * net_rate
      enddo

    enddo

    ! Loop over falloff-Troe reactions
    do ir = 1, nrc_troe

      ! Compute third-body effective concentration
      if (ni1_troe_tab(ns+1, ir)>0d0) then
        coM = 0.d0
        do is = 1, ns
          coM = coM + coi(is) * epsch_troe_tab(is, ir)
        enddo
      else 
        coM = 1.d0
      endif

      ! Compute forward and reverse rate-of-progress
      prod_fwd = 1.0d0
      prod_rev = 1.0d0
      do is = 1, ns
        if (ni1_troe_tab(is, ir) /= 0) prod_fwd = prod_fwd * coi(is)**ni1_troe_tab(is, ir)
        if (ni2_troe_tab(is, ir) /= 0) prod_rev = prod_rev * coi(is)**ni2_troe_tab(is, ir)
      enddo
      k = f_k_troe(ir,Tint,Tdiff,coM)
      rate_fwd = k(1) * prod_fwd
      rate_rev = k(2) * prod_rev
      net_rate = rate_fwd - rate_rev

      ! Sum up net production rate for each species
      do is = 1, ns
        deltani = ni2_troe_tab(is, ir) - ni1_troe_tab(is, ir)
        omegadot(is) = omegadot(is) + Wm_tab(is) * deltani * net_rate
      enddo

    enddo

    ! Loop over falloff-Lindemann reactions
    do ir = 1, nrc_lindemann

      ! Compute third-body effective concentration
      if (ni1_lind_tab(ns+1, ir)>0d0) then
        coM = 0.d0
        do is = 1, ns
          coM = coM + coi(is) * epsch_lind_tab(is, ir)
        enddo
      else 
        coM = 1.d0
      endif

      ! Compute forward and reverse rate-of-progress
      prod_fwd = 1.0d0
      prod_rev = 1.0d0
      do is = 1, ns
        if (ni1_lind_tab(is, ir) /= 0) prod_fwd = prod_fwd * coi(is)**ni1_lind_tab(is, ir)
        if (ni2_lind_tab(is, ir) /= 0) prod_rev = prod_rev * coi(is)**ni2_lind_tab(is, ir)
      enddo
      k = f_k_lindemann(ir,Tint,Tdiff,coM)
      rate_fwd = k(1) * prod_fwd
      rate_rev = k(2) * prod_rev
      net_rate = rate_fwd - rate_rev

      ! Sum up net production rate for each species
      do is = 1, ns
        deltani = ni2_lind_tab(is, ir) - ni1_lind_tab(is, ir)
        omegadot(is) = omegadot(is) + Wm_tab(is) * deltani * net_rate
      enddo

    enddo

  end subroutine general

end module FLINT_Lib_Chemistry_wdot
