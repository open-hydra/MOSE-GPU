! Source: https://cantera.org/3.1/reference/kinetics/rate-constants.html#sec-falloff-rate
module FLINT_Lib_Chemistry_falloff
  use FLINT_Lib_Chemistry_data
  implicit none
  private
  public :: f_k_troe
  public :: f_k_lindemann

contains

  ! Computes the Lindemann rates for a reaction given temperature, third-body concentration, and reaction index
  pure function f_k_lindemann(ireact, Tint, Tdiff, coM) result(rate)
    real(8), intent(in) :: coM
    integer, intent(in) :: ireact
    integer, intent(in) :: Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: rate(2)
    real(8) :: kinf, k0, kc, Fcent, Pr

    ! Compute the Troe rate using the precomputed tables
    kinf = f_kinf_lind(ireact, Tint, Tdiff)
    k0 = f_k0_lind(ireact, Tint, Tdiff)
    kc = f_kc_lind(ireact, Tint, Tdiff)

    ! Reduced pressure
    Pr = k0 * coM / kinf

    ! Forward rate
    rate(1) = kinf * ( Pr / (1d0+Pr) )
    ! Backward rate
    rate(2) = rate(1) / kc

  end function f_k_lindemann

  ! Computes the Troe rates for a reaction given temperature, third-body concentration, and reaction index
  pure function f_k_troe(ireact, Tint, Tdiff, coM) result(rate)
    real(8), intent(in) :: coM
    integer, intent(in) :: ireact
    integer, intent(in) :: Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: rate(2)
    real(8) :: kinf, k0, kc, Fcent, Pr

    ! Compute the Troe rate using the precomputed tables
    kinf = f_kinf_troe(ireact, Tint, Tdiff)
    k0 = f_k0_troe(ireact, Tint, Tdiff)
    kc = f_kc_troe(ireact, Tint, Tdiff)
    Fcent = f_Fcent(ireact, Tint, Tdiff)

    ! Reduced pressure
    Pr = k0 * coM / kinf

    ! Forward rate
    rate(1) = kinf * ( Pr / (1d0+Pr) ) * f_F(Pr, Fcent)
    ! Backward rate
    rate(2) = rate(1) / kc

  end function f_k_troe

  ! Computes the Troe falloff correction factor F
  pure function f_F(Pr, Fcent) result(F)
    implicit none

    real(8), intent(in)  :: Pr      ! Reduced pressure
    real(8), intent(in)  :: Fcent   ! Known center factor
    real(8)              :: F       ! Final Troe correction

    real(8) :: logFcent, logPr, c, n, f1

    logFcent = log10(Fcent)
    logPr = log10(Pr)

    c = -0.4d0 - 0.67d0 * logFcent
    n = 0.75d0 - 1.27d0 * logFcent

    f1 = (logPr + c) / (n - 0.14d0* (logPr + c))

    F = 10d0**(logFcent / (1d0 + f1*f1))

  end function f_F


end module FLINT_Lib_Chemistry_falloff