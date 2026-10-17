module FLINT_Lib_Chemistry_data
  implicit none

  !> Fixed upper bounds for the device chemistry work arrays. The GPU build
  !> cannot use runtime-sized (ns-dependent) automatics — each would be a
  !> per-thread per-call device-heap malloc through a serializing allocator.
  !> Sized for the mechanisms allowed on the device (mech_dev guard):
  !> ONERA-7 (ns=7) is the largest. Bump if a bigger device mechanism is added.
  !> DECOUPLED device-path array sizing. The species count is a thermo/composition
  !> property INDEPENDENT of the reacting mechanism: a FROZEN run carries species but
  !> no mechanism and may have MORE species than any device mechanism.
  !>   NSMAX  = max species on the device THERMO/transport path (frozen-capable; bound
  !>            by the MOSE device cap NPRIM_DEV). Sizes coi/cpi/Wilke-fi. Bump for a
  !>            frozen case with more species (independent of any mechanism).
  !>   NSCHEM = max species of a device-supported reacting MECHANISM (ONERA-7=7 largest).
  !>            Sizes the chemistry rhs/jac + the radau5 ODE solver ONLY (reactive runs);
  !>            kept tight for per-thread memory. Bump only when a bigger mech is added.
  integer, parameter, public :: NSMAX  = 8
  integer, parameter, public :: NSCHEM = 7
  integer, parameter, public :: NZMAX  = NSCHEM+1  !> max ODE system size = mech species + 1

  integer                              :: nrc
  integer, dimension(:), allocatable   :: rxn_type ! 0 -> Arrhenius, 1 -> Troe, 2 -> Lindemann
  ! Arrhenius
  integer                              :: nrc_arrh
  real(8), dimension(:,:), allocatable :: kf_tab
  real(8), dimension(:,:), allocatable :: kb_tab
  ! Device copies for the OpenACC chemistry port (uploaded once by
  ! flint_acc_upload_tables; referenced from acc-routine device code).
  !$acc declare create(kf_tab, kb_tab)
  ! Falloff-Troe
  integer                              :: nrc_troe
  real(8), dimension(:,:), allocatable :: kinf_troe_tab
  real(8), dimension(:,:), allocatable :: k0_troe_tab
  real(8), dimension(:,:), allocatable :: kc_troe_tab
  real(8), dimension(:,:), allocatable :: Fcent_tab
  ! Falloff-Lindemann
  integer                              :: nrc_lindemann
  real(8), dimension(:,:), allocatable :: kinf_lind_tab
  real(8), dimension(:,:), allocatable :: k0_lind_tab
  real(8), dimension(:,:), allocatable :: kc_lind_tab
  ! Arrhenius (for the general loop only)
  real(8), dimension(:,:), allocatable :: ni1_arrh_tab
  real(8), dimension(:,:), allocatable :: ni2_arrh_tab
  real(8), dimension(:,:), allocatable :: epsch_arrh_tab
  ! Falloff-Troe (for the general loop only)
  real(8), dimension(:,:), allocatable :: ni1_troe_tab
  real(8), dimension(:,:), allocatable :: ni2_troe_tab
  real(8), dimension(:,:), allocatable :: epsch_troe_tab
  ! Falloff-Lindemann (for the general loop only)
  real(8), dimension(:,:), allocatable :: ni1_lind_tab
  real(8), dimension(:,:), allocatable :: ni2_lind_tab
  real(8), dimension(:,:), allocatable :: epsch_lind_tab

contains

  pure function comp_ch_tabT(ireact,tab,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: tab(:,:), Tdiff
    ! Local
    real(8) :: a, b
    real(8) :: result
      
    a = tab(Tint(1),ireact)      ! int(T)   <- Tint(1)
    b = tab(Tint(2),ireact)      ! int(T)+1 <- Tint(2)
    result = a+(b-a)*Tdiff

  end function comp_ch_tabT

  pure function f_kf(ireact,Tint,Tdiff) result(result)
    !$acc routine seq
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = kf_tab(Tint(1),ireact)
    b = kf_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_kf

  pure function f_kb(ireact,Tint,Tdiff) result(result)
    !$acc routine seq
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = kb_tab(Tint(1),ireact)
    b = kb_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_kb

  pure function f_kc_lind(ireact,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = kc_lind_tab(Tint(1),ireact)
    b = kc_lind_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_kc_lind

  pure function f_kinf_lind(ireact,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = kinf_lind_tab(Tint(1),ireact)
    b = kinf_lind_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_kinf_lind

  pure function f_k0_lind(ireact,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = k0_lind_tab(Tint(1),ireact)
    b = k0_lind_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_k0_lind

  pure function f_kc_troe(ireact,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = kc_troe_tab(Tint(1),ireact)
    b = kc_troe_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_kc_troe

  pure function f_kinf_troe(ireact,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = kinf_troe_tab(Tint(1),ireact)
    b = kinf_troe_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_kinf_troe

  pure function f_k0_troe(ireact,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = k0_troe_tab(Tint(1),ireact)
    b = k0_troe_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_k0_troe

  pure function f_Fcent(ireact,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = Fcent_tab(Tint(1),ireact)
    b = Fcent_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_Fcent

  ! Free all chemistry data
  subroutine free_chemistry_data()
    implicit none
    if (allocated(rxn_type)) deallocate(rxn_type)
    if (allocated(kf_tab)) deallocate(kf_tab)
    if (allocated(kb_tab)) deallocate(kb_tab)
    if (allocated(kinf_lind_tab)) deallocate(kinf_lind_tab)
    if (allocated(k0_lind_tab)) deallocate(k0_lind_tab)
    if (allocated(kc_lind_tab)) deallocate(kc_lind_tab)
    if (allocated(kinf_troe_tab)) deallocate(kinf_troe_tab)
    if (allocated(k0_troe_tab)) deallocate(k0_troe_tab)
    if (allocated(kc_troe_tab)) deallocate(kc_troe_tab)
    if (allocated(Fcent_tab)) deallocate(Fcent_tab)
    if (allocated(ni1_arrh_tab)) deallocate(ni1_arrh_tab)
    if (allocated(ni2_arrh_tab)) deallocate(ni2_arrh_tab)
    if (allocated(epsch_arrh_tab)) deallocate(epsch_arrh_tab)
    if (allocated(ni1_lind_tab)) deallocate(ni1_lind_tab)
    if (allocated(ni2_lind_tab)) deallocate(ni2_lind_tab)
    if (allocated(epsch_lind_tab)) deallocate(epsch_lind_tab)
    if (allocated(ni1_troe_tab)) deallocate(ni1_troe_tab)
    if (allocated(ni2_troe_tab)) deallocate(ni2_troe_tab)
    if (allocated(epsch_troe_tab)) deallocate(epsch_troe_tab)
  end subroutine free_chemistry_data

end module FLINT_Lib_Chemistry_data