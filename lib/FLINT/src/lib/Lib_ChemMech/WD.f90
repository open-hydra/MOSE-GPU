  ! WD: Global Westbrook-Dryer mechanism
  ! 5 species & 3 reactions
  module WD_mod
    implicit none
    contains
  subroutine WD(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    real(8), intent(inout) :: roi(ns)
    real(8), intent(in)    :: temp 
    real(8), intent(out)   :: omegadot(ns)
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(ns), Tdiff
    real(8) :: prod1,prod2,prod3

    do is = 1, ns
      roi(is) = max(roi(is), 0.d0)
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
    enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! species: [CH4, O2, CO2, H2O, CO]

    ! CH4 + 1.5 O2 => CO + 2 H2O
    prod1 = comp_ch_tabT(1,kf_tab,Tint,Tdiff)*(coi(1)**0.70)*(coi(2)**0.80)

    ! CO + 0.5 O2 + H2O => CO2 + H2O
    prod2 = comp_ch_tabT(2,kf_tab,Tint,Tdiff)*(coi(4)**0.5)*coi(5)*(coi(2)**0.25)

    ! CO2 => CO + 0.5 O2
    prod3 = comp_ch_tabT(3,kf_tab,Tint,Tdiff)*coi(3)
     
    ! Chemical Source Terms
    omegadot = 0d0
    omegadot(1)=Wm_tab(1)*(-prod1)
    omegadot(2)=Wm_tab(2)*(-1.5*prod1-0.5*prod2+0.5*prod3)
    omegadot(3)=Wm_tab(3)*(prod2-prod3)
    omegadot(4)=Wm_tab(4)*(2*prod1)
    omegadot(5)=Wm_tab(5)*(prod1-prod2+prod3)

  end subroutine WD

  subroutine Andersen(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    real(8), intent(inout) :: roi(ns)
    real(8), intent(in)    :: temp 
    real(8), intent(out)   :: omegadot(ns)
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(ns), Tdiff
    real(8) :: prod1,prod2,prod3,prod4,prod5,prod6

    do is = 1, ns
      roi(is) = max(roi(is), 0.d0)
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
    enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1
    ! species: [CH4, O2, CO2, H2O, CO]

    ! CH4 + 1.5 O2 => CO + 2 H2O
    prod1 = comp_ch_tabT(1,kf_tab,Tint,Tdiff)*(coi(1)**0.70)*(coi(2)**0.80)

    ! CO + 0.5 O2 + H2O => CO2 + H2O
    prod2 = comp_ch_tabT(2,kf_tab,Tint,Tdiff)*(coi(4)**0.5)*coi(5)*(coi(2)**0.25)

    ! CO2 => CO + 0.5 O2
    prod3 = comp_ch_tabT(3,kf_tab,Tint,Tdiff)*(coi(3)**1.25)
     
    ! Chemical Source Terms
    omegadot = 0d0
    omegadot(1)=Wm_tab(1)*(-prod1)
    omegadot(2)=Wm_tab(2)*(-1.5*prod1-0.5*prod2+0.5*prod3)
    omegadot(3)=Wm_tab(3)*(prod2-prod3)
    omegadot(4)=Wm_tab(4)*(2*prod1)
    omegadot(5)=Wm_tab(5)*(prod1-prod2+prod3)

  end subroutine Andersen

  subroutine OSK(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    real(8), intent(inout) :: roi(ns)
    real(8), intent(in)    :: temp 
    real(8), intent(out)   :: omegadot(ns)
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(ns), Tdiff
    real(8) :: prod1,prod2,prod3,prod4,prod5,prod6

    do is = 1, ns
      roi(is) = max(roi(is), 0.d0)
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
    enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1
    ! species: [O2, CH4, H2O, CO2]

    ! CH4 + 2 O2 => CO2 + 2 H2O
    prod1 = comp_ch_tabT(1,kf_tab,Tint,Tdiff)*(coi(2)**0.7)*(coi(1)**0.8)

     
    ! Chemical Source Terms
    omegadot = 0d0
    omegadot(1)=Wm_tab(1)*(-2*prod1)
    omegadot(2)=Wm_tab(2)*(-prod1)
    omegadot(3)=Wm_tab(3)*(2*prod1)
    omegadot(4)=Wm_tab(4)*(prod1)

    !omegadot(5)=Wm_tab(5)*(prod1-prod2+prod3)

  end subroutine OSK
  


end module WD_mod
