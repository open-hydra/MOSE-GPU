module coronetti_mod
  implicit none
contains
  ! Coronetti for butadiene combustion (DOI:10.2514/1.B34760)
  ! Order of species: O2, C4H6, H2O, CO, CO2, H2, O, H, OH
  subroutine Coronetti(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    real(8), intent(inout) :: roi(ns)
    real(8), intent(in) :: temp 
    real(8), intent(out) :: omegadot(ns)
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(ns), Tdiff
    real(8) :: prod1, prod2, prod3, prod4, prod5, prod6
    real(8) :: kf, kb, CP1, CP2, xi, lin1, lin2
    real(8), parameter :: limitH2 = 1d-10, limitO2 = 1d-10, sigma = 23d0, tau = 17d0

    do is = 1, ns
      roi(is) = max(roi(is), 0.d0)
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
    enddo
    
    ! Preliminary
    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! irreversible reaction: C4H6 + 2 O2 -> 4 CO + 3 H2 
    kf = comp_ch_tabT(1,kf_tab,Tint,Tdiff)
    CP1 = (coi(2)**0.5)*(coi(1)**1.25)
    if (( coi(1) < limitO2 ).or.(coi(2)<limitO2) )then  
      prod1 = 0.d0
    else
      prod1 = kf*CP1
    end if

    ! irreversible reaction: C4H6 + 4 H2O -> 4 CO + 7 H2
    prod2 = comp_ch_tabT(2,kf_tab,Tint,Tdiff)*coi(2)*coi(3)
      
    ! reversible reaction: CO + H2O <--> CO2 + H2
    kf = comp_ch_tabT(3,kf_tab,Tint,Tdiff)
    kb = comp_ch_tabT(3,kb_tab,Tint,Tdiff)
    CP1 = coi(4)*coi(3)
    CP2 = coi(5)*coi(6)
    prod3 = kf*CP1 - kb*CP2

    ! reversible reaction: H2 + 1/2 O2 <--> H2O
    kf = comp_ch_tabT(4,kf_tab,Tint,Tdiff)
    kb = comp_ch_tabT(4,kb_tab,Tint,Tdiff)
    CP1 = (coi(6)**0.25)*(coi(1)**1.50)
    CP2 = (coi(3))*(coi(1))*(coi(6)**(-0.75))
    if (( coi(6) < limitH2 ).or.(coi(1)< limitH2)) then
      CP1 = 0.d0 
    end if
    if ((coi(3)<limitH2).or.(coi(1)<limitH2).or.(coi(6)<limitH2)) then
      CP2 = 0.d0 
    endif
    prod4 = kf*CP1 - kb*CP2
   
    ! reversible reaction: O2 <--> 2 O
    prod5 = comp_ch_tabT(5,kf_tab,Tint,Tdiff)*coi(1) &      
          - comp_ch_tabT(5,kb_tab,Tint,Tdiff)*coi(7)**2
      
    ! reversible reaction: H2O <--> OH + H
    prod6 = comp_ch_tabT(6,kf_tab,Tint,Tdiff)*coi(3) &
          - comp_ch_tabT(6,kb_tab,Tint,Tdiff)*coi(9)*coi(8)
      
    omegadot(1)=Wm_tab(1)*(-2.d0*prod1-0.5d0*prod4-prod5)         ! O2 
    omegadot(2)=Wm_tab(2)*(-prod1-prod2)                          ! C4H6
    omegadot(3)=Wm_tab(3)*(-4.d0*prod2-prod3+prod4-prod6)         ! H2O 
    omegadot(4)=Wm_tab(4)*(4.d0*prod1+4.d0*prod2-prod3)           ! CO
    omegadot(5)=Wm_tab(5)*(prod3)                                 ! CO2
    omegadot(6)=Wm_tab(6)*(3.d0*prod1+7.d0*prod2+prod3-prod4)     ! H2
    omegadot(7)=Wm_tab(7)*(2.d0*prod5)                            ! O
    omegadot(8)=Wm_tab(8)*(prod6)                                 ! H 
    omegadot(9)=Wm_tab(9)*(prod6)                                 ! OH

  end subroutine Coronetti
end module coronetti_mod
