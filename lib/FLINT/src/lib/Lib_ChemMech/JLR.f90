module JLRs_mod
  implicit none
contains
  ! John Lindstedt with Recombination Reaction Mechanism
  ! 9 species & 7 reactions
  subroutine JLR(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    real(8), intent(inout) :: roi(ns)
    real(8), intent(in)    :: temp 
    real(8), intent(out)   :: omegadot(ns)
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(ns), Tdiff
    real(8) :: prod1,prod2,prod3,prod4,prod5,prod6,prod7

   do is = 1, ns
     coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
     ! Loop done in order to avoid numerical issues in omegadot evaluation
     ! Very low coi could produce finite prods
     if (coi(is).lt.1d-10) coi(is) = 0.0d0
   enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! 0.5 CH4 + 1.25 O2 --> CO + 2 H2 - 0.5 CH4 + 0.75 O2
    prod1 = comp_ch_tabT(1,kf_tab,Tint,Tdiff)*(coi(2)**0.50)*(coi(1)**1.25)

    ! CH4 + H2O --> CO + 3 H2
    prod2 = comp_ch_tabT(2,kf_tab,Tint,Tdiff)*(coi(2)*coi(3))

    ! CO + H2O <--> CO2 + H2
    prod3 = comp_ch_tabT(3,kf_tab,Tint,Tdiff)*coi(4)*coi(3)- &
            comp_ch_tabT(3,kb_tab,Tint,Tdiff)*coi(5)*coi(6)

    ! 1/4 H2 + 3/2 O2 <--> H2O + O2 - 3/4 H2
    if (coi(6) < 1.d-10) then
      prod4 = (comp_ch_tabT(4,kf_tab,Tint,Tdiff)*(coi(6)**0.25)*(coi(1)**1.50))
    else
      prod4 = comp_ch_tabT(4,kf_tab,Tint,Tdiff)*(coi(6)**0.25)*(coi(1)**1.50) - & 
              comp_ch_tabT(4,kb_tab,Tint,Tdiff)*(coi(3))*(coi(1))*(coi(6)**(-0.75))
    endif

    ! O2 <--> 2O
    prod5 = comp_ch_tabT(5,kf_tab,Tint,Tdiff)*coi(1)-comp_ch_tabT(5,kb_tab,Tint,Tdiff)*coi(8)**2

    ! H2O <--> H + OH
    prod6 = comp_ch_tabT(6,kf_tab,Tint,Tdiff)*coi(3)-comp_ch_tabT(6,kb_tab,Tint,Tdiff)*coi(7)*coi(9)

    ! OH + H2 <--> H + H2O
    prod7 = comp_ch_tabT(7,kf_tab,Tint,Tdiff)*coi(9)*coi(6)- &
            comp_ch_tabT(7,kb_tab,Tint,Tdiff)*coi(7)*coi(3)
     
    ! Chemical Source Terms
    omegadot = 0d0
    omegadot(1)=Wm_tab(1)*(-0.5*prod1-0.5*prod4-prod5)        !O2    [kg/(m3*s)]
    omegadot(2)=Wm_tab(2)*(-prod1-prod2)                      !CH4   [kg/(m3*s)]
    omegadot(3)=Wm_tab(3)*(-prod2-prod3+prod4-prod6+prod7)    !H2O   [kg/(m3*s)]
    omegadot(4)=Wm_tab(4)*(-prod3+prod1+prod2)                !CO    [kg/(m3*s)]
    omegadot(5)=Wm_tab(5)*(prod3)                             !CO2   [kg/(m3*s)]
    omegadot(6)=Wm_tab(6)*(2*prod1+3*prod2+prod3-prod4-prod7) !H2    [kg/(m3*s)]
    omegadot(7)=Wm_tab(7)*(prod6+prod7)                       !H     [kg/(m3*s)]
    omegadot(8)=Wm_tab(8)*(2*prod5)                           !O     [kg/(m3*s)]
    omegadot(9)=Wm_tab(9)*(prod6-prod7)                       !OH    [kg/(m3*s)]

  end subroutine JLR

  ! Frassoldati: John Lindstedt with Recombination Reaction Mechanism
  ! 9 species & 6 reactions
  subroutine Frassoldati(roi,temp,omegadot)
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
     coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
     ! Loop done in order to avoid numerical issues in omegadot evaluation
     ! Very low coi could produce finite prods
     if (coi(is).lt.1d-10) coi(is) = 0.0d0
   enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! 0.5 CH4 + 1.25 O2 --> CO + 2 H2 - 0.5 CH4 + 0.75 O2
    prod1 = comp_ch_tabT(1,kf_tab,Tint,Tdiff)*(coi(2)**0.50)*(coi(1)**1.30)

    ! CH4 + H2O --> CO + 3 H2
    prod2 = comp_ch_tabT(2,kf_tab,Tint,Tdiff)*(coi(2)*coi(3))

    ! CO + H2O <--> CO2 + H2
    prod3 = comp_ch_tabT(3,kf_tab,Tint,Tdiff)*coi(4)*coi(3)- &
            comp_ch_tabT(3,kb_tab,Tint,Tdiff)*coi(5)*coi(6)

    ! 1/4 H2 + 3/2 O2 <--> H2O + O2 - 3/4 H2
    if (coi(6) < 1.d-10) then
      prod4 = (comp_ch_tabT(4,kf_tab,Tint,Tdiff)*(coi(6)**0.3)*(coi(1)**1.55))
    else
      prod4 = comp_ch_tabT(4,kf_tab,Tint,Tdiff)*(coi(6)**0.3)*(coi(1)**1.55) - & 
              comp_ch_tabT(4,kb_tab,Tint,Tdiff)*(coi(3))*(coi(1))*(coi(6)**(-0.75))
    endif

    ! O2 <--> 2O
    prod5 = comp_ch_tabT(5,kf_tab,Tint,Tdiff)*coi(1)-comp_ch_tabT(5,kb_tab,Tint,Tdiff)*coi(8)**2

    ! H2O <--> H + OH
    prod6 = comp_ch_tabT(6,kf_tab,Tint,Tdiff)*coi(3)-comp_ch_tabT(6,kb_tab,Tint,Tdiff)*coi(7)*coi(9)
     
    ! Chemical Source Terms
    omegadot = 0d0
    omegadot(1)=Wm_tab(1)*(-0.5*prod1-0.5*prod4-prod5)        !O2    [kg/(m3*s)]
    omegadot(2)=Wm_tab(2)*(-prod1-prod2)                      !CH4   [kg/(m3*s)]
    omegadot(3)=Wm_tab(3)*(-prod2-prod3+prod4-prod6)          !H2O   [kg/(m3*s)]
    omegadot(4)=Wm_tab(4)*(-prod3+prod1+prod2)                !CO    [kg/(m3*s)]
    omegadot(5)=Wm_tab(5)*(prod3)                             !CO2   [kg/(m3*s)]
    omegadot(6)=Wm_tab(6)*(2*prod1+3*prod2+prod3-prod4)       !H2    [kg/(m3*s)]
    omegadot(7)=Wm_tab(7)*(prod6)                             !H     [kg/(m3*s)]
    omegadot(8)=Wm_tab(8)*(2*prod5)                           !O     [kg/(m3*s)]
    omegadot(9)=Wm_tab(9)*(prod6)                             !OH    [kg/(m3*s)]

  end subroutine Frassoldati

  !> John Lindstedt with Recombination Reaction Mechanism adapted for
  !> kerosene decomposition products combustion 10 species & 8 reactions
  !> SPECIES: O2, C2H4, H2O, CO, CO2, H2, H, O, OH, RP-1 surrogate (C12H24)
  subroutine CKJLR10sp(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    real(8), intent(inout)  :: roi(ns)
    real(8), intent(in)  :: temp 
    real(8), intent(out) :: omegadot(ns)
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(ns), Tdiff
    real(8) :: prod1, prod2, prod3, prod4, prod5, prod6, prod7, prod8

    do is = 1, ns
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
      ! Loop done in order to avoid numerical issues in omegadot evaluation
      ! Very low coi could produce finite prods
      if (coi(is).lt.1d-10) coi(is) = 0.0d0
    enddo

    ! Preliminary
    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    !> 0.5 C2H4 + 1.25 O2 --> 2 CO + 2 H2 - 0.5 C2H4 + 0.25 O2
    prod1 = comp_ch_tabT(1,kf_tab,Tint,Tdiff)*(coi(2)**0.50)*(coi(1)**1.25)

    !> C2H4 + H2O --> 2 CO + 4 H2 - H2O
    prod2 = comp_ch_tabT(2,kf_tab,Tint,Tdiff)*coi(2)*coi(3)

    !> CO + H2O <--> CO2 + H2
    prod3 = comp_ch_tabT(3,kf_tab,Tint,Tdiff)*coi(4)*coi(3)-&
            comp_ch_tabT(3,kb_tab,Tint,Tdiff)*coi(5)*coi(6)

    !> 1/4 H2 + 3/2 O2 <--> H2O + O2 - 3/4 H2
    if (coi(6) < 1.d-10) then
      prod4 = comp_ch_tabT(4,kf_tab,Tint,Tdiff)*(coi(6)**0.25)*(coi(1)**1.50)
    else
      prod4 = comp_ch_tabT(4,kf_tab,Tint,Tdiff)*(coi(6)**0.25)*(coi(1)**1.50) - &
              comp_ch_tabT(4,kb_tab,Tint,Tdiff)*coi(3)*coi(1)*(coi(6)**(-0.75))
    endif

    !> O2 <--> 2O
    prod5 = comp_ch_tabT(5,kf_tab,Tint,Tdiff)*coi(1)-comp_ch_tabT(5,kb_tab,Tint,Tdiff)*coi(8)**2

    !> H2O <--> H + OH
    prod6 = comp_ch_tabT(6,kf_tab,Tint,Tdiff)*coi(3)-comp_ch_tabT(6,kb_tab,Tint,Tdiff)*coi(7)*coi(9)

    !> OH + H2 <--> H + H2O
    prod7 = comp_ch_tabT(7,kf_tab,Tint,Tdiff)*coi(9)*coi(6)- &
            comp_ch_tabT(7,kb_tab,Tint,Tdiff)*coi(7)*coi(3)

    !> C12H24 --> 6C2H4 
    prod8 = comp_ch_tabT(8,kf_tab,Tint,Tdiff)*coi(10)

    !> Chemical Source Terms
    omegadot = 0d0
    omegadot(1) =Wm_tab(1)*(-prod1-0.5*prod4-prod5)             !O2     [kg/(m3*s)]  
    omegadot(2) =Wm_tab(2)*(-prod1-prod2+6*prod8)               !C2H4   [kg/(m3*s)]  
    omegadot(3) =Wm_tab(3)*(-2*prod2-prod3+prod4-prod6+prod7)   !H2O    [kg/(m3*s)]  
    omegadot(4) =Wm_tab(4)*(2*prod1+2*prod2-prod3)              !CO     [kg/(m3*s)]  
    omegadot(5) =Wm_tab(5)*(prod3)                              !CO2    [kg/(m3*s)]  
    omegadot(6) =Wm_tab(6)*(2*prod1+4*prod2+prod3-prod4-prod7)  !H2     [kg/(m3*s)]  
    omegadot(7) =Wm_tab(7)*(prod6+prod7)                        !H      [kg/(m3*s)]  
    omegadot(8) =Wm_tab(8)*(2*prod5)                            !O      [kg/(m3*s)]  
    omegadot(9) =Wm_tab(9)*(prod6-prod7)                        !OH     [kg/(m3*s)]      
    omegadot(10)=Wm_tab(10)*(-prod8)                            !C12H24 [kg/(m3*s)]
    
  end subroutine CKJLR10sp

end module JLRs_mod