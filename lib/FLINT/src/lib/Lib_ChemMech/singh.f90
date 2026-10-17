module singh_mod
  implicit none
contains

  ! Singh (1994) global mechanism for ethylene combustion.
  ! 10 reversible reactions for 9 species.
  subroutine singh(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    real(8), intent(inout) :: roi(ns)
    real(8), intent(in) :: temp 
    real(8), intent(out) :: omegadot(ns)
    ! Local
    integer :: is, T_i, Tint(2)
    integer, parameter :: iH=1, iC2H4=2, iOH=3, iCO=4, iCO2=5
    integer, parameter :: iH2=6, iH2O=7, iO2=8, iO=9
    real(8) :: coi(ns), Tdiff
    real(8) :: M, prodf(11), prodb(11), prod(10)

    do is = 1, ns
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
      if (coi(is).lt.1d-12) coi(is) = 0d0
    enddo
    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! 3rd body eff: 2.5 H2, 16.0 H2O, 1.0 else
    M = sum(coi(iH:iCO2)) + 2.5d0*coi(iH2) + 16d0*coi(iH2O) + coi(iO2) + coi(iO)
    ! C2H4 + O2 <-> 2CO+ 2H2
    prodf(1)=comp_ch_tabT(1,kf_tab,Tint,Tdiff)*coi(iC2H4)*coi(iO2)
    prodb(1)=comp_ch_tabT(1,kb_tab,Tint,Tdiff)*coi(iCO)**2*coi(iH2)**2
    prod(1)=prodf(1)-prodb(1)
    ! CO + O + M <-> CO2 + M
    prodf(2)=comp_ch_tabT(2,kf_tab,Tint,Tdiff)*coi(iCO)*coi(iO)*M
    prodb(2)=comp_ch_tabT(2,kb_tab,Tint,Tdiff)*coi(iCO2)*M
    prod(2)=prodf(2)-prodb(2)
    !CO + OH <-> CO2 + H
    prodf(3)=comp_ch_tabT(3,kf_tab,Tint,Tdiff)*coi(iCO)*coi(iOH)
    prodb(3)=comp_ch_tabT(3,kb_tab,Tint,Tdiff)*coi(iCO2)*coi(iH)
    prod(3)=prodf(3)-prodb(3)
    !H2 + O2 <-> OH+OH
    prodf(4)=comp_ch_tabT(4,kf_tab,Tint,Tdiff)*coi(iH2)*coi(iO2)
    prodb(4)=comp_ch_tabT(4,kb_tab,Tint,Tdiff)*coi(iOH)**2
    prod(4)=prodf(4)-prodb(4)
    !H + O2 <-> OH + O
    prodf(5)=comp_ch_tabT(5,kf_tab,Tint,Tdiff)*coi(iH)*coi(iO2)
    prodb(5)=comp_ch_tabT(5,kb_tab,Tint,Tdiff)*coi(iOH)*coi(iO)
    prod(5)=prodf(5)-prodb(5)
    !OH +H2 <-> H2O + H
    prodf(6)=comp_ch_tabT(6,kf_tab,Tint,Tdiff)*coi(iOH)*coi(iH2)
    prodb(6)=comp_ch_tabT(6,kb_tab,Tint,Tdiff)*coi(iH2O)*coi(iH)
    prod(6)=prodf(6)-prodb(6)
    !O + H2 <-> OH + H
    prodf(7)=comp_ch_tabT(7,kf_tab,Tint,Tdiff)*coi(iO)*coi(iH2)
    prodb(7)=comp_ch_tabT(7,kb_tab,Tint,Tdiff)*coi(iOH)*coi(iH)
    prod(7)=prodf(7)-prodb(7)
    !OH + OH <-> H2O + O
    prodf(8)=comp_ch_tabT(8,kf_tab,Tint,Tdiff)*coi(iOH)**2
    prodb(8)=comp_ch_tabT(8,kb_tab,Tint,Tdiff)*coi(iH2O)*coi(iO)
    prod(8)=prodf(8)-prodb(8)
    !H + H + M <-> H2 + M
    prodf(9)=comp_ch_tabT(9,kf_tab,Tint,Tdiff)*coi(iH)**2*M
    prodb(9)=comp_ch_tabT(9,kb_tab,Tint,Tdiff)*coi(iH2)*M
    prod(9)=prodf(9)-prodb(9)
    !H + OH <-> H2O + M
    prodf(10)=comp_ch_tabT(10,kf_tab,Tint,Tdiff)*coi(iH)*coi(iOH)*M
    prodb(10)=comp_ch_tabT(10,kb_tab,Tint,Tdiff)*coi(iH2O)*M
    prod(10)=prodf(10)-prodb(10)

    ! H
    omegadot(iH)=wm_tab(iH)*(prod(3)-prod(5)+prod(6)+prod(7)-2*prod(9)-prod(10))
    ! C2H4
    omegadot(iC2H4)=wm_tab(iC2H4)*(-prod(1))
    ! OH
    omegadot(iOH)=wm_tab(iOH)*(-prod(3)+2*prod(4)+prod(5)-prod(6)+ &
                                prod(7)-2*prod(8)-prod(10))
    ! CO
    omegadot(iCO)=wm_tab(iCO)*(2*prod(1)-prod(2)-prod(3))
    ! CO2
    omegadot(iCO2)=wm_tab(iCO2)*(prod(2)+prod(3))
    ! H2
    omegadot(iH2)=wm_tab(iH2)*(2*prod(1)-prod(4)-prod(6)-prod(7)+prod(9))
    ! H2O
    omegadot(iH2O)=wm_tab(iH2O)*(prod(6)+prod(8)+prod(10))
    ! O2
    omegadot(iO2)=wm_tab(iO2)*(-prod(1)-prod(4)-prod(5))
    ! O
    omegadot(iO)=wm_tab(iO)*(-prod(2)+prod(5)-prod(7)+prod(8))

  end subroutine singh


  ! Singh (1994) + paraffin cracking
  subroutine Singh_WC32(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    real(8), intent(inout) :: roi(ns)
    real(8), intent(in) :: temp 
    real(8), intent(out) :: omegadot(ns)
    ! Local
    integer :: is, T_i, Tint(2)
    integer, parameter :: iH=1, iC2H4=2, iOH=3, iCO=4, iCO2=5
    integer, parameter :: iH2=6, iH2O=7, iO2=8, iO=9, iC32H66=10
    real(8) :: coi(ns), Tdiff
    real(8) :: M, prodf(11), prodb(11), prod(11)

    do is = 1, ns
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
      if (coi(is).lt.1d-12) coi(is) = 0d0
    enddo
    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! 3rd body eff: 2.5 H2, 16.0 H2O, 1.0 else
    M = sum(coi(iH:iCO2)) + 2.5d0*coi(iH2) + 16d0*coi(iH2O) + coi(iO2) + coi(iO)
    ! C2H4 + O2 <-> 2CO+ 2H2
    prodf(1)=comp_ch_tabT(1,kf_tab,Tint,Tdiff)*coi(iC2H4)*coi(iO2)
    prodb(1)=comp_ch_tabT(1,kb_tab,Tint,Tdiff)*coi(iCO)**2*coi(iH2)**2
    prod(1)=prodf(1)-prodb(1)
    ! CO + O + M <-> CO2 + M
    prodf(2)=comp_ch_tabT(2,kf_tab,Tint,Tdiff)*coi(iCO)*coi(iO)*M
    prodb(2)=comp_ch_tabT(2,kb_tab,Tint,Tdiff)*coi(iCO2)*M
    prod(2)=prodf(2)-prodb(2)
    !CO + OH <-> CO2 + H
    prodf(3)=comp_ch_tabT(3,kf_tab,Tint,Tdiff)*coi(iCO)*coi(iOH)
    prodb(3)=comp_ch_tabT(3,kb_tab,Tint,Tdiff)*coi(iCO2)*coi(iH)
    prod(3)=prodf(3)-prodb(3)
    !H2 + O2 <-> OH+OH
    prodf(4)=comp_ch_tabT(4,kf_tab,Tint,Tdiff)*coi(iH2)*coi(iO2)
    prodb(4)=comp_ch_tabT(4,kb_tab,Tint,Tdiff)*coi(iOH)**2
    prod(4)=prodf(4)-prodb(4)
    !H + O2 <-> OH + O
    prodf(5)=comp_ch_tabT(5,kf_tab,Tint,Tdiff)*coi(iH)*coi(iO2)
    prodb(5)=comp_ch_tabT(5,kb_tab,Tint,Tdiff)*coi(iOH)*coi(iO)
    prod(5)=prodf(5)-prodb(5)
    !OH +H2 <-> H2O + H
    prodf(6)=comp_ch_tabT(6,kf_tab,Tint,Tdiff)*coi(iOH)*coi(iH2)
    prodb(6)=comp_ch_tabT(6,kb_tab,Tint,Tdiff)*coi(iH2O)*coi(iH)
    prod(6)=prodf(6)-prodb(6)
    !O + H2 <-> OH + H
    prodf(7)=comp_ch_tabT(7,kf_tab,Tint,Tdiff)*coi(iO)*coi(iH2)
    prodb(7)=comp_ch_tabT(7,kb_tab,Tint,Tdiff)*coi(iOH)*coi(iH)
    prod(7)=prodf(7)-prodb(7)
    !OH + OH <-> H2O + O
    prodf(8)=comp_ch_tabT(8,kf_tab,Tint,Tdiff)*coi(iOH)**2
    prodb(8)=comp_ch_tabT(8,kb_tab,Tint,Tdiff)*coi(iH2O)*coi(iO)
    prod(8)=prodf(8)-prodb(8)
    !H + H + M <-> H2 + M
    prodf(9)=comp_ch_tabT(9,kf_tab,Tint,Tdiff)*coi(iH)**2*M
    prodb(9)=comp_ch_tabT(9,kb_tab,Tint,Tdiff)*coi(iH2)*M
    prod(9)=prodf(9)-prodb(9)
    !H + OH <-> H2O + M
    prodf(10)=comp_ch_tabT(10,kf_tab,Tint,Tdiff)*coi(iH)*coi(iOH)*M
    prodb(10)=comp_ch_tabT(10,kb_tab,Tint,Tdiff)*coi(iH2O)*M
    prod(10)=prodf(10)-prodb(10)
    !C32H66 --> 16C2H4 + H2
    if (coi(iC32H66)>1d-10) then
      prodf(11)=comp_ch_tabT(11,kf_tab,Tint,Tdiff)*(coi(iC32H66))
    else
      prodf(11)=0d0
    endif
    prodb(11)=0d0
    prod(11)=prodf(11)

    ! H
    omegadot(iH)=wm_tab(iH)*(prod(3)-prod(5)+prod(6)+prod(7)-2*prod(9)-prod(10))
    ! C2H4
    omegadot(iC2H4)=wm_tab(iC2H4)*(-prod(1)+16*prod(11))
    ! OH
    omegadot(iOH)=wm_tab(iOH)*(-prod(3)+2*prod(4)+prod(5)-prod(6)+ &
                                prod(7)-2*prod(8)-prod(10))
    ! CO
    omegadot(iCO)=wm_tab(iCO)*(2*prod(1)-prod(2)-prod(3))
    ! CO2
    omegadot(iCO2)=wm_tab(iCO2)*(prod(2)+prod(3))
    ! H2
    omegadot(iH2)=wm_tab(iH2)*(2*prod(1)-prod(4)-prod(6)-prod(7)+prod(9)+prod(11))
    ! H2O
    omegadot(iH2O)=wm_tab(iH2O)*(prod(6)+prod(8)+prod(10))
    ! O2
    omegadot(iO2)=wm_tab(iO2)*(-prod(1)-prod(4)-prod(5))
    ! O
    omegadot(iO)=wm_tab(iO)*(-prod(2)+prod(5)-prod(7)+prod(8))
    ! C32H66
    omegadot(iC32H66)=wm_tab(iC32H66)*(-prod(11))

  end subroutine Singh_WC32


  ! Variant for propene combustion
  subroutine singhC3H6(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    real(8), intent(inout) :: roi(ns)
    real(8), intent(in) :: temp 
    real(8), intent(out) :: omegadot(ns)
    ! Local
    integer :: is, T_i, Tint(2)
    integer, parameter :: iH=1, iC3H6=2, iOH=3, iCO=4, iCO2=5
    integer, parameter :: iH2=6, iH2O=7, iO2=8, iO=9
    real(8) :: coi(ns), Tdiff
    real(8) :: M, prodf(11), prodb(11), prod(10)

    do is = 1, ns
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
      if (coi(is).lt.1d-12) coi(is) = 0d0
    enddo
    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! 3rd body eff: 2.5 H2, 16.0 H2O, 1.0 else
    M = sum(coi(iH:iCO2)) + 2.5d0*coi(iH2) + 16d0*coi(iH2O) + coi(iO2) + coi(iO)
    ! C3H6 + 1.5 O2 => 3 CO + 3 H2
    prodf(1)=comp_ch_tabT(1,kf_tab,Tint,Tdiff)*coi(iC3H6)*coi(iO2)**1.5d0
    prodb(1)=comp_ch_tabT(1,kb_tab,Tint,Tdiff)*coi(iCO)**3*coi(iH2)**3
    prod(1)=prodf(1)-prodb(1)
    ! CO + O + M <-> CO2 + M
    prodf(2)=comp_ch_tabT(2,kf_tab,Tint,Tdiff)*coi(iCO)*coi(iO)*M
    prodb(2)=comp_ch_tabT(2,kb_tab,Tint,Tdiff)*coi(iCO2)*M
    prod(2)=prodf(2)-prodb(2)
    !CO + OH <-> CO2 + H
    prodf(3)=comp_ch_tabT(3,kf_tab,Tint,Tdiff)*coi(iCO)*coi(iOH)
    prodb(3)=comp_ch_tabT(3,kb_tab,Tint,Tdiff)*coi(iCO2)*coi(iH)
    prod(3)=prodf(3)-prodb(3)
    !H2 + O2 <-> OH+OH
    prodf(4)=comp_ch_tabT(4,kf_tab,Tint,Tdiff)*coi(iH2)*coi(iO2)
    prodb(4)=comp_ch_tabT(4,kb_tab,Tint,Tdiff)*coi(iOH)**2
    prod(4)=prodf(4)-prodb(4)
    !H + O2 <-> OH + O
    prodf(5)=comp_ch_tabT(5,kf_tab,Tint,Tdiff)*coi(iH)*coi(iO2)
    prodb(5)=comp_ch_tabT(5,kb_tab,Tint,Tdiff)*coi(iOH)*coi(iO)
    prod(5)=prodf(5)-prodb(5)
    !OH +H2 <-> H2O + H
    prodf(6)=comp_ch_tabT(6,kf_tab,Tint,Tdiff)*coi(iOH)*coi(iH2)
    prodb(6)=comp_ch_tabT(6,kb_tab,Tint,Tdiff)*coi(iH2O)*coi(iH)
    prod(6)=prodf(6)-prodb(6)
    !O + H2 <-> OH + H
    prodf(7)=comp_ch_tabT(7,kf_tab,Tint,Tdiff)*coi(iO)*coi(iH2)
    prodb(7)=comp_ch_tabT(7,kb_tab,Tint,Tdiff)*coi(iOH)*coi(iH)
    prod(7)=prodf(7)-prodb(7)
    !OH + OH <-> H2O + O
    prodf(8)=comp_ch_tabT(8,kf_tab,Tint,Tdiff)*coi(iOH)**2
    prodb(8)=comp_ch_tabT(8,kb_tab,Tint,Tdiff)*coi(iH2O)*coi(iO)
    prod(8)=prodf(8)-prodb(8)
    !H + H + M <-> H2 + M
    prodf(9)=comp_ch_tabT(9,kf_tab,Tint,Tdiff)*coi(iH)**2*M
    prodb(9)=comp_ch_tabT(9,kb_tab,Tint,Tdiff)*coi(iH2)*M
    prod(9)=prodf(9)-prodb(9)
    !H + OH <-> H2O + M
    prodf(10)=comp_ch_tabT(10,kf_tab,Tint,Tdiff)*coi(iH)*coi(iOH)*M
    prodb(10)=comp_ch_tabT(10,kb_tab,Tint,Tdiff)*coi(iH2O)*M
    prod(10)=prodf(10)-prodb(10)

    ! H
    omegadot(iH)=wm_tab(iH)*(prod(3)-prod(5)+prod(6)+prod(7)-2*prod(9)-prod(10))
    ! C3H6
    omegadot(iC3H6)=wm_tab(iC3H6)*(-prod(1))
    ! OH
    omegadot(iOH)=wm_tab(iOH)*(-prod(3)+2*prod(4)+prod(5)-prod(6)+ &
                                prod(7)-2*prod(8)-prod(10))
    ! CO
    omegadot(iCO)=wm_tab(iCO)*(3*prod(1)-prod(2)-prod(3))
    ! CO2
    omegadot(iCO2)=wm_tab(iCO2)*(prod(2)+prod(3))
    ! H2
    omegadot(iH2)=wm_tab(iH2)*(3*prod(1)-prod(4)-prod(6)-prod(7)+prod(9))
    ! H2O
    omegadot(iH2O)=wm_tab(iH2O)*(prod(6)+prod(8)+prod(10))
    ! O2
    omegadot(iO2)=wm_tab(iO2)*(-1.5*prod(1)-prod(4)-prod(5))
    ! O
    omegadot(iO)=wm_tab(iO)*(-prod(2)+prod(5)-prod(7)+prod(8))

  end subroutine singhC3H6

end module singh_mod