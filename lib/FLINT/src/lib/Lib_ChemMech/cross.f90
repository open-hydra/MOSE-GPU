module cross_mod
implicit none
contains
subroutine cross(roi,temp,omegadot)
use FLINT_Lib_Thermodynamic
use FLINT_Lib_Chemistry_data
use FLINT_Lib_Chemistry_falloff
implicit none
real(8), intent(inout)  :: roi(ns)
real(8), intent(in)  :: temp
real(8), intent(out) :: omegadot(ns) 

real(8) :: coi(ns), Tdiff 
real(8) :: M !< Third body
integer :: is, T_i, Tint(2)
real(8) :: prodf(1:33), prodb(1:33)
real(8) :: k(2) !< Troe rate coefficients


do is = 1, ns 
 coi(is)=roi(is)/Wm_tab(is) ! kmol/m^3
enddo 
T_i = int(temp) 
Tdiff  = temp-T_i 
Tint(1) = T_i 
Tint(2) = T_i + 1 
! reac n. 1: CL2 + H <=> CL + HCL
prodf(1)=f_kf(1,Tint,Tdiff)*(coi(9))*(coi(11))
prodb(1)=f_kb(1,Tint,Tdiff)*(coi(8))*(coi(16))
! reac n. 2: H2 + OH <=> H + H2O
prodf(2)=f_kf(2,Tint,Tdiff)*(coi(12))*(coi(19))
prodb(2)=f_kb(2,Tint,Tdiff)*(coi(11))*(coi(13))
! reac n. 3: C2H2 + OH <=> CH3 + CO
prodf(3)=f_kf(3,Tint,Tdiff)*(coi(1))*(coi(19))
prodb(3)=f_kb(3,Tint,Tdiff)*(coi(4))*(coi(6))
! reac n. 4: CH4 + CL <=> CH3 + HCL
prodf(4)=f_kf(4,Tint,Tdiff)*(coi(5))*(coi(8))
prodb(4)=f_kb(4,Tint,Tdiff)*(coi(4))*(coi(16))
! reac n. 5: CH4 + H <=> CH3 + H2
prodf(5)=f_kf(5,Tint,Tdiff)*(coi(5))*(coi(11))
prodb(5)=f_kb(5,Tint,Tdiff)*(coi(4))*(coi(12))
! reac n. 6: CH4 + OH <=> CH3 + H2O
prodf(6)=f_kf(6,Tint,Tdiff)*(coi(5))*(coi(19))
prodb(6)=f_kb(6,Tint,Tdiff)*(coi(4))*(coi(13))
! reac n. 7: CH3 + H + M <=> CH4 + M
M=coi(1)+coi(2)+coi(3)+coi(4)+coi(5)*2+coi(6)*1.5+coi(7)*2+coi(8)+coi(9)+coi(10)+coi(11)+coi(12)*2+coi(13)*6+coi(14)+coi(15)+coi(16)+coi(17)+coi(18)+coi(19)
prodf(7)=f_kf(7,Tint,Tdiff)*(coi(4))*(coi(11))*M
prodb(7)=f_kb(7,Tint,Tdiff)*(coi(5))*M
! reac n. 8: HCO + M <=> CO + H + M
M=coi(1)+coi(2)+coi(3)+coi(4)+coi(5)*2+coi(6)*1.5+coi(7)*2+coi(8)+coi(9)+coi(10)+coi(11)+coi(12)*2+coi(14)+coi(15)+coi(16)+coi(17)+coi(18)+coi(19)
prodf(8)=f_kf(8,Tint,Tdiff)*(coi(15))*M
prodb(8)=f_kb(8,Tint,Tdiff)*(coi(6))*(coi(11))*M
! reac n. 9: CO + OH <=> CO2 + H
prodf(9)=f_kf(9,Tint,Tdiff)*(coi(6))*(coi(19))
prodb(9)=f_kb(9,Tint,Tdiff)*(coi(7))*(coi(11))
! reac n. 10: CLO + CO <=> CL + CO2
prodf(10)=f_kf(10,Tint,Tdiff)*(coi(10))*(coi(6))
prodb(10)=f_kb(10,Tint,Tdiff)*(coi(8))*(coi(7))
! reac n. 11: C2H2 + O <=> CH2 + CO
prodf(11)=f_kf(11,Tint,Tdiff)*(coi(1))*(coi(18))
prodb(11)=f_kb(11,Tint,Tdiff)*(coi(3))*(coi(6))
! reac n. 12: CH2 + H2 <=> CH3 + H
prodf(12)=f_kf(12,Tint,Tdiff)*(coi(3))*(coi(12))
prodb(12)=f_kb(12,Tint,Tdiff)*(coi(4))*(coi(11))
! reac n. 13: CH2 + H + M <=> CH3 + M
M=coi(1)+coi(2)+coi(3)+coi(4)+coi(5)*2+coi(6)*1.5+coi(7)*2+coi(8)+coi(9)+coi(10)+coi(11)+coi(12)*2+coi(13)*6+coi(14)+coi(15)+coi(16)+coi(17)+coi(18)+coi(19)
prodf(13)=f_kf(13,Tint,Tdiff)*(coi(3))*(coi(11))*M
prodb(13)=f_kb(13,Tint,Tdiff)*(coi(4))*M
! reac n. 14: CH4 + O <=> CH3 + OH
prodf(14)=f_kf(14,Tint,Tdiff)*(coi(5))*(coi(18))
prodb(14)=f_kb(14,Tint,Tdiff)*(coi(4))*(coi(19))
! reac n. 15: CH3 + O <=> CH2O + H
prodf(15)=f_kf(15,Tint,Tdiff)*(coi(4))*(coi(18))
prodb(15)=f_kb(15,Tint,Tdiff)*(coi(2))*(coi(11))
! reac n. 16: CH2 + O <=> H + HCO
prodf(16)=f_kf(16,Tint,Tdiff)*(coi(3))*(coi(18))
prodb(16)=f_kb(16,Tint,Tdiff)*(coi(11))*(coi(15))
! reac n. 17: CH3 + OH <=> CH2 + H2O
prodf(17)=f_kf(17,Tint,Tdiff)*(coi(4))*(coi(19))
prodb(17)=f_kb(17,Tint,Tdiff)*(coi(3))*(coi(13))
! reac n. 18: CH2 + OH <=> CH2O + H
prodf(18)=f_kf(18,Tint,Tdiff)*(coi(3))*(coi(19))
prodb(18)=f_kb(18,Tint,Tdiff)*(coi(2))*(coi(11))
! reac n. 19: CL + H + M <=> HCL + M
M=coi(1)+coi(2)+coi(3)+coi(4)+coi(5)*2+coi(6)*1.5+coi(7)*2+coi(8)+coi(9)+coi(10)+coi(11)+coi(12)*2+coi(13)+coi(14)+coi(15)+coi(16)+coi(17)+coi(18)+coi(19)
prodf(19)=f_kf(19,Tint,Tdiff)*(coi(8))*(coi(11))*M
prodb(19)=f_kb(19,Tint,Tdiff)*(coi(16))*M
! reac n. 20: 2 CL + M <=> CL2 + M
M=coi(1)+coi(2)+coi(3)+coi(4)+coi(5)*2+coi(6)*1.5+coi(7)*2+coi(8)+coi(9)+coi(10)+coi(11)+coi(12)*2+coi(13)+coi(14)+coi(15)+coi(16)+coi(17)+coi(18)+coi(19)
prodf(20)=f_kf(20,Tint,Tdiff)*(coi(8)*coi(8))*M
prodb(20)=f_kb(20,Tint,Tdiff)*(coi(9))*M
! reac n. 21: CH2O + CL <=> HCL + HCO
prodf(21)=f_kf(21,Tint,Tdiff)*(coi(2))*(coi(8))
prodb(21)=f_kb(21,Tint,Tdiff)*(coi(16))*(coi(15))
! reac n. 22: CH3 + CLO <=> CH2O + CL + H
prodf(22)=f_kf(22,Tint,Tdiff)*(coi(4))*(coi(10))
prodb(22)=f_kb(22,Tint,Tdiff)*(coi(2))*(coi(8))*(coi(11))
! reac n. 23: CH3 + CLO <=> CH2O + HCL
prodf(23)=f_kf(23,Tint,Tdiff)*(coi(4))*(coi(10))
prodb(23)=f_kb(23,Tint,Tdiff)*(coi(2))*(coi(16))
! reac n. 24: H + HCL <=> CL + H2
prodf(24)=f_kf(24,Tint,Tdiff)*(coi(11))*(coi(16))
prodb(24)=f_kb(24,Tint,Tdiff)*(coi(8))*(coi(12))
! reac n. 25: HCL + O <=> CL + OH
prodf(25)=f_kf(25,Tint,Tdiff)*(coi(16))*(coi(18))
prodb(25)=f_kb(25,Tint,Tdiff)*(coi(8))*(coi(19))
! reac n. 26: 2 H + M <=> H2 + M
M=coi(1)+coi(2)+coi(3)+coi(4)+coi(5)*2+coi(6)+coi(7)+coi(8)+coi(9)+coi(10)+coi(11)+coi(12)+coi(13)+coi(14)+coi(15)+coi(16)+coi(17)+coi(18)+coi(19)
prodf(26)=f_kf(26,Tint,Tdiff)*(coi(11)*coi(11))*M
prodb(26)=f_kb(26,Tint,Tdiff)*(coi(12))*M
! reac n. 27: 2 H + H2 <=> H2 + H2
M=coi(12)
prodf(27)=f_kf(27,Tint,Tdiff)*(coi(11)*coi(11))*M
prodb(27)=f_kb(27,Tint,Tdiff)*(coi(12))*M
! reac n. 28: 2 H + H2O <=> H2 + H2O
M=coi(13)
prodf(28)=f_kf(28,Tint,Tdiff)*(coi(11)*coi(11))*M
prodb(28)=f_kb(28,Tint,Tdiff)*(coi(12))*M
! reac n. 29: H + HCO <=> CO + H2
prodf(29)=f_kf(29,Tint,Tdiff)*(coi(11))*(coi(15))
prodb(29)=f_kb(29,Tint,Tdiff)*(coi(6))*(coi(12))
! reac n. 30: 2 OH <=> H2O + O
prodf(30)=f_kf(30,Tint,Tdiff)*(coi(19)*coi(19))
prodb(30)=f_kb(30,Tint,Tdiff)*(coi(13))*(coi(18))
! reac n. 31: CH2 + CH4 <=> 2 CH3
prodf(31)=f_kf(31,Tint,Tdiff)*(coi(3))*(coi(5))
prodb(31)=f_kb(31,Tint,Tdiff)*(coi(4)*coi(4))
! reac n. 32: CH3 + HCO <=> CH4 + CO
prodf(32)=f_kf(32,Tint,Tdiff)*(coi(4))*(coi(15))
prodb(32)=f_kb(32,Tint,Tdiff)*(coi(5))*(coi(6))
! reac n. 33: CH2O + CH3 <=> CH4 + HCO
prodf(33)=f_kf(33,Tint,Tdiff)*(coi(2))*(coi(4))
prodb(33)=f_kb(33,Tint,Tdiff)*(coi(5))*(coi(15))
! species source terms
omegadot(1)=Wm_tab(1)*(+(0.0-1.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(11)-prodb(11)))
omegadot(2)=Wm_tab(2)*(+(1.0-0.0)*(prodf(15)-prodb(15))+(1.0-0.0)*(prodf(18)-prodb(18))+(0.0-1.0)*(prodf(21)-prodb(21))+(1.0-0.0)*(prodf(22)-prodb(22))+(1.0-0.0)*(prodf(23)-prodb(23))+(0.0-1.0)*(prodf(33)-prodb(33)))
omegadot(3)=Wm_tab(3)*(+(1.0-0.0)*(prodf(11)-prodb(11))+(0.0-1.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(13)-prodb(13))+(0.0-1.0)*(prodf(16)-prodb(16))+(1.0-0.0)*(prodf(17)-prodb(17))+(0.0-1.0)*(prodf(18)-prodb(18))+(0.0-1.0)*(prodf(31)-prodb(31)))
omegadot(4)=Wm_tab(4)*(+(1.0-0.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(5)-prodb(5))+(1.0-0.0)*(prodf(6)-prodb(6))+(0.0-1.0)*(prodf(7)-prodb(7))+(1.0-0.0)*(prodf(12)-prodb(12))+(1.0-0.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(14)-prodb(14))+(0.0-1.0)*(prodf(15)-prodb(15))+(0.0-1.0)*(prodf(17)-prodb(17))+(0.0-1.0)*(prodf(22)-prodb(22))+(0.0-1.0)*(prodf(23)-prodb(23))+(2.0-0.0)*(prodf(31)-prodb(31))+(0.0-1.0)*(prodf(32)-prodb(32))+(0.0-1.0)*(prodf(33)-prodb(33)))
omegadot(5)=Wm_tab(5)*(+(0.0-1.0)*(prodf(4)-prodb(4))+(0.0-1.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(7)-prodb(7))+(0.0-1.0)*(prodf(14)-prodb(14))+(0.0-1.0)*(prodf(31)-prodb(31))+(1.0-0.0)*(prodf(32)-prodb(32))+(1.0-0.0)*(prodf(33)-prodb(33)))
omegadot(6)=Wm_tab(6)*(+(1.0-0.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(8)-prodb(8))+(0.0-1.0)*(prodf(9)-prodb(9))+(0.0-1.0)*(prodf(10)-prodb(10))+(1.0-0.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(29)-prodb(29))+(1.0-0.0)*(prodf(32)-prodb(32)))
omegadot(7)=Wm_tab(7)*(+(1.0-0.0)*(prodf(9)-prodb(9))+(1.0-0.0)*(prodf(10)-prodb(10)))
omegadot(8)=Wm_tab(8)*(+(1.0-0.0)*(prodf(1)-prodb(1))+(0.0-1.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(10)-prodb(10))+(0.0-1.0)*(prodf(19)-prodb(19))+(0.0-2.0)*(prodf(20)-prodb(20))+(0.0-1.0)*(prodf(21)-prodb(21))+(1.0-0.0)*(prodf(22)-prodb(22))+(1.0-0.0)*(prodf(24)-prodb(24))+(1.0-0.0)*(prodf(25)-prodb(25)))
omegadot(9)=Wm_tab(9)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(20)-prodb(20)))
omegadot(10)=Wm_tab(10)*(+(0.0-1.0)*(prodf(10)-prodb(10))+(0.0-1.0)*(prodf(22)-prodb(22))+(0.0-1.0)*(prodf(23)-prodb(23)))
omegadot(11)=Wm_tab(11)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(7)-prodb(7))+(1.0-0.0)*(prodf(8)-prodb(8))+(1.0-0.0)*(prodf(9)-prodb(9))+(1.0-0.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(15)-prodb(15))+(1.0-0.0)*(prodf(16)-prodb(16))+(1.0-0.0)*(prodf(18)-prodb(18))+(0.0-1.0)*(prodf(19)-prodb(19))+(1.0-0.0)*(prodf(22)-prodb(22))+(0.0-1.0)*(prodf(24)-prodb(24))+(0.0-2.0)*(prodf(26)-prodb(26))+(0.0-2.0)*(prodf(27)-prodb(27))+(0.0-2.0)*(prodf(28)-prodb(28))+(0.0-1.0)*(prodf(29)-prodb(29)))
omegadot(12)=Wm_tab(12)*(+(0.0-1.0)*(prodf(2)-prodb(2))+(1.0-0.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(12)-prodb(12))+(1.0-0.0)*(prodf(24)-prodb(24))+(1.0-0.0)*(prodf(26)-prodb(26))+(1.0-0.0)*(prodf(27)-prodb(27))+(1.0-0.0)*(prodf(28)-prodb(28))+(1.0-0.0)*(prodf(29)-prodb(29)))
omegadot(13)=Wm_tab(13)*(+(1.0-0.0)*(prodf(2)-prodb(2))+(1.0-0.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(17)-prodb(17))+(1.0-0.0)*(prodf(30)-prodb(30)))
omegadot(14)=0.d0
omegadot(15)=Wm_tab(15)*(+(0.0-1.0)*(prodf(8)-prodb(8))+(1.0-0.0)*(prodf(16)-prodb(16))+(1.0-0.0)*(prodf(21)-prodb(21))+(0.0-1.0)*(prodf(29)-prodb(29))+(0.0-1.0)*(prodf(32)-prodb(32))+(1.0-0.0)*(prodf(33)-prodb(33)))
omegadot(16)=Wm_tab(16)*(+(1.0-0.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(19)-prodb(19))+(1.0-0.0)*(prodf(21)-prodb(21))+(1.0-0.0)*(prodf(23)-prodb(23))+(0.0-1.0)*(prodf(24)-prodb(24))+(0.0-1.0)*(prodf(25)-prodb(25)))
omegadot(17)=0.d0
omegadot(18)=Wm_tab(18)*(+(0.0-1.0)*(prodf(11)-prodb(11))+(0.0-1.0)*(prodf(14)-prodb(14))+(0.0-1.0)*(prodf(15)-prodb(15))+(0.0-1.0)*(prodf(16)-prodb(16))+(0.0-1.0)*(prodf(25)-prodb(25))+(1.0-0.0)*(prodf(30)-prodb(30)))
omegadot(19)=Wm_tab(19)*(+(0.0-1.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(6)-prodb(6))+(0.0-1.0)*(prodf(9)-prodb(9))+(1.0-0.0)*(prodf(14)-prodb(14))+(0.0-1.0)*(prodf(17)-prodb(17))+(0.0-1.0)*(prodf(18)-prodb(18))+(1.0-0.0)*(prodf(25)-prodb(25))+(0.0-2.0)*(prodf(30)-prodb(30)))
end subroutine cross
end module cross_mod
