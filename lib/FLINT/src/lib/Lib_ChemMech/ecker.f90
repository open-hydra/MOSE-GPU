module ecker_mod
implicit none
contains
subroutine ecker(roi,temp,omegadot)
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
real(8) :: prodf(1:28), prodb(1:28)
real(8) :: k(2) !< Troe rate coefficients


do is = 1, ns 
 coi(is)=roi(is)/Wm_tab(is) ! kmol/m^3
enddo 
T_i = int(temp) 
Tdiff  = temp-T_i 
Tint(1) = T_i 
Tint(2) = T_i + 1 
! reac n. 1: H2 + O2 <=> H + HO2
prodf(1)=f_kf(1,Tint,Tdiff)*(coi(1))*(coi(2))
prodb(1)=f_kb(1,Tint,Tdiff)*(coi(4))*(coi(3))
! reac n. 2: H + O2 <=> O + OH
prodf(2)=f_kf(2,Tint,Tdiff)*(coi(4))*(coi(2))
prodb(2)=f_kb(2,Tint,Tdiff)*(coi(5))*(coi(6))
! reac n. 3: H2 + O <=> H + OH
prodf(3)=f_kf(3,Tint,Tdiff)*(coi(1))*(coi(5))
prodb(3)=f_kb(3,Tint,Tdiff)*(coi(4))*(coi(6))
! reac n. 4: H2 + OH <=> H + H2O
prodf(4)=f_kf(4,Tint,Tdiff)*(coi(1))*(coi(6))
prodb(4)=f_kb(4,Tint,Tdiff)*(coi(4))*(coi(7))
! reac n. 5: 2 OH <=> H2O + O
prodf(5)=f_kf(5,Tint,Tdiff)*(coi(6)*coi(6))
prodb(5)=f_kb(5,Tint,Tdiff)*(coi(7))*(coi(5))
! reac n. 6: H + OH + M <=> H2O + M
M=coi(1)+coi(2)+coi(3)+coi(4)+coi(5)+coi(6)+coi(7)*6+coi(8)+coi(9)+coi(10)+coi(11)+coi(12)+coi(13)+coi(14)
prodf(6)=f_kf(6,Tint,Tdiff)*(coi(4))*(coi(6))*M
prodb(6)=f_kb(6,Tint,Tdiff)*(coi(7))*M
! reac n. 7: 2 H + M <=> H2 + M
M=coi(1)*2+coi(2)+coi(3)+coi(4)+coi(5)+coi(6)+coi(7)*6+coi(8)+coi(9)+coi(10)+coi(11)+coi(12)+coi(13)+coi(14)
prodf(7)=f_kf(7,Tint,Tdiff)*(coi(4)*coi(4))*M
prodb(7)=f_kb(7,Tint,Tdiff)*(coi(1))*M
! reac n. 8: H + O + M <=> OH + M
M=coi(1)+coi(2)+coi(3)+coi(4)+coi(5)+coi(6)+coi(7)*5+coi(8)+coi(9)+coi(10)+coi(11)+coi(12)+coi(13)+coi(14)
prodf(8)=f_kf(8,Tint,Tdiff)*(coi(4))*(coi(5))*M
prodb(8)=f_kb(8,Tint,Tdiff)*(coi(6))*M
! reac n. 9: H + O2 + M <=> HO2 + M
M=coi(1)*2+coi(2)+coi(3)+coi(4)+coi(5)+coi(6)+coi(7)*16+coi(8)+coi(9)+coi(10)+coi(11)+coi(12)+coi(13)+coi(14)
prodf(9)=f_kf(9,Tint,Tdiff)*(coi(4))*(coi(2))*M
prodb(9)=f_kb(9,Tint,Tdiff)*(coi(3))*M
! reac n. 10: H + HO2 <=> 2 OH
prodf(10)=f_kf(10,Tint,Tdiff)*(coi(4))*(coi(3))
prodb(10)=f_kb(10,Tint,Tdiff)*(coi(6)*coi(6))
! reac n. 11: H2 + HO2 <=> H2O + OH
prodf(11)=f_kf(11,Tint,Tdiff)*(coi(1))*(coi(3))
prodb(11)=f_kb(11,Tint,Tdiff)*(coi(7))*(coi(6))
! reac n. 12: HO2 + O <=> O2 + OH
prodf(12)=f_kf(12,Tint,Tdiff)*(coi(3))*(coi(5))
prodb(12)=f_kb(12,Tint,Tdiff)*(coi(2))*(coi(6))
! reac n. 13: HO2 + OH <=> H2O + O2
prodf(13)=f_kf(13,Tint,Tdiff)*(coi(3))*(coi(6))
prodb(13)=f_kb(13,Tint,Tdiff)*(coi(7))*(coi(2))
! reac n. 14: 2 HO2 <=> H2O2 + O2
prodf(14)=f_kf(14,Tint,Tdiff)*(coi(3)*coi(3))
prodb(14)=f_kb(14,Tint,Tdiff)*(coi(8))*(coi(2))
! reac n. 15: H + H2O2 <=> H2 + HO2
prodf(15)=f_kf(15,Tint,Tdiff)*(coi(4))*(coi(8))
prodb(15)=f_kb(15,Tint,Tdiff)*(coi(1))*(coi(3))
! reac n. 16: H2O2 + O <=> HO2 + OH
prodf(16)=f_kf(16,Tint,Tdiff)*(coi(8))*(coi(5))
prodb(16)=f_kb(16,Tint,Tdiff)*(coi(3))*(coi(6))
! reac n. 17: H2O2 + OH <=> H2O + HO2
prodf(17)=f_kf(17,Tint,Tdiff)*(coi(8))*(coi(6))
prodb(17)=f_kb(17,Tint,Tdiff)*(coi(7))*(coi(3))
! reac n. 18: H2O2 + M <=> 2 OH + M
M=coi(1)+coi(2)+coi(3)+coi(4)+coi(5)+coi(6)+coi(7)*15+coi(8)+coi(9)+coi(10)+coi(11)+coi(12)+coi(13)+coi(14)
prodf(18)=f_kf(18,Tint,Tdiff)*(coi(8))*M
prodb(18)=f_kb(18,Tint,Tdiff)*(coi(6)*coi(6))*M
! reac n. 19: 2 O + M <=> O2 + M
M=sum(coi(1:14))
prodf(19)=f_kf(19,Tint,Tdiff)*(coi(5)*coi(5))*M
prodb(19)=f_kb(19,Tint,Tdiff)*(coi(2))*M
! reac n. 20: CO + OH <=> CO2 + H
prodf(20)=f_kf(20,Tint,Tdiff)*(coi(9))*(coi(6))
prodb(20)=f_kb(20,Tint,Tdiff)*(coi(10))*(coi(4))
! reac n. 21: CO + O2 <=> CO2 + O
prodf(21)=f_kf(21,Tint,Tdiff)*(coi(9))*(coi(2))
prodb(21)=f_kb(21,Tint,Tdiff)*(coi(10))*(coi(5))
! reac n. 22: CO + O + M <=> CO2 + M
M=coi(1)*2.5+coi(2)+coi(3)+coi(4)+coi(5)+coi(6)+coi(7)*12+coi(8)+coi(9)*1.9+coi(10)*3.8+coi(11)+coi(12)+coi(13)+coi(14)
prodf(22)=f_kf(22,Tint,Tdiff)*(coi(9))*(coi(5))*M
prodb(22)=f_kb(22,Tint,Tdiff)*(coi(10))*M
! reac n. 23: H + HCL <=> CL + H2
prodf(23)=f_kf(23,Tint,Tdiff)*(coi(4))*(coi(11))
prodb(23)=f_kb(23,Tint,Tdiff)*(coi(13))*(coi(1))
! reac n. 24: CL2 + H <=> CL + HCL
prodf(24)=f_kf(24,Tint,Tdiff)*(coi(12))*(coi(4))
prodb(24)=f_kb(24,Tint,Tdiff)*(coi(13))*(coi(11))
! reac n. 25: HCL + OH <=> CL + H2O
prodf(25)=f_kf(25,Tint,Tdiff)*(coi(11))*(coi(6))
prodb(25)=f_kb(25,Tint,Tdiff)*(coi(13))*(coi(7))
! reac n. 26: HCL + O <=> CL + OH
prodf(26)=f_kf(26,Tint,Tdiff)*(coi(11))*(coi(5))
prodb(26)=f_kb(26,Tint,Tdiff)*(coi(13))*(coi(6))
! reac n. 27: 2 CL + M <=> CL2 + M
M=sum(coi(1:14))
prodf(27)=f_kf(27,Tint,Tdiff)*(coi(13)*coi(13))*M
prodb(27)=f_kb(27,Tint,Tdiff)*(coi(12))*M
! reac n. 28: CL + H + M <=> HCL + M
M=sum(coi(1:14))
prodf(28)=f_kf(28,Tint,Tdiff)*(coi(13))*(coi(4))*M
prodb(28)=f_kb(28,Tint,Tdiff)*(coi(11))*M
! species source terms
omegadot(1)=Wm_tab(1)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(0.0-1.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(7)-prodb(7))+(0.0-1.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(15)-prodb(15))+(1.0-0.0)*(prodf(23)-prodb(23)))
omegadot(2)=Wm_tab(2)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(0.0-1.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(9)-prodb(9))+(1.0-0.0)*(prodf(12)-prodb(12))+(1.0-0.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(14)-prodb(14))+(1.0-0.0)*(prodf(19)-prodb(19))+(0.0-1.0)*(prodf(21)-prodb(21)))
omegadot(3)=Wm_tab(3)*(+(1.0-0.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(9)-prodb(9))+(0.0-1.0)*(prodf(10)-prodb(10))+(0.0-1.0)*(prodf(11)-prodb(11))+(0.0-1.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(13)-prodb(13))+(0.0-2.0)*(prodf(14)-prodb(14))+(1.0-0.0)*(prodf(15)-prodb(15))+(1.0-0.0)*(prodf(16)-prodb(16))+(1.0-0.0)*(prodf(17)-prodb(17)))
omegadot(4)=Wm_tab(4)*(+(1.0-0.0)*(prodf(1)-prodb(1))+(0.0-1.0)*(prodf(2)-prodb(2))+(1.0-0.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(4)-prodb(4))+(0.0-1.0)*(prodf(6)-prodb(6))+(0.0-2.0)*(prodf(7)-prodb(7))+(0.0-1.0)*(prodf(8)-prodb(8))+(0.0-1.0)*(prodf(9)-prodb(9))+(0.0-1.0)*(prodf(10)-prodb(10))+(0.0-1.0)*(prodf(15)-prodb(15))+(1.0-0.0)*(prodf(20)-prodb(20))+(0.0-1.0)*(prodf(23)-prodb(23))+(0.0-1.0)*(prodf(24)-prodb(24))+(0.0-1.0)*(prodf(28)-prodb(28)))
omegadot(5)=Wm_tab(5)*(+(1.0-0.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(8)-prodb(8))+(0.0-1.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(16)-prodb(16))+(0.0-2.0)*(prodf(19)-prodb(19))+(1.0-0.0)*(prodf(21)-prodb(21))+(0.0-1.0)*(prodf(22)-prodb(22))+(0.0-1.0)*(prodf(26)-prodb(26)))
omegadot(6)=Wm_tab(6)*(+(1.0-0.0)*(prodf(2)-prodb(2))+(1.0-0.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(4)-prodb(4))+(0.0-2.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(8)-prodb(8))+(2.0-0.0)*(prodf(10)-prodb(10))+(1.0-0.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(16)-prodb(16))+(0.0-1.0)*(prodf(17)-prodb(17))+(2.0-0.0)*(prodf(18)-prodb(18))+(0.0-1.0)*(prodf(20)-prodb(20))+(0.0-1.0)*(prodf(25)-prodb(25))+(1.0-0.0)*(prodf(26)-prodb(26)))
omegadot(7)=Wm_tab(7)*(+(1.0-0.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(5)-prodb(5))+(1.0-0.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(17)-prodb(17))+(1.0-0.0)*(prodf(25)-prodb(25)))
omegadot(8)=Wm_tab(8)*(+(1.0-0.0)*(prodf(14)-prodb(14))+(0.0-1.0)*(prodf(15)-prodb(15))+(0.0-1.0)*(prodf(16)-prodb(16))+(0.0-1.0)*(prodf(17)-prodb(17))+(0.0-1.0)*(prodf(18)-prodb(18)))
omegadot(9)=Wm_tab(9)*(+(0.0-1.0)*(prodf(20)-prodb(20))+(0.0-1.0)*(prodf(21)-prodb(21))+(0.0-1.0)*(prodf(22)-prodb(22)))
omegadot(10)=Wm_tab(10)*(+(1.0-0.0)*(prodf(20)-prodb(20))+(1.0-0.0)*(prodf(21)-prodb(21))+(1.0-0.0)*(prodf(22)-prodb(22)))
omegadot(11)=Wm_tab(11)*(+(0.0-1.0)*(prodf(23)-prodb(23))+(1.0-0.0)*(prodf(24)-prodb(24))+(0.0-1.0)*(prodf(25)-prodb(25))+(0.0-1.0)*(prodf(26)-prodb(26))+(1.0-0.0)*(prodf(28)-prodb(28)))
omegadot(12)=Wm_tab(12)*(+(0.0-1.0)*(prodf(24)-prodb(24))+(1.0-0.0)*(prodf(27)-prodb(27)))
omegadot(13)=Wm_tab(13)*(+(1.0-0.0)*(prodf(23)-prodb(23))+(1.0-0.0)*(prodf(24)-prodb(24))+(1.0-0.0)*(prodf(25)-prodb(25))+(1.0-0.0)*(prodf(26)-prodb(26))+(0.0-2.0)*(prodf(27)-prodb(27))+(0.0-1.0)*(prodf(28)-prodb(28)))
omegadot(14)=0.d0
end subroutine ecker
end module ecker_mod
