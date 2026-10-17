module TSRCDF13_mod
implicit none
contains
subroutine TSRCDF13(roi,temp,omegadot)
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
real(8) :: prodf(1:46), prodb(1:46)
real(8) :: k(2) !< Troe rate coefficients


do is = 1, ns 
 coi(is)=roi(is)/Wm_tab(is) ! kmol/m^3
enddo 
T_i = int(temp) 
Tdiff  = temp-T_i 
Tint(1) = T_i 
Tint(2) = T_i + 1 
! reac n. 1: 2 O + M <=> O2 + M
M=coi(1)+coi(2)+coi(3)*2.4+coi(4)*3.6+coi(5)*1.75+coi(6)*15.4+coi(7)*2+coi(8)+coi(9)+coi(10)+coi(11)+coi(12)+coi(13)
prodf(1)=f_kf(1,Tint,Tdiff)*(coi(13)*coi(13))*M
prodb(1)=f_kb(1,Tint,Tdiff)*(coi(10))*M
! reac n. 2: H + O + M <=> OH + M
M=coi(1)+coi(2)+coi(3)*2+coi(4)*2+coi(5)*1.5+coi(6)*6+coi(7)*2+coi(8)+coi(9)+coi(10)+coi(11)+coi(12)+coi(13)
prodf(2)=f_kf(2,Tint,Tdiff)*(coi(8))*(coi(13))*M
prodb(2)=f_kb(2,Tint,Tdiff)*(coi(1))*M
! reac n. 3: H2 + O <=> H + OH
prodf(3)=f_kf(3,Tint,Tdiff)*(coi(3))*(coi(13))
prodb(3)=f_kb(3,Tint,Tdiff)*(coi(8))*(coi(1))
! reac n. 4: HO2 + O <=> O2 + OH
prodf(4)=f_kf(4,Tint,Tdiff)*(coi(2))*(coi(13))
prodb(4)=f_kb(4,Tint,Tdiff)*(coi(10))*(coi(1))
! reac n. 5: CH3 + O <=> CH2O + H
prodf(5)=f_kf(5,Tint,Tdiff)*(coi(9))*(coi(13))
prodb(5)=f_kb(5,Tint,Tdiff)*(coi(11))*(coi(8))
! reac n. 6: CH4 + O <=> CH3 + OH
prodf(6)=f_kf(6,Tint,Tdiff)*(coi(7))*(coi(13))
prodb(6)=f_kb(6,Tint,Tdiff)*(coi(9))*(coi(1))
! reac n. 7: CO + O + M <=> CO2 + M
M=coi(1)+coi(2)+coi(3)*2+coi(4)*3.5+coi(5)*1.5+coi(6)*6+coi(7)*2+coi(8)+coi(9)+coi(10)*6+coi(11)+coi(12)+coi(13)
prodf(7)=f_kf(7,Tint,Tdiff)*(coi(5))*(coi(13))*M
prodb(7)=f_kb(7,Tint,Tdiff)*(coi(4))*M
! reac n. 8: HCO + O <=> CO + OH
prodf(8)=f_kf(8,Tint,Tdiff)*(coi(12))*(coi(13))
prodb(8)=f_kb(8,Tint,Tdiff)*(coi(5))*(coi(1))
! reac n. 9: HCO + O <=> CO2 + H
prodf(9)=f_kf(9,Tint,Tdiff)*(coi(12))*(coi(13))
prodb(9)=f_kb(9,Tint,Tdiff)*(coi(4))*(coi(8))
! reac n. 10: CH2O + O <=> HCO + OH
prodf(10)=f_kf(10,Tint,Tdiff)*(coi(11))*(coi(13))
prodb(10)=f_kb(10,Tint,Tdiff)*(coi(12))*(coi(1))
! reac n. 11: CO + O2 <=> CO2 + O
prodf(11)=f_kf(11,Tint,Tdiff)*(coi(5))*(coi(10))
prodb(11)=f_kb(11,Tint,Tdiff)*(coi(4))*(coi(13))
! reac n. 12: CH2O + O2 <=> HCO + HO2
prodf(12)=f_kf(12,Tint,Tdiff)*(coi(11))*(coi(10))
prodb(12)=f_kb(12,Tint,Tdiff)*(coi(12))*(coi(2))
! reac n. 13: H + O2 + M <=> HO2 + M
M=coi(1)+coi(2)+coi(3)+coi(4)*1.5+coi(5)*0.75+coi(7)+coi(8)+coi(9)+coi(11)+coi(12)+coi(13)
prodf(13)=f_kf(13,Tint,Tdiff)*(coi(8))*(coi(10))*M
prodb(13)=f_kb(13,Tint,Tdiff)*(coi(2))*M
! reac n. 14: H + O2 + O2 <=> HO2 + O2
M=coi(10)
prodf(14)=f_kf(14,Tint,Tdiff)*(coi(8))*(coi(10))*M
prodb(14)=f_kb(14,Tint,Tdiff)*(coi(2))*M
! reac n. 15: H + O2 + H2O <=> HO2 + H2O
M=coi(6)
prodf(15)=f_kf(15,Tint,Tdiff)*(coi(8))*(coi(10))*M
prodb(15)=f_kb(15,Tint,Tdiff)*(coi(2))*M
! reac n. 16: H + O2 <=> O + OH
prodf(16)=f_kf(16,Tint,Tdiff)*(coi(8))*(coi(10))
prodb(16)=f_kb(16,Tint,Tdiff)*(coi(13))*(coi(1))
! reac n. 17: 2 H + M <=> H2 + M
M=coi(1)+coi(2)+coi(5)+coi(7)*2+coi(8)+coi(9)+coi(10)+coi(11)+coi(12)+coi(13)
prodf(17)=f_kf(17,Tint,Tdiff)*(coi(8)*coi(8))*M
prodb(17)=f_kb(17,Tint,Tdiff)*(coi(3))*M
! reac n. 18: 2 H + H2 <=> H2 + H2
M=coi(3)
prodf(18)=f_kf(18,Tint,Tdiff)*(coi(8)*coi(8))*M
prodb(18)=f_kb(18,Tint,Tdiff)*(coi(3))*M
! reac n. 19: 2 H + H2O <=> H2 + H2O
M=coi(6)
prodf(19)=f_kf(19,Tint,Tdiff)*(coi(8)*coi(8))*M
prodb(19)=f_kb(19,Tint,Tdiff)*(coi(3))*M
! reac n. 20: 2 H + CO2 <=> H2 + CO2
M=coi(4)
prodf(20)=f_kf(20,Tint,Tdiff)*(coi(8)*coi(8))*M
prodb(20)=f_kb(20,Tint,Tdiff)*(coi(3))*M
! reac n. 21: H + OH + M <=> H2O + M
M=coi(1)+coi(2)+coi(3)*0.73+coi(4)+coi(5)+coi(6)*3.65+coi(7)*2+coi(8)+coi(9)+coi(10)+coi(11)+coi(12)+coi(13)
prodf(21)=f_kf(21,Tint,Tdiff)*(coi(8))*(coi(1))*M
prodb(21)=f_kb(21,Tint,Tdiff)*(coi(6))*M
! reac n. 22: H + HO2 <=> H2O + O
prodf(22)=f_kf(22,Tint,Tdiff)*(coi(8))*(coi(2))
prodb(22)=f_kb(22,Tint,Tdiff)*(coi(6))*(coi(13))
! reac n. 23: H + HO2 <=> H2 + O2
prodf(23)=f_kf(23,Tint,Tdiff)*(coi(8))*(coi(2))
prodb(23)=f_kb(23,Tint,Tdiff)*(coi(3))*(coi(10))
! reac n. 24: H + HO2 <=> 2 OH
prodf(24)=f_kf(24,Tint,Tdiff)*(coi(8))*(coi(2))
prodb(24)=f_kb(24,Tint,Tdiff)*(coi(1)*coi(1))
! reac n. 25: CH3 + H (+M) <=> CH4 (+M)
M=coi(1)+coi(2)+coi(3)*2+coi(4)*2+coi(5)*1.5+coi(6)*6+coi(7)*2+coi(8)+coi(9)+coi(10)+coi(11)+coi(12)+coi(13)
k = f_k_troe(1,Tint,Tdiff,M)
prodf(25)=k(1)*(coi(9))*(coi(8))
prodb(25)=k(2)*(coi(7))
! reac n. 26: CH4 + H <=> CH3 + H2
prodf(26)=f_kf(25,Tint,Tdiff)*(coi(7))*(coi(8))
prodb(26)=f_kb(25,Tint,Tdiff)*(coi(9))*(coi(3))
! reac n. 27: H + HCO (+M) <=> CH2O (+M)
M=coi(1)+coi(2)+coi(3)*2+coi(4)*2+coi(5)*1.5+coi(6)*6+coi(7)*2+coi(8)+coi(9)+coi(10)+coi(11)+coi(12)+coi(13)
k = f_k_troe(2,Tint,Tdiff,M)
prodf(27)=k(1)*(coi(8))*(coi(12))
prodb(27)=k(2)*(coi(11))
! reac n. 28: H + HCO <=> CO + H2
prodf(28)=f_kf(26,Tint,Tdiff)*(coi(8))*(coi(12))
prodb(28)=f_kb(26,Tint,Tdiff)*(coi(5))*(coi(3))
! reac n. 29: CH2O + H <=> H2 + HCO
prodf(29)=f_kf(27,Tint,Tdiff)*(coi(11))*(coi(8))
prodb(29)=f_kb(27,Tint,Tdiff)*(coi(3))*(coi(12))
! reac n. 30: CO + H2 (+M) <=> CH2O (+M)
M=coi(1)+coi(2)+coi(3)*2+coi(4)*2+coi(5)*1.5+coi(6)*6+coi(7)*2+coi(8)+coi(9)+coi(10)+coi(11)+coi(12)+coi(13)
k = f_k_troe(3,Tint,Tdiff,M)
prodf(30)=k(1)*(coi(5))*(coi(3))
prodb(30)=k(2)*(coi(11))
! reac n. 31: H2 + OH <=> H + H2O
prodf(31)=f_kf(28,Tint,Tdiff)*(coi(3))*(coi(1))
prodb(31)=f_kb(28,Tint,Tdiff)*(coi(8))*(coi(6))
! reac n. 32: 2 OH <=> H2O + O
prodf(32)=f_kf(29,Tint,Tdiff)*(coi(1)*coi(1))
prodb(32)=f_kb(29,Tint,Tdiff)*(coi(6))*(coi(13))
! reac n. 33: HO2 + OH <=> H2O + O2
prodf(33)=f_kf(30,Tint,Tdiff)*(coi(2))*(coi(1))
prodb(33)=f_kb(30,Tint,Tdiff)*(coi(6))*(coi(10))
! reac n. 34: CH4 + OH <=> CH3 + H2O
prodf(34)=f_kf(31,Tint,Tdiff)*(coi(7))*(coi(1))
prodb(34)=f_kb(31,Tint,Tdiff)*(coi(9))*(coi(6))
! reac n. 35: CO + OH <=> CO2 + H
prodf(35)=f_kf(32,Tint,Tdiff)*(coi(5))*(coi(1))
prodb(35)=f_kb(32,Tint,Tdiff)*(coi(4))*(coi(8))
! reac n. 36: HCO + OH <=> CO + H2O
prodf(36)=f_kf(33,Tint,Tdiff)*(coi(12))*(coi(1))
prodb(36)=f_kb(33,Tint,Tdiff)*(coi(5))*(coi(6))
! reac n. 37: CH2O + OH <=> H2O + HCO
prodf(37)=f_kf(34,Tint,Tdiff)*(coi(11))*(coi(1))
prodb(37)=f_kb(34,Tint,Tdiff)*(coi(6))*(coi(12))
! reac n. 38: CH3 + HO2 <=> CH4 + O2
prodf(38)=f_kf(35,Tint,Tdiff)*(coi(9))*(coi(2))
prodb(38)=f_kb(35,Tint,Tdiff)*(coi(7))*(coi(10))
! reac n. 39: CO + HO2 <=> CO2 + OH
prodf(39)=f_kf(36,Tint,Tdiff)*(coi(5))*(coi(2))
prodb(39)=f_kb(36,Tint,Tdiff)*(coi(4))*(coi(1))
! reac n. 40: CH3 + O2 <=> CH2O + OH
prodf(40)=f_kf(37,Tint,Tdiff)*(coi(9))*(coi(10))
prodb(40)=f_kb(37,Tint,Tdiff)*(coi(11))*(coi(1))
! reac n. 41: CH3 + HCO <=> CH4 + CO
prodf(41)=f_kf(38,Tint,Tdiff)*(coi(9))*(coi(12))
prodb(41)=f_kb(38,Tint,Tdiff)*(coi(7))*(coi(5))
! reac n. 42: CH2O + CH3 <=> CH4 + HCO
prodf(42)=f_kf(39,Tint,Tdiff)*(coi(11))*(coi(9))
prodb(42)=f_kb(39,Tint,Tdiff)*(coi(7))*(coi(12))
! reac n. 43: HCO + H2O <=> CO + H + H2O
M=coi(6)
prodf(43)=f_kf(40,Tint,Tdiff)*(coi(12))*M
prodb(43)=f_kb(40,Tint,Tdiff)*(coi(5))*(coi(8))*M
! reac n. 44: HCO + M <=> CO + H + M
M=coi(1)+coi(2)+coi(3)*2+coi(4)*2+coi(5)*1.5+coi(7)*2+coi(8)+coi(9)+coi(10)+coi(11)+coi(12)+coi(13)
prodf(44)=f_kf(41,Tint,Tdiff)*(coi(12))*M
prodb(44)=f_kb(41,Tint,Tdiff)*(coi(5))*(coi(8))*M
! reac n. 45: HCO + O2 <=> CO + HO2
prodf(45)=f_kf(42,Tint,Tdiff)*(coi(12))*(coi(10))
prodb(45)=f_kb(42,Tint,Tdiff)*(coi(5))*(coi(2))
! reac n. 46: CH3 + OH <=> CH2O + H2
prodf(46)=f_kf(43,Tint,Tdiff)*(coi(9))*(coi(1))
prodb(46)=f_kb(43,Tint,Tdiff)*(coi(11))*(coi(3))
! species source terms
omegadot(1)=Wm_tab(1)*(+(1.0-0.0)*(prodf(2)-prodb(2))+(1.0-0.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(8)-prodb(8))+(1.0-0.0)*(prodf(10)-prodb(10))+(1.0-0.0)*(prodf(16)-prodb(16))+(0.0-1.0)*(prodf(21)-prodb(21))+(2.0-0.0)*(prodf(24)-prodb(24))+(0.0-1.0)*(prodf(31)-prodb(31))+(0.0-2.0)*(prodf(32)-prodb(32))+(0.0-1.0)*(prodf(33)-prodb(33))+(0.0-1.0)*(prodf(34)-prodb(34))+(0.0-1.0)*(prodf(35)-prodb(35))+(0.0-1.0)*(prodf(36)-prodb(36))+(0.0-1.0)*(prodf(37)-prodb(37))+(1.0-0.0)*(prodf(39)-prodb(39))+(1.0-0.0)*(prodf(40)-prodb(40))+(0.0-1.0)*(prodf(46)-prodb(46)))
omegadot(2)=Wm_tab(2)*(+(0.0-1.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(12)-prodb(12))+(1.0-0.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(14)-prodb(14))+(1.0-0.0)*(prodf(15)-prodb(15))+(0.0-1.0)*(prodf(22)-prodb(22))+(0.0-1.0)*(prodf(23)-prodb(23))+(0.0-1.0)*(prodf(24)-prodb(24))+(0.0-1.0)*(prodf(33)-prodb(33))+(0.0-1.0)*(prodf(38)-prodb(38))+(0.0-1.0)*(prodf(39)-prodb(39))+(1.0-0.0)*(prodf(45)-prodb(45)))
omegadot(3)=Wm_tab(3)*(+(0.0-1.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(17)-prodb(17))+(1.0-0.0)*(prodf(18)-prodb(18))+(1.0-0.0)*(prodf(19)-prodb(19))+(1.0-0.0)*(prodf(20)-prodb(20))+(1.0-0.0)*(prodf(23)-prodb(23))+(1.0-0.0)*(prodf(26)-prodb(26))+(1.0-0.0)*(prodf(28)-prodb(28))+(1.0-0.0)*(prodf(29)-prodb(29))+(0.0-1.0)*(prodf(30)-prodb(30))+(0.0-1.0)*(prodf(31)-prodb(31))+(1.0-0.0)*(prodf(46)-prodb(46)))
omegadot(4)=Wm_tab(4)*(+(1.0-0.0)*(prodf(7)-prodb(7))+(1.0-0.0)*(prodf(9)-prodb(9))+(1.0-0.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(35)-prodb(35))+(1.0-0.0)*(prodf(39)-prodb(39)))
omegadot(5)=Wm_tab(5)*(+(0.0-1.0)*(prodf(7)-prodb(7))+(1.0-0.0)*(prodf(8)-prodb(8))+(0.0-1.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(28)-prodb(28))+(0.0-1.0)*(prodf(30)-prodb(30))+(0.0-1.0)*(prodf(35)-prodb(35))+(1.0-0.0)*(prodf(36)-prodb(36))+(0.0-1.0)*(prodf(39)-prodb(39))+(1.0-0.0)*(prodf(41)-prodb(41))+(1.0-0.0)*(prodf(43)-prodb(43))+(1.0-0.0)*(prodf(44)-prodb(44))+(1.0-0.0)*(prodf(45)-prodb(45)))
omegadot(6)=Wm_tab(6)*(+(1.0-0.0)*(prodf(21)-prodb(21))+(1.0-0.0)*(prodf(22)-prodb(22))+(1.0-0.0)*(prodf(31)-prodb(31))+(1.0-0.0)*(prodf(32)-prodb(32))+(1.0-0.0)*(prodf(33)-prodb(33))+(1.0-0.0)*(prodf(34)-prodb(34))+(1.0-0.0)*(prodf(36)-prodb(36))+(1.0-0.0)*(prodf(37)-prodb(37)))
omegadot(7)=Wm_tab(7)*(+(0.0-1.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(25)-prodb(25))+(0.0-1.0)*(prodf(26)-prodb(26))+(0.0-1.0)*(prodf(34)-prodb(34))+(1.0-0.0)*(prodf(38)-prodb(38))+(1.0-0.0)*(prodf(41)-prodb(41))+(1.0-0.0)*(prodf(42)-prodb(42)))
omegadot(8)=Wm_tab(8)*(+(0.0-1.0)*(prodf(2)-prodb(2))+(1.0-0.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(5)-prodb(5))+(1.0-0.0)*(prodf(9)-prodb(9))+(0.0-1.0)*(prodf(13)-prodb(13))+(0.0-1.0)*(prodf(14)-prodb(14))+(0.0-1.0)*(prodf(15)-prodb(15))+(0.0-1.0)*(prodf(16)-prodb(16))+(0.0-2.0)*(prodf(17)-prodb(17))+(0.0-2.0)*(prodf(18)-prodb(18))+(0.0-2.0)*(prodf(19)-prodb(19))+(0.0-2.0)*(prodf(20)-prodb(20))+(0.0-1.0)*(prodf(21)-prodb(21))+(0.0-1.0)*(prodf(22)-prodb(22))+(0.0-1.0)*(prodf(23)-prodb(23))+(0.0-1.0)*(prodf(24)-prodb(24))+(0.0-1.0)*(prodf(25)-prodb(25))+(0.0-1.0)*(prodf(26)-prodb(26))+(0.0-1.0)*(prodf(27)-prodb(27))+(0.0-1.0)*(prodf(28)-prodb(28))+(0.0-1.0)*(prodf(29)-prodb(29))+(1.0-0.0)*(prodf(31)-prodb(31))+(1.0-0.0)*(prodf(35)-prodb(35))+(1.0-0.0)*(prodf(43)-prodb(43))+(1.0-0.0)*(prodf(44)-prodb(44)))
omegadot(9)=Wm_tab(9)*(+(0.0-1.0)*(prodf(5)-prodb(5))+(1.0-0.0)*(prodf(6)-prodb(6))+(0.0-1.0)*(prodf(25)-prodb(25))+(1.0-0.0)*(prodf(26)-prodb(26))+(1.0-0.0)*(prodf(34)-prodb(34))+(0.0-1.0)*(prodf(38)-prodb(38))+(0.0-1.0)*(prodf(40)-prodb(40))+(0.0-1.0)*(prodf(41)-prodb(41))+(0.0-1.0)*(prodf(42)-prodb(42))+(0.0-1.0)*(prodf(46)-prodb(46)))
omegadot(10)=Wm_tab(10)*(+(1.0-0.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(4)-prodb(4))+(0.0-1.0)*(prodf(11)-prodb(11))+(0.0-1.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(13)-prodb(13))+(0.0-1.0)*(prodf(14)-prodb(14))+(0.0-1.0)*(prodf(15)-prodb(15))+(0.0-1.0)*(prodf(16)-prodb(16))+(1.0-0.0)*(prodf(23)-prodb(23))+(1.0-0.0)*(prodf(33)-prodb(33))+(1.0-0.0)*(prodf(38)-prodb(38))+(0.0-1.0)*(prodf(40)-prodb(40))+(0.0-1.0)*(prodf(45)-prodb(45)))
omegadot(11)=Wm_tab(11)*(+(1.0-0.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(10)-prodb(10))+(0.0-1.0)*(prodf(12)-prodb(12))+(1.0-0.0)*(prodf(27)-prodb(27))+(0.0-1.0)*(prodf(29)-prodb(29))+(1.0-0.0)*(prodf(30)-prodb(30))+(0.0-1.0)*(prodf(37)-prodb(37))+(1.0-0.0)*(prodf(40)-prodb(40))+(0.0-1.0)*(prodf(42)-prodb(42))+(1.0-0.0)*(prodf(46)-prodb(46)))
omegadot(12)=Wm_tab(12)*(+(0.0-1.0)*(prodf(8)-prodb(8))+(0.0-1.0)*(prodf(9)-prodb(9))+(1.0-0.0)*(prodf(10)-prodb(10))+(1.0-0.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(27)-prodb(27))+(0.0-1.0)*(prodf(28)-prodb(28))+(1.0-0.0)*(prodf(29)-prodb(29))+(0.0-1.0)*(prodf(36)-prodb(36))+(1.0-0.0)*(prodf(37)-prodb(37))+(0.0-1.0)*(prodf(41)-prodb(41))+(1.0-0.0)*(prodf(42)-prodb(42))+(0.0-1.0)*(prodf(43)-prodb(43))+(0.0-1.0)*(prodf(44)-prodb(44))+(0.0-1.0)*(prodf(45)-prodb(45)))
omegadot(13)=Wm_tab(13)*(+(0.0-2.0)*(prodf(1)-prodb(1))+(0.0-1.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(4)-prodb(4))+(0.0-1.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(6)-prodb(6))+(0.0-1.0)*(prodf(7)-prodb(7))+(0.0-1.0)*(prodf(8)-prodb(8))+(0.0-1.0)*(prodf(9)-prodb(9))+(0.0-1.0)*(prodf(10)-prodb(10))+(1.0-0.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(16)-prodb(16))+(1.0-0.0)*(prodf(22)-prodb(22))+(1.0-0.0)*(prodf(32)-prodb(32)))
end subroutine TSRCDF13
end module TSRCDF13_mod
