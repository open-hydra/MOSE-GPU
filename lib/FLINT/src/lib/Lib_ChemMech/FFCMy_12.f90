module FFCMy_12_mod
implicit none
contains
subroutine FFCMy_12(roi,temp,omegadot)
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
real(8) :: prodf(1:38), prodb(1:38)
real(8) :: k(2) !< Falloff rate coefficients


do is = 1, ns 
 coi(is)=roi(is)/Wm_tab(is) ! kmol/m^3
enddo 
T_i = int(temp) 
Tdiff  = temp-T_i 
Tint(1) = T_i 
Tint(2) = T_i + 1 
! reac n. 1: H + O2 <=> O + OH
prodf(1)=f_kf(1,Tint,Tdiff)*(coi(2))*(coi(3))
prodb(1)=f_kb(1,Tint,Tdiff)*(coi(4))*(coi(5))
! reac n. 2: H2 + O <=> H + OH
prodf(2)=f_kf(2,Tint,Tdiff)*(coi(1))*(coi(4))
prodb(2)=f_kb(2,Tint,Tdiff)*(coi(2))*(coi(5))
! reac n. 3: H2 + O <=> H + OH
prodf(3)=f_kf(3,Tint,Tdiff)*(coi(1))*(coi(4))
prodb(3)=f_kb(3,Tint,Tdiff)*(coi(2))*(coi(5))
! reac n. 4: H2 + OH <=> H + H2O
prodf(4)=f_kf(4,Tint,Tdiff)*(coi(1))*(coi(5))
prodb(4)=f_kb(4,Tint,Tdiff)*(coi(2))*(coi(7))
! reac n. 5: 2 OH <=> H2O + O
prodf(5)=f_kf(5,Tint,Tdiff)*(coi(5)*coi(5))
prodb(5)=f_kb(5,Tint,Tdiff)*(coi(7))*(coi(4))
! reac n. 6: H2 + M <=> 2 H + M
M=coi(1)*2.5+coi(2)+coi(3)+coi(4)+coi(5)+coi(6)+coi(7)*12+coi(8)+coi(9)*2+coi(10)*1.9+coi(11)*3.8+coi(12)*2.5+coi(13)
prodf(6)=f_kf(6,Tint,Tdiff)*(coi(1))*M
prodb(6)=f_kb(6,Tint,Tdiff)*(coi(2)*coi(2))*M
! reac n. 7: H + O + M <=> OH + M
M=coi(1)*2.5+coi(2)+coi(3)+coi(4)+coi(5)+coi(6)+coi(7)*12+coi(8)+coi(9)*2+coi(10)*1.9+coi(11)*3.8+coi(12)*2.5+coi(13)
prodf(7)=f_kf(7,Tint,Tdiff)*(coi(2))*(coi(4))*M
prodb(7)=f_kb(7,Tint,Tdiff)*(coi(5))*M
! reac n. 8: H2O + M <=> H + OH + M
M=coi(1)*3+coi(2)+coi(3)*1.5+coi(4)+coi(5)+coi(6)+coi(8)+coi(9)*7+coi(10)*1.9+coi(11)*3.8+coi(12)*2.5+coi(13)*2
prodf(8)=f_kf(8,Tint,Tdiff)*(coi(7))*M
prodb(8)=f_kb(8,Tint,Tdiff)*(coi(2))*(coi(5))*M
! reac n. 9: H2O + H2O <=> H + OH + H2O
M=coi(7)
prodf(9)=f_kf(9,Tint,Tdiff)*(coi(7))*M
prodb(9)=f_kb(9,Tint,Tdiff)*(coi(2))*(coi(5))*M
! reac n. 10: H + O2 (+M) <=> HO2 (+M)
M=coi(1)*2+coi(2)+coi(3)*0.78+coi(4)+coi(5)+coi(6)+coi(7)*14+coi(8)+coi(9)*2+coi(10)*1.9+coi(11)*3.8+coi(12)*2.5+coi(13)
k = f_k_troe(1,Tint,Tdiff,M)
prodf(10)=k(1)*(coi(2))*(coi(3))
prodb(10)=k(2)*(coi(6))
! reac n. 11: H + HO2 <=> H2 + O2
prodf(11)=f_kf(10,Tint,Tdiff)*(coi(2))*(coi(6))
prodb(11)=f_kb(10,Tint,Tdiff)*(coi(1))*(coi(3))
! reac n. 12: H + HO2 <=> 2 OH
prodf(12)=f_kf(11,Tint,Tdiff)*(coi(2))*(coi(6))
prodb(12)=f_kb(11,Tint,Tdiff)*(coi(5)*coi(5))
! reac n. 13: H + HO2 <=> H2O + O
prodf(13)=f_kf(12,Tint,Tdiff)*(coi(2))*(coi(6))
prodb(13)=f_kb(12,Tint,Tdiff)*(coi(7))*(coi(4))
! reac n. 14: HO2 + O <=> O2 + OH
prodf(14)=f_kf(13,Tint,Tdiff)*(coi(6))*(coi(4))
prodb(14)=f_kb(13,Tint,Tdiff)*(coi(3))*(coi(5))
! reac n. 15: HO2 + OH <=> H2O + O2
prodf(15)=f_kf(14,Tint,Tdiff)*(coi(6))*(coi(5))
prodb(15)=f_kb(14,Tint,Tdiff)*(coi(7))*(coi(3))
! reac n. 16: HO2 + OH <=> H2O + O2
prodf(16)=f_kf(15,Tint,Tdiff)*(coi(6))*(coi(5))
prodb(16)=f_kb(15,Tint,Tdiff)*(coi(7))*(coi(3))
! reac n. 17: CO + O (+M) <=> CO2 (+M)
M=coi(1)*2.5+coi(2)+coi(3)+coi(4)+coi(5)+coi(6)+coi(7)*12+coi(8)+coi(9)*2+coi(10)*1.9+coi(11)*3.8+coi(12)*2.5+coi(13)
k = f_k_lindemann(1,Tint,Tdiff,M)
prodf(17)=k(1)*(coi(10))*(coi(4))
prodb(17)=k(2)*(coi(11))
! reac n. 18: CO + O2 <=> CO2 + O
prodf(18)=f_kf(16,Tint,Tdiff)*(coi(10))*(coi(3))
prodb(18)=f_kb(16,Tint,Tdiff)*(coi(11))*(coi(4))
! reac n. 19: CO + OH <=> CO2 + H
prodf(19)=f_kf(17,Tint,Tdiff)*(coi(10))*(coi(5))
prodb(19)=f_kb(17,Tint,Tdiff)*(coi(11))*(coi(2))
! reac n. 20: CO + OH <=> CO2 + H
prodf(20)=f_kf(18,Tint,Tdiff)*(coi(10))*(coi(5))
prodb(20)=f_kb(18,Tint,Tdiff)*(coi(11))*(coi(2))
! reac n. 21: CO + HO2 <=> CO2 + OH
prodf(21)=f_kf(19,Tint,Tdiff)*(coi(10))*(coi(6))
prodb(21)=f_kb(19,Tint,Tdiff)*(coi(11))*(coi(5))
! reac n. 22: CH4 + H <=> CH3 + H2
prodf(22)=f_kf(20,Tint,Tdiff)*(coi(9))*(coi(2))
prodb(22)=f_kb(20,Tint,Tdiff)*(coi(8))*(coi(1))
! reac n. 23: CH4 + O <=> CH3 + OH
prodf(23)=f_kf(21,Tint,Tdiff)*(coi(9))*(coi(4))
prodb(23)=f_kb(21,Tint,Tdiff)*(coi(8))*(coi(5))
! reac n. 24: CH4 + OH <=> CH3 + H2O
prodf(24)=f_kf(22,Tint,Tdiff)*(coi(9))*(coi(5))
prodb(24)=f_kb(22,Tint,Tdiff)*(coi(8))*(coi(7))
! reac n. 25: CH3 + H (+M) <=> CH4 (+M)
M=coi(1)*2+coi(2)+coi(3)+coi(4)+coi(5)+coi(6)+coi(7)*6+coi(8)+coi(9)*2+coi(10)*1.5+coi(11)*2+coi(12)*2.5+coi(13)
k = f_k_troe(2,Tint,Tdiff,M)
prodf(25)=k(1)*(coi(8))*(coi(2))
prodb(25)=k(2)*(coi(9))
! reac n. 26: CH3 + O <=> CH2O + H
prodf(26)=f_kf(23,Tint,Tdiff)*(coi(8))*(coi(4))
prodb(26)=f_kb(23,Tint,Tdiff)*(coi(12))*(coi(2))
! reac n. 27: CH3 + O => CO + H + H2
prodf(27)=f_kf(24,Tint,Tdiff)*(coi(8))*(coi(4))
prodb(27)=f_kb(24,Tint,Tdiff)*(coi(10))*(coi(2))*(coi(1))
! reac n. 28: CH3 + HO2 <=> CH4 + O2
prodf(28)=f_kf(25,Tint,Tdiff)*(coi(8))*(coi(6))
prodb(28)=f_kb(25,Tint,Tdiff)*(coi(9))*(coi(3))
! reac n. 29: CH3 + HO2 => CH2O + H + OH
prodf(29)=f_kf(26,Tint,Tdiff)*(coi(8))*(coi(6))
prodb(29)=f_kb(26,Tint,Tdiff)*(coi(12))*(coi(2))*(coi(5))
! reac n. 30: CH3 + O2 => CH2O + H + O
prodf(30)=f_kf(27,Tint,Tdiff)*(coi(8))*(coi(3))
prodb(30)=f_kb(27,Tint,Tdiff)*(coi(12))*(coi(2))*(coi(4))
! reac n. 31: CH3 + O2 <=> CH2O + OH
prodf(31)=f_kf(28,Tint,Tdiff)*(coi(8))*(coi(3))
prodb(31)=f_kb(28,Tint,Tdiff)*(coi(12))*(coi(5))
! reac n. 32: CH2O + CH3 => CH4 + CO + H
prodf(32)=f_kf(29,Tint,Tdiff)*(coi(12))*(coi(8))
prodb(32)=f_kb(29,Tint,Tdiff)*(coi(9))*(coi(10))*(coi(2))
! reac n. 33: CH2O (+M) <=> CO + H2 (+M)
M=coi(1)*2+coi(2)+coi(3)+coi(4)+coi(5)+coi(6)+coi(7)*6+coi(8)+coi(9)*2+coi(10)*1.5+coi(11)*2+coi(12)*2.5+coi(13)
k = f_k_troe(3,Tint,Tdiff,M)
prodf(33)=k(1)*(coi(12))
prodb(33)=k(2)*(coi(10))*(coi(1))
! reac n. 34: CH2O + H => CO + H2 + H
M=coi(2)
prodf(34)=f_kf(30,Tint,Tdiff)*(coi(12))*M
prodb(34)=f_kb(30,Tint,Tdiff)*(coi(10))*(coi(1))*M
! reac n. 35: CH2O + H => CO + H2 + H
M=coi(2)
prodf(35)=f_kf(31,Tint,Tdiff)*(coi(12))*M
prodb(35)=f_kb(31,Tint,Tdiff)*(coi(10))*(coi(1))*M
! reac n. 36: CH2O + O => CO + H + OH
prodf(36)=f_kf(32,Tint,Tdiff)*(coi(12))*(coi(4))
prodb(36)=f_kb(32,Tint,Tdiff)*(coi(10))*(coi(2))*(coi(5))
! reac n. 37: CH2O + OH => CO + H + H2O
prodf(37)=f_kf(33,Tint,Tdiff)*(coi(12))*(coi(5))
prodb(37)=f_kb(33,Tint,Tdiff)*(coi(10))*(coi(2))*(coi(7))
! reac n. 38: CH2O + O2 => CO + H + HO2
prodf(38)=f_kf(34,Tint,Tdiff)*(coi(12))*(coi(3))
prodb(38)=f_kb(34,Tint,Tdiff)*(coi(10))*(coi(2))*(coi(6))
! species source terms
omegadot(1)=Wm_tab(1)*(+(0.0-1.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(4)-prodb(4))+(0.0-1.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(22)-prodb(22))+(1.0-0.0)*(prodf(27)-prodb(27))+(1.0-0.0)*(prodf(33)-prodb(33))+(1.0-0.0)*(prodf(34)-prodb(34))+(1.0-0.0)*(prodf(35)-prodb(35)))
omegadot(2)=Wm_tab(2)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(2)-prodb(2))+(1.0-0.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(4)-prodb(4))+(2.0-0.0)*(prodf(6)-prodb(6))+(0.0-1.0)*(prodf(7)-prodb(7))+(1.0-0.0)*(prodf(8)-prodb(8))+(1.0-0.0)*(prodf(9)-prodb(9))+(0.0-1.0)*(prodf(10)-prodb(10))+(0.0-1.0)*(prodf(11)-prodb(11))+(0.0-1.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(19)-prodb(19))+(1.0-0.0)*(prodf(20)-prodb(20))+(0.0-1.0)*(prodf(22)-prodb(22))+(0.0-1.0)*(prodf(25)-prodb(25))+(1.0-0.0)*(prodf(26)-prodb(26))+(1.0-0.0)*(prodf(27)-prodb(27))+(1.0-0.0)*(prodf(29)-prodb(29))+(1.0-0.0)*(prodf(30)-prodb(30))+(1.0-0.0)*(prodf(32)-prodb(32))+(1.0-0.0)*(prodf(36)-prodb(36))+(1.0-0.0)*(prodf(37)-prodb(37))+(1.0-0.0)*(prodf(38)-prodb(38)))
omegadot(3)=Wm_tab(3)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(0.0-1.0)*(prodf(10)-prodb(10))+(1.0-0.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(14)-prodb(14))+(1.0-0.0)*(prodf(15)-prodb(15))+(1.0-0.0)*(prodf(16)-prodb(16))+(0.0-1.0)*(prodf(18)-prodb(18))+(1.0-0.0)*(prodf(28)-prodb(28))+(0.0-1.0)*(prodf(30)-prodb(30))+(0.0-1.0)*(prodf(31)-prodb(31))+(0.0-1.0)*(prodf(38)-prodb(38)))
omegadot(4)=Wm_tab(4)*(+(1.0-0.0)*(prodf(1)-prodb(1))+(0.0-1.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(7)-prodb(7))+(1.0-0.0)*(prodf(13)-prodb(13))+(0.0-1.0)*(prodf(14)-prodb(14))+(0.0-1.0)*(prodf(17)-prodb(17))+(1.0-0.0)*(prodf(18)-prodb(18))+(0.0-1.0)*(prodf(23)-prodb(23))+(0.0-1.0)*(prodf(26)-prodb(26))+(0.0-1.0)*(prodf(27)-prodb(27))+(1.0-0.0)*(prodf(30)-prodb(30))+(0.0-1.0)*(prodf(36)-prodb(36)))
omegadot(5)=Wm_tab(5)*(+(1.0-0.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(2)-prodb(2))+(1.0-0.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(4)-prodb(4))+(0.0-2.0)*(prodf(5)-prodb(5))+(1.0-0.0)*(prodf(7)-prodb(7))+(1.0-0.0)*(prodf(8)-prodb(8))+(1.0-0.0)*(prodf(9)-prodb(9))+(2.0-0.0)*(prodf(12)-prodb(12))+(1.0-0.0)*(prodf(14)-prodb(14))+(0.0-1.0)*(prodf(15)-prodb(15))+(0.0-1.0)*(prodf(16)-prodb(16))+(0.0-1.0)*(prodf(19)-prodb(19))+(0.0-1.0)*(prodf(20)-prodb(20))+(1.0-0.0)*(prodf(21)-prodb(21))+(1.0-0.0)*(prodf(23)-prodb(23))+(0.0-1.0)*(prodf(24)-prodb(24))+(1.0-0.0)*(prodf(29)-prodb(29))+(1.0-0.0)*(prodf(31)-prodb(31))+(1.0-0.0)*(prodf(36)-prodb(36))+(0.0-1.0)*(prodf(37)-prodb(37)))
omegadot(6)=Wm_tab(6)*(+(1.0-0.0)*(prodf(10)-prodb(10))+(0.0-1.0)*(prodf(11)-prodb(11))+(0.0-1.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(13)-prodb(13))+(0.0-1.0)*(prodf(14)-prodb(14))+(0.0-1.0)*(prodf(15)-prodb(15))+(0.0-1.0)*(prodf(16)-prodb(16))+(0.0-1.0)*(prodf(21)-prodb(21))+(0.0-1.0)*(prodf(28)-prodb(28))+(0.0-1.0)*(prodf(29)-prodb(29))+(1.0-0.0)*(prodf(38)-prodb(38)))
omegadot(7)=Wm_tab(7)*(+(1.0-0.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(8)-prodb(8))+(0.0-1.0)*(prodf(9)-prodb(9))+(1.0-0.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(15)-prodb(15))+(1.0-0.0)*(prodf(16)-prodb(16))+(1.0-0.0)*(prodf(24)-prodb(24))+(1.0-0.0)*(prodf(37)-prodb(37)))
omegadot(8)=Wm_tab(8)*(+(1.0-0.0)*(prodf(22)-prodb(22))+(1.0-0.0)*(prodf(23)-prodb(23))+(1.0-0.0)*(prodf(24)-prodb(24))+(0.0-1.0)*(prodf(25)-prodb(25))+(0.0-1.0)*(prodf(26)-prodb(26))+(0.0-1.0)*(prodf(27)-prodb(27))+(0.0-1.0)*(prodf(28)-prodb(28))+(0.0-1.0)*(prodf(29)-prodb(29))+(0.0-1.0)*(prodf(30)-prodb(30))+(0.0-1.0)*(prodf(31)-prodb(31))+(0.0-1.0)*(prodf(32)-prodb(32)))
omegadot(9)=Wm_tab(9)*(+(0.0-1.0)*(prodf(22)-prodb(22))+(0.0-1.0)*(prodf(23)-prodb(23))+(0.0-1.0)*(prodf(24)-prodb(24))+(1.0-0.0)*(prodf(25)-prodb(25))+(1.0-0.0)*(prodf(28)-prodb(28))+(1.0-0.0)*(prodf(32)-prodb(32)))
omegadot(10)=Wm_tab(10)*(+(0.0-1.0)*(prodf(17)-prodb(17))+(0.0-1.0)*(prodf(18)-prodb(18))+(0.0-1.0)*(prodf(19)-prodb(19))+(0.0-1.0)*(prodf(20)-prodb(20))+(0.0-1.0)*(prodf(21)-prodb(21))+(1.0-0.0)*(prodf(27)-prodb(27))+(1.0-0.0)*(prodf(32)-prodb(32))+(1.0-0.0)*(prodf(33)-prodb(33))+(1.0-0.0)*(prodf(34)-prodb(34))+(1.0-0.0)*(prodf(35)-prodb(35))+(1.0-0.0)*(prodf(36)-prodb(36))+(1.0-0.0)*(prodf(37)-prodb(37))+(1.0-0.0)*(prodf(38)-prodb(38)))
omegadot(11)=Wm_tab(11)*(+(1.0-0.0)*(prodf(17)-prodb(17))+(1.0-0.0)*(prodf(18)-prodb(18))+(1.0-0.0)*(prodf(19)-prodb(19))+(1.0-0.0)*(prodf(20)-prodb(20))+(1.0-0.0)*(prodf(21)-prodb(21)))
omegadot(12)=Wm_tab(12)*(+(1.0-0.0)*(prodf(26)-prodb(26))+(1.0-0.0)*(prodf(29)-prodb(29))+(1.0-0.0)*(prodf(30)-prodb(30))+(1.0-0.0)*(prodf(31)-prodb(31))+(0.0-1.0)*(prodf(32)-prodb(32))+(0.0-1.0)*(prodf(33)-prodb(33))+(0.0-1.0)*(prodf(34)-prodb(34))+(0.0-1.0)*(prodf(35)-prodb(35))+(0.0-1.0)*(prodf(36)-prodb(36))+(0.0-1.0)*(prodf(37)-prodb(37))+(0.0-1.0)*(prodf(38)-prodb(38)))
omegadot(13)=0.d0
end subroutine FFCMy_12
end module FFCMy_12_mod
