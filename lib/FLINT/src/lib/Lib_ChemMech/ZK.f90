module ZK_mod
implicit none
contains
subroutine ZK(roi,temp,omegadot)
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
real(8) :: prodf(1:51), prodb(1:51)
real(8) :: k(2) !< Troe rate coefficients


do is = 1, ns 
 coi(is)=roi(is)/Wm_tab(is) ! kmol/m^3
enddo 
T_i = int(temp) 
Tdiff  = temp-T_i 
Tint(1) = T_i 
Tint(2) = T_i + 1 
! reac n. 1: CH2O + O2 <=> HCO + HO2
prodf(1)=f_kf(1,Tint,Tdiff)*(coi(19))*(coi(6))
prodb(1)=f_kb(1,Tint,Tdiff)*(coi(18))*(coi(9))
! reac n. 2: H + O2 + M <=> HO2 + M
M=coi(3)+coi(4)+coi(5)+coi(7)+coi(8)*6+coi(9)+coi(10)+coi(11)*0.75+coi(12)*1.5+coi(13)+coi(14)+coi(15)+coi(16)+coi(17)+coi(18)+coi(19)+coi(20)+coi(21)+coi(22)+coi(23)+coi(24)*1.5+coi(25)
prodf(2)=f_kf(2,Tint,Tdiff)*(coi(4))*(coi(6))*M
prodb(2)=f_kb(2,Tint,Tdiff)*(coi(9))*M
! reac n. 3: H + O2 + O2 <=> HO2 + O2
M=coi(6)
prodf(3)=f_kf(3,Tint,Tdiff)*(coi(4))*(coi(6))*M
prodb(3)=f_kb(3,Tint,Tdiff)*(coi(9))*M
! reac n. 4: CH2O + H (+M) <=> CH3O (+M)
M=coi(1)+coi(2)+coi(3)*2+coi(4)+coi(5)+coi(6)+coi(7)+coi(8)*6+coi(9)+coi(10)+coi(11)*1.5+coi(12)*2+coi(13)+coi(14)+coi(15)+coi(16)+coi(17)*2+coi(18)+coi(19)+coi(20)+coi(21)+coi(22)+coi(23)+coi(24)*3+coi(25)
k = f_k_troe(1,Tint,Tdiff,M)
prodf(4)=k(1)*(coi(19))*(coi(4))
prodb(4)=k(2)*(coi(20))
! reac n. 5: 2 OH (+M) <=> H2O2 (+M)
M=coi(1)*0.7+coi(2)+coi(3)*2+coi(4)+coi(5)+coi(6)+coi(7)+coi(8)*6+coi(9)+coi(10)+coi(11)*1.5+coi(12)*2+coi(13)+coi(14)+coi(15)+coi(16)+coi(17)*2+coi(18)+coi(19)+coi(20)+coi(21)+coi(22)+coi(23)+coi(24)*3+coi(25)
k = f_k_troe(2,Tint,Tdiff,M)
prodf(5)=k(1)*(coi(7)*coi(7))
prodb(5)=k(2)*(coi(10))
! reac n. 6: HO2 + OH <=> H2O + O2
prodf(6)=f_kf(4,Tint,Tdiff)*(coi(9))*(coi(7))
prodb(6)=f_kb(4,Tint,Tdiff)*(coi(8))*(coi(6))
! reac n. 7: H2O2 + OH <=> H2O + HO2
prodf(7)=f_kf(5,Tint,Tdiff)*(coi(10))*(coi(7))
prodb(7)=f_kb(5,Tint,Tdiff)*(coi(8))*(coi(9))
! reac n. 8: H2O2 + OH <=> H2O + HO2
prodf(8)=f_kf(6,Tint,Tdiff)*(coi(10))*(coi(7))
prodb(8)=f_kb(6,Tint,Tdiff)*(coi(8))*(coi(9))
! reac n. 9: CH4 + OH <=> CH3 + H2O
prodf(9)=f_kf(7,Tint,Tdiff)*(coi(17))*(coi(7))
prodb(9)=f_kb(7,Tint,Tdiff)*(coi(16))*(coi(8))
! reac n. 10: 2 HO2 <=> H2O2 + O2
prodf(10)=f_kf(8,Tint,Tdiff)*(coi(9)*coi(9))
prodb(10)=f_kb(8,Tint,Tdiff)*(coi(10))*(coi(6))
! reac n. 11: 2 HO2 <=> H2O2 + O2
prodf(11)=f_kf(9,Tint,Tdiff)*(coi(9)*coi(9))
prodb(11)=f_kb(9,Tint,Tdiff)*(coi(10))*(coi(6))
! reac n. 12: CH3 + HO2 <=> CH4 + O2
prodf(12)=f_kf(10,Tint,Tdiff)*(coi(16))*(coi(9))
prodb(12)=f_kb(10,Tint,Tdiff)*(coi(17))*(coi(6))
! reac n. 13: CH3 + HO2 <=> CH3O + OH
prodf(13)=f_kf(11,Tint,Tdiff)*(coi(16))*(coi(9))
prodb(13)=f_kb(11,Tint,Tdiff)*(coi(20))*(coi(7))
! reac n. 14: CO + HO2 <=> CO2 + OH
prodf(14)=f_kf(12,Tint,Tdiff)*(coi(11))*(coi(9))
prodb(14)=f_kb(12,Tint,Tdiff)*(coi(12))*(coi(7))
! reac n. 15: CH2O + HO2 <=> H2O2 + HCO
prodf(15)=f_kf(13,Tint,Tdiff)*(coi(19))*(coi(9))
prodb(15)=f_kb(13,Tint,Tdiff)*(coi(10))*(coi(18))
! reac n. 16: CH3 + O2 <=> CH3O + O
prodf(16)=f_kf(14,Tint,Tdiff)*(coi(16))*(coi(6))
prodb(16)=f_kb(14,Tint,Tdiff)*(coi(20))*(coi(5))
! reac n. 17: CH3 + O2 <=> CH2O + OH
prodf(17)=f_kf(15,Tint,Tdiff)*(coi(16))*(coi(6))
prodb(17)=f_kb(15,Tint,Tdiff)*(coi(19))*(coi(7))
! reac n. 18: CH3 + H2O2 <=> CH4 + HO2
prodf(18)=f_kf(16,Tint,Tdiff)*(coi(16))*(coi(10))
prodb(18)=f_kb(16,Tint,Tdiff)*(coi(17))*(coi(9))
! reac n. 19: CH2O + CH3 <=> CH4 + HCO
prodf(19)=f_kf(17,Tint,Tdiff)*(coi(19))*(coi(16))
prodb(19)=f_kb(17,Tint,Tdiff)*(coi(17))*(coi(18))
! reac n. 20: CH3O + HO2 <=> CH2O + H2O2
prodf(20)=f_kf(18,Tint,Tdiff)*(coi(20))*(coi(9))
prodb(20)=f_kb(18,Tint,Tdiff)*(coi(19))*(coi(10))
! reac n. 21: CH3 + CH3O2 <=> 2 CH3O
prodf(21)=f_kf(19,Tint,Tdiff)*(coi(16))*(coi(25))
prodb(21)=f_kb(19,Tint,Tdiff)*(coi(20)*coi(20))
! reac n. 22: CH3O + O2 <=> CH2O + HO2
prodf(22)=f_kf(20,Tint,Tdiff)*(coi(20))*(coi(6))
prodb(22)=f_kb(20,Tint,Tdiff)*(coi(19))*(coi(9))
! reac n. 23: CH3 + O2 <=> CH3O2
prodf(23)=f_kf(21,Tint,Tdiff)*(coi(16))*(coi(6))
prodb(23)=f_kb(21,Tint,Tdiff)*(coi(25))
! reac n. 24: CH3 + CH3O <=> CH2O + CH4
prodf(24)=f_kf(22,Tint,Tdiff)*(coi(16))*(coi(20))
prodb(24)=f_kb(22,Tint,Tdiff)*(coi(19))*(coi(17))
! reac n. 25: CH4 + O <=> CH3 + OH
prodf(25)=f_kf(23,Tint,Tdiff)*(coi(17))*(coi(5))
prodb(25)=f_kb(23,Tint,Tdiff)*(coi(16))*(coi(7))
! reac n. 26: H + O2 <=> O + OH
prodf(26)=f_kf(24,Tint,Tdiff)*(coi(4))*(coi(6))
prodb(26)=f_kb(24,Tint,Tdiff)*(coi(5))*(coi(7))
! reac n. 27: H + O2 + H2O <=> HO2 + H2O
M=coi(8)
prodf(27)=f_kf(25,Tint,Tdiff)*(coi(4))*(coi(6))*M
prodb(27)=f_kb(25,Tint,Tdiff)*(coi(9))*M
! reac n. 28: H2 + O <=> H + OH
prodf(28)=f_kf(26,Tint,Tdiff)*(coi(3))*(coi(5))
prodb(28)=f_kb(26,Tint,Tdiff)*(coi(4))*(coi(7))
! reac n. 29: CH3 + O <=> CH2O + H
prodf(29)=f_kf(27,Tint,Tdiff)*(coi(16))*(coi(5))
prodb(29)=f_kb(27,Tint,Tdiff)*(coi(19))*(coi(4))
! reac n. 30: CO + O + M <=> CO2 + M
M=coi(1)*0.5+coi(2)+coi(3)*2+coi(4)+coi(5)+coi(6)*6+coi(7)+coi(8)*6+coi(9)+coi(10)+coi(11)*1.5+coi(12)*3.5+coi(13)+coi(14)+coi(15)+coi(16)+coi(17)*2+coi(18)+coi(19)+coi(20)+coi(21)+coi(22)+coi(23)+coi(24)*3+coi(25)
prodf(30)=f_kf(28,Tint,Tdiff)*(coi(11))*(coi(5))*M
prodb(30)=f_kb(28,Tint,Tdiff)*(coi(12))*M
! reac n. 31: H + OH + M <=> H2O + M
M=coi(1)*0.38+coi(2)+coi(3)*0.73+coi(4)+coi(5)+coi(6)+coi(7)+coi(8)*3.65+coi(9)+coi(10)+coi(11)+coi(12)+coi(13)+coi(14)+coi(15)+coi(16)+coi(17)*2+coi(18)+coi(19)+coi(20)+coi(21)+coi(22)+coi(23)+coi(24)*3+coi(25)
prodf(31)=f_kf(29,Tint,Tdiff)*(coi(4))*(coi(7))*M
prodb(31)=f_kb(29,Tint,Tdiff)*(coi(8))*M
! reac n. 32: CH3 + H (+M) <=> CH4 (+M)
M=coi(1)+coi(2)+coi(3)*2+coi(4)+coi(5)+coi(6)+coi(7)+coi(8)*6+coi(9)+coi(10)+coi(11)*1.5+coi(12)*2+coi(13)+coi(14)+coi(15)+coi(16)+coi(17)*2+coi(18)+coi(19)+coi(20)+coi(21)+coi(22)+coi(23)+coi(24)*3+coi(25)
k = f_k_troe(3,Tint,Tdiff,M)
prodf(32)=k(1)*(coi(16))*(coi(4))
prodb(32)=k(2)*(coi(17))
! reac n. 33: H + HCO (+M) <=> CH2O (+M)
M=coi(1)+coi(2)+coi(3)*2+coi(4)+coi(5)+coi(6)+coi(7)+coi(8)*6+coi(9)+coi(10)+coi(11)*1.5+coi(12)*2+coi(13)+coi(14)+coi(15)+coi(16)+coi(17)*2+coi(18)+coi(19)+coi(20)+coi(21)+coi(22)+coi(23)+coi(24)*3+coi(25)
k = f_k_troe(4,Tint,Tdiff,M)
prodf(33)=k(1)*(coi(4))*(coi(18))
prodb(33)=k(2)*(coi(19))
! reac n. 34: C2H4 + H (+M) <=> C2H5 (+M)
M=coi(1)+coi(2)+coi(3)*2+coi(4)+coi(5)+coi(6)+coi(7)+coi(8)*6+coi(9)+coi(10)+coi(11)*1.5+coi(12)*2+coi(13)+coi(14)+coi(15)+coi(16)+coi(17)*2+coi(18)+coi(19)+coi(20)+coi(21)+coi(22)+coi(23)+coi(24)*3+coi(25)
k = f_k_troe(5,Tint,Tdiff,M)
prodf(34)=k(1)*(coi(22))*(coi(4))
prodb(34)=k(2)*(coi(23))
! reac n. 35: C2H4 + H <=> C2H3 + H2
prodf(35)=f_kf(30,Tint,Tdiff)*(coi(22))*(coi(4))
prodb(35)=f_kb(30,Tint,Tdiff)*(coi(21))*(coi(3))
! reac n. 36: C2H6 + H <=> C2H5 + H2
prodf(36)=f_kf(31,Tint,Tdiff)*(coi(24))*(coi(4))
prodb(36)=f_kb(31,Tint,Tdiff)*(coi(23))*(coi(3))
! reac n. 37: H2 + OH <=> H + H2O
prodf(37)=f_kf(32,Tint,Tdiff)*(coi(3))*(coi(7))
prodb(37)=f_kb(32,Tint,Tdiff)*(coi(4))*(coi(8))
! reac n. 38: CH2 + OH <=> CH2O + H
prodf(38)=f_kf(33,Tint,Tdiff)*(coi(15))*(coi(7))
prodb(38)=f_kb(33,Tint,Tdiff)*(coi(19))*(coi(4))
! reac n. 39: C2H6 + OH <=> C2H5 + H2O
prodf(39)=f_kf(34,Tint,Tdiff)*(coi(24))*(coi(7))
prodb(39)=f_kb(34,Tint,Tdiff)*(coi(23))*(coi(8))
! reac n. 40: HCO + O2 <=> CO + HO2
prodf(40)=f_kf(35,Tint,Tdiff)*(coi(18))*(coi(6))
prodb(40)=f_kb(35,Tint,Tdiff)*(coi(11))*(coi(9))
! reac n. 41: HCO + M <=> CO + H + M
M=coi(1)+coi(2)+coi(3)*2+coi(4)+coi(5)+coi(6)+coi(7)+coi(8)*12+coi(9)+coi(10)+coi(11)*1.5+coi(12)*2+coi(13)+coi(14)+coi(15)+coi(16)+coi(17)*2+coi(18)+coi(19)+coi(20)+coi(21)+coi(22)+coi(23)+coi(24)*3+coi(25)
prodf(41)=f_kf(36,Tint,Tdiff)*(coi(18))*M
prodb(41)=f_kb(36,Tint,Tdiff)*(coi(11))*(coi(4))*M
! reac n. 42: CH3 + OH <=> CH2O + H2
prodf(42)=f_kf(37,Tint,Tdiff)*(coi(16))*(coi(7))
prodb(42)=f_kb(37,Tint,Tdiff)*(coi(19))*(coi(3))
! reac n. 43: CH2 + CH3 <=> C2H4 + H
prodf(43)=f_kf(38,Tint,Tdiff)*(coi(15))*(coi(16))
prodb(43)=f_kb(38,Tint,Tdiff)*(coi(22))*(coi(4))
! reac n. 44: CO + O2 <=> CO2 + O
prodf(44)=f_kf(39,Tint,Tdiff)*(coi(11))*(coi(6))
prodb(44)=f_kb(39,Tint,Tdiff)*(coi(12))*(coi(5))
! reac n. 45: CO + OH <=> CO2 + H
prodf(45)=f_kf(40,Tint,Tdiff)*(coi(11))*(coi(7))
prodb(45)=f_kb(40,Tint,Tdiff)*(coi(12))*(coi(4))
! reac n. 46: CH2O + OH <=> H2O + HCO
prodf(46)=f_kf(41,Tint,Tdiff)*(coi(19))*(coi(7))
prodb(46)=f_kb(41,Tint,Tdiff)*(coi(8))*(coi(18))
! reac n. 47: CH2O + H <=> H2 + HCO
prodf(47)=f_kf(42,Tint,Tdiff)*(coi(19))*(coi(4))
prodb(47)=f_kb(42,Tint,Tdiff)*(coi(3))*(coi(18))
! reac n. 48: CH4 + H <=> CH3 + H2
prodf(48)=f_kf(43,Tint,Tdiff)*(coi(17))*(coi(4))
prodb(48)=f_kb(43,Tint,Tdiff)*(coi(16))*(coi(3))
! reac n. 49: 2 CH3 (+M) <=> C2H6 (+M)
M=coi(1)*0.7+coi(2)+coi(3)*2+coi(4)+coi(5)+coi(6)+coi(7)+coi(8)*6+coi(9)+coi(10)+coi(11)*1.5+coi(12)*2+coi(13)+coi(14)+coi(15)+coi(16)+coi(17)*2+coi(18)+coi(19)+coi(20)+coi(21)+coi(22)+coi(23)+coi(24)*3+coi(25)
k = f_k_troe(6,Tint,Tdiff,M)
prodf(49)=k(1)*(coi(16)*coi(16))
prodb(49)=k(2)*(coi(24))
! reac n. 50: H + O2 + N2 <=> HO2 + N2
M=coi(2)
prodf(50)=f_kf(44,Tint,Tdiff)*(coi(4))*(coi(6))*M
prodb(50)=f_kb(44,Tint,Tdiff)*(coi(9))*M
! reac n. 51: H + O2 + AR <=> HO2 + AR
M=coi(1)
prodf(51)=f_kf(45,Tint,Tdiff)*(coi(4))*(coi(6))*M
prodb(51)=f_kb(45,Tint,Tdiff)*(coi(9))*M
! species source terms
omegadot(1)=0.d0
omegadot(2)=0.d0
omegadot(3)=Wm_tab(3)*(+(0.0-1.0)*(prodf(28)-prodb(28))+(1.0-0.0)*(prodf(35)-prodb(35))+(1.0-0.0)*(prodf(36)-prodb(36))+(0.0-1.0)*(prodf(37)-prodb(37))+(1.0-0.0)*(prodf(42)-prodb(42))+(1.0-0.0)*(prodf(47)-prodb(47))+(1.0-0.0)*(prodf(48)-prodb(48)))
omegadot(4)=Wm_tab(4)*(+(0.0-1.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(4)-prodb(4))+(0.0-1.0)*(prodf(26)-prodb(26))+(0.0-1.0)*(prodf(27)-prodb(27))+(1.0-0.0)*(prodf(28)-prodb(28))+(1.0-0.0)*(prodf(29)-prodb(29))+(0.0-1.0)*(prodf(31)-prodb(31))+(0.0-1.0)*(prodf(32)-prodb(32))+(0.0-1.0)*(prodf(33)-prodb(33))+(0.0-1.0)*(prodf(34)-prodb(34))+(0.0-1.0)*(prodf(35)-prodb(35))+(0.0-1.0)*(prodf(36)-prodb(36))+(1.0-0.0)*(prodf(37)-prodb(37))+(1.0-0.0)*(prodf(38)-prodb(38))+(1.0-0.0)*(prodf(41)-prodb(41))+(1.0-0.0)*(prodf(43)-prodb(43))+(1.0-0.0)*(prodf(45)-prodb(45))+(0.0-1.0)*(prodf(47)-prodb(47))+(0.0-1.0)*(prodf(48)-prodb(48))+(0.0-1.0)*(prodf(50)-prodb(50))+(0.0-1.0)*(prodf(51)-prodb(51)))
omegadot(5)=Wm_tab(5)*(+(1.0-0.0)*(prodf(16)-prodb(16))+(0.0-1.0)*(prodf(25)-prodb(25))+(1.0-0.0)*(prodf(26)-prodb(26))+(0.0-1.0)*(prodf(28)-prodb(28))+(0.0-1.0)*(prodf(29)-prodb(29))+(0.0-1.0)*(prodf(30)-prodb(30))+(1.0-0.0)*(prodf(44)-prodb(44)))
omegadot(6)=Wm_tab(6)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(0.0-1.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(10)-prodb(10))+(1.0-0.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(16)-prodb(16))+(0.0-1.0)*(prodf(17)-prodb(17))+(0.0-1.0)*(prodf(22)-prodb(22))+(0.0-1.0)*(prodf(23)-prodb(23))+(0.0-1.0)*(prodf(26)-prodb(26))+(0.0-1.0)*(prodf(27)-prodb(27))+(0.0-1.0)*(prodf(40)-prodb(40))+(0.0-1.0)*(prodf(44)-prodb(44))+(0.0-1.0)*(prodf(50)-prodb(50))+(0.0-1.0)*(prodf(51)-prodb(51)))
omegadot(7)=Wm_tab(7)*(+(0.0-2.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(6)-prodb(6))+(0.0-1.0)*(prodf(7)-prodb(7))+(0.0-1.0)*(prodf(8)-prodb(8))+(0.0-1.0)*(prodf(9)-prodb(9))+(1.0-0.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(14)-prodb(14))+(1.0-0.0)*(prodf(17)-prodb(17))+(1.0-0.0)*(prodf(25)-prodb(25))+(1.0-0.0)*(prodf(26)-prodb(26))+(1.0-0.0)*(prodf(28)-prodb(28))+(0.0-1.0)*(prodf(31)-prodb(31))+(0.0-1.0)*(prodf(37)-prodb(37))+(0.0-1.0)*(prodf(38)-prodb(38))+(0.0-1.0)*(prodf(39)-prodb(39))+(0.0-1.0)*(prodf(42)-prodb(42))+(0.0-1.0)*(prodf(45)-prodb(45))+(0.0-1.0)*(prodf(46)-prodb(46)))
omegadot(8)=Wm_tab(8)*(+(1.0-0.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(7)-prodb(7))+(1.0-0.0)*(prodf(8)-prodb(8))+(1.0-0.0)*(prodf(9)-prodb(9))+(1.0-0.0)*(prodf(31)-prodb(31))+(1.0-0.0)*(prodf(37)-prodb(37))+(1.0-0.0)*(prodf(39)-prodb(39))+(1.0-0.0)*(prodf(46)-prodb(46)))
omegadot(9)=Wm_tab(9)*(+(1.0-0.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(2)-prodb(2))+(1.0-0.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(7)-prodb(7))+(1.0-0.0)*(prodf(8)-prodb(8))+(0.0-2.0)*(prodf(10)-prodb(10))+(0.0-2.0)*(prodf(11)-prodb(11))+(0.0-1.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(13)-prodb(13))+(0.0-1.0)*(prodf(14)-prodb(14))+(0.0-1.0)*(prodf(15)-prodb(15))+(1.0-0.0)*(prodf(18)-prodb(18))+(0.0-1.0)*(prodf(20)-prodb(20))+(1.0-0.0)*(prodf(22)-prodb(22))+(1.0-0.0)*(prodf(27)-prodb(27))+(1.0-0.0)*(prodf(40)-prodb(40))+(1.0-0.0)*(prodf(50)-prodb(50))+(1.0-0.0)*(prodf(51)-prodb(51)))
omegadot(10)=Wm_tab(10)*(+(1.0-0.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(7)-prodb(7))+(0.0-1.0)*(prodf(8)-prodb(8))+(1.0-0.0)*(prodf(10)-prodb(10))+(1.0-0.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(15)-prodb(15))+(0.0-1.0)*(prodf(18)-prodb(18))+(1.0-0.0)*(prodf(20)-prodb(20)))
omegadot(11)=Wm_tab(11)*(+(0.0-1.0)*(prodf(14)-prodb(14))+(0.0-1.0)*(prodf(30)-prodb(30))+(1.0-0.0)*(prodf(40)-prodb(40))+(1.0-0.0)*(prodf(41)-prodb(41))+(0.0-1.0)*(prodf(44)-prodb(44))+(0.0-1.0)*(prodf(45)-prodb(45)))
omegadot(12)=Wm_tab(12)*(+(1.0-0.0)*(prodf(14)-prodb(14))+(1.0-0.0)*(prodf(30)-prodb(30))+(1.0-0.0)*(prodf(44)-prodb(44))+(1.0-0.0)*(prodf(45)-prodb(45)))
omegadot(13)=0.d0
omegadot(14)=0.d0
omegadot(15)=Wm_tab(15)*(+(0.0-1.0)*(prodf(38)-prodb(38))+(0.0-1.0)*(prodf(43)-prodb(43)))
omegadot(16)=Wm_tab(16)*(+(1.0-0.0)*(prodf(9)-prodb(9))+(0.0-1.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(13)-prodb(13))+(0.0-1.0)*(prodf(16)-prodb(16))+(0.0-1.0)*(prodf(17)-prodb(17))+(0.0-1.0)*(prodf(18)-prodb(18))+(0.0-1.0)*(prodf(19)-prodb(19))+(0.0-1.0)*(prodf(21)-prodb(21))+(0.0-1.0)*(prodf(23)-prodb(23))+(0.0-1.0)*(prodf(24)-prodb(24))+(1.0-0.0)*(prodf(25)-prodb(25))+(0.0-1.0)*(prodf(29)-prodb(29))+(0.0-1.0)*(prodf(32)-prodb(32))+(0.0-1.0)*(prodf(42)-prodb(42))+(0.0-1.0)*(prodf(43)-prodb(43))+(1.0-0.0)*(prodf(48)-prodb(48))+(0.0-2.0)*(prodf(49)-prodb(49)))
omegadot(17)=Wm_tab(17)*(+(0.0-1.0)*(prodf(9)-prodb(9))+(1.0-0.0)*(prodf(12)-prodb(12))+(1.0-0.0)*(prodf(18)-prodb(18))+(1.0-0.0)*(prodf(19)-prodb(19))+(1.0-0.0)*(prodf(24)-prodb(24))+(0.0-1.0)*(prodf(25)-prodb(25))+(1.0-0.0)*(prodf(32)-prodb(32))+(0.0-1.0)*(prodf(48)-prodb(48)))
omegadot(18)=Wm_tab(18)*(+(1.0-0.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(15)-prodb(15))+(1.0-0.0)*(prodf(19)-prodb(19))+(0.0-1.0)*(prodf(33)-prodb(33))+(0.0-1.0)*(prodf(40)-prodb(40))+(0.0-1.0)*(prodf(41)-prodb(41))+(1.0-0.0)*(prodf(46)-prodb(46))+(1.0-0.0)*(prodf(47)-prodb(47)))
omegadot(19)=Wm_tab(19)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(0.0-1.0)*(prodf(4)-prodb(4))+(0.0-1.0)*(prodf(15)-prodb(15))+(1.0-0.0)*(prodf(17)-prodb(17))+(0.0-1.0)*(prodf(19)-prodb(19))+(1.0-0.0)*(prodf(20)-prodb(20))+(1.0-0.0)*(prodf(22)-prodb(22))+(1.0-0.0)*(prodf(24)-prodb(24))+(1.0-0.0)*(prodf(29)-prodb(29))+(1.0-0.0)*(prodf(33)-prodb(33))+(1.0-0.0)*(prodf(38)-prodb(38))+(1.0-0.0)*(prodf(42)-prodb(42))+(0.0-1.0)*(prodf(46)-prodb(46))+(0.0-1.0)*(prodf(47)-prodb(47)))
omegadot(20)=Wm_tab(20)*(+(1.0-0.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(16)-prodb(16))+(0.0-1.0)*(prodf(20)-prodb(20))+(2.0-0.0)*(prodf(21)-prodb(21))+(0.0-1.0)*(prodf(22)-prodb(22))+(0.0-1.0)*(prodf(24)-prodb(24)))
omegadot(21)=Wm_tab(21)*(+(1.0-0.0)*(prodf(35)-prodb(35)))
omegadot(22)=Wm_tab(22)*(+(0.0-1.0)*(prodf(34)-prodb(34))+(0.0-1.0)*(prodf(35)-prodb(35))+(1.0-0.0)*(prodf(43)-prodb(43)))
omegadot(23)=Wm_tab(23)*(+(1.0-0.0)*(prodf(34)-prodb(34))+(1.0-0.0)*(prodf(36)-prodb(36))+(1.0-0.0)*(prodf(39)-prodb(39)))
omegadot(24)=Wm_tab(24)*(+(0.0-1.0)*(prodf(36)-prodb(36))+(0.0-1.0)*(prodf(39)-prodb(39))+(1.0-0.0)*(prodf(49)-prodb(49)))
omegadot(25)=Wm_tab(25)*(+(0.0-1.0)*(prodf(21)-prodb(21))+(1.0-0.0)*(prodf(23)-prodb(23)))
end subroutine ZK
end module ZK_mod
