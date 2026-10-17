module ONERA7_mod
implicit none
contains
subroutine ONERA_7(roi,temp,omegadot)
!$acc routine seq
use FLINT_Lib_Thermodynamic
use FLINT_Lib_Chemistry_data
use FLINT_Lib_Chemistry_falloff
implicit none
real(8), intent(inout)  :: roi(ns)
real(8), intent(in)  :: temp
real(8), intent(out) :: omegadot(ns)

real(8) :: coi(NSCHEM), Tdiff   ! fixed-size local (device: no malloc); 1:ns used
real(8) :: M !< Third body
integer :: is, T_i, Tint(2)
real(8) :: prodf(1:14), prodb(1:14)
real(8) :: k(2) !< Falloff rate coefficients


do is = 1, ns 
 coi(is)=roi(is)/Wm_tab(is) ! kmol/m^3
enddo 
T_i = int(temp) 
Tdiff  = temp-T_i 
Tint(1) = T_i 
Tint(2) = T_i + 1 
! reac n. 1: H2 + O2 => 2 OH
prodf(1)=f_kf(1,Tint,Tdiff)*(coi(3))*(coi(1))
prodb(1)=f_kb(1,Tint,Tdiff)*(coi(6)*coi(6))
! reac n. 2: 2 OH => H2 + O2
prodf(2)=f_kf(2,Tint,Tdiff)*(coi(6)*coi(6))
prodb(2)=f_kb(2,Tint,Tdiff)*(coi(3))*(coi(1))
! reac n. 3: H + O2 => O + OH
prodf(3)=f_kf(3,Tint,Tdiff)*(coi(4))*(coi(1))
prodb(3)=f_kb(3,Tint,Tdiff)*(coi(5))*(coi(6))
! reac n. 4: O + OH => H + O2
prodf(4)=f_kf(4,Tint,Tdiff)*(coi(5))*(coi(6))
prodb(4)=f_kb(4,Tint,Tdiff)*(coi(4))*(coi(1))
! reac n. 5: H2 + OH => H + H2O
prodf(5)=f_kf(5,Tint,Tdiff)*(coi(3))*(coi(6))
prodb(5)=f_kb(5,Tint,Tdiff)*(coi(4))*(coi(2))
! reac n. 6: H + H2O => H2 + OH
prodf(6)=f_kf(6,Tint,Tdiff)*(coi(4))*(coi(2))
prodb(6)=f_kb(6,Tint,Tdiff)*(coi(3))*(coi(6))
! reac n. 7: H2 + O => H + OH
prodf(7)=f_kf(7,Tint,Tdiff)*(coi(3))*(coi(5))
prodb(7)=f_kb(7,Tint,Tdiff)*(coi(4))*(coi(6))
! reac n. 8: H + OH => H2 + O
prodf(8)=f_kf(8,Tint,Tdiff)*(coi(4))*(coi(6))
prodb(8)=f_kb(8,Tint,Tdiff)*(coi(3))*(coi(5))
! reac n. 9: 2 OH => H2O + O
prodf(9)=f_kf(9,Tint,Tdiff)*(coi(6)*coi(6))
prodb(9)=f_kb(9,Tint,Tdiff)*(coi(2))*(coi(5))
! reac n. 10: H2O + O => 2 OH
prodf(10)=f_kf(10,Tint,Tdiff)*(coi(2))*(coi(5))
prodb(10)=f_kb(10,Tint,Tdiff)*(coi(6)*coi(6))
! reac n. 11: H + OH + M => H2O + M
M=coi(1)+coi(2)*12+coi(3)*2.5+coi(4)+coi(5)+coi(6)+coi(7)
prodf(11)=f_kf(11,Tint,Tdiff)*(coi(4))*(coi(6))*M
prodb(11)=f_kb(11,Tint,Tdiff)*(coi(2))*M
! reac n. 12: H2O + M => H + OH + M
M=coi(1)+coi(2)*12+coi(3)*2.5+coi(4)+coi(5)+coi(6)+coi(7)
prodf(12)=f_kf(12,Tint,Tdiff)*(coi(2))*M
prodb(12)=f_kb(12,Tint,Tdiff)*(coi(4))*(coi(6))*M
! reac n. 13: 2 H + M => H2 + M
M=coi(1)+coi(2)*12+coi(3)*2.5+coi(4)+coi(5)+coi(6)+coi(7)
prodf(13)=f_kf(13,Tint,Tdiff)*(coi(4)*coi(4))*M
prodb(13)=f_kb(13,Tint,Tdiff)*(coi(3))*M
! reac n. 14: H2 + M => 2 H + M
M=coi(1)+coi(2)*12+coi(3)*2.5+coi(4)+coi(5)+coi(6)+coi(7)
prodf(14)=f_kf(14,Tint,Tdiff)*(coi(3))*M
prodb(14)=f_kb(14,Tint,Tdiff)*(coi(4)*coi(4))*M
! species source terms
omegadot(1)=Wm_tab(1)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(4)-prodb(4)))
omegadot(2)=Wm_tab(2)*(+(1.0-0.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(9)-prodb(9))+(0.0-1.0)*(prodf(10)-prodb(10))+(1.0-0.0)*(prodf(11)-prodb(11))+(0.0-1.0)*(prodf(12)-prodb(12)))
omegadot(3)=Wm_tab(3)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(5)-prodb(5))+(1.0-0.0)*(prodf(6)-prodb(6))+(0.0-1.0)*(prodf(7)-prodb(7))+(1.0-0.0)*(prodf(8)-prodb(8))+(1.0-0.0)*(prodf(13)-prodb(13))+(0.0-1.0)*(prodf(14)-prodb(14)))
omegadot(4)=Wm_tab(4)*(+(0.0-1.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(7)-prodb(7))+(0.0-1.0)*(prodf(8)-prodb(8))+(0.0-1.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(12)-prodb(12))+(0.0-2.0)*(prodf(13)-prodb(13))+(2.0-0.0)*(prodf(14)-prodb(14)))
omegadot(5)=Wm_tab(5)*(+(1.0-0.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(4)-prodb(4))+(0.0-1.0)*(prodf(7)-prodb(7))+(1.0-0.0)*(prodf(8)-prodb(8))+(1.0-0.0)*(prodf(9)-prodb(9))+(0.0-1.0)*(prodf(10)-prodb(10)))
omegadot(6)=Wm_tab(6)*(+(2.0-0.0)*(prodf(1)-prodb(1))+(0.0-2.0)*(prodf(2)-prodb(2))+(1.0-0.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(4)-prodb(4))+(0.0-1.0)*(prodf(5)-prodb(5))+(1.0-0.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(7)-prodb(7))+(0.0-1.0)*(prodf(8)-prodb(8))+(0.0-2.0)*(prodf(9)-prodb(9))+(2.0-0.0)*(prodf(10)-prodb(10))+(0.0-1.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(12)-prodb(12)))
omegadot(7)=0.d0
end subroutine ONERA_7


!------------------------------------------------------------------------------
! Analytical Jacobian of ONERA_7's mass-production rate w.r.t. species mass
! densities and temperature.
!
!   dwdr(i,j) = d omegadot(i) / d roi(j)   [1/s]
!   dwdT(i)   = d omegadot(i) / d T        [kg/(m^3 s K)]
!
! Convention matches ONERA_7 above:
!   - coi(s) = roi(s) / Wm_tab(s)                        [kmol/m^3]
!   - rates linearly interpolated in T from kf_tab/kb_tab tables, so
!     dkf/dT = kf_tab(T_i+1,r) - kf_tab(T_i,r) is a free table lookup
!   - third-body efficiencies for reactions 11..14:
!     epsM = [1, 12, 2.5, 1, 1, 1, 1] (matches the hand-written M expression
!     used in the RHS above; do NOT change without updating ONERA_7).
!------------------------------------------------------------------------------
subroutine ONERA_7_jac(roi, temp, dwdr, dwdT)
!$acc routine seq
use FLINT_Lib_Thermodynamic
use FLINT_Lib_Chemistry_data
implicit none
real(8), intent(in)  :: roi(ns), temp
! dwdr: fixed leading dimension NSCHEM to match jac_native/chemjac_if storage
! (sequence association); only (1:ns,1:ns) written.
real(8), intent(out) :: dwdr(NSCHEM,ns)
real(8), intent(out) :: dwdT(ns)

! Local (fixed-size: device builds avoid per-call heap malloc; 1:ns used)
real(8) :: coi(NSCHEM), Tdiff
integer :: T_i, j
real(8) :: kf_r(14), kb_r(14)
real(8) :: dkf_dT(14), dkb_dT(14)
real(8) :: M, dM_dc(NSCHEM)
real(8) :: dRf_dc(NSCHEM), dRb_dc(NSCHEM), dnet_dc(NSCHEM)
real(8) :: dRf_dT_r, dRb_dT_r, dnet_dT_r
real(8) :: dwdr_c(NSCHEM,NSCHEM)     ! wrt coi; converted to wrt roi at the end
real(8), parameter :: epsM(7) = [1.d0, 12.d0, 2.5d0, 1.d0, 1.d0, 1.d0, 1.d0]

! 1) coi and T setup ---------------------------------------------------------
do j = 1, ns
  coi(j) = roi(j)/Wm_tab(j)
enddo
T_i   = int(temp)
Tdiff = temp - T_i

! 2) Forward/backward rate constants and their T-derivatives -----------------
!    f_k(r) = kf_tab(T_i,r) + Tdiff*(kf_tab(T_i+1,r) - kf_tab(T_i,r))
!    so dk/dT = kf_tab(T_i+1,r) - kf_tab(T_i,r)  (linear interp slope)
do j = 1, 14
  kf_r(j)   =  kf_tab(T_i  ,j) + Tdiff*(kf_tab(T_i+1,j) - kf_tab(T_i,j))
  kb_r(j)   =  kb_tab(T_i  ,j) + Tdiff*(kb_tab(T_i+1,j) - kb_tab(T_i,j))
  dkf_dT(j) =  kf_tab(T_i+1,j) - kf_tab(T_i,j)
  dkb_dT(j) =  kb_tab(T_i+1,j) - kb_tab(T_i,j)
enddo

! 3) Third-body M and its derivatives ----------------------------------------
M = coi(1) + 12.d0*coi(2) + 2.5d0*coi(3) + coi(4) + coi(5) + coi(6) + coi(7)
dM_dc(:) = epsM(:)

! 4) Per-reaction accumulation -----------------------------------------------
dwdr_c = 0.d0
dwdT   = 0.d0

! ---- reac 1: H2(3) + O2(1) => 2 OH(6)        prodf = kf*c3*c1, prodb = kb*c6*c6
dRf_dc = 0.d0; dRb_dc = 0.d0
dRf_dc(1) = kf_r(1)*coi(3);   dRf_dc(3) = kf_r(1)*coi(1)
dRb_dc(6) = 2.d0*kb_r(1)*coi(6)
dRf_dT_r = dkf_dT(1)*coi(3)*coi(1)
dRb_dT_r = dkb_dT(1)*coi(6)*coi(6)
dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
! stoich: O2(1):-1, H2(3):-1, OH(6):+2
dwdr_c(1,:) = dwdr_c(1,:) -      Wm_tab(1)*dnet_dc
dwdr_c(3,:) = dwdr_c(3,:) -      Wm_tab(3)*dnet_dc
dwdr_c(6,:) = dwdr_c(6,:) + 2.d0*Wm_tab(6)*dnet_dc
dwdT(1) = dwdT(1) -      Wm_tab(1)*dnet_dT_r
dwdT(3) = dwdT(3) -      Wm_tab(3)*dnet_dT_r
dwdT(6) = dwdT(6) + 2.d0*Wm_tab(6)*dnet_dT_r

! ---- reac 2: 2 OH => H2 + O2                 prodf = kf*c6*c6, prodb = kb*c3*c1
dRf_dc = 0.d0; dRb_dc = 0.d0
dRf_dc(6) = 2.d0*kf_r(2)*coi(6)
dRb_dc(1) = kb_r(2)*coi(3);   dRb_dc(3) = kb_r(2)*coi(1)
dRf_dT_r = dkf_dT(2)*coi(6)*coi(6)
dRb_dT_r = dkb_dT(2)*coi(3)*coi(1)
dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
! stoich: O2(1):+1, H2(3):+1, OH(6):-2
dwdr_c(1,:) = dwdr_c(1,:) +      Wm_tab(1)*dnet_dc
dwdr_c(3,:) = dwdr_c(3,:) +      Wm_tab(3)*dnet_dc
dwdr_c(6,:) = dwdr_c(6,:) - 2.d0*Wm_tab(6)*dnet_dc
dwdT(1) = dwdT(1) +      Wm_tab(1)*dnet_dT_r
dwdT(3) = dwdT(3) +      Wm_tab(3)*dnet_dT_r
dwdT(6) = dwdT(6) - 2.d0*Wm_tab(6)*dnet_dT_r

! ---- reac 3: H(4) + O2(1) => O(5) + OH(6)    prodf = kf*c4*c1, prodb = kb*c5*c6
dRf_dc = 0.d0; dRb_dc = 0.d0
dRf_dc(1) = kf_r(3)*coi(4);   dRf_dc(4) = kf_r(3)*coi(1)
dRb_dc(5) = kb_r(3)*coi(6);   dRb_dc(6) = kb_r(3)*coi(5)
dRf_dT_r = dkf_dT(3)*coi(4)*coi(1)
dRb_dT_r = dkb_dT(3)*coi(5)*coi(6)
dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
! stoich: O2(1):-1, H(4):-1, O(5):+1, OH(6):+1
dwdr_c(1,:) = dwdr_c(1,:) - Wm_tab(1)*dnet_dc
dwdr_c(4,:) = dwdr_c(4,:) - Wm_tab(4)*dnet_dc
dwdr_c(5,:) = dwdr_c(5,:) + Wm_tab(5)*dnet_dc
dwdr_c(6,:) = dwdr_c(6,:) + Wm_tab(6)*dnet_dc
dwdT(1) = dwdT(1) - Wm_tab(1)*dnet_dT_r
dwdT(4) = dwdT(4) - Wm_tab(4)*dnet_dT_r
dwdT(5) = dwdT(5) + Wm_tab(5)*dnet_dT_r
dwdT(6) = dwdT(6) + Wm_tab(6)*dnet_dT_r

! ---- reac 4: O(5) + OH(6) => H(4) + O2(1)    prodf = kf*c5*c6, prodb = kb*c4*c1
dRf_dc = 0.d0; dRb_dc = 0.d0
dRf_dc(5) = kf_r(4)*coi(6);   dRf_dc(6) = kf_r(4)*coi(5)
dRb_dc(1) = kb_r(4)*coi(4);   dRb_dc(4) = kb_r(4)*coi(1)
dRf_dT_r = dkf_dT(4)*coi(5)*coi(6)
dRb_dT_r = dkb_dT(4)*coi(4)*coi(1)
dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
! stoich: O2(1):+1, H(4):+1, O(5):-1, OH(6):-1
dwdr_c(1,:) = dwdr_c(1,:) + Wm_tab(1)*dnet_dc
dwdr_c(4,:) = dwdr_c(4,:) + Wm_tab(4)*dnet_dc
dwdr_c(5,:) = dwdr_c(5,:) - Wm_tab(5)*dnet_dc
dwdr_c(6,:) = dwdr_c(6,:) - Wm_tab(6)*dnet_dc
dwdT(1) = dwdT(1) + Wm_tab(1)*dnet_dT_r
dwdT(4) = dwdT(4) + Wm_tab(4)*dnet_dT_r
dwdT(5) = dwdT(5) - Wm_tab(5)*dnet_dT_r
dwdT(6) = dwdT(6) - Wm_tab(6)*dnet_dT_r

! ---- reac 5: H2(3) + OH(6) => H(4) + H2O(2)  prodf = kf*c3*c6, prodb = kb*c4*c2
dRf_dc = 0.d0; dRb_dc = 0.d0
dRf_dc(3) = kf_r(5)*coi(6);   dRf_dc(6) = kf_r(5)*coi(3)
dRb_dc(2) = kb_r(5)*coi(4);   dRb_dc(4) = kb_r(5)*coi(2)
dRf_dT_r = dkf_dT(5)*coi(3)*coi(6)
dRb_dT_r = dkb_dT(5)*coi(4)*coi(2)
dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
! stoich: H2O(2):+1, H2(3):-1, H(4):+1, OH(6):-1
dwdr_c(2,:) = dwdr_c(2,:) + Wm_tab(2)*dnet_dc
dwdr_c(3,:) = dwdr_c(3,:) - Wm_tab(3)*dnet_dc
dwdr_c(4,:) = dwdr_c(4,:) + Wm_tab(4)*dnet_dc
dwdr_c(6,:) = dwdr_c(6,:) - Wm_tab(6)*dnet_dc
dwdT(2) = dwdT(2) + Wm_tab(2)*dnet_dT_r
dwdT(3) = dwdT(3) - Wm_tab(3)*dnet_dT_r
dwdT(4) = dwdT(4) + Wm_tab(4)*dnet_dT_r
dwdT(6) = dwdT(6) - Wm_tab(6)*dnet_dT_r

! ---- reac 6: H(4) + H2O(2) => H2(3) + OH(6)  prodf = kf*c4*c2, prodb = kb*c3*c6
dRf_dc = 0.d0; dRb_dc = 0.d0
dRf_dc(2) = kf_r(6)*coi(4);   dRf_dc(4) = kf_r(6)*coi(2)
dRb_dc(3) = kb_r(6)*coi(6);   dRb_dc(6) = kb_r(6)*coi(3)
dRf_dT_r = dkf_dT(6)*coi(4)*coi(2)
dRb_dT_r = dkb_dT(6)*coi(3)*coi(6)
dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
! stoich: H2O(2):-1, H2(3):+1, H(4):-1, OH(6):+1
dwdr_c(2,:) = dwdr_c(2,:) - Wm_tab(2)*dnet_dc
dwdr_c(3,:) = dwdr_c(3,:) + Wm_tab(3)*dnet_dc
dwdr_c(4,:) = dwdr_c(4,:) - Wm_tab(4)*dnet_dc
dwdr_c(6,:) = dwdr_c(6,:) + Wm_tab(6)*dnet_dc
dwdT(2) = dwdT(2) - Wm_tab(2)*dnet_dT_r
dwdT(3) = dwdT(3) + Wm_tab(3)*dnet_dT_r
dwdT(4) = dwdT(4) - Wm_tab(4)*dnet_dT_r
dwdT(6) = dwdT(6) + Wm_tab(6)*dnet_dT_r

! ---- reac 7: H2(3) + O(5) => H(4) + OH(6)    prodf = kf*c3*c5, prodb = kb*c4*c6
dRf_dc = 0.d0; dRb_dc = 0.d0
dRf_dc(3) = kf_r(7)*coi(5);   dRf_dc(5) = kf_r(7)*coi(3)
dRb_dc(4) = kb_r(7)*coi(6);   dRb_dc(6) = kb_r(7)*coi(4)
dRf_dT_r = dkf_dT(7)*coi(3)*coi(5)
dRb_dT_r = dkb_dT(7)*coi(4)*coi(6)
dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
! stoich: H2(3):-1, H(4):+1, O(5):-1, OH(6):+1
dwdr_c(3,:) = dwdr_c(3,:) - Wm_tab(3)*dnet_dc
dwdr_c(4,:) = dwdr_c(4,:) + Wm_tab(4)*dnet_dc
dwdr_c(5,:) = dwdr_c(5,:) - Wm_tab(5)*dnet_dc
dwdr_c(6,:) = dwdr_c(6,:) + Wm_tab(6)*dnet_dc
dwdT(3) = dwdT(3) - Wm_tab(3)*dnet_dT_r
dwdT(4) = dwdT(4) + Wm_tab(4)*dnet_dT_r
dwdT(5) = dwdT(5) - Wm_tab(5)*dnet_dT_r
dwdT(6) = dwdT(6) + Wm_tab(6)*dnet_dT_r

! ---- reac 8: H(4) + OH(6) => H2(3) + O(5)    prodf = kf*c4*c6, prodb = kb*c3*c5
dRf_dc = 0.d0; dRb_dc = 0.d0
dRf_dc(4) = kf_r(8)*coi(6);   dRf_dc(6) = kf_r(8)*coi(4)
dRb_dc(3) = kb_r(8)*coi(5);   dRb_dc(5) = kb_r(8)*coi(3)
dRf_dT_r = dkf_dT(8)*coi(4)*coi(6)
dRb_dT_r = dkb_dT(8)*coi(3)*coi(5)
dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
! stoich: H2(3):+1, H(4):-1, O(5):+1, OH(6):-1
dwdr_c(3,:) = dwdr_c(3,:) + Wm_tab(3)*dnet_dc
dwdr_c(4,:) = dwdr_c(4,:) - Wm_tab(4)*dnet_dc
dwdr_c(5,:) = dwdr_c(5,:) + Wm_tab(5)*dnet_dc
dwdr_c(6,:) = dwdr_c(6,:) - Wm_tab(6)*dnet_dc
dwdT(3) = dwdT(3) + Wm_tab(3)*dnet_dT_r
dwdT(4) = dwdT(4) - Wm_tab(4)*dnet_dT_r
dwdT(5) = dwdT(5) + Wm_tab(5)*dnet_dT_r
dwdT(6) = dwdT(6) - Wm_tab(6)*dnet_dT_r

! ---- reac 9: 2 OH(6) => H2O(2) + O(5)        prodf = kf*c6*c6, prodb = kb*c2*c5
dRf_dc = 0.d0; dRb_dc = 0.d0
dRf_dc(6) = 2.d0*kf_r(9)*coi(6)
dRb_dc(2) = kb_r(9)*coi(5);   dRb_dc(5) = kb_r(9)*coi(2)
dRf_dT_r = dkf_dT(9)*coi(6)*coi(6)
dRb_dT_r = dkb_dT(9)*coi(2)*coi(5)
dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
! stoich: H2O(2):+1, O(5):+1, OH(6):-2
dwdr_c(2,:) = dwdr_c(2,:) +      Wm_tab(2)*dnet_dc
dwdr_c(5,:) = dwdr_c(5,:) +      Wm_tab(5)*dnet_dc
dwdr_c(6,:) = dwdr_c(6,:) - 2.d0*Wm_tab(6)*dnet_dc
dwdT(2) = dwdT(2) +      Wm_tab(2)*dnet_dT_r
dwdT(5) = dwdT(5) +      Wm_tab(5)*dnet_dT_r
dwdT(6) = dwdT(6) - 2.d0*Wm_tab(6)*dnet_dT_r

! ---- reac 10: H2O(2) + O(5) => 2 OH(6)       prodf = kf*c2*c5, prodb = kb*c6*c6
dRf_dc = 0.d0; dRb_dc = 0.d0
dRf_dc(2) = kf_r(10)*coi(5);  dRf_dc(5) = kf_r(10)*coi(2)
dRb_dc(6) = 2.d0*kb_r(10)*coi(6)
dRf_dT_r = dkf_dT(10)*coi(2)*coi(5)
dRb_dT_r = dkb_dT(10)*coi(6)*coi(6)
dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
! stoich: H2O(2):-1, O(5):-1, OH(6):+2
dwdr_c(2,:) = dwdr_c(2,:) -      Wm_tab(2)*dnet_dc
dwdr_c(5,:) = dwdr_c(5,:) -      Wm_tab(5)*dnet_dc
dwdr_c(6,:) = dwdr_c(6,:) + 2.d0*Wm_tab(6)*dnet_dc
dwdT(2) = dwdT(2) -      Wm_tab(2)*dnet_dT_r
dwdT(5) = dwdT(5) -      Wm_tab(5)*dnet_dT_r
dwdT(6) = dwdT(6) + 2.d0*Wm_tab(6)*dnet_dT_r

! ---- reac 11: H(4) + OH(6) + M => H2O(2) + M
! prodf = kf*c4*c6*M
! d/dc(j) = kf*[(c6 if j==4) + (c4 if j==6)]*M + kf*c4*c6*dM/dc(j)
! prodb = kb*c2*M
! d/dc(j) = kb*[(1 if j==2)]*M + kb*c2*dM/dc(j)
dRf_dc(:) = kf_r(11)*coi(4)*coi(6)*dM_dc(:)
dRf_dc(4) = dRf_dc(4) + kf_r(11)*coi(6)*M
dRf_dc(6) = dRf_dc(6) + kf_r(11)*coi(4)*M
dRb_dc(:) = kb_r(11)*coi(2)*dM_dc(:)
dRb_dc(2) = dRb_dc(2) + kb_r(11)*M
dRf_dT_r = dkf_dT(11)*coi(4)*coi(6)*M
dRb_dT_r = dkb_dT(11)*coi(2)*M
dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
! stoich: H2O(2):+1, H(4):-1, OH(6):-1
dwdr_c(2,:) = dwdr_c(2,:) + Wm_tab(2)*dnet_dc
dwdr_c(4,:) = dwdr_c(4,:) - Wm_tab(4)*dnet_dc
dwdr_c(6,:) = dwdr_c(6,:) - Wm_tab(6)*dnet_dc
dwdT(2) = dwdT(2) + Wm_tab(2)*dnet_dT_r
dwdT(4) = dwdT(4) - Wm_tab(4)*dnet_dT_r
dwdT(6) = dwdT(6) - Wm_tab(6)*dnet_dT_r

! ---- reac 12: H2O(2) + M => H(4) + OH(6) + M
! prodf = kf*c2*M ;  prodb = kb*c4*c6*M
dRf_dc(:) = kf_r(12)*coi(2)*dM_dc(:)
dRf_dc(2) = dRf_dc(2) + kf_r(12)*M
dRb_dc(:) = kb_r(12)*coi(4)*coi(6)*dM_dc(:)
dRb_dc(4) = dRb_dc(4) + kb_r(12)*coi(6)*M
dRb_dc(6) = dRb_dc(6) + kb_r(12)*coi(4)*M
dRf_dT_r = dkf_dT(12)*coi(2)*M
dRb_dT_r = dkb_dT(12)*coi(4)*coi(6)*M
dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
! stoich: H2O(2):-1, H(4):+1, OH(6):+1
dwdr_c(2,:) = dwdr_c(2,:) - Wm_tab(2)*dnet_dc
dwdr_c(4,:) = dwdr_c(4,:) + Wm_tab(4)*dnet_dc
dwdr_c(6,:) = dwdr_c(6,:) + Wm_tab(6)*dnet_dc
dwdT(2) = dwdT(2) - Wm_tab(2)*dnet_dT_r
dwdT(4) = dwdT(4) + Wm_tab(4)*dnet_dT_r
dwdT(6) = dwdT(6) + Wm_tab(6)*dnet_dT_r

! ---- reac 13: 2 H(4) + M => H2(3) + M
! prodf = kf*c4*c4*M ;  prodb = kb*c3*M
dRf_dc(:) = kf_r(13)*coi(4)*coi(4)*dM_dc(:)
dRf_dc(4) = dRf_dc(4) + 2.d0*kf_r(13)*coi(4)*M
dRb_dc(:) = kb_r(13)*coi(3)*dM_dc(:)
dRb_dc(3) = dRb_dc(3) + kb_r(13)*M
dRf_dT_r = dkf_dT(13)*coi(4)*coi(4)*M
dRb_dT_r = dkb_dT(13)*coi(3)*M
dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
! stoich: H2(3):+1, H(4):-2
dwdr_c(3,:) = dwdr_c(3,:) +      Wm_tab(3)*dnet_dc
dwdr_c(4,:) = dwdr_c(4,:) - 2.d0*Wm_tab(4)*dnet_dc
dwdT(3) = dwdT(3) +      Wm_tab(3)*dnet_dT_r
dwdT(4) = dwdT(4) - 2.d0*Wm_tab(4)*dnet_dT_r

! ---- reac 14: H2(3) + M => 2 H(4) + M
! prodf = kf*c3*M ;  prodb = kb*c4*c4*M
dRf_dc(:) = kf_r(14)*coi(3)*dM_dc(:)
dRf_dc(3) = dRf_dc(3) + kf_r(14)*M
dRb_dc(:) = kb_r(14)*coi(4)*coi(4)*dM_dc(:)
dRb_dc(4) = dRb_dc(4) + 2.d0*kb_r(14)*coi(4)*M
dRf_dT_r = dkf_dT(14)*coi(3)*M
dRb_dT_r = dkb_dT(14)*coi(4)*coi(4)*M
dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
! stoich: H2(3):-1, H(4):+2
dwdr_c(3,:) = dwdr_c(3,:) -      Wm_tab(3)*dnet_dc
dwdr_c(4,:) = dwdr_c(4,:) + 2.d0*Wm_tab(4)*dnet_dc
dwdT(3) = dwdT(3) -      Wm_tab(3)*dnet_dT_r
dwdT(4) = dwdT(4) + 2.d0*Wm_tab(4)*dnet_dT_r

! 5) Convert ∂/∂coi(j) → ∂/∂roi(j): divide each column by Wm_tab(j) -----------
do j = 1, ns
  dwdr(:,j) = dwdr_c(:,j) / Wm_tab(j)
enddo

! N2 (species 7) is inert: row stays zero by construction.

end subroutine ONERA_7_jac

end module ONERA7_mod
