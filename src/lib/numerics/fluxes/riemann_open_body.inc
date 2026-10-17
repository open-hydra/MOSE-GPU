! ---------------------------------------------------------------------------------------------------------------------
! The MOSE_open Riemann solvers that the device kernels do not inline, one body for the two device twins:
!   riemann_open_d   (attributes(device), Mod_Convective_W_cuf: the CUF convective kernels)
!   riemann_open (!$acc routine seq, Lib_Riemann: the record BC kernels and the OpenACC BC routines)
! rs = 4 AUSM+        MOSE_open Lib_Riemann_AUSM riemann_AUSMp
!      5 AUSM+M       MOSE_open Lib_Riemann_AUSM riemann_AUSMp_M (Chen et al. 2020)
!      6 HLLE++       MOSE_open Lib_Riemann_HLL  riemann_HLLEpp (Tramel et al.)
!      7 HLLC+Tramel  MOSE_open Lib_Riemann_HLL  riemann_HLLCp_Tramel: riemann_HLLC and riemann_HLLE blended by the sensor
!      8 LLF          MOSE_open Lib_Riemann_LF   riemann_LLF (Rusanov): central flux plus spectral-radius dissipation
!      9 SLAU2        MOSE_open Lib_Riemann_SLAU riemann_SLAU2: SLAU with the Kitamura-Shima 2013 pressure flux
!     10 HLLC+Chen    MOSE_open Lib_Riemann_HLL  riemann_HLLCp_Chen: HLLC Batten plus the Chen low-Mach and shear terms
!                     (bet = the Chen sensor h of the face, mco the floor of the local Mach number)
!     11 exact        MOSE_open Lib_Riemann_Godunov riemann_exact: Newton iteration on the star pressure with the
!                     gamma of each state (f_gamma); the face enthalpy carries the normal velocity only, as there
!     12 LMRoe        MOSE_open Lib_Riemann_Roe  Riemann_LMRoe: Roe with the Rieper low-Mach factor on the acoustic
!                     jump. MOSE_open never sets u2Roe there; here u2Roe = |u_roe|^2, the value of Riemann_MiczekRoe
!     13 MiczekRoe    MOSE_open Lib_Riemann_Roe  Riemann_MiczekRoe: the Miczek et al. 2015 preconditioned Roe flux
!     14 HLLC-PC      MOSE_open Lib_Riemann_HLL  Riemann_HLLCprec: HLLC Batten with the preconditioned signal velocities
!                     (urL_, urR_ = the Ur of the two cells of the face, integration-variables = prec)
! (MOSE_open b58ad91, statement by statement). bet is the shock sensor of the face, i.e. MOSE_open blk%beta of the face
! left cell (Tramel: 1 smooth, 0 at a shock; Chen: the pressure ratio h), 1 where MOSE_open passes 1 (symmetry). mco is
! riemann-options-Mco (AUSM+M). The thermodynamics are of FLINT, written out on the tables: f_Rtot, H / H0
! (f_tab: T = p/(R rho) clamped to [RO_TLO, RO_THI-1]), f_gamma (f_cp with the clamped f_tabT_expr indices), roe_averages
! (unclamped in MOSE_open; clamped here to the table bounds, a no-op for a state on the table).
! Included with #include after the twin has defined, for its own tables: RO_HTAB(i,s), RO_CPTAB(i,s), RO_RITAB(s), RO_WMTAB(s)
! and RO_TLO / RO_THI (FLINT Tmin / Tmax, the clamp bounds), and RO_MCO (riemann-options-Mco): the tables and Mco are read in
! place, never passed (module device arrays passed from the CUF kernels made NVVM redefine pgf90_dev_common), and RO_INV, the
! attribute of the five control scalars: value in the CUF twin (passed by reference from the CUF kernels, the literal and
! parameter actuals made NVVM redefine pgf90_dev_common in builds with more than 8 species), intent(in) in the OpenACC twin.
! Dummies: rs, bet, mco (read through RO_MCO by the OpenACC twin only), runiv_, ns_, the two states, the normal, the five
! fluxes.
! ---------------------------------------------------------------------------------------------------------------------
    implicit none
    integer,  RO_INV :: rs, ns_
    real(8),  RO_INV :: bet, mco, runiv_
    real(8),  intent(in)  :: dl(ns_), dr(ns_)
    real(8),  intent(in)  :: ul, vl, wl, pl, al, dltot, ur, vr, wr, pr, ar, drtot, nx, ny, nz
    real(8),  intent(in)  :: urL_, urR_   ! HLLC-PC (14): the reference velocities of the two cells
    ! HLLC-PC (14)
    real(8) :: pc_hL, pc_hR, pc_eL, pc_eR, pc_EtL, pc_EtR, pc_den, pc_VrL, pc_VrR, pc_VrRoe, pc_alL, pc_alR, pc_alRoe
    real(8) :: pc_vLs, pc_vRs, pc_vRoes, pc_cLs, pc_cRs, pc_cRoes
    real(8),  intent(out) :: F_r, F_u, F_v, F_w, F_E
    ! AUSM-family split constants (MOSE_open: beta = 1/8, alpha = 3/16)
    real(8), parameter :: ab_ = 0.125d0, aa_ = 0.1875d0
    real(8), parameter :: pi_ = 3.141592653589793d0
    integer :: s, Ti_, Ti1_, Ti2_, Til, Tir
    real(8) :: unl, unr, H_L, H_R, rho_, rg_, T_, Tc_, Td_, hs_, Hs_L, Hs_R, cp_, gam_L, gam_R, gam
    real(8) :: vel_
    ! AUSM+ / AUSM+M
    real(8) :: a_half, ML, MR, MplusL, MminusR, pplusL, pminusR, Mhalf, phalf, mdotp, mdotm
    real(8) :: Vl2, Vr2, hnorm, cs, c12, f_, g_, fo, rho12, Mp, M12, mdot, ps, gfac, pux, puy, puz
    ! Roe averages (roe_averages) and the HLL family
    real(8) :: Rl, Rr, inv_dltot, inv_drtot, Tl, Tr, dTl, dTr, srL, srR, inv_sr, rho_roe, u_roe, v_roe, w_roe
    real(8) :: R_roe, T_roe, h0l, h0r, cv_roe, sum_ei, invW, hl, hr, el, er, d_roe, h0_roe, gam_roe, vel2, a_roe, un_roe
    real(8) :: S1, S4, inv_denom, Sstar, pstar, U1S, U2S, U3S, U4S, U5S, e0l, e0r, FrL, FrR
    real(8) :: FrE, FuE, FvE, FwE, FeE, FrC, FuC, FvC, FwC, FeC
    real(8) :: lam1_hlle, lam2_hlle, lam3_hlle, eps_, lam1_roe, lam2_roe, lam3_roe, lam1_pp, lam2_pp, lam3_pp
    real(8) :: Ful, Fvl, Fwl, Fel, Fur, Fvr, Fwr, Fer
    real(8) :: drho, dp, du, dv, dw, dun, t1x, t1y, t1z, t2x, t2y, t2z, normt
    real(8) :: alpha_c, alpha_s1, alpha_s2, alpha_m, alpha_p, vv_
    ! LLF (8)
    real(8) :: lf_a, e0l_, e0r_
    ! SLAU2 (9)
    real(8) :: s2_af, s2_ml, s2_mr, s2_mt, s2_chi, s2_vnm, s2_frho, s2_vnl, s2_vnr, s2_bl, s2_br, s2_ph
    ! HLLC+Chen (10)
    real(8) :: ch_phiL, ch_phiR, ch_pref, ch_Mloc, ch_fM, ch_g, ch_fstar, ch_SK, ch_dU, ch_dutx, ch_duty, ch_dutz
    real(8) :: ch_lm, ch_sh
    ! exact (11)
    integer :: gx_iter
    real(8) :: gx_gm1, gx_gm4, gx_d1, gx_d4, gx_e1, gx_e4, gx_un1, gx_un4, gx_ut1, gx_vt1, gx_wt1, gx_ut4, gx_vt4
    real(8) :: gx_wt4, gx_r2i, gx_r1iv, gx_z, gx_p2, gx_p3, gx_u0, gx_vel, gx_dpr, gx_dpl, gx_a2, gx_a3, gx_vw1
    real(8) :: gx_vw2, gx_vw3, gx_vw4, gx_dum, gx_swl, gx_swr, gx_af, gx_pf, gx_rhof, gx_uf, gx_vf, gx_wf, gx_x
    real(8) :: gx_mrel, gx_c, gx_rp, gx_dum1, gx_delta, gx_df, gx_u
    ! LMRoe (12)
    real(8) :: lm_n1, lm_n2, lm_n3, lm_nmag, lm_tmpx, lm_tmpy, lm_tmpz, lm_tmag, lm_unL, lm_unR, lm_ut1L, lm_ut1R
    real(8) :: lm_ut2L, lm_ut2R, lm_hL, lm_hR, lm_H0L, lm_H0R, lm_FL1, lm_FL2, lm_FL3, lm_FL4, lm_FL5
    real(8) :: lm_FR1, lm_FR2, lm_FR3, lm_FR4, lm_FR5, lm_unRoe, lm_utmag, lm_Ma, lm_phi, lm_dunmod
    real(8) :: lm_abs1, lm_abs2, lm_abs3, lm_abs4, lm_abs5, lm_a2, lm_inva, lm_inva2, lm_al1, lm_al2, lm_al3, lm_al4
    real(8) :: lm_al5, lm_u2, lm_dut1, lm_dut2
    real(8) :: lm_r1_1, lm_r1_2, lm_r1_3, lm_r1_4, lm_r1_5, lm_r2_1, lm_r2_2, lm_r2_3, lm_r2_4, lm_r2_5
    real(8) :: lm_r3_1, lm_r3_2, lm_r3_3, lm_r3_4, lm_r3_5, lm_r4_1, lm_r4_2, lm_r4_3, lm_r4_4, lm_r4_5
    real(8) :: lm_r5_1, lm_r5_2, lm_r5_3, lm_r5_4, lm_r5_5, lm_d1, lm_d2, lm_d3, lm_d4, lm_d5
    ! MiczekRoe (13)
    real(8) :: mz_n1, mz_n2, mz_n3, mz_nmag, mz_tmpx, mz_tmpy, mz_tmpz, mz_tmag, mz_t1x, mz_t1y, mz_t1z
    real(8) :: mz_t2x, mz_t2y, mz_t2z, mz_unL, mz_unR, mz_ut1L, mz_ut1R, mz_ut2L, mz_ut2R, mz_hL, mz_hR, mz_H0L, mz_H0R
    real(8) :: mz_FL1, mz_FL2, mz_FL3, mz_FL4, mz_FL5, mz_FR1, mz_FR2, mz_FR3, mz_FR4, mz_FR5, mz_sRL, mz_sRR, mz_den
    real(8) :: mz_u, mz_v, mz_w, mz_un, mz_H, mz_rho, mz_u2, mz_kL, mz_kR, mz_kRoe, mz_chiL, mz_chiR, mz_chiRoe
    real(8) :: mz_cc, mz_Mcut, mz_Mloc, mz_mu, mz_del, mz_del2, mz_om, mz_p11, mz_p12, mz_p13, mz_p22, mz_p23
    real(8) :: mz_p32, mz_p33, mz_le, mz_lp, mz_lm, mz_s, mz_s2, mz_va1, mz_va2, mz_va3, mz_vb1, mz_vb2, mz_vb3
    real(8) :: mz_w1, mz_w2, mz_w3, mz_x1, mz_x2, mz_x3, mz_z1, mz_z2, mz_z3, mz_c1, mz_c2, mz_c3, mz_g1, mz_g2
    real(8) :: mz_g3, mz_ddrho, mz_ddun, mz_ddp, mz_ddut1, mz_ddut2, mz_ddux, mz_dduy, mz_dduz, mz_drho, mz_dp
    real(8) :: mz_dun, mz_dut1, mz_dut2

    unl = ul*nx + vl*ny + wl*nz
    unr = ur*nx + vr*ny + wr*nz

    select case (rs)

    ! ================================================================== AUSM+ (riemann_AUSMp)
    case (4)
      ! H_L = H0(pl, dl, sqrt(ul**2 + vl**2 + wl**2))
      vel_ = sqrt(ul**2 + vl**2 + wl**2)
      rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dl(s); end do
      rg_ = 0.d0; do s = 1, ns_; rg_ = rg_ + RO_RITAB(s)*dl(s)/rho_; end do
      Hs_L = 0.d0
      do s = 1, ns_
        T_ = pl/(rg_*rho_); Tc_ = max(dble(RO_TLO), min(T_, dble(RO_THI-1))); Ti_ = idint(Tc_); Td_ = Tc_ - Ti_
        Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
        hs_ = RO_HTAB(Ti1_,s) + (RO_HTAB(Ti2_,s) - RO_HTAB(Ti1_,s))*Td_
        Hs_L = Hs_L + hs_*dl(s)/rho_
      end do
      H_L = Hs_L + 0.5d0*vel_*vel_
      vel_ = sqrt(ur**2 + vr**2 + wr**2)
      rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dr(s); end do
      rg_ = 0.d0; do s = 1, ns_; rg_ = rg_ + RO_RITAB(s)*dr(s)/rho_; end do
      Hs_R = 0.d0
      do s = 1, ns_
        T_ = pr/(rg_*rho_); Tc_ = max(dble(RO_TLO), min(T_, dble(RO_THI-1))); Ti_ = idint(Tc_); Td_ = Tc_ - Ti_
        Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
        hs_ = RO_HTAB(Ti1_,s) + (RO_HTAB(Ti2_,s) - RO_HTAB(Ti1_,s))*Td_
        Hs_R = Hs_R + hs_*dr(s)/rho_
      end do
      H_R = Hs_R + 0.5d0*vel_*vel_
      a_half = sqrt( al * ar )
      ML = unl / a_half
      MR = unr / a_half
      if (abs(ML) < 1.d0) then
        MplusL = 0.25d0*(ML + 1.d0)**2 + ab_*(ML**2 - 1.d0)**2
        pplusL = 0.25d0*(ML + 1.d0)**2 * (2.d0 - ML) + aa_*ML*(ML**2 - 1.d0)**2
      else
        MplusL = 0.5d0*(ML + abs(ML))
        pplusL = 0.5d0*(1.d0 + sign(1.d0, ML))
      endif
      if (abs(MR) < 1.d0) then
        MminusR = -0.25d0*(MR - 1.d0)**2 - ab_*(MR**2 - 1.d0)**2
        pminusR = 0.25d0*(MR - 1.d0)**2 * (2.d0 + MR) - aa_*MR*(MR**2 - 1.d0)**2
      else
        MminusR = 0.5d0*(MR - abs(MR))
        pminusR = 0.5d0*(1.d0 - sign(1.d0, MR))
      endif
      Mhalf = MplusL + MminusR
      phalf = pplusL*pl + pminusR*pr
      mdotp = dltot * a_half * 0.5d0 * (Mhalf + abs(Mhalf))
      mdotm = drtot * a_half * 0.5d0 * (Mhalf - abs(Mhalf))
      F_r = mdotp + mdotm
      F_u = mdotp * ul + mdotm * ur + nx * phalf
      F_v = mdotp * vl + mdotm * vr + ny * phalf
      F_w = mdotp * wl + mdotm * wr + nz * phalf
      F_E = mdotp * H_L + mdotm * H_R

    ! ================================================================== AUSM+M (riemann_AUSMp_M)
    case (5)
      Vl2 = ul*ul + vl*vl + wl*wl
      Vr2 = ur*ur + vr*vr + wr*wr
      ! left: H (static), H0, f_gamma(dl, pl, dltot, f_Rtot(dl))
      rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dl(s); end do
      rg_ = 0.d0; do s = 1, ns_; rg_ = rg_ + RO_RITAB(s)*dl(s)/rho_; end do
      Hs_L = 0.d0
      do s = 1, ns_
        T_ = pl/(rg_*rho_); Tc_ = max(dble(RO_TLO), min(T_, dble(RO_THI-1))); Ti_ = idint(Tc_); Td_ = Tc_ - Ti_
        Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
        hs_ = RO_HTAB(Ti1_,s) + (RO_HTAB(Ti2_,s) - RO_HTAB(Ti1_,s))*Td_
        Hs_L = Hs_L + hs_*dl(s)/rho_
      end do
      vel_ = sqrt(Vl2); H_L = Hs_L + 0.5d0*vel_*vel_
      T_ = pl/(rg_*dltot); Ti_ = idint(T_); Td_ = T_ - Ti_
      Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
      cp_ = 0.d0
      do s = 1, ns_
        cp_ = cp_ + dl(s)/dltot*( RO_CPTAB(Ti1_,s) + (RO_CPTAB(Ti2_,s) - RO_CPTAB(Ti1_,s))*Td_ )
      end do
      gam_L = cp_/(cp_ - rg_)
      ! right
      rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dr(s); end do
      rg_ = 0.d0; do s = 1, ns_; rg_ = rg_ + RO_RITAB(s)*dr(s)/rho_; end do
      Hs_R = 0.d0
      do s = 1, ns_
        T_ = pr/(rg_*rho_); Tc_ = max(dble(RO_TLO), min(T_, dble(RO_THI-1))); Ti_ = idint(Tc_); Td_ = Tc_ - Ti_
        Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
        hs_ = RO_HTAB(Ti1_,s) + (RO_HTAB(Ti2_,s) - RO_HTAB(Ti1_,s))*Td_
        Hs_R = Hs_R + hs_*dr(s)/rho_
      end do
      vel_ = sqrt(Vr2); H_R = Hs_R + 0.5d0*vel_*vel_
      T_ = pr/(rg_*drtot); Ti_ = idint(T_); Td_ = T_ - Ti_
      Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
      cp_ = 0.d0
      do s = 1, ns_
        cp_ = cp_ + dr(s)/drtot*( RO_CPTAB(Ti1_,s) + (RO_CPTAB(Ti2_,s) - RO_CPTAB(Ti1_,s))*Td_ )
      end do
      gam_R = cp_/(cp_ - rg_)
      ! AUSMPW+ numerical sound speed (Eq. 28-29)
      gam   = 0.5d0*(gam_L + gam_R)
      hnorm = 0.5d0*( Hs_L + Hs_R )
      cs    = sqrt( 2.d0*(gam - 1.d0)/(gam + 1.d0)*hnorm )
      if (unl + unr >= 0.d0) then
        c12 = cs*cs/max(abs(unl), cs)
      else
        c12 = cs*cs/max(abs(unr), cs)
      end if
      ML = unl/c12 ; MR = unr/c12
      if (abs(ML) < 1.d0) then
        MplusL = 0.25d0*(ML + 1.d0)**2 + ab_*(ML*ML - 1.d0)**2
        pplusL = 0.25d0*(ML + 1.d0)**2*(2.d0 - ML) + aa_*ML*(ML*ML - 1.d0)**2
      else
        MplusL = 0.5d0*(ML + abs(ML))
        pplusL = 0.5d0*(1.d0 + sign(1.d0, ML))
      end if
      if (abs(MR) < 1.d0) then
        MminusR = -0.25d0*(MR - 1.d0)**2 - ab_*(MR*MR - 1.d0)**2
        pminusR = 0.25d0*(MR - 1.d0)**2*(2.d0 + MR) - aa_*MR*(MR*MR - 1.d0)**2
      else
        MminusR = 0.5d0*(MR - abs(MR))
        pminusR = 0.5d0*(1.d0 - sign(1.d0, MR))
      end if
      ! Mach-number function f and shock gate g (the sensor h of the face: switch = bet)
      f_ = 0.5d0*(1.d0 - cos(pi_*min(1.d0, max(abs(ML), abs(MR)))))
      g_ = 0.5d0*(1.d0 + cos(pi_*bet))
      rho12 = 0.5d0*(dltot + drtot)
      Mp    = -0.5d0*(1.d0 - f_)*(pr - pl)/(rho12*c12*c12)*(1.d0 - g_)
      M12   = MplusL + MminusR + Mp
      mdot  = M12*c12*merge(dltot, drtot, M12 >= 0.d0)
      fo = min(1.d0, max(f_, RO_MCO))
      ps = 0.5d0*(pl + pr) + 0.5d0*(pplusL - pminusR)*(pl - pr) &
         + fo*(pplusL + pminusR - 1.d0)*0.5d0*(pl + pr)
      gfac = g_*gam*(pl + pr)/(2.d0*c12)*pplusL*pminusR
      pux  = -gfac*(ur - ul)
      puy  = -gfac*(vr - vl)
      puz  = -gfac*(wr - wl)
      F_r = mdot
      if (M12 >= 0.d0) then
        F_u = mdot*ul + ps*nx + pux
        F_v = mdot*vl + ps*ny + puy
        F_w = mdot*wl + ps*nz + puz
        F_E = mdot*H_L
      else
        F_u = mdot*ur + ps*nx + pux
        F_v = mdot*vr + ps*ny + puy
        F_w = mdot*wr + ps*nz + puz
        F_E = mdot*H_R
      end if


    ! ================================================================== LLF / Rusanov (riemann_LLF)
    case (8)
      ! A = MAX(ABS(uln-al),ABS(urn-ar),ABS(uln+al),ABS(urn+ar)); fluxes(p,r,...): F_r = sum(r)*un, F_E = F_r*H0(p,r,|u|);
      ! E0(p,r,|u|) = (H(p,r) - p/sum(r)) + |u|^2/2
      lf_a = MAX(ABS(unl-al),ABS(unr-ar),ABS(unl+al),ABS(unr+ar))
      rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dr(s); end do
      rg_ = 0.d0; do s = 1, ns_; rg_ = rg_ + RO_RITAB(s)*dr(s)/rho_; end do
      Hs_R = 0.d0
      do s = 1, ns_
        T_ = pr/(rg_*rho_); Tc_ = max(dble(RO_TLO), min(T_, dble(RO_THI-1))); Ti_ = idint(Tc_); Td_ = Tc_ - Ti_
        Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
        hs_ = RO_HTAB(Ti1_,s) + (RO_HTAB(Ti2_,s) - RO_HTAB(Ti1_,s))*Td_
        Hs_R = Hs_R + hs_*dr(s)/rho_
      end do
      vel_ = sqrt(ur*ur+vr*vr+wr*wr)
      FrR = rho_*(ur*nx+vr*ny+wr*nz)
      Fur = FrR*ur + pr*nx; Fvr = FrR*vr + pr*ny; Fwr = FrR*wr + pr*nz; Fer = FrR*(Hs_R + 0.5d0*vel_*vel_)
      e0r_ = (Hs_R - pr/rho_) + 0.5d0*vel_*vel_
      rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dl(s); end do
      rg_ = 0.d0; do s = 1, ns_; rg_ = rg_ + RO_RITAB(s)*dl(s)/rho_; end do
      Hs_L = 0.d0
      do s = 1, ns_
        T_ = pl/(rg_*rho_); Tc_ = max(dble(RO_TLO), min(T_, dble(RO_THI-1))); Ti_ = idint(Tc_); Td_ = Tc_ - Ti_
        Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
        hs_ = RO_HTAB(Ti1_,s) + (RO_HTAB(Ti2_,s) - RO_HTAB(Ti1_,s))*Td_
        Hs_L = Hs_L + hs_*dl(s)/rho_
      end do
      vel_ = sqrt(ul*ul+vl*vl+wl*wl)
      FrL = rho_*(ul*nx+vl*ny+wl*nz)
      Ful = FrL*ul + pl*nx; Fvl = FrL*vl + pl*ny; Fwl = FrL*wl + pl*nz; Fel = FrL*(Hs_L + 0.5d0*vel_*vel_)
      e0l_ = (Hs_L - pl/rho_) + 0.5d0*vel_*vel_
      F_r = 0.5*(FrL + FrR - lf_a*(drtot - dltot))
      F_u = 0.5*(Fur + Ful - lf_a*(drtot*ur-dltot*ul))
      F_v = 0.5*(Fvr + Fvl - lf_a*(drtot*vr-dltot*vl))
      F_w = 0.5*(Fwr + Fwl - lf_a*(drtot*wr-dltot*wl))
      F_E = 0.5*(Fel + Fer - lf_a*(drtot*e0r_ - dltot*e0l_))

    ! ================================================================== SLAU2 (riemann_SLAU_base + riemann_SLAU2)
    case (9)
      Vl2 = ul**2 + vl**2 + wl**2
      Vr2 = ur**2 + vr**2 + wr**2
      s2_af = 0.5d0*(al + ar)
      s2_ml = unl / s2_af
      s2_mr = unr / s2_af
      s2_mt = min(1.d0, sqrt(0.5d0*(Vl2 + Vr2)) / s2_af)
      s2_chi = (1.d0 - s2_mt)**2
      s2_vnm = (dltot*abs(unl) + drtot*abs(unr)) / (dltot + drtot)
      s2_frho = -max(min(s2_ml,0.d0), -1.d0) * min(max(s2_mr,0.d0), 1.d0)
      s2_vnl = (1.d0 - s2_frho)*s2_vnm + s2_frho*abs(unl)
      s2_vnr = (1.d0 - s2_frho)*s2_vnm + s2_frho*abs(unr)
      F_r = 0.5d0 * ( dltot*(unl + s2_vnl) + drtot*(unr - s2_vnr) - (s2_chi/s2_af)*(pr - pl) )
      F_r = max(-1.0d150, min(F_r, 1.0d150))
      if (abs(s2_ml) < 1.d0) then
        s2_bl = 0.25d0*(2.d0 - s2_ml)*(s2_ml + 1.d0)**2
      else if (s2_ml >= 0.d0) then
        s2_bl = 1.d0
      else
        s2_bl = 0.d0
      endif
      if (abs(s2_mr) < 1.d0) then
        s2_br = 0.25d0*(2.d0 + s2_mr)*(s2_mr - 1.d0)**2
      else if (s2_mr >= 0.d0) then
        s2_br = 0.d0
      else
        s2_br = 1.d0
      endif
      ! SLAU2 pressure flux (Kitamura and Shima 2013)
      s2_ph = 0.5d0*(pl + pr) + 0.5d0*(s2_bl - s2_br)*(pl - pr) &
            + sqrt(0.5d0*(Vl2 + Vr2))*(s2_bl + s2_br - 1.d0)*s2_af*0.5d0*(dltot + drtot)
      s2_ph = max(-1.0d150, min(s2_ph, 1.0d150))
      rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dl(s); end do
      rg_ = 0.d0; do s = 1, ns_; rg_ = rg_ + RO_RITAB(s)*dl(s)/rho_; end do
      Hs_L = 0.d0
      do s = 1, ns_
        T_ = pl/(rg_*rho_); Tc_ = max(dble(RO_TLO), min(T_, dble(RO_THI-1))); Ti_ = idint(Tc_); Td_ = Tc_ - Ti_
        Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
        hs_ = RO_HTAB(Ti1_,s) + (RO_HTAB(Ti2_,s) - RO_HTAB(Ti1_,s))*Td_
        Hs_L = Hs_L + hs_*dl(s)/rho_
      end do
      rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dr(s); end do
      rg_ = 0.d0; do s = 1, ns_; rg_ = rg_ + RO_RITAB(s)*dr(s)/rho_; end do
      Hs_R = 0.d0
      do s = 1, ns_
        T_ = pr/(rg_*rho_); Tc_ = max(dble(RO_TLO), min(T_, dble(RO_THI-1))); Ti_ = idint(Tc_); Td_ = Tc_ - Ti_
        Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
        hs_ = RO_HTAB(Ti1_,s) + (RO_HTAB(Ti2_,s) - RO_HTAB(Ti1_,s))*Td_
        Hs_R = Hs_R + hs_*dr(s)/rho_
      end do
      H_L = Hs_L + 0.5d0*sqrt(Vl2)*sqrt(Vl2)
      H_R = Hs_R + 0.5d0*sqrt(Vr2)*sqrt(Vr2)
      H_L = max(-1.0d150, min(H_L, 1.0d150))
      H_R = max(-1.0d150, min(H_R, 1.0d150))
      if (F_r >= 0.d0) then
        F_u = F_r*ul + nx*s2_ph
        F_v = F_r*vl + ny*s2_ph
        F_w = F_r*wl + nz*s2_ph
        F_E = F_r*H_L
      else
        F_u = F_r*ur + nx*s2_ph
        F_v = F_r*vr + ny*s2_ph
        F_w = F_r*wr + nz*s2_ph
        F_E = F_r*H_R
      endif

    ! ================================================================== exact (riemann_exact, Godunov)
    case (11)
      rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dl(s); end do
      rg_ = 0.d0; do s = 1, ns_; rg_ = rg_ + RO_RITAB(s)*dl(s)/rho_; end do
      T_ = pl/(rg_*dltot); Ti_ = idint(T_); Td_ = T_ - Ti_
      Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
      cp_ = 0.d0
      do s = 1, ns_
        cp_ = cp_ + dl(s)/dltot*( RO_CPTAB(Ti1_,s) + (RO_CPTAB(Ti2_,s) - RO_CPTAB(Ti1_,s))*Td_ )
      end do
      gam_L = cp_/(cp_ - rg_)
      rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dr(s); end do
      rg_ = 0.d0; do s = 1, ns_; rg_ = rg_ + RO_RITAB(s)*dr(s)/rho_; end do
      T_ = pr/(rg_*drtot); Ti_ = idint(T_); Td_ = T_ - Ti_
      Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
      cp_ = 0.d0
      do s = 1, ns_
        cp_ = cp_ + dr(s)/drtot*( RO_CPTAB(Ti1_,s) + (RO_CPTAB(Ti2_,s) - RO_CPTAB(Ti1_,s))*Td_ )
      end do
      gam_R = cp_/(cp_ - rg_)
      gx_gm1 = gam_L-1.d0
      gx_gm4 = gam_R-1.d0
      gx_d1 = 0.5d0*gx_gm1
      gx_d4 = 0.5d0*gx_gm4
      gx_e1 = gam_L/gx_d1
      gx_e4 = gam_R/gx_d4
      gx_un1 = ul*nx+vl*ny+wl*nz
      gx_ut1 = ul-gx_un1*nx; gx_vt1 = vl-gx_un1*ny; gx_wt1 = wl-gx_un1*nz
      gx_un4 = ur*nx+vr*ny+wr*nz
      gx_ut4 = ur-gx_un4*nx; gx_vt4 = vr-gx_un4*ny; gx_wt4 = wr-gx_un4*nz
      gx_r2i  = gx_un1 + (2.d0/gx_gm1)*al
      gx_r1iv = gx_un4 - (2.d0/gx_gm4)*ar
      if (pl<pr) then
        gx_dum = gx_gm4/(2.d0*gam_R)
      else
        gx_dum = gx_gm1/(2.d0*gam_L)
      end if
      gx_z = (gx_gm1/gx_gm4)*(ar/al)*((pl/pr)**gx_dum)
      gx_u0 = (gx_z*gx_r2i+gx_r1iv)/(1.d0+gx_z)
      gx_vel = gx_u0
      gx_p3 = 1d0
      gx_p2 = 2d0*gx_p3
      gx_iter = 0
      ! MOSE_open leaves the wave speeds of a branch never taken unset; here they start at 0 and the sound speeds of the
      ! states (they are read only on a branch the iteration took)
      gx_swl = 0d0; gx_swr = 0d0; gx_vw1 = 0d0; gx_vw2 = 0d0; gx_vw3 = 0d0; gx_vw4 = 0d0; gx_a2 = al; gx_a3 = ar
      do while (dabs(1-(gx_p2/gx_p3))>=1d-5 .and. gx_iter<=1000)
        if (gx_vel<gx_un1) then
          ! shockleft(un1,a1,p1,vel,a2,p2,dpr,swl,gam1)
          gx_x = (gam_L+1.d00)*(gx_vel-gx_un1)/(4.d00*al)
          gx_mrel = gx_x-dsqrt(1.d00+gx_x**2.d00)
          gx_c = gam_L*pl/al
          gx_p2 = pl+gx_c*(gx_vel-gx_un1)*gx_mrel
          gx_dpr = 2.d00*gx_c*gx_mrel**3.d00/(1.d00+gx_mrel**2.d00)
          gx_rp = gx_p2/pl
          gx_swl = gx_un1+al*gx_mrel
          gx_a2 = al*dsqrt((gam_L+1.d00+(gam_L-1.d00)*gx_rp)/(gam_L+1.d00+(gam_L-1.d00)/gx_rp))
        else
          ! rareleft(p1,a1,un1,vel,p2,a2,dpr,vw1,vw2,gam1)
          gx_delta = .5d0*(gam_L-1)
          gx_dum1 = 2.d00*gam_L/(gam_L-1.d0)
          gx_a2 = al-gx_delta*(gx_vel-gx_un1)
          gx_p2 = pl*((gx_a2/al)**gx_dum1)
          gx_dpr = -gam_L*gx_p2/gx_a2
          gx_vw1 = gx_un1-al
          gx_vw2 = gx_vel-gx_a2
        endif
        if (gx_vel>gx_un4) then
          ! shockright(un4,a4,p4,vel,a3,p3,dpl,swr,gam4)
          gx_x = (gam_R+1.d00)*(gx_vel-gx_un4)/(4.d00*ar)
          gx_mrel = gx_x+dsqrt(1.d00+gx_x**2.d00)
          gx_c = gam_R*pr/ar
          gx_p3 = pr+gx_c*(gx_vel-gx_un4)*gx_mrel
          gx_dpl = 2.d00*gx_c*gx_mrel**3.d00/(1.d00+gx_mrel**2.d00)
          gx_rp = gx_p3/pr
          gx_swr = gx_un4+ar*gx_mrel
          gx_a3 = ar*dsqrt((gam_R+1.d00+(gam_R-1.d00)*gx_rp)/(gam_R+1.d00+(gam_R-1.d00)/gx_rp))
        else
          ! rareright(p4,a4,un4,vel,p3,a3,dpl,vw3,vw4,gam4)
          gx_delta = .5d0*(gam_R-1)
          gx_dum1 = 2.d00*gam_R/(gam_R-1.d0)
          gx_a3 = ar+gx_delta*(gx_vel-gx_un4)
          gx_p3 = pr*((gx_a3/ar)**gx_dum1)
          gx_dpl = gam_R*gx_p3/gx_a3
          gx_vw3 = gx_vel+gx_a3
          gx_vw4 = gx_un4+ar
        endif
        gx_vel = gx_vel-((gx_p2-gx_p3)/(gx_dpr-gx_dpl))
        gx_iter = gx_iter+1
      end do
      ! the fluxes at the interface
      if (gx_vel>0.d0) then
        gx_af = gx_a2
        gx_pf = gx_p2
        if (gx_vel<=gx_un1) then
          if (gx_swl>0.d0) then
            gx_af = al
            gx_vel = gx_un1
            gx_pf = pl
          end if
        else
          if (gx_vw1>0.d0) then
            gx_vel = gx_un1
            gx_af = al
            gx_pf = pl
          else
            if (gx_vw2>0.d0) then
              gx_af = (gx_un1+al/gx_d1)/(1.d0+1.d0/gx_d1)
              gx_vel = gx_af
              gx_pf = pl*(gx_af/al)**gx_e1
            end if
          end if
        end if
        gx_rhof = gam_L*gx_pf/(gx_af**2.d0)
        gx_uf = gx_vel*nx+gx_ut1; gx_vf = gx_vel*ny+gx_vt1; gx_wf = gx_vel*nz+gx_wt1
        ! H(pf, df) with df = dl/dltot*rhof
        rho_ = 0.d0; do s = 1, ns_; gx_df = dl(s)/dltot*gx_rhof; rho_ = rho_ + gx_df; end do
        rg_ = 0.d0; do s = 1, ns_; gx_df = dl(s)/dltot*gx_rhof; rg_ = rg_ + RO_RITAB(s)*gx_df/rho_; end do
        Hs_L = 0.d0
        do s = 1, ns_
          gx_df = dl(s)/dltot*gx_rhof
          T_ = gx_pf/(rg_*rho_); Tc_ = max(dble(RO_TLO), min(T_, dble(RO_THI-1))); Ti_ = idint(Tc_); Td_ = Tc_ - Ti_
          Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
          hs_ = RO_HTAB(Ti1_,s) + (RO_HTAB(Ti2_,s) - RO_HTAB(Ti1_,s))*Td_
          Hs_L = Hs_L + hs_*gx_df/rho_
        end do
      else
        gx_af = gx_a3
        gx_pf = gx_p3
        if (gx_vel>=gx_un4) then
          if (gx_swr<0.d0) then
            gx_af = ar
            gx_vel = gx_un4
            gx_pf = pr
          end if
        else
          if (gx_vw4<0.d0) then
            gx_af = ar
            gx_vel = gx_un4
            gx_pf = pr
          else
            if (gx_vw3<0.d0) then
              gx_af = -(gx_un4-ar/gx_d4)/(1.d0+1.d0/gx_d4)
              gx_vel = -gx_af
              gx_pf = pr*(gx_af/ar)**gx_e4
            end if
          end if
        end if
        gx_rhof = gam_R*gx_pf/(gx_af**2.d0)
        gx_uf = gx_vel*nx+gx_ut4; gx_vf = gx_vel*ny+gx_vt4; gx_wf = gx_vel*nz+gx_wt4
        ! H(pf, df) with df = dr/drtot*rhof
        rho_ = 0.d0; do s = 1, ns_; gx_df = dr(s)/drtot*gx_rhof; rho_ = rho_ + gx_df; end do
        rg_ = 0.d0; do s = 1, ns_; gx_df = dr(s)/drtot*gx_rhof; rg_ = rg_ + RO_RITAB(s)*gx_df/rho_; end do
        Hs_L = 0.d0
        do s = 1, ns_
          gx_df = dr(s)/drtot*gx_rhof
          T_ = gx_pf/(rg_*rho_); Tc_ = max(dble(RO_TLO), min(T_, dble(RO_THI-1))); Ti_ = idint(Tc_); Td_ = Tc_ - Ti_
          Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
          hs_ = RO_HTAB(Ti1_,s) + (RO_HTAB(Ti2_,s) - RO_HTAB(Ti1_,s))*Td_
          Hs_L = Hs_L + hs_*gx_df/rho_
        end do
      end if
      F_r = gx_rhof*gx_vel
      F_u = F_r*gx_uf + gx_pf*nx
      F_v = F_r*gx_vf + gx_pf*ny
      F_w = F_r*gx_wf + gx_pf*nz
      ! F_E = F_r*H0(pf,df,sqrt(vel*vel)): the face enthalpy with the normal velocity only, as MOSE_open
      gx_u = sqrt(gx_vel*gx_vel)
      F_E = F_r*(Hs_L + 0.5d0*gx_u*gx_u)

    ! ================================================================== MiczekRoe (Riemann_MiczekRoe)
    case (13)
      mz_nmag = sqrt(nx*nx + ny*ny + nz*nz)
      if (mz_nmag <= 0.d0) then
        F_r = 0.d0; F_u = 0.d0; F_v = 0.d0; F_w = 0.d0; F_E = 0.d0
      else
      mz_n1 = nx/mz_nmag; mz_n2 = ny/mz_nmag; mz_n3 = nz/mz_nmag
      if (abs(mz_n1) < 0.9d0) then
        mz_tmpx = 1.d0; mz_tmpy = 0.d0; mz_tmpz = 0.d0
      else
        mz_tmpx = 0.d0; mz_tmpy = 1.d0; mz_tmpz = 0.d0
      end if
      mz_t1x = mz_tmpx - (mz_tmpx*mz_n1 + mz_tmpy*mz_n2 + mz_tmpz*mz_n3)*mz_n1
      mz_t1y = mz_tmpy - (mz_tmpx*mz_n1 + mz_tmpy*mz_n2 + mz_tmpz*mz_n3)*mz_n2
      mz_t1z = mz_tmpz - (mz_tmpx*mz_n1 + mz_tmpy*mz_n2 + mz_tmpz*mz_n3)*mz_n3
      mz_tmag = sqrt(mz_t1x*mz_t1x + mz_t1y*mz_t1y + mz_t1z*mz_t1z)
      if (mz_tmag <= 0.d0) then
        mz_t1x = -mz_n2; mz_t1y =  mz_n1; mz_t1z = 0.d0
        mz_tmag = sqrt(mz_t1x*mz_t1x + mz_t1y*mz_t1y + mz_t1z*mz_t1z)
      end if
      mz_t1x = mz_t1x/mz_tmag; mz_t1y = mz_t1y/mz_tmag; mz_t1z = mz_t1z/mz_tmag
      mz_t2x = mz_n2*mz_t1z - mz_n3*mz_t1y
      mz_t2y = mz_n3*mz_t1x - mz_n1*mz_t1z
      mz_t2z = mz_n1*mz_t1y - mz_n2*mz_t1x
      mz_unL  = ul*mz_n1 + vl*mz_n2 + wl*mz_n3
      mz_unR  = ur*mz_n1 + vr*mz_n2 + wr*mz_n3
      mz_ut1L = ul*mz_t1x + vl*mz_t1y + wl*mz_t1z
      mz_ut1R = ur*mz_t1x + vr*mz_t1y + wr*mz_t1z
      mz_ut2L = ul*mz_t2x + vl*mz_t2y + wl*mz_t2z
      mz_ut2R = ur*mz_t2x + vr*mz_t2y + wr*mz_t2z
      ! hL = H(pl,dl), hR = H(pr,dr)
      rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dl(s); end do
      rg_ = 0.d0; do s = 1, ns_; rg_ = rg_ + RO_RITAB(s)*dl(s)/rho_; end do
      mz_hL = 0.d0
      do s = 1, ns_
        T_ = pl/(rg_*rho_); Tc_ = max(dble(RO_TLO), min(T_, dble(RO_THI-1))); Ti_ = idint(Tc_); Td_ = Tc_ - Ti_
        Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
        hs_ = RO_HTAB(Ti1_,s) + (RO_HTAB(Ti2_,s) - RO_HTAB(Ti1_,s))*Td_
        mz_hL = mz_hL + hs_*dl(s)/rho_
      end do
      rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dr(s); end do
      rg_ = 0.d0; do s = 1, ns_; rg_ = rg_ + RO_RITAB(s)*dr(s)/rho_; end do
      mz_hR = 0.d0
      do s = 1, ns_
        T_ = pr/(rg_*rho_); Tc_ = max(dble(RO_TLO), min(T_, dble(RO_THI-1))); Ti_ = idint(Tc_); Td_ = Tc_ - Ti_
        Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
        hs_ = RO_HTAB(Ti1_,s) + (RO_HTAB(Ti2_,s) - RO_HTAB(Ti1_,s))*Td_
        mz_hR = mz_hR + hs_*dr(s)/rho_
      end do
      mz_H0L = mz_hL + 0.5d0*(ul*ul + vl*vl + wl*wl)
      mz_H0R = mz_hR + 0.5d0*(ur*ur + vr*vr + wr*wr)
      mz_FL1 = dltot*mz_unL
      mz_FL2 = dltot*mz_unL*ul + pl*mz_n1
      mz_FL3 = dltot*mz_unL*vl + pl*mz_n2
      mz_FL4 = dltot*mz_unL*wl + pl*mz_n3
      mz_FL5 = dltot*mz_unL*mz_H0L
      mz_FR1 = drtot*mz_unR
      mz_FR2 = drtot*mz_unR*ur + pr*mz_n1
      mz_FR3 = drtot*mz_unR*vr + pr*mz_n2
      mz_FR4 = drtot*mz_unR*wr + pr*mz_n3
      mz_FR5 = drtot*mz_unR*mz_H0R
      mz_sRL = sqrt(dltot)
      mz_sRR = sqrt(drtot)
      mz_den = mz_sRL + mz_sRR
      mz_u  = (mz_sRL*ul + mz_sRR*ur)/mz_den
      mz_v  = (mz_sRL*vl + mz_sRR*vr)/mz_den
      mz_w  = (mz_sRL*wl + mz_sRR*wr)/mz_den
      mz_H  = (mz_sRL*mz_H0L + mz_sRR*mz_H0R)/mz_den
      mz_un = mz_u*mz_n1 + mz_v*mz_n2 + mz_w*mz_n3
      mz_rho = mz_sRL*mz_sRR
      mz_u2  = mz_u*mz_u + mz_v*mz_v + mz_w*mz_w
      ! kappa = f_gamma - 1 of each state, chi = a^2 - kappa h
      rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dl(s); end do
      rg_ = 0.d0; do s = 1, ns_; rg_ = rg_ + RO_RITAB(s)*dl(s)/rho_; end do
      T_ = pl/(rg_*dltot); Ti_ = idint(T_); Td_ = T_ - Ti_
      Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
      cp_ = 0.d0
      do s = 1, ns_
        cp_ = cp_ + dl(s)/dltot*( RO_CPTAB(Ti1_,s) + (RO_CPTAB(Ti2_,s) - RO_CPTAB(Ti1_,s))*Td_ )
      end do
      gam_L = cp_/(cp_ - rg_)
      rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dr(s); end do
      rg_ = 0.d0; do s = 1, ns_; rg_ = rg_ + RO_RITAB(s)*dr(s)/rho_; end do
      T_ = pr/(rg_*drtot); Ti_ = idint(T_); Td_ = T_ - Ti_
      Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
      cp_ = 0.d0
      do s = 1, ns_
        cp_ = cp_ + dr(s)/drtot*( RO_CPTAB(Ti1_,s) + (RO_CPTAB(Ti2_,s) - RO_CPTAB(Ti1_,s))*Td_ )
      end do
      gam_R = cp_/(cp_ - rg_)
      mz_kL = gam_L - 1.d0
      mz_kR = gam_R - 1.d0
      mz_chiL = al*al - mz_kL*mz_hL
      mz_chiR = ar*ar - mz_kR*mz_hR
      mz_kRoe = (mz_sRL*mz_kL + mz_sRR*mz_kR)/mz_den
      mz_chiRoe = (mz_sRL*mz_chiL + mz_sRR*mz_chiR)/mz_den
      mz_cc = sqrt( max(1.d-30, mz_chiRoe + mz_kRoe*(mz_H - 0.5d0*mz_u2)) )
      mz_Mcut = max(1.d-6, RO_MCO)
      mz_Mloc = sqrt(mz_u2)/mz_cc
      mz_mu   = min(1.d0, max(mz_Mloc, mz_Mcut))
      mz_del  = 1.d0/mz_mu - 1.d0
      mz_del2 = mz_del*mz_del
      mz_om   = 1.d0/(1.d0 + mz_del2)
      mz_drho = drtot - dltot
      mz_dp   = pr - pl
      mz_dun  = mz_unR - mz_unL
      mz_dut1 = mz_ut1R - mz_ut1L
      mz_dut2 = mz_ut2R - mz_ut2L
      mz_p11 = mz_un
      mz_p12 = mz_rho*(1.d0 + mz_del*mz_un/mz_cc)
      mz_p13 = mz_del/mz_cc
      mz_p22 = mz_un - mz_del*mz_cc
      mz_p23 = (1.d0 - mz_del*mz_un/mz_cc)/mz_rho
      mz_p32 = mz_rho*mz_cc*(mz_cc + mz_del*mz_un)
      mz_p33 = mz_un + mz_cc*mz_del
      mz_s2 = (1.d0 + mz_del2)*mz_cc*mz_cc - mz_del2*mz_un*mz_un
      mz_s  = sqrt( max(1.d-30, mz_s2) )
      mz_s2 = mz_s*mz_s
      mz_le = mz_un
      mz_lp = mz_un + mz_s
      mz_lm = mz_un - mz_s
      ! pa3_shift(lam, x, r): r1 = (p11-lam) x1 + p12 x2 + p13 x3, r2 = (p22-lam) x2 + p23 x3, r3 = p32 x2 + (p33-lam) x3
      mz_va1 = (mz_p11 - mz_lm)*mz_drho + mz_p12*mz_dun + mz_p13*mz_dp
      mz_va2 = (mz_p22 - mz_lm)*mz_dun + mz_p23*mz_dp
      mz_va3 = mz_p32*mz_dun + (mz_p33 - mz_lm)*mz_dp
      mz_vb1 = (mz_p11 - mz_lp)*mz_drho + mz_p12*mz_dun + mz_p13*mz_dp
      mz_vb2 = (mz_p22 - mz_lp)*mz_dun + mz_p23*mz_dp
      mz_vb3 = mz_p32*mz_dun + (mz_p33 - mz_lp)*mz_dp
      mz_w1 = (mz_p11 - mz_lp)*mz_va1 + mz_p12*mz_va2 + mz_p13*mz_va3
      mz_w2 = (mz_p22 - mz_lp)*mz_va2 + mz_p23*mz_va3
      mz_w3 = mz_p32*mz_va2 + (mz_p33 - mz_lp)*mz_va3
      mz_x1 = (mz_p11 - mz_le)*mz_va1 + mz_p12*mz_va2 + mz_p13*mz_va3
      mz_x2 = (mz_p22 - mz_le)*mz_va2 + mz_p23*mz_va3
      mz_x3 = mz_p32*mz_va2 + (mz_p33 - mz_le)*mz_va3
      mz_z1 = (mz_p11 - mz_le)*mz_vb1 + mz_p12*mz_vb2 + mz_p13*mz_vb3
      mz_z2 = (mz_p22 - mz_le)*mz_vb2 + mz_p23*mz_vb3
      mz_z3 = mz_p32*mz_vb2 + (mz_p33 - mz_le)*mz_vb3
      mz_c1 = abs(mz_le)/(-mz_s2)
      mz_c2 = abs(mz_lp)/(2.d0*mz_s2)
      mz_c3 = abs(mz_lm)/(2.d0*mz_s2)
      mz_g1 = mz_c1*mz_w1 + mz_c2*mz_x1 + mz_c3*mz_z1
      mz_g2 = mz_c1*mz_w2 + mz_c2*mz_x2 + mz_c3*mz_z2
      mz_g3 = mz_c1*mz_w3 + mz_c2*mz_x3 + mz_c3*mz_z3
      mz_ddrho = mz_g1 - (mz_rho*mz_del/mz_cc)*mz_om*mz_g2 - (mz_del2/(mz_cc*mz_cc))*mz_om*mz_g3
      mz_ddun  = mz_om*( mz_g2 + (mz_del/(mz_rho*mz_cc))*mz_g3 )
      mz_ddp   = mz_om*( -mz_rho*mz_cc*mz_del*mz_g2 + mz_g3 )
      mz_ddut1 = abs(mz_un)*mz_dut1
      mz_ddut2 = abs(mz_un)*mz_dut2
      mz_ddux = mz_ddun*mz_n1 + mz_ddut1*mz_t1x + mz_ddut2*mz_t2x
      mz_dduy = mz_ddun*mz_n2 + mz_ddut1*mz_t1y + mz_ddut2*mz_t2y
      mz_dduz = mz_ddun*mz_n3 + mz_ddut1*mz_t1z + mz_ddut2*mz_t2z
      F_r = 0.5d0*(mz_FL1 + mz_FR1) - 0.5d0*( mz_ddrho )
      F_u = 0.5d0*(mz_FL2 + mz_FR2) - 0.5d0*( mz_u*mz_ddrho + mz_rho*mz_ddux )
      F_v = 0.5d0*(mz_FL3 + mz_FR3) - 0.5d0*( mz_v*mz_ddrho + mz_rho*mz_dduy )
      F_w = 0.5d0*(mz_FL4 + mz_FR4) - 0.5d0*( mz_w*mz_ddrho + mz_rho*mz_dduz )
      F_E = 0.5d0*(mz_FL5 + mz_FR5) - 0.5d0*( 0.5d0*mz_u2*mz_ddrho &
                                      + mz_rho*(mz_u*mz_ddux + mz_v*mz_dduy + mz_w*mz_dduz) &
                                      + mz_ddp/mz_kRoe )
      end if

    ! ================================================================== HLLE++ (6), HLLC+Tramel (7), HLLC+Chen (10)
    ! and LMRoe (12): the Roe averages of MOSE_open roe_averages first
    case (6, 7, 10, 12)
      ! roe_averages (MOSE_open Lib_Riemann_HLL): Rl = f_Rtot(dl), Rr = f_Rtot(dr)
      rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dl(s); end do
      Rl = 0.d0; do s = 1, ns_; Rl = Rl + RO_RITAB(s)*dl(s)/rho_; end do
      rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dr(s); end do
      Rr = 0.d0; do s = 1, ns_; Rr = Rr + RO_RITAB(s)*dr(s)/rho_; end do
      inv_dltot = 1.d0 / dltot
      inv_drtot = 1.d0 / drtot
      Tl = pl * inv_dltot / Rl
      Tr = pr * inv_drtot / Rr
      Til = int(Tl); Til = max(RO_TLO, min(Til, RO_THI-1))
      Tir = int(Tr); Tir = max(RO_TLO, min(Tir, RO_THI-1))
      dTl = Tl - Til
      dTr = Tr - Tir
      srL = sqrt(dltot)
      srR = sqrt(drtot)
      inv_sr = 1.d0 / (srL + srR)
      rho_roe = srL * srR
      u_roe = (srR*ur + srL*ul) * inv_sr
      v_roe = (srR*vr + srL*vl) * inv_sr
      w_roe = (srR*wr + srL*wl) * inv_sr
      R_roe = (srR*Rr + srL*Rl) * inv_sr
      T_roe = (srR*Tr + srL*Tl) * inv_sr
      h0l = 0.d0; h0r = 0.d0; cv_roe = 0.d0; sum_ei = 0.d0
      do s = 1, ns_
        invW = runiv_ / RO_WMTAB(s)
        hl = RO_HTAB(Til, s) + (RO_HTAB(Til+1, s) - RO_HTAB(Til, s)) * dTl
        hr = RO_HTAB(Tir, s) + (RO_HTAB(Tir+1, s) - RO_HTAB(Tir, s)) * dTr
        el = hl - invW * Tl
        er = hr - invW * Tr
        h0l = h0l + hl * dl(s) * inv_dltot
        h0r = h0r + hr * dr(s) * inv_drtot
        d_roe = (srR*dr(s)*inv_drtot + srL*dl(s)*inv_dltot) * inv_sr
        cv_roe = cv_roe + d_roe * ( &
            0.5d0 * ( &
              RO_CPTAB(Til, s) + (RO_CPTAB(Til+1, s)-RO_CPTAB(Til, s))*dTl + &
              RO_CPTAB(Tir, s) + (RO_CPTAB(Tir+1, s)-RO_CPTAB(Tir, s))*dTr ) &
            - invW )
        sum_ei = sum_ei + d_roe * (srR*er + srL*el) * inv_sr
      end do
      h0l = h0l + 0.5d0*(ul*ul + vl*vl + wl*wl)
      h0r = h0r + 0.5d0*(ur*ur + vr*vr + wr*wr)
      h0_roe = (srR*h0r + srL*h0l) * inv_sr
      gam_roe = 1.d0 + R_roe / cv_roe
      vel2 = u_roe*u_roe + v_roe*v_roe + w_roe*w_roe
      a_roe = sqrt( (gam_roe - 1.d0) * &
            ( h0_roe - 0.5d0*vel2 + cv_roe*T_roe - sum_ei ) )
      un_roe = u_roe*nx + v_roe*ny + w_roe*nz

      if (rs == 6) then
        ! ---------------------------------------------------------------- HLLE++ (riemann_HLLEpp)
        S1 = min(0.d0, unl-al, un_roe-a_roe)
        S4 = max(0.d0, unr+ar, un_roe+a_roe)
        inv_denom = 1d0 / (S4 - S1)
        lam1_hlle = abs( (S4+S1)*un_roe*inv_denom - 2d0*S4*S1*inv_denom )
        lam2_hlle = abs( (S4+S1)*(un_roe-a_roe)*inv_denom - 2d0*S4*S1*inv_denom )
        lam3_hlle = abs( (S4+S1)*(un_roe+a_roe)*inv_denom - 2d0*S4*S1*inv_denom )
        eps_ = 0.3d0*(abs(un_roe)+a_roe)
        lam1_roe = abs(un_roe)
        vv_ = un_roe - a_roe                                   ! hartenHyman(un_roe - a_roe, eps)
        lam2_roe = abs(vv_); if (abs(vv_) <= eps_) lam2_roe = (vv_*vv_ + eps_*eps_)/(2d0*eps_)
        vv_ = un_roe + a_roe                                   ! hartenHyman(un_roe + a_roe, eps)
        lam3_roe = abs(vv_); if (abs(vv_) <= eps_) lam3_roe = (vv_*vv_ + eps_*eps_)/(2d0*eps_)
        lam1_pp = bet*lam1_roe + (1.d0 - bet)*lam1_hlle
        lam2_pp = bet*lam2_roe + (1.d0 - bet)*lam2_hlle
        lam3_pp = bet*lam3_roe + (1.d0 - bet)*lam3_hlle
        ! physical fluxes (fluxes): F_r = sum(r)*un, F_E = F_r*H0(p, r, sqrt(u*u+v*v+w*w))
        rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dl(s); end do
        FrL = rho_*(ul*nx+vl*ny+wl*nz)
        Ful = FrL*ul + pl*nx; Fvl = FrL*vl + pl*ny; Fwl = FrL*wl + pl*nz
        rg_ = 0.d0; do s = 1, ns_; rg_ = rg_ + RO_RITAB(s)*dl(s)/rho_; end do
        Hs_L = 0.d0
        do s = 1, ns_
          T_ = pl/(rg_*rho_); Tc_ = max(dble(RO_TLO), min(T_, dble(RO_THI-1))); Ti_ = idint(Tc_); Td_ = Tc_ - Ti_
          Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
          hs_ = RO_HTAB(Ti1_,s) + (RO_HTAB(Ti2_,s) - RO_HTAB(Ti1_,s))*Td_
          Hs_L = Hs_L + hs_*dl(s)/rho_
        end do
        vel_ = sqrt(ul*ul+vl*vl+wl*wl)
        Fel = FrL*(Hs_L + 0.5d0*vel_*vel_)
        rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dr(s); end do
        FrR = rho_*(ur*nx+vr*ny+wr*nz)
        Fur = FrR*ur + pr*nx; Fvr = FrR*vr + pr*ny; Fwr = FrR*wr + pr*nz
        rg_ = 0.d0; do s = 1, ns_; rg_ = rg_ + RO_RITAB(s)*dr(s)/rho_; end do
        Hs_R = 0.d0
        do s = 1, ns_
          T_ = pr/(rg_*rho_); Tc_ = max(dble(RO_TLO), min(T_, dble(RO_THI-1))); Ti_ = idint(Tc_); Td_ = Tc_ - Ti_
          Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
          hs_ = RO_HTAB(Ti1_,s) + (RO_HTAB(Ti2_,s) - RO_HTAB(Ti1_,s))*Td_
          Hs_R = Hs_R + hs_*dr(s)/rho_
        end do
        vel_ = sqrt(ur*ur+vr*vr+wr*wr)
        Fer = FrR*(Hs_R + 0.5d0*vel_*vel_)
        F_r = 0.5d0*(FrL + FrR)
        F_u = 0.5d0*(Ful + Fur)
        F_v = 0.5d0*(Fvl + Fvr)
        F_w = 0.5d0*(Fwl + Fwr)
        F_E = 0.5d0*(Fel + Fer)
        drho = drtot - dltot
        dp   = pr - pl
        du   = ur - ul
        dv   = vr - vl
        dw   = wr - wl
        dun  = unr - unl
        if (abs(nx) < 0.9d0) then
          t1x = 0.d0; t1y = -nz; t1z = ny
        else
          t1x = -nz; t1y = 0.d0; t1z = nx
        endif
        normt = dsqrt(t1x*t1x+t1y*t1y+t1z*t1z)
        t1x = t1x/normt ; t1y = t1y/normt ; t1z = t1z/normt
        t2x = ny*t1z - nz*t1y
        t2y = nz*t1x - nx*t1z
        t2z = nx*t1y - ny*t1x
        alpha_m  = (dp - rho_roe*a_roe*dun)/(2d0*a_roe*a_roe)
        alpha_c  = drho - dp/(a_roe*a_roe)
        alpha_s1 = rho_roe*(du*t1x + dv*t1y + dw*t1z)
        alpha_s2 = rho_roe*(du*t2x + dv*t2y + dw*t2z)
        alpha_p  = (dp + rho_roe*a_roe*dun)/(2d0*a_roe*a_roe)
        F_r = F_r - 0.5d0 * ( alpha_c*lam1_pp + alpha_m*lam2_pp + alpha_p*lam3_pp )
        F_u = F_u - 0.5d0*( alpha_c * lam1_pp * u_roe )
        F_v = F_v - 0.5d0*( alpha_c * lam1_pp * v_roe )
        F_w = F_w - 0.5d0*( alpha_c * lam1_pp * w_roe )
        F_u = F_u - 0.5d0 * lam1_pp * ( alpha_s1*t1x + alpha_s2*t2x )
        F_v = F_v - 0.5d0 * lam1_pp * ( alpha_s1*t1y + alpha_s2*t2y )
        F_w = F_w - 0.5d0 * lam1_pp * ( alpha_s1*t1z + alpha_s2*t2z )
        F_u = F_u - 0.5d0 * ( alpha_m * lam2_pp * (u_roe - a_roe * nx) + alpha_p * lam3_pp * (u_roe + a_roe * nx) )
        F_v = F_v - 0.5d0 * ( alpha_m * lam2_pp * (v_roe - a_roe * ny) + alpha_p * lam3_pp * (v_roe + a_roe * ny) )
        F_w = F_w - 0.5d0 * ( alpha_m * lam2_pp * (w_roe - a_roe * nz) + alpha_p * lam3_pp * (w_roe + a_roe * nz) )
        F_E = F_E - 0.5d0 * lam1_pp * ( alpha_c*0.5d0*(u_roe*u_roe + v_roe*v_roe + w_roe*w_roe) )
        F_E = F_E - 0.5d0 * lam1_pp * ( alpha_s1*(u_roe*t1x+v_roe*t1y+w_roe*t1z) + alpha_s2*(u_roe*t2x+v_roe*t2y+w_roe*t2z) )
        F_E = F_E - 0.5d0 * ( alpha_m * lam2_pp * (h0_roe - a_roe * un_roe) + alpha_p * lam3_pp * (h0_roe + a_roe * un_roe) )
      else if (rs == 7) then
        ! ---------------------------------------------------------------- HLLC+Tramel (riemann_HLLCp_Tramel)
        S1 = min(0.d0, unl - al, un_roe - a_roe)
        S4 = max(0.d0, unr + ar, un_roe + a_roe)
        ! riemann_HLLE
        inv_denom = 1.d0 / (S4 - S1)
        FrL = dltot*unl
        FrR = drtot*unr
        FrE = ( S4*FrL - S1*FrR + S1*S4*(drtot-dltot) ) * inv_denom
        FuE = ( S4*(FrL*ul + pl*nx) - S1*(FrR*ur + pr*nx) &
              + S1*S4*(drtot*ur - dltot*ul) ) * inv_denom
        FvE = ( S4*(FrL*vl + pl*ny) - S1*(FrR*vr + pr*ny) &
              + S1*S4*(drtot*vr - dltot*vl) ) * inv_denom
        FwE = ( S4*(FrL*wl + pl*nz) - S1*(FrR*wr + pr*nz) &
              + S1*S4*(drtot*wr - dltot*wl) ) * inv_denom
        FeE = ( S4*(FrL*h0l) - S1*(FrR*h0r) &
              + S1*S4*( drtot*(h0r - pr/drtot) &
                        -dltot*(h0l - pl/dltot) ) ) * inv_denom
        ! riemann_HLLC
        e0l = h0l - pl / dltot
        e0r = h0r - pr / drtot
        inv_denom = 1.d0 / ( dltot*(S1-unl) - drtot*(S4-unr) )
        Sstar = ( pr - pl + dltot*unl*(S1-unl) - drtot*unr*(S4-unr) ) * inv_denom
        pstar = dltot*(unl-S1)*(unl-Sstar) + pl
        if (S1 >= 0.d0) then
          FrC = dltot * unl
          FuC = FrC*ul + pl*nx; FvC = FrC*vl + pl*ny; FwC = FrC*wl + pl*nz; FeC = FrC*h0l
        elseif (Sstar >= 0.d0) then
          FrC = dltot * unl
          FuC = FrC*ul + pl*nx; FvC = FrC*vl + pl*ny; FwC = FrC*wl + pl*nz; FeC = FrC*h0l
          inv_denom = 1.d0 / (S1 - Sstar)
          U1S = dltot*(S1-unl)*inv_denom
          U2S = ((S1-unl)*dltot*ul + (pstar-pl)*nx)*inv_denom
          U3S = ((S1-unl)*dltot*vl + (pstar-pl)*ny)*inv_denom
          U4S = ((S1-unl)*dltot*wl + (pstar-pl)*nz)*inv_denom
          U5S = ((S1-unl)*dltot*e0l - pl*unl + pstar*Sstar)*inv_denom
          FrC = FrC + S1*(U1S - dltot)
          FuC = FuC + S1*(U2S - dltot*ul)
          FvC = FvC + S1*(U3S - dltot*vl)
          FwC = FwC + S1*(U4S - dltot*wl)
          FeC = FeC + S1*(U5S - dltot*e0l)
        elseif (S4 > 0.d0) then
          FrC = drtot * unr
          FuC = FrC*ur + pr*nx; FvC = FrC*vr + pr*ny; FwC = FrC*wr + pr*nz; FeC = FrC*h0r
          inv_denom = 1.d0 / (S4 - Sstar)
          U1S = drtot*(S4-unr)*inv_denom
          U2S = ((S4-unr)*drtot*ur + (pstar-pr)*nx)*inv_denom
          U3S = ((S4-unr)*drtot*vr + (pstar-pr)*ny)*inv_denom
          U4S = ((S4-unr)*drtot*wr + (pstar-pr)*nz)*inv_denom
          U5S = ((S4-unr)*drtot*e0r - pr*unr + pstar*Sstar)*inv_denom
          FrC = FrC + S4*(U1S - drtot)
          FuC = FuC + S4*(U2S - drtot*ur)
          FvC = FvC + S4*(U3S - drtot*vr)
          FwC = FwC + S4*(U4S - drtot*wr)
          FeC = FeC + S4*(U5S - drtot*e0r)
        else
          FrC = drtot * unr
          FuC = FrC*ur + pr*nx; FvC = FrC*vr + pr*ny; FwC = FrC*wr + pr*nz; FeC = FrC*h0r
        endif
        if (bet == 0d0) then
          F_r = FrE; F_u = FuE; F_v = FvE; F_w = FwE; F_E = FeE
        elseif (bet >= 1d0) then
          F_r = FrC; F_u = FuC; F_v = FvC; F_w = FwC; F_E = FeC
        else
          F_r = bet*FrC + (1d0-bet)*FrE
          F_u = bet*FuC + (1d0-bet)*FuE
          F_v = bet*FvC + (1d0-bet)*FvE
          F_w = bet*FwC + (1d0-bet)*FwE
          F_E = bet*FeC + (1d0-bet)*FeE
        endif
      else if (rs == 10) then
        ! ---------------------------------------------------------------- HLLC+Chen (riemann_HLLCp_Chen)
        e0l = h0l - pl/dltot
        e0r = h0r - pr/drtot
        S1 = min(0.d0, unl - al, un_roe - a_roe)
        S4 = max(0.d0, unr + ar, un_roe + a_roe)
        inv_denom = 1.d0 / ( dltot*(S1-unl) - drtot*(S4-unr) )
        Sstar = ( pr - pl + dltot*unl*(S1-unl) - drtot*unr*(S4-unr) ) * inv_denom
        pstar = dltot*(unl-S1)*(unl-Sstar) + pl
        ! Chen all-speed corrections, active only in the subsonic star region (switch = bet = the sensor h)
        ch_phiL = dltot*(S1 - unl)
        ch_phiR = drtot*(S4 - unr)
        ch_pref = ch_phiL*ch_phiR / (ch_phiR - ch_phiL)
        ch_Mloc = min(1.d0, max( sqrt(ul*ul+vl*vl+wl*wl)/al, sqrt(ur*ur+vr*vr+wr*wr)/ar ))
        ch_Mloc = max(RO_MCO, ch_Mloc)
        ch_fM   = ch_Mloc*sqrt(4.d0 + (1.d0 - ch_Mloc*ch_Mloc)**2)/(1.d0 + ch_Mloc*ch_Mloc)
        ch_g    = 1.d0 - bet**ch_Mloc
        ch_fstar = ch_fM + ch_g*(1.d0 - ch_fM)
        ch_dU   = unr - unl
        ch_dutx = (ur-ul) - ch_dU*nx
        ch_duty = (vr-vl) - ch_dU*ny
        ch_dutz = (wr-wl) - ch_dU*nz
        if (S1 >= 0.d0) then
          F_r = dltot*unl
          F_u = F_r*ul + pl*nx
          F_v = F_r*vl + pl*ny
          F_w = F_r*wl + pl*nz
          F_E = F_r*h0l
        elseif (Sstar >= 0.d0) then
          F_r = dltot*unl
          F_u = F_r*ul + pl*nx
          F_v = F_r*vl + pl*ny
          F_w = F_r*wl + pl*nz
          F_E = F_r*h0l
          inv_denom = 1.d0 / (S1 - Sstar)
          U1S = dltot*(S1-unl)*inv_denom
          U2S = ((S1-unl)*dltot*ul + (pstar-pl)*nx)*inv_denom
          U3S = ((S1-unl)*dltot*vl + (pstar-pl)*ny)*inv_denom
          U4S = ((S1-unl)*dltot*wl + (pstar-pl)*nz)*inv_denom
          U5S = ((S1-unl)*dltot*e0l - pl*unl + pstar*Sstar)*inv_denom
          F_r = F_r + S1*(U1S - dltot)
          F_u = F_u + S1*(U2S - dltot*ul)
          F_v = F_v + S1*(U3S - dltot*vl)
          F_w = F_w + S1*(U4S - dltot*wl)
          F_E = F_E + S1*(U5S - dltot*e0l)
          ch_SK = S1
          ch_lm = ch_pref*(ch_fstar - 1.d0)*ch_dU
          ch_sh = ch_pref*( ch_SK/(ch_SK - Sstar) )*ch_g
          F_u = F_u + ch_lm*nx + ch_sh*ch_dutx
          F_v = F_v + ch_lm*ny + ch_sh*ch_duty
          F_w = F_w + ch_lm*nz + ch_sh*ch_dutz
          F_E = F_E + ch_lm*Sstar
        elseif (S4 > 0.d0) then
          F_r = drtot*unr
          F_u = F_r*ur + pr*nx
          F_v = F_r*vr + pr*ny
          F_w = F_r*wr + pr*nz
          F_E = F_r*h0r
          inv_denom = 1.d0 / (S4 - Sstar)
          U1S = drtot*(S4-unr)*inv_denom
          U2S = ((S4-unr)*drtot*ur + (pstar-pr)*nx)*inv_denom
          U3S = ((S4-unr)*drtot*vr + (pstar-pr)*ny)*inv_denom
          U4S = ((S4-unr)*drtot*wr + (pstar-pr)*nz)*inv_denom
          U5S = ((S4-unr)*drtot*e0r - pr*unr + pstar*Sstar)*inv_denom
          F_r = F_r + S4*(U1S - drtot)
          F_u = F_u + S4*(U2S - drtot*ur)
          F_v = F_v + S4*(U3S - drtot*vr)
          F_w = F_w + S4*(U4S - drtot*wr)
          F_E = F_E + S4*(U5S - drtot*e0r)
          ch_SK = S4
          ch_lm = ch_pref*(ch_fstar - 1.d0)*ch_dU
          ch_sh = ch_pref*( ch_SK/(ch_SK - Sstar) )*ch_g
          F_u = F_u + ch_lm*nx + ch_sh*ch_dutx
          F_v = F_v + ch_lm*ny + ch_sh*ch_duty
          F_w = F_w + ch_lm*nz + ch_sh*ch_dutz
          F_E = F_E + ch_lm*Sstar
        else
          F_r = drtot*unr
          F_u = F_r*ur + pr*nx
          F_v = F_r*vr + pr*ny
          F_w = F_r*wr + pr*nz
          F_E = F_r*h0r
        endif
      else
        ! ---------------------------------------------------------------- LMRoe (Riemann_LMRoe)
        lm_nmag = sqrt(nx*nx + ny*ny + nz*nz)
        if (lm_nmag <= 0.d0) then
          F_r = 0.d0; F_u = 0.d0; F_v = 0.d0; F_w = 0.d0; F_E = 0.d0
        else
        lm_n1 = nx/lm_nmag; lm_n2 = ny/lm_nmag; lm_n3 = nz/lm_nmag
        if (abs(lm_n1) < 0.9d0) then
          lm_tmpx = 1.d0; lm_tmpy = 0.d0; lm_tmpz = 0.d0
        else
          lm_tmpx = 0.d0; lm_tmpy = 1.d0; lm_tmpz = 0.d0
        end if
        t1x = lm_tmpx - (lm_tmpx*lm_n1 + lm_tmpy*lm_n2 + lm_tmpz*lm_n3)*lm_n1
        t1y = lm_tmpy - (lm_tmpx*lm_n1 + lm_tmpy*lm_n2 + lm_tmpz*lm_n3)*lm_n2
        t1z = lm_tmpz - (lm_tmpx*lm_n1 + lm_tmpy*lm_n2 + lm_tmpz*lm_n3)*lm_n3
        lm_tmag = sqrt(t1x*t1x + t1y*t1y + t1z*t1z)
        if (lm_tmag <= 0.d0) then
          t1x = -lm_n2; t1y =  lm_n1; t1z = 0.d0
          lm_tmag = sqrt(t1x*t1x + t1y*t1y + t1z*t1z)
        end if
        t1x = t1x/lm_tmag; t1y = t1y/lm_tmag; t1z = t1z/lm_tmag
        t2x = lm_n2*t1z - lm_n3*t1y
        t2y = lm_n3*t1x - lm_n1*t1z
        t2z = lm_n1*t1y - lm_n2*t1x
        lm_unL  = ul*lm_n1 + vl*lm_n2 + wl*lm_n3
        lm_unR  = ur*lm_n1 + vr*lm_n2 + wr*lm_n3
        lm_ut1L = ul*t1x + vl*t1y + wl*t1z
        lm_ut1R = ur*t1x + vr*t1y + wr*t1z
        lm_ut2L = ul*t2x + vl*t2y + wl*t2z
        lm_ut2R = ur*t2x + vr*t2y + wr*t2z
        ! hL = H(pl,dl), hR = H(pr,dr); H0 = h + |u|^2/2 enters the physical fluxes (computed before roe_averages)
        rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dl(s); end do
        rg_ = 0.d0; do s = 1, ns_; rg_ = rg_ + RO_RITAB(s)*dl(s)/rho_; end do
        lm_hL = 0.d0
        do s = 1, ns_
          T_ = pl/(rg_*rho_); Tc_ = max(dble(RO_TLO), min(T_, dble(RO_THI-1))); Ti_ = idint(Tc_); Td_ = Tc_ - Ti_
          Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
          hs_ = RO_HTAB(Ti1_,s) + (RO_HTAB(Ti2_,s) - RO_HTAB(Ti1_,s))*Td_
          lm_hL = lm_hL + hs_*dl(s)/rho_
        end do
        rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dr(s); end do
        rg_ = 0.d0; do s = 1, ns_; rg_ = rg_ + RO_RITAB(s)*dr(s)/rho_; end do
        lm_hR = 0.d0
        do s = 1, ns_
          T_ = pr/(rg_*rho_); Tc_ = max(dble(RO_TLO), min(T_, dble(RO_THI-1))); Ti_ = idint(Tc_); Td_ = Tc_ - Ti_
          Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
          hs_ = RO_HTAB(Ti1_,s) + (RO_HTAB(Ti2_,s) - RO_HTAB(Ti1_,s))*Td_
          lm_hR = lm_hR + hs_*dr(s)/rho_
        end do
        lm_H0L = lm_hL + 0.5d0*(ul*ul + vl*vl + wl*wl)
        lm_H0R = lm_hR + 0.5d0*(ur*ur + vr*vr + wr*wr)
        lm_FL1 = dltot*lm_unL
        lm_FL2 = dltot*lm_unL*ul + pl*lm_n1
        lm_FL3 = dltot*lm_unL*vl + pl*lm_n2
        lm_FL4 = dltot*lm_unL*wl + pl*lm_n3
        lm_FL5 = dltot*lm_unL*lm_H0L
        lm_FR1 = drtot*lm_unR
        lm_FR2 = drtot*lm_unR*ur + pr*lm_n1
        lm_FR3 = drtot*lm_unR*vr + pr*lm_n2
        lm_FR4 = drtot*lm_unR*wr + pr*lm_n3
        lm_FR5 = drtot*lm_unR*lm_H0R
        lm_unRoe = u_roe*lm_n1 + v_roe*lm_n2 + w_roe*lm_n3
        lm_u2 = u_roe*u_roe + v_roe*v_roe + w_roe*w_roe
        lm_utmag = sqrt( max(0.d0, lm_u2 - lm_unRoe*lm_unRoe ) )
        if (a_roe > 0.d0) then
          lm_Ma = (abs(lm_unRoe) + lm_utmag)/a_roe
        else
          lm_Ma = 1.d0
        end if
        lm_phi = max(RO_MCO, min(1.d0, lm_Ma))
        drho = drtot - dltot
        dp   = pr - pl
        dun  = lm_unR - lm_unL
        lm_dut1 = lm_ut1R - lm_ut1L
        lm_dut2 = lm_ut2R - lm_ut2L
        lm_dunmod = lm_phi * dun
        lm_abs1 = abs(lm_unRoe - a_roe); lm_abs2 = abs(lm_unRoe); lm_abs3 = abs(lm_unRoe); lm_abs4 = abs(lm_unRoe)
        lm_abs5 = abs(lm_unRoe + a_roe)
        lm_a2    = a_roe*a_roe
        lm_inva  = 1.d0/a_roe
        lm_inva2 = 1.d0/lm_a2
        lm_al4 = drho - dp*lm_inva2
        lm_al2 = rho_roe * lm_dut1
        lm_al3 = rho_roe * lm_dut2
        lm_al1 = 0.5d0*( dp*lm_inva2 - rho_roe*lm_dunmod*lm_inva )
        lm_al5 = 0.5d0*( dp*lm_inva2 + rho_roe*lm_dunmod*lm_inva )
        lm_r1_1 = 1.d0
        lm_r1_2 = u_roe - a_roe*lm_n1
        lm_r1_3 = v_roe - a_roe*lm_n2
        lm_r1_4 = w_roe - a_roe*lm_n3
        lm_r1_5 = h0_roe - a_roe*lm_unRoe
        lm_r5_1 = 1.d0
        lm_r5_2 = u_roe + a_roe*lm_n1
        lm_r5_3 = v_roe + a_roe*lm_n2
        lm_r5_4 = w_roe + a_roe*lm_n3
        lm_r5_5 = h0_roe + a_roe*lm_unRoe
        lm_r4_1 = 1.d0
        lm_r4_2 = u_roe
        lm_r4_3 = v_roe
        lm_r4_4 = w_roe
        lm_r4_5 = 0.5d0*lm_u2
        lm_r2_1 = 0.d0
        lm_r2_2 = t1x
        lm_r2_3 = t1y
        lm_r2_4 = t1z
        lm_r2_5 = u_roe*t1x + v_roe*t1y + w_roe*t1z
        lm_r3_1 = 0.d0
        lm_r3_2 = t2x
        lm_r3_3 = t2y
        lm_r3_4 = t2z
        lm_r3_5 = u_roe*t2x + v_roe*t2y + w_roe*t2z
        lm_d1 = lm_abs1*lm_al1*lm_r1_1 + lm_abs2*lm_al2*lm_r2_1 + lm_abs3*lm_al3*lm_r3_1 &
              + lm_abs4*lm_al4*lm_r4_1 + lm_abs5*lm_al5*lm_r5_1
        lm_d2 = lm_abs1*lm_al1*lm_r1_2 + lm_abs2*lm_al2*lm_r2_2 + lm_abs3*lm_al3*lm_r3_2 &
              + lm_abs4*lm_al4*lm_r4_2 + lm_abs5*lm_al5*lm_r5_2
        lm_d3 = lm_abs1*lm_al1*lm_r1_3 + lm_abs2*lm_al2*lm_r2_3 + lm_abs3*lm_al3*lm_r3_3 &
              + lm_abs4*lm_al4*lm_r4_3 + lm_abs5*lm_al5*lm_r5_3
        lm_d4 = lm_abs1*lm_al1*lm_r1_4 + lm_abs2*lm_al2*lm_r2_4 + lm_abs3*lm_al3*lm_r3_4 &
              + lm_abs4*lm_al4*lm_r4_4 + lm_abs5*lm_al5*lm_r5_4
        lm_d5 = lm_abs1*lm_al1*lm_r1_5 + lm_abs2*lm_al2*lm_r2_5 + lm_abs3*lm_al3*lm_r3_5 &
              + lm_abs4*lm_al4*lm_r4_5 + lm_abs5*lm_al5*lm_r5_5
        F_r = 0.5d0*(lm_FL1 + lm_FR1) - 0.5d0*lm_d1
        F_u = 0.5d0*(lm_FL2 + lm_FR2) - 0.5d0*lm_d2
        F_v = 0.5d0*(lm_FL3 + lm_FR3) - 0.5d0*lm_d3
        F_w = 0.5d0*(lm_FL4 + lm_FR4) - 0.5d0*lm_d4
        F_E = 0.5d0*(lm_FL5 + lm_FR5) - 0.5d0*lm_d5
        end if
      end if
    case (14)
      ! ---------------------------------------------------------------- HLLC-PC (Riemann_HLLCprec)
      ! HLLC Batten with the signal velocities of the preconditioned system: Vr = min(Ur, a) of the two cells, the
      ! averages with the plain sqrt(rho) weights (sound speed included), alpha = (1 - Vr^2/a^2)/2
      unl = ul*nx + vl*ny + wl*nz
      unr = ur*nx + vr*ny + wr*nz
      rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dl(s); end do
      rg_ = 0.d0; do s = 1, ns_; rg_ = rg_ + RO_RITAB(s)*dl(s)/rho_; end do
      pc_hL = 0.d0
      do s = 1, ns_
        T_ = pl/(rg_*rho_); Tc_ = max(dble(RO_TLO), min(T_, dble(RO_THI-1))); Ti_ = idint(Tc_); Td_ = Tc_ - Ti_
        Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
        hs_ = RO_HTAB(Ti1_,s) + (RO_HTAB(Ti2_,s) - RO_HTAB(Ti1_,s))*Td_
        pc_hL = pc_hL + hs_*dl(s)/rho_
      end do
      rho_ = 0.d0; do s = 1, ns_; rho_ = rho_ + dr(s); end do
      rg_ = 0.d0; do s = 1, ns_; rg_ = rg_ + RO_RITAB(s)*dr(s)/rho_; end do
      pc_hR = 0.d0
      do s = 1, ns_
        T_ = pr/(rg_*rho_); Tc_ = max(dble(RO_TLO), min(T_, dble(RO_THI-1))); Ti_ = idint(Tc_); Td_ = Tc_ - Ti_
        Ti1_ = max(RO_TLO, min(Ti_, RO_THI)); Ti2_ = max(RO_TLO, min(Ti_+1, RO_THI))
        hs_ = RO_HTAB(Ti1_,s) + (RO_HTAB(Ti2_,s) - RO_HTAB(Ti1_,s))*Td_
        pc_hR = pc_hR + hs_*dr(s)/rho_
      end do
      pc_eL = pc_hL - pl/dltot
      pc_eR = pc_hR - pr/drtot
      h0l = pc_hL + 0.5d0*(ul**2 + vl**2 + wl**2)
      h0r = pc_hR + 0.5d0*(ur**2 + vr**2 + wr**2)
      pc_EtL = pc_eL + 0.5d0*(ul**2 + vl**2 + wl**2)
      pc_EtR = pc_eR + 0.5d0*(ur**2 + vr**2 + wr**2)
      FrL = dltot * unl
      Ful = FrL * ul + pl * nx; Fvl = FrL * vl + pl * ny; Fwl = FrL * wl + pl * nz; Fel = FrL * h0l
      FrR = drtot * unr
      Fur = FrR * ur + pr * nx; Fvr = FrR * vr + pr * ny; Fwr = FrR * wr + pr * nz; Fer = FrR * h0r
      srL = sqrt(dltot)
      srR = sqrt(drtot)
      pc_den = srL + srR
      u_roe = (srL*ul + srR*ur) / pc_den
      v_roe = (srL*vl + srR*vr) / pc_den
      w_roe = (srL*wl + srR*wr) / pc_den
      a_roe = (srL*al + srR*ar) / pc_den
      un_roe = u_roe*nx + v_roe*ny + w_roe*nz
      pc_VrL   = min(urL_, al)
      pc_VrR   = min(urR_, ar)
      pc_VrRoe = (srL*pc_VrL + srR*pc_VrR) / pc_den
      pc_alL   = (1.0d0 - (pc_VrL  **2) / (al  **2)) / 2.0d0
      pc_alR   = (1.0d0 - (pc_VrR  **2) / (ar  **2)) / 2.0d0
      pc_alRoe = (1.0d0 - (pc_VrRoe**2) / (a_roe**2)) / 2.0d0
      pc_vLs   = unl    * (1.0d0 - pc_alL)
      pc_vRs   = unr    * (1.0d0 - pc_alR)
      pc_vRoes = un_roe * (1.0d0 - pc_alRoe)
      pc_cLs   = sqrt(max(pc_alL  **2 * unl   **2 + pc_VrL  **2, 0.0d0))
      pc_cRs   = sqrt(max(pc_alR  **2 * unr   **2 + pc_VrR  **2, 0.0d0))
      pc_cRoes = sqrt(max(pc_alRoe**2 * un_roe**2 + pc_VrRoe**2, 0.0d0))
      S1 = min(0.0d0, pc_vLs - pc_cLs, pc_vRoes - pc_cRoes)
      S4 = max(0.0d0, pc_vRs + pc_cRs, pc_vRoes + pc_cRoes)
      Sstar = ( pr - pl + dltot*unl*(S1 - unl) - drtot*unr*(S4 - unr) ) / ( dltot*(S1 - unl) - drtot*(S4 - unr) )
      pstar = dltot*(unl - S1)*(unl - Sstar) + pl
      if (S1 >= 0.0d0) then
        F_r = FrL; F_u = Ful; F_v = Fvl; F_w = Fwl; F_E = Fel
      else if (S4 <= 0.0d0) then
        F_r = FrR; F_u = Fur; F_v = Fvr; F_w = Fwr; F_E = Fer
      else if (Sstar >= 0.0d0) then
        U1S = dltot*(S1 - unl)/(S1 - Sstar)
        U2S = ((S1 - unl)*dltot*ul + (pstar - pl)*nx)/(S1 - Sstar)
        U3S = ((S1 - unl)*dltot*vl + (pstar - pl)*ny)/(S1 - Sstar)
        U4S = ((S1 - unl)*dltot*wl + (pstar - pl)*nz)/(S1 - Sstar)
        U5S = ((S1 - unl)*dltot*pc_EtL - pl*unl + pstar*Sstar)/(S1 - Sstar)
        F_r = FrL + S1*(U1S - dltot)
        F_u = Ful + S1*(U2S - dltot*ul)
        F_v = Fvl + S1*(U3S - dltot*vl)
        F_w = Fwl + S1*(U4S - dltot*wl)
        F_E = Fel + S1*(U5S - dltot*pc_EtL)
      else
        U1S = drtot*(S4 - unr)/(S4 - Sstar)
        U2S = ((S4 - unr)*drtot*ur + (pstar - pr)*nx)/(S4 - Sstar)
        U3S = ((S4 - unr)*drtot*vr + (pstar - pr)*ny)/(S4 - Sstar)
        U4S = ((S4 - unr)*drtot*wr + (pstar - pr)*nz)/(S4 - Sstar)
        U5S = ((S4 - unr)*drtot*pc_EtR - pr*unr + pstar*Sstar)/(S4 - Sstar)
        F_r = FrR + S4*(U1S - drtot)
        F_u = Fur + S4*(U2S - drtot*ur)
        F_v = Fvr + S4*(U3S - drtot*vr)
        F_w = Fwr + S4*(U4S - drtot*wr)
        F_E = Fer + S4*(U5S - drtot*pc_EtR)
      end if
    case default
      F_r = 0.d0; F_u = 0.d0; F_v = 0.d0; F_w = 0.d0; F_E = 0.d0   ! never: Assign_RS refuses the other ids
    end select
