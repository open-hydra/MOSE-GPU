    ! ===== inlined convective physics (conv_physics.inc) =====
    ! inputs:  pm1,p0,pp1,pp2 (4-cell stencil), dlc(-1:2), nrm(3), area, bet
    ! output:  res0(1:nprim) = +Flux  (== old rloc(:,0))
    ! --- Reconstruction (the limiter of the case, limv) ---
    l0 = 0.5d0*(dlc(-1)+dlc(0)); l1 = 0.5d0*(dlc(0)+dlc(1)); l2 = 0.5d0*(dlc(1)+dlc(2))
    lm = 0.5d0*dlc(0); lp = 0.5d0*dlc(1)
    if (SD) then; beta_order = bet; else; beta_order = 1d0; endif
    idl0 = 1d0/l0; idl1 = 1d0/l1; idl2 = 1d0/l2
    do ss = 1, nprim
      slope0 = (p0(ss)-pm1(ss))*idl0
      slope1 = (pp1(ss)-p0(ss))*idl1
      slope2 = (pp2(ss)-pp1(ss))*idl2
      slopel = rlimiter_d(limv, slope1, slope0)   ! the limiter of the case (lim_dev): minmod verbatim for limv = 1
      sloper = rlimiter_d(limv, slope2, slope1)
      Prim1(ss) = p0(ss)  + beta_order*slopel*lm
      Prim4(ss) = pp1(ss) - beta_order*sloper*lp
    end do
    do ss = 1, nsc
      Prim1(ss) = max(Prim1(ss), 1d-20)
      Prim4(ss) = max(Prim4(ss), 1d-20)
    end do
    ! --- co_rotot_Rtot (left=1, right=4) ---
    rotot1 = 0d0; Rtot1 = 0d0
    do ss = 1, nsc; rotot1 = rotot1 + Prim1(ss); end do
    do ss = 1, nsc; Rtot1 = Rtot1 + Ri_tab_d(ss)*Prim1(ss)/rotot1; end do
    rotot4 = 0d0; Rtot4 = 0d0
    do ss = 1, nsc; rotot4 = rotot4 + Prim4(ss); end do
    do ss = 1, nsc; Rtot4 = Rtot4 + Ri_tab_d(ss)*Prim4(ss)/rotot4; end do
    ! --- f_ss (sound speed a1) ---
    T1d = Prim1(np)/(Rtot1*rotot1)
    Ti1 = idint(T1d); if (Ti1<0) Ti1=0; if (Ti1>NTMAX-1) Ti1=NTMAX-1
    Tdiffd = T1d - Ti1; cpd = 0d0
    do ss = 1, nsc
      cpi = cp_tab_d(Ti1,ss) + (cp_tab_d(Ti1+1,ss)-cp_tab_d(Ti1,ss))*Tdiffd
      cpd = cpd + Prim1(ss)/rotot1*cpi
    end do
    gamd = cpd/(cpd-Rtot1); a1 = dsqrt(gamd*Rtot1*T1d)
    ! --- f_ss (sound speed a4) ---
    T4d = Prim4(np)/(Rtot4*rotot4)
    Ti4 = idint(T4d); if (Ti4<0) Ti4=0; if (Ti4>NTMAX-1) Ti4=NTMAX-1
    Tdiffd = T4d - Ti4; cpd = 0d0
    do ss = 1, nsc
      cpi = cp_tab_d(Ti4,ss) + (cp_tab_d(Ti4+1,ss)-cp_tab_d(Ti4,ss))*Tdiffd
      cpd = cpd + Prim4(ss)/rotot4*cpi
    end do
    gamd = cpd/(cpd-Rtot4); a4 = dsqrt(gamd*Rtot4*T4d)
    ! --- riemann_HLLC (roe_averages inline) ---
    dltot = rotot1; drtot = rotot4
    ul = Prim1(nu); vl = Prim1(nv); wl = Prim1(nw); pl = Prim1(np); al = a1
    ur = Prim4(nu); vr = Prim4(nv); wr = Prim4(nw); pr = Prim4(np); ar = a4
    nx = nrm(1); ny = nrm(2); nz = nrm(3)
    unl = ul*nx + vl*ny + wl*nz
    unr = ur*nx + vr*ny + wr*nz
    Rl = Rtot1; Rr = Rtot4
    inv_dltot = 1d0/dltot; inv_drtot = 1d0/drtot
    Tl = pl*inv_dltot/Rl; Tr = pr*inv_drtot/Rr
    Til = int(Tl); if (Til<0) Til=0; if (Til>NTMAX-1) Til=NTMAX-1
    Tir = int(Tr); if (Tir<0) Tir=0; if (Tir>NTMAX-1) Tir=NTMAX-1
    dTl = Tl - Til; dTr = Tr - Tir
    srL = sqrt(dltot); srR = sqrt(drtot); inv_sr = 1d0/(srL+srR)
    rho_roe = srL*srR
    u_roe = (srR*ur+srL*ul)*inv_sr; v_roe = (srR*vr+srL*vl)*inv_sr; w_roe = (srR*wr+srL*wl)*inv_sr
    R_roe = (srR*Rr+srL*Rl)*inv_sr; T_roe = (srR*Tr+srL*Tl)*inv_sr
    h0l = 0d0; h0r = 0d0; cv_roe = 0d0; sum_ei = 0d0
    do ss = 1, nsc
      invW = Runiv/wm_tab_d(ss)
      hl = h_tab_d(Til,ss) + (h_tab_d(Til+1,ss)-h_tab_d(Til,ss))*dTl
      hr = h_tab_d(Tir,ss) + (h_tab_d(Tir+1,ss)-h_tab_d(Tir,ss))*dTr
      el = hl - invW*Tl; er = hr - invW*Tr
      h0l = h0l + hl*Prim1(ss)*inv_dltot
      h0r = h0r + hr*Prim4(ss)*inv_drtot
      d_roe = (srR*Prim4(ss)*inv_drtot + srL*Prim1(ss)*inv_dltot)*inv_sr
      cv_roe = cv_roe + d_roe*( 0.5d0*( cp_tab_d(Til,ss)+(cp_tab_d(Til+1,ss)-cp_tab_d(Til,ss))*dTl + &
               cp_tab_d(Tir,ss)+(cp_tab_d(Tir+1,ss)-cp_tab_d(Tir,ss))*dTr ) - invW )
      sum_ei = sum_ei + d_roe*(srR*er+srL*el)*inv_sr
    end do
    h0l = h0l + 0.5d0*(ul*ul+vl*vl+wl*wl)
    h0r = h0r + 0.5d0*(ur*ur+vr*vr+wr*wr)
    h0_roe = (srR*h0r+srL*h0l)*inv_sr
    gam_roe = 1d0 + R_roe/cv_roe
    vel2 = u_roe*u_roe + v_roe*v_roe + w_roe*w_roe
    a_roe = sqrt( (gam_roe-1d0)*( h0_roe - 0.5d0*vel2 + cv_roe*T_roe - sum_ei ) )
    un_roe = u_roe*nx + v_roe*ny + w_roe*nz
    e0l = h0l - pl/dltot; e0r = h0r - pr/drtot
    S1 = min(0d0, unl-al, un_roe-a_roe)
    S4 = max(0d0, unr+ar, un_roe+a_roe)
    if (rsv == 1) then   ! HLLC Batten (verbatim)
    inv_denom = 1d0/( dltot*(S1-unl) - drtot*(S4-unr) )
      Sstar = ( pr-pl + dltot*unl*(S1-unl) - drtot*unr*(S4-unr) )*inv_denom
      pstar = dltot*(unl-S1)*(unl-Sstar) + pl
      if (S1 >= 0d0) then
        F_r = dltot*unl; F_u = F_r*ul+pl*nx; F_v = F_r*vl+pl*ny; F_w = F_r*wl+pl*nz; F_E = F_r*h0l
      elseif (Sstar >= 0d0) then
        F_r = dltot*unl; F_u = F_r*ul+pl*nx; F_v = F_r*vl+pl*ny; F_w = F_r*wl+pl*nz; F_E = F_r*h0l
        inv_denom = 1d0/(S1-Sstar)
        U1S = dltot*(S1-unl)*inv_denom
        U2S = ((S1-unl)*dltot*ul + (pstar-pl)*nx)*inv_denom
        U3S = ((S1-unl)*dltot*vl + (pstar-pl)*ny)*inv_denom
        U4S = ((S1-unl)*dltot*wl + (pstar-pl)*nz)*inv_denom
        U5S = ((S1-unl)*dltot*e0l - pl*unl + pstar*Sstar)*inv_denom
        F_r = F_r + S1*(U1S-dltot); F_u = F_u + S1*(U2S-dltot*ul); F_v = F_v + S1*(U3S-dltot*vl)
        F_w = F_w + S1*(U4S-dltot*wl); F_E = F_E + S1*(U5S-dltot*e0l)
      elseif (S4 > 0d0) then
        F_r = drtot*unr; F_u = F_r*ur+pr*nx; F_v = F_r*vr+pr*ny; F_w = F_r*wr+pr*nz; F_E = F_r*h0r
        inv_denom = 1d0/(S4-Sstar)
        U1S = drtot*(S4-unr)*inv_denom
        U2S = ((S4-unr)*drtot*ur + (pstar-pr)*nx)*inv_denom
        U3S = ((S4-unr)*drtot*vr + (pstar-pr)*ny)*inv_denom
        U4S = ((S4-unr)*drtot*wr + (pstar-pr)*nz)*inv_denom
        U5S = ((S4-unr)*drtot*e0r - pr*unr + pstar*Sstar)*inv_denom
        F_r = F_r + S4*(U1S-drtot); F_u = F_u + S4*(U2S-drtot*ur); F_v = F_v + S4*(U3S-drtot*vr)
        F_w = F_w + S4*(U4S-drtot*wr); F_E = F_E + S4*(U5S-drtot*e0r)
      else
        F_r = drtot*unr; F_u = F_r*ur+pr*nx; F_v = F_r*vr+pr*ny; F_w = F_r*wr+pr*nz; F_E = F_r*h0r
      endif
    else if (rsv == 2) then
      ! HLLE (Lib_Riemann_HLL riemann_HLLE): S1/S4 of the HLLC above, two-wave flux
      inv_denom = 1.d0 / (S4 - S1)
      frl_ = dltot*unl; frr_ = drtot*unr
      F_r = ( S4*frl_ - S1*frr_ + S1*S4*(drtot-dltot) ) * inv_denom
      F_u = ( S4*(frl_*ul + pl*nx) - S1*(frr_*ur + pr*nx) + S1*S4*(drtot*ur - dltot*ul) ) * inv_denom
      F_v = ( S4*(frl_*vl + pl*ny) - S1*(frr_*vr + pr*ny) + S1*S4*(drtot*vr - dltot*vl) ) * inv_denom
      F_w = ( S4*(frl_*wl + pl*nz) - S1*(frr_*wr + pr*nz) + S1*S4*(drtot*wr - dltot*wl) ) * inv_denom
      F_E = ( S4*(frl_*h0l) - S1*(frr_*h0r) + S1*S4*( drtot*(h0r - pr/drtot) - dltot*(h0l - pl/dltot) ) ) * inv_denom
    else if (rsv == 3) then
      ! SLAU (Lib_Riemann_SLAU riemann_SLAU_base + riemann_SLAU): h0l/h0r are the total enthalpies H0 of the host
      af_ = 0.5d0*(al + ar)
      ml_ = unl / af_; mr_ = unr / af_
      mt_ = min(1.d0, sqrt(0.5d0*((ul*ul + vl*vl + wl*wl) + (ur*ur + vr*vr + wr*wr))) / af_)
      chi_ = (1.d0 - mt_)**2
      vnm_ = (dltot*abs(unl) + drtot*abs(unr)) / (dltot + drtot)
      frho_ = -max(min(ml_,0.d0), -1.d0) * min(max(mr_,0.d0), 1.d0)
      vnl_ = (1.d0 - frho_)*vnm_ + frho_*abs(unl)
      vnr_ = (1.d0 - frho_)*vnm_ + frho_*abs(unr)
      F_r = 0.5d0 * ( dltot*(unl + vnl_) + drtot*(unr - vnr_) - (chi_/af_)*(pr - pl) )
      if (abs(ml_) < 1.d0) then
        bl_ = 0.25d0*(2.d0 - ml_)*(ml_ + 1.d0)**2
      else if (ml_ >= 0.d0) then
        bl_ = 1.d0
      else
        bl_ = 0.d0
      endif
      if (abs(mr_) < 1.d0) then
        br_ = 0.25d0*(2.d0 + mr_)*(mr_ - 1.d0)**2
      else if (mr_ >= 0.d0) then
        br_ = 0.d0
      else
        br_ = 1.d0
      endif
      ph_ = 0.5d0*(pl + pr) + 0.5d0*(bl_ - br_)*(pl - pr) + (1.d0 - chi_)*(bl_ + br_ - 1.d0)*0.5d0*(pl + pr)
      if (F_r >= 0.d0) then
        F_u = F_r*ul + nx*ph_; F_v = F_r*vl + ny*ph_; F_w = F_r*wl + nz*ph_; F_E = F_r*h0l
      else
        F_u = F_r*ur + nx*ph_; F_v = F_r*vr + ny*ph_; F_w = F_r*wr + nz*ph_; F_E = F_r*h0r
      endif
    else
      ! AUSM+ (4), AUSM+M (5), HLLE++ (6), HLLC+Tramel (7), LLF (8), SLAU2 (9), HLLC+Chen (10), exact (11), LMRoe (12),
      ! MiczekRoe (13), HLLC-PC (14) of MOSE_open (riemann_open_body.inc); bet = the sensor of the face; HLLC-PC takes the
      ! Ur of the two cells of the face (MOSE_open Mod_Fluxes: blk%Ur(i), blk%Ur(i+1))
      urL_ = 0d0; urR_ = 0d0
      if (rsv == 14) then
        urL_ = urf(c0); urR_ = urf(cp1)
      end if
      call riemann_open_d ( rsv, bet, 0d0, Runiv, nsc, &   ! Mco: read in place (RO_MCO = mco_d(1)), never passed
                            Prim1(1:nsc), ul, vl, wl, pl, al, dltot, Prim4(1:nsc), ur, vr, wr, pr, ar, drtot, nx, ny, nz, &
                            urL_, urR_, F_r, F_u, F_v, F_w, F_E )
    end if
    ! --- Flux assembly -> res0 ---
    su = sign(0.5d0, F_r); sp = 0.5d0+su; sm = su-0.5d0
    irotot1 = 1d0/rotot1; irotot4 = 1d0/rotot4
    do ss = 1, nsc
      res0(ss) = F_r*area*( sp*(Prim1(ss)*irotot1) - sm*(Prim4(ss)*irotot4) )
    end do
    res0(nu) = F_u*area; res0(nv) = F_v*area; res0(nw) = F_w*area; res0(np) = F_E*area
    if (nprim > np) then
      do ss = np+1, nprim
        res0(ss) = F_r*area*( sp*Prim1(ss)*irotot1 - sm*Prim4(ss)*irotot4 )
      end do
    end if
