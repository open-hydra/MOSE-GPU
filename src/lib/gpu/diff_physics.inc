    ! ===== inlined DIFFUSIVE physics (diff_physics.inc) =====
    ! inputs:  p1..p10 (10-cell stencil), rho_n10(10), T_n10(10), M1l/M2l(3,3),
    !          nrm(3), area, wd1, wd2, rho_f,T_f,mil_f,kl_f,cp_f, ia/ib/ic (axes)
    ! output:  res0(1:nprim) = +Flux  (== old Res2z)
    ! --- gradients (Diffusive_Flux_cached): normal pair on axis ia ---
    irho1 = 1d0/rho_n10(1); irho2 = 1d0/rho_n10(2)
    Grad(1:nsc, ia) = p2(1:nsc)*irho2 - p1(1:nsc)*irho1
    Grad(np, ia)    = T_n10(2) - T_n10(1)
    Grad(nu:nw, ia) = p2(nu:nw) - p1(nu:nw)
    waldis = 0d0
    if (nprim > np) then
      Grad(np+1:nprim, ia) = p2(np+1:nprim)*irho2 - p1(np+1:nprim)*irho1
      waldis = 0.5d0*(wd1 + wd2)
    end if
    ! --- tangential gradient on axis ib (cells 3,4,5,6) ---
    irt1=1d0/rho_n10(3); irt2=1d0/rho_n10(4); irt3=1d0/rho_n10(5); irt4=1d0/rho_n10(6)
    ! STREAM cells 3..6: load w(c3..c6,s) at use -> folds into the FMA and dies. Bit-identical
    ! (same address/value/order as the p3..p6 arrays); frees ~40 regs (no live p3..p6).
    do s = 1, nsc
      Grad(s, ib) = (w(c4,s)*irt2 - w(c3,s)*irt1 + w(c6,s)*irt4 - w(c5,s)*irt3)*0.25d0
    end do
    Grad(np, ib)    = (T_n10(4) - T_n10(3) + T_n10(6) - T_n10(5))*0.25d0
    do s = nu, nw
      Grad(s, ib) = (w(c4,s) - w(c3,s) + w(c6,s) - w(c5,s))*0.25d0
    end do
    if (nprim > np) then
      do s = np+1, nprim
        Grad(s, ib) = (w(c4,s)*irt2 - w(c3,s)*irt1 + w(c6,s)*irt4 - w(c5,s)*irt3)*0.25d0
      end do
    end if
    ! --- tangential gradient on axis ic (cells 7,8,9,10) ---
    irt1=1d0/rho_n10(7); irt2=1d0/rho_n10(8); irt3=1d0/rho_n10(9); irt4=1d0/rho_n10(10)
    ! STREAM cells 7..10: load-at-use, frees ~40 regs (no live p7..p10).
    do s = 1, nsc
      Grad(s, ic) = (w(c8,s)*irt2 - w(c7,s)*irt1 + w(c10,s)*irt4 - w(c9,s)*irt3)*0.25d0
    end do
    Grad(np, ic)    = (T_n10(8) - T_n10(7) + T_n10(10) - T_n10(9))*0.25d0
    do s = nu, nw
      Grad(s, ic) = (w(c8,s) - w(c7,s) + w(c10,s) - w(c9,s))*0.25d0
    end do
    if (nprim > np) then
      do s = np+1, nprim
        Grad(s, ic) = (w(c8,s)*irt2 - w(c7,s)*irt1 + w(c10,s)*irt4 - w(c9,s)*irt3)*0.25d0
      end do
    end if
    ! --- metric transform: Grad = Grad . (0.5*(M1+M2)) ---
    Mav = 0.5d0*(M1l + M2l)
    do iv = 1, nprim
      gg1 = Grad(iv,1); gg2 = Grad(iv,2); gg3 = Grad(iv,3)
      Grad(iv,1) = gg1*Mav(1,1) + gg2*Mav(2,1) + gg3*Mav(3,1)
      Grad(iv,2) = gg1*Mav(1,2) + gg2*Mav(2,2) + gg3*Mav(3,2)
      Grad(iv,3) = gg1*Mav(1,3) + gg2*Mav(2,3) + gg3*Mav(3,3)
    end do
    ! --- face-averaged primitive ---
    do iv = 1, nprim
      Primf(iv) = 0.5d0*(p1(iv) + p2(iv))
    end do
    ! ===== Compute_Diffusive_Flux_cached =====
    drho = rho_f; dT = T_f; dcp = cp_f; dmil = mil_f; dkl = kl_f
    Ti = int(dT); Tdf = dT - Ti
    do dd = 1, 3
      Gvel(1,dd) = Grad(nu,dd); Gvel(2,dd) = Grad(nv,dd); Gvel(3,dd) = Grad(nw,dd)
    end do
    mie = 0d0
    if (RANS .and. eddy_dev == 2) then
      ! Spalart_Eddy_Viscosity: mu_t = fv1(chi) * rho*nu-tilde, chi = rho*nu-tilde / mu_l (face values)
      chi_ = Primf(nt)/dmil
      mie  = chi_**3/(chi_**3 + cv1**3) * Primf(nt)
    else if (RANS .and. eddy_dev == 3) then
      ! Wilcox2006_Eddy_Viscosity: mu_t = rho*k / max(omega, Clim sqrt(2 Sbar:Sbar / beta*)), Sbar = S - div(u)/3 I
      dome  = Primf(nt+1)/drho
      wdiag = 1d0/3d0 * (Gvel(1,1) + Gvel(2,2) + Gvel(3,3))
      Sij(1,1)=Gvel(1,1)-wdiag; Sij(2,2)=Gvel(2,2)-wdiag; Sij(3,3)=Gvel(3,3)-wdiag
      Sij(1,2)=0.5d0*(Gvel(1,2)+Gvel(2,1)); Sij(1,3)=0.5d0*(Gvel(1,3)+Gvel(3,1)); Sij(2,3)=0.5d0*(Gvel(2,3)+Gvel(3,2))
      Sij(2,1)=Sij(1,2); Sij(3,1)=Sij(1,3); Sij(3,2)=Sij(2,3)
      mie = Primf(nt) / max( dome, wx_Clim*sqrt(2d0*sum(Sij**2)/wx_beta_star) )
    else if (RANS) then
      ! SST_Eddy_Viscosity (vorticity O computed but unused upstream -> skipped, bit-exact)
      dkap = Primf(nt)/drho; dome = Primf(nt+1)/drho
      arg2 = max( 2d0*sqrt(dkap)/(beta_star*dome*waldis), 5d2*dmil/(drho*dome*waldis**2) )
      dF2 = tanh(arg2**2)
      Sij(1,1)=Gvel(1,1); Sij(2,2)=Gvel(2,2); Sij(3,3)=Gvel(3,3)
      Sij(1,2)=0.5d0*(Gvel(1,2)+Gvel(2,1)); Sij(1,3)=0.5d0*(Gvel(1,3)+Gvel(3,1)); Sij(2,3)=0.5d0*(Gvel(2,3)+Gvel(3,2))
      Sij(2,1)=Sij(1,2); Sij(3,1)=Sij(1,3); Sij(3,2)=Sij(2,3)
      Smag = sqrt(2d0*sum(Sij**2))
      mie = drho*dkap*alpha_1 / max(alpha_1*dome, Smag*dF2)
    end if
#ifdef MOSE_MCDIFF_KERNEL
    ! DIFF_MC twin: Dm = the mixture-averaged D_k of the face (co_DS_expr_dev on the face state, written by the fill into
    ! column 5+s) + mu_t/(rho Sct), as MOSE_open with schmidt <= 0 (the face row of the kernel: tface(f,.) or tfg(t_,.)).
    do s = 1, nsc
      Dm(s) = DCUF_TF(DCUF_ROW, 5+s)
    end do
    if (mie > 0d0) then
      do s = 1, nsc
        Dm(s) = Dm(s) + mie/(drho*sct_dev)
      end do
    end if
#else
    Dm(1:nsc) = (dmil/sc_dev + mie/sct_dev)/drho
#endif
    kappaf = dkl + mie*dcp/prt_dev
    ! Stress vector (Std, or 2eq adds the SST k diagonal)
    dmi = dmil + mie
    dlam = negativetwothird*dmi
    Diver = Gvel(1,1) + Gvel(2,2) + Gvel(3,3)
    lamDiv = dlam*Diver
    if (stress_dev == 2) then
      rhokd = negativetwothird*Primf(nt)
    else
      rhokd = 0d0
    end if
    Txx = lamDiv + 2d0*dmi*Gvel(1,1) + rhokd
    Tyy = lamDiv + 2d0*dmi*Gvel(2,2) + rhokd
    Tzz = lamDiv + 2d0*dmi*Gvel(3,3) + rhokd
    Txy = dmi*(Gvel(1,2)+Gvel(2,1)); Txz = dmi*(Gvel(1,3)+Gvel(3,1)); Tyz = dmi*(Gvel(2,3)+Gvel(3,2))
    Strs(1) = Txx*nrm(1) + Txy*nrm(2) + Txz*nrm(3)
    Strs(2) = Txy*nrm(1) + Tyy*nrm(2) + Tyz*nrm(3)
    Strs(3) = Txz*nrm(1) + Tyz*nrm(2) + Tzz*nrm(3)
    ! species + enthalpy diffusion
    DiffH = 0d0
    do iv = 1, nprim
      Flx(iv) = 0d0
    end do
    do s = 1, nsc
      DmGY = Dm(s)*(Grad(s,1)*nrm(1) + Grad(s,2)*nrm(2) + Grad(s,3)*nrm(3))
      Flx(s) = drho*area*DmGY
      hh1 = h_tab_d(Ti, s); hh2 = h_tab_d(Ti+1, s)
      hh  = hh1 + Tdf*(hh2 - hh1)
      DiffH = DiffH + DmGY*hh
    end do
    sflx = 0d0
    do s = 1, nsc; sflx = sflx + Flx(s); end do
    do s = 1, nsc; Flx(s) = Flx(s) - sflx*Primf(s)/drho; end do
    Flx(nu) = area*Strs(1); Flx(nv) = area*Strs(2); Flx(nw) = area*Strs(3)
    Flx(np) = area*( (Strs(1)*Primf(nu) + Strs(2)*Primf(nv) + Strs(3)*Primf(nw)) &
            + kappaf*(Grad(np,1)*nrm(1) + Grad(np,2)*nrm(2) + Grad(np,3)*nrm(3)) + DiffH*drho )
    ! RANS diffusive flux: SA (Spalart_RANS_Diffusive_Flux), Wilcox2006 (Wilcox2006_RANS_Diffusive_Flux: rho k / omega
    ! with the unlimited omega) or SST (mit == mie, same inputs; F1 + sigma blending)
    if (RANS .and. eddy_dev == 2) then
      gndot   = Grad(nt,1)*nrm(1) + Grad(nt,2)*nrm(2) + Grad(nt,3)*nrm(3)
      Flx(nt) = (dmil + Primf(nt))*gndot*area/sa_sigma
    else if (RANS .and. eddy_dev == 3) then
      gkdot = Grad(nt,1)*nrm(1)   + Grad(nt,2)*nrm(2)   + Grad(nt,3)*nrm(3)
      gwdot = Grad(nt+1,1)*nrm(1) + Grad(nt+1,2)*nrm(2) + Grad(nt+1,3)*nrm(3)
      Flx(nt)   = (dmil + wx_sigma_k*Primf(nt)/dome)*gkdot*area
      Flx(nt+1) = (dmil + wx_sigma_w*Primf(nt)/dome)*gwdot*area
    else if (RANS) then
      gkdot = Grad(nt,1)*nrm(1)   + Grad(nt,2)*nrm(2)   + Grad(nt,3)*nrm(3)
      gwdot = Grad(nt+1,1)*nrm(1) + Grad(nt+1,2)*nrm(2) + Grad(nt+1,3)*nrm(3)
      gkw   = Grad(nt+1,1)*Grad(nt,1) + Grad(nt+1,2)*Grad(nt,2) + Grad(nt+1,3)*Grad(nt,3)
      CDw   = sigmad_e*drho/dome * max(gkw, 1d-20)
      zeta  = min( max( sqrt(dkap)/(C_mu*dome*waldis), 5d2*dmil/(drho*dome*waldis**2) ), &
                   4d0*sigmaw_e*drho*dkap/(CDw*waldis**2) )
      dF1   = tanh(zeta**4)
      sigk  = sigma_k1*dF1 + sigma_k2*(1d0-dF1)
      sigw  = sigma_w1*dF1 + sigma_w2*(1d0-dF1)
      Flx(nt)   = (dmil + sigk*mie)*gkdot*area
      Flx(nt+1) = (dmil + sigw*mie)*gwdot*area
    end if
    res0(1:nprim) = Flx(1:nprim)
