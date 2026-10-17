    ! ===== inlined-physics locals (conv_decls.inc) =====
    real(R8) :: l0, l1, l2, lm, lp, beta_order, idl0, idl1, idl2
    real(R8) :: slope0, slope1, slope2, slopel, sloper, xx, yy, mm
    real(R8) :: Prim1(NPRIM_DEV), Prim4(NPRIM_DEV)
    real(R8) :: rotot1, Rtot1, rotot4, Rtot4, a1, a4, irotot1, irotot4
    real(R8) :: T1d, T4d, Tdiffd, cpd, cpi, gamd
    integer  :: Ti1, Ti4, Til, Tir, ss
    real(R8) :: F_r, F_u, F_v, F_w, F_E, su, sp, sm
    real(R8) :: nx, ny, nz, unl, unr, un_roe
    real(R8) :: rho_roe, u_roe, v_roe, w_roe, a_roe, h0_roe, h0l, h0r
    real(R8) :: Rl, Rr, Tl, Tr, dTl, dTr, inv_dltot, inv_drtot, srL, srR, inv_sr
    real(R8) :: cv_roe, sum_ei, R_roe, T_roe, gam_roe, vel2, hl, hr, el, er, invW, d_roe
    real(R8) :: dltot, drtot, ul, vl, wl, ur, vr, wr, pl, pr, al, ar
    real(R8) :: S1, S4, Sstar, pstar, inv_denom, U1S, U2S, U3S, U4S, U5S, e0l, e0r
    real(R8) :: frl_, frr_, ml_, mr_, mt_, chi_, vnm_, frho_, vnl_, vnr_, bl_, br_, af_, ph_   ! HLLE / SLAU (rsv = 2 / 3)
    real(R8) :: urL_, urR_   ! HLLC-PC (14): the Ur of the two cells of the face
