! 5s-a: c0 and the three per-block values (dt, vol, T seed) come from the caller, so the same
! body serves the per-block kernel and a kernel fused over the level's concatenated cell list.
! One cell of the conservative Newstate (Newstate_Conservative_blk / _list): prim2cons of
! PO_flat and w, the RK update on R_flat, cons2prim with the T_c seed, the clamps, the var-last scatter.
      do s = 1, nprim
        p_l(s)  = P_flat(c0,s)
        po_l(s) = PO_flat(c0,s)
        res_l(s)= R_flat(c0,s)
      end do
      dt = dt_l

      ! prim2cons(PO) -> consold
      rho  = sum( po_l(1:nsc) )
      Rg   = f_Rtot( po_l(1:nsc) )
      vmag = 0d0; do s = nu, nw; vmag = vmag + po_l(s)*po_l(s); end do
      Tloc = po_l(np) / (Rg*rho)
      T_i = int(Tloc); if (T_i < tmin_) T_i = tmin_; if (T_i > tmax_-1) T_i = tmax_-1; Td = Tloc - T_i
      hmix = 0d0
      do s = 1, nsc
        h_s  = h_tab(T_i,s) + Td*( h_tab(T_i+1,s) - h_tab(T_i,s) )
        hmix = hmix + h_s * po_l(s)/rho
      end do
      E0 = hmix - po_l(np)/rho + 0.5d0*vmag
      consold(1:nsc) = po_l(1:nsc)
      consold(nu:nw) = rho * po_l(nu:nw)
      consold(np)    = rho * E0
      if (nprim>np) consold(np+1:nprim) = po_l(np+1:nprim)

      ! prim2cons(P) -> cons
      rho  = sum( p_l(1:nsc) )
      Rg   = f_Rtot( p_l(1:nsc) )
      vmag = 0d0; do s = nu, nw; vmag = vmag + p_l(s)*p_l(s); end do
      Tloc = p_l(np) / (Rg*rho)
      T_i = int(Tloc); if (T_i < tmin_) T_i = tmin_; if (T_i > tmax_-1) T_i = tmax_-1; Td = Tloc - T_i
      hmix = 0d0
      do s = 1, nsc
        h_s  = h_tab(T_i,s) + Td*( h_tab(T_i+1,s) - h_tab(T_i,s) )
        hmix = hmix + h_s * p_l(s)/rho
      end do
      E0 = hmix - p_l(np)/rho + 0.5d0*vmag
      cons(1:nsc) = p_l(1:nsc)
      cons(nu:nw) = rho * p_l(nu:nw)
      cons(np)    = rho * E0
      if (nprim>np) cons(np+1:nprim) = p_l(np+1:nprim)

      ! time residual + RK stage (conservative form)
      if (.not. ps_) res_l(1:nprim) = -res_l(1:nprim)/vol_l*dt*strang   ! IRS: already scaled and smoothed (Mod_IRS_W)
      cons(1:nprim)  = consold(1:nprim) &
                     + ( cons(1:nprim) - consold(1:nprim) + res_l(1:nprim) ) * rkc

      ! cons2prim (Newton-Raphson for temperature, T_c seed)
      rho  = sum( cons(1:nsc) )
      Rg   = f_Rtot( cons(1:nsc) )
      p_l(1:nsc) = cons(1:nsc)
      p_l(nu:nw) = cons(nu:nw)/rho
      vmag   = 0d0; do s = nu, nw; vmag = vmag + (cons(s)/rho)*(cons(s)/rho); end do
      energy = cons(np)/rho - 0.5d0*vmag
      diff_en = 1d3; it = 1; T_it = tc_l
      do while ( (abs(diff_en/energy) >= 1d-6) .and. (it <= 20) )
        T_i = int(T_it); if (T_i < tmin_) T_i = tmin_; if (T_i > tmax_-1) T_i = tmax_-1; Td = T_it - T_i   ! NaN / off-table: no wild index, the thermo guard stops next
        cv = -Rg; enth = 0d0
        do s = 1, nsc
          cp_s = cp_tab(T_i,s) + Td*( cp_tab(T_i+1,s) - cp_tab(T_i,s) )
          h_s  = h_tab(T_i,s)  + Td*( h_tab(T_i+1,s)  - h_tab(T_i,s)  )
          cv   = cv   + cons(s)/rho * cp_s
          enth = enth + cons(s)/rho * h_s
        end do
        p_it    = rho*Rg*T_it
        diff_en = abs( (enth - p_it/rho) - energy )
        T_it    = T_it + (energy - (enth - p_it/rho))/cv
        it      = it + 1
      end do
      p_l(np) = rho*Rg*T_it
      if (nprim>np) p_l(np+1:nprim) = cons(np+1:nprim)

      ! check_and_fix_state clamps (species, SST k/w)
      do s = 1, nsc
        if (p_l(s) < 1d-20) p_l(s) = 1d-20
      end do
      if (RANS) then   ! the floors of SST_/Wilcox2006_/Spalart_Enforce_Realizability on the nrans turbulence variables
        do s = nt, nt + nrans - 1
          if (p_l(s) < 1d-10) p_l(s) = 1d-10
        end do
      end if

      ! scatter the updated primitive state back var-last
      do s = 1, nprim
        P_flat(c0,s) = p_l(s)
      end do
