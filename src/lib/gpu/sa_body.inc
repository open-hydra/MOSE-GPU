! One cell of the Spalart-Allmaras source (SA_Source_Blk): the mirror of the host
! SA_Source_Blk (Mod_Spalart.f90): velocity gradient, vorticity, production (SA or SA-R), non-conservative
! diffusion from the nu-tilde gradient, destruction, the var-last R_flat write. Shared by the box, list and
! fused sweeps so that they match. Reads the var-last w through PS/PV/PSUM (see the macros in Mod_Spalart.f90).
      rho  = rho_c_f(WR(i,j,k))   ! the flat twins (the same doubles the per-block cache holds)
      mil  = mil_c_f(WR(i,j,k))
      nit  = PS(nt,i,j,k) / rho   ! nu-tilde
      chi  = PS(nt,i,j,k) / mil   ! Spalart chi = rho*nu-tilde / mu_l
      dist = WDist(i,j,k)

#ifdef MOSE_RC_KERNEL
      ! SA-RC: the gradient of Compute_Velocity_Gradient, 1d-40 included (MOSE_open: Gradvel = gradv(i,j,k)%c)
      do c = 1, 3
      do a = 1, 3
        Gradvel(a,c) = gv_f(WR(i,j,k), a+3*(c-1))
      end do
      end do
#else
      ! Velocity gradient (central differences) -> metric: Gradvel = matmul(Gtmp, Metric)
      Gtmp(1:3,1) = ( PV(i+1,j,k) - PV(i-1,j,k) ) * 0.5d0
      Gtmp(1:3,2) = ( PV(i,j+1,k) - PV(i,j-1,k) ) * 0.5d0
      Gtmp(1:3,3) = ( PV(i,j,k+1) - PV(i,j,k-1) ) * 0.5d0
      do c = 1, 3
      do a = 1, 3
        Gradvel(a,c) = Gtmp(a,1)*Metric(1,c,i,j,k) &
                     + Gtmp(a,2)*Metric(2,c,i,j,k) &
                     + Gtmp(a,3)*Metric(3,c,i,j,k)
      end do
      end do
#endif

      ! Vorticity (absolute vorticity in a rotating frame) and its magnitude
      vort = Vorticity_Vector ( Gradvel )
      if (rotating_frame) then
        vort(1) = vort(1) + 2.0d0 * ov1_
        vort(2) = vort(2) + 2.0d0 * ov2_
        vort(3) = vort(3) + 2.0d0 * ov3_
      end if
      Om = Sqrt ( vort(1)**2 + vort(2)**2 + vort(3)**2 )

      ! Production (SA 2003 "neg" branch of S-tilde as in the host), fv1/fv2 inlined
      fv1_ = chi**3 / ( chi**3 + cv1**3 )
      fv2_ = 1d0 - chi / ( 1d0 + chi * fv1_ )
      Sbar = nit / ( cK_sq * dist**2 ) * fv2_   ! cK_sq = cK**2: the list kernels' dummy ck(:) hides the constant cK
      if ( Sbar >= (-c2*Om) ) then
        Stilde = Om + Sbar
      else
        Stilde = Om + Om * (c2**2*Om + c3*Sbar) / ((c3 - 2*c2)*Om - Sbar)
      end if
      if ( Stilde == 0d0 ) Stilde = 1d-40   ! just to avoid a division by zero
#ifdef MOSE_RC_KERNEL
      ! SA-RC: the Spalart-Shur factor fr1 on the production (MOSE_open Lib_Spalart SA_Source_Blk)
      Sij = Strain_Tensor ( Gradvel )
      S2 = 0d0
      do c = 1, 3
      do a = 1, 3
        S2 = S2 + Sij(a,c)**2
      end do
      end do
      abs_Eij = Sqrt ( 2d0 * S2 )
      rstar  = rc1_f(WR(i,j,k))
      D2     = 0.5d0*(abs_Eij**2 + Om**2)
      rtilde = rc2_f(WR(i,j,k)) / D2**2
      fr1    = (1.d0 + cr1)*(2.d0*rstar/(1.d0 + rstar))*(1.d0 - cr3*atan(cr2*rtilde)) - cr1
      Production = fr1 * cb1 * Stilde * nit
#else
      if ( sar ) then
        Sij = Strain_Tensor ( Gradvel )
        S2 = 0d0
        do c = 1, 3
        do a = 1, 3
          S2 = S2 + Sij(a,c)**2
        end do
        end do
        abs_Eij = Sqrt ( 2d0 * S2 )
        Production = cb1 * ( Stilde + Crot * Min(0d0, abs_Eij - Om) ) * nit
      else
        Production = cb1 * Stilde * nit
      end if
#endif

      ! nu-tilde gradient (nu-tilde = rho*nu-tilde / sum of the species densities) -> metric
      sum_im = PSUM(i-1,j,k) ; sum_ip = PSUM(i+1,j,k)
      sum_jm = PSUM(i,j-1,k) ; sum_jp = PSUM(i,j+1,k)
      sum_km = PSUM(i,j,k-1) ; sum_kp = PSUM(i,j,k+1)
      Gn(1) = ( PS(nt,i+1,j,k)/sum_ip - PS(nt,i-1,j,k)/sum_im ) * 0.5d0
      Gn(2) = ( PS(nt,i,j+1,k)/sum_jp - PS(nt,i,j-1,k)/sum_jm ) * 0.5d0
      Gn(3) = ( PS(nt,i,j,k+1)/sum_kp - PS(nt,i,j,k-1)/sum_km ) * 0.5d0
      do c = 1, 3
        Gradnit(c) = Gn(1)*Metric(1,c,i,j,k) + Gn(2)*Metric(2,c,i,j,k) + Gn(3)*Metric(3,c,i,j,k)
      end do
      Gradnit2  = Gradnit(1)**2 + Gradnit(2)**2 + Gradnit(3)**2
      Diffusion = cb2 / sigma * Gradnit2   ! non-conservative diffusion term

      ! Destruction
      r_Fun = nit / ( Stilde * cK_sq * dist**2 )
      r_Fun = Min ( r_Fun, 10d0 )
      g_Fun = r_Fun + cw2 * ( r_Fun**6 - r_Fun )
      fw    = g_Fun * ( ( 1d0 + cw3**6 ) / ( g_Fun**6 + cw3**6 ) )**(1.d0/6.d0)
      Destruction = cw1 * fw * ( nit / dist )**2
#ifdef MOSE_SACOMP_KERNEL
      ! SAcomp: the Paciorri-Sabetta compressibility correction of the production (MOSE_open Lib_Spalart), the sound
      ! speed of the cell from its species densities, pressure, density and gas constant
      do s = 1, nsc
        rhoi_(s) = PS(s,i,j,k)
      end do
      sound_ = f_ss_dev ( rhoi_, PS(np,i,j,k), rho, Rg_c_f(WR(i,j,k)) )
      call SA_Compressibility_Correction ( nit, chi, sound_, Om, Production )
#endif

      ! Source multiplied by rho: the SA equation is integrated in conservative form
      Source = rho * ( Production + Diffusion - Destruction ) * Volume(i,j,k)
      if (pimp) then   ! Patankar treatment of the destruction term (MOSE_open Lib_Spalart), box sweep only
        fnu_ = 1d0 / ( 1d0 + dtl(i,j,k) * 2d0 * cw1 * fw * nit / dist**2 )
        Source = Source * fnu_
      end if
      R_flat(WR(i,j,k), nt) = R_flat(WR(i,j,k), nt) - Source
