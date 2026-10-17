! One cell of the SST source (SST_Blk): gradients, blending, production,
! cross diffusion, dissipation, the var-last R_flat write; shared so the box and the list sweeps match.
      ! Cell-center thermo from cache
      rho  = rho_c_f(WR(i,j,k))   ! the flat twins (the same doubles the per-block cache holds)
      mil  = mil_c_f(WR(i,j,k))
      kap  = PS(nt  ,i,j,k) / rho
      ome  = PS(nt+1,i,j,k) / rho
      dist = WDist(i,j,k)

#ifdef MOSE_RC_KERNEL
      ! SST-RC: the gradient of Compute_Velocity_Gradient, 1d-40 included (MOSE_open: Gradvel = gradv(i,j,k)%c)
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
      Divel = Gradvel(1,1) + Gradvel(2,2) + Gradvel(3,3)

      ! Gradient of (k,w) -> metric: Grad = matmul(Grw, Metric)
      sum_im = PSUM(i-1,j,k) ; sum_ip = PSUM(i+1,j,k)
      sum_jm = PSUM(i,j-1,k) ; sum_jp = PSUM(i,j+1,k)
      sum_km = PSUM(i,j,k-1) ; sum_kp = PSUM(i,j,k+1)
      Grw(1,1) = ( PS(nt  ,i+1,j,k)/sum_ip - PS(nt  ,i-1,j,k)/sum_im ) * 0.5d0
      Grw(2,1) = ( PS(nt+1,i+1,j,k)/sum_ip - PS(nt+1,i-1,j,k)/sum_im ) * 0.5d0
      Grw(1,2) = ( PS(nt  ,i,j+1,k)/sum_jp - PS(nt  ,i,j-1,k)/sum_jm ) * 0.5d0
      Grw(2,2) = ( PS(nt+1,i,j+1,k)/sum_jp - PS(nt+1,i,j-1,k)/sum_jm ) * 0.5d0
      Grw(1,3) = ( PS(nt  ,i,j,k+1)/sum_kp - PS(nt  ,i,j,k-1)/sum_km ) * 0.5d0
      Grw(2,3) = ( PS(nt+1,i,j,k+1)/sum_kp - PS(nt+1,i,j,k-1)/sum_km ) * 0.5d0
      do c = 1, 3
        Grad(1,c) = Grw(1,1)*Metric(1,c,i,j,k) + Grw(1,2)*Metric(2,c,i,j,k) + Grw(1,3)*Metric(3,c,i,j,k)
        Grad(2,c) = Grw(2,1)*Metric(1,c,i,j,k) + Grw(2,2)*Metric(2,c,i,j,k) + Grw(2,3)*Metric(3,c,i,j,k)
      end do

      ! Blending (gather Grad rows to contiguous temps for compute_F)
      gradk_c(1)=Grad(1,1); gradk_c(2)=Grad(1,2); gradk_c(3)=Grad(1,3)
      gradw_c(1)=Grad(2,1); gradw_c(2)=Grad(2,2); gradw_c(3)=Grad(2,3)
      F = compute_F ( rho=rho, w=ome, k=kap, mil=mil, gradk=gradk_c, gradw=gradw_c, d=dist )
      gamma = gamma_1*F(1) + gamma_2*(1d0-F(1))
      beta  = beta_1 *F(1) + beta_2 *(1d0-F(1))

      ! Strain magnitude S = sqrt(2*sum(Sij**2))
      Sij = Strain_Tensor ( Gradvel )
      S2 = 0d0
      do c = 1, 3
      do a = 1, 3
        S2 = S2 + Sij(a,c)**2
      end do
      end do
      S = Sqrt ( 2d0 * S2 )

      ! Vorticity magnitude (absolute vorticity in a rotating frame)
      vort = Vorticity_Vector ( Gradvel )
      if (rotating_frame) then
        vort(1) = vort(1) + 2.0d0 * ov1_
        vort(2) = vort(2) + 2.0d0 * ov2_
        vort(3) = vort(3) + 2.0d0 * ov3_
      end if
      Wij = Vorticity_Tensor ( vort )
      O2 = 0d0
      do c = 1, 3
      do a = 1, 3
        O2 = O2 + Wij(a,c)**2
      end do
      end do
      O = Sqrt ( 2d0 * O2 )

      mi_t = rho*kap*alpha_1 / Max ( alpha_1*ome, S*F(2) )

      ! Production (SST-2003m), with k-eqn limiter; sst-production = compressible (SST-2003e, MOSE_open Lib_SST): the
      ! limiter acts on the deviatoric part mu_t (S^2 - 2/3 divU^2), -2/3 rho k divU and -2/3 gamma rho omega divU outside it
      Prod(1) = mi_t * S**2
      if ( sst_comp_c(1) ) Prod(1) = Prod(1) - 2d0/3d0 * mi_t * Divel**2
#ifdef MOSE_RC_KERNEL
      ! Rotation/curvature correction, Smirnov & Menter (2009), TMR SST-RC (MOSE_open Lib_SST SST_Blk)
      rstar  = rc1_f(WR(i,j,k))
      Drc    = sqrt( max ( S**2, 0.09d0*ome**2) )
      rtilde = rc2_f(WR(i,j,k)) / ( max(O, 1d-30)*Drc**3 )
      frot   = (1.d0 + cr1)*(2.d0*rstar/(1.d0 + rstar))*(1.d0 - cr3*atan(cr2*rtilde)) - cr1
      fr1    = max ( min(frot, 1.25d0), 0d0 )
      Prod(1) = fr1*Prod(1)
#endif
      Prod(1) = Min( Prod(1), 10d0*beta_star*rho*ome*kap )
      Prod(2) = ( gamma*rho / mi_t ) * Prod(1)
      if ( sst_comp_c(1) ) then
        Prod(1) = Prod(1) - 2d0/3d0 * rho * kap * Divel
        Prod(2) = Prod(2) - 2d0/3d0 * gamma * rho * ome * Divel
      end if

      ! Cross diffusion
      dkDotdw = Grad(1,1)*Grad(2,1) + Grad(1,2)*Grad(2,2) + Grad(1,3)*Grad(2,3)
      Diff = 2d0*(1d0-F(1))*rho*sigma_w2/ome*dkDotdw

      ! Dissipation
      Diss(1) = beta_star*rho*kap*ome
      Diss(2) = beta*rho*ome**2

      ! Source terms
      Source(1) = Prod(1) - Diss(1)
      Source(2) = Prod(2) - Diss(2) + Diff

      ! SST k/omega source -> coalesced var-last R_flat (WR = flat cidx,
      ! the same macro the w-reads use). Completes the W-path residual writers on R_flat.
      if (pimp) then   ! Patankar treatment of the destruction terms (MOSE_open Lib_SST), box sweep only
        fk_ = 1d0 / ( 1d0 + dtl(i,j,k) * beta_star * ome )
        fw_ = 1d0 / ( 1d0 + dtl(i,j,k) * 2d0 * beta * ome )
        R_flat(WR(i,j,k), nt)   = R_flat(WR(i,j,k), nt)   - Source(1) * Volume(i,j,k) * fk_
        R_flat(WR(i,j,k), nt+1) = R_flat(WR(i,j,k), nt+1) - Source(2) * Volume(i,j,k) * fw_
      else
        fk_ = 1d0
        R_flat(WR(i,j,k), nt:nt+1) = R_flat(WR(i,j,k), nt:nt+1) - Source * Volume(i,j,k)
      end if
      ! k in the total energy: the k source leaves the mean-flow energy with the Patankar factor of the k equation,
      ! as MOSE_open Lib_SST SST_Blk (Res(np) + Source(1) * Volume * fk; fk = 1 without point-implicit): A10 of
      ! 2026-10-16. x * 1d0 is exact: the bits without point-implicit are those of before.
      if (k_energy_coupling) R_flat(WR(i,j,k), np) = R_flat(WR(i,j,k), np) + Source(1) * Volume(i,j,k) * fk_
