! One cell of the Wilcox (2006) k-omega source (Wilcox_Blk): MOSE_open b29130b Lib_Wilcox2006 Wilcox_Blk statement by
! statement on the var-last w (PS/PV/PSUM, the macros of Mod_Wilcox2006.f90), the var-last R_flat write.
      rho = PSUM(i,j,k)
      kap = PS(nt  ,i,j,k) / rho
      ome = PS(nt+1,i,j,k) / rho

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

      ! Strain tensor; Sb2 = the sum of the squares of its deviatoric part (MOSE_open Sum(Sij_bar**2), same order)
      Divel = Gradvel(1,1) + Gradvel(2,2) + Gradvel(3,3)
      Sij = Strain_Tensor ( Gradvel )
      diag = 1d0/3d0 * Divel
      Sb2 = 0d0
      do c = 1, 3
      do a = 1, 3
        if ( a == c ) then
          Sb2 = Sb2 + ( Sij(a,c) - diag )**2
        else
          Sb2 = Sb2 + Sij(a,c)**2
        end if
      end do
      end do

      ! Eddy viscosity with the stress limiter
      ome_hat = Max ( ome, Clim * Sqrt( 2d0 * Sb2 / beta_star ) )
      mi_t = rho * kap / ome_hat

      ! Production P = tau_ij du_i/dx_j, tau_ij = 2 mu_t S_ij - (2/3 mu_t div(u) + 2/3 rho k) delta_ij
      diag = 2d0/3d0 * mi_t * Divel + 2d0/3d0 * rho * kap
      Prod(1) = 0d0
      do a = 1, 3
      do c = 1, 3
        Tij_ = 2d0 * mi_t * Sij(a,c)
        if ( a == c ) Tij_ = Tij_ - diag
        Prod(1) = Prod(1) + Tij_ * Gradvel(a,c)
      end do
      end do
      Prod(2) = ( gamma * ome / kap ) * Prod(1)
      Prod(1) = Min( Prod(1), 20d0*beta_star*rho*ome*kap ) ! k-eqn production limiter

      ! Gradient of k and omega (k = rho*k / sum of the species densities) -> metric: Grad = matmul(Grw, Metric)
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

      ! Cross diffusion, only where grad k . grad omega > 0
      dkDotdw = Grad(1,1)*Grad(2,1) + Grad(1,2)*Grad(2,2) + Grad(1,3)*Grad(2,3)
      sigma_d = 1d0/8d0
      if ( dkDotdw <= 0d0 ) sigma_d = 0d0
      Diff = sigma_d * rho / ome * dkDotdw

      ! Dissipation: beta = beta_0 f_beta(X_w), X_w = |W_ij W_jk S^_ki / (beta* omega)^3| with S^ = S - div(u)/2 I on all
      ! nine entries (MOSE_open assigns the diagonal only: see the module header) and W of the absolute vorticity
      do c = 1, 3
      do a = 1, 3
        Sij_hat(a,c) = Sij(a,c)
      end do
      end do
      Sij_hat(1,1) = Sij(1,1) - 0.5d0 * Divel
      Sij_hat(2,2) = Sij(2,2) - 0.5d0 * Divel
      Sij_hat(3,3) = Sij(3,3) - 0.5d0 * Divel
      vort = Vorticity_Vector ( Gradvel )
      if ( rotating_frame ) then
        vort(1) = vort(1) + 2.0d0 * ov1_
        vort(2) = vort(2) + 2.0d0 * ov2_
        vort(3) = vort(3) + 2.0d0 * ov3_
      end if
      Wij = Vorticity_Tensor ( vort )
      Num = 0d0
      do a = 1, 3
      do c = 1, 3
      do kk = 1, 3
        Num = Num + Wij(a,c) * Wij(c,kk) * Sij_hat(kk,a)
      end do
      end do
      end do
      X_w = Abs ( Num / ( beta_star * ome )**3 )
      f_beta = ( 1d0 + 85d0 * X_w ) / ( 1d0 + 100d0 * X_w )
      beta = beta_0 * f_beta

      Diss(1) = beta_star * rho * kap * ome
      Diss(2) = beta * rho * ome**2

      ! Source terms
      Source(1) = Prod(1) - Diss(1)
      Source(2) = Prod(2) - Diss(2) + Diff

      ! Point-implicit (Patankar) treatment of the destruction terms
      fk = 1d0 ; fw = 1d0
      if ( pimp ) then
        fk = 1d0 / ( 1d0 + dtl(i,j,k) * beta_star * ome )
        fw = 1d0 / ( 1d0 + dtl(i,j,k) * 2d0 * beta * ome )
      end if

      R_flat(WR(i,j,k), nt)   = R_flat(WR(i,j,k), nt)   - Source(1) * Volume(i,j,k) * fk
      R_flat(WR(i,j,k), nt+1) = R_flat(WR(i,j,k), nt+1) - Source(2) * Volume(i,j,k) * fw
      ! k in the total energy: the k source leaves the mean-flow energy with the Patankar factor of the k equation
      if ( k_energy_coupling ) R_flat(WR(i,j,k), np) = R_flat(WR(i,j,k), np) + Source(1) * Volume(i,j,k) * fk
