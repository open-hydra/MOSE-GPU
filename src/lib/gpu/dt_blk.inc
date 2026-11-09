! ---------------------------------------------------------------------------------------------
! Compute_dt_blk, instantiated twice in Mod_dt.f90:
!   - MOSE_MCDIFF_KERNEL undefined -> Compute_dt_blk, the previous source token for token;
!   - MOSE_MCDIFF_KERNEL defined   -> Compute_dt_blk_mc, two more arguments (vnnd = [MOSE-parameters] vnn-diff, sctd = Sct)
!     and the species-diffusion limit A' = vnnd * dl_f(d)^2 / max_k(D_k + mu_t/(rho Sct)) after the VNN one (D7).
! Expected macro: the routine name (DTW_NAME).
! Both take dt-rule (additive) and the stiffness coefficient of explicit turbulence sources (stiffc): with additive = .false.
! the step is the previous one, token for token; with additive = .true. it is the rule of MOSE_open main (Mod_dt.f90,
! compute), below.
! ---------------------------------------------------------------------------------------------
  subroutine DTW_NAME ( b, m_f, dl_f, rho_c, Rg_c, mil_c, &
#ifdef MOSE_MCDIFF_KERNEL
                                dtlocal, n, cfl, vnn, rampa_fac, dtglobal, urs, nur, prec, additive, stiffc, vnnd, sctd )
#else
                                dtlocal, n, cfl, vnn, rampa_fac, dtglobal, urs, nur, prec, additive, stiffc )
#endif
    use MOSE_Field_m, only: dt_f, n_cache_f
    use MOSE_Mod_SST, only: SST_Eddy_Viscosity
    use MOSE_Mod_Spalart, only: Spalart_Eddy_Viscosity
    use MOSE_Mod_Wilcox2006, only: Wilcox2006_Eddy_Viscosity
    use MOSE_Lib_RANS, only: eddy_dev, sc_dev, sct_dev, prt_dev, prl_dev, mcdiff_dev
#ifdef MOSE_MCDIFF_KERNEL
    use FLINT_Lib_ThermoTransport_dev, only: f_ss_dev, co_DS_expr_dev
#else
    use FLINT_Lib_ThermoTransport_dev, only: f_ss_dev
#endif
    use FLINT_Lib_ThermoTransport_dev, only: f_cp_dev, co_k_mi_lam_Wilke_dev, co_DS_expr_dev
    integer,  intent(in) :: b, n(3)
    real(R8), intent(in) :: cfl, vnn, rampa_fac
    real(R8), intent(inout) :: dtglobal
    integer,  intent(in) :: nur
    real(R8), intent(in) :: urs(nur)   ! Lib_Preconditioning ur_f: Ur of the cells (prec), 1 element otherwise
    logical,  intent(in) :: prec       ! integration-variables = prec
    logical,  intent(in) :: additive   ! dt-rule = additive (MOSE_open main)
    real(R8), intent(in) :: stiffc     ! additive rule: 1/dt += stiffc omega (explicit SST sources), 0 otherwise
    real(R8), dimension(3, 3, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: m_f
    real(R8), dimension(3, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: dl_f
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: rho_c, Rg_c, mil_c
    real(R8), dimension(1:n(1), 1:n(2), 1:n(3)), intent(out) :: dtlocal
    ! Local (per cell)
    integer  :: i, j, k, d, s, c0, boff, pi, pj
    real(R8) :: Sound, dtc, dtd, lambda, vel_, mie, mi, Beta_, Alpha_, cst_
#ifdef MOSE_MCDIFF_KERNEL
    real(R8), intent(in) :: vnnd, sctd   ! DIFF_MC A': vnn-diff and the turbulent Schmidt number
    integer  :: T_i, Tint(2)
    real(R8) :: Tc_, Tdf_, Dmax_, Dk(NPRIM_DEV)
#endif
    ! nprim is a compile-time parameter: exact-size local.
    real(R8) :: pl(nprim), vel(3), versor(3), Gvel0(3,3), rkw(2)
    ! additive rule (MOSE_open main compute)
    real(R8), parameter :: SIGMA_SA_DT = 2d0/3d0   ! diffusion constant of Spalart-Allmaras, as MOSE_open Mod_dt
    integer  :: Tia(2)
    real(R8) :: ra_c, ra_v, gn_, c12_, c13_, c23_, ss_, qp_, mr_, Ta_, Tdfa_, cpa_, cva_, mila_, kla_, Dmx_, Dkm_
    real(R8) :: Dka(NPRIM_DEV)
    call field_block_index_params(b, boff, pi, pj)

    !$acc data copyin(m_f, dl_f, rho_c, Rg_c, mil_c) copyout(dtlocal) copy(dtglobal)
    !$acc parallel loop collapse(3) gang vector reduction(min:dtglobal) &
    !$acc&  present(P_flat, m_f, dl_f, rho_c, Rg_c, mil_c, dtlocal, dt_f, urs) &
#ifdef MOSE_MCDIFF_KERNEL
    !$acc&  private(pl, vel, versor, Gvel0, rkw, c0, s, d, T_i, Tint, Tc_, Tdf_, Dmax_, Dk, Tia, Dka)
#else
    !$acc&  private(pl, vel, versor, Gvel0, rkw, c0, s, d, Tia, Dka)
#endif
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)
      c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      do s = 1, nprim
        pl(s) = P_flat(c0,s)
      end do
      Sound    = f_ss_dev( pl(1:nsc), pl(np), rho_c(i,j,k), Rg_c(i,j,k) )
      vel(1:3) = pl(nu:nw)
      dtc      = 1d8
      if (prec) then   ! MOSE_open Compute_dt with preconditioning
        Beta_  = 1.0d0/Sound**2
        Alpha_ = 0.5d0 * ( 1d0 - Beta_ * urs(c0)**2 )
      end if
      if (additive) then
      ! dt-rule = additive, the rule of MOSE_open main (Mod_dt.f90, compute): the convective and diffusive rates add up,
      !   1/dt = sum_d lambda_d |grad xi_d| / cfl + D_max sum_d |grad xi_d|^2 / vnn  (+ stiffc omega),
      ! grad xi_d the row d of the cell metric; in 3-D the viscous rate takes the skew factor (2 + m + 1/m)/4 on cells
      ! whose three directions are skewed against each other; D_max the largest diffusivity of momentum, energy, species
      ! and nu~ (SA). The wedge term s_h of MOSE_open is not here: the host stops the rule on meshes with ndir /= 3.
        ra_c = 0d0
        ra_v = 0d0
        do d = 1, ndir
          gn_  = norm2( m_f(d,1:3,i,j,k) )
          vel_ = abs( dot_product( vel(1:3), m_f(d,1:3,i,j,k) ) ) / gn_
          if (prec) then
            cst_   = Sqrt ( Alpha_**2 * vel_**2 + urs(c0)**2 )
            lambda = vel_ * ( 1d0 - Alpha_ ) + cst_
          else
            lambda = vel_ + Sound
          end if
          ra_c = ra_c + lambda * gn_
          ra_v = ra_v + gn_**2
        end do
        if (VISC) then
          if (ndir == 3) then
            c12_ = dot_product( m_f(1,1:3,i,j,k), m_f(2,1:3,i,j,k) ) / ( norm2( m_f(1,1:3,i,j,k) ) * norm2( m_f(2,1:3,i,j,k) ) )
            c13_ = dot_product( m_f(1,1:3,i,j,k), m_f(3,1:3,i,j,k) ) / ( norm2( m_f(1,1:3,i,j,k) ) * norm2( m_f(3,1:3,i,j,k) ) )
            c23_ = dot_product( m_f(2,1:3,i,j,k), m_f(3,1:3,i,j,k) ) / ( norm2( m_f(2,1:3,i,j,k) ) * norm2( m_f(3,1:3,i,j,k) ) )
            ss_  = c12_**2 + c13_**2 + c23_**2
            qp_  = c12_ * c13_ * c23_
            if ( ss_ + 2d0 * qp_ > 1d0 ) then
              mr_  = 2d0 * sqrt( ss_ / 3d0 ) * &
                     cos( acos( max( -1d0, min( 1d0, 3d0 * sqrt( 3d0 ) * qp_ / ss_**1.5d0 ) ) ) / 3d0 )
              ra_v = ra_v * 0.25d0 * ( 2d0 + mr_ + 1d0 / mr_ )
            end if
          end if
          Ta_    = pl(np) / ( rho_c(i,j,k) * Rg_c(i,j,k) )
          Tia(1) = int(Ta_)
          Tia(2) = Tia(1) + 1
          Tdfa_  = Ta_ - Tia(1)
          cpa_   = f_cp_dev( pl(1:nsc), Tia, Tdfa_, rho_c(i,j,k) )
          cva_   = cpa_ - Rg_c(i,j,k)
          call co_k_mi_lam_Wilke_dev ( pl(1:nsc), rho_c(i,j,k), Tia, Tdfa_, mila_, kla_ )
          if ( prl_dev > 0d0 ) kla_ = mila_ * cpa_ / prl_dev
          mie = 0d0
          if (RANS) then   ! velocity gradient 0 and wall distance 1e-6, as MOSE_open
            Gvel0 = 0d0
            if (eddy_dev == 2) then
              call Spalart_Eddy_Viscosity ( mie, pl(nt:nt+nrans-1), mila_, rho_c(i,j,k), Gvel0, 1d-6 )
            else
              rkw(1) = pl(nt); rkw(2) = pl(nt+1)
              if (eddy_dev == 3) then
                call Wilcox2006_Eddy_Viscosity ( mie, rkw, mila_, rho_c(i,j,k), Gvel0, 1d-6 )
              else
                call SST_Eddy_Viscosity ( mie, rkw, mila_, rho_c(i,j,k), Gvel0, 1d-6 )
              end if
            end if
          end if
          Dmx_ = max( 4d0/3d0 * ( mila_ + mie ) / rho_c(i,j,k), &
                      ( kla_ + mie * cpa_ / prt_dev ) / ( rho_c(i,j,k) * cva_ ) )
          if ( nsc > 1 ) then
            if ( mcdiff_dev ) then
              call co_DS_expr_dev ( pl(1:nsc), rho_c(i,j,k), Tia, Tdfa_, pl(np), Dka(1:nsc) )
              Dkm_ = Dka(1)
              do s = 2, nsc
                Dkm_ = max( Dkm_, Dka(s) )
              end do
              Dmx_ = max( Dmx_, Dkm_ + mie / ( rho_c(i,j,k) * sct_dev ) )
            else
              Dmx_ = max( Dmx_, ( mila_ / sc_dev + mie / sct_dev ) / rho_c(i,j,k) )
            end if
          end if
          if ( RANS .and. eddy_dev == 2 ) Dmx_ = max( Dmx_, ( mila_ + pl(nt) ) / ( SIGMA_SA_DT * rho_c(i,j,k) ) )
          dtc = 1d0 / ( ra_c / cfl + Dmx_ * ra_v / vnn )
        else
          dtc = cfl / ra_c
        end if
        ! explicit SST sources (time-accurate, point-implicit = false): the bound of their stiffness, MOSE_open
        ! SST_Stiffness_Bound, so that dt d <= 1
        if ( stiffc > 0d0 ) dtc = 1d0 / ( 1d0 / dtc + stiffc * max( pl(nt+1) / rho_c(i,j,k), 0d0 ) )
      else
      ! CFL condition per computational direction
      do d = 1, ndir
        versor(1:3) = m_f(d,1:3,i,j,k) / norm2( m_f(d,1:3,i,j,k) )
        vel_   = abs( dot_product( vel(1:3), versor(1:3) ) )
        if (prec) then   ! the pseudo-acoustic speed c* of the preconditioned system
          cst_   = Sqrt ( Alpha_**2 * vel_**2 + urs(c0)**2 )
          vel_   = vel_ * ( 1d0 - Alpha_ )
          lambda = vel_ + cst_
        else
          lambda = vel_ + Sound
        end if
        dtd = dl_f(d,i,j,k) / lambda * cfl
        dtc    = min( dtd, dtc )
      end do
      ! VNN (viscous) condition
      if (VISC) then
        mie = 0d0
        if (RANS) then
          Gvel0 = 0d0
          if (eddy_dev == 2) then   ! SA: mu_t = fv1(chi)*rho*nu-tilde
            call Spalart_Eddy_Viscosity ( mie, pl(nt:nt+nrans-1), mil_c(i,j,k), rho_c(i,j,k), Gvel0, 1d-6 )
          else
            rkw(1) = pl(nt); rkw(2) = pl(nt+1)   ! k-omega slots (dead branch in an nrans=1 binary; pl is NPRIM_DEV-padded)
            if (eddy_dev == 3) then   ! Wilcox2006: the zero gradient leaves the stress limiter out, mu_t = rho k / omega
              call Wilcox2006_Eddy_Viscosity ( mie, rkw, mil_c(i,j,k), rho_c(i,j,k), Gvel0, 1d-6 )
            else
              call SST_Eddy_Viscosity ( mie, rkw, mil_c(i,j,k), rho_c(i,j,k), Gvel0, 1d-6 )
            end if
          end if
        end if
        mi = mie + mil_c(i,j,k)
        do d = 1, ndir
          dtd = ( rho_c(i,j,k) * dl_f(d,i,j,k)**2 * vnn ) / mi
          dtc = min( dtd, dtc )
        end do
#ifdef MOSE_MCDIFF_KERNEL
        ! DIFF_MC A' (D7): the species-diffusion limit vnnd * dl_f(d)^2 / max_k(D_k + mu_t/(rho Sct)), D_k of the cell state
        ! (co_DS_expr_dev, rescaled to the cell pressure), in every direction.
        Tc_  = pl(np) / ( rho_c(i,j,k) * Rg_c(i,j,k) )
        T_i  = int(Tc_)
        Tdf_ = Tc_ - T_i
        Tint(1) = T_i
        Tint(2) = T_i + 1
        call co_DS_expr_dev ( pl(1:nsc), rho_c(i,j,k), Tint, Tdf_, pl(np), Dk(1:nsc) )
        Dmax_ = Dk(1)
        do s = 2, nsc
          Dmax_ = max( Dmax_, Dk(s) )
        end do
        Dmax_ = Dmax_ + mie/(rho_c(i,j,k)*sctd)
        do d = 1, ndir
          dtd = ( vnnd * dl_f(d,i,j,k)**2 ) / Dmax_
          dtc = min( dtd, dtc )
        end do
#endif
      end if
      end if   ! dt-rule
      dtc = dtc * rampa_fac
      dtlocal(i,j,k) = dtc
#ifndef MOSE_DT_NOTWIN
      if ( n_cache_f > 1 ) dt_f(c0) = dtc
#endif
      dtglobal = min( dtc, dtglobal )
    end do; end do; end do
    !$acc end data

  end subroutine DTW_NAME
