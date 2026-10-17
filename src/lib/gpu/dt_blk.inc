! ---------------------------------------------------------------------------------------------
! Compute_dt_blk, instantiated twice in Mod_dt.f90:
!   - MOSE_MCDIFF_KERNEL undefined -> Compute_dt_blk, the previous source token for token;
!   - MOSE_MCDIFF_KERNEL defined   -> Compute_dt_blk_mc, two more arguments (vnnd = [MOSE-parameters] vnn-diff, sctd = Sct)
!     and the species-diffusion limit A' = vnnd * dl_f(d)^2 / max_k(D_k + mu_t/(rho Sct)) after the VNN one (D7).
! Expected macro: the routine name (DTW_NAME).
! ---------------------------------------------------------------------------------------------
  subroutine DTW_NAME ( b, m_f, dl_f, rho_c, Rg_c, mil_c, &
#ifdef MOSE_MCDIFF_KERNEL
                                dtlocal, n, cfl, vnn, rampa_fac, dtglobal, urs, nur, prec, vnnd, sctd )
#else
                                dtlocal, n, cfl, vnn, rampa_fac, dtglobal, urs, nur, prec )
#endif
    use MOSE_Field_m, only: dt_f, n_cache_f
    use MOSE_Mod_SST, only: SST_Eddy_Viscosity
    use MOSE_Mod_Spalart, only: Spalart_Eddy_Viscosity
    use MOSE_Mod_Wilcox2006, only: Wilcox2006_Eddy_Viscosity
    use MOSE_Lib_RANS, only: eddy_dev
#ifdef MOSE_MCDIFF_KERNEL
    use FLINT_Lib_ThermoTransport_dev, only: f_ss_dev, co_DS_expr_dev
#else
    use FLINT_Lib_ThermoTransport_dev, only: f_ss_dev
#endif
    integer,  intent(in) :: b, n(3)
    real(R8), intent(in) :: cfl, vnn, rampa_fac
    real(R8), intent(inout) :: dtglobal
    integer,  intent(in) :: nur
    real(R8), intent(in) :: urs(nur)   ! Lib_Preconditioning ur_f: Ur of the cells (prec), 1 element otherwise
    logical,  intent(in) :: prec       ! integration-variables = prec
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
    call field_block_index_params(b, boff, pi, pj)

    !$acc data copyin(m_f, dl_f, rho_c, Rg_c, mil_c) copyout(dtlocal) copy(dtglobal)
    !$acc parallel loop collapse(3) gang vector reduction(min:dtglobal) &
    !$acc&  present(P_flat, m_f, dl_f, rho_c, Rg_c, mil_c, dtlocal, dt_f, urs) &
#ifdef MOSE_MCDIFF_KERNEL
    !$acc&  private(pl, vel, versor, Gvel0, rkw, c0, s, d, T_i, Tint, Tc_, Tdf_, Dmax_, Dk)
#else
    !$acc&  private(pl, vel, versor, Gvel0, rkw, c0, s, d)
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
      dtc = dtc * rampa_fac
      dtlocal(i,j,k) = dtc
      if ( n_cache_f > 1 ) dt_f(c0) = dtc
      dtglobal = min( dtc, dtglobal )
    end do; end do; end do
    !$acc end data

  end subroutine DTW_NAME
