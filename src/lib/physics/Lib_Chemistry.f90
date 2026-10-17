module MOSE_Lib_Chemistry
  use OSLo

  implicit none
  logical :: reactive

! MOSE_ODE_RODAS: swap the device chemistry integrator radau5_dev -> rodas_dev
! (device RODAS4 Rosenbrock; ~half the per-thread matrix storage of RADAU5 -> the
! occupancy lever). Compile-time (mirrors the mechanism-binary pattern) so the
! unused integrator is NOT compiled into the cell kernel (keeps its register/stack
! footprint at the chosen method's). Same CT call signature for both.
# if defined (MOSE_ODE_RODAS)
#   define ODE_DEV rodas_dev
# else
#   define ODE_DEV radau5_dev
# endif

contains



  !----------------------------------------------------------------------------
  ! GPU cell loop: same per-cell algebra as the host loop above, but the
  ! integration is the devirtualized device RADAU5 (radau5_dev). The whole
  ! P block is moved H2D/D2H per call (correctness-first; residency later).
  !----------------------------------------------------------------------------

  !> Chemistry on the var-LAST flat SoA w. Per-cell twin of
  !> Chemistry_Newstate_acc; byte-identical algebra (cons-species + T from w -> radau5_dev
  !> -> back), only the access changes (P(s,i,j,k) -> w(c0,s) gathered into Z). Lets the
  !> reacted state live in w so the next Strang substage's halo on P_flat / flux read it without
  !> a P round-trip. c0 = cidx computed inline from host scalars (field_block_index_params).
  subroutine Chemistry_Newstate_blk ( b, dt, n, hr )
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic, only: f_Rtot
    use FLINT_Lib_Radau5_dev, only: radau5_dev, rodas_dev, rtolt_dev, atolt_dev, fnewt_dev, NSMX
    use MOSE_Field_m, only: P_flat, field_block_index_params, hr_f, n_cache_f   ! hr_f: the flat twin of hr_c
    implicit none
    integer, intent(in) :: b, n(3)
    real(kind=8), dimension(n(1),n(2),n(3)), intent(in) :: dt
    real(kind=8), dimension( 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc ), intent(inout) :: hr
    ! Local
    integer :: lz, nsc_l, np_l, i, j, k, s, err, nf, nj, nstp, na, nr, c0, boff, pi, pj
    real(8) :: Z(NSMX), rho, Tb
    real(8) :: fnw

    lz = nsc + 1; nsc_l = nsc; np_l = np
    ! MOSE_CHEM_CT: radau5_dev is specialized to a compile-time ODE size.
    if (lz /= MOSE_CHEM_NZ) error stop "MOSE_CHEM_CT: nsc+1 /= MOSE_CHEM_NZ (rebuild with matching MOSE_NSC)"
    fnw = fnewt_dev
    call field_block_index_params(b, boff, pi, pj)

    !$acc data copyin(dt) copy(hr)
    !$acc parallel loop collapse(3) gang vector present(P_flat, dt, hr, hr_f) private(Z, c0, Tb)
    do k = 1, n(3); do j = 1, n(2); do i = 1, n(1)
      c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      rho = 0d0
      do s = 1, nsc_l
        Z(s) = P_flat(c0,s); rho = rho + Z(s)
      end do
      Z(lz) = P_flat(c0,np_l) / (rho*f_Rtot(Z(1:nsc_l)))
      Tb = Z(lz)

      call ODE_DEV ( lz, 0d0, Z(1:lz), dt(i,j,k), rtolt_dev(1:lz), atolt_dev(1:lz), &
                        fnw, err, nf, nj, nstp, na, nr )

      do s = 1, nsc_l
        if (Z(s)<=1d-20) Z(s) = 1d-20
        P_flat(c0,s) = Z(s)
      end do
      P_flat(c0,np_l) = Z(lz)*rho*f_Rtot(Z(1:nsc_l))
      hr(i,j,k) = abs(Z(lz) - Tb) / Tb
      if ( n_cache_f > 1 ) hr_f(c0) = hr(i,j,k)   ! the flat twin (the same double)
    end do; end do; end do
    !$acc end data
  end subroutine Chemistry_Newstate_blk

end module MOSE_Lib_Chemistry
