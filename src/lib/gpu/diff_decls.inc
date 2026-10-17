    ! ===== inlined-diffusive locals (diff_decls.inc) =====
    real(R8) :: irho1, irho2, irt1, irt2, irt3, irt4, waldis
    real(R8) :: Grad(NPRIM_DEV,3), Primf(NPRIM_DEV), Flx(NPRIM_DEV), Dm(NPRIM_DEV)
    real(R8) :: Mav(3,3), gg1, gg2, gg3
    real(R8) :: drho, dT, dcp, dmil, dkl, Tdf
    real(R8) :: Gvel(3,3), mie, dkap, dome, arg2, dF2, Sij(3,3), Smag
    real(R8) :: kappaf, dmi, dlam, Diver, lamDiv, rhokd
    real(R8) :: Txx, Tyy, Tzz, Txy, Txz, Tyz, Strs(3)
    real(R8) :: DiffH, DmGY, hh1, hh2, hh, sflx
    real(R8) :: gkdot, gwdot, gkw, CDw, zeta, dF1, sigk, sigw
    real(R8) :: chi_, gndot   ! SA (eddy_dev = 2)
    real(R8) :: wdiag         ! Wilcox2006 (eddy_dev = 3)
    integer  :: iv, dd, Ti
