! OVERFLOW2.1 Riemann solvers
! Reference:
!   Tramel, R., Nichols, R., & Buning, P. (2009). Addition of improved shock-capturing schemes to OVERFLOW 2.1. In 19th AIAA Computational Fluid Dynamics (p. 3988).
! Note:
! - HLLC+  | The minimum allowed beta is set to 0.0 instead of 0.4 to allow pure HLLE flows
! - HLLE++ | The Roe lambda_1 is set to un_roe instead of max(un_roe,a_roe). 
!            The adopted forumlation preserves pure Eulerian shear flow, contrary to the reference one.

module MOSE_Lib_Riemann_HLL
  implicit none
  private
  public :: riemann_HLLE, riemann_HLLC

  real(kind=8), public :: min_beta=0.0d0

contains

   ! HLLEM
   ! Einfeldt, B., Munz, C.C., Roe, P.L., and Sjogreen, B., "On Godunov-type methods near low densities", J. Comput. Phys., 92 (1991), pp. 273–295.

   subroutine riemann_HLLC(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,dummy,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
      !$acc routine seq
      use MOSE_Global_m, only: nsc
      implicit none

      ! Inputs
      real(8), intent(in) :: dl(nsc), dr(nsc)
      real(8), intent(in) :: ul, vl, wl, ur, vr, wr
      real(8), intent(in) :: pl, pr, al, ar
      real(8), intent(in) :: dltot, drtot
      real(8), intent(in) :: nx, ny, nz
      real(8), intent(in) :: dummy
      ! Outputs
      real(8), intent(out) :: F_r, F_u, F_v, F_w, F_E
      ! Roe + flow variables
      real(8) :: rho_roe, u_roe, v_roe, w_roe, a_roe, un_roe, h0_roe
      real(8) :: h0l, h0r, unl, unr
      ! HLLC variables
      real(8) :: S1, S4, Sstar, pstar
      real(8) :: inv_denom
      real(8) :: U1S, U2S, U3S, U4S, U5S
      real(8) :: e0l, e0r

      unl = ul*nx + vl*ny + wl*nz
      unr = ur*nx + vr*ny + wr*nz

      !------------------------------------------------
      ! Roe averages
      call roe_averages(nsc, dl, dr, dltot, drtot, &
                     ul, vl, wl, ur, vr, wr, pl, pr, h0l, h0r, &
                     rho_roe, u_roe, v_roe, w_roe, a_roe, h0_roe)

      un_roe = u_roe*nx + v_roe*ny + w_roe*nz

      ! Total energy (conservative variable) from total enthalpy
      e0l = h0l - pl / dltot
      e0r = h0r - pr / drtot

      !------------------------------------------------
      ! Batten signal speeds

      S1 = min(0.d0, unl - al, un_roe - a_roe)
      S4 = max(0.d0, unr + ar, un_roe + a_roe)

      inv_denom = 1.d0 / ( dltot*(S1-unl) - drtot*(S4-unr) )

      Sstar = ( pr - pl + dltot*unl*(S1-unl) - drtot*unr*(S4-unr) ) * inv_denom
      pstar = dltot*(unl-S1)*(unl-Sstar) + pl

      !------------------------------------------------
      ! HLLC fluxes
      if (S1 >= 0.d0) then

         F_r = dltot * unl
         F_u = F_r*ul + pl*nx
         F_v = F_r*vl + pl*ny
         F_w = F_r*wl + pl*nz
         F_E = F_r*h0l

      elseif (Sstar >= 0.d0) then

         F_r = dltot * unl
         F_u = F_r*ul + pl*nx
         F_v = F_r*vl + pl*ny
         F_w = F_r*wl + pl*nz
         F_E = F_r*h0l

         inv_denom = 1.d0 / (S1 - Sstar)

         U1S = dltot*(S1-unl)*inv_denom
         U2S = ((S1-unl)*dltot*ul + (pstar-pl)*nx)*inv_denom
         U3S = ((S1-unl)*dltot*vl + (pstar-pl)*ny)*inv_denom
         U4S = ((S1-unl)*dltot*wl + (pstar-pl)*nz)*inv_denom
         U5S = ((S1-unl)*dltot*e0l - pl*unl + pstar*Sstar)*inv_denom

         F_r = F_r + S1*(U1S - dltot)
         F_u = F_u + S1*(U2S - dltot*ul)
         F_v = F_v + S1*(U3S - dltot*vl)
         F_w = F_w + S1*(U4S - dltot*wl)
         F_E = F_E + S1*(U5S - dltot*e0l)

      elseif (S4 > 0.d0) then

         F_r = drtot * unr
         F_u = F_r*ur + pr*nx
         F_v = F_r*vr + pr*ny
         F_w = F_r*wr + pr*nz
         F_E = F_r*h0r

         inv_denom = 1.d0 / (S4 - Sstar)

         U1S = drtot*(S4-unr)*inv_denom
         U2S = ((S4-unr)*drtot*ur + (pstar-pr)*nx)*inv_denom
         U3S = ((S4-unr)*drtot*vr + (pstar-pr)*ny)*inv_denom
         U4S = ((S4-unr)*drtot*wr + (pstar-pr)*nz)*inv_denom
         U5S = ((S4-unr)*drtot*e0r - pr*unr + pstar*Sstar)*inv_denom

         F_r = F_r + S4*(U1S - drtot)
         F_u = F_u + S4*(U2S - drtot*ur)
         F_v = F_v + S4*(U3S - drtot*vr)
         F_w = F_w + S4*(U4S - drtot*wr)
         F_E = F_E + S4*(U5S - drtot*e0r)

      else

         F_r = drtot * unr
         F_u = F_r*ur + pr*nx
         F_v = F_r*vr + pr*ny
         F_w = F_r*wr + pr*nz
         F_E = F_r*h0r

      endif

   end subroutine riemann_HLLC


   subroutine riemann_HLLE(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,dummy,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
      !$acc routine seq
      use MOSE_Global_m, only: nsc
      implicit none

      ! Inputs
      real(8), intent(in) :: dl(nsc), dr(nsc)
      real(8), intent(in) :: ul, vl, wl, ur, vr, wr
      real(8), intent(in) :: pl, pr, al, ar
      real(8), intent(in) :: dltot, drtot
      real(8), intent(in) :: nx, ny, nz
      real(8), intent(in) :: dummy
      ! Outputs
      real(8), intent(out) :: F_r, F_u, F_v, F_w, F_E
      ! Roe data
      real(8) :: rho_roe, u_roe, v_roe, w_roe, a_roe, un_roe, h0_roe
      real(8) :: h0l, h0r, unl, unr
      ! HLLE variables
      real(8) :: S1, S4, inv_denom
      real(8) :: FrL, FrR

      unl = ul*nx + vl*ny + wl*nz 
      unr = ur*nx + vr*ny + wr*nz

      !------------------------------------------------
      ! Roe averages
      call roe_averages(nsc, dl, dr, dltot, drtot, &
                     ul, vl, wl, ur, vr, wr, pl, pr, h0l, h0r, &
                     rho_roe, u_roe, v_roe, w_roe, a_roe, h0_roe)

      un_roe = u_roe*nx + v_roe*ny + w_roe*nz

      !------------------------------------------------
      ! Batten signal speeds

      S1 = min(0.d0, unl - al, un_roe - a_roe)
      S4 = max(0.d0, unr + ar, un_roe + a_roe)

      !------------------------------------------------
      ! HLL fluxes

      inv_denom = 1.d0 / (S4 - S1)

      FrL = dltot*unl
      FrR = drtot*unr

      F_r = ( S4*FrL - S1*FrR + S1*S4*(drtot-dltot) ) * inv_denom

      F_u = ( S4*(FrL*ul + pl*nx) - S1*(FrR*ur + pr*nx) &
            + S1*S4*(drtot*ur - dltot*ul) ) * inv_denom

      F_v = ( S4*(FrL*vl + pl*ny) - S1*(FrR*vr + pr*ny) &
            + S1*S4*(drtot*vr - dltot*vl) ) * inv_denom

      F_w = ( S4*(FrL*wl + pl*nz) - S1*(FrR*wr + pr*nz) &
            + S1*S4*(drtot*wr - dltot*wl) ) * inv_denom

      F_E = ( S4*(FrL*h0l) - S1*(FrR*h0r) &
            + S1*S4*( drtot*(h0r - pr/drtot) &
                     -dltot*(h0l - pl/dltot) ) ) * inv_denom

   end subroutine riemann_HLLE


  ! HLLE++

  ! HLLEM - HLLE 
  ! Hybrid solver working with a mulidimensional shock detector

  ! HLLC - HLLE 
  ! Hybrid solver working with a mulidimensional shock detector


pure subroutine roe_averages(nsc, dl, dr, dltot, drtot, ul, vl, wl, ur, vr, wr, pl, pr, h0l, h0r, rho_roe, u_roe, v_roe, w_roe, a_roe, h0_roe)
      !$acc routine seq
  use FLINT_Lib_Thermodynamic
  implicit none
  ! Input
  integer, intent(in) :: nsc
  real(8), intent(in) :: dl(nsc), dr(nsc)
  real(8), intent(in) :: dltot, drtot
  real(8), intent(in) :: ul, vl, wl, ur, vr, wr
  real(8), intent(in) :: pl, pr
  ! Output
  real(8), intent(out) :: h0l, h0r
  real(8), intent(out) :: rho_roe, u_roe, v_roe, w_roe, a_roe, h0_roe
  ! Local 
  integer :: s, Til, Tir
  real(8) :: Rl, Rr, Tl, Tr, dTl, dTr
  real(8) :: inv_dltot, inv_drtot
  real(8) :: srL, srR, inv_sr
  real(8) :: cv_roe, sum_ei
  real(8) :: R_roe, T_roe
  real(8) :: gam_roe, vel2
  real(8) :: hl, hr, el, er, invW, d_roe

  !------------------------------------------------
  ! Precompute thermodynamics
  Rl = f_Rtot(dl)
  Rr = f_Rtot(dr)

  inv_dltot = 1.d0 / dltot
  inv_drtot = 1.d0 / drtot

  Tl = pl * inv_dltot / Rl
  Tr = pr * inv_drtot / Rr

  Til = int(Tl)
  Tir = int(Tr)
  dTl = Tl - Til
  dTr = Tr - Tir

  srL = sqrt(dltot)
  srR = sqrt(drtot)
  inv_sr = 1.d0 / (srL + srR)

  !------------------------------------------------
  ! Roe flow averages
  rho_roe = srL * srR

  u_roe = (srR*ur + srL*ul) * inv_sr
  v_roe = (srR*vr + srL*vl) * inv_sr
  w_roe = (srR*wr + srL*wl) * inv_sr

  R_roe = (srR*Rr + srL*Rl) * inv_sr
  T_roe = (srR*Tr + srL*Tl) * inv_sr

  !------------------------------------------------
  ! Species loop
  h0l = 0.d0
  h0r = 0.d0
  cv_roe = 0.d0
  sum_ei = 0.d0

  do s = 1, nsc

    invW = Runiv / Wm_tab(s)

    hl = h_tab(Til, s) + (h_tab(Til+1, s) - h_tab(Til, s)) * dTl
    hr = h_tab(Tir, s) + (h_tab(Tir+1, s) - h_tab(Tir, s)) * dTr

    el = hl - invW * Tl
    er = hr - invW * Tr

    h0l = h0l + hl * dl(s) * inv_dltot
    h0r = h0r + hr * dr(s) * inv_drtot

    d_roe = (srR*dr(s)*inv_drtot + srL*dl(s)*inv_dltot) * inv_sr

    cv_roe = cv_roe + d_roe * ( &
        0.5d0 * ( &
          cp_tab(Til, s) + (cp_tab(Til+1, s)-cp_tab(Til, s))*dTl + &
          cp_tab(Tir, s) + (cp_tab(Tir+1, s)-cp_tab(Tir, s))*dTr ) &
        - invW )

    sum_ei = sum_ei + d_roe * (srR*er + srL*el) * inv_sr

  end do

  !------------------------------------------------
  ! Final scalars
  h0l = h0l + 0.5d0*(ul*ul + vl*vl + wl*wl)
  h0r = h0r + 0.5d0*(ur*ur + vr*vr + wr*wr)

  h0_roe = (srR*h0r + srL*h0l) * inv_sr

  gam_roe = 1.d0 + R_roe / cv_roe

  vel2 = u_roe*u_roe + v_roe*v_roe + w_roe*w_roe

  a_roe = sqrt( (gam_roe - 1.d0) * &
        ( h0_roe - 0.5d0*vel2 + cv_roe*T_roe - sum_ei ) )

end subroutine roe_averages


! Subroutine for computing the conservative fluxes from primitive variables.
pure subroutine fluxes(p,r,u,v,w,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
  use MOSE_Global_m, only: nsc
  use FLINT_Lib_Thermodynamic, only: H0
  implicit none
  real(8), intent(IN)::  p           ! Pressure
  real(8), intent(IN)::  r(nsc)      ! Density
  real(8), intent(IN)::  u,v,w       ! Velocity
  real(8), intent(IN)::  nx,ny,nz    ! Normals
  real(8), intent(OUT):: F_r         ! Flux of mass conservation
  real(8), intent(OUT):: F_u,F_v,F_w ! Flux of momentum conservation
  real(8), intent(OUT):: F_E         ! Flux of energy conservation
  F_r = sum(r)*(u*nx+v*ny+w*nz)
  F_u = F_r*u + p*nx
  F_v = F_r*v + p*ny
  F_w = F_r*w + p*nz
  F_E = F_r*H0(p,r,sqrt(u*u+v*v+w*w))
end subroutine fluxes




end module MOSE_Lib_Riemann_HLL