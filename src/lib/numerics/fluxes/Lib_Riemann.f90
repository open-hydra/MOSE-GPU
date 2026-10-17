module MOSE_Lib_Riemann

  implicit none
  private
  public :: Assign_RS
  public :: Riemann   ! the device Riemann dispatch (rs_dev = 1 HLLC Batten, 2 HLLE, 3 SLAU, 4-7 below), acc routine
  public :: riemann_open   ! MOSE_open AUSM+ (4), AUSM+M (5), HLLE++ (6), HLLC+Tramel (7), LLF (8), SLAU2 (9), HLLC+Chen (10),
                           ! exact (11), LMRoe (12), MiczekRoe (13): riemann_open_body.inc

  logical, public :: SD_riemann

  !> Devirtualized Riemann-scheme id for the GPU device build (procedure
  !> pointers can't run on the GPU). 1 = HLLC Batten (riemann_HLLC), the only
  !> scheme the device path supports; 0 = unsupported (device hard-stops).
  integer, public :: rs_dev = 0
  !$acc declare create(rs_dev)
  !> riemann-options-Mco of MOSE_open (the low-Mach cutoff of AUSM+M): [MOSE-scheme] riemann-options-Mco.
  real(kind=8), public :: mco_dev = 0.d0
  !$acc declare create(mco_dev)

contains

  subroutine Assign_RS(RS_word, Minf_, Mco_)
    use MOSE_Lib_Preconditioning, only: obj_prec
    implicit none
    character(len=*), intent(inout) :: RS_word
    real(kind=8), intent(in)        :: Minf_
    real(kind=8), intent(in), optional :: Mco_

    rs_dev = 0
    mco_dev = 0.d0
    if (present(Mco_)) mco_dev = Mco_

    if ( RS_word(len_trim(RS_word)-2:len_trim(RS_word))=='SD' ) then
      RS_word = trim(RS_word(1:len_trim(RS_word)-3))
      SD_riemann = .true.
    else
      SD_riemann = .false.
    endif

    select case (RS_word)
    case ('HLLC Batten','HLLC')
      RS_word = 'HLLC Batten'
      if (SD_riemann) error stop 'riemann-solver HLLC with the shock detector (HLLCSD) is not on the GPU'
      rs_dev = 1   ! riemann_HLLC / the inline HLLC of the kernels
    case ('HLLE')
      RS_word = 'HLLE'
      rs_dev = 2   ! riemann_HLLE / the inline HLLE of the kernels
    case ('SLAU')
      RS_word = 'SLAU'
      rs_dev = 3   ! riemann_SLAU / the inline SLAU of the kernels
    case ('AUSM+')
      RS_word = 'AUSM+'
      rs_dev = 4   ! MOSE_open riemann_AUSMp (riemann_open_body.inc)
    case ('AUSM+M')
      RS_word = 'AUSM+M'
      rs_dev = 5   ! MOSE_open riemann_AUSMp_M (Chen et al. 2020): needs riemann-options-Mco, sensor Chen by default
      if (mco_dev == 0.d0) error stop '[ERROR] AUSM+M: riemann-options-Mco must be given in [MOSE-scheme] (MOSE_open: Mco)'
    case ('HLLE++')
      RS_word = 'HLLE++'
      rs_dev = 6   ! MOSE_open riemann_HLLEpp: the sensor (Tramel by default, Assign_SD) blends Roe and HLLE waves
    case ('HLLC+Tramel')
      RS_word = 'HLLC+Tramel'
      rs_dev = 7   ! MOSE_open riemann_HLLCp_Tramel: HLLC and HLLE blended by the sensor (Tramel by default)
    case ('LLF', 'Rusanov')
      RS_word = 'LLF'
      rs_dev = 8   ! MOSE_open riemann_LLF
    case ('SLAU2')
      RS_word = 'SLAU2'
      rs_dev = 9   ! MOSE_open riemann_SLAU2
    case ('HLLC+Chen')
      RS_word = 'HLLC+Chen'
      rs_dev = 10  ! MOSE_open riemann_HLLCp_Chen: sensor Chen by default, Mco the floor of the local Mach
      if (mco_dev == 0.d0) error stop '[ERROR] HLLC+Chen: riemann-options-Mco must be given in [MOSE-scheme] (MOSE_open: Mco)'
    case ('exact')
      RS_word = 'exact'
      rs_dev = 11  ! MOSE_open riemann_exact (Godunov)
    case ('LMRoe')
      RS_word = 'LMRoe'
      rs_dev = 12  ! MOSE_open Riemann_LMRoe (low-Mach Roe), Mco the floor of the low-Mach factor
      if (mco_dev == 0.d0) error stop '[ERROR] LMRoe: riemann-options-Mco must be given in [MOSE-scheme] (MOSE_open: Mco)'
    case ('MiczekRoe')
      RS_word = 'MiczekRoe'
      rs_dev = 13  ! MOSE_open Riemann_MiczekRoe, Mco the cut-off Mach number (floor 1e-6)
    case ('HLLC-PC')
      RS_word = 'HLLC-PC'
      rs_dev = 14  ! MOSE_open Riemann_HLLCprec: HLLC Batten with the preconditioned signal velocities
      if ( .not. obj_prec%enabled ) error stop &
        '[ERROR] Preconditioned HLL solver selected. integration-variables must be set to prec (as MOSE_open)'
    case default
      write(*,'(3A)') ' [ERROR] riemann-solver ', trim(RS_word), ' is not on the GPU: HLLC, HLLE, SLAU, AUSM+, AUSM+M, HLLE++, &
                      &HLLC+Tramel, LLF, Rusanov, SLAU2, HLLC+Chen, exact, LMRoe, MiczekRoe, HLLC-PC (with prec)'
      error stop 'Assign_RS: riemann-solver not on the GPU'
    end select

    !> Push the ids to the device: the record BC kernels (Lib_BC, through Riemann) read rs_dev and mco_dev from the device
    !> copies. Both are created on the device by IO_Chemistry, which runs first
    !> (the pattern of lim_dev in Assign_Limiter and of the RANS ids in Set_RANS_Model).
    !$acc update device(rs_dev, mco_dev)
    !> The limiter id is set by Assign_Limiter (called before, IO_Ini): 0 = a limiter the device kernels do not have.
    block
      use MOSE_Lib_Limiters, only: lim_dev
      if (lim_dev == 0) error stop 'Assign_RS: the limiter is not on the GPU (the device kernels have Min-Mod, Van Leer, &
                                   &Van Albada, MC, Super-Bee and first order)'
    end block

  end subroutine Assign_RS


  !> The device Riemann solver of the case, for the OpenACC BC routines and the host side of the device build:
  !> rs_dev = 1 HLLC Batten (riemann_HLLC), 2 HLLE (riemann_HLLE), 3 SLAU (riemann_SLAU), 4..14 the MOSE_open solvers of
  !> riemann_open_body.inc (riemann_open). The CUF convective kernels run the same solvers (conv_physics.inc,
  !> riemann_open_d); the record BC kernels of Lib_BC call this routine.
  subroutine Riemann(rs, dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,dummy,urL_,urR_,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    !$acc routine seq
    use MOSE_Global_m, only: nsc
    use MOSE_Lib_Riemann_HLL, only: riemann_HLLC, riemann_HLLE
    use FLINT_Lib_Thermodynamic, only: Runiv
    implicit none
    integer, intent(in) :: rs
    real(8), intent(in) :: dl(nsc), dr(nsc)
    real(8), intent(in) :: ul, vl, wl, ur, vr, wr, pl, pr, al, ar, dltot, drtot, nx, ny, nz, dummy
    real(8), intent(in) :: urL_, urR_   ! MOSE_open url, urr: the Ur of the two sides (HLLC-PC only)
    real(8), intent(out) :: F_r, F_u, F_v, F_w, F_E
    select case (rs)
    case (2)
      call riemann_HLLE(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,dummy,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    case (3)
      call riemann_SLAU(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    case (4:14)  ! dummy = the face sensor beta (MOSE_open passes blk%beta of the cell, 1 at symmetries)
      call riemann_open ( rs, dummy, mco_dev, Runiv, nsc, &
                              dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,nx,ny,nz,urL_,urR_,F_r,F_u,F_v,F_w,F_E )
    case default
      call riemann_HLLC(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,dummy,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    end select
  end subroutine Riemann

  !> SLAU on the device (Lib_Riemann_SLAU riemann_SLAU_base + riemann_SLAU): the total enthalpies H0 = h(T) + V^2/2 of the
  !> two states from the FLINT tables (the linear-in-T interpolation of h_tab the device kernels use), T = p/(rho R).
  subroutine riemann_SLAU(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    !$acc routine seq
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: h_tab, Ri_tab
    implicit none
    real(8), intent(in) :: dl(nsc), dr(nsc)
    real(8), intent(in) :: ul, vl, wl, ur, vr, wr, pl, pr, al, ar, dltot, drtot, nx, ny, nz
    real(8), intent(out) :: F_r, F_u, F_v, F_w, F_E
    real(8) :: unl, unr, h0l, h0r, Rl, Rr, Tl, Tr, dTl, dTr
    real(8) :: frl_, frr_, ml_, mr_, mt_, chi_, vnm_, frho_, vnl_, vnr_, bl_, br_, af_, ph_
    integer :: s, Til, Tir
    unl = ul*nx + vl*ny + wl*nz
    unr = ur*nx + vr*ny + wr*nz
    Rl = 0.d0; Rr = 0.d0
    do s = 1, nsc
      Rl = Rl + Ri_tab(s)*dl(s)/dltot
      Rr = Rr + Ri_tab(s)*dr(s)/drtot
    end do
    Tl = pl/(dltot*Rl); Tr = pr/(drtot*Rr)
    Til = int(Tl); Tir = int(Tr); dTl = Tl - Til; dTr = Tr - Tir
    h0l = 0.d0; h0r = 0.d0
    do s = 1, nsc
      h0l = h0l + ( h_tab(Til,s) + (h_tab(Til+1,s) - h_tab(Til,s))*dTl ) * dl(s)/dltot
      h0r = h0r + ( h_tab(Tir,s) + (h_tab(Tir+1,s) - h_tab(Tir,s))*dTr ) * dr(s)/drtot
    end do
    h0l = h0l + 0.5d0*(ul*ul + vl*vl + wl*wl)
    h0r = h0r + 0.5d0*(ur*ur + vr*vr + wr*wr)
    ! SLAU (Lib_Riemann_SLAU riemann_SLAU_base + riemann_SLAU): h0l/h0r are the total enthalpies H0 of the host
    af_ = 0.5d0*(al + ar)
    ml_ = unl / af_; mr_ = unr / af_
    mt_ = min(1.d0, sqrt(0.5d0*((ul*ul + vl*vl + wl*wl) + (ur*ur + vr*vr + wr*wr))) / af_)
    chi_ = (1.d0 - mt_)**2
    vnm_ = (dltot*abs(unl) + drtot*abs(unr)) / (dltot + drtot)
    frho_ = -max(min(ml_,0.d0), -1.d0) * min(max(mr_,0.d0), 1.d0)
    vnl_ = (1.d0 - frho_)*vnm_ + frho_*abs(unl)
    vnr_ = (1.d0 - frho_)*vnm_ + frho_*abs(unr)
    F_r = 0.5d0 * ( dltot*(unl + vnl_) + drtot*(unr - vnr_) - (chi_/af_)*(pr - pl) )
    if (abs(ml_) < 1.d0) then
      bl_ = 0.25d0*(2.d0 - ml_)*(ml_ + 1.d0)**2
    else if (ml_ >= 0.d0) then
      bl_ = 1.d0
    else
      bl_ = 0.d0
    endif
    if (abs(mr_) < 1.d0) then
      br_ = 0.25d0*(2.d0 + mr_)*(mr_ - 1.d0)**2
    else if (mr_ >= 0.d0) then
      br_ = 0.d0
    else
      br_ = 1.d0
    endif
    ph_ = 0.5d0*(pl + pr) + 0.5d0*(bl_ - br_)*(pl - pr) + (1.d0 - chi_)*(bl_ + br_ - 1.d0)*0.5d0*(pl + pr)
    if (F_r >= 0.d0) then
      F_u = F_r*ul + nx*ph_; F_v = F_r*vl + ny*ph_; F_w = F_r*wl + nz*ph_; F_E = F_r*h0l
    else
      F_u = F_r*ur + nx*ph_; F_v = F_r*vr + ny*ph_; F_w = F_r*wr + nz*ph_; F_E = F_r*h0r
    endif
  end subroutine riemann_SLAU

  !> MOSE_open AUSM+ / AUSM+M / HLLE++ / HLLC+Tramel for the OpenACC side and the host: the shared body of
  !> riemann_open_body.inc (the CUF twin is Mod_Convective_W_CUF riemann_open_d).
  subroutine riemann_open ( rs, bet, mco, runiv_, ns_, &
                                dl, ul, vl, wl, pl, al, dltot, dr, ur, vr, wr, pr, ar, drtot, nx, ny, nz, urL_, urR_, &
                                F_r, F_u, F_v, F_w, F_E )
    !$acc routine seq
    use FLINT_Lib_Thermodynamic, only: h_tab, cp_tab, Ri_tab, Wm_tab, Tmin, Tmax
#define RO_HTAB(i,s) h_tab(i,s)
#define RO_CPTAB(i,s) cp_tab(i,s)
#define RO_RITAB(s) Ri_tab(s)
#define RO_WMTAB(s) Wm_tab(s)
#define RO_TLO Tmin
#define RO_THI Tmax
#define RO_MCO mco
#define RO_INV intent(in)
#include "riemann_open_body.inc"
#undef RO_HTAB
#undef RO_CPTAB
#undef RO_RITAB
#undef RO_WMTAB
#undef RO_TLO
#undef RO_THI
#undef RO_MCO
#undef RO_INV
  end subroutine riemann_open

end module MOSE_Lib_Riemann