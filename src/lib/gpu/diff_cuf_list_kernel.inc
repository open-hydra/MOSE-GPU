! ---------------------------------------------------------------------------------------------
! The diffusive pass 1 over the face list of one direction of a multirate level (Mod_Fluxes Fluxes_List,
! Lib_Multirate_Base mr_fl_type), instantiated twice in Mod_Diffusive.f90 as diff_cuf_kernels.inc (MOSE_MCDIFF_KERNEL
! undefined / defined; macro DCUF_KL = the kernel name). The face of list position e: left cell c1 = c(e), right cell
! c1 + sa(e), the transverse neighbours along the axes ib (stride sb) and ic (stride sc) -- the ten cells of
! diff_cuf_i/j/k --, the cell metrics, face metric and wall distances gathered by position, the face thermo tf(e, :)
! of Diffusive_Face_Thermo_Fill_List. Axes of direction d as diff_cuf_i/j/k: (1,2,3), (2,1,3), (3,1,2).
! ---------------------------------------------------------------------------------------------
  attributes(global) launch_bounds(128) subroutine DCUF_KL ( &
       w, dfh, ne, d, c, sa, sb, sc, fmg, mmg, yng, rho_c, T_c, &
       tf, ncell, &
       RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
    integer, value :: ncell, ne, d
    integer, value :: stress_dev, eddy_dev
    logical, value :: RANS
    real(R8), value :: sc_dev, sct_dev, prt_dev
    real(R8) :: w(ncell, nprim), dfh(ne, nprim)
    integer  :: c(ne), sa(ne), sb(ne), sc(ne)
    real(R8) :: fmg(4, ne), mmg(3, 3, 2, ne), yng(2, ne)
#ifdef MOSE_MCDIFF_KERNEL
    real(R8) :: tf(ne, 5+nsc)
#else
    real(R8) :: tf(ne, 5)
#endif
    real(R8) :: rho_c(ncell), T_c(ncell)
    integer :: e, s, ia, ib, ic
    real(R8) :: p1(NPRIM_DEV),p2(NPRIM_DEV),p3(NPRIM_DEV),p4(NPRIM_DEV),p5(NPRIM_DEV)
    real(R8) :: p6(NPRIM_DEV),p7(NPRIM_DEV),p8(NPRIM_DEV),p9(NPRIM_DEV),p10(NPRIM_DEV)
    real(R8) :: rho_n10(10), T_n10(10), M1l(3,3), M2l(3,3), nrm(3), res0(NPRIM_DEV)
    real(R8) :: area, wd1, wd2, rho_f, T_f, mil_f, kl_f, cp_f
    integer  :: c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
    include 'diff_decls.inc'
    e = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (e > ne) return
    if (d == 1) then
      ia=1; ib=2; ic=3
    else if (d == 2) then
      ia=2; ib=1; ic=3
    else
      ia=3; ib=1; ic=2
    end if
    c1 = c(e)
    c2 = c1 + sa(e)
    c3 = c1 - sb(e); c4 = c1 + sb(e); c5 = c2 - sb(e); c6 = c2 + sb(e)
    c7 = c1 - sc(e); c8 = c1 + sc(e); c9 = c2 - sc(e); c10 = c2 + sc(e)
    do s = 1, nprim
      p1(s)=w(c1,s); p2(s)=w(c2,s)   ! p3..p10 streamed load-at-use in diff_physics.inc (Grad ib/ic)
    end do
    rho_n10(1)=rho_c(c1);   rho_n10(2)=rho_c(c2)
    rho_n10(3)=rho_c(c3); rho_n10(4)=rho_c(c4)
    rho_n10(5)=rho_c(c5); rho_n10(6)=rho_c(c6)
    rho_n10(7)=rho_c(c7); rho_n10(8)=rho_c(c8)
    rho_n10(9)=rho_c(c9); rho_n10(10)=rho_c(c10)
    T_n10(1)=T_c(c1);   T_n10(2)=T_c(c2)
    T_n10(3)=T_c(c3); T_n10(4)=T_c(c4)
    T_n10(5)=T_c(c5); T_n10(6)=T_c(c6)
    T_n10(7)=T_c(c7); T_n10(8)=T_c(c8)
    T_n10(9)=T_c(c9); T_n10(10)=T_c(c10)
    M1l = mmg(:,:,1,e); M2l = mmg(:,:,2,e)
    nrm(1)=fmg(1,e); nrm(2)=fmg(2,e); nrm(3)=fmg(3,e)
    area = fmg(4,e); wd1 = yng(1,e); wd2 = yng(2,e)
    rho_f=tf(e,1); T_f=tf(e,2); mil_f=tf(e,3); kl_f=tf(e,4); cp_f=tf(e,5)
#define DCUF_TF tf
#define DCUF_ROW e
#include "diff_physics.inc"
#undef DCUF_TF
#undef DCUF_ROW
    do s = 1, nprim
      dfh(e, s) = res0(s)
    end do
  end subroutine DCUF_KL
