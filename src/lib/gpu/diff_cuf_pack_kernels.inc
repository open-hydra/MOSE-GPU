! ---------------------------------------------------------------------------------------------
! The block-pack twins of the three CUF pass-1 kernels of the diffusive flux (MOSE_Pack_m): one launch per direction over
! the faces of every block of a batch. CUDA block e = cb0 + blockIdx%x belongs to local block p = cbp(e) and covers the
! faces cbt0(e)+1 .. cbt0(e)+128 of the block's own face range (the thread -> face decode of diff_cuf_*, so the same face
! is computed by the same arithmetic). Geometry from the flat pools (Mp, ynp by the flat cell index; fm*p by the block's
! chunk offset), face thermo from the batch buffer tfp and fluxes into dfhp, both at the block's face offset tofh(d,p).
! Instantiated twice in Mod_Diffusive.f90 as diff_cuf_kernels.inc: MOSE_MCDIFF_KERNEL undefined -> DCUF_PKI/PKJ/PKK
! = diff_cuf_*_pack (tfp of 5 columns); defined -> diff_cuf_*_pack_mc (5+nsc columns, the face D_k of the fills).
! ---------------------------------------------------------------------------------------------
  attributes(global) launch_bounds(128) subroutine DCUF_PKI ( &
       w, dfhp, fm1p, Mp, ynp, rho_c, T_c, tfp, ncell, nfp, nfm, ntfl, cbp, cbt0, cb0, &
       tn1, tn2, tn3, tboff, tpi, tpj, tofm, tofh, nloc, gc, &
       RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
    integer, value :: ncell, nfp, nfm, ntfl, cb0, nloc, gc
    integer, value :: stress_dev, eddy_dev
    logical, value :: RANS
    real(R8), value :: sc_dev, sct_dev, prt_dev
    real(R8) :: w(ncell, nprim), dfhp(nfp, nprim)
    real(R8) :: fm1p(4, nfm)
#ifdef MOSE_MCDIFF_KERNEL
    real(R8) :: tfp(ntfl, 5+nsc)
#else
    real(R8) :: tfp(ntfl, 5)
#endif
    real(R8) :: Mp(3, 3, ncell), ynp(ncell)
    real(R8) :: rho_c(ncell), T_c(ncell)
    integer  :: cbp(*), cbt0(*)
    integer  :: tn1(nloc), tn2(nloc), tn3(nloc), tboff(nloc), tpi(nloc), tpj(nloc), tofm(nloc), tofh(3, nloc)
    integer :: i, j, k, s, f, ia, ib, ic, p, e, fi, n1, n2, n3, boff, pi, pj
    real(R8) :: p1(NPRIM_DEV),p2(NPRIM_DEV),p3(NPRIM_DEV),p4(NPRIM_DEV),p5(NPRIM_DEV)
    real(R8) :: p6(NPRIM_DEV),p7(NPRIM_DEV),p8(NPRIM_DEV),p9(NPRIM_DEV),p10(NPRIM_DEV)
    real(R8) :: rho_n10(10), T_n10(10), M1l(3,3), M2l(3,3), nrm(3), res0(NPRIM_DEV)
    real(R8) :: area, wd1, wd2, rho_f, T_f, mil_f, kl_f, cp_f
    integer  :: c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
    integer :: t_, ni_, nj_, rem_
    include 'diff_decls.inc'
    e = cb0 + blockIdx%x
    p = cbp(e)
    n1 = tn1(p); n2 = tn2(p); n3 = tn3(p); boff = tboff(p); pi = tpi(p); pj = tpj(p)
    t_  = cbt0(e) + threadIdx%x
    ni_ = n1 - 1
    nj_ = n2
    if (t_ > ni_*nj_*n3) return
    i    = mod(t_-1, ni_) + 1
    rem_ = (t_-1) / ni_
    j    = mod(rem_, nj_) + 1
    k    = rem_ / nj_ + 1
    ia=1; ib=2; ic=3
    c1 =boff+(i  +gc-1)+(j  +gc-1)*pi+(k  +gc-1)*pi*pj+1
    c2 =boff+(i+1+gc-1)+(j  +gc-1)*pi+(k  +gc-1)*pi*pj+1
    c3 =boff+(i  +gc-1)+(j-1+gc-1)*pi+(k  +gc-1)*pi*pj+1
    c4 =boff+(i  +gc-1)+(j+1+gc-1)*pi+(k  +gc-1)*pi*pj+1
    c5 =boff+(i+1+gc-1)+(j-1+gc-1)*pi+(k  +gc-1)*pi*pj+1
    c6 =boff+(i+1+gc-1)+(j+1+gc-1)*pi+(k  +gc-1)*pi*pj+1
    c7 =boff+(i  +gc-1)+(j  +gc-1)*pi+(k-1+gc-1)*pi*pj+1
    c8 =boff+(i  +gc-1)+(j  +gc-1)*pi+(k+1+gc-1)*pi*pj+1
    c9 =boff+(i+1+gc-1)+(j  +gc-1)*pi+(k-1+gc-1)*pi*pj+1
    c10=boff+(i+1+gc-1)+(j  +gc-1)*pi+(k+1+gc-1)*pi*pj+1
    do s = 1, nprim
      p1(s)=w(c1,s); p2(s)=w(c2,s)
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
    M1l = Mp(:,:,c1); M2l = Mp(:,:,c2)
    fi = tofm(p) + i + (n1+1)*((j-1) + n2*(k-1)) + 1
    nrm(1)=fm1p(1,fi); nrm(2)=fm1p(2,fi); nrm(3)=fm1p(3,fi)
    area = fm1p(4,fi); wd1 = ynp(c1); wd2 = ynp(c2)
    f = tofh(1,p) + (i-1) + (j-1)*(n1-1) + (k-1)*(n1-1)*n2 + 1
    rho_f=tfp(f,1); T_f=tfp(f,2); mil_f=tfp(f,3); kl_f=tfp(f,4); cp_f=tfp(f,5)
#define DCUF_TF tfp
#define DCUF_ROW f
#include "diff_physics.inc"
#undef DCUF_TF
#undef DCUF_ROW
    do s = 1, nprim
      dfhp(f, s) = res0(s)
    end do
  end subroutine DCUF_PKI

  attributes(global) launch_bounds(128) subroutine DCUF_PKJ ( &
       w, dfhp, fm2p, Mp, ynp, rho_c, T_c, tfp, ncell, nfp, nfm, ntfl, cbp, cbt0, cb0, &
       tn1, tn2, tn3, tboff, tpi, tpj, tofm, tofh, nloc, gc, &
       RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
    integer, value :: ncell, nfp, nfm, ntfl, cb0, nloc, gc
    integer, value :: stress_dev, eddy_dev
    logical, value :: RANS
    real(R8), value :: sc_dev, sct_dev, prt_dev
    real(R8) :: w(ncell, nprim), dfhp(nfp, nprim)
    real(R8) :: fm2p(4, nfm)
#ifdef MOSE_MCDIFF_KERNEL
    real(R8) :: tfp(ntfl, 5+nsc)
#else
    real(R8) :: tfp(ntfl, 5)
#endif
    real(R8) :: Mp(3, 3, ncell), ynp(ncell)
    real(R8) :: rho_c(ncell), T_c(ncell)
    integer  :: cbp(*), cbt0(*)
    integer  :: tn1(nloc), tn2(nloc), tn3(nloc), tboff(nloc), tpi(nloc), tpj(nloc), tofm(nloc), tofh(3, nloc)
    integer :: i, j, k, s, f, ia, ib, ic, p, e, fi, n1, n2, n3, boff, pi, pj
    real(R8) :: p1(NPRIM_DEV),p2(NPRIM_DEV),p3(NPRIM_DEV),p4(NPRIM_DEV),p5(NPRIM_DEV)
    real(R8) :: p6(NPRIM_DEV),p7(NPRIM_DEV),p8(NPRIM_DEV),p9(NPRIM_DEV),p10(NPRIM_DEV)
    real(R8) :: rho_n10(10), T_n10(10), M1l(3,3), M2l(3,3), nrm(3), res0(NPRIM_DEV)
    real(R8) :: area, wd1, wd2, rho_f, T_f, mil_f, kl_f, cp_f
    integer  :: c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
    integer :: t_, ni_, nj_, rem_
    include 'diff_decls.inc'
    e = cb0 + blockIdx%x
    p = cbp(e)
    n1 = tn1(p); n2 = tn2(p); n3 = tn3(p); boff = tboff(p); pi = tpi(p); pj = tpj(p)
    t_  = cbt0(e) + threadIdx%x
    ni_ = n1
    nj_ = n2 - 1
    if (t_ > ni_*nj_*n3) return
    i    = mod(t_-1, ni_) + 1
    rem_ = (t_-1) / ni_
    j    = mod(rem_, nj_) + 1
    k    = rem_ / nj_ + 1
    ia=2; ib=1; ic=3
    c1 =boff+(i  +gc-1)+(j  +gc-1)*pi+(k  +gc-1)*pi*pj+1
    c2 =boff+(i  +gc-1)+(j+1+gc-1)*pi+(k  +gc-1)*pi*pj+1
    c3 =boff+(i-1+gc-1)+(j  +gc-1)*pi+(k  +gc-1)*pi*pj+1
    c4 =boff+(i+1+gc-1)+(j  +gc-1)*pi+(k  +gc-1)*pi*pj+1
    c5 =boff+(i-1+gc-1)+(j+1+gc-1)*pi+(k  +gc-1)*pi*pj+1
    c6 =boff+(i+1+gc-1)+(j+1+gc-1)*pi+(k  +gc-1)*pi*pj+1
    c7 =boff+(i  +gc-1)+(j  +gc-1)*pi+(k-1+gc-1)*pi*pj+1
    c8 =boff+(i  +gc-1)+(j  +gc-1)*pi+(k+1+gc-1)*pi*pj+1
    c9 =boff+(i  +gc-1)+(j+1+gc-1)*pi+(k-1+gc-1)*pi*pj+1
    c10=boff+(i  +gc-1)+(j+1+gc-1)*pi+(k+1+gc-1)*pi*pj+1
    do s = 1, nprim
      p1(s)=w(c1,s); p2(s)=w(c2,s)
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
    M1l = Mp(:,:,c1); M2l = Mp(:,:,c2)
    fi = tofm(p) + (i-1) + n1*(j + (n2+1)*(k-1)) + 1
    nrm(1)=fm2p(1,fi); nrm(2)=fm2p(2,fi); nrm(3)=fm2p(3,fi)
    area = fm2p(4,fi); wd1 = ynp(c1); wd2 = ynp(c2)
    f = tofh(2,p) + (i-1) + (j-1)*n1 + (k-1)*n1*(n2-1) + 1
    rho_f=tfp(f,1); T_f=tfp(f,2); mil_f=tfp(f,3); kl_f=tfp(f,4); cp_f=tfp(f,5)
#define DCUF_TF tfp
#define DCUF_ROW f
#include "diff_physics.inc"
#undef DCUF_TF
#undef DCUF_ROW
    do s = 1, nprim
      dfhp(f, s) = res0(s)
    end do
  end subroutine DCUF_PKJ

  attributes(global) launch_bounds(128) subroutine DCUF_PKK ( &
       w, dfhp, fm3p, Mp, ynp, rho_c, T_c, tfp, ncell, nfp, nfm, ntfl, cbp, cbt0, cb0, &
       tn1, tn2, tn3, tboff, tpi, tpj, tofm, tofh, nloc, gc, &
       RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
    integer, value :: ncell, nfp, nfm, ntfl, cb0, nloc, gc
    integer, value :: stress_dev, eddy_dev
    logical, value :: RANS
    real(R8), value :: sc_dev, sct_dev, prt_dev
    real(R8) :: w(ncell, nprim), dfhp(nfp, nprim)
    real(R8) :: fm3p(4, nfm)
#ifdef MOSE_MCDIFF_KERNEL
    real(R8) :: tfp(ntfl, 5+nsc)
#else
    real(R8) :: tfp(ntfl, 5)
#endif
    real(R8) :: Mp(3, 3, ncell), ynp(ncell)
    real(R8) :: rho_c(ncell), T_c(ncell)
    integer  :: cbp(*), cbt0(*)
    integer  :: tn1(nloc), tn2(nloc), tn3(nloc), tboff(nloc), tpi(nloc), tpj(nloc), tofm(nloc), tofh(3, nloc)
    integer :: i, j, k, s, f, ia, ib, ic, p, e, fi, n1, n2, n3, boff, pi, pj
    real(R8) :: p1(NPRIM_DEV),p2(NPRIM_DEV),p3(NPRIM_DEV),p4(NPRIM_DEV),p5(NPRIM_DEV)
    real(R8) :: p6(NPRIM_DEV),p7(NPRIM_DEV),p8(NPRIM_DEV),p9(NPRIM_DEV),p10(NPRIM_DEV)
    real(R8) :: rho_n10(10), T_n10(10), M1l(3,3), M2l(3,3), nrm(3), res0(NPRIM_DEV)
    real(R8) :: area, wd1, wd2, rho_f, T_f, mil_f, kl_f, cp_f
    integer  :: c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
    integer :: t_, ni_, nj_, rem_
    include 'diff_decls.inc'
    e = cb0 + blockIdx%x
    p = cbp(e)
    n1 = tn1(p); n2 = tn2(p); n3 = tn3(p); boff = tboff(p); pi = tpi(p); pj = tpj(p)
    t_  = cbt0(e) + threadIdx%x
    ni_ = n1
    nj_ = n2
    if (t_ > ni_*nj_*(n3-1)) return
    i    = mod(t_-1, ni_) + 1
    rem_ = (t_-1) / ni_
    j    = mod(rem_, nj_) + 1
    k    = rem_ / nj_ + 1
    ia=3; ib=1; ic=2
    c1 =boff+(i  +gc-1)+(j  +gc-1)*pi+(k  +gc-1)*pi*pj+1
    c2 =boff+(i  +gc-1)+(j  +gc-1)*pi+(k+1+gc-1)*pi*pj+1
    c3 =boff+(i-1+gc-1)+(j  +gc-1)*pi+(k  +gc-1)*pi*pj+1
    c4 =boff+(i+1+gc-1)+(j  +gc-1)*pi+(k  +gc-1)*pi*pj+1
    c5 =boff+(i-1+gc-1)+(j  +gc-1)*pi+(k+1+gc-1)*pi*pj+1
    c6 =boff+(i+1+gc-1)+(j  +gc-1)*pi+(k+1+gc-1)*pi*pj+1
    c7 =boff+(i  +gc-1)+(j-1+gc-1)*pi+(k  +gc-1)*pi*pj+1
    c8 =boff+(i  +gc-1)+(j+1+gc-1)*pi+(k  +gc-1)*pi*pj+1
    c9 =boff+(i  +gc-1)+(j-1+gc-1)*pi+(k+1+gc-1)*pi*pj+1
    c10=boff+(i  +gc-1)+(j+1+gc-1)*pi+(k+1+gc-1)*pi*pj+1
    do s = 1, nprim
      p1(s)=w(c1,s); p2(s)=w(c2,s)
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
    M1l = Mp(:,:,c1); M2l = Mp(:,:,c2)
    fi = tofm(p) + (i-1) + n1*((j-1) + n2*k) + 1
    nrm(1)=fm3p(1,fi); nrm(2)=fm3p(2,fi); nrm(3)=fm3p(3,fi)
    area = fm3p(4,fi); wd1 = ynp(c1); wd2 = ynp(c2)
    f = tofh(3,p) + (i-1) + (j-1)*n1 + (k-1)*n1*n2 + 1
    rho_f=tfp(f,1); T_f=tfp(f,2); mil_f=tfp(f,3); kl_f=tfp(f,4); cp_f=tfp(f,5)
#define DCUF_TF tfp
#define DCUF_ROW f
#include "diff_physics.inc"
#undef DCUF_TF
#undef DCUF_ROW
    do s = 1, nprim
      dfhp(f, s) = res0(s)
    end do
  end subroutine DCUF_PKK
