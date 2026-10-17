! ---------------------------------------------------------------------------------------------
! The three CUF pass-1 kernels of the diffusive flux, instantiated twice in
! Mod_Diffusive.f90:
!   - MOSE_MCDIFF_KERNEL undefined -> the production kernels diff_cuf_* : the previous source token for token
!     (tface of 5 columns, Dm = (mu/Sc + mu_t/Sct)/rho in diff_physics.inc);
!   - MOSE_MCDIFF_KERNEL defined   -> the twins diff_cuf_*_mc : tface of 5+nsc columns and, in diff_physics.inc,
!     Dm(s) = the face D_k (column 5+s, written by the fills) + mu_t/(rho Sct), as MOSE_open with schmidt <= 0.
! The launchers pick the twin on mcdiff_dev. Expected macros: the three kernel names (DCUF_K + I, J, K); each
! kernel names its face-thermo array and row for diff_physics.inc right before including it.
! ---------------------------------------------------------------------------------------------
  attributes(global) launch_bounds(128) subroutine DCUF_KI ( &
       w, dfhat, fm1, M, yn, rho_c, T_c, &
       tface, &
       ncell, nfacemax, &
       boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3, &
       RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
    integer, value :: ncell, nfacemax, boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3
    integer, value :: stress_dev, eddy_dev
    logical, value :: RANS
    real(R8), value :: sc_dev, sct_dev, prt_dev
    real(R8) :: w(ncell, nprim), dfhat(nfacemax, nprim)
    real(R8) :: fm1(4, 0:n1, n2, n3)
#ifdef MOSE_MCDIFF_KERNEL
    real(R8) :: tface(nfacemax, 5+nsc)   ! + DIFF_MC twin: the face D_k in columns 6:5+nsc (fills)
#else
    real(R8) :: tface(nfacemax, 5)   ! face-exact rho,T,mil,kl,cp (Diffusive_Face_Thermo_Fill)
#endif
    real(R8) :: M(3, 3, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    real(R8) :: yn(1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    real(R8) :: rho_c(ncell), T_c(ncell)   ! the flat twins rho_c_f / T_c_f (the same doubles as the per-block caches)
    integer :: i, j, k, s, f, ia, ib, ic
    real(R8) :: p1(NPRIM_DEV),p2(NPRIM_DEV),p3(NPRIM_DEV),p4(NPRIM_DEV),p5(NPRIM_DEV)
    real(R8) :: p6(NPRIM_DEV),p7(NPRIM_DEV),p8(NPRIM_DEV),p9(NPRIM_DEV),p10(NPRIM_DEV)
    real(R8) :: rho_n10(10), T_n10(10), M1l(3,3), M2l(3,3), nrm(3), res0(NPRIM_DEV)
    real(R8) :: area, wd1, wd2, rho_f, T_f, mil_f, kl_f, cp_f
    integer  :: c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
    integer :: t_, ni_, nj_, rem_
    include 'diff_decls.inc'
    ! flat 1D mapping: in-kernel decode, i fastest -> coalescing preserved;
    ! kills the (32,2,2) tile's partial-warp waste on narrow blocks. Bit-exact.
    t_  = (blockIdx%x-1)*blockDim%x + threadIdx%x
    ni_ = h1 - o1 + 1
    nj_ = h2 - o2 + 1
    if (t_ > ni_*nj_*(h3-o3+1)) return
    i    = mod(t_-1, ni_) + o1
    rem_ = (t_-1) / ni_
    j    = mod(rem_, nj_) + o2
    k    = rem_ / nj_ + o3
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
    M1l = M(:,:,i,j,k); M2l = M(:,:,i+1,j,k)
    nrm(1)=fm1(1,i,j,k); nrm(2)=fm1(2,i,j,k); nrm(3)=fm1(3,i,j,k)
    area = fm1(4,i,j,k); wd1 = yn(i,j,k); wd2 = yn(i+1,j,k)
    f = (i-1) + (j-1)*(n1-1) + (k-1)*(n1-1)*n2 + 1
    rho_f=tface(f,1); T_f=tface(f,2); mil_f=tface(f,3); kl_f=tface(f,4); cp_f=tface(f,5)
#define DCUF_TF tface
#define DCUF_ROW f
#include "diff_physics.inc"
#undef DCUF_TF
#undef DCUF_ROW
    f = (i-1) + (j-1)*(n1-1) + (k-1)*(n1-1)*n2 + 1
    do s = 1, nprim
      dfhat(f, s) = res0(s)
    end do
  end subroutine DCUF_KI

  attributes(global) launch_bounds(128) subroutine DCUF_KJ ( &
       w, dfhat, fm2, M, yn, rho_c, T_c, &
       tface, &
       ncell, nfacemax, &
       boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3, &
       RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
    integer, value :: ncell, nfacemax, boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3
    integer, value :: stress_dev, eddy_dev
    logical, value :: RANS
    real(R8), value :: sc_dev, sct_dev, prt_dev
    real(R8) :: w(ncell, nprim), dfhat(nfacemax, nprim)
    real(R8) :: fm2(4, n1, 0:n2, n3)
#ifdef MOSE_MCDIFF_KERNEL
    real(R8) :: tface(nfacemax, 5+nsc)   ! + DIFF_MC twin: the face D_k in columns 6:5+nsc (fills)
#else
    real(R8) :: tface(nfacemax, 5)   ! face-exact rho,T,mil,kl,cp (Diffusive_Face_Thermo_Fill)
#endif
    real(R8) :: M(3, 3, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    real(R8) :: yn(1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    real(R8) :: rho_c(ncell), T_c(ncell)   ! the flat twins rho_c_f / T_c_f (the same doubles as the per-block caches)
    integer :: i, j, k, s, f, ia, ib, ic
    real(R8) :: p1(NPRIM_DEV),p2(NPRIM_DEV),p3(NPRIM_DEV),p4(NPRIM_DEV),p5(NPRIM_DEV)
    real(R8) :: p6(NPRIM_DEV),p7(NPRIM_DEV),p8(NPRIM_DEV),p9(NPRIM_DEV),p10(NPRIM_DEV)
    real(R8) :: rho_n10(10), T_n10(10), M1l(3,3), M2l(3,3), nrm(3), res0(NPRIM_DEV)
    real(R8) :: area, wd1, wd2, rho_f, T_f, mil_f, kl_f, cp_f
    integer  :: c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
    integer :: t_, ni_, nj_, rem_
    include 'diff_decls.inc'
    ! flat 1D mapping: in-kernel decode, i fastest -> coalescing preserved;
    ! kills the (32,2,2) tile's partial-warp waste on narrow blocks. Bit-exact.
    t_  = (blockIdx%x-1)*blockDim%x + threadIdx%x
    ni_ = h1 - o1 + 1
    nj_ = h2 - o2 + 1
    if (t_ > ni_*nj_*(h3-o3+1)) return
    i    = mod(t_-1, ni_) + o1
    rem_ = (t_-1) / ni_
    j    = mod(rem_, nj_) + o2
    k    = rem_ / nj_ + o3
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
    M1l = M(:,:,i,j,k); M2l = M(:,:,i,j+1,k)
    nrm(1)=fm2(1,i,j,k); nrm(2)=fm2(2,i,j,k); nrm(3)=fm2(3,i,j,k)
    area = fm2(4,i,j,k); wd1 = yn(i,j,k); wd2 = yn(i,j+1,k)
    f = (i-1) + (j-1)*n1 + (k-1)*n1*(n2-1) + 1
    rho_f=tface(f,1); T_f=tface(f,2); mil_f=tface(f,3); kl_f=tface(f,4); cp_f=tface(f,5)
#define DCUF_TF tface
#define DCUF_ROW f
#include "diff_physics.inc"
#undef DCUF_TF
#undef DCUF_ROW
    f = (i-1) + (j-1)*n1 + (k-1)*n1*(n2-1) + 1
    do s = 1, nprim
      dfhat(f, s) = res0(s)
    end do
  end subroutine DCUF_KJ

  attributes(global) launch_bounds(128) subroutine DCUF_KK ( &
       w, dfhat, fm3, M, yn, rho_c, T_c, &
       tface, &
       ncell, nfacemax, &
       boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3, &
       RANS, stress_dev, eddy_dev, sc_dev, sct_dev, prt_dev)
    integer, value :: ncell, nfacemax, boff, pi, pj, gc, n1, n2, n3, o1, o2, o3, h1, h2, h3
    integer, value :: stress_dev, eddy_dev
    logical, value :: RANS
    real(R8), value :: sc_dev, sct_dev, prt_dev
    real(R8) :: w(ncell, nprim), dfhat(nfacemax, nprim)
    real(R8) :: fm3(4, n1, n2, 0:n3)
#ifdef MOSE_MCDIFF_KERNEL
    real(R8) :: tface(nfacemax, 5+nsc)   ! + DIFF_MC twin: the face D_k in columns 6:5+nsc (fills)
#else
    real(R8) :: tface(nfacemax, 5)   ! face-exact rho,T,mil,kl,cp (Diffusive_Face_Thermo_Fill)
#endif
    real(R8) :: M(3, 3, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    real(R8) :: yn(1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    real(R8) :: rho_c(ncell), T_c(ncell)   ! the flat twins rho_c_f / T_c_f (the same doubles as the per-block caches)
    integer :: i, j, k, s, f, ia, ib, ic
    real(R8) :: p1(NPRIM_DEV),p2(NPRIM_DEV),p3(NPRIM_DEV),p4(NPRIM_DEV),p5(NPRIM_DEV)
    real(R8) :: p6(NPRIM_DEV),p7(NPRIM_DEV),p8(NPRIM_DEV),p9(NPRIM_DEV),p10(NPRIM_DEV)
    real(R8) :: rho_n10(10), T_n10(10), M1l(3,3), M2l(3,3), nrm(3), res0(NPRIM_DEV)
    real(R8) :: area, wd1, wd2, rho_f, T_f, mil_f, kl_f, cp_f
    integer  :: c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
    integer :: t_, ni_, nj_, rem_
    include 'diff_decls.inc'
    ! flat 1D mapping: in-kernel decode, i fastest -> coalescing preserved;
    ! kills the (32,2,2) tile's partial-warp waste on narrow blocks. Bit-exact.
    t_  = (blockIdx%x-1)*blockDim%x + threadIdx%x
    ni_ = h1 - o1 + 1
    nj_ = h2 - o2 + 1
    if (t_ > ni_*nj_*(h3-o3+1)) return
    i    = mod(t_-1, ni_) + o1
    rem_ = (t_-1) / ni_
    j    = mod(rem_, nj_) + o2
    k    = rem_ / nj_ + o3
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
    M1l = M(:,:,i,j,k); M2l = M(:,:,i,j,k+1)
    nrm(1)=fm3(1,i,j,k); nrm(2)=fm3(2,i,j,k); nrm(3)=fm3(3,i,j,k)
    area = fm3(4,i,j,k); wd1 = yn(i,j,k); wd2 = yn(i,j,k+1)
    f = (i-1) + (j-1)*n1 + (k-1)*n1*n2 + 1
    rho_f=tface(f,1); T_f=tface(f,2); mil_f=tface(f,3); kl_f=tface(f,4); cp_f=tface(f,5)
#define DCUF_TF tface
#define DCUF_ROW f
#include "diff_physics.inc"
#undef DCUF_TF
#undef DCUF_ROW
    f = (i-1) + (j-1)*n1 + (k-1)*n1*n2 + 1
    do s = 1, nprim
      dfhat(f, s) = res0(s)
    end do
  end subroutine DCUF_KK
