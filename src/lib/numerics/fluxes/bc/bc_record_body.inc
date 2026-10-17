! One boundary record of BC_Fluxes_blk (idx set by the caller's loop).
! Shared by both loops so the box and the list launches do the same work per record.
        Im = bc_bi(idx); Jm = bc_bj(idx); Km = bc_bk(idx); Fm = bc_bf(idx)
        dir   = (Fm + 1)/2
        modfm = mod(Fm, 2)
        gstep = 2*modfm - 1

        ! Boundary-face index + metric (fm<dir> includes the boundary face at 0:n)
        Fi = Im; Fj = Jm; Fk = Km
        Ii = Im; Ij = Jm; Ik = Km
        select case (dir)
          case (1)
            Fi = Im - modfm; Ii = Im + gstep
            normal(1) = fm1(1,Fi,Fj,Fk); normal(2) = fm1(2,Fi,Fj,Fk); normal(3) = fm1(3,Fi,Fj,Fk)
            area = fm1(4,Fi,Fj,Fk)
          case (2)
            Fj = Jm - modfm; Ij = Jm + gstep
            normal(1) = fm2(1,Fi,Fj,Fk); normal(2) = fm2(2,Fi,Fj,Fk); normal(3) = fm2(3,Fi,Fj,Fk)
            area = fm2(4,Fi,Fj,Fk)
          case (3)
            Fk = Km - modfm; Ik = Km + gstep
            normal(1) = fm3(1,Fi,Fj,Fk); normal(2) = fm3(2,Fi,Fj,Fk); normal(3) = fm3(3,Fi,Fj,Fk)
            area = fm3(4,Fi,Fj,Fk)
        end select
        Bound_Dl = 0.5d0 * dl_f(dir, Im, Jm, Km)
        Int_Dl   = 0.5d0 * dl_f(dir, Ii, Ij, Ik)
        ! 5r-b: the ghost cell and the four dl points the connection twins used to index themselves
        Ig_ = Im; Jg_ = Jm; Kg_ = Km
        if ( dir == 1 ) Ig_ = Im - gstep
        if ( dir == 2 ) Jg_ = Jm - gstep
        if ( dir == 3 ) Kg_ = Km - gstep
        o1_ = modfm*(-gstep*gc) + (1-modfm)*( gstep)
        o2_ = modfm*(-gstep)
        o3_ = (1-modfm)*(-gstep)
        o4_ = (1-modfm)*(-gstep*gc) + modfm*( gstep)
        if ( dir == 1 ) then
          dlc4_(1) = dl_f(1, Im+o1_, Jm, Km); dlc4_(2) = dl_f(1, Im+o2_, Jm, Km)
          dlc4_(3) = dl_f(1, Im+o3_, Jm, Km); dlc4_(4) = dl_f(1, Im+o4_, Jm, Km)
        else if ( dir == 2 ) then
          dlc4_(1) = dl_f(2, Im, Jm+o1_, Km); dlc4_(2) = dl_f(2, Im, Jm+o2_, Km)
          dlc4_(3) = dl_f(2, Im, Jm+o3_, Km); dlc4_(4) = dl_f(2, Im, Jm+o4_, Km)
        else
          dlc4_(1) = dl_f(3, Im, Jm, Km+o1_); dlc4_(2) = dl_f(3, Im, Jm, Km+o2_)
          dlc4_(3) = dl_f(3, Im, Jm, Km+o3_); dlc4_(4) = dl_f(3, Im, Jm, Km+o4_)
        end if

        ! 5r-c: the values the branches used to index out of the per-block arrays
        Mrec_ = 0.0d0; ynrec_ = 0.0d0; Mgho_ = 0.0d0; yngho_ = 0.0d0; beta_ = 1.0d0
        if ( bc_typ(idx) == 1 .or. bc_typ(idx) == 666 ) then
          beta_  = beta(Im,Jm,Km)
          Mrec_  = m_f(:,:,Im,Jm,Km);  ynrec_ = yn(Im,Jm,Km)
          Mgho_  = m_f(:,:,Ig_,Jg_,Kg_); yngho_ = yn(Ig_,Jg_,Kg_)
        else if ( bc_typ(idx) == 3 .or. bc_typ(idx) == 5 .or. bc_typ(idx) == 6 ) then
          Mrec_  = m_f(:,:,Im,Jm,Km);  ynrec_ = yn(Im,Jm,Km)
        end if
#include "bc_record_apply.inc"
