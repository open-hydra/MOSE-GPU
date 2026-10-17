module MOSE_Mod_dt
  use iso_fortran_env, only: I4 => int32, R8 => real64

  use MOSE_Global_m,   only: gc, nprim, NPRIM_DEV, nsc, nu, nw, np, nt, nrans, ndir, VISC, RANS
  use MOSE_Field_m,    only: P_flat, field_block_index_params
  implicit none

  public :: Compute_dt_blk, Compute_dt_blk_mc, Dt_Anatomy_blk

contains


  subroutine Compute_dt ( domain, config )
    use MOSE_Types_m
    use MOSE_Global_m
    use MOSE_Lib_RANS
    use MOSE_Mod_Thermo_Cache, only: Refresh_Thermo_Cache_Block
    use MOSE_Mod_MPI, only: is_local_block, mpi_allreduce_min_r8, mpi_rank_
    use MOSE_IO_Table, only: vnn_diff   ! DIFF_MC A' (D7): [MOSE-parameters] vnn-diff, read with schmidt <= 0 only
    use MOSE_Field_m, only: field_is_allocated, field_alloc, pack_P_to_flat, field_update_P_flat_device
    use MOSE_Mod_Thermo_Cache, only: Thermo_Cache_rrt_blk, Thermo_Cache_visc_blk
    use MOSE_Lib_Preconditioning, only: obj_prec, ur_f, Prec_Alloc, Update_Derived_Variables
    implicit none

    type(MOSE_domain_type), intent(inout) :: domain
    type(config_type), intent(inout)      :: config
    ! Local
    integer :: i, j, k, b
    real(R8) :: dtcell, dtglobal, rampa_fac
    logical, save :: w_init_done = .false.   ! One-time seed of resident P_flat from the IC P
    character(len=8)  :: an_env
    logical, save :: an_done = .false.       ! MOSE_DT_ANAT: one-shot per-direction dt anatomy
    integer, save :: an_ncall = 0
    integer :: an_want, an_u, an_t, an_st, ii, jj, kk, an_n, an_m
    real(R8), allocatable :: an(:,:,:,:)
    character(len=64) :: an_fn


    dtglobal = domain % dtglobal

    !> Device path: per-block thermo refresh (device) + per-cell dt kernel with
    !> a min-reduction (order-independent -> no FP-reorder concern). The ramp
    !> factor is evaluated on the host and passed in.
    rampa_fac = 1d0
    if ( domain%iter < config%rampa_iter ) rampa_fac = real(domain%iter,R8)/config%rampa_iter
    ! Timestep + its interior thermo-cache refresh on the var-LAST
    ! flat SoA w. Compute_dt runs at the TOP of Explicit_Step, BEFORE the P/cache/geometry
    ! enter-data residency and before any Fill_Ghost — host P is current here (synced by
    ! the prior step's f5_P_to_host, or the IC on iter 1), so pack directly from host P.
    ! NB no `update host(P)`: P may not be device-resident yet on iter 1 and that would
    ! fault (the  lesson). Cache is interior-only (ghost prims uninitialised before
    ! Fill_Ghost -> OOB in the Wilke int(T) dispatch). Caches/geometry use present_or data
    ! clauses inside the kernels so this works both pre- and post-residency.
    if (.not. field_is_allocated()) call field_alloc(domain)
    ! ONE-TIME seed of the resident P_flat from the host IC P. Thereafter w is the
    ! device-resident truth (maintained by Newstate/chemistry); Compute_dt reads it -> no
    ! per-step P -> P_flat host pack (that pack was a per-step host transpose at THOR scale).
    if (.not. w_init_done) then
      call pack_P_to_flat(domain); call field_update_P_flat_device(); w_init_done = .true.
    end if
    ! integration-variables = prec: Ur of the starting state before the first time step (MOSE_open Wrap_Setup), then after
    ! every state update (Mod_Explicit); ur_f has one element otherwise
    call Prec_Alloc ()
    if ( obj_prec%enabled ) call Update_Derived_Variables ( domain, .true. )
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      call Thermo_Cache_rrt_blk ( b, domain%blk(b)%dim, 1, domain%blk(b)%dim(1),   &
             1, domain%blk(b)%dim(2), 1, domain%blk(b)%dim(3),                        &
             domain%blk(b)%rho_c, domain%blk(b)%Rg_c, domain%blk(b)%T_c )
      if (VISC) &
        call Thermo_Cache_visc_blk ( b, domain%blk(b)%dim, domain%blk(b)%rho_c,     &
               domain%blk(b)%T_c, domain%blk(b)%mil_c )
      if ( mcdiff_dev .and. vnn_diff > 0d0 ) then   ! DIFF_MC A' (D7): + the species-diffusion limit
      call Compute_dt_blk_mc ( b, domain%blk(b)%m_f, domain%blk(b)%dl_f,            &
             domain%blk(b)%rho_c, domain%blk(b)%Rg_c, domain%blk(b)%mil_c,            &
             domain%blk(b)%dtlocal, domain%blk(b)%dim, config%cfl, config%vnn,        &
             rampa_fac, dtglobal, ur_f, size(ur_f), obj_prec%enabled, vnn_diff, sct_dev )
      else
      call Compute_dt_blk ( b, domain%blk(b)%m_f, domain%blk(b)%dl_f,               &
             domain%blk(b)%rho_c, domain%blk(b)%Rg_c, domain%blk(b)%mil_c,            &
             domain%blk(b)%dtlocal, domain%blk(b)%dim, config%cfl, config%vnn,        &
             rampa_fac, dtglobal, ur_f, size(ur_f), obj_prec%enabled )
      end if
    end do
    call mpi_allreduce_min_r8(dtglobal, domain % dtglobal)
    ! MOSE_DT_ANAT=<n>: at the n-th Compute_dt call, the per-direction ingredients of every cell's
    ! limit (diagnostic only). OUTPUT/dtanat_r<rank>.bin: per local block int32 b, n(3), then
    ! real32 an(1:6,:,:,:) = convective times dl/lambda (3) and viscous times rho dl^2/mu (3), then an(13) = hr_c.
    ! OUTPUT/dtanat_r<rank>.txt: every cell below 4 dt_min with geometry and state.
    if (.not. an_done) then
      an_ncall = an_ncall + 1
      call get_environment_variable('MOSE_DT_ANAT', an_env, status=an_st)
      an_want = 0
      if (an_st == 0 .or. an_st == -1) read(an_env,*,iostat=an_st) an_want
      if (an_st /= 0) an_want = 0
      if (an_want > 0 .and. an_ncall >= an_want) then
        an_done = .true.
        write(an_fn,'(A,I0,A)') 'OUTPUT/dtanat_r', mpi_rank_, '.bin'
        open(newunit=an_u, file=trim(an_fn), access='stream', form='unformatted', status='replace')
        write(an_fn,'(A,I0,A)') 'OUTPUT/dtanat_r', mpi_rank_, '.txt'
        open(newunit=an_t, file=trim(an_fn), status='replace')
        write(an_t,'(A,ES16.8,A,I0,A,F8.4,A,F8.4)') '# dt_min ', domain%dtglobal, '  iter ', domain%iter, &
          '  cfl ', config%cfl, '  vnn ', config%vnn
        write(an_t,'(A)') '# b i j k  dtloc/dtmin  cfl*tc_i cfl*tc_j cfl*tc_k (/dtmin)  vnn*tv_i vnn*tv_j vnn*tv_k (/dtmin)' // &
          '  dl_i dl_j dl_k  yn  T  a  |u_i| |u_j| |u_k|'
        ! DIFF_MC (schmidt <= 0): 6 more rows, the species-diffusion time td = dl^2/max_k(D_k + mu_t/(rho Sct)) of A' per
        ! direction, the species k of the max, its Y_k and the max; OUTPUT/dtanat_mc_r<rank>.bin: per local block int32 b,
        ! n(3), then real32 an(14:19,:,:,:); in the .txt six more columns (td/dtmin without the vnn-diff factor, k, Y_k, Dmax).
        an_n = 13
        if ( mcdiff_dev ) then
          an_n = 19
          write(an_t,'(A,F8.4)') '# DIFF_MC columns: td_i td_j td_k (/dtmin, no vnn-diff factor)  k_max  Y_kmax  Dmax;  vnn-diff ', &
            vnn_diff
          write(an_fn,'(A,I0,A)') 'OUTPUT/dtanat_mc_r', mpi_rank_, '.bin'
          open(newunit=an_m, file=trim(an_fn), access='stream', form='unformatted', status='replace')
        end if
        do b = 1, domain % nb
          if (.not. is_local_block(b)) cycle
          allocate ( an(an_n, domain%blk(b)%dim(1), domain%blk(b)%dim(2), domain%blk(b)%dim(3)) )
          call Dt_Anatomy_blk ( b, domain%blk(b)%m_f, domain%blk(b)%dl_f, domain%blk(b)%rho_c, &
                 domain%blk(b)%Rg_c, domain%blk(b)%mil_c, domain%blk(b)%dtlocal, domain%blk(b)%dim, an, &
                 domain%blk(b)%hr_c, an_n, sct_dev )
          write(an_u) int(b,I4), int(domain%blk(b)%dim(1:3),I4), real(an(1:6,:,:,:), kind=4), &
                      real(an(13,:,:,:), kind=4)
          if ( mcdiff_dev ) write(an_m) int(b,I4), int(domain%blk(b)%dim(1:3),I4), real(an(14:19,:,:,:), kind=4)
          do kk = 1, domain%blk(b)%dim(3); do jj = 1, domain%blk(b)%dim(2); do ii = 1, domain%blk(b)%dim(1)
            if ( an(12,ii,jj,kk) < 4.0_R8*domain%dtglobal ) then
              if ( mcdiff_dev ) then
              write(an_t,'(4I5,F10.4,6F10.3,3ES12.4,ES12.4,F9.1,F9.1,3F9.1,3F10.3,I4,2ES12.4)') b, ii, jj, kk, &
                an(12,ii,jj,kk)/domain%dtglobal, config%cfl*an(1:3,ii,jj,kk)/domain%dtglobal, &
                config%vnn*an(4:6,ii,jj,kk)/domain%dtglobal, domain%blk(b)%dl_f(1:3,ii,jj,kk), &
                domain%blk(b)%yn(ii,jj,kk), an(11,ii,jj,kk), an(10,ii,jj,kk), an(7:9,ii,jj,kk), &
                an(14:16,ii,jj,kk)/domain%dtglobal, nint(an(17,ii,jj,kk)), an(18:19,ii,jj,kk)
              else
              write(an_t,'(4I5,F10.4,6F10.3,3ES12.4,ES12.4,F9.1,F9.1,3F9.1)') b, ii, jj, kk, &
                an(12,ii,jj,kk)/domain%dtglobal, config%cfl*an(1:3,ii,jj,kk)/domain%dtglobal, &
                config%vnn*an(4:6,ii,jj,kk)/domain%dtglobal, domain%blk(b)%dl_f(1:3,ii,jj,kk), &
                domain%blk(b)%yn(ii,jj,kk), an(11,ii,jj,kk), an(10,ii,jj,kk), an(7:9,ii,jj,kk)
              end if
            end if
          end do; end do; end do
          deallocate ( an )
        end do
        close(an_u); close(an_t)
        if ( mcdiff_dev ) close(an_m)
      end if
    end if
    call Clip_dt_To_Timemax ( domain, config )
    return
    contains

      subroutine compute ( rhoi, vel, p, rans_, met, dl, rho, Rgas, mil, dtmin, cfl, vnn )
        use MOSE_Global_m
        use FLINT_Lib_Thermodynamic, only: f_ss, co_DS_expr

        implicit none
        real(R8), intent(in)  :: rhoi(nsc), vel(3), p, rans_(:)
        real(R8), intent(in)  :: met(3,3), dl(3), cfl, vnn
        real(R8), intent(in)  :: rho, Rgas, mil
        real(R8), intent(out) :: dtmin
        ! Local
        integer :: d
        real(R8) :: Sound, dt, versor(3), lambda, mie, mi
        real(R8) :: vel_, dummy(3,3)=0d0
        integer  :: Tint(2)
        real(R8) :: Tc_, Tdf_, Dk(nsc), Dmax_

        Sound = f_ss ( rhoi, p, rho, Rgas )

        dtmin = 1d8

        do d = 1, ndir

          ! CFL condition along d-direction
          versor = met(d,:) / norm2 ( met(d,:) )
          vel_ = abs ( dot_product ( vel, versor ) )

          lambda  = vel_ + Sound
          dt = dl(d) / lambda * cfl
          dtmin = min ( dt, dtmin )

        enddo

        if ( VISC ) then
          mie = 0.d0
          if (RANS) then
            ! Note: approximate mit for 2 equation models.
            ! velocity gradient assumed 0 and small distance from wall
            call Eddy_Viscosity ( mut=mie, rans_variables=rans_, mul=mil, rho=rho, &
                                  vel_gradient=dummy, walldist=1d-6 )
          endif
          mi = mie + mil

          do d = 1, ndir

            ! VNN condition along d-direction
            dt = ( rho * dl(d)**2 * vnn ) / mi
            dtmin = min ( dt, dtmin )

          end do
          ! DIFF_MC A' (D7): species-diffusion limit vnn-diff * dl^2 / max_k(D_k + mu_t/(rho Sct))
          if ( mcdiff_dev .and. vnn_diff > 0d0 ) then
            Tc_  = p / ( rho * Rgas )
            Tint(1) = int(Tc_)
            Tint(2) = Tint(1) + 1
            Tdf_ = Tc_ - Tint(1)
            call co_DS_expr ( rhoi, rho, Tint, Tdf_, p, Dk )
            Dmax_ = maxval(Dk) + mie/(rho*sct_dev)
            do d = 1, ndir
              dt = ( vnn_diff * dl(d)**2 ) / Dmax_
              dtmin = min ( dt, dtmin )
            end do
          end if
        endif

      end subroutine compute

  end subroutine Compute_dt


  !> time-clip: shorten the global dt so that the step lands exactly on
  !> time-threshold (same-time field comparisons against an oracle). Off by
  !> default; only in time-accurate runs. The 1e-12 overshoot makes the
  !> `time >= timemax` stop test robust to rounding.
  subroutine Clip_dt_To_Timemax ( domain, config )
    use MOSE_Types_m
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    type(config_type), intent(in)         :: config
    real(R8) :: rem
    if ( .not. ( config%time_clip .and. config%time_accurate ) ) return
    rem = config%timemax - domain%time
    if ( rem > 0d0 .and. domain%dtglobal > rem ) domain%dtglobal = rem * ( 1d0 + 1d-12 )
    if ( config%dump_clip ) then   ! dump-clip (A5): land exactly on the next time-based output instant as well
      rem = config%time_from_call + config%dump_dt - domain%time
      if ( rem > 0d0 .and. domain%dtglobal > rem ) domain%dtglobal = rem * ( 1d0 + 1d-12 )
    end if
  end subroutine Clip_dt_To_Timemax


  !----------------------------------------------------------------------------
  ! GPU per-cell local time step + min reduction (device twin of Compute_dt's
  ! inner loop). Flat geometry m_f/dl_f (== old m%c / dl%c); f_ss_dev sound
  ! speed; devirtualized SST eddy viscosity for the VNN viscous limit (same
  ! Gradvel=0 / dist=1e-6 approximation as the host). min() is associative so
  ! the reduction is order-independent (no FP-reorder drift).
  !----------------------------------------------------------------------------

  subroutine Set_Global_dt ( domain )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, k

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle

      ! dtlocal is device-resident. The uniform global-dt write must happen ON
      ! the device (a host write would never reach the resident dtlocal that device
      ! Newstate/chemistry read) — this device port is precisely what unblocks dtlocal
      ! residency. One dummy level (component in) keeps the present-match; copy(dtlocal)
      ! is present_or so it also works on iter 1 before enter-data makes it resident.
      call Set_Global_dt_blk ( domain % blk(b) % dtlocal, domain % dtglobal, domain % blk(b) % dim, b )

    end do

  end subroutine Set_Global_dt


  !> Device twin of Set_Global_dt's per-block uniform write dtlocal(:,:,:) = dtglobal.
  !> dtlocal is resident (enter data); copy() is present_or so it is correct on iter 1
  !> (before residency) and a no-op in-place thereafter. dtg is passed by value.
  subroutine Set_Global_dt_blk ( dtlocal, dtg, n, b )
    use MOSE_Global_m, only: gc
    use MOSE_Field_m,  only: dt_f, n_cache_f, field_block_index_params
    implicit none
    integer,  intent(in)    :: n(3), b
    real(R8), intent(in)    :: dtg
    real(R8), intent(inout) :: dtlocal(1:n(1), 1:n(2), 1:n(3))
    integer :: i, j, k, nb1_, nb2_, nb3_, boff_, pi_, pj_
    real(R8) :: dtg_
    nb1_ = n(1); nb2_ = n(2); nb3_ = n(3); dtg_ = dtg
    call field_block_index_params(b, boff_, pi_, pj_)
    ! Dummy scalars/arrays referenced in the region are
    ! implicitly copied per launch; locals ride the launch parameters.
    !$acc parallel loop collapse(3) gang vector copy(dtlocal) present(dt_f)
    do k = 1, nb3_
    do j = 1, nb2_
    do i = 1, nb1_
      dtlocal(i,j,k) = dtg_
      if ( n_cache_f > 1 ) dt_f( boff_ + (i+gc-1) + (j+gc-1)*pi_ + (k+gc-1)*pi_*pj_ + 1 ) = dtg_
    end do
    end do
    end do
  end subroutine Set_Global_dt_blk



  ! ---- from Mod_dt_W.f90 (MOSE_Mod_dt) ----

  !> The routine lives in dt_blk.inc, instantiated as Compute_dt_blk (the previous source) and
  !> Compute_dt_blk_mc (+ the species-diffusion limit A', chosen by Compute_dt with schmidt <= 0 and vnn-diff > 0).
#define DTW_NAME Compute_dt_blk
#include "dt_blk.inc"
#undef DTW_NAME
#define MOSE_MCDIFF_KERNEL
#define DTW_NAME Compute_dt_blk_mc
#include "dt_blk.inc"
#undef DTW_NAME
#undef MOSE_MCDIFF_KERNEL

  !> Diagnostic (MOSE_DT_ANAT, the non-preconditioned speeds also with prec): the ingredients of Compute_dt_blk's limit, per cell and per
  !> direction, WITHOUT the cfl/vnn factors:
  !>   an(1:3)  = dl_f(d)/lambda_d        convective time of direction d
  !>   an(4:6)  = rho dl_f(d)^2 / mu      viscous time of direction d (mu = laminar + SST, as the solver)
  !>   an(7:9)  = |u . versor_d|
  !>   an(10)   = sound speed, an(11) = T, an(12) = dtlocal as the solver computed it, an(13) = hr_c
  !>   DIFF_MC (nan = 19, schmidt <= 0): an(14:16) = dl_f(d)^2 / max_k(D_k + mu_t/(rho Sct)) species-diffusion time of
  !>   direction d (without vnn-diff), an(17) = the species k of the max, an(18) = its Y_k, an(19) = max_k D_k + mu_t/(rho Sct)
  subroutine Dt_Anatomy_blk ( b, m_f, dl_f, rho_c, Rg_c, mil_c, dtlocal, n, an, hr_c, nan, sctd )
    use MOSE_Mod_SST, only: SST_Eddy_Viscosity
    use MOSE_Mod_Spalart, only: Spalart_Eddy_Viscosity
    use MOSE_Mod_Wilcox2006, only: Wilcox2006_Eddy_Viscosity
    use MOSE_Lib_RANS, only: eddy_dev
    use FLINT_Lib_ThermoTransport_dev, only: f_ss_dev, co_DS_expr_dev
    integer,  intent(in) :: b, n(3), nan
    real(R8), intent(in) :: sctd
    real(R8), dimension(3, 3, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: m_f
    real(R8), dimension(3, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: dl_f
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: rho_c, Rg_c, mil_c
    real(R8), dimension(1:n(1), 1:n(2), 1:n(3)), intent(in) :: dtlocal
    real(R8), intent(out) :: an(nan, n(1), n(2), n(3))
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: hr_c
    integer  :: i, j, k, d, s, c0, boff, pi, pj, n1_, n2_, n3_
    real(R8) :: Sound, lambda, vel_, mie, mi
    real(R8) :: pl(nprim), vel(3), versor(3), Gvel0(3,3), rkw(2)
    integer  :: T_i, Tint(2), km_
    real(R8) :: Tc_, Tdf_, Dmax_, Dk(NPRIM_DEV)
    logical  :: mcd_
    call field_block_index_params(b, boff, pi, pj)
    n1_ = n(1); n2_ = n(2); n3_ = n(3)
    mcd_ = ( nan >= 20 )
    !$acc data copyin(m_f, dl_f, rho_c, Rg_c, mil_c, dtlocal, hr_c) copyout(an)
    !$acc parallel loop collapse(3) gang vector present(P_flat, m_f, dl_f, rho_c, Rg_c, mil_c, dtlocal, an, hr_c) &
    !$acc&  private(pl, vel, versor, Gvel0, rkw, c0, s, d, Sound, lambda, vel_, mie, mi, T_i, Tint, km_, Tc_, Tdf_, Dmax_, Dk)
    do k = 1, n3_
    do j = 1, n2_
    do i = 1, n1_
      c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      do s = 1, nprim
        pl(s) = P_flat(c0,s)
      end do
      Sound    = f_ss_dev( pl(1:nsc), pl(np), rho_c(i,j,k), Rg_c(i,j,k) )
      vel(1:3) = pl(nu:nw)
      mie = 0d0
      if (VISC .and. RANS) then
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
      do d = 1, 3
        versor(1:3) = m_f(d,1:3,i,j,k) / norm2( m_f(d,1:3,i,j,k) )
        vel_   = abs( dot_product( vel(1:3), versor(1:3) ) )
        lambda = vel_ + Sound
        an(d,i,j,k)   = dl_f(d,i,j,k) / lambda
        an(3+d,i,j,k) = rho_c(i,j,k) * dl_f(d,i,j,k)**2 / mi
        an(6+d,i,j,k) = vel_
      end do
      an(10,i,j,k) = Sound
      an(11,i,j,k) = pl(np) / ( rho_c(i,j,k) * Rg_c(i,j,k) )
      an(12,i,j,k) = dtlocal(i,j,k)
      an(13,i,j,k) = hr_c(i,j,k)
      if (mcd_) then   ! DIFF_MC: the species-diffusion time of A' and the species that sets it
        Tc_  = an(11,i,j,k)
        T_i  = int(Tc_)
        Tdf_ = Tc_ - T_i
        Tint(1) = T_i
        Tint(2) = T_i + 1
        call co_DS_expr_dev ( pl(1:nsc), rho_c(i,j,k), Tint, Tdf_, pl(np), Dk(1:nsc) )
        km_ = 1
        do s = 2, nsc
          if ( Dk(s) > Dk(km_) ) km_ = s
        end do
        Dmax_ = Dk(km_) + mie/(rho_c(i,j,k)*sctd)
        do d = 1, 3
          an(13+d,i,j,k) = dl_f(d,i,j,k)**2 / Dmax_
        end do
        an(17,i,j,k) = real(km_, R8)
        an(18,i,j,k) = pl(km_) / rho_c(i,j,k)
        an(19,i,j,k) = Dmax_
      end if
    end do; end do; end do
    !$acc end data
  end subroutine Dt_Anatomy_blk
end module MOSE_Mod_dt
