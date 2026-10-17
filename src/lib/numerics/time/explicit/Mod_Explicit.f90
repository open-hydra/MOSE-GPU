! P is kept device-resident to avoid a per-kernel P transfer. Toggle off only to A/B it.
!
! CONSTRAINT — device kernels reading resident P must be called DIRECTLY, with exactly
! ONE dummy-argument level (as Newstate_* and Chemistry_Newstate_acc are), under
! `!$omp single`, with the matching `data copy(P)` / `present(P)` clause. Routing a
! kernel through an extra wrapper level hands nvfortran a host address it cannot match
! to the resident P present-table entry: the kernel then operates on a throwaway buffer
! and leaves the resident device P stale. This fails SILENTLY — no error, just wrong
! state — so it is not something the compiler will catch for you.
module MOSE_Mod_Explicit
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none

  logical, save :: field_check_done      = .false.  ! Phase-2 field-map self-test runs once per run

contains



  !> Phase-2 field-map self-test (one-time, env MOSE_FIELD_CHECK). Allocates the flat
  !> var-last field and verifies cidx is a bijection over the padded multi-block cells and
  !> unpack(pack(P))==P, then reports + error_stops on failure. Host-only, off by default
  !> (zero impact on the committed path). De-risks the SoA migration before any kernel
  !> consumes w (plan risk register: cidx/blk_off is the new single point of failure).
  subroutine field_selftest_once ( domain )
    use MOSE_Types_m,  only: MOSE_domain_type
    use MOSE_Field_m,  only: field_alloc, field_is_allocated, field_roundtrip_check
    type(MOSE_domain_type), intent(inout) :: domain
    character(len=8) :: val
    integer  :: ln, st
    real(R8) :: dev
    logical  :: ok
    if (field_check_done) return
    field_check_done = .true.
    call get_environment_variable ('MOSE_FIELD_CHECK', val, ln, st)
    if (.not. (st == 0 .and. ln > 0)) return
    if (.not. field_is_allocated()) call field_alloc(domain)
    call field_roundtrip_check (domain, dev, ok)
    write(*,'(a,l1,a,es10.3)') ' [FIELD-CHECK] flat var-last index map: bijection+identity ok=', &
      ok, '  roundtrip max|dev|=', dev
    if (.not. ok) error stop 'FIELD-CHECK: flat var-last index map is NOT a valid bijection/identity'
  end subroutine field_selftest_once


  subroutine Explicit_Step ( domain, config, MGL, level, External_Function )
    use MOSE_Types_m
    use MOSE_Global_m
    use MOSE_Mod_dt
    use MOSE_Lib_Strang, only: N_strang
    use MOSE_Lib_RK, only: obj_rk
    use MOSE_Lib_Extrapolation
    use MOSE_Mod_Fluxes
    use MOSE_Mod_BC
    use MOSE_Lib_RANS
    use MOSE_Lib_Newstate
    use MOSE_Lib_Omega_Relax, only: obj_omega_relax, Omega_Relax_Source, Omega_Relax_Floor
    use MOSE_Lib_Preconditioning, only: obj_prec, Update_Derived_Variables
    use MOSE_Lib_Chemistry
    use MOSE_Mod_Diagnostic
    use MOSE_IO_Solution
    use MOSE_Lib_Soot
    use MOSE_Lib_RotatingFrame, only: RotatingFrame_Source_Terms
    use MOSE_Mod_Thermo_Cache, only: Refresh_Thermo_Cache
    use MOSE_Mod_MPI, only: is_local_block, mpi_reduce_sum_r8_array, &
                           mpi_is_root, mpi_bcast_logical, mpi_bcast_integer, mpi_size_
    use MOSE_IO_RANS, only: obj_rans
    use MOSE_Field_m, only: field_is_allocated, field_alloc, field_sync_P_to_flat, &
                            field_update_P_flat_host, unpack_flat_to_P, &
                            field_update_PO_flat_host, unpack_flat_to_PO, &
                            field_add_R_to_R_flat

    implicit none
    integer, intent(in)                   :: MGL
    type(MOSE_domain_type), intent(inout) :: domain(MGL)
    type(config_type), intent(inout)      :: config
    integer, intent(inout)                :: level
    external :: External_Function
    ! Local
    logical  :: endsim, iosim, endmg, need_host_P, res_needs_host_P
    integer  :: i_rk, i_strang, b
    real(R8) :: average(nres)


    config % change_MG = .false.
    domain(level) % iter = domain(level) % iter + 1
    config % iter_from_call = config % iter_from_call + 1
    config % iter_general   = config % iter_general + 1

    ! Phase-2: one-time self-test of the flat var-last
    ! field index map, gated on env MOSE_FIELD_CHECK (host-only, off by default).
    call field_selftest_once ( domain(level) )
  
    if (config % HYDRA_time_accurate) then
      call Set_Global_dt ( domain(level) )
      domain(level) % time = domain(level) % time + domain(level) % dtglobal
    else
      domain(level) % dtglobal = 1d5
      call Compute_dt ( domain(level), config )  ! Compute local and minimum time step
      if ( config % time_accurate ) then
        call Set_Global_dt ( domain(level) )  ! Time-accurate: apply global minimum time step
        domain(level) % time = domain(level) % time + domain(level) % dtglobal 
      endif
    endif
     

    !$omp parallel
    strang: do i_strang = 1, N_strang ! Strang splitting: L = L_rk(dt/2) * L_chem(dt) * L_rk(dt/2)

      call Copy_State ( domain(level) )   ! P -> PO; device kernel when _OPENACC (no host P needed)

      rk: do i_rk = 1, obj_rk%n_RK

        ! The halo on w: the flat-w fill on one rank, the multi-rank driver (inter-rank exchange + local fill) on
        ! more. w interior is maintained by Compute_dt's one-time init + Newstate + chemistry.
        !$omp single
        if (.not. field_is_allocated()) call field_alloc(domain(level))
        !$omp end single
        call Fill_Ghost_Cell ( domain(level) )
        ! At substage (1,1) the interior cache
        ! bit-duplicates what Compute_dt just computed on the SAME w (nothing writes the
        ! w interior between them; the halo writes ghosts only) -> refresh ONLY the
        ! ghost shell. Guards: HYDRA_time_accurate skips Compute_dt entirely -> full
        ! refresh; strang-2 substage 1 refreshes fully too (chemistry rewrote w).
        call Refresh_Thermo_Cache ( domain(level), shell_only = ( i_strang == 1 .and.   &
               i_rk == 1 .and. .not. config % HYDRA_time_accurate ) )
        call Fluxes ( domain(level) )                        ! Convective and Diffusive fluxes

        call BC_Fluxes ( domain(level) )                     ! Boundary fluxes
        ! NB no post-BC P re-sync: the device BC kernels take P intent(in) (Lib_BC.f90,
        ! every BC_*_dev) and write only r — they cannot dirty P ghosts. Nothing between
        ! Fill_Ghost_Cell and the next P-ghost reader (host Extrapolate_Wall) writes P, so
        ! the ghosts filled by Fill_Ghost_Cell are still current. (Removed the /
        ! "restore/re-clean" — it was justified by the debunked premise that device BC
        ! writes P ghosts; Lib_BC:303,305 are the HOST BC_Symmetry_Visc, not the device path.)

        call External_Function ( domain(level) )             ! External function (e.g. source terms)
        
        if (ROTATING_FRAME) then
          ! This host source term accumulates into the device-
          ! resident r -> sync the device-accumulated r down, add on host, push back.
          ! It reads the state of the stage (Coriolis: the relative velocity, energy: w . centrifugal), which lives
          ! in the resident P_flat: P is refreshed only at the end of a step (not every step), so the source read the
          ! velocities of the start of the step in every stage (pressure-centrifugal: 1e-3 in v, w after one RK3
          ! step). The host P is filled from w first, as at the end of a step.
          call field_update_P_flat_host ()
          call unpack_flat_to_P ( domain(level) )
          do b = 1, domain(level) % nb
            if (is_local_block(b)) then
              !$acc update host(domain(level)%blk(b)%r)
            end if
          end do
          call RotatingFrame_Source_Terms ( domain(level) ) ! Coriolis and centrifugal source terms
          do b = 1, domain(level) % nb
            if (is_local_block(b)) then
              !$acc update device(domain(level)%blk(b)%r)
            end if
          end do
        end if

        if (RANS) then
          !$omp single
          call Extrapolate_Wall ( domain(level) )   ! the wall ghost of the turbulence variables, on w
          !$omp end single
          call RANS_Source_Terms ( domain(level) ) ! RHS of turbulence equations
        endif

        if (soot) then
          ! Host soot source accumulates into device-resident r.
          do b = 1, domain(level) % nb
            if (is_local_block(b)) then
              !$acc update host(domain(level)%blk(b)%r)
            end if
          end do
          call soot_source_terms ( domain(level) ) ! RHS of soot equations
          do b = 1, domain(level) % nb
            if (is_local_block(b)) then
              !$acc update device(domain(level)%blk(b)%r)
            end if
          end do
        end if

        ! Host sources (RotatingFrame/soot) accumulate into
        ! the var-first r (synced to the device by their brackets above); fold r into the
        ! var-LAST R_flat so Newstate reads R_flat ALONE. Same single addition the old bridge read
        ! (r + R_flat) performed per cell -> bit-exact. Skipped entirely when no host source
        ! is active (production THOR/RDE: r is never written, never read).
        if (ROTATING_FRAME .or. soot) then
          !$omp single
          call field_add_R_to_R_flat ( domain(level) )
          !$omp end single
        end if
        ! omega-source-relax: the point-implicit SST destruction on the resident residual R_flat
        if ( obj_omega_relax%active ) then
          !$omp single
          call Omega_Relax_Source ( domain(level) )
          !$omp end single
        end if
        call RK_Newstate ( domain(level), i_rk )             ! State update
        if ( obj_omega_relax%active ) then   ! omega >= rho * omfloor on the updated state
          !$omp single
          call Omega_Relax_Floor ( domain(level) )
          !$omp end single
        end if
        if ( obj_prec%enabled ) then   ! Ur of the updated state (MOSE_open Update_Derived_Variables)
          !$omp single
          call Update_Derived_Variables ( domain(level), .false. )
          !$omp end single
        end if

      enddo rk

      ! Chemistry, only if i_strang=1 and N_strang=2
      if ( i_strang + N_strang == 3 ) then
        ! React on the resident var-LAST w in place. NO unpack P_flat -> P
        ! -- P is the resident truth only at the end-of-step IO/residual boundary now; the
        ! next Strang substage (Copy_State w->PO_flat, halo/flux/BC/cache/Newstate/SST on w) never
        ! reads P interior. Removes the per-step chemistry P_flat -> P host transpose.
        !$omp single
        do b = 1, domain(level) % nb
          if (.not. is_local_block(b)) cycle
          if (domain(level)%blk(b)%no_chem) cycle
          call Chemistry_Newstate_blk ( b, domain(level) % blk(b) % dtlocal, domain(level) % blk(b) % dim, domain(level) % blk(b) % hr_c )
        enddo
        !$omp end single
      endif

    enddo strang
    !$omp end parallel
    ! OVERHEAD FIX: host P is read ONLY by Compute_Residual (periodic) and I/O writes, so
    ! gate the whole-field D2H+transpose to those iters (was unconditional, ~40% of loop).
    ! residuotot-based termination only fires on a residual iter (residuotot fresh) -> covered.
    res_needs_host_P = ( mod (domain(level) % iter, obj_io%res_diter) == 0 )
    ! Device residual (unconditional): runs on the device (resident P_flat/PO_flat, Compute_Residual_blk)
    ! for the nres=7 RANS layout, so a residual-only iteration no longer forces the
    ! whole-field w D2H + host P_flat -> P scatter. Probes are refreshed point-wise from w
    ! instead (Mod_Probes). Other layouts (nres/=7) keep the host path.
    if (nres == 7) res_needs_host_P = .false.
    need_host_P =      ( config % iter_from_call == 1 ) &
                  .or. res_needs_host_P &
                  .or. ( mod (domain(level) % iter, obj_io%sol_diter) == 0 ) &
                  .or. ( mod (domain(level) % iter, obj_io%bck_diter) == 0 ) &
                  .or. ( domain(level) % time >= config % time_from_call + obj_io%sol_dtime ) &
                  .or. ( config % iter_from_call >= domain(level) % itermax ) &
                  .or. ( domain(level) % time >= config % timemax )
    ! w is the resident truth; P is not maintained in the loop. Refresh the host
    ! P (D2H w + scatter w->host P) for Compute_Residual + I/O, ONLY when needed.
    if (need_host_P) then
      call field_update_P_flat_host ()
      call unpack_flat_to_P ( domain(level) )
    end if

    ! Save residuals
    if ( (config % iter_from_call == 1) .or. (mod (domain(level) % iter, obj_io%res_diter) == 0) ) then
      config % residuotot = 0.d0
      if (nres == 7) then
      ! Device residual: on the resident P_flat/PO_flat (device reduction) -- no PO_flat D2H, no
      ! host PO scatter, no dtlocal sync (Compute_Residual never reads its dt argument).
      ! Same interior cells, same algebra; only the accumulation order differs (device
      ! reduction tree vs host loop) -> trailing-digit shifts in residual-history.dat;
      ! the FIELD is untouched (residual is diagnostic + restop termination only).
      do b = 1, domain(level) % nb
        if (.not. is_local_block(b)) cycle
        call Compute_Residual_blk ( b, domain(level)%blk(b)%dim, config%residuotot )
      enddo
      else
      ! PO is no longer maintained on device (Copy_State does w->PO_flat). Refresh the
      ! host PO from the resident P_flato (D2H PO_flat + scatter PO_flat->host PO); dtlocal is resident.
      call field_update_PO_flat_host ()
      call unpack_flat_to_PO ( domain(level) )
      do b = 1, domain(level) % nb
        if (is_local_block(b)) then
          !$acc update host(domain(level)%blk(b)%dtlocal)
        end if
      end do
      do b = 1, domain(level) % nb
        if (.not. is_local_block(b)) cycle
        call Compute_Residual ( new=domain(level)%blk(b)%p, &
                                old=domain(level)%blk(b)%po, &
                                dt=domain(level)%blk(b)%dtlocal, &
                                n=domain(level)%blk(b)%dim, &
                                average=average, &
                                total=config%residuotot )
      enddo
      end if
      call mpi_reduce_sum_r8_array(config%residuotot, nres)
      if (mpi_is_root) config % residuotot = sqrt ( config % residuotot ) ! L2 norm time derivative
    endif

    ! Determine simulation control flags on root, then broadcast to all ranks
    if (mpi_is_root) then
      if (level == 1) then
        endsim = ( config % iter_from_call >= domain(1) % itermax ) &
            .or. ( config % residuotot(1) <= config % restop ) &
            .or. ( domain(1) % time >= config % timemax )
      else
        endmg = ( config % iter_from_call >= domain(level) % itermax ) &
           .or. ( config % residuotot(1) <= config % restop ) &
           .or. ( domain(level) % time >= config % timemax )
      endif

      iosim  = ( mod (domain(level) % iter, obj_io%sol_diter) == 0) &
          .or. ( mod (domain(level) % iter, obj_io%bck_diter) == 0) &
          .or. ( domain(level) % time >= config % time_from_call + obj_io%sol_dtime )

      if ( endsim ) then
        config % TODO = 3
      elseif ( iosim .or. endmg ) then
        config % TODO = 2
        if (endmg .and. config % MGL > 1) then
          config % change_MG = .true.
        endif
      else
        config % TODO = 1
      endif
    end if

    ! Broadcast simulation control from root to all ranks
    call mpi_bcast_integer(config % TODO)
    call mpi_bcast_logical(config % change_MG)
      
  end subroutine Explicit_Step


  !> Device residency of one domain, once at setup: the flat state w allocated and seeded from P, then the geometry, the
  !> work arrays, the thermo caches, P and the boundary / connection tables mapped to the device (the kernels' copyin /
  !> copy clauses are no-ops from here on, present_or semantics). The chimera / manifold flat tables need the field arena,
  !> so they are built here too. Static mesh, one multigrid level.
  subroutine Setup_Device_Residency ( domain )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Mod_BC,  only: Setup_BC_Device_Flat
    use MOSE_Field_m, only: field_is_allocated, field_alloc, pack_P_to_flat, field_update_P_flat_device
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: b

    if (.not. field_is_allocated()) call field_alloc ( domain )
    call pack_P_to_flat ( domain )
    call field_update_P_flat_device ()
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      ! Geometry: written once at setup, read-only thereafter.
      !$acc enter data copyin(domain%blk(b)%fm1, domain%blk(b)%fm2,  &
      !$acc&                  domain%blk(b)%fm3, domain%blk(b)%dl_f, &
      !$acc&                  domain%blk(b)%m_f, domain%blk(b)%yn,   &
      !$acc&                  domain%blk(b)%vol)
      ! r (residual), beta (shock sensor), PO (RK base state), dtlocal (per-cell dt):
      ! written and read only by device kernels in the loop — r by Fluxes/BC/SST then
      ! Newstate; beta by Fluxes then BC; PO by Copy_State then Newstate; dtlocal by
      ! Compute_dt/Set_Global_dt then Newstate/chemistry. Their only host reader is
      ! I/O (Write_Diagnostic), which does its own `acc update host`.
      !
      ! r must be zeroed here: on the W path it is only a bridge for host-side sources
      ! (RotatingFrame/soot), so with no such source it is never written — yet
      ! Write_Diagnostic still reads it back and adds R_flat on the host. Not every
      ! allocation path zeroes it (BlockSplit does, others may not), so without this
      ! that read is undefined.
      domain%blk(b)%r = 0.0_R8
      !$acc enter data copyin(domain%blk(b)%r,    &
      !$acc&                  domain%blk(b)%beta,  &
      !$acc&                  domain%blk(b)%PO,    &
      !$acc&                  domain%blk(b)%dtlocal)
      ! Thermo caches (rho_c/T_c/Rg_c/mil_c): written only by the device thermo_cache_*
      ! kernels, read only by device kernels. No host stage touches them — output T is
      ! recomputed from P in IO_Solution — so they need no host sync at all.
      ! CAVEAT: that holds for SST. Other RANS models (SSGLRR/Spalart) read these
      ! caches on the host and would need syncs before this is safe for them.
      ! (kl_c/cp_c are NOT made resident: face-exact recomputes conductivity/cp at
      ! faces, so nothing on the device reads those caches.)
      !$acc enter data copyin(domain%blk(b)%rho_c, domain%blk(b)%T_c, &
      !$acc&                  domain%blk(b)%Rg_c,  domain%blk(b)%mil_c, &
      !$acc&                  domain%blk(b)%hr_c)
      ! P: see the residency constraint at the top of this file. The host syncs
      ! (f5_P_to_host/device) around Copy_State / Fill_Ghost_Cell are still required
      ! so the host stages see current P.
      !$acc enter data copyin(domain%blk(b)%P)
    end do
    ! Domain-level resident BC tables. These are STATIC device dispatch/data arrays
    ! sized by nbound, re-copied on EVERY device BC/ghost call before. At THOR scale
    ! (nbnd~1.8M) that was hundreds of GB/iter of pointless H2D (the GPU stalled ~40%
    ! util): bc_pg(nprim,6,nbnd)~850 MB was re-copied per CONNECTION PAIR per substage
    ! (the dominant cost), and bc_bi..bc_tvec per BLOCK per substage. Resident -> the
    ! present_or copyin clauses in BC_Fluxes_blk / Fill_Ghost_Cell_*_acc / the
    ! connection kernel all become no-ops. bc_p0 (time-varying inflow) and bc_pg
    ! (non dev_ghost host refresh) are kept current with `update device` in Mod_BC.
    if (allocated(domain%bc_bi)) then
      !$acc enter data copyin(domain%bc_bi, domain%bc_bj,           &
      !$acc&  domain%bc_bk, domain%bc_bf, domain%bc_typ,      &
      !$acc&  domain%bc_mach, domain%bc_t0, domain%bc_p0,     &
      !$acc&  domain%bc_relfac, domain%bc_alpha, domain%bc_beta, &
      !$acc&  domain%bc_pamb, domain%bc_ci, domain%bc_rans,   &
      !$acc&  domain%bc_tvec, domain%bc_tw, domain%bc_pg)
    end if
    if (domain%n_conn > 0) then
      !$acc enter data copyin(domain%conn_im, domain%conn_jm,   &
      !$acc&  domain%conn_km, domain%conn_fm, domain%conn_is, &
      !$acc&  domain%conn_js, domain%conn_ks, domain%conn_fs, &
      !$acc&  domain%conn_d11, domain%conn_d12, domain%conn_d21, &
      !$acc&  domain%conn_d22, domain%conn_pg, domain%conn_all, domain%conn_all_bs)
    end if
    ! Chimera (666) / manifold (1001) flat tables: global-flat-index based, so they
    ! need the field arena; builds + enter-data's its own arrays (no-op without
    ! such faces). Placed here so every later Fill_Ghost_Cell_local / BC kernel
    ! sees them resident.
    call Setup_BC_Device_Flat ( domain )
  end subroutine Setup_Device_Residency


  subroutine Copy_State ( domain )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Field_m, only: field_is_allocated, field_copy_P_flat_to_PO_flat

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: i, j, k, b

    ! PO resident. Copy P->PO on device (one dummy level, like the
    ! other _acc kernels, so the resident PO present-matches). copyin(P) is a no-op
    ! when P is resident. One
    ! `!$omp single` around the whole block loop (like Newstate/chemistry dispatch).
    !$omp single
    ! Snapshot the var-last RK base on device (PO_flat <- w). P is no longer maintained
    ! during the substage loop (w is the resident truth), so the P->PO copy is DROPPED here;
    ! the host PO is refreshed only end-of-step (from PO_flat) for Compute_Residual. w is current
    ! at the substage start (previous Newstate/chemistry write; Compute_dt's one-time init).
    if (field_is_allocated()) call field_copy_P_flat_to_PO_flat ()
    !$omp end single

  end subroutine Copy_State




end module MOSE_Mod_Explicit
