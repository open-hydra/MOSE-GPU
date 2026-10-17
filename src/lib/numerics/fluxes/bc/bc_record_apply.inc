! The BC applied to one record, with every value already in a scalar: shared verbatim by the
! per-block kernels (which derive them from their own arrays) and by the fused kernel (which reads
! them from the per-level gather), so the two paths cannot drift.
        select case (bc_typ(idx))
          case (1, 666)   ! connection & chimera (donor data already in the receiver ghosts / bc_pg;
                          ! chimera ghosts filled by Ghost_Chimera -- host BC_Fluxes routes 666
                          ! through the connection twins identically)
            call BC_Connection_Eul ( PWARG, r, beta_, dlc4_, nb1_, nb2_, nb3_, &
                                         Im, Jm, Km, Fm, normal, area )
            if (VISC) &
              call BC_Connection_Visc ( PWARG, r, Mrec_, Mgho_, &
                                            ynrec_, yngho_, bc_pg(:,:,idx), &
                                            nb1_, nb2_, nb3_, Im, Jm, Km, Fm, normal, area )
          case (2)   ! rotational periodicity (wedge seam, faces 5-6; reads block-own
                     ! opposite-k cells + rotates donor (v,w) by delthe — no ghost)
            call BC_Rotational_Periodic_Eul ( PWARG, r, dl_f, nb1_, nb2_, nb3_, &
                                                  Im, Jm, Km, Fm, normal, area, delthe )
          case (3)   ! symmetry
            call BC_Symmetry_Eul ( PWARG, r, nb1_, nb2_, nb3_, Im, Jm, Km, Fm, &
                                       normal, area, Bound_Dl, Int_Dl )
            if (VISC) &
              call BC_Symmetry_Visc ( PWARG, r, Mrec_, ynrec_, &
                                          nb1_, nb2_, nb3_, Im, Jm, Km, Fm, normal, area )
          case (5)   ! no-slip wall = Euler pressure (mirrors case 3) + viscous wall stress
            call BC_Symmetry_Eul ( PWARG, r, nb1_, nb2_, nb3_, Im, Jm, Km, Fm, &
                                       normal, area, Bound_Dl, Int_Dl )
            if (VISC) &
              call BC_Wall_Heat_dev ( PWARG, r, Mrec_, ynrec_, &
                                      nb1_, nb2_, nb3_, Im, Jm, Km, Fm, normal, area )
          case (6)   ! isothermal no-slip wall: Euler pressure + viscous stress + conduction
            call BC_Symmetry_Eul ( PWARG, r, nb1_, nb2_, nb3_, Im, Jm, Km, Fm, &
                                       normal, area, Bound_Dl, Int_Dl )
            if (VISC) &
              call BC_Wall_Temperature_dev ( PWARG, r, Mrec_, ynrec_, &
                                             nb1_, nb2_, nb3_, Im, Jm, Km, Fm, normal, area, bc_tw(idx) )
          case (4)   ! inflow/outflow. Falls back to symmetry on closed (outflow) inlet.
            call BC_Inflow_Outflow ( PWARG, r, nb1_, nb2_, nb3_, Im, Jm, Km, Fm, normal, area, &
                                         bc_tvec(:,idx), bc_mach(idx), bc_t0(idx), bc_p0(idx), &
                                         bc_relfac(idx), bc_alpha(idx), bc_beta(idx), bc_pamb(idx), &
                                         bc_ci(:,idx), bc_rans(1:max(nrans,1),idx), ierr )
            if (ierr == 1) &
              call BC_Symmetry_Eul ( PWARG, r, nb1_, nb2_, nb3_, Im, Jm, Km, Fm, &
                                         normal, area, Bound_Dl, Int_Dl )
          case (22)  ! nozzle injector (real device twin; pressure-regime selector).
            call BC_Inflow_Nozzle ( PWARG, r, nb1_, nb2_, nb3_, Im, Jm, Km, Fm, normal, area, &
                                        bc_tvec(:,idx), bc_mach(idx), bc_t0(idx), bc_p0(idx), &
                                        bc_relfac(idx), bc_alpha(idx), bc_beta(idx), bc_pamb(idx), &
                                        bc_ci(:,idx), bc_rans(1:max(nrans,1),idx), ierr )
            if (ierr == 1) &
              call BC_Symmetry_Eul ( PWARG, r, nb1_, nb2_, nb3_, Im, Jm, Km, Fm, &
                                         normal, area, Bound_Dl, Int_Dl )
          case (24)  ! back-pressure outlet (MOSE_open 406): symmetry flux on reversed flow, as the host
            call BC_Outlet ( PWARG, r, nb1_, nb2_, nb3_, Im, Jm, Km, Fm, normal, area, &
                                 bc_pamb(idx), bc_relfac(idx), ierr )
            if (ierr == 1) &
              call BC_Symmetry_Eul ( PWARG, r, nb1_, nb2_, nb3_, Im, Jm, Km, Fm, &
                                         normal, area, Bound_Dl, Int_Dl )
          case (11)  ! extrapolation
            call BC_Extrapolation ( PWARG, r, nb1_, nb2_, nb3_, Im, Jm, Km, Fm, normal, area )
          case (1001) ! manifold mean connection: mass-flux inflow (mach=-10) with the
                      ! per-substage (T0, g, pstat) scalars written into bc_t0/bc_p0/bc_pamb
                      ! by BC_Manifold_Update; static params preset in Setup_BC_Device.
                      ! Mirrors the host BC_Fluxes case(1001) -> BC_Inflow_Outflow call.
            call BC_Inflow_Outflow ( PWARG, r, nb1_, nb2_, nb3_, Im, Jm, Km, Fm, normal, area, &
                                         bc_tvec(:,idx), bc_mach(idx), bc_t0(idx), bc_p0(idx), &
                                         bc_relfac(idx), bc_alpha(idx), bc_beta(idx), bc_pamb(idx), &
                                         bc_ci(:,idx), bc_rans(1:max(nrans,1),idx), ierr )
            ! NB: host case(1001) has NO symmetry fallback on ierr==1 (adds no flux) -- mirror that.
        end select
