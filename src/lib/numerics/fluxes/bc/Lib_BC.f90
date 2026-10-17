module MOSE_Lib_BC
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Global_m
  use MOSE_Types_m
  use MOSE_Parameters_m

  implicit none

! The BC routines read the state from the flat field (P_flat, dummy `w` here) and write the residual into
! R_flat (dummy `R_flat`): PCELL(s,i,j,k) is the read, WRCELL(s,i,j,k) the write, at the flat row
! boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1 of the block. PWARG threads the two arrays and the
! block's offsets through every signature and call site: a module allocatable cannot be referenced inside an
! `acc routine seq` (CUDA illegal address), so they travel as explicit-shape arguments.
#  define PWARG w, wr, ncw, nrw, boff, pi, pj, wroff, ur, nur
#  define URCELL(i,j,k) ur( boff + ((i)+gc-1) + ((j)+gc-1)*pi + ((k)+gc-1)*pi*pj + 1 )
#  define PCELL(s,i,j,k) w( boff + ((i)+gc-1) + ((j)+gc-1)*pi + ((k)+gc-1)*pi*pj + 1, s )
#  define WRCELL(s,i,j,k) wr( wroff + ((i)+gc-1) + ((j)+gc-1)*pi + ((k)+gc-1)*pi*pj + 1, s )

contains


  subroutine Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )

    implicit none
    integer, intent(in) :: fm
    integer, intent(out) :: modfm, modfm1, modfm2, modfm3

    modFm  = mod(Fm,2)
    modfm1 = 1-modFm  
    modfm2 = 1-2*modFm
    modfm3 = 2*( 2*modFm-1 )

  end subroutine


  subroutine Face_Index ( face, dir, cell_i, cell_j, cell_k, face_i, face_j, face_k )
    
    implicit none
    integer, intent(in) :: face, cell_i, cell_j, cell_k
    integer, intent(out) :: dir, face_i, face_j, face_k

    select case ( face )
    
      case(1:2)
        dir = 1
        face_i = cell_i - mod ( face, 2 )
        face_j = cell_j
        face_k = cell_k
      
      case(3:4)
        dir = 2
        face_i = cell_i
        face_j = cell_j - mod ( face, 2 )
        face_k = cell_k
      
      case(5:6)
        dir = 3
        face_i = cell_i
        face_j = cell_j
        face_k = cell_k - mod ( face, 2 )
  
    end select

  end subroutine Face_Index








  subroutine BC_Wall_Heat ( Im, Jm, Km, Fm, Blk, Heat_Flux, Ovar, w_wall )
    use MOSE_Mod_Fluid
    use MOSE_Lib_RANS
    use FLINT_Lib_Thermodynamic

    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    real(R8), intent(in) :: Heat_Flux
    type(MOSE_block_type), intent(inout) :: Blk
    real(R8), optional, dimension(8), intent(inout) :: Ovar
    real(R8), optional, dimension(3), intent(in) :: w_wall  !> Wall velocity in rotating frame [m/s]
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Face_k, Ig, Jg, Kg
    real(R8) :: Normal(3), Area, Dist, M(3,3), Prim(nprim), rho, Rgas, temp, Gradient(nprim,3)
    real(R8) :: mil, kl, Stress(3), Flux(nprim), Prim_Wall(nprim), rho_wall, Twall, dl


    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a
    Dist = 1d-20 ! since the flux is computed at the wall
    M = Blk % M (Im,Jm,Km) % c

    ! boundary cell variables
    Prim = Blk % P(:,Im,Jm,Km)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    temp = EOS( p=prim(np), ro=rho, R=Rgas )

    ! Approximation: Twall==Tcell. Therefore, mil(wall) = mil(cell) and mit = 0
    call co_k_mi_lam_Wilke ( Prim(1:nsc), rho, temp, mil, kl )

    ! Initialization
    Gradient = 0d0

    ! Difference in Across-face direction: symmetry; tangential direction: null
    ! w_wall: velocity of the wall in the rotating frame (0 for rotating blades, -Omega x r for stationary walls)
    if ( present(w_wall) ) then
      Gradient(nu:nw,Dir) = ( Prim(nu:nw) - w_wall ) * modfm3  ! Velocity relative to wall
    else
      Gradient(nu:nw,Dir) = Prim(nu:nw) * modfm3               ! Standard no-slip (w_wall=0)
    end if
    Gradient(np,Dir) = -Heat_Flux / ( 2d0*kl*dot_product( M(Dir,:), normal ) ) * modfm2

    ! Wall variables
    Twall = temp - Gradient(np,Dir)
    rho_wall = Prim(np)/( Rgas * Twall ) ! EOS
    Prim_wall(1:nsc) = Prim(1:nsc) * rho_wall / rho

    if (RANS)  then
      dl = blk % yn(im,jm,km)
      call RANS_Set_Wall_Values( mil, Prim_Wall(nt:nprim), dl )
      Gradient(nt:nprim,Dir) = ( Prim(nt:nprim)/rho - Prim_Wall(nt:nprim)/rho_wall ) * modfm3 
    end if
    
    Gradient = matmul ( Gradient, M )

    ! Stress vector
    Stress = Stress_Vector ( Gradient(nu:nw,:), Normal, mil, 0d0, Prim(nt:) )

    ! Fluxes
    Flux = 0.0d0
    Flux(nu:nw) = Stress * Area
    Flux(np) = Area * Heat_Flux

    if (RANS) then
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), &
                                 rans_variables=Prim_Wall(nt:nprim), &
                                 vel_gradient=Gradient(nu:nw,:), &
                                 rans_gradient=Gradient(nt:nprim,:), &
                                 mul=mil, rho=rho_wall, &
                                 area=area, normal=normal, dist=dist )
      ! Ghost-cell extrapolation of RANS variables
      ig = im - guide(fm,1)
      jg = jm - guide(fm,2)
      kg = km - guide(fm,3)
      call RANS_Extrapolate_Wall ( Prim(nt:nprim), Prim_Wall(nt:nprim), &
                                   rho, rho_wall, Blk % P(nt:nprim,Ig,Jg,Kg) )
    endif

    ! Residual update
    Blk % r (:,Im,Jm,Km) = Blk % r (:,Im,Jm,Km) - modfm2 * Flux

    if (present(Ovar)) call Compute_Wall_Properties(stress=Stress, pw=prim(np), qw=Heat_Flux, &
                                                    y=blk%dl(im,jm,km)%c(dir)*0.5d0,          &
                                                    Tw=Twall, rhow=rho_wall, mu=mil, exit_array=Ovar)

  end subroutine BC_Wall_Heat


  subroutine BC_Wall_Temperature ( Im, Jm, Km, Fm, Blk, T_wall, Ovar, w_wall )
    use MOSE_Mod_Fluid
    use MOSE_Lib_RANS
    use FLINT_Lib_Thermodynamic

    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    real(R8), intent(in) :: T_wall
    type(MOSE_block_type), intent(inout) :: Blk
    real(R8), optional, dimension(8), intent(out) :: Ovar
    real(R8), optional, dimension(3), intent(in) :: w_wall  !> Wall velocity in rotating frame [m/s]
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Face_k, ig, jg, kg
    real(R8) :: Normal(3), Area, Dist, M(3,3), Prim(nprim), rho, Rgas, T, rho_wall, dl
    real(R8) :: Prim_wall(nprim), mil, kl, Gradient(nprim,3), Stress(3), Flux(nprim)

    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a
    Dist = 1d-20 ! since the flux is computed at the wall
    M = Blk % M (Im,Jm,Km) % c

    ! boundary cell variables
    Prim = Blk % P(:,Im,Jm,Km)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    T = Prim(np) / ( rho * Rgas )

    ! Initialization
    Gradient = 0d0

    ! Boundary layer: p(wall) = p(cell). Also, non catalytic wall, so R(wall) = R(cell) and ci(wall) = ci(cell), but not roi
    rho_wall = Prim(np) / ( Rgas * T_wall )
    Prim_wall(1:nsc) = Prim(1:nsc) * rho_wall / rho

    ! Laminar viscosity/conductivity at interface computed at Twall
    call co_k_mi_lam_Wilke ( Prim_wall(1:nsc), rho_wall, T_wall, mil, kl )

    ! Difference in Across-face direction: symmetry; tangential direction: null
    ! w_wall: velocity of the wall in the rotating frame (0 for rotating blades, -Omega x r for stationary walls)
    if ( present(w_wall) ) then
      Gradient(nu:nw,Dir) = ( Prim(nu:nw) - w_wall ) * modfm3  ! Velocity relative to wall
    else
      Gradient(nu:nw,Dir) = Prim(nu:nw) * modfm3               ! Standard no-slip (w_wall=0)
    end if
    Gradient(np,Dir) = ( T - T_wall ) * modfm3                  ! Temperature

    if (RANS)  then
      dl = blk % yn(im,jm,km)
      call RANS_Set_Wall_Values( mil, Prim_Wall(nt:nprim), dl )
      Gradient(nt:nprim,Dir) = ( Prim(nt:nprim)/rho - Prim_Wall(nt:nprim)/rho_wall ) * modfm3 
    end if

    ! Transformation
    Gradient = matmul ( Gradient, M )

    ! Stress vector
    Stress = Stress_Vector ( Gradient(nu:nw,:), Normal, mil, 0d0, Prim(nt:) )

    ! Fluxes
    Flux = 0.0d0 
    Flux(nu:nw) = Stress * Area
    Flux(np) = Area * kl * dot_product ( Gradient(np,:), Normal )

    if (RANS) then
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), &
                                 rans_variables=Prim_Wall(nt:nprim), &
                                 vel_gradient=Gradient(nu:nw,:), &
                                 rans_gradient=Gradient(nt:nprim,:), &
                                 mul=mil, rho=rho_wall, &
                                 area=area, normal=normal, dist=dist )
      ! Ghost-cell extrapolation of RANS variables
      ig = im - guide(fm,1)
      jg = jm - guide(fm,2)
      kg = km - guide(fm,3)
      call RANS_Extrapolate_Wall ( Prim(nt:nprim), Prim_Wall(nt:nprim), &
                                   rho, rho_wall, Blk % P(nt:nprim,Ig,Jg,Kg) )
    endif

    ! Residual update
    Blk % r (:,Im,Jm,Km) = Blk % r (:,Im,Jm,Km) - modfm2 * Flux

    if (present(Ovar)) call Compute_Wall_Properties(stress=Stress, pw=prim(np), Tw=T_wall,          &
                                                    rhow=rho_wall, mu=mil,                          &
                                                    y=blk%dl(im,jm,km)%c(dir)*0.5d0,                &
                                                    qw=kl * dot_product ( Gradient(np,:), Normal ), &
                                                    exit_array=Ovar)

  end subroutine BC_Wall_Temperature



  subroutine BC_Wall_Melting ( Im, Jm, Km, Fm, Blk, T_wall, q_rad, Ovar )
    use MOSE_Mod_Fluid
    use MOSE_Mod_GSI
    use MOSE_Lib_RANS
    use FLINT_Lib_Thermodynamic

    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    real(R8), intent(in) :: T_wall, q_rad
    type(MOSE_block_type), intent(inout) :: Blk
    real(R8), optional, dimension(8), intent(inout) :: Ovar
    ! Local
    integer :: s, modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Face_k, Ig, Jg, Kg
    real(R8) :: Normal(3), Area, Dist, dl, M(3,3), Prim(nprim), rho, Rgas, T, Gradient(nprim,3)
    real(R8) :: rho_wall, Prim_Wall(nprim), mil, kl, q_conv, omega(nsc), mdot
    real(R8) :: Stress(3), Blowing(3), Sum_Omegai_Hwi, Flux(nprim)


    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a
    Dist = 1d-20 ! since the flux is computed at the wall
    M = Blk % M (Im,Jm,Km) % c

    ! boundary cell variables
    Prim = Blk % P(:,Im,Jm,Km)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    T = Prim(np) / ( rho * Rgas )

    ! Initialization
    Gradient = 0d0

    ! Difference in Across-face direction: symmetry; tangential direction: null
    Gradient(np,Dir) = ( T - T_wall ) * modfm3                  ! Temperature

    Gradient(np,:) = matmul ( Gradient(np,:), M )               ! Temperature gradient
    rho_wall = Prim(np) / ( Rgas * T_wall )                     ! Boundary layer: p(wall) = p(cell), R(wall) = R(cell) and ci(wall) = ci(cell), but not roi
    Prim_wall(1:nsc) = Prim(1:nsc) * rho_wall / rho             ! Partial densities at wall

    call co_k_mi_lam_Wilke ( Prim_wall(1:nsc), rho_wall, T_wall, mil, kl ) ! transport

    if (RANS) then ! turbulence variables
      dl = blk % yn(im,jm,km)
      call RANS_Set_Wall_Values( mil, Prim_Wall(nt:nprim), dl )
      Gradient(nt:nprim,Dir) = ( Prim(nt:nprim)/rho - Prim_Wall(nt:nprim)/rho_wall ) * modfm3
    end if

    q_conv = kl * dot_product( Gradient(np,:), Normal)          ! Convective heat flux

    call GSI_Melting ( q_conv + q_rad, T_wall, omega, mdot )    ! Gas-surface interaction for paraffin mass flux

    Blowing = mdot / rho_wall * Normal                          ! Blowing velocity
    Gradient(nu:nw,Dir) = ( Prim(nu:nw) - Blowing ) * modfm3    ! Velocity gradient due to blowing in direction normal to the face

    Gradient(nu:nw,:) = matmul ( Gradient(nu:nw,:), M )         ! Velocity gradient transformation

    Stress = Stress_Vector ( Gradient(nu:nw,:), Normal, mil, 0d0, Prim(nt:) )

    Sum_Omegai_Hwi = 0.d0
    do s = 1, nsc
      Sum_Omegai_Hwi = Sum_Omegai_Hwi + omega(s) * f_tabT( T_wall,s,h_tab )
    end do

    ! Fluxes
    Flux = 0.0d0
    Flux(1:nsc) = - omega * Area
    Flux(nu:nw) = ( Stress - mdot * Blowing ) * Area
    Flux(np) = Area * ( q_conv + q_rad - Sum_Omegai_Hwi - mdot * 0.5d0 * sum( Blowing**2 ) + &
                        dot_product( Stress, Blowing ) ) ! check sign q_rad!!!

    if (RANS) then
      Gradient(nt:nprim,:) = matmul ( Gradient(nt:nprim,:), M )
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), &
                                 rans_variables=Prim_Wall(nt:nprim), &
                                 vel_gradient=Gradient(nu:nw,:), &
                                 rans_gradient=Gradient(nt:nprim,:), &
                                 mul=mil, rho=rho_wall, area=area, &
                                 normal=normal, dist=dist )
      ! Ghost-cell extrapolation of RANS variables
      ig = im - guide(fm,1)
      jg = jm - guide(fm,2)
      kg = km - guide(fm,3)
      call RANS_Extrapolate_Wall ( Prim(nt:nprim), Prim_Wall(nt:nprim), &
                                   rho, rho_wall, Blk % P(nt:nprim,Ig,Jg,Kg) )
    endif

    ! Residual update
    Blk % r (:,Im,Jm,Km) = Blk % r (:,Im,Jm,Km) - modfm2 * Flux

    if (present(Ovar)) call Compute_Wall_Properties(stress=Stress, pw=prim(np), qw=q_conv+q_rad, &
                                                    mdot=mdot, y=blk%dl(im,jm,km)%c(dir)*0.5d0, &
                                                    Tw=T_wall, rhow=rho_wall, mu=mil, exit_array=Ovar)

  end subroutine BC_Wall_Melting


  subroutine BC_Wall_Pyrolysis ( Im, Jm, Km, Fm, Blk, qrad, fuel, Ovar )
    use MOSE_Mod_Fluid
    use MOSE_Mod_GSI
    use MOSE_Lib_RANS
    use FLINT_Lib_Thermodynamic

    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm, fuel
    real(R8), intent(in) :: qrad
    type(MOSE_block_type), intent(inout) :: Blk
    real(R8), optional, dimension(8), intent(inout) :: Ovar
    ! Local
    integer :: s, modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Face_k, Ig, Jg, Kg, iter
    real(R8) :: Normal(3), Area, Dist, dl, M(3,3), Prim(nprim), rho, Rgas, T, Gradient(nprim,3)
    real(R8) :: rho_wall, Prim_Wall(nprim), mil, kl, omega(nsc), mdot
    real(R8) :: Stress(3), Blowing(3), Sum_Omegai_Hwi, Flux(nprim)
    real(R8) :: T_wall, T_wall_old, balance, balanceold, q_conv, q_pyro, q_rad


    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a
    M = Blk % M (Im,Jm,Km) % c

    ! boundary cell variables
    Prim = Blk % P(:,Im,Jm,Km)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    T = Prim(np) / ( rho * Rgas )

    ! Initialization
    Gradient = 0d0

    ! Newton-Raphson loop initialization
    T_wall = T ! First guess: T_wall=T_cell
    rho_wall = Prim(np) / ( Rgas * T_wall )         ! BL: p(wall)=p(cell), R(w)=R(c) and ci(w)=ci(c)
    Prim_wall(1:nsc) = Prim(1:nsc) * rho_wall / rho ! Partial densities at wall

    select case ( fuel )
      case ( 9 ) ! HDPE
        call GSI_HDPE(T_wall, omega, mdot)
      case ( 10 ) ! PP
        call GSI_PP(T_wall, omega, mdot)
      case ( 12 ) ! Carbon
        call GSI_carbon(Prim_wall(1:nsc), T_wall, omega, mdot)
      case ( 13 ) ! HTPB
        call GSI_HTPB(T_wall, omega, mdot)
    end select
    mdot = mdot * (-modfm2) ! adjust sign
    omega = omega * (-modfm2) ! adjust sign

    call co_k_mi_lam_Wilke ( Prim_wall(1:nsc), rho_wall, T_wall, mil, kl ) ! transport

    q_conv = 0d0 ! convective heat flux
    select case ( fuel ) ! pyrolysis heat flux
    case ( 9 ) ! HDPE
      q_pyro = mdot * ( Dh_HDPE + cp_HDPE* (T_wall - Ti_HDPE) ) 
    case ( 10 ) ! PP
      q_pyro = mdot * ( Dh_PP + cp_PP* (T_wall - Ti_PP) )
    case ( 12 ) ! Carbon
      q_pyro = 0d0
      do s = 1, nsc
        q_pyro = q_pyro + omega(s)*f_tabT( T_wall,s,h_tab )
      end do    
    case ( 13 ) ! HTPB
      q_pyro = mdot * ( Dh_HTPB + cp_HTPB* (T_wall - Ti_HTPB) )
    end select
    q_rad = qrad * modfm2 ! adjust sign for radiation
    balance = q_conv - q_pyro - q_rad

    do iter = 1, 20 ! start of loop
      ! first part of the iteration
      balanceold = balance
      T_wall_old = T_wall
      T_wall = T_wall * (1d0 + 1d-4) ! increment
      Gradient = 0d0
      Gradient(np,Dir) = ( T - T_wall ) * modfm3
      Gradient(np,:) = matmul ( Gradient(np,:), M )
      rho_wall = Prim(np) / ( Rgas * T_wall )
      Prim_wall(1:nsc) = Prim(1:nsc) * rho_wall / rho
      select case ( fuel )
      case ( 9 ) ! HDPE
        call GSI_HDPE(T_wall, omega, mdot)
      case ( 10 ) ! PP
        call GSI_PP(T_wall, omega, mdot)
      case ( 12 ) ! Carbon
        call GSI_carbon(Prim_wall(1:nsc), T_wall, omega, mdot)
      case ( 13 ) ! HTPB
        call GSI_HTPB(T_wall, omega, mdot)
      end select
      mdot = mdot * (-modfm2)
      omega = omega * (-modfm2)
      call co_k_mi_lam_Wilke ( Prim_wall(1:nsc), rho_wall, T_wall, mil, kl )
      q_conv = kl * dot_product( Gradient(np,:), normal )
      select case ( fuel ) ! pyrolysis heat flux
      case ( 9 ) ! HDPE
        q_pyro = mdot * ( Dh_HDPE + cp_HDPE* (T_wall - Ti_HDPE) ) 
      case ( 10 ) ! PP
        q_pyro = mdot * ( Dh_PP + cp_PP* (T_wall - Ti_PP) )
      case ( 12 ) ! Carbon
        q_pyro = 0d0
        do s = 1, nsc
          q_pyro = q_pyro + omega(s)*f_tabT( T_wall,s,h_tab )
        end do    
      case ( 13 ) ! HTPB
        q_pyro = mdot * ( Dh_HTPB + cp_HTPB* (T_wall - Ti_HTPB) )
      end select
      balance = q_conv - q_pyro - q_rad

      ! new wall temperature (Newton):
      T_wall = T_wall_old - balanceold * (1d-4*T_wall_old)/(balance-balanceold)

      !if ( Abs ( balance ) < 1d-5 ) exit ! stop criteria

      ! second part of the iteration
      Gradient = 0d0
      Gradient(np,Dir) = ( T - T_wall ) * modfm3
      Gradient(np,:) = matmul ( Gradient(np,:), M )
      rho_wall = Prim(np) / ( Rgas * T_wall )
      Prim_wall(1:nsc) = Prim(1:nsc) * rho_wall / rho
      select case ( fuel )
      case ( 9 ) ! HDPE
        call GSI_HDPE(T_wall, omega, mdot)
      case ( 10 ) ! PP
        call GSI_PP(T_wall, omega, mdot)
      case ( 12 ) ! Carbon
        call GSI_carbon(Prim_wall(1:nsc), T_wall, omega, mdot)
      case ( 13 ) ! HTPB
        call GSI_HTPB(T_wall, omega, mdot)
      end select      
      mdot = mdot * (-modfm2)
      omega = omega * (-modfm2)
      call co_k_mi_lam_Wilke ( Prim_wall(1:nsc), rho_wall, T_wall, mil, kl )
      q_conv = kl * dot_product( Gradient(np,:), normal )
      select case ( fuel ) ! pyrolysis heat flux
      case ( 9 ) ! HDPE
        q_pyro = mdot * ( Dh_HDPE + cp_HDPE* (T_wall - Ti_HDPE) ) 
      case ( 10 ) ! PP
        q_pyro = mdot * ( Dh_PP + cp_PP* (T_wall - Ti_PP) )
      case ( 12 ) ! Carbon
        q_pyro = 0d0
        do s = 1, nsc
          q_pyro = q_pyro + omega(s)*f_tabT( T_wall,s,h_tab )
        end do    
      case ( 13 ) ! HTPB
        q_pyro = mdot * ( Dh_HTPB + cp_HTPB* (T_wall - Ti_HTPB) )
      end select
      balance = q_conv - q_pyro - q_rad
    end do

    ! From here on T_wall is known
    Gradient = 0d0
    Gradient(np,Dir) = ( T - T_wall ) * modfm3
    Gradient(np,:) = matmul ( Gradient(np,:), M )
    q_conv = kl * dot_product(Gradient(np,:), normal)
    rho_wall = Prim(np) / ( Rgas * T_wall )
    Prim_wall(1:nsc) = Prim(1:nsc) * rho_wall / rho

    Blowing = mdot / rho_wall * Normal                          ! Blowing velocity
    Gradient(nu:nw,Dir) = ( Prim(nu:nw) - Blowing ) * modfm3    ! Velocity gradient due to blowing in direction normal to the face
    Gradient(nu:nw,:) = matmul ( Gradient(nu:nw,:), M )         ! Velocity gradient transformation
    Stress = Stress_Vector ( Gradient(nu:nw,:), Normal, mil, 0d0, Prim(nt:) )

    Sum_Omegai_Hwi = 0d0
    do s = 1, nsc
      Sum_Omegai_Hwi = Sum_Omegai_Hwi + omega(s) * f_tabT( T_wall,s,h_tab )
    end do

    ! Fluxes
    Flux = 0.0d0
    Flux(1:nsc) = - omega * Area
    Flux(nu:nw) = ( Stress - mdot * Blowing ) * Area
    Flux(np) = Area * ( q_conv - Sum_Omegai_Hwi - mdot * 0.5d0 * sum( Blowing**2 ) + &
                        dot_product( Stress, Blowing ) ) ! check sign q_rad!!!

    if (RANS) then
      dl = blk % yn(im,jm,km)
      call RANS_Set_Blowing_Wall ( rho=rho_wall, mil=mil, &
                                   rans_variables=Prim_Wall(nt:nprim), &
                                   tau=stress, mdot=mdot, dist=dl )
      Gradient(nt:nprim,Dir) = ( Prim(nt:nprim)/rho - Prim_Wall(nt:nprim)/rho_wall ) * modfm3
      Gradient(nt:nprim,:) = matmul ( Gradient(nt:nprim,:), M )
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), &
                                 rans_variables=Prim_Wall(nt:nprim), &
                                 vel_gradient=Gradient(nu:nw,:), &
                                 rans_gradient=Gradient(nt:nprim,:), &
                                 mul=mil, rho=rho_wall, &
                                 area=area, normal=normal, dist=1d-20 )
      ! Ghost-cell extrapolation of RANS variables
      ig = im - guide(fm,1)
      jg = jm - guide(fm,2)
      kg = km - guide(fm,3)
      call RANS_Extrapolate_Wall ( Prim(nt:nprim), Prim_Wall(nt:nprim), &
                                   rho, rho_wall, Blk % P(nt:nprim,Ig,Jg,Kg) )
    endif

    ! Residual update
    Blk % r (:,Im,Jm,Km) = Blk % r (:,Im,Jm,Km) - modfm2 * Flux

    if (present(Ovar)) call Compute_Wall_Properties(stress=Stress, pw=prim(np), qw=q_conv, &
                                                    mdot=mdot, y=blk%dl(im,jm,km)%c(dir)*0.5d0, &
                                                    Tw=T_wall, rhow=rho_wall, mu=mil, exit_array=Ovar)

  end subroutine BC_Wall_Pyrolysis














  
  subroutine Compute_Wall_Properties(stress,pw,Tw,rhow,mu,qw,mdot,y,exit_array)

    implicit none
    real(R8), dimension(:), intent(in)    :: stress
    real(R8), intent(in)                  :: pw
    real(R8), dimension(8), intent(inout) :: exit_array
    real(R8), optional :: Tw, qw, rhow, mu, y, mdot
    real(R8) :: tauw, ut


    ! y+
    tauw = sqrt( sum ( Stress**2 ) )
    ut = sqrt ( tauw/rhow )
    exit_array(1) = rhow*y*ut/mu
    ! Wall stress
    exit_array(2:4) = stress
    ! Wall pressure
    exit_array(5) = pw
    ! Wall temperature
    if (present(Tw)) exit_array(6) = Tw
    ! Heat flux
    if (present(qw)) exit_array(7) = qw
    ! Mass flux
    if (present(mdot)) exit_array(8) = mdot

  end subroutine Compute_Wall_Properties



  !============================================================================
  ! Device twins of the boundary-flux routines.
  !
  ! These mirror the host BC_Symmetry_Eul / BC_Symmetry_Visc / BC_Extrapolation
  ! exactly, but: (a) take PLAIN device-mapped arrays (the block P/r + flat SoA
  ! geometry fm/dl_f/m_f/yn) instead of the AoS derived type; (b) devirtualize
  ! the riemann (HLLC) / limiter (MINMOD) procedure pointers and use the _dev
  ! fixed-size thermo; (c) use NPRIM_DEV-padded locals (no per-call device
  ! malloc). They carry `acc routine seq` and are driven by BC_Fluxes_blk
  ! (Mod_BC). All geometry (face normal/area, cell dl, guide step, modfm) is
  ! precomputed by the caller and passed in, so the routines stay arithmetic.
  !
  ! Geometry-from-Fm arithmetic (no `guide`/`Face_Index` table on the device):
  !   dir    = (Fm+1)/2                 ! 1,2,3 for face pairs (1,2),(3,4),(5,6)
  !   modfm  = mod(Fm,2)                ! 1 on low faces (1,3,5), 0 on high
  !   modfm2 = 1 - 2*modfm              ! -1 low, +1 high  (host's modfm2)
  !   gstep  = 2*modfm - 1              ! +1 low, -1 high  (host's guide() step)
  ! Validated against the host routines on dettube (ROTATING_FRAME/RSM absent in
  ! the device subset; Un_wall=0, no RSM branch).
  !============================================================================

  !> Device twin of BC_Symmetry_Eul (no rotating frame). Updates r(nu:nw).
  subroutine BC_Symmetry_Eul ( PWARG, r, n1, n2, n3, Im, Jm, Km, Fm, &
                                   normal, area, Bound_Dl, Int_Dl )
    !$acc routine seq
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, f_ss_dev
    use MOSE_Lib_Riemann,              only: Riemann, rs_dev
    use MOSE_Lib_Limiters,             only: rlimiter_MINMOD, rlimiter, lim_dev

    implicit none
    integer,  intent(in)    :: n1, n2, n3, Im, Jm, Km, Fm
    real(R8), intent(in)    :: normal(3), area, Bound_Dl, Int_Dl
    integer,  intent(in)    :: ncw, nrw, boff, pi, pj, wroff
    integer,  intent(in)    :: nur
    real(R8), intent(in)    :: ur(nur)   ! Ur of the cells (prec; 1 element otherwise): read by HLLC-PC only
    real(R8), intent(in)    :: w(ncw, nprim)
    real(R8), intent(inout) :: wr(nrw, nprim)   ! Res var-last: coalesced residual WRITE target
    real(R8), intent(inout) :: r(nprim, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    ! Local (NPRIM_DEV fixed-size; only 1:nprim / 1:np used)
    integer  :: dir, modfm, modfm1, modfm2, gstep, Int_i, Int_j, Int_k, s
    real(R8) :: Bound_Prim(NPRIM_DEV), Int_Prim(NPRIM_DEV), Ghost_Prim(NPRIM_DEV)
    real(R8) :: Prim1(NPRIM_DEV), Prim3(NPRIM_DEV), Face_Prim(NPRIM_DEV), Slope(NPRIM_DEV)
    real(R8) :: Slope32, Slope21, Un, Velt(3), Dist12, Dist23, Dist
    real(R8) :: rotot, Rtot, Sound, Un1, F_r, F_u, F_v, F_w, F_E, U1(3), U4(3), urL_, urR_

    dir    = (Fm + 1)/2
    modfm  = mod(Fm, 2)
    modfm1 = 1 - modfm
    modfm2 = 1 - 2*modfm
    gstep  = 2*modfm - 1
    Int_i = Im; Int_j = Jm; Int_k = Km
    select case (dir)
      case (1); Int_i = Im + gstep
      case (2); Int_j = Jm + gstep
      case (3); Int_k = Km + gstep
    end select

    do s = 1, nprim
      Bound_Prim(s) = PCELL(s, Im, Jm, Km)
      Int_Prim(s)   = PCELL(s, Int_i, Int_j, Int_k)
    end do

    ! Boundary cell normal velocity + tangential component
    Un = Bound_Prim(nu)*normal(1) + Bound_Prim(nv)*normal(2) + Bound_Prim(nw)*normal(3)

    ! Ghost cell: reflect the normal velocity about the wall (Un_wall = 0)
    do s = 1, nprim
      Ghost_Prim(s) = Bound_Prim(s)
    end do
    Ghost_Prim(nu) = Bound_Prim(nu) - 2d0*Un*normal(1)
    Ghost_Prim(nv) = Bound_Prim(nv) - 2d0*Un*normal(2)
    Ghost_Prim(nw) = Bound_Prim(nw) - 2d0*Un*normal(3)

    ! Stencil + distances (point (1),(2)=boundary,(3))
    do s = 1, nprim
      Prim1(s) = modfm*Ghost_Prim(s) + modfm1*Int_Prim(s)
      Prim3(s) = modfm*Int_Prim(s)   + modfm1*Ghost_Prim(s)
    end do
    Dist12 = 2d0*Bound_Dl*modfm + (Bound_Dl + Int_Dl)*modfm1
    Dist23 = (Bound_Dl + Int_Dl)*modfm + 2d0*Bound_Dl*modfm1

    ! Limited reconstruction at the boundary face (components 1:np)
    do s = 1, np
      Slope32  = (Prim3(s) - Bound_Prim(s)) / Dist23
      Slope21  = (Bound_Prim(s) - Prim1(s)) / Dist12
      Slope(s) = rlimiter ( lim_dev,  Slope32, Slope21 )
    end do
    Dist = Bound_Dl * modfm2
    do s = 1, np
      Face_Prim(s) = Bound_Prim(s) + Slope(s)*Dist
    end do

    ! Riemann problem at the symmetry interface (contact surface at the wall)
    call co_rotot_Rtot_dev ( Face_Prim(1:nsc), rotot, Rtot )
    Sound = f_ss_dev ( Face_Prim(1:nsc), Face_Prim(np), rotot, Rtot )
    Un = Face_Prim(nu)*normal(1) + Face_Prim(nv)*normal(2) + Face_Prim(nw)*normal(3)
    Velt(1) = Face_Prim(nu) - Un*normal(1)
    Velt(2) = Face_Prim(nv) - Un*normal(2)
    Velt(3) = Face_Prim(nw) - Un*normal(3)
    Un1 = Un * modfm2
    U1(1) =  Un1*normal(1) + Velt(1); U1(2) =  Un1*normal(2) + Velt(2); U1(3) =  Un1*normal(3) + Velt(3)
    U4(1) = -Un1*normal(1) + Velt(1); U4(2) = -Un1*normal(2) + Velt(2); U4(3) = -Un1*normal(3) + Velt(3)
    urL_ = 0d0; urR_ = 0d0
    if (rs_dev == 14) then   ! HLLC-PC: the Ur of the cell on both sides (MOSE_open Lib_BC_Fluxes_Symmetry)
      urL_ = URCELL(Im,Jm,Km); urR_ = urL_
    end if
    call Riemann ( rs_dev,  Face_Prim(1:nsc), U1(1), U1(2), U1(3), Face_Prim(np), Sound, rotot, &
                        Face_Prim(1:nsc), U4(1), U4(2), U4(3), Face_Prim(np), Sound, rotot, &
                        1d0, urL_, urR_, normal(1), normal(2), normal(3), F_r, F_u, F_v, F_w, F_E )

    ! Residual update (momentum only)
    WRCELL(nu,Im,Jm,Km) = WRCELL(nu,Im,Jm,Km) + modfm2 * F_u * area
    WRCELL(nv,Im,Jm,Km) = WRCELL(nv,Im,Jm,Km) + modfm2 * F_v * area
    WRCELL(nw,Im,Jm,Km) = WRCELL(nw,Im,Jm,Km) + modfm2 * F_w * area

  end subroutine BC_Symmetry_Eul


  !> Device twin of BC_Symmetry_Visc (SST, Wilcox2006, SA; RSM branch absent). Updates
  !> r(nu:np). mil/mie/stress mirror the device interior-diffusion path
  !> (co_k_mi_lam_Wilke_dev for the Wilke laminar viscosity, devirt SST eddy
  !> viscosity + stress_dev tensor) so the BC and the volume fluxes are mutually
  !> consistent on the device.
  subroutine BC_Symmetry_Visc ( PWARG, r, Mrec, ynrec, n1, n2, n3, Im, Jm, Km, Fm, normal, area )
    !$acc routine seq
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, co_k_mi_lam_Wilke_dev
    use MOSE_Mod_SST,                  only: SST_Eddy_Viscosity
    use MOSE_Mod_Spalart,              only: Spalart_Eddy_Viscosity
    use MOSE_Mod_Wilcox2006,           only: Wilcox2006_Eddy_Viscosity
    use MOSE_Mod_Fluid,                only: Stress_Vector_2eq, Stress_Vector_Std
    use MOSE_Lib_RANS,                 only: stress_dev, eddy_dev

    implicit none
    integer,  intent(in)    :: n1, n2, n3, Im, Jm, Km, Fm
    real(R8), intent(in)    :: normal(3), area
    integer,  intent(in)    :: ncw, nrw, boff, pi, pj, wroff
    integer,  intent(in)    :: nur
    real(R8), intent(in)    :: ur(nur)   ! Ur of the cells (prec; 1 element otherwise): read by HLLC-PC only
    real(R8), intent(in)    :: w(ncw, nprim)
    real(R8), intent(inout) :: wr(nrw, nprim)   ! Res var-last: coalesced residual WRITE target
    real(R8), intent(inout) :: r(nprim, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    real(R8), intent(in)    :: Mrec(3,3)   !< 5r-b: m_f at this record's cell
    real(R8), intent(in)    :: ynrec       !< 5r-b: yn at the same cell
    ! Local
    integer  :: dir, modfm, modfm2, a, b, c, s, ii, kk, Tint(2)
    real(R8) :: M(3,3), Grad(3,3), G2(3,3), dist, rho, Rgas, mil, kl, mie, T, Tdiff
    real(R8) :: Vel1(3), Vel2(3), Vel3(3), Vel4(3), Vel5(3), Vel6(3)
    real(R8) :: Vel7(3), Vel8(3), Vel9(3), Vel10(3), Vel(3), Stress(3)
    real(R8) :: Prim_b(NPRIM_DEV), d1, d3, d4, d7, d8

    dir    = (Fm + 1)/2
    modfm  = mod(Fm, 2)
    modfm2 = 1 - 2*modfm

    do s = 1, nprim
      Prim_b(s) = PCELL(s, Im, Jm, Km)
    end do
    M    = Mrec
    dist = ynrec

    call co_rotot_Rtot_dev ( Prim_b(1:nsc), rho, Rgas )
    T       = Prim_b(np) / (rho*Rgas)
    Tint(1) = int(T); Tint(2) = Tint(1) + 1; Tdiff = T - Tint(1)
    call co_k_mi_lam_Wilke_dev ( Prim_b(1:nsc), rho, Tint, Tdiff, mil, kl )

    ! Velocity stencil + axis mapping (mirror host select case Fm)
    Vel1(1) = PCELL(nu,Im,Jm,Km); Vel1(2) = PCELL(nv,Im,Jm,Km); Vel1(3) = PCELL(nw,Im,Jm,Km)
    select case (Fm)
      case (1:2)
        Vel3(1)=PCELL(nu,Im,Jm-1,Km); Vel3(2)=PCELL(nv,Im,Jm-1,Km); Vel3(3)=PCELL(nw,Im,Jm-1,Km)
        Vel4(1)=PCELL(nu,Im,Jm+1,Km); Vel4(2)=PCELL(nv,Im,Jm+1,Km); Vel4(3)=PCELL(nw,Im,Jm+1,Km)
        Vel7(1)=PCELL(nu,Im,Jm,Km-1); Vel7(2)=PCELL(nv,Im,Jm,Km-1); Vel7(3)=PCELL(nw,Im,Jm,Km-1)
        Vel8(1)=PCELL(nu,Im,Jm,Km+1); Vel8(2)=PCELL(nv,Im,Jm,Km+1); Vel8(3)=PCELL(nw,Im,Jm,Km+1)
        a = 1; b = 2; c = 3
      case (3:4)
        Vel3(1)=PCELL(nu,Im-1,Jm,Km); Vel3(2)=PCELL(nv,Im-1,Jm,Km); Vel3(3)=PCELL(nw,Im-1,Jm,Km)
        Vel4(1)=PCELL(nu,Im+1,Jm,Km); Vel4(2)=PCELL(nv,Im+1,Jm,Km); Vel4(3)=PCELL(nw,Im+1,Jm,Km)
        Vel7(1)=PCELL(nu,Im,Jm,Km-1); Vel7(2)=PCELL(nv,Im,Jm,Km-1); Vel7(3)=PCELL(nw,Im,Jm,Km-1)
        Vel8(1)=PCELL(nu,Im,Jm,Km+1); Vel8(2)=PCELL(nv,Im,Jm,Km+1); Vel8(3)=PCELL(nw,Im,Jm,Km+1)
        a = 2; b = 1; c = 3
      case (5:6)
        Vel3(1)=PCELL(nu,Im-1,Jm,Km); Vel3(2)=PCELL(nv,Im-1,Jm,Km); Vel3(3)=PCELL(nw,Im-1,Jm,Km)
        Vel4(1)=PCELL(nu,Im+1,Jm,Km); Vel4(2)=PCELL(nv,Im+1,Jm,Km); Vel4(3)=PCELL(nw,Im+1,Jm,Km)
        Vel7(1)=PCELL(nu,Im,Jm-1,Km); Vel7(2)=PCELL(nv,Im,Jm-1,Km); Vel7(3)=PCELL(nw,Im,Jm-1,Km)
        Vel8(1)=PCELL(nu,Im,Jm+1,Km); Vel8(2)=PCELL(nv,Im,Jm+1,Km); Vel8(3)=PCELL(nw,Im,Jm+1,Km)
        a = 3; b = 1; c = 2
    end select

    ! Symmetry-reflected mirror velocities
    d1 = Vel1(1)*normal(1) + Vel1(2)*normal(2) + Vel1(3)*normal(3)
    d3 = Vel3(1)*normal(1) + Vel3(2)*normal(2) + Vel3(3)*normal(3)
    d4 = Vel4(1)*normal(1) + Vel4(2)*normal(2) + Vel4(3)*normal(3)
    d7 = Vel7(1)*normal(1) + Vel7(2)*normal(2) + Vel7(3)*normal(3)
    d8 = Vel8(1)*normal(1) + Vel8(2)*normal(2) + Vel8(3)*normal(3)
    do s = 1, 3
      Vel2(s)  = Vel1(s) - 2d0*d1*normal(s)
      Vel5(s)  = Vel3(s) - 2d0*d3*normal(s)
      Vel6(s)  = Vel4(s) - 2d0*d4*normal(s)
      Vel9(s)  = Vel7(s) - 2d0*d7*normal(s)
      Vel10(s) = Vel8(s) - 2d0*d8*normal(s)
    end do

    do s = 1, 3
      Grad(s,a) = (Vel2(s) - Vel1(s)) * modfm2
      Grad(s,b) = (Vel4(s) - Vel3(s) + Vel6(s) - Vel5(s)) * 0.25d0
      Grad(s,c) = (Vel8(s) - Vel7(s) + Vel10(s) - Vel9(s)) * 0.25d0
    end do

    ! G2 = matmul(Grad, M)  ->  G2(ii,kk) = sum_s Grad(ii,s)*M(s,kk)
    do kk = 1, 3
      do ii = 1, 3
        G2(ii,kk) = Grad(ii,1)*M(1,kk) + Grad(ii,2)*M(2,kk) + Grad(ii,3)*M(3,kk)
      end do
    end do

    mie = 0d0
    if (RANS) then
      if (eddy_dev == 2) then
        call Spalart_Eddy_Viscosity ( mie, Prim_b(nt:nt+nrans-1), mil, rho, G2, dist )
      else if (eddy_dev == 3) then
        call Wilcox2006_Eddy_Viscosity ( mie, Prim_b(nt:nprim), mil, rho, G2, dist )
      else
        call SST_Eddy_Viscosity ( mie, Prim_b(nt:nprim), mil, rho, G2, dist )
      end if
    end if

    if (stress_dev == 2) then
      Stress = Stress_Vector_2eq ( G2, normal, mil, mie, Prim_b(nt:nprim) )
    else
      Stress = Stress_Vector_Std ( G2, normal, mil, mie, Prim_b(nt:nprim) )
    end if

    Vel(1) = 0.5d0*(Vel1(1) + Vel2(1))
    Vel(2) = 0.5d0*(Vel1(2) + Vel2(2))
    Vel(3) = 0.5d0*(Vel1(3) + Vel2(3))

    ! Residual update (viscous momentum + energy)
    WRCELL(nu,Im,Jm,Km) = WRCELL(nu,Im,Jm,Km) - modfm2 * area * Stress(1)
    WRCELL(nv,Im,Jm,Km) = WRCELL(nv,Im,Jm,Km) - modfm2 * area * Stress(2)
    WRCELL(nw,Im,Jm,Km) = WRCELL(nw,Im,Jm,Km) - modfm2 * area * Stress(3)
    WRCELL(np,Im,Jm,Km) = WRCELL(np,Im,Jm,Km) - modfm2 * area * &
                     ( Vel(1)*Stress(1) + Vel(2)*Stress(2) + Vel(3)*Stress(3) )

  end subroutine BC_Symmetry_Visc


  !> Device twin of BC_Wall_Heat for the LAMINAR, ADIABATIC no-slip wall (type 5).
  !> Mirrors BC_Symmetry_Visc's prologue (PCELL read; rho/Rgas/T; mil/kl via the
  !> Wilke device path, the same form the host BC_Wall_Heat uses) but the wall physics
  !> is SIMPLER than symmetry:
  !>   - no-slip (wall velocity = 0): the only velocity gradient is the wall-normal one
  !>     from the boundary cell, Grad(s,dir) = Prim(nu:nw)*modfm3; tangential columns 0
  !>     (host: Gradient(nu:nw,Dir) = Prim(nu:nw)*modfm3, the rest 0).
  !>   - adiabatic (Heat_Flux = 0): Flux(np) = Area*Heat_Flux = 0 and the no-slip wall
  !>     has Vel = 0 so there is no viscous-work energy term -> energy residual untouched.
  !>   - mie = 0 at the wall MOMENTUM stress (host passes 0d0 to Stress_Vector regardless of
  !>     model; for SST stress_dev=2 -> Stress_Vector_2eq adds the -2/3*rho*k term from Prim_b).
  !>   - k-omega (SST, Wilcox2006) wall: adds the k/omega diffusive flux into r(nt:nprim) via
  !>     SST_ / Wilcox2006_RANS_Diffusive_Flux. Adiabatic -> Twall=T -> rho_wall=rho. The matching omega-wall
  !>     GHOST is set in Extrapolate_Wall_blk.
  !> Updates r(nu:nw) (+ r(nt:nprim) for SST). Write target mirrors BC_Symmetry_Visc.
  subroutine BC_Wall_Heat_dev ( PWARG, r, Mrec, ynrec, n1, n2, n3, Im, Jm, Km, Fm, normal, area )
    !$acc routine seq
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, co_k_mi_lam_Wilke_dev
    use MOSE_Mod_Fluid,                only: Stress_Vector_2eq, Stress_Vector_Std
    use MOSE_Mod_SST,                  only: SST_RANS_Diffusive_Flux
    use MOSE_Mod_Wilcox2006,           only: Wilcox2006_RANS_Diffusive_Flux
    use MOSE_Mod_Spalart,              only: sa_sigma => sigma
    use MOSE_Lib_RANS,                 only: stress_dev, eddy_dev

    implicit none
    integer,  intent(in)    :: n1, n2, n3, Im, Jm, Km, Fm
    real(R8), intent(in)    :: normal(3), area
    integer,  intent(in)    :: ncw, nrw, boff, pi, pj, wroff
    integer,  intent(in)    :: nur
    real(R8), intent(in)    :: ur(nur)   ! Ur of the cells (prec; 1 element otherwise): read by HLLC-PC only
    real(R8), intent(in)    :: w(ncw, nprim)
    real(R8), intent(inout) :: wr(nrw, nprim)   ! Res var-last: coalesced residual WRITE target
    real(R8), intent(inout) :: r(nprim, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    real(R8), intent(in)    :: Mrec(3,3)   !< 5r-b: m_f at this record's cell
    real(R8), intent(in)    :: ynrec       !< 5r-b: yn at the same cell
    ! Local
    integer  :: dir, modfm, modfm2, modfm3, s, ii, kk, rr, Tint(2)
    real(R8) :: M(3,3), Grad(3,3), G2(3,3), rho, Rgas, mil, kl, T, Tdiff
    real(R8) :: Vel1(3), Stress(3), Prim_b(NPRIM_DEV)
    real(R8) :: rkw_w(2), GradR(2,3), GradR2(2,3), Flux_r(2), dist_w   ! k-omega (SST, Wilcox2006: nRANS=2)

    dir    = (Fm + 1)/2
    modfm  = mod(Fm, 2)
    modfm2 = 1 - 2*modfm
    modfm3 = 2*(2*modfm - 1)            ! host Compute_Modfm: = -2*modfm2

    do s = 1, nprim
      Prim_b(s) = PCELL(s, Im, Jm, Km)
    end do
    M = Mrec

    call co_rotot_Rtot_dev ( Prim_b(1:nsc), rho, Rgas )
    T       = Prim_b(np) / (rho*Rgas)
    Tint(1) = int(T); Tint(2) = Tint(1) + 1; Tdiff = T - Tint(1)
    call co_k_mi_lam_Wilke_dev ( Prim_b(1:nsc), rho, Tint, Tdiff, mil, kl )

    ! No-slip wall: only the wall-normal velocity gradient from the boundary cell
    ! (wall velocity = 0). Tangential gradient columns stay zero.
    Vel1(1) = Prim_b(nu); Vel1(2) = Prim_b(nv); Vel1(3) = Prim_b(nw)
    Grad = 0d0
    do s = 1, 3
      Grad(s, dir) = Vel1(s) * modfm3
    end do

    ! G2 = matmul(Grad, M)  ->  G2(ii,kk) = sum_s Grad(ii,s)*M(s,kk)
    do kk = 1, 3
      do ii = 1, 3
        G2(ii,kk) = Grad(ii,1)*M(1,kk) + Grad(ii,2)*M(2,kk) + Grad(ii,3)*M(3,kk)
      end do
    end do

    ! Adiabatic laminar wall: mie = 0 (no eddy viscosity at the wall).
    if (stress_dev == 2) then
      Stress = Stress_Vector_2eq ( G2, normal, mil, 0d0, Prim_b(nt:nprim) )
    else
      Stress = Stress_Vector_Std ( G2, normal, mil, 0d0, Prim_b(nt:nprim) )
    end if

    ! Residual update (viscous momentum only; energy flux = Area*Heat_Flux = 0,
    ! and the no-slip wall carries no viscous-work term).
    WRCELL(nu,Im,Jm,Km) = WRCELL(nu,Im,Jm,Km) - modfm2 * area * Stress(1)
    WRCELL(nv,Im,Jm,Km) = WRCELL(nv,Im,Jm,Km) - modfm2 * area * Stress(2)
    WRCELL(nw,Im,Jm,Km) = WRCELL(nw,Im,Jm,Km) - modfm2 * area * Stress(3)

    ! --- k-omega (SST, Wilcox2006) wall: k/omega diffusive flux (host BC_Wall_Heat `if (RANS)` branch) ---
    ! k-omega (nRANS=2). Adiabatic no-slip wall -> Twall=T -> rho_wall=rho, so the
    ! wall-density factors cancel. The omega-wall GHOST is set in Extrapolate_Wall_blk.
    if (RANS .and. eddy_dev == 2) then
      ! SA wall (Spalart_Set_Wall_Values: nu-tilde_wall = 0; Spalart_RANS_Diffusive_Flux): the wall-normal gradient only
      GradR(1,1) = 0d0; GradR(1,2) = 0d0; GradR(1,3) = 0d0
      GradR(1,dir) = Prim_b(nt) / rho * modfm3
      do kk = 1, 3
        GradR2(1,kk) = GradR(1,dir) * M(dir,kk)
      end do
      Flux_r(1) = area * mil * ( GradR2(1,1)*normal(1) + GradR2(1,2)*normal(2) + GradR2(1,3)*normal(3) ) / sa_sigma
      WRCELL(nt,Im,Jm,Km) = WRCELL(nt,Im,Jm,Km) - modfm2 * Flux_r(1)
    else if (RANS) then
      dist_w   = ynrec
      rkw_w(1) = 0d0                            ! SST_ / Wilcox2006_Set_Wall_Values: k_wall = 0
      rkw_w(2) = omega_wall_c(1) * mil / (dist_w*dist_w)    !   omega_wall (rkw = [rho*k, rho*w])
      ! Wall-normal RANS gradient (host: (Prim_b(nt:)/rho - rkw_w/rho_wall)*modfm3, rho_wall=rho);
      ! tangential columns 0, so matmul(.,M) keeps only the `dir` column.
      do rr = 1, 2
        GradR(rr,1) = 0d0; GradR(rr,2) = 0d0; GradR(rr,3) = 0d0
        GradR(rr,dir) = (Prim_b(nt-1+rr) - rkw_w(rr)) / rho * modfm3
        do kk = 1, 3
          GradR2(rr,kk) = GradR(rr,dir) * M(dir,kk)
        end do
      end do
      if (eddy_dev == 3) then   ! Wilcox2006 (at k_wall = 0 both models give mu_l grad . n)
        call Wilcox2006_RANS_Diffusive_Flux ( Flux_r, rkw_w, mil, G2, GradR2, rho, area, normal, 1d-20 )
      else
        call SST_RANS_Diffusive_Flux ( Flux_r, rkw_w, mil, G2, GradR2, rho, area, normal, 1d-20 )
      end if
      WRCELL(nt,  Im,Jm,Km) = WRCELL(nt,  Im,Jm,Km) - modfm2 * Flux_r(1)
      WRCELL(nt+1,Im,Jm,Km) = WRCELL(nt+1,Im,Jm,Km) - modfm2 * Flux_r(2)
    end if

  end subroutine BC_Wall_Heat_dev


  !> Device twin of BC_Wall_Temperature: no-slip wall at a PRESCRIBED temperature
  !> (type 6). Structure mirrors BC_Wall_Heat_dev; the host deltas it carries over:
  !>   - rho_wall = p/(Rgas*T_wall) (boundary-layer p(wall)=p(cell), non-catalytic
  !>     composition), NOT rho: the wall-density factors no longer cancel.
  !>   - mil/kl from the Wilke device path at the WALL state (scaled species, T_wall)
  !>     -- the host computes them at (Prim_wall, rho_wall, T_wall).
  !>   - energy flux = Area * kl * dot(Grad_T_row . M, normal), with the wall-normal
  !>     temperature gradient row Grad(np,dir) = (T - T_wall)*modfm3 (host line-for-line).
  !>   - SST: rkw_w(2) = omega_wall_c(1)*mil/(dist^2) with mil at T_wall, and the RANS gradient keeps
  !>     both densities: (Prim_b/rho - rkw_w/rho_wall)*modfm3; flux called with rho_wall.
  !> The matching omega-wall GHOST (with the same rho/rho_wall factor) is set in
  !> Extrapolate_Wall_* (device BC kernels cannot write P/w ghosts).
  subroutine BC_Wall_Temperature_dev ( PWARG, r, Mrec, ynrec, n1, n2, n3, Im, Jm, Km, Fm, normal, area, BC_Tw )
    !$acc routine seq
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, co_k_mi_lam_Wilke_dev, f_cp_dev
    use MOSE_Mod_Fluid,                only: Stress_Vector_2eq, Stress_Vector_Std
    use MOSE_Mod_SST,                  only: SST_RANS_Diffusive_Flux
    use MOSE_Mod_Wilcox2006,           only: Wilcox2006_RANS_Diffusive_Flux
    use MOSE_Mod_Spalart,              only: sa_sigma => sigma
    use MOSE_Lib_RANS,                 only: stress_dev, eddy_dev, prl_dev

    implicit none
    integer,  intent(in)    :: n1, n2, n3, Im, Jm, Km, Fm
    real(R8), intent(in)    :: normal(3), area, BC_Tw
    integer,  intent(in)    :: ncw, nrw, boff, pi, pj, wroff
    integer,  intent(in)    :: nur
    real(R8), intent(in)    :: ur(nur)   ! Ur of the cells (prec; 1 element otherwise): read by HLLC-PC only
    real(R8), intent(in)    :: w(ncw, nprim)
    real(R8), intent(inout) :: wr(nrw, nprim)   ! Res var-last: coalesced residual WRITE target
    real(R8), intent(inout) :: r(nprim, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    real(R8), intent(in)    :: Mrec(3,3)   !< 5r-b: m_f at this record's cell
    real(R8), intent(in)    :: ynrec       !< 5r-b: yn at the same cell
    ! Local
    integer  :: dir, modfm, modfm2, modfm3, s, ii, kk, rr, Tint(2)
    real(R8) :: M(3,3), Grad(3,3), G2(3,3), rho, Rgas, mil, kl, T, Tdiff
    real(R8) :: Vel1(3), Stress(3), Prim_b(NPRIM_DEV), Prim_w(NPRIM_DEV)
    real(R8) :: rho_wall, Gnp_dir, Gnp(3), Flux_np
    real(R8) :: rkw_w(2), GradR(2,3), GradR2(2,3), Flux_r(2), dist_w   ! k-omega (SST, Wilcox2006: nRANS=2)

    dir    = (Fm + 1)/2
    modfm  = mod(Fm, 2)
    modfm2 = 1 - 2*modfm
    modfm3 = 2*(2*modfm - 1)            ! host Compute_Modfm: = -2*modfm2

    do s = 1, nprim
      Prim_b(s) = PCELL(s, Im, Jm, Km)
    end do
    M = Mrec

    call co_rotot_Rtot_dev ( Prim_b(1:nsc), rho, Rgas )
    T = Prim_b(np) / (rho*Rgas)

    ! Wall state: p(wall) = p(cell), non-catalytic (same mass fractions), T = BC_Tw.
    rho_wall = Prim_b(np) / (Rgas * BC_Tw)
    do s = 1, nsc
      Prim_w(s) = Prim_b(s) * rho_wall / rho
    end do
    Tint(1) = int(BC_Tw); Tint(2) = Tint(1) + 1; Tdiff = BC_Tw - Tint(1)
    call co_k_mi_lam_Wilke_dev ( Prim_w(1:nsc), rho_wall, Tint, Tdiff, mil, kl )
    if (prl_dev > 0d0) kl = mil * f_cp_dev ( Prim_w(1:nsc), Tint, Tdiff, rho_wall ) / prl_dev   ! [MOSE-Physics] prandtl

    ! No-slip wall: only the wall-normal velocity gradient from the boundary cell
    ! (wall velocity = 0). Tangential gradient columns stay zero.
    Vel1(1) = Prim_b(nu); Vel1(2) = Prim_b(nv); Vel1(3) = Prim_b(nw)
    Grad = 0d0
    do s = 1, 3
      Grad(s, dir) = Vel1(s) * modfm3
    end do

    ! G2 = matmul(Grad, M)  ->  G2(ii,kk) = sum_s Grad(ii,s)*M(s,kk)
    do kk = 1, 3
      do ii = 1, 3
        G2(ii,kk) = Grad(ii,1)*M(1,kk) + Grad(ii,2)*M(2,kk) + Grad(ii,3)*M(3,kk)
      end do
    end do

    ! Laminar wall: mie = 0 (no eddy viscosity at the wall); mil at the WALL state.
    if (stress_dev == 2) then
      Stress = Stress_Vector_2eq ( G2, normal, mil, 0d0, Prim_b(nt:nprim) )
    else
      Stress = Stress_Vector_Std ( G2, normal, mil, 0d0, Prim_b(nt:nprim) )
    end if

    ! Wall-normal temperature gradient row (host Gradient(np,Dir) = (T - T_wall)*modfm3;
    ! tangential columns 0 -> after matmul only the `dir` row of M survives).
    Gnp_dir = (T - BC_Tw) * modfm3
    do kk = 1, 3
      Gnp(kk) = Gnp_dir * M(dir,kk)
    end do
    Flux_np = area * kl * ( Gnp(1)*normal(1) + Gnp(2)*normal(2) + Gnp(3)*normal(3) )

    ! Residual update (viscous momentum + conductive energy; no viscous-work term:
    ! the no-slip wall velocity is 0).
    WRCELL(nu,Im,Jm,Km) = WRCELL(nu,Im,Jm,Km) - modfm2 * area * Stress(1)
    WRCELL(nv,Im,Jm,Km) = WRCELL(nv,Im,Jm,Km) - modfm2 * area * Stress(2)
    WRCELL(nw,Im,Jm,Km) = WRCELL(nw,Im,Jm,Km) - modfm2 * area * Stress(3)
    WRCELL(np,Im,Jm,Km) = WRCELL(np,Im,Jm,Km) - modfm2 * Flux_np

    ! --- k-omega (SST, Wilcox2006) wall: k/omega diffusive flux (host BC_Wall_Temperature `if (RANS)`) ---
    ! k-omega (nRANS=2). Unlike the adiabatic wall the densities do NOT cancel:
    ! GradR = (Prim_b/rho - rkw_w/rho_wall)*modfm3, and the flux is called with rho_wall.
    if (RANS .and. eddy_dev == 2) then
      ! SA wall (Spalart_Set_Wall_Values: nu-tilde_wall = 0; Spalart_RANS_Diffusive_Flux): the wall-normal gradient only
      GradR(1,1) = 0d0; GradR(1,2) = 0d0; GradR(1,3) = 0d0
      GradR(1,dir) = Prim_b(nt) / rho * modfm3
      do kk = 1, 3
        GradR2(1,kk) = GradR(1,dir) * M(dir,kk)
      end do
      Flux_r(1) = area * mil * ( GradR2(1,1)*normal(1) + GradR2(1,2)*normal(2) + GradR2(1,3)*normal(3) ) / sa_sigma
      WRCELL(nt,Im,Jm,Km) = WRCELL(nt,Im,Jm,Km) - modfm2 * Flux_r(1)
    else if (RANS) then
      dist_w   = ynrec
      rkw_w(1) = 0d0                            ! SST_ / Wilcox2006_Set_Wall_Values: k_wall = 0
      rkw_w(2) = omega_wall_c(1) * mil / (dist_w*dist_w)    !   omega_wall (rkw = [rho*k, rho*w]), mil at T_wall
      do rr = 1, 2
        GradR(rr,1) = 0d0; GradR(rr,2) = 0d0; GradR(rr,3) = 0d0
        GradR(rr,dir) = ( Prim_b(nt-1+rr)/rho - rkw_w(rr)/rho_wall ) * modfm3
        do kk = 1, 3
          GradR2(rr,kk) = GradR(rr,dir) * M(dir,kk)
        end do
      end do
      if (eddy_dev == 3) then   ! Wilcox2006 (at k_wall = 0 both models give mu_l grad . n)
        call Wilcox2006_RANS_Diffusive_Flux ( Flux_r, rkw_w, mil, G2, GradR2, rho_wall, area, normal, 1d-20 )
      else
        call SST_RANS_Diffusive_Flux ( Flux_r, rkw_w, mil, G2, GradR2, rho_wall, area, normal, 1d-20 )
      end if
      WRCELL(nt,  Im,Jm,Km) = WRCELL(nt,  Im,Jm,Km) - modfm2 * Flux_r(1)
      WRCELL(nt+1,Im,Jm,Km) = WRCELL(nt+1,Im,Jm,Km) - modfm2 * Flux_r(2)
    end if

  end subroutine BC_Wall_Temperature_dev


  !> Device twin of BC_Extrapolation. Updates r(:) (all components).
  subroutine BC_Extrapolation ( PWARG, r, n1, n2, n3, Im, Jm, Km, Fm, normal, area )
    !$acc routine seq
    use FLINT_Lib_ThermoTransport_dev, only: H0_dev

    implicit none
    integer,  intent(in)    :: n1, n2, n3, Im, Jm, Km, Fm
    real(R8), intent(in)    :: normal(3), area
    integer,  intent(in)    :: ncw, nrw, boff, pi, pj, wroff
    integer,  intent(in)    :: nur
    real(R8), intent(in)    :: ur(nur)   ! Ur of the cells (prec; 1 element otherwise): read by HLLC-PC only
    real(R8), intent(in)    :: w(ncw, nprim)
    real(R8), intent(inout) :: wr(nrw, nprim)   ! Res var-last: coalesced residual WRITE target
    real(R8), intent(inout) :: r(nprim, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    ! Local
    integer  :: s, FluxSign
    real(R8) :: Prim(NPRIM_DEV), Flux(NPRIM_DEV), Un, Fmass, vmag

    do s = 1, nprim
      Prim(s) = PCELL(s, Im, Jm, Km)
    end do

    Un = Prim(nu)*normal(1) + Prim(nv)*normal(2) + Prim(nw)*normal(3)

    Fmass = 0d0
    do s = 1, nsc
      Flux(s) = Prim(s) * Un * area
      Fmass   = Fmass + Flux(s)
    end do
    Flux(nu) = Fmass*Prim(nu) + Prim(np)*area*normal(1)
    Flux(nv) = Fmass*Prim(nv) + Prim(np)*area*normal(2)
    Flux(nw) = Fmass*Prim(nw) + Prim(np)*area*normal(3)
    vmag     = sqrt( Prim(nu)**2 + Prim(nv)**2 + Prim(nw)**2 )
    Flux(np) = Fmass * H0_dev ( Prim(np), Prim(1:nsc), vmag )
    if (nprim > np) then
      do s = np+1, nprim
        Flux(s) = Prim(s) * Un * area
      end do
    end if

    FluxSign = 1 - 2*mod(Fm, 2)
    do s = 1, nprim
      WRCELL(s,Im,Jm,Km) = WRCELL(s,Im,Jm,Km) + FluxSign * Flux(s)
    end do

  end subroutine BC_Extrapolation


  !> Device twin of BC_Outflow. Subsonic outflow to p_Amb (isentropic) or pure
  !> extrapolation; fills Flux(1:nprim)+Fmass. Mirror of Lib_BC BC_Outflow.
  subroutine BC_Outflow ( Prim, rho, Rgas, Sound, Un, Gamma, Riem, Normal, Area, &
                              p_Amb, BC_rel_fac, BC_Sign, Flux, Fmass )
    !$acc routine seq
    use FLINT_Lib_Thermodynamic, only: f_tabT_expr, h_tab

    implicit none
    real(R8), intent(in)  :: Prim(NPRIM_DEV), rho, Rgas, Sound, Un, Gamma, Riem
    real(R8), intent(in)  :: Normal(3), Area, p_Amb, BC_Sign, BC_rel_fac
    real(R8), intent(out) :: Flux(NPRIM_DEV), Fmass
    integer  :: s, T_i, Tint(2)
    real(R8) :: Vel_T(3), Mach_N, Face_Sound, Isentropic, Face_Un, Tdiff
    real(R8) :: Face_Prim(NPRIM_DEV), Face_rho, Face_Rgas, Face_Enthalpy, Face_T
    real(R8) :: Mach_abs, pstar, denom, delta, p_exit, vmag, vt2

    vmag     = sqrt( Prim(nu)**2 + Prim(nv)**2 + Prim(nw)**2 )
    Vel_T(1) = Prim(nu) - Un*Normal(1)
    Vel_T(2) = Prim(nv) - Un*Normal(2)
    Vel_T(3) = Prim(nw) - Un*Normal(3)
    Mach_abs = vmag / Sound
    Mach_N   = abs( Un / Sound )
    p_exit   = p_Amb

    if ( (Mach_N < 1d0) .and. (p_Amb > 0d0) ) then
      delta = 0.5d0 * (Gamma - 1d0)
      denom = 1d0 / (1d0 + delta)**(Gamma/(Gamma-1d0))
      ! p* from the Mach number normal to the face, as MOSE_open Outflow (Lib_BC_Fluxes_Outflow.f90:93); it was the
      ! total Mach number (A3 of the alignment of 2026-10-16)
      pstar = Prim(np) * denom * (1d0 + delta*Mach_N**2)**(Gamma/(Gamma-1d0))
      if (pstar >= p_Amb) p_exit = pstar
      Isentropic    = (Gamma - 1d0) / (2d0*Gamma)
      Face_Sound    = Sound * ( (p_exit/Prim(np))**Isentropic )
      Face_Un       = (Riem - Face_Sound) / (BC_Sign*0.5d0*(Gamma-1d0))
      Face_Prim(np) = p_exit
      Face_rho      = Gamma * Face_Prim(np) / Face_Sound**2
      do s = 1, nsc
        Face_Prim(s) = Prim(s) / rho * Face_rho
      end do
      Face_Prim(nu) = Vel_T(1) + Face_Un*Normal(1)
      Face_Prim(nv) = Vel_T(2) + Face_Un*Normal(2)
      Face_Prim(nw) = Vel_T(3) + Face_Un*Normal(3)
      if (nprim > np) then
        do s = np+1, nprim
          Face_Prim(s) = Prim(s) / rho * Face_rho
        end do
      end if
      Face_Rgas = Rgas
    else
      do s = 1, nprim
        Face_Prim(s) = Prim(s)
      end do
      Face_Un = Un; Face_Sound = Sound; Face_rho = rho; Face_Rgas = Rgas
    end if

    Face_Enthalpy = 0d0
    Face_T        = Face_Prim(np) / (Face_Rgas*Face_rho)
    T_i = idint(Face_T); Tdiff = Face_T - T_i; Tint(1) = T_i; Tint(2) = T_i + 1
    do s = 1, nprim
      Flux(s) = 0d0
    end do
    do s = 1, nsc
      Face_Enthalpy = Face_Enthalpy + Face_Prim(s)/Face_rho * f_tabT_expr( s, h_tab, Tint, Tdiff )
      Flux(s)       = Face_Prim(s) * Area * Face_Un
    end do
    Fmass = Face_rho * Area * Face_Un
    Flux(nu) = Face_Prim(nu)*Fmass + Face_Prim(np)*Area*Normal(1)
    Flux(nv) = Face_Prim(nv)*Fmass + Face_Prim(np)*Area*Normal(2)
    Flux(nw) = Face_Prim(nw)*Fmass + Face_Prim(np)*Area*Normal(3)
    vt2      = Face_Prim(nu)**2 + Face_Prim(nv)**2 + Face_Prim(nw)**2
    Flux(np) = Fmass * (0.5d0*vt2 + Face_Enthalpy)
    if (nprim > np) then
      do s = np+1, nprim
        Flux(s) = Face_Prim(s)/Face_rho * Fmass
      end do
    end if

  end subroutine BC_Outflow


  !> Device twin of BC_Supersonic_Inflow (devirt riemann_HLLC + rlimiter_MINMOD).
  subroutine BC_Supersonic_Inflow ( Prim2, Int_Prim, modfm, modfm1, modfm2, BC_Mach, BC_p0, BC_T0, &
                                        BC_ci, BC_RANS, BC_rel_fac, BC_alpha, BC_beta, Normal, Area, Flux, Fmass )
    !$acc routine seq
    use FLINT_Lib_Thermodynamic,       only: f_Rtot, Ri_tab
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, f_ss_dev
    use MOSE_Lib_Riemann,              only: Riemann, rs_dev
    use MOSE_Lib_Limiters,             only: rlimiter_MINMOD, rlimiter, lim_dev

    implicit none
    integer,  intent(in)  :: modfm, modfm1, modfm2
    real(R8), intent(in)  :: Prim2(NPRIM_DEV), Int_Prim(NPRIM_DEV), BC_Mach, BC_p0, BC_T0, BC_rel_fac
    real(R8), intent(in)  :: BC_ci(nsc), BC_RANS(nrans), BC_alpha, BC_beta, Normal(3), Area
    real(R8), intent(out) :: Flux(NPRIM_DEV), Fmass
    integer  :: s
    real(R8) :: Sup_Prim(NPRIM_DEV), Sup_T, Sup_Rgas, Sup_rho, Sup_Sound, Prim1(NPRIM_DEV), Prim3(NPRIM_DEV)
    real(R8) :: Slope, Face_Prim(NPRIM_DEV), Prim_L(NPRIM_DEV), Prim_R(NPRIM_DEV)
    real(R8) :: rho_L, rho_R, Rgas_L, Rgas_R, Sound_L, Sound_R, F_r, F_u, F_v, F_w, F_E
    real(R8) :: su, Sel_L, Sel_R, vmag2, flow_Rgas, flow_Sound, flow_Mach, BC_Mach_local, BC_p0_local

    flow_Rgas  = f_Rtot( Prim2(1:nsc) )
    vmag2      = Prim2(nu)**2 + Prim2(nv)**2 + Prim2(nw)**2
    flow_Sound = f_ss_dev ( Prim2(1:nsc), Prim2(np), sum(Prim2(1:nsc)), flow_Rgas )
    flow_Mach  = sqrt(vmag2) / flow_Sound

    BC_Mach_local = BC_rel_fac*BC_Mach + (1d0-BC_rel_fac)*flow_Mach
    BC_p0_local   = BC_rel_fac*BC_p0   + (1d0-BC_rel_fac)*Prim2(np)

    Sup_Prim(np) = BC_p0_local
    Sup_T        = BC_T0
    Sup_Rgas     = 0d0
    do s = 1, nsc
      Sup_Rgas = Sup_Rgas + BC_ci(s)*Ri_tab(s)
    end do
    Sup_rho  = Sup_Prim(np) / (Sup_Rgas*Sup_T)
    do s = 1, nsc
      Sup_Prim(s) = Sup_rho * BC_ci(s)
    end do
    Sup_Sound    = f_ss_dev ( Sup_Prim(1:nsc), Sup_Prim(np), Sup_rho, Sup_Rgas )
    Sup_Prim(nu) = BC_Mach_local * Sup_Sound * cos(BC_alpha) * cos(BC_beta)
    Sup_Prim(nv) = BC_Mach_local * Sup_Sound * sin(BC_alpha) * cos(BC_beta)
    Sup_Prim(nw) = BC_Mach_local * Sup_Sound * sin(BC_beta)
    if (RANS) then
      do s = 1, nrans
        Sup_Prim(nt-1+s) = Sup_rho * BC_RANS(s)
      end do
    end if

    do s = 1, nprim
      Prim1(s) = modfm*Sup_Prim(s) + modfm1*Int_Prim(s)
      Prim3(s) = modfm*Int_Prim(s) + modfm1*Sup_Prim(s)
      Slope    = rlimiter ( lim_dev,  Prim3(s)-Prim2(s), Prim2(s)-Prim1(s) )
      Face_Prim(s) = Prim2(s) + 0.5d0*Slope*modfm2
      Prim_L(s) = modfm*Sup_Prim(s)  + modfm1*Face_Prim(s)
      Prim_R(s) = modfm*Face_Prim(s) + modfm1*Sup_Prim(s)
    end do

    call co_rotot_Rtot_dev ( Prim_L(1:nsc), rho_L, Rgas_L )
    Sound_L = f_ss_dev ( Prim_L(1:nsc), Prim_L(np), rho_L, Rgas_L )
    call co_rotot_Rtot_dev ( Prim_R(1:nsc), rho_R, Rgas_R )
    Sound_R = f_ss_dev ( Prim_R(1:nsc), Prim_R(np), rho_R, Rgas_R )

    call Riemann ( rs_dev,  Prim_L(1:nsc), Prim_L(nu), Prim_L(nv), Prim_L(nw), Prim_L(np), Sound_L, rho_L, &
                        Prim_R(1:nsc), Prim_R(nu), Prim_R(nv), Prim_R(nw), Prim_R(np), Sound_R, rho_R, &
                        1d0, Sound_L, Sound_R, Normal(1), Normal(2), Normal(3), F_r, F_u, F_v, F_w, F_E )   ! Ur = a: no preconditioning (MOSE_open)

    su    = sign ( 1d0, F_r )
    Sel_L = (1d0 + su) / 2d0
    Sel_R = (su - 1d0) / 2d0

    Fmass = F_r * Area
    do s = 1, nprim
      Flux(s) = 0d0
    end do
    do s = 1, nsc
      Flux(s) = Fmass * ( Sel_L*Prim_L(s)/rho_L - Sel_R*Prim_R(s)/rho_R )
    end do
    Flux(nu) = F_u * Area; Flux(nv) = F_v * Area; Flux(nw) = F_w * Area; Flux(np) = F_E * Area
    if (RANS) then
      do s = nt, nprim
        Flux(s) = Fmass * ( Sel_L*Prim_L(s)/rho_L - Sel_R*Prim_R(s)/rho_R )
      end do
    end if

  end subroutine BC_Supersonic_Inflow


  !> Device twin of BC_Subsonic_Inflow. error=1 signals the caller to fall back
  !> to a reflective wall (symmetry) for this face/iteration.
  subroutine BC_Subsonic_Inflow ( Prim, Un, Sound, Gamma, Riem, BC_alpha, BC_beta, BC_Mach, &
                                      BC_ci, BC_RANS, BC_T0, BC_p0, BC_rel_fac, BC_Sign, t_Vec, &
                                      Normal, Area, Flux, Fmass, error )
    !$acc routine seq
    use FLINT_Lib_Thermodynamic, only: f_tabT_expr, h_tab, cp_tab, Ri_tab

    implicit none
    real(R8), intent(in)  :: Prim(NPRIM_DEV), Un, Sound, Gamma, Riem, BC_alpha, BC_beta, BC_Mach
    real(R8), intent(in)  :: BC_ci(nsc), BC_RANS(nrans), BC_T0, BC_p0, BC_rel_fac, BC_Sign
    real(R8), intent(in)  :: t_Vec(3), Normal(3), Area
    real(R8), intent(out) :: Flux(NPRIM_DEV), Fmass
    integer,  intent(out) :: error
    integer  :: s, T_i, Tint(2)
    integer  :: iter_discr
    real(R8) :: Tdiff
    real(R8) :: b_Vec(3), b_Mod, alpha, beta, XA, XB, XC, XD, XE, XF, check, E1, E2
    real(R8) :: Face_Rgas, Sound0, Coef_D, Coef_C, YA, YB, YC, Discr, Un_Sign, Un3, Sound3
    real(R8) :: p3, Ut3, Ub3, mdot, Rgas3, cp3, T3, rho3, BC_p0_try
    real(R8) :: Face_T, Face_Sound, Face_Prim(NPRIM_DEV), Face_rho, Face_Enthalpy, Face_Un, vt2

    error = 0

    ! b = n x t  (then normalised)
    b_Vec(1) = Normal(2)*t_Vec(3) - Normal(3)*t_Vec(2)
    b_Vec(2) = Normal(3)*t_Vec(1) - Normal(1)*t_Vec(3)
    b_Vec(3) = Normal(1)*t_Vec(2) - Normal(2)*t_Vec(1)
    b_Mod    = sqrt( b_Vec(1)**2 + b_Vec(2)**2 + b_Vec(3)**2 )
    b_Vec(1) = b_Vec(1)/b_Mod; b_Vec(2) = b_Vec(2)/b_Mod; b_Vec(3) = b_Vec(3)/b_Mod

    if ( (BC_alpha == 0d0) .and. (BC_beta == 0d0) ) then
      alpha = atan ( Normal(2)/(Normal(1) + 1d-20) )
      beta  = atan ( Normal(3)/(Normal(1) + 1d-20) )
    else
      alpha = BC_alpha; beta = BC_beta
    end if

    XA = Normal(2) - Normal(1)*tan(alpha)
    XB =  t_Vec(2) -  t_Vec(1)*tan(alpha)
    XC =  b_Vec(2) -  b_Vec(1)*tan(alpha)
    XD = Normal(3) - Normal(1)*tan(beta)
    XE =  t_Vec(3) -  t_Vec(1)*tan(beta)
    XF =  b_Vec(3) -  b_Vec(1)*tan(beta)
    check = (XC*XE - XF*XB)
    E1 = (XF*XA - XC*XD) / check
    E2 = (XD*XB - XE*XA) / check

    if ( BC_Mach >= 0d0 ) then
      Discr     = -1d0
      BC_p0_try = BC_p0
      iter_discr = 0
      do while (Discr < 0d0)
        iter_discr = iter_discr + 1
        Face_Rgas = 0d0
        do s = 1, nsc
          Face_Rgas = Face_Rgas + BC_ci(s)*Ri_tab(s)
        end do
        Sound0 = sqrt ( Gamma*Face_Rgas*BC_T0 )
        Coef_D = (1d0 + E1**2 + E2**2) * (Gamma - 1d0)/2d0
        Coef_C = (BC_p0_try/Prim(np))**((Gamma-1d0)/Gamma)
        Coef_C = Coef_C * ((Sound/Sound0)**2)
        YA = ((Gamma-1d0)/2d0)**2 + Coef_C*Coef_D
        YB = BC_Sign*Riem*((Gamma-1d0)/2d0)
        YC = Riem**2 - Coef_C*Sound0**2
        Discr = YB**2 - YA*YC
        if ( (Discr < 0d0) .and. (abs(Discr) < 1d-2) ) Discr = 0d0
        if ( (Discr < 0d0) .and. (abs(Discr) > 1d-2) ) then
          ! The imposed (p0, T0) cannot drive inflow against this cell pressure. As MOSE_open
          ! (Lib_BC_Fluxes_Inflow.f90:505-524): relax p0 towards the cell pressure, at most 1000 times, then take
          ! Discr = 0. With rel_fac = 1 the relaxation changes nothing and MOSE_open's 1000 tries end on the state of
          ! the first pass: taken at once (before, the face was closed: a symmetry for this step; A4 of 2026-10-16).
          if ( BC_rel_fac == 1d0 ) then
            Discr = 0d0
            exit
          end if
          BC_p0_try = BC_rel_fac*BC_p0_try + (1d0-BC_rel_fac)*Prim(np)
        end if
        if ( iter_discr > 1000 ) then
          Discr = 0d0
          exit
        end if
      end do

      Un_Sign = sign ( 1d0, Un )
      Un3     = (YB + Un_Sign*sqrt(Discr)) / YA
      Sound3  = sqrt ( Sound0**2 - Coef_D*Un3**2 )
      p3      = BC_p0_try * (Sound3/Sound0)**(2d0*Gamma/(Gamma-1d0))
      Ut3     = Un3*E1
      Ub3     = Un3*E2
      Face_Sound    = Sound3
      Face_Prim(np) = p3
      Face_rho      = Gamma*Face_Prim(np)/Face_Sound**2
      do s = 1, nsc
        Face_Prim(s) = BC_ci(s)*Face_rho
      end do
      Face_Prim(nu) = Un3*Normal(1) + Ut3*t_Vec(1) + Ub3*b_Vec(1)
      Face_Prim(nv) = Un3*Normal(2) + Ut3*t_Vec(2) + Ub3*b_Vec(2)
      Face_Prim(nw) = Un3*Normal(3) + Ut3*t_Vec(3) + Ub3*b_Vec(3)
      Face_Un       = Un3
      if (RANS) then
        do s = 1, nrans
          Face_Prim(nt-1+s) = BC_RANS(s)*Face_rho
        end do
      end if

    else  ! BC_Mach < 0 : assigned mass flux + temperature
      mdot  = BC_p0
      p3    = Prim(np)
      Rgas3 = 0d0
      cp3   = 0d0
      T_i = idint(BC_T0); Tdiff = BC_T0 - T_i; Tint(1) = T_i; Tint(2) = T_i + 1
      do s = 1, nsc
        cp3   = cp3   + BC_ci(s)*f_tabT_expr( s, cp_tab, Tint, Tdiff )
        Rgas3 = Rgas3 + BC_ci(s)*Ri_tab(s)
      end do

      if ( BC_Mach /= -15d0 ) then
        if ( BC_Mach == -10d0 ) then
          T3 = ( -cp3 + sqrt( cp3**2 + 2d0*cp3*BC_T0*(mdot*Rgas3/p3)**2 ) ) &
               * ( p3/(mdot*Rgas3) )**2
        else
          T3 = BC_T0
        end if
        rho3      = p3 / (Rgas3*T3)
        Un3       = mdot / rho3
        Ut3       = Un3*E1
        Ub3       = Un3*E2
        Face_Rgas = Rgas3
        Face_rho  = rho3
        do s = 1, nsc
          Face_Prim(s) = Face_rho*BC_ci(s)
        end do
        Face_Prim(nu) = Un3*Normal(1) + Ut3*t_Vec(1) + Ub3*b_Vec(1)
        Face_Prim(nv) = Un3*Normal(2) + Ut3*t_Vec(2) + Ub3*b_Vec(2)
        Face_Prim(nw) = Un3*Normal(3) + Ut3*t_Vec(3) + Ub3*b_Vec(3)
        Face_Prim(np) = p3
        Face_Un       = Un3
      else  ! HRE tangential injector
        Face_Rgas = Rgas3
        Face_rho  = p3/(Face_Rgas*BC_T0)
        do s = 1, nsc
          Face_Prim(s) = Face_rho*BC_ci(s)
        end do
        Face_Prim(nu) = mdot/Face_rho
        Face_Prim(nv) = 0d0
        Face_Prim(nw) = Face_Prim(nu)*0.25d0
        Face_Prim(np) = p3
        Face_Un       = Face_Prim(nu)*Normal(1) + Face_Prim(nv)*Normal(2) + Face_Prim(nw)*Normal(3)
      end if
      if (RANS) then
        do s = 1, nrans
          Face_Prim(nt-1+s) = BC_RANS(s)*Face_rho
        end do
      end if
    end if

    Face_Enthalpy = 0d0
    Face_T        = Face_Prim(np) / (Face_Rgas*Face_rho)
    T_i = idint(Face_T); Tdiff = Face_T - T_i; Tint(1) = T_i; Tint(2) = T_i + 1
    do s = 1, nprim
      Flux(s) = 0d0
    end do
    do s = 1, nsc
      Face_Enthalpy = Face_Enthalpy + Face_Prim(s)/Face_rho * f_tabT_expr( s, h_tab, Tint, Tdiff )
      Flux(s)       = Face_Prim(s) * Face_Un * Area
    end do
    Fmass    = Face_rho * Face_Un * Area
    Flux(nu) = Fmass*Face_Prim(nu) + Face_Prim(np)*Area*Normal(1)
    Flux(nv) = Fmass*Face_Prim(nv) + Face_Prim(np)*Area*Normal(2)
    Flux(nw) = Fmass*Face_Prim(nw) + Face_Prim(np)*Area*Normal(3)
    vt2      = Face_Prim(nu)**2 + Face_Prim(nv)**2 + Face_Prim(nw)**2
    Flux(np) = Fmass * (0.5d0*vt2 + Face_Enthalpy)
    if (RANS) then
      do s = nt, nprim
        Flux(s) = Face_Prim(s)/Face_rho * Fmass
      end do
    end if

  end subroutine BC_Subsonic_Inflow


  !> Device twin of BC_Inflow_Outflow (no rotating frame — guarded at setup).
  !> Returns error=1 when the inlet must close to a reflective wall; the caller
  !> (BC_Fluxes_blk) then applies BC_Symmetry_Eul. The mdot diagnostic
  !> is not written back on the device path (it does not enter the residual).
  subroutine BC_Inflow_Outflow ( PWARG, r, n1, n2, n3, Im, Jm, Km, Fm, Normal, Area, t_Vec, &
                                     BC_Mach, BC_T0, BC_p0, BC_rel_fac, BC_alpha, BC_beta, BC_pAmb, &
                                     BC_ci, BC_RANS, error )
    !$acc routine seq
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, f_ss_dev, f_cp_dev

    implicit none
    integer,  intent(in)    :: n1, n2, n3, Im, Jm, Km, Fm
    real(R8), intent(in)    :: Normal(3), Area, t_Vec(3)
    real(R8), intent(in)    :: BC_Mach, BC_T0, BC_p0, BC_rel_fac, BC_alpha, BC_beta, BC_pAmb
    real(R8), intent(in)    :: BC_ci(nsc), BC_RANS(nrans)
    integer,  intent(in)    :: ncw, nrw, boff, pi, pj, wroff
    integer,  intent(in)    :: nur
    real(R8), intent(in)    :: ur(nur)   ! Ur of the cells (prec; 1 element otherwise): read by HLLC-PC only
    real(R8), intent(in)    :: w(ncw, nprim)
    real(R8), intent(inout) :: wr(nrw, nprim)   ! Res var-last: coalesced residual WRITE target
    real(R8), intent(inout) :: r(nprim, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    integer,  intent(out)   :: error
    integer  :: modfm, modfm1, modfm2, gstep, Int_i, Int_j, Int_k, dir, s, Tint(2)
    real(R8) :: Bound_Prim(NPRIM_DEV), Int_Prim(NPRIM_DEV), Un, BC_Sign
    real(R8) :: Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem, InOut, Flux(NPRIM_DEV), Fmass
    real(R8) :: T, Tdiff, cp

    error  = 0
    dir    = (Fm + 1)/2
    modfm  = mod(Fm, 2)
    modfm1 = 1 - modfm
    modfm2 = 1 - 2*modfm
    gstep  = 2*modfm - 1
    Int_i = Im; Int_j = Jm; Int_k = Km
    select case (dir)
      case (1); Int_i = Im + gstep
      case (2); Int_j = Jm + gstep
      case (3); Int_k = Km + gstep
    end select

    do s = 1, nprim
      Bound_Prim(s) = PCELL(s, Im, Jm, Km)
      Int_Prim(s)   = PCELL(s, Int_i, Int_j, Int_k)
    end do
    BC_Sign = real(modfm2, R8)

    Un = Bound_Prim(nu)*Normal(1) + Bound_Prim(nv)*Normal(2) + Bound_Prim(nw)*Normal(3)
    call co_rotot_Rtot_dev ( Bound_Prim(1:nsc), Bound_rho, Bound_Rgas )
    Bound_Sound = f_ss_dev ( Bound_Prim(1:nsc), Bound_Prim(np), Bound_rho, Bound_Rgas )
    T       = Bound_Prim(np) / (Bound_rho*Bound_Rgas)
    Tint(1) = int(T); Tint(2) = Tint(1) + 1; Tdiff = T - Tint(1)
    cp          = f_cp_dev ( Bound_Prim(1:nsc), Tint, Tdiff, Bound_rho )
    Bound_Gamma = cp / (cp - Bound_Rgas)

    Riem  = Bound_Sound + BC_Sign*0.5d0*(Bound_Gamma-1d0)*Un
    InOut = BC_Sign*Un
    if (BC_Mach < 0d0) InOut = -1d0

    if (InOut >= 0d0) then
      if (BC_pAmb == -3.14d0) then
        error = 1
        return
      else
        call BC_Outflow ( Bound_Prim, Bound_rho, Bound_Rgas, Bound_Sound, Un, Bound_Gamma, &
                              Riem, Normal, Area, BC_pAmb, BC_rel_fac, BC_Sign, Flux, Fmass )
      end if
    else
      if (BC_Mach > 1d0) then
        call BC_Supersonic_Inflow ( Bound_Prim, Int_Prim, modfm, modfm1, modfm2, BC_Mach, BC_p0, BC_T0, &
                                        BC_ci, BC_RANS, BC_rel_fac, BC_alpha, BC_beta, Normal, Area, Flux, Fmass )
      else if (BC_Mach < 1d0) then
        call BC_Subsonic_Inflow ( Bound_Prim, Un, Bound_Sound, Bound_Gamma, Riem, BC_alpha, BC_beta, BC_Mach, &
                                      BC_ci, BC_RANS, BC_T0, BC_p0, BC_rel_fac, BC_Sign, t_Vec, Normal, Area, &
                                      Flux, Fmass, error )
        if (error == 1) return
      end if
    end if

    do s = 1, nprim
      WRCELL(s,Im,Jm,Km) = WRCELL(s,Im,Jm,Km) + modfm2 * Flux(s)
    end do

  end subroutine BC_Inflow_Outflow


  !> Device twin of BC_Outlet (MOSE_open BC 406 = internal type 24): back-pressure outlet, the relaxed exit
  !> pressure p_exit = rel_fac*p_back + (1-rel_fac)*p(cell) through BC_Outflow; a reversed (inward) flow
  !> closes the face (error = 1 -> the caller applies the symmetry flux, as the host BC_Fluxes does). Mirror of
  !> Lib_BC BC_Outlet; Fmass is not accumulated on the device (host bookkeeping only), as in BC_Inflow_Outflow.
  subroutine BC_Outlet ( PWARG, r, n1, n2, n3, Im, Jm, Km, Fm, Normal, Area, BC_pExit, BC_rel_fac, error )
    !$acc routine seq
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, f_ss_dev, f_cp_dev

    implicit none
    integer,  intent(in)    :: n1, n2, n3, Im, Jm, Km, Fm
    real(R8), intent(in)    :: Normal(3), Area, BC_pExit, BC_rel_fac
    integer,  intent(in)    :: ncw, nrw, boff, pi, pj, wroff
    integer,  intent(in)    :: nur
    real(R8), intent(in)    :: ur(nur)   ! Ur of the cells (prec; 1 element otherwise): read by HLLC-PC only
    real(R8), intent(in)    :: w(ncw, nprim)
    real(R8), intent(inout) :: wr(nrw, nprim)   ! Res var-last: coalesced residual WRITE target
    real(R8), intent(inout) :: r(nprim, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    integer,  intent(out)   :: error
    integer  :: modfm, modfm2, s, Tint(2)
    real(R8) :: Bound_Prim(NPRIM_DEV), Un, BC_Sign
    real(R8) :: Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem, p_exit, Flux(NPRIM_DEV), Fmass
    real(R8) :: T, Tdiff, cp

    error  = 0
    modfm  = mod(Fm, 2)
    modfm2 = 1 - 2*modfm
    do s = 1, nprim
      Bound_Prim(s) = PCELL(s, Im, Jm, Km)
    end do
    BC_Sign = real(modfm2, R8)                        ! -1 on faces 1,3,5 ; +1 on faces 2,4,6
    Un = Bound_Prim(nu)*Normal(1) + Bound_Prim(nv)*Normal(2) + Bound_Prim(nw)*Normal(3)

    ! Reversed (inward) flow: this BC has no inlet closure, close the face.
    if ( BC_Sign * Un < 0d0 ) then
      error = 1
      return
    end if

    call co_rotot_Rtot_dev ( Bound_Prim(1:nsc), Bound_rho, Bound_Rgas )
    Bound_Sound = f_ss_dev ( Bound_Prim(1:nsc), Bound_Prim(np), Bound_rho, Bound_Rgas )
    T       = Bound_Prim(np) / (Bound_rho*Bound_Rgas)
    Tint(1) = int(T); Tint(2) = Tint(1) + 1; Tdiff = T - Tint(1)
    cp          = f_cp_dev ( Bound_Prim(1:nsc), Tint, Tdiff, Bound_rho )
    Bound_Gamma = cp / (cp - Bound_Rgas)

    ! riemann invariant a+-delta*u directed towards the interface
    Riem = Bound_Sound + BC_Sign * 0.5d0*( Bound_Gamma - 1d0 ) * Un

    ! Relaxed exit pressure: partially non-reflecting outlet.
    p_exit = BC_rel_fac * BC_pExit + ( 1d0 - BC_rel_fac ) * Bound_Prim(np)

    call BC_Outflow ( Bound_Prim, Bound_rho, Bound_Rgas, Bound_Sound, Un, Bound_Gamma, &
                          Riem, Normal, Area, p_exit, BC_rel_fac, BC_Sign, Flux, Fmass )

    do s = 1, nprim
      WRCELL(s,Im,Jm,Km) = WRCELL(s,Im,Jm,Km) + modfm2 * Flux(s)
    end do

  end subroutine BC_Outlet


  !> Device twin of BC_Connection_Eul. The donor block's solution already sits in
  !> the receiver ghost cells (filled by the host Ghost_Connection before
  !> BC_Fluxes), so this reads only the receiver `P` (+ flat dl) — no donor access.
  !> 4-point stencil straddling the face; HLLC at the interface. Updates r(:).
  subroutine BC_Connection_Eul ( PWARG, r, beta_c, dlc4, n1, n2, n3, Im, Jm, Km, Fm, Normal, Area )
    !$acc routine seq
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, f_ss_dev
    use MOSE_Lib_Riemann,              only: Riemann, rs_dev
    use MOSE_Lib_Shock_Detector,       only: sd_dev, SD_beta
    use MOSE_Lib_Limiters,             only: rlimiter_MINMOD, SD_space_order, rlimiter, lim_dev

    implicit none
    integer,  intent(in)    :: n1, n2, n3, Im, Jm, Km, Fm
    real(R8), intent(in)    :: Normal(3), Area, beta_c
    integer,  intent(in)    :: ncw, nrw, boff, pi, pj, wroff
    integer,  intent(in)    :: nur
    real(R8), intent(in)    :: ur(nur)   ! Ur of the cells (prec; 1 element otherwise): read by HLLC-PC only
    real(R8), intent(in)    :: w(ncw, nprim)
    real(R8), intent(inout) :: wr(nrw, nprim)   ! Res var-last: coalesced residual WRITE target
    real(R8), intent(in)    :: dlc4(4)   !< 5r-b: dl_f(dir,.) at the four reconstruction points
    real(R8), intent(inout) :: r(nprim, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    integer  :: dir, modfm, modfm1, modfm2, gstep, s, o1, o2, o3, o4
    integer  :: i1,j1,k1, i2,j2,k2, i3,j3,k3, i4,j4,k4
    real(R8) :: Prim1(NPRIM_DEV), Prim2(NPRIM_DEV), Prim3(NPRIM_DEV), Prim4(NPRIM_DEV)
    real(R8) :: Prim_L(NPRIM_DEV), Prim_R(NPRIM_DEV)
    real(R8) :: Dl1, Dl2, Dl3, Dl4, Dist12, Dist23, Dist34, beta_, d21, d32, d43
    real(R8) :: rho_L, rho_R, Rtot_L, Rtot_R, a_L, a_R, F_r, F_u, F_v, F_w, F_E, Fmass, su, sel_L, sel_R, urL_, urR_

    dir    = (Fm + 1)/2
    modfm  = mod(Fm, 2)
    modfm1 = 1 - modfm
    modfm2 = 1 - 2*modfm
    gstep  = 2*modfm - 1

    ! dir-axis offsets of the 4 stencil points (mirror the host modfm/guide algebra)
    o1 = modfm*(-gstep*gc) + modfm1*( gstep)
    o2 = modfm*(-gstep)
    o3 = modfm1*(-gstep)
    o4 = modfm1*(-gstep*gc) + modfm*( gstep)

    i1=Im;j1=Jm;k1=Km; i2=Im;j2=Jm;k2=Km; i3=Im;j3=Jm;k3=Km; i4=Im;j4=Jm;k4=Km
    select case (dir)
      case (1); i1=Im+o1; i2=Im+o2; i3=Im+o3; i4=Im+o4
      case (2); j1=Jm+o1; j2=Jm+o2; j3=Jm+o3; j4=Jm+o4
      case (3); k1=Km+o1; k2=Km+o2; k3=Km+o3; k4=Km+o4
    end select

    do s = 1, nprim
      Prim1(s) = PCELL(s,i1,j1,k1); Prim2(s) = PCELL(s,i2,j2,k2)
      Prim3(s) = PCELL(s,i3,j3,k3); Prim4(s) = PCELL(s,i4,j4,k4)
    end do

    if (sd_dev > 0) then   ! the MOSE_open sensor of the record's cell (blk%beta(Im,Jm,Km)) from its six face neighbours
      beta_ = SD_beta ( sd_dev, PCELL(np,Im,Jm,Km), PCELL(np,Im-1,Jm,Km), PCELL(np,Im+1,Jm,Km), &
                            PCELL(np,Im,Jm-1,Km), PCELL(np,Im,Jm+1,Km), PCELL(np,Im,Jm,Km-1), PCELL(np,Im,Jm,Km+1) )
    else if (SD_space_order) then
      beta_ = beta_c
    else
      beta_ = 1d0
    end if

    Dl1 = dlc4(1)*0.5d0; Dl2 = dlc4(2)*0.5d0
    Dl3 = dlc4(3)*0.5d0; Dl4 = dlc4(4)*0.5d0
    Dist12 = Dl1 + Dl2; Dist23 = Dl2 + Dl3; Dist34 = Dl3 + Dl4

    do s = 1, nprim
      d21 = (Prim2(s) - Prim1(s)) / Dist12
      d32 = (Prim3(s) - Prim2(s)) / Dist23
      d43 = (Prim4(s) - Prim3(s)) / Dist34
      Prim_L(s) = Prim2(s) + beta_ * rlimiter ( lim_dev, d32, d21) * Dl2
      Prim_R(s) = Prim3(s) - beta_ * rlimiter ( lim_dev, d43, d32) * Dl3
    end do

    call co_rotot_Rtot_dev ( Prim_L(1:nsc), rho_L, Rtot_L )
    a_L = f_ss_dev ( Prim_L(1:nsc), Prim_L(np), rho_L, Rtot_L )
    call co_rotot_Rtot_dev ( Prim_R(1:nsc), rho_R, Rtot_R )
    a_R = f_ss_dev ( Prim_R(1:nsc), Prim_R(np), rho_R, Rtot_R )

    urL_ = 0d0; urR_ = 0d0
    if (rs_dev == 14) then   ! HLLC-PC: MOSE_open blk%Ur(I2,J2,K2), blk%Ur(I3,J3,K3); the ghost one is 0 (never filled)
      urL_ = URCELL(i2,j2,k2); urR_ = URCELL(i3,j3,k3)
    end if
    call Riemann ( rs_dev,  Prim_L(1:nsc), Prim_L(nu), Prim_L(nv), Prim_L(nw), Prim_L(np), a_L, rho_L, &
                        Prim_R(1:nsc), Prim_R(nu), Prim_R(nv), Prim_R(nw), Prim_R(np), a_R, rho_R, &
                        beta_, urL_, urR_, Normal(1), Normal(2), Normal(3), F_r, F_u, F_v, F_w, F_E )   ! MOSE_open: the sensor of the cell

    su    = sign ( 0.5d0, F_r )
    sel_L = 0.5d0 + su
    sel_R = su - 0.5d0
    Fmass = F_r * Area

    do s = 1, nsc
      WRCELL(s,Im,Jm,Km) = WRCELL(s,Im,Jm,Km) + modfm2 * Fmass*( sel_L*Prim_L(s)/rho_L - sel_R*Prim_R(s)/rho_R )
    end do
    WRCELL(nu,Im,Jm,Km) = WRCELL(nu,Im,Jm,Km) + modfm2 * F_u*Area
    WRCELL(nv,Im,Jm,Km) = WRCELL(nv,Im,Jm,Km) + modfm2 * F_v*Area
    WRCELL(nw,Im,Jm,Km) = WRCELL(nw,Im,Jm,Km) + modfm2 * F_w*Area
    WRCELL(np,Im,Jm,Km) = WRCELL(np,Im,Jm,Km) + modfm2 * F_E*Area
    if (nprim > np) then
      do s = np+1, nprim
        WRCELL(s,Im,Jm,Km) = WRCELL(s,Im,Jm,Km) + modfm2 * Fmass*( sel_L*Prim_L(s)/rho_L - sel_R*Prim_R(s)/rho_R )
      end do
    end if

  end subroutine BC_Connection_Eul


  !> Device helper: primitive -> "viscous" variables (Y_s, vel, T, RANS/rho).
  !> Mirror of the Visc_Variables contained subroutine in BC_Connection_Visc.
  subroutine Visc_Variables ( Prim, Visc )
    !$acc routine seq
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev

    implicit none
    real(R8), intent(in)  :: Prim(NPRIM_DEV)
    real(R8), intent(out) :: Visc(NPRIM_DEV)
    integer  :: s
    real(R8) :: rho, Rgas

    call co_rotot_Rtot_dev ( Prim(1:nsc), rho, Rgas )
    do s = 1, nsc
      Visc(s) = Prim(s) / rho
    end do
    Visc(nu) = Prim(nu); Visc(nv) = Prim(nv); Visc(nw) = Prim(nw)
    Visc(np) = Prim(np) / (rho*Rgas)
    if (nprim > np) then
      do s = np+1, nprim
        Visc(s) = Prim(s) / rho
      end do
    end if

  end subroutine Visc_Variables


  !> Device twin of BC_Connection_Visc (SST). The donor stencil arrives in Pg(:,1)
  !> (adjacent ghost) and Pg(:,3:6) (ghost-side tangential neighbours); receiver
  !> interior neighbours come from P. Two Compute_Diffusive_Flux calls (the shared
  !> device viscous kernel) on the local- and ghost-side gradients. Updates r(:).
  subroutine BC_Connection_Visc ( PWARG, r, Mrec, Mgho, ynrec, yngho, Pg, n1, n2, n3, Im, Jm, Km, Fm, Normal, Area )
    !$acc routine seq
    use MOSE_Mod_Diffusive, only: Compute_Diffusive_Flux

    implicit none
    integer,  intent(in)    :: n1, n2, n3, Im, Jm, Km, Fm
    real(R8), intent(in)    :: Normal(3), Area
    integer,  intent(in)    :: ncw, nrw, boff, pi, pj, wroff
    integer,  intent(in)    :: nur
    real(R8), intent(in)    :: ur(nur)   ! Ur of the cells (prec; 1 element otherwise): read by HLLC-PC only
    real(R8), intent(in)    :: w(ncw, nprim)
    real(R8), intent(inout) :: wr(nrw, nprim)   ! Res var-last: coalesced residual WRITE target
    real(R8), intent(inout) :: r(nprim, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    real(R8), intent(in)    :: Mrec(3,3), Mgho(3,3)   !< 5r-b: m_f at the record's cell and at the ghost
    real(R8), intent(in)    :: ynrec, yngho
    real(R8), intent(in)    :: Pg(nprim, 6)
    integer  :: dir, modfm2, gstep, Ig, Jg, Kg, s, kk
    real(R8) :: M(3,3), Waldis
    real(R8) :: Prim_loc(NPRIM_DEV), Prim_ghost(NPRIM_DEV), Prim_avg(NPRIM_DEV), tmp(NPRIM_DEV)
    real(R8) :: Visc_loc(NPRIM_DEV), Visc_ghost(NPRIM_DEV)
    real(R8) :: Visc_ip(NPRIM_DEV), Visc_im(NPRIM_DEV), Visc_jp(NPRIM_DEV), Visc_jm(NPRIM_DEV)
    real(R8) :: Visc_kp(NPRIM_DEV), Visc_km(NPRIM_DEV)
    real(R8) :: Visc_g3(NPRIM_DEV), Visc_g4(NPRIM_DEV), Visc_g5(NPRIM_DEV), Visc_g6(NPRIM_DEV)
    real(R8) :: Gloc(NPRIM_DEV,3), Gghost(NPRIM_DEV,3), G2loc(NPRIM_DEV,3), G2ghost(NPRIM_DEV,3)
    real(R8) :: Flux_loc(NPRIM_DEV), Flux_ghost(NPRIM_DEV)

    dir    = (Fm + 1)/2
    modfm2 = 1 - 2*mod(Fm,2)
    gstep  = 2*mod(Fm,2) - 1
    Ig = Im; Jg = Jm; Kg = Km
    select case (dir)
      case (1); Ig = Im - gstep
      case (2); Jg = Jm - gstep
      case (3); Kg = Km - gstep
    end select

    do s = 1, 3
      M(1,s) = 0.5d0*( Mrec(1,s) + Mgho(1,s) )
      M(2,s) = 0.5d0*( Mrec(2,s) + Mgho(2,s) )
      M(3,s) = 0.5d0*( Mrec(3,s) + Mgho(3,s) )
    end do
    Waldis = 0.5d0*( ynrec + yngho )

    do s = 1, nprim
      Prim_loc(s)   = PCELL(s,Im,Jm,Km)
      Prim_ghost(s) = Pg(s,1)
    end do

    call Visc_Variables ( Prim_loc,   Visc_loc )
    call Visc_Variables ( Prim_ghost, Visc_ghost )
    do s = 1, nprim; tmp(s) = PCELL(s,Im+1,Jm,Km); end do; call Visc_Variables ( tmp, Visc_ip )
    do s = 1, nprim; tmp(s) = PCELL(s,Im-1,Jm,Km); end do; call Visc_Variables ( tmp, Visc_im )
    do s = 1, nprim; tmp(s) = PCELL(s,Im,Jm+1,Km); end do; call Visc_Variables ( tmp, Visc_jp )
    do s = 1, nprim; tmp(s) = PCELL(s,Im,Jm-1,Km); end do; call Visc_Variables ( tmp, Visc_jm )
    do s = 1, nprim; tmp(s) = PCELL(s,Im,Jm,Km+1); end do; call Visc_Variables ( tmp, Visc_kp )
    do s = 1, nprim; tmp(s) = PCELL(s,Im,Jm,Km-1); end do; call Visc_Variables ( tmp, Visc_km )
    do s = 1, nprim; tmp(s) = Pg(s,3); end do; call Visc_Variables ( tmp, Visc_g3 )
    do s = 1, nprim; tmp(s) = Pg(s,4); end do; call Visc_Variables ( tmp, Visc_g4 )
    do s = 1, nprim; tmp(s) = Pg(s,5); end do; call Visc_Variables ( tmp, Visc_g5 )
    do s = 1, nprim; tmp(s) = Pg(s,6); end do; call Visc_Variables ( tmp, Visc_g6 )

    select case (Fm)
      case (1:2)
        do s = 1, nprim
          Gloc(s,1)   = ( Visc_ghost(s) - Visc_loc(s) ) * modfm2
          Gloc(s,2)   = ( Visc_jp(s) - Visc_jm(s) ) * 0.5d0
          Gloc(s,3)   = ( Visc_kp(s) - Visc_km(s) ) * 0.5d0
          Gghost(s,1) = ( Visc_loc(s) - Visc_ghost(s) ) * (-modfm2)
          Gghost(s,2) = ( Visc_g4(s) - Visc_g3(s) ) * 0.5d0
          Gghost(s,3) = ( Visc_g6(s) - Visc_g5(s) ) * 0.5d0
        end do
      case (3:4)
        do s = 1, nprim
          Gloc(s,1)   = ( Visc_ip(s) - Visc_im(s) ) * 0.5d0
          Gloc(s,2)   = ( Visc_ghost(s) - Visc_loc(s) ) * modfm2
          Gloc(s,3)   = ( Visc_kp(s) - Visc_km(s) ) * 0.5d0
          Gghost(s,1) = ( Visc_g4(s) - Visc_g3(s) ) * 0.5d0
          Gghost(s,2) = ( Visc_loc(s) - Visc_ghost(s) ) * (-modfm2)
          Gghost(s,3) = ( Visc_g6(s) - Visc_g5(s) ) * 0.5d0
        end do
      case (5:6)
        do s = 1, nprim
          Gloc(s,1)   = ( Visc_ip(s) - Visc_im(s) ) * 0.5d0
          Gloc(s,2)   = ( Visc_jp(s) - Visc_jm(s) ) * 0.5d0
          Gloc(s,3)   = ( Visc_ghost(s) - Visc_loc(s) ) * modfm2
          Gghost(s,1) = ( Visc_g4(s) - Visc_g3(s) ) * 0.5d0
          Gghost(s,2) = ( Visc_g6(s) - Visc_g5(s) ) * 0.5d0
          Gghost(s,3) = ( Visc_loc(s) - Visc_ghost(s) ) * (-modfm2)
        end do
    end select

    ! G = matmul(G, M)  ->  G2(s,kk) = sum_j G(s,j)*M(j,kk)
    do kk = 1, 3
      do s = 1, nprim
        G2loc(s,kk)   = Gloc(s,1)*M(1,kk)   + Gloc(s,2)*M(2,kk)   + Gloc(s,3)*M(3,kk)
        G2ghost(s,kk) = Gghost(s,1)*M(1,kk) + Gghost(s,2)*M(2,kk) + Gghost(s,3)*M(3,kk)
      end do
    end do

    do s = 1, nprim
      Prim_avg(s) = 0.5d0*( Prim_loc(s) + Prim_ghost(s) )
    end do

    call Compute_Diffusive_Flux ( Prim_avg, G2loc,   Area, Normal, Waldis, Flux_loc )
    call Compute_Diffusive_Flux ( Prim_avg, G2ghost, Area, Normal, Waldis, Flux_ghost )

    do s = 1, nprim
      WRCELL(s,Im,Jm,Km) = WRCELL(s,Im,Jm,Km) - 0.5d0 * modfm2 * ( Flux_loc(s) + Flux_ghost(s) )
    end do

  end subroutine BC_Connection_Visc


  !> Device twin of BC_Rotational_Periodic_Eul (3D wedge seam, faces 5-6).
  !> Self-contained: reads the block's OWN opposite-k stencil [Km-1,Km | 1,2],
  !> rotates the donor-side (v,w) by the wedge angle delthe, then a connection-
  !> style HLLC. The pair is driven only from face 6 (Fm/=6 -> no-op, mirrors the
  !> host `if (Fm==6)`): the flux is applied to r(Km) AND, momentum rotated back
  !> to the θ=0 frame, to r(Ks=1). Inviscid only (no host viscous periodic twin).
  subroutine BC_Rotational_Periodic_Eul ( PWARG, r, dl_f, n1, n2, n3, Im, Jm, Km, Fm, &
                                              normal, area, delthe )
    !$acc routine seq
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, f_ss_dev
    use MOSE_Lib_Riemann,              only: Riemann, rs_dev
    use MOSE_Lib_Limiters,             only: rlimiter_MINMOD, rlimiter, lim_dev

    implicit none
    integer,  intent(in)    :: n1, n2, n3, Im, Jm, Km, Fm
    real(R8), intent(in)    :: normal(3), area, delthe
    integer,  intent(in)    :: ncw, nrw, boff, pi, pj, wroff
    integer,  intent(in)    :: nur
    real(R8), intent(in)    :: ur(nur)   ! Ur of the cells (prec; 1 element otherwise): read by HLLC-PC only
    real(R8), intent(in)    :: w(ncw, nprim)
    real(R8), intent(inout) :: wr(nrw, nprim)   ! Res var-last: coalesced residual WRITE target
    real(R8), intent(in)    :: dl_f(3, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    real(R8), intent(inout) :: r(nprim, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    integer  :: Ks, K1, K2, K3, K4, s
    real(R8) :: Prim1(NPRIM_DEV), Prim2(NPRIM_DEV), Prim3(NPRIM_DEV), Prim4(NPRIM_DEV)
    real(R8) :: Prim_L(NPRIM_DEV), Prim_R(NPRIM_DEV)
    real(R8) :: Dl1, Dl2, Dl3, Dl4, Dist12, Dist23, Dist34, d21, d32, d43
    real(R8) :: rho_L, rho_R, Rtot_L, Rtot_R, a_L, a_R, F_r, F_u, F_v, F_w, F_E, urL_, urR_
    real(R8) :: Fmass, su, sel_L, sel_R, cd, sd, v_rot, w_rot, Fv2, Fw2, sp_s

    if (n3 == 1) then   ! Check_Mesh_Type: ndir /= 3 exactly when the block has one cell in k
      ! nk = 1 (2D planar or 2D axisymmetric): no Riemann flux across the k seam, only the terms of the face
      ! normal, the pressure on ny (axisymmetry) and the centrifugal terms -- MOSE_open BC_Rotational_Periodic_Eul
      ! (ndir /= 3) and the host twin above. Identically zero on a planar mesh (ny = 0, w = 0). The 3D stencil
      ! below on an nk = 1 block fluxed the cell against itself through the ghosts and fed w (2D Euler tests).
      rho_L = 0d0
      do s = 1, nsc
        rho_L = rho_L + PCELL(s,Im,Jm,Km)
      end do
      v_rot = PCELL(nv,Im,Jm,Km); w_rot = PCELL(nw,Im,Jm,Km)
      Fv2 =   v_rot * normal(3) + w_rot * normal(2)      ! Face_V
      Fw2 = - v_rot * normal(2) + w_rot * normal(3)      ! Face_W
      Fmass = rho_L * area * w_rot
      sd = real(1 - 2*mod(Fm,2), R8)                     ! FluxSign: -1 on face 5, +1 on face 6
      WRCELL(nv,Im,Jm,Km) = WRCELL(nv,Im,Jm,Km) + sd * ( Fv2 * Fmass + PCELL(np,Im,Jm,Km) * area * normal(2) )
      WRCELL(nw,Im,Jm,Km) = WRCELL(nw,Im,Jm,Km) + sd * ( Fw2 * Fmass )
      return
    end if
    if (Fm /= 6) return            ! pair driven from face 6 only (3D)

    Ks = 1
    K1 = Km - 1; K2 = Km; K3 = Ks; K4 = Ks + 1

    do s = 1, nprim
      Prim1(s) = PCELL(s,Im,Jm,K1); Prim2(s) = PCELL(s,Im,Jm,K2)
      Prim3(s) = PCELL(s,Im,Jm,K3); Prim4(s) = PCELL(s,Im,Jm,K4)
    end do

    ! Rotate the donor-side (K3,K4) velocity (v,w) by delthe (face-6 convention)
    cd = cos(delthe); sd = sin(delthe)
    w_rot =  Prim3(nw)*cd + Prim3(nv)*sd;  v_rot = -Prim3(nw)*sd + Prim3(nv)*cd
    Prim3(nv) = v_rot; Prim3(nw) = w_rot
    w_rot =  Prim4(nw)*cd + Prim4(nv)*sd;  v_rot = -Prim4(nw)*sd + Prim4(nv)*cd
    Prim4(nv) = v_rot; Prim4(nw) = w_rot

    Dl1 = dl_f(3,Im,Jm,K1)*0.5d0; Dl2 = dl_f(3,Im,Jm,K2)*0.5d0
    Dl3 = dl_f(3,Im,Jm,K3)*0.5d0; Dl4 = dl_f(3,Im,Jm,K4)*0.5d0
    Dist12 = Dl1+Dl2; Dist23 = Dl2+Dl3; Dist34 = Dl3+Dl4

    do s = 1, nprim
      d21 = (Prim2(s)-Prim1(s))/Dist12
      d32 = (Prim3(s)-Prim2(s))/Dist23
      d43 = (Prim4(s)-Prim3(s))/Dist34
      Prim_L(s) = Prim2(s) + rlimiter ( lim_dev, d32,d21)*Dl2
      Prim_R(s) = Prim3(s) - rlimiter ( lim_dev, d43,d32)*Dl3
    end do

    call co_rotot_Rtot_dev ( Prim_L(1:nsc), rho_L, Rtot_L )
    a_L = f_ss_dev ( Prim_L(1:nsc), Prim_L(np), rho_L, Rtot_L )
    call co_rotot_Rtot_dev ( Prim_R(1:nsc), rho_R, Rtot_R )
    a_R = f_ss_dev ( Prim_R(1:nsc), Prim_R(np), rho_R, Rtot_R )

    urL_ = 0d0; urR_ = 0d0
    if (rs_dev == 14) then   ! HLLC-PC: the two seam cells (MOSE_open blk%Ur(Im,Jm,K2), blk%Ur(Im,Jm,K3))
      urL_ = URCELL(Im,Jm,K2); urR_ = URCELL(Im,Jm,K3)
    end if
    call Riemann ( rs_dev,  Prim_L(1:nsc), Prim_L(nu), Prim_L(nv), Prim_L(nw), Prim_L(np), a_L, rho_L, &
                        Prim_R(1:nsc), Prim_R(nu), Prim_R(nv), Prim_R(nw), Prim_R(np), a_R, rho_R, &
                        1d0, urL_, urR_, normal(1), normal(2), normal(3), F_r, F_u, F_v, F_w, F_E )

    su = sign(0.5d0, F_r); sel_L = 0.5d0 + su; sel_R = su - 0.5d0
    Fmass = F_r * area

    ! Apply to r(Km) [face-6 frame]
    do s = 1, nsc
      sp_s = Fmass*( sel_L*Prim_L(s)/rho_L - sel_R*Prim_R(s)/rho_R )
      WRCELL(s,Im,Jm,Km) = WRCELL(s,Im,Jm,Km) + sp_s
    end do
    WRCELL(nu,Im,Jm,Km) = WRCELL(nu,Im,Jm,Km) + F_u*area
    WRCELL(nv,Im,Jm,Km) = WRCELL(nv,Im,Jm,Km) + F_v*area
    WRCELL(nw,Im,Jm,Km) = WRCELL(nw,Im,Jm,Km) + F_w*area
    WRCELL(np,Im,Jm,Km) = WRCELL(np,Im,Jm,Km) + F_E*area
    if (nprim>np) then
      do s = np+1, nprim
        WRCELL(s,Im,Jm,Km) = WRCELL(s,Im,Jm,Km) + Fmass*( sel_L*Prim_L(s)/rho_L - sel_R*Prim_R(s)/rho_R )
      end do
    end if

    ! Apply to r(Ks=1) [θ=0 frame: rotate momentum back, opposite sign]
    Fv2 = ( F_v*cd + F_w*sd )*area
    Fw2 = (-F_v*sd + F_w*cd )*area
    do s = 1, nsc
      sp_s = Fmass*( sel_L*Prim_L(s)/rho_L - sel_R*Prim_R(s)/rho_R )
      WRCELL(s,Im,Jm,Ks) = WRCELL(s,Im,Jm,Ks) - sp_s
    end do
    WRCELL(nu,Im,Jm,Ks) = WRCELL(nu,Im,Jm,Ks) - F_u*area
    WRCELL(nv,Im,Jm,Ks) = WRCELL(nv,Im,Jm,Ks) - Fv2
    WRCELL(nw,Im,Jm,Ks) = WRCELL(nw,Im,Jm,Ks) - Fw2
    WRCELL(np,Im,Jm,Ks) = WRCELL(np,Im,Jm,Ks) - F_E*area
    if (nprim>np) then
      do s = np+1, nprim
        WRCELL(s,Im,Jm,Ks) = WRCELL(s,Im,Jm,Ks) - Fmass*( sel_L*Prim_L(s)/rho_L - sel_R*Prim_R(s)/rho_R )
      end do
    end if

  end subroutine BC_Rotational_Periodic_Eul


  !> Device twin of Trapezoidal (isentropic cp/T integral). The host allocates a
  !> T_vec(n_points) where n_points spans the integer-Kelvin interval (can be
  !> hundreds-thousands) -> rewritten ARRAY-FREE: accumulate the trapezoidal
  !> segments on the fly over [T_lo, int_sx..int_dx, T_hi]. Bit-faithful order.
  subroutine Trapezoidal ( sup_T, BC_T0, BC_ci, Rgas, integral )
    !$acc routine seq
    use FLINT_Lib_Thermodynamic, only: cp_tab, f_tabT
    implicit none
    real(R8), intent(in)  :: sup_T, BC_T0, Rgas, BC_ci(nsc)
    real(R8), intent(out) :: integral
    integer  :: i, s, int_sx, int_dx
    real(R8) :: T_lo, T_hi, T_prev, T_cur, cp_prev, cp_cur

    if (sup_T > BC_T0) then
      T_lo = BC_T0; T_hi = sup_T
    else
      T_lo = sup_T; T_hi = BC_T0
    end if
    int_sx = ceiling(T_lo); int_dx = floor(T_hi)

    integral = 0d0
    T_prev = T_lo; cp_prev = 0d0
    do s = 1, nsc
      cp_prev = cp_prev + BC_ci(s) * f_tabT(T_prev, s, cp_tab)
    end do
    do i = int_sx, int_dx                 ! interior integer nodes
      T_cur = real(i, R8); cp_cur = 0d0
      do s = 1, nsc
        cp_cur = cp_cur + BC_ci(s) * f_tabT(T_cur, s, cp_tab)
      end do
      integral = integral + 0.5d0*(T_cur - T_prev)*( cp_prev/Rgas/T_prev + cp_cur/Rgas/T_cur )
      T_prev = T_cur; cp_prev = cp_cur
    end do
    T_cur = T_hi; cp_cur = 0d0             ! final segment to the high endpoint
    do s = 1, nsc
      cp_cur = cp_cur + BC_ci(s) * f_tabT(T_cur, s, cp_tab)
    end do
    integral = integral + 0.5d0*(T_cur - T_prev)*( cp_prev/Rgas/T_prev + cp_cur/Rgas/T_cur )

    if (sup_T > BC_T0) integral = -integral

  end subroutine Trapezoidal


  !> Device twin of BC_Inflow_Nozzle (type 22). Pressure-regime selector: outflow
  !> or p>=p0 -> error=1 (caller falls back to symmetry); else subsonic/supersonic
  !> nozzle. Full residual update r(:) += modfm2*Flux (mass+momentum+energy).
  subroutine BC_Inflow_Nozzle ( PWARG, r, n1, n2, n3, Im, Jm, Km, Fm, Normal, Area, t_Vec, &
                                    BC_Mach, BC_T0, BC_p0, BC_rel_fac, BC_alpha, BC_beta, BC_pAmb, &
                                    BC_ci, BC_RANS, error )
    !$acc routine seq
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, f_ss_dev, f_cp_dev
    implicit none
    integer,  intent(in)    :: n1, n2, n3, Im, Jm, Km, Fm
    real(R8), intent(in)    :: Normal(3), Area, t_Vec(3)
    real(R8), intent(in)    :: BC_Mach, BC_T0, BC_p0, BC_rel_fac, BC_alpha, BC_beta, BC_pAmb
    real(R8), intent(in)    :: BC_ci(nsc), BC_RANS(nrans)
    integer,  intent(in)    :: ncw, nrw, boff, pi, pj, wroff
    integer,  intent(in)    :: nur
    real(R8), intent(in)    :: ur(nur)   ! Ur of the cells (prec; 1 element otherwise): read by HLLC-PC only
    real(R8), intent(in)    :: w(ncw, nprim)
    real(R8), intent(inout) :: wr(nrw, nprim)   ! Res var-last: coalesced residual WRITE target
    real(R8), intent(inout) :: r(nprim, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
    integer,  intent(out)   :: error
    integer  :: modfm, modfm1, modfm2, gstep, Int_i, Int_j, Int_k, dir, s, Tint(2)
    real(R8) :: Bound_Prim(NPRIM_DEV), Int_Prim(NPRIM_DEV), Un, BC_Sign
    real(R8) :: Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem, InOut, Flux(NPRIM_DEV), Fmass
    real(R8) :: T, Tdiff, cp, pSub, pSup, g

    error  = 0
    dir    = (Fm + 1)/2
    modfm  = mod(Fm, 2)
    modfm1 = 1 - modfm
    modfm2 = 1 - 2*modfm
    gstep  = 2*modfm - 1
    Int_i = Im; Int_j = Jm; Int_k = Km
    select case (dir)
      case (1); Int_i = Im + gstep
      case (2); Int_j = Jm + gstep
      case (3); Int_k = Km + gstep
    end select
    do s = 1, nprim
      Bound_Prim(s) = PCELL(s, Im, Jm, Km)
      Int_Prim(s)   = PCELL(s, Int_i, Int_j, Int_k)
    end do

    pSub = BC_alpha * 1.d5
    pSup = BC_beta  * 1.d5
    g    = BC_pAmb
    BC_Sign = real(modfm2, R8)

    Un = Bound_Prim(nu)*Normal(1) + Bound_Prim(nv)*Normal(2) + Bound_Prim(nw)*Normal(3)
    call co_rotot_Rtot_dev ( Bound_Prim(1:nsc), Bound_rho, Bound_Rgas )
    Bound_Sound = f_ss_dev ( Bound_Prim(1:nsc), Bound_Prim(np), Bound_rho, Bound_Rgas )
    T       = Bound_Prim(np) / (Bound_rho*Bound_Rgas)
    Tint(1) = int(T); Tint(2) = Tint(1) + 1; Tdiff = T - Tint(1)
    cp          = f_cp_dev ( Bound_Prim(1:nsc), Tint, Tdiff, Bound_rho )
    Bound_Gamma = cp / (cp - Bound_Rgas)
    Riem  = Bound_Sound + BC_Sign*0.5d0*(Bound_Gamma-1d0)*Un
    InOut = BC_Sign*Un

    if ( InOut >= 0d0 .or. Bound_Prim(np) >= BC_p0 ) then
      error = 1
      return
    else if ( Bound_Prim(np) >= pSub ) then
      call BC_Nozzle_Subsonic ( Bound_Prim, Un, Bound_Sound, Bound_Gamma, Riem, 0d0, 0d0, &
                                    0d0, BC_ci, BC_RANS, BC_T0, BC_p0, BC_rel_fac, BC_Sign, t_Vec, &
                                    Normal, Area, Flux, Fmass )
    else if ( Bound_Prim(np) >= pSup ) then
      call BC_Nozzle_Subsonic ( Bound_Prim, Un, Bound_Sound, Bound_Gamma, Riem, 0d0, 0d0, &
                                    -10d0, BC_ci, BC_RANS, BC_T0, g, BC_rel_fac, BC_Sign, t_Vec, &
                                    Normal, Area, Flux, Fmass )
    else
      if (BC_Mach == 0d0) then
        call BC_Nozzle_Supersonic ( Bound_Prim, Int_Prim, modfm, modfm1, modfm2, g, BC_p0, &
                                        BC_T0, BC_ci, BC_RANS, 0d0, 0d0, Normal, Area, Flux, Fmass, 1d0 )
      else
        call BC_Nozzle_Supersonic ( Bound_Prim, Int_Prim, modfm, modfm1, modfm2, BC_Mach, BC_p0, &
                                        BC_T0, BC_ci, BC_RANS, 0d0, 0d0, Normal, Area, Flux, Fmass, 0d0 )
      end if
    end if

    do s = 1, nprim
      WRCELL(s,Im,Jm,Km) = WRCELL(s,Im,Jm,Km) + modfm2 * Flux(s)
    end do
  end subroutine BC_Inflow_Nozzle


  !> Device twin of BC_Nozzle_Supersonic. Newton solve for the sonic/mass-flux
  !> inflow state, isentropic p via Trapezoidal, then HLLC at the interface.
  subroutine BC_Nozzle_Supersonic ( Prim2, Int_Prim, modfm, modfm1, modfm2, g, BC_p0, BC_T0, &
                                        BC_ci, BC_RANS, alpha_app, beta_app, Normal, Area, Flux, Fmass, BC_supersonic )
    !$acc routine seq
    use FLINT_Lib_Thermodynamic,       only: f_tabT, h_tab, cp_tab, dcpi_tab, Ri_tab
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev, f_ss_dev
    use MOSE_Lib_Riemann,              only: Riemann, rs_dev
    use MOSE_Lib_Limiters,             only: rlimiter_MINMOD, rlimiter, lim_dev
    implicit none
    integer,  intent(in)  :: modfm, modfm1, modfm2
    real(R8), intent(in)  :: Prim2(NPRIM_DEV), Int_Prim(NPRIM_DEV), g, BC_p0, BC_T0, BC_supersonic
    real(R8), intent(in)  :: BC_ci(nsc), BC_RANS(nrans), alpha_app, beta_app, Normal(3), Area
    real(R8), intent(out) :: Flux(NPRIM_DEV), Fmass
    integer  :: s
    real(R8) :: Sup_Prim(NPRIM_DEV), Sup_T, Sup_Rgas, Sup_rho, Sup_Sound
    real(R8) :: Prim1(NPRIM_DEV), Prim3(NPRIM_DEV), Slope, Face_Prim(NPRIM_DEV)
    real(R8) :: Prim_L(NPRIM_DEV), Prim_R(NPRIM_DEV)
    real(R8) :: rho_L, rho_R, Rgas_L, Rgas_R, Sound_L, Sound_R, F_r, F_u, F_v, F_w, F_E
    real(R8) :: su, Sel_L, Sel_R
    real(R8) :: h0_, cpSup, dcpSup, hSup, uSup, Sup_h0_, duSup, Sup_p, dpSup, dh0_Sup, gamSup, dgamSup, daSup
    real(R8) :: fun_T, diff_h0_

    h0_ = 0d0; Sup_Rgas = 0d0
    do s = 1, nsc
      h0_      = h0_      + BC_ci(s) * f_tabT( BC_T0, s, h_tab )
      Sup_Rgas = Sup_Rgas + BC_ci(s) * Ri_tab(s)
    end do
    Sup_T = 0.5d0*BC_T0
    diff_h0_ = 1d0
    do while (diff_h0_ >= 1.d-9)
      hSup = 0d0; cpSup = 0d0; dcpSup = 0d0
      do s = 1, nsc
        hSup   = hSup   + BC_ci(s) * f_tabT( Sup_T, s, h_tab )
        cpSup  = cpSup  + BC_ci(s) * f_tabT( Sup_T, s, cp_tab )
        dcpSup = dcpSup + BC_ci(s) * f_tabT( Sup_T, s, dcpi_tab )
      end do
      if (BC_supersonic == 1.d0) then
        call Trapezoidal ( Sup_T, BC_T0, BC_ci, Sup_Rgas, fun_T )
        Sup_p   = BC_p0 / exp(fun_T)
        Sup_rho = Sup_p / (Sup_Rgas*Sup_T)
        uSup    = g / Sup_rho
        Sup_h0_ = hSup + 0.5d0*uSup**2
        diff_h0_ = abs(Sup_h0_-h0_)/abs(h0_)
        dpSup   = BC_p0/exp(fun_T) * cpSup/(Sup_Rgas*Sup_T)
        duSup   = g*Sup_Rgas*(Sup_p - dpSup*Sup_T)/Sup_p**2
        dh0_Sup = cpSup + uSup*duSup
      else
        gamSup    = cpSup/(cpSup-Sup_Rgas)
        Sup_Sound = sqrt(gamSup*Sup_Rgas*Sup_T)
        uSup      = Sup_Sound*g
        Sup_h0_   = hSup + 0.5d0*uSup**2
        diff_h0_  = abs(Sup_h0_-h0_)/abs(h0_)
        dgamSup   = -Sup_Rgas*dcpSup/(cpSup-Sup_Rgas)**2
        daSup     = sqrt(Sup_Rgas)*0.5d0/sqrt(gamSup*Sup_T)*(gamSup + Sup_T*dgamSup)
        duSup     = g*daSup
        dh0_Sup   = cpSup + uSup*duSup
      end if
      Sup_T = Sup_T - (Sup_h0_-h0_)/dh0_Sup
    end do
    call Trapezoidal ( Sup_T, BC_T0, BC_ci, Sup_Rgas, fun_T )
    Sup_Prim(np) = BC_p0/exp(fun_T)
    Sup_rho = Sup_Prim(np)/(Sup_Rgas*Sup_T)
    do s = 1, nsc; Sup_Prim(s) = Sup_rho*BC_ci(s); end do
    Sup_Sound = f_ss_dev ( Sup_Prim(1:nsc), Sup_Prim(np), Sup_rho, Sup_Rgas )
    if (BC_supersonic == 1.d0) then
      Sup_Prim(nu) = g/Sup_rho*Normal(1); Sup_Prim(nv) = g/Sup_rho*Normal(2); Sup_Prim(nw) = g/Sup_rho*Normal(3)
    else
      Sup_Prim(nu) = g*Sup_Sound*Normal(1); Sup_Prim(nv) = g*Sup_Sound*Normal(2); Sup_Prim(nw) = g*Sup_Sound*Normal(3)
    end if
    if (RANS) then
      do s = nt, nprim; Sup_Prim(s) = Sup_rho*BC_RANS(s-nt+1); end do
    end if
    do s = 1, nprim
      Prim1(s) = modfm*Sup_Prim(s) + modfm1*Int_Prim(s)
      Prim3(s) = modfm*Int_Prim(s) + modfm1*Sup_Prim(s)
    end do
    do s = 1, nprim
      Slope = rlimiter ( lim_dev,  Prim3(s)-Prim2(s), Prim2(s)-Prim1(s) )
      Face_Prim(s) = Prim2(s) + 0.5d0*Slope*modfm2
    end do
    do s = 1, nprim
      Prim_L(s) = modfm*Sup_Prim(s) + modfm1*Face_Prim(s)
      Prim_R(s) = modfm*Face_Prim(s) + modfm1*Sup_Prim(s)
    end do
    call co_rotot_Rtot_dev ( Prim_L(1:nsc), rho_L, Rgas_L )
    Sound_L = f_ss_dev ( Prim_L(1:nsc), Prim_L(np), rho_L, Rgas_L )
    call co_rotot_Rtot_dev ( Prim_R(1:nsc), rho_R, Rgas_R )
    Sound_R = f_ss_dev ( Prim_R(1:nsc), Prim_R(np), rho_R, Rgas_R )
    call Riemann ( rs_dev,  Prim_L(1:nsc), Prim_L(nu), Prim_L(nv), Prim_L(nw), Prim_L(np), Sound_L, rho_L, &
                        Prim_R(1:nsc), Prim_R(nu), Prim_R(nv), Prim_R(nw), Prim_R(np), Sound_R, rho_R, &
                        1d0, Sound_L, Sound_R, Normal(1), Normal(2), Normal(3), F_r, F_u, F_v, F_w, F_E )   ! Ur = a: no preconditioning (MOSE_open)
    su = sign(1d0, F_r); Sel_L = (1d0+su)/2d0; Sel_R = (su-1d0)/2d0
    do s = 1, nprim; Flux(s) = 0d0; end do
    Fmass = F_r*Area
    do s = 1, nsc; Flux(s) = Fmass*( Sel_L*Prim_L(s)/rho_L - Sel_R*Prim_R(s)/rho_R ); end do
    Flux(nu) = F_u*Area; Flux(nv) = F_v*Area; Flux(nw) = F_w*Area; Flux(np) = F_E*Area
    if (RANS) then
      do s = nt, nprim; Flux(s) = Fmass*( Sel_L*Prim_L(s)/rho_L - Sel_R*Prim_R(s)/rho_R ); end do
    end if
  end subroutine BC_Nozzle_Supersonic


  !> Device twin of BC_Nozzle_Subsonic. Triple-nested Newton across the three
  !> input regimes (BC_Mach>=0: assigned T0; ==-10: assigned h0; <0: assigned mdot).
  subroutine BC_Nozzle_Subsonic ( Prim, Un, Sound, Gamma, Riem, alpha_app, beta_app, BC_Mach, &
                                      BC_ci, BC_RANS, BC_T0, BC_p0, BC_rel_fac, BC_Sign, t_Vec, Normal, Area, Flux, Fmass )
    !$acc routine seq
    use FLINT_Lib_Thermodynamic,       only: f_tabT, h_tab, cp_tab, dcpi_tab, Ri_tab
    use FLINT_Lib_ThermoTransport_dev, only: co_rotot_Rtot_dev
    implicit none
    real(R8), intent(in)  :: Prim(NPRIM_DEV), Un, Sound, Gamma, Riem, alpha_app, beta_app, BC_Mach
    real(R8), intent(in)  :: BC_ci(nsc), BC_RANS(nrans), BC_T0, BC_p0, BC_rel_fac, BC_Sign, t_Vec(3), Normal(3), Area
    real(R8), intent(out) :: Flux(NPRIM_DEV), Fmass
    integer  :: s
    real(R8) :: b_Vec(3), b_Mod, alpha, beta, XA, XB, XC, XD, XE, XF, check, E1, E2
    real(R8) :: p3, Ut3, Ub3, mdot, Rgas3, cp3, T3, rho3, h3, h0_3, u3, h0_3_iter, p03, h3_iter
    real(R8) :: Face_T, Face_Prim(NPRIM_DEV), Face_rho, Face_Rgas, Face_Enthalpy, Face_Un, Un3
    real(R8) :: diff_h0_, da, du, dgam, ddel, dT, dp, dh0_, dp0
    real(R8) :: Rgas_acu, T_acu, h_acu, p_acu, u_acu, cp_acu, dcp_acu, gamma_acu, a_acu, T_Bound, Bound_ci(nsc), fun_T, diff_T
    real(R8) :: rho_Bound

    b_Vec(1) = Normal(2)*t_Vec(3) - Normal(3)*t_Vec(2)
    b_Vec(2) = Normal(3)*t_Vec(1) - Normal(1)*t_Vec(3)
    b_Vec(3) = Normal(1)*t_Vec(2) - Normal(2)*t_Vec(1)
    b_Mod = sqrt(b_Vec(1)**2 + b_Vec(2)**2 + b_Vec(3)**2)
    b_Vec(1) = b_Vec(1)/b_Mod; b_Vec(2) = b_Vec(2)/b_Mod; b_Vec(3) = b_Vec(3)/b_Mod
    if ( alpha_app == 0d0 .and. beta_app == 0d0 ) then
      alpha = atan( Normal(2)/(Normal(1)+1d-20) ); beta = atan( Normal(3)/(Normal(1)+1d-20) )
    else
      alpha = alpha_app; beta = beta_app
    end if
    XA = Normal(2) - Normal(1)*tan(alpha); XB = t_Vec(2) - t_Vec(1)*tan(alpha); XC = b_Vec(2) - b_Vec(1)*tan(alpha)
    XD = Normal(3) - Normal(1)*tan(beta);  XE = t_Vec(3) - t_Vec(1)*tan(beta);  XF = b_Vec(3) - b_Vec(1)*tan(beta)
    check = XC*XE - XF*XB
    E1 = (XF*XA - XC*XD)/check
    E2 = (XD*XB - XE*XA)/check

    mdot = BC_p0
    Rgas3 = 0d0
    do s = 1, nsc; Rgas3 = Rgas3 + BC_ci(s)*Ri_tab(s); end do
    call co_rotot_Rtot_dev ( Prim(1:nsc), rho_Bound, Rgas_acu )
    diff_h0_ = 1.d0
    h0_3 = 0d0
    if ( BC_Mach == -10d0 .or. BC_Mach >= 0.d0 ) then
      do s = 1, nsc; h0_3 = h0_3 + BC_ci(s)*f_tabT( BC_T0, s, h_tab ); end do
    end if
    T_Bound = Prim(np)/rho_Bound/Rgas_acu
    T_acu   = 1.1d0*T_Bound
    do s = 1, nsc; Bound_ci(s) = Prim(s)/rho_Bound; end do

    do while (diff_h0_ >= 1.d-6)
      h_acu = 0d0; cp_acu = 0d0; dcp_acu = 0d0
      do s = 1, nsc
        h_acu   = h_acu   + Bound_ci(s)*f_tabT( T_acu, s, h_tab )
        cp_acu  = cp_acu  + Bound_ci(s)*f_tabT( T_acu, s, cp_tab )
        dcp_acu = dcp_acu + Bound_ci(s)*f_tabT( T_acu, s, dcpi_tab )
      end do
      gamma_acu = cp_acu/(cp_acu-Rgas_acu)
      a_acu     = sqrt(gamma_acu*Rgas_acu*T_acu)
      call Trapezoidal ( T_Bound, T_acu, Bound_ci, Rgas_acu, fun_T )
      p_acu = Prim(np)*exp(fun_T)
      u_acu = (Riem - a_acu)/(BC_Sign*(gamma_acu-1d0)*0.5d0)
      p3 = p_acu; u3 = u_acu
      dp   = Prim(np)*exp(fun_T)*cp_acu/(Rgas_acu*T_acu)
      dgam = -Rgas_acu*dcp_acu/(cp_acu-Rgas_acu)**2
      ddel = dgam*0.5d0
      da   = sqrt(Rgas_acu)*0.5d0/sqrt(gamma_acu*T_acu)*(gamma_acu + T_acu*dgam)
      du   = (-da*BC_Sign*(gamma_acu-1d0)*0.5d0 - (Riem-a_acu)*(BC_Sign*ddel))/(BC_Sign*(gamma_acu-1d0)*0.5d0)**2
      if (BC_Mach >= 0d0) then
        diff_T = 1.d0; h3 = h0_3 - 0.5d0*u3**2; T3 = 0.9d0*BC_T0
        do while (diff_T >= 1.d-6)
          cp3 = 0d0; h3_iter = 0d0
          do s = 1, nsc
            h3_iter = h3_iter + BC_ci(s)*f_tabT( T3, s, h_tab )
            cp3     = cp3     + BC_ci(s)*f_tabT( T3, s, cp_tab )
          end do
          diff_T = abs(h3-h3_iter)/abs(h3)
          T3 = T3 - (h3_iter-h3)/cp3
        end do
        cp3 = 0d0
        do s = 1, nsc; cp3 = cp3 + BC_ci(s)*f_tabT( T3, s, cp_tab ); end do
        call Trapezoidal ( T3, BC_T0, BC_ci, Rgas3, fun_T )
        p03 = p3*exp(fun_T)
        diff_h0_ = abs(p03 - BC_p0)/BC_p0
        dp0 = dp*exp(fun_T) + p3*exp(fun_T)/(Rgas3*T3)*u3*du
        T_acu = T_acu - (p03 - BC_p0)/dp0
      else
        dT = (1d0/Rgas3/mdot)*(dp*u3 + du*p3)
        rho3 = mdot/u3; T3 = p3/(Rgas3*rho3)
        h3 = 0d0; cp3 = 0d0
        do s = 1, nsc
          h3  = h3  + BC_ci(s)*f_tabT( T3, s, h_tab )
          cp3 = cp3 + BC_ci(s)*f_tabT( T3, s, cp_tab )
        end do
        if ( BC_Mach == -10d0 ) then
          h0_3_iter = h3 + 0.5d0*u3**2
          diff_h0_  = abs(h0_3 - h0_3_iter)/abs(h0_3)
          dh0_      = cp3*dT + u3*du
          T_acu     = T_acu - (h0_3_iter - h0_3)/dh0_
        else
          diff_h0_ = abs(T3 - BC_T0)/BC_T0
          T_acu    = T_acu - (T3 - BC_T0)/dT
        end if
      end if
    end do

    ! final state at converged T_acu
    h_acu = 0d0; cp_acu = 0d0; dcp_acu = 0d0
    do s = 1, nsc
      h_acu  = h_acu  + Bound_ci(s)*f_tabT( T_acu, s, h_tab )
      cp_acu = cp_acu + Bound_ci(s)*f_tabT( T_acu, s, cp_tab )
    end do
    gamma_acu = cp_acu/(cp_acu-Rgas_acu)
    a_acu     = sqrt(gamma_acu*Rgas_acu*T_acu)
    call Trapezoidal ( T_Bound, T_acu, Bound_ci, Rgas_acu, fun_T )
    p_acu = Prim(np)*exp(fun_T)
    u_acu = (Riem - a_acu)/(BC_Sign*(gamma_acu-1d0)*0.5d0)
    p3 = p_acu; u3 = u_acu
    if (BC_Mach >= 0.d0) then
      h3 = h0_3 - 0.5d0*u3**2; T3 = 0.9d0*BC_T0; diff_T = 1.d0
      do while (diff_T >= 1.d-6)
        cp3 = 0d0; h3_iter = 0d0
        do s = 1, nsc
          h3_iter = h3_iter + BC_ci(s)*f_tabT( T3, s, h_tab )
          cp3     = cp3     + BC_ci(s)*f_tabT( T3, s, cp_tab )
        end do
        diff_T = abs(h3-h3_iter)/abs(h3)
        T3 = T3 - (h3_iter-h3)/cp3
      end do
      rho3 = p3/(Rgas3*T3)
    else
      rho3 = mdot/u3
    end if

    Un3 = u3; Ut3 = Un3*E1; Ub3 = Un3*E2
    Face_Rgas = Rgas3; Face_rho = rho3
    do s = 1, nsc; Face_Prim(s) = Face_rho*BC_ci(s); end do
    Face_Prim(nu) = Un3*Normal(1) + Ut3*t_Vec(1) + Ub3*b_Vec(1)
    Face_Prim(nv) = Un3*Normal(2) + Ut3*t_Vec(2) + Ub3*b_Vec(2)
    Face_Prim(nw) = Un3*Normal(3) + Ut3*t_Vec(3) + Ub3*b_Vec(3)
    Face_Prim(np) = p3
    Face_Un = Un3
    if (RANS) then
      do s = nt, nprim; Face_Prim(s) = BC_RANS(s-nt+1)*Face_rho; end do
    end if

    Face_Enthalpy = 0d0
    Face_T = Face_Prim(np)/(Face_Rgas*Face_rho)
    do s = 1, nprim; Flux(s) = 0d0; end do
    do s = 1, nsc
      Face_Enthalpy = Face_Enthalpy + Face_Prim(s)/Face_rho*f_tabT( Face_T, s, h_tab )
      Flux(s) = Face_Prim(s)*Face_Un*Area
    end do
    Fmass = Face_rho*Face_Un*Area
    Flux(nu) = Fmass*Face_Prim(nu) + Face_Prim(np)*Area*Normal(1)
    Flux(nv) = Fmass*Face_Prim(nv) + Face_Prim(np)*Area*Normal(2)
    Flux(nw) = Fmass*Face_Prim(nw) + Face_Prim(np)*Area*Normal(3)
    Flux(np) = Fmass*( 0.5d0*(Face_Prim(nu)**2 + Face_Prim(nv)**2 + Face_Prim(nw)**2) + Face_Enthalpy )
    if (RANS) then
      do s = nt, nprim; Flux(s) = Face_Prim(s)/Face_rho*Fmass; end do
    end if
  end subroutine BC_Nozzle_Subsonic


end module MOSE_Lib_BC
