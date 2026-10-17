
module FLINT_Lib_Thermodynamic
  implicit none
  public:: E0, H0
  public:: prim2cons, cons2prim
  !-----------------------------------------------------------------------------------------------------------------------------------

  character(len=128)      :: FLINT_phase_prefix=''

  integer                 :: ns                 ! Number of species
  integer                 :: ne                 ! Number of atomic elements
  integer                 :: Tmin,Tmax          ! Extreme temperatures in tables
  real(kind=8), parameter :: Runiv = 8314.46261815324d0  ! Universal gas constant [J/(kmol*K)]: exact SI value N_A*k_B (CODATA 2018)

  integer, parameter :: s_str_len = 20
  character(len=s_str_len), allocatable :: species_names(:)
  character(len=s_str_len), allocatable :: elements_names(:)

  real(kind=8), dimension(:,:), allocatable :: species_composition

  real(kind=8), dimension(:), allocatable   :: wm_tab, Ri_tab
  real(kind=8), dimension(:,:), allocatable :: cp_tab, dcpi_tab, h_tab, s_tab
  ! Device copies for the OpenACC chemistry port (uploaded once by
  ! flint_acc_upload_tables; referenced from acc-routine device code).
  !$acc declare create(ns, Tmin, Tmax, wm_tab, Ri_tab, cp_tab, dcpi_tab, h_tab)
  real(kind=8), dimension(:,:), allocatable :: mi_tab, k_tab

  ! Species-pair constants used by co_fiij (Wilke mixing rule). They depend
  ! only on molecular weights, so are precomputed once after wm_tab is loaded.
  real(kind=8), dimension(:,:), allocatable :: Mi_Mj_pow_m025 !> (Wm(i)/Wm(j))**(-0.25)
  real(kind=8), dimension(:,:), allocatable :: inv_sqrt8_1p   !> 1 / sqrt(8 * (1 + Wm(i)/Wm(j)))
  ! Transport tables for the device thermo-cache kernel (GPU Phase 2). Declared
  ! after all four arrays exist (declare must follow the declarations).
  !$acc declare create(mi_tab, k_tab, Mi_Mj_pow_m025, inv_sqrt8_1p)

  ! Binary diffusion coefficient table (T, unique species pair).
  ! Stored upper-triangular: ndij = ns*(ns-1)/2 columns, ordered (1,2),(1,3),..,(1,ns),(2,3),..
  ! Binary diffusivities scale as D ~ 1/p, so the table is built at a reference pressure
  ! dij_pref [Pa] (parsed from the diffusion.dat TITLE) and rescaled by dij_pref/p at use.
  real(kind=8), dimension(:,:), allocatable :: dij_tab
  integer :: ndij = 0
  real(kind=8) :: dij_pref = 101325.d0
  ! Device copy of the binary diffusion table for co_DS_expr_dev (uploaded by
  ! flint_acc_upload_thermo when read_idealgas_diffusion has allocated it).
  !$acc declare create(dij_tab, dij_pref)

contains

  !> Function for computing the total specific entalpy (per unit of mass)
  pure function H0(p,r,u) result(entalpy)
    implicit none
    real(kind=8), intent(IN):: p       !< Pressure.
    real(kind=8), intent(IN):: r(:)    !< Density (\f$\rho\f$).
    real(kind=8), intent(IN):: u       !< Module of velocity vector.
    real(kind=8)            :: entalpy !< Total specific entalpy (per unit of mass).
    entalpy = H(p,r)+0.5d0*u*u
  endfunction H0

  !> Function for computing the static specific entalpy (per unit of mass)
  pure function H(p,r) result(entalpy)
    implicit none
    real(kind=8), intent(IN):: p       !< Pressure.
    real(kind=8), intent(IN):: r(:)    !< Density (\f$\rho\f$).
    real(kind=8)            :: entalpy !< Static specific entalpy (per unit of mass).
    integer  :: s
    real(kind=8) :: rho, h_, Rgas
    entalpy = 0.d0
    rho = sum(r)
    Rgas = f_Rtot(r)
    do s = 1, ns
      h_ = f_tab(p,rho,Rgas,s,h_tab)
      entalpy = entalpy+h_*r(s)/rho
    enddo
  endfunction H


  !> Function for computing the total specific energy (per unit of mass)
  pure function E0(p,r,u) result(energy)
    real(kind=8), intent(IN):: p      !< Pressure.
    real(kind=8), intent(IN):: r(:)   !< Density (\f$\rho\f$).
    real(kind=8), intent(IN):: u      !< Module of velocity vector.
    real(kind=8)               energy !< Total specific energy (per unit of mass).
    energy = E(p,r)+0.5d0*u*u
  endfunction E0


  !> Function for computing the static specific energy (per unit of mass)
  pure function E(p,r) result(energy)
    implicit none
    real(kind=8), intent(IN):: p      !< Pressure.
    real(kind=8), intent(IN):: r(:)   !< Density (\f$\rho\f$).
    real(kind=8)               energy !< Static specific energy (per unit of mass).
    energy = H(p,r)-p/sum(r)
  endfunction E


  pure function S0i(i,T) result(s)
    implicit none
    integer, intent(IN)     :: i
    real(kind=8), intent(IN):: T
    real(kind=8)            :: s

    s = f_tabT(T,i,s_tab)

  end function S0i


  pure function H0i(i,T) result(h)
    implicit none
    integer, intent(IN)     :: i
    real(kind=8), intent(IN):: T
    real(kind=8)            :: h

    h = f_tabT(T,i,h_tab)

  end function H0i


  !> Computes the values of C_P_0, H_0 and S_0
  pure subroutine COMP_THERMO_QUANTS(temp,N_reac,C_P_0,H_0,S_0)
    use FLINT_CEA_data, only: Rr
    !! I/O
    real(8), intent(in)          :: temp
    integer, intent(in)           :: N_reac
    real(8), intent(out)         :: C_P_0(N_reac), H_0(N_reac), S_0(N_reac)
    !! internal
    integer                       :: i_reac

    do i_reac = 1, N_reac
      C_P_0(i_reac) = cpi(i_reac,temp) * wm_tab(i_reac) / (Rr)
      H_0(i_reac) = H0i(i_reac,temp) * wm_tab(i_reac) / (Rr*temp)
      S_0(i_reac) = S0i(i_reac,temp) * wm_tab(i_reac) / (Rr)
    end do

  end subroutine COMP_THERMO_QUANTS


  !> Function for the perfect gas equation of state
  pure function EOS(p,ro,R) result(T)
    implicit none
    real(kind=8), intent(IN):: p      !< Pressure.
    real(kind=8), intent(IN):: ro     !< Density (\f$\rho\f$).
    real(kind=8), intent(IN):: R      !< Gas constant.
    real(kind=8)            :: T      !< Temperature.
    T = p/(ro*R)
  endfunction EOS


  !> Function to compute the conservative variables starting from the primitive ones.
  pure function cons2prim(cons,T) result(prim)
    implicit none
    real(kind=8), intent(IN) :: cons(:)
    real(kind=8), intent(IN) :: T
    real(kind=8)             :: prim(size(cons))
    real(kind=8) :: Rgas, rho, rhoi(1:ns)
    real(kind=8) :: T_iter ,p_iter, cv_iter, energy_iter, diff_en
    real(kind=8) :: energy, vel
    integer :: iter_T, np, nui, nuf

    np = size(cons)
    nui = ns + 1
    nuf = np - 1

    Rgas = f_Rtot(cons(1:ns))
    rhoi = cons(1:ns)
    rho = sum(rhoi)
    prim(1:ns) = rhoi                       ! density
    prim(nui:nuf) = cons(nui:nuf)/rho        ! velocity

    ! Newton Raphson
    vel = norm2(prim(nui:nuf))
    energy = cons(np)/rho-0.5d0*vel*vel
    diff_en = 1.0d+3; iter_T = 1; T_iter = T
    do while ( (abs(diff_en/energy) >= 1d-6) .and. (iter_T <= 20) )
      p_iter = rho*Rgas*T_iter
      cv_iter = f_cp(rhoi,T_iter,rho)-Rgas
      energy_iter = E(p_iter,rhoi)
      diff_en = abs(energy_iter-energy)
      T_iter = T_iter+(energy-energy_iter)/(cv_iter)
      iter_T = iter_T+1
    enddo
    !if (iter_T==21) write(*,*)'Warning: max NR iter reached'
    prim(np) = rho*Rgas*T_iter               ! pressure

  end function cons2prim


  !> Subroutine to compute the primitive variables starting from the conserveative ones 
  pure function prim2cons(prim) result(cons)
    implicit none
    real(kind=8), intent(IN) :: prim(:)
    real(kind=8)             :: cons(size(prim))
    real(kind=8) :: Rgas, rho, vel
    integer      :: nui, nuf

    nui = ns + 1
    nuf = size(prim) - 1

    Rgas = f_Rtot(prim(1:ns))
    rho = sum(prim(1:ns))
    vel = norm2(prim(nui:nuf))
    cons(1:ns)   = prim(1:ns)                           ! density
    cons(nui:nuf) = rho*prim(nui:nuf)                     ! momentum
    cons(nuf+1)   = rho*E0(prim(nuf+1),prim(1:ns),vel)   ! energy

  end function prim2cons


  pure subroutine co_rotot_Rtot ( rhoi, rho, R )
    implicit none
    real(8), intent(in)  :: rhoi(ns)
    real(8), intent(out) :: rho, R
  
    rho = sum(rhoi)
    R =f_Rtot(rhoi)
  
  endsubroutine co_rotot_Rtot


  pure function f_molecularWeight(rhoi) result(result)
    !$acc routine seq
    implicit none
    real(8), intent(in)  :: rhoi(ns)
    real(8) :: result
    integer :: s
    real(8) :: rho

    rho = sum(rhoi)
    result = 0.d0
    do s = 1, ns
      result = result+(rhoi(s)/rho/Wm_tab(s))
    enddo
    result = 1.d0/result

  endfunction f_molecularWeight


  pure function f_Rtot(rhoi) result(result)
    !$acc routine seq
    implicit none
    real(8), intent(in)  :: rhoi(ns)
    real(8) :: result
    integer :: s
    real(8) :: rho

    rho = sum(rhoi)
    result = 0.d0
    do s = 1, ns
      result = result+Ri_tab(s)*rhoi(s)/rho
    enddo

  endfunction f_Rtot


  pure function cpi(i,T) result(result)
    implicit none
    integer, intent(in) :: i
    real(8), intent(in) :: T
    real(8) :: Tdiff
    integer :: T_i, Tint(2)
    real(8) :: result

    T_i = idint(T)
    Tdiff  = T-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    result = f_tabT_expr(i,cp_tab,Tint,Tdiff)

  endfunction cpi


  pure function f_cp(rhoi,T,rho) result(result)
    implicit none
    real(8), intent(in) :: rhoi(ns)
    real(8), intent(in) :: T, rho
    real(8) :: cpi(ns)
    real(8) :: Tdiff
    integer :: s, T_i, Tint(2)
    real(8) :: result

    T_i = idint(T)
    Tdiff  = T-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    result = 0.d0
    do s = 1, ns
      cpi(s) = f_tabT_expr(s,cp_tab,Tint,Tdiff)
      result = result+rhoi(s)/rho*cpi(s)
    enddo

  endfunction f_cp


  pure function f_cp_expr(rhoi,Tint,Tdiff,rho) result(result)
    implicit none
    real(8), intent(in) :: rhoi(ns)
    real(8), intent(in) :: Tdiff, rho
    integer, intent(in) :: Tint(2)
    real(8) :: cpi(ns)
    integer :: s
    real(8) :: result

    result = 0.d0
    do s = 1, ns
      cpi(s) = f_tabT_expr(s,cp_tab,Tint,Tdiff)
      result = result+rhoi(s)/rho*cpi(s)
    enddo

  endfunction f_cp_expr


  pure function f_ss(rhoi,p,rho,Rtot) result(result)
    implicit none
    real(8), intent(in) :: rhoi(ns)
    real(8), intent(in) :: p, Rtot, rho
    real(8) :: T, cp, gam
    real(8) :: result
    
    T = p/(Rtot*rho)
    cp = f_cp(rhoi,T,rho)
    gam = cp/(cp-Rtot)
    result = dsqrt(gam*Rtot*T)

  endfunction f_ss


  pure function f_gamma(rhoi,p,rho,Rtot) result(gam)
    implicit none
    real(8), intent(in)  :: rhoi(ns), p, rho, Rtot
    real(8) :: gam
    real(8) :: T, cp

    T = p/(Rtot*rho)
    cp = f_cp(rhoi,T,rho)
    gam = cp/(cp-Rtot)

  endfunction f_gamma


  pure function mass2molar(Y_species) result(x_species)
    implicit none
    real(8), intent(in) :: Y_species(ns)
    real(8)             :: x_species(ns)
    integer  :: j
    real(8)  :: n(ns)

    do j = 1, ns
       n(j) = Y_species(j) / wm_tab(j)
    end do
    x_species(:) = n(:) / sum(n)

  end function mass2molar


  pure function molar2mass(x_species) result(y_species)
    implicit none
    real(8), intent(in) :: x_species(ns)
    real(8)             :: y_species(ns)
    integer  :: j
    real(8)  :: n(ns)

    do j = 1, ns
       n(j) = x_species(j) * wm_tab(j)
    end do
    y_species(:) = n(:) / sum(n)

  end function molar2mass


  pure subroutine co_fiij(fi,mi)
    implicit none
    real(8), intent(in)  :: mi(ns)
    real(8), intent(out) :: fi(ns,ns)
    real(8) :: sqrt_mi_num(ns), inv_sqrt_mi_den(ns), ratio, one_plus
    integer :: i, j

    do i = 1, ns
      sqrt_mi_num(i)     = dsqrt(mi(i))
      inv_sqrt_mi_den(i) = 1.d0 / dsqrt(mi(i) + 1d-20)
    enddo

    do j = 1, ns
      do i = 1, ns
        ratio    = sqrt_mi_num(i) * inv_sqrt_mi_den(j) * Mi_Mj_pow_m025(i,j)
        one_plus = 1.d0 + ratio
        fi(i,j)  = one_plus * one_plus * inv_sqrt8_1p(i,j)
      enddo
    enddo

  endsubroutine co_fiij


  pure function f_laminarViscosity(rhoi,p,rho,Rtot) result(mil)
    implicit none
    real(8), intent(in)     :: rhoi(ns), p, rho, Rtot
    real(8)                 :: mil
    real(8), dimension(ns)  :: mil_den, Xi, mi_fiij
    real(8)                 :: fi(ns,ns)
    real(8)                 :: Wmtot, T, Tdiff
    integer                 :: s ,i ,j ,T_i, Tint(2)

    T = p/(Rtot*rho)
    T_i = idint(T)
    Tdiff = T-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    Wmtot = f_molecularWeight(rhoi)

    !calcolo delle frazioni molari Xi, viscosità laminare da tabella, mi(s)
    do s = 1, ns
      Xi(s) = (rhoi(s)*Wmtot)/(rho*Wm_tab(s))
      mi_fiij(s) = f_tabT_expr(s,mi_tab,Tint,Tdiff)
    enddo

    ! calcolo del denominatore della legge di Wilke
    call co_fiij(fi,mi_fiij)

    do i = 1, ns
      mil_den(i)=0.d0
      do j = 1, ns
        mil_den(i) = mil_den(i)+Xi(j)*fi(i,j)
      enddo
    enddo
    
    ! calcolo della vicosità laminare
    mil = 0.d0
    do s = 1, ns
      mil = mil+Xi(s)*mi_fiij(s)/mil_den(s)
    enddo

  endfunction f_laminarViscosity


  pure subroutine co_k_mi_lam_Wilke(rhoi,rho,T,milam,klam)
    implicit none
    real(8), intent(in)  :: rhoi(ns), rho, T
    real(8), intent(out) :: milam, klam
    real(8) :: Tdiff
    integer :: T_i, Tint(2)

    T_i = idint(T)
    Tdiff  = T-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    call co_k_mi_lam_Wilke_expr(rhoi,rho,Tint,Tdiff,milam,klam)

  endsubroutine co_k_mi_lam_Wilke


  pure subroutine co_k_mi_lam_Wilke_expr(rhoi,rho,Tint,Tdiff,milam,klam)
    implicit none
    integer, intent(in)  :: Tint(2)
    real(8), intent(in)  :: rhoi(ns), rho, Tdiff
    real(8), intent(out) :: milam, klam
    real(8), dimension(ns) :: lam_den, Xi
    real(8)                 :: fi(ns,ns), klam_i(ns), milam_i(ns)
    real(8)                 :: Wmtot, inv_lam_den
    integer                 :: s ,i ,j

    Wmtot = f_molecularWeight(rhoi)

    ! calcolo delle frazioni molari Xi, viscosità laminare da tabella, mi(s)
    do s = 1, ns
      Xi(s) = rhoi(s)*Wmtot/(rho*Wm_tab(s))
      milam_i(s) = f_tabT_expr(s,mi_tab,Tint,Tdiff)
      klam_i(s)  = f_tabT_expr(s,k_tab,Tint,Tdiff)
    enddo

    ! calcolo del denominatore della legge di Wilke (same for k and mi computations)
    call co_fiij(fi,milam_i)

    do i=1,ns
      lam_den(i)=0.d0
      do j=1,ns
        lam_den(i)=lam_den(i)+Xi(j)*fi(i,j)
      enddo
    enddo

    ! calcolo della vicosità laminare
    milam = 0.d0
    klam = 0.d0
    do s = 1, ns
      inv_lam_den = 1/lam_den(s)
      milam = milam+Xi(s)*milam_i(s)*inv_lam_den
      klam = klam+Xi(s)*klam_i(s)*inv_lam_den
    enddo

  endsubroutine co_k_mi_lam_Wilke_expr


  !> Mixture-averaged (Wilke/Wassiljewa) multicomponent diffusion coefficients.
  !> Thin wrapper: derives temperature from the EOS and forwards to co_DS_expr.
  pure subroutine co_DS(rhoi,p,rho,Rtot,Dm)
    implicit none
    real(8), intent(in)  :: rhoi(ns), p, rho, Rtot
    real(8), intent(out) :: Dm(ns)
    real(8) :: T, Tdiff
    integer :: T_i, Tint(2)

    T = p/(rho*Rtot)
    T_i = idint(T)
    Tdiff = T-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    call co_DS_expr(rhoi,rho,Tint,Tdiff,p,Dm)

  endsubroutine co_DS


  !> Mixture-averaged diffusion coefficients from the (Tint,Tdiff) bracket, so callers
  !> that already hold the temperature indices (e.g. the diffusive flux) avoid recomputing.
  !>   Dm(s) = (1 - X_s) / sum_{j/=s} X_j / D_sj ,  with a near-pure-cell guard.
  !> The (1 - Y_s) Chemkin/Kee/Cantera variant is kept as commented reference in the
  !> final loop; edit the source there to switch the active formulation.
  !>
  !> Single sweep over the unique species pairs in the upper-triangular storage order of
  !> dij_tab, (1,2),(1,3),..,(1,ns),(2,3),..  Each binary D_ij is interpolated and inverted
  !> exactly once and scattered to both partners, so the cost is O(ns^2/2) with no ns x ns
  !> matrix and half the divisions of a naive double loop. The T-bracket is clamped once.
  pure subroutine co_DS_expr(rhoi,rho,Tint,Tdiff,pres,Dm)
    implicit none
    integer, intent(in)  :: Tint(2)
    real(8), intent(in)  :: rhoi(ns), rho, Tdiff, pres
    real(8), intent(out) :: Dm(ns)
    real(8), dimension(ns) :: Xi, Dm_den, Dsum
    real(8) :: Wmtot, Dval, invD, d0, d1
    integer :: s, i, j, p, Ti1, Ti2

    if (ns < 2) then
      Dm = 0.d0   ! single species: molecular diffusion is undefined / irrelevant
      return
    endif

    Wmtot = f_molecularWeight(rhoi)
    do s = 1, ns
      Xi(s) = rhoi(s)*Wmtot/(rho*Wm_tab(s))
    enddo

    ! Clamp the temperature bracket once: it is shared by every species pair.
    Ti1 = max(Tmin, min(Tint(1), Tmax))
    Ti2 = max(Tmin, min(Tint(2), Tmax))

    ! One pass over unique pairs; p tracks the dij_tab column (same storage order).
    Dm_den = 0.d0
    Dsum   = 0.d0
    p = 0
    do i = 1, ns
      do j = i+1, ns
        p = p + 1
        d0 = dij_tab(Ti1,p)
        d1 = dij_tab(Ti2,p)
        Dval = d0 + (d1-d0)*Tdiff
        invD = 1.d0 / Dval
        Dm_den(i) = Dm_den(i) + Xi(j)*invD   ! Wilke denominator (D_ij shared by i and j)
        Dm_den(j) = Dm_den(j) + Xi(i)*invD
        Dsum(i)   = Dsum(i) + Dval           ! running sum for the near-pure fallback
        Dsum(j)   = Dsum(j) + Dval
      enddo
    enddo

    ! Mixture-averaged numerator. Two literature forms share the same Wilke denominator
    ! Dm_den(s) = sum_{j/=s} X_j/D_sj; they differ only in the numerator:
    !   (1 - X_s)  Curtiss-Hirschfelder, mole-fraction form  -- ACTIVE below
    !   (1 - Y_s)  Chemkin/Kee, mass-fraction form (== Cantera mix_diff_coeffs)
    ! with Y_s = rhoi(s)/rho. The two diverge by ~10-20% per species in light/heavy
    ! mixtures (light species diffuse faster under the (1-Y_s) form). To match the
    ! Chemkin/Cantera standard instead, swap the active lines for the commented ones.
    do s = 1, ns
      if (1.d0 - Xi(s) < 1.d-10) then
        ! near-pure cell: D_s is ill-defined (0/0) but barely matters since grad(Y_s)
        ! is tiny and the flux is closed by the mass-conservation correction.
        ! Fall back to the mean binary diffusivity of species s.
        Dm(s) = Dsum(s) / dble(ns-1)
      else
        Dm(s) = (1.d0 - Xi(s)) / Dm_den(s)             ! (1 - X_s)  Curtiss-Hirschfelder [active]
        ! Dm(s) = (1.d0 - rhoi(s)/rho) / Dm_den(s)     ! (1 - Y_s)  Chemkin/Kee == Cantera mix_diff_coeffs
      endif
    enddo

    ! Correct from the table reference pressure to the local pressure: D ~ 1/p exactly
    ! (kinetic theory). The factor is uniform, so D_s scales linearly with it.
    Dm(1:ns) = Dm(1:ns) * ( dij_pref / pres )

  endsubroutine co_DS_expr


  pure function f_tab(p,rho,Rtot,s,tab) result(result)
    implicit none
    real(8), intent(in) :: rho, p, Rtot, tab(Tmin:Tmax,ns)
    integer, intent(in) :: s
    real(8) :: result
    real(8) :: T
    T = p/(Rtot*rho)
    ! chiamata alla funzione di interpolazione delle proprietà di tabella
    ! le proprietà per ora dipendono solo dalla temperatura, ma in generale
    ! possono anche dipendere dalla densità, motivo per cui in input alla funzione
    ! è anche data la variabile "ro"
    result = f_tabT(T,s,tab)
  endfunction f_tab


  pure function f_tabT(T,sp,tab) result(result)
    !$acc routine seq
    implicit none
    real(8), intent(in) :: T, tab(Tmin:Tmax,ns)
    integer, intent(in) :: sp
    real(8) :: result
    integer :: Tint(2)
    real(8) :: Tdiff
    
    Tint(1) = idint(T); Tint(2) = Tint(1)+1
    Tdiff = T-Tint(1)
    result = f_tabT_expr(sp,tab,Tint,Tdiff)

  endfunction f_tabT


  pure function f_tabT_expr(sp,tab,Tint,Tdiff) result(result)
    !$acc routine seq
    implicit none
    real(8), intent(in) :: Tdiff, tab(Tmin:Tmax,ns)
    integer, intent(in) :: sp, Tint(2)
    real(8) :: result
    real(8) :: Vij,Viij

    Vij=tab(Tint(1),sp)       ! int(T)   <- Tint(1)
    Viij=tab(Tint(2),sp)      ! int(T)+1 <- Tint(2)
    result = Vij+(Viij-Vij)*Tdiff

  endfunction f_tabT_expr

endmodule FLINT_Lib_Thermodynamic