module MOSE_Lib_Soot
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Parameters_m

  implicit none
  
  real(R8), parameter :: Mc = 12.011_R8              ! kg/kmol
  real(R8), parameter :: rhos = 1800.0_R8
  integer :: i_C2H2, i_CO2, i_CO, i_H2, i_O2, i_CS, i_C6H6 ! indici specie utili

  type :: soot_type
    real(R8) :: rho_soot
    logical :: is_on, LL91, LIN
    character(llen) :: model_name
    contains
      procedure, pass(self) :: ini_setup => read_ini
  end type soot_type

  type(soot_type), public :: obj_soot
  private :: read_ini

contains


  subroutine read_ini ( self, fini )
    use Finer, only: file_ini
    use MOSE_Parameters_m, only: codename
    use MOSE_Global_m, only: nsoot

    implicit none
    class(soot_type) :: self
    type(file_ini) :: fini
    integer :: error
    ! === Per lettura delle specie ===
    character(len=256) :: line, species
    real(R8) :: mw
    integer :: ios, idx, unit
    character(len=*), parameter :: phasefile = 'INPUT/phase.txt'

    call fini%get(section_name=trim(codename)//'-SOOT', option_name='model', val=self%model_name, error=error)
    if (error/=0) then
      self%is_on = .false.
      nsoot = 0
      return
    else
      self%is_on = .true.
      if (index(self%model_name,'LL91')>0) then
        self%LL91 = .true.
        self%LIN = .false.
      elseif (index(self%model_name,'LIN')>0) then
        self%LL91 = .false.
        self%LIN = .true.
      else
        error stop ( ' Invalid soot model.' )
      end if
    end if

    ! inizializzo
    i_O2   = -1
    i_C2H2 = -1
    i_CO2  = -1
    i_CO   = -1
    i_H2   = -1
    i_CS   = -1
    i_C6H6 = -1

    ! === Scansione del file phase.txt ===
    open(newunit=unit, file=phasefile, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'WARNING: cannot open ', trim(phasefile)
      return
    end if

    ! salta prima riga ("ideal-gas phase")
    read(unit, '(A)', iostat=ios) line

    idx = 0
    do
      read(unit, *, iostat=ios) species, mw
      if (ios /= 0) exit
      idx = idx + 1

      select case (adjustl(trim(species)))
      case ('O2')
        i_O2 = idx
      case ('C2H2')
        i_C2H2 = idx
      case ('CO2')
        i_CO2 = idx
      case ('CO')
        i_CO = idx
      case ('H2')
        i_H2 = idx
      case ('C(gr)')
        i_CS = idx
      case ('C6H6')
        i_C6H6 = idx
      end select
    end do
    close(unit)

    ! === Output di diagnostica ===
    write(*,'(A)')    ' Soot model species indices'
    write(*,'(A,I0)') '  O2   : ', i_O2
    write(*,'(A,I0)') '  C2H2 : ', i_C2H2
    write(*,'(A,I0)') '  C6H6 : ', i_C6H6
    write(*,'(A,I0)') '  CO2  : ', i_CO2
    write(*,'(A,I0)') '  CO   : ', i_CO
    write(*,'(A,I0)') '  H2   : ', i_H2
    write(*,'(A,I0)') '  C(gr): ', i_CS
    write(*,*)
    
    if (i_CS==-1) then
      self%is_on = .true.
      nsoot = 2
    else
      self%is_on = .true.
      nsoot = 1
    end if
    
  end subroutine read_ini


  subroutine LL91_soot_rates(rho, T, C_C2H2, C_O2, Y_soot, Nvol, &
                             R1, R2, R3, R4)
    implicit none
    real(R8), intent(in)  :: rho, T, C_C2H2, C_O2, Y_soot, Nvol
    real(R8), intent(out) :: R1, R2, R3, R4

    ! Constants
    real(R8), parameter :: pi   = 3.141592653589793_R8
    real(R8), parameter :: NA   = 6.02214076e26_R8
    real(R8), parameter :: kB   = 1.380649e-23_R8
    integer,  parameter :: Cmin = 100
    real(R8), parameter :: Ca   = 9.0_R8

    ! Locals
    real(R8) :: C_mol, d_p, S_area, f_Area
    real(R8) :: k1, k2, k3, k4

    ! Molar conc. of soot carbon
    C_mol = rho*max(Y_soot,0.0_R8)/Mc  ! kmol/m^3

    ! Diameter from volume-based number density
    d_p = (6.0_R8*Mc*C_mol / (pi*rhos*Nvol))**(1.0_R8/3.0_R8)

    ! Surface area per volume
    S_area = pi * d_p**2 * Nvol
    f_Area = sqrt(pi*((6.0_R8*Mc)/(pi*rhos))**(2.0_R8/3.0_R8)) * C_mol**(1.0_R8/3.0_R8) * Nvol**(1.0_R8/6.0_R8)

    ! --- Rate constants (LL91) ---
    k1 = 0.1e5_R8 * exp(-21100.0_R8/T)
    k2 = 0.6e4_R8 * exp(-12100.0_R8/T)
    k3 = 0.1e5_R8 * sqrt(T) * exp(-19680.0_R8/T)
    k4 = 2.0_R8 * Ca * ((6.0_R8*Mc)/(pi*rhos))**(1.0_R8/6.0_R8) * sqrt((6.0_R8*kB*T)/rhos)

    ! --- Reaction rates ---
    R1 = k1 * C_C2H2                   ! nucleation
    R2 = k2 * C_C2H2 * f_Area          ! surface growth
    R3 = k3 * S_area * C_O2            ! oxidation
    R4 = (2.0_R8 * NA / Cmin) * R1 &   ! Nucleation-Coagulation [#/m^3/s]
         - k4 * C_mol**(1.0_R8/6.0_R8) * Nvol**(11.0_R8/6.0_R8) 

  end subroutine LL91_soot_rates

  subroutine LIN_soot_rates(rho, T, C_C2H2, C_O2, Y_soot, Nvol, &
                            R1, R2, R3, R4)
    implicit none
    real(R8), intent(in)  :: rho, T, C_C2H2, C_O2, Y_soot, Nvol
    real(R8), intent(out) :: R1, R2, R3, R4

    ! Constants
    real(R8), parameter :: pi   = 3.141592653589793_R8
    real(R8), parameter :: NA   = 6.02214076e26_R8
    real(R8), parameter :: kB   = 1.380649e-23_R8
    integer,  parameter :: Cmin = 100
    real(R8), parameter :: Ca   = 9.0_R8
    real(R8), parameter :: Mc   = 12.011_R8     ! [kg/kmol]
    real(R8), parameter :: rhos = 1800.0_R8     ! soot density [kg/m^3]

    ! Locals
    real(R8) :: C_mol, d_p, S_area
    real(R8) :: k1, k2, k3, k4

    ! Molar concentration of soot carbon
    C_mol = rho * max(Y_soot, 0.0_R8) / Mc   ! kmol/m^3

    ! Mean particle diameter from number density
    d_p = (6.0_R8 * Mc * C_mol / (pi * rhos * Nvol))**(1.0_R8/3.0_R8)

    ! Surface area per unit volume
    S_area = pi * d_p**2 * Nvol

    ! --- Rate constants ---
    ! Lindstedt (2005) nucleation
    k1 = 0.63e4_R8 * exp(-21100.0_R8 / T)
    ! Lindstedt (1994) surface growth (Bockhorn eq. 27.36)
    k2 = 0.1e-11_R8 * exp(-12100.0_R8 / T)
    ! LL91 oxidation
    k3 = 0.1e5_R8 * sqrt(T) * exp(-19680.0_R8 / T)
    ! LL91 coagulation
    k4 = 2.0_R8 * Ca * ((6.0_R8*Mc)/(pi*rhos))**(1.0_R8/6.0_R8) * sqrt((6.0_R8*kB*T)/rhos)

    ! --- Reaction rates ---
    R1 = k1 * C_C2H2                   ! nucleation [kmol/m^3/s]
    R2 = k2 * C_C2H2  * Nvol           ! surface growth [kmol/m^3/s]
    R3 = k3 * S_area * C_O2            ! Oxidation (LL91)
    R4 = (2.0_R8 * NA / Cmin) * R1 &   ! Nucleation-Coagulation [#/m^3/s]
         - k4 * C_mol**(1.0_R8/6.0_R8) * Nvol**(11.0_R8/6.0_R8) 

  end subroutine LIN_soot_rates


  ! subroutine HACA_soot_rates(rho, T, C_C2H2, C_O2, C_H, C_OH, C_H2, C_H2O, Y_soot, Nvol, &
  !                            R1, R2, R3, R4)
  !   implicit none
  !   real(R8), intent(in)  :: rho, T
  !   real(R8), intent(in)  :: C_C2H2, C_O2, C_H, C_OH, C_H2, C_H2O
  !   real(R8), intent(in)  :: Y_soot, Nvol
  !   real(R8), intent(out) :: R1, R2, R3, R4

  !   ! -------------------------------------------------------------------
  !   ! HACA soot model (Appel-Bockhorn-Frenklach, Combust. Flame 121:122–136, 2000)
  !   ! - Nucleation: Lindstedt (2005)
  !   ! - Growth:     HACA (Frenklach & Wang, 1990; Balthasar & Frenklach, 2005)
  !   ! - Oxidation:  HACA (OH + O2 contributions)
  !   ! - Coagulation: -
  !   ! -------------------------------------------------------------------

  !   ! Constants
  !   real(R8), parameter :: pi   = 3.141592653589793_R8
  !   real(R8), parameter :: NA   = 6.02214076e26_R8
  !   real(R8), parameter :: kB   = 1.380649e-23_R8
  !   integer,  parameter :: Cmin = 100
  !   real(R8), parameter :: Ca   = 9.0_R8
  !   real(R8), parameter :: Mc   = 12.011_R8       ! [kg/kmol]
  !   real(R8), parameter :: rhos = 1800.0_R8       ! [kg/m^3]
  !   real(R8), parameter :: Rg   = 8314.462618_R8  ! [J/kmol/K]

  !   ! Locals
  !   real(R8) :: C_mol
  !   real(R8) :: k1, k3, k4
  !   real(R8) :: RT_kcal, chi_CH, a_param, b_param, alpha, denom, chi_S
  !   real(R8) :: fR1, rR1, fR2, rR2, fR3, fR4, fR5
  !   real(R8) :: rxns_s, rgrowth, RO2, ROH

  !   ! Molar concentration of soot carbon
  !   C_mol = rho * max(Y_soot, 0.0_R8) / Mc   ! kmol/m^3

  !   ! --- Nucleation (Lindstedt 2005)
  !   k1 = 0.63e4_R8 * exp(-21100.0_R8 / T)
  !   R1 = k1 * C_C2H2   ! [kmol/m^3/s]

  !   ! -------------------------------------------------------------------
  !   ! --- HACA growth (Frenklach & Wang 1990, Appel et al. 2000)
  !   ! -------------------------------------------------------------------
  !   RT_kcal = 1.9872036e-3_R8 * T         ! kcal/mol
  !   chi_CH  = 2.3e15_R8                   ! C–H sites/cm^2
  !   a_param = 33.167_R8 - 0.0154_R8 * T
  !   b_param = -2.5786_R8 + 0.00112_R8 * T

  !   ! Raw rates [cm^3/mol/site/s]
  !   fR1 = 4.2e13_R8 * exp(-13.0_R8 / RT_kcal) * C_H    * 1.0e-3_R8
  !   rR1 = 3.9e12_R8 * exp(-11.0_R8 / RT_kcal) * C_H2   * 1.0e-3_R8
  !   fR2 = 1.0e10_R8 * T**0.734_R8 * exp(-1.43_R8 / RT_kcal) * C_OH   * 1.0e-3_R8
  !   rR2 = 3.68e8_R8 * T**1.139_R8 * exp(-17.1_R8 / RT_kcal) * C_H2O  * 1.0e-3_R8
  !   fR3 = 2.0e13_R8 * C_H  * 1.0e-3_R8
  !   fR4 = 8.00e7_R8 * T**1.56_R8 * exp(-3.8_R8 / RT_kcal) * C_C2H2 * 1.0e-3_R8
  !   fR5 = 2.2e12_R8 * exp(-7.5_R8 / RT_kcal) * C_O2 * 1.0e-3_R8

  !   ! Steady-state soot radical site concentration chi_S [sites/cm^2]
  !   denom = rR1 + rR2 + fR3 + fR4 + fR5
  !   chi_S = chi_CH * (fR1 + fR2) / denom

  !   ! Alpha (fraction of active surface sites)
  !   alpha = tanh(a_param / log10((rho*Y_soot)/Nvol) + b_param)
  !   if (alpha < 0.0_R8) alpha = 1.0_R8

  !   ! Growth rate per surface [rxns/s·m²]
  !   rxns_s = fR4 * alpha * chi_S * 1.0e4_R8     ! [1/cm²/s] -> [1/m²/s]
  !   rgrowth = 2.0_R8 * rxns_s * Mc / NA         ! [kg/m²/s]
  !   R2 = 2.0_R8 * rxns_s / Na                   ! [kmol/m³/s]  (proportional to Nvol)

  !   ! -------------------------------------------------------------------
  !   ! --- HACA oxidation (O2 + OH)
  !   ! -------------------------------------------------------------------
  !   RO2 = 2.0_R8 * rxns_s * Mc / NA             ! same structure as growth for O2 sites
  !   ROH = 0.13_R8 * C_OH * sqrt(Rg*T/(2.0_R8*pi*17.0_R8)) * Mc  ! simplified term
  !   R3 = (RO2 + ROH) * Nvol / Mc                ! convert to [kmol/m³/s]

  !   ! --- Nucleation (same as LL91)
  !   R4 = (2.0_R8 * NA / Cmin) * R1

  ! end subroutine HACA_soot_rates


  subroutine soot_source_terms ( domain )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle

      call soot_source_blk ( domain % blk(b) % P,            &
                             domain % blk(b) % r,            &
                             domain % blk(b) % vol,          &
                             domain % blk(b) % dim           )
    
    end do

  end subroutine soot_source_terms


  subroutine soot_source_blk ( prim, res, vol, n )
    use MOSE_General_Types_m
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic

    implicit none
    integer, intent(in) :: n(3)
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: res
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: vol
    ! Local
    integer :: i, j, k, s
    real(R8) :: rho, Rgas, T
    real(R8) :: Yc2h2, Yo2, C_c2h2, C_o2, Y_soot, Nvol
    real(R8) :: R1, R2, R3, R4, S_rhoYs, S_rhoN, omega(nsc), energy


    !$omp do collapse (3) private ( rho, Rgas, T, Yc2h2, Yo2, Y_soot, Nvol, C_c2h2, C_o2 ) &
    !$omp private ( R1, R2, R3, R4, S_rhoYs, S_rhoN, omega, energy, i, j, k, s )
    
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)

      call Co_rotot_Rtot ( Prim(1:nsc,i,j,k), rho, Rgas )
      T = prim(np,i,j,k) / (rho * Rgas)

      ! --- dati specie / soot ---
      Yc2h2  = prim(i_C2H2,i,j,k) / rho
      Yo2    = prim(i_O2  ,i,j,k) / rho
      if (nsoot==1) then
        Y_soot = prim(i_CS,i,j,k) / rho        ! soot mass fraction
        Nvol   = prim(nc  ,i,j,k)              ! number per volume [#/m^3]
      else
        Y_soot = prim(nc  ,i,j,k) / rho        ! soot mass fraction
        Nvol   = prim(nc+1,i,j,k)              ! number per volume [#/m^3]
      end if

      ! concentrazioni molari (kmol/m^3)
      C_c2h2 = rho * Yc2h2 / Wm_tab(i_C2H2)
      C_o2   = rho * Yo2   / Wm_tab(i_O2)

      ! --- reaction rate modello LL91 (R1,R2,R3 in kmol/m^3/s di C, R4 in #/m^3/s) ---
      if (obj_soot%LL91) then
        call LL91_soot_rates(rho, T, C_c2h2, C_o2, Y_soot, Nvol, R1, R2, R3, R4)
      elseif (obj_soot%LIN) then
        call LIN_soot_rates(rho, T, C_c2h2, C_o2, Y_soot, Nvol, R1, R2, R3, R4)
      endif

      ! --- sorgenti per i residui ---
      ! massa soot: d(rho*Ys)/dt = Mc*(2R1 + R2 - R3)  [kg/m^3/s]
      S_rhoYs = Mc * ( 2.0*R1 + R2 - R3 )

      ! numero particelle per volume: d(rhoN)/dt = R4  [#/m^3/s]
      S_rhoN  = R4

      ! Accumulo nei residui
      if (nsoot==2) then
        res(nc  ,i,j,k) = res(nc  ,i,j,k) - S_rhoYs * vol(i,j,k)
        res(nc+1,i,j,k) = res(nc+1,i,j,k) - S_rhoN  * vol(i,j,k)
      else
        omega = 0.0_R8
        ! mass sources, written as production (so negative for consumption)
        omega(i_C2H2) = -(R1+R2) * Wm_tab(i_C2H2)
        omega(i_O2  ) = -0.5*R3  * Wm_tab(i_O2  )
        omega(i_H2  ) = +(R1+R2) * Wm_tab(i_H2  )
        omega(i_CO  ) = +R3      * Wm_tab(i_CO  )
        omega(i_CS  ) = -sum(omega) ! enforcing mass conservation
        ! update species residuals
        res(1:nsc,i,j,k) = res(1:nsc,i,j,k) - omega(1:nsc)*vol(i,j,k)
        ! update number density residual
        res(nc   ,i,j,k) = res(nc   ,i,j,k) - S_rhoN*vol(i,j,k) ! number density
        ! Energy source term
        energy = 0.0_R8
        do s = 1, nsc
          energy = energy + omega(s)*f_tabT(T,s,h_tab)
        end do
        res(np,i,j,k) = res(np,i,j,k) - energy*vol(i,j,k)
      end if

    enddo ; enddo ; enddo ! (i, j, k) loop

  end subroutine soot_source_blk




end module MOSE_Lib_Soot
