!>@brief: Module for legacy gas-surface interaction subroutines.
module MOSE_Mod_GSI
  use MOSE_IO_Ini_Open, only: MOSE_Load_Ini   ! one ini loader: MOSE_open dialect translated in memory
  use MOSE_Global_m

  implicit none

  ! Paraffin/Melting
  real(kind=8) :: T_i = 298.15d0             ! initial temperature of paraffin
  real(kind=8) :: cp_wax = 1946.0310193884d0 ! data taken from NIST website (last accessed Fri Nov 1 2019): https://webbook.nist.gov/cgi/cbook.cgi?ID=C544854&Units=SI&Mask=6EF
  real(kind=8) :: Dh_wax = 1.698285789316d+5 ! data taken from NIST website (last accessed Fri Nov 1 2019): https://webbook.nist.gov/cgi/cbook.cgi?ID=C544854&Units=SI&Mask=6EF 
  real(kind=8) :: y_SEBS = 0d0               ! SEBS mass fraction in fuel
  integer :: iFUEL = 10                      ! standard set for singhWC32
  integer :: iSEBS = 2

  ! Carbon
  real(kind=8), parameter :: wmc  = 12.0107d0 ! Carbon molecular weight
  real(kind=8), parameter :: pref = 1.01325d5 ! Athmospheric reference pressure 

  ! HDPE
  real(kind=8) :: cp_HDPE = 1255.2d0
  real(kind=8) :: Dh_HDPE = 2.72d+6
  real(kind=8) :: A_HDPE = 4588.8d3
  real(kind=8) :: Ea_HDPE = 251039.73d3
  real(kind=8) :: Ti_HDPE = 298.15d0

  ! PP
  real(kind=8) :: cp_PP = 1700d0
  real(kind=8) :: Dh_PP = 2.4823d6
  real(kind=8) :: A_PP = 2.12d15
  real(kind=8) :: Ea_PP = 2.12d8
  real(kind=8) :: Ti_PP = 298.15d0

  ! HTPB
  real(kind=8) :: cp_HTPB = 1632d0
  real(kind=8) :: Dh_HTPB = 1.1d6
  real(kind=8) :: Ti_HTPB = 298.15d0
  real(kind=8) :: A1_HTPB = 3.965d0
  real(kind=8) :: Ea1_HTPB = 55.8564d6
  real(kind=8) :: A2_HTPB = 11.04d-3
  real(kind=8) :: Ea2_HTPB = 20.54344d6
  real(kind=8) :: Ts_HTPB = 722d0
  real(kind=8) :: rho_HTPB = 960d0

contains

  subroutine Setup_GSI()  ! TODO: adattarlo al nuoco IO_Ini
    use Finer, only: file_ini
    implicit none
    integer :: error
    type(file_ini) :: fini

    call MOSE_Load_Ini(fini)
    ! HDPE
    call fini%get(section_name='HDPE-prop', option_name='cp', val=cp_HDPE, error=error)
    call fini%get(section_name='HDPE-prop', option_name='Dh', val=Dh_HDPE, error=error)
    call fini%get(section_name='HDPE-prop', option_name='A', val=A_HDPE, error=error)
    call fini%get(section_name='HDPE-prop', option_name='Ea', val=Ea_HDPE, error=error)
    call fini%get(section_name='HDPE-prop', option_name='T_i', val=Ti_HDPE, error=error)
    ! PP
    call fini%get(section_name='PP-prop', option_name='cp', val=cp_PP, error=error)
    call fini%get(section_name='PP-prop', option_name='Dh', val=Dh_PP, error=error)
    call fini%get(section_name='PP-prop', option_name='A', val=A_PP, error=error)
    call fini%get(section_name='PP-prop', option_name='Ea', val=Ea_PP, error=error)
    call fini%get(section_name='PP-prop', option_name='T_i', val=Ti_PP, error=error)
    ! Paraffin
    call fini%get(section_name='melting-prop', option_name='cp', val=cp_wax, error=error)
    call fini%get(section_name='melting-prop', option_name='Dh', val=Dh_wax, error=error)
    call fini%get(section_name='melting-prop', option_name='T_i', val=T_i, error=error)
    call fini%get(section_name='melting-prop', option_name='y_SEBS', val=y_SEBS, error=error)
    call fini%get(section_name='melting-prop', option_name='iFuel', val=iFUEL, error=error)
    call fini%get(section_name='melting-prop', option_name='iSEBS', val=iSEBS, error=error)
    ! HTPB
    call fini%get(section_name='HTPB-prop', option_name='cp', val=cp_HTPB, error=error)
    call fini%get(section_name='HTPB-prop', option_name='Dh', val=Dh_HTPB, error=error)
    call fini%get(section_name='HTPB-prop', option_name='A1', val=A1_HTPB, error=error)
    call fini%get(section_name='HTPB-prop', option_name='A2', val=A2_HTPB, error=error)
    call fini%get(section_name='HTPB-prop', option_name='Ea1', val=Ea1_HTPB, error=error)
    call fini%get(section_name='HTPB-prop', option_name='Ea2', val=Ea2_HTPB, error=error)
    call fini%get(section_name='HTPB-prop', option_name='Ts', val=Ts_HTPB, error=error)
    call fini%get(section_name='HTPB-prop', option_name='T_i', val=Ti_HTPB, error=error)
    call fini%get(section_name='HTPB-prop', option_name='rho', val=rho_HTPB, error=error)

  end subroutine Setup_GSI


  subroutine GSI_Melting ( qw, Tw, omega, mdot )

    implicit none
    real(kind=8), intent(in)    :: qw, Tw
    real(kind=8), intent(out)   :: omega(nsc), mdot
    real(kind=8) :: cp_solid, Dh_melting, Dh_conduction, heat_absorbed, mdot_wax

    cp_solid   = cp_wax
    dh_melting = dh_wax

    ! energy balance
    dh_conduction = cp_solid*(Tw - T_i)  ! heat of conduction in solid (J/kg)  
    heat_absorbed = dh_melting + dh_conduction  

    ! mass flux computed as a result of the energy balance at the grain surface
    mdot = qw / heat_absorbed        ! (kg/(s m^2))

    mdot_wax = (1.0d0-y_SEBS)*mdot
    
    !% Species source term
    omega = 0.0d0                     
    ! Injection of paraffin as C32H66 and ethylene for the SEBS
    if (y_SEBS>0.0d0) omega(iSEBS)  = y_SEBS*mdot ! C2H4
    omega(iFUEL) = mdot_wax                       ! C32H66
    
  end subroutine GSI_Melting


  subroutine GSI_Carbon ( roi, Tw, omegadot, mdot )
    use FLINT_Lib_Thermodynamic, only: Runiv, Ri_tab, Wm_tab

    implicit none
    real(kind=8), intent(in)    ::  roi(nsc), Tw
    real(kind=8), intent(inout) :: omegadot(nsc), mdot
    real(kind=8) :: dotm_H2O, dotm_CO2, dotm_O2, dotm_OH, dotm_O
    real(kind=8) :: k5, k6, k7, k8, kH2O, kCO2, kOH, kO, Y_term
    real(kind=8) :: pi(nsc), RuTw, inv_sqrtTw
    integer, parameter :: iH=1, iC2H4=2, iOH=3, iCO=4, iCO2=5
    integer, parameter :: iH2=6, iH2O=7, iO2=8, iO=9

    omegadot = 0d0  ! initialization

    ! Constants
    RuTw = Runiv*Tw
    inv_sqrtTw = 1.0d0/sqrt(Tw)

    ! Partial pressure in athmospheres
    pi = roi*Tw*Ri_tab/pref

    ! Rates
    kH2O = 4.80d5 *exp(-288.00d6/RuTw)  ! H2O
    kCO2 = 9.00d3 *exp(-285.00d6/RuTw)  ! CO2
    kOH  = 3.61d2 *inv_sqrtTw
    kO   = 6.655d2*inv_sqrtTw

    k5   = 2.40d3 *exp(-125.60d6/RuTw)  ! O2a
    k6   = 2.13d1 *exp(  17.17d6/RuTw)  ! O2b
    k7   = 5.35d-1*exp( -63.64d6/RuTw)  ! O2c
    k8   = 1.81d7 *exp(-406.10d6/RuTw)  ! O2d

    ! Compute ybig term
    Y_term = 1d0 / ( 1d0 + k8/(k7*pi(iO2)) )

    dotm_H2O = kH2O * sqrt(pi(iO2))
    dotm_CO2 = kCO2 * sqrt(pi(iCO2))
    dotm_O2  = Y_term * (k5*pi(iO2))/(1.0d0 + k6*pi(iO2)) + k7*pi(iO2)*(1d0-Y_term)
    dotm_OH  = kOH*pi(iOH)
    dotm_O   = kO*pi(iO)
    ! Compute source terms
    omegadot(iO2)  = -(Wm_tab(iO2)/wmc)*(0.5d0*dotm_O2)                                         ! O2
    omegadot(iH2O) = -(Wm_tab(iH2O)/wmc)*(dotm_H2O)                                             ! H2O
    omegadot(iCO)  =  (Wm_tab(iCO)/wmc)*(2.d0*dotm_CO2 + dotm_H2O + dotm_O2 + dotm_OH + dotm_O) ! CO
    omegadot(iCO2) = -(Wm_tab(iCO2)/wmc)*(dotm_CO2)                                             ! CO2
    omegadot(iH2)  =  (Wm_tab(iH2)/wmc)*(dotm_H2O)                                              ! H2
    omegadot(iO)   = -(Wm_tab(iO)/wmc)*(dotm_O)                                                 ! O
    omegadot(iH)   =  (Wm_tab(iH)/wmc)*(dotm_OH)                                                ! H
    omegadot(iOH)  = -(Wm_tab(iOH)/wmc)*(dotm_OH)                                               ! OH

    ! output of dotm
    mdot = sum ( omegadot )

  end subroutine GSI_carbon


  subroutine GSI_HDPE ( Tw, omegadot, mdot )
    use FLINT_Lib_Thermodynamic, only: Runiv

    implicit none
    real(kind=8), intent(in) :: Tw
    real(kind=8), intent(inout) :: omegadot(nsc), mdot

    omegadot = 0d0
    mdot = A_HDPE * exp(-Ea_HDPE/(2*Runiv*Tw))
    omegadot(2) = mdot ! ethylene, C2H4

  end subroutine GSI_HDPE


  subroutine GSI_HTPB ( Tw, omegadot, mdot )
    use FLINT_Lib_Thermodynamic, only: Runiv

    implicit none
    real(kind=8), intent(in) :: Tw
    real(kind=8), intent(inout) :: omegadot(nsc), mdot

    omegadot = 0d0
    if (Tw<=Ts_HTPB) then
      mdot = rho_HTPB * A1_HTPB * exp(-Ea1_HTPB/(Runiv*Tw))
    else
      mdot = rho_HTPB * A2_HTPB * exp(-Ea2_HTPB/(Runiv*Tw))
    end if
    omegadot(2) = mdot ! C4H6

  end subroutine GSI_HTPB


  subroutine GSI_PP ( Tw, omegadot, mdot )
    use FLINT_Lib_Thermodynamic, only: Runiv

    implicit none
    real(kind=8), intent(in) :: Tw
    real(kind=8), intent(inout) :: omegadot(nsc), mdot
    real(kind=8) :: Arrh, RR

    omegadot = 0d0
    Arrh = Ea_PP/(Runiv*Tw)
    RR = sqrt ( 4.60517d0 * ( 1d0 - Ti_PP/Tw + Dh_PP/(cp_PP*Tw) ) -Dh_PP/(cp_PP*Tw) )
    mdot = 910d0*sqrt( A_PP * exp(-Arrh) * 1d0/Arrh * 0.2d0/(910d0*cp_PP) * 1/RR )
    omegadot(2) = mdot ! propylene, C3H6

  end subroutine GSI_PP


end module MOSE_Mod_GSI
