module MOSE_IO_Ini
  use MOSE_IO_Ini_Open, only: MOSE_Load_Ini   ! one ini loader: MOSE_open dialect translated in memory
  
  implicit none

contains


  subroutine Read_Inifile ( config )
    use Finer, only: file_ini
    use MOSE_Types_m
    use MOSE_Parameters_m
    use MOSE_IO_Table, only: obj_tab
    use MOSE_IO_Chemistry, only: obj_chem
    use MOSE_IO_RANS, only: obj_rans
    use MOSE_IO_Solution
    use MOSE_Lib_RK, only: obj_rk
    use MOSE_Lib_IRS, only: obj_irs
    use MOSE_Lib_Omega_Relax, only: obj_omega_relax
    use MOSE_Lib_Preconditioning, only: obj_prec
    use MOSE_Lib_Extrapolation, only: wall_ghost_key, Assign_Extrap_Ghost, Assign_RANS_Wall_Ghost
    use MOSE_Lib_Newstate, only: Assign_State_Check
    use MOSE_Lib_Limiters, only: Assign_Limiter, SD_space_order
    use MOSE_Lib_Riemann, only: Assign_RS, SD_riemann
    use MOSE_Lib_Shock_Detector, only: SD, Assign_SD, sd_dev
    use MOSE_Lib_Newstate, only: Assign_Integration_Variables
    use MOSE_Mod_RANS, only: Set_RANS_Model
    use MOSE_Lib_Soot
    use MOSE_Lib_RotatingFrame, only: obj_rot

    implicit none
    type(config_type), intent(inout) :: config
    ! Local
    type(file_ini) :: fini
    integer :: error
    character(len=32) :: egk   ! [MOSE-scheme] extrapolation-ghost
    character(len=32) :: sdk   ! [MOSE-scheme] shock-detector


    ! Load input.ini
    call MOSE_Load_Ini(fini)

    ! Read input.ini
    call config%ini_setup(fini=fini)
    call obj_tab%ini_setup(fini=fini)
    call obj_chem%ini_setup(fini=fini)
    call obj_rans%ini_setup(fini=fini)
    call obj_io%ini_setup(fini=fini)
    call obj_rk%ini_setup(fini=fini)
    call obj_irs%ini_setup(fini=fini)
    call obj_prec%ini_setup(fini=fini)   ! before Assign_RS: HLLC-PC needs prec
    call obj_omega_relax%ini_setup(fini=fini)
    call Stop_On_Removed_Keys ( fini )
    ! The four defaults are MOSE_open's (Lib_Ghost, Ghost_Wall_Extrapolation, Check_And_Fix_State), whatever the dialect
    call fini%get(section_name=trim(codename)//'-scheme', option_name='wall-ghost', val=wall_ghost_key, error=error)
    if (error/=0) wall_ghost_key = 'cubic'
    egk = 'cubic'
    call fini%get(section_name=trim(codename)//'-scheme', option_name='extrapolation-ghost', val=egk, error=error)
    if (error/=0) egk = 'cubic'
    call Assign_Extrap_Ghost ( egk )   ! 400 faces: cubic (MOSE_open, default) or zero-gradient
    egk = 'copy'
    call fini%get(section_name=trim(codename)//'-scheme', option_name='rans-wall-ghost', val=egk, error=error)
    if (error/=0) egk = 'copy'
    call Assign_RANS_Wall_Ghost ( egk )   ! wall ghost of the RANS variables for their sources: copy (MOSE_open, default) or wall-value
    egk = 't'
    call fini%get(section_name=trim(codename)//'-scheme', option_name='state-check', val=egk, error=error)
    if (error/=0) egk = 't'
    call Assign_State_Check ( egk )   ! stop on p < 0 or NaN after every update (MOSE_open Check_And_Fix_State)
    call obj_soot%ini_setup(fini=fini)
    call obj_rot%ini_setup(fini=fini)
    
    ! Initialization
    call Input_Solution_Setup( config%newrun )
    call Assign_Limiter( config%limiter_name )
    call Assign_RS ( config%riemann_name, config%Minf, config%Mco )
    ! the MOSE_open shock sensor (Tramel / Chen): the key, else the one of the Riemann solver (Assign_Setup of MOSE_open)
    sdk = ''
    call fini%get(section_name=trim(codename)//'-scheme', option_name='shock-detector', val=sdk, error=error)
    if (error/=0) sdk = ''
    call Assign_SD ( sdk, config%riemann_name )
    if ( SD_riemann .or. SD_space_order) SD = .true.
    if ( trim(config%solver_name) == 'explicit' ) then
      call Assign_Integration_Variables( config%intvar_name )
    endif
    call Set_RANS_Model()

  end subroutine Read_Inifile


  subroutine Read_Inifile_Runtime ( config )
    use Finer, only: file_ini
    use MOSE_Types_m
    use MOSE_Parameters_m
    use MOSE_IO_Solution
    use MOSE_Lib_RK, only: obj_rk
    use MOSE_Lib_IRS, only: obj_irs
    use MOSE_Lib_Omega_Relax, only: obj_omega_relax
    use MOSE_Lib_Extrapolation, only: wall_ghost_key, Assign_Extrap_Ghost, Assign_RANS_Wall_Ghost
    use MOSE_Lib_Newstate, only: Assign_State_Check
    use MOSE_Lib_Limiters, only: Assign_Limiter, SD_space_order
    use MOSE_Lib_Riemann, only: Assign_RS, SD_riemann
    use MOSE_Lib_Shock_Detector, only: SD, Assign_SD, sd_dev
    use MOSE_Lib_Newstate, only: Assign_Integration_Variables
    
    implicit none
    type(config_type), intent(inout) :: config
    ! Local
    type(file_ini) :: fini
    integer :: error
    character(len=32) :: egk   ! [MOSE-scheme] extrapolation-ghost
    character(len=32) :: sdk   ! [MOSE-scheme] shock-detector

    ! Load input.ini
    call MOSE_Load_Ini(fini)

    ! Read input.ini
    call config%ini_setup(fini=fini)
    call obj_io%ini_setup(fini=fini)
    call obj_rk%ini_setup(fini=fini)
    call obj_irs%ini_setup(fini=fini)
    call obj_omega_relax%ini_setup(fini=fini)
    call Stop_On_Removed_Keys ( fini )
    ! The four defaults are MOSE_open's (Lib_Ghost, Ghost_Wall_Extrapolation, Check_And_Fix_State), whatever the dialect
    call fini%get(section_name=trim(codename)//'-scheme', option_name='wall-ghost', val=wall_ghost_key, error=error)
    if (error/=0) wall_ghost_key = 'cubic'
    egk = 'cubic'
    call fini%get(section_name=trim(codename)//'-scheme', option_name='extrapolation-ghost', val=egk, error=error)
    if (error/=0) egk = 'cubic'
    call Assign_Extrap_Ghost ( egk )   ! 400 faces: cubic (MOSE_open, default) or zero-gradient
    egk = 'copy'
    call fini%get(section_name=trim(codename)//'-scheme', option_name='rans-wall-ghost', val=egk, error=error)
    if (error/=0) egk = 'copy'
    call Assign_RANS_Wall_Ghost ( egk )   ! wall ghost of the RANS variables for their sources: copy (MOSE_open, default) or wall-value
    egk = 't'
    call fini%get(section_name=trim(codename)//'-scheme', option_name='state-check', val=egk, error=error)
    if (error/=0) egk = 't'
    call Assign_State_Check ( egk )   ! stop on p < 0 or NaN after every update (MOSE_open Check_And_Fix_State)

    ! Initialization
    call Assign_Limiter( config%limiter_name )
    call Assign_RS ( config%riemann_name, config%Minf, config%Mco )
    ! the MOSE_open shock sensor (Tramel / Chen): the key, else the one of the Riemann solver (Assign_Setup of MOSE_open)
    sdk = ''
    call fini%get(section_name=trim(codename)//'-scheme', option_name='shock-detector', val=sdk, error=error)
    if (error/=0) sdk = ''
    call Assign_SD ( sdk, config%riemann_name )
    if ( SD_riemann .or. SD_space_order) SD = .true.
    if ( trim(config%solver_name) == 'explicit' ) then
      call Assign_Integration_Variables( config%intvar_name )
    endif

  end subroutine Read_Inifile_Runtime


  !> Keys of the features this version does not have: multirate (branch multirate of DIMA-THRUST/MOSE), the
  !> line-implicit operator and the IMEX integrator (last in 85e3b63 of branch mose-launcher). A run that asks for
  !> one stops, rather than running the explicit scheme in silence; = f is accepted.
  subroutine Stop_On_Removed_Keys ( fini )
    use Finer, only: file_ini
    use MOSE_Parameters_m, only: codename
    implicit none
    type(file_ini), intent(inout) :: fini
    character(len=32) :: v
    integer :: error
    ! (the key names are literal: the launcher lists the keys the code reads from the source)
    v = 'f'
    call fini%get(section_name=trim(codename)//'-scheme', option_name='multirate', val=v, error=error)
    if ( error == 0 ) call stop_unless_off ( 'multirate', v, 'multirate (branch multirate of DIMA-THRUST/MOSE)' )
    v = 'f'
    call fini%get(section_name=trim(codename)//'-scheme', option_name='j-implicit', val=v, error=error)
    if ( error == 0 ) call stop_unless_off ( 'j-implicit', v, 'line-implicit operator (last in 85e3b63)' )
    v = 'f'
    call fini%get(section_name=trim(codename)//'-scheme', option_name='imex', val=v, error=error)
    if ( error == 0 ) call stop_unless_off ( 'imex', v, 'IMEX integrator (last in 85e3b63)' )

  contains

    subroutine stop_unless_off ( key, val, what )
      character(len=*), intent(in) :: key, what
      character(len=*), intent(inout) :: val
      integer :: i
      do i = 1, len_trim(val)
        if ( val(i:i) >= 'A' .and. val(i:i) <= 'Z' ) val(i:i) = achar( iachar(val(i:i)) + 32 )
      end do
      select case ( trim(val) )
      case ('f', '.false.', 'false')
      case default
        write(*,'(7A)') ' ERROR: [', trim(codename), '-scheme] ', key, ' = ', trim(val), &
                        ' -- this version of the code has no ' // what // ': remove the key'
        error stop 'a removed feature is requested in input.ini (see the ERROR line)'
      end select
    end subroutine stop_unless_off

  end subroutine Stop_On_Removed_Keys

end module MOSE_IO_Ini
