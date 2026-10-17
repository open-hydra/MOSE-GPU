module MOSE_IO_chemistry
  use MOSE_IO_Ini_Open, only: open_force_frozen
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Parameters_m

  implicit none

  type :: chem_type
    logical :: reactive
    character(llen) :: mechanism_name, ode_name
    character(llen) :: ratefile, stoichfile
    integer, allocatable :: no_chem_list(:)          ! Blocks without chemistry
    contains
      procedure, pass(self) :: ini_setup => read_ini
  end type chem_type

  type(chem_type), public :: obj_chem
  private :: read_ini

contains

  subroutine read_ini(self, fini)
    use Finer, only: file_ini
    use MOSE_Lib_Strang, only: set_StrangSplitting
    use FLINT_Load_chemistry
    use FLINT_Lib_Chemistry_wdot, only: Assign_Mechanism, chemistry_jacobian
    use FLINT_Lib_Chemistry_rhs,  only: use_analytical_jacobian, validate_analytical_jacobian
    use FLINT_Lib_Radau5_dev,     only: flint_acc_upload_tables, flint_acc_upload_thermo, set_radau5_dev_tols
    use FLINT_Lib_Thermodynamic,  only: ns, Tmin, Tmax, dij_tab, dij_pref
    use MOSE_Lib_Limiters,        only: SD_space_order, lim_dev
    use MOSE_Lib_RANS,            only: eddy_dev, ransdiff_dev, stress_dev, sc_dev, sct_dev, prt_dev, prl_dev, mcdiff_dev
    use MOSE_Lib_Riemann,         only: rs_dev, mco_dev
    use MOSE_Lib_Shock_Detector,  only: sd_dev
    use MOSE_Global_m, only : nsc, RANS, VISC, nu, nv, nw, np, nt, nrans, nprim, gc
    use oslo, only: setup_odesolver
    use strings, only: parse
    implicit none
    ! Local
    type(file_ini) :: fini
    class(chem_type) :: self
    integer :: max_num_steps_odesolve, neq, iopt(3)
    real(R8), allocatable :: RT(:), AT(:)
    real(R8) :: rtol_sp, rtol_pre, atol_sp, atol_pre
    integer :: error_mechanism, error
    character(llen) :: nochem_str
    integer :: ios
    character(llen) :: args(100)
    integer :: i, nochem_blocks

    self%reactive = .false.

    if (open_force_frozen) then   ! [MOSE-Physics] chemistry = frozen (MOSE_open dialect): the mechanism file is not read
      error_mechanism = 1
    else
      error_mechanism = read_chemistry(mech_name = self%mechanism_name)
    end if

    ! if (error_mechanism/=0) then 
    !   call fini%get(section_name='Chemistry', option_name='mechanism-name', val=self%mechanism_name, error=error_mechanism)
    !   call read_chemistry_file_legacy()
    ! endif

    ! ODE solver setup.
    if (error_mechanism==0) then

      self%reactive = .true.
      call Assign_Mechanism(self%mechanism_name)

      ! ODE solver selection.
      call fini%get(section_name= trim(codename)//'-ODE', option_name='solver', val=self%ode_name, error=error)
      if (error/=0) then
        write(*,*) ' Warning: MOSE-ODE not specified.'
        self%ode_name='H-radau5'
      end if
      ! ODE solver parameters.
      iopt = 0
      neq = nsc+1
      allocate(RT(neq),AT(neq))

      call fini%get(section_name= trim(codename)//'-ODE', option_name='max-steps-ode', val=max_num_steps_odesolve, error=error)
      if (error/=0) iopt(1)=100000

      call fini%get(section_name= trim(codename)//'-ODE', option_name='relative-tol-species', val=rtol_sp, error=error)
      if (error/=0) rtol_sp=1d-5
      RT(1:nsc) = rtol_sp
      
      call fini%get(section_name= trim(codename)//'-ODE', option_name='relative-tol-pressure', val=rtol_pre, error=error)
      if (error/=0) rtol_pre=1d-5
      RT(neq) = rtol_pre

      call fini%get(section_name= trim(codename)//'-ODE', option_name='absolute-tol-species', val=atol_sp, error=error)
      if (error/=0) atol_sp=1d-5
      AT(1:nsc) = atol_sp
      
      call fini%get(section_name= trim(codename)//'-ODE', option_name='absolute-tol-pressure', val=atol_pre, error=error)
      if (error/=0) atol_pre=1d-5
      AT(neq) = atol_pre

      call setup_odesolver(N=neq,solver=self%ode_name,RT=RT,AT=AT,iopt=iopt)
      !> GPU chemistry (OpenACC build): mirror the tolerances for the device
      !> cell loop and upload the read-only thermo/rate tables once.
      call set_radau5_dev_tols(neq, RT, AT)
      call flint_acc_upload_tables()

      !> Analytical-Jacobian opt-in. Default off; honored only when the active
      !> mechanism registered an analytical Jacobian.
      call fini%get(section_name= trim(codename)//'-ODE', option_name='analytical-jacobian', &
                    val=use_analytical_jacobian, error=error)
      if (error/=0) use_analytical_jacobian = .false.
      if (use_analytical_jacobian) then
        if (.not. associated(chemistry_jacobian)) then
          write(*,'(A)') ' [JAC] analytical-jacobian=true but mechanism "'//trim(self%mechanism_name)//&
                         '" has no analytical Jacobian. Falling back to FD.'
          use_analytical_jacobian = .false.
        else
          write(*,'(A)') ' [JAC] analytical Jacobian enabled (IJAC=1)'
          call validate_analytical_jacobian()
        end if
      end if

    else
      !> FROZEN run (no mechanism read -> error_mechanism/=0): chemistry stays
      !> off, but the device convective/diffusive/thermo kernels still need the
      !> read-only thermo tables resident. flint_acc_upload_tables lives in the
      !> reactive branch above, so without this the device sees ns=0 / Tmin=Tmax=0
      !> -> rotot=sum(rhoi(1:0))=0 -> Tl=p/rho/R=Inf -> int(Tl) OOB-indexes h_tab
      !> in the convective kernel (the frozen multi-block ILLEGAL_ADDRESS). Push
      !> thermo only; the rate tables kf_tab/kb_tab are unallocated here. No-op on
      !> the CPU build (the !$acc directives are the whole body).
      !> The frozen path is the FIRST device data op for FLINT thermo (no chemistry
      !> kernels ran to materialize the declare-create thermo SCALARS, unlike the
      !> reactive branch above). nvfortran honors `declare create` for module ARRAYS
      !> at startup but NOT scalars, so flint_acc_upload_thermo's
      !> `update device(ns,Tmin,Tmax)` faults "data not found on device". Create them
      !> here first; frozen-only -> no double-present clash on the reactive path.
      !> Same fix for the limiter scalars (Assign_Limiter runs after this point and
      !> its `update device(SD_space_order,lim_dev)` would otherwise fault on frozen;
      !> the values set there are pushed by that later update -- here we only create).
      !$acc enter data copyin(ns, Tmin, Tmax, SD_space_order, lim_dev, &
      !$acc&                  eddy_dev, ransdiff_dev, stress_dev, sc_dev, sct_dev, prt_dev, prl_dev, RANS, VISC, &
      !$acc&                  nsc, nu, nv, nw, np, nt, nrans, nprim, gc, rs_dev, mcdiff_dev, mco_dev, sd_dev)
      !> DIFF_MC: same for the scalar of the binary-diffusion table (dij_tab itself is an allocatable).
      if (allocated(dij_tab)) then
        !$acc enter data copyin(dij_pref)
      end if
      call flint_acc_upload_thermo()

    endif

    call set_strangsplitting(self%reactive)
        
    nochem_str = ''
    args= ''
    call fini%get(section_name='MOSE-nochem', option_name='blocks', val=nochem_str, error=error)
    if (allocated(self%no_chem_list)) deallocate(self%no_chem_list)

    nochem_str = trim(adjustl(nochem_str))
    if (error==0 .and. len_trim(nochem_str) > 0) then
      call parse(nochem_str,' ',args)
      nochem_blocks = count(args /= '')
      allocate(self%no_chem_list(nochem_blocks))
      do i = 1, nochem_blocks
        read(args(i),'(I3)',iostat=ios) self%no_chem_list(i)
        if (ios /= 0) then
          write(*,*) ' Warning: MOSE-nochem/blocks parsing failed. String: ', nochem_str
          deallocate(self%no_chem_list)
          allocate(self%no_chem_list(0))
        end if
      end do
    else
      allocate(self%no_chem_list(0))  ! lista vuota
    end if
    
  end subroutine read_ini


  ! subroutine read_chemistry_file_legacy( file1, file2 )
  !   use IR_Precision
  !   use MOSE_Lib_ChemMech, only: kf_tab, kb_tab, ni1_tab, ni2_tab, epsch_tab
  !   use MOSE_Global_m, only: nsc, nrc
  !   implicit none
  !   character(len=*), intent(in), optional :: file1, file2
  !   ! Local
  !   integer :: idum, unitfile
  !   integer :: i, j
  !   integer :: Tmax_, ios
  !   character(len=16) :: try
  !   logical :: is_present

  !   if (.not.present(file1) .and. .not.present(file2)) then
  !     do i = 1, 99
  !       try = 'MISCELA/'//trim(str(.true.,i))
  !       inquire(file=trim(try)//'sp/wm.dat',exist=is_present)
  !       if (is_present) then
  !         obj_chem%ratefile = trim(try)//'sp/rate.dat'
  !         obj_chem%stoichfile = trim(try)//'sp/stoich.dat'
  !         exit
  !       endif
  !     enddo
  !   endif

  !   if (present(file1)) then
  !     open(newunit=unitfile,file=file1,form='formatted',iostat=ios, status='old', action='read')
  !   else
  !     open(newunit=unitfile,file=obj_chem%ratefile,form='formatted',iostat=ios, status='old', action='read')
  !   endif
  !   if (ios/=0) return
  !   read(unitfile,*) nrc, Tmax_
  !   allocate(kf_tab(1:Tmax_,1:nrc))
  !   allocate(kb_tab(1:Tmax_,1:nrc))
  !   do j = 1, nrc; do i = 1, Tmax_
  !       read(unitfile,*) idum,idum,kf_tab(i,j),kb_tab(i,j)
  !   enddo; enddo 
  !   close(unitfile)

  !   if (present(file2)) then
  !     open(newunit=unitfile,file=file2,form='formatted',iostat=ios)
  !   else
  !     open(newunit=unitfile,file=obj_chem%stoichfile,form='formatted',iostat=ios)
  !   endif
  !   read(unitfile,*) nrc, nsc
  !   allocate(ni1_tab(1:nsc+1,1:nrc))
  !   allocate(ni2_tab(1:nsc+1,1:nrc))
  !   allocate(epsch_tab(1:nsc+1,1:nrc))
  !   do j = 1, nrc; do i = 1, nsc+1
  !       read(unitfile,*)idum,idum,ni1_tab(i,j),ni2_tab(i,j),epsch_tab(i,j)
  !   enddo; enddo 
  !   close(unitfile)

  ! end subroutine read_chemistry_file_legacy


end module MOSE_IO_chemistry