module FLINT_Load_ThermoTransport
  use FLINT_Lib_Thermodynamic
  implicit none

contains

  function read_idealgas_thermo(folder) result(ios)
    use strings, only: parse
    use Lib_Tecplot
    use Lib_ORION_data
    implicit none
    character(len=*), intent(in), optional :: folder
    ! Local
    integer           :: ios, i, unitfile, start, dummy1, dummy23
    character(256)    :: wholestring, args(2)
    character(512)    :: wmfile, thermofile(2)
    type(ORION_data)  :: orion

    if (present(folder)) then
      wmfile = trim(folder)//'/'//trim(FLINT_phase_prefix)//'phase.txt'
      thermofile(1) = trim(folder)//'/'//trim(FLINT_phase_prefix)//'thermo.dat'
      thermofile(2) = trim(folder)//'/'//trim(FLINT_phase_prefix)//'thermo.szplt'
    else
      wmfile = 'INPUT/'//trim(FLINT_phase_prefix)//'phase.txt'
      thermofile(1) = 'INPUT/'//trim(FLINT_phase_prefix)//'thermo.dat'
      thermofile(2) = 'INPUT/'//trim(FLINT_phase_prefix)//'thermo.szplt'
    endif

    ! File 1: phase
    open(newunit=unitFile,file=trim(wmfile),status='old',iostat=ios)
    if (ios/=0) then
      write(*,*) '[ERROR] phase.txt type file not found'
      return
    endif
    ios = 0; ns = -1
    read(unitFile,*)
    do while(ios==0)
      read(unitFile,'(A)',iostat=ios) wholestring
      ns = ns + 1
    enddo
    allocate(wm_tab(1:ns))
    allocate(Ri_tab(1:ns))
    allocate(species_names(1:ns))
    rewind(unitFile)
    read(unitFile,*)
    do i = 1, ns
      read(unitFile,'(A)') wholestring
      call parse(wholestring,' ',args)
      species_names(i) = trim(args(1))
      read(args(2),*) wm_tab(i)
    end do
    close(unitFile)

    Ri_tab = Runiv/wm_tab

    ! Precompute species-pair constants used by co_fiij (Wilke mixing rule)
    allocate(Mi_Mj_pow_m025(1:ns, 1:ns))
    allocate(inv_sqrt8_1p  (1:ns, 1:ns))
    block
      integer :: si, sj
      real(kind=8) :: ratio
      do sj = 1, ns
        do si = 1, ns
          ratio = wm_tab(si) / wm_tab(sj)
          Mi_Mj_pow_m025(si,sj) = ratio**(-0.25d0)
          inv_sqrt8_1p  (si,sj) = 1.d0 / dsqrt(8.d0 * (1.d0 + ratio))
        end do
      end do
    end block

    ! File 2: thermo
    ios = tec_read_points_multivars(orion,4,trim(thermofile(1)))
    if (ios/=0) then
      ios = tec_read_structured_multiblock(orion=orion, filename=trim(thermofile(2)))
    endif
    if (ios/=0) then
      write(*,*) '[ERROR] Reading thermo file'
      return
    endif
    dummy1  = lbound(orion%block(1)%mesh, dim=2)
    dummy23 = lbound(orion%block(1)%mesh, dim=3)
    Tmin = nint(orion%block(1)%mesh(1,dummy1,dummy23,dummy23))
    Tmax = Tmin + ubound(orion%block(1)%mesh, dim=2) - dummy1
    start = Tmin
    if (Tmin==1) Tmin = 0
    allocate(h_tab(Tmin:Tmax, 1:ns))
    allocate(s_tab, cp_tab, dcpi_tab, mold=h_tab)
    dummy23 = lbound(orion%block(1)%vars, dim=3)
    do i = 1, ns
      cp_tab(start:Tmax,i) = orion%block(i)%vars(1,:,dummy23,dummy23)
      h_tab(start:Tmax,i) = orion%block(i)%vars(2,:,dummy23,dummy23)
      s_tab(start:Tmax,i) = orion%block(i)%vars(3,:,dummy23,dummy23)
      dcpi_tab(start:Tmax,i) = orion%block(i)%vars(4,:,dummy23,dummy23)
    enddo
    if (Tmin == 0) then
      cp_tab(0,:) = cp_tab(1,:)
      h_tab(0,:) = cp_tab(1,:)*0.d0
      s_tab(0,:) = s_tab(1,:)
      dcpi_tab(0,:) = dcpi_tab(1,:)
    endif

    ios = read_idealgas_composition(folder)

  end function read_idealgas_thermo


  function read_idealgas_composition(folder) result(ios)
    use strings, only: parse
    implicit none
    character(len=*), intent(in), optional :: folder
    ! Local
    integer           :: ios, i, unitfile
    character(512)    :: file

    if (present(folder)) then
      file = trim(folder)//'/'//trim(FLINT_phase_prefix)//'composition.txt'
    else
      file = 'INPUT/'//trim(FLINT_phase_prefix)//'composition.txt'
    endif

    ! File 1: phase
    open(newunit=unitFile,file=trim(file),status='old',iostat=ios)
    if (ios/=0) then
      write(*,*) '[WARNING] Composition data not found'
      ios = 0
      return
    endif
    read(unitfile,*,iostat=ios) ne
    allocate(elements_names(1:ne))
    allocate(species_composition(1:ne,1:ns))
    do i = 1, ne
      read(unitFile,'(A)',iostat=ios) elements_names(i)
    enddo
    read(unitFile,*,iostat=ios) ! blank
    read(unitFile,*,iostat=ios) ! composition word
    do i = 1, ns
      read(unitFile,*,iostat=ios) species_composition(:,i)
    end do
    close(unitFile)

  end function read_idealgas_composition


  function read_idealgas_transport(folder) result(ios)
    use strings, only: parse
    use Lib_Tecplot
    use Lib_ORION_data
    implicit none
    character(len=*), intent(in), optional :: folder
    ! Local
    integer           :: ios, i, start, dummy1, dummy23
    logical           :: exists
    character(512)    :: transfile(2)
    type(ORION_data)  :: orion

    ios = 1

    if (present(folder)) then
      transfile(1) = trim(folder)//'/'//trim(FLINT_phase_prefix)//'transport.dat'
      transfile(2) = trim(folder)//'/'//trim(FLINT_phase_prefix)//'transport.szplt'
    else
      transfile(1) = 'INPUT/'//trim(FLINT_phase_prefix)//'transport.dat'
      transfile(2) = 'INPUT/'//trim(FLINT_phase_prefix)//'transport.szplt'
    endif

    inquire(file=trim(transfile(1)),exist=exists)
    if (exists) ios = tec_read_points_multivars(orion,2,trim(transfile(1)))
    inquire(file=trim(transfile(2)),exist=exists)
    if (exists) ios = tec_read_structured_multiblock(orion=orion, filename=trim(transfile(2)))
    if (ios/=0) return
    dummy1  = lbound(orion%block(1)%mesh, dim=2)
    dummy23 = lbound(orion%block(1)%mesh, dim=3)
    Tmin = nint(orion%block(1)%mesh(1,dummy1,dummy23,dummy23))
    Tmax = Tmin + ubound(orion%block(1)%mesh, dim=2) - dummy1
    start = Tmin
    if (Tmin==1) Tmin = 0
    allocate(mi_tab(Tmin:Tmax, 1:ns))
    allocate(k_tab(Tmin:Tmax, 1:ns))
    dummy23 = lbound(orion%block(1)%vars, dim=3)
    do i = 1, ns
      mi_tab(start:Tmax,i) = orion%block(i)%vars(1,:,dummy23,dummy23)
      k_tab(start:Tmax,i) = orion%block(i)%vars(2,:,dummy23,dummy23)
    enddo
    if (Tmin==0) then
      mi_tab(0,:) = mi_tab(1,:)
      k_tab(0,:) = k_tab(1,:)
    endif

  end function read_idealgas_transport


  ! Reads binary diffusion coefficients D_ij(T) from <prefix>diffusion.dat.
  ! One ZONE per unique unordered species pair, upper-triangular order
  ! (1,2),(1,3),..,(1,ns),(2,3),..  Each zone holds a single variable (Dij) vs T,
  ! on the same temperature grid as the thermo/transport tables.
  ! Populates dij_tab(Tmin:Tmax, 1:ndij) with ndij = ns*(ns-1)/2.
  function read_idealgas_diffusion(folder) result(ios)
    use Lib_Tecplot
    use Lib_ORION_data
    implicit none
    character(len=*), intent(in), optional :: folder
    ! Local
    integer           :: ios, i, start, dummy1, dummy23
    integer           :: Tmin_dummy, Tmax_dummy
    integer           :: u, ios2, ip
    logical           :: exists, attempted
    character(512)    :: difffile(2), titleline
    type(ORION_data)  :: orion

    ios = 1
    attempted = .false.
    dij_pref = 101325.d0   ! default reference pressure (1 atm) if TITLE has no Pref tag

    ! Single-species mixtures have no pairs: nothing to load.
    ndij = ns*(ns-1)/2
    if (ndij < 1) return

    if (present(folder)) then
      difffile(1) = trim(folder)//'/'//trim(FLINT_phase_prefix)//'diffusion.dat'
      difffile(2) = trim(folder)//'/'//trim(FLINT_phase_prefix)//'diffusion.szplt'
    else
      difffile(1) = 'INPUT/'//trim(FLINT_phase_prefix)//'diffusion.dat'
      difffile(2) = 'INPUT/'//trim(FLINT_phase_prefix)//'diffusion.szplt'
    endif

    inquire(file=trim(difffile(1)),exist=exists)
    if (exists) then
      attempted = .true.
      ! Parse the reference pressure from the TITLE line: "... (Pref=<value> Pa)".
      open(newunit=u, file=trim(difffile(1)), status='old', iostat=ios2)
      if (ios2==0) then
        read(u,'(A)',iostat=ios2) titleline
        if (ios2==0) then
          ip = index(titleline,'Pref=')
          if (ip>0) read(titleline(ip+5:),*,iostat=ios2) dij_pref
        endif
        close(u)
      endif
      ios = tec_read_points_multivars(orion,1,trim(difffile(1)))
    endif
    inquire(file=trim(difffile(2)),exist=exists)
    if (exists) then
      attempted = .true.
      ios = tec_read_structured_multiblock(orion=orion, filename=trim(difffile(2)))
    endif
    if (ios/=0) then
      if (attempted) ios = 2
      return
    endif

    if (size(orion%block) /= ndij) then
      ios = 4
      return
    endif

    dummy1  = lbound(orion%block(1)%mesh, dim=2)
    dummy23 = lbound(orion%block(1)%mesh, dim=3)
    Tmin_dummy = nint(orion%block(1)%mesh(1,dummy1,dummy23,dummy23))
    Tmax_dummy = Tmin_dummy + ubound(orion%block(1)%mesh, dim=2) - dummy1

    if (Tmax_dummy /= Tmax) then
      ios = 3
      return
    endif

    start = Tmin_dummy
    if (Tmin_dummy==1) Tmin_dummy = 0
    allocate(dij_tab(Tmin_dummy:Tmax, 1:ndij))
    dummy23 = lbound(orion%block(1)%vars, dim=3)
    do i = 1, ndij
      dij_tab(start:Tmax,i) = orion%block(i)%vars(1,:,dummy23,dummy23)
    enddo
    if (Tmin_dummy==0) dij_tab(0,:) = dij_tab(1,:)

  end function read_idealgas_diffusion

end module FLINT_Load_ThermoTransport
