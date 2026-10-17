module FLINT_Load_Chemistry
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none

contains

  function read_chemistry( folder, mech_name ) result(ios)
    use Lib_ORION_Data
    use Lib_Tecplot
    use FLINT_Lib_Chemistry_data
    use FLINT_Lib_Thermodynamic, only: ns, FLINT_phase_prefix
    implicit none
    character(len=*), intent(in), optional :: folder
    character(len=*), intent(out), optional :: mech_name
    integer :: idum, unitfile
    type(ORION_data)  :: orion
    integer :: i, j, ios, j0, j1, j2
    integer :: Ti1, Ti2, dummy1, dummy23
    character(len=32):: chardum

    nrc_arrh = 0
    nrc_troe = 0
    nrc_lindemann = 0

    !! Info
    if (present(folder)) then
      open(newunit=unitfile,file=trim(folder)//'/'//trim(FLINT_phase_prefix)//'chemistry-info.txt',form='formatted',status='old',action='read',iostat=ios)
    else
      open(newunit=unitfile,file='INPUT/'//trim(FLINT_phase_prefix)//'chemistry-info.txt',form='formatted',status='old',action='read',iostat=ios)
    endif
    if (ios/=0) then
      write(*,*) '[WARNING] Chemistry data not found'
      return
    endif
    ! Read mechanism name
    read(unitfile,*) mech_name
    read(unitfile,*)
    read(unitfile,'(A17,I4)') chardum, nrc
    do j = 1, 4; read(unitfile,*); enddo
    ! Read reaction type
    allocate(rxn_type(1:nrc))
    do j = 1, nrc
      read(unitfile,*) idum, chardum
      if (index(trim(chardum),'Troe')>0) then
        nrc_troe = nrc_troe + 1
        rxn_type(j) = 1
      elseif (index(trim(chardum),'Lindemann')>0) then
        nrc_lindemann = nrc_lindemann + 1
        rxn_type(j) = 2
      else
        nrc_arrh = nrc_arrh + 1
        rxn_type(j) = 0
      endif
    enddo
    ! Read info for general loop
    allocate(ni1_arrh_tab(1:ns+1,1:nrc_arrh))
    allocate(ni2_arrh_tab(1:ns+1,1:nrc_arrh))
    allocate(epsch_arrh_tab(1:ns+1,1:nrc_arrh))
    if (nrc_troe>0) then
      allocate(ni1_troe_tab(1:ns+1,1:nrc_troe))
      allocate(ni2_troe_tab(1:ns+1,1:nrc_troe))
      allocate(epsch_troe_tab(1:ns+1,1:nrc_troe))
    endif
    if (nrc_lindemann>0) then
      allocate(ni1_lind_tab(1:ns+1,1:nrc_lindemann))
      allocate(ni2_lind_tab(1:ns+1,1:nrc_lindemann))
      allocate(epsch_lind_tab(1:ns+1,1:nrc_lindemann))
    endif
    ! Read reaction info
    read(unitfile,*)
    read(unitfile,*)
    j0 = 0; j1 = 0; j2 = 0
    do j = 1, nrc
      if (rxn_type(j)==0) then
        j0 = j0+1
        do i = 1, ns+1
          read(unitfile,*)idum,chardum,ni1_arrh_tab(i,j0),ni2_arrh_tab(i,j0),epsch_arrh_tab(i,j0)
        enddo
      elseif (rxn_type(j)==1) then
        j1 = j1+1
        do i = 1, ns+1
          read(unitfile,*)idum,chardum,ni1_troe_tab(i,j1),ni2_troe_tab(i,j1),epsch_troe_tab(i,j1)
        enddo
      elseif (rxn_type(j)==2) then
        j2 = j2+1
        do i = 1, ns+1
          read(unitfile,*)idum,chardum,ni1_lind_tab(i,j2),ni2_lind_tab(i,j2),epsch_lind_tab(i,j2)
        enddo
      endif
    enddo 
    close(unitfile)

    !! Rate Arrhenius
    if (present(folder)) then
      ios = tec_read_points_multivars(orion,2,trim(folder)//'/'//trim(FLINT_phase_prefix)//'chemistry-Arrhenius.dat')
      if (ios/=0) ios = tec_read_structured_multiblock(orion=orion, filename=trim(folder)//'/'//trim(FLINT_phase_prefix)//'chemistry-Arrhenius.szplt')
    else
      ios = tec_read_points_multivars(orion,2,trim('INPUT/')//trim(FLINT_phase_prefix)//'chemistry-Arrhenius.dat')
      if (ios/=0) ios = tec_read_structured_multiblock(orion=orion, filename=trim('INPUT/')//trim(FLINT_phase_prefix)//'chemistry-Arrhenius.szplt')
    endif
    if (ios/=0) then
      write(*,*) '[ERROR] chemistry-Arrhenius file not found'
      stop
    endif
    dummy1  = lbound(orion%block(1)%mesh, dim=2)
    dummy23 = lbound(orion%block(1)%mesh, dim=3)
    Ti1 = nint(orion%block(1)%mesh(1,dummy1,dummy23,dummy23))
    Ti2 = Ti1 + ubound(orion%block(1)%mesh, dim=2) - dummy1
    allocate(kf_tab(Ti1:Ti2, 1:nrc_arrh))
    allocate(kb_tab(Ti1:Ti2, 1:nrc_arrh))
    dummy23 = lbound(orion%block(1)%vars, dim=3)
    do i = 1, nrc_arrh
      kf_tab(Ti1:Ti2,i) = orion%block(i)%vars(1,:,dummy23,dummy23)
      kb_tab(Ti1:Ti2,i) = orion%block(i)%vars(2,:,dummy23,dummy23)
    enddo

    !! Rate Troe
    if (nrc_troe/=0) then
      deallocate(orion%block)
      if (present(folder)) then
        ios = tec_read_points_multivars(orion,4,trim(folder)//'/'//trim(FLINT_phase_prefix)//'chemistry-Troe.dat')
        if (ios/=0) ios = tec_read_structured_multiblock(orion=orion, filename=trim(folder)//'/'//trim(FLINT_phase_prefix)//'chemistry-Troe.szplt')
      else
        ios = tec_read_points_multivars(orion,4,trim('INPUT/')//trim(FLINT_phase_prefix)//'chemistry-Troe.dat')
        if (ios/=0) ios = tec_read_structured_multiblock(orion=orion, filename=trim('INPUT/')//trim(FLINT_phase_prefix)//'chemistry-Troe.szplt')
      endif
      if (ios/=0) then
        write(*,*) '[ERROR] chemistry-Troe file not found'
        stop
      endif
      allocate(Fcent_tab(Ti1:Ti2, 1:nrc_troe))
      allocate(k0_troe_tab, kinf_troe_tab, kc_troe_tab, mold=Fcent_tab)
      do i = 1, nrc_troe
        kinf_troe_tab(Ti1:Ti2,i) = orion%block(i)%vars(1,:,dummy23,dummy23)
        k0_troe_tab(Ti1:Ti2,i)   = orion%block(i)%vars(2,:,dummy23,dummy23)
        kc_troe_tab(Ti1:Ti2,i)   = orion%block(i)%vars(3,:,dummy23,dummy23)
        Fcent_tab(Ti1:Ti2,i) = orion%block(i)%vars(4,:,dummy23,dummy23)
      enddo
    endif
  
    !! Rate Lindemann
    if (nrc_lindemann/=0) then
      deallocate(orion%block)
      if (present(folder)) then
        ios = tec_read_points_multivars(orion,3,trim(folder)//'/'//trim(FLINT_phase_prefix)//'chemistry-Lindemann.dat')
        if (ios/=0) ios = tec_read_structured_multiblock(orion=orion, filename=trim(folder)//'/'//trim(FLINT_phase_prefix)//'chemistry-Lindemann.szplt')
      else
        ios = tec_read_points_multivars(orion,3,trim('INPUT/')//trim(FLINT_phase_prefix)//'chemistry-Lindemann.dat')
        if (ios/=0) ios = tec_read_structured_multiblock(orion=orion, filename=trim('INPUT/')//trim(FLINT_phase_prefix)//'chemistry-Lindemann.szplt')
      endif
      if (ios/=0) then
        write(*,*) '[ERROR] chemistry-Lindemann file not found'
        stop
      endif
      allocate(kinf_lind_tab(Ti1:Ti2, 1:nrc_lindemann))
      allocate(k0_lind_tab, kc_lind_tab, mold=kinf_lind_tab)
      do i = 1, nrc_lindemann
        kinf_lind_tab(Ti1:Ti2,i) = orion%block(i)%vars(1,:,dummy23,dummy23)
        k0_lind_tab(Ti1:Ti2,i)   = orion%block(i)%vars(2,:,dummy23,dummy23)
        kc_lind_tab(Ti1:Ti2,i)   = orion%block(i)%vars(3,:,dummy23,dummy23)
      enddo
    endif

  end function read_chemistry


end module FLINT_Load_Chemistry