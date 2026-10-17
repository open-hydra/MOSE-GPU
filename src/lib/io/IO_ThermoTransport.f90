module MOSE_IO_Table
  use MOSE_Parameters_m
  use iso_fortran_env, only: R8_tab => real64

  implicit none

  type :: tab_type
  contains
      procedure, nopass :: ini_setup => read_ini
  end type tab_type

  type(tab_type) :: obj_tab
  private :: read_ini

  !> DIFF_MC (2026-09-26): laminar species-diffusion closure, key [MOSE-Physics] schmidt (mirror of MOSE_open).
  !>   key absent  -> constant laminar Schmidt from [MOSE-RANS] Sc (default 0.7): the production path, unchanged;
  !>   schmidt > 0 -> constant laminar Schmidt = schmidt (replaces [MOSE-RANS] Sc, as in the reference);
  !>   schmidt <= 0 -> mixture-averaged multicomponent D_k (FLINT co_DS_expr) from INPUT/<prefix>diffusion.dat.
  !> vnn-diff ([MOSE-parameters], default 0.8) is read only with schmidt <= 0: dt_D = vnn-diff * dl^2 / max_k D_k
  !> (D_k at the cell + mu_t/(rho Sct)); vnn-diff = 0 -> no such term.
  logical,       public :: mcdiff        = .false.
  logical,       public :: schmidt_given = .false.
  real(R8_tab),  public :: schmidt_val   = 0.7d0
  real(R8_tab),  public :: prandtl_val   = 0.0d0   !> [MOSE-Physics] prandtl: > 0 fixes the laminar conductivity kl = mu cp / Pr (MOSE_open semantics), 0 = transport table (Wilke)
  real(R8_tab),  public :: vnn_diff      = 0.8d0
  private :: read_diffusion_table

contains

  subroutine read_ini(fini)
    use Finer, only: file_ini
    use MOSE_Global_m, only: nsc, MOSE_phase_prefix
    use FLINT_Load_ThermoTransport
    use FLINT_Lib_Thermodynamic, only: FLINT_phase_prefix

    implicit none
    type(file_ini)  :: fini
    ! Local
    character(llen) :: folder
    integer :: error

    call fini%get(section_name=trim(codename)//'-IO', option_name='idealgas-folder', val=folder, error=error)
    if (error/=0) folder = 'INPUT/'

    FLINT_phase_prefix = MOSE_phase_prefix
    error = read_idealgas_thermo( folder )
    if (error /= 0) then
      write(*,*) 'Error reading ideal gas thermo data, using legacy tables'
      call read_tables_legacy()
    end if

    error = read_idealgas_transport( folder )
    if (error /= 0) then
      write(*,*) '[WARNING] Transport data not found'
    end if

    ! DIFF_MC: laminar species-diffusion closure. Key absent -> nothing below runs (production path unchanged).
    call fini%get(section_name=trim(codename)//'-Physics', option_name='schmidt', val=schmidt_val, error=error)
    schmidt_given = ( error == 0 )
    call fini%get(section_name=trim(codename)//'-Physics', option_name='prandtl', val=prandtl_val, error=error)
    if ( error /= 0 ) prandtl_val = 0.0d0
    if ( prandtl_val > 0.0d0 ) write(*,'(A,F8.4,A)') ' [Physics] prandtl = ', prandtl_val, ': laminar conductivity kl = mu cp / Pr (as MOSE_open)'
    if ( schmidt_given .and. schmidt_val <= 0.0d0 ) then
      mcdiff = .true.
      call fini%get(section_name=trim(codename)//'-parameters', option_name='vnn-diff', val=vnn_diff, error=error)
      if ( error /= 0 ) vnn_diff = 0.8d0
      if ( vnn_diff < 0.0d0 ) error stop '[DIFF_MC] vnn-diff must be >= 0 (0 = no species-diffusion time-step limit)'
      call read_diffusion_table ( folder )
    end if

    if (ns /= nsc) error stop 'this binary was built for a different species count than the mechanism (one binary per mechanism)'

  end subroutine read_ini

  !> DIFF_MC D3: read the binary-diffusion table (FLINT read_idealgas_diffusion) and STOP with a message on every failure
  !> (never run co_DS_expr on an unallocated or mismatched table). Checks, in order: exactly one table file exists (.dat or
  !> .szplt: with both, FLINT would read the two into one object); the reader succeeds (code 2 = unreadable, 4 = not ns(ns-1)/2 zones, 3 = Tmax differs from the thermo table); the table
  !> starts at or below the thermo Tmin (co_DS_expr clamps T to [Tmin, Tmax] of the thermo and indexes dij_tab there);
  !> and, for the ASCII table, the zone names are the pairs of the phase.txt species in upper-triangular order (the reader
  !> itself does not look at the names: a table for another species order would be read silently).
  subroutine read_diffusion_table ( folder )
    use FLINT_Load_ThermoTransport, only: read_idealgas_diffusion
    use FLINT_Lib_Thermodynamic,    only: ns, Tmin, Tmax, dij_tab, dij_pref, ndij, species_names, FLINT_phase_prefix
    implicit none
    character(*), intent(in) :: folder
    integer :: ios, u, i, j, p, q
    character(512) :: fdat, fszplt, line
    character(128) :: zname, expect
    logical :: exdat, exszplt

    fdat   = trim(folder)//'/'//trim(FLINT_phase_prefix)//'diffusion.dat'
    fszplt = trim(folder)//'/'//trim(FLINT_phase_prefix)//'diffusion.szplt'
    inquire(file=trim(fdat), exist=exdat)
    inquire(file=trim(fszplt), exist=exszplt)
    if ( .not. exdat .and. .not. exszplt ) then
      write(*,'(A)') ' [DIFF_MC] ERROR: schmidt <= 0 (multicomponent diffusion) needs '//trim(fdat)//' (or .szplt): not found'
      error stop '[DIFF_MC] binary-diffusion table missing'
    end if
    if ( exdat .and. exszplt ) then
      write(*,'(A)') ' [DIFF_MC] ERROR: both '//trim(fdat)//' and '//trim(fszplt)//' exist: keep one'
      error stop '[DIFF_MC] binary-diffusion table ambiguous (.dat and .szplt)'
    end if
    ios = read_idealgas_diffusion( folder )
    select case (ios)
    case (0)
    case (2)
      write(*,'(A)') ' [DIFF_MC] ERROR: the binary-diffusion table exists but could not be read: '//trim(fdat)
      error stop '[DIFF_MC] binary-diffusion table unreadable (code 2)'
    case (3)
      write(*,'(A,I0,A)') ' [DIFF_MC] ERROR: the binary-diffusion table does not end at the thermo Tmax = ', Tmax, &
                          ' K (regenerate it on the same T grid as thermo.dat)'
      error stop '[DIFF_MC] binary-diffusion table Tmax mismatch (code 3)'
    case (4)
      write(*,'(A,I0,A)') ' [DIFF_MC] ERROR: the binary-diffusion table does not have ns(ns-1)/2 = ', ns*(ns-1)/2, &
                          ' zones (one per species pair)'
      error stop '[DIFF_MC] binary-diffusion table zone count (code 4)'
    case default
      write(*,'(A,I0)') ' [DIFF_MC] ERROR: read_idealgas_diffusion failed with code ', ios
      error stop '[DIFF_MC] binary-diffusion table read failure'
    end select
    if ( .not. allocated(dij_tab) ) error stop '[DIFF_MC] binary-diffusion table not allocated after a successful read'
    if ( lbound(dij_tab,1) > Tmin .or. ubound(dij_tab,1) /= Tmax ) then
      write(*,'(A,I0,A,I0,A,I0,A,I0,A)') ' [DIFF_MC] ERROR: the binary-diffusion table covers T = ', lbound(dij_tab,1), '..', &
            ubound(dij_tab,1), ' K, the thermo table ', Tmin, '..', Tmax, ' K'
      error stop '[DIFF_MC] binary-diffusion table does not cover the thermo T range (Tmin)'
    end if
    if ( exdat ) then
      open(newunit=u, file=trim(fdat), status='old', action='read', iostat=ios)
      if ( ios /= 0 ) error stop '[DIFF_MC] cannot reopen the binary-diffusion table for the zone-name check'
      p = 0
      do
        read(u,'(A)',iostat=ios) line
        if ( ios /= 0 ) exit
        line = adjustl(line)
        if ( line(1:4) /= 'ZONE' ) cycle
        q = index(line, '"')
        zname = line(q+1:)
        q = index(zname, '"')
        if ( q > 0 ) zname = zname(1:q-1)
        p = p + 1
        expect = '(none: more zones than species pairs)'
        q = 0
        outer: do i = 1, ns
          do j = i+1, ns
            q = q + 1
            if ( q == p ) then
              expect = trim(species_names(i))//'-'//trim(species_names(j))
              exit outer
            end if
          end do
        end do outer
        if ( p > ns*(ns-1)/2 .or. trim(zname) /= trim(expect) ) then
          write(*,'(A,I0,A)') ' [DIFF_MC] ERROR: zone ', p, ' of '//trim(fdat)//' is "'//trim(zname)//'", expected "'// &
                              trim(expect)//'" (pairs of phase.txt in upper-triangular order)'
          error stop '[DIFF_MC] binary-diffusion table species order'
        end if
      end do
      close(u)
    else
      write(*,'(A)') ' [DIFF_MC] note: .szplt table, zone-name check skipped (read with the .dat for the check)'
    end if
    write(*,'(A)') ' [DIFF_MC] schmidt <= 0: mixture-averaged multicomponent diffusion (FLINT co_DS_expr), table '// &
                   trim(merge(fdat, fszplt, exdat))
    write(*,'(A,I0,A,I0,A,I0,A,ES12.5,A,F6.3)') ' [DIFF_MC]   pairs ', ndij, ', T ', lbound(dij_tab,1), '..', ubound(dij_tab,1), &
                   ' K, Pref ', dij_pref, ' Pa, vnn-diff ', vnn_diff
  end subroutine read_diffusion_table

  subroutine read_tables_legacy()
    use FLINT_Lib_Thermodynamic
    implicit none
    integer :: i, j
    integer :: ios, unitfile
    integer :: dummy_check_line(2)
    logical :: present_limits, is_present
    character(80) :: try
    character(200):: wmfile
    character(200):: tabmsfile
    character(2)  :: n

    do i = 1, 99
      write(n,'(I2)') i
      try = 'MISCELA/'//trim(adjustl(n))
      inquire(file=trim(try)//'sp/wm.dat',exist=is_present)
      if (is_present) then
        wmfile = trim(try)//'sp/wm.dat'
        tabmsfile = trim(try)//'sp/tabellams.dat'
        exit
      endif
    enddo

    ! Open file for molecular weights (usually wm.dat)
    open(newunit=unitfile, file=wmfile, form='formatted')
    ios = 0
    do while(ios==0)
      read(unitfile,*,iostat=ios)
      ns = ns+1
    end do
    ns = ns-1 ! Setting the number of species from the lines of wm.dat
    allocate(wm_tab(ns))
    rewind(unitfile)
    do i = 1, ns
      read(unitfile,*) wm_tab(i)
    end do
    close(unitfile)
    ! write(*,*) '  Molecular weights file: ', wmfile
    ! write(*,*) '  Number of species  -->  ', ns

    ! Gas constants for the species
    Ri_tab = Runiv/wm_tab

    ! Open file for thermodynamic and trasport properties tabmsfile (usually tabellams.dat)
    open(newunit=unitfile, file=tabmsfile, form='formatted')
    read(unitfile,'(2I8)',iostat=ios)  dummy_check_line
    if (sum(dummy_check_line)/=2) then
      present_limits = .true.
    else
      present_limits = .false.
      rewind(unitfile)
    endif
    ! Counting the number of lines for Tmax
    ios = 0
    do while(ios==0)
      read(unitfile,*,iostat=ios)
      Tmax = Tmax+1
    end do
    Tmax = (Tmax-1)/ns ! Setting the maximum temperature
    ! write(*,*) '  Temperature limit for the simulation  -->  ', Tmax
    rewind(unitfile)

    Tmin = 0
    allocate(h_tab(Tmin:Tmax, 1:ns))
    allocate(s_tab, cp_tab, dcpi_tab, mi_tab, k_tab, mold=h_tab)

    if (present_limits) read(unitfile,*)
    do j = 1, ns
      do i = 1, Tmax
        read(unitfile,*) dummy_check_line, h_tab(i,j), s_tab(i,j), cp_tab(i,j), dcpi_tab(i,j), mi_tab(i,j), k_tab(i,j)
      enddo
      s_tab(0,j) = s_tab(1,j)
      cp_tab(0,j) = cp_tab(1,j)
      dcpi_tab(0,j) = dcpi_tab(1,j)
      h_tab(0,j) = cp_tab(1,j)*0.d0
      s_tab(0,j) = s_tab(1,j)
      mi_tab(0,j) = mi_tab(1,j)
      k_tab(0,j) = k_tab(1,j)
    enddo
    close(unitfile)
    write(*,*) '  Legacy table read'

  end subroutine read_tables_legacy


end module MOSE_IO_Table
