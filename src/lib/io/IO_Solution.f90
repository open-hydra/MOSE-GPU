module MOSE_IO_Solution
  use MOSE_Parameters_m

  implicit none

  type :: IO_type
    integer             :: sol_diter, bck_diter, res_diter !> iteration interval to save the solution
    real(8)             :: sol_dtime, bck_dtime            !> time interval to save the solution
    integer             :: sol_field_every = 1             !> sol-field-every: 3D field at every N-th time-based dump (wall file at every one)
    integer             :: n_tdump = 0                     !> time-based dumps taken so far
    logical             :: sol_dtime_grid = .false.        !> sol-dtime-grid: time-based dumps on the fixed grid t0 + k*sol-dtime (no drift)
    logical             :: sol_overwrite, bck_overwrite    !> switch to overwrite the solution
    character(clen)     :: sol_format(2), bck_format(2), ini_format(2)   !> solution format (native,tecplot,vtk) and formatting (ascii,raw)
    logical             :: write_thermo, write_transport, write_composition
    character(len=llen) :: sol_variables
    character(len=hlen) :: Ovarnames
    character(len=llen) :: ORANSname
    integer             :: Onvar
    real(8)             :: IOtime
    character(llen)     :: meshfile
    character(llen)     :: nameinit
    contains
      procedure, pass(self) :: ini_setup => read_ini
  end type IO_type

  type(IO_type) :: obj_io
  ! NOTE: no "private :: read_ini" here — nvfortran rejects a PRIVATE
  ! statement naming a module procedure that is also a type-bound binding.

  !> Concrete procedure pointing to one of the subroutine realizations
  procedure(r_solution_if), pointer, public :: read_ic
  procedure(w_solution_if), pointer, public :: write_bck, write_solution
  
  !> Abstract interface relative to the finite-rate reactions source procedure
  abstract interface
    subroutine r_solution_if ( IOfield, newrun )
      use Lib_ORION_data

      implicit none
      type(ORION_data), intent(inout) :: IOfield
      logical, intent(in)             :: newrun

    end subroutine r_solution_if


    subroutine w_solution_if ( domain, IOfield, file )
      use IR_Precision
      use MOSE_Types_m
      use Lib_ORION_data
      use MOSE_Parameters_m
      use MOSE_Global_m

      implicit none
      type(MOSE_domain_type), intent(inout)  :: domain
      type(ORION_data), intent(inout)        :: IOfield
      character(llen), intent(in)            :: file

    end subroutine w_solution_if
  end interface

contains


  subroutine read_ini(self, fini)
    use Finer, only: file_ini
    use MOSE_Parameters_m, only: llen, codename
    use MOSE_IO_Wall
    use strings, only: parse

    implicit none
    class(IO_type), intent(inout) :: self
    type(file_ini) :: fini
    ! Local
    character(llen)  :: wholestring
    integer :: error, errini, errbck, errsol

    ! IO
    ! Two different output data are written:
    ! - bck : backup solution, only readable file formats may be employed. Only primitive variables are required
    ! - sol : real solution, not strictly readable file formats may be used. 
    !         Primitive and auxiliary variables are saved together with wall and diagnostic data
    ! If only the sol-format is specified, the sol field is used as bck
    
    ! Defaults
    self%ini_format(1) = 'tecplot'
    self%bck_format(1) = 'native' 
    self%sol_format(1) = 'tecplot'
    self%ini_format(2) = 'ascii'
    self%bck_format(2) = 'binary'
    self%sol_format(2) = 'ascii'
    self%bck_diter = huge(1)
    self%bck_dtime = huge(1.d0)
    self%sol_diter = huge(1)
    self%sol_dtime = huge(1.d0)
    self%bck_overwrite = .true.
    self%sol_overwrite = .true.
    
    call fini%get(section_name=trim(codename)//'-IO', option_name='ini-format', val=wholestring, error=errini)
    if (errini==0) call parse(wholestring,' ', self%ini_format)
    call fini%get(section_name=trim(codename)//'-IO', option_name='bck-format', val=wholestring, error=errbck)
    if (errbck==0) call parse(wholestring,' ', self%bck_format)
    call fini%get(section_name=trim(codename)//'-IO', option_name='sol-format', val=wholestring, error=errsol)
    if (errsol==0) call parse(wholestring,' ', self%sol_format)
    if (errbck/=0 .and. errsol==0) self%bck_format = self%sol_format
    call fini%get(section_name=trim(codename)//'-IO', option_name='bck-diter', val=self%bck_diter, error=error)
    call fini%get(section_name=trim(codename)//'-IO', option_name='bck-dtime', val=self%bck_dtime, error=error)
    call fini%get(section_name=trim(codename)//'-IO', option_name='sol-diter', val=self%sol_diter, error=error)
    call fini%get(section_name=trim(codename)//'-IO', option_name='sol-dtime', val=self%sol_dtime, error=error)
    call fini%get(section_name=trim(codename)//'-IO', option_name='sol-field-every', val=self%sol_field_every, error=error)
    if (error/=0) self%sol_field_every = 1
    self%sol_field_every = max(1, self%sol_field_every)
    call fini%get(section_name=trim(codename)//'-IO', option_name='sol-dtime-grid', val=self%sol_dtime_grid, error=error)
    if (error/=0) self%sol_dtime_grid = .false.
    call fini%get(section_name=trim(codename)//'-IO', option_name='bck-overwrite', val=self%bck_overwrite, error=error)
    call fini%get(section_name=trim(codename)//'-IO', option_name='sol-overwrite', val=self%sol_overwrite, error=error)
    ! Define variables to write
    self%sol_variables = 'thermo'
    call fini%get(section_name=trim(codename)//'-IO', option_name='sol-variables', val=self%sol_variables, error=error)
    wall_variables = 'mech thermal'
    call fini%get(section_name=trim(codename)//'-IO', option_name='wall-variables', val=wall_variables, error=error)

    ! Residual history file
    self%res_diter = 1
    call fini%get(section_name=trim(codename)//'-IO', option_name='diter-res', val=self%res_diter, error=error)

    ! Mesh file name
    call fini%get(section_name=trim(codename)//'-IO', option_name='mesh-file', val=self%meshfile, error=error)
    if (error/=0) self%meshfile = 'MESH/mesh.tec'

  end subroutine read_ini


  subroutine Input_solution_setup( newrun ) 
    use IR_Precision
    use MOSE_Global_m

    implicit none
    logical, intent(in) :: newrun
    ! Local
    logical         :: present
    integer         :: i
    character(llen) :: try
    character(8)    :: extension   !> must fit '.szplt' (6 chars); trim() used below

    
    if (newrun) then
      select case (obj_io%ini_format(1))
        case('native')
          extension = '.bin'
          read_ic  => read_native
        case('vtk')
          extension = '.vtm'
          read_ic  => read_vtk_tec
        case('tecplot')
        if (obj_io%ini_format(2)=='ascii') then
          extension = '.tec'
          else
          extension = '.szplt'
          endif
          read_ic  => read_vtk_tec
      end select
      obj_io%nameinit = 'INPUT/'//trim(MOSE_phase_prefix)//'ic'//trim(extension)
      inquire(file=obj_io%nameinit, exist=present)

    else
      select case (obj_io%bck_format(1))
      case('native')
        extension = '.bin'
        read_ic => read_native
      case('vtk')
        extension = '.vtm'
        read_ic => read_vtk_tec
      case('tecplot')
       if (obj_io%bck_format(2)=='ascii') then  ! restart branch: was wrongly reading ini_format(2)
          extension = '.tec'
        else
          extension = '.szplt'
        endif
        read_ic => read_vtk_tec
      end select
      obj_io%nameinit = 'OUTPUT/'//trim(MOSE_phase_prefix)//'field'//trim(extension)
      inquire(file=obj_io%nameinit, exist=present)

      if (.not.present) then
        i = 0
        do
          i = i+1
          try = 'OUTPUT/'//trim(MOSE_phase_prefix)//'field'//trim(str(.true.,i))//trim(extension)
          inquire(file=try,exist=present)
          if (present) then
            obj_io%nameinit = try
          else
            exit
          endif
        enddo
      endif

    endif

    select case (obj_io%bck_format(1))
      case('native')
        write_bck => write_native
      case('vtk')
        write_bck => write_vtk_tec
      case('tecplot')
        write_bck => write_vtk_tec
    end select

  end subroutine Input_solution_setup


  !> Using the native format, the mesh data and the field data are stored into two different files:
  !> - MESH/mesh.tec
  !> - INPUT/field.bin or OUTPUT/field<n>.bin
  !> As a consequence, the mesh is read via ORION tecplot library whereas the field via internal commands
  subroutine read_native ( IOfield, newrun )
    use Lib_Tecplot
    use Lib_ORION_data
    use MOSE_Global_m

    implicit none
    type(ORION_data), intent(inout) :: IOfield
    logical, intent(in)             :: newrun
    ! Local
    integer :: unitfile
    integer :: i, j, k, b, s, nb, nvar, ios


    !% Read the mesh file and store the geometric data into an orion-field
    IOfield%tec%format = 'ascii'
    ios = tec_read_structured_multiblock(orion=IOfield,filename=obj_io%meshfile)
    nb = size ( IOfield%block )

    !% Allocate orion-field vars
    
    nvar = nsc+4+nrans
    do b = 1, nb
      allocate(IOfield%block(b)%vars(1:nvar,1:IOfield%block(b)%Ni, &
                                            1:IOfield%block(b)%Nj, &
                                            1:IOfield%block(b)%Nk) )
    enddo

    !% Read the field data and store them into an orion-field
    open ( newunit=unitfile, form='unformatted', file=obj_io%nameinit, status='old', convert='big_endian' )
    read(unitfile) nb
    do b = 1, nb; read(unitfile); end do
    do b = 1, nb; read(unitfile); end do
    do b = 1, nb; read(unitfile); end do
    ! Species density
    do b = 1, nb
      read(unitfile) ( ( ( ( IOfield%block(b)%vars(s,i,j,k), s=1,nsc ), &
                             i=1,IOfield%block(b)%Ni ), &
                             j=1,IOfield%block(b)%Nj ), &
                             k=1,IOfield%block(b)%Nk )
    end do
    ! U, V, W, p
    do s = nsc+1, nsc+4
      do b = 1, nb
        read(unitfile) ( ( ( IOfield%block(b)%vars(s,i,j,k), &
                             i=1,IOfield%block(b)%Ni ), &
                             j=1,IOfield%block(b)%Nj ), &
                             k=1,IOfield%block(b)%Nk )
      end do
    enddo
    ! RANS variables
    if ( nrans > 0 ) then
      do b = 1, nb
        read(unitfile) ( ( ( ( IOfield%block(b)%vars(s,i,j,k), s=nsc+5,nvar ), &
                               i=1,IOfield%block(b)%Ni ), &
                               j=1,IOfield%block(b)%Nj ), &
                               k=1,IOfield%block(b)%Nk )
      end do
    end if
    ! Time
    read(unitfile) obj_io%IOtime
    close(unitfile)

  end subroutine read_native


  subroutine read_vtk_tec ( IOfield, newrun )
    use Lib_VTK
    use Lib_Tecplot
    use Lib_ORION_data
    use MOSE_Mod_MPI, only: mpi_abort_all

    implicit none
    type(ORION_data), intent(inout) :: IOfield
    logical, intent(in)             :: newrun
    ! Local
    integer:: error
    character(clen) :: format_1, format_2


    if ( newrun ) then
      format_1 = obj_io%ini_format(1)
      format_2 = obj_io%ini_format(2)
    else
      format_1 = obj_io%bck_format(1)
      format_2 = obj_io%bck_format(2)
    endif

    obj_io%IOtime = 0.d0

    select case(trim(format_1))
    case('tecplot')
      IOfield%tec%format = trim(format_2)
      error = tec_read_structured_multiblock(orion=IOfield,filename=trim(obj_io%nameinit))
    case('vtk')
      IOfield%tec%format = trim(format_2)
      error = vtk_read_structured_multiblock(orion=IOfield,vtmpath=obj_io%nameinit(1:len(trim(obj_io%nameinit))-4),vtspath='INPUT/vtk/field',time=obj_io%IOtime)
    end select

    if (error/=0) call mpi_abort_all("[ERROR] reading input file "//trim(obj_io%nameinit))

  end subroutine read_vtk_tec


  !> Output setup
  subroutine Output_Solution_Setup ( IOfield, MGL )
    use Lib_ORION_data
    use IR_Precision
    use MOSE_Global_m

    implicit none
    integer, intent(in)             :: MGL
    type(ORION_data), intent(inout) :: IOfield(MGL)
    ! Local
    integer :: b, i, m


    ! IO Variables specification
    obj_io%Ovarnames=' '
    obj_io%ORANSname=' '
    if ( nrans == 1 ) then
      obj_io%ORANSname=' "mi_tilde"'
    elseif ( nrans == 2 ) then
      obj_io%ORANSname=' "kappa" "omega"'
    elseif ( nrans == 7 ) then
      obj_io%ORANSname=" ru'u' rv'v' rw'w' ru'v' ru'w' rv'w' omega"
    end if
    do i = 1, nsc
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//'"rho('//trim(str(.true.,i))//')"'
    enddo
    obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "u" "v" "w" "p"'
    !Soot Variables
      if ( nsoot == 1 ) then
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "rho*Np"'
    end if
    if ( nsoot == 2 ) then
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "rho*Yp" "rho*Np"'
    end if
    !Passive scalars
    if (npass>0) then
      do i = 1, npass
        obj_io%Ovarnames = trim(obj_io%Ovarnames)//'"Pass('//trim(str(.true.,i))//')"'
      enddo
    endif 
   
     obj_io%Ovarnames = trim(obj_io%Ovarnames)//obj_io%ORANSname

    ! Auxiliary variables
    ! T is printed by default
    obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "T"'
    obj_io%Onvar = nprim + 1 ! T
    ! Thermodynamic properties
    obj_io%write_thermo = .false.
    if (index(obj_io%sol_variables,'thermo')>0) then
      obj_io%write_thermo = .true.
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "g" "R"'
      obj_io%Onvar = obj_io%Onvar + 2
    endif
    ! Transport properties
    obj_io%write_transport = .false.
    if (index(obj_io%sol_variables,'transport')>0) then
      obj_io%write_transport = .true.
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "mil" "kl"'
      obj_io%Onvar = obj_io%Onvar + 2
    endif
    if (obj_io%write_transport .and. (nrans>0)) then
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "mit"'
      obj_io%Onvar = obj_io%Onvar + 1
    endif
    ! Composition
    ! write_composition = .false.
    ! if (index(sol_variables,'compos')>0) then
    !   write_composition = .true.
    !   Ovarnames = trim(Ovarnames)//' "mu" "kl"'
    !   Onvar = Onvar + 2
    ! endif

    do m = 1, MGL
      ! VTK specification
      IOfield(m)%vtk%format = obj_io%sol_format(2)
      IOfield(m)%vtk%node = .false.

      ! Tecplot specification
      IOfield(m)%tec%node = .false.
      IOfield(m)%tec%format = obj_io%sol_format(2)
      IOfield(m)%tec%bc = .false.
       if (obj_io%sol_format(2)=='ascii') then
        IOfield(m)%tec%extension = '.tec'
      else
        IOfield(m)%tec%extension = '.szplt'
      endif
    end do

    ! Concretize the sol subroutine
    write_solution => write_vtk_tec

    do m = 1, MGL
      do b = 1, Size(IOfield(m)%block)
        IOfield(m)%block(b)%name = 'Block'//trim(str(.true.,b))
      enddo
    enddo

    ! Reallocate IOfield vars if the number of variables read into the backup file is different from the solution one
    ! The reallocation is performed after reading the ICs and before the first solution is written!
    do m = 1, MGL
      if ( obj_io%Onvar /= Size(IOfield(m)%block(1)%vars,1) ) then
        do b = 1, Size ( IOfield(m)%block )
          IOfield(m)%block(b)%name = 'Block'//trim(str(.true.,b))
          deallocate ( IOfield(m)%block(b)%vars )
          allocate( IOfield(m)%block(b)%vars(1:obj_io%Onvar, &
                                             1:IOfield(m)%block(b)%Ni, &
                                             1:IOfield(m)%block(b)%Nj, &
                                             1:IOfield(m)%block(b)%Nk) )
        enddo
      endif
    enddo

  end subroutine Output_Solution_Setup


  !> Directly write the domain data into a binary file using a native format
  subroutine write_native ( domain, IOfield, file )
    use IR_Precision
    use MOSE_Types_m
    use Lib_ORION_data
    use MOSE_Parameters_m
    use MOSE_Global_m
    use MOSE_Mod_MPI, only: mpi_is_root
    use MOSE_Mod_GhostExchange, only: gather_P_to_root, mpi_io_barrier

    implicit none
    type(MOSE_domain_type), intent(inout)  :: domain
    type(ORION_data), intent(inout)        :: IOfield
    character(llen), intent(in)            :: file
    ! Local
    character(len=llen) :: path
    integer :: unitfile, i, j, k, b, s, nb

    call gather_P_to_root(domain)

    if (mpi_is_root) then
      path = 'OUTPUT/'
      open ( newunit=unitfile, form='unformatted', file=trim(path)//trim(file)//'.bin', convert='big_endian' )

      nb = domain%nb
      write(unitfile) nb
      do b = 1, nb; write(unitfile) domain%blk(b)%dim(1); end do
      do b = 1, nb; write(unitfile) domain%blk(b)%dim(2); end do
      do b = 1, nb; write(unitfile) domain%blk(b)%dim(3); end do

      do b = 1, nb
        write(unitfile) ( ( ( ( domain%blk(b)%P(s,i,j,k), s=1,nsc ), i=1,domain%blk(b)%dim(1) ), j=1,domain%blk(b)%dim(2) ), k=1,domain%blk(b)%dim(3) )
      end do
      do s = nsc+1, np
        do b = 1, nb
          write(unitfile) ( ( ( domain%blk(b)%P(s,i,j,k), i=1,domain%blk(b)%dim(1) ), j=1,domain%blk(b)%dim(2) ), k=1,domain%blk(b)%dim(3) )
        end do
      enddo
      if ( nRANS > 0 ) then
        do b = 1, nb
          write(unitfile) ( ( ( ( domain%blk(b)%P(s,i,j,k), s=nt,nprim ), i=1,domain%blk(b)%dim(1) ), j=1,domain%blk(b)%dim(2) ), k=1,domain%blk(b)%dim(3) )
        end do
      endif

      write(unitfile) domain%time

      close(unitfile)
    end if

    ! Synchronize all ranks after I/O
    call mpi_io_barrier()

  end subroutine write_native


  !> Update the orion-field data and write them accordingly to the chosen format (vtk,tecplot)
  subroutine write_vtk_tec ( domain, IOfield, file )
    use IR_Precision
    use MOSE_Types_m
    use Lib_ORION_data
    use MOSE_Parameters_m
    use MOSE_Global_m
    use Lib_VTK
    use Lib_Tecplot
    use MOSE_Lib_RANS
    use FLINT_Lib_Thermodynamic
    use MOSE_Mod_MPI, only: mpi_is_root
    use MOSE_Mod_GhostExchange, only: gather_P_to_root, mpi_io_barrier

    implicit none
    type(MOSE_domain_type), intent(inout)  :: domain
    type(ORION_data), intent(inout)        :: IOfield
    character(llen), intent(in)            :: file
    ! Local
    character(len=llen) :: path
    character(len=llen) :: localpath_vtk
    integer             :: E_IO, b, i, j, k, pv
    real(8)             :: p, rho, Rtot, rhoi(nsc), T, mu, kl, mut
    real(8)             :: vel_gradient(3,3)

    call gather_P_to_root(domain)

    if (mpi_is_root) then
      path = 'OUTPUT/'

      ! Update IOfield variables with domain primitives and other variables
      do b = 1, size(IOfield%block)
        IOfield%block(b)%vars(1:nprim,:,:,:) = domain%blk(b)%P(:,1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
        ! Auxiliary variables
        do k = 1, IOfield%block(b)%Nk ; do j = 1, IOfield%block(b)%Nj ; do i = 1, IOfield%block(b)%Ni
          rhoi = domain%blk(b)%P(1:nsc,i,j,k)
          p    = domain%blk(b)%P(np,i,j,k)
          rho  = sum ( rhoi )
          Rtot = f_Rtot( rhoi )
          T = eos(p=p,ro=rho,R=Rtot)
          pv = 1; IOfield%block(b)%vars(nprim+1,i,j,k) = T
          if (obj_io%write_thermo) then
            pv = pv + 1; IOfield%block(b)%vars(nprim+pv,i,j,k) = f_gamma ( rhoi, p, rho, Rtot )
            pv = pv + 1; IOfield%block(b)%vars(nprim+pv,i,j,k) = Rtot
          endif
          if (obj_io%write_transport) then
            call co_k_mi_lam_Wilke(rhoi,rho,T,mu,kl)
            pv = pv +1; IOfield%block(b)%vars(nprim+pv,i,j,k) = mu
            pv = pv +1; IOfield%block(b)%vars(nprim+pv,i,j,k) = kl
            if (nrans>0) then
              vel_gradient(:,1) = ( domain%blk(b)%P(nu:nw,i+1,j,k) - domain%blk(b)%P(nu:nw,i-1,j,k) )/2d0
              vel_gradient(:,2) = ( domain%blk(b)%P(nu:nw,i,j+1,k) - domain%blk(b)%P(nu:nw,i,j-1,k) )/2d0
              vel_gradient(:,3) = ( domain%blk(b)%P(nu:nw,i,j,k+1) - domain%blk(b)%P(nu:nw,i,j,k-1) )/2d0
              vel_gradient = matmul ( vel_gradient + 1d-40, domain%blk(b)%m(i,j,k)%c )
              call Eddy_Viscosity ( mut=mut, rans_variables=domain%blk(b)%p(nt:,i,j,k), &
                                    mul=mu, rho=rho, vel_gradient=vel_gradient, &
                                    walldist=domain%blk(b)%yn(i,j,k))
              pv = pv +1; IOfield%block(b)%vars(nprim+pv,i,j,k) = mut
            end if
          endif
        enddo; enddo; enddo
      enddo

      ! Write the IOfield accordingly to the solution format
      select case(obj_io%sol_format(1))
      case('vtk')
        IOfield%vtk%format = obj_io%sol_format(2)
        localpath_vtk = trim(path)//'vtk/'
        call execute_command_line('mkdir -p '//trim(localpath_vtk))
        E_IO = vtk_write_structured_multiblock(orion=IOfield,vtspath=trim(localpath_vtk)//trim(file), &
                                                             vtmpath=trim(path)//trim(file),varnames=obj_io%Ovarnames,time=domain%time)
      case('tecplot')
        IOfield%tec%format = obj_io%sol_format(2)
        E_IO = tec_write_structured_multiblock(orion=IOfield,varnames=obj_io%Ovarnames,filename=trim(path)//trim(file)//trim(IOfield%tec%extension))
      end select
    end if

    ! Synchronize all ranks after I/O
    call mpi_io_barrier()

  end subroutine write_vtk_tec


end module MOSE_IO_Solution
