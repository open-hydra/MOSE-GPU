module MOSE_IO_BC
    use MOSE_Types_m

  implicit none


contains

  subroutine Setup_BC ( domain, MGL )
    use MOSE_IO_Wall
    use MOSE_Mod_BC_Q2D, only: Setup_Q2D_BC_Data

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain(MGL)
    integer, intent(in)                   :: MGL
    ! Local
    integer :: m

    do m = 1, MGL
      call Allocate_BC ( domain(m) )
      
      if (allocated(viscous_flag)) deallocate(viscous_flag)
      allocate ( viscous_flag( 6 * domain(m) % nb, 6 ) )
      viscous_flag = .false.

      if (allocated(coupling_flag)) deallocate(coupling_flag)
      allocate ( coupling_flag( 6 * domain(m) % nb, 6 ) )    
      coupling_flag = .false.
    
      call Read_BCfile ( domain(m) % bc, domain(m) % nbound, domain(m) % n_bf, m )

      if (m == 1) call Setup_Q2D_BC_Data ( domain(m) )

      if ((any(viscous_flag)) .and. (.not.VISC)) then
        write(*,*) ' [ERROR] Viscous wall BCs assigned with Eulerian model'
        stop
      end if
    
    end do

  end subroutine Setup_BC


  subroutine Allocate_BC ( domain )
    
    implicit none
    class(MOSE_domain_type) :: domain
    ! Local
    integer :: b, n, ni, nj, nk

    n = 0

    do b = 1, domain % nb

      ni = domain % blk(b) % dim(1)
      nj = domain % blk(b) % dim(2)
      nk = domain % blk(b) % dim(3)

      n = n + 2*nj*nk + 2*ni*nk + 2*nj*ni

    end do

    domain % nbound = n
    allocate ( domain % bc ( domain % nbound ) )
    allocate ( domain % n_bf ( domain % nb, 6 ) )
    
  end subroutine Allocate_BC


  subroutine Read_BCfile ( bc, n, n_bf, level )
    use MOSE_Global_m
    use MOSE_IO_Wall
    use MOSE_Mod_GSI, only : setup_GSI
    use IR_Precision
    use MOSE_Mod_MPI
    use MOSE_IO_BC_Open, only: BC_DIALECT_LEGACY, BC_DIALECT_OPEN, Detect_BC_Dialect, &
                               Skip_BC_Record, Read_BC_Record_Open
    implicit none
    type(MOSE_bc_type), dimension(n), intent(inout) :: bc
    integer, intent(in)                             :: n, level
    integer, dimension(1:,1:), intent(inout)        :: n_bf
    ! Local
    integer :: cc, i, s, nconnect, nwall, ngsi, nSRM, nio, nsym, nper, nchimera, nnozzle, ncoupled, next, nmanifold, nKaffs, nq2d
    integer :: di(5), ti
    integer :: unitfile, ios, cios, n_proof, ip
    integer :: dialect
    character(len=32) :: p0file
    character(len=256) :: q2d_line


    ! Open file
    if (level == 1) then
      open(newunit=unitfile,file='INPUT/'//trim(MOSE_phase_prefix)//'bc.txt',status='old',iostat=ios,action='read')
    else
      open(newunit=unitfile,file='INPUT/'//trim(MOSE_phase_prefix)//'bc'//trim(str(.true.,level))//'.txt',status='old',iostat=ios,action='read')
    endif
    if (ios/=0) error stop 'BC file not found'

    ! Counters for specific BC types
    cios = 0
    nconnect = 0
    nwall = 0
    ngsi = 0
    nSRM = 0
    nio = 0
    nsym = 0
    nper = 0
    nchimera = 0
    nnozzle = 0
    ncoupled = 0
    nmanifold = 0
    next = 0
    nKaffs = 0
    nq2d = 0

    ! Counter for number of cells per face in each block
    n_bf = 0

    ! Which dialect is this file written in? Legacy MOSE ids, or the 3-digit
    ! MOSE_open/ATLAS ids (translated record by record below).
    dialect = Detect_BC_Dialect ( unitfile )
    if (mpi_is_root .and. level == 1 .and. dialect == BC_DIALECT_OPEN) then
      write(*,*)
      write(*,'(A)') ' Boundary condition file in MOSE_open (ATLAS) format'
    endif

    ! Cheak BC file consistency
    ios = 0; n_proof = -1
    do while (ios==0)
      read( unitfile,*,iostat=ios ) di(1), di(2), di(3), di(4), di(5), ti
      call Skip_BC_Record ( unitfile, ti, dialect, ios )
      n_proof = n_proof + 1
    enddo

    if (n_proof /= n) then
      write(*,'(A84)') "[ERROR] Boundary conditions mesh is different than the one of the initial conditions"
      write(*,*) "        Number of assigned BCs in the file = ",n_proof
      write(*,*) "        Number of boundary faces in the IC = ",n
      stop
    endif

    ! Read file
    rewind(unitfile)
    do i = 1, n
       
      ! First line is equal for every BC type
      read( unitfile,*,iostat=ios ) bc(i)%b, bc(i)%i, bc(i)%j, bc(i)%k, bc(i)%f, bc(i)%type
      if (ios/=0) write(*,*) '  Error in BC file'

      ! n_bf update
      n_bf( bc(i) % b, bc(i) % f ) = n_bf( bc(i) % b, bc(i) % f ) + 1

      ! Second line is BC type-dependent. In the MOSE_open dialect the whole
      ! record is read and translated to the legacy id/fields first, so only the
      ! counters (summary print + WALL/CHIMERA/GSI flags) are left to update.
      if (dialect == BC_DIALECT_OPEN) then

        call Read_BC_Record_Open ( unitfile, bc(i), ios )

        select case( bc(i)%type )
          case(1);            nconnect = nconnect + 1
          case(2);            nper     = nper     + 1
          case(3);            nsym     = nsym     + 1
          case(4,24,25);      nio      = nio      + 1
          case(5,6);          nwall    = nwall    + 1
          case(8,9,10,12,13); nwall    = nwall    + 1; ngsi = ngsi + 1
          case(11);           next     = next     + 1
          case(14);           nwall    = nwall    + 1; nSRM = nSRM + 1
          case(22);           nnozzle  = nnozzle  + 1
          case(666);          nchimera = nchimera + 1
          case(667);          nq2d     = nq2d     + 1
          case(1000);         ncoupled = ncoupled + 1
          case(1001);         nmanifold = nmanifold + 1
        end select

        if (bc(i)%p0time%exists) cios = 1        ! drives the 'Time-dependent BC' notice

        cycle

      endif

      select case( bc(i)%type )
        case(1)
          nconnect = nconnect + 1
          read( unitfile,*,iostat=ios ) &
          bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22
          allocate ( bc(i) % Pg (nprim, 6) )
        
        case(2)
          nper = nper + 1
        
        case(3)
          nsym = nsym + 1
        
        case(667)
          nq2d = nq2d + 1
          allocate ( bc(i) % Pg (nprim, 6) )
          read(unitfile, '(A)', iostat=ios) q2d_line
          q2d_line = adjustl(q2d_line)
          ip = index(trim(q2d_line), ' ')
          if (ip > 0) then
            bc(i)%q2d_file = q2d_line(1:ip-1)
            if (index(q2d_line(ip:), 'periodic') > 0) bc(i)%q2d_periodic = .true.
          else
            bc(i)%q2d_file = trim(q2d_line)
          endif

        case(4)
          nio = nio + 1
          allocate( bc(i) % ci(1 : nsc+nrans) )
          read( Unitfile,*,iostat=ios ) &
          bc(i)%mach, bc(i)%T0, p0file, bc(i)%alpha, bc(i)%beta, bc(i)%pamb, bc(i)%rel_fac, (bc(i)%ci(s), s = 1, nsc+nrans)
          ! Time-varying properties: if the read string is not convertible to a real, the BC is time-dependent
          read(p0file,*,iostat=cios) bc(i)%p0
          if (cios/=0) call bc(i)%p0time%initialize(file=p0file,bar=.true.)

        case(5)
          nwall = nwall + 1
          viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          read(unitfile,*,iostat=ios) bc(i)%qw
        
        case(6)
          nwall = nwall + 1
          viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          read(unitfile,*,iostat=ios) bc(i)%Tw
        
        case(8)
          nwall = nwall + 1
          ngsi = ngsi + 1
          viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          read(unitfile,*,iostat=ios) bc(i)%Tw, bc(i)%qrad

        case(9,10,12,13)
          nwall = nwall + 1
          ngsi = ngsi + 1
          viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          read(unitfile,*,iostat=ios) bc(i)%qrad

        case(11)
          next = next + 1

        case(14)
          nwall = nwall + 1
          nSRM = nSRM + 1
          allocate( bc(i) % ci(1 : nsc) )
          read(unitfile,*,iostat=ios) bc(i)%aCoeff, bc(i)%n, bc(i)%pRef, bc(i)%Taf, bc(i)%haf, (bc(i)%ci(s), s = 1, nsc)

        case(22)
          nnozzle = nnozzle + 1
          allocate( bc(i) % ci(1 : nsc+nrans) )
          read( unitfile,*,iostat=ios ) &
          bc(i)%mach, bc(i)%T0, bc(i)%p0, bc(i)%alpha, bc(i)%beta, bc(i)%pamb, bc(i)%rel_fac, (bc(i)%ci(s), s = 1, nsc+nrans)

        case(666)
          nchimera = nchimera + 1
          read( unitfile,*,iostat=ios ) (bc(i)%ni(cc),cc=1,2)
          allocate(bc(i)%donorID(1:sum(bc(i)%ni),1:4))
          allocate(bc(i)%volume_fraction(1:sum(bc(i)%ni)))
          do s = 1, bc(i)%ni(1)
            read( unitfile,*,iostat=ios ) bc(i)%donorID(s,1:4), bc(i)%volume_fraction(s)
          enddo
          do s = bc(i)%ni(1)+1, bc(i)%ni(1)+bc(i)%ni(2)
            read( unitfile,*,iostat=ios ) bc(i)%donorID(s,1:4), bc(i)%volume_fraction(s)
          enddo
          allocate ( bc(i) % Pg (nprim, 6) )

        case(1000)
          ncoupled = ncoupled + 1
          coupling_flag( bc(i)%b , bc(i)%f ) = .true.
          read( unitfile,*,iostat=ios ) &
          bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22
          allocate(bc(i)%ext_flux(nprim))
          bc(i)%ext_flux = 0.0
          allocate ( bc(i) % Pg (1, 6) )
        
        case(1001)
          nmanifold = nmanifold + 1
          read( unitfile,*,iostat=ios ) bc(i)%bs, bc(i)%fs
          
        case(999)
          nKaffs = nKaffs + 1
           read(unitfile,*,iostat=ios) bc(i)%bs, bc(i)%fs

      end select

    end do
    ! Set global flags on ALL ranks (not just root) to avoid MPI deadlocks
    if (level == 1) then
      if (nwall > 0) WALL = .true.
      if (nchimera > 0) CHIMERA = .true.
      if (ngsi > 0) call Setup_GSI
    endif

    if (mpi_is_root) then
    if (level == 1) then
    write(*,*)
    write(*,'(A)') ' Boundary Conditions:'
    write(*,'(A,I8)') '  - Connection      --> ', nconnect
    write(*,'(A,I8)') '  - Viscous wall    --> ', nwall
    if (ngsi > 0) then
      write(*,'(A,I8)') '  - GSI (inc. wall) --> ', ngsi
    endif
    write(*,'(A,I8)') '  - Inflow/outflow  --> ', nio
    write(*,'(A,I8)') '  - Symmetry        --> ', nsym
    write(*,'(A,I8)') '  - Periodicity     --> ', nper
    write(*,'(A,I8)') '  - Extrapolation   --> ', next
    if (nchimera > 0) then
      write(*,'(A,I8)') '  - Chimera         --> ', nchimera
    endif
    if (ncoupled > 0) then
      write(*,'(A,I8)') '  - Coupled wall    --> ', ncoupled
    endif
    if (nSRM > 0) then
        write(*,'(A,I8)') '  - SRM wall        --> ', nSRM
    endif
    if (nkaffs > 0) then
       write(*,'(A,I8)') '  - Coupled 1D-MultiD    --> ', nkaffs
    endif
    if (nq2d > 0) then
       write(*,'(A,I8)') '  - Q2D Mapped      --> ', nq2d
    endif
  endif
    if (cios/=0) write(*,*) '  Time-dependent BC'
   endif
    close( Unitfile )

  end subroutine Read_BCfile


end module MOSE_IO_BC
