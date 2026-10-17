!> @ingroup Library
!> @{
!> @defgroup Lib_TecplotLibrary Lib_Tecplot
!> @}

!> Pure Fortran (2003+) library to write data conforming the Tecplot standard
!> @ingroup Lib_TecplotLibrary
module Lib_Tecplot
  use, intrinsic:: ISO_FORTRAN_ENV, only: stdout => OUTPUT_UNIT, stderr => ERROR_UNIT ! Standard output/error logical units.
  use IR_Precision
  use Lib_ORION_data
  implicit none
# if defined(TECIO)
  include "tecio.f90"
# endif
  private

  public:: tec_write_structured_multiblock
  public:: tec_read_structured_multiblock
  public:: tec_write_points_multivars
  public:: tec_read_points_multivars


contains

  ! Subroutine for computing the dimensions of the domain to be post-processed. 
  subroutine compute_dimensions(node,bc,Nx,Ny,Nz,gc,ni1,ni2,nj1,nj2,nk1,nk2,ci1,ci2,cj1,cj2,ck1,ck2)
    implicit none
    logical, intent(in) :: node 
    logical, intent(in) :: bc 
    integer, intent(IN) :: Nx, Ny, Nz, gc          ! Block-level data.
    integer, intent(OUT):: ni1,ni2,nj1,nj2,nk1,nk2 ! Bounds of dimensions of node-centered data.
    integer, intent(OUT):: ci1,ci2,cj1,cj2,ck1,ck2 ! Bounds of dimensions of cell-centered data.
    if (node) then
      ni1 = 0 ; ni2 = Nx
      nj1 = 0 ; nj2 = Ny
      nk1 = 0 ; nk2 = Nz

      ci1 = 1 ; ci2 = Nx
      cj1 = 1 ; cj2 = Ny
      ck1 = 1 ; ck2 = Nz
    else
      if (bc) then
        ni1 = 0 - gc ; ni2 = Nx + gc
        nj1 = 0 - gc ; nj2 = Ny + gc
        nk1 = 0 - gc ; nk2 = Nz + gc

        ci1 = 1 - gc ; ci2 = Nz + gc
        cj1 = 1 - gc ; cj2 = Ny + gc
        ck1 = 1 - gc ; ck2 = Nz + gc
      else
        ni1 = 0      ; ni2 = Nx
        nj1 = 0      ; nj2 = Ny
        nk1 = 0      ; nk2 = Nz

        ci1 = 1      ; ci2 = Nx
        cj1 = 1      ; cj2 = Ny
        ck1 = 1      ; ck2 = Nz
      endif
    endif
    ! If 2D domain -> nk = 0; ck = 1 
    if (nk1+nk2==0) then
      ck2 = 1; ck1 = 1
    endif
  endsubroutine compute_dimensions
  !> @}

  !> @ingroup Lib_PostProcessingPublicProcedure
  !> @{
  !> Function for writing ORION block data to Tecplot file.
  function tec_write_structured_multiblock(orion,varnames,filename,Nvars) result(err)
    implicit none
    type(orion_data), intent(in)              :: orion
    character(len=*), intent(in), optional    :: varnames
    character(len=*), intent(in)              :: filename !< File name of the output file.
    integer, intent(in), optional             :: Nvars    !< Input number of variables saved.
    logical :: meshonly, time_accurate
    integer :: err
# if defined(TECIO)
    integer, external::               tecini142,    &     ! |
                                      tecauxstr142, &     ! |
                                      teczne142,    &     ! | Tecplot external functions.
                                      tecdatd142,   &     ! |
                                      tecend142           ! |
# endif
    character(1), parameter:: tecendrec = char(0) !< End-character for binary-record end.
    character(500)::          tecvarname          !< Variables name for tecplot header file.
    character(500)::          teczoneheader       !< Tecplot string of zone header.
    character(500)::          tecvarform          !< Format for variables for tecplot file.
    integer, allocatable::    tecvarloc(:)        !< Tecplot array of variables location.
    character(500)::          tecvarlocstr        !< Tecplot string of variables location.
    integer, allocatable::    tecnull(:)          !< Tecplot null array.
    integer::                 tecunit             !< Free logic unit of tecplot file.
    integer::                 Debug      = 0
    integer::                 VIsDouble  = 0
    integer::                 FileType   = 0
    integer::                 fileFormat ! 0 == PLT, 1 == SZPLT
    integer::                 Nvar                !< Internal number of variables saved.
    integer::                 ndir                !< Internal number of dimensions.
    integer::                 Nblocks             !< Number of blocks.
    integer::                 b                   !< Counter.
    integer::                 NvarTot
    integer::                 gc
    real(R8P)::               time_
  
    ! Preliminary operations
    FileType   = 0
    Debug      = 1
    VIsDouble  = 0
    if (allocated(orion%block(1)%vars) .and. .not.present(Nvars)) then
      meshonly = .false.
      Nvar = size(orion%block(1)%vars,1)
    elseif (allocated(orion%block(1)%vars) .and. present(Nvars)) then
      meshonly = .false.
      Nvar = Nvars
    elseif (.not.allocated(orion%block(1)%vars)) then
      meshonly = .true.
      Nvar = 0
    endif
    Nblocks = size(orion%block)
    ndir = size(orion%block(1)%mesh, 1)
    gc = 1
    ! allocating dynamic arrays
    allocate(tecvarloc(1:ndir+Nvar))
    allocate(tecnull(1:ndir+Nvar))
    ! initializing tecplot variables
    call tec_init()
    ! time
    time_accurate = .true.
    time_ = orion%solutiontime
    if (time_<0) then
      time_ = abs(time_)
      time_accurate = .false.
    endif
    ! initializing tecplot file
    select case(orion%tec%format)
    case('binary')
#   if defined(TECIO)
      if (index(filename,'.plt')>0) then
        fileFormat = 0
      elseif (index(filename,'.szplt')>0) then
        fileFormat = 1
      else
        write(stderr,'(A)')'Tecplot binary file must end with ".plt" or ".szplt"'
        return
      endif
      err = tecini142(tecendrec,trim(tecvarname)//tecendrec,trim(filename)//tecendrec,'.'//tecendrec,fileFormat,FileType,Debug,VIsDouble)
      err = tecauxstr142("Time"//tecendrec,trim(str(n=time_))//tecendrec)
#   else
      stop "You can not write in binary formato without compiling against TecIO"
#   endif
    case('ascii')
      open(newunit=tecunit,file=trim(filename))
      write(tecunit,'(A)',iostat=err)trim(tecvarname)
    end select
    ! writing data blocks
    do b=1,Nblocks
      err = tec_blk_data(b = b)
    enddo
    ! finalizing tecplot file
    select case(orion%tec%format)
    case('binary')
# if defined(TECIO)
      err = tecend142()
# endif
    case('ascii')
      close(tecunit)
    end select
    ! deallocating dynamic arrays
    deallocate(tecvarloc)
    deallocate(tecnull)

  contains
    subroutine tec_init()
      implicit none
      integer  :: s ! Counter.
      if (.not.meshonly) then
        NvarTot = ndir + Nvar
      else
        NvarTot = ndir
      endif
      select case(orion%tec%format)
      case('binary')
#     if defined(TECIO)
        ! header variables names
        tecvarname = 'x'
        if (ndir>=2) tecvarname = trim(tecvarname)//' y'
        if (ndir==3) tecvarname = trim(tecvarname)//' z'
        if (.not.meshonly) then
          if (present(varnames)) then
            tecvarname = trim(tecvarname)//' '//trim(varnames)
          else
            do s = 1, Nvar
              tecvarname = trim(tecvarname)//' "var'//trim(str(.true.,s))//'"'
            enddo
          endif
        endif
        ! variables location
        if (orion%tec%node) then
          tecvarloc = 1
        else
          tecvarloc(1:ndir) = 1 ; tecvarloc(ndir+1:ndir+Nvar)= 0
        endif
        ! null array
        tecnull = 0
#     else
        stop "You can not write in binary format without compiling against TecIO"
#     endif
      case('ascii')
        ! header variables names
        tecvarname = ' VARIABLES ="x"'
        if (ndir>=2) tecvarname = trim(tecvarname)//' "y"'
        if (ndir==3) tecvarname = trim(tecvarname)//' "z"'
        if (.not.meshonly) then
          if (present(varnames)) then
            tecvarname = trim(tecvarname)//' '//trim(varnames)
          else
            do s = 1, Nvar
              tecvarname = trim(tecvarname)//' "var'//trim(str(.true.,s))//'"'
            enddo
          endif
        endif
        ! variables output format
        if (.not.meshonly) then
          write(tecvarform,'(A)')'('//trim(str(no_sign=.true.,n=Nvar))//'('//FR_P//',1X))'
        else
          write(tecvarform,'(A)')'('//trim(str(no_sign=.true.,n=3))//'('//FR_P//',1X))'
        endif
        ! variables location
        if (.not.meshonly) then
          if (orion%tec%node) then
            tecvarlocstr = ', VARLOCATION=([1-'//trim(str(.true.,NvarTot))//']=NODAL)'
          else
            if (Nvar==1) then
              tecvarlocstr = ', VARLOCATION=([1-'//trim(str(.true.,ndir))//']=NODAL,['//trim(str(.true.,ndir+1))//']=CELLCENTERED)'
            else
              tecvarlocstr = ', VARLOCATION=([1-'//trim(str(.true.,ndir))//']=NODAL,['//trim(str(.true.,ndir+1))//'-'//trim(str(.true.,NvarTot))//']=CELLCENTERED)'
            endif
          endif
        else
          tecvarlocstr = ', VARLOCATION=([1-'//trim(str(.true.,ndir))//']=NODAL)'
        endif
      end select
    endsubroutine tec_init

    function tec_blk_data(b) result(err)
      implicit none
      integer, intent(IN):: b           ! Block number.
      integer:: err                     ! Error trapping flag: 0 no errors, >0 error occurs.
      integer:: ni1,ni2,nj1,nj2,nk1,nk2 ! Bounds of dimensions of node-centered data.
      integer:: ci1,ci2,cj1,cj2,ck1,ck2 ! Bounds of dimensions of cell-centered data.
      integer:: nnode,ncell             ! Number of nodes and cells.
      integer:: i,j,k,s                 ! Counters.
      integer:: Nx, Ny, Nz
      integer:: start

      Nx = size(orion%block(b)%mesh,2)-1; Ny = size(orion%block(b)%mesh,3)-1; Nz = size(orion%block(b)%mesh,4)-1
      ! initialize the zone dimensions
      call compute_dimensions(node=orion%tec%node,bc=orion%tec%bc,             &
                              Nx=Nx,Ny=Ny,Nz=Nz,gc=gc,                         &
                              ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2, &
                              ci1=ci1,ci2=ci2,cj1=cj1,cj2=cj2,ck1=ck1,ck2=ck2)
      nnode = (ni2-ni1+1)*(nj2-nj1+1)*(nk2-nk1+1)
      ncell = (ci2-ci1+1)*(cj2-cj1+1)*(ck2-ck1+1)
      ! writing the block data
      select case(orion%tec%format)
      case('binary')
#     if defined(TECIO)
        err = teczne142(trim(orion%block(b)%name)//tecendrec,         &
                        0,                                            &
                        ni2-ni1+1,                                    &
                        nj2-nj1+1,                                    &
                        nk2-nk1+1,                                    &
                        0,                                            &
                        0,                                            &
                        0,                                            &
                        time_,                                        &
                        0,                                            &
                        0,                                            &
                        1,                                            & !1=>block,0=>point
                        0,                                            &
                        0,                                            &
                        0,                                            &
                        0,                                            &
                        0,                                            &
                        tecnull(1:nvar),                              &
                        tecvarloc(1:nvar),                            &
                        tecnull(1:nvar),                              &
                        0)
        err=tec_dat(N=nnode,dat=orion%block(b)%mesh(1,ni1:ni2,nj1:nj2,nk1:nk2))
        if (ndir>=2) &
          err=tec_dat(N=nnode,dat=orion%block(b)%mesh(2,ni1:ni2,nj1:nj2,nk1:nk2))
        if (ndir==3) &
          err=tec_dat(N=nnode,dat=orion%block(b)%mesh(3,ni1:ni2,nj1:nj2,nk1:nk2))
        if (.not.meshonly) then
          start = 1
          if (orion%tec%node) start = 0
          if (ni2==0) ni2 = 1
          if (nj2==0) nj2 = 1
          if (nk2==0) nk2 = 1
          ! Force values if 2D
          if (ndir==2) then
            nk1 = 1-start; nk2 = 1
          endif
          do s=1,Nvar
            err=tec_dat(N=ncell,dat=orion%block(b)%vars(s,ni1+start:ni2,nj1+start:nj2,nk1+start:nk2))
          enddo
        endif
#     endif
      case('ascii')
        ! tecplot zone header
        teczoneheader = ' ZONE  T = '//trim(orion%block(b)%name)//            &
                        ', I='//trim(str(no_sign=.true.,n=ni2-ni1+1))//       &
                        ', J='//trim(str(no_sign=.true.,n=nj2-nj1+1))//       &
                        ', K='//trim(str(no_sign=.true.,n=nk2-nk1+1))//       &
                        ', DATAPACKING=BLOCK'//adjustl(trim(tecvarlocstr))
        
        teczoneheader = trim(teczoneheader)//', SOLUTIONTIME='//trim(str(no_sign=.true.,n=time_))
        if (.not. time_accurate) teczoneheader = trim(teczoneheader)//', STRANDID = 0'
        write(tecunit,'(A)',iostat=err)trim(teczoneheader)
        write(tecunit,FR_P,iostat=err)(((orion%block(b)%mesh(1,i,j,k),i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
        if (ndir>1) &
          write(tecunit,FR_P,iostat=err)(((orion%block(b)%mesh(2,i,j,k),i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
        if (ndir>2) &
          write(tecunit,FR_P,iostat=err)(((orion%block(b)%mesh(3,i,j,k),i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
        if (.not.meshonly) then
          start = 1
          if (orion%tec%node) start = 0
          if (ni2==0) ni2 = 1
          if (nj2==0) nj2 = 1
          if (nk2==0) nk2 = 1
          ! Force values if 2D
          if (ndir==2) then
            nk1 = 1-start; nk2 = 1
          endif
          do s=1,Nvar
            write(tecunit,FR_P,iostat=err)(((orion%block(b)%vars(s,i,j,k),i=ni1+start,ni2),j=nj1+start,nj2),k=nk1+start,nk2)
          enddo
        endif
      end select
    endfunction tec_blk_data

#   if defined(TECIO)
    ! Function interface for using "tecdat" function.
    function tec_dat(N,dat) result(err) 
      implicit none
      integer, intent(IN) :: N        ! Number of data to save.
      real(8), intent(IN) :: dat(1:N) ! Data to save.
      integer :: err                  ! Error trapping flag: 0 no errors, >0 error occurs.
      err = tecdatd142(N,dat)
    endfunction tec_dat
#   endif

  endfunction tec_write_structured_multiblock


  function tec_write_points_multivars(orion,varnames,filename,Nvars) result(err)
    implicit none
    type(orion_data), intent(in)              :: orion
    character(len=*), intent(in), optional    :: varnames
    character(len=*), intent(in)              :: filename !< File name of the output file.
    integer, intent(in), optional             :: Nvars    !< Input number of variables saved.
    integer :: err
    character(500)::          tecvarname          !< Variables name for tecplot header file.
    character(500)::          teczoneheader       !< Tecplot string of zone header.
    character(500)::          tecvarform          !< Format for variables for tecplot file.
    integer::                 tecunit             !< Free logic unit of tecplot file.
    integer::                 Nvar                !< Internal number of variables saved.
    integer::                 Nblocks             !< Number of blocks.
    integer::                 b                   !< Counter.
  
    ! Preliminary operations
    if (allocated(orion%block(1)%vars) .and. .not.present(Nvars)) then
      Nvar = size(orion%block(1)%vars,1)
    elseif (allocated(orion%block(1)%vars) .and. present(Nvars)) then
      Nvar = Nvars
    elseif (.not.allocated(orion%block(1)%vars)) then
      err = 1
      return
    endif
    Nblocks = size(orion%block)

    ! initializing tecplot variables
    call tec_init()

    ! initializing tecplot file
    open(newunit=tecunit,file=trim(filename))
    write(tecunit,'(A)',iostat=err)trim(tecvarname)
    ! writing data blocks
    do b=1,Nblocks
      err = tec_blk_data(b = b)
    enddo
    ! finalizing tecplot file
    close(tecunit)

  contains
    subroutine tec_init()
      implicit none
      integer  :: s ! Counter.
      ! header variables names
      tecvarname = ' VARIABLES ='
      if (present(varnames)) then
        tecvarname = trim(tecvarname)//' '//trim(varnames)
      else
        do s = 1, Nvar
          tecvarname = trim(tecvarname)//' "var'//trim(str(.true.,s))//'"'
        enddo
      endif
      ! variables output format
      write(tecvarform,'(A)')'('//trim(str(no_sign=.true.,n=Nvar))//'('//FR_P//',1X))'
    endsubroutine tec_init

    function tec_blk_data(b) result(err)
      implicit none
      integer, intent(IN):: b           ! Block number.
      integer:: err                     ! Error trapping flag: 0 no errors, >0 error occurs.
      integer:: i
      integer:: Nx

      Nx = size(orion%block(b)%mesh,2)
      ! tecplot zone header
      teczoneheader = ' ZONE  T = '//trim(orion%block(b)%name)//     &
                      ', I='//trim(str(no_sign=.true.,n=Nx))//       &
                      ', F=POINT'
      write(tecunit,'(A)',iostat=err)trim(teczoneheader)
      do i = 1, Nx
        write(tecunit,*,iostat=err) orion%block(b)%mesh(1,i,1,1), orion%block(b)%vars(1:nvar,i,1,1)
      enddo
    endfunction tec_blk_data

  endfunction tec_write_points_multivars


  !> Function for reading ORION block data from Tecplot file.
  function tec_read_structured_multiblock(orion,filename) result(err)
    implicit none
    type(orion_data), intent(inout)                        :: orion
    character(len=*), intent(in)                           :: filename
    integer :: err

    if (index(filename,'.dat')>0 .or. index(filename,'.tec')>0) then
      err = tec_read_ascii(orion,filename)
    elseif (index(filename,'.szplt')>0) then
#     if defined(TECIO)
      err = tec_read_szplt(orion,filename)
#     else
      stop "You can not read in binary format without compiling against TecIO"
#     endif
    endif

  end function tec_read_structured_multiblock


  !> Function for reading ORION block data from Tecplot ascii file.
  function tec_read_ascii(orion,filename) result(err)
    use, intrinsic :: iso_fortran_env, only : iostat_end
    use strings, only: getvals, parse
    implicit none
    type(orion_data), intent(inout)                        :: orion
    character(len=*), intent(in)                           :: filename
    real(R8P) :: dummy_float, solutiontime
    logical :: meshonly
    integer :: err, end
    integer :: tecunit, ios, ios_prev, ios2
    integer :: i, j, k, d, b, b2, ck1, ck2
    integer, allocatable :: Ni(:), Nj(:), Nk(:), nskip(:)
    integer :: Nblocks, nlines, nvar, ndir
    character(500) :: line
    character(100) :: args(20), subargs(2)

    meshonly = .false.

    ! Open file
    open(newunit=tecunit,file=trim(filename),status='old',action='read',iostat=err)
    if (err/=0) return

    ios = 0
    do while(ios==0)
      read(tecunit,'(A)',iostat=ios) line
      call read_variables(line, orion%varnames)
      if (allocated(orion%varnames)) exit
    enddo

    rewind(tecunit)
    
    ! Count blocks and allocate Ni, Nj, Nk, nskip
    ios = 0; Nblocks = 0; nlines = -1
    do while(ios==0)
      read(tecunit,'(A)',iostat=ios) line
      nlines = nlines+1
      if (index(line,"ZONE")>0 .and. index(line,"ZONETYPE")==0) Nblocks = Nblocks+1
      if (index(line,"Zone")>0) Nblocks = Nblocks+1
    enddo
    allocate(Ni(Nblocks));allocate(Nj(Nblocks));allocate(Nk(Nblocks));allocate(Nskip(Nblocks));

    rewind(tecunit)

    ! Read blocks dimensions
    ios = 0; b = 1; ios2 = 0; nskip = 0; ios_prev = 0; b2 = 1
    do while(ios==0)
      read(tecunit,'(A)',iostat=ios) line
      if (index(line,'I=')>0) then
        call parse(line,',',args)
        do i = 1, 6!size(args)
          if (index(args(i),'I=')>0) then
            call parse(args(i),'=',subargs)
            read(subargs(2),'(I5)') Ni(b)
          endif
          if (index(args(i),'J=')>0) then
            call parse(args(i),'=',subargs)
            read(subargs(2),'(I5)') Nj(b)
          endif
          if (index(args(i),'K=')>0) then
            call parse(args(i),'=',subargs)
            read(subargs(2),'(I5)') Nk(b)
          endif
        enddo
        b = b+1
      endif
       ! Count not-floating lines
      read(line,*,iostat=ios2) dummy_float
      if ( (ios2==0 .and. index(line,'DATA')>0) .or. ios2/=0 ) then
        nskip(b2) = nskip(b2)+1
        ios2 = 1
      endif
      if (ios2==0 .and. ios_prev/=0) b2 = b2+1
      ios_prev = ios2
    enddo
    ! Allocate data
    allocate(orion%block(1:Nblocks))
    do b = 1, Nblocks
      orion%block(b)%Ni = Ni(b)-1
      orion%block(b)%Nj = Nj(b)-1
      orion%block(b)%Nk = Nk(b)-1
    enddo
    rewind(tecunit)

    ! Find solutiontime
    ios = 0; solutiontime = -10._R8P;
    do j=1,10
      read(tecunit,'(A)',iostat=ios) line
      if (index(line,'SOLUTIONTIME=')>0) then
        call parse(line,',',args)
        do i = 1, size(args)
          if (index(args(i),'SOLUTIONTIME=')>0) then
            call parse(args(i),'=',subargs)
            read(subargs(2),FR_P) solutiontime
          endif
        enddo
      endif
      if (solutiontime>0.0_R8P) exit
    enddo
    orion%solutiontime = solutiontime
    rewind(tecunit)

    if (Nk(1)==1) then
      ndir = 2
    else
      ndir = 3
    endif

    ! Count variables
    Nvar = nlines-sum(nskip)
    do b = 1, Nblocks
      nvar = nvar-ndir*((orion%block(b)%Ni+1)*(orion%block(b)%Nj+1)*(orion%block(b)%Nk+1))
    enddo
    nvar = nvar/(sum(orion%block(:)%Ni*orion%block(:)%Nj*max(orion%block(:)%Nk,1)))
    if (nvar==0) meshonly = .true.

    ! Read all
    end = 0; if (orion%tec%node) end = 1
    do b = 1, Nblocks
      allocate(orion%block(b)%mesh(1:ndir,0:orion%block(b)%Ni,0:orion%block(b)%Nj,0:orion%block(b)%Nk))
      if (.not.meshonly) allocate(orion%block(b)%vars(1:nvar,1:orion%block(b)%Ni,1:orion%block(b)%Nj,1:max(orion%block(b)%Nk,1)))
      call skip(tecunit,nskip(b))
      do d = 1, ndir
        do k = 0, orion%block(b)%Nk; do j = 0, orion%block(b)%Nj; do i = 0, orion%block(b)%Ni
              read(tecunit,*,iostat=err) orion%block(b)%mesh(d,i,j,k)
        enddo; enddo; enddo
      enddo
      if (ndir==2) then
        ck1 = 1; ck2 = 1
      else
        ck1 = 1; ck2 = orion%block(b)%Nk + end
      endif
      do d = 1, nvar
        do k = ck1, ck2; do j = 1, orion%block(b)%Nj+end; do i = 1, orion%block(b)%Ni+end
              read(tecunit,*,iostat=err) orion%block(b)%vars(d,i,j,k)
        enddo; enddo; enddo
      enddo
    enddo

    close(tecunit)

  end function tec_read_ascii


  !> Function for reading ORION block data from Tecplot file.
  function tec_read_points_multivars(orion,nvar,filename) result(err)
    use, intrinsic :: iso_fortran_env, only : iostat_end
    use strings, only: getvals, parse
    implicit none
    type(orion_data), intent(inout)  :: orion
    character(len=*), intent(in)     :: filename
    integer, intent(in)              :: nvar
    real(R8P) :: dummy_float
    integer :: err
    integer :: tecunit, ios, ios_prev
    integer :: i, j, k, b
    integer :: Nzones, nlines
    integer, allocatable :: nskip(:)
    character(500) :: line
    character(100) :: args(20), subargs(2)

    ! Open file
    open(newunit=tecunit,file=trim(filename),iostat=err,action='read',status='old')
    if (err/=0) return
    
    ! Count blocks and allocate data
    ios = 0; Nzones = 0; nlines = -1
    do while(ios==0)
      read(tecunit,'(A)',iostat=ios) line
      nlines = nlines+1
      if (index(line,"ZONE")>0 .and. index(line,"ZONETYPE")==0) then
        Nzones = Nzones+1
      elseif (index(line,"Zone")>0) then
        Nzones = Nzones+1
      elseif (index(line,"ZONE T")>0) then
        Nzones = Nzones+1
      endif
    enddo
    allocate(orion%block(1:Nzones))
    rewind(tecunit)

    ! Read blocks size
    ios = 0; b = 0
    do
      do while (index(line,'I=')==0 .and. ios/=iostat_end)
        read(tecunit,'(A)',iostat=ios) line
      enddo
      if (ios==iostat_end) exit
      b = b+1
      call parse(line,',',args)
      do i = 1, 2
        if (index(args(i),'I=')>0) then
          call parse(args(i),'=',subargs)
          read(subargs(2),'(I8)') orion%block(b)%Ni
          orion%block(b)%Ni = orion%block(b)%Ni
          orion%block(b)%Nj = 1
          orion%block(b)%Nk = 1
        endif
      enddo
      line = 'here we go'
    enddo
    rewind(tecunit)

    ! Count not-floating lines
    allocate(nskip(Nzones))
    nskip = 0; ios = 0; b = 1; ios_prev = 0
    do
      read(tecunit,'(A)',iostat=ios) line
      if (ios==iostat_end) exit
      read(line,*,iostat=ios) dummy_float
      if ((ios==0 .and. index(line,'DATA')>0) .or. ios/=0) then
        nskip(b) = nskip(b)+1
        ios = 1
      endif
      if (ios==0 .and. ios_prev/=0) b = b+1
      ios_prev = ios
    enddo
    rewind(tecunit)

    ! Read all
    do b = 1, Nzones
      allocate(orion%block(b)%mesh(1:1,1:orion%block(b)%Ni,1:orion%block(b)%Nj,1:orion%block(b)%Nk))
      allocate(orion%block(b)%vars(1:nvar,1:orion%block(b)%Ni,1:orion%block(b)%Nj,1:orion%block(b)%Nk))
      call skip(tecunit,nskip(b))
      do k = 1, orion%block(b)%Nk; do j = 1, orion%block(b)%Nj; do i = 1, orion%block(b)%Ni
            read(tecunit,*,iostat=err) orion%block(b)%mesh(1,i,j,k), orion%block(b)%vars(1:nvar,i,j,k)
      enddo; enddo; enddo
    enddo

    close(tecunit)

  end function tec_read_points_multivars


  subroutine read_variables(line,variables)
    implicit none
    character(len=*), intent(inout) :: line
    character(len=32), allocatable, intent(out) :: variables(:)
    character(len=32) :: variables_(150)
    integer :: nvar, start, end_pos

    nvar = 0

    ! Find the position of "VARIABLES"
    if (index(line, 'VARIABLES') > 0) then
        ! Remove the "VARIABLES =" part from the line
        line = trim(adjustl(line(index(line, '=')+1:)))

        ! Loop to extract each variable name
        do while (len_trim(line) > 0)
            ! Find the start and end positions of the variable name
            start = index(line, '"') + 1
            end_pos = index(line(start+1:), '"') + start

            ! Store the variable name in the array
            nvar = nvar + 1
            variables_(nvar) = trim(line(start:end_pos-1))

            ! Remove the extracted variable from the line
            if (index(line(end_pos+1:), ',') > 0 .or. index(line(end_pos+1:), ' ') > 0) then
                line = trim(adjustl(line(end_pos+2:)))
            else
                line = ''
            end if
        end do
    else
      return
    end if
    
    allocate(character(32)::variables(1:nvar))
    variables = variables_(1:nvar)

  end subroutine read_variables


  ! Convenience to extract a returned C string to a FOrtran string
  subroutine copyCharArrayToString(charArray, length, string)
    use iso_c_binding, only : C_NULL_CHAR, c_ptr, c_f_pointer
    implicit none
    type(c_ptr) :: charArray
    integer length
    character(*) string

    character, pointer :: charPointer(:)
    integer i

    call c_f_pointer(charArray, charPointer, [length])

    string = ' '
    do i = 1, length
        string(i:i) = charPointer(i)
    enddo
    string(length+1:length+1) = C_NULL_CHAR

  end

# if defined(TECIO)
  function tec_read_szplt(orion,filename) result(i)
    use iso_c_binding
    implicit none
    type(orion_data), intent(inout)                        :: orion
    character(len=*), intent(in)                           :: filename

    integer i, j, k, cnt
    character(256) inputFileName
    character(256) dataSetTitle, zoneTitle
    character(1024) varNames
    character, pointer :: stringPtr(:)
    integer nameLen, strLen
    integer(c_int8_t), allocatable :: int8Values(:)
    integer(c_int16_t), allocatable :: int16Values(:)
    integer(c_int32_t) :: numVars, var
    integer(c_int32_t) :: fileType
    integer(c_int32_t) :: vstart
    integer(c_int32_t) :: inputZone, numZones
    integer(c_int32_t) :: zoneType
    integer(c_int32_t) :: strandID
    integer(c_int32_t), allocatable :: int32Values(:)
    integer(c_int64_t) :: numValues
    integer(c_int64_t) :: iMax, jMax, kMax, ndir
    integer(c_int32_t), allocatable :: varTypes(:)
    integer(c_int32_t), allocatable :: passiveVarList(:)
    integer(c_int32_t), allocatable :: valueLocation(:)
    integer(c_int32_t), allocatable :: shareVarFromZone(:)
    real(c_float), allocatable :: floatValues(:)
    real(c_double) :: solutionTime
    real(c_double), allocatable :: doubleValues(:)
    type(c_ptr) :: inputFileHandle = C_NULL_PTR
    type(c_ptr) :: stringCPtr = C_NULL_PTR
    logical :: onlyNode

    inputFileName = trim(filename) // C_NULL_CHAR

    ! Open the input file for reading
    i = tecFileReaderOpen(inputFileName, inputFileHandle)

    if (i/=0) return

    ! Read info about the data set
    i = tecDataSetGetTitle(inputFileHandle, stringCPtr)
    call copyCharArrayToString(stringCPtr, &
        tecStringLength(stringCPtr), dataSetTitle)
    call tecStringFree(stringCPtr)
    i = tecDataSetGetNumVars(inputFileHandle, numVars)

    strLen = 0
    do var = 1, numVars
        i = tecVarGetName(inputFileHandle, var, stringCPtr)
        nameLen = tecStringLength(stringCPtr)
        call c_f_pointer(stringCPtr, stringPtr, [nameLen])
        if (var .gt. 1) then
            strLen = strLen + 1
            varNames(strLen : strLen) = ','
        endif
        do j = 1, nameLen
            varNames(strLen + j : strLen + j) = stringPtr(j)
        enddo
        strLen = strLen + nameLen
        call tecStringFree(stringCPtr)
    enddo
    varNames(strLen + 1 : strlen + 1) = C_NULL_CHAR

    i = tecFileGetType(inputFileHandle, fileType)
    i = tecDataSetGetNumZones(inputFileHandle, numZones)

    allocate(orion%block(1:numZones))

    ! Zones
    do inputZone = 1, numZones
        i = tecZoneGetType(inputFileHandle, inputZone, zoneType)
        if (zoneType == 6 .or. zoneType == 7) &
            stop "Unsupported inputZone type."
          
        ! Retrieve info about the inputZone
        i = tecZoneGetTitle(inputFileHandle, inputZone, stringCPtr)
        call copyCharArrayToString(stringCPtr, &
            tecStringLength(stringCPtr), zoneTitle)
        call tecStringFree(stringCPtr)

        orion%block(inputZone)%name = zoneTitle
        
        i = tecZoneGetIJK(inputFileHandle, inputZone, &
            iMax, jMax, kMax)

        orion%block(inputZone)%Ni = iMax-1
        orion%block(inputZone)%Nj = jMax-1
        orion%block(inputZone)%Nk = kMax-1

        if (jMax>1 .and. kMax>1) then
          ndir = 3
        elseif (jMax>1 .and. kMax==1) then
          ndir = 2
        elseif (jMax==1 .and. kMax==1) then
          ndir = 1
        endif

        allocate(valueLocation(numVars))
        do var = 1, numVars
            i = tecZoneVarGetValueLocation(inputFileHandle, inputZone, &
                var, valueLocation(var))
        enddo

        onlyNode = .false.
        do var = ndir+1, numVars
          if (valueLocation(var)==1) onlyNode = .true.
        enddo

        allocate(orion%block(inputZone)%mesh(1:ndir,0:iMax-1,0:jMax-1,0:kMax-1))
        if (onlyNode) then
          vstart = 0
          allocate(orion%block(inputZone)%vars(1:numVars-ndir,vstart:max(vstart,iMax-1),vstart:max(vstart,jMax-1),vstart:max(vstart,kMax-1)))
        else
          vstart = 1
          allocate(orion%block(inputZone)%vars(1:numVars-ndir,vstart:max(vstart,iMax-1),vstart:max(vstart,jMax-1),vstart:max(vstart,kMax-1)))
        endif
        
        allocate(varTypes(numVars))
        allocate(passiveVarList(numVars))
        allocate(shareVarFromZone(numVars))
        do var = 1, numVars
            i = tecZoneVarGetType(inputFileHandle, inputZone, &
                var, varTypes(var))
            i = tecZoneVarIsPassive(inputFileHandle, inputZone, &
                var, passiveVarList(var))
            i = tecZoneVarGetSharedZone(inputFileHandle, inputZone, &
                var, shareVarFromZone(var))
        enddo

        ! i = tecZoneConnectivityGetSharedZone(inputFileHandle, &
        !     inputZone, shareConnectivityFromZone)
        ! i = tecZoneFaceNbrGetMode(inputFileHandle, inputZone, &
        !     faceNeighborMode)
        ! if (faceNeighborMode > 4) faceNeighborMode = 1
        ! i = tecZoneFaceNbrGetNumConnections(inputfileHandle, &
        !     inputZone, numFaceConnections)

        i = tecZoneGetSolutionTime(inputFileHandle, inputZone, &
            solutionTime)
        orion%solutiontime = solutionTime      ! fixed 2026-09-06: was read and dropped -> szplt restarts began at t=0
        i = tecZoneGetStrandID(inputFileHandle, inputZone, strandID)
        ! if (solutionTime /= 0.0 .or. strandID /= 0) &
        !     i = tecZoneSetUnsteadyOptions(outputFileHandle, &
        !         outputZone, solutionTime, strandID)

        ! Read and write inputZone data
        do var = 1, numVars
            if (passiveVarList(var) == 0 .and. &
                  shareVarFromZone(var) == 0) then            
                i = tecZoneVarGetNumValues(inputFileHandle, &
                          inputZone, var, numValues)
                select case (varTypes(var))
                case (1) ! float
                    allocate(floatValues(numValues))
                    i = tecZoneVarGetFloatValues(inputFileHandle, &
                        inputZone, var, 1_c_int64_t, numValues, &
                        floatValues)
                    if (valueLocation(var)==1 .and. var<=ndir) then
                      cnt = 1
                      do k = 0, kMax-1; do j = 0, jMax-1; do i = 0, iMax-1
                            orion%block(inputZone)%mesh(var,i,j,k) = floatValues(cnt)
                            cnt = cnt + 1
                      enddo; enddo; enddo
                    else
                      cnt = 1
                      do k = vstart, max(vstart,kMax-1); do j = vstart, max(vstart,jMax-1); do i = vstart, max(vstart,iMax-1)
                            orion%block(inputZone)%vars(var-ndir,i,j,k) = floatValues(cnt)
                            cnt = cnt + 1
                      enddo; enddo; enddo
                    endif
                    deallocate(floatValues)
                case (2) ! double
                    allocate(doubleValues(numValues))
                    i = tecZoneVarGetDoubleValues(inputFileHandle, &
                        inputZone, var, 1_c_int64_t, numValues, &
                        doubleValues)
                    if (valueLocation(var)==1 .and. var<=ndir) then   ! fixed 2026-09-06: was floatValues (unallocated) and no var<=ndir
                      cnt = 1
                      do k = 0, kMax-1; do j = 0, jMax-1; do i = 0, iMax-1
                            orion%block(inputZone)%mesh(var,i,j,k) = doubleValues(cnt)
                            cnt = cnt + 1
                      enddo; enddo; enddo
                    else
                      cnt = 1
                      do k = vstart, max(vstart,kMax-1); do j = vstart, max(vstart,jMax-1); do i = vstart, max(vstart,iMax-1)
                            orion%block(inputZone)%vars(var-ndir,i,j,k) = doubleValues(cnt)
                            cnt = cnt + 1
                      enddo; enddo; enddo
                    endif
                    deallocate(doubleValues)
                case (3) ! int32_t
                    allocate(int32Values(numValues))
                    i = tecZoneVarGetInt32Values(inputFileHandle, &
                        inputZone, var, 1_c_int64_t, numValues, &
                        int32Values)
                    if (valueLocation(var)==1) then
                      cnt = 1
                      do k = 0, kMax-1; do j = 0, jMax-1; do i = 0, iMax-1
                            orion%block(inputZone)%mesh(var,i,j,k) = floatValues(cnt)
                            cnt = cnt + 1
                      enddo; enddo; enddo
                    else
                      cnt = 1
                       do k = vstart, max(vstart,kMax-1); do j = vstart, max(vstart,jMax-1); do i = vstart, max(vstart,iMax-1)
                            orion%block(inputZone)%vars(var-ndir,i,j,k) = floatValues(cnt)
                            cnt = cnt + 1
                      enddo; enddo; enddo
                    endif
                    deallocate(int32Values)
                case (4) ! int16_t
                    allocate(int16Values(numValues))
                    i = tecZoneVarGetInt16Values(inputFileHandle, &
                        inputZone, var, 1_c_int64_t, numValues, &
                        int16Values)
                    if (valueLocation(var)==1) then
                      cnt = 1
                      do k = 0, kMax-1; do j = 0, jMax-1; do i = 0, iMax-1
                            orion%block(inputZone)%mesh(var,i,j,k) = floatValues(cnt)
                            cnt = cnt + 1
                      enddo; enddo; enddo
                    else
                      cnt = 1
                       do k = vstart, max(vstart,kMax-1); do j = vstart, max(vstart,jMax-1); do i = vstart, max(vstart,iMax-1)
                            orion%block(inputZone)%vars(var-ndir,i,j,k) = floatValues(cnt)
                            cnt = cnt + 1
                      enddo; enddo; enddo
                    endif
                    deallocate(int16Values)
                case (5) ! uint8_t
                    allocate(int8Values(numValues))
                    i = tecZoneVarGetUInt8Values(inputFileHandle, &
                        inputZone, var, 1_c_int64_t, numValues, &
                        int8Values)
                    if (valueLocation(var)==1) then
                      cnt = 1
                      do k = 0, max(1,kMax-1); do j = 0, jMax-1; do i = 0, iMax-1
                            orion%block(inputZone)%mesh(var,i,j,k) = floatValues(cnt)
                            cnt = cnt + 1
                      enddo; enddo; enddo
                    else
                      cnt = 1
                       do k = vstart, max(vstart,kMax-1); do j = vstart, max(vstart,jMax-1); do i = vstart, max(vstart,iMax-1)
                            orion%block(inputZone)%vars(var-ndir,i,j,k) = floatValues(cnt)
                            cnt = cnt + 1
                      enddo; enddo; enddo
                    endif
                    deallocate(int8Values)
                endselect
            endif
        enddo
        
        deallocate(varTypes)
        deallocate(passiveVarList)
        deallocate(valueLocation)
        deallocate(shareVarFromZone)
                  
    enddo ! inputZone loop

    ! Close old and new files
    i = tecFileReaderClose(inputFileHandle)

  end function tec_read_szplt
# endif


  subroutine skip(u,n)
    implicit none
    integer, intent(in) :: u, n
    integer :: i
    do i = 1, n; read(u,*); enddo
  end subroutine skip

  end module Lib_Tecplot
