module functions
  use IR_Precision
  use Lib_ORION_data
  use Lib_VTK
  use Lib_Tecplot
  use Lib_PLOT3D
  implicit none


contains


  pure function extract_path(inpath) result(outpath)
    implicit none
    integer :: index_start
    character(len=256), intent(in)    :: inpath
    character(len=256)                :: outpath

    ! Find the last occurrence of '/' in the string
    index_start = len_trim(inpath)
    do while (index_start > 0)
      if (inpath(index_start:index_start) == '/') exit
      index_start = index_start - 1
    end do
    ! Extract the folder path
    outpath = inpath(1:index_start)
    if (outpath=='') outpath = './'

  end function extract_path


  subroutine read_vts(orion,files)
    implicit none
    type(orion_data), intent(inout) :: orion
    character(len=*), intent(in)    :: files(:)
    real(R8P), allocatable          :: x(:),y(:),z(:) ! Input geo arrays
    real(R8P), allocatable          :: v(:)           ! Input var arrays
    character(256), allocatable     :: varnames(:)
    integer(I4P)                    :: Nblocks,nx1,nx2,ny1,ny2,nz1,nz2,nn,nc
    integer(I4P)                    :: b, i, j, k, n, p, err, start

    Nblocks = size(files)

    allocate(orion%block(1:Nblocks))

    ! Read VTS file
    do b = 1, Nblocks
      
      call read_variables_name(files(b),varnames,orion%tec%node)
      ! Geometry
      err = VTK_INI_XML_READ(input_format=trim(orion%vtk%format),filename=trim(files(b)),mesh_topology='StructuredGrid',&
                              nx1=nx1,nx2=nx2,ny1=ny1,ny2=ny2,nz1=nz1,nz2=nz2)
      err = VTK_GEO_XML_READ(nx1=nx1,nx2=nx2,ny1=ny1,ny2=ny2,nz1=nz1,nz2=nz2,NN=nn,X=x,Y=y,Z=z)
      allocate(orion%block(b)%mesh(1:3,nx1:nx2,ny1:ny2,nz1:nz2))
      n = 0
      do k = nz1, nz2; do j = ny1, ny2; do i = nx1, nx2
            n = n + 1
            orion%block(b)%mesh(1,i,j,k) = x(n)
            orion%block(b)%mesh(2,i,j,k) = y(n)
            orion%block(b)%mesh(3,i,j,k) = z(n)
      enddo; enddo; enddo
      ! Variables field
      do p = 1, size(varnames)-1
        if (allocated(v)) deallocate(v)
        if (orion%tec%node) then
          err = VTK_VAR_XML_READ(var_location='node', varname=trim(varnames(p+1)), NC_NN=nn, NCOMP=nc, var=v)
          write(*,'(A)') ' - Data location : cell nodes'
          start = 0
        else
          err = VTK_VAR_XML_READ(var_location='cell', varname=trim(varnames(p+1)), NC_NN=nn, NCOMP=nc, var=v)
          write(*,'(A)') ' - Data location : cell centers'
          start = 1
        endif
        allocate(orion%block(b)%vars(1:size(varnames)-1,nx1+start:nx2,ny1+start:ny2,nz1+start:nz2))
        n = 0
        do k = nz1+start, nz2; do j = ny1+start, ny2; do i = nx1+start, nx2
              n = n + 1
              orion%block(b)%vars(p,i,j,k) = v(n)
        enddo; enddo; enddo
      enddo
    enddo
    err = VTK_END_XML_READ()

  end subroutine read_vts


  subroutine read_vtm(orion,file)
    implicit none
    type(orion_data), intent(inout) :: orion
    character(len=*), intent(in)    :: file
    integer :: E_IO

    E_IO = vtk_read_structured_multiblock(orion=orion,vtmpath=trim(file),vtspath='')
    if (E_IO/=0) then
      write(*,'(A)')' ERROR read VTK files'
      stop
    endif

  end subroutine read_vtm


  subroutine write_vtm(orion,file,varname_scalar)
    use strings, only: parse
    implicit none
    type(orion_data), intent(inout) :: orion
    character(len=*), intent(in)    :: file
    character(len=*), intent(inout) :: varname_scalar
    integer :: E_IO
    character(32)  :: name
    character(256) :: path

    name = file(1:len_trim(file)-4)
    path = extract_path(file)
    E_IO = vtk_write_structured_multiblock(orion=orion,vtspath=trim(path),vtmpath=trim(name),varnames=varname_scalar)
    if (E_IO/=0) then
      write(*,'(A)')' ERROR write VTK files'
      stop
    endif

  end subroutine write_vtm


  subroutine read_dat(orion,file)
    implicit none
    type(orion_data), intent(inout) :: orion
    character(len=*), intent(in)    :: file
    integer :: E_IO, ntry

    do ntry = 1, 313
      E_IO = tec_read_points_multivars(orion=orion,filename=trim(file),nvar=2)
      if (E_IO/=0) then
        if (allocated(orion%block)) deallocate(orion%block)
      else
        exit
      endif
    enddo

    if (E_IO/=0) then
      write(*,'(A)')' ERROR read Tecplot file'
      stop
    endif

  end subroutine read_dat


  subroutine read_tec(orion,file)
    implicit none
    type(orion_data), intent(inout) :: orion
    character(len=*), intent(in)    :: file
    integer :: E_IO

    ! Dispatch by extension: a single shared %tec%format cannot describe both an
    ! ascii .tec input and a binary .szplt output in the same run.
    if (index(trim(file),'.szplt')>0 .or. index(trim(file),'.plt')>0) then
      orion%tec%format = 'binary'
    else
      orion%tec%format = 'ascii'
    endif
    E_IO = tec_read_structured_multiblock(orion=orion,filename=trim(file))
    if (E_IO/=0) then
      write(*,'(A)')' ERROR read Tecplot file'
      stop
    endif

  end subroutine read_tec


  subroutine write_tec(orion,file,varname_scalar)
    implicit none
    type(orion_data), intent(inout) :: orion
    character(len=*), intent(in)    :: file,varname_scalar
    integer :: E_IO

    if (index(trim(file),'.szplt')>0 .or. index(trim(file),'.plt')>0) then
      orion%tec%format = 'binary'
    else
      orion%tec%format = 'ascii'
    endif
    E_IO = tec_write_structured_multiblock(orion=orion,varnames=varname_scalar,filename=trim(file))
    if (E_IO/=0) then
      write(*,'(A)')' ERROR write Tecplot file'
      stop
    endif

  end subroutine write_tec


  subroutine write_dat(orion,file,varname_scalar)
    implicit none
    type(orion_data), intent(inout) :: orion
    character(len=*), intent(in)    :: file,varname_scalar
    integer :: E_IO

    E_IO = tec_write_points_multivars(orion=orion,varnames=varname_scalar,filename=trim(file))
    if (E_IO/=0) then
      write(*,'(A)')' ERROR write Tecplot file'
      stop
    endif

  end subroutine write_dat


  subroutine read_p3d(orion,file)
    implicit none
    type(orion_data), intent(inout) :: orion
    character(len=*), intent(in)    :: file
    integer :: E_IO

    E_IO = p3d_read_multiblock(orion=orion,filename=trim(file))
    if (E_IO/=0) then
      write(*,'(A)')' ERROR write PLOT3D file'
      stop
    endif

  end subroutine read_p3d


  subroutine write_p3d(orion,file)
    implicit none
    type(orion_data), intent(inout) :: orion
    character(len=*), intent(in)    :: file
    integer :: E_IO

    E_IO = p3d_write_multiblock(orion=orion,filename=trim(file))
    if (E_IO/=0) then
      write(*,'(A)')' ERROR write PLOT3D file'
      stop
    endif

  end subroutine write_p3d

end module functions



program main
  use functions
  use Lib_ORION_data
  implicit none
  type(orion_data)             :: data
  integer                      :: b, st
  character(256)               :: infile, outfile, varname_scalar

  data%tec%format = 'binary'
  data%vtk%format = 'binary'
  data%p3d%format = 'ascii'

  call command_line_argument()

  write(*,*)
  write(*,'(2A)')' - Input file  : ',trim(infile)
  write(*,'(2A)')' - Output file : ',trim(outfile)
  write(*,*)

  if (index(trim(infile),'.dat')>0) call read_dat(data,infile)
  if (index(trim(infile),'.tec')>0) call read_tec(data,infile)
  if (index(trim(infile),'.szplt')>0) call read_tec(data,infile)
  if (index(trim(infile),'.p3d')>0) call read_p3d(data,infile)
  if (index(trim(infile),'.vtm')>0) call read_vtm(data,infile)

  varname_scalar = ''
  st = 1
  if (index(trim(outfile),'.dat')>0) st = 0
  do b = st, size(data%block(1)%vars,1)
    varname_scalar = trim(varname_scalar)//' '//'var'//trim(str(.true.,b))
  enddo
  do b = 1, size(data%block)
    data%block(b)%name = 'Block'//trim(str(.true.,b))
  enddo

  if (index(trim(outfile),'.dat')>0) call write_dat(data,outfile,varname_scalar)
  if (index(trim(outfile),'.tec')>0) call write_tec(data,outfile,varname_scalar)
  if (index(trim(outfile),'.szplt')>0) call write_tec(data,outfile,varname_scalar)
  if (index(trim(outfile),'.p3d')>0) call write_p3d(data,outfile)
  if (index(trim(outfile),'.vtm')>0) call write_vtm(data,outfile,varname_scalar)

  write(*,*)
  write(*,'(A)')' Done!'
  

contains


  subroutine command_line_argument()
    implicit none
    character(99):: arg
    integer :: arg_count, i

    ! Get the number of command-line arguments
    arg_count = COMMAND_ARGUMENT_COUNT()

    if (arg_count==0) then
      call print_help
      stop
    endif

    ! Loop through each command-line argument
    do i = 1, arg_count
      ! Get the i-th command-line argument
      call GET_COMMAND_ARGUMENT(i, arg)

      ! Check for different command-line options
      if (arg == '-h' .or. arg == '--help') then
        ! Display help information
        call print_help()
        stop

      elseif (index(arg, '--out-file=') > 0) then
        outfile = trim(adjustl(arg(12:)))

      elseif (index(arg, '--in-file=') > 0) then
        infile = trim(adjustl(arg(11:)))

      elseif (index(arg, '--out-format=') > 0) then
        if (index(outfile,'.p3d')>0) data%p3d%format = trim(adjustl(arg(14:)))
        if (index(outfile,'.vtm')>0) data%vtk%format = trim(adjustl(arg(14:)))

      elseif (index(arg, '--in-format=') > 0) then
        if (index(infile,'.p3d')>0) data%p3d%format = trim(adjustl(arg(14:)))
        if (index(infile,'.vtm')>0) data%vtk%format = trim(adjustl(arg(14:)))

      end if

      if (index(infile,'.dat')>0) data%tec%format = 'ascii'
      if (index(infile,'.tec')>0) data%tec%format = 'ascii'
      if (index(infile,'.szplt')>0) data%tec%format = 'binary'
      if (index(outfile,'.dat')>0) data%tec%format = 'ascii'
      if (index(outfile,'.tec')>0) data%tec%format = 'ascii'
      if (index(outfile,'.szplt')>0) data%tec%format = 'binary'

    end do

  end subroutine


  subroutine print_help
    implicit none
    write(*,*) "File format converter. Supported formats: PLOT3D, Tecplot, VTK"
    write(*,*)
    write(*,*) "Usage: ORION [input]"
    write(*,*) "Input:"
    write(*,*) "  -h, --help              Display this help message"
    write(*,*) "  --in-file=<file>        Choose the input file path"
    write(*,*) "  --in-format=<format>    Set the input format (e.g., binary,ascii)"
    write(*,*) "  --out-file=<file>       Choose the output file path"
    write(*,*) "  --out-format=<format>   Set the output format (e.g., binary,ascii)"
  end subroutine print_help

endprogram main
