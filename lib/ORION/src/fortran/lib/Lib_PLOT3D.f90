module Lib_PLOT3D
  use, intrinsic:: ISO_FORTRAN_ENV, only: stdout => OUTPUT_UNIT, stderr => ERROR_UNIT
  use IR_Precision
  use Lib_ORION_data
  implicit none
  private

  public:: p3d_read_multiblock
  public:: p3d_write_multiblock


contains


  function p3d_write_multiblock(orion,filename) result(err)
    implicit none
    type(orion_data), intent(in)              :: orion
    character(len=*), intent(in)              :: filename
    integer :: err
    integer :: unit, b, d, k, j, i, Nblocks, ndir
    
    Nblocks = size(orion%block)

    ! initializing file
    select case(orion%p3d%format)
    case('binary')
      stop "You can not write in binary P3D format"
    case('ascii')
      if (index(filename,'.')==0) then
        open(newunit=unit,file=trim(filename)//".p3d",iostat=err)
      else
        open(newunit=unit,file=trim(filename),iostat=err)
      endif
    end select

    if (orion%block(1)%Nk==0) then
      ndir = 2
    else
      ndir = 3
    endif

    ! writing data blocks
    write(unit,*,iostat=err) Nblocks
    
    do b = 1, Nblocks
      if (orion%block(b)%Nk==0) then
        write(unit,*,iostat=err) orion%block(b)%Ni+1, orion%block(b)%Nj+1
      else
        write(unit,*,iostat=err) orion%block(b)%Ni+1, orion%block(b)%Nj+1, orion%block(b)%Nk+1
      endif
    enddo

    do b = 1, Nblocks
      do d = 1, ndir
        do k = 0, orion%block(b)%Nk
          do j = 0, orion%block(b)%Nj
            do i = 0, orion%block(b)%Ni
              write(unit,*,iostat=err) orion%block(b)%mesh(d,i,j,k)
            enddo
          enddo
        enddo
      enddo
    enddo

    close(unit)

  endfunction p3d_write_multiblock


  function p3d_read_multiblock(orion,filename) result(err)
    use, intrinsic :: iso_fortran_env, only : iostat_end
    use strings, only: getvals, parse
    implicit none
    type(orion_data), intent(inout)              :: orion
    character(len=*), intent(in)                 :: filename
    integer :: err, dum
    integer :: unit, b, d, k, j, i, Nblocks, ndir
    character(len=256) :: line

    open(newunit=unit,file=trim(filename),iostat=err,form='formatted',status='old')
    if (err /= 0) return

    read(unit,*) Nblocks
    allocate(orion%block(1:Nblocks))

    ! Check ndir
    ndir = 0
    read(unit,*) line
    do i = 1, 3
      read(line, *, iostat=err) (dum, ndir=1,ndir)
        if (err /= 0) exit
        ndir = ndir + 1
    end do
    ndir = ndir - 1
    rewind(unit)
    read(unit,*)

    if (ndir==2) then
      orion%block(:)%Nk = 0
      do b = 1, Nblocks
        read(unit,*) orion%block(b)%Ni, orion%block(b)%Nj
        orion%block(b)%Ni = orion%block(b)%Ni-1
        orion%block(b)%Nj = orion%block(b)%Nj-1
      enddo
    else
      do b = 1, Nblocks
        read(unit,*) orion%block(b)%Ni, orion%block(b)%Nj, orion%block(b)%Nk
        orion%block(b)%Ni = orion%block(b)%Ni-1
        orion%block(b)%Nj = orion%block(b)%Nj-1
        orion%block(b)%Nk = orion%block(b)%Nk-1
      enddo
    endif

    do b = 1, Nblocks
      allocate(orion%block(b)%mesh(1:ndir,0:orion%block(b)%Ni,0:orion%block(b)%Nj,0:orion%block(b)%Nk))
      do d = 1, ndir
        do k = 0, orion%block(b)%Nk; do j = 0, orion%block(b)%Nj; do i = 0, orion%block(b)%Ni
              read(unit,*,iostat=err) orion%block(b)%mesh(d,i,j,k)
        enddo; enddo; enddo
      enddo
    enddo

    close(unit)

  endfunction p3d_read_multiblock


  end module Lib_PLOT3D