  program vtk_write_multiblock
  USE IR_Precision
  USE Lib_VTK, only: vtk_write_structured_multiblock
  use Lib_ORION_data
  USE, intrinsic:: ISO_FORTRAN_ENV, only: stdout=>OUTPUT_UNIT, stderr=>ERROR_UNIT
  !---------------------------------------------------------------------------------------------------------------------------------
  !> Procedure for testing multi-blocks writing.
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(orion_data)       :: orion
  integer(I4P), parameter:: nx1=0_I4P,nx2=90_I4P,ny1=0_I4P,ny2=50_I4P
  integer(I4P)           :: nz1 ,nz2
  integer(I4P)           :: i,j,k,b,E_IO
  character(len=10)      :: varnames
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! 3D test

  nz1 = 0_I4P; nz2 = 50_I4P
  allocate(orion%block(1:4))
  do b = 1, 4
    orion%block(b)%name = 'B3D'//trim(str(.true.,b))
    orion%block(b)%Ni = nx2
    orion%block(b)%Nj = ny2
    orion%block(b)%Nk = nz2
    allocate(orion%block(b)%mesh(1:3,nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(orion%block(b)%vars(1:2,nx1+1:nx2,ny1+1:ny2,nz1+1:nz2))
    do k=nz1,nz2
      do j=ny1,ny2
        do i=nx1,nx2
        orion%block(b)%mesh(1,i,j,k) = i*1._R8P
        orion%block(b)%mesh(2,i,j,k) = j*1._R8P
        orion%block(b)%mesh(3,i,j,k) = k*1._R8P
        enddo
      enddo
    enddo
    do k=nz1+1,nz2
      do j=ny1+1,ny2
        do i=nx1+1,nx2
        orion%block(b)%vars(1,i,j,k) = (i*j*k)*1._R8P
        orion%block(b)%vars(2,i,j,k) = -(i*j*k)*1._R8P
        enddo
      enddo
    enddo
    if (b>1) then
      orion%block(b)%mesh(1,:,:,:) = orion%block(b-1)%mesh(1,:,:,:) + nx2*1._R8P
    endif
  enddo

  varnames='var1 var2'
  orion%vtk%format='raw'
  E_IO = vtk_write_structured_multiblock(orion=orion,vtspath='',vtmpath='field3D',varnames=varnames)

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Pure 2D test
  
  deallocate(orion%block)
  allocate(orion%block(1:4))
  do b = 1, 4
    orion%block(b)%name = 'B2D'//trim(str(.true.,b))
    orion%block(b)%Ni = nx2
    orion%block(b)%Nj = ny2
    orion%block(b)%Nk = 1_I4P
    allocate(orion%block(b)%mesh(1:2,nx1:nx2,ny1:ny2,0:0))
    allocate(orion%block(b)%vars(1:2,nx1+1:nx2,ny1+1:ny2,1:1))
    do j=ny1,ny2
      do i=nx1,nx2
      orion%block(b)%mesh(1,i,j,0) = i*1._R8P
      orion%block(b)%mesh(2,i,j,0) = j*1._R8P
      enddo
    enddo
    do j=ny1+1,ny2
      do i=nx1+1,nx2
      orion%block(b)%vars(1,i,j,1) = (i*j*k)*1._R8P
      orion%block(b)%vars(2,i,j,1) = -(i*j*k)*1._R8P
      enddo
    enddo
    if (b>1) then
      orion%block(b)%mesh(1,:,:,:) = orion%block(b-1)%mesh(1,:,:,:) + nx2*1._R8P
    endif
  enddo

  varnames='var1 var2'
  orion%vtk%format='raw'
  E_IO = vtk_write_structured_multiblock(orion=orion,vtspath='',vtmpath='field2D',varnames=varnames)

  end program vtk_write_multiblock