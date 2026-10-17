program p3d_write
  use Lib_PLOT3D
  use Lib_ORION_data
  implicit none
  integer, parameter :: nx1=0                              !< X lower bound extent.
  integer, parameter :: nx2=10                             !< X upper bound extent.
  integer, parameter :: ny1=0                              !< Y lower bound extent.
  integer, parameter :: ny2=10                             !< Y upper bound extent.
  integer, parameter :: nz1=0                              !< Z lower bound extent.
  integer, parameter :: nz2=10                             !< Z upper bound extent.
  type(orion_data)   :: data
  integer            :: error                              !< Status error.
  integer            :: i, j, k                            !< Counter.

  allocate(data%block(1:2))

  data%block(1)%Ni = nx2
  data%block(1)%Nj = ny2
  data%block(1)%Nk = nz2

  data%block(2)%Ni = nx2+10
  data%block(2)%Nj = ny2
  data%block(2)%Nk = nz2

  allocate(data%block(1)%mesh(1:3,nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(data%block(2)%mesh(1:3,nx1:nx2+10,ny1:ny2,nz1:nz2))

  do k=nz1, nz2
    do j=ny1, ny2
      do i=nx1, nx2
        data%block(1)%mesh(1,i,j,k) = dble(i)
        data%block(1)%mesh(2,i,j,k) = dble(j)
        data%block(1)%mesh(3,i,j,k) = dble(k)
      enddo
    enddo
  enddo

  do k=nz1, nz2
    do j=ny1, ny2
      do i=nx1, nx2+10
        data%block(2)%mesh(1,i,j,k) = data%block(1)%mesh(1,nx2,j,k)+dble(i)
        data%block(2)%mesh(2,i,j,k) = dble(j)
        data%block(2)%mesh(3,i,j,k) = dble(k)
      enddo
    enddo
  enddo

  data%p3d%format = 'ascii'
  error = p3d_write_multiblock(orion=data,filename='mesh')

  ! 2D
  deallocate(data%block)
  allocate(data%block(1:2))

  data%block(1)%Ni = nx2
  data%block(1)%Nj = ny2
  data%block(1)%Nk = 0

  data%block(2)%Ni = nx2+10
  data%block(2)%Nj = ny2
  data%block(2)%Nk = 0

  allocate(data%block(1)%mesh(1:2,nx1:nx2,ny1:ny2,0:0))
  allocate(data%block(2)%mesh(1:2,nx1:nx2+10,ny1:ny2,0:0))

  do j=ny1, ny2
    do i=nx1, nx2
      data%block(1)%mesh(1,i,j,0) = dble(i)
      data%block(1)%mesh(2,i,j,0) = dble(j)
    enddo
  enddo

  do j=ny1, ny2
    do i=nx1, nx2+10
      data%block(2)%mesh(1,i,j,0) = data%block(1)%mesh(1,nx2,j,0)+dble(i)
      data%block(2)%mesh(2,i,j,0) = dble(j)
    enddo
  enddo

  data%p3d%format = 'ascii'
  error = p3d_write_multiblock(orion=data,filename='mesh2D')

end program p3d_write