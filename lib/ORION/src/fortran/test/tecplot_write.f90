program tecplot_write_multiblock
  use Lib_Tecplot
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

  data%tec%node = .false.
  data%tec%bc = .false.

  allocate(data%block(1:2))

  data%block(1)%name = 'blocco-A'
  data%block(2)%name = 'blocco-B'

  allocate(data%block(1)%mesh(1:3,nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(data%block(1)%vars(1:2,nx1+1:nx2,ny1+1:ny2,nz1+1:nz2))
  allocate(data%block(2)%mesh(1:3,nx1:nx2+10,ny1:ny2,nz1:nz2))
  allocate(data%block(2)%vars(1:2,nx1+1:nx2+10,ny1+1:ny2,nz1+1:nz2))

  do k=nz1, nz2
    do j=ny1, ny2
      do i=nx1, nx2
        data%block(1)%mesh(1,i,j,k) = dble(i)
        data%block(1)%mesh(2,i,j,k) = dble(j)
        data%block(1)%mesh(3,i,j,k) = dble(k)
      enddo
    enddo
  enddo
  do k=nz1+1, nz2
    do j=ny1+1, ny2
      do i=nx1+1, nx2
        data%block(1)%vars(1,i,j,k) = dble(i*j*k)
        data%block(1)%vars(2,i,j,k) = -dble(i*j*k)
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
  do k=nz1+1, nz2
    do j=ny1+1, ny2
      do i=nx1+1, nx2+10
        data%block(2)%vars(1,i,j,k) = dble(i*j*k)
        data%block(2)%vars(2,i,j,k) = -dble(i*j*k)
      enddo
    enddo
  enddo

  data%tec%format = 'binary'
  error = tec_write_structured_multiblock(orion=data,varnames='variable1 variable2',filename='tecfile.szplt')

  data%tec%format = 'binary'
  error = tec_write_structured_multiblock(orion=data,varnames='variable1 variable2',filename='tecfile.plt')

  data%tec%format = 'ascii'
  error = tec_write_structured_multiblock(orion=data,varnames='"variable1"',filename='tecfile.tec',Nvars=1)

  deallocate(data%block)

  allocate(data%block(1:2))

  data%block(1)%name = 'blocco-A'
  data%block(2)%name = 'blocco-B'

  allocate(data%block(1)%mesh(1:2,nx1:nx2,ny1:ny2,0:0))
  allocate(data%block(1)%vars(1:2,nx1+1:nx2,ny1+1:ny2,1:1))
  allocate(data%block(2)%mesh(1:2,nx1:nx2+10,ny1:ny2,0:0))
  allocate(data%block(2)%vars(1:2,nx1+1:nx2+10,ny1+1:ny2,1:1))

  do j=ny1, ny2
    do i=nx1, nx2
      data%block(1)%mesh(1,i,j,0) = dble(i)
      data%block(1)%mesh(2,i,j,0) = dble(j)
    enddo
  enddo

  do j=ny1+1, ny2
    do i=nx1+1, nx2
      data%block(1)%vars(1,i,j,1) = dble(i*j)
      data%block(1)%vars(2,i,j,1) = -dble(i*j)
    enddo
  enddo

  do j=ny1, ny2
    do i=nx1, nx2+10
      data%block(2)%mesh(1,i,j,0) = data%block(1)%mesh(1,nx2,j,0)+dble(i)
      data%block(2)%mesh(2,i,j,0) = dble(j)
    enddo
  enddo
  do j=ny1+1, ny2
    do i=nx1+1, nx2+10
      data%block(2)%vars(1,i,j,1) = dble(i*j)
      data%block(2)%vars(2,i,j,1) = -dble(i*j)
    enddo
  enddo

  data%tec%format = 'binary'
  error = tec_write_structured_multiblock(orion=data,varnames='variable1 variable2',filename='tec2D.szplt')

  data%tec%format = 'ascii'
  error = tec_write_structured_multiblock(orion=data,varnames='"variable1"',filename='tec2D.tec',Nvars=1)

end program tecplot_write_multiblock