program p3d_read
  use Lib_PLOT3D
  use Lib_ORION_data
  implicit none
  type(orion_data) :: data
  integer          :: error, i, nb

  data%p3d%format = 'ascii'
  error = p3d_read_multiblock(orion=data,filename='mesh.p3d')

  nb = size(data%block)
  write(*,*) '3D domain'
  write(*,*) 'Blocks number = ', nb
  do i = 1, nb
    write(*,*) 'Block size    = ', data%block(i)%Ni, data%block(i)%Nj, data%block(i)%Nk
  enddo

  deallocate(data%block)

  write(*,*)
  data%p3d%format = 'ascii'
  error = p3d_read_multiblock(orion=data,filename='mesh2D.p3d')

  nb = size(data%block)
  write(*,*) '2D domain'
  write(*,*) 'Blocks number = ', nb
  do i = 1, nb
    write(*,*) 'Block size    = ', data%block(i)%Ni, data%block(i)%Nj, data%block(i)%Nk
  enddo

end program p3d_read