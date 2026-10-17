program tecplot_read_multiblock
  use Lib_Tecplot
  use Lib_ORION_data
  implicit none
  type(orion_data) :: data
  integer          :: error, i, nb

  data%tec%node = .false.
  data%tec%bc = .false.

  write(*,*) '----------------------------------'
  write(*,*) ' TEC ASCII FILES '
  write(*,*) '----------------------------------'

  data%tec%format = 'ascii'
  error = tec_read_structured_multiblock(orion=data,filename='tecfile.tec')
  nb = size(data%block)

  write(*,*) '3D domain'
  write(*,*) 'Blocks number = ', nb
  do i = 1, nb
    write(*,*) 'Block size    = ', data%block(i)%Ni, data%block(i)%Nj, data%block(i)%Nk
  enddo

  deallocate(data%block)

  error = tec_read_structured_multiblock(orion=data,filename='tec2D.tec')
  nb = size(data%block)

  write(*,*) '2D domain'
  write(*,*) 'Blocks number = ', nb
  do i = 1, nb
    write(*,*) 'Block size    = ', data%block(i)%Ni, data%block(i)%Nj, data%block(i)%Nk
  enddo

  deallocate(data%block)

  write(*,*) '----------------------------------'
  write(*,*) ' SZPLT BINARY FILES '
  write(*,*) '----------------------------------'

  data%tec%format = 'binary'
  error = tec_read_structured_multiblock(orion=data,filename='tecfile.szplt')
  nb = size(data%block)

  write(*,*) '3D domain'
  write(*,*) 'Blocks number = ', nb
  do i = 1, nb
    write(*,*) 'Block name    = ', trim(data%block(i)%name) 
    write(*,*) 'Block size    = ', data%block(i)%Ni, data%block(i)%Nj, data%block(i)%Nk
  enddo

  deallocate(data%block)

  error = tec_read_structured_multiblock(orion=data,filename='tec2D.szplt')
  nb = size(data%block)

  write(*,*) '2D domain'
  write(*,*) 'Blocks number = ', nb
  do i = 1, nb
    write(*,*) 'Block size    = ', data%block(i)%Ni, data%block(i)%Nj, data%block(i)%Nk
  enddo

end program tecplot_read_multiblock