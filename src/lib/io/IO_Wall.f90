!>@brief: .
module MOSE_IO_Wall
  use Lib_ORION_data
  use MOSE_Parameters_m
  use MOSE_Global_m

  implicit none
  logical, allocatable :: viscous_flag(:,:)
  logical, allocatable :: coupling_flag(:,:)
  integer              :: Onwall
  type(ORION_data)     :: IOwall
  character(len=llen)  :: Owallnames
  character(len=llen)  :: wall_variables
  logical              :: write_mechanical, write_thermal, write_mass, write_composition

contains


  subroutine Initialize_Wall_File ( domain, extension )
    use MOSE_Types_m
    use IR_Precision
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    character(len=*), intent(in)       :: extension
    ! Local
    integer :: b, f, bb, i, j, k

    IOwall%tec%extension = extension

    write_mechanical = .false.
    write_thermal = .false.
    write_mass = .false.
    write_composition = .false.

    Onwall = 1
    Owallnames = '"y+"'
    if (index(wall_variables,'mech')>0) then
      Onwall = Onwall + 4
      Owallnames = trim(Owallnames)//' "tauX" "tauY" "tauZ" "pw"'
      write_mechanical = .true.
    endif
    if (index(wall_variables,'thermal')>0) then
      Onwall = Onwall + 2
      Owallnames = trim(Owallnames)//' "Tw" "qw"'
      write_thermal = .true.
    endif
    if (index(wall_variables,'mass')>0) then
      Onwall = Onwall + 1
      Owallnames = trim(Owallnames)//' "mdot"'
      write_mass = .true.
    endif

    allocate(IOwall%block(1:count(viscous_flag)))

    bb = 0
    do b = 1, domain % nb; do f = 1, 6
      if (viscous_flag( b , f )) then
        bb = bb + 1
        IOwall%block(bb)%name = 'B'//trim(str(.true.,b))//'F'//trim(str(.true.,f))
        select case(f)
        case(1)
          allocate(IOwall%block(bb)%mesh(1:3,0:0,0:domain%blk(b)%dim(2),0:domain%blk(b)%dim(3)))
          allocate(IOwall%block(bb)%vars(1:Onwall,1,1:domain%blk(b)%dim(2),1:domain%blk(b)%dim(3)))
          do k = 0, domain%blk(b)%dim(3); do j = 0, domain%blk(b)%dim(2)
              IOwall % block(bb) % mesh(1,0,j,k) = domain % blk(b) % node(0,j,k) % c(1)
              IOwall % block(bb) % mesh(2,0,j,k) = domain % blk(b) % node(0,j,k) % c(2)
              IOwall % block(bb) % mesh(3,0,j,k) = domain % blk(b) % node(0,j,k) % c(3)
              ! Section to fool compiler...
              if (IOwall%block(bb)%mesh(1,0,j,k)==0d0 .and. &
                  IOwall%block(bb)%mesh(2,0,j,k)==0d0 .and. &
                  IOwall%block(bb)%mesh(3,0,j,k)==0d0) &
              print*, IOwall % block(bb) % mesh(:,0,j,k)
          enddo; enddo
        case(2)
          allocate(IOwall%block(bb)%mesh(1:3,0:0,0:domain%blk(b)%dim(2),0:domain%blk(b)%dim(3)))
          allocate(IOwall%block(bb)%vars(1:Onwall,1,1:domain%blk(b)%dim(2),1:domain%blk(b)%dim(3)))
          do k = 0, domain%blk(b)%dim(3); do j = 0, domain%blk(b)%dim(2)
              IOwall % block(bb) % mesh(:,0,j,k) = domain % blk(b) % node(domain%blk(b)%dim(1),j,k) % c
          enddo; enddo
        case(3)
          allocate(IOwall%block(bb)%mesh(1:3,0:domain%blk(b)%dim(1),0:0,0:domain%blk(b)%dim(3)))
          allocate(IOwall%block(bb)%vars(1:Onwall,1:domain%blk(b)%dim(1),1,1:domain%blk(b)%dim(3)))
          do k = 0, domain%blk(b)%dim(3); do i = 0, domain%blk(b)%dim(1)
              IOwall % block(bb) % mesh(:,i,0,k) = domain % blk(b) % node(i,0,k) % c
          enddo; enddo
        case(4)
          allocate(IOwall%block(bb)%mesh(1:3,0:domain%blk(b)%dim(1),0:0,0:domain%blk(b)%dim(3)))
          allocate(IOwall%block(bb)%vars(1:Onwall,1:domain%blk(b)%dim(1),1,1:domain%blk(b)%dim(3)))
          do k = 0, domain%blk(b)%dim(3); do i = 0, domain%blk(b)%dim(1)
              IOwall % block(bb) % mesh(:,i,0,k) = domain % blk(b) % node(i,domain%blk(b)%dim(2),k) % c
          enddo; enddo
        case(5)
          allocate(IOwall%block(bb)%mesh(1:3,0:domain%blk(b)%dim(1),0:domain%blk(b)%dim(2),0:0))
          allocate(IOwall%block(bb)%vars(1:Onwall,1:domain%blk(b)%dim(1),1:domain%blk(b)%dim(2),1))
          do j = 0, domain%blk(b)%dim(2); do i = 0, domain%blk(b)%dim(1)
              IOwall % block(bb) % mesh(:,i,j,0) = domain % blk(b) % node(i,j,0) % c
          enddo; enddo
        case(6)
          allocate(IOwall%block(bb)%mesh(1:3,0:domain%blk(b)%dim(1),0:domain%blk(b)%dim(2),0:0))
          allocate(IOwall%block(bb)%vars(1:Onwall,1:domain%blk(b)%dim(1),1:domain%blk(b)%dim(2),1))
          do j = 0, domain%blk(b)%dim(2); do i = 0, domain%blk(b)%dim(1)
              IOwall % block(bb) % mesh(:,i,j,0) = domain % blk(b) % node(i,j,domain%blk(b)%dim(3)) % c
          enddo; enddo
        end select
      endif
    enddo; enddo

  end subroutine Initialize_Wall_File


  !> Update the orion-field data and write them accordingly to the chosen format (vtk,tecplot)
  subroutine Write_Wall_Solution( domain, file, format )
    use Lib_VTK
    use Lib_Tecplot
    use MOSE_Types_m
    use MOSE_Parameters_m
    use MOSE_Mod_MPI, only: mpi_is_root
    use MOSE_Mod_GhostExchange, only: mpi_io_barrier

    implicit none
    type(MOSE_domain_type), intent(inout)  :: domain
    character(llen), intent(in)            :: file
    character(clen), intent(in)            :: format(2)
    ! Local
    character(len=llen)           :: path
    character(len=llen)           :: localpath_vtk
    integer                       :: E_IO

    ! Each rank computes wall properties for its local blocks and reduces to root.
    ! (collective — all ranks must call)
    call Extract_Properties_From_BC( domain )

    if (mpi_is_root) then
      path = 'OUTPUT/'

      ! Write the IOwall accordingly to the solution format
      select case(format(1))
      case('vtk')
        IOwall%vtk%format = format(2)
        localpath_vtk = trim(path)//'vtk/'
        call execute_command_line('mkdir -p '//trim(localpath_vtk))
        E_IO = vtk_write_structured_multiblock(orion=IOwall,vtspath=trim(localpath_vtk)//trim(file), &
                                                             vtmpath=trim(path)//trim(file),varnames=Owallnames,time=domain%time)
      case('tecplot')
        IOwall%tec%format = format(2)
        E_IO = tec_write_structured_multiblock(orion=IOwall,varnames=Owallnames,filename=trim(path)//trim(file)//trim(IOwall%tec%extension))
      end select
    end if
    ! Synchronize all ranks after I/O
    call mpi_io_barrier()

  end subroutine Write_Wall_Solution


  subroutine Extract_Properties_From_BC ( domain )
    use MOSE_Types_m
    use MOSE_Global_m
    use MOSE_Lib_BC
    use MOSE_Mod_MPI, only: mpi_is_root, is_local_block, mpi_reduce_sum_r8_array

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    real(kind=8), dimension(8) :: Ovar
    integer :: f, lower, upper, i, b, bb, pv, Bm, Im, Jm, Km, Fm, Is, Js, Ks
    integer :: n_wall_total, wall_i
    real(R8), allocatable :: wall_buf(:)

    ! Count total BC cells across all viscous faces
    n_wall_total = 0
    do b = 1, domain%nb
      do f = 1, 6
        if (viscous_flag(b,f)) n_wall_total = n_wall_total + domain%n_bf(b,f)
      end do
    end do

    allocate(wall_buf(8 * n_wall_total))
    wall_buf = 0d0

    ! Each rank computes Ovar only for its local blocks
    upper = 0; wall_i = 0
    do b = 1, domain%nb
      do f = 1, 6
        lower = upper + 1                ! Update lower bound
        upper = upper + domain%n_bf(b,f) ! Upper bound: add number of cells on face f of block b
 
        if (.not. viscous_flag(b,f)) cycle

        do i = lower, upper
          wall_i = wall_i + 1
          Bm = domain%bc(i)%b
          if (.not. is_local_block(Bm)) cycle

          Im = domain % bc(i) % i
          Jm = domain % bc(i) % j
          Km = domain % bc(i) % k
          Fm = domain % bc(i) % f
          Ovar = 0d0
          select case ( domain % bc(i) % type )

            case (5) ! wall q
              call BC_Wall_Heat ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % qw, Ovar )

            case (6) ! wall T
              call BC_Wall_Temperature ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Tw, Ovar )

            case (8) ! melting fuel - paraffin
              call BC_Wall_Melting ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Tw, domain % bc(i) % qrad, Ovar )

            case (9,10,12,13) ! pyrolyzing fuel - HDPE, CC, ecc
              call BC_Wall_Pyrolysis ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % qrad, domain%bc(i) % type, Ovar )

          end select
          wall_buf(8*(wall_i-1)+1 : 8*wall_i) = Ovar
        enddo

      enddo ! face
    enddo ! block

    ! Reduce wall_buf to root. Each cell is owned by exactly one rank, so SUM gives
    ! correct values. On root, wall_buf is updated in-place by mpi_reduce_sum_r8_array.
    call mpi_reduce_sum_r8_array(wall_buf, 8*n_wall_total)

    ! Root fills IOwall from the gathered buffer
    if (mpi_is_root) then
      upper = 0; wall_i = 0; bb = 0
      do b = 1, domain%nb
        do f = 1, 6
          lower = upper + 1
          upper = upper + domain%n_bf(b,f)
          if (.not. viscous_flag(b,f)) cycle
          bb = bb + 1

          do i = lower, upper
            wall_i = wall_i + 1
            Ovar = wall_buf(8*(wall_i-1)+1 : 8*wall_i)

            select case(f)
              case(1,2); Is = 1               ; Js = domain%bc(i)%j ; Ks = domain%bc(i)%k
              case(3,4); Is = domain%bc(i)% i  ; Js = 1              ; Ks = domain%bc(i)%k
              case(5,6); Is = domain%bc(i)% i  ; Js = domain%bc(i)%j ; Ks = 1
            end select
            IOwall%block(bb)%vars(1, Is, Js, Ks) = Ovar(1)
            pv = 1
            if (write_mechanical) then
              IOwall%block(bb)%vars(pv+1:pv+4, Is, Js, Ks) = Ovar(2:5)
              pv = pv + 4
            endif
            if (write_thermal) then
              IOwall%block(bb)%vars(pv+1:pv+2, Is, Js, Ks) = Ovar(6:7)
              pv = pv + 2
            endif
            if (write_mass) then
              IOwall%block(bb)%vars(pv+1, Is, Js, Ks) = Ovar(8)
              pv = pv + 2
            endif
          enddo

        enddo ! face
      enddo ! block
    end if

    deallocate(wall_buf)

  end subroutine Extract_Properties_From_BC
  
end module MOSE_IO_Wall