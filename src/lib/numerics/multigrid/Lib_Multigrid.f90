module MOSE_Lib_Multigrid
  use MOSE_IO_Ini_Open, only: MOSE_Load_Ini   ! one ini loader: MOSE_open dialect translated in memory

  implicit none

contains


  subroutine Setup_Multigrid ( simulation )
    use MOSE_Types_m
    use MOSE_Parameters_m
    use Finer, only: file_ini
    use MOSE_Mod_Allocate_Data
  
    implicit none
    type(MOSE_simulation_type), intent(inout) :: simulation
    ! Local
    type(file_ini)   :: fini
    integer          :: b, m, i, j, k, nb, rap, error
    character(clen)  :: option


    call Check_Multigrid ( simulation%domain(1), simulation%config%MGL )

    nb = simulation%domain(1) % nb

    do m = 2, simulation%config%MGL
      
      allocate ( simulation%domain(m) % Blk(nb) )
      simulation%domain(m) % nb = nb
      simulation%domain(m) % time = simulation%domain(1) % time

      rap = 2
      do b = 1, nb
        simulation%domain(m) % Blk(b) % Dim = simulation%domain(m-1) % Blk(b) % Dim / rap
        simulation%domain(m) % Blk(b) % Dim(3) = Max ( 1, simulation%domain(m) % Blk(b) % Dim(3) ) ! 2D case
        call Allocate_Block ( simulation%domain(m) % Blk(b), simulation%domain(m) % Blk(b) % Dim, simulation%config%MGL )
      enddo

      call Coarse_Grid ( simulation%domain(m-1), simulation%domain(m) )

      call Coarse_IOfield ( simulation%IOfield(m-1), simulation%IOfield(m) )

      do b = 1, nb
        do k = 0, simulation%domain(m)%Blk(b)%Dim(3)
        do j = 0, simulation%domain(m)%Blk(b)%Dim(2)
        do i = 0, simulation%domain(m)%Blk(b)%Dim(1)
          simulation%IOfield(m)%block(b)%mesh(:,i,j,k) = simulation%domain(m) % Blk(b) % Node(i,j,k) % c
        enddo; enddo; enddo
      enddo

    enddo

    ! Read number of iterations for each domain
    if (.not. simulation%config%HYDRA_MG) then
      call MOSE_Load_Ini(fini)
      do m = 1, simulation%config%MGL
        write (option,'(A5,I0,A5)') 'level', m, '-iter'
        call fini%get(section_name=trim(codename)//'-multigrid', option_name=option, &
                      val=simulation%domain(m)%itermax, error=error)
        if ( error/=0 .and. m>1 ) simulation%domain(m)%itermax = 0 ! for grid1 itermax from -parameters section
      enddo
    endif

  end subroutine Setup_Multigrid


  subroutine Check_Multigrid ( domain, MGL )
    use MOSE_Types_m

    implicit none
    integer, intent(in)                   :: MGL
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, d, rap, check


    rap = 2**(MGL-1)

    do b = 1, domain % nb

      check = 0

      do d = 1, 3

        if ( Mod ( domain % Blk(b) % Dim(d), rap ) == 0 ) check = check + 1
        
        if (d == 3) then
          if ( domain % Blk(b) % Dim(d) == 1 ) check = check + 1
        endif

      enddo

      if ( check < 3 ) then
        write(*,*) ' Error in Check_Multigrid, block: ', b
        stop
      endif

    enddo

  end subroutine Check_Multigrid


  subroutine Coarse_Grid ( Fine, Coarse )
    use MOSE_Types_m

    implicit none
    type(MOSE_domain_type), intent(inout) :: Fine, Coarse
    ! Local
    integer :: b


    do b = 1, Coarse % nb ! Loop over blocks

      call Coarse_Grid_Blk ( Fine % Blk(b) % node, Coarse % Blk(b) % node, &
                             Fine % Blk(b) % dim,  Coarse % Blk(b) % dim )

    enddo

  end subroutine Coarse_Grid


  subroutine Coarse_Grid_Blk ( fNode, cNode, fDim, cDim )
    use MOSE_General_Types_m
    use MOSE_Global_m

    implicit none
    integer, intent(in), dimension(3) :: fDim, cDim
    type(MOSE_vector_3D_type), intent(in),  dimension(0:fDim(1),0:fDim(2),0:fDim(3)) :: fNode
    type(MOSE_vector_3D_type), intent(out), dimension(0:cDim(1),0:cDim(2),0:cDim(3)) :: cNode
    ! Local
    integer :: i, j, k, i2, j2, k2


    !$omp do collapse (2)
    do k = 0, fDim(3), 2-Mod(fDim(3),2)
    do j = 0, fDim(2), 2
    do i = 0, fDim(1), 2
      
      i2 = i / 2
      j2 = j / 2
      k2 = k / 2

      if ( fDim(3) == 1 ) then ! 2D
        k2 = k
      end if

      cNode(i2,j2,k2) % c = fNode(i,j,k) % c
    
    enddo; enddo; enddo
  
  end subroutine Coarse_Grid_Blk


  subroutine Coarse_IOfield ( fIOfield, cIOfield )
    use Lib_ORION_data
    use MOSE_Global_m

    implicit none
    type(ORION_data), intent(in)    :: fIOfield
    type(ORION_data), intent(inout) :: cIOfield
    ! Local
    integer :: b, nb, Ni, Nj, Nk, Onvar

    Onvar = nprim + 3  ! T, gamma, R

    nb = size( fIOfield % block )
    allocate ( cIOfield % block( nb ) )
    
    do b = 1, nb
      Ni = fIOfield % block(b)%Ni /2
      Nj = fIOfield % block(b)%Nj /2
      Nk = fIOfield % block(b)%Nk /2
      Nk = Max ( 1, Nk ) ! 2D case
      cIOfield % block(b) % Ni = Ni
      cIOfield % block(b) % Nj = Nj
      cIOfield % block(b) % Nk = Nk
      allocate( cIOfield%block(b)%mesh(1:3,0:Ni,0:Nj,0:Nk) )
      allocate( cIOfield%block(b)%vars(1:Onvar,Ni,Nj,Nk) )
    enddo

  end subroutine Coarse_IOfield

end module MOSE_Lib_Multigrid