!> @ingroup Library
!> @{
!> @defgroup Lib_ORION_dataLibrary Lib_ORION_data
!> @}

!> Module to store the ORION main data structure
!> @ingroup Lib_ORION_dataLibrary
module Lib_ORION_data
  use IR_Precision
  implicit none
  save
  private

  public :: copyORION

  type :: Type_tec_Format
    character(6):: extension = '.tec'    !> File extension.
    character(6):: format    = 'binary'  !> Binary or ascii file.
    logical     :: node      = .false.   !> Node or cell data location.
    logical     :: bc        = .false.   !> Saving or not boundary conditions cells.
  endtype Type_tec_Format

  type :: Type_vtk_Format
    character(6):: format    = 'binary'  !> Binary or ascii file.
    logical     :: node      = .false.   !> Node or cell data location.
  endtype Type_vtk_Format

  type :: Type_p3d_Format
    character(6):: format    = 'binary'  !> Binary or ascii file.
  endtype Type_p3d_Format

  type :: obj_block
    character(len=128) :: name
    integer :: Ni, Nj, Nk
    real(R8P), dimension(:,:,:,:), allocatable :: mesh
    real(R8P), dimension(:,:,:,:), allocatable :: vars
  end type obj_block

  type, public :: orion_data
    character(len=32), allocatable :: varnames(:)
    real(R8P) :: solutiontime
    type(obj_block), allocatable :: block(:)
    type(Type_tec_Format) :: tec
    type(Type_vtk_Format) :: vtk
    type(Type_p3d_Format) :: p3d
  endtype orion_data

contains

  subroutine copyORION(a,b)
    type(orion_data), intent(in)  :: a
    type(orion_data), intent(out) :: b
    integer :: i

    ! Copy non-allocatable members
    b%tec%format = a%tec%format
    b%vtk%format = a%vtk%format
    b%p3d%format = a%p3d%format

    ! Handle allocation of block array
    if (allocated(b%block)) then
      if (size(b%block) /= size(a%block)) then
        deallocate(b%block)
        allocate(b%block(size(a%block)))
      end if
    else
      allocate(b%block(size(a%block)))
    end if
  
    if (allocated(a%varnames)) then
      ! Handle allocation of varnames
      if (allocated(b%varnames)) then
        if (size(b%varnames) /= size(a%varnames)) then
          deallocate(b%varnames)
          allocate(b%varnames(size(a%varnames)))
        end if
      else
        allocate(b%varnames(size(a%varnames)))
      end if
      b%varnames = a%varnames
    endif
 

    ! Copy each block
    do i = 1, size(a%block)
      b%block(i)%name = a%block(i)%name
      b%block(i)%Ni   = a%block(i)%Ni
      b%block(i)%Nj   = a%block(i)%Nj
      b%block(i)%Nk   = a%block(i)%Nk

      ! Handle mesh allocation
      if (allocated(b%block(i)%mesh)) then
        if (any(shape(b%block(i)%mesh) /= shape(a%block(i)%mesh))) then
          deallocate(b%block(i)%mesh)
          allocate(b%block(i)%mesh, mold=a%block(i)%mesh)
        end if
      else
        allocate(b%block(i)%mesh, mold=a%block(i)%mesh)
      end if
      b%block(i)%mesh = a%block(i)%mesh

      if (allocated(a%block(i)%vars)) then
        ! Handle vars allocation
        if (allocated(b%block(i)%vars)) then
          if (any(shape(b%block(i)%vars) /= shape(a%block(i)%vars))) then
            deallocate(b%block(i)%vars)
            allocate(b%block(i)%vars, mold=a%block(i)%vars)
          end if
        else
          allocate(b%block(i)%vars, mold=a%block(i)%vars)
        end if
        b%block(i)%vars = a%block(i)%vars
      endif

    end do

  end subroutine copyORION

end module Lib_ORION_data