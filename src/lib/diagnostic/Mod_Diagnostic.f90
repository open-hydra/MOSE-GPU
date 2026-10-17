module MOSE_Mod_Diagnostic
  use MOSE_Parameters_m

  implicit none
  integer         :: unitRES             ! Residuals file unit
  character(llen) :: unitRES_format      ! Format for residual_history_file
  character(llen) :: Dvarnames = '"rho" "rhou" "rhov" "rhow" "rhoe" "dt" "beta"'

contains


  subroutine Compute_Residual ( new, old, dt, n, average, total )
    use MOSE_Global_m

    implicit none
    integer, intent(in)         :: n(3)
    real(kind=8), intent(in)    :: new ( nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc )
    real(kind=8), intent(in)    :: old ( nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc )
    real(kind=8), intent(in)    :: dt ( n(1), n(2), n(3) )
    real(kind=8), intent(out)   :: average(nres)
    real(kind=8), intent(inout) :: total(nres)
    ! Local
    integer :: i, j, k
    real(kind=8) :: resn(nres), residuo(nres)


    residuo = 0d0
    
    !$omp parallel
    !$omp do private ( resn ) reduction ( + : residuo ) collapse (3)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)
      
      resn(1) = abs ( sum(new(1:nsc,i,j,k) - old(1:nsc,i,j,k)) )  ! delta(rho)
      resn(2:nres) = abs ( new(nu:nprim,i,j,k) - old(nu:nprim,i,j,k) ) ! delta(...)
      residuo = residuo + resn*resn  ! sum ( r(i)**2 )

    enddo; enddo; enddo
    !$omp end parallel

    ! Average residual of this block
    average = sqrt ( residuo / float ( n(1)*n(2)*n(3) ) )
    ! Overall residual
    total = total + residuo

  end subroutine Compute_Residual


  !> Device residual (GPU, unconditional): device twin of Compute_Residual on the RESIDENT var-last P_flat/PO_flat
  !> (new state / old state of the last Strang half -- exactly what unpack_flat_to_P /
  !> unpack_flat_to_PO would have produced on the host). Same interior cells, same algebra;
  !> only the ACCESS (P/PO(:,i,j,k) -> P_flat/PO_flat(c0,:)) and the accumulation order (device
  !> reduction tree vs host loop) change -> residual-history.dat may shift in trailing
  !> digits; the field itself is untouched (the residual is diagnostic + restop only).
  !> Removes the residual-only host round-trip (P_flat + PO_flat D2H ~3.2 GB/iter at THOR scale
  !> + two host scatters) that the DEFAULT diter-res=1 forces every iteration.
  !>
  !> Hard-wired to the nres=7 layout (1 rho-sum + nu..nprim = 6 vars, i.e. RANS runs;
  !> nres = nprim-nsc+1): scalar reductions r1..r7 (OpenACC ARRAY reductions are not
  !> portable under nvfortran). The caller guards on nres==7 and keeps the host path
  !> otherwise. No dt argument: Compute_Residual never reads its dt.
  subroutine Compute_Residual_blk ( b, n, total )
    use MOSE_Global_m, only: gc, nsc, nu, nv, nw, np, nt, nprim, nres
    use MOSE_Field_m,  only: P_flat, PO_flat, field_block_index_params

    implicit none
    integer, intent(in)         :: b, n(3)
    real(kind=8), intent(inout) :: total(nres)
    ! Local
    integer      :: i, j, k, s, c0, boff, pi, pj
    real(kind=8) :: d, r1, r2, r3, r4, r5, r6, r7
    integer :: nb1_, nb2_, nb3_   ! grain-fix
    nb1_ = n(1); nb2_ = n(2); nb3_ = n(3)

    call field_block_index_params(b, boff, pi, pj)
    r1=0d0; r2=0d0; r3=0d0; r4=0d0; r5=0d0; r6=0d0; r7=0d0

    !$acc parallel loop collapse(3) gang vector present(P_flat, PO_flat) private(c0, d, s) &
    !$acc&  reduction(+:r1,r2,r3,r4,r5,r6,r7)
    do k = 1, nb3_
    do j = 1, nb2_
    do i = 1, nb1_
      c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      d = 0d0
      do s = 1, nsc
        d = d + ( P_flat(c0,s) - PO_flat(c0,s) )          ! delta(rho) = sum of species deltas
      end do
      r1 = r1 + d*d
      d = P_flat(c0,nu)    - PO_flat(c0,nu);    r2 = r2 + d*d
      d = P_flat(c0,nv)    - PO_flat(c0,nv);    r3 = r3 + d*d
      d = P_flat(c0,nw)    - PO_flat(c0,nw);    r4 = r4 + d*d
      d = P_flat(c0,np)    - PO_flat(c0,np);    r5 = r5 + d*d
      d = P_flat(c0,nt)    - PO_flat(c0,nt);    r6 = r6 + d*d
      d = P_flat(c0,nprim) - PO_flat(c0,nprim); r7 = r7 + d*d
    enddo; enddo; enddo

    total(1) = total(1) + r1
    total(2) = total(2) + r2
    total(3) = total(3) + r3
    total(4) = total(4) + r4
    total(5) = total(5) + r5
    total(6) = total(6) + r6
    total(7) = total(7) + r7

  end subroutine Compute_Residual_blk


  !> Update the orion-field data and write them accordingly to the chosen format (vtk,tecplot)
  subroutine Write_Diagnostic( domain, IOfield, file, format )
    use Lib_VTK
    use Lib_Tecplot
    use MOSE_Types_m
    use Lib_ORION_data
    use MOSE_Global_m
    use MOSE_Mod_MPI, only: mpi_is_root, is_local_block
    use MOSE_Mod_GhostExchange, only: mpi_io_barrier, gather_diagnostic_to_root

    implicit none
    type(MOSE_domain_type), intent(inout)  :: domain
    type(ORION_data), intent(inout)        :: IOfield
    character(llen), intent(in)            :: file
    character(clen), intent(in), optional  :: format(2)
    ! Local
    character(len=llen) :: path
    character(len=llen) :: localpath_vtk
    integer             :: E_IO, b, i, j, k

    ! r/beta/dtlocal are device-resident across the solve. Refresh the
    ! host copies before they are read here (and by the collective gather, which early-
    ! returns for a single rank so the read below is the actual host consumer).
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$acc update host(domain%blk(b)%r, domain%blk(b)%beta, domain%blk(b)%dtlocal)
    end do
    block
      use MOSE_Field_m, only: field_update_R_flat_host, unpack_flat_to_R_host
      ! The COMPLETE last-substage residual lives in the resident
      ! var-LAST R_flat (device writers accumulate into it; host bridge sources were folded
      ! in by field_add_R_to_R_flat). Rebuild the var-first r for this diagnostic on the
      ! host: D2H R_flat once, assign into blk%r interior. Also fixes the Stage-2/3 window
      ! where this output was missing the conv/diff/BC/SST content.
      call field_update_R_flat_host ()
      call unpack_flat_to_R_host ( domain )
    end block

    ! Gather R, dtlocal, beta from all ranks to root (collective)
    call gather_diagnostic_to_root(domain)

    if (mpi_is_root) then

      path = 'OUTPUT/'

      ! Update IOfield variables with domain residuals
      do b = 1, size(IOfield%block)
        IOfield%block(b)%vars(1,:,:,:) = sum(domain%blk(b)%r(1:nsc,1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk), dim=1)
        IOfield%block(b)%vars(2:5,:,:,:) = domain%blk(b)%r(nu:np,1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
        if ( nrans > 0 ) then ! turbulence variables
          IOfield%block(b)%vars(6:6+nrans-1,:,:,:) = domain%blk(b)%r(nt:nprim,1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
        end if
        ! Auxiliary variables: dt
        do k = 1, IOfield%block(b)%Nk ; do j = 1, IOfield%block(b)%Nj ; do i = 1, IOfield%block(b)%Ni
              IOfield%block(b)%vars(nres+1,i,j,k) = domain%blk(b)%dtlocal(i,j,k)
        enddo; enddo; enddo
        ! Auxiliary variables: beta
        do k = 1, IOfield%block(b)%Nk ; do j = 1, IOfield%block(b)%Nj ; do i = 1, IOfield%block(b)%Ni
              IOfield%block(b)%vars(nres+2,i,j,k) = domain%blk(b)%beta(i,j,k)
        enddo; enddo; enddo
      enddo

      ! Write the IOfield accordingly to the solution format
      select case(format(1))
      case('vtk')
        IOfield%vtk%format = format(2)
        localpath_vtk = trim(path)//'/vtk'
        call execute_command_line('mkdir -p '//trim(localpath_vtk))
        E_IO = vtk_write_structured_multiblock(orion=IOfield,vtspath=trim(localpath_vtk)//trim(file), &
                                                             vtmpath=trim(path)//trim(file),varnames=Dvarnames,time=domain%time)
      case('tecplot')
        IOfield%tec%format = format(2)
        E_IO = tec_write_structured_multiblock(Nvars=nres+2,orion=IOfield,varnames=Dvarnames,filename=trim(path)//trim(file)//'.tec')
      end select
    end if

    ! Synchronize all ranks after I/O
    call mpi_io_barrier()

  end subroutine Write_Diagnostic


end module MOSE_Mod_Diagnostic
