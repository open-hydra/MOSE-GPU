module MOSE_Lib_Interpolation

  implicit none

contains


  subroutine Restriction ( Fine, Coarse )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block

    implicit none
    type(MOSE_domain_type), intent(inout) :: Fine, Coarse
    ! Local
    integer :: b


    do b = 1, Coarse % nb ! Loop over blocks
      if (.not. is_local_block(b)) cycle

      call fine2coarse_prim ( Fine % Blk(b) % P,   Coarse % Blk(b) % P  , &
                              Fine % Blk(b) % vol, Coarse % Blk(b) % vol, &
                              Fine % Blk(b) % dim, Coarse % Blk(b) % dim )

    enddo

  end subroutine Restriction


  subroutine Prolongation ( Fine, Coarse )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block

    implicit none
    type(MOSE_domain_type), intent(inout) :: Fine, Coarse
    ! Local
    integer :: b


    do b = 1, Coarse % nb ! Loop over blocks
      if (.not. is_local_block(b)) cycle

      call coarse2fine_prim ( Fine % Blk(b) % P, Coarse % Blk(b) % P, &
                              Fine % Blk(b) % dim, Coarse % Blk(b) % dim )

    enddo

  end subroutine Prolongation


  subroutine fine2coarse_prim ( fPrim, cPrim, fVol, cVol, fDim, cDim )
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic

    implicit none
    integer, intent(in), dimension(3) :: fDim, cDim
    real(kind=8), intent(in), dimension(nprim,1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fPrim
    real(kind=8), intent(out),dimension(nprim,1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cPrim
    real(kind=8), intent(in), dimension(1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fVol
    real(kind=8), intent(in), dimension(1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cVol
    ! Local
    integer :: i, j, k, i2, j2, k2, i2d, j2d, k2d
    real(kind=8) :: guess(nprim), rho, Rgas, T0
    real(kind=8), dimension(nprim,fDim(1),fDim(2),fDim(3)) :: fCons
    real(kind=8), dimension(nprim,cDim(1),cDim(2),cDim(3)) :: cCons


    !$omp do collapse (2)
    do k = 1, fDim(3)
    do j = 1, fDim(2)
    do i = 1, fDim(1)
      fCons(np+1:nprim,i,j,k) = fPrim(np+1:nprim,i,j,k)
      fCons(1:np,i,j,k) = prim2cons ( fPrim(1:np,i,j,k) )
    end do ; end do ; end do

    !$omp do collapse (2)
    do k = 1, cDim(3)
    do j = 1, cDim(2)
    do i = 1, cDim(1)

      i2 = 2*i
      j2 = 2*j
      k2 = 2*k
      i2d = i2-1
      j2d = j2-1
      k2d = k2-1

      if ( fDim(3) == 1 ) then ! 2D
        cCons(:,i,j,k) = fCons(:,i2d,j2d,1) * fVol(i2d,j2d,1)  &
                       + fCons(:,i2, j2d,1) * fVol(i2, j2d,1)  &
                       + fCons(:,i2d,j2, 1) * fVol(i2d,j2, 1)  &
                       + fCons(:,i2, j2, 1) * fVol(i2, j2, 1)
      else
        cCons(:,i,j,k) = fCons(:,i2d,j2d,k2d) * fVol(i2d,j2d,k2d)  &
                       + fCons(:,i2, j2d,k2d) * fVol(i2, j2d,k2d)  &
                       + fCons(:,i2d,j2, k2d) * fVol(i2d,j2, k2d)  &
                       + fCons(:,i2, j2, k2d) * fVol(i2, j2, k2d)  &
                       + fCons(:,i2d,j2d,k2 ) * fVol(i2d,j2d,k2 )  &
                       + fCons(:,i2, j2d,k2 ) * fVol(i2, j2d,k2 )  &
                       + fCons(:,i2d,j2, k2 ) * fVol(i2d,j2, k2 )  &
                       + fCons(:,i2, j2, k2 ) * fVol(i2, j2, k2 )
      end if ! 2D

    enddo; enddo; enddo

    !$omp do collapse (2)
    do k = 1, cDim(3)
    do j = 1, cDim(2)
    do i = 1, cDim(1)
      cCons(:,i,j,k) = cCons(:,i,j,k) / cVol(i,j,k)
      ! Guess for Newton-Raphson with one fine cell
      i2 = 2*i
      j2 = 2*j
      k2 = 2*k
      if ( fDim(3) == 1 ) k2 = 1
      guess = fPrim(:, i2, j2, k2)
      call co_rotot_Rtot( guess(:), rho, Rgas )
      T0 = guess(np) / ( rho * Rgas )
      cPrim(np+1:nprim,i,j,k) = cCons(np+1:nprim,i,j,k) ! turbulence variables
      cPrim(1:np,i,j,k) = cons2prim ( cCons(1:np,i,j,k), T0 ) ! ( rho U p )
    enddo; enddo; enddo

  end subroutine fine2coarse_prim


  subroutine coarse2fine_prim ( fPrim, cPrim, fDim, cDim )
    use MOSE_Global_m

    implicit none
    integer, intent(in), dimension(3) :: fDim, cDim
    real(kind=8), intent(out), dimension(nprim,1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fPrim
    real(kind=8), intent(in),  dimension(nprim,1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cPrim
    ! Local
    integer :: i, j, k, rap, i2, j2, k2, i2d, j2d, k2d, im, jm, km, ip, jp, kp
    integer :: ii, jj, kk, counter, mask(3), id(6)
    real(kind=8) :: a1, a2, a3, a4, coeffs(8), interp(nprim)


    a1 = 27d0/64d0
    a2 = 9d0/64d0
    a3 = 3d0/64d0
    a4 = 1d0/64d0
    coeffs(1:8) = [ a1, a2, a2, a2, a3, a3, a3, a4 ]
    rap = 2
    
    if ( fDim(3) == 1 ) then ! 2D
      a1 = 9.d0/16.d0
      a2 = 3.d0/16.d0
      a3 = 1.d0/16.d0
      a4 = 0d0
      coeffs(1:8) = [ a1, a2, a2, a4, a3, a4, a4, a4 ]
    end if

    !$omp do collapse (2)
    do k = 1, cDim(3)
    do j = 1, cDim(2)            
    do i = 1, cDim(1)
            
      i2 = rap*i
      j2 = rap*j
      k2 = rap*k
      i2d = i2-(rap-1)
      j2d = j2-(rap-1)
      k2d = k2-(rap-1)
          
      im = Max (1,i-1)
      jm = Max (1,j-1)
      km = Max (1,k-1)

      ip = Min (cDim(1), i+1)
      jp = Min (cDim(2), j+1)
      kp = Min (cDim(3), k+1)

      if ( fDim(3) == 1 ) then
        k2 = 1
        k2d = 1 
      end if ! 2D

      id(1:6) = [ im, jm, km, ip, jp, kp ]
          
      counter = 1
          
      do kk = k2d, k2          
      do jj = j2d, j2            
      do ii = i2d, i2

        if ( fDim(3) == 1 ) then ! 2D
          if (counter==1) mask(1:3) = [1,2,3]
          if (counter==2) mask(1:3) = [4,2,3]
          if (counter==3) mask(1:3) = [1,5,3]
          if (counter==4) mask(1:3) = [4,5,3]
        else
          if (counter==1) mask(1:3) = [1,2,3]
          if (counter==2) mask(1:3) = [4,2,3]
          if (counter==3) mask(1:3) = [1,5,3]
          if (counter==4) mask(1:3) = [4,5,3]
          if (counter==5) mask(1:3) = [1,2,6]
          if (counter==6) mask(1:3) = [4,2,6]
          if (counter==7) mask(1:3) = [1,5,6]
          if (counter==8) mask(1:3) = [4,5,6]
        endif  

        interp = coeffs(1) * cPrim(:,i,j,k) &
               + coeffs(2) * cPrim(:,id(mask(1)),j,k) &
               + coeffs(3) * cPrim(:,i,id(mask(2)),k)    &
               + coeffs(4) * cPrim(:,i,j,id(mask(3)))    &
               + coeffs(5) * cPrim(:,id(mask(1)),id(mask(2)),k)   &
               + coeffs(6) * cPrim(:,id(mask(1)),j,id(mask(3)))   &
               + coeffs(7) * cPrim(:,i,id(mask(2)),id(mask(3)))   &
               + coeffs(8) * cPrim(:,id(mask(1)),id(mask(2)),id(mask(3)))

        fPrim(:,ii,jj,kk) = interp

        counter = counter + 1          
        
      enddo; enddo; enddo
    
    enddo; enddo; enddo
      
  end subroutine coarse2fine_prim

  
end module MOSE_Lib_Interpolation