#ifdef FULL_ALGEBRA
module ls_solver_sdirk
      implicit none
      save
      integer :: nvar
      type LSdata
        double precision, allocatable :: fjac(:,:),e(:,:)
        integer, allocatable :: ip(:)
      end type LSdata
contains
      subroutine lss_jac(t,y,jac,data)
      implicit none
      type(LSdata) :: data
      double precision :: t, y(nvar)
      external :: jac
      call jac(nvar,t,y,data%fjac)
      end subroutine lss_jac

      subroutine lss_decomp(data,hgamma,ising)
      implicit none
      type(LSdata) :: data
      integer :: ising, i, j
      double precision :: hgamma
      ! prepare left hand side
      do j=1,nvar
        do i=1,nvar
           data%e(i,j) = -data%fjac(i,j)
        end do
        data%e(j,j) = data%e(j,j) + hgamma
      end do

      call dgetrf(nvar,nvar,data%e,nvar,data%ip,ising)
      end subroutine lss_decomp

      subroutine lss_solve(data,rhs)
      implicit none
      type(LSdata) :: data
      double precision ::rhs(nvar)
      integer :: ising
      call dgetrs( 'n', nvar, 1, data%e, nvar, data%ip, rhs, nvar, ising)
      end subroutine lss_solve
 
      subroutine lss_init(data,n,nn)
      implicit none
      type(LSdata) :: data
      integer :: n,state
      integer :: nn ! necessary only for standardize the subroutine
      nvar = n
      allocate(data%ip(nvar),data%fjac(nvar,nvar),data%e(nvar,nvar),STAT = state)
      if(state .ne. 0)stop 'Allocation error for ip, fjac or e in lapack_init'
      end subroutine lss_init
   
      subroutine lss_free(data)
      implicit none
      type(LSdata) :: data
      integer :: state
      deallocate(data%ip,data%fjac,data%e,STAT=state)
      if(state .ne.0) stop 'Deallocation error for fjac and e in lapack_free'
      end subroutine lss_free
end module ls_solver_sdirk
#endif

! NOT WORKING WITH SPARSE_UMF OR SPARSE_LU

#ifdef SPARSE_UMF
module umf_sdirk
    implicit none 
    save 
    integer :: nvar, nnz
!   ax(nnz), ai(nnz), ap(nvar+1)
    double precision ,allocatable:: fjac(:,:), ax(:)
    integer, allocatable :: ai(:), ap(:)
    double precision ::  control(20), info(90)
    integer*8 :: symbolic, numeric
contains
    subroutine umf_decomp(ising,hgamma)
    !convert the matrix e from coordinate format to column compressed format
    integer :: i,j,idx,found,ising
    double precision :: hgamma
!   free previously used memory
    if(symbolic .ne. 0) then
      call umf4fsym(symbolic)
      call umf4fnum(numeric)
    end if

    idx = 1
    ! c index starts from 0, while fortran one starts from 1
    do j = 1,nvar
      found = 0
      do i = 1, nvar
        if ( fjac(i,j) .ne. 0d0) then
!          ax(idx) = e(i,j)
          if( i .eq. j ) then
            ax(idx) = -fjac(i,j) + hgamma
          else 
            ax(idx) = -fjac(i,j)
          end if      
          ! write the row in the row indices vector ai
          ai(idx) = i - 1
          if(found .eq. 0) then
            found = idx
          end if
          idx = idx + 1
        end if
      end do
      !must also update the column pointers vector ap
      ap(j) = found - 1
    end do
    !last element in ap must be nnza+1
    ap(nvar+1) = nnz
   
    call umf4def(control)
    control(1) = 0
    call umf4sym(nvar, nvar, ap, ai, ax, symbolic, control, info)
    if (info(1) .lt. 0) then
      print *, 'error occurred in umf4sym: ', info(1)
      stop
    endif

    call umf4num(ap, ai, ax, symbolic, numeric, control, info)
    if (info(1) .lt. 0) then
      print *, 'error occurred in umf4num: ', info(1)
      stop
    endif

    ising = 0
    end subroutine umf_decomp

    subroutine umf_solve(rhs)
    ! solve ax=b, without iterative refinement
    integer :: sys
    double precision :: x(nvar),rhs(nvar)
    sys = 0
    call umf4sol(sys, x, rhs, numeric, control, info)

    if (info(1) .lt. 0) then
      print *, 'error occurred in umf4sol: ', info(1)
      stop
    endif
    ! free the numeric factorization, to be continued
    rhs(:) = x(:)
    end subroutine umf_solve
      

    subroutine umf_init(n,nn)
    integer :: n, nn
    integer :: state
    nvar = n
    nnz  = nn   
    allocate(ax(nnz),ai(nnz),ap(nvar+1),fjac(nvar,nvar),STAT = state)
    if(state .ne. 0) stop 'Allocation error in umf_init'
    symbolic = 0
    numeric = 0
    end subroutine umf_init
 
    subroutine umf_free
    integer :: state
    deallocate(ax,ai,ap,fjac,STAT=state)
    call umf4fsym(symbolic)
    call umf4fnum(numeric)
    if(state .ne. 0) stop 'Deallocation error in umf_free'
    end subroutine umf_free

end module umf_sdirk
#endif

#ifdef SPARSE_LU
module superlu_sdirk
    implicit none
    save
    integer :: nvar, nnz
    double precision, allocatable :: fjac(:,:)
! ax(nnz),ai(nnz),ap(nvar+1)
    double precision, allocatable :: ax(:), b(:)
    integer, allocatable :: ai(:), ap(:)
    integer :: info, iopt, ldb, nrhs, tr
    integer*8 factors
contains
    subroutine superlu_decomp(ising,hgamma)
    !convert the matrix e from coordinate format to column compressed format
    integer :: i,j,idx,found,ising
    double precision :: hgamma
!   free previouslu used memory
    nrhs = 1
    ldb = nvar
    tr = 0
    if(factors .ne. 0) then 
      iopt = 3
      call c_fortran_dgssv( tr, iopt, nvar, nnz, nrhs, ax, ai, ap,&
                       b, ldb, factors, info )
    end if

    idx = 1
    do j = 1,nvar
      found = 0
      do i = 1, nvar
        if ( fjac(i,j) .ne. 0d0) then
!          ax(idx) = e(i,j)
          if( i .eq. j ) then
            ax(idx) = -fjac(i,j) + hgamma
          else
            ax(idx) = -fjac(i,j)
          end if
          ! write the row in the row indices vector ai
          ai(idx) = i 
          if(found .eq. 0) then
            found = idx
          end if
          idx = idx + 1
        end if
      end do
      !must also update the column pointers vector ap
      ap(j) = found 
    end do
    !last element in ap must be nnza+1
    ap(nvar+1) = nnz + 1
!   factorize the matrix. The factors are stored in *factors* handle.
    iopt = 1
    call c_fortran_dgssv( tr, iopt, nvar, nnz, nrhs, ax, ai, ap,&
                          b, ldb, factors, info )
    if(info .ne. 0) then
      write(*,*) 'INFO from failed factorization = ', info
    endif
    ising = info
    end subroutine superlu_decomp

    subroutine superlu_solve(rhs)
    double precision :: rhs(nvar)
!   solve the system using the existing factors.
    iopt = 2
    call c_fortran_dgssv( tr, iopt, nvar, nnz, nrhs, ax, ai, ap,&
                          rhs, ldb, factors, info )
    if(info .ne. 0) then
      write(*,*) 'INFO from failed triangular solve = ', info
    endif

    end subroutine superlu_solve


    subroutine superlu_init(n,nn)
    integer :: n,nn
    integer :: state
    nvar = n
    nnz  = nn
    factors = 0
    allocate(ax(nnz),ai(nnz),ap(nvar+1),b(nvar),fjac(nvar,nvar),STAT=state)
    if(state .ne. 0) stop 'Allocation error in superlu_init'
    end subroutine superlu_init

    subroutine superlu_free
    integer :: state
    deallocate(ax,ai,ap,b,fjac,STAT=state)
    if(state .ne. 0) stop 'Deallocation error in superlu_free'
    iopt = 3
    call c_fortran_dgssv( tr, iopt, nvar, nnz, nrhs, ax, ai, ap,&
                       b, ldb, factors, info )
    end subroutine superlu_free

end module superlu_sdirk
#endif
