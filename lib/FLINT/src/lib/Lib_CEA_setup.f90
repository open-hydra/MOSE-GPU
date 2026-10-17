module FLINT_CEA_setup

  real(kind=8), dimension(:), allocatable   :: el_weight
contains

  subroutine CEA_initialize_global()
    use FLINT_CEA_data
    use FLINT_Lib_Thermodynamic
    implicit none
    integer :: i, k, satomic

    Npt=1
    Nonly = 0
    Nomit = 0
    Trace = 0
    Short = .false.
    Massf = .false.
    Np = 1
    Trace = 0.
    Lsave = 0
    R = Rr/4184.D0
    R = Rr/1000.0
    Tp = .false.
    HP = .true.
    Vol = .true.
    Ions = .false.
    Nt = 1
    Size = 0.

    ng = ns
    nc = 0
    ngc = ng + nc
    Ngp1 = Ng + 1
    Nreac = ns
    nlm = ne
    prod(1:ns) = species_names
    Elmt(1:ne) = elements_names

    Mw = 0
    Mw(1:ns) = wm_tab

    !call build_species_composition(species_names)
    call assign_elemental_weight

    Jx = 0
    Deln = 0.d0
    A = 0.d0
    satomic = 0

    do i = 1, ngc
      do k = 1, ne
        A(k, i) = species_composition(k,i)
      end do

      ! Post-processing for atomic species
      if (nc == 0) then
        if (sum(species_composition(:,i)) == 1.0) then
          satomic = satomic + 1
          do k = 1, ne
            if (Elmt(k) == species_names(i)) then
              Jx(k) = i
              Jcm(k) = i
              exit
            end if
          end do
        end if
      end if
    end do

    ! FIND MISSING ELEMENTS (IF ANY) FOR COMPONENTS
    Nspx = Ngc
    if ( satomic<ne ) then
      do k = 1, ne
        if ( Jx(k)==0 ) then
          Nspx = Nspx + 1
          A(k,Nspx) = 1.
          Jx(k) = Nspx
          Jcm(k) = Nspx
        endif
      enddo
    endif

    Size = 18.420681D0

  end subroutine CEA_initialize_global


  subroutine CEA_initialize_local(Tt, Pecwt, Vv, B0P, Hsub0)
    use FLINT_CEA_data
    use FLINT_Lib_Thermodynamic
    implicit none
    real(8), dimension(:), intent(in) :: Pecwt
    real(8), intent(in) :: Tt
    real(8), intent(out) :: Vv
    real(8), intent(out) :: B0p(:,:)
    real(8), intent(out) :: Hsub0
    integer :: j, kr, n, jj
    real(8) :: bigb, dbi(nlm)
    real(8) :: Enth(ns)
    real(8) :: dat(ne)
    logical :: gaseous
    real(8) :: rm, wp, hpp

    Vv = 1d05/sum(Pecwt)

    ! IF OXIDANT, KR = 1
    ! IF FUEL, KR = 2
    kr = 2

    Wp = 0.d0
    B0p(:,2) = 0.d0
    hpp = 0.d0

    do n = 1, Nreac

      ! STORE ATOMIC SYMBOLS IN ELMT ARRAY.
      ! CALCULATE MOLECULAR WEIGHT.
      rm = 0.D0
      dat = 0
      do jj = 1, ne
        rm = rm + species_composition(jj,n)*el_weight(jj)
        dat(jj) = dat(jj) + species_composition(jj,n)
      enddo

      if (Pecwt(n)<=1d-15) cycle

      if (index('(L)',species_names(n))>0) then
        gaseous = .false.
      else
        gaseous = .true.
      end if
      ! ADD CONTRIBUTIONS TO WP, HPP, AND B0P(I,K)
      Wp = Wp + Pecwt(n)

      Enth(n) = H0i(n,Tt)* wm_tab(n) / Rr
      if ( Vol .and. gaseous ) Enth(n) = Enth(n) - Tt
      Hpp = Hpp + Enth(n)*Pecwt(n)/rm

      do j = 1,Nlm
        B0p(j,kr) = dat(j)*Pecwt(n)/rm + B0p(j,kr)
      enddo
    enddo

    Hsub0 = hpp

    if ( Wp/=0.d0 ) then
      do j = 1,Nlm
        B0p(j,kr) = B0p(j,kr)/Wp
      enddo
    endif

    dbi = DABS(b0p(1:nlm,2))
    bigb = maxval(dbi)
    Bcheck = bigb*.000001D0

  end subroutine CEA_initialize_local

  subroutine assign_elemental_weight()
    use FLINT_Lib_Thermodynamic
    use FLINT_CEA_data
    implicit none
    integer :: i, j

    allocate(el_weight(1:ne))
    
    do j = 1, ne
      do i = 1, 100
        if (trim(Symbol(i))==trim(elements_names(j))) el_weight(j) = Atmwt(i)
      enddo
    enddo

   end subroutine assign_elemental_weight

end module FLINT_CEA_setup