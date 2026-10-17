module MOSE_Mod_Probes
  use MOSE_IO_Ini_Open, only: MOSE_Load_Ini   ! one ini loader: MOSE_open dialect translated in memory
  use MOSE_Types_m
  use MOSE_Global_m
  use MOSE_General_Types_m, only: MOSE_vector_3D_type
  use MOSE_Parameters_m

  implicit none

  type :: real_ptr
    real(8), pointer :: p
  end type real_ptr

  type :: obj_probe
    type(MOSE_vector_3D_type)        :: location
    integer, dimension(4)            :: ilocation
    character(len=llen)              :: file
    integer                          :: nvar
    character(len=clen), allocatable :: names(:)
    type(real_ptr), allocatable      :: variables(:)
    real(8), pointer                 :: P(:)
    real(8)                          :: T, M
    real(8)                          :: dtime
    integer                          :: ntime
    integer                          :: diter
    integer                          :: unit
  contains
    private 
    procedure, pass(self) :: Place
    procedure, pass(self) :: Update
  end type obj_probe

  integer :: nprobes
  type(obj_probe), allocatable, target :: probe(:)

contains


  subroutine Setup_Probes( domain )
    use IR_precision
    use Finer, only: file_ini
    use strings, only: parse
    use MOSE_Mod_MPI, only: is_local_block

    implicit none
    type(MOSE_domain_type), intent(in)  :: domain
    ! Local
    type(file_ini)  :: fini
    integer         :: i, error
    character(llen) :: name, wholestring
    logical         :: restart


    call MOSE_Load_Ini(fini)
    ! Count the probes
    i = 0
    do
      i = i+1
      call fini%get(section_name=trim(codename)//'-probes', option_name='probe'//trim(str(.true.,i)), val=name, error=error)
      if (error/=0) exit
    enddo

    ! Allow both restart and newrun. If restart=true (newrun=false) append data to the old probes files
    restart = .false.
    ! Try restart
    call fini%get(section_name=trim(codename)//'-probes', option_name='restart', val=restart, error=error)
    ! If not found try newrun
    if (error/=0) then
      call fini%get(section_name=trim(codename)//'-probes', option_name='newrun', val=restart, error=error)
      if (error==0) restart = .not.restart
    endif

    ! Allocate the probes
    nprobes = i-1
    allocate(probe(nprobes))

    do i = 1, nprobes
      call fini%get(section_name=trim(codename)//'-probes', option_name='probe'//trim(str(.true.,i)), val=name, error=error)
      probe(i)%file = 'OUTPUT/'//trim(name)//'.txt'

      ! Vars name
      probe(i)%nvar = fini%count_values(section_name=name, option_name='variables')
      allocate(probe(i)%names(1:probe(i)%nvar))
      call fini%get(section_name=name, option_name='variables', val=wholestring, error=error)
      call parse(wholestring,' ',probe(i)%names)

      ! Frequency
      probe(i)%dtime = huge(1d0); probe(i)%diter = huge(1); probe(i)%ntime = 1
      call fini%get(section_name=name, option_name='dtime', val=probe(i)%dtime, error=error)
      call fini%get(section_name=name, option_name='diter', val=probe(i)%diter, error=error)

      ! Location
      probe(i)%location%c = 0.d0
      probe(i)%ilocation = 1
      call fini%get(section_name=name, option_name='index-position', val=probe(i)%ilocation, error=error)
      if (error/=0) then
        call fini%get(section_name=name, option_name='position', val=probe(i)%location%c, error=error)
        if (error/=0) then
          write(*,*) ' Missing probe position input'; stop
        endif
        call probe(i)%Place(domain)
      endif

      ! Only the rank owning the probe's block sets up pointers and opens the file
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle

      call Assign_Variables(probe(i), domain)

      if (restart) then
        open(newunit=probe(i)%unit,file=trim(probe(i)%file),status='OLD',iostat=error)
        if (error/=0) stop "[ERROR] You restarted from an old solution but the probes files were not found."
        error = 0
        do while ( error == 0 ); read (probe(i)%unit, *, iostat=error); enddo
        backspace (probe(i)%unit)
      else
        open(newunit=probe(i)%unit,file=trim(probe(i)%file))
      endif

    enddo

  end subroutine Setup_Probes


  !> Find the cell indexes related to the probe starting from space coords
  subroutine Place (self, domain)

    implicit none
    class(obj_probe)                   :: self
    type(MOSE_domain_type), intent(in) :: domain
    ! Local
    integer :: i, j, k, b
    real(8) :: d0, d


    d0 = huge(1d0)
    do b = 1, domain%nb
      do k = 0, domain%blk(b)%dim(3); do j = 0, domain%blk(b)%dim(2); do i = 1, domain%blk(b)%dim(1)
        d = norm2(self%location%c-domain%blk(b)%node(i,j,k)%c)
        if (d<d0) then
          d0 = d
          self%ilocation = [b , i, j, k]
        endif
      enddo; enddo; enddo
    enddo

  end subroutine Place


  subroutine Assign_Variables(probe, domain)
    use IR_precision

    implicit none
    type(obj_probe), intent(inout), target     :: probe
    type(MOSE_domain_type), intent(in), target :: domain
    ! Local
    integer :: v,s,b,i,j,k


    b = probe%ilocation(1)
    i = probe%ilocation(2)
    j = probe%ilocation(3)
    k = probe%ilocation(4)

    probe%P => domain%blk(b)%P(:,i,j,k)
    allocate(probe%variables(1:probe%nvar))

    do v = 1, probe%nvar
      do s = 1, nsc
        if (probe%names(v)=='rho('//trim(str(.true.,s))//')') probe%variables(v)%p => probe%P(s)
      enddo
      if (probe%names(v)=='u') probe%variables(v)%p => probe%P(nu)
      if (probe%names(v)=='v') probe%variables(v)%p => probe%P(nv)
      if (probe%names(v)=='w') probe%variables(v)%p => probe%P(nw)
      if (probe%names(v)=='p') probe%variables(v)%p => probe%P(np)
      if (probe%names(v)=='T') probe%variables(v)%p => probe%T
      if (probe%names(v)=='M') probe%variables(v)%p => probe%M
    enddo

  end subroutine Assign_Variables


  subroutine Update(self)
    use FLINT_Lib_Thermodynamic
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none
    class(obj_probe), intent(inout) :: self
    ! Local
    real(8) :: Rtot, speed_of_sound, rho, T_est

    ! Guard the thermo-table lookups: the host f_tabT path has NO off-table
    ! clamp (the device twins do), so a probed cell whose state is transiently
    ! NaN / off-table turns int(T) into a wild table index -> SIGSEGV. That
    ! killed production job 296241 (Tw300 leg, 2026-09-01) from inside probe
    ! post-processing. On a bad state freeze the derived quantities instead;
    ! the probed primitive (e.g. p) is still written by the caller.
    rho  = sum(self%P(1:nsc))
    Rtot = f_Rtot( self%P(1:nsc) )
    T_est = self%P(np) / (rho * Rtot)
    if ( .not. ieee_is_finite(T_est) .or. &
         T_est <= real(Tmin,8) .or. T_est >= real(Tmax,8) ) then
      self%T = T_est
      self%M = -1d0
      return
    end if

    speed_of_sound = f_ss(self%P(1:nsc),self%P(np),rho,Rtot)
    self%T = EOS(p=self%P(np),ro=rho,R=Rtot)
    self%M = norm2(self%P(nu:nw))/speed_of_sound

  end subroutine Update


  !> With the device residual, residual-only iterations no longer refresh
  !> the whole host P, so the blk%P columns the probe pointers alias would go STALE between
  !> I/O iterations. Refresh exactly the probed cells from the resident var-last w (strided
  !> micro-D2H, ~nprim tiny copies per probe -- negligible). This also fixes the pre-existing
  !> probe staleness whenever diter-res is set high: probes now always sample the current w.
  !> w interior is current at end-of-iteration (Newstate/chemistry write w); probes sample
  !> interior cells (Place clamps to interior indices).
  subroutine Probes_Refresh_From_flat()
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Field_m, only: P_flat, field_block_index_params, field_is_allocated

    implicit none
    integer :: n, b, i, j, k, c0, boff, pi, pj

    if (.not. field_is_allocated()) return   ! pre-solve call: host P still holds the IC
    do n = 1, nprobes
      if (.not. is_local_block(probe(n)%ilocation(1))) cycle
      b = probe(n)%ilocation(1)
      i = probe(n)%ilocation(2)
      j = probe(n)%ilocation(3)
      k = probe(n)%ilocation(4)
      call field_block_index_params(b, boff, pi, pj)
      c0 = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      !$acc update host(P_flat(c0:c0, 1:nprim))
      probe(n)%P(1:nprim) = P_flat(c0, 1:nprim)
    end do

  end subroutine Probes_Refresh_From_flat


  subroutine Write_Probes_Data( iter, time )
    use MOSE_Mod_MPI, only: is_local_block

    implicit none
    integer, intent(in) :: iter
    real(8), intent(in) :: time
    ! Local
    integer :: i, v


    call Probes_Refresh_From_flat()
    do i = 1, nprobes
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle
      if (mod(iter,probe(i)%diter)==0) then
        call probe(i)%update
        write(probe(i)%unit,*) iter, (probe(i)%variables(v)%p,v=1,probe(i)%nvar)
      elseif (time >= probe(i)%dtime*probe(i)%ntime) then
        call probe(i)%update
        write(probe(i)%unit,*) time, (probe(i)%variables(v)%p,v=1,probe(i)%nvar)
        probe(i)%ntime = probe(i)%ntime+1
      endif
    enddo

  end subroutine Write_Probes_Data


end module MOSE_Mod_Probes
