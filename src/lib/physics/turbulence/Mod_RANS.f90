module MOSE_Mod_RANS
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  
contains


  subroutine Set_RANS_Model ()
    use MOSE_Global_m
    use MOSE_IO_Table, only: mcdiff, schmidt_given, schmidt_val, prandtl_val
    use MOSE_IO_RANS
    use MOSE_Lib_RANS
    use MOSE_Mod_Fluid
    use MOSE_Mod_Spalart
    use MOSE_Mod_SSGLRR
    use MOSE_Mod_SST
    use MOSE_Mod_Wilcox2006
    use MOSE_Mod_QCR

    implicit none

    ! QCR2000, the one correction without a device twin: refused for every model
    if ( obj_rans%QCR2000 )     error stop 'turbulence-model with QCR2000 is not on the GPU'

    ! Setting RANS or NS model
    if ( index ( trim(obj_rans%rans_name), 'laminar' ) > 0 ) then
      if (0 /= nrans) error stop 'turbulence model: this binary was built for another number of RANS variables (one binary per model)'
      VISC = .true.
      RANS = .false.
      Eddy_Viscosity => null()
      RANS_Diffusive_Flux => null()
      Stress_Vector => Stress_Vector_Std
      stress_dev = 1   ! laminar -> Std (device build: viscous, no turbulence)
      RANS_Enforce_Realizability => null()

    elseif ( index ( trim(obj_rans%rans_name), 'SA' ) > 0 ) then
      if (1 /= nrans) error stop 'turbulence model: this binary was built for another number of RANS variables (one binary per model)'
      VISC = .true.
      RANS = .true.
      Eddy_Viscosity => Spalart_Eddy_Viscosity
      RANS_Diffusive_Flux => Spalart_RANS_Diffusive_Flux
      if ( obj_rans%QCR2000 ) then
        Stress_Vector => Stress_Vector_QCR2000
      else
        Stress_Vector => Stress_Vector_Std
      end if
      RANS_Source_Terms => Spalart_Source_Terms
      RANS_Set_Wall_Values => Spalart_Set_Wall_Values
      RANS_Set_Blowing_Wall => Spalart_Set_Blowing_Wall
      RANS_Extrapolate_Wall => Spalart_Extrapolate_Wall
      RANS_Enforce_Realizability => Spalart_Enforce_Realizability
      eddy_dev = 2; ransdiff_dev = 2   ! device SA: diff_physics.inc (eddy viscosity, nu-tilde flux), sa_body.inc (source)

    elseif ( index ( trim(obj_rans%rans_name), 'Wilcox2006' ) > 0 ) then
      if (2 /= nrans) error stop 'turbulence model: this binary was built for another number of RANS variables (one binary per model)'
      VISC = .true.
      RANS = .true.
      Eddy_Viscosity => Wilcox2006_Eddy_Viscosity
      RANS_Diffusive_Flux => Wilcox2006_RANS_Diffusive_Flux
      eddy_dev = 3; ransdiff_dev = 3   ! device Wilcox2006: diff_physics.inc (eddy viscosity, k/omega flux), wilcox_body.inc (source)
      if ( obj_rans%k_energy_coupling ) then
        Stress_Vector => Stress_Vector_2eq
        stress_dev = 2
      else
        Stress_Vector => Stress_Vector_Std
        stress_dev = 1
      end if
      RANS_Source_Terms => Wilcox2006_Source_Terms
      RANS_Set_Wall_Values => Wilcox2006_Set_Wall_Values
      if (obj_rans%blowing_corr) then
        RANS_Set_Blowing_Wall => SST_Blowing_Correction
      else
        RANS_Set_Blowing_Wall => SST_Blowing_noCorrection
      end if
      RANS_Extrapolate_Wall => Wilcox2006_Extrapolate_Wall
      RANS_Enforce_Realizability => Wilcox2006_Enforce_Realizability
    
    elseif ( index ( trim(obj_rans%rans_name), 'SST' ) > 0 ) then
      if (2 /= nrans) error stop 'turbulence model: this binary was built for another number of RANS variables (one binary per model)'
      VISC = .true.
      RANS = .true.
      Eddy_Viscosity => SST_Eddy_Viscosity
      RANS_Diffusive_Flux => SST_RANS_Diffusive_Flux
      eddy_dev = 1; ransdiff_dev = 1   ! device-supported SST diffusive
      if (obj_rans%k_energy_coupling) then
        Stress_Vector => Stress_Vector_2eq
        stress_dev = 2
      else
        Stress_Vector => Stress_Vector_Std
        stress_dev = 1
      end if
      RANS_Source_Terms => SST_Source_Terms
      RANS_Set_Wall_Values => SST_Set_Wall_Values
      if (obj_rans%blowing_corr) then
        RANS_Set_Blowing_Wall => SST_Blowing_Correction
      else
        RANS_Set_Blowing_Wall => SST_Blowing_noCorrection
      end if
      RANS_Extrapolate_Wall => SST_Extrapolate_Wall
      RANS_Enforce_Realizability => SST_Enforce_Realizability

    elseif ( index ( trim(obj_rans%rans_name), 'SSGLRR' ) > 0 ) then
      if (7 /= nrans) error stop 'turbulence model: this binary was built for another number of RANS variables (one binary per model)'
      VISC = .true.
      RANS = .true.
      obj_rans%RSM = .true.
      error stop 'turbulence-model SSGLRR is not on the GPU: use SST, Wilcox2006, SA, SA-R or laminar'
      Eddy_Viscosity => SSGLRR_Eddy_Viscosity
      if ( index ( trim(obj_rans%rans_name), '-SD' ) > 0 ) then
        RANS_Diffusive_Flux => SSGLRR_SD_RANS_Diffusive_Flux
      else
        RANS_Diffusive_Flux => SSGLRR_RANS_Diffusive_Flux
      end if
      Stress_Vector => Stress_Vector_RSM
      RANS_Source_Terms => SSGLRR_Source_Terms
      RANS_Set_Wall_Values => SSGLRR_Set_Wall_Values
      RANS_Extrapolate_Wall => SSGLRR_Extrapolate_Wall
      RANS_Enforce_Realizability => SSGLRR_Enforce_Realizability

    else
      if (0 /= nrans) error stop 'turbulence model: this binary was built for another number of RANS variables (one binary per model)'
      VISC = .false.
      RANS = .false.
      Eddy_Viscosity => null()
      RANS_Diffusive_Flux => null()
      Stress_Vector => null()
      obj_rans%rans_name = '  Inviscid flow'

    endif

    if ( RANS .and. obj_rans%point_implicit ) write(*,'(A)') ' [RANS] point-implicit treatment of the turbulence destruction terms (Patankar), as MOSE_open'

    ! Push the devirtualized turbulence ids and the diffusion non-dim numbers to
    ! the device. `acc declare create` only reserves device storage — the values
    ! set above (and read from obj_rans, which lives only on the host) must be
    ! synced once here, after the model is fully selected.
    ! DIFF_MC: [MOSE-Physics] schmidt > 0 replaces [MOSE-RANS] Sc (reference semantics); <= 0 selects multicomponent D_k.
    if ( schmidt_given .and. schmidt_val > 0d0 ) then
      write(*,'(A,F8.4,A,F8.4,A)') ' [DIFF_MC] [MOSE-Physics] schmidt = ', schmidt_val, ' replaces [MOSE-RANS] Sc = ', &
                                   obj_rans%Sc, ' (constant laminar Schmidt)'
      obj_rans%Sc = schmidt_val
    end if
    mcdiff_dev = mcdiff
    if ( mcdiff .and. .not. VISC ) write(*,'(A)') ' [DIFF_MC] note: schmidt <= 0 has no effect on an inviscid run'
    sc_dev = obj_rans%Sc; sct_dev = obj_rans%Sct; prt_dev = obj_rans%Prt; prl_dev = prandtl_val
    if ( RANS .and. eddy_dev == 0 ) error stop 'turbulence model without a device twin (eddy_dev = 0): SST, Wilcox2006, SA and SA-R run on the GPU'
    !$acc update device(eddy_dev, ransdiff_dev, stress_dev, sc_dev, sct_dev, prt_dev, prl_dev, RANS, VISC, mcdiff_dev)

  end subroutine Set_RANS_Model


end module MOSE_Mod_RANS