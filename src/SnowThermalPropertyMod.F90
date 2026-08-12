module SnowThermalPropertyMod

!!! Compute snowpack thermal conductivity and volumetric specific heat

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
#ifdef NOAHMP_ACC_COLUMNS
  use NoahmpAccDeviceMathShimMod, only : pow => acc_powf
#else
  ! pow is a C-ism supplied only by the shim; on the host path it
  ! comes from NoahmpMathHostMod, where it is the exact x**y.
  use NoahmpMathHostMod, only : pow
#endif

  implicit none

contains

  subroutine SnowThermalProperty(noahmp)
#ifdef NOAHMP_ACC_COLUMNS
!$acc routine seq
#endif

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: CSNOW
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: LoopInd                        ! loop index
#ifdef NOAHMP_ACC_COLUMNS
    real(kind=kind_noahmp), dimension(-NoahmpAccMaxSnowLayer+1:0) :: SnowDensBulk  ! bulk density of snow [kg/m3]
#else
    real(kind=kind_noahmp), allocatable, dimension(:) :: SnowDensBulk              ! bulk density of snow [kg/m3]
#endif

! --------------------------------------------------------------------
    associate(                                                                       &
              NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg        ,& ! in,  actual number of snow layers (negative)
              NumSnowLayerMax        => noahmp%config%domain%NumSnowLayerMax        ,& ! in,  maximum number of snow layers
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in,  thickness of snow/soil layers [m]
              OptSnowThermConduct    => noahmp%config%nmlist%OptSnowThermConduct    ,& ! in,  options for snow thermal conductivity schemes
              SnowIce                => noahmp%water%state%SnowIce                  ,& ! in,  snow layer ice [mm]
              SnowLiqWater           => noahmp%water%state%SnowLiqWater             ,& ! in,  snow layer liquid water [mm]
              SnowIceVol             => noahmp%water%state%SnowIceVol               ,& ! out, partial volume of snow ice [m3/m3]
              SnowLiqWaterVol        => noahmp%water%state%SnowLiqWaterVol          ,& ! out, partial volume of snow liquid water [m3/m3]
              SnowEffPorosity        => noahmp%water%state%SnowEffPorosity          ,& ! out, snow effective porosity [m3/m3]
              HeatCapacVolSnow       => noahmp%energy%state%HeatCapacVolSnow        ,& ! out, snow layer volumetric specific heat [J/m3/K]
              ThermConductSnow       => noahmp%energy%state%ThermConductSnow         & ! out, snow layer thermal conductivity [W/m/K]
             )
! ----------------------------------------------------------------------

    ! initialization
#ifndef NOAHMP_ACC_COLUMNS
    if (.not. allocated(SnowDensBulk)) allocate(SnowDensBulk(-NumSnowLayerMax+1:0))
#endif
    SnowDensBulk = 0.0

    !  effective porosity of snow
    do LoopInd = NumSnowLayerNeg+1, 0
       SnowIceVol(LoopInd)      = min(1.0, SnowIce(LoopInd)/(ThicknessSnowSoilLayer(LoopInd)*ConstDensityIce))
       SnowEffPorosity(LoopInd) = 1.0 - SnowIceVol(LoopInd)
       SnowLiqWaterVol(LoopInd) = min(SnowEffPorosity(LoopInd), &
                                      SnowLiqWater(LoopInd)/(ThicknessSnowSoilLayer(LoopInd)*ConstDensityWater))
    enddo

    ! thermal capacity of snow
    do LoopInd = NumSnowLayerNeg+1, 0
       SnowDensBulk(LoopInd)     = (SnowIce(LoopInd) + SnowLiqWater(LoopInd)) / ThicknessSnowSoilLayer(LoopInd)
       HeatCapacVolSnow(LoopInd) = ConstHeatCapacIce*SnowIceVol(LoopInd) + ConstHeatCapacWater*SnowLiqWaterVol(LoopInd)
      !HeatCapacVolSnow(LoopInd) = 0.525e06  ! constant
    enddo

    ! thermal conductivity of snow
    do LoopInd = NumSnowLayerNeg+1, 0
       if (OptSnowThermConduct == 1) &
          ThermConductSnow(LoopInd) = 3.2217e-6 * SnowDensBulk(LoopInd)*SnowDensBulk(LoopInd)     ! Stieglitz(yen,1965)
       if (OptSnowThermConduct == 2) &
          ThermConductSnow(LoopInd) = 2e-2 + 2.5e-6*SnowDensBulk(LoopInd)*SnowDensBulk(LoopInd)   ! Anderson, 1976
       if (OptSnowThermConduct == 3) &
          ThermConductSnow(LoopInd) = 0.35                                                        ! constant
       if (OptSnowThermConduct == 4) &
          ThermConductSnow(LoopInd) = 2.576e-6 * SnowDensBulk(LoopInd)*SnowDensBulk(LoopInd) + 0.074 ! Verseghy (1991)
       if (OptSnowThermConduct == 5) &
          ThermConductSnow(LoopInd) = 2.22 * pow(SnowDensBulk(LoopInd)/1000.0, 1.88)              ! Douvill(Yen, 1981)
    enddo

    ! deallocate local arrays to avoid memory leaks
#ifndef NOAHMP_ACC_COLUMNS
    deallocate(SnowDensBulk)
#endif

    end associate

  end subroutine SnowThermalProperty

end module SnowThermalPropertyMod
