module RunoffSurfaceXinAnJiangMod

!!! Compute surface infiltration rate and surface runoff based on XinAnJiang runoff scheme
!!! Reference: Knoben, W. J., et al., (2019): Modular Assessment of Rainfall-Runoff Models 
!!! Toolbox (MARRMoT) v1.2 an open-source, extendable framework providing implementations 
!!! of 46 conceptual hydrologic models as continuous state-space formulations.

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  ! The five x**y below all have RUNTIME exponents (shape parameters), and
  ! CCE 19.0.0 cannot link x**y with a runtime exponent inside
  ! !$acc routine seq -- it emits _HEXP/_HLOG, which nvlink cannot resolve.
  ! Route them through the device shim, as SoilHydraulicPropertyMod already
  ! does. Every base here is >= 0 by construction (SoilWaterTmp and
  ! SoilWaterFree are min()-clamped to their maxima, and the branch guarding
  ! the 0.5-TensionWatDistrInfl base requires it to be positive), so
  ! acc_powf's x<=0 -> 0 behaviour is never reached in a way that differs
  ! from Fortran's 0**y = 0.
#ifdef NOAHMP_ACC_COLUMNS
  use NoahmpAccDeviceMathShimMod, only : pow => acc_powf
#else
  ! pow is a C-ism supplied only by the shim; on the host path it
  ! comes from NoahmpMathHostMod, where it is the exact x**y.
  use NoahmpMathHostMod, only : pow
#endif

  implicit none

contains

  subroutine RunoffSurfaceXinAnJiang(noahmp, TimeStep)
#ifdef NOAHMP_ACC_COLUMNS
!$acc routine seq
#endif

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: COMPUTE_XAJ_SURFRUNOFF
! Original code: Prasanth Valayamkunnath <prasanth@ucar.edu>
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! ----------------------------------------------------------------------------------------

    implicit none

! IN & OUT variables
    type(noahmp_type)     , intent(inout) :: noahmp
    real(kind=kind_noahmp), intent(in)    :: TimeStep            ! timestep (may not be the same as model timestep)

! local variable
    integer                               :: LoopInd             ! do-loop index
    real(kind=kind_noahmp)                :: SoilWaterTmp        ! temporary soil water [m]
    real(kind=kind_noahmp)                :: SoilWaterMax        ! maximum soil water [m]
    real(kind=kind_noahmp)                :: SoilWaterFree       ! free soil water [m]
    real(kind=kind_noahmp)                :: SoilWaterFreeMax    ! maximum free soil water [m]
    real(kind=kind_noahmp)                :: RunoffSfcImp        ! impervious surface runoff [m]
    real(kind=kind_noahmp)                :: RunoffSfcPerv       ! pervious surface runoff [m]

! --------------------------------------------------------------------
    associate(                                                                 &
              NumSoilLayer         => noahmp%config%domain%NumSoilLayer       ,& ! in,  number of soil layers
              DepthSoilLayer       => noahmp%config%domain%DepthSoilLayer     ,& ! in,  depth [m] of layer-bottom from soil surface
              SoilMoisture         => noahmp%water%state%SoilMoisture         ,& ! in,  total soil moisture [m3/m3]
              SoilImpervFrac       => noahmp%water%state%SoilImpervFrac       ,& ! in,  fraction of imperviousness due to frozen soil
              SoilSfcInflowMean    => noahmp%water%flux%SoilSfcInflowMean     ,& ! in,  mean water input on soil surface [m/s]
              SoilMoistureSat      => noahmp%water%param%SoilMoistureSat      ,& ! in,  saturated value of soil moisture [m3/m3]
              SoilMoistureFieldCap => noahmp%water%param%SoilMoistureFieldCap ,& ! in,  reference soil moisture (field capacity) [m3/m3]
              TensionWatDistrInfl  => noahmp%water%param%TensionWatDistrInfl  ,& ! in,  Tension water distribution inflection parameter
              TensionWatDistrShp   => noahmp%water%param%TensionWatDistrShp   ,& ! in,  Tension water distribution shape parameter
              FreeWatDistrShp      => noahmp%water%param%FreeWatDistrShp      ,& ! in,  Free water distribution shape parameter
              RunoffSurface        => noahmp%water%flux%RunoffSurface         ,& ! out, surface runoff [m/s]
              InfilRateSfc         => noahmp%water%flux%InfilRateSfc           & ! out, infiltration rate at surface [m/s]
             )
! ----------------------------------------------------------------------

    ! initialization 
    SoilWaterTmp     = 0.0
    SoilWaterMax     = 0.0
    SoilWaterFree    = 0.0
    SoilWaterFreeMax = 0.0
    RunoffSfcImp     = 0.0
    RunoffSfcPerv    = 0.0
    RunoffSurface    = 0.0
    InfilRateSfc     = 0.0

    do LoopInd = 1, NumSoilLayer-2
       if ( (SoilMoisture(LoopInd)-SoilMoistureFieldCap(LoopInd)) > 0.0 ) then   ! soil moisture greater than field capacity
          SoilWaterFree = SoilWaterFree + &
                          (SoilMoisture(LoopInd)-SoilMoistureFieldCap(LoopInd)) * (-1.0) * DepthSoilLayer(LoopInd)
          SoilWaterTmp  = SoilWaterTmp + SoilMoistureFieldCap(LoopInd) * (-1.0) * DepthSoilLayer(LoopInd)
       else
          SoilWaterTmp  = SoilWaterTmp + SoilMoisture(LoopInd) * (-1.0) * DepthSoilLayer(LoopInd)
       endif
       SoilWaterMax     = SoilWaterMax + SoilMoistureFieldCap(LoopInd) * (-1.0) * DepthSoilLayer(LoopInd)
       SoilWaterFreeMax = SoilWaterFreeMax + &
                          (SoilMoistureSat(LoopInd)-SoilMoistureFieldCap(LoopInd)) * (-1.0) * DepthSoilLayer(LoopInd)
    enddo
    SoilWaterTmp  = min(SoilWaterTmp, SoilWaterMax)      ! tension water [m] 
    SoilWaterFree = min(SoilWaterFree, SoilWaterFreeMax) ! free water [m]

    ! impervious surface runoff R_IMP    
    RunoffSfcImp = SoilImpervFrac(1) * SoilSfcInflowMean * TimeStep

    ! solve pervious surface runoff (m) based on Eq. (310)
    if ( (SoilWaterTmp/SoilWaterMax) <= (0.5-TensionWatDistrInfl) ) then
       RunoffSfcPerv = (1.0-SoilImpervFrac(1)) * SoilSfcInflowMean * TimeStep * &
                       pow(0.5-TensionWatDistrInfl, 1.0-TensionWatDistrShp) * &
                       pow(SoilWaterTmp/SoilWaterMax, TensionWatDistrShp)
    else
       RunoffSfcPerv = (1.0-SoilImpervFrac(1)) * SoilSfcInflowMean * TimeStep * &
                       (1.0-(pow(0.5+TensionWatDistrInfl, 1.0-TensionWatDistrShp) * &
                       pow(1.0-(SoilWaterTmp/SoilWaterMax), TensionWatDistrShp)))
    endif

    ! estimate surface runoff based on Eq. (313)
    if ( SoilSfcInflowMean == 0.0 ) then
      RunoffSurface = 0.0
    else
      RunoffSurface = RunoffSfcPerv * (1.0-pow(1.0-(SoilWaterFree/SoilWaterFreeMax), FreeWatDistrShp)) + RunoffSfcImp
    endif
    RunoffSurface = RunoffSurface / TimeStep
    RunoffSurface = max(0.0,RunoffSurface)
    RunoffSurface = min(SoilSfcInflowMean, RunoffSurface)
    InfilRateSfc  = SoilSfcInflowMean - RunoffSurface

    end associate

  end subroutine RunoffSurfaceXinAnJiang

end module RunoffSurfaceXinAnJiangMod
