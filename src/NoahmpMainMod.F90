module NoahmpMainMod

!!! Main NoahMP module including all column model processes
!!! atmos forcing -> canopy intercept -> precip heat advect -> main energy -> main water -> main biogeochemistry -> balance check

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use AtmosForcingMod,            only : ProcessAtmosForcing
  use GeneralInitMod,             only : GeneralInit
  use PhenologyMainMod,           only : PhenologyMain
  use IrrigationPrepareMod,       only : IrrigationPrepare
  use IrrigationSprinklerMod,     only : IrrigationSprinkler
  use CanopyWaterInterceptMod,    only : CanopyWaterIntercept
  use PrecipitationHeatAdvectMod, only : PrecipitationHeatAdvect
  use EnergyMainMod,              only : EnergyMain
  use WaterMainMod,               only : WaterMain
  use BiochemNatureVegMainMod,    only : BiochemNatureVegMain
  use BiochemCropMainMod,         only : BiochemCropMain
  use BalanceErrorCheckMod,       only : BalanceWaterInit, BalanceWaterCheck, BalanceEnergyCheck 
 
  implicit none

contains

  subroutine NoahmpMain(noahmp)
#ifdef NOAHMP_ACC_COLUMNS
!$acc routine seq
#endif

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: NOAHMP_SFLX
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! --------------------------------------------------------------------
    associate(                                                                     &
              FlagDynamicVeg         => noahmp%config%domain%FlagDynamicVeg       ,& ! in,    flag to activate dynamic vegetation model
              FlagDynamicCrop        => noahmp%config%domain%FlagDynamicCrop      ,& ! in,    flag to activate dynamic crop model
              OptCropModel           => noahmp%config%nmlist%OptCropModel         ,& ! in,    option for crop model
              IrrigationAmtSprinkler => noahmp%water%state%IrrigationAmtSprinkler ,& ! inout, irrigation water amount [m] for sprinkler
              FlagCropland           => noahmp%config%domain%FlagCropland          & ! out,   flag to identify croplands
             )
! ----------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! Atmospheric forcing processing
    !--------------------------------------------------------------------- 

    call ProcessAtmosForcing(noahmp)

    !---------------------------------------------------------------------
    ! General initialization to prepare key variables
    !--------------------------------------------------------------------- 

    call GeneralInit(noahmp)

    !---------------------------------------------------------------------
    ! Prepare for water balance check
    !--------------------------------------------------------------------- 

#ifndef NOAHMP_ACC_COLUMNS
    call BalanceWaterInit(noahmp)
#endif

    !---------------------------------------------------------------------
    ! Phenology
    !--------------------------------------------------------------------- 

    call PhenologyMain(noahmp)

    !---------------------------------------------------------------------
    ! Irrigation prepare including trigger
    !--------------------------------------------------------------------- 

#ifndef NOAHMP_ACC_COLUMNS
    call IrrigationPrepare(noahmp)

    !---------------------------------------------------------------------
    ! Sprinkler irrigation
    !--------------------------------------------------------------------- 

    ! call sprinkler irrigation before canopy process to have canopy interception
    if ( (FlagCropland .eqv. .true.) .and. (IrrigationAmtSprinkler > 0.0) ) &
       call IrrigationSprinkler(noahmp)
#endif

    !---------------------------------------------------------------------
    ! Canopy water interception and precip heat advection
    !--------------------------------------------------------------------- 

    call CanopyWaterIntercept(noahmp)
    call PrecipitationHeatAdvect(noahmp)

    !---------------------------------------------------------------------
    ! Energy processes
    !--------------------------------------------------------------------- 

    call EnergyMain(noahmp)

    !---------------------------------------------------------------------
    ! Water processes
    !--------------------------------------------------------------------- 

    call WaterMain(noahmp)

    !---------------------------------------------------------------------
    ! Biochem processes (crop and carbon)
    !--------------------------------------------------------------------- 

    ! for generic vegetation
    if ( FlagDynamicVeg .eqv. .true. ) then
#ifndef NOAHMP_ACC_COLUMNS
       call BiochemNatureVegMain(noahmp)
#endif
    endif
   
    ! for explicit crop treatment
    if ( (OptCropModel == 1) .and. (FlagDynamicCrop .eqv. .true.) ) then
#ifndef NOAHMP_ACC_COLUMNS
       call BiochemCropMain(noahmp)
#endif
    endif

    !---------------------------------------------------------------------
    ! Error check for energy and water balance
    !--------------------------------------------------------------------- 

#ifndef NOAHMP_ACC_COLUMNS
    call BalanceWaterCheck(noahmp)
    call BalanceEnergyCheck(noahmp) 
#endif

    !---------------------------------------------------------------------
    ! End of all NoahMP column processes
    !--------------------------------------------------------------------- 

    end associate

  end subroutine NoahmpMain

end module NoahmpMainMod
