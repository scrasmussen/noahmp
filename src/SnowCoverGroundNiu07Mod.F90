module SnowCoverGroundNiu07Mod

!!! Compute ground snow cover fraction based on Niu and Yang (2007, JGR) scheme

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  ! Measured on CCE 19.0.0 with -h acc -target-accel=nvidia80, which device
  ! intrinsics actually LINK inside !$acc routine seq:
  !
  !     sqrt              LINKS  (maps to a hardware instruction)
  !     x**y  const expo  links  -- only because -O2 constant-folds it
  !     x**y  runtime exp FAILS  (_HEXP, _HLOG)
  !     exp / log / log10 FAILS  (_HEXP, _HLOG, _HLOG10)
  !     tanh / tan / atan / acos / cos   FAILS
  !
  ! So sqrt is the only intrinsic that can replace a shim here. tanh and pow
  ! both have to stay on the shim: SnowMeltFac is a runtime value, so `**`
  ! below would emit _HEXP/_HLOG and fail at nvlink.
  !
  ! NOTE both shims are crude. acc_powf is exp(y*log(x)) through two more
  ! approximations (and returns 0 for x<=0); acc_tanhf is (e-1)/(e+1) with
  ! e=acc_expf(2x), which cancels badly for small x -- exactly the
  ! shallow-snow regime FSNO is most sensitive to. These remain the prime
  ! suspects for the ~392k one-signed FSNO differences vs the CPU, and the
  ! fix is to make the shims accurate rather than to remove them.
#ifdef NOAHMP_ACC_COLUMNS
  use NoahmpAccDeviceMathShimMod, only : tanh => acc_tanhf, pow => acc_powf
#else
  ! pow is a C-ism supplied only by the shim; on the host path it
  ! comes from NoahmpMathHostMod, where it is the exact x**y.
  use NoahmpMathHostMod, only : pow
#endif

  implicit none

contains

  subroutine SnowCoverGroundNiu07(noahmp)
#ifdef NOAHMP_ACC_COLUMNS
!$acc routine seq
#endif

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: None (embedded in ENERGY subroutine)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    real(kind=kind_noahmp)           :: SnowDensBulk   ! bulk density of snow [Kg/m3]
    real(kind=kind_noahmp)           :: MeltFac        ! melting factor for snow cover frac

! --------------------------------------------------------------------
    associate(                                                     &
              SnowMeltFac     => noahmp%water%param%SnowMeltFac     ,& ! in,  snowmelt m parameter
              SnowCoverFac    => noahmp%water%param%SnowCoverFac    ,& ! in,  snow cover factor [m]
              SnowCoverFracMax=> noahmp%water%param%SnowCoverFracMax,& ! in,  legacy SCAMAX cap (only used when NOAHMP_LEGACY_PHYSICS)
              SnowDepth       => noahmp%water%state%SnowDepth       ,& ! in,  snow depth [m]
              SnowWaterEquiv  => noahmp%water%state%SnowWaterEquiv  ,& ! in,  snow water equivalent [mm]
              SnowCoverFrac   => noahmp%water%state%SnowCoverFrac    & ! out, snow cover fraction
             )
! ----------------------------------------------------------------------

    SnowCoverFrac = 0.0
    if ( SnowDepth > 0.0 ) then
         SnowDensBulk  = SnowWaterEquiv / SnowDepth
         MeltFac       = pow(SnowDensBulk / 100.0, SnowMeltFac)
        !SnowCoverFrac = tanh( SnowDepth /(2.5 * Z0 * MeltFac))
         SnowCoverFrac = tanh( SnowDepth /(SnowCoverFac * MeltFac)) ! C.He: bring hard-coded 2.5*z0 to MPTABLE
#ifdef NOAHMP_LEGACY_PHYSICS
         ! Legacy NoahMP capped FSNO at parameters%SCAMAX (a per-cell field from
         ! SPATIAL_SOIL); the He-et-al refactor dropped that cap. Restore it for
         ! bit-reproduction. Bridge sets SnowCoverFracMax = SCAMAX_2D(I,J).
         SnowCoverFrac = SnowCoverFracMax * SnowCoverFrac
#endif
    endif

    end associate

  end subroutine SnowCoverGroundNiu07

end module SnowCoverGroundNiu07Mod
