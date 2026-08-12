# CMakeLists.txt for the NoahMP submodule library
# Builds the refactored NoahMP physics and selected HRLDAS driver modules
# as a static library for use by WRF-Hydro.

file(GLOB NOAHMP_SUBMODULE_PHYSICS_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.F90
)

file(GLOB NOAHMP_SUBMODULE_UTILITY_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/utility/*.F90
)

# Select specific driver files needed for WRF-Hydro integration.
# Excludes NoahmpDriverMainMod (depends on module_ra_gfdleta),
# NoahmpReadNamelistMod (standalone namelist reader),
# NoahmpInitMainMod/NoahmpSnowInitMod/NoahmpGroundwaterInitMod (standalone init),
# and GroundWaterMmfMod (WRF-Hydro uses its own groundwater module).
set(NOAHMP_SUBMODULE_DRIVER_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/drivers/hrldas/NoahmpIOVarType.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/drivers/hrldas/NoahmpIOVarInitMod.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/drivers/hrldas/NoahmpReadTableMod.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/drivers/hrldas/PedoTransferSR2006Mod.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/drivers/hrldas/ConfigVarInTransferMod.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/drivers/hrldas/ConfigVarOutTransferMod.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/drivers/hrldas/ForcingVarInTransferMod.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/drivers/hrldas/ForcingVarOutTransferMod.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/drivers/hrldas/EnergyVarInTransferMod.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/drivers/hrldas/EnergyVarOutTransferMod.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/drivers/hrldas/WaterVarInTransferMod.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/drivers/hrldas/WaterVarOutTransferMod.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/drivers/hrldas/BiochemVarInTransferMod.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/drivers/hrldas/BiochemVarOutTransferMod.F90
)

add_library(noahmp_submodule STATIC
    ${NOAHMP_SUBMODULE_PHYSICS_SOURCES}
    ${NOAHMP_SUBMODULE_UTILITY_SOURCES}
    ${NOAHMP_SUBMODULE_DRIVER_SOURCES}
)

# Optional legacy-physics build: restores the older canopy-energy formulation
# (no canopy heat capacity term in the N-R solver, and the leaf boundary-layer
# resistance floor of 20 s/m with no upper cap) so output matches the built-in
# NoahMP used by the reference WRF-Hydro build.
option(NOAHMP_LEGACY_PHYSICS "Use legacy NoahMP physics (match built-in WRF-Hydro NoahMP)" OFF)
if(NOAHMP_LEGACY_PHYSICS)
    target_compile_definitions(noahmp_submodule PUBLIC NOAHMP_LEGACY_PHYSICS)
endif()

# Optional double-precision build: promotes kind_noahmp from 4 to 8 in
# noahmp/utility/Machine.F90, so all real(kind=kind_noahmp) variables inside
# the submodule are computed at FP64. The bridge driver and the rest of
# WRF-Hydro remain at FP32; the bridge handles conversions at the boundary.
# Use this to test whether the FP32 reorder-noise threshold-flip residual
# disappears when the submodule's internal arithmetic is FP64.
option(NOAHMP_DOUBLE_PREC "Compute NoahMP submodule physics at FP64" OFF)
if(NOAHMP_DOUBLE_PREC)
    target_compile_definitions(noahmp_submodule PUBLIC DOUBLE_PREC)
endif()

# Experimental OpenACC path for the refactored NoahMP column physics.  Default
# to the top-level GPU build setting, while still allowing an explicit override
# while this path is being brought up.
set(NOAHMP_ACC_COLUMNS_DEFAULT OFF)
if(OPENACC_GPU)
    set(NOAHMP_ACC_COLUMNS_DEFAULT ON)
endif()
option(NOAHMP_ACC_COLUMNS "Batch refactored NoahMP columns for OpenACC offload" ${NOAHMP_ACC_COLUMNS_DEFAULT})
if(NOAHMP_ACC_COLUMNS)
    target_compile_definitions(noahmp_submodule PUBLIC NOAHMP_ACC_COLUMNS)
    if(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
        # The top-level GPU build adds -h acc globally.  For this experimental
        # column path, keep unrelated NoahMP files host-only so their diagnostics,
        # glacier paths, irrigation code, and unused math do not enter nvlink.
        target_compile_options(noahmp_submodule PRIVATE "SHELL:-h noacc")

        set(NOAHMP_ACC_COLUMN_SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/src/NoahmpAccDeviceMathShimMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/NoahmpMainMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosForcingMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/GeneralInitMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/PhenologyMainMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/CanopyWaterInterceptMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/PrecipitationHeatAdvectMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/EnergyMainMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/WaterMainMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SnowCoverGroundNiu07Mod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/GroundRoughnessPropertyMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/GroundThermalPropertyMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SnowThermalPropertyMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SoilThermalPropertyMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SurfaceAlbedoMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SnowAgingBatsMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SnowAlbedoBatsMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/GroundAlbedoMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/CanopyRadiationTwoStreamMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SurfaceRadiationMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SurfaceEmissivityMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SoilWaterTranspirationMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/ResistanceGroundEvaporationMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/PsychrometricVariableMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SurfaceEnergyFluxVegetatedMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/VaporPressureSaturationMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/ResistanceAboveCanopyMostMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/ResistanceLeafToGroundMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/ResistanceCanopyStomataBallBerryMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SurfaceEnergyFluxBareGroundMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/ResistanceBareGroundMostMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SoilSnowTemperatureMainMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SoilSnowThermalDiffusionMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SoilSnowTemperatureSolverMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SoilSnowWaterPhaseChangeMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/CanopyHydrologyMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SnowWaterMainMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SnowfallBelowCanopyMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SnowpackCompactionMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SnowLayerCombineMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SnowLayerWaterComboMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SnowLayerDivideMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SnowpackHydrologyMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SoilWaterMainMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/RunoffSurfaceFreeDrainMod.F90
            # Added so RUNOFF_OPTION=7 works on the device. Without it the two
            # RunoffSurfaceXinAnJiang calls in SoilWaterMainMod sit inside
            # #ifndef NOAHMP_ACC_COLUMNS and are compiled out, so the GPU never
            # computes InfilRateSfc or RunoffSurface at all -- a systematic
            # one-signed dry bias in SOILSAT_TOP (mean 0.007, max 0.278).
            ${CMAKE_CURRENT_SOURCE_DIR}/src/RunoffSurfaceXinAnJiangMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SoilWaterDiffusionRichardsMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SoilMoistureSolverMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/RunoffSubSurfaceDrainageMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SoilHydraulicPropertyMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/MatrixSolverTriDiagonalMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SoilWaterSupercoolNiu06Mod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SoilWaterSupercoolKoren99Mod.F90
        )

        # -------------------------------------------------------------------
        # -O0 for the OpenACC sources.
        #
        # HARDCODED here rather than exposed as a cache option: re-running
        # cmake (x.sh does) kept resetting the cached value, and a mis-set
        # config produces a run that looks perfectly healthy while being
        # byte-identical to the previous one -- an easy way to waste a cycle.
        #
        # Why -O0 at all: SnowWaterMainMod and SoilWaterMainMod already
        # required it to work, i.e. Cray's OpenACC codegen at -O2 has
        # demonstrably produced wrong results in this code before.
        #
        # The "-O0 on ALL 52 sources" experiment has since been RUN and came
        # back NEGATIVE: exactly 1 of 906,542 differing LDASOUT values moved,
        # so the optimizer is exonerated and the GPU-vs-CPU bias is
        # source-level (it is the math shims). Hence this is back to the
        # original two files.
        #
        # DO NOT put -O0 on the whole list again without reverting the sqrt
        # change: -O0 stops the compiler inlining the hardware sqrt
        # instruction, so it emits a call to 'sqrtf' which has no device
        # implementation and nvlink fails. Measured on CCE 19.0.0:
        #     fp1 -O2 -> sqrt LINKS      fp1 -O0 -> sqrt FAILS (sqrtf)
        # Neither of the two files below calls sqrt, so they are safe at -O0.
        #
        # BISECTION (if ever needed again): replace the list below with an
        # explicit subset, halving each time. The configure step prints how
        # many sources got -O0; check that line, or
        # `grep -c -- "-O0" <builddir>/.../noahmp_submodule.dir/flags.make`,
        # before trusting a run.
        # -------------------------------------------------------------------
        set(NOAHMP_ACC_O0_SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SnowWaterMainMod.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/SoilWaterMainMod.F90)

        set_source_files_properties(
            ${NOAHMP_ACC_COLUMN_SOURCES}
            # fp1, NOT fp0. Measured on CCE 19.0.0: at -h fp0 the compiler
            # refuses the hardware sqrt instruction and emits a call to the
            # precise libm 'sqrtf', which has no device implementation, so
            # nvlink fails. At fp1/fp2 it inlines the instruction and links.
            # fp1 also matches the host, which the top-level CMakeLists
            # compiles at -h fp1 -- pinning these files at fp0 was a mismatch
            # with the CPU baseline, not a stricter-is-better choice.
            # (Only sqrt is recovered this way; exp/log/tanh/atan/acos/cos/tan
            # and x**y with a runtime exponent fail at EVERY fp level, which
            # is why NoahmpAccDeviceMathShimMod still has to exist.)
            PROPERTIES COMPILE_OPTIONS "-h;acc;-target-accel=nvidia80;-h;fp1"
        )
        set_source_files_properties(
            ${NOAHMP_ACC_O0_SOURCES}
            PROPERTIES COMPILE_OPTIONS "-h;acc;-target-accel=nvidia80;-h;fp1;-O0"
        )
        list(LENGTH NOAHMP_ACC_O0_SOURCES _noahmp_acc_o0_count)
        list(LENGTH NOAHMP_ACC_COLUMN_SOURCES _noahmp_acc_total)
        message("-- NoahMP ACC: -O0 applied to ${_noahmp_acc_o0_count} of ${_noahmp_acc_total} OpenACC sources")
    endif()
endif()

# Ensure free-form line length is unlimited for gfortran
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_compile_options(noahmp_submodule PRIVATE -ffree-line-length-none)
endif()
