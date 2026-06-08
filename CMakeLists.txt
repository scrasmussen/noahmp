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

# Ensure free-form line length is unlimited for gfortran
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_compile_options(noahmp_submodule PRIVATE -ffree-line-length-none)
endif()
