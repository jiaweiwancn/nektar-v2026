cmake_minimum_required(VERSION 3.16)
project(WJH27Feasibility CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(Nektar++ REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
get_filename_component(H27_SOURCE_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../../../.." ABSOLUTE)
add_definitions(${NEKTAR++_DEFINITIONS} -DENABLE_NEKTAR_EXCEPTIONS)
include_directories(${H27_SOURCE_ROOT}/library ${H27_SOURCE_ROOT}/solvers)
include_directories(SYSTEM ${NEKTAR++_INCLUDE_DIRS} ${NEKTAR++_TP_INCLUDE_DIRS})
# Build the exact historical adapter and pure WJ directional kernel afresh.
# This probe is a constitutive diagnostic, not a coupled EquationSystem.
add_subdirectory(../h26 h26-kernel EXCLUDE_FROM_ALL)
add_library(WJNativeMixedFlow STATIC
    ${H27_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/WJNativeMixedFlow.cpp
    ${H27_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/WJNativeWallV2.cpp
    ${H27_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/WJNativeRules.cpp
    ${H27_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/WJNativeScalar.cpp
    ${H27_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/WJNativeReference.cpp
    ${H27_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/WJNativeAdmissibility.cpp
    ${H27_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/WJNativeCoupledMMS.cpp
    ${H27_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/WJNativeCoupling.cpp
    ${H27_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/WJNativeFlowTangent.cpp
    ${H27_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/WJNativeTimeNewton.cpp
    ${H27_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/WJNativeConstitutive.cpp
    ${H27_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/WJNativePerturbation.cpp
    ${H27_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/WJNativeSteadyKrylov.cpp
    ${H27_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/WJNativeInitialFields.cpp
    # The isolated prefix predates these legacy factory dependencies. Build
    # their exact capsule sources rather than dropping historical factories.
    ${H27_SOURCE_ROOT}/library/SolverUtils/Turbulence/RANSRobustVariant.cpp
    ${H27_SOURCE_ROOT}/library/SolverUtils/Turbulence/RANSWallInnerStress.cpp
    ${H27_SOURCE_ROOT}/library/SolverUtils/Turbulence/EquilibriumAWFProfiles.cpp
    ${H27_SOURCE_ROOT}/library/SolverUtils/Turbulence/EquilibriumAWFTransition.cpp
    ${H27_SOURCE_ROOT}/library/SolverUtils/Turbulence/RANSWallGeometry.cpp)
target_compile_definitions(WJNativeMixedFlow PRIVATE SOLVER_UTILS_EXPORTS)
target_link_libraries(WJNativeMixedFlow PUBLIC WJNativeFlow)
target_compile_options(WJNativeMixedFlow PRIVATE
    -O3 -ffp-contract=off -fno-fast-math -Wall -Wextra -Wno-sign-compare)
add_executable(WJH27PrincipalProbe principal.cpp)
target_link_libraries(WJH27PrincipalProbe PRIVATE WJNativeMixedFlow)
target_compile_options(WJH27PrincipalProbe PRIVATE
    -O3 -ffp-contract=off -fno-fast-math -Wall -Wextra)
add_executable(WJH27MixedProbe mixed_probe.cpp)
target_link_libraries(WJH27MixedProbe PRIVATE WJNativeMixedFlow)
target_compile_options(WJH27MixedProbe PRIVATE
    -O3 -ffp-contract=off -fno-fast-math -Wall -Wextra -Wno-sign-compare)
add_executable(WJH27ScalarProbe scalar_probe.cpp)
target_link_libraries(WJH27ScalarProbe PRIVATE WJNativeMixedFlow)
target_compile_options(WJH27ScalarProbe PRIVATE
    -O3 -ffp-contract=off -fno-fast-math -Wall -Wextra -Wno-sign-compare)
include(source_identity.cmake)
wj_h27_source_identity("${H27_SOURCE_ROOT}" wj_source_identity)
target_compile_definitions(WJNativeMixedFlow PUBLIC WJ_H27_SOURCE_ID="${wj_source_identity}")
foreach(probe reference admissibility coupled_mms constitutive steady_krylov initial_fields perturbation flow_tangent)
    if(probe STREQUAL "reference")
        set(name WJH27ReferenceProbe)
    elseif(probe STREQUAL "admissibility")
        set(name WJH27AdmissibilityProbe)
    elseif(probe STREQUAL "coupled_mms")
        set(name WJH27CoupledMMSProbe)
    elseif(probe STREQUAL "constitutive")
        set(name WJH27ConstitutiveProbe)
    elseif(probe STREQUAL "steady_krylov")
        set(name WJH27SteadyKrylovProbe)
    elseif(probe STREQUAL "initial_fields")
        set(name WJH27InitialFieldsProbe)
    elseif(probe STREQUAL "perturbation")
        set(name WJH27PerturbationProbe)
    else()
        set(name WJH27FlowTangentProbe)
    endif()
    add_executable(${name} "${probe}_probe.cpp")
    target_link_libraries(${name} PRIVATE WJNativeMixedFlow)
    target_compile_options(${name} PRIVATE -O3 -ffp-contract=off -fno-fast-math -Wall -Wextra)
endforeach()
add_executable(WJH27TraceGeometryProbe trace_geometry_probe.cpp)
target_link_libraries(WJH27TraceGeometryProbe PRIVATE WJNativeMixedFlow)
target_compile_options(WJH27TraceGeometryProbe PRIVATE
    -O3 -ffp-contract=off -fno-fast-math -Wall -Wextra)
add_executable(WJH27CutoffPanelProbe cutoff_panel_probe.cpp)
target_link_libraries(WJH27CutoffPanelProbe PRIVATE WJNativeMixedFlow)
target_compile_options(WJH27CutoffPanelProbe PRIVATE
    -O3 -ffp-contract=off -fno-fast-math -Wall -Wextra)
add_executable(WJH27CutoffPanelsUnit cutoff_panels_unit.cpp)
target_compile_options(WJH27CutoffPanelsUnit PRIVATE
    -O3 -ffp-contract=off -fno-fast-math -Wall -Wextra)
# The transport diagnostic has its own main and requires a fresh steady
# physical session. It is deliberately absent from the production factory.
file(SHA256 "${CMAKE_CURRENT_LIST_DIR}/inner_k_fixture.cpp" inner_k_main)
string(SHA256 inner_k_source "WJInnerKFixture|${wj_source_identity}|${inner_k_main}")
add_executable(WJH27InnerKFixture inner_k_fixture.cpp)
target_link_libraries(WJH27InnerKFixture PRIVATE WJNativeMixedFlow)
target_compile_definitions(WJH27InnerKFixture PRIVATE WJ_INNER_K_FIXTURE_SOURCE_ID="${inner_k_source}")
target_compile_options(WJH27InnerKFixture PRIVATE
    -O3 -ffp-contract=off -fno-fast-math -Wall -Wextra)
# The real production main, Standard Driver and complete existing factory
# source inventory. No probe is relabelled as the application executable.
set(h27_solver_root "${H27_SOURCE_ROOT}/solvers/IncNavierStokesSolver")
file(READ "${h27_solver_root}/CMakeLists.txt" h27_solver_cmake)
string(REGEX MATCH "LIBRARY_SOURCES[^)]*" h27_library_block "${h27_solver_cmake}")
string(REGEX MATCHALL "\\./[A-Za-z0-9_/]+\\.cpp" h27_solver_sources "${h27_library_block}")
list(TRANSFORM h27_solver_sources PREPEND "${h27_solver_root}/")
add_executable(IncNavierStokesSolver "${h27_solver_root}/IncNavierStokesSolver.cpp"
    "${h27_solver_root}/EquationSystems/WJNativeCoupled.cpp" ${h27_solver_sources})
target_link_libraries(IncNavierStokesSolver PRIVATE WJNativeMixedFlow)
target_compile_options(IncNavierStokesSolver PRIVATE -O3 -ffp-contract=off -fno-fast-math -Wno-sign-compare)

add_executable(WJH27CoupledReturnProbe coupled_return_probe.cpp)
target_link_libraries(WJH27CoupledReturnProbe PRIVATE WJNativeMixedFlow)
target_compile_options(WJH27CoupledReturnProbe PRIVATE
    -O3 -ffp-contract=off -fno-fast-math -Wall -Wextra)
